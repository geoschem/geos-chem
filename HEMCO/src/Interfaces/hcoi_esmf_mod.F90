!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_esmf_mod
!
! !DESCRIPTION: Module HCOI\_ESMF\_MOD is the HEMCO-ESMF interface module.
!\\
!\\
! !INTERFACE:
!
MODULE HCOI_ESMF_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_Types_Mod

#if defined (ESMF_)
#include "MAPL_Generic.h"
  USE ESMF
  USE MAPL_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!

  ! ESMF environment only:
  PUBLIC :: HCO_SetServices
  PUBLIC :: HCO_SetExtState_ESMF
  PUBLIC :: HCO_Imp2Ext
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Diagn2Exp
  PRIVATE :: HCO_Imp2Ext2R
  PRIVATE :: HCO_Imp2Ext2S
  PRIVATE :: HCO_Imp2Ext2I
  PRIVATE :: HCO_Imp2Ext3R
  PRIVATE :: HCO_Imp2Ext3S
!
! !REVISION HISTORY:
!  10 Oct 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE HCO_Imp2Ext
     MODULE PROCEDURE HCO_Imp2Ext2R
     MODULE PROCEDURE HCO_Imp2Ext2S
     MODULE PROCEDURE HCO_Imp2Ext2I
     MODULE PROCEDURE HCO_Imp2Ext3R
     MODULE PROCEDURE HCO_Imp2Ext3S
  END INTERFACE HCO_Imp2Ext

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_SetServices
!
! !DESCRIPTION: Subroutine HCO\_SetServices registers all required HEMCO
! data so that it can be imported through the ESMF import state.
! This routine determines all required HEMCO input fields from the HEMCO
! configuration file. Note that each file needs an equivalent ESMF-style
! entry in the registry file (typically ExtData.rc). Otherwise, ESMF won't
! read these files and HEMCO will fail when attempting to get pointers to
! these data arrays.
!\\
!\\
! The field names provided in ExtData.rc must match the names in the HEMCO
! configuration file! Also, all time settings (average and update interval)
! and data units need to be properly specified in ExtData.rc.
! For now, ExtData.rc and HEMCO configuration file need to be synchronized
! manually. The pyHEMCO interface will automate this process!
!\\
!\\
! This routine also prepares an emissions export field for every species
! found in the HEMCO configuration file. These export fields will only
! be filled if specified so in the MAPL History registry.
! The corresponding HEMCO diagnostics must be created separately via
! Diagn\_Create (e.g. in hcoi\_gc\_diagn\_mod.F90).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_SetServices( am_I_Root,  GC, HcoConfig, &
                                  ConfigFile, RC )
!
! !USES:
!
      USE HCO_TYPES_MOD,    ONLY : ListCont
      USE HCO_DATACONT_MOD, ONLY : ListCont_NextCont
      USE HCO_CONFIG_MOD,   ONLY : Config_ReadFile
      USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
      USE HCO_CONFIG_MOD,   ONLY : Config_GetnSpecies
      USE HCO_CONFIG_MOD,   ONLY : Config_GetSpecNames
      USE HCO_DIAGN_MOD,    ONLY : DiagnFileOpen
      USE HCO_DIAGN_MOD,    ONLY : DiagnFileGetNext
      USE HCO_DIAGN_MOD,    ONLY : DiagnFileClose
!
! !ARGUMENTS:
!
      LOGICAL,             INTENT(IN   )             :: am_I_Root
      TYPE(ESMF_GridComp), INTENT(INOUT)             :: GC
      TYPE(ConfigObj),     POINTER                   :: HcoConfig
      CHARACTER(LEN=*),    INTENT(IN   )             :: ConfigFile
      INTEGER,             INTENT(  OUT)             :: RC
!
! !REVISION HISTORY:
!  29 Aug 2013 - C. Keller - Initial version.
!  10 Sep 2015 - C. Keller - Added RESTART=MAPL_RestartSkip.
!  21 Feb 2016 - C. Keller - Update to v2.0, added default diagnostics (optional)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  16 Mar 2018 - E. Lundgren - Expand log write to specify reading HEMCO
!                              diagnostic config file and adding HEMCO exports
!  16 Jul 2018 - E. Lundgren - Move verbose to within DoUse blocks
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                    :: LUN, ExtNr, Cat, Hier, SpaceDim
      INTEGER                    :: DIMS, VLOC
      INTEGER                    :: I, FLAG, nSpc, nDiagn
      INTEGER                    :: DefaultDim
      LOGICAL                    :: EOF
      LOGICAL                    :: FOUND, DefaultSet
      CHARACTER(LEN=31)          :: cName, SpcName, OutUnit
      CHARACTER(LEN=63)          :: DefaultSNAME, DefaultLNAME, DefaultUnit
      CHARACTER(LEN=63)          :: SNAME, UnitName
      CHARACTER(LEN=127)         :: LNAME
      CHARACTER(LEN=63), POINTER :: Spc(:)
      TYPE(ListCont),    POINTER :: CurrCont

      ! ================================================================
      ! HCO_SetServices begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('HCO_SetServices (HCOI_ESMF_MOD.F90)')

      ! Init
      Spc      => NULL()
      CurrCont => NULL()

      ! ---------------------------------------------------------------------
      ! Read file into buffer
      ! ---------------------------------------------------------------------

      CALL Config_ReadFile( am_I_Root, HcoConfig, TRIM(ConfigFile), 0, STATUS )
      ASSERT_(STATUS==HCO_SUCCESS)

      ! ---------------------------------------------------------------------
      ! Set services for all import fields
      ! ---------------------------------------------------------------------

      ! Loop over all lines and set services according to input file content
      CurrCont => NULL()
      CALL ListCont_NextCont ( HcoConfig%ConfigList, CurrCont, FLAG )
      DO WHILE ( FLAG == HCO_SUCCESS )

         ! Skip containers that are not defined
         IF ( .NOT. ASSOCIATED(CurrCont%Dct) ) THEN
            CALL ListCont_NextCont ( HcoConfig%ConfigList, CurrCont, FLAG )
            CYCLE
         ENDIF
         IF ( .NOT. ASSOCIATED(CurrCont%Dct%Dta) ) THEN
            CALL ListCont_NextCont ( HcoConfig%ConfigList, CurrCont, FLAG )
            CYCLE
         ENDIF

         ! Add arrays to import spec. Distinguish between 2D and 3D arrays.
         ! Note that we can ignore the time reading interval here, as this
         ! is automatically determined by ESMF based upon the registry file
         ! content!.
         ! Ignore containers with ncRead flag disabled. These are typically
         ! scalar fields directly read from the configuration file.
         IF ( .NOT. CurrCont%Dct%Dta%ncRead ) THEN

         ! Multiple data containers can use the same source data. In this
         ! case we only need to import the data once. The second, third, etc.
         ! containters registerd for the same source data have been assigned
         ! lower DtaHome values (in hco_config_mod.F90), so skip this container
         ! if flag is not -999 (=default).
         ELSEIF ( CurrCont%Dct%DtaHome /= -999 ) THEN

         ! Import 2D data
         ELSEIF ( CurrCont%Dct%Dta%SpaceDim == 2 ) THEN

            CALL MAPL_AddImportSpec(GC,                      &
               SHORT_NAME = TRIM(CurrCont%Dct%Dta%ncFile),   &
               LONG_NAME  = TRIM(CurrCont%Dct%Dta%ncFile),   &
               UNITS      = TRIM(CurrCont%Dct%Dta%OrigUnit), &
               DIMS       = MAPL_DimsHorzOnly,               &
               VLOCATION  = MAPL_VLocationNone,              &
               RESTART    = MAPL_RestartSkip,                &
               RC         = STATUS                            )

            ! Error trap
            IF ( STATUS /= ESMF_SUCCESS ) THEN
               WRITE(*,*) '2D import error: ', TRIM(CurrCont%Dct%Dta%ncFile)
               VERIFY_(STATUS)
            ENDIF

         ! Import 3D data: Assume central location in vertical dimension!
         ELSEIF ( CurrCont%Dct%Dta%SpaceDim == 3 ) THEN

            CALL MAPL_AddImportSpec(GC,                      &
               SHORT_NAME = TRIM(CurrCont%Dct%Dta%ncFile),   &
               LONG_NAME  = TRIM(CurrCont%Dct%Dta%ncFile),   &
               UNITS      = TRIM(CurrCont%Dct%Dta%OrigUnit), &
               DIMS       = MAPL_DimsHorzVert,               &
               VLOCATION  = MAPL_VLocationCenter,            &
               RESTART    = MAPL_RestartSkip,                &
               RC         = STATUS                            )

            ! Error trap
            IF ( STATUS /= ESMF_SUCCESS ) THEN
               WRITE(*,*) '3D import error: ', TRIM(CurrCont%Dct%Dta%ncFile)
               VERIFY_(STATUS)
            ENDIF

         ! Return w/ error if not 2D or 3D data
         ELSE
            ASSERT_(.FALSE.)
         ENDIF

         ! Advance to next container
         CALL ListCont_NextCont ( HcoConfig%ConfigList, CurrCont, FLAG )

      ENDDO

      ! Free pointer
      CurrCont => NULL()

      ! ---------------------------------------------------------------------
      ! Try to open diagnostics definition file
      ! ---------------------------------------------------------------------
      CALL DiagnFileOpen( HcoConfig, LUN, RC )
      ASSERT_(RC == HCO_SUCCESS )

      ! ---------------------------------------------------------------------
      ! If DiagnFile is found, prepare a diagnostics export for every entry
      ! ---------------------------------------------------------------------

      IF ( LUN > 0 ) THEN

         IF ( am_I_Root ) WRITE(*,*) 'Reading HEMCO configuration file: ', &
                                     TRIM(HcoConfig%ConfigFileName)
         DO

            ! Get next line
            CALL DiagnFileGetNext( HcoConfig, LUN,     cName,       &
                                   SpcName,   ExtNr,   Cat,   Hier, &
                                   SpaceDim,  OutUnit, EOF,   RC,   &
                                   lName=lName, UnitName=UnitName )
            IF ( RC /= HCO_SUCCESS ) RETURN

            ! Leave here if end of file
            IF ( EOF ) EXIT

            ! Define vertical dimension
            IF ( SpaceDim == 3 ) THEN
               DIMS = MAPL_DimsHorzVert
               VLOC = MAPL_VLocationCenter
            ELSE
               DIMS = MAPL_DimsHorzOnly
               VLOC = MAPL_VLocationNone
            ENDIF

            ! Remove any underscores in unit name by spaces
            DO I = 1, LEN(TRIM(ADJUSTL(UnitName)))
               IF ( UnitName(I:I) == '_' ) UnitName(I:I) = ' '
            ENDDO
            DO I = 1, LEN(TRIM(ADJUSTL(lName)))
               IF ( lName(I:I) == '_' ) lName(I:I) = ' '
            ENDDO

            ! Add to export state
            CALL MAPL_AddExportSpec(GC,       &
               SHORT_NAME         = cName,    &
               LONG_NAME          = lName,    &
               UNITS              = UnitName, &
               DIMS               = DIMS,     &
               VLOCATION          = VLOC,     &
               RC                 = STATUS     )
            IF ( STATUS /= ESMF_SUCCESS ) THEN
               WRITE(*,*) 'Cannot add to export: ',TRIM(SNAME)
               ASSERT_(.FALSE.)
            ELSE
               IF ( am_I_Root ) WRITE(*,*) 'adding HEMCO export: ', TRIM(cName)
            ENDIF

         ENDDO

         ! Close file
         CALL DiagnFileClose ( LUN )
      ENDIF

      ! ---------------------------------------------------------------------
      ! Eventually prepare a diagnostics export for every potential HEMCO
      ! species. This is optional and controlled by HEMCO setting
      ! DefaultDiagnSet.
      ! ---------------------------------------------------------------------
      CALL GetExtOpt( HcoConfig, -999, 'DefaultDiagnOn', &
                      OptValBool=DefaultSet, FOUND=FOUND, RC=RC )
      IF ( .NOT. FOUND ) DefaultSet = .FALSE.
      IF ( DefaultSet ) THEN

         ! Search for default diagnostics variable prefix
         CALL GetExtOpt( HcoConfig, -999, 'DefaultDiagnSname', &
                         OptValChar=DefaultSNAME, FOUND=FOUND, RC=RC )
         IF ( .NOT. FOUND ) DefaultSNAME = 'HEMCO_EMIS_'

         CALL GetExtOpt( HcoConfig, -999, 'DefaultDiagnLname', &
                         OptValChar=DefaultLNAME, FOUND=FOUND, RC=RC )
         IF ( .NOT. FOUND ) DefaultLNAME = 'HEMCO_emissions_of_species_'

         ! Search for default diagnostics dimension
         CALL GetExtOpt( HcoConfig, -999, 'DefaultDiagnDim', &
                         OptValInt=DefaultDim, FOUND=FOUND, RC=RC )
         IF ( .NOT. FOUND ) DefaultDim = 3
         DefaultDim = MAX(MIN(DefaultDim,3),2)

         ! Get units
         CALL GetExtOpt( HcoConfig, -999, 'DefaultDiagnUnit', &
                         OptValChar=DefaultUnit, FOUND=FOUND, RC=RC )
         IF ( .NOT. FOUND ) DefaultUnit = 'kg m-2 s-1'

         ! Get # of species and species names
         nSpc = Config_GetnSpecies( HcoConfig )
         CALL Config_GetSpecNames( HcoConfig, Spc, nSpc, RC )
         ASSERT_(RC == HCO_SUCCESS)

         ! Loop over all species and add to export state
         DO I = 1, nSpc
            SNAME = TRIM(DefaultSNAME)//TRIM(Spc(I))
            LNAME = TRIM(DefaultLNAME)//TRIM(Spc(I))
            CALL Diagn2Exp( GC, SNAME, LNAME, DefaultUnit, DefaultDim, __RC__ )
         ENDDO
      ENDIF

      ! ---------------------------------------------------------------------
      ! Cleanup
      ! ---------------------------------------------------------------------
      IF ( ASSOCIATED(Spc) ) DEALLOCATE(Spc)

      ! Return success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCO_SetServices
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn2Exp
!
! !DESCRIPTION: Subroutine Diagn2Exp is a helper routine to add a potential
! HEMCO diagnostics to the Export state.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Diagn2Exp( GC, SNAME, LNAME, UNITS, NDIM, RC )
!
! !USES:
!
!
! !ARGUMENTS:
!
      TYPE(ESMF_GridComp), INTENT(INOUT)   :: GC
      CHARACTER(LEN=*),    INTENT(IN   )   :: SNAME
      CHARACTER(LEN=*),    INTENT(IN   )   :: LNAME
      CHARACTER(LEN=*),    INTENT(IN   )   :: UNITS
      INTEGER,             INTENT(IN   )   :: NDIM
      INTEGER,             INTENT(  OUT)   :: RC
!
! !REVISION HISTORY:
!  05 Jan 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: DIMS, VLOC

      ! ================================================================
      ! Diagn2Exp begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('Diagn2Exp (HCOI_ESMF_MOD.F90)')

      ! Set horizontal variables
      IF ( NDIM == 3 ) THEN
         DIMS = MAPL_DimsHorzVert
         VLOC = MAPL_VLocationCenter
      ELSE
         DIMS = MAPL_DimsHorzOnly
         VLOC = MAPL_VLocationNone
      ENDIF

      ! Add to export state
      CALL MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = SNAME,                             &
         LONG_NAME          = LNAME,                            &
         UNITS              = UNITS,                            &
         DIMS               = DIMS,                             &
         VLOCATION          = VLOC,                             &
         RC                 = STATUS                             )
      IF ( STATUS /= ESMF_SUCCESS ) THEN
         WRITE(*,*) 'Cannot add to export: ',TRIM(SNAME)
         ASSERT_(.FALSE.)
      ENDIF

      ! Return w/ success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE Diagn2Exp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_SetExtState_ESMF
!
! !DESCRIPTION: Subroutine HCO\_SetExtState\_ESMF tries to populate some
! fields of the ExtState object from the ESMF import state.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_SetExtState_ESMF( HcoState, ExtState, RC )
!
! !USES:
!
      USE HCO_STATE_MOD,   ONLY : Hco_State
      USE HCOX_STATE_MOD,  ONLY : Ext_State
!
! !ARGUMENTS:
!
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(Ext_State),     POINTER         :: ExtState
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  06 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! ================================================================
      ! HCO_SetExtState_ESMF begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('HCO_SetExtState_ESMF (HCOI_ESMF_MOD.F90)')

      ! Get pointers to fields
      CALL HCO_Imp2Ext ( HcoState, ExtState%BYNCY, 'BYNCY', __RC__ )

#if defined( MODEL_GEOS )
      ! Get pointers to fields
      CALL HCO_Imp2Ext ( HcoState, ExtState%LFR, 'LFR', __RC__ )
#endif

      ! Return success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_SetExtState_ESMF
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext2S
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext copies fields from the import state to
! the HEMCO ExtState object.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext2S( HcoState, ExtDat, FldName, RC )
!
! !USES:
!
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
      USE HCO_STATE_MOD,   ONLY : Hco_State
      USE HCOX_STATE_MOD,  ONLY : ExtDat_2S
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_2S),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  08 Feb 2016 - C. Keller - Initial version
!  11 Apr 2017 - C. Keller - It's now ok if field not found.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG
      REAL,             POINTER    :: Ptr2D(:,:)   => NULL()
      INTEGER                      :: STAT

      ! ================================================================
      ! HCO_Imp2Ext2S begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('HCO_Imp2Ext2S (HCOI_ESMF_MOD.F90)')

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN
         ASSERT_( ASSOCIATED(HcoState%IMPORT) )
         CALL MAPL_GetPointer( HcoState%IMPORT, Ptr2D, TRIM(FldName), __RC__ )
         CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
         ASSERT_(STAT==HCO_SUCCESS)
         ExtDat%Arr%Val = 0.0
         IF ( ASSOCIATED( Ptr2D ) ) THEN
            WHERE( Ptr2D /= MAPL_UNDEF )
               ExtDat%Arr%Val = Ptr2D
            END WHERE
         ENDIF
         Ptr2D => NULL()

         ! Verbose
         IF ( HCO_IsVerb(HcoState%Config%Err,2) .AND.HcoState%amIRoot ) THEN
            CALL HCO_MSG('Passed from import to ExtState: '//TRIM(FldName))
         ENDIF

      ENDIF ! DoUse

      ! Return success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCO_Imp2Ext2S
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext3S
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext copies fields from the import state to
! the HEMCO ExtState object.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext3S( HcoState, ExtDat, FldName, RC )
!
! !USES:
!
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
      USE HCO_STATE_MOD,   ONLY : Hco_State
      USE HCOX_STATE_MOD,  ONLY : ExtDat_3S
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_3S),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  08 Feb 2016 - C. Keller - Initial version
!  11 Apr 2017 - C. Keller - It's now ok if field not found.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG
      REAL,             POINTER    :: Ptr3D(:,:,:)   => NULL()
      INTEGER                      :: L, NZ, OFF, STAT

      ! ================================================================
      ! HCO_Imp2Ext3S begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('HCO_Imp2Ext3S (HCOI_ESMF_MOD.F90)')

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN

         ASSERT_( ASSOCIATED(HcoState%IMPORT) )
         CALL MAPL_GetPointer( HcoState%IMPORT, Ptr3D, TRIM(FldName), __RC__ )
         ASSERT_( ASSOCIATED(Ptr3D) )

         ! Make sure the array in ExtDat is allocated and has the right size
         NZ = SIZE(Ptr3D,3)
         CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ, STAT )
         ASSERT_(STAT==HCO_SUCCESS)

         ! Pass field to ExtDat
         OFF = LBOUND(Ptr3D,3)
         DO L = 1, NZ
            WHERE ( Ptr3D(:,:,NZ-L+OFF) == MAPL_UNDEF )
               ExtDat%Arr%Val(:,:,L) = 0.0
            ELSEWHERE
               ExtDat%Arr%Val(:,:,L) = Ptr3D(:,:,NZ-L+OFF)
            END WHERE
         ENDDO
         Ptr3D => NULL()

         ! Verbose
         IF ( HCO_IsVerb(HcoState%Config%Err,2) .AND. HcoState%amIRoot ) THEN
            CALL HCO_MSG('Passed from import to ExtState: '//TRIM(FldName))
         ENDIF

      ENDIF ! DoUse


      ! Return success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCO_Imp2Ext3S
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext2R
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext copies fields from the import state to
! the HEMCO ExtState object.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext2R( HcoState, ExtDat, FldName, RC, Fld )
!
! !USES:
!
      USE HCO_STATE_MOD,   ONLY : Hco_State
      USE HCOX_STATE_MOD,  ONLY : ExtDat_2R
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_2R),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
      REAL(hp), OPTIONAL,  INTENT(IN)      :: Fld(HcoState%NX,HcoState%NY)
!
! !REVISION HISTORY:
!  08 Feb 2016 - C. Keller - Initial version
!  11 Apr 2017 - C. Keller - It's now ok if field not found.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG
      INTEGER                      :: STAT
      REAL,             POINTER    :: Ptr2D(:,:)   => NULL()
      LOGICAL                      :: Filled

      ! ================================================================
      ! HCO_Imp2Ext2R begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('HCO_Imp2Ext2R (HCOI_ESMF_MOD.F90)')

      ! Init
      Filled = .FALSE.

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN

         ASSERT_( ASSOCIATED(HcoState%IMPORT) )
         CALL MAPL_GetPointer( HcoState%IMPORT, Ptr2D, TRIM(FldName), NotFoundOk=.TRUE., __RC__ )

         CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
         ASSERT_(STAT==HCO_SUCCESS)

         ExtDat%Arr%Val = 0.0
         IF ( ASSOCIATED( Ptr2D ) ) THEN
            WHERE( Ptr2D /= MAPL_UNDEF )
               ExtDat%Arr%Val = Ptr2D
            END WHERE
            Filled = .TRUE.
         ELSE
            IF ( PRESENT(Fld) ) THEN
               ExtDat%Arr%Val = Fld
               Filled = .TRUE.
            ENDIF
         ENDIF
         Ptr2D => NULL()

         ! Error check
         IF ( .NOT. Filled ) THEN
            CALL HCO_ERROR(HcoState%Config%Err,'Cannot fill '//TRIM(FldName),RC)
            ASSERT_(.FALSE.)
         ENDIF

         ! Verbose
         IF ( HCO_IsVerb(HcoState%Config%Err,2) .AND. HcoState%amIRoot ) THEN
            CALL HCO_MSG(HcoState%Config%Err,'Passed from import to ExtState: '//TRIM(FldName))
         ENDIF

      ENDIF ! DoUse

      ! Return success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCO_Imp2Ext2R
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext3R
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext copies fields from the import state to
! the HEMCO ExtState object.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext3R( HcoState, ExtDat, FldName, RC )
!
! !USES:
!
      USE HCO_STATE_MOD,   ONLY : Hco_State
      USE HCOX_STATE_MOD,  ONLY : ExtDat_3R
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_3R),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  08 Feb 2016 - C. Keller - Initial version
!  11 Apr 2017 - C. Keller - It's now ok if field not found.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG
      INTEGER                      :: L, NZ, OFF, STAT
      REAL,             POINTER    :: Ptr3D(:,:,:) => NULL()

      ! ================================================================
      ! HCO_Imp2Ext3R begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('HCO_Imp2Ext3R (HCOI_ESMF_MOD.F90)')

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN

         ASSERT_( ASSOCIATED(HcoState%IMPORT) )
         CALL MAPL_GetPointer( HcoState%IMPORT, Ptr3D, TRIM(FldName), __RC__ )
         ASSERT_( ASSOCIATED(Ptr3D) )

         ! Make sure the array in ExtDat is allocated and has the right size
         NZ = SIZE(Ptr3D,3)
         CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ, STAT )
         ASSERT_(STAT==HCO_SUCCESS)

         ! Pass field to ExtDat
         OFF = LBOUND(Ptr3D,3)
         DO L = 1, NZ
            WHERE ( Ptr3D(:,:,NZ-L+OFF) == MAPL_UNDEF )
               ExtDat%Arr%Val(:,:,L) = 0.0
            ELSEWHERE
               ExtDat%Arr%Val(:,:,L) = Ptr3D(:,:,NZ-L+OFF)
            END WHERE
         ENDDO
         Ptr3D => NULL()

         ! Verbose
         IF ( HCO_IsVerb(HcoState%Config%Err,2) .AND. HcoState%amIRoot ) THEN
            CALL HCO_MSG('Passed from import to ExtState: '//TRIM(FldName))
         ENDIF

      ENDIF ! DoUse

      ! Return success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCO_Imp2Ext3R
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext2I
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext copies fields from the import state to
! the HEMCO ExtState object.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext2I( HcoState, ExtDat, FldName, RC )
!
! !USES:
!
      USE HCO_STATE_MOD,   ONLY : Hco_State
      USE HCOX_STATE_MOD,  ONLY : ExtDat_2I
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_2I),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  08 Feb 2016 - C. Keller - Initial version
!  11 Apr 2017 - C. Keller - It's now ok if field not found.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG
      INTEGER                      :: STAT
      REAL,             POINTER    :: Ptr2D(:,:)   => NULL()

      ! ================================================================
      ! HCO_Imp2Ext2I begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('HCO_Imp2Ext2I (HCOI_ESMF_MOD.F90)')

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN

         ASSERT_( ASSOCIATED(HcoState%IMPORT) )
         CALL MAPL_GetPointer( HcoState%IMPORT, Ptr2D, TRIM(FldName), __RC__ )

         CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
         ASSERT_(STAT==HCO_SUCCESS)

         ExtDat%Arr%Val = 0.0
         IF ( ASSOCIATED( Ptr2D ) ) THEN
            WHERE( Ptr2D /= MAPL_UNDEF )
               ExtDat%Arr%Val = Ptr2D
            END WHERE
         ENDIF
         Ptr2D => NULL()

         ! Verbose
         IF ( HCO_IsVerb(HcoState%Config%Err,2) .AND. HcoState%amIRoot ) THEN
            CALL HCO_MSG('Passed from import to ExtState: '//TRIM(FldName))
         ENDIF

      ENDIF ! DoUse

      ! Return success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCO_Imp2Ext2I
!EOC
#endif
END MODULE HCOI_ESMF_MOD
