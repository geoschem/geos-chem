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
!
! !PRIVATE MEMBER FUNCTIONS:
!      
  PRIVATE :: Diagn2Exp
!
! !REVISION HISTORY:
!  10 Oct 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
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
      SUBROUTINE HCO_SetServices( am_I_Root, GC, ConfigFile, RC )
!
! !USES:
!
      USE HCO_DATACONT_MOD, ONLY : ListCont
      USE HCO_CONFIG_MOD,   ONLY : Config_ReadFile, GetNextCont
      USE HCO_CONFIG_MOD,   ONLY : Config_ScalIDinUse
      USE HCO_CONFIG_MOD,   ONLY : Config_GetnSpecies
      USE HCO_CONFIG_MOD,   ONLY : Config_GetSpecNames
      USE HCO_DIAGN_MOD,    ONLY : DiagnFileOpen
      USE HCO_DIAGN_MOD,    ONLY : DiagnFileGetNext
      USE HCO_DIAGN_MOD,    ONLY : DiagnFileClose
!
! !ARGUMENTS:
!
      LOGICAL,             INTENT(IN   )   :: am_I_Root
      TYPE(ESMF_GridComp), INTENT(INOUT)   :: GC
      CHARACTER(LEN=*),    INTENT(IN   )   :: ConfigFile
      INTEGER,             INTENT(  OUT)   :: RC
!
! !REVISION HISTORY:
!  29 Aug 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                    :: LUN, ExtNr, Cat, Hier, SpaceDim
      INTEGER                    :: DIMS, VLOC
      INTEGER                    :: I, FLAG, nSpc, nDiagn
      LOGICAL                    :: EOF
      CHARACTER(LEN=31)          :: cName, SpcName, OutUnit 
      CHARACTER(LEN=63)          :: SNAME, LNAME, UNITS
      CHARACTER(LEN=63), POINTER :: Spc(:)   => NULL() 
      TYPE(ListCont),    POINTER :: CurrCont => NULL()

      ! ================================================================
      ! HCO_SetServices begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('HCO_SetServices (HCOI_ESMF_MOD.F90)')

      ! ---------------------------------------------------------------------
      ! Read file into buffer
      ! ---------------------------------------------------------------------

      CALL Config_ReadFile( am_I_Root, TRIM(ConfigFile), 0, STATUS )
      ASSERT_(STATUS==HCO_SUCCESS)

      ! ---------------------------------------------------------------------
      ! Set services for all import fields
      ! ---------------------------------------------------------------------

      ! Loop over all lines and set services according to input file content
      CurrCont => NULL()
      CALL GetNextCont ( CurrCont, FLAG )
      DO WHILE ( FLAG == HCO_SUCCESS )

         ! Skip containers that are not defined
         IF ( .NOT. ASSOCIATED(CurrCont%Dct) ) THEN
            CALL GetNextCont ( CurrCont, FLAG )
            CYCLE
         ENDIF
         IF ( .NOT. ASSOCIATED(CurrCont%Dct%Dta) ) THEN
            CALL GetNextCont ( CurrCont, FLAG )
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
         CALL GetNextCont ( CurrCont, FLAG )

      ENDDO

      ! Free pointer
      CurrCont => NULL()

      ! ---------------------------------------------------------------------
      ! Try to open diagnostics definition file 
      ! ---------------------------------------------------------------------
      CALL DiagnFileOpen( am_I_Root, LUN, RC )
      ASSERT_(RC == HCO_SUCCESS )

      ! ---------------------------------------------------------------------
      ! If DiagnFile is found, prepare a diagnostics export for every entry
      ! ---------------------------------------------------------------------

      IF ( LUN > 0 ) THEN

         DO 

            ! Get next line
            CALL DiagnFileGetNext( am_I_Root, LUN, cName, &
               SpcName, ExtNr, Cat, Hier, SpaceDim, OutUnit, EOF, RC ) 
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
  
            ! Add to export state 
            CALL MAPL_AddExportSpec(GC,       &
               SHORT_NAME         = cName,    &
               LONG_NAME          = cName,    &
               UNITS              = OutUnit,  &
               DIMS               = DIMS,     &
               VLOCATION          = VLOC,     &
               RC                 = STATUS     )
            IF ( STATUS /= ESMF_SUCCESS ) THEN
               WRITE(*,*) 'Cannot add to export: ',TRIM(SNAME)
               ASSERT_(.FALSE.)
            ENDIF

         ENDDO

         ! Close file
         CALL DiagnFileClose ( LUN )

      ! ---------------------------------------------------------------------
      ! If DiagnFile is is not found, prepare a diagnostics export for every
      ! potential HEMCO species
      ! ---------------------------------------------------------------------
      ELSE
         nSpc = Config_GetnSpecies( )
         CALL Config_GetSpecNames( Spc, nSpc, RC )
         ASSERT_(RC == HCO_SUCCESS) 

         ! All units in kg/m2/s
         UNITS = 'kg m-2 s-1'

         ! Loop over all species and add to export state
         DO I = 1, nSpc
            SNAME = 'EMIS_'//TRIM(Spc(I))
            LNAME = 'HEMCO_emissions_'//TRIM(Spc(I))
            CALL Diagn2Exp( GC, SNAME, LNAME, UNITS, 2, __RC__ )
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
      SUBROUTINE HCO_SetExtState_ESMF( am_I_Root, HcoState, ExtState, RC )
!
! !USES:
!
      USE HCO_STATE_MOD,   ONLY : Hco_State
      USE HCOX_STATE_MOD,  ONLY : Ext_State
!
! !ARGUMENTS:
!
      LOGICAL,             INTENT(IN   )   :: am_I_Root
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
      INTEGER                      :: STAT
      TYPE(ESMF_STATE), POINTER    :: IMPORT       => NULL()
      REAL,             POINTER    :: Ptr3D(:,:,:) => NULL()
      REAL,             POINTER    :: Ptr2D(:,:)   => NULL()

      ! ================================================================
      ! HCO_SetExtState_ESMF begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defines Iam and STATUS)
      __Iam__('HCO_SetExtState_ESMF (HCOI_ESMF_MOD.F90)')

      ! Assume failure until otherwise
      RC = HCO_FAIL

      ! Point to ESMF IMPORT object
      IMPORT => HcoState%IMPORT
      IF ( .NOT. ASSOCIATED(IMPORT) ) RETURN 

      ! Get pointers to fields
      CALL MAPL_GetPointer( IMPORT, Ptr3D, 'BYNCY', notFoundOK=.TRUE., __RC__ )
      IF ( ASSOCIATED(Ptr3D) ) THEN
         ExtState%BYNCY%Arr%Val => Ptr3D(:,:,HcoState%NZ:1:-1)
      ENDIF
      Ptr3D => NULL()

      ! Not needed at the moment
!      CALL MAPL_GetPointer( IMPORT, Ptr3D, 'RCCODE', notFoundOK=.TRUE., __RC__ )
!      IF ( ASSOCIATED(Ptr3D) ) THEN
!         ExtState%RCCODE%Arr%Val => Ptr3D(:,:,HcoState%NZ:1:-1)
!      ENDIF
!      Ptr3D => NULL()

!      CALL MAPL_GetPointer( IMPORT, Ptr2D, 'CNV_TOPP', notFoundOK=.TRUE., __RC__ )
!      IF ( ASSOCIATED(Ptr2D) ) THEN
!         ExtState%CNV_TOPP%Arr%Val => Ptr2D
!      ENDIF
      Ptr2D => NULL()

      ! Return success
      RC = HCO_SUCCESS 

      END SUBROUTINE HCO_SetExtState_ESMF 
!EOC
#endif
END MODULE HCOI_ESMF_MOD
