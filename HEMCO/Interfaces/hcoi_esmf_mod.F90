!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  PUBLIC :: HCO_SetServices 
  PUBLIC :: HCOI_ESMF_DIAGNCREATE
!
! !REVISION HISTORY:
!  10 Oct 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
! The corresponding HEMCO diagnostics can be created using subroutine
! HCOI\_ESMF\_DIAGNCREATE.
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
      INTEGER                    :: I, FLAG, nSpc
      CHARACTER(LEN=63), POINTER :: Spc(:)   => NULL() 
      TYPE(ListCont),    POINTER :: CurrCont => NULL()

      ! ================================================================
      ! HCO_SetServices begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defined Iam and STATUS)
      __Iam__('HCO_SetServices (HCOI_ESMF_MOD.F90)')

      ! ---------------------------------------------------------------------
      ! Read file into buffer
      ! ---------------------------------------------------------------------

      CALL Config_ReadFile( am_I_Root, TRIM(ConfigFile), STATUS )
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
      ! Prepare a diagnostics export for every potential HEMCO species
      ! ---------------------------------------------------------------------
      nSpc = Config_GetnSpecies( )
      CALL Config_GetSpecNames( Spc, nSpc, RC )
      ASSERT_(RC == HCO_SUCCESS) 

      ! Loop over all species
      DO I = 1, nSpc

         ! Add to export state
         call MAPL_AddExportSpec(GC,                               &
            SHORT_NAME         = 'EMIS_'//TRIM(Spc(I)),            &
            LONG_NAME          = 'HEMCO_emissions_'//TRIM(Spc(I)), &
            UNITS              = 'kg m-2 s-1',                     &
            DIMS               = MAPL_DimsHorzOnly,                &
            VLOCATION          = MAPL_VLocationNone,               &
            RC                 = STATUS                             )
         IF ( STATUS /= ESMF_SUCCESS ) THEN
            WRITE(*,*) 'Cannot add to export: ',TRIM(Spc(I))//'_EMIS'
            ASSERT_(.FALSE.)
         ENDIF
      ENDDO

      ! Cleanup
      IF ( ASSOCIATED(Spc) ) DEALLOCATE(Spc)

      ! Return success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCO_SetServices
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCOI_ESMF_DIAGNCREATE
!
! !DESCRIPTION: Subroutine HCOI\_ESMF\_DIAGNCREATE creates an emissions 
! diagnostics for every HEMCO species defined in the HEMCO state object.
!\\
!\\
! These diagnostics are expected to be written into the export states 
! defined when setting the HEMCO services (HCO\_SetServices). They represent
! the total species emissions. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOI_ESMF_DIAGNCREATE( am_I_Root, HcoState, RC )
!
! !USES:
!
      USE HCO_DIAGN_MOD, ONLY : DIAGN_CREATE
      USE HCO_STATE_MOD, ONLY : HCO_STATE
!
! !ARGUMENTS:
!
      LOGICAL,             INTENT(IN   )   :: am_I_Root
      TYPE(HCO_STATE),     POINTER         :: HcoState
      INTEGER,             INTENT(  OUT)   :: RC
!
! !REVISION HISTORY:
!  11 Nov 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                    :: I, N, HcoID
      CHARACTER(LEN=31)          :: Spc

      ! ================================================================
      ! HCOI_ESMF_DIAGNCREATE begins here
      ! ================================================================

      DO I = 1, HcoState%nSpc

         ! Get HEMCO species information
         Spc   = HcoState%Spc(I)%SpcName
         HcoID = HcoState%Spc(I)%HcoID
         IF ( HcoID <= 0 ) CYCLE
         
         ! Define diagnostics
         CALL Diagn_Create( am_I_Root,                     &
                            HcoState,                      &
                            cName    = 'EMIS_'//TRIM(Spc), &
                            ExtNr    = -1,                 &
                            Cat      = -1,                 &
                            Hier     = -1,                 &
                            HcoID    = HcoID,              &
                            SpaceDim = 2,                  &
                            LevIDx   = -1,                 &
                            OutUnit  = 'kg/m2/s',          &
                            WriteFreq = 'Manual',          &
                            AutoFill  = 1,                 &
                            cID       = N,                 & 
                            RC        = RC                  ) 
        IF ( RC /= HCO_SUCCESS ) RETURN

      ENDDO

      ! Return success
      RC = HCO_SUCCESS

      END SUBROUTINE HCOI_ESMF_DIAGNCREATE
!EOC
#endif
END MODULE HCOI_ESMF_MOD
