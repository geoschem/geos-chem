#if defined (ESMF_)
#include "MAPL_Generic.h"
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
  USE ESMF
  USE MAPL_Mod

  USE HCO_ERROR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!      
  PUBLIC :: HCO_SetServices 
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
! The field names provided in ExtData.rc must match the names in the HEMCO 
! configuration file! Also, all time settings (average and update interval) 
! and data units need to be properly specified in ExtData.rc.
! TODO: For now, ExtData.rc and HEMCO configuration file have to be 
! synchronized manually. The pyHEMCO interface will automate this process! 
!\\
!\\
! Since this routine is called from outside of the HEMCO environment, use the 
! MAPL specific error codes!
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
      INTEGER                    :: FLAG
      TYPE(ListCont), POINTER    :: CurrCont => NULL()

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

      ! Return success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCO_SetServices
!EOC
END MODULE HCOI_ESMF_MOD
#endif
