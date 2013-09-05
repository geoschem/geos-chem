!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_test_utils
!
! !DESCRIPTION: Module GIGC\_Test\_Utils contains debugging code for the
!  ESMF interface to the Grid-Independent GEOS-Chem (aka "GIGC").
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_Test_Utils
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Dump_Config
!
! !REMARKS:
!  Additional debugging routines can be added here as necessary
!
! !REVISION HISTORY: 
!  18 Oct 2012 - M. Long     - Initial version
!  18 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_test_utils.F90
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
! !IROUTINE: gigc_dump_config
!
! !DESCRIPTION: Prints out information about GEOS-Chem options.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Dump_Config( am_I_Root )
!
! !USES:
!
    USE LOGICAL_MOD
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN) :: am_I_root
! 
! !REVISION HISTORY: 
!  18 Oct 2012 - M. Long     - Initial version
!  18 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Dump_Config
!EOP
!------------------------------------------------------------------------------
!BOC
   
    IF ( am_I_Root ) THEN 

       !%%%% Aerosols %%%%
       WRITE(*,*) 'LCARB: ',    LCARB         ! Use carbon aerosol tracers?
       WRITE(*,*) 'LCRYST: ',   LCRYST        ! Use Crystalline aerosols?
       WRITE(*,*) 'LDEAD: ',    LDEAD         ! Use the DEAD/Zender dust?
       WRITE(*,*) 'LDUST: ',    LDUST         ! Use dust aerosol tracers?
       WRITE(*,*) 'LSULF: ',    LSULF         ! Use sulfate aerosol tracers?
       WRITE(*,*) 'LSOA: ',     LSOA          ! Use SOA tracers?
       WRITE(*,*) 'LSSALT: ',   LSSALT        ! Use sea-salt aerosol tracers?
       WRITE(*,*) 'LDICARB: ',  LDICARB       ! Use dicarbonyl chemistry
                                   
       !%%%% Chemistry %%%%stry %%%
       WRITE(*,*) 'LCHEM: ',    LCHEM         ! Use chemistry?
       WRITE(*,*) 'LKPP: ',     LKPP          ! Use KPP solver>?
                                   
       WRITE(*,*) 'LDRYD: ',    LDRYD         ! Use dry deposition?
                                   
       !%%%% Variable Tropoble Troppause %%%%
       WRITE(*,*) 'LVARTROP: ', LVARTROP      ! Use dynamic tropopause option?

    ENDIF

  END SUBROUTINE GIGC_Dump_Config
!EOC
END MODULE GIGC_Test_Utils
