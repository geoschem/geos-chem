!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_test_utils
!
! !DESCRIPTION: Module GC\_TEST\_UTILS contains debugging code.
!\\
!\\
! !INTERFACE:
!
MODULE GC_TEST_UTILS

  PUBLIC
!
! !REMARKS:

! !REVISION HISTORY: 
!  18 Oct 2012 - M. Long     - Initial version
!  18 Oct 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dump_gc_config
!
! !DESCRIPTION: Prints out information about GEOS-Chem options.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DUMP_GC_CONFIG( am_I_Root )
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

  END SUBROUTINE DUMP_GC_CONFIG
!EOC
END MODULE GC_TEST_UTILS
