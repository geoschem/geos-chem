MODULE GC_TEST_UTILS

  PUBLIC

CONTAINS

  SUBROUTINE DUMP_GC_CONFIG
    USE LOGICAL_MOD
    USE GC_Value_Mod, ONLY : MASTERPROC
    
    IF (MASTERPROC) THEN
       !%%%% Aerosols %%%%
       WRITE(*,*) 'LCARB: ',    LCARB          ! Use carbon aerosol tracers?
       WRITE(*,*) 'LCRYST: ',   LCRYST         ! Use Crystalline aerosols? (not implemented)
       WRITE(*,*) 'LDEAD: ',    LDEAD          ! Use the DEAD/Zender dust algorithm?
       WRITE(*,*) 'LDUST: ',    LDUST          ! Use dust aerosol tracers?
       WRITE(*,*) 'LSULF: ',    LSULF          ! Use sulfate aerosol tracers?
       WRITE(*,*) 'LSOA: ',     LSOA           ! Use secondary organic aerosol tracers?
       WRITE(*,*) 'LSSALT: ',   LSSALT         ! Use sea-salt aerosol tracers?
       WRITE(*,*) 'LDICARB: ',  LDICARB        ! Use dicarbonyl chemistry mechanism?
                                   
       !%%%% Chemistry %%%%stry %%%
       WRITE(*,*) 'LCHEM: ',    LCHEM          ! Use chemistry?
       WRITE(*,*) 'LKPP: ',     LKPP           ! Use KPP solver instead of SMVGEAR?
                                   
       WRITE(*,*) 'LDRYD: ',    LDRYD          ! Use dry deposition?
                                   
       !%%%% Variable Tropoble Troppause %%%%
       WRITE(*,*) 'LVARTROP: ', LVARTROP       ! Use dynamic tropopause option?

    ENDIF

  END SUBROUTINE DUMP_GC_CONFIG

END MODULE GC_TEST_UTILS
