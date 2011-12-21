!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: emissions_mod
!
! !DESCRIPTION: Module EMISSIONS\_MOD is used to call the proper emissions 
!  subroutines for the various GEOS-Chem simulations.
!\\
!\\
! !INTERFACE: 
!
      MODULE EMISSIONS_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: DO_EMISSIONS
!
! !PUBLIC MEMBER DATA:
!
      !FP_ISOP (6/2009)
      PUBLIC :: ISOP_SCALING,NOx_SCALING
!
! !REVISION HISTORY:
!  11 Feb 2003 - R. Yantosca - Initial version
!  (1 ) Now references DEBUG_MSG from "error_mod.f"
!  (2 ) Now references "Kr85_mod.f" (jsw, bmy, 8/20/03)
!  (3 ) Now references "carbon_mod.f" and "dust_mod.f" (rjp, tdf, bmy, 4/2/04)
!  (4 ) Now references "seasalt_mod.f" (rjp, bmy, bec, 4/20/04)
!  (5 ) Now references "logical_mod" & "tracer_mod.f" (bmy, 7/20/04)
!  (6 ) Now references "epa_nei_mod.f" and "time_mod.f" (bmy, 11/5/04)
!  (7 ) Now references "emissions_mod.f" (bmy, 12/7/04)
!  (8 ) Now calls EMISSSULFATE if LCRYST=T.  Also read EPA/NEI emissions for 
!        the offline aerosol simulation. (bmy, 1/11/05)
!  (9 ) Remove code for the obsolete CO-OH param simulation (bmy, 6/24/05)
!  (10) Now references "co2_mod.f" (pns, bmy, 7/25/05)
!  (11) Now references "emep_mod.f" (bdf, bmy, 10/1/05)
!  (12) Now references "gfed2_biomass_mod.f" (bmy, 3/30/06)
!  (13) Now references "bravo_mod.f" (rjp, kfb, bmy, 6/26/06)
!  (14) Now references "edgar_mod.f" (avd, bmy, 7/6/06)
!  (15) Now references "streets_anthro_mod.f" (yxw, bmy, 8/18/06)
!  (16) Now references "h2_hd_mod.f" (lyj, phs, 9/18/07)
!  (17) Now calls EMISSDR for tagged CO simulation (jaf, mak, bmy, 2/14/08)
!  (18) Now references "cac_anthro_mod.f" (amv, phs, 03/11/08)
!  (19) Now references "vistas_anthro_mod.f" (amv, 12/02/08)
!  (20) Bug fixe : add specific calls for Streets for the grid 0.5x0.666.
!        (dan, ccc, 3/11/09)
!  18 Dec 2009 - Aaron van D - Added emissions for nested grids @ 0.5 x 0.666
!  26 Fev 2010 - Fabien P.   - Add scaling for isoprene and Nox emissions
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !FP_ISOP. For scaling Isoprene and NOx emissions.
      REAL*8              :: ISOP_SCALING,NOx_SCALING

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_emissions
!
! !DESCRIPTION: Subroutine DO\_EMISSIONS is the driver routine which calls 
!  the appropriate emissions subroutine for the various GEOS-CHEM simulations. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_EMISSIONS
!
! !USES:
!
! NBIOMAX has been moved to CMN_SIZE by FP (hotp 7/31/09)
!      USE BIOMASS_MOD,            ONLY : NBIOMAX
      USE BIOMASS_MOD,            ONLY : COMPUTE_BIOMASS_EMISSIONS
      USE ARCTAS_SHIP_EMISS_MOD,  ONLY : EMISS_ARCTAS_SHIP
      USE BRAVO_MOD,              ONLY : EMISS_BRAVO
      USE C2H6_MOD,               ONLY : EMISSC2H6
      USE CAC_ANTHRO_MOD,         ONLY : EMISS_CAC_ANTHRO
      USE CAC_ANTHRO_MOD,         ONLY : EMISS_CAC_ANTHRO_05x0666
      USE CARBON_MOD,             ONLY : EMISSCARBON
      USE CH3I_MOD,               ONLY : EMISSCH3I
      USE CO2_MOD,                ONLY : EMISSCO2
      USE DUST_MOD,               ONLY : EMISSDUST
      USE EDGAR_MOD,              ONLY : EMISS_EDGAR
      USE EMEP_MOD,               ONLY : EMISS_EMEP
      USE EMEP_MOD,               ONLY : EMISS_EMEP_05x0666
      USE EPA_NEI_MOD,            ONLY : EMISS_EPA_NEI
      USE ERROR_MOD,              ONLY : DEBUG_MSG
      USE GLOBAL_CH4_MOD,         ONLY : EMISSCH4
      USE H2_HD_MOD,              ONLY : EMISS_H2_HD
      USE HCN_CH3CN_MOD,          ONLY : EMISS_HCN_CH3CN
      USE ICOADS_SHIP_MOD,        ONLY : EMISS_ICOADS_SHIP
      USE LOGICAL_MOD                
      USE MERCURY_MOD,            ONLY : EMISSMERCURY
      USE NEI2005_ANTHRO_MOD,     ONLY : EMISS_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD,     ONLY : EMISS_NEI2005_ANTHRO_05x0666
      USE RnPbBe_MOD,             ONLY : EMISSRnPbBe
      USE SEASALT_MOD,            ONLY : EMISSSEASALT
      USE STREETS_ANTHRO_MOD,     ONLY : EMISS_STREETS_ANTHRO
      USE STREETS_ANTHRO_MOD,     ONLY : EMISS_STREETS_ANTHRO_05x0666
      USE SULFATE_MOD,            ONLY : EMISSSULFATE 
      USE TIME_MOD,               ONLY : GET_MONTH,       GET_YEAR
      USE TIME_MOD,               ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
      USE TRACER_MOD                 
      USE TAGGED_CO_MOD,          ONLY : EMISS_TAGGED_CO
      USE VISTAS_ANTHRO_MOD,      ONLY : EMISS_VISTAS_ANTHRO
      USE TRACERID_MOD,           ONLY : IDTNK1, IDTSF1, IDTSS1 !(win ,7/17/09)
      USE TRACERID_MOD,           ONLY : IDTECIL1, IDTECOB1, IDTDUST1 !(win ,7/17/09)
      USE TRACERID_MOD,           ONLY : IDTOCIL1, IDTOCOB1 !(win ,7/17/09)

#     include "CMN_SIZE"               ! Size parameters
#     include "CMN_O3"                 ! FSCLYR
!
! !REVISION HISTORY: 
!  (1 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (2 ) Now calls Kr85 emissions if NSRCX == 12 (jsw, bmy, 8/20/03)
!  (3 ) Now calls EMISSCARBON and EMISSDUST for carbon aerosol and dust
!        aerosol chemistry (rjp, tdf, bmy, 4/2/04)
!  (4 ) Now calls EMISSSEASALT for seasalt emissions (rjp, bec, bmy, 4/20/04)
!  (5 ) Now use inquiry functions from "tracer_mod.f".  Now references
!        "logical_mod.f" (bmy, 7/20/04)
!  (6 ) Now references ITS_A_NEW_MONTH from "time_mod.f".  Now references
!        EMISS_EPA_NEI from "epa_nei_mod.f" (bmy, 11/5/04)
!  (7 ) Now calls EMISSMERCURY from "mercury_mod.f" (eck, bmy, 12/7/04)
!  (8 ) Now calls EMISSSULFATE if LCRYST=T.  Also read EPA/NEI emissions for
!        the offline sulfate simulation.  Also call EMISS_EPA_NEI for the
!        tagged CO simulation. (cas, bmy, stu, 1/10/05).
!  (9 ) Now call EMISSSEASALT before EMISSSULFATE (bec, bmy, 4/13/05)
!  (10) Now call EMISS_HCN_CH3CN from "hcn_ch3cn_mod.f".   Also remove all 
!        references to the obsolete CO-OH param simulation. (xyp, bmy, 6/23/05)
!  (11) Now call EMISSCO2 from "co2_mod.f" (pns, bmy, 7/25/05)
!  (12) Now references EMISS_EMEP from "emep_mod.f" (bdf, bmy, 11/1/05)
!  (13) Now call GFED2_COMPUTE_BIOMASS to read 1x1 biomass emissions and
!        regrid to the model resolution once per month. (bmy, 3/30/06)
!  (14) Now references EMISS_BRAVO from "bravo_mod.f" (rjp, kfb, bmy, 6/26/06)
!  (15) Now references EMISS_EDGAR from "edgar_mod.f" (avd, bmy, 7/6/06)
!  (16) Now references EMISS_STREETS_ANTHRO from "streets_anthro_mod.f"
!        (yxw, bmy, 8/17/06)
!  (17) Now calls EMISSDR for tagged CO simulation (jaf, mak, bmy, 2/18/08)
!  (18) Now references EMISS_CAC_ANTHRO from "cac_anthro_mod.f"
!        (amv, phs, 3/11/08)
!  (19) Now references EMISS_ARCTAS_SHIP from "arctas_ship_emiss_mod.f"
!        (phs, 5/12/08)
!  (20) Now references EMISS_VISTAS_ANTHR from "vistas_anthro_mod.f". Call
!        EMEP, and Streets every month (amv, 12/2/08)
!  (21) Now references EMISS_NEI2005_ANTHRO from "nei2005_anthro_mod.f"
!        (amv, 10/19/09)
!  (22) Reference to TRACERID_MOd for IDTDUST1 for calling EMISSDUST (Win, 7/17/09)
!  18 Dec 2009 - Aaron van D - Added emissions for nested grids @ 0.5 x 0.666
!  08 Feb 2010 - NBIOMAX is now in CMN_SIZE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                     :: MONTH, YEAR
      ! BIOMASS not used (hotp 8/3/09)
      !REAL*8                      :: BIOMASS(IIPAR,JJPAR,NBIOMAX)

      !=================================================================
      ! DO_EMISSIONS begins here!
      !=================================================================

      ! Get year and month
      MONTH = GET_MONTH()

      ! check if emissions year differs from met field year
      IF ( FSCALYR < 0 ) THEN
         YEAR = GET_YEAR()
      ELSE
         YEAR = FSCALYR
      ENDIF


      ! Get biomass burning emissions for use below
      IF ( LBIOMASS ) THEN
!#### HARDWIRE SWITCH : use YEAR instead of GET_YEAR() to use the same
!#### base year as anthropogenic emissions
!         CALL COMPUTE_BIOMASS_EMISSIONS( YEAR, MONTH )
         CALL COMPUTE_BIOMASS_EMISSIONS( GET_YEAR(), MONTH )
      ENDIF
         
      ! Test by simulation type
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !--------------------
         ! NOx-Ox-HC-aerosol
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666      !(dan)
#else
            CALL EMISS_STREETS_ANTHRO
#endif
         ENDIF

         ! Read EDGAR emissions once per month to get, at least 
         ! the NOx diurnal scale factors, and the EDGAR emissions 
         ! if necessary (amv, phs, 3/11/08)
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EDGAR( YEAR, MONTH )
         ENDIF

         ! Read EPA/NEI99 (USA) emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) CALL EMISS_EPA_NEI

         ! Read VISTAS (USA) emissions once per month
         IF ( LVISTAS .and. ITS_A_NEW_MONTH() ) 
     &        CALL EMISS_VISTAS_ANTHRO

         ! Read BRAVO (Mexico) emissions once per year
         IF ( LBRAVO .and. ITS_A_NEW_YEAR()  ) CALL EMISS_BRAVO

         ! Read EMEP (Europe) emissions once per year
         IF ( LEMEP  .and. ITS_A_NEW_MONTH()  ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_EMEP_05x0666
#else
            CALL EMISS_EMEP
#endif
         ENDIF

         ! Read CAC (Canada) emissions
         IF ( LCAC .and. ITS_A_NEW_YEAR() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_CAC_ANTHRO_05x0666
#else
            CALL EMISS_CAC_ANTHRO
#endif
         ENDIF

         ! Read NEI2005 (USA) emissions
         IF ( LNEI05 .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_NEI2005_ANTHRO_05x0666
#else
            CALL EMISS_NEI2005_ANTHRO
#endif
         ENDIF

         ! Read SO2 ARCTAS emissions
         IF ( LARCSHIP .AND. ITS_A_NEW_YEAR() )
     $        CALL EMISS_ARCTAS_SHIP( YEAR )

         ! Read ICOADS ship emissions once per month (cklee, 7/09/09)
         IF ( LICOADSSHIP .and. ITS_A_NEW_MONTH() ) 
     &        CALL EMISS_ICOADS_SHIP

         ! NOx-Ox-HC (w/ or w/o aerosols)
         CALL EMISSDR

         ! Emissions for various aerosol types
!prior to TOMAS addition (win, 7/17/09)
c$$$         IF ( LSSALT            ) CALL EMISSSEASALT
c$$$         IF ( LSULF .or. LCRYST ) CALL EMISSSULFATE
c$$$         IF ( LCARB             ) CALL EMISSCARBON
c$$$         IF ( LDUST             ) CALL EMISSDUST
         IF ( LSSALT .or. IDTSS1 > 0  ) CALL EMISSSEASALT
         IF ( LSULF .or. LCRYST .or. IDTSF1 > 0 ) CALL EMISSSULFATE
         IF ( LCARB .or. ( IDTECIL1 > 0 .and.   IDTECOB1 > 0 .and. 
     &        IDTOCIL1 > 0 .and. IDTOCOB1 > 0 ) ) CALL EMISSCARBON
         IF ( LDUST .or. IDTDUST1 > 0  ) CALL EMISSDUST

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN
         
         !--------------------
         ! Offline aerosol
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666      !(dan)
#else
            CALL EMISS_STREETS_ANTHRO
#endif
         ENDIF

         ! Read CAC emissions
         IF ( LCAC .and. ITS_A_NEW_YEAR() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_CAC_ANTHRO_05x0666
#else
            CALL EMISS_CAC_ANTHRO
#endif
         ENDIF

         ! Read EDGAR emissions once per month
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EDGAR( YEAR, MONTH )
         ENDIF

         ! Read EPA/NEI99 emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) CALL EMISS_EPA_NEI

         ! Read NEI2005 emissions once per month
         IF ( LNEI05 .and. ITS_A_NEW_MONTH() ) THEN
#if    defined(GRID05x0666)
            CALL EMISS_NEI2005_ANTHRO_05x0666
#else
            CALL EMISS_NEI2005_ANTHRO
#endif
         ENDIF

         ! Read BRAVO (Mexico) emissions once per year
         IF ( LBRAVO .and. ITS_A_NEW_YEAR()  ) CALL EMISS_BRAVO

         ! Read EMEP (Europe) emissions once per year
         IF ( LEMEP  .and. ITS_A_NEW_YEAR()  ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_EMEP_05x0666
#else       
            CALL EMISS_EMEP
#endif
         ENDIF

         ! Read SO2 ARCTAS emissions
         IF ( LARCSHIP .AND. ITS_A_NEW_YEAR() )
     $        CALL EMISS_ARCTAS_SHIP( YEAR )

         ! Read ICOADS ship emissions once per month !(cklee, 7/09/09)
         IF ( LICOADSSHIP .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_ICOADS_SHIP
         ENDIF

         ! Emissions for various aerosol types
         IF ( LSSALT            ) CALL EMISSSEASALT
         IF ( LSULF .or. LCRYST ) CALL EMISSSULFATE
         IF ( LCARB             ) CALL EMISSCARBON
         IF ( LDUST             ) CALL EMISSDUST

      ELSE IF ( ITS_A_RnPbBe_SIM() ) THEN
         
         !--------------------
         ! Rn-Pb-Be
         !--------------------
         CALL EMISSRnPbBe

      ELSE IF ( ITS_A_CH3I_SIM() ) THEN

         !--------------------
         ! CH3I
         !--------------------

         ! Emit CH3I
         CALL EMISSCH3I

      ELSE IF ( ITS_A_HCN_SIM() ) THEN

         !--------------------
         ! HCN - CH3CN
         !--------------------
         CALL EMISS_HCN_CH3CN( N_TRACERS, STT )

      ELSE IF ( ITS_A_TAGCO_SIM() ) THEN

         !--------------------
         ! Tagged CO
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         ! Bug fix: call every month now (pdk, phs, 3/17/09)
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666      !(dan)
#else
            CALL EMISS_STREETS_ANTHRO
#endif
         ENDIF

         ! Read CAC emissions
         IF ( LCAC .and. ITS_A_NEW_YEAR() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_CAC_ANTHRO_05x0666
#else
            CALL EMISS_CAC_ANTHRO
#endif
         ENDIF

         ! Read EDGAR emissions once per month
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EDGAR( YEAR, MONTH )
         ENDIF

         ! Read EPA (USA) emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) CALL EMISS_EPA_NEI

         ! Read NEI2005 (USA) emissions once per year
         IF ( LNEI05 .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_NEI2005_ANTHRO_05x0666
#else
            CALL EMISS_NEI2005_ANTHRO
#endif
         ENDIF

         ! Read BRAVO (Mexico) emissions once per year
         IF ( LBRAVO .and. ITS_A_NEW_YEAR()  ) CALL EMISS_BRAVO

         ! Read EPA (Europe) emissions once per year
         IF ( LEMEP  .and. ITS_A_NEW_YEAR()  ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_EMEP_05x0666
#else
            CALL EMISS_EMEP
#endif
         ENDIF

         ! Read ICOADS ship emissions once per month (cklee, 7/09/09)
         IF ( LICOADSSHIP .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_ICOADS_SHIP
         ENDIF

         ! Now call EMISSDR for Tagged CO fossil fuel emissions, 
         ! so that we get the same emissions for Tagged CO as 
         ! we do for the full-chemistry (jaf, mak, bmy, 2/14/08)
         CALL EMISSDR 
         
         ! Emit tagged CO
         CALL EMISS_TAGGED_CO

      ELSE IF ( ITS_A_C2H6_SIM() ) THEN

         !--------------------
         ! C2H6
         !--------------------

         ! Emit C2H6
         CALL EMISSC2H6

      ELSE IF ( ITS_A_CH4_SIM() ) THEN

         !--------------------
         ! CH4
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         ! Bug fix: call every month now (phs, 3/17/09)
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666      !(dan)
#else
            CALL EMISS_STREETS_ANTHRO
#endif
         ENDIF

         ! Emit CH4
         CALL EMISSCH4

      ELSE IF ( ITS_A_MERCURY_SIM() ) THEN

         !--------------------
         ! Mercury
         !--------------------
         CALL EMISSMERCURY

      ELSE IF ( ITS_A_CO2_SIM() ) THEN

         !--------------------
         ! CO2
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         ! Bug fix: call every month now (phs, 3/17/09)         
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666      !(dan)
#else
            CALL EMISS_STREETS_ANTHRO
#endif
         ENDIF

         ! Read CO2 ARCTAS SHIP emissions
         IF ( LARCSHIP .and. ITS_A_NEW_YEAR() )
     $        CALL EMISS_ARCTAS_SHIP( YEAR )

         ! Emit CO2
         CALL EMISSCO2

      ELSE IF ( ITS_A_H2HD_SIM() ) THEN

         !--------------------
         ! Offline H2/HD 
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         ! Bug fix: call every month now (phs, 3/17/09)
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666      !(dan)
#else
            CALL EMISS_STREETS_ANTHRO
#endif
         ENDIF

         ! Read EDGAR emissions once per month
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EDGAR( YEAR, MONTH )
         ENDIF

         ! Read CAC (Canada) emissions
         IF ( LCAC .and. ITS_A_NEW_YEAR() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_CAC_ANTHRO_05x0666
#else
            CALL EMISS_CAC_ANTHRO
#endif
         ENDIF

         ! Read EPA (USA) emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) CALL EMISS_EPA_NEI

         ! Read NEI2005 (USA) emissions
         IF ( LNEI05 .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_NEI2005_ANTHRO_05x0666
#else
            CALL EMISS_NEI2005_ANTHRO
#endif
         ENDIF

         ! Read BRAVO (Mexico) emissions once per year
         IF ( LBRAVO .and. ITS_A_NEW_YEAR()  ) CALL EMISS_BRAVO

         ! Read EPA (Europe) emissions once per year
         IF ( LEMEP  .and. ITS_A_NEW_YEAR()  ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_EMEP_05x0666
#else
            CALL EMISS_EMEP
#endif
         ENDIF

         ! Read ICOADS ship emissions once per month !(cklee, 7/09/09)
         IF ( LICOADSSHIP .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_ICOADS_SHIP
         ENDIF

         ! Emit H2/HD
         CALL EMISS_H2_HD

      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG ( '### DO_EMISSIONS: a EMISSIONS' )

      END SUBROUTINE DO_EMISSIONS
!EOC
      END MODULE EMISSIONS_MOD
