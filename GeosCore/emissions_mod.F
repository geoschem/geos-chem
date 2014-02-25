!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  01 Feb 2011 - C Friedman  - Added POP emissions
!  07 Feb 2011 - R. Yantosca - Now use EPA/NEI99 biofuel emissions when
!                              EPA/NEI05 anthro emissions are selected.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE DO_EMISSIONS( am_I_Root, Input_Opt, 
     &                         State_Met, State_Chm, RC )
!
! !USES:
!
      USE ARCTAS_SHIP_EMISS_MOD,  ONLY : EMISS_ARCTAS_SHIP
      USE BIOMASS_MOD,            ONLY : COMPUTE_BIOMASS_EMISSIONS
      USE BRAVO_MOD,              ONLY : EMISS_BRAVO
      USE C2H6_MOD,               ONLY : EMISSC2H6
      USE CAC_ANTHRO_MOD,         ONLY : EMISS_CAC_ANTHRO
      USE CAC_ANTHRO_MOD,         ONLY : EMISS_CAC_ANTHRO_05x0666
      USE CARBON_MOD,             ONLY : EMISSCARBON
      USE CH3I_MOD,               ONLY : EMISSCH3I
      USE CMN_O3_MOD
      USE CMN_SIZE_MOD
      USE CO2_MOD,                ONLY : EMISSCO2
      USE DUST_MOD,               ONLY : EMISSDUST
      USE EDGAR_MOD,              ONLY : EMISS_EDGAR
      USE EMEP_MOD,               ONLY : EMISS_EMEP
      USE EMEP_MOD,               ONLY : EMISS_EMEP_05x0666
      USE EPA_NEI_MOD,            ONLY : EMISS_EPA_NEI
      USE ERROR_MOD,              ONLY : DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,     ONLY : OptInput
      USE GIGC_State_Chm_Mod,     ONLY : ChmState
      USE GIGC_State_Met_Mod,     ONLY : MetState
      USE GLOBAL_CH4_MOD,         ONLY : EMISSCH4
      USE H2_HD_MOD,              ONLY : EMISS_H2_HD
      USE HCN_CH3CN_MOD,          ONLY : EMISS_HCN_CH3CN
      USE ICOADS_SHIP_MOD,        ONLY : EMISS_ICOADS_SHIP
      USE MERCURY_MOD,            ONLY : EMISSMERCURY
      USE NEI2005_ANTHRO_MOD,     ONLY : EMISS_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD,     ONLY : EMISS_NEI2005_ANTHRO_05x0666
      USE PARANOX_MOD,            ONLY : READ_PARANOX_LUT
      USE POPS_MOD,               ONLY : EMISSPOPS  !(clf, 2/1/2011)
      USE RCP_MOD,                ONLY : LOAD_RCP_EMISSIONS !cdh
      USE RETRO_MOD,              ONLY : EMISS_RETRO
      USE RnPbBe_MOD,             ONLY : EMISSRnPbBe
      USE SEASALT_MOD,            ONLY : EMISSSEASALT
      USE SSA_BROMINE_MOD,        ONLY : EMIT_Br2
      USE STREETS_ANTHRO_MOD,     ONLY : EMISS_STREETS_ANTHRO
      USE STREETS_ANTHRO_MOD,     ONLY : EMISS_STREETS_ANTHRO_05x0666
      USE SULFATE_MOD,            ONLY : EMISSSULFATE 
      USE TAGGED_CO_MOD,          ONLY : EMISS_TAGGED_CO
      USE TIME_MOD,               ONLY : GET_MONTH
      USE TIME_MOD,               ONLY : GET_YEAR
      USE TIME_MOD,               ONLY : ITS_A_NEW_MONTH
      USE TIME_MOD,               ONLY : ITS_A_NEW_YEAR
      USE TRACERID_MOD,           ONLY : IDTSO2, IDTC2H6
      USE VISTAS_ANTHRO_MOD,      ONLY : EMISS_VISTAS_ANTHRO

#if   defined( TOMAS )
      USE TRACERID_MOD,           ONLY : IDTSS1           ! (win, 1/25/10)
      USE TRACERID_MOD,           ONLY : IDTSF1           ! (win, 1/25/10)
      USE TRACERID_MOD,           ONLY : IDTECIL1         ! (win, 1/25/10)
      USE TRACERID_MOD,           ONLY : IDTOCIL1         ! (win, 1/25/10)
      USE TRACERID_MOD,           ONLY : IDTECOB1         ! (win, 1/25/10)
      USE TRACERID_MOD,           ONLY : IDTOCOB1         ! (win, 1/25/10)
      USE TRACERID_MOD,           ONLY : IDTDUST1         ! (win, 1/25/10)
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
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
!  01 Feb 2011 - CFriedman   - Added emissions for POPs
!  07 Feb 2011 - R. Yantosca - Use NEI99 biofuels when useing NEI05 anthro
!  17 Aug 2011 - R. Yantosca - Added call to RETRO anthro emissions
!  07 Feb 2012 - M. Payer    - Added call to read PARANOX look up tables
!  22 Mar 2012 - M. Payer    - Added call to EMISSC2H6
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  19 Oct 2012 - R. Yantosca - Rename CHEM_STATE argument to State_Chm
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Now pass all met fields via State_Met object
!  26 Mar 2013 - R. Yantosca - Now use logical fields from Input_Opt
!  27 Mar 2013 - R. Yantosca - Now remove DEVEL tags; pass objects to routines
!  22 Jul 2013 - M. Sulprizio- Now copy LRCP, LRCPSHIP, LRCPAIR from Input_Opt
!  15 Oct 2013 - R. Yantosca - Remove calls to PARANOX for simulations other
!                              than full-chemistry (they don't carry NO, O3)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: MONTH
      INTEGER :: YEAR
      REAL*8  :: SSA_Br2(IIPAR,JJPAR)

      ! For fields from Input_Opt
      LOGICAL :: ITS_A_FULLCHEM_SIM
      LOGICAL :: ITS_AN_AEROSOL_SIM
      LOGICAL :: ITS_A_RnPbBe_SIM
      LOGICAL :: ITS_A_CH3I_SIM
      LOGICAL :: ITS_A_HCN_SIM
      LOGICAL :: ITS_A_TAGCO_SIM
      LOGICAL :: ITS_A_C2H6_SIM
      LOGICAL :: ITS_A_CH4_SIM
      LOGICAL :: ITS_A_MERCURY_SIM
      LOGICAL :: ITS_A_POPS_SIM
      LOGICAL :: ITS_A_CO2_SIM
      LOGICAL :: ITS_A_H2HD_SIM
      LOGICAL :: LBIOMASS
      LOGICAL :: LSTREETS
      LOGICAL :: LRETRO
      LOGICAL :: LNEI99
      LOGICAL :: LVISTAS
      LOGICAL :: LBRAVO
      LOGICAL :: LEMEP
      LOGICAL :: LCAC
      LOGICAL :: LNEI05
      LOGICAL :: LARCSHIP
      LOGICAL :: LEDGARSHIP
      LOGICAL :: LEMEPSHIP
      LOGICAL :: LICOADSSHIP
      LOGICAL :: LSSALT
      LOGICAL :: LCARB
      LOGICAL :: LSULF
      LOGICAL :: LCRYST
      LOGICAL :: LDUST
      LOGICAL :: LSSABr2
      LOGICAL :: LPRT
      LOGICAL :: LRCP
      LOGICAL :: LRCPSHIP
      LOGICAL :: LRCPAIR
      INTEGER :: N_TRACERS

      !=================================================================
      ! DO_EMISSIONS begins here!
      !=================================================================

      ! Assume success
      RC                 = GIGC_SUCCESS

      ! Copy values from Input_Opt
      ITS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
      ITS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM
      ITS_A_RnPbBe_SIM   = Input_Opt%ITS_A_RnPbBe_SIM
      ITS_A_CH3I_SIM     = Input_Opt%ITS_A_CH3I_SIM
      ITS_A_HCN_SIM      = Input_Opt%ITS_A_HCN_SIM
      ITS_A_TAGCO_SIM    = Input_Opt%ITS_A_TAGCO_SIM
      ITS_A_C2H6_SIM     = Input_Opt%ITS_A_C2H6_SIM
      ITS_A_CH4_SIM      = Input_Opt%ITS_A_CH4_SIM
      ITS_A_MERCURY_SIM  = Input_Opt%ITS_A_MERCURY_SIM
      ITS_A_POPS_SIM     = Input_Opt%ITS_A_POPS_SIM
      ITS_A_CO2_SIM      = Input_Opt%ITS_A_CO2_SIM
      ITS_A_H2HD_SIM     = Input_Opt%ITS_A_H2HD_SIM
      LBIOMASS           = Input_Opt%LBIOMASS
      LSTREETS           = Input_Opt%LSTREETS
      LRETRO             = Input_Opt%LRETRO
      LNEI99             = Input_Opt%LNEI99
      LVISTAS            = Input_Opt%LVISTAS
      LBRAVO             = Input_Opt%LBRAVO
      LEMEP              = Input_Opt%LEMEP
      LCAC               = Input_Opt%LCAC
      LNEI05             = Input_Opt%LNEI05
      LARCSHIP           = Input_Opt%LARCSHIP
      LEDGARSHIP         = Input_Opt%LEDGARSHIP
      LEMEPSHIP          = Input_Opt%LEMEPSHIP
      LICOADSSHIP        = Input_Opt%LICOADSSHIP
      LSSALT             = Input_Opt%LSSALT
      LCARB              = Input_Opt%LCARB
      LSULF              = Input_Opt%LSULF
      LCRYST             = Input_Opt%LCRYST
      LDUST              = Input_Opt%LDUST
      LSSABr2            = Input_Opt%LSSABr2
      LPRT               = Input_Opt%LPRT
      LRCP               = Input_Opt%LRCP
      LRCPSHIP           = Input_Opt%LRCPSHIP
      LRCPAIR            = Input_Opt%LRCPAIR
      N_TRACERS          = Input_Opt%N_TRACERS

#if defined( TOMAS )
      ! Extra error checks for TOMAS simulation
      LSSALT             = ( LSSALT .and. ( IDTSS1   > 0 ) )
      LSULF              = ( LSULF  .and. ( IDTSF1   > 0 ) )
      LCARB              = ( LCARB  .and. ( IDTECIL1 > 0 ) 
     &                              .and. ( IDTECOB1 > 0 ) 
     &                              .and. ( IDTOCIL1 > 0 ) 
     &                              .and. ( IDTOCOB1 > 0 ) )
#endif

      ! Get year and month
      MONTH              = GET_MONTH()

      ! check if emissions year differs from met field year
      IF ( FSCALYR < 0 ) THEN
         YEAR            = GET_YEAR()
      ELSE
         YEAR            = FSCALYR
      ENDIF

      ! Get biomass burning emissions for use below
      IF ( LBIOMASS ) THEN
!----------------------------------------------------------------------------
!#### HARDWIRE SWITCH : use YEAR instead of GET_YEAR() to use the same
!#### base year as anthropogenic emissions
!         CALL COMPUTE_BIOMASS_EMISSIONS( am_I_Root, YEAR,      
!                                         MONTH,     Input_Opt, 
!                                         State_Chm, RC         )  
!----------------------------------------------------------------------------
         CALL COMPUTE_BIOMASS_EMISSIONS
     &      ( am_I_Root, GET_YEAR(), MONTH,     
     &        Input_Opt, State_Chm, RC      )  
      ENDIF
         
      ! Test by simulation type
      IF ( ITS_A_FULLCHEM_SIM ) THEN

         !--------------------
         ! NOx-Ox-HC-aerosol
         ! (aka "fullchem")
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
#else
            CALL EMISS_STREETS_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
#endif
         ENDIF

         ! Read EDGAR emissions once per month to get, at least 
         ! the NOx diurnal scale factors, and the EDGAR emissions 
         ! if necessary (amv, phs, 3/11/08)
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EDGAR
     &         ( am_I_Root, YEAR, MONTH, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read RETRO emissions once per month (wfr, 3/8/11)
         IF ( LRETRO .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_RETRO( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read RCP emissions once per month (cdh, 10/14/11)
         ! We read all emissions (for land, ship, air) even if only
         ! one is needed
         IF ( (LRCP .OR. LRCPSHIP .OR. LRCPAIR) .and. 
     &        ITS_A_NEW_MONTH() ) CALL LOAD_RCP_EMISSIONS( Input_Opt )

         ! Read EPA/NEI99 (USA) emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EPA_NEI( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read VISTAS (USA) emissions once per month
         IF ( LVISTAS .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_VISTAS_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read BRAVO (Mexico) emissions once per year
         IF ( LBRAVO .and. ITS_A_NEW_YEAR() ) THEN
            CALL EMISS_BRAVO( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read EMEP (Europe) emissions once per year
         IF ( LEMEP  .and. ITS_A_NEW_MONTH()  ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_EMEP_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! nested
#else
            CALL EMISS_EMEP
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read CAC (Canada) emissions
         IF ( LCAC .and. ITS_A_NEW_MONTH() ) THEN
#if   defined( GRID05x0666 )
            CALL EMISS_CAC_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! nested
#else
            CALL EMISS_CAC_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read NEI2005 (USA) emissions
         IF ( LNEI05 .and. ITS_A_NEW_MONTH() ) THEN
#if   defined( GRID05x0666 )
            CALL EMISS_EPA_NEI
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI99 biofuel, nested
            CALL EMISS_NEI2005_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI05 anthro,  nested
#else
            CALL EMISS_EPA_NEI
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI99 biofuel, global
            CALL EMISS_NEI2005_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI05 anthro,  global
#endif
         ENDIF

         ! Read SO2 ARCTAS emissions
         IF ( LARCSHIP .AND. ITS_A_NEW_YEAR() ) THEN
            CALL EMISS_ARCTAS_SHIP
     &         ( am_I_Root, YEAR, Input_Opt, State_Chm, RC )
         ENDIF
            
         ! Read ICOADS ship emissions once per month (cklee, 7/09/09)
         IF ( LICOADSSHIP .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_ICOADS_SHIP
     &         ( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read look up tables for PARANOX ship plume model (mpayer, 2/7/12)
         IF ( ( LICOADSSHIP   .or.  LEDGARSHIP .or. 
     &          LEMEPSHIP  ) .and. ITS_A_NEW_MONTH() )  THEN
            CALL READ_PARANOX_LUT
         ENDIF

         ! C2H6 emissions using RETRO are too low, use Yaping Xiao's
         ! offline C2H6 emissions instead (mpayer, 3/22/12)
         CALL EMISSC2H6
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

         ! NOx-Ox-HC (w/ or w/o aerosols)
         CALL EMISSDR
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

         ! Seasalt aerosol emission
         IF ( LSSALT ) THEN 
            CALL EMISSSEASALT
     &         ( am_I_Root, Input_Opt, State_Met, 
     &           State_Chm, RC,        ssa_br2    )
         ENDIF

         ! Sulfate aerosol emission
         IF ( LSULF .or. LCRYST ) THEN
            CALL EMISSSULFATE
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF

         ! Carbonaceous/SOA aerosol emission
         IF ( LCARB ) THEN
            CALL EMISSCARBON 
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF

         ! Mineral dust aerosol emission
         IF ( LDUST ) THEN
            CALL EMISSDUST
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF

         ! jpp, 3/2/10:
         ! Now distribute the Br2 emitted from SSA, through
         ! the full mixed layer.
         IF ( LSSABr2 ) THEN
            CALL EMIT_BR2( SSA_Br2 )
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM ) THEN
         
         !--------------------
         ! Offline aerosol
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! nested
#else
            CALL EMISS_STREETS_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read CAC emissions
         IF ( LCAC .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_CAC_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! nested
#else
            CALL EMISS_CAC_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read EDGAR emissions once per month
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EDGAR
     &         ( am_I_Root, YEAR,      MONTH,
     &           Input_Opt, State_Chm, RC     )
         ENDIF

         ! Read EPA/NEI99 emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EPA_NEI( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read NEI2005 emissions once per month
         IF ( LNEI05 .and. ITS_A_NEW_MONTH() ) THEN
#if   defined( GRID05x0666 )
            CALL EMISS_EPA_NEI
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI99 biofuel, nested
            CALL EMISS_NEI2005_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI05 anthro,  nested
#else
            CALL EMISS_EPA_NEI 
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI99 biofuel, global
            CALL EMISS_NEI2005_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI05 anthro,  global
#endif
         ENDIF

         ! Read BRAVO (Mexico) emissions once per year
         IF ( LBRAVO .and. ITS_A_NEW_YEAR() ) THEN
            CALL EMISS_BRAVO( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read EMEP (Europe) emissions once per year
         IF ( LEMEP  .and. ITS_A_NEW_YEAR()  ) THEN
#if   defined( GRID05x0666 )
            CALL EMISS_EMEP_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! nested
#else       
            CALL EMISS_EMEP
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read SO2 ARCTAS emissions
         IF ( LARCSHIP .AND. ITS_A_NEW_YEAR() ) THEN
            CALL EMISS_ARCTAS_SHIP
     &         ( am_I_Root, YEAR, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read ICOADS ship emissions once per month !(cklee, 7/09/09)
         IF ( LICOADSSHIP .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_ICOADS_SHIP
     &         ( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

!-----------------------------------------------------------------------------
! Prior to 10/15/13:
! PARANOX is not needed for offline aerosol simulations, which do not
! carry NO or O3.  Remove this subroutine call. (bmy, 10/15/13)
!         ! Read look up tables for PARANOX ship plume model (mpayer, 2/7/12)
!        IF ( ( LICOADSSHIP .or.  LEDGARSHIP .or. 
!     &          LEMEPSHIP ) .and. ITS_A_NEW_MONTH() ) THEN
!            CALL READ_PARANOX_LUT
!         ENDIF
!-----------------------------------------------------------------------------

         ! Seasalt aerosol emission
         IF ( LSSALT ) THEN 
            CALL EMISSSEASALT
     &         ( am_I_Root, Input_Opt, State_Met, 
     &           State_Chm, RC,        ssa_br2    )
         ENDIF

         ! Sulfate aerosol emission
         IF ( LSULF .or. LCRYST ) THEN
            CALL EMISSSULFATE
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF

         ! Carbonaceous/SOA aerosol emission
         IF ( LCARB ) THEN
            CALL EMISSCARBON 
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF

         ! Mineral dust aerosol emission
         IF ( LDUST ) THEN
            CALL EMISSDUST
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF

      ELSE IF ( ITS_A_RnPbBe_SIM ) THEN
         
         !--------------------
         ! Rn-Pb-Be
         !--------------------
         CALL EMISSRnPbBe
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ELSE IF ( ITS_A_CH3I_SIM ) THEN

         !--------------------
         ! CH3I
         !--------------------

         ! Emit CH3I
         CALL EMISSCH3I
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ELSE IF ( ITS_A_HCN_SIM ) THEN

         !--------------------
         ! HCN - CH3CN
         !--------------------
         CALL EMISS_HCN_CH3CN
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ELSE IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN

         !--------------------
         ! Tagged CO
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         ! Bug fix: call every month now (pdk, phs, 3/17/09)
         IF ( Input_Opt%LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! nested
#else
            CALL EMISS_STREETS_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read CAC emissions
         IF ( LCAC .and. ITS_A_NEW_MONTH() ) THEN
#if   defined( GRID05x0666 )
            CALL EMISS_CAC_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! nested
#else
            CALL EMISS_CAC_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read EDGAR emissions once per month
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EDGAR( am_I_Root, YEAR,      MONTH,
     &                        Input_Opt, State_Chm, RC     )
         ENDIF

         ! Read EPA (USA) emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EPA_NEI( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read NEI2005 (USA) emissions once per year
         IF ( LNEI05 .and. ITS_A_NEW_MONTH() ) THEN
#if   defined( GRID05x0666 )
            CALL EMISS_EPA_NEI
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI99 biofuel, nested
            CALL EMISS_NEI2005_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI05 anthro,  nested
#else
            CALL EMISS_EPA_NEI
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI99 biofuel, global
            CALL EMISS_NEI2005_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI05 anthro,  global
#endif
         ENDIF

         ! Read BRAVO (Mexico) emissions once per year
         IF ( LBRAVO .and. ITS_A_NEW_YEAR() ) THEN
            CALL EMISS_BRAVO( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read EPA (Europe) emissions once per year
         IF ( LEMEP  .and. ITS_A_NEW_YEAR()  ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_EMEP_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! nested
#else
            CALL EMISS_EMEP
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read ICOADS ship emissions once per month (cklee, 7/09/09)
         IF ( LICOADSSHIP .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_ICOADS_SHIP( am_I_Root, Input_Opt,
     &                              State_Chm, RC         )
         ENDIF

!-----------------------------------------------------------------------------
! Prior to 10/15/13:
! NOTE: PARANOX does not apply to tagged CO simulations, which do not carry
! NO or O3 tracers.  Disable this subroutine call (bmy, 10/15/13)
!         ! Read look up tables for PARANOX ship plume model (mpayer, 2/7/12)
!         IF ( ( LICOADSSHIP .or. LEDGARSHIP .or. 
!     &          LEMEPSHIP ) .and. ITS_A_NEW_MONTH() ) THEN
!            CALL READ_PARANOX_LUT
!         ENDIF
!-----------------------------------------------------------------------------

         ! Now call EMISSDR for Tagged CO fossil fuel emissions, 
         ! so that we get the same emissions for Tagged CO as 
         ! we do for the full-chemistry (jaf, mak, bmy, 2/14/08)
         CALL EMISSDR 
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         
         ! Emit tagged CO
         CALL EMISS_TAGGED_CO
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ELSE IF ( ITS_A_C2H6_SIM ) THEN

         !--------------------
         ! C2H6
         !--------------------

         ! Emit C2H6
         CALL EMISSC2H6
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ELSE IF ( ITS_A_CH4_SIM ) THEN

         !--------------------
         ! CH4
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         ! Bug fix: call every month now (phs, 3/17/09)
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! nested
#else
            CALL EMISS_STREETS_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Emit CH4
         CALL EMISSCH4
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ELSE IF ( ITS_A_MERCURY_SIM ) THEN

         !--------------------
         ! Mercury
         !--------------------
         CALL EMISSMERCURY
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ELSE IF ( ITS_A_POPS_SIM ) THEN

         !--------------------
         ! POPS
         !--------------------
         CALL EMISSPOPS
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ELSE IF ( ITS_A_CO2_SIM ) THEN

         !--------------------
         ! CO2
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         ! Bug fix: call every month now (phs, 3/17/09)         
         IF ( Input_Opt%LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! nested
#else
            CALL EMISS_STREETS_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read CO2 ARCTAS emissions
         IF ( LARCSHIP .AND. ITS_A_NEW_YEAR() ) THEN
            CALL EMISS_ARCTAS_SHIP
     &         ( am_I_Root, YEAR, Input_Opt, State_Chm, RC )
         ENDIF

         ! Emit CO2
         CALL EMISSCO2
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ELSE IF ( ITS_A_H2HD_SIM ) THEN

         !--------------------
         ! Offline H2/HD 
         !--------------------

         ! Read David Streets' emisisons over China / SE ASia
         ! Bug fix: call every month now (phs, 3/17/09)
         IF ( LSTREETS .and. ITS_A_NEW_MONTH() ) THEN
#if   defined(GRID05x0666)
            CALL EMISS_STREETS_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! nested
#else
            CALL EMISS_STREETS_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read EDGAR emissions once per month
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EDGAR
     &         ( am_I_Root, YEAR, MONTH, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read CAC (Canada) emissions
         IF ( LCAC .and. ITS_A_NEW_MONTH() ) THEN
#if   defined( GRID05x0666 )
            CALL EMISS_CAC_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! nested
#else
            CALL EMISS_CAC_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read EPA (USA) emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_EPA_NEI( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read NEI2005 (USA) emissions
         IF ( Input_Opt%LNEI05 .and. ITS_A_NEW_MONTH() ) THEN
#if   defined( GRID05x0666 )
            CALL EMISS_EPA_NEI
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI99 biofuel, nested
            CALL EMISS_NEI2005_ANTHRO_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI05 anthro,  nested
#else
            CALL EMISS_EPA_NEI
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI99 biofuel, global
            CALL EMISS_NEI2005_ANTHRO
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! NEI05 anthro,  global
#endif
         ENDIF

         ! Read BRAVO (Mexico) emissions once per year
         IF ( LBRAVO .and. ITS_A_NEW_YEAR() ) THEN
            CALL EMISS_BRAVO( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

         ! Read EMEP (Europe) emissions once per year
         IF ( LEMEP  .and. ITS_A_NEW_YEAR()  ) THEN
#if   defined( GRID05x0666 )
            CALL EMISS_EMEP_05x0666
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! nested
#else
            CALL EMISS_EMEP
     &         ( am_I_Root, Input_Opt, State_Chm, RC ) ! global
#endif
         ENDIF

         ! Read ICOADS ship emissions once per month !(cklee, 7/09/09)
         IF ( LICOADSSHIP .and. ITS_A_NEW_MONTH() ) THEN
            CALL EMISS_ICOADS_SHIP
     &         ( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF

!-----------------------------------------------------------------------------
! Prior to 10/15/13:
! NOTE: PARANOX does not apply to tagged CO simulations, which do not carry
! NO or O3 tracers.  Disable this subroutine call (bmy, 10/15/13)
!         ! Read look up tables for PARANOX ship plume model (mpayer, 2/7/12)
!         IF ( ( LICOADSSHIP .or. LEDGARSHIP .or. 
!     &          LEMEPSHIP ) .and. ITS_A_NEW_MONTH() ) THEN
!            CALL READ_PARANOX_LUT
!         ENDIF
!-----------------------------------------------------------------------------

         ! Emit H2/HD
         CALL EMISS_H2_HD
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG ( '### DO_EMISSIONS: a EMISSIONS' )

      END SUBROUTINE DO_EMISSIONS
!EOC
      END MODULE EMISSIONS_MOD
