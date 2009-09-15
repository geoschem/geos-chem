! $Id: logical_mod.f,v 1.22 2009/09/15 15:51:47 phs Exp $
      MODULE LOGICAL_MOD
!
!******************************************************************************
!  Module LOGICAL_MOD contains all of the logical switches used by GEOS-CHEM.
!  (bmy, 7/9/04, 9/24/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) LAIRNOX   (LOGICAL) : ON/OFF switch for aircraft NOx emissions
!  (2 ) LATEQ     (LOGICAL) : --
!  (3 ) LAVHRRLAI (LOGICAL) : ON/OFF switch for reading AVHRR-derived LAI data
!  (4 ) LBIONOX   (LOGICAL) : ON/OFF switch for biomass burning emissions
!  (5 ) LBBSEA    (LOGICAL) : ON/OFF switch for seasonal biomass emissions 
!  (6 ) LCARB     (LOGICAL) : ON/OFF switch for ONLINE CARBONACEOUS AEROSOLS
!  (7 ) LCHEM     (LOGICAL) : ON/OFF switch for CHEMISTRY
!  (8 ) LCONV     (LOGICAL) : ON/OFF switch for CLOUD CONVECTION
!  (9 ) LDBUG     (LOGICAL) : --
!  (10) LDEAD     (LOGICAL) : Toggles DEAD (=T) or GOCART (=F) dust emissions
!  (11) LDIAG     (LOGICAL) : -- 
!  (12) LDRYD     (LOGICAL) : ON/OFF switch for DRY DEPOSITION 
!  (13) LDUST     (LOGICAL) : ON/OFF switch for online DUST MOBILIZATION
!  (14) LEMBED    (LOGICAL) : ON/OFF switch for EMBEDDED CHEMISTRY
!  (15) LEMEP     (LOGICAL) : ON/OFF switch for EMEP EUROPEAN EMISSIONS
!  (16) LEMIS     (LOGICAL) : ON/OFF switch for EMISSIONS     
!  (17) LFFNOX    (LOGICAL) : ON/OFF switch for FOSSIL FUEL NOx        
!  (18) LFILL     (LOGICAL) : Argument for TPCORE (transport)
!  (19) LFOSSIL   (LOGICAL) : ON/OFF switch for ANTHROPOGENIC EMISSIONS
!  (20) LLIGHTNOX (LOGICAL) : ON/OFF switch for LIGHTNING NOx EMISSIONS
!  (21) LCTH      (LOGICAL) : ON/OFF switch for CTH LIGHTNING PARAMETERIZATION
!  (22) LMFLUX    (LOGICAL) : ON/OFF switch for MFLUX LIGHTNING PARAMETERIZ'N
!  (23) LPRECON   (LOGICAL) : ON/OFF switch for PRECON LIGHTNING PARAMETERIZ'N
!  (24) LMEGAN    (LOGICAL): ON/OFF switch for MEGAN BIOGENIC EMISSIONS for ISOP
!  (24) LMEGANMONO (LOGICAL): ON/OFF switch for MEGAN BIOGENIC EMISSIONS 
!                             for MONO and MBO.
!  (25) LMFCT     (LOGICAL) : Argument for TPCORE (transport)
!  (26) LMONOT    (LOGICAL) : Scales acetone to monoterpene emission
!  (27) LNEI99    (LOGICAL) : Toggles on EPA/NEI 99 emissions over cont. USA
!  (28) LPRT      (LOGICAL) : Toggles ND70 debug output (via DEBUG_MSG)
!  (29) LSHIPSO2  (LOGICAL) : ON/OFF switch for SO2 EMISSIONS FROM SHIP EXHAUST
!  (30) LSOA      (LOGICAL) : ON/OFF switch for SECONDARY ORGANIC AEROSOLS
!  (31) LSOILNOX  (LOGICAL) : ON/OFF switch for SOIL NOx EMISSIONS
!  (32) LSPLIT    (LOGICAL) : Splits
!  (33) LSSALT    (LOGICAL) : ON/OFF switch for online SEA SALT AEROSOLS
!  (34) LSTDRUN   (LOGICAL) : ON/OFF switch to save init/final masses std runs
!  (35) LSULF     (LOGICAL) : ON/OFF switch for online SULFATE AEROSOLS
!  (36) LSVGLB    (LOGICAL) : ON/OFF switch for SAVING A RESTART FILE
!  (37) LTPFV     (LOGICAL) : If =T, will use fvDAS TPCORE for GEOS-3 winds
!  (38) LTRAN     (LOGICAL) : ON/OFF switch for TRANSPORT  
!  (39) LTOMSAI   (LOGICAL) : ON/OFF switch for scaling biomass w/ TOMS AI data
!  (40) LTURB     (LOGICAL) : ON/OFF switch for PBL MIXING
!  (41) LUNZIP    (LOGICAL) : ON/OFF switch for unzipping met field data 
!  (42) LUPBD     (LOGICAL) : ON/OFF switch for STRATOSPHERIC O3, NOy BC's
!  (43) LWAIT     (LOGICAL) : ON/OFF switch for unzipping met fields in fg
!  (44) LWETD     (LOGICAL) : ON/OFF switch for WET DEPOSITION 
!  (45) LWINDO    (LOGICAL) : ON/OFF switch for WINDOW TRANSPORT (usually 1x1)
!  (46) LWOODCO   (LOGICAL) : ON/OFF switch for BIOFUEL EMISSIONS
!  (47) LDYNOCEAN (LOGICAL) : ON/OFF switch for OCEAN MERCURY MODULE
!  (48) LGFED2BB  (LOGICAL) : ON/OFF switch for GFED2 BIOMASS BURNING 
!  (49) LBRAVO    (LOGICAL) : ON/OFF switch for BRAVO EMISSIONS
!  (50) LEDGAR    (LOGICAL) : ON/OFF switch for EDGAR emissions
!  (51) LEDGARNOx (LOGICAL) : ON/OFF switch for EDGAR NOx emissions
!  (52) LEDGARCO  (LOGICAL) : ON/OFF switch for EDGAR CO emissions
!  (53) LEDGARSHIP(LOGICAL) : ON/OFF switch for EDGAR ship exhaust emissions
!  (54) LEDGARSOx (LOGICAL) : ON/OFF switch for EDGAR SOx emissions
!  (55) LSTREETS  (LOGICAL) : ON/OFF switch for David Streets' emissions
!  (56) LVARTROP  (LOGICAL) : ON/OFF switch for Variable Tropopause
!  (57) LOTDREG   (LOGICAL) : ON/OFF switch for OTD-LIS regional redistribution
!  (57) LOTDLOC   (LOGICAL) : ON/OFF switch for OTD-LIS local    redistribution
!  (58) LOTDSCALE (LOGICAL) : ON/OFF switch for scaling to OTD-LIS climatology
!  (59) LCAC      (LOGICAL) : ON/OFF switch for CAC Canadian anthro emissions
!  (60) LARCSHIP  (LOGICAL) : ON/OFF switch for ARCTAS ship SO2 emissions
!  (61) LEMEPSHIP (LOGICAL) : ON/OFF switch for EMEP ship emissions
!  (62) LVISTAS   (LOGICAL) : ON/OFF switch for VISTAS NOX anthro emissions
!  (63) L8DAYBB   (LOGICAL) : ON/OFF switch for 8-day GFED BB emissions
!  (64) L3HRBB    (LOGICAL) : ON/OFF switch for 3-hr GFED BB emissions
!  (65) LSYNOPBB  (LOGICAL) : ON/OFF switch for synoptic GFED BB emissions
!  (66) LICARTT   (LOGICAL) : ON/OFF switch for modified NEI99-EPA
!
!  (67) LSVCSPEC  (LOGICAL) : ON/OFF switch for using CSPEC restart values     
!
!  (68) LDICARB   (LOGICAL) : ON/OFF switch for SOG condensation 
!                             onto OC aerosols
!  (69) LCOOKE    (LOGICAL) : ON/OFF switch for overwritting OC/BC emissions
!                             from BOND with COOKE data over North America
!
!  SWITCHES FOR CH4 SIMULATION
!  ---------------------------
!  (70) LGAO      (LOGICAL) : ON/OFF switch for using Gas and Oil emissions
!  (71) LCOL      (LOGICAL) : ON/OFF switch for using Coal emissions
!  (72) LLIV      (LOGICAL) : ON/OFF switch for using Livestock emissions
!  (73) LWAST     (LOGICAL) : ON/OFF switch for using Waste emissions
!  (74) LRICE     (LOGICAL) : ON/OFF switch for using Rice emissions
!  (75) LOTANT    (LOGICAL) : ON/OFF switch for using Ot. Anthropogenic emiss.
!  (76) LWETL     (LOGICAL) : ON/OFF switch for using Wetland emissions
!  (77) LSOABS    (LOGICAL) : ON/OFF switch for using Soil Absorption
!  (78) LOTNAT    (LOGICAL) : ON/OFF switch for using Ot. Natural emissions
!  (79) LCH4BUD   (LOGICAL) : ON/OFF switch for computing CH4 budget
!
!
!  NOTES:
!  (1 ) Added LNEI99 switch to toggle EPA/NEI emissions (bmy, 11/5/04)
!  (2 ) Added LAVHRRLAI switch to toggle AVHRR LAI fields (bmy, 12/20/04)
!  (3 ) Added LMEGAN switch to toggle MEGAN biogenics (tmf, bmy, 10/20/05)
!  (4 ) Added LEMEP switch to toggle EMEP anthro emissions (bdf, bmy, 11/1/05)
!  (5 ) Added LDYNOCEAN switch for online ocean Hg model (bmy, 2/24/06)
!  (6 ) Added LGFED2BB switch for GFED2 BIOMASS BURNING (bmy, 4/5/06)
!  (7 ) Added LCTH, LMFLUX, LPRECON for lightning options (ltm, bmy, 5/5/06)
!  (8 ) Added LFUTURE (swu, bmy, 5/30/06)
!  (9 ) Added LBRAVO (rjp, kfb, bmy, 6/26/06)
!  (10) Added LEDGAR, LEDGARNOx, LEDGARCO, LEDGARSHIP, LEDGARSOx switches
!        for EDGAR emissions (avd, bmy, 7/6/06)
!  (11) Added LSTREETS for David Streets' emissions (bmy, 8/17/06)
!  (12) Added LVARTROP for variable tropopause (phs, 8/21/06)
!  (13) Added LOTDREG, LOTDLOC for regional or local OTD-LIS redistribution
!        of lightning flashes. (bmy, 1/31/07)
!  (14) Added LOTDSCALE (ltm, bmy, 9/24/07)
!  (15) Added LCAC, LARCSHIP, LEMEPSHIP (amv, phs, 3/8/08)
!  (16) Added LVISTAS (amv, 11/24/08)
!  (17) Added L8DAYBB, L3HRBB and LSYNOPBB for 8-day and 3-hr GFED BB
!        emissions (yc, phs, 02/12/07)
!  (18) Added LICARTT to account for Hudman corrections to EPA/NEI99
!        (phs, 1/26/09)
!  (19) Added LSVCSPEC (dkh, 02/12/09) 
!  (20) Added LMEGANMONO (ccc, tmf, 3/2/09)
!  (21) Added LDICARB (ccc, tmf, 3/10/09)
!  (22) Add LNLPBL, LARPBLH and LDEPBCK for non-local PBL scheme (lin, 5/29/09)
!  (23) Added LCOOKE (phs, 5/18/09)
!  (24) Added LKPP (phs, 5/28/09)
!  (25) Added LICOADSSHIP (cklee, 06/30/09)
!  (26) Added switches for CH4 emissions and CH4 budget (kjw, 8/18/09)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Aerosols
      LOGICAL :: LATEQ
      LOGICAL :: LCARB
      LOGICAL :: LCRYST
      LOGICAL :: LCOOKE
      LOGICAL :: LDEAD
      LOGICAL :: LDUST
      LOGICAL :: LSULF
      LOGICAL :: LSOA
      LOGICAL :: LSSALT
      LOGICAL :: LDICARB

      ! Chemistry
      LOGICAL :: LCHEM  
      LOGICAL :: LEMBED
      LOGICAL :: LKPP

      ! Cloud convection
      LOGICAL :: LCONV

      ! Diagnostics
      LOGICAL :: LDBUG
      LOGICAL :: LDIAG  
      LOGICAL :: LPRT
      LOGICAL :: LSTDRUN

      ! Dry deposition
      LOGICAL :: LDRYD   

      ! Emissions
      LOGICAL :: LAIRNOX
      LOGICAL :: LANTHRO 
      LOGICAL :: LBBSEA
      LOGICAL :: LBIONOX    ! <-- deprecated: replace w/ LBIOMASS soon
      LOGICAL :: LBIOMASS 
      LOGICAL :: LBIOFUEL
      LOGICAL :: LBIOGENIC
      LOGICAL :: LCAC
      LOGICAL :: LBRAVO
      LOGICAL :: LEDGAR
      LOGICAL :: LEDGARNOx
      LOGICAL :: LEDGARCO
      LOGICAL :: LEDGARSHIP
      LOGICAL :: LEDGARSOx
      LOGICAL :: LEMEP
      LOGICAL :: LEMIS  
      LOGICAL :: LFFNOX                 
      LOGICAL :: LFOSSIL    ! <-- deprecated: replace w/ LANTHRO soon
      LOGICAL :: LSTREETS
      LOGICAL :: LICARTT
      LOGICAL :: LICOADSSHIP    !(cklee, 6/30/09)
      LOGICAL :: LLIGHTNOX
      LOGICAL :: LOTDREG
      LOGICAL :: LOTDLOC
      LOGICAL :: LOTDSCALE  ! (ltm, 9/24/07)
      LOGICAL :: LCTH
      LOGICAL :: LMFLUX
      LOGICAL :: LPRECON
      LOGICAL :: LMEGAN
      LOGICAL :: LMEGANMONO
      LOGICAL :: LMONOT
      LOGICAL :: LNEI99    
      LOGICAL :: LSHIPSO2
      LOGICAL :: LSOILNOX
      LOGICAL :: LTOMSAI
      LOGICAL :: LWOODCO    ! <-- deprecated: replace w/ LBIOFUEL soon
      LOGICAL :: LAVHRRLAI
      LOGICAL :: LGFED2BB
      LOGICAL :: LFUTURE
      LOGICAL :: LARCSHIP 
      LOGICAL :: LEMEPSHIP
      LOGICAL :: LVISTAS
      LOGICAL :: L8DAYBB
      LOGICAL :: L3HRBB
      LOGICAL :: LSYNOPBB

      ! Transport and strat BC's
      LOGICAL :: LFILL
      LOGICAL :: LMFCT
      LOGICAL :: LTRAN   
      LOGICAL :: LTPFV
      LOGICAL :: LUPBD
      LOGICAL :: LWINDO

      ! Met fields
      LOGICAL :: LUNZIP
      LOGICAL :: LWAIT

      ! PBL mixing
      LOGICAL :: LTURB
      LOGICAL :: LNLPBL
      LOGICAL :: LARPBLH ! --> pblh_ar in vdiff_mod.f
      LOGICAL :: LDEPBCK ! --> drydep_back_cons in vdiff_mod.f

      ! Restart file
      LOGICAL :: LSVGLB
      LOGICAL :: LSVCSPEC
  

      ! Tagged simulations
      LOGICAL :: LSPLIT 

      ! Variable Tropopause
      LOGICAL :: LVARTROP

      ! Wet convection
      LOGICAL :: LWETD  

      ! Dynamic ocean mercury model
      LOGICAL :: LDYNOCEAN

      ! CH4 emissions
      LOGICAL :: LGAO
      LOGICAL :: LCOL
      LOGICAL :: LLIV
      LOGICAL :: LWAST
      LOGICAL :: LRICE
      LOGICAL :: LOTANT
      LOGICAL :: LWETL
      LOGICAL :: LSOABS
      LOGICAL :: LOTNAT
      LOGICAL :: LBFCH4
      LOGICAL :: LBMCH4
      LOGICAL :: LCH4BUD

      ! End of module
      END MODULE LOGICAL_MOD
