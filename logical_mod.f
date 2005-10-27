! $Id: logical_mod.f,v 1.5 2005/10/27 14:00:00 bmy Exp $
      MODULE LOGICAL_MOD
!
!******************************************************************************
!  Module LOGICAL_MOD contains all of the logical switches used by GEOS-CHEM.
!  (bmy, 7/9/04, 10/20/05)
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
!  (15) LEMIS     (LOGICAL) : ON/OFF switch for EMISSIONS     
!  (16) LFFNOX    (LOGICAL) : ON/OFF switch for FOSSIL FUEL NOx        
!  (17) LFILL     (LOGICAL) : Argument for TPCORE (transport)
!  (18) LFOSSIL   (LOGICAL) : ON/OFF switch for ANTHROPOGENIC EMISSIONS
!  (19) LLIGHTNOX (LOGICAL) : ON/OFF switch for LIGHTNING NOx EMISSIONS
!  (20) LMEGAN    (LOGICAL) : ON/OFF switch for MEGAN BIOGENIC EMISSIONS
!  (21) LMFCT     (LOGICAL) : Argument for TPCORE (transport)
!  (22) LMONOT    (LOGICAL) : Scales acetone to monoterpene emission
!  (23) LNEI99    (LOGICAL) : Toggles on EPA/NEI 99 emissions over cont. USA
!  (24) LPRT      (LOGICAL) : Toggles ND70 debug output (via DEBUG_MSG)
!  (25) LSHIPSO2  (LOGICAL) : ON/OFF switch for SO2 EMISSIONS FROM SHIP EXHAUST
!  (26) LSOA      (LOGICAL) : ON/OFF switch for SECONDARY ORGANIC AEROSOLS
!  (27) LSOILNOX  (LOGICAL) : ON/OFF switch for SOIL NOx EMISSIONS
!  (28) LSPLIT    (LOGICAL) : Splits
!  (29) LSSALT    (LOGICAL) : ON/OFF switch for online SEA SALT AEROSOLS
!  (30) LSTDRUN   (LOGICAL) : ON/OFF switch to save init/final masses std runs
!  (31) LSULF     (LOGICAL) : ON/OFF switch for online SULFATE AEROSOLS
!  (32) LSVGLB    (LOGICAL) : ON/OFF switch for SAVING A RESTART FILE
!  (33) LTPFV     (LOGICAL) : If =T, will use fvDAS TPCORE for GEOS-3 winds
!  (34) LTRAN     (LOGICAL) : ON/OFF switch for TRANSPORT  
!  (35) LTOMSAI   (LOGICAL) : ON/OFF switch for scaling biomass w/ TOMS AI data
!  (36) LTURB     (LOGICAL) : ON/OFF switch for PBL MIXING
!  (37) LUNZIP    (LOGICAL) : ON/OFF switch for unzipping met field data 
!  (38) LUPBD     (LOGICAL) : ON/OFF switch for STRATOSPHERIC O3, NOy BC's
!  (39) LWAIT     (LOGICAL) : ON/OFF switch for unzipping met fields in fg
!  (40) LWETD     (LOGICAL) : ON/OFF switch for WET DEPOSITION 
!  (41) LWINDO    (LOGICAL) : ON/OFF switch for WINDOW TRANSPORT (usually 1x1)
!  (42) LWOODCO   (LOGICAL) : ON/OFF switch for BIOFUEL EMISSIONS
!
!  NOTES:
!  (1 ) Added LNEI99 switch to toggle EPA/NEI emissions (bmy, 11/5/04)
!  (2 ) Added LAVHRRLAI switch to toggle AVHRR LAI fields (bmy, 12/20/04)
!  (3 ) Added LMEGAN switch to toggle MEGAN biogenics (tmf, bmy, 10/20/05)
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
      LOGICAL :: LDEAD
      LOGICAL :: LDUST
      LOGICAL :: LSULF
      LOGICAL :: LSOA
      LOGICAL :: LSSALT

      ! Chemistry
      LOGICAL :: LCHEM  
      LOGICAL :: LEMBED

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
      LOGICAL :: LBIONOX    ! <-- deprecated: replace w/ LBIOMASS soon
      LOGICAL :: LBIOMASS 
      LOGICAL :: LBIOFUEL
      LOGICAL :: LBIOGENIC
      LOGICAL :: LBBSEA
      LOGICAL :: LEMIS  
      LOGICAL :: LFFNOX                 
      LOGICAL :: LFOSSIL    ! <-- deprecated: replace w/ LANTHRO soon
      LOGICAL :: LLIGHTNOX
      LOGICAL :: LMEGAN
      LOGICAL :: LMONOT
      LOGICAL :: LNEI99    
      LOGICAL :: LSHIPSO2
      LOGICAL :: LSOILNOX
      LOGICAL :: LTOMSAI
      LOGICAL :: LWOODCO    ! <-- deprecated: replace w/ LBIOFUEL soon
      LOGICAL :: LAVHRRLAI

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

      ! Restart file
      LOGICAL :: LSVGLB

      ! Tagged simulations
      LOGICAL :: LSPLIT 

      ! Wet convection
      LOGICAL :: LWETD  

      ! End of module
      END MODULE LOGICAL_MOD
