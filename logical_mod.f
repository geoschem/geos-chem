! $Id: logical_mod.f,v 1.2 2004/12/02 21:48:38 bmy Exp $
      MODULE LOGICAL_MOD
!
!******************************************************************************
!  Module LOGICAL_MOD contains all of the logical switches used by GEOS-CHEM.
!  (bmy, 7/9/04, 11/5/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) LAIRNOX   (LOGICAL) : ON/OFF switch for aircraft NOx emissions
!  (2 ) LATEQ     (LOGICAL) : --
!  (3 ) LBIONOX   (LOGICAL) : ON/OFF switch for biomass burning emissions
!  (4 ) LBBSEA    (LOGICAL) : ON/OFF switch for seasonal biomass emissions 
!  (5 ) LCARB     (LOGICAL) : ON/OFF switch for ONLINE CARBONACEOUS AEROSOLS
!  (6 ) LCHEM     (LOGICAL) : ON/OFF switch for CHEMISTRY
!  (7 ) LCONV     (LOGICAL) : ON/OFF switch for CLOUD CONVECTION
!  (8 ) LDBUG     (LOGICAL) : --
!  (9 ) LDEAD     (LOGICAL) : Toggles DEAD (=T) or GOCART (=F) dust emissions
!  (10) LDIAG     (LOGICAL) : -- 
!  (11) LDRYD     (LOGICAL) : ON/OFF switch for DRY DEPOSITION 
!  (12) LDUST     (LOGICAL) : ON/OFF switch for online DUST MOBILIZATION
!  (13) LEMBED    (LOGICAL) : ON/OFF switch for EMBEDDED CHEMISTRY
!  (14) LEMIS     (LOGICAL) : ON/OFF switch for EMISSIONS     
!  (15) LFFNOX    (LOGICAL) : ON/OFF switch for FOSSIL FUEL NOx        
!  (16) LFILL     (LOGICAL) : Argument for TPCORE (transport)
!  (17) LFOSSIL   (LOGICAL) : ON/OFF switch for ANTHROPOGENIC EMISSIONS
!  (18) LLIGHTNOX (LOGICAL) : ON/OFF switch for LIGHTNING NOx EMISSIONS
!  (19) LMFCT     (LOGICAL) : Argument for TPCORE (transport)
!  (20) LMONOT    (LOGICAL) : Scales acetone to monoterpene emission
!  (21) LNEI99    (LOGICAL) : Toggles on EPA/NEI 99 emissions over cont. USA
!  (22) LPRT      (LOGICAL) : Toggles ND70 debug output (via DEBUG_MSG)
!  (23) LSHIPSO2  (LOGICAL) : ON/OFF switch for SO2 EMISSIONS FROM SHIP EXHAUST
!  (24) LSOA      (LOGICAL) : ON/OFF switch for SECONDARY ORGANIC AEROSOLS
!  (25) LSOILNOX  (LOGICAL) : ON/OFF switch for SOIL NOx EMISSIONS
!  (26) LSPLIT    (LOGICAL) : Splits
!  (27) LSSALT    (LOGICAL) : ON/OFF switch for online SEA SALT AEROSOLS
!  (28) LSTDRUN   (LOGICAL) : ON/OFF switch to save init/final masses std runs
!  (29) LSULF     (LOGICAL) : ON/OFF switch for online SULFATE AEROSOLS
!  (30) LSVGLB    (LOGICAL) : ON/OFF switch for SAVING A RESTART FILE
!  (31) LTPFV     (LOGICAL) : If =T, will use fvDAS TPCORE for GEOS-3 winds
!  (32) LTRAN     (LOGICAL) : ON/OFF switch for TRANSPORT  
!  (33) LTOMSAI   (LOGICAL) : ON/OFF switch for scaling biomass w/ TOMS AI data
!  (34) LTURB     (LOGICAL) : ON/OFF switch for PBL MIXING
!  (35) LUNZIP    (LOGICAL) : ON/OFF switch for unzipping met field data 
!  (36) LUPBD     (LOGICAL) : ON/OFF switch for STRATOSPHERIC O3, NOy BC's
!  (37) LWAIT     (LOGICAL) : ON/OFF switch for unzipping met fields in fg
!  (38) LWETD     (LOGICAL) : ON/OFF switch for WET DEPOSITION 
!  (39) LWINDO    (LOGICAL) : ON/OFF switch for WINDOW TRANSPORT (usually 1x1)
!  (40) LWOODCO   (LOGICAL) : ON/OFF switch for BIOFUEL EMISSIONS
!
!  NOTES:
!  (1 ) Added LNEI99 switch to toggle EPA/NEI emissions (bmy, 11/5/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Aerosols
      LOGICAL :: LATEQ
      LOGICAL :: LCARB
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
      LOGICAL :: LMONOT
      LOGICAL :: LNEI99    
      LOGICAL :: LSHIPSO2
      LOGICAL :: LSOILNOX
      LOGICAL :: LTOMSAI
      LOGICAL :: LWOODCO    ! <-- deprecated: replace w/ LBIOFUEL soon

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
