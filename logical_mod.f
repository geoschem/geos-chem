! $Id: logical_mod.f,v 1.1 2004/09/21 18:04:15 bmy Exp $
      MODULE LOGICAL_MOD
!
!******************************************************************************
!  Module LOGICAL_MOD contains all of the logical switches used by GEOS-CHEM.
!  (bmy, 7/9/04)
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
!  (21) LPRT      (LOGICAL) : Toggles ND70 debug output (via DEBUG_MSG)
!  (22) LSHIPSO2  (LOGICAL) : ON/OFF switch for SO2 EMISSIONS FROM SHIP EXHAUST
!  (23) LSOA      (LOGICAL) : ON/OFF switch for SECONDARY ORGANIC AEROSOLS
!  (24) LSOILNOX  (LOGICAL) : ON/OFF switch for SOIL NOx EMISSIONS
!  (25) LSPLIT    (LOGICAL) : Splits
!  (26) LSSALT    (LOGICAL) : ON/OFF switch for online SEA SALT AEROSOLS
!  (27) LSTDRUN   (LOGICAL) : ON/OFF switch to save init/final masses std runs
!  (28) LSULF     (LOGICAL) : ON/OFF switch for online SULFATE AEROSOLS
!  (29) LSVGLB    (LOGICAL) : ON/OFF switch for SAVING A RESTART FILE
!  (30) LTPFV     (LOGICAL) : If =T, will use fvDAS TPCORE for GEOS-3 winds
!  (31) LTRAN     (LOGICAL) : ON/OFF switch for TRANSPORT  
!  (32) LTOMSAI   (LOGICAL) : ON/OFF switch for scaling biomass w/ TOMS AI data
!  (33) LTURB     (LOGICAL) : ON/OFF switch for PBL MIXING
!  (34) LUNZIP    (LOGICAL) : ON/OFF switch for unzipping met field data 
!  (35) LUPBD     (LOGICAL) : ON/OFF switch for STRATOSPHERIC O3, NOy BC's
!  (36) LWAIT     (LOGICAL) : ON/OFF switch for unzipping met fields in fg
!  (37) LWETD     (LOGICAL) : ON/OFF switch for WET DEPOSITION 
!  (38) LWINDO    (LOGICAL) : ON/OFF switch for WINDOW TRANSPORT (usually 1x1)
!  (39) LWOODCO   (LOGICAL) : ON/OFF switch for BIOFUEL EMISSIONS
!
!  NOTES:
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
