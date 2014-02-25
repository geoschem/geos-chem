!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: isoropiaii_mod
!
! !DESCRIPTION: Module ISOROPIAII\_MOD contains the routines that provide
!  the interface between ISORROPIA II and GEOS-Chem.
!\\
!\\
!  The actual ISORROPIA II code which performs Na-SO4-NH3-NO3-Cl-(Ca-K-Mg) 
!  aerosol thermodynamic equilibrium is in \texttt{isoropiaIIcode.f}.
!\\
!\\
! !INTERFACE:
!
      MODULE ISOROPIAII_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_ISOROPIAII
      PUBLIC  :: DO_ISOROPIAII
      PUBLIC  :: GET_GNO3
      PUBLIC  :: GET_ISRINFO
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GET_HNO3
      PRIVATE :: INIT_ISOROPIAII
      PRIVATE :: SAFELOG10
      PRIVATE :: SET_HNO3
!
! !REMARKS:
!  Original Author:
!  *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
!  *** GEORGIA INSTITUTE OF TECHNOLOGY
!  *** WRITTEN BY ATHANASIOS NENES
!  *** UPDATED BY CHRISTOS FOUNTOUKIS
!                                                                             .
!  Original v1.3 isoropia implementation into GEOS-Chem by
!  Becky Alexander and Bob Yantosca (bec, bmy, 4/12/05, 11/2/05)
!                                                                             .
!  For Ca,K,Mg = 0, ISOROPIA II performs exactly like ISOROPIAv1.7
!  Ca, K, Mg, Na from dust is not currently considered
!                                                                             .
!  To implement ISOROPIA II into GEOS-Chem:
!    * cleanup_isoropiaII needs to be called from cleanup.f
!    * DO_ISOROPIA needs to be replaced with DO_ISOROPIAII in chemistry_mod.f
!    * Change ISOROPIA to ISOROPIAII in sulfate_mod.f
!    * add isoropiaII_mod.f, isoropiaIIcode.f, and irspia.inc to Makefile
!                                                                             .
!  ISOROPIA II implementation notes by Havala O.T. Pye:
!  (1) The original isoropia code from T.Nenes is left as unmodified as 
!       possible. Original isoropia code can be found in isoropiaIIcode.f
!       and common blocks can be found in isrpia.inc. For future upgrades
!       to isoropia, replace isrpia.inc and isoropiaIIcode.f with the new
!       version of isoropia and modify the call to ISOROPIA in this module.
!       Please let the original author know of any changes made to ISOROPIA.
!  (2) As of Nov 2007, routines using non-zero Ca, K, and Mg do not always 
!       conserve mass. Ca, K, and Mg are set to zero. 
!                                                                             .
!  NOTE: ISORROPIA is Greek for "equilibrium", in case you were wondering.
!
! !REVISION HISTORY: 
!  06 Jul 2007 - H. O. T. Pye - Initial version
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!  21 Apr 2010 - R. Yantosca  - Bug fix in DO_ISOROPIAII for offline aerosol
!  16 Feb 2011 - R. Yantosca  - Add modifications for APM from G. Luo
!  20 Aug 2013 - R. Yantosca  - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Array for offline HNO3 (for relaxation of M.M.)
      REAL*8,  ALLOCATABLE :: HNO3_sav(:,:,:)

      ! Array for offline use in sulfate_mod (SEASALT_CHEM) 
      REAL*8,  ALLOCATABLE :: GAS_HNO3(:,:,:)

      ! AEROPH: Save information related to aerosol pH (hotp 8/11/09)
      REAL*8,  ALLOCATABLE :: PH_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: HPLUS_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: WATER_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: SULRAT_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: NARAT_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: ACIDPUR_SAV(:,:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_isoropiaii
!
! !DESCRIPTION: Subroutine DO\_ISOROPIAII is the interface between the 
!  GEOS-Chem model and the aerosol thermodynamical equilibrium routine 
!  ISORROPIA II.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_ISOROPIAII( am_I_Root, Input_Opt, 
     &                          State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_HNO3_MOD,    ONLY : GET_GLOBAL_HNO3
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTHNO3
      USE TRACERID_MOD,       ONLY : IDTNIT
      USE TRACERID_MOD,       ONLY : IDTNH4
      USE TRACERID_MOD,       ONLY : IDTNH3
      USE TRACERID_MOD,       ONLY : IDTSALA
      USE TRACERID_MOD,       ONLY : IDTSO4
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
#if   defined( APM )          
      USE APM_INIT_MOD,       ONLY : NSO4
      USE APM_INIT_MOD,       ONLY : IDTSO4BIN1
      USE APM_INIT_MOD,       ONLY : IDTCTSEA
      USE APM_INIT_MOD,       ONLY : IDTCTBCOC
      USE APM_INIT_MOD,       ONLY : IDTCTDST
      USE APM_INIT_MOD,       ONLY : IDTCTSO4
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Original isoropia v1.3 implmentation: (rjp, bec, bmy, 12/17/01, 8/22/05)
! 
! !REVISION HISTORY: 
!  24 Aug 2007 - H. O. T. Pye - Initial version, in ISORROPIA II
!  18 Dec 2009 - H. O. T. Pye - Added division checks
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!  21 Apr 2010 - E. Sofen     - Prevent out-of-bounds errors for offline
!                               aerosol simulations where HNO3 is undefined
!  23 Jul 2010 - R. Yantosca  - Bug fix: corrected typo in ND42 diag section
!  30 Jul 2012 - R. Yantosca  - Now accept am_I_Root as an argument when
!                               running with the traditional driver main.F
!  31 Jul 2012 - R. Yantosca  - Now loop from 1..LLPAR for GIGC compatibilit
!  14 Nov 2012 - R. Yantosca  - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer     - Replaced all met field arrays with State_Met
!                               derived type object
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Array dimensions
      INTEGER, PARAMETER       :: NOTHERA  =  9
      INTEGER, PARAMETER       :: NCTRLA   =  2
      INTEGER, PARAMETER       :: NCOMPA   =  8
      INTEGER, PARAMETER       :: NIONSA   = 10
      INTEGER, PARAMETER       :: NGASAQA  =  3
      INTEGER, PARAMETER       :: NSLDSA   = 19

      ! Concentration lower limit [mole/m3]
      REAL*8,  PARAMETER       :: CONMIN = 1.0d-30
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER                  :: I,    J,    L,    N
      REAL*8                   :: ANO3, GNO3, RHI,  TEMPI
      REAL*8                   :: TCA,  TMG,  TK,   HNO3_DEN
      REAL*8                   :: TNA,  TCL,  TNH3, TNH4
      REAL*8                   :: TNIT, TNO3, TSO4, VOL
      REAL*8                   :: AERLIQ(NIONSA+NGASAQA+2)
      REAL*8                   :: AERSLD(NSLDSA) 
      REAL*8                   :: GAS(NGASAQA) 
      REAL*8                   :: OTHER(NOTHERA)
      REAL*8                   :: WI(NCOMPA)    
      REAL*8                   :: WT(NCOMPA)
      REAL*8                   :: CNTRL(NCTRLA)
      CHARACTER(LEN=255)       :: X 
      CHARACTER(LEN=15)        :: SCASI  
#if   defined( APM )
      REAL*8                   :: TSO4COAT
#endif

      !Temporary variables to check if division is safe
      REAL*8                   :: NUM_SAV, DEN_SAV
              
      ! AEROPH: Temporary variable for pH (hotp 8/11/09)
      REAL*8                   :: HPLUSTEMP

      ! debug variables
      INTEGER                  :: Itemp, Jtemp, Ltemp
      LOGICAL, SAVE            :: FIRSTCHECK = .TRUE.

      LOGICAL                  :: IT_IS_AN_AEROSOL_SIM
      LOGICAL                  :: IT_IS_A_FULLCHEM_SIM
      LOGICAL                  :: LPRT
                  
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER          :: STT(:,:,:,:)

      !=================================================================
      ! DO_ISOROPIAII begins here!
      !=================================================================

      ! Assume success
      RC                   = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LPRT                 = Input_Opt%LPRT  
      IT_IS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM
      IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM

      ! Location string
      X = 'DO_ISOROPIAII (isoropiaII_mod.f)'

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Make sure certain tracers are defined
         IF ( IDTSO4  == 0 ) CALL ERROR_STOP( 'IDTSO4 is undefined!', X)
         IF ( IDTNH3  == 0 ) CALL ERROR_STOP( 'IDTNH3 is undefined!', X)
         IF ( IDTNH4  == 0 ) CALL ERROR_STOP( 'IDTNH4 is undefined!', X)
         IF ( IDTNIT  == 0 ) CALL ERROR_STOP( 'IDTNIT is undefined!', X)
         IF ( IDTSALA == 0 ) CALL ERROR_STOP( 'IDTSALA is undefined!',X)

         ! Initialize arrays
         CALL INIT_ISOROPIAII

         ! Reset first-time flag
         FIRST = .FALSE. 
      ENDIF

      !=================================================================
      ! Check to see if we have to read in monthly mean HNO3
      !=================================================================
      IF ( IDTHNO3 == 0 ) THEN

         IF ( IT_IS_A_FULLCHEM_SIM ) THEN

            ! Coupled simulation: stop w/ error since we need HNO3
            CALL ERROR_STOP( 'IDTHNO3 is not defined!', X )
 
         ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

            ! Offline simulation: read monthly mean HNO3
            IF ( ITS_A_NEW_MONTH() ) THEN
               CALL GET_GLOBAL_HNO3( GET_MONTH() )
            ENDIF

            ! Initialize for each timestep (bec, bmy, 4/15/05)
            GAS_HNO3 = 0d0

         ELSE

            ! Otherwise stop w/ error
            CALL ERROR_STOP( 'Invalid simulation type!', X )

         ENDIF
      ENDIF

      ! AEROPH: Initialize arrays all the way up to LLPAR for
      ! aeroph. Arrays go up to LLPAR due to ND42 use (hotp 8/11/09)
      PH_SAV      = 0d0
      HPLUS_SAV   = 0d0
      WATER_SAV   = 0d0
      SULRAT_SAV  = 0d0
      NARAT_SAV   = 0d0
      ACIDPUR_SAV = 0d0

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      !=================================================================
      ! Loop over grid boxes and call ISOROPIA (see comments in the 
      ! ISOROPIA routine ISOROPIAIICODE.f which describes
      ! the input/output args)
      !=================================================================

      ! AEROPH: add HPLUSTEMP as private (hotp 8/11/09)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,      L,       N,      WI,   WT,  GAS,  TEMPI )
!$OMP+PRIVATE( RHI,  VOL,    TSO4,    TNH3,   TNA,  TCL, ANO3, GNO3  )
!$OMP+PRIVATE( TCA,  TMG,    TK,      CNTRL,  SCASI                  )
!$OMP+PRIVATE( TNO3, AERLIQ, AERSLD,  OTHER,  TNH4, TNIT             )
!$OMP+PRIVATE( HPLUSTEMP,    NUM_SAV, DEN_SAV, HNO3_DEN              )
#if   defined( APM )
!$OMP+PRIVATE( TSO4COAT             )
#endif
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip strat boxes 
         IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) CYCLE

         ! Initialize WI, WT
         DO N = 1, NCOMPA
            WI(N) = 0d0
            WT(N) = 0d0
         ENDDO

         ! Initialize GAS
         DO N = 1, NGASAQA
            GAS(N) = 0d0
         ENDDO

         ! Temperature [K]
         TEMPI    = State_Met%T(I,J,L)

         ! Relative humidity [unitless]
         RHI      = State_Met%RH(I,J,L) * 1.d-2

         ! Force RH in the range 0.01 - 0.98
         RHI      = MAX( 0.01d0, RHI )
         RHI      = MIN( 0.98d0, RHI )

         ! Volume of grid box [m3] 
         VOL      = State_Met%AIRVOL(I,J,L)

         !---------------------------------
         ! Compute quantities for ISOROPIA
         !---------------------------------

#if   defined( APM )
         
         ! %%% WITH APM %%%

         !GanLuo+ Deduce SULFLO form TSO4, 
         ! Add those sulfate coated on primary particles
         TSO4     = (SUM(STT(I,J,L,IDTSO4BIN1:(IDTSO4BIN1+NSO4-1)))
     &            ) * 1.d3 / ( 96.d0 * VOL )

         TSO4COAT = (STT(I,J,L,(IDTCTBCOC-1+1)) +
     &               STT(I,J,L,(IDTCTBCOC-1+2)) +
     &               STT(I,J,L,IDTCTDST)+STT(I,J,L,IDTCTSEA)) *
     &               1.d3 / ( 96.d0 * VOL )
         TSO4 = TSO4 + TSO4COAT
#else

         !%%% WITHOUT APM %%%

         ! Total SO4 [mole/m3]
         TSO4     = STT(I,J,L,IDTSO4) * 1.d3 / ( 96.d0 * VOL )

#endif
         ! Total NH3 [mole/m3] 
         TNH3     = STT(I,J,L,IDTNH4) * 1.d3 / ( 18.d0 * VOL )  +
     &              STT(I,J,L,IDTNH3) * 1.d3 / ( 17.d0 * VOL )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: The error-trap statement above will halt execution if IDTSALA is
!%%% undefined.  Therefore this IF statement is superfluous.  Comment out 
!%%% for clarity.  (hotp, bmy, 2/1/10)
!%%%
!%%%         IF ( IDTSALA > 0 ) THEN
            
            ! Total Na+ (30.61% by weight of seasalt) [mole/m3]
            TNA      = STT(I,J,L,IDTSALA) * 0.3061d0 * 1.d3 /
     &                                    ( 22.99d0  * VOL  )

            ! Total Cl- (55.04% by weight of seasalt) [mole/m3]
            TCL      = STT(I,J,L,IDTSALA) * 0.5504d0 * 1.d3 /
     &                                    ( 35.45d0  * VOL  )

!==============================================================================
!=== NOTE: As of 11/2007, ISORROPIAII does not conserve mass when Ca,K,Mg are 
!=== non-zero. If you would like to consider Ca, K, Mg from seasalt and dust, 
!=== isoropiaIIcode.f ISRP4F routines must be debugged.  (hotp, bmy, 2/1/10)
!===
!===            ! Total Ca2+ (1.16% by weight of seasalt) [mole/m3]
!===            TCA      = STT(I,J,L,IDTSALA) * 0.0116d0 * 1.d3 /
!===     &                                 ( 40.08d0  * VOL  )
!===
!===            ! Total K+   (1.1% by weight of seasalt)  [mole/m3]
!===            TK       = STT(I,J,L,IDTSALA) * 0.0110d0 * 1.d3 /
!===     &                                 ( 39.102d0 * VOL  )
!===
!===            ! Total Mg+  (3.69% by weight of seasalt) [mole/m3]
!===            TMG      = STT(I,J,L,IDTSALA) * 0.0369d0 * 1.d3 /
!===     &                                 ( 24.312d0 * VOL  )

            ! Set Ca, K, Mg to zero for time being (hotp, bmy, 2/1/10)
            TCA      = 0d0
            TK       = 0d0
            TMG      = 0d0
!==============================================================================
!%%%         ELSE
!%%%
!%%%            ! no seasalt, set to zero
!%%%            TNA = 0.d0
!%%%            TCL = 0.d0
!%%%            TCA = 0.d0
!%%%            TK  = 0.d0
!%%%            TMG = 0.d0
!%%%
!%%%         ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Compute gas-phase NO3
         IF ( IDTHNO3 > 0 ) THEN
            
            !---------------------
            ! COUPLED SIMULATION
            !---------------------

            ! Compute gas-phase HNO3 [mole/m3] from HNO3 tracer
            GNO3  = STT(I,J,L,IDTHNO3)
            GNO3  = MAX( GNO3 * 1.d3 / ( 63.d0 * VOL ), CONMIN )

            ! Aerosol-phase NO3 [mole/m3]
            ANO3     = STT(I,J,L,IDTNIT) * 1.d3 / ( 62.d0 * VOL )

            ! Total NO3 [mole/m3]
            TNO3    = GNO3 + ANO3

         ELSE

            !---------------------
            ! OFFLINE SIMULATION
            !---------------------

            ! Convert total inorganic NO3 from [ug/m3] to [mole/m3].
            ! GET_HNO3, lets HNO3 conc's evolve, but relaxes to 
            ! monthly mean values every 3h.
            TNO3  = GET_HNO3( I, J, L, State_Met ) * 1.d-6 / 63.d0

         ENDIF

         !---------------------------------
         ! Call ISOROPIA
         !---------------------------------

         ! set type of ISOROPIA call
         ! Forward problem, do not change this value
         ! 0d0 represents forward problem
         CNTRL(1) = 0.0d0

         ! Metastable for now
         ! 1d0 represents metastable problem
         CNTRL(2) = 1.0d0

         ! Insert concentrations [mole/m3] into WI & prevent underflow
         WI(1)    = MAX( TNA,  CONMIN )
         WI(2)    = MAX( TSO4, CONMIN )
         WI(3)    = MAX( TNH3, CONMIN )
         WI(4)    = MAX( TNO3, CONMIN )
         WI(5)    = MAX( TCL,  CONMIN )
         WI(6)    = MAX( TCA,  CONMIN )
         WI(7)    = MAX( TK,   CONMIN )
         WI(8)    = MAX( TMG,  CONMIN )

         ! Perform aerosol thermodynamic equilibrium 
         ! ISOROPIA can be found in ISOROPIAIICODE.f
         ! inputs are WI, RHI, TEMPI, CNTRL
         CALL ISOROPIA (WI,    RHI,  TEMPI,  CNTRL, 
     &                  WT,    GAS,  AERLIQ, AERSLD,
     &                  SCASI, OTHER)
        
         !---------------------------------
         ! Save back into tracer array
         !---------------------------------
#if   defined( APM )
         ! Nothing
#else
         ! Convert ISOROPIA output from [mole/m3] to [kg]
         TSO4 = MAX( 96.d-3 * VOL *   WT(2),            CONMIN )
#endif
         TNH3 = MAX( 17.d-3 * VOL *   GAS(1),           CONMIN )
         TNH4 = MAX( 18.d-3 * VOL * ( WT(3) - GAS(1) ), CONMIN )
         TNIT = MAX( 62.d-3 * VOL * ( WT(4) - GAS(2) ), CONMIN )

         ! Save tracers back into STT array [kg]
         ! no longer save TSO4 back into STT. SO4 is all aerosol phase
         ! (hotp 11/7/07)
         ! STT(I,J,L,IDTSO4) = TSO4
         STT(I,J,L,IDTNH3) = TNH3
         STT(I,J,L,IDTNH4) = TNH4
         STT(I,J,L,IDTNIT) = TNIT

         ! Special handling for HNO3 [kg]
         IF ( IDTHNO3 > 0 ) THEN
            
            !---------------------
            ! COUPLED SIMULATION
            !---------------------

            ! HNO3 [mole/m3] is in GAS(2); convert & store in STT [kg]
            STT(I,J,L,IDTHNO3) = MAX( 63.d-3 * VOL * GAS(2), CONMIN )

            ! Save for use in DEN_SAV expression below (sofen, 4/21/10)
            HNO3_DEN           = STT(I,J,L,IDTHNO3)

         ELSE

            !---------------------
            ! OFFLINE SIMULATION:
            !---------------------

            ! Convert total inorganic nitrate from [mole/m3] to [ug/m3] 
            ! and save for next time
            ! WT(4) is in [mole/m3] -- unit conv is necessary!
            CALL SET_HNO3( I, J, L, 63.d6 * WT(4) )

            ! Save for use in sulfate_mod (SEASALT_CHEM) for offline
            ! aerosol simulations (bec, 4/15/05)
            GAS_HNO3(I,J,L) = GAS(2)

            ! Save for use in DEN_SAV expression below (sofen, 4/21/10)
            HNO3_DEN        = GAS(2) * VOL * 63d-3

         ENDIF

         !-------------------------
         ! ND42 diagnostic arrays
         !-------------------------

         ! AEROPH: get pH related info to SAV arrays (hotp 8/11/09)
         ! HPLUSTEMP is H+ in mol/L water, AERLIQ1 is H, AERLIQ8 is H2O
         ! in mol/m3 air --> convert to mol/L water
         IF ( AERLIQ(8) < 1d-32 ) THEN
            ! Aerosol is dry so HPLUSTEMP and PH_SAV are undefined
            ! We force HPLUSTEMP to 1d20 and PH_SAV to -999d0.
            ! (hotp, ccc, 12/18/09)
            HPLUSTEMP       = 1d20
            PH_SAV(I,J,L)   = -999d0
         ELSE
            HPLUSTEMP       = AERLIQ(1) / AERLIQ(8) * 1d3/18d0

            ! Use SAFELOG10 to prevent NAN
            PH_SAV(I,J,L)   = -1d0 * SAFELOG10( HPLUSTEMP )
         ENDIF

         ! Additional Info
         HPLUS_SAV(I,J,L)   = AERLIQ(1)
         WATER_SAV(I,J,L)   = AERLIQ(8)
         SULRAT_SAV(I,J,L)  = OTHER(2)
         NARAT_SAV(I,J,L)   = OTHER(4)

         NUM_SAV            = ( STT(I,J,L,IDTNH3) /17d0         +
     &                          STT(I,J,L,IDTNH4) /18d0         +
     &                          STT(I,J,L,IDTSALA)*0.3061d0/23.0d0 )

         DEN_SAV            = ( STT(I,J,L,IDTSO4)  / 96d0   * 2d0     +
     &                          STT(I,J,L,IDTNIT)  / 62d0             + 
     &                          HNO3_DEN           / 63d0             +
     &                          STT(I,J,L,IDTSALA) * 0.55d0 / 35.45d0 )

         ! Value if DEN_SAV and NUM_SAV too small.
         ACIDPUR_SAV(I,J,L) = SAFE_DIV(NUM_SAV, DEN_SAV,
     &                                 0d0,       
     &                                 999d0)
      ENDDO
      ENDDO
      ENDDO 
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### ISOROPIAII: a AERO_THERMO' )

      ! Return to calling program
      END SUBROUTINE DO_ISOROPIAII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: safelog10
!
! !DESCRIPTION: Calculates the LOG (base 10) of a number X.  Returns a minimum
!  value if X is too small, in order to avoid NaN or Infinity problems.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SAFELOG10( X ) RESULT ( SAFLOG )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: X        ! Argument for LOG10 function
!
! !RETURN VALUE:
!
      REAL*8             :: SAFLOG   ! LOG10 output -- 
! 
! !REVISION HISTORY: 
!  11 Aug 2009 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( X <= 1d-20 ) THEN
          SAFLOG = -1d0*20d0   ! if X<0, make pH 20
      ELSE
          SAFLOG = LOG10(X)
      ENDIF

      END FUNCTION SAFELOG10
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_isrinfo
!
! !DESCRIPTION: Subroutine GET\_ISRINFO returns information related to 
!  aerosol pH.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_ISRINFO( I, J, L, N ) RESULT ( RETURNVALUE )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I   ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: J   ! GEOS-Chem latitude index
      INTEGER, INTENT(IN) :: L   ! GEOS-Chem level index
      INTEGER, INTENT(IN) :: N   ! Flag for which information is desired
!
! !RETURN VALUE:
!
      REAL*8              :: RETURNVALUE
! 
! !REVISION HISTORY: 
!  11 Aug 2009 - H. O. T. Pye - Initial version
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF     ( N == 1 ) THEN
         RETURNVALUE = PH_SAV( I, J, L )
      ELSEIF ( N == 2 ) THEN
         RETURNVALUE = HPLUS_SAV( I, J, L ) 
      ELSEIF ( N == 3 ) THEN
         RETURNVALUE = WATER_SAV( I, J, L ) 
      ELSEIF ( N == 4 ) THEN
         RETURNVALUE = SULRAT_SAV( I, J, L ) 
      ELSEIF ( N == 5 ) THEN
         RETURNVALUE = NARAT_SAV( I, J, L ) 
      ELSEIF ( N == 6 ) THEN
         RETURNVALUE = ACIDPUR_SAV( I, J, L ) 
      ELSE
         ! return large value to indicate problem
         RETURNVALUE = 99999d0
         !FP_ISOP
         WRITE(*,*) 'VALUE NOT DEFINED IN GET_ISRINFO'
      ENDIF

      END FUNCTION GET_ISRINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_hno3
!
! !DESCRIPTION: Subroutine GET\_HNO3 allows the HNO3 concentrations to evolve 
!  with time, but relaxes back to the monthly mean concentrations every 3 
!  hours.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_HNO3( I, J, L, State_Met ) RESULT ( HNO3_UGM3 )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_HNO3_MOD,    ONLY : GET_HNO3_UGM3
      USE TIME_MOD,           ONLY : GET_ELAPSED_MIN
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I           ! GEOS-Chem longitude index
      INTEGER,        INTENT(IN)  :: J           ! GEOS-Chem latitude index
      INTEGER,        INTENT(IN)  :: L           ! GEOS-Chem level index
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
! 
! !REVISION HISTORY: 
!  16 Dec 2002 - R. Yantosca  - Initial version, in ISORROPIA I
!  24 Mar 2003 - R. Yantosca  - Now use function GET_ELAPSED_MIN() from the 
!                               new "time_mod.f" to get the elapsed minutes 
!                               since the start of run.
!  06 Jul 2007 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8 :: HNO3_UGM3

      !=================================================================
      ! GET_HNO3 begins here!
      !=================================================================

      ! Relax to monthly mean HNO3 concentrations every 3 hours
      ! Otherwise just return the concentration in HNO3_sav
      IF ( MOD( GET_ELAPSED_MIN(), 180 ) == 0 ) THEN
         HNO3_UGM3 = GET_HNO3_UGM3( I, J, L, State_Met )
      ELSE
         HNO3_UGM3 = HNO3_sav(I,J,L)
      ENDIF

      ! Return to calling program
      END FUNCTION GET_HNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_hno3
!
! !DESCRIPTION: Subroutine SET\_HNO3 stores the modified HNO3 value back 
!  into the HNO3\_sav array for the next timestep.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_HNO3( I, J, L, HNO3_UGM3 )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I           ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: J           ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: L           ! GEOS-Chem longitude index
      REAL*8,  INTENT(IN) :: HNO3_UGM3   ! HNO3 concentration [ug/m3]
! 
! !REVISION HISTORY: 
!  16 Dec 2002 - R. Yantosca  - Initial version, in ISORROPIA I
!  06 Jul 2007 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      HNO3_sav(I,J,L) = HNO3_UGM3

      END SUBROUTINE SET_HNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gno3
!
! !DESCRIPTION: Function GET\_GNO3 returns the gas-phase HNO3 [v/v] for 
!  calculation of sea-salt chemistry in sulfate\_mod (SEASALT\_CHEM). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GNO3( I, J, L, HNO3_kg, State_Met )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: I          ! GEOS-Chem longitude index
      INTEGER,        INTENT(IN)    :: J          ! GEOS-Chem latitude index
      INTEGER,        INTENT(IN)    :: L          ! GEOS-Chem level index
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8,          INTENT(OUT)  :: HNO3_kg    ! Gas-phase HNO3 [kg]
! 
! !REVISION HISTORY: 
!  15 Apr 2005 - B. Alexander - Initial version, in ISORROPIA I
!  06 Jul 2007 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!  09 Nov 2012 - M. Payer     - Replaced all met field arrays with State_Met
!                               derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Zero variables
      HNO3_kg  = 0.D0

      ! convert from [mole/m3] to [kg]
      HNO3_kg = GAS_HNO3(I,J,L) * 63.d-3 * State_Met%AIRVOL(I,J,L) 

      ! Return to calling program
      END SUBROUTINE GET_GNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_isoropiaII
!
! !DESCRIPTION: Subroutine INIT\_ISOROPIAII initializes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ISOROPIAII
!
! !USES:
!
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE CMN_SIZE_MOD                   ! Size parameters
! 
! !REVISION HISTORY: 
!  06 Jul 2007 - H. O. T. Pye - Initial version
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!  31 Jul 2012 - R. Yantosca  - Declare diag arrays w/ LLPAR (not LLTROP)
!EOP
!------------------------------------------------------------------------------
!BOC   
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_ISOROPIAII begins here!
      !=================================================================

      ALLOCATE( HNO3_sav( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HNO3_sav' )
      HNO3_sav = 0d0

      ALLOCATE( GAS_HNO3( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GAS_HNO3' )
      GAS_HNO3 = 0d0

      ! AEROPH: diagnostic info (hotp 8/11/09)
      ! Allocate up to LLPAR, but zero above LLTROP
      ALLOCATE( PH_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PH_SAV' )
      PH_SAV = 0d0

      ALLOCATE( HPLUS_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0  ) CALL ALLOC_ERR( 'HPLUS_SAV' )
      HPLUS_SAV = 0d0
      
      ALLOCATE( WATER_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WATER_SAV' )
      WATER_SAV = 0d0

      ALLOCATE( SULRAT_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SULRAT_SAV' )
      SULRAT_SAV = 0d0

      ALLOCATE( NARAT_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NARAT_SAV' )
      NARAT_SAV = 0d0

      ALLOCATE( ACIDPUR_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ACIDPUR_SAV' )
      ACIDPUR_SAV = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_ISOROPIAII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_isoropiaII
!
! !DESCRIPTION: Subroutine CLEANUP\_ISOROPIAII deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ISOROPIAII
! 
! !REVISION HISTORY: 
!  06 Jul 2007 - H. O. T. Pye - Initial version
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( ALLOCATED( HNO3_sav    ) ) DEALLOCATE( HNO3_sav )
      IF ( ALLOCATED( GAS_HNO3    ) ) DEALLOCATE( GAS_HNO3 )
      ! AEROPH: Deallocate arrays for pH (hotp 8/11/09)
      IF ( ALLOCATED( PH_SAV      ) ) DEALLOCATE( PH_SAV     )
      IF ( ALLOCATED( HPLUS_SAV   ) ) DEALLOCATE( HPLUS_SAV  )
      IF ( ALLOCATED( WATER_SAV   ) ) DEALLOCATE( WATER_SAV  )
      IF ( ALLOCATED( SULRAT_SAV  ) ) DEALLOCATE( SULRAT_SAV )
      IF ( ALLOCATED( NARAT_SAV   ) ) DEALLOCATE( NARAT_SAV  )
      IF ( ALLOCATED( ACIDPUR_SAV ) ) DEALLOCATE( ACIDPUR_SAV)

      END SUBROUTINE CLEANUP_ISOROPIAII
!EOC
      END MODULE ISOROPIAII_MOD
