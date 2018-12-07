!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: isoropiaii_mod.F
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
      USE HCO_ERROR_MOD    ! For real precisions (hp)
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_ISOROPIAII
      PUBLIC  :: DO_ISOROPIAII
      PUBLIC  :: GET_GNO3
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
!  23 Sep 2014 - M. Sulprizio - Get global HNO3 for offline aerosol simulation
!                               from HEMCO 
!  14 Nov 2014 - M. Yannetti  - Added PRECISION_MOD
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Array for offline HNO3 (for relaxation of M.M.)
      REAL(fp), ALLOCATABLE :: HNO3_sav(:,:,:)

      ! Array for offline use in sulfate_mod (SEASALT_CHEM) 
      REAL(fp), ALLOCATABLE :: GAS_HNO3(:,:,:)

      ! HEMCO pointers
      REAL(sp), POINTER     :: HNO3(:,:,:) => NULL()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Add a C-preprocessor switch to skip calling ISORROPIA if the pressure
!%%% and/or temperature lie outside of the range that will produce a stable
!%%% solution.  This will eliminate the random noise observed in the
!%%% ISORROPIA output.
!%%%
!%%% Leaving this feature deactivated will replicate the prior behavior in 
!%%% v11-01 and earlier GEOS-Chem versions.  This will become the default
!%%% setting in a future version, but give the user the choice to activate
!%%% or deactivate this for now.
!%%% 
!%%%  -- Seb Eastham and Bob Yantosca (1/25/17)
!%%%
!#define SKIP_IF_P_AND_T_ARE_OUT_OF_RANGE 1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
      SUBROUTINE DO_ISOROPIAII( am_I_Root, Input_Opt,  State_Met,
     &                          State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD 
      USE ErrCode_Mod
      USE ERROR_MOD,            ONLY : DEBUG_MSG
      USE ERROR_MOD,            ONLY : ERROR_STOP
      USE ERROR_MOD,            ONLY : SAFE_DIV
      USE HCO_INTERFACE_MOD,    ONLY : HcoState
      USE HCO_EMISLIST_MOD,     ONLY : HCO_GetPtr
      USE Input_Opt_Mod,        ONLY : OptInput
      USE State_Chm_Mod,        ONLY : ChmState
      USE State_Diag_Mod,       ONLY : DgnState
      USE State_Chm_Mod,        ONLY : Ind_
      USE State_Met_Mod,        ONLY : MetState
      USE TIME_MOD,             ONLY : GET_MONTH
      USE TIME_MOD,             ONLY : ITS_A_NEW_MONTH
      USE IsorropiaII_Main_Mod, ONLY : Isoropia
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD,         ONLY : LD73,   ND73
      USE DIAG_MOD,             ONLY : AD73
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
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
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
!  25 Mar 2013 - M. Payer     - Now pass State_Chm object via the arg list
!  23 Sep 2014 - M. Sulprizio - Now get HNO3 for offline aerosol sim from HEMCO
!  16 Jun 2016 - K. Travis    - Now define species ID's with the Ind_ function 
!  17 Jun 2016 - R. Yantosca  - Now only define species ID's on the first call;
!                               we don't want to do it on each loop iteration
!  30 Jun 2016 - R. Yantosca  - Remove instances of STT.  Now get the advected
!                               species ID from State_Chm%Map_Advect.
!  10 Aug 2016 - R. Yantosca  - Remove temporary tracer-removal code
!  25 Jan 2017 - S. Eastham   - Can skip calling ISORROPIA if P and T are not 
!                               in the proper range to yield a valid solution.
!                               (Activate this with the #ifdef block above.)
!  10 Mar 2017 - C. Keller    - Add print out statement after initialization. 
!  02 Nov 2017 - R. Yantosca  - Now accept State_Diag as an argument
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
      REAL(fp),  PARAMETER       :: CONMIN = 1.0e-30_fp
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER, SAVE            :: id_HNO3, id_NH3,  id_NH4
      INTEGER, SAVE            :: id_NIT,  id_SALA, id_SO4

      ! Scalars
      INTEGER                  :: I,    J,    L,    N
      REAL(fp)                 :: ANO3, GNO3
      REAL(f8)                 :: RHI,  TEMPI, P_Pa
      REAL(fp)                 :: TCA,  TMG,  TK,   HNO3_DEN
      REAL(fp)                 :: TNA,  TCL,  TNH3, TNH4
      REAL(fp)                 :: TNIT, TNO3, TSO4, VOL
      REAL(f8)                 :: AERLIQ(NIONSA+NGASAQA+2)
      REAL(f8)                 :: AERSLD(NSLDSA) 
      REAL(f8)                 :: GAS(NGASAQA) 
      REAL(f8)                 :: OTHER(NOTHERA)
      REAL(f8)                 :: WI(NCOMPA)    
      REAL(f8)                 :: WT(NCOMPA)
      REAL(f8)                 :: CNTRL(NCTRLA)

      ! Strings
      CHARACTER(LEN=15)        :: SCASI  
      CHARACTER(LEN=255)       :: ErrMsg
      CHARACTER(LEN=255)       :: ThisLoc

      !Temporary variables to check if division is safe
      REAL(fp)                 :: NUM_SAV, DEN_SAV
              
      ! AEROPH: Temporary variable for pH (hotp 8/11/09)
      REAL(fp)                 :: HPLUSTEMP
      ! Temporary variable for SO4--
      REAL(fp)                 :: SULFTEMP
      ! Temporary variable for HSO4-
      REAL(fp)                 :: BISULTEMP
      ! Temporary variable for NO3-
      REAL(fp)                 :: NITRTEMP


      ! debug variables
      INTEGER                  :: Itemp, Jtemp, Ltemp
      LOGICAL, SAVE            :: FIRSTCHECK = .TRUE.

      LOGICAL                  :: IT_IS_AN_AEROSOL_SIM
      LOGICAL                  :: IT_IS_A_FULLCHEM_SIM
      LOGICAL                  :: LPRT
                  
      ! Pointers
      REAL(fp), POINTER        :: Spc(:,:,:,:)

      ! For getting fields from HEMCO 
      LOGICAL                  :: aIR

      ! Are we out of the range of valid inputs?
      Logical                  :: OutOfBounds

      !=================================================================
      ! DO_ISOROPIAII begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = 
     &  ' -> at DO_ISOROPIAII (in module GeosCore/isoropiaII_mod.F)'

      ! am I root? 
      aIR                  = am_I_Root

      ! Copy fields from INPUT_OPT to local variables for use below
      LPRT                 = Input_Opt%LPRT  
      IT_IS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM
      IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Make sure certain tracers are defined
         id_HNO3 = Ind_('HNO3') 
         id_NH3  = Ind_('NH3' )
         id_NH4  = Ind_('NH4' )
         id_NIT  = Ind_('NIT' )
         id_SALA = Ind_('SALA')
         id_SO4  = Ind_('SO4' )

         ! Make sure certain tracers are defined
         IF ( id_SO4 <= 0 ) THEN
            ErrMsg = 'SO4 is an undefined species!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         IF ( id_NH3 <= 0 ) THEN
            ErrMsg = 'NH3 is an undefined species!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         IF ( id_NH4 <= 0 ) THEN
            ErrMsg = 'NH4 is an undefined species!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         IF ( id_NIT <= 0 ) THEN
            ErrMsg = 'NIT is an undefined species!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         IF ( id_SALA <= 0 ) THEN
            ErrMsg = 'SALA is an undefined species!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Initialize arrays
         CALL INIT_ISOROPIAII()

         ! Check to see if we need to get HNO3 from HEMCO
         IF ( id_HNO3 <= 0 ) THEN

            IF ( IT_IS_A_FULLCHEM_SIM ) THEN

               ! Coupled simulation: stop w/ error since we need HNO3
               ErrMsg = 'HNO3 is an undefined species!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
 
            ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

               ! Offline simulation: get HNO3 from HEMCO (mps, 9/23/14)
               CALL HCO_GetPtr( aIR, HcoState, 'GLOBAL_HNO3', HNO3, RC )
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 
     &              'Cannot get pointer to HEMCO field GLOBAL_HNO3!'
                  CALL GC_Error( ErrMsg, RC, ThisLoc )
                  RETURN
               ENDIF

            ELSE

               ! ISORROPIA is only valid for full-chem or aerosol-only sims
               ErrMsg = 'Invalid simulation type!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN

            ENDIF
         ENDIF

         ! Print out
         IF ( am_I_Root ) THEN
            WRITE( 6, 100 ) REPEAT( '=', 79 )
            WRITE( 6, 110 )
            WRITE( 6, 100 ) REPEAT( '=', 79 )
         ENDIF

100      FORMAT( a                                            )
110      FORMAT( 'Successfully initialized ISORROPIA code II' )

         ! Reset first-time flag
         FIRST = .FALSE. 
      ENDIF

      ! Initialize for each timestep (bec, bmy, 4/15/05)
      IF ( IT_IS_AN_AEROSOL_SIM ) THEN
         GAS_HNO3 = 0e+0_fp
      ENDIF

      ! Point to chemical species array [kg]
      Spc => State_Chm%Species

      !=================================================================
      ! Loop over grid boxes and call ISOROPIA (see comments in the 
      ! ISOROPIA routine ISOROPIAIICODE.f which describes
      ! the input/output args)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,        L,          N,         WI      )
!$OMP+PRIVATE( WT,      GAS,      TEMPI,      RHI,       VOL     )    
!$OMP+PRIVATE( TSO4,    TNH3,     TNA,        TCL,       ANO3    )
!$OMP+PRIVATE( GNO3,    TCA,      TMG,        TK,        CNTRL   )
!$OMP+PRIVATE( SCASI,   P_PA,     TNO3,       AERLIQ,    AERSLD  )
!$OMP+PRIVATE( OTHER,   TNH4,     TNIT,       HPLUSTEMP, NUM_SAV )
!$OMP+PRIVATE( DEN_SAV, HNO3_DEN, OutOfBounds                    )
!$OMP+PRIVATE( SULFTEMP, BISULTEMP, NITRTEMP                     )
!$OMP+SCHEDULE( DYNAMIC, 1 )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only applying ISORROPIA II in troposphere
         IF ( State_Met%InStratMeso(I,J,L) ) CYCLE

         ! Initialize WI, WT
         DO N = 1, NCOMPA
            WI(N) = 0e+0_fp
            WT(N) = 0e+0_fp
         ENDDO

         ! Initialize GAS
         DO N = 1, NGASAQA
            GAS(N) = 0e+0_fp
         ENDDO

         ! Temperature [K]
         TEMPI    = State_Met%T(I,J,L)

         ! Pressure [Pa]
         P_Pa    = State_Met%PMid(I,J,L)*100.0e+0_f8

         ! Relative humidity [unitless]
         RHI      = State_Met%RH(I,J,L) * 1.e-2_fp

         ! Force RH in the range 0.01 - 0.98
         RHI      = MAX( 0.01e+0_fp, RHI )
         RHI      = MIN( 0.98e+0_fp, RHI )

         ! Volume of grid box [m3] 
         VOL      = State_Met%AIRVOL(I,J,L)

         !---------------------------------
         ! Compute quantities for ISOROPIA
         !---------------------------------

         ! Total SO4 [mole/m3]
         TSO4 = Spc(I,J,L,id_SO4) * 1.e+3_fp / ( 96.e+0_fp * VOL )

         ! Total NH3 [mole/m3] 
         TNH3 = Spc(I,J,L,id_NH4) * 1.e+3_fp / (18.e+0_fp * VOL ) +
     &          Spc(I,J,L,id_NH3) * 1.e+3_fp / (17.e+0_fp * VOL )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: The error-trap statement above will halt execution if id_SALA is
!%%% undefined.  Therefore this IF statement is superfluous.  Comment out 
!%%% for clarity.  (hotp, bmy, 2/1/10)
!%%%
!%%%         IF ( id_SALA > 0 ) THEN
            
            ! Total Na+ (30.61% by weight of seasalt) [mole/m3]
            TNA = Spc(I,J,L,id_SALA) * 0.3061e+0_fp * 1.e+3_fp /
     &        ( 22.99e+0_fp  * VOL  )

            ! Total Cl- (55.04% by weight of seasalt) [mole/m3]
            TCL = Spc(I,J,L,id_SALA) * 0.5504e+0_fp * 1.e+3_fp /
     &           ( 35.45e+0_fp  * VOL  )

!==============================================================================
!=== NOTE: As of 11/2007, ISORROPIAII does not conserve mass when Ca,K,Mg are 
!=== non-zero. If you would like to consider Ca, K, Mg from seasalt and dust, 
!=== isoropiaIIcode.f ISRP4F routines must be debugged.  (hotp, bmy, 2/1/10)
!===
!===            ! Total Ca2+ (1.16% by weight of seasalt) [mole/m3]
!===            TCA      = Spc(I,J,L,id_SALA) * 0.0116e+0_fp * 1.d3 /
!===     &                                 ( 40.08e+0_fp  * VOL  )
!===
!===            ! Total K+   (1.1% by weight of seasalt)  [mole/m3]
!===            TK       = Spc(I,J,L,id_SALA) * 0.0110e+0_fp * 1.d3 /
!===     &                                 ( 39.102e+0_fp * VOL  )
!===
!===            ! Total Mg+  (3.69% by weight of seasalt) [mole/m3]
!===            TMG      = Spc(I,J,L,id_SALA) * 0.0369e+0_fp * 1.d3 /
!===     &                                 ( 24.312e+0_fp * VOL  )

            ! Set Ca, K, Mg to zero for time being (hotp, bmy, 2/1/10)
            TCA      = 0e+0_fp
            TK       = 0e+0_fp
            TMG      = 0e+0_fp
!==============================================================================
!%%%         ELSE
!%%%
!%%%            ! no seasalt, set to zero
!%%%            TNA = 0.e+0_fp
!%%%            TCL = 0.e+0_fp
!%%%            TCA = 0.e+0_fp
!%%%            TK  = 0.e+0_fp
!%%%            TMG = 0.e+0_fp
!%%%
!%%%         ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Compute gas-phase NO3
         IF ( id_HNO3 > 0 ) THEN
            
            !---------------------
            ! COUPLED SIMULATION
            !---------------------

            ! Compute gas-phase HNO3 [mole/m3] from HNO3 tracer
            GNO3 = Spc(I,J,L,id_HNO3)
            GNO3 = MAX( GNO3 * 1.e+3_fp / ( 63.e+0_fp * VOL ), CONMIN )

            ! Aerosol-phase NO3 [mole/m3]
            ANO3 = Spc(I,J,L,id_NIT) * 1.e+3_fp /
     &           (62.e+0_fp * VOL )

            ! Total NO3 [mole/m3]
            TNO3    = GNO3 + ANO3

         ELSE

            !---------------------
            ! OFFLINE SIMULATION
            !---------------------

            ! Convert total inorganic NO3 from [ug/m3] to [mole/m3].
            ! GET_HNO3, lets HNO3 conc's evolve, but relaxes to 
            ! monthly mean values every 3h.
            TNO3  = GET_HNO3( I, J, L, State_Met ) * 1.e-6_fp / 
     &              63.e+0_fp

            ANO3 = 0.0e+0_fp
            GNO3 = TNO3

         ENDIF

         !---------------------------------
         ! Call ISOROPIA
         !---------------------------------

         ! set type of ISOROPIA call
         ! Forward problem, do not change this value
         ! 0e+0_fp represents forward problem
         CNTRL(1) = 0.0e+0_fp

         ! Metastable for now
         ! 1e+0_fp represents metastable problem
         CNTRL(2) = 1.0e+0_fp

         ! Insert concentrations [mole/m3] into WI & prevent underflow
         WI(1)    = MAX( TNA,  CONMIN )
         WI(2)    = MAX( TSO4, CONMIN )
         WI(3)    = MAX( TNH3, CONMIN )
         WI(4)    = MAX( TNO3, CONMIN )
         WI(5)    = MAX( TCL,  CONMIN )
         WI(6)    = MAX( TCA,  CONMIN )
         WI(7)    = MAX( TK,   CONMIN )
         WI(8)    = MAX( TMG,  CONMIN )

#if defined( SKIP_IF_P_AND_T_ARE_OUT_OF_RANGE )

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !%%% If the C-preprocessor switch is activated then check if
         !%%% pressure and temperature are in the range that will result
         !%%% in a stable solution.  If not, then we will skip calling
         !%%% ISORROPIA to avoid random noise in the output.
         !%%%
         !%%% NOTE: Turning this feature on will result in differences
         !%%% with respect to prior GEOS-Chem versions.  So we'll give
         !%%% the user the option to activate it or not.  At some point
         !%%% in the future this will become the default setting.
         !%%%
         !%%%  -- Seb Eastham and Bob Yantosca (1/25/17)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! SDE 2017-01-18: Recommendation from Shannon Capps
         ! Skip equilibrium if T < 250 K or P < 200 hPa
         OutOfBounds = ((P_Pa.lt.200.0e+2_f8).or.(TEMPI.lt.250.0e+0_f8))

#else

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !%%% Always call ISORROPIA, regardless of the values of pressure
         !%%% and temperature.  This will match the prior behavior of
         !%%% when comparing to v11-01 and earlier versions.
         !%%%
         !%%%  -- Seb Eastham and Bob Yantosca (1/25/17)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Never skip calling ISORROPIA
         OutOfBounds = .FALSE.

#endif

         IF ( OutOfBounds ) THEN

            ! %%% Skip equilibrium %%%
            ! Just keep things the way they are 
            ! and spoof the other outputs
            WT(:)     = WI(:)
            AERLIQ(:) = 0.0e+0_f8
            AERSLD(:) = 0.0e+0_f8
            GAS(:)    = 0.0e+0_f8
            OTHER(:)  = 0.0e+0_f8

            ! Separate NH3 and NH4
            TNH3 = Spc(I,J,L,id_NH3) * 1.e+3_fp / (17.e+0_fp * VOL )
            TNH4 = Spc(I,J,L,id_NH4) * 1.e+3_fp / (18.e+0_fp * VOL )

         ELSE

            ! %%% Perform aerosol thermodynamic equilibrium %%%
            ! ISOROPIA can be found in ISOROPIAIICODE.F
            ! inputs are WI, RHI, TEMPI, CNTRL
            CALL ISOROPIA( WI,    RHI,  TEMPI,  CNTRL, 
     &                     WT,    GAS,  AERLIQ, AERSLD,
     &                     SCASI, OTHER                 )

            ! Retrieve concentrations in mol/m3
            TSO4 = WT(2)
            TNH3 = GAS(1)
            TNH4 = WT(3) - GAS(1)
            GNO3 = GAS(2)
            TNO3 = WT(4)
            ANO3 = TNO3 - GNO3

         ENDIF
        
         !---------------------------------
         ! Save back into tracer array
         !---------------------------------
         ! Convert ISOROPIA output from [mole/m3] to [kg]
         TSO4 = MAX( 96.e-3_fp * VOL * TSO4, CONMIN )
         TNH3 = MAX( 17.e-3_fp * VOL * TNH3, CONMIN )
         TNH4 = MAX( 18.e-3_fp * VOL * TNH4, CONMIN )
         TNIT = MAX( 62.e-3_fp * VOL * ANO3, CONMIN )

         ! Save tracers back into Spc array [kg]
         ! no longer save TSO4 back into Spc. SO4 is all aerosol phase
         ! (hotp 11/7/07)
         ! Spc(I,J,L,id_SO4) = TSO4
         Spc(I,J,L,id_NH3) = TNH3
         Spc(I,J,L,id_NH4) = TNH4
         Spc(I,J,L,id_NIT) = TNIT

         ! Special handling for HNO3 [kg]
         IF ( id_HNO3 > 0 ) THEN
            
            !---------------------
            ! COUPLED SIMULATION
            !---------------------

            ! HNO3 [mole/m3] is in GAS(2); convert & store in Spc [kg]
            Spc(I,J,L,id_HNO3)
     &           = MAX( 63.e-3_fp * VOL * GNO3, CONMIN )

            ! Save for use in DEN_SAV expression below (sofen, 4/21/10)
            HNO3_DEN           = Spc(I,J,L,id_HNO3)

         ELSE

            !---------------------
            ! OFFLINE SIMULATION:
            !---------------------

            ! Convert total inorganic nitrate from [mole/m3] to [ug/m3] 
            ! and save for next time
            ! WT(4) is in [mole/m3] -- unit conv is necessary!
            CALL SET_HNO3( I, J, L, 63.e+6_f8 * TNO3 )

            ! Save for use in sulfate_mod (SEASALT_CHEM) for offline
            ! aerosol simulations (bec, 4/15/05)
            GAS_HNO3(I,J,L) = GNO3

            ! Save for use in DEN_SAV expression below (sofen, 4/21/10)
            HNO3_DEN        = GNO3 * VOL * 63e-3_fp

         ENDIF

         !-------------------------
         ! ND42 diagnostic arrays
         !-------------------------

         ! AEROPH: get pH related info to SAV arrays (hotp 8/11/09)
         ! HPLUSTEMP is H+ in mol/L water, AERLIQ1 is H, AERLIQ8 is H2O
         ! in mol/m3 air --> convert to mol/L water
         IF ( AERLIQ(8) < 1e-18_fp ) THEN
            ! Aerosol is dry so HPLUSTEMP and PH_SAV are undefined
            ! We force HPLUSTEMP to 1d20 and PH_SAV to -999e+0_fp.
            ! (hotp, ccc, 12/18/09)
            HPLUSTEMP       = 1e+20_fp
            SULFTEMP        = 1e-30_fp
            BISULTEMP       = 1e-30_fp
            NITRTEMP        = 1e-30_fp
            State_Chm%pHSav(I,J,L) = -999e+0_fp
         ELSE
            HPLUSTEMP       = AERLIQ(1) / AERLIQ(8) * 1e+3_fp/18e+0_fp
            SULFTEMP        = AERLIQ(5) / AERLIQ(8) * 1e+3_fp/18e+0_fp
            BISULTEMP       = AERLIQ(6) / AERLIQ(8) * 1e+3_fp/18e+0_fp
            NITRTEMP        = AERLIQ(7) / AERLIQ(8) * 1e+3_fp/18e+0_fp

            ! Use SAFELOG10 to prevent NAN
            State_Chm%pHSav(I,J,L) = -1e+0_fp * SAFELOG10( HPLUSTEMP )
         ENDIF

         ! Additional Info
         State_Chm%HplusSav(I,J,L)  = HPLUSTEMP
         State_Chm%WaterSav(I,J,L)  = AERLIQ(8) * 18e+6_fp ! mol/m3 -> ug/m3
         State_Chm%SulratSav(I,J,L) = SULFTEMP
         State_Chm%NaRatSav(I,J,L)  = NITRTEMP
         State_Chm%BisulSav(I,J,L)  = BISULTEMP

         NUM_SAV    = ( Spc(I,J,L,id_NH3) /17e+0_fp         +
     &                  Spc(I,J,L,id_NH4) /18e+0_fp         +
     &                  Spc(I,J,L,id_SALA)*0.3061e+0_fp/23.0e+0_fp )

         DEN_SAV    = ( Spc(I,J,L,id_SO4)  / 96e+0_fp  * 2e+0_fp  +
     &                  Spc(I,J,L,id_NIT)  / 62e+0_fp             + 
     &                  HNO3_DEN           / 63e+0_fp             +
     &                  Spc(I,J,L,id_SALA) *0.55e+0_fp / 35.45e+0_fp)

         ! Value if DEN_SAV and NUM_SAV too small.
         State_Chm%AcidPurSav(I,J,L) = SAFE_DIV(NUM_SAV, DEN_SAV,
     &                                          0e+0_fp, 999e+0_fp)

         !==============================================================
         ! DIAGNOSTICS: Archive ISORROPIA H+ and pH values
         !==============================================================

#if defined( BPCH_DIAG )
         !--------------------------------------------------------------
         ! ND73: bpch diagnostic (eam, 04/2015)
         !--------------------------------------------------------------
         IF ( ND73 > 0 ) THEN
            IF ( L <= LD73 ) THEN
               AD73(I,J,L,1) = AD73(I,J,L,1)+State_Chm%pHSav(   I,J,L)
               AD73(I,J,L,2) = AD73(I,J,L,2)+State_Chm%HplusSav(I,J,L)
               AD73(I,J,L,3) = AD73(I,J,L,3)+State_Chm%WaterSav(I,J,L)
            ENDIF
         ENDIF
#endif

      ENDDO
      ENDDO
      ENDDO 
!$OMP END PARALLEL DO

      ! Free pointer
      Spc => NULL()

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### ISOROPIAII: a DO_ISOROPIAII' )

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
      REAL(fp), INTENT(IN) :: X        ! Argument for LOG10 function
!
! !RETURN VALUE:
!
      REAL(fp)             :: SAFLOG   ! LOG10 output -- 
! 
! !REVISION HISTORY: 
!  11 Aug 2009 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( X <= 1e-20_fp ) THEN
          SAFLOG = -1e+0_fp*20e+0_fp   ! if X<0, make pH 20
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
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : AIRMW
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
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
!  23 Sep 2014 - M. Sulprizio - Now get HNO3 for offline aerosol sims from HEMCO
!  06 Jan 2016 - E. Lundgren  - Use global physical parameters
!  06 Feb 2018 - E. Lundgren  - Change GET_ELAPSED_MIN to GET_ELAPSED_SEC to
!                               match new timestep unit of seconds
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: HNO3_UGM3

      !=================================================================
      ! GET_HNO3 begins here!
      !=================================================================

      ! Relax to monthly mean HNO3 concentrations every 3 hours
      ! Otherwise just return the concentration in HNO3_sav
      IF ( MOD( GET_ELAPSED_SEC(), 10800 ) == 0 ) THEN
         ! HNO3 is in v/v (from HEMCO), convert to ug/m3
         ! First convert HNO3 from [v/v] to [kg]
         HNO3_UGM3 = HNO3( I, J, L ) * State_Met%AD(I,J,L) /
     &               ( AIRMW / 63e+0_fp )

         ! Then convert HNO3 from [kg] to [ug/m3]
         HNO3_UGM3 = HNO3_UGM3 * 1.e+9_fp / State_Met%AIRVOL(I,J,L)
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
      REAL(f8),  INTENT(IN) :: HNO3_UGM3   ! HNO3 concentration [ug/m3]
! 
! !REVISION HISTORY: 
!  16 Dec 2002 - R. Yantosca  - Initial version, in ISORROPIA I
!  06 Jul 2007 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!  11 Jan 2015 - M. Yannetti  - Some precision edits, for safety.
!  19 Sep 2016 - R. Yantosca  - Bug fix: #ifdef's should be lower case
!  19 Sep 2016 - R. Yantosca  - Bug fix: elif should be else
!EOP
!------------------------------------------------------------------------------
!BOC
#if !defined( USE_REAL_8 )
      HNO3_sav(I,J,L) = SNGL(HNO3_UGM3) ! if we are not using real*8
#else
      HNO3_sav(I,J,L) = HNO3_UGM3       ! if we are using real*8
#endif

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
      USE State_Met_Mod,      ONLY : MetState
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
      REAL(fp),          INTENT(OUT)  :: HNO3_kg    ! Gas-phase HNO3 [kg]
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
      HNO3_kg  = 0.e+0_fp

      ! convert from [mole/m3] to [kg]
      HNO3_kg = GAS_HNO3(I,J,L) * 63.e-3_fp * State_Met%AIRVOL(I,J,L) 

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
      USE CMN_SIZE_MOD      ! Size parameters
      USE ERROR_MOD,   ONLY : ALLOC_ERR
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
      HNO3_sav = 0e+0_fp

      ALLOCATE( GAS_HNO3( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GAS_HNO3' )
      GAS_HNO3 = 0e+0_fp

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

      ! Free pointers
      IF ( ASSOCIATED( HNO3       ) ) HNO3 => NULL()

      END SUBROUTINE CLEANUP_ISOROPIAII
!EOC
      END MODULE ISOROPIAII_MOD
