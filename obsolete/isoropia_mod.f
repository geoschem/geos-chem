! $Id: isoropia_mod.f,v 1.1 2010/02/02 16:57:45 bmy Exp $
      MODULE ISOROPIA_MOD
!
!******************************************************************************
!  Module ISOROPIA_MOD contains the routines from the ISORROPIA package, which
!  performs aerosol thermodynamical equilibrium (bec, bmy, 4/12/05, 11/2/05)
!
!  NOTE: ISORROPIA is Greek for "equilibrium", in case you were wondering. :-)
!
!  Original Author: 
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  Module Variables (NOTE: also see header file "isoropia.h"):
!  ============================================================================
!  (1  ) IACALC           : ?     
!  (2  ) IMAX             : Longitude dimension for BNC* arrays
!  (3  ) IPROB            : IPROB=0 calls the forward solver
!  (4  ) JMAX             : Latitude dimension for BNC* arrays
!  (5  ) MAXIT            : Max # of internal iterations 
!  (6  ) METSTBL          : METSTBL=1 no liquid phase (solid+gas only)
!  (7  ) NCOMP            : ?
!  (8  ) NDIV             : Number of  
!  (9  ) NERRMX           : ? 
!  (10 ) NGASAQ           : ?
!  (11 ) NIONS            : Number of ions 
!  (12 ) NOTHER           : ?
!  (13 ) NPAIR            : ?
!  (14 ) NRHS             : ? 
!  (15 ) NSLDS            : Number of solids 
!  (16 ) NSWEEP           : ?
!  (17 ) NSO4S            : ? 
!  (18 ) NASRD            : ? 
!  (19 ) NZSR             : ? 
!  (20 ) WFTYP            : ?
!  (21 ) EPS              : Tolerance value (1e-3)
!  (22 ) EPSACT           : Another tolerance  (5e-2)
!  (23 ) GREAT            : A very large number (1e10) 
!  (24 ) ONE              : Constant set to value of 1.0
!  (25 ) R                : Gas constant (82.0567d-6
!  (26 ) TINY             : A very small number (1e-20)
!  (27 ) TINY2            : Another very small number (1e-11)
!  (28 ) ZERO             : Constant set to value of 0.0
!  (29 ) BNC198           : Array to hold values for initializaton
!  (30 ) BNC223           : Array to hold values for initializaton
!  (31 ) BNC248           : Array to hold values for initializaton
!  (32 ) BNC273           : Array to hold values for initializaton  
!  (33 ) BNC298           : Array to hold values for initializaton   
!  (34 ) BNC323           : Array to hold values for initializaton   
!  (35 ) ASRAT            : ? 
!  (36 ) ASSO4            : ?
!  (37 ) HNO3_sav         : Array for offline HNO3 (for relaxation to M.M.)
!
!  Module Routines:
!  ============================================================================
!  (1  ) DO_ISOROPIA      : Interface between GEOS-CHEM and ISORROPIA     
!  (2  ) GET_GNO3         : Passes gas-phase HNO3 to "sulfate_mod.f"
!  (3  ) GET_HNO3         : Gets offline HNO3 conc. (relaxed to M.M. ev. 3 hrs)
!  (4  ) SET_HNO3         : Stores offline HNO3 conc. for next timestep
!  (5  ) ISOROPIA         : Calls other ISORROPIA routines
!  (6  ) INIT2            : Initialization routine
!  (7  ) INIT3            : Initialization routine
!  (8  ) GETASR           : Internal ISORROPIA routine
!  (9  ) CALCNA           : Internal ISORROPIA routine
!  (10 ) CALCNH3          : Internal ISORROPIA routine
!  (11 ) CALCNIAQ         : Internal ISORROPIA routine
!  (12 ) CALCNIAQ2        : Internal ISORROPIA routine
!  (13 ) CALCMR           : Internal ISORROPIA routine
!  (14 ) CALCMDRH         : Internal ISORROPIA routine
!  (15 ) CALCMDRP         : Internal ISORROPIA routine
!  (16 ) CALCHS4          : Internal ISORROPIA routine
!  (17 ) CALCPH           : Internal ISORROPIA routine
!  (18 ) CALCACT          : Internal ISORROPIA routine
!  (19 ) G                : Internal ISORROPIA routine
!  (20 ) RSTGAM           : Internal ISORROPIA routine
!  (21 ) KMFUL            : Internal ISORROPIA routine
!  (22 ) MKBI             : Internal ISORROPIA routine
!  (23 ) KMTAB            : Internal ISORROPIA routine
!  (24 ) READ_KMC         : Internal ISORROPIA routine
!  (25 ) READ_BINARY      : Internal ISORROPIA routine
!  (26 ) KM198            : Internal ISORROPIA routine
!  (27 ) KM223            : Internal ISORROPIA routine
!  (28 ) KM248            : Internal ISORROPIA routine
!  (29 ) KM273            : Internal ISORROPIA routine
!  (30 ) KM298            : Internal ISORROPIA routine
!  (31 ) KM323            : Internal ISORROPIA routine
!  (32 ) INIT_KMC         : Internal ISORROPIA routine
!  (33 ) CHRBLN           : Internal ISORROPIA routine
!  (34 ) SHFTRGHT         : Internal ISORROPIA routine
!  (35 ) RPLSTR           : Internal ISORROPIA routine
!  (36 ) PUSHEND          : Internal ISORROPIA routine
!  (37 ) APPENDEXT        : Internal ISORROPIA routine
!  (38 ) POLY3            : Internal ISORROPIA routine
!  (39 ) POLY3B           : Internal ISORROPIA routine
!  (40 ) FUNC             : Internal ISORROPIA routine
!  (41 ) EX10             : Internal ISORROPIA routine
!  (42 ) PUSHERR          : Internal ISORROPIA routine
!  (43 ) ISERRINF         : Internal ISORROPIA routine
!  (44 ) ERRSTAT          : Internal ISORROPIA routine
!  (45 ) ISRP2F           : Internal ISORROPIA routine
!  (46 ) ISRP3F           : Internal ISORROPIA routine
!  (47 ) CALCB4           : Internal ISORROPIA routine
!  (48 ) CALCB3           : Internal ISORROPIA routine
!  (49 ) CALCB3A          : Internal ISORROPIA routine
!  (50 ) FUNCB3A          : Internal ISORROPIA routine
!  (51 ) CALCB3B          : Internal ISORROPIA routine
!  (52 ) CALCB2           : Internal ISORROPIA routine
!  (53 ) CALCB2A          : Internal ISORROPIA routine
!  (54 ) CALCB2A2         : Internal ISORROPIA routine
!  (55 ) CALCB2B          : Internal ISORROPIA routine
!  (56 ) FUNCB2B          : Internal ISORROPIA routine
!  (57 ) CALCB1           : Internal ISORROPIA routine
!  (58 ) CALCB1A          : Internal ISORROPIA routine
!  (59 ) CALCB1B          : Internal ISORROPIA routine
!  (60 ) CALCC2           : Internal ISORROPIA routine
!  (61 ) CALCC1           : Internal ISORROPIA routine
!  (62 ) FUNCC1           : Internal ISORROPIA routine
!  (63 ) CALCD3           : Internal ISORROPIA routine
!  (64 ) FUNCD3           : Internal ISORROPIA routine
!  (65 ) CALCD2           : Internal ISORROPIA routine
!  (66 ) FUNCD2           : Internal ISORROPIA routine
!  (67 ) CALCD1           : Internal ISORROPIA routine
!  (68 ) CALCD1A          : Internal ISORROPIA routine
!  (69 ) CALCG1           : Internal ISORROPIA routine
!  (70 ) CALCG1A          : Internal ISORROPIA routine
!  (71 ) CALCG2           : Internal ISORROPIA routine
!  (72 ) CALCG2A          : Internal ISORROPIA routine
!  (73 ) FUNCG2A          : Internal ISORROPIA routine
!  (74 ) CALCG3           : Internal ISORROPIA routine
!  (75 ) CALCG3A          : Internal ISORROPIA routine
!  (76 ) FUNCG3A          : Internal ISORROPIA routine
!  (77 ) CALCG4           : Internal ISORROPIA routine
!  (78 ) FUNCG4A          : Internal ISORROPIA routine
!  (79 ) CALCG5           : Internal ISORROPIA routine
!  (80 ) FUNCG5A          : Internal ISORROPIA routine
!  (81 ) CALCH1           : Internal ISORROPIA routine
!  (82 ) CALCH1A          : Internal ISORROPIA routine
!  (83 ) CALCH2           : Internal ISORROPIA routine
!  (84 ) CALCH2A          : Internal ISORROPIA routine
!  (85 ) FUNCH2A          : Internal ISORROPIA routine
!  (86 ) CALCH3           : Internal ISORROPIA routine
!  (87 ) FUNCH3A          : Internal ISORROPIA routine
!  (88 ) CALCH4           : Internal ISORROPIA routine
!  (89 ) FUNCH4A          : Internal ISORROPIA routine
!  (90 ) CALCH5           : Internal ISORROPIA routine
!  (91 ) FUNCH5A          : Internal ISORROPIA routine
!  (92 ) CALCH6           : Internal ISORROPIA routine
!  (93 ) FUNCH6A          : Internal ISORROPIA routine
!  (94 ) CALCI1           : Internal ISORROPIA routine
!  (95 ) CALCI1A          : Internal ISORROPIA routine
!  (96 ) CALCI2           : Internal ISORROPIA routine
!  (97 ) CALCI2A          : Internal ISORROPIA routine
!  (98 ) FUNCI2A          : Internal ISORROPIA routine
!  (99 ) CALCI3           : Internal ISORROPIA routine
!  (100) CALCI3A          : Internal ISORROPIA routine
!  (101) FUNCI3A          : Internal ISORROPIA routine
!  (102) FUNCI3B          : Internal ISORROPIA routine
!  (103) CALCI4           : Internal ISORROPIA routine
!  (104) FUNCI4A          : Internal ISORROPIA routine
!  (105) CALCI5           : Internal ISORROPIA routine
!  (106) FUNCI5A          : Internal ISORROPIA routine
!  (107) CALCI6           : Internal ISORROPIA routine
!  (108) CALCJ1           : Internal ISORROPIA routine
!  (109) FUNCJ1           : Internal ISORROPIA routine
!  (110) CALCJ2           : Internal ISORROPIA routine
!  (111) FUNCJ2           : Internal ISORROPIA routine
!  (112) CALCJ3           : Internal ISORROPIA routine
!  (113) CALCNHA          : Internal ISORROPIA routine
!  (114) CALCHA           : Internal ISORROPIA routine
!  (115) INIT_ISOROPIA    : Allocates and zeroes all module arrays
!  (116) CLEANUP_ISOROPIA : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by isoropia_mod.f
!  ============================================================================
!  (1  ) bpch2_mod.f      : Module w/ routines for binary punch file I/O
!  (2  ) dao_mod.f        : Module w/ arrays for DAO met fields
!  (3  ) diag_mod.f       : Module w/ GEOS-CHEM diagnostic arrays
!  (4  ) directory_mod.f  : Module w/ GEOS-CHEM data & met field dirs
!  (5  ) error_mod.f      : Module w/ I/O error and NaN check routines
!  (6  ) grid_mod.f       : Module w/ horizontal grid information
!  (7  ) logical_mod.f    : Module w/ GEOS-CHEM logical switches
!  (8  ) time_mod.f       : Module w/ routines for computing time & date
!  (9  ) tracerid_mod.f   : Module w/ pointers to tracers & emissions
!  (10 ) transfer_mod.f   : Module w/ routines to cast & resize arrays
!  (11 ) tropopause_mod.f : Module w/ routines to read ann mean tropopause
!  
!  NOTES:
!  (1 ) Now references "tropopause_mod.f" (bmy, 8/22/05)
!  (2 ) Remove duplicate declaration of variables, and declaration of local
!        variables that were in "isoropia.h".  This is necessary to compile
!        with the LINUX/IFORT (v9) compiler. (bmy, 11/1/05)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "isoropia_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_ISOROPIA
      PUBLIC :: DO_ISOROPIA
      PUBLIC :: GET_GNO3

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! INTEGER parameters
      INTEGER, PARAMETER   :: IACALC  = 1          
      INTEGER, PARAMETER   :: IMAX   = 741
      INTEGER, PARAMETER   :: IPROB   = 0
      INTEGER, PARAMETER   :: JMAX    = 13
      INTEGER, PARAMETER   :: MAXIT   = 100
      INTEGER, PARAMETER   :: METSTBL = 1
      INTEGER, PARAMETER   :: NCOMP   = 5
      INTEGER, PARAMETER   :: NDIV    = 5
      INTEGER, PARAMETER   :: NERRMX  = 25
      INTEGER, PARAMETER   :: NGASAQ  = 3
      INTEGER, PARAMETER   :: NIONS   = 7
      INTEGER, PARAMETER   :: NOTHER  = 6
      INTEGER, PARAMETER   :: NPAIR   = 13
      INTEGER, PARAMETER   :: NRHS    = 20
      INTEGER, PARAMETER   :: NSLDS   = 9
      INTEGER, PARAMETER   :: NSWEEP  = 4
      INTEGER, PARAMETER   :: NSO4S   = 14
      INTEGER, PARAMETER   :: NASRD   = NSO4S * NRHS
      INTEGER, PARAMETER   :: NZSR    = 100
      INTEGER, PARAMETER   :: WFTYP   = 2 

      ! REAL*8 parameters
      REAL*8,  PARAMETER   :: EPS     = 1.d-3
      REAL*8,  PARAMETER   :: EPSACT  = 5.d-2        
      REAL*8,  PARAMETER   :: GREAT   = 1.d10
      REAL*8,  PARAMETER   :: ONE     = 1.0d0
      REAL*8,  PARAMETER   :: R       = 82.0567d-6
      REAL*8,  PARAMETER   :: TINY    = 1d-20
      REAL*8,  PARAMETER   :: TINY2   = 1.d-11
      REAL*8,  PARAMETER   :: ZERO    = 0.0d0

      ! Allocatable arrays
      REAL*4,  ALLOCATABLE :: BNC198(:,:)
      REAL*4,  ALLOCATABLE :: BNC223(:,:)
      REAL*4,  ALLOCATABLE :: BNC248(:,:)      
      REAL*4,  ALLOCATABLE :: BNC273(:,:)
      REAL*4,  ALLOCATABLE :: BNC298(:,:)
      REAL*4,  ALLOCATABLE :: BNC323(:,:)
      REAL*8,  ALLOCATABLE :: ASRAT(:)
      REAL*8,  ALLOCATABLE :: ASSO4(:)
      REAL*8,  ALLOCATABLE :: HNO3_sav(:,:,:)
      REAL*8,  ALLOCATABLE :: GAS_HNO3(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow after the CONTAINS statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_ISOROPIA
!
!******************************************************************************
!  Subroutine DO_ISOROPIA is the interface between the GEOS-CHEM model
!  and the aerosol thermodynamical equilibrium routine in "rpmares.f"
!  (rjp, bec, bmy, 12/17/01, 8/22/05)
! 
!  NOTES: 
!  (1 ) Aerosol concentrations are all in ug/m^3
!  (2 ) Now references ITS_IN_THE_STRAT from "tropopause_mod.f".  Now remove
!        reference to CMN, it's obsolete. (bmy, 8/22/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,         ONLY : AIRVOL, RH, T
      USE ERROR_MOD,       ONLY : DEBUG_MSG,       ERROR_STOP 
      USE GLOBAL_HNO3_MOD, ONLY : GET_GLOBAL_HNO3, GET_HNO3_UGM3
      USE LOGICAL_MOD,     ONLY : LPRT
      USE TIME_MOD,        ONLY : GET_MONTH,       ITS_A_NEW_MONTH
      USE TRACER_MOD
      USE TRACERID_MOD
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT

#     include "CMN_SIZE"        ! Size aprameters
#     include "isoropia.h"      ! ISOROPIA common blocks

      ! Local variables
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER                  :: I,    J,    L,    N
      REAL*8                   :: ANO3, GNO3, RHI,  TEMPI
      REAL*8                   :: TNa,  TCl,  TNH3, TNH4
      REAL*8                   :: TNIT, TNO3, TSO4, VOL
      REAL*8                   :: AERLIQ(NIONS+NGASAQ+2)
      REAL*8                   :: AERSLD(NSLDS) 
      REAL*8                   :: GAS(NGASAQ) 
      REAL*8                   :: OTHER(NOTHER)
      REAL*8                   :: WI(NCOMP)    
      REAL*8                   :: WT(NCOMP)
      CHARACTER(LEN=255)       :: X 

      ! Concentration lower limit [mole/m3]
      REAL*8,  PARAMETER       :: CONMIN = 1.0d-30
                  
      !=================================================================
      ! DO_ISOROPIA begins here!
      !=================================================================

      ! Location string
      X = 'DO_ISOROPIA (isoropia_mod.f)'

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Make sure certain tracers are defined
         IF ( IDTSO4  == 0 ) CALL ERROR_STOP( 'IDTSO4 is undefined!', X)
         IF ( IDTNH3  == 0 ) CALL ERROR_STOP( 'IDTNH3 is undefined!', X)
         IF ( IDTNH4  == 0 ) CALL ERROR_STOP( 'IDTNH4 is undefined!', X)
         IF ( IDTNIT  == 0 ) CALL ERROR_STOP( 'IDTNIT is undefined!', X)
         IF ( IDTSALA == 0 ) CALL ERROR_STOP( 'IDTSALA is undefined!',X)

         ! Initialize arrays
         CALL INIT_ISOROPIA

         ! Reset first-time flag
         FIRST = .FALSE. 
      ENDIF

      !=================================================================
      ! Check to see if we have to read in monthly mean HNO3
      !=================================================================
      IF ( IDTHNO3 == 0 ) THEN

         IF ( ITS_A_FULLCHEM_SIM() ) THEN

            ! Coupled simulation: stop w/ error since we need HNO3
            CALL ERROR_STOP( 'IDTHNO3 is not defined!', X )
 
         ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

            ! Offline simulation: read monthly mean HNO3
            IF ( ITS_A_NEW_MONTH() ) THEN
               CALL GET_GLOBAL_HNO3( GET_MONTH() )
            ENDIF

            ! Initialize for each timestep (bec, bmy, 4/15/05)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLTROP
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               GAS_HNO3(I,J,L) = 0d0
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE

            ! Otherwise stop w/ error
            CALL ERROR_STOP( 'Invalid simulation type!', X )

         ENDIF
      ENDIF

      !=================================================================
      ! Loop over grid boxes and call ISOROPIA (see comments in the 
      ! ISOROPIA routine below which describe the input/output args)
      !=================================================================
       
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,      L,      N,     WI,   WT,  GAS,  TEMPI )
!$OMP+PRIVATE( RHI,  VOL,    TSO4,   TNH3,  TNa,  TCl, ANO3, GNO3  )
!$OMP+PRIVATE( TNO3, AERLIQ, AERSLD, OTHER, TNH4, TNIT             )
!$OMP+SCHEDULE( DYNAMIC )  
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip strat boxes 
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! coordinates
         I_SAV = I
         J_SAV = J
         L_SAV = L

         ! Initialize WI, WT
         DO N = 1, NCOMP
            WI(N) = 0d0
            WT(N) = 0d0
         ENDDO

         ! Initialize GAS
         DO N = 1, NGASAQ
            GAS(N) = 0d0
         ENDDO

         ! Temperature [K]
         TEMPI    = T(I,J,L)

         ! Relative humidity [unitless]
         RHI      = RH(I,J,L) * 1.d-2

         ! Force RH in the range 0.01 - 0.95
	 RHI      = MAX( 0.01d0, RHI )
	 RHI      = MIN( 0.95d0, RHI )

         ! Volume of grid box [m3] 
         VOL      = AIRVOL(I,J,L)

         !---------------------------------
         ! Compute quantities for ISOROPIA
         !---------------------------------
                                
         ! Total SO4 [mole/m3]
         TSO4     = STT(I,J,L,IDTSO4) * 1.d3 / ( 96.d0 * VOL )

         ! Total NH3 [mole/m3] 
         TNH3     = STT(I,J,L,IDTNH4) * 1.d3 / ( 18.d0 * VOL )  +
     &              STT(I,J,L,IDTNH3) * 1.d3 / ( 17.d0 * VOL )

         IF ( IDTSALA > 0 ) THEN
            
            ! Total Na+ (30.61% by weight of seasalt) [mole/m3]
            TNa      = STT(I,J,L,IDTSALA) * 0.3061d0 * 1.d3 /
     &                                 ( 22.99d0  * VOL  )
	    ! Total Cl- (55.04% by weight of seasalt) [mole/m3]
            TCl      = STT(I,J,L,IDTSALA) * 0.5504d0 * 1.d3 /
     &                                 ( 35.45d0  * VOL  )
	 ELSE

	    ! no seasalt, set to zero
            TNa = ZERO
            TCl = ZERO

         ENDIF

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
            TNO3  = GET_HNO3( I,J,L ) * 1.d-6 / 63.d0

         ENDIF

         !---------------------------------
         ! Call ISOROPIA
         !---------------------------------

         ! Insert concentrations [mole/m3] into WI & prevent underflow
         WI(1)    = MAX( TNA,  CONMIN )
         WI(2)    = MAX( TSO4, CONMIN )
         WI(3)    = MAX( TNH3, CONMIN )
         WI(4)    = MAX( TNO3, CONMIN )
	 WI(5)    = MAX( TCL,  CONMIN )

         ! Perform aerosol thermodynamic equilibrium 
         CALL ISOROPIA( WI,  RHI,    TEMPI,  WT, 
     &                  GAS, AERLIQ, AERSLD, OTHER )

         !---------------------------------
         ! Save back into tracer array
         !---------------------------------
         
         ! Convert ISOROPIA output from [mole/m3] to [kg]
         TSO4 = MAX( 96.d-3 * VOL *   WT(2),            CONMIN )
         TNH3 = MAX( 17.d-3 * VOL *   GAS(1),           CONMIN )
         TNH4 = MAX( 18.d-3 * VOL * ( WT(3) - GAS(1) ), CONMIN )
         TNIT = MAX( 62.d-3 * VOL * ( WT(4) - GAS(2) ), CONMIN )

         ! Save tracers back into STT array [kg]
         STT(I,J,L,IDTSO4) = TSO4
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

	 ENDIF
          
      ENDDO
      ENDDO
      ENDDO 
!$OMP END PARALLEL DO 

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### ISOROPIA: a AERO_THERMO' )

      ! Return to calling program
      END SUBROUTINE DO_ISOROPIA

!------------------------------------------------------------------------------

      FUNCTION GET_HNO3( I, J, L ) RESULT ( HNO3_UGM3 )
!
!******************************************************************************
!  Subroutine GET_HNO3 allows the HNO3 concentrations to evolve with time,
!  but relaxes back to the monthly mean concentrations every 3 hours.
!  (bmy, 12/16/02, 3/24/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box lon, lat, alt indices
!
!  Function Value: 
!  ============================================================================
!  (1  ) HNO3_UGM3 (REAL*8 ) : HNO3 concentration in ug/m3
!
!  NOTES:
!  (1 ) Now use function GET_ELAPSED_MIN() from the new "time_mod.f" to
!        get the elapsed minutes since the start of run. (bmy, 3/24/03)
!******************************************************************************
!
      ! References to F90 modules
      USE GLOBAL_HNO3_MOD, ONLY : GET_HNO3_UGM3
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      REAL*8              :: HNO3_UGM3

      !=================================================================
      ! GET_HNO3 begins here!
      !=================================================================

      ! Relax to monthly mean HNO3 concentrations every 3 hours
      ! Otherwise just return the concentration in HNO3_sav
      IF ( MOD( GET_ELAPSED_MIN(), 180 ) == 0 ) THEN
         HNO3_UGM3 = GET_HNO3_UGM3( I, J, L )
      ELSE
         HNO3_UGM3 = HNO3_sav(I,J,L)
      ENDIF

      ! Return to calling program
      END FUNCTION GET_HNO3

!------------------------------------------------------------------------------

      SUBROUTINE SET_HNO3( I, J, L, HNO3_UGM3 )
!
!******************************************************************************
!  Subroutine SET_HNO3 stores the modified HNO3 value back into the HNO3_sav
!  array for the next timestep. (bmy, 12/16/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box lon, lat, alt indices
!  (4  ) HNO3_UGM3 (REAL*8 ) : HNO3 concentration in ug/m3
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      REAL*8,  INTENT(IN) :: HNO3_UGM3
      
      !=================================================================
      ! SET_HNO3 begins here!
      !=================================================================
      HNO3_sav(I,J,L) = HNO3_UGM3

      ! Return to calling program
      END SUBROUTINE SET_HNO3

!------------------------------------------------------------------------------

      SUBROUTINE GET_GNO3( I, J, L, HNO3_kg )
!
!******************************************************************************
!  Function GET_GNO3 returns the gas-phase HNO3 [v/v] for calculation of 
!  sea-salt chemistry in sulfate_mod (SEASALT_CHEM) (bec, bmy, 4/15/05)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) I, J, L
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) HNO3_kg [kg]
!
!  NOTES:
!******************************************************************************
!
      USE DAO_MOD,      ONLY : AIRVOL, AD

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L 
      
      ! Return value
      REAL*8, INTENT(OUT)  :: HNO3_kg ! [kg]
 
      !=================================================================
      ! GET_GNO3 begins here!
      !=================================================================

      ! Zero variables
      HNO3_kg  = 0.D0

      ! convert from [mole/m3] to [kg]
      HNO3_kg = GAS_HNO3(I,J,L) * 63.d-3 * AIRVOL(I,J,L) 

      ! Return to calling program
      END SUBROUTINE GET_GNO3

!------------------------------------------------------------------------------

      SUBROUTINE ISOROPIA( WI, RHI, TEMPI,  
     &                     WT, GAS, AERLIQ, AERSLD, OTHER )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISOROPIA
! *** THIS SUBROUTINE IS THE MASTER ROUTINE FOR THE ISORROPIA
!     THERMODYNAMIC EQUILIBRIUM AEROSOL MODEL (VERSION 1.1 and above)
!
! ======================== ARGUMENTS / USAGE ===========================
!
!  INPUT:
!  1. [WI] 
!     DOUBLE PRECISION array of length [5].
!     Concentrations, expressed in moles/m3. Depending on the type of
!     problem solved (specified in CNTRL(1)), WI contains either 
!     GAS+AEROSOL or AEROSOL only concentratios.
!     WI(1) - sodium
!     WI(2) - sulfate
!     WI(3) - ammonium
!     WI(4) - nitrate
!     WI(5) - chloride
!
!  2. [RHI] 
!     DOUBLE PRECISION variable.  
!     Ambient relative humidity expressed on a (0,1) scale.
!
!  3. [TEMPI]
!     DOUBLE PRECISION variable. 
!     Ambient temperature expressed in Kelvins. 
!
!  4. [CNTRL]
!     DOUBLE PRECISION array of length [2].
!     Parameters that control the type of problem solved.
!
!     CNTRL(1): Defines the type of problem solved.
!     0 - Forward problem is solved. In this case, array WI contains 
!         GAS and AEROSOL concentrations together.
!     1 - Reverse problem is solved. In this case, array WI contains
!         AEROSOL concentrations only.
!
!     CNTRL(2): Defines the state of the aerosol
!     0 - The aerosol can have both solid+liquid phases (deliquescent)
!     1 - The aerosol is in only liquid state (metastable aerosol)
!
!  OUTPUT:
!  1. [WT] 
!     DOUBLE PRECISION array of length [5].
!     Total concentrations (GAS+AEROSOL) of species, expressed in moles/m3. 
!     If the foreward probelm is solved (CNTRL(1)=0), array WT is 
!     identical to array WI.
!     WT(1) - total sodium
!     WT(2) - total sulfate
!     WT(3) - total ammonium
!     WT(4) - total nitrate
!     WT(5) - total chloride
!
!  2. [GAS]
!     DOUBLE PRECISION array of length [03]. 
!     Gaseous species concentrations, expressed in moles/m3. 
!     GAS(1) - NH3
!     GAS(2) - HNO3
!     GAS(3) - HCl 
!
!  3. [AERLIQ]
!     DOUBLE PRECISION array of length [11]. 
!     Liquid aerosol species concentrations, expressed in moles/m3. 
!     AERLIQ(01) - H+(aq)          
!     AERLIQ(02) - Na+(aq)         
!     AERLIQ(03) - NH4+(aq)
!     AERLIQ(04) - Cl-(aq)         
!     AERLIQ(05) - SO4--(aq)       
!     AERLIQ(06) - HSO4-(aq)       
!     AERLIQ(07) - NO3-(aq)        
!     AERLIQ(08) - H2O             
!     AERLIQ(09) - NH3(aq) (undissociated)
!     AERLIQ(10) - HNCl(aq) (undissociated)
!     AERLIQ(11) - HNO3(aq) (undissociated)
!     AERLIQ(12) - OH-(aq)
!
!  4. [AERSLD]
!     DOUBLE PRECISION array of length [09]. 
!     Solid aerosol species concentrations, expressed in moles/m3. 
!     AERSLD(01) - NaNO3(s)
!     AERSLD(02) - NH4NO3(s)
!     AERSLD(03) - NaCl(s)         
!     AERSLD(04) - NH4Cl(s)
!     AERSLD(05) - Na2SO4(s)       
!     AERSLD(06) - (NH4)2SO4(s)
!     AERSLD(07) - NaHSO4(s)
!     AERSLD(08) - NH4HSO4(s)
!     AERSLD(09) - (NH4)4H(SO4)2(s)
!
!  5. [SCASI]
!     CHARACTER*15 variable.
!     Returns the subcase which the input corresponds to.
!
!  6. [OTHER]
!     DOUBLE PRECISION array of length [6].
!     Returns solution information.
!
!     OTHER(1): Shows if aerosol water exists.
!     0 - Aerosol is WET
!     1 - Aerosol is DRY
!
!     OTHER(2): Aerosol Sulfate ratio, defined as (in moles/m3) :
!               (total ammonia + total Na) / (total sulfate)
!
!     OTHER(3): Sulfate ratio based on aerosol properties that defines 
!               a sulfate poor system:
!               (aerosol ammonia + aerosol Na) / (aerosol sulfate)
!           
!     OTHER(4): Aerosol sodium ratio, defined as (in moles/m3) :
!               (total Na) / (total sulfate)
!      
!     OTHER(5): Ionic strength of the aqueous aerosol (if it exists).
!      
!     OTHER(6): Total number of calls to the activity coefficient 
!               calculation subroutine.
! 
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters
#     include "isoropia.h"  ! ISOROPIA common blocks

      ! Arguments
      REAL*8,  INTENT(IN)  :: RHI, TEMPI
      REAL*8,  INTENT(IN)  :: WI(NCOMP)
      REAL*8,  INTENT(OUT) :: WT(NCOMP)
      REAL*8,  INTENT(OUT) :: GAS(NGASAQ)
      REAL*8,  INTENT(OUT) :: AERSLD(NSLDS) 
      REAL*8,  INTENT(OUT) :: AERLIQ(NIONS+NGASAQ+2)
      REAL*8,  INTENT(OUT) :: OTHER(NOTHER)

      ! Local variables
      INTEGER              :: K

      !=================================================================
      ! ISOROPIA begins here!
      !=================================================================
	
      ! We always solve the FORWARD PROBLEM
      IF ( WI(1) < TINY .AND. WI(5) < TINY ) THEN 

         ! Call ISRP2F if Na+, Cl- are both zero
         CALL ISRP2F( WI, RHI, TEMPI )

      ELSE

         ! Otherwise call ISRP3F for nonzero Na+, Cl-
         CALL ISRP3F (WI, RHI, TEMPI)

      ENDIF

      !=================================================================
      ! Save results to arrays (units = mole/m3)
      !=================================================================

      ! Gaseous aerosol species
      GAS(1) = GNH3                
      GAS(2) = GHNO3
      GAS(3) = GHCL

      ! Liquid aerosol species
      ! This gives AERLIQ(1,3,5-7) H, NH4, SO4, HSO4, NO3
      DO K = 1, NIONS              
         AERLIQ(K) = MOLAL(K)
      ENDDO

      ! This gives AERLIQ(9-11) NH3, HNCl, HNO3 aqeuous
      ! but GASAQ(1-3) is NH3, HNO3, HCl these are always zero
      DO K = 1, NGASAQ
         AERLIQ(NIONS+1+K)   = GASAQ(K)
      ENDDO

      ! This gives AERLIQ(8) H2O
      AERLIQ(NIONS+1)        = WATER*1.0D3/18.0D0

      ! This gives AERLIQ(12) OH-  this is always zero
      AERLIQ(NIONS+NGASAQ+2) = COH

      ! Solid aerosol species 
      AERSLD(1) = CNaNO3           
      AERSLD(2) = CNH4NO3
      AERSLD(3) = CNaCL
      AERSLD(4) = CNH4CL
      AERSLD(5) = CNa2SO4
      AERSLD(6) = CNH42S4
      AERSLD(7) = CNaHSO4
      AERSLD(8) = CNH4HS4
      AERSLD(9) = CLC

      ! Dry flag
      IF( WATER <= TINY ) THEN    
         OTHER(1) = 1.d0
      ELSE
         OTHER(1) = 0.d0
      ENDIF

      ! Other stuff
      OTHER(2) = SULRAT         
      OTHER(3) = SULRATW ! always 2 
      OTHER(4) = SODRAT  ! always zero
      OTHER(5) = IONIC   ! always zero
      OTHER(6) = ICLACT  ! activity counter
	
      ! Total gas+aerosol phase
      WT(1) = WI(1)             
      WT(2) = WI(2)
      WT(3) = WI(3) 
      WT(4) = WI(4)
      WT(5) = WI(5)

      ! For reverse problem
      IF ( IPROB > 0 .AND. WATER > TINY ) THEN 
         WT(3) = WT(3) + GNH3 
         WT(4) = WT(4) + GHNO3
         WT(5) = WT(5) + GHCL
      ENDIF

      ! Return to calling program
      END SUBROUTINE ISOROPIA

!------------------------------------------------------------------------------

      SUBROUTINE INIT2( WI, RHI, TEMPI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE INIT2
! *** THIS SUBROUTINE INITIALIZES ALL GLOBAL VARIABLES FOR AMMONIUM,
!     NITRATE, SULFATE AEROSOL SYSTEMS (SUBROUTINE ISRP2)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters
#     include "isoropia.h"  ! ISOROPIA common blocks

      ! Arguments
      REAL*8, INTENT(IN)   :: WI(NCOMP)
      REAL*8, INTENT(IN)   :: RHI, TEMPI

      ! Local variables
      INTEGER              :: K,  IRH
      REAL*8               :: IC, GII, GI0, XX 
      REAL*8 , PARAMETER   :: LN10=2.3025851d0
      REAL*8               :: T0, T0T, COEF, TCF, G130, G13I

      !=================================================================
      ! INIT2 begins here!
      !=================================================================

      ! Save input variables in common block
      IF (IPROB.EQ.0) THEN                 

         ! Forward calculation
         DO K=1,NCOMP
            W(K) = MAX(WI(K), TINY)
         ENDDO

      ELSE
         
         ! Reverse calculation
         DO 15 K=1,NCOMP                   
            WAER(K) = MAX(WI(K), TINY)
            W(K)    = ZERO
15       CONTINUE
      ENDIF

      RHB     = RHI
      TEMP    = TEMPI

      !=================================================================
      ! Calculate equilibrium constants
      !=================================================================

      XK1  = 1.015d-2   ! HSO4(aq)         <==> H(aq)     + SO4(aq)
      XK21 = 57.639d0   ! NH3(g)           <==> NH3(aq)
      XK22 = 1.805d-5   ! NH3(aq)          <==> NH4(aq)   + OH(aq)
      XK4  = 2.511d6    ! HNO3(g)          <==> H(aq)     + NO3(aq) ! ISORR
      !XK4  = 3.638e6    ! HNO3(g)          <==> H(aq)     + NO3(aq) ! SEQUIL
      XK41 = 2.100d5    ! HNO3(g)          <==> HNO3(aq)
      XK7  = 1.817d0    ! (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
      XK10 = 5.746d-17  ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! ISORR
      !XK10 = 2.985e-17  ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! SEQUIL
      XK12 = 1.382d2    ! NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
      XK13 = 29.268d0   ! (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
      XKW  = 1.010d-14  ! H2O              <==> H(aq)     + OH(aq)

      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0   = 298.15D0
         T0T  = T0/TEMP
         COEF = 1.0+LOG(T0T)-T0T
         XK1  = XK1 *EXP(  8.85d0*(T0T-1.0d0) + 25.14d0*COEF)
         XK21 = XK21*EXP( 13.79d0*(T0T-1.0d0) -  5.393d0*COEF)
         XK22 = XK22*EXP( -1.5d0*(T0T-1.0d0) + 26.92d0*COEF)
         XK4  = XK4 *EXP( 29.17d0*(T0T-1.0d0) + 16.830d0*COEF) !ISORR
         !XK4  = XK4 *EXP( 29.47*(T0T-1.0) + 16.840*COEF) ! SEQUIL
         XK41 = XK41*EXP( 29.17d0*(T0T-1.0d0) + 16.83d0*COEF)
         XK7  = XK7 *EXP( -2.65d0*(T0T-1.0) + 38.570*COEF)
         XK10 = XK10*EXP(-74.38*(T0T-1.0d0) +  6.12d0*COEF) ! ISORR
         !XK10 = XK10*EXP(-75.11*(T0T-1.0) + 13.460*COEF) ! SEQUIL
         XK12 = XK12*EXP( -2.87d0*(T0T-1.0d0) + 15.83d0*COEF)
         XK13 = XK13*EXP( -5.19d0*(T0T-1.0d0) + 54.4d0*COEF)
         XKW  = XKW *EXP(-22.52d0*(T0T-1.0d0) + 26.92d0*COEF)
      ENDIF

      XK2  = XK21*XK22       
      XK42 = XK4/XK41

      !=================================================================
      ! Calculate deliquescence relative humidities (unicomponent)
      !=================================================================
      DRH2SO4  = ZERO
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRNH4NO3 = 0.6183D0
      DRLC     = 0.6900D0

      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15D0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNH4NO3 = DRNH4NO3*EXP(852.d0*TCF)
         DRNH42S4 = DRNH42S4*EXP( 80.d0*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.d0*TCF) 
         DRLC     = DRLC    *EXP(186.d0*TCF) 
      ENDIF

      !=================================================================
      ! Calculate mutual deliquescence relative humidities
      !=================================================================
      DRMLCAB = 0.3780D0              ! (NH4)3H(SO4)2 & NH4HSO4 
      DRMLCAS = 0.6900D0              ! (NH4)3H(SO4)2 & (NH4)2SO4 
      DRMASAN = 0.6000D0              ! (NH4)2SO4     & NH4NO3

      ! comment out for time being
      !IF (INT(TEMP) .NE. 298) THEN    ! For the time being
      !   T0       = 298.15d0
      !   TCF      = 1.0/TEMP - 1.0/T0
      !   DRMLCAB  = DRMLCAB*EXP( 507.506*TCF) 
      !   DRMLCAS  = DRMLCAS*EXP( 133.865*TCF) 
      !   DRMASAN  = DRMASAN*EXP(1269.068*TCF)
      !ENDIF

      !=================================================================
      ! Liquid phase initialization
      !=================================================================
      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY

      DO K = 1, NPAIR
         MOLALR(K) = ZERO
         GAMA(K)   = 0.1d0
         GAMIN(K)  = GREAT
         GAMOU(K)  = GREAT
         M0(K)     = 1d5
      ENDDO

      DO K = 1, NPAIR
         GAMA(K) = 0.1d0
      ENDDO

      DO K = 1, NIONS
         MOLAL(K) = ZERO
      ENDDO
      COH = ZERO

      DO K = 1, NGASAQ
         GASAQ(K) = ZERO
      ENDDO

      !=================================================================
      ! Solid phase initialization
      !=================================================================
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CNaCL   = ZERO
      CNa2SO4 = ZERO
      CNaNO3  = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNaHSO4 = ZERO
      CLC     = ZERO

      !=================================================================
      ! Gas phase
      !=================================================================
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO

      !=================================================================
      ! Calculate ZSR parameters
      !=================================================================
      IRH     = MIN (INT(RHB*NZSR+0.5d0),NZSR)  ! Position in ZSR arrays
      IRH     = MAX (IRH, 1)

      ! NACl
      M0(01) = AWSC(IRH)      
      IF (M0(01) .LT. 100.0d0) THEN
         IC = M0(01)
         CALL KMTAB(IC,298.0d0,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,REAL(TEMP),GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,   GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(01) = M0(01)*EXP(LN10*(GI0-GII))
      ENDIF

      ! (NA)2SO4
      M0(02) = AWSS(IRH)      
      IF (M0(02) .LT. 100.0d0) THEN
         IC = 3.0d0*M0(02)
         CALL KMTAB(IC,298.0d0,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,REAL(TEMP),XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,   XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(02) = M0(02)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NANO3
      M0(03) = AWSN(IRH)      
      IF (M0(03) .LT. 100.0d0) THEN
         IC = M0(03)
         CALL KMTAB(IC,298.0d0,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,REAL(TEMP),XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,   XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(03) = M0(03)*EXP(LN10*(GI0-GII))
      ENDIF

      ! (NH4)2SO4
      M0(04) = AWAS(IRH)      
      IF (M0(04) .LT. 100.0d0) THEN
         IC = 3.0d0*M0(04)
         CALL KMTAB(IC,298.0d0,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,   XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(04) = M0(04)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NH4NO3
      M0(05) = AWAN(IRH)      
      IF (M0(05) .LT. 100.0d0) THEN
         IC     = M0(05)
         CALL KMTAB(IC,298.0d0,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,   XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
         M0(05) = M0(05)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NH4CL
      M0(06) = AWAC(IRH)      
      IF (M0(06) .LT. 100.0d0) THEN
         IC = M0(06)
         CALL KMTAB(IC,298.0d0,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,  XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
         M0(06) = M0(06)*EXP(LN10*(GI0-GII))
      ENDIF

      ! 2H-SO4
      M0(07) = AWSA(IRH)      
      IF (M0(07) .LT. 100.0d0) THEN
         IC = 3.0d0*M0(07)
         CALL KMTAB(IC,298.0d0,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,   XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
         M0(07) = M0(07)*EXP(LN10*(GI0-GII))
      ENDIF

      ! H-HSO4
      M0(08) = AWSA(IRH)      
      !### These are redundant, because M0(8) is not used
      !IF (M0(08) .LT. 100.0) THEN     
      !   IC = M0(08)
      !   CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
      !   CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
      !   CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX)
      !   M0(08) = M0(08)*EXP(LN10*(GI0-GII))
      !ENDIF

      ! NH4HSO4
      M0(09) = AWAB(IRH)      
      IF (M0(09) .LT. 100.0d0) THEN
         IC = M0(09)
         CALL KMTAB(IC,298.0d0,XX,XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX)
!         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
         CALL KMTAB(IC,TEMP,   XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
         M0(09) = M0(09)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NAHSO4
      M0(12) = AWSB(IRH)      
      IF (M0(12) .LT. 100.0d0) THEN
         IC = M0(12)
         CALL KMTAB(IC,298.0d0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GI0)
!         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
         CALL KMTAB(IC,TEMP,  XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
         M0(12) = M0(12)*EXP(LN10*(GI0-GII))
      ENDIF

      ! (NH4)3H(SO4)2
      M0(13) = AWLC(IRH)      
      IF (M0(13) .LT. 100.0d0) THEN
         IC     = 4.0d0*M0(13)
         CALL KMTAB(IC,298.0d0,XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G130   = 0.2d0*(3.0d0*GI0+2.0d0*GII)
!         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         CALL KMTAB(IC,TEMP,   XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G13I   = 0.2d0*(3.0d0*GI0+2.0d0*GII)
!         M0(13) = M0(13)*EXP(LN10*SNGL(G130-G13I))
         M0(13) = M0(13)*EXP(LN10*(G130-G13I))
      ENDIF

      !=================================================================
      ! OTHER INITIALIZATIONS
      !=================================================================
      ICLACT  =  0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
#if   defined( IBM_XLF ) 
      SCASE   =  '\?\?'
#else
      SCASE   =  '??'
#endif
      SULRATW =  2.D0
      SODRAT  =  ZERO
      NOFER   =  0
      STKOFL  = .FALSE.

      DO K = 1, NERRMX
         ERRSTK(K) = -999
         ERRMSG(K) = 'MESSAGE N/A'
      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT2

!------------------------------------------------------------------------------

      SUBROUTINE INIT3( WI, RHI, TEMPI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE INIT3
! *** THIS SUBROUTINE INITIALIZES ALL GLOBAL VARIABLES FOR AMMONIUM,
!     SODIUM, CHLORIDE, NITRATE, SULFATE AEROSOL SYSTEMS (SUBROUTINE 
!     ISRP3)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================

#     include "CMN_SIZE" ! IIPAR, JPAR, LLTROP
#     include "isoropia.h" 

      ! Arguments
      REAL*8,  INTENT(IN)    :: WI(NCOMP) 
      REAL*8,  INTENT(IN)    :: RHI, TEMPI

      ! Local variables ??
      INTEGER           :: K,  IRH
      REAL*8            :: IC, GII, GI0, XX, T0T
!      REAL*8            :: TEMP, RH
      REAL*8            :: TCF, G130, G13I, COEF
      REAL*8, PARAMETER :: LN10 = 2.3025851D0
      REAL*8, PARAMETER :: T0   = 298.15D0

      !=================================================================
      ! INIT3 begins here!
      !=================================================================

      ! Save input variables in common block
      IF (IPROB.EQ.0) THEN   
              
         ! Forward calculation
         DO K=1,NCOMP
            W(K) = MAX(WI(K), TINY)
         ENDDO

      ELSE

         ! Reverse calculation
         DO K=1,NCOMP                   
            WAER(K) = MAX(WI(K), TINY)
            W(K)    = ZERO
         ENDDO
      ENDIF

      RHB     = RHI
      TEMP    = TEMPI

      !=================================================================
      ! Calculate equilibrium constants
      !=================================================================

      XK1  = 1.015D-2  ! HSO4(aq)         <==> H(aq)     + SO4(aq)
      XK21 = 57.639D0  ! NH3(g)           <==> NH3(aq)
      XK22 = 1.805D-5  ! NH3(aq)          <==> NH4(aq)   + OH(aq)
      XK3  = 1.971D6   ! HCL(g)           <==> H(aq)     + CL(aq)
      XK31 = 2.500d3   ! HCL(g)           <==> HCL(aq)
      XK4  = 2.511d6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! ISORR
      !XK4  = 3.638e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! SEQUIL
      XK41 = 2.100d5   ! HNO3(g)          <==> HNO3(aq)
      XK5  = 0.4799D0  ! NA2SO4(s)        <==> 2*NA(aq)  + SO4(aq)
      XK6  = 1.086D-16 ! NH4CL(s)         <==> NH3(g)    + HCL(g)
      XK7  = 1.817D0   ! (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
      XK8  = 37.661D0  ! NACL(s)          <==> NA(aq)    + CL(aq)
      XK10 = 5.746D-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! ISORR
      !XK10 = 2.985e-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! SEQUIL
      XK11 = 2.413D4   ! NAHSO4(s)        <==> NA(aq)    + HSO4(aq)
      XK12 = 1.382D2   ! NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
      XK13 = 29.268D0  ! (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
      XK14 = 22.05D0   ! NH4CL(s)         <==> NH4(aq)   + CL(aq)
      XKW  = 1.010D-14 ! H2O              <==> H(aq)     + OH(aq)
      XK9  = 11.977D0  ! NANO3(s)         <==> NA(aq)    + NO3(aq)

      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0T = T0/TEMP
         COEF= 1.0d0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85d0*(T0T-1.0d0) + 25.14d0*COEF)
         XK21= XK21*EXP( 13.79d0*(T0T-1.0d0) -  5.393d0*COEF)
         XK22= XK22*EXP( -1.5d0*(T0T-1.0d0) + 26.92d0*COEF)
         XK3 = XK3 *EXP( 30.2d0*(T0T-1.0d0) + 19.91d0*COEF)
         XK31= XK31*EXP( 30.2d0*(T0T-1.0d0) + 19.91d0*COEF)
         XK4 = XK4 *EXP( 29.17d0*(T0T-1.0d0) + 16.83d0*COEF) !ISORR
         !XK4 = XK4 *EXP( 29.47*(T0T-1.0) + 16.840*COEF) ! SEQUIL
         XK41= XK41*EXP( 29.17d0*(T0T-1.0d0) + 16.83d0*COEF)
         XK5 = XK5 *EXP(  0.98d0*(T0T-1.0d0) + 39.5d0*COEF)
         XK6 = XK6 *EXP(-71.0d0*(T0T-1.0d0) +  2.4d0*COEF)
         XK7 = XK7 *EXP( -2.65d0*(T0T-1.0d0) + 38.57d0*COEF)
         XK8 = XK8 *EXP( -1.56d0*(T0T-1.0d0) + 16.9d0*COEF)
         XK9 = XK9 *EXP( -8.22d0*(T0T-1.0d0) + 16.01d0*COEF)
         XK10= XK10*EXP(-74.38d0*(T0T-1.0d0) +  6.12d0*COEF) ! ISORR
         !XK10= XK10*EXP(-75.11*(T0T-1.0) + 13.460*COEF) ! SEQUIL
         XK11= XK11*EXP(  0.79d0*(T0T-1.0d0) + 14.746d0*COEF)
         XK12= XK12*EXP( -2.87d0*(T0T-1.0d0) + 15.83d0*COEF)
         XK13= XK13*EXP( -5.19d0*(T0T-1.0d0) + 54.4d0*COEF)
         XK14= XK14*EXP( 24.55d0*(T0T-1.0d0) + 16.9d0*COEF)
         XKW = XKW *EXP(-22.52d0*(T0T-1.0d0) + 26.92d0*COEF)
      ENDIF
      XK2  = XK21*XK22       
      XK42 = XK4/XK41
      XK32 = XK3/XK31

      !=================================================================
      ! Calculate deliquescence relative humidities (unicomponent)
      !=================================================================
      DRH2SO4  = ZERO
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRLC     = 0.6900D0
      DRNACL   = 0.7528D0
      DRNANO3  = 0.7379D0
      DRNH4CL  = 0.7710D0
      DRNH4NO3 = 0.6183D0
      DRNA2SO4 = 0.9300D0
      DRNAHSO4 = 0.5200D0

      IF (INT(TEMP) .NE. 298) THEN
         TCF      = 1.0d0/TEMP - 1.0d0/T0
         DRNACL   = DRNACL  *EXP( 25.d0*TCF)
         DRNANO3  = DRNANO3 *EXP(304.d0*TCF)
         DRNA2SO4 = DRNA2SO4*EXP( 80.d0*TCF)
         DRNH4NO3 = DRNH4NO3*EXP(852.d0*TCF)
         DRNH42S4 = DRNH42S4*EXP( 80.d0*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.d0*TCF) 
         DRLC     = DRLC    *EXP(186.d0*TCF)
         DRNH4CL  = DRNH4Cl *EXP(239.d0*TCF)
         DRNAHSO4 = DRNAHSO4*EXP(-45.d0*TCF) 
      ENDIF

      !=================================================================
      ! Calculate mutual deliquescence relative humidities
      !=================================================================

      DRMLCAB = 0.378D0    ! (NH4)3H(SO4)2 & NH4HSO4 
      DRMLCAS = 0.690D0    ! (NH4)3H(SO4)2 & (NH4)2SO4 
      DRMASAN = 0.600D0    ! (NH4)2SO4     & NH4NO3
      DRMG1   = 0.460D0    ! (NH4)2SO4, NH4NO3, NA2SO4, NH4CL
      DRMG2   = 0.691D0    ! (NH4)2SO4, NA2SO4, NH4CL
      DRMG3   = 0.697D0    ! (NH4)2SO4, NA2SO4
      DRMH1   = 0.240D0    ! NA2SO4, NANO3, NACL, NH4NO3, NH4CL
      DRMH2   = 0.596D0    ! NA2SO4, NANO3, NACL, NH4CL
      DRMI1   = 0.240D0    ! LC, NAHSO4, NH4HSO4, NA2SO4, (NH4)2SO4
      DRMI2   = 0.363D0    ! LC, NAHSO4, NA2SO4, (NH4)2SO4  - NO DATA -
      DRMI3   = 0.610D0    ! LC, NA2SO4, (NH4)2SO4 
      DRMQ1   = 0.494D0    ! (NH4)2SO4, NH4NO3, NA2SO4
      DRMR1   = 0.663D0    ! NA2SO4, NANO3, NACL
      DRMR2   = 0.735D0    ! NA2SO4, NACL
      DRMR3   = 0.673D0    ! NANO3, NACL
      DRMR4   = 0.694D0    ! NA2SO4, NACL, NH4CL
      DRMR5   = 0.731D0    ! NA2SO4, NH4CL
      DRMR6   = 0.596D0    ! NA2SO4, NANO3, NH4CL
      DRMR7   = 0.380D0    ! NA2SO4, NANO3, NACL, NH4NO3
      DRMR8   = 0.380D0    ! NA2SO4, NACL, NH4NO3
      DRMR9   = 0.494D0    ! NA2SO4, NH4NO3
      DRMR10  = 0.476D0    ! NA2SO4, NANO3, NH4NO3
      DRMR11  = 0.340D0    ! NA2SO4, NACL, NH4NO3, NH4CL
      DRMR12  = 0.460D0    ! NA2SO4, NH4NO3, NH4CL
      DRMR13  = 0.438D0    ! NA2SO4, NANO3, NH4NO3, NH4CL
CCC      IF (INT(TEMP) .NE. 298) THEN
CCC         T0       = 298.15d0
CCC         TCF      = 1.0/TEMP - 1.0/T0
CCC         DRMLCAB  = DRMLCAB*EXP( 507.506*TCF) 
CCC         DRMLCAS  = DRMLCAS*EXP( 133.865*TCF) 
CCC         DRMASAN  = DRMASAN*EXP(1269.068*TCF)
CCC         DRMG1    = DRMG1  *EXP( 572.207*TCF)
CCC         DRMG2    = DRMG2  *EXP(  58.166*TCF)
CCC         DRMG3    = DRMG3  *EXP(  22.253*TCF)
CCC         DRMH1    = DRMH1  *EXP(2116.542*TCF)
CCC         DRMH2    = DRMH2  *EXP( 650.549*TCF)
CCC         DRMI1    = DRMI1  *EXP( 565.743*TCF)
CCC         DRMI2    = DRMI2  *EXP(  91.745*TCF)
CCC         DRMI3    = DRMI3  *EXP( 161.272*TCF)
CCC         DRMQ1    = DRMQ1  *EXP(1616.621*TCF)
CCC         DRMR1    = DRMR1  *EXP( 292.564*TCF)
CCC         DRMR2    = DRMR2  *EXP(  14.587*TCF)
CCC         DRMR3    = DRMR3  *EXP( 307.907*TCF)
CCC         DRMR4    = DRMR4  *EXP(  97.605*TCF)
CCC         DRMR5    = DRMR5  *EXP(  98.523*TCF)
CCC         DRMR6    = DRMR6  *EXP( 465.500*TCF)
CCC         DRMR7    = DRMR7  *EXP( 324.425*TCF)
CCC         DRMR8    = DRMR8  *EXP(2660.184*TCF)
CCC         DRMR9    = DRMR9  *EXP(1617.178*TCF)
CCC         DRMR10   = DRMR10 *EXP(1745.226*TCF)
CCC         DRMR11   = DRMR11 *EXP(3691.328*TCF)
CCC         DRMR12   = DRMR12 *EXP(1836.842*TCF)
CCC         DRMR13   = DRMR13 *EXP(1967.938*TCF)
CCC      ENDIF

      !=================================================================
      ! Liquid phase initialization
      !=================================================================
      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY

      DO K=1,NPAIR
         MOLALR(K)=ZERO
         GAMA(K)  =0.1d0
         GAMIN(K) =GREAT
         GAMOU(K) =GREAT
         M0(K)    =1d5
       ENDDO

      DO K=1,NPAIR
         GAMA(K) = 0.1d0
      ENDDO

      DO K=1,NIONS
         MOLAL(K)=ZERO
      ENDDO
      COH = ZERO

      DO K=1,NGASAQ
         GASAQ(K)=ZERO
      ENDDO

      !=================================================================
      ! Solid phase initialization
      !=================================================================
      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO

      !=================================================================
      ! Gas phase
      !=================================================================
      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO

      !=================================================================
      ! Calculate ZSR parameters
      !=================================================================
      IRH    = MIN (INT(RHB*NZSR+0.5d0),NZSR)  ! Position in ZSR arrays
      IRH    = MAX (IRH, 1)

      ! NACl
      IF (M0(01) .LT. 100.0d0) THEN
      M0(01) = AWSC(IRH)
         IC = M0(01)
         CALL KMTAB(IC,298.d0,    GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,SNGL(TEMP),GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(01) = M0(01)*EXP(LN10*(GI0-GII))
      ENDIF

      ! (NA)2SO4
      M0(02) = AWSS(IRH)
      IF (M0(02) .LT. 100.0d0) THEN
         IC = 3.0d0*M0(02)
         CALL KMTAB(IC,298.d0,    XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,SNGL(TEMP),XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,      XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(02) = M0(02)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NANO3 
      M0(03) = AWSN(IRH) 
      IF (M0(03) .LT. 100.0d0) THEN
         IC = M0(03)
         CALL KMTAB(IC,298.d0,    XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,SNGL(TEMP),XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,      XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(03) = M0(03)*EXP(LN10*(GI0-GII))
      ENDIF

      ! (NH4)2SO4
      M0(04) = AWAS(IRH)
      IF (M0(04) .LT. 100.0d0) THEN
         IC = 3.0d0*M0(04)
         CALL KMTAB(IC,298.d0,    XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,      XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(04) = M0(04)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NH4NO3
      M0(05) = AWAN(IRH)
      IF (M0(05) .LT. 100.0d0) THEN
         IC     = M0(05)
         CALL KMTAB(IC,298.d0,    XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,      XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
         M0(05) = M0(05)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NH4CL
      M0(06) = AWAC(IRH)
      IF (M0(06) .LT. 100.0d0) THEN
         IC = M0(06)
         CALL KMTAB(IC,298.d0,    XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,      XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
         M0(06) = M0(06)*EXP(LN10*(GI0-GII))
      ENDIF

      ! 2H-SO4
      M0(07) = AWSA(IRH)
      IF (M0(07) .LT. 100.0d0) THEN
         IC = 3.0d0*M0(07)
         CALL KMTAB(IC,298.d0,    XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX)
!         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,TEMP,      XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
         M0(07) = M0(07)*EXP(LN10*(GI0-GII))
      ENDIF

      ! H-HSO4
      M0(08) = AWSA(IRH)
CCC      IF (M0(08) .LT. 100.0) THEN     ! These are redundant, because M0(8) is not used
CCC         IC = M0(08)
CCC         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
CCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
CCCCCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX)
CCC         M0(08) = M0(08)*EXP(LN10*(GI0-GII))
CCC      ENDIF

      ! NH4HSO4
      M0(09) = AWAB(IRH)
      IF (M0(09) .LT. 100.0d0) THEN
         IC = M0(09)
         CALL KMTAB(IC,298.d0,    XX,XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX)
!         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
         CALL KMTAB(IC,TEMP,      XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
         M0(09) = M0(09)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NAHSO4
      M0(12) = AWSB(IRH)
      IF (M0(12) .LT. 100.0d0) THEN
         IC = M0(12)
         CALL KMTAB(IC,298.d0,    XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GI0)
!         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
         CALL KMTAB(IC,TEMP,      XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
         M0(12) = M0(12)*EXP(LN10*(GI0-GII))
      ENDIF

      ! (NH4)3H(SO4)2
      M0(13) = AWLC(IRH)
      IF (M0(13) .LT. 100.0d0) THEN
         IC     = 4.0d0*M0(13)
         CALL KMTAB(IC,298.d0,    XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G130   = 0.2d0*(3.0d0*GI0+2.0d0*GII)
!         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         CALL KMTAB(IC,TEMP,      XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G13I   = 0.2d0*(3.0d0*GI0+2.0d0*GII)
!         M0(13) = M0(13)*EXP(LN10*SNGL(G130-G13I))
         M0(13) = M0(13)*EXP(LN10*(G130-G13I))
      ENDIF

      !=================================================================
      ! OTHER INITIALIZATIONS
      !=================================================================
      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = '??'
      SULRATW = 2.D0
      NOFER   = 0
      STKOFL  =.FALSE.

      DO K=1,NERRMX
         ERRSTK(K) =-999
         ERRMSG(K) = 'MESSAGE N/A'
      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT3

!-----------------------------------------------------------------------------

      FUNCTION GETASR( SO4I, RHI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION GETASR
! *** CALCULATES THE LIMITING NH4+/SO4 RATIO OF A SULFATE POOR SYSTEM
!     (i.e. SULFATE RATIO = 2.0) FOR GIVEN SO4 LEVEL AND RH
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: SO4I, RHI

      ! Local variables
      INTEGER            :: IA1,   INDS,  INDR, INDSL
      INTEGER            :: INDSH, IPOSL, IPOSH
      REAL*8             :: GETASR, RAT,  A1, WF

            
      !=================================================================
      ! GETASR begins here!
      !=================================================================

      !### SOLVE USING FULL COMPUTATIONS, NOT LOOK-UP TABLES
      !W(2) = WAER(2)
      !W(3) = WAER(2)*2.0001D0
      !CALL CALCA2
      !SULRATW = MOLAL(3)/WAER(2)
      !CALL INIT1 (WI, RHI, TEMPI) ! Re-initialize COMMON BLOCK

      ! CALCULATE INDICES
      RAT    = SO4I/1.d-9    
      A1     = INT(LOG(RAT))                   ! Magnitude of RAT
      IA1    = INT(RAT/2.5d0/10.0d0**A1)

      INDS   = 4.0d0*A1 + MIN(IA1,4)
      INDS   = MIN(MAX(0, INDS), NSO4S-1) + 1     ! SO4 component of IPOS

      INDR   = INT(99.0d0-RHI*100.0d0) + 1
      INDR   = MIN(MAX(1, INDR), NRHS)            ! RH component of IPOS

      ! GET VALUE AND RETURN
      INDSL  = INDS
      INDSH  = MIN(INDSL+1, NSO4S)
      IPOSL  = (INDSL-1)*NRHS + INDR              ! Low position in array
      IPOSH  = (INDSH-1)*NRHS + INDR              ! High position in array

      WF     = (SO4I-ASSO4(INDSL))/(ASSO4(INDSH)-ASSO4(INDSL) + 1d-7)
      WF     = MIN(MAX(WF, 0.0d0), 1.0d0)

      GETASR = WF*ASRAT(IPOSH) + (1.0d0-WF)*ASRAT(IPOSL)

      ! Return to calling program
      END FUNCTION GETASR

!------------------------------------------------------------------------------

      SUBROUTINE CALCNA
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNA
! *** CALCULATES NITRATES SPECIATION
!
!     NITRIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES, THAT 
!     DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. THE NITRIC
!     ACID DISSOLVED IS CALCULATED FROM THE HNO3(G) -> (H+) + (NO3-) 
!     EQUILIBRIUM, USING THE (H+) FROM THE SULFATES.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      REAL*8 :: KAPA, X, DELT, ALFA, DIAK

      !### comment out for now
      !CHARACTER ERRINF*40

      !=================================================================
      ! CALCNA begins here!
      !=================================================================

      ! CALCULATE HNO3 DISSOLUTION
      X    = W(4) 
      DELT = 0.0d0
      IF (WATER.GT.TINY) THEN
         KAPA = MOLAL(1)
         ALFA = XK4*R*TEMP*(WATER/GAMA(10))**2.0d0
         DIAK = SQRT( (KAPA+ALFA)**2.0d0 + 4.0d0*ALFA*X)
         DELT = 0.5d0*(-(KAPA+ALFA) + DIAK)

         !### Comment out for now
         !IF (DELT/KAPA.GT.0.1d0) THEN
         !   WRITE (ERRINF,'(1PE10.3)') DELT/KAPA*100.0
         !   CALL PUSHERR (0019, ERRINF)    ! WARNING ERROR: NO SOLUTION
         !ENDIF
      ENDIF

      ! CALCULATE HNO3 SPECIATION IN THE GAS PHASE 
      GHNO3    = MAX(X-DELT, 0.0d0)  ! GAS HNO3

      ! CALCULATE HNO3 SPECIATION IN THE LIQUID PHASE
      MOLAL(7) = DELT                ! NO3-
      MOLAL(1) = MOLAL(1) + DELT     ! H+ 

      ! Return to calling program
      END SUBROUTINE CALCNA

!------------------------------------------------------------------------------

      SUBROUTINE CALCNH3
!
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNH3
! *** CALCULATES AMMONIA IN GAS PHASE
!
!     AMMONIA IN THE GAS PHASE IS ASSUMED A MINOR SPECIES, THAT 
!     DOES NOT SIGNIFICANTLY PERTURB THE AEROSOL EQUILIBRIUM. 
!     AMMONIA GAS IS CALCULATED FROM THE NH3(g) + (H+)(l) <==> (NH4+)(l)
!     EQUILIBRIUM, USING (H+), (NH4+) FROM THE AEROSOL SOLUTION.
!
!     THIS IS THE VERSION USED BY THE DIRECT PROBLEM
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
#     include "isoropia.h" 

      ! Local variables 
      REAL*8 :: A1, CHI1, CHI2, BB, CC, DIAK, PSI

      !=================================================================
      ! CALCNH3 begins here!
      !=================================================================

      ! Is there a liquid phase?
      IF ( WATER .LE. TINY ) RETURN

      ! Calculate NH3 sublimation
      A1   =  (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0d0
      CHI1 =  MOLAL(3)
      CHI2 =  MOLAL(1)

      BB   = (CHI2 + ONE/A1)          ! a=1; b!=1; c!=1 
      CC   = -CHI1/A1             
      DIAK =  SQRT(BB*BB - 4.D0*CC)   ! Always > 0
      PSI  =  0.5d0*(-BB + DIAK)        ! One positive root
      PSI  =  MAX(TINY, MIN(PSI,CHI1))! Constrict in acceptible range

      ! Calculate NH3 speciation in the gas phase
      GNH3     = PSI                 ! GAS HNO3

      ! Calculate NH3 affect in the liquid phase
      MOLAL(3) = CHI1 - PSI          ! NH4+
      MOLAL(1) = CHI2 + PSI          ! H+ 
 
      ! Return to calling program
      END SUBROUTINE CALCNH3

!------------------------------------------------------------------------------

      SUBROUTINE CALCNIAQ( NO3I, HI, DELT )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNIAQ
!
!     THIS SUBROUTINE CALCULATES THE HNO3(aq) GENERATED FROM (H,NO3-).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8 :: NO3I, HI, DELT

      ! Local variables
      REAL*8 :: A42, OM1, OM2, BB, CC, DD, DEL1, DEL2

      !=================================================================
      ! CALCNIAQ begins here!
      !=================================================================

      ! equilibrium constants
      A42  = XK42*WATER/(GAMA(10))**2. ! gama(hno3) assumed 1

      ! Find root
      OM1  =  NO3I          
      OM2  =  HI
      BB   = -(OM1+OM2+A42)
      CC   =  OM1*OM2
      DD   =  SQRT(BB*BB-4.D0*CC)

      DEL1 =  0.5D0*(-BB - DD)
      DEL2 =  0.5D0*(-BB + DD)

      ! Get appropriate root.
      IF (DEL1.LT.ZERO .OR. DEL1.GT.HI .OR. DEL1.GT.NO3I) THEN
         DELT = ZERO
      ELSE
         DELT = DEL1
         RETURN
      ENDIF

      IF (DEL2.LT.ZERO .OR. DEL2.GT.NO3I .OR. DEL2.GT.HI) THEN
         DELT = ZERO
      ELSE
         DELT = DEL2
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCNIAQ

!------------------------------------------------------------------------------

      SUBROUTINE CALCNIAQ2( GGNO3, NO3I, HI, NO3AQ )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNIAQ2
!
!     THIS SUBROUTINE CALCULATES THE UNDISSOCIATED HNO3(aq)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8 :: NO3I, NO3AQ, GGNO3, HI 

      ! Local variables
      REAL*8 :: A42, AKW, ALF1, ALF2, ALF3, BB, CC, DEL1

      !=================================================================
      ! CALCNIAQ2 begins here!
      !=================================================================

      ! Equilibrium constants
      A42  = XK42*WATER/(GAMA(10))**2.d0 ! GAMA(HNO3) assumed 1
      AKW  = XKW *RHB*WATER*WATER

      ! Find root
      ALF1  = NO3I - GGNO3
      ALF2  = GGNO3
      ALF3  = HI

      BB    = ALF3 + ALF1 + A42
      CC    = ALF3*ALF1 - A42*ALF2
      DEL1  = 0.5d0*(-BB + SQRT(BB*BB-4.D0*CC))

      ! Correct concentrations
      NO3I  = ALF1 + DEL1
      HI    = ALF3 + DEL1
      IF (HI.LE.TINY) HI = SQRT(AKW)   ! If solution is neutral.
      NO3AQ = ALF2 - DEL1

      ! Return to calling program
      END SUBROUTINE CALCNIAQ2

!------------------------------------------------------------------------------

      SUBROUTINE CALCMR
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCMR
! *** THIS SUBROUTINE CALCULATES:
!     1. ION PAIR CONCENTRATIONS (FROM [MOLAR] ARRAY)
!     2. WATER CONTENT OF LIQUID AEROSOL PHASE (FROM ZSR CORRELATION)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      REAL*8           :: SO4I, HSO4I, AML5, TOTS4, FRNH4, FRNO3, FRCL
      INTEGER          :: K
      CHARACTER(LEN=1) :: SC
      REAL*8      :: CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCMR begins here!
      !=================================================================

      ! CALCULATE ION PAIR CONCENTRATIONS ACCORDING TO SPECIFIC CASE
      SC =SCASE(1:1)                   ! SULRAT & SODRAT case

      !=================================================================
      ! NH4-SO4 SYSTEM ; SULFATE POOR CASE
      !=================================================================
      IF (SC.EQ.'A') THEN      
         ! (NH4)2SO4 - CORRECT FOR SO4 TO HSO4
         MOLALR(4) = MOLAL(5)+MOLAL(6) 

      !=================================================================
      ! NH4-SO4 SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'B') THEN
	 ! CORRECT FOR HSO4 DISSOCIATION 
         SO4I  = MOLAL(5)-MOLAL(1)     
         HSO4I = MOLAL(6)+MOLAL(1)              
         IF (SO4I.LT.HSO4I) THEN                
            MOLALR(13) = SO4I                   ! [LC] = [SO4]       
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)  ! NH4HSO4
         ELSE                                   
            MOLALR(13) = HSO4I                  ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)  ! (NH4)2SO4
         ENDIF

      !=================================================================
      ! NH4-SO4 SYSTEM ; SULFATE RICH CASE ; FREE ACID 
      !=================================================================
      ELSE IF (SC.EQ.'C') THEN
         MOLALR(4) = MOLAL(3)                     ! NH4HSO4
         MOLALR(7) = MAX(W(2)-W(3), ZERO)  ! H2SO4

      !=================================================================
      ! NH4-SO4 SYSTEM ; SULFATE RICH CASE ; FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'D') THEN      
         MOLALR(4) = MOLAL(5) + MOLAL(6)  ! (NH4)2SO4
         AML5      = MOLAL(3)-2.D0*MOLALR(4)     ! "free" NH4
         MOLALR(5) = MAX(MIN(AML5,MOLAL(7)), ZERO)
						! NH4NO3 = MIN("free", NO3)

      !=================================================================
      ! NH4-SO4-NO3 SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'E') THEN      
         SO4I  = MAX(MOLAL(5)-MOLAL(1),ZERO)
						! FROM HSO4 DISSOCIATION 
         HSO4I = MOLAL(6)+MOLAL(1)              
         IF (SO4I.LT.HSO4I) THEN                
            MOLALR(13) = SO4I                     ! [LC] = [SO4] 
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)    ! NH4HSO4
         ELSE                                   
            MOLALR(13) = HSO4I                    ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)    ! (NH4)2SO4
         ENDIF

      !=================================================================
      ! NH4-SO4-NO3 SYSTEM ; SULFATE RICH CASE ; FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'F') THEN      
         MOLALR(4) = MOLAL(3)                   ! NH4HSO4
         MOLALR(7) = MAX(MOLAL(5)+MOLAL(6)-
     &                                 MOLAL(3),ZERO)  ! H2SO4

      !=================================================================
      ! NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE POOR ; SODIUM POOR CASE
      !=================================================================
      ELSE IF (SC.EQ.'G') THEN      
         MOLALR(2) = 0.5d0*MOLAL(2)            ! NA2SO4
         TOTS4     = MOLAL(5)+MOLAL(6)         ! Total SO4
         MOLALR(4) = MAX(TOTS4 - MOLALR(2), ZERO) ! (NH4)2SO4
         FRNH4     = MAX(MOLAL(3) - 2.D0*MOLALR(4), ZERO)
         MOLALR(5) = MIN(MOLAL(7),FRNH4)       ! NH4NO3
         FRNH4     = MAX(FRNH4 - MOLALR(5), ZERO)
         MOLALR(6) = MIN(MOLAL(4), FRNH4)        ! NH4CL

      !=================================================================
      ! NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE POOR ; SODIUM RICH CASE
      ! RETREIVE DISSOLVED SALTS DIRECTLY FROM COMMON BLOCK /SOLUT/
      !=================================================================
      ELSE IF (SC.EQ.'H') THEN      
         MOLALR(1) = PSI7                                  ! NACL 
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(3) = PSI8                                  ! NANO3
         MOLALR(4) = ZERO                                  ! (NH4)2SO4
         FRNO3     = MAX(MOLAL(7) - MOLALR(3), ZERO) ! "FREE" NO3
         FRCL      = MAX(MOLAL(4) - MOLALR(1), ZERO) ! "FREE" CL
         MOLALR(5) = MIN(MOLAL(3),FRNO3)             ! NH4NO3
         FRNH4     = MAX(MOLAL(3) - MOLALR(5), ZERO)  ! "FREE" NH3
         MOLALR(6) = MIN(FRCL, FRNH4)                      ! NH4CL

      !=================================================================
      ! NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
      ! RETREIVE DISSOLVED SALTS DIRECTLY FROM COMMON BLOCK /SOLUT/
      !=================================================================
      ELSE IF (SC.EQ.'I') THEN      
         MOLALR(04) = PSI5                                 ! (NH4)2SO4
         MOLALR(02) = PSI4                                 ! NA2SO4
         MOLALR(09) = PSI1                                 ! NH4HSO4
         MOLALR(12) = PSI3                                 ! NAHSO4
         MOLALR(13) = PSI2                                 ! LC

      !=================================================================
      ! NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE RICH CASE ; FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'J') THEN      
         MOLALR(09) = MOLAL(3)                             ! NH4HSO4
         MOLALR(12) = MOLAL(2)                             ! NAHSO4
         MOLALR(07) = MOLAL(5)+MOLAL(6)-
     &                      MOLAL(3)-MOLAL(2)  ! H2SO4
         MOLALR(07) = MAX(MOLALR(07),ZERO)

      !=================================================================
      ! ----- REVERSE PROBLEMS ----- 
      !
      ! NH4-SO4-NO3 SYSTEM ; SULFATE POOR CASE
      !=================================================================
      ELSE IF (SC.EQ.'N') THEN      
         MOLALR(4) = MOLAL(5) + MOLAL(6) ! (NH4)2SO4
         AML5      = WAER(3)-2.D0*MOLALR(4)    ! "free" NH4
         MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO) ! NH4NO3 = MIN("free", NO3)

      !=================================================================
      ! NH4-SO4-NO3-NA-CL SYSTEM ; SULFATE POOR, SODIUM POOR CASE
      !=================================================================
      ELSE IF (SC.EQ.'Q') THEN      
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(4) = PSI6                                  ! (NH4)2SO4
         MOLALR(5) = PSI5                                  ! NH4NO3
         MOLALR(6) = PSI4                                  ! NH4CL

      !=================================================================
      ! NH4-SO4-NO3-NA-CL SYSTEM ; SULFATE POOR, SODIUM RICH CASE
      !=================================================================
      ELSE IF (SC.EQ.'R') THEN      
         MOLALR(1) = PSI3                                  ! NACL 
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(3) = PSI2                                  ! NANO3
         MOLALR(4) = ZERO                                  ! (NH4)2SO4
         MOLALR(5) = PSI5                                  ! NH4NO3
         MOLALR(6) = PSI4                                  ! NH4CL

      !=================================================================
      ! UNKNOWN CASE
      !=================================================================
      ELSE
         CALL PUSHERR (1001, ' ') ! FATAL ERROR: CASE NOT SUPPORTED 
      ENDIF

      !=================================================================
      ! CALCULATE WATER CONTENT ; ZSR CORRELATION
      !=================================================================
      WATER = ZERO
      DO K = 1, NPAIR
         WATER = WATER + MOLALR(K)/M0(K)
      ENDDO
      WATER = MAX(WATER, TINY)

      ! Return to calling program
      END SUBROUTINE CALCMR

!------------------------------------------------------------------------------

      SUBROUTINE CALCMDRH( RHI, RHDRY, RHLIQ, 
     &                     DRYCASE, LIQCASE )

!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCMDRH
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE 'DRY' SOLUTION (SUBROUTINE DRYCASE) AND THE
!     'SATURATED LIQUID' SOLUTION (SUBROUTINE LIQCASE).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8  :: RHI, RHDRY, RHLIQ
      EXTERNAL   DRYCASE, LIQCASE

      ! Local variables
      REAL*8  :: WF, ONEMWF
      REAL*8  :: CNH42SO, CNH4HSO, CLCO, CNH4N3O, CNH4CLO, CNA2SO 
      REAL*8  :: CNAHSO, CNANO, CNACLO, GNH3O, GHNO3O, GHCLO   
      REAL*8  :: DAMSUL, DSOSUL, DAMBIS, DSOBIS, DLC, DAMNIT
      REAL*8  :: DAMCHL, DSONIT, DSOCHL, DAMG, DHAG, DNAG
      INTEGER :: K

      !=================================================================
      ! CALCMDRH begins here!
      !=================================================================

      ! Find weight factor
      IF (WFTYP.EQ.0) THEN
         WF = ONE
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
      ENDIF
      ONEMWF  = ONE - WF

      !=================================================================
      ! FIND FIRST SECTION ; DRY ONE
      !=================================================================

      ! Find solution at mdrh by weighting dry & liquid solutions.
      CALL DRYCASE
      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  ! DRY AEROSOL

      CNH42SO = CNH42S4                  ! FIRST (DRY) SOLUTION
      CNH4HSO = CNH4HS4
      CLCO    = CLC 
      CNH4N3O = CNH4NO3
      CNH4CLO = CNH4CL
      CNA2SO  = CNA2SO4
      CNAHSO  = CNAHSO4
      CNANO   = CNANO3
      CNACLO  = CNACL
      GNH3O   = GNH3
      GHNO3O  = GHNO3
      GHCLO   = GHCL

      !=================================================================
      ! FIND SECOND SECTION ; DRY & LIQUID
      !=================================================================
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CLC     = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CNAHSO4 = ZERO
      CNANO3  = ZERO
      CNACL   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
      CALL LIQCASE                   ! SECOND (LIQUID) SOLUTION

      !=================================================================
      ! ADJUST THINGS FOR THE CASE THAT THE 
      ! LIQUID SUB PREDICTS DRY AEROSOL
      !=================================================================
      IF ( WATER .LE. TINY ) THEN

         ! Aqueous phase
         DO K = 1, NIONS
            MOLAL(K)= ZERO           
         ENDDO

         WATER   = ZERO

         ! Solid phase
         CNH42S4 = CNH42SO           
         CNA2SO4 = CNA2SO
         CNAHSO4 = CNAHSO
         CNH4HS4 = CNH4HSO
         CLC     = CLCO
         CNH4NO3 = CNH4N3O
         CNANO3  = CNANO
         CNACL   = CNACLO                                                  
         CNH4CL  = CNH4CLO 

         ! Gas Phase
         GNH3    = GNH3O             
         GHNO3   = GHNO3O
         GHCL    = GHCLO

         GOTO 200
      ENDIF

      !=================================================================
      ! FIND SALT DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
      !=================================================================
      DAMSUL  = CNH42SO - CNH42S4
      DSOSUL  = CNA2SO  - CNA2SO4
      DAMBIS  = CNH4HSO - CNH4HS4
      DSOBIS  = CNAHSO  - CNAHSO4
      DLC     = CLCO    - CLC
      DAMNIT  = CNH4N3O - CNH4NO3
      DAMCHL  = CNH4CLO - CNH4CL
      DSONIT  = CNANO   - CNANO3
      DSOCHL  = CNACLO  - CNACL

      !=================================================================
      ! FIND GAS DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
      !=================================================================
      DAMG    = GNH3O   - GNH3 
      DHAG    = GHCLO   - GHCL
      DNAG    = GHNO3O  - GHNO3

      !=================================================================
      ! FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
      !=================================================================

      ! LIQUID
      MOLAL(1)= ONEMWF*MOLAL(1)                                 ! H+
      MOLAL(2)= ONEMWF*(2.D0*DSOSUL + DSOBIS + DSONIT + DSOCHL) ! NA+
      MOLAL(3)= ONEMWF*(2.D0*DAMSUL + DAMG   + DAMBIS + DAMCHL +
     &                  3.D0*DLC    + DAMNIT )                  ! NH4+
      MOLAL(4)= ONEMWF*( DAMCHL + DSOCHL + DHAG)            ! CL-
      MOLAL(5)= ONEMWF*( DAMSUL + DSOSUL + DLC)             ! SO4--
      MOLAL(6)= ONEMWF*( MOLAL(6) + DSOBIS + DAMBIS + DLC) ! HSO4-
      MOLAL(7)= ONEMWF*(     DAMNIT + DSONIT + DNAG)            ! NO3-
      WATER   = ONEMWF*WATER

      ! SOLID
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
      CNACL   = WF*CNACLO  + ONEMWF*CNACL
      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL

      ! GAS
      GNH3    = WF*GNH3O   + ONEMWF*GNH3
      GHNO3   = WF*GHNO3O  + ONEMWF*GHNO3
      GHCL    = WF*GHCLO   + ONEMWF*GHCL

      ! Return to calling program
 200  CONTINUE
      END SUBROUTINE CALCMDRH

!------------------------------------------------------------------------------

      SUBROUTINE CALCMDRP( RHI, RHDRY, RHLIQ, 
     &                     DRYCASE, LIQCASE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCMDRP
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE 'DRY' SOLUTION (SUBROUTINE DRYCASE) AND THE
!     'SATURATED LIQUID' SOLUTION (SUBROUTINE LIQCASE).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8  :: RHI, RHDRY, RHLIQ
      EXTERNAL   DRYCASE, LIQCASE

      ! Local variables
      REAL*8  :: WF, ONEMWF
      REAL*8  :: CNH42SO, CNH4HSO, CLCO, CNH4N3O, CNH4CLO, CNA2SO 
      REAL*8  :: CNAHSO, CNANO, CNACLO
      REAL*8  :: DAMBIS, DSOBIS, DLC, A8
      REAL*8  :: HIEQ, HIEN, A2, A3, A4
      INTEGER :: K

      !=================================================================
      ! CALCMDRP begins here!
      !=================================================================

      ! FIND WEIGHT FACTOR
      IF (WFTYP.EQ.0) THEN
         WF = ONE
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
      ENDIF
      ONEMWF  = ONE - WF

      !=================================================================
      ! FIND FIRST SECTION ; DRY ONE 
      !=================================================================
      CALL DRYCASE
      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  ! DRY AEROSOL

      CNH42SO = CNH42S4              ! FIRST (DRY) SOLUTION
      CNH4HSO = CNH4HS4
      CLCO    = CLC 
      CNH4N3O = CNH4NO3
      CNH4CLO = CNH4CL
      CNA2SO  = CNA2SO4
      CNAHSO  = CNAHSO4
      CNANO   = CNANO3
      CNACLO  = CNACL

      !=================================================================
      ! FIND SECOND SECTION ; DRY & LIQUID 
      !=================================================================
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CLC     = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CNAHSO4 = ZERO
      CNANO3  = ZERO
      CNACL   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
      CALL LIQCASE                   ! SECOND (LIQUID) SOLUTION
 
      !=================================================================
      ! ADJUST THINGS FOR THE CASE THAT THE 
      ! LIQUID SUB PREDICTS DRY AEROSOL
      !=================================================================
      IF ( WATER .LE. TINY ) THEN
         WATER = ZERO

         DO K = 1, NIONS
            MOLAL(K)= ZERO
         ENDDO

         CALL DRYCASE
         GOTO 200
      ENDIF

      !=================================================================
      ! FIND SALT DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
      !=================================================================
      DAMBIS  = CNH4HSO - CNH4HS4
      DSOBIS  = CNAHSO  - CNAHSO4
      DLC     = CLCO    - CLC

      !=================================================================
      ! FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
      !=================================================================

      ! SOLID
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
      CNACL   = WF*CNACLO  + ONEMWF*CNACL
      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL

      ! LIQUID
      WATER   = ONEMWF*WATER

      MOLAL(2)= WAER(1) - 2.D0*CNA2SO4 - CNAHSO4 - CNANO3 -     
     &                         CNACL                            ! NA+
      MOLAL(3)= WAER(3) - 2.D0*CNH42S4 - CNH4HS4 - CNH4CL - 
     &                    3.D0*CLC     - CNH4NO3                ! NH4+
      MOLAL(4)= WAER(5) - CNACL - CNH4CL            ! CL-
      MOLAL(7)= WAER(4) - CNANO3 - CNH4NO3          ! NO3-
      MOLAL(6)= ONEMWF*(MOLAL(6) + DSOBIS + DAMBIS + DLC) ! HSO4-
      MOLAL(5)= WAER(2) - MOLAL(6) - CLC - 
     &                    CNH42S4 - CNA2SO4    ! SO4--

      A8      = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
      IF (MOLAL(5).LE.TINY) THEN
         HIEQ = SQRT(XKW *RHB*WATER*WATER)  ! Neutral solution
      ELSE
         HIEQ = A8*MOLAL(6)/MOLAL(5)          
      ENDIF
      HIEN    = MOLAL(4) + MOLAL(7) + 
     &          MOLAL(6) + 2.D0*MOLAL(5) -
     &          MOLAL(2) - MOLAL(3)
      MOLAL(1)= MAX (HIEQ, HIEN)                  ! H+

      ! GAS (ACTIVITY COEFS FROM LIQUID SOLUTION)
      A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

      GNH3    = MOLAL(3)/MAX(MOLAL(1),TINY)/A2
      GHNO3   = MOLAL(1)*MOLAL(7)/A3
      GHCL    = MOLAL(1)*MOLAL(4)/A4

      ! Return to calling program
 200  CONTINUE
      END SUBROUTINE CALCMDRP

!------------------------------------------------------------------------------

      SUBROUTINE CALCHS4( HI, SO4I, HSO4I, DELTA )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCHS4
! *** THIS SUBROUTINE CALCULATES THE HSO4 GENERATED FROM (H,SO4).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8 :: HI, SO4I, HSO4I, DELTA

      ! Local variables
      REAL*8 :: A8, BB, CC, DD, SQDD, DELTA1, DELTA2
      !CHARACTER ERRINF*40

      !=================================================================
      ! CALCHS4 begins here!
      !=================================================================

      ! IF TOO LITTLE WATER, DONT SOLVE
      IF (WATER.LE.1d1*TINY) THEN
         DELTA = ZERO 
         RETURN
      ENDIF

      !=================================================================
      ! CALCULATE HSO4 SPECIATION
      !=================================================================
      A8 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

      BB =-(HI + SO4I + A8)
      CC = HI*SO4I - HSO4I*A8
      DD = BB*BB - 4.D0*CC

      IF ( DD .GE. ZERO ) THEN
         SQDD   = SQRT(DD)
         DELTA1 = 0.5d0*(-BB + SQDD)
         DELTA2 = 0.5d0*(-BB - SQDD)
         IF (HSO4I.LE.TINY) THEN
            DELTA = DELTA2
         ELSEIF( HI*SO4I .GE. A8*HSO4I ) THEN
            DELTA = DELTA2
         ELSEIF( HI*SO4I .LT. A8*HSO4I ) THEN
            DELTA = DELTA1
         ELSE
            DELTA = ZERO
         ENDIF
      ELSE
         DELTA  = ZERO
      ENDIF

      !-----------------------------------------------------------------------
      ! Comment out for now
      !!=================================================================
      !!COMPARE DELTA TO TOTAL H+ ; ESTIMATE EFFECT OF HSO4 
      !!=================================================================
      !
      !HYD = MAX(HI, MOLAL(1))
      !IF (HYD.GT.TINY) THEN
      !   IF (DELTA/HYD.GT.0.1d0) THEN
      !      WRITE (ERRINF,'(1PE10.3)') DELTA/HYD*100.0
      !      CALL PUSHERR (0020, ERRINF)
      !   ENDIF
      !ENDIF
      !-----------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE CALCHS4

!------------------------------------------------------------------------------

      SUBROUTINE CALCPH( GG, HI, OHI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCPH
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8 :: GG, HI, OHI

      ! Local variables
      REAL*8 :: AKW, CN, BB, DD, CC

      AKW  = XKW *RHB*WATER*WATER
      CN   = SQRT(AKW)

      !=================================================================
      ! CALCPH begins here!
      !=================================================================

      ! GG = (negative charge) - (positive charge)
      ! H+ in excess
      IF ( GG .GT. TINY ) THEN                        
         BB  = -GG
         CC  = -AKW
         DD  =  BB*BB - 4.D0*CC
         HI  =  MAX(0.5D0*(-BB + SQRT(DD)),CN)
         OHI =  AKW/HI

      ! OH- in excess
      ELSE                                        
         BB  =  GG
         CC  = -AKW
         DD  =  BB*BB - 4.D0*CC
         OHI =  MAX(0.5D0*(-BB + SQRT(DD)),CN)
         HI  =  AKW/OHI

      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCPH

!------------------------------------------------------------------------------

      SUBROUTINE CALCACT
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCACT
! *** CALCULATES MULTI-COMPONENET ACTIVITY COEFFICIENTS FROM BROMLEYS
!     METHOD. THE BINARY ACTIVITY COEFFICIENTS ARE CALCULATED BY 
!     KUSIK-MEISNER RELATION (SUBROUTINE KMTAB or SUBROUTINE KMFUL). 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      REAL*8 , PARAMETER :: URF=0.5d0
!      REAL*8             :: EX10
      REAL*8             :: G0(3,4),ZPL,ZMI,AGAMA,SION,H,CH,F1(3),F2(4)
!      REAL*8             :: TEMP
!      INTEGER            :: I, J
      INTEGER            :: K, M
      REAL*8             :: MPL, XIJ, YJI, ERROU, ERRIN

      !=================================================================
      ! CALCACT begins here!
      !=================================================================
      ! This is an F77-style statement function, which is obsolete
      ! under F90.  We moved this into FUNCTION G below.
      ! (bec, bmy, 3/3/05)
      !G(I,J)= (F1(I)/Z(I) + F2(J)/Z(J+3)) / (Z(I)+Z(J+3)) - H

      !=================================================================
      ! SAVE ACTIVITIES IN OLD ARRAY
      !=================================================================
      IF ( FRST ) THEN               

         ! Outer loop
         DO K = 1, NPAIR
            GAMOU(K) = GAMA(K)
         ENDDO
      ENDIF

      ! Inner loop
      DO K = 1, NPAIR              
         GAMIN(K) = GAMA(K)
      ENDDO

      !=================================================================
      ! CALCULATE IONIC ACTIVITY OF SOLUTION
      !=================================================================
      IONIC = 0.0d0

      DO K = 1, NIONS
         IONIC = IONIC + MOLAL(K)*Z(K)*Z(K)
      ENDDO

      IONIC = MAX(MIN(0.5d0*IONIC/WATER,20.d0), TINY)

      !=================================================================
      ! CALCULATE BINARY ACTIVITY COEFFICIENTS
      !
      ! G0(1,1)=G11; G0(1,2)=G07; G0(1,3)=G08; G0(1,4)=G10; 
      ! G0(2,1)=G01; G0(2,2)=G02; G0(2,3)=G12; G0(2,4)=G03;
      ! G0(3,1)=G06; G0(3,2)=G04; G0(3,3)=G09; G0(3,4)=G05
      !=================================================================
      IF ( IACALC .EQ. 0 ) THEN              
         ! K.M.; FULL
!         CALL KMFUL (IONIC, REAL(TEMP),G0(2,1),G0(2,2),G0(2,4),
!     &               G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),
!     &               G0(1,4),G0(1,1),G0(2,3))
         CALL KMFUL (IONIC, TEMP,G0(2,1),G0(2,2),G0(2,4),
     &               G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),
     &               G0(1,4),G0(1,1),G0(2,3))
      ELSE                               
         ! K.M.; TABULATED
!         CALL KMTAB (IONIC, REAL(TEMP),G0(2,1),G0(2,2),G0(2,4),
!     &               G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),
!     &               G0(1,4),G0(1,1),G0(2,3))
         CALL KMTAB (IONIC, TEMP,G0(2,1),G0(2,2),G0(2,4),
     &               G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),
     &               G0(1,4),G0(1,1),G0(2,3))
      ENDIF

      !=================================================================
      ! CALCULATE MULTICOMPONENT ACTIVITY COEFFICIENTS
      !=================================================================
      AGAMA = 0.511d0*(298.0d0/TEMP)**1.5d0    ! Debye Huckel const. at T
      SION  = SQRT(IONIC)
      H     = AGAMA*SION/(1+SION)

      DO K = 1,3
         F1(K)=0.0d0
         F2(K)=0.0d0
      ENDDO
         
      F2(4)=0.0d0

      DO K = 1,3
         ZPL = Z(K)
         MPL = MOLAL(K)/WATER
         DO M = 1,4
            ZMI   = Z(M+3)
            CH    = 0.25d0*(ZPL+ZMI)*(ZPL+ZMI)/IONIC
            XIJ   = CH*MPL
            YJI   = CH*MOLAL(M+3)/WATER
!            F1(I) = F1(I) + SNGL(YJI*(G0(I,J) + ZPL*ZMI*H))
            F1(K) = F1(K) + (YJI*(G0(K,M) + ZPL*ZMI*H))
!            F2(J) = F2(J) + SNGL(XIJ*(G0(I,J) + ZPL*ZMI*H))
            F2(M) = F2(M) + (XIJ*(G0(K,M) + ZPL*ZMI*H))
         ENDDO
      ENDDO

      !=================================================================
      ! LOG10 OF ACTIVITY COEFFICIENTS
      !=================================================================
      GAMA(01) = G(H,2,1,F1(2),F2(1))*ZZ(01)                     ! NACL
      GAMA(02) = G(H,2,2,F1(2),F2(2))*ZZ(02)                     ! NA2SO4
      GAMA(03) = G(H,2,4,F1(2),F2(4))*ZZ(03)                     ! NANO3
      GAMA(04) = G(H,3,2,F1(3),F2(2))*ZZ(04)                     ! (NH4)2SO4
      GAMA(05) = G(H,3,4,F1(3),F2(4))*ZZ(05)                     ! NH4NO3
      GAMA(06) = G(H,3,1,F1(3),F2(1))*ZZ(06)                     ! NH4CL
      GAMA(07) = G(H,1,2,F1(1),F2(2))*ZZ(07)                     ! 2H-SO4
      GAMA(08) = G(H,1,3,F1(1),F2(3))*ZZ(08)                     ! H-HSO4
      GAMA(09) = G(H,3,3,F1(3),F2(3))*ZZ(09)                     ! NH4HSO4
      GAMA(10) = G(H,1,4,F1(1),F2(4))*ZZ(10)                     ! HNO3
      GAMA(11) = G(H,1,1,F1(1),F2(1))*ZZ(11)                     ! HCL
      GAMA(12) = G(H,2,3,F1(2),F2(3))*ZZ(12)                     ! NAHSO4
      GAMA(13) = 0.2d0*(3.d0*GAMA(04)+2.d0*GAMA(09))  ! LC ; SCAPE
      !GAMA(13) = 0.50*(GAMA(04)+GAMA(09))          ! LC ; SEQUILIB
      !GAMA(13) = 0.25*(3.0*GAMA(04)+GAMA(07))      ! LC ; AIM

      !=================================================================
      ! CONVERT LOG (GAMA) COEFFICIENTS TO GAMA
      !=================================================================
      DO K = 1, NPAIR
         !GAMA(I)=MAX(-5.0d0, MIN(GAMA(I),5.0d0) )   ! F77 LIBRARY ROUTINE
         !GAMA(I)=10.0**GAMA(I)
!         GAMA(I) = EX10(SNGL(GAMA(I)), 5.d0)          ! CUTOFF SET TO [-5,5]
         GAMA(K) = EX10(GAMA(K), 5.d0)          ! CUTOFF SET TO [-5,5]
         GAMA(K) = GAMIN(K)*(1.d0-URF) + URF*GAMA(K)  ! Under-relax GAMA's
      ENDDO

      !=================================================================
      ! SETUP ACTIVITY CALCULATION FLAGS
      !=================================================================

      ! OUTER CALCULATION LOOP ; ONLY IF FRST=.TRUE.
      IF ( FRST ) THEN          

         ! Convergence criterion
         ERROU = ZERO                    

         DO K = 1,NPAIR
            ERROU = MAX(ERROU, ABS((GAMOU(K)-GAMA(K))/GAMOU(K)))
         ENDDO

         ! Setup flags
         CALAOU = ERROU .GE. EPSACT      
         FRST   =.FALSE.
      ENDIF

      ! INNER CALCULATION LOOP ; ALWAYS
      ! Convergence criterion
      ERRIN = ZERO                       

      DO K = 1, NPAIR
         ERRIN = MAX (ERRIN, ABS((GAMIN(K)-GAMA(K))/GAMIN(K)))
      ENDDO

      CALAIN = ERRIN .GE. EPSACT

      ! Increment ACTIVITY call counter
      ICLACT = ICLACT + 1                

      ! Return to calling program
      END SUBROUTINE CALCACT

!------------------------------------------------------------------------------

      FUNCTION G( H, I, J, F1, F2 ) RESULT( VALUE )

#     include "isoropia.h" 

      INTEGER, INTENT(IN) :: I, J
      REAL*8,  INTENT(IN) :: H, F1, F2
      REAL*8              :: VALUE

      !=================================================================
      ! CALCACT begins here!
      !=================================================================
      VALUE = (F1/Z(I) + F2/Z(J+3)) / (Z(I)+Z(J+3)) - H 

      ! Return to calling program
      END FUNCTION G

!------------------------------------------------------------------------------

      SUBROUTINE RSTGAM
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE RSTGAM
! *** RESETS ACTIVITY COEFFICIENT ARRAYS TO DEFAULT VALUE OF 0.1
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER :: K

      DO K = 1, NPAIR
         GAMA(K) = 0.1
      ENDDO

      ! Return to calling program
      END SUBROUTINE RSTGAM

!------------------------------------------------------------------------------

      SUBROUTINE KMFUL( IONIC, TEMP, G01, G02, G03, G04, G05, 
     &                  G06,   G07,  G08, G09, G10, G11, G12 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE KMFUL
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments 
      REAL*8 :: IONIC, TEMP
      REAL*8 :: G01, G02, G03, G04, G05, G06
      REAL*8 :: G07, G08, G09, G10, G11, G12

      ! Local variables
      REAL*8 :: Z01, Z02, Z03, Z04, Z05, Z06
      REAL*8 :: Z07, Z08, Z09, Z10, Z11, SION
      REAL*8 :: TI,  CF1, CF2

      !=================================================================
      ! KMFUL begins here!
      !=================================================================
      
      ! Initialize
      Z01  = 1.d0
      Z02  = 2.d0
      Z03  = 1.d0
      Z04  = 2.d0
      Z05  = 1.d0
      Z06  = 1.d0
      Z07  = 2.d0
      Z08  = 1.d0
      Z10  = 1.d0
      Z11  = 1.d0
      SION = SQRT(IONIC)

      !=================================================================
      ! Coefficients at 25 oC
      !=================================================================
      CALL MKBI( 2.230d0, IONIC, SION, Z01, G01 )
      CALL MKBI( -0.19d0, IONIC, SION, Z02, G02 )
      CALL MKBI( -0.39d0, IONIC, SION, Z03, G03 )
      CALL MKBI( -0.25d0, IONIC, SION, Z04, G04 )
      CALL MKBI( -1.15d0, IONIC, SION, Z05, G05 )
      CALL MKBI( 0.820d0, IONIC, SION, Z06, G06 )
      CALL MKBI( -.100d0, IONIC, SION, Z07, G07 )
      CALL MKBI( 8.000d0, IONIC, SION, Z08, G08 )
      CALL MKBI( 2.600d0, IONIC, SION, Z10, G10 )
      CALL MKBI( 6.000d0, IONIC, SION, Z11, G11 )

      !=================================================================
      ! Correct for T other than 298 K
      !=================================================================
      TI  = TEMP-273.0d0
      IF (ABS(TI) .GT. 1.d0) THEN
         CF1 = 1.125d0-0.005d0*TI
         CF2 = (0.125d0-0.005d0*TI)*(0.039d0*IONIC**0.92d0-0.41d0*SION/
     &         (1.d0+SION))
         G01 = CF1*G01 - CF2*Z01
         G02 = CF1*G02 - CF2*Z02
         G03 = CF1*G03 - CF2*Z03
         G04 = CF1*G04 - CF2*Z04
         G05 = CF1*G05 - CF2*Z05
         G06 = CF1*G06 - CF2*Z06
         G07 = CF1*G07 - CF2*Z07
         G08 = CF1*G08 - CF2*Z08
         G10 = CF1*G10 - CF2*Z10
         G11 = CF1*G11 - CF2*Z11
      ENDIF

      G09 = G06 + G08 - G11
      G12 = G01 + G08 - G11

      ! Return to calling program
      END SUBROUTINE KMFUL

!------------------------------------------------------------------------------

      SUBROUTINE MKBI( Q, IONIC, SION, ZIP, BI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE MKBI
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!
!******************************************************************************
! 
      ! Arguments
      REAL*8 :: IONIC, Q, SION, ZIP, BI

      ! Local variables
      REAL*8 :: B, C, XX

      !=================================================================
      ! MKBI begins here!
      !=================================================================
      B = .75d0-.065d0*Q
      C = 1.d0

      IF (IONIC.LT.6.d0) C=1.d0+.055d0*Q*EXP(-.023d0*IONIC*IONIC*IONIC)

      XX = -0.5107d0*SION/(1.d0+C*SION)
      BI = (1.d0+B*(1.d0+.1d0*IONIC)**Q-B)
      BI = ZIP*LOG(BI) + ZIP*XX

      ! Return to calling program
      END SUBROUTINE MKBI

!------------------------------------------------------------------------------

      SUBROUTINE KMTAB( IN,  TEMP, G01, G02, G03, G04, G05,
     &                  G06, G07,  G08, G09, G10, G11, G12 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE KMTAB
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!     THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!     LOOKUP TABLES. THE IONIC ACTIVITY 'IONIC' IS INPUT, AND THE ARRAY
!     'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: IN, Temp
      REAL*8, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*8, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IND

      !=================================================================
      ! KMC_TAB begins here!
      !=================================================================

      ! Find temperature range
      IND = NINT((TEMP-198.d0)/25.d0) + 1
      IND = MIN(MAX(IND,1),6)

      !=================================================================
      ! Call appropriate routine
      !=================================================================
      IF ( IND .EQ. 1 ) THEN
         CALL KM198(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSEIF ( IND .EQ. 2 ) THEN
         CALL KM223(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSEIF ( IND .EQ. 3 ) THEN
         CALL KM248(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSEIF ( IND .EQ. 4 ) THEN
         CALL KM273(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSEIF ( IND .EQ. 5 ) THEN
         CALL KM298(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSE
         CALL KM323(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ENDIF

      ! Return to calling program
      END SUBROUTINE KMTAB
 
!------------------------------------------------------------------------------

      SUBROUTINE READ_KMC
!
!******************************************************************************
!  Subroutine READ_KMC reads data from binary files in order to initialize
!  parameters for the ISORROPIA routines.  This is necessary since the 
!  original code uses big BLOCK DATA statements which don't compile well on 
!  SGI. (rjp, bdf, bmy, 9/23/02, 7/28/05) 
!
!  NOTES:
!  (1 ) Now read files from "sulfate_sim_200508/isorropia".  Also remove
!        reference to obsolete "CMN_SETUP" (bmy, 7/28/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR

      ! Local variables
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_KMC begins here!
      !=================================================================

      ! Initialize arrays 
      !CALL INIT_KMC

      ! Read data at 198 K
      FILENAME = TRIM( DATA_DIR ) // 
     &           'sulfate_sim_200508/isorropia/kmc198.bin'
      CALL READ_BINARY( FILENAME, BNC198 )

      ! Read data at 223 K
      FILENAME = TRIM( DATA_DIR ) // 
     &           'sulfate_sim_200508/isorropia/kmc223.bin'
      CALL READ_BINARY( FILENAME, BNC223 ) 

      ! Read data at 248 K
      FILENAME = TRIM( DATA_DIR ) // 
     &           'sulfate_sim_200508/isorropia/kmc248.bin'
      CALL READ_BINARY( FILENAME, BNC248 )

      ! Read data at 273 K
      FILENAME = TRIM( DATA_DIR ) // 
     &           'sulfate_sim_200508/isorropia/kmc273.bin'
      CALL READ_BINARY( FILENAME, BNC273 )

      ! Read data at 298 K
      FILENAME = TRIM( DATA_DIR ) // 
     &           'sulfate_sim_200508/isorropia/kmc298.bin'
      CALL READ_BINARY( FILENAME, BNC298 )

      ! Read data at 323 K
      FILENAME = TRIM( DATA_DIR ) // 
     &           'sulfate_sim_200508/isorropia/kmc323.bin'
      CALL READ_BINARY( FILENAME, BNC323 )

      ! Return to calling program
      END SUBROUTINE READ_KMC

!------------------------------------------------------------------------------
      
      SUBROUTINE READ_BINARY( FILENAME, ARRAY )
!
!******************************************************************************
!  Subroutine READ_BINARY reads a binary file from disk. (rjp, bmy, 9/23/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHAR*255) : Name of the binary file to be read

!  Arguments as Output:
!  ============================================================================
!  (2 ) ARRAY    (REAL*4  ) : Data array
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_FILE, IOERROR
      
      ! Arguments
      CHARACTER(LEN=255), INTENT(IN)  :: FILENAME
      REAL*4,             INTENT(OUT) :: ARRAY(IMAX,JMAX)

      ! Local variables
      INTEGER                         :: I, J, IOS
      
      !=================================================================
      ! READ_BINARY begins here!
      !=================================================================

      ! Open the file
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD',
     &               FORM='UNFORMATTED',    IOSTAT=IOS )

      ! Error check
      IF ( IU_FILE /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_kmc:1' )
      
      ! Read data
      DO J = 1, JMAX
         READ( IU_FILE, IOSTAT=IOS ) ( ARRAY(I,J), I=1,IMAX )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_kmc:2' )       
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_BINARY

!-----------------------------------------------------------------------------

      SUBROUTINE KM198( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!*****************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM198
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
! 
!      TEMPERATURE IS 198K
! 
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!*****************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: IN
      REAL*8, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*8, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IPOS

      !=================================================================
      ! KMC323 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.3d2) THEN
         IPOS = NINT( 0.2d2*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.2d1*IN- 0.3d2)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC198(IPOS,1 )
      G02 = BNC198(IPOS,2 )
      G03 = BNC198(IPOS,3 )
      G04 = BNC198(IPOS,4 )
      G05 = BNC198(IPOS,5 )
      G06 = BNC198(IPOS,6 )
      G07 = BNC198(IPOS,7 )
      G08 = BNC198(IPOS,8 )
      G09 = BNC198(IPOS,9 )
      G10 = BNC198(IPOS,10)
      G11 = BNC198(IPOS,11)
      G12 = BNC198(IPOS,12)

      ! Return to calling program
      END SUBROUTINE KM198

!------------------------------------------------------------------------------

      SUBROUTINE KM223( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!******************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM223
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
! 
!      TEMPERATURE IS 223K
! 
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
! 
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: IN
      REAL*8, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*8, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variable
      INTEGER             :: IPOS

      !=================================================================
      ! KMC223 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.3d2) THEN
         IPOS = NINT( 0.2d2*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.2d1*IN- 0.3d2)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC223(IPOS,1)
      G02 = BNC223(IPOS,2)
      G03 = BNC223(IPOS,3)
      G04 = BNC223(IPOS,4)
      G05 = BNC223(IPOS,5)
      G06 = BNC223(IPOS,6)
      G07 = BNC223(IPOS,7)
      G08 = BNC223(IPOS,8)
      G09 = BNC223(IPOS,9)
      G10 = BNC223(IPOS,10)
      G11 = BNC223(IPOS,11)
      G12 = BNC223(IPOS,12)

      ! Return point ; End of subroutine
      END SUBROUTINE KM223

!------------------------------------------------------------------------------

      SUBROUTINE KM248( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!******************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM248
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!  
!      TEMPERATURE IS 248K
!  
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!******************************************************************************
!
      REAL*8, INTENT(IN)  :: IN
      REAL*8, INTENT(OUT) :: G01, G02, G03, G04, G05, G06 
      REAL*8, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER :: IPOS

      !=================================================================
      ! KMC248 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.3d2) THEN
         IPOS = NINT( 0.2d2*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.2d1*IN- 0.3d2)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC248(IPOS,1)
      G02 = BNC248(IPOS,2)
      G03 = BNC248(IPOS,3)
      G04 = BNC248(IPOS,4)
      G05 = BNC248(IPOS,5)
      G06 = BNC248(IPOS,6)
      G07 = BNC248(IPOS,7)
      G08 = BNC248(IPOS,8)
      G09 = BNC248(IPOS,9)
      G10 = BNC248(IPOS,10)
      G11 = BNC248(IPOS,11)
      G12 = BNC248(IPOS,12)

      ! Return point ; End of subroutine
      END SUBROUTINE KM248
 
!------------------------------------------------------------------------------

      SUBROUTINE KM273( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!*****************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM273
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!  
!      TEMPERATURE IS 273K
!
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!*****************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: IN
      REAL*8, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*8, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IPOS

      !=================================================================
      ! KMC273 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.3d2) THEN
         IPOS = NINT( 0.2d2*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.2d1*IN- 0.3d2)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC273(IPOS,1 )
      G02 = BNC273(IPOS,2 )
      G03 = BNC273(IPOS,3 )
      G04 = BNC273(IPOS,4 )
      G05 = BNC273(IPOS,5 ) 
      G06 = BNC273(IPOS,6 )
      G07 = BNC273(IPOS,7 )
      G08 = BNC273(IPOS,8 )
      G09 = BNC273(IPOS,9 )
      G10 = BNC273(IPOS,10)
      G11 = BNC273(IPOS,11)
      G12 = BNC273(IPOS,12)
      
      ! Return to calling program
      END SUBROUTINE KM273

!------------------------------------------------------------------------------

      SUBROUTINE KM298( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!******************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM298
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!
!      TEMPERATURE IS 298K
!
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: IN
      REAL*8, INTENT(OUT) :: G01, G02, G03, G04, G05, G06 
      REAL*8, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IPOS

      !=================================================================
      ! KMC298 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.3d02) THEN
         IPOS = NINT( 0.2d2*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.2d1*IN- 0.3d2)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC298(IPOS,1 )
      G02 = BNC298(IPOS,2 )
      G03 = BNC298(IPOS,3 )
      G04 = BNC298(IPOS,4 )
      G05 = BNC298(IPOS,5 )
      G06 = BNC298(IPOS,6 )
      G07 = BNC298(IPOS,7 )
      G08 = BNC298(IPOS,8 )
      G09 = BNC298(IPOS,9 )
      G10 = BNC298(IPOS,10)
      G11 = BNC298(IPOS,11)
      G12 = BNC298(IPOS,12)

      ! Return to calling program
      END SUBROUTINE KM298
 
!------------------------------------------------------------------------------

      SUBROUTINE KM323( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!******************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM323
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!  
!      TEMPERATURE IS 323K
!  
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!  
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!******************************************************************************
!
      REAL*8, INTENT(IN)  :: IN
      REAL*8, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*8, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IPOS

      !=================================================================
      ! KMC323 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.3d2) THEN
         ipos = NINT( 0.2d2*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.2d1*IN- 0.3d2)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC323(IPOS,1 )
      G02 = BNC323(IPOS,2 )
      G03 = BNC323(IPOS,3 )
      G04 = BNC323(IPOS,4 )
      G05 = BNC323(IPOS,5 )
      G06 = BNC323(IPOS,6 )
      G07 = BNC323(IPOS,7 )
      G08 = BNC323(IPOS,8 ) 
      G09 = BNC323(IPOS,9 )
      G10 = BNC323(IPOS,10)
      G11 = BNC323(IPOS,11)
      G12 = BNC323(IPOS,12)

      ! Return to calling program
      END SUBROUTINE KM323

!------------------------------------------------------------------------------

      SUBROUTINE INIT_KMC
!
!******************************************************************************
!  Subroutine CLEANUP_KMC initializes and zeroes all module arrays. 
!  (rjp, bmy, 9/23/02)
!******************************************************************************
!


      ! Return to calling program
      END SUBROUTINE INIT_KMC

!------------------------------------------------------------------------------

      SUBROUTINE CHRBLN( STR, IBLK )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE CHRBLN
!  Purpose        : Position of last non-blank character in a string
!  Author         : Athanasios Nenes
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  STR        is the CHARACTER variable containing the string examined
!  IBLK       is a INTEGER variable containing the position of last non 
!             blank character. If string is all spaces (ie '   '), then
!             the value returned is 1.
!
!  EXAMPLE:
!             STR = 'TEST1.DAT     '
!             CALL CHRBLN (STR, IBLK)
!          
!  after execution of this code segment, "IBLK" has the value "9", which
!  is the position of the last non-blank character of "STR".
!
!***********************************************************************
!
      ! Arguments
      CHARACTER(LEN=*) :: STR
      INTEGER          :: IBLK

      ! Local variables
      INTEGER          :: K, ILEN
      
      !=================================================================
      ! CHRBLN begins here!
      !=================================================================
      IBLK = 1                       ! Substring pointer (default=1)
      ILEN = LEN(STR)                ! Length of string 
      DO K = ILEN, 1, -1              
         IF (STR(K:K).NE.' ' .AND. STR(K:K).NE.CHAR(0)) THEN
            IBLK = K
            RETURN
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE CHRBLN

!------------------------------------------------------------------------------

      SUBROUTINE SHFTRGHT( CHR )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE SHFTRGHT
!  Purpose        : RIGHT-JUSTIFICATION FUNCTION ON A STRING
!  Author         : Athanasios Nenes
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  STRING     is the CHARACTER variable with the string to be justified
!
!  EXAMPLE:
!             STRING    = 'AAAA    '
!             CALL SHFTRGHT (STRING)
!          
!  after execution of this code segment, STRING contains the value
!  '    AAAA'.
!
!*************************************************************************
!
      ! Arguments
      CHARACTER(LEN=*) :: CHR

      ! Local variables
      INTEGER          :: K, I1, I2
      
      !=================================================================
      ! SHFTRGHT begins here!
      !=================================================================
      I1  = LEN(CHR)             ! Total length of string
      CALL CHRBLN(CHR,I2)        ! Position of last non-blank character
      IF (I2.EQ.I1) RETURN

      DO K = I2, 1, -1            ! Shift characters
         CHR(I1+K-I2:I1+K-I2) = CHR(K:K)
         CHR(K:K) = ' '
      ENDDO

      ! Return to calling program
      END SUBROUTINE SHFTRGHT

!------------------------------------------------------------------------------

      SUBROUTINE RPLSTR( STRING, OLD, NEW, IERR )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE RPLSTR
!  Purpose        : REPLACE CHARACTERS OCCURING IN A STRING
!  Author         : Athanasios Nenes
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  STRING     is the CHARACTER variable with the string to be edited
!  OLD        is the old character which is to be replaced
!  NEW        is the new character which OLD is to be replaced with
!  IERR       is 0 if everything went well, is 1 if 'NEW' contains 'OLD'.
!             In this case, this is invalid, and no change is done.
!
!  EXAMPLE:
!             STRING    = 'AAAA'
!             OLD       = 'A'
!             NEW       = 'B' 
!             CALL RPLSTR (STRING, OLD, NEW)
!          
!  after execution of this code segment, STRING contains the value
!  'BBBB'.
!
!*************************************************************************
!
      ! Arguments
      CHARACTER(LEN=*) :: STRING, OLD, NEW
      INTEGER          :: IERR

      ! Local variables
      INTEGER          :: ILO, IP

      !=================================================================
      ! RPLSTR begins here!
      !=================================================================
      ILO = LEN(OLD)

      ! CHECK AND SEE IF 'NEW' CONTAINS 'OLD', WHICH CANNOT
      IP = INDEX(NEW,OLD)
      IF (IP.NE.0) THEN
         IERR = 1
         RETURN
      ELSE
         IERR = 0
      ENDIF

      ! PROCEED WITH REPLACING
10    IP = INDEX(STRING,OLD)      ! SEE IF 'OLD' EXISTS IN 'STRING'
      IF (IP.EQ.0) RETURN         ! 'OLD' DOES NOT EXIST ; RETURN
      STRING(IP:IP+ILO-1) = NEW   ! REPLACE SUBSTRING 'OLD' WITH 'NEW'
      GOTO 10                     ! GO FOR NEW OCCURANCE OF 'OLD'

      ! Return to calling program
      END SUBROUTINE RPLSTR
        
!------------------------------------------------------------------------------

      SUBROUTINE INPTD( VAR, DEF, PROMPT, PRFMT, IERR )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE INPTD
!  Purpose        : Prompts user for a value (DOUBLE). A default value
!                   is provided, so if user presses <Enter>, the default
!                   is used. 
!  Author         : Athanasios Nenes
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  VAR        is the DOUBLE PRECISION variable which value is to be saved 
!  DEF        is a DOUBLE PRECISION variable, with the default value of VAR. !       
!  PROMPT     is a CHARACTER varible containing the prompt string.     
!  PRFMT      is a CHARACTER variable containing the FORMAT specifier
!             for the default value DEF.
!  IERR       is an INTEGER error flag, and has the values:
!             0 - No error detected.
!             1 - Invalid FORMAT and/or Invalid default value.
!             2 - Bad value specified by user
!
!  EXAMPLE:
!             CALL INPTD (VAR, 1.0D0, 'Give value for A ', '*', Ierr)
!          
!  after execution of this code segment, the user is prompted for the
!  value of variable VAR. If <Enter> is pressed (ie no value is specified)
!  then 1.0 is assigned to VAR. The default value is displayed in free-
!  format. The error status is specified by variable Ierr
!
!***********************************************************************
!
      ! Arguments
      CHARACTER(LEN=*)   :: PROMPT, PRFMT
      CHARACTER(LEN=128) :: BUFFER

      ! Local variables
      REAL*8             :: DEF,  VAR
      INTEGER            :: IERR, IEND

      !=================================================================
      ! INPTD begins here!
      !=================================================================
      IERR = 0

      ! Write default value to work buffer
      WRITE (BUFFER, FMT=PRFMT, ERR=10) DEF
      CALL CHRBLN (BUFFER, IEND)

      ! Prompt user for input and read it
      WRITE (*,*) PROMPT,' [',BUFFER(1:IEND),']: '
      READ  (*, '(A)', ERR=20, END=20) BUFFER
      CALL CHRBLN (BUFFER,IEND)

      ! Read data or set default ?
      IF ( IEND .EQ. 1 .AND. BUFFER(1:1) .EQ. ' ' ) THEN
         VAR = DEF
      ELSE
         READ( BUFFER, *, ERR=20, END=20 ) VAR
      ENDIF

      ! Return
30    RETURN

      !=================================================================`
      ! ERROR HANDLER
      !=================================================================
10    IERR = 1       ! Bad FORMAT and/or bad default value
      GOTO 30

20    IERR = 2       ! Bad number given by user
      GOTO 30

      END SUBROUTINE INPTD

!------------------------------------------------------------------------------

      SUBROUTINE PUSHEND( IUNIT )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE Pushend 
!  Purpose        : Positions the pointer of a sequential file at its end
!                   Simulates the ACCESS='APPEND' clause of a F77L OPEN
!                   statement with Standard Fortran commands.
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  Iunit      is a INTEGER variable, the file unit which the file is 
!             connected to.
!
!  EXAMPLE:
!             CALL PUSHEND (10)
!          
!  after execution of this code segment, the pointer of unit 10 is 
!  pushed to its end.
!
!***********************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: IUNIT

      ! Local variablese
      LOGICAL             :: OPNED

      !=================================================================
      ! PUSHEND begins here!
      !=================================================================

      ! INQUIRE IF Iunit CONNECTED TO FIL
      INQUIRE (UNIT=Iunit, OPENED=OPNED)
      IF (.NOT.OPNED) GOTO 25

      ! Iunit CONNECTED, PUSH POINTER TO END
10    READ (Iunit,'()', ERR=20, END=20)
      GOTO 10

      ! Return point
20    BACKSPACE( IUNIT )
25    RETURN

      ! Return to calling program
      END SUBROUTINE PUSHEND

!------------------------------------------------------------------------------

      SUBROUTINE APPENDEXT( FILENAME, DEFEXT, OVERWRITE )
!
!******************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE APPENDEXT
!  Purpose        : Fix extension in file name string
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  Filename   is the CHARACTER variable with the file name
!  Defext     is the CHARACTER variable with extension (including '.',
!             ex. '.DAT')
!  Overwrite  is a LOGICAL value, .TRUE. overwrites any existing extension
!             in "Filename" with "Defext", .FALSE. puts "Defext" only if 
!             there is no extension in "Filename".
!
!  EXAMPLE:
!             FILENAME1 = 'TEST.DAT'
!             FILENAME2 = 'TEST.DAT'
!             CALL APPENDEXT (FILENAME1, '.TXT', .FALSE.)
!             CALL APPENDEXT (FILENAME2, '.TXT', .TRUE. )
!          
!  after execution of this code segment, "FILENAME1" has the value 
!  'TEST.DAT', while "FILENAME2" has the value 'TEST.TXT'
!
!******************************************************************************
!
      ! Arguments 
      CHARACTER(LEN=*) :: Filename, Defext
      LOGICAL          :: Overwrite

      ! Local variables
      INTEGER          :: IDOT, IEND

      !=================================================================
      ! APPENDEXT begins here!
      !=================================================================
      CALL CHRBLN (Filename, Iend)

      ! Filename empty
      IF (Filename(1:1).EQ.' ' .AND. Iend.EQ.1) RETURN  

      ! Append extension ?
      Idot = INDEX (Filename, '.')                      
      IF (Idot.EQ.0) Filename = Filename(1:Iend)//Defext
      IF (Overwrite .AND. Idot.NE.0)
     &              Filename = Filename(:Idot-1)//Defext

      ! Return to calling program
      END SUBROUTINE APPENDEXT

!------------------------------------------------------------------------------

      SUBROUTINE POLY3( A1, A2, A3, ROOT, ISLV )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE POLY3
! *** FINDS THE REAL ROOTS OF THE THIRD ORDER ALGEBRAIC EQUATION:
!     X**3 + A1*X**2 + A2*X + A3 = 0.0
!     THE EQUATION IS SOLVED ANALYTICALLY.
!
!     PARAMETERS A1, A2, A3 ARE SPECIFIED BY THE USER. THE MINIMUM
!     NONEGATIVE ROOT IS RETURNED IN VARIABLE 'ROOT'. IF NO ROOT IS 
!     FOUND (WHICH IS GREATER THAN ZERO), ROOT HAS THE VALUE 1D30.
!     AND THE FLAG ISLV HAS A VALUE GREATER THAN ZERO.
!
!     SOLUTION FORMULA IS FOUND IN PAGE 32 OF:
!     MATHEMATICAL HANDBOOK OF FORMULAS AND TABLES
!     SCHAUM'S OUTLINE SERIES
!     MURRAY SPIEGER, McGRAW-HILL, NEW YORK, 1968
!     (GREEK TRANSLATION: BY SOTIRIOS PERSIDES, ESPI, ATHENS, 1976)
!
!     A SPECIAL CASE IS CONSIDERED SEPERATELY ; WHEN A3 = 0, THEN
!     ONE ROOT IS X=0.0, AND THE OTHER TWO FROM THE SOLUTION OF THE
!     QUADRATIC EQUATION X**2 + A1*X + A2 = 0.0
!     THIS SPECIAL CASE IS CONSIDERED BECAUSE THE ANALYTICAL FORMULA 
!     DOES NOT YIELD ACCURATE RESULTS (DUE TO NUMERICAL ROUNDOFF ERRORS)
!
! *** COPYRIGHT 1996-98, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
!#     include "isoropia.h" 

      ! Arguments
      REAL*8            :: A1, A2, A3, ROOT
      INTEGER           :: ISLV

      ! Local variables
      REAL*8, PARAMETER :: EXPON = 1.D0 / 3.D0
      REAL*8, PARAMETER :: ZERO  = 0.D0
      REAL*8, PARAMETER :: THET1 = 120.D0 / 180.D0 
      REAL*8, PARAMETER :: THET2 = 240.D0 / 180.D0
      REAL*8, PARAMETER :: PI    = 3.14159265358932d0
      REAL*8, PARAMETER :: EPS   = 1D-50
      REAL*8            :: R
      REAL*8            :: X(3), D, Q, THET, COEF, SSIG, S, TSIG, T, SQD
      INTEGER           :: I, IX
      

      !=================================================================
      ! POLY3 begins here!
      ! 
      ! SPECIAL CASE : QUADRATIC*X EQUATION 
      !=================================================================
      IF (ABS(A3).LE.EPS) THEN 
         ISLV = 1
         IX   = 1
         X(1) = ZERO
         D    = A1*A1-4.D0*A2
         IF (D.GE.ZERO) THEN
            IX   = 3
            SQD  = SQRT(D)
            X(2) = 0.5d0*(-A1+SQD)
            X(3) = 0.5d0*(-A1-SQD)
         ENDIF
      ELSE

      !=================================================================
      ! NORMAL CASE : CUBIC EQUATION
      !=================================================================

         ! DEFINE PARAMETERS Q, R, S, T, D 
         ISLV= 1
         Q   = (3.D0*A2 - A1*A1)/9.D0
         R   = (9.D0*A1*A2 - 27.D0*A3 - 2.D0*A1*A1*A1)/54.D0
         D   = Q*Q*Q + R*R
         
         !==============================================================
         ! CALCULATE ROOTS
         !==============================================================

         ! D < 0, THREE REAL ROOTS
         IF (D.LT.-EPS) THEN        ! D < -EPS  : D < ZERO
            IX   = 3
            THET = EXPON*ACOS(R/SQRT(-Q*Q*Q))
            COEF = 2.D0*SQRT(-Q)
            X(1) = COEF*COS(THET)            - EXPON*A1
            X(2) = COEF*COS(THET + THET1*PI) - EXPON*A1
            X(3) = COEF*COS(THET + THET2*PI) - EXPON*A1

         ! D = 0, THREE REAL (ONE DOUBLE) ROOTS
         ELSE IF (D.LE.EPS) THEN    ! -EPS <= D <= EPS  : D = ZERO
            IX   = 2
            SSIG = SIGN (1.D0, R)
            S    = SSIG*(ABS(R))**EXPON
            X(1) = 2.D0*S  - EXPON*A1
            X(2) =     -S  - EXPON*A1

         ! D > 0, ONE REAL ROOT
         ELSE                       ! D > EPS  : D > ZERO
            IX   = 1
            SQD  = SQRT(D)
            SSIG = SIGN (1.D0, R+SQD)       ! TRANSFER SIGN TO SSIG
            TSIG = SIGN (1.D0, R-SQD)
            S    = SSIG*(ABS(R+SQD))**EXPON ! EXPONENTIATE ABS() 
            T    = TSIG*(ABS(R-SQD))**EXPON
            X(1) = S + T - EXPON*A1
         ENDIF
      ENDIF

      !=================================================================
      ! SELECT APPROPRIATE ROOT
      !=================================================================
      ROOT = 1.D30
      DO I = 1, IX
         IF (X(I).GT.ZERO) THEN
            ROOT = MIN (ROOT, X(I))
            ISLV = 0
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE POLY3

!------------------------------------------------------------------------------

      SUBROUTINE POLY3B( A1, A2, A3, RTLW, RTHI, ROOT, ISLV )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE POLY3B
! *** FINDS A REAL ROOT OF THE THIRD ORDER ALGEBRAIC EQUATION:
!     X**3 + A1*X**2 + A2*X + A3 = 0.0
!     THE EQUATION IS SOLVED NUMERICALLY (BISECTION).
!
!     PARAMETERS A1, A2, A3 ARE SPECIFIED BY THE USER. THE MINIMUM
!     NONEGATIVE ROOT IS RETURNED IN VARIABLE 'ROOT'. IF NO ROOT IS 
!     FOUND (WHICH IS GREATER THAN ZERO), ROOT HAS THE VALUE 1D30.
!     AND THE FLAG ISLV HAS A VALUE GREATER THAN ZERO.
!
!     RTLW, RTHI DEFINE THE INTERVAL WHICH THE ROOT IS LOOKED FOR.
!
!  NOTES:
!  (1 ) Split off internal function FUNC and made it a module routine
!        since it is not possible to call internal routines from w/in
!        a parallelized region.  Now we must call FUNC(A1,A2,A3,X) 
!        instead of FUNC(X). (bmy, 3/10/05)
!******************************************************************************
!
      ! Arguments
      INTEGER            :: ISLV
      REAL*8             :: A1, A2, A3, RTLW, RTHI, ROOT
 
      ! Local variables
      INTEGER, PARAMETER :: MAXIT=100
      INTEGER, PARAMETER :: NDIV=5
      INTEGER            :: I, K1, K2
      REAL*8,  PARAMETER :: ZERO=0.0D0
      REAL*8,  PARAMETER :: EPS=1.0D-15
      REAL*8             :: X1, Y1, DX, X2, Y2, Y3, X3

      !=================================================================
      ! POLY3B begins here!
      !=================================================================
       
      ! Initial values for bisection
      X1   = RTLW
      Y1   = FUNC( A1, A2, A3, X1 )

      ! Is low a root?
      IF ( ABS(Y1) .LE. EPS ) THEN     
         ROOT = RTLW
         GOTO 50
      ENDIF

      !=================================================================
      ! Root tracking ; for the range of HI and LO
      !=================================================================
      DX = (RTHI-RTLW)/FLOAT(NDIV)

      DO I = 1, NDIV
         X2 = X1+DX
         Y2 = FUNC( A1, A2, A3, X2 )

         ! (Y1*Y2.LT.ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y2) .LT. ZERO ) GOTO 20 
         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================
      ! No subdivision with solution found
      !=================================================================
      IF ( ABS( Y2) .LT. EPS ) THEN   ! X2 is a root
         ROOT = X2
      ELSE
         ROOT = 1.d30
         ISLV = 1
      ENDIF
      GOTO 50

      !=================================================================
      ! BISECTION
      !=================================================================
 20   CONTINUE

      DO I = 1, MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNC( A1, A2, A3, X3 )

         ! (Y1*Y3 .LE. ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y3) .LE. ZERO ) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS( X2 - X1 ) .LE. EPS * X1 ) GOTO 40
      ENDDO

      !=================================================================
      ! Converged ; RETURN
      !=================================================================
 40   CONTINUE
      X3   = 0.5d0*(X1+X2)
      Y3   = FUNC( A1, A2, A3, X3 )
      ROOT = X3
      ISLV = 0

      ! Return to calling program
 50   CONTINUE

      ! Return to calling program
      END SUBROUTINE POLY3B
      
!------------------------------------------------------------------------------

      FUNCTION FUNC( A1, A2, A3, X ) RESULT ( VALUE )
!
!******************************************************************************
!  This was split off from POLY3B (bmy, 3/10/05)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: A1, A2, A3, X

      ! Function value
      REAL*8             :: VALUE

      ! This replaces the old statement function
      VALUE = X**3.d0 + A1*X**2.d0 + A2*X + A3

      ! Eventually replace the polynomial w/ this function below
      !VALUE = A3 + ( X * ( A2 + X * ( A1 + X )))

      ! Return to POLY3B
      END FUNCTION FUNC

!------------------------------------------------------------------------------

      FUNCTION EX10( X, K ) RESULT( VALUE )
! 
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION EX10
! *** 10^X FUNCTION ; ALTERNATE OF LIBRARY ROUTINE ; USED BECAUSE IT IS
!     MUCH FASTER BUT WITHOUT GREAT LOSS IN ACCURACY. , 
!     MAXIMUM ERROR IS 2%, EXECUTION TIME IS 42% OF THE LIBRARY ROUTINE 
!     (ON A 80286/80287 MACHINE, using Lahey FORTRAN 77 v.3.0).
!
!     EXPONENT RANGE IS BETWEEN -K AND K (K IS THE REAL ARGUMENT 'K')
!     MAX VALUE FOR K: 9.999
!     IF X < -K, X IS SET TO -K, IF X > K, X IS SET TO K
!
!     THE EXPONENT IS CALCULATED BY THE PRODUCT ADEC*AINT, WHERE ADEC
!     IS THE MANTISSA AND AINT IS THE MAGNITUDE (EXPONENT). BOTH 
!     MANTISSA AND MAGNITUDE ARE PRE-CALCULATED AND STORED IN LOOKUP
!     TABLES ; THIS LEADS TO THE INCREASED SPEED.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: X, K

      ! Return value
      REAL*8             :: VALUE

      ! Local variables
!      REAL*4             :: EX10, Y
      REAL*8             :: Y 
      INTEGER            :: K1, K2

      ! For integer part
      REAL*4, SAVE       :: AINT10(20) = (/
     & 0.1000d-08, 0.1000d-07, 0.1000d-06, 0.1000d-05, 0.1000d-04,
     & 0.1000d-03, 0.1000d-02, 0.1000d-01, 0.1000d+00, 0.1000d+01,
     & 0.1000d+02, 0.1000d+03, 0.1000d+04, 0.1000d+05, 0.1000d+06,
     & 0.1000d+07, 0.1000d+08, 0.1000d+09, 0.1000d+10, 0.1000d+11/)

      ! For decimal part
      REAL*4, SAVE       :: ADEC10(200) = (/
     & 0.1023d+00, 0.1047d+00, 0.1072d+00, 0.1096d+00, 0.1122d+00,
     & 0.1148d+00, 0.1175d+00, 0.1202d+00, 0.1230d+00, 0.1259d+00,
     & 0.1288d+00, 0.1318d+00, 0.1349d+00, 0.1380d+00, 0.1413d+00,
     & 0.1445d+00, 0.1479d+00, 0.1514d+00, 0.1549d+00, 0.1585d+00,
     & 0.1622d+00, 0.1660d+00, 0.1698d+00, 0.1738d+00, 0.1778d+00,
     & 0.1820d+00, 0.1862d+00, 0.1905d+00, 0.1950d+00, 0.1995d+00,
     & 0.2042d+00, 0.2089d+00, 0.2138d+00, 0.2188d+00, 0.2239d+00,
     & 0.2291d+00, 0.2344d+00, 0.2399d+00, 0.2455d+00, 0.2512d+00,
     & 0.2570d+00, 0.2630d+00, 0.2692d+00, 0.2754d+00, 0.2818d+00,
     & 0.2884d+00, 0.2951d+00, 0.3020d+00, 0.3090d+00, 0.3162d+00,
     & 0.3236d+00, 0.3311d+00, 0.3388d+00, 0.3467d+00, 0.3548d+00,
     & 0.3631d+00, 0.3715d+00, 0.3802d+00, 0.3890d+00, 0.3981d+00,
     & 0.4074d+00, 0.4169d+00, 0.4266d+00, 0.4365d+00, 0.4467d+00,
     & 0.4571d+00, 0.4677d+00, 0.4786d+00, 0.4898d+00, 0.5012d+00,
     & 0.5129d+00, 0.5248d+00, 0.5370d+00, 0.5495d+00, 0.5623d+00,
     & 0.5754d+00, 0.5888d+00, 0.6026d+00, 0.6166d+00, 0.6310d+00,
     & 0.6457d+00, 0.6607d+00, 0.6761d+00, 0.6918d+00, 0.7079d+00,
     & 0.7244d+00, 0.7413d+00, 0.7586d+00, 0.7762d+00, 0.7943d+00,
     & 0.8128d+00, 0.8318d+00, 0.8511d+00, 0.8710d+00, 0.8913d+00,
     & 0.9120d+00, 0.9333d+00, 0.9550d+00, 0.9772d+00, 0.1000d+01,
     & 0.1023d+01, 0.1047d+01, 0.1072d+01, 0.1096d+01, 0.1122d+01,
     & 0.1148d+01, 0.1175d+01, 0.1202d+01, 0.1230d+01, 0.1259d+01,
     & 0.1288d+01, 0.1318d+01, 0.1349d+01, 0.1380d+01, 0.1413d+01,
     & 0.1445d+01, 0.1479d+01, 0.1514d+01, 0.1549d+01, 0.1585d+01,
     & 0.1622d+01, 0.1660d+01, 0.1698d+01, 0.1738d+01, 0.1778d+01,
     & 0.1820d+01, 0.1862d+01, 0.1905d+01, 0.1950d+01, 0.1995d+01,
     & 0.2042d+01, 0.2089d+01, 0.2138d+01, 0.2188d+01, 0.2239d+01,
     & 0.2291d+01, 0.2344d+01, 0.2399d+01, 0.2455d+01, 0.2512d+01,
     & 0.2570d+01, 0.2630d+01, 0.2692d+01, 0.2754d+01, 0.2818d+01,
     & 0.2884d+01, 0.2951d+01, 0.3020d+01, 0.3090d+01, 0.3162d+01,
     & 0.3236d+01, 0.3311d+01, 0.3388d+01, 0.3467d+01, 0.3548d+01,
     & 0.3631d+01, 0.3715d+01, 0.3802d+01, 0.3890d+01, 0.3981d+01,
     & 0.4074d+01, 0.4169d+01, 0.4266d+01, 0.4365d+01, 0.4467d+01,
     & 0.4571d+01, 0.4677d+01, 0.4786d+01, 0.4898d+01, 0.5012d+01,
     & 0.5129d+01, 0.5248d+01, 0.5370d+01, 0.5495d+01, 0.5623d+01,
     & 0.5754d+01, 0.5888d+01, 0.6026d+01, 0.6166d+01, 0.6310d+01,
     & 0.6457d+01, 0.6607d+01, 0.6761d+01, 0.6918d+01, 0.7079d+01,
     & 0.7244d+01, 0.7413d+01, 0.7586d+01, 0.7762d+01, 0.7943d+01,
     & 0.8128d+01, 0.8318d+01, 0.8511d+01, 0.8710d+01, 0.8913d+01,
     & 0.9120d+01, 0.9333d+01, 0.9550d+01, 0.9772d+01, 0.1000d+02
     & /)

      !=================================================================
      ! EX10 begins here!
      !=================================================================
      
      ! Limit X to [-K, K] range:  MIN: -9.999, MAX: 9.999
      Y    = MAX(-K, MIN(X,K))   

      ! Get integer and decimal part
      K1   = INT(Y)
      K2   = INT(100*(Y-K1))

      ! Calculate EXP function
      VALUE = AINT10(K1+10) * ADEC10(K2+100)

      ! Return to calling program
      END FUNCTION EX10

!------------------------------------------------------------------------------

      SUBROUTINE PUSHERR( IERR, ERRINF )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE PUSHERR
! *** THIS SUBROUTINE SAVES AN ERROR MESSAGE IN THE ERROR STACK
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      INTEGER,          INTENT(IN) :: IERR
      CHARACTER(LEN=*), INTENT(IN) :: ERRINF

      ! Local variables
      

      !=================================================================
      ! PUSHERR begins here!
      !=================================================================

      ! SAVE ERROR CODE IF THERE IS ANY SPACE
      IF ( NOFER .LT. NERRMX ) THEN   
         NOFER         = NOFER + 1 
         ERRSTK(NOFER) = IERR
         ERRMSG(NOFER) = ERRINF   
         STKOFL        =.FALSE.

      ! STACK OVERFLOW
      ELSE
         STKOFL        =.TRUE.      
      ENDIF

      ! Return to calling program
      END SUBROUTINE PUSHERR

!------------------------------------------------------------------------------

      SUBROUTINE ISERRINF( ERRSTKI, ERRMSGI, NOFERI, STKOFLI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISERRINF
! *** THIS SUBROUTINE OBTAINS A COPY OF THE ERROR STACK (& MESSAGES) 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      CHARACTER(LEN=40), INTENT(OUT) :: ERRMSGI(NERRMX)
      INTEGER,           INTENT(OUT) :: ERRSTKI(NERRMX)
      INTEGER,           INTENT(OUT) :: NOFERI
      LOGICAL,           INTENT(OUT) :: STKOFLI

      ! Local variables
      INTEGER           :: I

      !=================================================================
      ! ISERRINF begins here!
      ! 
      ! OBTAIN WHOLE ERROR STACK
      !=================================================================
      DO I = 1, NOFER              ! Error messages & codes
         ERRSTKI(I) = ERRSTK(I)
         ERRMSGI(I) = ERRMSG(I)
      ENDDO

      STKOFLI = STKOFL
      NOFERI  = NOFER

      ! Return to calling program
      END SUBROUTINE ISERRINF
      
!------------------------------------------------------------------------------

      SUBROUTINE ERRSTAT( IO, IERR, ERRINF )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE ERRSTAT
! *** THIS SUBROUTINE REPORTS ERROR MESSAGES TO UNIT 'IO'
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      INTEGER           :: IO, IERR
      CHARACTER(LEN=*)  :: ERRINF

      ! Local variables
      INTEGER           :: IOK, IEND
      CHARACTER(LEN=4)  :: CER 
      CHARACTER(LEN=29) :: NCIS = 'NO CONVERGENCE IN SUBROUTINE '
      CHARACTER(LEN=27) :: NCIF = 'NO CONVERGENCE IN FUNCTION ' 
      CHARACTER(LEN=26) :: NSIS = 'NO SOLUTION IN SUBROUTINE '  
      CHARACTER(LEN=24) :: NSIF = 'NO SOLUTION IN FUNCTION '  

      !=================================================================
      ! ERRSTAT begins here!
      !
      ! WRITE ERROR IN CHARACTER
      !=================================================================
      WRITE (CER,'(I4)') IERR
      CALL RPLSTR (CER, ' ', '0',IOK)   ! REPLACE BLANKS WITH ZEROS
      CALL CHRBLN (ERRINF, IEND)        ! LAST POSITION OF ERRINF CHAR

      !=================================================================
      ! WRITE ERROR TYPE (FATAL, WARNING) 
      !=================================================================    
      IF ( IERR .EQ. 0 ) THEN
         WRITE(IO,1000) 'NO ERRORS DETECTED '
         GOTO 10

      ELSE IF ( IERR .LT. 0 ) THEN
         WRITE(IO,1000) 'ERROR STACK EXHAUSTED '
         GOTO 10

      ELSE IF ( IERR .GT. 1000 ) THEN
         WRITE(IO,1100) 'FATAL',CER

      ELSE
         WRITE(IO,1100) 'WARNING',CER

      ENDIF

      !=================================================================
      ! WRITE ERROR MESSAGE
      !================================================================= 

      ! FATAL MESSAGES
      IF ( IERR .EQ. 1001 ) THEN 
         CALL CHRBLN (SCASE, IEND)
         WRITE(IO,1000) 'CASE NOT SUPPORTED IN CALCMR ['//SCASE(1:IEND)
     &                  //']'

      ELSE IF ( IERR .EQ. 1002 ) THEN 
         CALL CHRBLN (SCASE, IEND)
         WRITE(IO,1000) 'CASE NOT SUPPORTED ['//SCASE(1:IEND)//']'

      ! WARNING MESSAGES
      ELSE IF ( IERR .EQ. 0001 ) THEN 
         WRITE(IO,1000) NSIS,ERRINF

      ELSE IF ( IERR .EQ. 0002 ) THEN 
         WRITE(IO,1000) NCIS,ERRINF

      ELSEIF ( IERR .EQ. 0003 ) THEN 
         WRITE(IO,1000) NSIF,ERRINF

      ELSE IF ( IERR .EQ. 0004 ) THEN 
         WRITE(IO,1000) NCIF,ERRINF

      ELSE IF ( IERR .EQ. 0019 ) THEN
         WRITE(IO,1000) 'HNO3(aq) AFFECTS H+, WHICH '//
     &                  'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE(IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSE IF ( IERR .EQ. 0020 ) THEN
         IF (W(4).GT.TINY .AND. W(5).GT.TINY) THEN
            WRITE(IO,1000) 'HSO4-SO4 EQUILIBRIUM MIGHT AFFECT HNO3,'
     &                    //'HCL DISSOLUTION'
         ELSE
            WRITE(IO,1000) 'HSO4-SO4 EQUILIBRIUM MIGHT AFFECT NH3 '
     &                    //'DISSOLUTION'
         ENDIF
         WRITE(IO,1000) 'DIRECT DECREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSE IF ( IERR .EQ. 0021 ) THEN
         WRITE(IO,1000) 'HNO3(aq),HCL(aq) AFFECT H+, WHICH '//
     &                  'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE(IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSE IF ( IERR .EQ. 0022 ) THEN
         WRITE(IO,1000) 'HCL(g) EQUILIBRIUM YIELDS NONPHYSICAL '//
     &                  'DISSOLUTION'
         WRITE(IO,1000) 'A TINY AMOUNT [',ERRINF(1:IEND),'] IS '//
     &                  'ASSUMED TO BE DISSOLVED'

      ELSE IF ( IERR .EQ. 0033 ) THEN
         WRITE(IO,1000) 'HCL(aq) AFFECTS H+, WHICH '//
     &                  'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE(IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSEIF ( IERR .EQ. 0050 ) THEN
         WRITE(IO,1000) 'TOO MUCH SODIUM GIVEN AS INPUT.'
         WRITE(IO,1000) 'REDUCED TO COMPLETELY NEUTRALIZE SO4,Cl,NO3.'
         WRITE(IO,1000) 'EXCESS SODIUM IS IGNORED.'

      ELSE
         WRITE(IO,1000) 'NO DIAGNOSTIC MESSAGE AVAILABLE'
      ENDIF

      ! Return
10    RETURN

      ! FORMAT statements
1000  FORMAT (1X,A:A:A:A:A)
1100  FORMAT (1X,A,' ERROR [',A4,']:')

      ! Return to calling program
      END SUBROUTINE ERRSTAT

!------------------------------------------------------------------------------

      SUBROUTINE ISRP2F( WI, RHI, TEMPI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISRP2F
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF 
!     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM. 
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!     THE AMBIENT RELATIVE HUMIDITY.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "CMN_SIZE" ! IIPAR, JPAR, LLTROP
#     include "isoropia.h" 

      REAL*8,  INTENT(IN) :: WI(NCOMP), RHI, TEMPI

      !=================================================================
      ! ISRP2F begins here!
      !=================================================================
      
      ! Allocate all module arrays and some common blocks
      CALL INIT2( WI, RHI, TEMPI )

      ! CALCULATE SULFATE RATIO
      SULRAT = W(3)/W(2)

      !=================================================================
      ! FIND CALCULATION REGIME FROM (SULRAT,RH)
      !
      ! SULFATE POOR 
      !=================================================================
      IF ( 2.0 .LE. SULRAT ) THEN                

         ! Only liquid (metastable)
         IF ( METSTBL .EQ. 1 ) THEN
            SCASE = 'D3'
            CALL CALCD3     
         ELSE

            ! NH42SO4,NH4NO3       ; case D1
            IF ( RHB .LT. DRNH4NO3 ) THEN    
               SCASE = 'D1'
               CALL CALCD1      
     
            ! NH42S4               ; case D2
            ELSE IF ( DRNH4NO3 .LE. RHB .AND. RHB .LT. DRNH42S4 ) THEN         
               SCASE = 'D2'
               CALL CALCD2      
     
            ! Only liquid          ; case D3
            ELSE IF ( DRNH42S4 .LE. RHB ) THEN
               SCASE = 'D3'
               CALL CALCD3      
            ENDIF
         ENDIF

      !=================================================================
      ! SULFATE RICH (NO ACID)
      !
      ! FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
      ! THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
      ! SUBROUTINES CALCB? ARE CALLED, AND THEN THE NITRIC ACID IS 
      ! DISSOLVED FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
      !=================================================================
      ELSE IF ( 1.0 .LE. SULRAT .AND. SULRAT .LT. 2.0 ) THEN 

         IF ( METSTBL .EQ. 1 ) THEN
            SCASE = 'B4'
            CALL CALCB4         ! Only liquid (metastable)
            SCASE = 'E4'
         ELSE
            
            ! NH4HSO4,LC,NH42SO4   ; case E1
            IF ( RHB .LT. DRNH4HS4 ) THEN         
               SCASE = 'B1'
               CALL CALCB1   
               SCASE = 'E1'

            ! LC,NH42S4            ; case E2
            ELSE IF ( DRNH4HS4 .LE. RHB .AND. RHB .LT. DRLC ) THEN         
               SCASE = 'B2'
               CALL CALCB2      
               SCASE = 'E2'
               
            ! NH42S4               ; case E3
            ELSE IF ( DRLC .LE. RHB .AND. RHB .LT. DRNH42S4 ) THEN         
               SCASE = 'B3'
               CALL CALCB3      
               SCASE = 'E3'
               
            ! Only liquid          ; case E4
            ELSE IF ( DRNH42S4 .LE. RHB ) THEN         
               SCASE = 'B4'
               CALL CALCB4      
               SCASE = 'E4'
            ENDIF
         ENDIF
         
         ! HNO3(g) DISSOLUTION
         CALL CALCNA            

      !=================================================================
      ! SULFATE RICH (FREE ACID)
      !
      ! FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
      ! THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM
      ! SUBROUTINE CALCC? IS CALLED, AND THEN THE NITRIC ACID IS 
      ! DISSOLVED FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
      !=================================================================
      ELSE IF ( SULRAT .LT. 1.0 ) THEN             

         ! Only liquid (metastable)
         IF ( METSTBL .EQ. 1 ) THEN
            SCASE = 'C2'
            CALL CALCC2            
            SCASE = 'F2'
         ELSE

            ! NH4HSO4              ; case F1
            IF ( RHB .LT. DRNH4HS4 ) THEN         
               SCASE = 'C1'
               CALL CALCC1         
               SCASE = 'F1'

            ! Only liquid          ; case F2
            ELSE IF ( DRNH4HS4 .LE. RHB ) THEN         
               SCASE = 'C2'
               CALL CALCC2              
               SCASE = 'F2'
            ENDIF
         ENDIF

         ! HNO3(g) DISSOLUTION
         CALL CALCNA 
       
      ENDIF

      ! Return to calling program 
      END SUBROUTINE ISRP2F

!-----------------------------------------------------------------------------

      SUBROUTINE ISRP3F( WI, RHI, TEMPI )

!*****************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISRP3F
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
!     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM. 
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM 
!     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
!
! *** COPYRIGHT 1996-2000 UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!*****************************************************************************

#     include "CMN_SIZE" ! IIPAR, JPAR, LLTROP
#     include "isoropia.h" 

      REAL*8 , INTENT(IN) :: WI(NCOMP), RHI, TEMPI
      REAL*8              :: REST, WIN(NCOMP)
      INTEGER             :: K

      !=================================================================
      ! ISRP3F begins here!
      !=================================================================

      DO K = 1, NCOMP
 	WIN(K) = WI(K)
      ENDDO

      ! Adjust for too little ammonium and chloride 
      WIN(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
      WIN(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3

      ! Adjust for too little sodium, sulfate, and nitrate combined
      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
         WIN(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
         WIN(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
      ENDIF

      ! Allocate all module arrays and some common blocks
      CALL INIT3 (WIN, RHI, TEMPI)

      ! Check if too much sodium, adjust and issue an error message
      REST = 2.D0*W(2) + W(4) + W(5) 
      IF (W(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
         W(1) = (ONE-1D-6)*REST         ! Adjust Na amount
         CALL PUSHERR (0050, 'ISRP3F')  ! Warning error: Na adjusted
      ENDIF

      ! Calculate sulfate and sodium ratios
      SULRAT = (W(1)+W(3))/W(2)
      SODRAT = W(1)/W(2)

      ! Find calculation regime from SULRAT, RH
      ! Sulfate poor; sodium poor
      IF ( 2.0 .LE. SULRAT .AND. SODRAT .LT. 2.0 ) THEN                
      IF( METSTBL .EQ. 1 ) THEN
         SCASE = 'G5'
         CALL CALCG5                 ! Only liquid (metastable)
      ELSE

         IF ( RHB .LT. DRNH4NO3 ) THEN    
            SCASE = 'G1'
            CALL CALCG1              ! NH42SO4,NH4NO3,NH4CL,NA2SO4

         ELSEIF ( DRNH4NO3 .LE. RHB .AND. RHB .LT. DRNH4CL ) THEN         
            SCASE = 'G2'
            CALL CALCG2              ! NH42SO4,NH4CL,NA2SO4

         ELSEIF ( DRNH4CL .LE. RHB  .AND. RHB .LT. DRNH42S4 ) THEN         
            SCASE = 'G3'
            CALL CALCG3              ! NH42SO4,NA2SO4
 
        ELSEIF ( DRNH42S4 .LE. RHB  .AND. RHB .LT. DRNA2SO4 ) THEN         
            SCASE = 'G4'
            CALL CALCG4              ! NA2SO4

         ELSEIF ( DRNA2SO4 .LE. RHB ) THEN         
            SCASE = 'G5'
            CALL CALCG5             ! Only liquid
         ENDIF
      ENDIF

      ! Sulfate poor, sodium rich
      ELSE IF ( SULRAT .GE. 2.0 .AND. SODRAT .GE. 2.0 ) THEN                

      IF( METSTBL .EQ. 1 ) THEN
         SCASE = 'H6'
         CALL CALCH6                 ! Only liquid (metastable)
      ELSE

         IF ( RHB .LT. DRNH4NO3 ) THEN    
            SCASE = 'H1'
            CALL CALCH1              ! NH4NO3,NH4CL,NA2SO4,NACL,NANO3

         ELSEIF ( DRNH4NO3 .LE. RHB .AND. RHB .LT. DRNANO3 ) THEN         
            SCASE = 'H2'
            CALL CALCH2             ! NH4CL,NA2SO4,NACL,NANO3

         ELSEIF ( DRNANO3 .LE. RHB  .AND. RHB .LT. DRNACL ) THEN         
            SCASE = 'H3'
            CALL CALCH3             ! NH4CL,NA2SO4,NACL

         ELSEIF ( DRNACL .LE. RHB   .AND. RHB .LT. DRNH4Cl ) THEN         
            SCASE = 'H4'
            CALL CALCH4             ! NH4CL,NA2SO4

         ELSEIF ( DRNH4Cl .LE. RHB .AND. RHB .LT. DRNA2SO4 ) THEN         
            SCASE = 'H5'
            CALL CALCH5              ! NA2SO4

         ELSEIF ( DRNA2SO4 .LE. RHB ) THEN         
            SCASE = 'H6'
            CALL CALCH6              ! NO SOLID
         ENDIF
      ENDIF

      ! Sulfate rich (no acid)
      ELSEIF ( 1.0 .LE. SULRAT .AND. SULRAT .LT. 2.0 ) THEN 

      IF( METSTBL .EQ. 1 ) THEN
         SCASE = 'I6'
         CALL CALCI6                 ! Only liquid (metastable)
      ELSE

         IF ( RHB .LT. DRNH4HS4 ) THEN         
            SCASE = 'I1'
            CALL CALCI1              ! NA2SO4,(NH4)2SO4,NAHSO4,NH4HSO4,LC

         ELSEIF ( DRNH4HS4 .LE. RHB .AND. RHB .LT. DRNAHSO4 ) THEN         
            SCASE = 'I2'
            CALL CALCI2              ! NA2SO4,(NH4)2SO4,NAHSO4,LC

         ELSEIF ( DRNAHSO4 .LE. RHB .AND. RHB .LT. DRLC ) THEN         
            SCASE = 'I3'
            CALL CALCI3              ! NA2SO4,(NH4)2SO4,LC

         ELSEIF ( DRLC .LE. RHB     .AND. RHB .LT. DRNH42S4 ) THEN         
            SCASE = 'I4'
            CALL CALCI4              ! NA2SO4,(NH4)2SO4

         ELSEIF ( DRNH42S4 .LE. RHB .AND. RHB .LT. DRNA2SO4 ) THEN         
            SCASE = 'I5'
            CALL CALCI5              ! NA2SO4

         ELSEIF ( DRNA2SO4 .LE. RHB ) THEN         
            SCASE = 'I6'
            CALL CALCI6              ! NO SOLIDS
         ENDIF
      ENDIF
                                    
      CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl       
      CALL CALCNH3                !                NH3 

      ! Sulfate rich (free acid)
      ELSEIF ( SULRAT .LT. 1.0 ) THEN             

      IF( METSTBL .EQ. 1 ) THEN
         SCASE = 'J3'
         CALL CALCJ3                 ! Only liquid (metastable)
      ELSE

         IF ( RHB .LT. DRNH4HS4 ) THEN         
            SCASE = 'J1'
            CALL CALCJ1              ! NH4HSO4,NAHSO4

         ELSEIF ( DRNH4HS4 .LE. RHB .AND. RHB .LT. DRNAHSO4 ) THEN         
            SCASE = 'J2'
            CALL CALCJ2              ! NAHSO4

         ELSEIF ( DRNAHSO4 .LE. RHB ) THEN         
            SCASE = 'J3'
            CALL CALCJ3              
         ENDIF
      ENDIF
                                    
      CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl       
      CALL CALCNH3                !                NH3 
      ENDIF

      ! Return to calling program 
      END SUBROUTINE ISRP3F

!------------------------------------------------------------------------------

      SUBROUTINE CALCB4
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB4
! *** CASE B4 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
!
!     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
!     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
!     AND THAT CALCULATED FROM ELECTRONEUTRALITY.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER :: K
      REAL*8  :: AK1, BET, GAM, BB, CC, DD

      !=================================================================
      ! CALCB4 begins here!
      !=================================================================

      ! Set flags to solve equations
      FRST       = .TRUE.
      CALAIN     = .TRUE.
      CALAOU     = .TRUE.

      !=================================================================
      ! CALCULATE WATER CONTENT
      !=================================================================

      ! Get dry salt content, and use for water.
      CALL CALCB1A       

      MOLALR(13) = CLC       
      MOLALR(9)  = CNH4HS4   
      MOLALR(4)  = CNH42S4   
      CLC        = ZERO
      CNH4HS4    = ZERO
      CNH42S4    = ZERO
      WATER      = MOLALR(13)/M0(13)+MOLALR(9)/
     &             M0(9)+MOLALR(4)/M0(4)

      ! NH4I
      MOLAL(3)   = W(3)   

      DO K = 1, NSWEEP
         AK1   = XK1*((GAMA(8)/GAMA(7))**2.d0)*(WATER/GAMA(7))
         BET   = W(2)
         GAM   = MOLAL(3)
         BB    = BET + AK1 - GAM
         CC    =-AK1*BET
         DD    = BB*BB - 4.D0*CC

         ! Speciation & water content
         MOLAL (5) = MIN(0.5d0*(-BB + SQRT(DD)), W(2)) ! SO4I
         MOLAL (6) = MAX(TINY,MIN(W(2)-MOLAL(5),
     &                     W(2)))         ! HSO4I
         MOLAL (1) = MAX(TINY,MIN(AK1*MOLAL(6)/
     &                     MOLAL(5),W(2))) ! HI
         CALL CALCMR              ! Water content

         ! Calculate activities or terminate internal loop 
         IF ( .NOT. CALAIN ) GOTO 30
         CALL CALCACT

      ENDDO

      ! Return to calling program
30    CONTINUE
      END SUBROUTINE CALCB4

!------------------------------------------------------------------------------

      SUBROUTINE CALCB3
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB3
! *** CASE B3 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!     3. SOLIDS POSSIBLE: (NH4)2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8 :: X, Y, TLC, TNH42S4, TNH4HS4

      !=================================================================
      ! CALCB3 begins here!
      !=================================================================

      ! Calculate equivalent amount of HSO4 and SO4
      X = MAX(2*W(2)-W(3), ZERO)   ! Equivalent NH4HSO4
      Y = MAX(W(3)  -W(2), ZERO)   ! Equivalent NH42SO4

      ! Calculate species according to relative abundance of HSO4
      IF ( X .LT. Y ) THEN               ! LC is the MIN (x,y)
         SCASE   = 'B3 ; SUBCASE 1'
         TLC     = X
         TNH42S4 = Y-X
         CALL CALCB3A (TLC,TNH42S4)      ! LC + (NH4)2SO4 
      ELSE
         SCASE   = 'B3 ; SUBCASE 2'
         TLC     = Y
         TNH4HS4 = X-Y
         CALL CALCB3B (TLC,TNH4HS4)      ! LC + NH4HSO4
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB3

!------------------------------------------------------------------------------

      SUBROUTINE CALCB3A( TLC, TNH42S4 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB3A
! *** CASE B3 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!     3. SOLIDS POSSIBLE: (NH4)2SO4
!
!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
!     AMOUNT OF SOLID (NH4)2SO4 DISSOLVED IN THE LIQUID PHASE.
!     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB3A CALCULATES THE
!     AMOUNT OF H+ PRODUCED (BASED ON THE SO4 RELEASED INTO THE
!     SOLUTION). THE SOLUBILITY PRODUCT OF (NH4)2SO4 IS USED AS THE 
!     OBJECTIVE FUNCTION.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 


      ! Arguments
      REAL*8  :: TLC, TNH42S4

      ! Local variables
      INTEGER :: K
!      REAL*8  :: FUNCB3A
      REAL*8  :: ZLO, ZHI, Z1, Z2, Z3, Y1, Y2, Y3, YLO, YHI, DZ, ZK

      !=================================================================
      ! CALCB3A begins here!
      !=================================================================
      CALAOU = .TRUE.         ! Outer loop activity calculation flag
      ZLO    = ZERO           ! MIN DISSOLVED (NH4)2SO4
      ZHI    = TNH42S4        ! MAX DISSOLVED (NH4)2SO4

      !=================================================================
      ! INITIAL VALUES FOR BISECTION (DISSOLVED (NH4)2SO4 
      !=================================================================
      Z1 = ZLO
      Y1 = FUNCB3A (Z1, TLC, TNH42S4)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DZ = (ZHI-ZLO)/FLOAT(NDIV)
      DO K = 1, NDIV
         Z2 = Z1+DZ
         Y2 = FUNCB3A (Z2, TLC, TNH42S4)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         Z1 = Z2
         Y1 = Y2
      ENDDO

      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION FOUND 
      !=================================================================
      YHI= Y1                    ! Save Y-value at HI position
      IF (ABS(Y2) .LT. EPS) THEN ! X2 IS A SOLUTION 
         RETURN

      ! { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
      ELSE IF ( YLO .LT. ZERO .AND. YHI. LT. ZERO ) THEN
         Z1 = ZHI
         Z2 = ZHI
         GOTO 40

      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
      ELSE IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         Z1 = ZLO
         Z2 = ZLO
         GOTO 40

      ! Warning error: no solution
      ELSE
         CALL PUSHERR( 0001, 'CALCB3A' ) 
         RETURN

      ENDIF

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
 20   CONTINUE

      DO K = 1, MAXIT
         Z3 = 0.5d0*(Z1+Z2)
         Y3 = FUNCB3A (Z3, TLC, TNH42S4)

         ! (Y1*Y3 .LE. ZERO)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            Z2    = Z3
         ELSE
            Y1    = Y3
            Z1    = Z3
         ENDIF
         IF (ABS(Z2-Z1) .LE. EPS*Z1) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR( 0002, 'CALCB3A' )

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
40    CONTINUE
      ZK = 0.5d0*(Z1+Z2)
      Y3 = FUNCB3A (ZK, TLC, TNH42S4)

      ! Return to calling program
      END SUBROUTINE CALCB3A

!------------------------------------------------------------------------------

!      FUNCTION FUNCB3A( ZK, Y, X ) 
      FUNCTION FUNCB3A( ZK, Y, X ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCB3A
! *** CASE B3 ; SUBCASE 1
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B3
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA3.
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: ZK, Y, X

      ! Return value
      REAL*8  :: VALUE

      ! Local variables
      INTEGER :: K
      REAL*8  :: KK, GRAT1, DD!, FUNCB3A
      
      !=================================================================
      ! FUNCB3A begins here!
      !=================================================================
      FRST   = .TRUE.
      CALAIN = .TRUE.

      ! Solve equations ; with iterations for activity coef.
      DO K = 1, NSWEEP
         GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
         DD    = SQRT( (ZK+GRAT1+Y)**2.d0 + 4.d0*Y*GRAT1)
         KK    = 0.5d0*(-(ZK+GRAT1+Y) + DD )

         ! Speciation & water content
         MOLAL (1) = KK                ! HI
         MOLAL (5) = KK+ZK+Y           ! SO4I
         MOLAL (6) = MAX (Y-KK, TINY)  ! HSO4I
         MOLAL (3) = 3.d0*Y+2*ZK        ! NH4I
         CNH42S4   = X-ZK              ! Solid (NH4)2SO4
         CALL CALCMR                   ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 30
         ENDIF
      ENDDO

      !=================================================================
      ! calculate objective function
      !=================================================================
30    CONTINUE
      !FUNCB3A= ( SO4I*NH4I**2.0 )/( XK7*(WATER/GAMA(4))**3.0 )      
!      FUNCB3A= MOLAL(5)*MOLAL(3)**2.0
!      FUNCB3A= FUNCB3A/(XK7*(WATER/GAMA(4))**3.0) - ONE
      VALUE= MOLAL(5)*MOLAL(3)**2.d0
      VALUE= VALUE/(XK7*(WATER/GAMA(4))**3.d0) - ONE
  
      ! Return to calling program
      END FUNCTION FUNCB3A

!------------------------------------------------------------------------------

      SUBROUTINE CALCB3B( Y, X )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB3B
! *** CASE B3 ; SUBCASE 2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. LIQUID PHASE ONLY IS POSSIBLE
!
!     SPECIATION CALCULATIONS IS BASED ON THE HSO4 <--> SO4 EQUILIBRIUM. 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: Y, X

      ! Local variables
      INTEGER            :: K
      REAL*8             :: KK, GRAT1, DD

      !=================================================================
      ! CALCB3B begins here!
      !=================================================================

      ! Outer loop activity calculation flag
      CALAOU = .FALSE.        
      FRST   = .FALSE.
      CALAIN = .TRUE.

      ! Solve equations ; with iterations for activity coef.
      DO K = 1, NSWEEP
         GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
         DD    = SQRT( (GRAT1+Y)**2.d0 + 4.d0*(X+Y)*GRAT1)
         KK    = 0.5d0*(-(GRAT1+Y) + DD )

         ! Speciation & water content
         MOLAL (1) = KK                   ! HI
         MOLAL (5) = Y+KK                 ! SO4I
         MOLAL (6) = MAX (X+Y-KK, TINY)   ! HSO4I
         MOLAL (3) = 3.d0*Y+X              ! NH4I
         CALL CALCMR                      ! Water content

         ! Calculate activities or terminate internal loop
         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT     
      ENDDO

      ! Return to calling program
30    CONTINUE
      END SUBROUTINE CALCB3B

!------------------------------------------------------------------------------

      SUBROUTINE CALCB2
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2
! *** CASE B2 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : LC, (NH4)2SO4
!
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON THE SULFATE RATIO:
!     1. WHEN BOTH LC AND (NH4)2SO4 ARE POSSIBLE (SUBROUTINE CALCB2A)
!     2. WHEN ONLY LC IS POSSIBLE (SUBROUTINE CALCB2B).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      REAL*8 :: X, Y

      !=================================================================
      ! CALCB2 begins here!
      !=================================================================

      ! CALCULATE EQUIVALENT AMOUNT OF HSO4 AND SO4
      X = MAX(2*W(2)-W(3), TINY)   ! Equivalent NH4HSO4
      Y = MAX(W(3)  -W(2), TINY)   ! Equivalent NH42SO4

      ! CALCULATE SPECIES ACCORDING TO RELATIVE ABUNDANCE OF HSO4 *********
      IF (X.LE.Y) THEN             ! LC is the MIN (x,y)
         SCASE = 'B2 ; SUBCASE 1'
         CALL CALCB2A (X,Y-X)      ! LC + (NH4)2SO4 POSSIBLE
      ELSE
         SCASE = 'B2 ; SUBCASE 2'
         CALL CALCB2B (Y,X-Y)      ! LC ONLY POSSIBLE
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB2

!------------------------------------------------------------------------------

      SUBROUTINE CALCB2A( TLC, TNH42S4 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2
! *** CASE B2 ; SUBCASE A. 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. SOLID PHASE ONLY POSSIBLE
!     3. SOLIDS POSSIBLE: LC, (NH4)2SO4
!
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE 
!
!     FOR SOLID CALCULATIONS, A MATERIAL BALANCE BASED ON THE STOICHIMETRIC
!     PROPORTION OF AMMONIUM AND SULFATE IS DONE TO CALCULATE THE AMOUNT 
!     OF LC AND (NH4)2SO4 IN THE SOLID PHASE.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8,  INTENT(IN)  :: TLC, TNH42S4


      !=================================================================
      ! CALCB2A begins here!
      !=================================================================

      ! Regime depends upon the ambient relative humidity
      IF ( RHB .LT. DRMLCAS ) THEN    

         ! Solids possible only
         SCASE   = 'B2 ; SUBCASE A1'    
         CLC     = TLC
         CNH42S4 = TNH42S4
         SCASE   = 'B2 ; SUBCASE A1'
      ELSE

         ! Liquid & solid phase possible
         SCASE = 'B2 ; SUBCASE A2'
         CALL CALCB2A2( TLC, TNH42S4 )   
         SCASE = 'B2 ; SUBCASE A2'
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB2A

!------------------------------------------------------------------------------

      SUBROUTINE CALCB2A2 ( TLC, TNH42S4)
! 
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2A2
! *** CASE B2 ; SUBCASE A2. 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. SOLID PHASE ONLY POSSIBLE
!     3. SOLIDS POSSIBLE: LC, (NH4)2SO4
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE CALCB2A1) AND THE
!     THE SOLID WITH LIQUID PHASE (SUBROUTINE CALCB3).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: TLC, TNH42S4 

      ! Local variables
      REAL*8  :: WF, ONEMWF, CLCO, CNH42SO

      !=================================================================
      ! CALCB2A2 begins here!
      !=================================================================

      ! Find weight factor
      IF ( WFTYP .EQ. 0 ) THEN
         WF = ZERO
      ELSE IF ( WFTYP .EQ. 1 ) THEN
         WF = 0.5D0
      ELSE
         WF = ( DRLC - RHB ) / ( DRLC - DRMLCAS )
      ENDIF

      ONEMWF  = ONE - WF

      !=================================================================
      ! FIND FIRST SECTION ; DRY ONE
      !=================================================================      
      CLCO     = TLC                     ! FIRST (DRY) SOLUTION   
      CNH42SO  = TNH42S4

      !=================================================================
      ! FIND SECOND SECTION ; DRY & LIQUID 
      !=================================================================
      CLC     = ZERO
      CNH42S4 = ZERO
      CALL CALCB3                      ! SECOND (LIQUID) SOLUTION

      !=================================================================
      ! FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS
      !=================================================================
      MOLAL(1)= ONEMWF*MOLAL(1)                             ! H+
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + 3.D0*(CLCO-CLC)) ! NH4+
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               ! SO4--
      MOLAL(6)= ONEMWF*(CLCO-CLC)                                 ! HSO4-

      WATER   = ONEMWF*WATER

      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4

      ! Return to calling program
      END SUBROUTINE CALCB2A2

!------------------------------------------------------------------------------

      SUBROUTINE CALCB2B( TLC, TNH4HS4 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2
! *** CASE B2 ; SUBCASE B 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!     3. SOLIDS POSSIBLE: LC
!
!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
!     AMOUNT OF SOLID LC DISSOLVED IN THE LIQUID PHASE.
!     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB2A CALCULATES THE
!     AMOUNT OF H+ PRODUCED (BASED ON THE HSO4, SO4 RELEASED INTO THE
!     SOLUTION). THE SOLUBILITY PRODUCT OF LC IS USED AS THE OBJECTIVE 
!     FUNCTION.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 


      ! Arguments
      REAL*8  :: TLC, TNH4HS4

      ! 
      INTEGER :: K
!      REAL*8  :: FUNCB2B
      REAL*8  :: ZLO, ZHI, X1, X2, X3, Y1, Y2, Y3, YLO, YHI, DX

      !=================================================================
      ! CALCB2B begins here!
      !=================================================================
      CALAOU = .TRUE.       ! Outer loop activity calculation flag
      ZLO    = ZERO
      ZHI    = TLC          ! High limit: all of it in liquid phase

      !=================================================================
      ! INITIAL VALUES FOR BISECTION 
      !=================================================================      
      X1 = ZHI
      Y1 = FUNCB2B (X1,TNH4HS4,TLC)
      IF (ABS(Y1).LE.EPS) RETURN
      YHI= Y1                        ! Save Y-value at Hi position

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO 
      !=================================================================
      DX = (ZHI-ZLO)/NDIV
      DO K = 1, NDIV
         X2 = X1-DX
         Y2 = FUNCB2B (X2,TNH4HS4,TLC)

         ! (Y1*Y2.LT.ZERO)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  

         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================      
      ! NO SUBDIVISION WITH SOLUTION FOUND 
      !=================================================================
      YLO= Y1                        ! Save Y-value at LO position
      IF ( ABS(Y2) .LT. EPS ) THEN   ! X2 IS A SOLUTION 
         RETURN

      !=================================================================
      ! { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
      !=================================================================
      ELSE IF ( YLO .LT. ZERO .AND. YHI .LT. ZERO ) THEN
         X1 = ZHI
         X2 = ZHI
         GOTO 40

      !=================================================================
      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
      !=================================================================
      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         X1 = ZLO
         X2 = ZLO
         GOTO 40

      ELSE
         ! Warning error: no solution
         CALL PUSHERR( 0001, 'CALCB2B' )    
         RETURN

      ENDIF

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
20    CONTINUE

      DO K = 1, MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCB2B (X3,TNH4HS4,TLC)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR( 0002, 'CALCB2B' )    

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
40    CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCB2B (X3,TNH4HS4,TLC)

      ! Return to calling program
      END SUBROUTINE CALCB2B

!------------------------------------------------------------------------------

      FUNCTION FUNCB2B( X, TNH4HS4, TLC ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCB2B
! *** CASE B2 ; 
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B2 ; SUBCASE 2
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCB2B.
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: X, TNH4HS4, TLC

      ! Return value
      REAL*8             :: VALUE

      ! Local variables
      INTEGER :: K
      REAL*8  :: GRAT2, PARM, DELTA, OMEGA!, FUNCB2B

      !=================================================================
      ! FUNCB2B begins here!
      !=================================================================
      FRST   = .TRUE.
      CALAIN = .TRUE.

      ! Solve equation
      DO K = 1, NSWEEP
         GRAT2 = XK1*WATER*(GAMA(8)/GAMA(7))**2.d0/GAMA(7)
         PARM  = X+GRAT2
         DELTA = PARM*PARM + 4.d0*(X+TNH4HS4)*GRAT2 ! Diakrinousa
         OMEGA = 0.5d0*(-PARM + SQRT(DELTA))         ! Thetiki riza (ie:H+>0)

         ! Speciation & water content
         MOLAL (1) = OMEGA                         ! HI
         MOLAL (3) = 3.d0*X+TNH4HS4                 ! NH4I
         MOLAL (5) = X+OMEGA                       ! SO4I
         MOLAL (6) = MAX (X+TNH4HS4-OMEGA, TINY)   ! HSO4I
         CLC       = MAX(TLC-X,ZERO)               ! Solid LC
         CALL CALCMR                               ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 30
         ENDIF
      ENDDO
      
      !=================================================================
      ! Calculate objective function
      !=================================================================
30    CONTINUE
      !FUNCB2B= ( NH4I**3.*SO4I*HSO4I )/( XK13*(WATER/GAMA(13))**5. )
!      FUNCB2B= (MOLAL(3)**3.)*MOLAL(5)*MOLAL(6)
!      FUNCB2B= FUNCB2B/(XK13*(WATER/GAMA(13))**5.) - ONE
      VALUE= (MOLAL(3)**3.d0)*MOLAL(5)*MOLAL(6)
      VALUE= VALUE/(XK13*(WATER/GAMA(13))**5.d0) - ONE

      ! Return to calling program
      END FUNCTION FUNCB2B

!------------------------------------------------------------------------------

      SUBROUTINE CALCB1
! 
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB1
! *** CASE B1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : LC, (NH4)2SO4, NH4HSO4
!
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCB1A)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      !=================================================================
      ! CALCB1 begins here!
      !=================================================================

      ! Regime depends upon the ambient relative humidity
      IF ( RHB .LT. DRMLCAB ) THEN  

         ! Solid phase only possible
         SCASE = 'B1 ; SUBCASE 1'  
         CALL CALCB1A              
         SCASE = 'B1 ; SUBCASE 1'

      ELSE

         ! Liquid & solid phase possible
         SCASE = 'B1 ; SUBCASE 2'
         CALL CALCB1B              
         SCASE = 'B1 ; SUBCASE 2'

      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB1

!------------------------------------------------------------------------------

      SUBROUTINE CALCB1A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB1A
! *** CASE B1 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH
!     2. THERE IS NO LIQUID PHASE
!     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
!                         BUT NOT BOTH)
!
!     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE AMOUNT OF LC
!     IS CALCULATED FROM THE (NH4)2SO4 AND NH4HSO4 WHICH IS LEAST
!     ABUNDANT (STOICHIMETRICALLY). THE REMAINING EXCESS OF SALT 
!     IS MIXED WITH THE LC.  
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 


      ! Local variables
      REAL*8 :: X, Y

      !=================================================================
      ! CALCB1 begins here! 
      !=================================================================
      X = 2.d0*W(2)-W(3)       ! Equivalent NH4HSO4
      Y = W(3)-W(2)         ! Equivalent (NH4)2SO4

      ! Calculate composition
      IF (X.LE.Y) THEN      ! LC is the MIN (x,y)
         CLC     = X        ! NH4HSO4 >= (NH4)2S04
         CNH4HS4 = ZERO
         CNH42S4 = Y-X
      ELSE
         CLC     = Y        ! NH4HSO4 <  (NH4)2S04
         CNH4HS4 = X-Y
         CNH42S4 = ZERO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB1A

!------------------------------------------------------------------------------

      SUBROUTINE CALCB1B
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB1B
! *** CASE B1 ; SUBCASE 2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
!                         BUT NOT BOTH)
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE CALCB1A) AND THE
!     THE SOLID WITH LIQUID PHASE (SUBROUTINE CALCB2).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************

#     include "isoropia.h" 

      ! Local variables
      REAL*8 :: WF, ONEMWF, CLCO, CNH42SO, CNH4HSO

      !=================================================================
      ! CALCB1B begins here!
      !=================================================================

      ! Find weight factor
      IF (WFTYP.EQ.0) THEN
         WF = ZERO
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (DRNH4HS4-RHB)/(DRNH4HS4-DRMLCAB)
      ENDIF
      ONEMWF  = ONE - WF

      !=================================================================
      ! FIND FIRST SECTION ; DRY ONE
      !=================================================================
      CALL CALCB1A
      CLCO     = CLC               ! FIRST (DRY) SOLUTION
      CNH42SO  = CNH42S4
      CNH4HSO  = CNH4HS4

      !=================================================================
      ! FIND SECOND SECTION ; DRY & LIQUID
      !=================================================================
      CLC     = ZERO
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CALL CALCB2                  ! SECOND (LIQUID) SOLUTION

      !=================================================================
      ! FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
      !=================================================================
      MOLAL(1)= ONEMWF*MOLAL(1)                        ! H+
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + (CNH4HSO-CNH4HS4)  
     &                + 3.D0*(CLCO-CLC))                          ! NH4+
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               ! SO4--
      MOLAL(6)= ONEMWF*(CNH4HSO-CNH4HS4 + CLCO-CLC)               ! HSO4-

      WATER   = ONEMWF*WATER

      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4

      ! Return to calling program
      END SUBROUTINE CALCB1B

!------------------------------------------------------------------------------

      SUBROUTINE CALCC2
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCC2
! *** CASE C2 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS ONLY A LIQUID PHASE
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER :: K
      REAL*8  :: LAMDA, KAPA, PSI, PARM, BB, CC

      !=================================================================
      ! CALCC2 begins here!
      !=================================================================
      CALAOU =.TRUE.         ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
      LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
      PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION

      ! Solve equations
      DO K = 1, NSWEEP
         PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
         BB    = PSI+PARM
         CC    =-PARM*(LAMDA+PSI)
         KAPA  = 0.5d0*(-BB+SQRT(BB*BB-4.d0*CC))

         ! Speciation & water content
         MOLAL(1) = PSI+KAPA                               ! HI
         MOLAL(3) = LAMDA                                  ! NH4I
         MOLAL(5) = KAPA                                   ! SO4I
         MOLAL(6) = MAX(LAMDA+PSI-KAPA, TINY)              ! HSO4I
         CH2SO4   = MAX(MOLAL(5)+MOLAL(6)-
     &                  MOLAL(3), ZERO)  ! Free H2SO4
         CALL CALCMR                                ! Water content

         ! Calculate activities or terminate internal loop
         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT   

      ENDDO

      ! Return to calling program
30    CONTINUE
      END SUBROUTINE CALCC2

!------------------------------------------------------------------------------

      SUBROUTINE CALCC1
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCC1
! *** CASE C1 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE: NH4HSO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************!

#     include "isoropia.h"     

! Local variables
      REAL*8  :: KLO, KHI, X1, X2, X3, Y1, Y2, Y3
!      REAL*8  :: FUNCC1, YLO, YHI, DX
      REAL*8  :: YLO, YHI, DX
      INTEGER :: K

      !=================================================================
      ! CALCC1 begins here!
      !=================================================================
      CALAOU = .TRUE.    ! Outer loop activity calculation flag
      KLO    = TINY    
      KHI    = W(3)

      !=================================================================
      ! INITIAL VALUES FOR BISECTION
      !=================================================================
      X1 = KLO
      Y1 = FUNCC1(X1)
      IF ( ABS(Y1) .LE. EPS ) GOTO 50
      YLO = Y1

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX = ( KHI - KLO )/ FLOAT( NDIV )

      DO K = 1, NDIV
         X2 = X1+DX
         Y2 = FUNCC1 (X2)

         ! (Y1*Y2 .LT. ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y2) .LT. ZERO ) GOTO 20 

         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION FOUND 
      !=================================================================

      ! Save Y-value at HI position
      YHI= Y2                 

      ! X2 IS A SOLUTION 
      IF ( ABS(Y2) .LT. EPS ) THEN   
         GOTO 50

      !=================================================================
      ! { YLO, YHI } < 0.0  SOLUTION IS ALWAYS UNDERSATURATED WITH NH4HS04
      !=================================================================
      ELSE IF ( YLO .LT. ZERO .AND. YHI .LT. ZERO ) THEN
         GOTO 50

      !=================================================================
      ! { YLO, YHI } > 0.0 SOLUTION IS ALWAYS SUPERSATURATED WITH NH4HS04
      !=================================================================
      ELSE IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         X1 = KLO
         X2 = KLO
         GOTO 40
      ELSE
         CALL PUSHERR (0001, 'CALCC1')    ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF

      !=================================================================
      ! CPERFORM BISECTION OF DISSOLVED NH4HSO4
      !=================================================================
20    CONTINUE

      DO K = 1, MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCC1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR( 0002, 'CALCC1' )    

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCC1 (X3)

      ! Return to calling program
 50   CONTINUE
      END SUBROUTINE CALCC1

!------------------------------------------------------------------------------

      FUNCTION FUNCC1( KAPA ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCC1
! *** CASE C1 ; 
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE C1
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCC1.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CANREGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: KAPA

      ! Return value
      REAL*8             :: VALUE

      ! Local variables
      INTEGER            :: K
      REAL*8             :: LAMDA, PSI, PAR1, PAR2, BB, CC!, FUNCC1

      !=================================================================
      ! FUNCC1 begins here!
      !=================================================================
      FRST   = .TRUE.
      CALAIN = .TRUE.

      !==================================================================
      ! Solve equations
      !==================================================================
      PSI = W(2)-W(3)
      DO K = 1, NSWEEP
         PAR1  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
         PAR2  = XK12*(WATER/GAMA(9))**2.d0
         BB    = PSI + PAR1
         CC    =-PAR1*(PSI+KAPA)
         LAMDA = 0.5d0*(-BB+SQRT(BB*BB-4*CC))

         ! Save concentrations in molal array
         MOLAL(1) = PSI+LAMDA                    ! HI
         MOLAL(3) = KAPA                         ! NH4I
         MOLAL(5) = LAMDA                        ! SO4I
         MOLAL(6) = MAX (ZERO, PSI+KAPA-LAMDA)   ! HSO4I
         CNH4HS4  = MAX(W(3)-MOLAL(3), ZERO)  ! Solid NH4HSO4
         CH2SO4   = MAX(PSI, ZERO)               ! Free H2SO4
         CALL CALCMR                             ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 30
         ENDIF
      ENDDO

      !==================================================================
      ! Calculate zero function
      !==================================================================
 30   CONTINUE
      !###FUNCC1= (NH4I*HSO4I/PAR2) - ONE
!      FUNCC1= (MOLAL(3)*MOLAL(6)/PAR2) - ONE
      VALUE= (MOLAL(3)*MOLAL(6)/PAR2) - ONE

      ! Return to calling program
      END FUNCTION FUNCC1

!------------------------------------------------------------------------------

      SUBROUTINE CALCD3
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCD3
! *** CASE D3
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. THERE IS OLNY A LIQUID PHASE
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER :: K
      REAL*8  :: PSI4LO, PSI4HI, X1, X2, X3, Y1, Y2, Y3, YLO, YHI, DX
      REAL*8  :: P4, YY, DELTA
!      REAL*8  :: FUNCD3
      REAL*8      :: CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )
 
      !=================================================================
      ! CALCD3 begins here!
      !=================================================================

      ! Find dry composition
      CALL CALCD1A

      ! Setup parameters
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3

      ! Assign initial PSI's
      PSI1 = CNH4NO3               
      PSI2 = CHI2
      PSI3 = ZERO   
      PSI4 = ZERO  

      MOLAL(5) = ZERO
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1

      ! Initial water
      CALL CALCMR              

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit

      !=================================================================
      ! INITIAL VALUES FOR BISECTION
      !=================================================================
 60   CONTINUE
      X1 = PSI4LO
      Y1 = FUNCD3( X1 )
      IF ( ABS(Y1) .LE. EPS ) RETURN
      YLO= Y1                 ! Save Y-value at HI position

      !=================================================================     
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX = ( PSI4HI - PSI4LO ) / FLOAT( NDIV )

      DO K = 1, NDIV
         X2 = X1+DX
         Y2 = FUNCD3 (X2)

         ! (Y1*Y2.LT.ZERO)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  

         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================     
      ! NO SUBDIVISION WITH SOLUTION FOUND 
      !=================================================================
      YHI= Y1                        ! Save Y-value at Hi position
      IF ( ABS(Y2) .LT. EPS ) THEN   ! X2 IS A SOLUTION 
         RETURN

      !=================================================================
      ! { YLO, YHI } < 0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
      !
      ! Physically I dont know when this might happen, but I have put 
      ! this branch in for completeness. I assume there is no solution; 
      ! all NO3 goes to the gas phase. (rjp)
      !=================================================================
      ELSE IF ( YLO .LT. ZERO .AND. YHI .LT. ZERO ) THEN
         P4 = TINY ! PSI4LO ! CHI4
         YY = FUNCD3(P4)
         GOTO 50

      !=================================================================
      ! { YLO, YHI } > 0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
      !
      ! This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate 
      ! evaporates and goes to the gas phase ; so I redefine the LO and 
      ! HI limits of PSI4 and proceed again with root tracking. (rjp)
      !=================================================================
      ELSE IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN

         PSI4HI = PSI4LO

         ! No solution; some NH3 evaporates
         PSI4LO = PSI4LO - 0.1d0*(PSI1+PSI2) 

         IF ( PSI4LO .LT. -( PSI1 + PSI2 ) ) THEN

            ! Warning error: no solution
            CALL PUSHERR( 0001, 'CALCD3' )  
            RETURN
         ELSE

            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  ! Initial water
            GOTO 60                      ! Redo root tracking
         ENDIF
      ENDIF

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
20    CONTINUE

      DO K = 1, MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCD3 (X3)

         ! (Y1*Y3 .LE. ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y3) .LE. ZERO ) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS( X2 - X1 ) .LE. EPS * ABS( X1 ) ) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR (0002, 'CALCD3')    

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================    
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCD3( X3 )

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN
      !================================================================= 
 50   CONTINUE

      IF ( MOLAL(1) .GT. TINY ) THEN
         CALL CALCHS4( MOLAL(1), MOLAL(5), 
     &                 ZERO, DELTA )
         MOLAL(1) = MOLAL(1) - DELTA            ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA            ! SO4  EFFECT
         MOLAL(6) = DELTA                       ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCD3

!------------------------------------------------------------------------------

      FUNCTION FUNCD3( P4 ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCD3
! *** CASE D3 
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ; 
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD3.
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: P4

      ! Return value
      REAL*8             :: VALUE

      ! Local variables
      INTEGER :: K
      REAL*8  :: BB, DENM, ABB, AHI!, FUNCD3
      REAL*8      :: CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCD3 begins here!
      !=================================================================
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4

      !=================================================================
      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF
      !=================================================================
      DO K = 1, NSWEEP
         A2   = XK7*(WATER/GAMA(4))**3.d0
         A3   = XK4*R*TEMP*(WATER/GAMA(10))**2.d0
         A4   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A7   = XKW *RHB*WATER*WATER

         PSI3 = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3 = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4) 
         PSI3 = MIN(MAX(PSI3, ZERO), CHI3)

         BB   = PSI4 - PSI3
!###old         AHI  = 0.5*(-BB + SQRT(BB*BB + 4.d0*A7)) ! This is correct also
!###        AHI  =2.0*A7/(BB+SQRT(BB*BB + 4.d0*A7)) ! Avoid overflow when HI->0

         DENM = BB+SQRT(BB*BB + 4.d0*A7)
         IF (DENM.LE.TINY) THEN       ! Avoid overflow when HI->0
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.d0*A7/ABB ! Taylor expansion of SQRT
         ENDIF
         AHI = 2.d0*A7/DENM

         ! SPECIATION & WATER CONTENT
         MOLAL (1) = AHI                             ! HI
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2         ! NH4I
         MOLAL (5) = PSI2                            ! SO4I
         MOLAL (6) = ZERO                            ! HSO4I
         MOLAL (7) = PSI3 + PSI1                     ! NO3I
         CNH42S4   = CHI2 - PSI2                     ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                            ! Solid NH4NO3
         GHNO3     = CHI3 - PSI3                     ! Gas HNO3
         GNH3      = CHI4 - PSI4                     ! Gas NH3
         CALL CALCMR                                 ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND.CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
      ENDDO

      !=================================================================
      ! Calculate objective function
      !=================================================================
20    CONTINUE
      VALUE= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 

      ! Return to calling program
      END FUNCTION FUNCD3

!------------------------------------------------------------------------------

      SUBROUTINE CALCD2
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCD2
! *** CASE D2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER :: K
      REAL*8  :: PSI4LO, PSI4HI, X1, X2, X3, Y1, Y2, Y3, YLO, YHI, DX
!      REAL*8  :: FUNCD2, P4, YY, DELTA
      REAL*8  :: P4, YY, DELTA
      REAL*8      :: CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCD2 begins here!
      !=================================================================

      ! Find dry composition
      CALL CALCD1A

      ! Setup parameters
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3

      PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
      PSI2 = CNH42S4
      PSI3 = ZERO   
      PSI4 = ZERO  

      MOLAL(5) = ZERO
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1
      CALL CALCMR                  ! Initial water

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit

      !=================================================================
      ! Initial values for bisection
      !=================================================================
 60   CONTINUE
      X1 = PSI4LO
      Y1 = FUNCD2 (X1)
      IF ( ABS( Y1 ) .LE. EPS ) RETURN

      ! Save Y-value at HI position
      YLO = Y1                  

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX   = (PSI4HI-PSI4LO)/FLOAT(NDIV)

      DO K = 1, NDIV
         X2 = X1+DX
         Y2 = FUNCD2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) THEN

            ! This is done, in case if Y(PSI4LO)>0, 
            ! but Y(PSI4LO+DX) < 0 (i.e.undersat)
            IF (Y1 .LE. Y2) GOTO 20 ! (Y1*Y2.LT.ZERO)
         ENDIF
         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================
      ! No subdivision with solution found
      !=================================================================
      YHI= Y1                          ! Save Y-value at Hi position
      IF ( ABS( Y2 ) .LT. EPS ) THEN   ! X2 IS A SOLUTION 
         RETURN

      !=================================================================
      ! { YLO, YHI } < 0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
      !
      ! Physically I dont know when this might happen, but I have put 
      ! this branch in for completeness. I assume there is no solution; 
      ! all NO3 goes to the gas phase. (rjp)
      !=================================================================
      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = TINY ! PSI4LO ! CHI4
         YY = FUNCD2(P4)
         GOTO 50

      !=================================================================
      ! { YLO, YHI } > 0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
      ! 
      ! This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate 
      ! evaporatesand goes to the gas phase ; so I redefine the LO and 
      ! HI limits of PSI4 and proceed again with root tracking. (rjp)
      !=================================================================
      ELSE IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         PSI4HI = PSI4LO

         ! No solution; some NH3 evaporates
         PSI4LO = PSI4LO - 0.1d0*(PSI1+PSI2) 

         IF ( PSI4LO .LT. -( PSI1 + PSI2 ) ) THEN

            ! Warning error: no solution
            CALL PUSHERR( 0001, 'CALCD2' )  
            RETURN
         ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  ! Initial water
            GOTO 60                      ! Redo root tracking
         ENDIF
      ENDIF

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
 20   CONTINUE

      DO K = 1, MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCD2 (X3)

         ! (Y1*Y3 .LE. ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y3) .LE. ZERO ) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR (0002, 'CALCD2')    

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
 40   CONTINUE

      ! 0.5*(X1+X2)  ! Get "low" side, it's acidic soln.
      X3 = MIN(X1,X2)           
      Y3 = FUNCD2 (X3)

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN
      !=================================================================
 50   CONTINUE
      
      IF ( MOLAL(1) .GT. TINY ) THEN
         CALL CALCHS4( MOLAL(1), MOLAL(5), 
     &                 ZERO, DELTA )
         MOLAL(1) = MOLAL(1) - DELTA          ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA         ! SO4  EFFECT
         MOLAL(6) = DELTA                          ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCD2

!------------------------------------------------------------------------------
      
      FUNCTION FUNCD2( P4 ) RESULT( VALUE )
! 
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCD2
! *** CASE D2 
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D2 ; 
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD2.
!
!******************************************************************************
!
#     include "isoropia.h" 

      REAL*8, INTENT(IN)  :: P4

      ! Return value
      REAL*8             :: VALUE

      ! Local variables
      INTEGER            :: K, ISLV
      REAL*8             :: PSI14, BB, DENM, ABB, AHI!, FUNCD2
      REAL*8      :: CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCD2 begins here!
      !=================================================================

      ! Reset activity coefficients to 0.1
      CALL RSTGAM       
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4
      PSI2   = CHI2

      !=================================================================       
      ! Solve equations ; with iterations for activity coef.
      !=================================================================
      DO K = 1, NSWEEP
         A2  = XK7*(WATER/GAMA(4))**3.d0
         A3  = XK4*R*TEMP*(WATER/GAMA(10))**2.d0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A7  = XKW *RHB*WATER*WATER

         IF ( CHI2 .GT. TINY .AND. WATER .GT. TINY ) THEN
            PSI14 = PSI1+PSI4

            CALL POLY3 (PSI14,0.25*PSI14**2.,-A2/4.D0, PSI2, ISLV)  ! PSI2
            IF ( ISLV.EQ.0 ) THEN
                PSI2 = MIN (PSI2, CHI2)
            ELSE
                PSI2 = ZERO
            ENDIF
         ENDIF

         PSI3  = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3  = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4) 
         !### PSI3  = MIN(MAX(PSI3, ZERO), CHI3)

         ! (BB > 0, acidic solution, <0 alkaline)
         BB   = PSI4 - PSI3 

         !==============================================================
         ! Do not change computation scheme for H+, 
         ! all others did not work well
         !==============================================================
         DENM = BB+SQRT(BB*BB + 4.d0*A7)

         ! Avoid overflow when HI->0
         IF ( DENM .LE. TINY ) THEN       
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.d0*A7/ABB ! Taylor expansion of SQRT
         ENDIF
         AHI = 2.d0*A7/DENM

         ! Speciation & water content
         MOLAL (1) = AHI                              ! HI
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2          ! NH4
         MOLAL (5) = PSI2                             ! SO4
         MOLAL (6) = ZERO                             ! HSO4
         MOLAL (7) = PSI3 + PSI1                      ! NO3
         CNH42S4   = CHI2 - PSI2                      ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                             ! Solid NH4NO3
         GHNO3     = CHI3 - PSI3                      ! Gas HNO3
         GNH3      = CHI4 - PSI4                      ! Gas NH3
         CALL CALCMR                                  ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
      ENDDO

      !=================================================================
      ! Calculate objective function
      !=================================================================
 20   CONTINUE
      VALUE= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 

      ! Return to calling program
      END FUNCTION FUNCD2

!------------------------------------------------------------------------------

      SUBROUTINE CALCD1
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCD1
! *** CASE D1 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!
!     THERE ARE TWO REGIMES DEFINED BY RELATIVE HUMIDITY:
!     1. RH < MDRH ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCD1A)
!     2. RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 



      !=================================================================
      ! CALCD1 begins here!
      !=================================================================
      IF ( RHB .LT. DRMASAN ) THEN    

         ! Solid phase only possible
         SCASE = 'D1 ; SUBCASE 1'   
         CALL CALCD1A            
         SCASE = 'D1 ; SUBCASE 1'

      ELSE

         ! Liquid & solid phase possible
         SCASE = 'D1 ; SUBCASE 2'   
         CALL CALCMDRH (RHB, DRMASAN, DRNH4NO3, 
     &                  CALCD1A, CALCD2)
         SCASE = 'D1 ; SUBCASE 2'

      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCD1

!------------------------------------------------------------------------------

      SUBROUTINE CALCD1A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCD1A
! *** CASE D1 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!
!     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
!     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
!     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
!     THE SOLID PHASE.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" ! CNH42S4, CNH4NO3, GNH3, GHNO3

      ! Local variables
      REAL*8  :: PARM, X, PS, OM, OMPS, DIAK, ZE

      !=================================================================
      ! CALCD1A begins here!
      !=================================================================
      PARM    = XK10/(R*TEMP)/(R*TEMP)

      ! Calculate NH4NO3 that volatizes
      CNH42S4 = W(2)                                    
      X       = MAX(ZERO, MIN(W(3)-2.d0*CNH42S4, W(4))) ! MAX NH4NO3
      PS      = MAX(W(3) - X - 2.d0*CNH42S4, ZERO)
      OM      = MAX(W(4) - X, ZERO)

      OMPS    = OM+PS
      DIAK    = SQRT(OMPS*OMPS + 4.d0*PARM)              ! DIAKRINOUSA
      ZE      = MIN(X, 0.5d0*(-OMPS + DIAK))              ! THETIKI RIZA

      ! Speciation
      CNH4NO3 = X  - ZE    ! Solid NH4NO3
      GNH3    = PS + ZE    ! Gas NH3
      GHNO3   = OM + ZE    ! Gas HNO3

      ! Return to calling program
      END SUBROUTINE CALCD1A

!------------------------------------------------------------------------------

      SUBROUTINE CALCG1
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCG1
! *** CASE G1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4CL, NA2SO4
!
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCG1A)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 


      !=================================================================
      ! CALCG1 begins here!
      !=================================================================

      ! Regime depends upon the ambient relative humidity
      IF (RHB.LT.DRMG1) THEN    
         SCASE = 'G1 ; SUBCASE 1'  
         CALL CALCG1A              ! SOLID PHASE ONLY POSSIBLE
         SCASE = 'G1 ; SUBCASE 1'
      ELSE
         SCASE = 'G1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
         CALL CALCMDRH (RHB, DRMG1, DRNH4NO3, CALCG1A, 
     &                  CALCG2A)
         SCASE = 'G1 ; SUBCASE 2'
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCG1

!------------------------------------------------------------------------------

      SUBROUTINE CALCG1A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCG1A
! *** CASE G1 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!
!     SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
!     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
!     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
!     THE SOLID PHASE.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      REAL*8  :: LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2
      REAL*8  :: ALF, BET, GAM, RTSQ, THETA1, THETA2, BB, CC
      REAL*8  :: DD, DQDD, DD1, DD2, SQDD1, SQDD2, SQDD
      REAL*8  :: A1, A2

      !=================================================================
      ! CALCG1A begins here!
      !=================================================================

      ! Calculate non-volatile solids
      CNA2SO4 = 0.5d0*W(1)
      CNH42S4 = W(2) - CNA2SO4

      ! Calculate volatile species
      ALF     = W(3) - 2.d0*CNH42S4
      BET     = W(5)
      GAM     = W(4)

      RTSQ    = R*TEMP*R*TEMP
      A1      = XK6/RTSQ
      A2      = XK10/RTSQ

      THETA1  = GAM - BET*(A2/A1)
      THETA2  = A2/A1

      ! Quadratic equation solution
      BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
      CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
      DD      = BB*BB - 4.0D0*CC
      IF (DD.LT.ZERO) GOTO 100   ! Solve each reaction seperately

      ! Two roots for KAPA, check and see if any valid
      SQDD    = SQRT(DD)
      KAPA1   = 0.5D0*(-BB+SQDD)
      KAPA2   = 0.5D0*(-BB-SQDD)
      LAMDA1  = THETA1 + THETA2*KAPA1
      LAMDA2  = THETA1 + THETA2*KAPA2

      IF (KAPA1.GE.ZERO .AND. LAMDA1.GE.ZERO) THEN
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND.
     &       BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF

      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND. 
     &       BET-KAPA2.GE.ZERO .AND. GAM-LAMDA2.GE.ZERO) THEN
             KAPA = KAPA2
             LAMDA= LAMDA2
             GOTO 200
         ENDIF
      ENDIF

      ! Separate solution of NH4Cl and NH4NO3 equilibria
100   KAPA  = ZERO
      LAMDA = ZERO
      DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
      DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)

      ! NH4Cl equilibrium
      IF (DD1.GE.ZERO) THEN
         SQDD1 = SQRT(DD1)
         KAPA1 = 0.5D0*(ALF+BET + SQDD1)
         KAPA2 = 0.5D0*(ALF+BET - SQDD1)

         IF (KAPA1.GE.ZERO .AND. KAPA1.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA1 
         ELSE IF (KAPA2.GE.ZERO .AND. KAPA2.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA2
         ELSE
            KAPA = ZERO
         ENDIF
      ENDIF

      ! NH4NO3 equilibrium
      IF (DD2.GE.ZERO) THEN
         SQDD2 = SQRT(DD2)
         LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
         LAMDA2= 0.5D0*(ALF+GAM - SQDD2)

         IF (LAMDA1.GE.ZERO .AND. LAMDA1.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1 
         ELSE IF (LAMDA2.GE.ZERO .AND. LAMDA2.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
         ELSE
            LAMDA = ZERO
         ENDIF
      ENDIF

      ! If both KAPA and LAMDA are > 0, then apply the existance criterion
      IF (KAPA.GT.ZERO .AND. LAMDA.GT.ZERO) THEN
         IF (BET .LT. LAMDA/THETA1) THEN
            KAPA = ZERO
         ELSE
            LAMDA= ZERO
         ENDIF
      ENDIF

      ! Calculate composition of volatile species
200   CONTINUE
      CNH4NO3 = LAMDA
      CNH4CL  = KAPA

      GNH3    = MAX(ALF - KAPA - LAMDA, ZERO)
      GHNO3   = MAX(GAM - LAMDA, ZERO)
      GHCL    = MAX(BET - KAPA, ZERO)

      ! Return to calling program
      END SUBROUTINE CALCG1A

!------------------------------------------------------------------------------

      SUBROUTINE CALCG2
!
!******************************************************************************
! *** ISORROPIA CODE
! *** SUBROUTINE CALCG2
! *** CASE G2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. LIQUID & SOLID PHASE ARE BOTH POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      INTEGER  :: K

      !=================================================================
      ! CALCG2 begins here!
      !=================================================================

      ! Regime depends on the existance of nitrates 
      IF (W(4).GT.TINY) THEN        ! NO3 EXISTS, WATER POSSIBLE
         SCASE = 'G2 ; SUBCASE 1'  
         CALL CALCG2A
         SCASE = 'G2 ; SUBCASE 1' 
      ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
         SCASE = 'G1 ; SUBCASE 1'  
         CALL CALCG1A
         SCASE = 'G1 ; SUBCASE 1'  
      ENDIF

      ! Regime depends on the existance of water and of the RH 
      IF (WATER.LE.TINY) THEN
         IF (RHB.LT.DRMG2) THEN             ! ONLY SOLIDS 
            WATER = TINY
            DO 10 K=1,NIONS
               MOLAL(K) = ZERO
10          CONTINUE
            CALL CALCG1A
            SCASE = 'G2 ; SUBCASE 2'  
         ELSE
            IF (W(5).GT. TINY) THEN
               SCASE = 'G2 ; SUBCASE 3'    ! MDRH (NH4CL, NA2SO4, NH42S4)  
               CALL CALCMDRH (RHB, DRMG2, DRNH4CL, CALCG1A, 
     &                        CALCG3A)
               SCASE = 'G2 ; SUBCASE 3'  
            ENDIF
            IF (WATER.LE.TINY .AND. RHB.GE.DRMG3) THEN
               SCASE = 'G2 ; SUBCASE 4'    ! MDRH (NA2SO4, NH42S4)
               CALL CALCMDRH (RHB, DRMG3, DRNH42S4, 
     &                        CALCG1A, CALCG4)
               SCASE = 'G2 ; SUBCASE 4'  
            ELSE
               WATER = TINY
               DO 20 K=1,NIONS
                  MOLAL(K) = ZERO
20             CONTINUE
               CALL CALCG1A
               SCASE = 'G2 ; SUBCASE 2'  
            ENDIF
         ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCG2

!------------------------------------------------------------------------------

      SUBROUTINE CALCG2A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCG2A
! *** CASE G2 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" ! ISOROPIA common blocks

      ! Local variables
      INTEGER ::     ISLV, K
      REAL*8  ::     PSI6LO, PSI6HI, X1, Y1, DX, X2, Y2
      REAL*8  ::     X3, Y3, DELTA

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
!$OMP THREADPRIVATE( /CASEG/ )

      !=================================================================
      ! CALCG2A begins here!
      !=================================================================

      CALAOU = .TRUE.   
      CHI1   = 0.5d0*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)

      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY

      WATER  = TINY

      ! Initial values for bisection
      X1 = PSI6LO
      Y1 = FUNCG2A (X1)
      IF (CHI6.LE.TINY) GOTO 50  
      ! Comment out
      !IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
      !IF (WATER .LE. TINY) GOTO 50               ! No water

      !=================================================================
      ! Root tracking ; for the range of HI and LO
      !=================================================================
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCG2A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      !=================================================================
      ! No subdivision with solution; If ABS(Y2)<EPS solution is assumed
      !=================================================================
      IF (ABS(Y2) .GT. EPS) WATER = TINY
      GOTO 50

      !=================================================================
      ! Perform bisection
      !=================================================================
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCG2A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCG2A') ! WARNING ERROR: NO CONVERGENCE

      !=================================================================
      ! Converged ; Return 
      !=================================================================
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      IF (X3.LE.TINY2) THEN   ! PRACTICALLY NO NITRATES, SO DRY SOLUTION
         WATER = TINY
      ELSE
         Y3 = FUNCG2A (X3)
      ENDIF

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN 
      !=================================================================
 50   CONTINUE

      ! Na2SO4 dissolution
      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        ! PSI1
         CALL POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ELSE
         PSI1 = ZERO
      ENDIF
      MOLAL(2) = 2.0D0*PSI1               ! Na+  EFFECT
      MOLAL(5) = MOLAL(5) + PSI1          ! SO4  EFFECT
      CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   ! NA2SO4(s) depletion

      ! HSO4 equilibrium
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), 
     &                 ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA     ! H+   AFFECT
         MOLAL(5) = MOLAL(5) - DELTA     ! SO4  AFFECT
         MOLAL(6) = DELTA                ! HSO4 AFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCG2A

!------------------------------------------------------------------------------

      FUNCTION FUNCG2A( X ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCG2A
! *** CASE G2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" ! ISOROPIA common blocks

      ! Arguments
      REAL*8, INTENT(IN)  :: X

      ! Local variables
      INTEGER ::     K, ISLV
      REAL*8  ::     VALUE
      REAL*8  ::     DENO, PSI20, SMIN, HI, OHI, DELT
      REAL*8  ::     BB, CC, DD, PSI31, PSI32

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
!$OMP THREADPRIVATE( /CASEG/ )

      !=================================================================
      ! FUNCG2A begins here!
      !=================================================================
      PSI6   = X
      PSI2   = CHI2
      PSI3   = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE. 

      !=================================================================
      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      !=================================================================
      DO 10 K=1,NSWEEP

         A1  = XK5 *(WATER/GAMA(2))**3.d0
         A2  = XK7 *(WATER/GAMA(4))**3.d0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.d0
         A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.d0
     
         DENO = MAX(CHI6-PSI6-PSI3, ZERO)
         PSI5 = CHI5/((A6/A5)*(DENO/PSI6) + ONE)

         PSI4 = MIN(PSI5+PSI6,CHI4)

         IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN     
            CALL POLY3 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
            IF (ISLV.EQ.0) PSI2 = MIN (PSI20, CHI2)
         ENDIF

         !==============================================================
         ! Save concentrations in MOLAL array
         !==============================================================
         MOLAL (2) = ZERO              ! NA
         MOLAL (3) = 2.d0*PSI2 + PSI4  ! NH4I
         MOLAL (4) = PSI6              ! CLI
         MOLAL (5) = PSI2              ! SO4I
         MOLAL (6) = ZERO              ! HSO4
         MOLAL (7) = PSI5              ! NO3I
         ! Comment out
         !MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
         SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+
     &               MOLAL(4)-MOLAL(2)-MOLAL(3)
         CALL CALCPH (SMIN, HI, OHI)
         MOLAL (1) = HI

         !==============================================================
         ! Calculate gas / solid species (liquid in MOLAL already)
         !==============================================================
         GNH3      = MAX(CHI4 - PSI4, TINY)
         GHNO3     = MAX(CHI5 - PSI5, TINY)
         GHCL      = MAX(CHI6 - PSI6, TINY)

         CNH42S4   = MAX(CHI2 - PSI2, ZERO)
         CNH4NO3   = ZERO

         !==============================================================
         ! NH4Cl(s) calculations
         !==============================================================
         A3   = XK6 /(R*TEMP*R*TEMP)
         IF (GNH3*GHCL.GT.A3) THEN
            DELT = MIN(GNH3, GHCL)
            BB = -(GNH3+GHCL)
            CC = GNH3*GHCL-A3
            DD = BB*BB - 4.D0*CC
            PSI31 = 0.5D0*(-BB + SQRT(DD))
            PSI32 = 0.5D0*(-BB - SQRT(DD))
            IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
               PSI3 = PSI31
            ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
               PSI3 = PSI32
            ELSE
               PSI3 = ZERO
            ENDIF
         ELSE
            PSI3 = ZERO
         ENDIF

         !=============================================================
         ! Calculate gas / solid species (liquid in MOLAL already)
         !=============================================================
         GNH3    = MAX(GNH3 - PSI3, TINY)
         GHCL    = MAX(GHCL - PSI3, TINY)
         CNH4CL  = PSI3

         !==============================================================
         ! Calculate MOLALR array, water and activities
         !==============================================================
         CALL CALCMR

         !==============================================================
         ! Calculate activities or terminate internal loop 
         !==============================================================
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      !=================================================================
      ! Calculate function value for outer loop 
      !=================================================================
 20   CONTINUE
      IF (CHI4.LE.TINY) THEN
         VALUE = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
      ELSE
         VALUE = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
      ENDIF

      ! Return to calling program
      END FUNCTION FUNCG2A

!------------------------------------------------------------------------------

      SUBROUTINE CALCG3
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCG3
! *** CASE G3
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. LIQUID & SOLID PHASE ARE BOTH POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h"  ! ISOROPIA common blocks

      ! Local variables
      INTEGER  :: K

      !=================================================================
      ! CALCG3 begins here!
      !=================================================================

      ! Regime depends on the existance of water and of the RH
      IF (W(4).GT.TINY .AND. W(5).GT.TINY) THEN 

         ! NO3,CL EXIST, WATER POSSIBLE
         SCASE = 'G3 ; SUBCASE 1'  
         CALL CALCG3A
         SCASE = 'G3 ; SUBCASE 1' 

      ELSE                               

         ! NO3, CL NON EXISTANT
         SCASE = 'G1 ; SUBCASE 1'  
         CALL CALCG1A
         SCASE = 'G1 ; SUBCASE 1'  

      ENDIF

      IF (WATER.LE.TINY) THEN

         IF (RHB.LT.DRMG3) THEN          

            ! ONLY SOLIDS 
            WATER = TINY
            DO 10 K=1,NIONS
               MOLAL(K) = ZERO
 10         CONTINUE
            CALL CALCG1A
            SCASE = 'G3 ; SUBCASE 2'  
            RETURN

         ELSE

            ! MDRH REGION (NA2SO4, NH42S4)  
            SCASE = 'G3 ; SUBCASE 3' 
            CALL CALCMDRH (RHB, DRMG3, DRNH42S4, CALCG1A, CALCG4)
            SCASE = 'G3 ; SUBCASE 3'  
 
         ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCG3

!------------------------------------------------------------------------------

      SUBROUTINE CALCG3A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCG3A
! *** CASE G3 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K, ISLV
      REAL*8  ::     PSI6LO, PSI6HI, X1, Y1, DX
      REAL*8  ::     X2, Y2, X3, Y3, DELTA

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
!$OMP THREADPRIVATE( /CASEG/ )

      !=================================================================
      ! CALCG3A begins here!
      !=================================================================

      CALAOU = .TRUE.   
      CHI1   = 0.5d0*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
 
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      WATER  = TINY

      !=================================================================
      ! INITIAL VALUES FOR BISECTION
      !=================================================================
      X1 = PSI6LO
      Y1 = FUNCG3A (X1)
      IF (CHI6.LE.TINY) GOTO 50  
      ! comment out
      !IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY .OR. WATER .LE. TINY) GOTO 50
      !IF (WATER .LE. TINY) RETURN                    ! No water

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2  = X1+DX 
         Y2  = FUNCG3A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1  = X2
         Y1  = Y2
 10   CONTINUE

      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
      !=================================================================
      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG3A (PSI6LO)
      GOTO 50

      !=================================================================
      ! PERFORM BISECTION 
      !=================================================================
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCG3A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCG3A') ! WARNING ERROR: NO CONVERGENCE

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCG3A (X3)
      
      !=================================================================
      ! FINAL CALCULATIONS
      !=================================================================
 50   CONTINUE

      ! Na2SO4 DISSOLUTION
      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN ! PSI1
         CALL POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI1 = MIN (PSI1, CHI1)
         ELSE
            PSI1 = ZERO
         ENDIF
      ELSE
         PSI1 = ZERO
      ENDIF
      MOLAL(2) = 2.0D0*PSI1     ! Na+  EFFECT
      MOLAL(5) = MOLAL(5) + PSI1 ! SO4  EFFECT
      CNA2SO4  = MAX(CHI1 - PSI1, ZERO) ! NA2SO4(s) depletion

      ! HSO4 equilibrium
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA         ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA         ! SO4  EFFECT
         MOLAL(6) = DELTA                    ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCG3A

!------------------------------------------------------------------------------

      FUNCTION FUNCG3A( X ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCG3A
! *** CASE G3 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h"  ! ISOROPIA common blocks

      ! Arguments
      REAL*8, INTENT(IN) :: X

      ! Local variables
      INTEGER ::     K, ISLV
      REAL*8  ::     BB, CC, DD, PSI20, SMIN, HI, OHI
      REAL*8  ::     VALUE

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
!$OMP THREADPRIVATE( /CASEG/ )

      !=================================================================
      ! FUNCG3A begins here!
      !=================================================================
      PSI6   = X
      PSI2   = CHI2
      FRST   = .TRUE.
      CALAIN = .TRUE. 

      !=================================================================       
      ! Solve equations ; with iterations for activity coef.
      !=================================================================
      DO 10 K=1,NSWEEP

      A1  = XK5 *(WATER/GAMA(2))**3.d0
      A2  = XK7 *(WATER/GAMA(4))**3.d0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.d0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.d0

      !=================================================================       
      ! Calculate dissociation quantities
      !=================================================================
      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF

      !IF(CHI4.GT.TINY) THEN
      IF(W(2).GT.TINY) THEN       
     
         ! Accounts for NH3 evaporation
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)  

         ! Patch proposed by Uma Shankar, 19/11/01
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF

      IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN     
         CALL POLY3 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
         IF (ISLV.EQ.0) PSI2 = MIN (PSI20, CHI2)
      ENDIF

      !=================================================================
      ! CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY)
      !=================================================================
      MOLAL (2) = ZERO                                ! Na
      MOLAL (3) = 2.d0*PSI2 + PSI4                    ! NH4I
      MOLAL (4) = PSI6                                ! CLI
      MOLAL (5) = PSI2                                ! SO4I
      MOLAL (6) = ZERO                                ! HSO4
      MOLAL (7) = PSI5                                ! NO3I

      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+
     &            MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
      GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3
      GHCL      = MAX(CHI6 - PSI6, TINY)              ! Gas HCl

      CNH42S4   = CHI2 - PSI2                         ! Solid (NH4)2SO4
      CNH4NO3   = ZERO                                ! Solid NH4NO3
      CNH4CL    = ZERO                                ! Solid NH4Cl

      CALL CALCMR                                     ! Water content

      !=================================================================
      ! Calculate activities or terminate internal loop
      !=================================================================
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT     
      ELSE
         GOTO 20
      ENDIF
 10   CONTINUE

      !=================================================================
      ! Calculate function value for outer loop
      !=================================================================
 20   CONTINUE
      VALUE = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      ! Return to calling program
      END FUNCTION FUNCG3A

!------------------------------------------------------------------------------

      SUBROUTINE CALCG4
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCG4
! *** CASE G4
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h"  ! ISOROPIA common blocks

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     PSI6LO, PSI6HI, X1, Y1, DX
      REAL*8  ::     X2, Y2, X3, Y3, DELTA

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
!$OMP THREADPRIVATE( /CASEG/ )

      !=================================================================
      ! CALCG4 begins here!
      !=================================================================

      CALAOU = .TRUE.   
      CHI1   = 0.5d0*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)

      PSI2   = CHI2
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      WATER  = CHI2/M0(4) + CHI1/M0(2)

      !=================================================================
      ! Initial values for bisection 
      !=================================================================
      X1 = PSI6LO
      Y1 = FUNCG4A (X1)
      IF (CHI6.LE.TINY) GOTO 50  
      ! Comment out
      !IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY .OR. WATER .LE. TINY) GOTO 50
      !IF (WATER .LE. TINY) RETURN                    ! No water

      !=================================================================
      ! Root tracking ; for the range of HI and LO
      !=================================================================
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2  = X1+DX
         Y2  = FUNCG4A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1  = X2
         Y1  = Y2
 10   CONTINUE

      !=================================================================  
      ! No subdivision with solution; IF ABS(Y2)<EPS solution is assumed
      !=================================================================  
      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG4A (PSI6LO)
      GOTO 50

      !=================================================================  
      ! Perform bisection 
      !=================================================================  
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCG4A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCG4') ! WARNING ERROR: NO CONVERGENCE

      !================================================================= 
      ! Converged ; return 
      !================================================================= 
40    CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCG4A (X3)

      !================================================================= 
      ! CALCULATE HSO4 SPECIATION AND RETURN
      !================================================================= 
 50   CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA         ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA         ! SO4  EFFECT
         MOLAL(6) = DELTA                    ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCG4

!------------------------------------------------------------------------------

      FUNCTION FUNCG4A( X ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCG4A
! *** CASE G4
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h"  ! ISOROPIA common blocks

      ! Arguments
      REAL*8, INTENT(IN) :: X

      ! Local variables
      INTEGER ::     K, ISLV
      REAL*8  ::     BB, CC, DD, CLI, SO4I, HI, OHI
      REAL*8  ::     NAI, NH4I, NO3I, VALUE

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
!$OMP THREADPRIVATE( /CASEG/ )

      !=================================================================
      ! FUNCG4A begins here!
      !=================================================================
      PSI6   = X
      PSI1   = CHI1
      FRST   = .TRUE.
      CALAIN = .TRUE. 

      !=================================================================
      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      !=================================================================
      DO 10 K=1,NSWEEP

         A1  = XK5 *(WATER/GAMA(2))**3.d0
         A2  = XK7 *(WATER/GAMA(4))**3.d0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.d0
         A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.d0

         ! CALCULATE DISSOCIATION QUANTITIES
         IF (CHI5.GE.TINY) THEN
            PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
         ELSE
            PSI5 = TINY
         ENDIF

         ! comment out
         !IF(CHI4.GT.TINY) THEN
         IF(W(2).GT.TINY) THEN       

            ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO) 

            ! Patch proposed by Uma shankar, 19/11/2001
            PSI4 =0.5d0*(-BB - SQRT(DD))
         ELSE
            PSI4 = TINY
         ENDIF

         !==============================================================
         ! CALCULATE CONCENTRATIONS
         !==============================================================
         NH4I = 2.d0*PSI2 + PSI4
         CLI  = PSI6
         SO4I = PSI2 + PSI1
         NO3I = PSI5
         NAI  = 2.0D0*PSI1  

         CALL CALCPH(2.d0*SO4I+NO3I+CLI-NAI-NH4I, HI, OHI)

         ! Na2SO4 DISSOLUTION
         IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN ! PSI1
            CALL POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
            IF (ISLV.EQ.0) THEN
               PSI1 = MIN (PSI1, CHI1)
            ELSE
               PSI1 = ZERO
            ENDIF
         ELSE
            PSI1 = ZERO
         ENDIF

         ! SAVE CONCENTRATIONS IN MOLAL ARRAY 
         MOLAL (1) = HI
         MOLAL (2) = NAI
         MOLAL (3) = NH4I
         MOLAL (4) = CLI
         MOLAL (5) = SO4I
         MOLAL (6) = ZERO
         MOLAL (7) = NO3I

         !==============================================================
         ! CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY)
         !==============================================================
         GNH3      = MAX(CHI4 - PSI4, TINY)
         GHNO3     = MAX(CHI5 - PSI5, TINY)
         GHCL      = MAX(CHI6 - PSI6, TINY)
         
         CNH42S4   = ZERO
         CNH4NO3   = ZERO
         CNH4CL    = ZERO

         !==============================================================
         ! CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES
         !==============================================================
         CALL CALCMR

         !==============================================================
         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         !==============================================================
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      !=================================================================
      ! CALCULATE FUNCTION VALUE FOR OUTER LOOP 
      !=================================================================
 20   CONTINUE
      VALUE = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      ! Return to calling program
      END FUNCTION FUNCG4A

!------------------------------------------------------------------------------

      SUBROUTINE CALCG5
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCG5
! *** CASE G5
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER  :: K
      REAL*8  ::     PSI6LO, PSI6HI, X1, Y1, DX
      REAL*8  ::     X2, Y2, X3, Y3, DELTA

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
!$OMP THREADPRIVATE( /CASEG/ )

      !=================================================================
      ! CALCG5 begins here!
      !=================================================================

      CALAOU = .TRUE.   
      CHI1   = 0.5d0*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
 
      PSI1   = CHI1
      PSI2   = CHI2
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      WATER  = CHI2/M0(4) + CHI1/M0(2)

      ! Initial values for bisection
      X1 = PSI6LO
      Y1 = FUNCG5A (X1)
      IF (CHI6.LE.TINY) GOTO 50  
      ! comment out
      !IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
      !IF (WATER .LE. TINY) RETURN                    ! No water

      !=================================================================
      ! Root tracking ; for the range of HI and LO
      !=================================================================
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCG5A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE

      !=================================================================
      ! No subdivision with solution; IF ABS(Y2)<EPS solution is assumed
      !=================================================================
      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG5A (PSI6HI)
      GOTO 50

      !=================================================================
      ! Perform bisection
      !=================================================================
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCG5A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE

      !=================================================================
      ! Converged ; return 
      !=================================================================
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCG5A (X3)

      !=================================================================
      ! Calculate HSO4 speciation and return
      !=================================================================
 50   CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  ! If quadrat.called
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA              ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA              ! SO4  EFFECT
         MOLAL(6) = DELTA                         ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCG5

!------------------------------------------------------------------------------

      FUNCTION FUNCG5A( X ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCG5A
! *** CASE G5
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h"  

      ! Arguments
      REAL*8, INTENT(IN)  :: X

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     AKK, BB, CC, DD, SMIN, HI, OHI
      REAL*8  ::     VALUE

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
!$OMP THREADPRIVATE( /CASEG/ )

      !=================================================================
      ! FUNCG5A begins here!
      !=================================================================
      PSI6   = X
      FRST   = .TRUE.
      CALAIN = .TRUE. 

      !=================================================================       
      ! Solve equations ; with iterations for activity coef.
      !=================================================================
      DO 10 K=1,NSWEEP
         A1  = XK5 *(WATER/GAMA(2))**3.d0
         A2  = XK7 *(WATER/GAMA(4))**3.d0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.d0
         A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.d0
         AKK = A4*A6

         !==============================================================       
         ! Calculate dissociation quantities
         !==============================================================
         IF (CHI5.GE.TINY) THEN
            PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
         ELSE
            PSI5 = TINY
         ENDIF

         ! Comment out
         !IF(CHI4.GT.TINY) THEN
         IF(W(2).GT.TINY) THEN  
            ! Accounts for NH3 evaporation
            BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
            CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
            DD   = MAX(BB*BB-4.d0*CC,ZERO) 

            ! Patch proposed by Uma Shankar, 19/11/01
            PSI4 =0.5d0*(-BB - SQRT(DD))
         ELSE
            PSI4 = TINY
         ENDIF

         ! Speciation & water content
         MOLAL (2) = 2.0D0*PSI1                 ! NAI
         MOLAL (3) = 2.d0*PSI2 + PSI4           ! NH4I
         MOLAL (4) = PSI6                       ! CLI
         MOLAL (5) = PSI2 + PSI1                ! SO4I
         MOLAL (6) = ZERO
         MOLAL (7) = PSI5                       ! NO3I

         SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+
     &               MOLAL(4)-MOLAL(2)-MOLAL(3)
         CALL CALCPH (SMIN, HI, OHI)
         MOLAL (1) = HI
 
         GNH3      = MAX(CHI4 - PSI4, TINY)     ! Gas NH3
         GHNO3     = MAX(CHI5 - PSI5, TINY)     ! Gas HNO3
         GHCL      = MAX(CHI6 - PSI6, TINY)     ! Gas HCl

         CNH42S4   = ZERO                       ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                       ! Solid NH4NO3
         CNH4CL    = ZERO                       ! Solid NH4Cl

         CALL CALCMR                            ! Water content

         ! Calculate activities or terminate internal loop
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      !=================================================================
      ! Calculate objective function
      !=================================================================
20    CONTINUE
      VALUE = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE

      ! Return to calling program
      END FUNCTION FUNCG5A

!------------------------------------------------------------------------------

      SUBROUTINE CALCH1
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCH1
! *** CASE H1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4
!
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCH1A)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!==============================================================================
!
#     include "isoropia.h" 



      !=================================================================
      ! CALCH1 begins here!
      !=================================================================

      IF ( RHB .LT. DRMH1 ) THEN    
         SCASE = 'H1 ; SUBCASE 1'  
         CALL CALCH1A            ! SOLID PHASE ONLY POSSIBLE
         SCASE = 'H1 ; SUBCASE 1'
      ELSE
         SCASE = 'H1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
         CALL CALCMDRH (RHB, DRMH1, DRNH4NO3, 
     &                  CALCH1A, CALCH2A)
         SCASE = 'H1 ; SUBCASE 2'
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCH1 

!-----------------------------------------------------------------------------

      SUBROUTINE CALCH1A
!
!*****************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCH1A
! *** CASE H1 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NANO3, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!*****************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      REAL*8  :: LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2, NAFR, NO3FR
      REAL*8  :: CLFR, ALF, BET, GAM, RTSQ, THETA1, THETA2
      REAL*8  :: BB, CC, DD, SQDD, DD1, DD2, SQDD1, SQDD2
      REAL*8  :: A1, A2

      !=================================================================
      ! CALCH1A begins here!
      !=================================================================

      ! CALCULATE NON VOLATILE SOLIDS 
      CNA2SO4 = W(2)
      CNH42S4 = ZERO
      NAFR    = MAX (W(1)-2*CNA2SO4, ZERO)
      CNANO3  = MIN (NAFR, W(4))
      NO3FR   = MAX (W(4)-CNANO3, ZERO)
      CNACL   = MIN (MAX(NAFR-CNANO3, ZERO), W(5))
      CLFR    = MAX (W(5)-CNACL, ZERO)

      ! CALCULATE VOLATILE SPECIES 
      ALF     = W(3)                     ! FREE NH3
      BET     = CLFR                     ! FREE CL
      GAM     = NO3FR                    ! FREE NO3

      RTSQ    = R*TEMP*R*TEMP
      A1      = XK6/RTSQ
      A2      = XK10/RTSQ

      THETA1  = GAM - BET*(A2/A1)
      THETA2  = A2/A1

      ! QUADRATIC EQUATION SOLUTION
      BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
      CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
      DD      = BB*BB - 4.0D0*CC
      IF (DD.LT.ZERO) GOTO 100   ! Solve each reaction seperately

      ! TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID
      SQDD    = SQRT(DD)
      KAPA1   = 0.5D0*(-BB+SQDD)
      KAPA2   = 0.5D0*(-BB-SQDD)
      LAMDA1  = THETA1 + THETA2*KAPA1
      LAMDA2  = THETA1 + THETA2*KAPA2

      IF (KAPA1.GE.ZERO .AND. LAMDA1.GE.ZERO) THEN
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND.
     &       BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF

      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND. 
     &       BET-KAPA2.GE.ZERO .AND. GAM-LAMDA2.GE.ZERO) THEN
             KAPA = KAPA2
             LAMDA= LAMDA2
             GOTO 200
         ENDIF
      ENDIF

      ! SEPARATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA 
100   CONTINUE
      KAPA  = ZERO
      LAMDA = ZERO
      DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
      DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)

      ! NH4CL EQUILIBRIUM
      IF (DD1.GE.ZERO) THEN
         SQDD1 = SQRT(DD1)
         KAPA1 = 0.5D0*(ALF+BET + SQDD1)
         KAPA2 = 0.5D0*(ALF+BET - SQDD1)

         IF (KAPA1.GE.ZERO .AND. KAPA1.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA1 
         ELSE IF (KAPA2.GE.ZERO .AND. KAPA2.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA2
         ELSE
            KAPA = ZERO
         ENDIF
      ENDIF

      ! NH4NO3 EQUILIBRIUM
      IF (DD2.GE.ZERO) THEN
         SQDD2 = SQRT(DD2)
         LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
         LAMDA2= 0.5D0*(ALF+GAM - SQDD2)

         IF (LAMDA1.GE.ZERO .AND. LAMDA1.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1 
         ELSE IF (LAMDA2.GE.ZERO .AND. LAMDA2.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
         ELSE
            LAMDA = ZERO
         ENDIF
      ENDIF

      ! IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION
      IF (KAPA.GT.ZERO .AND. LAMDA.GT.ZERO) THEN
         IF (BET .LT. LAMDA/THETA1) THEN
            KAPA = ZERO
         ELSE
            LAMDA= ZERO
         ENDIF
      ENDIF

      ! CALCULATE COMPOSITION OF VOLATILE SPECIES
 200  CONTINUE
      CNH4NO3 = LAMDA
      CNH4CL  = KAPA

      GNH3    = ALF - KAPA - LAMDA
      GHNO3   = GAM - LAMDA
      GHCL    = BET - KAPA

      ! Return to calling program
      END SUBROUTINE CALCH1A 

!------------------------------------------------------------------------------

      SUBROUTINE CALCH2
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCH2
! *** CASE H2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : NH4Cl, NA2SO4, NANO3, NACL
!
!     THERE ARE THREE REGIMES IN THIS CASE:
!     1. NH4NO3(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCH2A)
!     2. NH4NO3(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY 
!     3. NH4NO3(s) NOT POSSIBLE, AND RH >= MDRH. (MDRH REGION)
!
!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES H1A, H2B
!     RESPECTIVELY (BECAUSE MDRH POINTS COINCIDE).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      !=================================================================
      ! CALCH2 begins here!
      !=================================================================

      ! REGIME DEPENDS ON THE EXISTANCE OF NITRATES
      IF ( W(4) .GT. TINY ) THEN        
 
         ! NO3 EXISTS, WATER POSSIBLE
         SCASE = 'H2 ; SUBCASE 1'  
         CALL CALCH2A                                 
         SCASE = 'H2 ; SUBCASE 1'  

      ELSE                          

         ! NO3 NON EXISTANT, WATER NOT POSSIBLE
         SCASE = 'H2 ; SUBCASE 1'  
         CALL CALCH1A
         SCASE = 'H2 ; SUBCASE 1'  

      ENDIF

      IF ( WATER .LE. TINY .AND. RHB .LT. DRMH2 ) THEN      

         ! DRY AEROSOL
         SCASE = 'H2 ; SUBCASE 2'  

      ELSE IF ( WATER .LE. TINY .AND. RHB .GE. DRMH2 ) THEN  

         ! MDRH OF H2
         SCASE = 'H2 ; SUBCASE 3'
         CALL CALCMDRH (RHB, DRMH2, DRNANO3, 
     &                  CALCH1A, CALCH3)
         SCASE = 'H2 ; SUBCASE 3'

      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCH2 

!------------------------------------------------------------------------------

      SUBROUTINE CALCH2A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCH2A
! *** CASE H2 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::      K
      REAL*8  ::     FRNA, PSI6LO, PSI6HI, X1, Y1, DX
      REAL*8  ::     X2, Y2, X3, Y3, DELTA

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCH2A begins here!
      !=================================================================

      ! SETUP PARAMETERS 
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)

      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      !=================================================================
      ! INITIAL VALUES FOR BISECTION 
      !=================================================================
      X1 = PSI6LO
      Y1 = FUNCH2A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO 
      !=================================================================
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH2A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
      !=================================================================
      IF (Y2 .GT. EPS) Y2 = FUNCH2A (PSI6LO)
      GOTO 50

      !=================================================================
      ! PERFORM BISECTION 
      !=================================================================
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCH2A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCH2A')    ! WARNING ERROR: NO CONVERGENCE

      !=================================================================
      ! CONVERGED ; RETURN 
      !=================================================================
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCH2A (X3)

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN 
      !=================================================================
 50   CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA           ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA           ! SO4  EFFECT
         MOLAL(6) = DELTA                      ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCH2A 

!------------------------------------------------------------------------------

      FUNCTION FUNCH2A( X ) RESULT( VALUE )
!
!******************************************************************************
! 
!  *** ISORROPIA CODE
!  *** SUBROUTINE FUNCH2A
!  *** CASE H2 ; SUBCASE 1
! 
!      THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!      1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!      2. THERE IS BOTH A LIQUID & SOLID PHASE
!      3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
! 
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
! 
#     include "isoropia.h" 
 
      ! Arguments
      REAL*8, INTENT(IN)  :: X

      ! Local variables
      INTEGER ::     K, ISLV
      REAL*8  ::     VALUE
      REAL*8  ::     A64, A9, BB, CC, DD, DIAK, AA
      REAL*8  ::     SMIN, HI, OHI, DELT, PSI31, PSI32

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCH2A begins here!
      !=================================================================

      ! SETUP PARAMETERS
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A1  = XK5 *(WATER/GAMA(2))**3.d0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.d0
         A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.d0
         A7  = XK8 *(WATER/GAMA(1))**2.d0
         A8  = XK9 *(WATER/GAMA(3))**2.d0
         A64 = (XK3*XK2/XKW)*(GAMA(10)/GAMA(5)/GAMA(11))**2.d0
         A64 = A64*(R*TEMP*WATER)**2.d0
         A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

         ! CALCULATE DISSOCIATION QUANTITIES
         PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
         PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
         PSI5 = MAX(PSI5, TINY)

         IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN 
            ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
         ELSE
            PSI4 = TINY
         ENDIF

         IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN 
            ! NACL DISSOLUTION
            DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
            PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
         ENDIF

         IF (CHI8.GT.TINY .AND. WATER.GT.TINY) THEN     
            ! NANO3 DISSOLUTION
            DIAK = (PSI7-PSI5)**2.D0 + 4.D0*A8
            PSI8 = 0.5D0*( -(PSI7+PSI5) + SQRT(DIAK) )
            PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
         ENDIF

         IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN     
            ! NA2SO4 DISSOLUTION
            AA = PSI7+PSI8
            BB = AA*AA
            CC =-A1/4.D0
            CALL POLY3 (AA, BB, CC, PSI1, ISLV)
            IF (ISLV.EQ.0) THEN
               PSI1 = MIN (PSI1, CHI1)
            ELSE
               PSI1 = ZERO
            ENDIF
         ENDIF

         ! CALCULATE SPECIATION 
         MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                 ! NAI
         MOLAL (3) = PSI4                                    ! NH4I
         MOLAL (4) = PSI6 + PSI7                             ! CLI
         MOLAL (5) = PSI2 + PSI1                             ! SO4I
         MOLAL (6) = ZERO                                    ! HSO4I
         MOLAL (7) = PSI5 + PSI8                             ! NO3I

         SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+
     &               MOLAL(4)-MOLAL(2)-MOLAL(3)
         CALL CALCPH (SMIN, HI, OHI)
         MOLAL (1) = HI

         GNH3      = MAX(CHI4 - PSI4, TINY)
         GHNO3     = MAX(CHI5 - PSI5, TINY)
         GHCL      = MAX(CHI6 - PSI6, TINY)

         CNH42S4   = ZERO
         CNH4NO3   = ZERO
         CNACL     = MAX(CHI7 - PSI7, ZERO)
         CNANO3    = MAX(CHI8 - PSI8, ZERO)
         CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 

         ! NH4Cl(s) calculations
         A3   = XK6 /(R*TEMP*R*TEMP)
         DELT = MIN(GNH3, GHCL)
         BB = -(GNH3+GHCL)
         CC = GNH3*GHCL-A3
         DD = BB*BB - 4.D0*CC
         PSI31 = 0.5D0*(-BB + SQRT(DD))
         PSI32 = 0.5D0*(-BB - SQRT(DD))
         IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
            PSI3 = PSI31
         ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
            PSI3 = PSI32
         ELSE
            PSI3 = ZERO
         ENDIF

         ! CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY)
         GNH3    = MAX(GNH3 - PSI3, TINY)
         GHCL    = MAX(GHCL - PSI3, TINY)
         CNH4CL  = PSI3

         ! Water content
         CALL CALCMR 

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE FUNCTION VALUE FOR OUTER LOOP
 20   CONTINUE
      VALUE = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A64 - ONE

      ! Return to calling program
      END FUNCTION FUNCH2A

!-----------------------------------------------------------------------

      SUBROUTINE CALCH3

!***********************************************************************
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCH3
C *** CASE H3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     FRNA, PSI6LO, PSI6HI, X1, Y1, DX
      REAL*8  ::     X2, Y2, X3, Y3, DELTA

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCH3 begins here!
      !=================================================================

      ! REGIME DEPENDS ON THE EXISTANCE OF NITRATES 
      IF ( W(4) .LE. TINY ) THEN    ! NO3 NOT EXIST, WATER NOT POSSIBLE
         SCASE = 'H3'  
         CALL CALCH1A
         SCASE = 'H3'  
         RETURN
      ENDIF
C
      ! SETUP PARAMETERS
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)

      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      !=================================================================
      ! INITIAL VALUES FOR BISECTION
      !=================================================================
      X1 = PSI6LO
      Y1 = FUNCH3A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH3A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20 ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
      !=================================================================
      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH3A (PSI6LO)
      GOTO 50

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCH3A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCH3')    ! WARNING ERROR: NO CONVERGENCE

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCH3A (X3)

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN
      !=================================================================
 50   CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA         ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA         ! SO4  EFFECT
         MOLAL(6) = DELTA                    ! HSO4 EFFECT
      ENDIF
      
      ! Return to calling program
      END SUBROUTINE CALCH3 

!------------------------------------------------------------------------------

      FUNCTION FUNCH3A( X ) RESULT( VALUE )
!
!******************************************************************************
!
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCH3A
! *** CASE H3
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN) :: X

      ! Local variables
      INTEGER ::     K, ISLV
      REAL*8  ::     VALUE
      REAL*8  ::     A9, BB, CC, DD, DIAK, AA, SMIN, HI, OHI
      REAL*8  ::     DELT, PSI31, PSI32
 
      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCH3A begins here!
      !=================================================================

      ! SETUP PARAMETERS
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 

      !=================================================================
      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      !=================================================================
      DO 10 K=1,NSWEEP

         A1  = XK5 *(WATER/GAMA(2))**3.d0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.d0
         A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.d0
         A7  = XK8 *(WATER/GAMA(1))**2.d0
         A8  = XK9 *(WATER/GAMA(3))**2.d0
         A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

         ! CALCULATE DISSOCIATION QUANTITIES
         PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
         PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
         PSI5 = MAX(PSI5, TINY)

         IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN  
            ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
         ELSE
            PSI4 = TINY
         ENDIF

         IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN 
            ! NACL DISSOLUTION
            DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
            PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
            PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
         ENDIF

         IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN 
            ! NA2SO4 DISSOLUTION
            AA = PSI7+PSI8
            BB = AA*AA
            CC =-A1/4.D0
            CALL POLY3 (AA, BB, CC, PSI1, ISLV)
            IF (ISLV.EQ.0) THEN
               PSI1 = MIN (PSI1, CHI1)
            ELSE
               PSI1 = ZERO
            ENDIF
         ENDIF

         !==============================================================
         ! CALCULATE SPECIATION
         !==============================================================
         MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1             ! NAI
         MOLAL (3) = PSI4                                ! NH4I
         MOLAL (4) = PSI6 + PSI7                         ! CLI
         MOLAL (5) = PSI2 + PSI1                         ! SO4I
         MOLAL (6) = ZERO
         MOLAL (7) = PSI5 + PSI8                         ! NO3I
         SMIN      = 2.d0*MOLAL(5)+MOLAL(7)
     &             + MOLAL(4)-MOLAL(2)-MOLAL(3)
         CALL CALCPH (SMIN, HI, OHI)
         MOLAL (1) = HI

         !==============================================================
         ! CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY)
         !==============================================================
         GNH3      = MAX(CHI4 - PSI4, TINY)
         GHNO3     = MAX(CHI5 - PSI5, TINY)
         GHCL      = MAX(CHI6 - PSI6, TINY)

         CNH42S4   = ZERO
         CNH4NO3   = ZERO
         CNACL     = MAX(CHI7 - PSI7, ZERO)
         CNANO3    = MAX(CHI8 - PSI8, ZERO)
         CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 

         ! NH4Cl(s) calculations
         A3   = XK6 /(R*TEMP*R*TEMP)
         DELT = MIN(GNH3, GHCL)
         BB = -(GNH3+GHCL)
         CC = GNH3*GHCL-A3
         DD = BB*BB - 4.D0*CC
         PSI31 = 0.5D0*(-BB + SQRT(DD))
         PSI32 = 0.5D0*(-BB - SQRT(DD))
         IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
            PSI3 = PSI31
         ELSEIF ( DELT-PSI32 .GT. ZERO .AND. PSI32 .GT. ZERO ) THEN
            PSI3 = PSI32
         ELSE
            PSI3 = ZERO
         ENDIF

         !==============================================================
         ! CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY)
         !==============================================================
         GNH3    = MAX( GNH3 - PSI3, TINY )
         GHCL    = MAX( GHCL - PSI3, TINY )
         CNH4CL  = PSI3
     
         ! Water content
         CALL CALCMR                                 

         !==============================================================
         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         !==============================================================
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      !=================================================================
      ! CALCULATE FUNCTION VALUE FOR OUTER LOOP 
      !=================================================================
20    CONTINUE
      VALUE = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      ! Return to calling program
      END FUNCTION FUNCH3A 

!------------------------------------------------------------------------------

      SUBROUTINE CALCH4
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCH4
! *** CASE H4
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     FRNA, PSI6LO, PSI6HI, X1, Y1, DX
      REAL*8  ::     X2, Y2, X3, Y3, DELTA

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )


      !=================================================================
      ! CALCH4 begins here!
      !=================================================================

      ! REGIME DEPENDS ON THE EXISTANCE OF NITRATES 
      IF ( W(4) .LE. TINY .AND. 
     &     W(5) .LE. TINY ) THEN  
         SCASE = 'H4'  
         CALL CALCH1A
         SCASE = 'H4'  
         RETURN
      ENDIF

      ! SETUP PARAMETERS
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
C
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      !=================================================================
      ! INITIAL VALUES FOR BISECTION
      !=================================================================
      X1 = PSI6LO
      Y1 = FUNCH4A (X1)
      IF ( ABS(Y1) .LE. EPS .OR. CHI6 .LE. TINY ) GOTO 50  

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH4A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE
      
      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
      !=================================================================
      IF (ABS(Y2) .GT. EPS) Y2 = FUNCH4A (PSI6LO)
      GOTO 50

      !=================================================================
      ! PERFORM BISECTION 
      !=================================================================
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCH4A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCH4') ! WARNING ERROR: NO CONVERGENCE

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCH4A (X3)

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN 
      !=================================================================
 50   CONTINUE
      IF ( MOLAL(1) .GT. TINY .AND. 
     &     MOLAL(5) .GT. TINY ) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA          ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA          ! SO4  EFFECT
         MOLAL(6) = DELTA                     ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCH4 

!------------------------------------------------------------------------------

      FUNCTION FUNCH4A( X ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCH4A
! *** CASE H4
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: X

      ! Local variables
      INTEGER ::     K, ISLV
      REAL*8  ::     VALUE
      REAL*8  ::     A9, BB, CC, DD, AA, SMIN, HI, OHI, DELT
      REAL*8  ::     PSI31, PSI32

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )


      !=================================================================
      ! FUNCH4A begins here!
      !=================================================================

      ! SETUP PARAMETERS
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A1  = XK5 *(WATER/GAMA(2))**3.d0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.d0
         A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.d0
         A7  = XK8 *(WATER/GAMA(1))**2.d0
         A8  = XK9 *(WATER/GAMA(3))**2.d0
         A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
      
         ! CALCULATE DISSOCIATION QUANTITIES       
         PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
         PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
         PSI5 = MAX(PSI5, TINY)

         IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN  
            ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
         ELSE
            PSI4 = TINY
         ENDIF

         IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN 
            ! NA2SO4 DISSOLUTION
            AA = PSI7+PSI8
            BB = AA*AA
            CC =-A1/4.D0
            CALL POLY3 (AA, BB, CC, PSI1, ISLV)
            IF (ISLV.EQ.0) THEN
               PSI1 = MIN (PSI1, CHI1)
            ELSE
               PSI1 = ZERO
            ENDIF
         ENDIF

         ! CALCULATE SPECIATION 
         MOLAL(2) = PSI8 + PSI7 + 2.D0*PSI1                ! NAI
         MOLAL(3) = PSI4                                   ! NH4I
         MOLAL(4) = PSI6 + PSI7                            ! CLI
         MOLAL(5) = PSI2 + PSI1                            ! SO4I
         MOLAL(6) = ZERO
         MOLAL(7) = PSI5 + PSI8                            ! NO3I

         SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+
     &               MOLAL(4)-MOLAL(2)-MOLAL(3)
         CALL CALCPH (SMIN, HI, OHI)
         MOLAL(1) = HI

         ! CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY)
         GNH3      = MAX(CHI4 - PSI4, TINY)
         GHNO3     = MAX(CHI5 - PSI5, TINY)
         GHCL      = MAX(CHI6 - PSI6, TINY)

         CNH42S4   = ZERO
         CNH4NO3   = ZERO
         CNACL     = MAX(CHI7 - PSI7, ZERO)
         CNANO3    = MAX(CHI8 - PSI8, ZERO)
         CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 

         ! NH4Cl(s) calculations
         A3   = XK6 /(R*TEMP*R*TEMP)
         DELT = MIN(GNH3, GHCL)
         BB = -(GNH3+GHCL)
         CC = GNH3*GHCL-A3
         DD = BB*BB - 4.D0*CC
         PSI31 = 0.5D0*(-BB + SQRT(DD))
         PSI32 = 0.5D0*(-BB - SQRT(DD))
         IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
            PSI3 = PSI31
         ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
            PSI3 = PSI32
         ELSE
            PSI3 = ZERO
         ENDIF

         ! CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY)
         GNH3    = MAX(GNH3 - PSI3, TINY)
         GHCL    = MAX(GHCL - PSI3, TINY)
         CNH4CL  = PSI3

         ! Water content
         CALL CALCMR                           
         
         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE FUNCTION VALUE FOR OUTER LOOP
 20   CONTINUE
      VALUE = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      ! Return to calling program
      END FUNCTION FUNCH4A 

!------------------------------------------------------------------------------

      SUBROUTINE CALCH5
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCH5
! *** CASE H5
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::      K
      REAL*8  ::     FRNA, PSI6LO, PSI6HI, X1, Y1, DX
      REAL*8  ::     X2, Y2, X3, Y3, DELTA

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCH5 begins here!
      !=================================================================

      ! REGIME DEPENDS ON THE EXISTANCE OF NITRATES 
      IF ( W(4) .LE. TINY .AND. 
     &     W(5) .LE. TINY ) THEN  
         SCASE = 'H5'  
         CALL CALCH1A
         SCASE = 'H5'  
         RETURN
      ENDIF

      ! SETUP PARAMETERS
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)

      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      !=================================================================
      ! INITIAL VALUES FOR BISECTION
      !=================================================================
      X1 = PSI6LO
      Y1 = FUNCH5A (X1)
      IF ( ABS(Y1) .LE. EPS .OR. CHI6 .LE. TINY ) GOTO 50  

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH5A (X2)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y2) .LT. ZERO ) GOTO 20 ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
      !=================================================================
      IF ( ABS(Y2) .GT. EPS) Y2 = FUNCH5A (PSI6LO)
      GOTO 50

      !=================================================================
      ! PERFORM BISECTION 
      !=================================================================
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCH5A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS(X2-X1) .LE. EPS*X1 ) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCH5') ! WARNING ERROR: NO CONVERGENCE

      !=================================================================
      ! CONVERGED ; RETURN 
      !=================================================================
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCH5A (X3)

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN
      !=================================================================
 50   CONTINUE
      IF ( MOLAL(1) .GT. TINY .AND. 
     &     MOLAL(5) .GT. TINY ) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA               ! H+   EFECT
         MOLAL(5) = MOLAL(5) - DELTA               ! SO4  EFFECT
         MOLAL(6) = DELTA                          ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCH5 

!------------------------------------------------------------------------------

      FUNCTION FUNCH5A( X ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCH5A
! *** CASE H5
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NONE
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN) :: X

      ! Local variables
      REAL*8  ::     VALUE
      INTEGER ::     K, ISLV
      REAL*8  ::     A9, BB, CC, DD, AA, SMIN, HI, OHI

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCH5A begins here!
      !=================================================================

      ! SETUP PARAMETERS 
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A1  = XK5 *(WATER/GAMA(2))**3.d0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.d0
         A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.d0
         A7  = XK8 *(WATER/GAMA(1))**2.d0
         A8  = XK9 *(WATER/GAMA(3))**2.d0
         A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

         ! CALCULATE DISSOCIATION QUANTITIES
         PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
         PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
         PSI5 = MAX(PSI5, TINY)

         ! First try 3rd order soln
         IF ( CHI1 .GT. TINY .AND. WATER .GT. TINY ) THEN  
            BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
         ELSE
            PSI4 = TINY
         ENDIF

         IF ( CHI1 .GT. TINY .AND. WATER .GT. TINY ) THEN     
            ! NA2SO4 DISSOLUTION
            AA = PSI7+PSI8
            BB = AA*AA
            CC =-A1/4.D0
            CALL POLY3 (AA, BB, CC, PSI1, ISLV)
            IF (ISLV.EQ.0) THEN
               PSI1 = MIN (PSI1, CHI1)
            ELSE
               PSI1 = ZERO
            ENDIF
         ENDIF

         ! CALCULATE SPECIATION 
         MOLAL(2) = PSI8 + PSI7 + 2.D0*PSI1                ! NAI
         MOLAL(3) = PSI4                                   ! NH4I
         MOLAL(4) = PSI6 + PSI7                            ! CLI
         MOLAL(5) = PSI2 + PSI1                            ! SO4I
         MOLAL(6) = ZERO
         MOLAL(7) = PSI5 + PSI8                            ! NO3I

         SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+
     &               MOLAL(4)-MOLAL(2)-MOLAL(3)
         CALL CALCPH (SMIN, HI, OHI)
         MOLAL(1) = HI

         GNH3      = MAX(CHI4 - PSI4, TINY)
         GHNO3     = MAX(CHI5 - PSI5, TINY)
         GHCL      = MAX(CHI6 - PSI6, TINY)

         CNH42S4   = ZERO
         CNH4NO3   = ZERO
         CNACL     = MAX(CHI7 - PSI7, ZERO)
         CNANO3    = MAX(CHI8 - PSI8, ZERO)
         CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 

         ! Water content
         CALL CALCMR                               

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE FUNCTION VALUE FOR OUTER LOOP
 20   CONTINUE
      VALUE = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      ! Return to calling program
      END FUNCTION FUNCH5A 

!------------------------------------------------------------------------------

      SUBROUTINE CALCH6
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCH6
! *** CASE H6
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     FRNA, PSI6LO, PSI6HI, X1, Y1, DX
      REAL*8  ::     X2, Y2, X3, Y3, DELTA

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCH6 begins here!
      !=================================================================

      ! SETUP PARAMETERS 
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)

      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      !=================================================================
      ! INITIAL VALUES FOR BISECTION
      !=================================================================
      X1 = PSI6LO
      Y1 = FUNCH6A (X1)
      IF ( ABS(Y1) .LE. EPS .OR. CHI6 .LE. TINY ) GOTO 50  

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO 
      !=================================================================
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH6A (X2)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y2) .LT. ZERO ) GOTO 20 ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
      !=================================================================
      IF ( ABS(Y2) .GT. EPS ) Y2 = FUNCH6A (PSI6LO)
      GOTO 50

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCH6A (X3)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO ) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS(X2-X1) .LE. EPS*X1 ) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCH6') ! WARNING ERROR: NO CONVERGENCE

      !=================================================================
      ! CONVERGED ; RETURN 
      !=================================================================
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCH6A (X3)

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN
      !=================================================================
 50   CONTINUE
      IF ( MOLAL(1) .GT. TINY .AND. 
     &     MOLAL(5) .GT. TINY ) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5),  ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA              ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA              ! SO4  EFFECT
         MOLAL(6) = DELTA                         ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCH6 

!------------------------------------------------------------------------------

      FUNCTION FUNCH6A( X ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCH6A
! *** CASE H6
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: X

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     VALUE
      REAL*8  ::     A9, BB, CC, DD, SMIN, HI, OHI

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCH6A begins here!
      !=================================================================

      ! SETUP PARAMETERS 
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A1  = XK5 *(WATER/GAMA(2))**3.d0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.d0
         A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.d0
         A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.d0
         A7  = XK8 *(WATER/GAMA(1))**2.d0
         A8  = XK9 *(WATER/GAMA(3))**2.d0
         A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

         ! CALCULATE DISSOCIATION QUANTITIES
         PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
         PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
         PSI5 = MAX(PSI5, TINY)

	 ! First try 3rd order soln
         IF ( CHI1 .GT. TINY .AND. WATER .GT. TINY ) THEN  
            BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
         ELSE
            PSI4 = TINY
         ENDIF

         ! CALCULATE SPECIATION
         MOLAL(2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
         MOLAL(3) = PSI4                                  ! NH4I
         MOLAL(4) = PSI6 + PSI7                           ! CLI
         MOLAL(5) = PSI2 + PSI1                           ! SO4I
         MOLAL(6) = ZERO                                  ! HSO4I
         MOLAL(7) = PSI5 + PSI8                           ! NO3I

         SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+
     &               MOLAL(4)-MOLAL(2)-MOLAL(3)
         CALL CALCPH (SMIN, HI, OHI)
         MOLAL(1) = HI

         GNH3      = MAX(CHI4 - PSI4, TINY)
         GHNO3     = MAX(CHI5 - PSI5, TINY)
         GHCL      = MAX(CHI6 - PSI6, TINY)
         
         CNH42S4   = ZERO
         CNH4NO3   = ZERO
         CNACL     = MAX(CHI7 - PSI7, ZERO)
         CNANO3    = MAX(CHI8 - PSI8, ZERO)
         CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 
         
         ! Water content
         CALL CALCMR            

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP 
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE FUNCTION VALUE FOR OUTER LOOP
 20   CONTINUE
      VALUE = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE

      ! Return to calling program
      END FUNCTION FUNCH6A 

!------------------------------------------------------------------------------

      SUBROUTINE CALCI1
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCI1
! *** CASE I1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4
!
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCI1A)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!**************************************************************************
!
#     include "isoropia.h" 

      !=================================================================
      ! CALCI1 begins here!
      !=================================================================

      ! REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY 
      IF (RHB.LT.DRMI1) THEN    

         ! SOLID PHASE ONLY POSSIBLE
         SCASE = 'I1 ; SUBCASE 1'  
         CALL CALCI1A              ! SOLID PHASE ONLY POSSIBLE
         SCASE = 'I1 ; SUBCASE 1'

      ELSE

         ! LIQUID & SOLID PHASE POSSIBLE
         SCASE = 'I1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
         CALL CALCMDRH (RHB, DRMI1, DRNH4HS4, 
     &                  CALCI1A, CALCI2A)
         SCASE = 'I1 ; SUBCASE 2'

      ENDIF

      ! AMMONIA IN GAS PHASE
      CALL CALCNH3

      ! Return to calling program
      END SUBROUTINE CALCI1 

!------------------------------------------------------------------------------

      SUBROUTINE CALCI1A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCI1A
! *** CASE I1 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      REAL*8  :: FRSO4, FRNH4

      !=================================================================
      ! CALCI1A begins here!
      !=================================================================

      ! CALCULATE NON VOLATILE SOLIDS
      CNA2SO4 = 0.5D0*W(1)
      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CNH42S4 = ZERO
      FRSO4   = MAX(W(2)-CNA2SO4, ZERO)

      CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
      FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
      FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)

      IF ( FRSO4 .LE. TINY ) THEN
         CLC     = MAX(CLC - FRNH4, ZERO)
         CNH42S4 = 2.D0*FRNH4

      ELSEIF ( FRNH4 .LE. TINY ) THEN
         CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
         CLC     = MAX(CLC-FRSO4, ZERO)
         IF (CNA2SO4.GT.TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
         ENDIF
      ENDIF

      ! CALCULATE GAS SPECIES
      GHNO3 = W(4)
      GHCL  = W(5)
      GNH3  = ZERO

       ! Return to calling program
      END SUBROUTINE CALCI1A 

!------------------------------------------------------------------------------

      SUBROUTINE CALCI2
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCI2
! *** CASE I2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
!
!     THERE ARE THREE REGIMES IN THIS CASE:
!     1. NH4HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCI2A)
!     2. NH4HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY 
!     3. NH4HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL 
!
!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES I1A, I2B
!     RESPECTIVELY
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local varaibles
      INTEGER :: K

      !=================================================================
      ! CALCI2 begins here!
      !=================================================================

      ! FIND DRY COMPOSITION 
      CALL CALCI1A

      ! REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH
      IF ( CNH4HS4 .GT. TINY ) THEN
         SCASE = 'I2 ; SUBCASE 1'  
         CALL CALCI2A                       
         SCASE = 'I2 ; SUBCASE 1'  
      ENDIF

      IF ( WATER .LE. TINY ) THEN
         IF ( RHB .LT. DRMI2 ) THEN         

            ! SOLID SOLUTION ONLY
            WATER = TINY
            DO 10 K=1,NIONS
               MOLAL(K) = ZERO
 10         CONTINUE
            CALL CALCI1A
            SCASE = 'I2 ; SUBCASE 2'  

         ELSE IF ( RHB .GE. DRMI2 ) THEN 

            ! MDRH OF I2
            SCASE = 'I2 ; SUBCASE 3'
            CALL CALCMDRH (RHB, DRMI2, DRNAHSO4, CALCI1A, CALCI3A)
            SCASE = 'I2 ; SUBCASE 3'
         ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCI2 

!------------------------------------------------------------------------------

      SUBROUTINE CALCI2A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCI2A
! *** CASE I2 ; SUBCASE A
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     PSI2LO, PSI2HI, X1, Y1, YHI, DX
      REAL*8  ::     X2, Y2, X3, Y3

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCI2A begins here!
      !=================================================================

      ! FIND DRY COMPOSITION
      ! Needed when called from CALCMDRH 
      CALL CALCI1A    

      ! SETUP PARAMETERS
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = ZERO   
      PSI3 = ZERO   
      PSI4 = ZERO  
      PSI5 = ZERO

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI2LO = ZERO                ! Low  limit
      PSI2HI = CHI2                ! High limit

      ! INITIAL VALUES FOR BISECTION
      X1 = PSI2HI
      Y1 = FUNCI2A (X1)
      YHI= Y1                      ! Save Y-value at HI position

      ! YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
      IF ( YHI .LT. EPS ) GOTO 50

      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = FUNCI2A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20 ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC  
      IF (Y2.GT.EPS) Y2 = FUNCI3A (ZERO)
      GOTO 50

      ! PERFORM BISECTION 
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCI2A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCI2A') ! WARNING ERROR: NO CONVERGENCE

      ! CONVERGED ; RETURN 
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCI2A (X3)

 50   CONTINUE
      RETURN
      
      ! Return to calling program
      END SUBROUTINE CALCI2A 

!------------------------------------------------------------------------------

      FUNCTION FUNCI2A( P2 )  RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCI2A
! *** CASE I2 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN) :: P2

      ! Local variables
      INTEGER ::     K, ISLV
      REAL*8  ::     VALUE
      REAL*8  ::     AA, BB, CC, DD

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCI2A begins here!
      !=================================================================

      ! SETUP PARAMETERS
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI2   = P2                  ! Save PSI2 in COMMON BLOCK
      PSI3   = CHI3
      PSI4   = CHI4
      PSI5   = CHI5
      PSI6   = ZERO

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.C
      DO 10 K=1,NSWEEP

         A3 = XK11*(WATER/GAMA(9))**2.d0
         A4 = XK5 *(WATER/GAMA(2))**3.d0
         A5 = XK7 *(WATER/GAMA(4))**3.d0
         A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
         A7 = SQRT(A4/A5)

         ! CALCULATE DISSOCIATION QUANTITIES
         IF ( CHI5 .GT. TINY .AND. WATER .GT. TINY ) THEN     
            PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7 
            PSI5 = MAX(MIN (PSI5, CHI5), TINY)
         ENDIF

         IF ( CHI4 .GT. TINY .AND. WATER .GT. TINY ) THEN     
            AA   = PSI2+PSI5+PSI6+PSI3
            BB   = PSI3*AA
            CC   = 0.25D0*(PSI3*PSI3*(PSI2+PSI5+PSI6)-A4)
            CALL POLY3 (AA, BB, CC, PSI4, ISLV)
            IF (ISLV.EQ.0) THEN
               PSI4 = MIN (PSI4, CHI4)
            ELSE
               PSI4 = ZERO
            ENDIF
         ENDIF

         IF (CHI3.GT.TINY .AND. WATER.GT.TINY) THEN     
            AA   = 2.D0*PSI4 + PSI2 + PSI1 - PSI6
            BB   = 2.D0*PSI4*(PSI2 + PSI1 - PSI6) - A3
            CC   = ZERO
            CALL POLY3 (AA, BB, CC, PSI3, ISLV)
            IF (ISLV.EQ.0) THEN
               PSI3 = MIN (PSI3, CHI3)
            ELSE
               PSI3 = ZERO
            ENDIF
         ENDIF
         
         BB   = PSI2 + PSI4 + PSI5 + A6 ! PSI6
         CC   =-A6*(PSI2 + PSI3 + PSI1)
         DD   = BB*BB - 4.D0*CC
         PSI6 = 0.5D0*(-BB + SQRT(DD))

         ! CALCULATE SPECIATION
         MOLAL(1) = PSI6                           ! HI
         MOLAL(2) = 2.D0*PSI4 + PSI3               ! NAI
         MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1   ! NH4I
         MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6      ! SO4I
         MOLAL(6) = PSI2 + PSI3 + PSI1 - PSI6      ! HSO4I
         CLC       = CHI2 - PSI2
         CNAHSO4   = CHI3 - PSI3
         CNA2SO4   = CHI4 - PSI4
         CNH42S4   = CHI5 - PSI5
         CNH4HS4   = ZERO
         CALL CALCMR                               ! Water content

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE FUNCTION VALUE FOR OUTER LOOP
 20   CONTINUE
      A2    = XK13*(WATER/GAMA(13))**5.d0
      VALUE = MOLAL(5)*MOLAL(6)*
     &        MOLAL(3)**3.D0/A2 - ONE

      ! Return to calling program
      END FUNCTION FUNCI2A 

!------------------------------------------------------------------------------

      SUBROUTINE CALCI3
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCI3
! *** CASE I3
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
!
!     THERE ARE THREE REGIMES IN THIS CASE:
!     1.(NA,NH4)HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCI3A)
!     2.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY 
!     3.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL 
!
!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES I1A, I2B
!     RESPECTIVELY
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER :: K

      !=================================================================
      ! CALCI3 begins here!
      !=================================================================

      ! FIND DRY COMPOSITION
      CALL CALCI1A

      ! REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RHC
      IF ( CNH4HS4 .GT. TINY .OR. CNAHSO4 .GT. TINY ) THEN
         ! FULL SOLUTION
         SCASE = 'I3 ; SUBCASE 1'  
         CALL CALCI3A                     
         SCASE = 'I3 ; SUBCASE 1'  
      ENDIF

      IF ( WATER .LE. TINY ) THEN

         IF ( RHB .LT. DRMI3 ) THEN         

            ! SOLID SOLUTION
            WATER = TINY
            DO 10 K=1,NIONS
               MOLAL(K) = ZERO
 10         CONTINUE
            CALL CALCI1A
            SCASE = 'I3 ; SUBCASE 2'  
            
         ELSE IF ( RHB .GE. DRMI3 ) THEN     

            ! MDRH OF I3
            SCASE = 'I3 ; SUBCASE 3'
            CALL CALCMDRH (RHB, DRMI3, DRLC, CALCI1A, CALCI4)
            SCASE = 'I3 ; SUBCASE 3'
         ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCI3 

!--------------------------------------------------------------------------

      SUBROUTINE CALCI3A

!***************************************************************************
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCI3A
C *** CASE I3 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC
C
C *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
#     include "isoropia.h" 


      INTEGER  :: K
      REAL*8  :: PSI2LO, PSI2HI, X1, Y1, YHI, DX
      REAL*8  :: X2, Y2, X3, Y3
      REAL*8  :: CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCI3A begins here!
      !=================================================================

      ! FIND DRY COMPOSITION 
      ! Needed when called from CALCMDRH
      CALL CALCI1A         

      ! SETUP PARAMETERS
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = ZERO   
      PSI3 = CNAHSO4
      PSI4 = ZERO  
      PSI5 = ZERO

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI2LO = ZERO                ! Low  limit
      PSI2HI = CHI2                ! High limit

      ! INITIAL VALUES FOR BISECTION
      X1 = PSI2HI
      Y1 = FUNCI3A (X1)
      YHI= Y1                      ! Save Y-value at HI position

      ! YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
      IF (YHI.LT.EPS) GOTO 50
      
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = FUNCI3A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC  
      IF (Y2.GT.EPS) Y2 = FUNCI3A (ZERO)
      GOTO 50

      ! PERFORM BISECTION
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCI3A (X3)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO ) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS(X2-X1) .LE. EPS*X1 ) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCI3A') ! WARNING ERROR: NO CONVERGENCE

      ! CONVERGED ; RETURN 
40    CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCI3A (X3)

50    CONTINUE
      RETURN

      ! Return to calling program
      END SUBROUTINE CALCI3A 

!------------------------------------------------------------------------------

      FUNCTION FUNCI3A( P2 ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCI3A
! *** CASE I3 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: P2 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     VALUE
      REAL*8  ::     PSI4LO, PSI4HI, X1, Y1, YHI, DX
      REAL*8  ::     BB, CC, DD, X2, Y2, X3, Y3

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCI3A begins here!
      !=================================================================

      ! SETUP PARAMETERS
      PSI2   = P2                  ! Save PSI2 in COMMON BLOCK
      PSI4LO = ZERO                ! Low  limit for PSI4
      PSI4HI = CHI4                ! High limit for PSI4

      ! IF NH3 =0, CALL FUNCI3B FOR Y4=0 
      IF ( CHI4 .LE. TINY ) THEN
         VALUE = FUNCI3B (ZERO)
         GOTO 50
      ENDIF

      ! INITIAL VALUES FOR BISECTION
      X1 = PSI4HI
      Y1 = FUNCI3B (X1)
      IF ( ABS(Y1) .LE. EPS ) GOTO 50
      YHI= Y1                         ! Save Y-value at HI position

      ! YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4
      IF ( YHI .LT. ZERO ) GOTO 50

      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = MAX(X1-DX, PSI4LO)
         Y2 = FUNCI3B (X2)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y2) .LT. ZERO ) GOTO 20 ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NA2SO4
      IF ( Y2 .GT. EPS ) Y2 = FUNCI3B (PSI4LO)
      GOTO 50

      ! PERFORM BISECTION
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCI3B (X3)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO ) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS(X2-X1) .LE. EPS*X1 ) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0004, 'FUNCI3A') ! WARNING ERROR: NO CONVERGENCE

      ! INNER LOOP CONVERGED
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCI3B (X3)

      ! CALCULATE FUNCTION VALUE FOR OUTER LOOP
 50   CONTINUE
      A2    = XK13*(WATER/GAMA(13))**5.d0
      VALUE = MOLAL(5)*MOLAL(6)*
     &        MOLAL(3)**3.D0/A2 - ONE

      ! Return to calling program
      END FUNCTION FUNCI3A 

!------------------------------------------------------------------------------

      FUNCTION FUNCI3B( P4 ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCI3B
! *** CASE I3 ; SUBCASE 2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC
!
!     SOLUTION IS SAVED IN COMMON BLOCK /CASE/
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: P4

      ! Function value
      INTEGER ::     K
      REAL*8  ::     VALUE
      REAL*8  ::     BB, CC, DD
      
      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCI3B begins here!
      !=================================================================

      ! SETUP PARAMETERS
      PSI4   = P4   
      FRST   = .TRUE.
      CALAIN = .TRUE.

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. 
      DO 10 K=1,NSWEEP

         A4 = XK5*(WATER/GAMA(2))**3.d0
         A5 = XK7*(WATER/GAMA(4))**3.d0
         A6 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
         A7 = SQRT(A4/A5)

         !CALCULATE DISSOCIATION QUANTITIES
         BB   = PSI2 + PSI4 + PSI5 + A6 ! PSI6
         CC   =-A6*(PSI2 + PSI3 + PSI1)
         DD   = BB*BB - 4.D0*CC
         PSI6 = 0.5D0*(-BB + SQRT(DD))

         PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7 
         PSI5 = MIN (PSI5, CHI5)

         ! CALCULATE SPECIATION
         MOLAL(1) = PSI6                                  ! HI
         MOLAL(2) = 2.D0*PSI4 + PSI3                      ! NAI
         MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1          ! NH4I
         MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6             ! SO4I
         MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 - PSI6, TINY)  ! HSO4I
         CLC      = MAX(CHI2 - PSI2, ZERO)
         CNAHSO4  = ZERO
         CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
         CNH42S4  = MAX(CHI5 - PSI5, ZERO)
         CNH4HS4  = ZERO
         CALL CALCMR                                      ! Water content

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF ( FRST. AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE OBJECTIVE FUNCTION
 20   CONTINUE
      A4    = XK5 *(WATER/GAMA(2))**3.d0    
      VALUE = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE

      ! Return to calling program
      END FUNCTION FUNCI3B 

!------------------------------------------------------------------------------

      SUBROUTINE CALCI4
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCI4
! *** CASE I4
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     PSI4LO, PSI4HI, Y1, X1, YHI, DX
      REAL*8  ::     X2, Y2, YLO, Y3, X3

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCI4 begins here!
      !=================================================================

      ! FIND DRY COMPOSITION
      CALL CALCI1A

      ! SETUP PARAMETERS
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = CLC   
      PSI3 = CNAHSO4
      PSI4 = ZERO  
      PSI5 = ZERO

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = ZERO                ! Low  limit
      PSI4HI = CHI4                ! High limit

      ! IF NA2SO4(S) =0, CALL FUNCI4B FOR Y4=0
      IF ( CHI4 .LE. TINY ) THEN
         Y1 = FUNCI4A (ZERO)
         GOTO 50
      ENDIF

      ! INITIAL VALUES FOR BISECTION
      X1 = PSI4HI
      Y1 = FUNCI4A (X1)
      YHI= Y1                      ! Save Y-value at HI position

      ! YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4
      IF ( ABS(Y1) .LE. EPS .OR. YHI .LT. ZERO ) GOTO 50

      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LOC
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1-DX
         Y2 = FUNCI4A (X2)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y2) .LT. ZERO ) GOTO 20 ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH4CL  
      YLO= Y1                   ! Save Y-value at Hi position
      IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         Y3 = FUNCI4A (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN ! X2 IS A SOLUTION 
         GOTO 50
      ELSE
         CALL PUSHERR (0001, 'CALCI4') ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF

      ! PERFORM BISECTION
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCI4A (X3)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO ) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS(X2-X1) .LE. EPS*X1 ) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCI4') ! WARNING ERROR: NO CONVERGENCE

      ! CONVERGED ; RETURN
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCI4A (X3)

 50   CONTINUE
      RETURN
      
      ! Return to calling program
      END SUBROUTINE CALCI4 
         
!------------------------------------------------------------------------------

      FUNCTION FUNCI4A( P4 ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCI4A
! *** CASE I4
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: P4

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     VALUE
      REAL*8  ::     BB, CC, DD

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCI4A begins here!
      !=================================================================

      ! SETUP PARAMETERS
      PSI4   = P4     ! PSI3 already assigned in FUNCI4A
      FRST   = .TRUE.
      CALAIN = .TRUE.

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A4 = XK5 *(WATER/GAMA(2))**3.d0
         A5 = XK7 *(WATER/GAMA(4))**3.d0
         A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
         A7 = SQRT(A4/A5)

         ! CALCULATE DISSOCIATION QUANTITIES
         BB   = PSI2 + PSI4 + PSI5 + A6 ! PSI6
         CC   =-A6*(PSI2 + PSI3 + PSI1)
         DD   = BB*BB - 4.D0*CC
         PSI6 = 0.5D0*(-BB + SQRT(DD))
     
         PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7 
         PSI5 = MIN (PSI5, CHI5)

         ! CALCULATE SPECIATION
         MOLAL(1) = PSI6                            ! HI
         MOLAL(2) = 2.D0*PSI4 + PSI3                ! NAI
         MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1    ! NH4I
         MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6       ! SO4I
         MOLAL(6) = PSI2 + PSI3 + PSI1 - PSI6       ! HSO4I
         CLC       = ZERO
         CNAHSO4   = ZERO
         CNA2SO4   = CHI4 - PSI4
         CNH42S4   = CHI5 - PSI5
         CNH4HS4   = ZERO
         CALL CALCMR                                ! Water content

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE OBJECTIVE FUNCTION
 20   CONTINUE
      A4    = XK5 *(WATER/GAMA(2))**3.d0    
      VALUE = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE

      ! Return to calling program
      END FUNCTION FUNCI4A 

!------------------------------------------------------------------------------

      SUBROUTINE CALCI5
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCI5
! *** CASE I5
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     PSI4LO, PSI4HI, Y1, X1, YHI, DX
      REAL*8  ::     X2, Y2, YLO, Y3, X3

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCI5 begins here!
      !=================================================================

      ! FIND DRY COMPOSITION
      CALL CALCI1A
      
      ! SETUP PARAMETERS
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = CLC   
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = CNH42S4

      CALAOU =.TRUE.               ! Outer loop activity calculation flag
      PSI4LO = ZERO                ! Low  limit
      PSI4HI = CHI4                ! High limit

      ! IF NA2SO4(S) =0, CALL FUNCI5B FOR Y4=0
      IF ( CHI4 .LE. TINY ) THEN
         Y1 = FUNCI5A (ZERO)
         GOTO 50
      ENDIF

      ! INITIAL VALUES FOR BISECTION
      X1 = PSI4HI
      Y1 = FUNCI5A (X1)
      YHI= Y1                      ! Save Y-value at HI position

      ! YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4
      IF ( ABS(Y1) .LE. EPS .OR. YHI .LT. ZERO ) GOTO 50

      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1-DX
         Y2 = FUNCI5A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH4CL  
      YLO= Y1                      ! Save Y-value at Hi position
      IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         Y3 = FUNCI5A (ZERO)
         GOTO 50
      ELSE IF ( ABS(Y2) .LT. EPS ) THEN ! X2 IS A SOLUTION 
         GOTO 50
      ELSE
         CALL PUSHERR (0001, 'CALCI5') ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF

      ! PERFORM BISECTION
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCI5A (X3)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO ) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS(X2-X1) .LE. EPS*X1 ) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCI5') ! WARNING ERROR: NO CONVERGENCE

      ! CONVERGED ; RETURN
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCI5A (X3)
      
 50   CONTINUE
      RETURN

      ! Return to calling program
      END SUBROUTINE CALCI5 

!------------------------------------------------------------------------------

      FUNCTION FUNCI5A( P4 ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCI5A
! *** CASE I5
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: P4

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     VALUE
      REAL*8  ::     BB, CC, DD

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! FUNCI5A begins here!
      !=================================================================

      ! SETUP PARAMETERS
      PSI4   = P4     ! PSI3 already assigned in FUNCI5A
      FRST   = .TRUE.
      CALAIN = .TRUE.

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A4 = XK5 *(WATER/GAMA(2))**3.d0
         A5 = XK7 *(WATER/GAMA(4))**3.d0
         A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

         ! CALCULATE DISSOCIATION QUANTITIES
         BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
         CC   =-A6*(PSI2 + PSI3 + PSI1)
         DD   = BB*BB - 4.D0*CC
         PSI6 = 0.5D0*(-BB + SQRT(DD))

         ! CALCULATE SPECIATION
         MOLAL(1) = PSI6                            ! HI
         MOLAL(2) = 2.D0*PSI4 + PSI3                ! NAI
         MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1    ! NH4I
         MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6       ! SO4I
         MOLAL(6) = PSI2 + PSI3 + PSI1 - PSI6       ! HSO4I
         CLC       = ZERO
         CNAHSO4   = ZERO
         CNA2SO4   = CHI4 - PSI4
         CNH42S4   = ZERO
         CNH4HS4   = ZERO
         CALL CALCMR                                ! Water content

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE OBJECTIVE FUNCTION
 20   CONTINUE
      A4     = XK5 *(WATER/GAMA(2))**3.d0    
      VALUE = MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE

      ! Return to calling program
      END FUNCTION FUNCI5A 

!------------------------------------------------------------------------------

      SUBROUTINE CALCI6
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCI6
! *** CASE I6
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. SOLID & LIQUID AEROSOL POSSIBLE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     BB, CC, DD

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
!$OMP THREADPRIVATE( /SOLUT/ )

      !=================================================================
      ! CALCI6 begins here!
      !=================================================================
      
      ! FIND DRY COMPOSITION
      CALL CALCI1A

      ! SETUP PARAMETERS
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = CLC   
      PSI3 = CNAHSO4
      PSI4 = CNA2SO4
      PSI5 = CNH42S4

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

         ! CALCULATE DISSOCIATION QUANTITIES
         BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
         CC   =-A6*(PSI2 + PSI3 + PSI1)
         DD   = BB*BB - 4.D0*CC
         PSI6 = 0.5D0*(-BB + SQRT(DD))

         ! CALCULATE SPECIATION
         MOLAL(1) = PSI6                                    ! HI
         MOLAL(2) = 2.D0*PSI4 + PSI3                        ! NAI
         MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1            ! NH4I
         MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6               ! SO4I
         MOLAL(6) = PSI2 + PSI3 + PSI1 - PSI6               ! HSO4I
         CLC       = ZERO
         CNAHSO4   = ZERO
         CNA2SO4   = CHI4 - PSI4
         CNH42S4   = ZERO
         CNH4HS4   = ZERO
         CALL CALCMR                                        ! Water content

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! Exit
20    CONTINUE
      RETURN

      ! Return to calling program
      END SUBROUTINE CALCI6 

!------------------------------------------------------------------------------

      SUBROUTINE CALCJ1
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCJ1
! *** CASE J1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     PSI1LO, PSI1HI, X1,  Y1, YHI, DX
      REAL*8  ::     X2,     Y2,     YLO, Y3, X3

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3
      COMMON /CASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3
!$OMP THREADPRIVATE( /CASEJ/ )

      !=================================================================
      ! CALCJ1 begins here!
      !=================================================================

      ! SETUP PARAMETERS
      CALAOU =.TRUE.               ! Outer loop activity calculation flag
      CHI1   = W(1)                ! Total NA initially as NaHSO4
      CHI2   = W(3)                ! Total NH4 initially as NH4HSO4

      PSI1LO = TINY                ! Low  limit
      PSI1HI = CHI1                ! High limit

      ! INITIAL VALUES FOR BISECTION 
      X1 = PSI1HI
      Y1 = FUNCJ1 (X1)
      YHI= Y1                      ! Save Y-value at HI position

      ! YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH42SO4 ****
      IF ( ABS(Y1) .LE. EPS .OR. YHI .LT. ZERO ) GOTO 50

      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1-DX
         Y2 = FUNCJ1 (X2)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y2) .LT. ZERO ) GOTO 20 ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH42SO4
      YLO= Y1                      ! Save Y-value at Hi position
      IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         Y3 = FUNCJ1 (ZERO)
         GOTO 50
      ELSE IF ( ABS(Y2) .LT. EPS ) THEN   ! X2 IS A SOLUTION 
         GOTO 50
      ELSE
         CALL PUSHERR (0001, 'CALCJ1')    ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF

      ! PERFORM BISECTION 
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCJ1 (X3)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO ) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS(X2-X1) .LE. EPS*X1 ) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCJ1') ! WARNING ERROR: NO CONVERGENCE

      ! CONVERGED ; RETURN
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCJ1 (X3)
 
 50   CONTINUE
      RETURN

      ! Return to calling program
      END SUBROUTINE CALCJ1 

!------------------------------------------------------------------------------

      FUNCTION FUNCJ1( P1 ) RESULT( VALUE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCJ1
! *** CASE J1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h"

      ! Arguments
      REAL*8, INTENT(IN)  :: P1

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     VALUE
      REAL*8  ::     BB, CC, DD

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3
      COMMON /CASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3
!$OMP THREADPRIVATE( /CASEJ/ )

      !=================================================================
      ! FUNCJ1 begins here!
      !=================================================================

      ! SETUP PARAMETERS      
      FRST   = .TRUE.
      CALAIN = .TRUE.

      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
      PSI1   = P1

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A1 = XK11 *(WATER/GAMA(12))**2.d0
         A2 = XK12 *(WATER/GAMA(09))**2.d0
         A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

         PSI2 = 0.5d0*(-(LAMDA+PSI1) + SQRT((LAMDA+PSI1)**2.D0+4.D0*A2)) ! PSI2
         PSI2 = MIN (PSI2, CHI2)

         BB   = A3+LAMDA                      ! KAPA
         CC   =-A3*(LAMDA + PSI2 + PSI1)
         DD   = BB*BB-4.D0*CC
         KAPA = 0.5D0*(-BB+SQRT(DD))    

         ! SAVE CONCENTRATIONS IN MOLAL ARRAY
         MOLAL(1) = LAMDA + KAPA                  ! HI
         MOLAL(2) = PSI1                          ! NAI
         MOLAL(3) = PSI2                          ! NH4I
         MOLAL(4) = ZERO
         MOLAL(5) = KAPA                          ! SO4I
         MOLAL(6) = LAMDA + PSI1 + PSI2 - KAPA    ! HSO4I
         MOLAL(7) = ZERO
         CALL CALCMR                              ! Water content

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE OBJECTIVE FUNCTION
 20   CONTINUE
      VALUE = MOLAL(2)*MOLAL(6)/A1 - ONE

      ! Return to calling program
      END FUNCTION FUNCJ1 

!------------------------------------------------------------------------------

      SUBROUTINE CALCJ2
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCJ2
! *** CASE J2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : NAHSO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     PSI1LO, PSI1HI, X1, Y1, DX
      REAL*8  ::     X2, Y2, YL0, Y3, X3, YHI, YLO

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3
      COMMON /CASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3
!$OMP THREADPRIVATE( /CASEJ/ )

      !=================================================================
      ! CALCJ2 begins here!
      !=================================================================

      ! SETUP PARAMETERS
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      CHI1   = W(1)                ! NA TOTAL
      CHI2   = W(3)                ! NH4 TOTAL
      PSI1LO = TINY                ! Low  limit
      PSI1HI = CHI1                ! High limit

      ! INITIAL VALUES FOR BISECTION
      X1 = PSI1HI
      Y1 = FUNCJ2 (X1)
      YHI= Y1                      ! Save Y-value at HI position

      ! YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH42SO4
      IF ( ABS(Y1) .LE. EPS .OR. YHI .LT. ZERO ) GOTO 50

      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
      DO 10 K=1,NDIV
         X2 = X1-DX
         Y2 = FUNCJ2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
 10   CONTINUE

      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH42SO4
      YLO= Y1                      ! Save Y-value at Hi position
      IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         Y3 = FUNCJ2 (ZERO)
         GOTO 50
      ELSE IF ( ABS(Y2) .LT. EPS ) THEN ! X2 IS A SOLUTION 
         GOTO 50
      ELSE
         CALL PUSHERR (0001, 'CALCJ2') ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF

      ! PERFORM BISECTION
 20   CONTINUE
      DO 30 K=1,MAXIT
         X3 = 0.5d0*(X1+X2)
         Y3 = FUNCJ2 (X3)
         IF ( SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO ) THEN ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS(X2-X1) .LE. EPS*X1 ) GOTO 40
 30   CONTINUE
      CALL PUSHERR (0002, 'CALCJ2') ! WARNING ERROR: NO CONVERGENCE

      ! CONVERGED ; RETURN
 40   CONTINUE
      X3 = 0.5d0*(X1+X2)
      Y3 = FUNCJ2 (X3)
      
 50   CONTINUE
      RETURN

      ! Return to calling program
      END SUBROUTINE CALCJ2 

!------------------------------------------------------------------------------

      FUNCTION FUNCJ2( P1 ) RESULT( VALUE )
! 
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE FUNCJ2
! *** CASE J2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Arguments
      REAL*8, INTENT(IN)  :: P1

      ! Local variables
      INTEGER ::     K
      REAL*8  ::     VALUE
      REAL*8  ::     BB, CC, DD

      ! Local common blocks
      REAL*8  ::     CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3
      COMMON /CASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3
!$OMP THREADPRIVATE( /CASEJ/ )

      !=================================================================
      ! FUNCJ2 begins here!
      !=================================================================

      ! SETUP PARAMETERS
      FRST   = .TRUE.
      CALAIN = .TRUE.

      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
      PSI1   = P1
      PSI2   = CHI2                           ! ALL NH4HSO4 DELIQUESCED

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A1 = XK11 *(WATER/GAMA(12))**2.d0
         A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

         ! CALCULATE DISSOCIATION QUANTITIES
         BB   = A3+LAMDA                          ! KAPA
         CC   =-A3*(LAMDA + PSI1 + PSI2)
         DD   = BB*BB-4.D0*CC
         KAPA = 0.5D0*(-BB+SQRT(DD))

         ! CALCULATE SPECIATION
         MOLAL(1) = LAMDA + KAPA                  ! HI
         MOLAL(2) = PSI1                          ! NAI
         MOLAL(3) = PSI2                          ! NH4I
         MOLAL(4) = ZERO                          ! CLI
         MOLAL(5) = KAPA                          ! SO4I
         MOLAL(6) = LAMDA + PSI1 + PSI2 - KAPA    ! HSO4I
         MOLAL(7) = ZERO                          ! NO3I
         CALL CALCMR                              ! Water content

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE

      ! CALCULATE OBJECTIVE FUNCTION
 20   CONTINUE
      VALUE = MOLAL(2)*MOLAL(6)/A1 - ONE

      ! Return to calling program
      END FUNCTION FUNCJ2 

!------------------------------------------------------------------------------

      SUBROUTINE CALCJ3
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCJ3
! *** CASE J3
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS ONLY A LIQUID PHASE
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      ! Local variables
      INTEGER :: K
      REAL*8  :: LAMDA, KAPA
      REAL*8  :: BB, CC, DD
      REAL*8  :: CHI1, CHI2, PSI1, PSI2, A3

      !=================================================================
      ! CALCJ3 begins here!
      !=================================================================

      ! SETUP PARAMETERS 
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.

      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
      CHI1   = W(1)                           ! NA TOTAL as NaHSO4
      CHI2   = W(3)                           ! NH4 TOTAL as NH4HSO4
      PSI1   = CHI1
      PSI2   = CHI2                           ! ALL NH4HSO4 DELIQUESCED

      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF.
      DO 10 K=1,NSWEEP

         A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0

         ! CALCULATE DISSOCIATION QUANTITIES
         BB   = A3+LAMDA                        ! KAPA
         CC   =-A3*(LAMDA + PSI1 + PSI2)
         DD   = BB*BB-4.D0*CC
         KAPA = 0.5D0*(-BB+SQRT(DD))

         ! CALCULATE SPECIATION
         MOLAL(1) = LAMDA + KAPA                 ! HI
         MOLAL(2) = PSI1                         ! NAI
         MOLAL(3) = PSI2                         ! NH4I
         MOLAL(4) = ZERO                         ! CLI
         MOLAL(5) = KAPA                         ! SO4I
         MOLAL(6) = LAMDA + PSI1 + PSI2 - KAPA   ! HSO4I
         MOLAL(7) = ZERO                         ! NO3I
         CALL CALCMR                             ! Water content

         ! CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOPC
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 50
         ENDIF
 10   CONTINUE

      ! Exit
 50   CONTINUE
      RETURN

      ! Return to calling program
      END SUBROUTINE CALCJ3 

!------------------------------------------------------------------------------

      SUBROUTINE CALCNHA
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNHA
!
!     THIS SUBROUTINE CALCULATES THE DISSOLUTION OF HCL, HNO3 AT
!     THE PRESENCE OF (H,SO4). HCL, HNO3 ARE CONSIDERED MINOR SPECIES,
!     THAT DO NOT SIGNIFICANTLY AFFECT THE EQUILIBRIUM POINT.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 


      INTEGER  :: ISLV
      REAL*8   :: DELCL, DELNO, OMEGA, C1, C2, C3
      REAL*8   :: M1, M2, M3
      REAL*8   :: A3, A4, CHI3, CHI4
      CHARACTER ERRINF*40 

      !=================================================================
      ! CALCNHA begins here!
      !=================================================================    

      ! SPECIAL CASE; WATER=ZERO 
      IF ( WATER .LE. TINY ) THEN
         GOTO 55

      ! SPECIAL CASE; HCL=HNO3=ZERO
      ELSE IF ( W(5) .LE. TINY .AND. W(4) .LE. TINY ) THEN
         GOTO 60

      ! SPECIAL CASE; HCL=ZERO
      ELSE IF ( W(5) .LE. TINY ) THEN
         CALL CALCNA              ! CALL HNO3 DISSOLUTION ROUTINE
         GOTO 60

      ! SPECIAL CASE; HNO3=ZERO 
      ELSE IF ( W(4) .LE. TINY ) THEN
         CALL CALCHA               ! CALL HCL DISSOLUTION ROUTINE
         GOTO 60
      ENDIF

      ! CALCULATE EQUILIBRIUM CONSTANTSC
      A3 = XK4*R*TEMP*(WATER/GAMA(10))**2.d0   ! HNO3
      A4 = XK3*R*TEMP*(WATER/GAMA(11))**2.d0   ! HCL

      ! CALCULATE CUBIC EQUATION COEFFICIENTS
      DELCL = ZERO
      DELNO = ZERO

      OMEGA = MOLAL(1)       ! H+
      CHI3  = W(4)           ! HNO3
      CHI4  = W(5)           ! HCL

      C1    = A3*CHI3
      C2    = A4*CHI4
      C3    = A3 - A4

      M1    = (C1 + C2 + (OMEGA+A4)*C3)/C3
      M2    = ((OMEGA+A4)*C2 - A4*C3*CHI4)/C3
      M3    =-A4*C2*CHI4/C3

      ! CALCULATE ROOTS
      CALL POLY3 (M1, M2, M3, DELCL, ISLV) ! HCL DISSOLUTION
      IF ( ISLV .NE. 0 ) THEN
         DELCL = TINY           ! TINY AMOUNTS OF HCL ASSUMED WHEN NO ROOT 
         WRITE (ERRINF,'(1PE7.1)') TINY
         CALL PUSHERR (0022, ERRINF) ! WARNING ERROR: NO SOLUTION
      ENDIF
      DELCL = MIN(DELCL, CHI4)

      DELNO = C1*DELCL/(C2 + C3*DELCL)  
      DELNO = MIN(DELNO, CHI3)

      IF ( DELCL .LT. ZERO .OR. DELNO .LT. ZERO .OR.
     &   DELCL .GT. CHI4 .OR. DELNO .GT. CHI3       ) THEN
         DELCL = TINY  ! TINY AMOUNTS OF HCL ASSUMED WHEN NO ROOT 
         DELNO = TINY
         WRITE (ERRINF,'(1PE7.1)') TINY
         CALL PUSHERR (0022, ERRINF)    ! WARNING ERROR: NO SOLUTION
      ENDIF

! Comment out
!C
!C *** COMPARE DELTA TO TOTAL H+ ; ESTIMATE EFFECT TO HSO4 ***************
!C
!      IF ((DELCL+DELNO)/MOLAL(1).GT.0.1d0) THEN
!         WRITE (ERRINF,'(1PE10.3)') (DELCL+DELNO)/MOLAL(1)*100.0
!         CALL PUSHERR (0021, ERRINF)   
!      ENDIF

      ! EFFECT ON LIQUID PHASE
 50   CONTINUE
      MOLAL(1) = MOLAL(1) + (DELNO+DELCL) ! H+   CHANGE
      MOLAL(4) = MOLAL(4) + DELCL          ! CL-  CHANGE
      MOLAL(7) = MOLAL(7) + DELNO          ! NO3- CHANGE

      ! EFFECT ON GAS PHASE
 55   CONTINUE
      GHCL     = MAX(W(5) - MOLAL(4), TINY)
      GHNO3    = MAX(W(4) - MOLAL(7), TINY)
      
      ! Exit
 60   CONTINUE
      RETURN

      ! Return to calling program
      END SUBROUTINE CALCNHA 
!------------------------------------------------------------------------------

      SUBROUTINE CALCHA
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCHA
! *** CALCULATES CHLORIDES SPECIATION
!
!     HYDROCHLORIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES,  
!     AND DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. THE 
!     HYDROCHLORIC ACID DISSOLVED IS CALCULATED FROM THE 
!     HCL(G) <-> (H+) + (CL-) 
!     EQUILIBRIUM, USING THE (H+) FROM THE SULFATES.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
#     include "isoropia.h" 

      REAL*8 :: X, DELT, ALFA, DIAK
      REAL*8 :: KAPA

      !=================================================================
      ! CALCHA begins here!
      !=================================================================

      ! CALCULATE HCL DISSOLUTION
      X    = W(5) 
      DELT = 0.0d0
      IF ( WATER .GT. TINY ) THEN
         KAPA = MOLAL(1)
         ALFA = XK3*R*TEMP*(WATER/GAMA(11))**2.d0
         DIAK = SQRT( (KAPA+ALFA)**2.d0 + 4.0*ALFA*X)
         DELT = 0.5d0*(-(KAPA+ALFA) + DIAK)
         ! Comment out
         !IF (DELT/KAPA.GT.0.1d0) THEN
         !   WRITE (ERRINF,'(1PE10.3)') DELT/KAPA*100.0
         !   CALL PUSHERR (0033, ERRINF)    
         !ENDIF
      ENDIF

      ! CALCULATE HCL SPECIATION IN THE GAS PHASE
      GHCL     = MAX(X-DELT, 0.0d0)  ! GAS HCL

      ! CALCULATE HCL SPECIATION IN THE LIQUID PHASE
      MOLAL(4) = DELT                ! CL-
      MOLAL(1) = MOLAL(1) + DELT     ! H+ 

      ! Return to calling program
      END SUBROUTINE CALCHA 

!------------------------------------------------------------------------------

      SUBROUTINE INIT_ISOROPIA
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE SETPARM
! *** THIS SUBROUTINE REDEFINES THE SOLUTION PARAMETERS OF ISORROPIA
!
! ======================== ARGUMENTS / USAGE ===========================
!
! *** NOTE: IF NEGATIVE VALUES ARE GIVEN FOR A PARAMETER, IT IS
!     IGNORED AND THE CURRENT VALUE IS USED INSTEAD.
! 
!  INPUT:
!  1. [WFTYPI] 
!     INTEGER variable.
!     Defines the type of weighting algorithm for the solution in Mutual 
!     Deliquescence Regions (MDR's):
!     0 - MDR's are assumed dry. This is equivalent to the approach 
!         used by SEQUILIB.
!     1 - The solution is assumed "half" dry and "half" wet throughout
!         the MDR.
!     2 - The solution is a relative-humidity weighted mean of the
!         dry and wet solutions (as defined in Nenes et al., 1998)
!
!  2. [IACALCI] 
!     INTEGER variable.
!     Method of activity coefficient calculation:
!     0 - Calculate coefficients during runtime
!     1 - Use precalculated tables
! 
!  3. [EPSI] 
!     DOUBLE PRECITION variable.
!     Defines the convergence criterion for all iterative processes
!     in ISORROPIA, except those for activity coefficient calculations
!     (EPSACTI controls that).
!
!  4. [MAXITI]
!     INTEGER variable.
!     Defines the maximum number of iterations for all iterative 
!     processes in ISORROPIA, except for activity coefficient calculations 
!     (NSWEEPI controls that).
!
!  5. [NSWEEPI]
!     INTEGER variable.
!     Defines the maximum number of iterations for activity coefficient 
!     calculations.
! 
!  6. [EPSACTI] 
!     DOUBLE PRECISION variable.
!     Defines the convergence criterion for activity coefficient 
!     calculations.
! 
!  7. [NDIV] 
!     INTEGER variable.
!     Defines the number of subdivisions needed for the initial root
!     tracking for the bisection method. Usually this parameter should 
!     not be altered, but is included for completeness.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"
#     include "isoropia.h" 

      ! Local variables
      INTEGER :: AS 

      !=================================================================
      ! INIT_ISOROPIA begins here!
      !=================================================================
      ALLOCATE( ASRAT( NASRD ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'ASRAT' )
      ASRAT = 0d0

      ALLOCATE( ASSO4( NSO4S ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'ASSO4' )
      ASSO4 = 0d0

      ALLOCATE( BNC198( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC198' )
      BNC198 = 0e0

      ALLOCATE( BNC223( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC223' )
      BNC223 = 0e0

      ALLOCATE( BNC248( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC248' )
      BNC248 = 0e0

      ALLOCATE( BNC273( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC273' )
      BNC273 = 0e0
     
      ALLOCATE( BNC298( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC298' )
      BNC298 = 0e0

      ALLOCATE( BNC323( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC323' )
      BNC323 = 0e0

      ALLOCATE( HNO3_sav( IIPAR, JJPAR, LLTROP ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HNO3_sav' )
      HNO3_sav = 0d0

      ALLOCATE( GAS_HNO3( IIPAR, JJPAR, LLTROP ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GAS_HNO3' )
      GAS_HNO3 = 0d0

      ! 
      ZZ      = (/1,2,1,2,1,1,2,1,1,1,1,1,2/)
      Z       = (/1,1,1,1,2,1,1/)

      !=================================================================
      ! ZSR RELATIONSHIP PARAMETERS
      !=================================================================

      ! AWAS = ammonium sulfate
      AWAS = (/
     & 100.d0,100.d0,100.d0,100.d0,100.d0,100.d0,
     $ 100.d0,100.d0,100.d0,100.d0,100.d0,
     & 100.d0,100.d0,100.d0,100.d0,100.d0,100.d0,
     & 100.d0,100.d0,100.d0,100.d0,100.d0,
     & 100.d0,100.d0,100.d0,100.d0,100.d0,100.d0,
     & 100.d0,100.d0,100.d0,100.d0,100.d0,
     & 30.d0, 30.d0, 30.d0, 29.54d0, 
     & 28.25d0, 27.06d0, 25.94d0,
     & 24.89d0, 23.9d0, 22.97d0, 22.1d0, 21.27d0, 20.48d0,
     & 19.73d0, 19.02d0, 18.34d0, 17.69d0,
     & 17.07d0, 16.48d0, 15.91d0, 15.37d0, 14.85d0,
     & 14.34d0, 13.86d0, 13.39d0, 12.94d0, 12.5d0,
     & 12.08d0, 11.67d0, 11.27d0, 10.88d0, 10.51d0, 10.14d0,
     &  9.79d0, 9.44d0, 9.1d0, 8.78d0,
     &  8.45d0, 8.14d0, 7.83d0, 7.53d0, 7.23d0,
     &  6.94d0, 6.65d0, 6.36d0, 6.08d0, 5.81d0,
     &  5.53d0, 5.26d0, 4.99d0, 4.72d0, 4.46d0, 
     &  4.19d0, 3.92d0, 3.65d0, 3.38d0, 3.11d0,
     &  2.83d0, 2.54d0, 2.25d0, 1.95d0, 1.63d0,
     &  1.31d0, 0.97d0, 0.63d0, 0.3d0, 0.001d0/)

      ! AWSN= sodium nitrate
      AWSN = (/
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,685.59d0,
     & 451.d0,336.46d0,268.48d0,223.41d0,191.28d0,
     & 167.2d0,148.46d0,133.44d0,121.12d0,110.83d0,
     & 102.09d0,94.57d0,88.03d0,82.29d0,77.2d0,
     & 72.65d0,68.56d0,64.87d0,61.51d0,58.44d0,
     & 55.62d0,53.03d0,50.63d0,48.4d0,46.32d0,44.39d0,
     & 42.57d0,40.87d0,39.27d0,37.76d0,
     & 36.33d0,34.98d0,33.7d0,32.48d0,31.32d0,
     & 30.21d0,29.16d0,28.14d0,27.18d0,26.25d0,
     & 25.35d0,24.5d0,23.67d0,22.87d0,22.11d0,21.36d0,
     & 20.65d0,19.95d0,19.28d0,18.62d0,
     & 17.99d0,17.37d0,16.77d0,16.18d0,15.61d0,15.05d0,
     & 14.51d0,13.98d0,13.45d0,12.94d0,
     & 12.44d0,11.94d0,11.46d0,10.98d0,10.51d0,10.04d0,
     &  9.58d0, 9.12d0, 8.67d0, 8.22d0,
     &  7.77d0, 7.32d0, 6.88d0, 6.43d0, 5.98d0, 5.53d0,
     &  5.07d0, 4.61d0, 4.15d0, 3.69d0,
     &  3.22d0, 2.76d0, 2.31d0, 1.87d0, 1.47d0, 1.1d0,
     &  0.77d0, 0.48d0, 0.23d0, 0.001d0/)

      ! AWSC = sodium chloride
      AWSC =(/
     &  100.d0, 100.d0, 100.d0, 100.d0, 100.d0,
     &  100.d0, 100.d0, 100.d0, 100.d0, 100.d0,
     &  100.d0, 100.d0, 100.d0, 100.d0, 100.d0,
     &  100.d0, 100.d0, 100.d0, 100.d0,16.34d0,
     & 16.28d0,16.22d0,16.15d0,16.09d0,16.02d0,
     & 15.95d0,15.88d0,15.8d0,15.72d0,15.64d0,
     & 15.55d0,15.45d0,15.36d0,15.25d0,15.14d0,
     & 15.02d0,14.89d0,14.75d0,14.6d0,14.43d0,
     & 14.25d0,14.04d0,13.81d0,13.55d0,13.25d0,
     & 12.92d0,12.56d0,12.19d0,11.82d0,11.47d0,
     & 11.13d0,10.82d0,10.53d0,10.26d0,10.d0, 
     &  9.76d0, 9.53d0, 9.3d0, 9.09d0, 8.88d0,
     &  8.67d0, 8.48d0, 8.28d0, 8.09d0, 7.9d0, 
     &  7.72d0, 7.54d0, 7.36d0, 7.17d0, 6.99d0,
     &  6.81d0, 6.63d0, 6.45d0, 6.27d0, 6.09d0,
     &  5.91d0, 5.72d0, 5.53d0, 5.34d0, 5.14d0,
     &  4.94d0, 4.74d0, 4.53d0, 4.31d0, 4.09d0,
     &  3.86d0, 3.62d0, 3.37d0, 3.12d0, 2.85d0,
     &  2.58d0, 2.3d0, 2.01d0, 1.72d0, 1.44d0,
     &  1.16d0, 0.89d0, 0.64d0, 0.4d0, 0.18d0/)

      ! AWAC = ammonium chloride
      AWAC = (/
     &  100.d0, 100.d0, 100.d0, 100.d0, 100.d0,
     &  100.d0, 100.d0, 100.d0, 100.d0, 100.d0,
     &  100.d0, 100.d0, 100.d0, 100.d0, 100.d0,
     &  100.d0, 100.d0, 100.d0, 100.d0,31.45d0,
     & 31.3d0,31.14d0,30.98d0,30.82d0,30.65d0,
     & 30.48d0,30.3d0,30.11d0,29.92d0,29.71d0,
     & 29.5d0,29.29d0,29.0d06,28.82d0,28.57d0,
     & 28.3d0,28.03d0,27.78d0,27.78d0,27.77d0,
     & 27.77d0,27.43d0,27.07d0,26.67d0,26.21d0,
     & 25.73d0,25.18d0,24.56d0,23.84d0,23.01d0,
     & 22.05d0,20.97d0,19.85d0,18.77d0,17.78d0,
     & 16.89d0,16.1d0,15.39d0,14.74d0,14.14d0,
     & 13.59d0,13.06d0,12.56d0,12.09d0,11.65d0,
     & 11.22d0,10.81d0,10.42d0,10.03d0, 9.66d0,
     &  9.3d0, 8.94d0, 8.59d0, 8.25d0, 7.92d0,
     &  7.59d0, 7.27d0, 6.95d0, 6.63d0, 6.32d0,
     &  6.01d0, 5.7d0, 5.39d0, 5.08d0, 4.78d0,
     &  4.47d0, 4.17d0, 3.86d0, 3.56d0, 3.25d0,
     &  2.94d0, 2.62d0, 2.3d0, 1.98d0, 1.65d0,
     &  1.32d0, 0.97d0, 0.62d0, 0.26d0, 0.13d0/)

      ! AWSS = sodium sulfate
      AWSS = (/
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,
     & 14.3d0,14.3d0,14.3d0,14.3d0,14.3d0,14.3d0,
     & 14.3d0,14.3d0,14.3d0,14.3d0,
     & 14.3d0,14.3d0,14.3d0,14.3d0,14.3d0,
     & 14.3d0,14.3d0,14.3d0,14.3d0,14.3d0,
     & 14.3d0,14.3d0,14.3d0,14.21d0,12.53d0,11.47d0,
     & 10.66d0,10.01d0, 9.46d0, 8.99d0, 8.57d0,
     &  8.19d0, 7.85d0, 7.54d0, 7.25d0, 6.98d0,
     &  6.74d0, 6.5d0, 6.29d0, 6.08d0, 5.88d0, 5.7d0,
     &  5.52d0, 5.36d0, 5.2d0, 5.04d0,
     &  4.9d0, 4.75d0, 4.54d0, 4.34d0, 4.14d0, 3.93d0,
     &  3.71d0, 3.49d0, 3.26d0, 3.02d0,
     &  2.76d0, 2.49d0, 2.2d0, 1.89d0, 1.55d0, 1.18d0,
     &  0.82d0, 0.49d0, 0.22d0, 0.001d0/)
 
      ! AWAB = ammonium bisulfate
      AWAB = (/356.45d0,296.51d0,253.21d0,220.47d0,194.85d0,
     & 174.24d0,157.31d0,143.16d0,131.15d0,120.82d0,
     & 111.86d0,103.99d0,97.04d0,90.86d0,85.31d0,
     & 80.31d0,75.78d0,71.66d0,67.9d0,64.44d0,
     &  61.25d0,58.31d0,55.58d0,53.04d0,50.68d0,
     &  48.47d0,46.4d0,44.46d0,42.63d0,40.91d0,
     &  39.29d0,37.75d0,36.3d0,34.92d0,33.61d0,
     &  32.36d0,31.18d0,30.04d0,28.96d0,27.93d0,
     &  26.94d0,25.99d0,25.08d0,24.21d0,23.37d0,
     &  22.57d0,21.79d0,21.05d0,20.32d0,19.63d0,
     &  18.96d0,18.31d0,17.68d0,17.07d0,16.49d0,
     &  15.92d0,15.36d0,14.83d0,14.31d0,13.8d0,
     &  13.31d0,12.83d0,12.36d0,11.91d0,11.46d0,
     &  11.03d0,10.61d0,10.2d0, 9.8d0, 9.41d0,
     &   9.02d0, 8.64d0, 8.28d0, 7.91d0, 7.56d0,
     &   7.21d0, 6.87d0, 6.54d0, 6.21d0, 5.88d0,
     &   5.56d0, 5.25d0, 4.94d0, 4.63d0, 4.33d0,
     &   4.03d0, 3.73d0, 3.44d0, 3.14d0, 2.85d0,
     &   2.57d0, 2.28d0, 1.99d0, 1.71d0, 1.42d0,
     &   1.14d0, 0.86d0, 0.57d0, 0.29d0, 0.001d0/)

      ! AWSA = sulfuric acid
      AWSA = (/
     & 34.d0,33.56d0,29.22d0,26.55d0,24.61d0,
     & 23.11d0,21.89d0,20.87d0,19.99d0,
     & 19.21d0,18.51d0,17.87d0,17.29d0,
     & 16.76d0,16.26d0,15.8d0,15.37d0,14.95d0,14.56d0,
     & 14.2d0,13.85d0,13.53d0,13.22d0,12.93d0,
     & 12.66d0,12.4d0,12.14d0,11.9d0,11.67d0,
     & 11.44d0,11.22d0,11.01d0,10.8d0,10.6d0,
     & 10.4d0,10.2d0,10.01d0,9.83d0,9.65d0,9.47d0,
     & 9.3d0,9.13d0,8.96d0,8.81d0,8.64d0,8.48d0,
     & 8.33d0,8.17d0,8.02d0,7.87d0,7.72d0,7.58d0,
     & 7.44d0,7.3d0,7.16d0,7.02d0,6.88d0,6.75d0,
     & 6.61d0,6.48d0,6.35d0,6.21d0,6.08d0,5.95d0,
     & 5.82d0,5.69d0,5.56d0,5.44d0,5.31d0,5.18d0,
     & 5.05d0,4.92d0,4.79d0,4.66d0,4.53d0,4.4d0,
     & 4.27d0,4.14d0,4.d0,3.87d0,3.73d0,3.6d0,
     & 3.46d0,3.31d0,3.17d0,3.02d0,2.87d0,2.72d0,
     & 2.56d0,2.4d0,2.23d0,2.05d0,1.87d0,1.68d0,
     & 1.48d0,1.27d0,1.05d0,0.807d0,0.552d0,0.281d0/)

      ! AWLC = (NH4)3H(SO4)2
      AWLC = (/
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,
     & 17.d0,16.5d0,15.94d0,15.31d0,14.71d0,14.14d0,
     & 13.6d0,13.08d0,12.59d0,12.12d0,11.68d0,
     & 11.25d0,10.84d0,10.44d0,10.07d0, 9.71d0,
     &  9.36d0, 9.02d0, 8.7d0, 8.39d0, 8.09d0, 7.8d0,
     &  7.52d0, 7.25d0, 6.99d0, 6.73d0,
     &  6.49d0, 6.25d0, 6.02d0, 5.79d0, 5.57d0,
     &  5.36d0, 5.15d0, 4.95d0, 4.76d0, 4.56d0,
     &  4.38d0, 4.2d0, 4.02d0, 3.84d0, 3.67d0,
     &  3.51d0, 3.34d0, 3.18d0, 3.02d0, 2.87d0,
     &  2.72d0, 2.57d0, 2.42d0, 2.28d0, 2.13d0,
     &  1.99d0, 1.85d0, 1.71d0, 1.57d0, 1.43d0,
     &  1.3d0, 1.16d0, 1.02d0, 0.89d0, 0.75d0,
     &  0.61d0, 0.46d0, 0.32d0, 0.16d0, 0.001d0/)

      ! AWAN = ammonium nitrate
      AWAN = (/
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,
     & 1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,1.d5,
     &       97.17d0,92.28d0,87.66d0,83.15d0,
     & 78.87d0,74.84d0,70.98d0,67.46d0,64.11d0,
     & 60.98d0,58.07d0,55.37d0,52.85d0,50.43d0,
     & 48.24d0,46.19d0,44.26d0,42.4d0,40.7d0,
     & 39.1d0,37.54d0,36.1d0,34.69d0,33.35d0,
     & 32.11d0,30.89d0,29.71d0,28.58d0,27.46d0,
     & 26.42d0,25.37d0,24.33d0,23.89d0,22.42d0,
     & 21.48d0,20.56d0,19.65d0,18.76d0,17.91d0,
     & 17.05d0,16.23d0,15.4d0,14.61d0,13.82d0,
     & 13.03d0,12.3d0,11.55d0,10.83d0,10.14d0,
     &  9.44d0, 8.79d0, 8.13d0, 7.51d0, 6.91d0,
     &  6.32d0, 5.75d0, 5.18d0, 4.65d0, 4.14d0,
     &  3.65d0, 3.16d0, 2.71d0, 2.26d0, 1.83d0,
     &  1.42d0, 1.03d0, 0.66d0, 0.3d0, 0.001d0/)

      ! AWSB = sodium bisulfate
      AWSB = (/ 173.72d0,156.88d0,142.8d0,130.85d0,120.57d0,
     &          111.64d0,103.8d0,96.88d0,90.71d0,85.18d0,
     & 80.2d0,75.69d0,71.58d0,67.82d0,64.37d0,61.19d0,
     & 58.26d0,55.53d0,53.d0,50.64d0,
     & 48.44d0,46.37d0,44.44d0,42.61d0,40.9d0,39.27d0,
     & 37.74d0,36.29d0,34.91d0,33.61d0,
     & 32.36d0,31.18d0,30.05d0,28.97d0,27.94d0,
     & 26.95d0,26.d0,25.1d0,24.23d0,23.39d0,
     & 22.59d0,21.81d0,21.07d0,20.35d0,19.65d0,
     & 18.98d0,18.34d0,17.71d0,17.11d0,16.52d0,
     & 15.95d0,15.4d0,14.87d0,14.35d0,13.85d0,13.36d0,
     & 12.88d0,12.42d0,11.97d0,11.53d0,
     & 11.1d0,10.69d0,10.28d0, 9.88d0, 9.49d0,
     &  9.12d0, 8.75d0, 8.38d0, 8.03d0, 7.68d0,
     &  7.34d0, 7.01d0, 6.69d0, 6.37d0, 6.06d0,
     &  5.75d0, 5.45d0, 5.15d0, 4.86d0, 4.58d0,
     &  4.3d0, 4.02d0, 3.76d0, 3.49d0, 3.23d0,
     &  2.98d0, 2.73d0, 2.48d0, 2.24d0, 2.01d0,
     &  1.78d0, 1.56d0, 1.34d0, 1.13d0, 0.92d0,
     &  0.73d0, 0.53d0, 0.35d0, 0.17d0, 0.001d0/)

      ASSO4 = (/ 1.0d-9, 2.5d-9, 5.0d-9, 7.5d-9, 1.0d-8,
     &           2.5d-8, 5.0d-8, 7.5d-8, 1.0d-7, 2.5d-7, 
     &           5.0d-7, 7.5d-7, 1.0d-6, 5.0d-6 /)

      ASRAT = (/
     & 1.020464d0,  0.9998130d0, 0.9960167d0, 
     & 0.9984423d0, 1.004004d0,  1.010885d0,  
     & 1.018356d0,  1.026726d0,  1.034268d0,
     & 1.043846d0,  1.052933d0,  1.062230d0,  
     & 1.062213d0,  1.080050d0,  1.088350d0,
     & 1.096603d0,  1.104289d0,  1.111745d0,  
     & 1.094662d0,  1.121594d0,  1.268909d0,
     & 1.242444d0,  1.233815d0,  1.232088d0,
     & 1.234020d0,  1.238068d0,  1.243455d0,
     & 1.250636d0,  1.258734d0,  1.267543d0,
     & 1.276948d0,  1.286642d0,  1.293337d0, 
     & 1.305592d0,  1.314726d0,  1.323463d0,  
     & 1.333258d0,  1.343604d0,  1.344793d0, 
     & 1.355571d0,  1.431463d0,  1.405204d0,
     & 1.395791d0,  1.393190d0,  1.394403d0, 
     & 1.398107d0,  1.403811d0,  1.411744d0,  
     & 1.420560d0,  1.429990d0,  1.439742d0, 
     & 1.449507d0,  1.458986d0,  1.468403d0, 
     & 1.477394d0,  1.487373d0,  1.495385d0,
     & 1.503854d0,  1.512281d0,  1.520394d0,
     & 1.514464d0,  1.489699d0,  1.480686d0, 
     & 1.478187d0,  1.479446d0,  1.483310d0, 
     & 1.489316d0,  1.497517d0,  1.506501d0, 
     & 1.515816d0,  1.524724d0,  1.533950d0,
     & 1.542758d0,  1.551730d0,  1.559587d0, 
     & 1.568343d0,  1.575610d0,  1.583140d0,  
     & 1.590440d0,  1.596481d0,  1.567743d0, 
     & 1.544426d0,  1.535928d0,  1.533645d0, 
     & 1.535016d0,  1.539003d0,  1.545124d0, 
     & 1.553283d0,  1.561886d0,  1.570530d0,
     & 1.579234d0,  1.587813d0,  1.595956d0, 
     & 1.603901d0,  1.611349d0,  1.618833d0, 
     & 1.625819d0,  1.632543d0,  1.639032d0, 
     & 1.645276d0,  1.707390d0,  1.689553d0,  
     & 1.683198d0,  1.681810d0,  1.683490d0, 
     & 1.687477d0,  1.693148d0,  1.700084d0,  
     & 1.706917d0,  1.713507d0,  1.719952d0, 
     & 1.726190d0,  1.731985d0,  1.737544d0, 
     & 1.742673d0,  1.747756d0,  1.752431d0, 
     & 1.756890d0,  1.761141d0,  1.765190d0,
     & 1.785657d0,  1.771851d0,  1.767063d0, 
     & 1.766229d0,  1.767901d0,  1.771455d0, 
     & 1.776223d0,  1.781769d0,  1.787065d0, 
     & 1.792081d0,  1.796922d0,  1.801561d0,  
     & 1.805832d0,  1.809896d0,  1.813622d0, 
     & 1.817292d0,  1.820651d0,  1.823841d0,  
     & 1.826871d0,  1.829745d0,  1.822215d0, 
     & 1.810497d0,  1.806496d0,  1.805898d0, 
     & 1.807480d0,  1.810684d0,  1.814860d0, 
     & 1.819613d0,  1.824093d0,  1.828306d0,
     & 1.832352d0,  1.836209d0,  1.839748d0, 
     & 1.843105d0,  1.846175d0,  1.849192d0, 
     & 1.851948d0,  1.854574d0,  1.857038d0, 
     & 1.859387d0,  1.844588d0,  1.834208d0,  
     & 1.830701d0,  1.830233d0,  1.831727d0, 
     & 1.834665d0,  1.838429d0,  1.842658d0,
     & 1.846615d0,  1.850321d0,  1.853869d0, 
     & 1.857243d0,  1.860332d0,  1.863257d0, 
     & 1.865928d0,  1.868550d0,  1.870942d0, 
     & 1.873208d0,  1.875355d0,  1.877389d0,
     & 1.899556d0,  1.892637d0,  1.890367d0,
     & 1.890165d0,  1.891317d0,  1.893436d0, 
     & 1.896036d0,  1.898872d0,  1.901485d0, 
     & 1.903908d0,  1.906212d0,  1.908391d0,  
     & 1.910375d0,  1.912248d0,  1.913952d0, 
     & 1.915621d0,  1.917140d0,  1.918576d0,  
     & 1.919934d0,  1.921220d0,  1.928264d0, 
     & 1.923245d0,  1.921625d0,  1.921523d0, 
     & 1.922421d0,  1.924016d0,  1.925931d0, 
     & 1.927991d0,  1.929875d0,  1.931614d0,
     & 1.933262d0,  1.934816d0,  1.936229d0, 
     & 1.937560d0,  1.938769d0,  1.939951d0, 
     & 1.941026d0,  1.942042d0,  1.943003d0, 
     & 1.943911d0,  1.941205d0,  1.937060d0,  
     & 1.935734d0,  1.935666d0,  1.936430d0, 
     & 1.937769d0,  1.939359d0,  1.941061d0,
     & 1.942612d0,  1.944041d0,  1.945393d0, 
     & 1.946666d0,  1.947823d0,  1.948911d0, 
     & 1.949900d0,  1.950866d0,  1.951744d0, 
     & 1.952574d0,  1.953358d0,  1.954099d0,
     & 1.948985d0,  1.945372d0,  1.944221d0, 
     & 1.944171d0,  1.944850d0,  1.946027d0, 
     & 1.947419d0,  1.948902d0,  1.950251d0, 
     & 1.951494d0,  1.952668d0,  1.953773d0,  
     & 1.954776d0,  1.955719d0,  1.956576d0, 
     & 1.957413d0,  1.958174d0,  1.958892d0,  
     & 1.959571d0,  1.960213d0,  1.977193d0, 
     & 1.975540d0,  1.975023d0,  1.975015d0, 
     & 1.975346d0,  1.975903d0,  1.976547d0, 
     & 1.977225d0,  1.977838d0,  1.978401d0,
     & 1.978930d0,  1.979428d0,  1.979879d0, 
     & 1.980302d0,  1.980686d0,  1.981060d0, 
     & 1.981401d0,  1.981722d0,  1.982025d0,  1.982312d0 /)

      ! Return to calling program
      END SUBROUTINE INIT_ISOROPIA 

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_ISOROPIA
!
!******************************************************************************
!  Subroutine CLEANUP_KMC deallocates all module arrays. (rjp, bmy, 9/23/02)
!******************************************************************************
!
      IF ( ALLOCATED( BNC198   ) ) DEALLOCATE( BNC198   )
      IF ( ALLOCATED( BNC223   ) ) DEALLOCATE( BNC223   )
      IF ( ALLOCATED( BNC248   ) ) DEALLOCATE( BNC248   )
      IF ( ALLOCATED( BNC273   ) ) DEALLOCATE( BNC273   )
      IF ( ALLOCATED( BNC298   ) ) DEALLOCATE( BNC298   )
      IF ( ALLOCATED( BNC323   ) ) DEALLOCATE( BNC323   )
      IF ( ALLOCATED( HNO3_sav ) ) DEALLOCATE( HNO3_sav )
      IF ( ALLOCATED( GAS_HNO3 ) ) DEALLOCATE( GAS_HNO3 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_ISOROPIA

!------------------------------------------------------------------------------

      END MODULE ISOROPIA_MOD
