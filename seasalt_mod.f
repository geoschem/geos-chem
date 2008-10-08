! $Id: seasalt_mod.f,v 1.14 2008/10/08 18:30:31 bmy Exp $
      MODULE SEASALT_MOD
!
!******************************************************************************
!  Module SEASALT_MOD contains arrays and routines for performing either a
!  coupled chemistry/aerosol run or an offline seasalt aerosol simulation.
!  Original code taken from Mian Chin's GOCART model and modified accordingly.
!  (bec, rjp, bmy, 6/22/00, 7/18/08)
!
!  Seasalt aerosol species: (1) Accumulation mode (usually 0.1 -  0.5 um)
!                           (2) Coarse mode       (usually 0.5 - 10.0 um)
!
!  NOTE: You can change the bin sizes for accumulation mode and coarse
!        mode seasalt in the "input.geos" file in v7-yy-zz and higher.
!
!  Module Variables:
!  ============================================================================
!  (1 ) DRYSALA  (INTEGER) : Drydep index for accumulation mode sea salt
!  (2 ) DRYSALC  (INTEGER) : Drydep index for coarse mode sea salt
!  (3 ) NSALT    (INTEGER) : Number of sea salt tracers
!  (4 ) IDDEP    (INTEGER) : Drydep index array for sea salt tracers
!  (5 ) REDGE    (REAL*8 ) : Array for edges of seasalt radius bins
!  (6 ) RMID     (REAL*8 ) : Array for centers of seasalt radius bins
!  (7 ) SRC      (REAL*8 ) : Array for baseline seasalt emission/bin [kg/m2]
!  (7 ) SRC_N    (REAL*8 ) : Array for baseline seasalt emission/bin [#/m2]
!  (8 ) SS_DEN   (REAL*8 ) : Sea salt density [kg/m3]
!  (9 ) ALK_EMIS (REAL*8 ) : Array for alkalinity [kg]
!  (10) N_DENS   (REAL*8 ) : Number density of seasalt emissions [#/m3]
!
!  Module Routines:
!  ============================================================================
!  (1 ) CHEMSEASALT        : Driver routine for sea salt loss processes
!  (2 ) WET SETTLING       : Routine which performs wet settling of sea salt
!  (3 ) DRY_DEPOSITION     : Routine which performs dry deposition of sea salt
!  (4 ) EMISSSEASALT       : Driver routine for sea salt emissions
!  (5 ) SRCSALT            : Updates surface mixing ratio for sea salt
!  (6 ) GET_ALK            : Gets the alkalinity of seasalt emissions
!  (6 ) INIT_SEASALT       : Allocates all module arrays
!  (7 ) CLEANUP_SEASALT    : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by "seasalt_mod.f":
!  ============================================================================
!  (1 ) dao_mod.f          : Module w/ arrays for GMAO met fields
!  (2 ) diag_mod.f         : Module w/ GEOS-CHEM diagnostic arrays
!  (3 ) drydep_mod.f       : Module w/ GEOS-CHEM drydep routines
!  (4 ) error_mod.f        : Module w/ I/O error and NaN check routines
!  (5 ) grid_mod.f         : Module w/ horizontal grid information
!  (6 ) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (7 ) pbl_mix_mod.f      : Module w/ routines for PBL height & mixing
!  (8 ) pressure_mod.f     : Module w/ routines to compute P(I,J,L)
!  (9 ) time_mod.f         : Module w/ routines to compute date & time
!  (10) tracer_mod.f       : Module w/ GEOS-CHEM tracer array STT etc.
!  (11) tracerid_mod.f     : Module w/ pointers to tracers & emissions 
!
!  References:
!  ============================================================================
!  (1 ) Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
!  (2 ) Gong, S., L. Barrie, and J.-P. Blanchet, "Modeling sea-salt
!        aerosols in the atmosphere. 1. Model development", J. Geophys. Res.,
!        v. 102, 3805-3818, 1997.
!
!  NOTES:  
!  (1 ) Now references "logical_mod.f" and "tracer_mod.f".  Comment out 
!        SS_SIZE, this has been replaced by SALA_REDGE_um and SALC_REDGE_um
!        from "tracer_mod.f".  Increased NR_MAX to 200. (bmy, 7/20/04)
!  (2 ) Added error check in EMISSSEASALT (bmy, 1/20/05)
!  (3 ) Now references "pbl_mix_mod.f" (bmy, 2/22/05)
!  (4 ) Added routine GET_ALK to account for alkalinity. (bec, bmy, 4/13/05)
!  (5 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (6 ) Now only call dry deposition routine if LDRYD=T (bec, bmy, 5/23/06)
!  (7 ) Remove unused variables from GET_ALK.  Also fixed variable declaration
!        bug in WET_SETTLING. (bec, bmy, 9/5/06)
!  (8 ) Extra error check for low RH in WET_SETTLING (phs, 6/11/08)
!  (9 ) Bug fix to remove a double-substitution in GET_ALK (bec, bmy, 7/18/08)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "seasalt_mod.f"
      !=================================================================

      ! Make everyting PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CHEMSEASALT
      PUBLIC :: EMISSSEASALT
      PUBLIC :: CLEANUP_SEASALT
      PUBLIC :: GET_ALK

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      INTEGER, PARAMETER   :: NSALT = 2
      INTEGER, PARAMETER   :: NR_MAX = 200
      INTEGER              :: DRYSALA, DRYSALC

      ! Arrays
      INTEGER              :: IDDEP(NSALT)
      REAL*8,  ALLOCATABLE :: REDGE(:,:)   
      REAL*8,  ALLOCATABLE :: RMID(:,:)
      REAL*8,  ALLOCATABLE :: SRC(:,:)   
      REAL*8,  ALLOCATABLE :: SRC_N(:,:)   
      REAL*8,  ALLOCATABLE :: ALK_EMIS(:,:,:,:)
      REAL*8,  ALLOCATABLE :: N_DENS(:,:,:,:)    
      REAL*8               :: SS_DEN(NSALT)    = (/ 2200.d0, 2200.d0 /)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CHEMSEASALT
!
!******************************************************************************
!  Subroutine CHEMSEASALT is the interface between the GEOS-CHEM main program 
!  and the seasalt chemistry routines that mostly calculates seasalt dry 
!  deposition (rjp, bmy, 1/24/02, 5/23/06)
!
!  NOTES:
!  (1 ) Now reference STT from "tracer_mod.f".  Now references LPRT from
!        "logical_mod.f" (bmy, 7/20/04)
!  (2 ) Now only call DRY_DEPOSITION if LDRYD=T (bec, bmy, 5/23/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPNAME, NUMDEP
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE LOGICAL_MOD,  ONLY : LPRT,    LDRYD
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSALA, IDTSALC

#     include "CMN_SIZE"     ! Size parameters 

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: N

      !=================================================================
      ! CHEMSEASALT begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Initialize (if necessary)
         CALL INIT_SEASALT

         ! Find drydep species in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE ( 'SALA' )
                  DRYSALA = N
               CASE ( 'SALC' )
                  DRYSALC = N
               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO

         ! Store in IDDEP array
         IDDEP(1) = DRYSALA
         IDDEP(2) = DRYSALC

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Maybe someday we should merge these two separate calculations
      ! into one (rjp, 4/3/04)
      !=================================================================

      !-------------------
      ! Accumulation mode
      !-------------------
      CALL WET_SETTLING( STT(:,:,:,IDTSALA), 1 )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSEASALT: WET_SET, Accum' )

      IF ( LDRYD ) THEN
         CALL DRY_DEPOSITION( STT(:,:,:,IDTSALA), 1 )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSEASALT: DRY_DEP, Accum' )
      ENDIF

      !-------------------
      ! Coarse mode
      !-------------------
      CALL WET_SETTLING( STT(:,:,:,IDTSALC), 2 )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSEASALT: WET_SET, Coarse' )

      IF ( LDRYD ) THEN
         CALL DRY_DEPOSITION( STT(:,:,:,IDTSALC), 2 )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSEASALT: DRY_DEP, Coarse')
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEMSEASALT

!------------------------------------------------------------------------------

      SUBROUTINE WET_SETTLING( TC, N )
!
!******************************************************************************
!  Subroutine WET_SETTLING performs wet settling of sea salt.
!  (bec, rjp, bmy, 4/20/04, 6/11/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Sea salt tracer [kg]
!  (2 ) N  (INTEGER) : N=1 is accum mode; N=2 is coarse mode
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Contains modified tracer
!
!  NOTES:
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (3 ) Bug fix: DTCHEM has to be REAL*8, not integer. (bmy, 9/7/06)
!  (4 ) Now limit relative humidity to [tiny(real*8),0.99] range for DLOG
!         argument (phs, 5/1/08)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : T, BXHEIGHT, RH
      USE DIAG_MOD,      ONLY : AD44
      USE DRYDEP_MOD,    ONLY : DEPSAV
      USE PRESSURE_MOD,  ONLY : GET_PCENTER
      USE TRACER_MOD,    ONLY : SALA_REDGE_um, SALC_REDGE_um, XNUMOL
      USE TRACERID_MOD,  ONLY : IDTSALA,       IDTSALC
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE GRID_MOD,      ONLY : GET_AREA_CM2

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_GCTM"      ! g0
#     include "CMN_DIAG"      ! ND44

      ! Argumetns
      INTEGER, INTENT(IN)    :: N
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,      J,     L
      REAL*8                 :: DELZ,   DELZ1, REFF,     DEN
      REAL*8                 :: P,      DP,    PDP,      TEMP        
      REAL*8                 :: CONST,  SLIP,  VISC,     FAC1
      REAL*8                 :: FAC2,   FLUX,  AREA_CM2, RHB
      REAL*8                 :: RCM,    RWET,  RATIO_R,  RHO
      REAL*8                 :: TOT1,   TOT2,  DTCHEM
      REAL*8                 :: VTS(LLPAR)  
      REAL*8                 :: TC0(LLPAR)
      
      ! Parameters
      REAL*8,  PARAMETER     :: C1 =  0.7674d0 
      REAL*8,  PARAMETER     :: C2 =  3.079d0 
      REAL*8,  PARAMETER     :: C3 =  2.573d-11
      REAL*8,  PARAMETER     :: C4 = -1.424d0

      !=================================================================
      ! WET_SETTLING begins here!
      !=================================================================
      
      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Sea salt density [kg/m3]
      DEN  = SS_DEN( N )

      ! Seasalt effective radius (i.e. midpt of radius bin) [m]
      SELECT CASE ( N )

         ! Accum mode
         CASE( 1 )
            REFF = 0.5d-6 * ( SALA_REDGE_um(1) + SALA_REDGE_um(2) )

         ! Coarse mode
         CASE( 2 ) 
            REFF = 0.5d-6 * ( SALC_REDGE_um(1) + SALC_REDGE_um(2) )
            
      END SELECT

      ! Sea salt radius [cm]
      RCM  = REFF * 100d0  

      ! Exponential factors
      FAC1 = C1 * ( RCM**C2 )
      FAC2 = C3 * ( RCM**C4 )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,     L,    VTS,  P,        TEMP, RHB,  RWET ) 
!$OMP+PRIVATE( RATIO_R, RHO,   DP,   PDP,  CONST,    SLIP, VISC, TC0  )
!$OMP+PRIVATE( DELZ,    DELZ1, TOT1, TOT2, AREA_CM2, FLUX             )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR       

         ! Initialize 
         DO L = 1, LLPAR
            VTS(L) = 0d0
         ENDDO

         ! Loop over levels
         DO L = 1, LLPAR

            ! Pressure at center of the level [kPa]
            P       = GET_PCENTER(I,J,L) * 0.1d0

            ! Temperature [K]
            TEMP    = T(I,J,L)

            ! Cap RH at 0.99 
            RHB     = MIN( 0.99d0, RH(I,J,L) * 1d-2 )

            ! Safety check (phs, 5/1/08)
            RHB     = MAX( TINY(RHB), RHB           )

            ! Aerosol growth with relative humidity in radius [m] 
            ! (Gerber, 1985)
            RWET    = 0.01d0*(FAC1/(FAC2-DLOG(RHB))+RCM**3.d0)**0.33d0

            ! Ratio dry over wet radii at the cubic power
            RATIO_R = ( REFF / RWET )**3.d0

            ! Density of the wet aerosol (kg/m3)
            RHO     = RATIO_R * DEN + ( 1.d0 - RATIO_R ) * 1000.d0

            ! Dp = particle diameter [um]
            DP      = 2.d0 * RWET * 1.d6        

            ! PdP = P * dP [hPa * um]
            PDp     = P * Dp

            ! Constant
            CONST   = 2.d0 * RHO * RWET**2 * g0 / 9.d0

            !===========================================================
            ! NOTE: Slip correction factor calculations following 
            ! Seinfeld, pp464 which is thought to be more accurate 
            ! but more computation required. (rjp, 1/24/02)
            !
            ! # air molecule number density
            ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
            !
            ! # gas mean free path
            ! lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
            !
            ! # Slip correction
            ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp     
            !     &     / (2. * lamda))) / Dp
            !
            ! NOTE: Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
            ! which produces slip correction factore with small error
            ! compared to the above with less computation.
            !===========================================================  
          
            ! Slip correction factor (as function of P*dp)
            Slip = 1.d0+(15.60d0 + 7.0d0 * EXP(-0.059d0 * PDp)) / PDp

            ! Viscosity [Pa*s] of air as a function of temperature 
            VISC = 1.458d-6 * (Temp)**(1.5d0) / ( Temp + 110.4d0 )

            ! Settling velocity [m/s]
            VTS(L) = CONST * Slip / VISC
         ENDDO

         ! Method is to solve bidiagonal matrix which is
         ! implicit and first order accurate in z (rjp, 1/24/02)

         ! Save initial tracer concentration in column
         DO L = 1, LLPAR
            TC0(L) = TC(I,J,L)
         ENDDO

         ! We know the boundary condition at the model top
         L    = LLTROP
         DELZ = BXHEIGHT(I,J,L)

         TC(I,J,L) = TC(I,J,L) / ( 1.d0 + DTCHEM * VTS(L) / DELZ )

         DO L = LLTROP-1, 1, -1
            DELZ  = BXHEIGHT(I,J,L)
            DELZ1 = BXHEIGHT(I,J,L+1)
            TC(I,J,L) = 1.d0 / ( 1.d0 + DTCHEM * VTS(L) / DELZ )
     &                * ( TC(I,J,L) + DTCHEM * VTS(L+1) / DELZ1
     &                *  TC(I,J,L+1) )
         ENDDO
         
         !==============================================================
         ! ND44 diagnostic: sea salt loss [molec/cm2/s]
         !==============================================================
         IF ( ND44 > 0 ) THEN

            ! Initialize
            TOT1 = 0d0
            TOT2 = 0d0
            
            ! Compute column totals of TCO(:) and TC(I,J,:,N)
            DO L = 1, LLPAR
               TOT1 = TOT1 + TC0(L)
               TOT2 = TOT2 + TC(I,J,L)
            ENDDO

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            ! Convert sea salt flux from [kg/s] to [molec/cm2/s]
            FLUX     = ( TOT1 - TOT2 ) / DTCHEM
            FLUX     = FLUX * XNUMOL(IDTSALA) / AREA_CM2 
   
            ! Store in AD44 array
            AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE WET_SETTLING

!------------------------------------------------------------------------------

      SUBROUTINE DRY_DEPOSITION( TC, N )
!
!******************************************************************************
!  Subroutine DRY_DEPOSITION computes the loss of sea salt by dry deposition
!  at the surface, using an implicit method. (bec, rjp, bmy, 4/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Sea salt tracer [kg]
!  (2 ) N  (INTEGER) : N=1 is accum mode; N=2 is coarse mode
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Contains modified tracer
!
!  NOTES:
!  (1 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV 
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTSALA,   IDTSALC
      USE TIME_MOD,     ONLY : GET_MONTH, GET_TS_CHEM
      USE GRID_MOD,     ONLY : GET_AREA_CM2

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_GCTM"     ! g0
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      INTEGER, INTENT(IN)    :: N
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,        J,     L,      DTCHEM
      REAL*8                 :: OLD,      NEW,   G,      REFF
      REAL*8                 :: DIAM,     U_TS0, REYNOL, ALPHA 
      REAL*8                 :: BETA,     GAMMA, DENS,   FLUX 
      REAL*8                 :: AREA_CM2, TOT1,  TOT2

      ! Parameters
      REAL*8,  PARAMETER     :: RHOA = 1.25d-3

      !=================================================================
      ! DRY_DEPOSITION begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, OLD, NEW, FLUX )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Old tracer concentration [kg]
            OLD  = TC(I,J,1)

            ! New tracer concentration [kg]
            NEW  = OLD * EXP( -DEPSAV(I,J,IDDEP(N)) * DTCHEM  )

            !===========================================================
            ! ND44 diagnostic: sea salt drydep loss [molec/cm2/s]
            !===========================================================
            IF ( ND44 > 0 ) THEN
               
               ! Convert drydep loss from [kg/s] to [molec/cm2/s]
               FLUX = ( OLD - NEW ) / DTCHEM 
               FLUX = FLUX * XNUMOL(IDTSALA) / AREA_CM2 
            
               ! Store in AD44
               AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
            ENDIF

            ! Update tracer array
            TC(I,J,1) = NEW 
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DRY_DEPOSITION

!------------------------------------------------------------------------------

      SUBROUTINE EMISSSEASALT
!
!******************************************************************************
!  Subroutine EMISSSEASALT is the interface between the GEOS-CHEM model
!  and the SEASALT emissions routines in "seasalt_mod.f".
!  (bec, rjp, bmy, 3/24/03, 2/22/05)
!
!  NOTES:
!  (1 ) Now references LPRT from "logical_mod.f" and STT from "tracer_mod.f".
!        (bmy, 7/20/04)
!  (2 ) Now make sure IDTSALA, IDTSALC are nonzero before calling SRCSALT.
!        (bmy, 1/26/05)
!  (3 ) Remove reference to header file "CMN" (bmy, 2/22/05)
!  (4 ) Now call INIT_SEASALT on the first timestep.  Also initialize ALK_EMIS
!        and N_DENS on each timestep. (bec, bmy, 4/13/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE LOGICAL_MOD,  ONLY : LPRT
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSALA, IDTSALC

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: I, J, L, N

      !=================================================================
      ! EMISSSEASALT begins here! 
      !=================================================================

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### in EMISSEASALT' )

      ! Allocate all module arrays (bec, bmy, 4/13/05)
      IF ( FIRST ) THEN  
         CALL INIT_SEASALT
         FIRST = .FALSE.
      ENDIF

      ! Initialize for each timestep (bec, bmy, 4/13/05)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, NSALT
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ALK_EMIS(I,J,L,N) = 0d0
         N_DENS(I,J,L,N)   = 0d0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Accumulation mode emissions
      IF ( IDTSALA > 0 ) THEN
         CALL SRCSALT( STT(:,:,:,IDTSALA), 1 )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSEASALT: Accum' )
      ENDIF

      ! Coarse mode emissions
      IF ( IDTSALC > 0 ) THEN
         CALL SRCSALT( STT(:,:,:,IDTSALC), 2 )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSEASALT: Coarse' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE EMISSSEASALT

!-----------------------------------------------------------------------------

      SUBROUTINE SRCSALT( TC, N )
!
!******************************************************************************
!  Subroutine SRCSALT updates the surface mixing ratio of dry sea salt
!  aerosols for NSALT size bins.  The generation of sea salt aerosols
!  has been parameterized following Monahan et al. [1986] parameterization
!  as described by Gong et al. [1997].  (bec, rjp, bmy, 4/20/04, 10/25/05)
! 
!  Contact: Becky Alexander (bec@io.harvard.edu) or 
!           Rokjin Park     (rjp@io.harvard.edu)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Sea salt tracer array [v/v]
!  (2 ) N  (INTEGER) : N=1 denotes accumulation mode; N=2 denotes coarse mode
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Contains modified sea salt concentration [v/v]
!
!  References:
!  ============================================================================
!  (1 ) Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
!  (2 ) Gong, S., L. Barrie, and J.-P. Blanchet, "Modeling sea-salt
!        aerosols in the atmosphere. 1. Model development", J. Geophys. Res.,
!        v. 102, 3805-3818, 1997.
! 
!  NOTES:
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_TOP_L from "pbl_mix_mod.f".
!        Removed reference to header file CMN.  Removed reference to 
!        "pressure_mod.f".  (bmy, 2/22/05)
!  (3 ) Now also compute alkalinity and number density of seasalt emissions.
!        (bec, bmy, 4/13/05)
!  (4 ) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : PBL, AD, IS_WATER, AIRVOL
      USE DIAG_MOD,      ONLY : AD08
      USE ERROR_MOD,     ONLY : DEBUG_MSG, ERROR_STOP
      USE GRID_MOD,      ONLY : GET_AREA_M2
      USE PBL_MIX_MOD,   ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,      ONLY : GET_TS_EMIS
      USE TRACER_MOD,    ONLY : SALA_REDGE_um, SALC_REDGE_um, XNUMOL

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DIAG"      ! ND44, ND08
#     include "CMN_GCTM"      ! PI

      ! Arguments
      INTEGER, INTENT(IN)    :: N
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables 
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I,     J,      L
      INTEGER                :: R,     NR,     NTOP
      REAL*8                 :: W10M,  DTEMIS, R0
      REAL*8                 :: R1,    CONST,  CONST_N
      REAL*8                 :: FEMIS, A_M2
      REAL*8                 :: SALT(IIPAR,JJPAR)
      REAL*8                 :: SALT_N(IIPAR,JJPAR)

      ! Increment of radius for Emission integration (um)
      REAL*8, PARAMETER      :: DR    = 5.d-2
      REAL*8, PARAMETER      :: BETHA = 1.d0

      ! External functions
      REAL*8,  EXTERNAL      :: SFCWINDSQR

      !=================================================================
      ! SRCSALT begins here!
      !=================================================================
      
      ! Emission timestep [s]
      DTEMIS = GET_TS_EMIS() * 60d0

      ! Constant [volume * time * other stuff??] 
      CONST   = 4d0/3d0 * PI * DR * DTEMIS * 1.d-18 * 1.373d0

      ! Constant for converting to [#/m2] (bec, bmy, 4/13/05)
      CONST_N = DTEMIS * DR * 1.373d0

      ! Lower and upper limit of size bin N [um]
      SELECT CASE( N ) 
       
         ! Accum mode
         CASE( 1 )
            R0 = SALA_REDGE_um(1)
            R1 = SALA_REDGE_um(2)
          
         ! Coarse mode
         CASE( 2 )
            R0 = SALC_REDGE_um(1)
            R1 = SALC_REDGE_um(2)
            
      END SELECT

      ! Number of radius size bins
      NR = INT( ( ( R1 - R0 ) / DR ) + 0.5d0 ) 

      ! Error check
      IF ( NR > NR_MAX ) THEN
         CALL ERROR_STOP( 'Too many bins!', 'SRCSALT (seasalt_mod.f)' )
      ENDIF

      ! Initialize source
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         SALT(I,J)   = 0d0
         SALT_N(I,J) = 0d0
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Define edges and midpoints of each incrmental radius bin
      ! This only has to be done once per sea salt type
      !=================================================================
      IF ( FIRST ) THEN 

         ! Lower edge of 0th bin
         REDGE(0,N) = R0
      
         ! Loop over the # of radius bins
         DO R = 1, NR

            ! Midpoint of IRth bin
            RMID(R,N)  = REDGE(R-1,N) + ( DR / 2d0 )

            ! Upper edge of IRth bin
            REDGE(R,N) = REDGE(R-1,N) + DR 

            ! Sea salt base source [kg/m2]
            SRC(R,N)  = CONST * SS_DEN( N ) 
     &           * ( 1.d0 + 0.057d0*( BETHA * RMID(R,N) )**1.05d0 )
     &           * 10d0**( 1.19d0*
     &                  EXP(-((0.38d0-LOG(BETHA*RMID(R,N)))/0.65d0)**2))
     &           / BETHA**2         

            ! Sea salt base source [#/m2] (bec, bmy, 4/13/05)
            SRC_N(R,N) = CONST_N * (1.d0/RMID(R,N)**3)
     &           * (1.d0+0.057d0*(BETHA*RMID(R,N))**1.05d0)
     &           * 10d0**(1.19d0*EXP(-((0.38d0-LOG(BETHA*RMID(R,N)))
     &           /0.65d0)**2))/ BETHA**2  

!### Debug
!###           WRITE( 6, 100 ) R,REDGE(R-1,N),RMID(R,N),REDGE(R,N),SRC(R,N)
!### 100        FORMAT( 'IR, R0, RMID, R1: ', i3, 3f11.4,2x,es13.6 )
         ENDDO
      
         ! Reset only after N=NSALT
         IF ( FIRST .and.  N == NSALT ) FIRST = .FALSE.
      ENDIF
    
      !=================================================================
      ! Emission is integrated over a given size range for each bin
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, R, A_M2, W10M )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO J = 1, JJPAR
            
         ! Grid box surface area [m2]
         A_M2 = GET_AREA_M2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Test if this is a water box
            IF ( IS_WATER(I,J) ) THEN

               ! Wind speed at 10 m altitude [m/s]
               W10M = SQRT( SFCWINDSQR(I,J) )

               ! Loop over size bins
               DO R = 1, NR

                  ! Update seasalt source into SALT [kg]
                  SALT(I,J)   = SALT(I,J) + 
     &                          ( SRC(R,N)   * A_M2 * W10M**3.41d0 )

                  ! Update seasalt source into SALT_N [#] 
                  ! (bec, bmy, 4/13/05)
                  SALT_N(I,J) = SALT_N(I,J) +
     &                          ( SRC_N(R,N) * A_M2 * W10M**3.41d0 )

               ENDDO
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Now partition seasalt emissions through boundary layer
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, L, FEMIS )
!$OMP+SCHEDULE( DYNAMIC ) 
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Layer in which the PBL top occurs
         NTOP = CEILING( GET_PBL_TOP_L( I, J ) )
        
         ! Loop thru the boundary layer
         DO L = 1, NTOP

            ! Fraction of the PBL spanned by box (I,J,L) [unitless]
            FEMIS             = GET_FRAC_OF_PBL( I, J, L )

            ! Add seasalt emissions into box (I,J,L) [kg]
            TC(I,J,L)         = TC(I,J,L) + ( FEMIS * SALT(I,J) )

	    ! Alkalinity [kg] (bec, bmy, 4/13/05)
            ALK_EMIS(I,J,L,N) = SALT(I,J)

	    ! Number density [#/m3] (bec, bmy, 4/13/05)
	    N_DENS(I,J,L,N)   = SALT_N(I,J) / AIRVOL(I,J,L)              

         ENDDO

         ! ND08 diagnostic: sea salt emissions [kg]
         IF ( ND08 > 0 ) THEN
            AD08(I,J,N) = AD08(I,J,N) + SALT(I,J)
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      ! Return to calling program
      END SUBROUTINE SRCSALT

!------------------------------------------------------------------------------

      SUBROUTINE GET_ALK( I, J, L, ALK1, ALK2, Kt1, Kt2, Kt1N, Kt2N )
!
!******************************************************************************
!  Function GET_ALK returns the seasalt alkalinity emitted at each timestep to
!  sulfate_mod.f for chemistry on seasalt aerosols.
!  (bec, 12/7/04, 9/5/06)
! 
!  Arguments as Input:
!  ============================================================================
!
!  NOTES:
!  (1 ) Becky Alexander says we can remove AREA1, AREA2 (bec, bmy, 9/5/06)
!  (2 ) Bug fix to remove a double-substitution.  Replace code lines for 
!        TERM{123}A, TERM{123}B, TERM{123}AN, TERM{123}BN. (bec, bmy, 7/18/08)
!******************************************************************************
!
      USE DAO_MOD,      ONLY : AD, RH
      USE ERROR_MOD,    ONLY : IT_IS_NAN
      USE TRACER_MOD,   ONLY : SALA_REDGE_um, SALC_REDGE_um

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L 
      
      ! Return value
      REAL*8, INTENT(OUT)  :: ALK1, ALK2 ! [kg]
      REAL*8, INTENT(OUT)  :: Kt1, Kt2, Kt1N, Kt2N ! [s-1]
 
      REAL*8,  PARAMETER :: PI = 3.14159265
      REAL*8             :: N1, N2, Kt
      REAL*8             :: HGF, ALK
      REAL*8             :: RAD1, RAD2, RAD3
      REAL*8             :: term1a, term2a, term3a
      REAL*8             :: term1b, term2b, term3b
      REAL*8             :: term1aN, term2aN, term3aN
      REAL*8             :: term1bN, term2bN, term3bN
      REAL*8             :: const1, const2, const1N, const2N
      REAL*8             :: a1, a2, b1, b2, a1N, a2N, b1N, b2N
      REAL*8,  PARAMETER :: MINDAT = 1.d-20
      INTEGER            :: IRH
      REAL*8,  PARAMETER   :: gamma_SO2 = 0.11d0 !from Worsnop et al. (1989)
      REAL*8,  PARAMETER   :: gamma_HNO3 = 0.2d0 !from JPL [2001] 
      REAL*8,  PARAMETER   :: Dg = 0.2d0 !gas phase diffusion coeff. [cm2/s]
      REAL*8,  PARAMETER   :: v = 3.0d4 !cm/s

      LOGICAL, SAVE :: FIRST = .TRUE.
 
      !=================================================================
      ! GET_ALK begins here!
      !=================================================================

      ! Zero variables
      ALK1  = 0.D0
      ALK2  = 0.D0
      KT1   = 0.D0
      KT2   = 0.D0
      KT1N  = 0.D0
      KT2N  = 0.D0
      N1    = 0.D0
      N2    = 0.D0

      ! [kg] use this when not transporting alk
      ALK1  = ALK_EMIS(I,J,L,1)
      ALK2  = ALK_EMIS(I,J,L,2)

      !-----------------------------------------------------------------------
      ! NOTE: If you want to transport alkalinity then uncomment this section
      ! (bec, bmy, 4/13/05)
      ! 
      !! alkalinity [v/v] to [kg] use this when transporting alk
      !! or using Liao et al [2004] assumption of a continuous supply of
      ! alkalinity based on Laskin et al. [2003]
      !ALK1 = STT(I,J,L,IDTSALA) * AD(I,J,L)/TCVV(IDTSALA)
      !ALK2 = STT(I,J,L,IDTSALC) * AD(I,J,L)/TCVV(IDTSALC)
      !-----------------------------------------------------------------------

      ! Conversion from [m-3] --> [cm-3]
      N1 = N_DENS(I,J,L,1) * 1.d-6
      N2 = N_DENS(I,J,L,2) * 1.d-6

      ALK = ALK1 + ALK2

      ! If there is any alkalinity ...
      IF ( ALK > MINDAT ) THEN

         ! set humidity index IRH as a percent
         IRH = RH(I,J,L)
         IRH = MAX(  1, IRH )
         IRH = MIN( 99, IRH )

         ! hygroscopic growth factor for sea-salt from Chin et al. (2002)
         IF ( IRH < 100 ) HGF = 2.2d0
         IF ( IRH < 99  ) HGF = 1.9d0
         IF ( IRH < 95  ) HGF = 1.8d0
         IF ( IRH < 90  ) HGF = 1.6d0
         IF ( IRH < 80  ) HGF = 1.5d0
         IF ( IRH < 70  ) HGF = 1.4d0
         IF ( IRH < 50  ) HGF = 1.0d0
	
         ! radius of sea-salt aerosol size bins [cm] accounting for 
         ! hygroscopic growth
         RAD1 = SALA_REDGE_um(1) * HGF * 1.d-4 
         RAD2 = SALA_REDGE_um(2) * HGF * 1.d-4 
         RAD3 = SALC_REDGE_um(2) * HGF * 1.d-4 

         !----------------------------------
         ! SO2 uptake onto fine particles 
         !----------------------------------

	 ! calculate gas-to-particle rate constant for uptake of 
	 ! SO2 onto fine sea-salt aerosols [Jacob, 2000] analytical solution
         CONST1 = 4.D0/(V*GAMMA_SO2)
         A1     = (RAD1/DG)+CONST1
         B1     = (RAD2/DG)+CONST1
!-----------------------------------------------------------------------------
! Prior to 7/18/08:
! Becky Alexander's fix to remove double-substitution (bec, bmy, 7/18/08)
!         TERM1A = (((B1/DG)**2)+(2.0D0*CONST1*B1/DG)+(CONST1**2)) -
!     &            (((A1/DG)**2)+(2.0D0*CONST1*A1/DG)+(CONST1**2))
!         TERM2A = 2.D0*CONST1*(((B1/DG)+CONST1)-((A1/DG)+CONST1))
!         TERM3A = (CONST1**2)*(LOG((B1/DG)+CONST1) -
!     &            LOG((A1/DG)+CONST1))
!         KT1    = 4.D0*PI*N1*(DG**2)*(TERM1A - TERM2A + TERM3A)
!-----------------------------------------------------------------------------
         TERM1A = ((B1**2)/2.0d0) - ((A1**2)/2.0d0)
         TERM2A = 2.D0*CONST1*(B1-A1)
         TERM3A = (CONST1**2)*LOG(B1/A1)
         KT1    = 4.D0*PI*N1*(DG**3)*(TERM1A - TERM2A + TERM3A)

         !----------------------------------
         ! SO2 uptake onto coarse particles 
         !----------------------------------
         
	 ! calculate gas-to-particle rate constant for uptake of 
	 ! SO2 onto coarse sea-salt aerosols [Jacob, 2000] analytical solution
         CONST2 = 4.D0/(V*GAMMA_SO2)
         A2     = (RAD2/DG)+CONST2
         B2     = (RAD3/DG)+CONST2
!------------------------------------------------------------------------------
! Prior to 7/18/08:
! Becky Alexander's fix to remove double-substitution (bec, bmy, 7/18/08)
! Remove these lines:
!         TERM1B = (((B2/DG)**2)+(2.0D0*CONST2*B2/DG)+(CONST2**2)) -
!     &            (((A2/DG)**2)+(2.0D0*CONST2*A2/DG)+(CONST2**2))
!         TERM2B = 2.D0*CONST2*(((B2/DG)+CONST2)-((A2/DG)+CONST2))
!         TERM3B = (CONST2**2)*(LOG((B2/DG)+CONST2) -
!     &             LOG((A2/DG)+CONST2))
!         KT2    = 4.D0*PI*N2*(DG**2)*(TERM1B - TERM2B + TERM3B)
!------------------------------------------------------------------------------
         TERM1B = ((B2**2)/2.0d0) - ((A2**2)/2.0d0)
         TERM2B = 2.D0*CONST2*(B2-A2)
         TERM3B = (CONST2**2)*LOG(B2/A2)
         KT2    = 4.D0*PI*N2*(DG**3)*(TERM1B - TERM2B + TERM3B)
         KT     = KT1 + KT2

         !----------------------------------
         ! HNO3 uptake onto fine particles 
         !----------------------------------

         ! calculate gas-to-particle rate constant for uptake of 
         ! HNO3 onto fine sea-salt aerosols [Jacob, 2000] analytical solution
         CONST1N = 4.D0/(V*GAMMA_HNO3)
         A1N     = (RAD1/DG)+CONST1N
         B1N     = (RAD2/DG)+CONST1N
!-----------------------------------------------------------------------------
! Prior to 7/18/08:
! Becky Alexander's fix to remove double-substitution (bec, bmy, 7/18/08)
! Remove these lines:
!         TERM1AN = (((B1N/DG)**2)+(2.0D0*CONST1N*B1N/DG)+(CONST1N**2)) -
!     &             (((A1N/DG)**2)+(2.0D0*CONST1N*A1N/DG)+(CONST1N**2))
!         TERM2AN = 2.D0*CONST1N*(((B1N/DG)+CONST1N)-((A1N/DG)+CONST1N))
!         TERM3AN = (CONST1N**2)*(LOG((B1N/DG)+CONST1N) -
!     &             LOG((A1N/DG)+CONST1N))
!         KT1N    = 4.D0*PI*N1*(DG**2)*(TERM1AN - TERM2AN + TERM3AN)
!-----------------------------------------------------------------------------
         TERM1AN = ((B1N**2)/2.0d0) - ((A1N**2)/2.0d0)
         TERM2AN = 2.D0*CONST1N*(B1N-A1N)
         TERM3AN = (CONST1N**2)*LOG(B1N/A1N)
         KT1N    = 4.D0*PI*N1*(DG**3)*(TERM1AN - TERM2AN + TERM3AN)

         !----------------------------------
         ! HNO3 uptake onto coarse particles 
         !----------------------------------

	 ! calculate gas-to-particle rate constant for uptake of 
	 ! HNO3 onto coarse sea-salt aerosols [Jacob, 2000] analytical solution
         CONST2N = 4.D0/(V*GAMMA_HNO3)
         A2N     = (RAD2/DG)+CONST2N
         B2N     = (RAD3/DG)+CONST2N
!-----------------------------------------------------------------------------
! Prior to 7/18/08:
! Becky Alexander's fix to remove double-substitution (bec, bmy, 7/18/08)
! Remove these lines:
!         TERM1BN = (((B2N/DG)**2)+(2.0D0*CONST2N*B2N/DG)+(CONST2N**2)) -
!     &             (((A2N/DG)**2)+(2.0D0*CONST2N*A2N/DG)+(CONST2N**2))
!         TERM2BN = 2.D0*CONST2N*(((B2N/DG)+CONST2N)-((A2N/DG)+CONST2N))
!         TERM3BN = (CONST2N**2)*(LOG((B2N/DG)+CONST2N) -
!     &             LOG((A2N/DG)+CONST2N))
!         KT2N    = 4.D0*PI*N2*(DG**2)*(TERM1BN - TERM2BN + TERM3BN)
!-----------------------------------------------------------------------------
         TERM1BN = ((B2N**2)/2.0d0) - ((A2N**2)/2.0d0)
         TERM2BN = 2.D0*CONST2N*(B2N-A2N)
         TERM3BN = (CONST2N**2)*LOG(B2N/A2N)
         KT2N    = 4.D0*PI*N2*(DG**3)*(TERM1BN - TERM2BN + TERM3BN)


      ELSE

         ! If no alkalinity, set everything to zero
         KT1  = 0.D0
         KT1N = 0.D0
         KT2  = 0.D0
         KT2N = 0.D0

      ENDIF

      ! Return to calling program
      END SUBROUTINE GET_ALK

!------------------------------------------------------------------------------

      SUBROUTINE INIT_SEASALT
!
!******************************************************************************
!  Subroutine INIT_SEASALT initializes and zeroes all module arrays 
!  (bmy, 4/26/04, 4/13/05)
!
!  NOTES:
!  (1 ) Now exit if we have allocated arrays before.  Now also allocate 
!        ALK_EMIS & N_DENS.  Now reference CMN_SIZE. (bec, bmy, 4/13/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS

      !=================================================================
      ! INIT_SEASALT begins here!
      !=================================================================

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      ! Allocate arrays
      ALLOCATE( REDGE( 0:NR_MAX, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'REDGE' )
      REDGE = 0d0

      ALLOCATE( RMID( NR_MAX, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RMID' )
      RMID = 0d0

      ALLOCATE( SRC( NR_MAX, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SRC' )
      SRC = 0d0

      ALLOCATE( SRC_N( NR_MAX, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SRC_N' )
      SRC_N = 0d0

      ALLOCATE( ALK_EMIS( IIPAR, JJPAR, LLTROP, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALK_EMIS' )
      ALK_EMIS = 0d0

      ALLOCATE( N_DENS( IIPAR, JJPAR, LLTROP, NSALT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'N_DENS' )
      N_DENS = 0d0

      ! Reset IS_INIT
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_SEASALT

!----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_SEASALT
!
!******************************************************************************
!  Subroutine INIT_SEASALT deallocates all module arrays 
!  (bmy, 4/26/04, 4/13/05)
!
!  NOTES:
!  (1 ) Now deallocates ALK_EMIS, N_DENS, SRC_N (bec, bmy, 4/13/05)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_SEASALT begins here!
      !=================================================================
      IF ( ALLOCATED( REDGE    ) ) DEALLOCATE( REDGE    )
      IF ( ALLOCATED( RMID     ) ) DEALLOCATE( RMID     )
      IF ( ALLOCATED( SRC      ) ) DEALLOCATE( SRC      )
      IF ( ALLOCATED( SRC_N    ) ) DEALLOCATE( SRC_N    )
      IF ( ALLOCATED( ALK_EMIS ) ) DEALLOCATE( ALK_EMIS )
      IF ( ALLOCATED( N_DENS   ) ) DEALLOCATE( N_DENS   )      

      ! Return to calling program
      END SUBROUTINE CLEANUP_SEASALT

!------------------------------------------------------------------------------

      ! End of module
      END MODULE SEASALT_MOD
