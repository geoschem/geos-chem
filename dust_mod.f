! $Id: dust_mod.f,v 1.17 2007/11/16 18:47:38 bmy Exp $
      MODULE DUST_MOD
!
!******************************************************************************
!  Module DUST_MOD contains routines for computing dust aerosol emissions,
!  chemistry, and optical depths. (rjp, tdf, bmy, 4/14/04, 1/23/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) FRAC_S   (REAL*8 ) : Fraction of each size classes (GINOUX only)
!  (2 ) DUSTREFF (REAL*8 ) : Dust particle radii [m]
!  (3 ) DUSTDEN  (REAL*8 ) : Soil density [kg/m3]
!  (4 ) IDDEP    (INTEGER) : Dust ID flags for drydep
!  (5 ) DRYDST1  (INTEGER) : Index for DST1 in drydep array
!  (6 ) DRYDST2  (INTEGER) : Index for DST2 in drydep array
!  (7 ) DRYDST3  (INTEGER) : Index for DST3 in drydep array
!  (8 ) DRYDST4  (INTEGER) : Index for DST4 in drydep array
!
!  Module Routines:
!  ============================================================================
!  (1 ) CHEMDUST           : Driver routine for dust chemistry
!  (2 ) DRY_SETTLING       : Routine which performs dust settling
!  (3 ) DRY_DEPOSITION     : Routine which performs dust dry deposition
!  (4 ) EMISSDUST          : Driver routine for dust emission
!  (5 ) SRC_DUST_DEAD      : Dust emissions according to DEAD   source function
!  (6 ) SRC_DUST_GINOUX    : Dust emissions according to GINOUX source function
!  (7 ) RDUST_ONLINE       : Computes dust optical depths (online dust)
!  (8 ) RDUST_OFFLINE      : Computes dust optical depths (monthly mean dust)
!  (9 ) INIT_DUST          : Allocates & initializes all module variables
!  (10) GET_DUST_ALK       : Returns dust alkalinity, and Het. Rate coeffs.
!                                                         tdf 4/15/08
!  (11) CLEANUP_DUST       : Deallocates all module variables
!
!  GEOS-CHEM modules referenced by "dust_mod.f"
!  ============================================================================
!  (1 ) dao_mod.f          : Module containing arrays for DAO met fields
!  (2 ) diag_mod.f         : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) directory_mod.f    : Module containing GEOS-CHEM data & met field dirs
!  (4 ) drydep_mod.f       : Module containing GEOS-CHEM drydep routines
!  (5 ) dust_dead_mod.f    : Module containing Zender's DEAD dust routines
!  (6 ) error_mod.f        : Module containing I/O error and NaN check routines
!  (7 ) file_mod.f         : Contains file unit numbers and error checks
!  (8 ) grid_mod.f         : Module containing horizontal grid information
!  (9 ) logical_mod.f      : Module containing GEOS-CHEM logical switches
!  (10) pressure_mod.f     : Module containing routines to compute P(I,J,L)
!  (11) time_mod.f         : Module containing routines for computing time/date
!  (12) tracer_mod.f       : Module containing GEOS-CHEM tracer array STT etc.
!  (13) tracerid_mod.f     : Module containing pointers to tracers & emissions
!  
!  NOTES:
!  (1 ) Bug fix in SRC_DUST_DEAD (bmy, 4/14/04)
!  (2 ) Now references "logical_mod.f", "directory_mod.f", and "tracer_mod.f"
!        Added comments. (bmy, 7/2/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) Bug fix in snow height computation (bmy, 11/18/05)
!  (6 ) Now only do drydep if LDRYD=T (bmy, 5/23/06)
!  (7 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (8 ) Updated output print statement in SRC_DUST_DEAD (bmy, 1/23/07)
!  (9 ) Modifications for GEOS-5 (bmy, 1/24/07)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "dust_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE                :: DRYDST1, DRYDST2,  DRYDST3, DRYDST4 
!tdf
      PRIVATE                :: DRYDAL1, DRYDAL2,  DRYDAL3, DRYDAL4
!tdf
      PRIVATE                :: DRYSO4D1, DRYSO4D2,  DRYSO4D3, DRYSO4D4
!tdf following variables now used also in sulfate_mod;   tdf 4/30/08
!tdf      PRIVATE                :: DUSTDEN, DUSTREFF, FRAC_S,  IDDEP
      PRIVATE                :: FRAC_S,  IDDEP

      ! PRIVATE module routines
      PRIVATE                :: DRY_SETTLING
      PRIVATE                :: DRY_DEPOSITION
      PRIVATE                :: SRC_DUST_DEAD
      PRIVATE                :: SRC_DUST_GINOUX

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      INTEGER                :: DRYDST1, DRYDST2, DRYDST3, DRYDST4
!tdf
      INTEGER                :: DRYDAL1, DRYDAL2,  DRYDAL3, DRYDAL4
      INTEGER                :: DRYSO4D1, DRYSO4D2,  DRYSO4D3, DRYSO4D4
      INTEGER, ALLOCATABLE   :: IDDEP(:)
!tdf
      INTEGER, ALLOCATABLE   :: JDDEP(:)
      REAL*8,  ALLOCATABLE   :: FRAC_S(:)
      REAL*8,  ALLOCATABLE   :: DUSTREFF(:)
      REAL*8,  ALLOCATABLE   :: DUSTDEN(:)
      
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!-----------------------------------------------------------------------

      SUBROUTINE CHEMDUST
!
!******************************************************************************
!  Subroutine CHEMDUST is the interface between the GEOS-CHEM main program and
!  the dust chemistry routines that mostly calculates dust dry deposition.
!  (tdf, bmy, 3/30/04, 5/23/06)
!
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" and LDUST from "logical_mod.f"
!        (bmy, 7/20/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now only do dry deposition if LDRYD = T (bmy, 5/23/06)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LDRYD,   LDUST
      USE DRYDEP_MOD,   ONLY : DEPNAME, NUMDEP
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTDST1, IDTDST2, IDTDST3, IDTDST4
!tdf
      ! Include dust alkalinity tracer IDs               tdf 04/10/08
      USE TRACERID_MOD, ONLY : IDTDAL1, IDTDAL2, IDTDAL3, IDTDAL4
!tdf
      ! Include dust sulfate tracer IDs                  tdf 06/16/08
      USE TRACERID_MOD, ONLY : IDTSO4D1, IDTSO4D2, IDTSO4D3, IDTSO4D4
!tdf
      ! Include dust nitrate tracer IDs                  tdf 06/16/08
      USE TRACERID_MOD, ONLY : IDTNITD1, IDTNITD2, IDTNITD3, IDTNITD4
!tdf
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE LOGICAL_MOD, ONLY : LPRT
!tdf
#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: N
!tdf
      INTEGER               :: IDISP

      !=================================================================
      ! CHEMDUST begins here!
      !=================================================================

      ! Execute on first call only
      IF ( FIRST ) THEN
 
         ! Stop w/ error if dust tracer flags are undefined
         IF ( IDTDST1 + IDTDST2 + IDTDST3 + IDTDST4 == 0 ) THEN
            IF ( LDUST ) THEN 
               CALL ERROR_STOP( 
     &              'LDUST=T but dust tracers are undefined!',
     &              'EMISSDUST ("dust_mod.f")' )
            ENDIF
         ENDIF

         ! Allocate arrays (if necessary)
         CALL INIT_DUST

         ! Find drydep species in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE ( 'DST1' )
                  DRYDST1 = N
               CASE ( 'DST2' )
                  DRYDST2 = N
               CASE ( 'DST3' )
                  DRYDST3 = N
               CASE ( 'DST4' )
                  DRYDST4 = N
!tdf
         ! Include dust alkalinity tracers                tdf 04/10/08
               CASE ( 'DSTAL1' )
                  DRYDAL1 = N
               CASE ( 'DSTAL2' )
                  DRYDAL2 = N
               CASE ( 'DSTAL3' )
                  DRYDAL3 = N
               CASE ( 'DSTAL4' )
                  DRYDAL4 = N
!tdf
         ! Include dust sulfate tracers                  tdf 06/16/08
               CASE ( 'SO4D1' )
                  DRYSO4D1 = N
               CASE ( 'SO4D2' )
                  DRYSO4D2 = N
               CASE ( 'SO4D3' )
                  DRYSO4D3 = N
               CASE ( 'SO4D4' )
                  DRYSO4D4 = N

               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO

         ! This may lead to out of bounds errors
         IDDEP(1) = DRYDST1
         IDDEP(2) = DRYDST2
         IDDEP(3) = DRYDST3
         IDDEP(4) = DRYDST4
 
         ! Include dust alkalinity tracers                tdf 04/10/08
         IDDEP(5) = DRYDAL1
         IDDEP(6) = DRYDAL2
         IDDEP(7) = DRYDAL3
         IDDEP(8) = DRYDAL4

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Do dust settling & deposition
      !=================================================================

      ! Dust settling
      IDISP = 0
      CALL DRY_SETTLING(   STT(:,:,:,IDTDST1:IDTDST4), 
     +                     IDTDST1, IDISP )
!tdf
      ! Include dust alkalinity settling                 tdf 04/10/08
      IDISP = 4
      CALL DRY_SETTLING(   STT(:,:,:,IDTDAL1:IDTDAL4), 
     +                     IDTDAL1, IDISP )

      ! Note: DRY SETTLING of DSTNIT, DSTSO4 occurs in GRAV_SETTLING
      !       called from CHEMSULFATE in SULFATE_MOD.    tdf 04/10/08

      ! Dust deposition
!tdf
      IF ( LDRYD ) THEN
!tdf
      IDISP = 0
      CALL DRY_DEPOSITION( STT(:,:,:,IDTDST1:IDTDST4),
     +                     IDTDST1, IDISP )
!tdf
      IDISP = 4
      CALL DRY_DEPOSITION( STT(:,:,:,IDTDAL1:IDTDAL4),
     +                     IDTDAL1, IDISP )

      ! Note: DRY DEPOSITION of DSTNIT, DSTSO4 occurs in CHEM_NIT,
      !       and CHEM_SO4, in SULFATE_MOD.              tdf 04/10/08

      ENDIF
!tdf
      IF ( LPRT ) CALL DEBUG_MSG ( '### CHEMDUST: COMPLETED')
      ! Return to calling program
      END SUBROUTINE CHEMDUST

!------------------------------------------------------------------------------

      SUBROUTINE DRY_SETTLING( TC, INDEX, IDISP )
!
!******************************************************************************
!  Subroutine DRY_SETTLING computes the dry settling of dust tracers.
!  (tdf, bmy, 3/30/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Dust tracer array 
!tdf
!  (2 ) INDEX (INTEGER) : Tracer index for XNUMOL
!  (3 ) IDISP : displacement index for IIDEP
!
!  NOTES
!  (1 ) Updated comments, cosmetic changes (bmy, 3/30/04)
!  (2 ) Remove reference to CMN, it's not needed (bmy, 7/20/04)
!  (3 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
! 
      USE DAO_MOD,      ONLY : T, BXHEIGHT
      USE DIAG_MOD,     ONLY : AD44
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
! tdf
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE LOGICAL_MOD, ONLY : LPRT
! tdf
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_GCTM"     ! g0
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      REAL*8, INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

      ! Local variables
      INTEGER               :: I, J, L, N
!tdf
      INTEGER               :: INDEX,    IDISP
      REAL*8                :: DT_SETTL, DELZ,  DELZ1
      REAL*8                :: REFF,     DEN,   CONST   
      REAL*8                :: NUM,      LAMDA, FLUX
      REAL*8                :: AREA_CM2, TC0(LLPAR)
      REAL*8                :: TOT1,     TOT2

      ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa      
      REAL*8                :: P 

      ! Diameter of aerosol [um]
      REAL*8                :: Dp

      ! Pressure * DP
      REAL*8                :: PDp 

      ! Temperature (K)    
      REAL*8                :: TEMP        

      ! Slip correction factor
      REAL*8                :: Slip        

      ! Viscosity of air (Pa s)
      REAL*8                :: Visc   

      ! Settling velocity of particle (m/s)
      REAL*8                :: VTS(LLPAR)  
      
      ! Parameters
      REAL*8,  PARAMETER    :: C1 =  0.7674D0
      REAL*8,  PARAMETER    :: C2 =  3.079d0 
      REAL*8,  PARAMETER    :: C3 =  2.573D-11
      REAL*8,  PARAMETER    :: C4 = -1.424d0

      !=================================================================
      ! DRY_SETTLING begins here!
      !=================================================================

      ! Dust settling timestep [s]
      DT_SETTL = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,        L,    N,     DEN,  REFF, DP    )
!$OMP+PRIVATE( CONST, AREA_CM2, VTS,  TEMP,  P,    PDP,  SLIP  )
!$OMP+PRIVATE( VISC,  TC0,      DELZ, DELZ1, TOT1, TOT2, FLUX  )

      ! Loop over dust bins
      DO N = 1, NDSTBIN

         ! Initialize
         DEN   = DUSTDEN(N)
         REFF  = DUSTREFF(N)
         DP    = 2D0 * REFF * 1.D6              ! Dp [um] = particle diameter
         CONST = 2D0 * DEN * REFF**2 * G0 / 9D0
         print *, 'Inside DRY_SETTLING, N, DEN, REFF ',N,DEN, REFF

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2(J)

            ! Loop over longitudes
            DO I = 1, IIPAR
            
               ! Initialize settling velocity
               DO L = 1, LLPAR
                  VTS(L) = 0d0
               ENDDO

               ! Loop over levels
               DO L = 1, LLPAR

                  ! Get P [kPa], T [K], and P*DP
                  P    = GET_PCENTER(I,J,L) * 0.1d0
                  TEMP = T(I,J,L)
                  PDP  = P * DP

                  !=====================================================
                  ! # air molecule number density
                  ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
                  !
                  ! # gas mean free path
                  ! lamda = 1.d6 / 
                  !     &   ( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
                  !
                  ! # Slip correction
                  ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * 
                  !      &  exp( -1.1 * Dp / (2. * lamda))) / Dp
                  !=====================================================
                  ! NOTE, Slip correction factor calculations following 
                  !       Seinfeld, pp464 which is thought to be more 
                  !       accurate but more computation required.
                  !=====================================================

                  ! Slip correction factor as function of (P*dp)
                  SLIP = 1d0 + 
     &                   ( 15.60d0 + 7.0d0 * EXP(-0.059d0*PDP) ) / PDP
            
                  !=====================================================
                  ! NOTE, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
                  ! which produce slip correction factor with small 
                  ! error compared to the above with less computation.
                  !=====================================================

                  ! Viscosity [Pa s] of air as a function of temp (K)
                  VISC = 1.458d-6 * (TEMP)**(1.5d0) / ( TEMP + 110.4d0 )

                  ! Settling velocity [m/s]
                  VTS(L) = CONST * SLIP / VISC

               ENDDO

               ! Method is to solve bidiagonal matrix 
               ! which is implicit and first order accurate in Z
               DO L = 1, LLPAR
                  TC0(L) = TC(I,J,L,N)
               ENDDO

               ! We know the boundary condition at the model top
               L           = LLTROP
               DELZ        = BXHEIGHT(I,J,L)
               TC(I,J,L,N) = TC(I,J,L,N) / 
     &                       ( 1.d0 + DT_SETTL * VTS(L) / DELZ )

               DO L = LLTROP-1, 1, -1
                  DELZ        = BXHEIGHT(I,J,L)
                  DELZ1       = BXHEIGHT(I,J,L+1)
                  TC(I,J,L,N) = 1.d0 / 
     &                          ( 1.d0 + DT_SETTL * VTS(L)   / DELZ )
     &                 * (TC(I,J,L,N)  + DT_SETTL * VTS(L+1) / DELZ1
     &                 *  TC(I,J,L+1,N) )
               ENDDO

               !========================================================      
               ! ND44: Dry deposition diagnostic [#/cm2/s]
               !========================================================
               IF ( ND44 > 0 ) THEN

                  ! Initialize
                  TOT1 = 0d0
                  TOT2 = 0d0
            
                  ! Compute column totals of TCO(:) and TC(I,J,:,N)
                  DO L = 1, LLPAR
                     TOT1 = TOT1 + TC0(L)
                     TOT2 = TOT2 + TC(I,J,L,N)
                  ENDDO

                  ! Convert dust flux from [kg/s] to [#/cm2/s]
                  FLUX = ( TOT1 - TOT2 ) / DT_SETTL  
                  FLUX = FLUX * XNUMOL(INDEX) / AREA_CM2 

                  ! Save in AD44
!tdf                  AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
      AD44(I,J,IDDEP(N+IDISP),1) = AD44(I,J,IDDEP(N+IDISP),1) + FLUX
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
!tdf
      IF ( LPRT ) CALL DEBUG_MSG ( '### DRY_SETTLING: COMPLETED')

      ! Return to calling program
      END SUBROUTINE DRY_SETTLING

!------------------------------------------------------------------------------

      SUBROUTINE DRY_DEPOSITION ( TC, INDEX, IDISP )
!
!******************************************************************************
!  Subroutine DRY_DEPOSITION computes the loss of dust due to dry deposition
!  at the surface using an implicit method. (tdf, bmy, 3/30/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Dust tracer array   
!tdf
!  (2 ) INDEX (INTEGER) : Tracer index for XNUMOL
!  (3 ) IDISP : displacement index for IIDEP
!
!  NOTES: 
!  (1 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV 
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
! tdf
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE LOGICAL_MOD, ONLY : LPRT
! tdf
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      REAL*8, INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)
      INTEGER               :: INDEX
      INTEGER               :: IDISP
!tdf
      ! local variables
      INTEGER               :: I,   J,   L,      N
      REAL*8                :: OLD, NEW, DTCHEM, FLUX, AREA_CM2
!tdf
      REAL*8                :: ACID

      !=================================================================
      ! DRY_DEPOSITION begins here!
      !=================================================================

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over dust bins
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, OLD, NEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over dust bins
      DO N = 1, NDSTBIN    

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2(J)

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Original dust concentration at surface
               OLD = TC(I,J,1,N)
!tdf
               ! Dust left after dry deposition
!tdf : Deliberately use same freq. (DEPSAV) for dust alkalinity as for dust
               NEW = OLD * EXP( -DEPSAV(I,J,IDDEP(N)) * DTCHEM  )
!tdf
               !========================================================
               ! ND44 diagnostic: dust drydep loss [#/cm2/s]
               !========================================================
               IF ( ND44 > 0 ) THEN

                  ! Convert drydep flux from [kg/s] to [#/cm2/s]
                  FLUX     = ( OLD - NEW ) / DTCHEM 
                  FLUX     = FLUX * XNUMOL(INDEX) / AREA_CM2 
               
                  ! Store in AD44
!tdf           AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
      AD44(I,J,IDDEP(N+IDISP),1) = AD44(I,J,IDDEP(N+IDISP),1) + FLUX
               ENDIF

               ! Save back into STT
               TC(I,J,1,N) = NEW 
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
ctdf
      IF ( LPRT ) CALL DEBUG_MSG ( '### DRY_DEPOSITION: COMPLETED')

      ! Return to calling program
      END SUBROUTINE DRY_DEPOSITION

!------------------------------------------------------------------------------

      SUBROUTINE EMISSDUST
!
!******************************************************************************
!  Subroutine EMISSDUST is the driver routine for the dust emission
!  module.  You may call either the GINOUX or the DEAD dust source 
!  function. (tdf, bmy, 3/30/04, 10/3/05)
!
!  NOTES:
!  (1 ) Now reference LDEAD, LDUST, LPRT from "logical_mod.f".  Now reference!
!        STT from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F(0 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP, DEBUG_MSG
      USE LOGICAL_MOD,  ONLY : LDEAD, LDUST, LPRT
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTDST1, IDTDST2, IDTDST3, IDTDST4
!tdf
      ! Include dust alkalinity                         tdf 04/10/08
      USE TRACERID_MOD, ONLY : IDTDAL1, IDTDAL2, IDTDAL3, IDTDAL4

#     include "CMN_SIZE"     ! Size parameters
      
      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      
      !=================================================================
      ! EMISSDUST begins here!
      !=================================================================

      ! Execute on first-call only
      IF ( FIRST ) THEN

         ! Return if dust ID flags are not defined
         IF ( IDTDST1 + IDTDST2 + IDTDST3 + IDTDST4 == 0 ) THEN
            IF ( LDUST ) THEN 
               CALL ERROR_STOP( 
     &              'LDUST=T but dust tracers are undefined!',
     &              'EMISSDUST ("dust_mod.f")' )
            ENDIF
         ENDIF
          
         ! Allocate module arrays
         CALL INIT_DUST
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSDUST: a INIT_DUST' )

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Call appropriate emissions routine
      !=================================================================
      IF ( LDEAD ) THEN

         ! Use Zender's DEAD dust source function
         ! Include DUST Alkalinity SOURCE          tdf 04/10/08
!tdf         CALL SRC_DUST_DEAD( STT(:,:,:,IDTDST1:IDTDST4) )
         CALL SRC_DUST_DEAD( STT(:,:,:,IDTDST1:IDTDST4),
     &                       STT(:,:,:,IDTDAL1:IDTDAL4) )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_DEAD' )
         
      ELSE

         ! Use Paul Ginoux's dust source function
!tdf Include DUST Alkalinity SOURCE              tdf 04/10/08
!tdf         CALL SRC_DUST_GINOUX( STT(:,:,:,IDTDST1:IDTDST4) )
         CALL SRC_DUST_GINOUX( STT(:,:,:,IDTDST1:IDTDST4),
     &                         STT(:,:,:,IDTDAL1:IDTDAL4) )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_GINOUX')

      ENDIF     
      ! Return to calling program
      END SUBROUTINE EMISSDUST

!------------------------------------------------------------------------------

!tdf      SUBROUTINE SRC_DUST_DEAD( TC )
      SUBROUTINE SRC_DUST_DEAD( TC, TALK )
!
!******************************************************************************
!  Subroutine SRC_DUST_DEAD is the DEAD model dust emission scheme, 
!  alternative to Ginoux scheme.  Increments the TC array with emissions 
!  from the DEAD model.  (tdf, bmy, 4/8/04, 1/23/07)
!
!  Input:
!         SRCE_FUNK Source function                               (-)
!         for 1: Sand, 2: Silt, 3: Clay
!         DUSTDEN   Dust density                                  (kg/m3)
!         DUSTREFF  Effective radius                              (um)
!         AD        Air mass for each grid box                    (kg)
!         NTDT      Time step                                     (s)
!         W10M      Velocity at the anemometer level (10meters)   (m/s)
!         GWET      Surface wetness                               (-)
!
!  Parameters used in GEOS-CHEM
!
!  Longitude: IIPAR 
!  Latitude : JJPAR 
!  Levels   : LLPAR  =  20 (GEOS-1), 26 (GEOS-strat), 30 (GEOS-terra)
!  Size bins: NDSTBIN  =   4
!
!  Dust properties used in GOCART
!
!  Size classes: 01-1, 1-1.8, 1.8-3, 3-6 (um)
!  Radius: 0.7, 1.5, 2.5, 4  (um)
!  Density: 2500, 2650, 2650, 2650 (kg/m3)
!
!  NOTES:
!  (1 ) Added OpenMP parallelization, added comments (bmy, 4/8/04)
!  (2 ) Bug fix: DSRC needs to be held PRIVATE (bmy, 4/14/04)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Bug fix: It should be SNOW/1d3 not SNOW*1d3 (tdf, bmy, 11/18/05)
!  (6 ) Updated output statement (bmy, 1/23/07)
!  (7 ) Use SNOMAS (m H2O) for GEOS-5 (bmy, 1/24/07)
!  (8 ) Now included dust alkalinity source            (tdf 04/10/08)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : BXHEIGHT,     GWETTOP,   LWI
      USE DAO_MOD,       ONLY : SNOW,         SPHU,      T    
      USE DAO_MOD,       ONLY : TS,           UWND,      VWND
      USE DAO_MOD,       ONLY : SNOMAS
      USE DUST_DEAD_MOD, ONLY : GET_TIME_INVARIANT_DATA, GET_ORO
      USE DUST_DEAD_MOD, ONLY : GET_MONTHLY_DATA,        DST_MBL
      USE DIAG_MOD,      ONLY : AD06
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE GRID_MOD,      ONLY : GET_YMID_R
      USE PRESSURE_MOD,  ONLY : GET_PEDGE,       GET_PCENTER 
      USE TIME_MOD,      ONLY : GET_TS_EMIS,     GET_MONTH   
      USE TIME_MOD,      ONLY : GET_DAY_OF_YEAR, ITS_A_NEW_MONTH
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
!tdf acidic uptake on dust
      USE LOGICAL_MOD,   ONLY : LDSTUP

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DIAG"      ! ND06
#     include "CMN_GCTM"      ! g0

      !----------------
      ! Arguments
      !----------------
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

!tdf Include dust Alkalinity                            tdf 4/10/2K8
      REAL*8,  INTENT(INOUT) :: TALK(IIPAR,JJPAR,LLPAR,NDSTBIN)

      !-----------------
      ! Local variables
      !-----------------

      ! Scalars
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I,      J,      L,       N
      INTEGER                :: M,      IOS,    INC,     LAT_IDX
      INTEGER                :: NDB,    NSTEP
      REAL*8                 :: W10M,   DEN,    DIAM,    U_TS0 
      REAL*8                 :: U_TS,   SRCE_P, Reynol,  YMID_R
      REAL*8                 :: ALPHA,  BETA,   GAMMA,   CW
      REAL*8                 :: DTSRCE, XTAU,   P1,      P2
      REAL*8                 :: DOY             
      CHARACTER(LEN=255)     :: FILENAME

      ! Arrays
      INTEGER                :: OROGRAPHY(IIPAR,JJPAR)
      REAL*8                 :: PSLON(IIPAR)         ! surface pressure
      REAL*8                 :: PTHICK(IIPAR)        ! delta P (L=1)
      REAL*8                 :: PMID(IIPAR)          ! mid layer P (L=1)
      REAL*8                 :: TLON(IIPAR)          ! temperature (L=1)
      REAL*8                 :: THLON(IIPAR)         ! pot. temp. (L=1)
      REAL*8                 :: ULON(IIPAR)          ! U-wind (L=1)
      REAL*8                 :: VLON(IIPAR)          ! V-wind (L=1)
      REAL*8                 :: BHT2(IIPAR)          ! half box height (L=1)
      REAL*8                 :: Q_H2O(IIPAR)         ! specific humidity (L=1)
      REAL*8                 :: ORO(IIPAR)           ! "orography" 
      REAL*8                 :: SNW_HGT_LQD(IIPAR)   ! equivalent snow ht.
      REAL*8                 :: DSRC(IIPAR,NDSTBIN)  ! dust mixing ratio incr.

      !----------------
      ! Parameters
      !----------------
      REAL*8, PARAMETER      :: Ch_dust = 9.375d-10
      REAL*8, PARAMETER      :: G       = g0 * 1.D2
      REAL*8, PARAMETER      :: RHOA    = 1.25D-3
      REAL*8, PARAMETER      :: CP      = 1004.16d0
      REAL*8, PARAMETER      :: RGAS    = 8314.3d0 / 28.97d0
      REAL*8, PARAMETER      :: AKAP    = RGAS     / CP
      REAL*8, PARAMETER      :: P1000   = 1000d0

      ! External functions
      REAL*8,  EXTERNAL       :: SFCWINDSQR

      !=================================================================
      ! SRC_DUST_DEAD begins here!
      !=================================================================      

      ! DTSRCE is the emission timestep in seconds
      DTSRCE = GET_TS_EMIS() * 60d0

      ! DOY is the day of year (0-365 or 0-366)
      DOY    = DBLE( GET_DAY_OF_YEAR() )

      !=================================================================
      ! Read data fields for the DEAD model from disk
      !=================================================================      
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   )
         WRITE( 6, 110   )
         WRITE( 6, 120   )
         WRITE( 6, 130   ) 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       
         ! FORMAT strings
 100     FORMAT( 'D E A D   D U S T   M O B I L I Z A T I O N'         )
 110     FORMAT( 'Routines from DEAD model by Charlie Zender et al'    )
 120     FORMAT( 'Modified for GEOS-CHEM by D. Fairlie and R. Yantosca')
 130     FORMAT( 'Last Modification Date: 1/23/07'                     )

         ! Read fields for DEAD that are time-invariant
         CALL GET_TIME_INVARIANT_DATA
         
         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Read monthly data for DEAD
      IF ( ITS_A_NEW_MONTH() ) THEN
         CALL GET_MONTHLY_DATA
      ENDIF

      !=================================================================
      ! Call dust mobilization scheme
      !=================================================================

      ! Make OROGRAPHY array from GEOS-CHEM LWI
      CALL GET_ORO( OROGRAPHY )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,      P1,   P2,   PTHICK, PMID, TLON        )
!$OMP+PRIVATE( THLON, ULON,   VLON, BHT2, Q_H2O,  ORO,  SNW_HGT_LQD )
!$OMP+PRIVATE( N,     YMID_R, DSRC                                  )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Pressure [hPa] at bottom and top edge of level 1
            P1             = GET_PEDGE(I,J,1)
            P2             = GET_PEDGE(I,J,2)

            ! Pressure thickness of 1st layer [Pa]
            PTHICK(I)      = ( P1 - P2 ) * 100d0

            ! Pressure at midpt of surface layer [Pa]
            PMID(I)        = GET_PCENTER(I,J,1) * 100d0 

            ! Temperature [K] at surface layer
            TLON(I)        = T(I,J,1)

            ! Potential temperature [K] at surface layer
            THLON(I)       = TLON(I) * ( P1000 / PMID(I) )**AKAP

            ! U and V winds at surface [m/s]
            ULON(I)        = UWND(I,J,1)
            VLON(I)        = VWND(I,J,1)
            
            ! Half box height at surface [m]
            BHT2(I)        = BXHEIGHT(I,J,1) / 2.d0

            ! Specific humidity at surface [kg H2O/kg air]
            Q_H2O(I)       = SPHU(I,J,1) / 1000.d0         

            ! Orography at surface
            ! Ocean is 0; land is 1; ice is 2
            ORO(I)         = OROGRAPHY(I,J) 

            ! Snow height [m H2O]
#if   defined( GEOS_5  )
            SNW_HGT_LQD(I) = SNOMAS(I,J)
#else
            SNW_HGT_LQD(I) = SNOW(I,J) / 1000d0
#endif
            ! Dust tracer and increments
            DO N = 1, NDSTBIN
               DSRC(I,N)   = 0.0d0
            ENDDO
         ENDDO 

         !==============================================================
         ! Call dust mobilization driver (DST_MBL) for latitude J
         !==============================================================
         
         ! Latitude in RADIANS
         YMID_R = GET_YMID_R(J)

         ! Call DEAD dust mobilization
         CALL DST_MBL( DOY,    BHT2,  J,     YMID_R, ORO,    
     &                 PTHICK, PMID,  Q_H2O, DSRC,   SNW_HGT_LQD, 
     &                 DTSRCE, TLON,  THLON, VLON,   ULON,
     &                 FIRST,  J )

         ! Update
         DO N = 1, NDSTBIN
         DO I = 1, IIPAR

            ! Add dust emissions into tracer array [kg]
            TC(I,J,1,N) = TC(I,J,1,N) + DSRC(I,N) 
!tdf
            IF ( LDSTUP ) THEN
            ! Include DUST Alkalinity SOURCE, assuming an alkalinity
!tdf with 3% Ca, there's also 1% equ. Mg, makes 4%
            ! of 4% by weight [kg].                  !tdf 05/10/08
              TALK(I,J,1,N) = TALK(I,J,1,N) + 0.04 * DSRC(I,N)
!tdf with 1.5% Ca, there's also 0.5% equ. Mg, makes 2%
!tdf            ! of 2% by weight [kg].                  !tdf 05/10/08
!tdf              TALK(I,J,1,N) = TALK(I,J,1,N) + 0.02 * DSRC(I,N)
            ENDIF
 
            ! ND19 diagnostics [kg]
            IF ( ND06 > 0 ) THEN
               AD06(I,J,N) = AD06(I,J,N) + DSRC(I,N)
               IF ( LDSTUP ) THEN
!tdf dust alkalninity diagnostic
!tdf with 3% Ca, there's also 1% equ. Mg, makes 4%
            ! of 4% by weight [kg].                  !tdf 05/10/08
               AD06(I,J,N+NDSTBIN) = AD06(I,J,N+NDSTBIN)
     &                             + 0.04 * DSRC(I,N)
            ! of 2% by weight [kg].                  !tdf 05/10/08
!tdf     &                             + 0.02 * DSRC(I,N)
               ENDIF
            ENDIF
         ENDDO
         ENDDO

      ENDDO  
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRC_DUST_DEAD

!------------------------------------------------------------------------------

!tdf      SUBROUTINE SRC_DUST_GINOUX( TC )
      SUBROUTINE SRC_DUST_GINOUX( TC, TALK )
!
!******************************************************************************
!  Paul GINOUX dust source function 
!  (Added to GEOS-CHEM, tdf, bmy, 4/8/04, 7/20/04)
!
!  This subroutine updates the surface mixing ratio of dust aerosols for
!  NDSTBIN size bins. The uplifting of dust depends in space on the source 
!  function, and in time and space on the soil moisture and surface
!  wind speed (10 meters). Dust is uplifted if the wind speed is greater
!  than a threshold velocity which is calculated with the formula of
!  Marticorena et al.  (JGR, v.102, p 23277-23287, 1997).
!  To run this subroutine you need the source function which can be
!  obtained by contacting Paul Ginoux at ginoux@rondo.gsfc.nasa.gov
!  If you are not using GEOS DAS met fields, you will most likely need
!  to adapt the adjusting parameter.
! 
!  Contact: Paul Ginoux (ginoux@rondo.gsfc.nasa.gov) 
! 
! 
!  Input:
!    SRCE_FUNK Source function                               (-)
!              for 1: Sand, 2: Silt, 3: Clay
!
!    DUSTDEN   Dust density                                  (kg/m3)
!    DUSTREFF  Effective radius                              (um)
!    AD        Air mass for each grid box                    (kg)
!    NTDT      Time step                                     (s)
!    W10m      Velocity at the anemometer level (10meters)   (m/s)
!    GWET      Surface wetness                               (-)
!       
! 
!  Parameters used in GEOS-CHEM
! 
!  Longitude: IIPAR 
!  Latitude : JJPAR 
!  Levels   : LLPAR  =  20 (GEOS-1), 26 (GEOS-strat), 30 (GEOS-terra)
!  Size bins: NDSTBIN  =   4
! 
!  Dust properties used in GOCART
! 
!  Size classes: 01-1, 1-1.8, 1.8-3, 3-6 (um)
!  Radius: 0.7, 1.5, 2.5, 4  (um)
!  Density: 2500, 2650, 2650, 2650 (kg/m3)
! 
!  References:
!  ============================================================================
!  (1 ) Ginoux, P., M. Chin, I. Tegen, J. Prospero, B. Hoben, O. Dubovik,
!        and S.-J. Lin, "Sources and distributions of dust aerosols simulated
!        with the GOCART model", J. Geophys. Res., 2001
!  (2 ) Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
!
!  NOTES:
!  (1 ) Added OpenMP parallelization (bmy, 4/8/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now included dust alkalinity source            (tdf 04/10/08)
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DAO_MOD,       ONLY : GWETTOP   
      USE DIAG_MOD,      ONLY : AD06
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE TIME_MOD,      ONLY : GET_TS_EMIS
      USE GRID_MOD,      ONLY : GET_AREA_M2
!tdf acidic uptake on dust
      USE LOGICAL_MOD,   ONLY : LDSTUP

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! ND19, LD13 (for now)
#     include "CMN_GCTM"  ! g0

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)

!tdf Include dust Alkalinity                            tdf 4/10/2K8
      REAL*8,  INTENT(INOUT) :: TALK(IIPAR,JJPAR,LLPAR,NDSTBIN)

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, N, M, IOS
      INTEGER                :: IPOINT(NDSTBIN) = (/3, 2, 2, 2/)
      REAL*4                 :: ARRAY(IIPAR,JJPAR,3)
      REAL*8, SAVE           :: SRCE_FUNC(IIPAR,JJPAR,3)
      REAL*8                 :: W10M,   DEN,  DIAM,   U_TS0, U_TS
      REAL*8                 :: SRCE_P, DSRC, REYNOL, ALPHA, BETA
      REAL*8                 :: GAMMA,  CW,   DTSRCE, AREA_M2
      CHARACTER(LEN=255)     :: FILENAME

      ! Transfer coeff for type natural source  (kg*s2/m5)
      REAL*8, PARAMETER      :: CH_DUST  = 9.375d-10
      REAL*8, PARAMETER      :: G        = G0 * 1.d2
      REAL*8, PARAMETER      :: RHOA     = 1.25d-3

      ! External functions
      REAL*8,  EXTERNAL      :: SFCWINDSQR

      !=================================================================
      ! SRC_DUST_GINOUX begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      !=================================================================
      ! Read dust source function
      !=================================================================
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   )
         WRITE( 6, 110   )
         WRITE( 6, 120   )
         WRITE( 6, 130   ) 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       
         ! FORMAT strings
 100     FORMAT( 'G I N O U X   D U S T   M O B I L I Z A T I O N'     )
 110     FORMAT( 'Routines originally by Paul Ginoux, GSFC'            )
 120     FORMAT( 'Modified for GEOS-CHEM by D. Fairlie and R. Yantosca')
 130     FORMAT( 'Last Modification Date: 4/6/04'                      )

         ! Filename
         FILENAME = TRIM( DATA_DIR )        // 
     &              'dust_200203/NSP.dust.' // GET_RES_EXT()

         ! Open file
         OPEN( 65, FILE=FILENAME,       STATUS='OLD', 
     &             FORM= 'UNFORMATTED', IOSTAT=IOS )
         IF ( IOS > 0 ) CALL IOERROR( IOS, 65, 'SRC_DUST_GINOUX:1' )

         ! Read data
         READ( 65, IOSTAT=IOS ) ARRAY
         IF ( IOS > 0 ) CALL IOERROR( IOS, 65, 'SRC_DUST_GINOUX:2' )

         ! Close file
         CLOSE( 65 )

         ! Cast to REAL*8
         SRCE_FUNC = ARRAY

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF
        
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,      J,     M,    N,      DEN,   DIAM    )
!$OMP+PRIVATE( REYNOL, ALPHA, BETA, GAMMA,  U_TS0, AREA_M2 )
!$OMP+PRIVATE( CW,     U_TS,  W10M, SRCE_P, DSRC           )

      ! Loop over size bins
      DO N = 1, NDSTBIN

         !==============================================================
         ! Threshold velocity as a function of the dust density and the 
         ! diameter from Bagnold (1941), valid for particles larger 
         ! than 10 um.
         ! 
         ! u_ts0 = 6.5*sqrt(dustden(n)*g0*2.*dustreff(n))
         !
         ! Threshold velocity from Marticorena and Bergametti
         ! Convert units to fit dimensional parameters
         !==============================================================
         DEN    = DUSTDEN(N) * 1.d-3                 ! [g/cm3]
         DIAM   = 2d0 * DUSTREFF(N) * 1.d2           ! [cm in diameter]
         REYNOL = 1331.d0 * DIAM**(1.56d0) + 0.38d0  ! [Reynolds number]
         ALPHA  = DEN * G * DIAM / RHOA
         BETA   = 1d0 + ( 6.d-3 / ( DEN * G * DIAM**(2.5d0) ) )
         GAMMA  = ( 1.928d0 * REYNOL**(0.092d0) ) - 1.d0
         
         !==============================================================
         ! I think the 129.d-5 is to put U_TS in m/sec instead of cm/sec
         ! This is a threshold friction velocity!       from M&B
         ! i.e. Ginoux uses the Gillette and Passi formulation
         ! but has substituted Bagnold's Ut with M&B's U*t.
         ! This appears to be a problem.  (tdf, 4/2/04)
         !==============================================================

         ! [m/s] 
         U_TS0  = 129.d-5 * SQRT( ALPHA ) * SQRT( BETA ) / SQRT( GAMMA )
         M      = IPOINT(N)

         ! Loop over latitudes
         DO J = 1, JJPAR
            
            ! Get grid box surface area [m2]
            AREA_M2 = GET_AREA_M2(J)

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Fraction of emerged surfaces 
               ! (subtract lakes, coastal ocean,...)
               CW = 1.d0

               ! Case of surface dry enough to erode
               IF ( GWETTOP(I,J) < 0.2d0 ) THEN
               
                  U_TS = U_TS0*( 1.2d0 + 
     &                           0.2d0*LOG10( MAX(1.d-3,GWETTOP(I,J))))
                  U_TS = MAX( 0.d0, U_TS )

               ELSE

                  ! Case of wet surface, no erosion
                  U_TS = 100.d0

               ENDIF

               ! 10m wind speed [m/s]
               W10M   = SQRT( SFCWINDSQR(I,J) )

               ! Units are m2
               SRCE_P = FRAC_S(N) * SRCE_FUNC(I,J,M) * AREA_M2

               ! Dust source increment [kg]
               DSRC   = CW     * CH_DUST * SRCE_P * W10M**2
     &                * ( W10M - U_TS )  * DTSRCE  

               ! Not less than zero
               IF ( DSRC < 0.d0 ) DSRC = 0.d0

               ! Dust SOURCE at first model level [kg].
               TC(I,J,1,N) = TC(I,J,1,N) + DSRC 
!tdf
               IF ( LDSTUP ) THEN
               ! Include DUST Alkalinity SOURCE, assuming an alkalinity
!tdf               ! of 6.8% by weight [kg].                  !tdf 04/10/08
!tdf               TALK(I,J,1,N) = TALK(I,J,1,N) + 0.068 * DSRC
!tdf 3% Ca + equ 1% Mg = 4% alkalinity
               ! of 4% by weight [kg].                  !tdf 05/10/08
                 TALK(I,J,1,N) = TALK(I,J,1,N) + 0.04 * DSRC
!tdf 1.5% Ca + equ 0.5% Mg = 2% alkalinity
!tdf               ! of 2% by weight [kg].                  !tdf 05/10/08
!tdf               TALK(I,J,1,N) = TALK(I,J,1,N) + 0.02 * DSRC
               ENDIF

               !========================================================
               ! ND06 diagnostics: dust emissions [kg/timestep]
               !========================================================
               IF ( ND06 > 0 ) THEN
                  AD06(I,J,N) = AD06(I,J,N) + DSRC
!tdf dust alkalninity diagnostic
                  IF ( LDSTUP ) THEN
                    AD06(I,J,N+NDSTBIN) = AD06(I,J,N+NDSTBIN)
     &                                  + 0.04 * DSRC
!tdf     &                                  + 0.02 * DSRC
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRC_DUST_GINOUX

!------------------------------------------------------------------------------

      SUBROUTINE RDUST_ONLINE( DUST )
!
!******************************************************************************
!  Subroutine RDUST reads global mineral dust concentrations as determined 
!  by P. Ginoux.  Calculates dust optical depth at each level for the
!  FAST-J routine "set_prof.f". (rvm, rjp, tdf, bmy, 4/1/04, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DUST (REAL*8) : Dust from soils [kg/m3]
!
!  NOTES:
!  (1 ) Bundled into "dust_mod.f" (bmy, 4/1/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f".  Now parallelize over
!        the L-dimension for ND21 diagnostics. (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,    ONLY : ERADIUS, IXSAVE, IYSAVE 
      USE COMODE_MOD,    ONLY : IZSAVE,  JLOP,   TAREA
      USE DAO_MOD,       ONLY : BXHEIGHT
      USE DIAG_MOD,      ONLY : AD21
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      IMPLICIT NONE

#     include "cmn_fj.h"   ! LPAR, CMN_SIZE
#     include "jv_cmn.h"   ! ODMDUST, QAA, RAA
#     include "CMN_DIAG"   ! ND21, LD21
#     include "comode.h"   ! NTTLOOP

      ! Arguments
      REAL*8, INTENT(IN) :: DUST(IIPAR,JJPAR,LLPAR,NDUST)

      ! Local variables
      INTEGER            :: I, J, JLOOP, L, N
      REAL*8             :: MSDENS(NDUST), XTAU

      !=================================================================
      ! RDUST_ONLINE begins here!
      !=================================================================

      ! Dust density 
      MSDENS(1) = 2500.0d0
      MSDENS(2) = 2500.0d0
      MSDENS(3) = 2500.0d0
      MSDENS(4) = 2500.0d0
      MSDENS(5) = 2650.0d0
      MSDENS(6) = 2650.0d0
      MSDENS(7) = 2650.0d0

      !=================================================================     
      ! Convert concentration [kg/m3] to optical depth [unitless].
      !
      ! ODMDUST = ( 0.75 * BXHEIGHT * CONC * QAA ) / 
      !           ( MSDENS * RAA * 1e-6 )
      ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
      !
      !  Units ==> DUST     [ kg/m3    ]
      !            MSDENS   [ kg/m3    ]
      !            RAA      [ um       ]
      !            BXHEIGHT [ m        ]
      !            QAA      [ unitless ]
      !            ODMDUST  [ unitless ]
      !
      ! NOTES: 
      ! (1) Do the calculation at QAA(4,:) (i.e. 999 nm).          
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, NDUST
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ODMDUST(I,J,L,N) = 0.75d0        * BXHEIGHT(I,J,L) * 
     &                      DUST(I,J,L,N) * QAA(4,14+N)     / 
     &                     ( MSDENS(N) * RAA(4,14+N) * 1.0D-6 )
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !==============================================================
      ! Calculate Dust Surface Area
      !
      ! Units ==> DUST     [ kg dust/m^3 air    ]
      !           MSDENS   [ kg dust/m^3 dust   ]
      !           RAA      [ um                 ]
      !           TAREA    [ cm^2 dust/cm^3 air ]
      !           ERADIUS  [ cm                 ]
      !
      ! NOTE: first find volume of dust (cm3 dust/cm3 air), then 
      !       multiply by 3/radius to convert to surface area in cm2
      !  
      ! TAREA(:,1:NDUST) and ERADIUS(:,1:NDUST) are for 
      ! the NDUST FAST-J dust wavelength bins (read into DUST)
      !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L, N )
      DO N     = 1, NDUST
      DO JLOOP = 1, NTTLOOP

         ! Compute 3-D grid box indices
         I = IXSAVE(JLOOP)
         J = IYSAVE(JLOOP)
         L = IZSAVE(JLOOP)

         ERADIUS(JLOOP,N) = RAA(4,14+N) * 1.0D-4

         TAREA(JLOOP,N)   = 3.D0 / ERADIUS(JLOOP,N) *
     &                      DUST(I,J,L,N) / MSDENS(N)  
! tdf
         if (I .eq. 1 .and. J .eq. 63 .and. L .eq. 6) then
            print *,' RDUST_ONLINE: TAREA'
            print *,' IBIN,    TAREA'
            write (6,51) N, TAREA(JLOOP,N)
   51 format (' ',I4,E12.3,'  cm2/cm3')
            endif
! tdf
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! ND21 Diagnostic: 
      !
      ! Tracer #1: Cloud optical depths    (from "optdepth_mod.f")
      ! Tracer #2: Max Overlap Cld Frac    (from "optdepth_mod.f")
      ! Tracer #3: Random Overlap Cld Frac (from "optdepth_mod.f")
      ! Tracer #4: Dust optical depths at 400 nm (from all size bins)
      ! Tracer #5: Dust surface areas (from all size bins)
      !==============================================================
      IF ( ND21 > 0 ) THEN

         DO N = 1, NDUST

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L )
            DO L = 1, LD21
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               !--------------------------------------
               ! ND21 tracer #4: Dust optical depths
               !--------------------------------------
               AD21(I,J,L,4) = AD21(I,J,L,4) + 
     &              ( ODMDUST(I,J,L,N) * QAA(2,14+N) / QAA(4,14+N) )

               !--------------------------------------
               ! ND21 tracer #5: Dust surface areas
               !--------------------------------------
               IF ( L <= LLTROP ) THEN
               
                  ! Convert 3-D indices to 1-D index
                  ! JLOP is only defined in the tropopause
                  JLOOP = JLOP(I,J,L)
             
                  ! Add to AD21
                  IF ( JLOOP > 0 ) THEN
                     AD21(I,J,L,5) = AD21(I,J,L,5) + TAREA(JLOOP,N)
                  ENDIF
               ENDIF
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDDO
      ENDIF 

      ! Return to calling program
      END SUBROUTINE RDUST_ONLINE

!------------------------------------------------------------------------------

      SUBROUTINE RDUST_OFFLINE( THISMONTH, THISYEAR )
!
!******************************************************************************
!  Subroutine RDUST_OFFLINE reads global mineral dust concentrations as 
!  determined by P. Ginoux.  Calculates dust optical depth at each level for 
!  the FAST-J routine "set_prof.f". (rvm, bmy, 9/30/00, 8/4/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Number of the current month (1-12)
!  (2 ) THISYEAR  (INTEGER) : 4-digit year number (e.g. 1996, 2001)
!
!  NOTES:
!  (1 ) RDUST was patterned after rdaerosol.f (rvm, 9/30/00)
!  (2 ) Don't worry about rewinding the binary file...reading from
!        binary files is pretty fast.  And it's only done once a month.
!  (3 ) Now references punch file utility routines from F90 module
!        "bpch2_mod.f".  Also reference variable DATA_DIR from the
!         header file "CMN_SETUP". (bmy, 9/30/00) 
!  (4 ) Now selects proper GEOS-STRAT dust field for 1996 or 1997.
!        Also need to pass THISYEAR thru the arg list. (rvm, bmy, 11/21/00)
!  (5 ) CONC is now declared as REAL*8 (rvm, bmy, 12/15/00)
!  (6 ) Removed obsolete code from 12/15/00 (bmy, 12/21/00)
!  (7 ) CONC(IGLOB,JGLOB,LGLOB,NDUST) is now CONC(IIPAR,JJPAR,LLPAR,NDUST).
!        Now use routine TRANSFER_3D from "transfer_mod.f" to cast from REAL*4
!        to REAL*8 and also to convert from {IJL}GLOB to IIPAR,JJPAR,LLPAR 
!        space.  Use 3 arguments in call to GET_TAU0.  Updated comments.
!        (bmy, 9/26/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (9 ) Now reference ERADIUS, IXSAVE, IYSAVE, IZSAVE, TAREA from 
!        "comode_mod.f".  Compute ERADIUS and TAREA for the NDUST dust
!        size bins from FAST-J.  Renamed CONC to DUST to avoid conflicts.
!        Also reference NTTLOOP from "comode.h".  Also added parallel
!        DO-loops.  Also renamed MONTH and YEAR to THISMONTH and THISYEAR
!        to avoid conflicts w/ other variables. (bmy, 11/15/01)
!  (10) Bug fix: Make sure to use 1996 dust data for Dec 1995 for the
!        GEOS-STRAT met field dataset.  Set off CASE statement with an
!        #if defined( GEOS_STRAT ) block. (rvm, bmy, 1/2/02)
!  (11) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (12) Now report dust optical depths in ND21 diagnostic at 400 nm.  Now
!       report dust optical depths as one combined diagnostic field instead 
!        of 7 separate fields.  Now reference JLOP from "comode_mod.f".  
!        Now save aerosol surface areas as tracer #5 of the ND21 diagnostic.  
!        (rvm, bmy, 2/28/02)
!  (13) Remove declaration for TIME, since that is also defined in the
!        header file "comode.h" (bmy, 3/20/02)
!  (14) Now read mineral dust files directly from the DATA_DIR/dust_200203/
!        subdirectory (bmy, 4/2/02)
!  (15) Now reference BXHEIGHT from "dao_mod.f".  Also reference ERROR_STOP
!        from "error_mod.f". (bmy, 10/15/02)
!  (16) Now call READ_BPCH2 with QUIET=TRUE to suppress extra informational
!        output from being printed.  Added cosmetic changes. (bmy, 3/14/03)
!  (17) Since December 1997 dust data does not exist, use November 1997 dust
!        data as a proxy. (bnd, bmy, 6/30/03)
!  (18) Bundled into "dust_mod.f" and renamed to RDUST_OFFLINE. (bmy, 4/1/04)
!  (19) Now references DATA_DIR from "directory_mod.f".  Now parallelize over 
!        the L-dimension for ND21 diagnostic. (bmy, 7/20/04)
!  (20) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (21) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE COMODE_MOD,    ONLY : ERADIUS, IXSAVE, IYSAVE
      USE COMODE_MOD,    ONLY : IZSAVE,  JLOP,   TAREA
      USE DAO_MOD,       ONLY : BXHEIGHT
      USE DIAG_MOD,      ONLY : AD21
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      IMPLICIT NONE

#     include "cmn_fj.h"   ! LPAR, CMN_SIZE
#     include "jv_cmn.h"   ! ODMDUST, QAA, RAA
#     include "CMN_DIAG"   ! ND21, LD21
#     include "comode.h"   ! NTTLOOP

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH, THISYEAR

      ! Local variables
      INTEGER             :: I, J, JLOOP, L, N
      INTEGER, SAVE       :: MONTH_LAST = -999
      REAL*4              :: TEMP(IGLOB,JGLOB,LGLOB)
      REAL*8              :: DUST(IIPAR,JJPAR,LLPAR,NDUST)
      REAL*8              :: MSDENS(NDUST), XTAU
      CHARACTER (LEN=255) :: FILENAME

      !=================================================================
      ! RDUST begins here!
      !
      ! Read aerosol data from the binary punch file during the first 
      ! chemistry timestep and, after that, at the start of each month.
      !=================================================================
      IF ( THISMONTH /= MONTH_LAST ) THEN   
         
         ! Save the current month
         MONTH_LAST = THISMONTH

         ! Get TAU0 value used to index the punch file
         ! Use the "generic" year 1985
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
         
         ! Select proper dust file name for GEOS-1, GEOS-3, or GEOS-4
         FILENAME = TRIM( DATA_DIR ) // 'dust_200203/dust.' //
     &              GET_NAME_EXT()   // '.'                 // 
     &              GET_RES_EXT()

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - RDUST: Reading ', a )

         ! Read aerosol concentrations [kg/m3] for each 
         ! dust type from the binary punch file
         DO N = 1, NDUST 
            CALL READ_BPCH2( FILENAME, 'MDUST-$', N,     XTAU,
     &                       IGLOB,     JGLOB,    LGLOB, TEMP, 
     &                       QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, DUST(:,:,:,N) )
         ENDDO

         !==============================================================
         ! Convert concentration [kg/m3] to optical depth [unitless].
         !
         ! ODMDUST = ( 0.75 * BXHEIGHT * CONC * QAA ) / 
         !           ( MSDENS * RAA * 1e-6 )
         ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
         !
         !  Units ==> DUST     [ kg/m3    ]
         !            MSDENS   [ kg/m3    ]
         !            RAA      [ um       ]
         !            BXHEIGHT [ m        ]
         !            QAA      [ unitless ]
         !            ODMDUST  [ unitless ]
         !
         ! NOTES: 
         ! (1) Do the calculation at QAA(4,:) (i.e. 999 nm).          
         !==============================================================
         MSDENS(1) = 2500.0d0
         MSDENS(2) = 2500.0d0
         MSDENS(3) = 2500.0d0
         MSDENS(4) = 2500.0d0
         MSDENS(5) = 2650.0d0
         MSDENS(6) = 2650.0d0
         MSDENS(7) = 2650.0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, NDUST
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ODMDUST(I,J,L,N) = 0.75d0        * BXHEIGHT(I,J,L) * 
     &                         DUST(I,J,L,N) * QAA(4,14+N)     / 
     &                        ( MSDENS(N) * RAA(4,14+N) * 1.0D-6 )
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Echo information
         WRITE( 6, 110 ) 
 110     FORMAT( '     - RDUST: Finished computing optical depths' )

         !==============================================================
         ! Calculate Dust Surface Area
         !
         ! Units ==> DUST     [ kg dust/m^3 air    ]
         !           MSDENS   [ kg dust/m^3 dust   ]
         !           RAA      [ um                 ]
         !           TAREA    [ cm^2 dust/cm^3 air ]
         !           ERADIUS  [ cm                 ]
         !
         ! NOTE: first find volume of dust (cm3 dust/cm3 air), then 
         !       multiply by 3/radius to convert to surface area in cm2
         !  
         ! TAREA(:,1:NDUST) and ERADIUS(:,1:NDUST) are for 
         ! the NDUST FAST-J dust wavelength bins (read into DUST)
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L, N )
         DO N     = 1, NDUST
         DO JLOOP = 1, NTTLOOP

            ! Compute 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)

            ERADIUS(JLOOP,N) = RAA(4,14+N) * 1.0D-4

            TAREA(JLOOP,N)   = 3.D0 / ERADIUS(JLOOP,N) *
     &                         DUST(I,J,L,N) / MSDENS(N)  
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! ND21 Diagnostic: 
         !
         ! Tracer #1: Cloud optical depths    (from "optdepth_mod.f")
         ! Tracer #2: Max Overlap Cld Frac    (from "optdepth_mod.f")
         ! Tracer #3: Random Overlap Cld Frac (from "optdepth_mod.f")
         ! Tracer #4: Dust optical depths at 400 nm (from all size bins)
         ! Tracer #5: Dust surface areas (from all size bins)
         !==============================================================
         IF ( ND21 > 0 ) THEN
            DO N = 1, NDUST

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L ) 
               DO L = 1, LD21
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  !--------------------------------------
                  ! ND21 tracer #4: Dust optical depths
                  !--------------------------------------
                  AD21(I,J,L,4) = AD21(I,J,L,4) + 
     &               ( ODMDUST(I,J,L,N) * QAA(2,14+N) / QAA(4,14+N) )

                  !--------------------------------------
                  ! ND21 tracer #5: Dust surface areas
                  !--------------------------------------
                  IF ( L <= LLTROP ) THEN

                     ! Convert 3-D indices to 1-D index
                     ! JLOP is only defined in the tropopause
                     JLOOP = JLOP(I,J,L)
             
                     ! Add to AD21
                     IF ( JLOOP > 0 ) THEN
                        AD21(I,J,L,5) = AD21(I,J,L,5) + TAREA(JLOOP,N)
                     ENDIF
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ENDDO
         ENDIF 
      ENDIF

      ! Return to calling program
      END SUBROUTINE RDUST_OFFLINE

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DUST
!
!******************************************************************************
!  Subroutine INIT_DUST allocates all module arrays (bmy, 3/30/04)
! 
!  NOTES:
!  (1 ) Now references LDEAD from "logical_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_MOD, ONLY : LDEAD
      USE ERROR_MOD,   ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! Size parameters
      
      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS

      !=================================================================
      ! INIT_DUST begins here!
      !=================================================================

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

!tdf  Double size of IDDEP to account for dust alkalinity   tdf 04/10/08

      ! Drydep flags
!tdf      ALLOCATE( IDDEP( NDSTBIN ), STAT=AS )      
      ALLOCATE( IDDEP( NDSTBIN*2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IDDEP' )
      IDDEP = 0

!tdf  Introduce JDDEP to permit dry deposition frequencies for sulfate
!tdf  to be used when acid dust mass exceeds 1% of dust mass

      ! Drydep flags
      ALLOCATE( JDDEP( NDSTBIN*2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IDDEP' )
      JDDEP = 0

      ! Dust radii
      ALLOCATE( DUSTREFF( NDSTBIN ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DUSTREFF' )
      DUSTREFF(1:NDSTBIN) = (/ 0.73d-6, 1.4d-6,  2.4d-6,  4.5d-6  /)

      ! Dust density
      ALLOCATE( DUSTDEN( NDSTBIN ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DUSTREFF' )
      DUSTDEN(1:NDSTBIN) = (/ 2500.d0, 2650.d0, 2650.d0, 2650.d0 /)

      ! These only have to be allocated for the Ginoux source function
      IF ( .not. LDEAD ) THEN
         ALLOCATE( FRAC_S( NDSTBIN ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRAC_S' )
         FRAC_S(1:NDSTBIN)   = (/ 0.095d0, 0.3d0, 0.3d0, 0.3d0 /)
      ENDIF

      ! Reset flag
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_DUST

!-----------------------------------------------------------------------------

      SUBROUTINE GET_DUST_ALK( I, J, L, ALK_d, KTS, KTN, KTH )
!
!******************************************************************************
!  Function GET_DUST_ALK returns
!           (1) dust alkalinity,   ALK_d(NDSTBIN) [v/v]
!  NB seasalt alkalinity returned from GET_ALK is in [kg]   tdf 04/09/08
!           (2) rate coefficients, KTS(NDSTBIN), KTN(NDSTBIN), for uptake
!               of SO2 and HNO3 on dust for use in sulfate_mod.f for
!               chemistry on dust aerosols.
! tdf 2/13/2K9
!           (3) rate coefficients KTH(NDSTBIN) for uptake of H2SO4
! tdf 08/21/2K9: KTH is now a FRACTION of the size-weighted total
!               area of aerosols in the grid box
!
!  GET_DUST_ALK is analogous to GET_ALK for seasalt
!  (bec, 12/7/04; tdf 04/08/08)
!
!  Arguments as Input:
!  ============================================================================
!
!  NOTES:
!  Arguments as Output:
!  KTS (NDSTBIN): rate coefficients for uptake of SO2 on dust  [s-1]
!  KTN (NDSTBIN): rate coefficients for uptake of HNO3 on dust [s-1]
! tdf 2/13/2K9
!  KTH (NDSTBIN): rate coefficients for uptake of H2SO4 on dust [s-1]
! tdf 08/21/2K9: KTH is now a FRACTION of the size-weighted total
!               area of aerosols in the grid box
!
!  Variables used:
!  STT                    : GEOS-CHEM Tracer array
!  Note in the GET_ALK routine for seasalt, STT is in [v/v]
!  AD                     : Air mass in each grid box [kg]
!  AIRVOL                 : Air volume in each grid box [m^3]
!  TCVV                   : Molecular weight air / Molecular weight tracer
!  STT * AIR_DENS / TCVV  is then
!  [moles/mole] * [kg/m^3] * [kg/mole]/[kg/mole]  = [kg/m^3]
!  Need V, GAMMA_SO2, GAMMA_HNO3
!******************************************************************************
!
! tdf
      USE DAO_MOD,      ONLY : AD, RH, AIRVOL
      USE TRACER_MOD,   ONLY : STT, TCVV
      USE ERROR_MOD,    ONLY : IT_IS_NAN

      ! Include dust tracers                            tdf 04/10/08
      USE TRACERID_MOD, ONLY : IDTDST1, IDTDST2, IDTDST3, IDTDST4

      ! Include dust alkalinity tracers                 tdf 04/10/08
      USE TRACERID_MOD, ONLY : IDTDAL1, IDTDAL2, IDTDAL3, IDTDAL4

!tdf for TAREA
      USE COMODE_MOD,    ONLY : ERADIUS, IXSAVE, IYSAVE
      USE COMODE_MOD,    ONLY : IZSAVE,  JLOP,   TAREA
!tdf for TAREA

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L

      ! Return value
      REAL*8, INTENT(OUT)  :: ALK_d(NDSTBIN)                 ! [kg]
      REAL*8, INTENT(OUT)  :: KTS  (NDSTBIN)                 ! [s-1]
      REAL*8, INTENT(OUT)  :: KTN  (NDSTBIN)                 ! [s-1]
! tdf 08/21/2K9
      REAL*8, INTENT(OUT)  :: KTH  (NDSTBIN)                 ! [fraction]

      ! Local
! tdf Dust Surface Areas                               ! tdf 08/20/09
      REAL*8             :: AREA_d(NDSTBIN)            ! [cm^2/cm^3]
! tdf Dust Surface Areas within sub-bins 1-4 of BIN 1  ! tdf 08/20/09
      REAL*8             :: AREA_sd1(4)                ! [cm^2/cm^3]
! tdf Dust Effective Radius                            ! tdf 08/20/09
      REAL*8             :: RD_d(NDSTBIN)              ! [cm]
! tdf Dust Effective Radii for sub-bins 1-4 of BIN 1   ! tdf 08/20/09
      REAL*8             :: RD_sd1(4)                  ! [cm]
! tdf Dust size-weighted Surface Areas                 ! tdf 08/20/09
      REAL*8             :: DF_AREA_d(NDSTBIN)         ! [1/s]
! tdf Dust size-weighted Surface Areas for sub-bins 1-4! tdf 08/20/09
      REAL*8             :: DF_AREA_sd1(4)             ! [1/s]
!tdf
      REAL*8,  PARAMETER :: PI = 3.14159265
      REAL*8             :: N1, KT1, KT1N
      REAL*8             :: AREA, HGF
      REAL*8             :: CONST1, CONST2, CONST
      REAL*8             :: AIR_DENS
      REAL*8             :: A1, B1, A1N, B1N, T1, R1
      REAL*8             :: DD, RD, DM                 ! tdf 04/16/08
      REAL*8             :: TOTAL_AREA, DF_TOTAL_AREA  ! tdf 08/20/09
      REAL*8,  PARAMETER :: MINDAT = 1.d-20
!
      INTEGER            :: IRH
      INTEGER            :: IBIN, ISBIN
      INTEGER            :: JLOOP
! tdf
      REAL*8             :: DST_d (NDSTBIN), ALK, GAMS, GAMN
      REAL*8             :: SULF_AREA, BC_AREA, OC_AREA
      REAL*8             :: SULF_RAD, BC_RAD, OC_RAD
      REAL*8             :: SULF_FAC, BC_FAC, OC_FAC
      REAL*8             :: SSA_AREA, SSC_AREA
      REAL*8             :: SSA_RAD, SSC_RAD
      REAL*8             :: SSA_FAC, SSC_FAC
! tdf

!  Need this for dust
!      REAL*8,  PARAMETER   :: GAMMA_SO2 = 0.05d0  !(Song & Carmichael, 2001)
! 200 times smaller 8/28/2K9
      REAL*8,  PARAMETER   :: GAMMA_SO2 = 2.5d-4

!tdf V9 4/1/2K9 Applying Song et al.(2007) reduced value 
!tdf 
      REAL*8,  PARAMETER   :: GAMMA_H2SO4 = 1.d0  

!  Need this for dust
! tdf      REAL*8,  PARAMETER   :: GAMMA_HNO3 = 0.1d0 ! (Song & Carmichael, 2001)
! 200 times smaller    8/28/2K9
      REAL*8,  PARAMETER   :: GAMMA_HNO3 = 5.0d-4 

      REAL*8,  PARAMETER   :: DG = 0.2d0 !gas phase diffusion coeff. [cm2/s]
      REAL*8,  PARAMETER   :: v = 3.0d4 !cm/s

      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! GET_DUST_ALK begins here!
      !=================================================================

      ! Zero variables
      DO IBIN = 1, NDSTBIN
         ALK_d (IBIN) = 0.0D0
         KTS   (IBIN) = 0.0D0
         KTN   (IBIN) = 0.0D0
         KTH   (IBIN) = 0.0D0
         AREA_d (IBIN) = 0.D0
         RD_d (IBIN) = 0.D0
      END DO

      AIR_DENS = AD(I,J,L) / AIRVOL(I,J,L)        ! air density [kg/m3]

! Retrieve Dust Alkalinity [v/v] from STT array
      ALK_d(1) = STT(I,J,L,IDTDAL1)
      ALK_d(2) = STT(I,J,L,IDTDAL2)
      ALK_d(3) = STT(I,J,L,IDTDAL3)
      ALK_d(4) = STT(I,J,L,IDTDAL4)

! Dust [kg/m3] from STT, used to compute dust surface area
! Units: (moles/mole).(kg(air)/m3).(kg(dust)/mole)/(kg(air)/mole)
      DST_d(1) = STT(I,J,L,IDTDST1) * AIR_DENS/TCVV(IDTDST1)
      DST_d(2) = STT(I,J,L,IDTDST2) * AIR_DENS/TCVV(IDTDST2)
      DST_d(3) = STT(I,J,L,IDTDST3) * AIR_DENS/TCVV(IDTDST3)
      DST_d(4) = STT(I,J,L,IDTDST4) * AIR_DENS/TCVV(IDTDST4)

! tdf Now get aerosol surface area from TAREA (cm2/cm3)
      JLOOP = JLOP(I,J,L)
      SULF_AREA = TAREA(JLOOP,NDUST+1)
      BC_AREA   = TAREA(JLOOP,NDUST+2)
      OC_AREA   = TAREA(JLOOP,NDUST+3)
      SSA_AREA  = TAREA(JLOOP,NDUST+4)
      SSC_AREA  = TAREA(JLOOP,NDUST+5)

! tdf Now get aerosol effective radius from ERADIUS (cm)
      SULF_RAD  = ERADIUS(JLOOP,NDUST+1)
      BC_RAD    = ERADIUS(JLOOP,NDUST+2)
      OC_RAD    = ERADIUS(JLOOP,NDUST+3)
      SSA_RAD   = ERADIUS(JLOOP,NDUST+4)
      SSC_RAD   = ERADIUS(JLOOP,NDUST+5)

! tdf Quotients [s/cm] used to weight surface area for H2SO4 uptake
      SULF_FAC = (SULF_RAD / DG + 4./(V*GAMMA_H2SO4) )
      BC_FAC   = (  BC_RAD / DG + 4./(V*GAMMA_H2SO4) )
      OC_FAC   = (  OC_RAD / DG + 4./(V*GAMMA_H2SO4) )
      SSA_FAC  = ( SSA_RAD / DG + 4./(V*GAMMA_H2SO4) )
      SSC_FAC  = ( SSC_RAD / DG + 4./(V*GAMMA_H2SO4) )

!tdf Surface areas and effective radii for sub-bins 1-4 of dust bin 1
      DO ISBIN = 1, 4
         T1 = TAREA  (JLOOP,ISBIN)
         R1 = ERADIUS(JLOOP,ISBIN)
         AREA_sd1    (ISBIN) = T1
         RD_sd1      (ISBIN) = R1
!tdf surface area for sub bins 1-4 in bin 1, weighted by gas-phase
!tdf diffusion and collision limitations
!tdf used to compute proportionate uptake of H2SO4 only  [1/s]
         DF_AREA_sd1 (ISBIN) = T1
     +                       / (R1/DG + 4.0D0/(V*GAMMA_H2SO4))
      END DO

!-----------------------------------------------------------------------
! Very Simple Formulation: For each size bin (i)   ! tdf 8/20/09
! Dust Area density = 3 * Dust Mass density  / (REFF(i) * DUSTDEN)
! TAREA computed   in RDUST_ONLINE - Units: cm^2(dust) / cm^3(air) 
! ERADIUS computed in RDUST_ONLINE - Units: cm
! NB: I am now subdividing the submicron dust size bin 
!     using TAREA (JLKOOP,1->4), and ERADIUS (JLOOP,1->4).
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------
!  Find Dust surface area density in grid-box, AREA_d [cm^2/cm^3].
!  Also find the size-weighted surface area density, DF_AREA_d [1/s].
!  The latter represents the gas-phase diffusion and surface 
!  limited weighting and is used to determine the fraction of H2SO4 
!  taken up on dust, versus taken up on other aerosols.
!                                                 tdf 08/21/09
!-------------------------------------------------------------------------

! tdf Loop over size bins  (NDSTBIN = 4)
      DO IBIN = 1, NDSTBIN

! Dust Area density in grid box,      AREA_d [cm^2/cm^3]    tdf 8/21/09
! Dust weighted surface area density, DF_AREA_d [1/s]       tdf 8/21/09

         IF (IBIN .EQ. 1) THEN
      ! For Dust size bin 1, sum over the 4 size sub bins  tdf 8/21/09
            AREA_d   (IBIN) = AREA_sd1(1) + AREA_sd1(1)       ![cm^2/cm^3]
     +                      + AREA_sd1(3) + AREA_sd1(4) 
            DF_AREA_d(IBIN) = DF_AREA_sd1(1) + DF_AREA_sd1(3) ![1/s]
     +                      + DF_AREA_sd1(2) + DF_AREA_sd1(4)
         ELSE
            T1 = TAREA(JLOOP,3+IBIN)      ! [cm^2/cm^3]
            R1 = ERADIUS(JLOOP,3+IBIN)    ! [cm]
            AREA_d   (IBIN) = T1          ! [cm^2/cm^3]
            DF_AREA_d(IBIN) = T1          ! [1/s]
     +                       / (R1/DG + 4.0D0/(V*GAMMA_H2SO4))
         ENDIF

      END DO
! tdf
         if (I .eq. 1 .and. J .eq. 63 .and. L .eq. 6) then
           print *,' GET_DUST_ALK: AREA'
           print *,' SULF_AREA ',SULF_AREA,'cm2/cm3'
           print *,'   BC_AREA ',BC_AREA,'cm2/cm3'
           print *,'   OC_AREA ',OC_AREA,'cm2/cm3'
           print *,'  SSA_AREA ',SSA_AREA,'cm2/cm3'
           print *,'  SSC_AREA ',SSC_AREA,'cm2/cm3'
           print *,'  AREA_d(1) ',AREA_d(1),'cm2/cm3'
           print *,'  AREA_d(2) ',AREA_d(2),'cm2/cm3'
           print *,'  AREA_d(3) ',AREA_d(3),'cm2/cm3'
           print *,'  AREA_d(4) ',AREA_d(4),'cm2/cm3'
         endif
! tdf
! tdf total aerosol surface area  [cm^2/cm^3]
      TOTAL_AREA = SULF_AREA + BC_AREA + OC_AREA
     +           + SSA_AREA  + SSC_AREA
     +           + AREA_d(1) + AREA_d(2) + AREA_d(3) + AREA_d(4)

! tdf total surface area weighted by gas-phase diffusion limitation [1/s]
      DF_TOTAL_AREA = SULF_AREA / SULF_FAC + BC_AREA / BC_FAC
     +              + OC_AREA / OC_FAC + SSA_AREA / SSA_FAC
     +              + SSC_AREA / SSC_FAC
     +              + DF_AREA_d(1) 
     +              + DF_AREA_d(2) 
     +              + DF_AREA_d(3) 
     +              + DF_AREA_d(4) 

         if (I .eq. 1 .and. J .eq. 63 .and. L .eq. 6) then
           print *,' GET_DUST_ALK: DF_AREA'
           print *,' DF_SULF_AREA  ',SULF_AREA/SULF_FAC,'1/s'
           print *,'   DF_BC_AREA  ',BC_AREA/BC_FAC,'1/s'
           print *,'   DF_OC_AREA  ',OC_AREA/OC_FAC,'1/s'
           print *,'  DF_SSA_AREA  ',SSA_AREA/SSA_FAC,'1/s'
           print *,'  DF_SSC_AREA  ',SSC_AREA/SSC_FAC,'1/s'
           print *,'  DF_AREA_d(1) ',DF_AREA_d(1),'1/s'
           print *,'  DF_AREA_d(2) ',DF_AREA_d(2),'1/s'
           print *,'  DF_AREA_d(3) ',DF_AREA_d(3),'1/s'
           print *,'  DF_AREA_d(4) ',DF_AREA_d(4),'1/s'
         endif

! tdf Total Dust Alkalinity
      ALK = ALK_d(1) + ALK_d(2) + ALK_d(3) + ALK_d(4)  ! [v/v]

      ! set humidity index IRH as a percent
         IRH = RH(I,J,L)
         IRH = MAX(  1, IRH )
         IRH = MIN( 99, IRH )

      ! hygroscopic growth factor for dust: Set to NO GROWTH for now
      IF ( IRH < 100 ) HGF = 1.0d0

! tdf Loop over size bins (NDSTBIN = 4)

      DO IBIN = 1, NDSTBIN

      !----------------------------------
      ! SO2 uptake onto particles
      !----------------------------------

!tdf 2/11/2K9
!tdf Following relative uptake rates of Preszler-Prince et al.(2007)
         IF (IRH >= 90 ) THEN     
           GAMS = GAMMA_SO2 * 2.
         ELSE IF (IRH >= 84 ) THEN  
           GAMS = GAMMA_SO2 * (0.5  + 1.5*(IRH-84.)/(90.-84.))
         ELSE IF (IRH >= 76 ) THEN 
           GAMS = GAMMA_SO2 * (0.16 + 0.34*(IRH-76.)/(84.-76.))
         ELSE IF (IRH >= 33 ) THEN  
           GAMS = GAMMA_SO2 * (0.03 + 0.13*(IRH-33.)/(76.-33.))
         ELSE IF (IRH >= 20 ) THEN 
           GAMS = GAMMA_SO2*0.03
         ELSE                        ! 0.0 below 20%
           GAMS = GAMMA_SO2*0.0
         ENDIF

      ! Check for sufficient alkalinity          tdf 3/28/2K8
         IF ( ALK > MINDAT ) THEN

         ! calculate gas-to-particle rate constant for uptake of
         ! SO2 onto dust aerosols [Jacob, 2000] analytical solution
         ! Corrected based on discussions with Becky     tdf 07/14/08
! tdf
            KT1    = 0.d0
            IF (IBIN .EQ. 1) THEN
! tdf Sum over the 1-4 sub-bins for bin 1      ! tdf 08/21/2K9
              DO ISBIN = 1, 4
                 RD = RD_sd1 (ISBIN)           ! effective radius [cm]
                 AREA = AREA_sd1 (ISBIN)       ! Dust Surface Area [cm^2/cm^3]
                 CONST1 = 4.D0/(V*GAMS)        ! Collision [s/cm]
                 CONST2 = RD/DG                ! Diffusion [s/cm]
                 CONST  = CONST1 + CONST2
                 KT1    = KT1 + AREA / CONST   ! [cm^2/cm^3] * [cm/s] = [1/s]
              END DO
            ELSE
               RD     = RD_d (IBIN)            ! effective radius [cm]
               AREA   = AREA_d (IBIN)          ! Dust Surface Area [cm^2/cm^3]
               CONST1 = 4.D0/(V*GAMS)          ! Collision [s/cm]
               CONST2 = RD/DG                  ! Diffusion [s/cm]
               CONST  = CONST1 + CONST2
               KT1    = AREA / CONST           ! [cm^2/cm^3] * [cm/s] = [1/s]
            ENDIF

            KTS(IBIN) = KT1

         ELSE
         ! If no alkalinity, set rate coefficients to zero
!tdf
            KTS(IBIN)  = 0.D0

         ENDIF

      !----------------------------------
      ! HNO3 uptake onto particles
      !----------------------------------

!tdf 2/11/2K9
!tdf Following uptake coefficients of Liu et al.(2007)
         IF (IRH >= 80 ) THEN         
           GAMN = GAMMA_HNO3 * 2.1
         ELSE IF (IRH >= 70 ) THEN   
           GAMN = GAMMA_HNO3 * (1.3  + 0.7 * (IRH - 70.)/ 10.)
         ELSE IF (IRH >= 60 ) THEN 
           GAMN = GAMMA_HNO3 * (1.0  + 0.3 * (IRH - 60.)/  10.)
         ELSE IF (IRH >= 50 ) THEN
           GAMN = GAMMA_HNO3 * (0.7  + 0.3 * (IRH - 50.)/ 10.)
         ELSE IF (IRH >= 30 ) THEN  
           GAMN = GAMMA_HNO3 * (0.19 + 0.255 *(IRH - 30.)/ 10.)
         ELSE IF (IRH >= 10 ) THEN  
           GAMN = GAMMA_HNO3 * (0.03 + 0.08 *(IRH - 10.)/ 10.)
         ELSE                     
           GAMN = GAMMA_HNO3*0.0
         ENDIF
!tdf 2/11/2K9

      ! Check for sufficient alkalinity      tdf 3/28/2K8
         IF ( ALK > MINDAT ) THEN

         ! calculate gas-to-particle rate constant for uptake of
         ! HNO3 onto dust aerosols [Jacob, 2000] analytical solution
         ! Corrected based on discussions with Becky     tdf 07/14/08

            KT1    = 0.d0
            IF (IBIN .EQ. 1) THEN
! tdf Sum over the 1-4 sub-bins for bin 1      ! tdf 08/21/2K9
              DO ISBIN = 1, 4
                 RD = RD_sd1 (ISBIN)           ! effective radius [cm]
                 AREA = AREA_sd1 (ISBIN)       ! Dust Surface Area [cm^2/cm^3]
                 CONST1 = 4.D0/(V*GAMN)          ! Collision [s/cm]
                 CONST2 = RD/DG                  ! Diffusion [s/cm]
                 CONST  = CONST1 + CONST2
                 KT1    = KT1 + AREA / CONST     ! [cm^2/cm^3] * [cm/s] = [1/s]
              END DO
            ELSE
              CONST1 = 4.D0/(V*GAMN)    ! Collision [s/cm]
              CONST2 = RD/DG                  ! Diffusion [s/cm]
              CONST  = CONST1 + CONST2
              KT1    = AREA / CONST           ! [cm^2/cm^3] * [cm/s] = [1/s]
            ENDIF

            KTN(IBIN) = KT1

         ELSE
         ! If no alkalinity, set rate coefficients to zero
!tdf
            KTN(IBIN)  = 0.D0

         ENDIF

      !----------------------------------
      ! H2SO4 uptake onto particles
      !----------------------------------

      ! Uptake not limited by dust alkalinity      tdf 3/02/2K9

!tdf As of 08/20/09, we use AREA and size weighted uptake
!tdf where now KTH is a fractional uptake for each size bin
!tdf with respect to total aerosol surface area.

         KT1    = DF_AREA_d(IBIN) / DF_TOTAL_AREA   ! Fraction

         KTH(IBIN) = KT1

! tdf End Loop over size bins
      END DO

      ! Return arrays: ALK_d, KTS, KTN, KTH

      ! Return to calling program
      END SUBROUTINE GET_DUST_ALK

!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DUST
!
!******************************************************************************
!  Subroutine CLEANUP_DUST deallocates all module arrays (bmy, 3/30/04)
! 
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DUST begins here!
      !=================================================================
      IF ( ALLOCATED( IDDEP    ) ) DEALLOCATE( IDDEP    )      
      IF ( ALLOCATED( FRAC_S   ) ) DEALLOCATE( FRAC_S   )
      IF ( ALLOCATED( DUSTREFF ) ) DEALLOCATE( DUSTREFF )
      IF ( ALLOCATED( DUSTDEN  ) ) DEALLOCATE( DUSTDEN  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DUST

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DUST_MOD
