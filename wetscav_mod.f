! $Id: wetscav_mod.f,v 1.6 2004/04/13 14:52:33 bmy Exp $
      MODULE WETSCAV_MOD
!
!******************************************************************************
!  Module WETSCAV_MOD contains arrays for used in the wet scavenging of
!  tracer in cloud updrafts, rainout, and washout. (bmy, 2/28/00, 4/5/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) NSOLMAX (INTEGER) : Max # of soluble tracers       [unitless]
!  (2 ) NSOL    (INTEGER) : Actual # of soluble tracers    [unitless]
!  (3 ) IDWETD  (INTEGER) : Index array for WETDEP routine [unitless]
!  (4 ) Vud     (REAL*8 ) : Array for updraft velocity     [m/s]
!  (5 ) CLDLIQ  (REAL*8 ) : Array for cloud liquid water   [cm3 H2O/cm3 air]
!  (6 ) CLDICE  (REAL*8 ) : Array for cloud ice content    [cm3 ice/cm3 air]
!  (7 ) C_H2O   (REAL*8 ) : Array for Mixing ratio of ,      
!                            water, computed from Eice(T)  [v/v]
!  (8 ) PDOWN   (REAL*8 ) : Precip thru bottom of grid box [cm3 H2O/cm2 area/s]
!  (9 ) QQ      (REAL*8 ) : Rate of new precip formation   [cm3 H2O/cm3 air/s]
!  (10) EPSILON (REAL*8 ) : A very small positive number   [unitless]
!  (11) H2O2s   (REAL*8 ) : Array to save H2O2 for wetdep  [v/v]
!  (12) SO2s    (REAL*8 ) : Array to save SO2 for wetdep   [v/v]
!
!  Module Routines:
!  ============================================================================
!  (1 ) MAKE_QQ           : Constructs the QQ field (precipitable water)
!  (2 ) E_ICE             : Computes saturation vapor pressure for ice 
!  (3 ) COMPUTE_L2G       : Computes the ratio [v/v liquid] / [v/v gas] 
!  (4 ) COMPUTE_F         : Computes fraction of tracer lost in cloud updrafts
!  (5 ) GET_ISOL          : Returns correct index for ND38 diagnostic
!  (6 ) F_AEROSOL         : Internal routine to COMPUTE_F
!  (7 ) GET_ISOL          : Internal routine to COMPUTE_F
!  (8 ) RAINOUT           : Computes fraction of soluble tracer lost to rainout
!  (9 ) GET_RAINFRAC      : Internal routine to RAINOUT
!  (10) WASHOUT           : Computes fraction of soluble tracer lost to washout
!  (11) WASHFRAC_AEROSOL  : Internal routine to WASHOUT
!  (12) WASHFRAC_LIQ_GAS  : Internal routine to WASHOUT
!  (13) WETDEP            : Driver routine for computing wet deposition losses
!  (14) LS_K_RAIN         : Computes K_RAIN (for LS precipitation)
!  (15) LS_F_PRIME        : Computes F_PRIME (for LS precipitation)
!  (16) CONV_F_PRIME      : Computes F_PRIME (for convective precipitation)
!  (17) SAFETY            : Stops WETDEP w/ error msg if negative tracer found
!  (18) WETDEPID          : Initalizes the IDWETD array for routine WETDEP
!  (20) GET_WETDEP_NMAX   : Returns max # of soluble tracers per simulation
!  (21) GET_WETDEP_NSOL   : Returns actual # of soluble tracers per simulation
!  (21) GET_WETDEP_IDWETD : Returns CTM tracer # of for a given wetdep species 
!  (22) INIT_WETSCAV      : Initializes fields used for computing wetdep losses
!  (23) CLEANUP_WETSCAV   : Deallocates all allocatable module arrays
!
!  GEOS-CHEM modules referenced by wetscav_mod.f
!  ============================================================================
!  (1 ) dao_mod.f      : Module containing arrays for DAO met fields
!  (2 ) diag_mod.f     : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) error_mod.f    : Module containing NaN and other error check routines
!  (4 ) pressure_mod.f : Module containing routines to compute P(I,J,L)
!  (5 ) tracerid_mod.f : Module containing pointers to tracers and emissions
!
!  References:
!  ============================================================================
!  (1 ) Liu,H., D.J. Jacob, I. Bey and R.M. Yantosca, "Constraints from 210Pb 
!        and 7Be on wet deposition and transport in a global three-dimensional
!        chemical tracer model driven by assimilated meteorological fields", 
!        JGR, Vol 106, pp 12109-12128, 2001.
!  (2 ) D.J. Jacob, H. Liu, C. Mari, and R. M. Yantosca, "Harvard wet 
!        deposition scheme for GMI", Harvard Atmospheric Chemistry Modeling 
!        Group, March 2000.
!  (3 ) Chin, M., D.J. Jacob, G.M. Gardner, M.S. Foreman-Fowler, and P.A. 
!        Spiro, "A global three-dimensional model of tropospheric sulfate", 
!        J. Geophys. Res., 101, 18667-18690, 1996.
!  (4 ) Balkanski, Y  D.J. Jacob, G.M. Gardner, W.C. Graustein, and K.K.
!        Turekian, "Transport and Residence Times of Tropospheric Aerosols
!        from a Global Three-Dimensional Simulation of 210Pb", JGR, Vol 98, 
!        (D11) pp 20573-20586, 1993.  
!  (5 ) Giorgi, F, & W.L. Chaimedes, "Rainout Lifetimes of Highly Soluble
!        Aerosols and Gases as Inferred from Simulations With a General
!        Circulation Model", JGR, Vol 86 (D13) pp 14367-14376, 1986.  
!
!  NOTES:
!  (1 ) Now trap allocation errors with routine ALLOC_ERR. (bmy, 7/11/00)
!  (2 ) Moved routine MAKE_QQ here from "dao_mod.f" (bmy, 10/12/00)
!  (3 ) Reordered arguments in INIT_PRECIP (bmy, 10/12/00)
!  (4 ) Updated comments (bmy, 9/4/01)
!  (5 ) Bug fix in MAKE_QQ: BXHEIGHT is sized IIPAR,JJPAR,LLPAR (bmy, 10/4/01)
!  (6 ) Removed obsolete, commented-out code from 10/01 (bmy, 11/26/01)
!  (7 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (8 ) Now zero allocatable arrays (bmy, 8/5/02)
!  (9 ) Bug fix: ND39 diagnostic now closes the budget.  Also bundled several
!        standalone routines into this module.  Now references F90 module
!        "tracerid_mod.f".  Also set NSOLMAX=10 since we now have sulfate
!        tracers for wetdep.   Now prevent out-of-bounds errors in routine
!        WETDEP.  Added GET_WETDEP_NMAX function to return max # of soluble
!        tracers for allocating diagnostic arrays.  Added functions 
!        GET_WETDEP_NSOL and GET_WETDEP_IDWETD.  Now init H2O2s and SO2s
!        to the initial H2O2 and SO2 from STT.  Updated comments. 
!        (qli, bmy, 1/14/03)
!  (10) Improvements for SO2/SO4 scavenging (rjp, bmy, 3/23/03)
!  (11) Now references "time_mod.f".  Added driver routine DO_WETDEP to
!        remove cumbersome calling sequence from MAIN program.  Also declared
!        WETDEP and MAKE_QQ PRIVATE to this module. (bmy, 3/27/03)
!  (11) Add parallelization to routine WETDEP (bmy, 3/17/04)
!  (12) Added carbon and dust aerosol tracers (rjp, tdf, bmy, 4/5/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS
      !=================================================================

      ! PRIVATE routines
      PRIVATE :: COMPUTE_L2G, E_ICE,    RAINOUT
      PRIVATE :: WASHOUT,     WETDEPID, MAKE_QQ
      PRIVATE :: WETDEP,      GET_ISOL

      ! PRIVATE variables
      PRIVATE :: IDWETD, NSOLMAX, NSOL,  Vud, C_H2O
      PRIVATE :: CLDLIQ, CLDICE,  PDOWN, QQ,  EPSILON

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      !---------------------------------------------------------------------
      ! Prior to 4/5/04:
      ! Adjust NSOLMAX for carbon & dust aerosol tracers (rjp, bmy, 4/5/04)
      !INTEGER, PARAMETER   :: NSOLMAX = 10
      !---------------------------------------------------------------------
      INTEGER, PARAMETER   :: NSOLMAX = 18
      INTEGER              :: NSOL 
      INTEGER              :: IDWETD(NSOLMAX)
      REAL*8,  ALLOCATABLE :: Vud(:,:)
      REAL*8,  ALLOCATABLE :: C_H2O(:,:,:)
      REAL*8,  ALLOCATABLE :: CLDLIQ(:,:,:)
      REAL*8,  ALLOCATABLE :: CLDICE(:,:,:)
      REAL*8,  ALLOCATABLE :: PDOWN(:,:,:)
      REAL*8,  ALLOCATABLE :: QQ(:,:,:)
      REAL*8,  ALLOCATABLE :: H2O2s(:,:,:)
      REAL*8,  ALLOCATABLE :: SO2s(:,:,:)
      REAL*8,  PARAMETER   :: EPSILON = 1d-32

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_WETDEP
!
!******************************************************************************
!  Subroutine DO_WETDEP is a driver for the wet deposition code, called
!  from the MAIN program. (bmy, 3/27/03)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : DEBUG_MSG

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPRT for debug

      !=================================================================
      ! DO_WETDEP begins here!
      !=================================================================

      ! Wetdep by large-scale (stratiform) precip
      CALL MAKE_QQ( .TRUE. )
      CALL WETDEP(  .TRUE. )
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_WETDEP: after LS wetdep' )

      ! Wetdep by convective precip
      CALL MAKE_QQ( .FALSE. )
      CALL WETDEP(  .FALSE. )
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_WETDEP: after conv wetdep' )

      ! Return to calling program
      END SUBROUTINE DO_WETDEP

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_QQ( LS )
!
!*****************************************************************************
!  Subroutine MAKE_QQ computes the large-scale or convective precipitation
!  fields for use with wetdep.f. (hyl, bmy, 2/29/00, 11/8/02)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) LS       : = T for Large-scale precip, =F otherwise
!
!  DAO met fields from "dao_mod.f:"
!  ===========================================================================
!  (1 ) AIRDEN   : Density of air in grid box (I,J,L) [kg air/m^3]
!  (2 ) BXHEIGHT : Height of grid box (I,J,L) in [m]
!  (3 ) MOISTQ   : DAO field for change in specific    
!                  humidity due to moist processes    [kg H2O/kg air/s]
!  (4 ) PREACC   : DAO total accumulated precipitaton [mm/day]
!  (5 ) PRECON   : DAO convective precipitation       [mm/day]
!
!  References (see above for full citations):
!  ===========================================================================
!  (1 ) Liu et al, 2000
!  (2 ) Jacob et al, 2000
!
!  NOTES:
!  (1 ) Now we partition MOISTQ into large-scale and convective parts, using
!        total precipitation PREACC and convective precipitation PRECON (both
!        are vertical integral amounts). The precipitation field at altitudes
!        (PDOWN) is also made (hyl, djj, 10/17/98).
!  (2 ) MAKE_QQ is written in Fixed-Form Fortran 90. (bmy, 4/2/99)!
!  (3 ) AIRDEN, MOISTQ, QQ, and PDOWN are dimensioned (LLPAR,IIPAR,JJPAR) 
!       in order to maximize loop efficiency when processing an (I,J) 
!       column layer by layer. (bmy, 3/14/00)
!  (4 ) MOISTQ is originally [g H2O/kg air/day], and is converted in
!        READ_A6 to [kg H2O/kg air/s]. (bmy, 3/14/00)
!  (5 ) Now reference PREACC, PRECON from "dao_mod.f" instead of from
!        common block header file "CMN_PRECIP" (bmy, 6/26/00)
!  (6 ) Now pass BXHEIGHT as an argument.  Also added to "dao_mod.f". 
!        (bmy, 6/26/00)
!  (7 ) Moved from "dao_mod.f" to "wetscav_mod.f".  Also made PREACC
!        and PRECON into arguments. (bmy, 10/12/00)
!  (8 ) Updated comments (bmy, 9/4/01)
!  (9 ) BXHEIGHT is now sized (IIPAR,JJPAR,LLPAR) (bmy, 10/4/01)
!  (10) Removed obsolete, commented-out code from 10/01 (bmy, 11/26/01)
!  (11) Now reference met field arrays directly from "dao_mod.f" (bmy, 11/8/02)
!******************************************************************************
! 
      ! References to F90 modules
      USE DAO_MOD,   ONLY : AIRDEN, BXHEIGHT, MOISTQ, PREACC, PRECON
      USE ERROR_MOD, ONLY : ALLOC_ERR
      
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      LOGICAL, INTENT(IN)  :: LS

      ! Local variables
      INTEGER              :: I, J, L, AS
      REAL*8               :: PTEMP, FRAC
      LOGICAL              :: FIRST = .TRUE.

      !=================================================================
      ! MAKE_QQ begins here!
      !=================================================================
      IF ( FIRST ) THEN

         ! Allocate PDOWN on first call
         ALLOCATE( PDOWN( LLPAR, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PDOWN' )
         PDOWN = 0d0
      
         ! Allocate QQ on first call
         ALLOCATE( QQ( LLPAR, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'QQ' )
         QQ = 0d0
         
         ! Reset flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Loop over surface grid boxes
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !==============================================================
         ! If there is total precipitation in the (I,J) column, then:
         ! 
         ! (1) Compute FRAC, the large scale fraction (if LS = .TRUE.) 
         !     or convective fraction (if LS = .FALSE.) total 
         !     precipitation.  FRAC is computed from PREACC and PRECON.
         !
         ! (2) Compute QQ, the rate of formation of precipitation 
         !     [cm3 H2O/cm3 air/s].  From MOISTQ [kg H2O/kg air/s], 
         !     the unit conversion is: 
         !
         !     kg H2O   |   m^3 H2O   | AIRDEN kg air         m^3 H2O
         !  ------------+-------------+--------------- ==> -------------   
         !   kg air * s | 1000 kg H2O |    m^3 air          m^3 air * s
         !
         ! and
         !
         !         m^3 H2O                         cm^3 H2O
         !      -------------  is equivalent to  -------------- 
         !       m^3 air * s                      cm^3 air * s!
         !       
         ! since the same conversion factor (10^6 cm^3/m^3) is in both
         ! the numerator and the denominator.
         !
         ! Therefore, the equation for QQ is:
         !
         !   QQ(L,I,J) = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1000.0
         !     
         ! (3) Compute PDOWN, the column precipitation 
         !     [cm3 H2O/cm2 air/s], by multiplying QQ(L,I,J) by 
         !     BXHEIGHT(I,J,L) * 100 cm.  
         !
         ! (4) The reason why we do not force PTEMP to be positive is 
         !     that PREACC is the integral of the MOISTQ field.  MOISTQ 
         !     contains both negative (evap) and positive (precip) 
         !     values.  If we forced PTEMP to be positive, then we would
         !     be adding extra precipitation to PDOWN (hyl, bmy, 3/6/99).
         !==============================================================
         IF ( PREACC(I,J) > 0d0 ) THEN

            ! Large scale or convective fraction of precipitation
            IF ( LS ) THEN
               FRAC = ( PREACC(I,J) - PRECON(I,J) ) / PREACC(I,J) 
            ELSE
               FRAC = PRECON(I,J) / PREACC(I,J)
            ENDIF

            ! Start at the top of the atmosphere
            L = LLPAR

            ! Compute QQ and PDOWN.  Keep PTEMP for the next level
            QQ(L,I,J)    = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1d3
            PTEMP        = QQ(L,I,J) * BXHEIGHT(I,J,L) * 1d2
            PDOWN(L,I,J) = PTEMP

            ! PDOWN cannot be negative
            IF ( PDOWN(L,I,J) < 0d0 ) PDOWN(L,I,J) = 0.d0

            ! Loop down from LLPAR to the surface
            DO L = LLPAR-1, 1, -1
               
               ! Compute QQ and PDOWN.  Keep PTEMP for the next level.
               QQ(L,I,J)    = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1d3
               PDOWN(L,I,J) = PTEMP + QQ(L,I,J) * BXHEIGHT(I,J,L) * 1d2  
               PTEMP        = PDOWN(L,I,J)

               ! PDOWN cannot be negative
               IF ( PDOWN(L,I,J) < 0.0d0 ) PDOWN(L,I,J) = 0.d0
            ENDDO
  
         !==============================================================
         ! If there is no precipitation reaching the surface in the 
         ! (I,J) column, then assume any precipitation at altitude to 
         ! be large-scale.
         ! 
         ! (1) Assume the large scale fraction = 1d0, 
         !                convective fraction  = 0d0
         ! (2) Compute QQ as described above
         ! (3) Compute PDOWN as described above
         !==============================================================
         ELSE

            ! Assume large-scale precipitation!
            IF ( LS ) THEN
               FRAC = 1d0
            ELSE         
               FRAC = 0d0
            ENDIF

            ! Start at the top of the atmosphere
            L = LLPAR

            ! Compute QQ and PDOWN.  Keep PTEMP for the next level
            QQ(L,I,J)    = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1d3
            PTEMP        = QQ(L,I,J) * BXHEIGHT(I,J,L) * 1d2
            PDOWN(L,I,J) = PTEMP
           
            ! PDOWN cannot be negative
            IF( PDOWN(L,I,J) < 0d0 ) PDOWN(L,I,J) = 0.d0

            ! Loop down from LLPAR to the surface
            DO L = LLPAR-1, 1, -1
              
               ! Compute QQ and PDOWN.  Keep PTEMP for the next level
               QQ(L,I,J)    = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1d3
               PDOWN(L,I,J) = PTEMP + QQ(L,I,J) * BXHEIGHT(I,J,L) * 1d2 
               PTEMP        = PDOWN(L,I,J)

               ! PDOWN cannot be negative
               IF ( PDOWN(L,I,J) < 0.0d0 ) PDOWN(L,I,J) = 0.d0
            ENDDO
         ENDIF
      ENDDO  ! J
      ENDDO  ! I

      ! Return to calling program
      END SUBROUTINE MAKE_QQ

!------------------------------------------------------------------------------

      FUNCTION E_ICE( TC ) 
!
!******************************************************************************
!  Subroutine E_ICE computes Eice(T), the saturation vapor pressure of ice
!  at a given Celsius temperature.  (bmy, djj, 2/10/00, 11/8/02)
!  
!  Arguments as Input:
!  ===========================================================================
!  (1 ) TC (REAL*8) : Ambient temperature [C]
!
!  NOTES:
!  (1 ) The saturation vapor pressure of ice is generated by a polynomial
!        formulation from -50 C < T < 0 C.  From -120 C to -50 C, log10
!        interpolations are done in 10 C bins.  This should be good enough
!        for most purposes. (bmy, djj, 2/10/00) 
!  (2 ) Use the "D" exponent (e.g. 1.0d+00) to explicitly force double
!        precision for all floating point numbers (bmy, 2/10/00)
!  (3 ) Eliminate double declaration of E_ICE as REAL*8 (bmy, 7/16/01)
!  (4 ) Now references GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run.  Also updated comments
!        and made cosmetic changes.  Now bundled into "wetscav_mod.f". 
!        (bmy, 11/8/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
     
      ! Arguments 
      REAL*8, INTENT(IN) :: TC

      ! Local variables
      REAL*8             :: E_ICE

      ! Parameters for the polynomial (Lowe & Ficke, 1974)
      REAL*8, PARAMETER  :: A0 = 6.109177956d+00
      REAL*8, PARAMETER  :: A1 = 5.034698970d-01
      REAL*8, PARAMETER  :: A2 = 1.886013408d-02
      REAL*8, PARAMETER  :: A3 = 4.176223716d-04
      REAL*8, PARAMETER  :: A4 = 5.824720280d-06
      REAL*8, PARAMETER  :: A5 = 4.838803174d-08
      REAL*8, PARAMETER  :: A6 = 1.838826904d-10

      !=================================================================
      ! E_ICE begins here!
      !
      ! The polynomial formulation for Eice(T), -50 C <= T <= 0 C, 
      ! is as follows:
      !
      !   Eice(T) = a0 + T(a1 + T(a2 + T(a3 + T(a4 + T(a5 + a6T)))))
      !
      ! with coefficients a0...a6 given above by Lowe & Ficke, 1974.
      !
      ! The log10 interpolations for temperatures below -50 C are 
      ! done in 10 degree C bins.  The equation is as follows:
      !
      !  T0  = Temperature in C at lower bound of temperature bin
      !  T1  = Temperature in C at upper bound of temperature bin
      !  T   = Temperature for which Eice(T) is desired
      !  E0  = Eice(T0) = value of Eice at T0 (obtained from literature)
      !  E1  = Eice(T1) = value of Eice at T1 (obtained from literature)
      !  E   = Eice(T)  = value of Eice at desired temperature T
      !
      !                         {  T - T0                              }
      ! log10(E) =  log10(E0) + { -------- * [ log10(E1) - log10(E0) ] }
      !                         {    10                                }
      !
      !  The relevant values of T0, T1, E0, E1 for each bin are: 
      !
      !  Bin      T0       T1      E0 = Eice(T0)    E1 = Eice(T1)     
      ! ----------------------------------------------------------------
      !   1     -60 C    -50 C     1.0765e-2 mb     3.947e-2  mb
      !   2     -70 C    -60 C     2.577e-3  mb     1.0765e-2 mb
      !   3     -80 C    -70 C     5.307e-4  mb     2.577e-3  mb
      !   4     -90 C    -80 C     9.17e-5   mb     5.307e-4  mb
      !   5    -100 C    -90 C     1.29e-5   mb     9.17e-5   mb
      !   6    -110 C   -100 C     1.4e-6    mb     1.29e-5   mb
      !   7    -120 C   -110 C     1.0e-7    mb     1.4e-6    mb
      !
      ! For computational expediency, the following terms were also 
      ! computed for use in the code below:
      !  
      !  Bin           log10(E0)         [ log10(E1) - log10(E0) ]  
      ! ----------------------------------------------------------------
      !   1        -1.96798596584            5.6425309224e-01            
      !   2        -2.58888558145            6.2089961561e-01            
      !   3        -3.27515091237            6.8626533092e-01 
      !   4        -4.03763066433            7.6247975196e-01  
      !   5        -4.88941028970            8.5177962537e-01 
      !   6        -5.85387196432            9.6446167462e-01
      !   7        -7.00000000000            1.14612803568e+00
      !=================================================================

      ! Use a polynomial formulation for the saturation vapor 
      ! pressure over ice from -50 C to 0 C (Lowe & Ficke, 1974)
      IF ( TC >= -50d0 .and. TC <= 0d0 ) THEN
         E_ICE = A0 + TC * ( A1 + TC * ( A2 + TC * ( A3 + 
     &                TC * ( A4 + TC * ( A5 + A6 * TC ) ) ) ) )

      ! Bin 1 -- Log10 interpolation 
      ELSE IF ( TC >= -60d0 .and. TC < -50d0 ) THEN
         E_ICE = -1.96798596584d0 + 
     &           ( ( ( TC + 60d0 ) / 10d0 ) * 5.6425309224d-1 )

         E_ICE = 10d0**E_ICE

      ! Bin 2 -- Log10 interpolation
      ELSE IF ( TC >= -70d0 .and. TC < -60d0 ) THEN
         E_ICE = -2.58888558145d0 +
     &           ( ( ( TC + 70d0 ) / 10d0 ) * 6.2089961561d-1 ) 
                 
         E_ICE = 10d0**E_ICE

      ! Bin 3 -- Log10 interpolation
      ELSE IF ( TC >= -80d0 .and. TC < -70d0 ) THEN
         E_ICE = -3.27515091237d0 +
     &           ( ( ( TC + 80d0 ) / 10d0 ) * 6.8626533092d-1 ) 

         E_ICE = 10d0**E_ICE

      ! Bin 4 -- Log10 interpolation
      ELSE IF ( TC >= -90d0 .and. TC < -80d0 ) THEN
         E_ICE = -4.03763066433d0 +
     &           ( ( ( TC + 90d0 ) / 10d0 ) * 7.6247975196d-1 ) 

         E_ICE = 10d0**E_ICE

      ! Bin 5 -- Log10 interpolation
      ELSE IF ( TC >= -100d0 .and. TC < -90d0 ) THEN
         E_ICE =  -4.88941028970d0  + 
     &           ( ( ( TC + 100d0 ) / 10d0 ) * 8.5177962537d-1 ) 

         E_ICE = 10d0**E_ICE

      ! Bin 6 -- Log10 interpolation
      ELSE IF ( TC >= -110d0 .and. TC < -100d0 ) THEN
         E_ICE =  -5.85387196432d0 + 
     &           ( ( ( TC + 110d0 ) / 10d0 ) * 9.6446167462d-1 ) 

         E_ICE = 10d0**E_ICE

      ! Bin 7 -- Log10 interpolation
      ELSE IF ( TC >= -120d0 .and. TC < -110d0 ) THEN
         E_ICE = -7d0 + ( ( ( TC + 120d0 ) / 10d0 ) * 1.14612803568d0 ) 

         E_ICE = 10d0**E_ICE

      ! Below -120 C stop with an error message
      ELSE
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '( ''Temperature: '', f10.5, '' C'' )' ) TC
         WRITE( 6, '( ''T must be in the range -120 C <= T <= 0C'' )' )
         WRITE( 6, '( ''STOP in e_ice.f'' )' )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END FUNCTION E_ICE

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_L2G( Kstar298, H298_R, TK, H2OLIQ, L2G )
!
!******************************************************************************
!  Subroutine COMPUTE_L2G computes the ratio L2G = Cliq / Cgas, which is 
!  the mixing ratio of tracer in the liquid phase, divided by the mixing 
!  ratio of tracer in the gas phase.  (bmy, 2/23/00, 11/8/02)
!
!  The ratio Cliq / Cgas is obtained via Henry's law.  The appropriate 
!  values of Kstar298 and H298_R must be supplied for each tracer.  
!  (cf Jacob et al 2000, p. 3)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) Kstar298 (REAL*8) : Eff. Henry's law constant @ 298 K [moles/atm]
!  (2 ) H298_R   (REAL*8) : Henry's law coefficient           [K]
!  (3 ) TK       (REAL*8) : Temperature at grid box (I,J,L)   [K]
!  (4 ) H2OLIQ   (REAL*8) : Liquid water content at (I,J,L)   [cm3 H2O/cm3 air]
!
!  Arguments as Output:
!  ============================================================================
!  (5 ) L2G      (REAL*8) : Cliq/Cgas ratio for given tracer  [unitless]
!
!  References (see above for full citations):
!  ===========================================================================
!  (1 ) Jacob et al, 2000
!
!  NOTES:
!  (1 ) Bundled into "wetscav_mod.f" (bmy, 11/8/02)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: KStar298, H298_R, TK, H2OLIQ
      REAL*8, INTENT(OUT) :: L2G
      
      ! Local variables
      REAL*8              :: Kstar

      ! R = universal gas constant [atm/moles/K]
      REAL*8, PARAMETER   :: R = 8.32d-2

      ! INV_T0 = 1/298 K
      REAL*8, PARAMETER   :: INV_T0 = 1d0 / 298d0

      !=================================================================
      ! COMPUTE_L2G begins here!
      !=================================================================

      ! Get Kstar, the effective Henry's law constant for temperature TK
      Kstar = Kstar298 * EXP( -H298_R * ( ( 1d0 / TK ) - INV_T0 ) )

      ! Use Henry's Law to get the ratio:
      ! [ mixing ratio in liquid phase / mixing ratio in gas phase ]
      L2G   = Kstar * H2OLIQ * R * TK

      ! Return to calling program
      END SUBROUTINE COMPUTE_L2G

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_F( N, F, ISOL )
!
!******************************************************************************
!  Subroutine COMPUTE_F computes F, the fraction of soluble tracer lost by 
!  scavenging in convective cloud updrafts. (hyl, bmy, djj, 2/23/00, 4/5/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N    (INTEGER) : Tracer number
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) F    (REAL*8)  : Fraction of tracer scavenged in cloud updraft [0-1]
!  (3 ) ISOL (INTEGER) : Index number for ND38 diagnostic                     
!
!  References (see above for full citations):
!  ===========================================================================
!  (1 ) Jacob et al, 2000
!  (2 ) Chin et al, 1996
!
!  NOTES:
!  (1 ) Currently works computes scavenging fractions for either full
!        chemistry simulation (NSRCX == 3) or Rn-Pb-Be chemistry simulation
!        (NSRCX == 1).  Set the scavenging fraction to zero for other
!        simulations which do not carry soluble tracers. (bmy, 3/2/00)
!  (2 ) Need to call INIT_SCAV to initialize the Vud, C_H2O, CLDLIQ, 
!        and CLDICE fields once per timestep. (bmy, 2/23/00)
!  (3 ) For aerosols only: now apply Eq. 2 for all temperatures.  Also
!        use the distance between the grid box centers in Eq. 2.  Updated
!        comments and made some cosmetic changes (hyl, bmy, 6/18/01)
!  (4 ) Remove IREF, JREF -- these are obsolete.  T is now dimensioned
!        (IIPAR,JJPAR,LLPAR).  T(IREF,JREF,L) is now T(I,J,L). (bmy, 9/27/01)
!  (5 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (6 ) Fix 2 bugs for aerosol scavenging in Rn-Pb-Be simulation: 
!        (a) set F(:,:,1) = 0 since we don't do any scavenging there.  
!        (b) DO L = 2, LLPAR to avoid any subscript range out of bounds 
!        errors (rjp, hyl, bmy, 1/10/02)
!  (7 ) Now set F=0 in the first level for all tracers.  Also now
!        compute the distance between grid box centers and use that in 
!        in Eq. 10 from Jacob et al, 2000 to compute F. (hyl, bmy, 1/24/02)
!  (8 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (9 ) Now reference T from "dao_mod.f" instead of from "CMN".  Also reference
!        BXHEIGHT from "dao_mod.f" instead of from "CMN_NOX".  Now bundled
!        into "wetscav_mod.f".  Now references IDTHNO3, IDTH2O2, etc, from
!        F90 module "tracerid_mod.f".  Added internal routines F_AEROSOL
!        and GET_ISOL.  Rewritten so that we don't duplicate code for 
!        different chemistry simulations. (bmy, 1/17/03)
!  (10) Now compute F for SO2 in the same way for both fullchem and offline 
!        simulations (rjp, bmy, 3/23/03)
!  (11) Added slots for carbon aerosol & dust tracers.  Now modified internal
!        routine GET_ISOL so it's not hardwired anymore. (rjp, bmy, 4/5/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,     ONLY : BXHEIGHT, T
      USE TRACERID_MOD

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! NSRCX

      ! Arguments
      INTEGER, INTENT(IN)  :: N
      REAL*8,  INTENT(OUT) :: F(IIPAR,JJPAR,LLPAR)
      INTEGER, INTENT(OUT) :: ISOL

      ! Local variables 
      INTEGER              :: I, J, L, NN
      REAL*8               :: L2G, I2G, C_TOT, F_L, F_I, K, TMP, SO2LOSS

      ! Kc is the conversion rate from cloud condensate to precip [s^-1]
      REAL*8, PARAMETER    :: KC   = 5d-3

      ! CONV = 0.6 * SQRT( 1.9 ), used for the ice to gas ratio for H2O2
      REAL*8, PARAMETER    :: CONV = 8.27042925126d-1

      !=================================================================
      ! COMPUTE_F begins here!
      !
      ! For aerosol tracers, compute F with internal routine F_AEROSOL.
      !
      ! ISOL = tracer index for the ND38 diagnostic.  Values are:
      !
      ! Tracer   Rn-Pb-Be run   Fullchem run   Offline sulfate run
      ! ------   ------------   ------------   -------------------
      !  210Pb         1              -                -
      !  7Be           2              -                -
      !  HNO3          -              1                -
      !  H2O2          -              2                7 
      !  CH2O          -              3                -
      !  MP            -              4                -
      !  SO2           -              5                1
      !  SO4           -              6                2
      !  MSA           -              7                3
      !  NH3           -              8                4
      !  NH4           -              9                5
      !  NIT           -             10                6
      !=================================================================

      !----------------------------
      ! 210Pb (aerosol)
      !----------------------------
      IF ( N == IDTPB ) THEN
         CALL F_AEROSOL( F )
         !ISOL = 1
         ISOL = GET_ISOL( N )

      !----------------------------        
      ! 7Be (aerosol)
      !----------------------------
      ELSE IF ( N == IDTBE7 ) THEN
         CALL F_AEROSOL( F )
         !ISOL = 2
         ISOL = GET_ISOL( N )

      !----------------------------
      ! HNO3 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTHNO3 ) THEN
         CALL F_AEROSOL( F )
         !ISOL = 1
         ISOL = GET_ISOL( N ) 

      !----------------------------
      ! H2O2 (liquid & ice phases)
      !----------------------------
      ELSE IF ( N == IDTH2O2 ) THEN

         ! No scavenging at the surface
         F(:,:,1) = 0d0

         ! Start scavenging at level 2
         DO L = 2, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Compute ice to gas ratio for H2O2 by co-condensation
            ! (Eq. 9, Jacob et al, 2000)
            IF ( C_H2O(I,J,L) > 0d0 ) THEN 
               I2G = ( CLDICE(I,J,L) / C_H2O(I,J,L) ) * CONV
            ELSE
               I2G = 0d0
            ENDIF

            ! Compute liquid to gas ratio for H2O2, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            CALL COMPUTE_L2G( 8.3d4,    -7.4d3, 
     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )

            ! Fraction of H2O2 in liquid & ice phases
            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1d0 + L2G + I2G
            F_L   = L2G / C_TOT
            F_I   = I2G / C_TOT

            ! Compute the rate constant K.  The retention factor for 
            ! liquid H2O2 is 0.05 for 248 K < T < 268 K and 1.0 for 
            ! T >= 268 K. (Eq. 1, Jacob et al, 2000)
            IF ( T(I,J,L) >= 268d0 ) THEN
               K = KC * ( F_L + F_I )
 
            ELSE IF ( T(I,J,L) > 248d0  .and. T(I,J,L) < 268d0 ) THEN
               K = KC * ( ( 5d-2 * F_L ) + F_I ) 

            ELSE
               K = KC * F_I
                  
            ENDIF

            ! Distance between grid box centers [m]
            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) ) 
            
            ! Compute F, the fraction of scavenged H2O2.
            ! (Eq. 2, Jacob et al, 2000)
            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )

         ENDDO
         ENDDO
         ENDDO
            
         ! ND38 index
         !ISOL = GET_ISOL( 2, 7 )
         ISOL = GET_ISOL( N )

      !----------------------------
      ! CH2O (liquid phase only)
      !----------------------------
      ELSE IF ( N == IDTCH2O ) THEN 

         ! No scavenging at the surface
         F(:,:,1) = 0d0

         ! Apply scavenging in levels 2 and higher
         DO L = 2, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Compute liquid to gas ratio for CH2O, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            CALL COMPUTE_L2G( 3.0d3,    -7.2d3, 
     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )

            ! Fraction of CH2O in liquid phase 
            ! NOTE: CH2O does not exist in the ice phase!
            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1d0 + L2G
            F_L   = L2G / C_TOT

            ! Compute the rate constant K.  The retention factor 
            ! for liquid CH2O is 0.0 for T <= 248K and 0.02 for 
            ! 248 K < T < 268 K. (Eq. 1, Jacob et al, 2000)
            IF ( T(I,J,L) >= 268d0 ) THEN
               K = KC * F_L  

            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
               K = KC * ( 2d-2 * F_L ) 

            ELSE
               K = 0d0

            ENDIF

            ! Distance between grid box centers [m]
            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) ) 
               
            ! F is the fraction of CH2O scavenged out of the updraft
            ! (Eq. 2, Jacob et al, 2000)
            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )

         ENDDO
         ENDDO
         ENDDO

         ! ND38 index
         !ISOL = 3
         ISOL = GET_ISOL( N ) 
  
      !----------------------------
      ! CH3OOH (liquid phase only)
      !----------------------------
      ELSE IF ( N == IDTMP ) THEN

         ! No scavenging at the surface
         F(:,:,1) = 0d0

         ! Start scavenging at level 2
         DO L = 2, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Compute liquid to gas ratio for CH3OOH, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            CALL COMPUTE_L2G( 3.1d2,    -5.2d3, 
     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )

            ! Fraction of CH3OOH in liquid phase
            ! NOTE: CH3OOH does not exist in the ice phase!
            ! (Eq. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1d0 + L2G
            F_L   = L2G / C_TOT

            ! Compute the rate constant K.  The retention factor  
            ! for liquid CH3OOH is 0.0 for T <= 248 K and 0.02 for 
            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
            IF ( T(I,J,L) >= 268d0 ) THEN
               K = KC * F_L  

            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
               K = KC * ( 2d-2 * F_L ) 
                  
            ELSE
               K = 0d0

            ENDIF
               
            ! Distance between grid box centers [m]
            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) ) 

            ! F is the fraction of CH3OOH scavenged out of the updraft
            ! (Eq. 2, Jacob et al, 2000)
            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )

         ENDDO
         ENDDO
         ENDDO

         ! ND38 index
         !ISOL = 4
         ISOL = GET_ISOL( N )

      !---------------------------
      ! SO2 (aerosol)
      !---------------------------
      ELSE IF ( N == IDTSO2 ) THEN

         ! Compute fraction of SO2 scavenged
         CALL F_AEROSOL( F )
         !ISOL = GET_ISOL( 5, 1 )
         ISOL = GET_ISOL( N )

         !==============================================================
         ! Coupled full chemistry/aerosol simulation:
         ! Use the wet scavenging formula of Chin et al [1996], 
         ! such that a soluble fraction of SO2 is limited by the
         ! availability of H2O2 in the precipitating grid box. 
         ! Scavenge the soluble SO2 at the same rate as the sulfate.
         ! Update H2O2_sav and SO2_sav for use in RAINOUT, WASHOUT
         !==============================================================
         DO L = 2, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
                
            IF ( ( H2O2s(I,J,L) < SO2s(I,J,L) ) .AND. 
     &           ( SO2s(I,J,L)  > EPSILON     ) ) THEN

               ! Compute F
               F(I,J,L)     = F(I,J,L) * ( H2O2s(I,J,L)/SO2s(I,J,L) )
               H2O2s(I,J,L) = 0d0

            ELSE
                  
               ! Update H2O2
               H2O2s(I,J,L) = H2O2s(I,J,L) - SO2s(I,J,L)
               H2O2s(I,J,L) = MAX( H2O2s(I,J,L), EPSILON )

            ENDIF

            ! Update SO2
            SO2s(I,J,L)     = SO2s(I,J,L) * ( 1d0 - F(I,J,L) )
            SO2s(I,J,L)     = MAX( SO2s(I,J,L), EPSILON )
              
         ENDDO
         ENDDO
         ENDDO

      !----------------------------
      ! SO4 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTSO4 ) THEN
         CALL F_AEROSOL( F ) 
         !ISOL = GET_ISOL( 6, 2 )
         ISOL = GET_ISOL( N ) 

      !---------------------------
      ! MSA (aerosol)
      !---------------------------
      ELSE IF ( N == IDTMSA ) THEN
         CALL F_AEROSOL( F )
         !ISOL = GET_ISOL( 7, 3 )
         ISOL = GET_ISOL( N )

      !----------------------------
      ! NH3 (liquid & ice phases)
      !----------------------------
      ELSE IF ( N == IDTNH3 ) THEN
         
         ! No scavenging at surface
         F(:,:,1) = 0d0

         ! Start scavenging at level 2
         DO L = 2, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Compute ice to gas ratio for NH3 by co-condensation
            ! (Eq. 9, Jacob et al, 2000)
            IF ( C_H2O(I,J,L) > 0d0 ) THEN 
               I2G = ( CLDICE(I,J,L) / C_H2O(I,J,L) ) * CONV
            ELSE
               I2G = 0d0
            ENDIF

            ! Compute liquid to gas ratio for NH3, using
            ! the appropriate parameters for Henry's law
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            CALL COMPUTE_L2G( 3.3d6,    -4.1d3, 
     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )

            ! Fraction of NH3 in liquid & ice phases
            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
            C_TOT = 1d0 + L2G + I2G
            F_L   = L2G / C_TOT
            F_I   = I2G / C_TOT

            ! Compute the rate constant K.  The retention factor  
            ! for liquid NH3 is 0.0 for T <= 248 K and 0.05 for 
            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
            IF ( T(I,J,L) >= 268d0 ) THEN
               K = KC * ( F_L + F_I )
                  
            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
               K = KC * ( ( 5d-2 * F_L ) + F_I )
                  
            ELSE
               K = KC * F_I

            ENDIF
  
            ! F is the fraction of NH3 scavenged out of the updraft
            ! (Eq. 2, Jacob et al, 2000)
            F(I,J,L) = 1d0 - EXP( -K * BXHEIGHT(I,J,L) / Vud(I,J) )
         ENDDO
         ENDDO
         ENDDO
       
         ! ND38 index
         !ISOL = GET_ISOL( 8, 4 )
         ISOL = GET_ISOL( N )

      !----------------------------
      ! NH4 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTNH4 ) THEN
         CALL F_AEROSOL( F )
         !ISOL = GET_ISOL( 9, 5 )
         ISOL = GET_ISOL( N )
         
      !----------------------------
      ! NIT (aerosol)
      !----------------------------
      ELSE IF ( N == IDTNIT ) THEN
         CALL F_AEROSOL( F )
         !ISOL = GET_ISOL( 10, 6 )
         ISOL = GET_ISOL( N )

      !----------------------------
      ! BC HYDROPHILIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTBCPI ) THEN
         CALL F_AEROSOL( F )
         ISOL = GET_ISOL( N )

      !----------------------------
      ! OC HYDROPHILIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTOCPI ) THEN
         CALL F_AEROSOL( F )
         ISOL = GET_ISOL( N )

      !----------------------------
      ! BC HYDROPHOBIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTBCPO ) THEN

         ! Force not to be lost in convective updraft for now
         F    = 0d0
         ISOL = GET_ISOL( N )

      !----------------------------
      ! OC HYDROPHOBIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTOCPO ) THEN

         ! Force not to be lost in convective updraft for now
         F    = 0d0
         ISOL = GET_ISOL( N )

      !----------------------------
      ! DUST1 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST1 ) THEN
         CALL F_AEROSOL( F )
         ISOL = GET_ISOL( N )
         
      !----------------------------
      ! DUST2 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST2 ) THEN
         CALL F_AEROSOL( F )
         ISOL = GET_ISOL( N )
         
      !----------------------------
      ! DUST3 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST3 ) THEN
         CALL F_AEROSOL( F )
         ISOL = GET_ISOL( N )
         
      !----------------------------
      ! DUST4 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST4 ) THEN
         CALL F_AEROSOL( F )
         ISOL = GET_ISOL( N )

      !----------------------------
      ! Insoluble tracer, set F=0
      !----------------------------
      ELSE
         F(:,:,:) = 0d0
         ISOL     = 0

      ENDIF
         
      !=================================================================
      ! Internal subroutines -- they can "see" all variables 
      ! that have been declared in COMPUTE_F
      !=================================================================
      CONTAINS

      !-----------------------------------------------------------------------

      SUBROUTINE F_AEROSOL( F ) 

      !=================================================================
      ! Internal routine F_AEROSOL returns the fraction of aerosol
      ! scavenged in updrafts (bmy, 11/7/02)
      !=================================================================

      ! Arguments
      REAL*8, INTENT(OUT) :: F(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER             :: I, J, L
      REAL*8              :: TMP
     
      !=================================================================
      ! GET_F_AEROSOL begins here!
      !
      ! Aerosol tracers are 100% in the cloud condensate phase, so 
      ! we set K = Kc, and compute F accordingly (cf Jacob et al 2000 )    
      !=================================================================
      
      ! Turn off scavenging in the first level by setting F = 0
      F(:,:,1) = 0d0
  
      ! Apply scavenging in levels 2 and higher
      DO L = 2, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
                
         ! Distance between grid box centers [m]
         TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) ) 
       
         ! (Eq. 2, Jacob et al, 2000, with K = Kc)
         F(I,J,L) = 1d0 - EXP( -KC * TMP / Vud(I,J) )
            
      ENDDO
      ENDDO
      ENDDO

      ! Return to COMPUTE_F
      END SUBROUTINE F_AEROSOL 

!----------------------------------------------------------------------------
! Prior to 4/5/04:    
! Now make GET_ISOL a module routine (bmy, 4/5/04)
!      FUNCTION GET_ISOL( I3, I10 ) RESULT( I )
!       
!      !=================================================================
!      ! Internal function GET_ISOL returns the value of ISOL (tracer
!      ! index for ND38) for either coupled or offline sulfate runs.
!      !=================================================================
!       
!      ! Arguments
!      INTEGER, INTENT(IN) :: I3, I10
!       
!      ! Function value
!      INTEGER             :: I
!       
!      !=================================================================
!      ! GET_ISOL begins here!
!      !=================================================================
!      SELECT CASE ( NSRCX )
!         CASE ( 3 ) 
!            I = I3
!         CASE ( 10 )
!            I = I10 
!      END SELECT
!       
!      ! Return to COMPUTE_F
!      END FUNCTION GET_ISOL
!------------------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE COMPUTE_F

!------------------------------------------------------------------------------

      FUNCTION GET_ISOL( N_TEST ) RESULT( VALUE )
!
!******************************************************************************
!  Function GET_ISOL returns the value of ISOL (tracer index for ND38) for 
!  all simulation types.  (bmy, 4/5/04)
!
!
! NOTES:
! (1 ) Now initializes a lookup table for faster execution.  Now made into
!        an EXTERNAL function. (rjp, bmy, 4/5/04)
!******************************************************************************
!
#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! NTRACE

      ! Arguments
      INTEGER, INTENT(IN) :: N_TEST
      
      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER, SAVE       :: NSOL_INDEX(NNPAR)
      INTEGER             :: I, L, N

      ! Function value
      INTEGER             :: VALUE
 
      !=================================================================
      ! GET_ISOL begins here!
      !=================================================================

      ! Initialize lookup table on the first call
      IF ( FIRST ) THEN
      
         ! Initialize
         NSOL_INDEX(:) = 0
   
         ! Loop over tracers
         DO N = 1, NTRACE
            
            ! Loop over soluble tracers
            DO L = 1, NSOL

               ! Test if tracer N is among the soluble tracers
               IF ( IDWETD(L) == N ) THEN

                  ! Save location into the lookup table
                  NSOL_INDEX(N) = L 

                  ! Go to next N
                  GOTO 100
               ENDIF
            ENDDO

 100        CONTINUE
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return value
      VALUE = NSOL_INDEX(N_TEST)

      ! Return to COMPUTE_F
      END FUNCTION GET_ISOL

!------------------------------------------------------------------------------

      SUBROUTINE RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )
!
!******************************************************************************
!  Subroutine RAINOUT computes RAINFRAC, the fraction of soluble tracer
!  lost to rainout events in precipitation. (djj, bmy, 2/28/00, 4/5/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : Grid box lon-lat-alt indices
!  (4  ) N        (INTEGER) : Tracer number
!  (5  ) K_RAIN   (REAL*8 ) : Rainout rate constant for tracer N [s^-1]
!  (6  ) DT       (REAL*8 ) : Timestep for rainout event         [s]
!  (7  ) F        (REAL*8 ) : Fraction of grid box precipitating [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (8  ) RAINFRAC (REAL*8)  : Fraction of tracer lost to rainout [unitless]
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Jacob et al, 2000
!  (2 ) Chin et al, 1996
!
!  NOTES:
!  (1 ) Currently works for either full chemistry simulation (NSRCX == 3) 
!        or Rn-Pb-Be chemistry simulation (NSRCX == 1).  Other simulations
!        do not carry soluble tracer, so set RAINFRAC = 0. (bmy, 2/28/00)
!  (2 ) Need to call INIT_SCAV to initialize the Vud, C_H2O, CLDLIQ, 
!        and CLDICE fields once per dynamic timestep. (bmy, 2/28/00)
!  (3 ) K_RAIN, the rainout rate constant, and F, the areal fraction of the 
!        grid box undergoing precipitiation, are computed according to 
!        Giorgi & Chaimedes, as described in Jacob et al, 2000.
!  (4 ) Now no longer suppress scavenging of HNO3 and aerosol below 258K.
!        Updated comments, cosmetic changes.  Now set TK = T(I,J,L) since
!        T is now sized (IIPAR,JJPAR,LLPAR) in "CMN". (djj, hyl, bmy, 1/24/02)
!  (5 ) Eliminated obsolete code (bmy, 2/27/02)
!  (6 ) Now reference T from "dao_mod.f".  Updated comments.  Now bundled 
!        into "wetscav_mod.f". Now refererences "tracerid_mod.f".  Also 
!        removed reference to CMN since we don't need NSRCX. (bmy, 11/8/02)
!  (7 ) Now updated for carbon & dust aerosol tracers (rjp, bmy, 4/5/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,     ONLY : T
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE TRACERID_MOD

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L, N
      REAL*8,  INTENT(IN)  :: K_RAIN, DT, F
      REAL*8,  INTENT(OUT) :: RAINFRAC

      ! Local variables 
      REAL*8               :: L2G, I2G, C_TOT, F_L, F_I, K, TK, SO2LOSS

      ! CONV = 0.6 * SQRT( 1.9 ), used for the ice to gas ratio for H2O2
      REAL*8, PARAMETER    :: CONV = 8.27042925126d-1

      !=================================================================
      ! RAINOUT begins here!
      !
      ! For aerosols, set K = K_RAIN and compute RAINFRAC according
      ! to Eq. 10 of Jacob et al 2000.  Call function GET_RAINFRAC.
      !=================================================================

      ! Save the local temperature in TK for convenience
      TK = T(I,J,L)

      !----------------------------
      ! 210Pb (aerosol)
      !----------------------------
      IF ( N == IDTPB ) THEN 
         RAINFRAC = GET_RAINFRAC( K_RAIN )
         
      !----------------------------
      ! 7Be (aerosol)
      !----------------------------
      ELSE IF ( N == IDTBE7 ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )
         
      !----------------------------
      ! HNO3 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTHNO3 ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !----------------------------
      ! H2O2 (liquid & ice phases)
      !----------------------------
      ELSE IF ( N == IDTH2O2 ) THEN

         ! Compute ice to gas ratio for H2O2 by co-condensation
         ! (Eq. 9, Jacob et al, 2000)
         IF ( C_H2O(I,J,L) > 0d0 ) THEN 
            I2G = ( CLDICE(I,J,L) / C_H2O(I,J,L) ) * CONV
         ELSE
            I2G = 0d0
         ENDIF

         ! Compute liquid to gas ratio for H2O2, using
         ! the appropriate parameters for Henry's law
         ! (Eqs. 7, 8 and Table 1, Jacob et al, 2000)
         CALL COMPUTE_L2G( 8.3d4, -7.4d3, TK, CLDLIQ(I,J,L), L2G )

         ! Fraction of H2O2 in liquid & ice phases
         ! (Eqs. 4, 5, 6, Jacob et al, 2000)
         C_TOT = 1d0 + L2G + I2G
         F_L   = L2G / C_TOT
         F_I   = I2G / C_TOT

         ! Compute the rate constant K.  The retention factor  
         ! for liquid H2O2 is 0.05 for 248 K < T < 268 K, and 
         ! 1.0 for T >= 268 K.  (Eq. 1, Jacob et al, 2000)
         IF ( TK >= 268d0 ) THEN
            K = K_RAIN * ( F_L + F_I ) 
               
         ELSE IF ( TK > 248d0  .and. TK < 268d0 ) THEN
            K = K_RAIN * ( ( 5d-2 * F_L ) + F_I ) 

         ELSE 
            K = K_RAIN * F_I

         ENDIF

         ! Compute RAINFRAC, the fraction of rained-out H2O2
         ! (Eq. 10, Jacob et al, 2000)
         RAINFRAC = GET_RAINFRAC( K )

      !----------------------------
      ! CH2O (liquid phase only)
      !----------------------------      
      ELSE IF ( N == IDTCH2O ) THEN 
                  
         ! Compute liquid to gas ratio for CH2O, using
         ! the appropriate parameters for Henry's law
         ! (Eqs. 7, 8 and Table 1, Jacob et al, 2000)
         CALL COMPUTE_L2G( 3.0d3, -7.2d3, TK, CLDLIQ(I,J,L), L2G )
            
         ! Fraction of CH2O in liquid phase 
         ! NOTE: CH2O does not exist in the ice phase!
         ! (Eqs. 4, 5, Jacob et al, 2000)
         C_TOT = 1d0 + L2G
         F_L   = L2G / C_TOT

         ! Compute the rate constant K.  The retention factor  
         ! for liquid CH2O is 0.02 for 248 K < T < 268 K, and 
         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
         IF ( TK >= 268d0 ) THEN
            K = K_RAIN * F_L 

         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
            K = K_RAIN * ( 2d-2 * F_L )
               
         ELSE
            K = 0d0

         ENDIF

         ! Compute RAINFRAC, the fraction of rained-out CH2O
         ! (Eq. 10, Jacob et al, 2000)
         RAINFRAC = GET_RAINFRAC( K )
  
      !----------------------------
      ! CH3OOH (liquid phase only)
      !----------------------------
      ELSE IF ( N == IDTMP ) THEN

         ! Compute liquid to gas ratio for CH3OOH, using
         ! the appropriate parameters for Henry's law
         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
         CALL COMPUTE_L2G( 3.1d2, -5.2d3, TK, CLDLIQ(I,J,L), L2G )

         ! Fraction of CH3OOH in liquid phase
         ! NOTE: CH3OOH does not exist in the ice phase!
         ! (Eqs. 4, 5, Jacob et al, 2000)
         C_TOT = 1d0 + L2G
         F_L   = L2G / C_TOT

         ! Compute the rate constant K.  The retention factor  
         ! for liquid CH3OOH is 0.02 for 248 K < T < 268 K, and 
         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
         IF ( TK >= 268d0 ) THEN
            K = K_RAIN * F_L  

         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
            K = K_RAIN * ( 2d-2 * F_L )  
            
         ELSE
            K = 0d0
               
         ENDIF
  
         ! Compute RAINFRAC, the fraction of rained-out CH3OOH
         ! (Eq. 10, Jacob et al, 2000)
         !RAINFRAC = F * ( 1 - EXP( -K * DT ) )
         RAINFRAC = GET_RAINFRAC( K ) 

      !----------------------------
      ! SO2
      !----------------------------
      ELSE IF ( N == IDTSO2 ) THEN

         !==============================================================
         ! NOTE: SO2 and H2O2 are in [v/v] and here RAINFRAC contains 
         ! the amount of SO2 lost due to rainout normalized by the
         ! total SO2 -- so that in WETDEP routine mulitiplying SO2 in 
         ! [kg] will produce correct amount.  Need to verify this. 
         ! (rjp, 01/16/02)
         !==============================================================

         ! Treat SO2 as an aerosol
         RAINFRAC = GET_RAINFRAC( K_RAIN )

         ! Update SO2 and H2O2
         IF ( SO2s(I,J,L) > EPSILON ) THEN 
         
            ! Limit RAINFRAC 
            SO2LOSS      = MIN( SO2s(I,J,L), H2O2s(I,J,L) )
            RAINFRAC     = SO2LOSS * RAINFRAC / SO2s(I,J,L)
            RAINFRAC     = MAX( RAINFRAC, 0d0 )
         
            ! Update saved H2O2 concentration
            H2O2s(I,J,L) = H2O2s(I,J,L) - ( SO2s(I,J,L) * RAINFRAC )
            H2O2s(I,J,L) = MAX( H2O2s(I,J,L), EPSILON )
         
         ELSE
            RAINFRAC = 0D0
         
         ENDIF
         
         ! Update saved SO2 concentration
         SO2s(I,J,L) = SO2s(I,J,L) * ( 1.D0 - RAINFRAC )
         SO2s(I,J,L) = MAX( SO2s(I,J,L), EPSILON )         

      !----------------------------
      ! SO4 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTSO4 ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !----------------------------
      ! MSA (aerosol)
      !----------------------------
      ELSE IF ( N == IDTMSA ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !----------------------------
      ! NH3 (liquid & ice phases)
      !----------------------------
      ELSE IF ( N == IDTNH3 ) THEN

         ! Compute ice to gas ratio for NH3 by co-condensation
         ! (Eq. 9, Jacob et al, 2000)
         IF ( C_H2O(I,J,L) > 0d0 ) THEN 
            I2G = ( CLDICE(I,J,L) / C_H2O(I,J,L) ) * CONV
         ELSE
            I2G = 0d0
         ENDIF
                  
         ! Compute liquid to gas ratio for NH3, using
         ! the appropriate parameters for Henry's law
         ! (Seinfeld and Pandis, p343 eq. 6.8)
         ! PH    = 4.5  ! Assumed PH for typical cloud drop
         ! Hstar = 1.054d11 * (10.**(-PH)) == 3.3d6
         CALL COMPUTE_L2G( 3.3d6, -4.1d3, TK, CLDLIQ(I,J,L), L2G )

         ! Fraction of NH3 in liquid & ice phases
         ! (Eqs. 4, 5, 6, Jacob et al, 2000)
         C_TOT = 1d0 + L2G + I2G
         F_L   = L2G / C_TOT
         F_I   = I2G / C_TOT

         ! Compute the rate constant K.  The retention factor  
         ! for liquid NH3 is 0.05 for 248 K < T < 268 K, and 
         ! 1.0 for T >= 268 K.  (Eq. 1, Jacob et al, 2000)
         IF ( TK >= 268d0 ) THEN
            K = K_RAIN * ( F_L + F_I ) 

         ELSE IF ( TK > 248d0  .and. TK < 268d0 ) THEN
            K = K_RAIN * ( ( 5d-2 * F_L ) + F_I ) 

         ELSE 
            K = K_RAIN * F_I

         ENDIF

         ! Compute RAINFRAC, the fraction of rained-out NH3
         ! (Eq. 10, Jacob et al, 2000)
         RAINFRAC = GET_RAINFRAC( K )

      !----------------------------
      ! NH4 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTNH4 ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !----------------------------
      ! NIT (aerosol)
      !----------------------------
      ELSE IF ( N == IDTNIT ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !----------------------------
      ! BC HYDROPHILIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTBCPI ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !----------------------------
      ! OC HYDROPHILIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTOCPI ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !-----------------------------
      ! BC HYDROPHOBIC (aerosol) 
      !-----------------------------
      ELSE IF ( N == IDTBCPO ) THEN

         ! No rainout 
         RAINFRAC = 0.0D0                  

      !-----------------------------
      ! OC HYDROPHOBIC (aerosol) 
      !-----------------------------
      ELSE IF ( N == IDTOCPO ) THEN

         ! No rainout
         RAINFRAC = 0.0D0 

      !----------------------------
      ! DUST1 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST1 ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !----------------------------
      ! DUST2 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST2 ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !----------------------------
      ! DUST3 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST3 ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )

      !----------------------------
      ! DUST4 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST4 ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN )      

      !----------------------------
      ! ERROR: insoluble tracer!
      !----------------------------  
      ELSE
         CALL ERROR_STOP( 'Invalid tracer!', 'RAINOUT (wetscav_mod.f)' )

      ENDIF
      
      !=================================================================
      ! Internal subroutines -- they can "see" all 
      ! variables that have been declared in RAINOUT
      !================================================================= 
      CONTAINS

      !------------------------------------------------------------------------

      FUNCTION GET_RAINFRAC( K ) RESULT( RAINFRAC )

      !=================================================================
      ! Internal function GET_RAINFRAC computes the fraction of tracer
      ! lost to rainout according to Jacob et al 2000. (bmy, 11/8/02)
      !=================================================================

      ! Arguments
      REAL*8, INTENT(IN) :: K
      
      ! Local variables
      REAL*8             :: RAINFRAC

      !=================================================================
      ! GET_RAINFRAC begins here!
      !=================================================================

      ! (Eq. 10, Jacob et al, 2000 ) 
      RAINFRAC = F * ( 1 - EXP( -K * DT ) )

      ! Return to RAINOUT
      END FUNCTION GET_RAINFRAC

      !------------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE RAINOUT

!------------------------------------------------------------------------------

      SUBROUTINE WASHOUT( I, J, L, N, PP, DT, F, WASHFRAC, AER )
!
!******************************************************************************
!  Subroutine WASHOUT computes WASHFRAC, the fraction of soluble tracer
!  lost to washout events in precipitation. (djj, bmy, 2/28/00, 4/5/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : Grid box lon-lat-alt indices [unitless]
!  (4  ) N        (INTEGER) : Tracer number                [unitless]
!  (5  ) PP       (REAL*8 ) : Precip rate thru at bottom    
!                             of grid box (I,J,L)          [cm3 H2O/cm2 air/s]
!  (6  ) DT       (REAL*8 ) : Timestep for rainout event   [s]
!  (7  ) F        (REAL*8 ) : Fraction of grid box 
!                             that is precipitating        [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (8  ) WASHFRAC (REAL*8)  : Fraction of tracer lost to rainout [unitless]
!  (9  ) AER      (LOGICAL) : = T if the tracer is an aerosol, =F otherwise
!
!  Reference (see above for full citations):
!  ============================================================================
!  (1  ) Jacob et al, 2000
!
!  NOTES:
!  (1 ) Currently works for either full chemistry simulation (NSRCX == 3) 
!        or Rn-Pb-Be chemistry simulation (NSRCX == 1).  Other simulations
!        do not carry soluble tracers, so set WASHFRAC = 0. 
!  (2 ) K_WASH, the rainout rate constant, and F, the areal fraction of the 
!        grid box undergoing precipitiation, are computed according to 
!        Giorgi & Chaimedes, as described in Jacob et al, 2000.
!  (3 ) Washout is only done for T >= 268 K, when the cloud condensate is
!        in the liquid phase. 
!  (4 ) T(I+I0,J+J0,L) is now T(I,J,L).  Removed IREF, JREF -- these are 
!        obsolete.  Updated comments. (bmy, 9/27/01)
!  (5 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (6 ) Now reference BXHEIGHT, T from "dao_mod.f".  Also remove reference
!        to "CMN_NOX".  Updated comments.  Now bundled into "wetscav_mod.f".
!        Now also references "tracerid_mod.f".  Added internal routines
!        WASHFRAC_AEROSOL and WASHFRAC_LIQ_GAS.  Also removed reference to
!        CMN since we don't need to use NSRCX here. (bmy, 11/6/02)
!  (7 ) Updated for carbon aerosol and dust tracers (rjp, bmy, 4/5/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : BXHEIGHT, T
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACERID_MOD

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L, N
      REAL*8,  INTENT(IN)  :: PP, DT, F
      REAL*8,  INTENT(OUT) :: WASHFRAC
      LOGICAL, INTENT(OUT) :: AER

      ! Local variables 
      REAL*8               :: L2G,  LP, DZ, TK, WASHFRAC_F_14, SO2LOSS

      ! First order washout rate constant for HNO3, aerosols = 1 cm^-1
      REAL*8, PARAMETER    :: K_WASH = 1d0

      !=================================================================
      ! WASHOUT begins here!
      !
      ! Call either WASHFRAC_AEROSOL or WASHFRAC_LIQ_GAS to compute the
      ! fraction of tracer lost to washout according to Jacob et al 2000
      !=================================================================

      ! TK is Kelvin temperature 
      TK = T(I,J,L)

      ! DZ is the height of the grid box in cm
      DZ = BXHEIGHT(I,J,L) * 1d2
      
      ! Initialize some variables
      L2G           = 0d0
      LP            = 0d0
      WASHFRAC_F_14 = 0d0

      !----------------------------
      ! 210Pb (aerosol)
      !----------------------------
      IF ( N == IDTPB ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! 7Be (aerosol)
      !----------------------------
      ELSE IF ( N == IDTBE7 ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! HNO3 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTHNO3 ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! H2O2 (liquid & gas phases)
      !----------------------------
      ELSE IF ( N == IDTH2O2 ) THEN
         AER      = .FALSE.
         WASHFRAC = WASHFRAC_LIQ_GAS( 8.3d4, -7.4d3 )

      !----------------------------
      ! CH2O (liquid & gas phases)
      !----------------------------
      ELSE IF ( N == IDTCH2O ) THEN 
         AER      = .FALSE.
         WASHFRAC = WASHFRAC_LIQ_GAS( 3.0d3, -7.2d3 )

      !----------------------------
      ! MP (liquid & gas phases)
      !----------------------------
      ELSE IF ( N == IDTMP ) THEN
         AER      = .FALSE.
         WASHFRAC = WASHFRAC_LIQ_GAS( 3.1d2, -5.2d3 )

      !----------------------------
      ! SO2 (aerosol treatment)
      !----------------------------
      ELSE IF ( N == IDTSO2 ) THEN

         !==============================================================
         ! NOTE: Even though SO2 is not an aerosol we treat it as SO4 in
         ! wet scavenging.  When evaporation occurs, it returns to SO4.
         !==============================================================
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

         !==============================================================
         ! Use the wet-scavenging following [Chin et al, 1996] such 
         ! that a soluble fraction of SO2 is limited by the availability 
         ! of H2O2 in the precipitating grid box.  Then scavenge the 
         ! soluble SO2 at the same rate as sulfate.
         !==============================================================
         IF ( TK >= 268d0 .AND. SO2s(I,J,L) > EPSILON ) THEN
         
            ! Adjust WASHFRAC
            SO2LOSS  = MIN( SO2s(I,J,L), H2O2s(I,J,L) )
            WASHFRAC = SO2LOSS * WASHFRAC / SO2s(I,J,L)
            WASHFRAC = MAX( WASHFRAC, 0d0 )
         
            ! Update saved H2O2 concentration
            IF ( SO2s(I,J,L) < H2O2s(I,J,L) ) THEN
               H2O2s(I,J,L) = H2O2s(I,J,L) - SO2s(I,J,L)
               H2O2s(I,J,L) = MAX( H2O2s(I,J,L), EPSILON )
            ELSE
               H2O2s(I,J,L) = 0.d0
            ENDIF
         
         ELSE
            WASHFRAC = 0d0
         
         ENDIF
         
         ! Update saved SO2 concentration 
         SO2s(I,J,L) = SO2s(I,J,L) * ( 1d0 - WASHFRAC )
         SO2s(I,J,L) = MAX( SO2s(I,J,L), EPSILON )
           
      !----------------------------
      ! SO4 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTSO4 ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! MSA (aerosol)
      !----------------------------
      ELSE IF ( N == IDTMSA ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! NH3 (liquid & gas phases)
      !----------------------------
      ELSE IF ( N == IDTNH3 ) THEN
         AER      = .FALSE.
         WASHFRAC = WASHFRAC_LIQ_GAS( 3.3d6, -4.1d3 )

      !----------------------------
      ! NH4 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTNH4 ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! NIT (aerosol)
      !----------------------------
      ELSE IF ( N == IDTNIT ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! BC HYDROPHILIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTBCPI ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! OC HYDROPHILIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTOCPI ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! BC HYDROPHOBIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTBCPO ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL() 

      !----------------------------
      ! OC HYDROPHOBIC (aerosol)
      !----------------------------
      ELSE IF ( N == IDTOCPO ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! DUST1 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST1 ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! DUST2 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST2 ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! DUST3 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST3 ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! DUST4 (aerosol)
      !----------------------------
      ELSE IF ( N == IDTDST4 ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL()

      !----------------------------
      ! ERROR: Insoluble tracer
      !----------------------------
      ELSE 
         CALL ERROR_STOP( 'Invalid tracer!', 'WASHOUT (wetscav_mod.f)' )

      ENDIF

      !=================================================================
      ! Internal subroutines -- they can "see" all 
      ! variables that have been declared in RAINOUT
      !================================================================= 
      CONTAINS

      !------------------------------------------------------------------------

      FUNCTION WASHFRAC_AEROSOL() RESULT( WASHFRAC )

      !=================================================================
      ! Internal function WASHFRAC_AEROSOL returns the fraction of
      ! soluble aerosol tracer lost to washout (bmy, 11/8/02)
      !=================================================================

      ! Function value
      REAL*8 :: WASHFRAC

      ! Washout only happens at or above 268 K
      IF ( TK >= 268d0 ) THEN
         WASHFRAC = F * ( 1d0 - EXP( -K_WASH * ( PP / F ) * DT ) )
      ELSE
         WASHFRAC = 0d0
      ENDIF

      ! Return to WASHFRAC
      END FUNCTION WASHFRAC_AEROSOL

      !------------------------------------------------------------------------

      FUNCTION WASHFRAC_LIQ_GAS( Kstar298, H298_R ) RESULT( WASHFRAC )

      !=================================================================
      ! Internal function WASHFRAC_LIQ_GAS returns the fraction of
      ! soluble liquid/gas phase tracer lost to washout (bmy, 11/8/02)
      !
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) Kstar298 : Eff. Henry's law constant @ 298 K [moles/atm]
      ! (2 ) H298_R   : Henry's law coefficient           [K]
      !=================================================================

      ! Arguments
      REAL*8, INTENT(IN) :: Kstar298, H298_R

      ! Local variables
      REAL*8             :: LP, WASHFRAC, WASHFRAC_F_14

      !=================================================================
      ! WASHFRAC_LIQ_GAS begins here!
      !=================================================================

      ! Suppress washout below 268 K
      IF ( TK >= 268d0 ) THEN

         ! Rainwater content in the grid box (Eq. 17, Jacob et al, 2000)
         LP = ( PP * DT ) / ( F * DZ ) 

         ! Compute liquid to gas ratio for H2O2, using the appropriate 
         ! parameters for Henry's law -- also use rainwater content Lp
         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
         CALL COMPUTE_L2G( Kstar298, H298_R, TK, LP, L2G )

         ! Washout fraction from Henry's law (Eq. 16, Jacob et al, 2000)
         WASHFRAC = L2G / ( 1d0 + L2G )

         ! Washout fraction / F from Eq. 14, Jacob et al, 2000
         WASHFRAC_F_14 = 1d0 - EXP( -K_WASH * ( PP / F ) * DT )

         ! Do not let the Henry's law washout fraction exceed
         ! ( washout fraction / F ) from Eq. 14 -- this is a cap
         IF ( WASHFRAC > WASHFRAC_F_14 ) WASHFRAC = WASHFRAC_F_14
            
      ELSE
         WASHFRAC = 0d0
            
      ENDIF

      ! Return to WASHOUT
      END FUNCTION WASHFRAC_LIQ_GAS

      !------------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE WASHOUT

!------------------------------------------------------------------------------

      SUBROUTINE WETDEP( LS )
!
!******************************************************************************
!  Subroutine WETDEP computes the downward mass flux of tracer due to washout 
!  and rainout of aerosols and soluble tracers in a column.  The timestep is 
!  the dynamic timestep. (hyl, bey, bmy, djj, 4/2/99, 3/17/04)
!
!  The precip fields through the bottom of each level are indexed as follows:
!
!       Layer          GISS-CTM II         GEOS-CTM
!
!      ------------------------------------------------- Top of Atm.
!        LM            PSSW4(I,J,LM-1)   PDOWN(LM,I,J)
!                          |                  |
!      ====================V==================V========= Max Extent 
!        LM-1          PSSW4(I,J,LM)     PDOWN(LM-1,I,J)  of Clouds
!                          |                  |
!      --------------------V------------------V---------
!                         ...                ...             
!
!      -------------------------------------------------
!        4             PSSW4(I,J,3)      PDOWN(4,I,J)
!                          |                  |
!      --------------------V------------------V----------
!        3             PSSW4(I,J,2)      PDOWN(3,I,J)
!                          |                  |
!      --------------------V------------------V--------- Cloud base
!        2             PSSW4(I,J,1)      PDOWN(2,I,J) 
!                          |                  |
!      -  -  -  -  -  -  - V -  -   -   -   - V -  -  - 
!        1                               PDOWN(1,I,J) 
!                                             |
!      =======================================V========= Ground
!
!  From the diagram, we have the following for layer L:
!     
!  GISS-CTM:
!  (a) Precip coming in  thru top    of layer L = PSSW4(I,J,L  )
!  (b) Precip going  out thru bottom of layer L = PSSW4(I,J,L-1)
!
!  GEOS-CHEM
!  (a) Precip coming in  thru top    of layer L = PDOWN(L+1,I,J)
!  (b) Precip going  out thru bottom of layer L = PDOWN(L,  I,J) 
!
!  Thus: Precip coming in:  PSSW4(I,J,L  ) is analogous to PDOWN(L+1,I,J).
!        Precip going  out: PSSW4(I,J,L-1) is analogous to PDOWN(L,I,J  ).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LS    : =T for Large-Scale precipitation; =F otherwise   
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Jacob et al, 2000
!  (2 ) Balkanski et al, 1993
!  (3 ) Giorgi & Chaimedes, 1986
!
!  NOTES: 
!  (1 ) WETDEP should be called twice, once with LS = .TRUE. and once
!        with LS = .FALSE.  This will handle both large-scale and
!        convective precipitation. (bmy, 2/28/00)
!  (2 ) Call subroutine MAKE_QQ to construct the QQ and PDOWN precipitation
!        fields before calling WETDEP. (bmy, 2/28/00)
!  (3 ) Since we are working with an (I,J) column, the ordering of the
!        loops goes J - I - L - N.  Dimension arrays DSTT, PDOWN, QQ
!        to take advantage of this optimal configuration (bmy, 2/28/00)
!  (4 ) Use double-precision exponents to force REAL*8 accuracy
!        (e.g. 1d0, bmy, 2/28/00)
!  (5 ) Diagnostics ND16, ND17, ND18, and ND39 use allocatable arrays 
!        from "diag_mod.f"  (bmy, bey, 3/14/00)
!  (6 ) WETDEP only processes soluble tracers and/or aerosols, as are
!        defined in the NSOL and IDWETD arrays (bmy, 3/14/00)
!  (7 ) Add kludge to prevent wet deposition in the stratosphere (bmy, 6/21/00)
!  (8 ) Removed obsolete code from 10/27/00 (bmy, 12/21/00)
!  (9 ) Remove IREF, JREF -- they are obsolete (bmy, 9/27/01)
!  (10) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (11) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (12) Now reference BXHEIGHT from "dao_mod.f".  Also references routine
!        GEOS_CHEM_STOP from "error_mod.f".  Also fix ND39 diagnostic so that
!        the budget of tracer lost to wetdep is closed.  Now bundled into
!        "wetscav_mod.f".  Now only save to AD16, AD17, AD18, AD39 if L<=LD16,
!        L<=LD17, L<=LD18, and L<=LD39 respectively; this avoids out-of-bounds
!        array errors. Updated comments, cosmetic changes. (qli, bmy, 11/26/02)
!  (13) References IDTSO2, IDTSO4 from "tracerid_mod.f". SO2 in sulfate 
!        chemistry is wet-scavenged on the raindrop and converted to SO4 by 
!        aqueous chem. If evaporation occurs then SO2 comes back as SO4.
!        (rjp, bmy, 3/23/03)  
!  (14) Now use function GET_TS_DYN() from "time_mod.f" (bmy, 3/27/03)
!  (15) Now parallelize over outermost J-loop.  Also move internal routines
!        LS_K_RAIN, LS_F_PRIME, CONV_F_PRIME, and SAFETY to the module, since
!        we cannot call internal routines from w/in a parallel loop. 
!        (bmy, 3/18/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : BXHEIGHT
      USE DIAG_MOD,     ONLY : AD16, AD17, AD18, CT16, CT17, CT18, AD39 
      USE ERROR_MOD,    ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,     ONLY : GET_TS_DYN
      USE TRACERID_MOD, ONLY : IDTSO2, IDTSO4
      
      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT
#     include "CMN_DIAG"  ! Diagnostic arrays and switches 

      ! Arguments
      LOGICAL, INTENT(IN) :: LS

      ! Local Variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL             :: AER

      INTEGER             :: I, IDX, J, L, N, NN
      
      REAL*8              :: Q,     QDOWN,  DT,        DT_OVER_TAU
      REAL*8              :: K,     K_MIN,  K_RAIN,    RAINFRAC
      REAL*8              :: F,     FTOP,   F_PRIME,   WASHFRAC
      REAL*8              :: LOST,  GAINED, MASS_WASH, MASS_NOWASH
      REAL*8              :: ALPHA, ALPHA2, WETLOSS,   L_PLUS_W
      REAL*8              :: XDSTT, TMP

      ! DSTT is the accumulator array of rained-out 
      ! soluble tracer for a given (I,J) column
      REAL*8              :: DSTT(NSOLMAX,LLPAR)

      !=================================================================
      ! WETDEP begins here!
      !
      ! (1)  I n i t i a l i z e   V a r i a b l e s
      !=================================================================

      ! Dynamic timestep [s]
      DT = GET_TS_DYN() * 60d0
      
      ! Select index for diagnostic arrays -- will archive either
      ! large-scale or convective rainout/washout fractions
      IF ( LS ) THEN
         IDX = 1
      ELSE
         IDX = 2
      ENDIF

      !=================================================================
      ! (2)  L o o p   O v e r   (I, J)   S u r f a c e   B o x e s
      !
      ! Process rainout / washout by columns.
      !=================================================================

#if   !defined( LINUX_EFC ) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,           FTOP,      DSTT,     ALPHA    )
!$OMP+PRIVATE( ALPHA2,  F,           F_PRIME,   GAINED,   K_RAIN   )
!$OMP+PRIVATE( LOST,    MASS_NOWASH, MASS_WASH, RAINFRAC, WASHFRAC )
!$OMP+PRIVATE( WETLOSS, L,           Q,         NN,       N        )
!$OMP+PRIVATE( XDSTT,   QDOWN,       AER,       TMP                )
!$OMP+SCHEDULE( DYNAMIC )
#endif
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Zero arrays for this column
         FTOP      = 0d0
         DSTT(:,:) = 0d0

         !==============================================================
         ! (3)  R a i n o u t   F r o m   T o p   L a y e r  (L = LLPAR) 
         !
         ! Assume that rainout is happening in the top layer if 
         ! QQ(LLPAR,I,J) > 0.  In other words, if any precipitation 
         ! forms in grid box (I,J,LLPAR), assume that all of it falls 
         ! down to lower levels.
         !
         ! Soluble gases/aerosols are incorporated into the raindrops 
         ! and are completely removed from grid box (I,J,LLPAR).  There 
         ! is no evaporation and "resuspension" of aerosols during a 
         ! rainout event.
         !
         ! For large-scale (a.k.a. stratiform) precipitation, the first 
         ! order rate constant for rainout in the grid box (I,J,L=LLPAR) 
         ! (cf. Eq. 12, Jacob et al, 2000) is given by:
         !
         !                        Q        
         !    K_RAIN = K_MIN + -------    [units: s^-1]
         !                      L + W    
         !          
         ! and the areal fraction of grid box (I,J,L=LLPAR) that 
         ! is actually experiencing large-scale precipitation 
         ! (cf. Eq. 11, Jacob et al, 2000) is given by: 
         ! 
         !                  Q               
         !    F' =  -------------------   [unitless]
         !           K_RAIN * ( L + W )    
         !
         ! Where:
         !
         !    K_MIN  = minimum value for K_RAIN         
         !           = 1.0e-4 [s^-1]
         !
         !    L + W  = condensed water content in cloud 
         !           = 1.5e-6 [cm3 H2O/cm3 air]
         !
         !    Q = QQ = rate of precipitation formation 
         !             [ cm3 H2O / cm3 air / s ]
         !
         ! For convective precipitation, K_RAIN = 5.0e-3 [s^-1], and the
         ! expression for F' (cf. Eq. 13, Jacob et al, 2000) becomes:
         !
         !                                  { DT        }
         !                    FMAX * Q * MIN{ --- , 1.0 }
         !                                  { TAU       }
         !  F' = ------------------------------------------------------
         !               { DT        }
         !        Q * MIN{ --- , 1.0 }  +  FMAX * K_RAIN * ( L + W )
         !               { TAU       } 
         !
         ! Where:
         !
         !    Q = QQ = rate of precipitation formation 
         !             [cm3 H2O/cm3 air/s]
         !
         !    FMAX   = maximum value for F' 
         !           = 0.3
         !
         !    DT     = dynamic time step from the CTM [s]
         !
         !    TAU    = duration of rainout event 
         !           = 1800 s (30 min)
         !
         !    L + W  = condensed water content in cloud 
         !           = 2.0e-6 [cm3 H2O/cm3 air]
         !
         ! K_RAIN and F' are needed to compute the fraction of tracer
         ! in grid box (I,J,L=LLPAR) lost to rainout.  This is done in 
         ! module routine RAINOUT.
         !==============================================================

         ! Zero variables for this level
         ALPHA       = 0d0
         ALPHA2      = 0d0
         F           = 0d0
         F_PRIME     = 0d0
         GAINED      = 0d0
         K_RAIN      = 0d0
         LOST        = 0d0
         MASS_NOWASH = 0d0
         MASS_WASH   = 0d0
         RAINFRAC    = 0d0
         WASHFRAC    = 0d0
         WETLOSS     = 0d0

         ! Start at the top of the atmosphere
         L = LLPAR

         ! If precip forms at (I,J,L), assume it all rains out
         IF ( QQ(L,I,J) > 0d0 ) THEN

            ! Q is the new precip that is forming within grid box (I,J,L)
            Q = QQ(L,I,J)

            ! Compute K_RAIN and F' for either large-scale or convective
            ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000) 
            IF ( LS ) THEN
               K_RAIN  = LS_K_RAIN( Q )
               F_PRIME = LS_F_PRIME( Q, K_RAIN )
            ELSE
               K_RAIN  = 1.5d-3
               F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
            ENDIF
            
            ! Set F = F', since there is no FTOP at L = LLPAR
            F = F_PRIME

            ! Only compute rainout if F > 0. 
            ! This helps to eliminate unnecessary CPU cycles.
            IF ( F > 0d0 ) THEN 

               ! ND16 diagnostic...save LS and Conv fractions
               IF ( ND16 > 0 .and. L <= LD16 ) THEN
                  AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
                  CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
               ENDIF

               ! ND17 diagnostic...increment counter
               IF ( ND17 > 0 .and. L <= LD17 ) THEN
                  CT17(I,J,L,IDX) = CT17(I,J,L,IDX) + 1
               ENDIF

               ! Loop over soluble tracers and/or aerosol tracers    
               DO NN = 1, NSOL
                  N = IDWETD(NN)

                  ! Call subroutine RAINOUT to compute the fraction
                  ! of tracer lost to rainout in grid box (I,J,L=LLPAR)
                  CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )

                  ! WETLOSS is the amount of soluble tracer 
                  ! lost to rainout in grid box (I,J,L=LLPAR)
                  WETLOSS = STT(I,J,L,N) * RAINFRAC

                  ! Remove rainout losses in grid box (I,J,L=LLPAR) from STT
                  STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

                  ! DSTT is an accumulator array for rained-out tracers.  
                  ! The tracers in DSTT are in the liquid phase and will 
                  ! precipitate to the levels below until a washout occurs.
                  ! Initialize DSTT at (I,J,L=LLPAR) with WETLOSS.
                  DSTT(NN,L) = WETLOSS

                  ! ND17 diagnostic...LS and conv rainout fractions [unitless]
                  IF ( ND17 > 0 .and. L <= LD17 ) THEN
                     AD17(I,J,L,NN,IDX) = 
     &                    AD17(I,J,L,NN,IDX) + RAINFRAC / F 
                  ENDIF

                  ! ND39 diag - save rainout losses in [kg/s]
                  IF ( ND39 > 0 .and. L <= LD39 ) THEN
                     XDSTT          = WETLOSS / DT
                     AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
                  ENDIF

                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     !------------------------------------------------------
                     ! Prior to 3/18/04:
                     ! SAFETY is now a MODULE routine, so we need to pass
                     ! all of the variables for debug output (bmy, 3/18/04)
                     !CALL SAFETY( I, J, L, N, 4 )
                     !------------------------------------------------------
                     CALL SAFETY( I, J, L, N, 4, 
     &                            LS,         PDOWN(L,I,J), QQ(L,I,J),   
     &                            ALPHA,      ALPHA2,       RAINFRAC,    
     &                            WASHFRAC,   MASS_WASH,    MASS_NOWASH, 
     &                            WETLOSS,    GAINED,       LOST,        
     &                            DSTT(NN,:), STT(I,J,:,N) )
                  ENDIF
               ENDDO
            ENDIF

            ! Save FTOP for the next lower level 
            FTOP = F
         ENDIF

         !==============================================================
         ! (4)  R a i n o u t   i n   t h e   M i d d l e   L e v e l s
         ! 
         ! Rainout occurs when there is more precipitation in grid box 
         ! (I,J,L) than in grid box (I,J,L+1).  In other words, rainout 
         ! occurs when the amount of rain falling through the bottom of 
         ! grid box (I,J,L) is more than the amount of rain coming in 
         ! through the top of grid box (I,J,L). 
         !
         ! Thus ( PDOWN(L,I,J) > 0 and QQ(L,I,J) > 0 ) is the 
         ! criterion for Rainout.
         !
         ! Soluble gases/aerosols are incorporated into the raindrops 
         ! and are completely removed from grid box (I,J,L).  There is 
         ! no evaporation and "resuspension" of aerosols during a 
         ! rainout event.
         !
         ! Compute K_RAIN and F' for grid box (I,J,L) exactly as 
         ! described above in Section (4).  K_RAIN and F' depend on 
         ! whether we have large-scale or convective precipitation.
         !
         ! F' is the areal fraction of grid box (I,J,L) that is 
         ! precipitating.  However, the effective area of precipitation
         ! that layer L sees (cf. Eqs. 11-13, Jacob et al, 2000) is 
         ! given by:
         !
         !                   F = MAX( F', FTOP )
         !
         ! where FTOP = F' at grid box (I,J,L+1), that is, for the grid
         ! box immediately above the current grid box.  
         !
         ! Therefore, the effective area of precipitation in grid box
         ! (I,J,L) depends on the area of precipitation in the grid 
         ! boxes above it.
         !
         ! Having computed K_RAIN and F for grid box (I,J,L), call 
         ! routine RAINOUT to compute the fraction of tracer lost to 
         ! rainout conditions.
         !==============================================================
         DO L = LLPAR-1, 2, -1

            ! Zero variables for each level
            ALPHA       = 0d0
            ALPHA2      = 0d0
            F           = 0d0
            F_PRIME     = 0d0
            GAINED      = 0d0
            K_RAIN      = 0d0
            LOST        = 0d0
            MASS_NOWASH = 0d0
            MASS_WASH   = 0d0
            RAINFRAC    = 0d0
            WASHFRAC    = 0d0
            WETLOSS     = 0d0

            ! Rainout criteria
            IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) > 0d0 ) THEN

               ! Q is the new precip that is forming within grid box (I,J,L)
               Q = QQ(L,I,J)

               ! Compute K_RAIN and F' for either large-scale or convective
               ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000) 
               IF ( LS ) THEN
                  K_RAIN  = LS_K_RAIN( Q )
                  F_PRIME = LS_F_PRIME( Q, K_RAIN )
               ELSE
                  K_RAIN  = 1.5d-3
                  F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
               ENDIF

               ! F is the effective area of precip seen by grid box (I,J,L) 
               F = MAX( F_PRIME, FTOP )

               ! Only compute rainout if F > 0. 
               ! This helps to eliminate unnecessary CPU cycles. 
               IF ( F > 0d0 ) THEN

                  ! ND16 diagnostic...save F 
                  IF ( ND16 > 0 .and. L <= LD16 ) THEN
                     AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
                     CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1 
                  ENDIF

                  ! ND17 diagnostic...increment counter
                  IF ( ND17 > 0 .and. L <= LD17 ) THEN
                     CT17(I,J,L,IDX) = CT17(I,J,L,IDX) + 1
                  ENDIF

                  ! Loop over soluble tracers and/or aerosol tracers    
                  DO NN = 1, NSOL
                     N = IDWETD(NN)

                     ! Call subroutine RAINOUT to comptue the fraction
                     ! of tracer lost to rainout in grid box (I,J,L) 
                     CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )

                     ! WETLOSS is the amount of tracer in grid box 
                     ! (I,J,L) that is lost to rainout.
                     WETLOSS = STT(I,J,L,N) * RAINFRAC

                     ! Subtract the rainout loss in grid box (I,J,L) from STT
                     STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

                     ! Add to DSTT the tracer lost to rainout in grid box 
                     ! (I,J,L) plus the tracer lost to rainout from grid box 
                     ! (I,J,L+1), which has by now precipitated down into 
                     ! grid box (I,J,L).  DSTT will continue to accumulate 
                     ! rained out tracer in this manner until a washout 
                     ! event occurs.
                     DSTT(NN,L) = DSTT(NN,L+1) + WETLOSS

                     ! ND17 diagnostic...rainout fractions [unitless]
                     IF ( ND17 > 0 .and. L <= LD17 ) THEN
                        AD17(I,J,L,NN,IDX) = 
     &                       AD17(I,J,L,NN,IDX) + RAINFRAC / F
                     ENDIF

                     ! ND39 diag -- save rainout losses in [kg/s]
                     IF ( ND39 > 0 .and. L <= LD39 ) THEN
                        XDSTT          = WETLOSS / DT
                        AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
                     ENDIF

                     ! Negative tracer...call subroutine SAFETY
                     IF ( STT(I,J,L,N) < 0d0 ) THEN
                        !-----------------------------------------------------
                        ! Prior to 3/18/04:
                        ! SAFETY is now a MODULE routine, so we need to pass
                        ! all of the variables for debug output (bmy, 3/18/04)
                        !CALL SAFETY( I, J, L, N, 5 )
                        !-----------------------------------------------------
                        CALL SAFETY( I, J, L, N, 5, 
     &                               LS,           PDOWN(L,I,J),  
     &                               QQ(L,I,J),    ALPHA,        
     &                               ALPHA2,       RAINFRAC,     
     &                               WASHFRAC,     MASS_WASH,    
     &                               MASS_NOWASH,  WETLOSS,      
     &                               GAINED,       LOST,         
     &                               DSTT(NN,:),   STT(I,J,:,N) )
                     ENDIF
                  ENDDO
               ENDIF

               ! Save FTOP for next level
               FTOP = F 

            !==============================================================
            ! (5)  W a s h o u t   i n   t h e   m i d d l e   l e v e l s
            ! 
            ! Washout occurs when we have evaporation (or no precipitation 
            ! at all) at grid box (I,J,L), but have rain coming down from 
            ! grid box (I,J,L+1).
            !
            ! Thus PDOWN(L,I,J) > 0 and QQ(L,I,J) <= 0 is the criterion 
            ! for Washout.  Also recall that QQ(L,I,J) < 0 denotes 
            ! evaporation and not precipitation.
            !
            ! A fraction ALPHA of the raindrops falling down from grid 
            ! box (I,J,L+1) to grid box (I,J,L) will evaporate along the 
            ! way.  ALPHA is given by:
            !  
            !            precip leaving (I,J,L+1) - precip leaving (I,J,L)
            !  ALPHA = ---------------------------------------------------
            !                     precip leaving (I,J,L+1)
            !
            !
            !                    -QQ(L,I,J) * DZ(I,J,L)
            !        =         --------------------------
            !                        PDOWN(L+1,I,J)
            !
            ! We assume that a fraction ALPHA2 = 0.5 * ALPHA of the 
            ! previously rained-out aerosols and HNO3 coming down from 
            ! level (I,J,L+1) will evaporate and re-enter the atmosphere 
            ! in the gas phase in grid box (I,J,L).  This process is 
            ! called "resuspension".  
            !
            ! For non-aerosol species, the amount of previously rained 
            ! out mass coming down from grid box (I,J,L+1) to grid box 
            ! (I,J,L) is figured into the total mass available for 
            ! washout in grid box (I,J,L).  We therefore do not have to
            ! use the fraction ALPHA2 to compute the resuspension.
            !
            ! NOTE from Hongyu Liu about ALPHA (hyl, 2/29/00)
            ! =============================================================
            ! If our QQ field was perfect, the evaporated amount in grid 
            ! box (I,J,L) would be at most the total rain amount coming 
            ! from above (i.e. PDOWN(I,J,L+1) ). But this is not true for 
            ! the MOISTQ field we are using.  Sometimes the evaporation in 
            ! grid box (I,J,L) can be more than the rain amount from above.  
            ! The reason is our "evaporation" also includes the effect of 
            ! cloud detrainment.  For now we cannot find a way to 
            ! distinguish betweeen the two. We then decided to release 
            ! aerosols in both the detrained air and the evaporated air. 
            !
            ! Therefore, we should use this term in the numerator:
            ! 
            !                -QQ(I,J,L) * BXHEIGHT(I,J,L) 
            !
            ! instead of the term:
            ! 
            !                PDOWN(L+1)-PDOWN(L)
            !
            ! Recall that in make_qq.f we have restricted PDOWN to 
            ! positive values, otherwise, QQ would be equal to 
            ! PDOWN(L+1)-PDOWN(L).           
            !==============================================================
            ELSE IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) <= 0d0 ) THEN

               ! QDOWN is the precip leaving thru the bottom of box (I,J,L)
               ! Q     is the new precip that is forming within box (I,J,L)
               QDOWN = PDOWN(L,I,J)
               Q     = QQ(L,I,J)

               ! Since no precipitation is forming within grid box (I,J,L),
               ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
               F = FTOP
 
               ! Only compute washout if F > 0.
               ! This helps to eliminate needless CPU cycles.
               IF ( F > 0d0 ) THEN

                  ! ND16 diagnostic...save F (fraction of grid box raining)
                  IF ( ND16 > 0d0 .and. L <= LD16 ) THEN
                     AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
                     CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
                  ENDIF

                  ! ND18 diagnostic...increment counter
                  IF ( ND18 > 0 .and. L <= LD18 ) THEN
                     CT18(I,J,L,IDX) = CT18(I,J,L,IDX) + 1
                  ENDIF

                  ! Loop over soluble tracers and/or aerosol tracers    
                  DO NN = 1, NSOL
                     N  = IDWETD(NN)

                     ! Call WASHOUT to compute the fraction of 
                     ! tracer lost to washout in grid box (I,J,L)
                     CALL WASHOUT( I,     J,  L, N, 
     &                             QDOWN, DT, F, WASHFRAC, AER )
                  
                     !=====================================================
                     ! Washout of aerosol tracers -- 
                     ! this is modeled as a kinetic process
                     !=====================================================
                     IF ( AER ) THEN

                        ! ALPHA is the fraction of the raindrops that 
                        ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
                        ALPHA = ( ABS( Q ) * BXHEIGHT(I,J,L) * 100d0 ) / 
     &                            PDOWN(L+1,I,J) 

                        ! ALPHA2 is the fraction of the rained-out aerosols
                        ! that gets resuspended in grid box (I,J,L)
                        ALPHA2 = 0.5d0 * ALPHA

                        ! GAINED is the rained out aerosol coming down from 
                        ! grid box (I,J,L+1) that will evaporate and re-enter 
                        ! the atmosphere in the gas phase in grid box (I,J,L).
                        GAINED = DSTT(NN,L+1) * ALPHA2

                        ! Amount of aerosol lost to washout in grid box
                        ! (qli, bmy, 10/29/02)
                        WETLOSS = STT(I,J,L,N) * WASHFRAC - GAINED

                        ! Remove washout losses in grid box (I,J,L) from STT.
                        ! Add the aerosol that was reevaporated in (I,J,L).
                        ! SO2 in sulfate chemistry is wet-scavenged on the
                        ! raindrop and converted to SO4 by aqeuous chem.
                        ! If evaporation occurs then SO2 comes back as SO4
                        ! (rjp, bmy, 3/23/03)
                        IF ( N == IDTSO2 ) THEN
                            STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4) 
     &                                        + GAINED * 96D0 / 64D0

                            STT(I,J,L,N)      = STT(I,J,L,N) 
     &                                        - ( WETLOSS + GAINED )
                        ELSE
                            STT(I,J,L,N)      = STT(I,J,L,N) - WETLOSS
                        ENDIF

                        ! LOST is the rained out aerosol coming down from
                        ! grid box (I,J,L+1) that will remain in the liquid
                        ! phase in grid box (I,J,L) and will NOT re-evaporate.
                        LOST = DSTT(NN,L+1) - GAINED

                        ! Add the washed out tracer from grid box (I,J,L) to 
                        ! DSTT.  Also add the amount of tracer coming down
                        ! from grid box (I,J,L+1) that does NOT re-evaporate.
                        DSTT(NN,L) = DSTT(NN,L+1) + WETLOSS

                        ! ND18 diagnostic...divide washout fraction by F
                        IF ( ND18 > 0 .and. L <= LD18 ) THEN
                           AD18(I,J,L,NN,IDX) = 
     &                          AD18(I,J,L,NN,IDX) + WASHFRAC / F
                        ENDIF

                     !=====================================================
                     ! Washout of non-aerosol tracers
                     ! This is modeled as an equilibrium process
                     !=====================================================
                     ELSE
                  
                        ! MASS_NOWASH is the amount of non-aerosol tracer in 
                        ! grid box (I,J,L) that is NOT available for washout.
                        MASS_NOWASH = ( 1d0 - F ) * STT(I,J,L,N)
                     
                        ! MASS_WASH is the total amount of non-aerosol tracer
                        ! that is available for washout in grid box (I,J,L).
                        ! It consists of the mass in the precipitating
                        ! part of box (I,J,L), plus the previously rained-out
                        ! tracer coming down from grid box (I,J,L+1).
                        ! (Eq. 15, Jacob et al, 2000).
                        MASS_WASH = ( F * STT(I,J,L,N) ) + DSTT(NN,L+1)

                        ! WETLOSS is the amount of tracer mass in 
                        ! grid box (I,J,L) that is lost to washout.
                        ! (Eq. 16, Jacob et al, 2000)
                        WETLOSS = MASS_WASH * WASHFRAC - DSTT(NN,L+1)

                        ! The tracer left in grid box (I,J,L) is what was
                        ! in originally in the non-precipitating fraction 
                        ! of the box, plus MASS_WASH, less WETLOSS. 
                        STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS  
                  
                        ! Add washout losses in grid box (I,J,L) to DSTT 
                        DSTT(NN,L) = DSTT(NN,L+1) + WETLOSS

                        ! ND18 diagnostic...we don't have to divide the
                        ! washout fraction by F since this is accounted for.
                        IF ( ND18 > 0 .and. L <= LD18 ) THEN
                           AD18(I,J,L,NN,IDX) = 
     &                          AD18(I,J,L,NN,IDX) + WASHFRAC
                        ENDIF
                     ENDIF

                     ! ND39 diag -- save rainout losses in [kg/s]
                     IF ( ND39 > 0 .and. L <= LD39 ) THEN
                        XDSTT          = WETLOSS / DT
                        AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
                     ENDIF
  
                     ! Negative tracer...call subroutine SAFETY
                     IF ( STT(I,J,L,N) < 0d0 ) THEN
                        !-----------------------------------------------------
                        ! Prior to 3/18/04:
                        ! SAFETY is now a MODULE routine, so we need to pass
                        ! all of the variables for debug output (bmy, 3/18/04)
                        !CALL SAFETY( I, J, L, N, 6 )
                        !-----------------------------------------------------
                        CALL SAFETY( I, J, L, N, 6, 
     &                               LS,           PDOWN(L,I,J), 
     &                               QQ(L,I,J),    ALPHA,        
     &                               ALPHA2,       RAINFRAC,     
     &                               WASHFRAC,     MASS_WASH,    
     &                               MASS_NOWASH,  WETLOSS,      
     &                               GAINED,       LOST,         
     &                               DSTT(NN,:),   STT(I,J,:,N) )
                     ENDIF
                  ENDDO  
               ENDIF

               ! Save FTOP for next level
               FTOP = F   

            !===========================================================
            ! (6)  N o   D o w n w a r d   P r e c i p i t a t i o n 
            !
            ! If there is no precipitation leaving grid box (I,J,L), 
            ! then  set F, the effective area of precipitation in grid 
            ! box (I,J,L), to zero.
            !
            ! Also, all of the previously rained-out tracer that is now 
            ! coming down from grid box (I,J,L+1) will evaporate and 
            ! re-enter the atmosphere in the gas phase in grid box 
            ! (I,J,L).  This is called "resuspension".
            !===========================================================
            ELSE IF ( ABS( PDOWN(L,I,J) ) < 1d-30 ) THEN

               ! No precipitation at grid box (I,J,L), thus F = 0
               F = 0d0

               ! Loop over soluble tracers and/or aerosol tracers            
               DO NN = 1, NSOL
                  N = IDWETD(NN) 

                  ! WETLOSS is the amount of tracer in grid box (I,J,L) 
                  ! that is lost to rainout. (qli, bmy, 10/29/02)
                  WETLOSS = -DSTT(NN,L+1)

                  ! All of the rained-out tracer coming from grid box
                  ! (I,J,L+1) goes back into the gas phase at (I,J,L)
                  ! In evap, SO2 comes back as SO4 (rjp, bmy, 3/23/03)
                  IF ( N == IDTSO2 ) THEN
                     STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4) 
     &                                 - ( WETLOSS * 96d0 / 64d0 )
                  ELSE
                     STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
                  ENDIF

                  ! There is nothing rained out/washed out in grid box
                  ! (I,J,L), so set DSTT at grid box (I,J,L) to zero.
		  DSTT(NN,L) = 0d0
                  
                  ! ND39 diag -- save rainout losses in [kg/s]
                  IF ( ND39 > 0 .and. L <= LD39 ) THEN
                     XDSTT          = WETLOSS / DT
                     AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
                  ENDIF

                  ! Negative tracer...call subroutine SAFETY
                  !--------------------------------------------------------
                  ! Prior to 3/18/04:
                  ! SAFETY is now a MODULE routine, so we need to pass
                  ! all of the variables for debug output (bmy, 3/18/04)
                  !IF ( STT(I,J,L,N) < 0d0 ) CALL SAFETY( I, J, L, N, 7 )
                  !--------------------------------------------------------
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 7, 
     &                            LS,         PDOWN(L,I,J), QQ(L,I,J),   
     &                            ALPHA,      ALPHA2,       RAINFRAC,    
     &                            WASHFRAC,   MASS_WASH,    MASS_NOWASH, 
     &                            WETLOSS,    GAINED,       LOST,        
     &                            DSTT(NN,:), STT(I,J,:,N) )
                  ENDIF
               ENDDO

               ! Save FTOP for next level
               FTOP = F
            ENDIF 
         ENDDO                  ! L

         !==============================================================
         ! (8)  W a s h o u t   i n   L e v e l   1
         !
         ! Assume all of the tracer precipitating down from grid box 
         ! (I,J,L=2) to grid box (I,J,L=1) gets washed out in grid box 
         ! (I,J,L=1).
         !==============================================================

         ! Zero variables for this level
         ALPHA       = 0d0
         ALPHA2      = 0d0
         F           = 0d0
         F_PRIME     = 0d0
         GAINED      = 0d0
         K_RAIN      = 0d0
         LOST        = 0d0
         MASS_NOWASH = 0d0
         MASS_WASH   = 0d0
         RAINFRAC    = 0d0
         WASHFRAC    = 0d0
         WETLOSS     = 0d0
         
         ! We are at the surface, set L = 1
         L = 1

         ! Washout at level 1 criteria
         IF ( PDOWN(L+1,I,J) > 0d0 ) THEN

            ! QDOWN is the precip leaving thru the bottom of box (I,J,L+1)
            QDOWN = PDOWN(L+1,I,J)

            ! Since no precipitation is forming within grid box (I,J,L),
            ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
            F = FTOP

            ! Only compute washout if F > 0.
            ! This helps to eliminate unnecessary CPU cycles.
            IF ( F > 0d0 ) THEN

               ! ND16 diagnostic...save F 
               IF ( ND16 > 0 .and. L <= LD16 ) THEN
                  AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
                  CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
               ENDIF

               ! ND18 diagnostic...increment counter
               IF ( ND18 > 0 .and. L <= LD18 ) THEN
                  CT18(I,J,L,IDX) = CT18(I,J,L,IDX) + 1
               ENDIF

               ! Loop over soluble tracers and/or aerosol tracers
               DO NN = 1, NSOL
                  N = IDWETD(NN)

                  ! Call WASHOUT to compute the fraction of tracer 
                  ! in grid box (I,J,L) that is lost to washout.  
                  CALL WASHOUT( I,     J,  L, N, 
     &                          QDOWN, DT, F, WASHFRAC, AER )

                  ! NOTE: for HNO3 and aerosols, there is an F factor
                  ! already present in WASHFRAC.  For other soluble
                  ! gases, we need to multiply by the F (hyl, bmy, 10/27/00)
                  IF ( AER ) THEN
                     WETLOSS = STT(I,J,L,N) * WASHFRAC
                  ELSE
                     WETLOSS = STT(I,J,L,N) * WASHFRAC * F
                  ENDIF

                  ! Subtract WETLOSS from STT
                  STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS     
              
                  ! ND18 diagnostic...LS and conv washout fractions [unitless]
                  IF ( ND18 > 0 .and. L <= LD18 ) THEN

                     ! Only divide WASHFRAC by F for aerosols, since
                     ! for non-aerosols this is already accounted for
                     IF ( AER ) THEN
                        TMP = WASHFRAC / F
                     ELSE
                        TMP = WASHFRAC
                     ENDIF
                     
                     AD18(I,J,L,NN,IDX) = AD18(I,J,L,NN,IDX) + TMP
                  ENDIF

                  ! ND39 diag -- save washout loss in [kg/s]
                  IF ( ND39 > 0 .and. L <= LD39 ) THEN
                     XDSTT          = WETLOSS / DT
                     AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
                  ENDIF

                  !-----------------------------------------------------
                  ! Dirty kludge to prevent wet deposition from removing 
                  ! stuff from stratospheric boxes -- this can cause 
                  ! negative tracer (rvm, bmy, 6/21/00)
                  !
                  IF ( STT(I,J,L,N) < 0d0 .and. L > 23 ) THEN
                      WRITE ( 6, 101 ) I, J, L, N, 7
 101                  FORMAT( 'WETDEP - STT < 0 at ', 3i4, 
     &                        ' for tracer ', i4, 'in area ', i4 )
                      PRINT*, 'STT:', STT(I,J,:,N)
                      STT(I,J,L,N) = 0d0
                  ENDIF
                  !-----------------------------------------------------

                  ! Negative tracer...call subroutine SAFETY
                  !--------------------------------------------------------
                  ! Prior to 3/18/04:
                  ! SAFETY is now a MODULE routine, so we need to pass
                  ! all of the variables for debug output (bmy, 3/18/04)
                  !IF ( STT(I,J,L,N) < 0d0 ) CALL SAFETY( I, J, L, N, 8 )
                  !--------------------------------------------------------
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 8, 
     &                            LS,         PDOWN(L,I,J), QQ(L,I,J),   
     &                            ALPHA,      ALPHA2,       RAINFRAC,    
     &                            WASHFRAC,   MASS_WASH,    MASS_NOWASH, 
     &                            WETLOSS,    GAINED,       LOST,        
     &                            DSTT(NN,:), STT(I,J,:,N) )
                  ENDIF
               ENDDO
            ENDIF    
         ENDIF
      ENDDO	
      ENDDO	
#if   !defined( LINUX_EFC )
!$OMP END PARALLEL DO
#endif

      ! Return to calling program
      END SUBROUTINE WETDEP

!------------------------------------------------------------------------------

      FUNCTION LS_K_RAIN( Q ) RESULT( K_RAIN )
!
!******************************************************************************
!  Function LS_K_RAIN computes K_RAIN, the first order rainout rate constant 
!  for large-scale (a.k.a. stratiform) precipitation (bmy, 3/18/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) Q      (REAL*8) : Rate of precip formation [cm3 H2O/cm3 air/s]
!
!  Function value:
!  ============================================================================
!  (2 ) K_RAIN (REAL*8) : 1st order rainout rate constant [s-1]
!
!  NOTES:
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: Q
      
      ! Function value
      REAL*8             :: K_RAIN

      !==================================================================
      ! LS_K_RAIN begins here!
      !==================================================================

      ! Compute rainout rate constant K in s^-1 (Eq. 12, Jacob et al, 2000).
      ! 1.0d-4 = K_MIN, a minimum value for K_RAIN 
      ! 1.5d-6 = L + W, the condensed water content (liq + ice) in the cloud
      K_RAIN = 1.0d-4 + ( Q / 1.5d-6 ) 
      
      ! Return to WETDEP
      END FUNCTION LS_K_RAIN

!------------------------------------------------------------------------------

      FUNCTION LS_F_PRIME( Q, K_RAIN ) RESULT( F_PRIME )
!
!******************************************************************************
!  Function LS_F_PRIME computes F', the fraction of the grid box that is 
!  precipitating during large scale (a.k.a. stratiform) precipitation.
!  (bmy, 3/18/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) Q       (REAL*8) : Rate of precip formation [cm3 H2O/cm3 air/s]
!  (2 ) K_RAIN  (REAL*8) : 1st order rainout rate constant [s-1]
!
!  Function value:
!  ============================================================================
!  (3 ) F_PRIME (REAL*8) : Fraction of grid box undergoing LS precip [unitless]
!
!  NOTES:
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: Q, K_RAIN

      ! Function value
      REAL*8             :: F_PRIME

      !=================================================================
      ! LS_F_PRIME begins here!
      !=================================================================

      ! Compute F', the area of the grid box undergoing precipitation
      ! 1.5d-6 = L + W, the condensed water content [cm3 H2O/cm3 air]
      F_PRIME = Q / ( K_RAIN * 1.5d-6 )

      ! Return to WETDEP
      END FUNCTION LS_F_PRIME

!------------------------------------------------------------------------------

      FUNCTION CONV_F_PRIME( Q, K_RAIN, DT ) RESULT( F_PRIME )
!
!******************************************************************************
!  Function CONV_F_PRIME computes F', the fraction of the grid box that is 
!  precipitating during convective precipitation. (bmy, 3/18/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) Q       (REAL*8) : Rate of precip formation        [cm3 H2O/cm3 air/s]
!  (2 ) K_RAIN  (REAL*8) : 1st order rainout rate constant [s-1]
!  (3 ) DT      (REAL*8) : Wet deposition timestep         [s]
!
!  Function value:
!  ============================================================================
!  (4 ) F_PRIME (REAL*8) : Frac. of grid box undergoing CONV precip [unitless]
!
!  NOTES:
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: Q, K_RAIN, DT
      
      ! Local variables
      REAL*8             :: TIME

      ! Function value
      REAL*8             :: F_PRIME

      !=================================================================
      ! CONV_F_PRIME begins here!
      !=================================================================
      
      ! Assume the rainout event happens in 30 minutes (1800 s)
      ! Compute the minimum of DT / 1800s and 1.0
      TIME = MIN( DT / 1800d0, 1d0 )

      ! Compute F' for convective precipitation (Eq. 13, Jacob et al, 2000)
      ! 0.3  = FMAX, the maximum value of F' for convective precip
      ! 2d-6 = L + W, the condensed water content [cm3 H2O/cm3 air]
      F_PRIME = ( 0.3d0 * Q * TIME ) / 
     &          ( ( Q * TIME ) + ( 0.3d0 * K_RAIN * 2d-6 ) )

      ! Return to WETDEP
      END FUNCTION CONV_F_PRIME

!------------------------------------------------------------------------------

      SUBROUTINE SAFETY( I,         J,           L,        N,     
     &                   A,         LS,          PDOWN,    QQ,       
     &                   ALPHA,     ALPHA2,      RAINFRAC, WASHFRAC, 
     &                   MASS_WASH, MASS_NOWASH, WETLOSS,  GAINED,   
     &                   LOST,      DSTT,        STT )
!
!******************************************************************************
!  Subroutine SAFETY stops the run with debug output and an error message 
!  if negative tracers are found. (bmy, 3/18/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) Q       (REAL*8) : Rate of precip formation        [cm3 H2O/cm3 air/s]
!  (2 ) K_RAIN  (REAL*8) : 1st order rainout rate constant [s-1]
!  (3 ) DT      (REAL*8) : Wet deposition timestep         [s]
!
!  Function value:
!  ============================================================================
!  (4 ) F_PRIME (REAL*8) : Frac. of grid box undergoing CONV precip [unitless]
!
!  NOTES:
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

#     include "CMN_SIZE"

      ! Arguments
      LOGICAL, INTENT(IN) :: LS
      INTEGER, INTENT(IN) :: I, J, L, N, A
      REAL*8,  INTENT(IN) :: PDOWN,    QQ,       ALPHA,     ALPHA2
      REAL*8,  INTENT(IN) :: RAINFRAC, WASHFRAC, MASS_WASH, MASS_NOWASH
      REAL*8,  INTENT(IN) :: WETLOSS,  GAINED,   LOST,      DSTT(LLPAR)
      REAL*8,  INTENT(IN) :: STT(LLPAR)

      !=================================================================
      ! SAFETY begins here!
      !=================================================================
      
      ! Print line
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Write error message and stop the run
      WRITE ( 6, 100 ) I, J, L, N, A
 100  FORMAT( 'WETDEP - STT < 0 at ', 3i4, ' for tracer ', i4, 
     &        ' in area ', i4 )

      PRINT*, 'LS          : ', LS
      PRINT*, 'PDOWN       : ', PDOWN
      PRINT*, 'QQ          : ', QQ
      PRINT*, 'ALPHA       : ', ALPHA
      PRINT*, 'ALPHA2      : ', ALPHA2
      PRINT*, 'RAINFRAC    : ', RAINFRAC
      PRINT*, 'WASHFRAC    : ', WASHFRAC
      PRINT*, 'MASS_WASH   : ', MASS_WASH
      PRINT*, 'MASS_NOWASH : ', MASS_NOWASH
      PRINT*, 'WETLOSS     : ', WETLOSS
      PRINT*, 'GAINED      : ', GAINED
      PRINT*, 'LOST        : ', LOST
      PRINT*, 'DSTT(NN,:)  : ', DSTT(:)
      PRINT*, 'STT(I,J,:N) : ', STT(:)

      ! Print line
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Deallocate memory and stop
      CALL GEOS_CHEM_STOP

      ! Return to WETDEP
      END SUBROUTINE SAFETY

!------------------------------------------------------------------------------

      SUBROUTINE WETDEPID
!
!******************************************************************************
!  Subroutine WETDEPID sets up the index array of soluble tracers used in
!  the WETDEP routine above (bmy, 11/8/02, 4/5/04)
! 
!  NOTES:
!  (1 ) Now references "tracerid_mod.f".  Also references "CMN" in order to
!        pass variables NSRCX and NTRACE. (bmy, 11/8/02)
!  (2 ) Updated for carbon aerosol & dust tracers (rjp, bmy, 4/5/04)
!******************************************************************************
!
      ! References To F90 modules
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE TRACERID_MOD

#     include "CMN_SIZE" ! Size parameters
#     include "CMN"      ! NSRCX, NTRACE, TCNAME

      ! Local variables
      INTEGER            :: N, NN

      !=================================================================
      ! WETDEPID begins here!
      !=================================================================

      ! Zero NSOL
      NSOL = 0

      ! Sort soluble tracers into IDWETD
      DO N = 1, NTRACE
         IF ( N == IDTPB ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTPB

         ELSE IF ( N == IDTBE7 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTBE7

         ELSE IF ( N == IDTHNO3 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTHNO3

         ELSE IF ( N == IDTH2O2 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTH2O2

         ELSE IF ( N == IDTCH2O ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTCH2O

         ELSE IF ( N == IDTMP ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTMP

         ELSE IF ( N == IDTSO2 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTSO2

         ELSE IF ( N == IDTSO4 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTSO4

         ELSE IF ( N == IDTMSA ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTMSA

         ELSE IF ( N == IDTNH3 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTNH3

         ELSE IF ( N == IDTNH4 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTNH4
            
         ELSE IF ( N == IDTNIT ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTNIT
            
         ELSE IF ( N == IDTBCPI ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTBCPI

         ELSE IF ( N == IDTOCPI ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTOCPI

         ELSE IF ( N == IDTBCPO ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTBCPO

         ELSE IF ( N == IDTOCPO ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTOCPO

         ELSE IF ( N == IDTDST1 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTDST1

         ELSE IF ( N == IDTDST2 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTDST2

         ELSE IF ( N == IDTDST3 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTDST3

         ELSE IF ( N == IDTDST4 ) THEN
            NSOL         = NSOL + 1
            IDWETD(NSOL) = IDTDST4

         ENDIF
      ENDDO

      ! Error check: Make sure that NSOL is less than NSOLMAX
      IF ( NSOL > NSOLMAX ) THEN
         CALL ERROR_STOP( 'NSOL > NSOLMAX!', 'WETDEPID (wetscav_mod.f)')
      ENDIF
      
      ! Also check to see if NSOL is larger than the maximum
      ! number of soluble tracers for a particular simulation
      IF ( NSOL > GET_WETDEP_NMAX( NSRCX ) ) THEN
         CALL ERROR_STOP( 'NSOL > NMAX', 'WETDEPID (wetscav_mod.f)')
      ENDIF

      !=================================================================
      ! Echo list of soluble tracers to the screen
      !=================================================================
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'WETDEPID: List of soluble tracers: '
      WRITE( 6, '(a)'   ) '  #   Name  Tracer Molec.'
      WRITE( 6, '(a)'   ) '            Number Weight'
      WRITE( 6, '(a)'   ) REPEAT( '-', 26 )

      DO NN = 1, NSOL
         N = IDWETD(NN)
         WRITE( 6, '(i3,3x,a4,3x,i3,3x,f6.3)' )
     &        NN, TRIM( TCNAME(N) ), N, TCMASS(N)*1d-3
      ENDDO

      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE WETDEPID

!------------------------------------------------------------------------------

      FUNCTION GET_WETDEP_NMAX( NSRCX ) RESULT ( NMAX )
!
!******************************************************************************
!  Function GET_WETDEP_NMAX returns the maximum number of soluble tracers
!  for a given type of simulation.  Primarily used for allocation of 
!  diagnostic arrays. (bmy, 12/2/02, 4/5/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NSRCX (INTEGER) : GEOS-CHEM simulation type flag
!
!  NOTES:
!  (1 ) Modified to include carbon & dust aerosol tracers (rjp, bmy, 4/5/04)
!******************************************************************************
!
#     include "CMN_SIZE"   ! Size Parameters
#     include "CMN_SETUP"  ! LSULF

      ! Arguments
      INTEGER, INTENT(IN) :: NSRCX 
      
      ! Function value
      INTEGER             :: NMAX

      !=================================================================
      ! GET_WETDEP_NMAX begins here!
      !
      ! NOTE: If you add tracers to a simulation, update as necessary
      !=================================================================
      SELECT CASE ( NSRCX )

         ! Rn-Pb-Be: 2 soluble tracers
         CASE ( 1 )
            NMAX = 2

         ! Fullchem simulation
         CASE ( 3 )
            !--------------------------------------------------------------
            ! Prior to 4/5/04:
            ! Now account for carbon & dust aerosol tracers too
            !IF ( LSULF ) THEN
            !   NMAX = 10       ! Coupled simulation: 10 soluble tracers
            !ELSE 
            !   NMAX = 4        ! Regular simulation:  4 soluble tracers 
            !ENDIF
            !--------------------------------------------------------------
            NMAX = 4                            ! HNO3, H2O2
            IF ( LSULF ) NMAX = NMAX + 6        ! SO2, SO4, MSA, NH3, NH4, NIT
            IF ( LCARB ) NMAX = NMAX + 4        ! BCPI, BCPO, OCPI, OCPO
            IF ( LDUST ) NMAX = NMAX + NDSTBIN  ! plus # of dust bins

         ! Offline sulfate simulation: 7 soluble tracers
         CASE ( 10 ) 
            !--------------------------------------------------------------
            ! Prior to 4/5/04:
            ! Now account for carbon & dust aerosol tracers too
            !NMAX = 7
            !--------------------------------------------------------------
            NMAX = 0
            IF ( LSULF ) NMAX = NMAX + 7        ! add 7 sulfur species
            IF ( LCARB ) NMAX = NMAX + 4        ! add carbonaceous aerosols
            IF ( LDUST ) NMAX = NMAX + NDSTBIN  ! Add number of dust bins

         ! Default: No soluble tracers
         CASE DEFAULT
            NMAX = 0

      END SELECT

      ! Return to calling program
      END FUNCTION GET_WETDEP_NMAX

!------------------------------------------------------------------------------

      FUNCTION GET_WETDEP_NSOL() RESULT( N_SOLUBLE )
!
!******************************************************************************
!  Function GET_WETDEP_NSOL returns NSOL (# of soluble tracers) to a calling
!  program outside WETSCAV_MOD.  This is so that we can keep NSOL declared
!  as a PRIVATE variable. (bmy, 1/10/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: N_SOLUBLE
      
      !=================================================================
      ! GET_WETDEP_NSOL begins here!
      !=================================================================

      ! Get the # of soluble tracers
      N_SOLUBLE = NSOL
     
      ! Return to calling program
      END FUNCTION GET_WETDEP_NSOL

!------------------------------------------------------------------------------

      FUNCTION GET_WETDEP_IDWETD( NWET ) RESULT( N )
!
!******************************************************************************
!  Function GET_WETDEP_IDWETD returns the tracer number of wet deposition 
!  species  NWET.  This is meant to be called outside of WETSCAV_MOD so that 
!  IDWETD can be kept as a PRIVATE variable. (bmy, 1/10/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NWET (INTEGER) : Wet deposition species N
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NWET
      
      ! Function value
      INTEGER             :: N

      !=================================================================
      ! GET_WETDEP_IDWETD begins here!
      !=================================================================

      ! Make sure NWET is valid
      IF ( NWET < 1 .or. NWET > NSOLMAX ) THEN
         CALL ERROR_STOP( 'Invalid value of NWET!', 
     &                    'GET_N_WETDEP (wetscav_mod.f)' )
      ENDIF

      ! Get the tracer # for wet deposition species N
      N = IDWETD(NWET)
     
      ! Return to calling program
      END FUNCTION GET_WETDEP_IDWETD

!------------------------------------------------------------------------------

      SUBROUTINE INIT_WETSCAV
!
!******************************************************************************
!  Subroutine INIT_WETSCAV initializes updraft velocity, cloud liquid water
!  content, cloud ice content, and mixing ratio of water fields, which
!  are used in the wet scavenging routines. (bmy, 2/23/00, 1/27/03)
!
!  NOTES:
!  (1 ) References "e_ice.f" -- routine to compute Eice(T).
!  (2 ) Vud, CLDLIQ, CLDICE, C_H2O are all independent of tracer, so we
!        can compute them once per timestep, before calling the cloud 
!        convection and wet deposition routines.
!  (3 ) Set C_H2O = 0 below -120 Celsius.  E_ICE(T) has a lower limit of
!        -120 Celsius, so temperatures lower than this will cause a stop
!        with an error message. (bmy, 6/15/00)
!  (4 ) Replace {IJL}GLOB with IIPAR,JJPAR,LLPAR.  Also rename PW to P.
!        Remove IREF, JREF, these are obsolete.  Now reference IS_WATER
!        from "dao_mod.f" to determine water boxes. 
!  (5 ) Removed obsolete code from 9/01.  Updated comments and made
!        cosmetic changes. (bmy, 10/24/01)
!  (6 ) Now use routine GET_PCENTER from "pressure_mod.f" to compute the
!        pressure at the midpoint of grid box (I,J,L).  Also removed P and
!        SIG from the argument list (dsa, bdf, bmy, 8/20/02)
!  (7 ) Now reference T from "dao_mod.f".  Updated comments.  Now allocate
!        Vud, C_H2O, CLDLIQ and CLDICE here on the first call.  Now references
!        ALLOC_ERR from "error_mod.f".  Now set H2O2s and SO2s to the initial
!        values from for the first call to COMPUTE_F .  Now call WETDEPID
!        on the first call to initialize the wetdep index array. (bmy, 1/27/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : T, IS_WATER
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TRACERID_MOD, ONLY : IDTH2O2, IDTSO2

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT

      ! Local variables
      INTEGER             :: I, J, L, AS
      REAL*8              :: PL, TC, TK
      LOGICAL, SAVE       :: FIRST = .TRUE.      
            
      !=================================================================
      ! INIT_WETSCAV begins here!
      !=================================================================
      IF ( FIRST ) THEN

         ! Allocate Vud on first call
         ALLOCATE( Vud( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'Vud' )
         Vud = 0d0

         ! Allocate C_H2O on first call
         ALLOCATE( C_H2O( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'C_H2O' )
         C_H2O = 0d0

         ! Allocate CLDLIQ on first call
         ALLOCATE( CLDLIQ( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDLIQ' )
         CLDLIQ = 0d0

         ! Allocate CLDICE on first call
         ALLOCATE( CLDICE( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDICE' )
         CLDICE = 0d0

         ! Allocate H2O2s for wet deposition 
         ALLOCATE( H2O2s( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'H2O2s' )

         ! Set H2O2s to the initial H2O2 from STT, so that we will have
         ! nonzero values for the first call to COMPUTE_F (bmy, 1/14/03)
         IF ( IDTH2O2 > 0 ) THEN
            H2O2s = STT(:,:,:,IDTH2O2)
         ELSE
            H2O2s = 0d0
         ENDIF

         ! Allocate SO2s for wet deposition
         ALLOCATE( SO2s( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2s' )

         ! Set SO2s to the initial SO2 from STT, so that we will have
         ! nonzero values for the first call to COMPUTE_F (bmy, 1/14/03)
         IF ( IDTSO2 > 0 ) THEN
            SO2s = STT(:,:,:,IDTSO2)
         ELSE
            SO2s = 0d0
         ENDIF

         ! Index soluble tracers for cloud scavenging and wetdep
         CALL WETDEPID

         ! Reset flag
         FIRST = .FALSE. 
      ENDIF

      !=================================================================
      ! Compute Vud, CLDLIQ, CLDICE, C_H2O, following Jacob et al, 2000.
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TK, TC, PL )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Compute Temp [K], Temp [C], and Pressure [hPa]
         TK = T(I,J,L)
         TC = TK - 273.15d0
         PL = GET_PCENTER(I,J,L)

         !==============================================================
         ! Compute Vud -- 5 m/s over oceans, 10 m/s over land
         ! Assume Vud is the same at all altitudes; the array can be 2-D
         !==============================================================
         IF ( L == 1 ) THEN
            IF ( IS_WATER( I, J ) ) THEN
               Vud(I,J) = 5d0
            ELSE
               Vud(I,J) = 10d0
            ENDIF
         ENDIF

         !==============================================================
         ! CLDLIQ, the cloud liquid water content [cm3 H2O/cm3 air], 
         ! is a function of the local Kelvin temperature:
         !  
         !    CLDLIQ = 2e-6                    [     T >= 268 K    ]
         !    CLDLIQ = 2e-6 * ((T - 248) / 20) [ 248 K < T < 268 K ]
         !    CLDLIQ = 0                       [     T <= 248 K    ]
         !==============================================================
         IF ( TK >= 268d0 ) THEN
            CLDLIQ(I,J,L) = 2d-6

         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
            CLDLIQ(I,J,L) = 2d-6 * ( ( TK - 248d0 ) / 20d0 )

         ELSE
            CLDLIQ(I,J,L) = 0d0
            
         ENDIF
           
         !=============================================================
         ! CLDICE, the cloud ice content [cm3 ice/cm3 air] is given by:
         !
         !    CLDICE = 2e-6 - CLDLIQ
         !=============================================================
         CLDICE(I,J,L) = 2d-6 - CLDLIQ(I,J,L)

         !=============================================================
         ! C_H2O is the mixing ratio of water computed from Eice(Tc) 
         ! (the sat. vapor pressure of ice at the local Celsius temp.)
         !
         ! C_H2O is given by Dalton's Law as:
         !
         !       C_H2O = Eice( Tc(I,J,L) ) / P(I,J,L)
         !
         ! where P(L) = pressure in grid box (I,J,L)
         !
         ! and   Tc(I,J,L) is the Celsius temp. of grid box (I,J,L).
         !
         ! and   Eice( Tc(I,J,L) ) is the saturation vapor pressure 
         !       of ice at temperature Tc(I,J,L) -- computed in 
         !       routine E_ICE above.
         !
         ! NOTE: Do not call E_ICE for T > 0C, since C_H2O only needs 
         !       to be defined for T <= 0 anyway.
         !==============================================================
         IF ( TC > 0d0 .or. TC < -120d0 ) THEN
            C_H2O(I,J,L) = 0d0
         ELSE
            C_H2O(I,J,L) = E_ICE( TC ) / PL
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE INIT_WETSCAV

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_WETSCAV

      !=================================================================
      ! Subroutine CLEANUP_WETSCAV deallocates arrays for 
      ! wet scavenging / wet deposition
      !=================================================================
      IF ( ALLOCATED( Vud    ) ) DEALLOCATE( Vud    )
      IF ( ALLOCATED( C_H2O  ) ) DEALLOCATE( C_H2O  )
      IF ( ALLOCATED( CLDLIQ ) ) DEALLOCATE( CLDLIQ )
      IF ( ALLOCATED( CLDICE ) ) DEALLOCATE( CLDICE )
      IF ( ALLOCATED( PDOWN  ) ) DEALLOCATE( PDOWN  )
      IF ( ALLOCATED( QQ     ) ) DEALLOCATE( QQ     )
      IF ( ALLOCATED( H2O2s  ) ) DEALLOCATE( H2O2s  )
      IF ( ALLOCATED( SO2s   ) ) DEALLOCATE( SO2s   )

      ! Return to calling program
      END SUBROUTINE CLEANUP_WETSCAV

!-----------------------------------------------------------------------------

      END MODULE WETSCAV_MOD
