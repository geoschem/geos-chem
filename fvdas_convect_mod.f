! $Id: fvdas_convect_mod.f,v 1.10 2006/03/24 20:22:47 bmy Exp $
      MODULE FVDAS_CONVECT_MOD
!
!******************************************************************************
!  Module FVDAS_CONVECT_MOD contains routines (originally from NCAR) which 
!  perform shallow and deep convection for the GEOS-4/fvDAS met fields.  
!  These routines account for shallow and deep convection, plus updrafts 
!  and downdrafts.  (pjr, dsa, bmy, 6/26/03, 2/1/06)
!  
!  Module Variables:
!  ============================================================================
!  (1 ) RLXCLM   (LOGICAL) : Logical to relax column versus cloud triplet
!  (2 ) LIMCNV   (INTEGER) : Maximum CTM level for HACK convection
!  (3 ) CMFTAU   (REAL*8 ) : Characteristic adjustment time scale for HACK [s]
!  (4 ) EPS      (REAL*8 ) : A very small number [unitless]
!  (5 ) GRAV     (REAL*8 ) : Gravitational constant [m/s2]
!  (6 ) SMALLEST (REAL*8 ) : The smallest double-precision number 
!  (7 ) TINYNUM  (REAL*8 ) : 2 times the machine epsilon for dble-precision
!  (8 ) TINYALT  (REAL*8 ) : arbitrary small num used in transport estimates
!
!  Module Routines:
!  ============================================================================
!  (1 ) INIT_FVDAS_CONVECT : Initializes fvDAS convection scheme
!  (2 ) FVDAS_CONVECT      : fvDAS convection routine, called from MAIN 
!  (3 ) HACK_CONV          : HACK convection scheme routine
!  (4 ) ARCCONVTRAN        : Sets up fields for ZHANG/MCFARLANE convection
!  (5 ) CONVTRAN           : ZHANG/MCFARLANE convection scheme routine
!  (6 ) WHENFGT            : Test function -- not sure what this does?
!
!  GEOS-CHEM modules referenced by fvdas_convect_mod.f:
!  ============================================================================
!  (1 ) pressure_mod.f     : Module containing routines to compute P(I,J,L)
!
!  NOTES: 
!  (1 ) Contains new updates for GEOS-4/fvDAS convection.  Also added OpenMP
!        parallel loop over latitudes in FVDAS_CONVECT. (swu, bmy, 1/21/04)
!  (2 ) Now prevent FTMP, QTMP arrays from being held PRIVATE w/in the
!        parallel loop in routine DO_CONVECTION (bmy, 7/20/04)
!  (3 ) Now pass wet-scavenged Hg2 to "ocean_mercury_mod.f" (sas, bmy, 1/21/05)
!  (4 ) Rewrote parallel loops to avoid problems w/ OpenMP.  Also modified
!        for updated Hg simulation. (cdh, bmy, 2/1/06)
!******************************************************************************
!
      IMPLICIT NONE
      
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "fvdas_convect_mod.f"
      !=================================================================

      ! Declare everything PRIVATE ...
      PRIVATE
      
      ! ... except these routines
      PUBLIC :: INIT_FVDAS_CONVECT
      PUBLIC :: FVDAS_CONVECT

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Variables
      INTEGER            :: LIMCNV      
      
      ! Constants
      LOGICAL, PARAMETER :: RLXCLM   = .TRUE.
      REAL*8,  PARAMETER :: CMFTAU   = 3600.d0
      REAL*8,  PARAMETER :: EPS      = 1.0d-13   
      REAL*8,  PARAMETER :: GRAV     = 9.8d0
      REAL*8,  PARAMETER :: SMALLEST = TINY(1D0)
      REAL*8,  PARAMETER :: TINYALT  = 1.0d-36       
      REAL*8,  PARAMETER :: TINYNUM  = 2*EPSILON(1D0)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE INIT_FVDAS_CONVECT
!
!******************************************************************************
!  Subroutine INIT_FVDAS_CONVECT initializes the HACK and 
!  ZHANG/MCFARLANE convection schemes for GEOS-4/fvDAS met fields. 
!  (dsa, swu, bmy, 6/26/03, 12/17/03)
!
!  NOTES:
!  (1 ) Now compute HYPI in a more efficient way (bmy, 12/17/03)
!******************************************************************************
!
      ! References to F90 modules
      USE PRESSURE_MOD, ONLY : GET_PEDGE

#     include "CMN_SIZE"  ! Size parameters
      
      ! Local variables
      INTEGER            :: I, J, L, L2
      REAL*8             :: HYPI(LLPAR+1)

      !=================================================================
      ! INIT_FVDAS_CONVECT begins here!
      !
      ! Find the model level that roughly corresponds to 40 hPa and
      ! only let convection take place below that level (LIMCNV)
      !=================================================================
      
      ! Take I, J at midpt of region 
      ! (For global grids, this should be the equatorial box)
      I = IIPAR / 2
      J = JJPAR / 2

      ! Construct array of pressure edges [hPa] for column (I,J) 
      DO L = 1, LLPAR+1
         L2       = (LLPAR+1) - L + 1
         HYPI(L2) = GET_PEDGE(I,J,L)
      ENDDO

      ! Limit convection to regions below 40 hPa
      IF ( HYPI(1) >= 40d0 ) THEN
         LIMCNV = 1
      ELSE
         DO L = 1, LLPAR
            IF ( HYPI(L) < 40d0 .AND. HYPI(L+1) >= 40d0 ) THEN
               LIMCNV = L
               GOTO 10
            ENDIF
         ENDDO
         LIMCNV = LLPAR + 1
      ENDIF

      ! Exit loop
 10   CONTINUE

      !=================================================================
      ! Echo output
      !=================================================================

      WRITE( 6, 100 )  LIMCNV, HYPI(LIMCNV) 
 100  FORMAT( '     - GEOS-4 convection is capped at L = ', i3, 
     &       ', or approx ', f6.1, ' hPa' )

      ! Return to calling program
      END SUBROUTINE INIT_FVDAS_CONVECT

!------------------------------------------------------------------------------

      SUBROUTINE FVDAS_CONVECT( TDT,   NTRACE, Q,     RPDEL, ETA, 
     &                          BETA,  MU,     MD,    EU,    DP,    
     &                          NSTEP, FRACIS, TCVV,  INDEXSOL )
!
!******************************************************************************
!  Subroutine FVDAS_CONVECT is the convection driver routine for GEOS-4/fvDAS
!  met fields.  It calls both HACK and ZHANG/MCFARLANE convection schemes.
!  (pjr, dsa, bmy, 6/26/03, 12/13/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TDT     (REAL*8 ) : 2 * delta-T                          [s]
!  (2 ) NTRACE  (INTEGER) : Number of tracers to transport       [unitless]
!  (3 ) Q       (REAL*8 ) : Array of transported tracers         [v/v]
!  (4 ) RPDEL   (REAL*8 ) : 1 / DP                               [1/hPa]
!  (5 ) ETA     (REAL*8 ) : GMAO Hack convective mass flux       [kg/m2/s]
!  (6 ) BETA    (REAL*8 ) : GMAO Hack overshoot parameter        [unitless]
!  (7 ) MU      (REAL*8 ) : GMAO updraft mass flux   (ZMMU)      [Pa/s]
!  (8 ) MD      (REAL*8 ) : GMAO downdraft mass flux (ZMMD)      [Pa/s]
!  (9 ) EU      (REAL*8 ) : GMAO updraft entrainment (ZMEU)      [Pa/s]
!  (10) DP      (REAL*8 ) : Delta-pressure between level edges   [Pa]
!  (11) NSTEP   (INTEGER) : Time step index                      [unitless]
!  (12) FRACIS  (REAL*8 ) : Fraction of tracer that is insoluble [unitless]
!  (13) TCVV    (REAL*8 ) : Array of Molwt(AIR)/molwt(Tracer)    [unitless]
!  (14) INDEXSOL(INTEGER) : Index array of soluble tracers       [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) Q      (REAL*8 ) : Modified tracer array                 [v/v]
! 
!  Important Local Variables:
!  ============================================================================
!  (1 ) LENGATH(INTEGER)  : Number of lons where deep conv. happens at lat=J
!  (2 ) IDEEP  (INTEGER)  : Lon indices where deep convection happens at lat=J
!  (3 ) JT     (INTEGER)  : Cloud top layer for columns undergoing conv.
!  (4 ) MX     (INTEGER)  : Cloud bottom layer for columns undergoing conv.
!  (5 ) DSUBCLD(REAL*8 )  : Delta pressure from cloud base to sfc
!  (6 ) DU     (REAL*8 )  : Mass detraining from updraft (lon-alt array)
!  (7 ) ED     (REAL*8 )  : Mass entraining from downdraft (lon-alt array)
!  (8 ) DPG    (REAL*8 )  : gathered .01*dp (lon-alt array)
!  (8 ) EUG    (REAL*8 )  : gathered eu (lon-alt array) 
!  (9 ) MUG    (REAL*8 )  : gathered mu (lon-alt array)   
!  (10) MDG    (REAL*8 )  : gathered md (lon-alt array)
!
!  NOTES:
!  (1 ) Added TCVV and INDEXSOL to the arg list and in the call to CONVTRAN.  
!        Now perform convection in a loop over NSTEP iterations.  Added
!        an OpenMP parallel loop over latitude.   Removed IL1G and IL2G,
!        since these are no longer needed in this routine.  Now put NTRACE 
!        before Q on the arg list. (bmy, 1/21/04)
!  (2 ) Handle parallel loops differently for Intel Fortran Compilers, since
!        for some reason the code dies if large arrays (QTMP, FTMP) are held
!        PRIVATE in parallel loops. (bmy, 7/20/04)
!  (3 ) Added LINUX_IFORT switch for Intel v8/v9 compilers (bmy, 10/18/05)
!  (4 ) Rewrote parallel loops so that we pass entire arrays to the various
!        subroutines instead of array slices such as (:,J,:).  This can cause
!        problems with OpenMP for some compilers. (bmy, 12/13/05)
!******************************************************************************

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NSTEP, NTRACE             
      INTEGER, INTENT(IN)    :: INDEXSOL(NTRACE) 
      REAL*8,  INTENT(IN)    :: TDT                
      REAL*8,  INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NTRACE)
      REAL*8,  INTENT(IN)    :: RPDEL(IIPAR,JJPAR,LLPAR)  
      REAL*8,  INTENT(IN)    :: ETA(IIPAR,JJPAR,LLPAR)    
      REAL*8,  INTENT(IN)    :: BETA(IIPAR,JJPAR,LLPAR)   
      REAL*8,  INTENT(IN)    :: MU(IIPAR,JJPAR,LLPAR)     
      REAL*8,  INTENT(IN)    :: MD(IIPAR,JJPAR,LLPAR)     
      REAL*8,  INTENT(IN)    :: EU(IIPAR,JJPAR,LLPAR)     
      REAL*8,  INTENT(IN)    :: DP(IIPAR,JJPAR,LLPAR)     
      REAL*8,  INTENT(IN)    :: FRACIS(IIPAR,JJPAR,LLPAR,NTRACE) 
      REAL*8,  INTENT(IN)    :: TCVV(NTRACE)

      ! Local variables
      INTEGER                :: I, J, L, N, LENGATH, ISTEP
      INTEGER                :: JT(IIPAR)
      INTEGER                :: MX(IIPAR)
      INTEGER                :: IDEEP(IIPAR)
      REAL*8                 :: DSUBCLD(IIPAR)
      REAL*8                 :: DPG(IIPAR,LLPAR)
      REAL*8                 :: DUG(IIPAR,LLPAR)
      REAL*8                 :: EDG(IIPAR,LLPAR)
      REAL*8                 :: EUG(IIPAR,LLPAR)
      REAL*8                 :: MDG(IIPAR,LLPAR)
      REAL*8                 :: MUG(IIPAR,LLPAR)

      !=================================================================
      ! FVDAS_CONVECT begins here!
      !=================================================================

      ! Loop over latitudes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J,       MUG, MDG, DUG,   EUG,     EDG,  DPG )
!$OMP+PRIVATE( DSUBCLD, JT,  MX,  IDEEP, LENGATH, ISTEP     )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Gather mass flux arrays, compute mass fluxes, and determine top
         ! and bottom of Z&M convection.  LENGATH = # of longitudes in the
         ! band I=1,IIPAR where deep convection happens at latitude J.
         CALL ARCONVTRAN( J,   DP,  MU,    MD,      
     &                    EU,  MUG, MDG,   DUG, 
     &                    EUG, EDG, DPG,   DSUBCLD, 
     &                    JT,  MX,  IDEEP, LENGATH )

         ! Loop over internal convection timestep
         DO ISTEP = 1, NSTEP 
               
            !-----------------------------------
            ! ZHANG/MCFARLANE (deep) convection 
            !-----------------------------------

            ! Only call CONVTRAN where convection happens
            ! (i.e. at latitudes where LENGATH > 0)
            IF ( LENGATH > 0 ) THEN
               CALL CONVTRAN( J,     NTRACE,    Q,      MUG,  MDG,      
     &                        DUG,   EUG,       EDG,    DPG,  DSUBCLD,  
     &                        JT,    MX,        IDEEP,  1,    LENGATH, 
     &                        NSTEP, 0.5D0*TDT, FRACIS, TCVV, INDEXSOL )
            ENDIF

            !-----------------------------------
            ! HACK (shallow) convection                   
            !-----------------------------------
            CALL HACK_CONV( J, TDT, RPDEL, ETA, BETA, NTRACE, Q ) 
 
         ENDDO 

      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE FVDAS_CONVECT

!------------------------------------------------------------------------------

      SUBROUTINE HACK_CONV( J, TDT, RPDEL, ETA, BETA, NTRACE, Q )
!
!******************************************************************************
!  Subroutine HACK_CONV computes the convective mass flux adjustment to all 
!  tracers using the convective mass fluxes and overshoot parameters for the 
!  Hack scheme. (pjr, dsa, bmy, 6/26/03, 12/13/05)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) J      (INTEGER) : GEOS-CHEM Latitude index               [unitless]
!  (2 ) TDT    (REAL*8)  : 2 delta-t                              [s       ]
!  (3 ) RPDEL  (REAL*8)  : Reciprocal of pressure-thickness array [1/hPa   ]
!  (4 ) ETA    (REAL*8)  : GMAO Hack convective mass flux (HKETA) [kg/m2/s ]
!  (5 ) BETA   (REAL*8)  : GMAO Hack overshoot parameter (HKBETA) [unitless]
!  (6 ) NTRACE (INTEGER) : Number of tracers in the Q array       [unitless]
!  (7 ) Q      (REAL*8)  : Tracer concentrations                  [v/v     ]   
!  
!  Arguments as Output:
!  ============================================================================
!  (7 ) Q      (REAL*8)  : Modified tracer concentrations         [v/v     ]
!
!  Important Local Variables:
!  ============================================================================
!  (1 ) INDX1  (INTEGER) : Longitude indices for condition true
!  (2 ) ADJFAC (REAL*8 ) : Adjustment factor (relaxation related)
!  (3 ) BOTFLX (REAL*8 ) : Bottom constituent mixing ratio flux
!  (4 ) CMRC   (REAL*8 ) : constituent mix rat ("in-cloud")
!  (5 ) CMRH   (REAL*8 ) : interface constituent mixing ratio 
!  (6 ) DCMR1  (REAL*8 ) : Q convective change (lower lvl)
!  (7 ) DCMR2  (REAL*8 ) : Q convective change (mid level)
!  (8 ) DCMR3  (REAL*8 ) : Q convective change (upper lvl)
!  (9 ) EFAC1  (REAL*8 ) : Ratio q to convectively induced chg (btm lvl)
!  (10) EFAC2  (REAL*8 ) : Ratio q to convectively induced chg (mid lvl)
!  (11) EFAC3  (REAL*8 ) : Ratio q to convectively induced chg (top lvl)
!  (12) ETAGDT (REAL*8 ) : ETA * GRAV * DT
!  (13) TOPFLX (REAL*8 ) : Top constituent mixing ratio flux
!
!  NOTES:
!  (1 ) Updated comments.  Added NTRACE as an argument.  Now also force 
!        double-precision with the "D" exponents.  (bmy, 6/26/03)
!  (2 ) Now pass J via the arg list.  Now dimension RPDEL, ETA, BETA, and Q
!        with and make all input arrays dimensioned
!        with (IIPAR,JJPAR,LLPAR,...) to avoid seg fault error in OpenMP
!        on certain platforms.
!******************************************************************************
!
#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: J, NTRACE
      REAL*8,  INTENT(IN)    :: TDT
      REAL*8,  INTENT(IN)    :: RPDEL(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: ETA(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: BETA(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NTRACE)

      ! Local variables
      INTEGER                :: I, II, K, LEN1, M
      INTEGER                :: INDX1(IIPAR)
      REAL*8                 :: ADJFAC, BOTFLX, TOPFLX              
      REAL*8                 :: EFAC1,  EFAC2,  EFAC3
      REAL*8                 :: CMRC(IIPAR)         
      REAL*8                 :: CMRH(IIPAR,LLPAR+1)   
      REAL*8                 :: DCMR1(IIPAR)        
      REAL*8                 :: DCMR2(IIPAR)        
      REAL*8                 :: DCMR3(IIPAR)        
      REAL*8                 :: ETAGDT(IIPAR)      

      !=================================================================
      ! HACK_CONV begins here!
      !
      ! Ensure that characteristic adjustment time scale (cmftau) 
      ! assumed in estimate of eta isn't smaller than model time scale 
      ! (tdt).  The time over which the convection is assumed to act 
      ! (the adjustment time scale) can be applied with each application 
      ! of the three-level cloud model, or applied to the column 
      ! tendencies after a "hard" adjustment (i.e., on a 2-delta t 
      ! time scale) is evaluated
      !=================================================================
      IF ( RLXCLM ) THEN
         ADJFAC = TDT / ( MAX( TDT, CMFTAU ) )
      ELSE
         ADJFAC = 1.0D0
      ENDIF

      !=================================================================
      ! Begin moist convective mass flux adjustment procedure. 
      ! The formalism ensures that negative cloud liquid water can 
      ! never occur.
      !=================================================================
      DO 70 K = LLPAR-1, LIMCNV+1, -1
         LEN1 = 0
         DO I = 1, IIPAR
            IF ( ETA(I,J,K) /= 0.0 ) THEN
               ETAGDT(I)   = ETA(I,J,K) * GRAV * TDT *0.01d0  ![hPa]
               LEN1        = LEN1 + 1
               INDX1(LEN1) = I
            ELSE
               ETAGDT(I)   = 0.0d0
            ENDIF
         ENDDO
         
         ! Skip to next level
         IF ( LEN1 <= 0 ) GOTO 70

         !==============================================================
         ! Next, convectively modify passive constituents.  For now, 
         ! when applying relaxation time scale to thermal fields after 
         ! entire column has undergone convective overturning, 
         ! constituents will be mixed using a "relaxed" value of the mass
         ! flux determined above.  Although this will be inconsistent 
         ! with the treatment of the thermal fields, it's computationally 
         ! much cheaper, no more-or-less justifiable, and consistent with 
         ! how the history tape mass fluxes would be used in an off-line 
         ! mode (i.e., using an off-line transport model)
         !==============================================================
         DO 50 M = 1, NTRACE
            DO 40 II = 1, LEN1
               I = INDX1(II)

               ! If any of the reported values of the constituent is 
               ! negative in the three adjacent levels, nothing will 
               ! be done to the profile.  Skip to next longitude.
               IF ( ( Q(I,J,K+1,M) < 0.0 )  .OR. 
     &              ( Q(I,J,K,M)   < 0.0 )  .OR.
     &              ( Q(I,J,K-1,M) < 0.0 ) ) GOTO 40

               ! Specify constituent interface values (linear interpolation)
               CMRH(I,K  ) = 0.5d0 *( Q(I,J,K-1,M) + Q(I,J,K  ,M) )
               CMRH(I,K+1) = 0.5d0 *( Q(I,J,K  ,M) + Q(I,J,K+1,M) )
              
               CMRC(I) = Q(I,J,K+1,M)

               ! Determine fluxes, flux divergence => changes due to 
               ! convection.  Logic must be included to avoid producing 
               ! negative values. A bit messy since there are no a priori 
               ! assumptions about profiles.  Tendency is modified (reduced) 
               ! when pending disaster detected.
               BOTFLX   = ETAGDT(I)*(CMRC(I) - CMRH(I,K+1))*ADJFAC
               TOPFLX   = BETA(I,J,K)*ETAGDT(I)*
     &                    (CMRC(I)-CMRH(I,K))*ADJFAC
               DCMR1(I) = -BOTFLX*RPDEL(I,J,K+1)
               EFAC1    = 1.0d0
               EFAC2    = 1.0d0
               EFAC3    = 1.0d0
               
               IF ( Q(I,J,K+1,M)+DCMR1(I) < 0.0 ) THEN
                  EFAC1 = MAX(TINYALT,ABS(Q(I,J,K+1,M)/DCMR1(I)) - EPS)
               ENDIF

               IF ( EFAC1 == TINYALT .OR. EFAC1 > 1.0 ) EFAC1 = 0.0D0
               DCMR1(I) = -EFAC1*BOTFLX*RPDEL(I,J,K+1)
               DCMR2(I) = (EFAC1*BOTFLX - TOPFLX)*RPDEL(I,J,K)
               
               IF ( Q(I,J,K,M)+DCMR2(I) < 0.0 ) THEN
                  EFAC2 = MAX(TINYALT,ABS(Q(I,J,K,M)/DCMR2(I)) - EPS)
               ENDIF
               
               IF ( EFAC2 == TINYALT .OR. EFAC2 > 1.0 ) EFAC2 = 0.0D0
               DCMR2(I) = (EFAC1*BOTFLX - EFAC2*TOPFLX)*RPDEL(I,J,K)
               DCMR3(I) = EFAC2*TOPFLX*RPDEL(I,J,K-1)

               IF ( Q(I,J,K-1,M)+DCMR3(I) < 0.0 ) THEN
                  EFAC3 = MAX(TINYALT,ABS(Q(I,J,K-1,M)/DCMR3(I)) - EPS)
               ENDIF

               IF ( EFAC3 == TINYALT .OR. EFAC3 > 1.0 ) EFAC3 = 0.0D0
               EFAC3    = MIN(EFAC2,EFAC3)
               DCMR2(I) = (EFAC1*BOTFLX - EFAC3*TOPFLX)*RPDEL(I,J,K)
               DCMR3(I) = EFAC3*TOPFLX*RPDEL(I,J,K-1)
               
               Q(I,J,K+1,M) = Q(I,J,K+1,M) + DCMR1(I)
               Q(I,J,K  ,M) = Q(I,J,K  ,M) + DCMR2(I)
               Q(I,J,K-1,M) = Q(I,J,K-1,M) + DCMR3(I)
 40         CONTINUE
 50      CONTINUE
 70   CONTINUE
      
      ! Return to calling program
      END SUBROUTINE HACK_CONV

!------------------------------------------------------------------------------

      SUBROUTINE ARCONVTRAN( J,   DP,  MU,    MD,  
     &                       EU,  MUG, MDG,   DUG, 
     &                       EUG, EDG, DPG,   DSUBCLD, 
     &                       JTG, JBG, IDEEP, LENGATH )
!
!******************************************************************************
!  Subroutine ARCONVTRAN sets up the convective transport using archived mass
!  fluxes from the ZHANG/MCFARLANE convection scheme.  The setup involves:
!    (1) Gather mass flux arrays.
!    (2) Calc the mass fluxes that are determined by mass balance.
!    (3) Determine top and bottom of convection.
!  (pjr, dsa, swu, bmy, 6/26/03, 12/13/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J       (INTEGER) : GEOS-CHEM latitude index          [unitless]
!  (2 ) DP      (REAL*8 ) : Delta pressure between interfaces [Pa      ]      
!  (3 ) MU      (REAL*8 ) : Mass flux up                      [kg/m2/s ] 
!  (4 ) MD      (REAL*8 ) : Mass flux down                    [kg/m2/s ] 
!  (5 ) EU      (REAL*8 ) : Mass entraining from updraft      [1/s     ]     
!
!  Arguments as Output:
!  ============================================================================
!  (6 ) MUG     (REAL*8 ) : Gathered mu (lon-alt array)
!  (7 ) MDG     (REAL*8 ) : Gathered md (lon-alt array)
!  (8 ) DUG     (REAL*8 ) : Mass detraining from updraft (lon-alt array)
!  (9 ) EUG     (REAL*8 ) : Gathered eu (lon-alt array)
!  (10) EDG     (REAL*8 ) : Mass entraining from downdraft (lon-alt array)
!  (11) DPG     (REAL*8 ) : Gathered .01*dp (lon-alt array)
!  (12) DSUBCLD (REAL*8 ) : Delta pressure from cloud base to sfc (lon-alt arr)
!  (13) JTG     (INTEGER) : Cloud top layer for columns undergoing conv.
!  (14) JBG     (INTEGER) : Cloud bottom layer for columns undergoing conv.
!  (15) IDEEP   (INTEGER) : Index of longitudes where deep conv. happens
!  (16) LENGATH (INTEGER) : Length of gathered arrays
! 
!  NOTES:
!  (1 ) Removed NSTEP from arg list; it's not used.  Also zero arrays in order
!        to prevent them from being filled with compiler junk for latitudes
!        where no convection occurs at all. (bmy, 1/21/04)
!  (2 ) Now dimension DP, MU, MD, EU as (IIPAR,JJPAR,LLPAR) to avoid seg fault
!        error in OpenMP.  Also now pass the GEOS-CHEM latitude index J via
!        the argument list. (bmy, 12/13/05)
!******************************************************************************
!
#     include "CMN_SIZE"   ! Size parameters
      
      ! Arguments
      INTEGER, INTENT(IN)  :: J 
      REAL*8,  INTENT(IN)  :: DP(IIPAR,JJPAR,LLPAR) 
      REAL*8,  INTENT(IN)  :: MU(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)  :: MD(IIPAR,JJPAR,LLPAR) 
      REAL*8,  INTENT(IN)  :: EU(IIPAR,JJPAR,LLPAR) 
      REAL*8,  INTENT(OUT) :: MUG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: MDG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DUG(IIPAR,LLPAR)      
      REAL*8,  INTENT(OUT) :: EUG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: EDG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DPG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DSUBCLD(IIPAR)   
      INTEGER, INTENT(OUT) :: JTG(IIPAR)
      INTEGER, INTENT(OUT) :: JBG(IIPAR)
      INTEGER, INTENT(OUT) :: IDEEP(IIPAR)
      INTEGER, INTENT(OUT) :: LENGATH

      ! Local variables
      INTEGER              :: I, K, LENPOS 
      INTEGER              :: INDEX(IIPAR)
      REAL*8               :: SUM(IIPAR)
      REAL*8               :: RDPG(IIPAR,LLPAR)      

      !=================================================================
      ! ARCONVTRAN begins here!
      !=================================================================

      ! Initialize arrays
      DPG     = 0d0
      DSUBCLD = 0d0
      DUG     = 0d0
      EDG     = 0d0
      EUG     = 0d0
      JTG     = LLPAR
      JBG     = 1
      MDG     = 0d0
      MUG     = 0d0
      RDPG    = 0d0
      SUM     = 0d0
      
      !=================================================================
      ! First test if convection exists in the lon band I=1,IIPAR
      !=================================================================      

      ! Sum all upward mass fluxes in the longitude band I=1,IIPAR
      DO K = 1, LLPAR
      DO I = 1, IIPAR
         SUM(I) = SUM(I) + MU(I,J,K)
      ENDDO
      ENDDO

      ! IDEEP is the index of longitudes where SUM( up mass flux ) > 0
      ! LENGATH is the # of values where SUM( up mass flux ) > 0
      CALL WHENFGT( IIPAR, SUM, 1, 0d0, IDEEP, LENGATH )
      
      ! Return if there is no convection the longitude band
      IF ( LENGATH == 0 ) RETURN

      !=================================================================
      ! Gather input mass fluxes in places where there is convection
      !=================================================================
      DO K = 1, LLPAR
      DO I = 1, LENGATH

         ! Convert Pa->hPa
         DPG(I,K)  = 0.01d0 * DP(IDEEP(I),J,K)  
         RDPG(I,K) = 1.d0 / DPG(I,K)

         ! Convert Pa/s->hPa/s
         MUG(I,K)  = MU(IDEEP(I),J,K) *  0.01d0           
         MDG(I,K)  = MD(IDEEP(I),J,K) *  0.01d0

         ! Convert Pa/s->1/s
         EUG(I,K)  = EU(IDEEP(I),J,K) *  0.01d0 * RDPG(I,K) 
      ENDDO
      ENDDO

      !=================================================================
      ! Calc DU and ED in places where there is convection
      !=================================================================
      DO K = 1, LLPAR-1
      DO I = 1, LENGATH
         DUG(I,K) = EUG(I,K) - ( MUG(I,K) - MUG(I,K+1) ) * RDPG(I,K)
         EDG(I,K) = ( MDG(I,K) - MDG(I,K+1) ) * RDPG(I,K)
      ENDDO
      ENDDO

      DO I = 1, LENGATH
         DUG(I,LLPAR) = EUG(I,LLPAR) - MUG(I,LLPAR) * RDPG(I,LLPAR)
         EDG(I,LLPAR) = 0.0d0
      ENDDO

      DO K = 1, LLPAR
      DO I = 1, LENGATH
         IF ( DUG(I,K) < 1.d-7*EUG(I,K) ) DUG(I,K) = 0.0d0
      ENDDO
      ENDDO

      !=================================================================
      ! Find top and bottom layers with updrafts.
      !=================================================================
      DO I = 1, LENGATH
         JTG(I) = LLPAR
         JBG(I) = 1
      ENDDO

      ! Loop over altitudes
      DO K = 2, LLPAR
         
         ! Find places in the gathered array where MUG > 0
         CALL WHENFGT( LENGATH, MUG(:,K), 1, 0D0, INDEX, LENPOS )
         
         ! Compute top & bottom layers
         DO I = 1, LENPOS    
            JTG(INDEX(I)) = MIN( K-1, JTG(INDEX(I)) )
            JBG(INDEX(I)) = MAX( K,   JBG(INDEX(I)) )
         ENDDO
      ENDDO

      !=================================================================
      ! Calc delta p between srfc and cloud base.
      !=================================================================
      DO I = 1, LENGATH
         DSUBCLD(I) = DPG(I,LLPAR)
      ENDDO

      DO K = LLPAR-1, 2, -1
      DO I = 1, LENGATH
         IF ( JBG(I) <= K ) THEN
            DSUBCLD(I) = DSUBCLD(I) + DPG(I,K)
         ENDIF
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ARCONVTRAN

!------------------------------------------------------------------------------

      SUBROUTINE CONVTRAN( J,     NTRACE, Q,      MU,   MD,       
     &                     DU,    EU,     ED,     DP,   DSUBCLD,  
     &                     JT,    MX,     IDEEP,  IL1G, IL2G,     
     &                     NSTEP, DELT,   FRACIS, TCVV, INDEXSOL )
!
!******************************************************************************
!  Subroutine CONVTRAN applies the convective transport of trace species
!  (assuming moist mixing ratio) using the ZHANG/MCFARLANE convection scheme. 
!  (pjr, dsa, bmy, 6/26/03, 12/13/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J        (INTEGER) : GEOS-CHEM latitude index              [unitless]
!  (2 ) NTRACE   (INTEGER) : Number of tracers to transport        [unitless]
!  (3 ) Q        (REAL*8 ) : Tracer conc. including moisture       [v/v     ]
!  (4 ) MU       (REAL*8 ) : Mass flux up                          [hPa/s   ]
!  (5 ) MD       (REAL*8 ) : Mass flux down                        [hPa/s   ]
!  (6 ) DU       (REAL*8 ) : Mass detraining from updraft          [1/s     ]
!  (7 ) EU       (REAL*8 ) : Mass entraining from updraft          [1/s     ]
!  (8 ) ED       (REAL*8 ) : Mass entraining from downdraft        [1/s     ]
!  (9 ) DP       (REAL*8 ) : Delta pressure between interfaces
!  (10) DSUBCLD  (REAL*8 ) : Delta pressure from cloud base to sfc
!  (11) JT       (INTEGER) : Index of cloud top for each column
!  (12) MX       (INTEGER) : Index of cloud top for each column
!  (13) IDEEP    (INTEGER) : Gathering array
!  (14) IL1G     (INTEGER) : Gathered min lon indices over which to operate
!  (15) IL2G     (INTEGER) : Gathered max lon indices over which to operate
!  (16) NSTEP    (INTEGER) : Time step index
!  (17) DELT     (REAL*8 ) : Time step
!  (18) FRACIS   (REAL*8 ) : Fraction of tracer that is insoluble
!  (19) TCVV     (REAL*8 ) : Ratio of air mass / tracer mass 
!  (20) INDEXSOL (INTEGER) : Index array of soluble tracer numbers
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) Q       (REAL*8 ) : Contains modified tracer mixing ratios [v/v]
!
!  Important Local Variables:
!  ============================================================================
!  (1 ) CABV    (REAL*8 ) : Mixing ratio of constituent above
!  (2 ) CBEL    (REAL*8 ) : Mix ratio of constituent beloqw
!  (3 ) CDIFR   (REAL*8 ) : Normalized diff between cabv and cbel
!  (4 ) CHAT    (REAL*8 ) : Mix ratio in env at interfaces
!  (5 ) CMIX    (REAL*8 ) : Gathered tracer array 
!  (6 ) COND    (REAL*8 ) : Mix ratio in downdraft at interfaces
!  (7 ) CONU    (REAL*8 ) : Mix ratio in updraft at interfaces
!  (8 ) DCONDT  (REAL*8 ) : Gathered tend array 
!  (9 ) FISG    (REAL*8 ) : gathered insoluble fraction of tracer
!  (10) KBM     (INTEGER) : Highest altitude index of cloud base [unitless]
!  (11) KTM     (INTEGER) : Hightet altitude index of cloud top  [unitless]
!  (12) MBSTH   (REAL*8 ) : Threshold for mass fluxes
!  (13) SMALL   (REAL*8 ) : A small number
!
!  NOTES:
!  (1 ) Added references to "diag_mod.f", "grid_mod.f", and "CMN_DIAG.  
!        Also added TCVV and INDEXSOL as arguments.  Now only save LD38
!        levels of the ND38 diagnostic.  Now place NTRACE before Q in the
!        arg list. (swu, bmy, 1/21/04)
!  (2 ) Now pass Hg2 that is wet scavenged to "ocean_mercury_mod.f" for
!        computation of mercury fluxes (sas, bmy, 1/21/05)
!  (3 ) Now dimension Q and FRACIS of size (IIPAR,JJPAR,LLPAR,NTRACE), in 
!        order to avoid seg faults with OpenMP.  Also renamed GEOS-CHEM 
!        latitude index LATI_INDEX to J.  Now references ITS_A_MERCURY_SIM 
!        from "tracer_mod.f". Now references IS_Hg2 from "tracerid_mod.f.
!        Now do not call ADD_Hg2_WD if we are not using the dynamic online
!        ocean model.  Now references LDYNOCEAN from "logical_mod.f".
!        (cdh, bmy, 2/27/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,          ONLY : AD38 
      USE GRID_MOD,          ONLY : GET_AREA_M2
      USE LOGICAL_MOD,       ONLY : LDYNOCEAN
      USE OCEAN_MERCURY_MOD, ONLY : ADD_Hg2_WD
      USE TRACER_MOD,        ONLY : ITS_A_MERCURY_SIM
      USE TRACERID_MOD,      ONLY : IS_Hg2

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_DIAG"          ! ND38, LD38

      ! Arguments
      INTEGER, INTENT(IN)        :: J
      INTEGER, INTENT(IN)        :: NTRACE  
      REAL*8,  INTENT(INOUT)     :: Q(IIPAR,JJPAR,LLPAR,NTRACE)  
      REAL*8,  INTENT(IN)        :: MU(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)        :: MD(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)        :: DU(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)        :: EU(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)        :: ED(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)        :: DP(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)        :: DSUBCLD(IIPAR)      
      INTEGER, INTENT(IN)        :: JT(IIPAR)          
      INTEGER, INTENT(IN)        :: MX(IIPAR)          
      INTEGER, INTENT(IN)        :: IDEEP(IIPAR)       
      INTEGER, INTENT(IN)        :: IL1G               
      INTEGER, INTENT(IN)        :: IL2G 
      INTEGER, INTENT(IN)        :: NSTEP               
      REAL*8,  INTENT(IN)        :: DELT                
      REAL*8,  INTENT(IN)        :: FRACIS(IIPAR,JJPAR,LLPAR,NTRACE) 
      REAL*8,  INTENT(IN)        :: TCVV(NTRACE)
      INTEGER, INTENT(IN)        :: INDEXSOL(NTRACE)

      ! Local variables
      LOGICAL                    :: IS_Hg
      INTEGER                    :: II,      JJ,      LL,      NN
      INTEGER                    :: I,       K,       KBM,     KK     
      INTEGER                    :: KKP1,    KM1,     KP1,     KTM     
      INTEGER                    :: M,       ISTEP
      REAL*8                     :: CABV,    CBEL,    CDIFR,   CD2
      REAL*8                     :: DENOM,   SMALL,   MBSTH,   MUPDUDP
      REAL*8                     :: MINC,    MAXC,    QN,      FLUXIN
      REAL*8                     :: FLUXOUT, NETFLUX, AREA_M2, WET_Hg2     
      REAL*8                     :: CHAT(IIPAR,LLPAR)     
      REAL*8                     :: COND(IIPAR,LLPAR)     
      REAL*8                     :: CMIX(IIPAR,LLPAR)     
      REAL*8                     :: FISG(IIPAR,LLPAR)     
      REAL*8                     :: CONU(IIPAR,LLPAR)     
      REAL*8                     :: DCONDT(IIPAR,LLPAR)   

      !=================================================================
      ! CONVTRAN begins here!
      !=================================================================

      ! Is this a mercury simulation with dynamic ocean model?
      IS_Hg = ( ITS_A_MERCURY_SIM() .and. LDYNOCEAN )

      ! A small number
      SMALL = 1.d-36

      ! Threshold below which we treat the mass fluxes as zero (in mb/s)
      MBSTH = 1.d-15

      !=================================================================
      ! Find the highest level top and bottom levels of convection
      !=================================================================
      KTM = LLPAR
      KBM = LLPAR
      DO I = IL1G, IL2G
         KTM = MIN( KTM, JT(I) )
         KBM = MIN( KBM, MX(I) )
      ENDDO

      ! Loop ever each tracer
      DO M = 1, NTRACE

         ! Gather up the tracer and set tend to zero
         DO K = 1,    LLPAR
         DO I = IL1G, IL2G
            CMIX(I,K) = Q(IDEEP(I),J,K,M)
            IF ( CMIX(I,K) < 4.d0*SMALLEST ) CMIX(I,K) = 0D0
            FISG(I,K) = FRACIS(IDEEP(I),J,K,M)
         ENDDO
         ENDDO

         !==============================================================
         ! From now on work only with gathered data
         ! Interpolate environment tracer values to interfaces
         !==============================================================
         DO K = 1, LLPAR
            KM1 = MAX( 1, K-1 )

            DO I = IL1G, IL2G
               MINC = MIN( CMIX(I,KM1), CMIX(I,K) )
               MAXC = MAX( CMIX(I,KM1), CMIX(I,K) )

               IF ( MINC < 0d0 ) THEN 
                  CDIFR = 0.d0
               ELSE
                  CDIFR = ABS( CMIX(I,K)-CMIX(I,KM1) ) / MAX(MAXC,SMALL)
               ENDIF
               
               !------------------------------------------------------------
               ! The following 2 variables are actually NOT used
               ! (swu, 12/17/03)
               !DENOM = MAX( MAXC, SMALL )       
               !CD2   = ABS( CMIX(I,K) - CMIX(I,KM1) ) / DENOM
               !------------------------------------------------------------

               IF ( CDIFR > 1.d-6 ) THEN

                  ! If the two layers differ significantly.
                  ! use a geometric averaging procedure
                  CABV = MAX( CMIX(I,KM1), MAXC*TINYNUM, SMALLEST )
                  CBEL = MAX( CMIX(I,K),   MAXC*TINYNUM, SMALLEST )

                  CHAT(I,K) = LOG( CABV / CBEL)
     &                       /   ( CABV - CBEL)
     &                       *     CABV * CBEL

               ELSE             

                  ! Small diff, so just arithmetic mean
                  CHAT(I,K) = 0.5d0 * ( CMIX(I,K) + CMIX(I,KM1) )
               ENDIF

               ! Provisional up and down draft values
               CONU(I,K) = CHAT(I,K)
               COND(I,K) = CHAT(I,K)

               ! Provisional tends
               DCONDT(I,K) = 0.d0
            ENDDO
         ENDDO

         !==============================================================
         ! Do levels adjacent to top and bottom
         !==============================================================
         K   = 2
         KM1 = 1
         KK  = LLPAR 

         DO I = IL1G, IL2G
            MUPDUDP = MU(I,KK) + DU(I,KK) * DP(I,KK)

            IF ( MUPDUDP > MBSTH ) THEN
               CONU(I,KK) = ( EU(I,KK)*CMIX(I,KK)*DP(I,KK) ) 
     &                 /MUPDUDP  
            ENDIF

            IF ( MD(I,K) < -MBSTH ) THEN
               COND(I,K) =  (-ED(I,KM1)*CMIX(I,KM1)*DP(I,KM1))
     &                      /MD(I,K) 
            ENDIF
         ENDDO

         !==============================================================
         ! Updraft from bottom to top
         !==============================================================
         DO KK = LLPAR-1,1,-1
            KKP1 = MIN( LLPAR, KK+1 )

            DO I = IL1G,IL2G
               MUPDUDP = MU(I,KK) + DU(I,KK) * DP(I,KK)
                IF ( MUPDUDP > MBSTH ) THEN
                  CONU(I,KK) = (MU(I,KKP1)*CONU(I,KKP1) *FISG(I,KK)
     &                         +EU(I,KK)*CMIX(I,KK)*DP(I,KK)
     &                         )/MUPDUDP 
               ENDIF
            ENDDO
         ENDDO

         !==============================================================
         ! Downdraft from top to bottom
         !==============================================================
         DO K = 3, LLPAR
            KM1 = MAX( 1, K-1 )

            DO I = IL1G, IL2G
               IF ( MD(I,K) < -MBSTH ) THEN
                  COND(I,K) =  (  MD(I,KM1)*COND(I,KM1) 
     $                           -ED(I,KM1)*CMIX(I,KM1)
     $                         *DP(I,KM1))/MD(I,K)
               ENDIF
            ENDDO
         ENDDO

         DO K = KTM, LLPAR
            KM1 = MAX( 1,     K-1 )
            KP1 = MIN( LLPAR, K+1 )

            DO I = IL1G, IL2G

               ! Version 3 limit fluxes outside convection to mass in 
               ! appropriate layer.  These limiters are probably only safe
               ! for positive definite quantitities.  It assumes that mu 
               ! and md already satify a courant number limit of 1

               FLUXIN =  MU(I,KP1)* CONU(I,KP1) * FISG(I,K)
     $                + (MU(I,K)+MD(I,K)) * CMIX(I,KM1) 
     $                -  MD(I,K)  * COND(I,K)
   
               FLUXOUT = MU(I,K)   * CONU(I,K)     
     $                 +(MU(I,KP1)+MD(I,KP1)) * CMIX(I,K)
     $                 - MD(I,KP1) * COND(I,KP1) 

!------------------------------------------------------------------------------
! !!!!!!! backup: also works OK !!!!!!!!! (swu, 12/17/03)
!              FLUXIN =  MU(I,KP1)* CONU(I,KP1) 
!     $                +  MU(I,K)  * 0.5d0*(CHAT(I,K)+CMIX(I,KM1)) 
!     $                -  MD(I,K)  * COND(I,K)   
!     $                -  MD(I,KP1)* 0.5d0*(CHAT(I,KP1)+CMIX(I,KP1))
!
!               FLUXOUT = MU(I,K)   * CONU(I,K)     
!     $                 + MU(I,KP1) * 0.5d0*(CHAT(I,KP1)+CMIX(I,K))
!     $                 - MD(I,KP1) * COND(I,KP1) 
!     $                 - MD(I,K)   * 0.5d0*(CHAT(I,K)+CMIX(I,K))
!
!               FLUXIN =  MU(I,KP1)* CONU(I,KP1) 
!     $                +  MU(I,K)  * CHAT(I,K)
!     $                -  MD(I,K)  * COND(I,K)   
!     $                -  MD(I,KP1)* CHAT(I,KP1)
!
!               FLUXOUT = MU(I,K)   * CONU(I,K)     
!     $                 + MU(I,KP1) * CHAT(I,KP1)
!     $                 - MD(I,KP1) * COND(I,KP1) 
!     $                 - MD(I,K)   * CHAT(I,K)
!------------------------------------------------------------------------------

               !========================================================
               ! ND38 Diagnostic: loss of soluble tracer [kg/s] to
               ! convective rainout ("WETDCV-$") (swu, bmy, 12/17/03) 
               !========================================================

               ! Soluble tracer index
               NN = INDEXSOL(M)

               ! Only save to ND38 if it's turned on, if there are soluble 
               ! tracers, and if we are below the LD38 level limit
               IF ( ND38 > 0 .and. NN > 0 ) THEN

                  ! GEOS-CHEM lon, lat, alt indices
                  II = IDEEP(I)
                  JJ = J
                  LL = LLPAR - K + 1
                  
                  ! Only save up to LD38 vertical levels
                  IF ( LL <= LD38 ) THEN
                  
                     ! Grid box surface area [m2] 
                     AREA_M2 = GET_AREA_M2( JJ )  

                     ! Save loss in [kg/s]
                     AD38(II,JJ,LL,NN) = AD38(II,JJ,LL,NN) +
     &                    MU(I,KP1)   * AREA_M2     * 100d0         / 
     &                    GRAV        * CONU(I,KP1) * (1-FISG(I,K)) / 
     &                    TCVV(M)     / FLOAT(NSTEP)
                  ENDIF
               ENDIF

               !========================================================
               ! Pass the amount of Hg2 lost in wet scavenging [kg] 
               ! to "ocean_mercury_mod.f" w/ ADD_Hg2_WET.  
               !
               ! NOTE: DELT is already divided by NSTEP (as passed from
               ! the calling program) so we don't have to divide by
               ! it here, as is done above for ND38. (sas, bmy, 1/21/05)
               !
               ! ALSO NOTE: Do not do this computation if we are not
               ! using the online dynamic ocean (i.e. if LDYNOCEAN=F).
               ! (bmy, 2/27/06)
               !========================================================
               IF ( IS_Hg .and. IS_Hg2( M ) ) THEN

                  ! GEOS-CHEM lon & lat indices
                  II      = IDEEP(I)
                  JJ      = J

                  ! Grid box surface area [m2] 
                  AREA_M2 = GET_AREA_M2( JJ )

                  ! Hg2 wet-scavenged out of the column [kg]
                  WET_Hg2 = MU(I,KP1) * AREA_M2     * 100d0         / 
     &                      GRAV      * CONU(I,KP1) * (1-FISG(I,K)) /
     &                      TCVV(M)   * DELT 

                  ! Pass to "ocean_mercury_mod.f"
                  CALL ADD_Hg2_WD( II, J, M, WET_Hg2 )
               ENDIF

               NETFLUX = FLUXIN - FLUXOUT
               
               IF ( ABS(NETFLUX) < MAX(FLUXIN,FLUXOUT)*TINYNUM) THEN
                  NETFLUX = 0.D0
               ENDIF

               DCONDT(I,K) = NETFLUX / DP(I,K)
            ENDDO 
         ENDDO    

         DO K = KBM, LLPAR             
            KM1 = MAX( 1, K-1 )

            DO I = IL1G, IL2G

               IF ( K == MX(I) ) THEN

                  FLUXIN  =(MU(I,K)+MD(I,K))* CMIX(I,KM1)              
     $                    - MD(I,K)*COND(I,K)

                  FLUXOUT = MU(I,K)*CONU(I,K) 

!----------------------------------------------------------------------------
! !!!!!! BACK UP; also works well !!!!!!!! (swu, 12/17/03)
!                  FLUXIN  = MU(I,K)*0.5d0*(CHAT(I,K)+CMIX(I,KM1))
!     $                    - MD(I,K)*COND(I,K)
!
!                  FLUXOUT = MU(I,K)*CONU(I,K) 
!     $                    - MD(I,K)*0.5d0*(CHAT(I,K)+CMIX(I,K))
!----------------------------------------------------------------------------

                  NETFLUX = FLUXIN - FLUXOUT

                  IF (ABS(NETFLUX).LT.MAX(FLUXIN,FLUXOUT)*TINYNUM) THEN
                     NETFLUX = 0.d0
                  ENDIF

                  DCONDT(I,K) = NETFLUX / DP(I,K)

               ELSE IF ( K > MX(I) ) THEN

                  DCONDT(I,K) = 0.D0

               ENDIF

            ENDDO 
         ENDDO    

         !==============================================================
         ! Update and scatter data back to full arrays
         !==============================================================
         DO K = 1, LLPAR
            KP1 = MIN( LLPAR, K+1 )
            DO I = IL1G, IL2G    
            
               QN = CMIX(I,K) + DCONDT(I,K) * DELT 

               ! Do not make Q negative!!! (swu, 12/17/03)
               IF ( QN < 0d0 ) THEN
                  QN = 0d0
               ENDIF            

               Q(IDEEP(I),J,K,M) = QN
            ENDDO 
         ENDDO    

      ENDDO  !M ; End of tracer loop

      ! Return to calling program
      END SUBROUTINE CONVTRAN

!-----------------------------------------------------------------------------

      SUBROUTINE WHENFGT( N, ARRAY, INC, TARGET, INDEX, NVAL )
!
!******************************************************************************
!  Subroutine WHENFGT examines a 1-D vector and returns both an index array
!  of elements and the number of elements which are greater than a certain 
!  target value.  This routine came with the fvDAS convection code, we just
!  cleaned it up and added comments. (swu, bmy, 1/21/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N      (INTEGER) : Number of elements in ARRAY
!  (2 ) ARRAY  (REAL*8 ) : 1-D vector to be examined
!  (3 ) INC    (INTEGER) : Increment for stepping thru ARRAY
!  (4 ) TARGET (REAL*8 ) : Value that ARRAY will be tested against
!
!  Arguments as Output:
!  ============================================================================
!  (5 ) INDEX  (INTEGER) : Array  of places where ARRAY(I) > TARGET
!  (6 ) NVAL   (INTEGER) : Number of places where ARRAY(I) > TARGET
!
!  NOTES:
!  (1 ) Updated comments.  Now use F90 style declarations.  (bmy, 1/21/04)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)  :: N, INC
      REAL*8,  INTENT(IN)  :: ARRAY(N), TARGET
      INTEGER, INTENT(OUT) :: INDEX(N), NVAL

      ! Local variables
      INTEGER              :: I, INA

      !=================================================================
      ! WHENFGT begins here!
      !=================================================================

      ! Initialize
      INA      = 1
      NVAL     = 0
      INDEX(:) = 0

      ! Loop thru the array
      DO I = 1, N

         ! If the current element of ARRAY is greater than TARGET,
         ! then increment NVAL and save the element # in INDEX
         IF ( ARRAY(INA) > TARGET ) THEN
	    NVAL        = NVAL + 1
	    INDEX(NVAL) = I
         ENDIF

         ! Skip ahead by INC elements
         INA = INA + INC
      ENDDO

      ! Return to calling program
      END SUBROUTINE WHENFGT

!------------------------------------------------------------------------------

      ! End of module
      END MODULE FVDAS_CONVECT_MOD
