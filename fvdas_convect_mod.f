! $Id: fvdas_convect_mod.f,v 1.1 2003/06/30 20:26:09 bmy Exp $
      MODULE FVDAS_CONVECT_MOD
!
!******************************************************************************
!  Module FVDAS_CONVECT_MOD contains routines (originally from NCAR) which 
!  perform shallow and deep convection for the GEOS-4/fvDAS met fields.  
!  These routines account for shallow and deep convection, plus updrafts 
!  and downdrafts.  (pjr, dsa, bmy, 6/26/03)
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
!  (6 ) WHENFGT            : Test funtion
!
!  GEOS-CHEM modules referenced by fvdas_convect_mod.f:
!  ============================================================================
!  (1 ) pressure_mod.f     : Module containing routines to compute P(I,J,L)
!
!  NOTES:  
!******************************************************************************
!
      IMPLICIT NONE
      
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "fvdas_convect_mod.f"
      !=================================================================

      ! Declare everything PRIVATE ...
      PRIVATE
      
      ! ... except routines INIT_FVDAS_CONVECT and FVDAS_CONVECT
      PUBLIC :: INIT_FVDAS_CONVECT, FVDAS_CONVECT

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
!  Subroutine INIT_FVDAS_CONVECT initializes the HACK and ZHANG/MCFARLANE
!  convection schemes for GEOS-4/fvDAS met fields. (dsa, bmy, 6/26/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) HYPI    (REAL*8) : Reference pressures at interfaces [hPa] 
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE PRESSURE_MOD, ONLY : GET_PEDGE

#     include "CMN_SIZE"  ! Size parameters
      
      ! Local variables
      INTEGER            :: I, J, L
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

      ! Construct array of pressure edges at (I,J)
      DO L = 1, LLPAR+1
         HYPI(L) = GET_PEDGE(I,J,L)
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
 100  FORMAT( '       - fvDAS convection is capped at L = ', i3, 
     &       ', or approx ', f6.1, ' hPa' )

      ! Return to calling program
      END SUBROUTINE INIT_FVDAS_CONVECT

!------------------------------------------------------------------------------

      SUBROUTINE FVDAS_CONVECT( TDT, Q,  RPDEL, ETA, BETA,  NTRACE, 
     &                          MU,  MD, EU,    DP,  NSTEP, FRACIS )
!
!******************************************************************************
!  Subroutine FVDAS_CONVECT is the convection driver routine for GEOS-4/fvDAS
!  met fields.  It calls both HACK and ZHANG/MCFARLANE convection schemes.
!  (pjr, dsa, bmy, 6/26/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TDT    (REAL*8 ) : 2 * delta-T                          [s]
!  (2 ) Q      (REAL*8 ) : Array of transported tracers         [v/v]
!  (3 ) RPDEL  (REAL*8 ) : 1./pde                               [1/hPa]
!  (4 ) ETA    (REAL*8 ) : GMAO Hack convective mass flux       [kg/m2/s]
!  (5 ) BETA   (REAL*8 ) : GMAO Hack overshoot parameter        [unitless]
!  (6 ) NTRACE (INTEGER) : Number of tracers to transport       [unitless]
!  (7 ) MU     (REAL*8 ) : GMAO updraft mass flux   (ZMMU)      [ ]
!  (8 ) MD     (REAL*8 ) : GMAO downdraft mass flux (ZMMD)      [ ]
!  (9 ) EU     (REAL*8 ) : GMAO updraft entrainment (ZMEU)      [ ]
!  (10) DP     (REAL*8 ) : Delta-pressure between level edges   [hPa]
!  (11) NSTEP  (INTEGER) : Time step index                      [unitless]
!  (12) FRACIS (REAL*8 ) : Fraction of tracer that is insoluble [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) STT    (REAL*8 ) : Modified tracer array              [v/v]
! 
!  Important Local Variables:
!  ============================================================================
!  (1 ) IDEEP  (INTEGER)  : Gathering array
!  (2 ) IL1G   (INTEGER)  : Gathered min lon indices over which to operate
!  (3 ) IL2G   (INTEGER)  : Gathered max lon indices over which to operate
!  (4 ) JT     (INTEGER)  : Index of cloud top for each column
!  (5 ) LENGATH(INTEGER)  : ??       
!  (6 ) DSUBCLD(REAL*8 )  : Delta pressure from cloud base to sfc
!  (7 ) DPG    (REAL*8 )  : gathered .01*dp
!  (8 ) DU     (REAL*8 )  : Mass detraining from updraft
!  (9 ) ED     (REAL*8 )  : Mass entraining from downdraft
!  (10) EUG    (REAL*8 )  : gathered eu
!  (11) MUG    (REAL*8 )  : gathered mu
!  (12) MDG    (REAL*8 )  : gathered md
!  (13) MX     (INTEGER)  : Index of cloud top for each column
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Arguments
      INTEGER, INTENT(IN)    :: NSTEP, NTRACE             
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

      ! Local variables
      INTEGER                :: J
      INTEGER                :: JT(IIPAR)
      INTEGER                :: MX(IIPAR)
      INTEGER                :: IDEEP(IIPAR)
      INTEGER                :: IL1G=1
      INTEGER                :: IL2G=JJPAR
      INTEGER                :: LENGATH
      REAL*8                 :: QTMP(IIPAR,LLPAR,NTRACE)
      REAL*8                 :: DU(IIPAR,LLPAR)
      REAL*8                 :: ED(IIPAR,LLPAR)
      REAL*8                 :: DSUBCLD(IIPAR)
      REAL*8                 :: MUG(IIPAR,LLPAR)
      REAL*8                 :: MDG(IIPAR,LLPAR)
      REAL*8                 :: EUG(IIPAR,LLPAR)
      REAL*8                 :: DPG(IIPAR,LLPAR)
      
      !=================================================================
      ! FVDAS_CONVECT begins here!
      !=================================================================

      ! Loop over latitudes
      DO J = 1, JJPAR
         
         ! Save latitude slice of STT into Q
         QTMP(:,:,:) = Q(:,J,:,:)

         !----------------------------
         ! ZHANG/MCFARLANE convection
         !----------------------------
         CALL ARCONVTRAN( NSTEP,       DP(:,J,:),    MU(:,J,:),
     &                    MD(:,J,:),   EU(:,J,:),    MUG,   
     &                    MDG,         DU,           EUG, 
     &                    ED,          DPG,          DSUBCLD, 
     &                    JT,          MX,           IDEEP, 
     &                    LENGATH )

         CALL CONVTRAN(   QTMP,        NTRACE,       MU(:,J,:),  
     &                    MD(:,J,:),   DU,           EU(:,J,:),         
     &                    ED,          DP(:,J,:),    DSUBCLD,    
     &                    JT,          MX,           IDEEP, 
     &                    IL1G,        IL2G,         NSTEP, 
     &                    0.5D0*TDT,   FRACIS(:,J,:,:) )

         
         !----------------------------
         ! HACK convection
         !----------------------------
         CALL HACK_CONV(  TDT,         RPDEL(:,J,:), ETA(:,J,:), 
     &                    BETA(:,J,:), NTRACE,       QTMP )

         ! Store latitude slice back into STT
         Q(:,J,:,:) = QTMP(:,:,:)
      ENDDO

      ! Return to calling program
      END SUBROUTINE FVDAS_CONVECT

!------------------------------------------------------------------------------

      SUBROUTINE HACK_CONV( TDT, RPDEL, ETA, BETA, NTRACE, Q )
!
!******************************************************************************
!  Subroutine HACK_CONV computes the convective mass flux adjustment to all 
!  tracers using the convective mass fluxes and overshoot parameters for the 
!  Hack scheme. (pjr, dsa, bmy, 6/26/03)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) TDT    (REAL*8)  : 2 delta-t                              [s)
!  (2 ) RPDEL  (REAL*8)  : Reciprocal of pressure-thickness array [1/hPa]
!  (3 ) ETA    (REAL*8)  : GMAO Hack convective mass flux (HKETA) [kg/m2/s]
!  (4 ) BETA   (REAL*8)  : GMAO Hack overshoot parameter (HKBETA) [unitless]
!  (5 ) NTRACE (INTEGER) : Number of tracers in the Q array       [unitless]
!  (6 ) Q      (REAL*8)  : Tracer concentrations                  [v/v]    
!  
!  Arguments as Output:
!  ============================================================================
!  (6 ) Q     (REAL*8)   : Modified tracer concentrations         [v/v]       
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
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACE
      REAL*8,  INTENT(IN)    :: TDT
      REAL*8,  INTENT(IN)    :: RPDEL(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: ETA(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: BETA(IIPAR,LLPAR)
      REAL*8,  INTENT(INOUT) :: Q(IIPAR,LLPAR,NTRACE)

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
            IF ( ETA(I,K) /= 0.0 ) THEN
               ETAGDT(I)   = ETA(I,K) * GRAV * TDT
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
               IF ( ( Q(I,K+1,M) < 0.0 )  .OR. 
     &              ( Q(I,K,M)   < 0.0 )  .OR.
     &              ( Q(I,K-1,M) < 0.0 ) ) GOTO 40

               ! Specify constituent interface values (linear interpolation)
               CMRH(I,K  ) = 0.5d0 *( Q(I,K-1,M) + Q(I,K  ,M) )
               CMRH(I,K+1) = 0.5d0 *( Q(I,K  ,M) + Q(I,K+1,M) )
               
               CMRC(I) = Q(I,K+1,M)

               ! Determine fluxes, flux divergence => changes due to 
               ! convection.  Logic must be included to avoid producing 
               ! negative values. A bit messy since there are no a priori 
               ! assumptions about profiles.  Tendency is modified (reduced) 
               ! when pending disaster detected.
               BOTFLX   = ETAGDT(I)*(CMRC(I) - CMRH(I,K+1))*ADJFAC
               TOPFLX   = BETA(I,K)*ETAGDT(I)*(CMRC(I)-CMRH(I,K))*ADJFAC
               DCMR1(I) = -BOTFLX*RPDEL(I,K+1)
               EFAC1    = 1.0d0
               EFAC2    = 1.0d0
               EFAC3    = 1.0d0
               
               IF ( Q(I,K+1,M)+DCMR1(I) < 0.0 ) THEN
                  EFAC1 = MAX(TINYALT,ABS(Q(I,K+1,M)/DCMR1(I)) - EPS)
               ENDIF

               IF ( EFAC1 == TINYALT .OR. EFAC1 > 1.0 ) EFAC1 = 0.0D0
               DCMR1(I) = -EFAC1*BOTFLX*RPDEL(I,K+1)
               DCMR2(I) = (EFAC1*BOTFLX - TOPFLX)*RPDEL(I,K)
               
               IF ( Q(I,K,M)+DCMR2(I) < 0.0 ) THEN
                  EFAC2 = MAX(TINYALT,ABS(Q(I,K  ,M)/DCMR2(I)) - EPS)
               ENDIF
               
               IF ( EFAC2 == TINYALT .OR. EFAC2 > 1.0 ) EFAC2 = 0.0D0
               DCMR2(I) = (EFAC1*BOTFLX - EFAC2*TOPFLX)*RPDEL(I,K)
               DCMR3(I) = EFAC2*TOPFLX*RPDEL(I,K-1)

               IF ( Q(I,K-1,M)+DCMR3(I) < 0.0 ) THEN
                  EFAC3 = MAX(TINYALT,ABS(Q(I,K-1,M)/DCMR3(I)) - EPS)
               ENDIF

               IF ( EFAC3 == TINYALT .OR. EFAC3 > 1.0 ) EFAC3 = 0.0D0
               EFAC3    = MIN(EFAC2,EFAC3)
               DCMR2(I) = (EFAC1*BOTFLX - EFAC3*TOPFLX)*RPDEL(I,K)
               DCMR3(I) = EFAC3*TOPFLX*RPDEL(I,K-1)
               
               Q(I,K+1,M) = Q(I,K+1,M) + DCMR1(I)
               Q(I,K  ,M) = Q(I,K  ,M) + DCMR2(I)
               Q(I,K-1,M) = Q(I,K-1,M) + DCMR3(I)
 40         CONTINUE
 50      CONTINUE
 70   CONTINUE
      
      ! Return to calling program
      END SUBROUTINE HACK_CONV

!------------------------------------------------------------------------------

      SUBROUTINE ARCONVTRAN( NSTEP, DP,  MU,    MD, 
     &                       EU,    MUG, MDG,   DUG, 
     &                       EUG,   EDG, DPG,   DSUBCLD, 
     &                       JTG,   JBG, IDEEP, LENGATH )
!
!******************************************************************************
!  Subroutine ARCONVTRAN sets up the convective transport using archived mass
!  fluxes from the ZHANG/MCFARLANE convection scheme.  The setup involves:
!    (1) Gather mass flux arrays.
!    (2) Calc the mass fluxes that are determined by mass balance.
!    (3) Determine top and bottom of convection.
!  (pjr, dsa, bmy, 6/26/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NSTEP   (INTEGER) : Time step index
!  (2 ) DP      (REAL*8 ) : Delta pressure between interfaces [Pa] ?? [hPa]
!  (3 ) MU      (REAL*8 ) : Mass flux up                      [kg/m2/s]
!  (4 ) MD      (REAL*8 ) : Mass flux down                    [kg/m2/s]
!  (5 ) EU      (REAL*8 ) : Mass entraining from updraft      [1/s]
!
!  Arguments as Output:
!  ============================================================================
!  (6 ) MUG     (REAL*8 ) : Gathered mu
!  (7 ) MDG     (REAL*8 ) : Gathered md
!  (8 ) DUG     (REAL*8 ) : Mass detraining from updraft (gathered)
!  (9 ) EUG     (REAL*8 ) : Gathered eu
!  (10) EDG     (REAL*8 ) : Mass entraining from downdraft (gathered)
!  (11) DPG     (REAL*8 ) : Gathered .01*dp
!  (12) DSUBCLD (REAL*8 ) : Delta pressure from cloud base to sfc (gathered)
!  (13) JTG     (INTEGER) : ??
!  (14) JBG     (INTEGER) : ??
!  (15) IDEEP   (INTEGER) : ??
!  (16) LENGATH (INTEGER) : ??
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE" ! Size parameters
      
      ! Arguments
      INTEGER, INTENT(IN)  :: NSTEP
      INTEGER, INTENT(OUT) :: JTG(IIPAR)
      INTEGER, INTENT(OUT) :: JBG(IIPAR)
      INTEGER, INTENT(OUT) :: IDEEP(IIPAR)
      INTEGER, INTENT(OUT) :: LENGATH
      REAL*8,  INTENT(IN)  :: DP(IIPAR,LLPAR) 
      REAL*8,  INTENT(IN)  :: MU(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)  :: MD(IIPAR,LLPAR) 
      REAL*8,  INTENT(IN)  :: EU(IIPAR,LLPAR) 
      REAL*8,  INTENT(OUT) :: MUG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: MDG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DUG(IIPAR,LLPAR)      
      REAL*8,  INTENT(OUT) :: EUG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: EDG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DPG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DSUBCLD(IIPAR)   

      ! Local variables
      INTEGER              :: I, K, LENPOS 
      INTEGER              :: INDEX(IIPAR)
      REAL*8               :: SUM(IIPAR)
      REAL*8               :: RDPG(IIPAR,LLPAR)      

      !=================================================================
      ! ARCONVTRAN begins here!
      !=================================================================

      ! Gathered array contains all columns with a updraft.
      DO I = 1, IIPAR
         SUM(I) = 0.d0
      ENDDO

      DO K = 1, LLPAR
      DO I = 1, IIPAR
         SUM(I) = SUM(I) + MU(I,K)
      ENDDO
      ENDDO
      CALL WHENFGT( IIPAR, SUM, 1, 0D0, IDEEP, LENGATH )

      ! Return if LENGATH is zero
      IF ( LENGATH == 0 ) return

      !=================================================================
      ! Gather input mass fluxes
      !=================================================================
      DO K = 1, LLPAR
      DO I = 1, LENGATH
         DPG(I,K)  = 0.01d0 * DP(IDEEP(I),K)        ! convert Pa -> hPa
         RDPG(I,K) = 1.d0 / DPG(I,K)
         MUG(I,K)  = MU(IDEEP(I),K) * GRAV * 0.01d0 ! convert kg/m2/s -> hPa/s
         MDG(I,K)  = MD(IDEEP(I),K) * GRAV * 0.01d0
         EUG(I,K)  = EU(IDEEP(I),K)
      ENDDO
      ENDDO

      !=================================================================
      ! Calc DU and ED
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

      DO K = 2, LLPAR
         CALL WHENFGT( LENGATH, MUG(1,K), 1, 0D0, INDEX, LENPOS )
         DO I = 1, LENPOS
            JTG(INDEX(I)) = MIN( K-1, JTG(INDEX(I)) )
            JBG(INDEX(I)) = MAX( K, JBG(INDEX(I)) )
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

      SUBROUTINE CONVTRAN( Q,    NTRACE, MU,   MD,      DU,
     &                     EU,   ED,     DP,   DSUBCLD, JT,   
     &                     MX,   IDEEP,  IL1G, IL2G,    NSTEP,   
     &                     DELT, FRACIS )
!
!******************************************************************************
!  Subroutine CONVTRAN applies the convective transport of trace species
!  (assuming moist mixing ratio) using the ZHANG/MCFARLANE convection scheme. 
!  (pjr, dsa, bmy, 6/26/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) Q       (REAL*8 ) : Tracer concentrations including moisture [v/v]
!  (2 ) NTRACE  (INTEGER) : Number of tracers to transport           [unitless]
!  (3 ) MU      (REAL*8 ) : Mass flux up             
!  (4 ) MD      (REAL*8 ) : Mass flux down
!  (5 ) DU      (REAL*8 ) : Mass detraining from updraft
!  (6 ) EU      (REAL*8 ) : Mass entraining from updraft
!  (7 ) ED      (REAL*8 ) : Mass entraining from downdraft
!  (8 ) DP      (REAL*8 ) : Delta pressure between interfaces
!  (9 ) DSUBCLD (REAL*8 ) : Delta pressure from cloud base to sfc
!  (10) JT      (INTEGER) : Index of cloud top for each column
!  (11) MX      (INTEGER) : Index of cloud top for each column
!  (12) IDEEP   (INTEGER) : Gathering array
!  (13) IL1G    (INTEGER) : Gathered min lon indices over which to operate
!  (14) IL2G    (INTEGER) : Gathered max lon indices over which to operate
!  (15) NSTEP   (INTEGER) : Time step index
!  (16) DELT    (REAL*8 ) : Time step
!  (17) FRACIS  (REAL*8 ) : Fraction of tracer that is insoluble
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) Q       (REAL*8 ) : Contains modified tracer mixing ratios [v/v]
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
!******************************************************************************
!
#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACE             
      INTEGER, INTENT(IN)    :: JT(IIPAR)          
      INTEGER, INTENT(IN)    :: MX(IIPAR)          
      INTEGER, INTENT(IN)    :: IDEEP(IIPAR)       
      INTEGER, INTENT(IN)    :: IL1G               
      INTEGER, INTENT(IN)    :: IL2G               
      INTEGER, INTENT(IN)    :: NSTEP               
      REAL*8,  INTENT(INOUT) :: Q(IIPAR,LLPAR,NTRACE)  
      REAL*8,  INTENT(IN)    :: MU(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: MD(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: DU(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: EU(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: ED(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: DP(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: DSUBCLD(IIPAR)      
      REAL*8,  INTENT(IN)    :: DELT                
      REAL*8,  INTENT(IN)    :: FRACIS(IIPAR,LLPAR,NTRACE) 

      ! Local variables
      INTEGER                :: I,     K,      KBM,     KK,     KKP1
      INTEGER                :: KM1,   KP1,    KTM,     M
      REAL*8                 :: CABV,  CBEL,   CDIFR,   CD2,    DENOM
      REAL*8                 :: SMALL, MBSTH,  MUPDUDP, MINC,   MAXC
      REAL*8                 :: QN,    FLUXIN, FLUXOUT, NETFLUX             
      REAL*8                 :: CHAT(IIPAR,LLPAR)     
      REAL*8                 :: COND(IIPAR,LLPAR)     
      REAL*8                 :: CMIX(IIPAR,LLPAR)     
      REAL*8                 :: FISG(IIPAR,LLPAR)     
      REAL*8                 :: CONU(IIPAR,LLPAR)     
      REAL*8                 :: DCONDT(IIPAR,LLPAR)   

      !=================================================================
      ! CONVTRAN begins here!
      !=================================================================

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
            CMIX(I,K) = Q(IDEEP(I),K,M)
            IF ( CMIX(I,K) < 4.*SMALLEST ) CMIX(I,K) = 0D0
            FISG(I,K) = FRACIS(IDEEP(I),K,M)
         ENDDO
         ENDDO

         !==============================================================
         ! From now on work only with gathered data
         ! Interpolate environment tracer values to interfaces
         !==============================================================
         DO K = 1, LLPAR
            KM1 = MAX(1,K-1)

            DO I = IL1G, IL2G
               MINC = MIN( CMIX(I,KM1), CMIX(I,K) )
               MAXC = MAX( CMIX(I,KM1), CMIX(I,K) )

               IF ( MINC < 0 ) THEN
                  CDIFR = 0.D0
               ELSE
                  CDIFR = ABS( CMIX(I,K)-CMIX(I,KM1) ) / MAX(MAXC,SMALL)
               ENDIF

               DENOM = MAX( MAXC, SMALL )
               CD2   = ABS( CMIX(I,K) - CMIX(I,KM1) ) / DENOM

               IF ( CDIFR > 1.D-6 ) THEN

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
               CONU(I,KK) = ( 
     &                       +EU(I,KK)*FISG(I,KK)*CMIX(I,KK)*DP(I,KK)
     &                      )/MUPDUDP
            ENDIF

            IF ( MD(I,K) < -MBSTH ) THEN
               COND(I,K) =  (  
     &                     -ED(I,KM1)*FISG(I,KM1)*CMIX(I,KM1)*DP(I,KM1)
     &                      )/MD(I,K)
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
                  CONU(I,KK) = (  MU(I,KKP1)*CONU(I,KKP1) 
     &                         +EU(I,KK)*FISG(I,KK)*CMIX(I,KK)*DP(I,KK)
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
     $                           -ED(I,KM1)*FISG(I,KM1)*CMIX(I,KM1)
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
               FLUXIN =  MU(I,KP1)*CONU(I,KP1) 
     $                +  MU(I,K)*MIN(CHAT(I,K),CMIX(I,KM1)) 
     $                - (MD(I,K)  *COND(I,K)   
     $                +  MD(I,KP1)*MIN(CHAT(I,KP1),CMIX(I,KP1)))

               FLUXOUT = MU(I,K)*CONU(I,K)     
     $                 + MU(I,KP1)*MIN(CHAT(I,KP1),CMIX(I,K))
     $                 -(MD(I,KP1)*COND(I,KP1) 
     $                 + MD(I,K)*MIN(CHAT(I,K),CMIX(I,K)))

               NETFLUX = FLUXIN - FLUXOUT

               IF ( ABS(NETFLUX) < MAX(FLUXIN,FLUXOUT)*TINYNUM) THEN
                  NETFLUX = 0.D0
               ENDIF

               DCONDT(I,K) = NETFLUX / DP(I,K)
            enddo
         enddo

         DO K = KBM, LLPAR             
            KM1 = MAX( 1, K-1 )

            DO I = IL1G, IL2G

               IF ( K == MX(I) ) THEN

                  FLUXIN  = MU(I,K)*MIN(CHAT(I,K),CMIX(I,KM1))
     $                    - MD(I,K)*COND(I,K)

                  FLUXOUT = MU(I,K)*CONU(I,K) 
     $                    - MD(I,K)*MIN(CHAT(I,K),CMIX(I,K))

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
               QN              = CMIX(I,K) + DCONDT(I,K) * 2.D0 * DELT
               Q(IDEEP(I),K,M) = QN
            ENDDO
         ENDDO

      ! End of tracer loop
      ENDDO

      ! Return to calling program
      END SUBROUTINE CONVTRAN

!-----------------------------------------------------------------------------

      SUBROUTINE WHENFGT( N, ARRAY, INC, TARGET, INDEX, NVAL )
!
!******************************************************************************
!  Subroutine WHENFGT is a
!
!  Arguments as Input:
!  ============================================================================
!  
!******************************************************************************
!
      ! Arguments
      INTEGER :: INDEX(*), NVAL, INC, N
      REAL*8  :: ARRAY(*), TARGET

      ! Local variables
      INTEGER :: I, INA

      !=================================================================
      ! WHENFGT begins here!
      !=================================================================
      INA  = 1
      NVAL = 0

      IF ( INC < 0 ) INA = (-INC)*(N-1)+1

      DO I = 1, N
         IF ( ARRAY(INA) > TARGET ) THEN
	    NVAL        = NVAL+1
	    INDEX(NVAL) = I
         ENDIF
         INA = INA + INC
      ENDDO

      ! Return to calling program
      END SUBROUTINE WHENFGT

!------------------------------------------------------------------------------

      END MODULE FVDAS_CONVECT_MOD
