! $Id: gcap_convect_mod.f,v 1.1 2009/09/16 14:06:26 bmy Exp $
      MODULE GCAP_CONVECT_MOD
!
!******************************************************************************
!  Module GCAP_CONVECT_MOD contains routines (originally from GISS) which
!  perform shallow and deep convection for the GCAP met fields.  This module
!  was based on FVDAS_CONVECT_MOD. (swu, bmy, 6/9/05, 12/19/06)
!  
!  Module Variables:
!  ============================================================================
!  (1 ) GRAV     (REAL*8 ) : Gravitational constant [m/s2]
!  (2 ) SMALLEST (REAL*8 ) : A very small double-precision number 
!  (3 ) TINYNUM  (REAL*8 ) : 2 times the SMALLEST
!
!  Module Routines:
!  ============================================================================
!  (1 ) INIT_GCAP_CONVECT  : Initializes fvDAS convection scheme
!  (2 ) GCAP_CONVECT       : GCAP/GISS convection driver
!  (4 ) ARCCONVTRAN        : Sets up fields for ZHANG/MCFARLANE convection
!  (5 ) CONVTRAN           : ZHANG/MCFARLANE convection scheme routine
!  (6 ) WHENFGT            : Test funtion
!
!  GEOS-CHEM modules referenced by fvdas_convect_mod.f:
!  ============================================================================
!  (1 ) pressure_mod.f     : Module containing routines to compute P(I,J,L)
!
!  NOTES:  
!  (1 ) Rewrote parallel loops to avoid problems w/ OpenMP (bmy, 12/13/05)
!  (2 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (3 ) More bug fixes for SUN 4100 platform.  Make SMALLEST = 1d-60 to avoid 
!        problems (bmy, 12/19/06)
!  02 Mar 2012 - R. Yantosca - Now reference the new grid_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!  
      IMPLICIT NONE
      
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "gcap_convect_mod.f"
      !=================================================================

      ! Declare everything PRIVATE ...
      PRIVATE
      
      ! ... except these routines
      PUBLIC :: GCAP_CONVECT

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Constants
      REAL*8,  PARAMETER :: GRAV     = 9.8d0
      REAL*8,  PARAMETER :: SMALLEST = 1d-60
      REAL*8,  PARAMETER :: TINYNUM  = 2*SMALLEST

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GCAP_CONVECT( TDT,    Q,      NTRACE,   DP,   
     &                         NSTEP,  FRACIS, TCVV,     INDEXSOL, 
     &                         UPDE,   DNDE,   ENTRAIN,  DETRAINE, 
     &                         UPDN,   DNDN,   DETRAINN ) 
!
!******************************************************************************
!  Subroutine GCAP_CONVECT is the convection driver routine for GEOS-4/fvDAS
!  met fields.  It calls the ZHANG/MCFARLANE convection scheme.
!  (swu, bmy, 6/9/05, 12/19/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TDT    (REAL*8 ) : 2 * delta-T                          [s]
!  (2 ) Q      (REAL*8 ) : Array of transported tracers         [v/v]
!  (3 ) RPDEL  (REAL*8 ) : 1./pde                               [1/hPa]
!  (4 ) ETA    (REAL*8 ) : GMAO Hack convective mass flux       [kg/m2/s]
!  (5 ) BETA   (REAL*8 ) : GMAO Hack overshoot parameter        [unitless]
!  (6 ) NTRACE (INTEGER) : Number of tracers to transport       [unitless]
!  (7 ) MU     (REAL*8 ) : GMAO updraft mass flux   (ZMMU)      [ ]pa/s
!  (8 ) MD     (REAL*8 ) : GMAO downdraft mass flux (ZMMD)      [ ]pa/s
!  (9 ) EU     (REAL*8 ) : GMAO updraft entrainment (ZMEU)      [ ]pa/s
!  (10) DP     (REAL*8 ) : Delta-pressure between level edges   [hPa]pa
!  (11) NSTEP  (INTEGER) : Time step index                      [unitless]
!  (12) FRACIS (REAL*8 ) : Fraction of tracer that is insoluble [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) Q      (REAL*8 ) : Modified tracer array              [v/v]
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
!  (8 ) MX     (INTEGER)  : Index of cloud top for each column
!
!  NOTES:
!  (1 ) Rewrote parallel loops so that we pass entire arrays to the various
!        subroutines instead of array slices such as (:,J,:).  This can cause
!        problems with OpenMP for some compilers. (bmy, 12/13/05)
!  (2 ) Now don't call CONVTRAN if LENGATH=0 (bmy, 12/19/06)
!  04 Nov 2013 - M. Sulprizio- Eliminate array temporaries by accepting
!                              pointer-valued arguments via assumed-shape dummy
!                              arguments
!******************************************************************************
!

      USE CMN_SIZE_MOD      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACE 
      REAL*8,  INTENT(IN)    :: TDT                
      REAL*8,  INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NTRACE)          
      REAL*8,  INTENT(IN)    :: DP(IIPAR,JJPAR,LLPAR)     
      INTEGER, INTENT(IN)    :: NSTEP
      REAL*8,  INTENT(IN)    :: FRACIS(IIPAR,JJPAR,LLPAR,NTRACE) 
      REAL*8,  INTENT(IN)    :: TCVV(NTRACE)
      INTEGER, INTENT(IN)    :: INDEXSOL(NTRACE) 
!------------------------------------------------------------------------------
! Prior to 11/4/13:
! Eliminate array temporary by accepting pointer variables into this routine
! via assumed-shape arguments. (mpayer, 11/4/13)
!      REAL*8,  INTENT(IN)    :: UPDE(IIPAR,JJPAR,LLPAR)     
!      REAL*8,  INTENT(IN)    :: DNDE(IIPAR,JJPAR,LLPAR)     
!      REAL*8,  INTENT(IN)    :: ENTRAIN(IIPAR,JJPAR,LLPAR)
!      REAL*8,  INTENT(IN)    :: DETRAINE(IIPAR,JJPAR,LLPAR)
!      REAL*8,  INTENT(IN)    :: UPDN(IIPAR,JJPAR,LLPAR)     
!      REAL*8,  INTENT(IN)    :: DNDN(IIPAR,JJPAR,LLPAR)     
!      REAL*8,  INTENT(IN)    :: DETRAINN(IIPAR,JJPAR,LLPAR)
!------------------------------------------------------------------------------
      REAL*8,  INTENT(IN)    :: UPDE(:,:,:)
      REAL*8,  INTENT(IN)    :: DNDE(:,:,:)
      REAL*8,  INTENT(IN)    :: ENTRAIN(:,:,:)
      REAL*8,  INTENT(IN)    :: DETRAINE(:,:,:)
      REAL*8,  INTENT(IN)    :: UPDN(:,:,:)
      REAL*8,  INTENT(IN)    :: DNDN(:,:,:)
      REAL*8,  INTENT(IN)    :: DETRAINN(:,:,:)

      ! Local variables
      INTEGER                :: I, J, L, N, LENGATH, ISTEP
      INTEGER                :: JT(IIPAR)
      INTEGER                :: MX(IIPAR)
      INTEGER                :: IDEEP(IIPAR)
      INTEGER                :: IL1G=1
      INTEGER                :: IL2G
      REAL*8                 :: DPG(IIPAR,LLPAR)
      REAL*8                 :: ED(IIPAR,LLPAR)
      REAL*8                 :: UPDEG(IIPAR,LLPAR)
      REAL*8                 :: DNDEG(IIPAR,LLPAR)
      REAL*8                 :: ENTRAING(IIPAR,LLPAR)
      REAL*8                 :: DETRAINEG(IIPAR,LLPAR)
      REAL*8                 :: TOTALDNDEG(IIPAR,LLPAR)
      REAL*8                 :: UPDNG(IIPAR,LLPAR)
      REAL*8                 :: DNDNG(IIPAR,LLPAR)
      REAL*8                 :: DETRAINNG(IIPAR,LLPAR)
      REAL*8                 :: TOTALDNDNG(IIPAR,LLPAR)
      REAL*8                 :: ENTRAINN(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: ENTRAINNG(IIPAR,LLPAR)

      !=================================================================
      ! GCAP_CONVECT begins here!
      !=================================================================

      ! Fake entrainment in non-entraining plumes (swu, bmy, 6/9/05)
      ENTRAINN = 0d0
      IL2G=JJPAR

      ! Loop over latitudes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J,         ISTEP,      UPDEG,  DNDEG, DETRAINEG )
!$OMP+PRIVATE( ENTRAING,  TOTALDNDEG, DPG,    JT,    MX        )
!$OMP+PRIVATE( IDEEP,     LENGATH,    UPDNG,  DNDNG, DETRAINNG )
!$OMP+PRIVATE( ENTRAINNG, TOTALDNDNG                           )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         
         !----------------------------
         ! Entraining convection
         !---------------------------- 

         ! Set up convection fields
         CALL ARCONVTRAN(  J,          NSTEP,     DP,        UPDE,  
     &                     DNDE,       DETRAINE,  ENTRAIN,   UPDEG, 
     &                     DNDEG,      DETRAINEG, ENTRAING,  TOTALDNDEG, 
     &                     DPG,        JT,        MX,        IDEEP,    
     &                     LENGATH )

         ! Internal convection steps
         DO ISTEP = 1, NSTEP 

            ! If there are nonzero fields at this J, do the convection
            IF ( LENGATH > 0 ) THEN
               CALL CONVTRAN( J,          NTRACE, Q,         
     &                        UPDEG,      DNDEG,  DETRAINEG, ENTRAING,  
     &                        TOTALDNDEG, DPG,    JT,        MX,        
     &                        IDEEP,      1,      LENGATH,   NSTEP,     
     &                        0.5D0*TDT,  FRACIS, TCVV,      INDEXSOL ) 
            ENDIF

         ENDDO 

         !----------------------------
         ! Non-entraining convection
         !---------------------------- 

         ! Set up convection fields
         CALL ARCONVTRAN(  J,          NSTEP,     DP,        UPDN,  
     &                     DNDN,       DETRAINN,  ENTRAINN,  UPDNG, 
     &                     DNDNG,      DETRAINNG, ENTRAINNG, TOTALDNDNG, 
     &                     DPG,        JT,        MX,        IDEEP,  
     &                     LENGATH )

         ! Loop over internal convection timesteps
         DO ISTEP = 1, NSTEP  

            !  If there are nonzero fields at this J, do the convection
            IF ( LENGATH > 0 ) THEN
               CALL CONVTRAN( J,          NTRACE, Q,          
     &                        UPDNG,      DNDNG,  DETRAINNG, ENTRAINNG, 
     &                        TOTALDNDNG, DPG,    JT,        MX,         
     &                        IDEEP,      1,      LENGATH,   NSTEP, 
     &                        0.5D0*TDT,  FRACIS, TCVV,      INDEXSOL ) 
            ENDIF

         ENDDO

      ENDDO 
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE GCAP_CONVECT

!------------------------------------------------------------------------------

      SUBROUTINE ARCONVTRAN( J,     NSTEP,    DP,  MU,  MD, 
     &                       DU,    EU,       MUG, MDG, DUG,   
     &                       EUG,   TOTALMDG, DPG, JTG, JBG, 
     &                       IDEEP, LENGATH )
!
!******************************************************************************
!  Subroutine ARCONVTRAN sets up the convective transport using archived mass
!  fluxes from the ZHANG/MCFARLANE convection scheme.  The setup involves:
!    (1) Gather mass flux arrays.
!    (2) Calc the mass fluxes that are determined by mass balance.
!    (3) Determine top and bottom of convection.
!  (pjr, dsa, bmy, 6/26/03, 12/13/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J       (INTEGER) : GEOS-CHEM latitude index          [unitless]
!  (2 ) NSTEP   (INTEGER) : Time step index
!  (3 ) DP      (REAL*8 ) : Delta pressure between interfaces [Pa      ]     
!  (4 ) MU      (REAL*8 ) : Mass flux up                      [kg/m2/s ]
!  (5 ) MD      (REAL*8 ) : Mass flux down                    [kg/m2/s ]
!  (6 ) EU      (REAL*8 ) : Mass entraining from updraft      [1/s     ]    
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) MUG     (REAL*8 ) : Gathered mu                                Pa/s
!  (8 ) MDG     (REAL*8 ) : Gathered md                                Pa/s
!  (9 ) DUG     (REAL*8 ) : Mass detraining from updraft (gathered)    Pa/S
!  (10) EUG     (REAL*8 ) : Gathered eu                                Pa/S
!  (11) EDG     (REAL*8 ) : Mass entraining from downdraft (gathered)  Pa/s
!  (12) DPG     (REAL*8 ) : Gathered                                   Pa
!  (13) DSUBCLD (REAL*8 ) : Delta pressure from cloud base to sfc (gathered)
!  (14) JTG     (INTEGER) : Cloud top layer for columns undergoing conv.
!  (15) JBG     (INTEGER) : Cloud bottom layer for columns undergoing conv.
!  (16) IDEEP   (INTEGER) : Index of longitudes where deep conv. happens
!  (17) LENGATH (INTEGER) : Length of gathered arrays
!
!  NOTES:
!  (1 ) Now dimension DP, MU, MD, DU, EU as (IIPAR,JJPAR,LLPAR) to avoid seg 
!        fault error in OpenMP.  Also now pass the GEOS-CHEM latitude index J 
!        via the argument list.  Add comments. (bmy, 12/13/05)
!  04 Nov 2013 - M. Sulprizio- Eliminate array temporaries by accepting
!                              pointer-valued arguments via assumed-shape dummy
!                              arguments
!******************************************************************************
!
      USE CMN_SIZE_MOD    ! Size parameters
      
      ! Arguments
      INTEGER, INTENT(IN)  :: J
      INTEGER, INTENT(IN)  :: NSTEP
      REAL*8,  INTENT(IN)  :: DP(IIPAR,JJPAR,LLPAR) 
!------------------------------------------------------------------------------
! Prior to 11/4/13:
! Eliminate array temporary by accepting pointer variables into this routine
! via assumed-shape arguments. (mpayer, 11/4/13)
!      REAL*8,  INTENT(IN)  :: MU(IIPAR,JJPAR,LLPAR)
!      REAL*8,  INTENT(IN)  :: MD(IIPAR,JJPAR,LLPAR)
!      REAL*8,  INTENT(IN)  :: DU(IIPAR,JJPAR,LLPAR)
!      REAL*8,  INTENT(IN)  :: EU(IIPAR,JJPAR,LLPAR) 
!------------------------------------------------------------------------------
      REAL*8,  INTENT(IN)  :: MU(:,:,:)
      REAL*8,  INTENT(IN)  :: MD(:,:,:)
      REAL*8,  INTENT(IN)  :: DU(:,:,:)
      REAL*8,  INTENT(IN)  :: EU(:,:,:)
      REAL*8,  INTENT(OUT) :: MUG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: MDG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DUG(IIPAR,LLPAR)      
      REAL*8,  INTENT(OUT) :: EUG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: TOTALMDG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DPG(IIPAR,LLPAR)
      INTEGER, INTENT(OUT) :: JTG(IIPAR)
      INTEGER, INTENT(OUT) :: JBG(IIPAR)
      INTEGER, INTENT(OUT) :: IDEEP(IIPAR)
      INTEGER, INTENT(OUT) :: LENGATH

      ! Local variables
      INTEGER              :: I, K, LENPOS 
      INTEGER              :: INDEX(IIPAR)
      REAL*8               :: SUM(IIPAR)
      REAL*8               :: RDPG(IIPAR,LLPAR)      
      REAL*8               :: TOTALMD(IIPAR,LLPAR)

      !=================================================================
      ! ARCONVTRAN begins here!
      !=================================================================

      ! Gathered array contains all columns with a updraft.
      DO I = 1, IIPAR
         SUM(I) = 0.d0
      ENDDO

      DO K = 1, LLPAR
      DO I = 1, IIPAR
         SUM(I) = SUM(I) + MU(I,J,K)

         ! Calculate totalMD --- all the downdrafts coming downstairs
         IF ( K == 1 ) THEN 
            TOTALMD(I,K) = MD(I,J,K)
         ELSE 
            TOTALMD(I,K) = TOTALMD(I,K-1) + MD(I,J,K)
         ENDIF

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
         DPG(I,K)      = DP(IDEEP(I),J,K)    !Pa
         MUG(I,K)      = MU(IDEEP(I),J,K)    !Pa/s
         MDG(I,K)      = MD(IDEEP(I),J,K) 
         EUG(I,K)      = EU(IDEEP(I),J,K)    
         DUG(I,K)      = DU(IDEEP(I),J,K) 
         TOTALMDG(I,K) = TOTALMD(IDEEP(I),K) !!!=sum( MD(ideep(I),1:K) )
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
         
         CALL WHENFGT( LENGATH, MUG(:,K), 1, 0D0, INDEX, LENPOS )
         
         DO I = 1, LENPOS    
            JTG(INDEX(I)) = MIN( K-1, JTG(INDEX(I)) )
            JBG(INDEX(I)) = MAX( K,   JBG(INDEX(I)) )
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ARCONVTRAN

!------------------------------------------------------------------------------

      SUBROUTINE CONVTRAN( J,    NTRACE, Q,       MU,     MD,       
     &                     DU,   EU,     TOTALMD, DP,     JT,   
     &                     MX,   IDEEP,  IL1G,    IL2G,   NSTEP,         
     &                     DELT, FRACIS, TCVV,    INDEXSOL ) 
!
!******************************************************************************
!  Subroutine CONVTRAN applies the convective transport of trace species
!  (assuming moist mixing ratio) using the ZHANG/MCFARLANE convection scheme. 
!  (swu, bmy, 6/9/05, 12/19/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J        (INTEGER) : GEOS-CHEM latitude index
!  (2 ) NTRACE   (INTEGER) : Number of tracers to transport         [unitless]
!  (3 ) Q        (REAL*8 ) : Tracer conc. including moisture        [v/v     ]
!  (4 ) MU       (REAL*8 ) : Mass flux up                           [Pa/s    ]
!  (5 ) MD       (REAL*8 ) : Mass flux down                         [Pa/s    ]
!  (6 ) DU       (REAL*8 ) : Mass detraining from updraft           [Pa/s    ]
!  (7 ) EU       (REAL*8 ) : Mass entraining from updraft           [Pa/s    ]
!  (8 ) ED       (REAL*8 ) : Mass entraining from downdraft         [Pa/s    ]
!  (9 ) DP       (REAL*8 ) : Delta pressure between interfaces      [Pa      ]
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
!  (1 ) Now dimension Q and FRACIS of size (IIPAR,JJPAR,LLPAR,NTRACE), in 
!        order to avoid seg faults with OpenMP.  Also renamed GEOS-CHEM 
!        latitude index LATI_INDEX to J.  Added comments. (bmy, 12/13/05)
!  (2 ) Bug fix: avoid div by zero in formula for CHAT (bmy, 12/19/06)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD38, CONVFLUP
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE PRESSURE_MOD, ONLY : GET_PEDGE

      USE CMN_SIZE_MOD     ! Size parameters
      USE CMN_DIAG_MOD 
  
      ! Arguments
      INTEGER, INTENT(IN)    :: J
      INTEGER, INTENT(IN)    :: NTRACE             
      REAL*8,  INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NTRACE)  
      REAL*8,  INTENT(IN)    :: MU(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: MD(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: DU(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: EU(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: TOTALMD(IIPAR,LLPAR)      
      REAL*8,  INTENT(IN)    :: DP(IIPAR,LLPAR)      
      INTEGER, INTENT(IN)    :: JT(IIPAR)          
      INTEGER, INTENT(IN)    :: MX(IIPAR)          
      INTEGER, INTENT(IN)    :: IDEEP(IIPAR)       
      INTEGER, INTENT(IN)    :: IL1G               
      INTEGER, INTENT(IN)    :: IL2G               
      INTEGER, INTENT(IN)    :: NSTEP               
      REAL*8,  INTENT(IN)    :: DELT                
      REAL*8,  INTENT(IN)    :: FRACIS(IIPAR,JJPAR,LLPAR,NTRACE) 
      REAL*8,  INTENT(IN)    :: TCVV(NTRACE)
      INTEGER, INTENT(IN)    :: INDEXSOL(NTRACE)

      ! Local variables
      INTEGER                :: I,     K,      KBM,     KK,     KKP1
      INTEGER                :: KM1,   KP1,    KTM,     M,      ISTEP
      INTEGER                :: II,    JJ,     LL,      NN
      REAL*8                 :: CABV,  CBEL,   CDIFR,   CD2,    DENOM
      REAL*8                 :: SMALL, MBSTH,  MUPDUDP, MINC,   MAXC
      REAL*8                 :: QN,    FLUXIN, FLUXOUT, NETFLUX             
      REAL*8                 :: CHAT(IIPAR,LLPAR)     
      REAL*8                 :: COND(IIPAR,LLPAR)     
      REAL*8                 :: CMIX(IIPAR,LLPAR)     
      REAL*8                 :: FISG(IIPAR,LLPAR)     
      REAL*8                 :: CONU(IIPAR,LLPAR)     
      REAL*8                 :: DCONDT(IIPAR,LLPAR)   
      REAL*8                 :: AREA_M2,           DELTAP
      REAL*8                 :: TRC_BFCONVTRAN,    TRC_AFCONVTRAN
      REAL*8                 :: PLUMEIN, PLUMEOUT, PLUMECHANGE

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
            CMIX(I,K) = Q(IDEEP(I),J,K,M)
            FISG(I,K) = FRACIS(IDEEP(I),J,K,M)
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

               IF ( CDIFR > 1.D-6 ) THEN

                  ! If the two layers differ significantly.
                  ! use a geometric averaging procedure
                  CABV = MAX( CMIX(I,KM1), MAXC*TINYNUM, SMALLEST )
                  CBEL = MAX( CMIX(I,K),   MAXC*TINYNUM, SMALLEST )

                  ! If CABV-CBEL is zero then set CHAT=SMALLEST
                  ! so that we avoid div by zero (bmy, 12/19/06) 
                  IF ( ABS( CABV - CBEL ) > 0d0 ) THEN
                     CHAT(I,K) = LOG( CABV / CBEL )
     &                         /    ( CABV - CBEL )
     &                         *      CABV * CBEL
                  ELSE
                     CHAT(I,K) = SMALLEST
                  ENDIF

               ELSE             

                  ! If CDFIR <= 1d6, just use arithmetic mean
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
            PLUMEIN = MU(I,KK)

            IF ( PLUMEIN > MBSTH ) THEN
                CONU(I,KK) = CMIX(I,KK) 
            ENDIF

            IF ( MD(I,K) < -MBSTH ) THEN
                COND(I,K) = 0.5d0 * ( CMIX(I,KM1) + CONU(I,KM1) )
            ENDIF
         ENDDO

         !==============================================================
         ! Updraft from bottom to top
         !==============================================================
         DO KK = LLPAR-1,1,-1
            KKP1 = MIN( LLPAR, KK+1 )

            DO I = IL1G, IL2G
               PLUMEIN     = MU(I,KKP1) + EU(I,KK) 
               PLUMEOUT    = MU(I,KK) + DU(I,KK) - 0.5D0*MD(I,KK)
               PLUMECHANGE = PLUMEOUT - PLUMEIN

               IF ( PLUMECHANGE > MBSTH ) THEN
                  IF ( PLUMEOUT > MBSTH ) THEN
                     CONU(I,KK) = (MU(I,KKP1)*CONU(I,KKP1) *FISG(I,KK)
     &                          + EU(I,KK)*CMIX(I,KK)
     &                          + PLUMECHANGE*CMIX(I,KK)  )
     &                          / PLUMEOUT
                  ENDIF   
                  
               ELSE
                  IF ( PLUMEIN > MBSTH ) THEN
                     CONU(I,KK) = ( MU(I,KKP1)*CONU(I,KKP1) *FISG(I,KK)
     &                          + EU(I,KK)*CMIX(I,KK) )
     &                          / PLUMEIN
                  ENDIF
               ENDIF

               IF ( CONU(I,KK) < 0.0D0 ) THEN
                  WRITE(6,*) 'Warning! negative conu!!!', conu(I,KK)
                  CALL FLUSH(6)
               !ELSE IF ( CONU(I,KK) > 1.0e-10 ) THEN
               !    write(6,*) 'Warning! Too big conu!!!', conu(I,KK)
               !    call flush(6)
               ENDIF
            ENDDO
         ENDDO

         !==============================================================
         ! Downdraft from top to bottom
         !==============================================================
         DO K = 3, LLPAR
            KM1 = MAX( 1, K-1 )

            DO I = IL1G, IL2G

               IF ( TOTALMD(I,K) < -MBSTH ) THEN
                  IF ( MD(I,K) < -MBSTH ) THEN
                     COND(I,K) = ( TOTALMD(I,KM1)*COND(I,KM1) 
     $                  + 0.5D0 * MD(I,K) * ( CMIX(I,K)+CONU(I,K) ))
     $                  / TOTALMD(I,K)
                  ELSE
                     COND(I,K) = COND(I,KM1)
                 ENDIF
               ENDIF
               
               IF ( COND(I,K) < 0.0D0 ) THEN
                  WRITE(6,*) 'WARNING! negative cond!!!', cond(I,K)
                  CALL FLUSH(6)
                !ELSE IF ( COND(I,K) > 1.0e-10 ) THEN
                !   write(6,*) 'Warning! Too big cond!!!', cond(I,K)
                !   call flush(6)
               ENDIF
            ENDDO
         ENDDO

         DO K = 1, LLPAR
            KM1 = MAX( 1,     K-1 )
            KP1 = MIN( LLPAR, K+1 )
            
            DO I = IL1G, IL2G

               ! Version 3 limit fluxes outside convection to mass in 
               ! appropriate layer.  These limiters are probably only safe
               ! for positive definite quantitities.  It assumes that mu 
               ! and md already satify a courant number limit of 1

!               FLUXIN =  MU(I,KP1)* CONU(I,KP1) * FISG(I,K)
!     $                + (MU(I,K)+ totalMD(I,K)) * CMIX(I,KM1) 
!     $                -  totalMD(I,K)  * COND(I,K)
!   
!              FLUXOUT =  MU(I,K)   * CONU(I,K)     
!     $                + (MU(I,KP1)+ totalMD(I,KP1))*CMIX(I,K)
!     $                 - totalMD(I,KP1) * COND(I,KP1) 

               IF ( K == LLPAR ) THEN

                  FLUXIN  = MU(I,K)        * CMIX(I,KM1)              
     &                    - TOTALMD(I,KM1) * COND(I,KM1)

                  FLUXOUT = MU(I,K)        * CONU(I,K) 
     &                    - TOTALMD(I,KM1) * CMIX(I,K)

               ELSE
           
                  FLUXIN  =  MU(I,KP1)      * CONU(I,KP1) * FISG(I,K)
     &                    +  MU(I,K)        * CMIX(I,KM1) 
     &                    -  TOTALMD(I,KM1) * COND(I,KM1)
     &                    -  TOTALMD(I,K)   * CMIX(I,KP1) * FISG(I,K)
   
                  FLUXOUT = MU(I,K)        * CONU(I,K)     
     &                    + MU(I,KP1)      * CMIX(I,K)
     &                    - TOTALMD(I,K)   * COND(I,K) 
     &                    - TOTALMD(I,KM1) * CMIX(I,K)
               ENDIF

!!!!!!!!!!!!!!!!!!!backup: also works OK !!!!!!!!!!!!!!!!!!!!!!!
!!              FLUXIN =  MU(I,KP1)* CONU(I,KP1) 
!!     $                +  MU(I,K)  * 0.5d0*(CHAT(I,K)+CMIX(I,KM1)) 
!!     $                -  MD(I,K)  * COND(I,K)   
!!     $                -  MD(I,KP1)* 0.5d0*(CHAT(I,KP1)+CMIX(I,KP1))
!!
!!               FLUXOUT = MU(I,K)   * CONU(I,K)     
!!     $                 + MU(I,KP1) * 0.5d0*(CHAT(I,KP1)+CMIX(I,K))
!!     $                 - MD(I,KP1) * COND(I,KP1) 
!!     $                 - MD(I,K)   * 0.5d0*(CHAT(I,K)+CMIX(I,K))
!!
!!               FLUXIN =  MU(I,KP1)* CONU(I,KP1) 
!!     $                +  MU(I,K)  * CHAT(I,K)
!!     $                -  MD(I,K)  * COND(I,K)   
!!     $                -  MD(I,KP1)* CHAT(I,KP1)
!!
!!               FLUXOUT = MU(I,K)   * CONU(I,K)     
!!     $                 + MU(I,KP1) * CHAT(I,KP1)
!!     $                 - MD(I,KP1) * COND(I,KP1) 
!!     $                 - MD(I,K)   * CHAT(I,K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               !==================================================
               ! ND38 Diagnostic: loss of soluble tracer to wet
               ! scavenging in cloud updrafts [kg/s].  
               !==================================================
               NN = INDEXSOL(M)

               IF ( ND38 > 0 .and. NN > 0 ) THEN
 
                  ! Grid box indices 
                  II = IDEEP(I)
                  JJ = J
                  LL = LLPAR - K + 1

                  ! Grid box surface area [m2] 
                  AREA_M2 = GET_AREA_M2( II, JJ, LL ) 

                  ! Save into AD38 array [kg/s]
                  AD38(II,JJ,LL,NN) = AD38(II,JJ,LL,NN) 
     &                 +  MU(I,KP1)    * AREA_M2 / GRAV * CONU(I,KP1) 
     &                 * (1-FISG(I,K)) / TCVV(M) / FLOAT(NSTEP) 
     &                 -  TOTALMD(I,K) * AREA_M2 / GRAV * CMIX(I,KP1) 
     &                 * (1-FISG(I,K)) / TCVV(M) / FLOAT(NSTEP)
               ENDIF

               IF ( ND14 > 0 ) THEN 
                  II = IDEEP(I)
                  JJ = J
                  LL = LLPAR - K + 1
                  
                  ! Grid box surface area [m2]
                  AREA_M2 = GET_AREA_M2( II, JJ, LL ) 

                  CONVFLUP(II,JJ,LL,M) = CONVFLUP(II,JJ,LL,M)  
     &              + MU(I,K) * AREA_M2  * (CONU(I,K)-CMIX(I,KM1))
     &              / GRAV / TCVV(M) / FLOAT(NSTEP)
     &              - TOTALMD(I,KM1) * AREA_M2 * (CMIX(I,K)-COND(I,KM1)) 
     &              / GRAV / TCVV(M) / FLOAT(NSTEP)
                  
               ENDIF 

               NETFLUX = FLUXIN - FLUXOUT

               IF ( DP(I,K)< 0.0D0 ) THEN 
                  WRITE(6,*) 'WARNING! negative DP!!!', DP(I,K)
                  CALL FLUSH(6)
               ENDIF


               DCONDT(I,K)= NETFLUX/DP(I,K) !AD(IDEEP(I),lati_index,llpar+1-k)
            ENDDO               !I
         ENDDO                  !K


         DO K = KBM, LLPAR             
            KM1 = MAX( 1, K-1 )
            
            DO I = IL1G, IL2G

              !!!temp diag ATTENTION HERE!!!!

               IF ( K == (MX(I) + 100000) ) THEN
                  
                  FLUXIN  =(MU(I,K)+MD(I,K))* CMIX(I,KM1)              
     $                    - MD(I,K)*COND(I,K)

                  FLUXOUT = MU(I,K)*CONU(I,K) 

!!!!!!!!!!!!!!!!!!!!!!BACK UP; also works well !!!!!!!!!!!!!!!!!!!!!
!                  FLUXIN  = MU(I,K)*0.5d0*(CHAT(I,K)+CMIX(I,KM1))
!     $                    - MD(I,K)*COND(I,K)
!
!                  FLUXOUT = MU(I,K)*CONU(I,K) 
!     $                    - MD(I,K)*0.5d0*(CHAT(I,K)+CMIX(I,K))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                  NETFLUX = FLUXIN - FLUXOUT

                  IF (ABS(NETFLUX).LT.MAX(FLUXIN,FLUXOUT)*TINYNUM) THEN
                     NETFLUX = 0.d0
                  ENDIF
                  
                  DCONDT(I,K) = NETFLUX / DP(I,K)
                  
               ELSE IF ( K > MX(I) ) THEN

                  !!!!DCONDT(I,K) = 0.D0

               ENDIF

            ENDDO  !I
         ENDDO     !K

         !==============================================================
         ! Update and scatter data back to full arrays
         !==============================================================
         DO K = 1, LLPAR
            KP1 = MIN( LLPAR, K+1 )
            DO I = IL1G, IL2G    
            
               QN = CMIX(I,K) + DCONDT(I,K) * DELT 

               ! Do not make Q negative!!!
               IF ( QN < 0d0 ) then
                  QN = 0D0
               ENDIF            

               Q(IDEEP(I),J,K,M) = QN
            ENDDO   
         ENDDO      
         
      ENDDO   ! End of tracer loop

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

      ! End of module
      END MODULE GCAP_CONVECT_MOD
