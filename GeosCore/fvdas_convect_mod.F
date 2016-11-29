! $Id: fvdas_convect_mod.f,v 1.1 2009/09/16 14:06:28 bmy Exp $
      MODULE FVDAS_CONVECT_MOD
!
!******************************************************************************
!  Module FVDAS_CONVECT_MOD contains routines (originally from NCAR) which 
!  perform shallow and deep convection for the GEOS-4/fvDAS met fields.  
!  These routines account for shallow and deep convection, plus updrafts 
!  and downdrafts.  (pjr, dsa, bmy, 6/26/03, 12/19/06)
!  
!  Module Variables:
!  ============================================================================
!  (1 ) RLXCLM   (LOGICAL) : Logical to relax column versus cloud triplet
!  (2 ) LIMCNV   (INTEGER) : Maximum CTM level for HACK convection
!  (3 ) CMFTAU   (REAL(fp) ) : Characteristic adjustment time scale for HACK [s]
!  (4 ) EPS      (REAL(fp) ) : A very small number [unitless]
!  (5 ) SMALLEST (REAL(fp) ) : The smallest double-precision number 
!  (6 ) TINYNUM  (REAL(fp) ) : 2 times EPS
!  (7 ) TINYALT  (REAL(fp) ) : arbitrary small num used in transport estimates
!
!  Module Routines:
!  ============================================================================
!  (1 ) INIT_FVDAS_CONVECT : Initializes fvDAS convection scheme
!  (2 ) FVDAS_CONVECT      : fvDAS convection routine, called from MAIN 
!  (3 ) HACK_CONV          : HACK convection scheme routine
!  (4 ) ARCCONVTRAN        : Sets up fields for ZHANG/MCFARLANE convection
!  (5 ) CONVTRAN           : ZHANG/MCFARLANE convection scheme routine
!  (6 ) WHENFGT            : Returns index array of points > a reference value
!
!  NOTES: 
!  (1 ) Contains new updates for GEOS-4/fvDAS convection.  Also added OpenMP
!        parallel loop over latitudes in FVDAS_CONVECT. (swu, bmy, 1/21/04)
!  (2 ) Now prevent FTMP, QTMP arrays from being held PRIVATE w/in the
!        parallel loop in routine DO_CONVECTION (bmy, 7/20/04)
!  (3 ) Now pass wet-scavenged Hg2 to "ocean_mercury_mod.f" (sas, bmy, 1/21/05)
!  (4 ) Rewrote parallel loops to avoid problems w/ OpenMP.  Also modified
!        for updated Hg simulation. (cdh, bmy, 2/1/06)
!  (5 ) Rewrote DO loops in HACK_CONV for better optmization (bmy, 3/28/06)
!  (6 ) Split up Hg2 IF statement into 2 separate statements (bmy, 4/17/06)
!  (7 ) Minor fix in ND38 diagnostic: replace 1 w/ 1d0 (bmy, 5/24/06)
!  (8 ) Updated for ND14 diagnostic.  Now treat "negative" detrainment as 
!        entrainment, which will better conserve mixing ratio in convection.
!        (swu, bmy, 6/27/06)
!  (9 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (10) Bug fix in CONVTRAN to avoid div potential div by zero.  Make SMALLEST
!        = 1d-60 to avoid problems (bmy, 12/19/06)
!  01 Mar 2012 - R. Yantosca - Now references the new grid_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  28 Oct 2013 - R. Yantosca - Now avoid array temporaries
!  13 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!                              Remove dependency on pressure_mod.
!  27 Apr 2015 - M. Yannetti - Fixed SMALLEST
!  23 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v)
!  06 Jan 2016 - E. Lundgren - Use global parameter g0 instead of local GRAV
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!******************************************************************************
!
      USE PhysConstants, ONLY : g0
      USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

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
      REAL(fp),  PARAMETER :: CMFTAU   = 3600.e+0_fp
      REAL(fp),  PARAMETER :: EPS      = 1.0e-13_fp   
      !REAL(fp),  PARAMETER :: GRAV     = 9.8e+0_fp ! Use global g0
      REAL(fp),  PARAMETER :: SMALLEST = TINY(1.0e+0_fp)
      REAL(fp),  PARAMETER :: TINYALT  = 1.0e-36_fp       
      REAL(fp),  PARAMETER :: TINYNUM  = 2*SMALLEST

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE INIT_FVDAS_CONVECT( am_I_Root, State_Met )
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
      USE CMN_SIZE_MOD             ! Size parameters
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)         :: am_I_Root   ! Are we on the root CPU?
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L, L2
      REAL(fp)            :: HYPI(LLPAR+1)

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
         HYPI(L2) = State_Met%PEDGE(I,J,L)
      ENDDO

      ! Limit convection to regions below 40 hPa
      IF ( HYPI(1) >= 40e+0_fp ) THEN
         LIMCNV = 1
      ELSE
         DO L = 1, LLPAR
            IF ( HYPI(L) < 40e+0_fp .AND. HYPI(L+1) >= 40e+0_fp ) THEN
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
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 )  LIMCNV, HYPI(LIMCNV) 
 100     FORMAT( '     - GEOS-4 convection is capped at L = ', i3, 
     &          ', or approx ', f6.1, ' hPa' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_FVDAS_CONVECT

!------------------------------------------------------------------------------

      SUBROUTINE FVDAS_CONVECT( TDT,   NTRACE, Q,     RPDEL, ETA, 
     &                          BETA,  MU,     MD,    EU,    DP,    
     &                          NSTEP, FRACIS, INDEXSOL,
     &                          Input_Opt, State_Met, State_Chm )
!
!******************************************************************************
!  Subroutine FVDAS_CONVECT is the convection driver routine for GEOS-4/fvDAS
!  met fields.  It calls both HACK and ZHANG/MCFARLANE convection schemes.
!  (pjr, dsa, bmy, 6/26/03, 12/13/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TDT     (REAL(fp) ) : 2 * delta-T                          [s]
!  (2 ) NTRACE  (INTEGER  ) : Number of tracers to transport       [unitless]
!  (3 ) Q       (REAL(fp) ) : Array of transported tracers  [kg/kg total air]
!  (4 ) RPDEL   (REAL(fp) ) : 1 / DP                               [1/hPa]
!  (5 ) ETA     (REAL(fp) ) : GMAO Hack convective mass flux       [kg/m2/s]
!  (6 ) BETA    (REAL(fp) ) : GMAO Hack overshoot parameter        [unitless]
!  (7 ) MU      (REAL(fp) ) : GMAO updraft mass flux   (ZMMU)      [Pa/s]
!  (8 ) MD      (REAL(fp) ) : GMAO downdraft mass flux (ZMMD)      [Pa/s]
!  (9 ) EU      (REAL(fp) ) : GMAO updraft entrainment (ZMEU)      [Pa/s]
!  (10) DP      (REAL(fp) ) : Delta-pressure between level edges   [Pa]
!  (11) NSTEP   (INTEGER  ) : Time step index                      [unitless]
!  (12) FRACIS  (REAL(fp) ) : Fraction of tracer that is insoluble [unitless]
!  (13) INDEXSOL(INTEGER  ) : Index array of soluble tracers       [unitless]
!  (14) Input_opt(OPTINPUT) : Input options object
!  (15) State_Met(METSTATE) : Meteorological state object
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) Q       (REAL(fp) ) : Modified tracer array          [kg/kg total air]
! 
!  Important Local Variables:
!  ============================================================================
!  (1 ) LENGATH (INTEGER) : Number of lons where deep conv. happens at lat=J
!  (2 ) IDEEP   (INTEGER) : Lon indices where deep convection happens at lat=J
!  (3 ) JT      (INTEGER) : Cloud top layer for columns undergoing conv.
!  (4 ) MX      (INTEGER) : Cloud bottom layer for columns undergoing conv.
!  (5 ) DSUBCLD (REAL(fp) ) : Delta pressure from cloud base to sfc
!  (6 ) DU      (REAL(fp) ) : Mass detraining from updraft (lon-alt array)
!  (7 ) ED      (REAL(fp) ) : Mass entraining from downdraft (lon-alt array)
!  (8 ) DPG     (REAL(fp) ) : gathered .01*dp (lon-alt array)
!  (9 ) EUG     (REAL(fp) ) : gathered eu (lon-alt array) 
!  (10) MUG     (REAL(fp) ) : gathered mu (lon-alt array)   
!  (11) MDG     (REAL(fp) ) : gathered md (lon-alt array)
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
!  28 Oct 2013 - R. Yantosca - Now accept pointer-valued arguments into this
!                              routine via assumed-shape dummy arguments.
!                              This will eliminate array temporaries.
!  23 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v)
!******************************************************************************
!
! !USES
!
      USE CMN_SIZE_MOD
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Met_Mod,      ONLY : MetState
      USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)    :: NSTEP, NTRACE             
      INTEGER,        INTENT(IN)    :: INDEXSOL(NTRACE) 
      REAL(fp),       INTENT(IN)    :: TDT                
      REAL(fp),       INTENT(IN)    :: RPDEL(IIPAR,JJPAR,LLPAR)  
      REAL(fp),       INTENT(IN)    :: ETA (:,:,:)
      REAL(fp),       INTENT(IN)    :: BETA(:,:,:)
      REAL(fp),       INTENT(IN)    :: MU  (:,:,:)
      REAL(fp),       INTENT(IN)    :: MD  (:,:,:)
      REAL(fp),       INTENT(IN)    :: EU  (:,:,:)
      REAL(fp),       INTENT(IN)    :: DP(IIPAR,JJPAR,LLPAR)     
      REAL(fp),       INTENT(IN)    :: FRACIS(IIPAR,JJPAR,LLPAR,NTRACE) 
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NTRACE)
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J, L, N, LENGATH, ISTEP
      INTEGER                       :: JT(IIPAR)
      INTEGER                       :: MX(IIPAR)
      INTEGER                       :: IDEEP(IIPAR)
      REAL(fp)                      :: DSUBCLD(IIPAR)
      REAL(fp)                      :: DPG(IIPAR,LLPAR)
      REAL(fp)                      :: DUG(IIPAR,LLPAR)
      REAL(fp)                      :: EDG(IIPAR,LLPAR)
      REAL(fp)                      :: EUG(IIPAR,LLPAR)
      REAL(fp)                      :: MDG(IIPAR,LLPAR)
      REAL(fp)                      :: MUG(IIPAR,LLPAR)

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
     &                    DUG,   EUG,       EDG,    DPG,  DSUBCLD,  
     &                    JT,    MX,        IDEEP,  1,    LENGATH, 
     &                    NSTEP, 0.5e+0_fp*TDT, FRACIS, 
     &                    INDEXSOL, Input_Opt, State_Met, State_Chm )
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
!  Hack scheme. (pjr, dsa, bmy, 6/26/03, 3/28/06)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) J      (INTEGER)   : GEOS-CHEM Latitude index               [unitless]
!  (2 ) TDT    (REAL(fp) ) : 2 delta-t                              [s       ]
!  (3 ) RPDEL  (REAL(fp) ) : Reciprocal of pressure-thickness array [1/hPa   ]
!  (4 ) ETA    (REAL(fp) ) : GMAO Hack convective mass flux (HKETA) [kg/m2/s ]
!  (5 ) BETA   (REAL(fp) ) : GMAO Hack overshoot parameter (HKBETA) [unitless]
!  (6 ) NTRACE (INTEGER)   : Number of tracers in the Q array       [unitless]
!  (7 ) Q      (REAL(fp) ) : Tracer concentrations           [kg/kg total air]
!  
!  Arguments as Output:
!  ============================================================================
!  (7 ) Q      (REAL(fp) ) : Modified tracer concentrations  [kg/kg total air]
!
!  Important Local Variables:
!  ============================================================================
!  (1 ) ADJFAC (REAL(fp) ) : Adjustment factor (relaxation related)
!  (2 ) BOTFLX (REAL(fp) ) : Bottom constituent mixing ratio flux
!  (3 ) CMRC   (REAL(fp) ) : Constituent mixing ratio ("in-cloud")
!  (4 ) CMRH   (REAL(fp) ) : Interface constituent mixing ratio 
!  (5 ) DCMR1  (REAL(fp) ) : Q convective change (lower level)
!  (6 ) DCMR2  (REAL(fp) ) : Q convective change (mid level)
!  (7 ) DCMR3  (REAL(fp) ) : Q convective change (upper level)
!  (8 ) EFAC1  (REAL(fp) ) : Ratio q to convectively induced change (bot level)
!  (9 ) EFAC2  (REAL(fp) ) : Ratio q to convectively induced change (mid level)
!  (10) EFAC3  (REAL(fp) ) : Ratio q to convectively induced change (top level)
!  (11) ETAGDT (REAL(fp) ) : ETA * G0 * DT
!  (12) TOPFLX (REAL(fp) ) : Top constituent mixing ratio flux
!
!  NOTES:
!  (1 ) Updated comments.  Added NTRACE as an argument.  Now also force 
!        double-precision with the "D" exponents.  (bmy, 6/26/03)
!  (2 ) Now pass J via the arg list.  Now dimension RPDEL, ETA, BETA, and Q
!        with and make all input arrays dimensioned
!        with (IIPAR,JJPAR,LLPAR,...) to avoid seg fault error in OpenMP
!        on certain platforms.
!  (3 ) Rewrote DO loops and changed 1-D arrays into scalars in order to
!        improve optimization, particularly for the Intel IFORT v9 compiler.
!        (bmy, 3/28/06)
!******************************************************************************
!
      USE CMN_SIZE_MOD      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: J, NTRACE
      REAL(fp),  INTENT(IN)    :: TDT
      REAL(fp),  INTENT(IN)    :: RPDEL(IIPAR,JJPAR,LLPAR)
      REAL(fp),  INTENT(IN)    :: ETA (:,:,:)
      REAL(fp),  INTENT(IN)    :: BETA(:,:,:)
      REAL(fp),  INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NTRACE)

      ! Local variables
      INTEGER                :: I,      K,      M
      REAL(fp)               :: ADJFAC, BOTFLX, TOPFLX              
      REAL(fp)               :: EFAC1,  EFAC2,  EFAC3
      REAL(fp)               :: CMRC,   DCMR1,  DCMR2        
      REAL(fp)               :: DCMR3,  ETAGDT, CMRH(IIPAR,LLPAR+1)   

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
         ADJFAC = 1e+0_fp
      ENDIF

      !=================================================================
      ! Begin moist convective mass flux adjustment procedure. 
      ! The formalism ensures that negative cloud liquid water can 
      ! never occur.
      !
      ! Rewrote DO loops and changed 1-D arrays into scalars in order
      ! to optimization, esp. for Intel IFORT compiler. (bmy, 3/28/06)
      !=================================================================
      
      ! Loop over tracers
      DO M = 1, NTRACE

         ! Initialize
         CMRH(:,:) = 0e+0_fp

      ! Loop over levels and longitudes
      DO K = LLPAR-1, LIMCNV+1, -1
      DO I = 1,       IIPAR

         ! Initialize
         ETAGDT    = 0e+0_fp
         CMRC      = 0e+0_fp
         BOTFLX    = 0e+0_fp
         TOPFLX    = 0e+0_fp
         EFAC1     = 0e+0_fp
         EFAC2     = 0e+0_fp
         EFAC3     = 0e+0_fp
         DCMR1     = 0e+0_fp
         DCMR2     = 0e+0_fp
         DCMR3     = 0e+0_fp
         
         ! Only proceed for boxes with nonzero mass flux
         IF ( ETA(I,J,K) > 0e+0_fp ) THEN
            ETAGDT = ETA(I,J,K) * G0 * TDT * 0.01e+0_fp   ![hPa]
         ELSE
            CYCLE
         ENDIF
         
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

         ! If any of the reported values of the constituent is 
         ! negative in the three adjacent levels, nothing will 
         ! be done to the profile.  Skip to next longitude.
         IF ( ( Q(I,J,K+1,M) < 0e+0_fp )  .OR. 
     &        ( Q(I,J,K,  M) < 0e+0_fp )  .OR.
     &        ( Q(I,J,K-1,M) < 0e+0_fp ) ) CYCLE
         
         ! Specify constituent interface values (linear interpolation)
         CMRH(I,K  ) = 0.5e+0_fp *( Q(I,J,K-1,M) + Q(I,J,K  ,M) )
         CMRH(I,K+1) = 0.5e+0_fp *( Q(I,J,K  ,M) + Q(I,J,K+1,M) )
              
         ! In-cloud mixing ratio
         CMRC        = Q(I,J,K+1,M)

         ! Determine fluxes, flux divergence => changes due to convection.  
         ! Logic must be included to avoid producing negative values. 
         ! A bit messy since there are no a priori assumptions about profiles.
         ! Tendency is modified (reduced) when pending disaster detected.
         BOTFLX = ETAGDT * ( CMRC - CMRH(I,K+1) ) * ADJFAC
         TOPFLX = BETA(I,J,K) * ETAGDT * ( CMRC - CMRH(I,K) ) * ADJFAC
         DCMR1  = -BOTFLX * RPDEL(I,J,K+1)
         EFAC1  = 1.0e+0_fp
         EFAC2  = 1.0e+0_fp
         EFAC3  = 1.0e+0_fp
               
         ! K+1th level
         IF ( Q(I,J,K+1,M) + DCMR1 < 0e+0_fp ) THEN
            EFAC1 = MAX( TINYALT, ABS( Q(I,J,K+1,M) / DCMR1 ) - EPS )
         ENDIF

         IF ( EFAC1 == TINYALT .or. EFAC1 > 1e+0_fp ) EFAC1 = 0e+0_fp
         DCMR1 =  -EFAC1 * BOTFLX            * RPDEL(I,J,K+1)
         DCMR2 = ( EFAC1 * BOTFLX - TOPFLX ) * RPDEL(I,J,K)
               
         ! Kth level
         IF ( Q(I,J,K,M) + DCMR2 < 0e+0_fp ) THEN
            EFAC2 = MAX( TINYALT, ABS( Q(I,J,K,M) / DCMR2 ) - EPS )
         ENDIF
               
         IF ( EFAC2 == TINYALT .or. EFAC2 > 1e+0_fp ) EFAC2 = 0e+0_fp
         DCMR2 = ( EFAC1 * BOTFLX - EFAC2 * TOPFLX ) * RPDEL(I,J,K)
         DCMR3 =                    EFAC2 * TOPFLX   * RPDEL(I,J,K-1)
         
         ! K-1th level
         IF ( Q(I,J,K-1,M) + DCMR3 < 0e+0_fp ) THEN
            EFAC3 = MAX( TINYALT, ABS( Q(I,J,K-1,M) / DCMR3 ) - EPS )
         ENDIF

         IF ( EFAC3 == TINYALT .or. EFAC3 > 1e+0_fp ) EFAC3 = 0e+0_fp
         EFAC3 = MIN( EFAC2, EFAC3 )
         DCMR2 = ( EFAC1 * BOTFLX - EFAC3 * TOPFLX ) * RPDEL(I,J,K)
         DCMR3 =                    EFAC3 * TOPFLX   * RPDEL(I,J,K-1)
               
         ! Save back into tracer array (levels K+1, K, K-1)
         Q(I,J,K+1,M) = Q(I,J,K+1,M) + DCMR1
         Q(I,J,K  ,M) = Q(I,J,K  ,M) + DCMR2
         Q(I,J,K-1,M) = Q(I,J,K-1,M) + DCMR3
      ENDDO
      ENDDO
      ENDDO
      
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
!  (pjr, dsa, swu, bmy, 6/26/03, 6/27/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J       (INTEGER)   : GEOS-CHEM latitude index          [unitless]
!  (2 ) DP      (REAL(fp) ) : Delta pressure between interfaces [Pa      ]      
!  (3 ) MU      (REAL(fp) ) : Mass flux up                      [kg/m2/s ] 
!  (4 ) MD      (REAL(fp) ) : Mass flux down                    [kg/m2/s ] 
!  (5 ) EU      (REAL(fp) ) : Mass entraining from updraft      [1/s     ]     
!
!  Arguments as Output:
!  ============================================================================
!  (6 ) MUG     (REAL(fp) ) : Gathered mu (lon-alt array)
!  (7 ) MDG     (REAL(fp) ) : Gathered md (lon-alt array)
!  (8 ) DUG     (REAL(fp) ) : Mass detraining from updraft (lon-alt array)
!  (9 ) EUG     (REAL(fp) ) : Gathered eu (lon-alt array)
!  (10) EDG     (REAL(fp) ) : Mass entraining from downdraft (lon-alt array)
!  (11) DPG     (REAL(fp) ) : Gathered .01*dp (lon-alt array)
!  (12) DSUBCLD (REAL(fp) ) : Delta pressure from cloud base to sfc
!                              (lon-alt arr)
!  (13) JTG     (INTEGER)   : Cloud top layer for columns undergoing conv.
!  (14) JBG     (INTEGER)   : Cloud bottom layer for columns undergoing conv.
!  (15) IDEEP   (INTEGER)   : Index of longitudes where deep conv. happens
!  (16) LENGATH (INTEGER)   : Length of gathered arrays
! 
!  NOTES:
!  (1 ) Removed NSTEP from arg list; it's not used.  Also zero arrays in order
!        to prevent them from being filled with compiler junk for latitudes
!        where no convection occurs at all. (bmy, 1/21/04)
!  (2 ) Now dimension DP, MU, MD, EU as (IIPAR,JJPAR,LLPAR) to avoid seg fault
!        error in OpenMP.  Also now pass the GEOS-CHEM latitude index J via
!        the argument list. (bmy, 12/13/05)
!  (3 ) Now treat "negative detrainment" as entrainment, which will better 
!        conserve mixing ratio (swu, bmy, 6/27/06)
!  28 Oct 2013 - R. Yantosca - Eliminate array temporaries
!******************************************************************************
!
      USE CMN_SIZE_MOD    ! Size parameters
      
      ! Arguments
      INTEGER,   INTENT(IN)  :: J 
      REAL(fp),  INTENT(IN)  :: DP(IIPAR,JJPAR,LLPAR) 
      REAL(fp),  INTENT(IN)  :: MU (:,:,:)
      REAL(fp),  INTENT(IN)  :: MD (:,:,:)
      REAL(fp),  INTENT(IN)  :: EU (:,:,:)
      REAL(fp),  INTENT(OUT) :: MUG(IIPAR,LLPAR)
      REAL(fp),  INTENT(OUT) :: MDG(IIPAR,LLPAR)
      REAL(fp),  INTENT(OUT) :: DUG(IIPAR,LLPAR)      
      REAL(fp),  INTENT(OUT) :: EUG(IIPAR,LLPAR)
      REAL(fp),  INTENT(OUT) :: EDG(IIPAR,LLPAR)
      REAL(fp),  INTENT(OUT) :: DPG(IIPAR,LLPAR)
      REAL(fp),  INTENT(OUT) :: DSUBCLD(IIPAR)   
      INTEGER,   INTENT(OUT) :: JTG(IIPAR)
      INTEGER,   INTENT(OUT) :: JBG(IIPAR)
      INTEGER,   INTENT(OUT) :: IDEEP(IIPAR)
      INTEGER,   INTENT(OUT) :: LENGATH

      ! Local variables
      INTEGER                :: I, K, LENPOS 
      INTEGER                :: INDEX(IIPAR)
      REAL(fp)               :: SUM(IIPAR)
      REAL(fp)               :: RDPG(IIPAR,LLPAR)      

      !=================================================================
      ! ARCONVTRAN begins here!
      !=================================================================

      ! Initialize arrays
      DPG     = 0e+0_fp
      DSUBCLD = 0e+0_fp
      DUG     = 0e+0_fp
      EDG     = 0e+0_fp
      EUG     = 0e+0_fp
      JTG     = LLPAR
      JBG     = 1
      MDG     = 0e+0_fp
      MUG     = 0e+0_fp
      RDPG    = 0e+0_fp
      SUM     = 0e+0_fp
      
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
      CALL WHENFGT( IIPAR, SUM, 1, 0e+0_fp, IDEEP, LENGATH )
      
      ! Return if there is no convection the longitude band
      IF ( LENGATH == 0 ) RETURN

      !=================================================================
      ! Gather input mass fluxes in places where there is convection
      !=================================================================
      DO K = 1, LLPAR
      DO I = 1, LENGATH

         ! Convert Pa->hPa
         DPG(I,K)  = 0.01e+0_fp * DP(IDEEP(I),J,K)  
         RDPG(I,K) = 1.e+0_fp / DPG(I,K)

         ! Convert Pa/s->hPa/s
         MUG(I,K)  = MU(IDEEP(I),J,K) *  0.01e+0_fp           
         MDG(I,K)  = MD(IDEEP(I),J,K) *  0.01e+0_fp

         ! Convert Pa/s->1/s
         EUG(I,K)  = EU(IDEEP(I),J,K) *  0.01e+0_fp * RDPG(I,K) 
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
         EDG(I,LLPAR) = 0.0e+0_fp
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
         CALL WHENFGT( LENGATH, MUG(:,K), 1, 0e+0_fp, INDEX, LENPOS )
         
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
     &                     NSTEP, DELT,   FRACIS, INDEXSOL,
     &                     Input_Opt, State_Met, State_Chm )
!
!******************************************************************************
!  Subroutine CONVTRAN applies the convective transport of trace species
!  (assuming moist mixing ratio) using the ZHANG/MCFARLANE convection scheme. 
!  (pjr, dsa, bmy, 6/26/03, 12/19/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J        (INTEGER)   : GEOS-CHEM latitude index              [unitless]
!  (2 ) NTRACE   (INTEGER)   : Number of tracers to transport        [unitless]
!  (3 ) Q        (REAL(fp) ) : Tracer conc. including moisture [kg/kg total air]
!  (4 ) MU       (REAL(fp) ) : Mass flux up                          [hPa/s   ]
!  (5 ) MD       (REAL(fp) ) : Mass flux down                        [hPa/s   ]
!  (6 ) DU       (REAL(fp) ) : Mass detraining from updraft          [1/s     ]
!  (7 ) EU       (REAL(fp) ) : Mass entraining from updraft          [1/s     ]
!  (8 ) ED       (REAL(fp) ) : Mass entraining from downdraft        [1/s     ]
!  (9 ) DP       (REAL(fp) ) : Delta pressure between interfaces
!  (10) DSUBCLD  (REAL(fp) ) : Delta pressure from cloud base to sfc
!  (11) JT       (INTEGER)   : Index of cloud top for each column
!  (12) MX       (INTEGER)   : Index of cloud top for each column
!  (13) IDEEP    (INTEGER)   : Gathering array
!  (14) IL1G     (INTEGER)   : Gathered min lon indices over which to operate
!  (15) IL2G     (INTEGER)   : Gathered max lon indices over which to operate
!  (16) NSTEP    (INTEGER)   : Time step index
!  (17) DELT     (REAL(fp) ) : Time step
!  (18) FRACIS   (REAL(fp) ) : Fraction of tracer that is insoluble
!  (19) INDEXSOL (INTEGER)   : Index array of soluble tracer numbers
!  (20) Input_opt(OPTINPUT)  : Input options object
!  (21) State_Met(METSTATE)  : Meteorological state object
!  (21) State_Chm(CHMSTATE)  : Chemistry state object
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) Q       (REAL(fp) )  : Contains modified tracer mixing ratios 
!                              [kg/kg total air] 
!  Important Local Variables:
!  ============================================================================
!  (1 ) CABV    (REAL(fp) )  : Mixing ratio of constituent above
!  (2 ) CBEL    (REAL(fp) )  : Mix ratio of constituent beloqw
!  (3 ) CDIFR   (REAL(fp) )  : Normalized diff between cabv and cbel
!  (4 ) CHAT    (REAL(fp) )  : Mix ratio in env at interfaces
!  (5 ) CMIX    (REAL(fp) )  : Gathered tracer array 
!  (6 ) COND    (REAL(fp) )  : Mix ratio in downdraft at interfaces
!  (7 ) CONU    (REAL(fp) )  : Mix ratio in updraft at interfaces
!  (8 ) DCONDT  (REAL(fp) )  : Gathered tend array 
!  (9 ) FISG    (REAL(fp) )  : gathered insoluble fraction of tracer
!  (10) KBM     (INTEGER)    : Highest altitude index of cloud base [unitless]
!  (11) KTM     (INTEGER)    : Hightet altitude index of cloud top  [unitless]
!  (12) MBSTH   (REAL(fp) )  : Threshold for mass fluxes
!  (13) SMALL   (REAL(fp) )  : A small number
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
!  (4 ) Split Hg2 IF statement into 2 IF statements so as to avoid seg faults.
!        (bmy, 4/17/06)
!  (5 ) Replace 1 with 1d0 in ND38 diagnostic (bmy, 5/24/06)
!  (6 ) Updated for ND14 diagnostic (swu, bmy, 6/12/06)
!  (7 ) Now treat "negative detrainment" as entrainment, which will better 
!        conserve mixing ratio (swu, bmy, 6/27/06)
!  (8 ) Bug fix: avoid div by zero in formula for CHAT (bmy, 12/19/06)
!  (9 ) Updates for mercury simulation (ccc, 5/6/10)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  23 Jun 2014 - R. Yantosca - Now accept Input_Opt
!  23 Jun 2014 - R. Yantosca - Now remove references to logical_mod.F
!  25 Jun 2014 - R. Yantosca - Remove references to tracer_mod.F
!  27 Apr 2015 - E. Lundgren - Replace 100/G0 with global param g0_100
!  15 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v)
!  22 Apr 2016 - R. Yantosca - Now get IS_Hg2, Is_HgP from species database
!  25 Apr 2016 - R. Yantosca - Now pass Hg category # to ADD_Hg2_* functions
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!******************************************************************************
!
      ! References to F90 modules
      USE CMN_DIAG_MOD             ! ND38, LD38, ND14, LD14
      USE CMN_SIZE_MOD             ! Size parameters
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_WD, ADD_HgP_WD !eck/eds
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_SNOWPACK !CDH
#if defined( BPCH_DIAG )
      USE DIAG_MOD,           ONLY : AD38
#endif
      USE DIAG_MOD,           ONLY : CONVFLUP 
      USE GC_GRID_MOD,        ONLY : GET_AREA_M2
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants            ! Universal constants
      USE Species_Mod,        ONLY : Species
      USE State_Met_Mod,      ONLY : MetState
      USE State_Chm_Mod,      ONLY : ChmState

      ! Arguments
      INTEGER,        INTENT(IN)    :: J
      INTEGER,        INTENT(IN)    :: NTRACE  
      REAL(fp),       INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NTRACE)  
      REAL(fp),       INTENT(IN)    :: MU(IIPAR,LLPAR)      
      REAL(fp),       INTENT(IN)    :: MD(IIPAR,LLPAR)      
      REAL(fp),       INTENT(IN)    :: DU(IIPAR,LLPAR)      
      REAL(fp),       INTENT(IN)    :: EU(IIPAR,LLPAR)      
      REAL(fp),       INTENT(IN)    :: ED(IIPAR,LLPAR)      
      REAL(fp),       INTENT(IN)    :: DP(IIPAR,LLPAR)      
      REAL(fp),       INTENT(IN)    :: DSUBCLD(IIPAR)      
      INTEGER,        INTENT(IN)    :: JT(IIPAR)          
      INTEGER,        INTENT(IN)    :: MX(IIPAR)          
      INTEGER,        INTENT(IN)    :: IDEEP(IIPAR)       
      INTEGER,        INTENT(IN)    :: IL1G               
      INTEGER,        INTENT(IN)    :: IL2G 
      INTEGER,        INTENT(IN)    :: NSTEP               
      REAL(fp),       INTENT(IN)    :: DELT                
      REAL(fp),       INTENT(IN)    :: FRACIS(IIPAR,JJPAR,LLPAR,NTRACE) 
      INTEGER,        INTENT(IN)    :: INDEXSOL(NTRACE)
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object

      ! Local variables
      LOGICAL                      :: IS_Hg
      INTEGER                      :: II,      JJ,      LL,      NN
      INTEGER                      :: I,       K,       KBM,     KK     
      INTEGER                      :: KKP1,    KM1,     KP1,     KTM     
      INTEGER                      :: M,       ISTEP,   Hg_Cat
      REAL(fp)                     :: CABV,    CBEL,    CDIFR,   CD2
      REAL(fp)                     :: DENOM,   SMALL,   MBSTH,   MUPDUDP
      REAL(fp)                     :: MINC,    MAXC,    QN,      FLUXIN
      REAL(fp)                     :: D_NSTEP, FLUXOUT, NETFLUX, AREA_M2
      REAL(fp)                     :: WET_Hg2, WET_HgP, PLUMEIN !eck/eds
      REAL(fp)                     :: CHAT(IIPAR,LLPAR)     
      REAL(fp)                     :: COND(IIPAR,LLPAR)     
      REAL(fp)                     :: CMIX(IIPAR,LLPAR)     
      REAL(fp)                     :: FISG(IIPAR,LLPAR)     
      REAL(fp)                     :: CONU(IIPAR,LLPAR)     
      REAL(fp)                     :: DCONDT(IIPAR,LLPAR)   

      ! Pointers
      TYPE(Species), POINTER       :: SpcInfo

      !=================================================================
      ! CONVTRAN begins here!
      !=================================================================

      ! Is this a mercury simulation with dynamic ocean model?
      IS_Hg   = ( Input_Opt%ITS_A_MERCURY_SIM .and. 
     &            Input_Opt%LDYNOCEAN               )

      ! Initialize 
      SpcInfo => NULL()

      ! A small number
      SMALL   = 1.e-36_fp

      ! Threshold below which we treat the mass fluxes as zero (in mb/s)
      MBSTH   = 1.e-15_fp

      ! Convert NSTEP to REAL(fp) for use below
      D_NSTEP = NSTEP

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

         ! Point to the Species Database entry for tracer #M
         SpcInfo => State_Chm%SpcData(M)%Info

         ! Gather up the tracer and set tend to zero
         DO K = 1,    LLPAR
         DO I = IL1G, IL2G
            CMIX(I,K) = Q(IDEEP(I),J,K,M)
            IF ( CMIX(I,K) < 4.e+0_fp*SMALLEST ) CMIX(I,K) = 0e+0_fp
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

               IF ( MINC < 0e+0_fp ) THEN 
                  CDIFR = 0.e+0_fp
               ELSE
                  CDIFR = ABS( CMIX(I,K)-CMIX(I,KM1) ) / MAX(MAXC,SMALL)
               ENDIF
               
               !------------------------------------------------------------
               ! The following 2 variables are actually NOT used
               ! (swu, 12/17/03)
               !DENOM = MAX( MAXC, SMALL )       
               !CD2   = ABS( CMIX(I,K) - CMIX(I,KM1) ) / DENOM
               !------------------------------------------------------------

               IF ( CDIFR > 1.e-6_fp ) THEN

                  ! If the two layers differ significantly.
                  ! use a geometric averaging procedure
                  CABV = MAX( CMIX(I,KM1), MAXC*TINYNUM, SMALLEST )
                  CBEL = MAX( CMIX(I,K),   MAXC*TINYNUM, SMALLEST )

                  ! If CABV-CBEL is zero then set CHAT=SMALLEST
                  ! so that we avoid div by zero (bmy, 12/19/06)
                  IF ( ABS( CABV - CBEL ) > 0e+0_fp ) THEN
                     CHAT(I,K) = LOG( CABV / CBEL )
     &                         /    ( CABV - CBEL )
     &                         *      CABV * CBEL
                  ELSE
                     CHAT(I,K) = SMALLEST
                  ENDIF

               ELSE             

                  ! Small diff, so just arithmetic mean
                  CHAT(I,K) = 0.5e+0_fp * ( CMIX(I,K) + CMIX(I,KM1) )
               ENDIF

               ! Provisional up and down draft values
               CONU(I,K) = CHAT(I,K)
               COND(I,K) = CHAT(I,K)

               ! Provisional tends
               DCONDT(I,K) = 0.e+0_fp
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

            ! Layer LLPAR (ground layer) CLOUD does not have updraft
            ! entering, so assume tracer mixing ratio is same as
            ! the environment (swu, bmy, 6/27/06)
            IF ( MUPDUDP > MBSTH ) THEN
               CONU(I,KK) = CMIX(I,KK)
            ENDIF
            
            ! MD(I,2) is the downdraft entering from layer 2 from 
            ! layer 1 (model top layer); assumed to have the same
            ! mixing ratio as the environment (swu, bmy, 6/27/06)
            IF ( MD(I,K) < -MBSTH ) THEN
               COND(I,K) = CMIX(I,KM1)
            ENDIF 
         ENDDO

         !==============================================================
         ! Updraft from bottom to top
         !==============================================================
         DO KK = LLPAR-1,1,-1
            KKP1 = MIN( LLPAR, KK+1 )

            DO I = IL1G, IL2G

               ! Test for "negative detrainment"
               IF ( DU(I,KK) < 0e+0_fp ) THEN

                  !-----------------------------------------------------
                  ! If negative DU (detrainment) happens, which implies 
                  ! that the input metfields are not well constrained 
                  ! and EU is inaccurate, we apply the correction by 
                  ! treating the negative detrainment as extra 
                  ! entrainment. (swu, bmy, 06/27/06)
                  !-----------------------------------------------------
                  
                  ! Air mass flux going into layer KK of the CLOUD
                  PLUMEIN = MU(I,KKP1) + ( EU(I,KK) * DP(I,KK) )
     &                                 - ( DU(I,KK) * DP(I,KK) )

                  ! Compute concentration
                  IF ( PLUMEIN > MBSTH ) THEN
                     CONU(I,KK) = ( MU(I,KKP1)*CONU(I,KKP1)*FISG(I,KK)
     &                          +   EU(I,KK)  *CMIX(I,KK)  *DP(I,KK)
     &                          -   DU(I,KK)  *CMIX(I,KK)  *DP(I,KK)  )
     &                          /   PLUMEIN
                  ENDIF             

               ELSE    

                  !-----------------------------------------------------
                  ! Normal condition; so just mix up EU and MU   
                  !-----------------------------------------------------

                  ! Air mass flux going into layer KK of the CLOUD
                  PLUMEIN = MU(I,KKP1) + ( EU(I,KK) * DP(I,KK) )

                  ! Compute concentration
                  IF ( PLUMEIN > MBSTH ) THEN
                     CONU(I,KK) = ( MU(I,KKP1)*CONU(I,KKP1)*FISG(I,KK)
     &                          +   EU(I,KK)  *CMIX(I,KK)  *DP(I,KK)  )
     &                          / PLUMEIN
                  ENDIF
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
               ! ND14 Diagnostic: net upward flux of tracer [kg/s] in 
               ! cloud convection ("CV-FLX-$") (swu, bmy, 6/12/06) 
               ! 
               ! The ND14 diagnostic consists of 4 terms (T1..T4):
               ! -------------------------------------------------------
               ! 
               ! T1: + Mass flux going upward from level K --> K-1 
               !       (notice that the vertical levels are flipped)
               !
               ! T2: - Mass flux going downward from level K-1 --> K
               !       due to large scale subsidence
               !
               ! T3: - Mass flux going downward from level K-1 --> K 
               !       associated with the downdraft plume
               !
               ! T4: + Mass flux going up from level K --> K-1 due 
               !       to enviromental compensation for the downdraft. 
               !
               ! These terms are computed as follows:
               ! -------------------------------------------------------
               !
               ! AIRFLUX:  MU(I,K) * AREA_M2 * 100 / G0  
               !            = air mass (upward) flux in kg/s
               !
               ! T1:      +AIRFLUX * CONU(I,K) 
               !            = tracer mass upward flux [kg/s]
               !
               ! T2:      -AIRFLUX * CMIX(I,K-1)
               !            = subsidence of tracer [kg/s]
               !
               ! T3:      -AIRFLUX * CMIX(I,K)
               !            = downdraft flux of tracer [kg/s]
               !
               ! T4:      +AIRFLUX * COND(I,K-1)
               !            = compensatory upward tracer flux [kg/s] 
               !
               ! Where:
               !
               ! MU       = Grid box surface area         [hPa/s   ] 
               ! AREA_M2  = Mixing ratio in updraft       [m2      ]
               ! CONU     = Mixing ratio in updraft       [kg/kg   ]
               ! COND     = Mixing ratio in downdraft     [kg/kg   ]
               ! CMIX     = Gathered tracer array         [kg/kg   ]
               ! G0     = Acceleration due to gravity   [m/s2    ]
               ! D_NSTEP  = # of convection timesteps     [unitless]
               !
               ! Dividing by the number of time steps within each 
               ! convection step is simply accounting for the scale 
               ! factors (SCALECONV) in diag3.f.
               !========================================================

               ! Only save ND14 if it's turned on
               IF ( ND14 > 0 ) THEN 

                  ! GEOS-Chem lon, lat, alt indices
                  II = IDEEP(I)
                  JJ = J
                  LL = LLPAR - K + 1

                  ! Grid box surface area [m]
                  AREA_M2 = GET_AREA_M2( II, JJ, LL ) 

                  ! Only save from levels 1..LD14 
                  IF ( LL < LD14 ) THEN

                     ! Net upward convective flux [kg/s]
                     CONVFLUP(II,JJ,LL,M) = CONVFLUP(II,JJ,LL,M)  
     &                + (   MU(I,K)   * ( CONU(I,K) - CMIX(I,KM1) )       
     &                    - MD(I,KM1) * ( CMIX(I,K) - COND(I,KM1) ) )      
     &                  * AREA_M2 * g0_100 / D_NSTEP

                  ENDIF
               ENDIF 

               !========================================================
               ! ND38 Diagnostic: loss of soluble tracer [kg/s] to
               ! convective rainout ("WETDCV-$") (swu, bmy, 12/17/03)
               !
               ! The loss of soluble tracer is given by (cf ND14):
               !
               !   MU(I,K+1) * AREA_M2 * 100 / G0  
               !      = Air mass (upward) flux from level K+1 -> K
               !        (Note that vertical levels are reversed)
               !
               ! * CONU(I,K+1) * ( 1 - FISG(I,K) )
               !      = Tracer mass upward from level K+1 -> K [kg/s]
               !
               ! Where:
               !
               ! CONU(I,K+1)   = Tracer in mass flux from K+1 -> K
               ! 1 - FISG(I,K) = Fraction of tracer lost in convective
               !                  updraft going from level K+1 -> K
               !
               !========================================================

               ! Soluble tracer index
               NN = INDEXSOL(M)

               ! Only save to ND38 if it's turned on, if there are soluble 
               ! tracers, and if we are below the LD38 level limit
#if defined( BPCH_DIAG )
               IF ( ( ND38 > 0 .or. Input_Opt%LGTMM ) .and. 
     &              NN > 0 ) THEN

                  ! GEOS-CHEM lon, lat, alt indices
                  II = IDEEP(I)
                  JJ = J
                  LL = LLPAR - K + 1
                  
                  ! Only save up to LD38 vertical levels
                  IF ( LL <= LD38 ) THEN
                  
                     ! Grid box surface area [m2] 
                     AREA_M2 = GET_AREA_M2( II, JJ, LL )  

                     ! Save loss in [kg/s]
                     AD38(II,JJ,LL,NN) = AD38(II,JJ,LL,NN) +
     &                              +  MU(I,KP1) *  CONU(I,KP1) 
     &                                 * ( 1e+0_fp - FISG(I,K) )
     &                                 * AREA_M2 * g0_100 / D_NSTEP  

                  ENDIF
               ENDIF
#endif

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
               !               
               !========================================================

               ! If this is a Hg simulation ...
               IF ( IS_Hg ) THEN

                  ! ... and if this is one of the Hg2 tracers
                  if ( SpcInfo%Is_Hg2 ) THEN

                     ! GEOS-CHEM lon & lat indices
                     II      = IDEEP(I)
                     JJ      = J

                     ! Grid box surface area [m2] 
                     AREA_M2 = GET_AREA_M2( II, JJ, 1 )

                     ! Hg2 wet-scavenged out of the column [kg]
                     WET_Hg2 = MU(I,KP1) * CONU(I,KP1)
     &                          * ( 1e+0_fp - FISG(I,K) ) 
     &                          * AREA_M2 * g0_100 / DELT 

                     ! Category # for this Hg2 species
                     Hg_Cat = SpcInfo%Hg_Cat

                     ! Pass to "ocean_mercury_mod.f"
                     CALL ADD_Hg2_WD      ( II, J, Hg_Cat, WET_Hg2   )
                     CALL ADD_Hg2_SNOWPACK( II, J, Hg_Cat, WET_Hg2, 
     &                                                     State_Met )
                  ENDIF
                  !eck/eds
                  IF ( SpcInfo%Is_HgP ) THEN

                     ! GEOS-CHEM lon & lat indices
                     II      = IDEEP(I)
                     JJ      = J

                     ! Grid box surface area [m2] 
                     AREA_M2 = GET_AREA_M2( II, JJ, 1 )

                     ! Hg2 wet-scavenged out of the column [kg]
                     WET_HgP = MU(I,KP1) * CONU(I,KP1)
     &                         * ( 1e+0_fp - FISG(I,K) )
     &                          * AREA_M2 * g0_100 / DELT 

                     ! Category # for this Hg2 species
                     Hg_Cat = SpcInfo%Hg_Cat

                     ! Pass to "ocean_mercury_mod.f"
                     CALL ADD_HgP_WD      ( II, J, Hg_Cat, WET_HgP   )
                     CALL ADD_Hg2_SNOWPACK( II, J, Hg_Cat, WET_HgP,
     &                                                     State_Met )
                  ENDIF
               ENDIF

               NETFLUX = FLUXIN - FLUXOUT
               
               IF ( ABS(NETFLUX) < MAX(FLUXIN,FLUXOUT)*TINYNUM) THEN
                  NETFLUX = 0.e+0_fp
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
                     NETFLUX = 0.e+0_fp
                  ENDIF

                  DCONDT(I,K) = NETFLUX / DP(I,K)

               ELSE IF ( K > MX(I) ) THEN

                  DCONDT(I,K) = 0.e+0_fp

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
               IF ( QN < 0e+0_fp ) THEN
                  QN = 0e+0_fp
               ENDIF            

               Q(IDEEP(I),J,K,M) = QN
            ENDDO 
         ENDDO    

         ! Free pointer to species database
         SpcInfo => NULL()

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
!  (1 ) N      (INTEGER)   : Number of elements in ARRAY
!  (2 ) ARRAY  (REAL(fp) ) : 1-D vector to be examined
!  (3 ) INC    (INTEGER)   : Increment for stepping thru ARRAY
!  (4 ) TARGET (REAL(fp) ) : Value that ARRAY will be tested against
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
      REAL(fp),  INTENT(IN)  :: ARRAY(N), TARGET
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
