! $Id: convection_mod.f,v 1.5 2004/04/19 15:09:51 bmy Exp $
      MODULE CONVECTION_MOD
!
!******************************************************************************
!  Module CONVECTION_MOD contains routines which select the proper convection
!  code for GEOS-1, GEOS-STRAT, GEOS-3, or GEOS-4 met field data sets. 
!  (bmy, 6/28/03, 2/23/04)
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_CONVECTION       : Wrapper routine, chooses correct convection code
!  (2 ) NFCLDMX             : Convection routine for GEOS-1, GEOS-S, GEOS-3
!
!  GEOS-CHEM modules referenced by convection_mod.f
!  ============================================================================
!  (1 ) dao_mod.f           : Module containing arrays for DAO met fields    
!  (2 ) diag_mod.f          : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) fvdas_convect_mod.f : Module w/ convection code for fvDAS met fields
!  (4 ) grid_mod.f          : Module containing horizontal grid information
!  (5 ) pressure_mod.f      : Module containing routines to compute P(I,J,L)
!  (6 ) time_mod.f          : Module containing routines for computing time
!  (7 ) wetscav_mod.f       : Module containing routines for wetdep/scavenging
!
!  NOTES:
!  (1 ) Contains new updates for GEOS-4/fvDAS convection.  Also now references
!        "error_mod.f".  Now make F in routine NFCLDMX a 4-D array to avoid
!        memory problems on the Altix. (bmy, 1/27/04)
!  (2 ) Bug fix: Now pass NTRACE elements of TCVV to FVDAS_CONVECT in routine 
!        DO_CONVECTION (bmy, 2/23/04)   
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_CONVECTION
!
!******************************************************************************
!  Subroutine DO_CONVECTION is a wrapper for the GEOS-CHEM convection code.
!  It calls the proper convection routine for GEOS-1, GEOS-S, GEOS-3, or 
!  GEOS-4 met field types. (bmy, 6/26/03, 2/23/04)
!
!  NOTES:
!  (1 ) Now updated for GEOS-4/fvDAS convection.  Also placed calls to 
!         DEBUG_MSG for debugging msg output. (swu, bmy, 1/27/04)
!  (2 ) Only pass NTRACE elements of TCVV to FVDAS_CONVECT to prevent array
!        out-of-bounds error for less than NTRACE tracers (bmy, 2/23/04)
!******************************************************************************
!     
      ! References to F90 modules
      USE DAO_MOD,           ONLY : HKETA, HKBETA, ZMEU, ZMMU, ZMMD
      USE DIAG_MOD,          ONLY : AD37
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FVDAS_CONVECT_MOD, ONLY : INIT_FVDAS_CONVECT, FVDAS_CONVECT
      USE WETSCAV_MOD,       ONLY : COMPUTE_F
      USE TIME_MOD,          ONLY : GET_TS_CONV
      USE PRESSURE_MOD              

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT, TCVV, NTRACE, LPRT
#     include "CMN_DIAG"  ! ND37, LD37 

#if   defined( GEOS_4 )
      
      ! More local variables 
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I, ISOL, J, L, L2, N, NSTEP  
      INTEGER       :: INDEXSOL(NTRACE) 
      INTEGER       :: CONVDT    
      REAL*8        :: F(IIPAR,JJPAR,LLPAR,NTRACE)
      REAL*8        :: RPDEL(IIPAR,JJPAR,LLPAR)
      REAL*8        :: DP(IIPAR,JJPAR,LLPAR)
      REAL*8        :: TEMP_DIAG
      REAL*8        :: P1, P2, TDT   

      ! Declare temporary arrays to hold flipped met fields.  This means 
      ! that we won't have to reflip them after convection. (bmy, 12/12/03)
      REAL*8        :: HKBETA_F(IIPAR,JJPAR,LLPAR)
      REAL*8        :: HKETA_F(IIPAR,JJPAR,LLPAR)
      REAL*8        :: ZMEU_F(IIPAR,JJPAR,LLPAR)
      REAL*8        :: ZMMD_F(IIPAR,JJPAR,LLPAR)
      REAL*8        :: ZMMU_F(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! DO_CONVECTION begins here! (for GEOS-4!)
      !=================================================================
      
      ! Convection timestep [s]
      CONVDT = GET_TS_CONV() * 60d0 
       
      ! NSTEP is the # of internal convection timesteps.  According to
      ! notes in the old convection code, 300s works well. (swu, 12/12/03)
      NSTEP  = CONVDT / 300    
      NSTEP  = MAX( NSTEP, 1 ) 

      ! TIMESTEP*2; will be divided by 2 before passing to CONVTRAN 
      TDT    = DBLE( CONVDT ) *2.0D0 / DBLE( NSTEP )

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_FVDAS_CONVECT
         FIRST = .FALSE.
      ENDIF
         
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_CONVECTION: a INIT_FV' )

      !=================================================================
      ! Before calling convection, compute the fraction of insoluble
      ! tracer (Finsoluble) lost in updrafts.  Finsoluble = 1-Fsoluble.
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, ISOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, NTRACE

         ! Get fraction of tracer scavenged and the soluble tracer 
         ! index (ISOL). For non-soluble tracers, F=0 and ISOL=0.
         CALL COMPUTE_F( N, F(:,:,:,N), ISOL ) 
         
         ! Store ISOL in an array for later use
         INDEXSOL(N) = ISOL

         ! Loop over grid boxes
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! ND37 diagnostic: store fraction of tracer 
            ! lost in moist convective updrafts ("MC-FRC-$")
            IF ( ND37 > 0 .and. ISOL > 0 .and. L <= LD37 ) THEN
               AD37(I,J,L,ISOL) = AD37(I,J,L,ISOL) + F(I,J,L,N) 
            ENDIF

            ! GEOS-4 convection routines need the insoluble fraction
            F(I,J,L,N) = 1d0 - F(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
       
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_CONVECTION: a COMPUTE_F' )

      !=================================================================
      ! Compute pressure thickness arrays DP and RPDEL
      ! These arrays are indexed from atm top --> surface
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, L2, P1, P2 )
      DO L = 1, LLPAR

         ! L2 runs from the atm top down to the surface
         L2 = LLPAR - L + 1

         ! Loop over surface grid boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR
               
            ! Pressure at bottom and top edges of grid box [hPa]
            P1 = GET_PEDGE(I,J,L)
            P2 = GET_PEDGE(I,J,L+1)

            ! DP = Pressure difference between top & bottom edges [Pa]
            DP(I,J,L2) = ( P1 - P2 ) * 100.0d0

            ! RPDEL = reciprocal of DP [1/hPa]
            RPDEL(I,J,L2) = 100.0d0 / DP(I,J,L2) 
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_CONVECTION: a DP, RPDEL' )
   
      !=================================================================
      ! Flip arrays in the vertical and call FVDAS_CONVECT
      !=================================================================

      ! Flip 4-D arrays for FVDAS_CONVECT
      STT(:,:,1:LLPAR,1:NTRACE) = STT(:,:,LLPAR:1:-1,1:NTRACE)
      F  (:,:,1:LLPAR,1:NTRACE) = F  (:,:,LLPAR:1:-1,1:NTRACE)

      ! Flip met fields for FVDAS_CONVECT.  Store them in temporary
      ! local arrays so that we don't have to reflip them,  This
      ! is more efficient and saves clock cycles. (bmy, 12/16/03)
      HKBETA_F(:,:,1:LLPAR) = HKBETA(:,:,LLPAR:1:-1)
      HKETA_F (:,:,1:LLPAR) = HKETA (:,:,LLPAR:1:-1)
      ZMEU_F  (:,:,1:LLPAR) = ZMEU  (:,:,LLPAR:1:-1)
      ZMMD_F  (:,:,1:LLPAR) = ZMMD  (:,:,LLPAR:1:-1)
      ZMMU_F  (:,:,1:LLPAR) = ZMMU  (:,:,LLPAR:1:-1)

      ! Call the fvDAS convection routines (originally from NCAR!)
      CALL FVDAS_CONVECT( TDT,    NTRACE,  STT(:,:,:,1:NTRACE),    
     &                    RPDEL,  HKETA_F, HKBETA_F, 
     &                    ZMMU_F, ZMMD_F,  ZMEU_F,  
     &                    DP,     NSTEP,   F,      
     &                    TCVV(1:NTRACE),  INDEXSOL )

      ! Reflip STT array after FVDAS_CONVECT
      STT(:,:,1:LLPAR,1:NTRACE) = STT(:,:,LLPAR:1:-1,1:NTRACE)

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_CONVECTION: a FVDAS_CONVECT' )

#else

      !==================================================================
      ! DO_CONVECTION begins here! (for GEOS-3)
      !==================================================================

      ! Call the S-J Lin convection routine for GEOS-1, GEOS-S, GEOS-3
      CALL NFCLDMX( NTRACE, TCVV(1:NTRACE), STT(:,:,:,1:NTRACE) )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_CONVECTION: a NFCLDMX' )

#endif

      ! Return to calling program
      END SUBROUTINE DO_CONVECTION

!------------------------------------------------------------------------------

      SUBROUTINE NFCLDMX( NC, TCVV, Q )
!
!******************************************************************************
!  Subroutine NFCLDMX is S-J Lin's cumulus transport module for 3D GSFC-CTM,
!  modified for the GEOS-CHEM model.  The "NF" stands for "no flipping", and
!  denotes that you don't have to flip the tracer array Q in the main
!  program before passing it to NFCLDMX. (bmy, 2/12/97, 1/26/04)
!
!  NOTE: NFCLDMX can be used with GEOS-1, GEOS-STRAT, and GEOS-3 met fields.
!  For GEOS-4/fVDAS, you must use the routines in "fvdas_convect_mod.f"
!  (bmy, 6/26/03)
!
!  Arguments as input:
!  ==========================================================================
!  (1 ) NC    : TOTAL number of tracers (soluble + insoluble) [unitless]
!  (2 ) TCVV  : MW air (g/mol) / MW of tracer (g/mol)         [unitless]
!
!  Arguments as Input/Output:
!  ============================================================================
!  (3 )  Q    : Tracer concentration                          [v/v]
!
!  Reference:
!  ============================================================================
!  Lin, SJ.  "Description of the parameterization of cumulus transport
!     in the 3D Goddard Chemistry Transport Model, NASA/GSFC, 1996.
!
!  Vertical indexing:
!  ============================================================================
!  The indexing of the vertical sigma levels has been changed from 
!  SJ-Lin's original code:
!
!                 Old Method          New Method
!                  (SJ Lin)    
!
!               ------------------------------------- Top of Atm.
!                  k = 1               k = NLAY
!               ===================================== Max Extent 
!                  k = 2               k = NLAY-1      of Clouds
!               -------------------------------------
!              
!                   ...                 ...             
!
!               -------------------------------------
!                  k = NLAY-3          k = 4
!               -------------------------------------
!                  k = NLAY-2          k = 3
!               ------------------------------------- Cloud base
!                  k = NLAY-1          k = 2
!               -    -    -    -    -    -    -    - 
!                  k = NLAY            k = 1
!               ===================================== Ground
!
!      which means that:
!
!                 Old Method                      New Method
!                  (SJ Lin)
!
!            k-1      ^                      k+1      ^
!            ---------|---------             ---------|---------
!                     |                               |
!                  CLDMAS(k)                       CLDMAS(k)
!             
!                                 becomes
!            k      DTRN(k),                 k      DTRN(k),       
!                 QC(k), Q(k)                     QC(k), Q(k)   
!          
!                     ^                               ^
!            ---------|---------             ---------|---------
!                     |                               |   
!            k+1   CLDMAS(k+1)               k-1   CLDMAS(k-1)
!
!
!      i.e., the lowest level    used to be  NLAY  but is now  1
!            the level below k   used to be  k+1   but is now  k-1.
!            the level above k   used to be  k-1   but is now  k+1
!            the top of the atm. used to be  1     but is now  NLAY.
!
!  The old method required that the vertical dimensions of the CLDMAS, DTRN, 
!  and Q arrays had to be flipped before and after calling CLDMX.  Also, 
!  diagnostic arrays generated within CLDMX also had to be flipped.  The new 
!  indexing eliminates this requirement (and also saves on array operations).  
!
!  Major Modifications:
!  ============================================================================
!  Original Author:   Shian-Jiann Lin, Code 910.3, NASA/GSFC
!  Original Release:  12 February 1997
!                     Version 3, Detrainment and Entrainment are considered.
!                     The algorithm reduces to that of version 2 if Dtrn = 0.
! 
!  Modified By:       Bob Yantosca, for Harvard Atmospheric Sciences
!  Modified Release:  27 January 1998
!                     Version 3.11, contains features of V.3 but also 
!                     scavenges soluble tracer in wet convective updrafts.
!                     
!                     28 April 1998
!                     Version 3.12, now includes mass flux diagnostic
!         
!                     11 November 1999
!                     Added mass-flux diagnostics
!
!                     04 January 2000
!                     Updated scavenging constant AS2
!
!                     14 March 2000
!                     Added new wet scavenging code and diagnostics
!                     based on the GMI algorithm
!
!                     02 May 2000
!                     Added parallel loop over tracers
!
!  NOTES:              
!  (1 ) NFCLDMX is written in Fixed-Form Fortran 90.
!  (2 ) Added TCVV to the argument list.  Also cleaned up argument
!        and local variable declarations. (bey, bmy, 11/10/99)
!  (3 ) AD38 and CONVFLUP are now declared allocatable in "diag_mod.f". 
!        (bmy, 11/29/99)
!  (4 ) Bug fix for tagged CO tracer run (bey, bmy, 1/4/00)
!  (5 ) Add new routines for computing scavenging coefficients,
!        as well as adding the AD37 diagnostic array. (bmy, 3/14/00)
!  (6 ) Updated comments (bmy, 10/2/01)
!  (7 ) Now print a header to stdout on the first call, to confirm that 
!        NFCLDMX has been called (bmy, 4/15/02)
!  (8 ) Remove PZ from the arg list -- it isn't used! (bmy, 8/22/02)
!  (9 ) Fixed ND38 diagnostic so that it now reports correctly (must divide
!        by DNS).  Updatec comments, cosmetic changes. (bmy, 1/27/03)
!  (10) Bug fix: remove duplicate K from PRIVATE declaration (bmy, 3/23/03)
!  (11) Now removed all arguments except NC, TCVV, Q from the arg list -- the
!        other arguments can be supplied via F90 modules.  Now references
!        "dao_mod.f", "grid_mod.f", "pressure_mod.f", and "time_mod.f".
!        (bmy, 3/27/03)
!  (12) Bundled into "convection_mod.f" (bmy, 6/26/03)
!  (13) Make sure K does not go out of bounds in ND38 diagnostic.  Now make 
!        F a 4-D array in order to avoid memory problems on the Altix.  
!        (bmy, 1/27/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD,   CLDMAS, DTRN=>DTRAIN
      USE DIAG_MOD,     ONLY : AD37, AD38,   CONVFLUP
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE PRESSURE_MOD, ONLY : GET_BP
      USE TIME_MOD,     ONLY : GET_TS_CONV
      USE WETSCAV_MOD,  ONLY : COMPUTE_F

      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches & arrays

      ! Arguments
      INTEGER, INTENT(IN)    :: NC 
      REAL*8,  INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NC)
      REAL*8,  INTENT(IN)    :: TCVV(NC)

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, K, KTOP, L, N, NDT
      INTEGER                :: IC, ISTEP, JUMP, JS, JN, NS
      INTEGER                :: IMR, JNP, NLAY
      REAL*8,  SAVE          :: DSIG(LLPAR)
      REAL*8                 :: SDT, CMOUT, ENTRN, DQ, AREA_M2
      REAL*8                 :: T0, T1, T2, T3, T4, TSUM, DELQ
      REAL*8                 :: DTCSUM(IIPAR,JJPAR,LLPAR,NNPAR)

      ! F is the fraction of tracer lost to wet scavenging in updrafts
      REAL*8                 :: F(IIPAR,JJPAR,LLPAR,NC)

      ! Local Work arrays
      REAL*8                 :: BMASS(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: QB(IIPAR,JJPAR)
      REAL*8                 :: MB(IIPAR,JJPAR)
      REAL*8                 :: QC(IIPAR,JJPAR) 

      ! TINY = a very small number
      REAL*8, PARAMETER      :: TINY = 1d-14 

      ! ISOL is an index for the diagnostic arrays
      INTEGER                :: ISOL

      ! QC_PRES and QC_SCAV are the amounts of tracer 
      ! preserved against and lost to wet scavenging
      REAL*8                 :: QC_PRES, QC_SCAV 

      ! DNS is the double precision value for NS
      REAL*8                 :: DNS
     
      !=================================================================
      ! NFCLDMX begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'N F C L D M X  -- by S-J Lin'
         WRITE( 6, '(a)' ) 'Modified for GEOS-CHEM by Bob Yantosca'
         WRITE( 6, '(a)' ) 'Last Modification Date: 1/27/04'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! DSIG is the sigma-level thickness (NOTE: this assumes that
         ! we are using a pure-sigma grid.  Use new routine for fvDAS.)
         DO L = 1, LLPAR        
            DSIG(L) = GET_BP(L) - GET_BP(L+1)
         ENDDO

         ! Reset first time flag
         FIRST = .FALSE.
      ENDIF
      
      ! Define dimensions
      IMR  = IIPAR
      JNP  = JJPAR
      NLAY = LLPAR

      ! Convection timestep [s]
      NDT  = GET_TS_CONV() * 60d0

      !=================================================================
      ! Define active convective region, from J = JS(outh) to 
      ! J = JN(orth), and to level K = KTOP. 
      !
      ! Polar regions are too cold to have moist convection.
      ! (Dry convection should be done elsewhere.)
      !
      ! We initialize the ND14 diagnostic each time we start a new 
      ! time step loop.  Only initialize DTCSUM array if the ND14 
      ! diagnostic is turned on.  This saves a quite a bit of time. 
      ! (bmy, 12/15/99)       
      !=================================================================
      IF ( ND14 > 0 ) DTCSUM = 0d0

      KTOP = NLAY - 1
      JUMP = (JNP-1) / 20
      JS   = 1 + JUMP
      JN   = JNP - JS + 1

      !=================================================================
      ! Internal time step for convective mixing is 300 sec.
      ! Doug Rotman (LLNL) says that 450 sec works just as well.       
      !=================================================================
      NS  = NDT / 300
      NS  = MAX(NS,1)
      SDT = FLOAT(NDT) / FLOAT(NS)
      DNS = DBLE( NS )

!=============================================================================
!  BMASS has units of kg/m^2 and is equivalent to AD(I,J,L) / AREA_M2
!
!   Ps - Pt (mb)| P2 - P1 | 100 Pa |  s^2  | 1  |  1 kg        kg
!  -------------+---------+--------+-------+----+--------  =  -----
!               | Ps - Pt |   mb   | 9.8 m | Pa | m^2 s^2      m^2
!
!  This is done to keep BMASS in the same units as CLDMAS * SDT
!
!  We can parallelize over levels here.  The only quantities that need to 
!  be held local are the loop counters (I, IC, J, JREF, K). (bmy, 5/2/00)
!
!  Now use routine GET_AREA_M2 from "grid_mod.f" to get surface area of
!  grid boxes in m2. (bmy, 2/4/03)
!=============================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_M2, K )
!$OMP+SCHEDULE( DYNAMIC )
      DO K = 1, NLAY
         DO J = 1, JJPAR
            AREA_M2 = GET_AREA_M2( J )
            DO I = 1, IMR
               BMASS(I,J,K) = AD(I,J,K) / AREA_M2
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! (1)  T r a c e r   L o o p 
      !
      ! We now parallelize over tracers, since tracers are independent 
      ! of each other.  The parallel loop only takes effect if you 
      ! compile with the f90 "-mp" switch.  Otherwise the compiler will 
      ! interpret the parallel-processing directives as comments, and 
      ! the loop will execute on a single thread.
      !
      ! The following types of quantities must be held local for 
      ! parallelization:
      ! (1) Loop counters ( I, IC, ISTEP, J, K )
      ! (2) Scalars that are assigned values inside the tracer loop: 
      !     ( CMOUT, DELQ, ENTRN, ISOL, QC_PRES, etc. )
      ! (3) Arrays independent of tracer ( F, MB, QB, QC )
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( CMOUT, DELQ, ENTRN, I, IC, ISOL, ISTEP, J, AREA_M2, K  ) 
!$OMP+PRIVATE( MB, QB, QC, QC_PRES, QC_SCAV, T0, T1, T2, T3, T4, TSUM )
!$OMP+SCHEDULE( DYNAMIC )
      DO IC = 1, NC

         !==============================================================
         ! (2)  S c a v e n g i n g   i n   C l o u d   U p d r a f t s
         !
         ! Call COMPUTE_F to compute the fraction of tracer scavenged 
         ! in convective cloud updrafts.  COMPUTE_F works for both full 
         ! chemistry (NSRCX == 3) and Rn-Pb-Be chemistry (NSRCX == 1) 
         ! simulations.  It is best to compute the fraction of tracer 
         ! scavenged at this point, outside the internal time step loop.  
         ! This will avoid having to repeat the entire calculation of F 
         ! for NS times.
         !
         ! ISOL, which is returned from COMPUTE_F, is the tracer index 
         ! used for diagnostic arrays AD37 and AD38.  ISOL = 0 for all 
         ! non-soluble tracers.   
         !==============================================================
         CALL COMPUTE_F( IC, F(:,:,:,IC), ISOL ) 

         ! ND37 diagnostic -- store F only for soluble tracers
         IF ( ND37 > 0 .and. ISOL > 0 ) THEN
            DO K = 1, LD37
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               AD37(I,J,K,ISOL) = AD37(I,J,K,ISOL) + F(I,J,K,IC) 
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      
         !==============================================================
         ! (3)  I n t e r n a l   T i m e   S t e p   L o o p
         !
         ! The internal time step is currently 300 seconds.  
         !==============================================================
         DO ISTEP = 1, NS

            !===========================================================
            ! (4)  B e l o w   C l o u d   B a s e  (K < 3)
            !
            ! Loop over longitude and latitude (I,J), and consider what
            ! is below the cloud base (i.e. below the 2nd sigma level).  
            !===========================================================
            DO J = JS, JN
            DO I = 1, IMR

               !========================================================
               ! (4.1) If Cloud Mass Flux exists at (I,J,2), 
               !       then compute QB.
               !
               ! QB is "weighted average" mixing ratio below the cloud 
               ! base.  QB is used to compute QC, which is the mixing 
               ! ratio of the air that moved in cumulus transport up 
               ! to the next level.  MB is the total mass of air below 
               ! the cloud base.
               !========================================================

               IF ( CLDMAS(I,J,2) .gt. TINY ) THEN 
                  QB(I,J) = 
     &               ( Q(I,J,1,IC) * DSIG(1)   + 
     &                 Q(I,J,2,IC) * DSIG(2) ) / ( DSIG(1) + DSIG(2) )

                  MB(I,J) = BMASS(I,J,1) + BMASS(I,J,2)

!=============================================================================
!         Total mass of tracer below cloud base (i.e. MB  * QB ) +   
!         Subsidence into cloud base from above (i.e. SDT * C(2) * Q(3) ) 
!  QC =  -----------------------------------------------------------------
!             Total air mass below cloud base (i.e. MB + C(2)*Q(3) )
!=============================================================================
                  QC(I,J) = 
     &               ( MB(I,J)       * QB(I,J)       + 
     &                 CLDMAS(I,J,2) * Q(I,J,3,IC)   * SDT ) /
     &               ( MB(I,J)       + CLDMAS(I,J,2) * SDT )
     
                  !=====================================================
                  ! DQ = QB - QC is the total mass to be transported 
                  ! out of the cloud base.  Changes below cloud base 
                  ! are proportional to the background mass.
                  ! 
                  ! Subtract DQ from Q(*,*,K=1,*) and from Q(*,*,K=2,*), 
                  ! but do not make Q(*,*,K=1,*) or Q(*,*,K=2,*) < 0.
                  !=====================================================
!-----------------------------------------------------------------------------
!  Prior to 1/4/00:
!  This ensures additivity when using tagged CO tracers (bey, bmy, 1/4/00)
!                  DQ = QB(I,J) - QC(I,J)
!                  
!                  IF ( DQ .GT. Q(I,J,1,IC)   .OR. 
!     &                 DQ .GT. Q(I,J,2,IC) ) THEN
!                     Q(I,J,2,IC) = QC(I,J)
!                     Q(I,J,1,IC) = QC(I,J)
!                  ELSE
!                     Q(I,J,2,IC) = Q(I,J,2,IC) - DQ
!                     Q(I,J,1,IC) = Q(I,J,1,IC) - DQ
!                  ENDIF
!-----------------------------------------------------------------------------
                  Q(I,J,2,IC) = QC(I,J)
                  Q(I,J,1,IC) = QC(I,J)

               !========================================================
               ! If there is no Cloud mass flux, set QC = Q(K=3) 
               ! at this I,J location
               !========================================================
               ELSE
                  QC(I,J) = Q(I,J,3,IC)
               ENDIF
            ENDDO !I
            ENDDO !J

            !===========================================================
            ! (5)  A b o v e   C l o u d   B a s e
            !
            ! Loop over sigma levels K and longitude and latitude (I,J) 
            !
            ! Recall that K+1 is the level BELOW level K, and that K-1 
            ! is the level ABOVE level K.  This is the way that SJ Lin 
            ! indexes the sigma levels.
            !===========================================================
            DO K = 3, KTOP
            DO J = JS, JN

               ! Grid box surface area [m2]
               AREA_M2 = GET_AREA_M2( J )

               DO I = 1, IMR

!==============================================================================
!  (5.1)  M a s s   B a l a n c e   i n   C l o u d  ===>  QC(I,J)
!
!  QC_PRES = QC(I,J,K-1) * AP = amt of Qc preserved against wet scavenging
!  QC_SCAV = QC(I,J,K-1) * AS = amt of Qc lost to wet scavenging
!  CMOUT   = air mass flowing out of cloud at level K 
!  ENTRN   = Entrainment: air mass flowing into cloud at level K   
!
!  If ENTRN > 0 then compute the new value of QC(I,J):
!
!                CLDMAS(K-1) * QC_PRES  +  ENTRN(K) * Q(K)
!    QC(I,J) =  -------------------------------------------
!                      CLDMAS(I,J,K) + DTRN(I,J,K)
!
!            =   tracer mass coming in from below      (i.e. level K-1) + 
!                tracer mass coming in from this level (i.e. level K)
!               -----------------------------------------------------------
!                             total mass coming into cloud
!
!  Otherwise, preserve the previous value of QC(I,J).  This will ensure
!  that TERM1 - TERM2 is not a negative quantity (see below).
!  
!  Entrainment must be >= 0 (since we cannot have a negative flux of air
!  into the cloud).  This condition is strong enough to ensure that
!  CMOUT > 0 and will prevent floating-point exception.
!==============================================================================

                  IF ( CLDMAS(I,J,K-1) .gt. TINY ) THEN
                     CMOUT    = CLDMAS(I,J,K) + DTRN(I,J,K)
                     ENTRN    = CMOUT         - CLDMAS(I,J,K-1)
                     QC_PRES  = QC(I,J) * ( 1d0 - F(I,J,K,IC) )
                     QC_SCAV  = QC(I,J) * F(I,J,K,IC)

                     IF ( ENTRN .ge. 0 ) THEN
                        QC(I,J) = ( CLDMAS(I,J,K-1) * QC_PRES       + 
     &                              ENTRN           * Q(I,J,K,IC) ) / 
     &                              CMOUT
                     ENDIF   

!==============================================================================
!  (5.2)  M a s s   B a l a n c e   i n   L e v e l  ===>  Q(I,J,K,IC)
!
!  The cumulus transport above the cloud base is done as follows:
!     C_k-1  = cloud air mass flux from level k-1 to level k
!     C_k    = cloud air mass flux from level k   to level k+1
!     QC_k-1 = mixing ratio of tracer INSIDE CLOUD at level k-1
!     QC_k   = mixing ratio of tracer INSIDE CLOUD at level k
!     Q_k    = mixing ratio of tracer in level k
!     Q_k+1  = mixing ratio of tracer in level k+1
!                                
!                       |                    |
!      k+1     ^        |       Cloud        |3)      C_k * Q_k+1
!              |        |         ^          |            |
!      --------|--------+---------|----------+------------|--------
!              |        |         |          |            V
!      k     C_k        |2)   C_k * QC_k     |
!                       |                    |
!                       |                    |
!                       |         ^          |4)    C_k-1 * Q_k           
!              ^        |         |          |            |
!      --------|--------+---------|----------+------------|----------
!              |        |         |          |            |
!      k-1   C_k-1      |1) C_k-1 * QC_k-1   |            V
!                       |         * AP       |
!
!  There are 4 terms that contribute to mass flow in and out of level k:
!
!  1) C_k-1 * QC_PRES = tracer convected from level k-1 to level k 
!  2) C_k   * QC_k    = tracer convected from level k   to level k+1 
!  3) C_k   * Q_k+1   = tracer subsiding from level k+1 to level k 
!  4) C_k-1 * Q_k     = tracer subsiding from level k   to level k-1 
!
!  Therefore the change in tracer concentration is given by
!     DELQ = (Term 1) - (Term 2) + (Term 3) - (Term 4)
!
!  and Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ.  
!
!  The term T0 is the amount of tracer that is scavenged out of the box.
!  Compute that term here for the ND38 diagnostic below. (bmy, 1/27/03)
!==============================================================================
                     T0   =  CLDMAS(I,J,K-1) * QC_SCAV   
                     T1   =  CLDMAS(I,J,K-1) * QC_PRES
                     T2   = -CLDMAS(I,J,K  ) * QC(I,J       )
                     T3   =  CLDMAS(I,J,K  ) * Q (I,J,K+1,IC)
                     T4   = -CLDMAS(I,J,K-1) * Q (I,J,K,  IC)

                     TSUM = T1 + T2 + T3 + T4 

                     DELQ = ( SDT / BMASS(I,J,K) ) * TSUM

                     ! If DELQ > Q then do not make Q negative!!!
                     IF ( Q(I,J,K,IC) + DELQ < 0 ) THEN
                        DELQ = -Q(I,J,K,IC)
                     ENDIF

                     Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ
                                
                     !==================================================
                     ! ND14 Diagnostic: Upward mass flux due to wet 
                     ! convection.  The diagnostic ND14 works only for 
                     ! levels >= than 3.  This is ok for now since I 
                     ! want to have a closed budget with a box starting 
                     ! from the ground.
                     !
                     ! DTCSUM(I,J,K,IC) is the flux (kg/box/sec) in 
                     ! the box (I,J), for the tracer IC going out of 
                     ! the top of the layer K to the layer above (K+1) 
                     ! (bey, 11/10/99).                     
                     !==================================================
                     IF ( ND14 > 0 ) THEN
                        DTCSUM(I,J,K,IC) = DTCSUM(I,J,K,IC) +  
     &                       (-T2-T3) * AREA_M2 / TCVV(IC)  
                     ENDIF

                     !==================================================
                     ! ND38 Diagnostic: loss of soluble tracer to wet
                     ! scavenging in cloud updrafts [kg/s].  We must 
                     ! divide by DNS, the # of internal timesteps.
                     !==================================================
                     IF ( ND38 > 0 .and. ISOL > 0 .and. K <= LD38 ) THEN
                        AD38(I,J,K,ISOL) = AD38(I,J,K,ISOL) +
     &                       ( T0 * AREA_M2 / ( TCVV(IC) * DNS ) )
                     ENDIF
  
                  !==================================================
                  ! No cloud transport if cloud mass flux < TINY; 
                  ! Change Qc to q
                  !==================================================
                  ELSE                     
                     QC(I,J) = Q(I,J,K,IC)
                  ENDIF
               ENDDO  !I
            ENDDO     !J
            ENDDO     !K
         ENDDO        !NSTEP
      ENDDO           !IC
!$OMP END PARALLEL DO

      !=================================================================
      ! ND14 Diagnostic: Store into the CONVFLUP array.  
      ! Also divide by the number of internal timesteps (DNS).
      !================================================================= 
      IF ( ND14 > 0 ) THEN
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, IC, J, K ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO IC = 1,  NC
         DO K  = 3,  KTOP
         DO J  = JS, JN
         DO I  = 1,  IMR
            CONVFLUP(I,J,K,IC) = CONVFLUP(I,J,K,IC) +
     &                           ( DTCSUM(I,J,K,IC) / DNS )
         ENDDO
         ENDDO
         ENDDO
         ENDDO 
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE NFCLDMX

!------------------------------------------------------------------------------

      END MODULE CONVECTION_MOD
