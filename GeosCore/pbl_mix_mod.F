!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pbl_mix_mod
!
! !DESCRIPTION: Module PBL\_MIX\_MOD contains routines and variables used to 
!  compute the planetary boundary layer (PBL) height and to mix tracers 
!  underneath the PBL top.  
!\\
!\\
! !INTERFACE: 
!
      MODULE PBL_MIX_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_PBL_MIX
      PUBLIC  :: DO_PBL_MIX
      PUBLIC  :: GET_FRAC_OF_PBL
      PUBLIC  :: GET_FRAC_UNDER_PBLTOP
      PUBLIC  :: GET_PBL_MAX_L
      PUBLIC  :: GET_PBL_TOP_hPa
      PUBLIC  :: GET_PBL_TOP_L
      PUBLIC  :: GET_PBL_TOP_m
      PUBLIC  :: GET_PBL_THICK
      PUBLIC  :: INIT_PBL_MIX
      PUBLIC  :: COMPUTE_PBL_HEIGHT

#if defined ( DEVEL )
      PUBLIC :: PBL_TOP_L, PBL_TOP_m
#endif

!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: TURBDAY
!
! !REVISION HISTORY:
!  11 Feb 2005 - R. Yantosca - Initial version
!  (1 ) Now modified for GCAP and GEOS-5 met fields (bmy, 5/24/05)
!  (2 ) Remove reference to "CMN" and XTRA2. (bmy, 8/30/05)
!  (3 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (4 ) Add INIT_PBL_MIX and COMPUTE_PBL_HEIGHT as PUBLIC routines 
!        (lin, 5/29/09)
!  (5 ) Extend tracers for APM simulation (GanLuo, 2010)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES
      !
      ! IMIX        : Array for integer # of levels under PBL top
      ! FPBL        : Array for frac # of levels under PBL top
      ! F_OF_PBL    : Array for frac of box (I,J,L) w/in PBL
      ! F_UNDER_TOP : Array for frac of box (I,J,L) under PBL top
      ! PBL_TOP_hPa : Array for PBL top [hPa]
      ! PBL_TOP_L   : Array for PBL top [model levels]
      ! PBL_TOP_m   : Array for PBL top [m]
      ! PBL_THICK   : Array for PBL thickness [hPa]
      !=================================================================

      ! Scalars
      INTEGER              :: PBL_MAX_L
      
      ! Arrays
      INTEGER, ALLOCATABLE :: IMIX(:,:)
      REAL*8,  ALLOCATABLE :: FPBL(:,:)
      REAL*8,  ALLOCATABLE :: F_OF_PBL(:,:,:)
      REAL*8,  ALLOCATABLE :: F_UNDER_TOP(:,:,:)
      REAL*8,  ALLOCATABLE :: PBL_TOP_hPa(:,:)
      REAL*8,  ALLOCATABLE :: PBL_TOP_L(:,:)
      REAL*8,  ALLOCATABLE :: PBL_TOP_m(:,:)
      REAL*8,  ALLOCATABLE :: PBL_THICK(:,:)
      REAL*8,  ALLOCATABLE :: XTRA2(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_pbl_mix
!
! !DESCRIPTION: Subroutine DO\_PBL\_MIX is the driver routine for planetary 
!  boundary layer mixing.  The PBL layer height and related quantities are 
!  always computed.  Complete mixing of tracers underneath the PBL top is 
!  toggled by the DO\_TURBDAY switch.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_PBL_MIX( DO_TURBDAY, Input_Opt, 
     &                       State_Met,  State_Chm  )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
#if   defined( APM ) 
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: DO_TURBDAY  ! =T means call TURBDAY 
                                                   !    for full PBL mixing
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  07 Sep 2011 - G. Luo      - Add modifications for APM
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL, SAVE  :: FIRST = .TRUE.
      INTEGER        :: N_TOT_TRC
      INTEGER        :: N_TRACERS

      ! Arrays
      REAL*8         :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! DO_PBL_MIX begins here!
      !=================================================================

      ! Get fields from the Input_OPt object
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Pick the proper # of tracers
#if   defined( APM )  
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! If APM is used, then convect
#else                                    !  the G-C advected tracers +
      N_TOT_TRC = N_TRACERS              !  the APM aerosol tracers
#endif

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_PBL_MIX
         FIRST = .FALSE.
      ENDIF
         
      ! Compute PBL height and related quantities
      CALL COMPUTE_PBL_HEIGHT( State_Met )

      ! Do complete mixing of tracers in the PBL (if necessary)
      IF ( DO_TURBDAY ) CALL TURBDAY( N_TOT_TRC, STT, TCVV, State_Met )

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE DO_PBL_MIX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_pbl_height
!
! !DESCRIPTION: Subroutine COMPUTE\_PBL\_HEIGHT computes the PBL height and 
!  other related quantities.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COMPUTE_PBL_HEIGHT( State_Met )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PEDGE

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_GCTM_MOD             ! Scale height
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  (1 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (2 ) Remove reference to "CMN" and XTRA2 -- they're obsolete.  Also do not
!        force BLTOP, BLTHIK to minimum values for GEOS-STRAT met fields.
!        (bmy, 8/30/05)
!  (3 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,     J,      L,    LTOP
      REAL*8  :: BLTOP, BLTHIK, DELP
      REAL*8  :: P(0:LLPAR)

      !=================================================================
      ! COMPUTE_PBL_HEIGHT begins here!
      !=================================================================

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, P, BLTOP, BLTHIK, LTOP, DELP )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !----------------------------------------------
         ! Define pressure edges:
         ! P(L-1) = P at bottom edge of box (I,J,L)
         ! P(L  ) = P at top    edge of box (I,J,L)
         !----------------------------------------------

         ! Pressure at level edges [hPa]
         DO L = 0, LLPAR
            P(L) = GET_PEDGE(I,J,L+1)
         ENDDO
         
         !----------------------------------------------
         ! Find PBL top and thickness [hPa]
         !----------------------------------------------

         ! BLTOP = pressure at PBL top [hPa]
         ! Use barometric law since PBL is in [m]
         BLTOP  = P(0) * EXP( -State_Met%PBLH(I,J) / SCALE_HEIGHT )

         ! BLTHIK is PBL thickness [hPa]
         BLTHIK = P(0) - BLTOP

         !----------------------------------------------
         ! Find model level where BLTOP occurs
         !----------------------------------------------
         LTOP = 0

         ! Loop over levels
         DO L = 1, LLPAR

            ! Exit when we get to the PBL top level
            IF ( BLTOP > P(L) ) THEN
               LTOP = L
               EXIT
            ENDIF

         ENDDO 

         !----------------------------------------------
         ! Define various related quantities
         !----------------------------------------------

         ! IMIX(I,J)   is the level where the PBL top occurs at (I,J)
         ! IMIX(I,J)-1 is the number of whole levels below the PBL top
         IMIX(I,J)        = LTOP 

         ! Fraction of the IMIXth level underneath the PBL top
         FPBL(I,J)        = 1d0 - ( BLTOP     - P(LTOP) ) / 
     &                            ( P(LTOP-1) - P(LTOP) ) 

         ! PBL top [model layers]
         PBL_TOP_L(I,J)   = FLOAT( IMIX(I,J) - 1 ) + FPBL(I,J)

         ! PBL top [hPa]
         PBL_TOP_hPa(I,J) = BLTOP

         ! Zero PBL top [m] -- compute below
         PBL_TOP_m(I,J)   = 0d0

         ! PBL thickness [hPa]
         PBL_THICK(I,J)   = BLTHIK

         !==============================================================
         ! Loop up to the maximum tropopause level
         !==============================================================
         DO L = 1, LLTROP

            ! Thickness of grid box (I,J,L) [hPa]
            DELP = P(L-1) - P(L)
         
            IF ( L < IMIX(I,J) ) THEN

               !--------------------------------------------
               ! (I,J,L) lies completely below the PBL top
               !--------------------------------------------

               ! Fraction of grid box (I,J,L) w/in the PBL
               F_OF_PBL(I,J,L)    = DELP / BLTHIK

               ! Fraction of grid box (I,J,L) underneath PBL top
               F_UNDER_TOP(I,J,L) = 1d0

               ! PBL height [m]
               PBL_TOP_m(I,J)     = PBL_TOP_m(I,J) +
     &                              State_Met%BXHEIGHT(I,J,L)

            ELSE IF ( L == IMIX(I,J) ) THEN

               !--------------------------------------------
               ! (I,J,L) straddles the PBL top
               !--------------------------------------------

               ! Fraction of grid box (I,J,L) w/in the PBL
               F_OF_PBL(I,J,L)    = ( P(L-1) - BLTOP ) / BLTHIK

               ! Fraction of grid box (I,J,L) underneath PBL top
               F_UNDER_TOP(I,J,L) = FPBL(I,J) 
               
               ! PBL height [m] 
               PBL_TOP_m(I,J)     = PBL_TOP_m(I,J) + 
     &                            ( State_Met%BXHEIGHT(I,J,L) *
     &                              FPBL(I,J) )

            ELSE

               !--------------------------------------------
               ! (I,J,L) lies completely above the PBL top
               !--------------------------------------------

               ! Fraction of grid box (I,J,L) w/in the PBL
               F_OF_PBL(I,J,L)    = 0d0

               ! Fraction of grid box (I,J,L) underneath PBL top
               F_UNDER_TOP(I,J,L) = 0d0
               
            ENDIF 

!### Debug
!            IF ( I==23 .and. J==34 .and. L < 6 ) THEN 
!               PRINT*, '###--------------------------------------'
!               PRINT*, '### COMPUTE_PBL_HEIGHT'
!               PRINT*, '### I, J, L     : ', I, J, L
!               PRINT*, '### P(L-1)      : ', P(L-1)
!               PRINT*, '### P(L)        : ', P(L)
!               PRINT*, '### PBL(I,J)    : ', PBL(I,J)
!               PRINT*, '### F_OF_PBL    : ', F_OF_PBL(I,J,L)
!               PRINT*, '### F_UNDER_TOP : ', F_UNDER_TOP(I,J,L)
!               PRINT*, '### IMIX        : ', IMIX(I,J)
!               PRINT*, '### FPBL        : ', FPBL(I,J)
!               PRINT*, '### PBL_TOP_hPa : ', PBL_TOP_hPa(I,J)
!               PRINT*, '### PBL_TOP_L   : ', PBL_TOP_L(I,J)
!               PRINT*, '### DELP        : ', DELP
!               PRINT*, '### BLTHIK      : ', BLTHIK
!               PRINT*, '### BLTOP       : ', BLTOP
!               PRINT*, '### BXHEIGHT    : ', BXHEIGHT(I,J,L)
!               PRINT*, '### PBL_TOP_m   : ', PBL_TOP_m(I,J)
!               PRINT*, '### other way m : ', 
!     &               P(0) * EXP( -PBL_TOP_hPa(I,J) / SCALE_HEIGHT )
!            ENDIF

         ENDDO

         ! Error check
         IF ( ABS( SUM( F_OF_PBL(I,J,:) ) - 1.d0 ) > 1.d-3 ) THEN
            PRINT*, 'bad sum at: ', I, J
            CALL ERROR_STOP( 'Error in computing F_OF_PBL!',
     &                       'COMPUTE_PBL_HEIGHT ("pbl_mix_mod.f")' )
         ENDIF
      ENDDO
      ENDDO 
!$OMP END PARALLEL DO

      ! Model level where PBL top occurs
      PBL_MAX_L = MAXVAL( IMIX ) 

      END SUBROUTINE COMPUTE_PBL_HEIGHT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: turbday
!
! !DESCRIPTION: !  Subroutine TURBDAY executes the GEOS-Chem boundary layer 
!  mixing algorithm (full PBL mixing).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TURBDAY( NTRC, TC, TCVV, State_Met )              
!
! !USES:
!
      USE DIAG_MOD,           ONLY : TURBFLUP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CONV

      USE CMN_SIZE_MOD             ! Size parameters    
      USE CMN_DIAG_MOD             ! ND15
!
! !INPUT PARAMETERS: 
!
      ! Number of tracers used in computation
      INTEGER,        INTENT(IN)    :: NTRC             

      ! MW air (g/mol) / MW tracer (g/mol)     [ unitless ]
      REAL*8,         INTENT(IN)    :: TCVV(NTRC)

      ! Meteorology State object
      TYPE(MetState), INTENT(IN)    :: State_Met
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Tracer concentration [v/v]
      REAL*8,         INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NTRC) 
! 
! !REMARKS:
!  Original subroutine by Dale Allen, Univ of MD.
!
! !REVISION HISTORY: 
!  30 Jan 1998 - I. Bey, R. Yantosca - Initial version
!  (1 ) TURBDAY is written in Fixed-Form Fortran 90.  Also use F90
!        syntax for declarations (bmy, 4/1/99).
!  (2 ) New tracer concentrations are returned in TC.
!  (3 ) PS(I,J) is ACTUAL surface pressure and not Psurface - PTOP 
!  (4 ) Change in tracer in kg is now stored in DTC(I,J,L,N).  This makes
!        it easier to compute diagnostic quantities.  The new mixing ratio
!        is computed as TC(I,J,L,N) = TC(I,J,L,N) + DTC(I,J,L,N) / AD(I,J,L).
!  (5 ) XTRA2(*,*,5) is the height of the PBL in # of layers.  So if the
!        PBL top is located in the middle of the 3rd sigma layer at (I,J)
!        the value of XTRA2(I,J,5) would be 2.5.  The XTRA2 variable is
!        used by the HCTM drydep subroutines...it really is a historical
!        holdover.
!  (6 ) Restore the following NDxx diagnostics: (a) ND63 : Mass balance 
!        (CNVUPP) (b) ND15 : Mass change due to mixing in the boundary layer 
!  (7 ) Now pass TCVV and NCONV for the mass flux diagnostics.  Also
!        updated comments and cleaned up a few things. (bey, bmy, 11/10/99)
!  (8 ) Remove PTOP and XNUMOL from the arg list.  PTOP is now a parameter 
!        in "CMN_SIZE".  XNUMOL is no longer used in TURBDAY. (bmy, 2/10/00)
!  (9 ) Also removed obsolete ND63 diagnostics and updated comments.
!        (bmy, 4/12/00)
!  (10) Now use NTRC instead of NNPAR to dimension variables TC, TCVV, DTC, 
!        and DTCSUM (bmy, 10/17/00). 
!  (11) Removed obsolete code from 10/17/00 (bmy, 12/21/00)
!  (12) If the PBL depth is very small (or zero), then assume a PBL depth 
!        of 2 mb -- this prevents NaN's from propagating throughout the
!        code.  Also updated comments & made cosmetic changes. (bmy, 3/9/01)
!  (13) DTCSUM was declared twice but wasn't used.  Elminate declarations
!        to DTCSUM. (bmy, 7/16/01)
!  (14) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  Also updated comments. 
!        Also remove IREF, JREF and some debug output. (bmy, 9/25/01)
!  (15) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (16) Now takes in P=PS-PTOP instead of PS.  Redimension SIGE to 
!        (1:LLPAR+1).  
!  (17) Renamed PS to PZ so as not to conflict w/ the existing P variable.
!        Now pass P-PTOP thru PZ, in order to ensure that P and AD are
!        consistent w/ each other.  Added parallel DO-loops. Updated comments,
!        cosmetic changes.  Now print a header to stdout on the first call,
!        to confirm that TURBDAY has been called. (bmy, 4/11/02)
!  (18) Now use GET_PEDGE from "pressure_mod.f" to compute the pressure
!        at the bottom edge of grid box (I,J,L).  Deleted obsolete code from 
!        4/02.  Removed PZ, SIGE from the argument list, since we now compute
!        pressure from GET_PEDGE. (dsa, bdf, bmy, 8/22/02)  
!  (19)	Now reference AD, PBL from "dao_mod.f".  Now removed DXYP from the 
!        arg list, use GET_AREA_M2 from "grid_mod.f" instead.  Now removed 
!        NCONV, ALPHA_d, ALPHA_n from the arg list.  Now no longer reference 
!        SUNCOS.  Now set A(:,:)=1 day & nite; we assume full mixing all the 
!        time regardless of SUNCOS.  Updated comments, cosmetic changes.
!        (bmy, 2/11/03)
!  (20) Now can handle PBL field in meters for GEOS-4/fvDAS.  Also the
!        atmospheric scale height from CMN_GCTM. (bmy, 6/23/03)
!  (21) Now bundled into "pbl_mix_mod.f".  Broke off the part which computes
!        PBL height and related quantities into COMPUTE_PBL_HEIGHT. 
!        (bmy, 2/15/05)
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!   2 Mar 2012 - R. Yantosca - Remove reference to GET_AREA_M2
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I,  J,  L,     LTOP,    N
      REAL*8        :: AA, CC, CC_AA, AREA_M2, DTCONV
      REAL*8        :: A(IIPAR,JJPAR)
      REAL*8        :: DTC(IIPAR,JJPAR,LLPAR,NTRC)  

      ! Pointers
      REAL*8, POINTER :: AD(:,:,:)

      !=================================================================
      ! TURBDAY begins here!          
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'T U R B D A Y  -- by Dale Allen, U. Md.'
         WRITE( 6, '(a)' ) 'Modified for GEOS-CHEM by Bob Yantosca'
         WRITE( 6, '(a)' ) 'Last Modification Date: 2/4/03'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! Reset first time flag
         FIRST = .FALSE.
      ENDIF
      
      !=================================================================
      ! Do the boundary layer mixing
      !=================================================================

      ! Initialize pointers
      AD    => State_Met%AD

      ! Convection timestep [s]
      DTCONV = GET_TS_CONV() * 60d0

      ! Loop over Lat/Long grid boxes (I,J)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, AA, CC, CC_AA ) 
      DO J = 1, JJPAR
      DO I = 1, IIPAR
  
         ! We assume full mixing in the boundary layer, so the A
         ! coefficients are 1 everywhere, day & night (bmy, 2/11/03)
         A(I,J) = 1d0

         ! Calculate air mass within PBL at grid box (I,J,L)
         AA = 0.d0 
         DO L = 1, IMIX(I,J)-1 
            AA = AA + AD(I,J,L) 
         ENDDO       
               
         L  = IMIX(I,J) 
         AA = AA + AD(I,J,L) * FPBL(I,J)

         ! Loop over tracers 
         DO N = 1, NTRC
     
            !===========================================================
            ! Calculate tracer mass within PBL at grid box (I,J,L)
            !===========================================================
 
            ! Sum mass from (I,J,L) below the PBL top
            CC = 0.d0
            DO L = 1, IMIX(I,J)-1 
               CC = CC + AD(I,J,L) * TC(I,J,L,N)
            ENDDO       
               
            ! Then also sum mass from (I,J,L) which straddle the PBL top
            L     = IMIX(I,J) 
            CC    = CC + AD(I,J,L) * TC(I,J,L,N) * FPBL(I,J)
            
            ! CC/AA is the mean mixing ratio of tracer at 
            ! (I,J) from L=1 to L=LTOP
            CC_AA = CC / AA

            !========================================================
            ! TC(I,J,L,N) new  = TC(I,J,L,N) old + 
            !                    ( DTC(I,J,L,N) / AD(I,J,L) )
            !
            ! where
            !
            ! DTC(I,J,L,N) = [ alpha * (mean MR below PBL) * 
            !                  Airmass at (I,J,L) ] -
            !                [ alpha * TC(I,J,L,N) old     * 
            !                  Airmass at (I,J,L) ]
            !
            ! DTC is thus the change in mass (kg) due to BL mixing, 
            ! so DTC/AD is the change in (V/V) mixing ratio units.  
            !========================================================

            ! For grid boxes (I,J,L) which lie below the PBL top
            DO L = 1, IMIX(I,J)-1 
               DTC(I,J,L,N) = ( A(I,J) * CC_AA       * AD(I,J,L)  -
     &                          A(I,J) * TC(I,J,L,N) * AD(I,J,L) ) 

               TC(I,J,L,N) = TC(I,J,L,N) + DTC(I,J,L,N) / AD(I,J,L)
            ENDDO 

            ! For grid boxes (I,J,L) which straddle the PBL top
            L = IMIX(I,J) 

            DTC(I,J,L,N)  = 
     &           ( A(I,J) * FPBL(I,J)  * CC_AA       * AD(I,J,L) - 
     &             A(I,J) * FPBL(I,J)  * TC(I,J,L,N) * AD(I,J,L) ) 
               
            TC(I,J,L,N) = TC(I,J,L,N) + DTC(I,J,L,N) / AD(I,J,L)

            !=======================================================
            ! ND15 Diagnostic: 
            ! mass change due to mixing in the boundary layer
            !=======================================================
            IF ( ND15 > 0 ) THEN
               DO L = 1, IMIX(I,J)
                  TURBFLUP(I,J,L,N) = TURBFLUP(I,J,L,N) +
     &                 DTC(I,J,L,N) / ( TCVV(N) * DTCONV )
               ENDDO
            ENDIF
         ENDDO                  
      ENDDO                     
      ENDDO                     
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------------
!  Original code...leave here for reference (bmy, 11/10/99)
!                    TC(I,J,L,N) = 
!     &                ( A(I,J)     * AIRMAS(I,J,L) * CC/AA +
!     &                (1-A(I,J)) * TC(I,J,L,N)   * AIRMAS(I,J,L)) / 
!     &                AIRMAS(I,J,L) 
!
!                 TC(I,J,L,N) = 
!     &              ( A(I,J)        * FPBL(I,J)       * 
!     &                AIRMAS(I,J,L) * CC/AA           +
!     &               ( 1 - A(I,J)   * FPBL(I,J) )     *
!     &                TC(I,J,L,N)   * AIRMAS(I,J,L) ) / AIRMAS(I,J,L) 
!-----------------------------------------------------------------------------

      ! Free pointers
      NULLIFY( AD )

      END SUBROUTINE TURBDAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_frac_of_pbl
!
! !DESCRIPTION: Function GET\_FRAC\_OF\_PBL returns the fraction of grid box 
!  (I,J,L) that lies within the planetary boundary layer.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_FRAC_OF_PBL( I, J, L ) RESULT( FRAC )
!
! !USES:
!
      USE CMN_SIZE_MOD   ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I, J, L    ! Lon, lat, lev indices
!
! !RETURN VALUE:
!
      REAL*8              :: FRAC       ! Fraction of box (I,J,L) in the PBL
! 
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! GET_FRAC_OF_PBL begins here!
      !=================================================================
      IF ( L <= LLTROP ) THEN
         FRAC = F_OF_PBL(I,J,L)
      ELSE
         FRAC = 0d0
      ENDIF

      ! Return to calling program
      END FUNCTION GET_FRAC_OF_PBL
   !EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_frac_under_pbltop
!
! !DESCRIPTION: Function GET\_FRAC\_UNDER\_PBLTOP returns the fraction of 
!  grid box (I,J,L) that lies underneath the planetary boundary layer top.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_FRAC_UNDER_PBLTOP( I, J, L ) RESULT( FRAC )
!
! !USES:
!
      USE CMN_SIZE_MOD                 ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I, J, L   ! Lon, lat, level indices
!
! !RETURN VALUE:
!
      REAL*8              :: FRAC      ! Fraction of box (I,J,L) below PBL top
! 
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( L <= LLTROP ) THEN
         FRAC = F_UNDER_TOP(I,J,L)
      ELSE
         FRAC = 0d0
      ENDIF

      END FUNCTION GET_FRAC_UNDER_PBLTOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pbl_max_l
!
! !DESCRIPTION: Function GET\_PBL\_MAX\_L returns the model level at the 
!  highest part of the planetary boundary layer.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PBL_MAX_L() RESULT( TOP )
!
! !RETURN VALUE:
!
      INTEGER  :: TOP   ! Highest extent of PBL [model levels]
! 
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      TOP = PBL_MAX_L

      END FUNCTION GET_PBL_MAX_L
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pbl_top_hpa
!
! !DESCRIPTION: Function GET\_PBL\_TOP\_hPa returns the planetary boundary 
!  layer top [hPa] at a given GEOS-Chem surface location (I,J). 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PBL_TOP_hPa( I, J ) RESULT( TOP )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I, J   ! Lon and lat indices
!
! !RETURN VALUE:
!
      REAL*8              :: TOP    ! PBL top [hPa]
!
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      TOP = PBL_TOP_hPa(I,J)

      END FUNCTION GET_PBL_TOP_hPa
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pbl_top_l
!
! !DESCRIPTION: Function GET\_PBL\_TOP\_L returns the planetary boundary 
!  layer top [model levels] at a given GEOS-Chem surface location (I,J). 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PBL_TOP_L( I, J ) RESULT( TOP )
!
! !USES:
!

!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I, J   ! Lon and lat indices
!
! !RETURN VALUE:
!
      REAL*8              :: TOP    ! PBL top [model levels]
! 
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      TOP = PBL_TOP_L(I,J)

      END FUNCTION GET_PBL_TOP_L
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pbl_top_m
!
! !DESCRIPTION: Function GET\_PBL\_TOP\_m returns the planetary boundary 
!  layer top [m] at a given GEOS-CHEM surface location (I,J). 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PBL_TOP_m( I, J ) RESULT( TOP )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I, J   ! Lon and lat indices
!
! !RETURN VALUE:
!
      REAL*8              :: TOP    ! PBL top [m]
! 
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      TOP = PBL_TOP_m(I,J)

      END FUNCTION GET_PBL_TOP_m
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Function GET\_PBL\_THICK returns the thickness of the PBL
!  at a given surface location (I,J).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PBL_THICK( I, J ) RESULT( THICK )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I, J    ! Lon and lat indices
!
! !RETURN VALUE:
!
      REAL*8              :: THICK   ! PBL thickness [hPa]
! 
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THICK = PBL_THICK(I,J)

      END FUNCTION GET_PBL_THICK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pbl_mix
!
! !DESCRIPTION: Subroutine INIT\_PBL\_MIX allocates and zeroes module arrays
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_PBL_MIX
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE CMN_SIZE_MOD
! 
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_PBL_MIX begins here!
      !=================================================================

      ! Scalars
      PBL_MAX_L = 0

      ! Arrays
      ALLOCATE( IMIX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IMIX' )
      IMIX = 0

      ALLOCATE( FPBL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FPBL' )
      FPBL = 0d0

      ALLOCATE( F_OF_PBL( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'F_OF_PBL' )
      F_OF_PBL = 0d0

      ALLOCATE( F_UNDER_TOP( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'F_UNDER_TOP' )
      F_UNDER_TOP = 0d0

      ALLOCATE( PBL_TOP_hPa( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PBL_TOP_hPa' )
      PBL_TOP_hPa = 0d0

      ALLOCATE( PBL_TOP_L( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PBL_TOP_L' )
      PBL_TOP_L = 0d0

      ALLOCATE( PBL_TOP_m( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PBL_TOP_m' )
      PBL_TOP_m = 0d0

      ALLOCATE( PBL_THICK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PBL_THICK' )
      PBL_THICK = 0d0

      END SUBROUTINE INIT_PBL_MIX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_pbl_mix
!
! !DESCRIPTION: Subroutine CLEANUP\_PBL\_MIX allocates and zeroes 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_PBL_MIX
! 
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_PBL_MIX begins here!
      !=================================================================     
      IF ( ALLOCATED( IMIX        ) ) DEALLOCATE( IMIX        )
      IF ( ALLOCATED( FPBL        ) ) DEALLOCATE( FPBL        )
      IF ( ALLOCATED( F_OF_PBL    ) ) DEALLOCATE( F_OF_PBL    )
      IF ( ALLOCATED( F_UNDER_TOP ) ) DEALLOCATE( F_UNDER_TOP )
      IF ( ALLOCATED( PBL_TOP_hPa ) ) DEALLOCATE( PBL_TOP_hPa )
      IF ( ALLOCATED( PBL_TOP_L   ) ) DEALLOCATE( PBL_TOP_L   ) 
      IF ( ALLOCATED( PBL_TOP_m   ) ) DEALLOCATE( PBL_TOP_m   ) 
      IF ( ALLOCATED( PBL_THICK   ) ) DEALLOCATE( PBL_THICK   )

      END SUBROUTINE CLEANUP_PBL_MIX
!EOC
      END MODULE PBL_MIX_MOD
