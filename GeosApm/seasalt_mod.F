!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: seasalt_mod
!
! !DESCRIPTION: Module SEASALT\_MOD contains arrays and routines for performing
!  either a coupled chemistry/aerosol run or an offline seasalt aerosol
!  simulation. Original code taken from Mian Chin's GOCART model and modified
!  accordingly. (bec, rjp, bmy, 6/22/00, 11/23/09)
!\\
!\\
! !INTERFACE: 
!
      MODULE SEASALT_MOD
!
! !USES:
!

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CHEMSEASALT
      PUBLIC :: EMISSSEASALT
      PUBLIC :: CLEANUP_SEASALT
      PUBLIC :: GET_ALK
      PUBLIC :: INIT_SEASALT
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: SALT_V
      PUBLIC :: DMID
!
! !REMARKS:
!  Seasalt aerosol species: (1) Accumulation mode (usually 0.1 -  0.5 um)
!                           (2) Coarse mode       (usually 0.5 - 10.0 um)
!                                                                             .
!  NOTE: You can change the bin sizes for accumulation mode and coarse
!        mode seasalt in the "input.geos" file in v7-yy-zz and higher.
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
! !REVISION HISTORY:
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
!  (10) Save surface emissions separately (emis_save) for non-local scheme.
!        (ccc, 5/14/09)
!  (11) Bug fixes in GET_ALK and SRCSALT (bec, lyj, bmy, 11/23/09)
!  (12) Add size-resolved emission subroutine SRCSALT30 and reference to 
!        tomas_mod.f. (win, 7/17/09)
!  16 Feb 2011 - R. Yantosca - Extend seasalt bins for APM simulation (G. Luo)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  09 Feb 2012 - R. Yantosca - Moved APM routines to end of module, which 
!                              facilitates comparison with GeosCore code
!  16 Feb 2012 - R. Yantosca - Moved SRCSALT30 to end of module
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  04 Mar 2013 - R. Yantosca - Now call INIT_SULFATE from the init stage
!                              which facilitates connection to GEOS-5 GCM
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER   :: NSALT = 2    ! Number of sea salt tracers
      INTEGER, PARAMETER   :: NR_MAX = 200

#if   defined( APM )
      INTEGER, PARAMETER   :: NSALTBIN = 20
#endif
!
! !PRIVATE TYPES:
!
      !========================================================================
      ! Module Variables:
      !
      ! DRYSALA  : Drydep index for accumulation mode sea salt
      ! DRYSALC  : Drydep index for coarse mode sea salt
      ! IDDEP    : Drydep index array for sea salt tracers
      ! REDGE    : Array for edges of seasalt radius bins
      ! RMID     : Array for centers of seasalt radius bins
      ! SRC      : Array for baseline seasalt emission/bin [kg/m2]
      ! SRC_N    : Array for baseline seasalt emission/bin [#/m2]
      ! SS_DEN   : Sea salt density [kg/m3]
      ! ALK_EMIS : Array for alkalinity [kg]
      ! N_DENS   : Number density of seasalt emissions [#/m3]
      ! SALT_V   : Log-normal volum size distribution for sea salt
      !=======================================================================-

      ! Scalars
      INTEGER              :: DRYSALA
      INTEGER              :: DRYSALC
#if   defined( APM )
      INTEGER, SAVE        :: IDTSEA1, IDTSEA2
#endif

      ! Arrays
      INTEGER              :: IDDEP(NSALT)
      REAL*8,  ALLOCATABLE :: REDGE(:,:)   
      REAL*8,  ALLOCATABLE :: RMID(:,:)
      REAL*8,  ALLOCATABLE :: SRC(:,:)   
      REAL*8,  ALLOCATABLE :: SRC_N(:,:)   
      REAL*8,  ALLOCATABLE :: ALK_EMIS(:,:,:,:)
      REAL*8,  ALLOCATABLE :: N_DENS(:,:,:,:)    
      REAL*8,  ALLOCATABLE :: SALT_V(:)
      REAL*8,  ALLOCATABLE :: DMID(:)
      REAL*8               :: SS_DEN(NSALT)    = (/ 2200.d0, 2200.d0 /)

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
! !IROUTINE: chemseasalt
!
! !DESCRIPTION: Subroutine CHEMSEASALT is the interface between the GEOS-CHEM
!  main program and the seasalt chemistry routines that mostly calculates
!  seasalt dry deposition (rjp, bmy, 1/24/02, 5/23/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMSEASALT( am_I_Root, Input_Opt, 
     &                        State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DRYDEP_MOD,         ONLY : DEPNAME
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE DRYDEP_MOD,         ONLY : NUMDEP
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTSALA
      USE TRACERID_MOD,       ONLY : IDTSALC
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!      
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY:
!  (1 ) Now reference STT from "tracer_mod.f".  Now references LPRT from
!        "logical_mod.f" (bmy, 7/20/04)
!  (2 ) Now only call DRY_DEPOSITION if LDRYD=T (bec, bmy, 5/23/06)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  13 Nov 2012 - R. Yantosca - Now add Input_Opt, RC arguments for GIGC
!  15 Nov 2012 - M. Payer    - Now pass met fields via State_Met object
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Non-SAVEd variables
      LOGICAL        :: prtDebug
      LOGICAL        :: LDRYD
      LOGICAL        :: LNLPBL
      LOGICAL        :: LPRT
      INTEGER        :: N

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEMSEASALT begins here!
      !=================================================================

      ! Assume success
      RC       = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LDRYD    = Input_Opt%LDRYD
      LNLPBL   = Input_Opt%LNLPBL
      LPRT     = Input_Opt%LPRT

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT      => State_Chm%Tracers

      ! Do we have to print debug output?
      prtDebug = ( LPRT .and. am_I_Root )

      !=================================================================
      ! Maybe someday we should merge these two separate calculations
      ! into one (rjp, 4/3/04)
      !=================================================================

      !-------------------
      ! Accumulation mode
      !-------------------
      CALL WET_SETTLING( am_I_Root,          Input_Opt, State_Met,
     &                   STT(:,:,:,IDTSALA), 1,         RC         )


      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSEASALT: WET_SET, Accum' )
      ENDIF

      ! If LNLPBL (non local PBL mixing) is turned on, do sea salt
      ! dry deposition in vdiff as for all the other aerosols (jaegle 5/11/11)
      IF ( LDRYD .and. ( .not. LNLPBL ) ) THEN
         CALL DRY_DEPOSITION( am_I_Root,          Input_Opt, 
     &                        STT(:,:,:,IDTSALA), 1,         RC )

         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMSEASALT: DRY_DEP, Accum' )
         ENDIF
      ENDIF

      !-------------------
      ! Coarse mode
      !-------------------
      CALL WET_SETTLING( am_I_Root,          Input_Opt, State_Met,
     &                   STT(:,:,:,IDTSALC), 2,         RC         )

      IF ( prtDebug ) THEN 
         CALL DEBUG_MSG( '### CHEMSEASALT: WET_SET, Coarse' )
      ENDIF

      ! If LNLPBL (non local PBL mixing) is turned on, do sea salt
      ! dry deposition in vdiff as for all the other aerosols (jaegle 5/11/11)
      IF ( LDRYD .and. ( .not. LNLPBL ) ) THEN
         CALL DRY_DEPOSITION( am_I_Root,          Input_Opt,
     &                        STT(:,:,:,IDTSALC), 2,         RC )

         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMSEASALT: DRY_DEP, Coarse')
         ENDIF
      ENDIF

#if   defined( APM )
      !-------------------
      ! APM microphysics
      !-------------------
      CALL WET_SETTLINGBIN( STT(:,:,:,IDTSEA1:IDTSEA2),
     &                      State_Met, State_Chm )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSEASALT: WET_SET, Bins' )
      ENDIF

      IF ( LDRYD ) THEN
         CALL DRY_DEPOSITIONBIN( STT(:,:,:,IDTSEA1:IDTSEA2),
     &                           State_Chm )
         IF (  prtDebug ) THEN
            CALL DEBUG_MSG('### CHEMSEASALT: DRY_DEP, Bins')
         ENDIF
      ENDIF
#endif

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEMSEASALT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: wet_settling
!
! !DESCRIPTION: Subroutine WET\_SETTLING performs wet settling of sea salt.
!  (bec, rjp, bmy, 4/20/04, 6/11/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WET_SETTLING( am_I_Root, Input_Opt, State_Met, 
     &                         TC,        N,         RC         )
!
! !USES:
!
      USE CMN_GCTM_MOD
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TRACERID_MOD,       ONLY : IDTSALA
      USE TRACERID_MOD,       ONLY : IDTSALC
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE GRID_MOD,           ONLY : GET_AREA_CM2
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: N                     ! 1=accum mode; 
                                                             ! 2=coarse mode
      LOGICAL,        INTENT(IN)    :: am_I_Root             ! Root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
      TYPE(MetState), INTENT(IN)    :: State_Met             ! MeteorologyState
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,         INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! Sea salt [kg]
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC                    ! Success/failure 
!
! !REVISION HISTORY:
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (3 ) Bug fix: DTCHEM has to be REAL*8, not integer. (bmy, 9/7/06)
!  (4 ) Now limit relative humidity to [tiny(real*8),0.99] range for DLOG
!         argument (phs, 5/1/08)
!  (5 ) Update sea salt density calculation using Tang et al. (1997) (bec,
!       jaegle 5/11/11)
!  (6 ) Update hygroscopic growth for sea salt using Lewis and Schwartz (2006)
!       and and density calculation based on Tang et al. (1997) (bec, jaegle
!       5/11/11)
!  (7 ) Itegrate settling velocity over entire size distribution (jaegle
!        5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%XNUMOL(1:N_TRACERS)
!  12 Jun 2013 - R. Yantosca - Bug fix: SALT_MASS needs to be !OMP PRIVATE
!  12 Jun 2013 - R. Yantosca - Reformatted some comments for clarity
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I,      J,     L,        N_TRACERS
      REAL*8                 :: DELZ,   DELZ1, REFF,     DEN
      REAL*8                 :: P,      DP,    PDP,      TEMP        
      REAL*8                 :: CONST,  SLIP,  VISC,     FAC1
      REAL*8                 :: FAC2,   FLUX,  AREA_CM2, RHB
      ! replace RCM with RUM (radis in micron) jaegle 5/11/11
      REAL*8                 :: RUM,    RWET,  RATIO_R,  RHO
      REAL*8                 :: TOT1,   TOT2,  DTCHEM
      REAL*8                 :: VTS(LLPAR)  
      REAL*8                 :: TC0(LLPAR)
      ! New variables (jaegle 5/11/11)
      REAL*8                 :: SW
      REAL*8                 :: R0,       R1, NR, DEDGE, SALT_MASS
      REAL*8                 :: SALT_MASS_TOTAL, VTS_WEIGHT, DMIDW
      REAL*8                 :: WTP, RHO1 
      INTEGER                :: ID
      LOGICAL, SAVE          :: FIRST = .TRUE.
! !DEFINED PARAMETERS:
!      
      REAL*8,  PARAMETER     :: C1 =  0.7674d0 
      REAL*8,  PARAMETER     :: C2 =  3.079d0 
      REAL*8,  PARAMETER     :: C3 =  2.573d-11
      REAL*8,  PARAMETER     :: C4 = -1.424d0
      ! Parameters for polynomial coefficients to derive seawater
      ! density. From Tang et al. (1997) (jaegle 5/11/11)
      REAL*8,  PARAMETER     :: A1 =  7.93d-3
      REAL*8,  PARAMETER     :: A2 = -4.28d-5
      REAL*8,  PARAMETER     :: A3 =  2.52d-6
      REAL*8,  PARAMETER     :: A4 = -2.35d-8
      ! increment of radius for integration of settling velocity (um)
      REAL*8, PARAMETER      :: DR    = 5.d-2
      ! parameter for convergence
      REAL*8,  PARAMETER     :: EPSI = 1.0D-4
      ! parameters for assumed size distribution of acc and coarse mode
      ! sea salt aerosols (jaegle 5/11/11)
      ! geometric dry mean diameters (microns)
      REAL*8,  PARAMETER     ::   RG_A = 0.085d0
      REAL*8,  PARAMETER     ::   RG_C = 0.4d0
      ! sigma of the size distribution
      REAL*8,  PARAMETER     ::   SIG_A = 1.5d0
      REAL*8,  PARAMETER     ::   SIG_C = 1.8d0

      ! Local variables for Input_Opt quantities
      REAL*8                 :: SALA_REDGE_um(2)
      REAL*8                 :: SALC_REDGE_um(2)
      REAL*8                 :: XNUMOL(Input_Opt%N_TRACERS)

      !=================================================================
      ! WET_SETTLING begins here!
      !=================================================================


      ! Assume success
      RC            = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      SALA_REDGE_um = Input_Opt%SALA_REDGE_um 
      SALC_REDGE_um = Input_Opt%SALC_REDGE_um
      N_TRACERS     = Input_Opt%N_TRACERS
      XNUMOL        = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Chemistry timestep [s]
      DTCHEM        = GET_TS_CHEM() * 60d0

      ! Sea salt density [kg/m3]
      DEN            = SS_DEN( N )

      ! Seasalt effective radius (i.e. midpt of radius bin) [m]
      SELECT CASE ( N )

         ! Accum mode
         ! add R0 and R1 = edges if the sea salt size bins (jaegle 5/11/11)
         CASE( 1 )
            REFF = 0.5d-6 * ( SALA_REDGE_um(1) + SALA_REDGE_um(2) )
            R0 = SALA_REDGE_um(1)
            R1 = SALA_REDGE_um(2)

         ! Coarse mode
         CASE( 2 ) 
            REFF = 0.5d-6 * ( SALC_REDGE_um(1) + SALC_REDGE_um(2) )
            R0 = SALC_REDGE_um(1)
            R1 = SALC_REDGE_um(2)
            
      END SELECT

      ! Number of dry radius size bins between lowest radius (accumulation
      ! mode) and largest radii (coarse mode) (jaegle 5/11/11)
      
      NR = INT( ( ( SALC_REDGE_um(2) - SALA_REDGE_um(1) ) / DR ) 
     &                  + 0.5d0 )

      ! Error check
      IF ( NR > NR_MAX ) THEN
        CALL ERROR_STOP( 'Too many bins!', 'SRCSALT (seasalt_mod.f)')
      ENDIF

      !=================================================================
      ! Define the volume size distribution of sea-salt. This only has
      ! to be done once. We assume that sea-salt is the combination of a 
      ! coarse mode and accumulation model log-normal distribution 
      ! functions (jaegle 5/11/11)
      !=================================================================
      IF ( FIRST) THEN

        ! Lower edge of 0th bin
	DEDGE=SALA_REDGE_um(1) * 2d0

	! Loop over diameters
        DO ID = 1, NR
           ! Diameter of mid-point in microns
           DMID(ID)  = DEDGE + ( DR )

	   ! Calculate the dry volume size distribution as the sum of two 
           ! log-normal size distributions. The parameters for the size 
           ! distribution are based on Reid et al. and Quinn et al.
	   ! The scaling factors 13. and 0.8 for acc and coarse mode 
           ! aerosols are chosen to obtain a realistic distribution  
	   ! SALT_V (D) = dV/dln(D) [um3]
	   SALT_V(ID) = PI / 6d0* (DMID(ID)**3) * (
     &         13d0*exp(-0.5*( LOG(DMID(ID))-LOG(RG_A*2d0) )**2d0/
     &                   LOG(SIG_A)**2d0 )
     &         /( sqrt(2d0 * PI) * LOG(SIG_A) )  +
     &         0.8d0*exp(-0.5*( LOG(DMID(ID))-LOG(RG_C*2d0) )**2d0/
     &                   LOG(SIG_C)**2d0)
     &         /( sqrt(2d0 * PI) * LOG(SIG_C) )  )
	   ! update the next edge
	   DEDGE = DEDGE + DR*2d0
        ENDDO

        ! Reset after the first time
        IF ( FIRST ) FIRST = .FALSE.
      ENDIF


      ! Sea salt radius [cm]
      !RCM  = REFF * 100d0  
      ! The radius used in the Gerber formulation for hygroscopic growth
      ! of sea salt should be in microns (RUM) instead of cm (RCM). Replace RCM
      ! with RUM (jaegle 5/11/11)
      !RUM  = REFF * 1d6  

      ! Exponential factors
      !FAC1 = C1 * ( RCM**C2 )
      !FAC2 = C3 * ( RCM**C4 )
      ! Replace with RUM (jaegle 5/11/11)
      !FAC1 = C1 * ( RUM**C2 )
      !FAC2 = C3 * ( RUM**C4 )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,        J,     L,     VTS,             P          )     
!$OMP+PRIVATE( TEMP,     RHB,   RWET,  RATIO_R,         RHO        )
!$OMP+PRIVATE( DP,       PDP,   CONST, SLIP,            VISC       )
!$OMP+PRIVATE( TC0,      DELZ,  DELZ1, TOT1,            TOT2       )
!$OMP+PRIVATE( AREA_CM2, FLUX,  ID,    SALT_MASS_TOTAL, VTS_WEIGHT )
!$OMP+PRIVATE( DMIDW,    RHO1,  WTP,   SALT_MASS                   )
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
            TEMP    = State_Met%T(I,J,L)

            ! Cap RH at 0.99 
            RHB     = MIN( 0.99d0, State_Met%RH(I,J,L) * 1d-2 )

            ! Safety check (phs, 5/1/08)
            RHB     = MAX( TINY(RHB), RHB           )

            ! Aerosol growth with relative humidity in radius [m] 
            ! (Gerber, 1985)
            !RWET    = 0.01d0*(FAC1/(FAC2-DLOG(RHB))+RCM**3.d0)**0.33d0
	    ! Fix bugs in the Gerber formula:  a log10 (instead of ln) 
            ! should be used and the dry radius should be expressed in 
            ! micrometers (instead of cm) also add more significant 
	    ! digits to the exponent (should be 1/3) (jaegle 5/11/11)
            !RWET    = 1d-6*(FAC1/(FAC2-LOG10(RHB))+RUM**3.d0)**0.33333d0

            ! Use equation 5 in Lewis and Schwartz (2006) for sea 
            ! salt growth (bec, jaegle 5/11/11)
            RWET = REFF * (4.d0 / 3.7d0) *
     &              ( (2.d0 - RHB)/(1.d0 - RHB) )**(1.d0/3.d0)


            ! Ratio dry over wet radii at the cubic power
            RATIO_R = ( REFF / RWET )**3.d0

            ! Density of the wet aerosol (kg/m3)
            RHO     = RATIO_R * DEN + ( 1.d0 - RATIO_R ) * 1000.d0

            ! Above density calculation is chemically unsound because 
            ! it ignores chemical solvation.   Iteratively solve Tang et al., 
            ! 1997 equation 5 to calculate density of wet aerosol (kg/m3) 
            ! (bec, jaegle 5/11/11)
            RATIO_R = ( REFF / RWET )
            ! Assume an initial density of 1000 kg/m3
            RHO  = 1000.D0
            RHO1 = 0.d0 !initialize (bec, 6/21/10)
            DO WHILE ( ABS( RHO1-RHO ) .gt. EPSI )
                ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet) 
                WTP    = 100.d0 * DEN/RHO * RATIO_R**3.d0
                ! Then calculate density of wet aerosol using equation 5 
                ! in Tang et al., 1997 [kg/m3]
                RHO1   = ( 0.9971d0 + (A1 * WTP) + (A2 * WTP**2.d0) + 
     $               (A3 * WTP**3.d0) + (A4 * WTP**4.d0) ) * 1000.d0
                ! Now calculate new weight percent using above density 
                ! calculation
                WTP    = 100.d0 * DEN/RHO1 * RATIO_R**3.d0
                ! Now recalculate new wet density [kg/m3]
                RHO   = ( 0.9971d0 + (A1 * WTP) + (A2 * WTP**2.d0) + 
     $              (A3 * WTP**3.d0) + (A4 * WTP**4.d0) ) * 1000.d0
            ENDDO

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

            ! This settling velocity is for the mid-point of the size bin. 
            ! In the following we derive scaling factors to take into account
	    ! the strong dependence on radius of the settling velocity and the
	    ! mass size distribution:
	    !  VTS_WEIGHTED = total( M(k) x VTS(k)) / total( M(k) ) 
	    ! The settling velocity is a function of the radius squared 
            ! (see definition of CONST above) 
	    ! so VTS(k) = VTS * (RMID(k)/RWET)^2
            ! (jaegle 5/11/11)

	    SALT_MASS_TOTAL = 0d0
	    VTS_WEIGHT      = 0d0
	    DO ID = 1, NR
	       ! Calculate mass of wet aerosol (Dw = wet diameter, D = 
               ! dry diameter): dM/dlnDw = dV/dlnDw * RHO, we assume that 
               ! the density of sea-salt doesn't change much over the size 
               ! range.  and  dV/dlnDw = dV/dlnD * dlnD/dlnDw = 
               ! dV/dlnD * Dw/D = dV/dlnD * Rwet/Rdry
	       ! Further convert to dM/dDw = dM/dln(Dw) * dln(Dw)/Dw = 
               ! dM/dln(Dw)/Dw
	       ! Overall = dM/dDw = dV/dlnD * Rwet/Rdry * RHO /Rw
	       ! 
	       IF (DMID(ID) .ge. R0*2d0 .and. DMID(ID) .le. R1*2d0 ) THEN
	         DMIDW = DMID(ID) * RWET/REFF  ! wet radius [um]
	         SALT_MASS   = SALT_V(ID) * RWET/REFF * RHO / (DMIDW*0.5d0)
	         VTS_WEIGHT  = VTS_WEIGHT + 
     &              SALT_MASS * VTS(L) * (DMIDW/(RWET*1d6*2d0) )**2d0 *
     &                            (2d0 * DR *  RWET/REFF)
	         SALT_MASS_TOTAL=SALT_MASS_TOTAL+SALT_MASS *
     &                            (2d0 * DR *  RWET/REFF)
	       ENDIF

            ENDDO
            ! Calculate the weighted settling velocity:
            VTS(L) = VTS_WEIGHT/SALT_MASS_TOTAL
         ENDDO

         ! Method is to solve bidiagonal matrix which is
         ! implicit and first order accurate in z (rjp, 1/24/02)

         ! Save initial tracer concentration in column
         DO L = 1, LLPAR
            TC0(L) = TC(I,J,L)
         ENDDO

         ! We know the boundary condition at the model top
         L    = LLTROP
         DELZ = State_Met%BXHEIGHT(I,J,L)

         TC(I,J,L) = TC(I,J,L) / ( 1.d0 + DTCHEM * VTS(L) / DELZ )

         DO L = LLTROP-1, 1, -1
            DELZ  = State_Met%BXHEIGHT(I,J,L)
            DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
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
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! Convert sea salt flux from [kg/s] to [molec/cm2/s]
            FLUX     = ( TOT1 - TOT2 ) / DTCHEM
            FLUX     = FLUX * XNUMOL(IDTSALA) / AREA_CM2 
   
            ! Store in AD44 array
            AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE WET_SETTLING
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dry_deposition
!
! !DESCRIPTION: Subroutine DRY\_DEPOSITION computes the loss of sea salt by dry
!  deposition at the surface, using an implicit method (bec, rjp, bmy, 4/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DRY_DEPOSITION( am_I_Root, Input_Opt, TC, N, RC )
!
! !USES:
!
      USE CMN_GCTM_MOD
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACERID_MOD,       ONLY : IDTSALA
      USE TRACERID_MOD,       ONLY : IDTSALC
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE GRID_MOD,           ONLY : GET_AREA_CM2
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root             ! Is this root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
      INTEGER,        INTENT(IN)    :: N                     ! 1=accum; 2=coarse
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,         INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! Sea salt [kg]
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
! !REVISION HISTORY:
!  (1 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (2 ) Update to calculate the drydep throughout the entire PBL instead of
!       just at the surface. This is more in line with what is done in
!       dry_dep.f. This is only used if LNLPBL is turned off (or for GEOS-4 and
!       prior met fields). (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  20 Nov 2012 - R. Yantosca - Bug fix: use temp array to archive ND44 drydep
!                              then save back to AD44 outside parallel loop
!  20 Nov 2012 - R. Yantosca - Now loop from 1..LLPAR for GIGC.  Remove
!                              references to GET_PBL_MAX_L.
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%XNUMOL(1:N_TRACERS)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N_TRACERS
      INTEGER            :: I,        J,     L,      DTCHEM
      REAL*8             :: OLD,      NEW,   G,      REFF
      REAL*8             :: DIAM,     U_TS0, REYNOL, ALPHA 
      REAL*8             :: BETA,     GAMMA, DENS,   FLUX 
      REAL*8             :: AREA_CM2, TOT1,  TOT2
      REAL*8             :: F_UNDER_TOP, FREQ

      ! Arrays
      REAL*8             :: ND44_TMP(IIPAR,JJPAR,LLPAR)
      REAL*8             :: XNUMOL(Input_Opt%N_TRACERS)
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER :: RHOA = 1.25d-3

      !=================================================================
      ! DRY_DEPOSITION begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Chemistry timestep [s]
      DTCHEM    = GET_TS_CHEM() * 60d0

      ! Zero temporary array for drydep diagnostic
      IF ( ND44 > 0 ) ND44_TMP = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, OLD, NEW, FLUX )
!$OMP+PRIVATE( L, F_UNDER_TOP   , FREQ ) ! (jaegle 5/11/11)
!$OMP+SCHEDULE( DYNAMIC )

      ! Note, for grid-independent GEOS-Chem we have to loop over
      ! all levels and then just shunt if we are in the stratosphere
      ! (bmy, 11/20/12)
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( I, J, L )

         ! Fraction of box (I,J,L) under PBL top [unitless]
	 F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Only apply drydep to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0d0 ) THEN

            ! Sea salt dry deposition frequency [1/s] accounting
            ! for fraction of each grid box located beneath the PBL top
            FREQ = DEPSAV(I,J,IDDEP(N)) * F_UNDER_TOP

            ! Only apply drydep loss if FREQ is nonzero
            IF ( FREQ > 0d0 ) THEN

               ! Old tracer concentration [kg]
               OLD  = TC(I,J,L)

               ! New tracer concentration [kg]
               NEW  = OLD * EXP( -FREQ * DTCHEM )

               ! Old tracer concentration [kg]
               !OLD  = TC(I,J,1)

               ! New tracer concentration [kg]
               !NEW  = OLD * EXP( -DEPSAV(I,J,IDDEP(N)) * DTCHEM  )

               !===========================================================
               ! ND44 diagnostic: sea salt drydep loss [molec/cm2/s]
               !===========================================================
               IF ( ND44 > 0 ) THEN
            
                  ! Convert drydep loss from [kg/s] to [molec/cm2/s]
                  FLUX = ( OLD - NEW ) / DTCHEM 
                  FLUX = FLUX * XNUMOL(IDTSALA) / AREA_CM2 
         
                  ! Save in ND44 temporary array
                  ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
               ENDIF

               ! Update tracer array
               TC(I,J,L) = NEW 
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !====================================================================
      ! ND44 diagnostic: save into AD44 array, summing in the vertical
      !====================================================================
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD44(I,J,IDDEP(N),1) = SUM( ND44_TMP(I,J,:) ) 
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      END SUBROUTINE DRY_DEPOSITION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissseasalt
!
! !DESCRIPTION: Subroutine EMISSSEASALT is the interface between the GEOS-Chem
!  model and the SEASALT emissions routines in "seasalt\_mod.f".
!  (bec, rjp, bmy, 3/24/03, 2/22/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISSSEASALT( am_I_Root, Input_Opt, State_Met, 
     &                         State_Chm, RC,        SSA_Br2    )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTSALA, IDTSALC
      USE VDIFF_PRE_MOD,      ONLY : emis_save ! (Lin, 03/31/09)
#if   defined( TOMAS )
      USE TRACERID_MOD,       ONLY : IDTNK1 ! (win, 7/17/09)
      USE TRACERID_MOD,       ONLY : IDTSS1 ! (win, 7/17/09)
      USE TOMAS_MOD,          ONLY : IBINS  ! (win, 7/17/09)
#endif
#if   defined( APM )
      USE APM_INIT_MOD,       ONLY : IDTSEABIN1
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,           INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput),    INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState),    INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState),    INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,           INTENT(OUT)   :: RC          ! Success or failure?
      REAL*8,  OPTIONAL, INTENT(OUT)   :: SSA_Br2(IIPAR, JJPAR)
!
! !REVISION HISTORY: 
!  (1 ) Now references LPRT from "logical_mod.f" and STT from "tracer_mod.f".
!        (bmy, 7/20/04)
!  (2 ) Now make sure IDTSALA, IDTSALC are nonzero before calling SRCSALT.
!        (bmy, 1/26/05)
!  (3 ) Remove reference to header file "CMN" (bmy, 2/22/05)
!  (4 ) Now call INIT_SEASALT on the first timestep.  Also initialize ALK_EMIS
!        and N_DENS on each timestep. (bec, bmy, 4/13/05)
!  (5 ) Call SRCSALT30 for size-resolved sea-salt emission (win, 7/17/09)
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE   :: FIRST = .TRUE.
      INTEGER         :: I, J, L, N

      ! jpp, 3/2/10
      REAL*8          :: tmpbr2(IIPAR, JJPAR)

      ! For fields from Input_Opt
      LOGICAL         :: LNLPBL
      LOGICAL         :: LPRT

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! EMISSSEASALT begins here! 
      !=================================================================

      ! Assume success
      RC          =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LNLPBL      = Input_Opt%LNLPBL
      LPRT        = Input_Opt%LPRT

      ! Initialize Br2 emission arrays
      tmpbr2(:,:) = 0.d0
      IF ( PRESENT ( SSA_Br2 ) ) SSA_Br2(:,:) = 0.d0

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### in EMISSEASALT' )

#if   defined( APM )
         ! Define extra APM seasalt bins
         IDTSEA1 = IDTSEABIN1
         IDTSEA2 = IDTSEABIN1 + NSALTBIN - 1
#endif

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
         CALL SRCSALT( STT(:,:,:,IDTSALA), tmpbr2(:,:), 1, State_Met )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSEASALT: Accum' )
      ENDIF

      IF ( PRESENT ( SSA_Br2 ) ) THEN
         ! 1. store contributions of Br2 emissions
         !    from accumulation mode: (jpp, 3/2/10)
         SSA_Br2(:,:) = tmpbr2(:,:)

         ! reset the temporary variable
         tmpbr2(:,:)  = 0.d0
      ENDIF

      ! Coarse mode emissions
      IF ( IDTSALC > 0 ) THEN
         CALL SRCSALT( STT(:,:,:,IDTSALC), tmpbr2(:,:), 2, State_Met )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSEASALT: Coarse' )
      ENDIF

#if defined( DEVEL )
      State_Chm%TRAC_TEND(:,:,1,IDTSALA) = 
     &State_Chm%TRAC_TEND(:,:,1,IDTSALA) + ALK_EMIS(:,:,1,1)
      State_Chm%TRAC_TEND(:,:,1,IDTSALC) = 
     &State_Chm%TRAC_TEND(:,:,1,IDTSALC) + ALK_EMIS(:,:,1,2)
#endif

#if   defined( APM )
      ! Sea salt emissions for extra APM bins
      CALL SRCSALTBIN( STT(:,:,:,IDTSEA1:IDTSEA2), State_Met )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSEASALT: Bin' )
#endif

      IF ( PRESENT ( SSA_Br2 ) ) THEN
         ! 2. add in contributions of Br2 emissions
         !    from the coarse mode: (jpp, 3/2/10)
         SSA_Br2(:,:) = SSA_Br2(:,:) + tmpbr2(:,:)
      ENDIF

#if   defined( TOMAS )
      ! 30-bin Sea-salt emissions
      IF ( IDTSS1 > 0 ) THEN

         ! For now, stop the run if choose to use Lin's PBL scheme 
         ! when running TOMAS aerosols (win, 7/17/09)
         IF( LNLPBL ) THEN
            print *,'Currently the TOMAS size-resolved sea salt ',
     &           'emission code does not support the option to do',
     &           'non-local PBL scheme (convection menu)!  Try not ',
     &           'using the non-local PBL or modify SRCSALT30!'
            stop
         ENDIF
         
         CALL SRCSALT30( STT(:,:,:,IDTNK1:IDTNK1+IBINS-1), 
     &                   STT(:,:,:,IDTSS1:IDTSS1+IBINS-1 ),
     &                   State_Met )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSSEASALT: TOMAS' )
      ENDIF
#endif

      ! save surface emis for non-local PBL mixing (vdiff_mod.f) (units: kg)
      ! (Lin, 06/09/08) 
      IF ( LNLPBL ) THEN
!$OMP PARALLEL DO
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            emis_save(I,J,IDTSALA) = ALK_EMIS(I,J,1,1)
            emis_save(I,J,IDTSALC) = ALK_EMIS(I,J,1,2)
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE EMISSSEASALT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcsalt
!
! !DESCRIPTION: Subroutine SRCSALT is based on the sea salt source function of
!  Gong (2003) with the empirical sea surface temperature (SST) dependence of
!  Jaegle et al. (2011). This SST dependence was derived based on comparisons
!  to cruise observations of coarse mode sea salt mass concentrations.
! 
!  Contact: Lyatt Jaegle (jaegle@uw.edu)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SRCSALT( TC, SSA_Br2, N, State_Met )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_WATER
      USE DIAG_MOD,           ONLY : AD08
      USE ERROR_MOD,          ONLY : DEBUG_MSG, ERROR_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE SSA_BROMINE_MOD,    ONLY : EMISS_SSA_BROMINE
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TRACER_MOD,         ONLY : SALA_REDGE_um, SALC_REDGE_um
      USE TRACER_MOD,         ONLY : XNUMOL

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND44, ND08
      USE CMN_GCTM_MOD             ! PI

      !%%% NOTE: Keep this here as a placeholder, but we should
      !%%% eventually replace this with Input_Opt%LNLPBL
      USE LOGICAL_MOD,         ONLY : LNLPBL
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: N  ! N=1 is accum mode; N=2 is coarse mode
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,       INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! Sea salt [v/v]
!
! !OUTPUT PARAMETERS:
!
      !jpp, 3/2/10
      REAL*8,       INTENT(OUT)   :: SSA_Br2(IIPAR,JJPAR)  ! Br2 sea salt emissions
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
!  (2 ) Gong, S., L. Barrie, and J.-P. Blanchet, "Modeling sea-salt
!        aerosols in the atmosphere. 1. Model development", J. Geophys. Res.,
!        v. 102, 3805-3818, 1997.
!  (3 ) Gong, S. L., "A parameterization of sea-salt aerosol source function
!        for sub- and super-micron particles", Global Biogeochem.  Cy., 17(4), 
!        1097, doi:10.1029/2003GB002079, 2003.
!  (4 ) Jaegle, L., P.K. Quinn, T.S. Bates, B. Alexander, J.-T. Lin, "Global
!        distribution of sea salt aerosols: New constraints from in situ and 
!        remote sensing observations", Atmos. Chem. Phys., 11, 3137-3157, 
!        doi:10.5194/acp-11-3137-2011.
!
! !REVISION HISTORY: 
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_TOP_L from "pbl_mix_mod.f".
!        Removed reference to header file CMN.  Removed reference to 
!        "pressure_mod.f".  (bmy, 2/22/05)
!  (3 ) Now also compute alkalinity and number density of seasalt emissions.
!        (bec, bmy, 4/13/05)
!  (4 ) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) The source function is for wet aerosol radius (RH=80%, with a radius
!        twice the size of dry aerosols) so BETHA should be set to 2 
!        instead of 1.  Also now use LOG10 instead of LOG in the expressions
!        for the seasalt base source, since we need the logarithm to the base
!        10. (jaegle, bec, bmy, 11/23/09)
!  (6 ) Update to use the Gong (2003) source function (jaegle 5/11/11)
!  (7 ) Apply an empirical sea surface temperature dependence to Gong (2003)
!       (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I,     J,      L
      INTEGER                :: R,     NR,     NTOP
      REAL*8                 :: W10M,  DTEMIS, R0
      REAL*8                 :: R1,    CONST,  CONST_N
      REAL*8                 :: FEMIS, A_M2
      REAL*8                 :: SALT(IIPAR,JJPAR)
      REAL*8                 :: SALT_N(IIPAR,JJPAR)
      ! New variables (jaegle 5/11/11)
      REAL*8                 :: A, B, SST, SCALE
      ! jpp, 3/2/10
      REAL*8                 :: br2_kg_emiss, salt_current

      ! Increment of radius for Emission integration (um)
      REAL*8, PARAMETER      :: DR    = 5.d-2
      REAL*8, PARAMETER      :: BETHA = 2.d0

      ! External functions
      REAL*8,  EXTERNAL      :: SFCWINDSQR

      !=================================================================
      ! SRCSALT begins here!
      !=================================================================

      ! initialize the Br2 sea-salt emissions
      SSA_Br2(:,:) = 0.d0
      
      ! Emission timestep [s]
      DTEMIS = GET_TS_EMIS() * 60d0

      ! no longer used (jaegle 5/11/11)
      ! Constant [volume * time * other stuff??] 
      !CONST   = 4d0/3d0 * PI * DR * DTEMIS * 1.d-18 * 1.373d0

      !CONST_N = DTEMIS * DR * 1.373d0
      !  Constant for converting from [#/m2/s/um] to [#/m2]
      CONST_N = DTEMIS * (DR * BETHA)

      ! Lower and upper limit of size bin N [um]
      ! Note that these are dry size bins. In order to
      ! get wet (RH=80%) sizes, we need to multiply by
      ! BETHA.
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


            ! Sea salt base source [#/m2]. Note that the Gong formulation
            ! is for r80 (radius at 80% RH), so we need to multiply RMID
            ! by the scaling factor BETHA=2.
            A =  4.7*(1.+30.*(BETHA*RMID(R,N)))
     &                       **(-0.017*(BETHA*RMID(R,N))**(-1.44))
            B =  (0.433d0-LOG10(BETHA*RMID(R,N))) / 0.433d0
            SRC_N(R,N) = CONST_N * 1.373 * (1.d0/(BETHA*RMID(R,N))**(A))
     &           * (1.d0+0.057d0*(BETHA*RMID(R,N))**3.45d0)
     &           * 10d0**(1.607d0*EXP(-(B**2)))

            ! Sea salt base source [kg/m2]: multiply the number of particles
            ! by the dry volume multiplied by the dry density of sea-salt.
            SRC(R,N)  =   SRC_N(R,N) * 4d0/3d0 * PI * 1.d-18
     &           *  SS_DEN( N )  *  (RMID(R,N))**3


            !-----------------------------------------------------------
            ! IMPORTANT NOTE!
            !
            ! In mathematics, "LOG" means "log10". 
            ! In Fortran,     "LOG" means "ln" (and LOG10 is "log10").
            !
            ! The following equations require log to the base 10, so 
            ! we need to use the Fortran function LOG10 instead of LOG. 
            ! (jaegle, bmy, 11/23/09)
            !-----------------------------------------------------------

!            ! Old Monahan et al. (1986) formulation
!            ! Sea salt base source [kg/m2]
!            CONST_N = DTEMIS * (DR * BETHA)
!            SRC(R,N)  = CONST * SS_DEN( N )
!     &           * ( 1.d0 + 0.057d0*( BETHA * RMID(R,N) )**1.05d0 )
!     &           * 10d0**( 1.19d0*
!     &             EXP(-((0.38d0-LOG10(BETHA*RMID(R,N)))/0.65d0)**2))
!     &           / BETHA**2

!            ! Sea salt base source [#/m2] (bec, bmy, 4/13/05)
!            SRC_N(R,N) = CONST_N * (1.d0/RMID(R,N)**3)
!     &           * (1.d0+0.057d0*(BETHA*RMID(R,N))**1.05d0)
!     &           * 10d0**(1.19d0*EXP(-((0.38d0-LOG10(BETHA*RMID(R,N)))
!     &           /0.65d0)**2))/ BETHA**2

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
!$OMP+PRIVATE( SST , SCALE         ) !(jaegle 5/11/11)
!$OMP+PRIVATE( salt_current, br2_kg_emiss )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO J = 1, JJPAR
            
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Grid box surface area [m2]
            A_M2 = GET_AREA_M2( I, J, 1 )

            ! Test if this is a water box
            IF ( IS_WATER( I, J, State_Met ) ) THEN

               ! Wind speed at 10 m altitude [m/s]
               W10M = SQRT( SFCWINDSQR( I, J,
     &                                  State_Met%U10M,
     &                                  State_Met%V10M) )

               ! Sea surface temperature in Celcius (jaegle 5/11/11)
               SST = State_Met%TSKIN(I,J) - 273.15d0

               ! Limit SST to 0-30C range
               SST = MAX( SST , 0d0 )  ! limit to  0C
               SST = MIN( SST , 30d0 ) ! limit to 30C

               ! Empirical SST scaling factor (jaegle 5/11/11)
               SCALE = 0.329d0 + 0.0904d0*SST -
     &                 0.00717d0*SST**2d0 + 0.000207d0*SST**3d0

               ! Reset to using original Gong (2003) emissions (jaegle 6/30/11)
               !SCALE = 1.0d0

               ! The source function calculated with GEOS-4 2x2.5 wind speeds
               ! is too high compared to GEOS-5 at the same resolution. The 10m
               ! winds in GEOS-4 are too rapid. To correct this, apply a global
               ! scaling factor of 0.72 (jaegle 5/11/11)
#if   defined( GEOS_4 )
                SCALE = SCALE * 0.72d0
#endif

               ! Loop over size bins
               DO R = 1, NR

                  ! Update seasalt source into SALT [kg]
                  SALT(I,J)   = SALT(I,J) + 
     &               ( SCALE * SRC(R,N)   * A_M2 * W10M**3.41d0 )

                  ! Update seasalt source into SALT_N [#] 
                  ! (bec, bmy, 4/13/05)
                  SALT_N(I,J) = SALT_N(I,J) +
     &               ( SCALE * SRC_N(R,N) * A_M2 * W10M**3.41d0 )

                  ! jpp, 3/3/10: since the SALT arrays are integrations
                  !              I cannot use them for each independent
                  !              radius... that's why I'm getting too
                  !              much bromine. So I'm making a tmp
                  !              array to store only the current
                  !              Dry Radius bin. [kg]
                  salt_current = ( SRC(R,N)   * A_M2 * W10M**3.41d0 )

                  ! --------------------------------------------------
                  ! jpp, 3/2/10: Accounting for the bromine emissions
                  ! now. Store mass flux [kg] of Br2 based on how much
                  ! aerosol there is emitted in this box.
                  ! --------------------------------------------------
                  CALL EMISS_SSA_BROMINE( J, RMID(R, N), 
     &                  salt_current, br2_kg_emiss )

                  SSA_Br2(I,J) = SSA_Br2(I,J) + br2_kg_emiss

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

            ! Move surface emis to vdiff_mod.f for non-local PBL mixing
            ! --- emis is saved in emis_save(:,:,:) in SUBROUTINE EMISSSEASALT 
            ! (Lin, 06/09/08) 
            IF ( LNLPBL ) FEMIS = 0.D0

            ! Add seasalt emissions into box (I,J,L) [kg]
            TC(I,J,L)         = TC(I,J,L) + ( FEMIS * SALT(I,J) )

            ! Alkalinity [kg] (bec, bmy, 4/13/05)
            ALK_EMIS(I,J,L,N) = SALT(I,J)

            ! Number density [#/m3] (bec, bmy, 4/13/05)
            N_DENS(I,J,L,N)   = SALT_N(I,J) / State_Met%AIRVOL(I,J,L)

         ENDDO

         ! ND08 diagnostic: sea salt emissions [kg]
         IF ( ND08 > 0 ) THEN
            AD08(I,J,N) = AD08(I,J,N) + SALT(I,J)
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      END SUBROUTINE SRCSALT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_alk
!
! !DESCRIPTION: Subroutine GET\_ALK returns the seasalt alkalinity emitted at
!  each timestep to sulfate\_mod.f for chemistry on seasalt aerosols.
!  (bec, 12/7/04, 11/23/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ALK( I, J, L, ALK1, ALK2, Kt1, Kt2, Kt1N, Kt2N,
     &                    State_Met )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACER_MOD,         ONLY : SALA_REDGE_um, SALC_REDGE_um
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I, J, L 
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(OUT)  :: ALK1, ALK2           ! [kg]
      REAL*8, INTENT(OUT)  :: Kt1, Kt2, Kt1N, Kt2N ! [s-1]
!
! !REVISION HISTORY: 
!  (1 ) Becky Alexander says we can remove AREA1, AREA2 (bec, bmy, 9/5/06)
!  (2 ) Bug fix to remove a double-substitution.  Replace code lines for 
!        TERM{123}A, TERM{123}B, TERM{123}AN, TERM{123}BN. (bec, bmy, 7/18/08)
!  (3 ) Updated hygroscopic growth parameters (bec, bmy, 11/23/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
      REAL*8,  PARAMETER :: gamma_SO2 = 0.11d0 !from Worsnop et al. (1989)
      REAL*8,  PARAMETER :: gamma_HNO3 = 0.2d0 !from JPL [2001] 
      REAL*8,  PARAMETER :: Dg = 0.2d0 !gas phase diffusion coeff. [cm2/s]
      REAL*8,  PARAMETER :: v = 3.0d4  !cm/s

      LOGICAL, SAVE      :: FIRST = .TRUE.
 
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
      !ALK1 = STT(I,J,L,IDTSALA) * State_Met%AD(I,J,L)/TCVV(IDTSALA)
      !ALK2 = STT(I,J,L,IDTSALC) * State_Met%AD(I,J,L)/TCVV(IDTSALC)
      !-----------------------------------------------------------------------

      ! Conversion from [m-3] --> [cm-3]
      N1 = N_DENS(I,J,L,1) * 1.d-6
      N2 = N_DENS(I,J,L,2) * 1.d-6

      ALK = ALK1 + ALK2

      ! If there is any alkalinity ...
      IF ( ALK > MINDAT ) THEN

         ! set humidity index IRH as a percent
         IRH = State_Met%RH(I,J,L)
         IRH = MAX(  1, IRH )
         IRH = MIN( 99, IRH )

         ! Hygroscopic growth factor for sea-salt from Chin et al. (2002)
         ! Updated (bec, bmy, 11/23/09)
         IF ( IRH < 100 ) HGF = 4.8d0
         IF ( IRH < 99  ) HGF = 2.9d0
         IF ( IRH < 95  ) HGF = 2.4d0
         IF ( IRH < 90  ) HGF = 2.0d0
         IF ( IRH < 80  ) HGF = 1.8d0
         IF ( IRH < 70  ) HGF = 1.6d0
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

      END SUBROUTINE GET_ALK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_seasalt
!
! !DESCRIPTION: Subroutine INIT\_SEASALT initializes and zeroes all module
!  arrays (bmy, 4/26/04, 4/13/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_SEASALT( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root
      TYPE(OptInput), INTENT(IN)  :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  (1 ) Now exit if we have allocated arrays before.  Now also allocate 
!        ALK_EMIS & N_DENS.  Now reference CMN_SIZE. (bec, bmy, 4/13/05)
!  (2 ) Added SALT_V and DMID (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  04 Mar 2013 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC as arguments
!   4 Mar 2013 - R. Yantosca - Now test for DEPSALA, DEPSALC here

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS, N

      !=================================================================
      ! INIT_SEASALT begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

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

      ALLOCATE( SALT_V( NR_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SALT_V' )
      SALT_V = 0d0

      ALLOCATE( DMID( NR_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DMID' )
      DMID = 0d0

      ! Zero the IDDEP array
      IDDEP = 0

      ! Find drydep species in DEPSAV
      IF ( Input_Opt%LDRYD ) THEN 
         DO N = 1, Input_Opt%NUMDEP
            SELECT CASE ( TRIM( Input_Opt%DEPNAME(N) ) )
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
      ENDIF

      ! Reset IS_INIT
      IS_INIT = .TRUE.

      END SUBROUTINE INIT_SEASALT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_seasalt
!
! !DESCRIPTION: Subroutine CLEANUP\_SEASALT deallocates all module arrays 
!  (bmy, 4/26/04, 4/13/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_SEASALT
!
! !REVISION HISTORY:
!  (1 ) Now deallocates ALK_EMIS, N_DENS, SRC_N (bec, bmy, 4/13/05)
!  (2 ) Deallocated SALT_V and DMID (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
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
      IF ( ALLOCATED( SALT_V   ) ) DEALLOCATE( SALT_V   )
      IF ( ALLOCATED( DMID     ) ) DEALLOCATE( DMID     )

      END SUBROUTINE CLEANUP_SEASALT
!EOC
#if   defined( TOMAS )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcsalt30
!
! !DESCRIPTION: Subroutine SRCSALT30 emits sea-salt into the 30-bin sea-salt 
!  mass and aerosol number arrays.  Sea-salt emission parameterization of 
!  Clarke et al. [2006] (win, 7/17/09)
!\\
!\\
! !INTERFACE:

      SUBROUTINE SRCSALT30( TC1, TC2, State_Met )
!
! !USES:
!
      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND59
      USE DIAG_MOD,           ONLY : AD59_NUMB, AD59_SALT
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TOMAS_MOD,          ONLY : IBINS, Xk
      USE TRACERID_MOD             ! IDTNK1, IDTSS1
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      ! TC1 : Aerosol number tracer array [no.]
      ! TC2 (REAL*8 ) : Sea salt tracer array [kg]
      REAL*8,  INTENT(INOUT) :: TC1(IIPAR,JJPAR,LLPAR,IBINS)
      REAL*8,  INTENT(INOUT) :: TC2(IIPAR,JJPAR,LLPAR,IBINS)
! 
! !AUTHOR:
!  Contact: Win Trivitayanurak (win@cmu.edu)  
! 
!  Arguments as Input/Output:
!  ============================================================================
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Clarke, A.D., Owens, S., Zhou, J. " An ultrafine sea-salt flux from 
!        breaking waves: Implications for CCN in the remote marine atmosphere"
!        JGR, 2006
! 
! !REVISION HISTORY:
!  (1 ) Originally from emisnaN3clarke.f in GISS GCM-II' (win, 7/18/07)
!  (2 ) Now partition emission throughout the PBL (win, 7/18/07)
!  (3 ) Add COEF to adjust emission in a 1x1 nested-grid (win, 4/27/08)
!  16 Feb 2012 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - A_M2 is now a scalar
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I,      J,      L,     K
      INTEGER                :: NTOP
      REAL*4                 :: FOCEAN, W10M,   DTEMIS
      REAL*8                 :: F100,   W,      NUM
      REAL*8                 :: DBIN(IBINS),    A(IBINS)
      REAL*8                 :: A_M2,           FEMIS
       
      ! External functions
      REAL*8,  EXTERNAL      :: SFCWINDSQR

      ! Coefficient to adjust emission in 1x1 grid (win, 4/27/08)
      REAL*8                 :: COEF

#if  defined( TOMAS12 ) || defined( TOMAS15 )

      data Dbin /
# if  defined( TOMAS15 )
     & 0.0d0      , 0.0d0      ,  0.0d0      ,
# endif
     & 9.68859E-09, 1.53797E-08, 2.44137E-08, 3.87544E-08,
     & 6.15187E-08, 9.76549E-08, 1.55017E-07, 2.46075E-07,
     & 3.90620E-07, 6.20070E-07, 9.84300E-07, 3.12500E-06/

      data A /
# if  defined( TOMAS15 )
     &        0.0d0,       0.0d0,        0.0d0,
# endif
     &  4607513.229, 9309031.200, 12961629.010, 13602132.943,
     & 11441451.509, 9387934.311,  8559624.313,  7165322.549,
     &  4648135.263, 2447035.933,  3885009.997,  1006980.679/ ! make same Nk as 30 bins.

#else
!else we are using 30 or 40 bin model

      DATA Dbin /
# if  defined( TOMAS40 )
     & 0.0d0      , 0.0d0      ,  0.0d0      ,  0.0d0      ,  0.0d0   ,
     & 0.0d0      , 0.0d0      ,  0.0d0      ,  0.0d0      ,  0.0d0   ,
# endif
     & 9.68859E-09, 1.22069E-08, 1.53797E-08, 1.93772E-08, 2.44137E-08,
     & 3.07594E-08, 3.87544E-08, 4.88274E-08, 6.15187E-08, 7.75087E-08,
     & 9.76549E-08, 1.23037E-07, 1.55017E-07, 1.95310E-07, 2.46075E-07,
     & 3.10035E-07, 3.90620E-07, 4.92150E-07, 6.20070E-07, 7.81239E-07,
     & 9.84300E-07, 1.24014E-06, 1.56248E-06, 1.96860E-06, 2.48028E-06,
     & 3.12496E-06, 3.93720E-06, 4.96056E-06, 6.24991E-06, 7.87440E-06/
      DATA A /
# if  defined( TOMAS40 )
     & 0.0d0      , 0.0d0      ,  0.0d0      ,  0.0d0      ,  0.0d0   ,
     & 0.0d0      , 0.0d0      ,  0.0d0      ,  0.0d0      ,  0.0d0   ,
# endif
     & 1719793.975, 2887719.254, 4086059.079, 5222972.121, 6172287.155,
     & 6789341.855, 6954290.435, 6647842.508, 6030292.470, 5411159.039,
     & 4920485.633, 4467448.678, 4379031.834, 4180592.479, 3836983.331,
     & 3328339.218, 2675909.440, 1972225.823, 1384692.112, 1062343.821,
     & 913194.1118, 859176.8257, 812688.4300, 719215.3301, 580735.2991,
     & 418247.5535, 273217.6572, 183340.5653, 132174.9032,      0.0000/

#endif
     
      !=================================================================
      ! SRCSALT30 begins here!
      !=================================================================

      ! Depending on the grid resolution. 4x5 (default) doesn't need adjusting coeff
#if   defined( GRID4x5  )
      COEF = 1.d0
#elif defined( GRID2x25 )
      COEF = 1.d0
#elif defined( GRID1x1  ) && defined( NESTED_CH )
      COEF = 0.77d0
      ! Note: monthly emission over the China nested-grid domain is about
      !       1.3 times higher compared to the same domain in 4x5 resolution
      !       Thus applying 1/1.30 = 0.77 factor to correct the emission.
#else
      CALL ERROR_STOP('Adjust seasalt emiss coeff for your model res.?',
     &                'SRCSALT30: seasalt_mod.f')
#endif


      ! Emission timestep [s]
      DTEMIS = GET_TS_EMIS() * 60d0

      ! Loop over grid cells
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]
         A_M2  = GET_AREA_M2( I, J, 1 )
         
         ! Check if over ocean assuming only gridcells that are
         ! at least 50% water are oceans
         FOCEAN   = 1d0 - State_Met%FRCLND(I,J) 
         IF (FOCEAN > 0.5) THEN 

            ! Wind speed at 10 m altitude [m/s]
!            W10M = SQRT( SFCWINDSQR( I, J ) )
            W10M = SQRT( SFCWINDSQR( I, J ,
     &                               State_Met%U10M,
     &                               State_Met%V10M  ) )

            ! in ocean area - calc wind speed/eqm conc
            ! calculate the fraction of whitecap coverage
            W = 3.84E-6 * W10M ** (3.41)

            ! Loop over bins
            DO  K = 1, IBINS

               F100 = A(K)
               
               !===============================================================
               ! Calculate sea-salt emission 
               !===============================================================
               NUM = F100 * W * A_M2 * FOCEAN * DTEMIS * COEF

               !===============================================================
               ! Partition sea-salt emissions through boundary layer
               !===============================================================

               ! Layer in which the PBL top occurs
               NTOP = CEILING( GET_PBL_TOP_L( I, J ) )
             
               ! Loop thru the boundary layer
               DO L = 1, NTOP

                  ! Fraction of the PBL spanned by box (I,J,L) [unitless]
                  FEMIS = GET_FRAC_OF_PBL( I, J, L )

               !================================================================
               ! Add sea-salt number to the tracer array
               !================================================================

                  TC1(I,J,L,K) = TC1(I,J,L,K) + ( NUM * FEMIS )
                  TC2(I,J,L,K) = TC2(I,J,L,K) + 
     &                           NUM * SQRT( Xk(K) * Xk(K+1)) * FEMIS

!                  TC1(I,J,1,K) = TC1(I,J,1,K) + NUM
!                  TC2(I,J,1,K) = TC2(I,J,1,K) + 
!     &                           NUM * SQRT( Xk(K) * Xk(K+1))

               !===============================================================
               ! debug - sfarina used for finding the source of a NaN
               !===============================================================
!
!                  IF( IT_IS_NAN( TC1(I,J,L,K) ) .OR.
!     &                IT_IS_NAN( TC2(I,J,L,K) )  ) THEN
!                     PRINT *,'+++++++ Found NaN in SRCSALT30  +++++++'
!                     PRINT *,'Location (I,J,L):',I,J,L,'Bin',K
!                     CALL  ERROR_STOP('SRCSALT30 ','seasalt_mod.f')
!                  ENDIF

               ENDDO

               !==============================================================
               ! ND59 Diagnostic: Sea salt emission in [kg/box/timestep]       
               !==============================================================
               IF ( ND59 > 0) THEN
                  AD59_NUMB(I,J,1,k) = AD59_NUMB(I,J,1,k) + NUM
                  AD59_SALT(I,J,1,k) = AD59_SALT(I,J,1,k) + 
     &                                       NUM*sqrt(xk(k)*xk(k+1))
               ENDIF

            ENDDO
            
         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE SRCSALT30
!EOC
#endif
#if   defined( APM )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: wet_settlingbin
!
! !DESCRIPTION: Wet settling routine for APM microphysics (by Gan Luo)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WET_SETTLINGBIN( TC, State_Met, State_Chm )
!
! !USES:
!
      USE APM_INIT_MOD,       ONLY : RSALT
      USE APM_DRIV_MOD,       ONLY : GFTOT3D,DENWET3D
      USE APM_INIT_MOD,       ONLY : NCTSEA
      USE APM_INIT_MOD,       ONLY : IDTCTSEA
      USE DRYDEP_MOD,         ONLY : DEPSAVSIZE
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACER_MOD,         ONLY : N_TRACERS, N_APMTRA

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND44
      USE CMN_GCTM_MOD             ! g0
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NSALTBIN)
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY: 
!   9 Feb 2012 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  17 Oct 2013 - M. Sulprizio- Now pass STT via the State_Chm object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: N, IDTEMP
      INTEGER                :: I,      J,     L
      REAL*8                 :: DELZ,   DELZ1, REFF
      REAL*8                 :: P,      DP,    PDP,      TEMP        
      REAL*8                 :: CONST,  SLIP,  VISC
      REAL*8                 :: RWET, RHO, DTCHEM
      REAL*8                 :: TC0(LLPAR), VTS(LLPAR), NAREA, TAREA

      REAL*8,DIMENSION(IIPAR,JJPAR,LLPAR) :: MASS0, MASS1
      REAL*8,DIMENSION(IIPAR,JJPAR,LLPAR,NSALTBIN) :: MASS

      ! Make a pointer to the tracer array
      REAL*8, POINTER :: STT(:,:,:,:)

      ! Point to STT
      STT => State_Chm%TRACERS
      
      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      IF(NCTSEA>=1)THEN
        MASS0=0.D0
        MASS=0.D0
        MASS1=0.D0
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
        DO L=1,LLPAR
        DO J=1,JJPAR
        DO I=1,IIPAR
          MASS0(I,J,L)=SUM(TC(I,J,L,1:NSALTBIN))
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, VTS, P, TEMP, RWET ) 
!$OMP+PRIVATE( RHO, DP, PDP, CONST, SLIP, VISC )
!$OMP+PRIVATE( DELZ, DELZ1, REFF, NAREA, TAREA, TC0 )
!$OMP+SCHEDULE( DYNAMIC )
      DO N=1,NSALTBIN

      ! Seasalt effective radius (i.e. midpt of radius bin) [m]
      REFF = RSALT(N)

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
            TEMP    = State_Met%T(I,J,L)

            RWET    = REFF * GFTOT3D(I,J,L,2)

            ! Density of the wet aerosol (kg/m3)
            RHO     = DENWET3D(I,J,L,2) * 1000.D0 !luo change from g/cm3

            ! Dp = particle diameter [um]
            DP      = 2.d0 * RWET * 1.d6        

            ! PdP = P * dP [hPa * um]
            PDp     = P * Dp

            ! Constant
            CONST   = 2.d0 * RHO * RWET**2 * g0 / 9.d0

            ! Slip correction factor (as function of P*dp)
            Slip = 1.d0+(15.60d0 + 7.0d0 * EXP(-0.059d0 * PDp)) / PDp

            ! Viscosity [Pa*s] of air as a function of temperature 
            VISC = 1.458d-6 * (Temp)**(1.5d0) / ( Temp + 110.4d0 )

            ! Settling velocity [m/s]
            VTS(L) = CONST * Slip / VISC

         ENDDO

         ! Method is to solve bidiagonal matrix which is
         ! implicit and first order accurate in z (rjp, 1/24/02)

         ! We know the boundary condition at the model top
         DO L = 1, LLPAR
            TC0(L) = TC(I,J,L,N)
         ENDDO

         ! We know the boundary condition at the model top
         DO L = LLTROP, 2, -1
            DELZ        = State_Met%BXHEIGHT(I,J,L)
            TC(I,J,L,N) = TC0(L) /
     &      ( 1.D0 + DTCHEM * VTS(L) / DELZ )
            MASS(I,J,(L-1),N) = TC0(L) *
     &      (1.D0-1.D0/( 1.D0 + DTCHEM * VTS(L) / DELZ ))
         ENDDO

         IF(NCTSEA>=1)THEN
         DO L=1,LLPAR
            MASS1(I,J,L)=MASS1(I,J,L)+TC(I,J,L,N)
         ENDDO
         ENDIF

         DO L = LLTROP, 1, -1
            TC(I,J,L,N) = TC(I,J,L,N) + MASS(I,J,L,N)
         ENDDO

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(NCTSEA>=1)THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, IDTEMP )
!$OMP+SCHEDULE( DYNAMIC )
      DO N=1,NCTSEA
        IDTEMP=IDTCTSEA+N-1
      DO J=1,JJPAR
      DO I=1,IIPAR

      DO L=1,LLTROP-1
      IF(MASS0(I,J,L)>1.D-30)THEN
        STT(I,J,L,IDTEMP)=STT(I,J,L,IDTEMP)*MASS1(I,J,L)/MASS0(I,J,L)
      ENDIF
      IF(MASS0(I,J,L+1)>1.D-30)THEN
        STT(I,J,L,IDTEMP)=STT(I,J,L,IDTEMP)+
     &  STT(I,J,L+1,IDTEMP)*SUM(MASS(I,J,L,1:NSALTBIN))/MASS0(I,J,L+1)
      ENDIF
      ENDDO

      L=LLTROP
      IF(MASS0(I,J,L)>1.D-30)THEN
        STT(I,J,L,IDTEMP)=STT(I,J,L,IDTEMP)*MASS1(I,J,L)/MASS0(I,J,L)
      ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Clear the pointer
      NULLIFY( STT )

      END SUBROUTINE WET_SETTLINGBIN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dry_deposition bin
!
! !DESCRIPTION: Dry deposition routine for APM microphysics
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DRY_DEPOSITIONBIN( TC, State_Chm )
!
! !USES:
!
      USE DRYDEP_MOD,         ONLY : DEPSAVSIZE 
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TIME_MOD,           ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACER_MOD,         ONLY : N_TRACERS, N_APMTRA
      USE APM_INIT_MOD,       ONLY : NSO4,NCTSEA
      USE APM_INIT_MOD,       ONLY : IDTCTSEA

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_GCTM_MOD             ! g0
      USE CMN_DIAG_MOD             ! ND44
!
! !INPUT/OUTPUT ARGUMENTS:
!
      REAL*8,         INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NSALTBIN)
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY: 
!   9 Feb 2012 - R. Yantosca - Added ProTeX headers
!  17 Oct 2013 - M. Sulprizio- Now pass STT via the State_Chm object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8,  DIMENSION(IIPAR,JJPAR) :: MASS0, MASS

      ! Local variables
      INTEGER                :: N, IDTEMP
      INTEGER                :: I,        J,     L,      DTCHEM
      REAL*8                 :: OLD, NEW, NAREA, TAREA

      ! Make a pointer to the tracer array
      REAL*8, POINTER :: STT(:,:,:,:)

      ! Point to STT
      STT => State_Chm%TRACERS

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      IF(NCTSEA>=1)THEN
        MASS0=0.D0
        MASS=0.D0
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
!$OMP+SCHEDULE( DYNAMIC )
        DO J=1,JJPAR
        DO I=1,IIPAR
          MASS0(I,J)=SUM(TC(I,J,1,1:NSALTBIN))
        ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, OLD, NEW, NAREA, TAREA )
!$OMP+SCHEDULE( DYNAMIC )

      DO N=1,NSALTBIN

      ! Loop over latitudes
      DO J = 1, JJPAR
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Old tracer concentration [kg]
            OLD  = TC(I,J,1,N)

            ! New tracer concentration [kg]
            NEW  = OLD * EXP( -DEPSAVSIZE(I,J,(N+NSO4)) * DTCHEM  )

            ! Update tracer array
            TC(I,J,1,N) = NEW 

         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF(NCTSEA>=1)THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
!$OMP+SCHEDULE( DYNAMIC )
        DO J=1,JJPAR
        DO I=1,IIPAR
          MASS(I,J)=SUM(TC(I,J,1,1:NSALTBIN))
        ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

      IF(NCTSEA>=1)THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, IDTEMP)
!$OMP+SCHEDULE( DYNAMIC )
        DO N=1,NCTSEA
          IDTEMP=IDTCTSEA+N-1
        DO J=1,JJPAR
        DO I=1,IIPAR
        IF(MASS0(I,J)>1.D-30)THEN
          STT(I,J,1,IDTEMP)=STT(I,J,1,IDTEMP)
     &    * (MASS(I,J)/MASS0(I,J))
        ENDIF
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Clear the pointer
      NULLIFY( STT )

      END SUBROUTINE DRY_DEPOSITIONBIN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcsaltbin
!
! !DESCRIPTION: SRCSALT routine for APM microphysics
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SRCSALTBIN( TC, State_Met )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_WATER
      USE ERROR_MOD,          ONLY : DEBUG_MSG, ERROR_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE APM_INIT_MOD,       ONLY : DFMSALT9  ! kg m-2 s-1

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND44, ND08
      USE CMN_GCTM_MOD             ! PI
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NSALTBIN)
!
! !REVISION HISTORY: 
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: N
      INTEGER                :: I,     J,      L
      INTEGER                :: NTOP
      REAL*8                 :: W10M,  DTEMIS
      REAL*8                 :: FEMIS, A_M2
      REAL*8                 :: SALT(IIPAR,JJPAR)

      ! Increment of radius for Emission integration (um)
      REAL*8, PARAMETER      :: BETHA = 1.d0

      ! External functions
      REAL*8,  EXTERNAL      :: SFCWINDSQR

      !=================================================================
      ! SRCSALT begins here!
      !=================================================================

      ! Emission timestep [s]
      DTEMIS = GET_TS_EMIS() * 60d0

      DO N=1,NSALTBIN
      ! Initialize source
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         SALT(I,J)   = 0d0
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, A_M2, W10M )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Grid box surface area [m2]
            A_M2 = GET_AREA_M2( I, J, 1 )

            ! Test if this is a water box
            IF ( IS_WATER( I, J, State_Met ) ) THEN

               ! Wind speed at 10 m altitude [m/s]
               W10M = SQRT( SFCWINDSQR(I ,J,
     &                                 State_Met%U10M,
     &                                 State_Met%V10M ) )

               ! Loop over size bins

                  ! Update seasalt source into SALT [kg]
                  SALT(I,J)   = SALT(I,J) + 
                  ! DFMSALT9: Sea-salt mass flux dFM (kg m-2 s-1) at U10 = 9 m/s
     &            DFMSALT9(N)*(W10M/9.0)**3.41d0 *A_M2* DTEMIS

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
            TC(I,J,L,N)       = TC(I,J,L,N) + ( FEMIS * SALT(I,J) )

         ENDDO

      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      ENDDO

      END SUBROUTINE SRCSALTBIN          
!EOC
#endif

      END MODULE SEASALT_MOD
