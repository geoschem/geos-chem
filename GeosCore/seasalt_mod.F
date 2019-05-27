!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: seasalt_mod.F
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
      USE HCO_ERROR_MOD
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
      USE PHYSCONSTANTS

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CHEMSEASALT
      PUBLIC  :: CLEANUP_SEASALT
      PUBLIC  :: INIT_SEASALT
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: WET_SETTLING
      PRIVATE :: CHEM_MOPO
      PRIVATE :: CHEM_MOPI
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC  :: SALT_V
      PUBLIC  :: DMID
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
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  16 Feb 2012 - R. Yantosca - Moved SRCSALT30 to end of module
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  04 Mar 2013 - R. Yantosca - Now call INIT_SULFATE from the init stage
!                              which facilitates connection to GEOS-5 GCM
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  23 Jun 2014 - R. Yantosca - Remove code now made obsolete by HEMCO
!  03 Nov 2014 - C. Keller   - Moved GET_ALK to sulfate_mod.F. 
!  20 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  12 Jun 2015 - R. Yantosca - Removed DRY_DEPOSITION routine, because we
!                              now apply drydep in mixing_mod.F90
!  16 Jun 2016 - C. Miller   - Now define species ID's with Ind_ function
!  17 Jun 2016 - R. Yantosca - Add species ID's as module variables
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
      !========================================================================
      ! Module Variables:
      !
      ! NSALT    : # of sea salt tracers
      ! NR_MAX   : Number of size bins
      ! SMALLNUM : A small number (epsilon) for numerical cutoff
      ! SS_DEN   : Sea salt density [kg/m3]
      ! IDDEP    : Drydep index array for sea salt tracers
      ! REDGE    : Array for edges of seasalt radius bins
      ! DMID     : Array for centers of seasalt radius bins
      ! SRC      : Array for baseline seasalt emission/bin [kg/m2]
      ! SRC_N    : Array for baseline seasalt emission/bin [#/m2]
      ! ALK_EMIS : Array for alkalinity [kg]
      ! N_DENS   : Number density of seasalt emissions [#/m3]
      ! SALT_V   : Log-normal volum size distribution for sea salt
      ! 
      ! Note: sea salt emissions are now calculated in HEMCO following the
      ! original code (hcox_seasalt_mod.F). This has made some of the arrays
      ! obsolete. Alkalinity and number density are now calculated in 
      ! sulfate_mod.F through the HEMCO interface (ckeller, 11/03/2014).
      !=======================================================================-
!
! !DEFINED PARAMETERS:
!
      INTEGER,  PARAMETER   :: NSALT    = 2
      INTEGER,  PARAMETER   :: NR_MAX   = 200
      REAL(fp), PARAMETER   :: SMALLNUM = 1e-20_fp
!
! !PRIVATE TYPES:
!
      ! Arrays
      INTEGER               :: IDDEP (NSALT)
      REAL(fp)              :: SS_DEN(NSALT)

      ! Allocatable arrays
      REAL(fp), ALLOCATABLE :: SALT_V(:    )
      REAL(fp), ALLOCATABLE :: DMID  (:    )
      REAL(fp), ALLOCATABLE :: OCCONV(:,:,:)

      ! Species ID flags (formerly in tracerid_mod.F)
      INTEGER               :: id_MOPO
      INTEGER               :: id_MOPI
      INTEGER               :: id_NK1
      INTEGER               :: id_SALA
      INTEGER               :: id_SALC
      INTEGER               :: id_SS1

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
      SUBROUTINE CHEMSEASALT( am_I_Root, Input_Opt,  State_Met,
     &                        State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
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
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!      
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REMARKS:
!  Dry deposition is now handled in mixing_mod.F90.  We have removed
!  the calls to the DRY_DEPOSITION routine here. (bmy, 6/12/15)

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
!  12 Jun 2015 - R. Yantosca - Drydep is now handled in mixing_mod.F90,
!                              so we can greatly collapse this code
!  17 Jun 2016 - R. Yantosca - Add error checks before calling routines
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  10 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  20 Sep 2016 - R. Yantosca - Bug fix: rewrite invalid ( LMPOA > 0 ) test
!  02 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL           :: prtDebug

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! CHEMSEASALT begins here!
      !=================================================================

      ! Initialize
      RC       =  GC_SUCCESS
      Spc      => State_Chm%Species

      ! Do we have to print debug output?
      prtDebug = ( Input_Opt%LPRT .and. am_I_Root )

      !=================================================================
      ! Maybe someday we should merge these two separate calculations
      ! into one (rjp, 4/3/04)
      !=================================================================

      !=================================================================
      ! Accumulation mode wet settling
      !=================================================================
      IF ( id_SALA > 0 ) THEN
       CALL WET_SETTLING( am_I_Root,  Input_Opt,
     &                    State_Met,  State_Chm,
     &                    State_Diag, Spc(:,:,:,id_SALA),
     &                    1,          RC                  )

         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMSEASALT: WET_SET, Accum' )
         ENDIF
      ENDIF

      !=================================================================
      ! Coarse mode wet settling
      !=================================================================
      IF ( id_SALC > 0 ) THEN
        CALL WET_SETTLING( am_I_Root,  Input_Opt,
     &                     State_Met,  State_Chm,
     &                     State_Diag, Spc(:,:,:,id_SALC),
     &                     2,          RC                  )

         IF ( prtDebug ) THEN 
            CALL DEBUG_MSG( '### CHEMSEASALT: WET_SET, Coarse' )
         ENDIF
      ENDIF

      !=================================================================
      ! Do chemistry for marine organic aerosol tracers 
      !=================================================================
      IF ( Input_Opt%LMPOA ) THEN

         ! Chemistry for hydrophobic MOA
         IF ( id_MOPO > 0 ) THEN
            CALL CHEM_MOPO( am_I_Root, Input_Opt,
     &                      State_Chm, State_Diag, RC )
            IF ( prtDebug ) THEN
               CALL DEBUG_MSG( '### CHEMSEASALT: a CHEM_MOPO' )
            ENDIF
         ENDIF

         ! Chemistry for hydrophilic MOA
         IF ( id_MOPI > 0 ) THEN
            CALL CHEM_MOPI( am_I_Root, Input_Opt,
     &                      State_Chm, State_Diag, RC )
            IF ( prtDebug ) THEN
               CALL DEBUG_MSG( '### CHEMSEASALT: a CHEM_MOPI' )
            ENDIF
         ENDIF

      ENDIF

      ! Free pointer
      Spc => NULL()

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
      SUBROUTINE WET_SETTLING( am_I_Root, Input_Opt,  State_Met,
     &                         State_Chm, State_Diag, TC, 
     &                         N,         RC                     )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD44
#endif
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GC_GRID_MOD,        ONLY : GET_AREA_CM2
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: N           ! 1=accum mode;
                                                   ! 2=coarse mode
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! MeteorologyState
      TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
      REAL(fp),       INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! Sea salt [kg]
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success/failure
!
! !REMARKS:
!  TODO: Declare State_Chm as INTENT(INOUT) and make TC a local pointer
!  to State_Chm%SPECIES.  Pass in the ID field to index State_Chm.
!
! !REVISION HISTORY:
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (3 ) Bug fix: DTCHEM has to be REAL(fp), not integer. (bmy, 9/7/06)
!  (4 ) Now limit relative humidity to [tiny(real(fp)),0.99] range for DLOG
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
!  26 Feb 2015 - E. Lundgren - Replace GET_PCENTER with State_Met%PMID.
!                              Remove dependency on pressure_mod.
!  25 Jan 2016 - E. Lundgren - Update netcdf drydep flux diagnostic
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  31 May 2016 - E. Lundgren - Remove usage of Input_Opt%XNUMOL
!  22 Jun 2016 - M. Yannetti - Replaced TCVV with spec db and phys constant
!  20 Oct 2016 - R. Yantosca - Only print out debug output on the root CPU
!  02 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                  :: PrtDebug
      INTEGER                  :: I,      J,     L,        ND
      REAL(fp)                 :: DELZ,   DELZ1, REFF,     DEN
      REAL(fp)                 :: P,      DP,    PDP,      TEMP        
      REAL(fp)                 :: CONST,  SLIP,  VISC,     FAC1
      REAL(fp)                 :: FAC2,   FLUX,  AREA_CM2, RHB
      ! replace RCM with RUM (radis in micron) jaegle 5/11/11
      REAL(fp)                 :: RUM,    RWET,  RATIO_R
      REAL(fp)                 :: TOT1,   TOT2,  DTCHEM
      REAL(fp)                 :: VTS(LLPAR)  
      REAL(fp)                 :: TC0(LLPAR)
      ! New variables (jaegle 5/11/11)
      REAL(fp)                 :: SW
      REAL(fp)                 :: R0,       R1, NR, DEDGE, SALT_MASS
      REAL(fp)                 :: SALT_MASS_TOTAL, VTS_WEIGHT, DMIDW
      REAL(f8)                 :: WTP 
      INTEGER                  :: ID
      LOGICAL, SAVE            :: FIRST = .TRUE.
      REAL(f8)                 :: RHO, RHO1
      CHARACTER(LEN=255)       :: ErrMsg
      CHARACTER(LEN=255)       :: ThisLoc
!
! !DEFINED PARAMETERS:
!      
      REAL(fp),  PARAMETER     :: C1 =  0.7674e+0_fp 
      REAL(fp),  PARAMETER     :: C2 =  3.079e+0_fp 
      REAL(fp),  PARAMETER     :: C3 =  2.573e-11_fp
      REAL(fp),  PARAMETER     :: C4 = -1.424e+0_fp
      ! Parameters for polynomial coefficients to derive seawater
      ! density. From Tang et al. (1997) (jaegle 5/11/11)
      REAL(fp),  PARAMETER     :: A1 =  7.93e-3_fp
      REAL(fp),  PARAMETER     :: A2 = -4.28e-5_fp
      REAL(fp),  PARAMETER     :: A3 =  2.52e-6_fp
      REAL(fp),  PARAMETER     :: A4 = -2.35e-8_fp
      ! increment of radius for integration of settling velocity (um)
      REAL(fp), PARAMETER      :: DR    = 5.e-2_fp
      ! parameter for convergence
      REAL(f8),  PARAMETER     :: EPSI = 1.0e-4_f8
      ! parameters for assumed size distribution of acc and coarse mode
      ! sea salt aerosols (jaegle 5/11/11)
      ! geometric dry mean diameters (microns)
      REAL(fp),  PARAMETER     :: RG_A = 0.085e+0_fp
      REAL(fp),  PARAMETER     :: RG_C = 0.4e+0_fp
      ! sigma of the size distribution
      REAL(fp),  PARAMETER     :: SIG_A = 1.5e+0_fp
      REAL(fp),  PARAMETER     :: SIG_C = 1.8e+0_fp

      ! Local variables for Input_Opt quantities
      REAL(fp)                 :: SALA_REDGE_um(2)
      REAL(fp)                 :: SALC_REDGE_um(2)

      !=================================================================
      ! WET_SETTLING begins here!
      !=================================================================

      ! Initialize
      RC       = GC_SUCCESS
      prtDebug = ( Input_Opt%LPRT .and. am_I_Root )
      ErrMsg   = ''
      ThisLoc  = 
     &  ' -> at WET_SETTLING (in module GeosCore/seasalt_mod.F)'

      ! Copy fields from INPUT_OPT to local variables for use below
      SALA_REDGE_um = Input_Opt%SALA_REDGE_um 
      SALC_REDGE_um = Input_Opt%SALC_REDGE_um

      ! Chemistry timestep [s]
      DTCHEM        = GET_TS_CHEM()

      ! Sea salt density [kg/m3]
      DEN            = SS_DEN( N )

      ! Seasalt effective radius (i.e. midpt of radius bin) [m]
      SELECT CASE ( N )

         ! Accum mode
         ! add R0 and R1 = edges if the sea salt size bins (jaegle 5/11/11)
         CASE( 1 )
            REFF = 0.5e-6_fp * ( SALA_REDGE_um(1) + SALA_REDGE_um(2) )
            R0 = SALA_REDGE_um(1)
            R1 = SALA_REDGE_um(2)

         ! Coarse mode
         CASE( 2 ) 
            REFF = 0.5e-6_fp * ( SALC_REDGE_um(1) + SALC_REDGE_um(2) )
            R0 = SALC_REDGE_um(1)
            R1 = SALC_REDGE_um(2)
            
      END SELECT

      ! Number of dry radius size bins between lowest radius (accumulation
      ! mode) and largest radii (coarse mode) (jaegle 5/11/11)
      
      NR = INT( ( ( SALC_REDGE_um(2) - SALA_REDGE_um(1) ) / DR ) 
     &                  + 0.5e+0_fp )

      ! Trap potential errors
      IF ( NR > NR_MAX ) THEN
         ErrMsg = 'Too many bins!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Define the volume size distribution of sea-salt. This only has
      ! to be done once. We assume that sea-salt is the combination of a 
      ! coarse mode and accumulation model log-normal distribution 
      ! functions (jaegle 5/11/11)
      !=================================================================
      IF ( FIRST) THEN

        ! Lower edge of 0th bin
	DEDGE=SALA_REDGE_um(1) * 2e+0_fp

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
	   SALT_V(ID) = PI / 6e+0_fp* (DMID(ID)**3) * (
     &         13e+0_fp*exp(-0.5*( LOG(DMID(ID))-LOG(RG_A*2e+0_fp) )
     &         **2e+0_fp/LOG(SIG_A)**2e+0_fp )
     &         /( sqrt(2e+0_fp * PI) * LOG(SIG_A) )  +
     &         0.8e+0_fp*exp(-0.5*( LOG(DMID(ID))-LOG(RG_C*2e+0_fp) )
     &         **2e+0_fp/LOG(SIG_C)**2e+0_fp)
     &         /( sqrt(2e+0_fp * PI) * LOG(SIG_C) )  )
	   ! update the next edge
	   DEDGE = DEDGE + DR*2e+0_fp
        ENDDO

        ! Reset after the first time
        IF ( FIRST ) FIRST = .FALSE.
      ENDIF

      IF ( prtDebug ) CALL DEBUG_MSG('SEASALT: STARTING WET_SETTLING')

      ! Sea salt radius [cm]
      !RCM  = REFF * 100e+0_fp  
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
!$OMP+PRIVATE( DMIDW,    RHO1,  WTP,   SALT_MASS,       ND         )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR       

         ! Initialize 
         DO L = 1, LLPAR
            VTS(L) = 0e+0_fp
         ENDDO

         ! Loop over levels
         DO L = 1, LLPAR

            ! Pressure at center of the level [kPa]
            P       = State_Met%PMID(I,J,L) * 0.1e+0_fp

            ! Temperature [K]
            TEMP    = State_Met%T(I,J,L)

            ! Cap RH at 0.99 
            RHB     = MIN( 0.99e+0_fp, State_Met%RH(I,J,L) * 1e-2_fp )

            ! Safety check (phs, 5/1/08)
            RHB     = MAX( TINY(RHB), RHB           )

            ! Aerosol growth with relative humidity in radius [m] 
            ! (Gerber, 1985)
            !RWET    = 0.01e+0_fp*(FAC1/(FAC2-DLOG(RHB))+RCM**3.e+0_fp)**0.33e+0_fp
	    ! Fix bugs in the Gerber formula:  a log10 (instead of ln) 
            ! should be used and the dry radius should be expressed in 
            ! micrometers (instead of cm) also add more significant 
	    ! digits to the exponent (should be 1/3) (jaegle 5/11/11)
            !RWET    = 1d-6*(FAC1/(FAC2-LOG10(RHB))+RUM**3.e+0_fp)**0.33333e+0_fp

            ! Use equation 5 in Lewis and Schwartz (2006) for sea 
            ! salt growth (bec, jaegle 5/11/11)
            RWET = REFF * (4.e+0_fp / 3.7e+0_fp) *
     &       ( (2.e+0_fp - RHB)/(1.e+0_fp - RHB) )**(1.e+0_fp/3.e+0_fp)


            ! Ratio dry over wet radii at the cubic power
            RATIO_R = ( REFF / RWET )**3.e+0_fp

            ! Density of the wet aerosol (kg/m3)
            RHO = RATIO_R * DEN + ( 1.e+0_fp - RATIO_R ) * 1000.e+0_fp

            ! Above density calculation is chemically unsound because 
            ! it ignores chemical solvation.   Iteratively solve Tang et al., 
            ! 1997 equation 5 to calculate density of wet aerosol (kg/m3) 
            ! (bec, jaegle 5/11/11)
            RATIO_R = ( REFF / RWET )
            ! Assume an initial density of 1000 kg/m3
            RHO  = 1000.e+0_f8
            RHO1 = 0.e+0_f8 !initialize (bec, 6/21/10)
            DO WHILE ( ABS( RHO1-RHO ) .gt. EPSI )
                ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet) 
                WTP    = 100.e+0_f8 * DEN/RHO * RATIO_R**3.e+0_f8
                ! Then calculate density of wet aerosol using equation 5 
                ! in Tang et al., 1997 [kg/m3]
                RHO1 = (0.9971e+0_f8 + (A1 * WTP) + (A2 * WTP**2.e+0_f8)
     $               + (A3 * WTP**3.e+0_f8) + (A4 * WTP**4.e+0_f8) ) 
     $               * 1000.e+0_f8
                ! Now calculate new weight percent using above density 
                ! calculation
                WTP    = 100.e+0_f8 * DEN/RHO1 * RATIO_R**3.e+0_f8
                ! Now recalculate new wet density [kg/m3]
                RHO = (0.9971e+0_f8 + (A1 * WTP) + (A2 * WTP**2.e+0_f8)
     $              + (A3 * WTP**3.e+0_f8) + (A4 * WTP**4.e+0_f8) ) 
     $              * 1000.e+0_f8
            ENDDO


            ! Dp = particle diameter [um]
            DP      = 2.e+0_fp * RWET * 1.e+6_fp        

            ! PdP = P * dP [hPa * um]
            PDp     = P * Dp

            ! Constant
            CONST   = 2.e+0_fp * RHO * RWET**2 * g0 / 9.e+0_fp

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
            Slip = 1.e+0_fp+(15.60e+0_fp + 7.0e+0_fp 
     &            * EXP(-0.059e+0_fp * PDp)) / PDp

            ! Viscosity [Pa*s] of air as a function of temperature 
            VISC = 1.458e-6_fp * (Temp)**(1.5e+0_fp) 
     &            / ( Temp + 110.4e+0_fp )

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

	    SALT_MASS_TOTAL = 0e+0_fp
	    VTS_WEIGHT      = 0e+0_fp
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
	       IF ( DMID(ID) .ge. R0*2e+0_fp   .and. 
     &              DMID(ID) .le. R1*2e+0_fp ) THEN
	         DMIDW = DMID(ID) * RWET/REFF  ! wet radius [um]
	         SALT_MASS = SALT_V(ID) * RWET/REFF * RHO 
     &                      / (DMIDW*0.5e+0_fp)
	         VTS_WEIGHT  = VTS_WEIGHT + 
     &              SALT_MASS * VTS(L) * (DMIDW/(RWET*1d6*2e+0_fp) )
     &              ** 2e+0_fp * (2e+0_fp * DR *  RWET/REFF)
	         SALT_MASS_TOTAL=SALT_MASS_TOTAL+SALT_MASS *
     &                            (2e+0_fp * DR *  RWET/REFF)
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
         L    = LLCHEM
         DELZ = State_Met%BXHEIGHT(I,J,L)

         TC(I,J,L) = TC(I,J,L) / ( 1.e+0_fp + DTCHEM * VTS(L) / DELZ )

         DO L = LLCHEM-1, 1, -1
            DELZ  = State_Met%BXHEIGHT(I,J,L)
            DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
            TC(I,J,L) = 1.e+0_fp / ( 1.e+0_fp + DTCHEM * VTS(L) / DELZ )
     &                * ( TC(I,J,L) + DTCHEM * VTS(L+1) / DELZ1
     &                *  TC(I,J,L+1) )
         ENDDO

         !==============================================================
         ! DIAGNOSTIC: Sea salt loss to dry deposition [molec/cm2/s]
         !==============================================================
#if defined( BPCH_DIAG )
         IF ( ND44 > 0 ) THEN
            !-----------------------------------------------------------
            ! ND44 DIAGNOSTIC (bpch)
            ! Dry deposition flux loss [molec/cm2/s]
            !
            ! NOTE: Bpch diagnostics are being phased out.
            !-----------------------------------------------------------
            ! Initialize
            TOT1 = 0e+0_fp
            TOT2 = 0e+0_fp

            ! Compute column totals of TCO(:) and TC(I,J,:,N)
            DO L = 1, LLPAR
               TOT1 = TOT1 + TC0(L)
               TOT2 = TOT2 + TC(I,J,L)
            ENDDO

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! Convert sea salt flux from [kg/s] to [molec/cm2/s]
            FLUX     = ( TOT1 - TOT2 ) / DTCHEM
            FLUX     = FLUX * AVO / ( AIRMW / ( AIRMW 
     &                 / State_Chm%SpcData(id_SALA)%Info%emMW_g )
     &                 * 1.e-3_fp ) / AREA_CM2 

            ND             = IDDEP(N)
            AD44(I,J,ND,1) = AD44(I,J,ND,1) + FLUX
         ENDIF
#endif

         !-----------------------------------------------------------
         ! HISTORY (aka netCDF diagnostics)
         ! Dry deposition flux loss [molec/cm2/s]
         !
         ! NOTE: Eventually think about converting this
         ! diagnostic to more standard units [kg/m2/s]
         !-----------------------------------------------------------
         IF ( State_Diag%Archive_DryDepChm .OR.
     &        State_Diag%Archive_DryDep        ) THEN

            ! Initialize
            TOT1 = 0e+0_fp
            TOT2 = 0e+0_fp

            ! Compute column totals of TCO(:) and TC(I,J,:,N)
            DO L = 1, LLPAR
               TOT1 = TOT1 + TC0(L)
               TOT2 = TOT2 + TC(I,J,L)
            ENDDO

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! Convert sea salt flux from [kg/s] to [molec/cm2/s]
            FLUX     = ( TOT1 - TOT2 ) / DTCHEM
            FLUX     = FLUX * AVO / ( AIRMW / ( AIRMW 
     &                 / State_Chm%SpcData(id_SALA)%Info%emMW_g )
     &                 * 1.e-3_fp ) / AREA_CM2

            ! Drydep index
            ND = IDDEP(N)

            ! Drydep flux in chemistry only
            State_Diag%DryDepChm(I,J,ND) = FLUX
         ENDIF

      ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO

      IF ( prtDebug ) CALL DEBUG_MSG('SEASALT: ENDING WET_SETTLING')

      END SUBROUTINE WET_SETTLING
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_mopo
!
! !DESCRIPTION: Subroutine CHEM\_MOPO modifies hydrophobic marine organic
!  aerosol concentrations based on the conversion to hydrophilic marine 
!  organic aerosols.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_MOPO( am_I_Root, Input_Opt,
     &                      State_Chm, State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  10 Jul 2015 - E. Lundgren - Initial version (based on routine Chem_OCPO)
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  02 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER             :: I,      J,   L
      REAL(fp)            :: DTCHEM, KOC, TC0, CNEW, RKT, FREQ

      ! Pointers
      REAL(fp), POINTER   :: Spc(:,:,:,:)
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER :: OC_LIFE = 1.15e+0_fp

      !=================================================================
      ! CHEM_MOPO begins here!
      !=================================================================

      ! Assume success
      RC        = GC_SUCCESS

      ! Initialize
      KOC       = 1.e+0_fp / ( 86400e+0_fp * OC_LIFE )
      DTCHEM    = GET_TS_CHEM()
      OCCONV    = 0e+0_fp

      ! Set pointer to GEOS-Chem tracer array [kg]
      Spc      => State_Chm%Species
      
      !=================================================================
      ! For tracers with dry deposition, the loss rate of dry dep is 
      ! combined in chem loss term.
      !
      ! Conversion from hydrophobic to hydrophilic:  
      ! e-folding time 1.15 days 
      ! ----------------------------------------
      ! Use an e-folding time of 1.15 days or a convertion rate 
      ! of 1.0e-5 /sec. 
      !    Hydrophobic --> Hydrophilic,  k  = 1.0e-5          
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, FREQ, RKT, CNEW )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial OC [kg]
         TC0  = Spc(I,J,L,id_MOPO)

         ! Zero drydep freq
         ! ### NOTE: Remove this later, but need to make
         ! ### sure we don't incur numerical diffs (bmy, 6/12/15)
         FREQ = 0e+0_fp

         ! Amount of MOPO left after chemistry and drydep [kg]
         RKT  = ( KOC + FREQ ) * DTCHEM
         CNEW = TC0 * EXP( -RKT )

         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

         ! Amount of MOPO converted to MOPI [kg/timestep]
         OCCONV(I,J,L) = ( TC0 - CNEW ) * KOC / ( KOC + FREQ )

         ! Store modified OC concentration back in tracer array
         Spc(I,J,L,id_MOPO) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE CHEM_MOPO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_mopi
!
! !DESCRIPTION: Subroutine CHEM\_MOPI modifies hydrophilic marine organic
!  aerosol concentrations based on the conversion from hydrophobic marine 
!  organic aerosols.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_MOPI( am_I_Root, Input_Opt,
     &                      State_Chm, State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root             ! Root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
! !REMARKS:
! 
! !REVISION HISTORY:
!  10 Jul 2015 - E. Lundgren - Initial version (based on routine Chem_OCPI)
!  02 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER  :: I,   J,     L
      REAL(fp) :: TC0, CNEW, CCV

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! CHEM_MOPI begins here!
      !=================================================================

      ! Assume success
      RC        = GC_SUCCESS

      ! Set pointer to GEOS-Chem tracer array [kg]
      Spc      => State_Chm%Species
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, CCV, CNEW )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial H-philic OC [kg]
         TC0 = Spc(I,J,L,id_MOPI)

         ! H-philic OC that used to be H-phobic OC [kg]
         CCV = OCCONV(I,J,L)

         ! Add the amount of converted MOPO to MOPI
         CNEW = TC0 + CCV

         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

         ! Store modified concentration back in tracer array [kg]
         Spc(I,J,L,id_MOPI) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! Cleanup and quit
      !=================================================================

      ! Zero OCCONV array for next timestep
      OCCONV = 0e+0_fp

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE CHEM_MOPI
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
      SUBROUTINE INIT_SEASALT( am_I_Root, Input_Opt,
     &                         State_Chm, State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Diag_Mod,     ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  (1 ) Now exit if we have allocated arrays before.  Now also allocate 
!        ALK_EMIS & N_DENS.  Now reference CMN_SIZE. (bec, bmy, 4/13/05)
!  (2 ) Added SALT_V and DMID (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  04 Mar 2013 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC as arguments
!  04 Mar 2013 - R. Yantosca - Now test for DEPSALA, DEPSALC here
!  12 Jun 2015 - R. Yantosca - Remove obsolete, commented-out code
!  23 Sep 2015 - R. Yantosca - Now accept State_Chm as an argument so that
!                              we can use the species database
!  30 Sep 2015 - R. Yantosca - DD_A_Density is renamed to Density
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  17 Jun 2016 - R. Yantosca - Now define species ID's on the first call
!  13 Oct 2016 - R. Yantosca - Now allocate OCCONV for marine-POA simulations
!  02 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL, SAVE          :: IS_INIT = .FALSE.
      INTEGER                :: AS, N

      ! Objects
      TYPE(Species), POINTER :: SpcInfo

      !=================================================================
      ! INIT_SEASALT begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      ! Define species indices
      id_MOPI = Ind_('MOPI')
      id_MOPO = Ind_('MOPO')
      id_SALA = Ind_('SALA')
      id_SALC = Ind_('SALC')
      id_NK1  = Ind_('NK1' )
      id_SS1  = Ind_('SS1' )

      ! Initialize pointer
      SpcInfo => NULL()

      ALLOCATE( SALT_V( NR_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SALT_V' )
      SALT_V = 0e+0_fp

      ALLOCATE( DMID( NR_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DMID' )
      DMID = 0e+0_fp

      ! Allocate OCCONV only for marine-POA simulations (bmy, 10/13/16)
      IF ( Input_Opt%LMPOA ) THEN
         ALLOCATE( OCCONV( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCCONV' )
         OCCONV = 0e+0_fp
      ENDIF

      ! Zero the IDDEP array
      IDDEP  = 0
      SS_DEN = 0

      ! Find drydep species in DEPSAV
      IF ( Input_Opt%LDRYD ) THEN 

         ! Loop over all species
         DO N = 1, State_Chm%nSpecies
            
            ! Get info about this species from the species database
            SpcInfo => State_Chm%SpcData(N)%Info

            ! Assign parameters to each species
            SELECT CASE ( TRIM( SpcInfo%Name ) )
               CASE ( 'SALA' )
                  IDDEP(1)  = SpcInfo%DryDepID
                  SS_DEN(1) = SpcInfo%Density
               CASE ( 'SALC' )
                  IDDEP(2)  = SpcInfo%DryDepID
                  SS_DEN(2) = SpcInfo%Density
               CASE DEFAULT
                  ! Nothing
            END SELECT        

            ! Free pointer
            SpcInfo => NULL()
         ENDDO
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
!  13 Oct 2016 - R. Yantosca - Bug fix: now deallocate OCCONV array
!EOP
!------------------------------------------------------------------------------
!BOC
!

      !=================================================================
      ! CLEANUP_SEASALT begins here!
      !=================================================================
      IF ( ALLOCATED( SALT_V   ) ) DEALLOCATE( SALT_V   )
      IF ( ALLOCATED( DMID     ) ) DEALLOCATE( DMID     )
      IF ( ALLOCATED( OCCONV   ) ) DEALLOCATE( OCCONV   )

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
      USE DAO_MOD,            ONLY : IS_WATER
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE GC_GRID_MOD,        ONLY : GET_AREA_M2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TOMAS_MOD,          ONLY : IBINS, Xk
#if defined( BPCH_DIAG )
      USE DIAG_MOD,           ONLY : AD59_NUMB, AD59_SALT
#endif
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      ! TC1 : Aerosol number tracer array [no.]
      ! TC2 (REAL(fp) ) : Sea salt tracer array [kg]
      REAL(fp),  INTENT(INOUT) :: TC1(IIPAR,JJPAR,LLPAR,IBINS)
      REAL(fp),  INTENT(INOUT) :: TC2(IIPAR,JJPAR,LLPAR,IBINS)
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
!  10 Mar 2014 - J. Pierce   - Avoid sea salt emissions over ice
!  25 Jul 2014 - R. Yantosca - Remove reference to function SFCWINDSQR
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                  :: I,      J,      L,     K
      INTEGER                  :: NTOP
      REAL*4                   :: FOCEAN, W10M,   DTEMIS
      REAL(fp)                 :: F100,   W,      NUM
      REAL(fp)                 :: DBIN(IBINS),    A(IBINS)
      REAL(fp)                 :: A_M2,           FEMIS
      REAL(fp)                 :: SFCWINDSQR

      ! Coefficient to adjust emission in 1x1 grid (win, 4/27/08)
      REAL(fp)                 :: COEF

#if  defined( TOMAS12 ) || defined( TOMAS15 )

      data Dbin /
# if  defined( TOMAS15 )
     & 0.0e+0_fp      , 0.0e+0_fp      ,  0.0e+0_fp      ,
# endif
     & 9.68859E-09, 1.53797E-08, 2.44137E-08, 3.87544E-08,
     & 6.15187E-08, 9.76549E-08, 1.55017E-07, 2.46075E-07,
     & 3.90620E-07, 6.20070E-07, 9.84300E-07, 3.12500E-06/

      data A /
# if  defined( TOMAS15 )
     &        0.0e+0_fp,       0.0e+0_fp,        0.0e+0_fp,
# endif
     &  4607513.229, 9309031.200, 12961629.010, 13602132.943,
     & 11441451.509, 9387934.311,  8559624.313,  7165322.549,
     &  4648135.263, 2447035.933,  3885009.997,  1006980.679/ ! make same Nk as 30 bins.

#else
!else we are using 30 or 40 bin model

      DATA Dbin /
# if  defined( TOMAS40 )
     & 0.0e+0_fp, 0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp,
     & 0.0e+0_fp, 0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp,
# endif
     & 9.68859E-09, 1.22069E-08, 1.53797E-08, 1.93772E-08, 2.44137E-08,
     & 3.07594E-08, 3.87544E-08, 4.88274E-08, 6.15187E-08, 7.75087E-08,
     & 9.76549E-08, 1.23037E-07, 1.55017E-07, 1.95310E-07, 2.46075E-07,
     & 3.10035E-07, 3.90620E-07, 4.92150E-07, 6.20070E-07, 7.81239E-07,
     & 9.84300E-07, 1.24014E-06, 1.56248E-06, 1.96860E-06, 2.48028E-06,
     & 3.12496E-06, 3.93720E-06, 4.96056E-06, 6.24991E-06, 7.87440E-06/
      DATA A /
# if  defined( TOMAS40 )
     & 0.0e+0_fp, 0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp,
     & 0.0e+0_fp, 0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp,
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
      COEF = 1.e+0_fp
#elif defined( GRID2x25 )
      COEF = 1.e+0_fp
#else
      CALL ERROR_STOP('Adjust seasalt emiss coeff for your model res.?',
     &                'SRCSALT30: seasalt_mod.f')
#endif


      ! Emission timestep [s]
      DTEMIS = GET_TS_EMIS()

      ! Loop over grid cells
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]
         A_M2  = GET_AREA_M2( I, J, 1 )
         
         ! Check if over ocean assuming only gridcells that are
         ! at least 50% water are oceans (J. Pierce, 3/10/14)
         IF ( IS_WATER( I, J, State_Met ) ) THEN
            FOCEAN = 1e+0_fp - State_Met%FRCLND(I,J)
         ELSE
            FOCEAN = 0.e+0_fp
         ENDIF

         IF (FOCEAN > 0.5e+0_fp) THEN 

            ! Wind speed at 10 m altitude [m/s]
            SFCWINDSQR = State_Met%U10M(I,J)**2  
     &                 + State_Met%V10M(I,J)**2
            W10M       = SQRT( SFCWINDSQR )

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
#if defined( BPCH_DIAG )
               IF ( ND59 > 0) THEN
                  AD59_NUMB(I,J,1,k) = AD59_NUMB(I,J,1,k) + NUM
                  AD59_SALT(I,J,1,k) = AD59_SALT(I,J,1,k) + 
     &                                       NUM*sqrt(xk(k)*xk(k+1))
               ENDIF
#endif

            ENDDO
            
         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE SRCSALT30
!EOC
#endif
      END MODULE SEASALT_MOD
