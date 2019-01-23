!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_environment_mod.F90
!
! !DESCRIPTION: Module GC\_ENVIRONMENT\_MOD establishes the runtime 
!  environment for the GEOS-Chem.  It is designed to receive model parameter
!  and geophysical environment information and allocate memory based upon it.
!\\
!\\
!  It provides routines to do the following:
!
! \begin{itemize}
! \item Allocate geo-spatial arrays
! \item Initialize met. field derived type.
! \item Initialize Chemistry, Metorology, Emissions, and Physics States
! \end{itemize}
!
! !INTERFACE: 
!
MODULE GC_Environment_Mod
!
! !USES
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: GC_Allocate_All
  PUBLIC  :: GC_Init_StateObj
  PUBLIC  :: GC_Init_Grid
  PUBLIC  :: GC_Init_Extra
  PUBLIC  :: GC_Init_Regridding
!
! !PRIVATE MEMBER FUNCTIONS:
!
#if defined( TOMAS )
  PRIVATE :: INIT_TOMAS_MICROPHYSICS
#endif
!
! !REMARKS:
!  For consistency, we should probably move the met state initialization
!  to the same module where the met state derived type is contained.
!
! !REVISION HISTORY:
!  26 Jan 2012 - M. Long     - Created module file
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  19 Oct 2012 - R. Yantosca - Removed routine INIT_LOCAL_MET, this is now
!                              handled in Headers/gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_environment_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  28 Aug 2015 - R. Yantosca - Remove Get_nSchm_nSchmBry; stratospheric 
!                              chemistry fields are now read by HEMCO
!  16 Aug 2016 - M. Sulprizio- Rename from gigc_environment_mod.F90 to 
!                              gc_environment_mod.F90. The "gigc" nomenclature
!                              is no longer used.
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC        
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_allocate_all
!
! !DESCRIPTION: Subroutine GC\_ALLOCATE\_ALL allocates all LAT/LON 
!  ALLOCATABLE arrays for global use by the GEOS-Chem either as a standalone 
!  program or module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_Allocate_All( am_I_Root,       Input_Opt,       &
                              RC,              value_I_LO,      &
                              value_J_LO,      value_I_HI,      &
                              value_J_HI,      value_IM,        &
                              value_JM,        value_LM,        &
                              value_IM_WORLD,  value_JM_WORLD,  &
                              value_LM_WORLD,  value_LLSTRAT )
!
! !USES:
!
#if defined( BPCH_DIAG )
    USE CMN_DIAG_Mod,       ONLY : Init_CMN_DIAG
#endif
    USE CMN_FJX_MOD,        ONLY : Init_CMN_FJX
    USE CMN_O3_Mod,         ONLY : Init_CMN_O3
    USE CMN_SIZE_Mod,       ONLY : Init_CMN_SIZE
    USE ErrCode_Mod  
    USE Input_Opt_Mod
    USE VDIFF_PRE_Mod,      ONLY : Init_Vdiff_Pre

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root        ! Are we on the root CPU?
    INTEGER,        OPTIONAL      :: value_I_LO       ! Min local lon index
    INTEGER,        OPTIONAL      :: value_J_LO       ! Min local lat index
    INTEGER,        OPTIONAL      :: value_I_HI       ! Max local lon index
    INTEGER,        OPTIONAL      :: value_J_HI       ! Max local lat index
    INTEGER,        OPTIONAL      :: value_IM         ! Local # of lons
    INTEGER,        OPTIONAL      :: value_JM         ! Local # of lats
    INTEGER,        OPTIONAL      :: value_LM         ! Local # of levels
    INTEGER,        OPTIONAL      :: value_IM_WORLD   ! Global # of lons
    INTEGER,        OPTIONAL      :: value_JM_WORLD   ! Global # of lats
    INTEGER,        OPTIONAL      :: value_LM_WORLD   ! Global # of levels
    INTEGER,        OPTIONAL      :: value_LLSTRAT    ! # of strat. levels
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt        ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC               ! Success or failure?
!
! !REMARKS:
!  For error checking, return up to the main routine w/ an error code.
!  This can be improved upon later.
!
! !REVISION HISTORY: 
!  26 Jan 2012 - M. Long     - Initial version
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  17 Oct 2012 - R. Yantosca - Add am_I_Root, RC as arguments
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Allocate_All
!  30 Oct 2012 - R. Yantosca - Now pass am_I_Root, RC to SET_COMMSOIL_MOD
!  01 Nov 2012 - R. Yantosca - Now zero the fields of the Input Options object
!  16 Nov 2012 - R. Yantosca - Remove this routine from the #ifdef DEVEL block
!  27 Nov 2012 - R. Yantosca - Now pass Input_Opt to INIT_COMODE_LOOP
!  03 Dec 2012 - R. Yantosca - Now pass am_I_Root, RC to INIT_CMN_SIZE
!  03 Dec 2012 - R. Yantosca - Add optional arguments to accept dimension
!                              size information from the ESMF interface
!  13 Dec 2012 - R. Yantosca - Remove reference to obsolete CMN_DEP_mod.F
!  13 Dec 2012 - R. Yantosca - Remove reference to obsolete CMN_mod.F
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete CMN_NOX_mod.F
!  25 Jul 2014 - R. Yantosca - Remove reference to obsolete commsoil_mod.F90
!  25 Jul 2014 - R. Yantosca - Now call INIT_GET_NDEP (GeosCore/get_ndep_mod.F)
!  04 Aug 2015 - C. Keller   - Now pass LLTROP and LLSTRAT to INIT_CMN_SIZE.
!  17 Jun 2016 - R. Yantosca - Move call to INIT_GET_NDEP to GIGC_INIT_EXTRA
!                              which is called after species database init
!  30 Jun 2016 - M. Sulprizio- Remove call to INIT_COMODE_LOOP; it's obsolete
!  20 Dec 2017 - R. Yantosca - Return when encountering errors
!  29 Dec 2017 - C. Keller   - Now accept value of LLSTRAT from Input_Opt
!  14 Mar 2018 - E. Lundgren - Input_Opt parameter is IN only, not INOUT
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODEL_GEOS )
    ! Integers
    INTEGER            :: LLTROP
#endif

    ! Strings
    INTEGER            :: LLSTRAT
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    
    !=======================================================================
    ! GC_Allocate_All begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at GC_Allocate_All  (in module GeosCore/gc_environment_mod.F90)'

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ) || defined( ESMF_ )
    !-----------------------------------------------------------------------
    !          %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
    !
    ! Pass dimension sizes obtained from the ESMF interface to routine 
    ! INIT_CMN_SIZE via several optional arguments (i.e. "value_*").  This
    ! obviates the need to call INIT_CMN_SIZE from a higher level in the
    ! code (i.e. from GIGC_Chunk_Init in ESMF/gigc_chunk_mod.F90).
    ! (bmy, 12/3/12)
    !-----------------------------------------------------------------------

#if defined( MODEL_GEOS )
    LLTROP = 40
    ! 132 layers
    IF ( value_LM==132) LLTROP = 80
#endif

    ! Set dimensions in CMN_SIZE
    CALL Init_CMN_SIZE( am_I_Root      = am_I_Root,       &
                        Input_Opt      = Input_Opt,       &
                        RC             = RC,              &
                        value_I_LO     = value_I_LO,      &
                        value_J_LO     = value_J_LO,      &
                        value_I_HI     = value_I_HI,      &
                        value_J_HI     = value_J_HI,      &
                        value_IM       = value_IM,        &
                        value_JM       = value_JM,        &
                        value_LM       = value_LM,        &
                        value_IM_WORLD = value_IM_WORLD,  &
                        value_JM_WORLD = value_JM_WORLD,  &
                        value_LM_WORLD = value_LM_WORLD,  &
#if defined( MODEL_GEOS )
                        value_LLTROP   = LLTROP,          &
                        value_LLSTRAT  = value_LLSTRAT )
#else
                        value_LLTROP   = 40,              &
                        value_LLSTRAT  = 59            )
#endif

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_CMN_Size"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Exit upon error
    IF ( RC /= GC_SUCCESS ) RETURN

#else
    !-----------------------------------------------------------------------
    !          %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
    !
    ! Current practice in the standard GEOS-Chem is to set dimension sizes
    ! from parameters IGLOB, JGLOB, LGLOB in CMN_SIZE_mod.F.  Therefore,
    ! we do not need to call INIT_CMN_SIZE with optional parameters as is
    ! done when connecting to the ESMF interface.  
    !-----------------------------------------------------------------------

    ! Set dimensions in CMN_SIZE
    CALL Init_CMN_SIZE( am_I_Root, Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_CMN_SIZE"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#endif

#if defined( BPCH_DIAG )
    ! Set dimensions in CMN_DEP_mod.F and allocate arrays
    CALL Init_CMN_DIAG( am_I_Root, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_CMN_DIAG"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

    ! Initialize CMN_O3_mod.F
    CALL Init_CMN_O3( am_I_Root, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_CMN_O3"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Initialize vdiff_pre_mod.F90
    CALL Init_VDIFF_PRE( am_I_Root, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_Vdiff_Pre"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Initialize CMN_FJX_mod.F
    CALL Init_CMN_FJX( am_I_Root, Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_CMN_FJX"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
          
  END SUBROUTINE GC_Allocate_All
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_init_stateobj
!
! !DESCRIPTION: Subroutine GC\_INIT\_STATEOBJ initializes the top-level data 
!  structures that are either passed to/from GC or between GC components 
!  (emis->transport->chem->diagnostics->etc)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_Init_StateObj( am_I_Root,  Diag_List, Input_Opt, State_Chm, &
                               State_Diag, State_Met, RC                    )
!
! !USES:
!
    USE DiagList_Mod
    USE Diagnostics_Mod
    USE ErrCode_Mod
    USE Input_Opt_Mod
    USE State_Chm_Mod
    USE State_Met_Mod
    USE State_Diag_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(DgnList),  INTENT(IN)    :: Diag_List   ! Diagnostics list object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  Need to add better error checking, currently we just return upon error.
!
! !REVISION HISTORY: 
!  26 Jan 2012 - M. Long     - Initial version
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  16 Oct 2012 - R. Yantosca - Renamed LOCAL_MET argument to State_Met
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE  argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Call Init_Chemistry_State (in gc_type2_mod.F90,
!                              which was renamed from INIT_CHEMSTATE)
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_errcode_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference IGAS in Headers/comode_loop_mod.F
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Init_All
!  26 Oct 2012 - R. Yantosca - Now call Get_nSchm, nSchmBry to find out the
!                              number of strat chem species and Bry species
!  01 Nov 2012 - R. Yantosca - Now use LSCHEM from logical_mod.F
!  09 Nov 2012 - R. Yantosca - Now pass Input Options object for GIGC
!  26 Feb 2013 - R. Yantosca - Now pass Input_Opt to Init_GIGC_State_Chm
!  28 Aug 2015 - R. Yantosca - Remove strat-chem options from call to 
!                              Init_GIGC_State_Chm; this is done by HEMCO
!  25 Jan 2016 - R. Yantosca - Bug fix: Declare Input_Opt as INTENT(INOUT),
!                              to match the declaration in INIT_GIGC_STATE_CHM
!  28 Jan 2016 - M. Sulprizio- Remove NBIOMAX from call to Init_GIGC_State_Chm
!  30 Jun 2016 - M. Sulprizio- Remove nSpecies from call to Init_GIGC_State_Chm
!  26 Jun 2017 - R. Yantosca - Now call GC_ERROR to give better error feedback
!  05 Jul 2017 - R. Yantosca - Now initialize the Diagnostics State object
!  26 Sep 2017 - E. Lundgren - Pass diagnostics list object as argument
!  16 Nov 2017 - E. Lundgren - Do not pass IIPAR, JJPAR, LLPAR, NDUST, AER
!                              since available from CMN_Size_Mod in routines
!  01 Feb 2018 - E. Lundgren - Initialize new diagnostics_mod module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    ErrMsg  = ''
    ThisLoc = ' -> at GC_Init_StateObj (in GeosCore/gc_environment_mod.F90)'
    
    !=======================================================================
    ! Initialize the Meteorology State object
    !=======================================================================
    CALL Init_State_Met( am_I_Root   = am_I_Root,   & ! Root CPU (T/F?)
                         State_Met   = State_Met,   & ! Meteorology State
                         RC          = RC          )  ! Success or failure?

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Initialize the Chemistry State object
    !=======================================================================
    CALL Init_State_Chm(  am_I_Root  = am_I_Root,   &  ! Root CPU (Y/N)?
                          Input_Opt  = Input_Opt,   &  ! Input Options
                          State_Chm  = State_Chm,   &  ! Chemistry State
                          RC         = RC          )   ! Success or failure
    
    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Initialize the Diagnostics State object
    !=======================================================================
    CALL Init_State_Diag( am_I_Root  = am_I_Root,   &  ! Root CPU (Y/N)?
                          Input_Opt  = Input_Opt,   &  ! Input Options
                          State_Chm  = State_Chm,   &  ! Chemistry State
                          Diag_List  = Diag_List,   &  ! Diagnostic list obj 
                          State_Diag = State_Diag,  &  ! Chemistry State
                          RC         = RC          )   ! Success or failure

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_State_Diag"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
  END SUBROUTINE GC_Init_StateObj
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_init_grid
!
! !DESCRIPTION: Subroutine GC\_INIT\_GRID calls routines from
!  gc\_grid\_mod.F90 to initialize the horizontal grid parameters.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_Init_Grid( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE GC_GRID_MOD,        ONLY : COMPUTE_GRID
    USE GC_GRID_MOD,        ONLY : INIT_GRID
    USE GC_GRID_MOD,        ONLY : SET_XOFFSET, SET_YOFFSET
    USE Input_Opt_Mod,      ONLY : OptInput
    USE TRANSFER_MOD,       ONLY : INIT_TRANSFER
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  The module grid_mod.F90 has been modified to save grid parameters in 3D
!  format, which will facilitate interfacing GEOS-Chem to a GCM.
!                                                                             .
!  The module global_grid_mod.F90 contains several of the global grid arrays
!  (*_g) originally in grid_mod.F. These arrays are used in regridding GFED3
!  biomass emissions, which are available on a 0.5x0.5 global grid. The global
!  arrays may need to be used in the future for regridding other emissions for
!  nested grids.
!                                                                             .
! !REVISION HISTORY: 
!  01 Mar 2012 - R. Yantosca - Initial version
!  01 May 2012 - M. Payer    - Add call to COMPUTE_GLOBAL_GRID for nested grids
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  30 Nov 2012 - R. Yantosca - Accept external DLON, DLAT from ESMF interface
!  26 Feb 2013 - R. Yantosca - Now pass I_LO, J_LO to COMPUTE_GRID
!  28 Feb 2013 - R. Yantosca - Bug fix for GEOS-5 interface: Now call 
!                              Compute_Grid with 1..IIPAR, I..JJPAR
!  01 Jul 2013 - R. Yantosca - Don't use 1/2 sized polar boxes for GCAP
!  25 Jun 2014 - R. Yantosca - Now accept Input_Opt via the arg list
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete global_grid_mod
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!  26 Jan 2018 - M. Sulprizio- Moved to gc_environment_mod.F90 from input_mod.F
!                              and renamed from Initialize_Grid to GC_Init_Grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: JSP,    JNP

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! Define inputs depending on the grid that is selected
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
       ' -> at GC_Init_Grid (in module GeosCore/gc_environment_mod.F90)'

    ! Set global offsets for the horizontal grid
    CALL SET_XOFFSET( Input_Opt%NESTED_I0 )
    CALL SET_YOFFSET( Input_OPt%NESTED_J0 )

    ! Initialze quantities for transfer_mod.F
    CALL INIT_TRANSFER( 0, 0 )

#if   defined( GRID4x5 )

    !-----------------------------
    ! Global 4 x 5 grid
    !-----------------------------
    JSP           = 1                  ! Lat index of S Pole
    JNP           = JM_WORLD           ! Lat index of N pole
    DLON          = 5.0e+0_fp          ! Delta-longitude array [degrees]
    DLAT          = 4.0e+0_fp          ! Delta-latitude array [degrees]
    DLAT(:,JSP,:) = 2.0e+0_fp          ! Half-sized S Pole boxes
    DLAT(:,JNP,:) = 2.0e+0_fp          ! Half-sized N Pole boxes

#elif defined( GRID2x25 )

    !-----------------------------
    ! Global 2 x 2.5 grid
    !-----------------------------
    JSP           = 1                  ! Lat index of S Pole
    JNP           = JM_WORLD           ! Lat index of N pole
    DLON          = 2.5e+0_fp          ! Delta-longitude array [degrees]
    DLAT          = 2.0e+0_fp          ! Delta-latitude array [degrees]
    DLAT(:,JSP,:) = 1.0e+0_fp          ! Half-sized S Pole boxes
    DLAT(:,JNP,:) = 1.0e+0_fp          ! Half-sized N Pole boxes

#elif defined( GRID05x0625 )

    !-------------------=---------
    ! Nested 0.5 x 0.625 grids
    !-----------------------------
    JSP           = 0                  ! Lat index of S Pole
    JNP           = 0                  ! Lat index of N pole
    DLON          = 0.625e+0_fp        ! Delta-longitude array [degrees]
    DLAT          = 0.5e+0_fp          ! Delta-latitude array [degrees]

#elif defined( GRID025x03125 )

    !-----------------------------
    ! Nested 0.25 x 0.3125 grids
    !-----------------------------
    JSP           = 0                  ! Lat index of S Pole
    JNP           = 0                  ! Lat index of N pole
    DLON          = 0.3125e+0_fp       ! Delta-longitude array [degrees]
    DLAT          = 0.25e+0_fp         ! Delta-latitude array [degrees]

#elif defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )

    !-----------------------------
    ! Connecting to GCM via ESMF
    !-----------------------------
    JSP           = 1
    JNP           = JM_WORLD

    ! NOTE: DLON, DLAT are defined in routine GIGC_Get_Options in
    ! ESMF/gigc_initialization_mod.F90, so no need to define them
    ! here.  (bmy, 12/3/12)

#endif

    !=================================================================
    ! Initialize the horizontal grid
    !=================================================================

    ! Allocate module arrays in grid_mod.F90
    CALL Init_Grid   ( am_I_Root = am_I_Root, &   ! Are we on the root CPU?
                       Input_Opt = Input_Opt, &   ! Input Options object
                       IM        = IIPAR,     &   ! # of lons on this CPU
                       JM        = JJPAR,     &   ! # of lats on this CPU
                       LM        = LLPAR,     &   ! # of levs on this CPU
                       RC        = RC         )   ! Success or failure?
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Grid"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
    ! Compute the horiziontal grid properties
    CALL Compute_Grid( am_I_Root = am_I_Root, &   ! Are we on the root CPU?
                       I1        = 1,         &   ! Min lon index, this CPU
                       J1        = 1,         &   ! Min lat index, this CPU
                       L1        = 1,         &   ! Min lev index, this CPU
                       I2        = IIPAR,     &   ! Max lon index, this CPU
                       J2        = JJPAR,     &   ! Max lat index, this CPU
                       L2        = LLPAR,     &   ! Max lev index, this CPU
                       JSP       = JSP,       &   ! Lat index of South Pole
                       JNP       = JNP,       &   ! Lat index of North Pole
                       DLON      = DLON,      &   ! Delta-longitudes [deg]
                       DLAT      = DLAT,      &   ! Delta-latitudes  [deg]
                       I_LO      = I_LO,      &   ! Min global lon index
                       J_LO      = J_LO,      &   ! Min global lat index
                       RC        = RC         )   ! Success or failure?
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Compute_Grid"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

  END SUBROUTINE GC_Init_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_init_extra
!
! !DESCRIPTION: Suborutine GC\_INIT\_EXTRA initializes other GEOS-Chem 
!  modules that have not been initialized in either GC\_Allocate\_All or
!  GC\_Init\_all.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_Init_Extra( am_I_Root, Diag_List,  Input_Opt, &
                            State_Chm, State_Diag, RC         )
!
! !USES:
!
    USE Aerosol_Mod,        ONLY : Init_Aerosol
    USE Carbon_Mod,         ONLY : Init_Carbon
    USE CO2_Mod,            ONLY : Init_CO2
    USE C2H6_Mod,           ONLY : Init_C2H6
    USE Depo_Mercury_Mod,   ONLY : Init_Depo_Mercury
    USE Diag03_Mod,         ONLY : Init_Diag03
    USE Diag04_Mod,         ONLY : Init_Diag04
    USE Diag20_Mod,         ONLY : Init_Diag20
    USE Diag41_Mod,         ONLY : Init_Diag41
    USE Diag42_Mod,         ONLY : Init_Diag42
    USE Diag48_Mod,         ONLY : Init_Diag48
    USE Diag49_Mod,         ONLY : Init_Diag49
    USE Diag50_Mod,         ONLY : Init_Diag50
    USE Diag51_Mod,         ONLY : Init_Diag51
    USE Diag51b_Mod,        ONLY : Init_Diag51b
    USE Diag53_Mod,         ONLY : Init_Diag53
    USE Diag56_Mod,         ONLY : Init_Diag56
    USE Diag63_Mod,         ONLY : Init_Diag63
    USE Diag_OH_Mod,        ONLY : Init_Diag_OH
    USE DiagList_Mod,       ONLY : DgnList
    USE Drydep_Mod,         ONLY : Init_Drydep
    USE Dust_Mod,           ONLY : Init_Dust
    USE ErrCode_Mod
    USE Error_Mod,          ONLY : Debug_Msg
    USE Gamap_Mod,          ONLY : Do_Gamap
    USE Get_Ndep_Mod,       ONLY : Init_Get_Ndep
    USE Global_CH4_Mod,     ONLY : Init_Global_CH4
    USE Input_Mod,          ONLY : Do_Error_Checks
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Land_Mercury_Mod,   ONLY : Init_Land_Mercury
    USE Mercury_Mod,        ONLY : Init_Mercury
    USE Modis_Lai_Mod,      ONLY : Init_Modis_Lai
    USE Ocean_Mercury_Mod,  ONLY : Init_Ocean_Mercury
    USE POPs_Mod,           ONLY : Init_POPs
    USE Seasalt_Mod,        ONLY : Init_SeaSalt
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE Sulfate_Mod,        ONLY : Init_Sulfate
    USE Tagged_CO_Mod,      ONLY : Init_Tagged_CO
    USE Tagged_O3_Mod,      ONLY : Init_Tagged_O3
    USE Toms_Mod,           ONLY : Init_Toms
    USE TPCORE_BC_Mod,      ONLY : Init_TPCORE_BC
    USE Vdiff_Pre_Mod,      ONLY : Set_Vdiff_Values
    USE WetScav_Mod,        ONLY : Init_WetScav
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(DgnList ), INTENT(IN)    :: Diag_List   ! Diagnostics list object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Several of the INIT routines now called within GC_Init_Extra had 
!  originally been called from the Run method.  We now gather these INIT
!  routines here so that they may be called from the Initialization method.
!  This is necessary when connecting GEOS-Chem to the GEOS-5 GCM via ESMF.
!                                                                             .
! !REVISION HISTORY: 
!  04 Mar 2013 - R. Yantosca - Initial revision
!  05 Mar 2013 - R. Yantosca - Now call INIT_AEROSOL (GeosCore/aerosol_mod.F)
!  15 Mar 2013 - R. Yantosca - Now call INIT_LINOZ (GeosCore/linoz_mod.F)
!  29 Mar 2013 - R. Yantosca - Now call INIT_TROPOPAUSE (so that we can pass
!                              a LVARTROP from Input_Opt and not logical_mod.F)
!  10 Apr 2014 - R. Yantosca - Now call INIT_TAGGED_CO
!  10 Apr 2014 - R. Yantosca - Now call INIT_TAGGED_OX and INIT_GLOBAL_CH4
!  11 Apr 2014 - R. Yantosca - Now call INIT_C2H6 and INIT_HCN_CH3CN
!  14 Apr 2014 - R. Yantosca - Also call INIT_C2H6 if it's a fullchem sim
!                              since we read C2H6 emissions from c2h6_mod.F
!  25 Jun 2014 - R. Yantosca - Now call INIT_MODIS_LAI      
!  25 Jun 2014 - R. Yantosca - Now call SET_VDIFF_VALUES (vdiff_pre_mod.F90)
!  25 Aug 2014 - M. Sulprizio- Now call INIT_POPS
!  16 Mar 2015 - R. Yantosca - Now call INIT_TOMS here
!  28 Aug 2015 - R. Yantosca - Now initialize drydep & wetdep here, so that
!                              we can take advantage of the species database
!  03 Sep 2015 - R. Yantosca - Now call INIT_WETSCAV instead of WETDEPID
!  21 Sep 2015 - R. Yantosca - Now pass State_Chm to INIT_POPS
!  22 Sep 2015 - R. Yantosca - Bug fix: only call INIT_WETSCAV if convection,
!                              wetdep, or chemistry ist turned on.  This
!                              replicates the prior behavior,
!  23 Sep 2015 - R. Yantosca - Now pass State_Chm to INIT_SEASALT
!  23 Sep 2015 - R. Yantosca - Now pass State_Chm to INIT_SULFATE
!  25 Apr 2016 - R. Yantosca - Now call INIT_DIAG03 here
!  25 Apr 2016 - R. Yantosca - Now initialize all mercury modules from here
!  23 Jun 2016 - R. Yantosca - Now call INIT_DIAG_OH from here
!  16 Aug 2016 - M. Sulprizio- Rename from GIGC_Init_Extra to GC_Init_Extra.
!                              The "gigc" nomenclature is no longer used.
!  05 Oct 2017 - R. Yantosca - Now accept State_Diag as an argument and
!                              pass it to the various init routines
!  07 Nov 2017 - R. Yantosca - Now accept Diag_List as an argument
!  07 Nov 2017 - R. Yantosca - Return error condition to main level
!  26 Jan 2018 - M. Sulprizio- Moved to gc_environment_mod.F90 from input_mod.F
!  07 Aug 2018 - H.P. Lin    - Unify init routines to accept Input_Opt, State_Chm,
!                              State_Diag
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: prtDebug
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! GC_Init_Extra begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    prtDebug = ( am_I_Root .and. Input_Opt%LPRT )
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at GC_Init_Extra (in module GeosCore/gc_environment_mod.F90)'

    !=================================================================
    ! Do some error checking before initializing modules
    !=================================================================

    ! Make sure various input options are consistent with the
    ! species that are defined for the simulation
    CALL Do_Error_Checks( am_I_Root, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Do_Error_Checks"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Call setup routines for drydep
    !=================================================================
    IF ( Input_Opt%LDRYD ) THEN

       ! Setup for dry deposition
       CALL Init_DryDep( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Drydep"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
      
       ! Print extra info message for Hg simulation
       IF ( Input_Opt%ITS_A_MERCURY_SIM .and. Input_Opt%LSPLIT ) THEN
          WRITE ( 6, 120 )
          WRITE ( 6, 121 )
       ENDIF
    ENDIF

    ! FORMAT strings
120 FORMAT( /, 'All tagged Hg2 species have the same dep velocity ' &
               'as the total Hg2 species.' )
121 FORMAT( 'All tagged HgP species have the same dep velocity ' &
            'as the total HgP species.' )

    !=================================================================
    ! Call setup routines for wetdep
    !
    ! We need to initialize the wetdep module if either wet 
    ! deposition or convection is turned on, so that we can do the 
    ! large-scale and convective scavenging.  Also initialize the 
    ! wetdep module if both wetdep and convection are turned off, 
    ! but chemistry is turned on.  The INIT_WETSCAV routine will also 
    ! allocate the H2O2s and SO2s arrays that are referenced in the 
    ! convection code. (bmy, 9/23/15)
    !=================================================================
    IF ( Input_Opt%LCONV   .or. &
         Input_Opt%LWETD   .or. &    
         Input_Opt%LCHEM ) THEN
       CALL Init_WetScav( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Wetscav"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF
         
    !=================================================================
    ! Call setup routines from other F90 modules
    !=================================================================

    !-----------------------------------------------------------------
    ! Initialize the MODIS leaf area index module
    !-----------------------------------------------------------------
    CALL Init_Modis_LAI( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Modis_LAI"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Call SET_VDIFF_VALUES so that we can pass several values from 
    ! Input_Opt to the vdiff_mod.F90.  This replaces the functionality
    ! of logical_mod.F and tracer_mod.F..  This has to be called
    ! after the input.geos file has been read from disk.
    !-----------------------------------------------------------------
    CALL Set_Vdiff_Values( am_I_Root, Input_Opt, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Set_Vdiff_Values"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF


    !-----------------------------------------------------------------
    ! Initialize the GET_NDEP_MOD for soil NOx deposition (bmy, 6/17/16)
    !-----------------------------------------------------------------
    CALL Init_Get_Ndep( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Get_Ndep"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Initialize "carbon_mod.F"
    !-----------------------------------------------------------------
    IF ( Input_Opt%LCARB ) THEN
       CALL Init_Carbon( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in ""!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Initialize "dust_mod.F"
    !-----------------------------------------------------------------
    IF ( Input_Opt%LDUST ) then
       CALL Init_Dust( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Dust"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Initialize "seasalt_mod.F
    !-----------------------------------------------------------------
    IF ( Input_Opt%LSSALT ) THEN
       CALL Init_Seasalt( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Seasalt"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Initialize "sulfate_mod.F"
    !-----------------------------------------------------------------
    IF ( Input_Opt%LSULF ) THEN
       CALL Init_Sulfate( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Sulfate"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Initialize "aerosol_mod.F"
    !-----------------------------------------------------------------
    IF ( Input_Opt%LSULF .or. Input_Opt%LCARB    .or. &
         Input_Opt%LDUST .or. Input_Opt%LSSALT ) THEN
       CALL Init_Aerosol( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Aerosol"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Initialize "toms_mod.F"
    !-----------------------------------------------------------------
    CALL Init_Toms( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Toms"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Initialize specialty simulation modules here
    !=================================================================

    !-----------------------------------------------------------------
    ! CO2
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_CO2_SIM ) THEN
       CALL Init_CO2( am_I_Root, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_CO2"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! C2H6 -- note: simulation is obsolete, needs help!
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
         Input_Opt%ITS_A_C2H6_SIM     ) THEN
       CALL Init_C2H6( am_I_Root, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_C2H6"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! CH4
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_CH4_SIM ) THEN
       CALL Init_Global_Ch4( am_I_Root, Input_Opt, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Global_CH4"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Tagged CO
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_TAGCO_SIM   .or. &
         Input_Opt%ITS_A_H2HD_SIM  ) THEN
       CALL Init_Tagged_CO( am_I_Root, Input_Opt, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Tagged_CO"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Tagged O3
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
       CALL Init_Tagged_O3( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Tagged_O3"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! POPs
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_POPS_SIM ) THEN
       CALL Init_POPs( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_POPs"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Mercury
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN 

       ! Main mercury module
       CALL Init_Mercury( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Mercury"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Land mercury module
       CALL Init_Land_Mercury( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Land_Mercury"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Mercury deposition module
       CALL Init_Depo_Mercury ( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Depo_Mercury"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Ocean mercury module
       CALL Init_Ocean_Mercury( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Ocean_Mercury"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

#if defined( TOMAS )
    !-----------------------------------------------------------------
    ! TOMAS
    !-----------------------------------------------------------------

    ! Initialize the TOMAS microphysics package, if necessary
    CALL Init_Tomas_Microphysics( am_I_Root, Input_Opt, State_Chm, RC )    
#endif

    !-----------------------------------------------------------------
    ! Initialize tpcore_bc for nested grid simulations
    !-----------------------------------------------------------------
    CALL INIT_TPCORE_BC( am_I_Root, Input_Opt, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Tpcore_BC"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Call setup routines for certain GEOS-Chem diagnostics
    ! This allows us to use the species database object.
    !=================================================================

    ! Allocate and initialize variables
    CALL Ndxx_Setup( am_I_Root, Input_Opt, State_Chm, RC )

    ! Allocate diagnostic arrays
    CALL Init_Diag04
    CALL Init_Diag41
    CALL Init_Diag42( am_I_Root, Input_Opt, RC )
    IF ( Input_Opt%DO_ND48 ) CALL INIT_DIAG48 ( am_I_Root, Input_Opt, RC )
    IF ( Input_Opt%DO_ND49 ) CALL Init_Diag49 ( am_I_Root, Input_Opt, RC )
    IF ( Input_Opt%DO_ND50 ) CALL Init_Diag50 ( am_I_Root, Input_Opt, RC )
    IF ( Input_Opt%DO_ND51 ) CALL Init_Diag51 ( am_I_Root, Input_Opt, RC )
    IF ( Input_Opt%DO_ND51b) CALL Init_Diag51b( am_I_Root, Input_Opt, RC )
    CALL Init_Diag53
    CALL Init_Diag56
    IF ( Input_Opt%DO_ND63 ) CALL Init_Diag63 ( am_I_Root, Input_Opt, RC )

    ! Initialize the Hg diagnostics (bpch)
    CALL Init_Diag03( State_Chm )

    ! Initialize ND20 for tagged O3 simulation
    IF ( Input_Opt%DO_SAVE_O3 ) THEN
       CALL Init_Diag20( am_I_Root, Input_Opt, RC )
    ENDIF

    ! Enable Mean OH (or CH3CCl3) diag for runs which need it
    CALL Init_Diag_OH( am_I_Root, Input_Opt, RC )

#if !defined( ESMF_ )
    !--------------------------------------------------------------------
    ! Write out diaginfo.dat, tracerinfo.dat files for this simulation
    !
    ! NOTE: Do not do this for GCHP, because this will cause a file to
    ! be written out to disk for each core.  
    !
    ! ALSO NOTE: Eventually we will remove the ESMF_ C-preprocessor
    ! but we have to keep it for the time being (bmy, 4/11/18)
    !--------------------------------------------------------------------
    CALL Do_Gamap( am_I_Root, Input_Opt, State_Chm, RC )
#endif

    IF ( prtDebug ) CALL DEBUG_MSG( '### a GC_INIT_EXTRA' )

  END SUBROUTINE GC_Init_Extra
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_init_regridding
!
! !DESCRIPTION: Internal subroutine GC\_INIT\_REGRIDDING passes several
!  variables to regrid\_a2a\_mod.F90, where they are locally shadowed.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_Init_Regridding( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod
    USE ErrCode_Mod  
    USE Input_Opt_Mod,      ONLY : OptInput
    USE GC_Grid_Mod
    USE Regrid_A2A_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REMARKS:
!  This routine is a wrapper for Init_Map_A2A in regrid_a2a_mod.F90.  
!  Passing variables via Init_Map_A2A helps us to remove dependencies on
!  other GEOS-Chem routines from regrid_a2a_mod.F90.  This in turn 
!  facilitates the implementation of the HEMCO emissions package.
!
! !REVISION HISTORY: 
!  15 Jul 2014 - R. Yantosca - Initial version
!  05 Mar 2015 - R. Yantosca - Now read data w/r/t ExtData/HEMCO
!  26 Jan 2018 - M. Sulprizio- Moved to gc_environment_mod.F90 from main.F and
!                              renamed from Initialize_Regridding to
!                              GC_Init_Regridding
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!     
    ! Scalars
    INTEGER            :: I, J
    CHARACTER(LEN=255) :: DIR

    ! Arrays
    REAL(fp)           :: LONEDG(IIPAR+1       )  ! W longitude edges [deg]
    REAL(fp)           :: LATSIN(       JJPAR+1)  ! SIN(Lat edges)    [1  ]
    REAL(fp)           :: AREAS (IIPAR, JJPAR  )  ! Surface Areas     [m2 ]

    !================================================================
    ! GC_Init_Regridding begins here!
    !================================================================

    ! Assume success
    RC  = GC_SUCCESS

    ! Directory where netCDF ifles are found
    DIR = TRIM( Input_Opt%DATA_DIR ) // 'HEMCO/MAP_A2A/v2014-07/'

    ! Initialize longitudes [deg]
    DO I = 1, IIPAR+1
       LONEDG(I) = GET_XEDGE( I, 1, 1 )
    ENDDO

    ! Initialize sines of latitude [1]
    DO J = 1, JJPAR+1
       LATSIN(J) = GET_YSIN( 1, J, 1 )
    ENDDO

    ! Initialize surface areas [m2]
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       AREAS(I,J) = GET_AREA_M2( I, J, 1 )
    ENDDO
    ENDDO

    ! Pass to regrid_a2a_mod.F90, where these will be shadowed locally
    CALL Init_Map_A2A( IIPAR, JJPAR, LONEDG, LATSIN, AREAS, DIR )
      
  END SUBROUTINE GC_Init_Regridding
!EOC
#if defined( TOMAS )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tomas_microphysics
!
! !DESCRIPTION: Subroutine INIT\_TOMAS\_MICROPHYS will initialize the 
!  TOMAS microphysics package.  This replaces the former subroutine
!  READ\_MICROPHYSICS\_MENU.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_TOMAS_MICROPHYSICS( am_I_Root, Input_Opt, &
                                      State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE TOMAS_MOD,          ONLY : INIT_TOMAS
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REMARKS:
!  We now invoke TOMAS by compiling GEOS-Chem and setting either the TOMAS=yes
!  (30 bins, default) or TOMAS40=yes (40 bins, optional) switches.  The old
!  LTOMAS logical switch is now obsolete because all of the TOMAS code is
!  segregated from the rest of GEOS-Chem with #ifdef blocks.  Therefore,
!  we no longer need to read the microphysics menu, but we still need to
!  apply some error checks and then call INIT_TOMAS. (bmy, 4/23/13)
!                                                                             .
!  The Ind_() function now defines all species ID's.  It returns -1 if
!  a species cannot be found.  The prior behavior was to return 0 if a 
!  species wasn't found.  Therefore, in order to preserve the logic of the
!  error checks, we must force any -1's returned by Ind_() to 0's in
!  this subroutine.
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now read LNEI99 -- switch for EPA/NEI99 emissions (bmy, 11/5/04)
!  (2 ) Now read LAVHRRLAI-switch for using AVHRR-derived LAI (bmy, 12/20/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now read LMEGAN -- switch for MEGAN biogenics (tmf, bmy, 10/20/05)
!  (5 ) Now read LEMEP -- switch for EMEP emissions (bdf, bmy, 11/1/05)
!  (6 ) Now read LGFED2BB -- switch for GFED2 biomass emissions (bmy, 4/5/06)
!  (7 ) Now read LOTDLIS, LCTH, LMFLUX, LPRECON for lightning options 
!        (bmy, 5/10/06)
!  (8 ) Now read LBRAVO for BRAVO Mexican emissions (rjp, kfb, bmy, 6/26/06)
!  (9 ) Now read LEDGAR for EDGAR emissions (avd, bmy, 7/11/06)
!  (10) Now read LSTREETS for David Streets' emissions (bmy, 8/17/06)
!  (11) Kludge: Reset LMFLUX or LPRECON to LCTH, as the MFLUX and PRECON
!        lightning schemes have not yet been implemented.  Rename LOTDLIS
!        to LOTDREG.  Also read LOTDLOC for the OTD-LIS local redistribution
!        of lightning flashes (cf B. Sauvage).  Make sure LOTDREG and 
!        LOTDLOC are not both turned on at the same time. (bmy, 1/31/07)
!  (12) Add LOTDSCALE to the list of LNOx options (ltm, bmy, 9/24/07)
!  (13) Add new error traps for OTD-LIS options, dependent on met field type
!        (ltm, bmy, 11/29/07)
!  (14) Bug fix, create string variables for ERROR_STOP (bmy, 1/24/08)
!  (15) Now read LCAC for CAC emissions (amv, 1/09/2008)
!  (16) Now read LEDGARSHIP, LARCSHIP and LEMEPSHIP (phs, 12/5/08)
!  (17) Fixed typo in message for GEOS-3 (bmy, 10/30/08)
!  (18) Now read LVISTAS (amv, 12/2/08)
!  (19) Now read L8DAYBB, L3HRBB and LSYNOPBB for GFED2 8-days and 3hr
!        emissions, and LICARTT for corrected EPA (phs, yc, 12/17/08)
!  (20) Add a specific switch for MEGAN emissions for monoterpenes and MBO
!       (ccc, 2/2/09)
!  (21) Now read LICOADSSHIP (cklee, 6/30/09)
!  (22) Bug fix: for now, if LEMEPSHIP is turned on but LEMEP is turned off,
!        just turn off LEMEPSHIP and print a warning msg. (mak, bmy, 10/18/09)
!  (23) Now accounts for NEI2005 (amv, phs, 10/9/09)
!  (24) Included optional flag for using MODIS LAI data (mpb,2009).
!  (25) Included optional flag for using PCEEA model (mpb, 2009)
!  (26) Now force settings for EU, NA, CC nested grids (amv, bmy, 12/18/09)
!  (27) Now force MEGAN to use MODIS LAI (ccarouge, bmy, 2/24/10)
!  (28) Add separate switch for NOx fertilizer. (fp, 2/29/10)
!  (29) Add scaling for isoprene and NOx emissions. (fp, 2/29/10)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  23 Apr 2013 - R. Yantosca - Renamed to INIT_TOMAS_MICROPHYS
!  30 Jan 2014 - R. Yantosca - INIT_TOMAS accepts am_I_Root, Input_Opt, RC
!  16 Jun 2016 - K. Travis   - Now define species ID's with the Ind_ function 
!  16 Jun 2016 - E. Lundgren - INIT_TOMAS now accepts State_Chm
!  22 Jun 2016 - R. Yantosca - Force -1's returned by Ind_() to zeroes,
!                              in order to preserve the program logic
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, I

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_TOMAS_MICROPHYSICS begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = &
       ' -> at Init_Tomas_Microphysics (in GeosCore/input_mod.F)'

    ! Halt with error if we are trying to run TOMAS in a simulation
    ! that does not have any defined aerosols
    ! Turn off switches for simulations that don't use microphysics
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )  .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       ErrMsg = 'TOMAS needs to run with either a full-chemistry ' // &
                'or offline aerosol simulation!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Halt with error if none of the TOMAS aerosol species are defined
    I = MAX( Ind_('NK1'  ,'A'), 0 ) &
      + MAX( Ind_('SF1'  ,'A'), 0 ) &
      + MAX( Ind_('SS1'  ,'A'), 0 ) &
      + MAX( Ind_('ECOB1','A'), 0 ) &
      + MAX( Ind_('ECIL1','A'), 0 ) &
      + MAX( Ind_('OCOB1','A'), 0 ) &
      + MAX( Ind_('OCIL1','A'), 0 ) &
      + MAX( Ind_('DUST1','A'), 0 )

    IF ( I == 0 ) THEN 
       ErrMsg = 'None of the TOMAS aerosols are defined!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
      
    ! Halt with error if sulfate aerosols are not defined
    IF( Ind_('SF1') > 0 .and. ( .not. Input_Opt%LSULF ) ) THEN
       ErrMsg = 'Need LSULF on for TOMAS-Sulfate to work (for now)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Halt with error if carbonaceous aerosols are not defined
    I = MAX( Ind_('ECOB1','A'), 0 ) &
      + MAX( Ind_('ECIL1','A'), 0 ) & 
      + MAX( Ind_('OCOB1','A'), 0 ) &
      + MAX( Ind_('OCIL1','A'), 0 )

    IF ( I > 0 .and. (.not. Input_Opt%LCARB ) ) THEN
       ErrMsg = 'Need LCARB on for TOMAS-carb to work (for now)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Halt with error if dust aerosols are turned on.
    ! TOMAS defines its own dust aerosol species.
!    IF ( Ind_('DUST1') > 0 .AND. Input_Opt%LDUST ) THEN
!       MSG = 'Need to turn off LDUST for TOMAS Dust to work'
!       CALL ERROR_STOP( MSG, LOCATION )
!    ENDIF

    !=================================================================
    ! All error checks are satisfied, so initialize TOMAS
    !=================================================================
    CALL INIT_TOMAS( am_I_Root, Input_Opt, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_TOMAS"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE INIT_TOMAS_MICROPHYSICS
!EOC
#endif
END MODULE GC_Environment_Mod
