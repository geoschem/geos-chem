!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: depo_mercury_mod.F90
!
! !DESCRIPTION: Module DEPO\_MERCURY\_MOD contains routines to handle
!  deposition fluxes for mercury.
!
! !INTERFACE:
!
MODULE DEPO_MERCURY_MOD
!
! !USES:
!
  USE PRECISION_MOD   ! For Geos-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: ADD_Hg2_DD
  PUBLIC :: ADD_Hg2_WD
  PUBLIC :: ADD_HgP_DD
  PUBLIC :: ADD_HgP_WD
  PUBLIC :: ADD_HG2_SNOWPACK
  PUBLIC :: RESET_HG_DEP_ARRAYS
  PUBLIC :: CHECK_DIMENSIONS
#ifdef BPCH_DIAG
  PUBLIC :: READ_GTMM_RESTART
  PUBLIC :: MAKE_GTMM_RESTART
  PUBLIC :: UPDATE_DEP
#endif
  PUBLIC :: INIT_DEPO_MERCURY
  PUBLIC :: CLEANUP_DEPO_MERCURY
!
! !PUBLIC DATA MEMBERS:
!
  PUBLIC :: DD_HG2, DD_HGP, WD_HG2, WD_HGP
  PUBLIC :: HG2mth_wd, HG0mth_dd, HG2mth_dd
  PUBLIC :: LHGSNOW

  REAL(fp),  ALLOCATABLE :: DD_Hg2(:,:,:)
  REAL(fp),  ALLOCATABLE :: DD_HgP(:,:,:)
  REAL(fp),  ALLOCATABLE :: WD_Hg2(:,:,:)
  REAL(fp),  ALLOCATABLE :: WD_HgP(:,:,:)
  REAL(fp),  ALLOCATABLE :: HG0mth_dd(:,:)
  REAL(fp),  ALLOCATABLE :: HG2mth_dd(:,:)
  REAL(fp),  ALLOCATABLE :: HG2mth_wd(:,:)
  REAL(fp),  ALLOCATABLE :: Hg0dryGEOS(:,:)
  REAL(fp),  ALLOCATABLE :: HgIIdryGEOS(:,:)
  REAL(fp),  ALLOCATABLE :: HgIIwetGEOS(:,:)
!
! !PRIVATE DATA MEMBERS:
!
  CHARACTER(LEN=255)     :: GTMM_RST_FILE
  LOGICAL                :: LHGSNOW
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Scalars for Hg indexing
  INTEGER          :: N_Hg_CATS
  INTEGER          :: id_Hg0
  INTEGER          :: id_Hg2

  ! Pointers for Hg Indexing
  ! NOTE: Because these are SAVEd pointers (by virtue of being
  ! module variables) we can set these to NULL here. (bmy, 4/29/16)
  INTEGER, POINTER :: Hg0_Id_List (:) => NULL()
  INTEGER, POINTER :: Hg2_Id_List (:) => NULL()
  INTEGER, POINTER :: HgP_Id_List (:) => NULL()
  INTEGER, POINTER :: Hg0_CAT(:) => NULL()
  INTEGER, POINTER :: Hg2_CAT(:) => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_Hg2_dd
!
! !DESCRIPTION: Subroutine ADD\_Hg2\_DD computes the amount of Hg(II) dry
!  deposited out of the atmosphere into the column array DD\_Hg2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ADD_Hg2_DD( I, J, HG_CAT, DRY_Hg2 )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I, J      ! Grid box lon & lat indices
    INTEGER,  INTENT(IN) :: Hg_Cat    ! Hg category number
    REAL(fp), INTENT(IN) :: DRY_Hg2   ! Hg(II) dry deposited out of the
                                      !  atmosphere [kg]
!
! !REVISION HISTORY:
!  19 Jan 2005 - S. Strode, C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Store dry deposited Hg(II) into DD_Hg2 array
    IF ( Hg_Cat > 0 ) THEN
       DD_Hg2(I,J,Hg_Cat) = DD_Hg2(I,J,Hg_Cat) + DRY_Hg2
    ENDIF

  END SUBROUTINE ADD_Hg2_DD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_Hg2_wd
!
! !DESCRIPTION: Subroutine ADD\_Hg2\_WD computes the amount of Hg(II) wet
!  scavenged out of the atmosphere into the column array WD\_Hg2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ADD_Hg2_WD( I, J, Hg_Cat, WET_Hg2 )
!
! !INPUT PARAMETERS:
!
    INTEGER,   INTENT(IN) :: I, J      ! Grid box lon & lat indices
    INTEGER,   INTENT(IN) :: Hg_Cat    ! Hg category number
    REAL(fp),  INTENT(IN) :: WET_Hg2   ! Hg(II) scavenged out of the
                                       !  atmosphere [kg]
!
! !REVISION HISTORY:
!  19 Jan 2005 - S. Strode, C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Store wet deposited Hg(II) into WD_Hg2 array
    IF ( Hg_Cat > 0 ) THEN
       WD_Hg2(I,J,Hg_Cat) = WD_Hg2(I,J,Hg_Cat) + WET_Hg2
    ENDIF

  END SUBROUTINE ADD_Hg2_WD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_HgP_dd
!
! !DESCRIPTION: Subroutine ADD\_HgP\_DD computes the amount of HgP dry
!  deposited out of the atmosphere into the column array DD\_HgP.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ADD_HgP_DD( I, J, Hg_Cat, DRY_HgP )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I, J      ! Grid box lon & lat indices
    INTEGER,  INTENT(IN) :: Hg_Cat    ! Hg category number
    REAL(fp), INTENT(IN) :: DRY_HgP   ! HgP dry deposited out of the
                                      !  atmosphere [kg]
!
! !REVISION HISTORY:
!  19 Jan 2005 - S. Strode, C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Store dry deposited Hg(II) into DD_Hg2 array
    IF ( Hg_Cat > 0 ) THEN
       DD_HgP(I,J,Hg_Cat) = DD_HgP(I,J,Hg_Cat) + DRY_HgP
    ENDIF

  END SUBROUTINE ADD_HgP_DD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_HgP_wd
!
! !DESCRIPTION: Subroutine ADD\_HgP\_WD computes the amount of HgP wet
!  scavenged out of the atmosphere into the column array WD\_HgP.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ADD_HgP_WD( I, J, Hg_Cat, WET_HgP )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I, J      ! Grid box lon & lat indices
    INTEGER,  INTENT(IN) :: Hg_Cat    ! Hg category number
    REAL(fp), INTENT(IN) :: WET_HgP   ! HgP scavenged out of the
                                      !  atmosphere [kg]
!
! !REVISION HISTORY:
!  19 Jan 2005 - S. Strode, C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Store wet deposited HgP into WD_HgP array
    IF ( Hg_Cat > 0 ) THEN
       WD_HgP(I,J,Hg_Cat) = WD_HgP(I,J,Hg_Cat) + WET_HgP
    ENDIF

  END SUBROUTINE ADD_HgP_WD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_hg2_snowpack
!
! !DESCRIPTION: Subroutine ADD\_Hg2\_SNOWPACKS adds Hg2 deposition to snowpack.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ADD_HG2_SNOWPACK( I, J, Hg_Cat, DEP_Hg2, &
                               State_Met, State_Chm, State_Diag )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03, ND03
#endif
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Met_Mod,      ONLY : MetState
    USE Time_Mod,           ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J        ! Grid box lon & lat indices
    INTEGER,        INTENT(IN)    :: Hg_Cat      ! Hg category number
    REAL(fp),       INTENT(IN)    :: Dep_Hg2     ! Hg2 (or HgP) deposited
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry   State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !REVISION HISTORY:
!  02 Sep 2008 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: NN
    LOGICAL  :: IS_SNOW_OR_ICE
    REAL(fp) :: FRAC_SNOW_OR_ICE
    REAL(fp) :: FRAC_O
    REAL(fp) :: DT

    ! Pointers
    REAL(fp), POINTER :: SNOW_HG_OC(:,:,:)
    REAL(fp), POINTER :: SNOW_HG_LN(:,:,:)
    REAL(fp), POINTER :: SNOW_HG_STORED_OC(:,:,:)
    REAL(fp), POINTER :: SNOW_HG_STORED_LN(:,:,:)

    !=================================================================
    ! ADD_HG2_SNOWPACK begins here!
    !=================================================================

    ! Point to fields in State_Chm
    SNOW_HG_OC        => State_Chm%SnowHgOcean
    SNOW_HG_LN        => State_Chm%SnowHgLand
    SNOW_HG_STORED_OC => State_Chm%SnowHgOceanStored
    SNOW_HG_STORED_LN => State_Chm%SnowHgLandStored

    ! Store the Hg category number in NN
    NN = Hg_Cat

    ! Return if snowpack model is disabled
    IF ( .NOT. LHGSNOW ) RETURN

    ! Chemistry timestep [s]
    DT = GET_TS_CHEM()

    ! Don't let fraction be greater than 1
    FRAC_SNOW_OR_ICE = MIN( State_Met%FRSNO(I,J)     + &
                            State_Met%FRSEAICE(I,J)  + &
                            State_Met%FRLANDIC(I,J), 1e+0_fp)
    IS_SNOW_OR_ICE   = ( FRAC_SNOW_OR_ICE > 0e+0_fp )

    ! Ocean fraction (vs land)
    FRAC_O = State_Met%FROCEAN(I,J)

    ! Check if there is snow on the ground, or if this is sea ice
    IF ( IS_SNOW_OR_ICE ) THEN

       ! Add 60% of deposition to surface (i.e. assume 40% accumulates
       ! in surface) multiplied by the fraction of the box that is
       ! snow or ice (i.e. 1 for all met fields but MERRA and GEOS-5.7.x)
       SNOW_HG_OC(I,J,NN) = SNOW_HG_OC(I,J,NN) + &
                            FRAC_O * FRAC_SNOW_OR_ICE * &
                            MAX( 0.6e+0_fp*DEP_HG2, 0e+0_fp )

       ! Add remaining deposited Hg to reservoir for delivery to ocean
       ! when snow melts later (jaf, 6/17/11)
       ! This is Hg in snowpack over ocean that CANNOT be
       ! re-emitted to atmosphere
       SNOW_HG_STORED_OC(I,J,NN) = SNOW_HG_STORED_OC(I,J,NN) + &
                                   FRAC_O * FRAC_SNOW_OR_ICE * &
                                   MAX( 0.4e+0_fp*DEP_HG2, 0e+0_fp )

       ! This is Hg in snowpack over land that can potentially be
       ! re-emitted to atmosphere
       SNOW_HG_LN(I,J,NN) = SNOW_HG_LN(I,J,NN) + &
                            (1e+0_fp - FRAC_O) * FRAC_SNOW_OR_ICE * &
                            MAX( 0.6e+0_fp*DEP_HG2, 0e+0_fp )

       ! This is Hg in snowpack over land that CANNOT be
       ! re-emitted to atmosphere
       SNOW_HG_STORED_LN(I,J,NN) = SNOW_HG_STORED_LN(I,J,NN) + &
                                   ( 1e+0_fp - FRAC_O ) * &
                                   FRAC_SNOW_OR_ICE &
                                   * MAX( 0.4e+0_fp*DEP_HG2, 0e+0_fp )

#ifdef BPCH_DIAG
       ! Store diagnostic of TOTAL HgII/HgP deposition to snow/ice
       IF ( ND03 > 0 ) AD03(I,J,21,1) = AD03(I,J,21,1)   + &
                                        FRAC_SNOW_OR_ICE * &
                                        MAX(DEP_HG2, 0e+0_fp)
#endif

       !--------------------------------------------------------------
       ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
       !
       ! Store diagnostic of TOTAL HgII/HgP deposition to snow/ice
       ! NOTE: Units are now kg/s
       !--------------------------------------------------------------
       IF ( State_Diag%Archive_FluxHg2HgPfromAirToSnow ) THEN
          State_Diag%FluxHg2HgPfromAirToSnow(I,J) = &
             State_Diag%FluxHg2HgPfromAirToSnow(I,J) + &
             ( FRAC_SNOW_OR_ICE * MAX( DEP_HG2, 0e+0_fp ) ) !/ DT
       ENDIF

    ENDIF

    ! Free pointers
    SNOW_HG_OC        => NULL()
    SNOW_HG_LN        => NULL()
    SNOW_HG_STORED_OC => NULL()
    SNOW_HG_STORED_LN => NULL()

  END SUBROUTINE ADD_HG2_SNOWPACK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: reset_hg_dep_arrays
!
! !DESCRIPTION: Subroutine RESET\_Hg\_DEP\_ARRAYS resets the wet and dry
!  deposition arrays for Hg(II) and Hg(p) to zero. This allows us to call
!  OCEAN\_MERCURY\_FLUX and LAND\_MERCURY\_FLUX in any order in MERCURY\_MOD.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RESET_HG_DEP_ARRAYS
!
! !REVISION HISTORY:
!  02 Sep 2008 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Reset deposition arrays.
    DD_Hg2 = 0e+0_fp
    WD_Hg2 = 0e+0_fp
    DD_HgP = 0e+0_fp
    WD_HgP = 0e+0_fp

  END SUBROUTINE RESET_HG_DEP_ARRAYS
!EOC
#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: make_gtmm_restart
!
! !DESCRIPTION: MAKE\_GTMM\_RESTART writes a GTMM restart file with deposition
!  fluxes and store deposition fluxes for continuous runs.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MAKE_GTMM_RESTART( Input_Opt, State_Grid, NYMD, NHMS, TAU, RC )
!
! !USES:
!
    USE BPCH2_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,       ONLY : OptInput
    USE inquireMod,          ONLY : findFreeLUN
    USE TIME_MOD,            ONLY : EXPAND_DATE
    USE TIME_MOD,            ONLY : GET_CT_DYN, GET_CT_CHEM
    USE State_Grid_Mod,      ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)   :: NYMD        ! Year-Month-Date
    INTEGER,        INTENT(IN)   :: NHMS        ! and Hour-Min-Sec
    REAL(f8),       INTENT(IN)   :: TAU         ! TAU value: hrs since 1/1/85
    TYPE(OptInput), INTENT(IN)   :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)   :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)  :: RC          ! Success or failure?

! !REVISION HISTORY:
!  15 Sep 2009 - C. Carouge  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: HALFPOLAR, CENTER180, IU_FILE
    INTEGER               :: IFIRST,    JFIRST,    LFIRST
    INTEGER               :: N,         NN
    REAL(f8)              :: TS_DYN,    TS_CHEM
    REAL(f4)              :: LONRES,    LATRES
    REAL(f4)              :: ARRAY(State_Grid%NX,State_Grid%NY,1)
    CHARACTER(LEN=20)     :: MODELNAME
    CHARACTER(LEN=40)     :: CATEGORY,  UNIT,     RESERVED
    CHARACTER(LEN=255)    :: FILENAME

    !=================================================================
    ! MAKE_GTMM_RESTART begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Initialize values
    IFIRST    = State_Grid%XMinOffset + 1
    JFIRST    = State_Grid%YMinOffset + 1
    LFIRST    = 1
    HALFPOLAR = GET_HALFPOLAR()
    CENTER180 = 1
    LONRES    = State_Grid%DX
    LATRES    = State_Grid%DY
    MODELNAME = GET_MODELNAME( Input_Opt, State_Grid )
    CATEGORY  = 'DRYD-FLX'
    RESERVED  = ''
    UNIT      = 'molec/cm2/s'

    ! Find a free file LUN
    IU_FILE   = findFreeLUN()

    ! Expand date in filename
    FILENAME  = TRIM( Input_Opt%RUN_DIR ) // TRIM( GTMM_RST_FILE )
    CALL EXPAND_DATE( FILENAME, NYMD, NHMS )
    
    ! Echo info
    WRITE( 6, 100 ) TRIM( FILENAME ), IU_FILE
100 FORMAT( '     - MAKE_RESTART_FILE: Writing ', a, ' on unit ', i4 )

    ! Open BPCH file for output
    CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME )

    !---------------------------
    ! Total Hg(0) dry deposition
    !---------------------------
    ARRAY(:,:,1) = HG0mth_dd

    CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,   &
                HALFPOLAR, CENTER180, CATEGORY, N,        &
                UNIT,      TAU,       TAU,      RESERVED, &
                State_Grid%NX, State_Grid%NY, 1, IFIRST,  &
                JFIRST,    LFIRST,    ARRAY(:,:,1) )

    !---------------------------
    ! Hg(II) dry deposition
    !---------------------------
    ARRAY(:,:,1) = HG2mth_dd

    CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,   &
                HALFPOLAR, CENTER180, CATEGORY, N,        &
                UNIT,      TAU,       TAU,      RESERVED, &
                State_Grid%NX, State_Grid%NY, 1, IFIRST,  &
                JFIRST,    LFIRST,    ARRAY(:,:,1) )

    !---------------------------
    ! Hg(II) wet deposition
    !---------------------------
    CATEGORY  = 'WETDLS-$'
    ARRAY(:,:,1) = HG2mth_wd

    CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,   &
                HALFPOLAR, CENTER180, CATEGORY, N,        &
                UNIT,      TAU,       TAU,      RESERVED, &
                State_Grid%NX, State_Grid%NY, 1, IFIRST,  &
                JFIRST,    LFIRST,    ARRAY(:,:,1) )

    ! Close file
    CLOSE( IU_FILE )

  END SUBROUTINE MAKE_GTMM_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gtmm_restart
!
! !DESCRIPTION: Subroutine READ\_GTMM\_RESTART reads dry and wet deposition
!  for mercury from GTMM restart.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_GTMM_RESTART( Input_Opt, State_Chm, State_Grid, &
                                YYYYMMDD,    HHMMSS,              &
                                Hg0dryGEOS,  HgIIdryGEOS,         &
                                HgIIwetGEOS, RC           )
!
! !USES:
!
    USE BPCH2_MOD,            ONLY : OPEN_BPCH2_FOR_READ
    USE ErrCode_Mod
    USE ERROR_MOD,            ONLY : DEBUG_MSG
    USE FILE_MOD,             ONLY : IOERROR
    USE Input_Opt_Mod,        ONLY : OptInput
    USE inquireMod,           ONLY : findFreeLun
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Grid_Mod,       ONLY : GrdState
    USE TIME_MOD,             ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    INTEGER,        INTENT(IN)    :: YYYYMMDD, HHMMSS
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT)   :: RC          ! Success or failure?!
    REAL(fp), DIMENSION(State_Grid%NX, State_Grid%NY)   :: Hg0dryGEOS
    REAL(fp), DIMENSION(State_Grid%NX, State_Grid%NY)   :: HgIIdryGEOS
    REAL(fp), DIMENSION(State_Grid%NX, State_Grid%NY)   :: HgIIwetGEOS
!
! !REMAKRS:
!  NOTE: THIS ROUTINE WILL HAVE TO BE UPDATED TO CONVERT TO NETCDF FORMAT!!!
!
! !REVISION HISTORY:
!  15 Sep 2009 - C. Carouge  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: IOS, I, IU_FILE, J, L, NN, N_gtmm
    INTEGER               :: NCOUNT(State_Chm%nAdvect)
    REAL*4                :: FLUX(State_Grid%NX,State_Grid%NY)
    CHARACTER(LEN=255)    :: FILENAME

    ! For binary punch file, version 2.0
    INTEGER               :: NI,        NJ,      NL
    INTEGER               :: IFIRST,    JFIRST,  LFIRST
    INTEGER               :: NTRACER,   NSKIP
    INTEGER               :: HALFPOLAR, CENTER180
    REAL*4                :: LONRES,    LATRES
    REAL(fp)              :: ZTAU0,     ZTAU1
    CHARACTER(LEN=20)     :: MODELNAME
    CHARACTER(LEN=40)     :: CATEGORY
    CHARACTER(LEN=40)     :: UNIT
    CHARACTER(LEN=40)     :: RESERVED

    !=================================================================
    ! READ_GTMM_RESTART begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Find a free file LUN
    IU_FILE  = findFreeLUN()

    ! Copy input file name to a local variable
    FILENAME = TRIM( Input_Opt%RUN_DIR ) // TRIM( GTMM_RST_FILE )

    ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
    CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

    ! Echo some input to the screen
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, 100   )
       WRITE( 6, 110   ) TRIM( FILENAME ), IU_FILE
100    FORMAT( 'G T M M  H g   R E S T A R T   F I L E   I N P U T' )
110    FORMAT( /, 'READ_GTMM_RESTART: Reading ', a, ' on unit ', i6 )
    ENDIF

    ! Open the binary punch file for input
    CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

    !=================================================================
    ! Read concentrations -- store in the TRACER array
    !=================================================================
    DO
       READ( IU_FILE, IOSTAT=IOS ) &
            MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

       ! IOS < 0 is end-of-file, so exit
       IF ( IOS < 0 ) EXIT

       ! IOS > 0 is a real I/O error -- print error message
       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_rst:1' )

       READ( IU_FILE, IOSTAT=IOS ) &
            CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
            NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,   &
            NSKIP

       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_rst:2' )

       READ( IU_FILE, IOSTAT=IOS ) &
            ( ( FLUX(I,J), I=1,NI ), J=1,NJ )

       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_rst:3' )

       !==============================================================
       ! Assign data from the TRACER array to the STT array.
       !==============================================================

       ! Process dry deposition data
       IF ( CATEGORY(1:8) == 'DRYD-FLX' ) THEN

          ! Make sure array dimensions are of global size
          ! (NI=State_Grid%NX; NJ=State_Grid%NY, NL=State_Grid%NZ),
          ! or stop the run
          CALL CHECK_DIMENSIONS( NI, NJ, NL, State_Grid )

          ! Save into arrays
          IF ( ANY( Hg0_Id_List == NTRACER ) ) THEN

             !----------
             ! Hg(0)
             !----------

             ! Get the Hg category #
             NN              = Hg0_Cat(NTRACER)

             ! Store ocean Hg(0) in Hg0aq array
             Hg0dryGEOS(:,:) = FLUX(:,:)

             ! Increment NCOUNT
             NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1

          ELSE IF ( ANY( Hg2_Id_List == NTRACER ) ) THEN

             !----------
             ! Hg(II)
             !----------

             ! Get the Hg category #
             NN               = Hg2_Cat(NTRACER)

             ! Store ocean Hg(II) in Hg2_aq array
             HgIIdryGEOS(:,:) = FLUX(:,:)

             ! Increment NCOUNT
             NCOUNT(NTRACER)  = NCOUNT(NTRACER) + 1

          ENDIF
       ENDIF

       ! Process wet deposition data
       IF ( CATEGORY(1:8) == 'WETDLS-$' ) THEN

          ! Make sure array dimensions are of global size
          ! (NI=State_Grid%NX; NJ=State_Grid%NY, NL=State_Grid%NZ),
          ! or stop the run
          ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          ! %%%% KLUDGE: CHECK_DIMENSIONS only works for NL=1 !!!!
          ! And we are only interested by the surface flux...
          NL = 1
          ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          CALL CHECK_DIMENSIONS( NI, NJ, NL, State_Grid )

          IF ( ANY( Hg2_Id_List == NTRACER ) ) THEN

             !----------
             ! Hg(II)
             !----------

             ! Get the Hg category #
             NN               = Hg2_Cat(NTRACER)

             ! Store ocean Hg(II) in Hg2_aq array
             HgIIwetGEOS(:,:) = FLUX(:,:)

             ! Increment NCOUNT
             NCOUNT(NTRACER)  = NCOUNT(NTRACER) + 1

          ENDIF
       ENDIF
    ENDDO

    ! Close file
    CLOSE( IU_FILE )

  END SUBROUTINE READ_GTMM_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_dep
!
! !DESCRIPTION: Subroutine UPDATE\_DEP update the monthly average for wet and
!  dry deposition of Hg0 and Hg2 for mercury from GTMM restart.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE UPDATE_DEP( NN )
!
! !USES:
!
    USE TIME_MOD,     ONLY : GET_CT_DYN,  GET_CT_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER :: NN    ! Hg2 ID for wet deposition
!
! !REVISION HISTORY:
!  04 June 2010  - C. Carouge  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: N
    REAL(fp) :: SCALEDYN, SCALECHEM

    !=================================================================
    ! UPDATE_DEP begins here!
    !=================================================================

    ! counter variables
    SCALEDYN   = DBLE( GET_CT_DYN()  ) + 1e-32_fp
    SCALECHEM  = DBLE( GET_CT_CHEM() ) + 1e-32_fp

    !! Hg2 total wet deposition at the surface
    !HG2mth_wd = HG2mth_wd + ( SUM(AD38(:,:,:,NN), DIM=3) + &
    !                          SUM(AD39(:,:,:,NN), DIM=3) ) / SCALEDYN
    !
    !! Hg0 total dry deposition at the surface
    !N         = id_Hg0
    !HG0mth_dd = HG0mth_dd + AD44(:,:,N,1) / SCALECHEM
    !
    !! Hg2 total dry deposition at the surface
    !N         = id_Hg2
    !HG2mth_dd = HG2mth_dd + AD44(:,:,N,1) / SCALECHEM

  END SUBROUTINE UPDATE_DEP
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dimensions
!
! !DESCRIPTION: Subroutine CHECK\_DIMENSIONS makes sure that the dimensions of
!  the Hg restart file extend to cover the entire grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_DIMENSIONS( NI, NJ, NL, State_Grid )
!
! !USES:
!
    USE ERROR_MOD,      ONLY : GEOS_CHEM_STOP
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: NI, NJ, NL
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  30 Aug 2010 - S. Strode, C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! CHECK_DIMENSIONS begins here!
    !=================================================================

    ! Error check longitude dimension: NI must equal State_Grid%NX
    IF ( NI /= State_Grid%NX ) THEN
       WRITE( 6, 100 )
100    FORMAT( 'ERROR reading in Hg restart file', / &
               'Wrong number of longitudes encountered', / &
               'STOP in CHECK_DIMENSIONS ("depo_mercury_mod.F90")' )
       PRINT*, "NI: ", NI
       PRINT*, "State_Grid%NX: ", State_Grid%NX
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       CALL GEOS_CHEM_STOP
    ENDIF

    ! Error check latitude dimension: NJ must equal State_Grid%NY
    IF ( NJ /= State_Grid%NY ) THEN
       WRITE( 6, 110 )
110    FORMAT( 'ERROR reading in Hg restart file', / &
               'Wrong number of latitudes encountered', / &
               'STOP in CHECK_DIMENSIONS ("depo_mercury_mod.F90")' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       CALL GEOS_CHEM_STOP
    ENDIF

    ! Error check vertical dimension: NL must equal State_Grid%NZ
    IF ( NL /= 1 ) THEN
       WRITE( 6, 120 )
120    FORMAT( 'ERROR reading in Hg restart file', / &
               'Wrong number of vertical encountered', / &
               'STOP in CHECK_DIMENSIONS ("depo_mercury_mod.F90")' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       CALL GEOS_CHEM_STOP
    ENDIF

  END SUBROUTINE CHECK_DIMENSIONS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_depo_mercury
!
! !DESCRIPTION: Subroutine INIT\_DEPO\_MERCURY initialize deposition arrays
!  for mercury.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_DEPO_MERCURY( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
!  23 Apr 2010 - C. Carouge  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N, C

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! INIT_DEPO_MERCURY begins here!
    !=================================================================

    ! GTMM restart file name
    GTMM_RST_FILE = Input_Opt%GTMM_RST_FILE

    !================================================================
    ! Initialize local Hg indexing variables
    !================================================================

    ! Now get # of tagHg categories from State_Chm
    N_Hg_CATS     = State_Chm%N_Hg_CATS

    ! Hg species index corresponding to a given Hg category number
    Hg0_Id_List   => State_Chm%Hg0_Id_List
    Hg2_Id_List   => State_Chm%Hg2_Id_List
    HgP_Id_List   => State_Chm%HgP_Id_List

    ! Category numbers for each Hg0 tracer
    ALLOCATE( Hg0_Cat( State_Chm%nSpecies ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'Hg0_Cat' )
    Hg0_Cat = 0

    ! Category numbers for each Hg2 tracer
    ALLOCATE( Hg2_Cat( State_Chm%nSpecies ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'Hg0_Cat' )
    Hg2_Cat = 0

    ! Loop ovver all Hg species
    DO N = 1, State_Chm%nSpecies

       ! Point to the Species Database entry for tracer # N
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Define local id_Hg0 and id_Hg2 flags
       SELECT CASE( TRIM( SpcInfo%Name ) )
       CASE( 'Hg0' )
          id_Hg0 = SpcInfo%ModelId
       CASE( 'Hg2' )
          id_Hg2 = SpcInfo%ModelId
       CASE DEFAULT
          ! Do nothing
       END SELECT

       ! Store the Hg0 and Hg2 category numbers for each tracer
       ! for future lookup in READ_GTMM_RESTART (bmy, 4/26/16)
       IF ( SpcInfo%Is_Hg0 ) THEN
          Hg0_Cat(N) = SpcInfo%Hg_Cat
       ELSE iF ( SpcInfo%Is_Hg2 ) THEN
          Hg2_Cat(N) = SpcInfo%Hg_Cat
       ENDIF

       ! Free pointer
       SpcInfo => NULL()
    ENDDO

    !================================================================
    ! Allocate arrays
    !================================================================
    ALLOCATE( DD_Hg2( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'DD_Hg2' )
    DD_Hg2 = 0e+0_fp

    ALLOCATE( DD_HgP( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'DD_HgP' )
    DD_HgP = 0e+0_fp

    ALLOCATE( WD_Hg2( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'WD_Hg2' )
    WD_Hg2 = 0e+0_fp

    ALLOCATE( WD_HgP( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'WD_HgP' )
    WD_HgP = 0e+0_fp

    IF ( Input_Opt%LGTMM ) THEN
       ALLOCATE( HG0mth_dd( State_Grid%NX, State_Grid%NY ), STAT=RC )
       IF ( RC /= 0 ) CALL ALLOC_ERR( 'HG0mth_dd' )
       HG0mth_dd = 0e+0_fp

       ALLOCATE( HG2mth_dd( State_Grid%NX, State_Grid%NY ), STAT=RC )
       IF ( RC /= 0 ) CALL ALLOC_ERR( 'HG2mth_dd' )
       HG2mth_dd = 0e+0_fp

       ALLOCATE( HG2mth_wd( State_Grid%NX, State_Grid%NY ), STAT=RC )
       IF ( RC /= 0 ) CALL ALLOC_ERR( 'HG2mth_wd' )
       HG2mth_wd = 0e+0_fp
    ENDIF

  END SUBROUTINE INIT_DEPO_MERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_depo_mercury
!
! !DESCRIPTION: Subroutine CLEANUP\_DEPO\_MERCURY deallocate all arrays
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_DEPO_MERCURY
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
    IF ( ALLOCATED( DD_Hg2      ) ) DEALLOCATE( DD_Hg2      )
    IF ( ALLOCATED( DD_HgP      ) ) DEALLOCATE( DD_HgP      )
    IF ( ALLOCATED( WD_Hg2      ) ) DEALLOCATE( WD_Hg2      )
    IF ( ALLOCATED( WD_HgP      ) ) DEALLOCATE( WD_HgP      )
    IF ( ALLOCATED( HG0mth_dd   ) ) DEALLOCATE( HG0mth_dd   )
    IF ( ALLOCATED( HG2mth_dd   ) ) DEALLOCATE( HG2mth_dd   )
    IF ( ALLOCATED( HG2mth_wd   ) ) DEALLOCATE( HG2mth_wd   )

    ! Free pointers
    Hg0_Id_List => NULL()
    Hg2_Id_List => NULL()
    HgP_Id_List => NULL()

    IF ( ASSOCIATED( Hg0_Cat ) ) DEALLOCATE( Hg0_CAT )
    IF ( ASSOCIATED( Hg2_Cat ) ) DEALLOCATE( Hg2_CAT )

  END SUBROUTINE CLEANUP_DEPO_MERCURY
!EOC
END MODULE DEPO_MERCURY_MOD
