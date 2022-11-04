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
  PUBLIC :: INIT_DEPO_MERCURY
  PUBLIC :: CLEANUP_DEPO_MERCURY
!
! !PUBLIC DATA MEMBERS:
!
  PUBLIC :: DD_HG2, DD_HGP, WD_HG2, WD_HGP
  PUBLIC :: HG2mth_wd, HG0mth_dd, HG2mth_dd
  PUBLIC :: LHGSNOW

  REAL(fp),  ALLOCATABLE :: DD_Hg2(:,:)
  REAL(fp),  ALLOCATABLE :: DD_HgP(:,:)
  REAL(fp),  ALLOCATABLE :: WD_Hg2(:,:)
  REAL(fp),  ALLOCATABLE :: WD_HgP(:,:)
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
  INTEGER :: id_Hg0
  INTEGER :: id_Hg2

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
  SUBROUTINE ADD_Hg2_DD( I, J, DRY_Hg2 )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I, J      ! Grid box lon & lat indices
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
    DD_Hg2(I,J) = DD_Hg2(I,J) + DRY_Hg2

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
  SUBROUTINE ADD_Hg2_WD( I, J, WET_Hg2 )
!
! !INPUT PARAMETERS:
!
    INTEGER,   INTENT(IN) :: I, J      ! Grid box lon & lat indices
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
    WD_Hg2(I,J) = WD_Hg2(I,J) + WET_Hg2

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
  SUBROUTINE ADD_HgP_DD( I, J, DRY_HgP )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I, J      ! Grid box lon & lat indices
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
    DD_HgP(I,J) = DD_HgP(I,J) + DRY_HgP

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
  SUBROUTINE ADD_HgP_WD( I, J, WET_HgP )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I, J      ! Grid box lon & lat indices
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
    WD_HgP(I,J) = WD_HgP(I,J) + WET_HgP

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
  SUBROUTINE ADD_HG2_SNOWPACK( I, J, DEP_Hg2, State_Met, State_Chm, State_Diag )
!
! !USES:
!
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Met_Mod,      ONLY : MetState
    USE Time_Mod,           ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J        ! Grid box lon & lat indices
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
    LOGICAL  :: IS_SNOW_OR_ICE
    REAL(fp) :: FRAC_SNOW_OR_ICE
    REAL(fp) :: FRAC_O
    REAL(fp) :: DT

    ! Pointers
    REAL(fp), POINTER :: SNOW_HG_OC(:,:)
    REAL(fp), POINTER :: SNOW_HG_LN(:,:)
    REAL(fp), POINTER :: SNOW_HG_STORED_OC(:,:)
    REAL(fp), POINTER :: SNOW_HG_STORED_LN(:,:)

    !=================================================================
    ! ADD_HG2_SNOWPACK begins here!
    !=================================================================

    ! Point to fields in State_Chm
    SNOW_HG_OC        => State_Chm%SnowHgOcean
    SNOW_HG_LN        => State_Chm%SnowHgLand
    SNOW_HG_STORED_OC => State_Chm%SnowHgOceanStored
    SNOW_HG_STORED_LN => State_Chm%SnowHgLandStored

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
       SNOW_HG_OC(I,J) = SNOW_HG_OC(I,J) + &
                         FRAC_O * FRAC_SNOW_OR_ICE * &
                         MAX( 0.6e+0_fp*DEP_HG2, 0e+0_fp )

       ! Add remaining deposited Hg to reservoir for delivery to ocean
       ! when snow melts later (jaf, 6/17/11)
       ! This is Hg in snowpack over ocean that CANNOT be
       ! re-emitted to atmosphere
       SNOW_HG_STORED_OC(I,J) = SNOW_HG_STORED_OC(I,J) + &
                                FRAC_O * FRAC_SNOW_OR_ICE * &
                                MAX( 0.4e+0_fp*DEP_HG2, 0e+0_fp )

       ! This is Hg in snowpack over land that can potentially be
       ! re-emitted to atmosphere
       SNOW_HG_LN(I,J) = SNOW_HG_LN(I,J) + &
                         (1e+0_fp - FRAC_O) * FRAC_SNOW_OR_ICE * &
                         MAX( 0.6e+0_fp*DEP_HG2, 0e+0_fp )

       ! This is Hg in snowpack over land that CANNOT be
       ! re-emitted to atmosphere
       SNOW_HG_STORED_LN(I,J) = SNOW_HG_STORED_LN(I,J) + &
                                   ( 1e+0_fp - FRAC_O ) * &
                                   FRAC_SNOW_OR_ICE &
                                   * MAX( 0.4e+0_fp*DEP_HG2, 0e+0_fp )

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

    ! Assume success
    RC = GC_SUCCESS

    ! Exit if this is a dry-run
    IF ( Input_Opt%DryRun ) RETURN

    !================================================================
    ! Allocate arrays
    !================================================================
    ALLOCATE( DD_Hg2( State_Grid%NX, State_Grid%NY ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'DD_Hg2' )
    DD_Hg2 = 0e+0_fp

    ALLOCATE( DD_HgP( State_Grid%NX, State_Grid%NY ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'DD_HgP' )
    DD_HgP = 0e+0_fp

    ALLOCATE( WD_Hg2( State_Grid%NX, State_Grid%NY ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'WD_Hg2' )
    WD_Hg2 = 0e+0_fp

    ALLOCATE( WD_HgP( State_Grid%NX, State_Grid%NY ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'WD_HgP' )
    WD_HgP = 0e+0_fp

    IF ( Input_Opt%LGTMM ) THEN

       ! GTMM restart file name
       GTMM_RST_FILE = Input_Opt%GTMM_RST_FILE

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

  END SUBROUTINE CLEANUP_DEPO_MERCURY
!EOC
END MODULE DEPO_MERCURY_MOD
