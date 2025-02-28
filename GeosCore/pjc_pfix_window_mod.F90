!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pjc_pfix_window_mod.F90
!
! !DESCRIPTION: Module PJC\_PFIX\_WINDOW\_MOD contains routines which implements
!  the Philip Cameron-Smith pressure fixer.  Specially modified for GEOS-Chem
!  nested grid simulation. (yxw, dan, bmy, 11/6/08)
!\\
!\\
! !INTERFACE:
!
MODULE PJC_PFIX_WINDOW_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: CLEANUP_PJC_PFIX_WINDOW
  PUBLIC :: DO_PJC_PFIX_WINDOW
!
! !REVISION HISTORY:
!  (1 ) Adapted from "pjc_pfix_mod.f" (bmy, 11/6/08)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
  ! ======================================================================
  ! Module Variables:
  ! ======================================================================
  ! AI          : Vertical coord "A" for hybrid grid [hPa]
  ! BI          : Vertical coord "B" for hybrid grid [unitless]
  ! CLAT_FV     : Grid box center latitude [radians]
  ! COSE_FV     : COSINE of grid box edge latitudes [radians]
  ! COSP_FV     : COSINE of grid box center latitudes [radians]
  ! DAP         : Delta-A vertical coordinate [hPa]
  ! DBK         : Delta-B vertical coordinate [unitless]
  ! DLAT_FV     : Latitude extent of grid boxes [radians]
  ! ELAT_FV     : Grid box edge latitudes [radians]
  ! GEOFAC      : Geometric factor for N-S advection
  ! GW_FV       : Diff of SINE btw grid box lat edges [unitless]
  ! MCOR        : Grid box surface areas [m2]
  ! REL_AREA    : Relative surface area of grid box [fraction]
  ! RGW_FV      : Reciprocal of GW_FV [radians
  ! SINE_FV     : SINE of lat at grid box edges [unitless]
  ! GEOFAC_PC   : Geometric factor for N-S advection @ poles
  ! DLON_FV     : Longitude extent of a grid box [radians]
  ! LOC_PROC    : Local processor number
  ! PR_DIAG     : Flag for printing diagnostic message
  ! IMP_NBORDER : Used for ghost zones for MPI ???
  ! I1_GL       : ind of 1st  global lon       (no ghost zones)
  ! I2_GL       : ind of last global lon       (no ghost zones)
  ! JU1_GL      : ind of 1st  global "u" lat   (no ghost zones)
  ! JV1_GL      : ind of 1st  global "v" lat   (no ghost zones)
  ! J2_GL       : ind of last global "u&v" lat (no ghost zones)
  ! K1_GL       : ind of 1st  global alt       (no ghost zones)
  ! K2_GL       : ind of last global alt       (no ghost zones)
  ! ILO_GL      : I1_GL  - IMP_NBORDER        (has ghost zones)
  ! IHI_GL      : I2_GL  + IMP_NBORDER        (has ghost zones)
  ! JULO_GL     : JU1_GL - IMP_NBORDER        (has ghost zones)
  ! JVLO_GL     : JV1_GL - IMP_NBORDER        (has ghost zones)
  ! JHI_GL      : J2_GL  + IMP_NBORDER        (has ghost zones)
  ! I1          : ind of first local lon       (no ghost zones)
  ! I2          : ind of last  local lon       (no ghost zones)
  ! JU1         : ind of first local "u" lat   (no ghost zones)
  ! JV1         : ind of first local "v" lat   (no ghost zones)
  ! J2          : ind of last  local "u&v" lat (no ghost zones)
  ! K1          : index of first local alt     (no ghost zones)
  ! K2          : index of last  local alt     (no ghost zones)
  ! ILO         : I1  - IMP_NBORDER           (has ghost zones)
  ! IHI         : I2  + IMP_NBORDER           (has ghost zones)
  ! JULO        : JU1 - IMP_NBORDER           (has ghost zones)
  ! JVLO        : JV1 - IMP_NBORDER           (has ghost zones)
  ! JHI         : J2  + IMP_NBORDER           (has ghost zones)
  !========================================================================

  ! Allocatable arrays
  REAL(fp), ALLOCATABLE :: AI(:)
  REAL(fp), ALLOCATABLE :: BI(:)
  REAL(fp), ALLOCATABLE :: CLAT_FV(:)
  REAL(fp), ALLOCATABLE :: COSE_FV(:)
  REAL(fp), ALLOCATABLE :: COSP_FV(:)
  REAL(fp), ALLOCATABLE :: DAP(:)
  REAL(fp), ALLOCATABLE :: DBK(:)
  REAL(fp), ALLOCATABLE :: DLAT_FV(:)
  REAL(fp), ALLOCATABLE :: ELAT_FV(:)
  REAL(fp), ALLOCATABLE :: GEOFAC(:)
  REAL(fp), ALLOCATABLE :: GW_FV(:)
  REAL(fp), ALLOCATABLE :: MCOR(:,:)
  REAL(fp), ALLOCATABLE :: REL_AREA(:,:)
  REAL(fp), ALLOCATABLE :: RGW_FV(:)
  REAL(fp), ALLOCATABLE :: SINE_FV(:)

  ! Scalar variables
  LOGICAL               :: PR_DIAG
  INTEGER               :: LOC_PROC
  REAL(fp)              :: GEOFAC_PC
  REAL(fp)              :: DLON_FV

  ! Dimensions for GMI code (from "imp_dims")
  INTEGER               :: IMP_NBORDER
  INTEGER               :: I1_GL,  I2_GL,   JU1_GL,  JV1_GL
  INTEGER               :: J2_GL,  K1_GL,   K2_GL,   ILO_GL
  INTEGER               :: IHI_GL, JULO_GL, JVLO_GL, JHI_GL
  INTEGER               :: I1,     I2,      JU1,     JV1
  INTEGER               :: J2,     K1,      K2,      ILO
  INTEGER               :: IHI,    JULO,    JVLO,    JHI
  INTEGER               :: ILAT,   ILONG,   IVERT,   J1P
  INTEGER               :: J2P

  ! Dimensions for nested grids
  INTEGER               :: I1_W,     I2_W,      JU1_W
  INTEGER               :: J2_W,     J1P_W,     J2P_W
  INTEGER               :: BUFF_SIZE

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_pjc_pfix_window
!
! !DESCRIPTION: Subroutine DO_PJC_PFIX is the driver routine for the Philip
!  Cameron-Smith pressure fixer for the fvDAS transport scheme.
!  (bdf, bmy, 5/8/03, 3/5/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_PJC_PFIX_WINDOW( State_Grid, D_DYN, P1, P2, &
                                 UWND,  VWND, XMASS, YMASS )
!
! !USES:
!
    USE PhysConstants        ! Physical constants
    Use State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    REAL(fp),       INTENT(IN)  :: D_DYN       ! Dynamic timestep [s]
    REAL(fp),       INTENT(IN)  :: P1   (State_Grid%NX, &! True Psurf at middle
                                         State_Grid%NY)  ! dyn tstep [hPa]
    REAL(fp),       INTENT(IN)  :: P2   (State_Grid%NX, &! True Psurf at end of
                                         State_Grid%NY)  ! dyn timestep [hPa]
    REAL(fp),       INTENT(IN)  :: UWND (State_Grid%NX, &! Zonal wind [m/s]
                                         State_Grid%NY, &
                                         State_Grid%NZ)
    REAL(fp),       INTENT(IN)  :: VWND (State_Grid%NX, &! Meridional wind[m/s]
                                         State_Grid%NY, &
                                         State_Grid%NZ)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: XMASS(State_Grid%NX, &! E-W mass fluxes
                                         State_Grid%NY, &!  [kg/s]
                                         State_Grid%NZ)
    REAL(fp),       INTENT(OUT) :: YMASS(State_Grid%NX, &! N-S mass fluxes
                                         State_Grid%NY, &!  [kg/s]
                                         State_Grid%NZ)
!
! !AUTHOR:
!  Brendan Field and Bob Yantosca (5/8/03)
!
! !REMARKS:
!  We assume that the winds are on the A-GRID, since this is the input that
!  the GEOS-4/fvDAS transport scheme takes. (bdf, bmy, 5/8/03)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE        :: FIRST = .TRUE.
    INTEGER              :: I, J
    REAL(fp)             :: P2_TMP(State_Grid%NX,State_Grid%NY)
!
! !DEFINED PARAMETERS:
!
    LOGICAL, PARAMETER   :: INTERP_WINDS     = .TRUE.  ! winds are interp'd
    INTEGER, PARAMETER   :: MET_GRID_TYPE    = 0       ! A-GRID
    INTEGER, PARAMETER   :: ADVEC_CONSRV_OPT = 0       ! 2=floating pressure
    INTEGER, PARAMETER   :: PMET2_OPT        = 1       ! leave at 1
    INTEGER, PARAMETER   :: PRESS_FIX_OPT    = 1       ! Turn on P-Fixer

    !=================================================================
    ! DO_PJC_PFIX begins here!
    !=================================================================

    ! Initialize on first call
    IF ( FIRST ) THEN

       ! Initialize/allocate module variables
       CALL INIT_PJC_PFIX_WINDOW( State_Grid )

       ! Calculate advection surface-area factors
       CALL CALC_ADVECTION_FACTORS( MCOR, REL_AREA, GEOFAC, GEOFAC_PC)

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    ! Copy P2 into P2_TMP (yxw, bmy, 3/5/07)
    P2_TMP = P2

    ! Call PJC pressure fixer w/ the proper arguments
    ! NOTE: P1 and P2 are now "true" surface pressure, not PS-PTOP!!!
    CALL ADJUST_PRESS( State_Grid,                        &
                       'GEOS-CHEM',        INTERP_WINDS,  &
                       .TRUE.,             MET_GRID_TYPE, &
                       ADVEC_CONSRV_OPT,   PMET2_OPT,     &
                       PRESS_FIX_OPT,      D_DYN,         &
                       GEOFAC_PC,          GEOFAC,        &
                       COSE_FV,            COSP_FV,       &
                       REL_AREA,           DAP,           &
                       DBK,                P1,            &
                       P2_TMP,             P2_TMP,        &
                       UWND,               VWND,          &
                       XMASS,              YMASS )

  END SUBROUTINE DO_PJC_PFIX_WINDOW
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Pressure
!
! !DESCRIPTION: Subroutine Calc\_Pressure recalculates the new surface
!  pressure from the adjusted air masses XMASS and YMASS.  This is useful
!  for debugging purposes. (bdf, bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Pressure( State_Grid, XMASS, YMASS, &
                            RGW_FV, PS_NOW, PS_AFTER )
!
! !USES:
!
    USE CMN_SIZE_Mod,   ONLY : PTOP
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid             ! Grid State object
    REAL(fp),       INTENT(IN)  :: XMASS(State_Grid%NX,  &! E-W mass flux from
                                         State_Grid%NY,  &!  pressure fixer
                                         State_Grid%NZ)
    REAL(fp),       INTENT(IN)  :: YMASS(State_Grid%NX,  &! N-S mass flux from
                                         State_Grid%NY,  &!  pressure fixer
                                         State_Grid%NZ)
    REAL(fp),       INTENT(IN)  :: PS_NOW(State_Grid%NX, &! Surface P - PTOP
                                          State_Grid%NY)  !  at current time
    REAL(fp),       INTENT(IN)  :: RGW_FV(State_Grid%NY)  ! Latitude factor
                                                          ! 1/(SINE(J+1)-SIN(J))
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: PS_AFTER(State_Grid%NX,&! Sfc pressure - PTOP
                                            State_Grid%NY) ! adjusted by P-fixer
!
! !AUTHOR:
!   Brendan Field and Bob Yantosca (5/8/03)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I, J, L
    REAL(fp) :: DELP (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp) :: DELP1(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp) :: PE(State_Grid%NX,State_Grid%NZ+1,State_Grid%NY)

    !=================================================================
    ! CALC_PRESSURE begins here!
    !=================================================================
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       DELP1(I,J,L) = DAP(L) + ( DBK(L) * PS_NOW(I,J) )
    ENDDO
    ENDDO
    ENDDO

    DO L = 1, State_Grid%NZ
       DO J = 2, State_Grid%NY-1

          DO I =1, State_Grid%NX-1
             DELP(I,J,L) = DELP1(I,J,L) + &
                           XMASS(I,J,L) - XMASS(I+1,J,L) + &
                         ( YMASS(I,J,L) - YMASS(I,J+1,L) ) * RGW_FV(J)
          ENDDO

          DELP(State_Grid%NX,J,L) = DELP1(State_Grid%NX,J,L) + &
                XMASS(State_Grid%NX,J,L) - XMASS(1,J,L) + &
              ( YMASS(State_Grid%NX,J,L) - YMASS(State_Grid%NX,J+1,L) ) * &
                RGW_FV(J)
       ENDDO

       DO I = 1, State_Grid%NX
          DELP(I,1,L) = DELP1(I,1,L) - YMASS(I,2,L) * RGW_FV(1)
       ENDDO

       ! Compute average
       CALL XPAVG( DELP(1,1,L), State_Grid%NX )

       DO I = 1, State_Grid%NX
          DELP(I,State_Grid%NY,L) = DELP1(I,State_Grid%NY,L) + &
               YMASS(I,State_Grid%NY,L) * RGW_FV(State_Grid%NY)
       ENDDO

       ! Compute average
       CALL XPAVG( DELP(1,State_Grid%NY,L), State_Grid%NX )

    ENDDO

    !=================================================================
    ! Make the pressure
    !=================================================================
    DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          PE(I,1,J) = PTOP
       ENDDO

       DO L = 1,State_Grid%NZ
          DO I = 1,State_Grid%NX
             PE(I,L+1,J) = PE(I,L,J) + DELP(I,J,L)
          ENDDO
       ENDDO

       DO I = 1,State_Grid%NX
          PS_AFTER(I,J) = PE(I,State_Grid%NZ+1,J)
       ENDDO
    ENDDO

  END SUBROUTINE Calc_Pressure
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Advection_Factors
!
! !DESCRIPTION: Subroutine Calc\_Advection\_Factors calculates the relative
!   area of each grid box, and the geometrical factors used by this modified
!   version of TPCORE.  These geomoetrical DO assume that the space is
!   regularly gridded, but do not assume any link between the surface area
!   and the linear dimensions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Advection_Factors(mcor, rel_area, geofac, geofac_pc)
!
! !USES:
!
    USE PhysConstants  ! Physical constants
!
! !INPUT PARAMETERS:
!
    ! Area of grid box (m^2)
    REAL(fp), INTENT(IN)  :: mcor(i1_gl :i2_gl, ju1_gl:j2_gl)
!
! !OUTPUT PARAMETERS:
!
    ! relative surface area of grid box (fraction)
    REAL(fp), INTENT(OUT) :: rel_area(i1_gl :i2_gl, ju1_gl:j2_gl)

    ! Geometrical factor for meridional advection; geofac uses
    ! correct spherical geometry, and replaces acosp as the
    ! meridional geometrical factor in tpcore
    REAL(fp), INTENT(OUT) :: geofac(ju1_gl:j2_gl)

    ! Special geometrical factor (geofac) for Polar cap
    REAL(fp), INTENT(OUT) :: geofac_pc
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: ij

    ! Variables not used. (ccc, 8/3/10)
    !REAL(fp) :: dp           ! spacing in latitude (rad)
    !REAL(fp) :: ri2_gl
    !REAL(fp) :: rj2m1
    REAL(fp) :: total_area

    !----------------
    !Begin execution.
    !----------------

    ! Not used. (ccc, 8/3/10)
    !ri2_gl = i2_gl

    !---------------------------------
    !Set the relative area (rel_area).
    !---------------------------------

    total_area = Sum (mcor(:,:))

    rel_area(:,:) = mcor(:,:) / total_area

    !---------------------------------------------------------
    !Calculate geometrical factor for meridional advection.
    !Note that it is assumed that all grid boxes in a latitude
    !band are the same.
    !---------------------------------------------------------

    ! Not used for nested grids. (ccc, 8/3/10)
    !rj2m1 = j2_gl - 1
    !dp    = PI / 360D0

    ! The total area does not cover the full globe so use an other definition
    ! for the geometric factor. (lzh, ccc, 8/3/10)
    do ij = ju1_gl, j2_gl
       !geofac(ij) = dp / (2.0e+0_fp * rel_area(1,ij) * ri2_gl)
       geofac(ij) = 1.e+0_fp / COSP_FV(ij)
    end do

    ! geofac_pc used only for polar cap so no need. (ccc, 8/3/10)
    !geofac_pc = dp / (2.0e+0_fp * Sum (rel_area(1,ju1_gl:ju1_gl+1)) * ri2_gl)

    ! Make sure to return with a value
    geofac_pc = 0.0_fp

  END SUBROUTINE Calc_Advection_Factors
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adjust_Press
!
! !DESCRIPTION: Subroutine Adjust\_Press initializes and calls the
!  pressure fixer code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Adjust_Press(State_Grid, &
       metdata_name_org, do_timinterp_winds, new_met_rec,         &
       met_grid_type, advec_consrv_opt, pmet2_opt, press_fix_opt, &
       tdt, geofac_pc, geofac, cose, cosp, rel_area, dap, dbk,    &
       pctm1, pctm2, pmet2, uu, vv, xmass, ymass)
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object

    ! First  part of metdata_name, e.g., "NCAR"
    CHARACTER(LEN=*) :: metdata_name_org

    ! Time interpolate wind fields?
    LOGICAL :: do_timinterp_winds

    ! New met record?
    LOGICAL :: new_met_rec

    ! Met grid type, A or C
    INTEGER :: met_grid_type

    ! Advection_conserve option
    INTEGER :: advec_consrv_opt

    ! pmet2 option
    INTEGER :: pmet2_opt

    ! pressure fixer option
    INTEGER :: press_fix_opt

    ! Model time step [s]
    REAL(fp)  :: tdt

    ! Special geometrical factor (geofac) for Polar cap
    REAL(fp)  :: geofac_pc

    ! Geometrical factor for meridional advection; geofac uses
    ! correct spherical geometry, and replaces acosp as the
    ! meridional geometrical factor in tpcore
    REAL(fp)  :: geofac  (ju1_gl:j2_gl)

    ! Cosines of grid box edges and centers
    REAL(fp)  :: cose    (ju1_gl:j2_gl)
    REAL(fp)  :: cosp    (ju1_gl:j2_gl)

    ! Pressure difference across layer from (ai * pt) term [hPa]
    REAL(fp)  :: dap     (k1:k2)

    ! Difference in bi across layer - the dSigma term
    REAL(fp)  :: dbk     (k1:k2)

    ! Relative surface area of grid box (fraction)
    REAL(fp)  :: rel_area( i1_gl:i2_gl,   ju1_gl:j2_gl)

    ! Metfield surface pressure at t1+tdt [hPa]
    REAL(fp)  :: pmet2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

    ! CTM surface pressure at t1 [hPa]
    REAL(fp)  :: pctm1(ilo_gl:ihi_gl, julo_gl:jhi_gl)

    ! CTM surface pressure at t1+tdt [hPa]
    REAL(fp)  :: pctm2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

    ! Wind velocity, x direction at t1+tdt/2 [m/s]
    REAL(fp)  :: uu(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

    ! Wind velocity, y direction at t1+tdt/2 [m/s]
    REAL(fp)  :: vv(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Horizontal mass flux in E-W direction [hPa]
    REAL(fp)  :: xmass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

    ! Horizontal mass flux in N-S direction [hPa]
    REAL(fp)  :: ymass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    logical, save :: DO_ADJUST_PRESS_DIAG = .TRUE.

    !----------------------
    !Variable declarations.
    !----------------------

    logical, save :: first = .true.

    !--------------------------------------------------
    !dgpress   : global-pressure discrepancy
    !press_dev : RMS difference between pmet2 and pctm2
    !            (weighted by relative area)
    !--------------------------------------------------
    real(fp)  :: dgpress
    real(fp)  :: press_dev

    !-------------------------------------------------------------
    !dps : change of surface pressure from met field pressure (mb)
    !-------------------------------------------------------------
    real(fp)  :: dps(i1_gl:i2_gl, ju1_gl:j2_gl)

    !--------------------------------------------
    !dps_ctm : CTM surface pressure tendency (mb)
    !--------------------------------------------
    real(fp) :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)

    !---------------------------------------------------------------------
    !xmass_fixed : horizontal mass flux in E-W direction after fixing (mb)
    !ymass_fixed : horizontal mass flux in N-S direction after fixing (mb)
    !---------------------------------------------------------------------
    real(fp)  :: xmass_fixed(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
    real(fp)  :: ymass_fixed(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)

    !-------------
    !Dummy indexes
    !-------------
    !integer :: ij, il

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6, *) 'Adjust_Press called by ', loc_proc
    end if

    dps     = 0.0_fp    ! (bmy, 7/10/17)
    dps_ctm = 0.0e+0_fp

    dgpress =  Sum ( (pmet2(i1_gl:i2_gl, ju1_gl:j2_gl) - &
                      pctm1(i1_gl:i2_gl, ju1_gl:j2_gl)   ) &
                 * rel_area(i1_gl:i2_gl, ju1_gl:j2_gl)     )

    if (pmet2_opt == 1) then
       pmet2(:,:) = pmet2(:,:) - dgpress
    end if

    !### Debug
    !###if (DO_ADJUST_PRESS_DIAG) then
    !###  Write (6, *) 'Global mean surface pressure change (mb) = ',
    !###                dgpress
    !###end if

    !===================
    call Init_Press_Fix &
    !===================
         (State_Grid, &
          metdata_name_org, met_grid_type, tdt, geofac_pc, geofac, &
          cose, cosp, dap, dbk, dps, dps_ctm, rel_area, pctm1, pmet2, &
          uu, vv, xmass, ymass)

    if (press_fix_opt == 1) then

       !======================
       call Do_Press_Fix_Llnl &
       !======================
            (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area, &
             xmass, ymass, xmass_fixed, ymass_fixed )

       xmass(:,:,:) = xmass_fixed(:,:,:)
       ymass(:,:,:) = ymass_fixed(:,:,:)

    end if

    if ((advec_consrv_opt == 0) .or. &
        (advec_consrv_opt == 1)) then

       dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl) = &
            pmet2(i1_gl:i2_gl, ju1_gl:j2_gl) - &
            pctm1(i1_gl:i2_gl, ju1_gl:j2_gl)

       !-----------------------------------------------
       !else if (advec_consrv_opt == 2) then do nothing
       !-----------------------------------------------

    end if

    pctm2(i1_gl:i2_gl, ju1_gl:j2_gl) = &
         pctm1(i1_gl:i2_gl, ju1_gl:j2_gl) + &
         dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)

    if (DO_ADJUST_PRESS_DIAG) then

       !-------------------------------------------------------
       !Calculate the RMS pressure deviation (diagnostic only).
       !-------------------------------------------------------

       press_dev = &
            Sqrt (Sum (((pmet2(i1_gl:i2_gl,ju1_gl:j2_gl) - &
                         pctm2(i1_gl:i2_gl,ju1_gl:j2_gl))**2 * &
                      rel_area(i1_gl:i2_gl,ju1_gl:j2_gl))))

       !### Debug
       !###Write (6, *) 'RMS deviation between pmet2 & pctm2 (mb) = ', &
       !###             press_dev

    end if

  END SUBROUTINE Adjust_Press
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Press_Fix
!
! !DESCRIPTION: Subroutine Init\_Press\_Fix initializes the pressure fixer.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Press_Fix( State_Grid, &
       metdata_name_org, met_grid_type, tdt, geofac_pc, geofac,    &
       cose, cosp, dap, dbk, dps, dps_ctm, rel_area, pctm1, pmet2, &
       uu, vv, xmass, ymass)
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    ! Grid State object
    TYPE(GrdState), INTENT(IN) :: State_Grid

    ! First part of metdata_name, e.g., "NCAR"
    CHARACTER(LEN=*) :: metdata_name_org

    ! Met grid type, A or C
    INTEGER          :: met_grid_type

    ! Model Time step [s]
    REAL(fp)         :: tdt

    ! Special geometrical factor (geofac) for Polar cap
    REAL(fp)         :: geofac_pc

    ! Cosine of grid box edges and centers
    REAL(fp)         :: cose    (ju1_gl:j2_gl)
    REAL(fp)         :: cosp    (ju1_gl:j2_gl)

    ! Geometrical factor for meridional advection; geofac uses
    ! correct spherical geometry, and replaces acosp as the
    ! meridional geometrical factor in tpcore
    REAL(fp)         :: geofac  (ju1_gl:j2_gl)

    ! Pressure difference across layer from (ai * pt) term [hPa]
    REAL(fp)         :: dap     (k1:k2)

    ! Difference in bi across layer - the dSigma term
    REAL(fp)         :: dbk     (k1:k2)

    !rel_area : relative surface area of grid box (fraction)
    REAL(fp)         :: rel_area( i1_gl:i2_gl,   ju1_gl:j2_gl)

    ! Metfield surface pressure at t1 [hPa]
    REAL(fp)         ::  pmet2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

    ! CTM surface pressure at t1 [hPa]
    REAL(fp)         ::  pctm1(ilo_gl:ihi_gl, julo_gl:jhi_gl)

    ! CTM surface pressure at t1+tdt [hPa]
    REAL(fp)         ::  pctm2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

    ! Wind velocity, x direction at t1+tdt/2 [m/s]
    REAL(fp)         :: uu(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
    
    ! Wind velocity, y direction at t1+tdt/2 [m/s]
    REAL(fp)         :: vv(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
!
! !OUTPUT PARAMETERS:
!
    ! Horizontal mass flux in E-W direction [hPa]
    REAL(fp)         :: xmass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

    ! Horizontal mass flux in N-S direction [hPa]
    REAL(fp)         :: ymass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

    ! Change of surface pressure from met field pressure (mb)
    REAL(fp)         :: dps(i1_gl:i2_gl, ju1_gl:j2_gl)

    ! CTM surface pressure tendency (mb)
    REAL(fp)         :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !--------------------------------------------------------------
    !dpi   : divergence at a grid point; used to calculate vertical
    !        motion (mb)
    !--------------------------------------------------------------
    real(fp)  :: dpi(i1:i2, ju1:j2, k1:k2)

    !---------------------------------------------------------------------
    !crx   : Courant number in E-W direction
    !cry   : Courant number in N-S direction
    !delp1 : pressure thickness, the psudo-density in a hydrostatic system
    !        at t1 (mb)
    !delpm : pressure thickness, the psudo-density in a hydrostatic system
    !        at t1+tdt/2 (approximate) (mb)
    !pu    : pressure at edges in "u"  (mb)
    !---------------------------------------------------------------------
    real(fp)  :: crx  (ilo:ihi, julo:jhi, k1:k2)
    real(fp)  :: cry  (ilo:ihi, julo:jhi, k1:k2)
    real(fp)  :: delp1(ilo:ihi, julo:jhi, k1:k2)
    real(fp)  :: delpm(ilo:ihi, julo:jhi, k1:k2)
    real(fp)  :: pu   (ilo:ihi, julo:jhi, k1:k2)

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6,*) 'Init_Press_Fix called by ', loc_proc
    end if

    ! Need to initialize local arrays (bmy, 7/10/17)
    dpi   = 0.0_fp
    crx   = 0.0_fp
    cry   = 0.0_fp
    delp1 = 0.0_fp
    delpm = 0.0_fp
    pu    = 0.0_fp

    ! not treat poles (lzh, 07/20/2010)
    !!========================
    !call Average_Press_Poles &
    !!========================
    !  (rel_area, pctm1)
    !
    !!========================
    !call Average_Press_Poles &
    !!========================
    !  (rel_area, pmet2)

    !-------------------------------------------------------------------
    !We need to calculate pressures at t1+tdt/2.  One ought to use pctm2
    !in the call to Set_Press_Terms, but since we don't know it yet, we
    !are forced to use pmet2.  This should be good enough because it is
    !only used to convert the winds to the mass fluxes, which is done
    !crudely anyway and the mass fluxes will still get fixed OK.
    !-------------------------------------------------------------------

    dps(i1:i2,ju1:j2) = pmet2(i1:i2,ju1:j2) - pctm1(i1:i2,ju1:j2)

    !====================
    call Set_Press_Terms &
    !====================
         (dap, dbk, pctm1, pmet2, delp1, delpm, pu)

    !===================
    call Convert_Winds &
    !===================
         (State_Grid, met_grid_type, tdt, cosp, crx, cry, uu, vv)

    !=========================
    call Calc_Horiz_Mass_Flux &
    !=========================
         (State_Grid, cose, delpm, uu, vv, xmass, ymass, tdt, cosp)

    !====================
    call Calc_Divergence &
    !====================
         (.false., geofac_pc, geofac, dpi, xmass, ymass)

    dps_ctm(i1:i2,ju1:j2) = Sum (dpi(i1:i2,ju1:j2,:), dim=3)

  END SUBROUTINE Init_Press_Fix
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Press_Fix_Llnl
!
! !DESCRIPTION: Subroutine Do\_Press\_Fix\_Llnl fixes the mass fluxes to
!  match the met field pressure tendency.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Press_Fix_Llnl &
       (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area, &
       xmass, ymass, xmass_fixed, ymass_fixed)
!
! !INPUT PARAMETERS:
!
    ! Special geometrical factor (geofac) for Polar cap
    REAL(fp), INTENT(IN)  :: geofac_pc

    ! Geometrical factor for meridional advection; geofac uses
    ! correct spherical geometry, and replaces acosp as the
    !  meridional geometrical factor in tpcore
    REAL(fp), INTENT(IN)  :: geofac  (ju1_gl:j2_gl)

    ! Difference in bi across layer - the dSigma term
    REAL(fp), INTENT(IN)  :: dbk     (k1:k2)

    ! Change of surface pressure from met field pressure [hPa]
    REAL(fp), INTENT(IN)  :: dps     (i1:i2, ju1:j2)

    ! Relative surface area of grid box (fraction)
    REAL(fp), INTENT(IN)  :: rel_area(i1:i2, ju1:j2)

    ! Horizontal mass fluxes in E-W and N-S directions [hPa]
    REAL(fp), INTENT(IN)  :: xmass      (ilo:ihi, julo:jhi, k1:k2)
    REAL(fp), INTENT(IN)  :: ymass      (ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
    ! Sum over vertical of dpi calculated from original mass fluxes [hPa]
    REAL(fp), INTENT(OUT) :: dps_ctm (i1:i2, ju1:j2)

    ! Horizontal mass flux in E-W and N-S directions after fixing [hPa]
    REAL(fp), INTENT(OUT) :: xmass_fixed(ilo:ihi, julo:jhi, k1:k2)
    REAL(fp), INTENT(OUT) :: ymass_fixed(ilo:ihi, julo:jhi, k1:k2)
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: il, ij, ik

    REAL(fp)  :: dgpress
    REAL(fp)  :: fxmean
    REAL(fp)  :: ri2

    ! Arrays
    REAL(fp)  :: fxintegral(i1:i2+1)
    REAL(fp)  :: mmfd(ju1:j2)
    REAL(fp)  :: mmf (ju1:j2)
    REAL(fp)  :: ddps(i1:i2, ju1:j2)

    !------------------------------------------------------------------------
    !dpi : divergence at a grid point; used to calculate vertical motion (mb)
    !------------------------------------------------------------------------
    REAL(fp)  :: dpi(i1:i2, ju1:j2, k1:k2)
    REAL(fp)  :: xcolmass_fix(ilo:ihi, julo:jhi)
    REAL(fp)  :: xx

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6,*) 'Do_Press_Fix_Llnl called by ', loc_proc
    end if

    !------------------------------------------------------------
    ! Initialize variables (bmy, 7/10/17)
    !------------------------------------------------------------

    ! Arguments
    xmass_fixed  = xmass
    ymass_fixed  = ymass

    ! Local variables
    dgpress      = 0.0_fp
    fxmean       = 0.0_fp
    ri2          = i2_gl - 2 * BUFF_SIZE
    fxintegral   = 0.0_fp
    mmfd         = 0.0_fp
    mmf          = 0.0_fp
    ddps         = 0.0_fp
    dpi          = 0.0_fp
    xcolmass_fix = 0.0_fp
    xx           = 0.0_fp

    !------------------------------------------------------------
    !Calculate difference between GCM and LR predicted pressures.
    !------------------------------------------------------------

    ddps(:,:) = dps(:,:) - dps_ctm(:,:)

    !--------------------------------------
    ! Calculate global-pressure discrepancy.
    !--------------------------------------

    !dgpress = Sum (ddps(i1:i2,ju1:j2) * rel_area(i1:i2,ju1:j2))

    xx = sum(rel_area(i1_w:i2_w,j1p_w:j2p_w))
    dgpress = Sum (ddps(i1_w:i2_w,j1p_w:j2p_w) * &
         rel_area(i1_w:i2_w,j1p_w:j2p_w)) / xx

    !----------------------------------------------------------
    !Calculate mean meridional flux divergence (df/dy).
    !Note that mmfd is actually the zonal mean pressure change,
    !which is related to df/dy by geometrical factors.
    !----------------------------------------------------------

    !------------------------
    !Handle non-Pole regions.
    !------------------------

    ! Work on the inner window only (lzh, ccc, 8/3/10)
    !do ij = j1p, j2p
    !  mmfd(ij) = -(sum(ddps(:,ij)) / ri2 - dgpress)
    !end do

    do ij = j1p_w, j2p_w
       mmfd(ij) = -(sum(ddps(i1_w:i2_w,ij)) / ri2 - dgpress)
    end do

    ! No special case for poles, no poles. (ccc, 8/3/10)
    !!---------------------------------------------
    !!Handle poles.
    !!Note that polar boxes have all been averaged.
    !!---------------------------------------------
    !
    !mmfd(ju1)   = -(ddps(1,ju1)   - dgpress)
    !mmfd(ju1+1) = -(ddps(1,ju1+1) - dgpress)
    !mmfd(j2-1)  = -(ddps(1,j2-1)  - dgpress)
    !mmfd(j2)    = -(ddps(1,j2)    - dgpress)

    !---------------------------------------------
    !Calculate mean meridional fluxes (cos(e)*fy).
    !---------------------------------------------

    ! Use geofac, no polar cap. (ccc, 8/3/10)
    !mmf(j1p) = mmfd(ju1) / geofac_pc
    mmf(j1p_w) = mmfd(ju1_w) / geofac(j1p_w)

    ! Work on inner domain. (ccc, 8/3/10)
    !do ij = j1p, j2p
    do ij = j1p_w, j2p_w-1
       mmf(ij+1) = mmf(ij) + mmfd(ij) / geofac(ij)
    end do

    !------------------------------------------------------------
    !Fix latitude bands.
    !Note that we don't need to worry about geometry here because
    !all boxes in a latitude band are identical.
    !Note also that fxintegral(i2+1) should equal fxintegral(i1),
    !i.e., zero.
    !------------------------------------------------------------

    ! Work on inner domain (ccc, 8/3/10)
    !do ij = j1p, j2p
    do ij = j1p_w, j2p_w

       fxintegral(:) = 0.0e+0_fp

       !do il = i1, i2
       do il = i1_w, i2_w
          fxintegral(il+1) = fxintegral(il) - &
                             (ddps(il,ij) - dgpress) - mmfd(ij)
       end do

       fxmean = Sum (fxintegral(i1+1:i2+1)) / ri2
       !fxmean = Sum (fxintegral(i1_w+1:i2_w+1)) / ri2

       !do il = i1, i2
       do il = i1_w, i2_w
          xcolmass_fix(il,ij) = fxintegral(il) - fxmean
       end do

    end do

    !-------------------------------------
    !Distribute colmass_fix's in vertical.
    !-------------------------------------

    do ik = k1, k2
    !do ij = j1p, j2p
    !do il = i1, i2
    do ij = j1p_w, j2p_w
    do il = i1_w, i2_w

       xmass_fixed(il,ij,ik) = xmass(il,ij,ik) + xcolmass_fix(il,ij) * dbk(ik)

    end do
    end do
    end do

    ! Grid stops at j2p if nested domain (ccc, 8/3/10)
    !do ik = k1, k2
    !do ij = j1p, j2p+1
    !do il = i1, i2
    !
    !   ymass_fixed(il,ij,ik) = ymass(il,ij,ik) + mmf(ij) * dbk(ik)
    !
    !end do
    !end do
    !end do

    do ik = k1, k2
    !do ij = j1p, j2p
    !do il = i1, i2
    do ij = j1p_w, j2p_w
    do il = i1_w, i2_w

       ymass_fixed(il,ij,ik) = ymass(il,ij,ik) + mmf(ij) * dbk(ik)

    end do
    end do
    end do

    !====================
    call Calc_Divergence &
    !====================
         (.false., geofac_pc, geofac, dpi, xmass_fixed, ymass_fixed)

    dps_ctm(i1:i2,ju1:j2) = Sum (dpi(i1:i2,ju1:j2,:), dim=3)

  END SUBROUTINE Do_Press_Fix_Llnl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Convert_Winds
!
! !DESCRIPTION: Subroutine Convert\_Winds converts winds on A or C grid to
!  Courant \# on C grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Winds(State_Grid, igd, tdt, cosp, crx, cry, uu, vv)
!
! !USES:
!
    USE PhysConstants ! Re, PI
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    ! Grid State object
    TYPE(GrdState), INTENT(IN) :: State_Grid

    ! A or C grid
    INTEGER,  INTENT(IN)  :: igd

    ! Model time step [s]
    REAL(fp), INTENT(IN)  :: tdt

    ! Cosine of grid box centers
    REAL(fp), INTENT(IN)  :: cosp(ju1_gl:j2_gl)

    ! Wind velocity in E-W (UU) and N-S (VV) directions at t1+tdt/2 [m/s]
    REAL(fp), INTENT(IN)  :: uu  (ilo:ihi, julo:jhi, k1:k2)
    REAL(fp), INTENT(IN)  :: vv  (ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
    ! Courant number in E-W (CRX) and N-S (CRY) directions
    REAL(fp), INTENT(OUT) :: crx (ilo:ihi, julo:jhi, k1:k2)
    REAL(fp), INTENT(OUT) :: cry (ilo:ihi, julo:jhi, k1:k2)
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REMARKS:
!  Use GEOS-CHEM physical constants Re, PI to be consistent with other
!  usage everywhere (bmy, 5/5/03)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    logical, save :: first = .true.

    integer :: il, ij

    !-------------------------------
    !dl : spacing in longitude (rad)
    !dp : spacing in latitude  (rad)
    !-------------------------------
    real(fp)  :: dl
    real(fp)  :: dp

    real(fp)  :: ri2
    real(fp)  :: rj2m1

    !------------------------
    !dtdy  : dt/dy      (s/m)
    !dtdy5 : 0.5 * dtdy (s/m)
    !------------------------
    real(fp), save :: dtdy
    real(fp), save :: dtdy5

    !------------------------
    !dtdx  : dt/dx      (s/m)
    !dtdx5 : 0.5 * dtdx (s/m)
    !------------------------
    real(fp), allocatable, save :: dtdx (:)
    real(fp), allocatable, save :: dtdx5(:)

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6, *) 'Convert_Winds called by ', loc_proc
    end if

    !==========
    if (first) then
    !==========

       first = .false.

       Allocate (dtdx (ju1_gl:j2_gl))
       Allocate (dtdx5(ju1_gl:j2_gl))
       dtdx = 0.0_fp; dtdx5 = 0.0_fp

       ri2   = i2_gl
       rj2m1 = j2_gl - 1

       IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
          dl    = 2.0e+0_fp * PI / 1152.0_fp
          dp    = PI /720e+0_fp
       ELSEIF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
          dl    = 2.0e+0_fp * PI / 576.0_fp
          dp    = PI /360e+0_fp
       ENDIF

       dtdy  = tdt / (Re * dp)
       dtdy5 = 0.5_fp * dtdy

       !-----lzh----------
       !dtdx (ju1_gl) = 0.0e+0_fp
       !dtdx5(ju1_gl) = 0.0e+0_fp
       !
       !do ij = ju1_gl + 1, j2_gl - 1
       !
       !  dtdx (ij) = tdt / (dl * Re * cosp(ij))
       !  dtdx5(ij) = 0.5e+0_fp * dtdx(ij)
       !
       !end do
       !
       !dtdx (j2_gl)  = 0.0e+0_fp
       !dtdx5(j2_gl)  = 0.0e+0_fp

       !-----------------------------------------------
       ! for nested NA or EA (lzh, 07/20/2010)
       do ij = ju1_gl, j2_gl
          dtdx (ij) = tdt / (dl * Re * cosp(ij))
          dtdx5(ij) = 0.5e+0_fp * dtdx(ij)
       end do
       !-----------------------------------------------

    end if

    !=============
    if (igd == 0) then  ! A grid.
    !=============

       do ij = ju1+1, j2-1
          do il = i1+1, i2
             crx(il,ij,:) = dtdx5(ij) * &
                            (uu(il,ij,:) + uu(il-1,ij,  :))
          end do
          ! No periodicity (ccc, 8/3/10)
          !crx(1,ij,:) = dtdx5(ij) * &
          !              (uu(1,ij,:) + uu(i2,ij,  :))
       end do

       do ij = ju1+1, j2
          do il = i1, i2
             cry(il,ij,:) = dtdy5 * &
                            (vv(il,ij,:) + vv(il,  ij-1,:))
          end do
       end do

    !====
    else  ! C grid.
    !====

       ! No ghost zones. (ccc, 8/3/10)
       !do ij = ju1, j2
       !  do il = i1, i2
       do ij = ju1+1, j2
          do il = i1+1, i2

             crx(il,ij,:) = dtdx(ij) * uu(il-1,ij,  :)

             cry(il,ij,:) = dtdy     * vv(il,  ij-1,:)

          end do
       end do

    end if

  END SUBROUTINE Convert_Winds
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Horiz_Mass_Flux
!
! !DESCRIPTION: Subroutine Calc\_Horiz\_Mass\_Flux calculates the horizontal
!  mass flux for non-GISS met data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Horiz_Mass_Flux &
    (State_Grid, cose, delpm, uu, vv, xmass, ymass, tdt, cosp)
!
! !USES:
!
    USE PhysConstants ! Re, Pi
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    ! Grid State object
    TYPE(GrdState), INTENT(IN) :: State_Grid

    ! Timestep [s]
    REAL(fp)  :: tdt

    ! Cosine of grid box edges
    REAL(fp)  :: cose (ju1_gl:j2_gl)

    ! Cosine of grid box centers
    REAL(fp)  :: cosp (ju1_gl:j2_gl)

    ! Pressure thickness, the pseudo-density in a
    ! hdrostatic system  at t1+tdt/2 (approximate) [hPa]
    REAL(fp)  :: delpm(ilo:ihi, julo:jhi, k1:k2)

    ! E-W (UU) and N-S (VV) winds [m/s]
    REAL(fp)  :: uu  (ilo:ihi, julo:jhi, k1:k2)
    REAL(fp)  :: vv  (ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
    ! Horizontal mass flux in E-W and N-S directions [hPa]
    REAL(fp)  :: xmass(ilo:ihi, julo:jhi, k1:k2)
    REAL(fp)  :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO
!   John Tannahill, LLNL (jrt@llnl.gov)
!
! !REMARKS:
!   Use GEOS-CHEM physical constants Re, PI to be consistent with other
!   usage everywhere (bmy, 5/5/03)

! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer   :: ij
    integer   :: il
    integer   :: jst, jend
    real(fp)  :: dl
    real(fp)  :: dp

    real(fp)  :: ri2
    real(fp)  :: rj2m1
    real(fp)  :: factx
    real(fp)  :: facty

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6,*) 'Calc_Horiz_Mass_Flux called by ', loc_proc
    end if

    ri2   = i2_gl
    rj2m1 = j2_gl - 1

    IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
       dl    = 2.0e+0_fp * PI / 1152.0_fp
       dp    = PI /720.0_fp
    ELSEIF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
       dl    = 2.0e+0_fp * PI / 576.0_fp
       dp    = PI /360.0_fp
    ENDIF

    facty  = 0.5_fp * tdt / (Re * dp)

    !-----------------------------------
    !Calculate E-W horizontal mass flux.
    !-----------------------------------

    do ij = ju1, j2

       factx = 0.5_fp * tdt / (dl * Re * cosp(ij))

       do il = i1+1, i2
          xmass(il,ij,:) = factx * &
                           (uu(il,ij,:) * delpm(il,ij,:) + &
                            uu(il-1,ij,:) * delpm(il-1,ij,:))
       end do

       ! No periodicity. (ccc, 8/3/10)
       !xmass(i1,ij,:) = factx * &
       !  (uu(i1,ij,:) * delpm(i1,ij,:)+ &
       !   uu(i2,ij,:) * delpm(i2,ij,:))

    end do

    !-----------------------------------
    !Calculate N-S horizontal mass flux.
    !-----------------------------------

    do ij = ju1+1, j2

       ymass(i1:i2,ij,:) = facty * &
            cose(ij) * (vv(i1:i2,ij,:)*delpm(i1:i2,ij,:) + &
            vv(i1:i2,ij-1,:)*delpm(i1:i2,ij-1,:))

    end do

  END SUBROUTINE Calc_Horiz_Mass_Flux
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Divergence
!
! !DESCRIPTION: Subroutine Calc\_Divergence calculates the divergence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Divergence &
       (do_reduction, geofac_pc, geofac, dpi, xmass, ymass)
!
! !INPUT PARAMETERS:
!
    ! Set to F if called on Master; set to T if called by Slaves
    ! (NOTE: this doesn't seem to be used!)
    LOGICAL,  INTENT(IN)  :: do_reduction

    ! Special geometrical factor (geofac) for Polar cap
    REAL(fp), INTENT(IN)  :: geofac_pc

    ! geometrical factor for meridional advection; geofac uses
    ! correct spherical geometry, and replaces acosp as the
    ! meridional geometrical factor in tpcore
    REAL(fp), INTENT(IN)  :: geofac(ju1_gl:j2_gl)

    ! horizontal mass fluxes in E-W and N-S directions [hPa]
    REAL(fp), INTENT(IN)  :: xmass (ilo:ihi, julo:jhi, k1:k2)
    REAL(fp), INTENT(IN)  :: ymass (ilo:ihi, julo:jhi, k1:k2)
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Divergence at a grid point; used to calculate vertical motion [hPa]
    REAL(fp), INTENT(OUT) :: dpi   (i1:i2, ju1:j2, k1:k2)
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REVISION HISTORY:
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL(fp).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: il, ij
    !integer :: jst, jend

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6,*) 'Calc_Divergence called by ', loc_proc
    end if

    ! Initialize DPI, it's an output
    dpi = 0.0_fp

    !-------------------------
    !Calculate N-S divergence.
    !-------------------------

    ! No polar cap. (ccc, 8/3/10)
    !do ij = j1p, j2p
    !
    !  dpi(i1:i2,ij,:) =
    !    (ymass(i1:i2,ij,:) - ymass(i1:i2,ij+1,:)) * geofac(ij)
    !
    !end do
    !
    !if(j1p.ne.2) then
    !  dpi(:,2,:) = 0.
    !  dpi(:,j2-1,:) = 0.
    !endif

    !do ij = j1p_w, j2p_w
    do ij = j1p, j2p-1

       dpi(i1:i2,ij,:) = (ymass(i1:i2,ij,:) - ymass(i1:i2,ij+1,:)) * geofac(ij)

    end do

    !-----lzh-----------------------
    !!===========================
    !call Do_Divergence_Pole_Sum &
    !!===========================
    !  (do_reduction, geofac_pc, dpi, ymass)
    ! comment out for nested NA (lzh, 07/20/2010)
    !dpi(:,1,:) = 0.  ! (lzh, 07/20/2010)
    !qpi(:,j2,:) = 0.

    !-------------------------
    !Calculate E-W divergence.
    !-------------------------

    do ij = j1p,j2p
       do il = i1, i2-1
          dpi(il,ij,:) = dpi(il,ij,:) + xmass(il,ij,:) - xmass(il+1,ij,:)
       end do
       ! No periodicity. (ccc, 8/3/10)
       !dpi(i2,ij,:) = dpi(i2,ij,:) + xmass(i2,ij,:) - xmass(1,ij,:)
    end do

  END SUBROUTINE Calc_Divergence
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Press_Terms
!
! !DESCRIPTION: Subroutine Set\_Press\_Terms sets the pressure terms.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Press_Terms &
       (dap, dbk, pres1, pres2, delp1, delpm, pu)
!
! !INPUT PARAMETERS:
!
    ! Pressure difference across layer from (ai * pt) term [hPa]
    REAL(fp), INTENT(IN)  :: dap  (k1:k2)

    ! Difference in bi across layer - the dSigma term
    REAL(fp), INTENT(IN)  :: dbk  (k1:k2)

    ! Surface pressure at t1 [hPa]
    REAL(fp), INTENT(IN)  :: pres1(ilo:ihi, julo:jhi)

    ! Surface pressure at t1+tdt [hPa]
    REAL(fp), INTENT(IN)  :: pres2(ilo:ihi, julo:jhi)
!
! !OUTPUT PARAMETERS:
!
    ! Pressure thickness, the psudo-density in a
    ! hydrostatic system at t1 [hPa]
    REAL(fp), INTENT(OUT) :: delp1(ilo:ihi, julo:jhi, k1:k2)

    ! Pressure thickness, the psudo-density in a
    ! hydrostatic system at t1+tdt/2 (approximate)  [hPa]
    REAL(fp), INTENT(OUT) :: delpm(ilo:ihi, julo:jhi, k1:k2)

    ! Pressure at edges in "u" [hPa]
    REAL(fp), INTENT(OUT) :: pu   (ilo:ihi, julo:jhi, k1:k2)
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: il, ij, ik
    integer :: jst, jend

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6,*) 'Set_Press_Terms called by ', loc_proc
    end if

    do ik = k1, k2

       delp1(:,:,ik) = dap(ik) + (dbk(ik) * pres1(:,:))

       delpm(:,:,ik) = dap(ik) + &
            (dbk(ik) * 0.5e+0_fp * (pres1(:,:) + pres2(:,:)))

    end do

    do ij = ju1, j2
       do il = i1+1, i2
          pu(il,ij,:) = 0.5e+0_fp * (delpm(il,ij,:) + delpm(il-1,ij,:))
       end do

       ! No periodicity. (ccc, 8/3/10)
       !pu(i1,ij,:) = 0.5e+0_fp * (delpm(i1,ij,:) + delpm(i2,ij,:))

    end do

  END SUBROUTINE Set_Press_Terms
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Xpavg
!
! !description: Subroutine Xpavg replaces each element of a vector with
!  the average of the entire array. (bmy, 5/7/03)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE XPAVG( P, IM )
!
! !USES:
!
    USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
    ! Dimension of P
    INTEGER,  INTENT(IN)    :: IM
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! 1-D vector to be averaged
    REAL(fp), INTENT(INOUT) :: P(IM)
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL                   :: AVG

    !=================================================================
    ! XPAVG begins here!
    !=================================================================

    ! Error check IM
    IF ( IM == 0 ) THEN
       CALL ERROR_STOP( 'Div by zero!', 'XPAVG ("pjc_pfix_mod.F90")' )
    ENDIF

    ! Take avg of entire P array
    AVG  = SUM( P ) / DBLE( IM )

    ! Store average value in all elements of P
    P(:) = AVG

  END SUBROUTINE XPAVG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pjc_pfix
!
! !DESCRIPTION: INIT\_PJC\_PFIX_WINDOW allocates and initializes module arrays
!  and variables.  GMI dimension variables will be used for compatibility with
!  the Phil Cameron-Smith P-fixer. (bdf, bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_PJC_PFIX_WINDOW( State_Grid )
!
! !USES:
!
    USE ERROR_MOD,    ONLY   : ALLOC_ERR,   ERROR_STOP
    USE PhysConstants        ! Re, PI, etc...
    USE PRESSURE_MOD, ONLY   : GET_AP,      GET_BP
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN ) :: State_Grid  ! Grid State object
!
! !AUTHOR:
!   Brendan Field and Bob Yantosca (5/8/03)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS, I, J, L

    !=================================================================
    ! INIT_PJC_PFIX_WINDOW begins here!
    !
    ! Initialize dimensions for GMI pressure-fixer code
    !=================================================================
    IMP_NBORDER = 0
    I1_GL       = 1
    I2_GL       = State_Grid%NX
    JU1_GL      = 1
    JV1_GL      = 1
    J2_GL       = State_Grid%NY
    K1_GL       = 1
    K2_GL       = State_Grid%NZ
    ILO_GL      = I1_GL  - IMP_NBORDER
    IHI_GL      = I2_GL  + IMP_NBORDER
    JULO_GL     = JU1_GL - IMP_NBORDER
    JVLO_GL     = JV1_GL - IMP_NBORDER
    JHI_GL      = J2_GL  + IMP_NBORDER
    I1          = I1_GL
    I2          = I2_GL
    JU1         = JU1_GL
    JV1         = JV1_GL
    J2          = J2_GL
    K1          = K1_GL
    K2          = K2_GL
    ILO         = ILO_GL
    IHI         = IHI_GL
    JULO        = JULO_GL
    JVLO        = JVLO_GL
    JHI         = JHI_GL
    ! No polar cap. (ccc, 8/3/10)
    !J1P         = 3
    J1P         = 1
    J2P         = J2_GL - J1P + 1
    ! Used only to check dimensions
    ILAT        = J2_GL - JU1_GL + 1
    ILONG       = I2_GL -  I1_GL + 1
    IVERT       = K2_GL -  K1_GL + 1

    ! To add a buffer zone to calculate p-fixer for nested grid
    ! simulations. The p-fixer is not calculated for the edge boxes.
    ! (lzh, ccc, 8/3/10)
    BUFF_SIZE     = 2
    I1_W          = I1_GL + BUFF_SIZE
    I2_W          = I2_GL - BUFF_SIZE
    JU1_W         = JU1_GL + BUFF_SIZE
    J2_W          = J2_GL - BUFF_SIZE
    J1P_W         = 1 + BUFF_SIZE
    J2P_W         = J2_GL - J1P_W + 1

    ! Error check longitude
    IF ( ILONG /= State_Grid%NX ) THEN
       CALL ERROR_STOP( 'Invalid longitude dimension ILONG!', &
                        'INIT_PJC_FIX_WINDOW' )
    ENDIF

    ! Error check latitude
    IF ( ILAT /= State_Grid%NY ) THEN
       CALL ERROR_STOP( 'Invalid latitude dimension ILAT!', &
                        'INIT_PJC_FIX_WINDOW' )
    ENDIF

    ! Error check altitude
    IF ( IVERT /= State_Grid%NZ ) THEN
       CALL ERROR_STOP( 'Invalid altitude dimension IVERT!', &
                        'INIT_PJC_FIX_WINDOW' )
    ENDIF

    !=================================================================
    ! Allocate module arrays (use dimensions from GMI code)
    !=================================================================
    ALLOCATE( AI( K1_GL-1:K2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AI' )

    ALLOCATE( BI( K1_GL-1:K2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'BI' )

    ALLOCATE( DAP( K1_GL:K2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'DAP' )

    ALLOCATE( DBK( K1_GL:K2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'DBK' )

    ALLOCATE( CLAT_FV( JU1_GL:J2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLAT_FV' )

    ALLOCATE( COSE_FV( JU1_GL:J2_GL+1 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSE_FV' )

    ALLOCATE( COSP_FV( JU1_GL:J2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSP_FV' )

    ALLOCATE( DLAT_FV( JU1_GL:J2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'DLAT_FV' )

    ALLOCATE( ELAT_FV( JU1_GL:J2_GL+1 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'ELAT_FV' )

    ALLOCATE( GEOFAC( JU1_GL:J2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'GEOFAC' )

    ALLOCATE( GW_FV( JU1_GL:J2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'GW_FV' )

    ALLOCATE( MCOR( I1_GL:I2_GL, JU1_GL:J2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'MCOR' )

    ALLOCATE( REL_AREA( I1_GL:I2_GL, JU1_GL:J2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'REL_AREA' )

    ALLOCATE( RGW_FV( JU1_GL:J2_GL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'RGW_FV' )

    ALLOCATE( SINE_FV( JU1_GL:J2_GL+1 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SINE_FV' )

    !=================================================================
    ! Initialize arrays and variables
    !=================================================================

    ! Grid box surface areas [m2]
    DO J = JU1_GL, J2_GL
    DO I =  I1_GL, I2_GL
       MCOR(I,J) = State_Grid%Area_M2(I,J)
    ENDDO
    ENDDO

    ! Hybrid grid vertical coords: Ai [hPa] and Bi [unitless]
    DO L = K1_GL-1, K2_GL
       AI(L) = GET_AP( L+1 )
       BI(L) = GET_BP( L+1 )
    ENDDO

    ! Delta A [hPa] and Delta B [unitless]
    DO L = K1_GL, K2_GL
       !-------------------------------------------------------------
       ! NOTE:, this was the original code.  But since AI is already
       ! in hPa, we shouldn't need to multiply by PTOP again.  This
       ! should only matter for the fvDAS fields.  Also, DBK needs
       ! to be positive (bmy, 5/8/03)
       !DAP(L) = ( AI(L) - AI(L-1) ) * PTOP
       !DBK(L) = BI(L) - BI(L-1)
       !-------------------------------------------------------------
       DAP(L) = AI(L-1) - AI(L)
       DBK(L) = BI(L-1) - BI(L)
    ENDDO

    ! Grid box center latitudes [radians]
    DO J = JU1_GL, J2_GL
       CLAT_FV(J) = State_Grid%YMid_R(1,J)
    ENDDO

    ! Longitude spacing
    IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
       DLON_FV    = 2.0_fp * PI / 1152.0_fp
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
       DLON_FV    = 2.0_fp * PI / 576.0_fp
    ENDIF
    
    ! Latitude edge at south pole [radians]
    !ELAT_FV(1) = -0.5e+0_fp * PI
    IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
       ELAT_FV(1) = CLAT_FV(1) - 0.125_fp * PI / 180.0_fp
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
       ELAT_FV(1) = CLAT_FV(1) - 0.25_fp  * PI / 180.0_fp
    ENDIF

    ! SIN and COS of lat edge at south pole [unitless]
    !SINE_FV(1) = -1.e+0_fp
    !COSE_FV(1) =  0.e+0_fp
    ! for nested NA or EA (lzh, 07/20/2010)
    SINE_FV(1) =  SIN( ELAT_FV(1) )
    COSE_FV(1) =  COS( ELAT_FV(1) )

    ! Latitude edges [radians] (w/ SIN & COS) at intermediate latitudes
    DO J = JU1_GL+1, J2_GL  !2, State_Grid%NY
       ELAT_FV(J) = 0.5e+0_fp * ( CLAT_FV(J-1) + CLAT_FV(J) )
       SINE_FV(J) = SIN( ELAT_FV(J) )
       COSE_FV(J) = COS( ELAT_FV(J) )
    ENDDO

    ! Latitude edge at North Pole [radians]
    !ELAT_FV(J2_GL+1) = 0.5e+0_fp * PI
    ! for nested NA or EA (lzh, 07/20/2010)
    IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
       ELAT_FV(J2_GL+1) = CLAT_FV(J2_GL)+0.125_fp * PI / 180.0_fp
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
       ELAT_FV(J2_GL+1) = CLAT_FV(J2_GL)+0.25_fp  * PI / 180.0_fp
    ENDIF

    ! SIN of lat edge at North Pole
    !SINE_FV(J2_GL+1) = 1.e+0_fp
    ! for nested NA or EA (lzh, 07/20/2010)
    SINE_FV(J2_GL+1) =  SIN( ELAT_FV(J2_GL+1) )
    COSE_FV(J2_GL+1) =  COS( ELAT_FV(J2_GL+1) )

    ! Latitude extent of South polar box [radians]
    !DLAT_FV(1) = 2.e+0_fp * ( ELAT_FV(2) - ELAT_FV(1) )
    ! comment out for nested NA or EA (lzh, 07/20/2010)

    ! Latitude extent of boxes at intermediate latitudes [radians]
    !DO J = JU1_GL+1, J2_GL-1  ! 2, State_Grid%NY-1
    ! for nested NA or EA (lzh, 07/20/2010)
    DO J = JU1_GL, J2_GL
       DLAT_FV(J) = ELAT_FV(J+1) - ELAT_FV(J)
    ENDDO

    ! Latitude extent of North polar box [radians]
    !DLAT_FV(J2_GL) = 2.e+0_fp * ( ELAT_FV(J2_GL+1) - ELAT_FV(J2_GL) )
    ! comment out for nested NA or EA (lzh, 07/20/2010)

    ! Other stuff
    DO J = JU1_GL, J2_GL
       GW_FV(J)   = SINE_FV(J+1) - SINE_FV(J)
       COSP_FV(J) = GW_FV(J)     / DLAT_FV(J)
       RGW_FV(J)  = 1.e+0_fp         / GW_FV(J)
    ENDDO

  END SUBROUTINE INIT_PJC_PFIX_WINDOW
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_pjc_pfix_window
!
! !DESCRIPTION: Subroutine CLEANUP\_PJC\_PFIX\_WINDOW deallocates all module
!  arrays (bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_PJC_PFIX_WINDOW
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( ALLOCATED( AI       ) ) DEALLOCATE( AI       )
    IF ( ALLOCATED( BI       ) ) DEALLOCATE( BI       )
    IF ( ALLOCATED( CLAT_FV  ) ) DEALLOCATE( CLAT_FV  )
    IF ( ALLOCATED( COSE_FV  ) ) DEALLOCATE( COSE_FV  )
    IF ( ALLOCATED( COSP_FV  ) ) DEALLOCATE( COSP_FV  )
    IF ( ALLOCATED( DAP      ) ) DEALLOCATE( DAP      )
    IF ( ALLOCATED( DBK      ) ) DEALLOCATE( DBK      )
    IF ( ALLOCATED( DLAT_FV  ) ) DEALLOCATE( DLAT_FV  )
    IF ( ALLOCATED( ELAT_FV  ) ) DEALLOCATE( ELAT_FV  )
    IF ( ALLOCATED( GEOFAC   ) ) DEALLOCATE( GEOFAC   )
    IF ( ALLOCATED( GW_FV    ) ) DEALLOCATE( GW_FV    )
    IF ( ALLOCATED( MCOR     ) ) DEALLOCATE( MCOR     )
    IF ( ALLOCATED( REL_AREA ) ) DEALLOCATE( REL_AREA )
    IF ( ALLOCATED( RGW_FV   ) ) DEALLOCATE( RGW_FV   )
    IF ( ALLOCATED( SINE_FV  ) ) DEALLOCATE( SINE_FV )

  END SUBROUTINE CLEANUP_PJC_PFIX_WINDOW
!EOC
END MODULE PJC_PFIX_WINDOW_MOD
