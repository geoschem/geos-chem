!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pjc_pfix_mod.F90
!
! !DESCRIPTION: Module Pjc\_Pfix\_Mod contains routines which implements the
!  Philip Cameron-Smith pressure fixer for the new fvDAS transport
!  scheme. (bdf, bmy, 5/8/03, 10/27/03)
!\\
!\\
! !INTERFACE:
!
MODULE PJC_PFIX_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Do_Pjc_Pfix
  PUBLIC  :: Cleanup_Pjc_Pfix
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Calc_Pressure
  PRIVATE :: Calc_Advection_Factors
  PRIVATE :: Adjust_Press
  PRIVATE :: Init_Press_Fix
  PRIVATE :: Do_Press_Fix_LLNL
  PRIVATE :: Average_Press_Poles
  PRIVATE :: Convert_Winds
  PRIVATE :: Calc_Horiz_Mass_Flux
  PRIVATE :: Calc_Divergence
  PRIVATE :: Set_Press_Terms
  PRIVATE :: Do_Divergence_Pole_Sum
  PRIVATE :: Xpavg
  PRIVATE :: Init_Pjc_Pfix
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!  Brendan Field and Bob Yantosca (5/8/03)
!  Modified for new GMI TPCORE by Claire Carouge (ccarouge@seas.harvard.edu)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
!
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
  LOGICAL             :: PR_DIAG
  INTEGER             :: LOC_PROC
  REAL(fp)            :: GEOFAC_PC
  REAL(fp)            :: DLON_FV

  ! Dimensions for GMI code (from "imp_dims")
  INTEGER             :: IMP_NBORDER
  INTEGER             :: I1_GL,  I2_GL,   JU1_GL,  JV1_GL
  INTEGER             :: J2_GL,  K1_GL,   K2_GL,   ILO_GL
  INTEGER             :: IHI_GL, JULO_GL, JVLO_GL, JHI_GL
  INTEGER             :: I1,     I2,      JU1,     JV1
  INTEGER             :: J2,     K1,      K2,      ILO
  INTEGER             :: IHI,    JULO,    JVLO,    JHI
  INTEGER             :: ILAT,   ILONG,   IVERT,   J1P
  INTEGER             :: J2P

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
! !IROUTINE: Do_Pjc_Pfix
!
! !DESCRIPTION: Subroutine Do\_Pjc\_Pfix is the driver routine for the Philip
!  Cameron-Smith pressure fixer for the fvDAS transport scheme.
!  (bdf, bmy, 5/8/03, 3/5/07)
!\\
!\\
!  We assume that the winds are on the A-GRID, since this is the input that
!  the fvDAS transport scheme takes. (bdf, bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Pjc_Pfix( State_Grid, D_DYN, P1, P2, UWND, VWND, XMASS, YMASS )
!
! !USES:
!
    USE PhysConstants   ! Physical constants
    Use State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    REAL(fp),       INTENT(IN)  :: D_DYN       ! Dynamic timestep [s]
    REAL(fp),       INTENT(IN)  :: P1(:,:)     ! True PSurface at middle of
                                               ! dynamic timestep [hPa]
    REAL(fp),       INTENT(IN)  :: P2(:,:)     ! True PSurface at end of
                                               ! dynamic timestep [hPa]
    REAL(fp),       INTENT(IN)  :: UWND(State_Grid%NX, & ! Zonal wind [m/s]
                                        State_Grid%NY, &
                                        State_Grid%NZ)
    REAL(fp),       INTENT(IN)  :: VWND(State_Grid%NX, & ! Meridional wind [m/s]
                                        State_Grid%NY, &
                                        State_Grid%NZ)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: XMASS(State_Grid%NX, & ! E-W mass fluxes
                                         State_Grid%NY, & !  [mixing ratio]
                                         State_Grid%NZ)

    REAL(fp),       INTENT(OUT) :: YMASS(State_Grid%NX, & ! N-S mass fluxes
                                         State_Grid%NY, & !  [mixing ratio]
                                         State_Grid%NZ)
!
! !AUTHOR:
!  Brendan Field and Bob Yantosca (5/8/03)
!
! !REMARKS:
!  (1 ) Now P1 and P2 are "true" surface pressures, and not PS-PTOP.  If using
!        this P-fixer w/ GEOS-3 winds, pass true surface pressure to this
!        routine. (bmy, 10/27/03)
!  (2 ) Now define P2_TMP array for passing to ADJUST_PRESS (yxw, bmy, 3/5/07)
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
    INTEGER              :: I, J, K
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
       CALL INIT_PJC_PFIX( State_Grid )

       ! Calculate advection surface-area factors
       CALL CALC_ADVECTION_FACTORS( MCOR, REL_AREA, GEOFAC, GEOFAC_PC)

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    ! Copy P2 into P2_TMP (yxw, bmy, 3/5/07)
    P2_TMP = P2

    ! Call PJC pressure fixer w/ the proper arguments
    ! NOTE: P1 and P2 are now "true" surface pressure, not PS-PTOP!!!
    CALL ADJUST_PRESS( 'GEOS-CHEM',        INTERP_WINDS,  &
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

  END SUBROUTINE Do_Pjc_Pfix
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
    TYPE(GrdState), INTENT(IN)  :: State_Grid      ! Grid State object
    REAL(fp),       INTENT(IN)  :: XMASS(State_Grid%NX, & ! E-W mass flux from
                                         State_Grid%NY, & !  pressure fixer
                                         State_Grid%NZ)
    REAL(fp),       INTENT(IN)  :: YMASS(State_Grid%NX, & ! N-S mass flux from
                                         State_Grid%NY, & !  pressure fixer
                                         State_Grid%NZ)
    REAL(fp),       INTENT(IN)  :: PS_NOW(State_Grid%NX, & ! Sfc pressure - PTOP
                                          State_Grid%NY)   !  at current time
    REAL(fp),       INTENT(IN)  :: RGW_FV(State_Grid%NY)   ! Latitude factor
                                                        ! 1/(SINE(J+1) - SIN(J))
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
    INTEGER :: ij
    REAL(fp)  :: dp           ! spacing in latitude (rad)
    REAL(fp)  :: ri2_gl
    REAL(fp)  :: rj2m1
    REAL(fp)  :: total_area

    !----------------
    !Begin execution.
    !----------------

    ri2_gl = i2_gl

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

    rj2m1 = j2_gl - 1
    dp    = PI / rj2m1

    do ij = ju1_gl, j2_gl
       geofac(ij) = dp / (2.0e+0_fp * rel_area(1,ij) * ri2_gl)
    end do

    geofac_pc = dp / (2.0e+0_fp * Sum (rel_area(1,ju1_gl:ju1_gl+1)) * ri2_gl)

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
  SUBROUTINE Adjust_Press &
        (metdata_name_org, do_timinterp_winds, new_met_rec,         &
         met_grid_type, advec_consrv_opt, pmet2_opt, press_fix_opt, &
         tdt, geofac_pc, geofac, cose, cosp, rel_area, dap, dbk,    &
         pctm1, pctm2, pmet2, uu, vv, xmass, ymass)
!
! !INPUT PARAMETERS:
!
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
    !dps : change of surface pressure from met field pressure [hPa]
    !-------------------------------------------------------------
    real(fp)  :: dps(i1_gl:i2_gl, ju1_gl:j2_gl)

    !--------------------------------------------
    !dps_ctm : CTM surface pressure tendency [hPa]
    !--------------------------------------------
    real(fp) :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)

    !---------------------------------------------------------------------
    !xmass_fixed : horizontal mass flux in E-W direction after fixing [hPa]
    !ymass_fixed : horizontal mass flux in N-S direction after fixing [hPa]
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

    dps_ctm(:,:) = 0.0e+0_fp

    dgpress =  Sum ( (pmet2(i1_gl:i2_gl, ju1_gl:j2_gl) - &
                      pctm1(i1_gl:i2_gl, ju1_gl:j2_gl)   ) &
                 * rel_area(i1_gl:i2_gl, ju1_gl:j2_gl)     )

    if (pmet2_opt == 1) then
       pmet2(:,:) = pmet2(:,:) - dgpress
    end if

    !### Debug
    !###if (DO_ADJUST_PRESS_DIAG) then
    !###  Write (6, *) 'Global mean surface pressure change [hPa] = ',
    !###                dgpress
    !###end if

    !===================
    call Init_Press_Fix &
    !===================
         (metdata_name_org, met_grid_type, tdt, geofac_pc, geofac, &
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
       !###Write (6, *) 'RMS deviation between pmet2 & pctm2 [hPa] = ',
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
  SUBROUTINE Init_Press_Fix &
       (metdata_name_org, met_grid_type, tdt, geofac_pc, geofac, &
       cose, cosp, dap, dbk, dps, dps_ctm, rel_area, pctm1, pmet2, &
       uu, vv, xmass, ymass)
!
! !INPUT PARAMETERS:
!
    ! First part of metdata_name, e.g., "NCAR"
    CHARACTER(LEN=*) :: metdata_name_org

    ! Met grid type, A or C
    INTEGER          :: met_grid_type

    ! Model Time step [s]
    REAL(fp) :: tdt

    ! Special geometrical factor (geofac) for Polar cap
    REAL(fp)         :: geofac_pc

    ! Cosine of grid box edges and centers
    REAL(fp)         :: cose(ju1_gl:j2_gl)
    REAL(fp)         :: cosp(ju1_gl:j2_gl)

    ! Geometrical factor for meridional advection; geofac uses
    ! correct spherical geometry, and replaces acosp as the
    ! meridional geometrical factor in tpcore
    REAL(fp)         :: geofac(ju1_gl:j2_gl)

    ! Pressure difference across layer from (ai * pt) term [hPa]
    REAL(fp)         :: dap(k1:k2)

    ! Difference in bi across layer - the dSigma term
    REAL(fp)         :: dbk(k1:k2)

    ! relative surface area of grid box (fraction)
    REAL(fp)         :: rel_area( i1_gl:i2_gl, ju1_gl:j2_gl)

    ! Metfield surface pressure at t1 [hPa]
    REAL(fp)         :: pmet2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

    ! CTM surface pressure at t1 [hPa]
    REAL(fp)         :: pctm1(ilo_gl:ihi_gl, julo_gl:jhi_gl)

    ! CTM surface pressure at t1+tdt [hPa]
    REAL(fp)         :: pctm2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

    ! Wind velocity, x direction at t1+tdt/2 [m/s]
    REAL(fp)         :: uu(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

    ! Wind velocity, y direction at t1+tdt/2 [m/s]
    REAL(fp)         :: vv(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
!
! !OUTPUT PARAMETERS:
!
    ! Horizontal mass flux in E-W direction [hPa]
    REAL(fp)  :: xmass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

    ! Horizontal mass flux in N-S direction [hPa]
    REAL(fp)  :: ymass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

    ! Change of surface pressure from met field pressure [hPa]
    REAL(fp)  :: dps(i1_gl:i2_gl, ju1_gl:j2_gl)

    ! CTM surface pressure tendency [hPa]
    REAL(fp)  :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)
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
    !        motion [hPa]
    !--------------------------------------------------------------
    real(fp)  :: dpi(i1:i2, ju1:j2, k1:k2)

    !---------------------------------------------------------------------
    !crx   : Courant number in E-W direction
    !cry   : Courant number in N-S direction
    !delp1 : pressure thickness, the psudo-density in a hydrostatic system
    !        at t1 [hPa]
    !delpm : pressure thickness, the psudo-density in a hydrostatic system
    !        at t1+tdt/2 (approximate) [hPa]
    !pu    : pressure at edges in "u"  [hPa]
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

    !========================
    call Average_Press_Poles &
    !========================
         (rel_area, pctm1)

    !========================
    call Average_Press_Poles &
    !========================
         (rel_area, pmet2)

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
         (met_grid_type, tdt, cosp, crx, cry, uu, vv)

    !=========================
    call Calc_Horiz_Mass_Flux &
    !=========================
         (cose, delpm, uu, vv, xmass, ymass, tdt, cosp)

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
    REAL(fp), INTENT(IN)   :: geofac_pc

    ! Geometrical factor for meridional advection; geofac uses
    ! correct spherical geometry, and replaces acosp as the
    !  meridional geometrical factor in tpcore
    REAL(fp), INTENT(IN)   :: geofac(ju1_gl:j2_gl)

    ! Difference in bi across layer - the dSigma term
    REAL(fp), INTENT(IN)   :: dbk(k1:k2)

    ! Change of surface pressure from met field pressure [hPa]
    REAL(fp), INTENT(IN)   :: dps(i1:i2, ju1:j2)

    ! Relative surface area of grid box (fraction)
    REAL(fp), INTENT(IN)   :: rel_area(i1:i2, ju1:j2)

    ! Horizontal mass fluxes in E-W and N-S directions [hPa]
    REAL(fp), INTENT(IN)   :: xmass(ilo:ihi, julo:jhi, k1:k2)
    REAL(fp), INTENT(IN)   :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
    ! Sum over vertical of dpi calculated from original mass fluxes [hPa]
    REAL(fp),  INTENT(OUT) :: dps_ctm(i1:i2, ju1:j2)

    ! Horizontal mass flux in E-W and N-S directions after fixing [hPa]
    REAL(fp),  INTENT(OUT) :: xmass_fixed(ilo:ihi, julo:jhi, k1:k2)
    REAL(fp),  INTENT(OUT) :: ymass_fixed(ilo:ihi, julo:jhi, k1:k2)
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
    !dpi: divergence at a grid point; used to calculate vertical motion [hPa]
    !------------------------------------------------------------------------
    real(fp)  :: dpi(i1:i2, ju1:j2, k1:k2)

    real(fp)  :: xcolmass_fix(ilo:ihi, julo:jhi)

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6,*) 'Do_Press_Fix_Llnl called by ', loc_proc
    end if

    ri2 = i2_gl

    mmfd(:) = 0.0e+0_fp

    xcolmass_fix(:,:)   = 0.0e+0_fp

    xmass_fixed (:,:,:) = xmass(:,:,:)
    ymass_fixed (:,:,:) = ymass(:,:,:)

    !------------------------------------------------------------
    !Calculate difference between GCM and LR predicted pressures.
    !------------------------------------------------------------

    ddps(:,:) = dps(:,:) - dps_ctm(:,:)

    !--------------------------------------
    ! Calculate global-pressure discrepancy.
    !--------------------------------------

    dgpress = Sum (ddps(i1:i2,ju1:j2) * rel_area(i1:i2,ju1:j2))

    !----------------------------------------------------------
    !Calculate mean meridional flux divergence (df/dy).
    !Note that mmfd is actually the zonal mean pressure change,
    !which is related to df/dy by geometrical factors.
    !----------------------------------------------------------

    !------------------------
    !Handle non-Pole regions.
    !------------------------

    do ij = j1p, j2p
       mmfd(ij) = -(sum(ddps(:,ij)) / ri2 - dgpress)
    end do

    !---------------------------------------------
    !Handle poles.
    !Note that polar boxes have all been averaged.
    !---------------------------------------------

    mmfd(ju1)   = -(ddps(1,ju1)   - dgpress)
    mmfd(ju1+1) = -(ddps(1,ju1+1) - dgpress)
    mmfd(j2-1)  = -(ddps(1,j2-1)  - dgpress)
    mmfd(j2)    = -(ddps(1,j2)    - dgpress)

    !---------------------------------------------
    !Calculate mean meridional fluxes (cos(e)*fy).
    !---------------------------------------------

    mmf(j1p) = mmfd(ju1) / geofac_pc

    do ij = j1p, j2p
       mmf(ij+1) = mmf(ij) + mmfd(ij) / geofac(ij)
    end do

    !------------------------------------------------------------
    !Fix latitude bands.
    !Note that we don't need to worry about geometry here because
    !all boxes in a latitude band are identical.
    !Note also that fxintegral(i2+1) should equal fxintegral(i1),
    !i.e., zero.
    !------------------------------------------------------------
    do ij = j1p, j2p

       fxintegral(:) = 0.0e+0_fp

       do il = i1, i2
          fxintegral(il+1) = fxintegral(il) - &
                             (ddps(il,ij) - dgpress) - mmfd(ij)
       end do

       fxmean = Sum (fxintegral(i1+1:i2+1)) / ri2

       do il = i1, i2
          xcolmass_fix(il,ij) = fxintegral(il) - fxmean
       end do

    end do

    !-------------------------------------
    !Distribute colmass_fix's in vertical.
    !-------------------------------------

    do ik = k1, k2
    do ij = j1p, j2p
    do il = i1, i2

       xmass_fixed(il,ij,ik) = xmass(il,ij,ik) + xcolmass_fix(il,ij) * dbk(ik)

    end do
    end do
    end do

    do ik = k1, k2
    do ij = j1p, j2p+1
    do il = i1, i2

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
! !IROUTINE: Average_Press_Poles
!
! !DESCRIPTION: Subroutine Average\_Press\_Poles averages pressure at the
!  Poles when the Polar cap is enlarged.  It makes the last two latitudes
!  equal.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Average_Press_Poles(rel_area, press)
!
! !INPUT PARAMETERS:
!
    ! Relative surface area of grid box (fraction)
    REAL(fp), INTENT(IN)    :: rel_area(i1:i2, ju1:j2)
!
! !OUTPUT PARAMETERS:
!
    ! Surface pressure [hPa]
    REAL(fp), INTENT(INOUT) :: press   (ilo:ihi, julo:jhi)
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)  :: meanp

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6,*) 'Average_Press_Poles called by ', loc_proc
    end if

    meanp = Sum (rel_area(i1:i2,ju1:ju1+1) * press(i1:i2,ju1:ju1+1)) / &
            Sum (rel_area(i1:i2,ju1:ju1+1))

    press(i1:i2,ju1:ju1+1) = meanp

    meanp = Sum (rel_area(i1:i2,j2-1:j2) * press(i1:i2,j2-1:j2)) / &
            Sum (rel_area(i1:i2,j2-1:j2))

    press(i1:i2,j2-1:j2) = meanp

  END SUBROUTINE Average_Press_Poles
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
  SUBROUTINE Convert_Winds(igd, tdt, cosp, crx, cry, uu, vv)
!
! !USES:
!
    USE PhysConstants ! Re, PI
!
! !INPUT PARAMETERS:
!
    ! A or C grid
    INTEGER, INTENT(IN)  :: igd

    ! Model time step [s]
    REAL(fp),  INTENT(IN)  :: tdt

    ! Cosine of grid box centers
    REAL(fp),  INTENT(IN)  :: cosp(ju1_gl:j2_gl)

    ! Wind velocity in E-W (UU) and N-S (VV) directions at t1+tdt/2 [m/s]
    REAL(fp),  INTENT(IN)  :: uu  (ilo:ihi, julo:jhi, k1:k2)
    REAL(fp),  INTENT(IN)  :: vv  (ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
    ! Courant number in E-W (CRX) and N-S (CRY) directions
    REAL(fp),  INTENT(OUT) :: crx (ilo:ihi, julo:jhi, k1:k2)
    REAL(fp),  INTENT(OUT) :: cry (ilo:ihi, julo:jhi, k1:k2)
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
       dtdx = 0.0e+0_fp; dtdx5 = 0.0e+0_fp

       ri2   = i2_gl
       rj2m1 = j2_gl - 1

       dl    = 2.0e+0_fp * PI / ri2
       dp    = PI / rj2m1

       dtdy  = tdt / (Re * dp)
       dtdy5 = 0.5e+0_fp * dtdy


       dtdx (ju1_gl) = 0.0e+0_fp
       dtdx5(ju1_gl) = 0.0e+0_fp

       do ij = ju1_gl + 1, j2_gl - 1

          dtdx (ij) = tdt / (dl * Re * cosp(ij))
          dtdx5(ij) = 0.5e+0_fp * dtdx(ij)

       end do

       dtdx (j2_gl)  = 0.0e+0_fp
       dtdx5(j2_gl)  = 0.0e+0_fp

    end if

    !=============
    if (igd == 0) then  ! A grid.
    !=============

       do ij = ju1+1, j2-1
          do il = i1+1, i2
             crx(il,ij,:) = dtdx5(ij) * &
                            (uu(il,ij,:) + uu(il-1,ij,  :))
          end do
          crx(1,ij,:) = dtdx5(ij) * &
                        (uu(1,ij,:) + uu(i2,ij,  :))
       end do

       do ij = ju1+1, j2
          do il = i1, i2
             cry(il,ij,:) = dtdy5 * (vv(il,ij,:) + vv(il,  ij-1,:))
          end do
       end do

    !====
    else  ! C grid.
    !====

       do ij = ju1, j2
          do il = i1, i2

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
       (cose, delpm, uu, vv, xmass, ymass, tdt, cosp)
!
! !USES:
!
    USE PhysConstants ! Re, Pi
!
! !INPUT PARAMETERS:
!
    ! Timestep [s]
    REAL(fp), INTENT(IN)   :: tdt

    ! Cosine of grid box edges
    REAL(fp), INTENT(IN)   :: cose (ju1_gl:j2_gl)

    ! Cosine of grid box centers
    REAL(fp), INTENT(IN)   :: cosp (ju1_gl:j2_gl)

    ! Pressure thickness, the pseudo-density in a
    ! hdrostatic system  at t1+tdt/2 (approximate) [hPa]
    REAL(fp), INTENT(IN)   :: delpm(ilo:ihi, julo:jhi, k1:k2)

    ! E-W (UU) and N-S (VV) winds [m/s]
    REAL(fp), INTENT(IN)   :: uu  (ilo:ihi, julo:jhi, k1:k2)
    REAL(fp), INTENT(IN)   :: vv  (ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
    ! Horizontal mass flux in E-W and N-S directions [hPa]
    REAL(fp), INTENT(OUT)  :: xmass(ilo:ihi, julo:jhi, k1:k2)
    REAL(fp), INTENT(OUT)  :: ymass(ilo:ihi, julo:jhi, k1:k2)
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
    INTEGER   :: ij
    INTEGER   :: il
    INTEGER   :: jst, jend
    REAL(fp)  :: dl
    REAL(fp)  :: dp

    REAL(fp)  :: ri2
    REAL(fp)  :: rj2m1
    REAL(fp)  :: factx
    REAL(fp)  :: facty

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6,*) 'Calc_Horiz_Mass_Flux called by ', loc_proc
    end if

    ri2   = i2_gl
    rj2m1 = j2_gl - 1

    dl    = 2.0e+0_fp * PI / ri2
    dp    = PI / rj2m1

    facty  = 0.5e+0_fp * tdt / (Re * dp)

    !-----------------------------------
    !Calculate E-W horizontal mass flux.
    !-----------------------------------

    do ij = ju1, j2

       factx = 0.5e+0_fp * tdt / (dl * Re * cosp(ij))

       do il = i1+1, i2
          xmass(il,ij,:) = factx * &
                           (uu(il  ,ij,:) * delpm(il  ,ij,:) + &
                            uu(il-1,ij,:) * delpm(il-1,ij,:))
       end do

       xmass(i1,ij,:) = factx * &
                        (uu(i1,ij,:) * delpm(i1,ij,:) + &
                         uu(i2,ij,:) * delpm(i2,ij,:))

    end do

    !-----------------------------------
    !Calculate N-S horizontal mass flux.
    !-----------------------------------

    do ij = ju1+1, j2

       ymass(i1:i2,ij,:) = facty * &
            cose(ij) * (vv(i1:i2,ij,:)*delpm(i1:i2,ij,:) + &
            vv(i1:i2,ij-1,:)*delpm(i1:i2,ij-1,:))

    end do

    ymass(i1:i2,ju1,:) = facty * &
         cose(ju1) * (vv(i1:i2,ju1,:)*delpm(i1:i2,ju1,:))

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
    LOGICAL, INTENT(IN)    :: do_reduction

    ! Special geometrical factor (geofac) for Polar cap
    REAL(fp),  INTENT(IN)    :: geofac_pc

    ! geometrical factor for meridional advection; geofac uses
    ! correct spherical geometry, and replaces acosp as the
    ! meridional geometrical factor in tpcore
    REAL(fp),  INTENT(IN)    :: geofac(ju1_gl:j2_gl)

    ! horizontal mass fluxes in E-W and N-S directions [hPa]
    REAL(fp),  INTENT(IN)    :: xmass (ilo:ihi, julo:jhi, k1:k2)
    REAL(fp),  INTENT(IN)    :: ymass (ilo:ihi, julo:jhi, k1:k2)
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Divergence at a grid point; used to calculate vertical motion [hPa]
    REAL(fp),  INTENT(INOUT) :: dpi   (i1:i2, ju1:j2, k1:k2)
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
    integer :: jst, jend

    !----------------
    !Begin execution.
    !----------------

    if (pr_diag) then
       Write (6,*) 'Calc_Divergence called by ', loc_proc
    end if

    !-------------------------
    !Calculate N-S divergence.
    !-------------------------

    do ij = j1p, j2p

       dpi(i1:i2,ij,:) = (ymass(i1:i2,ij,:) - ymass(i1:i2,ij+1,:)) * geofac(ij)

    end do

    !-------------------------
    !Calculate E-W divergence.
    !-------------------------

    do ij = j1p,j2p
       do il = i1, i2-1
          dpi(il,ij,:) = dpi(il,ij,:) + xmass(il,ij,:) - xmass(il+1,ij,:)
       end do
       dpi(i2,ij,:) = dpi(i2,ij,:) + xmass(i2,ij,:) - xmass(1,ij,:)
    end do

    !===========================
    call Do_Divergence_Pole_Sum &
    !===========================
         (do_reduction, geofac_pc, dpi, ymass)

    ! Added this IF statemetn (ccarouge, 12/3/08)
    if (j1p /= ju1_gl+1) then

       !--------------------------------------------
       !Polar cap enlarged:  copy dpi to polar ring.
       !--------------------------------------------

       if (ju1 == ju1_gl) then

          dpi(:,ju1+1,:) = dpi(:,ju1,:)

       end if

       if (j2 == j2_gl) then

          dpi(:,j2-1,:)  = dpi(:,j2,:)

       end if

    end if

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

       pu(i1,ij,:) = 0.5e+0_fp * (delpm(i1,ij,:) + delpm(i2,ij,:))

    end do

  END SUBROUTINE Set_Press_Terms
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Divergence_Pole_Sum
!
! !DESCRIPTION: Do\_Divergence\_Pole\_Sum sets the divergence at the Poles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Divergence_Pole_Sum(do_reduction, geofac_pc, dpi, ymass)
!
! !INPUT PARAMETERS:
!
    ! Set to T if called on Master; set to F if called by Slaves
    ! (NOTE: This does not seem to be used!)
    LOGICAL :: do_reduction

    ! Special geometrical factor (geofac) for Polar cap
    REAL(fp)  :: geofac_pc

    ! horizontal mass flux in N-S direction [hPa]
    REAL(fp)  :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
    ! Divergence at a grid point; used to calculate vertical motion [hPa]
    REAL(fp)  :: dpi  ( i1:i2,   ju1:j2,  k1:k2)
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
    !----------------------
    !Variable declarations.
    !----------------------

    integer :: il, ik

    real(fp)  :: ri2

    real(fp)  :: mean_np(k1:k2)
    real(fp)  :: mean_sp(k1:k2)
    real(fp)  :: sumnp  (k1:k2)
    real(fp)  :: sumsp  (k1:k2)

    !----------------
    !Begin execution.
    !----------------

    ri2 = i2_gl

    !==================
    if (ju1 == ju1_gl) then
    !==================

       do ik = k1, k2

          sumsp(ik) = 0.0e+0_fp

          do il = i1, i2

             sumsp(ik) = sumsp(ik) + ymass(il,j1p,ik)

          end do

       end do

       do ik = k1, k2

          mean_sp(ik) = -sumsp(ik) / ri2 * geofac_pc

          do il = i1, i2

             dpi(il,ju1,ik) = mean_sp(ik)

          end do

       end do

    !======
    end if
    !======

    !================
    if (j2 == j2_gl) then
    !================

       do ik = k1, k2

          sumnp(ik) = 0.0e+0_fp

          do il = i1, i2

             sumnp(ik) = sumnp(ik) + ymass(il,j2p+1,ik)

          end do

       end do

       do ik = k1, k2

          mean_np(ik) = sumnp(ik) / ri2 * geofac_pc

          do il = i1, i2

             dpi(il,j2,ik) = mean_np(ik)

          end do

       end do

    !======
    end if
    !======

  END SUBROUTINE Do_Divergence_Pole_Sum
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
  SUBROUTINE Xpavg( P, IM )
!
! !USES:
!
    USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
    ! Dimension of P
    INTEGER, INTENT(IN)    :: IM
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! 1-D vector to be averaged
    REAL(fp),  INTENT(INOUT) :: P(IM)
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
    REAL(fp)                 :: AVG

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

  END SUBROUTINE Xpavg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Pjc_Pfix
!
! !DESCRIPTION: Subroutine Init\_Pjc\_Pfix allocates and initializes module
!  arrays and variables.  GMI dimension variables will be used for
!  compatibility with the Phil Cameron-Smith P-fixer. (bdf, bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Pjc_Pfix( State_Grid )
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ALLOC_ERR,   ERROR_STOP
    USE PRESSURE_MOD,   ONLY : GET_AP,      GET_BP
    USE PhysConstants        ! Re, PI, etc...
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
    ! INIT_PJC_PFIX begins here!
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
    ILAT        = J2_GL - JU1_GL + 1
    ILONG       = I2_GL -  I1_GL + 1
    IVERT       = K2_GL -  K1_GL + 1
    J1P         = 3
    J2P         = J2_GL - J1P + 1

    ! Error check longitude
    IF ( ILONG /= State_Grid%NX ) THEN
       CALL ERROR_STOP( 'Invalid longitude dimension ILONG!', &
                        'INIT_PJC_FIX ("pjc_pfix_mod.F90")' )
    ENDIF

    ! Error check latitude
    IF ( ILAT /= State_Grid%NY ) THEN
       CALL ERROR_STOP( 'Invalid latitude dimension ILAT!', &
                        'INIT_PJC_FIX ("pjc_pfix_mod.F90")' )
    ENDIF

    ! Error check altitude
    IF ( IVERT /= State_Grid%NZ ) THEN
       CALL ERROR_STOP( 'Invalid altitude dimension IVERT!', &
                        'INIT_PJC_FIX ("pjc_pfix_mod.F90")' )
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
    DLON_FV    = 2.e+0_fp * PI / DBLE( I2_GL )

    ! Latitude edge at south pole [radians]
    ELAT_FV(1) = -0.5e+0_fp * PI

    ! SIN and COS of lat edge at south pole [unitless]
    SINE_FV(1) = -1.e+0_fp
    COSE_FV(1) =  0.e+0_fp

    ! Latitude edges [radians] (w/ SIN & COS) at intermediate latitudes
    DO J = JU1_GL+1, J2_GL  !2, State_Grid%NY
       ELAT_FV(J) = 0.5e+0_fp * ( CLAT_FV(J-1) + CLAT_FV(J) )
       SINE_FV(J) = SIN( ELAT_FV(J) )
       COSE_FV(J) = COS( ELAT_FV(J) )
    ENDDO

    ! Latitude edge at North Pole [radians]
    ELAT_FV(J2_GL+1) = 0.5e+0_fp * PI

    ! SIN of lat edge at North Pole
    SINE_FV(J2_GL+1) = 1.e+0_fp

    ! Latitude extent of South polar box [radians]
    DLAT_FV(1) = 2.e+0_fp * ( ELAT_FV(2) - ELAT_FV(1) )

    ! Latitude extent of boxes at intermediate latitudes [radians]
    DO J = JU1_GL+1, J2_GL-1  ! 2, State_Grid%NY-1
       DLAT_FV(J) = ELAT_FV(J+1) - ELAT_FV(J)
    ENDDO

    ! Latitude extent of North polar box [radians]
    DLAT_FV(J2_GL) = 2.e+0_fp * ( ELAT_FV(J2_GL+1) - ELAT_FV(J2_GL) )

    ! Other stuff
    DO J = JU1_GL, J2_GL
       GW_FV(J)   = SINE_FV(J+1) - SINE_FV(J)
       COSP_FV(J) = GW_FV(J)     / DLAT_FV(J)
       RGW_FV(J)  = 1.e+0_fp         / GW_FV(J)
    ENDDO

  END SUBROUTINE Init_Pjc_Pfix
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Pjc_Pfix
!
! !DESCRIPTION: Subroutine Cleanup\_Pjc\_Pfix deallocates all module arrays
!  (bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Pjc_Pfix
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

  END SUBROUTINE Cleanup_Pjc_Pfix
!EOC
END MODULE Pjc_Pfix_Mod
