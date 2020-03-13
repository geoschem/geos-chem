!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: calc_met_mod.F90
!
! !DESCRIPTION: Module CALC\_MET\_MOD (formerly DAO\_MOD) contains
!  subroutines that compute, interpolate, or otherwise process met field data.
!\\
!\\
! !INTERFACE:
!
MODULE CALC_MET_MOD
!
! !USES:
!
  USE PhysConstants          ! Physical constants
  USE PRECISION_MOD          ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: AVGPOLE
  PUBLIC  :: AIRQNT
  PUBLIC  :: GET_COSINE_SZA
  PUBLIC  :: GET_OBK
  PUBLIC  :: INTERP
  PUBLIC  :: SET_DRY_SURFACE_PRESSURE
  PUBLIC  :: Set_Met_AgeOfAir
#if defined( ESMF_ ) || defined( EXTERNAL_GRID )
  PUBLIC  :: GIGC_Cap_Tropopause_Prs
#endif
!
! !REVISION HISTORY:
!  26 Jun 2010 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: avgpole
!
! !DESCRIPTION: Subroutine AVGPOLE computes average quantity near polar caps,
!  defined by (J = 1, 2) and (J = State\_Grid%NY-1, State\_Grid%NY).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AVGPOLE( State_Grid, Z )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Quantity to be averaged over the pole (usually PS)
    REAL(fp), INTENT(INOUT) :: Z(State_Grid%NX,State_Grid%NY)
!
! !REVISION HISTORY:
!  30 Jan 1998 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J
    REAL(fp)  :: TOTAL_Z1, TOTAL_Z2, TOTAL_Z3, TOTAL_Z4
    REAL(fp)  :: TOTAL_A1, TOTAL_A2, TOTAL_A3, TOTAL_A4

    !=================================================================
    ! AVGPOLE begins here!
    !=================================================================

    IF ( State_Grid%NestedGrid ) THEN
       ! Nested grids typically do not extend to poles. Add better check
       ! here is users wish to define custom nested grids near the poles.
       RETURN
    ENDIF

    TOTAL_Z1 = 0.
    TOTAL_Z2 = 0.
    TOTAL_Z3 = 0.
    TOTAL_Z4 = 0.
    TOTAL_A1 = 0.
    TOTAL_A2 = 0.
    TOTAL_A3 = 0.
    TOTAL_A4 = 0.

    DO I = 1, State_Grid%NX
       TOTAL_Z1 = TOTAL_Z1 + State_Grid%Area_M2(I,1) * Z(I,1)

       TOTAL_Z2 = TOTAL_Z2 + State_Grid%Area_M2(I,2) * Z(I,2)

       TOTAL_Z3 = TOTAL_Z3 + State_Grid%Area_M2(I,State_Grid%NY-1) &
                  * Z(I,State_Grid%NY-1)

       TOTAL_Z4 = TOTAL_Z4 + State_Grid%Area_M2(I,State_Grid%NY) &
                  * Z(I,State_Grid%NY)

       TOTAL_A1 = TOTAL_A1 + State_Grid%Area_M2(I,1)
       TOTAL_A2 = TOTAL_A2 + State_Grid%Area_M2(I,2)
       TOTAL_A3 = TOTAL_A3 + State_Grid%Area_M2(I,State_Grid%NY-1)
       TOTAL_A4 = TOTAL_A4 + State_Grid%Area_M2(I,State_Grid%NY)
    ENDDO

    DO I = 1, State_Grid%NX
       Z(I,1) = (TOTAL_Z1 + TOTAL_Z2) / (TOTAL_A1 + TOTAL_A2)
       Z(I,2) = Z(I,1)
       Z(I,State_Grid%NY-1) = (TOTAL_Z3 + TOTAL_Z4) / (TOTAL_A3 + TOTAL_A4)
       Z(I,State_Grid%NY  ) = Z(I,State_Grid%NY-1)
    ENDDO

  END SUBROUTINE AVGPOLE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: airqnt
!
! !DESCRIPTION: Subroutine AIRQNT sets several members of State\_Met, the
!  meteorology object of derived type MetState, and optionally updates
!  the tracer concentrations to conserve tracer mass when air quantities
!  change.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AIRQNT( Input_Opt, State_Chm, State_Grid, State_Met, &
                     RC, Update_Mixing_Ratio )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE PhysConstants,  ONLY : AIRMW, AVO
    USE Pressure_Mod
    USE Time_Mod,       ONLY : Get_LocalTime
    USE Time_Mod,       ONLY : Get_LocalTime_In_Sec
    USE Time_Mod,       ONLY : Get_Ts_Dyn
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)           :: Input_Opt  ! Input Options object
    TYPE(GrdState), INTENT(IN)           :: State_Grid ! Grid State object
    LOGICAL,        INTENT(IN), OPTIONAL :: Update_Mixing_Ratio ! Default is yes
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT)        :: State_Met ! Meteorology State obj
    TYPE(ChmState), INTENT(INOUT)        :: State_Chm ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)          :: RC        ! Success or failure?
!
! !REMARKS:
!  Met fields updated by AIRQNT:
!  ========================================================================
!  (1)  PEDGE     (REAL(fp)) : Moist air pressure at grid box bottom      [hPa]
!  (2)  PEDGE_DRY (REAL(fp)) : Dry air partial pressure at box bottom     [hPa]
!  (3)  PMID      (REAL(fp)) : Moist air pressure at grid box centroid    [hPa]
!  (4)  PMID_DRY  (REAL(fp)) : Dry air partial pressure at box centroid   [hPa]
!  (5)  PMEAN     (REAL(fp)) : Altitude-weighted mean moist air pressure  [hPa]
!  (6)  PMEAN_DRY (REAL(fp)) : Alt-weighted mean dry air partial pressure [hPa]
!  (7)  DELP      (REAL(fp)) : Delta-P extent of grid box                 [hPa]
!                              (Same for both moist and dry air since we
!                              assume constant water vapor pressure
!                              across box)
!  (8)  AIRDEN    (REAL(fp)) : Mean grid box dry air density            [kg/m^3]
!                              (defined as total dry air mass/box vol)
!  (9)  MAIRDEN   (REAL(fp)) : Mean grid box moist air density          [kg/m^3]
!                              (defined as total moist air mass/box vol)
!  (10) AD        (REAL(fp)) : Total dry air mass in grid box             [kg]
!  (11) ADMOIST   (REAL(fp)) : Total moist air mass in grid box           [kg]
!  (12) BXHEIGHT  (REAL(fp)) : Vertical height of grid box                [m]
!  (13) AIRVOL    (REAL(fp)) : Volume of grid box                         [m^3]
!  (14) MOISTMW   (REAL(fp)) : Molecular weight of moist air in box     [g/mol]
!
! !REVISION HISTORY:
!  30 Jan 1998 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: Dt_Sec
    INTEGER             :: I,         J,          L
    INTEGER             :: L_CG,      L_TP,       N
    REAL(fp)            :: PEdge_Top, Esat,       Ev_mid,    Ev_edge
    REAL(fp)            :: Ev_mean,   PMEAN,      PMEAN_DRY, EsatA
    REAL(fp)            :: EsatB,     EsatC,      EsatD
    REAL(fp)            :: SPHU_kgkg, AVGW_moist, H,         FRAC
    REAL(fp)            :: Pb,        Pt
    LOGICAL             :: UpdtMR

    ! Arrays
    LOGICAL             :: IsLocNoon (State_Grid%NX,State_Grid%NY)
    REAL(f8)            :: LocTime   (State_Grid%NX,State_Grid%NY)
    INTEGER             :: LocTimeSec(State_Grid%NX,State_Grid%NY)

    ! Strings
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
!
! DEFINED PARAMETERS:
!
    ! Conversion factors
    REAL(fp), PARAMETER :: TOFFSET = -2.7315e+2_fp   ! K   -> deg C
    REAL(fp), PARAMETER :: PCONV   = 1.00e+2_fp      ! hPa -> Pa
    REAL(fp), PARAMETER :: RHCONV  = 1.00e-2_fp      ! %   -> fraction
    REAL(fp), PARAMETER :: MCONV   = 1.00e-3_fp      ! g   -> kg

    ! Empirical parameters for water vapor saturation pressure
    ! (Source: Nordquist, 1973. "Numerical Approximiations of
    !  Selected Meteorological Parameters Related to Cloud Physics"
    !  Text quality clarifications from Stipanuk, 1973. "Algorithms
    !  for Generating a Skew-T, Log P Diagram and Computing Selected
    !  Meteorological Quantities")
    REAL(fp), PARAMETER :: ESATP1  =  2.3832241e+1_fp
    REAL(fp), PARAMETER :: ESATP2  = -5.02808e+0_fp
    REAL(fp), PARAMETER :: ESATP3  =  8.1328e-3_fp
    REAL(fp), PARAMETER :: ESATP4  =  3.49149e+0_fp
    REAL(fp), PARAMETER :: ESATP5  = -1.3028844e+3_fp
    REAL(fp), PARAMETER :: ESATP6  = -1.3816e-7_fp
    REAL(fp), PARAMETER :: ESATP7  =  1.1344e+1_fp
    REAL(fp), PARAMETER :: ESATP8  = -3.03998e-2_fp
    REAL(fp), PARAMETER :: ESATP9  = -2.949076e+3_fp

    !=================================================================
    ! AIRQNT begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at AIRQNT (in module GeosCore/dao_mod.F)'
    Dt_Sec   = Get_Ts_Dyn()

    ! Shadow variable for mixing ratio update
    UpdtMR = .TRUE.
    IF ( PRESENT(update_mixing_ratio) ) UpdtMR = update_mixing_ratio

    ! Pre-compute local solar time = UTC + Lon/15
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J   )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Local time only depends on longitude, but longitude is a
       ! function of (I,J) for cubed-sphere grids.  Therefore, use
       ! (I,J) in the call to GET_LOCALTIME.  Obtain the local time
       ! both in hours and in seconds (bmy, 4/18/19)
       LocTime   (I,J) = Get_LocalTime       ( I, J, 1, State_Grid )
       LocTimeSec(I,J) = Get_LocalTime_In_Sec( I, J, 1, State_Grid )

       ! Pick the boxes that are closest to local noon (12hr = 43200 s).
       ! Use local time in seconds, which avoids roundoff issues.
       IsLocNoon(I,J) = ( LocTimeSec(I,J)          <= 43200  .and. &
                          LocTimeSec(I,J) + Dt_Sec >= 43200 )


       ! Land: LWI=1 and ALBEDO less than 69.5%
       State_Met%IsLand(I,J) = ( NINT( State_Met%LWI(I,J) ) == 1 .and. &
                               State_Met%ALBD(I,J)  <  0.695e+0_fp )

       ! Water: LWI=0 and ALBEDO less than 69.5%
       State_Met%IsWater(I,J) = ( NINT( State_Met%LWI(I,J) ) == 0 .and. &
                                State_Met%ALBD(I,J)  <  0.695e+0_fp )

       ! Ice: LWI=2 or ALBEDO > 69.5%
       State_Met%IsIce(I,J) = ( NINT( State_Met%LWI(I,J) ) == 2 .or. &
                              State_Met%ALBD(I,J)  >= 0.695e+0_fp )

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !=============================================================
    ! Update air quantities
    !=============================================================

    !$OMP PARALLEL DO                                       &
    !$OMP DEFAULT( SHARED                                 ) &
    !$OMP PRIVATE( I,       J,         L,       Pedge_Top ) &
    !$OMP PRIVATE( EsatA,   EsatB,     EsatC,   EsatD     ) &
    !$OMP PRIVATE( Esat,    SPHU_kgkg, Ev_mid,  Ev_edge   ) &
    !$OMP PRIVATE( PMean,   Ev_mean,   PMean_Dry          )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !=============================================================
       ! Set wet air pressures [hPa]
       ! (lower edge, delta, centroid, and spatially-weighted mean)
       !=============================================================

       ! Pressure at bottom edge of grid box [hPa]
       State_Met%PEDGE(I,J,L) = GET_PEDGE(I,J,L)

       ! Pressure at top edge of grid box [hPa]
       PEdge_Top = GET_PEDGE(I,J,L+1)

       ! Pressure at bottom edge of grid box [hPa](level State_Grid%NZ+1 only)
       IF ( L == State_Grid%NZ ) THEN
          State_Met%PEDGE(I,J,L+1) = PEdge_Top
       ENDIF

       ! Pressure difference between top & bottom edges [hPa]
       State_Met%DELP(I,J,L) = State_Met%PEDGE(I,J,L) - PEdge_Top

       ! Arithmetic average pressure in grid box [hPa] defined as
       ! ( PEDGE(L) + PEDGE(L+1) ) / 2. This represents the grid box
       ! mass centroid pressure. Use in the ideal gas law yields
       ! a local air density at the centroid.
       State_Met%PMID(I,J,L) = GET_PCENTER( I, J, L )

       !=============================================================
       ! Calculate water vapor saturation pressure [hPa]
       ! from temperature
       !=============================================================

       ! Prepare the intermediate terms
       EsatA = ESATP1 + ESATP2 * log10( State_Met%T(I,J,L) )
       EsatB = ESATP3 * 10**( ESATP4 + ESATP5 / State_Met%T(I,J,L) )
       EsatC = ESATP6 * 10**( ESATP7 + ESATP8 * State_Met%T(I,J,L) )
       EsatD = ESATP9 / State_Met%T(I,J,L)

       ! Saturation water vapor pressure [hPa]
       Esat  = 10**( EsatA + EsatB + EsatC + EsatD )

       !=============================================================
       ! Set AVGW, the water vapor volume mixing ratio from
       ! specific humidity [mol H2O / mol dry air]
       !=============================================================
       !
       ! vol H2O       dry air    H2O         kg H2O      kg wet air
       ! ----------  = molec wt / molec wt * ---------- * ---------
       ! vol dry air   [g/mol]    [g/mol]    kg wet air   kg dry air
       !
       !      thus AVGW = AIRMW * SPHU_kgkg / ( H2OMW * (1-SPHU_kgkg) )
       !
       ! where,
       !        SPHU_kgkg = specific humidity in [kg/kg]
       !

       ! Convert specific humidity from [g/kg] to [kg/kg]
       SPHU_kgkg = State_Met%SPHU(I,J,L) * MCONV

       ! Water vapor volume mixing ratio [v/v dry air]
       State_Met%AVGW(I,J,L) = AIRMW * SPHU_kgkg / &
                               ( H2OMW * (1.0e+0_fp - SPHU_kgkg ) )

       !=============================================================
       ! Calculate water vapor partial pressures [hPa] from relative
       ! humidity
       !=============================================================

       ! At vertical midpoint of grid box
       Ev_mid  = State_Met%RH(I,J,L) * RHCONV * Esat

       ! At bottom edge of grid box
       Ev_edge = State_Met%PEDGE(I,J,L) * Ev_mid / State_Met%PMID(I,J,L)

       !=============================================================
       ! Set grid box height [m]
       !=============================================================
       !
       ! BXHEIGHT is the height (Delta-Z) of grid box (I,J,L)
       ! in meters. It is calculated using the hypsometric equation.
       ! A virtual temperature is calculated to enable use of
       ! of moist air pressure with dry air molecular weight.
       !
       !              Gas constant   virtual
       !              for dry air  * grid box
       !              [J/K/mol]      temp [K]       bottom edge P
       ! height [m] = ----------------------- * ln(---------------)
       !                     g [m/s2]                top edge P
       !
       !  where,
       !
       !                         Grid box temperature [K]
       !    Virtual  = -------------------------------------------
       !    Temp [K]       H20 partial pressure          MW_H2O
       !               1 - -------------------- * ( 1 - --------- )
       !                    moist air  pressure         MW_dryair
       !
       ! Source: Wallace and Hobbes "Atmospheric Science: An
       !         Introductory Survey"
       !
       ! Assume constant temperature and moisture across grid box.
       !

       ! Grid box potential temperature [K]
       ! NOTE: Due to the parallelization, we cannot assume that
       ! State_Met%PEDGE(I,J,1) has been defined.  So always call
       ! GET_PEDGE(I,J,1) to return the proper surface pressure
       ! (bmy, 2/23/18)
       State_Met%THETA(I,J,L)  = State_Met%T(I,J,L) * &
                                 ( GET_PEDGE( I, J, 1 ) / &
                                 State_Met%PMID(I,J,L) )**0.286

       ! Grid box virtual temperature [K]
       State_Met%TV(I,J,L) = State_Met%T(I,J,L) / (1 - Ev_edge / &
                             State_Met%PEDGE(I,J,L) * ( 1 - H2OMW / AIRMW ) )

       ! Grid box box height [m]
       State_Met%BXHEIGHT(I,J,L) = Rdg0 * State_Met%TV(I,J,L) * &
                                   LOG( State_Met%PEDGE(I,J,L) / PEdge_Top )

       !==============================================================
       ! Set grid box volume [m3]
       !==============================================================
       !
       ! AIRVOL is the volume of grid box (I,J,L) in meters^3
       !
       ! AREA_M2 is the Delta-X * Delta-Y surface area of grid
       ! boxes (I,J,L=1), that is, at the earth's surface.
       !
       ! Since the thickness of the atmosphere is much smaller
       ! than the radius of the earth, we can make the "thin
       ! atmosphere" approximation, namely:
       !
       !               (Rearth + h) ~ Rearth
       !
       ! Therefore, the Delta-X * Delta-Y surface area of grid
       ! boxes that are above the earth's surface will be
       ! approx. the same as AREA_M2.  Thus we are justified
       ! in using AREA_M2 for grid boxes (I, J, L > 1)
       !
       State_Met%AIRVOL(I,J,L) = State_Met%BXHEIGHT(I,J,L) * &
                                 State_Grid%AREA_M2(I,J)

       !==============================================================
       ! Set grid box dry air partial pressures at grid box edges
       ! and at grid box centroid (arithmetic mean pressure level)
       ! [hPa]. Assume constant humidity across grid box.
       !==============================================================

       ! Partial pressure of dry air at lower edge of grid box [hPa]
       State_Met%PEDGE_DRY(I,J,L) = State_Met%PEDGE(I,J,L) - Ev_edge

       ! Set dry air partial pressure for level State_Grid%NZ+1 lower edge
       IF ( L == State_Grid%NZ ) THEN
          State_Met%PEDGE_DRY(I,J,L+1) = Pedge_Top - Ev_edge
       ENDIF

       ! Partial pressure of dry air at box centroid [hPa]
       State_Met%PMID_DRY(I,J,L) = State_Met%PMID(I,J,L) - Ev_mid

       ! Set previous dry P difference to current dry P difference
       ! prior to updating with new met values
       State_Met%DP_DRY_PREV(I,J,L) = State_Met%DELP_DRY(I,J,L)

       ! Update dry pressure difference as calculated from the
       ! dry surface pressure
       State_Met%DELP_DRY(I,J,L) = GET_DELP_DRY(I,J,L)

       !==============================================================
       ! Set mass of dry air in grid box [kg]
       !==============================================================

       ! AD = mass of dry air in grid box (I,J,L)
       !
       ! given by:
       !
       !        Dry air pressure   Grid box   100 [Pa]   1 [s2]
       ! Mass = difference       * surface  * -------- * -----
       !        across grid        area       1 [hPa]    g [m]
       !        box [hPa]          [m2]
       !
       ! NOTES:
       ! Dry air pressure difference across grid box is calculate
       ! from the surface dry pressure and GMAO A and B parameters
       ! (see GeosUtil/pressure_mod.F). The dry surface pressure
       ! that is used is calculated from the GMAO wet surface pressure
       ! (or time-interpolated value) and the A and B parameters
       ! (see SET_DRY_SURFACE_PRESSURE). It is not derived from the
       ! wet edge pressures. This distinction is important for
       ! mass conservation.
       State_Met%AD(I,J,L) = ( State_Met%DELP_DRY(I,J,L) * G0_100 ) * &
                             State_Grid%AREA_M2(I,J)

       !==============================================================
       ! Set grid box dry air partial pressures at grid box
       ! altitude-weighted mean pressure [hPa]
       ! Assume constant humidity across grid box.
       !==============================================================

       ! Mean altitude-weighted pressures in grid box [hPa] defined as
       ! average P(z) over z. Use in the ideal gas law yields a
       ! mean density equivalent to total mass per volume in grid box.
       PMEAN = State_Met%DELP(I,J,L) / log( State_Met%PEDGE(I,J,L) / PEdge_Top )
       Ev_mean   = PMEAN * Ev_mid / State_Met%PMID(I,J,L)
       PMEAN_DRY = PMEAN - Ev_mean

       ! NOTE: Try the below definition in the future to change the
       ! AIRDEN equation to use the ideal gas law, thereby removing
       ! area-independence of AIRDEN
       !State_Met%PMEAN_DRY( I,J,L ) = State_Met%DELP_DRY(I,J,L) / &
       !                               log( State_Met%PEDGE(I,J,L) / &
       !                                    PEdge_Top )

       !==============================================================
       ! Set grid box densities
       !==============================================================

       ! Set grid box dry air density [kg/m3]
       !
       ! NOTE: Air density is derived from dry surface pressure
       ! following implementation of the moisture fix, and therefore
       ! its equation must be dry air mass divided by volume.
       ! This is because the level pressures derived from the dry
       ! surface pressure may be used for extracting mass per level,
       ! but not a representative pressure for that level. The ideal
       ! gas law requires a representative level pressure. Eventually
       ! a representative pressure must be derived that, when used in
       ! the ideal gas law, will replicate the air density defined as
       ! mass/volume. Moving to use of the ideal gas law is necessary
       ! for grid-independence (ewl, 9/16/16)
       State_Met%AIRDEN(I,J,L) = State_Met%AD(I,J,L) / State_Met%AIRVOL(I,J,L)

       ! Set grid box dry air number density [molec/cm3]
       State_Met%AIRNUMDEN(I,J,L) = State_Met%AIRDEN(I,J,L) * 1e-3_fp * &
                                    AVO / AIRMW

       ! Set grid box moist air density [kg/m3] using the ideal gas law
       ! and pressures derived from GMAO pressure
       !
       !  MAIRDEN = density of moist air [kg/m^3],
       !  given by:
       !
       !            Partial      Molec     Partial       Molec
       !            pressure   * wt of   + pressure of * wt of
       !            of dry air   air       water vapor   water
       ! Moist      [hPa]        [g/mol]   [hPa]         [g/mol]
       ! Air     =  ------------------------------------------------
       ! Density    Universal gas constant * Temp * 1000  * 0.01
       !                   [J/K/mol]         [K]   [g/kg]   [hPa/Pa]
       !
       ! NOTE: MAIRDEN is used in wetscav_mod only
       State_Met%MAIRDEN(I,J,L) = ( PMEAN_DRY * AIRMW + Ev_mean * H2OMW ) * &
                                  PCONV * MCONV / RSTARG / State_Met%T(I,J,L)

       !==============================================================
       ! Define the various query fields of State_Met
       !
       ! NOTE: For convenience, we set State_Met%InPbl in routine
       ! COMPUTE_PBL_HEIGHT (in module GeosCore/pbl_mix_mod.F).
       !
       ! NOTE: For certain queries we test against level numbers,
       ! (e.g. LLSTRAT, LLCHEM), but should really test level
       ! pressure edges, so that this algorithm will be robust if
       ! we switch to different met fields or interface with
       ! different ESMs.  Add this at a later time. (bmy, 1/8/18)
       !==============================================================

       ! Is this grid box within the troposphere?
       State_Met%InTroposphere(I,J,L) = &
            ( State_Met%PEDGE(I,J,L) > State_Met%TROPP(I,J) )

       ! Is this grid box within the stratosphere or mesosphere?
       State_Met%InStratMeso(I,J,L) = &
            ( .not. State_Met%InTroposphere(I,J,L) )

       ! Is this grid box within the stratosphere (but not mesosphere)?
       State_Met%InStratosphere(I,J,L) = &
            ( L <= State_Grid%MaxStratLev .and. State_Met%InStratMeso(I,J,L) )

       ! Is grid box (I,J,L) within the chemistry grid?
       IF ( L > State_Grid%MaxChemLev ) THEN

          ! Chemistry is not done higher than the mesopause
          State_Met%InChemGrid(I,J,L) = .FALSE.

       ELSE

          IF ( Input_Opt%LUCX ) THEN

             ! UCX mechanisms: Chemistry grid goes up to stratopause
             State_Met%InChemGrid(I,J,L) = ( L <= State_Grid%MaxChemLev )

          ELSE

             ! "tropchem" mechanisms: Chemistry grid goes up to tropopause
             State_Met%InChemGrid(I,J,L) = State_Met%InTroposphere(I,J,L)

          ENDIF

       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !=================================================================
    ! Compute more tropopause and chemistry grid quantities.  This
    ! has to be done after the State_Met%InChemGrid etc. fields have
    ! been totally initialized.
    !
    ! Also compute if it is near local noon in a grid box.
    ! This will be used for the J-value diagnostics.
    !=================================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L_CG, L_TP, H, Pb, Pt, FRAC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !--------------------------------------------------------------
       ! Local solar time and box that is nearest to local noon
       !--------------------------------------------------------------
       State_Met%LocalSolarTime(I,J) = LocTime(I,J)
       State_Met%IsLocalNoon(I,J)    = IsLocNoon(I,J)

       !--------------------------------------------------------------
       ! Compute highest level of chemistry grid and the troposphere
       !
       ! NOTE: The COUNT intrinsic function returns the number of
       ! locations where the column State_Met%InChemGrid(I,J,:) is
       ! TRUE.  This is also equivalent to the maximum level of the
       ! chemistry grid.  Ditto for State_Met%InTroposphere.
       !--------------------------------------------------------------

       ! Highest level in the column at (I,J)
       ! that is fully within the chemistry grid
       L_CG  = COUNT( State_Met%InChemGrid(I,J,:) )

       ! Highest level in the column at (I,J)
       ! that is fully within the tropopause
       L_TP  = COUNT( State_Met%InTroposphere(I,J,:) )

       !--------------------------------------------------------------
       ! Compute tropopause height
       !--------------------------------------------------------------

       ! Get height (from surface to top edge) of all boxes that lie
       ! totally w/in the troposphere.  NOTE: Grid box (I,J,L_TP-1)
       ! is the highest purely tropospheric grid box in the column.
       !
       ! BMY COMMENT: L_TP may be the highest purely tropospheric
       ! grid box in the column, not L_TP-1. That might have been
       ! due to historical baggage.  In any case, this is only used
       ! to compute State_Met%TropHt which is a purely diagnostic
       ! output, so it won't make a big difference even if it is
       ! incorrect.  Check into this later on. (bmy, 1/17/18)
       H = SUM( State_Met%BXHEIGHT(I,J,1:L_TP-1) )

       ! Get the pressures [hPa] at the bottom and top edges
       ! of the grid box in which the tropopause occurs
       Pb = State_Met%PEDGE(I,J,L_TP  )
       Pt = State_Met%PEDGE(I,J,L_TP+1)

       ! FRAC is the fraction of the grid box (I,J,L_TP)
       ! that lies totally within the troposphere
       FRAC = ( Pb - State_Met%TROPP(I,J) ) / ( Pb - Pt )

       ! Add to H the height [m] of the purely tropospheric
       ! fraction of grid box (I,J,L_TP)
       H = H + ( FRAC * State_Met%BXHEIGHT(I,J,L_TP) )

       !--------------------------------------------------------------
       ! Save in the relevant fields of State_Met
       !--------------------------------------------------------------
       State_Met%ChemGridLev(I,J) = L_CG           ! Max chemgrid level
       State_Met%TropLev    (I,J) = L_TP           ! Max tropopause level
       State_Met%TropHt     (I,J) = H / 1000.0_fp  ! Troposphere ht [km]

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !=================================================================
    ! Update species concentrations with updated air quantities
    ! if update mixing ratio flag indicates to do so
    !=================================================================

    ! Update the mixing ratio with new air quantities such that species
    ! mass is conserved if (1) update_mixing_ratio flag is not passed,
    ! or (2) update_mixing_ratio flag is passed as true.
    ! NOTE: The only places where mixing ratio is not currently updated
    ! following air quantity change is during GEOS-Chem initialization and
    ! in transport after the pressure fixer is applied
    IF ( UpdtMR ) THEN
       !IF ( .not. PRESENT( update_mixing_ratio ) .or. update_mixing_ratio ) THEN

       ! The concentration update formula works only for dry mixing ratios
       ! (kg/kg or v/v) so check if units are correct
       IF ( State_Chm%Spc_units == 'kg/kg dry' .or. &
            State_Chm%Spc_units == 'v/v dry' ) THEN

          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, L, N )
          DO N = 1, State_Chm%nSpecies
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N) * &
                                          State_Met%DP_DRY_PREV(I,J,L) / &
                                          State_Met%DELP_DRY(I,J,L)
          ENDDO
          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ELSE
          ErrMsg = 'Incorrect species units: ' // TRIM( State_Chm%Spc_Units )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE AIRQNT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: interp
!
! !DESCRIPTION: Subroutine INTERP linearly interpolates GEOS-Chem I6 and I3
!  fields (winds, surface pressure, temperature, surface albedo, specific
!  humidity etc.) to the current dynamic timestep.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INTERP( NTIME0, NTIME1, NTDT, Input_Opt, State_Grid, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: NTIME0     ! Elapsed time [s] at
                                                !  start of outer time step
    INTEGER,        INTENT(IN)    :: NTIME1     ! Elapsed time [s] at
                                                !  current time
    INTEGER,        INTENT(IN)    :: NTDT       ! Dynamic timestep [s]
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid ! State Grid object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !REMARKS:
!  Different met fields are archived at I6 (instantaneous 6-hr) and
!  I3 (instantaneous 3-hr) time resolution depending on the specific product.
!
! !REVISION HISTORY:
!  30 Jan 1998 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER   :: I,        J,        L
    REAL(fp)  :: D_NTIME0, D_NTIME1, D_NDT
    REAL(fp)  :: D_NTDT,   TM,       TC2
    REAL(fp)  :: YSOUTH,   YNORTH

    !=================================================================
    ! Initialization
    !=================================================================

    ! Convert time variables from FLOAT to DBLE
    D_NTIME0 = DBLE( NTIME0 )
    D_NTIME1 = DBLE( NTIME1 )
    D_NTDT   = DBLE( NTDT   )

    D_NDT    = 10800e+0_fp ! For 3-hr instantaneous fields

    ! Fraction of timestep elapsed at mid point of this dyn timestep
    TM  = ( D_NTIME1 + D_NTDT/2e+0_fp - D_NTIME0 ) / D_NDT

    ! Fraction of timestep elapsed at the end of this dyn timestep
    TC2 = ( D_NTIME1 + D_NTDT - D_NTIME0 ) / D_NDT

    ! For I3 fields, need to reset fraction after 3 hours (ewl, 5/12/2015)
    IF ( TM > 1.0e+0_fp ) THEN
       TM  = TM  - 1.0e+0_fp
       TC2 = TC2 - 1.0e+0_fp
    ENDIF

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, YSOUTH, YNORTH )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !----------------------------------------------------
       ! Interpolate 2D variables
       !----------------------------------------------------
       IF ( L == 1 ) THEN

          ! North & south edges of box
          YSOUTH = State_Grid%YEdge(I,J  )
          YNORTH = State_Grid%YEdge(I,J+1)

          ! Interpolate wet pressure [hPa] to the end of the dynamic timestep
          State_Met%PSC2_WET(I,J)  = State_Met%PS1_WET(I,J) + &
                                   ( State_Met%PS2_WET(I,J) - &
                                     State_Met%PS1_WET(I,J) ) * TC2

          ! Do the same for dry pressure [hPa] (ewl, 5/4/16)
          State_Met%PSC2_DRY(I,J)  = State_Met%PS1_DRY(I,J) + &
                                   ( State_Met%PS2_DRY(I,J) - &
                                     State_Met%PS1_DRY(I,J) ) * TC2

          ! Even though TROPP is a 3-hour average field, we
          ! we still need to make sure to cap TROPP in the
          ! polar regions (if the entire box is outside 60S-60N)
          ! so that we don't do chemistry at an abnormally high
          ! altitude.  Set TROPP in the polar regions to 200 hPa.
          ! (jal, phs, bmy, 9/18/07)
          IF ( YSOUTH >= 60.0_fp .OR. YNORTH <= -60.0_fp ) THEN
             State_Met%TROPP(I,J) = MAX( State_Met%TROPP(I,J), 200e+0_fp )
          ENDIF
       ENDIF

       !----------------------------------------------------
       ! Interpolate 3D variables
       !----------------------------------------------------

       ! Interpolate temperature
       State_Met%T(I,J,L) = State_Met%TMPU1(I,J,L) + &
            ( State_Met%TMPU2(I,J,L) - State_Met%TMPU1(I,J,L) ) * TM

       ! Interpolate specific humidity
       State_Met%SPHU(I,J,L) = State_Met%SPHU1(I,J,L) + &
            ( State_Met%SPHU2(I,J,L) - State_Met%SPHU1(I,J,L) ) * TM

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE INTERP
!EOC
#if defined( ESMF_ ) || defined( EXTERNAL_GRID )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_cap_tropopause_prs
!
! !DESCRIPTION: Subroutine GIGC\_CAP\_TROPOPAUSE\_PRS caps the tropopause
!  pressure in polar latitudes to 200 hPa, so that we don't end up doing
!  troposheric chemistry too high over the poles.  This is done in the
!  standalone GEOS-Chem, and we also need to apply this when running
!  GEOS-Chem within the GEOS-5 GCM.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Cap_Tropopause_Prs( Input_Opt, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt     ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid    ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met     ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure
!
! !REMARKS:
!  Jennifer Logan (see correspondence below) suggested that we should cap the
!  variable tropopause at 200hPa in near-polar regions (90-60S and 60-90N),
!  to avoid the problem with anomalously high tropopause heights at high
!  latitudes. This fix was standardized in GEOS-Chem v7-04-13.
!                                                                             .
!  Jennifer Logan wrote:
!     I think we should restrict the tropopause at latitudes > 60 deg. to
!     pressures greater than 200 mb (about 11 km). From Fig. 3 in Seidel and
!     Randel, there are tropopause (TP) heights as high as 13.5 km in the
!     Antarctic (median height is ~9.8 km, 250 mb), but I don't think we want
!     to be doing trop. chem there. The median TP pressure at ~80 N is ~300 mb,
!     compared to ~250 mb at 70-85 S. The extratropical TP heights are higher
!     (lower pressure) in the SH than in the NH according to Fig. 3.
!     This approach is also very easy to explain in a paper.
!
! !REVISION HISTORY:
!  14 Mar 2013 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I,       J
    REAL*8  :: YSOUTH,  YNORTH

    ! Assume success
    RC = GC_SUCCESS

      ! Return if option not set
#ifdef MODEL_GEOS
    IF ( .NOT. Input_Opt%LCAPTROP ) RETURN
#endif

    ! Loop over grid boxes on this PET
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! North & south edges of box
       YSOUTH = State_Grid%YEdge(I,J  )
       YNORTH = State_Grid%YEdge(I,J+1)

       ! Cap tropopause height at 200 hPa polewards of 60N and 60S
       IF ( YSOUTH >= 60d0 .or. YNORTH <= -60d0 ) THEN
          State_Met%TROPP(I,J) = MAX( State_Met%TROPP(I,J), 200d0 )
       ENDIF

    ENDDO
    ENDDO

  END SUBROUTINE GIGC_Cap_Tropopause_Prs
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_dry_surface_pressure
!
! !DESCRIPTION: Subroutine SET\_DRY\_SURFACE\_PRESSURE sets the dry
!  surface pressures PS1\_DRY or PS2\_DRY by removing the water vapor from
!  the column constructed with MET pressure fields PS1\_WET or PS2\_WET.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, PS_ID)
!
! !USES:
!
    USE ERROR_MOD,            ONLY : CHECK_VALUE
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE PRESSURE_MOD,         ONLY : GET_AP, GET_BP
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)     :: PS_ID       ! 1 = PS1, 2 = PS2
    TYPE(GrdState), INTENT(IN)     :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT)  :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  This subroutine is an adaptation of the GEOS-Chem moisture fix
!  implemented by Meemong Lee (JPL) in the adjoint model. Like PS1_WET
!  and PS2_WET, from which PS1_DRY and PS2_DRY are derived, these values
!  are interpolated within routine INTERP to derive instantaneous PSC2_DRY.
!  Note that PSC2_WET and PSC2_DRY are not used to fetch pressures and
!  calculate air quantities until after advection. Also note that the
!  dry surface pressure calculated in this routine may be used to
!  calculate delta dry pressures across a level by utilitizing GMAO
!  parameters Ap and Bp but should not be used to extract dry pressure
!  edge values as a height proxy.
!
! !REVISION HISTORY:
!  03 May 2002 - E. Lundgren - Initial version (M. Lee, JPL)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L
    REAL(fp)           :: PEDGE_BOT, PEDGE_TOP

    ! Pointers
    REAL(fp), POINTER  :: PS_WET(:,:) => NULL()
    REAL(fp), POINTER  :: PS_DRY(:,:) => NULL()
    REAL(fp), POINTER  :: SPHU(:,:,:) => NULL()

    !=================================================================
    ! SET_DRY_SURFACE_PRESSURE begins here!
    !=================================================================

    ! Set pointer to the appropriate humidity
    IF ( PS_ID == 1 ) THEN
       SPHU   => State_Met%SPHU1
    ELSE IF ( PS_ID == 2 ) THEN
       SPHU   => State_Met%SPHU2
    ENDIF

    ! Set pointers to State_Met pressure fields.
    IF ( PS_ID == 1 ) THEN
       PS_WET => State_Met%PS1_WET
       PS_DRY => State_Met%PS1_DRY
    ELSE IF ( PS_ID == 2 ) THEN
       PS_WET => State_Met%PS2_WET
       PS_DRY => State_Met%PS2_DRY
    ELSE
       ! Throw an error (to be added)
    ENDIF

    ! Reset dry surface pressure to TOA value
    PS_DRY = GET_AP(State_Grid%NZ+1)

    ! Calculate dry surface pressure from GMAO wet pressure as the
    ! column sum of wet delta pressures with humidity removed
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, PEDGE_BOT, PEDGE_TOP )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
    DO L = 1, State_Grid%NZ
       PEDGE_BOT   = GET_AP(L) + GET_BP(L) * PS_WET(I,J)
       PEDGE_TOP   = GET_AP(L+1) + GET_BP(L+1) * PS_WET(I,J)
       PS_DRY(I,J) = PS_DRY(I,J) + ( PEDGE_BOT - PEDGE_TOP ) * &
                     ( 1.e+0_fp - SPHU(I,J,L) * 1.0e-3_fp )
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! If dry pressure is negative, set equal to moist pressure
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       IF ( PS_DRY(I,J) < 0.e+0_fp) THEN
          PS_DRY(I,J) = PS_WET(I,J)
       ENDIF
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Nullify pointers
    PS_WET  => NULL()
    PS_DRY  => NULL()
    SPHU    => NULL()

  END SUBROUTINE SET_DRY_SURFACE_PRESSURE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_obk
!
! !DESCRIPTION: Function GET\_OBK returns the Monin-Obhukov length at a grid
!  box (I,J).
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_OBK( I, J, State_Met ) RESULT( OBK )
!
! !USES:
!
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN) :: I           ! Longitude index
    INTEGER,         INTENT(IN) :: J           ! Latitude  index
    TYPE(MetState),  INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
    REAL(fp)                    :: OBK         ! Monin-Obhukhov length
!
! !REMARKS:
!
!
! !REVISION HISTORY:
!  25 May 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! For all GEOS met fields:
    !
    ! The direct computation of the Monin-Obhukov length is:
    !
    !            - Air density * Cp * T(surface air) * Ustar^3
    !    OBK =  -----------------------------------------------
    !              Kappa       * g  * Sensible Heat flux
    !
    ! Cp    = 1000 J / kg / K = specific heat of air at constant P
    ! Kappa = 0.4             = Von Karman's constant
    !
    !
    !  Also test the denominator in order to prevent div by zero.
    !=================================================================

    ! Local variables
    REAL(fp)            :: NUM, DEN

    ! Parameters
    REAL(fp), PARAMETER :: KAPPA = 0.4e+0_fp
    REAL(fp), PARAMETER :: CP    = 1000.0e+0_fp

    ! Numerator
    NUM = - State_Met%AIRDEN(I,J,1) * CP                   * &
            State_Met%TS(I,J)       * State_Met%USTAR(I,J) * &
            State_Met%USTAR(I,J)    * State_Met%USTAR(I,J)

    ! Denominator
    DEN =  KAPPA * g0 * State_Met%HFLUX(I,J)

    ! Prevent div by zero
    IF ( ABS( DEN ) > 0e+0_fp ) THEN
       OBK = NUM / DEN
    ELSE
       OBK = 1.0e+5_fp
    ENDIF

  END FUNCTION GET_OBK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_cosine_sza
!
! !DESCRIPTION: Routine GET\_COSINE\_SZA is a driver for calling the COSSZA
!  routine from dao\_mod.F.  This routine calls COSSZA twice.  The first call
!  computes the sun angles at the current time and midpoint of the current
!  chemistry time step.  The second call computes the sun angles 5 hours
!  prior to the current time (for the PARANOX ship emissions plume model).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_COSINE_SZA( Input_Opt, State_Grid, State_Met, RC  )
!
! USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE JULDAY_MOD,         ONLY : JULDAY
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_DAY_OF_YEAR
    USE TIME_MOD,           ONLY : GET_DAY
    USE TIME_MOD,           ONLY : GET_GMT
    USE TIME_MOD,           ONLY : GET_HOUR
    USE TIME_MOD,           ONLY : GET_MINUTE
    USE TIME_MOD,           ONLY : GET_SECOND
    USE TIME_MOD,           ONLY : GET_MONTH
    USE TIME_MOD,           ONLY : GET_YEAR
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  07 Feb 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete histor
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARAIBLES
!
    ! Scalars
    INTEGER :: DOY, YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
    REAL*8  :: DDAY, JD,  GMT

    !=================================================================
    ! Get cosine(SZA) at the current time (SUNCOS) and at the
    ! midpoint of the chemistry timestep (SUNCOS_MID)
    !=================================================================

    ! Assume success
    RC     = GC_SUCCESS

    ! Current time
    DOY    = GET_DAY_OF_YEAR()                       ! Current day of year
    YEAR   = GET_YEAR()                              ! Current year
    MONTH  = GET_MONTH()                             ! Current month
    DAY    = GET_DAY()                               ! Current day of month
    HOUR   = GET_HOUR()                              ! Current GMT hour
    MINUTE = GET_MINUTE()                            ! Current GMT minutes
    SECOND = GET_SECOND()                            ! Current GMT seconds
    GMT    = GET_GMT()                               ! Current GMT
    DDAY   = DAY + ( HOUR/24d0 ) + &                 ! Current decimal day
             ( MINUTE/1440d0 ) + ( SECOND/3600d0 )
    JD     = JULDAY( YEAR, MONTH, DDAY )             ! Current Julian date

    ! Compute cosine(SZA) quantities for the current time
    CALL COSSZA( DOY, HOUR, State_Grid, State_Met )

  END SUBROUTINE GET_COSINE_SZA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cossza
!
! !DESCRIPTION: COSSZA computes the cosine of the solar zenith angle, given
!  the day of the year and GMT hour.  The cosine of the solar zenith
!  angle is returned at both the current time and at the midpoint of the
!  chemistry timestep (i.e. for the centralized chemistry timestep option).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE COSSZA( DOY, GMT_HOUR, State_Grid, State_Met )
!
! !USES:
!
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MINUTE, GET_SECOND
    USE TIME_MOD,           ONLY : GET_LOCALTIME
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: DOY          ! Day of the year
    INTEGER,        INTENT(IN)    :: GMT_HOUR     ! Hour of day
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met    ! Meteorology State
!
! !REMARKS:
!  Hour angle (AHR) is a function of longitude.  AHR is zero at solar noon,
!  and increases by 15 deg for every hour before or after solar noon.  Hour
!  angle can be thought of as the time in hours since the sun last passed
!  the meridian (i.e. the time since the last local noon).
!                                                                             .
!  The cosine of the solar zenith angle (SZA) is given by:
!                                                                             .
!     cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
!                                                                             .
!  where LAT = the latitude angle,
!        DEC = the solar declination angle,
!        AHR = the hour angle, all in radians.
!                                                                             .
!  If SUNCOS < 0, then the sun is below the horizon, and therefore does not
!  contribute to any solar heating.
!
! !REVISION HISTORY:
!  21 Jan 1998 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: I,        J
    INTEGER              :: SECOND,   MINUTE,   TS_SUN,   FACTOR
    REAL(fp)             :: TIMLOC,   DEC,      C_DEC
    REAL(fp)             :: S_DEC,    R,        YMID_R
    REAL(fp)             :: C_YMID_R, S_YMID_R
    REAL(fp)             :: AHR,      SUNCOS,   SUNCOS_MID
    REAL(f8)             :: GMT_CUR,  GMT_MID
!
! !DEFINED PARAMETERS:
!
    ! Coefficients for solar declination angle
    REAL(fp),  PARAMETER :: A0 = 0.006918e+0_fp
    REAL(fp),  PARAMETER :: A1 = 0.399912e+0_fp
    REAL(fp),  PARAMETER :: A2 = 0.006758e+0_fp
    REAL(fp),  PARAMETER :: A3 = 0.002697e+0_fp
    REAL(fp),  PARAMETER :: B1 = 0.070257e+0_fp
    REAL(fp),  PARAMETER :: B2 = 0.000907e+0_fp
    REAL(fp),  PARAMETER :: B3 = 0.000148e+0_fp

    !=================================================================
    ! Initialization
    !=================================================================

    ! Quantities for central chemistry timestep
    TS_SUN   = GET_TS_CHEM()                         ! Chemistry interval
    SECOND   = GET_SECOND()                          ! Current seconds
    MINUTE   = GET_MINUTE()                          ! Current minutes
    FACTOR   = ( MINUTE * 60 + SECOND ) / TS_SUN     ! Multiplying factor

    ! GMT at the current time
    GMT_CUR  = DBLE( GMT_HOUR          ) &
             + ( DBLE( TS_SUN * FACTOR ) / 3600e+0_fp )

    ! GMT at the midpoint of the chemistry time interval
    GMT_MID  = ( DBLE( GMT_HOUR        )        ) &
             + ( DBLE( TS_SUN * FACTOR ) / 3600e+0_fp ) &
             + ( DBLE( TS_SUN / 2      ) / 3600e+0_fp )

    ! Path length of earth's orbit traversed since Jan 1 [radians]
    R        = ( 2e+0_fp * PI / 365e+0_fp ) * DBLE( DOY - 1 )

    ! Solar declination angle (low precision formula) [radians]
    DEC      = A0 - A1*COS(         R ) + B1*SIN(         R ) &
                  - A2*COS( 2e+0_fp*R ) + B2*SIN( 2e+0_fp*R ) &
                  - A3*COS( 3e+0_fp*R ) + B3*SIN( 3e+0_fp*R )

    ! Pre-compute sin & cos of DEC outside of DO loops (for efficiency)
    S_DEC    = SIN( DEC )
    C_DEC    = COS( DEC )

    !=================================================================
    ! Compute cosine of solar zenith angle
    !=================================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,      J,   YMID_R, S_YMID_R,  C_YMID_R ) &
    !$OMP PRIVATE( TIMLOC, AHR, SUNCOS, SUNCOS_MID          )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Latitude of grid box [radians]
       YMID_R     = State_Grid%YMid_R(I,J)

       ! Pre-compute sin & cos of DEC outside of I loop (for efficiency)
       S_YMID_R   = SIN( YMID_R )
       C_YMID_R   = COS( YMID_R )

       !==============================================================
       ! Compute cosine of SZA at the current GMT time
       !==============================================================

       ! Local time at box (I,J) [hours]
       TIMLOC     = GET_LOCALTIME( I, J, 1, State_Grid, GMT=GMT_CUR )

       ! Hour angle at box (I,J) [radians]
       AHR        = ABS( TIMLOC - 12e+0_fp ) * 15e+0_fp * PI_180

       ! Cosine of solar zenith angle at box (I,J) [unitless]
       SUNCOS     = ( S_YMID_R * S_DEC              ) &
                  + ( C_YMID_R * C_DEC * COS( AHR ) )

       !==============================================================
       ! Compute cosine of SZA at the midpoint of the chem timestep
       ! Required for photolysis, chemistry, emissions, drydep
       !==============================================================

       ! Local time [hours] at box (I,J) at the midpt of the chem timestep
       TIMLOC     = GET_LOCALTIME( I, J, 1, State_Grid, GMT=GMT_MID )

       ! Hour angle at box (I,J) [radians]
       AHR        = ABS( TIMLOC - 12e+0_fp ) * 15e+0_fp * PI_180

       ! Corresponding cosine( SZA ) at box (I,J) [unitless]
       SUNCOS_MID = ( S_YMID_R * S_DEC              ) &
                  + ( C_YMID_R * C_DEC * COS( AHR ) )

       !==============================================================
       ! Copy data into fields of the Meteorology State object
       !==============================================================

       ! COS(SZA) at the current time
       State_Met%SUNCOS    (I,J) = SUNCOS

       ! COS(SZA) @ the midpoint of the current chemistry timestep
       State_Met%SUNCOSmid (I,J) = SUNCOS_MID

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE COSSZA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_met_ageofair
!
! !DESCRIPTION: Subroutine Set\_Met\_AgeOfAir adds the time step (in seconds)
!  to every grid box every time step with a total sink at the surface every
!  time step to reproduce GMI tracer mechanism.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Met_AgeOfAir( State_Grid, State_Met )
!
! !USES:
!
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_DYN
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  21 Dec 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER        :: I, J, L
    REAL(fp), SAVE :: TimeStep

    !=================================================================
    ! SET_MET_AGEOFAIR begins here!
    !=================================================================

    ! Get timestep [s]
    TimeStep = GET_TS_DYN()

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L )
    DO L=1,State_Grid%NZ
    DO J=1,State_Grid%NY
    DO I=1,State_Grid%NX

       IF ( L == 1 ) THEN
          ! Set the surface to a sink
          State_Met%AgeOfAir(I,J,L) = 0
       ELSE
          ! Otherwise add time step [s]
          State_Met%AgeOfAir(I,J,L) = State_Met%AgeOfAir(I,J,L) + TimeStep
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Set_Met_AgeOfAir
!EOC
END MODULE CALC_MET_MOD
