!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: land_mercury_mod.F90
!
! !DESCRIPTION: Module LAND\_MERCURY\_MOD contains variables and routines for
! the land emissions for the GEOS-Chem mercury simulation.
!\\
!\\
! !INTERFACE:
!
MODULE LAND_MERCURY_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD    ! For real precisions (hp)
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp, f4, f8)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: BIOMASSHG
  PUBLIC :: VEGEMIS
  PUBLIC :: SOILEMIS
  PUBLIC :: LAND_MERCURY_FLUX
  !============================================================================
  ! Disable GTMM until it is brought up-to-date (mps, 3/10/19)
  !PUBLIC :: GTMM_DR
  !============================================================================
  PUBLIC :: SNOWPACK_MERCURY_FLUX
  PUBLIC :: INIT_LAND_MERCURY
  PUBLIC :: CLEANUP_LAND_MERCURY
!
! !REVISION HISTORY:
!  02 Jun 2010 - N. E. Selin, C. Carouge - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Keep a shadow variable of N_HG_CATS for backwards compatibility
  INTEGER                :: N_HG_CATS

  ! Plant transpiration rate [m/s]
  REAL(fp),  ALLOCATABLE :: TRANSP(:,:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: land_mercury_flux
!
! !DESCRIPTION: Subroutine LAND\_MERCURY\_FLUX calculates emissions of Hg(0)
!  from prompt recycling of previously deposited mercury to land, in [kg/s].
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LAND_MERCURY_FLUX( LFLUX, LHGSNOW, State_Grid, State_Met )
!
! !USES:
!
    USE DEPO_MERCURY_MOD,   ONLY : WD_HGP, WD_HG2, DD_HGP, DD_HG2
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: LHGSNOW     ! Use Hg0 from snow?
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    ! Hg0 flux [kg/s]
    REAL(fp),       INTENT(OUT) :: LFLUX(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
!
! !REVISION HISTORY:
!  30 Aug 2010 - N. E. Selin, C. Holmes, B. Corbitt - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: DTSRCE, REEMFRAC
    REAL(fp)              :: FRAC_SNOW_OR_ICE, FRAC_SNOWFREE_LAND
    LOGICAL               :: IS_LAND_OR_ICE
    REAL(fp), PARAMETER   :: SEC_PER_YR = 365.25e+0_fp * 86400e+0_fp
    INTEGER               :: I,      J,      NN

    !=================================================================
    ! LAND_MERCURY_FLUX begins here!
    !=================================================================

    ! Emission timestep [s]
    DTSRCE = GET_TS_EMIS()

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,              J,                NN                 ) &
    !$OMP PRIVATE( REEMFRAC,       FRAC_SNOW_OR_ICE, FRAC_SNOWFREE_LAND ) &
    !$OMP PRIVATE( IS_LAND_OR_ICE                                       )
    DO J  = 1, State_Grid%NY
    DO I  = 1, State_Grid%NX
    DO NN = 1, N_Hg_CATS

       ! Distinguish between ice/snow and snow-free land
       FRAC_SNOW_OR_ICE    = MIN( State_Met%FRSNO(I,J)     + &
                                  State_Met%FRSEAICE(I,J)  + &
                                  State_Met%FRLANDIC(I,J), 1e+0_fp)
       FRAC_SNOWFREE_LAND  = MAX( State_Met%FRLAND(I,J)    - &
                                  State_Met%FRSNO(I,J),    0e+0_fp)

       IS_LAND_OR_ICE      = (( FRAC_SNOWFREE_LAND > 0e+0_fp ) .OR. &
                              ( FRAC_SNOW_OR_ICE   > 0e+0_fp ))

       ! If snow or ice on the ground, reemission fraction is 0.6,
       ! otherwise 0.2
       IF ( IS_LAND_OR_ICE ) THEN

          IF (LHGSNOW) THEN
             REEMFRAC = 0.2e+0_fp * FRAC_SNOWFREE_LAND
          ELSE
             REEMFRAC = 0.2e+0_fp * FRAC_SNOWFREE_LAND + &
                        0.6e+0_fp * FRAC_SNOW_OR_ICE
          ENDIF

          ! Mass of emitted Hg(0), kg
          LFLUX(I,J,NN) = ( WD_HgP(I,J,NN) + WD_Hg2(I,J,NN) + &
                            DD_HgP(I,J,NN) + DD_Hg2(I,J,NN) ) * REEMFRAC

          ! Emission rate of Hg(0). Convert kg /timestep -> kg/s
          LFLUX(I,J,NN) = LFLUX(I,J,NN) / DTSRCE

       ELSE

          ! No flux from non-land surfaces (water, sea ice)
          LFLUX(I,J,NN) = 0e+0_fp

       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE LAND_MERCURY_FLUX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: biomasshg
!
! !DESCRIPTION: Subroutine BIOMASSHG is the subroutine for Hg(0) emissions
!  from biomass burning. These emissions are active only for present day
!  simulations and not for preindustrial simulations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE BIOMASSHG( Input_Opt, EHg0_bb, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_ERROR_MOD
    USE HCO_STATE_MOD,      ONLY : HCO_STATE
    USE HCO_INTERFACE_MOD,  ONLY : GetHcoDiagn
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),           INTENT(IN)  :: Input_Opt   ! Input Options
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), DIMENSION(:,:), INTENT(OUT) :: EHg0_bb
    INTEGER,                  INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Emissions are based on an inventory of CO emissions from biomass burning
!  (Duncan et al. J Geophys Res 2003), multiplied by a Hg/CO ratio in BB
!  plumes from Franz Slemr (Poster, EGU 2006).
!                                                                             .
!  Slemr surveyed emission factors from measurements worldwide. Although his
!  best estimate was 1.5e-7 mol Hg/ mol CO, we chose the highest value
!  (2.1e-7 mol Hg/ mol CO) in the range because the simulations shown in
!  Selin et al. (GBC 2008) required large Hg(0) emissions to sustain
!  reasonable atmospheric Hg(0) concentrations. (eck, 11/13/2008)
!
! !REVISION HISTORY:
!  30 Jul 2008 - N. E. Selin, C. Holmes, B. Corbitt - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J

    ! Pointers
    REAL(sp), POINTER  :: Ptr2D(:,:)

    ! Strings
    CHARACTER(LEN=63)  :: DgnName
    CHARACTER(LEN=255) :: ThisLoc
    CHARACTER(LEN=512) :: ErrMsg

    !=================================================================
    ! BIOMASSHG begins here!
    !=================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    ErrMsg  =  ''
    ThisLoc =  ' -> at BIOMASSHG (in GeosCore/land_mercury_mod.F90)'
    Ptr2d   => NULL()

    ! Do biomass Hg emissions if biomass burning is on and it is a
    ! present-day simulation (i.e. not preindustrial)
    IF ( .NOT. Input_Opt%LPREINDHG ) THEN

       IF ( .NOT. ASSOCIATED(HcoState) ) THEN
          ErrMsg = 'HcoState object is not associated!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get diagnostics from HEMCO. The HG0 biomass diagnostics
       ! is defined in hcoi_gc_diagn_mod.F90 and becomes updated
       ! by HEMCO automatically. Output unit is as specified when
       ! defining diagnostics (kg/m2/s).
       DgnName = 'BIOMASS_HG0'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       EHg0_bb = Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)

    ELSE

       ! No emissions for preindustrial period, or when BB is turned off.
       EHg0_bb = 0e+0_fp

    ENDIF

    ! Free pointer
    Ptr2D    => NULL()

  END SUBROUTINE BIOMASSHG
!EOC
!-----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vegemis
!
! !DESCRIPTION: Subroutine VEGEMIS is the subroutine for Hg(0) emissions from
!  vegetation by evapotranspiration.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE VEGEMIS( Input_Opt, State_Met, LVEGEMIS, EHg0_dist, EHg0_vg, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH, ITS_A_NEW_MONTH
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),           INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(MetState),           INTENT(IN)  :: State_Met  ! Met State object
    LOGICAL,                  INTENT(IN)  :: LVEGEMIS   !
    REAL(fp), DIMENSION(:,:), INTENT(IN)  :: EHg0_dist  !
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), DIMENSION(:,:), INTENT(OUT) :: EHg0_vg
    INTEGER,                  INTENT(OUT) :: RC         ! Success or failure?
!
! !REMARKS:
!  Vegetation emissions are proportional to the evapotranspiration rate and the
!  soil water mercury content. We assume a constant concentration of mercury
!  in soil matter, based on the preindustrial and present-day simulations
!  described in Selin et al. (GBC 2008) and in SOILEMIS subroutine. From the
!  soil matter Hg concentration, we calculate a soil water Hg concentration in
!  equilibrium (Allison and Allison, 2005).
!  NASA provides a climatology of evapotranspiration based on a water budget
!  model (Mintz and Walker, 1993).
!
! Calculate vegetation emissions following Xu et al (1999)
!    Fc = Ec Cw
!
!    Fc is Hg0 flux (ng m-2 s-1)
!    Ec is canopy transpiration (m s-1)
!    Cw is conc of Hg0 in surface soil water (ng m-3)
!
! Calculate Cw from the Allison and Allison (2005) equilibrium formula
!    Cw = Cs / Kd
!
!    Cs is the concentration of Hg is surface soil solids, ng/g
!    Kd is the equilibrium constant = [sorbed]/[dissolved]
!       log Kd = 3.8 L/kg -> Kd = 6310 L /kg = 6.31D-3 m3/g
!
! We assume a global mean Cs = 45 ng/g for the preindustrial period. In
! iterative simulations we redistribute this according to the deposition
! pattern while maintining the global mean. The scaling factor, EHg0_dist,
! also accounts for the anthropogenic enhancement of soil Hg in the present
! day.

! !REVISION HISTORY:
!  30 Aug 2010 - N. Eckley, C. Holmes, B. Corbitt - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)        :: DRYSOIL_HG, SOILWATER_HG, VEG_EMIS
    INTEGER         :: I, J

    ! Soil Hg sorption to dissolution ratio, m3/g
    REAL(fp), PARAMETER  :: Kd = 6.31e-3_fp

    ! Preindustrial global mean soil Hg concentration, ng Hg /g dry soil
    REAL(fp), PARAMETER  :: DRYSOIL_PREIND_HG = 45e+0_fp

    !=================================================================
    ! VEGEMIS begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! No emissions through transpiration if we use Bess' GCAP emissions
    ! Bug fix: VEGEMIS shouldn't be tied to GCAP emissions
    ! (jaf, eds, 4/1/11)
    EHg0_vg = 0e+0_fp

  END SUBROUTINE VEGEMIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soilemis
!
! !DESCRIPTION: Subroutine SOILEMIS is the subroutine for Hg(0) emissions
!  from soils.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SOILEMIS( EHg0_dist, EHg0_so, State_Grid, State_Met )
!
! !USES:
!
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH, ITS_A_NEW_MONTH
!
! !INPUT PARAMETERS:
!
    REAL(fp), DIMENSION(:,:), INTENT(IN)  :: EHg0_dist
    TYPE(GrdState),           INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState),           INTENT(IN)  :: State_Met   ! Met State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), DIMENSION(:,:), INTENT(OUT) :: EHg0_so
!
! !REMARKS:
!  Soil emissions are a function of solar radiation at ground level
!  (accounting for attenuation by leaf canopy) and surface temperature.
!  The radiation dependence from Zhang et al. (2000) is multiplied by the
!  temperature dependence from Poissant and Casimir (1998).
!  Finally, this emission factor is multiplied by the soil mercury
!  concentration and scaled to meet the global emission total.
!
!  Comments on soil Hg concentration:
!  ----------------------------------
!  We chose the preindustrial value of 45 ng Hg /g dry soil as the mean of
!  the range quoted in Selin et al. (GBC 2008): 20-70 ng/g (Andersson, 1967;
!  Shacklette et al., 1971; Richardson et al., 2003; Frescholtz and Gustin,
!  2004). Present-day soil concentrations are thought to be 15% greater than
!  preindustrial (Mason and Sheu 2002), but such a difference is much less
!  than the range of concentrations found today, so not well constrained.
!  We calculate the present-day soil Hg distribution by adding a global mean
!  6.75 ng/g (=0.15 * 45 ng/g) according to present-day Hg deposition.
!  (eck, 11/13/08)
!
! !REVISION HISTORY:
!  30 Aug 2010 - N. Eckley, B. Corbitt - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL             :: IS_SNOWFREE_LAND
    INTEGER             :: I, J
    REAL(fp)            :: SOIL_EMIS, DIMLIGHT, TAUZ, LIGHTFRAC
    REAL(fp)            :: AREA_M2, DRYSOIL_HG
    REAL(fp)            :: SUNCOSVALUE
    REAL(fp)            :: FRAC_SNOWFREE_LAND
    REAL(fp)            :: SOIL_EMIS_FAC
!
! !DEFINED PARAMETERS:
!
    ! Preindustrial global mean soil Hg concentration, ng Hg /g dry soil
    REAL(fp), PARAMETER  :: DRYSOIL_PREIND_HG = 45e+0_fp

    ! Scaling factor for emissions, g soil /m2 /h
    ! (This parameter is beta in Eq 3 of Selin et al., GBC 2008.
    ! The value in paper is actually DRYSOIL_PREIND_HG * SOIL_EMIS_FAC
    ! and the stated units are incorrect. The paper should have stated
    ! beta = 1.5D15 / 45D0 = 3.3D13)
    ! This parameter is tuned in the preindustrial simulation
    ! so that total deposition to soil equals total emission from soil,
    ! while also requiring global mean soil Hg concentration of 45 ng/g
    ! Update for v11 GEOS-FP at 4x5 (J Fisher 3/2016)
    ! Assume same factor for all resolutions for now and therefore
    ! comment out resolution dependence. Users of 2x2.5 or nested
    ! resolutions are encouraged to test and update these as needed (J
    ! Fisher 4/2016)
    !IF      ( TRIM(State_Grid%GridRes) == '4.0x5.0') THEN
       SOIL_EMIS_FAC = 1.6e-2_fp * 0.9688e+0_fp
    !ELSE IF ( TRIM(State_Grid%GridRes) == '2.0x2.5' ) THEN
    !   SOIL_EMIS_FAC=1.6e-2_fp
    !ELSE IF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
    !   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !   ! bmy, 8/12/15, add pre-defined SOIL_EMIS_FAC for 05x0667 simulation
    !   ! This is a non-physical value and would need to be changed for an
    !   ! actual mercury simulation
    !   !
    !   ! %%% NOTE: Bob Y. used same value as for GRID05x0666  %%%
    !   ! %%% Team Hg will have to supply a better value later %%%
    !   SOIL_EMIS_FAC = 1.6e-2_fp  ! yzh
    !   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !ELSE IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
    !   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !   ! lhu, 1/5/2012, add pre-defined SOIL_EMIS_FAC for 05x0667 simulation
    !   ! This is a non-physical value and would need to be changed for an
    !   ! actual mercury simulation
    !   SOIL_EMIS_FAC = -9.0e+99_fp
    !   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !ELSE
    !   ! This sets a default value for GCM & ESMF COUPLING
    !   SOIL_EMIS_FAC = 2.4e-2_fp*.5742e+0_fp ! for sunlight function
    !ENDIF

    !=================================================================
    ! SOILEMIS begins here!
    !=================================================================

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,          J,    SOIL_EMIS                       ) &
    !$OMP PRIVATE( DRYSOIL_HG, TAUZ, LIGHTFRAC, AREA_M2, SUNCOSVALUE ) &
    !$OMP PRIVATE( IS_SNOWFREE_LAND, FRAC_SNOWFREE_LAND              )
    DO J=1, State_Grid%NY
    DO I=1, State_Grid%NX

       FRAC_SNOWFREE_LAND = MAX( State_Met%FRLAND(I,J) - &
                                 State_Met%FRSNO(I,J), 0e+0_fp )
       IS_SNOWFREE_LAND   = ( FRAC_SNOWFREE_LAND > 0e+0_fp )

       IF ( IS_SNOWFREE_LAND ) THEN

          ! attenuate solar radiation based on function of leaf area index
          ! Jacob and Wofsy 1990 equations 8 & 9
          TAUZ = State_Met%MODISLAI(I,J) * 0.5e+0_fp

          ! For very low and below-horizon solar zenith angles, use
          ! same attenuation as for SZA=85 degrees
          SUNCOSVALUE = MAX( State_Met%SUNCOS(I,J), 0.09e+0_fp )

          ! fraction of light reaching the surface is
          ! attenuated based on LAI
          LIGHTFRAC = EXP( -TAUZ / SUNCOSVALUE )

          ! Dry soil Hg concentration, ng Hg /g soil
          DRYSOIL_HG = DRYSOIL_PREIND_HG * EHg0_dist(I,J)

          ! Soil emissions, ng /m2 /h
          SOIL_EMIS = EXP( 0.0011 * State_Met%SWGDN(I,J) * LIGHTFRAC ) * &
                      DRYSOIL_HG * SOIL_EMIS_FAC

          ! Grid box surface area [m2]
          AREA_M2   = State_Grid%Area_M2(I,J)

          ! convert soilnat from ng /m2 /h -> kg /gridbox /s
          EHg0_so(I,J) = SOIL_EMIS * AREA_M2 * 1e-12_fp / &
                         ( 60e+0_fp * 60e+0_fp )

          ! Multiply by fractional land area
          EHg0_so(I,J) = EHg0_so(I,J) * FRAC_SNOWFREE_LAND

       ELSE

          ! no soil emissions from water and ice
          EHg0_so(I,J) = 0e+0_fp

       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE SOILEMIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: snowpack_mercury_flux
!
! !DESCRIPTION: Subroutine SNOWPACK\_MERCURY\_FLUX calculates emission of
!  Hg(0) from snow and ice.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SNOWPACK_MERCURY_FLUX( FLUX, LHGSNOW, State_Chm, &
                                    State_Grid, State_Met )
!
! !USES:
!
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: LHGSNOW     ! Use Hg from snow?
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    ! Hg0 flux [kg/s]
    REAL(fp),       INTENT(OUT) :: FLUX(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
!
! !REMARKS:
!  Emissions are a linear function of Hg mass stored in the snowpack. The
!  Hg lifetime in snow is assumed to be 180 d when T< 270K and 7 d when T>270K
!
!     E = k * SNOW_HG     : k = 6D-8 if T<270K, 1.6D-6 otherwise
!
!  These time constants reflect the time scales of emission observed in the
!  Arctic and in field studies. Holmes et al 2010
!
!  Formulation from Holmes et al. 2010 is now obsolete. Instead, emissions
!  from snow are tied to solar radiation, not temperature. Effective rate
!  constant is in the mid-range of values estimated by Durnford and Dastoor
!  (2011) and consistent with surface air Hg0 observations, as described in
!  Fisher et al. (2011, in review).
!
! !REVISION HISTORY:
!  15 Sep 2009 - C. Holmes, S. Carouge - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER       :: I, J, NN, MONTH
    REAL(fp)      :: DTSRCE, SNOW_HG_OC_NEW, K_EMIT, SWRAD
    REAL(fp)      :: SNOW_HG_LN_NEW, FLUX_TMP

    ! For solar flux dependent photo-reduction (jaf, 11/18/11)
    ! This is in units of m2/W/s so that when multiplied by SWRAD we
    ! end up with units of 1/s. This corresponds to k ~0.001 1/h
    ! when SWRAD = 100 W/m2
    REAL(fp), PARAMETER :: K_EMIT_0 = 2.5e-9_fp

    ! Pointers
    REAL(fp), POINTER :: SNOW_HG_OC(:,:,:)
    REAL(fp), POINTER :: SNOW_HG_LN(:,:,:)

    !=================================================================
    ! SNOWPACK_MERCURY_FLUX begins here!
    !=================================================================

    ! Initialize
    FLUX = 0e+0_fp

    ! Return to calling program if snowpack model is disabled
    IF ( .NOT. LHGSNOW ) RETURN

    ! Point to fields in State_Chm
    SNOW_HG_OC => State_Chm%SnowHgOcean
    SNOW_HG_LN => State_Chm%SnowHgLand

    ! Emission timestep [s]
    DTSRCE = GET_TS_EMIS()

    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED )   &
    !$OMP PRIVATE( I, J, NN ) &
    !$OMP PRIVATE( SNOW_HG_OC_NEW, K_EMIT, SWRAD, SNOW_HG_LN_NEW ) &
    !$OMP PRIVATE( FLUX_TMP )
    DO J  = 1, State_Grid%NY
    DO I  = 1, State_Grid%NX

       ! If the sun is set, then no emissions, go to next box
       IF ( State_Met%SUNCOS(I,J) < 0e+0_fp ) CYCLE

       ! Zero variables
       K_EMIT = 0e+0_fp
       SWRAD  = 0e+0_fp

       ! Set solar radiation @ ground
       SWRAD = State_Met%SWGDN(I,J)

       ! Compute K_EMIT
       K_EMIT = K_EMIT_0 * SWRAD

       ! SWRAD can be very small at edge of sunlit zone, probably due to
       ! averaging. This can make K_EMIT unrealistically small (e.g.:
       ! 1d-13 * 1d-9 ~ 1d-22). In IF statement below, we only use & save
       ! K_EMIT if it is >= K_EMIT_0 (jaf, 5/19/11)
       IF ( K_EMIT >= K_EMIT_0 ) THEN

          ! Loop over # of categories
          DO NN = 1, N_Hg_CATS

             ! Zero temporary reservoir on each iteration
             FLUX_TMP = 0D0

             ! Check if there is Hg that could be emitted
             IF ( SNOW_HG_OC(I,J,NN) > 0e+0_fp ) THEN

                ! New mass of snow in Hg
                SNOW_HG_OC_NEW = SNOW_HG_OC(I,J,NN) * EXP(-K_EMIT*DTSRCE)

                FLUX_TMP = MAX(SNOW_HG_OC(I,J,NN) - SNOW_HG_OC_NEW, 0e+0_fp)

                SNOW_HG_OC(I,J,NN) = SNOW_HG_OC_NEW

             ENDIF !SNOW_HG_OC > 0

             ! Check if there is Hg in land snow that could be emitted
             IF ( SNOW_HG_LN(I,J,NN) > 0e+0_fp ) THEN

                ! New mass of snow in Hg
                SNOW_HG_LN_NEW = SNOW_HG_LN(I,J,NN) * EXP(-K_EMIT*DTSRCE)

                FLUX_TMP = FLUX_TMP + MAX(SNOW_HG_LN(I,J,NN)-SNOW_HG_LN_NEW,0D0)

                SNOW_HG_LN(I,J,NN) = SNOW_HG_LN_NEW

             ENDIF !SNOW_HG_LN > 0

             ! Convert mass -> flux
             FLUX(I,J,NN) = FLUX_TMP / DTSRCE

          ENDDO !Loop over NN
       ENDIF !K_EMIT >= K_EMIT_0

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers
    SNOW_HG_OC => NULL()
    SNOW_HG_LN => NULL()

  END SUBROUTINE SNOWPACK_MERCURY_FLUX
!EOC
!==============================================================================
! Disable GTMM for now. This code needs to be brought-up-to-date so that
! input data are in netCDF and read in via HEMCO (mps, 3/10/19)
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: gtmm_dr
!!
!! !DESCRIPTION: GTMM\_DR is a driver to call GTMM from GEOS-Chem.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE GTMM_DR( Input_Opt, State_Grid, State_Met, State_Chm, Hg0gtm, RC )
!!
!! !USES:
!!
!    USE BPCH2_MOD
!    USE DEPO_MERCURY_MOD,   ONLY : CHECK_DIMENSIONS
!    USE DEPO_MERCURY_MOD,   ONLY : WD_Hg2, WD_HgP, DD_HgP, DD_Hg2
!    USE DEPO_MERCURY_MOD,   ONLY : READ_GTMM_RESTART
!    USE ErrCode_Mod
!    USE FILE_MOD,           ONLY : IOERROR
!    USE Input_Opt_Mod,      ONLY : OptInput
!    USE inquireMod,         ONLY : findFreeLun
!    USE State_Chm_Mod,      ONLY : ChmState
!    USE State_Grid_Mod,     ONLY : GrdState
!    USE State_Met_Mod,      ONLY : MetState
!    USE TIME_MOD,           ONLY : EXPAND_DATE, YMD_EXTRACT
!    USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
!!
!! !INPUT PARAMETERS:
!!
!    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!!
!! !OUTPUT PARAMETERS:
!!
!    ! Emission of Hg0 calculated by GTMM for the month [kg/s]
!    REAL(fp),       INTENT(OUT)   :: Hg0gtm(State_Grid%NX,State_Grid%NY)
!
!    ! Success or failure?
!    INTEGER,        INTENT(OUT)   :: RC
!!
!! !REMARKS:
!!  ##########################################################################
!!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!!  #####    WILL EVENTUALLY BE READ IN FROM netCDF FILES VIA HEMCO!     #####
!!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!!  ##########################################################################
!!
!! !REVISION HISTORY:
!!  15 Sep 2009 - C. Carouge  - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER   :: YEAR    ! Current year
!    INTEGER   :: MONTH   ! Current month
!    INTEGER   :: DAY
!    INTEGER   :: NYMD, NHMS
!
!    REAL(fp)  :: TSURF(State_Grid%NX, State_Grid%NY)    ! Ground temperature
!    REAL(fp)  :: PRECIP(State_Grid%NX, State_Grid%NY)   ! Total precipitation
!                                                        ! for the month
!    REAL(fp)  :: SOLAR_W(State_Grid%NX, State_Grid%NY)  ! Solar radiation for
!                                                        ! the month
!
!    REAL*4,   :: TRACER(State_Grid%NX, State_Grid%NY)   ! Temporary array
!
!    ! Monthly average deposition arrays
!    REAL(fp)  :: Hg0mth_dry(State_Grid%NX, State_Grid%NY)
!    REAL(fp)  :: Hg2mth_dry(State_Grid%NX, State_Grid%NY)
!    REAL(fp)  :: Hg2mth_wet(State_Grid%NX, State_Grid%NY)
!
!    INTEGER               :: IOS, I, J, L, IU_FILE
!
!    CHARACTER(LEN=255)    :: FILENAME
!
!    ! For binary punch file, version 2.0
!    INTEGER               :: NI,        NJ,      NL
!    INTEGER               :: IFIRST,    JFIRST,  LFIRST
!    INTEGER               :: HALFPOLAR, CENTER180
!    INTEGER               :: NTRACER,   NSKIP
!    REAL*4                :: LONRES,    LATRES
!    REAL(fp)                :: ZTAU0,     ZTAU1
!    CHARACTER(LEN=20)     :: MODELNAME
!    CHARACTER(LEN=40)     :: CATEGORY
!    CHARACTER(LEN=40)     :: UNIT
!    CHARACTER(LEN=40)     :: RESERVED
!
!    !=================================================================
!    ! GTMM_DR begins here!
!    !=================================================================
!
!    ! Assume success
!    RC         = GC_SUCCESS
!
!    ! Initialise arrays
!    NYMD       = GET_NYMD()
!    NHMS       = GET_NHMS()
!
!    TSURF      = 0e+0_fp
!    PRECIP     = 0e+0_fp
!    SOLAR_W    = 0e+0_fp
!    Hg0gtm     = 0e+0_fp
!
!    ! Reset deposition arrays.
!    Hg0mth_dry = 0e+0_fp
!    Hg2mth_dry = 0e+0_fp
!    Hg2mth_wet = 0e+0_fp
!
!    CALL YMD_EXTRACT( NYMD, YEAR, MONTH, DAY )
!
!    !=================================================================
!    ! Read monthly meteorology fields
!    !=================================================================
!
!    ! Find a free file LUN
!    IU_FILE  = findFreeLun()
!
!    ! File name
!    FILENAME = TRIM( Input_Opt%DATA_DIR ) // &
!               TRIM( Input_Opt%RES_DIR  ) // 'mercury_201007/' // &
!               'mean_metfields/' // GET_NAME_EXT() // &
!               '/mean_YYYYMM.bpch'
!
!    ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!    CALL EXPAND_DATE( FILENAME, NYMD, NHMS )
!
!    ! Echo some input to the screen
!    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!    WRITE( 6, 100   )
!    WRITE( 6, 110   ) TRIM( FILENAME ), IU_FILE
!100 FORMAT( 'G T M M  H g   M E T   F I L E   I N P U T' )
!110 FORMAT( /, 'GTMM_DR: Reading ', a, ' on unit ', i6  )
!
!    ! Open the binary punch file for input
!    CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )
!
!    !-----------------------------------------------------------------
!    ! Read concentrations -- store in the TRACER array
!    !-----------------------------------------------------------------
!    DO
!       READ( IU_FILE, IOSTAT=IOS ) &
!             MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!       ! IOS < 0 is end-of-file, so exit
!       IF ( IOS < 0 ) EXIT
!
!       ! IOS > 0 is a real I/O error -- print error message
!       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_dr:1' )
!
!       READ( IU_FILE, IOSTAT=IOS ) &
!             CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED, &
!             NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,  NSKIP
!
!       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_dr:2' )
!
!       READ( IU_FILE, IOSTAT=IOS ) &
!           ( ( TRACER(I,J), I=1,NI ), J=1,NJ )
!
!       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_dr:3' )
!
!       !--------------------------------------------------------------
!       ! Assign data from the TRACER array to the arrays.
!       !--------------------------------------------------------------
!
!       ! Process dry deposition data
!       IF ( CATEGORY(1:8) == 'GMAO-2D' ) THEN
!
!          ! Make sure array dimensions are of global size
!          ! (NI=State_Grid%NX; NJ=State_Grid%NY, NL=1), or stop the run
!          CALL CHECK_DIMENSIONS( NI, NJ, NL, State_Grid )
!
!          ! Save into arrays
!          IF ( NTRACER == 54 .OR. NTRACER == 59 ) THEN
!
!             !----------
!             ! Surface temperature
!             !----------
!
!             ! Store surface temperature in TSURF array
!             TSURF(:,:)   = TRACER(:,State_Grid%NY:1:-1)
!
!          ELSE IF ( NTRACER == 26 .OR. NTRACER == 29 ) THEN
!
!             !----------
!             ! Total precipitation
!             !----------
!
!             ! Store precipitation in PRECIP array
!             PRECIP(:,:)   = TRACER(:,State_Grid%NY:1:-1)
!
!          ELSE IF ( NTRACER == 37 .OR. NTRACER == 51 ) THEN
!
!             !----------
!             ! Solar radiation
!             !----------
!
!             ! Store solar radiation in SOLAR_W array
!             SOLAR_W(:,:) = TRACER(:,State_Grid%NY:1:-1)
!
!          ENDIF
!       ENDIF
!
!    ENDDO
!
!    ! Close file
!    CLOSE( IU_FILE )
!
!    !=================================================================
!    ! Read GTMM restart file to get data from previous month
!    !=================================================================
!    CALL READ_GTMM_RESTART( Input_Opt,  State_Chm, State_Grid, &
!                            NYMD,       NHMS,                  &
!                            Hg0mth_dry, Hg2mth_dry, Hg2mth_wet, RC  )
!
!    !=================================================================
!    ! Call GTMM model
!    !=================================================================
!#if defined( GTMM_Hg )
!    CALL GTMM_coupled( YEAR,    &
!                       MONTH,   &
!                       Hg0mth_dry(:,State_Grid%NY:1:-1), &
!                       Hg2mth_dry(:,State_Grid%NY:1:-1), &
!                       Hg2mth_wet(:,State_Grid%NY:1:-1), &
!                       TSURF,   &
!                       PRECIP,  &
!                       SOLAR_W, &
!                       Hg0gtm(:,State_Grid%NY:1:-1) )
!#endif
!
!    DO J = 1, State_Grid%NY
!    DO I = 1, State_Grid%NX
!       IF ( .NOT. State_Met%IsLand(I,J) ) Hg0gtm(I,J) = 0e+0_fp
!    ENDDO
!    ENDDO
!
!  END SUBROUTINE GTMM_DR
!!EOC
!==============================================================================
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_land_mercury
!
! !DESCRIPTION: Subroutine INIT\_LAND\_MERCURY allocates and zeroes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_LAND_MERCURY( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
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
!  14 Sep 2009 - C. Carouge  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! INIT_MERCURY begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Save State_Chm%N_Hg_CATS in a shadow variable
    N_Hg_CATS = State_Chm%N_Hg_CATS

  END SUBROUTINE INIT_LAND_MERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_land_mercury
!
! !DESCRIPTION: Subroutine CLEANUP\_LAND\_MERCURY deallocates all module
!  arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_LAND_MERCURY
!
! !REVISION HISTORY:
!  14 Sep 2009 - C. Carouge  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

  END SUBROUTINE CLEANUP_LAND_MERCURY
!EOC
END MODULE LAND_MERCURY_MOD
