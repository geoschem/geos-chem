!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: co2_mod.F90
!
! !DESCRIPTION: Module CO2\_MOD contains variables and routines used for the
!  CO2 simulation.  A tagged CO2 simulation capability has now been added.
!\\
!\\
!  References:
!
!  \begin{itemize}
!  \item Andres, R.J, G. Marland, I. Fung, and E. Matthews, \emph{A 1x1
!        distribution of carbon dioxide emissions from fossil fuel
!        consumption and cement manufacture}, \underline{Glob. Biogeochem.
!        Cycles}, \textbf{10}, 419-429, 1996.
!  \item Corbett and Koehler (2003) \emph{Updated emissions from ocean
!        shipping}, \underline{J. Geophys. Res.}, \textbf{108}, D20, 4650.
!  \item Corbett and Koehler (2004) \emph{Considering alternative input
!        parameters in an activity-based ship fuel consumption and emissions
!        model: Reply ...} \underline{J. Geophys. Res.}, D23303.
!  \item Endresen et al. (2007) \emph{A historical reconstruction of ships
!        fuel consumption and emissions}, \underline{J. Geophys. Res.}
!        \textbf{112}, D12301.
!  \item Kim et al. (2005) \emph{System for assessing Aviation's Global
!        Emissions (SAGE) Version 1.5 global Aviation Emissions Inventories
!        for 2000-2004}
!  \item Kim et al. (2007) \emph{System for assessing Aviation's Global
!        Emissions (SAGE) Part 1: Model description and inventory results}
!  \item LeQuere et al. (2009) \emph{Trends in the sources and sinks of carbon
!        dioxide}, \underline{Nature Geoscience}, doi:10.1038/ngeo689.
!  \item Olsen and Randerson (2004), \emph{Differences between surface and
!        column atmospheric CO2 and implications for carbon cycle research},
!        \underline{J. Geophys. Res.}, \textbf{109}, D02301,
!  \item Potter et al. (1993), \emph{Terrestrial Ecosystem Production:
!        A process model based on global satellite and surface data},
!        \underline{Glob. Biogeochem. Cycles}, \textbf{7}(4), 811-841.
!  \item Randerson, J.T, M.V. Thompson, T.J.Conway, I.Y. Fung, and C.B. Field,
!        \emph{The contribution of terrestrial sources and sinks to trends
!        in the seasonal cycle of atmospheric carbon dioxide},
!        \underline{Glob. Biogeochem. Cycles},\textbf{11}, 535-560, 1997.
!  \item Suntharalingam et al. (2005) \emph{Infulence of reduced carbon
!        emissions and oxidation on the distribution of atmospheric CO2:
!        Implications for inversion analysis}, BGC, 19, GB4003.
!  \item Takahashi, T, R. Feely, R. Weiss, R. Wanninkof, D. Chipman,
!        S. Sutherland, and T. Takahashi (1997), \emph{Global air-sea flux
!        of CO2: An estimate based on measurements of sea-air pCO2 difference},
!        \underline{Proceedings of the National Academy of Sciences},
!        \textbf{94}, 8292-8299.
!  \item Takahashi, T, et al. (2009), \emph{Climatological mean and decadal
!        change in surface ocean pCO2, and net sea-air CO2 flux over the
!        global oceans}, \textbf{Deep-Sea Research II},
!        doi:10.1016/jdsr2/2008.12.009.
!  \item Yevich, R. and J. A. Logan, \emph{An assesment of biofuel use and
!        burning of agricultural waste in the developing world},
!        \underline{Glob. Biogeochem. Cycles}, \textbf{17}, 1095,
!        doi:10.1029/2002GB001952, 2003.
!  \item Sausen, R. and Schumann, U. "Estimates of the Climate Response to
!        Aircraft CO2 and NOx Emissions Scenarios", Climate Change,
!        44: 27-58, 2000
!  \item Wilkersen, J.T. et al. \emph{Analysis of emission data from global
!        commercial Aviation: 2004 and 2006}, \underline{Atmos. chem. Phys.
!        Disc.}, \textbf{10}, 2945-2983, 2010.
!  \end{itemize}
!
! !INTERFACE:
!
MODULE CO2_MOD
!
! !USES:
!
  USE PhysConstants       ! Physical constants
  USE inquireMod,    ONLY : findFreeLUN
  USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CLEANUP_CO2
  PUBLIC  :: INIT_CO2
  PUBLIC  :: EMISSCO2
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DEF_BIOSPH_CO2_REGIONS_F
  PRIVATE :: DEF_OCEAN_CO2_REGIONS_F
  PRIVATE :: DEF_FOSSIL_CO2_REGIONS_F
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  WARNING! Tagged CO2 simulation only work for 2 x 2.5 grid! %%%
!  %%%  Someone will have to make this more general later on...    %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                             .
! !REVISION HISTORY:
!  16 Aug 2005 - P. Suntharalingam - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER, ALLOCATABLE :: FOSSIL_REGION(:,:)
  INTEGER, ALLOCATABLE :: BIOSPH_REGION(:,:)
  INTEGER, ALLOCATABLE :: OCEAN_REGION(:,:)
!
! !DEFINED PARAMETERS:
!
  ! FMOL_CO2     - kg CO2 / mole CO2
  REAL(fp),  PARAMETER   :: FMOL_CO2   = 44e-3_fp

  ! FMOL_C       - kg C   / mole C
  REAL(fp),  PARAMETER   :: FMOL_C     = 12e-3_fp

  ! XNUMOL_CO2   - molecules CO2 / kg CO2
  REAL(fp),  PARAMETER   :: XNUMOL_CO2 = AVO / FMOL_CO2

  ! XNUMOL_C     - molecules C / kg C
  REAL(fp),  PARAMETER   :: XNUMOL_C = AVO / FMOL_C

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissco2
!
! !DESCRIPTION: Subroutine EMISSCO2 is the driver routine for CO2 emissions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSCO2( Input_Opt,  State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE HCO_STATE_MOD,      ONLY : HCO_GetHcoID
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  The initial condition for CO2 has to be at least 50 ppm or higher or else
!  the balanced biosphere fluxes will make STT negative. (pns, bmy, 8/16/05)
!
! !REVISION HISTORY:
!  16 Aug 2005 - P. Suntharalingam - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                    :: I, J, L, N, NA
    INTEGER                    :: nAdvect
    REAL(fp)                   :: A_CM2, DTSRCE, E_CO2

    ! SAVEd scalars
    LOGICAL,         SAVE      :: FIRST      = .TRUE.

    ! Strings
    CHARACTER(LEN=255)         :: ThisLoc
    CHARACTER(LEN=512)         :: ErrMsg

    ! Pointers
    REAL(f4),        POINTER   :: CO2_COPROD(:,:,:  )
    REAL(fp),        POINTER   :: Spc       (:,:,:,:)
!
! !DEFINED PARAMETERS:
!
    REAL(fp),        PARAMETER :: CM2PERM2 = 1.d4
    REAL(fp),        PARAMETER :: CM3PERM3 = 1.d6

    !=================================================================
    ! EMISSCO2 begins here!
    !=================================================================

    ! Initialize
    RC          = GC_SUCCESS
    ErrMsg      = ''
    ThisLoc     = ' -> at EMISSCO2 (in module GeosCore/co2_mod.F)'
    CO2_COPROD  => NULL()
    Spc         => NULL()

    ! Import emissions from HEMCO (through HEMCO state)
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       ErrMsg = 'The "HcoState" object is not defined!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Emission timestep
    DTSRCE = HcoState%TS_EMIS

    ! Number of advected species
    nAdvect = State_Chm%nAdvect

    !=================================================================
    ! Species ID setup and error checks (first-time only)
    !=================================================================
    IF ( FIRST ) THEN

       !--------------------------------------------------------------
       ! Error check: For now, the emission grid must be
       ! on the simulation grid.
       !--------------------------------------------------------------
       IF ( HcoState%NX /= State_Grid%NX .OR. &
            Hcostate%NY /= State_Grid%NY .OR. &
            Hcostate%NZ /= State_Grid%NZ     ) THEN
          ErrMsg = 'The HEMCO grid not same as the sim. grid!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Set first-time flag to false
       FIRST = .FALSE.
    ENDIF

    !=================================================================
    ! Process emissions and save diagnostics
    ! #10 CO2 production from CO oxidation
    !
    ! NOTE: Emissions for all other tagged species are now handled
    ! directly via HEMCO.  Therefore, we can move the IF statement
    ! out of the DO loop, so that the DO loop will only execute
    ! if Input_Opt%LCHEMCO2 is TRUE.  We can also move the L-loop
    ! to the outermost loop, which is more efficient.
    !=================================================================
    IF ( Input_Opt%LCHEMCO2 ) THEN

       ! Point to chemical species array [kg/kg dry air]
       Spc => State_Chm%Species

       ! %%% NOTE: This might be able to be done just in init %%%
       ! %%% as the HEMCO pointer will evolve in time %%%%%%%%%%%
       ! Get a pointer to the CO2 production from CO oxidation
       ! This is now listed in the NON-EMISSIONS DATA section of
       ! the HEMCO configuration file. (bmy, 4/17/15)
       CALL HCO_GetPtr( HcoState, 'CO2_COPROD', CO2_COPROD, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'CO2 production is not defined in HEMCO!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Loop over all grid boxes
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, E_CO2, N )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Production is in [kg C/m3], convert to [molec/cm2/s]
          E_CO2 = CO2_COPROD(I,J,L) &                      ! kg/m3
                  / CM3PERM3        &                      ! => kg/cm3
                  * XNUMOL_C        &                      ! => molec/cm3
                  / DTSRCE          &                      ! => molec/cm3/s
                  *State_Met%BXHEIGHT(I,J,L) * 100         ! => molec/cm2/s

          !==========================================================
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Save production of CO2 from CO oxidation [kg/m2/s]
          !==========================================================
          IF ( State_Diag%Archive_ProdCO2fromCO ) THEN
             State_Diag%ProdCO2fromCO(I,J,L) = E_CO2        & ! molec/cm2/s
                                               / XNUMOL_CO2 & ! => kg/cm2/s
                                               * 1e4_fp       ! => kg/m2/s

          ENDIF

          ! Convert emissions from [molec/cm2/s] to [kg/kg dry air]
          ! (ewl, 9/11/15)
          E_CO2  =  E_CO2 * DTSRCE * CM2PERM2 / &
                    ( XNUMOL_CO2 * State_Met%DELP(I,J,L) &
                    * G0_100 * ( 1.0e+0_fp &
                    - State_Met%SPHU(I,J,L) * 1.0e-3_fp ) )

          ! Add to Species #1: Total CO2 [kg/kg]
          Spc(I,J,L,1) = Spc(I,J,L,1) + E_CO2

          ! Add to Species #10: Chemical Source of CO2 [kg/kg]
          IF ( nAdvect > 9 ) THEN
             Spc(I,J,L,10) = Spc(I,J,L,10) + E_CO2
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! Free pointers
    Spc        => NULL()
    CO2_COPROD => NULL()

  END SUBROUTINE EMISSCO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: def_biosph_co2_regions_f
!
! !DESCRIPTION: Subroutine DEF\_BIOSPH\_CO2\_REGIONS defines the land
!  biospheric and ocean CO2 exchange regions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DEF_BIOSPH_CO2_REGIONS_F( State_Grid, REGION )
!
! !USES:
!
    USE FILE_MOD,       ONLY : IOERROR
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: REGION(State_Grid%NX,State_Grid%NY)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  WARNING! Tagged CO2 simulation only work for 2 x 2.5 grid!  %%%
!  %%%  Someone will have to make this more general later on...     %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, IOS, IU_FILE
    INTEGER                :: TMP(State_Grid%NX,State_Grid%NY)
    INTEGER                :: LAND_REG(State_Grid%NX,State_Grid%NY)
    CHARACTER(LEN=255)     :: LANDFILE
    CHARACTER(LEN=144)     :: ROW
    CHARACTER(LEN=1)       :: CHAR1(State_Grid%NX,State_Grid%NY)

    !=================================================================
    ! Reading LAND BIOSPHERE REGIONS
    !=================================================================

    LANDFILE  = 'Regions_land.dat'

    WRITE(*,*) ' '
100 FORMAT( '     - READ_REGIONS: Reading ', a )
    WRITE( 6, 100 ) TRIM( LANDFILE )

    ! Initialize ARRAY
    LAND_REG = 0

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open file
    OPEN( IU_FILE, FILE = TRIM( LANDFILE ), FORM='FORMATTED', IOSTAT=IOS )
    IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:1' )

    ! Read data
    DO J = 1, State_Grid%NY
       IF (State_Grid%NX ==  72) READ( IU_FILE, '(72A)', IOSTAT=IOS ) ROW
       IF (State_Grid%NX == 144) READ( IU_FILE,'(144A)', IOSTAT=IOS ) ROW
       WRITE (*,'(A)') ROW

       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:2' )

       DO I = 1, State_Grid%NX
          CHAR1(I,J) = ROW(I:I)
          IF (CHAR1(I,J) == ' ') CHAR1(I,J) = '0'
          READ (CHAR1(I,J),'(I1)') TMP(I,J)
       ENDDO
    ENDDO

    ! Close file
    CLOSE( IU_FILE )

    ! Flip array in the North-South Direction
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       LAND_REG(I,J) = TMP(I,State_Grid%NY-J+1)
    ENDDO
    ENDDO
    WRITE(*,*) ' '

    !=================================================================
    ! Loop over entire globe -- multiprocessor
    !=================================================================

    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       !-----------------------------------------------------------------------
       ! Species #13 -- Canadian Tundra
       IF (LAND_REG(I,J) == 1 .and. I > 5 .and. I <= 60) THEN
          REGION(I,J) = 13
       !----------------------------------------------------------------------
       ! Species #14 -- NA Boreal Forest
       ELSE IF (LAND_REG(I,J) == 2 .and. I <= 60) THEN
          REGION(I,J) = 14
       !-----------------------------------------------------------------------
       ! Species #15 -- Western US/Mexico
       ELSE IF (LAND_REG(I,J) == 3 .and. I <= 60) THEN
          REGION(I,J) = 15
       !-----------------------------------------------------------------------
       ! Species #16 -- Central NA Agricultural
       ELSE IF (LAND_REG(I,J) == 4 .and. I <= 60) THEN
          REGION(I,J) = 16
       !-----------------------------------------------------------------------
       ! Species #17 -- NA Mixed Forest
       ELSE IF (LAND_REG(I,J) == 5 .and. I <= 60) THEN
          REGION(I,J) = 17
       !-----------------------------------------------------------------------
       ! Species #18 -- Central America and Caribbean
       ELSE IF (LAND_REG(I,J) == 6 .and. I <= 60) THEN
          REGION(I,J) = 18
       !-----------------------------------------------------------------------
       ! Species #19 -- SA Tropical Rain Forest
       ELSE IF (LAND_REG(I,J) == 7 .and. I <= 60) THEN
          REGION(I,J) = 19
       !-----------------------------------------------------------------------
       ! Species #20 -- SA Coast and Mountains
       ELSE IF (LAND_REG(I,J) == 8 .and. I <= 60) THEN
          REGION(I,J) = 20
       !-----------------------------------------------------------------------
       ! Species #21 -- SA Wooded Grasslands
       ELSE IF (LAND_REG(I,J) == 9 .and. I <= 60) THEN
          REGION(I,J) = 21
       !-----------------------------------------------------------------------
       ! Species #22 -- Eurasian Tundra
       ELSE IF (LAND_REG(I,J) == 1 .and. (I>60 .or. I<=5)) THEN
          REGION(I,J) = 22
       !-----------------------------------------------------------------------
       ! Species #23 -- Eurasian Boreal Coniferous Forest
       ELSE IF (LAND_REG(I,J) == 2 .and. I > 60 .and. J > 65) THEN
          REGION(I,J) = 23
       !-----------------------------------------------------------------------
       ! Species #24 -- Eurasian Boreal Deciduous Forest
       ELSE IF (LAND_REG(I,J) == 5 .and. I > 60 .and. J > 65) THEN
          REGION(I,J) = 24
       !-----------------------------------------------------------------------
       ! Species #25 -- South and Central Europe
       ELSE IF (LAND_REG(I,J) == 6 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 25
       !-----------------------------------------------------------------------
       ! Species #26 -- Central Asian Grasslands
       ELSE IF (LAND_REG(I,J) == 4 .and. I > 60 .and. J > 46) THEN
          REGION(I,J) = 26
       !-----------------------------------------------------------------------
       ! Species #27 -- Central Asian Desert
       ELSE IF (LAND_REG(I,J) == 8 .and. I >100 .and. I <118) THEN
          REGION(I,J) = 27
       !-----------------------------------------------------------------------
       ! Species #28 -- East Asia Mainland
       ELSE IF (LAND_REG(I,J) == 3 .and. I > 100) THEN
          REGION(I,J) = 28
       !-----------------------------------------------------------------------
       ! Species #29 -- Japan
       ELSE IF (LAND_REG(I,J) == 9 .and. I > 100) THEN
          REGION(I,J) = 29
       !-----------------------------------------------------------------------
       ! Species #30 -- Northern African Desert
       ELSE IF (LAND_REG(I,J) == 8 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 30
       !-----------------------------------------------------------------------
       ! Species #31 -- Northern Africa Grasslands
       ELSE IF (LAND_REG(I,J) == 3 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 31
       !-----------------------------------------------------------------------
       ! Species #32 -- Africa Tropical Forest
       ELSE IF (LAND_REG(I,J) == 7 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 32
       !-----------------------------------------------------------------------
       ! Species #33 -- Southern Africa Grasslands
       ELSE IF (LAND_REG(I,J) == 4 .and. I > 60 .and. J < 50) THEN
          REGION(I,J) = 33
       !-----------------------------------------------------------------------
       ! Species #34 -- Southern African Desert
       ELSE IF (LAND_REG(I,J) == 9 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 34
       !-----------------------------------------------------------------------
       ! Species #35 -- Middle East
       ELSE IF (LAND_REG(I,J) == 2 .and. J > 40 .and. J < 65) THEN
          REGION(I,J) = 35
       !-----------------------------------------------------------------------
       ! Species #36 -- India and bordering countries
       ELSE IF (LAND_REG(I,J) == 5 .and. I > 60 .and. J < 65) THEN
          REGION(I,J) = 36
       !-----------------------------------------------------------------------
       ! Species #37 -- Maritime Asia (Indonesia, Malaysia, New Guinea, etc.)
       ELSE IF (LAND_REG(I,J) == 7 .and. I > 100) THEN
          REGION(I,J) = 37
       !-----------------------------------------------------------------------
       ! Species #38 -- Australian Forest/Grassland
       ELSE IF (LAND_REG(I,J) == 6 .and. I > 100) THEN
          REGION(I,J) = 38
       !-----------------------------------------------------------------------
       ! Species #39 -- Australian Desert
       ELSE IF (LAND_REG(I,J) == 8 .and. I >116 .and. J < 46) THEN
          REGION(I,J) = 39
       !-----------------------------------------------------------------------
       ! Species #40 -- New Zealand
       ELSE IF (LAND_REG(I,J) == 2 .and. I > 120) THEN
          REGION(I,J) = 40
       !-----------------------------------------------------------------------
       ! Species #52 -- CO2 from everywhere else (Remote Islands & Ice Caps)
       ELSE
          REGION(I,J) = 52
       !-----------------------------------------------------------------------
       ENDIF
    ENDDO
    ENDDO

  END SUBROUTINE DEF_BIOSPH_CO2_REGIONS_F
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: def_ocean_co2_regions_f
!
! !DESCRIPTION: Subroutine DEF\_OCEAN\_CO2\_REGIONS defines CO2 regions
!  for ocean exchange.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DEF_OCEAN_CO2_REGIONS_F( State_Grid, REGION )
!
! !USES:
!
    USE FILE_MOD,       ONLY : IOERROR
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: REGION(State_Grid%NX,State_Grid%NY)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  WARNING! Tagged CO2 simulation only work for 2 x 2.5 grid!  %%%
!  %%%  Someone will have to make this more general later on...     %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    INTEGER               :: I, J, IU_FILE, IOS
    INTEGER               :: TMP(State_Grid%NX,State_Grid%NY)
    INTEGER               :: OCEAN_REG(State_Grid%NX,State_Grid%NY)
    CHARACTER(LEN=255)    :: OCEANFILE
    CHARACTER(LEN=144)    :: ROW
    CHARACTER(LEN=1)      :: CHAR1(State_Grid%NX,State_Grid%NY)

    !=================================================================
    ! DEF_CO2_OCEAN_REGIONS begins here!
    !=================================================================

    OCEANFILE = 'Regions_ocean.dat'

    WRITE( 6, 100 ) TRIM( OCEANFILE )
100 FORMAT( '     - READ_REGIONS: Reading ', a )
    WRITE(*,*) ' '

    ! Initialize ARRAYS
    OCEAN_REG = 0

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open file
    OPEN( IU_FILE, FILE = TRIM( OCEANFILE ), FORM='FORMATTED', IOSTAT=IOS )
    IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:1' )

    ! Read data
    DO J = 1, State_Grid%NY
       IF (State_Grid%NX ==  72) READ( IU_FILE, '(72A)', IOSTAT=IOS ) ROW
       IF (State_Grid%NX == 144) READ( IU_FILE,'(144A)', IOSTAT=IOS ) ROW
       WRITE (*,'(A)') ROW

       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:2' )

       DO I = 1, State_Grid%NX
          CHAR1(I,J) = ROW(I:I)
          IF (CHAR1(I,J) == ' ') CHAR1(I,J) = '0'
          READ (CHAR1(I,J),'(I1)') TMP(I,J)
       ENDDO
    ENDDO

    ! Close file
    CLOSE( IU_FILE )

    ! Flip array in the North-South Direction
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       OCEAN_REG(I,J) = TMP(I,State_Grid%NY-J+1)
    ENDDO
    ENDDO
    WRITE(*,*) ' '

    !=================================================================
    ! Loop over entire globe -- multiprocessor
    !=================================================================

    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       !-----------------------------------------------------------------------
       ! Species #41 -- Arctic Ocean
       IF       (OCEAN_REG(I,J) == 5 .and. J > 60) THEN
          REGION(I,J) = 41
       !-----------------------------------------------------------------------
       ! Species #42 -- North Pacific
       ELSE IF  (OCEAN_REG(I,J) == 1) THEN
          REGION(I,J) = 42
       !-----------------------------------------------------------------------
       ! Region #43 -- Tropical West Pacific
       ELSE IF  (OCEAN_REG(I,J) == 2) THEN
          REGION(I,J) = 43
       !-----------------------------------------------------------------------
       ! Species #44 -- Tropical East Pacific
       ELSE IF  (OCEAN_REG(I,J) == 3) THEN
          REGION(I,J) = 44
       !-----------------------------------------------------------------------
       ! Species #45-- South Pacific
       ELSE IF  (OCEAN_REG(I,J) == 4) THEN
          REGION(I,J) = 45
       !-----------------------------------------------------------------------
       ! Species #46 -- North Atlantic
       ELSE IF  (OCEAN_REG(I,J) == 6 .and. J > 45) THEN
          REGION(I,J) = 46
       !-----------------------------------------------------------------------
       ! Species #47 -- Tropical Atlantic
       ELSE IF  (OCEAN_REG(I,J) == 7) THEN
          REGION(I,J) = 47
       !-----------------------------------------------------------------------
       ! Species #48 -- South Atlantic
       ELSE IF  (OCEAN_REG(I,J) == 8) THEN
          REGION(I,J) = 48
       !-----------------------------------------------------------------------
       ! Species #49 -- Tropical Indian Ocean
       ELSE IF  (OCEAN_REG(I,J) == 5 .and. J < 60) THEN
          REGION(I,J) = 49
       !-----------------------------------------------------------------------
       ! Species #50 -- Southern Indian Ocean
       ELSE IF  (OCEAN_REG(I,J) == 6 .and. J < 45) THEN
          REGION(I,J) = 50
       !-----------------------------------------------------------------------
       ! Species #51 -- Southern (Antacrtic) Ocean
       ELSE IF  (OCEAN_REG(I,J) == 9) THEN
          REGION(I,J) = 51
       !-----------------------------------------------------------------------
       ! Species #52 -- CO2 from everywhere else (Remote Islands & Ice Caps)
       ELSE
          REGION(I,J) = 52
       !-----------------------------------------------------------------------
       ENDIF
    ENDDO
    ENDDO

  END SUBROUTINE DEF_OCEAN_CO2_REGIONS_F
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: def_fossil_co2_regions_f
!
! !DESCRIPTION:  Subroutine DEF\_FOSSIL\_CO2\_REGIONS defines CO2 regions
!  for anthropogenic emissions
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DEF_FOSSIL_CO2_REGIONS_F( State_Grid, REGION )
!
! !USES:
!
    USE FILE_MOD,       ONLY : IOERROR
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: REGION(State_Grid%NX,State_Grid%NY)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  WARNING! Tagged CO2 simulation only work for 2 x 2.5 grid!  %%%
!  %%%  Someone will have to make this more general later on...     %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, IU_FILE, IOS
    INTEGER                :: TMP(State_Grid%NX,State_Grid%NY)
    INTEGER                :: REG_CODE(State_Grid%NX,State_Grid%NY)
    CHARACTER(LEN=255)     :: FILENAME
    CHARACTER(LEN=144)     :: ROW
    CHARACTER(LEN=1)       :: CHAR1(State_Grid%NX,State_Grid%NY)

    !=================================================================
    ! DEF_CO2_FOSSIL_REGIONS begins here!
    !=================================================================

    FILENAME  = 'Regions_land.dat'

    WRITE( 6, 100 ) TRIM( FILENAME )
100 FORMAT( '     - READ_REGIONS: Reading ', a )

    ! Initialize ARRAYS
    REG_CODE = 0

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open file
    OPEN( IU_FILE, FILE = TRIM( FILENAME ), FORM='FORMATTED', IOSTAT=IOS )
    IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:1' )

    ! Read data
    DO J = 1, State_Grid%NY
       IF (State_Grid%NX ==  72) READ( IU_FILE, '(72A)', IOSTAT=IOS ) ROW
       IF (State_Grid%NX == 144) READ( IU_FILE,'(144A)', IOSTAT=IOS ) ROW
       WRITE (*,'(A)') ROW

       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:2' )

       DO I = 1, State_Grid%NX
          CHAR1(I,J) = ROW(I:I)
          IF (CHAR1(I,J) == ' ') CHAR1(I,J) = '0'
          READ (CHAR1(I,J),'(I1)') TMP(I,J)
       ENDDO
    ENDDO

    ! Close file
    CLOSE( IU_FILE )

    ! Flip array in the North-South Direction
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       REG_CODE(I,J) = TMP(I,State_Grid%NY-J+1)
    ENDDO
    ENDDO

    !=================================================================
    ! Loop over entire globe -- multiprocessor
    !=================================================================
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       !-----------------------------------------------------------------------
       ! Species #13 -- Canadian Tundra
       IF (REG_CODE(I,J) == 1 .and. I > 5 .and. I <= 60) THEN
          REGION(I,J) = 13
       !-----------------------------------------------------------------------
       ! Species #14 -- NA Boreal Forest
       ELSE IF (REG_CODE(I,J) == 2 .and. I <= 60) THEN
          REGION(I,J) = 14
       !-----------------------------------------------------------------------
       ! Species #15 -- Western US/Mexico
       ELSE IF (REG_CODE(I,J) == 3 .and. I <= 60) THEN
          REGION(I,J) = 15
       !-----------------------------------------------------------------------
       ! Species #16 -- Central NA Agricultural
       ELSE IF (REG_CODE(I,J) == 4 .and. I <= 60) THEN
          REGION(I,J) = 16
       !-----------------------------------------------------------------------
       ! Species #17 -- NA Mixed Forest
       ELSE IF (REG_CODE(I,J) == 5 .and. I <= 60) THEN
          REGION(I,J) = 17
       !-----------------------------------------------------------------------
       ! Species #18 -- Central America and Caribbean
       ELSE IF (REG_CODE(I,J) == 6 .and. I <= 60) THEN
          REGION(I,J) = 18
       !-----------------------------------------------------------------------
       ! Species #19 -- SA Tropical Rain Forest
       ELSE IF (REG_CODE(I,J) == 7 .and. I <= 60) THEN
          REGION(I,J) = 19
       !-----------------------------------------------------------------------
       ! Species #20 -- SA Coast and Mountains
       ELSE IF (REG_CODE(I,J) == 8 .and. I <= 60) THEN
          REGION(I,J) = 20
       !-----------------------------------------------------------------------
       ! Species #21 -- SA Wooded Grasslands
       ELSE IF (REG_CODE(I,J) == 9 .and. I <= 60) THEN
          REGION(I,J) = 21
       !-----------------------------------------------------------------------
       ! Species #22 -- Eurasian Tundra
       ELSE IF (REG_CODE(I,J) == 1 .and. (I>60 .or. I<=5)) THEN
          REGION(I,J) = 22
       !-----------------------------------------------------------------------
       ! Species #23 -- Eurasian Boreal Coniferous Forest
       ELSE IF (REG_CODE(I,J) == 2 .and. I > 60 .and. J > 55) THEN
          REGION(I,J) = 23
       !-----------------------------------------------------------------------
       ! Species #24 -- Eurasian Boreal Deciduous Forest
       ELSE IF (REG_CODE(I,J) == 5 .and. I > 60 .and. J > 64) THEN
          REGION(I,J) = 24
       !-----------------------------------------------------------------------
       ! Species #25 -- South and Central Europe
       ELSE IF (REG_CODE(I,J) == 6 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 25
       !-----------------------------------------------------------------------
       ! Species #26 -- Central Asian Grasslands
       ELSE IF (REG_CODE(I,J) == 4 .and. I > 60 .and. J > 46) THEN
          REGION(I,J) = 26
       !-----------------------------------------------------------------------
       ! Species #27 -- Central Asian Desert
       ELSE IF (REG_CODE(I,J) == 8 .and. I >100 .and. I <118) THEN
          REGION(I,J) = 27
       !-----------------------------------------------------------------------
       ! Species #28 -- East Asia Mainland
       ELSE IF (REG_CODE(I,J) == 3 .and. I > 100) THEN
          REGION(I,J) = 28
       !-----------------------------------------------------------------------
       ! Species #29 -- Japan
       ELSE IF (REG_CODE(I,J) == 9 .and. I > 100) THEN
          REGION(I,J) = 29
       !-----------------------------------------------------------------------
       ! Species #30 -- Northern African Desert
       ELSE IF (REG_CODE(I,J) == 8 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 30
       !-----------------------------------------------------------------------
       ! Species #31 -- Northern Africa Grasslands
       ELSE IF (REG_CODE(I,J) == 3 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 31
       !-----------------------------------------------------------------------
       ! Species #32 -- Africa Tropical Forest
       ELSE IF (REG_CODE(I,J) == 7 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 32
       !-----------------------------------------------------------------------
       ! Species #33 -- Southern Africa Grasslands
       ELSE IF (REG_CODE(I,J) == 4 .and. I > 60 .and. J < 50) THEN
          REGION(I,J) = 33
       !-----------------------------------------------------------------------
       ! Species #34 -- Southern African Desert
       ELSE IF (REG_CODE(I,J) == 9 .and. I > 60 .and. I <100) THEN
          REGION(I,J) = 34
       !-----------------------------------------------------------------------
       ! Species #35 -- Middle East
       ELSE IF (REG_CODE(I,J) == 2 .and. J > 40 .and. J < 60) THEN
          REGION(I,J) = 35
       !-----------------------------------------------------------------------
       ! Species #36 -- India and bordering countries
       ELSE IF (REG_CODE(I,J) == 5 .and. I > 60 .and. J < 64) THEN
          REGION(I,J) = 36
       !-----------------------------------------------------------------------
       ! Species #37 -- Maritime Asia (Indonesia, Malaysia, New Guinea, etc.)
       ELSE IF (REG_CODE(I,J) == 7 .and. I > 100) THEN
          REGION(I,J) = 37
       !-----------------------------------------------------------------------
       ! Species #38 -- Australian Forest/Grassland
       ELSE IF (REG_CODE(I,J) == 6 .and. I > 100) THEN
          REGION(I,J) = 38
       !-----------------------------------------------------------------------
       ! Species #39 -- Australian Desert
       ELSE IF (REG_CODE(I,J) == 8 .and. I > 116 .and. J <45) THEN
          REGION(I,J) = 39
       !-----------------------------------------------------------------------
       ! Species #40 -- New Zealand
       ELSE IF (REG_CODE(I,J) == 2 .and. I > 120) THEN
          REGION(I,J) = 40
       !-----------------------------------------------------------------------
       ! Species #52 -- CO2 from everywhere else (Remote Islands & Ice Caps)
       ELSE
          REGION(I,J) = 52
       !-----------------------------------------------------------------------
       ENDIF
    ENDDO
    ENDDO

  END SUBROUTINE DEF_FOSSIL_CO2_REGIONS_F
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_co2
!
! !DESCRIPTION: Subroutine INIT\_CO2 allocates memory to module arrays and
!  reads in annual mean emissions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_CO2( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE    :: IS_INIT = .FALSE.
    INTEGER          :: AS

    ! For values from Input_Opt
    LOGICAL          :: LFOSSIL
    LOGICAL          :: LCHEMCO2
    LOGICAL          :: LBIODIURNAL
    LOGICAL          :: LBIONETCLIM
    LOGICAL          :: LOCEAN
    LOGICAL          :: LSHIP
    LOGICAL          :: LPLANE
    LOGICAL          :: LFFBKGRD
    LOGICAL          :: LBIOSPHTAG,  LFOSSILTAG,  LSHIPTAG
    LOGICAL          :: LPLANETAG

    !=================================================================
    ! INIT_CO2 begins here!
    !=================================================================

    ! Return success
    RC          = GC_SUCCESS

    ! Exit if we have already intialized
    IF ( IS_INIT ) RETURN

    ! Copy values from Input_Opt
    LFOSSIL     = Input_Opt%LFOSSIL
    LCHEMCO2    = Input_Opt%LCHEMCO2
    LBIODIURNAL = Input_Opt%LBIODIURNAL
    LBIONETCLIM = Input_Opt%LBIONETCLIM
    LOCEAN      = Input_Opt%LOCEAN
    LSHIP       = Input_Opt%LSHIP
    LPLANE      = Input_Opt%LPLANE
    LFFBKGRD    = Input_Opt%LFFBKGRD
    LBIOSPHTAG  = Input_Opt%LBIOSPHTAG
    LFOSSILTAG  = Input_Opt%LFOSSILTAG
    LSHIPTAG    = Input_Opt%LSHIPTAG
    LPLANETAG   = Input_Opt%LPLANETAG

    ! Array for Fossil Fuel regions
    ALLOCATE( FOSSIL_REGION( State_Grid%NX, State_Grid%NY ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'FOSSIL_REGION' )
    FOSSIL_REGION = 0

    ! Array for Biospheric regions
    ALLOCATE( BIOSPH_REGION( State_Grid%NX, State_Grid%NY ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOSPH_REGION' )
    BIOSPH_REGION = 0

    ! Array for Ocean Regions
    ALLOCATE( OCEAN_REGION( State_Grid%NX, State_Grid%NY ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCEAN_REGION' )
    OCEAN_REGION = 0

    !=================================================================
    ! Set up regions
    !=================================================================

    ! Set up tagged regions for balanced biosphere & ocean
    IF ( Input_Opt%LBIOSPHTAG ) THEN
       CALL DEF_BIOSPH_CO2_REGIONS_F( State_Grid, BIOSPH_REGION )
       CALL DEF_OCEAN_CO2_REGIONS_F( State_Grid, OCEAN_REGION )
    ENDIF

    ! Set up tagged regions for fossil fuel
    IF ( Input_Opt%LFOSSILTAG ) THEN
       CALL DEF_FOSSIL_CO2_REGIONS_F( State_Grid, FOSSIL_REGION )
    ENDIF

    ! Reset IS_INIT flag
    IS_INIT = .TRUE.

  END SUBROUTINE INIT_CO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_co2
!
! !DESCRIPTION: Subroutine CLEANUP\_CO2 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_CO2
!
! !REVISION HISTORY:
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CLEANUP_CO2 begins here!
    !=================================================================
    IF ( ALLOCATED( FOSSIL_REGION ) ) DEALLOCATE( FOSSIL_REGION )
    IF ( ALLOCATED( BIOSPH_REGION ) ) DEALLOCATE( BIOSPH_REGION )
    IF ( ALLOCATED( OCEAN_REGION  ) ) DEALLOCATE( OCEAN_REGION  )

  END SUBROUTINE CLEANUP_CO2
!EOC
END MODULE CO2_MOD
