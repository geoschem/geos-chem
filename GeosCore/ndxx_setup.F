!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ndxx_setup.F
!
! !DESCRIPTION: Subroutine NDXX\_SETUP dynamically allocates memory for
!  certain diagnostic arrays that  are declared allocatable in "diag\_mod.f".
!\\
!\\
!  This allows us to reduce the amount of memory that needs to be declared
!  globally.  We only allocate memory for arrays if the corresponding
!  diagnostic is turned on.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NDXX_SETUP( am_I_Root,  Input_Opt, State_Chm,
     &                       State_Grid, RC )
!
! !USES:
!
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Grid_Mod,     ONLY : GrdState
      USE State_Met_Mod,      ONLY : MetState
#ifdef BPCH_DIAG
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD
      USE DIAG_OH_MOD,        ONLY : INIT_DIAG_OH
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR,   ERROR_STOP
      USE PLANEFLIGHT_MOD,    ONLY : SETUP_PLANEFLIGHT
#endif
#ifdef TOMAS
      USE TOMAS_MOD,          ONLY : IBINS, ICOMP, IDIAG   !(win, 7/9/09)
#endif

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)     :: am_I_Root   ! Are we on the root CPU?
      TYPE(GrdState), INTENT(IN)     :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT)  :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)    :: RC          ! Success or failure?
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE WILL BE A STUB UNLESS GEOS-Chem IS COMPILED    %%%
!  %%%  WITH THE BPCH_DIAG=y OPTION. (bmy, 10/4/19)                      %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! !REVISION HISTORY:
!  16 Jun 1998 - I. Bey, R. Yantosca - Initial version
!  See the Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef BPCH_DIAG
      INTEGER :: NMAX, AS, NEMISS, LMAX, nAdvect, NFAM
      LOGICAL :: IT_IS_A_FULLCHEM_SIM
      LOGICAL :: IT_IS_A_MERCURY_SIM
      LOGICAL :: IT_IS_A_TAGO3_SIM
      LOGICAL :: LDUST
      LOGICAL :: LCARB
      LOGICAL :: LSSALT
      LOGICAL :: LDRYD
      LOGICAL :: LGTMM

      !=================================================================
      ! NDXX_SETUP begins here!
      !=================================================================
      LD15 = 1
      LD16 = 1
      LD17 = 1
      LD18 = 1
      LD19 = 1
      LD21 = 1
      LD22 = 1
      LD24 = 1
      LD25 = 1
      LD26 = 1
      LD31 = 1
      LD37 = 1
      LD38 = 1
      LD39 = 1
      LD43 = 1
      LD45 = 1
      LD47 = 1
      LD52 = 1
      LD54 = 1
      LD66 = 1
      LD68 = 1

      !=================================================================
      ! ND01: Rn, Pb, Be emissions
      !=================================================================
      IF ( ND01 > 0 ) THEN
         LD01 = MIN( ND01, State_Grid%NZ )
         ALLOCATE( AD01( State_Grid%NX, State_Grid%NY, LD01, 2 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD01' )
      ENDIF

      !=================================================================
      ! ND02: Rn, Pb, Be decay
      !=================================================================
      IF ( ND02 > 0 ) THEN
         LD02 = MIN( ND02, State_Grid%NZ )
         ALLOCATE( AD02( State_Grid%NX, State_Grid%NY, LD02, nAdvect ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD02' )
      ENDIF

      !=================================================================
      ! ND04: CO2 source - see ??
      !
      ! ND05: Sulfate Prod/loss
      !=================================================================
      IF ( ND05 > 0 ) THEN

         ! Now size ND05 for all vertical levels (bmy, 3/9/18)
         LD05 = MIN( ND05, State_Grid%NZ )

         ALLOCATE( AD05( State_Grid%NX, State_Grid%NY, LD05, PD05 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD05' )
      ENDIF

      !=================================================================
      ! ND07: Carbonaceous aerosols emissions and chemical conversion
      !=================================================================
      IF ( ND07 > 0 .and. LCARB ) THEN
         LD07 = MIN( ND07, State_Grid%NZ )

         ALLOCATE( AD07_BC( State_Grid%NX, State_Grid%NY, LD07 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07_BC' )

         ALLOCATE( AD07_OC( State_Grid%NX, State_Grid%NY, LD07 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07_OC' )

         ! SOAupdate (hotp 6/15/09)
         ! semivolpoa2: increase dimension for POA (hotp 3/3/09)
         ! semivolpoa4: increase dimension for OPOA(hotp 3/27/09)
         ! semivolpoa4: add dimension for POG ox to OPOG (hotp 3/28/09)
         ! NAPSOA: increase AD07_HC dimension by one to 9 for production
         ! of SOA from IVOCs (hotp 7/22/09)
         ALLOCATE( AD07_HC( State_Grid%NX, State_Grid%NY, LD07, 9 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07_HC' )
      ENDIF

      !=================================================================
      ! ND12: Distribution of emissions in boundary layer [fraction]
      !       --> uses AD12 array (allocatable)
      !=================================================================
      LD12 = MIN( ND12, State_Grid%NZ )

      IF ( ND12 > 0 ) THEN
         ALLOCATE( AD12( State_Grid%NX, State_Grid%NY, LD12 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PL24H' )
      ENDIF

      !=================================================================
      ! ND14: Upward flux of from wet conv [kg/s]
      !       --> uses CONVFLUP array (allocatable)
      !=================================================================
      IF ( ND14 > 0 ) THEN
         LD14 = MIN( ND14, State_Grid%NZ )
         NMAX = nAdvect

         ALLOCATE( CONVFLUP( State_Grid%NX, State_Grid%NY,
     &                       State_Grid%NZ, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CONVFLUP' )
      ENDIF

      !=================================================================
      ! ND15: Mass change from BL-mixing [kg/s]
      !       --> uses TURBFLUP array (allocatable)
      !=================================================================
      IF ( ND15 > 0 ) THEN
         LD15 = MIN( ND15, State_Grid%NZ )
         NMAX = nAdvect

         ALLOCATE( TURBFLUP( State_Grid%NX, State_Grid%NY,
     &                       State_Grid%NZ, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TURBFLUP' )
      ENDIF

      !=================================================================
      ! ND16: Fraction of grid box experiencing large-scale and
      !       convective precipitation --> uses AD16 array (allocatable)
      !=================================================================
      IF ( ND16 > 0 ) THEN
         LD16 = MIN( ND16, State_Grid%NZ )

         ! Store LS Precip fraction only
         ALLOCATE( AD16( State_Grid%NX, State_Grid%NY, LD16, 1 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD16' )

         ! Counter array for AD16
         ALLOCATE( CT16( State_Grid%NX, State_Grid%NY, LD16, 1 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CT16' )
      ENDIF

      !=================================================================
      ! ND17: Fraction of tracer lost to rainout (in both large-scale
      !       and conv precipitation) --> uses AD17 array (allocatable)
      !=================================================================
      IF ( ND17 > 0 ) THEN

         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep

         ! Turn off ND17 if there are no wetdep species
         IF ( NMAX == 0 ) THEN
            ND17 = 0
            IF ( am_I_Root ) THEN
               print*, '### No wetdep species, so turn off ND17'
            ENDIF
         ENDIF

         ! Number of levels
         LD17 = MIN( ND17, State_Grid%NZ )

         ! Store both LS precip only
         ALLOCATE( AD17( State_Grid%NX, State_Grid%NY, LD17, NMAX, 1 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD17' )

         ! Counter array for AD17
         ALLOCATE( CT17( State_Grid%NX, State_Grid%NY, LD17, 1 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CT17' )
      ENDIF

      !=================================================================
      ! ND18: Fraction of tracer lost to washout (in both large-scale
      !       and convective precipitation) --> uses AD18 array (alloc.)
      !=================================================================
      IF ( ND18 > 0 ) THEN

         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep

         ! Turn off ND17 if there are no wetdep species
         IF ( NMAX == 0 ) THEN
            ND18 = 0
            IF ( am_I_Root ) THEN
               print*, '### No wetdep species, so turn off ND18'
            ENDIF
         ENDIF

         ! Number of levels
         LD18 = MIN( ND18, State_Grid%NZ )

         ! Store both LS precip fraction only
         ALLOCATE( AD18( State_Grid%NX, State_Grid%NY, LD18, NMAX, 1 ),
     &             STAT=AS)
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD18' )

         ! Counter array for AD17
         ALLOCATE( CT18( State_Grid%NX, State_Grid%NY, LD18, 1 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CT18' )
      ENDIF

      !=================================================================
      ! ND19: CH4 loss by OH
      !=================================================================
      IF ( ND19 > 0 ) THEN
         LD19 = MIN( ND19, State_Grid%NZ )

         ALLOCATE( AD19( State_Grid%NX, State_Grid%NY, LD19, PD19 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD19' )
      ENDIF

      !=================================================================
      ! ND20: Save O3 P-L losses to disk for single-tracer O3 run
      !       in the PL24H array.  Also turn on ND65, since the P-L
      !       rates are computed by ND65.
      !=================================================================
      IF ( ND20 > 0 ) THEN
         IF ( ND65 == 0 ) THEN
            CALL ERROR_STOP( 'ND65 must be turned on for ND20 output!',
     &                       'ndxx_setup.f'  )
         ENDIF
      ENDIF
#endif
#if defined( BPCH_DIAG )
      !=================================================================
      ! ND21: Optical depths and cloud fractions [unitless]
      !       --> uses AD21 array (allocatable)
      !=================================================================
      IF ( Input_Opt%ND21 > 0 ) THEN
         ! Now size diagnostic arrays w/ all vertical levels (bmy, 3/9/18)
         Input_Opt%LD21 = MIN( Input_Opt%ND21, State_Grid%MaxChemLev )
         LD21 = Input_Opt%LD21

         ! For regular
         ALLOCATE( AD21( State_Grid%NX, State_Grid%NY, Input_Opt%LD21,
     &                   PD21 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD21' )
      ENDIF
#endif
#if defined( BPCH_DIAG )

      !=================================================================
      ! ND22: J-value diagnostics [s^-1]
      !       --> uses AD22 array (allocatable)
      !=================================================================
      IF ( ND22 > 0 ) THEN

         ! Now size J-value diagnostic for all vertical levels (bmy, 3/9/18)
         LD22 = MIN( ND22, State_Grid%NZ  )

         ! Accumulating diagnostic array
         ALLOCATE( AD22( State_Grid%NX, State_Grid%NY, LD22, PD22 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD22' )

         ! Number of times where LT is between HR1_JV and HR2_JV
         ALLOCATE( CTJV( State_Grid%NX, State_Grid%NY ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTJV' )
      ENDIF

      ! Locations where LT is between HR1_JV and HR2_JV
      ! NOTE: LTJV is currently used for netCDF diagnostics in flexchem_mod.F90
      ! regardless of the ND22 setting, so we always need to allocate it here
      ! (mps, 1/5/18)
      ALLOCATE( LTJV( State_Grid%NX, State_Grid%NY ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTJV' )

      !=================================================================
      ! ND27: Flux of Ox across the annual mean tropopause [kg/s]
      !       ND27 will also turn on ND24, ND25, ND26 diagnostics
      !=================================================================
      IF ( ND27 > 0 ) THEN
         ND24 = State_Grid%NZ
         ND25 = State_Grid%NZ
         ND26 = State_Grid%NZ
      ENDIF

      !=================================================================
      ! ND24: Eastward mass flux from transport [kg/s]
      !       --> uses MASSFLEW array (allocatable)
      !=================================================================
      IF ( ND24 > 0 ) THEN
         LD24 = MIN( ND24, State_Grid%NZ )
         NMAX = nAdvect
         ALLOCATE( MASSFLEW( State_Grid%NX, State_Grid%NY,
     &                       State_Grid%NZ, NMAX ), STAT=AS)
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLEW' )
      ELSE
         ALLOCATE( MASSFLEW( 1, 1, 1, 1 ), STAT=AS)
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLEW' )
      ENDIF

      !=================================================================
      ! ND25: Northward mass flux from transport [kg/s]
      !       --> uses MASSFLNS array (allocatable)
      !=================================================================
      IF ( ND25 > 0 ) THEN
         LD25 = MIN( ND25, State_Grid%NZ )
         NMAX = nAdvect
         ALLOCATE( MASSFLNS( State_Grid%NX, State_Grid%NY,
     &                       State_Grid%NZ, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLNS' )
      ELSE
         ALLOCATE( MASSFLNS( 1, 1, 1, 1 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLNS' )
      ENDIF

      !=================================================================
      ! ND26: Vertical mass flux from transport [kg/s]
      !       --> uses MASSFLUP array (allocatable)
      !=================================================================
      IF ( ND26 > 0 ) THEN
         LD26 = MIN( ND26, State_Grid%NZ )
         NMAX = nAdvect
         ALLOCATE( MASSFLUP( State_Grid%NX, State_Grid%NY,
     &                       State_Grid%NZ, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLUP' )
      ELSE
         ALLOCATE( MASSFLUP( 1, 1, 1, 1 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLUP' )
      ENDIF

      !=================================================================
      ! ND30: Land/water/ice flags
      !=================================================================
      IF ( ND30 > 0 ) THEN
         ALLOCATE( AD30( State_Grid%NX, State_Grid%NY ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD30' )
      ENDIF

      !=================================================================
      ! ND31: 3-D pressure edges [hPa] --> Uses AD31 array (allocatable)
      !=================================================================
      IF ( ND31 > 0 ) THEN
         LD31 = MIN( ND31, State_Grid%NZ+1 )
         ALLOCATE( AD31( State_Grid%NX, State_Grid%NY, LD31 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD31' )
      ENDIF

      !=================================================================
      ! ND33: Column sum of tracer [kg]
      !       --> uses AD33 array (allocatable)!
      !=================================================================
      IF ( ND33 > 0 ) THEN
         ALLOCATE( AD33( State_Grid%NX, State_Grid%NY, nAdvect ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD33' )
      ENDIF

      !=================================================================
      ! ND35: Tracer at 500 mb [v/v] (this is ~ level 9 for GEOS-CHEM)
      !       --> uses AD35 array (allocatable)
      !=================================================================
      IF ( ND35 > 0 ) THEN
         ALLOCATE( AD35( State_Grid%NX, State_Grid%NY, nAdvect ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD35' )
      ENDIF

      !=================================================================
      ! ND37: Fraction of tracer scavenged in cloud updrafts
      !       --> Uses AD37 array (allocatable)
      !=================================================================
      IF ( ND37 > 0 ) THEN

         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep

         ! Turn off ND37 if there are no wetdep species
         IF ( NMAX == 0 ) THEN
            ND37 = 0
            IF ( am_I_Root ) THEN
               print*, '### No wetdep species, so turn off ND37'
            ENDIF
         ENDIF

         ! Number of levels
         LD37 = MIN( ND37, State_Grid%NZ )

         ! Allocate array accordingly
         IF ( NMAX > 0 ) THEN
            ALLOCATE( AD37( State_Grid%NX, State_Grid%NY, LD37, NMAX ),
     &                STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD37' )
         ENDIF
      ENDIF

      !=================================================================
      ! ND38: Rainout of tracer (nfcldmx.f)
      !       --> uses AD38 array (allocatable)
      !=================================================================
      ! Allocate AD38 for GTMM. (ccc, 10/29/09)
      IF ( ND38 > 0 .OR. LGTMM ) THEN

         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep
         LD38 = MIN( ND38, State_Grid%NZ )

         ! If GTMM on and LD38 off, ND38 should be State_Grid%NZ,(ccc, 10/29/09)
         IF ( ND38 == 0 ) LD38 = State_Grid%NZ

         ! Allocate AD38 array accordingly
         IF ( NMAX > 0 ) THEN
            ALLOCATE( AD38( State_Grid%NX, State_Grid%NY, LD38, NMAX ),
     &                STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD38' )
         ENDIF
      ENDIF

      !=================================================================
      ! ND39: Rainout of tracer (wetdep.f)
      !       --> uses AD39 array (allocatable)
      !=================================================================
      ! Allocate AD39 for LGTMM. (ccc, 10/29/09)
      IF ( ND39 > 0 .OR. LGTMM ) THEN

         LD39 = MIN( ND39, State_Grid%NZ )
         ! If GTMM on and ND38 off, ND38 should be State_Grid%NZ,(ccc, 10/29/09)
         IF ( ND39 == 0 ) LD39 = State_Grid%NZ

         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep

         ! Allocate AD39 array accordingly
         IF ( NMAX > 0 ) THEN
            ALLOCATE( AD39( State_Grid%NX, State_Grid%NY, LD39, NMAX ),
     &                STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD39' )
         ENDIF
      ENDIF

      !=================================================================
      ! ND42: SOA concentration
      !
      ! ND43: Chemical diagnostics: OH [molec/cm3/s] and HO2 [v/v]
      !       --> uses AD43 array (allocatable)
      !=================================================================
      IF ( ND43 > 0 ) THEN

         ! Now size diagnostic arrays with all vertical levels (bmy, 3/9/18)
         LD43 = MIN( ND43, State_Grid%NZ  )

         ! Accumulating diagnostic array
         ALLOCATE( AD43( State_Grid%NX, State_Grid%NY, LD43, PD43 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD43' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTOH( State_Grid%NX, State_Grid%NY, LD43 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTOH' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTHO2( State_Grid%NX, State_Grid%NY, LD43 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTHO2' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTO1D( State_Grid%NX, State_Grid%NY, LD43 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTO1D' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTO3P( State_Grid%NX, State_Grid%NY, LD43 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTO3P' )
      ENDIF

      ! The following are also used for netcdf diagnostics and therefore
      ! must come out of the BPCH_DIAG ifdef and not be predicated on
      ! an NDXX flag. Additional work should be done to only allocated
      ! these arrays if either the bpch or netcdf diagnostic is on
      ! (ewl, 1/5/2018)

      ! Locations where LT is between HR1_OH and HR2_OH
      ALLOCATE( LTOH( State_Grid%NX, State_Grid%NY ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTOH' )

      ! Locations where LT is between HR1_OH and HR2_OH
      ALLOCATE( LTHO2( State_Grid%NX, State_Grid%NY ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTHO2' )


      ! Assume successful return
      RC = GC_SUCCESS

      nAdvect = State_Chm%nAdvect
      LDUST   = Input_Opt%LDUST
      LCARB   = Input_Opt%LCARB

#ifdef TOMAS
      !=================================================================
      ! ND44: Drydep fluxes [s-1] and drydep velocities [cm/s]
      !       --> uses AD44 arrays (allocatable)
      !=================================================================
      IF ( .not. LDRYD ) ND44 = 0

      IF ( ND44 > 0 ) THEN
         ! add space in diag array for TOMAS aerosol mass (win, 7/14/09)
         ! Now use State_Chm%nDryDep for # dry depositing species
         ! (ewl, 10/14/15)
         IF ( Ind_('NK1') > 1 ) THEN
            NMAX = State_Chm%nDryDep + ( ICOMP - IDIAG )* IBINS
         ELSE
            NMAX = State_Chm%nDryDep
         ENDIF
         ! Allocate AD44 array
         ALLOCATE( AD44( State_Grid%NX, State_Grid%NY, NMAX, 2 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD44' )


      !=================================================================
      ! ND59: Size-resolved primary aerosol emissions      !(win, 7/9/09)
      !         Emissions to number, sulfate, sea-salt, carb, dust
      !      ----> save 3-D (I,J,1) or up to (I,J,2)
      !=================================================================
      IF ( ND59 > 0 ) THEN
         LD59 = MIN( ND59, State_Grid%NZ )

         ! Number emission
         ALLOCATE( AD59_NUMB( State_Grid%NX, State_Grid%NY, 2, IBINS ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_NUMB' )

         ALLOCATE( AD59_SULF( State_Grid%NX, State_Grid%NY, 2, IBINS ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_SULF' )

         ALLOCATE( AD59_SALT( State_Grid%NX, State_Grid%NY, 2, IBINS ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_SALT' )

         ALLOCATE( AD59_ECIL( State_Grid%NX, State_Grid%NY, 2, IBINS ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_ECIL' )

         ALLOCATE( AD59_ECOB( State_Grid%NX, State_Grid%NY, 2, IBINS ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_ECOB' )

         ALLOCATE( AD59_OCIL( State_Grid%NX, State_Grid%NY, 2, IBINS ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_OCIL' )

         ALLOCATE( AD59_OCOB( State_Grid%NX, State_Grid%NY, 2, IBINS ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_OCOB' )

         ALLOCATE( AD59_DUST( State_Grid%NX, State_Grid%NY, 2, IBINS ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_DUST' )
      ENDIF

      !=================================================================
      ! ND60: TOMAS microphysical process rates (condensation,
      !           coagulation, nucleation, aqueous oxidation,
      !           error-fudging)
      !       ---> save 2-D (J,L) for 30-bin of each aerosol species
      !=================================================================
      IF ( ND60 > 0 ) THEN
         LD60 = MIN( ND60, State_Grid%NZ )

         ! Now the array dimension is IBINS*(ICOMP-IDIAG+1) because
         ! we need it for all prognostic mass species + 1 for number
         ! IDIAG = # of diagnostic species.  (win, 9/27/08)
         !Condensation rate
         ALLOCATE(AD60_COND(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)),
     &            STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_COND' )

         !Coagulation rate
         ALLOCATE(AD60_COAG(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)),
     &            STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_COAG' )

         !Nucleation rate
         ALLOCATE(AD60_NUCL(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)),
     &            STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_NUCL' )

         !Aqueous oxidation rate
         ALLOCATE(AD60_AQOX(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)),
     &            STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_AQOX' )

         !Accumulated error-fudging
         ALLOCATE(AD60_ERROR(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1))
     &            ,STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_ERROR' )

         !SOA Condensation rate
         ALLOCATE(AD60_SOA(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)),
     &            STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_SOA' )
      ENDIF

      !=================================================================
      ! ND61: 3-D TOMAS process rate diagnostic
      !       --> Uses AD61 array (allocatable)
      !  NOTE: ND61 is used for 10-nm particle formation
      !        rate and cluster-size nucleation rate.  So the array
      !        is declared for (NX,NY,LD61,2) (win, 10/6/08)
      !=================================================================
      IF ( ND61 > 0 ) THEN
         LD61 = MIN( ND61, State_Grid%NZ )

         ALLOCATE( AD61( State_Grid%NX, State_Grid%NY, LD61, PD61 ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD61' )
         ALLOCATE( AD61_INST( State_Grid%NX, State_Grid%NY, LD61, PD61),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD61_INST' )
      ENDIF
#endif

#ifdef RRTMG
      !=================================================================
      ! ND72: Radiative output diagnostics (TOASW, SRFSW, TOALW, SRFLW,
      !       GTOASW, GSRFSW, GTOALW, GSRFLW, ATOASW, ASRFSW, ATOALW,
      !       ASRFLW) for both clear and all sky (prefix CLR and ALL)
      !       --> uses AD72 array (allocatable)
      !=================================================================
      IF ( ND72 > 0 ) THEN
         ALLOCATE( AD72( State_Grid%NX, State_Grid%NY, PD72 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD72' )
      ENDIF

      Input_Opt%LD59  = LD59
      Input_Opt%LD60  = LD60
      Input_Opt%LD72  = LD72
#endif

#endif
      END SUBROUTINE NDXX_SETUP
!EOC
