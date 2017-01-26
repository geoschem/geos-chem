!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: co2_mod.F
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
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now references biomass_mod.f (bmy, 9/27/06)
!  (3 ) Tagged CO2 capability developed (dbj)
!  (4 ) Implemented monthly and annual fossil fuel inventories 
!        (R.Nassar 2009-03-10)
!  (5 ) Implemented CO2 emissions from shipping and aviation (R.Nassar 2010)
!  (6 ) Implemented monthly CO2 chemical production and surface correction 
!        (R.Nassar 2010)   
!  25 Feb 2011 - R. Nassar  - Now read updated CDIAC CO2 emissions data
!  07 Sep 2011 - P. Kasibhatla - Modified to include GFED3
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  20 Jun 2014 - R. Yantosca - Remove obsolete emissions code; we now use HEMCO
!  07 Nov 2014 - M. Yannetti - Added PRECISION_MOD
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

      ! XNUMOL_CO2   - molecules CO2 / kg CO2 
      REAL(fp),  PARAMETER   :: XNUMOL_CO2 = AVO / FMOL_CO2

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
      SUBROUTINE EMISSCO2( am_I_Root, Input_Opt,
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE DIAG04_MOD,         ONLY : AD04, ND04
      USE DIAG04_MOD,         ONLY : AD04_plane,    AD04_chem
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
      USE HCO_STATE_MOD,      ONLY : HCO_GetHcoID
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!

      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  The initial condition for CO2 has to be at least 50 ppm or higher or else
!  the balanced biosphere fluxes will make STT negative. (pns, bmy, 8/16/05)
! 
! 
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam - Initial version
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) We now get CO2 biomass emissions from biomass_mod.f.  This allows us 
!        to use either GFED2 or default Duncan et al biomass emissions. 
!        (bmy, 9/27/06)
!  (3 ) Tagged tracer capability added. This requires the editable region 
!        files Regions_land.dat and Regions_ocean.dat in the run directory
!        (rnassar,dbj, 2009)
!  (4 ) New tracers for emissions from international and domestic shipping,
!        international and domestic aviation, and the chemical CO2 source
!        from the oxidation of CO, CH4, and other organics (rnassar,dbj, 2009)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  27 Aug 2014 - C. Keller   - Now interface with HEMCO.
!  13 Apr 2015 - R. Nassar   - Some bug fixes for HEMCO implementation, 
!                              e.g. unit fixes for aviation and chem source
!  11 Sep 2015 - E. Lundgren - Tracer units are now kg/kg dry air (prev kg)
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                    :: I, J, L, N, NA
      INTEGER                    :: nAdvect
      LOGICAL                    :: LFOSSIL
      LOGICAL                    :: LCHEMCO2
      LOGICAL                    :: LBIODIURNAL
      LOGICAL                    :: LBIONETCLIM
      LOGICAL                    :: LOCEAN
      LOGICAL                    :: LSHIP
      LOGICAL                    :: LPLANE
      LOGICAL                    :: LFFBKGRD
      LOGICAL                    :: LBIOSPHTAG
      LOGICAL                    :: LFOSSILTAG
      LOGICAL                    :: LSHIPTAG
      LOGICAL                    :: LPLANETAG
      REAL(fp)                   :: A_CM2, DTSRCE, E_CO2

      ! SAVEd scalars
      LOGICAL,         SAVE      :: FIRST      = .TRUE.
      LOGICAL,         SAVE      :: BIOFUEL_ON = .TRUE.
      LOGICAL,         SAVE      :: BIOMASS_ON = .TRUE.
      INTEGER,         SAVE      :: IDff       = -1
      INTEGER,         SAVE      :: IDoc       = -1
      INTEGER,         SAVE      :: IDbal      = -1
      INTEGER,         SAVE      :: IDbb       = -1
      INTEGER,         SAVE      :: IDbf       = -1
      INTEGER,         SAVE      :: IDnte      = -1
      INTEGER,         SAVE      :: IDse       = -1
      INTEGER,         SAVE      :: IDav       = -1
      INTEGER,         SAVE      :: IDch       = -1
      INTEGER,         SAVE      :: IDcorr     = -1

      ! Strings
      CHARACTER(LEN=255)         :: LOC = 'EMISSCO2 (co2_mod.F)'

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

      ! Return success
      RC          = GC_SUCCESS

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
  
      ! Initialize pointers
      CO2_COPROD  => NULL()
      Spc         => NULL()

      ! Import emissions from HEMCO (through HEMCO state)
      IF ( .NOT. ASSOCIATED(HcoState) ) THEN
         CALL ERROR_STOP ( 'HcoState not defined!', LOC ) 
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
         ! Query HEMCO for all of the species IDs
         !--------------------------------------------------------------
         IDff   = HCO_GetHcoID( 'CO2ff',   HcoState )
         IDoc   = HCO_GetHcoID( 'CO2oc',   HcoState )
         IDbal  = HCO_GetHcoID( 'CO2bal',  HcoState )
         IDbb   = HCO_GetHcoID( 'CO2bb',   HcoState )
         IDbf   = HCO_GetHcoID( 'CO2bf',   HcoState )
         IDnte  = HCO_GetHcoID( 'CO2nte',  HcoState )
         IDse   = HCO_GetHcoID( 'CO2se',   HcoState )
         IDav   = HCO_GetHcoID( 'CO2av',   HcoState )
         IDch   = HCO_GetHcoID( 'CO2ch',   HcoState )
         IDcorr = HCO_GetHcoID( 'CO2corr', HcoState )

         !--------------------------------------------------------------
         ! Add error checks.  Check to see that the species index
         ! and the corresponding field in the HEMCO state are defined.
         ! Move this outside the parallel loop to avoid problems when
         ! exiting while w/in a parallel region. (bmy, 4/17/15)
         !--------------------------------------------------------------
         IF ( LFOSSIL ) THEN
            IF ( IDff < 1 ) THEN
               CALL ERROR_STOP( 'CO2ff not defined', LOC )
            ENDIF
            IF ( .not. ASSOCIATED(HcoState%Spc(IDff)%Emis%Val) ) THEN
               CALL ERROR_STOP('Fossil emissions not defined', LOC )
            ENDIF
         ENDIF

         IF ( LOCEAN ) THEN
            IF ( IDoc < 1 ) THEN
               CALL ERROR_STOP( 'CO2oc not defined', LOC )
            ENDIF
            IF ( .not. ASSOCIATED(HcoState%Spc(IDoc)%Emis%Val) ) THEN
               CALL ERROR_STOP( 'Oceanic emissions not defined', LOC )
            ENDIF
         ENDIF

         IF ( LBIODIURNAL ) THEN
            IF ( IDbal < 1 ) THEN
               CALL ERROR_STOP( 'CO2bal not defined' ,LOC )
            ENDIF
            IF ( .not. ASSOCIATED(HcoState%Spc(IDbal)%Emis%Val) ) THEN
               CALL ERROR_STOP( 'Balanced biosphere not defined', LOC )
            ENDIF
         ENDIF

         IF ( LBIONETCLIM ) THEN
            IF ( IDnte < 1 ) THEN
               CALL ERROR_STOP(' CO2nte not defined', LOC )
            ENDIF
            IF ( .not. ASSOCIATED(HcoState%Spc(IDnte)%Emis%Val) ) THEN
               CALL ERROR_STOP( 'Net terr. exchange not defined', LOC )
            ENDIF
         ENDIF

         IF ( LSHIP ) THEN
            IF ( IDse < 1 ) THEN
               CALL ERROR_STOP( 'CO2se not defined', LOC )
            ENDIF
            IF ( .not. ASSOCIATED(HcoState%Spc(IDse)%Emis%Val) ) THEN
               CALL ERROR_STOP('Ship emissions not defined',LOC)
            ENDIF
         ENDIF

         IF ( LPLANE ) THEN
            IF ( IDav < 1 ) THEN
               CALL ERROR_STOP( 'CO2av not defined', LOC )
            ENDIF
            IF ( .not. ASSOCIATED(HcoState%Spc(IDav)%Emis%Val) ) THEN
               CALL ERROR_STOP( 'Aviation emissions not defined', LOC )
            ENDIF
         ENDIF

         IF ( LCHEMCO2 ) THEN
            IF ( IDch < 1 ) THEN
               CALL ERROR_STOP( 'CO2ch not defined', LOC )
            ENDIF
            IF ( IDcorr < 1 ) THEN
               CALL ERROR_STOP( 'CO2corr not defined', LOC )
            ENDIF
            IF ( .not. ASSOCIATED(HcoState%Spc(IDcorr)%Emis%Val) ) THEN
               CALL ERROR_STOP( 'Surface correction not defined', LOC )
            ENDIF
         ENDIF

         !--------------------------------------------------------------
         ! Also test if the biofuel & biomass species are on
         ! (i.e. if the pointers in the HEMCO state are associated).
         ! We can do this just once, outside of the parallel loop
         !--------------------------------------------------------------
         BIOFUEL_ON = ASSOCIATED( HcoState%Spc(IDbf)%Emis%Val )
         BIOMASS_ON = ASSOCIATED( HcoState%Spc(IDbb)%Emis%Val )

         !--------------------------------------------------------------
         ! Error check: For now, the emission grid must be
         ! on the simulation grid.
         !--------------------------------------------------------------
         IF ( HcoState%NX /= IIPAR .OR. 
     &        Hcostate%NY /= JJPAR .OR.
     &        Hcostate%NZ /= LLPAR     ) THEN
            CALL ERROR_STOP( 'HEMCO grid not same as sim. grid!', LOC )
         ENDIF 

         ! Set first-time flag to false
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read in monthly, daily or variable emissions fields
      !=================================================================      

      ! Point to the chemical species ([kg/kg] dry air)
      Spc  => State_Chm%Species

      ! Get a pointer to the CO2 production from CO oxidation
      ! This is now listed in the NON-EMISSIONS DATA section of
      ! the HEMCO configuration file. (bmy, 4/17/15)
      CALL HCO_GetPtr( am_I_Root,   HcoState, 
     &                'CO2_COPROD', CO2_COPROD, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'CO2 production not defined', LOC )
      ENDIF

      !=================================================================
      ! Process emissions and save diagnostics
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J, I, A_CM2, E_CO2, N )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

!         ! Grid box surface area [cm2]
!         A_CM2 = HcoState%Grid%AREA_M2%Val( I, J ) * CM2PERM2

         !-------------------------------------------
         ! #1: Total CO2
         ! #2: CO2 from fossil fuel emissions 
         !-------------------------------------------
         IF ( LFOSSIL ) THEN

            ! Get emissions [kg/m2/s] and convert to [molec/cm2/s]
            E_CO2 = HcoState%Spc(IDff)%Emis%Val(I,J,1) 
     &            / CM2PERM2 
     &            * XNUMOL_CO2

            ! ND04 diag: Fossil Fuel CO2 [molec/cm2/s] 
            IF ( ND04 > 0 ) THEN
               AD04(I,J,1) = AD04(I,J,1) + E_CO2
            ENDIF

! old
!            ! Convert from [molec/cm2/s] to [kg]
!            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
! new
            ! Get emissions [kg/m2/s] and convert to [kg/kg dry air]
            ! using method in ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod
            ! (ewl, 9/11/15)
            E_CO2  =  HcoState%Spc(IDff)%Emis%Val(I,J,1)
     &                * DTSRCE / ( State_Met%DELP(I,J,1) 
     &                * G0_100 * ( 1.0e+0_fp
     &                - State_Met%SPHU(I,J,1) * 1.0e-3_fp ) )        
! end new (ewl)

            ! Add to Species #1: Total CO2 [kg/kg]
            Spc(I,J,1,1)   = Spc(I,J,1,1) + E_CO2

	     ! Add to Species #2: Fossil CO2 [kg/kg]
            IF ( nAdvect > 1 ) THEN
               Spc(I,J,1,2)   = Spc(I,J,1,2) + E_CO2

               ! Split Fossil Fuel CO2 into geographic regions
               IF ( LFOSSILTAG ) THEN
                  N            = FOSSIL_REGION(I,J)
                  Spc(I,J,1,N) = Spc(I,J,1,N) + E_CO2
               ENDIF
            ENDIF

	     ! Add to Species #12: Background with Fossil CO2 [kg/kg]
            IF ( LFFBKGRD ) THEN
               Spc(I,J,1,12)   = Spc(I,J,1,12) + E_CO2
            ENDIF

	   ! Note: To define the background as including fossil fuels here 
           ! (instead of during data processing) for tagged-CO2 runs which 
           ! will be used for estimating natural bio/ocean fluxes etc., and 
           ! accepting the FF inventory
           ! use Spc(I,J,1,12) = Spc(I,J,1,12) + E_CO2.
	   ! Shipping and aviation emissions can be included in a similar way.
	   
         ENDIF

         !-------------------------------------------
         ! #3: CO2 from ocean exchange
         !-------------------------------------------
         IF ( LOCEAN ) THEN

            ! Get emissions [kg/m2/s] and convert to [molec/cm2/s]
            E_CO2 = HcoState%Spc(IDoc)%Emis%Val(I,J,1) 
     &            / CM2PERM2 
     &            * XNUMOL_CO2

            ! ND04 diag: Ocean CO2 [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,2) = AD04(I,J,2) + E_CO2
            ENDIF

! old
!            ! Convert from [molec/cm2/s] to [kg]
!            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
! new
            ! Get emissions [kg/m2/s] and convert to [kg/kg dry air]
            ! using method in ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod
            ! (ewl, 9/11/15)
            E_CO2  =  HcoState%Spc(IDoc)%Emis%Val(I,J,1)
     &                * DTSRCE / ( State_Met%DELP(I,J,1) 
     &                * G0_100 * ( 1.0e+0_fp
     &                - State_Met%SPHU(I,J,1) * 1.0e-3_fp ) )        
! end new (ewl)

            ! Add to Species #1: Total CO2 [kg/kg]
            Spc(I,J,1,1)   = Spc(I,J,1,1) + E_CO2

            ! Add to Species #3: Ocean CO2 [kg/kg]
            IF ( nAdvect > 2 ) THEN
               Spc(I,J,1,3)   = Spc(I,J,1,3) + E_CO2   

               ! Split ocean CO2 into geographic regions
               IF ( LBIOSPHTAG ) THEN
                  N            = OCEAN_REGION(I,J)
                  Spc(I,J,1,N) = Spc(I,J,1,N) + E_CO2
               ENDIF
            ENDIF
         ENDIF

         !-------------------------------------------
         ! #4: CO2 from balanced biosphere emissions
         !-------------------------------------------
         IF ( LBIODIURNAL ) THEN

            ! Get emissions [kg/m2/s] and convert to [molec/cm2/s]
            E_CO2 = HcoState%Spc(IDbal)%Emis%Val(I,J,1) 
     &            / CM2PERM2 
     &            * XNUMOL_CO2

            ! ND04 diag: Bal Bio CO2 [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,3) = AD04(I,J,3) + E_CO2
            ENDIF

! old
!           ! Convert from [molec/cm2/s] to [kg CO2]
!            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 
! new
            ! Get emissions [kg/m2/s] and convert to [kg/kg dry air]
            ! using method in ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod
            ! (ewl, 9/11/15)
            E_CO2  =  HcoState%Spc(IDbal)%Emis%Val(I,J,1)
     &                * DTSRCE / ( State_Met%DELP(I,J,1) 
     &                * G0_100 * ( 1.0e+0_fp
     &                - State_Met%SPHU(I,J,1) * 1.0e-3_fp ) )        
! end new (ewl)
 
            ! Add to Species #1 -- total CO2 [kg/kg]
            Spc(I,J,1,1)   = Spc(I,J,1,1) + E_CO2

            ! Add to Species #4 -- Bal Bio CO2 [kg/kg]
            IF ( nAdvect > 3 ) THEN
               Spc(I,J,1,4)   = Spc(I,J,1,4) + E_CO2
                
               	! Split biospheric CO2 exchange into geographic regions
               IF ( LBIOSPHTAG ) THEN
                  N            = BIOSPH_REGION(I,J)
                  Spc(I,J,1,N) = Spc(I,J,1,N) + E_CO2
               ENDIF
            ENDIF	
         ENDIF

         !-------------------------------------------
         ! #5: CO2 from biomass burning emissions
         !-------------------------------------------
         IF ( IDbb > 0 ) THEN

            ! Emission array will only be defined if biomass burning
            ! is enabled in HEMCO configuration file.
            IF ( BIOMASS_ON ) THEN

               ! Get emissions [kg/m2/s] and convert to [molec/cm2/s]
               E_CO2 = HcoState%Spc(IDbb)%Emis%Val(I,J,1) 
     &               / CM2PERM2 
     &               * XNUMOL_CO2
   
               ! ND04 diag: Biomass burning CO2 [molec/cm2/s]
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,4) = AD04(I,J,4) + E_CO2
               ENDIF

! old   
!               ! Convert from [molec/cm2/s] to [kg]
!               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 
! new
               ! Get emissions [kg/m2/s] and convert to [kg/kg dry air]
               ! using method in ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod
               ! (ewl, 9/11/15)
               E_CO2  =  HcoState%Spc(IDbb)%Emis%Val(I,J,1)
     &                   * DTSRCE / ( State_Met%DELP(I,J,1) 
     &                   * G0_100 * ( 1.0e+0_fp
     &                   - State_Met%SPHU(I,J,1) * 1.0e-3_fp ) )        
! end new (ewl)   
               ! Add to Species #1: Total CO2 [kg/kg]
               Spc(I,J,1,1)   = Spc(I,J,1,1) + E_CO2
               
               ! Add to Species #5: Biomass burning CO2 [kg/kg]
               IF ( nAdvect > 4 ) THEN
                  Spc(I,J,1,5)   = Spc(I,J,1,5) + E_CO2
   
                  ! Split Bioburn CO2 into geographic regions
                  IF ( LBIOSPHTAG ) THEN
                     N            = BIOSPH_REGION(I,J)
                     Spc(I,J,1,N) = Spc(I,J,1,N) + E_CO2
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         !-------------------------------------------
         ! #6: CO2 from biofuel emissions
         !-------------------------------------------
         IF ( IDbf > 0 ) THEN

            ! Emission array will only be defined if biofuel emissions
            ! are specified in HEMCO configuration file. 
            IF ( BIOFUEL_ON ) THEN
   
               ! Get emissions [kg/m2/s] and convert to [molec/cm2/s]
               E_CO2 = HcoState%Spc(IDbf)%Emis%Val(I,J,1) 
     &               / CM2PERM2 
     &               * XNUMOL_CO2
   
               ! ND04 diag: Biofuel CO2 [molec/cm2/s] 
               IF ( ND04 > 0 ) THEN
                  AD04(I,J,5) = AD04(I,J,5) + E_CO2
               ENDIF

! old   
!               ! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
!               E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
! new
               ! Get emissions [kg/m2/s] and convert to [kg/kg dry air]
               ! using method in ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod
               ! (ewl, 9/11/15)
               E_CO2  =  HcoState%Spc(IDbf)%Emis%Val(I,J,1)
     &                   * DTSRCE / ( State_Met%DELP(I,J,1) 
     &                   * G0_100 * ( 1.0e+0_fp
     &                   - State_Met%SPHU(I,J,1) * 1.0e-3_fp ) )        
! end new (ewl)
   
               ! Add to Species #1: Total CO2 [kg/kg]
               Spc(I,J,1,1)   = Spc(I,J,1,1) + E_CO2
   
               ! Add to Species #6: Biofuel CO2 [kg/kg]
               IF (nAdvect > 5) THEN
                  Spc(I,J,1,6)   = Spc(I,J,1,6) + E_CO2
   
                  ! Split BF CO2 into geographic regions
                  IF ( LBIOSPHTAG ) THEN
                     N            = BIOSPH_REGION(I,J)
                     Spc(I,J,1,N) = Spc(I,J,1,N) + E_CO2
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
   
         !-------------------------------------------
         ! #7: CO2 from net terrestrial exchange
         !-------------------------------------------
         IF ( LBIONETCLIM ) THEN

            ! Get emissions [kg/m2/s] and convert to [molec/cm2/s]
            E_CO2 = HcoState%Spc(IDnte)%Emis%Val(I,J,1) 
     &            / CM2PERM2 
     &            * XNUMOL_CO2

            ! ND04 diag: net terrestrial exchange [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,6) = AD04(I,J,6) + E_CO2
            ENDIF

! old
!            ! Convert from [molec/cm2/s] to [kg]
!            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
! new
            ! Get emissions [kg/m2/s] and convert to [kg/kg dry air]
            ! using method in ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod
            ! (ewl, 9/11/15)
            E_CO2  =  HcoState%Spc(IDnte)%Emis%Val(I,J,1)
     &                * DTSRCE / ( State_Met%DELP(I,J,1) 
     &                * G0_100 * ( 1.0e+0_fp
     &                - State_Met%SPHU(I,J,1) * 1.0e-3_fp ) )        
! end new (ewl)

            ! Add to Species #1: Total CO2 [kg CO2]
            Spc(I,J,1,1)   = Spc(I,J,1,1) + E_CO2

            ! Add to Species #7: Net Terr exchange CO2 [kg]
            IF ( nAdvect > 6 ) THEN
               Spc(I,J,1,7)   = Spc(I,J,1,7) + E_CO2

               ! Split Net Terr Exch CO2 into geographic regions
               IF ( LBIOSPHTAG ) THEN
                  N            = BIOSPH_REGION(I,J)
                  Spc(I,J,1,N) = Spc(I,J,1,N) + E_CO2
               ENDIF
            ENDIF	
         ENDIF

         !-------------------------------------------
         ! #8: CO2 from ship emissions
         !-------------------------------------------
         IF ( LSHIP ) THEN

            ! Get emissions [kg/m2/s] and convert to [molec/cm2/s]
            E_CO2 = HcoState%Spc(IDse)%Emis%Val(I,J,1) 
     &            / CM2PERM2 
     &            * XNUMOL_CO2

            ! ND04 diag: Ship CO2 [molec/cm2/s] 
            IF ( ND04 > 0 ) THEN
               AD04(I,J,7) = AD04(I,J,7) + E_CO2
            ENDIF

! old
!            ! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
!            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
! new
            ! Get emissions [kg/m2/s] and convert to [kg/kg dry air]
            ! using method in ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod
            ! (ewl, 9/11/15)
            E_CO2  =  HcoState%Spc(IDse)%Emis%Val(I,J,1)
     &                * DTSRCE / ( State_Met%DELP(I,J,1) 
     &                * G0_100 * ( 1.0e+0_fp
     &                - State_Met%SPHU(I,J,1) * 1.0e-3_fp ) )        
! end new (ewl)

            ! Add to Species #1: Total CO2 [kg/kg]
            Spc(I,J,1,1)   = Spc(I,J,1,1) + E_CO2

            ! Add to Species #8: Ship CO2 [kg/kg]
            IF (nAdvect > 7) THEN
               Spc(I,J,1,8)   = Spc(I,J,1,8) + E_CO2
 
              	! Tagged species for global ship emissions
               IF ( LSHIPTAG ) THEN
                  Spc(I,J,1,53) = Spc(I,J,1,53) + E_CO2
               ENDIF
            ENDIF

	   ! Add to Species #12: Background with Fossil CO2 [kg/kg]
           !-------------------------------------------
	   !IF ( LFFBKGRD ) THEN
           !   Spc(I,J,1,12)   = Spc(I,J,1,12) + E_CO2
	   !ENDIF
           !-------------------------------------------
	   ! Uncomment to include ship CO2 emissions in the background

         ENDIF
	   
         !-------------------------------------------
         ! #9: CO2 from aviation emissions
         !-------------------------------------------
         IF ( LPLANE ) THEN

            DO L = 1, LLPAR

               ! Get emissions [kg/m2/s] and convert to [molec/cm2/s]
               E_CO2 = HcoState%Spc(IDav)%Emis%Val(I,J,L) 
     &               / CM2PERM2 
     &               * XNUMOL_CO2

               ! ND04 diag: Aviation CO2 [molec/cm2/s]  
               ! !Unit change in v10 
               IF ( ND04 > 0 ) THEN
                  AD04_plane(I,J,L) = AD04_plane(I,J,L) + E_CO2
               ENDIF

! old
!               ! Convert E_CO2 from [molec CO2/cm2/s] to [kg]   
!               ! Unit change in v10
!               E_CO2 = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
! new
               ! Get emissions [kg/m2/s] and convert to [kg/kg dry air]
               ! using method in ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod
               ! (ewl, 9/11/15)
               E_CO2  =  HcoState%Spc(IDav)%Emis%Val(I,J,L)
     &                   * DTSRCE / ( State_Met%DELP(I,J,L) 
     &                   * G0_100 * ( 1.0e+0_fp
     &                   - State_Met%SPHU(I,J,L) * 1.0e-3_fp ) )        
! end new (ewl)

               ! Add to Species #1: Total CO2 [kg/kg]
               Spc(I,J,L,1)   = Spc(I,J,L,1) + E_CO2

               ! Add to Species #9: Aviation CO2 [kg/kg]
               IF ( nAdvect > 8 ) THEN
                  Spc(I,J,L,9) = Spc(I,J,L,9) + E_CO2
 
              	  ! Tagged species for global aviation emissions
                  IF ( LPLANETAG ) THEN
                     Spc(I,J,L,54) = Spc(I,J,L,54) + E_CO2
                  ENDIF
               ENDIF
	
	      ! Add to Species #12: Background with Fossil CO2 [kg]
              !-------------------------------------------
	      !IF (LFFBKGRD) THEN
	      !   Spc(I,J,L,12)   = Spc(I,J,L,12) + E_CO2
	      !ENDIF
              !-------------------------------------------
	      ! Uncomment to include aviation CO2 emissions in the background

            ENDDO
         ENDIF
 
         !-------------------------------------------
         ! #10 CO2 production from CO oxidation
         !-------------------------------------------
         IF ( LCHEMCO2 ) THEN

            DO L = 1, LLPAR
         
               ! Production is in [kg/m3], convert to [molec/cm2/s]
               E_CO2 = CO2_COPROD(I,J,L)
     &               / CM3PERM3
     &               * XNUMOL_CO2
     &               / DTSRCE
     &               * 100 !per area (cm2) instead of per volume (cm3)

	       ! ND04 diag: CO2 chem source from CO loss (3-D) [molec/cm2/s] 
               !Unit change in v10 
               IF ( ND04 > 0 ) THEN
                  AD04_chem(I,J,L) = AD04_chem(I,J,L) + E_CO2
               ENDIF

! old
!               ! Convert E_CO2 from [molec CO2/cm2/s] to [kg]   
!               ! Unit change in v10
!               E_CO2 = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
! new
               ! Convert emissions from [molec/cm2/s] to [kg/kg dry air]
               ! (ewl, 9/11/15)
               E_CO2  =  E_CO2 * DTSRCE * CM2PERM2 /
     &                   ( XNUMOL_CO2 * State_Met%DELP(I,J,L) 
     &                   * G0_100 * ( 1.0e+0_fp
     &                   - State_Met%SPHU(I,J,L) * 1.0e-3_fp ) )        
! end new (ewl)

               ! Add to Species #1: Total CO2 [kg/kg]
               Spc(I,J,L,1)   = Spc(I,J,L,1) + E_CO2

               ! Add to Species #10: Chemical Source of CO2 [kg/kg]
               IF (nAdvect > 9) THEN
                  Spc(I,J,L,10)   = Spc(I,J,L,10) + E_CO2
               ENDIF
		
            ENDDO
         ENDIF

         !-------------------------------------------
         ! #11 CO2 surface correction for CO oxidation
         !-------------------------------------------
         IF ( LCHEMCO2 ) THEN
            
            ! Get emissions [kg/m2/s] and convert to [molec/cm2/s]
            ! Emissions are only defined at surface, so set L=1.
            E_CO2 = HcoState%Spc(IDcorr)%Emis%Val(I,J,1) 
     &            / CM2PERM2 
     &            * XNUMOL_CO2

            ! ND04 diag: CO2 chem source surface correction [molec/cm2/s] 
            IF ( ND04 > 0 ) THEN
               AD04(I,J,10) = AD04(I,J,10) - E_CO2 ! SUBTRACT
            ENDIF

! old
!            ! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
!            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
! new
            ! Get emissions [kg/m2/s] and convert to [kg/kg dry air]
            ! using method in ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod
            ! (ewl, 9/11/15)
            E_CO2  =  HcoState%Spc(IDcorr)%Emis%Val(I,J,1)
     &                * DTSRCE / ( State_Met%DELP(I,J,1) 
     &                * G0_100 * ( 1.0e+0_fp
     &                - State_Met%SPHU(I,J,1) * 1.0e-3_fp ) )        
! end new (ewl)

            ! Subtract from Species #1: Total CO2 [kg/kg]
            Spc(I,J,1,1)   = Spc(I,J,1,1) - E_CO2

            ! Subtract from Species #11: Chem Source Surf Correction [kg/kg]
            IF ( nAdvect > 10 ) THEN
               Spc(I,J,1,11)   = Spc(I,J,1,11) - E_CO2
            ENDIF
		
         ENDIF

         !-------------------------------------------
         ! #12: Background CO2
         !-------------------------------------------
         ! Background CO2 without fossil fuels is obtained by setting 
         ! species 12 in the restart file, equal to species 1 at the start of 
         ! a run.  No sources or sinks (chemical or surface) act on this 
         ! species if LFFBKGRD == FALSE, it is simply the advected intial CO2.
	 ! To include fossil fuels in the background, see the first loop.

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

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
      SUBROUTINE DEF_BIOSPH_CO2_REGIONS_F( REGION )
!
! !USES:
!
      USE CMN_SIZE_MOD        ! Size parameters
      USE FILE_MOD,      ONLY : IOERROR
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT)   :: REGION(IIPAR,JJPAR)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  WARNING! Tagged CO2 simulation only work for 2 x 2.5 grid!  %%%
!  %%%  Someone will have to make this more general later on...     %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  13 Apr 2013 - R. Nassar - minor fixes to region boundary indices
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, IOS, IU_FILE
      INTEGER                :: TMP(IIPAR,JJPAR)
      INTEGER                :: LAND_REG(IIPAR,JJPAR)
      CHARACTER(LEN=255)     :: LANDFILE
      CHARACTER(LEN=144)     :: ROW
      CHARACTER(LEN=1)       :: CHAR1(IIPAR,JJPAR)

      !=================================================================
      ! Reading LAND BIOSPHERE REGIONS
      !=================================================================
      
      LANDFILE  = 'Regions_land.dat'

      WRITE(*,*) ' '
 100  FORMAT( '     - READ_REGIONS: Reading ', a ) 
      WRITE( 6, 100 ) TRIM( LANDFILE )

      ! Initialize ARRAY
      LAND_REG = 0

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open file
      OPEN( IU_FILE, FILE = TRIM( LANDFILE ), 
     &      FORM='FORMATTED', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:1' )      

      ! Read data
      DO J = 1, JJPAR
         IF (IIPAR ==  72) READ( IU_FILE, '(72A)', IOSTAT=IOS ) ROW
         IF (IIPAR == 144) READ( IU_FILE,'(144A)', IOSTAT=IOS ) ROW
         WRITE (*,'(A)') ROW

         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:2' )

         DO I = 1, IIPAR
            CHAR1(I,J) = ROW(I:I)
            IF (CHAR1(I,J) == ' ') CHAR1(I,J) = '0'
            READ (CHAR1(I,J),'(I1)') TMP(I,J)
         ENDDO
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Flip array in the North-South Direction
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         LAND_REG(I,J) = TMP(I,JJPAR-J+1)
      ENDDO
      ENDDO
      WRITE(*,*) ' '

      !=================================================================
      ! Loop over entire globe -- multiprocessor
      !=================================================================

      DO J = 1, JJPAR
      DO I = 1, IIPAR
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
      SUBROUTINE DEF_OCEAN_CO2_REGIONS_F( REGION )
!
! !USES:
!
      USE CMN_SIZE_MOD        ! Size parameters
      USE FILE_MOD,      ONLY : IOERROR
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT)  :: REGION(IIPAR,JJPAR)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  WARNING! Tagged CO2 simulation only work for 2 x 2.5 grid!  %%%
!  %%%  Someone will have to make this more general later on...     %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER               :: I, J, IU_FILE, IOS
      INTEGER               :: TMP(IIPAR,JJPAR), OCEAN_REG(IIPAR,JJPAR)
      CHARACTER(LEN=255)    :: OCEANFILE
      CHARACTER(LEN=144)    :: ROW
      CHARACTER(LEN=1)      :: CHAR1(IIPAR,JJPAR)

      !=================================================================
      ! DEF_CO2_OCEAN_REGIONS begins here!
      !=================================================================
      
      OCEANFILE = 'Regions_ocean.dat'

      WRITE( 6, 100 ) TRIM( OCEANFILE )
 100  FORMAT( '     - READ_REGIONS: Reading ', a ) 
      WRITE(*,*) ' '

      ! Initialize ARRAYS
      OCEAN_REG = 0

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open file
      OPEN( IU_FILE, FILE = TRIM( OCEANFILE ), 
     &      FORM='FORMATTED', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:1' )      

      ! Read data
      DO J = 1, JJPAR
         IF (IIPAR ==  72) READ( IU_FILE, '(72A)', IOSTAT=IOS ) ROW
         IF (IIPAR == 144) READ( IU_FILE,'(144A)', IOSTAT=IOS ) ROW
         WRITE (*,'(A)') ROW

         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:2' )

         DO I = 1, IIPAR
            CHAR1(I,J) = ROW(I:I)
            IF (CHAR1(I,J) == ' ') CHAR1(I,J) = '0'
            READ (CHAR1(I,J),'(I1)') TMP(I,J)
         ENDDO
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Flip array in the North-South Direction
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         OCEAN_REG(I,J) = TMP(I,JJPAR-J+1)
      ENDDO
      ENDDO
      WRITE(*,*) ' '

      !=================================================================
      ! Loop over entire globe -- multiprocessor
      !=================================================================

      DO J = 1, JJPAR
      DO I = 1, IIPAR
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
      SUBROUTINE DEF_FOSSIL_CO2_REGIONS_F( REGION )
!
! !USES:
!
      USE CMN_SIZE_MOD        ! Size parameters
      USE FILE_MOD,      ONLY : IOERROR
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT)   :: REGION(IIPAR,JJPAR)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  WARNING! Tagged CO2 simulation only work for 2 x 2.5 grid!  %%%
!  %%%  Someone will have to make this more general later on...     %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
! !REVISION HISTORY: 
!  18 May 2010 - R. Nassar, D. Jones - Initial version
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, IU_FILE, IOS
      INTEGER                :: TMP(IIPAR,JJPAR), REG_CODE(IIPAR,JJPAR)
      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=144)     :: ROW
      CHARACTER(LEN=1)       :: CHAR1(IIPAR,JJPAR)

      !=================================================================
      ! DEF_CO2_FOSSIL_REGIONS begins here!
      !=================================================================
      
      FILENAME  = 'Regions_land.dat'

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_REGIONS: Reading ', a ) 

      ! Initialize ARRAYS
      REG_CODE = 0

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open file
      OPEN( IU_FILE, FILE = TRIM( FILENAME ), 
     &      FORM='FORMATTED', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:1' )      

      ! Read data
      DO J = 1, JJPAR
         IF (IIPAR ==  72) READ( IU_FILE, '(72A)', IOSTAT=IOS ) ROW
         IF (IIPAR == 144) READ( IU_FILE,'(144A)', IOSTAT=IOS ) ROW
         WRITE (*,'(A)') ROW
                
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_regions:2' )

         DO I = 1, IIPAR
            CHAR1(I,J) = ROW(I:I)
            IF (CHAR1(I,J) == ' ') CHAR1(I,J) = '0'
            READ (CHAR1(I,J),'(I1)') TMP(I,J)
         ENDDO
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Flip array in the North-South Direction
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         REG_CODE(I,J) = TMP(I,JJPAR-J+1)
      ENDDO
      ENDDO

      !=================================================================
      ! Loop over entire globe -- multiprocessor
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
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
      SUBROUTINE INIT_CO2( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      ! References to F90 modules
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE Input_Opt_Mod,      ONLY : OptInput
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
! !REVISION HISTORY: 
!  16 Aug 2005 - P. Suntharalingam   - Initial version
!  18 May 2010 - R. Nassar, D. Jones - Updated 
!  25 Mar 2013 - R. Yantosca - Now use logical fields from Input_Opt
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
      ALLOCATE( FOSSIL_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FOSSIL_REGION' )
      FOSSIL_REGION = 0

      ! Array for Biospheric regions
      ALLOCATE( BIOSPH_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOSPH_REGION' )
      BIOSPH_REGION = 0

      ! Array for Ocean Regions
      ALLOCATE( OCEAN_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCEAN_REGION' )
      OCEAN_REGION = 0

      !=================================================================
      ! Set up regions 
      !=================================================================

      ! Set up tagged regions for balanced biosphere & ocean
      IF ( Input_Opt%LBIOSPHTAG ) THEN
         CALL DEF_BIOSPH_CO2_REGIONS_F( BIOSPH_REGION )
         CALL DEF_OCEAN_CO2_REGIONS_F( OCEAN_REGION )
      ENDIF

      ! Set up tagged regions for fossil fuel
      IF ( Input_Opt%LFOSSILTAG ) THEN
         CALL DEF_FOSSIL_CO2_REGIONS_F( FOSSIL_REGION )
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
!  18 May 2010 - R. Nassar, D. Jones - Updated 
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
