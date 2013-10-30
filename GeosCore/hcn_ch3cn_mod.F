! $Id: hcn_ch3cn_mod.f,v 1.2 2010/03/15 19:33:23 ccarouge Exp $
      MODULE HCN_CH3CN_MOD
!
!******************************************************************************
!  Module HCN_CH3CN_MOD contains variables and routines that are used for the 
!  geographically tagged HCN/CH3CN simulation. (qli, xyp, bmy, 8/16/05,9/27/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) HCN_BB_REGION     : Array to denote tagged HCN biomass tracers
!  (2 ) HCN_DF_REGION     : Array to denote tagged HCN fossil fuel tracers
!  (3 ) CH3CN_BB_REGION   : Array to denote tagged CH3CN biomass tracers
!  (4 ) CH3CN_DF_REGION   : Array to denote tagged CH3CN fossil fuel tracers
!  (5 ) EMIS_CO_df        : Array for CO from domestic fossil fuel
!  (6 ) HCN_INDEX         : Index array for HCN tracers
!  (7 ) CH3CN_INDEX       : Index array for CH3CN tracers
!  (8 ) SCNR89            : Weekday/weekend scenarios for fossil fuel scaling
!  (9 ) TODH              : Time of day scale factor for hydrocarbon emissions
!  (10) TODN              : Time of day scale factor for NOx emissions
!  (11) TODB              : Time of day scale factor for biogenic emissions
!
!  Module Routines:
!  ============================================================================
!  (1 ) DEFINE_BB_REGIONS : Defines geographic regions for biomass burn
!  (2 ) DEFINE_DF_REGIONS : Defines geographic regions for fossil fuel
!  (3 ) EMISS_HCN_CH3CN   : Emits into geographically "tagged" tracers
!  (4 ) CHEM_HCN_CH3CN    : Does chemistry for "tagged" tracers
!  (5 ) INIT_HCN_CH3CN    : Allocates and initializes module arrays
!  (6 ) CLEANUP_HCN_CH3CN : Deallocates module arrays
!
!  GEOS-Chem modules referenced by hcn_ch3cn_mod.f
!  ============================================================================
!  (1 ) biomass_mod.f     : Module w/ routines to read biomass emissions
!  (2 ) bpch2_mod.f       : Module w/ routines for binary punch file I/O
!  (3 ) dao_mod.f         : Module w/ arrays for DAO met fields!
!  (4 ) diag_mod.f        : Module w/ GEOS-Chem diagnostic arrays
!  (5 ) directory_mod.f   : Module w/ GEOS-Chem data & met field dirs
!  (6 ) geia_mod.f        : Module w/ routines to read anthro emissions
!  (7 ) global_oh_mod.f   : Module w/ routines to read 3-D OH field
!  (8 ) grid_mod.f        : Module w/ horizontal grid information
!  (9 ) global_oh_mod.f   : Module w/ routines to read 3-D OH field
!  (10) logical_mod.f     : Module w/ GEOS-Chem logical switches
!  (11) pbl_mix_mod.f     : Module w/ routines for PBL height & mixing
!  (12) time_mod.f        : Module w/ routines for computing time & date
!  (13) tracerid_mod.f    : Module w/ pointers to tracers & emissions
!  (14) transfer_mod.f    : Module w/ routines to cast & resize arrays
!  
! 
!  Tagged HCN/CH3CN tracers:
!  ============================================================================
!  (1 ) Total HCN
!  (2 ) Total CH3CN
!  (3 ) HCN from Asian biomass burning
!  (4 ) HCN from elsewhere biomass burning 
!  (5 ) HCN from Asian domestic fossil fuel 
!  (6 ) HCN from elsewhere domestic fossil fuel
!  (7 ) CH3CN from Asian biomass burning
!  (8 ) CH3CN from elsewhere biomass burning 
!  (9 ) CH3CN from Asian domestic fossil fuel 
!  (10) CH3CN from elsewhere domestic fossil fuel
!
!  References:
!  ============================================================================
!  (1 ) Li, Q.B., D.J. Jacob, R.M. Yantosca, C.L. Heald, H.B. Singh, M. Koike, 
!        Y.Zhao, G.W. Sachse, and D.G. Streets, "A Global 3-D Model Evaluation 
!        of the Atmospheric Budgets of HCN and CH3CN: Constraints From 
!        Aircraft Measurements Over the Western Pacific", J. Geophys. Res., 
!        108(D21), 2003
!  (2 ) Nightingale et al [2000a], J. Geophys. Res, 14, 373-387
!  (3 ) Nightingale et al [2000b], Geophys. Res. Lett, 27, 2117-2120
!
!  NOTES:
!  (1 ) Now use Nightingale et al [2000b] formulation for KL (bmy, 8/16/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Remove duplicate variable declarations for Linux IFORT v9 compiler
!        (bmy, 11/2/05)
!  (4 ) Now modified for new "biomass_mod.f" (bmy, 4/5/06)
!  (5 ) BIOMASS(:,:,IDBCO) from "biomass_mod.f" is now in units of 
!        [molec CO/cm2/s].  Adjust unit conversion accordingly. (bmy, 9/27/06)
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "tagged_hcn_ch3cn_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CHEM_HCN_CH3CN
      PUBLIC :: CLEANUP_HCN_CH3CN
      PUBLIC :: EMISS_HCN_CH3CN

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      REAL*8,  PARAMETER   :: MAIR         = 28.96d-3           ! kg/mol
      REAL*8,  PARAMETER   :: MHCN         = 27d-3              ! kg/mol
      REAL*8,  PARAMETER   :: MCH3CN       = 41d-3              ! kg/mol
      REAL*8,  PARAMETER   :: XNUMOL_AIR   = 6.022d23 / MAIR    ! molec/kg
      REAL*8,  PARAMETER   :: XNUMOL_HCN   = 6.022d23 / MHCN    ! molec/kg
      REAL*8,  PARAMETER   :: XNUMOL_CH3CN = 6.022d23 / MCH3CN  ! molec/kg

      ! Allocatable arrays
      INTEGER, ALLOCATABLE :: HCN_REG_bb(:,:)
      INTEGER, ALLOCATABLE :: HCN_REG_df(:,:)
      INTEGER, ALLOCATABLE :: CH3CN_REG_bb(:,:)
      INTEGER, ALLOCATABLE :: CH3CN_REG_df(:,:)
      REAL*8,  ALLOCATABLE :: EMIS_CO_df(:,:)

      ! Fixed-size arrays
      INTEGER              :: HCN_INDEX(5)
      INTEGER              :: CH3CN_INDEX(5)
      REAL*8               :: SCNR89(3,3)
      REAL*8               :: TODH(6)
      REAL*8               :: TODN(6)
      REAL*8               :: TODB(6)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS     

!------------------------------------------------------------------------------

      SUBROUTINE DEFINE_BB_REGIONS
!
!******************************************************************************
!  Subroutine DEFINE_BB_REGIONS defines the geographic regions for biomass 
!  burning emissions for the tagged HCN/CH3CN simulation. (xyp, bmy, 6/30/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) REGION (INTEGER) : Array of Fossil Fuel CO regions
! 
!  NOTES:
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90 
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90 
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID, GET_YMID

      USE CMN_SIZE_MOD    ! Size parameters

      ! Local variables
      INTEGER              :: I, J
      REAL*8               :: X, Y

      !=================================================================
      ! DEFINE_BB_REGIONS begins here!
      !=================================================================

      ! Loop over latitudes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, X, Y )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Longitude and latitude [degrees]
         X = GET_XMID( I, J, 1 )
         Y = GET_YMID( I, J, 1 )

         ! Region #3: SE Asian BB HCN (1st sub-box)
         IF      ( ( X >= 72.5 .AND. X < 127.5 )  .AND.
     &             ( Y >=  8.0 .AND. Y <  28.0 ) ) THEN
            HCN_REG_bb(I,J) = 3

         ! Region #3: SE Asian HCN BB (2nd sub-box)
         ELSE IF ( ( X >= 72.5 .AND. X < 152.5 )  .AND.
     &             ( Y >= 28.0 .AND. Y <  48.0 ) ) THEN
            HCN_REG_bb(I,J) = 3
  
         ! Region #4: HCN BB from elsewhere
         ELSE
            HCN_REG_bb(I,J) = 4

         ENDIF

         ! CH3CN tracer #'s are HCN tagged tracers + 4
         CH3CN_REG_bb(I,J)  = HCN_REG_bb(I,J) + 4

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DEFINE_BB_REGIONS

!------------------------------------------------------------------------------

      SUBROUTINE DEFINE_DF_REGIONS
!
!******************************************************************************
!  Subroutine DEFINE_DF_REGIONS defines the geographic regions for domestic 
!  fossil fuel emissions for the HCN/CH3CN simulation. (xyp, bmy, 6/30/05)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) REGION (INTEGER) : Array of Fossil Fuel regions 
!
!  NOTES
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID, GET_YMID

      USE CMN_SIZE_MOD    ! Size parameters

      ! Local variables
      INTEGER              :: I, J
      REAL*8               :: X, Y

      !=================================================================
      ! DEFINE_DF_REGIONS begins here!
      !=================================================================

      ! Loop over latitudes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, X, Y )
      DO J = 1, JJPAR
         
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Longitude and latitude [degrees]
            X = GET_XMID( I, J, 1 )
            Y = GET_YMID( I, J, 1 )         
         
            ! Region #5: HCN Asian DF (1st sub-box)
            IF      ( ( X >= 72.5 .AND. X < 127.5 )  .AND.
     &                ( Y >=  8.0 .AND. Y <  28.0 ) ) THEN
               HCN_REG_df(I,J) = 5

            ! Region #5: HCN Asian DF (2nd sub-box)
            ELSE IF ( ( X >= 72.5 .AND. X < 152.5 )  .AND.
     &             ( Y >= 28.0 .AND. Y <  48.0 ) ) THEN
               HCN_REG_df(I,J) = 5
   
            ! Region #6: HCN DF from elsewhere
            ELSE
               HCN_REG_df(I,J) = 6
               
            ENDIF

            ! CH3CN tracer #'s are HCN tagged tracers + 4
            CH3CN_REG_df(I,J)  = HCN_REG_df(I,J) + 4

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DEFINE_DF_REGIONS

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_HCN_CH3CN( am_I_Root, Input_Opt,
     &                            State_Met, State_Chm, RC )
!
!******************************************************************************
!  Subroutine EMISS_HCN_CH3CN reads in CO emissions and scale them to get
!  HCN/CH3CN emissions for the tagged HCN/CH3CN run. (bmy, 8/16/05, 9/27/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N_TRACERS (INTEGER) : Number of tracers
!  (2 ) STT       (REAL*8 ) : Tracer array [kg]
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now modified for new "biomass_mod.f" (bmy, 4/5/06)
!  (3 ) BIOMASS(:,:,IDBCO) from "biomass_mod.f" is now in units of 
!        [molec CO/cm2/s].  Adjust unit conversion accordingly. (bmy, 9/27/06)
!  (4 ) IDBCO is now defined in TRACERID_MOD. (fp, 6/2009)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!******************************************************************************
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE BIOMASS_MOD,   ONLY : BIOMASS!,         IDBCO
      USE TRACERID_MOD,  ONLY : IDBCO
      USE GEIA_MOD,      ONLY : GET_DAY_INDEX,   GET_IHOUR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE DIAG_MOD,      ONLY : AD09_em
      USE PBL_MIX_MOD,   ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,      ONLY : GET_MONTH,       GET_TS_CHEM,  GET_TAU 
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
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I,      J,      L,        N, NTAU
      INTEGER                :: IHOUR,  INDEX,  MONTH,    PBL_MAX 
      REAL*8                 :: ACM2,   E_CObb, E_COdf,   SFAC89 
      REAL*8                 :: HCN_bb, HCN_df, CH3CN_bb, CH3CN_df
      REAL*8                 :: DTSRCE, FRAC

      ! Emission ratios for HCN/CH3CN from biomass burning 
      ! and domestic fossil fuel
      REAL*8,  PARAMETER     :: EHCN_bb   = 0.27d-2
      REAL*8,  PARAMETER     :: EHCN_df   = 1.60d-2
      REAL*8,  PARAMETER     :: ECH3CN_bb = 0.20d-2
      REAL*8,  PARAMETER     :: ECH3CN_df = 0.25d-2

      ! For fields from Input_Opt
      LOGICAL                :: LSPLIT

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! EMISS_TAGGED_HCN_CH3CN begins here!
      !=================================================================

      ! Assume success
      RC      = GIGC_SUCCESS

      ! Copy values from Input_Opt
      LSPLIT  = Input_Opt%LSPLIT

      ! DTSRCE is the number of seconds per emission timestep
      DTSRCE  = GET_TS_CHEM() * 60d0

      ! Get the highest extent of the PBL [levels]
      PBL_MAX = GET_PBL_MAX_L()

      ! Get the current month
      MONTH   = GET_MONTH()

      ! Current TAU value (integer)
      NTAU    = GET_TAU()

      ! First-time initialization
      IF ( FIRST ) THEN 
         CALL INIT_HCN_CH3CN
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Process biomass burning/domestic fossil fuel HCN/CH3CN emissions
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J, ACM2, E_CObb, INDEX,  SFAC89,   E_COdf   )
!$OMP+PRIVATE( IHOUR, N, L,    HCN_bb, HCN_df, CH3CN_bb, CH3CN_df )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid area in [cm2]
         ACM2 = GET_AREA_CM2( I, J, 1 )

         !-----------------------------------------------------------------
         ! (1) Process biomass burning HCN/CH3CN emissions
         !-----------------------------------------------------------------

         ! Get CO biomass burning [molec CO/cm2/s]
         E_CObb = BIOMASS(I,J,IDBCO)

         ! ND09: biomass burning HCN/CH3CN emissions [molec/cm2/s]
         IF ( ND09 > 0 ) THEN
            AD09_em(I,J,1) = AD09_em(I,J,1) + ( EHCN_bb   * E_CObb )
            AD09_em(I,J,2) = AD09_em(I,J,2) + ( ECH3CN_bb * E_CObb )
         ENDIF

         ! Convert [molec CO/cm2/s] to [mole/grid box]: 1/6.022d23 = 1.66d-24
         E_CObb = E_CObb * 1.66d-24 * ACM2 * DTSRCE

         !-----------------------------------------------------------------
         ! (2) Process domestic fossil fuel HCN/CH3CN emissions
         !-----------------------------------------------------------------

         ! SFAC89 is the Weekday/Saturday/Sunday scale factor
         INDEX     = GET_DAY_INDEX( NTAU )
         SFAC89    = SCNR89( 2, INDEX ) 

         ! E_COdf is DF CO emissions in [molec CO/cm2/s]
         ! Scale E_COdf by the day-of-week scale factor SFAC89
         E_COdf    = EMIS_CO_df(I,J) * SFAC89

         ! Scale E_COdf by the time-of-day scale factor TODH
         ! IHOUR is the index for the time-of-day scale factor TODH
         IHOUR     = GET_IHOUR( I )
         E_COdf    = E_COdf * TODH(IHOUR)

         ! Enhance E_COdf by 18.5% to account for oxidation 
         ! from anthropogenic VOC's (bnd, bmy, 6/8/01)
         E_COdf    = E_COdf * 1.185d0
            
         ! Get HCN domestic fossil fuel region # (either =5 or =6)
         N         = HCN_REG_df(I,J)

         ! To achieve the best fit to the observed HCN-CH3CN-CO correlations 
         ! in the boundary layer, we have to double the residential coal 
         ! burning source from Asia. This leads us to reduce the residential 
         ! coal burning source from the rest of the world by a factor of eight
         ! to achieve a best fit to the observed vertical distributions of HCN
         ! and CH3CN. [According to Li et al 2003.] (xyp, 6/22/05)
         IF ( N == 5 ) THEN
            E_COdf = E_COdf * 2.1d0   ! Asian domestic fossil fuel
         ELSE
            E_COdf = E_COdf / 8.0d0   ! Elsewhere domestic fossil fuel
         ENDIF

         ! ND09: domestic fossil fuel HCN/CH3CN emissions [molec/cm2/s]
         IF ( ND09 > 0 ) THEN
            AD09_em(I,J,3) = AD09_em(I,J,3) + ( EHCN_df   * E_COdf )
            AD09_em(I,J,4) = AD09_em(I,J,4) + ( ECH3CN_df * E_COdf )
         ENDIF

         ! Convert [molec CO/cm2/s] to [mole/grid box]: 1/6.022d23 = 1.66d-24
         E_COdf    = E_COdf * 1.66d-24 * ACM2 * DTSRCE       

         !-----------------------------------------------------------------
         ! (3) Partition emissions throughout the boundary layer
         !-----------------------------------------------------------------

         ! Loop up to the highest PBL level
         DO L = 1, PBL_MAX
            
            ! Fraction of the PBL occupied by this layer
            FRAC            = GET_FRAC_OF_PBL( I, J, L )

            ! HCN biomass and domestic fossil fuel emissions
            HCN_bb          = FRAC * MHCN * EHCN_bb * E_CObb
            HCN_df          = FRAC * MHCN * EHCN_df * E_COdf

            ! CH3CN biomass and domestic fossil fuel emissions
            CH3CN_bb        = FRAC * MCH3CN * ECH3CN_bb * E_CObb
            CH3CN_df        = FRAC * MCH3CN * ECH3CN_df * E_COdf

            ! Add total HCN emissions (BB+DF) into STT
            STT(I,J,L,1)    = STT(I,J,L,1) + ( HCN_bb   + HCN_df  )

            ! Add total CH3CN emissions (BB+DF) into STT
            STT(I,J,L,2)    = STT(I,J,L,2) + ( CH3CN_bb + CH3CN_df )

            ! If we are using tagged tracers ...
            IF ( LSPLIT ) THEN
               
               ! Add emissions into tagged HCN biomass tracers
               N            = HCN_REG_bb(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + HCN_bb

               ! Add emissions into tagged HCN dom. fossil tracers
               N            = HCN_REG_df(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + HCN_df

               ! Add emissions into tagged CH3CN biomass tracers
               N            = CH3CN_REG_bb(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + CH3CN_bb

               ! Add emissions into tagged CH3CN dom. fossil tracers
               N            = CH3CN_REG_df(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + CH3CN_df

            ENDIF
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      ! Return to calling program
      END SUBROUTINE EMISS_HCN_CH3CN

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_HCN_CH3CN( am_I_Root, Input_Opt,
     &                           State_Met, State_Chm, RC )
!
!******************************************************************************
!  Subroutine CHEM_HCN_CH3CN computes the loss of HCN and CH3CN due to 
!  reaction with OH and ocean uptake. (xyp, bmy, 8/16/05, 11/2/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FIRSTCHEM (LOGICAL) : = T if this is the first call to this routine
!
!  NOTES:
!  (1 ) Now use Nightingale et al [2000b] formulation for KL (bmy, 8/16/05)
!  (2 ) Bug fix: remove duplicate declaration of KTMP (bmy, 11/2/05)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  23 Oct 2013 - R. Yantosca - Now pass objects to GET_GLOBAL_OH routine
!******************************************************************************
! 
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD09, AD09_em
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_OH_MOD,      ONLY : OH, GET_GLOBAL_OH
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_MONTH
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
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
! !LOCAL VARIABLES:
!      
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I,      J,    L,    N,     NN,  N_MAX
      REAL*8                 :: K0,     K1,   KTMP, KRATE, TMP, DTCHEM
      REAL*8                 :: H,      U,    TC,   SC,    KL,  KG
      REAL*8                 :: KKG,    CL,   SR,   CG,    FLUX
      REAL*8                 :: ACM2,   AMT_LOST,   OCEAN_HCN
      REAL*8                 :: FOCEAN, OCEAN_CH3CN

      ! Undersaturation ratios for HCN/CH3CN in seawater
      REAL*8, PARAMETER      :: ALPHA_HCN   = 0.21d0
      REAL*8, PARAMETER      :: ALPHA_CH3CN = 0.12d0

      ! Coefficients for fitting the Schmdit number for HCN in seawater
      REAL*8, PARAMETER      :: A0          = 2008.917d0
      REAL*8, PARAMETER      :: A1          =  -83.235d0
      REAL*8, PARAMETER      :: A2          =    1.348d0
      REAL*8, PARAMETER      :: A3          =   -0.009d0
      
      ! Coefficients for fitting the Schmdit number for CH3CN in seawater
      REAL*8, PARAMETER      :: B0          = 2745.722d0
      REAL*8, PARAMETER      :: B1          = -113.763d0
      REAL*8, PARAMETER      :: B2          =    1.843d0
      REAL*8, PARAMETER      :: B3          =   -0.012d0

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER        :: STT(:,:,:,:)

      ! External functions
      REAL*8, EXTERNAL       :: BOXVL
     
      !=================================================================
      ! CHEM_HCN_CH3CN begins here! 
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! First-time initialization (if not already done)
      IF ( FIRST ) THEN
         CALL INIT_HCN_CH3CN
         FIRST = .FALSE.
      ENDIF

      ! Read offline OH fields once per month
      IF ( ITS_A_NEW_MONTH() ) THEN 
!----------------------------------------------------------------------------
! Prior to 10/23/13:
!         CALL GET_GLOBAL_OH( GET_MONTH(), State_Met )
!----------------------------------------------------------------------------
         CALL GET_GLOBAL_OH( am_I_Root, Input_Opt, 
     &                       State_Met, GET_MONTH(), RC )

      ENDIF
     
      ! Compute number of tracers to process
      IF ( Input_Opt%LSPLIT ) THEN
         N_MAX = 5
      ELSE
         N_MAX = 1
      ENDIF

      !=================================================================
      ! Do HCN and CH3CN chemistry
      !=================================================================

      ! Chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, K0, K1, TMP, KTMP, KRATE, NN, N, AMT_LOST )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !------------------------------------------------------------------
         ! (1) HCN loss via reaction with OH
         !------------------------------------------------------------------

         K0    = 7.4d-33 
         K1    = 9.0d-15    * ( State_Met%T(I,J,L)  / 300d0 ) ** 3.2d0
         TMP   = K0  /  K1  *   State_Met%AD(I,J,L) *
     &           XNUMOL_AIR / BOXVL(I,J,L,State_Met)

         ! K: [cm3/molec/s]
         KTMP  = K1 * TMP / ( 1d0 + TMP )      
     &         * EXP ( -0.511d0 / ( 1d0 + LOG10( TMP ) ** 2d0 ) )

         ! Rate constant for rxn w/ OH [units??]
         KRATE = KTMP * OH(I,J,L) * DTCHEM

         ! Subtract lost HCN from STT array
         DO NN = 1, N_MAX 

            ! Get the pr
            N = HCN_INDEX(NN)

            ! Compute the amount of tracer that is lost to OH
            AMT_LOST     = KRATE * STT(I,J,L,N)

            ! Remove lost tracer from STT array (avoid negatives!)
            STT(I,J,L,N) = MAX( STT(I,J,L,N) - AMT_LOST, 0d0 )
            
            ! ND09 diagnostic: HCN/CH3CN loss via OH [kg]
            IF ( ND09 > 0 ) THEN
               AD09(I,J,L,N) = AD09(I,J,L,N) + AMT_LOST
            ENDIF
         ENDDO

         !------------------------------------------------------------------
         ! (2) CH3CN loss via reaction with OH
         !------------------------------------------------------------------

         ! K: [cm3/molec/s]
         KTMP  = 7.8d-13 * EXP( -1050d0 / State_Met%T(I,J,L) )
         KRATE = KTMP * OH(I,J,L) * DTCHEM

         ! Subtract lost CH3CN tracer from STT
         DO NN = 1, N_MAX 

            ! Get the proper tracer number
            N = CH3CN_INDEX(NN)

            ! Compute the amount of tracer that is lost to OH
            AMT_LOST     = KRATE * STT(I,J,L,N)

            ! Remove lost CH3CN tracer from STT array (avoid negatives!)
            STT(I,J,L,N) = MAX( STT(I,J,L,N) - AMT_LOST, 0d0 )
            
            ! ND09 diagnostic: CH3CN loss via OH [kg]
            IF ( ND09 > 0 ) THEN
               AD09(I,J,L,N) = AD09(I,J,L,N) + AMT_LOST
            ENDIF
         ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! HCN and CH3CN ocean uptake
      !=================================================================

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,   J,  FOCEAN, OCEAN_HCN, OCEAN_CH3CN, ACM2     )
!$OMP+PRIVATE( U,   TC, H,      SC,        KL,          KG       )
!$OMP+PRIVATE( KKG, NN, N,      CG,        FLUX,        AMT_LOST )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Fraction of a grid box that is ocean
         FOCEAN          = 1d0 - State_Met%FRCLND(I,J) 

         ! Initialize HCN and CH3CN [kg] lost into the ocean
         OCEAN_HCN       = 0d0
         OCEAN_CH3CN     = 0d0

         ! Make sure there is > 50% ocean (not ice) in the grid box
         IF ( FOCEAN > 0.5d0 .AND. State_Met%ALBD(I,J) <= 0.4d0 ) THEN

            ! Grid box area in [cm2]
            ACM2         = GET_AREA_CM2( I, J, 1 ) 
            
            ! Wind speed [m/s] at 10m above the surface 
            U            = SQRT( State_Met%U10M(I,J)**2 + 
     &                           State_Met%V10M(I,J)**2 )

            ! Surface temperature [C]
            TC           = State_Met%TS(I,J) - 273.15d0  

            !-----------------------------------------------------------
            ! (1) HCN ocean uptake
            !-----------------------------------------------------------

            ! Henry's law constant for HCN [unitless]
            H            = 7.93d4 * EXP( -5000d0 / State_Met%TS(I,J) ) 
            
            ! SC is Schmidt # for HCN in seawater [unitless]
            SC           = A0 + TC * ( A1 + TC * ( A2 + TC * ( A3 )))

            ! KL: conductance for mass transfer in liquid phase 
            ! (Nightingale 2000b), which has unit of [cm/h]
            KL           = ( 0.24d0*U*U + 0.061d0*U ) * SQRT( 600d0/SC ) 

            ! KG: conductance for mass transfer in gas phase (Asher 1997)
            ! Convert from m/s to cm/h by multiplying 360000
            KG           = ( 15.3d0 + 940.6d0 * U ) 

            ! KKG: transfer velocity on a gas phase basis (Liss & Slater 1974)
            ! Convert from [cm/h] to [cm/s] by dividing 3600
            KKG          = 2.78d-4 * KL * KG / ( KL + H * KG )

            ! Loop over HCN tagged tracers
            DO NN = 1, N_MAX
               
               ! Get HCN tagged tracer number
               N            = HCN_INDEX(NN)

               ! Bulk concentration of HCN in gas phase [kg/cm3]
               CG           = STT(I,J,1,N) / BOXVL(I,J,1,State_Met)

               ! Air-to-sea flux of HCN [kg/cm2/s]
               FLUX         = ALPHA_HCN * KKG * CG     

               ! Amount of tagged tracer lost to ocean [kg]
               AMT_LOST     = FLUX * FOCEAN * ACM2 * DTCHEM

               ! Save total HCN lost to ocean for ND09 diag [molec/cm2/s]
               IF ( N == 1 ) THEN
                  OCEAN_HCN = AMT_LOST * XNUMOL_HCN / ( ACM2 * DTCHEM )
               ENDIF

               ! Subtract ocean loss from STT array [kg/box/step]
               STT(I,J,1,N) = MAX( STT(I,J,1,N) - AMT_LOST, 0d0 )

            ENDDO

            !-----------------------------------------------------------
            ! (2) CH3CN ocean uptake
            !-----------------------------------------------------------

            ! Henry's law constant for CH3CN [unitless]
            H            = 861.7d0 * EXP( -4100d0 / State_Met%TS(I,J) ) 

            ! SC is Schmidt # for HCN in seawater [unitless]
            SC           = B0 + TC * ( B1 + TC * ( B2 + TC * ( B3 )))

            ! KL: conductance for mass transfer in liquid phase
            ! (Wanninkhof 1992), which has units of [cm/h]
            KL           = ( 0.222d0 * U * U  + 0.333d0 * U )
     &                   * ( SC / 600d0 )**( -0.5d0 )

            ! KG: conductance for mass transfer in gas phase (Asher 1997)
            ! Convert from m/s to cm/h by mutiplying by 360000
            KG           = ( 12.4d0 + 763.3d0 * U ) 

            ! KKG: transfer velocity on a gas phase basis (Liss & Slater 1974)
            ! Convert from [cm/h] to [cm/s] by dividing by 3600
            KKG          = 2.78d-4 * KL * KG / ( KL + H * KG )

            ! Loop over CH3HCN tagged tracers
            DO NN = 1, N_MAX
               
               ! Get CH3CN tagged tracer number
               N              = CH3CN_INDEX(NN)

               ! Bulk concentration of CH3CN in gas phase [kg/cm3]
               CG             = STT(I,J,1,N) / BOXVL(I,J,1,State_Met)

               ! Air-to-sea flux of HCN [kg/cm2/s]
               FLUX           = ALPHA_HCN * KKG * CG     

               ! Amount of tagged tracer lost to ocean [kg]
               AMT_LOST       = FLUX * FOCEAN * ACM2 * DTCHEM

               ! Save total HCN lost to ocean for ND09 diag [molec/cm2/s]
               IF ( N == 2 ) THEN
                  OCEAN_CH3CN = AMT_LOST * XNUMOL_CH3CN / (ACM2*DTCHEM) 
               ENDIF

               ! Subtract ocean loss from STT array [kg/box/step]
               STT(I,J,1,N)   = MAX( STT(I,J,1,N) - AMT_LOST, 0d0 )

            ENDDO
         ENDIF

         !--------------------------------------------------------------
         ! ND10 diag: Save HCN and CH3CN ocean uptake in [molec/cm2/s]
         !--------------------------------------------------------------
         IF ( ND09 > 0 ) THEN
            AD09_em(I,J,5) = AD09_em(I,J,5) + OCEAN_HCN 
            AD09_em(I,J,6) = AD09_em(I,J,6) + OCEAN_CH3CN 
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      ! Return to calling program
      END SUBROUTINE CHEM_HCN_CH3CN

!------------------------------------------------------------------------------

      SUBROUTINE READ_EMISSIONS
!
!******************************************************************************
!  Subroutine READ_EMISSIONS reads the domestic fossil fuel emissions from
!  disk. (bmy, 6/29/05, 10/3/05)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_CO   (REAL*4) : GEIA anthro CO   (no seasonality, 1 level )
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE GEIA_MOD,      ONLY : READ_TODX
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD   ! Size parameters

      ! Local variables
      REAL*4              :: ARRAY(IIPAR,JJPAR,1)
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_EMISSIONS begins here!
      !=================================================================

      ! Define the binary punch file name
      FILENAME = TRIM( DATA_DIR )                         //
     &           'HCN_200507/domfos_CO_for_TRACEP.'       // 
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT() 
      
      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_EMISSIONS: Reading ', a )

      ! Read time-of-day and day-of-week scale factors for GEIA emissions
      CALL READ_TODX( TODN, TODH, TODB, SCNR89 )

      ! Read CO (tracer #4): aseasonal
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    4,  
     &                 0d0,       IIPAR,        JJPAR,     
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMIS_CO_df )

      ! Return to calling program
      END SUBROUTINE READ_EMISSIONS

!------------------------------------------------------------------------------

      SUBROUTINE INIT_HCN_CH3CN
!
!******************************************************************************
!  Subroutine INIT_TAGGED_HCN_CH3CN allocates memory to module arrays.
!  (bmy, 6/29/05)
! 
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD  ! Size parameters

      ! Local variables
      LOGICAL, SAVE      :: IS_INIT = .FALSE.
      INTEGER            :: AS
      
      !=================================================================
      ! INIT_TAGGED_CO begins here!
      !=================================================================

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      ! Allocate arrays
      ALLOCATE( HCN_REG_bb( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HCN_REG_bb' )         

      ALLOCATE( HCN_REG_df( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HCN_REG_df' )

      ALLOCATE( CH3CN_REG_bb( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH3CN_REG_bb' )         

      ALLOCATE( CH3CN_REG_df( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH3CN_REG_df' )
      
      ALLOCATE( EMIS_CO_df( IIPAR, JJPAR ), STAT=as )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMIS_CO_df' )

      ! Define geographic regions for biomass burning
      CALL DEFINE_BB_REGIONS

      ! Define geographic regions for domestic fossil fuel burning
      CALL DEFINE_DF_REGIONS

      ! Read domestic fossil fuel emissions
      CALL READ_EMISSIONS

      ! Index of HCN tracers
      HCN_INDEX(:)   = (/ 1, 3, 4, 5, 6  /)

      ! Index of CH3CN tracers
      CH3CN_INDEX(:) = (/ 2, 7, 8, 9, 10 /)

      ! Set flag
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_HCN_CH3CN

!------------------------------------------------------------------------------
  
      SUBROUTINE CLEANUP_HCN_CH3CN
!
!******************************************************************************
!  Subroutine CLEANUP_HCN_CH3CN deallocates memory from previously
!  allocated module arrays (bmy, 6/23/05)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_HCN_CH3CN begins here!
      !=================================================================
      IF ( ALLOCATED( HCN_REG_bb    ) ) DEALLOCATE( HCN_REG_bb   )
      IF ( ALLOCATED( HCN_REG_df    ) ) DEALLOCATE( HCN_REG_df   )
      IF ( ALLOCATED( CH3CN_REG_bb  ) ) DEALLOCATE( CH3CN_REG_bb )
      IF ( ALLOCATED( CH3CN_REG_df  ) ) DEALLOCATE( CH3CN_REG_df )
      IF ( ALLOCATED( EMIS_CO_df    ) ) DEALLOCATE( EMIS_CO_df   )

      ! Return to calling program
      END SUBROUTINE CLEANUP_HCN_CH3CN

!------------------------------------------------------------------------------

      ! End of module
      END MODULE HCN_CH3CN_MOD
