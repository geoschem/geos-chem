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
!  14 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  11 Mar 2015 - R. Yantosca - Now get OH from the HEMCO list
!******************************************************************************
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

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
      PUBLIC :: INIT_HCN_CH3CN

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      REAL(fp), PARAMETER   :: MAIR         = 28.96e-3_fp           ! kg/mol
      REAL(fp), PARAMETER   :: MHCN         = 27e-3_fp              ! kg/mol
      REAL(fp), PARAMETER   :: MCH3CN       = 41e-3_fp              ! kg/mol
      REAL(fp), PARAMETER   :: XNUMOL_AIR   = 6.022e+23_fp / MAIR    ! molec/kg
      REAL(fp), PARAMETER   :: XNUMOL_HCN   = 6.022e+23_fp / MHCN    ! molec/kg
      REAL(fp), PARAMETER   :: XNUMOL_CH3CN = 6.022e+23_fp / MCH3CN  ! molec/kg

      ! Allocatable arrays
      INTEGER,  ALLOCATABLE :: HCN_REG_bb(:,:)
      INTEGER,  ALLOCATABLE :: HCN_REG_df(:,:)
      INTEGER,  ALLOCATABLE :: CH3CN_REG_bb(:,:)
      INTEGER,  ALLOCATABLE :: CH3CN_REG_df(:,:)
      REAL(fp), ALLOCATABLE :: EMIS_CO_df(:,:)

      ! Fixed-size arrays
      INTEGER               :: HCN_INDEX(5)
      INTEGER               :: CH3CN_INDEX(5)
      REAL(fp)              :: SCNR89(3,3)
      REAL(fp)              :: TODH(6)
      REAL(fp)              :: TODN(6)
      REAL(fp)              :: TODB(6)

      ! Pointers to fields in the HEMCO data structure.
      ! These need to be declared REAL(f4), aka REAL*4.
      REAL(f4), POINTER     :: OH(:,:,:)   => NULL()

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
      REAL(fp)               :: X, Y

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
      REAL(fp)               :: X, Y

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
!  11 Apr 2014 - R. Yantosca - Remove call to INIT_HCN_CH3CN
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!******************************************************************************
! 
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD09, AD09_em
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
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
      REAL(fp)                 :: K0,     K1,   KTMP, KRATE, TMP, DTCHEM
      REAL(fp)                 :: H,      U,    TC,   SC,    KL,  KG
      REAL(fp)                 :: KKG,    CL,   SR,   CG,    FLUX
      REAL(fp)                 :: ACM2,   AMT_LOST,   OCEAN_HCN
      REAL(fp)                 :: FOCEAN, OCEAN_CH3CN
      REAL(fp)                 :: BOXVL

      ! Undersaturation ratios for HCN/CH3CN in seawater
      REAL(fp), PARAMETER      :: ALPHA_HCN   = 0.21e+0_fp
      REAL(fp), PARAMETER      :: ALPHA_CH3CN = 0.12e+0_fp

      ! Coefficients for fitting the Schmdit number for HCN in seawater
      REAL(fp), PARAMETER      :: A0          = 2008.917e+0_fp
      REAL(fp), PARAMETER      :: A1          =  -83.235e+0_fp
      REAL(fp), PARAMETER      :: A2          =    1.348e+0_fp
      REAL(fp), PARAMETER      :: A3          =   -0.009e+0_fp
      
      ! Coefficients for fitting the Schmdit number for CH3CN in seawater
      REAL(fp), PARAMETER      :: B0          = 2745.722e+0_fp
      REAL(fp), PARAMETER      :: B1          = -113.763e+0_fp
      REAL(fp), PARAMETER      :: B2          =    1.843e+0_fp
      REAL(fp), PARAMETER      :: B3          =   -0.012e+0_fp

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER        :: STT(:,:,:,:)

      ! Strings
      CHARACTER(LEN=255) :: LOC = 'CHEM_HCN_CH3CN (hcn_ch3cn_mod.F)'

      !=================================================================
      ! CHEM_HCN_CH3CN begins here! 
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Get a pointer to the OH field from the HEMCO list (bmy, 3/11/15)
      IF ( ITS_A_NEW_MONTH() ) THEN 
         CALL HCO_GetPtr( am_I_Root, 'GLOBAL_OH', OH,   RC )
         IF ( RC /= GIGC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to GLOBAL_OH', LOC )
         ENDIF
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
      DTCHEM = GET_TS_CHEM() * 60e+0_fp

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, K0, K1, TMP, KTMP, KRATE, NN, N, AMT_LOST )
!$OMP+PRIVATE( BOXVL                                              )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box volume [cm3]
         BOXVL = State_Met%AIRVOL(I,J,L) * 1e+6_fp

         !------------------------------------------------------------------
         ! (1) HCN loss via reaction with OH
         !------------------------------------------------------------------

         K0    = 7.4e-33_fp 
         K1    = 9.0e-15_fp    * ( State_Met%T(I,J,L)  / 300e+0_fp ) 
     &           ** 3.2e+0_fp
         TMP   = K0  /  K1  *   State_Met%AD(I,J,L) *
     &           XNUMOL_AIR / BOXVL

         ! K: [cm3/molec/s]
         KTMP  = K1 * TMP / ( 1e+0_fp + TMP )      
     &         * EXP ( -0.511e+0_fp / ( 1e+0_fp + LOG10( TMP ) 
     &         ** 2e+0_fp ) )

         ! Rate constant for rxn w/ OH [units??]
         KRATE = KTMP * OH(I,J,L) * DTCHEM

         ! Subtract lost HCN from STT array
         DO NN = 1, N_MAX 

            ! Get the pr
            N = HCN_INDEX(NN)

            ! Compute the amount of tracer that is lost to OH
            AMT_LOST     = KRATE * STT(I,J,L,N)

            ! Remove lost tracer from STT array (avoid negatives!)
            STT(I,J,L,N) = MAX( STT(I,J,L,N) - AMT_LOST, 0e+0_fp )
            
            ! ND09 diagnostic: HCN/CH3CN loss via OH [kg]
            IF ( ND09 > 0 ) THEN
               AD09(I,J,L,N) = AD09(I,J,L,N) + AMT_LOST
            ENDIF
         ENDDO

         !------------------------------------------------------------------
         ! (2) CH3CN loss via reaction with OH
         !------------------------------------------------------------------

         ! K: [cm3/molec/s]
         KTMP  = 7.8e-13_fp * EXP( -1050e+0_fp / State_Met%T(I,J,L) )
         KRATE = KTMP * OH(I,J,L) * DTCHEM

         ! Subtract lost CH3CN tracer from STT
         DO NN = 1, N_MAX 

            ! Get the proper tracer number
            N = CH3CN_INDEX(NN)

            ! Compute the amount of tracer that is lost to OH
            AMT_LOST     = KRATE * STT(I,J,L,N)

            ! Remove lost CH3CN tracer from STT array (avoid negatives!)
            STT(I,J,L,N) = MAX( STT(I,J,L,N) - AMT_LOST, 0e+0_fp )
            
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
!$OMP+PRIVATE( BOXVL                                             )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box volume [cm3]
         BOXVL           = State_Met%AIRVOL(I,J,L) * 1e+6_fp

         ! Fraction of a grid box that is ocean
         FOCEAN          = 1e+0_fp - State_Met%FRCLND(I,J) 

         ! Initialize HCN and CH3CN [kg] lost into the ocean
         OCEAN_HCN       = 0e+0_fp
         OCEAN_CH3CN     = 0e+0_fp

         ! Make sure there is > 50% ocean (not ice) in the grid box
         IF ( FOCEAN > 0.5e+0_fp .AND. State_Met%ALBD(I,J) 
     &      <= 0.4e+0_fp ) THEN

            ! Grid box area in [cm2]
            ACM2         = GET_AREA_CM2( I, J, 1 ) 
            
            ! Wind speed [m/s] at 10m above the surface 
            U            = SQRT( State_Met%U10M(I,J)**2 + 
     &                           State_Met%V10M(I,J)**2 )

            ! Surface temperature [C]
            TC           = State_Met%TS(I,J) - 273.15e+0_fp  

            !-----------------------------------------------------------
            ! (1) HCN ocean uptake
            !-----------------------------------------------------------

            ! Henry's law constant for HCN [unitless]
            H = 7.93e+4_fp * EXP( -5000e+0_fp / State_Met%TS(I,J) ) 
            
            ! SC is Schmidt # for HCN in seawater [unitless]
            SC = A0 + TC * ( A1 + TC * ( A2 + TC * ( A3 )))

            ! KL: conductance for mass transfer in liquid phase 
            ! (Nightingale 2000b), which has unit of [cm/h]
            KL = (0.24e+0_fp*U*U + 0.061e+0_fp*U) * SQRT(600e+0_fp/SC) 

            ! KG: conductance for mass transfer in gas phase (Asher 1997)
            ! Convert from m/s to cm/h by multiplying 360000
            KG = ( 15.3e+0_fp + 940.6e+0_fp * U ) 

            ! KKG: transfer velocity on a gas phase basis (Liss & Slater 1974)
            ! Convert from [cm/h] to [cm/s] by dividing 3600
            KKG = 2.78e-4_fp * KL * KG / ( KL + H * KG )

            ! Loop over HCN tagged tracers
            DO NN = 1, N_MAX
               
               ! Get HCN tagged tracer number
               N            = HCN_INDEX(NN)

               ! Bulk concentration of HCN in gas phase [kg/cm3]
               CG           = STT(I,J,1,N) / BOXVL

               ! Air-to-sea flux of HCN [kg/cm2/s]
               FLUX         = ALPHA_HCN * KKG * CG     

               ! Amount of tagged tracer lost to ocean [kg]
               AMT_LOST     = FLUX * FOCEAN * ACM2 * DTCHEM

               ! Save total HCN lost to ocean for ND09 diag [molec/cm2/s]
               IF ( N == 1 ) THEN
                  OCEAN_HCN = AMT_LOST * XNUMOL_HCN / ( ACM2 * DTCHEM )
               ENDIF

               ! Subtract ocean loss from STT array [kg/box/step]
               STT(I,J,1,N) = MAX( STT(I,J,1,N) - AMT_LOST, 0e+0_fp )

            ENDDO

            !-----------------------------------------------------------
            ! (2) CH3CN ocean uptake
            !-----------------------------------------------------------

            ! Henry's law constant for CH3CN [unitless]
            H = 861.7e+0_fp * EXP( -4100e+0_fp / State_Met%TS(I,J) ) 

            ! SC is Schmidt # for HCN in seawater [unitless]
            SC = B0 + TC * ( B1 + TC * ( B2 + TC * ( B3 )))

            ! KL: conductance for mass transfer in liquid phase
            ! (Wanninkhof 1992), which has units of [cm/h]
            KL = ( 0.222e+0_fp * U * U  + 0.333e+0_fp * U )
     &         * ( SC / 600e+0_fp )**( -0.5e+0_fp )

            ! KG: conductance for mass transfer in gas phase (Asher 1997)
            ! Convert from m/s to cm/h by mutiplying by 360000
            KG = ( 12.4e+0_fp + 763.3e+0_fp * U ) 

            ! KKG: transfer velocity on a gas phase basis (Liss & Slater 1974)
            ! Convert from [cm/h] to [cm/s] by dividing by 3600
            KKG = 2.78e-4_fp * KL * KG / ( KL + H * KG )

            ! Loop over CH3HCN tagged tracers
            DO NN = 1, N_MAX
               
               ! Get CH3CN tagged tracer number
               N              = CH3CN_INDEX(NN)

               ! Bulk concentration of CH3CN in gas phase [kg/cm3]
               CG             = STT(I,J,1,N) / BOXVL

               ! Air-to-sea flux of HCN [kg/cm2/s]
               FLUX           = ALPHA_HCN * KKG * CG     

               ! Amount of tagged tracer lost to ocean [kg]
               AMT_LOST       = FLUX * FOCEAN * ACM2 * DTCHEM

               ! Save total HCN lost to ocean for ND09 diag [molec/cm2/s]
               IF ( N == 2 ) THEN
                  OCEAN_CH3CN = AMT_LOST * XNUMOL_CH3CN / (ACM2*DTCHEM) 
               ENDIF

               ! Subtract ocean loss from STT array [kg/box/step]
               STT(I,J,1,N)   = MAX( STT(I,J,1,N) - AMT_LOST, 0e+0_fp )

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

      SUBROUTINE INIT_HCN_CH3CN( am_I_Root, Input_Opt, RC )
!
!******************************************************************************
!  Subroutine INIT_TAGGED_HCN_CH3CN allocates memory to module arrays.
!  (bmy, 6/29/05)
! 
!  NOTES:
!  This routine is now called from GIGC_Init_Extra in input_mod.F. 
!******************************************************************************
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

      !=================================================================
      ! INIT_TAGGED_CO begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate arrays
      ALLOCATE( HCN_REG_bb( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'HCN_REG_bb' )         

      ALLOCATE( HCN_REG_df( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'HCN_REG_df' )

      ALLOCATE( CH3CN_REG_bb( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'CH3CN_REG_bb' )         

      ALLOCATE( CH3CN_REG_df( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'CH3CN_REG_df' )
      
      ALLOCATE( EMIS_CO_df( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMIS_CO_df' )

      ! Define geographic regions for biomass burning
      CALL DEFINE_BB_REGIONS

      ! Define geographic regions for domestic fossil fuel burning
      CALL DEFINE_DF_REGIONS

!      ! Read domestic fossil fuel emissions
!      CALL READ_EMISSIONS

      ! Index of HCN tracers
      HCN_INDEX(:)   = (/ 1, 3, 4, 5, 6  /)

      ! Index of CH3CN tracers
      CH3CN_INDEX(:) = (/ 2, 7, 8, 9, 10 /)

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

      ! Free the OH pointer
      OH => NULL()

      ! Return to calling program
      END SUBROUTINE CLEANUP_HCN_CH3CN

!------------------------------------------------------------------------------

      ! End of module
      END MODULE HCN_CH3CN_MOD
