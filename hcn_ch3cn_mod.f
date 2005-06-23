! $Id: hcn_ch3cn_mod.f,v 1.1 2005/06/23 19:32:56 bmy Exp $
      MODULE HCN_CH3CN_MOD
!
!******************************************************************************
!  Module HCN_CH3CN_MOD contains variables and routines that are used for the 
!  geographically tagged HCN/CH3CN simulation. (qli, xyp, bmy, 6/23/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) BB_REGION         : Array w/ geographic regions for biomass burning
!  (2 ) DF_REGION         : Array w/ geographic regions for domestic fossilfuel
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
!  GEOS-CHEM modules referenced by hcn_ch3cn_mod.f
!  ============================================================================
!  (1 ) diag_mod.f        : Module w/ GEOS-CHEM diagnostic arrays
!  (2 ) dao_mod.f         : Module w/ arrays for DAO met fields!
!  (3 ) biomass_mod.f     : Module w/ routines to read biomass emissions
!  (4 ) geia_mod          : Module w/ routines to read anthro emissions
!  (5 ) global_oh_mod.f   : Module w/ routines to read 3-D OH field
! 
!  Tagged HCN/CH3CN tracers:
!  ============================================================================
!  (1 ) Total HCN
!  (2 ) HCN from Asian biomass burning
!  (3 ) HCN from elsewhere biomass burning 
!  (4 ) HCN from Asian domestic fossil fuel 
!  (5 ) HCN from elsewhere domestic fossil fuel
!  (6 ) Total CH3CN
!  (7 ) CH3CN from Asian biomass burning
!  (8 ) CH3CN from elsewhere biomass burning 
!  (9 ) CH3CN from Asian domestic fossil fuel 
!  (10) CH3CN from elsewhere domestic fossil fuel
!
!  References:
!  ============================================================================
!
!  NOTES:
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

      ! Arrays
      INTEGER, ALLOCATABLE :: BB_REGION(:,:)
      INTEGER, ALLOCATABLE :: DF_REGION(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS     

!------------------------------------------------------------------------------

      SUBROUTINE DEFINE_BB_REGIONS( REGION )
!
!******************************************************************************
!  Subroutine DEFINE_BB_REGIONS defines the geographic regions 
!  for biomass burning emissions for the tagged HCN/CH3CN simulation. 
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) REGION (INTEGER) : Array of Fossil Fuel CO regions 
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, GET_XMID, GET_YMID

#     include "CMN_SIZE" 

      ! Arguments
      INTEGER, INTENT(OUT) :: REGION(IIPAR,JJPAR)

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

         ! Latitude [degrees]
         Y = GET_YMID( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Longitude [degrees]
            X = GET_XMID( I )

            ! Region #2/#7 -- SE Asian BB (1st sub-box)
            IF      ( ( X >= 72.5 .AND. X < 127.5 )  .AND.
     &                ( Y >=  8.0 .AND. Y <  28.0 ) ) THEN
               REGION(I,J) = 2

            ! Region #2/#7 -- SE Asian BB (2nd sub-box)
            ELSE IF ( ( X >= 72.5 .AND. X < 152.5 )  .AND.
     &                ( Y >= 28.0 .AND. Y <  48.0 ) ) THEN
               REGION(I,J) = 2
  
            ! Region #3 -- African BB
            ELSE IF ( ( X >= -17.5 .and. X < 65.0 )  .and.
     &                ( Y >= -36.5 .and. Y < 36.0 ) ) THEN
               REGION(I,J) = 3

            ! Region #4 -- South America
            ELSE IF ( ( X >= -100.0 .and. X < -30.0 )  .and.
     &                ( Y >= -60.0 .and. Y < 0.0 ) ) THEN
               REGION(I,J) = 4

            ! Region #5 -- BB from elsewhere
            ELSE
               REGION(I,J) = 5
               
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DEFINE_BB_REGIONS

!------------------------------------------------------------------------------

      SUBROUTINE DEFINE_DF_REGIONS( REGION )
!
!******************************************************************************
!  Subroutine DEFINE_DF_REGIONS defines the geographic regions
!  for domestic fossil fuel emissions for the tagged HCN/CH3CN simulation. 
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) REGION (INTEGER) : Array of Fossil Fuel regions 
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID, GET_YMID

#     include "CMN_SIZE"    ! Size parameters 

      ! Arguments
      INTEGER, INTENT(OUT) :: REGION(IIPAR,JJPAR)

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

         ! Latitude [degrees]
         Y = GET_YMID( J )         

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Longitude [degrees]
            X = GET_XMID( I )
         
            ! Region #4/#9 -- Asian DF (1st sub-box)
            IF      ( ( X >= 72.5 .AND. X < 127.5 )  .AND.
     &                ( Y >=  8.0 .AND. Y <  28.0 ) ) THEN
               REGION(I,J) = 6

            ! Region #4/#9 -- Asian DF (2nd sub-box)
            ELSE IF ( ( X >= 72.5 .AND. X < 152.5 )  .AND.
     &                ( Y >= 28.0 .AND. Y <  48.0 ) ) THEN
               REGION(I,J) = 6
   
            ! Region #5/#10 -- DF from elsewhere
            ELSE
               REGION(I,J) = 7
               
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DEFINE_DF_REGIONS

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_HCN_CH3CN( N_TRACERS, STT )
!
!******************************************************************************
!  Subroutine EMISS_HCN_CH3CN reads in CO emissions and scale them to 
!  get HCN/CH3CN emissions for the tagged HCN/CH3CN run.
!
!  Arguments as Input:
!  ============================================================================
!  (1) FIRSTEMISS (LOGICAL) : = T if this is the first call to this routine
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,   ONLY : BURNEMIS, BIOBURN
      USE GEIA_MOD
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      !USE DIAG_MOD,      ONLY : AD10
      USE LOGICAL_MOD,   ONLY : LSPLIT
      USE PBL_MIX_MOD,   ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX
      USE TIME_MOD,      ONLY : GET_TS_CHEM, GET_MONTH, ITS_A_NEW_MONTH
      
#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_O3"        ! SCNR89, TODH, TODB, TODN, EMISTRCO
!#     include "CMN_DIAG"      ! ND10
#     include "comtrid.h"     ! IDBCO

      ! Arguments
      INTEGER, INTENT(IN)    :: N_TRACERS
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: PBL_MAX 
      INTEGER                :: I,      J,      L,      N
      REAL*8                 :: ACM2,   E_CObb, E_COdf, SFAC89 
      REAL*8                 :: DTSRCE, HCN_bb, HCN_df, FRAC

      ! Emission ratios for HCN/CH3CN from biomass burning 
      ! and domestic fossil fuel
      REAL*8,  PARAMETER     :: EHCN_bb   = 0.27d-2
      REAL*8,  PARAMETER     :: EHCN_df   = 1.60d-2
      REAL*8,  PARAMETER     :: ECH3CN_bb = 0.20d-2
      REAL*8,  PARAMETER     :: ECH3CN_df = 0.25d-2

      ! External functions
      REAL*8, EXTERNAL       :: BOXVL

      !=================================================================
      ! EMISS_TAGGED_HCN_CH3CN begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN 

         ! Read time-of-day and day-of-week scale factors for GEIA emissions
         CALL READ_TODX( TODN, TODH, TODB, SCNR89 )

         ! Read domestic fossil fuel CO emissions from GEIA
         CALL READ_GEIA( E_RCO=EMISTRCO  )

         ! Allocate all module arrays
         CALL INIT_TAGGED_HCN_CH3CN

         ! Set first-time flag to false
         FIRST = .FALSE.
      ENDIF

      ! DTSRCE is the number of seconds per emission timestep
      DTSRCE = GET_TS_CHEM() * 60d0

      !=================================================================
      ! Process biomass burning/domestic fossil fuel HCN/CH3CN emissions
      !=================================================================
      CALL BIOBURN( MONTH )

!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I,      J,      L,     N,    ACM2,   E_CObb )
!!$OMP+PRIVATE( SFAC89, E_COdf, IHOUR, FRAC, HCN_bb, HCN_df ) 
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid area in [cm2]
         ACM2 = GET_AREA_CM2( J )

         !-----------------------------------------------------------------
         ! (1) Process biomass burning HCN/CH3CN emissions
         !-----------------------------------------------------------------

         ! Convert [molec CO/cm3/s] to [molec CO/cm2/s]
         E_CObb = BURNEMIS(IDBCO,I,J) * BOXVL(I,J,1) / ACM2

         !! Diagnostic ND10: biomass burning HCN/CH3CN emissions [molec/cm2/s]
         !IF ( ND10 > 0 ) THEN
         !   AD10(I,J,1) = AD10(I,J,1) + EHCN_BB * E_CObb
         !ENDIF

         ! Convert [molec CO/cm2/s] to [mole/grid box]: 1/6.022d23 = 1.66d-24
         E_CObb = E_CObb * 1.66d-24 * ACM2 * DTSRCE

         !-----------------------------------------------------------------
         ! (2) Process domestic fossil fuel HCN/CH3CN emissions
         !-----------------------------------------------------------------

         ! SFAC89 is the Weekday/Saturday/Sunday scale factor
         SFAC89 = SCNR89( 2, GET_DAY_INDEX( NTAU ) ) 

         ! E_COdf is DF CO emissions in [molec CO/cm2/s]
         ! Scale E_COdf by the day-of-week scale factor SFAC89
         E_COdf = EMISTRCO(I,J) * SFAC89

         ! Scale E_COdf by the time-of-day scale factor TODH
         ! IHOUR is the index for the time-of-day scale factor TODH
         IHOUR = GET_IHOUR( I, TOFDAY, NDYN, DISIZE )
         E_COdf  = E_COdf * TODH(IHOUR)

         ! Enhance E_COdf by 18.5% to account for oxidation 
         ! from anthropogenic VOC's (bnd, bmy, 6/8/01)
         E_COdf = E_COdf * 1.185d0
            
         ! Get domestic fossil fuel region #
         N = DF_REGION(I,J)

         ! To achieve the best fit to the observed HCN-CH3CN-CO correlations 
         ! in the boundary layer, we have to double the residential coal 
         ! burning source from Asia. This leads us to reduce the residential 
         ! coal burning source from the rest of the world by a factor of eight
         ! to achieve a best fit to the observed vertical distributions of HCN
         ! and CH3CN. (xyp, 6/22/05)
         IF ( N == 6 ) THEN
            E_COdf = E_COdf * 2.1d0   ! Asian domestic fossil fuel
         ELSE 
            E_COdf = E_COdf / 8.0d0   ! Elsewhere domestic fossil fuel
         ENDIF

         !! ND10: domestic fossil fuel HCN/CH3CN emissions [molec/cm2/s]
         !IF ( ND10 > 0 ) THEN
         !   AD10(I,J,2) = AD10(I,J,2) + EHCN_DF   * E_COdf
         !ENDIF

         ! Convert [molec CO/cm2/s] to [mole/grid box]: 1/6.022d23 = 1.66d-24
         E_COdf = E_COdf * 1.66d-24 * ACM2 * DTSRCE       

         !-----------------------------------------------------------------
         ! (3) Partition emissions throughout the boundary layer
         !-----------------------------------------------------------------

         ! Loop up to the highest PBL level
         DO L = 1, PBL_MAX
            
            ! Fraction of the PBL occupied by this layer
            FRAC            = GET_FRAC_OF_PBL( I, J, L )

            ! HCN biomass burning emissions
            HCN_bb          = FRAC * MHCN * EHCN_BB * E_CObb
         
            ! HCN domestic fossil fuel emissions
            HCN_df          = FRAC * MHCN * EHCN_DF * E_COdf

            ! Total HCN/CH3CN emissions (BB+DF)
            STT(I,J,L,1)    = STT(I,J,L,1) + HCN_bb + HCN_df

            ! If we are using tagged tracers ...
            IF ( LSPLIT ) THEN
               
               ! Add emissions into tagged biomass tracers
               N            = BB_REGION(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + HCN_bb

               ! Add emissions into tagged domestic fossil fuel tracers
               N            = DF_REGION(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + HCN_df

            ENDIF
         ENDDO
      ENDDO
      ENDDO
!!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMISS_TAGGED_HCN_CH3CN

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_HCN_CH3CN( N_TRACERS, STT )
!
!******************************************************************************
! HCN/CH3CN loss due to reaction with OH and ocean uptake.
!
!  Arguments as Input:
!  ============================================================================
!  (1) FIRSTCHEM (LOGICAL) : = T if this is the first call to this routine
!******************************************************************************
! 
      ! References to F90 modules
      USE DAO_MOD,       ONLY : ALBD, TS, U10M, V10M
      USE DIAG_MOD,      ONLY : AD09, AD10
      USE GLOBAL_OH_MOD, ONLY : OH, GET_GLOBAL_OH
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE TIME_MOD,      ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters 
#     include "CMN_DIAG"     ! ND10
#     include "CMN_DEP"      ! FRCLND 
#     include "CMN_SETUP"    ! LSPLIT

      ! Arguments
      ! Arguments
      INTEGER, INTENT(IN)    :: N_TRACERS
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
      
      ! Local variables
      INTEGER, SAVE          :: LAST_MONTH=-1

      ! Local variables
      INTEGER                :: I,    J,  K,    L,     N,      N_MAX
      REAL*8                 :: K0,   K1, TMP,  KRATE, DTCHEM, KTMP
      REAL*8                 :: H,    U,  TC,   SC,    KL,     KG
      REAL*8                 :: KKG,  CL, SR,   CG,    FLUX,   FOCEAN
      REAL*8                 :: ACM2, AMT_LOST, OCEAN_HCN, OCEAN_CH3CN

      ! Undersaturation ratios for HCN/CH3CN in seawater
      REAL*8, PARAMETER      :: ALPHA_HCN   = 0.21d0
      REAL*8, PARAMETER      :: ALPHA_CH3CN = 0.12d0

      ! Coefficients for fitting the Schmdit number for HCN in seawater
      REAL*8, PARAMETER      :: A0 = 2008.917d0
      REAL*8, PARAMETER      :: A1 =  -83.235d0
      REAL*8, PARAMETER      :: A2 =    1.348d0
      REAL*8, PARAMETER      :: A3 =   -0.009d0

      ! Coefficients for fitting the Schmdit number for CH3CN in seawater
      REAL*8, PARAMETER      :: B0 = 2745.722d0
      REAL*8, PARAMETER      :: B1 = -113.763d0
      REAL*8, PARAMETER      :: B2 =    1.843d0
      REAL*8, PARAMETER      :: B3 =   -0.012d0

      ! External functions
      REAL*8, EXTERNAL       :: BOXVL
      
      !=================================================================
      ! CHEM_HCN_CH3CN begins here! 
      !=================================================================

      ! First-time initialization (if not already done)
      IF ( FIRST ) THEN
         CALL INIT_HCN_CH3N
         FIRST = .FALSE.
      ENDIF

      ! Read offline OH fields once per month
      IF ( ITS_A_NEW_MONTH() ) THEN 
         CALL GET_GLOBAL_OH( GET_MONTH() )
      ENDIF
     
      ! Compute number of tracers to process
      IF ( LSPLIT ) THEN
         N_MAX = N_TRACERS
      ELSE
         N_MAX = 1
      ENDIF

      !=================================================================
      ! Do HCN and CH3CN chemistry
      !=================================================================

      ! Chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over grid boxes
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, K0, K1, TMP, KTMP, KRATE, N )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !------------------------------------------------------------------
         ! (1) HCN loss via reaction with OH
         !------------------------------------------------------------------

         K0    = 7.4d-33 
         K1    = 9.0d-15 * ( T(I,J,L) / 300d0 ) ** 3.2d0

         ! AD: air mass in kg
         TMP   = K0 / K1 * AD(I,J,L) * XNUMOL_AIR / BOXVL(I,J,L)

         ! K: [cm3/molec/s]
         KTMP  = K1 * TMP / ( 1d0 + TMP )      
     &         * EXP ( -0.511d0 / ( 1d0 + LOG10( TMP ) ** 2d0 ) )

         KRATE = KTMP * OH(I,J,L) * DTCHEM

         !--------------------------------------------------------------
         ! (2) Subtract lost tracer from STT array
         !--------------------------------------------------------------
         DO N = 1, N_MAX 

            ! Compute the amount of tracer that is lost to OH
            AMT_LOST     = KRATE * STT(I,J,L,N)

            ! Remove lost tracer from STT array (avoid negatives!)
            STT(I,J,L,N) = MAX( STT(I,J,L,N) - AMT_LOST, 0d0 )
            
            ! ND09 diagnostic: HCN/CH3CN loss via OH [kg]
            IF ( ND09 > 0 ) THEN
               AD09(I,J,L,N) = AD09(I,J,L,N) + AMT_LOST
            ENDIF
         ENDDO
      ENDDO
      ENDDO
      ENDDO
!!$OMP END PARALLEL DO

      !=================================================================
      ! HCN and CH3CN ocean uptake
      !=================================================================

      ! Loop over grid boxes
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I,    J,         ACM2, FOCEAN,   U,   TC  )
!!$OMP+PRIVATE( H,    SC,        KL,   KG,       KKG, CG  )
!!$OMP+PRIVATE( FLUX, OCEAN_HCN, N,    AMT_LOST, CL       )
!!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box area in [cm2]
         ACM2            = GET_AREA_CM2( J ) 
            
         ! Fraction of a grid box that is ocean
         FOCEAN          = 1d0 - FRCLND(I,J) 

         !--------------------------------------------------------------
         ! Only compute ocean sink if there is more than 50% ocean
         ! in the grid box, and if it is not ice (albedo > 0.4)
         ! (mje, rvm, bmy, 11/26/01)
         !--------------------------------------------------------------
         IF ( FOCEAN > 0.5d0 .AND. ALBD(I,J) <= 0.4d0 ) THEN

            ! Wind speed [m/s] at 10m above the surface 
            U            = SQRT( U10M(I,J)**2 + V10M(I,J)**2 )

            ! Surface temperature [C]
            TC           = TS(I,J) - 273.15d0  

            ! Dimensionless Henry's law constants for HCN
            H            = 7.93d4 * EXP( -5000d0 / TS(I,J) ) 
            
            ! SC is Schmidt # for HCN in seawater [unitless]
            SC           = A0 + TC * ( A1 + TC * ( A2 + TC * ( A3 )))

            ! KL: conductance for mass transfer in liquid phase 
            ! (Nightingale 2000), which has unit of [cm/h]
            KL           = ( 0.222d0 * U * U + 0.333d0 * U ) 
     &                   * ( SC / 600d0 )**( -0.5d0 )  

            ! KG: conductance for mass transfer in gas phase (Asher 1997)
            ! Convert from m/s to cm/h by multiplying 360000
            KG           = ( 15.3d0 + 940.6d0 * U ) 

            ! KKG: transfer velocity on a gas phase basis (Liss & Slater 1974)
            ! Convert from [cm/h] to [cm/s] by dividing 3600
            KKG          = 2.78d-4 * KL * KG / ( KL + H * KG )

            ! CG: bulk concentration of HCN in gas phase [kg/cm3]
            CG           = STT(I,J,1,1) / BOXVL(I,J,1)

            ! FLUX: air-to-sea flux of HCN [kg/cm2/s]
            FLUX         = ALPHA_HCN * KKG * CG     
            CL           = ( 1d0 - ALPHA_HCN ) * CG / H

            ! Amount of HCN lost to the ocean [kg]
            OCEAN_HCN    = FLUX * FOCEAN * ACM2 * DTCHEM 

            ! Subtract ocean loss from STT array [kg/box/step]
            STT(I,J,1,1) = MAX( STT(I,J,1,1) - OCEAN_HCN, 0d0 )

            ! If there is more than one tracer
            IF ( LSPLIT ) THEN

               ! Subtract ocean loss for tagged tracers
               DO N = 2, N_TRACERS

                  ! FLUX: air-to-sea flux of tagged tracers [kg/cm2/s]
                  FLUX         = ALPHA_HCN    * KKG * 
     *                           STT(I,J,1,N) / BOXVL(I,J,1)    

                  ! Amount of tagged tracer lost to ocean [kg]
                  AMT_LOST     = FLUX * FOCEAN * ACM2 * DTCHEM

                  ! Remove lost tagged tracer from STT array [kg]
                  STT(I,J,1,N) = MAX( STT(I,J,1,N) - AMT_LOST, 0d0 )
               ENDDO
            ENDIF

         !--------------------------------------------------------------
         ! If there is less than 50% water in the grid box, or  
         ! if there is ice on the ocean, then zero the ocean sink
         !--------------------------------------------------------------
         ELSE

            ! Set to zero
            OCEAN_HCN = 0d0
            CL        = 0d0

         ENDIF

         !--------------------------------------------------------------
         ! ND10 diag: Save HCN/CH3CN ocean uptake in [molec/cm2/s]
         !--------------------------------------------------------------
         IF ( ND10 > 0 ) THEN
            AD10(I,J,3) = AD10(I,J,3) + 
     &                    ( OCEAN_HCN * XNUMOL_HCN / ( ACM2 * DTCHEM ) )

            ! CL in kg/cm3
            ! how to deal with SR when less than 50% water in grid box???
            AD10(I,J,4) = AD10(I,J,4) + CL
         ENDIF
      ENDDO
      ENDDO
!!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_HCN_CH3CN

!------------------------------------------------------------------------------

      SUBROUTINE INIT_HCN_CH3CN
!
!******************************************************************************
!  Subroutine INIT_TAGGED_HCN_CH3CN allocates memory to module arrays.
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      LOGICAL, SAVE      :: IS_INIT = .FALSE.
      INTEGER            :: AS
      
      !=================================================================
      ! INIT_TAGGED_CO begins here!
      !=================================================================

      ! 
      IF ( IS_INIT ) RETURN

      ! Allocate BB_REGION -- array for biomass burning regions
      ALLOCATE( BB_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BB_REGION' )         

      ! Allocate DF_REGION -- array for fossil fuel regions
      ALLOCATE( DF_REGION( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DF_REGION' )
      
      ! Define geographic regions for biomass burning
      CALL DEFINE_BB_REGIONS( BB_REGION )

      ! Define geographic regions for domestic fossil fuel burning
      CALL DEFINE_DF_REGIONS( DF_REGION )      

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
      IF ( ALLOCATED( BB_REGION ) ) DEALLOCATE( BB_REGION )
      IF ( ALLOCATED( DF_REGION ) ) DEALLOCATE( DF_REGION )

      ! Return to calling program
      END SUBROUTINE CLEANUP_HCN_CH3CN

!------------------------------------------------------------------------------

      ! End of module
      END MODULE HCN_CH3CN_MOD
