! $Id: ocean_mercury_mod.f,v 1.3 2005/03/29 15:52:43 bmy Exp $
      MODULE OCEAN_MERCURY_MOD
!
!******************************************************************************
!  Module OCEAN_MERCURY_MOD contains variables and routines needed to compute
!  the oceanic flux of mercury.  Original code by Sarah Strode at UWA/Seattle.
!  (sas, bmy, 1/21/05, 2/24/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) DD_Hg2 (REAL*8) : Array for Hg(II) dry deposited to ocean   [kg      ]
!  (2 ) Hg0aq  (REAL*8) : Array for ocean mass of Hg(0)             [kg      ]
!  (3 ) Hg2aq  (REAL*8) : Array for ocean mass of Hg(II)            [kg      ]
!  (4 ) MLD    (REAL*8) : Array for monthly-mean mixed layer depths [cm      ]
!  (5 ) NPP    (REAL*8) : Array for monthly-mean net primary prod.  [unitless]
!  (6 ) RAD    (REAL*8) : Array for monthly-mean solar radiation    [W/m2    ]
!  (7 ) WD_Hg2 (REAL*8) : Array for Hg(II) dry deposited to ocean   [kg      ]
!
!  Module Routines:
!  ============================================================================
!  (1 ) ADD_Hg2_DD            : Archives Hg2 lost to drydep in DD_HG2
!  (2 ) ADD_Hg2_WD            : Archives Hg2 lost to wetdep in WD_HG2
!  (3 ) OCEAN_MERCURY_FLUX    : Routine to compute flux of oceanic mercury
!  (4 ) OCEAN_MERCURY_READ    : Routine to read MLD, NPP, RADSWG data from disk
!  (5 ) INIT_OCEAN_MERCURY    : Allocates and zeroes module arrays
!  (6 ) CLEANUP_OCEAN_MERCURY : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by ocean_mercury_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f           : Module w/ routines for binary pch file I/O
!  (2 ) dao_mod.f             : Module w/ arrays for DAO met fields
!  (3 ) diag03_mod.f          : Module w/ ND03 diagnostic arrays 
!  (2 ) file_mod.f            : Module w/ file unit numbers and error checks
!  (9 ) grid_mod.f            : Module w/ horizontal grid information
!  (10) logical_mod.f         : Module w/ GEOS-CHEM logical switches
!  (11) pressure_mod.f        : Module w/ routines to compute P(I,J,L)
!  (12) time_mod.f            : Module w/ routines to compute date & time
!  (13) tracer_mod.f          : Module w/ GEOS-CHEM tracer array STT etc.
!  (14) tracerid_mod.f        : Module w/ pointers to tracers & emissions
!  (15) transfer_mod.f        : Module w/ routines to cast & resize arrays
!
!  References:
!  ============================================================================
!  (1 ) Xu et al (1999). Formulation of bi-directional atmosphere-surface
!        exchanges of elemental mercury.  Atmospheric Environment 
!        33, 4345-4355.
!  (2 ) Nightingale et al (2000).  In situ evaluation of air-sea gas exchange
!        parameterizations using novel conservative and volatile tracers.  
!        Global Biogeochemical Cycles, 14, 373-387.
!  (3 ) Lin and Tau (2003).  A numerical modelling study on regional mercury 
!        budget for eastern North America.  Atmos. Chem. Phys. Discuss., 
!        3, 983-1015.  And other references therein.
!  (4 ) Poissant et al (2000).  Mercury water-air exchange over the upper St.
!        Lawrence River and Lake Ontario.  Environ. Sci. Technol., 34, 
!        3069-3078. And other references therein.
!  (5 ) Wangberg et al. (2001).  Estimates of air-sea exchange of mercury in 
!        the Baltic Sea.  Atmospheric Environment 35, 5477-5484.
!  (6 ) Clever, Johnson and Derrick (1985).  The Solubility of Mercury and some
!        sparingly soluble mercury salts in water and aqueous electrolyte
!        solutions.  J. Phys. Chem. Ref. Data, Vol. 14, No. 3, 1985.
!
!  Nomenclature: 
!  ============================================================================
!  (1 ) Hg(0)  a.k.a. Hg0 : Elemental mercury
!  (2 ) Hg(II) a.k.a. Hg2 : Divalent  mercury
!
!  NOTES:
!  (1 ) Modified ocean flux w/ Sarah's new Ks value (sas, bmy, 2/24/05)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "ocean_mercury_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: ADD_Hg2_DD
      PUBLIC :: ADD_Hg2_WD
      PUBLIC :: INIT_OCEAN_MERCURY
      PUBLIC :: CLEANUP_OCEAN_MERCURY
      PUBLIC :: OCEAN_MERCURY_FLUX

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      CHARACTER(LEN=255)  :: Hg0_Hg2_FILE = 'ocean_Hg0_Hg2.YYYYMMDDhh' 

      ! Arrays
      REAL*8, ALLOCATABLE :: DD_Hg2(:,:)
      REAL*8, ALLOCATABLE :: Hg0aq(:,:)
      REAL*8, ALLOCATABLE :: Hg2aq(:,:)
      REAL*8, ALLOCATABLE :: MLD(:,:)
      REAL*8, ALLOCATABLE :: NPP(:,:)
      REAL*8, ALLOCATABLE :: RAD(:,:)
      REAL*8, ALLOCATABLE :: WD_Hg2(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE ADD_Hg2_DD( I, J, DRY_Hg2 )
!
!******************************************************************************
!  Subroutine ADD_Hg2_WD computes the amount of Hg(II) dry deposited out of
!  the atmosphere into the column array DD_Hg2. (sas, bmy, 1/19/05)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) I       (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J       (INTEGER) : GEOS-CHEM latitude index
!  (3 ) DRY_Hg2 (REAL*8 ) : Hg(II) dry deposited out of the atmosphere [kg]
!
!  NOTES:
!******************************************************************************
!
      ! Arguments as input
      INTEGER, INTENT(IN) :: I, J
      REAL*8,  INTENT(IN) :: DRY_Hg2
 
      !=================================================================
      ! ADD_Hg2_DD begins here!
      !=================================================================
      DD_Hg2(I,J) = DD_Hg2(I,J) + DRY_Hg2

      ! Return to calling program
      END SUBROUTINE ADD_Hg2_DD

!------------------------------------------------------------------------------

      SUBROUTINE ADD_Hg2_WD( I, J, WET_Hg2 )
!
!******************************************************************************
!  Subroutine ADD_Hg2_WD computes the amount of Hg(II) wet scavenged out of
!  the atmosphere into the column array WD_Hg2. (sas, bmy, 1/19/05)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) I       (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J       (INTEGER) : GEOS-CHEM latitude index
!  (3 ) WET_Hg2 (REAL*8 ) : Hg(II) scavenged out of the atmosphere
!
!  NOTES:
!******************************************************************************
!
      ! Arguments as input
      INTEGER, INTENT(IN) :: I, J
      REAL*8,  INTENT(IN) :: WET_Hg2
 
      !=================================================================
      ! ADD_Hg2_WD begins here!
      !=================================================================
      WD_Hg2(I,J) = WD_Hg2(I,J) + WET_Hg2

      ! Return to calling program
      END SUBROUTINE ADD_Hg2_WD

!------------------------------------------------------------------------------

      SUBROUTINE OCEAN_MERCURY_FLUX( FLUX )
!
!******************************************************************************
!  Subroutine OCEAN_MERCURY_FLUX calculates annual monthly emissions from 
!  the ocean in [kg/s].  (sas, bmy, 1/19/05, 2/24/05)
!
!  Arguments as Output
!  ============================================================================
!  (1 ) FLUX (REAL*8) : Flux of Hg(0) from the ocean [kg/s]
!
!  NOTES:
!  (1 ) Change Ks to make ocean flux for 2001 = 2.03e6 kg/year.
!        (sas, bmy, 2/24/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AIRVOL, ALBD, TS, RADSWG
      USE DIAG03_MOD,   ONLY : AD03, ND03
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE TIME_MOD,     ONLY : GET_TS_EMIS, GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHg0, IDTHg2

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DEP"      ! FRCLND

      ! Arguments 
      REAL*8,  INTENT(OUT)  :: FLUX(IIPAR,JJPAR)

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: I,       J,        L,        THISMONTH
      REAL*8                :: AREA_M2, ARG1,     CHg0,     CHg0aq
      REAL*8                :: CHg2aq,  DTSRCE,   FRAC_L,   FRAC_O  
      REAL*8                :: H,       Hg2_LOST, Hg2_SUNK, K1
      REAL*8                :: Ks,      Kw,       MLDCM,    MLDCMC  
      REAL*8                :: TC,      TK,       TOTDEP,   Sc
      REAL*8                :: ScCO2,   SUNKHg,   USQ

      ! Conversion factor from [cm/h * ng/L] --> [kg/m2/s]
      REAL*8,  PARAMETER    :: TO_KGM2S = 1.0D-11 / 3600D0 

      ! External functions
      REAL*8,  EXTERNAL     :: SFCWINDSQR 
      
      !=================================================================
      ! OCEAN_MERCURY_FLUX begins here!
      !=================================================================

      !---------------------------
      ! Read monthly data
      !---------------------------
      IF ( ITS_A_NEW_MONTH() ) THEN

         ! Get current month
         THISMONTH = GET_MONTH()

         ! Get monthly MLD, NPP, radiation
         CALL OCEAN_MERCURY_READ( THISMONTH )

      ENDIF

      !=================================================================
      ! Compute flux of Hg0 from the ocean (notes by Sarah Strode):
      !
      ! Net flux is given by the equation:
      !
      !    F = Kw * ( Caq - Cg/H )                   [Xu et al, 1999]
      ! 
      ! Kw is the exchange parameter (piston velocity) [cm/h] given by:
      !
      !    Kw = 0.25 * u^2 / SQRT( Sc / ScCO2 )    [Nightingale, 2000]
      !
      ! u^2 is the square of the wind speed (10m above ground) [m/s].
      !
      ! Sc is the Schmidt number [unitless] for Hg(0):
      !
      !    Sc = nu/D = ( 0.017*exp(-0.025T) ) / ( 6.0*10^-7*T + 10^-5 )
      !       with T in deg. C         
      !       [Lin and Tao, 2003 and Poissant et al., 2000]
      !
      ! Caq = 1.5 pM is the surface water concentration       
      !    [Lamborg et al., 2002]
      !
      ! Convert Caq to ng/L via 1.5 * 10^-12 * atomicWeight(Hg) * 10^9
      !
      ! Cg is the gas-phase concentration
      !
      ! H is the dimensionless Henry coefficient for elemental mercury
      !
      !    H = exp(4633.3/T - 14.52) where T is sea temp. in Kelvin
      !       [Wangberg et al., 1999 and Clever et al, 1985]
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Define sinking term Ks [unitless] to get 2.03e6 kg/year 
      ! (sas, bmy, 3/8/05)
#if   defined( GEOS_4 )
      Ks     = 7.3d-7 * DTSRCE
#else
      Ks     = 9.2d-7 * DTSRCE
#endif

      ! We are at the surface level
      L      = 1

      ! Loop over latitudes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,        AREA_M2,  MLDCM,  MLDCMC, Kw     )
!$OMP+PRIVATE( K1,   SUNKHg,   TK,       TC,     FRAC_L, FRAC_O )
!$OMP+PRIVATE( H,    Sc,       ScCO2,    Usq,    CHg0,   TOTDEP )
!$OMP+PRIVATE( ARG1, Hg2_LOST, Hg2_SUNK, CHg0aq                 )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Grid box surface area [m2]
         AREA_M2 = GET_AREA_M2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Initialize values
            MLDCM     = MLD(I,J)
            MLDCMC    = MLDCM
            FLUX(I,J) = 0d0
            Kw        = 0d0
            K1        = 0d0
            SUNKHg    = 0d0
            TK        = 0d0
            TC        = 0d0

            ! Get fractions of land and ocean in the grid box [unitless]
            FRAC_L    = FRCLND(I,J)
            FRAC_O    = 1d0 - FRAC_L

            !===========================================================
            ! Make sure we are in an ocean box
            !===========================================================
            IF ( ALBD(I,J) <= 0.4d0 .and. FRAC_L < 0.8d0 ) THEN

               ! Calculate K1 based on NPP & RAD
               !--------------------------------------------------------
               ! Uncomment if using monthly-mean RADSWG fields
               !K1        = 5.3d-11 * DTSRCE * NPP(I,J) * RAD(I,J) 
               !--------------------------------------------------------
               ! Uncomment if using daily RADSWG fields
               K1        = 5.3d-11 * DTSRCE * NPP(I,J) * RADSWG(I,J) 
               !--------------------------------------------------------
               
               ! Surface air temperature in both [K] and [C]
               ! (Use as surrogate for SST, cap at freezing point)
               TK        = MAX( TS(I,J), 273.15d0 )
               TC        = TK - 273.15d0

               ! Henry's law constant [unitless]  
               H         = EXP( 4633.3d0 / TK - 14.52d0 )

               ! Schmidt # for Hg [unitless] 
               Sc        = ( 0.017d0 * EXP( -0.025d0 * TC ) ) / 
     &                     ( 6.0d-7  * TC + 1d-5            )
            
               ! Schmidt # of CO2 [unitless]
               ScCO2     = 644.7d0 + TC * ( -6.16d0 + TC * ( 0.11d0 ) ) 

               ! Square of surface (actually 10m) wind speed [m2/s2]
               Usq       = SFCWINDSQR(I,J)

               ! Piston velocity [cm/h], now from Nightingale
               Kw        = ( 0.25d0 * Usq ) / SQRT( Sc / ScCO2 )

               ! Gas phase Hg(0) concentration: convert [kg] -> [ng/L]
               CHg0      = STT(I,J,L,IDTHg0) * 1.0D9 / AIRVOL(I,J,L)        

               ! Cap mixed layer depth for Hg2 reduction
               MLDCMC    = MIN( MLDCMC, 2500d0 )

               ! Cap MLD for mixing at 100m
               MLDCM     = MIN( MLDCM, 1.0d4 )

               !========================================================
               ! Calculate oceanic Hg(0) and Hg(II)
               ! A certain fraction of Hg(II) will convert to Hg(0)
               !========================================================
               IF ( MLDCM > 0.99d0 ) THEN

                  !-----------------------------------------------------
                  ! Compute ocean mass and concentration of Hg(II)
                  !-----------------------------------------------------

                  ! Total Hg(II) deposited on ocean surface [kg]
                  TOTDEP       = ( WD_Hg2(I,J) + DD_Hg2(I,J) ) * FRAC_O

                  ! Add deposited Hg(II) to the Hg(II) ocean mass [kg]
                  Hg2aq(I,J)   = Hg2aq(I,J) + TOTDEP

                  ! Pre-compute argument for exponential
                  ARG1         = K1 * MLDCMC / MLDCM

                  ! Amount of Hg(II) ocean mass --> Hg(0) ocean mass
                  Hg2_LOST     = Hg2aq(I,J) * ( 1d0 - EXP(    -ARG1 ) ) 
                  
                  ! Amount of Hg(II) that sinks in the ocean [kg]
                  Hg2_SUNK     = Hg2aq(I,J) * ( 1d0 - EXP( -KS-ARG1 ) )
                  
                  ! Archive Hg(II) sinking loss for ND03 [kg/m2/s]
                  SUNKHg       = ( Ks / DTSRCE ) * Hg2aq(I,J) / AREA_M2

                  ! Error check
                  IF ( Hg2aq(I,J) < Hg2_SUNK ) THEN
                     WRITE(6, *) 'got neg 2',K1,Hg2_SUNK,Hg2aq(I,J)
                  ENDIF
                  
                  ! Hg(II) ocean mass after sinking and conversion [kg]
                  Hg2aq(I,J)   = Hg2aq(I,J) * EXP( -KS-ARG1 )

                  ! Hg(II) ocean concentration [ng/L]
                  CHg2aq       = ( Hg2aq(I,J) * 1d11   ) /
     &                           ( AREA_M2    * FRAC_O ) / MLDCM
                 
                  !-----------------------------------------------------
                  ! Compute ocean mass and concentration of Hg(0)
                  !-----------------------------------------------------

                  ! Add converted Hg(II) to the ocean mass of Hg(0) [kg]
                  Hg0aq(I,J)   = Hg0aq(I,J) + Hg2_LOST

                  ! Concentration of Hg(0) in the ocean [ng/L]
                  CHg0aq       = ( Hg0aq(I,J) * 1d11   ) / 
     &                           ( AREA_M2    * FRAC_O ) / MLDCM 

                  !-----------------------------------------------------
                  ! Compute flux of Hg(0) from the ocean to the air
                  !-----------------------------------------------------

                  ! Compute ocean flux of Hg0 [cm/h*ng/L]
                  FLUX(I,J)    = Kw *( CHg0aq - ( H * CHg0 ) ) 
                
                  ! Convert [cm/h*ng/L] --> [kg/m2/s] --> [kg/s]
                  ! Also account for ocean fraction of grid box
                  FLUX(I,J)    = FLUX(I,J) * TO_KGM2S * AREA_M2 * FRAC_O

                  ! Cap the flux w/ the available Hg(0) ocean mass
                  IF ( FLUX(I,J) * DTSRCE > Hg0aq(I,J) ) THEN
                     FLUX(I,J) = Hg0aq(I,J) / DTSRCE
                  ENDIF
               
                  ! Remove amt of Hg(0) that is leaving the ocean [kg]
                  Hg0aq(I,J)  = Hg0aq(I,J) - ( FLUX(I,J) * DTSRCE )

               ELSE

                  !-----------------------------------------------------
                  ! Grid box is mostly land; set Hg(0) flux to zero
                  !-----------------------------------------------------
                  FLUX(I,J)   = 0d0

               ENDIF 

               !========================================================
               ! ND03 diagnostics: archive aqueous Hg0 & Hg2 mass [kg]
               ! 2=Hg0 ocean mass; 7=Hg2 ocean mass; 8=Sunk Hg2; 10=Kw
               !========================================================
               IF ( ND03 > 0 ) THEN
                  AD03(I,J,2)  = AD03(I,J,2)  + Hg0aq(I,J)
                  AD03(I,J,7)  = AD03(I,J,7)  + Hg2aq(I,J)
                  AD03(I,J,8)  = AD03(I,J,8)  + SUNKHg
                  AD03(I,J,10) = AD03(I,J,10) + Kw
               ENDIF

            ENDIF

            ! Zero amts of deposited Hg2 for next timestep [kg]
            DD_Hg2(I,J) = 0d0
            WD_Hg2(I,J) = 0d0

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      !WRITE( 6, * ) 'Total Oceanic Hg Emission: ', SUM(FLUX)*DTSRCE

      ! Return to calling program
      END SUBROUTINE OCEAN_MERCURY_FLUX

!------------------------------------------------------------------------------

      SUBROUTINE OCEAN_MERCURY_READ( THISMONTH )
!
!******************************************************************************
!  Subroutine OCEAN_MERCURY_READ reads in the mixed layer depth, net primary 
!  productivity and radiation climatology for each month.  This is needed for
!  the ocean flux computation. (sas, bmy, 1/20/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Month to read fields (1-12)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: THISMONTH

      ! Local Variables
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: XTAU
      CHARACTER(LEN=255)     :: FILENAME
     
      !=================================================================
      ! MERCURY_READMONTH begins here!
      !=================================================================
     
!------------------------------------------------------------------------------
! Uncomment if using monthly mean RADSWG fields -- we need
! to create these files later on (bmy, 1/20/05)
!      !------------------------------
!      ! Radiation climatology [W/m2]
!      !------------------------------
!
!      ! Radiation file name
!      FILENAME = TRIM( DATA_DIR )              // 
!     &           'mercury_200501/radswg.geos.' // GET_RES_EXT()
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
!
!      ! TAU0 value (uses baseline year 1985)
!      XTAU = GET_TAU0( THISMONTH, 1, 1985 )
!
!      ! Read from disk
!      CALL READ_BPCH2( FILENAME, 'DAO-FLDS',    2,  
!     &                 XTAU,      IGLOB,        JGLOB,     
!     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
!
!      ! Resize and cast to REAL*8
!      CALL TRANSFER_2D( ARRAY(:,:,1), RAD )
!------------------------------------------------------------------------------
      
      !------------------------------
      ! Mixed layer depth [cm]
      !------------------------------

      ! MLD file name
      FILENAME = TRIM( DATA_DIR )           // 
     &           'mercury_200501/mld.geos.' // GET_RES_EXT()      

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MERCURY_READMONTH: Reading ', a )  

      ! TAU0 value (uses year 1985)
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read from disk; original units are [m]
      CALL READ_BPCH2( FILENAME, 'BXHGHT-$',    5,  
     &                 XTAU,      IGLOB,        JGLOB,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), MLD )

      ! Convert [m] to [cm]
      MLD = MLD * 100d0

      !------------------------------
      ! Net primary productivity 
      !------------------------------
 
      ! NPP file name
      FILENAME = TRIM( DATA_DIR )                 // 
     &           'mercury_200501/modis_npp.geos.' // GET_RES_EXT() 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU0 values (uses year 2003)
      XTAU = GET_TAU0( THISMONTH, 1, 2003 )
 
      ! Read data
      CALL READ_BPCH2( FILENAME, 'GLOB-NPP',    1,  
     &                 XTAU,      IGLOB,        JGLOB,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), NPP )
  
      ! Return to calling program
      END SUBROUTINE OCEAN_MERCURY_READ

!------------------------------------------------------------------------------

      SUBROUTINE INIT_OCEAN_MERCURY
!
!******************************************************************************
!  Subroutine INIT_OCEAN_MERCURY allocates and zeroes module arrays.  The
!  initial masses of oceanic Hg(0) and Hg(II) are also read from disk.
!  (sas, bmy, 1/19/05)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE TIME_MOD,     ONLY : EXPAND_DATE, GET_NYMD, GET_NHMS, GET_TAU
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      INTEGER               :: AS, NYMD, NHMS
      REAL*4                :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                :: XTAU
      CHARACTER(LEN=255)    :: FILENAME

      !=================================================================
      ! INIT_OCEAN_MERCURY begins here!
      !=================================================================
      ALLOCATE( DD_Hg2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DD_Hg2' )
      DD_Hg2 = 0d0

      ALLOCATE( Hg0aq( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Hg0aq' )
      Hg0aq = 0d0

      ALLOCATE( Hg2aq( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Hg2aq' )
      Hg2aq = 0d0

      ALLOCATE( MLD( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MLD' )
      MLD = 0d0

      ALLOCATE( NPP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NPP' )
      NPP = 0d0

      !-----------------------------------------------------------------
      ! Uncomment if using monthly-mean RADSWG fields
      !ALLOCATE( RAD( IIPAR, JJPAR ), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'RAD' )
      !RAD = 0d0
      !-----------------------------------------------------------------

      ALLOCATE( WD_Hg2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WD_Hg2' )
      WD_Hg2 = 0d0

      !=================================================================
      ! Read initial ocean mass of Hg(0) and Hg2
      !=================================================================

      ! Get time values
      NYMD = GET_NYMD()
      NHMS = GET_NHMS()
      XTAU = GET_TAU()

      ! Expand date in filename
      FILENAME = Hg0_Hg2_FILE
      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

      !----------------
      ! Oceanic Hg(0) 
      !----------------

      ! Read from disk
      CALL READ_BPCH2( FILENAME, 'HG-SRCE',     2,  
     &                 XTAU,      IGLOB,        JGLOB,     
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), Hg0aq )

      !----------------
      ! Oceanic Hg(II) 
      !----------------

      ! Read from disk
      CALL READ_BPCH2( FILENAME, 'HG-SRCE',     7,  
     &                 XTAU,      IGLOB,        JGLOB,     
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), Hg2aq )

      !----------------
      ! Print totals
      !----------------
      WRITE( 6, '(a)' ) 
      WRITE( 6, 111   ) SUM( Hg0aq ) * 1d-6
      WRITE( 6, 112   ) SUM( Hg2aq ) * 1d-6

      ! FORMAT strings
 111  FORMAT( 'Initial Oceanic Hg(0)  : ', f7.3, ' [Gg]' )
 112  FORMAT( 'Initial Oceanic Hg(II) : ', f7.3, ' [Gg]' )
 
      ! Return to calling program
      END SUBROUTINE INIT_OCEAN_MERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_OCEAN_MERCURY
!
!******************************************************************************
!  Subroutine CLEANUP_OCEAN_MERCURY deallocates all arrays.  The final oceanic
!  masses of Hg(0) and Hg(II) are also written to disk. (sas, bmy, 1/20/05)
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE FILE_MOD,   ONLY : IU_FILE
      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,   ONLY : EXPAND_DATE, GET_NYMD, GET_NHMS, 
     &                       GET_TAU,     GET_TAUb
      USE TRACER_MOD, ONLY : ITS_A_MERCURY_SIM

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER            :: NYMD,      NHMS
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: IFIRST,    JFIRST,   LFIRST
      REAL*4             :: LONRES,    LATRES,   ARRAY(IGLOB,JGLOB,1)
      REAL*8             :: TAU,       TAUb
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY,  UNIT,     RESERVED
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! CLEANUP_OCEAN_MERCURY begins here!
      !=================================================================

      ! Get TAU values at start of run and at current time
      TAUb = GET_TAUb()
      TAU  = GET_TAU()

      ! Save ocean mercury concentrations (if TAU > TAUb)
      IF ( ITS_A_MERCURY_SIM() .and. TAU > TAUb ) THEN 

         ! Initialize values
         IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
         JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
         LFIRST    = 1
         HALFPOLAR = 1
         CENTER180 = 1
         LONRES    = DISIZE
         LATRES    = DJSIZE
         MODELNAME = GET_MODELNAME()
         NYMD      = GET_NYMD()
         NHMS      = GET_NHMS()
         CATEGORY  = 'HG-SRCE'
         RESERVED  = ''
         UNIT      = 'kg'

         ! Expand date in filename
         FILENAME  = Hg0_Hg2_FILE
         CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

         ! Open BPCH file for output
         CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME )

         !-----------------------
         ! Final Oceanic Hg(0)
         !-----------------------

         ! Cast from REAL*8 to REAL*4
         ARRAY(:,:,1) = Hg0aq

         ! Save to disk
         CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 2,    
     &               UNIT,      TAU,       TAU,      RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !-----------------------
         ! Final Oceanic Hg(II)
         !-----------------------

         ! Cast from REAL*8 to REAL*4
         ARRAY(:,:,1) = Hg2aq

         ! Save to disk
         CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 7,    
     &               UNIT,      TAU,       TAU,      RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         ! Close file
         CLOSE( IU_FILE )
      ENDIF

      !=================================================================
      ! Deallocate arrays
      !=================================================================
      IF ( ALLOCATED( DD_Hg2  ) ) DEALLOCATE( DD_Hg2  )
      IF ( ALLOCATED( Hg0aq   ) ) DEALLOCATE( Hg0aq   )  
      IF ( ALLOCATED( Hg2aq   ) ) DEALLOCATE( Hg2aq   )
      IF ( ALLOCATED( MLD     ) ) DEALLOCATE( MLD     )
      IF ( ALLOCATED( NPP     ) ) DEALLOCATE( NPP     )
      !-----------------------------------------------------------------
      ! Uncomment if using monthly-mean RADSWG fields
      !IF ( ALLOCATED( RAD     ) ) DEALLOCATE( RAD     )
      !-----------------------------------------------------------------
      IF ( ALLOCATED( WD_Hg2  ) ) DEALLOCATE( WD_Hg2  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_OCEAN_MERCURY

!------------------------------------------------------------------------------
     
      ! End of module
      END MODULE OCEAN_MERCURY_MOD
