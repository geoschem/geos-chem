! $Id: ocean_mercury_mod.f,v 1.7 2006/03/24 20:22:55 bmy Exp $
      MODULE OCEAN_MERCURY_MOD
!
!******************************************************************************
!  Module OCEAN_MERCURY_MOD contains variables and routines needed to compute
!  the oceanic flux of mercury.  Original code by Sarah Strode at UWA/Seattle.
!  (sas, bmy, 1/21/05, 3/20/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) Hg_RST_FILE (CHAR   )  : Name of restart file with ocean tracers
!  (2 ) USE_CHECKS  (LOGICAL)  : Flag for turning on error-checking  
!  (3 ) MAX_RELERR  (REAL*8 )  : Max error for total-tag error check [unitless]
!  (4 ) MAX_ABSERR  (REAL*8 )  : Max abs error for total-tag err chk [unitless]
!  (5 ) MAX_FLXERR  (REAL*8 )  : Max error tol for flux error check  [unitless]
!  (6 ) CDEEP       (REAL*8 )  : Conc. of Hg0, Hg2, HgC below MLD    [pM      ]
!  (7 ) DD_Hg2      (REAL*8 )  : Array for Hg(II) dry dep'd to ocean [kg      ]
!  (8 ) dMLD        (REAL*8 )  : Array for Change in ocean MLD       [cm      ]
!  (9 ) Hg0aq       (REAL*8 )  : Array for ocean mass of Hg(0)       [kg      ]
!  (10) Hg2aq       (REAL*8 )  : Array for ocean mass of Hg(II)      [kg      ]
!  (11) HgC         (REAL*8 )  : Array for ocean mass of HgC         [kg      ]
!  (12) MLD         (REAL*8 )  : Array for instantaneous ocean MLD   [cm      ]
!  (13) MLDav       (REAL*8 )  : Array for monthly mean ocean MLD    [cm      ]
!  (14) newMLD      (REAL*8 )  : Array for next month's ocean MLD    [cm      ]
!  (15) NPP         (REAL*8 )  : Array for mean net primary prod.    [unitless]
!  (16) RAD         (REAL*8 )  : Array for mean solar radiation      [W/m2    ]
!  (17) UPVEL       (REAL*8 )  : Array for ocean upwelling velocity  [m/s     ]
!  (18) WD_Hg2      (REAL*8 )  : Array for Hg(II) wet dep'd to ocean [kg      ]
!
!  Module Routines:
!  ============================================================================
!  (1 ) ADD_Hg2_DD             : Archives Hg2 lost to drydep in DD_HG2
!  (2 ) ADD_Hg2_WD             : Archives Hg2 lost to wetdep in WD_HG2
!  (3 ) OCEAN_MERCURY_FLUX     : Routine to compute flux of oceanic mercury
!  (4 ) OCEAN_MERCURY_READ     : Routine to read MLD, NPP, RADSWG data fields
!  (5 ) GET_MLD_FOR_NEXT_MONTH : Routine to read MLD for the next month
!  (6 ) MLD_ADJUSTMENT         : Adjusts MLD 
!  (7 ) READ_OCEAN_Hg_RESTART  : Reads restart file with ocean Hg tracers
!  (8 ) CHECK_DIMENSIONS       : Checks dims of data blocks from restart file
!  (9 ) CHECK_DATA_BLOCKS      : Checks for missing/multiple data blocks
!  (10) MAKE_OCEAN_Hg_RESTART  : Writes new restart file with ocean Hg tracers
!  (11) CHECK_ATMOS_MERCURY    : Checks mass of total & tagged atm Hg0 & Hg2 
!  (12) CHECK_OCEAN_MERCURY    : Checks mass of total & tagged oc Hg0 & Hg2
!  (13) CHECK_OCEAN_FLUXES     : Checks mass of total & tagged DD & WD fluxes
!  (14) INIT_OCEAN_MERCURY     : Allocates and zeroes all module variables
!  (15) CLEANUP_OCEAN_MERCURY  : Deallocates all module variables
!
!  GEOS-CHEM modules referenced by ocean_mercury_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f            : Module w/ routines for binary pch file I/O
!  (2 ) dao_mod.f              : Module w/ arrays for DAO met fields
!  (3 ) diag03_mod.f           : Module w/ ND03 diagnostic arrays 
!  (2 ) file_mod.f             : Module w/ file unit numbers and error checks
!  (9 ) grid_mod.f             : Module w/ horizontal grid information
!  (10) logical_mod.f          : Module w/ GEOS-CHEM logical switches
!  (11) pressure_mod.f         : Module w/ routines to compute P(I,J,L)
!  (12) time_mod.f             : Module w/ routines to compute date & time
!  (13) tracer_mod.f           : Module w/ GEOS-CHEM tracer array STT etc.
!  (14) tracerid_mod.f         : Module w/ pointers to tracers & emissions
!  (15) transfer_mod.f         : Module w/ routines to cast & resize arrays
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
!  (3 ) HgC               : Colloidal mercury
!
!  NOTES:
!  (1 ) Modified ocean flux w/ Sarah's new Ks value (sas, bmy, 2/24/05)
!  (2 ) Now get HALFPOLAR for GCAP or GEOS grids (bmy, 6/28/05)
!  (3 ) Now can read data for both GCAP or GEOS grids (bmy, 8/16/05)
!  (4 ) Include updates from S. Strode and C. Holmes (cdh, sas, bmy, 3/20/06)
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
      PUBLIC :: READ_OCEAN_Hg_RESTART
      PUBLIC :: MAKE_OCEAN_Hg_RESTART

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      LOGICAL              :: USE_CHECKS
      CHARACTER(LEN=255)   :: Hg_RST_FILE

      ! Parameters
      REAL*4,  PARAMETER   :: MAX_RELERR =  5.0d-2
      REAL*4,  PARAMETER   :: MAX_ABSERR =  5.0d-3
      REAL*4,  PARAMETER   :: MAX_FLXERR =  5.0d-1 
      REAL*8,  PARAMETER   :: CDEEP(3)   =  (/ 6d-11, 5d-10, 5d-10 /) 

      ! Arrays
      REAL*8,  ALLOCATABLE :: DD_Hg2(:,:,:)
      REAL*8,  ALLOCATABLE :: dMLD(:,:)
      REAL*8,  ALLOCATABLE :: Hg0aq(:,:,:)
      REAL*8,  ALLOCATABLE :: Hg2aq(:,:,:)
      REAL*8,  ALLOCATABLE :: HgC(:,:)
      REAL*8,  ALLOCATABLE :: MLD(:,:)
      REAL*8,  ALLOCATABLE :: MLDav(:,:)
      REAL*8,  ALLOCATABLE :: newMLD(:,:)
      REAL*8,  ALLOCATABLE :: NPP(:,:)
      REAL*8,  ALLOCATABLE :: RAD(:,:)
      REAL*8,  ALLOCATABLE :: UPVEL(:,:)
      REAL*8,  ALLOCATABLE :: WD_Hg2(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE ADD_Hg2_DD( I, J, N, DRY_Hg2 )
!
!******************************************************************************
!  Subroutine ADD_Hg2_WD computes the amount of Hg(II) dry deposited 
!  out of the atmosphere into the column array DD_Hg2. 
!  (sas, cdh, bmy, 1/19/05, 1/9/06)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) I       (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J       (INTEGER) : GEOS-CHEM latitude  index
!  (3 ) N       (INTEGER) : GEOS-CHEM tracer    index
!  (4 ) DRY_Hg2 (REAL*8 ) : Hg(II) dry deposited out of the atmosphere [kg]
!
!  NOTES:
!  (1 ) DD_Hg2 is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_Hg2_CAT to return the Hg category #. (cdh, bmy, 1/9/06)
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_MOD,  ONLY : LDYNOCEAN
      USE TRACERID_MOD, ONLY : GET_Hg2_CAT

      ! Arguments as input
      INTEGER, INTENT(IN)   :: I, J, N
      REAL*8,  INTENT(IN)   :: DRY_Hg2
 
      ! Local variables
      INTEGER               :: NN

      !=================================================================
      ! ADD_Hg2_DD begins here!
      !=================================================================

      ! Get the index for DD_Hg2 based on the tracer number
      NN = GET_Hg2_CAT( N )

      ! Store dry deposited Hg(II) into DD_Hg2 array
      IF ( NN > 0 ) THEN
         DD_Hg2(I,J,NN) = DD_Hg2(I,J,NN) + DRY_Hg2
      ENDIF

      ! Return to calling program
      END SUBROUTINE ADD_Hg2_DD

!------------------------------------------------------------------------------

      SUBROUTINE ADD_Hg2_WD( I, J, N, WET_Hg2 )
!
!******************************************************************************
!  Subroutine ADD_Hg2_WD computes the amount of Hg(II) wet scavenged 
!  out of the atmosphere into the column array WD_Hg2. 
!  (sas, cdh, bmy, 1/19/05, 1/9/06)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) I       (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J       (INTEGER) : GEOS-CHEM latitude  index
!  (3 ) N       (INTEGER) : GEOS-CHEM tracer    index
!  (4 ) WET_Hg2 (REAL*8 ) : Hg(II) scavenged out of the atmosphere
!
!  NOTES:
!  (1 ) DD_Hg2 is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_Hg2_CAT to return the Hg category #. (cdh, bmy, 1/9/06)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : GET_Hg2_CAT

      ! Arguments as input
      INTEGER, INTENT(IN)   :: I, J, N
      REAL*8,  INTENT(IN)   :: WET_Hg2
 
      ! Local variables
      INTEGER               :: NN

      !=================================================================
      ! ADD_Hg2_WD begins here!
      !=================================================================

      ! Get Hg2 category number
      NN = GET_Hg2_CAT( N ) 

      ! Store wet deposited Hg(II) into WD_Hg2 array
      IF ( NN > 0 ) THEN
         WD_Hg2(I,J,NN) = WD_Hg2(I,J,NN) + WET_Hg2
      ENDIF

      ! Return to calling program
      END SUBROUTINE ADD_Hg2_WD

!------------------------------------------------------------------------------

      SUBROUTINE OCEAN_MERCURY_FLUX( FLUX )
!
!******************************************************************************
!  Subroutine OCEAN_MERCURY_FLUX calculates emissions of Hg(0) from 
!  the ocean in [kg/s].  (sas, bmy, 1/19/05, 1/9/06)
!
!  NOTE: The emitted flux may be negative when ocean conc. is very low. 
!    
!  Arguments as Output
!  ============================================================================
!  (1 ) FLUX (REAL*8) : Flux of Hg(0) from the ocean [kg/s]
!
!  NOTES:
!  (1 ) Change Ks to make ocean flux for 2001 = 2.03e6 kg/year.
!        (sas, bmy, 2/24/05)
!  (2 ) Rewritten to include Sarah Strode's latest ocean Hg model code.
!        (sas, cdh, bmy, 1/9/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : AIRVOL, ALBD, TS, RADSWG
      USE DIAG03_MOD,    ONLY : AD03, ND03
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE GRID_MOD,      ONLY : GET_AREA_M2
      USE LOGICAL_MOD,   ONLY : LSPLIT
      USE TIME_MOD,      ONLY : GET_TS_EMIS,     GET_MONTH 
      USE TIME_MOD,      ONLY : ITS_A_NEW_MONTH, ITS_MIDMONTH
      USE TRACER_MOD,    ONLY : STT,             TRACER_MW_KG
      USE TRACERID_MOD,  ONLY : ID_Hg_tot,       ID_Hg_oc
      USE TRACERID_MOD,  ONLY : ID_Hg0,          N_Hg_CATS
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DEP"       ! FRCLND

      ! Arguments 
      REAL*8,  INTENT(OUT)  :: FLUX(IIPAR,JJPAR,N_Hg_CATS)

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      CHARACTER(LEN=255)    :: FILENAME
      INTEGER               :: I,         J,        NN, C
      INTEGER               :: N,         N_tot_oc
      INTEGER               :: NEXTMONTH, THISMONTH

      REAL*8                :: A_M2,     DTSRCE,   MLDCM,    MLDCMC
      REAL*8                :: CHg0aq,   CHg0
      REAL*8                :: K1,       Kc,       Ks,       Kw
      REAL*8                :: Kcon,     Ksink,    TC,       TK       
      REAL*8                :: Sc,       ScCO2,    USQ,      MHg
      REAL*8                :: Hg2_RED,  Hg2_GONE, Hg2_CONV, HgC_SUNK
      REAL*8                :: FRAC_L,   FRAC_O,   H,        TOTDEP
      REAL*8                :: EF,       oldMLD,   XTAU
      REAL*8                :: E_CONV,   E_SINK,   E_RED
      REAL*8                :: DIFFUSION(3) 

      ! Conversion factor from [cm/h * ng/L] --> [kg/m2/s]
      REAL*8,  PARAMETER    :: TO_KGM2S = 1.0D-11 / 3600D0 

      ! External functions
      REAL*8,  EXTERNAL     :: SFCWINDSQR 
     
      !=================================================================
      ! OCEAN_MERCURY_FLUX begins here!
      !=================================================================

      ! Loop limit for use below
      IF ( LSPLIT ) THEN
         N_tot_oc = 2
      ELSE
         N_tot_oc = 1
      ENDIF

      ! Molecular weight of Hg (applicable to all tagged tracers)
      MHg = TRACER_MW_KG(ID_Hg_tot)

      !-----------------------------------------------
      ! Check tagged & total sums (if necessary)
      !-----------------------------------------------
      IF ( USE_CHECKS .and. LSPLIT ) THEN
         CALL CHECK_ATMOS_MERCURY( 'start of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_MERCURY( 'start of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_FLUXES ( 'start of OCEAN_MERCURY_FLUX' )
      ENDIF

      !-----------------------------------------------
      ! Read monthly NPP, RADSW, MLD, UPVEL data
      !-----------------------------------------------
      IF ( ITS_A_NEW_MONTH() ) THEN

         ! Get current month
         THISMONTH = GET_MONTH()

         ! Get monthly MLD, NPP, etc.
         CALL OCEAN_MERCURY_READ( THISMONTH )

      ENDIF    
     
      !-----------------------------------------------
      ! MLD and entrainment change in middle of month
      !-----------------------------------------------
      IF ( ITS_MIDMONTH() ) THEN

         ! Get current month
         THISMONTH = GET_MONTH()

         ! Read next month's MLD
         CALL GET_MLD_FOR_NEXT_MONTH( THISMONTH )
         
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

      !--------------------------------------------------------------
      ! Determine Ks (sinking term) [unitless]
      !--------------------------------------------------------------
      ! Prior to 2/24/05
      ! Change Ks to make ocean flux for 2001 = 2.03e6 kg/year.
      ! (sas, bmy, 2/24/05)
      !Ks     = 8.4d-7 * DTSRCE 
      !--------------------------------------------------------------
#if   defined( GEOS_4 )
      Ks = 7.7d-21 * DTSRCE 
#else 
      Ks = 5.2d-8  * DTSRCE
#endif  

      ! Hg2 --> colloidal conversion rate
      Kc = 3.1d-22 * DTSRCE    

      ! Diffused mass of (Hg0, Hg2, HgC) across thermocline [kg/m2/timestep]
      ! Based on a fixed gradient at the thermocline
      ! DIFFUSION  = (Diff. coeff.) * (Gradient) * (Hg molar mass) * DT
      DIFFUSION(1) = 5.0d-5 * 3.0d-12 * MHg * DTSRCE  
      DIFFUSION(2) = 5.0d-5 * 5.0d-12 * MHg * DTSRCE
      DIFFUSION(3) = 5.0d-5 * 5.0d-12 * MHg * DTSRCE

!### Debug output
!      ! Write Ocean masses to log file
!      WRITE( 6, 101) SUM( Hg0aq(:,:,ID_Hg_tot) ), 
!     &               SUM( Hg2aq(:,:,ID_Hg_tot) ),
!     &               SUM( HgC                  )

 101  FORMAT( '     - TOTAL OCEAN MASS: ',
     &        'Hg0 ', ES9.2, ' |  Hg2 ', ES9.2, ' |  HgC ', ES9.2 )

      ! Loop over latitudes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,      A_M2,     OLDMLD,   MLDCM    )
!$OMP+PRIVATE( MLDCMC,  Kw,     K1,       HgC_SUNK, Hg2_CONV )
!$OMP+PRIVATE( TK,      TC,     FRAC_L,   FRAC_O,   H        )
!$OMP+PRIVATE( Sc,      ScCO2,  EF,       Ksink,    Kcon     )
!$OMP+PRIVATE( Usq,     C,      NN,       E_RED,    E_CONV   )
!$OMP+PRIVATE( E_SINK,  TOTDEP, Hg2_RED,  Hg2_GONE, N        )
!$OMP+PRIVATE( CHg0aq,  CHg0                                 )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Grid box surface area [m2]
         A_M2 = GET_AREA_M2( J )

      ! Loop over longitudes
      DO I = 1, IIPAR

         ! Initialize values
         OLDMLD     = MLDav(I,J)
         MLDav(I,J) = MLDav(I,J) + dMLD(I,J) * DTSRCE
         MLDcm      = MLDav(I,J)
         MLDcmc     = MLDCM
         Kw         = 0d0
         K1         = 0d0
         HgC_SUNK   = 0d0
         Hg2_CONV   = 0d0
         TK         = 0d0
         TC         = 0d0

         ! Get fractions of land and ocean in the grid box [unitless]
         FRAC_L     = FRCLND(I,J)
         FRAC_O     = 1d0 - FRAC_L

         ! Change ocean mass due to mixed layer depth change
         ! Keep before next IF so that we adjust mass in ice-covered boxes 
         CALL MLD_ADJUSTMENT( I, J, OLDMLD*1d-2, MLDcm*1d-2 )

         !===========================================================
         ! Make sure we are in an ocean box
         !===========================================================
         IF ( ( ALBD(I,J) <= 0.4d0 ) .and. 
     &        ( FRAC_L    <  0.8d0 ) .and.
     &        ( MLDCM     > 0.99d0 )      ) THEN

            !--------------------------------------------------------
            ! Calculate K1 (reduction) based on NPP & RAD
            !--------------------------------------------------------

            ! For Daily RADSWG fields
            K1     = 3.1D-24 * DTSRCE * NPP(I,J) * RADSWG(I,J)
     &             * A_M2 * FRAC_O !for Reed ScHg               
               
            ! Surface air temperature in both [K] and [C]
            ! (Use as surrogate for SST, cap at freezing point)
            TK     = MAX( TS(I,J), 273.15d0 )
            TC     = TK - 273.15d0

            ! Henry's law constant (liquid->gas) [unitless] [L air/L water]  
            H      = EXP( 4633.3d0 / TK - 14.52d0 )

            ! Schmidt # for Hg [unitless] 
            Sc     = ( 0.017d0 * EXP( -0.025d0 * TC ) ) / 
     &               ( 7.4D-8 * sqrt(2.6*18.0)*TK / 14.8)  ! Reid
            
            ! Schmidt # of CO2 [unitless]
            ScCO2  = 644.7d0 + TC * ( -6.16d0 + TC * ( 0.11d0 ) ) 

            ! EF ratio for particle sinking based on Laws et al. 2000 
            EF     = MAX( (0.63 - 0.02 * TC), 0.0) ! keep export > 0
            Ksink  = Ks * EF * NPP(I,J) * A_M2 *FRAC_O

            ! Conversion rate Hg2 -> HgC [unitless]
            Kcon   = Kc * NPP(I,J) * A_M2 * FRAC_O

            ! Square of surface (actually 10m) wind speed [m2/s2]
            Usq    = SFCWINDSQR(I,J)

            ! Piston velocity [cm/h], now from Nightingale
            Kw     = ( 0.25d0 * Usq ) / SQRT( Sc / ScCO2 )

            ! Cap mixed layer depth for Hg2 reduction at 100m
            MLDcmc = MIN( MLDcmc, 1d4 )

            !-----------------------------------------------------------
            ! Physical transport for tracers, Part I:
            ! Diffusion from below thermocline
            !-----------------------------------------------------------

            ! Loop over total Hg (and ocean Hg if necessary)
            DO C = 1, N_tot_oc

               ! Get Hg category #
               IF ( C == 1 ) NN = ID_Hg_tot
               IF ( C == 2 ) NN = ID_Hg_oc

               ! Hg0
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN) 
     &                       + ( DIFFUSION(1) * A_M2 * FRAC_O )

               ! Hg2
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                       + ( DIFFUSION(2) * A_M2 * FRAC_O )
 
               ! Hg colloidal
               IF ( C == 1 ) THEN
                  HgC(I,J)   = HgC(I,J) 
     &                       + ( DIFFUSION(3) * A_M2 * FRAC_O )
               ENDIF

            ENDDO

            !-----------------------------------------------------------
            ! Physical transport for tracers, Part II:
            ! Upward current transport (Ekman pumping)
            ! Upward mass flux is:
            ! Mass = (Vol upwelling water) * (Conc. below thermocline)
            ! Mass = (VEL * AREA * TIME  ) * (C * Molar Mass )
            !-----------------------------------------------------------
            IF ( UPVEL(I,J) > 0d0 ) THEN
                 
               ! Loop over total Hg (and ocean Hg if necessary)
               DO C = 1, N_tot_oc
 
                  ! Get Hg category #
                  IF ( C == 1 ) NN = ID_Hg_tot
                  IF ( C == 2 ) NN = ID_Hg_oc

                  ! Hg0 
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeep(1) )

                  ! Hg2 
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeep(2) )

                  ! Hg colloidal
                  IF ( C == 1 ) THEN
                     HgC(I,J)   = HgC(I,J) + UPVEL(I,J) 
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeep(3) )
                  ENDIF
                  
               ENDDO
                  
            !----------------------------------------------------------
            ! Physical transport for TOTAL TRACERS, Part III:
            ! Downward current transport (Ekman pumping)
            ! Treated as a deposition velocity
            ! d(Mass)/dt = - VEL * Mass / BoxHeight
            !----------------------------------------------------------
            ELSE  

               ! Loop over all types of tagged tracers
               DO NN = 1, N_Hg_CATS

                  ! Hg0
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) 
     &                * ( 1d0 + UPVEL(I,J) * DTSRCE / ( MLDcm * 1d-2 ) )
                  
                  ! Hg2
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) 
     &                * ( 1d0 + UPVEL(I,J) * DTSRCE / ( MLDcm * 1d-2 ) )
               
                  ! Hg colloidal
                  IF ( NN == 1 ) THEN
                     HgC(I,J)  = HgC(I,J) 
     &                * ( 1d0 + UPVEL(I,J) * DTSRCE / ( MLDcm * 1d-2 ) )
                  ENDIF

               ENDDO

            ENDIF

            !===========================================================
            ! Calculate reduction, conversion, sinking, evasion
            !
            ! (1) Hg2 -> HgC and HgC sinks
            ! (2) Hg2 -> Hg0 and Hg0 evades
            !
            ! NOTE: N is the GEOS-CHEM tracer # (for STT)
            !       and NN is the Hg category # (for Hg0aq, Hg2aq, HgC)
            !===========================================================

            ! Loop over all Hg categories
            DO NN = 1, N_Hg_CATS

               ! Reset flux each timestep
               FLUX(I,J,NN)  = 0d0 

               !--------------------------------------------------------
               ! Precompute exponents
               !--------------------------------------------------------

               ! Exponent for reduction [unitless]
               E_RED         = EXP( -K1 * MLDCMC / MLDCM )

               ! Exponent for conversion of Hg(II) -> Hg(C) [unitless]
               E_CONV        = EXP( - Kcon )

               ! Exponent for sinking Hg(C) --> deep ocean [unitless]
               E_SINK        = EXP( - Ksink )
                  
               !--------------------------------------------------------
               ! Calculate new Hg(II) mass
               !--------------------------------------------------------

               ! Total Hg(II) deposited on ocean surface [kg]
               TOTDEP        = (WD_Hg2(I,J,NN) + DD_Hg2(I,J,NN))*FRAC_O 

               ! Add deposited Hg(II) to the Hg(II) ocean mass [kg]
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + TOTDEP        

               ! Mass of Hg(II)  -->  Hg(C) 
               Hg2_CONV      = Hg2aq(I,J,NN) * ( 1d0 - E_CONV ) 

               ! Mass of Hg(II)  -->  Hg(0) 
               Hg2_RED       = Hg2aq(I,J,NN) * ( 1d0 - E_RED ) 

               ! Amount of Hg(II) that is lost [kg]
               Hg2_GONE      = Hg2_CONV + Hg2_RED

               ! Cap Hg2_GONE with available Hg2
               IF ( Hg2_GONE > Hg2aq(I,J,NN) ) THEN 
                  Hg2_GONE   = MIN( Hg2_GONE, Hg2aq(I,J,NN) )
               ENDIF

               ! Hg(II) ocean mass after reduction and conversion [kg]
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN) - Hg2_GONE

               !--------------------------------------------------------
               ! Calculate new Hg(C) mass
               !--------------------------------------------------------
               IF ( NN == 1 ) THEN

                  ! HgC ocean mass after conversion
                  HgC(I,J)   = HgC(I,J) + Hg2_CONV
                     
                  ! Archive Hg(C) sinking loss for ND03 [kg]
                  HgC_SUNK   = HgC(I,J) * ( 1d0 - E_SINK )

                  ! HgC ocean mass after sinking [kg]
                  HgC(I,J)   = HgC(I,J) - HgC_SUNK

                  ! Store Hg2_CONV for total tracer only
                  IF ( ND03 > 0 ) THEN
                     AD03(I,J,12) = AD03(I,J,12) + Hg2_CONV 
                  ENDIF

               ENDIF

               !--------------------------------------------------------
               ! Calculate new Hg(0) mass
               !--------------------------------------------------------

               ! Hg0 tracer number (for STT)
               N             = ID_Hg0(NN)

               ! Add converted Hg(II) to the ocean mass of Hg(0) [kg]
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + Hg2_RED     

               !--------------------------------------------------------
               ! Calculate oceanic and gas-phase concentration of Hg(0)
               !--------------------------------------------------------
                  
               ! Concentration of Hg(0) in the ocean [ng/L]
               CHg0aq        = ( Hg0aq(I,J,NN) * 1d11   ) /   
     &                         ( A_M2          * FRAC_O ) / MLDcm 
               
               ! Gas phase Hg(0) concentration: convert [kg] -> [ng/L]
               CHg0          = STT(I,J,1,N) * 1.0D9 / AIRVOL(I,J,1)
               
               !--------------------------------------------------------
               ! Compute flux of Hg(0) from the ocean to the air
               !--------------------------------------------------------

               ! Compute ocean flux of Hg0 [cm/h*ng/L]
               FLUX(I,J,NN)  = Kw * ( CHg0aq - ( H * CHg0 ) ) 

               ! Convert [cm/h*ng/L] --> [kg/m2/s] --> [kg/s]
               ! Also account for ocean fraction of grid box
               FLUX(I,J,NN)  = FLUX(I,J,NN) * TO_KGM2S * A_M2 *FRAC_O 

               !--------------------------------------------------------
               ! Flux limited by ocean and atm Hg(0)
               !--------------------------------------------------------

               ! Cap the flux w/ the available Hg(0) ocean mass
               IF ( FLUX(I,J,NN) * DTSRCE > Hg0aq(I,J,NN) ) THEN 
                  FLUX(I,J,NN) = Hg0aq(I,J,NN) / DTSRCE        
               ENDIF
               
               ! Cap the neg flux w/ the available Hg(0) atm mass
               IF ( (-FLUX(I,J,NN) * DTSRCE ) > STT(I,J,1,N) ) THEN
                  FLUX(I,J,NN) = -STT(I,J,1,N) / DTSRCE       
               ENDIF
                
               !--------------------------------------------------------
               ! Remove amt of Hg(0) that is leaving the ocean [kg]
               !--------------------------------------------------------
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN) - ( FLUX(I,J,NN) * DTSRCE ) 
               
            ENDDO   

            !-----------------------------------------------------------
            ! ND03 diagnostics ("OCEAN-HG")
            !-----------------------------------------------------------
            IF ( ND03 > 0 ) THEN

               ! Aqueous Hg(0) mass [kg]
               AD03(I,J,2)  = AD03(I,J,2)  + Hg0aq(I,J,ID_Hg_tot) 

               ! Aqueous Hg(II) mass [kg] 
               AD03(I,J,7)  = AD03(I,J,7)  + Hg2aq(I,J,ID_Hg_tot) 

               ! Hg2 sunk deep into the ocean [kg]
               AD03(I,J,8)  = AD03(I,J,8)  + HgC_SUNK

               ! Kw (piston velocity) [cm/s]
               AD03(I,J,10) = AD03(I,J,10) + Kw

               ! Hg converted to colloidal [kg/m2/s]
               AD03(I,J,11) = AD03(I,J,11) + HgC(I,J) 
            ENDIF
            
         !==============================================================
         ! If we are not in an ocean box, set Hg(0) flux to zero
         !==============================================================
         ELSE

            DO NN = 1, N_Hg_CATS 
               FLUX(I,J,NN) = 0d0
            ENDDO               

         ENDIF 
      
         !==============================================================
         ! Zero amts of deposited Hg2 for next timestep [kg]  
         !==============================================================
         DO NN = 1, N_Hg_CATS  
            DD_Hg2(I,J,NN) = 0d0
            WD_Hg2(I,J,NN) = 0d0
         ENDDO                

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Check tagged & total sums (if necessary)
      !=================================================================
      IF ( USE_CHECKS .and. LSPLIT ) THEN
         CALL CHECK_ATMOS_MERCURY(  'end of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_MERCURY(  'end of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_FLUXES (  'end of OCEAN_MERCURY_FLUX' )
         CALL CHECK_FLUX_OUT( FLUX, 'end of OCEAN_MERCURY_FLUX' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE OCEAN_MERCURY_FLUX

!------------------------------------------------------------------------------

      SUBROUTINE OCEAN_MERCURY_READ( THISMONTH )
!
!******************************************************************************
!  Subroutine OCEAN_MERCURY_READ reads in the mixed layer depth, net primary 
!  productivity, upwelling and radiation climatology for each month.  
!  This is needed for the ocean flux computation. 
!  (sas, cdh, bmy, 1/20/05, 1/9/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Month to read fields (1-12)
!
!  NOTES:
!  (1 ) Modified for S. Strode's latest ocean Hg code.  Now read files
!        from DATA_DIR_1x1/mercury_200511. (sas, cdh, bmy, 1/9/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: THISMONTH

      ! Local Variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME
     
      !=================================================================
      ! OCEAN_MERCURY_READ begins here!
      !=================================================================
     
      !------------------------------
      ! Mixed layer depth [cm]
      !------------------------------

      ! MLD file name
      FILENAME = TRIM( DATA_DIR_1x1 )       // 
     &           'mercury_200511/mld.geos.' // GET_RES_EXT()      
     
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - OCEAN_MERCURY_READ: Reading ', a )  

      ! TAU0 value (uses year 1985)
      TAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read from disk; original units are [m]
      CALL READ_BPCH2( FILENAME, 'BXHGHT-$',    5,  
     &                 TAU,       IGLOB,        JGLOB,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), MLD )

      ! Convert [m] to [cm]
      MLD = MLD * 100d0

      ! First-time only: Set MDLav [cm] to MLD of first month
      IF ( FIRST ) THEN   
         MLDav = MLD   
         dMLD  = 0.0
         FIRST = .FALSE.
      ENDIF

      !--------------------------------
      ! Net primary productivity 
      !--------------------------------
 
      ! NPP file name
      FILENAME = TRIM( DATA_DIR_1x1 )             // 
     &           'mercury_200511/modis_npp.geos.' // GET_RES_EXT() 
     
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU0 values (uses year 2003)
      TAU = GET_TAU0( THISMONTH, 1, 2003 )
 
      ! Read data
      CALL READ_BPCH2( FILENAME, 'GLOB-NPP',    1,  
     &                 TAU,       IGLOB,        JGLOB,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), NPP )
  
      !---------------------------------
      ! Ekman upwelling velocity [cm/s]
      !---------------------------------

      ! NPP file name
      FILENAME = TRIM( DATA_DIR_1x1 )               // 
     &           'mercury_200511/ekman_upvel.geos.' // GET_RES_EXT() 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU0 value (uses year 1985)
      TAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read from disk; original units are [cm/s]
      CALL READ_BPCH2( FILENAME, 'EKMAN-V',     1,  
     &                 TAU,       IGLOB,        JGLOB,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), UPVEL )

      ! convert [cm/s] to [m/s]
      UPVEL = UPVEL * 1.D-2
  
      ! Return to calling program
      END SUBROUTINE OCEAN_MERCURY_READ

!------------------------------------------------------------------------------

      SUBROUTINE GET_MLD_FOR_NEXT_MONTH( THISMONTH )
!
!******************************************************************************
!  Subroutine GET_MLD_FOR_NEXT_MONTH reads the mixed-layer depth (MLD) 
!  values for the next month. (sas, cdh, bmy, 1/9/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Now read files from DATA_DIR_1x1/mercury_200511 (bmy, 1/9/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_TAU0, GET_RES_EXT, READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: THISMONTH

      ! Local variables
      INTEGER                :: I, J, NEXTMONTH
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! GET_MLD_FOR_NEXT_MONTH begins here!
      !=================================================================
      
      ! MLD file name
      FILENAME = TRIM( DATA_DIR_1x1 )       // 
     &           'mercury_200511/mld.geos.' // GET_RES_EXT()      

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - GET_MLD_FOR_NEXT_MONTH: Reading ', a )  
      
      ! Get the next month
      NEXTMONTH = MOD( THISMONTH, 12 ) +1

      ! TAU0 value for next month (uses year 1985)
      TAU       = GET_TAU0( NEXTMONTH, 1, 1985 )

      ! Read from disk; original units are [m]
      CALL READ_BPCH2( FILENAME, 'BXHGHT-$',    5,  
     &                 TAU,       IGLOB,        JGLOB,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), newMLD )

      ! Convert [m] to [cm]
      newMLD = newMLD * 100d0

      ! get rate of change of MLD; convert [cm/month] -> [cm/s] 
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         dMLD(I,J) = (newMLD(I,J) - MLD(I,J)) / ( 3.6d3 *24d0 * 30.5d0 )
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE GET_MLD_FOR_NEXT_MONTH

!------------------------------------------------------------------------------

      SUBROUTINE MLD_ADJUSTMENT( I, J, MLDold, MLDnew )
!
!******************************************************************************
!  Subroutine MLD_ADJUSTMENT entrains new water when mixed layer depth deepens
!  and conserves concentration (leaves mass behind) when mixed layer shoals.
!  (sas, cdh, bmy, 4/18/05, 1/9/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I      (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J      (INTEGER) : GEOS-CHEM latitude index
!  (3 ) MLDold (REAL*8 ) : Old ocean mixed layer depth [m]
!  (4 ) MLDnew (REAL*8 ) : New ocean mixed layer depth [m]
!
!  NOTES:
!******************************************************************************
!
      ! Reference to fortran90 modules
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACERID_MOD, ONLY : ID_Hg_tot, ID_Hg_oc, N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DEP"      ! FRCLND
      
      ! Arguments
      INTEGER, INTENT(IN)   :: I, J 
      REAL*8,  INTENT(IN)   :: MLDold, MLDnew  

      ! Local variables
      INTEGER               :: C,    NN,     N_tot_oc    
      REAL*8                :: A_M2, DELTAH, FRAC_O,  MHg

      !=================================================================
      ! MLD_ADJUSTMENT begins here!
      !=================================================================

      ! Loop limit for use below
      IF ( LSPLIT ) THEN
         N_tot_oc = 2
      ELSE
         N_tot_oc = 1
      ENDIF

      ! Grid box surface area [m2]
      A_M2   = GET_AREA_M2( J )

      ! Fraction of box that is ocean
      FRAC_O = 1d0 - FRCLND(I,J)
      
      ! Molecular weight of Hg (valid for all tagged tracers)
      MHg    = TRACER_MW_KG(ID_Hg_tot)

      ! Test if MLD increased
      IF ( MLDnew > MLDold ) THEN

         !==============================================================
         ! IF MIXED LAYER DEPTH HAS INCREASED:
         !
         ! Entrain water with a concentration specified by CDeep
         !
         ! Entrained Mass = ( Vol water entrained ) * CDeep * Molar mass
         !                = ( DELTAH * AREA * FRAC_O ) * CDeep * MHg
         !==============================================================

         ! Increase in MLD [m]
         DELTAH = MLDnew - MLDold

         ! Loop over total Hg (and ocean Hg if necessary)
         DO C = 1, N_tot_oc

            ! Get Hg category number
            IF ( C == 1 ) NN = ID_Hg_tot
            IF ( C == 2 ) NN = ID_Hg_oc
                        
            ! Hg0
            Hg0aq(I,J,NN) = Hg0aq(I,J,NN)
     &                    + ( DELTAH * CDeep(1) * MHg * A_M2 * FRAC_O )

            ! Hg2
            Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                    + ( DELTAH * CDeep(2) * MHg * A_M2 * FRAC_O )

            ! HgC
            IF ( C == 1 ) THEN
               HgC(I,J)   = HgC(I,J)          
     &                    + ( DELTAH * CDeep(3) * MHg * A_M2 * FRAC_O )
            ENDIF

         ENDDO
               
      ELSE 

         !==============================================================
         ! IF MIXED LAYER DEPTH HAS DECREASED:
         !
         ! Conserve concentration, but shed mass for ALL tracers.  
         ! Mass changes by same ratio as volume.
         !==============================================================

         ! Avoid dividing by zero
         IF ( MLDold > 0d0 ) THEN

            ! Update Hg0 and Hg2 categories
            DO NN = 1, N_Hg_CATS
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN) * ( MLDnew / MLDold )
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN) * ( MLDnew / MLDold )
            ENDDO
            
            ! Update colloidal Hg
            HgC(I,J) = HgC(I,J) * ( MLDnew / MLDold )
         
         ENDIF

      ENDIF
      
      ! Return to calling program
      END SUBROUTINE MLD_ADJUSTMENT

!------------------------------------------------------------------------------

      SUBROUTINE READ_OCEAN_Hg_RESTART( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_OCEAN_Hg_RESTART initializes GEOS-CHEM oceanic mercury 
!  tracer masses from a restart file. (sas, cdh, bmy, 1/9/06)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
! 
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : OPEN_BPCH2_FOR_READ
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE FILE_MOD,     ONLY : IU_FILE,     IOERROR
      USE LOGICAL_MOD,  ONLY : LSPLIT,      LPRT
      USE TIME_MOD,     ONLY : EXPAND_DATE
      USE TRACER_MOD,   ONLY : STT,         TRACER_NAME, TRACER_MW_G
      USE TRACERID_MOD, ONLY : GET_Hg0_CAT, GET_Hg2_CAT, N_Hg_CATS
      USE TRACERID_MOD, ONLY : ID_Hg0,      ID_Hg2

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)   :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER               :: I, IOS, J, L, NN, N_oc
      INTEGER               :: YEAR, MONTH, DAY
      INTEGER               :: NCOUNT(NNPAR) 
      REAL*4                :: Hg_OCEAN(IIPAR,JJPAR,1)
      CHARACTER(LEN=255)    :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER               :: NI,        NJ,      NL
      INTEGER               :: IFIRST,    JFIRST,  LFIRST
      INTEGER               :: NTRACER,   NSKIP
      INTEGER               :: HALFPOLAR, CENTER180
      REAL*4                :: LONRES,    LATRES
      REAL*8                :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)     :: MODELNAME
      CHARACTER(LEN=40)     :: CATEGORY
      CHARACTER(LEN=40)     :: UNIT     
      CHARACTER(LEN=40)     :: RESERVED

      !=================================================================
      ! READ_OCEAN_Hg_RESTART begins here!
      !=================================================================

      ! Initialize some variables
      NCOUNT(:)       = 0
      Hg_OCEAN(:,:,:) = 0e0

      ! Copy input file name to a local variable
      FILENAME        = TRIM( Hg_RST_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Echo some input to the screen
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100   ) 
      WRITE( 6, 110   ) TRIM( FILENAME )
 100  FORMAT( 'O C E A N   H g   R E S T A R T   F I L E   I N P U T' )
 110  FORMAT( /, 'READ_OCEAN_Hg_RESTART: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )
      
      ! Echo more output
      WRITE( 6, 120 )
 120  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
     &        /, '(in volume mixing ratio units: v/v)' )
      
      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      DO 
         READ( IU_FILE, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_oc_hg_rst:1' )

         READ( IU_FILE, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_oc_hg_rst:2' )

         READ( IU_FILE, IOSTAT=IOS ) 
     &        ( ( ( Hg_OCEAN(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_oc_hg_rst:3' )

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================
  
         ! Only process concentration data (i.e. mixing ratio)
         IF ( CATEGORY(1:8) == 'OCEAN-HG' ) THEN 

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            CALL CHECK_DIMENSIONS( NI, NJ, NL )

            ! Save into arrays
            IF ( ANY( ID_Hg0 == NTRACER ) ) THEN

               !----------
               ! Hg(0)
               !----------
               
               ! Get the Hg category #
               NN              = GET_Hg0_CAT( NTRACER )

               ! Store ocean Hg(0) in Hg0aq array
               Hg0aq(:,:,NN)   = Hg_OCEAN(:,:,1)
               
               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1

            ELSE IF ( ANY( ID_Hg2 == NTRACER ) ) THEN
               
               !----------
               ! Hg(II)
               !----------

               ! Get the Hg category #
               NN              = GET_Hg2_CAT( NTRACER )

               ! Store ocean Hg(II) in Hg2_aq array
               Hg2aq(:,:,NN)   = Hg_OCEAN(:,:,1)

               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1

            ELSE IF ( NTRACER == 3 ) THEN

               !----------
               ! Hg(C)
               !----------

               ! Colloidal Hg
               HgC(:,:)        = Hg_OCEAN(:,:,1)

               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1

            ENDIF
         ENDIF
      ENDDO

      ! Close file
      CLOSE( IU_FILE )      

      !=================================================================
      ! Examine data blocks, print totals, and return
      !=================================================================

      ! Tagged simulation has 17 ocean tracers; otherwise 3
      IF ( LSPLIT ) THEN
         N_oc = 17
      ELSE
         N_oc = 3
      ENDIF

      ! Check for missing or duplicate data blocks
      CALL CHECK_DATA_BLOCKS( N_oc, NCOUNT )

      !=================================================================
      ! Print totals
      !=================================================================

      ! Echo info
      WRITE( 6, 130 )

      ! Hg0
      DO NN = 1, N_Hg_CATS
         WRITE( 6, 140 ) ID_Hg0(NN), TRACER_NAME( Id_Hg0(NN) ), 
     &                   SUM( Hg0aq(:,:,NN) ),  'kg'
      ENDDO

      ! Hg2
      DO NN = 1, N_Hg_CATS
         WRITE( 6, 140 ) ID_Hg2(NN), TRACER_NAME( Id_Hg2(NN) ), 
     &                   SUM( Hg0aq(:,:,NN) ), 'kg'
      ENDDO

      ! HgC
      WRITE( 6, 140 ) 3, 'HgC       ', SUM( HgC ), 'kg'

      ! Format strings
 130  FORMAT( /, 'Total masses for each ocean tracer: ' ) 
 140  FORMAT( 'Tracer ', i3, ' (', a10, ') ', es12.5, 1x, a4)

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Make sure tagged & total tracers sum up
      IF ( USE_CHECKS .and. LSPLIT ) THEN
         CALL CHECK_OCEAN_MERCURY( 'end of READ_OCEAN_Hg_RESTART' )
      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_OCEAN_Hg_RST: read file' )

      ! Return to calling program
      END SUBROUTINE READ_OCEAN_Hg_RESTART

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_DIMENSIONS( NI, NJ, NL ) 
!
!******************************************************************************
!  Subroutine CHECK_DIMENSIONS makes sure that the dimensions of the Hg 
!  restart file extend to cover the entire grid. (sas, cdh, bmy, 12/16/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NI (INTEGER) : Number of longitudes read from restart file
!  (2 ) NJ (INTEGER) : Number of latitudes  read from restart file
!  (3 ) NL (INTEGER) : Numbef of levels     read from restart file
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NI, NJ, NL

#     include "CMN_SIZE"

      !=================================================================
      ! CHECK_DIMENSIONS begins here!
      !=================================================================

      ! Error check longitude dimension: NI must equal IIPAR
      IF ( NI /= IIPAR ) THEN
         WRITE( 6, 100 ) 
 100     FORMAT( 'ERROR reading in Hg restart file', /
     &           'Wrong number of longitudes encountered', /
     &           'STOP in CHECK_DIMENSIONS ("ocean_mercury_mod.f")' )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Error check latitude dimension: NJ must equal JJPAR
      IF ( NJ /= JJPAR ) THEN
         WRITE( 6, 110 ) 
 110     FORMAT( 'ERROR reading in Hg restart file', /
     &           'Wrong number of longitudes encountered', /
     &           'STOP in CHECK_DIMENSIONS ("ocean_mercury_mod.f")' )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF
      
      ! Error check vertical dimension: NL must equal LLPAR
      IF ( NL /= 1 ) THEN
         WRITE( 6, 120 ) 
 120     FORMAT( 'ERROR reading in Hg restart file', /
     &           'Wrong number of longitudes encountered', /
     &           'STOP in CHECK_DIMENSIONS ("ocean_mercury_mod.f")' )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK_DIMENSIONS

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )
!
!******************************************************************************
!  Subroutine CHECK_DATA_BLOCKS checks to see if we have multiple or 
!  missing data blocks for a given tracer. (sas, cdh, bmy, 1/9/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N_TRACERS (INTEGER) : Number of tracers
!  (2 ) NCOUNT    (INTEGER) : Ctr array - # of data blocks found per tracer
!
!  NOTES:
!******************************************************************************
!      
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: N_TRACERS, NCOUNT(NNPAR)
  
      ! Local variables
      INTEGER             :: N

      !=================================================================
      ! CHECK_DATA_BLOCKS begins here! 
      !=================================================================

      ! Loop over all tracers
      DO N = 1, N_TRACERS

         ! Stop if a tracer has more than one data block 
         IF ( NCOUNT(N) > 1 ) THEN 
            WRITE( 6, 100 ) N
            WRITE( 6, 120 ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP
         ENDIF
         
         ! Stop if a tracer has no data blocks 
         IF ( NCOUNT(N) == 0 ) THEN
            WRITE( 6, 110 ) N
            WRITE( 6, 120 ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP
         ENDIF
      ENDDO

      ! FORMAT statements
 100  FORMAT( 'More than one record found for tracer : ', i4 )
 110  FORMAT( 'No records found for tracer : ',           i4 ) 
 120  FORMAT( 'STOP in CHECK_DATA_BLOCKS (restart_mod.f)'    )

      ! Return to calling program
      END SUBROUTINE CHECK_DATA_BLOCKS

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_OCEAN_Hg_RESTART( NYMD, NHMS, TAU )
!
!******************************************************************************
!  Subroutine MAKE_OCEAN_Hg_RESTART writes an ocean mercury restart file.
!  (sas, cdh, bmy, 1/9/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE FILE_MOD,     ONLY : IU_FILE
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE TIME_MOD,     ONLY : EXPAND_DATE, GET_TAU
      USE TRACERID_MOD, ONLY : ID_Hg_tot,   ID_Hg0
      USE TRACERID_MOD, ONLY : ID_Hg2,      N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments           
      INTEGER, INTENT(IN)   :: NYMD, NHMS
      REAL*8,  INTENT(IN)   :: TAU

      ! Local variables
      INTEGER               :: HALFPOLAR, CENTER180
      INTEGER               :: IFIRST,    JFIRST,   LFIRST
      INTEGER               :: N,         NN
      REAL*4                :: LONRES,    LATRES,   ARRAY(IGLOB,JGLOB,1)
      CHARACTER(LEN=20)     :: MODELNAME
      CHARACTER(LEN=40)     :: CATEGORY,  UNIT,     RESERVED
      CHARACTER(LEN=255)    :: FILENAME

      !=================================================================
      ! MAKE_OCEAN_Hg_RESTART begins here!
      !=================================================================

      ! Initialize values
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LFIRST    = 1
      HALFPOLAR = GET_HALFPOLAR()
      CENTER180 = 1
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      CATEGORY  = 'OCEAN-HG'
      RESERVED  = ''
      UNIT      = 'kg'

      ! Expand date in filename
      FILENAME  = Hg_RST_FILE
      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_RESTART_FILE: Writing ', a )

      ! Open BPCH file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME )

      !---------------------------
      ! Total Hg(0) in ocean
      !---------------------------
      N            = ID_Hg0(Id_Hg_tot)
      ARRAY(:,:,1) = Hg0aq(:,:,ID_Hg_tot)

      CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, N, 
     &            UNIT,      TAU,       TAU,      RESERVED,
     &            IIPAR,     JJPAR,     1,        IFIRST,
     &            JFIRST,    LFIRST,    ARRAY(:,:,1) )

      !---------------------------
      ! Total Hg(II) in ocean
      !---------------------------
      N            = ID_Hg2(ID_Hg_tot)
      ARRAY(:,:,1) = Hg2aq(:,:,ID_Hg_tot)

      CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, N, 
     &            UNIT,      TAU,       TAU,      RESERVED,
     &            IIPAR,     JJPAR,     1,        IFIRST,
     &            JFIRST,    LFIRST,    ARRAY(:,:,1) )

      !---------------------------
      ! Total HgC in ocean
      !---------------------------
      N            = 3
      ARRAY(:,:,1) = HgC(:,:)

      CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, N, 
     &            UNIT,      TAU,       TAU,      RESERVED,
     &            IIPAR,     JJPAR,     1,        IFIRST,
     &            JFIRST,    LFIRST,    ARRAY(:,:,1) )

      ! Save tagged ocean tracers if present
      IF ( LSPLIT ) THEN

         !------------------------
         ! Tagged Hg(0) in ocean
         !------------------------
         DO NN = 2, N_Hg_CATS
            N            = ID_Hg0(NN)
            ARRAY(:,:,1) = Hg0aq(:,:,NN)

            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N, 
     &                  UNIT,      TAU,       TAU,      RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO

         !------------------------
         ! Tagged Hg(II) in ocean
         !------------------------
         DO NN = 2, N_Hg_CATS
            N            = ID_Hg2(NN)
            ARRAY(:,:,1) = Hg2aq(:,:,NN)

            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N, 
     &                  UNIT,      TAU,       TAU,      RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF

      ! Close file
      CLOSE( IU_FILE )

      ! Make sure tagged & total tracers sum up
      IF ( USE_CHECKS .and. LSPLIT ) THEN
         CALL CHECK_OCEAN_MERCURY( 'end of MAKE_OCEAN_Hg_RESTART' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE MAKE_OCEAN_Hg_RESTART

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_ATMOS_MERCURY( LOC )
!
!******************************************************************************
!  Subroutine CHECK_ATMOS_MERCURY tests whether the total and tagged tracers 
!  the GEOS-CHEM tracer array STT sum properly within each grid box.
!  (cdh, bmy, 2/24/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LOC (CHARACTER) : Name of routine where CHECK_ATMOS_MERCURY is called
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD,          ONLY : STT
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE TRACERID_MOD,        ONLY : ID_Hg0,    ID_Hg2,   ID_HgP
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS

#     include "CMN_SIZE"            ! Size parameters

      ! Arguments as Input
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL                      :: FLAG
      INTEGER                      :: I,       J,       L
      INTEGER                      :: N,       NN
      REAL*8                       :: Hg0_tot, Hg0_tag, RELERR0, ABSERR0      
      REAL*8                       :: Hg2_tot, Hg2_tag, RELERR2, ABSERR2
      REAL*8                       :: HgP_tot, HgP_tag, RELERRP, ABSERRP

      !=================================================================
      ! CHECK_ATMOS_MERCURY begins here!
      !=================================================================

      ! Set error flags
      FLAG = .FALSE.

      ! Loop over grid boxes
! OMP PARALLEL DO
! OMP+DEFAULT( SHARED )
! OMP+PRIVATE( I,       J,       L,       N,      NN            )
! OMP+PRIVATE( Hg0_tot, RELERR0, ABSERR0                        )
! OMP+PRIVATE( Hg2_tot, RELERR2, ABSERR2                        )
! OMP+PRIVATE( HgP_tot, RELERRP, ABSERRP                        )
! OMP+REDUCTION( +:     Hg0_tag, Hg2_tag, HgP_tag               )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize
         Hg0_tot = 0d0
         Hg0_tag = 0d0
         RELERR0 = 0d0
         ABSERR0 = 0d0
         Hg2_tot = 0d0
         Hg2_tag = 0d0
         RELERR2 = 0d0
         ABSERR2 = 0d0
         HgP_tot = 0d0
         Hgp_tag = 0d0
         RELERRP = 0d0
         ABSERRP = 0d0

         !--------
         ! Hg(0)
         !--------

         ! Total Hg(0)
         N       = ID_Hg0(ID_Hg_tot)
         Hg0_tot = STT(I,J,L,N)

         ! Sum of tagged Hg(0)
         DO NN = 2, N_Hg_CATS
            N       = ID_Hg0(NN) 
            Hg0_tag = Hg0_tag + STT(I,J,L,N)
         ENDDO

         ! Absolute error for Hg0
         ABSERR0 = ABS( Hg0_tot - Hg0_tag )

         ! Relative error for Hg0 (avoid div by zero)
         IF ( Hg0_tot > 0d0 ) THEN
            RELERR0 = ABS( ( Hg0_tot - Hg0_tag ) / Hg0_tot )
         ELSE
            RELERR0 = -999d0
         ENDIF

         !--------
         ! Hg(II)
         !--------

         ! Total Hg(II)
         N       = ID_Hg2(ID_Hg_tot)
         Hg2_tot = STT(I,J,L,N)

         ! Sum of tagged Hg(II)
         DO NN = 2, N_Hg_CATS
            N       = ID_Hg2(NN) 
            Hg2_tag = Hg2_tag + STT(I,J,L,N)
         ENDDO

         ! Absolute error for Hg2
         ABSERR2 = ABS( Hg2_tot - Hg2_tag )

         ! Relative error for Hg2 (avoid div by zero)
         IF ( Hg2_tot > 0d0 ) THEN
            RELERR2 = ABS( ( Hg2_tot - Hg2_tag ) / Hg2_tot )
         ELSE
            RELERR2 = -999d0
         ENDIF

         !--------
         ! HgP
         !--------

         ! Total Hg(P)
         N       = ID_HgP(ID_Hg_tot)
         HgP_tot = STT(I,J,L,N)

         ! Sum of tagged Hg(P)
         DO NN = 2, N_Hg_CATS
            N = ID_HgP(NN)
            IF ( N > 0 ) HgP_tag = HgP_tag + STT(I,J,L,N)
         ENDDO

         ! Absolute error for HgP
         ABSERRP = ABS( HgP_tot - HgP_tag )

         ! Relative error for HgP (avoid div by zero)
         IF ( HgP_tot > 0d0 ) THEN
            RELERRP = ABS( ( HgP_tot - HgP_tag ) / HgP_tot )
         ELSE
            RELERRP = -999d0
         ENDIF

         !----------------------------
         ! Hg(0) error is too large
         !----------------------------
         IF ( RELERR0 > MAX_RELERR .and. ABSERR0 > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 100 ) I, J, L, Hg0_tot, Hg0_tag, RELERR0, ABSERR0
! OMP END CRITICAL
         ENDIF

         !----------------------------
         ! Hg(0) error is too large
         !----------------------------
         IF ( RELERR2 > MAX_RELERR .and. ABSERR2 > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 110 ) I, J, L, Hg2_tot, Hg2_tag, RELERR2, ABSERR2
! OMP END CRITICAL
         ENDIF

         !----------------------------
         ! HgP error is too large
         !----------------------------
         IF ( RELERRP > MAX_RELERR .and. ABSERRP > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 120 ) I, J, L, HgP_tot, HgP_tag, RELERRP, ABSERRP
! OMP END CRITICAL
         ENDIF
      ENDDO
      ENDDO
      ENDDO
! OMP END PARALLEL DO

      ! FORMAT strings
 100  FORMAT( 'Hg0 error: ', 3i5, 4es13.6 )
 110  FORMAT( 'Hg2 error: ', 3i5, 4es13.6 )
 120  FORMAT( 'HgP error: ', 3i5, 4es13.6 )
 
      ! Stop if Hg0 and Hg2 errors are too large
      IF ( FLAG ) THEN
         CALL ERROR_STOP( 'Tagged Hg0, Hg2, HgP do not add up!', LOC )
      ENDIF

      ! Return to calling program 
      END SUBROUTINE CHECK_ATMOS_MERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_OCEAN_MERCURY( LOC )
!
!******************************************************************************
!  Subroutine CHECK_TAGGED_HG_OC tests whether tagged tracers in Hg0aq and
!  Hg2aq add properly within each grid box. (cdh, bmy, 2/24/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LOC (CHARACTER) : Name of routine where CHECK_OCEAN_MERCURY is called
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE LOGICAL_MOD,         ONLY : LSPLIT
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS

#     include "CMN_SIZE"            ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL, SAVE                :: FIRST = .TRUE.
      LOGICAL                      :: FLAG
      INTEGER                      :: I,       J
      REAL*8                       :: Hg0_tot, Hg0_tag, RELERR0, ABSERR0      
      REAL*8                       :: Hg2_tot, Hg2_tag, RELERR2, ABSERR2

      !=================================================================
      ! CHECK_OCEAN_MERCURY begins here!
      !=================================================================

      ! Set error condition flag
      FLAG = .FALSE.

      ! Loop over ocean surface boxes
! OMP PARALLEL DO
! OMP+DEFAULT( SHARED )
! OMP+PRIVATE( I, J, Hg0_tot, Hg0_tag, RELERR0, ABSERR0 ) 
! OMP+PRIVATE        Hg2_tot, Hg2_tag, RELERR2, ABSERR2 )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !--------------------------------------
         ! Relative and absolute errors for Hg0 
         !--------------------------------------
         Hg0_tot = Hg0aq(I,J,ID_Hg_tot)
         Hg0_tag = SUM( Hg0aq(I,J,2:N_Hg_CATS) )
         ABSERR0 = ABS( Hg0_tot - Hg0_tag )

         ! Avoid div by zero
         IF ( Hg0_tot > 0d0 ) THEN
            RELERR0 = ABS( ( Hg0_tot - Hg0_tag ) / Hg0_tot )
         ELSE
            RELERR0 = -999d0
         ENDIF

         !--------------------------------------
         ! Relative and absolute errors for Hg2
         !--------------------------------------
         Hg2_tot = Hg2aq(I,J,ID_Hg_tot)
         Hg2_tag = SUM( Hg2aq(I,J,2:N_Hg_CATS) )
         ABSERR2 = ABS( Hg2_tot - Hg2_tag )

         ! Avoid div by zero
         IF ( Hg2_tot > 0d0 ) THEN
            RELERR2 = ABS( ( Hg2_tot - Hg2_tag ) / Hg2_tot )
         ELSE
            RELERR2 = -999d0
         ENDIF

         !--------------------------------------
         ! Hg(0) error is too large
         !--------------------------------------
         IF ( RELERR0 > MAX_RELERR .and. ABSERR0 > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 100 ) I, J, Hg0_tot, Hg0_tag, RELERR0, ABSERR0
! OMP END CRITICAL
         ENDIF

         !--------------------------------------
         ! Hg(II) error is too large
         !--------------------------------------
         IF ( RELERR2 > MAX_RELERR .and. ABSERR2 > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 110 ) I, J, Hg2_tot, Hg2_tag, RELERR2, ABSERR2
! OMP END CRITICAL
         ENDIF
      ENDDO
      ENDDO
! OMP END PARALLEL DO

      ! FORMAT strings
 100  FORMAT( 'Hg0aq error: ', 2i5, 4es13.6 )
 110  FORMAT( 'Hg2aq error: ', 2i5, 4es13.6 )

      ! Stop if Hg0 and Hg2 errors are too large
      IF ( FLAG ) THEN
         CALL ERROR_STOP( 'Tagged Hg0aq, Hg2aq do not add up!', LOC )
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK_OCEAN_MERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_OCEAN_FLUXES( LOC )
!
!******************************************************************************
!  Subroutine CHECK_OCEAN_FLUXES tests whether the drydep and wetdep fluxes in
!  DD_Hg2 and WD_Hg2 sum together in each grid box. (cdh, bmy, 3/20/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LOC (CHARACTER) : Name of routine where CHECK_OCEAN_FLUXES is called
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE LOGICAL_MOD,         ONLY : LSPLIT
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS

#     include "CMN_SIZE"            ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL                      :: FLAG
      INTEGER                      :: I,         J
      REAL*8                       :: DD_tot,    DD_tag 
      REAL*8                       :: DD_RELERR, DD_ABSERR      
      REAL*8                       :: WD_tot,    WD_tag
      REAL*8                       :: WD_RELERR, WD_ABSERR

      !=================================================================
      ! CHECK_OCEAN_MERCURY begins here!
      !=================================================================

      ! Echo
      WRITE( 6, 100 )
 100  FORMAT( '     - In CHECK_OCEAN_FLUXES' )

      ! Set error condition flag
      FLAG = .FALSE.

      ! Loop over ocean surface boxes
! OMP PARALLEL DO
! OMP+DEFAULT( SHARED )
! OMP+PRIVATE( I, J, DD_tot, DD_tag, DD_RELERR, DD_ABSERR )
! OMP+PRIVATE(       WD_tot, WD_tag, WD_RELERR, WD_ABSERR )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------------------
         ! Absolute & relative errors for DD_Hg2
         !---------------------------------------
         DD_tot    = DD_Hg2(I,J,1)
         DD_tag    = SUM( DD_Hg2(I,J,2:N_Hg_CATS) )
         DD_ABSERR = ABS( DD_tot - DD_tag ) 

         ! Avoid div by zero
         IF ( DD_tot > 0d0 ) THEN
            DD_RELERR = ABS( ( DD_tot - DD_tag ) / DD_tot )
         ELSE
            DD_RELERR = -999d0
         ENDIF

         !---------------------------------------
         ! Absolute & relative errors for WD_Hg2
         !---------------------------------------
         WD_tot    = WD_Hg2(I,J,1)
         WD_tag    = SUM( WD_Hg2(I,J,2:N_Hg_CATS) )
         WD_ABSERR = ABS( WD_tot - WD_tag )

         ! Avoid div by zero
         IF ( WD_tot > 0d0 ) THEN
            WD_RELERR = ABS( ( WD_tot - WD_tag ) / WD_tot )
         ELSE
            WD_RELERR = -999d0
         ENDIF

         !---------------------------------------
         ! DD flux error is too large
         !---------------------------------------
         IF ( DD_RELERR > MAX_RELERR .and. DD_ABSERR > MAX_FLXERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 110 ) I, J, DD_tot, DD_tag, DD_RELERR, DD_ABSERR
! OMP END CRITICAL
         ENDIF

         !---------------------------------------
         ! WD flux error is too large
         !---------------------------------------
         IF ( WD_RELERR > MAX_RELERR .and. WD_ABSERR > MAX_FLXERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 120 ) I, J, WD_tot, WD_tag, WD_RELERR, WD_ABSERR
! OMP END CRITICAL
         ENDIF
      ENDDO
      ENDDO
! OMP END PARALLEL DO

      ! FORMAT strings
 110  FORMAT( 'DD_Hg2 flux error: ', 2i5, 4es13.6 )
 120  FORMAT( 'WD_Hg2 flux error: ', 2i5, 4es13.6 )

      ! Stop if Hg0 and Hg2 errors are too large
      IF ( FLAG ) THEN
         CALL ERROR_STOP( 'Tagged DD, WD fluxes do not add up!', LOC )
      ENDIf

      ! Return to calling program
      END SUBROUTINE CHECK_OCEAN_FLUXES

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_FLUX_OUT( FLUX, LOC )
!
!******************************************************************************
!  Subroutine CHECK_FLUX_OUT tests whether tagged quantities in FLUX sum 
!  together in each grid box. (cdh, bmy, 3/20/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FLUX (REAL*8)   : Flux array (output of OCEAN_MERCURY_FLUX)
!  (2 ) LOC (CHARACTER) : Name of routine where CHECK_FLUX_OUT is called
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE LOGICAL_MOD,         ONLY : LSPLIT
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS

#     include "CMN_SIZE"            ! Size parameters

      ! Arguments
      REAL*8,           INTENT(IN) :: FLUX(IIPAR,JJPAR,N_Hg_CATS)            
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL                      :: FLAG
      INTEGER                      :: I,          J
      REAL*8                       :: FLX_tot,    FLX_tag
      REAL*8                       :: FLX_RELERR, FLX_ABSERR

      !=================================================================
      ! CHECK_FLUX_OUT begins here!
      !=================================================================

      ! Echo
      WRITE( 6, 100 )
 100  FORMAT( '     - In CHECK_FLUX_OUT' )

      ! Set error condition flag
      FLAG = .FALSE.

      ! Loop over ocean surface boxes
! OMP PARALLEL DO
! OMP+DEFAULT( SHARED )
! OMP+PRIVATE( I, J, FLX_tot, FLX_tag, FLX_err )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !----------------------------------------
         ! Absolute & relative errors for FLX_Hg2
         !----------------------------------------
         FLX_tot    = FLUX(I,J,1)
         FLX_tag    = SUM( FLUX(I,J,2:N_Hg_CATS) )
         FLX_ABSERR = ABS( FLX_tot - FLX_tag )
         
         ! Avoid div by zero
         IF ( FLX_tot > 0d0 ) THEN
            FLX_RELERR = ABS( ( FLX_tot - FLX_tag ) / FLX_tot )
         ELSE
            FLX_RELERR = -999d0
         ENDIF

         !----------------------------
         ! Flux error is too large
         !----------------------------
         IF ( FLX_RELERR > MAX_RELERR  .and. 
     &        FLX_ABSERR > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 110 ) I, J, FLX_tot,    FLX_tag, 
     &                            FLX_RELERR, FLX_ABSERR
! OMP END CRITICAL
         ENDIF

      ENDDO
      ENDDO
! OMP END PARALLEL DO

      ! FORMAT strings
 110  FORMAT( 'FLX_Hg2 flux error: ', 2i5, 4es13.6 )
 
      ! Stop if Hg0 and Hg2 errors are too large
      IF ( FLAG ) THEN
         CALL ERROR_STOP( 'Tagged emission fluxes do not add up!', LOC )
      ENDIf

      ! Return to calling program
      END SUBROUTINE CHECK_FLUX_OUT

!------------------------------------------------------------------------------

      SUBROUTINE INIT_OCEAN_MERCURY( THIS_Hg_RST_FILE, THIS_USE_CHECKS )
!
!******************************************************************************
!  Subroutine INIT_OCEAN_MERCURY allocates and zeroes module arrays.  
!  (sas, cdh, bmy, 1/19/05, 2/27/06)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now just allocates arrays.  We have moved the reading of the ocean
!        Hg restart file to READ_OCEAN_Hg_RESTART.  Now make Hg0aq and Hg2aq
!        3-D arrays. Now pass Hg_RST_FILE and USE_CHECKS from "input_mod.f"
!        via the argument list. (cdh, sas, bmy, 2/27/06)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE TRACERID_MOD, ONLY : N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: THIS_Hg_RST_FILE
      LOGICAL,          INTENT(IN) :: THIS_USE_CHECKS

      ! Local variables
      INTEGER                      :: AS

      !=================================================================
      ! INIT_OCEAN_MERCURY begins here!
      !=================================================================

      ! Ocean Hg restart file name
      Hg_RST_FILE = THIS_Hg_RST_FILE
      
      ! Turn on error checks for tagged & total sums?
      USE_CHECKS  = THIS_USE_CHECKS
      
      ! Allocate arrays
      ALLOCATE( DD_Hg2( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DD_Hg2' )
      DD_Hg2 = 0d0

      ALLOCATE( dMLD( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'dMLD' )
      dMLD = 0d0

      ALLOCATE( Hg0aq( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Hg0aq' )
      Hg0aq = 0d0

      ALLOCATE( Hg2aq( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Hg2aq' )
      Hg2aq = 0d0

      ALLOCATE( HgC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HgC' )
      HgC = 0d0

      ALLOCATE( MLD( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MLD' )
      MLD = 0d0

      ALLOCATE( MLDav( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MLDav' )
      MLDav = 0d0

      ALLOCATE( newMLD( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'newMLD' )
      newMLD = 0d0

      ALLOCATE( NPP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NPP' )
      NPP = 0d0

      ALLOCATE( UPVEL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UPVEL' )
      UPVEL = 0d0

      ALLOCATE( WD_Hg2( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WD_Hg2' )
      WD_Hg2 = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_OCEAN_MERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_OCEAN_MERCURY
!
!******************************************************************************
!  Subroutine CLEANUP_OCEAN_MERCURY deallocates all arrays.  
!  (sas, cdh, bmy, 1/20/05, 1/9/06)
!  
!  NOTES:
!  (1 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (2 ) Now just deallocate arrays.  We have moved the writing of the Hg
!        restart file to MAKE_OCEAN_Hg_RESTART.  Now also deallocate HgC, dMLD
!        and MLDav arrays. (sas, cdh, bmy, 1/9/06)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_OCEAN_MERCURY begins here!
      !=================================================================
      IF ( ALLOCATED( DD_Hg2  ) ) DEALLOCATE( DD_Hg2  )
      IF ( ALLOCATED( dMLD    ) ) DEALLOCATE( dMLD    )
      IF ( ALLOCATED( Hg0aq   ) ) DEALLOCATE( Hg0aq   )  
      IF ( ALLOCATED( Hg2aq   ) ) DEALLOCATE( Hg2aq   )
      IF ( ALLOCATED( HgC     ) ) DEALLOCATE( HgC     )  
      IF ( ALLOCATED( MLD     ) ) DEALLOCATE( MLD     )
      IF ( ALLOCATED( MLDav   ) ) DEALLOCATE( MLDav   )
      IF ( ALLOCATED( newMLD  ) ) DEALLOCATE( newMLD  )
      IF ( ALLOCATED( NPP     ) ) DEALLOCATE( NPP     )
      IF ( ALLOCATED( UPVEL   ) ) DEALLOCATE( UPVEL   )
      IF ( ALLOCATED( WD_Hg2  ) ) DEALLOCATE( WD_Hg2  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_OCEAN_MERCURY

!------------------------------------------------------------------------------
     
      ! End of module
      END MODULE OCEAN_MERCURY_MOD
