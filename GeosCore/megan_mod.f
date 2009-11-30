! $Id: megan_mod.f,v 1.2 2009/11/30 19:57:56 ccarouge Exp $
      MODULE MEGAN_MOD
!
!******************************************************************************
!  Module MEGAN_MOD contains variables and routines specifying the 
!  algorithms that control the MEGAN inventory of biogenic emissions.
!  (dsa, tmf, bmy, 11/17/04, 11/6/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) T_DAY(R(MAXIJ,DAY_DIM)    : TS at each gridbox for last 8 3 hour periods
!  (2 ) T_15(MAXIJ,NUM_DAYS   )   : Holds average daily TS for past NUM_DAYS days
!  (3 ) T_15_AVG(MAXIJ)           : Is the average TS over past NUM_DAYS days 
!  (4 ) T_DAILY(MAXIJ)            : Daily average TS at each gridbox
!  (5 ) PARDR_DAY(R(MAXIJ,DAY_DIM): PARDR at each gridbox for last 8 3-hour periods
!  (6 ) PARDR_15(MAXIJ,NUM_DAYS)  : Holds average daily PARDR for past NUM_DAYS days
!  (7 ) PARDR_15_AVG(MAXIJ)       : Is the average PARDR over past NUM_DAYS days 
!  (8 ) PARDR_DAILY(MAXIJ)        : Daily average PARDR at each gridbox
!  (9 ) PARDF_DAY(R(MAXIJ,DAY_DIM): PARDF at each gridbox for last 8 3-hour periods
!  (10) PARDF_15(MAXIJ,NUM_DAYS)  : Holds average daily PARDF for past NUM_DAYS days
!  (11) PARDF_15_AVG(MAXIJ)       : Is the average PARDF over past NUM_DAYS days 
!  (12) PARDF_DAILY(MAXIJ)        : Daily average PARDF at each gridbox
!  (13) DAY_DIM                   : number 3h periods in a day
!  (14) NUM_DAYS                  : Number days in averaging periods
!  (15) AEF_ISOP(MAXIJ)           : Annual emission factor for isoprene
!  (16) AEF_MONOT(MAXIJ)          : Annual emission factor for (total) monoterpenes
!  (17) AEF_MBO(MAXIJ)            : Annual emission factor for methyl butenol
!  (18) AEF_OVOC(MAXIJ)           : Annual emission factor for other biogenic VOCs
!  (19) AEF_APINE(MAXIJ)          : Annual emission factor for alpha-pinene
!  (20) AEF_BPINE(MAXIJ)          : Annual emission factor for beta-pinene
!  (21) AEF_LIMON(MAXIJ)          : Annual emission factor for limonene
!  (22) AEF_SABIN(MAXIJ)          : Annual emission factor for sabine
!  (23) AEF_MYRCN(MAXIJ)          : Annual emission factor for myrcene
!  (24) AEF_CAREN(MAXIJ)          : Annual emission factor for 3-carene
!  (25) AEF_OCIMN(MAXIJ)          : Annual emission factor for ocimene
!  (26) AEF_SPARE(MAXIJ)          : Temporary array for holding monoterpenes

!  Module Routines:
!  ============================================================================
!  (1 ) GET_EMISOP_MEGAN  : Returns ISOPRENE Emission at a gridbox 
!  (2 ) GET_EMMONOT_MEGAN : Returns MONOTERPENE Emission at a gridbox 
!  (3 ) GET_EMMBO_MEGAN   : Returns METHYL BUTENOL Emission at a gridbox 
!  (4 ) GET_EMOVOC_MEGAN  : Returns other BVOC Emission at a gridbox 
!  (5 ) GET_AEF           : Reads archived Annual emission factor (AEF) 

!  (9 ) UPDATE_T_DAY      : Get TSKIN for each gridbox, update T_DAY, 
!                           call once every 3h
!  (10) UPDATE_T_15_AVG   : Update T_15 and T_15_AVG, call once a day
!  (11) INIT_MEGAN        : Allocate + Initialize module variables
!  (12) CLEANUP_MEGAN     : Deallocate module variables
!
!  GEOS-CHEM modules referenced by megan_mod.f
!  ============================================================================
!  (1 ) a3_read_mod.f     : Module for reading A-3 fields
!  (2 ) directory_mod.f   : Module w/ GEOS-CHEM data & met field dirs
!  (3 ) error_mod.f       : Module w/ I/O error and NaN check routines
!  (4 ) julday_mod.f      : Module w/ astronomical Julian date routines
!  (5 ) lai_mod.f         : Module w/ routines & arrays to read AVHRR LAI
!  (6 ) regrid_1x1_mod.f  : Module w/ routines to regrid 1x1 data  
!  (8 ) time_mod.f        : Module w/ routines for computing time & date
!
!
!
!  References:
!  ============================================================================
!  (1 ) Guenther, A., et al., A global model of natural volatile organic 
!        commpound emissions, J.Geophys. Res., 100, 8873-8892, 1995.
!  (2 ) Wang, Y., D. J. Jacob, and J. A. Logan, Global simulation of 
!        tropospheric O3-Nox-hydrocarbon chemistry: 1. Model formulation, J. 
!        Geophys. Res., 103, D9, 10713-10726, 1998.
!  (3 ) Guenther, A., B. Baugh, G. Brasseur, J. Greenberg, P. Harley, L. 
!        Klinger, D. Serca, and L. Vierling, Isoprene emission estimates and 
!        uncertanties for the Central African EXPRESSO study domain, J. 
!        Geophys. Res., 104, 30,625-30,639, 1999.
!  (4 ) Guenther, A. C., T. Pierce, B. Lamb, P. Harley, and R. Fall, Natural 
!        emissions of non-methane volatile organic compounds, carbon 
!        monoxide, and oxides of nitrogen from North America, Atmos. Environ.,
!        34, 2205-2230, 2000.
!  (5 ) Guenther, A., and C. Wiedinmyer, User's guide to Model of Emissions of 
!        Gases and Aerosols from Nature. http://cdp.ucar.edu. (Nov. 3, 2004) 
!  (6 ) Guenther, A., AEF for methyl butenol, personal commucation. (Nov, 2004)
!
!  NOTES:
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (6/2003).  Updated to 
!        latest algorithm and modified for the standard code by May Fu 
!        (11/2004).
!  (2 ) All emission are currently calculated using TS from DAO met field.
!        TS is the surface air temperature, which should be carefully 
!        distinguished from TSKIN. (tmf, 11/20/2004)
!  (3 ) In GEOS4, the TS used here are the T2M in the A3 files, read in 
!        'a3_read_mod.f'. 
!  (4 ) Bug fix: change #if block to also cover GCAP met fields (bmy, 12/6/05)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Bug fix: Skip Feb 29th if GCAP in INIT_MEGAN (phs, 9/18/07)
!  (7 ) Added routine GET_AEF_05x0666 to read hi-res AEF data for the GEOS-5
!        0.5 x 0.666 nested grid simulations (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!
      ! Error module
      USE ERROR_MOD,     ONLY : ERROR_STOP

      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "megan_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE :: DAY_DIM, NUM_DAYS
      PRIVATE :: T_DAY,     T_15,     T_15_AVG,     T_DAILY
      PRIVATE :: PARDR_DAY, PARDR_15, PARDR_15_AVG, PARDR_DAILY
      PRIVATE :: PARDF_DAY, PARDF_15, PARDF_15_AVG, PARDF_DAILY  

      ! PRIVATE module functions (mpb,2009)
      PRIVATE :: GET_GAMMA_LAI
      PRIVATE :: GET_GAMMA_LEAF_AGE
      PRIVATE :: GET_GAMMA_P
      PRIVATE :: GET_GAMMA_T_ISOP
      PRIVATE :: GET_GAMMA_T_NISOP
      PRIVATE :: GET_GAMMA_P_PECCA

      ! PUBLIC module functios (mpb,2009)
      PUBLIC  :: ACTIVITY_FACTORS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      INTEGER, PARAMETER  :: DAY_DIM        = 8
      INTEGER, PARAMETER  :: NUM_DAYS       = 10      ! (mpb,2009)
      REAL*8,  PARAMETER  :: WM2_TO_UMOLM2S = 4.766d0

      ! Some conversions factors (mpb,2009)
      REAL*8, PARAMETER   :: D2RAD =  3.14159d0 /   180.0d0
      REAL*8, PARAMETER   :: RAD2D =    180.0d0 / 3.14159d0 
      REAL*8, PARAMETER   :: PI    =  3.14159d0

      ! Past light & temperature conditions (mpb,2009)
      ! (1) Temperature at 2m (TS):
      REAL*8, ALLOCATABLE :: T_DAY(:,:,:)
      REAL*8, ALLOCATABLE :: T_15(:,:,:)
      REAL*8, ALLOCATABLE :: T_15_AVG(:,:)
      REAL*8, ALLOCATABLE :: T_DAILY(:,:)
      ! (2) PAR Direct:
      REAL*8, ALLOCATABLE :: PARDR_DAILY(:,:)
      REAL*8, ALLOCATABLE :: PARDR_DAY(:,:,:)
      REAL*8, ALLOCATABLE :: PARDR_15(:,:,:)
      REAL*8, ALLOCATABLE :: PARDR_15_AVG(:,:)
      ! (3) PAR Diffuse: 
      REAL*8, ALLOCATABLE :: PARDF_DAILY(:,:)
      REAL*8, ALLOCATABLE :: PARDF_DAY(:,:,:)
      REAL*8, ALLOCATABLE :: PARDF_15(:,:,:)
      REAL*8, ALLOCATABLE :: PARDF_15_AVG(:,:)

      ! Basal emission factors (mpb,2009)
      REAL*8, ALLOCATABLE :: AEF_ISOP(:,:)    
      REAL*8, ALLOCATABLE :: AEF_MONOT(:,:)
      REAL*8, ALLOCATABLE :: AEF_MBO(:,:)    
      REAL*8, ALLOCATABLE :: AEF_OVOC(:,:)    
      REAL*8, ALLOCATABLE :: AEF_APINE(:,:)
      REAL*8, ALLOCATABLE :: AEF_BPINE(:,:)
      REAL*8, ALLOCATABLE :: AEF_LIMON(:,:)
      REAL*8, ALLOCATABLE :: AEF_SABIN(:,:)    
      REAL*8, ALLOCATABLE :: AEF_MYRCN(:,:)
      REAL*8, ALLOCATABLE :: AEF_CAREN(:,:)    
      REAL*8, ALLOCATABLE :: AEF_OCIMN(:,:)      
      REAL*8, ALLOCATABLE :: AEF_SPARE(:,:) 

      ! Path to MEGAN emission factors
      CHARACTER(LEN=20)   :: MEGAN_SUBDIR = 'MEGAN_200909/'

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_EMISOP_MEGAN( I,  J,      SUNCOS, 
     &                           TS, Q_DIR, Q_DIFF, XNUMOL )
     &         RESULT( EMISOP )
!
!******************************************************************************
!  Subroutine GET_EMISOP_MEGAN computes ISOPRENE EMISSIONS in units of 
!  [atoms C/box] using the MEGAN inventory. (dsa, tmf, 9/03, 10/24/05)  
!  
!  Function GET_EMISOP_MEGAN is called from "emissdr.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I         (INTEGER ) : GEOS-CHEM longitude index
!  (2 ) J         (INTEGER ) : GEOS-CHEM latitude index
!  (3 ) SUNCOS    (REAL*8  ) : Cos( solar zenith angle )
!  (4 ) TS        (REAL*8  ) : Local surface air temperature [K]
!  (5 ) Q_DIR     (REAL*8  ) : Flux of direct PAR above canopy [W/m2]
!  (6 ) Q_DIFF    (REAL*8  ) : Flux of diffuser PAR above canopy [W/m2]
!  (7 ) XNUMOL    (REAL*8  ) : Number of atoms C / kg C 
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2000, 2004, 2006
!  (2 ) Wang,    et al, 1998
!  (3 ) Guenther et al, 2007, MEGAN v2.1 User mannual 
!
!  Notes:
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/20/04)
!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (3 ) Restructing of function & implementation of activity factors (mpb,2009)
!
!******************************************************************************
!
      ! References to F90 modules
      USE LAI_MOD,     ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M
      USE LOGICAL_MOD, ONLY : LPECCA 

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I,      J
      REAL*8,  INTENT(IN) :: SUNCOS, TS, XNUMOL, Q_DIR, Q_DIFF

      ! Local variables
      INTEGER             :: IJLOOP
      REAL*8              :: GAMMA_LAI
      REAL*8              :: GAMMA_LEAF_AGE
      REAL*8              :: GAMMA_P
      REAL*8              :: GAMMA_P_PECCA
      REAL*8              :: GAMMA_T
      REAL*8              :: GAMMA_SM
      REAL*8              :: D_BTW_M, Q_DIR_2, Q_DIFF_2
            
      ! Function return value
      REAL*8              :: EMISOP

      !=================================================================
      ! GET_EMISOP_MEGAN begins here!
      !================================================================= 

      ! Initialize return value & activity factors
      EMISOP          = 0.d0
      GAMMA_T         = 0.d0
      GAMMA_LAI       = 0.d0
      GAMMA_LEAF_AGE  = 0.d0
      GAMMA_P         = 0.d0
      GAMMA_SM        = 0.d0

      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )

      ! 1-D array index
      IJLOOP   = ( (J-1) * IIPAR ) + I
      
      ! Convert Q_DIR and Q_DIFF from (W/m2) to (micromol/m2/s)
      Q_DIR_2  = Q_DIR  * WM2_TO_UMOLM2S
      Q_DIFF_2 = Q_DIFF * WM2_TO_UMOLM2S

      !---------------------------------------------------
      ! Do only during day and over continent
      ! Only interested in terrestrial biosphere (pip)
      !---------------------------------------------------
      IF ( SUNCOS > 0d0 ) THEN

         IF ( ISOLAI(I,J) * AEF_ISOP(I,J) > 0d0 ) THEN

            IF ( LPECCA ) THEN 
                
               ! Activity factor for leaf area
               GAMMA_LAI = GET_GAMMA_LAI( MISOLAI(I,J) )


               ! Activity factor for light 
               GAMMA_P = GET_GAMMA_P_PECCA( I, J, Q_DIR_2, Q_DIFF_2,
     &                       PARDR_15_AVG(I,J), PARDF_15_AVG(I,J)  )

            ELSE 

               ! Using a canopy model set this to 1
               GAMMA_LAI = 1.0d0 

               ! Activity factor for light 
               GAMMA_P = GET_GAMMA_P( ISOLAI(I,J), SUNCOS,
     &                                 Q_DIR_2, Q_DIFF_2 )

            ENDIF

            ! Activity factor for leaf age 
            GAMMA_LEAF_AGE = GET_GAMMA_LEAF_AGE( MISOLAI(I,J), 
     &                                 PMISOLAI(I,J), D_BTW_M, 
     &                                'ISOP' , T_15_AVG(I,J) )

            ! Activity factor for temperature
            GAMMA_T = GET_GAMMA_T_ISOP( TS, T_15_AVG(I,J), 
     &                                      T_DAILY(I,J) )

            ! Activity factor for soil moisture
            GAMMA_SM = 1.d0

         ELSE

            ! If it's night or over ocean, set activity factors to zero
            GAMMA_T         = 0.d0
            GAMMA_LAI       = 0.d0
            GAMMA_LEAF_AGE  = 0.d0
            GAMMA_P         = 0.d0
            GAMMA_SM        = 0.d0

         ENDIF

         ! Isoprene emission is the product of all these
         EMISOP    = AEF_ISOP(I,J) * GAMMA_LAI * GAMMA_LEAF_AGE 
     &                             * GAMMA_T   * GAMMA_P  * GAMMA_SM  

         ! Convert from [kg/box] to [atoms C/box]
         EMISOP    = EMISOP * XNUMOL

      ENDIF

      ! return to calling program
      END FUNCTION GET_EMISOP_MEGAN

!------------------------------------------------------------------------------

      FUNCTION GET_EMMBO_MEGAN( I,  J,      SUNCOS, 
     &                          TS, Q_DIR, Q_DIFF, XNUMOL ) 
     &         RESULT( EMMBO )
!
!******************************************************************************
!  Subroutine GET_EMMBO_MEGAN computes METHYLBUTENOL EMISSIONS in units of 
!  [atoms C/box] using the MEGAN inventory. (dsa, tmf, bmy, 9/03, 10/24/05)  
!
!  Function GET_EMMBO_MEGAN is called from "emissdr.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I      (INTEGER ) : GEOS-CHEM longitude index
!  (2 ) J      (INTEGER ) : GEOS-CHEM latitude index
!  (3 ) TS     (REAL*8  ) : Local surface air temperature (K)
!  (4 ) XNUMOL (REAL*8  ) : Number of atoms C / kg C 
!  (5 ) Q_DIR  (REAL*8 )  : flux of direct PAR above canopy (W m-2)
!  (6 ) Q_DIFF (REAL*8 )  : flux of diffuser PAR above canopy (W m-2)
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2000, 2004, 2006
!  (2 ) Wang,    et al, 1998
!  (3 ) Guenther et al, 2007, MEGAN v2.1 User mannual 
!
!  NOTES:
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/20/04)
!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (3 ) Restructing of function & implementation of activity factors (mpb,2009)
!
!******************************************************************************
!
      ! References to F90 modules
      USE LAI_MOD,     ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M
      USE LOGICAL_MOD, ONLY : LPECCA 

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I,      J
      REAL*8,  INTENT(IN) :: SUNCOS, TS, XNUMOL, Q_DIR, Q_DIFF

      ! Local variable
      REAL*8              :: GAMMA_LAI
      REAL*8              :: GAMMA_LEAF_AGE
      REAL*8              :: GAMMA_P
      REAL*8              :: GAMMA_T
      REAL*8              :: GAMMA_SM
      REAL*8              :: D_BTW_M, Q_DIR_2, Q_DIFF_2
      REAL*8, PARAMETER   :: BETA = 0.09

      ! Function return value
      REAL*8              :: EMMBO
            
      !=================================================================
      ! GET_EMMBO_MEGAN begins here!
      !================================================================= 

      ! Initialize return value & activity factors
      EMMBO           = 0.d0
      GAMMA_T         = 0.d0
      GAMMA_LAI       = 0.d0
      GAMMA_LEAF_AGE  = 0.d0
      GAMMA_P         = 0.d0
      GAMMA_SM        = 0.d0

      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )

      ! Convert Q_DIR and Q_DIFF from [W/m2] to [umol/m2/s]
      Q_DIR_2  = Q_DIR  * WM2_TO_UMOLM2S
      Q_DIFF_2 = Q_DIFF * WM2_TO_UMOLM2S

      !-------------------------------------------------
      ! Do only during day and over continent
      ! Only interested in terrestrial biosphere (pip)
      !-------------------------------------------------
      IF ( SUNCOS > 0d0 ) THEN

         IF ( ISOLAI(I,J) * AEF_MBO(I,J) > 0d0 ) THEN

            IF ( LPECCA ) THEN 
                
               ! Activity factor for leaf area
               GAMMA_LAI = GET_GAMMA_LAI( MISOLAI(I,J) )

               ! Activity factor for light 
               GAMMA_P = GET_GAMMA_P_PECCA( I, J, Q_DIR_2, Q_DIFF_2,
     &                        PARDR_15_AVG(I,J), PARDF_15_AVG(I,J) )
            ELSE 

               ! Using a canopy model set this to 1
               GAMMA_LAI = 1.0d0 

               ! Activity factor for light 
               GAMMA_P = GET_GAMMA_P( ISOLAI(I,J), SUNCOS,
     &                                 Q_DIR_2, Q_DIFF_2 )
 
            ENDIF

            ! Activity factor for leaf age 
            GAMMA_LEAF_AGE = GET_GAMMA_LEAF_AGE( MISOLAI(I,J), 
     &                                 PMISOLAI(I,J), D_BTW_M, 
     &                                'MBOT' , T_15_AVG(I,J) )
         
            ! Activity factor for temperature
            GAMMA_T = GET_GAMMA_T_NISOP( TS, BETA )

            ! Activity factor for soil moisture
            GAMMA_SM = 1.d0

         ELSE

            ! If it's night or over ocean, set activity factors to zero
            GAMMA_T         = 0.d0
            GAMMA_LAI       = 0.d0
            GAMMA_LEAF_AGE  = 0.d0
            GAMMA_P         = 0.d0
            GAMMA_SM        = 0.d0

         ENDIF

         ! MBO emissions in [kg/box]
         EMMBO = AEF_MBO(I,J) * GAMMA_LAI * GAMMA_LEAF_AGE 
     &                        * GAMMA_T   * GAMMA_P * GAMMA_SM

         ! Convert from [atoms C/box] to [kg/box]
         EMMBO = EMMBO * XNUMOL

      ENDIF

      ! Return to calling program
      END FUNCTION GET_EMMBO_MEGAN

!------------------------------------------------------------------------------

      FUNCTION GET_EMMONOG_MEGAN( I , J , SUNCOS , TS , Q_DIR ,
     &                            Q_DIFF , XNUMOL , MONO_SPECIES ) 
     &            RESULT( EMMONOT )
!
!******************************************************************************
!  Subroutine GET_EMMONOG_MEGAN computes generic ('G') MONOTERPENE EMISSIONS for 
!  individual monoterpene species in units of [atoms C/box] using the new v2.1 
!  MEGAN inventory emission factor maps. (mpb,2008)  
!
!  Function GET_EMMONOG_MEGAN is called from "emissdr.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I      (INTEGER )               : GEOS-CHEM longitude index
!  (2 ) J      (INTEGER )               : GEOS-CHEM latitude index
!  (3 ) SUNCOS (REAL*8  )               : Cos( solar zenith angle )
!  (4 ) TS     (REAL*8  )               : Local surface air temperature (K)
!  (5 ) Q_DIR  (REAL*8  )               : flux of direct PAR above canopy (W m-2)
!  (6 ) Q_DIFF (REAL*8  )               : flux of diffuser PAR above canopy (W m-2) 
!  (7 ) XNUMOL (REAL*8  )               : Number of atoms C / kg C 
!  (8 ) MONO_SPECIES (CHARACTER(LEN=5)) : Monoterpene species
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2004, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 User Manual
!
!  Notes:
!  (1 ) Written by Michael Barkley (2008), based on old monoterpene code by dsa,tmf.
!  (2 ) Uses gamma factors instead of exchange factors, this includes
!        calling of a new temperature algorithm which use a beta factor.(mpb,2008)
!******************************************************************************
!
      ! References to F90 modules
      USE LAI_MOD,     ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M
      USE LOGICAL_MOD, ONLY : LPECCA 

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)           :: I,  J
      REAL*8,  INTENT(IN)           :: TS, XNUMOL, Q_DIR, Q_DIFF
      CHARACTER(LEN=5), INTENT(IN)  :: MONO_SPECIES
      REAL*8,  INTENT(IN)           :: SUNCOS

      ! Local variable
      INTEGER             :: IJLOOP
      REAL*8              :: GAMMA_LAI
      REAL*8              :: GAMMA_LEAF_AGE
      REAL*8              :: GAMMA_P
      REAL*8              :: GAMMA_T
      REAL*8              :: GAMMA_SM
      REAL*8              :: D_BTW_M
      REAL*8              :: LDF 
      REAL*8,PARAMETER    :: BETA = 0.09
      REAL*8              :: Q_DIR_2, Q_DIFF_2

      ! Function return value
      REAL*8              :: EMMONOT

      !=================================================================
      ! GET_EMMONOT_MEGAN begins here!
      !================================================================= 

      ! Initialize return value & activity factors
      EMMONOT         = 0.d0
      GAMMA_T         = 0.d0
      GAMMA_LAI       = 0.d0
      GAMMA_LEAF_AGE  = 0.d0
      GAMMA_P         = 0.d0
      GAMMA_SM        = 0.d0
      AEF_SPARE       = 0.d0

      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )

      ! Convert Q_DIR and Q_DIFF from [W/m2] to [umol/m2/s]
      Q_DIR_2  = Q_DIR  * WM2_TO_UMOLM2S
      Q_DIFF_2 = Q_DIFF * WM2_TO_UMOLM2S

      ! Need to determine which monoterpene AEFs 
      ! we need to use (mpb,2009)
      SELECT CASE( MONO_SPECIES) 
      CASE( 'APINE' )
         AEF_SPARE = AEF_APINE
         LDF       = 0.1
      CASE( 'BPINE' )
         AEF_SPARE = AEF_BPINE
         LDF       = 0.1
      CASE( 'LIMON' )
         AEF_SPARE = AEF_LIMON
         LDF       = 0.05
      CASE( 'SABIN' )
         AEF_SPARE = AEF_SABIN
         LDF       = 0.1
      CASE( 'MYRCN' )
         AEF_SPARE = AEF_MYRCN
         LDF       = 0.05
      CASE( 'CAREN' )
         AEF_SPARE = AEF_CAREN
         LDF       = 0.05
      CASE( 'OCIMN' )
         AEF_SPARE = AEF_OCIMN
         LDF       = 0.8
         CASE DEFAULT
           CALL ERROR_STOP( 'Invalid MONOTERPENE species', 
     &                      'GET_EMMONOG_MEGAN (megan_mod.f)' )
      END SELECT   

      !-----------------------------------------------------
      ! Only interested in terrestrial biosphere (pip)
      ! If (local LAI != 0 .AND. baseline emission !=0 ) 
      !-----------------------------------------------------
      IF ( ISOLAI(I,J) * AEF_SPARE(I,J) > 0d0 ) THEN

            ! Calculate gamma PAR only if sunlight conditions
            IF ( SUNCOS > 0d0 ) THEN

               IF ( LPECCA ) THEN 
               
                  GAMMA_P = GET_GAMMA_P_PECCA(I, J, Q_DIR_2, Q_DIFF_2,
     &                         PARDR_15_AVG(I,J),  PARDF_15_AVG(I,J) )
               ELSE 
  
                  GAMMA_P = GET_GAMMA_P( ISOLAI(I,J), SUNCOS,
     &                                    Q_DIR_2, Q_DIFF_2 )
               END IF 

            ELSE             
 
               ! If night
               GAMMA_P = 0.d0  
           
            END IF

            ! Activity factor for leaf area
            GAMMA_LAI = GET_GAMMA_LAI( MISOLAI(I,J) )

            ! Activity factor for leaf age 
            GAMMA_LEAF_AGE = GET_GAMMA_LEAF_AGE( MISOLAI(I,J), 
     &                                 PMISOLAI(I,J), D_BTW_M, 
     &                                 'MONO', T_15_AVG(I,J) )

            ! Activity factor for temperature
            GAMMA_T = GET_GAMMA_T_NISOP( TS, BETA )

            ! Activity factor for soil moisture
            GAMMA_SM = 1.d0

         ELSE

            ! set activity factors to zero
            GAMMA_T         = 0.d0
            GAMMA_LAI       = 0.d0
            GAMMA_LEAF_AGE  = 0.d0
            GAMMA_P         = 0.d0
            GAMMA_SM        = 0.d0

      END IF
    
      ! Monoterpene emission is the product of all these; must be 
      ! careful to distinguish between canopy & PECCA models.
      IF ( LPECCA ) THEN

         EMMONOT    = AEF_SPARE(I,J) * GAMMA_LEAF_AGE * GAMMA_T 
     &                               * GAMMA_SM       * GAMMA_LAI
     &                       * ( (1.d0 - LDF) + (LDF * GAMMA_P) )

      ELSE 

         EMMONOT    = AEF_SPARE(I,J) * GAMMA_LEAF_AGE * GAMMA_T 
     &                               * GAMMA_SM       
     &                * (  GAMMA_LAI * (1.d0 - LDF) + (LDF * GAMMA_P) )

      END IF 

      ! Convert from [kg/box] to [atoms C/box]
      EMMONOT  = EMMONOT * XNUMOL

      ! return to calling program
      END FUNCTION GET_EMMONOG_MEGAN

!------------------------------------------------------------------------------

      FUNCTION GET_EMMONOT_MEGAN( I , J , SUNCOS , TS , Q_DIR ,
     &                            Q_DIFF , XNUMOL ) 
     &            RESULT( EMMONOT )
!
!******************************************************************************
!  Subroutine GET_EMMONOT_MEGAN computes the TOTAL MONOTERPENE EMISSIONS in 
!  units of [atoms C/box] using the MEGAN v2.1 inventory. (mpb,2009)  
!
!  Function GET_EMMONOT_MEGAN is called from "emissdr.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I      (INTEGER ) : GEOS-CHEM longitude index
!  (2 ) J      (INTEGER ) : GEOS-CHEM latitude index
!  (3 ) SUNCOS (REAL*8  ) : Cos( solar zenith angle )
!  (4 ) TS     (REAL*8  ) : Local surface air temperature (K)
!  (5 ) Q_DIR  (REAL*8  ) : flux of direct PAR above canopy (W m-2)
!  (6 ) Q_DIFF (REAL*8  ) : flux of diffuser PAR above canopy (W m-2) 

!  (7 ) XNUMOL (REAL*8  ) : Number of atoms C / kg C 
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2000, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 User Manual
!
!  Notes:
!  (1 ) Original code by Michael Barkley (mpb,2009).

!******************************************************************************
!
      ! References to F90 modules
      USE LAI_MOD,    ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)           :: I,  J
      REAL*8,  INTENT(IN)           :: TS, XNUMOL, Q_DIR, Q_DIFF
      REAL*8,  INTENT(IN)           :: SUNCOS

      ! Function return value
      REAL*8              :: EMMONOT

      ! Local variables
      REAL*8                        :: MONO
      INTEGER                       :: K 
      INTEGER, PARAMETER            :: N = 7
      CHARACTER(LEN=5),DIMENSION(7) :: SPECIES      
      SPECIES = ( / 'APINE' , 'BPINE' , 'LIMON' , 'SABIN' , 'MYRCN' ,
     &            'CAREN' , 'OCIMN' / )

      !=================================================================
      ! GET_EMMONOT_MEGAN begins here!
      !================================================================= 

      ! Initialize
      EMMONOT = 0.d0

      DO K = 1 , N

         MONO = GET_EMMONOG_MEGAN( I , J , TS , SUNCOS , Q_DIR ,
     &                             Q_DIFF , XNUMOL , SPECIES(K) ) 

         EMMONOT = EMMONOT + MONO

      END DO

      ! Return to calling program
      END FUNCTION GET_EMMONOT_MEGAN

!------------------------------------------------------------------------------

      SUBROUTINE ACTIVITY_FACTORS( I , J , TS , SUNCOS , Q_DIR  ,
     &                             Q_DIFF , XNUMOL , SPECIES    , 
     &                             GAMMA_LAI , GAMMA_LEAF_AGE   , 
     &                             GAMMA_P , GAMMA_T , GAMMA_SM )

!
!******************************************************************************
!  Subroutine ACTIVITY_FACTORS computes the gamma activity factors which adjust
!  the emission factors to the current weather & vegetation conditions. 
!  Here they are calculated by (default) for isoprene. 
!
!  Function ACTIVITY_FACTORS is called from "emissdr.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I         (INTEGER )    : GEOS-CHEM longitude index
!  (2 ) J         (INTEGER )    : GEOS-CHEM latitude index
!  (3 ) TS        (REAL*8  )    : Local surface air temperature [K]
!  (4 ) SUNCOS    (REAL*8  )    : Cos( solar zenith angle )
!  (5 ) Q_DIR     (REAL*8  )    : Flux of direct PAR above canopy [W/m2]
!  (6 ) Q_DIFF    (REAL*8  )    : Flux of diffuser PAR above canopy [W/m2]
!  (7 ) XNUMOL    (REAL*8  )    : Number of atoms C / kg C 
!  (8 ) SPECIES   (CHAR, LEN=4) : Species (ISOP,MONO,MBOT); not used at present.
!
!  Notes:
!  (1 ) Original code written by Michael Barkley (mpb,2009).
! 
!******************************************************************************

      ! References to F90 modules
      USE LAI_MOD,     ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M
      USE LOGICAL_MOD, ONLY : LPECCA 

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)           :: I,      J
      REAL*8,  INTENT(IN)           :: SUNCOS, TS, XNUMOL, Q_DIR, Q_DIFF
      CHARACTER(LEN=4), INTENT(IN)  :: SPECIES

      ! Local variables
      INTEGER             :: IJLOOP
      REAL*8,  INTENT(OUT):: GAMMA_LAI
      REAL*8,  INTENT(OUT):: GAMMA_LEAF_AGE
      REAL*8,  INTENT(OUT):: GAMMA_P
      REAL*8,  INTENT(OUT):: GAMMA_T
      REAL*8,  INTENT(OUT):: GAMMA_SM
      REAL*8              :: D_BTW_M, Q_DIR_2, Q_DIFF_2
      REAL*8, PARAMETER   :: BETA = 0.09

      !=================================================================
      ! ACTIVITY_FACTORS begins here!
      !================================================================= 

      ! Initialize
      GAMMA_T         = 0.d0
      GAMMA_LAI       = 0.d0
      GAMMA_LEAF_AGE  = 0.d0
      GAMMA_P         = 0.d0
      GAMMA_SM        = 0.d0

      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )

      ! 1-D array index
      IJLOOP   = ( (J-1) * IIPAR ) + I
      
      ! Convert Q_DIR and Q_DIFF from (W/m2) to (micromol/m2/s)
      Q_DIR_2  = Q_DIR  * WM2_TO_UMOLM2S
      Q_DIFF_2 = Q_DIFF * WM2_TO_UMOLM2S

      !---------------------------------------------------
      ! Do only during day and over continent
      ! Only interested in terrestrial biosphere (pip)
      !---------------------------------------------------
      IF ( SUNCOS > 0d0 ) THEN

         IF ( ISOLAI(I,J) * AEF_ISOP(I,J) > 0d0 ) THEN

            IF ( LPECCA ) THEN 
                
               ! Activity factor for leaf area
               GAMMA_LAI = GET_GAMMA_LAI( MISOLAI(I,J) )

               ! Activity factor for light 
               GAMMA_P = GET_GAMMA_P_PECCA( I, J, Q_DIR_2, Q_DIFF_2,
     &                        PARDR_15_AVG(I,J), PARDF_15_AVG(I,J) )
            ELSE 

               ! Using a canopy model set this to 1
               GAMMA_LAI = 1.0d0 

               ! Activity factor for light 
               GAMMA_P = GET_GAMMA_P( ISOLAI(I,J), SUNCOS,
     &                                 Q_DIR_2, Q_DIFF_2 )
 
            ENDIF

            ! Activity factor for leaf age 
            GAMMA_LEAF_AGE = GET_GAMMA_LEAF_AGE( MISOLAI(I,J), 
     &                                 PMISOLAI(I,J), D_BTW_M, 
     &                                'ISOP' , T_15_AVG(I,J) )
         
            ! Activity factor for temperature
            GAMMA_T = GET_GAMMA_T_ISOP( TS, T_15_AVG(I,J), 
     &                                      T_DAILY(I,J) )

            ! Activity factor for soil moisture
            GAMMA_SM = 1.d0

         ELSE

            ! If it's night or over ocean, set activity factors to zero
            GAMMA_T         = 0.d0
            GAMMA_LAI       = 0.d0
            GAMMA_LEAF_AGE  = 0.d0
            GAMMA_P         = 0.d0
            GAMMA_SM        = 0.d0

         ENDIF

      ENDIF

      ! return to calling program
      END SUBROUTINE  ACTIVITY_FACTORS

!------------------------------------------------------------------------------

      FUNCTION GET_GAMMA_P_PECCA( I , J , Q_DIR_2, Q_DIFF_2 ,
     &                            PARDR_AVG_SIM ,  PARDF_AVG_SIM  )
     &              RESULT( GAMMA_P_PECCA )

!
!******************************************************************************
!  Computes the PECCA gamma activity factor with sensitivity to LIGHT 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I             (INTEGER) : Index of longitude box
!  (2 ) J             (INTEGER) : Index of latitude box
!  (3 ) Q_DIR_2       (REAL*8 ) : flux of direct PAR above canopy [umol m-2 s-1]
!  (4 ) Q_DIFF_2      (REAL*8 ) : flux of diffuse PAR above canopy [umol m-2 s-1]
!  (5 ) PARDR_AVG_SIM (REAL*8 ) : Avg. flux of direct PAR above canopy [W/m2]
!  (6 ) PARDR_AVG_SIM (REAL*8 ) : Avg. flux of diffuse PAR above canopy [W/m2]
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 user guide
!
!  NOTES:
!  (1 ) Here PAR*_AVG_SIM is the average light conditions over the simulation 
!       period. I've set this = 10 days to be consistent with temperature & as 
!       outlined in Guenther et al, 2006. (mpb,2009)
!  (2 ) Code was taken & adapted directly from the MEGAN v2.1 source code.
!       (mpb,2009)
!
!******************************************************************************
!
      USE TIME_MOD,   ONLY : GET_DAY_OF_YEAR
      USE TIME_MOD,   ONLY : GET_LOCALTIME
      USE GRID_MOD,   ONLY : GET_YMID  

#     include "CMN_GCTM"   ! Physical constants (why?!)

      ! Arguments as input
      INTEGER, INTENT(IN) :: I , J
      REAL*8,  INTENT(IN) :: PARDR_AVG_SIM
      REAL*8,  INTENT(IN) :: PARDF_AVG_SIM
      REAL*8,  INTENT(IN) :: Q_DIR_2, Q_DIFF_2

      ! Return value
      REAL*8              :: GAMMA_P_PECCA

      ! Local Variables
      REAL*8              :: LUT , LAT 
      REAL*8              :: mmPARDR_DAILY
      REAL*8              :: mmPARDF_DAILY
      REAL*8              :: PAC_DAILY, PAC_INSTANT, C_PPFD
      REAL*8              :: PTOA, PHI
      REAL*8              :: BETA,   SINbeta 
      INTEGER             :: DOY 
      REAL*8              :: AAA, BBB

      !-----------------------------------------------------------------
      ! Compute GAMMA_P_PECCA
      !-----------------------------------------------------------------  

      ! Initialize
      C_PPFD   = 0.0d0
      PTOA     = 0.0d0

      ! Convert past light conditions to micromol/m2/s 
      mmPARDR_DAILY   = PARDR_AVG_SIM  * WM2_TO_UMOLM2S
      mmPARDF_DAILY   = PARDF_AVG_SIM  * WM2_TO_UMOLM2S

      ! Work out the light at the top of the canopy.
      PAC_DAILY    = mmPARDR_DAILY + mmPARDF_DAILY
      PAC_INSTANT  = Q_DIR_2       +  Q_DIFF_2

      ! Get day of year, local-time and latitude
      DOY   = GET_DAY_OF_YEAR()
      LUT   = GET_LOCALTIME( I )
      LAT   = GET_YMID( J )

      ! Get solar elevation angle
      SINbeta      =  SOLAR_ANGLE( DOY , LUT , LAT  )
      BETA         =  ASIN( SINbeta ) * RAD2D       

      IF ( SINbeta .LE. 0.0d0 ) THEN

         GAMMA_P_PECCA = 0.0d0

      ELSEIF ( SINbeta .GT. 0.0d0 ) THEN       

         ! PPFD at top of atmosphere
         PTOA    = 3000.0d0 + 99.0d0 * 
     &             COS( 2.d0 * 3.14d0 *( DOY - 10 ) / 365 )

         ! Above canopy transmission
         PHI     = PAC_INSTANT / ( SINbeta * PTOA )

         ! Work out gamma P
         BBB     = 1.0d0 + 0.0005d0 *( PAC_DAILY - 400.0d0  ) 
         AAA     = ( 2.46d0 * BBB * PHI ) - ( 0.9d0 * PHI**2 )

         GAMMA_P_PECCA = SINbeta * AAA

      ENDIF

       ! Screen unforced errors. IF solar elevation angle is 
       ! less than 1 THEN gamma_p can not be greater than 0.1.
       IF ( BETA .LT. 1.0 .AND. GAMMA_P_PECCA .GT. 0.1) THEN
          GAMMA_P_PECCA  = 0.0
       ENDIF
   
      ! Prevent negative values
      GAMMA_P_PECCA = MAX( GAMMA_P_PECCA , 0d0 )

      ! return to calling program
      END FUNCTION GET_GAMMA_P_PECCA

!------------------------------------------------------------------------------

      FUNCTION SOLAR_ANGLE( DOY , SHOUR , LAT )
     &            RESULT( SINbeta )

!******************************************************************************
!  Computes the local solar angle for a given day of year, latitude and 
!  longitude (or local time) (mpb, 2009)
!
!  Function SOLAR_ANGLE is called from GAMMA_P_PECCA of "megan_mod.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DOY   (INTEGER )   : Day of year
!  (2 ) SHOUR (REAL*8  )   : Local time at longitude
!  (3 ) LAT   (REAL*8  )   : Latitude
!
!  ============================================================================
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN v2.1 user mannual 2007-09
!
!  Notes:
!  (1 ) This code was taken directly from the MEGAN v2.1 source code.(mpb, 2009)
!
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)           :: DOY 
      REAL*8,  INTENT(IN)           :: SHOUR
      REAL*8,  INTENT(IN)           :: LAT 

      ! Function return value
      REAL*8    :: SINbeta

      ! Local variables
      REAL*8    :: BETA                 ! solar elevation angle
      REAL*8    :: sindelta, cosdelta, A, B

      ! Calculation of sin beta 
      sindelta = -SIN( 0.40907d0 ) * 
     &            COS( 6.28d0 * ( DOY + 10 ) / 365 )

      cosdelta = (1-sindelta**2.)**0.5

      A = SIN( LAT * D2RAD ) * sindelta
      B = COS( LAT * D2RAD ) * cosdelta

      SINbeta = A + B * COS( 2.0d0 * PI * ( SHOUR-12 ) / 24 ) 

      END FUNCTION SOLAR_ANGLE 

!------------------------------------------------------------------------------

      FUNCTION GET_GAMMA_T_ISOP( T, PT_15 , PT_1 ) RESULT( GAMMA_T )
!
!******************************************************************************
!  Computes the temperature sensitivity for ISOPRENE ONLY (mpb, 2009)
!
!  Function GET_HEA_T is called from GET_EMISOP_MEGAN, GET_EMMBO_MEGAN 
!  and GET_EMMONOG_MEGAN of "megan_mod.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) T     (REAL*8 ) : Current leaf temperature, the surface air 
!                          temperature field (TS) is assumed equivalent to 
!                          the leaf temperature over forests.
!  (2 ) PT_15 (REAL*8 ) : Average leaf temperature over the past 15 days
!  (3 ) PT_1  (REAL*8 ) : Average leaf temperature over the past arbitray day(s).
!                         This is not used at present (but might be soon!).
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995
!  (2 ) Guenther et al, 2006
!  (3 ) Guenther et al, MEGAN v2.1 user mannual 2007-08
!
!  Notes:
!  (1 ) Includes the latest MEGAN v2.1 temperature algorithm (mpb, 2009).
!       Note, this temp-dependence is the same for the PECCA & hybrid models.
!******************************************************************************
!
      ! Arguments
      REAL*8,  INTENT(IN) :: T, PT_15 , PT_1

      ! Local Variables
      REAL*8              :: C_T,   CT1,   CT2
      REAL*8              :: E_OPT, T_OPT, X

      ! Ideal gas constant [J/mol/K]
      REAL*8, PARAMETER   :: R   = 8.314d-3

      ! Function return value
      REAL*8              :: GAMMA_T

      ! Normalization Factor 
      REAL*8, PARAMETER   ::  NORM = 0.99558166894501043d0    
                                   
      !=================================================================!
      !             ALWAYS CHECK THE ABOVE NORMALIZATION                !
      !=================================================================!
       E_OPT = 1.75d0 * EXP( 0.08d0 * ( PT_15  - 2.97d2 ) )     
       T_OPT = 3.13d2 + ( 6.0d-1 * ( PT_15 - 2.97d2 ) )
       CT1   = 80d0
       CT2   = 200d0

      ! Variable related to temperature 
      X     = ( 1.d0/T_OPT - 1.d0/T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including 
      ! effect of average temperature over previous 15 days, based on 
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) / 
     &        ( CT2 - CT1 * ( 1.d0 - EXP( CT2 * X ) ) )

      ! Hourly emission activity = C_T
      ! Prevent negative values
      GAMMA_T = MAX( C_T * NORM , 0d0 )

      ! Return to calling program
      END FUNCTION GET_GAMMA_T_ISOP

!------------------------------------------------------------------------------

      FUNCTION GET_GAMMA_T_NISOP( T , BETA ) RESULT( GAMMA_T )
!
!******************************************************************************
!  Computes the temperature activity factor (GAMMA_T) for BVOCs OTHER than 
!  isoprene. (mpb,2008)
!    
!  GAMMA_T =  exp[BETA*(T-Ts)]
!
!             where BETA   = temperature dependent parameter
!                   Ts     = standard temperature (normally 303K, 30C)
!
!  Function GET_GAMMA_T_NISOP is called by GET_EMMONOG_MEGAN and GET_EMMBO_MEGAN
!  of "megan_mod.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) T     (REAL*8 ) : Current leaf temperature, the surface air 
!                          temperature field (TS) is assumed equivalent to 
!                          the leaf temperature over forests.
!  (2 ) BETA  (REAL*8)  : Temperature factor per species (from MEGAN user 
!                          mannual). Beta = 0.09 for MBO and for monoterpene 
!                          sepecies (APINE, BPINE, LIMON, SABIN, MYRCN, CAREN 
!                          & OCIMN). Pass as an argument in case this changes.
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!
!  Notes:
!  (1 ) Original code by Michael Barkley (2009).
!       Note: If T = Ts  (i.e. standard conditions) then GAMMA_T = 1 
!******************************************************************************
!
      ! Arguments
      REAL*8,  INTENT(IN) :: T , BETA

      ! Standard reference temperature [K]
      REAL*8, PARAMETER   :: Ts = 303.0

      ! Function return value
      REAL*8              :: GAMMA_T

      !=================================================================
      ! GET_GAMMAT_NISOP begins here!
      !================================================================= 

      GAMMA_T = EXP( BETA * ( T - Ts ) )

      ! Return to calling program
      END FUNCTION GET_GAMMA_T_NISOP

!------------------------------------------------------------------------------

      FUNCTION GET_GAMMA_P( LAI, SUNCOS1, Q_DIR_2, Q_DIFF_2 ) 
     &              RESULT( GAMMA_P )
!
!******************************************************************************
!  *** REVAMPED FUNCTION ***
!  Computes the gamma activity factor with sensitivity to LIGHT (aka 'PAR')
!
!  C_PPFD: Effect of increasing PPFD up to a saturation point, where emission 
!          level off, based on Eq 4abc from Guenther et al. (1999)
!          In addition, a 5 layered canopy model based on Eqs 12-16 
!          from Guenther et al. (1995) is included to correct for light 
!          attenuation in the canopy.
!  
!  Function GET_GAMMA_P is called by the functions GET_EMISOP_MEGAN, GET_EMMBO_MEGAN 
!  GET_EMMONOG_MEGAN of "megan_mod.f"
!
!  Arguments as Input:
!  ============================================================================
!  (3 ) LAI  (REAL*8 )    : cumulative leaf area index above leaf
!  (4 ) SUNCOS1 (REAL*8 ) : Cosine of solar zenith angle
!  (5 ) Q_DIR_2 (REAL*8 ) : flux of direct PAR above canopy [umol m-2 s-1]
!  (6 ) Q_DIFF_2(REAL*8 ) : flux of diffuser PAR above canopy [umol m-2 s-1]
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995
!  (2 ) Wang     et al, 1998
!  (3 ) Guenther et al, 1999
!  (5 ) Guenther et al, 2004
!
!  NOTES:
!  (1 ) Original code by Dorian Abbot and by May Fu.
!  (2 ) This code was extracted from the previous GET_HEA_TL function. (mpb,2009)
!******************************************************************************
!
#     include "CMN_GCTM"   ! Physical constants

      ! Arguments
      REAL*8,  INTENT(IN) :: LAI, SUNCOS1, Q_DIR_2, Q_DIFF_2

      ! Local Variables

      !-----------------------------------------------------------------
      ! Canopy model variables (Eqs 12-16 from Guenther et al, 1995)
      !-----------------------------------------------------------------

      ! Return value
      REAL*8              :: GAMMA_P

      ! C_PPFD: Effect of increasing PPFD up to a saturation point, where 
      ! emissions level off, based on Eq 4abc from Guenther et al. (1999)
      ! In addition, a 5 layered canopy model based on Eqs 12-16 
      ! from Guenther et al. (1995) is included to correct for light 
      ! attenuation in the canopy.
      REAL*8              :: C_PPFD

      ! LAYERS: number of layers in canopy model
      INTEGER, PARAMETER  :: LAYERS = 5

      ! A: mean leaf-Sun angle, 60 degrees represents a 
      !    spherical leaf angle distribution
      REAL*8,  PARAMETER  :: A = 6.0d1 * PI / 1.8d2 

      ! WEIGHTGAUSS : weights for gaussian integration
      REAL*8,  PARAMETER  :: WEIGHTGAUSS(LAYERS) = (/ 1.1846d-1, 
     &                                                2.3931d-1, 
     &                                                2.8444d-1, 
     &                                                2.3931d-1, 
     &                                                1.1846d-1 /)

      ! DISTGAUSS: points to evaluate fcn for gaussian integration
      REAL*8,  PARAMETER  :: DISTGAUSS(LAYERS)   = (/ 4.6910d-2, 
     &                                                2.3075d-1, 
     &                                                5.0000d-1, 
     &                                                7.6924d-1, 
     &                                                9.5308d-1 /)

      ! SCAT: Scattering coefficient
      REAL*8, PARAMETER   :: SCAT = 2.0d-1

      ! REFLD: Reflection coefficient for diffuse light
      REAL*8, PARAMETER   :: REFLD = 5.7d-2

      ! CLUSTER: clustering coefficient (accounts for leaf clumping 
      ! influence on mean projected leaf area in the direction 
      ! of the sun's beam) use 0.85 for default
      REAL*8, PARAMETER   :: CLUSTER = 8.5d-1

      ! NORMAL_FACTOR : C_PPFD calculated with LAI = 5, and above canopy 
      ! total PPFD = 1500 umol/m2/s, and Q_DIFF = 0.2 * above canopy total 
      ! PPFD.  May Fu calculated this to be 1.8967d0
      ! Quickly checked by mpb in 2008 & found to be:
      REAL*8, PARAMETER   :: NORMAL_FACTOR = 1.80437817285271950d0

      ! F_SUN: Fraction leaves sunlit
      REAL*8              :: F_SUN

      ! F_SHADE: Fraction leaves shaded
      REAL*8              :: F_SHADE

      ! LAI_DEPTH: Cumulative LAI above current layer
      REAL*8              :: LAI_DEPTH

      ! Q_SUN: Flux density of PAR on sunny leaves [umol/m2/s]
      REAL*8              :: Q_SUN
 
      ! Q_SHADE: Flux density of PAR on shaded leaves [umol/m2/s]
      REAL*8              :: Q_SHADE, Q_SHADE_1, Q_SHADE_2

      ! KB: Extinction coefficient for black leaves for direct (beam)
      REAL*8              :: KB, KBP

      ! KD: Extinction coefficient for black leaves for diffuse light
      REAL*8              :: KD,   KDP
      REAL*8              :: REFLB

      ! C_P_SUN: C_PPFD at layer I for sunlit leaves
      REAL*8              :: C_P_SUN
      
      ! C_P_SHADE: C_PPFD at layer I for shaded leaves
      REAL*8              :: C_P_SHADE

      ! C_PPFD_I    : C_PPFD at layer I
      REAL*8              :: C_PPFD_I

      ! Empirical functions (Eq 4a, 4b, 4c from Guenther et al, 1999)
      REAL*8              :: ALPHA, CL

      ! ??
      REAL*8              :: P

      ! Index for layers
      INTEGER             :: I

      !-----------------------------------------------------------------
      ! Compute C_PPFD
      !-----------------------------------------------------------------  

      ! Initialize
      C_PPFD  = 0.d0

      ! 0.5 and 0.8 assume a spherical leaf-angle distribution
      KB      = 0.5d0 * CLUSTER / SUNCOS1
      KD      = 0.8d0 * CLUSTER
      P       = SQRT( 1.d0 - SCAT )
      REFLB   = 1.d0 - 
     &          EXP(-2.d0*KB * ( 1.d0-P ) / ( 1.d0+P ) / ( 1.d0+KB ) )
      KBP     = KB * P
      KDP     = KD * P

      ! 5-layer Gaussian integration over canopy
      DO I = 1, LAYERS 

         ! Cumulative LAI above layer I
         LAI_DEPTH  = LAI * DISTGAUSS( I )

         ! Fraction sun and shade leaves at layer I
         F_SUN      = EXP( -KB * LAI_DEPTH ) 
         F_SHADE    = 1.d0 - F_SUN

         ! For PAR on shaded leaves
         Q_SHADE_1  = Q_DIFF_2 * KDP * ( 1.d0 - REFLD ) *
     &                EXP( -KDP * LAI_DEPTH ) 

         ! For PAR on shaded leaves
         Q_SHADE_2  = Q_DIR_2 * ( KBP * ( 1.d0 - REFLB ) * 
     &                EXP( -KBP * LAI_DEPTH ) - KB * ( 1.d0 - SCAT ) * 
     &                EXP( -KB * LAI_DEPTH ) )

         ! PAR on shaded leaves
         Q_SHADE    = ( Q_SHADE_1 + Q_SHADE_2 ) / ( 1.d0 - SCAT )

         ! PAR on sunlit leaves
         Q_SUN      = Q_SHADE + KB * Q_DIR_2

         ! Update C_P_SUN and C_P_SHADE at layer I
         ! The following already accounts for canopy attenuation
         ! (Guenther et al, 1999)
         ALPHA      = 1.0d-3 + ( 8.5d-4 ) * LAI_DEPTH
         CL         = ( 1.42d0 ) * EXP( -( 3.0d-1 ) * LAI_DEPTH )

         C_P_SUN    = ALPHA * CL * Q_SUN / 
     &                SQRT( 1.d0 + ALPHA*ALPHA * Q_SUN*Q_SUN )

         C_P_SHADE  = ALPHA * CL * Q_SHADE / 
     &                SQRT( 1.d0 + ALPHA*ALPHA * Q_SHADE*Q_SHADE )

         ! Update C_PPFD_I at layer I
         C_PPFD_I   = F_SUN * C_P_SUN + F_SHADE * C_P_SHADE

         ! Add on to total C_PPFD
         C_PPFD     = C_PPFD + WEIGHTGAUSS( I ) * C_PPFD_I * LAI

      ENDDO

      ! Prevent negative values.
      GAMMA_P = MAX( C_PPFD / NORMAL_FACTOR, 0d0 )

      ! return to calling program
      END FUNCTION GET_GAMMA_P

!------------------------------------------------------------------------------

      FUNCTION GET_GAMMA_LEAF_AGE( CMLAI, PMLAI, T , SPECIES , TT ) 
     &            RESULT( GAMMA_LEAF_AGE )
!
!******************************************************************************
!  Computes the gamma exchange activity factor which is sensitive to leaf
!    age (= GAMMA_LEAF_AGE).
!  
!  Function GET_GAMMA_LEAF_AGE is called from GET_EMISOP_MEGAN, GET_EMMBO_MEGAN 
!  and GET_EMMONOG_MEGAN of "megan_mod.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) CMLAI   (REAL*8 )   : Current month leaf area index at gridbox 
!  (2 ) PMLAI   (REAL*8 )   : Last month leaf area index at gridbox 
!  (3 ) T       (REAL*8 )   : Number of days between current and previous LAI.
!  (4 ) SPECIES (CHARACTER) : BVOC species 
!  (5 ) TT      (REAL*8 )   : Daily average temperature
!
!  References (see above for full citations):
!  ============================================================================
!  (3 ) Guenther et al, 2006
!  (5 ) Guenther et al, MEGAN user mannual 2007-08
!
!
!  NOTES:
!  (1 ) Original code by Dorian Abbot (9/2003). Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) Update to publically released (as of 11/2004) MEGAN algorithm and 
!        modified for the standard code by May Fu (11/2004).
!  (3 ) Algorithm is based on the latest User's Guide (tmf, 11/19/04)
!  (4 ) Renamed & now includes specific relative emission activity factors for
!       each BVOC based on MEGAN v2.1 algorithm (mpb,2008)
!  (5 ) Now calculate TI (number of days after budbreak required to induce 
!       iso. em.) and TM (number of days after budbreak required to reach 
!       peak iso. em. rates) using the daily average temperature, instead 
!       of using fixed values (mpb,2008)
!       NOTE: Can create 20% increases in tropics (Guenther et al 2006)
!  (6 ) Implemented change for the calculation of FGRO if ( CMLAI > PMLAI ),
!       i.e. if LAI has increased with time, and used new values for 
!       all foilage fractions if ( CMLAI = PMLAI ). Also removed TG variable 
!       as not now needed. (mpb,2000)
!******************************************************************************

      ! Arguments
      REAL*8, INTENT(IN)            :: T
      REAL*8, INTENT(IN)            :: CMLAI, PMLAI
      CHARACTER(LEN=4), INTENT(IN)  :: SPECIES
      REAL*8, INTENT(IN)            :: TT

      ! Function return value
      REAL*8             :: GAMMA_LEAF_AGE      
          
      ! Local Variables

      ! M_AGE: leaf age factor
      REAL*8             :: M_AGE

      ! FNEW, ANEW: new foliage that emits negligible amounts of isoprene
      REAL*8             :: FNEW

      ! FGRO, AGRO: growing foliage that emits isoprene at low rates
      REAL*8             :: FGRO

      ! FMAT, AMAT: mature foliage that emits isoprene at peak rates
      REAL*8             :: FMAT

      ! FSEN, ASEN: senescing foliage that emits isoprene at reduced rates
      REAL*8             :: FSEN

      ! TI: number of days after budbreak required to induce iso. em.
      REAL*8             :: TI   

      ! TM: number of days after budbreak required to reach peak iso. em. rates
      REAL*8             :: TM

      ! Index for the relative emission activity (mpb,2008)
      INTEGER            :: AINDX
 
      ! Use species specific relative emission activity factors (mpb,2008)
      INTEGER, PARAMETER :: N_CAT = 4
      REAL*8             :: ANEW(N_CAT)
      REAL*8             :: AGRO(N_CAT)
      REAL*8             :: AMAT(N_CAT)
      REAL*8             :: ASEN(N_CAT)
      REAL*8             :: NORM_V(N_CAT)

      ! Constant Factors (not used)
      DATA    ANEW(  1),  AGRO(  1),  AMAT(  1),  ASEN(  1)
     &     /  1.0d0    ,  1.0d0    ,  1.0d0    ,  1.0d0      /
      ! Monoterpenes 
      DATA    ANEW(  2),  AGRO(  2),  AMAT(  2),  ASEN(  2)
     &     /  2.0d0   ,   1.8d0    ,  0.95d0   ,  1.0d0      /
      ! Isoprene and MBO 
      DATA    ANEW(  3),  AGRO(  3),  AMAT(  3),  ASEN(  3)
     &     /  0.05d0   ,  0.6d0    ,  1.125d0  ,  1.0d0      /
      ! Current set-up (v7-04-11) - used for other VOCs (OVOC)
      DATA    ANEW(  4),  AGRO(  4),  AMAT(  4),  ASEN(  4)
     &     /  0.01d0   ,  0.5d0    ,  1.00d0   ,  0.33d0     /

      ! Normalization factors
      DATA   NORM_V(1) , NORM_V(2) , NORM_V(3) , NORM_V(4) 
     &     /                 1.0d0 , ! Constant
     &       0.96153846153846145d0 , ! Monoterpenes
     &       0.94339622641509424d0 , ! Isoprene & MBO
     &         1.132502831257078d0 / ! Other VOCs

      !=================================================================
      ! GET_GAMMA_LEAF_AGE begins here!
      !================================================================= 
      
      !-----------------------
      ! Compute TI and TM 
      ! (mpb,2009)
      !-----------------------

      IF ( TT <= 303.d0 ) THEN
         TI = 5.0d0 + 0.7 * ( 300.0d0 - TT )
      ELSEIF ( TT >  303.d0 ) THEN
         TI = 2.9d0
      ENDIF   
      TM = 2.3d0 * TI

      !-----------------------
      ! Compute M_AGE
      !-----------------------

      IF ( CMLAI == PMLAI ) THEN !(i.e. LAI stays the same) 

         !New values           Old Vlaues (mpb, 2009)
         FMAT = 0.8d0            !1.d0  
         FNEW = 0.d0             !0.d0
         FGRO = 0.1d0            !0.d0
         FSEN = 0.1d0            !0.d0

      ELSE IF ( CMLAI > PMLAI ) THEN !(i.e. LAI has INcreased) 

         IF ( T > TI ) THEN
            FNEW = ( TI / T ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE
            FNEW = 1.d0 - ( PMLAI / CMLAI )
         ENDIF

         IF ( T > TM ) THEN
            FMAT = ( PMLAI / CMLAI ) +
     &             ( ( T - TM ) / T ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE 
            FMAT = ( PMLAI / CMLAI )
         ENDIF

         ! Implement new condition for FGRO (mpb,2009)
         FSEN = 0.d0
         FGRO = 1.d0 - FNEW - FMAT

      ELSE ! This is the case if  PMLAI > CMLAI (i.e. LAI has DEcreased) 

         FSEN = ( PMLAI - CMLAI ) / PMLAI
         FMAT = 1.d0 - FSEN
         FGRO = 0.d0
         FNEW = 0.d0

      ENDIF

      ! Choose relative emission activity (mpb,2009)
      SELECT CASE ( TRIM(SPECIES) )
      CASE ('CONSTANT')
         AINDX = 1
      CASE ('MONO')
         AINDX = 2
      CASE ( 'ISOP','MBOT' )
         AINDX = 3
      CASE ('OVOC')
         AINDX = 4
      CASE DEFAULT
           CALL ERROR_STOP( 'Invalid BVOC species', 
     &                      'GET_GAMMA_LEAF_AGE (megan_mod.f)' )
      END SELECT

      ! Age factor
      M_AGE = FNEW*ANEW(AINDX) + FGRO*AGRO(AINDX) + 
     &        FMAT*AMAT(AINDX) + FSEN*ASEN(AINDX)  

      !Note: I think this should be normalized.
      !But, in the megan code I've consistently download 
      !from NCAR it never is...

      ! Normalize & prevent negative values
      GAMMA_LEAF_AGE = MAX( M_AGE * NORM_V(AINDX) , 0d0 )

      ! return to calling program
      END FUNCTION GET_GAMMA_LEAF_AGE

!------------------------------------------------------------------------------

      FUNCTION GET_GAMMA_LAI( CMLAI ) RESULT( GAMMA_LAI )
!
!******************************************************************************
!  Computes the gamma exchange activity factor which is sensitive to leaf
!    area (= GAMMA_LAI).
!  
!  Function GET_GAMMA_LAI is called from GET_EMISOP_MEGAN, GET_EMMBO_MEGAN 
!  and GET_EMMONOG_MEGAN of "megan_mod.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) CMLAI   (REAL*8 )   : Current month leaf area index at gridbox 
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!
!  NOTES:
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) Update to publically released (as of 11/2004) MEGAN algorithm and 
!        modified for the standard code by May Fu (11/2004).
!  (3 ) Algorithm is based on the latest MEGAN v2.1 User's Guide (mpb,2009)
!******************************************************************************

      ! Arguments
      REAL*8, INTENT(IN)  :: CMLAI

      ! Return value
      REAL*8              :: GAMMA_LAI 

      !-----------------------
      ! Compute GAMMA_LAI
      !-----------------------
      GAMMA_LAI = 0.49d0 * CMLAI / SQRT( 1.d0 + 0.2d0 * CMLAI*CMLAI )

      END FUNCTION GET_GAMMA_LAI 

!------------------------------------------------------------------------------

      SUBROUTINE GET_AEF
!
!******************************************************************************
!  Subroutine GET_AEF reads Annual Emission Factor for all biogenic VOC 
!  species from disk.  Function GET_AEF is called from "main.f"
!  (tmf, bmy, 10/24/05)
!
!  Reference
!  ============================================================================
!  (5 ) Guenther et al, 2004 
!
!  Notes:
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) AEF detailed in the latest MEGAN User's Guide (tmf, 11/19/04)
!  (3 ) Bug fix (tmf, 11/30/04)
!  (4 ) Now reads 1x1 files and regrids to current resolution (bmy, 10/24/05)
!  (5 ) Uses new v2.1 emission factors maps for isoprene, MBO and 7 monoterpene
!        species, download in 2009. (mpb,2009)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_RES_EXT, READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
      USE TIME_MOD,       ONLY : GET_TS_EMIS
      USE GRID_MOD,       ONLY : GET_AREA_M2

#     include "CMN_SIZE"       ! Size parameters

      ! Local Variables
      INTEGER                 :: I, J, IJLOOP
      REAL*4                  :: ARRAY(I1x1,J1x1,1)
      REAL*8                  :: DTSRCE, AREA_M2, FACTOR
      CHARACTER(LEN=255)      :: FILENAME


      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! In MEGAN v2.1 the emission factors are expressed in units of:
      !
      ! ** micro-grams of compound per m2 per hour **
      !
      ! Previously this was:
      !
      ! ** micro-grams of CARBON per m2 per hour **    
      !
      ! We must therefore apply a conversion factor to change the new 
      ! 'compound' emissions into carbon emissions.
      REAL*8                  ::  ISOP2CARBON               
      REAL*8                  ::  MONO2CARBON
      REAL*8                  ::  MBO2CARBON
      !
      ! Note, where it is written as [ug C/m2/hr] think of C = compound
      ! (mpb,2008)
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !=================================================================
      ! GET_AEF begins here!
      !=================================================================

      ! Emission timestep [min]
      DTSRCE = GET_TS_EMIS()

      !-----------------------------------------------
      ! Determine the 'carbon conversion' factors
      ! = molar mass carbon / molar mass compound
      ! (mpb,2008)
      !-----------------------------------------------

      ISOP2CARBON =  60d0 /  68d0
      MONO2CARBON = 120d0 / 136d0
      MBO2CARBON  =  60d0 /  86d0
      
      !---------------------------------------------
      ! Read in ISOPRENE Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &   'MEGAN_200909/MEGANv2.1_AEF_ISOPRENE.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
100   FORMAT( '     - GET_AEF: Reading ', a )

      ! Read data at 1x1 resolution [ug C/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_ISOP ) 

      ! Loop over longitudes
      DO J = 1, JJPAR
       
         ! Grid box surface area [m2]
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0
         
 
         ! Loop over latitudes
         DO I = 1, IIPAR

            ! Convert AEF_ISOP To [kg C/box]
            AEF_ISOP(I,J) = AEF_ISOP(I,J) * FACTOR * ISOP2CARBON

         ENDDO
      ENDDO

      !---------------------------------------------
      ! Read in MONOTERPENE Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &  TRIM( MEGAN_SUBDIR ) // 'MEGANv2.1_AEF_MONOTERPENES.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_MONOT ) 

      ! Loop over longitudes
      DO J = 1, JJPAR

         ! Surface area 
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR  = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over latitudes
         DO I = 1, IIPAR

            ! Convert AEF_MONOT to [kg C/box]
            AEF_MONOT(I,J) = AEF_MONOT(I,J) * FACTOR * MONO2CARBON

         ENDDO
      ENDDO

      !---------------------------------------------
      ! Read in MBO Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &    TRIM( MEGAN_SUBDIR ) // 'MEGANv2.1_AEF_MBO.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_MBO ) 

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_MBO(I,J) = AEF_MBO(I,J) * FACTOR * MBO2CARBON
            
         ENDDO
      ENDDO

      !---------------------------------------------
      ! Read in other VOC Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &           'MEGAN_200510/MEGAN_AEF_MTP.geos.1x1' 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_OVOC ) 

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_OVOC(I,J) = AEF_OVOC(I,J) * FACTOR
            
         ENDDO
      ENDDO

      !--------------------------------------------------
      ! Read in other Monoterpene Annual Emission Factors
      !-------------------------------------------------

      ! 1) Alpha Pinene (APIN) 

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &   TRIM( MEGAN_SUBDIR ) // 'MEGANv2.1_AEF_ALPHA_PINENE.geos.1x1' 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_APINE ) 

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_APINE(I,J) = AEF_APINE(I,J) * FACTOR * MONO2CARBON
            
         ENDDO
      ENDDO

      ! 2) Beta Pinene (BPIN) 

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &     TRIM( MEGAN_SUBDIR ) // 'MEGANv2.1_AEF_BETA_PINENE.geos.1x1' 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_BPINE ) 

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_BPINE(I,J) = AEF_BPINE(I,J) * FACTOR * MONO2CARBON
            
         ENDDO
      ENDDO

      ! 3) Limonene (LIMON) 

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &   TRIM( MEGAN_SUBDIR ) // 'MEGANv2.1_AEF_LIMONENE.geos.1x1' 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_LIMON ) 

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_LIMON(I,J) = AEF_LIMON(I,J) * FACTOR * MONO2CARBON
            
         ENDDO
      ENDDO

      ! 4) Sabinene (SABIN) 

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &   TRIM( MEGAN_SUBDIR ) // 'MEGANv2.1_AEF_SABINENE.geos.1x1' 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_SABIN ) 

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_SABIN(I,J) = AEF_SABIN(I,J) * FACTOR * MONO2CARBON
            
         ENDDO
      ENDDO

      ! 5) Myrcene (MYRCN) 

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &   TRIM( MEGAN_SUBDIR ) // 'MEGANv2.1_AEF_MYRCENE.geos.1x1' 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_MYRCN ) 

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_MYRCN(I,J) = AEF_MYRCN(I,J) * FACTOR * MONO2CARBON
            
         ENDDO
      ENDDO

      ! 6) 3-Carene (CAREN) 

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &   TRIM( MEGAN_SUBDIR ) // 'MEGANv2.1_AEF_CARENE.geos.1x1' 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_CAREN ) 

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_CAREN(I,J) = AEF_CAREN(I,J) * FACTOR * MONO2CARBON
            
         ENDDO
      ENDDO

      ! 7) Ocimene (OCIMN) 

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &   TRIM( MEGAN_SUBDIR ) // 'MEGANv2.1_AEF_OCIMENE.geos.1x1' 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,   
     &                 0d0,       I1x1,      J1x1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_1x1( 'ug C/m2/hr', ARRAY, AEF_OCIMN ) 

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_OCIMN(I,J) = AEF_OCIMN(I,J) * FACTOR * MONO2CARBON
            
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE GET_AEF

!------------------------------------------------------------------------------

      SUBROUTINE GET_AEF_05x0666
!
!******************************************************************************
!  Subroutine GET_AEF reads Annual Emission Factor for all biogenic VOC 
!  species from disk.  Function GET_AEF is called from "main.f".  Specially
!  constructed to read 0.5 x 0.666 nested grid data for the GEOS-5 nested
!  grid simulations. (yxw, dan, bmy, 11/6/08)
!
!  Reference
!  ============================================================================
!  (5 ) Guenther et al, 2004 
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_RES_EXT, READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_05x0666
      USE TIME_MOD,       ONLY : GET_TS_EMIS
      USE GRID_MOD,       ONLY : GET_AREA_M2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR

#     include "CMN_SIZE"       ! Size parameters

      ! Local Variables
      INTEGER                 :: I, J, IJLOOP
      REAL*4                  :: ARRAY(I05x0666,J05x0666,1)
      REAL*8                  :: GEOS_05x0666(I05x0666,J05x0666,1)  !(dan)
      REAL*8                  :: DTSRCE, AREA_M2, FACTOR
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! GET_AEF begins here!
      !=================================================================

      ! Emission timestep [min]
      DTSRCE = GET_TS_EMIS()

      !---------------------------------------------
      ! Read in ISOPRENE Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR ) //
     &           'MEGAN_200510/MEGAN_AEF_ISOP.geos.05x0666'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
100   FORMAT( '     - GET_AEF: Reading ', a )

      ! Read data at 1x1 resolution [ug C/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,
     &                 0d0,       I05x0666,      J05x0666,
     &                 1,         ARRAY,     QUIET=.TRUE. )

         ! Cast to REAL*8 before regridding
         GEOS_05x0666(:,:,1) = ARRAY(:,:,1)

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_05x0666( 1, 'ug C/m2/hr', GEOS_05x0666, AEF_ISOP )

      ! Loop over longitudes
      DO J = 1, JJPAR

         ! Grid box surface area [m2]
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over latitudes
         DO I = 1, IIPAR

            ! Convert AEF_ISOP To [kg C/box]
            AEF_ISOP(I,J) = AEF_ISOP(I,J) * FACTOR

         ENDDO
      ENDDO

      !---------------------------------------------
      ! Read in MONOTERPENE Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR ) //
     &           'MEGAN_200510/MEGAN_AEF_MTP.geos.05x0666'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,
     &                 0d0,       I05x0666,      J05x0666,
     &                 1,         ARRAY,     QUIET=.TRUE. )

         ! Cast to REAL*8 before regridding
         GEOS_05x0666(:,:,1) = ARRAY(:,:,1)

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_05x0666( 1,'ug C/m2/hr', GEOS_05x0666, AEF_MONOT )

      ! Loop over longitudes
      DO J = 1, JJPAR

         ! Surface area 
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR  = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over latitudes
         DO I = 1, IIPAR

            ! Convert AEF_MONOT to [kg C/box]
            AEF_MONOT(I,J) = AEF_MONOT(I,J) * FACTOR

         ENDDO
      ENDDO

      !---------------------------------------------
      ! Read in MBO Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR ) //
     &           'MEGAN_200510/MEGAN_AEF_MBO.geos.05x0666'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,
     &                 0d0,       I05x0666,      J05x0666,
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      GEOS_05x0666(:,:,1) = ARRAY(:,:,1)


      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_05x0666( 1,'ug C/m2/hr',GEOS_05x0666, AEF_MBO )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_MBO(I,J) = AEF_MBO(I,J) * FACTOR

         ENDDO
      ENDDO

      !---------------------------------------------
      ! Read in other VOC Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR ) //
     &           'MEGAN_200510/MEGAN_AEF_MTP.geos.05x0666'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data at 1x1 resolution [ug C/m2/hr]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', 99,
     &                 0d0,       I05x0666,      J05x0666,
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      GEOS_05x0666(:,:,1) = ARRAY(:,:,1)

      ! Regrid from 1x1 to the current grid (also cast to REAL*8)
      CALL DO_REGRID_05x0666( 1,'ug C/m2/hr', GEOS_05x0666, AEF_OVOC )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid box surface area
         AREA_M2 = GET_AREA_M2( J )

         ! Conversion factor from [ug C/m2/hr] to [kg C/box]
         FACTOR = 1.d-9 / 3600.d0 * AREA_M2 * DTSRCE * 60.d0

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert AEF to [kg C/box/step]
            AEF_OVOC(I,J) = AEF_OVOC(I,J) * FACTOR

         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE GET_AEF_05x0666

!------------------------------------------------------------------------------

      SUBROUTINE UPDATE_T_DAY
!
!******************************************************************************
!  Subroutine UPDATE_T_DAY must be called every time the A-3 fields are 
!  updated. Each 3h TS value for each gridbox is moved up one spot in the 
!  matrix and the current value is put in the last spot. (dsa, 6/17/03)
!
!  NOTES:
!  (1 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (2 ) In GEOS4, TS are originally T2M in the A3 files, read in 
!        'a3_read_mod.f'. 
!******************************************************************************
!
      USE MEGANUT_MOD    ! We use all functions from the module
#     include "CMN_SIZE" ! Size parameters

      ! Local Variables
      INTEGER           :: I, J, D

!--- Moved all functions to module MEGANUT_MOD
!      ! External functions
!      REAL*8, EXTERNAL  :: XLTMMP
!      REAL*8, EXTERNAL  :: XLPARDR ! (mpb,2009)
!      REAL*8, EXTERNAL  :: XLPARDF ! (mpb,2009)


      !=================================================================
      ! UPDATE_T_DAY begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, D )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Move each day up
         DO D = DAY_DIM, 2, -1
            ! Need PAR as well as Temp  (mpb,2009)
            T_DAY(I,J,D)     = T_DAY(I,J,D-1)
            PARDR_DAY(I,J,D) = PARDR_DAY(I,J,D-1)
            PARDF_DAY(I,J,D) = PARDF_DAY(I,J,D-1)
         ENDDO
            
         ! Store 
         T_DAY(I,J,1)     = XLTMMP(I,J)
         PARDF_DAY(I,J,1) = XLPARDF(I,J)
         PARDR_DAY(I,J,1) = XLPARDR(I,J)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
          
      ! return to calling program
      END SUBROUTINE UPDATE_T_DAY

!------------------------------------------------------------------------------

      SUBROUTINE UPDATE_T_15_AVG
!
!******************************************************************************
!  Subroutine UPDATE_T_15_AVG should be called at the beginning of each day.
!  It loops through the gridboxes doing the following:
!     1. Average T_DAY over the 8 TS values during the day.  
!     2. Push the daily average TS values through T_15, throwing out the 
!        oldest and putting the newest (the T_DAY average) in the last spot 
!     3. Get T_15_AVG by averaging T_15 over the 15 day period. 
!  (dsa, tmf, bmy, 6/17/03, 10/24/05)
!
!  NOTES:
!  (1 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (2 ) In GEOS4, TS are originally T2M in the A3 files, read in 
!        'a3_read_mod.f'. 
!******************************************************************************
! 
      IMPLICIT NONE

#     include "CMN_SIZE" ! MAXIJ

      ! Local Variables
      INTEGER           :: I,     J,      D
      REAL*8            :: D_DIM, D_DAYS, TMP_T
      REAL*8            :: TMP_PARDR , TMP_PARDF ! (mpb,2009)

      !=================================================================
      ! UPDATE_T_15_AVG begins here!
      !=================================================================

      ! Convert to REAL*8
      D_DIM  = DBLE( DAY_DIM  )
      D_DAYS = DBLE( NUM_DAYS )

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, D, TMP_T, TMP_PARDR, TMP_PARDF )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Average T_DAY over the 8 TS values during the day.
         TMP_T = SUM( T_DAY(I,J,:) ) / D_DIM
         
         ! Do the same for light (mpb,2009)
         TMP_PARDR = SUM( PARDR_DAY(I,J,:) ) / D_DIM
         TMP_PARDF = SUM( PARDF_DAY(I,J,:) ) / D_DIM

         ! Push the daily average TS values through T_15,
         ! throwing out the oldest 
         DO D = NUM_DAYS, 2, -1
            T_15(I,J,D)     = T_15(I,J,D-1)
            PARDR_15(I,J,D) = PARDR_15(I,J,D-1)
            PARDF_15(I,J,D) = PARDF_15(I,J,D-1)
         ENDDO

         ! Put the newest daily average TS value in the first spot
         T_15(I,J,1) = TMP_T

         ! Get T_15_AVG by averaging T_15 over the 15 day period.
         T_15_AVG(I,J) = SUM( T_15(I,J,:) ) / D_DAYS 

         ! Assign daily average temperature to T_DAILY (mpb,2009)
         T_DAILY(I,J)  = TMP_T

         ! Repeat for PAR diffuse & direct (mpb,2009)

         PARDR_15(I,J,1)   = TMP_PARDR
         PARDR_15_AVG(I,J) = SUM( PARDR_15(I,J,:) ) / D_DAYS 
         PARDR_DAILY(I,J)  = TMP_PARDR

         PARDF_15(I,J,1)   = TMP_PARDF
         PARDF_15_AVG(I,J) = SUM( PARDF_15(I,J,:) ) / D_DAYS 
         PARDF_DAILY(I,J)  = TMP_PARDF


      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE UPDATE_T_15_AVG

!------------------------------------------------------------------------------

      SUBROUTINE INIT_MEGAN
!
!******************************************************************************
!  Subroutine INIT_MEGAN allocates and initializes the T_DAY, T_15,  
!  T_15_AVG, and AEF_* arrays. (dsa, tmf, bmy, 10/24/05, 11/6/08)
!
!  NOTES:
!  (1 ) Change the logic in the #if block for G4AHEAD. (bmy, 12/6/05)
!  (2 ) Bug fix: skip Feb 29th if GCAP (phs, 9/18/07)
!  (3 ) Now call GET_AEF_05x0666 for GEOS-5 nested grids (yxw,dan,bmy, 11/6/08)
!******************************************************************************
!
      ! References to F90 modules
      USE A3_READ_MOD
      USE FILE_MOD,    ONLY : IU_A3
      USE JULDAY_MOD,  ONLY : CALDATE
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LAI_MOD,     ONLY : INIT_LAI
      USE LOGICAL_MOD, ONLY : LUNZIP
      USE TIME_MOD,    ONLY : GET_FIRST_A3_TIME, GET_JD
      USE TIME_MOD,    ONLY : ITS_A_LEAPYEAR,    YMD_EXTRACT
      
#     include "CMN_SIZE"    ! Size parameters

      ! Local Variables
      LOGICAL              :: GCAP_LEAP
      INTEGER              :: AS
      INTEGER              :: DATE_T15b(2)
      INTEGER              :: NYMD_T15b, NHMS_T15b, NYMD_T15, NHMS_T15
      INTEGER              :: I,         J,         G4AHEAD
      INTEGER              :: THISYEAR,  THISMONTH, THISDAY, BACK_ONE
      REAL*8               :: JD_T15b,   JD_T15
      
      !=================================================================
      ! INIT_MEGAN begins here!
      !=================================================================

      ! Allocate arrays
      ALLOCATE( T_DAY( IIPAR, JJPAR, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'T_DAY' )
      T_DAY = 0d0

      ALLOCATE( T_15( IIPAR, JJPAR, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'T_15' )
      T_15 = 0d0

      ALLOCATE( T_15_AVG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'T_15_AVG' )
      T_15_AVG = 0d0

      ! Daily averaged temperature (mpb,2009)
      ALLOCATE( T_DAILY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'T_DAILY' )
      T_DAILY = 0d0

      ALLOCATE( AEF_ISOP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_ISOP' )
      AEF_ISOP = 0d0

      ALLOCATE( AEF_MONOT( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_MONOT' )
      AEF_MONOT = 0d0

      ALLOCATE( AEF_MBO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_MBO' )
      AEF_MBO = 0d0

      ALLOCATE( AEF_OVOC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_OVOC' )
      AEF_OVOC = 0d0

      ! New monoterpene species (mpb,2008)

      ALLOCATE( AEF_APINE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_APINE' )
      AEF_APINE = 0d0

      ALLOCATE( AEF_BPINE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_BPINE' )
      AEF_BPINE = 0d0

      ALLOCATE( AEF_LIMON( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_LIMON' )
      AEF_LIMON = 0d0

      ALLOCATE( AEF_SABIN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_SABIN' )
      AEF_SABIN = 0d0

      ALLOCATE( AEF_MYRCN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_MYRCN' )
      AEF_MYRCN = 0d0

      ALLOCATE( AEF_CAREN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_CAREN' )
      AEF_CAREN = 0d0

      ALLOCATE( AEF_OCIMN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_OCIMN' )
      AEF_OCIMN = 0d0

      ALLOCATE( AEF_SPARE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_SPARE' )
      AEF_SPARE = 0d0

      ! Allocate arrays for light (mpb,2009)

      ! -- Direct --
      ALLOCATE( PARDR_DAY( IIPAR, JJPAR, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDR_DAY' )
      T_DAY = 0d0

      ALLOCATE( PARDR_15( IIPAR, JJPAR, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDR_15' )
      T_15 = 0d0

      ALLOCATE( PARDR_15_AVG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDR_15_AVG' )
      T_15_AVG = 0d0

      ALLOCATE( PARDR_DAILY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDR_DAILY' )
      T_DAILY = 0d0

      ! -- Diffuse --
      ALLOCATE( PARDF_DAY( IIPAR, JJPAR, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDF_DAY' )
      T_DAY = 0d0

      ALLOCATE( PARDF_15( IIPAR, JJPAR, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDF_15' )
      T_15 = 0d0

      ALLOCATE( PARDF_15_AVG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDF_15_AVG' )
      T_15_AVG = 0d0

      ALLOCATE( PARDF_DAILY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDF_DAILY' )
      T_DAILY = 0d0


      ! Get annual emission factors for MEGAN inventory
#if   defined( GRID05x0666 )
      CALL GET_AEF_05x0666     ! GEOS-5 nested grids only
#else
      CALL GET_AEF             ! Global simulations
#endif 

      ! Initialize LAI arrays
      CALL INIT_LAI

      !=================================================================
      ! Read A3 fields for the 15 days before the start of the run
      ! This section has been added so that the previous 15 day temp.
      ! average can be calculated for biogenic emissions.  Do only if 
      ! MEGAN biogenic emissions must be calculated.  
      !=================================================================

      ! Get the first time for reading A-3 files
      DATE_T15b = GET_FIRST_A3_TIME()
      NYMD_T15b = DATE_T15b(1)
      NHMS_T15b = DATE_T15b(2)

      ! Astronomical Julian Date of the A3 file at start of run
      JD_T15b   = GET_JD( NYMD_T15b, NHMS_T15b )

#if   defined( GEOS_3 )

      ! For GEOS-1, GEOS-STRAT, GEOS-3, the A-3 fields are timestamped 
      ! by ending time: 00Z, 03Z, 06Z, 09Z, 12Z, 15Z, 18Z, 21Z.  
      G4AHEAD   = 0

#else

      ! For GEOS4, the A-3 fields are timestamped by the center of 
      ! the 3-hr period: 01:30Z, 04:30Z, 07:30Z, 10:30Z, 
      ! 13:30Z, 16:30Z, 19:30Z, 22:30Z
      G4AHEAD   = 13000

#endif
      
      !------------------------------------------------------
      ! GCAP: Need to test if it's leap year (phs, 9/18/07)
      !------------------------------------------------------

      ! Initialize
      THISDAY   = 0
      GCAP_LEAP = .FALSE.
      BACK_ONE  = 0

#if   defined( GCAP )

      ! Extract year, month, day from NYMD_T15b
      CALL YMD_EXTRACT( NYMD_T15b, THISYEAR, THISMONTH, THISDAY )

      ! If it's a leap year then set appropriate variables
      IF ( ITS_A_LEAPYEAR( THISYEAR, FORCE=.TRUE. )  .AND.
     &     THISMONTH == 3                            .AND.
     &     THISDAY   <  16  ) THEN 
         GCAP_LEAP = .TRUE.
         BACK_ONE  = 1
      ENDIF

#endif

      ! Remove any leftover A-3 files in temp dir (if necessary)
      IF ( LUNZIP ) THEN
         CALL UNZIP_A3_FIELDS( 'remove all' )
      ENDIF

      ! Loop over 15 days
      DO I = 15+BACK_ONE, 1, -1

         ! Skip February 29th for GCAP (phs, 9/18/07)
         IF ( GCAP_LEAP .AND. I == THISDAY ) CYCLE

         ! Julian day at start of each of the 15 days 
         JD_T15 = JD_T15b - DBLE( I ) * 1.d0

         ! Call CALDATE to compute the current YYYYMMDD and HHMMSS
         CALL CALDATE( JD_T15, NYMD_T15, NHMS_T15 )

         ! Unzip A-3 files for archving (if necessary)
         IF ( LUNZIP ) THEN
            CALL UNZIP_A3_FIELDS( 'unzip foreground', NYMD_T15 )
         ENDIF

         ! Loop over 3h periods during day
         DO J = 0, 7

            ! Open A-3 fields
            CALL OPEN_A3_FIELDS( NYMD_T15, 30000*J + G4AHEAD )

            ! Read A-3 fields from disk
            CALL GET_A3_FIELDS(  NYMD_T15, 30000*J + G4AHEAD ) 

            ! Update hourly temperatures
            CALL UPDATE_T_DAY
         ENDDO       

         ! Compute 15-day average temperatures
         CALL UPDATE_T_15_AVG

         ! Remove A-3 file from temp dir (if necessary)
         IF ( LUNZIP ) THEN
            CALL UNZIP_A3_FIELDS( 'remove date', NYMD_T15 )
         ENDIF

      ENDDO

      ! Close the A-3 file
      CLOSE( IU_A3 )

      ! Remove any leftover A-3 files in temp dir
      IF ( LUNZIP ) THEN
         CALL UNZIP_A3_FIELDS(  'remove all' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_MEGAN

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_MEGAN
!
!******************************************************************************
!  Subroutine CLEANUP_MEGAN deallocates all allocated arrays at the
!  end of a GEOS-CHEM model run. (dsa, tmf, bmy, 6/17/03, 10/24/05)
!
!  NOTES:
!******************************************************************************
!
    
      !=================================================================
      ! CLEANUP_MEGAN begins here!
      !=================================================================

      IF ( ALLOCATED( T_DAY         ) ) DEALLOCATE( T_DAY         )
      IF ( ALLOCATED( T_15          ) ) DEALLOCATE( T_15          )
      IF ( ALLOCATED( T_15_AVG      ) ) DEALLOCATE( T_15_AVG      )
      IF ( ALLOCATED( T_DAILY       ) ) DEALLOCATE( T_DAILY       )
      IF ( ALLOCATED( PARDR_DAY     ) ) DEALLOCATE( PARDR_DAY     )
      IF ( ALLOCATED( PARDR_15      ) ) DEALLOCATE( PARDR_15      )
      IF ( ALLOCATED( PARDR_15_AVG  ) ) DEALLOCATE( PARDR_15_AVG  )
      IF ( ALLOCATED( PARDR_DAILY   ) ) DEALLOCATE( PARDR_DAILY   )
      IF ( ALLOCATED( PARDF_DAY     ) ) DEALLOCATE( PARDF_DAY     )
      IF ( ALLOCATED( PARDF_15      ) ) DEALLOCATE( PARDF_15      )
      IF ( ALLOCATED( PARDF_15_AVG  ) ) DEALLOCATE( PARDF_15_AVG  )
      IF ( ALLOCATED( PARDF_DAILY   ) ) DEALLOCATE( PARDF_DAILY   )

      IF ( ALLOCATED( AEF_ISOP  ) ) DEALLOCATE( AEF_ISOP  )
      IF ( ALLOCATED( AEF_MONOT ) ) DEALLOCATE( AEF_MONOT )
      IF ( ALLOCATED( AEF_MBO   ) ) DEALLOCATE( AEF_MBO   )
      IF ( ALLOCATED( AEF_OVOC  ) ) DEALLOCATE( AEF_OVOC  )
      IF ( ALLOCATED( AEF_APINE ) ) DEALLOCATE( AEF_APINE )
      IF ( ALLOCATED( AEF_BPINE ) ) DEALLOCATE( AEF_BPINE )
      IF ( ALLOCATED( AEF_LIMON ) ) DEALLOCATE( AEF_LIMON )
      IF ( ALLOCATED( AEF_SABIN ) ) DEALLOCATE( AEF_SABIN )
      IF ( ALLOCATED( AEF_MYRCN ) ) DEALLOCATE( AEF_MYRCN )
      IF ( ALLOCATED( AEF_CAREN ) ) DEALLOCATE( AEF_CAREN )
      IF ( ALLOCATED( AEF_OCIMN ) ) DEALLOCATE( AEF_OCIMN )
      IF ( ALLOCATED( AEF_SPARE ) ) DEALLOCATE( AEF_SPARE )

      ! Return to calling program
      END SUBROUTINE CLEANUP_MEGAN

!------------------------------------------------------------------------------

      ! End of module
      END MODULE MEGAN_MOD
