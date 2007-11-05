! $Id: megan_mod.f,v 1.6 2007/11/05 16:16:22 bmy Exp $
      MODULE MEGAN_MOD
!
!******************************************************************************
!  Module MEGAN_MOD contains variables and routines specifying the 
!  algorithms that control the MEGAN inventory of biogenic emissions.
!  (dsa, tmf, bmy, 11/17/04, 9/18/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) T_DAY(R(MAXIJ,DAY_DIM): TS at each gridbox for last 8 3h periods
!  (2 ) T_15(MAXIJ,NUM_DAYS)  : Average daily TS over past NUM_DAYS days
!  (3 ) T_15_AVG (MAXIJ)      : 24h average TS over past NUM_DAYS days 
!  (4 ) DAY_DIM               : number 3h periods in a day
!  (5 ) NUM_DAYS              : Number days in averaging periods
!  (6 ) AEF_ISOP(MAXIJ)       : Annual emission factor for isoprene
!  (7 ) AEF_MONOT(MAXIJ)      : Annual emission factor for monoterpenes
!  (8 ) AEF_MBO(MAXIJ)        : Annual emission factor for methyl butenol
!  (9 ) AEF_OVOC(MAXIJ)       : Annual emission factor for other biogenic VOCs
!
!  Module Routines:
!  ============================================================================
!  (1 ) GET_EMISOP_MEGAN  : Returns ISOPRENE Emission at a gridbox 
!  (2 ) GET_EMMONOT_MEGAN : Returns MONOTERPENE Emission at a gridbox 
!  (3 ) GET_EMMBO_MEGAN   : Returns METHYL BUTENOL Emission at a gridbox 
!  (4 ) GET_EMOVOC_MEGAN  : Returns other BVOC Emission at a gridbox 
!  (5 ) GET_AEF           : Reads archived Annual emission factor (AEF) 
!  (6 ) GET_MEA           : Calculates monthly exchanged activity (MEA)
!  (7 ) GET_HEA_TL        : Calculates hourly exchange activity factor (HEA_TL)
!                            with sensitivity to BOTH TEMPERATURE and LIGHT
!                            (for ISOP and MBO)
!  (8 ) GET_HEA_T         : Calculates hourly exchange activity factor (HEA_T)
!                            with sensitivity to TEMPERATURE ONLY 
!                            (for all BVOC except ISOP and MBO)
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
!  Info about new MEGAN coefficients (tmf, 10/24/05)
!  ============================================================================
!  Updated coefficients in GET_HEA_T and GET_HEA_TL to account for the use of 
!  'surface temperature' instead of 'leaf temperature' (as in User Guide).
!  See Alex Guenther's email below:
!
!  > Hi Paul,
!  > I recall that your emissions were too high when you used the "skin"
!  > temperature from your model- and I have seen that the "skin" temperature
!  > from other models are much higher than canopy leaf temperatures.
!  > 
!  > I am not sure which coefficients you are using in the temperature
!  > algorithm but if you use
!  >    E_opt = 1.9, C_T1 = 76, C_T2 = 160, and T_opt = 316
!  > then you should use air temperature but if you are using
!  >    E_opt =1.9, C_T1a =95, C_T1b = 230, T_opt =312.5
!  > then you should be using leaf temperature to drive it.
!  >
!  > The top set of coefficients are based on results from a canopy scale
!  > model that calculates leaf temperature (i.e. it is a fit of air
!  > temperature vs canopy scale isoprene emission). So you can then say that
!  > you are accounting for the difference between air and leaf temperature.
!  > cheers,
!  > Alex
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "megan_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE :: T_DAY, T_15, T_15_AVG, DAY_DIM, NUM_DAYS

      ! PRIVATE module functions
      PRIVATE :: GET_MEA
      PRIVATE :: GET_HEA_TL
      PRIVATE :: GET_HEA_T

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      INTEGER, PARAMETER  :: DAY_DIM        = 8
      INTEGER, PARAMETER  :: NUM_DAYS       = 15
      REAL*8,  PARAMETER  :: WM2_TO_UMOLM2S = 4.766d0

      ! Arrays
      REAL*8, ALLOCATABLE :: T_DAY(:,:,:)
      REAL*8, ALLOCATABLE :: T_15(:,:,:)
      REAL*8, ALLOCATABLE :: T_15_AVG(:,:)
      REAL*8, ALLOCATABLE :: AEF_ISOP(:,:)    
      REAL*8, ALLOCATABLE :: AEF_MONOT(:,:)
      REAL*8, ALLOCATABLE :: AEF_MBO(:,:)    
      REAL*8, ALLOCATABLE :: AEF_OVOC(:,:)    

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_EMISOP_MEGAN( I,  J,      SUNCOS, 
     &                           TS, XNUMOL, Q_DIR, Q_DIFF )
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
!  (5 ) XNUMOL    (REAL*8  ) : Number of atoms C / kg C 
!  (6 ) Q_DIR     (REAL*8  ) : Flux of direct PAR above canopy [W/m2]
!  (7 ) Q_DIFF    (REAL*8  ) : Flux of diffuser PAR above canopy [W/m2]
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995
!  (2 ) Wang,    et al, 1998
!  (3 ) Guenther et al, 1999
!  (4 ) Guenther et al, 2000
!  (5 ) Guenther et al, 2004
!
!  Notes:
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/20/04)
!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE LAI_MOD,    ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I,      J
      REAL*8,  INTENT(IN) :: SUNCOS, TS, XNUMOL, Q_DIR, Q_DIFF

      ! Local variables
      INTEGER             :: IJLOOP
      REAL*8              :: MEA,     DEA,     HEA_TL 
      REAL*8              :: D_BTW_M, Q_DIR_2, Q_DIFF_2
            
      ! Function return value
      REAL*8              :: EMISOP

      !=================================================================
      ! GET_EMISOP_MEGAN begins here!
      !================================================================= 

      ! Initialize
      EMISOP   = 0.d0
      MEA      = 0.d0
      DEA      = 0.d0
      HEA_TL   = 0.d0

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

         ! If (local LAI != 0 .AND. baseline emission !=0 ) 
         IF ( ISOLAI(I,J) * AEF_ISOP(I,J) > 0d0 ) THEN

            ! Hourly exchange activitiy for temp & light
            HEA_TL = GET_HEA_TL( TS,          T_15_AVG(I,J), 
     &                           ISOLAI(I,J), SUNCOS,
     &                           Q_DIR_2,     Q_DIFF_2 )

            ! Daily exchange activity.  Alex Guenther advised us 
            ! to set DEA = 1.d0 for now (tmf, 10/24/05)
            DEA    = 1.d0

            ! Monthly exchange activity
            MEA    = GET_MEA( MISOLAI(I,J), PMISOLAI(I,J), D_BTW_M )

         ELSE

            ! If it's night or over ocean, set activity factors to zero
            HEA_TL = 0.d0
            DEA    = 0.d0
            MEA    = 0.d0

         ENDIF
    
         ! Isoprene emission is the product of all these
         EMISOP    = AEF_ISOP(I,J) * HEA_TL * DEA * MEA

         ! Convert from [kg/box] to [atoms C/box]
         EMISOP    = EMISOP * XNUMOL
         
      ENDIF

      ! return to calling program
      END FUNCTION GET_EMISOP_MEGAN

!------------------------------------------------------------------------------

      FUNCTION GET_EMMONOT_MEGAN( I, J, TS, XNUMOL ) RESULT( EMMONOT )
!
!******************************************************************************
!  Subroutine GET_EMMONOT_MEGAN computes MONOTERPENE EMISSIONS in units of 
!  [atoms C/box] using the MEGAN inventory. (dsa, tmf, 9/03, 11/20/04)  
!
!  Function GET_EMMONOT_MEGAN is called from "emissdr.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I      (INTEGER ) : GEOS-CHEM longitude index
!  (2 ) J      (INTEGER ) : GEOS-CHEM latitude index
!  (3 ) TS     (REAL*8  ) : Local surface air temperature (K)
!  (4 ) XNUMOL (REAL*8  ) : Number of atoms C / kg C 
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995
!  (2 ) Wang,    et al, 1998
!  (3 ) Guenther et al, 1999
!  (4 ) Guenther et al, 2000
!  (5 ) Guenther et al, 2004
!
!  Notes:
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/20/04)
!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE LAI_MOD,    ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I,  J
      REAL*8,  INTENT(IN) :: TS, XNUMOL

      ! Local variable
      INTEGER             :: IJLOOP
      REAL*8              :: MEA, DEA, HEA_T, D_BTW_M

      ! Function return value
      REAL*8              :: EMMONOT

      !=================================================================
      ! GET_EMMONOT_MEGAN begins here!
      !================================================================= 

      ! Initialize
      EMMONOT  = 0.d0
      MEA      = 0.d0
      DEA      = 0.d0
      HEA_T    = 0.d0

      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )

      !-----------------------------------------------------
      ! Only interested in terrestrial biosphere (pip)
      ! If (local LAI != 0 .AND. baseline emission !=0 ) 
      !-----------------------------------------------------
      IF ( ISOLAI(I,J) * AEF_MONOT(I,J) > 0d0 ) THEN

         ! Hourly exchange activity (temp only)
         HEA_T = GET_HEA_T( TS, T_15_AVG(I,J) ) 

         ! Daily exchange activity.  Alex Guenther advised us to 
         ! set DEA = 1.d0 for now (tmf, 10/24/05)
         DEA   = 1.d0

         ! Monthly exchange activity
         MEA   = GET_MEA( MISOLAI(I,J), PMISOLAI(I,J), D_BTW_M )

      ELSE

         ! Otherwise set all activities to zero
         HEA_T = 0.d0
         DEA   = 0.d0
         MEA   = 0.d0

      ENDIF    

      ! Monoterpene emissions [kg/box]
      EMMONOT  = AEF_MONOT(I,J) * HEA_T * DEA * MEA

      ! Convert from [kg/box] to [atoms C/box]
      EMMONOT  = EMMONOT * XNUMOL

      ! return to calling program
      END FUNCTION GET_EMMONOT_MEGAN

!------------------------------------------------------------------------------

      FUNCTION GET_EMMBO_MEGAN( I,  J,      SUNCOS, 
     &                          TS, XNUMOL, Q_DIR, Q_DIFF ) 
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
!  (1 ) Guenther et al, 1995
!  (2 ) Wang     et al, 1998
!  (3 ) Guenther et al, 1999
!  (4 ) Guenther et al, 2000
!  (5 ) Guenther et al, 2004
!
!  NOTES:
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/20/04)
!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE LAI_MOD,    ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I,      J
      REAL*8,  INTENT(IN) :: SUNCOS, TS, XNUMOL, Q_DIR, Q_DIFF

      ! Local variable
      REAL*8              :: MEA,     DEA,     HEA_TL
      REAL*8              :: D_BTW_M, Q_DIR_2, Q_DIFF_2

      ! Function return value
      REAL*8              :: EMMBO
            
      !=================================================================
      ! GET_EMMBO_MEGAN begins here!
      !================================================================= 

      ! Initialize
      EMMBO    = 0.d0
      MEA      = 0.d0
      DEA      = 0.d0
      HEA_TL   = 0.d0

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

         ! If local LAI > 0 and baseline emission > 0 ...
         IF ( ISOLAI(I,J) * AEF_MBO(I,J) > 0d0 ) THEN

            ! Hourly exchange activity (based on temp & light)
            HEA_TL    = GET_HEA_TL( TS,          T_15_AVG(I,J), 
     &                              ISOLAI(I,J), SUNCOS,
     &                              Q_DIR_2,     Q_DIFF_2 )

            ! Daily exchange activity.  Alex Guenther advised us 
            ! to set DEA = 1.d0 for now (tmf, 11/20/05)
            DEA    = 1.d0

            ! Monthly exchange activity
            MEA    = GET_MEA( MISOLAI(I,J), PMISOLAI(I,J), D_BTW_M )
       
         ELSE

            ! Otherwise set activities to zero
            HEA_TL = 0.d0
            DEA    = 0.d0
            MEA    = 0.d0

         ENDIF

         ! MBO emissions in [kg/box]
         EMMBO = AEF_MBO(I,J) * HEA_TL * DEA * MEA

         ! Convert from [atoms C/box] to [kg/box]
         EMMBO = EMMBO * XNUMOL
         
      ENDIF

      ! Return to calling program
      END FUNCTION GET_EMMBO_MEGAN

!------------------------------------------------------------------------------

      FUNCTION GET_EMOVOC_MEGAN( I, J, TS, XNUMOL ) RESULT( EMOVOC )
!
!******************************************************************************
!  Subroutine GET_EMOVOC_MEGAN computes other BVOC EMISSIONS in units of 
!  [atoms C/box] using the MEGAN inventory. (dsa, tmf, 9/03, 11/20/04)  
!
!  Function GET_EMOVOC_MEGAN is called from "emissdr.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I      (INTEGER ) : GEOS-CHEM longitude index
!  (2 ) J      (INTEGER ) : GEOS-CHEM latitude index
!  (2 ) TS     (REAL*8  ) : Local surface air temperature (K)
!  (3 ) XNUMOL (REAL*8  ) : Number of atoms C / kg C 
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995
!  (2 ) Wang     et al, 1998
!  (3 ) Guenther et al, 1999
!  (4 ) Guenther et al, 2000
!  (5 ) Guenther et al, 2004
!  (6 ) Guenther pers.comm, 2004
!
!  NOTES:
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/20/04)
!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE LAI_MOD,    ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M

#     include "CMN_SIZE"   ! MAXIJ

      ! Arguments
      INTEGER, INTENT(IN) :: I,  J
      REAL*8,  INTENT(IN) :: TS, XNUMOL

      ! Function return value
      REAL*8              :: EMOVOC

      ! Local variable
      REAL*8              :: MEA, DEA, HEA_T, D_BTW_M

      !=================================================================
      ! GET_EMOVOC_MEGAN begins here!
      !================================================================= 

      ! Initialize
      EMOVOC  = 0.d0
      MEA     = 0.d0
      DEA     = 0.d0
      HEA_T   = 0.d0

      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M = DBLE( DAYS_BTW_M )

      !--------------------------------------------------
      ! Only interested in terrestrial biosphere (pip)
      ! If (local LAI != 0 .AND. baseline emission !=0 ) 
      !--------------------------------------------------
      IF ( ISOLAI(I,J) * AEF_OVOC(I,J) > 0d0 ) THEN

         ! Hourly exchange activity (temp only)
         HEA_T = GET_HEA_T( TS, T_15_AVG(I,J ) ) 

         ! Daily exchange activity.  Alex Guenther advised us 
         ! to set DEA = 1.d0 for now (tmf, 10/24/05)
         DEA   = 1.d0

         ! Monthly exchange activity
         MEA   = GET_MEA( MISOLAI(I,J), PMISOLAI(I,J), D_BTW_M )
       
      ELSE
         
         ! Otherwise set activities to zero
         HEA_T = 0.d0
         DEA   = 0.d0
         MEA   = 0.d0

      ENDIF

      ! OVOC emissions [kg/box]
      EMOVOC = AEF_OVOC(I,J) * HEA_T * DEA * MEA

      ! Convert from [kg/box] to [atoms C/box]
      EMOVOC = EMOVOC * XNUMOL

      ! return to calling program
      END FUNCTION GET_EMOVOC_MEGAN

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
#     include "CMN_SIZE" ! Size parameters

      ! Local Variables
      INTEGER           :: I, J, D

      ! External functions
      REAL*8, EXTERNAL  :: XLTMMP

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
            T_DAY(I,J,D) = T_DAY(I,J,D-1)
         ENDDO
            
         ! Store 
         T_DAY(I,J,1) = XLTMMP(I,J)
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

      !=================================================================
      ! UPDATE_T_15_AVG begins here!
      !=================================================================

      ! Convert to REAL*8
      D_DIM  = DBLE( DAY_DIM  )
      D_DAYS = DBLE( NUM_DAYS )

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, D, TMP_T )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Average T_DAY over the 8 TS values during the day.
         TMP_T = SUM( T_DAY(I,J,:) ) / D_DIM
         
         ! Push the daily average TS values through T_15,
         ! throwing out the oldest 
         DO D = NUM_DAYS, 2, -1
            T_15(I,J,D) = T_15(I,J,D-1)
         ENDDO

         ! Put the newest daily average TS value in the first spot
         T_15(I,J,1) = TMP_T

         ! Get T_15_AVG by averaging T_15 over the 15 day period.
         T_15_AVG(I,J) = SUM( T_15(I,J,:) ) / D_DAYS 

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE UPDATE_T_15_AVG

!------------------------------------------------------------------------------
! NOTE: This subroutine is not called by anything, so commment out for now.
! (bmy, 10/24/05)
!
!      SUBROUTINE MEGAN_FACTORS( I,      J,     SUNCOS, TS, 
!     &                          XNUMOL, Q_DIR, Q_DIFF, MEA, 
!     &                          DEA,    HEA_T, HEA_TL )
!!
!!******************************************************************************
!!  Subroutine MEGAN_FACTORS computes the emission activity factors in 
!!   the MEGAN inventory. (dsa, tmf, bmy, 9/03, 10/24/05)  
!!
!!  Function GET_EMISOP_MEGAN is called from "emissdr.f"
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) I      (INTEGER ) : GEOS-CHEM longitude index
!!  (2 ) J      (INTEGER ) : GEOS-CHEM latitude index
!!  (3 ) SUNCOS (REAL*8  ) : 1-D array of cos( solar zenith angle )
!!  (4 ) TS     (REAL*8  ) : Local surface air temperature (K)
!!  (5 ) XNUMOL (REAL*8  ) : Number of atoms C / kg C 
!!  (6 ) Q_DIR  (REAL*8  ) : flux of direct PAR above canopy (W m-2)
!!  (7 ) Q_DIFF (REAL*8  ) : flux of diffuser PAR above canopy (W m-2)
!!
!!  Arguments as Output:
!!  ============================================================================
!!  (9 ) MEA    (REAL*8  ) : MEGAN monthly emission factor
!!  (10) DEA    (REAL*8  ) : MEGAN daily emission factor
!!  (11) HEA_T  (REAL*8  ) : MEGAN hourly emission factor as F(T)
!!  (12) HEA_TL (REAL*8  ) : MEGAN hourly emission factor as F(T,light)
!!
!!  References (see above for full citations):
!!  ============================================================================
!!  (1 ) Guenther et al, 1995
!!  (2 ) Wang     et al, 1998
!!  (3 ) Guenther et al, 1999
!!  (4 ) Guenther et al, 2000
!!  (5 ) Guenther et al, 2004
!!
!!  NOTES:
!!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!!        algorithm and modified for the standard code by May Fu (11/20/04)
!!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!!        met field. TS is the surface air temperature, which should be 
!!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!!
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE LAI_MOD,     ONLY : ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M
!
!#     include "CMN_SIZE"    ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: I,J
!      REAL*8,  INTENT(IN)  :: SUNCOS(MAXIJ)
!      REAL*8,  INTENT(IN)  :: TS,  XNUMOL, Q_DIR, Q_DIFF
!      REAL*8,  INTENT(OUT) :: MEA, DEA,    HEA_T, HEA_TL
!
!      ! Local variables
!      INTEGER              :: IJLOOP
!      REAL*8               :: D_BTW_M, Q_DIR_2, Q_DIFF_2
!            
!      !=================================================================
!      ! MEGAN_FACTORS begins here!
!      !================================================================= 
!
!      ! 1-D array index
!      IJLOOP   = ( (J-1) * IIPAR ) + I
!
!      ! Initialize
!      MEA      = 0d0
!      DEA      = 0d0
!      HEA_T    = 0d0
!      HEA_TL   = 0d0
!
!      ! Number of days between MISOLAI and PMISOLAI
!      D_BTW_M  = DBLE( DAYS_BTW_M )
!
!      ! cCnvert Q_DIR and Q_DIFF from [W/m2] to [umol/m2/s]
!      Q_DIR_2  = Q_DIR  * WM2_TO_UMOLM2S
!      Q_DIFF_2 = Q_DIFF * WM2_TO_UMOLM2S 
!
!      !---------------------------------------------------
!      ! If local LAI > 0 and baseline emission > 0 ...
!      !---------------------------------------------------
!      IF ( ISOLAI(I,J) * AEF_ISOP(I,J) > 0d0 ) THEN
!
!         ! Hourly exchange factor (temp only)
!         HEA_T  = GET_HEA_T( TS, T_15_AVG(I,J) ) 
!
!         ! For light sensitive species: do only during day
!         IF ( SUNCOS(IJLOOP) > 0d0 ) THEN
!
!            ! Hourly exchange factor (temp & light)
!            HEA_TL = GET_HEA_TL( TS,          T_15_AVG(I,J), 
!     &                           ISOLAI(I,J), SUNCOS(IJLOOP),
!     &                           Q_DIR_2,     Q_DIFF_2 )
!         ELSE 
!            
!            ! Otherwise set to zero
!            HEA_TL = 0.d0
!
!         ENDIF
!
!         ! Daily exchange activity.  Alex Guenther suggests 
!         ! that we set DEA = 1.d0 for now
!         DEA       = 1.d0
! 
!         ! Monthly exchange activity
!         MEA       = GET_MEA( MISOLAI(I,J), PMISOLAI(I,J), D_BTW_M )     
!
!      ELSE
!
!         ! Otherwise set activities to zero
!         HEA_T     = 0.d0
!         HEA_TL    = 0.d0
!         DEA       = 0.d0
!         MEA       = 0.d0
!
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE MEGAN_FACTORS
!
!------------------------------------------------------------------------------

      FUNCTION GET_HEA_TL( T,       PT_15,   LAI, 
     &                     SUNCOS1, Q_DIR_2, Q_DIFF_2 ) RESULT( HEA_TL )
!
!******************************************************************************
!  Computes the hourly exchange activity factor (HEA_TL) with sensitivity to 
!  both TEMPERATHRE and LIGHT for ISOP emission.  (tmf, 11/18/04, 10/24/05)
!
!  HEA_TL = C_T * C_PPFD
!     
!  C_T:    Effect of temperature on leaf isoprene emission, including effect 
!          of average temperature over previous 15 days, based on Eq 5abc 
!          from Guenther et al. (1999)
!
!  C_PPFD: Effect of increasing PPFD up to a saturation point, where emission 
!          level off, based on Eq 4abc from Guenther et al. (1999)
!          In addition, a 5 layered canopy model based on Eqs 12-16 
!          from Guenther et al. (1995) is included to correct for light 
!          attenuation in the canopy.
!  
!  Function GET_HEA_TL is called from GET_EMISOP_MEGAN of "megan_mod.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) T    (REAL*8 )    : Current leaf temperature, the surface air temp
!                           field (TS) is assumed equivalent to the leaf
!                           temperature over forests.
!  (2 ) PT_15 (REAL*8 )   : Average leaf tempearture over the past 15 days
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
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/2004)
!  (2 ) Algorithm is based on Guenther (1999) instead of the latest User's 
!        Guide since the User's Guide does not contain sensitivity to 
!        temperature in the past 15 days. (tmf, 11/19/04)
!  (3 ) Algorithm updated to Users' Guide as of Nov 3, 2004.  Sensitivity on 
!        temperature history is included in DEA from GET_DEA. (tmf, 11/19/04)
!  (4 ) Significant modifications have been made to the canopy model in 
!        accordence with the MEGAN isoprene emissions model from Guenther et 
!        al. (unpublished as of 6/9/03) (dsa, 6/5/03)
!  (5 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (6 ) Switch back to Guenther et al. (1999) algorithm and account for 
!        tempearture history in HEA.  This is accompanied by setting DEA = 1.
!        (Alex Guenther, personal communication, 11/25/04)
!******************************************************************************
!
#     include "CMN_GCTM"   ! Physical constants

      ! Arguments
      REAL*8,  INTENT(IN) :: T, PT_15
      REAL*8,  INTENT(IN) :: LAI, SUNCOS1, Q_DIR_2, Q_DIFF_2

      ! Function return value
      REAL*8              :: HEA_TL

      ! Local Variables

      !-----------------------------------------------------------------
      ! Based on Eq 5a, 5b, and 5c in Guenther et al. (1999)
      !-----------------------------------------------------------------

      ! C_T: Effect of temperature on leaf BVOC emission, including 
      ! effect of average temperature over previous 15 days, based on 
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      REAL*8              :: C_T

      ! E_OPT: maximum normalized emission capacity
      REAL*8              :: E_OPT

      ! T_OPT: temperature at which E_OPT occurs
      REAL*8              :: T_OPT

      ! X: variable related to temp.
      REAL*8              :: X

      ! CT1, CT2: energy of activation and deactivation, respectively 
      REAL*8              :: CT1, CT2

      ! R - ideal gas constant (J mol-1 K-1)
      REAL*8, PARAMETER   :: R = 8.314d-3


      !-----------------------------------------------------------------
      ! Canopy model variables (Eqs 12-16 from Guenther et al, 1995)
      !-----------------------------------------------------------------

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
      REAL*8, PARAMETER   :: NORMAL_FACTOR = 1.8967d0

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

      !=================================================================
      ! GET_HEA_TL begins here!
      !================================================================= 

      !-----------------------------------------------------------------
      ! Compute C_T
      !
      ! NOTES:
      ! (1) Algorithm is from Eq 5a, 5b, and 5c in Guenther et al, 1999
      ! (2) Alex Guenther suggests that we set DEA=1 for now
      ! (3) E_OPT and T_OPT depend on PT_15 according to Eq 5b and 5c
      !-----------------------------------------------------------------

      ! E_OPT: maximum normalized emission capacity
      E_OPT = 1.9d0 * EXP( 1.25d-1 * ( PT_15 - 3.01d2 ) )

      ! T_OPT: temperature at which E_OPT occurs
      T_OPT = 3.16d2 + ( 5.0d-1 * ( PT_15 - 3.01d2 ) )

      ! Energy of activation 
      CT1   = 76d0

      ! Energy of deactivation 
      CT2   = 160d0

      ! Variable related to temperature 
      X     = ( 1.d0/T_OPT - 1.d0/T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including 
      ! effect of average temperature over previous 15 days, based on 
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) / 
     &        ( CT2 - CT1 * ( 1.d0 - EXP( CT2 * X ) ) )

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

      ! Finally, HEA_TL.  Prevent negative values.
      HEA_TL = MAX( C_T * C_PPFD / NORMAL_FACTOR, 0d0 )

      ! return to calling program
      END FUNCTION GET_HEA_TL

!------------------------------------------------------------------------------

      FUNCTION GET_HEA_T( T, PT_15 ) RESULT( HEA_T )
!
!******************************************************************************
!  Computes the hourly exchange activity factor (HEA_T) with sensitivity 
!  to only TEMPERATHRE for all BVOC emission except ISOP. 
!  (tmf, bmy, 11/25/04, 10/24/05)
!
!  Function GET_HEA_T is called from GET_EMISOP_MEGAN of "megan_mod.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) T     (REAL*8 ) : Current leaf temperature, the surface air 
!                          temperature field (TS) is assumed equivalent to 
!                          the leaf temperature over forests.
!  (2 ) PT_15 (REAL*8 ) : Average leaf temperature over the past 15 days
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995
!  (2 ) Wang     et al, 1998
!  (3 ) Guenther et al, 1999
!  (5 ) Guenther et al, 2004
!
!  Notes:
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/2004)
!  (2 ) Algorithm is based on Guenther (1999) instead of the latest User's 
!        Guide, since the User's Guide does not contain sensitivity to !
!        temperature in the past 15 days. (tmf, 11/19/04)
!  (3 ) Algorithm updated to Users' Guide as of Nov 3, 2004.  Sensitivity on 
!        temperature history is included in DEA from GET_DEA. (tmf, 11/19/04)
!  (4 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (5 ) Switch back to Guenther et al. (1999) algorithm and account for !
!        tempearture history in HEA.  This is accompanied by setting DEA = 1. 
!        (Alex Guenther, personal communication, 11/25/04)
!  (6 ) For consistency, C_T is calculated according to Guenther et al. (1999) 
!        algorithm.  Alex Guenther suggested that we use HEA TYPE 2 algorithm 
!        in Eq 12 from User's guide.  However, this makes incorporation of 
!        temperature history more difficult.  tmf checked the difference 
!        between Eq 9 and 12 from User's guide and they are essentially the 
!        same.  Therefore we see no need to differentiate the C_T in HEA_T 
!        and HEA_TL.
!******************************************************************************
!
      ! Arguments
      REAL*8,  INTENT(IN) :: T, PT_15

      ! Local Variables
      REAL*8              :: C_T,   CT1,   CT2
      REAL*8              :: E_OPT, T_OPT, X

      ! Ideal gas constant [J/mol/K]
      REAL*8, PARAMETER   :: R   = 8.314d-3

      ! Function return value
      REAL*8              :: HEA_T

      !=================================================================
      ! GET_HEA_T begins here!
      !
      ! NOTES:
      ! (1) Algorithm is from Eq 5a, 5b, and 5c in Guenther et al, 1999
      ! (2) Alex Guenther suggests that we set DEA=1 for now
      ! (3) E_OPT and T_OPT depend on PT_15 according to Eq 5b and 5c
      !================================================================= 

      ! E_OPT: maximum normalized emission capacity
      E_OPT = 1.9d0 * EXP( 1.25d-1 * ( PT_15 - 3.01d2 ) )

      ! T_OPT: temperature at which E_OPT occurs
      T_OPT = 3.16d2 + ( 5.0d-1 * ( PT_15 - 3.01d2 ) )

      ! Energy of activation 
      CT1   = 76d0

      ! Energy of deactivation 
      CT2   = 160d0

      ! Variable related to temperature 
      X     = ( 1.d0/T_OPT - 1.d0/T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including 
      ! effect of average temperature over previous 15 days, based on 
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) / 
     &        ( CT2 - CT1 * ( 1.d0 - EXP( CT2 * X ) ) )

      ! Hourly emission activity = C_T
      ! Prevent negative values
      HEA_T = MAX( C_T, 0d0 )

      ! Return to calling program
      END FUNCTION GET_HEA_T

!------------------------------------------------------------------------------

      FUNCTION GET_MEA( CMLAI, PMLAI, T ) RESULT( MEA )
!
!******************************************************************************
!  Computes the monthly exchange activity factor (MEA). (tmf, 10/24/05)
!     MEA = M_LAI * M_AGE * M_H
!  
!  Function GET_MEA is called from GET_EMISOP_MEGAN of "megan_mod.f"
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) CMLAI (REAL*8 ) : Current month leaf area index at gridbox 
!  (2 ) PMLAI (REAL*8 ) : Last month leaf area index at gridbox 
!  (3 ) T     (REAL*8 ) : Number of days between current and previous LAI.
!
!  References (see above for full citations):
!  ============================================================================
!  (3 ) Guenther et al, 1999
!  (5 ) Guenther et al, 2004
!
!  NOTES:
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) Update to publically released (as of 11/2004) MEGAN algorithm and 
!        modified for the standard code by May Fu (11/2004).
!  (3 ) Algorithm is based on the latest User's Guide (tmf, 11/19/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: T
      REAL*8, INTENT(IN) :: CMLAI, PMLAI

      ! Function return value
      REAL*8             :: MEA        
          
      ! Local Variables

      ! M_LAI: leaf area factor
      REAL*8             :: M_LAI 

      ! M_AGE: leaf age factor
      REAL*8             :: M_AGE

      ! M_H: sensible heat flux factor      
      REAL*8             :: M_H

      ! FNEW, ANEW: new foliage that emits negligible amounts of isoprene
      REAL*8             :: FNEW
      REAL*8,  PARAMETER :: ANEW = 1.d-2

      ! FGRO, AGRO: growing foliage that emits isoprene at low rates
      REAL*8             :: FGRO
      REAL*8,  PARAMETER :: AGRO = 5.d-1

      ! FMAT, AMAT: mature foliage that emits isoprene at peak rates
      REAL*8             :: FMAT
      REAL*8,  PARAMETER :: AMAT = 1.d0

      ! FSEN, ASEN: senescing foliage that emits isoprene at reduced rates
      REAL*8             :: FSEN
      REAL*8,  PARAMETER :: ASEN = 3.3d-1

      ! TI: number of days after budbreak required to induce iso. em.
      REAL*8,  PARAMETER :: TI   = 12.d0

      ! TM: number of days after budbreak required to reach peak iso. em. rates
      REAL*8,  PARAMETER :: TM   = 28.d0

      ! Variable for storing T or TM
      REAL*8             :: TG

      !=================================================================
      ! GET_MEA begins here!
      !================================================================= 

      !-----------------------
      ! Compute M_LAI
      !-----------------------
      M_LAI = 0.49d0 * CMLAI / SQRT( 1.d0 + 0.2d0 * CMLAI*CMLAI )
      
      !-----------------------
      ! Compute M_AGE
      !-----------------------
      IF ( T > TM ) THEN 
         TG = TM 
      ELSE 
         TG = T
      ENDIF

      IF ( CMLAI == PMLAI ) THEN

         FMAT = 1.d0  
         FNEW = 0.d0
         FGRO = 0.d0
         FSEN = 0.d0

      ELSE IF ( CMLAI > PMLAI ) THEN

         FSEN = 0.d0

         IF ( T > TI ) THEN
            FNEW = ( TI / T ) * ( 1.d0 -  PMLAI / CMLAI )
            FGRO = ( ( TG - TI ) / T ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE
            FNEW = 1.d0 - ( PMLAI / CMLAI )
            FGRO = 0.d0
         ENDIF

         IF ( T > TM ) THEN
            FMAT = ( PMLAI / CMLAI ) +
     &             ( ( T - TM ) / T ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE 
            FMAT = ( PMLAI / CMLAI )
         ENDIF

      ELSE

         FSEN = ( PMLAI - CMLAI ) / PMLAI
         FMAT = 1.d0 - FSEN
         FGRO = 0.d0
         FNEW = 0.d0

      ENDIF

      ! Age factor
      M_AGE = FNEW*ANEW + FGRO*AGRO + FMAT*AMAT + FSEN*ASEN 

      !--------------------------------------------------------------
      ! Compute M_H = 1.d0 + 0.03d0 * (H_MONTH - H_SEASON)
      ! 
      ! where H_MONTH is the daytime sensible heat flux of the past 
      ! month and H_SEASON is the daytime sensible heat flux of the 
      ! entire growing season
      !--------------------------------------------------------------

      ! Disable M_H for now (tmf, 10/24/05)
      M_H = 1.d0     

      !-------------------------------------
      ! Compute MEA = M_LAI * M_AGE * M_H
      !-------------------------------------

      ! Prevent negative values
      MEA = MAX( M_LAI * M_AGE * M_H, 0d0 )

      ! return to calling program
      END FUNCTION GET_MEA

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

      !=================================================================
      ! GET_AEF begins here!
      !=================================================================

      ! Emission timestep [min]
      DTSRCE = GET_TS_EMIS()

      !---------------------------------------------
      ! Read in ISOPRENE Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &           'MEGAN_200510/MEGAN_AEF_ISOP.geos.1x1'

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
            AEF_ISOP(I,J) = AEF_ISOP(I,J) * FACTOR

         ENDDO
      ENDDO

      !---------------------------------------------
      ! Read in MONOTERPENE Annual Emission Factors
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
            AEF_MONOT(I,J) = AEF_MONOT(I,J) * FACTOR

         ENDDO
      ENDDO

      !---------------------------------------------
      ! Read in MBO Annual Emission Factors
      !---------------------------------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &           'MEGAN_200510/MEGAN_AEF_MBO.geos.1x1'

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
            AEF_MBO(I,J) = AEF_MBO(I,J) * FACTOR
            
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

      ! Return to calling program
      END SUBROUTINE GET_AEF

!------------------------------------------------------------------------------

      SUBROUTINE INIT_MEGAN
!
!******************************************************************************
!  Subroutine INIT_MEGAN allocates and initializes the T_DAY, T_15,  
!  T_15_AVG, and AEF_* arrays. (dsa, tmf, bmy, 10/24/05, 9/18/07)
!
!  NOTES:
!  (1 ) Change the logic in the #if block for G4AHEAD. (bmy, 12/6/05)
!  (2 ) Bug fix: skip Feb 29th if GCAP (phs, 9/18/07)
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

      ! Get annual emission factors for MEGAN inventory
      CALL GET_AEF

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
      !---------------------------------
      ! Prior to 9/18/07:
      ! Modified for GCAP (phs, 9/18/07)
      !DO I = 15, 1, -1
      !---------------------------------
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
      IF ( ALLOCATED( T_DAY     ) ) DEALLOCATE( T_DAY     )
      IF ( ALLOCATED( T_15      ) ) DEALLOCATE( T_15      )
      IF ( ALLOCATED( T_15_AVG  ) ) DEALLOCATE( T_15_AVG  )
      IF ( ALLOCATED( AEF_ISOP  ) ) DEALLOCATE( AEF_ISOP  )
      IF ( ALLOCATED( AEF_MONOT ) ) DEALLOCATE( AEF_MONOT )
      IF ( ALLOCATED( AEF_MBO   ) ) DEALLOCATE( AEF_MBO   )
      IF ( ALLOCATED( AEF_OVOC  ) ) DEALLOCATE( AEF_OVOC  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_MEGAN

!------------------------------------------------------------------------------

      ! End of module
      END MODULE MEGAN_MOD
