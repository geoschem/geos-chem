! $Id: megan_mod.f,v 1.4 2010/03/09 21:44:17 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: megan_mod
!
! !DESCRIPTION: Module MEGAN\_MOD contains variables and routines specifying 
!  the algorithms that control the MEGAN inventory of biogenic emissions.
!\\
!\\
!  References:
!  \begin{itemize}
!  \item Guenther, A., et al., \emph{A global model of natural volatile 
!        organic compound emissions}, \underline{J.Geophys. Res.}, 
!        \textbf{100}, 8873-8892, 1995.
!  \item Wang, Y., D. J. Jacob, and J. A. Logan, \emph{Global simulation of 
!        tropospheric O3-Nox-hydrocarbon chemistry: 1. Model formulation}, 
!        \underline{J. Geophys. Res.}, \textbf{103}, D9, 10713-10726, 1998.
!  \item Guenther, A., B. Baugh, G. Brasseur, J. Greenberg, P. Harley, L. 
!        Klinger, D. Serca, and L. Vierling, \emph{Isoprene emission estimates 
!        and uncertanties for the Central African EXPRESSO study domain}, 
!        \underline{J. Geophys. Res.}, \textbf{104}, 30,625-30,639, 1999.
!  \item Guenther, A. C., T. Pierce, B. Lamb, P. Harley, and R. Fall, 
!        \emph{Natural emissions of non-methane volatile organic compounds, 
!        carbon monoxide, and oxides of nitrogen from North America}, 
!        \underline{Atmos. Environ.}, \textbf{34}, 2205-2230, 2000.
!  \item Guenther, A., and C. Wiedinmyer, \emph{User's guide to Model of 
!        Emissions of Gases and Aerosols from Nature}. http://cdp.ucar.edu. 
!        (Nov. 3, 2004) 
!  \item Guenther, A., \emph{AEF for methyl butenol}, personal commucation. 
!        (Nov, 2004)
!  \end{itemize}
!
! !INTERFACE:
!
      MODULE MEGAN_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD                                 ! Size parameters
      USE CMN_GCTM_MOD                                 ! Physical constants
      USE ERROR_MOD                                    ! Error trapping

      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !DEFINED PARAMETERS:
! 
      ! Scalar
#if   defined( MERRA ) || defined( GEOS_57 )
      INTEGER, PARAMETER  :: DAY_DIM        = 24       ! # of 1-hr periods/day
#else
      INTEGER, PARAMETER  :: DAY_DIM        = 8        ! # of 3-hr periods/day
#endif
      INTEGER, PARAMETER  :: NUM_DAYS       = 10       ! # of days to avg 
      REAL*8,  PARAMETER  :: WM2_TO_UMOLM2S = 4.766d0  ! W/m2 -> umol/m2/s

      REAL*8, PARAMETER   :: D2RAD =  PI_180           ! Degrees to radians
      REAL*8, PARAMETER   :: RAD2D =  1d0 / PI_180     ! Radians to degrees

!
! !PRIVATE TYPES:
!
      ! Past light & temperature conditions (mpb,2009)
      ! (1) Temperature at 2m (TS):
      REAL*8, ALLOCATABLE :: T_DAILY(:,:)       ! Daily averaged sfc temp
      REAL*8, ALLOCATABLE :: T_DAY(:,:,:)       ! Holds 1 day of sfc temp data
      REAL*8, ALLOCATABLE :: T_15(:,:,:)        ! Holds 15 days of daily avg T
      REAL*8, ALLOCATABLE :: T_15_AVG(:,:)      ! Sfc temp avg'd over NUM_DAYS

      ! (2) PAR Direct: 
      REAL*8, ALLOCATABLE :: PARDR_DAILY(:,:)   ! Average daily PARDR
      REAL*8, ALLOCATABLE :: PARDR_DAY(:,:,:)   ! Holds 1 day of PARDR data
      REAL*8, ALLOCATABLE :: PARDR_15(:,:,:)    ! 10 days of daily avg'd PARDR
      REAL*8, ALLOCATABLE :: PARDR_15_AVG(:,:)  ! PARDR averaged over NUM_DAYS

      ! (3) PAR Diffuse: 
      REAL*8, ALLOCATABLE :: PARDF_DAILY(:,:)   ! Average daily PARDR
      REAL*8, ALLOCATABLE :: PARDF_DAY(:,:,:)   ! Holds 1-day of PARDR data
      REAL*8, ALLOCATABLE :: PARDF_15(:,:,:)    ! 10 days of daily avg'd PARDR 
      REAL*8, ALLOCATABLE :: PARDF_15_AVG(:,:)  ! PARDF averaged over NUM_DAYS

      ! Annual emission factor arrays (mpb,2009)
      REAL*8, ALLOCATABLE :: AEF_ISOP(:,:)      ! Isoprene, atomsC/box 
      REAL*8, ALLOCATABLE :: AEF_MBOX(:,:)      ! Methyl butenol, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_APIN(:,:)      ! Alpha-pinene, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_BPIN(:,:)      ! Beta-pinene, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_CARE(:,:)      ! 3-Carene, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_LIMO(:,:)      ! Limonene, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_MYRC(:,:)      ! Myrcene, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_OCIM(:,:)      ! Ocimen, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_SABI(:,:)      ! Sabine, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_OMON(:,:)      ! Other monoterpenes, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_MOHX(:,:)      ! Methanol, molec/box
      REAL*8, ALLOCATABLE :: AEF_ACET(:,:)      ! Acetone, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_ETOH(:,:)      ! Ethanol, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_CH2O(:,:)      ! Formaldehyde, molec/box  
      REAL*8, ALLOCATABLE :: AEF_ALD2(:,:)      ! Acetaldehyde, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_FAXX(:,:)      ! Formic acid, molec/box
      REAL*8, ALLOCATABLE :: AEF_AAXX(:,:)      ! Acetic acid, molec/box
      REAL*8, ALLOCATABLE :: AEF_C2H4(:,:)      ! Ethene, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_TOLU(:,:)      ! Toluene, atomsC/box
      REAL*8, ALLOCATABLE :: AEF_HCNX(:,:)      ! HCN, molec/ box
      REAL*8, ALLOCATABLE :: AEF_PRPE(:,:)      ! >= C3 alkenes, atomsC/box
      REAL*8, ALLOCATABLE :: NORM_FAC(:)        ! Normalization factor

      ! Path to MEGAN emission factors
      CHARACTER(LEN=20)   :: MEGAN_SUBDIR = 'MEGAN_201211/'
!
! !PUBLIC MEMBER FUNCTIONS: 
!
      PUBLIC  :: GET_MEGAN_EMISSIONS ! dbm, new MEGAN driver routine
                                     ! for all compounds (6/21/2012)
      PUBLIC  :: INIT_MEGAN        
      PUBLIC  :: CLEANUP_MEGAN
      PUBLIC  :: UPDATE_T_DAY      
      PUBLIC  :: UPDATE_T_15_AVG   
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GET_MEGAN_PARAMS
      PRIVATE :: GET_MEGAN_AEF
      PRIVATE :: GET_GAMMA_PAR_PCEEA
      PRIVATE :: GET_GAMMA_T_LI
      PRIVATE :: GET_GAMMA_T_LD
      PRIVATE :: GET_GAMMA_LAI
      PRIVATE :: GET_GAMMA_AGE
      PRIVATE :: GET_GAMMA_SM
      PRIVATE :: CALC_NORM_FAC
      PRIVATE :: SOLAR_ANGLE
      PRIVATE :: CALC_AEF

! 
! !REVISION HISTORY: 
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
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  09 Mar 2010 - R. Yantosca - Minor bug fix in GET_EMMONOT_MEGAN
!  17 Mar 2010 - H. Pye      - AEF_SPARE must be a scalar local variable
!                              in GET_EMMONOT_MEGAN for parallelization.
!  20 Aug 2010 - R. Yantosca - Move CMN_SIZE to top of module
!  20 Aug 2010 - R. Yantosca - Now set DAY_DIM = 24 for MERRA, since the
!                              surface temperature is now an hourly field.
!  22 Nov 2011 - R. Yantosca - Do not use erroneous AEF's for nested grids 
!  06 Dec 2011 - E. Fischer  - Added Acetone emissions 
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!  11 Apr 2012 - R. Yantosca - Replace lai_mod.F with modis_lai_mod.F90
!  01 Sep 2010 - R. Yantosca - Bug fix in INIT_MEGAN: now only read in 
!                              NUM_DAYS (instead of 15) days of sfc temp data
!  19 Jan 2013 - D. Millet - Streamlined computations into a single driver routine
!                            and updated emissions according to MEGAN 2.1 as
!                            described in:
!                            Guenther et al., The Model of Emissions of Gases
!                            and Aerosols from Nature version 2.1 (MEGAN2.1):
!                            an extended and updated framework for modeling 
!                            biogenic emissions, GMD, 5, 1471-1492, 2012.
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_emisop_megan
!
! !DESCRIPTION: Subroutine GET\_EMISOP\_MEGAN computes isoprene emissions in 
!  units of [atoms C/box] using the MEGAN inventory. (dbm, 12/2012)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MEGAN_EMISSIONS( I,  J,     SUNCOS, 
     &                              TS, Q_DIR, Q_DIFF, CMPD )
     &     RESULT( MEGAN_EMIS )
!
! !USES:
!
      USE MODIS_LAI_MOD, ONLY : ISOLAI     => GC_LAI       ! Daily LAI
      USE MODIS_LAI_MOD, ONLY : PMISOLAI   => GC_LAI_PM    ! Prev month's LAI
      USE MODIS_LAI_MOD, ONLY : MISOLAI    => GC_LAI_CM    ! Curr month's LAI
      USE MODIS_LAI_MOD, ONLY : DAYS_BTW_M => DAYS_BTW_MON ! LAI month interval
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      INTEGER, INTENT(IN)            :: I, J    ! GEOS-Chem lon & lat indices
      REAL*8,  INTENT(IN)            :: SUNCOS  ! Solar zenith angle [unitless]
      REAL*8,  INTENT(IN)            :: TS      ! Surface temperature [K]
      REAL*8,  INTENT(IN)            :: Q_DIR   ! Flux of direct PAR above canopy [W/m2]
      REAL*8,  INTENT(IN)            :: Q_DIFF  ! Flux of diffuse PAR above canopy [W/m2]
      CHARACTER(LEN=*),  INTENT(IN)  :: CMPD    ! Compound name (dbm, 6/21/2012)
!
! !RETURN VALUE:
!
      REAL*8                          :: MEGAN_EMIS ! VOC emission in atomsC/box 
                                                    ! or molec/box, depending on
                                                    ! units the compound is carried in
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2000, 2004, 2006
!  (2 ) Wang,    et al, 1998
!  (3 ) Guenther et al, 2007, MEGAN v2.1 User mannual 
!  (4 ) Guenther et al., 2012 GMD MEGANv2.1 description and associated code at
!                                 http://acd.ucar.edu/~guenther/MEGAN/
!
! !REVISION HISTORY: 
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/20/04)
!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (3 ) Restructing of function & implementation of activity factors (mpb,2009)
!  (4 ) Dbm, new driver routine for all MEGAN compounds (6/21/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8              :: GAMMA_LAI
      REAL*8              :: GAMMA_AGE
      REAL*8              :: GAMMA_PAR
      REAL*8              :: GAMMA_T_LD
      REAL*8              :: GAMMA_T_LI
      REAL*8              :: GAMMA_SM
      REAL*8              :: AEF
      REAL*8              :: D_BTW_M, Q_DIR_2, Q_DIFF_2
      REAL*8              :: BETA, LDF, CT1, CEO
      REAL*8              :: ANEW, AGRO, AMAT, AOLD
      LOGICAL             :: BIDIR      
      
      !=================================================================
      ! GET_MEGAN_EMISSIONS begins here!
      !================================================================= 
      
      ! Initialize parameters, gamma values, and return value
      MEGAN_EMIS = 0d0
      GAMMA_LAI  = 0d0
      GAMMA_AGE  = 0d0
      GAMMA_T_LD = 0d0
      GAMMA_T_LI = 0d0
      GAMMA_PAR  = 0d0
      GAMMA_SM   = 0d0
      BETA       = 0d0
      AEF        = 0d0
      LDF        = 0d0
      CT1        = 0d0
      CEO        = 0d0
      ANEW       = 0d0
      AGRO       = 0d0
      AMAT       = 0d0
      AOLD       = 0d0
      BIDIR      = .FALSE.     

      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )
      
      ! Convert Q_DIR and Q_DIFF from (W/m2) to (micromol/m2/s)
      Q_DIR_2  = Q_DIR  * WM2_TO_UMOLM2S
      Q_DIFF_2 = Q_DIFF * WM2_TO_UMOLM2S

      ! --------------------------------------------
      ! Get MEGAN parameters for this compound
      ! --------------------------------------------
      CALL GET_MEGAN_PARAMS ( CMPD, BETA, LDF,  CT1,  CEO, 
     &                        ANEW, AGRO, AMAT, AOLD, BIDIR )

      ! --------------------------------------------
      ! Get base emission factor for this compound and grid square
      ! Units: atomsC/box or molec/box
      ! --------------------------------------------
      AEF = GET_MEGAN_AEF ( I, J, CMPD )
     
      !-----------------------------------------------------
      ! Only interested in terrestrial biosphere
      ! If (local LAI != 0 .AND. baseline emission !=0 ) 
      !-----------------------------------------------------
      IF ( ISOLAI(I,J) * AEF > 0d0 ) THEN

         ! --------------------------------------------------
         ! GAMMA_par (light activity factor)
         ! --------------------------------------------------

         ! Calculate GAMMA PAR only during day
         IF ( SUNCOS > 0d0 ) THEN
            
            GAMMA_PAR = GET_GAMMA_PAR_PCEEA(I, J, 
     &           Q_DIR_2, Q_DIFF_2,
     &           PARDR_15_AVG(I,J),  PARDF_15_AVG(I,J) )

         ELSE             
 
            ! If night
            GAMMA_PAR = 0.d0  
            
         ENDIF 

         ! --------------------------------------------------
         ! GAMMA_T_LI (temperature activity factor for light-independent fraction)
         ! --------------------------------------------------
         GAMMA_T_LI = GET_GAMMA_T_LI( TS, BETA )

         ! --------------------------------------------------
         ! GAMMA_T_LD (temperature activity factor for light-dependent fraction)
         ! --------------------------------------------------
         GAMMA_T_LD = GET_GAMMA_T_LD( TS, T_15_AVG(I,J), 
     &        T_DAILY(I,J), CT1, CEO )

         ! --------------------------------------------------
         ! GAMMA_LAI (leaf area index activity factor)
         ! --------------------------------------------------
         GAMMA_LAI = GET_GAMMA_LAI( MISOLAI(I,J), BIDIR )


         ! --------------------------------------------------
         ! GAMMA_AGE (leaf age activity factor)
         ! --------------------------------------------------
         GAMMA_AGE = GET_GAMMA_AGE( MISOLAI(I, J), PMISOLAI(I, J), 
     &        D_BTW_M, T_15_AVG(I,J), ANEW, AGRO, AMAT, AOLD )

         ! --------------------------------------------------
         ! GAMMA_SM (soil moisture activity factor)
         ! --------------------------------------------------
         GAMMA_SM = GET_GAMMA_SM( I, J, CMPD )

      ELSE

         ! set activity factors to zero
         GAMMA_PAR  = 0d0
         GAMMA_T_LI = 0d0
         GAMMA_T_LD = 0d0
         GAMMA_LAI  = 0d0
         GAMMA_AGE  = 0d0
         GAMMA_SM   = 0d0

      END IF
    
      ! Emission is the product of all of these.
      ! Units here are atomsC/box or molec/box as appropriate for the compound.
      ! Normalization factor ensures product of GAMMA values is 1.0 under standard
      ! conditions. 
      MEGAN_EMIS = NORM_FAC(1) * AEF * GAMMA_AGE * GAMMA_SM * GAMMA_LAI
     &     * ((1d0 - LDF) * GAMMA_T_LI + (LDF * GAMMA_PAR * GAMMA_T_LD))

      ! return to calling program
      END FUNCTION GET_MEGAN_EMISSIONS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_MEGAN_PARAMS
!
! !DESCRIPTION: Subroutine GET\_MEGAN\_PARAMS returns the emission parameters for
!  each MEGAN compound needed to compute emissions. Called from GET_MEGAN_EMISSIONS.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MEGAN_PARAMS( CPD, BTA, LIDF, C_T1,
     &     C_EO, A_NEW, A_GRO, A_MAT, A_OLD, BI_DIR )
!
! !USES:
!
! !INPUT AND OUTPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN) :: CPD    ! Compound name
      REAL*8,        INTENT(INOUT) :: BTA    ! Beta coefficient for temperature activity 
                                       ! factor for light-independent fraction
      REAL*8,        INTENT(INOUT) :: LIDF   ! Light-dependent fraction of emissions
      REAL*8,        INTENT(INOUT) :: C_T1   ! CT1 parameter for temperature activity 
                                       ! factor for light-dependent fraction
      REAL*8,        INTENT(INOUT) :: C_EO   ! Ceo parameter for temperature activity 
                                       ! factor for light-dependent fraction
      REAL*8,        INTENT(INOUT) :: A_NEW  ! Relative emission factor (new leaves)
      REAL*8,        INTENT(INOUT) :: A_GRO  ! Relative emission factor (growing leaves)
      REAL*8,        INTENT(INOUT) :: A_MAT  ! Relative emission factor (mature leaves)
      REAL*8,        INTENT(INOUT) :: A_OLD  ! Relative emission factor (old leaves)
      LOGICAL,       INTENT(INOUT) :: BI_DIR ! Logical flag to indicate bidirectional exchange
!                                                                             .
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, (GMD 2012) and associated MEGANv2.1 source code 
! 
! !REVISION HISTORY: 
!  (1 ) Created by dbm 07/2012
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !-----------------------------------------------------------------
      ! GET_MEGAN_PARAMS
      !-----------------------------------------------------------------  

      ! Initialize values
      BTA    = 0d0           
      LIDF   = 0d0    
      C_T1   = 0d0    
      C_EO   = 0d0     
      A_NEW  = 0d0  
      A_GRO  = 0d0   
      A_MAT  = 0d0    
      A_OLD  = 0d0   
      BI_DIR = .FALSE. 

      ! ---------------------------
      ! Note that not all the above compounds are used in standard chemistry
      ! simulations, but they are provided here for future incorporation or 
      ! specialized applications. More compounds can be added as needed
      ! by adding the corresponding CPD name and the appropriate paramaters.
      ! (dbm, 01/2013)
      ! ---------------------------

      ! Isoprene, MBO
      IF ( TRIM(CPD) == 'ISOP' .OR.
     &     TRIM(CPD) == 'MBOX' ) THEN
         BTA    = 0.13d0  ! Not actually used for ISOP, MBO
         LIDF   = 1.d0    
         C_T1   = 95d0    
         C_EO   = 2d0     
         A_NEW  = 0.05d0  
         A_GRO  = 0.6d0   
         A_MAT  = 1.d0    
         A_OLD  = 0.9d0   
         BI_DIR = .FALSE. 

      ! Myrcene, sabinene, alpha-pinene
      ELSE IF ( TRIM(CPD) == 'MYRC' .OR. 
     &          TRIM(CPD) == 'SABI' .OR.
     &          TRIM(CPD) == 'APIN' ) THEN
         BTA    = 0.10d0  
         LIDF   = 0.6d0   
         C_T1   = 80d0    
         C_EO   = 1.83d0  
         A_NEW  = 2d0     
         A_GRO  = 1.8d0   
         A_MAT  = 1.d0    
         A_OLD  = 1.05d0  
         BI_DIR = .FALSE. 

      ! Limonene, 3-carene, beta-pinene
      ELSE IF ( TRIM(CPD) == 'LIMO' .OR. 
     &          TRIM(CPD) == 'CARE' .OR.
     &          TRIM(CPD) == 'BPIN' ) THEN
         BTA    = 0.10d0  
         LIDF   = 0.2d0   
         C_T1   = 80d0    
         C_EO   = 1.83d0  
         A_NEW  = 2d0     
         A_GRO  = 1.8d0   
         A_MAT  = 1.d0    
         A_OLD  = 1.05d0  
         BI_DIR = .FALSE. 
       
      ! t-beta-ocimene
      ELSE IF ( TRIM(CPD) == 'OCIM' ) THEN
         BTA    = 0.10d0  
         LIDF   = 0.8d0   
         C_T1   = 80d0    
         C_EO   = 1.83d0  
         A_NEW  = 2d0     
         A_GRO  = 1.8d0   
         A_MAT  = 1.d0    
         A_OLD  = 1.05d0  
         BI_DIR = .FALSE. 

      ! Other monoterpenes (lumped)
      ELSE IF ( TRIM(CPD) == 'OMON' ) THEN
         BTA    = 0.10d0  
         LIDF   = 0.4d0   
         C_T1   = 80d0    
         C_EO   = 1.83d0  
         A_NEW  = 2d0     
         A_GRO  = 1.8d0   
         A_MAT  = 1.d0    
         A_OLD  = 1.05d0  
         BI_DIR = .FALSE. 
         
      ! Methanol
      ELSE IF ( TRIM(CPD) == 'MOHX' ) THEN
         BTA    = 0.08d0  
         LIDF   = 0.8d0   
         C_T1   = 60d0    
         C_EO   = 1.6d0   
         A_NEW  = 3.5d0   
         A_GRO  = 3d0     
         A_MAT  = 1.d0    
         A_OLD  = 1.2d0   
         BI_DIR = .FALSE. 

      ! Acetone
      ELSE IF ( TRIM(CPD) == 'ACET' ) THEN
         BTA    = 0.1d0   
         LIDF   = 0.2d0   
         C_T1   = 80d0    
         C_EO   = 1.83d0  
         A_NEW  = 1d0     
         A_GRO  = 1d0     
         A_MAT  = 1d0     
         A_OLD  = 1d0     
         BI_DIR = .FALSE. 
         
      ! Bidirectional VOC: Ethanol, formaldehyde, acetaldehyde, formic acid, acetic acid
      ELSE IF ( TRIM(CPD) == 'ETOH' .OR.
     &          TRIM(CPD) == 'CH2O' .OR. 
     &          TRIM(CPD) == 'ALD2' .OR.
     &          TRIM(CPD) == 'FAXX' .OR.
     &          TRIM(CPD) == 'AAXX' ) THEN
         BTA    = 0.13d0  
         LIDF   = 0.8d0   
         C_T1   = 95d0    
         C_EO   = 2d0     
         A_NEW  = 1d0     
         A_GRO  = 1d0     
         A_MAT  = 1d0     
         A_OLD  = 1d0     
         BI_DIR = .TRUE.  
        
      ! Stress VOCs: ethene, toluene, HCN
      ! There are others species in this category but none are currently 
      ! used in GEOS-Chem
      ELSE IF ( TRIM(CPD) == 'C2H4' .OR.
     &          TRIM(CPD) == 'TOLU' .OR. 
     &          TRIM(CPD) == 'HCNX' ) THEN
         BTA    = 0.1d0   
         LIDF   = 0.8d0   
         C_T1   = 80d0    
         C_EO   = 1.83d0  
         A_NEW  = 1d0     
         A_GRO  = 1d0     
         A_MAT  = 1d0     
         A_OLD  = 1d0     
         BI_DIR = .FALSE. 
         
      ! Other VOCs: >C2 alkenes
      ! This includes propene, butene and very minor contribution from larger alkenes
      ELSE IF ( TRIM(CPD) == 'PRPE' ) THEN
         BTA    = 0.1d0   
         LIDF   = 0.2d0   
         C_T1   = 80d0    
         C_EO   = 1.83d0  
         A_NEW  = 1d0     
         A_GRO  = 1d0     
         A_MAT  = 1d0     
         A_OLD  = 1d0     
         BI_DIR = .FALSE. 
         
      ! Calls for any other MEGAN compounds (e.g. sesquiterpenes, CO, etc.) can
      ! be added following the above format based on the parameters in Guenther 2012 or the 
      ! MEGAN source code (dbm, 6/21/2012).
      ELSE

         CALL ERROR_STOP( 'Invalid compound name for MEGAN emissions', 
     &        ' (megan_mod.f)' )
         
      ENDIF
     
      END SUBROUTINE GET_MEGAN_PARAMS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_MEGAN_AEF
!
! !DESCRIPTION: Returns the appropriate AEF value for a given compound
!  and grid square.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MEGAN_AEF( I, J, CPD )
     &              RESULT( EMFAC )
!
! !USES:
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN) :: I,  J  ! Lon & lat indices 
      CHARACTER(LEN=*), INTENT(IN) :: CPD    ! Compound name
      
!
! !RETURN VALUE:
!
      REAL*8              :: EMFAC             ! MEGAN base emission factor
                                               ! (molec/box or atomsC/box)
                                               ! for grid cell (I,J)
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2012, MEGANv2.1 source code
! 
! !REVISION HISTORY: 
!  (1 ) Created 11/2012 by dbm
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !-----------------------------------------------------------------
      ! GET_MEGAN_AEF
      !-----------------------------------------------------------------  


      ! Find appropriate tracer
      SELECT CASE ( TRIM( CPD ) )

      CASE ( 'ISOP' )
         EMFAC = AEF_ISOP(I,J)
      CASE ( 'MBOX' )
         EMFAC = AEF_MBOX(I,J)
      CASE ( 'MYRC' )
         EMFAC = AEF_MYRC(I,J)
      CASE ( 'SABI' )
         EMFAC = AEF_SABI(I,J)
      CASE ( 'APIN' )
         EMFAC = AEF_APIN(I,J)
      CASE ( 'LIMO' )
         EMFAC = AEF_LIMO(I,J)
      CASE ( 'CARE' )
         EMFAC = AEF_CARE(I,J)
      CASE ( 'BPIN' )
         EMFAC = AEF_BPIN(I,J)
      CASE ( 'OCIM' )
         EMFAC = AEF_OCIM(I,J)
      CASE ( 'OMON' )
         EMFAC = AEF_OMON(I,J)
      CASE ( 'MOHX' )
         EMFAC = AEF_MOHX(I,J)
      CASE ( 'ACET' )
         EMFAC = AEF_ACET(I,J)
      CASE ( 'ETOH' )
         EMFAC = AEF_ETOH(I,J)
      CASE ( 'CH2O' )
         EMFAC = AEF_CH2O(I,J)
      CASE ( 'ALD2' )
         EMFAC = AEF_ALD2(I,J)
      CASE ( 'FAXX' )
         EMFAC = AEF_FAXX(I,J)
      CASE ( 'AAXX' )
         EMFAC = AEF_AAXX(I,J)
      CASE ( 'C2H4' )
         EMFAC = AEF_C2H4(I,J)
      CASE ( 'TOLU' )
         EMFAC = AEF_TOLU(I,J)
      CASE ( 'HCNX' )
         EMFAC = AEF_HCNX(I,J)
      CASE ( 'PRPE' )
         EMFAC = AEF_PRPE(I,J)
      CASE DEFAULT     
         CALL ERROR_STOP( 'Invalid compound name for MEGAN AEF', 
     &        ' (megan_mod.f)' )
      END SELECT
         
      ! return to calling program
      END FUNCTION GET_MEGAN_AEF
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_GAMMA_PAR_PCEEA
!
! !DESCRIPTION: Computes the PCEEA gamma activity factor with sensitivity 
!  to LIGHT.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_PAR_PCEEA( I , J , Q_DIR_2, Q_DIFF_2 ,
     &                            PARDR_AVG_SIM ,  PARDF_AVG_SIM  )
     &              RESULT( GAMMA_P_PCEEA )
!
! !USES:
!
      USE TIME_MOD,   ONLY : GET_DAY_OF_YEAR
      USE TIME_MOD,   ONLY : GET_LOCALTIME
      USE GRID_MOD,   ONLY : GET_YMID  
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I,  J             ! Lon & lat indices 
      REAL*8,  INTENT(IN) :: PARDR_AVG_SIM     ! Average direct PAR [W/m2]
      REAL*8,  INTENT(IN) :: PARDF_AVG_SIM     ! Average diffuse PAR [W/m2]
      REAL*8,  INTENT(IN) :: Q_DIR_2           ! Direct PAR [umol/m2/s]
      REAL*8,  INTENT(IN) :: Q_DIFF_2          ! Diffuse PAR [umol/m2/s]
!
! !RETURN VALUE:
!
      REAL*8              :: GAMMA_P_PCEEA     ! GAMMA factor for light
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 user guide
! 
! !REVISION HISTORY: 
!  (1 ) Here PAR*_AVG_SIM is the average light conditions over the simulation 
!       period. I've set this = 10 days to be consistent with temperature & as 
!       outlined in Guenther et al, 2006. (mpb,2009)
!  (2 ) Code was taken & adapted directly from the MEGAN v2.1 source code.
!       (mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8              :: LUT , LAT 
      REAL*8              :: mmPARDR_DAILY
      REAL*8              :: mmPARDF_DAILY
      REAL*8              :: PAC_DAILY, PAC_INSTANT, C_PPFD
      REAL*8              :: PTOA, PHI
      REAL*8              :: BETA,   SINbeta 
      INTEGER             :: DOY 
      REAL*8              :: AAA, BBB

      !-----------------------------------------------------------------
      ! Compute GAMMA_PAR_PCEEA
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
      LUT   = GET_LOCALTIME( I, J, 1 )
      LAT   = GET_YMID( I, J, 1 )

      ! Get solar elevation angle
      SINbeta      =  SOLAR_ANGLE( DOY , LUT , LAT  )
      BETA         =  ASIN( SINbeta ) * RAD2D       

      IF ( SINbeta .LE. 0.0d0 ) THEN

         GAMMA_P_PCEEA = 0.0d0

      ELSEIF ( SINbeta .GT. 0.0d0 ) THEN       

         ! PPFD at top of atmosphere
         PTOA    = 3000.0d0 + 99.0d0 * 
     &             COS( 2.d0 * 3.14d0 *( DOY - 10 ) / 365 )

         ! Above canopy transmission
         PHI     = PAC_INSTANT / ( SINbeta * PTOA )

         ! Work out gamma P
         BBB     = 1.0d0 + 0.0005d0 *( PAC_DAILY - 400.0d0  ) 
         AAA     = ( 2.46d0 * BBB * PHI ) - ( 0.9d0 * PHI**2 )

         GAMMA_P_PCEEA = SINbeta * AAA

      ENDIF

       ! Screen unforced errors. IF solar elevation angle is 
       ! less than 1 THEN gamma_p can not be greater than 0.1.
       IF ( BETA .LT. 1.0 .AND. GAMMA_P_PCEEA .GT. 0.1) THEN
          GAMMA_P_PCEEA  = 0.0
       ENDIF
   
      ! Prevent negative values
      GAMMA_P_PCEEA = MAX( GAMMA_P_PCEEA , 0d0 )

      ! return to calling program
      END FUNCTION GET_GAMMA_PAR_PCEEA
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gamma_t_li
!
! !DESCRIPTION: Function GET\_GAMMA\_T\_LI computes the temperature 
!  activity factor (GAMMA\_T_LI) for the light-independent fraction of emissions
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_T_LI( T, BETA ) RESULT( GAMMA_T_LI )
!
! !INPUT PARAMETERS: 
!
      ! Current leaf temperature [K], the surface air temperature field (TS) 
      ! is assumed equivalent to the leaf temperature over forests.
      REAL*8, INTENT(IN) :: T        

      ! Temperature factor per species
      REAL*8, INTENT(IN) :: BETA     
!
! !RETURN VALUE:
!
      REAL*8             :: gamma_T_LI  ! Activity factor for the light-independent
                                        ! fraction of emissions
!
! !REMARKS:
!  GAMMA_T =  exp[Beta*(T - T_Standard)]
!                                                                             .
!             where Beta   = temperature dependent parameter
!                   Ts     = standard temperature (normally 303K, 30C)
!                                                                             .
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
!
! !REVISION HISTORY: 
!  (1 ) Original code by Michael Barkley (2009).
!       Note: If T = Ts  (i.e. standard conditions) then GAMMA_T = 1 
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  (2 ) Modified to GAMMA_T_LI (dbm, 6/21/2012)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Standard reference temperature [K]
      REAL*8, PARAMETER   :: T_STANDARD = 303.d0

      !=================================================================
      ! GET_GAMMAT_LI
      !================================================================= 

      GAMMA_T_LI = EXP( BETA * ( T - T_STANDARD ) )

      END FUNCTION GET_GAMMA_T_LI
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_GAMMA_T_LD
!
! !DESCRIPTION: Function GET\_GAMMA\_T\_LD computes the temperature 
!  sensitivity for the light-dependent fraction of emissions.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_T_LD( T, PT_15, PT_1, CT1, CEO ) 
     &     RESULT( GAMMA_T_LD )
!
! !INPUT PARAMETERS: 
!

      ! Current leaf temperature, the surface air temperature field (TS) 
      ! is assumed equivalent to the leaf temperature over forests.
      REAL*8,  INTENT(IN) :: T 

       ! Average leaf temperature over the past 15 days
      REAL*8,  INTENT(IN) :: PT_15

      ! Average leaf temperature over the past arbitray day(s).
      ! This is not used at present
      REAL*8,  INTENT(IN) :: PT_1

      ! Compound-specific parameters for light-dependent temperature activity
      ! factor (dbm, 6/21/2012)
      REAL*8,  INTENT(IN) :: CT1, CEO
      
!
! !RETURN VALUE:
!
      ! Temperature activity factor for the light-dependent fraction of 
      ! emissions
      REAL*8              :: GAMMA_T_LD
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995
!  (2 ) Guenther et al, 2006
!  (3 ) Guenther et al, MEGAN v2.1 user mannual 2007-08
!  (4 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
! 
! 
! !REVISION HISTORY: 
!  (1 ) Includes the latest MEGAN v2.1 temperature algorithm (mpb, 2009).
!       Note, this temp-dependence is the same for the PCEEA & hybrid models.
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  (2 ) Modified to gamma_t_ld and to permit compound specific parameters
!       CT1 and Ceo (dbm, 6/21/2012)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8              :: C_T,   CT2
      REAL*8              :: E_OPT, T_OPT, X

      ! Ideal gas constant [J/mol/K] / 1000
      REAL*8, PARAMETER   :: R   = 8.314d-3

      !=================================================================
      ! GET_GAMMA_T_LD
      !================================================================= 
      E_OPT = CEO * EXP( 0.08d0 * ( PT_15  - 2.97d2 ) )           
      T_OPT = 3.13D2 + ( 6.0d-1 * ( PT_15 - 2.97d2 ) )
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
      GAMMA_T_LD = MAX( C_T , 0D0 )

      ! Return to calling program
      END FUNCTION GET_GAMMA_T_LD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Function GET\_GAMMA\_LAI computes the gamma exchange activity 
!  factor which is sensitive to leaf area (= GAMMA_LAI).  Called from
!  GET_EMISOP_MEGAN, GET_EMMBO_MEGAN, and GET_EMMONOG_MEGAN.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_LAI( CMLAI, BIDIREXCH ) RESULT( GAMMA_LAI )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN)  :: CMLAI       ! Current month's LAI [cm2/cm2]
      LOGICAL, INTENT(IN)  :: BIDIREXCH   ! Logical flag indicating whether
                                          ! the compound undergoes bidirectional
                                          ! exchange
!
! !RETURN VALUE:
!
      REAL*8              :: GAMMA_LAI
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) Update to publically released (as of 11/2004) MEGAN algorithm and 
!        modified for the standard code by May Fu (11/2004).
!  (3 ) Algorithm is based on the latest MEGAN v2.1 User's Guide (mpb,2009)
!  (4 ) Updated to treat bidirectional exchange compounds appropriately (dbm, 6/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      !-----------------------
      ! Compute GAMMA_LAI
      !-----------------------

      ! Formulation for birectional compounds is as described for 
      ! ALD2 in Millet et al., ACP 2010
      IF ( BIDIREXCH ) THEN
         
          IF ( CMLAI <= 6d0) THEN
               
            ! if lai less than 2:
            IF ( CMLAI <= 2d0 ) THEN
               GAMMA_LAI = 0.5d0 * CMLAI
               
            ! if between 2 and 6:
            ELSE
               GAMMA_LAI = 1d0 - 0.0625 * ( CMLAI - 2d0 )
            END IF

            ! keep at 0.75 for LAI > 6
         ELSE
            GAMMA_LAI = 0.75d0
         END IF

      ! For all other compounds use the standard gamma_lai formulation
      ELSE
         GAMMA_LAI = 0.49d0 * CMLAI / SQRT( 1.d0 + 0.2d0 * CMLAI*CMLAI )
      ENDIF

      END FUNCTION GET_GAMMA_LAI 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Function GET\_GAMMA\_AGE computes the gamma exchange 
!  activity factor which is sensitive to leaf age (= GAMMA\_AGE).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_AGE( CMLAI, PMLAI, DBTWN, TT,
     &                        AN, AG, AM, AO ) 
     &            RESULT( GAMMA_AGE )
!
! !INPUT PARAMETERS: 
!
      REAL*8,           INTENT(IN) :: CMLAI     ! Current month's LAI [cm2/cm2]
      REAL*8,           INTENT(IN) :: PMLAI     ! Previous months LAI [cm2/cm2]
      REAL*8,           INTENT(IN) :: DBTWN     ! Number of days between 
      REAL*8,           INTENT(IN) :: TT        ! Daily average temperature [K]
      REAL*8,           INTENT(IN) :: AN, AG    ! Relative emission factors for
      REAL*8,           INTENT(IN) :: AM, AO    ! new, growing, mature, old leaves
!
! !RETURN VALUE:
!
      REAL*8                       :: GAMMA_AGE ! Activity factor    
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
! 
! !REVISION HISTORY: 
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
!  (7 ) Updated to pass leaf age activity factors as arguments (dbm, 6/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      REAL*8             :: FNEW  ! Fraction of new leaves in canopy
      REAL*8             :: FGRO  ! Fraction of growing leaves
      REAL*8             :: FMAT  ! Fraction of mature leaves
      REAL*8             :: FOLD  ! Fraction of old leaves

      ! TI: number of days after budbreak required to induce emissions
      REAL*8             :: TI   

      ! TM: number of days after budbreak required to reach peak emissions
      REAL*8             :: TM

      !=================================================================
      ! GET_GAMMA_AGE begins here!
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
      ! Compute GAMMA_AGE
      !-----------------------

      IF ( CMLAI == PMLAI ) THEN !(i.e. LAI stays the same) 

         FNEW = 0.d0        
         FGRO = 0.1d0       
         FMAT = 0.8d0        
         FOLD = 0.1d0       

      ELSE IF ( CMLAI > PMLAI ) THEN !(i.e. LAI has increased) 

         ! Calculate Fnew
         IF ( DBTWN > TI ) THEN
            FNEW = ( TI / DBTWN ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE
            FNEW = 1.d0 - ( PMLAI / CMLAI )
         ENDIF

         ! Calculate FMAT
         IF ( DBtwn > TM ) THEN
            FMAT = ( PMLAI / CMLAI ) +
     &             (( DBTWN - TM ) / DBTWN ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE 
            FMAT = ( PMLAI / CMLAI )
         ENDIF

         ! Calculate Fgro and Fold
         FGRO = 1.d0 - FNEW - FMAT
         FOLD = 0.d0
         

      ELSE ! This is the case if  PMLAI > CMLAI (i.e. LAI has decreased) 

         FNEW = 0.d0
         FGRO = 0.d0
         FOLD = ( PMLAI - CMLAI ) / PMLAI
         FMAT = 1.d0 - FOLD
                  
      ENDIF
 
      ! Age factor
      GAMMA_AGE = FNEW * AN + FGRO * AG +
     &        FMAT * AM + FOLD * AO

      ! Prevent negative values
      GAMMA_AGE = MAX( GAMMA_AGE , 0D0 )

      ! return to calling program
      END FUNCTION GET_GAMMA_AGE

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Function GET\_GAMMA\_SM computes activity factor for soil moisture
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_SM( I, J, CMPD )
     &            RESULT( GAMMA_SM )
!
! !USES:
!
      USE DAO_MOD,     ONLY : GWETROOT ! For soil moisture
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN)  :: I, J     ! GEOS-Chem lon & lat indices
      CHARACTER(LEN=*), INTENT(IN)  :: CMPD     ! Compound name (dbm, 6/21/2012)
      
! !RETURN VALUE:
!
      REAL*8                        :: GAMMA_SM ! Activity factor    
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, ACP 2006
!  (2 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
! 
! !REVISION HISTORY: 
!  (1 ) Created by dbm (6/2012). We are not currently using a soil moisture effect
!       for isoprene. For all compounds other than acetaldehyde and ethanol,
!       gamma_sm =1 presently.
!------------------------------------------------------------------------------
!BOC
!

      !=================================================================
      ! GET_GAMMA_SM begins here!
      !================================================================= 
      
      ! By default gamma_sm is 1.0
      GAMMA_SM = 1d0

      ! GWETROOT is only defined for GEOS-5
#if   defined( GEOS_5 )

      IF ( TRIM( CMPD ) == 'ALD2' .OR. TRIM ( CMPD ) == 'ETOH' ) THEN
         
         ! GWETROOT = degree of saturation or wetness in the root-zone
         ! (top meter of soil). This is defined as the ratio of the volumetric soil
         ! moisture to the porosity. We use a soil moisture activity factor 
         ! for ALD2 to account for stimulation of emission by flooding.
         ! (Millet et al., ACP 2010)
         ! Constant value of 1.0 for GWETROOT = 0-0.9, increasing linearly to 3.0 at GWETROOT =1.
         GAMMA_SM = MAX( 20d0* GWETROOT(I,J) -17d0, 1.d0)

      ENDIF

#endif


      ! return to calling program
      END FUNCTION GET_GAMMA_SM

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_NORM_FAC
!
! !DESCRIPTION: Function CALC\_NORM\_FAC calculates the normalization factor 
!  needed to compute emissions. Called from GET_MEGAN_EMISSIONS.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_NORM_FAC

!
! !USES:
!
! !INPUT PARAMETERS: 
!
                                       ! None
!                                                                             .
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, (GMD 2012) and associated MEGANv2.1 source code 
! 
! !REVISION HISTORY: 
!  (1 ) Created by dbm 11/2012. We calculate only 1 normalization factor for all
!       compounds based on the isoprene gamma values. Formally there should be a 
!       different normalization factor for each compound, but we are following Alex
!       Guenther's approach here and the MEGAN source code.
!       "Hi Dylan, sorry for being so slow to get back to you.
!        Since the change is only a few percent or less, I didn't
!        bother to assign a different normalization factor to each
!        compound.  Since the MEGAN canopy environment model also 
!        has 8 different canopy types (tropical broadleaf tree, 
!        conifer tree, etc.) then to be correct we should have a 
!        different CCE for each canopy type for each compound class
!        (which would be 160 slightly different values of CCE)."
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      REAL*8             :: PAC_DAILY, PHI, BBB, AAA, GAMMA_P_STANDARD
      REAL*8             :: GAMMA_T_LI_STANDARD
      REAL*8             :: GAMMA_SM_STANDARD
      REAL*8             :: CMLAI, GAMMA_LAI_STANDARD
      REAL*8             :: GAMMA_AGE_STANDARD
      REAL*8             :: PT_15, T, R, CEO, CT1, E_OPT, T_OPT, CT2, X
      REAL*8             :: GAMMA_T_LD_STANDARD
      REAL*8             :: LDF, GAMMA_STANDARD

      !-----------------------------------------------------------------
      ! CALC_NORM_FAC
      !-----------------------------------------------------------------  

      ! -----------------
      ! GAMMA_P for standard conditions
      ! -----------------
      ! Based on Eq. 11b from Guenther et al., 2006
      ! Using standard conditions of phi = 0.6, solar angle of 60 deg,
      ! and P_daily = 400
      ! Note corrigendum for Eq. 11b in that paper, should be a
      ! minus sign before the 0.9.
      PAC_DAILY = 400d0
      PHI       = 0.6d0
      BBB       = 1.0d0 + 0.0005d0 *( PAC_DAILY - 400.0d0  ) 
      AAA       = ( 2.46d0 * BBB * PHI ) - ( 0.9d0 * PHI**2 )
      ! sin(60) = 0.866
      GAMMA_P_STANDARD = 0.866d0 * AAA

      ! -----------------
      ! GAMMA_T_LI for standard conditions
      ! -----------------
      ! gamma_t_li = EXP( Beta * ( T - T_Standard ) )
      ! This is 1.0 for T = T_Standard
      GAMMA_T_LI_STANDARD = 1d0

      ! -----------------
      ! GAMMA_SM for standard conditions
      ! -----------------
      ! Standard condition is soil moisture = 0.3 m^3/m^3
      ! GAMMA_SM = 1.0 for all compounds under this condition
      GAMMA_SM_STANDARD = 1d0

      ! -----------------
      ! GAMMA_LAI for standard conditions 
      ! -----------------
      ! Standard condition is LAI = 5
      CMLAI = 5d0
      GAMMA_LAI_STANDARD = 0.49d0 * 
     &     CMLAI / SQRT( 1.d0 + 0.2d0 * CMLAI*CMLAI )
     
      ! -----------------
      ! GAMMA_AGE for standard conditions
      ! -----------------
      ! Standard condition is 0% new, 10% growing, 80% mature, 10% old foliage
      ! Isoprene uses A_NEW  = 0.05d0, A_GRO  = 0.6d0, A_MAT  = 1.d0, A_OLD  = 0.9d0   
      GAMMA_AGE_STANDARD = 0.1d0*0.6d0 + 0.8d0*1d0 + 0.1d0*0.9d0

      ! -----------------
      ! GAMMA_T_LD for standard conditions
      ! -----------------
      ! Standard condition is
      ! PT_15 = average leaf temp over past 24-240 hours = 297K
      ! T = air temperature = 303K
      PT_15 = 297d0
      T     = 303d0
      R     = 8.314d-3
      ! parameters for isoprene
      CEO = 2d0
      CT1 = 95d0

      E_OPT = CEO * EXP( 0.08d0 * ( PT_15  - 2.97d2 ) )           
      T_OPT = 3.13d2 + ( 6.0d-1 * ( PT_15 - 2.97d2 ) )
      CT2   = 200d0

      ! Variable related to temperature 
      X     = ( 1.d0/T_OPT - 1.d0/T ) / R

      GAMMA_T_LD_STANDARD = E_OPT * CT2 * EXP( CT1 * X ) / 
     &     ( CT2 - CT1 * ( 1.d0 - EXP( CT2 * X ) ) )

      ! -----------------
      ! Overall GAMMA_STANDARD
      ! -----------------
      ! LDF = 1d0 for isoprene
      LDF = 1d0
      GAMMA_STANDARD = 
     &     GAMMA_AGE_STANDARD * GAMMA_SM_STANDARD * GAMMA_LAI_STANDARD
     &     * ((1d0 - LDF) * GAMMA_T_LI_STANDARD 
     &     + (LDF * GAMMA_P_STANDARD * GAMMA_T_LD_STANDARD))
      ! This ends up being 1.0101081.
      
      NORM_FAC = 1d0 / GAMMA_STANDARD

      END SUBROUTINE CALC_NORM_FAC
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Function SOLAR\_ANGLE computes the local solar angle for a 
!  given day of year, latitude and longitude (or local time).  Called from  
!  routine GAMMA\_PAR\_PCEEA.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SOLAR_ANGLE( DOY, SHOUR, LAT ) RESULT( SINbeta )
!
! !USES:
!

!
! !INPUT PARAMETERS: 
!
      ! Arguments
      INTEGER, INTENT(IN) :: DOY       ! Day of year
      REAL*8,  INTENT(IN) :: SHOUR     ! Local time 
      REAL*8,  INTENT(IN) :: LAT       ! Latitude
!
! !RETURN VALUE:
!
      REAL*8              :: SINbeta   ! Sin of the local solar angle
!
! !REMARKS:
!  References (see above for full citations):
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN v2.1 user mannual 2007-09
! 
! !REVISION HISTORY: 
!  (1 ) This code was taken directly from the MEGAN v2.1 source code.(mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_AEF
!
! !DESCRIPTION: Subroutine CALC\_AEF reads Emission Factors for all
!  biogenic VOC species from disk.  CALC\_AEF is called from 
!  "SUBROUTINE INIT_MEGAN".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_AEF
!
! !USES:
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_RES_EXT, READ_BPCH2
      USE TIME_MOD,       ONLY : GET_TS_EMIS
      USE GRID_MOD,       ONLY : GET_AREA_M2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR

!
! !REMARKS:
!  Reference: (5 ) Guenther et al, 2004 
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) AEF detailed in the latest MEGAN User's Guide (tmf, 11/19/04)
!  (3 ) Bug fix (tmf, 11/30/04)
!  (4 ) Now reads 1x1 files and regrids to current resolution (bmy, 10/24/05)
!  (5 ) Uses new v2.1 emission factors maps for isoprene, MBO and 7 monoterpene
!     species, download in 2009. (mpb,2009)
!  (6 ) Now use 2.1 emission factors for isoprene, MBO, and 7 monoterpenes. EFs for
!       other compounds are computed by reading in the PFT fractions and multiplying
!       the fractions by corresponding EF values. (dbm, 11/2012)
!  (7 ) Also, now read in the EF values already gridded to model resolution.
!       (dbm, 11/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                 :: I, J, P, ARR_IND
      REAL*4                  :: ARRAY   (IIPAR, JJPAR, 1)
      REAL*4                  :: ARRAY_16(IIPAR, JJPAR, 16)
      REAL*8                  :: DTSRCE, AREA_M2, FACTOR
      CHARACTER(LEN=255)      :: FILENAME
      REAL*8,  PARAMETER      :: AVO = 6.022D23
      REAL*8                  :: PFT_EF_OMON(15), PFT_EF_MOHX(15)
      REAL*8                  :: PFT_EF_ACET(15), PFT_EF_BIDR(15)
      REAL*8                  :: PFT_EF_STRS(15), PFT_EF_OTHR(15)
! --->
! dbm, compute EF maps for a-pinene and myrcene as well since there seems to 
! be an issue with the EF maps for these species provided on the MEGAN
! data portal
      REAL*8                  :: PFT_EF_APIN(15), PFT_EF_MYRC(15)
! <---
      REAL*8                  :: EM_FRAC_ALD2(15), EM_FRAC_ETOH(15)
      REAL*8                  :: EM_FRAC_FAXX(15), EM_FRAC_AAXX(15)
      REAL*8                  :: EM_FRAC_CH2O(15)
      
      !=================================================================
      ! CALC_AEF begins here!
      !=================================================================

      ! ---------------------------
      ! Note that not all these compounds are used in standard chemistry
      ! simulations, but they are provided here for future incorporation or 
      ! specialized applications. More compounds can be added as needed
      ! by adding the corresponding PFT-specific emission factors and 
      ! emission category fraction.
      ! (dbm, 01/2013)
      ! ---------------------------

      ! Emission timestep [min]
      DTSRCE = GET_TS_EMIS()
      
      !---------------------------------------------
      ! Initialize AEF arrays for safety
      !---------------------------------------------
      AEF_ISOP(:,:) = 0d0
      AEF_MBOX(:,:) = 0d0
      AEF_MYRC(:,:) = 0d0
      AEF_SABI(:,:) = 0d0
      AEF_APIN(:,:) = 0d0
      AEF_LIMO(:,:) = 0d0
      AEF_CARE(:,:) = 0d0
      AEF_BPIN(:,:) = 0d0
      AEF_OCIM(:,:) = 0d0
      AEF_OMON(:,:) = 0d0
      AEF_MOHX(:,:) = 0d0
      AEF_ACET(:,:) = 0d0
      AEF_ETOH(:,:) = 0d0
      AEF_CH2O(:,:) = 0d0
      AEF_ALD2(:,:) = 0d0
      AEF_FAXX(:,:) = 0d0
      AEF_AAXX(:,:) = 0d0
      AEF_C2H4(:,:) = 0d0
      AEF_TOLU(:,:) = 0d0
      AEF_HCNX(:,:) = 0d0
      AEF_PRPE(:,:) = 0d0



      !---------------------------------------------
      ! Read in ISOPRENE Emission Factors
      !---------------------------------------------

      ! File name
      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
     &     'MEGAN2.1_EF_ISOP_ug.m2.h_geos.' // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
100   FORMAT( '     - CALC_AEF: Reading ', a )

      ! Read data [ug compound/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
     &                 0d0,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8
      AEF_ISOP = ARRAY(:,:,1)

      !---------------------------------------------
      ! Read in MBO Emission Factors
      !---------------------------------------------

      ! File name
      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
     &     'MEGAN2.1_EF_MBO_ug.m2.h_geos.' // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [ug compound/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
     &                 0d0,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8
      AEF_MBOX = ARRAY(:,:,1)

      !---------------------------------------------
      ! Read in A-Pinene Emission Factors
      !---------------------------------------------

! ---> Now compute EF maps for a-pinene and myrcene (dbm, 12/2012)
c$$$      ! File name
c$$$      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
c$$$     &     'MEGAN2.1_EF_APIN_ug.m2.h_geos.' // GET_RES_EXT()
c$$$
c$$$      ! Echo info
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$
c$$$      ! Read data [ug compound/m2/hr] 
c$$$      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
c$$$     &                 0d0,       IIPAR,     JJPAR,     
c$$$     &                 1,         ARRAY,     QUIET=.TRUE. )
c$$$
c$$$      ! Cast to REAL*8
c$$$      AEF_APIN = ARRAY(:,:,1)
! <---

      !---------------------------------------------
      ! Read in B-Pinene Emission Factors
      !---------------------------------------------

      ! File name
      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
     &     'MEGAN2.1_EF_BPIN_ug.m2.h_geos.' // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [ug compound/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
     &                 0d0,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8
      AEF_BPIN = ARRAY(:,:,1)

      !---------------------------------------------
      ! Read in Carene Emission Factors
      !---------------------------------------------

      ! File name
      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
     &     'MEGAN2.1_EF_CARE_ug.m2.h_geos.' // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [ug compound/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
     &                 0d0,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8
      AEF_CARE = ARRAY(:,:,1)

      !---------------------------------------------
      ! Read in Limonene Emission Factors
      !---------------------------------------------

      ! File name
      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
     &     'MEGAN2.1_EF_LIMO_ug.m2.h_geos.' // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [ug compound/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
     &                 0d0,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8
      AEF_LIMO = ARRAY(:,:,1)

      !---------------------------------------------
      ! Read in Myrcene Emission Factors
      !---------------------------------------------

! ---> Now compute EF maps for a-pinene and myrcene (dbm, 12/2012)
      ! File name
c$$$      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
c$$$     &     'MEGAN2.1_EF_MYRC_ug.m2.h_geos.' // GET_RES_EXT()
c$$$
c$$$      ! Echo info
c$$$      WRITE( 6, 100 ) TRIM( FILENAME )
c$$$
c$$$      ! Read data [ug compound/m2/hr] 
c$$$      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
c$$$     &                 0d0,       IIPAR,     JJPAR,     
c$$$     &                 1,         ARRAY,     QUIET=.TRUE. )
c$$$
c$$$      ! Cast to REAL*8
c$$$      AEF_MYRC = ARRAY(:,:,1)
! <---

      !---------------------------------------------
      ! Read in Ocimene Emission Factors
      !---------------------------------------------

      ! File name
      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
     &     'MEGAN2.1_EF_OCIM_ug.m2.h_geos.' // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [ug compound/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
     &                 0d0,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8
      AEF_OCIM = ARRAY(:,:,1)

      !---------------------------------------------
      ! Read in Sabinene Emission Factors
      !---------------------------------------------

      ! File name
      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
     &     'MEGAN2.1_EF_SABI_ug.m2.h_geos.' // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [ug compound/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
     &                 0d0,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8
      AEF_SABI = ARRAY(:,:,1)

      !---------------------------------------------
      ! Read in PFT fractions
      !---------------------------------------------

      ! File name
      FILENAME = TRIM( DATA_DIR ) // TRIM( MEGAN_SUBDIR ) //
     &     'CLM4_PFT.geos.' // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [ug compound/m2/hr] 
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE', -1,   
     &                 0d0,       IIPAR,     JJPAR,     
     &                 1,         ARRAY_16,  QUIET=.TRUE. )

      ! CLM4 PFT coverage (%)
      ! PFT_BARE                = ARRAY_16(:,:,1 )
      ! PFT_NDLF_EVGN_TMPT_TREE = ARRAY_16(:,:,2 )
      ! PFT_NDLF_EVGN_BORL_TREE = ARRAY_16(:,:,3 )
      ! PFT_NDLF_DECD_BORL_TREE = ARRAY_16(:,:,4 )
      ! PFT_BDLF_EVGN_TROP_TREE = ARRAY_16(:,:,5 )
      ! PFT_BDLF_EVGN_TMPT_TREE = ARRAY_16(:,:,6 )
      ! PFT_BDLF_DECD_TROP_TREE = ARRAY_16(:,:,7 )
      ! PFT_BDLF_DECD_TMPT_TREE = ARRAY_16(:,:,8 )
      ! PFT_BDLF_DECD_BORL_TREE = ARRAY_16(:,:,9 )
      ! PFT_BDLF_EVGN_SHRB      = ARRAY_16(:,:,10)
      ! PFT_BDLF_DECD_TMPT_SHRB = ARRAY_16(:,:,11)
      ! PFT_BDLF_DECD_BORL_SHRB = ARRAY_16(:,:,12)
      ! PFT_C3_ARCT_GRSS        = ARRAY_16(:,:,13)
      ! PFT_C3_NARC_GRSS        = ARRAY_16(:,:,14)
      ! PFT_C4_GRSS             = ARRAY_16(:,:,15)
      ! PFT_CROP                = ARRAY_16(:,:,16)

      ! -------------------------------------------
      ! PFT-specific EFs from Guenther et al., 2012
      ! in ug compound/m2/h
      ! PFTs 1-15 in the table correspond to #2-16 
      ! (i.e., excluding bare ground #1) in the above array.
      ! -------------------------------------------
      ! Compound Class EF1 EF2 EF3 EF4 EF5 EF6 EF7 EF8 EF9 EF10 EF11 EF12 EF13 EF14 EF15
      ! --------------------------------------------------------------------------------
      ! Other Monoterp 180 180 170 150 150 150 150 150 110 200  110  5    5    5    5
      ! Methanol       900 900 900 500 900 500 900 900 900 900  900  500  500  500  900
      ! Acetone        240 240 240 240 240 240 240 240 240 240  240  80   80   80   80
      ! Bidirect VOC   500 500 500 500 500 500 500 500 500 500  500  80   80   80   80
      ! Stress VOC     300 300 300 300 300 300 300 300 300 300  300  300  300  300  300
      ! Other VOC      140 140 140 140 140 140 140 140 140 140  140  140  140  140  140  
      ! --------------------------------------------------------------------------------

!     One thing to note is these are net emissions to the canopy atmosphere but not net
!     emissions to the above canopy atmosphere since they don't account for within-canopy
!     deposition. Only an issue for OVOCs.
      
!                     EF1    EF2    EF3    EF4    EF5    
      PFT_EF_OMON = (/180d0, 180d0, 170d0, 150d0, 150d0, 
!                     EF6    EF7    EF8    EF9    EF10
     &                150d0, 150d0, 150d0, 110d0, 200d0, 
!                     EF11   EF12   EF13   EF14   EF15
     &                110d0, 5d0  , 5d0  , 5d0  , 5d0/)

!                     EF1    EF2    EF3    EF4    EF5    
      PFT_EF_MOHX = (/900d0, 900d0, 900d0, 500d0, 900d0, 
!                     EF6    EF7    EF8    EF9    EF10
     &                500d0, 900d0, 900d0, 900d0, 900d0, 
!                     EF11   EF12   EF13   EF14   EF15
     &                900d0, 500d0, 500d0, 500d0, 900d0/)

!                     EF1    EF2    EF3    EF4    EF5    
      PFT_EF_ACET = (/240d0, 240d0, 240d0, 240d0, 240d0, 
!                     EF6    EF7    EF8    EF9    EF10
     &                240d0, 240d0, 240d0, 240d0, 240d0, 
!                     EF11   EF12   EF13   EF14   EF15
     &                240d0, 80d0 , 80d0 , 80d0 , 80d0/)

!                     EF1    EF2    EF3    EF4    EF5    
      PFT_EF_BIDR = (/500d0, 500d0, 500d0, 500d0, 500d0, 
!                     EF6    EF7    EF8    EF9    EF10
     &                500d0, 500d0, 500d0, 500d0, 500d0, 
!                     EF11   EF12   EF13   EF14   EF15
     &                500d0, 80d0 , 80d0 , 80d0 , 80d0/)

!                     EF1    EF2    EF3    EF4    EF5    
      PFT_EF_STRS = (/300d0, 300d0, 300d0, 300d0, 300d0, 
!                     EF6    EF7    EF8    EF9    EF10
     &                300d0, 300d0, 300d0, 300d0, 300d0, 
!                     EF11   EF12   EF13   EF14   EF15
     &                300d0, 300d0, 300d0, 300d0, 300d0/)

!                     EF1    EF2    EF3    EF4    EF5    
      PFT_EF_OTHR = (/140d0, 140d0, 140d0, 140d0, 140d0, 
!                     EF6    EF7    EF8    EF9    EF10
     &                140d0, 140d0, 140d0, 140d0, 140d0, 
!                     EF11   EF12   EF13   EF14   EF15
     &                140d0, 140d0, 140d0, 140d0, 140d0/)

! ---> Now compute EFs for a-pinene and myrcene as well (dbm, 12/2012)
!                     EF1    EF2    EF3    EF4    EF5    
      PFT_EF_APIN = (/500d0, 500d0, 510d0, 600d0, 400d0, 
!                     EF6    EF7    EF8    EF9    EF10
     &                600d0, 400d0, 400d0, 200d0, 300d0, 
!                     EF11   EF12   EF13   EF14   EF15
     &                200d0, 2d0,   2d0,   2d0,   2d0/)

!                     EF1    EF2    EF3    EF4    EF5    
      PFT_EF_MYRC = (/70d0,  70d0,  60d0,  80d0,  30d0, 
!                     EF6    EF7    EF8    EF9    EF10
     &                80d0,  30d0,  30d0,  30d0,  50d0, 
!                     EF11   EF12   EF13   EF14   EF15
     &                30d0,  0.3d0, 0.3d0, 0.3d0, 0.3d0/)
! <---


      ! Other monoterpenes, methanol, acetone, MBO are each 100% of thier
      ! respective categories. The VOCs within the stress category each 
      ! account for a specific fraction of emissions across all PFTs
      ! (ethene 58%, toluene 3%, HCN 1.5%). The VOCs within the 
      ! other category also account for a given fraction of emissions
      ! across all PFTs (propene 48%, butene 24%, other alkenes 0.2%). But
      ! VOCs in the bidirectional category account for a different amount of
      ! the total flux for the different PFTs. So in this case we define a
      ! vector containing these fractions.
      
      ! Acetaldehyde: 40% of bidirectional category flux, except 25%
      ! for grasses and crops
      EM_FRAC_ALD2 = (/0.40d0, 0.40d0, 0.40d0, 0.40d0, 0.40d0,
     &                 0.40d0, 0.40d0, 0.40d0, 0.40d0, 0.40d0,
     &                 0.40d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0/)

      ! Ethanol: 40% of bidirectional category flux, except 25%
      ! for grasses and crops
      EM_FRAC_ETOH = (/0.40d0, 0.40d0, 0.40d0, 0.40d0, 0.40d0,
     &                 0.40d0, 0.40d0, 0.40d0, 0.40d0, 0.40d0,
     &                 0.40d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0/)

      ! Formic acid: 6% of bidirectional category flux, except 15%
      ! for grasses and crops
      EM_FRAC_FAXX = (/0.06d0, 0.06d0, 0.06d0, 0.06d0, 0.06d0,
     &                 0.06d0, 0.06d0, 0.06d0, 0.06d0, 0.06d0,
     &                 0.06d0, 0.15d0, 0.15d0, 0.15d0, 0.15d0/)

      ! Acetic acid: 6% of bidirectional category flux, except 15%
      ! for grasses and crops
      EM_FRAC_AAXX = (/0.06d0, 0.06d0, 0.06d0, 0.06d0, 0.06d0,
     &                 0.06d0, 0.06d0, 0.06d0, 0.06d0, 0.06d0,
     &                 0.06d0, 0.15d0, 0.15d0, 0.15d0, 0.15d0/)

      ! Formaldehyde: 8% of bidirectional category flux, except 20%
      ! for grasses and crops
      EM_FRAC_CH2O = (/0.08d0, 0.08d0, 0.08d0, 0.08d0, 0.08d0,
     &                 0.08d0, 0.08d0, 0.08d0, 0.08d0, 0.08d0,
     &                 0.08d0, 0.20d0, 0.20d0, 0.20d0, 0.20d0/)

      ! Loop through plant types
      DO P = 1, 15

         ! Add 1 to Array_16 index to skip bare ground
         ARR_IND = P + 1
         
! ---> Now compute EFs for a-pinene and myrcene as well (dbm, 12/2012)
         ! a-pinene: 100% of category
         AEF_APIN(:,:) = AEF_APIN(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * PFT_EF_APIN(P)

         ! a-pinene: 100% of category
         AEF_MYRC(:,:) = AEF_MYRC(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * PFT_EF_MYRC(P)
! <---

         ! Other monoterpenes: 100% of category
         AEF_OMON(:,:) = AEF_OMON(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * PFT_EF_OMON(P)

         ! Methanol: 100% of category
         AEF_MOHX(:,:) = AEF_MOHX(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * PFT_EF_MOHX(P)

         ! Acetone: 100% of category
         AEF_ACET(:,:) = AEF_ACET(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * PFT_EF_ACET(P)

         ! Ethanol: variable fraction of category
         AEF_ETOH(:,:) = AEF_ETOH(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * 
     &        PFT_EF_BIDR(P) * EM_FRAC_ETOH(P)

         ! Formaldehyde: variable fraction of category
         AEF_CH2O(:,:) = AEF_CH2O(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * 
     &        PFT_EF_BIDR(P) * EM_FRAC_CH2O(P)

         ! Acetaldehyde: variable fraction of category
         AEF_ALD2(:,:) = AEF_ALD2(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * 
     &        PFT_EF_BIDR(P) * EM_FRAC_ALD2(P)

         ! Formic acid: variable fraction of category
         AEF_FAXX(:,:) = AEF_FAXX(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * 
     &        PFT_EF_BIDR(P) * EM_FRAC_FAXX(P)

         ! Acetic acid: variable fraction of category
         AEF_AAXX(:,:) = AEF_AAXX(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * 
     &        PFT_EF_BIDR(P) * EM_FRAC_AAXX(P)

         ! Ethene: 58% of "stress" category for all PFTs
         AEF_C2H4(:,:) = AEF_C2H4(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * PFT_EF_STRS(P) * 0.58d0

         ! Toluene: 3% of "stress" category for all PFTs
         AEF_TOLU(:,:) = AEF_TOLU(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * PFT_EF_STRS(P) * 0.03d0

         ! HCN: 1.5% of "stress" category for all PFTs
         AEF_HCNX(:,:) = AEF_HCNX(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * PFT_EF_STRS(P) * 0.015d0

         ! Propene: 48% of "other" category for all PFTs
         ! Butene:  24% of "other" category for all PFTs
         ! Larger alkenes: 0.2% of "other" category for all PFTs
         ! Total: 72.2% 
         AEF_PRPE(:,:) = AEF_PRPE(:,:) + 
     &        ARRAY_16(:,:,ARR_IND)/100d0 * PFT_EF_OTHR(P) * 0.722d0
         
      ENDDO

      ! Loop over longitudes
      DO J = 1, JJPAR
       
         ! Grid box surface area [m2]
         AREA_M2 = GET_AREA_M2( I, J, 1 )

         ! Conversion factor from [ug compound/m2/hr] to [g compound/box]
         FACTOR = 1.d-6 / 60.d0 * AREA_M2 * DTSRCE 
          
         ! Loop over latitudes
         DO I = 1, IIPAR

            ! Convert AEF arrays To [atomsC / box]
            AEF_ISOP(I,J)= AEF_ISOP(I,J) * FACTOR* AVO * 5d0 /68.117d0
            AEF_MBOX(I,J)= AEF_MBOX(I,J) * FACTOR* AVO * 5d0 /86.132d0
            AEF_APIN(I,J)= AEF_APIN(I,J) * FACTOR* AVO * 10d0 /136.234d0
            AEF_BPIN(I,J)= AEF_BPIN(I,J) * FACTOR* AVO * 10d0 /136.234d0
            AEF_CARE(I,J)= AEF_CARE(I,J) * FACTOR* AVO * 10d0 /136.234d0
            AEF_LIMO(I,J)= AEF_LIMO(I,J) * FACTOR* AVO * 10d0 /136.234d0 
            AEF_MYRC(I,J)= AEF_MYRC(I,J) * FACTOR* AVO * 10d0 /136.234d0
            AEF_OCIM(I,J)= AEF_OCIM(I,J) * FACTOR* AVO * 10d0 /136.234d0 
            AEF_SABI(I,J)= AEF_SABI(I,J) * FACTOR* AVO * 10d0 /136.234d0 
            AEF_OMON(I,J)= AEF_OMON(I,J) * FACTOR* AVO * 10d0 /136.234d0
            AEF_ACET(I,J)= AEF_ACET(I,J) * FACTOR* AVO * 3d0 /58.079d0 
            AEF_ETOH(I,J)= AEF_ETOH(I,J) * FACTOR* AVO * 2d0 /46.068d0 
            AEF_ALD2(I,J)= AEF_ALD2(I,J) * FACTOR* AVO * 2d0 /44.053d0 
            AEF_C2H4(I,J)= AEF_C2H4(I,J) * FACTOR* AVO * 2d0 /28.053d0  
            AEF_TOLU(I,J)= AEF_TOLU(I,J) * FACTOR* AVO * 7d0 /92.138d0 
            AEF_PRPE(I,J)= AEF_PRPE(I,J) * FACTOR* AVO * 3d0 /42.080d0  

            ! Methanol, formaldehyde, formic acid, acetic acid, HCN are
            ! carried in molec, not atomsC
            AEF_MOHX(I,J) = AEF_MOHX(I,J) * FACTOR * AVO / 32.042d0
            AEF_CH2O(I,J) = AEF_CH2O(I,J) * FACTOR * AVO / 30.026d0 
            AEF_FAXX(I,J) = AEF_FAXX(I,J) * FACTOR * AVO / 46.025d0
            AEF_AAXX(I,J) = AEF_AAXX(I,J) * FACTOR * AVO / 60.052d0
            AEF_HCNX(I,J) = AEF_HCNX(I,J) * FACTOR * AVO / 27.025d0

         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CALC_AEF
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_t_day
!
! !DESCRIPTION: Subroutine UPDATE\_T\_DAY must be called every time the 
!  A-3 fields are updated. Each 3h TS value for each gridbox is moved up one 
!  spot in the matrix and the current value is put in the last spot.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE UPDATE_T_DAY
!
! !USES:
!
      USE MEGANUT_MOD     ! We use all functions from the module

! 
! !REVISION HISTORY: 
!  (1 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (2 ) In GEOS4, TS are originally T2M in the A3 files, read in 
!        'a3_read_mod.f'. 
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Subroutine UPDATE\_T\_15\_AVG should be called at the 
!  beginningof each day. It loops through the gridboxes doing the following:
!
!  \begin{enumerate}
!  \item Average T\_DAY over the 8 TS values during the day.  
!  \item Push the daily average TS values through T_15, throwing out the 
!        oldest and putting the newest (the T\_DAY average) in the last spot 
!  \item Get T\_15\_AVG by averaging T\_15 over the 15 day period. 
!  \end{enumerate}
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE UPDATE_T_15_AVG
!
! !USES:
!
      IMPLICIT NONE
! 
! !REVISION HISTORY: 
!  01 Oct 1995 - M. Prather  - Initial version
!  (1 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (2 ) In GEOS4, TS are originally T2M in the A3 files, read in 
!        'a3_read_mod.f'. 
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_megan
!
! !DESCRIPTION: Subroutine INIT\_MEGAN allocates and initializes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_MEGAN
!
! !USES:
!
      USE A3_READ_MOD
      USE GEOS57_READ_MOD,  ONLY : GEOS57_READ_A1
      USE inquireMod,       ONLY : findFreeLUN
      USE MERRA_A1_MOD
      USE JULDAY_MOD,  ONLY : CALDATE
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LUNZIP
      USE TIME_MOD,    ONLY : GET_FIRST_A3_TIME, GET_JD
      USE TIME_MOD,    ONLY : ITS_A_LEAPYEAR,    YMD_EXTRACT
      
! 
! !REVISION HISTORY: 
!  (1 ) Change the logic in the #if block for G4AHEAD. (bmy, 12/6/05)
!  (2 ) Bug fix: skip Feb 29th if GCAP (phs, 9/18/07)
!  (3 ) Now call GET_AEF_05x0666 for GEOS-5 nested grids (yxw,dan,bmy, 11/6/08)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  26 Aug 2010 - R. Yantosca - Now reference merra_a1_mod.f
!  01 Sep 2010 - R. Yantosca - Now read in NUM_DAYS of sfc temp data (this had
!                              been hardwired to 15 days previously)
!  07 Feb 2011 - R. Yantosca - Fix typos: make sure to zero out the proper 
!                              PARDF_* and PARDR_* arrays after allocation
!  22 Nov 2011 - R. Yantosca - Do not use erroneous AEF's for nested grids
!  08 Feb 2012 - R. Yantosca - Now read surface temperature for GEOS-5.7.x
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  11 Apr 2012 - R. Yantosca - Now remove the call to INIT_LAI; we shall 
!                              initialize the LAI arrays from main.F
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL              :: GCAP_LEAP
      INTEGER              :: AS,        IUNIT
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

      ALLOCATE( AEF_MBOX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_MBOX' )
      AEF_MBOX = 0d0

      ALLOCATE( AEF_APIN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_APIN' )
      AEF_APIN = 0d0

      ALLOCATE( AEF_BPIN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_BPIN' )
      AEF_BPIN = 0d0

      ALLOCATE( AEF_CARE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_CARE' )
      AEF_CARE = 0d0

      ALLOCATE( AEF_LIMO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_LIMO' )
      AEF_LIMO = 0d0

      ALLOCATE( AEF_MYRC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_MYRC' )
      AEF_MYRC = 0d0

      ALLOCATE( AEF_OCIM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_OCIM' )
      AEF_OCIM = 0d0

      ALLOCATE( AEF_SABI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_SABI' )
      AEF_SABI = 0d0

      ALLOCATE( AEF_OMON( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_OMON' )
      AEF_OMON = 0d0

      ALLOCATE( AEF_MOHX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_MOHX' )
      AEF_MOHX = 0d0

      ALLOCATE( AEF_ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_ACET' )
      AEF_ACET = 0d0

      ALLOCATE( AEF_ETOH( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_ETOH' )
      AEF_ETOH = 0d0

      ALLOCATE( AEF_CH2O( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_CH2O' )
      AEF_CH2O = 0d0

      ALLOCATE( AEF_ALD2( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_ALD2' )
      AEF_ALD2 = 0d0

      ALLOCATE( AEF_FAXX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_FAXX' )
      AEF_FAXX = 0d0

      ALLOCATE( AEF_AAXX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_AAXX' )
      AEF_AAXX = 0d0

      ALLOCATE( AEF_C2H4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_C2H4' )
      AEF_C2H4 = 0d0

      ALLOCATE( AEF_TOLU( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_TOLU' )
      AEF_TOLU = 0d0

      ALLOCATE( AEF_HCNX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_HCNX' )
      AEF_HCNX = 0d0

      ALLOCATE( AEF_PRPE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEF_PRPE' )
      AEF_PRPE = 0d0

      ALLOCATE( NORM_FAC( 1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NORM_FAC' )
      NORM_FAC = -99d0

      ! Allocate arrays for light (mpb,2009)

      ! -- Direct --
      ALLOCATE( PARDR_DAY( IIPAR, JJPAR, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDR_DAY' )
      PARDR_DAY = 0d0

      ALLOCATE( PARDR_15( IIPAR, JJPAR, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDR_15' )
      PARDR_15 = 0d0

      ALLOCATE( PARDR_15_AVG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDR_15_AVG' )
      PARDR_15_AVG = 0d0

      ALLOCATE( PARDR_DAILY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDR_DAILY' )
      PARDR_DAILY = 0d0

      ! -- Diffuse --
      ALLOCATE( PARDF_DAY( IIPAR, JJPAR, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDF_DAY' )
      PARDF_DAY = 0d0

      ALLOCATE( PARDF_15( IIPAR, JJPAR, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDF_15' )
      PARDF_15 = 0d0

      ALLOCATE( PARDF_15_AVG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDF_15_AVG' )
      PARDF_15_AVG = 0d0

      ALLOCATE( PARDF_DAILY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDF_DAILY' )
      PARDF_DAILY = 0d0


      ! Generate annual emission factors for MEGAN inventory
      CALL CALC_AEF             

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

#if   defined( MERRA ) || defined( GEOS_57 )

      ! For MERRA and GEOS-5.7.x, the A1 fields are timestamped
      ! on the half-hours: 00:30, 01:30, 02:30, ... 23:30
      G4AHEAD   = 003000

#else

      ! For GEOS4, the A-3 fields are timestamped by the center of 
      ! the 3-hr period: 01:30Z, 04:30Z, 07:30Z, 10:30Z, 
      ! 13:30Z, 16:30Z, 19:30Z, 22:30Z
      G4AHEAD   = 013000

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
      DO I = NUM_DAYS+BACK_ONE, 1, -1

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
         DO J = 0, DAY_DIM-1  

#if   defined( GEOS_57 )

            !-----------------------------
            ! GEOS-5.7.x met fields
            ! Sfc temp is hourly data
            !-----------------------------
            CALL GEOS57_READ_A1( NYMD_T15, 010000*J + G4AHEAD ) 

#elif   defined( MERRA )

            !-----------------------------
            ! MERRA met fields
            ! Sfc temp is hourly data
            !-----------------------------

            ! Open A1 fields
            CALL OPEN_MERRA_A1_FIELDS( NYMD_T15, 
     &                                 010000*J + G4AHEAD, 
     &                                 IUNIT=IUNIT                  )

            ! Read A1 fields from disk
            CALL GET_MERRA_A1_FIELDS(  NYMD_T15, 010000*J + G4AHEAD ) 

#else
            !-----------------------------
            ! All other met field types!
            ! Sfc temp is 3-hourly data
            !-----------------------------

            ! Open A3 fields
            CALL OPEN_A3_FIELDS( NYMD_T15, 
     &                           030000*J + G4AHEAD,  
     &                           IUNIT=IUNIT                  )

            ! Read A3 fields from disk
            CALL GET_A3_FIELDS(  NYMD_T15, 030000*J + G4AHEAD ) 

#endif

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

#if   !defined( GEOS_57 ) 
      ! Close any open files
      CLOSE( IUNIT )
#endif

      ! Remove any leftover A-3 files in temp dir
      IF ( LUNZIP ) THEN
         CALL UNZIP_A3_FIELDS(  'remove all' )
      ENDIF

      ! Calculate normalization factor (dbm, 11/2012)
      CALL CALC_NORM_FAC

      ! Return to calling program
      END SUBROUTINE INIT_MEGAN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Subroutine CLEANUP\_MEGAN deallocates all allocated arrays 
!  at the end of a GEOS-Chem model run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_MEGAN
! 
! !REVISION HISTORY: 
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_MEGAN begins here!
      !=================================================================

      IF ( ALLOCATED( T_DAILY       ) ) DEALLOCATE( T_DAILY       )
      IF ( ALLOCATED( T_DAY         ) ) DEALLOCATE( T_DAY         )
      IF ( ALLOCATED( T_15          ) ) DEALLOCATE( T_15          )
      IF ( ALLOCATED( T_15_AVG      ) ) DEALLOCATE( T_15_AVG      )

      IF ( ALLOCATED( PARDR_DAILY   ) ) DEALLOCATE( PARDR_DAILY   )
      IF ( ALLOCATED( PARDR_DAY     ) ) DEALLOCATE( PARDR_DAY     )
      IF ( ALLOCATED( PARDR_15      ) ) DEALLOCATE( PARDR_15      )
      IF ( ALLOCATED( PARDR_15_AVG  ) ) DEALLOCATE( PARDR_15_AVG  )

      IF ( ALLOCATED( PARDF_DAILY   ) ) DEALLOCATE( PARDF_DAILY   )
      IF ( ALLOCATED( PARDF_DAY     ) ) DEALLOCATE( PARDF_DAY     )
      IF ( ALLOCATED( PARDF_15      ) ) DEALLOCATE( PARDF_15      )
      IF ( ALLOCATED( PARDF_15_AVG  ) ) DEALLOCATE( PARDF_15_AVG  )

      IF ( ALLOCATED( AEF_ISOP  ) ) DEALLOCATE( AEF_ISOP  )
      IF ( ALLOCATED( AEF_MBOX  ) ) DEALLOCATE( AEF_MBOX  )
      IF ( ALLOCATED( AEF_APIN  ) ) DEALLOCATE( AEF_APIN  )
      IF ( ALLOCATED( AEF_BPIN  ) ) DEALLOCATE( AEF_BPIN  )
      IF ( ALLOCATED( AEF_CARE  ) ) DEALLOCATE( AEF_CARE  )
      IF ( ALLOCATED( AEF_LIMO  ) ) DEALLOCATE( AEF_LIMO  )
      IF ( ALLOCATED( AEF_MYRC  ) ) DEALLOCATE( AEF_MYRC  )
      IF ( ALLOCATED( AEF_OCIM  ) ) DEALLOCATE( AEF_OCIM  )
      IF ( ALLOCATED( AEF_SABI  ) ) DEALLOCATE( AEF_SABI  )
      IF ( ALLOCATED( AEF_OMON  ) ) DEALLOCATE( AEF_OMON  )
      IF ( ALLOCATED( AEF_MOHX  ) ) DEALLOCATE( AEF_MOHX  )
      IF ( ALLOCATED( AEF_ACET  ) ) DEALLOCATE( AEF_ACET  )
      IF ( ALLOCATED( AEF_ETOH  ) ) DEALLOCATE( AEF_ETOH  )
      IF ( ALLOCATED( AEF_CH2O  ) ) DEALLOCATE( AEF_CH2O  )
      IF ( ALLOCATED( AEF_ALD2  ) ) DEALLOCATE( AEF_ALD2  )
      IF ( ALLOCATED( AEF_FAXX  ) ) DEALLOCATE( AEF_FAXX  )
      IF ( ALLOCATED( AEF_AAXX  ) ) DEALLOCATE( AEF_AAXX  )
      IF ( ALLOCATED( AEF_C2H4  ) ) DEALLOCATE( AEF_C2H4  )
      IF ( ALLOCATED( AEF_TOLU  ) ) DEALLOCATE( AEF_TOLU  )
      IF ( ALLOCATED( AEF_HCNX  ) ) DEALLOCATE( AEF_HCNX  )
      IF ( ALLOCATED( AEF_PRPE  ) ) DEALLOCATE( AEF_PRPE  )
      IF ( ALLOCATED( NORM_FAC  ) ) DEALLOCATE( NORM_FAC  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_MEGAN
!EOC
      END MODULE MEGAN_MOD
