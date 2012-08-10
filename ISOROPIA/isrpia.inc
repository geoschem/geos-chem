!******************************************************************************
!  File "isrpia.inc" contains the common block declarations and some 
!  parameters for the ISORROPIAII code. This file replaces "isoropia.h" from
!  previous ISORROPIA implementations.
!  ( hotp 8/1/07 )
!  ( previous implementation: bec, bmy, 3/7/05, 6/28/06 )
!
!  Please contact the original authors before making modifications to
!  ISOROPIAII
!
!  Original Documentation:
!  *** ISORROPIA CODE II
!  *** INCLUDE FILE 'ISRPIA.INC'
!  *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!      AND VARIABLES. 
!
!  *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
!  *** GEORGIA INSTITUTE OF TECHNOLOGY
!  *** WRITTEN BY ATHANASIOS NENES
!  *** UPDATED BY CHRISTOS FOUNTOUKIS
!
!  NOTES: 
!  (1 ) Renamed common block GAS to GASCB to avoid confusion with array
!        GAS in isoropiaIIcode.f (hotp)
!  (2 ) Explicity declared all variables in this file (hotp)
!  (3 ) Separated common blocks so that each block has only one type of data (hotp)
!  (4 ) Declared all common blocks (except READ-ONLY blocks) as
!        THREADPRIVATE (hotp)
!  (5 ) Divided common blocks so that each has no more than 8 variables 
!        (this is to speed compilation) (hotp)
!  (6 ) Divided common blocks so that READ-ONLY variables are not mixed with
!        variables that are modified (hotp)
!  (7 ) Removed DRYF since it is not used (hotp)
!  (8 ) Bug fix: CRRAT was being truncated at column 72 (hotp, bmy, 6/1/10)
!  (9 ) IONIC should be REAL*8 for consistency with -r8 flag during
!        compiling of isoropiaIIcode routines (hotp, bmy, 6/28/10)
!******************************************************************************
!
      ! leave this implicit statement for isoropiaIIcode.f variables
      ! ideally all variables would be declared
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      ! INTEGER parameters
      INTEGER, PARAMETER :: NCOMP  =   8
      INTEGER, PARAMETER :: NIONS  =  10
      INTEGER, PARAMETER :: NGASAQ =   3
      INTEGER, PARAMETER :: NSLDS  =  19
      INTEGER, PARAMETER :: NPAIR  =  23
      INTEGER, PARAMETER :: NZSR   = 100
      INTEGER, PARAMETER :: NERRMX =  25

      !=================================================================
      ! Input variables
      !=================================================================

      INTEGER            :: IPROB, METSTBL
      COMMON /INPT/         IPROB, METSTBL
!$OMP THREADPRIVATE( /INPT/ )
       
      INTEGER            :: NADJ
      COMMON /INPT2/        NADJ
      ! INPT2 is read-only

      REAL*8             :: W(NCOMP), WAER(NCOMP), TEMP, RH
      COMMON /INPT3/        W,        WAER,        TEMP, RH
!$OMP THREADPRIVATE( /INPT3/ )

      !=================================================================
      ! Water activities of pure salt solutions
      !=================================================================

      ! /ZSR/ is read-only and doesn't have to be declared THREADPRIVATE
      ! block does not need to be split up to speed compilation
      REAL*8     :: AWAS(NZSR), AWSS(NZSR), AWAC(NZSR), AWSC(NZSR)
      REAL*8     :: AWAN(NZSR), AWSN(NZSR), AWSB(NZSR), AWAB(NZSR)
      REAL*8     :: AWSA(NZSR), AWLC(NZSR), AWCS(NZSR), AWCN(NZSR)
      REAL*8     :: AWCC(NZSR), AWPS(NZSR), AWPB(NZSR), AWPN(NZSR)
      REAL*8     :: AWPC(NZSR), AWMS(NZSR), AWMN(NZSR), AWMC(NZSR)
      COMMON /ZSR / AWAS,       AWSS,       AWAC,       AWSC,
     &              AWAN,       AWSN,       AWSB,       AWAB,
     &              AWSA,       AWLC,       AWCS,       AWCN,
     &              AWCC,       AWPS,       AWPB,       AWPN,
     &              AWPC,       AWMS,       AWMN,       AWMC
      ! ZSR is read-only

      !=================================================================
      ! Deliquescence relative humidities
      !=================================================================

      INTEGER    :: WFTYP
      COMMON /WFT/  WFTYP
      ! WFTYP is read-only

      REAL*8     :: DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3 
      REAL*8     :: DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
      REAL*8     :: DRCASO4,  DRCANO32, DRCACL2,  DRK2SO4,  DRKHSO4
      REAL*8     :: DRKNO3,   DRKCL,    DRMGSO4,  DRMGNO32, DRMGCL2
      COMMON /DRH / DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3 
      COMMON /DRH2/ DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
      COMMON /DRH3/ DRCASO4,  DRCANO32, DRCACL2,  DRK2SO4,  DRKHSO4
      COMMON /DRH4/ DRKNO3,   DRKCL,    DRMGSO4,  DRMGNO32, DRMGCL2
!$OMP THREADPRIVATE( /DRH/ )
!$OMP THREADPRIVATE( /DRH2/ )
!$OMP THREADPRIVATE( /DRH3/ )
!$OMP THREADPRIVATE( /DRH4/ )

      REAL*8     ::  DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2
      REAL*8     ::  DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2
      REAL*8     ::  DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3
      REAL*8     ::  DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8
      REAL*8     ::  DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13
      COMMON /MDRH/  DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2
      COMMON /MDRH2/ DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2
      COMMON /MDRH3/ DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3
      COMMON /MDRH4/ DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8
      COMMON /MDRH5/ DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13
!$OMP THREADPRIVATE( /MDRH/ )
!$OMP THREADPRIVATE( /MDRH2/ )
!$OMP THREADPRIVATE( /MDRH3/ )
!$OMP THREADPRIVATE( /MDRH4/ )
!$OMP THREADPRIVATE( /MDRH5/ )

      REAL*8     ::  DRMO1,    DRMO2,    DRMO3,    DRML1,    DRML2
      REAL*8     ::  DRML3,    DRMM1,    DRMM2,    DRMP1,    DRMP2
      REAL*8     ::  DRMP3,    DRMP4,    DRMP5,    DRMV1
      COMMON /MDRH6/ DRMO1,    DRMO2,    DRMO3,    DRML1,    DRML2
      COMMON /MDRH7/ DRML3,    DRMM1,    DRMM2,    DRMP1,    DRMP2
      COMMON /MDRH8/ DRMP3,    DRMP4,    DRMP5,    DRMV1
!$OMP THREADPRIVATE( /MDRH6/ )
!$OMP THREADPRIVATE( /MDRH7/ )
!$OMP THREADPRIVATE( /MDRH8/ )

      !=================================================================
      ! Variables for liquid aerosol phase
      !=================================================================

      ! /IONS/ size does not seem to slow compilation
      REAL*8     :: MOLAL(NIONS), MOLALR(NPAIR), M0(NPAIR)
      REAL*8     :: GAMA(NPAIR)  
      REAL*8     :: GAMOU(NPAIR), GAMIN(NPAIR),  GASAQ(NGASAQ)
      REAL*8     :: COH,          CHNO3,         CHCL
      REAL*8     :: WATER         
      COMMON /IONS/ MOLAL,        MOLALR,        M0,
     &              GAMA,         
     &              GAMOU,        GAMIN,         GASAQ,
     &              COH,          CHNO3,         CHCL,
     &              WATER    
!$OMP THREADPRIVATE( /IONS/ )

      REAL*8     ::  ZZ(NPAIR),   Z(NIONS),      EPSACT
      COMMON /IONS2/ ZZ,          Z,             EPSACT
      ! IONS2 is read-only
      
      !----------------------------------------------------------------------
      ! Prior to 6/28/10:
      ! IONIC should be REAL*8 for consistency with -r8 flag during
      ! compiling of isoropiaIIcode routines (hotp 6/23/10)
#if   defined( ESMF_ )
#if   defined( DEVEL )
      REAL*4     ::  IONIC
#endif
#else
      REAL*8     ::  IONIC
#endif
      !----------------------------------------------------------------------
      COMMON /IONS3/ IONIC
!$OMP THREADPRIVATE( /IONS3/ )

      ! DRYF removed since it is not used (hotp)
      LOGICAL    ::  CALAOU, CALAIN, FRST
      COMMON /IONS4/ CALAOU, CALAIN, FRST
!$OMP THREADPRIVATE( /IONS4/ )

      INTEGER    ::  IACALC
      COMMON /IONS5/ IACALC
      ! IONS5 is read-only
      !=================================================================
      ! Variables for solid aerosol phase
      !=================================================================

      REAL*8     ::  CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4 
      REAL*8     ::  CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC, CCASO4
      REAL*8     ::  CCANO32, CCACL2,  CK2SO4,  CKHSO4,  CKNO3, CKCL
      REAL*8     ::  CMGSO4,  CMGNO32, CMGCL2
      COMMON /SALT/  CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4 
      COMMON /SALT2/ CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC, CCASO4
      COMMON /SALT3/ CCANO32, CCACL2,  CK2SO4,  CKHSO4,  CKNO3, CKCL
      COMMON /SALT4/ CMGSO4,  CMGNO32, CMGCL2
!$OMP THREADPRIVATE( /SALT/ )
!$OMP THREADPRIVATE( /SALT2/ )
!$OMP THREADPRIVATE( /SALT3/ )
!$OMP THREADPRIVATE( /SALT4/ )

      !=================================================================
      ! Variables for gas phase
      !=================================================================

      REAL*8     ::  GNH3, GHNO3, GHCL 
      COMMON /GASCB/ GNH3, GHNO3, GHCL 
!$OMP THREADPRIVATE( /GASCB/ )

      !=================================================================
      ! Equilibrium constants
      !=================================================================

      REAL*8     :: XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10
      REAL*8     :: XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41
      REAL*8     :: XK42,XK15,XK16,XK17,XK18,XK19,XK20,XK23
      REAL*8     :: XK24,XK25
      COMMON /EQUK/  XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10
      COMMON /EQUK2/ XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41
      COMMON /EQUK3/ XK42,XK15,XK16,XK17,XK18,XK19,XK20,XK23
      COMMON /EQUK4/ XK24,XK25
      !C     &              , XK26, XK27
!$OMP THREADPRIVATE( /EQUK/ )
!$OMP THREADPRIVATE( /EQUK2/ )
!$OMP THREADPRIVATE( /EQUK3/ )
!$OMP THREADPRIVATE( /EQUK4/ )

      !=================================================================
      ! Molecular Weights
      !=================================================================

      REAL*8     :: R, IMW(NIONS), WMW(NCOMP), SMW(NPAIR)
      COMMON /OTHR/ R, IMW,        WMW,        SMW
      ! OTHR is read-only

      !=================================================================
      ! Solution/info variables
      !=================================================================

      CHARACTER(LEN=15) :: SCASE
      COMMON /CASE/        SCASE
!$OMP THREADPRIVATE( /CASE/ )

      REAL*8            :: SULRATW, SULRAT,  SODRAT,
     &                     SO4RAT,  CRNARAT, CRRAT
      COMMON /CASE2/       SULRATW, SULRAT,  SODRAT, 
     &                     SO4RAT,  CRNARAT, CRRAT
!$OMP THREADPRIVATE( /CASE2/ )

      REAL*8            :: EPS
      COMMON /SOLN/        EPS
      ! SOLN is read-only

      INTEGER           :: MAXIT, NSWEEP, NDIV, ICLACT
      COMMON /SOLN2/       MAXIT, NSWEEP, NDIV
      ! SOLN 2 is read-only
      
      COMMON /SOLN3/       ICLACT
!$OMP THREADPRIVATE( /SOLN3/ )

      !=================================================================
      ! Error system
      !=================================================================

      CHARACTER(LEN=40) :: ERRMSG(NERRMX)
      COMMON /EROR/        ERRMSG
!$OMP THREADPRIVATE( /EROR/ )

      INTEGER           :: ERRSTK(NERRMX), NOFER   
      COMMON /EROR2/       ERRSTK,         NOFER   
!$OMP THREADPRIVATE( /EROR2/ )

      LOGICAL           :: STKOFL   
      COMMON /EROR3/       STKOFL   
!$OMP THREADPRIVATE( /EROR3/ )

      !=================================================================
      ! Generic Variables
      !=================================================================

      CHARACTER(LEN=15) :: VERSION
      COMMON /CGEN/        VERSION
      ! CGEN is read-only

      REAL*8            :: GREAT, TINY, TINY2, ZERO, ONE
      COMMON /CGEN2/       GREAT, TINY, TINY2, ZERO, ONE
      ! CGEN2 is read-only

      !=================================================================
      ! coordinates for debugging
      !=================================================================
      INTEGER :: ICOOR, JCOOR, LCOOR
      COMMON /LOC/ ICOOR, JCOOR, LCOOR
!$OMP THREADPRIVATE( /LOC/ )
C
C *** END OF INCLUDE FILE **********************************************
C
