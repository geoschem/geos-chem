! $Id: isoropia.h,v 1.1 2010/02/25 21:07:02 bmy Exp $
!
!******************************************************************************
!  Header file "isoropia.h" contains common block declarations for the 
!  ISORROPIA code.  We had to keep these common blocks so that we could
!  declare these as THREADPRIVATE for the OpenMP parallel loop in
!  routine AERO_THERMO in "isoropia_mod.f". (bec, bmy, 3/7/05, 6/28/06) 
!
!  Assume all common blocks are THREADPRIVATE unless otherwise specified.
!  All parameters are specified in MODULE section of "isoropia_mod.f".
!
!  Original Documentation:
!  *** ISORROPIA CODE
!  *** INCLUDE FILE 'ISRPIA.INC'
!  *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!      AND VARIABLES. 
!
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  NOTES:
!  (1 ) Common blocks were split up so that each common block only holds one
!        type of data (INTEGER, LOGICAL, REAL*8).  This avoids problems with
!        mixing data types.  Also now explicitly declare all variables and
!        remove old F77-style IMPLICIT typing.  Also removed some variables
!        which are constants and declared them as PARAMETERS w/in the
!        "isoropia_mod.f". (bmy, 3/7/05)
!  (2 ) Rename the /GAS/ common block to /GAS1/, because there are variables
!        named GAS within "isoropia_mod.f".  This avoids namespace confusion, 
!        which may cause some compilers to choke. (psk, bmy, 6/28/06)
!******************************************************************************
!
      !=================================================================
      ! Input variables
      !=================================================================

      REAL*8             :: W(NCOMP), WAER(NCOMP), TEMP, RHB
      COMMON /INPT2/        W,        WAER,        TEMP, RHB

      !=================================================================
      ! Water activities of pure salt solutions
      !=================================================================
      
      ! /ZSR/ is read-only and doesn't have to be declared THREADPRIVATE
      REAL*8             :: AWAS(NZSR), AWSS(NZSR), AWAC(NZSR)
      REAL*8             :: AWSC(NZSR), AWAN(NZSR), AWSN(NZSR)
      REAL*8             :: AWSB(NZSR), AWAB(NZSR), AWSA(NZSR)
      REAL*8             :: AWLC(NZSR)
      COMMON /ZSR/          AWAS,       AWSS,       AWAC,       
     &                      AWSC,       AWAN,       AWSN,       
     &                      AWSB,       AWAB,       AWSA,       
     &                      AWLC

      !=================================================================
      ! Deliquescence relative humidities
      !=================================================================

      REAL*8             :: DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL   
      REAL*8             :: DRNANO3,  DRNA2SO4, DRNH4HS4, DRLC     
      REAL*8             :: DRNH4NO3, DRNH4CL
      COMMON /DRH1/         DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL, 
     &                      DRNANO3,  DRNA2SO4, DRNH4HS4, DRLC,     
     &                      DRNH4NO3, DRNH4CL

      REAL*8             :: DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1   
      REAL*8             :: DRMG2,    DRMG3,    DRMH1,    DRMH2    
      REAL*8             :: DRMI1,    DRMI2,    DRMI3,    DRMQ1    
      REAL*8             :: DRMR1,    DRMR2,    DRMR3,    DRMR4
      REAL*8             :: DRMR5,    DRMR6,    DRMR7,    DRMR8
      REAL*8             :: DRMR9,    DRMR10,   DRMR11,   DRMR12  
      REAL*8             :: DRMR13
      COMMON /DRH2/         DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,   
     &                      DRMG2,    DRMG3,    DRMH1,    DRMH2,    
     &                      DRMI1,    DRMI2,    DRMI3,    DRMQ1,    
     &                      DRMR1,    DRMR2,    DRMR3,    DRMR4,    
     &                      DRMR5,    DRMR6,    DRMR7,    DRMR8,
     &                      DRMR9,    DRMR10,   DRMR11,   DRMR12,  
     &                      DRMR13

      !=================================================================
      ! Variables for liquid aerosol phase
      !=================================================================

      LOGICAL            :: CALAIN, CALAOU, DRYF, FRST   
      COMMON /IONS0/        CALAIN, CALAOU, DRYF, FRST

      ! /IONS1/ is readonly and doesn't have to be declared THREADPRIVATE
      INTEGER            :: ZZ(NPAIR),    Z(NIONS)
      COMMON /IONS1/        ZZ,           Z

      REAL*8             :: MOLAL(NIONS), MOLALR(NPAIR), GAMA(NPAIR) 
      REAL*8             :: GAMOU(NPAIR), GAMIN(NPAIR),  M0(NPAIR)
      REAL*8             :: GASAQ(NGASAQ)
      COMMON /IONS2/        MOLAL,        MOLALR,        GAMA, 
     &                      GAMOU,        GAMIN,         M0,
     &                      GASAQ

      REAL*8             :: COH,          CHNO3,         CHCL
      REAL*8             :: WATER,        IONIC   
      COMMON /IONS3/        COH,          CHNO3,         CHCL,         
     &                      WATER,        IONIC         

      !=================================================================
      ! Variables for solid aerosol phase
      !=================================================================

      REAL*8             :: CH2SO4, CNH42S4, CNH4HS4, CNaCL,   CNa2SO4
      REAL*8             :: CNaNO3, CNH4NO3, CNH4CL,  CNaHSO4, CLC
      COMMON /SALT/         CH2SO4, CNH42S4, CNH4HS4, CNACL,   CNA2SO4, 
     &                      CNANO3, CNH4NO3, CNH4CL,  CNAHSO4, CLC

      !=================================================================
      ! Variables for gas phase
      !=================================================================

      REAL*8             :: GNH3, GHNO3, GHCL 
      COMMON /GAS1/         GNH3, GHNO3, GHCL 

      !=================================================================
      ! Equilibrium constants
      !=================================================================

      REAL*8             :: XK1,  XK2,  XK3,  XK4,  XK5,  XK6,  XK7
      REAL*8             :: XK8,  XK9,  XK10, XK11, XK12, XK13, XK14 
      REAL*8             :: XKW,  XK21, XK22, XK31, XK32, XK41, XK42
      COMMON /EQUK/         XK1,  XK2,  XK3,  XK4,  XK5,  XK6,  XK7, 
     &                      XK8,  XK9,  XK10, XK11, XK12, XK13, XK14,  
     &                      XKW,  XK21, XK22, XK31, XK32, XK41, XK42

      !=================================================================
      ! Solution/info variables
      !=================================================================

      CHARACTER(LEN=15)  :: SCASE
      COMMON /CASE0/        SCASE

      REAL*8             :: SULRATW, SULRAT, SODRAT
      COMMON /CASE1/        SULRATW, SULRAT, SODRAT

      INTEGER            :: ICLACT
      COMMON /SOLN2/        ICLACT

      !=================================================================
      ! Error system
      !=================================================================

      CHARACTER(LEN=40)  :: ERRMSG
      COMMON /EROR0/        ERRMSG(NERRMX)

      INTEGER            :: ERRSTK,         NOFER   
      COMMON /EROR1/        ERRSTK(NERRMX), NOFER  

      LOGICAL            :: STKOFL   
      COMMON /EROR2/        STKOFL

      !=================================================================
      ! Coordinates for debugging
      !=================================================================

      INTEGER            :: I_SAV, J_SAV, L_SAV 
      COMMON /COORD/        I_SAV, J_SAV, L_SAV

      !=================================================================
      ! THREADPRIVATE declarations for OpenMP parallelization
      !=================================================================
!$OMP THREADPRIVATE( /CASE0/ )
!$OMP THREADPRIVATE( /CASE1/ )
!$OMP THREADPRIVATE( /COORD/ )
!$OMP THREADPRIVATE( /DRH1/  )
!$OMP THREADPRIVATE( /DRH2/  )
!$OMP THREADPRIVATE( /EQUK/  )
!$OMP THREADPRIVATE( /EROR0/ )
!$OMP THREADPRIVATE( /EROR1/ )
!$OMP THREADPRIVATE( /EROR2/ )
!$OMP THREADPRIVATE( /GAS1/  )
!$OMP THREADPRIVATE( /INPT2/ )
!$OMP THREADPRIVATE( /IONS0/ )
!$OMP THREADPRIVATE( /IONS2/ )
!$OMP THREADPRIVATE( /IONS3/ )
!$OMP THREADPRIVATE( /SALT/  )
!$OMP THREADPRIVATE( /SOLN2/ )

