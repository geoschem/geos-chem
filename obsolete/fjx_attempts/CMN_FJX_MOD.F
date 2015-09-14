!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cmn_fjx_mod
!
! !DESCRIPTION: Module CMN\_FJX\_MOD contains parameters and global variables
!  used to interface between Harvard chemistry and UC-Irvine photolysis 
!  programs (Fast-J/Fast-JX), along with all Fast-J(X) global variables
!  and some physical constants for the GEOS-Chem chemistry code.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_FJX_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD,      ONLY : IIPAR, JJPAR, LLPAR
      USE CMN_SIZE_MOD,      ONLY : NDUST, NAER, NRH
#if defined( DEVEL )
      USE SMV_DIMENSION_MOD, ONLY : MAX_COLUMN
#endif

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      ! New (SDE 03/28/13)
      ! Index in RAA & QAA of 999 nm wavelength
      INTEGER, PARAMETER :: IND999 = 4

      ! Adapted from CMN_FJ_MOD
      ! max # of photolysis rxns = 4 + IPHOT (see comode.h)
      ! FP increased JPMAX since IPHOT was increased (hotp 7/31/09)
      !INTEGER, PARAMETER :: JPMAX = 79
      INTEGER, PARAMETER :: JPMAX = 89

      ! Adapted from JV_MIE_MOD
      ! NL    Maximum number of levels after insertion of extra Mie levels
      ! N__   Number of levels in Mie grid: 2*(2*lpar+2+jaddto(1))+3
      ! M__   Number of Gauss points used
!-----------------------------------------------------------------------
!     NL=1500 was too small for dicarbonyls, so we upped it to 2000.
!     Uncomment this line to restore the previous definition (phs, 10/9/09)
!      INTEGER, PARAMETER :: NL  = 1500
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: NL  = 2000
      INTEGER, PARAMETER :: N__ = 2*NL
      INTEGER, PARAMETER :: M__ = 4 

      ! Adapted from JV_CMN (SDE 03/27/13)
#if defined( DEVEL )
      ! NB  Number of levels in CTM plus one for above model top
      ! NC  Number of levels in the fundamental Fast-J grid
      ! NS  Maximum number of species which require J-values calculating 
      ! NW  Maximum number of wavelength bins that can be used
      ! NP  Maximum number of aerosol/cloud types that can be used
      ! MX  Number of aerosol/cloud types supplied from CTM
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! NOTE: MAX_COLUMN is set to 47L for GRIDREDUCED or 72L otherwise
      ! This is kludge to let us test the DEVEL code (bmy, mlong, 8/10/12)
      INTEGER, PARAMETER :: NB      = MAX_COLUMN+1
#else
      INTEGER, PARAMETER :: NB      = LLPAR+1
#endif
      INTEGER, PARAMETER :: NC      = 2*NB
      INTEGER, PARAMETER :: NS      = 51
      INTEGER, PARAMETER :: NW      = 15
      INTEGER, PARAMETER :: NP      = 56
      INTEGER, PARAMETER :: MX      = 35

      REAL*8,  PARAMETER :: RAD     = 6375.d5   
      REAL*8,  PARAMETER :: ZZHT    = 5.d5      
      REAL*8,  PARAMETER :: dtaumax = 1.d0
      REAL*8,  PARAMETER :: dtausub = 1.d0
      REAL*8,  PARAMETER :: dsubdiv = 10.d0
      REAL*8,  PARAMETER :: szamax  = 98.0d0
!
! !PUBLIC DATA MEMBERS:
!
      ! Character variables
      CHARACTER*20 TITLEA(NP)
      CHARACTER*78 TITLE0
      CHARACTER*7  TITLEJ(3,NS), jlabel(JPMAX)

      !-----------------------------------------------------------------------
      ! These common blocks MUST NOT be held local (bmy, 5/2/00)

      REAL*8  :: WBIN(NW+1),WL(NW),FL(NW),QO2(NW,3),QO3(NW,3)       
      REAL*8  :: Q1D(NW,3),QQQ(NW,2,NS-3),QRAYL(NW),TQQ(3,NS)       
      REAL*8  :: WAA(4,NP),QAA(4,NP)                                 
      REAL*8  :: PAA(8,4,NP),RAA(4,NP),SSA(4,NP),QBC(NW)
            
      INTEGER :: NJVAL,NW1,NW2,NAA,NLBATM                                 

      REAL*8  :: WAA_AOD(NP),QAA_AOD(NP),PAA_AOD(8,NP)             
      REAL*8  :: RAA_AOD(NP),SSA_AOD(NP)

      REAL*8  :: TREF(51,18,12),OREF(51,18,12),BREF(51)           
#if defined( DEVEL )
      REAL*8, ALLOCATABLE  :: ODMDUST(:,:,:,:)
      REAL*8, ALLOCATABLE  :: ODAER(:,:,:,:)
#else
      REAL*8  :: ODMDUST(IIPAR,JJPAR,LLPAR,NDUST   )
      REAL*8  :: ODAER  (IIPAR,JJPAR,LLPAR,NAER*NRH)
#endif

      REAL*8  :: JFACTA(JPMAX),ZPDEP(NW,7)
      INTEGER :: NPDEP,JPDEP(NS),JIND(JPMAX)

      INTEGER :: MIEDX(MX)

      !-----------------------------------------------------------------------
      ! These common blocks MUST be held local for the parallelization 
      ! (bmy, 5/2/00)
      REAL*8   :: TJ(NB),PJ(NB+1),DM(NB),DO3(NB),DBC(NB),Z(NB)    
      REAL*8   :: AER(MX,NB),AMF(NB,NB),RFLECT,SZA,U0,TANHT
      REAL*8   :: zj(NB,JPMAX)
      REAL*8   :: FFF(NW,NB),VALJ(NS)
      INTEGER  :: jadsub(NC)

!$OMP THREADPRIVATE( TJ, PJ, DM, DO3, DBC, Z)           
!$OMP THREADPRIVATE( AER, AMF, RFLECT, SZA, U0, TANHT)
!$OMP THREADPRIVATE( zj )
!$OMP THREADPRIVATE( FFF, VALJ )
!$OMP THREADPRIVATE( jadsub )
#else
      ! NB  Number of levels in CTM plus one for above model top
      ! NC  Number of levels in the fundamental Fast-J grid
      ! NS  Maximum number of species which require J-values calculating 
      ! NW  Maximum number of wavelength bins that can be used
      ! NP  Maximum number of aerosol/cloud types that can be used
      ! MX  Number of aerosol/cloud types supplied from CTM
!----------------------------------------------------------------------------
! Prior to 8/10/12:
! Replace LPAR with LLPAR (bmy, 8/10/12)
!      INTEGER, PARAMETER :: NB      = LPAR+1
!----------------------------------------------------------------------------
      INTEGER, PARAMETER :: NB      = LLPAR+1
      INTEGER, PARAMETER :: NC      = 2*NB
      INTEGER, PARAMETER :: NS      = 51
      INTEGER, PARAMETER :: NW      = 15
      INTEGER, PARAMETER :: NP      = 56
      INTEGER, PARAMETER :: MX      = 35

      REAL*8,  PARAMETER :: RAD     = 6375.d5   
      REAL*8,  PARAMETER :: ZZHT    = 5.d5      
      REAL*8,  PARAMETER :: dtaumax = 1.d0
      REAL*8,  PARAMETER :: dtausub = 1.d0
      REAL*8,  PARAMETER :: dsubdiv = 10.d0
      REAL*8,  PARAMETER :: szamax  = 98.0d0
!
! !PUBLIC DATA MEMBERS:
!
      ! Character variables
      CHARACTER*20 TITLEA(NP)
      CHARACTER*78 TITLE0
      CHARACTER*7  TITLEJ(3,NS), jlabel(JPMAX)

      !-----------------------------------------------------------------------
      ! These common blocks MUST NOT be held local (bmy, 5/2/00)

      REAL*8  :: WBIN(NW+1),WL(NW),FL(NW),QO2(NW,3),QO3(NW,3)       
      REAL*8  :: Q1D(NW,3),QQQ(NW,2,NS-3),QRAYL(NW),TQQ(3,NS)       
      REAL*8  :: WAA(4,NP),QAA(4,NP)                                 
      REAL*8  :: PAA(8,4,NP),RAA(4,NP),SSA(4,NP),QBC(NW)
            
      INTEGER :: NJVAL,NW1,NW2,NAA,NLBATM                                 

      REAL*8  :: WAA_AOD(NP),QAA_AOD(NP),PAA_AOD(8,NP)             
      REAL*8  :: RAA_AOD(NP),SSA_AOD(NP)

      REAL*8  :: TREF(51,18,12),OREF(51,18,12),BREF(51)           
!-----------------------------------------------------------
! Prior to 8/10/12:
! Replace IPAR, JPAR, LPAR with IIPAr, JJPAR, LLPAR
!      REAL*8  :: ODMDUST(IPAR,JPAR,LPAR,NDUST)
!      REAL*8  :: ODAER(IPAR,JPAR,LPAR,NAER*NRH)
!-----------------------------------------------------------
      REAL*8  :: ODMDUST(IIPAR,JJPAR,LLPAR,NDUST)
      REAL*8  :: ODAER  (IIPAR,JJPAR,LLPAR,NAER*NRH)

      REAL*8  :: jfacta(JPMAX),zpdep(NW,7)
      INTEGER :: npdep,jpdep(NS),jind(JPMAX)

      INTEGER :: MIEDX(MX)

      !-----------------------------------------------------------------------
      ! Split off GLYX-chemistry specific arrays into separate common blocks
      ! (ccarouge, bmy, 8/20/09)
      INTEGER :: PDEPF(7)

      REAL*8  :: MGLYPDEP(NW, 3)

      !-----------------------------------------------------------------------
      ! These common blocks MUST be held local for the parallelization 
      ! (bmy, 5/2/00)
      REAL*8   :: TJ(NB),PJ(NB+1),DM(NB),DO3(NB),DBC(NB),Z(NB)    
      REAL*8   :: AER(MX,NB),AMF(NB,NB),RFLECT,SZA,U0,TANHT
!----------------------------------------------
! Prior to 8/10/12:
! Replace LPAR with LLPAR (bmy, 8/10/12)
!      REAL*8   :: zj(LPAR,JPMAX)
!      REAL*8   :: FFF(NW,lpar),VALJ(NS)
!----------------------------------------------
      REAL*8   :: zj(LLPAR,JPMAX)
      REAL*8   :: FFF(NW,LLPAR),VALJ(NS)
      INTEGER  :: jadsub(NC)

!$OMP THREADPRIVATE( TJ, PJ, DM, DO3, DBC, Z)           
!$OMP THREADPRIVATE( AER, AMF, RFLECT, SZA, U0, TANHT)
!$OMP THREADPRIVATE( zj )
!$OMP THREADPRIVATE( FFF, VALJ )
!$OMP THREADPRIVATE( jadsub )
#endif

      ! Adapted from CMN_FJ_MOD
      ! Variables for number of layers and number of photolysis rxns
      INTEGER            :: JPNL, JPPJ       

      ! Branches for photolysis species
      INTEGER            :: BRANCH(JPMAX)

      ! Names of photolysis species
      ! FP increased length of RNAMES for species indistinguishable
      ! with only 4 characters (hotp 7/31/09)
      ! used in jv_index and rd_js.f
      !CHARACTER (LEN=4)  :: RNAMES
      CHARACTER (LEN=7)  :: RNAMES(JPMAX)

      ! Mapping array from Harvard species names to UCI species names
      INTEGER            :: RINDEX(JPMAX)

      ! Output J values
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      REAL*8, ALLOCATABLE :: ZPJ(:,:,:,:)
#else
      REAL*8             :: ZPJ(LLPAR,JPMAX,IIPAR,JJPAR)
#endif

      ! Adapted from JV_MIE_MOD
      ! Arrays
      REAL*8  :: A(M__),           B(M__,M__),   C1(M__)
      REAL*8  :: H(M__),           AA(M__,M__),  CC(M__,M__)
      REAL*8  :: S(M__,M__),       W(M__,M__),   U1(M__,M__)
      REAL*8  :: V1(M__),          WT(M__),      EMU(M__)
      REAL*8  :: PM(M__,2*M__),    PM0(2*M__),   POMEGA(2*M__,N__)
      REAL*8  :: ZTAU(N__),        FZ(N__),      FJ(N__)
      REAL*8  :: DD(M__,M__,N__),  RR(M__,N__)
      REAL*8  :: ZREFL,            ZFLUX

      ! Scalars
      REAL*8  :: RADIUS,           ZU0
      INTEGER :: ND,               N
      INTEGER :: M,                MFIT

      !=================================================================
      ! Declare the following global variables as THREADPRIVATE for the
      ! OpenMP parallelization on all platforms (bmy, 3/23/03)
      !=================================================================
!$OMP THREADPRIVATE( A,B,C1,H,AA,CC,S,W,U1,V1,WT,EMU,PM,PM0,POMEGA )
!$OMP THREADPRIVATE( ZTAU,FZ,FJ,DD,RR,ZREFL,ZFLUX,RADIUS,ZU0 )
!$OMP THREADPRIVATE( ND,N,M,MFIT )
!
! !REMARKS:
!  This file is an amalgam of the original GEOS-Chem files CMN_FJ_MOD,
!  JV_CMN_MOD and JV_MIE_MOD.
!
!  Based on code from Oliver Wild (9 Jul 1999)
!
!  Keep increasing NL if necessary.  This will avoid the "too many levels
!  in photolysis" error.
!
!  NOTES for CTM Interface (bmy, 10/27/99, 3/23/03)
!  =====================================================================
!  (1) Change JPNL and JPPJ from parameters to variables, which are 
!      set in "inphot.f".  This allows the user to switch the number 
!      of levels at run-time via the CTM inputs. 
!                                                                             .
!  (2) Now make RAD, ZZHT, DTAUMAX, DTAUSUB, DSUBDIV, SZAMAX into
!      parameters instead of holding them in common blocks.  
!                                                                             .
!  (3) Create new common blocks /WLLOC/ and /JVLOC/ to hold certain
!      quantities -Xlocal for parallel code (ppm, 4/98, bmy, 9/21/99)
!                                                                             .
!  (4) The common blocks that must be held -Xlocal are:
!         /ATMOS/, /JVSUB/, /WLLOC/, /JVLOC/ 
!                                                                             .
!  (4a) Declare the above commons THREADPRIVATE for the Compaq
!       Alpha platform (bmy, 7/10/01)
!                                                                             .
!  (5) Break MIEDX off from the WLLOC common block, since it must
!      not be declared LOCAL for the parallelization. (bmy, 5/2/00)
!                                                                             .
!  (6) For including aerosol optical depths: (rvm, bmy, 9/30/00)
!      (a) Increase MX from 3 to 10 .  
!      (c) Add ODMDUST(IIPAR,JJPAR,LLPAR,NDUST) to common block /CLIM/
!                                                                             .
!  (7) Move NDUST to CMN_SIZE to avoid conflicts (bmy, 11/15/01)
!                                                                             .
!  (8) For updating aerosol optical depths again (rvm, bmy, 2/27/02):
!      (a) Change NP from 21 to 56
!      (b) Change MX from 10 to 35
!      (c) Add ODAER(IIPAR,JJPAR,LLPAR,NAER*NRH) to common block /CLIM/
!                                                                             .
!  (9) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!       compilation.  Also added & continuation characters in column 73
!       to allow header files to be included in F90 freeform files.
!       Also changed comment character from "C" to "!" to allow this
!       file to be inlined into freeform source code. (bmy, 6/25/02)
!                                                                             .
!  (10) Renamed cpp switch from DEC_COMPAQ to COMPAQ.  Also declare common
!        blocks ATMOS, JVLOC, WLLOC, JVSUB as !$OMP THREADPRIVATE for
!        all platforms. (bmy, 3/23/03)
!                                                                             .
!  (11) Added new pressure denpendencies algorithm parameters 
!         for MGLY. (tmf, 1/7/09)
!                                                                             .
!  (12) Added 'pdepf' as pressure dependancy function selector. (tmf, 1/31/06)
!                                                                             .
!  (13) Split off PDEPF and MGLYPDEP into separate common blocks to avoid
!        warnings on IFORT 9 (ccarouge, bmy, 8/20/09)
!                                                                             .
!  (14) Add new optical variables for AOD calculation (clh, 05/06/10)
!
! !REVISION HISTORY:
!  (1 ) Uses Fortran 90 declarations for parameters and variables
!  (2 ) Pass CTM size parameters and preprocessor switches via CMN_SIZE.
!  (3 ) Update JPMAX for new chemistry mechanism (amf, bmy, 4/20/00)
!  (4 ) Return JPMAX to original setting (bmy, 9/25/00)
!  (5 ) Return JPMAX to 55 for peroxy recycling (again) (bmy, 12/20/00)
!  (6 ) Now need to use the window parameters IIPAR,JJPAR,LLPAR (bmy, 9/25/01)
!  (7 ) Changed RCS ID tags to by adding a ! comment character to allow
!        freeform compilation.  Also added & continuation characters in 
!        column 73 to allow header files to be included in F90 freeform files.
!        Also changed comment character from "C" to "!", to allow this
!        file to be inlined into freeform source code. (bmy, 6/25/02)
!  (8 ) Now declare common blocks /MIEBLK/ and /MINDEX/ as THREADPRIVATE for
!        all platforms (bmy, 3/23/03)
!  (9 ) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!        compilation. (bmy, 6/25/02)
!  (10) Replaced ESIG array with ETAA and ETAB arrays for the hybrid
!        pressure formulation.  Also deleted PREST, since we don't need that
!        anymore. (bmy, 8/23/02)
!  (11) Removed ETAA and ETAB arrays.  We now compute PJ directly from the 
!        GET_PEDGE routine. (bmy, 10/30/07)
!  (12) Set NL to 1000 to avoid SMVGEAR crash with GEOS-5.2.0 (09/01/08)
!  (13) Increase photolysis rxns JPMAX = 79 (tmf, 1/7/09)
!  (14) Increase photolysis rxns JPMAX = 89 for Isoprene (fp, 2/2/10)
!  (15) Increase species name length. (fp, 2/2/10)
!  10 Aug 2012 - R. Yantosca - Replace IPAR, JPAR, LPAR w/ IIPAR, JJPAR, LLPAR
!  27 Mar 2013 - S. D. Eastham - Combine CMN_FJ_MOD, JV_CMN and JV_MIE_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
#if defined( DEVEL )
      CONTAINS

      SUBROUTINE SET_CMN_FJX_MOD

      IMPLICIT NONE

      INTEGER AS

      ALLOCATE(ZPJ(LLPAR,JPMAX,IIPAR,JJPAR), STAT=AS)

      ALLOCATE(
!     & TJ(NB),PJ(NB+1),DM(NB),DO3(NB),DBC(NB),Z(NB),
!     & AER(MX,NB),AMF(NB,NB), zj(LLPAR,JPMAX),
!     & FFF(NW,LLPAR), jadsub(NC),
     & ODMDUST(IIPAR,JJPAR,LLPAR,NDUST),
     & ODAER(IIPAR,JJPAR,LLPAR,NAER*NRH),
     & STAT=AS)

      END SUBROUTINE SET_CMN_FJX_MOD
#endif
      END MODULE CMN_FJX_MOD
!EOC
