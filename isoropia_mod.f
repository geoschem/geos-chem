! $Id: isoropia_mod.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      MODULE ISOROPIA_MOD
!
!******************************************************************************
!
!  *** ISORROPIA CODE
!  *** INCLUDE FILE 'ISRPIA.INC'
!  *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!      AND VARIABLES. 
! 
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "sulfate_mod.f"
      !=================================================================

      ! Make everything private!
      PRIVATE

      ! Make the following routines public
      PUBLIC :: AERO_THERMO, ISOROPIA

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      INTEGER,           PARAMETER   :: IMAX   = 741
      INTEGER,           PARAMETER   :: JMAX   = 13
      INTEGER,           PARAMETER   :: NCOMP  = 5
      INTEGER,           PARAMETER   :: NIONS  = 7
      INTEGER,           PARAMETER   :: NGASAQ = 3
      INTEGER,           PARAMETER   :: NSLDS  = 9  
      INTEGER,           PARAMETER   :: NPAIR  = 13
      INTEGER,           PARAMETER   :: NZSR   = 100
      INTEGER,           PARAMETER   :: NERRMX = 25
      INTEGER,           PARAMETER   :: NSO4S  = 14
      INTEGER,           PARAMETER   :: NRHS   = 20
      INTEGER,           PARAMETER   :: NASRD  = NSO4S * NRHS

      ! Allocatable arrays
      REAL*8,            ALLOCATABLE :: ASRAT(:)
      REAL*8,            ALLOCATABLE :: ASSO4(:)
      REAL*8,            ALLOCATABLE :: AWAB(:)
      REAL*8,            ALLOCATABLE :: AWAC(:) 
      REAL*8,            ALLOCATABLE :: AWAN(:) 
      REAL*8,            ALLOCATABLE :: AWAS(:)
      REAL*8,            ALLOCATABLE :: AWLC(:)
      REAL*8,            ALLOCATABLE :: AWSA(:)
      REAL*8,            ALLOCATABLE :: AWSB(:)
      REAL*8,            ALLOCATABLE :: AWSC(:)
      REAL*8,            ALLOCATABLE :: AWSS(:)
      REAL*8,            ALLOCATABLE :: AWSN(:)
      REAL*4,            ALLOCATABLE :: BNC198(:,:)
      REAL*4,            ALLOCATABLE :: BNC223(:,:)
      REAL*4,            ALLOCATABLE :: BNC248(:,:)      
      REAL*4,            ALLOCATABLE :: BNC273(:,:)
      REAL*4,            ALLOCATABLE :: BNC298(:,:)
      REAL*4,            ALLOCATABLE :: BNC323(:,:)
      INTEGER,           ALLOCATABLE :: ERRSTK(:)
      CHARACTER(LEN=40), ALLOCATABLE :: ERRMSG(:)
      REAL*8,            ALLOCATABLE :: GAMA(:)
      REAL*8,            ALLOCATABLE :: GAMIN(:)
      REAL*8,            ALLOCATABLE :: GAMOU(:)
      REAL*8,            ALLOCATABLE :: GASAQ(:)
      REAL*8,            ALLOCATABLE :: IMW(:)
      REAL*8,            ALLOCATABLE :: M0(:)
      REAL*8,            ALLOCATABLE :: MOLAL(:)
      REAL*8,            ALLOCATABLE :: MOLALR(:)
      REAL*8,            ALLOCATABLE :: SMW(:)
      REAL*8,            ALLOCATABLE :: W(:)
      REAL*8,            ALLOCATABLE :: WMW(:)
      REAL*8,            ALLOCATABLE :: WAER(:)
      REAL*8,            ALLOCATABLE :: Z(:)
      REAL*8,            ALLOCATABLE :: ZZ(:)

      ! Scalar variables
      LOGICAL                        :: CALAIN,   CALAOU
      LOGICAL                        :: DRYF,     FRST,     STKOFL   
      INTEGER                        :: IACALC,   ICLACT,   IPROB 
      INTEGER                        :: MAXIT,    METSTBL,  NDIV
      INTEGER                        :: NSWEEP,   WFTYP,    NOFER
      REAL*4                         :: IONIC
      REAL*8                         :: A1,       A2,       A3
      REAL*8                         :: A4,       A5,       A6
      REAL*8                         :: A7,       A8,       ONE
      REAL*8                         :: CHI1,     CHI2,     CHI3
      REAL*8                         :: CHI4,     CHI5,     CHI6
      REAL*8                         :: CHI7,     CHI8,     GREAT
      REAL*8                         :: CH2SO4,   CNH42S4,  CNH4HS4 
      REAL*8                         :: CNACL,    CNA2SO4,  CNANO3 
      REAL*8                         :: CNH4NO3,  CNH4CL,   CNAHSO4 
      REAL*8                         :: CLC,      CHCL,     CHNO3    
      REAL*8                         :: COH,      EPS,      EPSACT
      REAL*8                         :: GNH3,     GHNO3,    GHCL 
      REAL*8                         :: DRH2SO4,  DRNH42S4, DRNAHSO4
      REAL*8                         :: DRNACL,   DRNANO3,  DRNA2SO4 
      REAL*8                         :: DRNH4HS4, DRLC,     DRNH4NO3
      REAL*8                         :: DRNH4CL,  DRMLCAB,  DRMLCAS
      REAL*8                         :: DRMASAN,  DRMG1,    DRMG2
      REAL*8                         :: DRMG3,    DRMH1,    DRMH2
      REAL*8                         :: DRMI1,    DRMI2     DRMI3
      REAL*8                         :: DRMQ1,    DRMR1,    DRMR2
      REAL*8                         :: DRMR3,    DRMR4,    DRMR5
      REAL*8                         :: DRMR6,    DRMR7,    DRMR8
      REAL*8                         :: DRMR9,    DRMR10,   DRMR11
      REAL*8                         :: DRMR12,   DRMR13  
      REAL*8                         :: PSI1,     PSI2,     PSI3
      REAL*8                         :: PSI4,     PSI5,     PSI6
      REAL*8                         :: PSI7,     PSI8,     R
      REAL*8                         :: RH,       SULRATW,  SULRAT
      REAL*8                         :: SODRAT,   TEMP,     TINY
      REAL*8                         :: TINY2,    WATER,    ZERO
      REAL*8                         :: XK1,      XK2,      XK3
      REAL*8                         :: XK4,      XK5,      XK6
      REAL*8                         :: XK7,      XK8,      XK9
      REAL*8                         :: XK10,     XK11,     XK12
      REAL*8                         :: XK13,     XK14,     XKW
      REAL*8                         :: XK21,     XK22,     XK31
      REAL*8                         :: XK32,     XK41,     XK42
      CHARACTER(LEN=14)              :: VERSION
      CHARACTER(LEN=15)              :: SCASE

      !=================================================================
      ! MODULE ROUTINES -- follow after the CONTAINS statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE AERO_THERMO( HNO3, STT )
!
!***********************************************************************
!  Subroutine AERO_THERMO is the interface between the GEOS-CHEM model
!  and the aerosol thermodynamical equilibrium routine in "rpmares.f"
!  [rjp, 12/17/01]
! 
!  Arguments as Input:
!  =====================================================================
!  (2 ) TEMP    (REAL*8  ) : Air temperature                  [K]
!  (4 ) RH      (REAL*8  ) : Fractional relative humidity
!  (8 ) VOL     (REAL*8  ) : Volume  of air  in a grid box    [m^3]
!  (  ) HNO3    (REAL*8  ) : HNO3 contains monthly HNO3 conc. [kg]
!  (13) STT     (REAL*8  ) : Array for tracer concentrations  [kg]  
!
!  Arguments as Output:
!  =====================================================================
!
!  (13) STT     (REAL*8 ) : STT contains updated concentrations
!
!  NOTES: Aerosol concentrations are all in ug/m^3
!***********************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : AIRVOL, RH, T

      IMPLICIT NONE

#     include "CMN_SIZE"

      ! Kg/box
      REAL*8,  INTENT(INOUT) :: HNO3(IIPAR,JJPAR,LLPAR)     
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,NNPAR)

      ! Local variables for RPMARES
      LOGICAL, SAVE      :: FIRST = .TRUE.
      INTEGER            :: I, J, L
      REAL*8             :: SO4     ! total sulfate aerosol for input ( ug/m^3)
      REAL*8             :: ASO4    ! sulfate aerosol(SO4--)
      REAL*8             :: ANO3    ! nitrate aerosol
      REAL*8             :: AH2O    ! water content in sulfate aerosol
      REAL*8             :: ANH4    ! ammonium in aerosol
      REAL*8             :: GNH3    ! ammonia (GAS)
      REAL*8             :: GNO3    ! nitric acid (gas)
      REAL*8             :: AHSO4   ! bisulfate (AEROSOL)
      REAL*8             :: WI(5),  CNTRL(2),   RHI
      REAL*8             :: TEMPI,  VOL,        WT(5)
      REAL*8             :: GAS(3), AERLIQ(12), AERSLD(9), OTHER(6)
      REAL*8             :: TSO4, TNH3, TNO3
      CHARACTER(LEN=15)  :: SCASE

      ! concentration lower limit [ ug/m**3 ]
      REAL*8,  PARAMETER :: CONMIN = 1.0D-30
                  
      !=================================================================
      ! AERO_THERMO begins here!
      !=================================================================

      ! Initialize module arrays and variables
      IF ( FIRST ) THEN
         CALL INIT_ISOROPIA
         FIRST = .FALSE.
      ENDIF
       
      ! 0 = (SOLID+LIQUID), 1 = METASTABLE(ONLY LIQUID)
      CNTRL(1) = 0D0 
      CNTRL(2) = 1D0         

      !=================================================================
      ! Loop over grid boxes and call ISOROPIA (see comments in the 
      ! ISOROPIA routine below which describe the input/output args)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TEMPI, RHI, VOL, OTHER, SCASE )
!$OMP+PRIVATE( WI, TSO4, TNH3, TNO3, WT, GAS, AERSLD, AERLIQ )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Temperature [K]
         TEMPI  = T(I,J,L)

         ! Relative humidity [unitless]
         RHI    = RH(I,J,L) * 1.d-2

         ! Volume of grid box [m3] 
         VOL    = AIRVOL(I,J,L)

         ! Convert tracers from [kg] to [mole/m3]
         TSO4   = STT(I,J,L,3) * 1.d3 / ( 96.D0 * VOL )
         TNH3   = STT(I,J,L,6) * 1.D3 / ( 18.D0 * VOL )  +
     &            STT(I,J,L,5) * 1.D3 / ( 17.D0 * VOL )
         TNO3   = HNO3(I,J,L)  * 1.D3 / ( 63.D0 * VOL )

         ! Insert into WI array
         WI     = 0D0
         WI(2)  = MAX( TSO4, CONMIN )
         WI(3)  = MAX( TNH3, CONMIN )
         WI(4)  = MAX( TNO3, CONMIN )

         ! Unit for isoropia input is mole/m^3
         CALL ISOROPIA( WI, RHI, TEMPI,  CNTRL,
     &                  WT, GAS, AERLIQ, AERSLD, SCASE, OTHER )

         ! Convert output from ISOROPIA back to [kg/box] and store in STT
         STT(I,J,L,3) = MAX( 96.D-3 * VOL * WT(2),              CONMIN )
         STT(I,J,L,5) = MAX( 17.D-3 * VOL * GAS(1),             CONMIN )
         STT(I,J,L,6) = MAX( 18.D-3 * VOL * ( WT(3) - GAS(1) ), CONMIN )
         STT(I,J,L,7) = MAX( 62.D-3 * VOL * ( WT(4) - GAS(2) ), CONMIN )
         !HNO3(I,J,L ) = MAX(63.D-3*AIRVOL* GAS(2), CONMIN)
         !WH2O(I,J,L ) = 18.D-3 * AIRVOL *   AERLIQ(8)
 
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO 

      ! Return to calling program
      END SUBROUTINE AERO_THERMO

!------------------------------------------------------------------------------

      SUBROUTINE ISOROPIA( WI, RHI, TEMPI,  CNTRL, 
     &                     WT, GAS, AERLIQ, AERSLD, SCASI, OTHER )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISOROPIA
! *** THIS SUBROUTINE IS THE MASTER ROUTINE FOR THE ISORROPIA
!     THERMODYNAMIC EQUILIBRIUM AEROSOL MODEL (VERSION 1.1 and above)
!
! ======================== ARGUMENTS / USAGE ===========================
!
!  INPUT:
!  1. [WI] 
!     DOUBLE PRECISION array of length [5].
!     Concentrations, expressed in moles/m3. Depending on the type of
!     problem solved (specified in CNTRL(1)), WI contains either 
!     GAS+AEROSOL or AEROSOL only concentratios.
!     WI(1) - sodium
!     WI(2) - sulfate
!     WI(3) - ammonium
!     WI(4) - nitrate
!     WI(5) - chloride
!
!  2. [RHI] 
!     DOUBLE PRECISION variable.  
!     Ambient relative humidity expressed on a (0,1) scale.
!
!  3. [TEMPI]
!     DOUBLE PRECISION variable. 
!     Ambient temperature expressed in Kelvins. 
!
!  4. [CNTRL]
!     DOUBLE PRECISION array of length [2].
!     Parameters that control the type of problem solved.
!
!     CNTRL(1): Defines the type of problem solved.
!     0 - Forward problem is solved. In this case, array WI contains 
!         GAS and AEROSOL concentrations together.
!     1 - Reverse problem is solved. In this case, array WI contains
!         AEROSOL concentrations only.
!
!     CNTRL(2): Defines the state of the aerosol
!     0 - The aerosol can have both solid+liquid phases (deliquescent)
!     1 - The aerosol is in only liquid state (metastable aerosol)
!
!  OUTPUT:
!  1. [WT] 
!     DOUBLE PRECISION array of length [5].
!     Total concentrations (GAS+AEROSOL) of species, expressed in moles/m3. 
!     If the foreward probelm is solved (CNTRL(1)=0), array WT is 
!     identical to array WI.
!     WT(1) - total sodium
!     WT(2) - total sulfate
!     WT(3) - total ammonium
!     WT(4) - total nitrate
!     WT(5) - total chloride
!
!  2. [GAS]
!     DOUBLE PRECISION array of length [03]. 
!     Gaseous species concentrations, expressed in moles/m3. 
!     GAS(1) - NH3
!     GAS(2) - HNO3
!     GAS(3) - HCl 
!
!  3. [AERLIQ]
!     DOUBLE PRECISION array of length [11]. 
!     Liquid aerosol species concentrations, expressed in moles/m3. 
!     AERLIQ(01) - H+(aq)          
!     AERLIQ(02) - Na+(aq)         
!     AERLIQ(03) - NH4+(aq)
!     AERLIQ(04) - Cl-(aq)         
!     AERLIQ(05) - SO4--(aq)       
!     AERLIQ(06) - HSO4-(aq)       
!     AERLIQ(07) - NO3-(aq)        
!     AERLIQ(08) - H2O             
!     AERLIQ(09) - NH3(aq) (undissociated)
!     AERLIQ(10) - HNCl(aq) (undissociated)
!     AERLIQ(11) - HNO3(aq) (undissociated)
!     AERLIQ(12) - OH-(aq)
!
!  4. [AERSLD]
!     DOUBLE PRECISION array of length [09]. 
!     Solid aerosol species concentrations, expressed in moles/m3. 
!     AERSLD(01) - NaNO3(s)
!     AERSLD(02) - NH4NO3(s)
!     AERSLD(03) - NaCl(s)         
!     AERSLD(04) - NH4Cl(s)
!     AERSLD(05) - Na2SO4(s)       
!     AERSLD(06) - (NH4)2SO4(s)
!     AERSLD(07) - NaHSO4(s)
!     AERSLD(08) - NH4HSO4(s)
!     AERSLD(09) - (NH4)4H(SO4)2(s)
!
!  5. [SCASI]
!     CHARACTER*15 variable.
!     Returns the subcase which the input corresponds to.
!
!  6. [OTHER]
!     DOUBLE PRECISION array of length [6].
!     Returns solution information.
!
!     OTHER(1): Shows if aerosol water exists.
!     0 - Aerosol is WET
!     1 - Aerosol is DRY
!
!     OTHER(2): Aerosol Sulfate ratio, defined as (in moles/m3) :
!               (total ammonia + total Na) / (total sulfate)
!
!     OTHER(3): Sulfate ratio based on aerosol properties that defines 
!               a sulfate poor system:
!               (aerosol ammonia + aerosol Na) / (aerosol sulfate)
!           
!     OTHER(4): Aerosol sodium ratio, defined as (in moles/m3) :
!               (total Na) / (total sulfate)
!      
!     OTHER(5): Ionic strength of the aqueous aerosol (if it exists).
!      
!     OTHER(6): Total number of calls to the activity coefficient 
!               calculation subroutine.
! 
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Parameters
      INTEGER, PARAMETER    :: NCTRL=2, NOTHER=6

      ! Arguments
      REAL*8, INTENT(IN)    :: RHI, TEMPI
      REAL*8, INTENT(IN)    :: WI(NCOMP)
      REAL*8, INTENT(IN)    :: CNTRL(NCTRL)
      REAL*8, INTENT(OUT)   :: WT(NCOMP)
      REAL*8, INTENT(OUT)   :: GAS(NGASAQ)
      REAL*8, INTENT(OUT)   :: AERSLD(NSLDS) 
      REAL*8, INTENT(OUT)   :: AERLIQ(NIONS+NGASAQ+2)
      REAL*8, INTENT(OUT)   :: OTHER(NOTHER)

      ! Local variables
      INTEGER               :: I
      CHARACTER(LEN=15)     :: SCASI

      !=================================================================
      ! ISOROPIA begins here!
      !=================================================================

      ! PROBLEM TYPE (0=FOREWARD, 1=REVERSE)
      IPROB   = NINT(CNTRL(1))

      ! AEROSOL STATE (0=SOLID+LIQUID, 1=METASTABLE
      METSTBL = NINT(CNTRL(2))

      !=================================================================
      ! SOLVE FORWARD PROBLEM
      !=================================================================
      CALL ISRP2F( WI, RHI, TEMPI )

      ! IF METASTABLE AND NO WATER - RESOLVE AS NORMAL
      !IF (WATER.LE.TINY .AND. METSTBL.EQ.1) THEN
      !   METSTBL = 0
      !   GOTO 50
      !ENDIF

      !=================================================================
      ! Save results to arrays (units = mole/m3)
      !=================================================================

      ! Gaseous aerosol species
      GAS(1) = GNH3                
      GAS(2) = GHNO3
      GAS(3) = GHCL

      ! Liquid aerosol species
      DO I = 1, NIONS              
         AERLIQ(I) = MOLAL(I)
      ENDDO

      DO I = 1, NGASAQ
         AERLIQ(NIONS+1+I) = GASAQ(I)
      ENDDO

      AERLIQ(NIONS+1)        = WATER*1.0D3/18.0D0
      AERLIQ(NIONS+NGASAQ+2) = COH

      ! Solid aerosol species
      AERSLD(1) = CNANO3           
      AERSLD(2) = CNH4NO3
      AERSLD(3) = CNACL
      AERSLD(4) = CNH4CL
      AERSLD(5) = CNA2SO4
      AERSLD(6) = CNH42S4
      AERSLD(7) = CNAHSO4
      AERSLD(8) = CNH4HS4
      AERSLD(9) = CLC

      ! Dry flag
      IF(WATER.LE.TINY) THEN    
         OTHER(1) = 1.d0
      ELSE
         OTHER(1) = 0.d0
      ENDIF

      ! Other stuff
      OTHER(2) = SULRAT         
      OTHER(3) = SULRATW
      OTHER(4) = SODRAT
      OTHER(5) = IONIC
      OTHER(6) = ICLACT

      SCASI = SCASE

      ! Total gas+aerosol phase
      WT(1) = WI(1)             
      WT(2) = WI(2)
      WT(3) = WI(3) 
      WT(4) = WI(4)
      WT(5) = WI(5)

      IF (IPROB.GT.0 .AND. WATER.GT.TINY) THEN 
         WT(3) = WT(3) + GNH3 
         WT(4) = WT(4) + GHNO3
         WT(5) = WT(5) + GHCL
      ENDIF

      ! Return to calling program
      END SUBROUTINE ISOROPIA

!------------------------------------------------------------------------------

      SUBROUTINE INIT2( WI, RHI, TEMPI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE INIT2
! *** THIS SUBROUTINE INITIALIZES ALL GLOBAL VARIABLES FOR AMMONIUM,
!     NITRATE, SULFATE AEROSOL SYSTEMS (SUBROUTINE ISRP2)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8          :: WI(NCOMP), RHI, TEMPI

      ! Local variables ??
      INTEGER         :: I,  IRH
      REAL            :: IC, GII, GI0, XX
      REAL, PARAMETER :: LN10=2.3025851
      REAL*8          :: T0, T0T, COEF, TCF, G130, G13I

      !=================================================================
      ! INIT2 begins here!
      !=================================================================

      ! Save input variables in common block
      IF (IPROB.EQ.0) THEN                 

         ! Forward calculation
         DO I=1,NCOMP
            W(I) = MAX(WI(I), TINY)
         ENDDO

      ELSE
         
         ! Reverse calculation
         DO 15 I=1,NCOMP                   
            WAER(I) = MAX(WI(I), TINY)
            W(I)    = ZERO
15       CONTINUE
      ENDIF

      RH      = RHI
      TEMP    = TEMPI

      !=================================================================
      ! Calculate equilibrium constants
      !=================================================================

      XK1  = 1.015e-2  ! HSO4(aq)         <==> H(aq)     + SO4(aq)
      XK21 = 57.639    ! NH3(g)           <==> NH3(aq)
      XK22 = 1.805e-5  ! NH3(aq)          <==> NH4(aq)   + OH(aq)
      XK4  = 2.511e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! ISORR
      !XK4  = 3.638e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! SEQUIL
      XK41 = 2.100e5   ! HNO3(g)          <==> HNO3(aq)
      XK7  = 1.817     ! (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
      XK10 = 5.746e-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! ISORR
      !XK10 = 2.985e-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! SEQUIL
      XK12 = 1.382e2   ! NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
      XK13 = 29.268    ! (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
      XKW  = 1.010e-14 ! H2O              <==> H(aq)     + OH(aq)

      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0  = 298.15D0
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK4 = XK4 *EXP( 29.17*(T0T-1.0) + 16.830*COEF) !ISORR
         !XK4 = XK4 *EXP( 29.47*(T0T-1.0) + 16.840*COEF) ! SEQUIL
         XK41= XK41*EXP( 29.17*(T0T-1.0) + 16.830*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK10= XK10*EXP(-74.38*(T0T-1.0) +  6.120*COEF) ! ISORR
         !XK10= XK10*EXP(-75.11*(T0T-1.0) + 13.460*COEF) ! SEQUIL
         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF

      XK2  = XK21*XK22       
      XK42 = XK4/XK41

      !=================================================================
      ! Calculate deliquescence relative humidities (unicomponent)
      !=================================================================
      DRH2SO4  = ZERO
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRNH4NO3 = 0.6183D0
      DRLC     = 0.6900D0

      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15D0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNH4NO3 = DRNH4NO3*EXP(852.*TCF)
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF) 
         DRLC     = DRLC    *EXP(186.*TCF) 
      ENDIF

      !=================================================================
      ! Calculate mutual deliquescence relative humidities
      !=================================================================
      DRMLCAB = 0.3780D0              ! (NH4)3H(SO4)2 & NH4HSO4 
      DRMLCAS = 0.6900D0              ! (NH4)3H(SO4)2 & (NH4)2SO4 
      DRMASAN = 0.6000D0              ! (NH4)2SO4     & NH4NO3

      ! comment out for time being
      !IF (INT(TEMP) .NE. 298) THEN    ! For the time being
      !   T0       = 298.15d0
      !   TCF      = 1.0/TEMP - 1.0/T0
      !   DRMLCAB  = DRMLCAB*EXP( 507.506*TCF) 
      !   DRMLCAS  = DRMLCAS*EXP( 133.865*TCF) 
      !   DRMASAN  = DRMASAN*EXP(1269.068*TCF)
      !ENDIF

      !=================================================================
      ! Liquid phase initialization
      !=================================================================
      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY

      DO I = 1, NPAIR
         MOLALR(I) = ZERO
         GAMA(I)   = 0.1
         GAMIN(I)  = GREAT
         GAMOU(I)  = GREAT
         M0(I)     = 1d5
      ENDDO

      DO I = 1, NPAIR
         GAMA(I) = 0.1d0
      ENDDO

      DO I = 1, NIONS
         MOLAL(I) = ZERO
      ENDDO
      COH = ZERO

      DO I = 1, NGASAQ
         GASAQ(I) = ZERO
      ENDDO

      !=================================================================
      ! Solid phase initialization
      !=================================================================
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CNACL   = ZERO
      CNA2SO4 = ZERO
      CNANO3  = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNAHSO4 = ZERO
      CLC     = ZERO

      !=================================================================
      ! Gas phase
      !=================================================================
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO

      !=================================================================
      ! Calculate ZSR parameters
      !=================================================================
      IRH     = MIN (INT(RH*NZSR+0.5),NZSR)  ! Position in ZSR arrays
      IRH     = MAX (IRH, 1)

      ! NACl
      M0(01) = AWSC(IRH)      
      IF (M0(01) .LT. 100.0) THEN
         IC = M0(01)
         CALL KMTAB(IC,298.0,     GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,REAL(TEMP),GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(01) = M0(01)*EXP(LN10*(GI0-GII))
      ENDIF

      ! (NA)2SO4
      M0(02) = AWSS(IRH)      
      IF (M0(02) .LT. 100.0) THEN
         IC = 3.0*M0(02)
         CALL KMTAB(IC,298.0,     XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,REAL(TEMP),XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(02) = M0(02)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NANO3
      M0(03) = AWSN(IRH)      
      IF (M0(03) .LT. 100.0) THEN
         IC = M0(03)
         CALL KMTAB(IC,298.0,     XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,REAL(TEMP),XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(03) = M0(03)*EXP(LN10*(GI0-GII))
      ENDIF

      ! (NH4)2SO4
      M0(04) = AWAS(IRH)      
      IF (M0(04) .LT. 100.0) THEN
         IC = 3.0*M0(04)
         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
         M0(04) = M0(04)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NH4NO3
      M0(05) = AWAN(IRH)      
      IF (M0(05) .LT. 100.0) THEN
         IC     = M0(05)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
         M0(05) = M0(05)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NH4CL
      M0(06) = AWAC(IRH)      
      IF (M0(06) .LT. 100.0) THEN
         IC = M0(06)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
         M0(06) = M0(06)*EXP(LN10*(GI0-GII))
      ENDIF

      ! 2H-SO4
      M0(07) = AWSA(IRH)      
      IF (M0(07) .LT. 100.0) THEN
         IC = 3.0*M0(07)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX)
         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
         M0(07) = M0(07)*EXP(LN10*(GI0-GII))
      ENDIF

      ! H-HSO4
      M0(08) = AWSA(IRH)      
      !### These are redundant, because M0(8) is not used
      !IF (M0(08) .LT. 100.0) THEN     
      !   IC = M0(08)
      !   CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
      !   CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
      !   CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX)
      !   M0(08) = M0(08)*EXP(LN10*(GI0-GII))
      !ENDIF

      ! NH4HSO4
      M0(09) = AWAB(IRH)      
      IF (M0(09) .LT. 100.0) THEN
         IC = M0(09)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX)
         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
         M0(09) = M0(09)*EXP(LN10*(GI0-GII))
      ENDIF

      ! NAHSO4
      M0(12) = AWSB(IRH)      
      IF (M0(12) .LT. 100.0) THEN
         IC = M0(12)
         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GI0)
         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
         M0(12) = M0(12)*EXP(LN10*(GI0-GII))
      ENDIF

      ! (NH4)3H(SO4)2
      M0(13) = AWLC(IRH)      
      IF (M0(13) .LT. 100.0) THEN
         IC     = 4.0*M0(13)
         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G130   = 0.2*(3.0*GI0+2.0*GII)
         CALL KMTAB(IC,REAL(TEMP),XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
         G13I   = 0.2*(3.0*GI0+2.0*GII)
         M0(13) = M0(13)*EXP(LN10*SNGL(G130-G13I))
      ENDIF

      !=================================================================
      ! OTHER INITIALIZATIONS
      !=================================================================
      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = '??'
      SULRATW = 2.D0
      SODRAT  = ZERO
      NOFER   = 0
      STKOFL  =.FALSE.

      DO I = 1, NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT2

!-----------------------------------------------------------------------------

      FUNCTION GETASR( SO4I, RHI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION GETASR
! *** CALCULATES THE LIMITING NH4+/SO4 RATIO OF A SULFATE POOR SYSTEM
!     (i.e. SULFATE RATIO = 2.0) FOR GIVEN SO4 LEVEL AND RH
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: SO4I, RHI

      ! Local variables
      INTEGER            :: IA1,   INDS,  INDR, INDSL
      INTEGER            :: INDSH, IPOSL, IPOSH  
      REAL*8             :: GETASR, RAT,  A1, WF
            
      !=================================================================
      ! GETASR begins here!
      !=================================================================

      !### SOLVE USING FULL COMPUTATIONS, NOT LOOK-UP TABLES
      !W(2) = WAER(2)
      !W(3) = WAER(2)*2.0001D0
      !CALL CALCA2
      !SULRATW = MOLAL(3)/WAER(2)
      !CALL INIT1 (WI, RHI, TEMPI) ! Re-initialize COMMON BLOCK

      ! CALCULATE INDICES
      RAT    = SO4I/1.E-9    
      A1     = INT(ALOG10(RAT))                   ! Magnitude of RAT
      IA1    = INT(RAT/2.5/10.0**A1)

      INDS   = 4.0*A1 + MIN(IA1,4)
      INDS   = MIN(MAX(0, INDS), NSO4S-1) + 1     ! SO4 component of IPOS

      INDR   = INT(99.0-RHI*100.0) + 1
      INDR   = MIN(MAX(1, INDR), NRHS)            ! RH component of IPOS

      ! GET VALUE AND RETURN
      INDSL  = INDS
      INDSH  = MIN(INDSL+1, NSO4S)
      IPOSL  = (INDSL-1)*NRHS + INDR              ! Low position in array
      IPOSH  = (INDSH-1)*NRHS + INDR              ! High position in array

      WF     = (SO4I-ASSO4(INDSL))/(ASSO4(INDSH)-ASSO4(INDSL) + 1e-7)
      WF     = MIN(MAX(WF, 0.0), 1.0)

      GETASR = WF*ASRAT(IPOSH) + (1.0-WF)*ASRAT(IPOSL)

      ! Return to calling program
      END FUNCTION GETASR

!------------------------------------------------------------------------------

      SUBROUTINE CALCNA
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNA
! *** CALCULATES NITRATES SPECIATION
!
!     NITRIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES, THAT 
!     DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. THE NITRIC
!     ACID DISSOLVED IS CALCULATED FROM THE HNO3(G) -> (H+) + (NO3-) 
!     EQUILIBRIUM, USING THE (H+) FROM THE SULFATES.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      REAL*8 :: KAPA, X, DELT, ALFA, DIAK

      !### comment out for now
      !CHARACTER ERRINF*40

      !=================================================================
      ! CALCNA begins here!
      !=================================================================

      ! CALCULATE HNO3 DISSOLUTION
      X    = W(4) 
      DELT = 0.0d0
      IF (WATER.GT.TINY) THEN
         KAPA = MOLAL(1)
         ALFA = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         DIAK = SQRT( (KAPA+ALFA)**2.0 + 4.0*ALFA*X)
         DELT = 0.5*(-(KAPA+ALFA) + DIAK)

         !### Comment out for now
         !IF (DELT/KAPA.GT.0.1d0) THEN
         !   WRITE (ERRINF,'(1PE10.3)') DELT/KAPA*100.0
         !   CALL PUSHERR (0019, ERRINF)    ! WARNING ERROR: NO SOLUTION
         !ENDIF
      ENDIF

      ! CALCULATE HNO3 SPECIATION IN THE GAS PHASE 
      GHNO3    = MAX(X-DELT, 0.0d0)  ! GAS HNO3

      ! CALCULATE HNO3 SPECIATION IN THE LIQUID PHASE
      MOLAL(7) = DELT                ! NO3-
      MOLAL(1) = MOLAL(1) + DELT     ! H+ 

      ! Return to calling program
      END SUBROUTINE CALCNA

!------------------------------------------------------------------------------

      SUBROUTINE CALCNH3
!
!=======================================================================
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNH3
! *** CALCULATES AMMONIA IN GAS PHASE
!
!     AMMONIA IN THE GAS PHASE IS ASSUMED A MINOR SPECIES, THAT 
!     DOES NOT SIGNIFICANTLY PERTURB THE AEROSOL EQUILIBRIUM. 
!     AMMONIA GAS IS CALCULATED FROM THE NH3(g) + (H+)(l) <==> (NH4+)(l)
!     EQUILIBRIUM, USING (H+), (NH4+) FROM THE AEROSOL SOLUTION.
!
!     THIS IS THE VERSION USED BY THE DIRECT PROBLEM
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      ! Local variables
      REAL*8 :: A1, CHI1, CHI2, BB, CC, DIAK, PSI, XK42

      !=================================================================
      ! CALCNH3 begins here!
      !=================================================================

      ! Is there a liquid phase?
      IF ( WATER .LE. TINY ) RETURN

      ! Calculate NH3 sublimation
      A1   =  (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      CHI1 =  MOLAL(3)
      CHI2 =  MOLAL(1)

      BB   = (CHI2 + ONE/A1)          ! a=1; b!=1; c!=1 
      CC   = -CHI1/A1             
      DIAK =  SQRT(BB*BB - 4.D0*CC)   ! Always > 0
      PSI  =  0.5*(-BB + DIAK)        ! One positive root
      PSI  =  MAX(TINY, MIN(PSI,CHI1))! Constrict in acceptible range

      ! Calculate NH3 speciation in the gas phase
      GNH3     = PSI                 ! GAS HNO3

      ! Calculate NH3 affect in the liquid phase
      MOLAL(3) = CHI1 - PSI          ! NH4+
      MOLAL(1) = CHI2 + PSI          ! H+ 
 
      ! Return to calling program
      END SUBROUTINE CALCNH3

!------------------------------------------------------------------------------

      SUBROUTINE CALCNIAQ( NO3I, HI, DELT )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNIAQ
!
!     THIS SUBROUTINE CALCULATES THE HNO3(aq) GENERATED FROM (H,NO3-).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8 :: NO3I, HI, DELT

      ! Local variables
      REAL*8 :: A42, OM1, OM2, BB, CC, DD, DEL1, DEL2

      !=================================================================
      ! CALCNIAQ begins here!
      !=================================================================

      ! equilibrium constants
      A42  = XK42*WATER/(GAMA(10))**2. ! gama(hno3) assumed 1

      ! Find root
      OM1  =  NO3I          
      OM2  =  HI
      BB   = -(OM1+OM2+A42)
      CC   =  OM1*OM2
      DD   =  SQRT(BB*BB-4.D0*CC)

      DEL1 =  0.5D0*(-BB - DD)
      DEL2 =  0.5D0*(-BB + DD)

      ! Get appropriate root.
      IF (DEL1.LT.ZERO .OR. DEL1.GT.HI .OR. DEL1.GT.NO3I) THEN
         DELT = ZERO
      ELSE
         DELT = DEL1
         RETURN
      ENDIF

      IF (DEL2.LT.ZERO .OR. DEL2.GT.NO3I .OR. DEL2.GT.HI) THEN
         DELT = ZERO
      ELSE
         DELT = DEL2
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCNIAQ

!------------------------------------------------------------------------------

      SUBROUTINE CALCNIAQ2( GGNO3, NO3I, HI, NO3AQ )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCNIAQ2
!
!     THIS SUBROUTINE CALCULATES THE UNDISSOCIATED HNO3(aq)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8 :: NO3I, NO3AQ, GGNO3, HI 

      ! Local variables
      REAL*8 :: A42, AKW, ALF1, ALF2, ALF3, BB, CC, DEL1

      !=================================================================
      ! CALCNIAQ2 begins here!
      !=================================================================

      ! Equilibrium constants
      A42  = XK42*WATER/(GAMA(10))**2. ! GAMA(HNO3) assumed 1
      AKW  = XKW *RH*WATER*WATER

      ! Find root
      ALF1  = NO3I - GGNO3
      ALF2  = GGNO3
      ALF3  = HI

      BB    = ALF3 + ALF1 + A42
      CC    = ALF3*ALF1 - A42*ALF2
      DEL1  = 0.5*(-BB + SQRT(BB*BB-4.D0*CC))

      ! Correct concentrations
      NO3I  = ALF1 + DEL1
      HI    = ALF3 + DEL1
      IF (HI.LE.TINY) HI = SQRT(AKW)   ! If solution is neutral.
      NO3AQ = ALF2 - DEL1

      ! Return to calling program
      END SUBROUTINE CALCNIAQ2

!------------------------------------------------------------------------------

      SUBROUTINE CALCMR
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCMR
! *** THIS SUBROUTINE CALCULATES:
!     1. ION PAIR CONCENTRATIONS (FROM [MOLAR] ARRAY)
!     2. WATER CONTENT OF LIQUID AEROSOL PHASE (FROM ZSR CORRELATION)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      REAL*8           :: SO4I, HSO4I, AML5, TOTS4, FRNH4, FRNO3, FRCL
      INTEGER          :: I
      CHARACTER(LEN=1) :: SC

      !=================================================================
      ! CALCMR begins here!
      !=================================================================

      ! CALCULATE ION PAIR CONCENTRATIONS ACCORDING TO SPECIFIC CASE
      SC =SCASE(1:1)                   ! SULRAT & SODRAT case

      !=================================================================
      ! NH4-SO4 SYSTEM ; SULFATE POOR CASE
      !=================================================================
      IF (SC.EQ.'A') THEN      
         MOLALR(4) = MOLAL(5)+MOLAL(6) ! (NH4)2SO4 - CORRECT FOR SO4 TO HSO4

      !=================================================================
      ! NH4-SO4 SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'B') THEN
         SO4I  = MOLAL(5)-MOLAL(1)     ! CORRECT FOR HSO4 DISSOCIATION 
         HSO4I = MOLAL(6)+MOLAL(1)              
         IF (SO4I.LT.HSO4I) THEN                
            MOLALR(13) = SO4I                   ! [LC] = [SO4]       
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)  ! NH4HSO4
         ELSE                                   
            MOLALR(13) = HSO4I                  ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)  ! (NH4)2SO4
         ENDIF

      !=================================================================
      ! NH4-SO4 SYSTEM ; SULFATE RICH CASE ; FREE ACID 
      !=================================================================
      ELSE IF (SC.EQ.'C') THEN
         MOLALR(4) = MOLAL(3)                     ! NH4HSO4
         MOLALR(7) = MAX(W(2)-W(3), ZERO)         ! H2SO4

      !=================================================================
      ! NH4-SO4 SYSTEM ; SULFATE RICH CASE ; FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'D') THEN      
         MOLALR(4) = MOLAL(5) + MOLAL(6)          ! (NH4)2SO4
         AML5      = MOLAL(3)-2.D0*MOLALR(4)      ! "free" NH4
         MOLALR(5) = MAX(MIN(AML5,MOLAL(7)), ZERO)! NH4NO3 = MIN("free", NO3)

      !=================================================================
      ! NH4-SO4-NO3 SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'E') THEN      
         SO4I  = MAX(MOLAL(5)-MOLAL(1),ZERO)      ! FROM HSO4 DISSOCIATION 
         HSO4I = MOLAL(6)+MOLAL(1)              
         IF (SO4I.LT.HSO4I) THEN                
            MOLALR(13) = SO4I                     ! [LC] = [SO4] 
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)    ! NH4HSO4
         ELSE                                   
            MOLALR(13) = HSO4I                    ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)    ! (NH4)2SO4
         ENDIF

      !=================================================================
      ! NH4-SO4-NO3 SYSTEM ; SULFATE RICH CASE ; FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'F') THEN      
         MOLALR(4) = MOLAL(3)                              ! NH4HSO4
         MOLALR(7) = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3),ZERO)  ! H2SO4

      !=================================================================
      ! NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE POOR ; SODIUM POOR CASE
      !=================================================================
      ELSE IF (SC.EQ.'G') THEN      
         MOLALR(2) = 0.5*MOLAL(2)                          ! NA2SO4
         TOTS4     = MOLAL(5)+MOLAL(6)                     ! Total SO4
         MOLALR(4) = MAX(TOTS4 - MOLALR(2), ZERO)          ! (NH4)2SO4
         FRNH4     = MAX(MOLAL(3) - 2.D0*MOLALR(4), ZERO)
         MOLALR(5) = MIN(MOLAL(7),FRNH4)                   ! NH4NO3
         FRNH4     = MAX(FRNH4 - MOLALR(5), ZERO)
         MOLALR(6) = MIN(MOLAL(4), FRNH4)                  ! NH4CL

      !=================================================================
      ! NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE POOR ; SODIUM RICH CASE
      ! RETREIVE DISSOLVED SALTS DIRECTLY FROM COMMON BLOCK /SOLUT/
      !=================================================================
      ELSE IF (SC.EQ.'H') THEN      
         MOLALR(1) = PSI7                                  ! NACL 
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(3) = PSI8                                  ! NANO3
         MOLALR(4) = ZERO                                  ! (NH4)2SO4
         FRNO3     = MAX(MOLAL(7) - MOLALR(3), ZERO)       ! "FREE" NO3
         FRCL      = MAX(MOLAL(4) - MOLALR(1), ZERO)       ! "FREE" CL
         MOLALR(5) = MIN(MOLAL(3),FRNO3)                   ! NH4NO3
         FRNH4     = MAX(MOLAL(3) - MOLALR(5), ZERO)       ! "FREE" NH3
         MOLALR(6) = MIN(FRCL, FRNH4)                      ! NH4CL

      !=================================================================
      ! NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
      ! RETREIVE DISSOLVED SALTS DIRECTLY FROM COMMON BLOCK /SOLUT/
      !=================================================================
      ELSE IF (SC.EQ.'I') THEN      
         MOLALR(04) = PSI5                                 ! (NH4)2SO4
         MOLALR(02) = PSI4                                 ! NA2SO4
         MOLALR(09) = PSI1                                 ! NH4HSO4
         MOLALR(12) = PSI3                                 ! NAHSO4
         MOLALR(13) = PSI2                                 ! LC

      !=================================================================
      ! NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE RICH CASE ; FREE ACID
      !=================================================================
      ELSE IF (SC.EQ.'J') THEN      
         MOLALR(09) = MOLAL(3)                             ! NH4HSO4
         MOLALR(12) = MOLAL(2)                             ! NAHSO4
         MOLALR(07) = MOLAL(5)+MOLAL(6)-MOLAL(3)-MOLAL(2)  ! H2SO4
         MOLALR(07) = MAX(MOLALR(07),ZERO)

      !=================================================================
      ! ----- REVERSE PROBLEMS ----- 
      !
      ! NH4-SO4-NO3 SYSTEM ; SULFATE POOR CASE
      !=================================================================
      ELSE IF (SC.EQ.'N') THEN      
         MOLALR(4) = MOLAL(5) + MOLAL(6)          ! (NH4)2SO4
         AML5      = WAER(3)-2.D0*MOLALR(4)       ! "free" NH4
         MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO) ! NH4NO3 = MIN("free", NO3)

      !=================================================================
      ! NH4-SO4-NO3-NA-CL SYSTEM ; SULFATE POOR, SODIUM POOR CASE
      !=================================================================
      ELSE IF (SC.EQ.'Q') THEN      
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(4) = PSI6                                  ! (NH4)2SO4
         MOLALR(5) = PSI5                                  ! NH4NO3
         MOLALR(6) = PSI4                                  ! NH4CL

      !=================================================================
      ! NH4-SO4-NO3-NA-CL SYSTEM ; SULFATE POOR, SODIUM RICH CASE
      !=================================================================
      ELSE IF (SC.EQ.'R') THEN      
         MOLALR(1) = PSI3                                  ! NACL 
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(3) = PSI2                                  ! NANO3
         MOLALR(4) = ZERO                                  ! (NH4)2SO4
         MOLALR(5) = PSI5                                  ! NH4NO3
         MOLALR(6) = PSI4                                  ! NH4CL

      !=================================================================
      ! UNKNOWN CASE
      !=================================================================
      ELSE
         CALL PUSHERR (1001, ' ') ! FATAL ERROR: CASE NOT SUPPORTED 
      ENDIF

      !=================================================================
      ! CALCULATE WATER CONTENT ; ZSR CORRELATION
      !=================================================================
      WATER = ZERO
      DO I = 1, NPAIR
         WATER = WATER + MOLALR(I)/M0(I)
      ENDDO
      WATER = MAX(WATER, TINY)

      ! Return to calling program
      END SUBROUTINE CALCMR

!------------------------------------------------------------------------------

      SUBROUTINE CALCMDRH( RHI, RHDRY, RHLIQ, DRYCASE, LIQCASE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCMDRH
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE 'DRY' SOLUTION (SUBROUTINE DRYCASE) AND THE
!     'SATURATED LIQUID' SOLUTION (SUBROUTINE LIQCASE).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8  :: RHI, RHDRY, RHLIQ
      EXTERNAL   DRYCASE, LIQCASE

      ! Local variables
      REAL*8  :: WF, ONEMWF
      REAL*8  :: CNH42SO, CNH4HSO, CLCO, CNH4N3O, CNH4CLO, CNA2SO 
      REAL*8  :: CNAHSO, CNANO, CNACLO, GNH3O, GHNO3O, GHCLO   
      REAL*8  :: DAMSUL, DSOSUL, DAMBIS, DSOBIS, DLC, DAMNIT
      REAL*8  :: DAMCHL, DSONIT, DSOCHL, DAMG, DHAG, DNAG
      INTEGER :: I

      !=================================================================
      ! CALCMDRH begins here!
      !=================================================================

      ! Find weight factor
      IF (WFTYP.EQ.0) THEN
         WF = ONE
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
      ENDIF
      ONEMWF  = ONE - WF

      !=================================================================
      ! FIND FIRST SECTION ; DRY ONE
      !=================================================================

      ! Find solution at mdrh by weighting dry & liquid solutions.
      CALL DRYCASE
      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  ! DRY AEROSOL

      CNH42SO = CNH42S4                  ! FIRST (DRY) SOLUTION
      CNH4HSO = CNH4HS4
      CLCO    = CLC 
      CNH4N3O = CNH4NO3
      CNH4CLO = CNH4CL
      CNA2SO  = CNA2SO4
      CNAHSO  = CNAHSO4
      CNANO   = CNANO3
      CNACLO  = CNACL
      GNH3O   = GNH3
      GHNO3O  = GHNO3
      GHCLO   = GHCL

      !=================================================================
      ! FIND SECOND SECTION ; DRY & LIQUID
      !=================================================================
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CLC     = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CNAHSO4 = ZERO
      CNANO3  = ZERO
      CNACL   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
      CALL LIQCASE                   ! SECOND (LIQUID) SOLUTION

      !=================================================================
      ! ADJUST THINGS FOR THE CASE THAT THE 
      ! LIQUID SUB PREDICTS DRY AEROSOL
      !=================================================================
      IF ( WATER .LE. TINY ) THEN

         ! Aqueous phase
         DO I = 1, NIONS
            MOLAL(I)= ZERO           
         ENDDO

         WATER   = ZERO

         ! Solid phase
         CNH42S4 = CNH42SO           
         CNA2SO4 = CNA2SO
         CNAHSO4 = CNAHSO
         CNH4HS4 = CNH4HSO
         CLC     = CLCO
         CNH4NO3 = CNH4N3O
         CNANO3  = CNANO
         CNACL   = CNACLO                                                  
         CNH4CL  = CNH4CLO 

         ! Gas Phase
         GNH3    = GNH3O             
         GHNO3   = GHNO3O
         GHCL    = GHCLO

         GOTO 200
      ENDIF

      !=================================================================
      ! FIND SALT DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
      !=================================================================
      DAMSUL  = CNH42SO - CNH42S4
      DSOSUL  = CNA2SO  - CNA2SO4
      DAMBIS  = CNH4HSO - CNH4HS4
      DSOBIS  = CNAHSO  - CNAHSO4
      DLC     = CLCO    - CLC
      DAMNIT  = CNH4N3O - CNH4NO3
      DAMCHL  = CNH4CLO - CNH4CL
      DSONIT  = CNANO   - CNANO3
      DSOCHL  = CNACLO  - CNACL

      !=================================================================
      ! FIND GAS DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
      !=================================================================
      DAMG    = GNH3O   - GNH3 
      DHAG    = GHCLO   - GHCL
      DNAG    = GHNO3O  - GHNO3

      !=================================================================
      ! FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
      !=================================================================

      ! LIQUID
      MOLAL(1)= ONEMWF*MOLAL(1)                                 ! H+
      MOLAL(2)= ONEMWF*(2.D0*DSOSUL + DSOBIS + DSONIT + DSOCHL) ! NA+
      MOLAL(3)= ONEMWF*(2.D0*DAMSUL + DAMG   + DAMBIS + DAMCHL +
     &                  3.D0*DLC    + DAMNIT )                  ! NH4+
      MOLAL(4)= ONEMWF*(     DAMCHL + DSOCHL + DHAG)            ! CL-
      MOLAL(5)= ONEMWF*(     DAMSUL + DSOSUL + DLC)             ! SO4--
      MOLAL(6)= ONEMWF*(   MOLAL(6) + DSOBIS + DAMBIS + DLC)    ! HSO4-
      MOLAL(7)= ONEMWF*(     DAMNIT + DSONIT + DNAG)            ! NO3-
      WATER   = ONEMWF*WATER

      ! SOLID
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
      CNACL   = WF*CNACLO  + ONEMWF*CNACL
      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL

      ! GAS
      GNH3    = WF*GNH3O   + ONEMWF*GNH3
      GHNO3   = WF*GHNO3O  + ONEMWF*GHNO3
      GHCL    = WF*GHCLO   + ONEMWF*GHCL

      ! Return to calling program
 200  CONTINUE
      END SUBROUTINE CALCMDRH

!------------------------------------------------------------------------------

      SUBROUTINE CALCMDRP( RHI, RHDRY, RHLIQ, DRYCASE, LIQCASE )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCMDRP
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE 'DRY' SOLUTION (SUBROUTINE DRYCASE) AND THE
!     'SATURATED LIQUID' SOLUTION (SUBROUTINE LIQCASE).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8  :: RHI, RHDRY, RHLIQ
      EXTERNAL   DRYCASE, LIQCASE

      ! Local variables
      REAL*8  :: WF, ONEMWF
      REAL*8  :: CNH42SO, CNH4HSO, CLCO, CNH4N3O, CNH4CLO, CNA2SO 
      REAL*8  :: CNAHSO, CNANO, CNACLO
      REAL*8  :: DAMBIS, DSOBIS, DLC, A8
      REAL*8  :: HIEQ, HIEN, A2, A3, A4
      INTEGER :: I

      !=================================================================
      ! CALCMDRP begins here!
      !=================================================================

      ! FIND WEIGHT FACTOR
      IF (WFTYP.EQ.0) THEN
         WF = ONE
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
      ENDIF
      ONEMWF  = ONE - WF

      !=================================================================
      ! FIND FIRST SECTION ; DRY ONE 
      !=================================================================
      CALL DRYCASE
      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  ! DRY AEROSOL

      CNH42SO = CNH42S4              ! FIRST (DRY) SOLUTION
      CNH4HSO = CNH4HS4
      CLCO    = CLC 
      CNH4N3O = CNH4NO3
      CNH4CLO = CNH4CL
      CNA2SO  = CNA2SO4
      CNAHSO  = CNAHSO4
      CNANO   = CNANO3
      CNACLO  = CNACL

      !=================================================================
      ! FIND SECOND SECTION ; DRY & LIQUID 
      !=================================================================
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CLC     = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CNAHSO4 = ZERO
      CNANO3  = ZERO
      CNACL   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
      CALL LIQCASE                   ! SECOND (LIQUID) SOLUTION
 
      !=================================================================
      ! ADJUST THINGS FOR THE CASE THAT THE 
      ! LIQUID SUB PREDICTS DRY AEROSOL
      !=================================================================
      IF ( WATER .LE. TINY ) THEN
         WATER = ZERO

         DO I = 1, NIONS
            MOLAL(I)= ZERO
         ENDDO

         CALL DRYCASE
         GOTO 200
      ENDIF

      !=================================================================
      ! FIND SALT DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
      !=================================================================
      DAMBIS  = CNH4HSO - CNH4HS4
      DSOBIS  = CNAHSO  - CNAHSO4
      DLC     = CLCO    - CLC

      !=================================================================
      ! FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
      !=================================================================

      ! SOLID
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
      CNACL   = WF*CNACLO  + ONEMWF*CNACL
      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL

      ! LIQUID
      WATER   = ONEMWF*WATER

      MOLAL(2)= WAER(1) - 2.D0*CNA2SO4 - CNAHSO4 - CNANO3 -     
     &                         CNACL                            ! NA+
      MOLAL(3)= WAER(3) - 2.D0*CNH42S4 - CNH4HS4 - CNH4CL - 
     &                    3.D0*CLC     - CNH4NO3                ! NH4+
      MOLAL(4)= WAER(5) - CNACL - CNH4CL                        ! CL-
      MOLAL(7)= WAER(4) - CNANO3 - CNH4NO3                      ! NO3-
      MOLAL(6)= ONEMWF*(MOLAL(6) + DSOBIS + DAMBIS + DLC)       ! HSO4-
      MOLAL(5)= WAER(2) - MOLAL(6) - CLC - CNH42S4 - CNA2SO4    ! SO4--

      A8      = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
      IF (MOLAL(5).LE.TINY) THEN
         HIEQ = SQRT(XKW *RH*WATER*WATER)  ! Neutral solution
      ELSE
         HIEQ = A8*MOLAL(6)/MOLAL(5)          
      ENDIF
      HIEN    = MOLAL(4) + MOLAL(7) + MOLAL(6) + 2.D0*MOLAL(5) -
     &          MOLAL(2) - MOLAL(3)
      MOLAL(1)= MAX (HIEQ, HIEN)                                ! H+

      ! GAS (ACTIVITY COEFS FROM LIQUID SOLUTION)
      A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-

      GNH3    = MOLAL(3)/MAX(MOLAL(1),TINY)/A2
      GHNO3   = MOLAL(1)*MOLAL(7)/A3
      GHCL    = MOLAL(1)*MOLAL(4)/A4

      ! Return to calling program
 200  CONTINUE
      END SUBROUTINE CALCMDRP

!------------------------------------------------------------------------------

      SUBROUTINE CALCHS4( HI, SO4I, HSO4I, DELTA )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCHS4
! *** THIS SUBROUTINE CALCULATES THE HSO4 GENERATED FROM (H,SO4).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8 :: HI, SO4I, HSO4I, DELTA

      ! Local variables
      REAL*8 :: A8, BB, CC, DD, SQDD, DELTA1, DELTA2
      !CHARACTER ERRINF*40

      !=================================================================
      ! CALCHS4 begins here!
      !=================================================================

      ! IF TOO LITTLE WATER, DONT SOLVE
      IF (WATER.LE.1d1*TINY) THEN
         DELTA = ZERO 
         RETURN
      ENDIF

      !=================================================================
      ! CALCULATE HSO4 SPECIATION
      !=================================================================
      A8 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.

      BB =-(HI + SO4I + A8)
      CC = HI*SO4I - HSO4I*A8
      DD = BB*BB - 4.D0*CC

      IF ( DD .GE. ZERO ) THEN
         SQDD   = SQRT(DD)
         DELTA1 = 0.5*(-BB + SQDD)
         DELTA2 = 0.5*(-BB - SQDD)
         IF (HSO4I.LE.TINY) THEN
            DELTA = DELTA2
         ELSEIF( HI*SO4I .GE. A8*HSO4I ) THEN
            DELTA = DELTA2
         ELSEIF( HI*SO4I .LT. A8*HSO4I ) THEN
            DELTA = DELTA1
         ELSE
            DELTA = ZERO
         ENDIF
      ELSE
         DELTA  = ZERO
      ENDIF

      !-----------------------------------------------------------------------
      ! Comment out for now
      !!=================================================================
      !!COMPARE DELTA TO TOTAL H+ ; ESTIMATE EFFECT OF HSO4 
      !!=================================================================
      !
      !HYD = MAX(HI, MOLAL(1))
      !IF (HYD.GT.TINY) THEN
      !   IF (DELTA/HYD.GT.0.1d0) THEN
      !      WRITE (ERRINF,'(1PE10.3)') DELTA/HYD*100.0
      !      CALL PUSHERR (0020, ERRINF)
      !   ENDIF
      !ENDIF
      !-----------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE CALCHS4

!------------------------------------------------------------------------------

      SUBROUTINE CALCPH( GG, HI, OHI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCPH
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8 :: GG, HI, OHI

      ! Local variables
      REAL*8 :: AKW, CN, BB, DD, CC

      AKW  = XKW *RH*WATER*WATER
      CN   = SQRT(AKW)

      !=================================================================
      ! CALCPH begins here!
      !=================================================================

      ! GG = (negative charge) - (positive charge)
      ! H+ in excess
      IF ( GG .GT. TINY ) THEN                        
         BB  = -GG
         CC  = -AKW
         DD  =  BB*BB - 4.D0*CC
         HI  =  MAX(0.5D0*(-BB + SQRT(DD)),CN)
         OHI =  AKW/HI

      ! OH- in excess
      ELSE                                        
         BB  =  GG
         CC  = -AKW
         DD  =  BB*BB - 4.D0*CC
         OHI =  MAX(0.5D0*(-BB + SQRT(DD)),CN)
         HI  =  AKW/OHI

      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCPH

!------------------------------------------------------------------------------

      SUBROUTINE CALCACT
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCACT
! *** CALCULATES MULTI-COMPONENET ACTIVITY COEFFICIENTS FROM BROMLEYS
!     METHOD. THE BINARY ACTIVITY COEFFICIENTS ARE CALCULATED BY 
!     KUSIK-MEISNER RELATION (SUBROUTINE KMTAB or SUBROUTINE KMFUL). 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!******************************************************************************
!
      ! Local variables
      REAL, PARAMETER :: URF=0.5
      REAL            :: EX10
      REAL            :: G0(3,4),ZPL,ZMI,AGAMA,SION,H,CH,F1(3),F2(4)
      REAL            :: G
      INTEGER         :: I, J
      REAL*8          :: MPL, XIJ, YJI, ERROU, ERRIN

      !=================================================================
      ! CALCACT begins here!
      !=================================================================
      G(I,J)= (F1(I)/Z(I) + F2(J)/Z(J+3)) / (Z(I)+Z(J+3)) - H

      !=================================================================
      ! SAVE ACTIVITIES IN OLD ARRAY
      !=================================================================
      IF ( FRST ) THEN               

         ! Outer loop
         DO I = 1, NPAIR
            GAMOU(I) = GAMA(I)
         ENDDO
      ENDIF

      ! Inner loop
      DO I = 1, NPAIR              
         GAMIN(I) = GAMA(I)
      ENDDO

      !=================================================================
      ! CALCULATE IONIC ACTIVITY OF SOLUTION
      !=================================================================
      IONIC = 0.0

      DO I = 1, NIONS
         IONIC = IONIC + MOLAL(I)*Z(I)*Z(I)
      ENDDO

      IONIC = MAX(MIN(0.5*IONIC/WATER,20.d0), TINY)

      !=================================================================
      ! CALCULATE BINARY ACTIVITY COEFFICIENTS
      !
      ! G0(1,1)=G11; G0(1,2)=G07; G0(1,3)=G08; G0(1,4)=G10; 
      ! G0(2,1)=G01; G0(2,2)=G02; G0(2,3)=G12; G0(2,4)=G03;
      ! G0(3,1)=G06; G0(3,2)=G04; G0(3,3)=G09; G0(3,4)=G05
      !=================================================================
      IF ( IACALC .EQ. 0 ) THEN              
         ! K.M.; FULL
         CALL KMFUL (IONIC, REAL(TEMP),G0(2,1),G0(2,2),G0(2,4),
     &               G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),
     &               G0(1,4),G0(1,1),G0(2,3))
      ELSE                               
         ! K.M.; TABULATED
         CALL KMTAB (IONIC, REAL(TEMP),G0(2,1),G0(2,2),G0(2,4),
     &               G0(3,2),G0(3,4),G0(3,1),G0(1,2),G0(1,3),G0(3,3),
     &               G0(1,4),G0(1,1),G0(2,3))
      ENDIF

      !=================================================================
      ! CALCULATE MULTICOMPONENT ACTIVITY COEFFICIENTS
      !=================================================================
      AGAMA = 0.511*(298.0/TEMP)**1.5    ! Debye Huckel const. at T
      SION  = SQRT(IONIC)
      H     = AGAMA*SION/(1+SION)

      DO I = 1,3
         F1(I)=0.0
         F2(I)=0.0
      ENDDO
         
      F2(4)=0.0

      DO I = 1,3
         ZPL = Z(I)
         MPL = MOLAL(I)/WATER
         DO J = 1,4
            ZMI   = Z(J+3)
            CH    = 0.25*(ZPL+ZMI)*(ZPL+ZMI)/IONIC
            XIJ   = CH*MPL
            YJI   = CH*MOLAL(J+3)/WATER
            F1(I) = F1(I) + SNGL(YJI*(G0(I,J) + ZPL*ZMI*H))
            F2(J) = F2(J) + SNGL(XIJ*(G0(I,J) + ZPL*ZMI*H))
         ENDDO
      ENDDO

      !=================================================================
      ! LOG10 OF ACTIVITY COEFFICIENTS
      !=================================================================
      GAMA(01) = G(2,1)*ZZ(01)                     ! NACL
      GAMA(02) = G(2,2)*ZZ(02)                     ! NA2SO4
      GAMA(03) = G(2,4)*ZZ(03)                     ! NANO3
      GAMA(04) = G(3,2)*ZZ(04)                     ! (NH4)2SO4
      GAMA(05) = G(3,4)*ZZ(05)                     ! NH4NO3
      GAMA(06) = G(3,1)*ZZ(06)                     ! NH4CL
      GAMA(07) = G(1,2)*ZZ(07)                     ! 2H-SO4
      GAMA(08) = G(1,3)*ZZ(08)                     ! H-HSO4
      GAMA(09) = G(3,3)*ZZ(09)                     ! NH4HSO4
      GAMA(10) = G(1,4)*ZZ(10)                     ! HNO3
      GAMA(11) = G(1,1)*ZZ(11)                     ! HCL
      GAMA(12) = G(2,3)*ZZ(12)                     ! NAHSO4
      GAMA(13) = 0.20*(3.0*GAMA(04)+2.0*GAMA(09))  ! LC ; SCAPE
      !GAMA(13) = 0.50*(GAMA(04)+GAMA(09))          ! LC ; SEQUILIB
      !GAMA(13) = 0.25*(3.0*GAMA(04)+GAMA(07))      ! LC ; AIM

      !=================================================================
      ! CONVERT LOG (GAMA) COEFFICIENTS TO GAMA
      !=================================================================
      DO I = 1, NPAIR
         !GAMA(I)=MAX(-5.0d0, MIN(GAMA(I),5.0d0) )   ! F77 LIBRARY ROUTINE
         !GAMA(I)=10.0**GAMA(I)
         GAMA(I) = EX10(SNGL(GAMA(I)), 5.0)          ! CUTOFF SET TO [-5,5]
         GAMA(I) = GAMIN(I)*(1.0-URF) + URF*GAMA(I)  ! Under-relax GAMA's
      ENDDO

      !=================================================================
      ! SETUP ACTIVITY CALCULATION FLAGS
      !=================================================================

      ! OUTER CALCULATION LOOP ; ONLY IF FRST=.TRUE.
      IF ( FRST ) THEN          

         ! Convergence criterion
         ERROU = ZERO                    

         DO I = 1,NPAIR
            ERROU = MAX(ERROU, ABS((GAMOU(I)-GAMA(I))/GAMOU(I)))
         ENDDO

         ! Setup flags
         CALAOU = ERROU .GE. EPSACT      
         FRST   =.FALSE.
      ENDIF

      ! INNER CALCULATION LOOP ; ALWAYS
      ! Convergence criterion
      ERRIN = ZERO                       

      DO I = 1, NPAIR
         ERRIN = MAX (ERRIN, ABS((GAMIN(I)-GAMA(I))/GAMIN(I)))
      ENDDO

      CALAIN = ERRIN .GE. EPSACT

      ! Increment ACTIVITY call counter
      ICLACT = ICLACT + 1                

      ! Return to calling program
      END SUBROUTINE CALCACT

!------------------------------------------------------------------------------

      SUBROUTINE RSTGAM
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE RSTGAM
! *** RESETS ACTIVITY COEFFICIENT ARRAYS TO DEFAULT VALUE OF 0.1
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      INTEGER :: I

      DO I = 1, NPAIR
         GAMA(I) = 0.1
      ENDDO

      ! Return to calling program
      END SUBROUTINE RSTGAM

!------------------------------------------------------------------------------

      SUBROUTINE KMFUL( IONIC, TEMP, G01, G02, G03, G04, G05, 
     &                  G06,   G07,  G08, G09, G10, G11, G12 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE KMFUL
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments 
      REAL*4 :: IONIC, TEMP
      REAL*4 :: G01, G02, G03, G04, G05, G06
      REAL*4 :: G07, G08, G09, G10, G11, G12

      ! Local variables
      REAL*4 :: Z01, Z02, Z03, Z04, Z05, Z06
      REAL*4 :: Z07, Z08, Z09, Z10, Z11, SION
      REAL*4 :: TI,  CF1, CF2

      !=================================================================
      ! KMFUL begins here!
      !=================================================================
      
      ! Initialize
      Z01  = 1
      Z02  = 2
      Z03  = 1
      Z04  = 2
      Z05  = 1
      Z06  = 1
      Z07  = 2
      Z08  = 1
      Z10  = 1
      Z11  = 1
      SION = SQRT(IONIC)

      !=================================================================
      ! Coefficients at 25 oC
      !=================================================================
      CALL MKBI( 2.230, IONIC, SION, Z01, G01 )
      CALL MKBI( -0.19, IONIC, SION, Z02, G02 )
      CALL MKBI( -0.39, IONIC, SION, Z03, G03 )
      CALL MKBI( -0.25, IONIC, SION, Z04, G04 )
      CALL MKBI( -1.15, IONIC, SION, Z05, G05 )
      CALL MKBI( 0.820, IONIC, SION, Z06, G06 )
      CALL MKBI( -.100, IONIC, SION, Z07, G07 )
      CALL MKBI( 8.000, IONIC, SION, Z08, G08 )
      CALL MKBI( 2.600, IONIC, SION, Z10, G10 )
      CALL MKBI( 6.000, IONIC, SION, Z11, G11 )

      !=================================================================
      ! Correct for T other than 298 K
      !=================================================================
      TI  = TEMP-273.0
      IF (ABS(TI) .GT. 1.0) THEN
         CF1 = 1.125-0.005*TI
         CF2 = (0.125-0.005*TI)*(0.039*IONIC**0.92-0.41*SION/(1.+SION))
         G01 = CF1*G01 - CF2*Z01
         G02 = CF1*G02 - CF2*Z02
         G03 = CF1*G03 - CF2*Z03
         G04 = CF1*G04 - CF2*Z04
         G05 = CF1*G05 - CF2*Z05
         G06 = CF1*G06 - CF2*Z06
         G07 = CF1*G07 - CF2*Z07
         G08 = CF1*G08 - CF2*Z08
         G10 = CF1*G10 - CF2*Z10
         G11 = CF1*G11 - CF2*Z11
      ENDIF

      G09 = G06 + G08 - G11
      G12 = G01 + G08 - G11

      ! Return to calling program
      END SUBROUTINE KMFUL

!------------------------------------------------------------------------------

      SUBROUTINE MKBI( Q, IONIC, SION, ZIP, BI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE MKBI
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!
!******************************************************************************
! 
      ! Arguments
      REAL*4 :: IONIC, Q, SION, ZIP, BI

      ! Local variables
      REAL*4 :: B, C, XX, BI

      !=================================================================
      ! MKBI begins here!
      !=================================================================
      B = .75-.065*Q
      C = 1.0

      IF (IONIC.LT.6.0) C=1.+.055*Q*EXP(-.023*IONIC*IONIC*IONIC)

      XX = -0.5107*SION/(1.+C*SION)
      BI = (1.+B*(1.+.1*IONIC)**Q-B)
      BI = ZIP*ALOG10(BI) + ZIP*XX

      ! Return to calling program
      END SUBROUTINE MKBI

!------------------------------------------------------------------------------

      SUBROUTINE KMTAB( IN,  TEMP, G01, G02, G03, G04, G05,
     &                  G06, G07,  G08, G09, G10, G11, G12 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE KMTAB
! *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!     THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!     LOOKUP TABLES. THE IONIC ACTIVITY 'IONIC' IS INPUT, AND THE ARRAY
!     'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*4, INTENT(IN)  :: IN, Temp
      REAL*4, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*4, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IND

      !=================================================================
      ! KMC_TAB begins here!
      !=================================================================

      ! Find temperature range
      IND = NINT((TEMP-198.0)/25.0) + 1
      IND = MIN(MAX(IND,1),6)

      !=================================================================
      ! Call appropriate routine
      !=================================================================
      IF ( IND .EQ. 1 ) THEN
         CALL KM198(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSEIF ( IND .EQ. 2 ) THEN
         CALL KM223(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSEIF ( IND .EQ. 3 ) THEN
         CALL KM248(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSEIF ( IND .EQ. 4 ) THEN
         CALL KM273(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSEIF ( IND .EQ. 5 ) THEN
         CALL KM298(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ELSE
         CALL KM323(IN,G01,G02,G03,G04,G05,G06,G07,G08,G09,G10,G11,G12)

      ENDIF

      ! Return to calling program
      END SUBROUTINE KMTAB
 
!------------------------------------------------------------------------------

      SUBROUTINE READ_KMC
!
!******************************************************************************
!  Subroutine READ_KMC reads data from binary files in order to initialize
!  parameters for the ISORROPIA routines.  This is necessary since the 
!  original code uses big BLOCK DATA statements which don't compile well on 
!  SGI. (rjp, bdf, bmy, 9/23/02) 
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SETUP" ! DATA_DIR


      ! Local variables
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_KMC begins here!
      !=================================================================

      ! Initialize arrays 
      !CALL INIT_KMC

      ! Read data at 198 K
      FILENAME = TRIM( DATA_DIR ) // 'sulfate_sim_200210/kmc198.bin'
      CALL READ_BINARY( FILENAME, BNC198 )

      ! Read data at 223 K
      FILENAME = TRIM( DATA_DIR ) // 'sulfate_sim_200210/kmc223.bin'
      CALL READ_BINARY( FILENAME, BNC223 ) 

      ! Read data at 248 K
      FILENAME = TRIM( DATA_DIR ) // 'sulfate_sim_200210/kmc248.bin'
      CALL READ_BINARY( FILENAME, BNC248 )

      ! Read data at 273 K
      FILENAME = TRIM( DATA_DIR ) // 'sulfate_sim_200210/kmc273.bin'
      CALL READ_BINARY( FILENAME, BNC273 )

      ! Read data at 298 K
      FILENAME = TRIM( DATA_DIR ) // 'sulfate_sim_200210/kmc298.bin'
      CALL READ_BINARY( FILENAME, BNC298 )

      ! Read data at 323 K
      FILENAME = TRIM( DATA_DIR ) // 'sulfate_sim_200210/kmc323.bin'
      CALL READ_BINARY( FILENAME, BNC323 )

      ! Return to calling program
      END SUBROUTINE READ_KMC

!------------------------------------------------------------------------------
      
      SUBROUTINE READ_BINARY( FILENAME, ARRAY )
!
!******************************************************************************
!  Subroutine READ_BINARY reads a binary file from disk. (rjp, bmy, 9/23/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHAR*255) : Name of the binary file to be read

!  Arguments as Output:
!  ============================================================================
!  (2 ) ARRAY    (REAL*4  ) : Data array
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_FILE, IOERROR
      
      ! Arguments
      CHARACTER(LEN=255), INTENT(IN)  :: FILENAME
      REAL*4,             INTENT(OUT) :: ARRAY(IMAX,JMAX)

      ! Local variables
      INTEGER                         :: I, J, IOS
      
      !=================================================================
      ! READ_BINARY begins here!
      !=================================================================

      ! Open the file
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD',
     &               FORM='UNFORMATTED',    IOSTAT=IOS )

      ! Error check
      IF ( IU_FILE /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_kmc:1' )
      
      ! Read data
      DO J = 1, JMAX
         READ( IU_FILE, IOSTAT=IOS ) ( ARRAY(I,J), I=1,IMAX )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_kmc:2' )       
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_BINARY

!-----------------------------------------------------------------------------

      SUBROUTINE KM198( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!*****************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM198
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
! 
!      TEMPERATURE IS 198K
! 
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!*****************************************************************************
!
      ! Arguments
      REAL*4, INTENT(IN)  :: IN
      REAL*4, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*4, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IPOS

      !=================================================================
      ! KMC323 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.300000E+02) THEN
         IPOS = NINT( 0.200000E+02*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC198(IPOS,1 )
      G02 = BNC198(IPOS,2 )
      G03 = BNC198(IPOS,3 )
      G04 = BNC198(IPOS,4 )
      G05 = BNC198(IPOS,5 )
      G06 = BNC198(IPOS,6 )
      G07 = BNC198(IPOS,7 )
      G08 = BNC198(IPOS,8 )
      G09 = BNC198(IPOS,9 )
      G10 = BNC198(IPOS,10)
      G11 = BNC198(IPOS,11)
      G12 = BNC198(IPOS,12)

      ! Return to calling program
      END SUBROUTINE KM198

!------------------------------------------------------------------------------

      SUBROUTINE KM223( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!******************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM223
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
! 
!      TEMPERATURE IS 223K
! 
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
! 
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!******************************************************************************
!
      ! Arguments
      REAL*4, INTENT(IN)  :: IN
      REAL*4, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*4, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variable
      INTEGER             :: IPOS

      !=================================================================
      ! KMC223 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.300000E+02) THEN
         IPOS = NINT( 0.200000E+02*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC223(IPOS,1)
      G02 = BNC223(IPOS,2)
      G03 = BNC223(IPOS,3)
      G04 = BNC223(IPOS,4)
      G05 = BNC223(IPOS,5)
      G06 = BNC223(IPOS,6)
      G07 = BNC223(IPOS,7)
      G08 = BNC223(IPOS,8)
      G09 = BNC223(IPOS,9)
      G10 = BNC223(IPOS,10)
      G11 = BNC223(IPOS,11)
      G12 = BNC223(IPOS,12)

      ! Return point ; End of subroutine
      END SUBROUTINE KM223

!------------------------------------------------------------------------------

      SUBROUTINE KM248( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!******************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM248
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!  
!      TEMPERATURE IS 248K
!  
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!******************************************************************************
!
      REAL*4, INTENT(IN)  :: IN
      REAL*4, INTENT(OUT) :: G01, G02, G03, G04, G05, G06 
      REAL*4, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER :: IPOS

      !=================================================================
      ! KMC248 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.300000E+02) THEN
         IPOS = NINT( 0.200000E+02*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC248(IPOS,1)
      G02 = BNC248(IPOS,2)
      G03 = BNC248(IPOS,3)
      G04 = BNC248(IPOS,4)
      G05 = BNC248(IPOS,5)
      G06 = BNC248(IPOS,6)
      G07 = BNC248(IPOS,7)
      G08 = BNC248(IPOS,8)
      G09 = BNC248(IPOS,9)
      G10 = BNC248(IPOS,10)
      G11 = BNC248(IPOS,11)
      G12 = BNC248(IPOS,12)

      ! Return point ; End of subroutine
      END SUBROUTINE KM248
 
!------------------------------------------------------------------------------

      SUBROUTINE KM273( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!*****************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM273
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!  
!      TEMPERATURE IS 273K
!
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!*****************************************************************************
!
      ! Arguments
      REAL*4, INTENT(IN)  :: IN
      REAL*4, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*4, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IPOS

      !=================================================================
      ! KMC273 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.300000E+02) THEN
         IPOS = NINT( 0.200000E+02*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC273(IPOS,1 )
      G02 = BNC273(IPOS,2 )
      G03 = BNC273(IPOS,3 )
      G04 = BNC273(IPOS,4 )
      G05 = BNC273(IPOS,5 ) 
      G06 = BNC273(IPOS,6 )
      G07 = BNC273(IPOS,7 )
      G08 = BNC273(IPOS,8 )
      G09 = BNC273(IPOS,9 )
      G10 = BNC273(IPOS,10)
      G11 = BNC273(IPOS,11)
      G12 = BNC273(IPOS,12)
      
      ! Return to calling program
      END SUBROUTINE KM273

!------------------------------------------------------------------------------

      SUBROUTINE KM298( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!******************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM298
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!
!      TEMPERATURE IS 298K
!
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!******************************************************************************
!
      ! Arguments
      REAL*4, INTENT(IN)  :: IN
      REAL*4, INTENT(OUT) :: G01, G02, G03, G04, G05, G06 
      REAL*4, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IPOS

      !=================================================================
      ! KMC298 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.300000E+02) THEN
         IPOS = NINT( 0.200000E+02*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC298(IPOS,1 )
      G02 = BNC298(IPOS,2 )
      G03 = BNC298(IPOS,3 )
      G04 = BNC298(IPOS,4 )
      G05 = BNC298(IPOS,5 )
      G06 = BNC298(IPOS,6 )
      G07 = BNC298(IPOS,7 )
      G08 = BNC298(IPOS,8 )
      G09 = BNC298(IPOS,9 )
      G10 = BNC298(IPOS,10)
      G11 = BNC298(IPOS,11)
      G12 = BNC298(IPOS,12)

      ! Return to calling program
      END SUBROUTINE KM298
 
!------------------------------------------------------------------------------

      SUBROUTINE KM323( IN, G01, G02, G03, G04, G05, G06,
     &                      G07, G08, G09, G10, G11, G12 )
!
!******************************************************************************
!  *** ISORROPIA CODE
!  *** SUBROUTINE KM323
!  *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
!      THE COMPUTATIONS HAVE BEEN PERFORMED AND THE RESULTS ARE STORED IN
!      LOOKUP TABLES. THE IONIC ACTIVITY 'IN' IS INPUT, AND THE ARRAY
!      'BINARR' IS RETURNED WITH THE BINARY COEFFICIENTS. 
!  
!      TEMPERATURE IS 323K
!  
!  *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!  *** WRITTEN BY ATHANASIOS NENES
!  
!  NOTES:
!  (1 ) Adapted for GEOS-CHEM (rjp, bdf, bmy, 9/28/02)
!******************************************************************************
!
      REAL*4, INTENT(IN)  :: IN
      REAL*4, INTENT(OUT) :: G01, G02, G03, G04, G05, G06
      REAL*4, INTENT(OUT) :: G07, G08, G09, G10, G11, G12

      ! Local variables
      INTEGER             :: IPOS

      !=================================================================
      ! KMC323 begins here!
      !=================================================================

      ! Find position in arrays for bincoef
      IF (IN.LE. 0.300000E+02) THEN
         ipos = NINT( 0.200000E+02*IN) + 1
      ELSE
         IPOS = 600+NINT( 0.200000E+01*IN- 0.300000E+02)
      ENDIF

      IPOS = MIN( IPOS, IMAX )

      ! Assign values to return array
      G01 = BNC323(IPOS,1 )
      G02 = BNC323(IPOS,2 )
      G03 = BNC323(IPOS,3 )
      G04 = BNC323(IPOS,4 )
      G05 = BNC323(IPOS,5 )
      G06 = BNC323(IPOS,6 )
      G07 = BNC323(IPOS,7 )
      G08 = BNC323(IPOS,8 ) 
      G09 = BNC323(IPOS,9 )
      G10 = BNC323(IPOS,10)
      G11 = BNC323(IPOS,11)
      G12 = BNC323(IPOS,12)

      ! Return to calling program
      END SUBROUTINE KM323

!------------------------------------------------------------------------------

      SUBROUTINE INIT_KMC
!
!******************************************************************************
!  Subroutine CLEANUP_KMC initializes and zeroes all module arrays. 
!  (rjp, bmy, 9/23/02)
!******************************************************************************
!


      ! Return to calling program
      END SUBROUTINE INIT_KMC

!------------------------------------------------------------------------------

      SUBROUTINE CHRBLN( STR, IBLK )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE CHRBLN
!  Purpose        : Position of last non-blank character in a string
!  Author         : Athanasios Nenes
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  STR        is the CHARACTER variable containing the string examined
!  IBLK       is a INTEGER variable containing the position of last non 
!             blank character. If string is all spaces (ie '   '), then
!             the value returned is 1.
!
!  EXAMPLE:
!             STR = 'TEST1.DAT     '
!             CALL CHRBLN (STR, IBLK)
!          
!  after execution of this code segment, "IBLK" has the value "9", which
!  is the position of the last non-blank character of "STR".
!
!***********************************************************************
!
      ! Arguments
      CHARACTER(LEN=*) :: STR
      INTEGER          :: IBLK

      ! Local variables
      INTEGER          :: I, ILEN
      
      !=================================================================
      ! CHRBLN begins here!
      !=================================================================
      IBLK = 1                       ! Substring pointer (default=1)
      ILEN = LEN(STR)                ! Length of string 
      DO I = ILEN, 1, -1              
         IF (STR(i:i).NE.' ' .AND. STR(i:i).NE.CHAR(0)) THEN
            IBLK = i
            RETURN
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE CHRBLN

!------------------------------------------------------------------------------

      SUBROUTINE SHFTRGHT( CHR )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE SHFTRGHT
!  Purpose        : RIGHT-JUSTIFICATION FUNCTION ON A STRING
!  Author         : Athanasios Nenes
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  STRING     is the CHARACTER variable with the string to be justified
!
!  EXAMPLE:
!             STRING    = 'AAAA    '
!             CALL SHFTRGHT (STRING)
!          
!  after execution of this code segment, STRING contains the value
!  '    AAAA'.
!
!*************************************************************************
!
      ! Arguments
      CHARACTER(LEN=*) :: CHR

      ! Local variables
      INTEGER          :: I, I1, I2
      
      !=================================================================
      ! SHFTRGHT begins here!
      !=================================================================
      I1  = LEN(CHR)             ! Total length of string
      CALL CHRBLN(CHR,I2)        ! Position of last non-blank character
      IF (I2.EQ.I1) RETURN

      DO I = I2, 1, -1            ! Shift characters
         CHR(I1+I-I2:I1+I-I2) = CHR(I:I)
         CHR(I:I) = ' '
      ENDDO

      ! Return to calling program
      END SUBROUTINE SHFTRGHT

!------------------------------------------------------------------------------

      SUBROUTINE RPLSTR( STRING, OLD, NEW, IERR )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE RPLSTR
!  Purpose        : REPLACE CHARACTERS OCCURING IN A STRING
!  Author         : Athanasios Nenes
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  STRING     is the CHARACTER variable with the string to be edited
!  OLD        is the old character which is to be replaced
!  NEW        is the new character which OLD is to be replaced with
!  IERR       is 0 if everything went well, is 1 if 'NEW' contains 'OLD'.
!             In this case, this is invalid, and no change is done.
!
!  EXAMPLE:
!             STRING    = 'AAAA'
!             OLD       = 'A'
!             NEW       = 'B' 
!             CALL RPLSTR (STRING, OLD, NEW)
!          
!  after execution of this code segment, STRING contains the value
!  'BBBB'.
!
!*************************************************************************
!
      ! Arguments
      CHARACTER(LEN=*) :: STRING, OLD, NEW
      INTEGER          :: IERR

      ! Local variables
      INTEGER          :: ILO, IP

      !=================================================================
      ! RPLSTR begins here!
      !=================================================================
      ILO = LEN(OLD)

      ! CHECK AND SEE IF 'NEW' CONTAINS 'OLD', WHICH CANNOT
      IP = INDEX(NEW,OLD)
      IF (IP.NE.0) THEN
         IERR = 1
         RETURN
      ELSE
         IERR = 0
      ENDIF

      ! PROCEED WITH REPLACING
10    IP = INDEX(STRING,OLD)      ! SEE IF 'OLD' EXISTS IN 'STRING'
      IF (IP.EQ.0) RETURN         ! 'OLD' DOES NOT EXIST ; RETURN
      STRING(IP:IP+ILO-1) = NEW   ! REPLACE SUBSTRING 'OLD' WITH 'NEW'
      GOTO 10                     ! GO FOR NEW OCCURANCE OF 'OLD'

      ! Return to calling program
      END SUBROUTINE RPLSTR
        
!------------------------------------------------------------------------------

      SUBROUTINE INPTD( VAR, DEF, PROMPT, PRFMT, IERR )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE INPTD
!  Purpose        : Prompts user for a value (DOUBLE). A default value
!                   is provided, so if user presses <Enter>, the default
!                   is used. 
!  Author         : Athanasios Nenes
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  VAR        is the DOUBLE PRECISION variable which value is to be saved 
!  DEF        is a DOUBLE PRECISION variable, with the default value of VAR. !       
!  PROMPT     is a CHARACTER varible containing the prompt string.     
!  PRFMT      is a CHARACTER variable containing the FORMAT specifier
!             for the default value DEF.
!  IERR       is an INTEGER error flag, and has the values:
!             0 - No error detected.
!             1 - Invalid FORMAT and/or Invalid default value.
!             2 - Bad value specified by user
!
!  EXAMPLE:
!             CALL INPTD (VAR, 1.0D0, 'Give value for A ', '*', Ierr)
!          
!  after execution of this code segment, the user is prompted for the
!  value of variable VAR. If <Enter> is pressed (ie no value is specified)
!  then 1.0 is assigned to VAR. The default value is displayed in free-
!  format. The error status is specified by variable Ierr
!
!***********************************************************************
!
      ! Arguments
      CHARACTER(LEN=*)   :: PROMPT, PRFMT
      CHARACTER(LEN=128) :: BUFFER

      ! Local variables
      REAL*8             :: DEF,  VAR
      INTEGER            :: IERR, IEND

      !=================================================================
      ! INPTD begins here!
      !=================================================================
      IERR = 0

      ! Write default value to work buffer
      WRITE (BUFFER, FMT=PRFMT, ERR=10) DEF
      CALL CHRBLN (BUFFER, IEND)

      ! Prompt user for input and read it
      WRITE (*,*) PROMPT,' [',BUFFER(1:IEND),']: '
      READ  (*, '(A)', ERR=20, END=20) BUFFER
      CALL CHRBLN (BUFFER,IEND)

      ! Read data or set default ?
      IF ( IEND .EQ. 1 .AND. BUFFER(1:1) .EQ. ' ' ) THEN
         VAR = DEF
      ELSE
         READ( BUFFER, *, ERR=20, END=20 ) VAR
      ENDIF

      ! Return
30    RETURN

      !=================================================================`
      ! ERROR HANDLER
      !=================================================================
10    IERR = 1       ! Bad FORMAT and/or bad default value
      GOTO 30

20    IERR = 2       ! Bad number given by user
      GOTO 30

      END SUBROUTINE INPTD

!------------------------------------------------------------------------------

      SUBROUTINE PUSHEND( IUNIT )
!
!*************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE Pushend 
!  Purpose        : Positions the pointer of a sequential file at its end
!                   Simulates the ACCESS='APPEND' clause of a F77L OPEN
!                   statement with Standard Fortran commands.
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  Iunit      is a INTEGER variable, the file unit which the file is 
!             connected to.
!
!  EXAMPLE:
!             CALL PUSHEND (10)
!          
!  after execution of this code segment, the pointer of unit 10 is 
!  pushed to its end.
!
!***********************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: IUNIT

      ! Local variablese
      LOGICAL             :: OPNED

      !=================================================================
      ! PUSHEND begins here!
      !=================================================================

      ! INQUIRE IF Iunit CONNECTED TO FIL
      INQUIRE (UNIT=Iunit, OPENED=OPNED)
      IF (.NOT.OPNED) GOTO 25

      ! Iunit CONNECTED, PUSH POINTER TO END
10    READ (Iunit,'()', ERR=20, END=20)
      GOTO 10

      ! Return point
20    BACKSPACE( IUNIT )
25    RETURN

      ! Return to calling program
      END SUBROUTINE PUSHEND

!------------------------------------------------------------------------------

      SUBROUTINE APPENDEXT( FILENAME, DEFEXT, OVERWRITE )
!
!******************************************************************************
!
!  TOOLBOX LIBRARY v.1.0 (May 1995)
!
!  Program unit   : SUBROUTINE APPENDEXT
!  Purpose        : Fix extension in file name string
!
!  ======================= ARGUMENTS / USAGE =============================
!
!  Filename   is the CHARACTER variable with the file name
!  Defext     is the CHARACTER variable with extension (including '.',
!             ex. '.DAT')
!  Overwrite  is a LOGICAL value, .TRUE. overwrites any existing extension
!             in "Filename" with "Defext", .FALSE. puts "Defext" only if 
!             there is no extension in "Filename".
!
!  EXAMPLE:
!             FILENAME1 = 'TEST.DAT'
!             FILENAME2 = 'TEST.DAT'
!             CALL APPENDEXT (FILENAME1, '.TXT', .FALSE.)
!             CALL APPENDEXT (FILENAME2, '.TXT', .TRUE. )
!          
!  after execution of this code segment, "FILENAME1" has the value 
!  'TEST.DAT', while "FILENAME2" has the value 'TEST.TXT'
!
!******************************************************************************
!
      ! Arguments 
      CHARACTER(LEN=*) :: Filename, Defext
      LOGICAL          :: Overwrite

      ! Local variables
      INTEGER          :: IDOT, IEND

      !=================================================================
      ! APPENDEXT begins here!
      !=================================================================
      CALL CHRBLN (Filename, Iend)

      ! Filename empty
      IF (Filename(1:1).EQ.' ' .AND. Iend.EQ.1) RETURN  

      ! Append extension ?
      Idot = INDEX (Filename, '.')                      
      IF (Idot.EQ.0) Filename = Filename(1:Iend)//Defext
      IF (Overwrite .AND. Idot.NE.0)
     &              Filename = Filename(:Idot-1)//Defext

      ! Return to calling program
      END SUBROUTINE APPENDEXT

!------------------------------------------------------------------------------

      SUBROUTINE POLY3( A1, A2, A3, ROOT, ISLV )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE POLY3
! *** FINDS THE REAL ROOTS OF THE THIRD ORDER ALGEBRAIC EQUATION:
!     X**3 + A1*X**2 + A2*X + A3 = 0.0
!     THE EQUATION IS SOLVED ANALYTICALLY.
!
!     PARAMETERS A1, A2, A3 ARE SPECIFIED BY THE USER. THE MINIMUM
!     NONEGATIVE ROOT IS RETURNED IN VARIABLE 'ROOT'. IF NO ROOT IS 
!     FOUND (WHICH IS GREATER THAN ZERO), ROOT HAS THE VALUE 1D30.
!     AND THE FLAG ISLV HAS A VALUE GREATER THAN ZERO.
!
!     SOLUTION FORMULA IS FOUND IN PAGE 32 OF:
!     MATHEMATICAL HANDBOOK OF FORMULAS AND TABLES
!     SCHAUM'S OUTLINE SERIES
!     MURRAY SPIEGER, McGRAW-HILL, NEW YORK, 1968
!     (GREEK TRANSLATION: BY SOTIRIOS PERSIDES, ESPI, ATHENS, 1976)
!
!     A SPECIAL CASE IS CONSIDERED SEPERATELY ; WHEN A3 = 0, THEN
!     ONE ROOT IS X=0.0, AND THE OTHER TWO FROM THE SOLUTION OF THE
!     QUADRATIC EQUATION X**2 + A1*X + A2 = 0.0
!     THIS SPECIAL CASE IS CONSIDERED BECAUSE THE ANALYTICAL FORMULA 
!     DOES NOT YIELD ACCURATE RESULTS (DUE TO NUMERICAL ROUNDOFF ERRORS)
!
! *** COPYRIGHT 1996-98, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8            :: A1, A2, A3, ROOT
      INTEGER           :: ISLV

      ! Local variables
      REAL*8, PARAMETER :: EXPON = 1.D0 / 3.D0
      REAL*8, PARAMETER :: ZERO  = 0.D0
      REAL*8, PARAMETER :: THET1 = 120.D0 / 180.D0 
      REAL*8, PARAMETER :: THET2 = 240.D0 / 180.D0
      REAL*8, PARAMETER :: PI    = 3.14159265358932d0
      REAL*8, PARAMETER :: EPS   = 1D-50
      REAL*8            :: X(3), D, Q, THET, COEF, SSIG, S, TSIG, T, SQD
      INTEGER           :: I, IX
      

      !=================================================================
      ! POLY3 begins here!
      ! 
      ! SPECIAL CASE : QUADRATIC*X EQUATION 
      !=================================================================
      IF (ABS(A3).LE.EPS) THEN 
         ISLV = 1
         IX   = 1
         X(1) = ZERO
         D    = A1*A1-4.D0*A2
         IF (D.GE.ZERO) THEN
            IX   = 3
            SQD  = SQRT(D)
            X(2) = 0.5*(-A1+SQD)
            X(3) = 0.5*(-A1-SQD)
         ENDIF
      ELSE

      !=================================================================
      ! NORMAL CASE : CUBIC EQUATION
      !=================================================================

         ! DEFINE PARAMETERS Q, R, S, T, D 
         ISLV= 1
         Q   = (3.D0*A2 - A1*A1)/9.D0
         R   = (9.D0*A1*A2 - 27.D0*A3 - 2.D0*A1*A1*A1)/54.D0
         D   = Q*Q*Q + R*R
         
         !==============================================================
         ! CALCULATE ROOTS
         !==============================================================

         ! D < 0, THREE REAL ROOTS
         IF (D.LT.-EPS) THEN        ! D < -EPS  : D < ZERO
            IX   = 3
            THET = EXPON*ACOS(R/SQRT(-Q*Q*Q))
            COEF = 2.D0*SQRT(-Q)
            X(1) = COEF*COS(THET)            - EXPON*A1
            X(2) = COEF*COS(THET + THET1*PI) - EXPON*A1
            X(3) = COEF*COS(THET + THET2*PI) - EXPON*A1

         ! D = 0, THREE REAL (ONE DOUBLE) ROOTS
         ELSE IF (D.LE.EPS) THEN    ! -EPS <= D <= EPS  : D = ZERO
            IX   = 2
            SSIG = SIGN (1.D0, R)
            S    = SSIG*(ABS(R))**EXPON
            X(1) = 2.D0*S  - EXPON*A1
            X(2) =     -S  - EXPON*A1

         ! D > 0, ONE REAL ROOT
         ELSE                       ! D > EPS  : D > ZERO
            IX   = 1
            SQD  = SQRT(D)
            SSIG = SIGN (1.D0, R+SQD)       ! TRANSFER SIGN TO SSIG
            TSIG = SIGN (1.D0, R-SQD)
            S    = SSIG*(ABS(R+SQD))**EXPON ! EXPONENTIATE ABS() 
            T    = TSIG*(ABS(R-SQD))**EXPON
            X(1) = S + T - EXPON*A1
         ENDIF
      ENDIF

      !=================================================================
      ! SELECT APPROPRIATE ROOT
      !=================================================================
      ROOT = 1.D30
      DO I = 1, IX
         IF (X(I).GT.ZERO) THEN
            ROOT = MIN (ROOT, X(I))
            ISLV = 0
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE POLY3

!------------------------------------------------------------------------------

      SUBROUTINE POLY3B( A1, A2, A3, RTLW, RTHI, ROOT, ISLV )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE POLY3B
! *** FINDS A REAL ROOT OF THE THIRD ORDER ALGEBRAIC EQUATION:
!     X**3 + A1*X**2 + A2*X + A3 = 0.0
!     THE EQUATION IS SOLVED NUMERICALLY (BISECTION).
!
!     PARAMETERS A1, A2, A3 ARE SPECIFIED BY THE USER. THE MINIMUM
!     NONEGATIVE ROOT IS RETURNED IN VARIABLE 'ROOT'. IF NO ROOT IS 
!     FOUND (WHICH IS GREATER THAN ZERO), ROOT HAS THE VALUE 1D30.
!     AND THE FLAG ISLV HAS A VALUE GREATER THAN ZERO.
!
!     RTLW, RTHI DEFINE THE INTERVAL WHICH THE ROOT IS LOOKED FOR.
!
!******************************************************************************
!
      ! Arguments
      INTEGER            :: ISLV
      REAL*8             :: A1, A2, A3, RTLW, RTHI, ROOT

      ! Local variables
      INTEGER, PARAMETER :: MAXIT=100
      INTEGER, PARAMETER :: NDIV=5
      INTEGER            :: ISLV, I, K1, K2
      REAL*8,  PARAMETER :: ZERO=0.D0
      REAL*8,  PARAMETER :: EPS=1D-15
      REAL*8             :: EX10, X1, Y1, DX, X2, Y2, Y3, X3

      ! Statement function
      !FUNC(X) = X**3.d0 + A1*X**2.0 + A2*X + A3

      !=================================================================
      ! POLY3B begins here!
      !=================================================================
       
      ! Initial values for bisection
      X1   = RTLW
      Y1   = FUNC(X1)

      ! Is low a root?
      IF ( ABS(Y1) .LE. EPS ) THEN     
         ROOT = RTLW
         GOTO 50
      ENDIF

      !=================================================================
      ! Root tracking ; for the range of HI and LO
      !=================================================================
      DX = (RTHI-RTLW)/FLOAT(NDIV)

      DO I = 1, NDIV
         X2 = X1+DX
         Y2 = FUNC (X2)

         ! (Y1*Y2.LT.ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y2) .LT. ZERO ) GOTO 20 
         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================
      ! No subdivision with solution found
      !=================================================================
      IF ( ABS( Y2) .LT. EPS ) THEN   ! X2 is a root
         ROOT = X2
      ELSE
         ROOT = 1.d30
         ISLV = 1
      ENDIF
      GOTO 50

      !=================================================================
      ! BISECTION
      !=================================================================
20    CONTINUE

      DO I = 1, MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNC (X3)

         ! (Y1*Y3 .LE. ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y3) .LE. ZERO ) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS( X2 - X1 ) .LE. EPS * X1 ) GOTO 40
      ENDDO

      !=================================================================
      ! Converged ; RETURN
      !=================================================================
40    CONTINUE
      X3   = 0.5*(X1+X2)
      Y3   = FUNC (X3)
      ROOT = X3
      ISLV = 0

      ! Return to calling program
50    CONTINUE

      CONTAINS

      !-----------------------------------------------------------------------

      FUNCTION FUNC( X ) RESULT ( VALUE )

      REAL*8 :: X, VALUE

      ! This replaces the old statement function
      VALUE = X**3.d0 + A1*X**2.0 + A2*X + A3

      ! Return to POLY3B
      END FUNCTION FUNC

      !-----------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE POLY3B
      
!------------------------------------------------------------------------------

      FUNCTION EX10( X, K )
! 
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION EX10
! *** 10^X FUNCTION ; ALTERNATE OF LIBRARY ROUTINE ; USED BECAUSE IT IS
!     MUCH FASTER BUT WITHOUT GREAT LOSS IN ACCURACY. , 
!     MAXIMUM ERROR IS 2%, EXECUTION TIME IS 42% OF THE LIBRARY ROUTINE 
!     (ON A 80286/80287 MACHINE, using Lahey FORTRAN 77 v.3.0).
!
!     EXPONENT RANGE IS BETWEEN -K AND K (K IS THE REAL ARGUMENT 'K')
!     MAX VALUE FOR K: 9.999
!     IF X < -K, X IS SET TO -K, IF X > K, X IS SET TO K
!
!     THE EXPONENT IS CALCULATED BY THE PRODUCT ADEC*AINT, WHERE ADEC
!     IS THE MANTISSA AND AINT IS THE MAGNITUDE (EXPONENT). BOTH 
!     MANTISSA AND MAGNITUDE ARE PRE-CALCULATED AND STORED IN LOOKUP
!     TABLES ; THIS LEADS TO THE INCREASED SPEED.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*4, INTENT(IN) :: X, K

      ! Local variables
      REAL*4             :: EX10, Y
      INTEGER            :: K1, K2

      ! For integer part
      REAL*4, SAVE       :: AINT10(20) = (/
     & 0.1000E-08, 0.1000E-07, 0.1000E-06, 0.1000E-05, 0.1000E-04,
     & 0.1000E-03, 0.1000E-02, 0.1000E-01, 0.1000E+00, 0.1000E+01,
     & 0.1000E+02, 0.1000E+03, 0.1000E+04, 0.1000E+05, 0.1000E+06,
     & 0.1000E+07, 0.1000E+08, 0.1000E+09, 0.1000E+10, 0.1000E+11/)

      ! For decimal part
      REAL*4, SAVE       :: ADEC10(200) = (/
     & 0.1023E+00, 0.1047E+00, 0.1072E+00, 0.1096E+00, 0.1122E+00,
     & 0.1148E+00, 0.1175E+00, 0.1202E+00, 0.1230E+00, 0.1259E+00,
     & 0.1288E+00, 0.1318E+00, 0.1349E+00, 0.1380E+00, 0.1413E+00,
     & 0.1445E+00, 0.1479E+00, 0.1514E+00, 0.1549E+00, 0.1585E+00,
     & 0.1622E+00, 0.1660E+00, 0.1698E+00, 0.1738E+00, 0.1778E+00,
     & 0.1820E+00, 0.1862E+00, 0.1905E+00, 0.1950E+00, 0.1995E+00,
     & 0.2042E+00, 0.2089E+00, 0.2138E+00, 0.2188E+00, 0.2239E+00,
     & 0.2291E+00, 0.2344E+00, 0.2399E+00, 0.2455E+00, 0.2512E+00,
     & 0.2570E+00, 0.2630E+00, 0.2692E+00, 0.2754E+00, 0.2818E+00,
     & 0.2884E+00, 0.2951E+00, 0.3020E+00, 0.3090E+00, 0.3162E+00,
     & 0.3236E+00, 0.3311E+00, 0.3388E+00, 0.3467E+00, 0.3548E+00,
     & 0.3631E+00, 0.3715E+00, 0.3802E+00, 0.3890E+00, 0.3981E+00,
     & 0.4074E+00, 0.4169E+00, 0.4266E+00, 0.4365E+00, 0.4467E+00,
     & 0.4571E+00, 0.4677E+00, 0.4786E+00, 0.4898E+00, 0.5012E+00,
     & 0.5129E+00, 0.5248E+00, 0.5370E+00, 0.5495E+00, 0.5623E+00,
     & 0.5754E+00, 0.5888E+00, 0.6026E+00, 0.6166E+00, 0.6310E+00,
     & 0.6457E+00, 0.6607E+00, 0.6761E+00, 0.6918E+00, 0.7079E+00,
     & 0.7244E+00, 0.7413E+00, 0.7586E+00, 0.7762E+00, 0.7943E+00,
     & 0.8128E+00, 0.8318E+00, 0.8511E+00, 0.8710E+00, 0.8913E+00,
     & 0.9120E+00, 0.9333E+00, 0.9550E+00, 0.9772E+00, 0.1000E+01,
     & 0.1023E+01, 0.1047E+01, 0.1072E+01, 0.1096E+01, 0.1122E+01,
     & 0.1148E+01, 0.1175E+01, 0.1202E+01, 0.1230E+01, 0.1259E+01,
     & 0.1288E+01, 0.1318E+01, 0.1349E+01, 0.1380E+01, 0.1413E+01,
     & 0.1445E+01, 0.1479E+01, 0.1514E+01, 0.1549E+01, 0.1585E+01,
     & 0.1622E+01, 0.1660E+01, 0.1698E+01, 0.1738E+01, 0.1778E+01,
     & 0.1820E+01, 0.1862E+01, 0.1905E+01, 0.1950E+01, 0.1995E+01,
     & 0.2042E+01, 0.2089E+01, 0.2138E+01, 0.2188E+01, 0.2239E+01,
     & 0.2291E+01, 0.2344E+01, 0.2399E+01, 0.2455E+01, 0.2512E+01,
     & 0.2570E+01, 0.2630E+01, 0.2692E+01, 0.2754E+01, 0.2818E+01,
     & 0.2884E+01, 0.2951E+01, 0.3020E+01, 0.3090E+01, 0.3162E+01,
     & 0.3236E+01, 0.3311E+01, 0.3388E+01, 0.3467E+01, 0.3548E+01,
     & 0.3631E+01, 0.3715E+01, 0.3802E+01, 0.3890E+01, 0.3981E+01,
     & 0.4074E+01, 0.4169E+01, 0.4266E+01, 0.4365E+01, 0.4467E+01,
     & 0.4571E+01, 0.4677E+01, 0.4786E+01, 0.4898E+01, 0.5012E+01,
     & 0.5129E+01, 0.5248E+01, 0.5370E+01, 0.5495E+01, 0.5623E+01,
     & 0.5754E+01, 0.5888E+01, 0.6026E+01, 0.6166E+01, 0.6310E+01,
     & 0.6457E+01, 0.6607E+01, 0.6761E+01, 0.6918E+01, 0.7079E+01,
     & 0.7244E+01, 0.7413E+01, 0.7586E+01, 0.7762E+01, 0.7943E+01,
     & 0.8128E+01, 0.8318E+01, 0.8511E+01, 0.8710E+01, 0.8913E+01,
     & 0.9120E+01, 0.9333E+01, 0.9550E+01, 0.9772E+01, 0.1000E+02
     & /)

      !=================================================================
      ! EX10 begins here!
      !=================================================================
      
      ! Limit X to [-K, K] range:  MIN: -9.999, MAX: 9.999
      Y    = MAX(-K, MIN(X,K))   

      ! Get integer and decimal part
      K1   = INT(Y)
      K2   = INT(100*(Y-K1))

      ! Calculate EXP function
      EX10 = AINT10(K1+10) * ADEC10(K2+100)

      ! Return to calling program
      END FUNCTION EX10

!------------------------------------------------------------------------------

      SUBROUTINE PUSHERR( IERR, ERRINF )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE PUSHERR
! *** THIS SUBROUTINE SAVES AN ERROR MESSAGE IN THE ERROR STACK
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      INTEGER,          INTENT(IN) :: IERR
      CHARACTER(LEN=*), INTENT(IN) :: ERRINF

      ! Local variables
      

      !=================================================================
      ! PUSHERR begins here!
      !=================================================================

      ! SAVE ERROR CODE IF THERE IS ANY SPACE
      IF ( NOFER .LT. NERRMX ) THEN   
         NOFER         = NOFER + 1 
         ERRSTK(NOFER) = IERR
         ERRMSG(NOFER) = ERRINF   
         STKOFL        =.FALSE.

      ! STACK OVERFLOW
      ELSE
         STKOFL        =.TRUE.      
      ENDIF

      ! Return to calling program
      END SUBROUTINE PUSHERR

!------------------------------------------------------------------------------

      SUBROUTINE ISERRINF( ERRSTKI, ERRMSGI, NOFERI, STKOFLI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISERRINF
! *** THIS SUBROUTINE OBTAINS A COPY OF THE ERROR STACK (& MESSAGES) 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      CHARACTER(LEN=40), INTENT(OUT) :: ERRMSGI(NERRMX)
      INTEGER,           INTENT(OUT) :: ERRSTKI(NERRMX)
      INTEGER,           INTENT(OUT) :: NOFERI
      LOGICAL,           INTENT(OUT) :: STKOFLI

      ! Local variables
      INTEGER           :: I

      !=================================================================
      ! ISERRINF begins here!
      ! 
      ! OBTAIN WHOLE ERROR STACK
      !=================================================================
      DO I = 1, NOFER              ! Error messages & codes
         ERRSTKI(I) = ERRSTK(I)
         ERRMSGI(I) = ERRMSG(I)
      ENDDO

      STKOFLI = STKOFL
      NOFERI  = NOFER

      ! Return to calling program
      END SUBROUTINE ISERRINF
      
!------------------------------------------------------------------------------

      SUBROUTINE ERRSTAT( IO, IERR, ERRINF )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE ERRSTAT
! *** THIS SUBROUTINE REPORTS ERROR MESSAGES TO UNIT 'IO'
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      INTEGER           :: IO, IERR
      CHARACTER(LEN=*)  :: ERRINF

      ! Local variables
      INTEGER           :: IOK, IEND
      CHARACTER(LEN=4)  :: CER 
      CHARACTER(LEN=29) :: NCIS = 'NO CONVERGENCE IN SUBROUTINE '
      CHARACTER(LEN=27) :: NCIF = 'NO CONVERGENCE IN FUNCTION ' 
      CHARACTER(LEN=26) :: NSIS = 'NO SOLUTION IN SUBROUTINE '  
      CHARACTER(LEN=24) :: NSIF = 'NO SOLUTION IN FUNCTION '  

      !=================================================================
      ! ERRSTAT begins here!
      !
      ! WRITE ERROR IN CHARACTER
      !=================================================================
      WRITE (CER,'(I4)') IERR
      CALL RPLSTR (CER, ' ', '0',IOK)   ! REPLACE BLANKS WITH ZEROS
      CALL CHRBLN (ERRINF, IEND)        ! LAST POSITION OF ERRINF CHAR

      !=================================================================
      ! WRITE ERROR TYPE (FATAL, WARNING) 
      !=================================================================    
      IF ( IERR .EQ. 0 ) THEN
         WRITE(IO,1000) 'NO ERRORS DETECTED '
         GOTO 10

      ELSE IF ( IERR .LT. 0 ) THEN
         WRITE(IO,1000) 'ERROR STACK EXHAUSTED '
         GOTO 10

      ELSE IF ( IERR .GT. 1000 ) THEN
         WRITE(IO,1100) 'FATAL',CER

      ELSE
         WRITE(IO,1100) 'WARNING',CER

      ENDIF

      !=================================================================
      ! WRITE ERROR MESSAGE
      !================================================================= 

      ! FATAL MESSAGES
      IF ( IERR .EQ. 1001 ) THEN 
         CALL CHRBLN (SCASE, IEND)
         WRITE(IO,1000) 'CASE NOT SUPPORTED IN CALCMR ['//SCASE(1:IEND)
     &                  //']'

      ELSE IF ( IERR .EQ. 1002 ) THEN 
         CALL CHRBLN (SCASE, IEND)
         WRITE(IO,1000) 'CASE NOT SUPPORTED ['//SCASE(1:IEND)//']'

      ! WARNING MESSAGES
      ELSE IF ( IERR .EQ. 0001 ) THEN 
         WRITE(IO,1000) NSIS,ERRINF

      ELSE IF ( IERR .EQ. 0002 ) THEN 
         WRITE(IO,1000) NCIS,ERRINF

      ELSEIF ( IERR .EQ. 0003 ) THEN 
         WRITE(IO,1000) NSIF,ERRINF

      ELSE IF ( IERR .EQ. 0004 ) THEN 
         WRITE(IO,1000) NCIF,ERRINF

      ELSE IF ( IERR .EQ. 0019 ) THEN
         WRITE(IO,1000) 'HNO3(aq) AFFECTS H+, WHICH '//
     &                  'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE(IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSE IF ( IERR .EQ. 0020 ) THEN
         IF (W(4).GT.TINY .AND. W(5).GT.TINY) THEN
            WRITE(IO,1000) 'HSO4-SO4 EQUILIBRIUM MIGHT AFFECT HNO3,'
     &                    //'HCL DISSOLUTION'
         ELSE
            WRITE(IO,1000) 'HSO4-SO4 EQUILIBRIUM MIGHT AFFECT NH3 '
     &                    //'DISSOLUTION'
         ENDIF
         WRITE(IO,1000) 'DIRECT DECREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSE IF ( IERR .EQ. 0021 ) THEN
         WRITE(IO,1000) 'HNO3(aq),HCL(aq) AFFECT H+, WHICH '//
     &                  'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE(IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSE IF ( IERR .EQ. 0022 ) THEN
         WRITE(IO,1000) 'HCL(g) EQUILIBRIUM YIELDS NONPHYSICAL '//
     &                  'DISSOLUTION'
         WRITE(IO,1000) 'A TINY AMOUNT [',ERRINF(1:IEND),'] IS '//
     &                  'ASSUMED TO BE DISSOLVED'

      ELSE IF ( IERR .EQ. 0033 ) THEN
         WRITE(IO,1000) 'HCL(aq) AFFECTS H+, WHICH '//
     &                  'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE(IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'

      ELSEIF ( IERR .EQ. 0050 ) THEN
         WRITE(IO,1000) 'TOO MUCH SODIUM GIVEN AS INPUT.'
         WRITE(IO,1000) 'REDUCED TO COMPLETELY NEUTRALIZE SO4,Cl,NO3.'
         WRITE(IO,1000) 'EXCESS SODIUM IS IGNORED.'

      ELSE
         WRITE(IO,1000) 'NO DIAGNOSTIC MESSAGE AVAILABLE'
      ENDIF

      ! Return
10    RETURN

      ! FORMAT statements
1000  FORMAT (1X,A:A:A:A:A)
1100  FORMAT (1X,A,' ERROR [',A4,']:')

      ! Return to calling program
      END SUBROUTINE ERRSTAT

!------------------------------------------------------------------------------

      SUBROUTINE ISRP2F( WI, RHI, TEMPI )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE ISRP2F
! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF 
!     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM. 
!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!     THE AMBIENT RELATIVE HUMIDITY.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8  :: WI(NCOMP), RHI, TEMPI

      !=================================================================
      ! ISRP2F begins here!
      !=================================================================
      
      ! Allocate all module arrays and some common blocks
      CALL INIT2( WI, RHI, TEMPI )

      ! CALCULATE SULFATE RATIO
      SULRAT = W(3)/W(2)

      !=================================================================
      ! FIND CALCULATION REGIME FROM (SULRAT,RH)
      !
      ! SULFATE POOR 
      !=================================================================
      IF ( 2.0 .LE. SULRAT ) THEN                

         ! Only liquid (metastable)
         IF ( METSTBL .EQ. 1 ) THEN
            SCASE = 'D3'
            CALL CALCD3         
         ELSE

            ! NH42SO4,NH4NO3       ; case D1
            IF ( RH .LT. DRNH4NO3 ) THEN    
               SCASE = 'D1'
               CALL CALCD1      
     
            ! NH42S4               ; case D2
            ELSE IF ( DRNH4NO3 .LE. RH .AND. RH .LT. DRNH42S4 ) THEN         
               SCASE = 'D2'
               CALL CALCD2      
     
            ! Only liquid          ; case D3
            ELSE IF ( DRNH42S4 .LE. RH ) THEN
               SCASE = 'D3'
               CALL CALCD3      
            ENDIF
         ENDIF

      !=================================================================
      ! SULFATE RICH (NO ACID)
      !
      ! FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
      ! THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
      ! SUBROUTINES CALCB? ARE CALLED, AND THEN THE NITRIC ACID IS 
      ! DISSOLVED FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
      !=================================================================
      ELSE IF ( 1.0 .LE. SULRAT .AND. SULRAT .LT. 2.0 ) THEN 

         IF ( METSTBL .EQ. 1 ) THEN
            SCASE = 'B4'
            CALL CALCB4         ! Only liquid (metastable)
            SCASE = 'E4'
         ELSE
            
            ! NH4HSO4,LC,NH42SO4   ; case E1
            IF ( RH .LT. DRNH4HS4 ) THEN         
               SCASE = 'B1'
               CALL CALCB1      
               SCASE = 'E1'

            ! LC,NH42S4            ; case E2
            ELSE IF ( DRNH4HS4 .LE. RH .AND. RH .LT. DRLC ) THEN         
               SCASE = 'B2'
               CALL CALCB2      
               SCASE = 'E2'
               
            ! NH42S4               ; case E3
            ELSE IF ( DRLC .LE. RH .AND. RH .LT. DRNH42S4 ) THEN         
               SCASE = 'B3'
               CALL CALCB3      
               SCASE = 'E3'
               
            ! Only liquid          ; case E4
            ELSE IF ( DRNH42S4 .LE. RH ) THEN         
               SCASE = 'B4'
               CALL CALCB4      
               SCASE = 'E4'
            ENDIF
         ENDIF
         
         ! HNO3(g) DISSOLUTION
         CALL CALCNA            

      !=================================================================
      ! SULFATE RICH (FREE ACID)
      !
      ! FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
      ! THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM
      ! SUBROUTINE CALCC? IS CALLED, AND THEN THE NITRIC ACID IS 
      ! DISSOLVED FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
      !=================================================================
      ELSE IF ( SULRAT .LT. 1.0 ) THEN             

         ! Only liquid (metastable)
         IF ( METSTBL .EQ. 1 ) THEN
            SCASE = 'C2'
            CALL CALCC2            
            SCASE = 'F2'

         ELSE

            ! NH4HSO4              ; case F1
            IF ( RH .LT. DRNH4HS4 ) THEN         
               SCASE = 'C1'
               CALL CALCC1         
               SCASE = 'F1'

            ! Only liquid          ; case F2
            ELSE IF ( DRNH4HS4 .LE. RH ) THEN         
               SCASE = 'C2'
               CALL CALCC2              
               SCASE = 'F2'
            ENDIF
         ENDIF

         ! HNO3(g) DISSOLUTION
         CALL CALCNA            
      ENDIF

      ! Return to calling program 
      END SUBROUTINE ISRP2F

!------------------------------------------------------------------------------

      SUBROUTINE CALCB4
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB4
! *** CASE B4 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
!
!     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
!     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
!     AND THAT CALCULATED FROM ELECTRONEUTRALITY.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      INTEGER :: I
      REAL*8  :: AK1, BET, GAM, BB, CC, DD

      !=================================================================
      ! CALCB4 begins here!
      !=================================================================

      ! Set flags to solve equations
      FRST       = .TRUE.
      CALAIN     = .TRUE.
      CALAOU     = .TRUE.

      !=================================================================
      ! CALCULATE WATER CONTENT
      !=================================================================

      ! Get dry salt content, and use for water.
      CALL CALCB1A         

      MOLALR(13) = CLC       
      MOLALR(9)  = CNH4HS4   
      MOLALR(4)  = CNH42S4   
      CLC        = ZERO
      CNH4HS4    = ZERO
      CNH42S4    = ZERO
      WATER      = MOLALR(13)/M0(13)+MOLALR(9)/M0(9)+MOLALR(4)/M0(4)

      ! NH4I
      MOLAL(3)   = W(3)   

      DO I = 1, NSWEEP
         AK1   = XK1*((GAMA(8)/GAMA(7))**2.)*(WATER/GAMA(7))
         BET   = W(2)
         GAM   = MOLAL(3)
         BB    = BET + AK1 - GAM
         CC    =-AK1*BET
         DD    = BB*BB - 4.D0*CC

         ! Speciation & water content
         MOLAL (5) = MIN(0.5*(-BB + SQRT(DD)), W(2))           ! SO4I
         MOLAL (6) = MAX(TINY,MIN(W(2)-MOLAL(5),W(2)))         ! HSO4I
         MOLAL (1) = MAX(TINY,MIN(AK1*MOLAL(6)/MOLAL(5),W(2))) ! HI
         CALL CALCMR                                           ! Water content

         ! Calculate activities or terminate internal loop 
         IF ( .NOT. CALAIN ) GOTO 30
         CALL CALCACT
      ENDDO

      ! Return to calling program
30    CONTINUE
      END SUBROUTINE CALCB4

!------------------------------------------------------------------------------

      SUBROUTINE CALCB3
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB3
! *** CASE B3 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!     3. SOLIDS POSSIBLE: (NH4)2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8 :: X, Y, TLC, TNH42S4, TNH4HS4

      !=================================================================
      ! CALCB3 begins here!
      !=================================================================

      ! Calculate equivalent amount of HSO4 and SO4
      X = MAX(2*W(2)-W(3), ZERO)   ! Equivalent NH4HSO4
      Y = MAX(W(3)  -W(2), ZERO)   ! Equivalent NH42SO4

      ! Calculate species according to relative abundance of HSO4
      IF ( X .LT. Y ) THEN               ! LC is the MIN (x,y)
         SCASE   = 'B3 ; SUBCASE 1'
         TLC     = X
         TNH42S4 = Y-X
         CALL CALCB3A (TLC,TNH42S4)      ! LC + (NH4)2SO4 
      ELSE
         SCASE   = 'B3 ; SUBCASE 2'
         TLC     = Y
         TNH4HS4 = X-Y
         CALL CALCB3B (TLC,TNH4HS4)      ! LC + NH4HSO4
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB3

!------------------------------------------------------------------------------

      SUBROUTINE CALCB3A( TLC, TNH42S4 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB3A
! *** CASE B3 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!     3. SOLIDS POSSIBLE: (NH4)2SO4
!
!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
!     AMOUNT OF SOLID (NH4)2SO4 DISSOLVED IN THE LIQUID PHASE.
!     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB3A CALCULATES THE
!     AMOUNT OF H+ PRODUCED (BASED ON THE SO4 RELEASED INTO THE
!     SOLUTION). THE SOLUBILITY PRODUCT OF (NH4)2SO4 IS USED AS THE 
!     OBJECTIVE FUNCTION.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8  :: TLC, TNH42S4

      ! Local variables
      INTEGER :: I
      REAL*8  :: FUNCB3A
      REAL*8  :: ZLO, ZHI, Z1, Z2, Z3, Y1, Y2, Y3, YLO, YHI, DZ, ZK

      !=================================================================
      ! CALCB3A begins here!
      !=================================================================
      CALAOU = .TRUE.         ! Outer loop activity calculation flag
      ZLO    = ZERO           ! MIN DISSOLVED (NH4)2SO4
      ZHI    = TNH42S4        ! MAX DISSOLVED (NH4)2SO4

      !=================================================================
      ! INITIAL VALUES FOR BISECTION (DISSOLVED (NH4)2SO4 
      !=================================================================
      Z1 = ZLO
      Y1 = FUNCB3A (Z1, TLC, TNH42S4)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DZ = (ZHI-ZLO)/FLOAT(NDIV)
      DO I = 1, NDIV
         Z2 = Z1+DZ
         Y2 = FUNCB3A (Z2, TLC, TNH42S4)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         Z1 = Z2
         Y1 = Y2
      ENDDO

      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION FOUND 
      !=================================================================
      YHI= Y1                    ! Save Y-value at HI position
      IF (ABS(Y2) .LT. EPS) THEN ! X2 IS A SOLUTION 
         RETURN

      ! { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
      ELSE IF ( YLO .LT. ZERO .AND. YHI. LT. ZERO ) THEN
         Z1 = ZHI
         Z2 = ZHI
         GOTO 40

      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
      ELSE IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         Z1 = ZLO
         Z2 = ZLO
         GOTO 40

      ! Warning error: no solution
      ELSE
         CALL PUSHERR( 0001, 'CALCB3A' ) 
         RETURN

      ENDIF

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
 20   CONTINUE

      DO I = 1, MAXIT
         Z3 = 0.5*(Z1+Z2)
         Y3 = FUNCB3A (Z3, TLC, TNH42S4)

         ! (Y1*Y3 .LE. ZERO)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  
            Y2    = Y3
            Z2    = Z3
         ELSE
            Y1    = Y3
            Z1    = Z3
         ENDIF
         IF (ABS(Z2-Z1) .LE. EPS*Z1) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR( 0002, 'CALCB3A' )

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
40    CONTINUE
      ZK = 0.5*(Z1+Z2)
      Y3 = FUNCB3A (ZK, TLC, TNH42S4)

      ! Return to calling program
      END SUBROUTINE CALCB3A

!------------------------------------------------------------------------------

      FUNCTION FUNCB3A( ZK, Y, X )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCB3A
! *** CASE B3 ; SUBCASE 1
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B3
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA3.
!
!******************************************************************************
!
      ! Arguments
      REAL*8  :: ZK, Y, X

      ! Local variables
      INTEGER :: I
      REAL*8  :: KK, GRAT1, DD, FUNCB3A
      
      !=================================================================
      ! FUNCB3A begins here!
      !=================================================================
      FRST   = .TRUE.
      CALAIN = .TRUE.

      ! Solve equations ; with iterations for activity coef.
      DO I = 1, NSWEEP
         GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         DD    = SQRT( (ZK+GRAT1+Y)**2. + 4.0*Y*GRAT1)
         KK    = 0.5*(-(ZK+GRAT1+Y) + DD )

         ! Speciation & water content
         MOLAL (1) = KK                ! HI
         MOLAL (5) = KK+ZK+Y           ! SO4I
         MOLAL (6) = MAX (Y-KK, TINY)  ! HSO4I
         MOLAL (3) = 3.0*Y+2*ZK        ! NH4I
         CNH42S4   = X-ZK              ! Solid (NH4)2SO4
         CALL CALCMR                   ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 30
         ENDIF
      ENDDO

      !=================================================================
      ! calculate objective function
      !=================================================================
30    CONTINUE
      !FUNCB3A= ( SO4I*NH4I**2.0 )/( XK7*(WATER/GAMA(4))**3.0 )      
      FUNCB3A= MOLAL(5)*MOLAL(3)**2.0
      FUNCB3A= FUNCB3A/(XK7*(WATER/GAMA(4))**3.0) - ONE
  
      ! Return to calling program
      END FUNCTION FUNCB3A

!------------------------------------------------------------------------------

      SUBROUTINE CALCB3B( Y, X )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB3B
! *** CASE B3 ; SUBCASE 2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. LIQUID PHASE ONLY IS POSSIBLE
!
!     SPECIATION CALCULATIONS IS BASED ON THE HSO4 <--> SO4 EQUILIBRIUM. 
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: Y, X

      ! Local variables
      INTEGER            :: I
      REAL*8             :: KK, GRAT1, DD

      !=================================================================
      ! CALCB3B begins here!
      !=================================================================

      ! Outer loop activity calculation flag
      CALAOU = .FALSE.        
      FRST   = .FALSE.
      CALAIN = .TRUE.

      ! Solve equations ; with iterations for activity coef.
      DO I = 1, NSWEEP
         GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         DD    = SQRT( (GRAT1+Y)**2. + 4.0*(X+Y)*GRAT1)
         KK    = 0.5*(-(GRAT1+Y) + DD )

         ! Speciation & water content
         MOLAL (1) = KK                   ! HI
         MOLAL (5) = Y+KK                 ! SO4I
         MOLAL (6) = MAX (X+Y-KK, TINY)   ! HSO4I
         MOLAL (3) = 3.0*Y+X              ! NH4I
         CALL CALCMR                      ! Water content

         ! Calculate activities or terminate internal loop
         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT     
      ENDDO

      ! Return to calling program
30    CONTINUE
      END SUBROUTINE CALCB3B

!------------------------------------------------------------------------------

      SUBROUTINE CALCB2
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2
! *** CASE B2 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : LC, (NH4)2SO4
!
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON THE SULFATE RATIO:
!     1. WHEN BOTH LC AND (NH4)2SO4 ARE POSSIBLE (SUBROUTINE CALCB2A)
!     2. WHEN ONLY LC IS POSSIBLE (SUBROUTINE CALCB2B).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      REAL*8 :: X, Y

      !=================================================================
      ! CALCB2 begins here!
      !=================================================================

      ! CALCULATE EQUIVALENT AMOUNT OF HSO4 AND SO4
      X = MAX(2*W(2)-W(3), TINY)   ! Equivalent NH4HSO4
      Y = MAX(W(3)  -W(2), TINY)   ! Equivalent NH42SO4

      ! CALCULATE SPECIES ACCORDING TO RELATIVE ABUNDANCE OF HSO4 *********
      IF (X.LE.Y) THEN             ! LC is the MIN (x,y)
         SCASE = 'B2 ; SUBCASE 1'
         CALL CALCB2A (X,Y-X)      ! LC + (NH4)2SO4 POSSIBLE
      ELSE
         SCASE = 'B2 ; SUBCASE 2'
         CALL CALCB2B (Y,X-Y)      ! LC ONLY POSSIBLE
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB2

!------------------------------------------------------------------------------

      SUBROUTINE CALCB2A( TLC, TNH42S4 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2
! *** CASE B2 ; SUBCASE A. 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. SOLID PHASE ONLY POSSIBLE
!     3. SOLIDS POSSIBLE: LC, (NH4)2SO4
!
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE 
!
!     FOR SOLID CALCULATIONS, A MATERIAL BALANCE BASED ON THE STOICHIMETRIC
!     PROPORTION OF AMMONIUM AND SULFATE IS DONE TO CALCULATE THE AMOUNT 
!     OF LC AND (NH4)2SO4 IN THE SOLID PHASE.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: TLC, TNH42S4

      !=================================================================
      ! CALCB2A begins here!
      !=================================================================

      ! Regime depends upon the ambient relative humidity
      IF ( RH .LT. DRMLCAS ) THEN    

         ! Solids possible only
         SCASE   = 'B2 ; SUBCASE A1'    
         CLC     = TLC
         CNH42S4 = TNH42S4
         SCASE   = 'B2 ; SUBCASE A1'
      ELSE

         ! Liquid & solid phase possible
         SCASE = 'B2 ; SUBCASE A2'
         CALL CALCB2A2( TLC, TNH42S4 )   
         SCASE = 'B2 ; SUBCASE A2'
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB2A

!------------------------------------------------------------------------------

      SUBROUTINE CALCB2A2 (TLC, TNH42S4)
! 
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2A2
! *** CASE B2 ; SUBCASE A2. 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. SOLID PHASE ONLY POSSIBLE
!     3. SOLIDS POSSIBLE: LC, (NH4)2SO4
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE CALCB2A1) AND THE
!     THE SOLID WITH LIQUID PHASE (SUBROUTINE CALCB3).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: TLC, TNH42S4

      ! Local variables
      REAL*8  :: WF, ONEMWF, CLCO, CNH42SO

      !=================================================================
      ! CALCB2A2 begins here!
      !=================================================================

      ! Find weight factor
      IF ( WFTYP .EQ. 0 ) THEN
         WF = ZERO
      ELSE IF ( WFTYP .EQ. 1 ) THEN
         WF = 0.5D0
      ELSE
         WF = ( DRLC - RH ) / ( DRLC - DRMLCAS )
      ENDIF

      ONEMWF  = ONE - WF

      !=================================================================
      ! FIND FIRST SECTION ; DRY ONE
      !=================================================================      
      CLCO     = TLC                     ! FIRST (DRY) SOLUTION   
      CNH42SO  = TNH42S4

      !=================================================================
      ! FIND SECOND SECTION ; DRY & LIQUID 
      !=================================================================
      CLC     = ZERO
      CNH42S4 = ZERO
      CALL CALCB3                        ! SECOND (LIQUID) SOLUTION

      !=================================================================
      ! FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS
      !=================================================================
      MOLAL(1)= ONEMWF*MOLAL(1)                                   ! H+
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + 3.D0*(CLCO-CLC)) ! NH4+
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               ! SO4--
      MOLAL(6)= ONEMWF*(CLCO-CLC)                                 ! HSO4-

      WATER   = ONEMWF*WATER

      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4

      ! Return to calling program
      END SUBROUTINE CALCB2A2

!------------------------------------------------------------------------------

      SUBROUTINE CALCB2B( TLC, TNH4HS4 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB2
! *** CASE B2 ; SUBCASE B 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!     3. SOLIDS POSSIBLE: LC
!
!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
!     AMOUNT OF SOLID LC DISSOLVED IN THE LIQUID PHASE.
!     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB2A CALCULATES THE
!     AMOUNT OF H+ PRODUCED (BASED ON THE HSO4, SO4 RELEASED INTO THE
!     SOLUTION). THE SOLUBILITY PRODUCT OF LC IS USED AS THE OBJECTIVE 
!     FUNCTION.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!

      ! Arguments
      REAL*8  :: TLC, TNH4HS4

      ! 
      INTEGER :: I
      REAL*8  :: FUNCB2B
      REAL*8  :: ZLO, ZHI, X1, X2, X3, Y1, Y2, Y3, YLO, YHI, DX

      !=================================================================
      ! CALCB2B begins here!
      !=================================================================
      CALAOU = .TRUE.       ! Outer loop activity calculation flag
      ZLO    = ZERO
      ZHI    = TLC          ! High limit: all of it in liquid phase

      !=================================================================
      ! INITIAL VALUES FOR BISECTION 
      !=================================================================      
      X1 = ZHI
      Y1 = FUNCB2B (X1,TNH4HS4,TLC)
      IF (ABS(Y1).LE.EPS) RETURN
      YHI= Y1                        ! Save Y-value at Hi position

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO 
      !=================================================================
      DX = (ZHI-ZLO)/NDIV
      DO I = 1, NDIV
         X2 = X1-DX
         Y2 = FUNCB2B (X2,TNH4HS4,TLC)

         ! (Y1*Y2.LT.ZERO)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  

         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================      
      ! NO SUBDIVISION WITH SOLUTION FOUND 
      !=================================================================
      YLO= Y1                        ! Save Y-value at LO position
      IF ( ABS(Y2) .LT. EPS ) THEN   ! X2 IS A SOLUTION 
         RETURN

      !=================================================================
      ! { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
      !=================================================================
      ELSE IF ( YLO .LT. ZERO .AND. YHI .LT. ZERO ) THEN
         X1 = ZHI
         X2 = ZHI
         GOTO 40

      !=================================================================
      ! { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
      !=================================================================
      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         X1 = ZLO
         X2 = ZLO
         GOTO 40

      ELSE
         ! Warning error: no solution
         CALL PUSHERR( 0001, 'CALCB2B' )    
         RETURN

      ENDIF

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
20    CONTINUE

      DO I = 1, MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCB2B (X3,TNH4HS4,TLC)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR( 0002, 'CALCB2B' )    

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
40    CONTINUE
      X3 = 0.5*(X1+X2)
      Y3 = FUNCB2B (X3,TNH4HS4,TLC)

      ! Return to calling program
      END SUBROUTINE CALCB2B

!------------------------------------------------------------------------------

      FUNCTION FUNCB2B( X, TNH4HS4, TLC )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCB2B
! *** CASE B2 ; 
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B2 ; SUBCASE 2
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCB2B.
!
!******************************************************************************
!
      ! Arguments
      REAL*8  :: X, TNH4HS4, TLC

      ! Local variables
      INTEGER :: I
      REAL*8  :: GRAT2, PARM, DELTA, OMEGA, FUNCB2B

      !=================================================================
      ! FUNCB2B begins here!
      !=================================================================
      FRST   = .TRUE.
      CALAIN = .TRUE.

      ! Solve equation
      DO I = 1, NSWEEP
         GRAT2 = XK1*WATER*(GAMA(8)/GAMA(7))**2./GAMA(7)
         PARM  = X+GRAT2
         DELTA = PARM*PARM + 4.0*(X+TNH4HS4)*GRAT2 ! Diakrinousa
         OMEGA = 0.5*(-PARM + SQRT(DELTA))         ! Thetiki riza (ie:H+>0)

         ! Speciation & water content
         MOLAL (1) = OMEGA                         ! HI
         MOLAL (3) = 3.0*X+TNH4HS4                 ! NH4I
         MOLAL (5) = X+OMEGA                       ! SO4I
         MOLAL (6) = MAX (X+TNH4HS4-OMEGA, TINY)   ! HSO4I
         CLC       = MAX(TLC-X,ZERO)               ! Solid LC
         CALL CALCMR                               ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 30
         ENDIF
      ENDDO
      
      !=================================================================
      ! Calculate objective function
      !=================================================================
30    CONTINUE
      !FUNCB2B= ( NH4I**3.*SO4I*HSO4I )/( XK13*(WATER/GAMA(13))**5. )
      FUNCB2B= (MOLAL(3)**3.)*MOLAL(5)*MOLAL(6)
      FUNCB2B= FUNCB2B/(XK13*(WATER/GAMA(13))**5.) - ONE

      ! Return to calling program
      END FUNCTION FUNCB2B

!------------------------------------------------------------------------------

      SUBROUTINE CALCB1
! 
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB1
! *** CASE B1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : LC, (NH4)2SO4, NH4HSO4
!
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCB1A)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      !=================================================================
      ! CALCB1 begins here!
      !=================================================================

      ! Regime depends upon the ambient relative humidity
      IF ( RH .LT. DRMLCAB ) THEN  

         ! Solid phase only possible
         SCASE = 'B1 ; SUBCASE 1'  
         CALL CALCB1A              
         SCASE = 'B1 ; SUBCASE 1'

      ELSE

         ! Liquid & solid phase possible
         SCASE = 'B1 ; SUBCASE 2'
         CALL CALCB1B              
         SCASE = 'B1 ; SUBCASE 2'

      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB1

!------------------------------------------------------------------------------

      SUBROUTINE CALCB1A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB1A
! *** CASE B1 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH
!     2. THERE IS NO LIQUID PHASE
!     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
!                         BUT NOT BOTH)
!
!     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE AMOUNT OF LC
!     IS CALCULATED FROM THE (NH4)2SO4 AND NH4HSO4 WHICH IS LEAST
!     ABUNDANT (STOICHIMETRICALLY). THE REMAINING EXCESS OF SALT 
!     IS MIXED WITH THE LC.  
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      REAL*8 :: X, Y, CLC

      !=================================================================
      ! CALCB1 begins here! 
      !=================================================================
      X = 2*W(2)-W(3)       ! Equivalent NH4HSO4
      Y = W(3)-W(2)         ! Equivalent (NH4)2SO4

      ! Calculate composition
      IF (X.LE.Y) THEN      ! LC is the MIN (x,y)
         CLC     = X        ! NH4HSO4 >= (NH4)2S04
         CNH4HS4 = ZERO
         CNH42S4 = Y-X
      ELSE
         CLC     = Y        ! NH4HSO4 <  (NH4)2S04
         CNH4HS4 = X-Y
         CNH42S4 = ZERO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCB1A

!------------------------------------------------------------------------------

      SUBROUTINE CALCB1B
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCB1B
! *** CASE B1 ; SUBCASE 2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
!                         BUT NOT BOTH)
!
!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE CALCB1A) AND THE
!     THE SOLID WITH LIQUID PHASE (SUBROUTINE CALCB2).
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************

      ! Local variables
      REAL*8 :: WF, ONEMWF, CLCO, CNH42SO, CNH4HSO

      !=================================================================
      ! CALCB1B begins here!
      !=================================================================

      ! Find weight factor
      IF (WFTYP.EQ.0) THEN
         WF = ZERO
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (DRNH4HS4-RH)/(DRNH4HS4-DRMLCAB)
      ENDIF
      ONEMWF  = ONE - WF

      !=================================================================
      ! FIND FIRST SECTION ; DRY ONE
      !=================================================================
      CALL CALCB1A
      CLCO     = CLC               ! FIRST (DRY) SOLUTION
      CNH42SO  = CNH42S4
      CNH4HSO  = CNH4HS4

      !=================================================================
      ! FIND SECOND SECTION ; DRY & LIQUID
      !=================================================================
      CLC     = ZERO
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CALL CALCB2                  ! SECOND (LIQUID) SOLUTION

      !=================================================================
      ! FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
      !=================================================================
      MOLAL(1)= ONEMWF*MOLAL(1)                                   ! H+
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + (CNH4HSO-CNH4HS4)  
     &                + 3.D0*(CLCO-CLC))                          ! NH4+
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               ! SO4--
      MOLAL(6)= ONEMWF*(CNH4HSO-CNH4HS4 + CLCO-CLC)               ! HSO4-

      WATER   = ONEMWF*WATER

      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4

      ! Return to calling program
      END SUBROUTINE CALCB1B

!------------------------------------------------------------------------------

      SUBROUTINE CALCC2
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCC2
! *** CASE C2 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS ONLY A LIQUID PHASE
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      INTEGER :: I
      REAL*8  :: LAMDA, KAPA, PSI, PARM, BB, CC

      !=================================================================
      ! CALCC2 begins here!
      !=================================================================
      CALAOU =.TRUE.         ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
      LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
      PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION

      ! Solve equations
      DO I = 1, NSWEEP
         PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         BB    = PSI+PARM
         CC    =-PARM*(LAMDA+PSI)
         KAPA  = 0.5*(-BB+SQRT(BB*BB-4.0*CC))

         ! Speciation & water content
         MOLAL(1) = PSI+KAPA                               ! HI
         MOLAL(3) = LAMDA                                  ! NH4I
         MOLAL(5) = KAPA                                   ! SO4I
         MOLAL(6) = MAX(LAMDA+PSI-KAPA, TINY)              ! HSO4I
         CH2SO4   = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3), ZERO)  ! Free H2SO4
         CALL CALCMR                                       ! Water content

         ! Calculate activities or terminate internal loop
         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT     
      ENDDO

      ! Return to calling program
30    CONTINUE
      END SUBROUTINE CALCC2

!------------------------------------------------------------------------------

      SUBROUTINE CALCC1
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCC1
! *** CASE C1 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE: NH4HSO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************!
      ! Local variables
      REAL*8  :: KLO, KHI, X1, X2, X3, Y1, Y2, Y3
      REAL*8  :: FUNCC1, YLO, YHI, DX
      INTEGER :: I

      !=================================================================
      ! CALCC1 begins here!
      !=================================================================
      CALAOU = .TRUE.    ! Outer loop activity calculation flag
      KLO    = TINY    
      KHI    = W(3)

      !=================================================================
      ! INITIAL VALUES FOR BISECTION
      !=================================================================
      X1 = KLO
      Y1 = FUNCC1(X1)
      IF ( ABS(Y1) .LE. EPS ) GOTO 50
      YLO = Y1

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX = ( KHI - KLO )/ FLOAT( NDIV )

      DO I = 1, NDIV
         X2 = X1+DX
         Y2 = FUNCC1 (X2)

         ! (Y1*Y2 .LT. ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y2) .LT. ZERO ) GOTO 20 

         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================
      ! NO SUBDIVISION WITH SOLUTION FOUND 
      !=================================================================

      ! Save Y-value at HI position
      YHI= Y2                 

      ! X2 IS A SOLUTION 
      IF ( ABS(Y2) .LT. EPS ) THEN   
         GOTO 50

      !=================================================================
      ! { YLO, YHI } < 0.0  SOLUTION IS ALWAYS UNDERSATURATED WITH NH4HS04
      !=================================================================
      ELSE IF ( YLO .LT. ZERO .AND. YHI .LT. ZERO ) THEN
         GOTO 50

      !=================================================================
      ! { YLO, YHI } > 0.0 SOLUTION IS ALWAYS SUPERSATURATED WITH NH4HS04
      !=================================================================
      ELSE IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         X1 = KLO
         X2 = KLO
         GOTO 40
      ELSE
         CALL PUSHERR (0001, 'CALCC1')    ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF

      !=================================================================
      ! CPERFORM BISECTION OF DISSOLVED NH4HSO4
      !=================================================================
20    CONTINUE

      DO I = 1, MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCC1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR( 0002, 'CALCC1' )    

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
 40   CONTINUE
      X3 = 0.5*(X1+X2)
      Y3 = FUNCC1 (X3)

      ! Return to calling program
 50   CONTINUE
      END SUBROUTINE CALCC1

!------------------------------------------------------------------------------

      FUNCTION FUNCC1( KAPA )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCC1
! *** CASE C1 ; 
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE C1
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCC1.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CANREGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: KAPA

      ! Local variables
      INTEGER            :: I
      REAL*8             :: LAMDA, PSI, PAR1, PAR2, BB, CC, FUNCC1

      !=================================================================
      ! FUNCC1 begins here!
      !=================================================================
      FRST   = .TRUE.
      CALAIN = .TRUE.

      !==================================================================
      ! Solve equations
      !==================================================================
      PSI = W(2)-W(3)
      DO I = 1, NSWEEP
         PAR1  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
         PAR2  = XK12*(WATER/GAMA(9))**2.0
         BB    = PSI + PAR1
         CC    =-PAR1*(PSI+KAPA)
         LAMDA = 0.5*(-BB+SQRT(BB*BB-4*CC))

         ! Save concentrations in molal array
         MOLAL(1) = PSI+LAMDA                    ! HI
         MOLAL(3) = KAPA                         ! NH4I
         MOLAL(5) = LAMDA                        ! SO4I
         MOLAL(6) = MAX (ZERO, PSI+KAPA-LAMDA)   ! HSO4I
         CNH4HS4  = MAX(W(3)-MOLAL(3), ZERO)     ! Solid NH4HSO4
         CH2SO4   = MAX(PSI, ZERO)               ! Free H2SO4
         CALL CALCMR                             ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 30
         ENDIF
      ENDDO

      !==================================================================
      ! Calculate zero function
      !==================================================================
 30   CONTINUE
      !###FUNCC1= (NH4I*HSO4I/PAR2) - ONE
      FUNCC1= (MOLAL(3)*MOLAL(6)/PAR2) - ONE

      ! Return to calling program
      END FUNCTION FUNCC1

!------------------------------------------------------------------------------

      SUBROUTINE CALCD3
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCD3
! *** CASE D3
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. THERE IS OLNY A LIQUID PHASE
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      INTEGER :: I
      REAL*8  :: PSI4LO, PSI4HI, X1, X2, X3, Y1, Y2, Y3, YLO, YHI, DX
      REAL*8  :: P4, YY, DELTA
      REAL*8  :: FUNCD3
 
      !=================================================================
      ! CALCD3 begins here!
      !=================================================================

      ! Find dry composition
      CALL CALCD1A

      ! Setup parameters
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3

      ! Assign initial PSI's
      PSI1 = CNH4NO3               
      PSI2 = CHI2
      PSI3 = ZERO   
      PSI4 = ZERO  

      MOLAL(5) = ZERO
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1

      ! Initial water
      CALL CALCMR                  

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit

      !=================================================================
      ! INITIAL VALUES FOR BISECTION
      !=================================================================
60    CONTINUE
      X1 = PSI4LO
      Y1 = FUNCD3( X1 )
      IF ( ABS(Y1) .LE. EPS ) RETURN
      YLO= Y1                 ! Save Y-value at HI position

      !=================================================================     
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX = ( PSI4HI - PSI4LO ) / FLOAT( NDIV )

      DO I = 1, NDIV
         X2 = X1+DX
         Y2 = FUNCD3 (X2)

         ! (Y1*Y2.LT.ZERO)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  

         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================     
      ! NO SUBDIVISION WITH SOLUTION FOUND 
      !=================================================================
      YHI= Y1                        ! Save Y-value at Hi position
      IF ( ABS(Y2) .LT. EPS ) THEN   ! X2 IS A SOLUTION 
         RETURN

      !=================================================================
      ! { YLO, YHI } < 0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
      !
      ! Physically I dont know when this might happen, but I have put 
      ! this branch in for completeness. I assume there is no solution; 
      ! all NO3 goes to the gas phase. (rjp)
      !=================================================================
      ELSE IF ( YLO .LT. ZERO .AND. YHI .LT. ZERO ) THEN
         P4 = TINY ! PSI4LO ! CHI4
         YY = FUNCD3(P4)
         GOTO 50

      !=================================================================
      ! { YLO, YHI } > 0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
      !
      ! This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate 
      ! evaporates and goes to the gas phase ; so I redefine the LO and 
      ! HI limits of PSI4 and proceed again with root tracking. (rjp)
      !=================================================================
      ELSE IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN

         PSI4HI = PSI4LO

         ! No solution; some NH3 evaporates
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) 

         IF ( PSI4LO .LT. -( PSI1 + PSI2 ) ) THEN

            ! Warning error: no solution
            CALL PUSHERR( 0001, 'CALCD3' )  
            RETURN
         ELSE

            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  ! Initial water
            GOTO 60                      ! Redo root tracking
         ENDIF
      ENDIF

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
20    CONTINUE

      DO I = 1, MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCD3 (X3)

         ! (Y1*Y3 .LE. ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y3) .LE. ZERO ) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ( ABS( X2 - X1 ) .LE. EPS * ABS( X1 ) ) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR (0002, 'CALCD3')    

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================    
40    CONTINUE
      X3 = 0.5*(X1+X2)
      Y3 = FUNCD3( X3 )

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN
      !================================================================= 
50    CONTINUE

      IF ( MOLAL(1) .GT. TINY ) THEN
         CALL CALCHS4( MOLAL(1), MOLAL(5), ZERO, DELTA )
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCD3

!------------------------------------------------------------------------------

      FUNCTION FUNCD3( P4 )
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCD3
! *** CASE D3 
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ; 
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD3.
!
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: P4

      ! Local variables
      INTEGER :: I
      REAL*8  :: BB, DENM, ABB, AHI, FUNCD3

      !=================================================================
      ! FUNCD3 begins here!
      !=================================================================
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4

      !=================================================================
      ! SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF
      !=================================================================
      DO I = 1, NSWEEP
         A2   = XK7*(WATER/GAMA(4))**3.0
         A3   = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A7   = XKW *RH*WATER*WATER

         PSI3 = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3 = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4) 
         PSI3 = MIN(MAX(PSI3, ZERO), CHI3)

         BB   = PSI4 - PSI3
!###old         AHI  = 0.5*(-BB + SQRT(BB*BB + 4.d0*A7)) ! This is correct also
!###        AHI  =2.0*A7/(BB+SQRT(BB*BB + 4.d0*A7)) ! Avoid overflow when HI->0

         DENM = BB+SQRT(BB*BB + 4.d0*A7)
         IF (DENM.LE.TINY) THEN       ! Avoid overflow when HI->0
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.0*A7/ABB ! Taylor expansion of SQRT
         ENDIF
         AHI = 2.0*A7/DENM

         ! SPECIATION & WATER CONTENT
         MOLAL (1) = AHI                             ! HI
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2         ! NH4I
         MOLAL (5) = PSI2                            ! SO4I
         MOLAL (6) = ZERO                            ! HSO4I
         MOLAL (7) = PSI3 + PSI1                     ! NO3I
         CNH42S4   = CHI2 - PSI2                     ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                            ! Solid NH4NO3
         GHNO3     = CHI3 - PSI3                     ! Gas HNO3
         GNH3      = CHI4 - PSI4                     ! Gas NH3
         CALL CALCMR                                 ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND.CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
      ENDDO

      !=================================================================
      ! Calculate objective function
      !=================================================================
20    CONTINUE
      !### FUNCD3= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE 
      FUNCD3= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 
      RETURN

      ! Return to calling program
      END FUNCTION FUNCD3

!------------------------------------------------------------------------------

      SUBROUTINE CALCD2
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCD2
! *** CASE D2
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!     3. SOLIDS POSSIBLE : (NH4)2SO4
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      INTEGER :: I
      REAL*8  :: PSI4LO, PSI4HI, X1, X2, X3, Y1, Y2, Y3, YLO, YHI, DX
      REAL*8  :: FUNCD2, P4, YY, DELTA

      !=================================================================
      ! CALCD2 begins here!
      !=================================================================

      ! Find dry composition
      CALL CALCD1A

      ! Setup parameters
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3

      PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
      PSI2 = CNH42S4
      PSI3 = ZERO   
      PSI4 = ZERO  

      MOLAL(5) = ZERO
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1
      CALL CALCMR                  ! Initial water

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit

      !=================================================================
      ! Initial values for bisection
      !=================================================================
 60   CONTINUE
      X1 = PSI4LO
      Y1 = FUNCD2 (X1)
      IF ( ABS( Y1 ) .LE. EPS ) RETURN

      ! Save Y-value at HI position
      YLO = Y1                  

      !=================================================================
      ! ROOT TRACKING ; FOR THE RANGE OF HI AND LO
      !=================================================================
      DX   = (PSI4HI-PSI4LO)/FLOAT(NDIV)

      DO I = 1, NDIV
         X2 = X1+DX
         Y2 = FUNCD2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) THEN

            ! This is done, in case if Y(PSI4LO)>0, 
            ! but Y(PSI4LO+DX) < 0 (i.e.undersat)
            IF (Y1 .LE. Y2) GOTO 20 ! (Y1*Y2.LT.ZERO)
         ENDIF
         X1 = X2
         Y1 = Y2
      ENDDO

      !=================================================================
      ! No subdivision with solution found
      !=================================================================
      YHI= Y1                          ! Save Y-value at Hi position
      IF ( ABS( Y2 ) .LT. EPS ) THEN   ! X2 IS A SOLUTION 
         RETURN

      !=================================================================
      ! { YLO, YHI } < 0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
      !
      ! Physically I dont know when this might happen, but I have put 
      ! this branch in for completeness. I assume there is no solution; 
      ! all NO3 goes to the gas phase. (rjp)
      !=================================================================
      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = TINY ! PSI4LO ! CHI4
         YY = FUNCD2(P4)
         GOTO 50

      !=================================================================
      ! { YLO, YHI } > 0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
      ! 
      ! This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate 
      ! evaporatesand goes to the gas phase ; so I redefine the LO and 
      ! HI limits of PSI4 and proceed again with root tracking. (rjp)
      !=================================================================
      ELSE IF ( YLO .GT. ZERO .AND. YHI .GT. ZERO ) THEN
         PSI4HI = PSI4LO

         ! No solution; some NH3 evaporates
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) 

         IF ( PSI4LO .LT. -( PSI1 + PSI2 ) ) THEN

            ! Warning error: no solution
            CALL PUSHERR( 0001, 'CALCD2' )  
            RETURN
         ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  ! Initial water
            GOTO 60                      ! Redo root tracking
         ENDIF
      ENDIF

      !=================================================================
      ! PERFORM BISECTION
      !=================================================================
 20   CONTINUE

      DO I = 1, MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNCD2 (X3)

         ! (Y1*Y3 .LE. ZERO)
         IF ( SIGN(1.d0,Y1) * SIGN(1.d0,Y3) .LE. ZERO ) THEN  
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
      ENDDO

      ! Warning error: no convergence
      CALL PUSHERR (0002, 'CALCD2')    

      !=================================================================
      ! CONVERGED ; RETURN
      !=================================================================
 40   CONTINUE

      ! 0.5*(X1+X2)  ! Get "low" side, it's acidic soln.
      X3 = MIN(X1,X2)           
      Y3 = FUNCD2 (X3)

      !=================================================================
      ! CALCULATE HSO4 SPECIATION AND RETURN
      !=================================================================
 50   CONTINUE
      
      IF ( MOLAL(1) .GT. TINY ) THEN
         CALL CALCHS4( MOLAL(1), MOLAL(5), ZERO, DELTA )
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCD2

!------------------------------------------------------------------------------
      
      FUNCTION FUNCD2( P4 )
! 
!******************************************************************************
!
! *** ISORROPIA CODE
! *** FUNCTION FUNCD2
! *** CASE D2 
!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D2 ; 
!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD2.
!
!******************************************************************************
!
      REAL*8, INTENT(IN) :: P4

      ! Local variables
      INTEGER            :: I, ISLV
      REAL*8             :: PSI14, BB, DENM, ABB, AHI, FUNCD2

      !=================================================================
      ! FUNCD2 begins here!
      !=================================================================

      ! Reset activity coefficients to 0.1
      CALL RSTGAM       
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4
      PSI2   = CHI2

      !=================================================================       
      ! Solve equations ; with iterations for activity coef.
      !=================================================================
      DO I = 1, NSWEEP
         A2  = XK7*(WATER/GAMA(4))**3.0
         A3  = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A7  = XKW *RH*WATER*WATER

         IF ( CHI2 .GT. TINY .AND. WATER .GT. TINY ) THEN
            PSI14 = PSI1+PSI4

            CALL POLY3 (PSI14,0.25*PSI14**2.,-A2/4.D0, PSI2, ISLV)  ! PSI2
            IF ( ISLV.EQ.0 ) THEN
                PSI2 = MIN (PSI2, CHI2)
            ELSE
                PSI2 = ZERO
            ENDIF
         ENDIF

         PSI3  = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3  = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4) 
         !### PSI3  = MIN(MAX(PSI3, ZERO), CHI3)

         ! (BB > 0, acidic solution, <0 alkaline)
         BB   = PSI4 - PSI3 

         !==============================================================
         ! Do not change computation scheme for H+, 
         ! all others did not work well
         !==============================================================
         DENM = BB+SQRT(BB*BB + 4.d0*A7)

         ! Avoid overflow when HI->0
         IF ( DENM .LE. TINY ) THEN       
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.d0*A7/ABB ! Taylor expansion of SQRT
         ENDIF
         AHI = 2.d0*A7/DENM

         ! Speciation & water content
         MOLAL (1) = AHI                              ! HI
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2          ! NH4
         MOLAL (5) = PSI2                             ! SO4
         MOLAL (6) = ZERO                             ! HSO4
         MOLAL (7) = PSI3 + PSI1                      ! NO3
         CNH42S4   = CHI2 - PSI2                      ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                             ! Solid NH4NO3
         GHNO3     = CHI3 - PSI3                      ! Gas HNO3
         GNH3      = CHI4 - PSI4                      ! Gas NH3
         CALL CALCMR                                  ! Water content

         ! Calculate activities or terminate internal loop
         IF ( FRST .AND. CALAOU .OR. .NOT. FRST .AND. CALAIN ) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
      ENDDO

      !=================================================================
      ! Calculate objective function
      !=================================================================
 20   CONTINUE

      !### FUNCD2= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE 
      FUNCD2= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 

      ! Return to calling program
      END FUNCTION FUNCD2

!------------------------------------------------------------------------------

      SUBROUTINE CALCD1
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCD1
! *** CASE D1 
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!
!     THERE ARE TWO REGIMES DEFINED BY RELATIVE HUMIDITY:
!     1. RH < MDRH ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCD1A)
!     2. RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      !=================================================================
      ! CALCD1 begins here!
      !=================================================================
      IF ( RH .LT. DRMASAN ) THEN    

         ! Solid phase only possible
         SCASE = 'D1 ; SUBCASE 1'   
         CALL CALCD1A            
         SCASE = 'D1 ; SUBCASE 1'

      ELSE

         ! Liquid & solid phase possible
         SCASE = 'D1 ; SUBCASE 2'   
         CALL CALCMDRH (RH, DRMASAN, DRNH4NO3, CALCD1A, CALCD2)
         SCASE = 'D1 ; SUBCASE 2'

      ENDIF

      ! Return to calling program
      END SUBROUTINE CALCD1

!------------------------------------------------------------------------------

      SUBROUTINE CALCD1A
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE CALCD1A
! *** CASE D1 ; SUBCASE 1
!
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!
!     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
!     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
!     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
!     THE SOLID PHASE.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!******************************************************************************
!
      ! Local variables
      REAL*8  :: PARM, X, PS, OM, OMPS, DIAK, ZE

      !=================================================================
      ! CALCD1A begins here!
      !=================================================================
      PARM    = XK10/(R*TEMP)/(R*TEMP)

      ! Calculate NH4NO3 that volatizes
      CNH42S4 = W(2)                                    
      X       = MAX(ZERO, MIN(W(3)-2.0*CNH42S4, W(4)))  ! MAX NH4NO3
      PS      = MAX(W(3) - X - 2.0*CNH42S4, ZERO)
      OM      = MAX(W(4) - X, ZERO)

      OMPS    = OM+PS
      DIAK    = SQRT(OMPS*OMPS + 4.0*PARM)              ! DIAKRINOUSA
      ZE      = MIN(X, 0.5*(-OMPS + DIAK))              ! THETIKI RIZA

      ! Speciation
      CNH4NO3 = X  - ZE    ! Solid NH4NO3
      GNH3    = PS + ZE    ! Gas NH3
      GHNO3   = OM + ZE    ! Gas HNO3

      ! Return to calling program
      END SUBROUTINE CALCD1A

!------------------------------------------------------------------------------

      SUBROUTINE INIT_ISOROPIA
!
!******************************************************************************
!
! *** ISORROPIA CODE
! *** SUBROUTINE SETPARM
! *** THIS SUBROUTINE REDEFINES THE SOLUTION PARAMETERS OF ISORROPIA
!
! ======================== ARGUMENTS / USAGE ===========================
!
! *** NOTE: IF NEGATIVE VALUES ARE GIVEN FOR A PARAMETER, IT IS
!     IGNORED AND THE CURRENT VALUE IS USED INSTEAD.
! 
!  INPUT:
!  1. [WFTYPI] 
!     INTEGER variable.
!     Defines the type of weighting algorithm for the solution in Mutual 
!     Deliquescence Regions (MDR's):
!     0 - MDR's are assumed dry. This is equivalent to the approach 
!         used by SEQUILIB.
!     1 - The solution is assumed "half" dry and "half" wet throughout
!         the MDR.
!     2 - The solution is a relative-humidity weighted mean of the
!         dry and wet solutions (as defined in Nenes et al., 1998)
!
!  2. [IACALCI] 
!     INTEGER variable.
!     Method of activity coefficient calculation:
!     0 - Calculate coefficients during runtime
!     1 - Use precalculated tables
! 
!  3. [EPSI] 
!     DOUBLE PRECITION variable.
!     Defines the convergence criterion for all iterative processes
!     in ISORROPIA, except those for activity coefficient calculations
!     (EPSACTI controls that).
!
!  4. [MAXITI]
!     INTEGER variable.
!     Defines the maximum number of iterations for all iterative 
!     processes in ISORROPIA, except for activity coefficient calculations 
!     (NSWEEPI controls that).
!
!  5. [NSWEEPI]
!     INTEGER variable.
!     Defines the maximum number of iterations for activity coefficient 
!     calculations.
! 
!  6. [EPSACTI] 
!     DOUBLE PRECISION variable.
!     Defines the convergence criterion for activity coefficient 
!     calculations.
! 
!  7. [NDIV] 
!     INTEGER variable.
!     Defines the number of subdivisions needed for the initial root
!     tracking for the bisection method. Usually this parameter should 
!     not be altered, but is included for completeness.
!
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
!
!
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

      ! Local variables
      INTEGER :: AS 

      !=================================================================
      ! INIT_ISOROPIA begins here!
      !=================================================================
      ALLOCATE( ASRAT( NASRD ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'ASRAT' )
      ASRAT = 0d0

      ALLOCATE( ASSO4( NSO4S ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'ASSO4' )
      ASSO4 = 0d0

      ALLOCATE( AWAB( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWAB' ) 
      AWAB = 0d0

      ALLOCATE( AWAC( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWAC' ) 
      AWAC = 0d0

      ALLOCATE( AWAN( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWAN' ) 
      AWAN = 0d0

      ALLOCATE( AWAS( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWAS' ) 
      AWAS = 0d0

      ALLOCATE( AWLC( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWLC' ) 
      AWLC = 0d0

      ALLOCATE( AWSA( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWSA' ) 
      AWSA = 0d0
       
      ALLOCATE( AWSB( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWSB' ) 
      AWSB = 0d0  

      ALLOCATE( AWSC( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWSC' ) 
      AWSC = 0d0

      ALLOCATE( AWSS( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWSS' ) 
      AWSS = 0d0

      ALLOCATE( AWSN( NZSR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AWSN' ) 
      AWSN = 0d0

      ALLOCATE( BNC198( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC198' )
      BNC198 = 0e0

      ALLOCATE( BNC223( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC223' )
      BNC223 = 0e0

      ALLOCATE( BNC248( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC248' )
      BNC248 = 0e0

      ALLOCATE( BNC273( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC273' )
      BNC273 = 0e0
     
      ALLOCATE( BNC298( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC298' )
      BNC298 = 0e0

      ALLOCATE( BNC323( IMAX, JMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNC323' )
      BNC323 = 0e0

      ALLOCATE( ERRMSG( NERRMX ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'ERRSTK' )
      ERRMSG = ''

      ALLOCATE( ERRSTK( NERRMX ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'ERRSTK' )
      ERRSTK = 0d0

      ALLOCATE( GAMA( NPAIR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GAMA' ) 
      GAMA = 0d0

      ALLOCATE( GAMIN( NPAIR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GAMIN' ) 
      GAMIN = 0d0

      ALLOCATE( GAMOU( NPAIR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GAMOU' ) 
      GAMOU = 0d0

      ALLOCATE( GASAQ( NGASAQ ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GASAQ' ) 
      GASAQ = 0d0

      ALLOCATE( M0( NPAIR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'M0' ) 
      M0 = 0d0

      ALLOCATE( MOLAL( NIONS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOLAL' ) 
      MOLAL = 0d0

      ALLOCATE( MOLALR( NPAIR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOLALR' ) 
      MOLALR = 0d0

      ALLOCATE( IMW( NIONS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IMW' ) 
      IMW = 0d0

      ALLOCATE( SMW( NPAIR ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'WMW' )
      SMW = 0d0

      ALLOCATE( W( NCOMP ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'W' )
      W = 0d0

      ALLOCATE( WAER( NCOMP ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'W' )
      WAER = 0d0      

      ALLOCATE( WMW( NCOMP ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'WMW' )
      WMW = 0d0  

      ALLOCATE( Z( NIONS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Z' ) 
      Z = 0d0

      ALLOCATE( ZZ( NPAIR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZZ' ) 
      ZZ = 0d0

      !=================================================================
      ! Initialize module variables
      !=================================================================
      TEMP    = 298.0d0
      R       = 82.0567d-6
      RH      = 0.9d0
      EPS     = 1d-3
      MAXIT   = 100
      TINY    = 1d-20
      GREAT   = 1d10
      ZERO    = 0.0d0
      ONE     = 1.0d0
      NSWEEP  = 4
      TINY2   = 1d-11
      NDIV    = 5
      MOLAL   = 0.0d0 
      MOLALR  = 0.0d0 
      GAMA    = 0.1d0
      GAMOU   = 1d10
      GAMIN   = 1d10   
      CALAIN  = .TRUE.
      CALAOU  = .TRUE.
      EPSACT  = 5D-2        
      ICLACT  = 0
      IACALC  = 1          
      WFTYP   = 2 
      ERRSTK  = 0   
      ERRMSG  = ' '   
      NOFER   = 0 
      STKOFL  = .FALSE. 
      IPROB   = 0
      METSTBL = 0 
      VERSION = '1.5 (12/12/01)'


      SMW     = (/ 58.5,142.,85.0,132.,80.0,53.5,98.0,98.0,115.,63.0,
     &             36.5,120.,247./)

      IMW     = (/ 1.0,23.0,18.0,35.5,96.0,97.0,63.0/)
      WMW     = (/ 23.0,98.0,17.0,63.0,36.5/)
      ZZ      = (/1,2,1,2,1,1,2,1,1,1,1,1,2/)
      Z       = (/1,1,1,1,2,1,1/)

      !=================================================================
      ! ZSR RELATIONSHIP PARAMETERS
      !=================================================================

      ! AWAS = ammonium sulfate
      AWAS = (/
     & 100.,100.,100.,100.,100.,100.,100.,100.,100.,100.,100.,
     & 100.,100.,100.,100.,100.,100.,100.,100.,100.,100.,100.,
     & 100.,100.,100.,100.,100.,100.,100.,100.,100.,100.,100.,
     & 30.,30.,30.,29.,54.,28.25,27.06,25.94,
     & 24.89,23.90,22.97,22.10,21.27,20.48,19.73,19.02,18.34,17.69,
     & 17.07,16.48,15.91,15.37,14.85,14.34,13.86,13.39,12.94,12.50,
     & 12.08,11.67,11.27,10.88,10.51,10.14, 9.79, 9.44, 9.10, 8.78,
     &  8.45, 8.14, 7.83, 7.53, 7.23, 6.94, 6.65, 6.36, 6.08, 5.81,
     &  5.53, 5.26, 4.99, 4.72, 4.46, 4.19, 3.92, 3.65, 3.38, 3.11,
     &  2.83, 2.54, 2.25, 1.95, 1.63, 1.31, 0.97, 0.63, 0.30, 0.001/)

      ! AWSN= sodium nitrate
      AWSN = (/
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,685.59,
     & 451.00,336.46,268.48,223.41,191.28,
     & 167.20,148.46,133.44,121.12,110.83,
     & 102.09,94.57,88.03,82.29,77.20,72.65,68.56,64.87,61.51,58.44,
     & 55.62,53.03,50.63,48.40,46.32,44.39,42.57,40.87,39.27,37.76,
     & 36.33,34.98,33.70,32.48,31.32,30.21,29.16,28.14,27.18,26.25,
     & 25.35,24.50,23.67,22.87,22.11,21.36,20.65,19.95,19.28,18.62,
     & 17.99,17.37,16.77,16.18,15.61,15.05,14.51,13.98,13.45,12.94,
     & 12.44,11.94,11.46,10.98,10.51,10.04, 9.58, 9.12, 8.67, 8.22,
     &  7.77, 7.32, 6.88, 6.43, 5.98, 5.53, 5.07, 4.61, 4.15, 3.69,
     &  3.22, 2.76, 2.31, 1.87, 1.47, 1.10, 0.77, 0.48, 0.23, 0.001/)

      ! AWSC = sodium chloride
      AWSC =(/
     &  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,
     &  100., 100., 100., 100., 100., 100., 100., 100., 100.,16.34,
     & 16.28,16.22,16.15,16.09,16.02,15.95,15.88,15.80,15.72,15.64,
     & 15.55,15.45,15.36,15.25,15.14,15.02,14.89,14.75,14.60,14.43,
     & 14.25,14.04,13.81,13.55,13.25,12.92,12.56,12.19,11.82,11.47,
     & 11.13,10.82,10.53,10.26,10.00, 9.76, 9.53, 9.30, 9.09, 8.88,
     &  8.67, 8.48, 8.28, 8.09, 7.90, 7.72, 7.54, 7.36, 7.17, 6.99,
     &  6.81, 6.63, 6.45, 6.27, 6.09, 5.91, 5.72, 5.53, 5.34, 5.14,
     &  4.94, 4.74, 4.53, 4.31, 4.09, 3.86, 3.62, 3.37, 3.12, 2.85,
     &  2.58, 2.30, 2.01, 1.72, 1.44, 1.16, 0.89, 0.64, 0.40, 0.18/)

      ! AWAC = ammonium chloride
      AWAC = (/
     &  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,
     &  100., 100., 100., 100., 100., 100., 100., 100., 100.,31.45,
     & 31.30,31.14,30.98,30.82,30.65,30.48,30.30,30.11,29.92,29.71,
     & 29.50,29.29,29.06,28.82,28.57,28.30,28.03,27.78,27.78,27.77,
     & 27.77,27.43,27.07,26.67,26.21,25.73,25.18,24.56,23.84,23.01,
     & 22.05,20.97,19.85,18.77,17.78,16.89,16.10,15.39,14.74,14.14,
     & 13.59,13.06,12.56,12.09,11.65,11.22,10.81,10.42,10.03, 9.66,
     &  9.30, 8.94, 8.59, 8.25, 7.92, 7.59, 7.27, 6.95, 6.63, 6.32,
     &  6.01, 5.70, 5.39, 5.08, 4.78, 4.47, 4.17, 3.86, 3.56, 3.25,
     &  2.94, 2.62, 2.30, 1.98, 1.65, 1.32, 0.97, 0.62, 0.26, 0.13/)

      ! AWSS = sodium sulfate
      AWSS = (/
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,
     & 14.30,14.30,14.30,14.30,14.30,14.30,14.30,14.30,14.30,14.30,
     & 14.30,14.30,14.30,14.30,14.30,14.30,14.30,14.30,14.30,14.30,
     & 14.30,14.30,14.30,14.21,12.53,11.47,
     & 10.66,10.01, 9.46, 8.99, 8.57, 8.19, 7.85, 7.54, 7.25, 6.98,
     &  6.74, 6.50, 6.29, 6.08, 5.88, 5.70, 5.52, 5.36, 5.20, 5.04,
     &  4.90, 4.75, 4.54, 4.34, 4.14, 3.93, 3.71, 3.49, 3.26, 3.02,
     &  2.76, 2.49, 2.20, 1.89, 1.55, 1.18, 0.82, 0.49, 0.22, 0.001/)

      ! AWAB = ammonium bisulfate
      AWAB = (/356.45,296.51,253.21,220.47,194.85,
     & 174.24,157.31,143.16,131.15,120.82,
     & 111.86,103.99,97.04,90.86,85.31,80.31,75.78,71.66,67.90,64.44,
     &  61.25,58.31,55.58,53.04,50.68,48.47,46.40,44.46,42.63,40.91,
     &  39.29,37.75,36.30,34.92,33.61,32.36,31.18,30.04,28.96,27.93,
     &  26.94,25.99,25.08,24.21,23.37,22.57,21.79,21.05,20.32,19.63,
     &  18.96,18.31,17.68,17.07,16.49,15.92,15.36,14.83,14.31,13.80,
     &  13.31,12.83,12.36,11.91,11.46,11.03,10.61,10.20, 9.80, 9.41,
     &   9.02, 8.64, 8.28, 7.91, 7.56, 7.21, 6.87, 6.54, 6.21, 5.88,
     &   5.56, 5.25, 4.94, 4.63, 4.33, 4.03, 3.73, 3.44, 3.14, 2.85,
     &   2.57, 2.28, 1.99, 1.71, 1.42, 1.14, 0.86, 0.57, 0.29, 0.001/)

      ! AWSA = sulfuric acid
      AWSA = (/
     & 34.0,33.56,29.22,26.55,24.61,23.11,21.89,20.87,19.99,
     & 19.21,18.51,17.87,17.29,16.76,16.26,15.8,15.37,14.95,14.56,
     & 14.20,13.85,13.53,13.22,12.93,12.66,12.40,12.14,11.90,11.67,
     & 11.44,11.22,11.01,10.8,10.60,10.4,10.2,10.01,9.83,9.65,9.47,
     & 9.3,9.13,8.96,8.81,8.64,8.48,8.33,8.17,8.02,7.87,7.72,7.58,
     & 7.44,7.30,7.16,7.02,6.88,6.75,6.61,6.48,6.35,6.21,6.08,5.95,
     & 5.82,5.69,5.56,5.44,5.31,5.18,5.05,4.92,4.79,4.66,4.53,4.40,
     & 4.27,4.14,4.,3.87,3.73,3.6,3.46,3.31,3.17,3.02,2.87,2.72,
     & 2.56,2.4,2.23,2.05,1.87,1.68,1.48,1.27,1.05,0.807,0.552,0.281/)

      ! AWLC = (NH4)3H(SO4)2
      AWLC = (/
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,
     & 17.0,16.5,15.94,15.31,14.71,14.14,
     & 13.60,13.08,12.59,12.12,11.68,11.25,10.84,10.44,10.07, 9.71,
     &  9.36, 9.02, 8.70, 8.39, 8.09, 7.80, 7.52, 7.25, 6.99, 6.73,
     &  6.49, 6.25, 6.02, 5.79, 5.57, 5.36, 5.15, 4.95, 4.76, 4.56,
     &  4.38, 4.20, 4.02, 3.84, 3.67, 3.51, 3.34, 3.18, 3.02, 2.87,
     &  2.72, 2.57, 2.42, 2.28, 2.13, 1.99, 1.85, 1.71, 1.57, 1.43,
     &  1.30, 1.16, 1.02, 0.89, 0.75, 0.61, 0.46, 0.32, 0.16, 0.001/)

      ! AWAN = ammonium nitrate
      AWAN = (/
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,
     & 1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,1.e5,
     &       97.17,92.28,87.66,83.15,78.87,74.84,70.98,67.46,64.11,
     & 60.98,58.07,55.37,52.85,50.43,48.24,46.19,44.26,42.40,40.70,
     & 39.10,37.54,36.10,34.69,33.35,32.11,30.89,29.71,28.58,27.46,
     & 26.42,25.37,24.33,23.89,22.42,21.48,20.56,19.65,18.76,17.91,
     & 17.05,16.23,15.40,14.61,13.82,13.03,12.30,11.55,10.83,10.14,
     &  9.44, 8.79, 8.13, 7.51, 6.91, 6.32, 5.75, 5.18, 4.65, 4.14,
     &  3.65, 3.16, 2.71, 2.26, 1.83, 1.42, 1.03, 0.66, 0.30, 0.001/)

      ! AWSB = sodium bisulfate
      AWSB = (/ 173.72,156.88,142.80,130.85,120.57,
     &          111.64,103.80,96.88,90.71,85.18,
     & 80.20,75.69,71.58,67.82,64.37,61.19,58.26,55.53,53.00,50.64,
     & 48.44,46.37,44.44,42.61,40.90,39.27,37.74,36.29,34.91,33.61,
     & 32.36,31.18,30.05,28.97,27.94,26.95,26.00,25.10,24.23,23.39,
     & 22.59,21.81,21.07,20.35,19.65,18.98,18.34,17.71,17.11,16.52,
     & 15.95,15.40,14.87,14.35,13.85,13.36,12.88,12.42,11.97,11.53,
     & 11.10,10.69,10.28, 9.88, 9.49, 9.12, 8.75, 8.38, 8.03, 7.68,
     &  7.34, 7.01, 6.69, 6.37, 6.06, 5.75, 5.45, 5.15, 4.86, 4.58,
     &  4.30, 4.02, 3.76, 3.49, 3.23, 2.98, 2.73, 2.48, 2.24, 2.01,
     &  1.78, 1.56, 1.34, 1.13, 0.92, 0.73, 0.53, 0.35, 0.17, 0.001/)

      ASSO4 = (/ 1.0E-9, 2.5E-9, 5.0E-9, 7.5E-9, 1.0E-8,
     &           2.5E-8, 5.0E-8, 7.5E-8, 1.0E-7, 2.5E-7, 
     &           5.0E-7, 7.5E-7, 1.0E-6, 5.0E-6 /)

      ASRAT = (/
     & 1.020464, 0.9998130, 0.9960167, 0.9984423, 1.004004, 1.010885,  
     & 1.018356, 1.026726,  1.034268,  1.043846,  1.052933, 1.062230,  
     & 1.062213, 1.080050,  1.088350,  1.096603,  1.104289, 1.111745,  
     & 1.094662, 1.121594,  1.268909,  1.242444,  1.233815, 1.232088,
     & 1.234020, 1.238068,  1.243455,  1.250636,  1.258734, 1.267543,
     & 1.276948, 1.286642,  1.293337,  1.305592,  1.314726, 1.323463,  
     & 1.333258, 1.343604,  1.344793,  1.355571,  1.431463, 1.405204,
     & 1.395791, 1.393190,  1.394403,  1.398107,  1.403811, 1.411744,  
     & 1.420560, 1.429990,  1.439742,  1.449507,  1.458986, 1.468403, 
     & 1.477394, 1.487373,  1.495385,  1.503854,  1.512281, 1.520394,
     & 1.514464, 1.489699,  1.480686,  1.478187,  1.479446, 1.483310, 
     & 1.489316, 1.497517,  1.506501,  1.515816,  1.524724, 1.533950,
     & 1.542758, 1.551730,  1.559587,  1.568343,  1.575610, 1.583140,  
     & 1.590440, 1.596481,  1.567743,  1.544426,  1.535928, 1.533645, 
     & 1.535016, 1.539003,  1.545124,  1.553283,  1.561886, 1.570530,
     & 1.579234, 1.587813,  1.595956,  1.603901,  1.611349, 1.618833, 
     & 1.625819, 1.632543,  1.639032,  1.645276,  1.707390, 1.689553,  
     & 1.683198, 1.681810,  1.683490,  1.687477,  1.693148, 1.700084,  
     & 1.706917, 1.713507,  1.719952,  1.726190,  1.731985, 1.737544, 
     & 1.742673, 1.747756,  1.752431,  1.756890,  1.761141, 1.765190,
     & 1.785657, 1.771851,  1.767063,  1.766229,  1.767901, 1.771455, 
     & 1.776223, 1.781769,  1.787065,  1.792081,  1.796922, 1.801561,  
     & 1.805832, 1.809896,  1.813622,  1.817292,  1.820651, 1.823841,  
     & 1.826871, 1.829745,  1.822215,  1.810497,  1.806496, 1.805898, 
     & 1.807480, 1.810684,  1.814860,  1.819613,  1.824093, 1.828306,
     & 1.832352, 1.836209,  1.839748,  1.843105,  1.846175, 1.849192, 
     & 1.851948, 1.854574,  1.857038,  1.859387,  1.844588, 1.834208,  
     & 1.830701, 1.830233,  1.831727,  1.834665,  1.838429, 1.842658,
     & 1.846615, 1.850321,  1.853869,  1.857243,  1.860332, 1.863257, 
     & 1.865928, 1.868550,  1.870942,  1.873208,  1.875355, 1.877389,
     & 1.899556, 1.892637,  1.890367,  1.890165,  1.891317, 1.893436, 
     & 1.896036, 1.898872,  1.901485,  1.903908,  1.906212, 1.908391,  
     & 1.910375, 1.912248,  1.913952,  1.915621,  1.917140, 1.918576,  
     & 1.919934, 1.921220,  1.928264,  1.923245,  1.921625, 1.921523, 
     & 1.922421, 1.924016,  1.925931,  1.927991,  1.929875, 1.931614,
     & 1.933262, 1.934816,  1.936229,  1.937560,  1.938769, 1.939951, 
     & 1.941026, 1.942042,  1.943003,  1.943911,  1.941205, 1.937060,  
     & 1.935734, 1.935666,  1.936430,  1.937769,  1.939359, 1.941061,
     & 1.942612, 1.944041,  1.945393,  1.946666,  1.947823, 1.948911, 
     & 1.949900, 1.950866,  1.951744,  1.952574,  1.953358, 1.954099,
     & 1.948985, 1.945372,  1.944221,  1.944171,  1.944850, 1.946027, 
     & 1.947419, 1.948902,  1.950251,  1.951494,  1.952668, 1.953773,  
     & 1.954776, 1.955719,  1.956576,  1.957413,  1.958174, 1.958892,  
     & 1.959571, 1.960213,  1.977193,  1.975540,  1.975023, 1.975015, 
     & 1.975346, 1.975903,  1.976547,  1.977225,  1.977838, 1.978401,
     & 1.978930, 1.979428,  1.979879,  1.980302,  1.980686, 1.981060, 
     & 1.981401, 1.981722,  1.982025,  1.982312 /)

      ! Return to calling program
      END SUBROUTINE INIT_ISOROPIA 

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_ISOROPIA
!
!******************************************************************************
!  Subroutine CLEANUP_KMC deallocates all module arrays. (rjp, bmy, 9/23/02)
!******************************************************************************
!
      IF ( ALLOCATED( BNC198 ) ) DEALLOCATE( BNC198 )
      IF ( ALLOCATED( BNC223 ) ) DEALLOCATE( BNC223 )
      IF ( ALLOCATED( BNC248 ) ) DEALLOCATE( BNC248 )
      IF ( ALLOCATED( BNC273 ) ) DEALLOCATE( BNC273 )
      IF ( ALLOCATED( BNC298 ) ) DEALLOCATE( BNC298 )
      IF ( ALLOCATED( BNC323 ) ) DEALLOCATE( BNC323 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_ISOROPIA

!------------------------------------------------------------------------------

      END MODULE ISOROPIA_MOD
