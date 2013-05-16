!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: gmistrat_mod
!     
! !DESCRIPTION: Module GMISTRAT\_MOD contains routines and variables which
!  are associated with the addition of full stratospheric chemistry to
!  GEOS-Chem (based on the NASA GMI implementation
!\\   
!\\   
! !INTERFACE: 
!
      MODULE GMISTRAT_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN
      USE ERROR_MOD,  ONLY : DEBUG_MSG
      USE LOGICAL_MOD,ONLY : LPRT
      USE CMN_SIZE_MOD

      IMPLICIT NONE
#     include "define.h"
      PRIVATE

!
! !PUBLIC DATA MEMBERS:
!
      CHARACTER(LEN=255), PUBLIC :: TRAC_FILE_ROOT ! Root directory of 2D data
      CHARACTER(LEN=255), PUBLIC :: MONTREAL_FILE_ROOT ! Root directoy for WHO data 
                          PUBLIC :: T_STS ! Max temperature of STS formation (K)
      REAL*8,             PUBLIC :: T_NAT_SUPERCOOL ! NAT supercooling (K)
      REAL*8,             PUBLIC :: P_ICE_SUPERSAT ! Ice supersaturation (-)
                          PUBLIC :: KHETI_SLA ! Reaction cofactors for SLA
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: SET_INITIAL_MIXRATIOS
      PUBLIC  :: EMISS_BASIC
      PUBLIC  :: SETTLE_STRAT_AER
      PUBLIC  :: SO4_PHOTFRAC
      PUBLIC  :: CALC_STRAT_AER
      PUBLIC  :: GET_STRAT_OPT
      PUBLIC  :: KG_STRAT_AER
      PUBLIC  :: RHO_STRAT_AER
      PUBLIC  :: CHECK_NAT
      PUBLIC  :: MAKE_PSC_FILE
      PUBLIC  :: INIT_GMISTRAT
      PUBLIC  :: CLEANUP_GMISTRAT
!
! PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: APPLY_2DTRAC
      PRIVATE :: READ_SFC
      PRIVATE :: GET_MONTREAL
      PRIVATE :: TERNARY
      PRIVATE :: READ_PSC_FILE
      PRIVATE :: CARSLAW_DENSITY
      PRIVATE :: CALC_H2SO4_GAS
      PRIVATE :: CALC_SLA_GAMMA
      PRIVATE :: MOLEC_SPEED
!
! !REVISION HISTORY: 
!  26 Mar 2013 - S. D. Eastham - Initial version
!  04 Apr 2013 - S. D. Eastham - Rolled several routines into module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! MODULE PARAMETERS
      !
      ! AER_NLEVS       : Number of levels in AER data
      ! T_STS           : Maximum temperature of STS formation (K)
      ! I_SLA           : Index of liquid aerosols
      ! I_SPA           : Index of particulate PSCs
      ! AVOGADRO        : Avogadro's number (#/mole)
      !
      !=================================================================

      INTEGER, PARAMETER            :: AER_NLEVS=51
      REAL*8, PARAMETER             :: T_STS=240.0d0
      INTEGER, PARAMETER            :: I_SLA=1
      INTEGER, PARAMETER            :: I_SPA=2
      REAL*8,  PARAMETER            :: AVOGADRO=6.022d23

!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Scalars
      !
      ! TRAC_IDX           : Tracer index for output
      ! N_EMIT             : Number of tracers to emit
      ! MAX_EMIT           : Limit on tracer emission count
      ! SFC_MONTH          : Current month
      ! SLA_VA             : SLA volume-area conversion
      ! SLA_RR             : SLA effective-liquid radius conversion
      ! SLA_VR             : SLA volume-effective radius conversion
      ! NATMW              : Molar mass of NAT (kg/kmol)
      ! ICEMEW             : Molar mass of ice (kg/kmol)
      ! DENSNAT            : Density of pure NAT (kg/m3)
      ! DENSICE            : Density of pure ice (kg/m3)
      ! ISR_ClNO3          : ClNO3 MW (inverse sqrt) (kg/kmol)^-0.5 
      ! ISR_BrNO3          : BrNO3 MW (inverse sqrt) (kg/kmol)^-0.5 
      ! ISR_N2O5           : N2O5 MW (inverse sqrt) (kg/kmol)^-0.5
      ! ISR_HOCl           : HOCl MW (inverse sqrt) (kg/kmol)^-0.5
      ! ISR_HOBr           : HOBr MW (inverse sqrt) (kg/kmol)^-0.5
      !
      ! Arrays
      !
      ! AER_MR             : AER 2D mixing ratios (v/v)
      ! AER_PLEVS          : Pressure levels of AER data (hPa)
      ! IDT_EMIT           : IDs of emitted tracers
      ! GRID_EMIT          : Surface CFC mixing ratios
      ! SFC_MONTH          : Month of last read-in
      ! RAD_AER            : Strat. aerosol radius (cm)
      ! KG_AER             : Aerosol mass (kg/box)
      ! SAD_AER            : Aerosol surface area density (cm2/cm3)
      ! NDENS_AER          : Aerosol number density (#/m3)
      ! RHO_AER            : Aerosol mass density (kg/m3 aerosol)
      ! STATE_PSC          : Gridbox PSC type (see Kirner et al)
      ! KHETI_SLA          : Liquid aerosol reaction cofactors
      ! AERFRAC            : Mass fraction of species in liquid aerosols
      ! AERFRACIND         : Indices of liquid aerosol species
      !
      ! Strings
      !
      ! TRAC_FILE          : Input filename
      ! MONTH_STR          : Current month
      ! TRAC_SET           : Tracer name to be overwritten
      ! TRAC_EMIT          ! Tracer names for emission
      !
      ! Logicals
      ! 
      ! TRAC_ADD            : Add to (rather than overwrite) tracer
      ! 
      !=================================================================

      ! Scalars
      INTEGER                           :: TRAC_IDX
      INTEGER                           :: N_EMIT
      INTEGER, PARAMETER                :: MAX_EMIT=50
      INTEGER                           :: SFC_MONTH
      REAL*8                            :: SLA_VA
      REAL*8                            :: SLA_RR
      REAL*8                            :: SLA_VR
      REAL*8, PARAMETER                 :: NATMW   = 117.0
      REAL*8, PARAMETER                 :: ICEMW   = 18.0
      REAL*8, PARAMETER                 :: DENSNAT = 1626.d0
      REAL*8, PARAMETER                 :: DENSICE = 990.0d0
      REAL*8, PARAMETER                 :: ISR_ClNO3=1.d0/sqrt(97.46d0)
      REAL*8, PARAMETER                 :: ISR_BrNO3=1.d0/sqrt(141.9d0)
      REAL*8, PARAMETER                 :: ISR_N2O5 =1.d0/sqrt(108.0d0)
      REAL*8, PARAMETER                 :: ISR_HOCl =1.d0/sqrt(52.46d0)
      REAL*8, PARAMETER                 :: ISR_HOBr =1.d0/sqrt(96.91d0)

      ! Arrays
      REAL*8,DIMENSION(:,:),ALLOCATABLE     :: AER_MR
      REAL*8,DIMENSION(:),ALLOCATABLE       :: AER_PLEVS
      INTEGER,DIMENSION(MAX_EMIT)           :: IDT_EMIT
      REAL*8,DIMENSION(:,:),ALLOCATABLE     :: GRID_EMIT
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: RAD_AER
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: KG_AER
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: SAD_AER
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: NDENS_AER
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: RHO_AER
      INTEGER,DIMENSION(:,:,:),ALLOCATABLE  :: STATE_PSC
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: KHETI_SLA
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: AERFRAC
      INTEGER,DIMENSION(:),ALLOCATABLE      :: AERFRACIND
 
      ! Strings
      CHARACTER(LEN=255)                :: TRAC_FILE
      CHARACTER(LEN=2)                  :: MONTH_STR
      CHARACTER(LEN=255)                :: TRAC_SET 
      CHARACTER*20,DIMENSION(50)        :: TRAC_EMIT

      ! Logicals
      LOGICAL                           :: TRAC_ADD

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_initial_mixratios
!
! !DESCRIPTION: Subroutine SET\_INITIAL\_MIXRATIOS is a public interface.
!  Concentrations of each of the following tracers are read from
!  2D estimates made using the AER 2D model:
!                                                                             .
!     (# ) TRC                 = Constituents         < Controlled by
! ----------------------------------------------------------------------------
!     (1 ) CH4                 = CH4                  < LSETCH4
!     (2 ) N2O                 = N2O                  < LSETN2O
!     (3 ) OCS                 = OCS                  < LSETOCS
!     (4 ) SO4                 = H2SO4                < LSETH2SO4
!     (5 ) CFC                 = CFC-113/114/115      < LSETCFC
!     (6 ) HCFC                = HCFC-22/141b/142b    < LSETCFC
!     (7 ) CFC11               = CFC-11               < LSETCFC
!     (8 ) CFC12               = CFC-12               < LSETCFC
!     (9 ) H1202               = Halon 1202           < LSETCFC
!     (10) H1211               = Halon 1211           < LSETCFC
!     (11) H1301               = Halon 1301           < LSETCFC
!     (12) H2402               = Halon 2402           < LSETCFC
!     (13) Cl2                 = Cl2                  < LSETCL
!     (14) ClOx                = Cl + ClO             < LSETCL
!     (15) CCl4                = CCl4                 < LSETCL
!     (16) CH3Cl               = CH3Cl                < LSETCL
!     (17) CH3CCl3             = CH3CCl3              < LSETCL
!     (18) HCl                 = HCl                  < LSETCL
!     (19) HOCl                = HOCl                 < LSETCL
!     (20) Cl2O2               = Cl2O2                < LSETCL
!     (21) ClNO2               = ClNO2                < LSETCL
!     (22) ClONO2              = ClNO3                < LSETCL
!     (23) OClO                = OClO                 < LSETCL
!     (24) ClOO                = ClOO                 < LSETCL
!     (25) BrCl                = BrCl                 < (LSETCL || LSETBR || LSETBRSTRAT)
!     (26) Br2                 = Br2                  < LSETBR || LSETBRSTRAT
!     (27) Br                  = Br                   < LSETBR || LSETBRSTRAT 
!     (28) BrO                 = BrO                  < LSETBR || LSETBRSTRAT
!     (29) HOBr                = HOBr                 < LSETBR || LSETBRSTRAT
!     (30) HBr                 = HBr                  < LSETBR || LSETBRSTRAT
!     (31) BrNO2               = BrNO2                < LSETBR || LSETBRSTRAT
!     (32) BrNO3               = BrNO3                < LSETBR || LSETBRSTRAT
!     (33) CHBr3               = CHBr3                < LSETBR || LSETBRSTRAT
!     (34) CH2Br2              = CH2Br2               < LSETBR || LSETBRSTRAT
!     (35) CH3Br               = CH3Br                < LSETBR || LSETBRSTRAT
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_INITIAL_MIXRATIOS( am_I_Root )
!
! !USES:
!
      USE LOGICAL_MOD,          ONLY : LSETCH4,  LSETOCS,  LSETCFC
      USE LOGICAL_MOD,          ONLY : LSETCL,   LSETBR,   LSETN2O
      USE LOGICAL_MOD,          ONLY : LSETH2SO4,LGMISTRAT, LSETBRSTRAT
      USE TRACERID_MOD
      USE TRACER_MOD,           ONLY : N_TRACERS, STT,TCVV
      USE DAO_MOD,              ONLY : AD
      USE TIME_MOD,             ONLY : GET_MONTH
      USE ERROR_MOD,            ONLY : ALLOC_ERR,ERROR_STOP
      USE GRID_MOD,             ONLY : GET_YMID
      USE TIME_MOD,             ONLY : GET_YEAR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_YEAR
      USE LOGICAL_MOD,          ONLY : LFUTURE
!
! !REVISION HISTORY: 
!  26 Mar 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: AS, N, I, J, L
      REAL*8             :: PLAST, PCURR, PNEXT
      REAL*8             :: YLAT
      CHARACTER(LEN=255) :: MSG
      LOGICAL            :: USE2DDATA
      REAL*8             :: C3090S, C0030S, C0030N, C3090N
      INTEGER            :: CH4_YEAR
      LOGICAL,INTENT(IN) :: am_I_Root
      LOGICAL,PARAMETER  :: STRAT2DCH4=.TRUE.

      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100   )
      WRITE( 6, 110   )
      WRITE( 6, 120   )
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! FORMAT strings
 100  FORMAT( 'T R A C E R   I N I T I A L I Z A T I O N'   )
 110  FORMAT( 'Routine written by SEBASTIAN D. EASTHAM' )
 120  FORMAT( 'Last modified: 03/26/13' )

      ! Get current month string
      WRITE( MONTH_STR, '(i0.2)' ) GET_MONTH()

      ! Are we using 2D data?
      USE2DDATA = (LSETOCS.or.LSETCFC.or.LSETCL.or.LSETBR.or.LSETH2SO4)
      USE2DDATA = (USE2DDATA.or.STRAT2DCH4.or.LSETBRSTRAT)
      USE2DDATA = (USE2DDATA.or.LSETN2O)

#if !defined(GRID4x5) && !defined(GRID2x25)
      IF (LGMISTRAT.and.USE2DDATA) THEN
       MSG = 'Zonal means preprocessed only for 2x2.5 and 4x5 grid.'
       CALL ERROR_STOP( MSG, 'SET_INITIAL_MIXRATIOS (gmistrat_mod.f)!' )
      ELSE
       TRAC_FILE_ROOT = ''
      ENDIF
#elif defined(GRID2x25)
      TRAC_FILE_ROOT = TRIM(TRAC_FILE_ROOT)
     &      // '/Grid2x25/InitCFC_'
#elif defined(GRID4x5)
      TRAC_FILE_ROOT = TRIM(TRAC_FILE_ROOT)
     &      // '/Grid4x5/InitCFC_'
#endif

      IF (USE2DDATA) THEN
         ! Allocate array of input pressure levels
         ALLOCATE( AER_PLEVS( AER_NLEVS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AER_PLEVS' )

         ! Set input pressure levels (hPa)
         AER_PLEVS = (/ 0.2200d+00, 0.2600d+00, 0.3100d+00, 0.3600d+00,
     &                  0.4300d+00, 0.5100d+00, 0.6000d+00, 0.7100d+00,
     &                  0.8400d+00, 0.9900d+00, 1.1700d+00, 1.3800d+00,
     &                  1.6300d+00, 1.9300d+00, 2.2800d+00, 2.6900d+00,
     &                  3.1800d+00, 3.7600d+00, 4.4400d+00, 5.2500d+00,
     &                  6.2000d+00, 7.3200d+00, 8.6500d+00, 1.0220d+01,
     &                  1.2070d+01, 1.4260d+01, 1.6850d+01, 1.9910d+01,
     &                  2.3520d+01, 2.7780d+01, 3.2820d+01, 3.8770d+01,
     &                  4.5810d+01, 5.4110d+01, 6.3930d+01, 7.5220d+01,
     &                  8.9220d+01, 1.0540d+02, 1.2451d+02, 1.4710d+02,
     &                  1.7377d+02, 2.0529d+02, 2.4252d+02, 2.8650d+02,
     &                  3.3847d+02, 3.9985d+02, 4.7237d+02, 5.5804d+02,
     &                  6.5924d+02, 7.7880d+02, 9.2004d+02 /)
 
         ALLOCATE( AER_MR( IGLOB, AER_NLEVS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AER_MR' )
         AER_MR = 0d0

         IF ( am_I_Root ) THEN
            WRITE(6,*) 'Setting initial values for:'
            WRITE(6,*) '|---------|---------|---------|'
         ENDIF

         ! Bromine species (except Halons)
         IF (LSETBR.or.LSETBRSTRAT) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|Br2      |Br       |BrO      |'
               WRITE(6,*) '|HOBr     |HBr      |BrNO2    |'
               WRITE(6,*) '|BrONO2   |CHBr3    |CH2Br2   |'
               WRITE(6,*) '|CH3Br    |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTBr2
            TRAC_SET  = 'BR2'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
         
            TRAC_IDX = IDTBr
            TRAC_SET  = 'BR'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
         
            TRAC_IDX = IDTBrO
            TRAC_SET  = 'BRO'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
         
            TRAC_IDX = IDTHOBr
            TRAC_SET  = 'HOBR'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
         
            TRAC_IDX = IDTHBr
            TRAC_SET  = 'HBR'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
     
            TRAC_IDX = IDTBrNO2
            TRAC_SET  = 'BRNO2'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
     
            TRAC_IDX = IDTBrNO3
            TRAC_SET  = 'BRNO3'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
     
            TRAC_IDX = IDTCHBr3 
            TRAC_SET  = 'CHBR3'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
     
            TRAC_IDX = IDTCH2Br2
            TRAC_SET  = 'CH2BR2'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
     
            TRAC_IDX = IDTCH3Br
            TRAC_SET  = 'CH3BR'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
         ENDIF

         ! Br-Cl (special case)
         IF (LSETCL.or.LSETBR.or.LSETBRSTRAT) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|BrCl     |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTBrCl
            TRAC_SET  = 'BRCL'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.)
     
         ENDIF

         ! Chlorine species (except CFCS, HCFCs and Halons)
         IF (LSETCL) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|Cl2      |ClO      |Cl2O2    |'
               WRITE(6,*) '|Cl       |ClONO2   |HCl      |'
               WRITE(6,*) '|HOCl     |OClO     |CH3Cl    |'
               WRITE(6,*) '|CH3CCl3  |CCl4     |ClOO     |'
               WRITE(6,*) '|         |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTCl
            TRAC_SET  = 'CL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTClO
            TRAC_SET  = 'CLO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTClNO2
            TRAC_SET  = 'CLNO2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTClNO3
            TRAC_SET  = 'CLNO3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTHOCl
            TRAC_SET  = 'HOCL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTOClO
            TRAC_SET  = 'OCLO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTClOO
            TRAC_SET  = 'CLOO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCl2
            TRAC_SET  = 'CL2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTCl2O2
            TRAC_SET  = 'CL2O2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTHCl
            TRAC_SET  = 'HCL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
     
            ! Handle chloroalkanes      
            TRAC_IDX = IDTCH3Cl
            TRAC_SET  = 'CH3CL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCH3CCl3
            TRAC_SET  = 'CH3CCL3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCCl4
            TRAC_SET  = 'CCL4'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
         ENDIf
         
         ! CFCs, HCFCs and Halons
         IF (LSETCFC) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|H-1211   |H-1301   |H-2402   |'
               WRITE(6,*) '|HCFC-123 |HCFC-141b|HCFC-142b|'
               WRITE(6,*) '|HCFC-22  |R-11     |R-12     |'
               WRITE(6,*) '|R-113    |R-114    |R-115    |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            ! H-1211, 1301 and 2402
            TRAC_IDX = IDTH1211
            TRAC_SET  = 'CBRCLF2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTH1301
            TRAC_SET  = 'CBRF3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTH2402
            TRAC_SET  = 'C2BR2F4'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
     
            ! CFC  11, 12, and [113/114/115]
            TRAC_IDX = IDTCFC11
            TRAC_SET  = 'CFCL3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTCFC12
            TRAC_SET  = 'CF2CL2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTCFCX
            TRAC_SET  = 'C2CL3F3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_SET  = 'C2CL2F4'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.)
            TRAC_SET  = 'C2CLF5'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.)
     
            ! HCFC 22, 123, 141b and 142b
            TRAC_IDX = IDTHCFC22
            TRAC_SET  = 'CHCLF2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_IDX = IDTHCFCX
            TRAC_SET  = 'C2HF3CL2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
            TRAC_SET  = 'C2H3FCL2'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.)
            TRAC_SET  = 'C2H3F2CL'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.)
         ENDIF

         IF (LSETOCS) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|OCS      |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTOCS
            TRAC_SET  = 'OCS'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
         ENDIF

         IF (LSETH2SO4) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|H2SO4    |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTSO4
            ! Actually H2SO4 aerosol
            TRAC_SET  = 'SO4'
            CALL APPLY_2DTRAC(.FALSE.,.TRUE.,.FALSE.)
         ENDIF
      ENDIF ! USE2DDATA

      IF (LSETN2O) THEN
         IF ( am_I_Root ) THEN
            WRITE(6,*) '|N2O      |         |         |'
            WRITE(6,*) '|---------|---------|---------|'
            IF ( LPRT ) CALL DEBUG_MSG( '### GMISTRAT: Forcing 2D N2O' )
         ENDIF

         TRAC_IDX = IDTN2O
         TRAC_SET  = 'N2O'
         CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.)
      ENDIF

      ! CH4 taken from the existing GEOS-Chem 3D input data,
      ! rather than from 2D data

      IF (LSETCH4) THEN
         IF ( am_I_Root ) THEN
            WRITE(6,*) '|CH4      |         |         |'
            WRITE(6,*) '|---------|---------|---------|'
         ENDIF

         ! Adapted from GEOS-Chem routine

         ! Pick the CH4 concentration [ppbv] for the proper lat bin
         ! CH4 values are read in "chemdr.f"
         IF ( LFUTURE ) THEN
            CH4_YEAR = GET_FUTURE_YEAR()
         ELSE
            CH4_YEAR = GET_YEAR()
         ENDIF

         ! Get CH4 [ppbv] in 4 latitude bins for each year
         CALL GET_GLOBAL_CH4( CH4_YEAR, .TRUE., C3090S,
     &                        C0030S,   C0030N, C3090N, .FALSE. )

         ! Convert from [ppbv] to [kg CH4/mol air]
         C3090S = C3090S * (1.d-9) / TCVV(IDTCH4)
         C0030S = C0030S * (1.d-9) / TCVV(IDTCH4)
         C0030N = C0030N * (1.d-9) / TCVV(IDTCH4)
         C3090N = C3090N * (1.d-9) / TCVV(IDTCH4)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, YLAT )
         DO I=1,IIPAR   
         DO J=1,JJPAR   
         DO L=1,LLPAR   
            ! Grid-box latitude index
            YLAT  = GET_YMID( I,J,L )
            IF ( YLAT < -30d0 ) THEN
               STT(I,J,L,IDTCH4) = C3090S*AD(I,J,L)
            ELSE IF ( YLAT >= -30d0 .and. YLAT < 0d0  ) THEN
               STT(I,J,L,IDTCH4) = C0030S*AD(I,J,L)
            ELSE IF ( YLAT >=   0d0 .and. YLAT < 30d0 ) THEN
               STT(I,J,L,IDTCH4) = C0030N*AD(I,J,L)
            ELSE
               STT(I,J,L,IDTCH4) = C3090N*AD(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Set stratospheric CH4 from AER data?
         IF ( STRAT2DCH4 ) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'Setting stratospheric CH4 from 2D data.'
            ENDIF
            ! Code to set strat using AER 2D data
            TRAC_IDX = IDTCH4
            TRAC_SET  = 'CH4'
            CALL APPLY_2DTRAC(.FALSE.,.TRUE.,.FALSE.)
         ENDIF
      ENDIF

      ! Clean up
      IF ( ALLOCATED( AER_MR    ) ) DEALLOCATE( AER_MR    )
      IF ( ALLOCATED( AER_PLEVS ) ) DEALLOCATE( AER_PLEVS )
 
      END SUBROUTINE SET_INITIAL_MIXRATIOS
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apply_2dtrac
!
! !DESCRIPTION: Subroutine APPLY\_2DTRAC reads in and applies 2D data 
!  from AER 2D model output to the 3D GEOS-Chem grid. Zonally-averaged
!  mixing ratios vary by month, but no interpolation is performed.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APPLY_2DTRAC( TRAC_ADD, STRAT_ONLY, SFC_ONLY )
!
! !USES:
!
      USE FILE_MOD,      ONLY : IOERROR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE PRESSURE_MOD,  ONLY : GET_PCENTER
      USE TRACER_MOD,    ONLY : STT, TCVV
      USE CHEMGRID_MOD,  ONLY : ITS_IN_THE_STRATMESO
      USE DAO_MOD,       ONLY : AD
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)              :: TRAC_ADD     ! Add (instead of overwrite) 
      LOGICAL, INTENT(IN)              :: STRAT_ONLY   ! Only write to stratosphere 
      LOGICAL, INTENT(IN)              :: SFC_ONLY     ! Only write to surface 
!
! !REVISION HISTORY: 
!  26 Mar 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: IOS
      INTEGER :: ILEV
      INTEGER :: I, J, L
      REAL*8  :: PCENTER
      REAL*8  :: TRACMR,TRACMASS
      INTEGER :: VERTCOUNT
      LOGICAL :: FOUNDLEV,EXTRAP
      REAL*8  :: TRACMULT
      LOGICAL :: SET_CELL
      INTEGER :: IU_FILE
      INTEGER :: L_MAX

      !=================================================================
      ! APPLY_2DTRAC begins here!
      !=================================================================

      TRAC_FILE = TRIM(TRAC_FILE_ROOT)//TRIM(TRAC_SET )
      TRAC_FILE = TRIM(TRAC_FILE)//'_'//TRIM(MONTH_STR)//'.dat'

      ! Error check
      IF ((STRAT_ONLY).and.(SFC_ONLY)) THEN
         CALL ERROR_STOP('Cannot have both SFC_ONLY and STRAT_ONLY!',
     &                           'APPLY_2DTRAC(gmistrat_mod.F)')
      ENDIF

      ! Get a free LUN
      IU_FILE = findFreeLUN()

      IOS = 1
      OPEN( IU_FILE, FILE=TRIM( TRAC_FILE ),STATUS='OLD',IOSTAT=IOS) 
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) 'GMISTRAT_MOD: Could not read ', TRIM( TRAC_FILE )
         CALL IOERROR( IOS, IU_FILE,'GMISTRAT_MOD:APPLY_2DTRAC')
      ENDIF
   
      ! Read in mixing ratios
      DO ILEV = 1,51
         READ(IU_FILE, 110, IOSTAT=IOS ) AER_MR(:,ILEV)
      ENDDO

#if defined(GRID2x25)
110   FORMAT(46E10.3)
#elif defined(GRID4x5)
110   FORMAT(91E10.3)
#endif

      CLOSE(IU_FILE)    
 
      IF (SFC_ONLY) THEN
         L_MAX = 1
      ELSE
         L_MAX = LLPAR
      ENDIF

      ! Scan through STT, element by element
      DO J = 1,JJPAR
         DO I = 1,IIPAR
            ! Vertcount is the layer count for the input, where layer 1
            ! is at the *top* of the atmosphere
            VERTCOUNT = 51
            EXTRAP = .TRUE.
            DO L = 1, L_MAX
               IF (STRAT_ONLY) THEN
                  ! Only want to set non-trop cells
                  SET_CELL = ITS_IN_THE_STRATMESO(I,J,L)
               ELSE
                  SET_CELL = .TRUE.
               ENDIF
               ! Pressure at center of cell
               PCENTER = GET_PCENTER(I,J,L)
               FOUNDLEV = (PCENTER.gt.AER_PLEVS(VERTCOUNT))
               DO WHILE (.not. FOUNDLEV)
                  IF (VERTCOUNT.eq.1) THEN
                    ! At top layer; use it anyway
                    FOUNDLEV = .TRUE.
                    EXTRAP = .TRUE.
                  ELSE
                    VERTCOUNT = VERTCOUNT - 1
                    FOUNDLEV = (PCENTER.gt.AER_PLEVS(VERTCOUNT))
                    EXTRAP = .FALSE.
                  ENDIF
               ENDDO

               ! Interpolate mixing ratio
               IF (EXTRAP) THEN
                  ! Just take outside value if at edges
                  TRACMR = AER_MR(J,VERTCOUNT)
               ELSE
                  ! Interpolate by pressure
                  TRACMR = (AER_PLEVS(VERTCOUNT+1)-PCENTER)
                  TRACMR = TRACMR/(AER_PLEVS(VERTCOUNT+1)-
     &                               AER_PLEVS(VERTCOUNT))
                  TRACMR = TRACMR * 
     &               (AER_MR(J,VERTCOUNT+1)-AER_MR(J,VERTCOUNT))
                  TRACMR = AER_MR(J,VERTCOUNT+1) - TRACMR
               ENDIF
               ! Send to STT (currently kg/box)
               IF (SET_CELL) THEN
                  ! Convert from [ppbv] to [kg]
                  TRACMASS = TRACMR * AD(I,J,L) / TCVV(TRAC_IDX)
                  IF (TRAC_ADD) THEN
                     STT(I,J,L,TRAC_IDX) = STT(I,J,L,TRAC_IDX)+TRACMASS
                  ELSE
                     STT(I,J,L,TRAC_IDX) = TRACMASS
                  ENDIF
                ENDIF
            ENDDO
         ENDDO
      ENDDO
 
      ! Return to calling program
      END SUBROUTINE APPLY_2DTRAC
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_basic
!
! !DESCRIPTION: Subroutine EMISS\_BASIC sets surface mixing ratios of N2O,
!  OCS and ozone-depleting substances covered by the Montreal protocol.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_BASIC
!
! !USES:
!
      USE TIME_MOD,     ONLY : GET_MONTH
      USE TRACER_MOD,   ONLY : STT, TCVV
      USE DAO_MOD,      ONLY : AD
      USE LOGICAL_MOD,  ONLY : LBASICEMIS
!
! !INPUT PARAMETERS:
!
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, N
      REAL*8             :: N2O_MULT,OCS_MULT,CFC_MULT,HCFC_MULT
      REAL*8             :: DT_YR

      !=================================================================
      ! EMISS_BASIC begins here!
      !=================================================================

       ! If not correct hour block and day read in new fields
       IF (SFC_MONTH/=GET_MONTH()) THEN
          CALL READ_SFC
          SFC_MONTH    = GET_MONTH()
       END IF
  
      ! Loop over surface grid boxes
      DO J = 1, JJPAR
         DO I = 1, IIPAR
            ! Set kg/box using vol/vol (see CONVERT_UNITS in DAO_MOD)
            DO N=1,N_EMIT
               IF (IDT_EMIT(N) .ne. 0) THEN
                  STT(I,J,1,IDT_EMIT(N))=GRID_EMIT(J,N)*AD(I,J,1)/
     &                                 TCVV(IDT_EMIT(N))
               ENDIF
            ENDDO
         ENDDO ! I
      ENDDO ! J

      IF ( LPRT ) CALL DEBUG_MSG( '### GMISTRAT: Emissions complete' )

      END SUBROUTINE EMISS_BASIC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_sfc
!
! !DESCRIPTION: Subroutine READ\_SFC fills out the surface emissions for some
!  species introduced as part of the unified chemistry upgrade.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_SFC
!
! !USES:
!
      USE TIME_MOD,       ONLY : GET_MONTH,GET_YEAR
      USE TRACERID_MOD,   ONLY : IDTHCFCX
!
! !INPUT PARAMETERS:
!
!
! !OUTPUT VARIABLES:
!
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8,PARAMETER :: OCS_MR = 5.00d-10 ! OCS is fixed to 500 pptv
      INTEGER          :: N, TARG_LINE
      CHARACTER(LEN=20):: LOCAL_NAME

      !=================================================================
      ! READ_SFC begins here!
      !=================================================================

      ! Re-zero GRID_EMIT
      GRID_EMIT = 0d0

      ! Determine which line needs to be read in
      TARG_LINE = (12*(GET_YEAR()-1959)) + GET_MONTH()

      ! Loop over tracer IDs pegged for emissions
      DO N=1,N_EMIT
         LOCAL_NAME = TRIM(TRAC_EMIT(N))
         SELECT CASE (LOCAL_NAME)
            CASE ('CFCX')
               ! Add 3 CFCs
               LOCAL_NAME = TRIM('CFC113')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('CFC114')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('CFC115')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
            CASE ('HCFCX')
               ! Add 2 HCFCs
               LOCAL_NAME = TRIM('HCFC141b')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('HCFC142b')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               ! Also need HCFC123 - no Montreal data
               ! Use surface-level MR from 2D data instead
               TRAC_IDX = IDTHCFCX
               TRAC_SET  = 'C2HF3CL2'
               CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.TRUE.)
            CASE ('H2402')
               ! Limited data available for H2402
               LOCAL_NAME = TRIM('C2BR2F4')
               CALL GET_MONTREAL(LOCAL_NAME,GET_YEAR()-1949,N)
            CASE ('OCS')
               ! Simple boundary condition
               GRID_EMIT(:,N) = OCS_MR
            CASE DEFAULT
               ! Standard 1:1
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
         END SELECT
      ENDDO

      END SUBROUTINE READ_SFC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_montreal
!
! !DESCRIPTION: Subroutine GET\_MONTREAL sets mean surface mixing ratio of a 
!  given species based on forecasts from WHO estimates. Values are
!  interpolated between two monthly averages.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MONTREAL (TRAC_NAME,TARG_LINE,N)
!
! !USES:
!
      USE FILE_MOD,  ONLY : IOERROR
      USE ERROR_MOD, ONLY : ERROR_STOP, ALLOC_ERR
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=20),INTENT(IN)   :: TRAC_NAME
      INTEGER,INTENT(IN)             :: TARG_LINE
      INTEGER,INTENT(IN)             :: N
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: FILENAME
      REAL*8                         :: CURRLINE
      REAL*8                         :: SCALEFAC
      REAL*8                         :: TEMPMR,TEMPYR
      INTEGER                        :: ILINE
      REAL*8                         :: MEAN_MIXRATIO
      INTEGER                        :: IU_FILE, IOS, AS

      !=================================================================
      ! GET_MONTREAL begins here!
      !=================================================================

      FILENAME = TRIM(MONTREAL_FILE_ROOT) // TRIM(TRAC_NAME) // 
     &             '_boundary'

      ! Get a free LUN
      IU_FILE = findFreeLUN()

      IOS=1

      OPEN( IU_FILE, FILE=TRIM( FILENAME ), IOSTAT=IOS )
      IF ( IOS /= 0 ) THEN 
         WRITE(6,*) 'GMISTRAT: FAILED TO FIND ', TRIM(FILENAME)
          CALL IOERROR( IOS, IU_FILE, 'GMISTRAT_MOD:GET_MONTREAL')
      ENDIF

      ! First line not important
      READ(IU_FILE,*,IOSTAT=IOS)
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) 'GMISTRAT: READ ERROR IN ', TRIM(FILENAME)
          CALL IOERROR( IOS, IU_FILE, 'GMISTRAT_MOD:GET_MONTREAL')
      ENDIF

      ! Second line contains scaling factor
      READ(IU_FILE,101,IOSTAT=IOS) SCALEFAC

101   FORMAT(E7.1)

      ! Third and fourth line also to be ignored
      DO ILINE=1,2
         READ(IU_FILE,*,IOSTAT=IOS)
         IF ( IOS /= 0 ) THEN
            WRITE(6,*) 'GMISTRAT: READ ERROR IN ', TRIM(FILENAME)
             CALL IOERROR( IOS, IU_FILE, 'GMISTRAT_MOD:GET_MONTREAL')
         ENDIF
      ENDDO

      ! Skip to line specified
      DO ILINE=1,(TARG_LINE-1)
         READ(IU_FILE,*,IOSTAT=IOS) ! Skip line
         IF ( IOS /= 0 ) THEN 
            WRITE(6,*) 'GMISTRAT: READ ERROR IN ', TRIM(FILENAME)
             CALL IOERROR( IOS, IU_FILE, 'GMISTRAT_MOD:GET_MONTREAL')
         ENDIF
      ENDDO

      ! Read relevant line
      IF (TRAC_NAME.eq.'C2BR2F4') THEN
         WRITE(6,*) 'C2Br2F4 encountered - reading...'
         READ(IU_FILE,103,IOSTAT=IOS) TEMPYR,TEMPMR
      ELSE
         READ(IU_FILE,102,IOSTAT=IOS) TEMPYR,TEMPMR
      ENDIF
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) 'GMISTRAT: READ ERROR IN ', TRIM(FILENAME)
          CALL IOERROR( IOS, IU_FILE, 'GMISTRAT_MOD:GET_MONTREAL')
      ENDIF

      CLOSE(IU_FILE)

      MEAN_MIXRATIO = TEMPMR * SCALEFAC

      GRID_EMIT(:,N) = GRID_EMIT(:,N) + MEAN_MIXRATIO

102   FORMAT(F7.2,F11.2)
      ! C2Br2F4 only
103   FORMAT(F7.2,F8.2)

      ! Return to calling program
      END SUBROUTINE GET_MONTREAL
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: settle_strat_aer
!
! !DESCRIPTION: Subroutine SETTLE\_STRAT\_AER performs gravitational settling
!  of stratospheric aerosols. It is copied largely from GRAV\_SETTLING in
!  sulfate\_mod.F. All of this is ignored if APM is active.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETTLE_STRAT_AER
!
! !USES:
!
      USE DAO_MOD,       ONLY : T, BXHEIGHT, AD
      USE PRESSURE_MOD,  ONLY : GET_PCENTER, GET_PEDGE
      USE TRACER_MOD,    ONLY : STT, XNUMOL, TRACER_MW_G
      USE TRACERID_MOD,  ONLY : IDTNIT,IDTBCPI
      USE TIME_MOD,      ONLY : GET_ELAPSED_SEC, GET_TS_CHEM
      USE COMODE_MOD,    ONLY : JLOP,WERADIUS,ERADIUS
      USE LOGICAL_MOD,   ONLY : LGMISTRAT,LGRAVSTRAT
      USE CHEMGRID_MOD,  ONLY : ITS_IN_THE_STRATMESO
      USE CHEMGRID_MOD,  ONLY : ITS_IN_THE_TROP,ITS_IN_THE_CHEMGRID

      USE CMN_SIZE_MOD        ! Size parameters
      USE CMN_DIAG_MOD        ! ND44
      USE CMN_GCTM_MOD        ! g0
      USE CMN_FJX_MOD,   ONLY : RAA, IND999
      USE ERROR_MOD,     ONLY : IT_IS_NAN,ERROR_STOP     ! Test for NaN
!
! !REVISION HISTORY: 
!  11 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Common
      INTEGER                :: I,      J,     L,        DTCHEM
      REAL*8                 :: DELZ,   DELZ1
      REAL*8                 :: P,      TEMP,  DP,       PDP     
      REAL*8                 :: CONST,  SLIP,  VISC,     FAC1
      REAL*8                 :: FAC2,   FLUX,  AREA_CM2, RHB
      REAL*8                 :: RCM
      REAL*8                 :: TOT1,   TOT2
      INTEGER                :: KLOOP, K
      LOGICAL                :: NATCOL
     
      ! Specific to each class
      REAL*8                 :: RWET(2),CONST_V(2)
      REAL*8                 :: RHO(2),RATIO_R(2),REFF(2)
      REAL*8                 :: VTS(LLPAR,2)
 
      ! Used for old Seinfeld & Pandis slip factor calc
      REAL*8                 :: sp_Lambda, sp_Num

      ! Parameters
#if defined( APM )
      REAL*8, PARAMETER      :: BCDEN = 1800.d0
#else
      REAL*8, PARAMETER      :: BCDEN = 1000.d0 ! density (kg/m3)
#endif

      ! Indexing
      INTEGER, PARAMETER     :: IBC  = 1
      INTEGER, PARAMETER     :: ILIQ = 2
      INTEGER, PARAMETER     :: NSETTLE = 2
      INTEGER                :: IAERO
      LOGICAL                :: RUNCALC

      ! Partitioning
      REAL*8, DIMENSION(3,2) :: PHASEMASS
      REAL*8                 :: SEDMASS
      INTEGER                :: IDTCURRENT

      ! NAT only
      REAL*8                 :: VFALLMAX
      REAL*8,DIMENSION(LLPAR):: VNAT
      INTEGER                :: MAXALT, MINALT
      REAL*8                 :: BXMIN, SEDSTEP, TEMPREAL, P_0, P_ABOVE
      REAL*8                 :: PSED0, PSEDABOVE, INVAIR_0, INVAIR_ABOVE
      REAL*8                 :: XPSC_ABOVE, XPSC_0, XPSC_BELOW
      REAL*8                 :: XNAT_ABOVE, XNAT_0, XNAT_BELOW
      REAL*8                 :: XNO3_ABOVE, XNO3_0, XNO3_BELOW
      REAL*8                 :: XICE_ABOVE, XICE_0, XICE_BELOW
      REAL*8                 :: SEDPSC, SEDNAT, SEDNO3, SEDICE, SEDQUANT
      REAL*8                 :: BELOWGRAD, ABOVEGRAD
      INTEGER                :: LOCALPROFILE, NUMSEDSTEPS, STARTPT, ISED

      !=================================================================
      ! SETTLE_STRAT_AER begins here!
      !=================================================================

      ! Return if gravitational settling disabled
      IF (.not. LGRAVSTRAT) RETURN

      ! Return if it's the start of the run
      IF ( GET_ELAPSED_SEC() == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! First settle liquid aerosols (SLA) using scheme found
      ! elsewhere in GEOS-Chem
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,     L,    VTS,  P,        TEMP, RHB,  RWET ) 
!$OMP+PRIVATE( RATIO_R, RHO,   DP,   PDP,  CONST,    SLIP, VISC )
!$OMP+PRIVATE( DELZ,    DELZ1, TOT1, TOT2, AREA_CM2, FLUX, IAERO      )
!$OMP+PRIVATE( SP_NUM, SP_LAMBDA, RUNCALC )
!$OMP+PRIVATE( CONST_V, K, PHASEMASS, IDTCURRENT, SEDMASS )
!$OMP+PRIVATE( BXMIN, NATCOL, MINALT, MAXALT, VFALLMAX, SEDSTEP )
!$OMP+PRIVATE( LOCALPROFILE, NUMSEDSTEPS, TEMPREAL, STARTPT, ISED )
!$OMP+PRIVATE( PSED0, P_0, P_ABOVE, ABOVEGRAD, PSEDABOVE )
!$OMP+PRIVATE( BELOWGRAD, INVAIR_0, INVAIR_ABOVE )
!$OMP+PRIVATE( SEDNO3, SEDNAT, SEDICE, SEDPSC, SEDQUANT )
!$OMP+PRIVATE( XNO3_ABOVE, XNAT_ABOVE, XICE_ABOVE, XPSC_ABOVE )
!$OMP+PRIVATE( XNO3_0, XNAT_0, XICE_0, XPSC_0, VNAT )
!$OMP+PRIVATE( XNO3_BELOW, XNAT_BELOW, XICE_BELOW, XPSC_BELOW )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR       

         ! Initialize 
         DO L = 1, LLPAR
            VTS(L,:) = 0d0
            VNAT(L) = 0d0
         ENDDO
         NATCOL = .FALSE.

         ! Arbitrary limits
         BXMIN = 1.d20
         MINALT = LLPAR+1
         MAXALT = 0
         VFALLMAX = 0d0

         ! Loop over levels
         DO L = 1, LLPAR

            ! Temperature [K]
            TEMP    = T(I,J,L)
 
            RUNCALC = ITS_IN_THE_STRATMESO(I,J,L)

            IF (RUNCALC) THEN
               ! Need to translate for BC radii
               IF (ITS_IN_THE_CHEMGRID(I,J,L)) THEN
                  KLOOP = JLOP(I,J,L)  
                  RWET(IBC) = WERADIUS(KLOOP,2)*1.d-2
               ELSE
                  ! Use defaults, assume dry (!)
                  RWET(IBC) = RAA(IND999,29)
               ENDIF

               ! Taken from aerosol_mod (MSDENS(2))
               RHO(IBC) = BCDEN

               ! Get aerosol properties
               RWET(ILIQ) = RAD_AER(I,J,L,I_SLA)*1.d-2
               RHO(ILIQ)  = RHO_AER(I,J,L,I_SLA)

               ! Do we need to sediment NAT?
               IF (NDENS_AER(I,J,L,I_SPA).gt.TINY(0d0)) THEN
                  NATCOL = .TRUE.
                  BXMIN = MIN(BXMIN,BXHEIGHT(I,J,L))
                  MINALT = MIN(L,MINALT)
                  MAXALT = MAX(L,MAXALT)
               ENDIF
            ENDIF

            IF (.not.RUNCALC) THEN
               VTS(L,:) = 0d0
               VNAT(L) = 0d0
            ELSE
               ! Pressure at center of the level [kPa]
               P       = GET_PCENTER(I,J,L) * 0.1d0

               ! Calculate common variables first
               sp_Num = P * 1.d3 * 6.023d23 / (8.314 * Temp)
               sp_Lambda = 1.d6/(1.41421*sp_Num*3.141592*(3.7d-10)**2)

               ! Viscosity [Pa*s] of air as a function of temperature 
               VISC = 1.458d-6 * (Temp)**(1.5d0) / ( Temp + 110.4d0 )

               DO IAERO=1,2
                  IF (RWET(IAERO).le.TINY(0d0)) THEN
                     VTS(L,IAERO) = 0d0
                  ELSE
                     ! Dp = particle diameter [um]
                     DP = 2.d0 * RWET(IAERO) * 1.d6        

                     ! PdP = P * dP [hPa * um]
                     PDp     = P * Dp

                     ! Constant
                     CONST = 2.d0 * RHO(IAERO) * RWET(IAERO)**2 * g0 / 9.d0

                     !===========================================================
                     ! NOTE: Slip correction factor calculations following 
                     ! Seinfeld, pp464 which is thought to be more accurate 
                     ! but more computation required. (rjp, 1/24/02)
                     !
                     ! # air molecule number density
                     ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
                     !
                     ! # gas mean free path
                     ! lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
                     !
                     ! # Slip correction
                     ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp     
                     !     &     / (2. * lamda))) / Dp
                     !
                     ! NOTE: Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
                     ! which produces slip correction factore with small error
                     ! compared to the above with less computation.
                     !===========================================================  
                   
                     ! Slip correction factor (as function of P*dp)
                     ! Slip = 1.d0+(15.60d0 + 7.0d0 *
     &               !           EXP(-0.059d0 * PDp)) / PDp
                     ! Reverting to Seinfeld and Pandis
                     Slip = 1. + 2. * sp_Lambda * (1.257 + 0.4*exp(-1.1*
     &                     Dp / (2. * sp_Lambda))) / Dp
   
                     ! Settling velocity [m/s]
                     VTS(L,IAERO) = CONST * Slip / VISC 
                  ENDIF ! RWET
               ENDDO ! IAERO
               
               ! Now solid PSC particles
               IF (NATCOL) THEN
                  ! sp_Num: Air molecule#/m3
                  VNAT(L) = CALC_FALLVEL(RHO_AER(I,J,L,I_SPA),
     &                          RAD_AER(I,J,L,I_SPA),TEMP,P)
                  IF (VNAT(L).gt.VFALLMAX) VFALLMAX = VNAT(L)
               ELSE
                  VNAT(L) = 0.d0
               ENDIF
            ENDIF ! RUNCALC
         ENDDO

         ! First apply simpler SLA sedimentation scheme
         ! Handle model top condition
         L = LLPAR
         DELZ = BXHEIGHT(I,J,L)
         
         DO IAERO=1,2
            CONST_V(IAERO) = 1.d0/(1.d0+DTCHEM*VTS(L,IAERO)/DELZ)
         ENDDO

         ! Zero arrays
         PHASEMASS(:,:) = 0d0

         ! Only want to sediment fraction of species currently
         ! in the aerosol
         DO K = 1,6
            ! Only process transported tracers
            IDTCURRENT = AERFRACIND(K)
            IF (IDTCURRENT.ne.0) THEN
               ! Calculate local phase partitioning
               ! Total upper gridbox mass
               PHASEMASS(3,2) = STT(I,J,L,IDTCURRENT)
               ! Aerosol-phase upper gridbox mass
               PHASEMASS(2,2) = AERFRAC(I,J,L,K)*PHASEMASS(3,2)
               ! Gas-phase upper gridbox mass
               PHASEMASS(1,2) = PHASEMASS(3,2) - PHASEMASS(2,2)

               ! Calculate total sedimented mass
               SEDMASS = PHASEMASS(2,2) * (1.d0-CONST_V(ILIQ))

               ! Remove from upper gridbox
               PHASEMASS(2,2) = PHASEMASS(2,2) - SEDMASS
               PHASEMASS(3,2) = PHASEMASS(1,2) + PHASEMASS(2,2)

               ! Recalculate phase fractions
               IF (PHASEMASS(3,2).gt.TINY(1d0)) THEN
                  AERFRAC(I,J,L,K) = PHASEMASS(2,2)/PHASEMASS(3,2)
               ELSE
                  AERFRAC(I,J,L,K) = 0d0
               ENDIF

               ! Store result
               STT(I,J,L,IDTCURRENT) = PHASEMASS(3,2)
            ENDIF
         ENDDO
         STT(I,J,L,IDTBCPI)= STT(I,J,L,IDTBCPI)* CONST_V(IBC) 

         DO L = LLPAR-1,1,-1
            IF (ITS_IN_THE_TROP(I,J,L+1)) CYCLE
            DELZ  = BXHEIGHT(I,J,L)
            DELZ1 = BXHEIGHT(I,J,L+1)

            DO K=1,6
               IDTCURRENT = AERFRACIND(K)
               IF (IDTCURRENT.ne.0) THEN
                  ! Total upper gridbox mass
                  PHASEMASS(3,2) = STT(I,J,L+1,IDTCURRENT)
                  ! Aerosol-phase upper gridbox mass
                  PHASEMASS(2,2) = AERFRAC(I,J,L+1,K)*PHASEMASS(3,2)
                  ! Gas-phase upper gridbox mass
                  PHASEMASS(1,2) = PHASEMASS(3,2) - PHASEMASS(2,2)
                  
                  ! Total lower gridbox mass
                  PHASEMASS(3,1) = STT(I,J,L,IDTCURRENT)
                  ! Aerosol-phase lower gridbox mass
                  PHASEMASS(2,1) = AERFRAC(I,J,L,K)*PHASEMASS(3,1)
                  ! Gas-phase lower gridbox mass
                  PHASEMASS(1,1) = PHASEMASS(3,1) - PHASEMASS(2,1)
                  
                  ! New lower gridbox mass
                  PHASEMASS(2,1) = 1.d0/(1.d0+DTCHEM*VTS(L,ILIQ)/DELZ)
     &               * (PHASEMASS(2,1)+DTCHEM*VTS(L+1,ILIQ)/DELZ1
     &               * PHASEMASS(2,2))

                  ! Calculate new total mass in lower gridbox
                  PHASEMASS(3,1) = PHASEMASS(2,1) + PHASEMASS(1,1)

                  ! Recalculate phase fraction
                  IF (PHASEMASS(3,1).gt.TINY(1d0)) THEN
                     AERFRAC(I,J,L,K) = PHASEMASS(2,1)/PHASEMASS(3,1)
                  ELSE
                     AERFRAC(I,J,L,K) = 0d0
                  ENDIF

                  ! Store result
                  STT(I,J,L,IDTCURRENT) = PHASEMASS(3,1)
               ENDIF
            ENDDO
            STT(I,J,L,IDTBCPI) = 1.d0/(1.d0+DTCHEM*VTS(L,IBC)/DELZ)
     &              * (STT(I,J,L,IDTBCPI)+DTCHEM*VTS(L+1,IBC)/DELZ1
     &              * STT(I,J,L+1,IDTBCPI)) 
         ENDDO

         ! Now perform trapezoidal scheme for particulates
         ! Calculate maximum allowable timestep (seconds)
         SEDSTEP = BXMIN/VFALLMAX
         IF(DTCHEM.le.SEDSTEP)THEN
            NUMSEDSTEPS = 1
            SEDSTEP = DTCHEM
         ELSE
            ! Will need to run iteratively
            ! Calculate minimum necessary number of steps, limiting to
            ! 10 steps if excessive
            TEMPREAL = DTCHEM/SEDSTEP 
            NUMSEDSTEPS = CEILING(TEMPREAL)
            NUMSEDSTEPS = MIN(10,NUMSEDSTEPS)
            SEDSTEP = DTCHEM/(NUMSEDSTEPS*1.d0)
            VFALLMAX = BXMIN/SEDSTEP
         ENDIF

         SEDSTEPLOOP: DO ISED=1,NUMSEDSTEPS
            STARTPT = MAX(1,MINALT-1)
            ! XPSC is the number of molecules tied up in solid particles
            ! per m3 in a grid box
            L = STARTPT
            XNO3_0 = STT(I,J,L,IDTNIT) * AIRMW / 
     &                      (TRACER_MW_G(IDTNIT)*AD(I,J,L))
            XNAT_0 = XNO3_0 * 4.d0
            XICE_0 = (KG_AER(I,J,L,I_SPA)-STT(I,J,L,IDTNIT))
     &                    *AIRMW/(ICEMW*AD(I,J,L))
            XPSC_0 = XNAT_0 + XICE_0
            P_0 = 100.0d0 * GET_PCENTER(I,J,L)
            INVAIR_0 = AIRMW/AD(I,J,L)
            IF (L .ne. LLPAR) THEN
               INVAIR_ABOVE = AIRMW/AD(I,J,L+1)
               XNO3_ABOVE = STT(I,J,L+1,IDTNIT) * INVAIR_ABOVE / 
     &                             (TRACER_MW_G(IDTNIT))
               ! NAT = HNO3.3H2O = 4 molecules
               XNAT_ABOVE = XNO3_ABOVE * 4.d0
               XICE_ABOVE = (KG_AER(I,J,L+1,I_SPA)-STT(I,J,L+1,IDTNIT))
     &                                    *INVAIR_ABOVE/ICEMW
               XPSC_ABOVE = XNAT_ABOVE + XICE_ABOVE
               P_ABOVE = 100.0d0 * GET_PCENTER(I,J,L+1)
               PSEDABOVE = VNAT(L+1)*SEDSTEP*1.d2*
     &            (GET_PEDGE(I,J,L)-GET_PEDGE(I,J,L+1))/BXHEIGHT(I,J,L)
!               PSEDABOVE=g0 * (1.0d3*AIRMW) * P_ABOVE * 
!     &            VNAT(L+1) / (Rd * T(I,J,L+1))*SEDSTEP
            ELSE
               INVAIR_ABOVE = 0.0d0
               XNO3_ABOVE = 0.0d0
               XNAT_ABOVE = 0.0d0
               XICE_ABOVE = 0.0d0
               XPSC_ABOVE = 0.0d0
               P_ABOVE = GET_PEDGE(I,J,LLPAR+1) * 100.0d0
               PSEDABOVE = 0.0d0
            ENDIF
            ABOVEGRAD = (XPSC_ABOVE-XPSC_0)/(P_0-P_ABOVE)
            SED_LLOOP: DO L=STARTPT,MAXALT
               ! Actually calculating sedimentation for the box above this
               ! one. By the time we get around to processing a given
               ! box, it is actually the (i-1)th box
               ! We are therefore concerned with calculating
               ! sedimentation *out* of the box above us and *into* this
               ! one
               XNO3_BELOW = XNO3_0
               XNAT_BELOW = XNAT_0
               XICE_BELOW = XICE_0
               XPSC_BELOW = XPSC_0
               XNO3_0 = XNO3_ABOVE
               XNAT_0 = XNAT_ABOVE
               XICE_0 = XICE_ABOVE
               XPSC_0 = XPSC_ABOVE
               P_0 = P_ABOVE
               PSED0 = PSEDABOVE
               INVAIR_0 = INVAIR_ABOVE
               
               IF (L.le.(LLPAR-2)) THEN
                  INVAIR_ABOVE = AIRMW/AD(I,J,L+2)
                  XNO3_ABOVE = STT(I,J,L+2,IDTNIT) * INVAIR_ABOVE / 
     &                             TRACER_MW_G(IDTNIT)
                  XNAT_ABOVE = XNO3_ABOVE * 4.d0
                  XICE_ABOVE = (KG_AER(I,J,L+2,I_SPA)-
     &                STT(I,J,L+2,IDTNIT))/ICEMW
                  XPSC_ABOVE = XNAT_ABOVE + XICE_ABOVE
                  P_ABOVE = GET_PCENTER(I,J,L+2) * 100.0d0
                  PSEDABOVE=VNAT(L+2)*SEDSTEP*1.d2*(GET_PEDGE(I,J,L+1)-
     &               GET_PEDGE(I,J,L+2))/BXHEIGHT(I,J,L+2)
!                  PSEDABOVE = g0 * (1.0d3*AIRMW) * P_ABOVE *
!     &               VNAT(L+2)/(Rd * T(I,J,L+2))*SEDSTEP
               ELSE
                  INVAIR_ABOVE = 0.0d0
                  XNO3_ABOVE = 0.0d0
                  XNAT_ABOVE = 0.0d0
                  XICE_ABOVE = 0.0d0
                  XPSC_ABOVE = 0.0d0
                  P_ABOVE = GET_PEDGE(I,J,LLPAR) * 100.0d0
                  PSEDABOVE = 0.0d0
               ENDIF
               ! Note reversal of pressure values as pressure falls with
               ! height. If a value is positive, PSC particle substance
               ! mixing ratios are increasing with height, and vice versa
               BELOWGRAD = ABOVEGRAD
               ABOVEGRAD = (XPSC_ABOVE-XPSC_0)/(P_0-P_ABOVE)

               ! NB: Order of following conditionals is important! Could
               ! theoretically do without localprofile as a variable, but
               ! should be sure to check for the possible (but phenomenally
               ! unlikely) condition that the gradients are both zero
               IF ((ABOVEGRAD .eq. 0d0) .and. (BELOWGRAD .eq. 0d0)) THEN
                  ! Extremelely unlikely case, but use standard scheme
                  ! for flat vertical profile
                  LOCALPROFILE = 0
               ELSEIF ((ABOVEGRAD .le. 0d0) .and.
     &                  (BELOWGRAD .le. 0d0)) THEN
                  ! Consistently decreasing with height - above a peak
                  LOCALPROFILE = -1
               ELSEIF ((ABOVEGRAD .ge. 0d0) .and. 
     &                  (BELOWGRAD .ge. 0d0)) THEN
                  ! Consistently increasing with height - below a peak
                  LOCALPROFILE = +1
               ELSE
                  ! Local minmum or maximum
                  LOCALPROFILE = 0
               ENDIF

               IF (LOCALPROFILE.eq.0) THEN
                  ! Standard sedimentation scheme
                  SEDQUANT = XPSC_0 * PSED0
               ELSEIF (LOCALPROFILE.eq.1) THEN
                  ! Currently below a peak
                  ! Use gradient of bottom and centre boxes
                  IF (XPSC_0 .le. (0.5d0 * XPSC_ABOVE)) THEN
                     ! Interpret as nearing peak; use lower gradient
                     SEDQUANT = (XPSC_0+XPSC_BELOW)*0.5d0*PSED0
                  ELSE
                     SEDQUANT = (XPSC_ABOVE+XPSC_0)*0.50*PSED0
                  ENDIF
               ELSE
                  ! Above a peak
                  SEDQUANT = (XPSC_0+XPSC_BELOW)*0.5d0*PSED0
               ENDIF
               ! Divide sedimenting quantity by the pressure difference
               ! across the box being sedimented from
               ! Note conversion from hPa to Pa for denominator, ie
               ! multiply numerator by (1/100)
               IF (L .ne. LLPAR) THEN
                  SEDQUANT = (1.d-2) * SEDQUANT /
     &               (GET_PEDGE(I,J,L+1)-GET_PEDGE(I,J,L+2))
               ELSE
                  ! This shouldn't be possible?
                  CALL ERROR_STOP('Unknown sedimentation error',
     &                            'gmistrat_mod.F')
               ENDIF
         
               ! Apply limits so that sedimented quantity is:
               ! Greater than or equal to zero
               ! Less than or equal to the total available quantity

               ! Note that we are calculating sedimentation using
               ! the total ice and NAT, but are only actually
               ! transporting the local NO3
               SEDQUANT = MAX(0.0d0,SEDQUANT)
               ! X???_0 refers to the box losing mass (L+1)
               SEDNO3 = SEDQUANT * (XNO3_0/XPSC_0)
               SEDICE = SEDQUANT * (XICE_0/XPSC_0)
               SEDNAT = SEDNO3 * 4.d0
               
               ! Convert v/v to kg/box
               SEDNO3 = SEDNO3*TRACER_MW_G(IDTNIT)/INVAIR_0
               SEDICE = SEDICE*ICEMW/INVAIR_0
               SEDNAT = SEDNAT*NATMW/INVAIR_0
               
               SEDNO3 = MIN(SEDNO3,STT(I,J,L+1,IDTNIT))
               STT(I,J,L,IDTNIT) = STT(I,J,L,IDTNIT) + SEDNO3
               STT(I,J,L+1,IDTNIT)=STT(I,J,L+1,IDTNIT)-SEDNO3

               ! Now correct aerosol totals
               SEDPSC = SEDNAT + SEDICE
               SEDPSC = MIN(SEDPSC,KG_AER(I,J,L+1,I_SPA))
               KG_AER(I,J,L,I_SPA) = KG_AER(I,J,L,I_SPA) + SEDPSC
               KG_AER(I,J,L+1,I_SPA)=KG_AER(I,J,L+1,I_SPA)-SEDPSC
               
            ENDDO SED_LLOOP
         ENDDO SEDSTEPLOOP
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE SETTLE_STRAT_AER
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_h2so4_gas
!
! !DESCRIPTION: Subroutine CALC\_H2SO4\_GAS calculates the fraction of strat.
!  SO4 aerosol which can be considered to be gaseous H2SO4.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_H2SO4_GAS
!
! !USES:
!
      USE DAO_MOD,      ONLY : T, AD
      USE TRACER_MOD,   ONLY : STT, TRACER_MW_G
      USE TRACERID_MOD, ONLY : IDTSO4
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE CHEMGRID_MOD, ONLY : ITS_IN_THE_TROP

      USE CMN_GCTM_MOD, ONLY : AIRMW
!
! !REVISION HISTORY: 
!  11 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST=.TRUE.
      REAL*8,PARAMETER   :: GF_THRESHOLD = 0.0d0
      REAL*8,PARAMETER   :: GF_RANGE     = 1.0d-8
      REAL*8,PARAMETER   :: GF_DELTAHBYR = 10156.d0
      REAL*8,PARAMETER   :: GF_T0        = 360.d0
      REAL*8,PARAMETER   :: GF_TC        = 905.d0
      REAL*8,SAVE        :: GF_LOGP0
      REAL*8,SAVE        :: GF_BFACTOR
      REAL*8,SAVE        :: GF_ATMCONV
      REAL*8,SAVE        :: GF_INVT0
      REAL*8             :: GF_INVT,GF_LOGPSULFATE,GF_CFACTOR
      REAL*8             :: GF_AFACTOR
      REAL*8             :: GF_PP,GF_PVAP,GF_DIFF

      INTEGER            :: I, J, L
      REAL*8             :: PCENTER, TCENTER, H2SO4SUM, INVAIR

      !=================================================================
      ! CALC_H2SO4_GAS begins here!
      !=================================================================

      IF (FIRST) THEN
         FIRST = .FALSE.
         ! Calculate H2SO4 gas phase prefactors
         GF_INVT0 = 1.d0/GF_T0
         GF_LOGP0 = (-1.d0*GF_DELTAHBYR*GF_INVT0) + 16.259d0
         GF_BFACTOR = 0.38d0/(GF_TC - GF_T0)
         GF_ATMCONV = LOG(101325.d0)
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PCENTER, TCENTER, INVAIR ) 
!$OMP+PRIVATE( H2SO4SUM, GF_DIFF )
!$OMP+PRIVATE( GF_PP, GF_INVT, GF_CFACTOR, GF_LOGPSULFATE, GF_PVAP )
!$OMP+SCHEDULE( DYNAMIC )

      DO J=1, JJPAR
      DO I=1, IIPAR
      LLOOP: DO L=1, LLPAR
         ! Only interested in low-pressure boxes
         PCENTER = GET_PCENTER(I,J,L)
         TCENTER = T(I,J,L)
         INVAIR  = AIRMW / AD(I,J,L)
         IF (PCENTER.ge.1.d2) THEN
            AERFRAC(I,J,L,1) = 1d0
            CYCLE LLOOP
         ELSEIF (ITS_IN_THE_TROP(I,J,L)) THEN
            ! Don't want to interfere with tropospheric aerosols
            AERFRAC(I,J,L,1) = 1d0
            CYCLE LLOOP
         ELSE
            H2SO4SUM = STT(I,J,L,IDTSO4  )*INVAIR/TRACER_MW_G(IDTSO4  )
            ! Use approximation from Kumala (1990)
            GF_PP = H2SO4SUM*PCENTER
            GF_INVT = 1./TCENTER
            GF_CFACTOR = 1.d0+(LOG(GF_T0*GF_INVT))-(GF_T0*GF_INVT)
            GF_LOGPSULFATE = GF_LOGP0 + (GF_DELTAHBYR*(GF_INVT0 -
     &               GF_INVT + (GF_BFACTOR*GF_CFACTOR)))
            GF_LOGPSULFATE = GF_LOGPSULFATE + GF_ATMCONV
            GF_PVAP = 1.d-2 * EXP(GF_LOGPSULFATE)
            GF_DIFF = (GF_PVAP+GF_THRESHOLD) - GF_PP
            IF (GF_DIFF .lt. 0) THEN
               AERFRAC(I,J,L,1) = 1.d0
            ELSEIF (GF_DIFF .lt. GF_RANGE) THEN
               AERFRAC(I,J,L,1) = 1.d0-(GF_DIFF/GF_RANGE)
            ELSE
               AERFRAC(I,J,L,1) = 0d0
            ENDIF
         ENDIF
      ENDDO LLOOP
      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      IF ( LPRT ) CALL DEBUG_MSG( '### GMISTRAT: H2SO4 partitioned' )

      END SUBROUTINE CALC_H2SO4_GAS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SO4_PHOTFRAC
!
! !DESCRIPTION: FUNCTION SO4\_PHOTFRAC returns the fraction of H2SO4 which 
!  is available for photolysis.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION SO4_PHOTFRAC(I,J,L)
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: I,J,L      ! Location indices
!
! !OUTPUT VARIABLES:
!
!      REAL*8, INTENT(OUT)         :: PHOTFRAC   ! Gaseous fraction of H2SO4
!
! !REVISION HISTORY: 
!  11 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! SO4_PHOTFRAC begins here!
      !=================================================================

      SO4_PHOTFRAC = 1.d0 - AERFRAC(I,J,L,1)

      END FUNCTION SO4_PHOTFRAC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_strat_aer
!
! !DESCRIPTION: Subroutine CALC\_STRAT\_AER calculates aerosol properties    
!  stratosphere using the thermodynamic parameterization described in
!  Kirner et al. (`Simulation of polar stratospheric clouds in the
!  chemistry-climate-model EMAC via the submodel PSC', Geosci. Mod. Dev.,
!  4, 169-182, 2011).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_STRAT_AER
!
! !USES:
!
      USE DAO_MOD,          ONLY : T, AD, CLDF, QI, AIRVOL, AVGW
      USE PRESSURE_MOD,     ONLY : GET_PCENTER
      USE LOGICAL_MOD,      ONLY : LPRT, LHOMNUCNAT, LSOLIDPSC
      USE ERROR_MOD,        ONLY : DEBUG_MSG
      USE GRID_MOD,         ONLY : GET_YEDGE
      USE TRACER_MOD,       ONLY : STT, TRACER_MW_G
      USE TRACERID_MOD,     ONLY : IDTSO4
      USE TRACERID_MOD,     ONLY : IDTHNO3, IDTNIT, IDTBrNO3, IDTHBr
      USE TRACERID_MOD,     ONLY : IDTHOBr, IDTHCl, IDTClNO3, IDTHOCl
      USE CHEMGRID_MOD,     ONLY : ITS_IN_THE_STRATMESO
      USE CHEMGRID_MOD,     ONLY : ITS_IN_THE_CHEMGRID
      USE COMODE_MOD,       ONLY : CSPEC,JLOP
      USE ERROR_MOD,        ONLY : ERROR_STOP

      USE CMN_GCTM_MOD,     ONLY : AIRMW, PI
!
! !REVISION HISTORY: 
!  13 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Limits on PSC formation
      REAL*8,PARAMETER      :: PSC_MAXLAT =  45.0d0
      REAL*8,PARAMETER      :: PSC_MINLAT = -55.0d0
      REAL*8,PARAMETER      :: PSC_PMAX   =  18.0d3
      REAL*8,PARAMETER      :: PSC_PMIN   =   5.0d2

      ! Saturation and nucleation properties of HNO3
      REAL*8, PARAMETER     :: TSATHNO3_A = -2.7836d0
      REAL*8, PARAMETER     :: TSATHNO3_B = -0.00088d0
      REAL*8, PARAMETER     :: TSATHNO3_C = 38.9855d0
      REAL*8, PARAMETER     :: TSATHNO3_D = -11397.0d0
      REAL*8, PARAMETER     :: TSATHNO3_E = 0.009179d0

      ! Saturation and nucleation properties of water
      REAL*8, PARAMETER     :: TSATH2O_A  = -2663.5d0
      REAL*8, PARAMETER     :: TSATH2O_B  = 12.537d0

      ! Peak pressure at which NAT can form homogeneously
      REAL*8, PARAMETER     :: P_MAXNAT = 1.40d4 ! Pa

      ! Universal gas constant
      REAL*8, PARAMETER     :: R_UNIV = 8.31448d0

      ! Maximum temperature for PSC formation (K)
      REAL*8, PARAMETER     :: T_MAX = 215.0d0

      ! Limits on NAT/ice formation
      REAL*8, PARAMETER     :: MIN_RAD = 1.0d-7 ! m
      REAL*8, PARAMETER     :: MAX_NDENS=42.0d3 ! #/m3

      ! Local conditions
      REAL*8                :: TCENTER, PCENTER, PCENTER_PA, DENAIR
      REAL*8                :: INVAIR, TINV, TOFFSET

      ! Gridbox mixing ratios and partial pressures
      REAL*8                :: HNO3SUM, H2OSUM
      REAL*8                :: HNO3PP,  H2OPP
      REAL*8                :: PSATHNO3_SUPERCOOL
      REAL*8                :: PSATH2O_SUPERSAT
      REAL*8                :: PSATHNO3, PSATH2O
      REAL*8                :: H2SO4SUM
      REAL*8                :: ClNO3SUM, HClSUM, HOClSUM
      REAL*8                :: BrNO3SUM, HBrSUM, HOBrSUM

      ! Gridbox aerosol and phase data
      REAL*8                :: HNO3_BOX_G, HNO3_BOX_L, HNO3_BOX_S
      REAL*8                :: H2O_BOX_G,  H2O_BOX_L,  H2O_BOX_S
      REAL*8                :: H2SO4_BOX_G,H2SO4_BOX_L
      REAL*8                :: HCl_BOX_G,  HCl_BOX_L
      REAL*8                :: HOCl_BOX_G, HOCl_BOX_L
      REAL*8                :: HBr_BOX_G,  HBr_BOX_L
      REAL*8                :: HOBr_BOX_G, HOBr_BOX_L
      REAL*8                :: HNO3GASFRAC, HClGASFRAC, HOClGASFRAC
      REAL*8                :: HBrGASFRAC, HOBrGASFRAC
      REAL*8                :: VOL_NAT, VOL_ICE, VOL_SLA
      REAL*8                :: RAD_AER_BOX,RHO_AER_BOX
      REAL*8                :: KG_AER_BOX,NDENS_AER_BOX,SAD_AER_BOX
      REAL*8                :: KG_NAT, KG_ICE, KG_NO3

      ! SLA weight fractions
      REAL*8                :: W_H2SO4, W_H2O, W_HNO3
      REAL*8                :: W_HCl, W_HOCl, W_HBr, W_HOBr

      ! Reaction prefactors
      REAL*8                :: KHET_COMMON
      REAL*8                :: KHET_SPECIFIC

      ! Grid box location
      REAL*8                :: BOX_LAT_S, BOX_LAT_N
      LOGICAL               :: IS_VALID, IS_POLAR, IS_STRAT

      ! Local properties
      REAL*8, DIMENSION(11) :: GAMMA_BOX
      INTEGER               :: STATE_LOCAL

      ! Loop variables
      INTEGER               :: I, J, L, K, JLOOP

      !=================================================================
      ! CALC_STRAT_AER begins here!
      !=================================================================

      IF (LPRT) THEN
         CALL DEBUG_MSG( '### GMISTRAT: start CALC_STRAT_AER' )
      ENDIF

      ! Partition H2SO4 before proceeding
      CALL CALC_H2SO4_GAS

      ! Loop over latitude boxes first
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,J,L,K,IS_POLAR,IS_STRAT,PCENTER,PCENTER_PA )
!$OMP+PRIVATE( DENAIR, INVAIR, PSATHNO3, PSATH2O, JLOOP )
!$OMP+PRIVATE( VOL_NAT, VOL_ICE, VOL_SLA, PSATHNO3_SUPERCOOL )
!$OMP+PRIVATE( TCENTER, TINV, IS_VALID, BOX_LAT_S, BOX_LAT_N )
!$OMP+PRIVATE( RAD_AER_BOX, RHO_AER_BOX, KG_AER_BOX, NDENS_AER_BOX )
!$OMP+PRIVATE( SAD_AER_BOX, KG_NAT, KG_ICE, KG_NO3,PSATH2O_SUPERSAT )
!$OMP+PRIVATE( H2OSUM,H2OPP,H2O_BOX_S,H2O_BOX_L,H2O_BOX_G, H2SO4SUM )
!$OMP+PRIVATE( HNO3SUM,HNO3PP,HNO3_BOX_S,HNO3_BOX_L,HNO3_BOX_G )
!$OMP+PRIVATE( BrNO3SUM,HBrSUM,HOBrSUM,ClNO3SUM,HClSUM,HOClSUM )
!$OMP+PRIVATE( STATE_LOCAL,GAMMA_BOX,HBrGASFRAC,HOBrGASFRAC )
!$OMP+PRIVATE( HNO3GASFRAC,HClGASFRAC,HOClGASFRAC,TOFFSET )
!$OMP+PRIVATE( W_H2SO4,W_H2O,W_HCl,W_HOCl,W_HBr,W_HOBr,W_HNO3 )
!$OMP+PRIVATE( HCl_BOX_G,HCl_BOX_L,HOCl_BOX_G,HOCl_BOX_L,H2SO4_BOX_G )
!$OMP+PRIVATE( HBr_BOX_G,HBr_BOX_L,HOBr_BOX_G,HOBr_BOX_L,H2SO4_BOX_L )
!$OMP+PRIVATE( KHET_COMMON, KHET_SPECIFIC )
!$OMP+SCHEDULE( DYNAMIC )
      DO J=1,JJPAR
         BOX_LAT_S = GET_YEDGE(1,J,1)
         BOX_LAT_N = GET_YEDGE(1,J+1,1)
         IS_POLAR = ((BOX_LAT_S.lt.PSC_MINLAT).or.
     &               (BOX_LAT_N.gt.PSC_MAXLAT))
         ! Loop over longitudinal boxes
         DO I=1,IIPAR
            LLOOP: DO L=1,LLPAR

               ! Get local conditions
               PCENTER = GET_PCENTER(I,J,L) 
               PCENTER_PA = PCENTER*1.d2
               TCENTER = T(I,J,L)
               TINV = 1.d0/TCENTER

               ! Apply other limits from Kirner et al.
               IS_STRAT = ITS_IN_THE_STRATMESO(I,J,L)
               IS_VALID = (IS_POLAR).and.(IS_STRAT).and.       
     &                        (.not.((PCENTER_PA.lt.PSC_PMIN).or.
     &                               (PCENTER_PA.gt.PSC_PMAX)))

               ! Initialize variables
               VOL_NAT      = 0.0d0
               KG_NAT       = 0.0d0
               VOL_ICE      = 0.0d0
               KG_ICE       = 0.0d0
 
               RAD_AER_BOX  = 0.0d0
               RHO_AER_BOX  = 0.0d0
               KG_AER_BOX   = 0.0d0
               NDENS_AER_BOX= 0.0d0

               H2O_BOX_S    = 0.0d0
               H2O_BOX_L    = 0.0d0
               H2O_BOX_G    = 0.0d0
               HNO3_BOX_S   = 0.0d0
               HNO3_BOX_L   = 0.0d0
               HNO3_BOX_G   = 0.0d0
               H2SO4_BOX_L  = 0.0d0
               H2SO4_BOX_G  = 0.0d0
               HCl_BOX_L    = 0.0d0
               HCl_BOX_G    = 0.0d0
               HOCl_BOX_L   = 0.0d0
               HOCl_BOX_G   = 0.0d0
               HBr_BOX_L    = 0.0d0
               HBr_BOX_G    = 0.0d0
               HOBr_BOX_L   = 0.0d0
               HOBr_BOX_G   = 0.0d0

               STATE_LOCAL  = STATE_PSC(I,J,L)

               ! Calculate local air density
               DENAIR = AVOGADRO * PCENTER_PA / (TCENTER*R_UNIV)
               INVAIR = AIRMW / AD(I,J,L)

               ! Get available NO3 mass
               KG_NO3 = (STT(I,J,L,IDTHNO3)*
     &                   TRACER_MW_G(IDTNIT)/TRACER_MW_G(IDTHNO3)) + 
     &                   STT(I,J,L,IDTNIT)

               ! Calculate HNO3 mixing ratio
               HNO3SUM = STT(I,J,L,IDTHNO3)*INVAIR/
     &                    TRACER_MW_G(IDTHNO3)
               HNO3SUM = HNO3SUM + STT(I,J,L,IDTNIT)*INVAIR/
     &                    TRACER_MW_G(IDTNIT)

               ! Revert to static H2O for now
               H2OSUM = AVGW(I,J,L)

               ! Calculate partial pressures (Pa)
               HNO3PP = PCENTER_PA * HNO3SUM
               H2OPP  = PCENTER_PA * H2OSUM

               IF (.not.IS_VALID) THEN
                  ! No PSCs (SSA only)
                  STATE_LOCAL = 0
               ELSEIF (TCENTER.gt.T_MAX) THEN
                  ! STS only
                  STATE_LOCAL = 1
               ELSE
                  ! Calculate saturation pressures
                  PSATHNO3 = 10.d0**(((TSATHNO3_A+(TCENTER*
     &            TSATHNO3_B))*LOG10(H2OPP*760.0d0/101325.0d0))+
     &            (TSATHNO3_C+(TSATHNO3_D*TINV)+(TSATHNO3_E*
     &            TCENTER)))*101325.0d0/760.0d0

                  PSATH2O = 10.0d0**((TSATH2O_A*TINV)+TSATH2O_B)

                  ! Supersaturation requirement for ice
                  PSATH2O_SUPERSAT = PSATH2O * P_ICE_SUPERSAT

                  ! If homogeneous NAT nucleation allowed, calculate
                  ! threshold saturation pressure
                  IF (LHOMNUCNAT) THEN
                     ! Calculate as if temperature is (T+T_NAT_SUPERCOOL)
                     TOFFSET = TCENTER + T_NAT_SUPERCOOL
                     PSATHNO3_SUPERCOOL = 10.d0**(((TSATHNO3_A+(TOFFSET*
     &                  TSATHNO3_B))*LOG10(H2OPP*760.0d0/101325.0d0))+
     &                  (TSATHNO3_C+(TSATHNO3_D/TOFFSET)+(TSATHNO3_E*
     &                  TOFFSET)))*101325.0d0/760.0d0
                  ELSE
                     ! Make homogeneous nucleation impossible
                     TOFFSET = 280.0d0
                     PSATHNO3_SUPERCOOL = HNO3PP + 1d0
                  ENDIF
                  
                  ! Use local ice mass ratio from GEOS-5 data
                  KG_ICE = QI(I,J,L)*CLDF(L,I,J)*AD(I,J,L)
                  H2O_BOX_S = QI(I,J,L)*CLDF(L,I,J)*AIRMW/ICEMW
                  VOL_ICE = H2O_BOX_S * DENAIR * (1.d-3) * ICEMW /
     &                             (DENSICE*AVOGADRO) ! m3 ice/m3 air

                  ! Ice exists/possible?
                  IF ((H2OPP.gt.PSATH2O_SUPERSAT).or.
     &                   (H2O_BOX_S.gt.TINY(0d0))) THEN
                     STATE_LOCAL = 3
                     HNO3_BOX_S = 1d0
                  ELSE
                     ! If ice not possible could still have NAT
                     IF (STATE_LOCAL.gt.1) THEN
                        ! Had PSCs here last time
                        IF (((HNO3PP.gt.PSATHNO3_SUPERCOOL).and.
     &                       (PSATHNO3.gt.PSATH2O)).or.
     &                       (PCENTER.lt.P_MAXNAT)) THEN
                           STATE_LOCAL = 2
                           HNO3_BOX_S = 1d0
                        ELSE
                           ! STS only
                           STATE_LOCAL = 1
                        ENDIF 
                     ELSE
                        ! STS only
                        STATE_LOCAL = 1
                     ENDIF
                  ENDIF
               ENDIF

               ! Store state
               STATE_PSC(I,J,L) = STATE_LOCAL

               ! Only continue if we want online solid PSCs
               IF (LSOLIDPSC) THEN

                  ! Calculate NAT if relevant
                  IF ((HNO3_BOX_S.gt.TINY(1d0)).and.
     &                               (STATE_LOCAL.ge.2)) THEN
                     HNO3_BOX_S = (HNO3PP-PSATHNO3)/PCENTER_PA
                     HNO3_BOX_S = MAX(0d0,HNO3_BOX_S)

                     ! Calculate m3 NAT/m3 air
                     ! HNO3_BOX_S is the number of moles of HNO3
                     ! which will be frozen into HNO3.3H2O (NAT)
                     ! Therefore volume calculation must be done
                     ! with care!
                     VOL_NAT = HNO3_BOX_S * DENAIR * (1.d-3) *
     &                  NATMW / (DENSNAT*AVOGADRO)
                     KG_NAT = HNO3_BOX_S * NATMW *
     &                        AD(I,J,L) / AIRMW
                  ENDIF

                  ! Calculate particle properties
                  IF (STATE_LOCAL.lt.2) THEN
                     ! Zero all!
                     KG_AER_BOX   = 0d0
                     RAD_AER_BOX  = 0d0
                     RHO_AER_BOX  = DENSICE
                     NDENS_AER_BOX= 0d0
                  ELSE
                     KG_AER_BOX = KG_NAT + KG_ICE
                     RAD_AER_BOX = MIN_RAD
                     NDENS_AER_BOX = (3.0d0*(VOL_NAT+VOL_ICE)/
     &                  (4.0d0*PI*(RAD_AER_BOX**3.0d0)))
                     IF (NDENS_AER_BOX.gt.MAX_NDENS) THEN
                        NDENS_AER_BOX = MAX_NDENS
                        RAD_AER_BOX = (3.0d0*(VOL_NAT+VOL_ICE)/
     &                     (4.0d0*PI*MAX_NDENS))**(1.d0/3.d0)
                     ENDIF
                     RHO_AER_BOX = 
     &                  ((VOL_ICE*DENSICE)+(VOL_NAT*DENSNAT))/
     &                                    (VOL_ICE+VOL_NAT)
                  ENDIF

                  ! Calculate SAD (cm2/cm3)
                  SAD_AER_BOX = 4.0d-2 * RAD_AER_BOX *
     &               RAD_AER_BOX * NDENS_AER_BOX * PI
               ELSE
                  ! Solid PSCs not active
                  RAD_AER_BOX = 0d0
                  RHO_AER_BOX = 1000d0
                  KG_AER_BOX = 0d0
                  NDENS_AER_BOX = 0d0
                  SAD_AER_BOX = 0d0
                  HNO3_BOX_S = 0d0
                  H2O_BOX_S = 0d0
               ENDIF

               ! Store in outer arrays
               RAD_AER(I,J,L,I_SPA)  = RAD_AER_BOX*1.d2 ! cm
               RHO_AER(I,J,L,I_SPA)  = RHO_AER_BOX ! kg/m3
               KG_AER(I,J,L,I_SPA)   = KG_AER_BOX  ! kg
               NDENS_AER(I,J,L,I_SPA)= NDENS_AER_BOX!#/m3
               SAD_AER(I,J,L,I_SPA)  = SAD_AER_BOX ! cm2/cm3

               ! Repartition NIT and HNO3 in strat/meso
               IF (LSOLIDPSC.and.IS_STRAT) THEN
                  ! Convert NAT from kg NAT to kg NO3
                  STT(I,J,L,IDTNIT) = KG_NAT*TRACER_MW_G(IDTNIT)/NATMW
                  ! Remove (kg NO3 as NAT) from total kg NO3
                  ! then convert to kg HNO3
                  STT(I,J,L,IDTHNO3) = (KG_NO3-STT(I,J,L,IDTNIT))
     &                  *TRACER_MW_G(IDTHNO3)/TRACER_MW_G(IDTNIT) 
               ENDIF

               ! Now start liquid aerosol consideration
               ! Start by assuming all non-solid H2O/HNO3 is gaseous
               HNO3_BOX_G = HNO3SUM - HNO3_BOX_S
               HNO3_BOX_L = 0d0
               H2O_BOX_G = H2OSUM - H2O_BOX_S
               H2O_BOX_L = 0d0

               ! Calculate mixing ratios of other relevant species
               H2SO4SUM = STT(I,J,L,IDTSO4  )*INVAIR/
     &                  TRACER_MW_G(IDTSO4  )
               BrNO3SUM = STT(I,J,L,IDTBrNO3)*INVAIR/
     &                  TRACER_MW_G(IDTBrNO3)
               ClNO3SUM = STT(I,J,L,IDTClNO3)*INVAIR/
     &                  TRACER_MW_G(IDTClNO3)
               HOClSUM  = STT(I,J,L,IDTHOCl )*INVAIR/
     &                  TRACER_MW_G(IDTHOCl )
               HClSUM   = STT(I,J,L,IDTHCl  )*INVAIR/
     &                  TRACER_MW_G(IDTHCl  )
               HOBrSUM  = STT(I,J,L,IDTHOBr )*INVAIR/
     &                  TRACER_MW_G(IDTHOBr )
               HBrSUM   = STT(I,J,L,IDTHBr  )*INVAIR/
     &                  TRACER_MW_G(IDTHBr  )

               ! H2SO4 gas fraction calculated earlier throughout grid
               ! Consider gaseoues H2SO4 to be unavailable for SLA
               H2SO4_BOX_L = H2SO4SUM * AERFRAC(I,J,L,1)
               H2SO4_BOX_G = H2SO4SUM - H2SO4_BOX_L

               ! Zero local properties
               RHO_AER_BOX   = 1000d0
               RAD_AER_BOX   = 0d0 
               KG_AER_BOX    = 0d0
               NDENS_AER_BOX = 0d0
               SAD_AER_BOX   = 0d0
               VOL_SLA = 0d0

               IF (.not.IS_STRAT) THEN
                  ! Use JPL 10-06/Oslo CTM data, where available,
                  ! for conventional sulfates/H2SO4
                  GAMMA_BOX(1)  = 0.1d0
                  GAMMA_BOX(2)  = 0.0d0
                  GAMMA_BOX(3)  = 0.0d0
                  GAMMA_BOX(4)  = 0.0d0
                  GAMMA_BOX(5)  = 0.3d0
                  GAMMA_BOX(6)  = 0.4d0
                  GAMMA_BOX(7)  = 0.9d0
                  GAMMA_BOX(8)  = 0.0d0
                  GAMMA_BOX(9)  = 0.0d0
                  GAMMA_BOX(10) = 0.2d0
                  GAMMA_BOX(11) = 0.0d0
               ELSEIF (H2SO4_BOX_L.lt.1d-15) THEN
                  ! No aerosol to speak of
                  DO K=1,11
                     GAMMA_BOX(K) = 0.0d0
                  ENDDO
               ELSE
                  IF (STATE_LOCAL.eq.0) THEN
                     ! Allow binary H2SO4.nH2O only
                     CALL TERNARY (PCENTER,TCENTER,H2OSUM,H2SO4_BOX_L,
     &                  0.d0   ,HClSUM,HOClSUM,HBrSUM,HOBrSUM,
     &                  W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr,
     &                  HNO3GASFRAC,HClGASFRAC,HOClGASFRAC,
     &                  HBrGASFRAC,HOBrGASFRAC,VOL_SLA,RHO_AER_BOX)

                        ! For safety's sake, zero out HNO3 uptake
                        HNO3GASFRAC = 1.d0
                        W_H2O = W_H2O + W_HNO3
                        W_HNO3 = 0.d0
                        HNO3_BOX_G = HNO3SUM - HNO3_BOX_S
                        HNO3_BOX_L = 0.d0
                  ELSE
                     ! As per Buchholz, use only non-NAT HNO3 for STS
                     HNO3_BOX_G = HNO3SUM - HNO3_BOX_S
                     CALL TERNARY (PCENTER,TCENTER,H2OSUM,H2SO4_BOX_L,
     &                  HNO3_BOX_G,HClSUM,HOClSUM,HBrSUM,HOBrSUM,
     &                  W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr,
     &                  HNO3GASFRAC,HClGASFRAC,HOClGASFRAC,
     &                  HBrGASFRAC,HOBrGASFRAC,VOL_SLA,RHO_AER_BOX)
        
                        ! Partition HNO3 here for safety
                        HNO3_BOX_G = HNO3_BOX_G*HNO3GASFRAC
                        HNO3_BOX_L = HNO3SUM - (HNO3_BOX_G+HNO3_BOX_S)
                  ENDIF

                  ! Partition minor species
                  HCl_BOX_G  = HClSUM *HClGASFRAC
                  HCl_BOX_L  = HClSUM -HCl_BOX_G
                  HOCl_BOX_G = HOClSUM*HOClGASFRAC
                  HOCl_BOX_L = HOClSUM-HOCl_BOX_G
                  HBr_BOX_G  = HBrSUM *HBrGASFRAC
                  HBr_BOX_L  = HBrSUM -HBr_BOX_G
                  HOBr_BOX_G = HOBrSUM*HOBrGASFRAC
                  HOBr_BOX_L = HOBrSUM-HOBr_BOX_G

                  ! Calculate SLA parameters (Grainger 1995)
                  SAD_AER_BOX = SLA_VA*(VOL_SLA**0.751d0)        ! cm2/cm3
                  RAD_AER_BOX = SLA_VR*SLA_RR*(VOL_SLA**0.249d0) ! m
                  KG_AER_BOX  = RHO_AER_BOX*VOL_SLA*AIRVOL(I,J,L)! kg

                  IF (VOL_SLA.gt.1.d-30) THEN
                     ! Approximate particles as spherical for calculation
                     ! of aerosol number density
                     NDENS_AER_BOX = VOL_SLA*3.d0/
     &                       (4.d0*PI*(RAD_AER_BOX**3.d0))

                     ! DENAIR in #/m3 - convert to #/cm3
                     ! RHO_AER_BOX in kg/m3 - convert to g/cm3
                     ! RAD_AER_BOX in m - convert to cm
                     CALL  CALC_SLA_GAMMA(DENAIR*1.d-6,TCENTER,PCENTER,
     &                   W_H2SO4,H2OSUM,HClSUM,HBrSUM,HOBrSUM,
     &                   ClNO3SUM,BrNO3SUM,RHO_AER_BOX*1.d-3,
     &                   RAD_AER_BOX*1.d2,GAMMA_BOX)
                  ELSE
                     ! Ignore SLA
                     DO K=1,11
                        GAMMA_BOX(K) = 0.0d0
                     ENDDO
                  ENDIF
               ENDIF

               ! Send properties to larger array
               ! Convert sticking coefficients into
               ! premultiplying factors (Kirner)
               KHET_COMMON = 0.25d0*MOLEC_SPEED(TCENTER,1d0)

               ! N2O5 + H2O/HCl
               KHET_SPECIFIC= KHET_COMMON*ISR_N2O5
               KHETI_SLA(I,J,L,1)  = GAMMA_BOX(1 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,2)  = GAMMA_BOX(2 )*KHET_SPECIFIC

               ! ClNO3 + H2O/HCl/HBr
               KHET_SPECIFIC= KHET_COMMON*ISR_ClNO3
               KHETI_SLA(I,J,L,3)  = GAMMA_BOX(3 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,4)  = GAMMA_BOX(4 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,5)  = GAMMA_BOX(5 )*KHET_SPECIFIC
 
               ! BrNO3 + H2O/HCl
               KHET_SPECIFIC= KHET_COMMON*ISR_BrNO3
               KHETI_SLA(I,J,L,6)  = GAMMA_BOX(6 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,7)  = GAMMA_BOX(7 )*KHET_SPECIFIC

               ! HOCl + HCl/HBr
               KHET_SPECIFIC= KHET_COMMON*ISR_HOCl
               KHETI_SLA(I,J,L,8)  = GAMMA_BOX(8 )*KHET_SPECIFIC
               KHETI_SLA(I,J,L,9)  = GAMMA_BOX(9 )*KHET_SPECIFIC

               ! HOBr + HBr/HCl
               KHET_SPECIFIC= KHET_COMMON*ISR_HOBr
               KHETI_SLA(I,J,L,10) = GAMMA_BOX(10)*KHET_SPECIFIC
               KHETI_SLA(I,J,L,11) = GAMMA_BOX(11)*KHET_SPECIFIC

               RAD_AER(I,J,L,I_SLA)  = RAD_AER_BOX*1.d2 ! cm
               RHO_AER(I,J,L,I_SLA)  = RHO_AER_BOX      ! kg/m3
               KG_AER(I,J,L,I_SLA)   = KG_AER_BOX       ! kg
               NDENS_AER(I,J,L,I_SLA)= NDENS_AER_BOX    ! #/m3
               SAD_AER(I,J,L,I_SLA)  = SAD_AER_BOX      ! cm2/cm3
            ENDDO LLOOP ! L
         ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO

      IF (LPRT) THEN
         CALL DEBUG_MSG( '### GMISTRAT: after CALC_STRAT_AER' )
         WRITE(6,'(a)') 'Strat. aerosol diagnostics summary:'
         WRITE(6,'(a,1x,E12.4)') 'Peak particle # density (#/m3): ',
     &                               MAXVAL(NDENS_AER(:,:,:,I_SPA))
         WRITE(6,'(a,1x,E12.4)') 'Peak particle density (kg/m3) : ',
     &                               MAXVAL(RHO_AER(:,:,:,I_SPA))
         WRITE(6,'(a,1x,E12.4)') 'Peak particle SAD (cm2/cm3)   : ',
     &                               MAXVAL(SAD_AER(:,:,:,I_SPA))
         WRITE(6,'(a,1x,E12.4)') 'Peak particle radius (cm)     : ',
     &                               MAXVAL(RAD_AER(:,:,:,I_SPA))
         WRITE(6,'(a,1x,E12.4)') 'Peak SSA/STS # density (#/m3) : ',
     &                               MAXVAL(NDENS_AER(:,:,:,I_SLA))
         WRITE(6,'(a,1x,E12.4)') 'Peak SSA/STS density (kg/m3)  : ',
     &                               MAXVAL(RHO_AER(:,:,:,I_SLA))
         WRITE(6,'(a,1x,E12.4)') 'Peak SSA/STS SAD (cm2/cm3)    : ',
     &                               MAXVAL(SAD_AER(:,:,:,I_SLA))
         WRITE(6,'(a,1x,E12.4)') 'Peak SSA/STS radius (cm)      : ',
     &                               MAXVAL(RAD_AER(:,:,:,I_SLA))
      ENDIF

      END SUBROUTINE CALC_STRAT_AER
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kg_strat_aer
!
! !DESCRIPTION: Function KG\_STRAT\_AER returns the calculated mass of a  
!  stratospheric aerosol. The routine is essentially just an
!  interface to allow external routines to "see" the arrays.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION KG_STRAT_AER (I,J,L,IAER)
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: I,J,L      ! Grid indices
      INTEGER,INTENT(IN)          :: IAER       ! Aerosol index
                                                ! 1 = SSA (pure H2SO4)
                                                ! 2 = STS
                                                ! 3 = Solid PSC
!
! !REVISION HISTORY: 
!  18 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! KG_STRAT_AER begins here!
      !=================================================================

      KG_STRAT_AER = KG_AER(I,J,L,IAER)

      END FUNCTION KG_STRAT_AER
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rho_strat_aer
!
! !DESCRIPTION: Function RHO\_STRAT\_AER returns the calculated
!  stratospheric aerosol mass density.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION RHO_STRAT_AER (I,J,L,IAER)
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: I,J,L      ! Grid indices
      INTEGER,INTENT(IN)          :: IAER       ! Aerosol index:
                                                ! 1 = Liquid aerosol
                                                ! 2 = Solid aerosol
!
! !REVISION HISTORY: 
!  18 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! RHO_STRAT_AER begins here!
      !=================================================================

      RHO_STRAT_AER = RHO_AER(I,J,L,IAER)

      END FUNCTION RHO_STRAT_AER
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_strat_opt
!
! !DESCRIPTION: Subroutine GET\_STRAT\_OPT returns local optical properties 
!  for a given stratospheric aerosol. The routine is essentially just an
!  interface to allow external routines to "see" the arrays. However, local
!  aerosol radius is adjusted from liquid to effective radius for aerosol
!  optical depth calculations with liquid aerosols.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_STRAT_OPT (I,J,L,IAER,RAER,REFF,SAD)
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: I,J,L      ! Grid indices
      INTEGER,INTENT(IN)          :: IAER       ! Aerosol index
                                                ! 1 = SSA (pure H2SO4)
                                                ! 2 = STS
                                                ! 3 = Solid PSC
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT)         :: REFF       ! Effective radius (cm)
      REAL*8, INTENT(OUT)         :: RAER       ! Physical radius (cm)
      REAL*8, INTENT(OUT)         :: SAD        ! Surface area density (cm2/cm3)
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  17 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! GET_STRAT_OPT begins here!
      !=================================================================

      SAD  = SAD_AER(I,J,L,IAER)
      ! For SLA, convert liquid radius to effective optical radius
      RAER = RAD_AER(I,J,L,IAER)
      IF (IAER.eq.I_SLA) THEN
         REFF = RAER/SLA_RR
      ELSE
         REFF = RAER
      ENDIF

      END SUBROUTINE GET_STRAT_OPT
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_nat
!
! !DESCRIPTION: Subroutine CHECK\_NAT determines whether the solid PSC is 
!  composed of ice or NAT (needed for heterogeneous chemistry), or indeed
!  if there is any direct PSC calculation at all. This is important for
!  determining whether to use the JPP or Kirner scheme for ice cloud radii.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_NAT (I,J,L,IS_NAT,IS_PSC)
!
! !USES:
!
      USE TRACERID_MOD,      ONLY : IDTNIT
      USE TRACER_MOD,        ONLY : STT
      USE CHEMGRID_MOD,      ONLY : ITS_IN_THE_STRAT
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: I,J,L      ! Grid indices
!
! !OUTPUT VARIABLES:
!
      LOGICAL,INTENT(OUT)         :: IS_NAT     ! Is surface NAT?
      LOGICAL,INTENT(OUT)         :: IS_PSC     ! Are there solid PSCs?
!
! !REVISION HISTORY: 
!  17 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! CHECK_NAT begins here!
      !=================================================================

      IS_PSC = ((STATE_PSC(I,J,L).ge.2).and.(ITS_IN_THE_STRAT(I,J,L)))
      IS_NAT = (IS_PSC).and.(STT(I,J,L,IDTNIT).gt.TINY(1d0))

      END SUBROUTINE CHECK_NAT
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ternary
!
! !DESCRIPTION: Subroutine TERNARY calculates the composition of SSA/STS 
!  aerosols using a paramaterization from Carslaw et al. "A Thermodynamic
!  Model of the System HCl-HNO3-H2SO4-H2O, Including Solubilities of HBr,
!  from <200 to 328 K". The bulk of this code was taken directly from the
!  Global Modeling Initiative implementation by David Considine.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TERNARY (PCENTER_IN,TCENTER_IN,H2OSUM_IN,H2SO4SUM,
     &                 HNO3SUM,HClSUM,HOClSUM,HBRSum,HOBrSUM,
     &                 W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr,
     &                 HNO3GASFRAC,HClGASFRAC,HOClGASFRAC,
     &                 HBrGASFRAC,HOBrGASFRAC,SLA_VOL,SLA_RHO)
!
! !USES:
!
      USE CMN_GCTM_MOD,      ONLY : PI
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN)  :: PCENTER_IN   ! Pressure (hPa)
      REAL*8, INTENT(IN)  :: TCENTER_IN   ! Temperature (K)
      REAL*8, INTENT(IN)  :: H2OSUM_IN    ! Total H2O mixing ratio
      REAL*8, INTENT(IN)  :: H2SO4SUM     ! Liquid H2SO4 mixing ratio
      REAL*8, INTENT(IN)  :: HNO3SUM      ! Total HNO3 mixing ratio
      REAL*8, INTENT(IN)  :: HClSUM       ! Total HCl mixing ratio
      REAL*8, INTENT(IN)  :: HOClSUM      ! Total HOCl mixing ratio
      REAL*8, INTENT(IN)  :: HBrSUM       ! Total HBr mixing ratio
      REAL*8, INTENT(IN)  :: HOBrSUM      ! Total HOBr mixing ratio
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT) :: W_H2SO4      ! kg H2SO4/kg SLA
      REAL*8, INTENT(OUT) :: W_H2O        ! kg H2O  /kg SLA
      REAL*8, INTENT(OUT) :: W_HNO3       ! kg HNO3 /kg SLA
      REAL*8, INTENT(OUT) :: W_HCl        ! kg HCl  /kg SLA
      REAL*8, INTENT(OUT) :: W_HOCl       ! kg HOCl /kg SLA
      REAL*8, INTENT(OUT) :: W_HBr        ! kg HBr  /kg SLA
      REAL*8, INTENT(OUT) :: W_HOBr       ! kg HOBr /kg SLA
      REAL*8, INTENT(OUT) :: HNO3GASFRAC  ! Gas fraction HNO3     
      REAL*8, INTENT(OUT) :: HClGASFRAC   ! Gas fraction HCl     
      REAL*8, INTENT(OUT) :: HOClGASFRAC  ! Gas fraction HOCl     
      REAL*8, INTENT(OUT) :: HBrGASFRAC   ! Gas fraction HBr     
      REAL*8, INTENT(OUT) :: HOBrGASFRAC  ! Gas fraction HOBr     
      REAL*8, INTENT(OUT) :: SLA_VOL      ! Aerosol volume (m3/m3)
      REAL*8, INTENT(OUT) :: SLA_RHO      ! Aer. mass density (kg/m3)
!
! !REVISION HISTORY: 
!  19 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Derived inputs
      REAL*8 :: H2OSUM
      REAL*8 :: TCENTER
      REAL*8 :: PCENTER

      ! Partial pressures
      REAL*8 :: PATMH2O
      REAL*8 :: PATMHNO3
      REAL*8 :: PATMHCl
      REAL*8 :: PATMHOCl
      REAL*8 :: PATMHBr
      REAL*8 :: PATMHOBr

      ! Molar densities (mol/m3)
      REAL*8 :: MOLDENS_H2SO4

      ! Mass totals
      REAL*8 :: M_H2SO4,M_HNO3
      REAL*8 :: M_HCl,M_HOCl
      REAL*8 :: M_HBr,M_HOBr

      ! Binary solutions denoted with BIN
      ! Mole fractions
      REAL*8 :: X_H2SO4_BIN
      REAL*8 :: X_HNO3_BIN

      ! Mass fractions
      REAL*8 :: M_H2SO4_BIN
      REAL*8 :: M_HNO3_BIN

      ! Effective Henry's Law coefficients
      REAL*8 :: H_H2SO4_BIN
      REAL*8 :: H_HNO3_BIN
      REAL*8 :: H_HCL, H_HOCL
      REAL*8 :: H_HBr, H_HOBr

      ! Frost point
      REAL*8 :: T_ICE

      ! Equilibrium vapor pressure
      REAL*8 :: PVAP_HNO3
      REAL*8 :: PVAP_HCl
      REAL*8 :: PVAP_HBr
      REAL*8 :: PVAP_HOBr

      ! R in m3.atm/(mol K)
      REAL*8, PARAMETER :: R_ATM = 8.205d-5

      ! Transitional variables
      REAL*8 :: DENSITY
      REAL*8 :: TEMPERATURE
      REAL*8 :: PR
      REAL*8 :: TR
      REAL*8 :: TT

      ! Coefficients (q, k) for calculation of H* for H2SO4 and HNO3
      REAL*8,DIMENSION(10) :: QN,QS
      REAL*8,DIMENSION(7)  :: KN,KS

      ! Derived parameters
      REAL*8 :: A,B,C,PHI
      
! ======================================================================
        DATA QN/14.5734d0,0.0615994d0,-1.14895d0,0.691693d0,-0.098863d0,
     &  0.0051579d0,0.123472d0,-0.115574d0,0.0110113d0,0.0097914d0/
        DATA QS/14.4700d0,0.0638795d0,-3.29597d0,1.778224d0,-0.223244d0,
     &  0.0086486d0,0.536695d0,-0.335164d0,0.0265153d0,0.0157550d0/
        DATA KN/-39.136d0,6358.4d0,83.29d0,-17650.0d0,198.53d0,
     &  -11948.d0,-28.469d0/
        DATA KS/-21.661d0,2724.2d0,51.81d0,-15732.0d0,47.004d0,
     &       -6969.0d0,-4.6183d0/
! ======================================================================

      SAVE QN, QS, KN, KS

      !=================================================================
      ! TERNARY begins here!
      !=================================================================

      ! Routine only valid for certain limits
      H2OSUM = MAX(H2OSUM_IN,5.0d-7)
      PCENTER = MAX(PCENTER_IN,5.0d0)
      TCENTER = TCENTER_IN

      ! Calculate partial pressure of H2O & HNO3
      ! PCENTER is in hPa
      PATMH2O  = H2OSUM  * PCENTER / 1013.0d0
      PATMHNO3 = HNO3SUM * PCENTER / 1013.0d0
      PATMHCl  = HClSUM  * PCENTER / 1013.0d0
      PATMHOCl = HOClSUM * PCENTER / 1013.0d0
      PATMHBr  = HBrSUM  * PCENTER / 1013.0d0
      PATMHOBr = HOBrSUM * PCENTER / 1013.0d0

      ! Moles of H2SO4 per m3 air
      MOLDENS_H2SO4 = 100.d0*PCENTER*H2SO4SUM/(8.314d0*TCENTER)

      ! Nucleation temperature of ice
      T_ICE = 2668.70d0/
     &   (10.4310d0-(LOG(PATMH2O)+LOG(760.0d0))/LOG(10.0d0))

      ! Pressure relation
      PR = LOG(PATMH2O)+18.4d0
      
      ! Therefore if temperature lower, set to T_ICE-3
      IF (TCENTER .lt. (T_ICE-3.0d0)) THEN
         TCENTER = (T_ICE-3.0d0)
      ENDIF

      IF (TCENTER .lt. 185.0d0) THEN
         TCENTER = 185.0d0
      ENDIF

      ! ??
      TT = TCENTER * R_ATM * MOLDENS_H2SO4
      
      ! Temperature relation
      TR = 1.0d4/TCENTER-43.4782608d0

      ! Determine H2SO4/H2O pure solution concentration
      ! Mole fraction of H2SO4 in binary solution
      X_H2SO4_BIN = 1.0d0/(2.0d0*(KS(3)+KS(4)/TCENTER))*
     &   (-KS(1)-KS(2)/TCENTER-((KS(1)+KS(2)/TCENTER)**
     &   2.0d0-4.0d0*(KS(3)+KS(4)/TCENTER)*(KS(5)+KS(6)/
     &   TCENTER+KS(7)*LOG(TCENTER)-LOG(PATMH2O)))**0.5d0)

      ! Molality (mol H2SO4/kg H2O) in binary solution
      M_H2SO4_BIN = 55.51d0*X_H2SO4_BIN/(1.0d0-X_H2SO4_BIN)

      IF ((TCENTER.le.215.0d0).AND.(PATMHNO3.gt.0.0d0)) THEN
         ! Determine HNO3/H2SO4/H2O solution composition
         H_H2SO4_BIN = EXP(QS(1)+QS(2)*TR**2+(QS(3)+QS(4)*TR+
     &      QS(5)*TR**2+QS(6)*TR**3)*PR+(QS(7)+QS(8)*TR+
     &      QS(9)*TR**2)*PR**2+QS(10)*TR*PR**3)
         X_HNO3_BIN=1.0d0/(2.0d0*(KN(3)+KN(4)/TCENTER))*
     &      (-KN(1)-KN(2)/TCENTER-((KN(1)+KN(2)/TCENTER)**
     &      2-4.0d0*(KN(3)+KN(4)/TCENTER)*(KN(5)+
     &      KN(6)/TCENTER+KN(7)*LOG(TCENTER)-LOG(PATMH2O)
     &      ))**0.5d0)
         M_HNO3_BIN=55.51d0*X_HNO3_BIN/(1.0d0-X_HNO3_BIN)
         H_HNO3_BIN=EXP(QN(1)+QN(2)*TR**2+(QN(3)+QN(4)*TR+QN(5)*
     &      TR**2+QN(6)*TR**3)*PR+(QN(7)+QN(8)*TR+QN(9)*TR**2)*
     &      PR**2+QN(10)*TR*PR**3)
         A=(TT*H_HNO3_BIN*M_HNO3_BIN**2-TT*H_H2SO4_BIN*M_HNO3_BIN*
     &      M_H2SO4_BIN-2.0d0*M_HNO3_BIN**2*M_H2SO4_BIN+
     &      M_HNO3_BIN*M_H2SO4_BIN**2+H_HNO3_BIN*M_HNO3_BIN*
     &      M_H2SO4_BIN*PATMHNO3-H_H2SO4_BIN*M_H2SO4_BIN**2*
     &      PATMHNO3)/(M_HNO3_BIN**2-M_HNO3_BIN*M_H2SO4_BIN)
         B=M_H2SO4_BIN*(-2.0d0*TT*H_HNO3_BIN*M_HNO3_BIN+TT*
     &      H_H2SO4_BIN*M_H2SO4_BIN+M_HNO3_BIN*M_H2SO4_BIN-
     &      H_HNO3_BIN*M_H2SO4_BIN*PATMHNO3)/(M_HNO3_BIN-
     &      M_H2SO4_BIN)
         C=(TT*H_HNO3_BIN*M_HNO3_BIN*M_H2SO4_BIN**2)/
     &      (M_HNO3_BIN-M_H2SO4_BIN)
         PHI=ATAN(SQRT(4.0d0*(A**2-3.0d0*B)**3-(-2.0d0*A**3+
     &      9.0d0*A*B-27.0d0*C)**2)/(-2.0d0*A**3+9.0d0*A*B-
     &      27.0d0*C))
         IF (PHI.lt.0.d0) THEN
            PHI = PHI + PI
         ENDIF
         M_H2SO4=-1.0d0/3.0d0*(A+2.0d0*SQRT(A**2-3.0d0*B)*
     &      COS((PI+PHI)/3.0d0))
         M_HNO3=M_HNO3_BIN*(1.0d0-M_H2SO4/M_H2SO4_BIN)
         W_H2SO4 = M_H2SO4*0.098076d0/(1.0d0+M_H2SO4*
     &      0.098076d0+M_HNO3*0.063012d0)

         ! Check for low H2SO4
         IF (M_H2SO4 .lt. 0d0) THEN
             M_H2SO4 = 0.0d0
             M_HNO3 = M_HNO3_BIN
             W_H2SO4 = 0.0d0
         ENDIF

         PVAP_HNO3=M_HNO3/(H_HNO3_BIN*M_HNO3/(M_HNO3+
     &      M_H2SO4)+H_H2SO4_BIN*M_H2SO4/(M_HNO3+M_H2SO4))
         W_HNO3 = M_HNO3*0.063012d0/(1.0d0+M_H2SO4*
     &      0.098076d0+M_HNO3*0.063012d0)
         HNO3GASFRAC=(1.0d0-(PATMHNO3-PVAP_HNO3)/PATMHNO3)
         M_HNO3=MAX(M_HNO3,0.0d0)
      ELSE
         ! Solution is pure H2SO4/H2O
         M_H2SO4 = M_H2SO4_BIN
         M_HNO3 = 0.0d0
         W_H2SO4 = M_H2SO4_BIN*0.098076d0/(1.0d0+M_H2SO4_BIN*
     &      0.098076d0)
         W_HNO3 = 0.0d0
         PVAP_HNO3 = 0.0d0
         HNO3GASFRAC=1.0d0
      ENDIF

      ! Handle HCl (Luo et al., Vapor pressures of
      ! H2SO4/HNO3/HCl/HBr/H2O solutions to low stratospheric
      ! temperatures, 1995)
      IF (PATMHCL.gt.0.0d0) THEN
         H_HCL = EXP(-(21.0d0+46.610d0*W_HNO3+4.0690d0*W_H2SO4-
     &      4.8370d0*SQRT(W_HNO3)+2.1860d0*SQRT(W_H2SO4)-63.0d0*
     &      W_HNO3**2-40.170d0*W_HNO3*W_H2SO4-1.5710d0*W_H2SO4**2)-
     &      1.0d0/TCENTER*(-7437.0d0-8327.80d0*W_HNO3+1300.90d0*
     &      W_H2SO4+1087.20d0*SQRT(W_HNO3)-242.710d0*SQRT(W_H2SO4)+
     &      18749.0d0*W_HNO3**2+18500.0d0*W_HNO3*W_H2SO4+
     &      5632.0d0*W_H2SO4**2)-LOG(W_HNO3+0.610d0*W_H2SO4)-
     &      LOG(36.461d0/(1000.0d0+98.076d0*M_H2SO4+63.012d0*
     &      M_HNO3)))*1013.0d0
         M_HCl = (1.0d0/R_ATM/TCENTER*PATMHCL)/
     &      (MOLDENS_H2SO4/M_H2SO4 + 1.0d0/R_ATM/TCENTER/H_HCL)
         W_HCL=M_HCL*36.461d0/(1.0d3+98.076d0*M_H2SO4+
     &      63.012d0*M_HNO3)
         PVAP_HCl = M_HCl/H_HCl
         HClGASFRAC=1.0d0-(PATMHCL-PVAP_HCL)/PATMHCL
      ELSE
         W_HCl=0.0d0
         HClGASFRAC=1.0d0
      ENDIF

      ! Now HOCl
      IF (PATMHOCl.gt.0.0d0) THEN
         H_HOCl=EXP(6.49460d0-(-0.041070d0+54.56d0/TCENTER)*
     &      (M_H2SO4+M_HNO3)-5862.0d0*(1.0d0/298.15d0-1.0d0/
     &      TCENTER))
         M_HOCl=(1.0d0/R_ATM/TCENTER*PATMHOCl)/
     &      (MOLDENS_H2SO4/M_H2SO4 + 1.0d0/R_ATM/TCENTER/H_HOCL)
         W_HOCL=M_HOCL*52.46d0/(1.0d3+98.076d0*M_H2SO4+
     &      63.012d0*M_HNO3)
         ! Realistically expect no gas phase removal
         HOCLGASFRAC=1.0d0 
      ELSE
         W_HOCl=0.0d0
         HOClGASFRAC=1.0d0
      ENDIF
 
      ! Now HBr (Luo et al., Vapor pressures of
      ! H2SO4/HNO3/HCl/HBr/H2O solutions to low stratospheric
      ! temperatures, 1995)
      IF (PATMHBr.gt.0.0d0) THEN
         H_HBr = EXP(-(17.83d0+1.02d0*W_HNO3-1.08d0*W_H2SO4+
     &      3.9d0*SQRT(W_HNO3)+4.38d0*SQRT(W_H2SO4)-8.87d0*
     &      W_HNO3**2-17.0d0*W_HNO3*W_H2SO4+3.73d0*W_H2SO4**2)-
     &      1.0d0/TCENTER*(-8220.5d0-362.76d0*W_HNO3+658.93d0*
     &      W_H2SO4-914.0d0*SQRT(W_HNO3)-955.3d0*SQRT(W_H2SO4)+
     &      9976.6d0*W_HNO3**2+19778.5d0*W_HNO3*W_H2SO4+
     &      7680.0d0*W_H2SO4**2)-LOG(W_HNO3+0.410d0*W_H2SO4)-
     &      LOG(36.461d0/(1000.0d0+98.076d0*M_H2SO4+63.012d0*
     &      M_HNO3)))*1013.0d0
         M_HBr = (1.0d0/R_ATM/TCENTER*PATMHBr)/
     &      (MOLDENS_H2SO4/M_H2SO4 + 1.0d0/R_ATM/TCENTER/H_HBr)
         W_HBr=M_HBr*80.91d0/(1.0d3+98.076d0*M_H2SO4+
     &      63.012d0*M_HNO3)
         PVAP_HBr = M_HBr/H_HBr
         HBrGASFRAC=1.0d0-(PATMHBr-PVAP_HBr)/PATMHBr
      ELSE
         W_HBr=0.0d0
         HBrGASFRAC=1.0d0
      ENDIF

      ! Finally HOBr (Hanson and Ravishankara, Heterogeneous
      ! chemistry of Bromine species in sulfuric acid under
      ! stratospheric conditions, 1995)
      IF (PATMHOBr.gt.0.0d0) THEN
         ! Hanson and Ravishankara state that the volume-based
         ! Henry's Law coefficient for HOBr in H2SO4 is 10^6 M/atm.
         ! The molality-based Henry's law constant, H_HOBr, is
         ! therefore:
         H_HOBr = (1.0d6) * MOLDENS_H2SO4 / M_H2SO4
         M_HOBr = (1.0d0/R_ATM/TCENTER*PATMHOBr)/
     &      (MOLDENS_H2SO4/M_H2SO4 + 1.0d0/R_ATM/TCENTER/H_HOBr)
         W_HOBr = M_HOBr*96.911d0/(1.0d3+98.076d0*M_H2SO4+
     &      63.012d0*M_HNO3)
         PVAP_HOBr = M_HOBr/H_HOBr
         HOBrGASFRAC=1.0d0-(PATMHOBr-PVAP_HOBr)/PATMHOBr
      ELSE
         W_HOBr=0d0
         HOBrGASFRAC=1.0d0
      ENDIF

      ! Take W_H2O as remainder
      W_H2O = 1.d0-(W_H2SO4+W_HNO3+W_HCl+W_HOCl+W_HBr+W_HOBr)

      ! Aerosol mass density in kg/m3 aerosol
      SLA_RHO = CARSLAW_DENSITY(M_H2SO4,M_HNO3,TCENTER)

      ! Aerosol volume in m3/m3 air
      SLA_VOL = (MOLDENS_H2SO4*98.076d0/W_H2SO4/SLA_RHO)*1.d-3

      END SUBROUTINE TERNARY
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carslaw_density
!
! !DESCRIPTION: Function CARSLAW\_DENSITY determines the density of a
!  sol'n through a relationship from Carslaw et al.. Result is in kg/m3.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION CARSLAW_DENSITY(CS,CN,T)
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN)          :: CS         ! H2SO4 molality (mol H2SO4/kg solvent)
      REAL*8, INTENT(IN)          :: CN         ! HNO3 molality (mol HNO3/kg solvent)
      REAL*8, INTENT(IN)          :: T          ! Temperature (K)
!
! !REVISION HISTORY: 
!  19 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                :: DENSS,DENSN

      !=================================================================
      ! CARSLAW_DENSITY begins here!
      !=================================================================

        DENSS=1000.0d0+123.64d0*CS-5.6d-4*CS*T**2  
     &       -29.54d0*CS**1.5d0 + 1.814d-4*CS**1.5d0*T**2  
     &       + 2.343d0*CS**2  -1.487d-3*CS**2*T  
     &       -1.324d-5*CS**2*T**2

        DENSN=1000.0d0+85.107d0*CN-5.043d-4*CN*T**2  
     &       -18.96d0*CN**1.5d0 + 1.427d-4*CN**1.5d0*T**2  
     &       + 1.458d0*CN**2  -1.198d-3*CN**2*T  
     &       -9.703d-6*CN**2*T**2

        CARSLAW_DENSITY=1.0d0/((1.0d0/DENSS*CS/(CS+CN)  
     &  +1.0d0/DENSN*CN/(CS+CN)))
      RETURN

      END FUNCTION CARSLAW_DENSITY
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_fallvel
!
! !DESCRIPTION: Function CALC\_FALLVEL calculates the terminal velocity of a
!  solid particle.
!\\
!\\
! !INTERFACE:
!
      FUNCTION CALC_FALLVEL(DENSITY,RADIUS,TCENTER,PCENTER) RESULT(VEL)
!
! !USES:
!
      USE CMN_GCTM_MOD,    ONLY : PI, g0
!
! !INPUT PARAMETERS:
!
      REAL*8,INTENT(IN)   :: RADIUS  ! Particle radius (cm)
      REAL*8,INTENT(IN)   :: DENSITY ! Particle density (kg/m3)
      REAL*8,INTENT(IN)   :: TCENTER ! Local temperature (K)
      REAL*8,INTENT(IN)   :: PCENTER ! Local pressure (kPa)
!
! !OUTPUT VARIABLES:
!
      REAL*8              :: VEL ! Fall velocity (m/s)
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  11 Aug 2012 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8              :: Vy ! Intermediate velocity (m/s)
      REAL*8,PARAMETER    :: eta=6.45d-8 ! Constant (kg/(msK))
      REAL*8              :: val_x ! Dimensionless variable
      REAL*8,DIMENSION(3) :: alpha ! Auxiliary variables
      REAL*8              :: PR ! Pressure times radius

      !=================================================================
      ! CALC_FALLVEL begins here!
      !=================================================================

      DATA ALPHA/1.49d-5,5.02d-6,2.64d-5/

      ! Sanity check
      IF ((RADIUS.le.0.d0).or.(DENSITY.le.0.d0)) THEN
         VEL=0.d0
      ELSE
         ! PCENTER (kPa -> Pa) = *1.d3
         ! RADIUS  (cm  -> m ) = *1.d-2
         ! Therefore multiply PR by 10
         PR = PCENTER * RADIUS * 10d0
         VAL_X = -1.0d0*PR/(ALPHA(3)*TCENTER)
         VAL_X = ALPHA(2)*TCENTER*EXP(VAL_X)/PR
         VAL_X = 1.0d0 + VAL_X + (ALPHA(1)*TCENTER/PR)
         Vy = g0*DENSITY*RADIUS*RADIUS*(1.d-4)/(4.5*ETA*TCENTER)
         VEL = 0.893d0 * Vy * VAL_X
      ENDIF

      END FUNCTION CALC_FALLVEL
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cacl_sla_gamma
!
! !DESCRIPTION: Subroutine CALC\_SLA\_GAMMA calculates 11 different sticking 
!  coefficients on the surface of local stratospheric liquid aerosols, 
!  relevant to each of the 11 reactions listed in Kirner's paper.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_SLA_GAMMA( NDENS, T, P, WT_FRC, H2OSUM, HClSUM,
     &                           HBrSUM, HOBrSUM, ClNO3SUM, BrNO3SUM,
     &                           RHO, ARAD, RXNGAMMA )
!
! !USES:
!
      USE CMN_GCTM_MOD,       ONLY : PI
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN)  :: NDENS    ! Air number density (molec/cm3)
      REAL*8, INTENT(IN)  :: T        ! Temperature (K)
      REAL*8, INTENT(IN)  :: P        ! Pressure (hPa)
      REAL*8, INTENT(IN)  :: WT_FRC   ! Weight fraction of H2SO4 (kg/kg)
      REAL*8, INTENT(IN)  :: H2OSUM   ! H2O mixing ratio
      REAL*8, INTENT(IN)  :: HClSUM   ! HCl mixing ratio
      REAL*8, INTENT(IN)  :: HBrSUM   ! HBr mixing ratio
      REAL*8, INTENT(IN)  :: HOBrSUM  ! HOBr mixing ratio
      REAL*8, INTENT(IN)  :: ClNO3SUM ! ClNO3 mixing ratio
      REAL*8, INTENT(IN)  :: BrNO3SUM ! BrNO3 mixing ratio
      REAL*8, INTENT(IN)  :: RHO      ! STS density (g/cm3)
      REAL*8, INTENT(IN)  :: ARAD     ! SLA radius (cm)
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT) :: RXNGAMMA(11) ! Premultiplying factors
!
! !REVISION HISTORY: 
!  10 Oct 2012 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8               :: WT      ! Weight percentage H2SO4 (100*kg/kg)
      REAL*8               :: H2OPP   ! Partial pressure of H2O (hPa)
      REAL*8               :: HClPP   ! Partial pressure of HCl (atm)
      REAL*8               :: HBrPP   ! Partial pressure of HBr (atm)
      REAL*8               :: HOBrPP  ! Partial pressure of HOBr (atm)
      REAL*8               :: ClNO3PP ! Partial p. of ClONO2 (atm)
      REAL*8               :: BrNO3PP ! Partial p. of BrONO2 (atm)
      REAL*8               :: PSATH2O ! Water vapor sat. pressure (hPa)
      REAL*8               :: ACTH2O  ! Activity of water
      REAL*8               :: MOLAL   ! Molality of H2SO4 (mol H2SO4/kg solvent)
      REAL*8, DIMENSION(3) :: Z       ! Parameters for H2SO4 sol'n
      REAL*8               :: M_H2SO4 ! Mass of H2SO4

      ! HOBr parameters
      REAL*8               :: c_HOBr
      REAL*8               :: SHOBr
      REAL*8               :: HHOBr
      REAL*8               :: DHOBr
      REAL*8               :: kHOBr_HCl
      REAL*8               :: GHOBrrxn
      REAL*8               :: lHOBr
      REAL*8               :: fHOBr
      REAL*8               :: gHOBr_HCl

      ! HOCl parameters
      REAL*8               :: c_HOCl
      REAL*8               :: SHOCl
      REAL*8               :: HHOCl
      REAL*8               :: DHOCl
      REAL*8               :: kHOCl_HCl
      REAL*8               :: GHOClrxn
      REAL*8               :: lHOCl
      REAL*8               :: fHOCl
      REAL*8               :: gHOCl_HCl

      ! ClNO3 parameters
      REAL*8               :: c_ClNO3
      REAL*8               :: SClNO3
      REAL*8               :: HClNO3
      REAL*8               :: DClNO3
      REAL*8               :: GClNO3rxn
      REAL*8               :: lClNO3
      REAL*8               :: fClNO3
      REAL*8               :: gClNO3
      REAL*8               :: gClNO3_HCl
      REAL*8               :: gClNO3_H2O
 
      ! N2O5 parameters
      REAL*8, DIMENSION(3) :: AK

      ! Other parameters
      REAL*8 :: kH2O,kH,khdr,GbH2O,HHCl,MHCl,kHCl,GbHCl,Gs,FHCl,Gsp
      REAL*8 :: GbHClp, Gb, khydr, kII, k_dl

      ! Interim variables
      REAL*8               :: X,A,H,T_THRESHOLD,aH
      REAL*8, PARAMETER    :: MAX_T_DIFF = 6.0d0

      ! Control whether to run HCl calculations
      LOGICAL              :: HClOK

      !=================================================================
      ! CALC_SLA_GAMMA begins here!
      !=================================================================

      PSATH2O = EXP(18.452406985d0-3505.1578807d0/T-330918.55082d0/(T*T)
     &     +12725068.262d0/(T*T*T))       ! Saturation pressure of H2O
      H2OPP = H2OSUM * P                  ! Partial pressure of H2O
      ACTH2O = MAX((H2OPP/PSATH2O),1.0d0) ! Water activity

      ! Calculate molality of solution
      WT = 100.0d0*WT_FRC ! Convert from fraction to %
      MOLAL = 1000.0d0 * (WT/98.0d0/(100.0-WT))

      ! Parameters for H2SO4 solution
      !----------------------------------------------------------
      ! The solution density is calculated earlier, including
      ! contributions from HNO3. This code treats it as a binary
      ! solution - so far this is just a kludge. Need to update
      ! all this code to acknowledge the presence of at least
      ! HNO3 (e.g. X is still calculated based on pure H2O
      ! solvent!)
      !----------------------------------------------------------
      !Z(1) =   0.12364d0-5.6d-7*T*T
      !Z(2) =   -0.02954d0+1.814d-7*T*T
      !Z(3) =   2.343d-3-1.487d-6*T-1.324d-8*T*T
      !RHO  =   1.0d0+Z(1)*MOLAL+Z(2)*MOLAL**1.5+Z(3)*MOLAL*MOLAL
      !----------------------------------------------------------
      M_H2SO4 =   RHO*WT/9.8 ! Molality (mol H2SO4/kg solvent)
      X    =   WT/(WT+(100.-WT)*98./18.)
      A    =   169.5+5.18*WT-0.0825*WT*WT+3.27d-3*WT*WT*WT
      T_THRESHOLD = 144.11+0.166*WT-0.015*WT*WT+2.18d-4*WT*WT*WT
      IF ((T-T_THRESHOLD).gt.MAX_T_DIFF) THEN
         H = A*T**(-1.43)*EXP(448./(T-T_THRESHOLD))
      ELSE
         H = A*T**(-1.43)*EXP(448./MAX_T_DIFF)
      ENDIF

      aH   =   EXP(60.51-0.095*WT+0.0077*WT*WT-1.61e-5*WT*WT*WT
     &     -(1.76+2.52d-4*WT*WT)*SQRT(T) +
     &     (-805.89+253.05*WT**0.076)/SQRT(T))

      HClPP = HClSUM*P/1013.25 ! Note atm, not hPa
      ClNO3PP = ClNO3SUM*P/1013.25
      BrNO3PP = BrNO3SUM*P/1013.25
 
      ! Should we bother running HCl calculations?
      HClOK = (HClPP .gt. 1.d-30)

      ! Reaction 1. N2O5 + H2O (hydrolysis of N2O5)
      AK(1)=-25.5265-0.133188*WT+0.0093084*WT**2-9.0194E-5*WT**3
      AK(2)=9283.76+115.345*WT-5.19258*WT**2+0.0483464*WT**3
      AK(3)=-851801-22191.2*WT+766.916*WT**2-6.85427*WT**3
      RXNGAMMA(1)=exp(AK(1)+AK(2)/T+AK(3)/T**2)
      
      ! Reaction 2. N2O5 + HCl
      ! JPL 10-06 suggests near-zero gamma
      RXNGAMMA(2) = TINY(1d0)

      ! Reactions 3/4. ClNO3 + H2O/HCl
      c_ClNO3     =  1474.d0*SQRT(T)
      SClNO3      =  0.306d0+24.d0/T
      HClNO3      =  1.6d-6*EXP(4710.D0/T)*EXP(-SClNO3*M_H2SO4)
      DClNO3      =  5d-8*T/h
      kH2O        =  1.95d10*EXP(-2800.d0/T)
      kH          =  1.22d12*EXP(-6200.d0/T)
      khydr       =  kH2O*ACTH2O + kH*aH*ACTH2O
      GbH2O       =  4.D0*HClNO3*0.082D0*T*SQRT(DClNO3*khydr)/c_ClNO3
      HHCl        =  (0.094D0-0.61D0*X+1.2D0*X*X)*EXP(-8.68D0+
     &                                   (8515.D0-10718.D0*X**0.7)/T)
      MHCl        =  HHCl *HClPP
      kHCl        =  7.9D11*aH*DClNO3*MHCl
      lClNO3      =  SQRT(DClNO3/(khydr+kHCl))
      fClNO3      =  1.D0/tanh(arad/lClNO3)- lClNO3/arad
      GClNO3rxn   =  fClNO3*GbH2O *SQRT(1.D0+kHCl/khydr)
      GbHCl       =  GClNO3rxn* kHCl/(kHCl+ khydr)
      Gs          =  66.12D0*EXP(-1374.D0/T)*HClNO3*MHCl
      FHCl        =  1.D0/(1.D0+0.612D0*(Gs+GbHCl)* ClNO3PP/ HClPP)
      Gsp         =  FHCl*Gs
      GbHClp      =  FHCl*GbHCl
      Gb          =  GbHClp  + GClNO3rxn* khydr/( kHCl+ khydr)
      gClNO3      =  1.D0/(1.D0+1.D0/(Gsp + Gb))
      gClNO3_HCl  =  gClNO3 *(Gsp + GbHClp)/(Gsp + Gb)
      gClNO3_H2O  =  gClNO3 - gClNO3_HCl

      RXNGAMMA(3) =  gClNO3_H2O
      IF (HClOK) THEN
         RXNGAMMA(4) =  gClNO3_HCl
      ELSE
         RXNGAMMA(4) = TINY(1d0)
      ENDIF

      ! Reaction 5. ClNO3 + HBr
      ! Not present in JPL 10-06 for H2SO4
      RXNGAMMA(5) = TINY(1d0)

      ! Reaction 6. BrNO3 + H2O
!      RXNGAMMA(6) = 1.0/(1.0/0.88+exp(-17.832+0.245*WT))
      RXNGAMMA(6) = 1.D0/(1.D0/0.80D0+1.D0/(exp(29.2D0-0.4D0*WT )+0.11))
      
      ! Reaction 7. BrNO3 + HCl
      RXNGAMMA(7) = 0.9d0 ! JPL 10-06

      ! Reaction 8. HOCl + HCl
      IF (HClOK) THEN
         c_HOCl    =  MOLEC_SPEED(T,52.46D0)
         SHOCl     =  0.0776D0+59.18D0/T
         HHOCl     =  1.91D-6*EXP(5862.4D0/T)*EXP(-SHOCl*M_H2SO4)
         DHOCl     =  6.4d-8*T/H
         kHOCl_HCl =  1.25d9*aH*DHOCl*MHCl
         lHOCl     =  SQRT(DHOCl/kHOCl_HCl)
         if (lHOCl.gt.(1.d3*arad)) then
            ! Limiting rate
            fHOCl = arad/(3.d0*lHOCl)
         else
            fHOCl     =  1.D0/tanh(arad/lHOCl)- lHOCl/arad
         endif
         GHOClrxn  =  4.D0*HHOCl*0.082D0*T*sqrt(DHOCl*kHOCl_HCl)/c_HOCl
         IF (fHOCl.eq.0.) THEN
            gHOCl_HCl =  TINY(1d0)
         ELSE
            gHOCl_HCl =  1.D0/(1.D0+1.D0/(fHOCl*GHOClrxn*FHCl))
         ENDIF
      ELSE
         gHOCl_HCl = TINY(1d0)
      ENDIF
   
      RXNGAMMA(8) = gHOCl_HCl

      ! Reaction 9. HOCl + HBr
      ! Not yet implemented for STS; JPL 10-06 suggests complex
      ! relationship, not yet sufficiently well understood or
      ! parameterized for the purposes of simulation. Ignore for now
      RXNGAMMA(9) = TINY(1d0)

      ! Reaction 10. HOBr + HCl
      IF (HClOK) THEN
         HOBrPP = HOBrSUM*P/1013.25 ! Note atm, not hPa
         c_HOBr    =  MOLEC_SPEED(T,96.91D0)
         SHOBr     =  0.0776D0+59.18D0/T
!         HHOBr     =  30.D0
         HHOBR     = exp(-9.86D0+5427.D0/T)
         DHOBr     =  1.E-8
         kII       = exp(154.d0-1.63d0*WT)*exp(-(3.85d4-478.d0*WT)/T)
         k_dl      = 7.5D14*(DHOBr*arad*1.e7)
         IF (kII.gt.k_dl) kII=k_dl
         kHOBr_HCl =  kII*HHOBr*HOBrPP
         GHOBrrxn  =  4.D0*HHCl*0.082D0*T*sqrt(DHOBr*kHOBr_HCl)/c_HOBr
         lHOBr     =  sqrt(DHOBr/kHOBr_HCl)
         if (lHOBr.gt.(1.d3*arad)) then
            ! Limiting rate
            fHOBr = arad/(3.d0*lHOBr)
         else
            fHOBr     =  1.D0/tanh(arad/lHOBr)- lHOBr/arad
         endif
         gHOBr_HCl =  1.D0/(1.D0+1.D0/(fHOBr*GHOBrrxn))
         RXNGAMMA(10) = gHOBr_HCl
      ELSE
         RXNGAMMA(10) = TINY(1d0)
      ENDIF

      ! Reaction 11. HOBr + HBr
      ! Data from JPL limited; ignore for now
      RXNGAMMA(11) = TINY(1d0)

      ! Return to calling program
      END SUBROUTINE CALC_SLA_GAMMA
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: molec_speed
!
! !DESCRIPTION: Function MOLEC\_SPEED calculates the mean velocity of gas
!  phase particles based on temperature and molecular mass.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION MOLEC_SPEED(T,MOLMASS)
!
! !USES:
!
      USE CMN_GCTM_MOD, ONLY : PI
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN)          :: T       ! Temperature (K)
      REAL*8, INTENT(IN)          :: MOLMASS ! Molecular mass (g/mol)
!
! !REVISION HISTORY: 
!  10 Oct 2012 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! MOLEC_SPEED begins here!
      !=================================================================

      MOLEC_SPEED=SQRT(8.0d0*8.31451d7*T/(PI*MOLMASS))

      END FUNCTION MOLEC_SPEED
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: make_psc_file
!
! !DESCRIPTION: Subroutine MAKE\_PSC\_FILE writes PSC state information
!  into a checkpoint file (binary punch file format).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MAKE_PSC_FILE( YYYYMMDD, HHMMSS )
!
! !USES:
!
      USE BPCH2_MOD
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE FILE_MOD,    ONLY : IOERROR
      USE GRID_MOD,    ONLY : GET_XOFFSET
      USE GRID_MOD,    ONLY : GET_YOFFSET
      USE LOGICAL_MOD, ONLY : LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE
      USE TIME_MOD,    ONLY : GET_TAU

      USE CMN_SIZE_MOD                   ! Size parameters
      USE COMODE_LOOP_MOD                ! IGAS
      USE CMN_MOD                        ! TAU , NSRCX, LSOILNOX
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: YYYYMMDD   ! YYYY/MM/DD GMT date
      INTEGER, INTENT(IN)  :: HHMMSS     ! hh:mm:ss   GMT time
! 
! !REVISION HISTORY: 
!  28 Apr 2013 - S. D. Eastham - Initial version, based on MAKE_CSPEC_FILE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: IU_RST
      INTEGER              :: I,    I0, IOS, J,  J0, L, N, JLOOP
      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH
      CHARACTER(LEN=255)   :: FILENAME
      CHARACTER(LEN=255)   :: OUTPUT_CHECKPT_FILE


      ! Temporary storage arrays for checkpointed variables
      REAL*4               :: TMP(ILONG, ILAT, IPVERT)

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      ! make HALFPOLAR variable (hotp 2/25/09)
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 


      !=================================================================
      ! MAKE_PSC_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_CHECKPT_FILE = 'restart.psc.YYYYMMDDhh'

      ! Clear some arrays 
      ! use minimum value instead of zero hotp 2/25/09
      TMP(:,:,:)   = 0


      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Checkpoint File: ' // 
     &           'Instantaneous PSC state (unitless)'
      CATEGORY = 'IJ-PSC-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      ! get value of HALFPOLAR hotp 2/25/09
      HALFPOLAR = GET_HALFPOLAR()

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      ! Find a free file LUN
      IU_RST   = findFreeLUN()

      ! Copy the output checkpoint file name into a local variable
      FILENAME = TRIM( OUTPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to filename
      !FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME ), IU_RST
 100  FORMAT( '     - MAKE_PSC_FILE: Writing ', a, ' on unit ', i4 )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each checkpointed quantity to the checkpoint file
      !=================================================================

      ! Checkpt additional values for full chem simulation
      UNIT = '-'
          
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, IPVERT
         DO J = 1, ILAT
         DO I = 1, ILONG
        
              TMP(I,J,L) = REAL(STATE_PSC(I,J,L))+0.1d0

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
 
      ! Save the species name to the RESERVED slot of the bpch file
      RESERVED = 'STATE_PSC'

      ! Write the data block to the CSPEC checkpoint file
      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               ILONG,     ILAT,      IPVERT,    I0+1,
     &               J0+1,      1,         TMP(:,:,:) )
 
      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_PSC_FILE: wrote file' )

      END SUBROUTINE MAKE_PSC_FILE
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_psc_file
!
! !DESCRIPTION: Subroutine READ\_PSC\_FILE initializes PSC state information
!  from a checkpoint file (binary punch file format).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_PSC_FILE( am_I_Root, YYYYMMDD,  
     &                          HHMMSS,    IT_EXISTS ) 
!
! !USES:
!
      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE FILE_MOD,    ONLY : IOERROR
      USE FILE_MOD,    ONLY : FILE_EXISTS
      USE LOGICAL_MOD, ONLY : LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE

      USE CMN_SIZE_MOD                  ! Size parameters
      USE COMODE_LOOP_MOD               ! ITLOOP, IGAS
      USE CMN_MOD                        ! TAU , NSRCX, LSOILNOX
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN) :: am_I_Root  ! Is this the root CPU?
      INTEGER, INTENT(IN) :: YYYYMMDD   ! YYYY/MM/DD GMT date
      INTEGER, INTENT(IN) :: HHMMSS     ! hh:mm:ss   GMT time
! 
! !REVISION HISTORY: 
!  28 Apr 2013 - S. D. Eastham - Initial version, based on READ_CSPEC_FILE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, IOS, IU_RST, J, L, N, NN, NTL
      REAL*4              :: TMP(ILONG,ILAT,IPVERT)
      LOGICAL             :: IT_EXISTS 

      REAL*8              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: INPUT_CHECKPT_FILE
      CHARACTER(LEN=255)  :: MSG
     

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_PSC_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_CHECKPT_FILE = 'restart.psc.YYYYMMDDhh'

      ! Initialize some variables
      TMP(:,:,:) = 0

      !=================================================================
      ! Open checkpoint file and read top-of-file header
      !=================================================================

      ! Find a free file LUN
      IU_RST   = findFreeLUN()
     
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to name
      !FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )

      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), IU_RST
      ENDIF
 100  FORMAT( '     - READ_PSC_FILE: Reading ', a, ' on unit ', i4 )
 
      ! Check to see if PSC restart file exists
      IT_EXISTS = FILE_EXISTS( FILENAME )
      IF ( .not. IT_EXISTS ) THEN 
         RETURN
      ENDIF 

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )

      ! Read the values of STATE_PSC
      READ( IU_RST, IOSTAT=IOS )
     &    MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS < 0 is end-of-file, so exit
      IF ( IOS < 0 ) GOTO 555

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_psc_file:13' )

      ! Read data block header
      READ( IU_RST, IOSTAT=IOS )
     &      CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &      NTL,      NN,       NL,   IFIRST, JFIRST, LFIRST,
     &      NSKIP

      ! Error check
      IF ( IOS /= 0 ) 
     &   CALL IOERROR(IOS,IU_RST,'read_psc_file:14' )

      ! Read data block
      READ( IU_RST, IOSTAT=IOS )
     &    ( ( ( TMP(I,J,L), I= 1, NTL), J=1,NN ), L = 1, NL)

      ! Error check
      IF ( IOS /= 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_psc_file:16' )

      !==============================================================
      ! Assign data from the TMP array to STATE_PSC
      !==============================================================

      ! Only process PSC data 
      IF ( CATEGORY(1:8) == 'IJ-PSC-$' .and.
     &     NTL           == ILONG      .and. 
     &     NN            == ILAT       .and. 
     &     NL            == IPVERT            ) THEN

         DO I=1,NTL
         DO J=1,NN
         DO L=1,NL
         STATE_PSC(I,J,L) = INT(FLOOR(TMP(I,J,L)))
         ENDDO
         ENDDO
         ENDDO

      ELSE
         CALL ERROR_STOP(' Restart data is not correct ', 
     &                   ' reading STATE_PSC, gmistrat_mod')

      ENDIF

 555  CONTINUE

      ! Close file
      CLOSE( IU_RST )      

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_PSC_FILE: read file' )

      END SUBROUTINE READ_PSC_FILE
!
!EOC
!!------------------------------------------------------------------------------
!!               MIT Laboratory for Aviation and the Environment               !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: srtn
!!
!! !DESCRIPTION: Subroutine SRTN is a reference subroutine, with more than one
!!  line of description.
!!\\
!!\\
!! !INTERFACE:
!!
!      SUBROUTINE SRTN (INARG1,IOARG2,OARG3)
!!
!! !USES:
!!
!      USE THIS_MOD, ONLY : THIS_FN
!!
!! !INPUT PARAMETERS:
!!
!      REAL*8, INTENT(IN)          :: INARG1     ! Description of INARG1
!      REAL*8, INTENT(INOUT)       :: IOARG2     ! Description of IOARG2
!!
!! !OUTPUT VARIABLES:
!!
!      REAL*8, INTENT(OUT)         :: OARG3      ! Description of OARG3
!!
!! !REMARKS:
!! (1) A remark
!! 
!! !REVISION HISTORY: 
!!  28 Mar 2013 - S. D. Eastham - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      INTEGER, PARAMETER    :: IPARAM1=2
!
!      !=================================================================
!      ! SRTN begins here!
!      !=================================================================
!
!      OARG3 = INARG1 + IOARG2
!      IOARG2 = OARG3/2.d0*IPARAM1
!
!      END SUBROUTINE SRTN
!!
!!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gmistrat
!
! !DESCRIPTION: Subroutine INIT\_GMISTRAT initializes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GMISTRAT( am_I_Root, YYYYMMDD, HHMMSS )
!
! !USES:
!
      USE TRACERID_MOD
      USE LOGICAL_MOD,  ONLY : LN2OEMIS, LCFCEMIS, LCLEMIS,  LBREMIS
      USE LOGICAL_MOD,  ONLY : LOCSEMIS
      USE ERROR_MOD,    ONLY : ALLOC_ERR
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root  ! Is this the root CPU?
      INTEGER, INTENT(IN) :: YYYYMMDD   ! YYYY/MM/DD GMT date
      INTEGER, INTENT(IN) :: HHMMSS     ! hh:mm:ss   GMT time
!
! !OUTPUT VARIABLES:
!
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N, AS
      LOGICAL :: IT_EXISTS

      !=================================================================
      ! INIT_GMISTRAT begins here!
      !=================================================================

      WRITE( 6,'(a)') REPEAT( '=', 79 )
      WRITE( 6,'(a)') 'U N I F I E D   C H E M I S T R Y'
      WRITE( 6,'(a)') 'Routines written by SEBASTIAN D. EASTHAM'
      WRITE( 6,'(a)') 'Last Modification Date: 05/09/13'
      WRITE( 6,'(a)') REPEAT( '=', 79 )

      ! Calculate conversion factors for SLA
      ! Factor to convert volume (m3 SLA/m3 air) to
      ! surface area density (cm2 SLA/cm3 air)
      SLA_VA = (8.406d-8)*(10.d0**(12.d0*0.751d0))

      ! Factor to convert effective radius to 
      ! liquid radius (unitless)
      SLA_RR = EXP(-0.173d0)

      ! Factor to convert volume (m3/m3) to effective
      ! radius (m)
      SLA_VR = (0.357d-6)*(10.d0**(12.d0*0.249))

      ! Set "current" month to impossible value
      SFC_MONTH = 0

      ! Determine which tracers are to be emitted
      N_EMIT = 0
      IDT_EMIT = 0

      IF (LCFCEMIS) THEN
         ! R-11/12/113/114/115
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFC11
         TRAC_EMIT(N_EMIT) = 'CFC11'

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFC12
         TRAC_EMIT(N_EMIT) = 'CFC12'

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFCX
         TRAC_EMIT(N_EMIT) = 'CFCX'

         ! HCFC-22/123/141b/142b
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTHCFC22
         TRAC_EMIT(N_EMIT) = 'HCFC22'

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTHCFCX
         TRAC_EMIT(N_EMIT) = 'HCFCX'

         ! Halons (H-1301/2402/1211)
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH1301
         TRAC_EMIT(N_EMIT) = 'H1301'

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH1211
         TRAC_EMIT(N_EMIT) = 'H1211'

         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH2402
         TRAC_EMIT(N_EMIT) = 'H2402'
      ENDIF

      IF (LCLEMIS) THEN
         ! CCl4
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCCl4
         TRAC_EMIT(N_EMIT) = 'CCL4'

         ! CH3Cl
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3Cl
         TRAC_EMIT(N_EMIT) = 'CH3CL'

         ! CH3CCl3
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3CCl3
         TRAC_EMIT(N_EMIT) = 'CH3CCL3'
      ENDIF

      IF (LOCSEMIS) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTOCS
         TRAC_EMIT(N_EMIT) = 'OCS'
      ENDIF

      IF (LN2OEMIS) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTN2O
         TRAC_EMIT(N_EMIT) = 'N2O'
      ENDIF 

      IF (LBREMIS) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3Br
         TRAC_EMIT(N_EMIT) = 'CH3BR'
      ENDIF

      IF (N_EMIT.gt.0) THEN
         ALLOCATE( GRID_EMIT( JJPAR, N_EMIT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'GRID_EMIT' )
         GRID_EMIT = 0d0
      ENDIF

      ! Clear other tracer names for safety
      DO N=(N_EMIT+1),MAX_EMIT
         TRAC_EMIT(N) = ''
         IDT_EMIT(N) = 0
      ENDDO

      ! Need to put a trailing slash on the file root
      MONTREAL_FILE_ROOT = TRIM(MONTREAL_FILE_ROOT) // '/'

      ! Initialize PSC variables
      ALLOCATE( RAD_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RAD_AER' )
      RAD_AER = 0d0
      
      ALLOCATE( KG_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'KG_AER' )
      KG_AER = 0d0
      
      ALLOCATE( SAD_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SAD_AER' )
      SAD_AER = 0d0
      
      ALLOCATE( NDENS_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NDENS_AER' )
      NDENS_AER = 0d0
      
      ALLOCATE( RHO_AER( IIPAR, JJPAR, LLPAR, NSTRATAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RHO_AER' )
      RHO_AER = 0d0
      
      ALLOCATE( STATE_PSC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'STATE_PSC' )
      STATE_PSC = 0d0
      
      ! Mass fraction of tracers contained in liquid aerosol
      ! Indices: 1 - SO4
      !          2 - HNO3
      !          3 - HCl
      !          4 - HOCl
      !          5 - HBr
      !          6 - HOBr
      ALLOCATE( AERFRAC( IIPAR, JJPAR, LLPAR, 6 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AERFRAC' )
      AERFRAC = 0d0

      ALLOCATE( AERFRACIND( 6 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AERFRACIND' )
      AERFRACIND(1) = IDTSO4
      AERFRACIND(2) = IDTHNO3
      AERFRACIND(3) = IDTHCl
      AERFRACIND(4) = IDTHOCl
      AERFRACIND(5) = IDTHBr
      AERFRACIND(6) = IDTHOBr
    
      ! Sticking coefficients for PSC reactions on SLA 
      ALLOCATE( KHETI_SLA( IIPAR, JJPAR, LLPAR, 11 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'KHETI_SLA' )
      KHETI_SLA = 0d0
     
      ! Read in PSC state
      CALL READ_PSC_FILE( am_I_Root, YYYYMMDD,  
     &                    HHMMSS,    IT_EXISTS )

      IF ((am_I_Root).and.(.not.IT_EXISTS)) THEN
         WRITE(6,*) 
     &  '    - GMISTRAT: PSC restart not found, initialize PSC-free'
      ENDIF

      END SUBROUTINE INIT_GMISTRAT
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_gmistrat
!
! !DESCRIPTION: Subroutine CLEANUP\_GMISTRAT deallocates module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GMISTRAT
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!  13 Apr 2013 - S. D. Eastham - Added PSC arrays
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_GMISTRAT begins here!
      !=================================================================

      IF ( ALLOCATED( GRID_EMIT  ) ) DEALLOCATE( GRID_EMIT  )
      IF ( ALLOCATED( RAD_AER    ) ) DEALLOCATE( RAD_AER    )
      IF ( ALLOCATED( SAD_AER    ) ) DEALLOCATE( SAD_AER    )
      IF ( ALLOCATED( KG_AER     ) ) DEALLOCATE( KG_AER     )
      IF ( ALLOCATED( RHO_AER    ) ) DEALLOCATE( RHO_AER    )
      IF ( ALLOCATED( STATE_PSC  ) ) DEALLOCATE( STATE_PSC  )
      IF ( ALLOCATED( NDENS_AER  ) ) DEALLOCATE( NDENS_AER  )
      IF ( ALLOCATED( AERFRAC    ) ) DEALLOCATE( AERFRAC    )
      IF ( ALLOCATED( AERFRACIND ) ) DEALLOCATE( AERFRACIND )
      IF ( ALLOCATED( KHETI_SLA  ) ) DEALLOCATE( KHETI_SLA  )

      END SUBROUTINE CLEANUP_GMISTRAT
!
!EOC
      END MODULE GMISTRAT_MOD
