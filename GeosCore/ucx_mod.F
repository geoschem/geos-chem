!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: ucx_mod
!     
! !DESCRIPTION: Module UCX\_MOD contains routines and variables which
!  are associated with the addition of full stratospheric chemistry to
!  GEOS-Chem (based on the NASA GMI implementation, forming the Unified
!  Chemistry eXtension (UCX).
!\\   
!\\   
! !INTERFACE: 
!
      MODULE UCX_MOD
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
      CHARACTER(LEN=255), PUBLIC :: AVG_FILE_ROOT ! Directory for diurnal avg data
      CHARACTER(LEN=255), PUBLIC :: NOON_FILE_ROOT ! Directory for noontime data
      CHARACTER(LEN=255), PUBLIC :: MONTREAL_FILE_ROOT ! Root directoy for WHO data 
                          PUBLIC :: T_STS ! Max temperature of STS formation (K)
      REAL*8,             PUBLIC :: T_NAT_SUPERCOOL ! NAT supercooling (K)
      REAL*8,             PUBLIC :: P_ICE_SUPERSAT ! Ice supersaturation (-)
                          PUBLIC :: KHETI_SLA ! Reaction cofactors for SLA
                          PUBLIC :: NDENS_AER ! See below
      INTEGER,            PUBLIC :: CFCYEAR ! Year for CFC emissions
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: SET_INITIAL_MIXRATIOS
      PUBLIC  :: SET_H2O_TRAC
      PUBLIC  :: SET_CLOCK_TRAC
      PUBLIC  :: EMISS_BASIC
      PUBLIC  :: SETTLE_STRAT_AER
      PUBLIC  :: SO4_PHOTFRAC
      PUBLIC  :: UCX_NOX
      PUBLIC  :: UCX_H2SO4PHOT
      PUBLIC  :: CALC_STRAT_AER
      PUBLIC  :: GET_STRAT_OPT
      PUBLIC  :: KG_STRAT_AER
      PUBLIC  :: RHO_STRAT_AER
      PUBLIC  :: CHECK_NAT
      PUBLIC  :: MAKE_PSC_FILE
      PUBLIC  :: INIT_UCX
      PUBLIC  :: CLEANUP_UCX
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
      ! R_UNIV          : Universal gas constant (J/mol K)
      ! INITMR_BASIS    : Year for which the initializing mixing ratios
      !                   were calculated (needed for future-scaling)
      !
      !=================================================================

      INTEGER, PARAMETER            :: AER_NLEVS=51
      REAL*8, PARAMETER             :: T_STS=240.0d0
      INTEGER, PARAMETER            :: I_SLA=1
      INTEGER, PARAMETER            :: I_SPA=2
      REAL*8,  PARAMETER            :: AVOGADRO=6.022d23
      REAL*8, PARAMETER             :: R_UNIV = 8.31448d0
      INTEGER, PARAMETER            :: INITMR_BASIS = 2005

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
      ! OFFSETYEAR         : Number of year-ends since simulation start
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
      ! CLOCKMR            : Current base MR for the clock tracer
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
      ! NOX_O              : Monthly mean noontime O3P/O1D for NOx calcs
      ! NOX_J              : Monthly mean noontime J-rates for NOx calcs
      ! SO4_TOPPHOT        : Photolysis rate at the top of the chemgrid (1/s)
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
      ! STRAT2DCH4          : Use 2D estimates of stratospheric CH4?
      ! 
      !=================================================================

      ! Scalars
      INTEGER                           :: TRAC_IDX
      INTEGER                           :: N_EMIT
      INTEGER, PARAMETER                :: MAX_EMIT=50
      INTEGER                           :: SFC_MONTH
      INTEGER                           :: OFFSETYEAR
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
      REAL*8                            :: CLOCKMR

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
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: NOX_O
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: NOX_J
      REAL*8,DIMENSION(:,:),ALLOCATABLE     :: SO4_TOPPHOT
 
      ! Strings
      CHARACTER(LEN=255)                :: TRAC_FILE
      CHARACTER(LEN=2)                  :: MONTH_STR
      CHARACTER(LEN=255)                :: TRAC_SET 
      CHARACTER*20,DIMENSION(50)        :: TRAC_EMIT

      ! Logicals
      LOGICAL                           :: TRAC_ADD
      LOGICAL,PARAMETER                 :: STRAT2DCH4=.TRUE.

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
      USE LOGICAL_MOD,          ONLY : LSETH2SO4,LUCX,     LSETBRSTRAT
      USE LOGICAL_MOD,          ONLY : LSETNOYSTRAT
      USE LOGICAL_MOD,          ONLY : LBASICEMIS
      USE TRACERID_MOD
      USE TRACER_MOD,           ONLY : N_TRACERS, STT,TCVV
      USE DAO_MOD,              ONLY : AD
      USE TIME_MOD,             ONLY : GET_MONTH
      USE ERROR_MOD,            ONLY : ALLOC_ERR,ERROR_STOP
      USE GRID_MOD,             ONLY : GET_YMID
      USE TIME_MOD,             ONLY : GET_YEAR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_YEAR
      USE LOGICAL_MOD,          ONLY : LFUTURE, LFUTURECFC
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
      REAL*8             :: SCALEFAC
      INTEGER            :: FUTURE_LINE
      INTEGER            :: BASIS_LINE

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
      USE2DDATA = (USE2DDATA.or.LSETNOYSTRAT)
      USE2DDATA = (USE2DDATA.or.LSETN2O.or.LBASICEMIS)

      IF (USE2DDATA) THEN

         ! Calculate future scaling factor
         FUTURE_LINE = (12*(CFCYEAR+OFFSETYEAR-1959)) + GET_MONTH()
         BASIS_LINE = (12*(INITMR_BASIS-1959)) + GET_MONTH()

         ALLOCATE( AER_MR( JGLOB, AER_NLEVS ), STAT=AS )
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

            ! Scaling data only available for CH3Br
            SCALEFAC = 1.d0

            TRAC_IDX = IDTBr2
            TRAC_SET  = 'BR2'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
         
            TRAC_IDX = IDTBr
            TRAC_SET  = 'BR'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
         
            TRAC_IDX = IDTBrO
            TRAC_SET  = 'BRO'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
         
            TRAC_IDX = IDTHOBr
            TRAC_SET  = 'HOBR'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
         
            TRAC_IDX = IDTHBr
            TRAC_SET  = 'HBR'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
     
            TRAC_IDX = IDTBrNO2
            TRAC_SET  = 'BRNO2'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
     
            TRAC_IDX = IDTBrNO3
            TRAC_SET  = 'BRNO3'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
     
            TRAC_IDX = IDTCHBr3 
            TRAC_SET  = 'CHBR3'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
     
            TRAC_IDX = IDTCH2Br2
            TRAC_SET  = 'CH2BR2'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
     
            TRAC_IDX = IDTCH3Br
            TRAC_SET  = 'CH3BR'
            SCALEFAC = GET_MONTREAL('CH3Br',FUTURE_LINE)/
     &                    GET_MONTREAL('CH3Br',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
         ENDIF

         ! Stratospheric NOy
         IF (LSETNOYSTRAT) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|NOx      |HNO3     |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF

            ! No Montreal scaling
            SCALEFAC = 1.d0

            TRAC_IDX = IDTHNO3
            TRAC_SET = 'HNO3'
            CALL APPLY_2DTRAC(.FALSE.,.TRUE.,.FALSE.,SCALEFAC)

            TRAC_IDX = IDTNOx
            TRAC_SET = 'NOX'
            CALL APPLY_2DTRAC(.FALSE.,.TRUE.,.FALSE.,SCALEFAC)
         ENDIF

         ! Br-Cl (special case)
         IF (LSETCL.or.LSETBR.or.LSETBRSTRAT) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|BrCl     |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF

            ! No Montreal scaling
            SCALEFAC = 1.d0

            TRAC_IDX = IDTBrCl
            TRAC_SET  = 'BRCL'
            CALL APPLY_2DTRAC(.FALSE.,LSETBRSTRAT,.FALSE.,SCALEFAC)
     
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

            ! No Montreal scaling for non-long-lived-species
            SCALEFAC = 1d0

            TRAC_IDX = IDTCl
            TRAC_SET  = 'CL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_IDX = IDTClO
            TRAC_SET  = 'CLO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_IDX = IDTClNO2
            TRAC_SET  = 'CLNO2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_IDX = IDTClNO3
            TRAC_SET  = 'CLNO3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_IDX = IDTHOCl
            TRAC_SET  = 'HOCL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_IDX = IDTOClO
            TRAC_SET  = 'OCLO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_IDX = IDTClOO
            TRAC_SET  = 'CLOO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
     
            TRAC_IDX = IDTCl2
            TRAC_SET  = 'CL2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_IDX = IDTCl2O2
            TRAC_SET  = 'CL2O2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
     
            TRAC_IDX = IDTHCl
            TRAC_SET  = 'HCL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
     
            ! Handle chloroalkanes      
            TRAC_IDX = IDTCH3Cl
            TRAC_SET  = 'CH3CL'
            SCALEFAC = GET_MONTREAL('CH3Cl',FUTURE_LINE)/
     &                    GET_MONTREAL('CH3Cl',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
     
            TRAC_IDX = IDTCH3CCl3
            TRAC_SET  = 'CH3CCL3'
            SCALEFAC = GET_MONTREAL('CH3CCl3',FUTURE_LINE)/
     &                    GET_MONTREAL('CH3CCl3',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
     
            TRAC_IDX = IDTCCl4
            TRAC_SET  = 'CCL4'
            SCALEFAC = GET_MONTREAL('CCl4',FUTURE_LINE)/
     &                    GET_MONTREAL('CCl4',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
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
            SCALEFAC = GET_MONTREAL('H1211',FUTURE_LINE)/
     &                    GET_MONTREAL('H1211',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_IDX = IDTH1301
            TRAC_SET  = 'CBRF3'
            SCALEFAC = GET_MONTREAL('H1301',FUTURE_LINE)/
     &                    GET_MONTREAL('H1301',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_IDX = IDTH2402
            TRAC_SET  = 'C2BR2F4'

            SCALEFAC=GET_MONTREAL('C2BR2F4',CFCYEAR+OFFSETYEAR-1949)/
     &               GET_MONTREAL('C2BR2F4',INITMR_BASIS-1949)

            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
     
            ! CFC  11, 12, and [113/114/115]
            TRAC_IDX = IDTCFC11
            TRAC_SET  = 'CFCL3'
            SCALEFAC = GET_MONTREAL('CFC11',FUTURE_LINE)/
     &                    GET_MONTREAL('CFC11',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)

            TRAC_IDX = IDTCFC12
            TRAC_SET  = 'CF2CL2'
            SCALEFAC = GET_MONTREAL('CFC12',FUTURE_LINE)/
     &                    GET_MONTREAL('CFC12',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)

            TRAC_IDX = IDTCFCX
            TRAC_SET  = 'C2CL3F3'
            SCALEFAC = GET_MONTREAL('CFC113',FUTURE_LINE)/
     &                    GET_MONTREAL('CFC113',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)

            TRAC_SET  = 'C2CL2F4'
            SCALEFAC = GET_MONTREAL('CFC114',FUTURE_LINE)/
     &                    GET_MONTREAL('CFC114',BASIS_LINE)
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_SET  = 'C2CLF5'
            SCALEFAC = GET_MONTREAL('CFC115',FUTURE_LINE)/
     &                    GET_MONTREAL('CFC115',BASIS_LINE)
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.,SCALEFAC)
     
            ! HCFC 22, 123, 141b and 142b
            TRAC_IDX = IDTHCFC22
            TRAC_SET  = 'CHCLF2'
            SCALEFAC = GET_MONTREAL('HCFC22',FUTURE_LINE)/
     &                    GET_MONTREAL('HCFC22',BASIS_LINE)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)

            TRAC_IDX = IDTHCFCX
            TRAC_SET  = 'C2HF3CL2'
            ! Take scaling factor for HCFC22 (no data for 123)
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
            TRAC_SET  = 'C2H3FCL2'
            SCALEFAC = GET_MONTREAL('HCFC141b',FUTURE_LINE)/
     &                    GET_MONTREAL('HCFC141b',BASIS_LINE)
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.,SCALEFAC)

            TRAC_SET  = 'C2H3F2CL'
            SCALEFAC = GET_MONTREAL('HCFC142b',FUTURE_LINE)/
     &                    GET_MONTREAL('HCFC142b',BASIS_LINE)
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.,.FALSE.,SCALEFAC)
         ENDIF

         IF (LSETOCS) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|OCS      |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTOCS
            TRAC_SET  = 'OCS'

            ! No Montreal scaling
            SCALEFAC = 1.d0

            CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
         ENDIF

         IF (LSETH2SO4) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|H2SO4    |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            ! No Montreal scaling
            SCALEFAC = 1.d0

            TRAC_IDX = IDTSO4
            ! Actually H2SO4 aerosol
            TRAC_SET  = 'SO4'
            CALL APPLY_2DTRAC(.FALSE.,.TRUE.,.FALSE.,SCALEFAC)
         ENDIF
      ENDIF ! USE2DDATA

      IF (LSETN2O) THEN
         IF ( am_I_Root ) THEN
            WRITE(6,*) '|N2O      |         |         |'
            WRITE(6,*) '|---------|---------|---------|'
            IF ( LPRT ) CALL DEBUG_MSG( '### UCX: Forcing 2D N2O' )
         ENDIF

         TRAC_IDX = IDTN2O
         TRAC_SET  = 'N2O'
         SCALEFAC = GET_MONTREAL('N2O',FUTURE_LINE)/
     &                 GET_MONTREAL('N2O',BASIS_LINE)
         CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.FALSE.,SCALEFAC)
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
            ! No Montreal scaling
            SCALEFAC = 1.d0

            ! Code to set strat using AER 2D data
            TRAC_IDX = IDTCH4
            TRAC_SET  = 'CH4'
            CALL APPLY_2DTRAC(.FALSE.,.TRUE.,.FALSE.,SCALEFAC)
         ENDIF
      ENDIF
 
      END SUBROUTINE SET_INITIAL_MIXRATIOS
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ucx_nox
!
! !DESCRIPTION: Subroutine UCX\_NOX calculates NOx and N2O loss  
!  rates above the chemistry grid, based on estimates of j-rates from 
!  a 2D model and simple photochemical assumptiones.              
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE UCX_NOX
!
! !USES:
!
      USE ERROR_MOD,     ONLY : ERROR_STOP, DEBUG_MSG
      USE LOGICAL_MOD,   ONLY : LPRT
      USE TRACER_MOD,    ONLY : STT, TCVV
      USE TRACERID_MOD,  ONLY : IDTNOX, IDTN2O, IDTOX
      USE CHEMGRID_MOD,  ONLY : ITS_IN_THE_CHEMGRID
      USE DAO_MOD,       ONLY : AD, T, BXHEIGHT
      USE PRESSURE_MOD,  ONLY : GET_PCENTER
      USE TIME_MOD,      ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_MONTH
      USE TIME_MOD,      ONLY : ITS_A_LEAPYEAR, GET_HOUR, GET_LOCALTIME
      USE TIME_MOD,      ONLY : GET_MINUTE
      USE CMN_GCTM_MOD,  ONLY : PI, RE
      USE GRID_MOD,      ONLY : GET_YMID
      !-----SDE DEBUG
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE FILE_MOD,      ONLY : IOERROR
      USE TIME_MOD,      ONLY : GET_YEAR
      !-----SDE DEBUG
!
! !REVISION HISTORY: 
!  16 Jul 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! For indexing the O1D/O3P concentration array
      INTEGER, PARAMETER :: O3PIDX = 1
      INTEGER, PARAMETER :: O1DIDX = 2

      ! For indexing the J-rate coefficient array
      INTEGER, PARAMETER :: JNOIDX = 1
      INTEGER, PARAMETER :: JNO2IDX= 2
      INTEGER, PARAMETER :: JNO3IDX= 3
      INTEGER, PARAMETER :: JN2OIDX= 4

      ! Reaction rates and indices
      REAL*8,DIMENSION(12):: RRATE
      INTEGER, PARAMETER :: k_JNO = 6
      INTEGER, PARAMETER :: k_JNO2 = 4
      INTEGER, PARAMETER :: k_JNO3 = 5
      INTEGER, PARAMETER :: k_JN2O = 12

      ! Intermediate variables
      REAL*8 :: LOCALNOX, LOCALN2O, LOCALO3, LOCALO1D, LOCALO3P
      REAL*8 :: FRACNO, FRACNO2, FRACNO3, FRACN
      REAL*8 :: NOXRATE, N2ORATE, KGNOX, KGN2O
      REAL*8 :: NO_ALPHA, NO_BETA, NO_GAMMA, NO_EPSILON
      REAL*8 :: NO_QA, NO_QX, NO_QC, DNOX, DN2O
      REAL*8 :: MESONOX_DELTA
      REAL*8 :: MESON2O_DELTA

      ! Local air number density (molec/cm3) and box mass (kg)
      REAL*8 :: NDAIR
      REAL*8 :: XAIR

      ! Local temperature (K) and inverted T (1/K)
      REAL*8 :: T3K, TINV

      ! Chemistry timestep in seconds
      REAL*8 :: DTCHEM

      ! Timestep in hours
      REAL*8 :: DTCHEM_HR

      ! Timing information
      INTEGER       :: NDAY
      INTEGER, SAVE :: LASTMONTH=0
      !-----SDE DEBUG
      INTEGER, SAVE :: LASTDAY=0
      !-----SDE DEBUG

      ! Grid indexing
      INTEGER :: I,J,L
      REAL*8  :: MIDLAT, ZMID, ZDEL, ZBASE
 
      ! Local daylight fraction
      REAL*8             :: DAYFRAC
      REAL*8             :: SOL_ME, SOL_FE, SOL_LS, SOL_DEC
      REAL*8             :: SOL_ECC2, SOL_B, SOL_C, SOL_TEST
      REAL*8             :: TSTART, TEND, TMIN, TMAX, GMTSTART
      REAL*8, PARAMETER  :: SOL_ECC=0.016
      INTEGER            :: NDYEAR, AS
      LOGICAL            :: CYCLEBOX
      CHARACTER(LEN=255) :: DBGMSG

      !-----SDE DEBUG
      INTEGER            :: IU_FILE, IOS
      CHARACTER(LEN=255) :: STRATNOXFILE, NOXSTRING, STRATNOXSTAT
      CHARACTER(LEN=255) :: STRATNOXACCESS
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: NOXARR
      LOGICAL, PARAMETER :: WRITENOX = .FALSE.
      !-----SDE DEBUG

      !=================================================================
      ! UCX_NOX begins here!
      !=================================================================

      ! Retrieve monthly mean data if necessary
      IF (LASTMONTH.ne.GET_MONTH()) THEN
         LASTMONTH = GET_MONTH()
         CALL GET_NOXCOEFF( LASTMONTH )
      ENDIF

      ! Get chemistry step length in seconds
      DTCHEM = GET_TS_CHEM()*60.d0
      DTCHEM_HR = DTCHEM/3600.d0

      IF (ITS_A_LEAPYEAR()) THEN
         NDYEAR = 366
      ELSE
         NDYEAR = 365
      ENDIF

      NDAY = GET_DAY_OF_YEAR()

      ! Calculate timestep start and end (GMT)
      GMTSTART = REAL(GET_HOUR()) + (REAL(GET_MINUTE())/60.d0)
      GMTSTART = GMTSTART - (DTCHEM_HR/2.d0)

      ! Eccentricity^2
      SOL_ECC2 = SOL_ECC*SOL_ECC

      ! Calculate solar declination
      SOL_ME = 2.d0*PI*(REAL(NDAY)/REAL(NDYEAR)-0.009589)
      SOL_FE = SOL_ME + (2.*SOL_ECC*SIN(SOL_ME)) + 
     &         (1.25d0*SOL_ECC2*SIN(2.d0*SOL_ME))
      SOL_LS = SOL_FE - (77.240d0*PI/180.d0)
      SOL_DEC = ASIN(SIN(23.44*PI/180.d0)*SIN(SOL_LS))

      ! Reset NOx/N2O mass counters
      MESONOX_DELTA = 0d0
      MESON2O_DELTA = 0d0

      !-----SDE DEBUG
      IF (WRITENOX) THEN
         ALLOCATE( NOXARR( IIPAR,JJPAR,LLPAR,8 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXARR' )
         NOXARR=0d0
      ENDIF
      !-----SDE DEBUG

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, DAYFRAC, RRATE               )
!$OMP+PRIVATE( T3K, TINV, XAIR,  KGNOX, KGN2O        )
!$OMP+PRIVATE( FRACNO, FRACNO2, FRACNO3, LOCALNOX )
!$OMP+PRIVATE( LOCALN2O, LOCALO1D, LOCALO3P, LOCALO3 )
!$OMP+PRIVATE( NOXRATE, N2ORATE, FRACN, NO_ALPHA    )
!$OMP+PRIVATE( NO_BETA, NO_GAMMA, NDAIR, MIDLAT      )
!$OMP+PRIVATE( SOL_B, SOL_C, SOL_TEST, NO_EPSILON    )
!$OMP+PRIVATE( ZMID, ZBASE, ZDEL, CYCLEBOX           )
!$OMP+PRIVATE( TSTART, TEND, TMIN, TMAX, DNOX, DN2O  )
!$OMP+PRIVATE( NO_QA, NO_QX, NO_QC                   )
      DO J=1,JJPAR
      DO I=1,IIPAR
      ! Grid box longitude/latitude [degrees]
      MIDLAT = PI*GET_YMID( I, J, 1 )/180.d0
      SOL_B = TAN(SOL_DEC)*TAN(MIDLAT)
      SOL_C = COS(SOL_DEC)*COS(MIDLAT)
      !SOLFUNC  = 1.d0-(0.034d0*cos(dble(NDAY-172)*2.d0*PI/365.d0))
      TSTART = GET_LOCALTIME(I,J,1,GMTSTART)
      TEND = TSTART + DTCHEM_HR
      !DAYFRAC = SOLFUNC
      ZBASE = 0d0
      DO L=1,LLPAR
         ZDEL = BXHEIGHT(I,J,L)
         ZMID = ZBASE + (ZDEL/2.d0)
         ZBASE = ZBASE + ZDEL
         CYCLEBOX = .FALSE.
         ! Are the pre-calculated J-rates zero (ie in darkness)?
         CYCLEBOX = (ANY(NOX_J(I,J,L,:).lt.TINY(1d0)))
         IF (ITS_IN_THE_CHEMGRID(I,J,L).or.CYCLEBOX) THEN
!            !-----SDE DEBUG
!            WRITE(IU_FILE,222) I,J,L,0.d0,0.d0,0.d0,0.d0,0.d0,
!     &              0.d0,0.d0,0.d0,GMTSTART,DTCHEM_HR
!            ! Change this to CYCLE
!            CYCLEBOX = .TRUE.
!            !-----SDE DEBUG
            CYCLE
         ENDIF
         ! Calculate daylight fraction
         SOL_TEST = (-SOL_B) - SQRT(ZMID*(ZMID+2.d0*RE)/
     &     ((ZMID+RE)*(ZMID+RE)))/SOL_C
         IF (SOL_TEST.ge.0.9998d0) THEN
            DAYFRAC = 0.d0
            CYCLEBOX=.TRUE.
            TMIN = 0d0
            TMAX = 0d0
         ELSEIF (SOL_TEST.le.-0.998d0) THEN
            DAYFRAC = 1.d0
            TMIN = TSTART
            TMAX = TEND
            CYCLEBOX=.FALSE.
         ELSE
            DAYFRAC = ACOS(SOL_TEST)/PI
            ! Determine how much of the current timestep is within daylight
            TMIN = 12.d0*(1.d0-DAYFRAC)
            TMAX = 12.d0*(1.d0+DAYFRAC)
            IF ((TMIN.ge.TEND).or.(TMAX.le.TSTART)) THEN
               ! In darkness
               DAYFRAC = 0.d0
               CYCLEBOX = .TRUE.
            ELSEIF ((TSTART.gt.TMIN).and.(TEND.lt.TMAX)) THEN
               ! Full daylight
               DAYFRAC = 1.d0
            ELSE
               ! Partial daylight
               TMIN = MAX(TSTART,TMIN)
               TMAX = MIN(TEND,TMAX)
               DAYFRAC = (TMAX-TMIN)/DTCHEM_HR
            ENDIF
         ENDIF
!         !-----SDE DEBUG
!         WRITE(6,'(3(I0.3,1x),5(F10.4,1x),L1,1x,3(E10.4,1x))') I, J, L,
!     &      TSTART,TEND,TMIN,TMAX,DAYFRAC,CYCLEBOX,SOL_B,SOL_C,SOL_TEST
!         !-----SDE DEBUG
         IF (CYCLEBOX) THEN
!            !-----SDE DEBUG
!            WRITE(IU_FILE,222) I,J,L,0.d0,0.d0,0.d0,0.d0,0.d0,
!     &              0.d0,0.d0,1.d0,GMTSTART,DTCHEM_HR
!            !-----SDE DEBUG
            CYCLE
         ENDIF
         ! Reset reaction rates
         RRATE = 0d0

         ! Retrieve air mass (kg) and local temperature (K)
         T3K = T(I,J,L)
         TINV = 1.d0/T3K

         ! Calculate air number density (molec/cm3)
         NDAIR = GET_PCENTER(I,J,L)*1.d-4*AVOGADRO/(R_UNIV*T3K)
         XAIR = NDAIR/AD(I,J,L)

         ! Calculate local NOx, N2O and O3 number density
         KGNOX = STT(I,J,L,IDTNOX)
         KGN2O = STT(I,J,L,IDTN2O)
         LOCALNOX = KGNOX*TCVV(IDTNOX)*XAIR
         LOCALN2O = KGN2O*TCVV(IDTN2O)*XAIR
         LOCALO3  = STT(I,J,L,IDTOX)*TCVV(IDTOX)*XAIR

         ! Retrieve local O3P/O1D mixing ratios and relevant
         ! j-rates from interpolated 2D arrays
         LOCALO3P = NOX_O(I,J,L,1)*NDAIR*DAYFRAC
         LOCALO1D = NOX_O(I,J,L,2)*NDAIR*DAYFRAC

         ! Calculate reaction rates
         ! 1:  NO2 + O -> NO + O2
         RRATE(1)  = 5.1d-12*dexp(210.d0*TINV)
         ! 2:  NO + O3 -> NO2 + O2
         RRATE(2)  = 3.0d-12*dexp(-1500.d0*TINV)
         ! 3:  NO2 + O3 -> NO3 + O2
         RRATE(3)  = 1.2d-13*dexp(-2450.d0*TINV)
         ! 4:  NO2 + hv -> NO + O1D
         RRATE(k_JNO2) = NOX_J(I,J,L,JNO2IDX)*DAYFRAC
         ! 5:  NO3 + hv -> NO2 + O
         RRATE(k_JNO3) = NOX_J(I,J,L,JNO3IDX)*DAYFRAC
         ! 6:  NO + hv -> N + O
         RRATE(k_JNO ) = NOX_J(I,J,L,JNOIDX)*DAYFRAC
         ! 7:  N + NO2 -> N2O + O
         RRATE(7)  = 5.8d-12*dexp(220.d0*TINV)
         ! 8:  N + NO -> N2 + O
         RRATE(8)  = 2.1d-11*dexp(100.d0*TINV)
         ! 9:  N + O2 -> NO + O
         RRATE(9)  = 1.5d-11*dexp(-3600.d0*TINV)
         ! 10:  N2O + O1D -> N2 + O2
         RRATE(10) = 4.63d-11*dexp(20.d0*TINV)
         ! 11:  N2O + O1D -> 2NO
         RRATE(11) = 7.25d-11*dexp(20.d0*TINV)
         ! 12:  N2O + hv -> N2 + O1D
         RRATE(k_JN2O) = NOX_J(I,J,L,JN2OIDX)*DAYFRAC

         ! Partition NOx into N, NO, NO2 and NO3 based on PSSA
         ! Two cases: Daytime/nighttime
         NO_ALPHA =          RRATE(k_JNO)/
     &                  (RRATE(9)*0.21d0*NDAIR)
         NO_BETA  = (RRATE(k_JNO2)+(RRATE(1)*LOCALO3P))/
     &                     (RRATE(2)*LOCALO3)
         NO_GAMMA =        (RRATE(3)*LOCALO3)/
     &                         RRATE(k_JNO3)

         ! Calculate the partition fractions
         FRACNO2 = 1.d0/(1.d0+NO_GAMMA+(NO_BETA*(1.d0+NO_ALPHA)))
         FRACNO3 = NO_GAMMA * FRACNO2
         FRACNO  = NO_BETA  * FRACNO2
         FRACN   = NO_ALPHA * FRACNO

         ! Estimate net production rates for NOx and N2O
         NOXRATE = 2.d0*((RRATE(11)*LOCALN2O*LOCALO1D) -
     &    (((RRATE(7)*FRACN*FRACNO2)+(RRATE(8)*FRACN*FRACNO))
     &                                   *LOCALNOX*LOCALNOX))
         N2ORATE = (RRATE(7)*FRACN*FRACNO2*LOCALNOX*LOCALNOX) - 
     &   ((RRATE(k_JN2O)+((RRATE(10)+RRATE(11))*LOCALO1D))*LOCALN2O)
        
         ! Calculate fractional loss and apply to tracer
         ! Explicit Euler method (fast)
         DNOX = NOXRATE * DTCHEM / (TCVV(IDTNOX)*XAIR)
         IF ((DNOX*-1d0).gt.KGNOX) THEN
            DNOX = -1d0*KGNOX
            KGNOX = 0d0
         ELSE
            KGNOX = KGNOX + DNOX
         ENDIF
         STT(I,J,L,IDTNOX) = KGNOX

         MESONOX_DELTA = MESONOX_DELTA + DNOX

         DN2O = N2ORATE * DTCHEM / (TCVV(IDTN2O)*XAIR)
         IF ((DN2O*-1d0).gt.KGN2O) THEN
            DN2O = -1d0*KGN2O
            KGN2O = 0d0
         ELSE
            KGN2O = KGN2O + DN2O
         ENDIF
         STT(I,J,L,IDTN2O) = KGN2O

         MESON2O_DELTA = MESON2O_DELTA + DN2O
         
         !-----SDE DEBUG
         IF (WRITENOX) THEN
            NOXARR(I,J,L,1) = N2ORATE/NDAIR  ! (v/v)/s
            NOXARR(I,J,L,2) = NOXRATE/NDAIR  ! (v/v)/s
            NOXARR(I,J,L,3) = LOCALN2O/NDAIR ! v/v
            NOXARR(I,J,L,4) = LOCALNOX/NDAIR ! v/v
            NOXARR(I,J,L,5) = FRACN
            NOXARR(I,J,L,6) = FRACNO
            NOXARR(I,J,L,7) = FRACNO2
            NOXARR(I,J,L,8) = FRACNO3
         ENDIF
         !-----SDE DEBUG

!         ! Implicit Euler method (stable)
!         NO_ALPHA  = 2.d0*RRATE(11)*LOCALO1D
!         NO_GAMMA  = RRATE(7)*FRACN*FRACNO2   
!         NO_BETA   = 2.d0*(NO_GAMMA + (RRATE(8)*FRACN*FRACNO))
!         NO_EPSILON= RRATE(k_JN2O) + (LOCALO1D*(RRATE(10)+RRATE(11)))
!         ! Temporary storage variable
!         NO_QX     = DTCHEM*NO_ALPHA/(1.d0+(NO_EPSILON*DTCHEM))
!         ! Calculate quadratic coefficients (B = 1)
!         NO_QA     = DTCHEM*(NO_BETA - (NO_QX*NO_GAMMA))
!         NO_QC     = -1.d0*(LOCALNOX + (LOCALN2O*NO_QX))
!         ! Calculate new molecular concentrations
!         LOCALNOX  = (-1.d0 + SQRT(1.d0 - (4.d0*NO_QA*NO_QC)))/
!     &                              (2.d0*NO_QA)
!         LOCALN2O  = (LOCALN2O+(DTCHEM*LOCALNOX*LOCALNOX*NO_GAMMA))/
!     &               (1.d0 + (NO_EPSILON*DTCHEM))
!         KGNOX     = LOCALNOX/(TCVV(IDTNOX)*XAIR)
!         KGN2O     = LOCALN2O/(TCVV(IDTN2O)*XAIR)
!         STT(I,J,L,IDTNOX) = KGNOX
!         STT(I,J,L,IDTN2O) = KGN2O
      ENDDO ! J
      ENDDO ! I
      ENDDO ! L
!$OMP END PARALLEL DO

      IF (LPRT) THEN
         ! Print mean NOx tendency in mesosphere
         DBGMSG = ' ### UCX_NOX: Mesospheric NOx processed'
         CALL DEBUG_MSG(TRIM(DBGMSG))
         WRITE(DBGMSG,'(a,2(1x,F10.5),a)') ' ### Time (start/step):',
     &       GMTSTART, DTCHEM_HR, ' hours'
         CALL DEBUG_MSG(TRIM(DBGMSG))
         WRITE(DBGMSG,'(a,1x,E10.4,1x,a)') ' ### NOx direct delta: ',
     &       MESONOX_DELTA, ' kg'
         CALL DEBUG_MSG(TRIM(DBGMSG))
         WRITE(DBGMSG,'(a,1x,E10.4,1x,a)') ' ### N2O direct delta: ',
     &       MESON2O_DELTA, ' kg'
         CALL DEBUG_MSG(TRIM(DBGMSG))
      ENDIF

      !-----SDE DEBUG
      ! Write full set of diagnostics to file
      IF (WRITENOX) THEN
         IF (LASTDAY.ne.NDAY) THEN
            STRATNOXSTAT = 'REPLACE'
            STRATNOXACCESS='SEQUENTIAL'
         ELSE
            STRATNOXSTAT = 'OLD'
            STRATNOXACCESS='APPEND'
         ENDIF
         LASTDAY = NDAY

         ! Get a free LUN
         IU_FILE = findFreeLUN()
         WRITE(STRATNOXFILE,'(a,I0.4,I0.3)') 'STRATNOX_',GET_YEAR(),NDAY
         IOS = 1
         OPEN( IU_FILE, FILE=TRIM( STRATNOXFILE ),
     &          ACCESS=TRIM( STRATNOXACCESS ),
     &          STATUS=TRIM(STRATNOXSTAT),IOSTAT=IOS) 
         IF ( IOS /= 0 ) THEN
           WRITE(6,*) 'UCX_MOD: Could not read ',TRIM(STRATNOXFILE)
           CALL IOERROR( IOS, IU_FILE,'UCX_MOD:UCX_NOX')
         ENDIF
         DO I=1,IIPAR
         DO J=1,JJPAR
         DO L=LLCHEM,LLPAR
         IF (NOXARR(I,J,L,7).gt.0d0) THEN
            WRITE(IU_FILE, 222) I,J,L,NOXARR(I,J,L,1),NOXARR(I,J,L,2),
     &                                NOXARR(I,J,L,3),NOXARR(I,J,L,4),
     &                                NOXARR(I,J,L,5),NOXARR(I,J,L,6),
     &                                NOXARR(I,J,L,7),NOXARR(I,J,L,8),
     &                                GMTSTART,DTCHEM_HR
         ENDIF
         ENDDO
         ENDDO
         ENDDO
222      FORMAT (3(I0.3,1x), 8(E14.5,1x), F10.5, 1x, F10.5)
         CLOSE(IU_FILE)
         DEALLOCATE(NOXARR)
      ENDIF
      !-----SDE DEBUG

      END SUBROUTINE UCX_NOX
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_noxcoeff
!
! !DESCRIPTION: Subroutine GET\_NOXCOEFF reads in O1D and O3P mixing  
!  ratios along with NO, NO2, NO3 and N2O J-rates from 2D data, 
!  interpolating onto the 3D grid and storing in NOX_O and NOX_J.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_NOXCOEFF( TARG_MONTH )
!
! !USES:
!
      USE FILE_MOD,      ONLY : IOERROR
      USE ERROR_MOD,     ONLY : ERROR_STOP, ALLOC_ERR
      USE PRESSURE_MOD,  ONLY : GET_PCENTER
      USE TRACER_MOD,    ONLY : STT, TCVV
      USE CHEMGRID_MOD,  ONLY : ITS_IN_THE_CHEMGRID
      USE DAO_MOD,       ONLY : AD
      !-----SDE DEBUG
      USE CHEMGRID_MOD,  ONLY : GET_CHEMGRID_LEVEL
      USE FAST_JX_MOD,   ONLY : FJXFUNC
      !-----SDE DEBUG
!
! !INPUT PARAMETERS: 
!
      INTEGER,INTENT(IN) :: TARG_MONTH
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
      INTEGER :: I, J, L, ITRAC
      REAL*8  :: PCENTER
      REAL*8  :: CURRVAL
      INTEGER :: VERTCOUNT, AS
      LOGICAL :: FOUNDLEV,EXTRAP
      INTEGER :: IU_FILE
      LOGICAL :: ISRATE
      CHARACTER(LEN=255) :: NOX_FILE
      CHARACTER(LEN=255) :: TARG_TRAC
      CHARACTER(LEN=255) :: DBGMSG
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE :: NOXDATA2D
      INTEGER :: LSTART
      !-----SDE DEBUG
      REAL*8,DIMENSION(4) :: JTEMP
      !-----SDE DEBUG

      !=================================================================
      ! GET_NOXCOEFF begins here!
      !=================================================================

      ! Allocate and zero 2D data array
      ALLOCATE( NOXDATA2D( JJPAR, 51, 6 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXDATA2D' )
      NOXDATA2D = 0d0

      ! Clear interpolated arrays
      NOX_O = 0d0
      NOX_J = 0d0

      !-----SDE DEBUG
      ! j-value data not yet available
      !DO ITRAC = 1,2
      !-----SDE DEBUG
      DO ITRAC = 1,6
         SELECT CASE (ITRAC)
            CASE ( 1 )
               TARG_TRAC = 'O'
            CASE ( 2 )
               TARG_TRAC = 'O1D'
            CASE ( 3 )
               TARG_TRAC = 'JNO'
            CASE ( 4 )
               TARG_TRAC = 'JNO2'
            CASE ( 5 )
               TARG_TRAC = 'JNO3'
            CASE ( 6 )
               TARG_TRAC = 'JN2O'
         END SELECT
         WRITE(NOX_FILE,'(a,a,a,I0.2,a)') TRIM(NOON_FILE_ROOT),
     &    TRIM(TARG_TRAC), '_', TARG_MONTH, '.dat'
!         NOX_FILE = TRIM(NOON_FILE_ROOT)//TRIM(TARG_TRAC)
!         NOX_FILE = TRIM(NOX_FILE)//'_'//TRIM(TARG_MONTH)//'.dat'

         ! Get a free LUN
         IU_FILE = findFreeLUN()

         IOS = 1
         OPEN( IU_FILE, FILE=TRIM( NOX_FILE ),STATUS='OLD',IOSTAT=IOS) 
         IF ( IOS /= 0 ) THEN
            WRITE(6,*) 'UCX_MOD: Could not read ', TRIM(NOX_FILE)
            CALL IOERROR( IOS, IU_FILE,'UCX_MOD:GET_NOXCOEFF')
         ENDIF
  
        IF (LPRT) THEN
            WRITE(DBGMSG,'(a,a)') ' ### UCX: Reading ', 
     &        TRIM( NOX_FILE )
            CALL DEBUG_MSG( DBGMSG )
         ENDIF
 
      ! Read in data
         DO ILEV = 1,51
            READ(IU_FILE, 110, IOSTAT=IOS ) NOXDATA2D(:,ILEV,ITRAC)
            IF ( IOS /= 0 ) THEN
               WRITE(6,'(a,a,I4,a,1x,a)') 'UCX_MOD: Error reading '
     & , 'line ', ILEV, ' in file ', TRIM( NOX_FILE )
               CALL IOERROR( IOS, IU_FILE,'UCX_MOD:GET_NOXCOEFF')
            ENDIF
         ENDDO

#if defined(GRID2x25)
110   FORMAT(91E10.3)
#elif defined(GRID4x5)
110   FORMAT(46E10.3)
#endif

         CLOSE(IU_FILE)
      ENDDO 


      ! Scan through target array, element by element
      DO J = 1,JJPAR
      DO I = 1,IIPAR
         ! Vertcount is the layer count for the input, where layer 1
         ! is at the *top* of the atmosphere
         VERTCOUNT = 51
         EXTRAP = .TRUE.
         LSTART = GET_CHEMGRID_LEVEL(I,J)

         DO L = LSTART, LLPAR
            IF (ITS_IN_THE_CHEMGRID(I,J,L)) CYCLE
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

            DO ITRAC=1,6
               ISRATE = (ITRAC.gt.2)
               ! Interpolate data
               IF (EXTRAP) THEN
                  ! Just take outside value if at edges
                  CURRVAL = NOXDATA2D(J,VERTCOUNT,ITRAC)
               ELSE
                  ! Interpolate by pressure
                  CURRVAL = (AER_PLEVS(VERTCOUNT+1)-PCENTER)
                  CURRVAL = CURRVAL/(AER_PLEVS(VERTCOUNT+1)-
     &                               AER_PLEVS(VERTCOUNT))
                  CURRVAL = CURRVAL * 
     &               (NOXDATA2D(J,VERTCOUNT+1,ITRAC)-
     &                NOXDATA2D(J,VERTCOUNT,ITRAC))
                  CURRVAL = NOXDATA2D(J,VERTCOUNT+1,ITRAC) - CURRVAL
               ENDIF
               IF (.not.ISRATE) THEN
                  ! Reading in mixing ratios (v/v)
                  NOX_O(I,J,L,ITRAC) = CURRVAL
               ELSE
                  ! J-rate (no conversion necessary)
                  NOX_J(I,J,L,ITRAC-2) = CURRVAL
               ENDIF
            ENDDO ! ITRAC
         ENDDO ! L
      ENDDO ! I
      ENDDO ! J

      DEALLOCATE( NOXDATA2D )
 
      END SUBROUTINE GET_NOXCOEFF
 
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
      SUBROUTINE APPLY_2DTRAC( TRAC_ADD, STRAT_ONLY, SFC_ONLY, SCALEFAC)
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
      REAL*8, INTENT(IN)               :: SCALEFAC     ! Future scaling factor
      
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
      CHARACTER(LEN=255) :: DBGMSG

      !=================================================================
      ! APPLY_2DTRAC begins here!
      !=================================================================

      TRAC_FILE = TRIM(AVG_FILE_ROOT)//TRIM(TRAC_SET )
      TRAC_FILE = TRIM(TRAC_FILE)//'_'//TRIM(MONTH_STR)//'.dat'

      ! Error check
      IF ((STRAT_ONLY).and.(SFC_ONLY)) THEN
         CALL ERROR_STOP('Cannot have both SFC_ONLY and STRAT_ONLY!',
     &                           'APPLY_2DTRAC(UCX_mod.F)')
      ENDIF

      ! Get a free LUN
      IU_FILE = findFreeLUN()

      IOS = 1
      OPEN( IU_FILE, FILE=TRIM( TRAC_FILE ),STATUS='OLD',IOSTAT=IOS) 
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) 'UCX_MOD: Could not read ', TRIM( TRAC_FILE )
         CALL IOERROR( IOS, IU_FILE,'UCX_MOD:APPLY_2DTRAC')
      ENDIF
  
      IF (LPRT) THEN
         WRITE(DBGMSG,'(a,a)') ' ### UCX: Reading ', 
     &        TRIM( TRAC_FILE )
         CALL DEBUG_MSG( DBGMSG )
      ENDIF
 
      ! Read in mixing ratios
      DO ILEV = 1,51
         READ(IU_FILE, 110, IOSTAT=IOS ) AER_MR(:,ILEV)
         IF ( IOS /= 0 ) THEN
            WRITE(6,'(a,I4,a,1x,a)') 'UCX_MOD: Error reading line '
     & , ILEV, ' in file ', TRIM( TRAC_FILE )
            CALL IOERROR( IOS, IU_FILE,'UCX_MOD:APPLY_2DTRAC')
         ENDIF
      ENDDO

#if defined(GRID2x25)
110   FORMAT(91E10.3)
#elif defined(GRID4x5)
110   FORMAT(46E10.3)
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
                  TRACMASS = SCALEFAC*TRACMR*AD(I,J,L)/TCVV(TRAC_IDX)
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
      USE LOGICAL_MOD,  ONLY : LFUTURECFC
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
      CHARACTER(LEN=255) :: OUTMSG

      !=================================================================
      ! EMISS_BASIC begins here!
      !=================================================================

       ! If not correct hour block and day read in new fields
       IF (SFC_MONTH/=GET_MONTH()) THEN
          IF (GET_MONTH().lt.SFC_MONTH) THEN
             OFFSETYEAR = OFFSETYEAR + 1
          ENDIF
          SFC_MONTH = GET_MONTH()
          IF (LPRT) THEN
             WRITE(OUTMSG,'(a,a,1x,I4)') ' ### UCX: ',
     &          'Using Montreal reference year', (CFCYEAR+OFFSETYEAR)
             CALL DEBUG_MSG( OUTMSG )
          ENDIF
          CALL READ_SFC
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

      IF ( LPRT ) CALL DEBUG_MSG( '### UCX: Emissions complete' )

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
      USE TRACERID_MOD,   ONLY : IDTHCFCX
      USE LOGICAL_MOD,    ONLY : LFUTURECFC
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
      REAL*8,PARAMETER :: H2_MR = 5.00d-7  ! H2 is fixed to 500 ppbv
      INTEGER          :: N, TARG_LINE
      CHARACTER(LEN=20):: LOCAL_NAME

      !=================================================================
      ! READ_SFC begins here!
      !=================================================================

      ! Re-zero GRID_EMIT
      GRID_EMIT = 0d0

      ! Determine which line needs to be read in
      TARG_LINE = (12*(CFCYEAR+OFFSETYEAR-1959)) + SFC_MONTH

      ! Loop over tracer IDs pegged for emissions
      DO N=1,N_EMIT
         LOCAL_NAME = TRIM(TRAC_EMIT(N))
         SELECT CASE (LOCAL_NAME)
            CASE ('CFCX')
               ! Add 3 CFCs
               LOCAL_NAME = TRIM('CFC113')
               CALL SET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('CFC114')
               CALL SET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('CFC115')
               CALL SET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
            CASE ('HCFCX')
               ! Add 2 HCFCs
               LOCAL_NAME = TRIM('HCFC141b')
               CALL SET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('HCFC142b')
               CALL SET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               ! Also need HCFC123 - no Montreal data
               ! Use surface-level MR from 2D data instead
               !TRAC_IDX = IDTHCFCX
               !TRAC_SET  = 'C2HF3CL2'
               !CALL APPLY_2DTRAC(.FALSE.,.FALSE.,.TRUE.)
               ! Ignore HCFC123 for the moment
            CASE ('H2402')
               ! Limited data available for H2402
               LOCAL_NAME = TRIM('C2BR2F4')
               CALL SET_MONTREAL(LOCAL_NAME,
     &                      CFCYEAR+OFFSETYEAR-1949,N)
            CASE ('OCS')
               ! Simple boundary condition
               GRID_EMIT(:,N) = OCS_MR
            CASE ('H2')
               ! Simple boundary condition
               GRID_EMIT(:,N) = H2_MR
            CASE DEFAULT
               ! Standard 1:1
               CALL SET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
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
! !IROUTINE: set_montreal
!
! !DESCRIPTION: Subroutine SET\_MONTREAL calls GET\_MONTREAL for the surface 
!  mixing ratio for a specific year, and applies it to GRID_EMIT.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_MONTREAL (TRAC_NAME,TARG_LINE,N)
!
! !USES:
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=20),INTENT(IN)   :: TRAC_NAME
      INTEGER,INTENT(IN)             :: TARG_LINE
      INTEGER,INTENT(IN)             :: N
!
! !REVISION HISTORY: 
!  10 Oct 2013 - S. D. Eastham - Wrapper built from edges of original GET
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! SET_MONTREAL begins here!
      !=================================================================

      GRID_EMIT(:,N) = GRID_EMIT(:,N)+GET_MONTREAL(TRAC_NAME,TARG_LINE)

      ! Return to calling routine
      END SUBROUTINE SET_MONTREAL
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_montreal
!
! !DESCRIPTION: Function GET\_MONTREAL gets mean surface mixing ratio of a 
!  given species based on forecasts from WHO estimates. Values are
!  interpolated between two monthly averages.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MONTREAL (TRAC_NAME,TARG_LINE) RESULT (MEAN_MIXRATIO)
!
! !USES:
!
      USE FILE_MOD,  ONLY : IOERROR
      USE ERROR_MOD, ONLY : ERROR_STOP, ALLOC_ERR
!
! !INPUT PARAMETERS:
!
      REAL*8                         :: MEAN_MIXRATIO
      CHARACTER*(*),INTENT(IN)       :: TRAC_NAME
      INTEGER,INTENT(IN)             :: TARG_LINE
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!  10 Oct 2013 - S. D. Eastham - Split into two routines (GET and SET)
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
      INTEGER                        :: IU_FILE, IOS, AS

      !=================================================================
      ! SET_MONTREAL begins here!
      !=================================================================

      FILENAME = TRIM(MONTREAL_FILE_ROOT) // TRIM(TRAC_NAME) // 
     &             '_boundary'

      ! Get a free LUN
      IU_FILE = findFreeLUN()

      IOS=1

      OPEN( IU_FILE, FILE=TRIM( FILENAME ), IOSTAT=IOS )
      IF ( IOS /= 0 ) THEN 
         WRITE(6,*) 'UCX: FAILED TO FIND ', TRIM(FILENAME)
          CALL IOERROR( IOS, IU_FILE, 'UCX_MOD:SET_MONTREAL')
      ENDIF

      ! First line not important
      READ(IU_FILE,*,IOSTAT=IOS)
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) 'UCX: READ ERROR IN ', TRIM(FILENAME)
          CALL IOERROR( IOS, IU_FILE, 'UCX_MOD:SET_MONTREAL')
      ENDIF

      ! Second line contains scaling factor
      READ(IU_FILE,101,IOSTAT=IOS) SCALEFAC

101   FORMAT(E7.1)

      ! Third and fourth line also to be ignored
      DO ILINE=1,2
         READ(IU_FILE,*,IOSTAT=IOS)
         IF ( IOS /= 0 ) THEN
            WRITE(6,*) 'UCX: READ ERROR IN ', TRIM(FILENAME)
             CALL IOERROR( IOS, IU_FILE, 'UCX_MOD:SET_MONTREAL')
         ENDIF
      ENDDO

      ! Skip to line specified
      DO ILINE=1,(TARG_LINE-1)
         READ(IU_FILE,*,IOSTAT=IOS) ! Skip line
         IF ( IOS /= 0 ) THEN 
            WRITE(6,*) 'UCX: READ ERROR IN ', TRIM(FILENAME)
             CALL IOERROR( IOS, IU_FILE, 'UCX_MOD:SET_MONTREAL')
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
         WRITE(6,*) 'UCX: READ ERROR IN ', TRIM(FILENAME)
          CALL IOERROR( IOS, IU_FILE, 'UCX_MOD:SET_MONTREAL')
      ENDIF

      CLOSE(IU_FILE)

      MEAN_MIXRATIO = TEMPMR * SCALEFAC

102   FORMAT(F7.2,F11.2)
      ! C2Br2F4 only
103   FORMAT(F7.2,F8.2)

      ! Return to calling program
      END FUNCTION GET_MONTREAL
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
      USE TRACERID_MOD,  ONLY : IDTNIT,IDTBCPI,IDTH2O,IDTHNO3
      USE TIME_MOD,      ONLY : GET_ELAPSED_SEC, GET_TS_CHEM
      USE COMODE_MOD,    ONLY : JLOP,WERADIUS,ERADIUS
      USE LOGICAL_MOD,   ONLY : LUCX,LGRAVSTRAT
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
      REAL*8                 :: SEDPSC, SEDNAT, SEDICE, SEDQUANT
      REAL*8                 :: SEDH2O, SEDNO3
      REAL*8                 :: BELOWGRAD, ABOVEGRAD
      INTEGER                :: LOCALPROFILE, NUMSEDSTEPS, STARTPT, ISED

      ! SDE 2013-07-07: Debug variables
      REAL*8                 :: D_XH2O, D_PPA,  D_TK,   D_XHNO3,D_XNAT
      REAL*8                 :: D_ICE1, D_ICE0, D_H2O1, D_H2O0
      REAL*8                 :: D_NAT1, D_NAT0, D_SEDQ, D_PSED, D_HNO3
      REAL*8                 :: D_RMAX, D_VMAX, MAXLOSS, FRACLOSS

      !=================================================================
      ! SETTLE_STRAT_AER begins here!
      !=================================================================

      ! Return if gravitational settling disabled
      IF (.not. LGRAVSTRAT) RETURN

      ! Return if it's the start of the run
      IF ( GET_ELAPSED_SEC() == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! For debugging
      MAXLOSS= 0d0
      D_XHNO3= 0d0
      D_XNAT = 0d0
      D_PPA  = 0d0
      D_TK   = 0d0
      D_XH2O = 0d0
      D_SEDQ = 0d0
      D_PSED = 0d0
      D_ICE0 = 0d0
      D_ICE1 = 0d0
      D_H2O0 = 0d0
      D_H2O1 = 0d0
      D_NAT0 = 0d0
      D_NAT1 = 0d0
      D_HNO3 = 0d0
      D_RMAX = 0d0
      D_VMAX = 0d0

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
!$OMP+PRIVATE( SEDNO3, SEDNAT, SEDICE, SEDPSC, SEDH2O, SEDQUANT )
!$OMP+PRIVATE( XNO3_ABOVE, XNAT_ABOVE, XICE_ABOVE, XPSC_ABOVE )
!$OMP+PRIVATE( XNO3_0, XNAT_0, XICE_0, XPSC_0, VNAT )
!$OMP+PRIVATE( XNO3_BELOW, XNAT_BELOW, XICE_BELOW, XPSC_BELOW )
!$OMP+PRIVATE( FRACLOSS )
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
 
            ! Pressure at center of the level [kPa]
            P       = GET_PCENTER(I,J,L) * 0.1d0

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
         DO K = 1,7
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

         !-----SDE DEBUG
!         WRITE(6,'(a,1x,2I4,1x,E10.2)') 'Pre-STS H2O sum.. for I, J: ',
!     &                      I,J, SUM(STT(I,J,:,IDTH2O))
         !-----SDE DEBUG

         DO L = LLPAR-1,1,-1
            IF (ITS_IN_THE_TROP(I,J,L+1)) CYCLE
            DELZ  = BXHEIGHT(I,J,L)
            DELZ1 = BXHEIGHT(I,J,L+1)

            DO K=1,7
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
         IF (VFALLMAX.gt.TINY(1d0)) THEN
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

            !-----SDE DEBUG
!            WRITE(6,'(a,1x,2I4,1x,E10.2)') 'Post-STS H2O sum. for I, J: ',
!     &                         I,J, SUM(STT(I,J,:,IDTH2O))
!            WRITE(6,'(a,1x,2I4,1x,E10.2)') 'Max fall velocity for I, J: ',
!     &                         I,J,VFALLMAX
            !-----SDE DEBUG

            SEDSTEPLOOP: DO ISED=1,NUMSEDSTEPS
               STARTPT = MAX(1,MINALT-1)
               ! XPSC is the number of molecules tied up in solid particles
               ! per m3 in a grid box
               L = STARTPT
               IF (STATE_PSC(I,J,L).lt.2) THEN
                  XNO3_0 = 0d0
                  XNAT_0 = 0d0
                  XICE_0 = 0d0
               ELSE
                  XNO3_0 = STT(I,J,L,IDTNIT) * AIRMW / 
     &                         (TRACER_MW_G(IDTNIT)*AD(I,J,L))
                  XNAT_0 = XNO3_0! * 4.d0
                  XICE_0 = (KG_AER(I,J,L,I_SPA)-
     &                 ((NATMW/TRACER_MW_G(IDTNIT))*
     &                STT(I,J,L,IDTNIT)))*AIRMW/(ICEMW*AD(I,J,L))
               ENDIF
               !-----SDE DEBUG
               IF (XICE_0 .lt. 0d0) THEN
                  WRITE(6,'(a)') 'ERROR! L,XNO3 XNAT XICE KG AER:'
                  WRITE(6,'(I3,1x,4(E10.4,1x))') L,XNO3_0,XNAT_0,XICE_0,
     &                  KG_AER(I,J,L,I_SPA)
                  CALL ERROR_STOP('Negative mixing ratio!',
     &                               'UCX_mod.F')
               ENDIF
               !-----SDE DEBUG
               XPSC_0 = XNAT_0 + XICE_0
               P_0 = 100.0d0 * GET_PCENTER(I,J,L)
               INVAIR_0 = AIRMW/AD(I,J,L)
               IF (L .lt. (LLPAR-1)) THEN
                  INVAIR_ABOVE = AIRMW/AD(I,J,L+1)
                  IF (STATE_PSC(I,J,L+1).lt.2) THEN
                     XNO3_ABOVE = 0d0
                     XNAT_ABOVE = 0d0
                     XICE_ABOVE = 0d0
                  ELSE
                     XNO3_ABOVE = STT(I,J,L+1,IDTNIT) * INVAIR_ABOVE / 
     &                                   (TRACER_MW_G(IDTNIT))
                     ! NAT = HNO3.3H2O = 4 molecules
                     XNAT_ABOVE = XNO3_ABOVE! * 4.d0
!                     XICE_ABOVE = (KG_AER(I,J,L+1,I_SPA)-STT(I,J,L+1,IDTNIT))
!     &                                          *INVAIR_ABOVE/ICEMW
                     XICE_ABOVE = (KG_AER(I,J,L+1,I_SPA)-
     &                   ((NATMW/TRACER_MW_G(IDTNIT))*
     &                   STT(I,J,L+1,IDTNIT)))*INVAIR_ABOVE/ICEMW
                  ENDIF
                  XPSC_ABOVE = XNAT_ABOVE + XICE_ABOVE
                  P_ABOVE = 100.0d0 * GET_PCENTER(I,J,L+1)
                  PSEDABOVE = VNAT(L+1)*SEDSTEP*1.d2*
     &         (GET_PEDGE(I,J,L+1)-GET_PEDGE(I,J,L+2))/BXHEIGHT(I,J,L+1)
!                  PSEDABOVE=g0 * (1.0d3*AIRMW) * P_ABOVE * 
!     &               VNAT(L+1) / (Rd * T(I,J,L+1))*SEDSTEP
               ELSE
                  INVAIR_ABOVE = 1.0d0
                  XNO3_ABOVE = 0.0d0
                  XNAT_ABOVE = 0.0d0
                  XICE_ABOVE = 0.0d0
                  XPSC_ABOVE = 0.0d0
                  P_ABOVE = GET_PEDGE(I,J,L+2) * 100.0d0
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
                  
                  IF (L.lt.(LLPAR-1)) THEN
                     INVAIR_ABOVE = AIRMW/AD(I,J,L+2)
                     IF (STATE_PSC(I,J,L+2).ge.2) THEN
                        XNO3_ABOVE = STT(I,J,L+2,IDTNIT)*INVAIR_ABOVE/ 
     &                                TRACER_MW_G(IDTNIT)
                        XNAT_ABOVE = XNO3_ABOVE! * 4.d0
                        XICE_ABOVE = (KG_AER(I,J,L+2,I_SPA)-
     &                      ((NATMW/TRACER_MW_G(IDTNIT))*
     &                      STT(I,J,L+2,IDTNIT)))*INVAIR_ABOVE/ICEMW
                     ELSE
                        XNO3_ABOVE = 0d0
                        XNAT_ABOVE = 0d0
                        XICE_ABOVE = 0d0
                     ENDIF
                     XPSC_ABOVE = XNAT_ABOVE + XICE_ABOVE
                     P_ABOVE = GET_PCENTER(I,J,L+2) * 100.0d0
                     PSEDABOVE=VNAT(L+2)*SEDSTEP*1.d2*
     &                  (GET_PEDGE(I,J,L+2)-GET_PEDGE(I,J,L+3))/
     &                  BXHEIGHT(I,J,L+2)
!                     PSEDABOVE = g0 * (1.0d3*AIRMW) * P_ABOVE *
!     &                  VNAT(L+2)/(Rd * T(I,J,L+2))*SEDSTEP
                  ELSE
                     INVAIR_ABOVE = 0.0d0
                     XNO3_ABOVE = 0.0d0
                     XNAT_ABOVE = 0.0d0
                     XICE_ABOVE = 0.0d0
                     XPSC_ABOVE = 0.0d0
                     P_ABOVE = GET_PEDGE(I,J,L+2) * 100.0d0
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
                  IF ((ABOVEGRAD .eq. 0d0) .and. 
     &                            (BELOWGRAD .eq. 0d0)) THEN
                     ! Extremelely unlikely case, but use standard scheme
                     ! for flat vertical profile
                     LOCALPROFILE = 0
                  ELSEIF ((ABOVEGRAD .le. 0d0) .and.
     &                     (BELOWGRAD .le. 0d0)) THEN
                     ! Consistently decreasing with height - above a peak
                     LOCALPROFILE = -1
                  ELSEIF ((ABOVEGRAD .ge. 0d0) .and. 
     &                     (BELOWGRAD .ge. 0d0)) THEN
                     ! Consistently increasing with height - below a peak
                     LOCALPROFILE = +1
                  ELSE
                     ! Local minmum or maximum
                     LOCALPROFILE = 0
                  ENDIF

                  IF (LOCALPROFILE.eq.0) THEN
                     ! Standard sedimentation scheme
                     SEDQUANT = XPSC_0 * PSED0
                     SEDNO3   = XNO3_0 * PSED0
                     SEDICE   = XICE_0 * PSED0
                  ELSEIF (LOCALPROFILE.eq.1) THEN
                     ! Currently below a peak
                     ! Use gradient of bottom and centre boxes
                     IF (XPSC_0 .le. (0.5d0 * XPSC_ABOVE)) THEN
                        ! Interpret as nearing peak; use lower gradient
                        SEDQUANT = (XPSC_0+XPSC_BELOW)*0.5d0*PSED0
                        SEDNO3   = (XNO3_BELOW+XNO3_0)*0.50*PSED0
                        SEDICE   = (XICE_BELOW+XICE_0)*0.50*PSED0
                     ELSE
                        SEDQUANT = (XPSC_ABOVE+XPSC_0)*0.50*PSED0
                        SEDNO3   = (XNO3_ABOVE+XNO3_0)*0.50*PSED0
                        SEDICE   = (XICE_ABOVE+XICE_0)*0.50*PSED0
                     ENDIF
                  ELSE
                     ! Above a peak
                     SEDQUANT = (XPSC_0+XPSC_BELOW)*0.5d0*PSED0
                     SEDNO3   = (XNO3_0+XNO3_BELOW)*0.5d0*PSED0
                     SEDICE   = (XICE_0+XICE_BELOW)*0.5d0*PSED0
                  ENDIF
                  ! Divide sedimenting quantity by the pressure difference
                  ! across the box being sedimented from
                  ! Note conversion from hPa to Pa for denominator, ie
                  ! multiply numerator by (1/100)
                  IF (L .ne. LLPAR) THEN
                     SEDQUANT = (1.d-2) * SEDQUANT /
     &                  (GET_PEDGE(I,J,L+1)-GET_PEDGE(I,J,L+2))
                     SEDNO3 = (1.d-2) * SEDNO3 / 
     &                  (GET_PEDGE(I,J,L+1)-GET_PEDGE(I,J,L+2))
                     SEDICE = (1.d-2) * SEDICE / 
     &                  (GET_PEDGE(I,J,L+1)-GET_PEDGE(I,J,L+2))
                  ELSE
                     ! This shouldn't be possible?
                     CALL ERROR_STOP('Unknown sedimentation error',
     &                               'UCX_mod.F')
                  ENDIF
            
                  ! Apply limits so that sedimented quantity is:
                  ! Greater than or equal to zero
                  ! Less than or equal to the total available quantity

                  ! Note that we are calculating sedimentation using
                  ! the total ice and NAT, but are only actually
                  ! transporting the local NO3
                  SEDQUANT = MAX(0.0d0,SEDQUANT)
                  SEDNO3 = MAX(0.0d0,SEDNO3)
                  SEDICE = MAX(0.0d0,SEDICE)
                  
                  !-----SDE DEBUG
                  ! Convert v/v to kg/box
                  SEDNO3 = SEDNO3*TRACER_MW_G(IDTNIT)/INVAIR_0
                  SEDICE = SEDICE*ICEMW/INVAIR_0
                  
                  SEDNO3 = MIN(SEDNO3,STT(I,J,L+1,IDTNIT))
                  STT(I,J,L,IDTNIT) = STT(I,J,L,IDTNIT) + SEDNO3
                  STT(I,J,L+1,IDTNIT)=STT(I,J,L+1,IDTNIT)-SEDNO3

                  ! Settle the ice out too
                  SEDH2O = MIN(SEDICE,STT(I,J,L+1,IDTH2O))
                  STT(I,J,L,IDTH2O) = STT(I,J,L,IDTH2O) + SEDH2O
                  STT(I,J,L+1,IDTH2O)=STT(I,J,L+1,IDTH2O)-SEDH2O

                  ! Now correct aerosol totals
                  SEDNAT = SEDNO3*NATMW/TRACER_MW_G(IDTNIT)
                  SEDPSC = SEDNAT + SEDICE
                  SEDPSC = MIN(SEDPSC,KG_AER(I,J,L+1,I_SPA))
                  KG_AER(I,J,L,I_SPA) = KG_AER(I,J,L,I_SPA) + SEDPSC
                  KG_AER(I,J,L+1,I_SPA)=KG_AER(I,J,L+1,I_SPA)-SEDPSC
                  !-----SDE DEBUG
                  
                  ! Store for debug
                  FRACLOSS = SEDH2O/STT(I,J,L+1,IDTH2O)
                  IF (FRACLOSS .gt. MAXLOSS) THEN
                     MAXLOSS = FRACLOSS
                     D_ICE0 = XICE_0*ICEMW/INVAIR_0
                     D_ICE1 = D_ICE0 - SEDICE
                     D_H2O1 = STT(I,J,L+1,IDTH2O)
                     D_H2O0 = D_H2O1 + SEDH2O
                     D_NAT1 = STT(I,J,L+1,IDTNIT)*NATMW/
     &                                   TRACER_MW_G(IDTNIT)
                     D_NAT0 = D_NAT1 + SEDNAT
                     D_HNO3 = STT(I,J,L+1,IDTHNO3)
                     D_RMAX = RAD_AER(I,J,L+1,I_SPA)
                     D_VMAX = VNAT(L+1)
                     D_SEDQ = SEDQUANT
                     D_PSED = PSED0
                     D_XH2O = D_H2O0*INVAIR_0/ICEMW
                     D_XHNO3= D_HNO3*INVAIR_0/TRACER_MW_G(IDTHNO3)
                     D_XNAT = D_NAT0*INVAIR_0/NATMW
                     D_PPA  = P_0
                     D_TK   = TEMP
                  ENDIF

               ENDDO SED_LLOOP
            ENDDO SEDSTEPLOOP
         !-----SDE DEBUG
!         WRITE(6,'(a,1x,2I4,1x,E10.2)') 'Post-PSC H2O sum. for I, J: ',
!     &                      I,J, SUM(STT(I,J,:,IDTH2O))
         !-----SDE DEBUG
         ENDIF ! VFALLMAX > 0

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF (LPRT) THEN
         CALL DEBUG_MSG( '### UCX: after SETTLE_STRAT_AER' )
         WRITE(6,'(a)') 'Max ice loss statistics:'
         WRITE(6,'(a,1x,2E12.4)') 'Source ice (before/after) (kg): ',
     &                         D_ICE0, D_ICE1
         WRITE(6,'(a,1x,2E12.4)') 'Source H2O (before/after) (kg): ',
     &                         D_H2O0, D_H2O1
         WRITE(6,'(a,1x,1E12.4)') 'Source HNO3 (total) (kg)      : ',
     &                         D_HNO3
         WRITE(6,'(a,1x,2E12.4)') 'Source NAT (before/after) (kg): ',
     &                         D_NAT0, D_NAT1
         WRITE(6,'(a,1x,3E12.4)') 'Source MR (H2O-HNO3-NAT) (v/v): ',
     &                         D_XH2O, D_XHNO3, D_XNAT
         WRITE(6,'(a,1x,2E12.4)') 'Source T (K) and P (Pa)       : ',
     &                         D_TK,   D_PPA
         WRITE(6,'(a,1x,2E12.4)') 'Particle radius (cm) & V (m/s): ',
     &                         D_RMAX, D_VMAX
         WRITE(6,'(a,1x,2E12.4)') 'Sedquant and PSed variables   : ',
     &                         D_SEDQ, D_PSED
      ENDIF

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

      IF ( LPRT ) CALL DEBUG_MSG( '### UCX: H2SO4 partitioned' )

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
      USE LOGICAL_MOD,      ONLY : LACTIVEH2O
      USE ERROR_MOD,        ONLY : DEBUG_MSG
      USE GRID_MOD,         ONLY : GET_YEDGE
      USE TRACER_MOD,       ONLY : STT, TRACER_MW_G
      USE TRACERID_MOD,     ONLY : IDTSO4, IDTH2O
      USE TRACERID_MOD,     ONLY : IDTHNO3, IDTNIT, IDTBrNO3, IDTHBr
      USE TRACERID_MOD,     ONLY : IDTHOBr, IDTHCl, IDTClNO3, IDTHOCl
      USE CHEMGRID_MOD,     ONLY : ITS_IN_THE_STRATMESO
      USE CHEMGRID_MOD,     ONLY : ITS_IN_THE_CHEMGRID
      USE COMODE_MOD,       ONLY : CSPEC,JLOP
      USE ERROR_MOD,        ONLY : ERROR_STOP, IS_SAFE_DIV
      !-----SDE DEBUG
!      USE ERROR_MOD,        ONLY : IT_IS_NAN
      !-----SDE DEBUG

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

      ! Allow PSC formation outsize Kirner limits?
      LOGICAL,PARAMETER     :: PSC_FULL   =.FALSE.

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

      ! SDE 2013-07-07: Debug variables
      REAL*8                 :: D_XH2O, D_PPA,  D_TK,   D_XHNO3,D_XNAT
      REAL*8                 :: D_MAXN, D_FRACN,D_PP,   D_PS

      !=================================================================
      ! CALC_STRAT_AER begins here!
      !=================================================================

      IF (LPRT) THEN
         CALL DEBUG_MSG( '### UCX: start CALC_STRAT_AER' )
      ENDIF

      ! Partition H2SO4 before proceeding
      CALL CALC_H2SO4_GAS

      D_XH2O = 0d0
      D_PPA  = 0d0
      D_TK   = 0d0
      D_XHNO3= 0d0
      D_XNAT = 0d0
      D_MAXN = 0d0
      D_PP   = 0d0
      D_PS   = 0d0
      
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
!$OMP+PRIVATE( KHET_COMMON, KHET_SPECIFIC, D_FRACN )
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
               IS_VALID = (IS_VALID.or.PSC_FULL)

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

               !-----SDE DEBUG
               ! Revert to static H2O for now
               !H2OSUM = AVGW(I,J,L)
               H2OSUM = STT(I,J,L,IDTH2O)*INVAIR/
     &                   TRACER_MW_G(IDTH2O)
               !-----SDE DEBUG

               ! Calculate partial pressures (Pa)
               HNO3PP = PCENTER_PA * HNO3SUM
               H2OPP  = PCENTER_PA * H2OSUM

               IF (.not.IS_VALID) THEN
                  ! No PSCs (SSA only)
                  STATE_LOCAL = 0
               ELSEIF (TCENTER.gt.T_MAX) THEN
                  ! STS/SSA only
                  STATE_LOCAL = 1
               ELSE
                  ! Calculate saturation pressures
                  PSATHNO3 = (10.d0**(((TSATHNO3_A+(TCENTER*
     &            TSATHNO3_B))*LOG10(H2OPP*760.0d0/101325.0d0))+
     &            (TSATHNO3_C+(TSATHNO3_D*TINV)+(TSATHNO3_E*
     &            TCENTER))))*101325.0d0/760.0d0

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
                  
                  IF (LACTIVEH2O) THEN
                     ! Only interested in sign
                     IF (STATE_LOCAL .gt. 1) THEN
                        H2O_BOX_S = H2OPP-PSATH2O
                     ELSE
                        H2O_BOX_S = H2OPP-PSATH2O_SUPERSAT
                     ENDIF
                  ELSE
                     ! Use local ice mass ratio from GEOS-5 data
                     ! Note that we are using this only for its sign!
                     H2O_BOX_S = QI(I,J,L)
                  ENDIF

                  ! Ice exists/possible?
                  IF (H2O_BOX_S.gt.TINY(0d0)) THEN
                     STATE_LOCAL = 3
                     HNO3_BOX_S = HNO3PP-PSATHNO3
                  ELSE
                     ! If ice not possible could still have NAT
                     H2O_BOX_S = 0d0
                     ! 1. Homogeneous nucleation
                     IF ((LHOMNUCNAT).and.(STATE_LOCAL.eq.1)) THEN
                        HNO3_BOX_S = HNO3PP-PSATHNO3_SUPERCOOL
                     ENDIF
                     ! 2. Box formerly contained ice or NAT
                     IF (STATE_LOCAL.eq.2) THEN
                        HNO3_BOX_S = HNO3PP-PSATHNO3
                     ENDIF
                     IF (HNO3_BOX_S.gt.TINY(1d0)) THEN
                        STATE_LOCAL = 2
                     ELSE
                        STATE_LOCAL = 1
                     ENDIF
                  ENDIF
               ENDIF

               ! Store state
               STATE_PSC(I,J,L) = STATE_LOCAL

               ! Only continue if we want online solid PSCs
               IF (LSOLIDPSC) THEN

                  IF (STATE_LOCAL.eq.3) THEN
                     ! Form ice PSCs
                     IF (LACTIVEH2O) THEN
                        H2O_BOX_S = (H2OPP-PSATH2O)/PCENTER_PA
                        H2O_BOX_S = MAX(0d0,H2O_BOX_S)
                        KG_ICE = H2O_BOX_S*ICEMW*AD(I,J,L)/AIRMW
                     ELSE
                        H2O_BOX_S = QI(I,J,L)*CLDF(L,I,J)*AIRMW/ICEMW
                        KG_ICE = QI(I,J,L)*CLDF(L,I,J)*AD(I,J,L)
                     ENDIF
                     VOL_ICE = H2O_BOX_S * DENAIR * (1.d-3) * ICEMW /
     &                             (DENSICE*AVOGADRO) ! m3 ice/m3 air 
                  ELSE
                     VOL_ICE = 0d0
                     H2O_BOX_S = 0d0
                     KG_ICE = 0d0
                  ENDIF

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

                     ! Debug
                     IF (HNO3SUM.gt.TINY(0d0)) THEN
                        D_FRACN = HNO3_BOX_S/HNO3SUM
                     ELSE
                        D_FRACN = 0d0
                     ENDIF
!
!                     IF (D_FRACN.gt.D_MAXN) THEN
!                        D_MAXN = D_FRACN
                     IF (HNO3_BOX_S.gt.D_XNAT) THEN
                        D_MAXN = D_FRACN
                        D_XH2O = H2OSUM
                        D_XHNO3= HNO3SUM
                        D_XNAT = HNO3_BOX_S
                        D_PPA  = PCENTER_PA
                        D_PP   = HNO3PP
                        D_PS   = PSATHNO3
                        D_TK   = TCENTER
                     ENDIF
                     
                  ELSE
                     HNO3_BOX_S = 0d0
                     VOL_NAT = 0d0
                     KG_NAT = 0d0
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
                     !----SDE DEBUG
!                     IF (IT_IS_NAN(RHO_AER_BOX)) THEN
!                        WRITE(6,'(a,1x,4I4)') 'BAD RHO! I,J,L,STATE:',
!     &                                I,J,L,STATE_LOCAL
!                        WRITE(6,'(a,1x,3E14.4)') 'VOL/DENS/KG ICE : ',
!     &                                 VOL_ICE,DENSICE,KG_NAT
!                        WRITE(6,'(a,1x,3E14.4)') 'VOL/DENS/KG NAT : ',
!     &                                 VOL_NAT,DENSNAT,KG_NAT
!                        CALL ERROR_STOP('HARDKILL','HARDKILL')
!                     ENDIF
                     !----SDE DEBUG
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
               VOL_SLA       = 0d0
               W_H2O         = 0d0
               W_H2SO4       = 1d0

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

                  !-----SDE DEBUG
!                  IF (RAD_AER_BOX.gt.1.d-7) THEN
!                  WRITE(6,'(a,1x,3I3)') 'Aerosol properties at: ', I,J,L
!                  WRITE(6,'(a,1x,4(E10.4,1x))') 'P T H2O H2SO4    : ', 
!     &                          PCENTER,TCENTER,H2OSUM,H2SO4_BOX_L
!                  WRITE(6,'(a,1x,4(E10.4,1x))') 'Ws Wh V RHO      : ', 
!     &                          W_H2SO4,W_H2O,VOL_SLA,RHO_AER_BOX
!                  WRITE(6,'(a,1x,2(E10.4,1x))') 'A(um2/cm3) r(um) : ', 
!     &                          SAD_AER_BOX*1.d8,RAD_AER_BOX*1.d6
!                  ENDIF
                  !-----SDE DEBUG

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

               ! Store liquid fractions
               ! Liquid H2O is removed from the sum, then it is assumed
               ! that the pre-calculated solid H2O is taken out of this
               ! liquid total
               H2O_BOX_L = (98.09d0/18.02d0)*H2SO4_BOX_L*(W_H2O/W_H2SO4)
               H2O_BOX_L = MAX(0d0,MIN(H2O_BOX_L-H2O_BOX_S,H2OSUM))
               H2O_BOX_G = MAX(0d0,H2O_BOX_G - (H2O_BOX_L+H2O_BOX_S)) 
  
               ! If very low number density, ignore settling
               AERFRAC(I,J,L,2) = 0d0 
               AERFRAC(I,J,L,3) = 0d0 
               AERFRAC(I,J,L,4) = 0d0 
               AERFRAC(I,J,L,5) = 0d0 
               AERFRAC(I,J,L,6) = 0d0 
               AERFRAC(I,J,L,7) = 0d0
 
               IF ((HNO3SUM.gt.1d0).and.
     &            (IS_SAFE_DIV(HNO3_BOX_L,HNO3SUM))) THEN
                  AERFRAC(I,J,L,2) = HNO3_BOX_L/HNO3SUM
               ENDIF

               IF ((HClSUM.gt.1d0).and.
     &            (IS_SAFE_DIV(HCl_BOX_L,HClSUM))) THEN
                  AERFRAC(I,J,L,3) = HCl_BOX_L/HClSUM
               ENDIF

               IF ((HOClSUM.gt.1d0).and.
     &            (IS_SAFE_DIV(HOCl_BOX_L,HOClSUM))) THEN
                  AERFRAC(I,J,L,4) = HOCl_BOX_L/HOClSUM
               ENDIF

               IF ((HBrSUM.gt.1d0).and.
     &            (IS_SAFE_DIV(HBr_BOX_L,HBrSUM))) THEN
                  AERFRAC(I,J,L,5) = HBr_BOX_L/HBrSUM
               ENDIF

               IF ((HOBrSUM.gt.1d0).and.
     &            (IS_SAFE_DIV(HOBr_BOX_L,HOBrSUM))) THEN
                  AERFRAC(I,J,L,6) = HOBr_BOX_L/HOBrSUM
               ENDIF

               IF ((H2OSUM.gt.1d0).and.
     &            (IS_SAFE_DIV(H2O_BOX_L,H2OSUM))) THEN
                  AERFRAC(I,J,L,7) = H2O_BOX_L/H2OSUM
               ENDIF

               ! Send properties to larger array
               ! Convert sticking coefficients into
               ! premultiplying factors ((Kirner)
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

!               !-----SDE DEBUG
!               WRITE(6,'(a,3(1x,I3))') '---------I J L       : ',
!     &                      I, J, L
!               WRITE(6,'(a,2(1x,E10.4))') ' T (K) P (Pa)     : ',
!     &                      TCENTER, PCENTER_PA
!               WRITE(6,'(a,3(1x,E10.4))') ' H2SO4 HNO3 H2O   : ',
!     &                      H2SO4SUM, HNO3SUM, H2OSUM
!               WRITE(6,'(a,3(1x,E10.4))') ' HCl H2SO4-l WT%  : ',
!     &                      HClSUM, H2SO4_BOX_L, 1.d2*W_H2SO4
!               WRITE(6,'(a,2(1x,E10.4))') ' r (m) SAD (u2/c3): ',
!     &                      RAD_AER_BOX, SAD_AER_BOX*1.d8
!               WRITE(6,'(a,1(1x,E10.4))') ' RHO (g/cm3)      : ',
!     &                      RHO_AER_BOX
!               WRITE(6,'(a,3(1x,E10.4))') ' GHCl GH2O GHOCl  : ',
!     &                      GAMMA_BOX(3), GAMMA_BOX(4), GAMMA_BOX(8)
!               !-----SDE DEBUG

            ENDDO LLOOP ! L
         ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO

!      !-----SDE DEBUG
!      CALL ERROR_STOP('HARDKILL','HARDKILL')
!      !-----SDE DEBUG

      IF (LPRT) THEN
         CALL DEBUG_MSG( '### UCX: after CALC_STRAT_AER' )
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
         WRITE(6,'(a,1x,E12.4)') 'NAT/NO3 fraction in max-NAT box  : ',
     &                               D_MAXN                       
         WRITE(6,'(a,3(1x,E12.4))') 'Local MR (H2O-HNO3-NAT) (v/v) : ',
     &                               D_XH2O,D_XHNO3,D_XNAT        
         WRITE(6,'(a,3(1x,E12.4))') 'Pressure (box-HNO3-sat) (Pa)  : ',
     &                               D_PPA,D_PP,D_PS        
         WRITE(6,'(a,1x,E12.4)') 'Local temperature (K)         : ',
     &                               D_TK                    
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
                                                ! 1 = Liquid aerosols
                                                ! 2 = Solid PSC
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
      SUBROUTINE CHECK_NAT (I,J,L,IS_NAT,IS_PSC,IS_STRAT)
!
! !USES:
!
      USE TRACERID_MOD,      ONLY : IDTNIT
      USE TRACER_MOD,        ONLY : STT
      USE CHEMGRID_MOD,      ONLY : ITS_IN_THE_STRAT
      USE LOGICAL_MOD,       ONLY : LPSCCHEM
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: I,J,L      ! Grid indices
!
! !OUTPUT VARIABLES:
!
      LOGICAL,INTENT(OUT)         :: IS_NAT     ! Is surface NAT?
      LOGICAL,INTENT(OUT)         :: IS_PSC     ! Are there solid PSCs?
      LOGICAL,INTENT(OUT)         :: IS_STRAT   ! Are we in the strat?
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

      IS_STRAT = ITS_IN_THE_STRAT(I,J,L)
      IS_PSC = ((LPSCCHEM).and.(STATE_PSC(I,J,L).ge.2).and.(IS_STRAT))
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

      ! Temporary - for debug
      USE ERROR_MOD,     ONLY : IT_IS_NAN,ERROR_STOP     ! Test for NaN
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
      
      ! Debug message
      CHARACTER(LEN=255)   :: DBGMSG

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

      ! Carslaw only valid for 2e-5 < PPH2O < 2e-3 (hPa)
      PATMH2O = MAX(PATMH2O,1.9738465d-8)
      PATMH2O = MIN(PATMH2O,1.9738465e-6)

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

      ! SDE DEBUG
      IF (W_H2SO4.gt.1.0d0) THEN
         WRITE(DBGMSG,'(a,E10.4)') 'PATMH2O: ', PATMH2O
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,E10.4)') 'TCENTER: ', TCENTER
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,E10.4)') 'TCENTER_IN: ', TCENTER_IN
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,E10.4)') 'MOLDENS_H2SO4: ', MOLDENS_H2SO4
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,E10.4)') 'H2SO4SUM: ', H2SO4SUM
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,E10.4)') 'PATMH2O: ', PATMH2O
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,E10.4)') 'X_H2SO4_BIN: ', X_H2SO4_BIN
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,E10.4)') 'M_H2SO4_BIN: ', M_H2SO4_BIN
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,E10.4)') 'M_H2SO4: ', M_H2SO4
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,E10.4)') 'W_H2SO4:    ', W_H2SO4
         CALL DEBUG_MSG( DBGMSG )
         CALL ERROR_STOP('BAD WEIGHT FRAC','UCX_mod')
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
      USE ERROR_MOD,       ONLY : ERROR_STOP
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

      ! Velocities should be of the order of 0.1 m/s
      IF (VEL.gt.10.d0) THEN
         CALL ERROR_STOP(' Excessive fall velocity? ', 
     &                   ' CALC_FALLVEL, UCX_mod')
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
      ! Temporary - for debug
      USE ERROR_MOD,     ONLY : IT_IS_NAN,ERROR_STOP     ! Test for NaN
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

      ! Control whether to run calculations
      LOGICAL              :: HClOK, HOBrOK

      ! Debug variables
      INTEGER              :: I
      CHARACTER(LEN=255)   :: DBGMSG

      !=================================================================
      ! CALC_SLA_GAMMA begins here!
      !=================================================================

      PSATH2O = EXP(18.452406985d0-3505.1578807d0/T-330918.55082d0/(T*T)
     &     +12725068.262d0/(T*T*T))       ! Saturation pressure of H2O
      H2OPP = H2OSUM * P                  ! Partial pressure of H2O
      ACTH2O = MAX((H2OPP/PSATH2O),1.0d0) ! Water activity

      ! Calculate molality of solution
      !WT = MIN(100.0d0,100.0d0*WT_FRC) ! Convert from fraction to %
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

      HClPP = HClSUM*P/1013.25d0 ! Note atm, not hPa
      ClNO3PP = ClNO3SUM*P/1013.25d0
      BrNO3PP = BrNO3SUM*P/1013.25d0
      HOBrPP = HOBrSUM*P/1013.25 ! Note atm, not hPa
 
      ! Should we bother running calculations?
      HClOK = (HClPP .gt. 1.d-30)
      HOBrOK = (HOBrPP .gt. 1.d-30)

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
      if (lClNO3.gt.(1.d5*arad)) then
         ! Limiting rate
         fClNO3   =  arad/(3.d0*lClNO3)
      else
         fClNO3   =  1.D0/tanh(arad/lClNO3)- lClNO3/arad
      endif
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

      IF (HClOK) THEN
         RXNGAMMA(3) =  gClNO3_H2O
         RXNGAMMA(4) =  gClNO3_HCl
      ELSE
         RXNGAMMA(3) = TINY(1d0)
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
         if (lHOCl.gt.(1.d5*arad)) then
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
      IF ((HClOK).and.(HOBrOK)) THEN
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

      ! SDE 2013-10-18: DEBUG
      DO I=1,11
         IF (IT_IS_NAN(RXNGAMMA(I))) THEN
            WRITE(DBGMSG,'(a,I2)') 'RXNGAMMA NaN: ', I
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'NDENS: ', NDENS
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'T: ', T
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'P: ', P
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'WT_FRC: ', WT_FRC
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'H2OSUM: ', H2OSUM
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'HClSUM: ', HClSUM
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'HBrSUM: ', HBrSUM
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'HOBrSUM: ', HOBrSUM
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'ClNO3SUM: ', ClNO3SUM
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'BrNO3SUM: ', BrNO3SUM
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'RHO: ', RHO
            CALL DEBUG_MSG( DBGMSG )
            WRITE(DBGMSG,'(a,E10.4)') 'ARAD: ', ARAD
            CALL DEBUG_MSG( DBGMSG )
            CALL ERROR_STOP('BAD GAMMA','UCX_mod')
         ENDIF
      ENDDO

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
!               MIT Laboratory for Aviation and the Environment               !
!               MIT Laboratory for Aviation and the Environment               !
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
!               MIT Laboratory for Aviation and the Environment               !
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
     &                   ' reading STATE_PSC, UCX_mod')

      ENDIF

 555  CONTINUE

      ! Close file
      CLOSE( IU_RST )      

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_PSC_FILE: read file' )

      END SUBROUTINE READ_PSC_FILE
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_clock_trac
!
! !DESCRIPTION: Subroutine SET\_CLOCK\_TRAC sets the clock tracer mixing 
!  ratio within the bottom five grid levels, increasing by a fixed rate of
!  0.5 ppbv/day
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CLOCK_TRAC(STEPLEN)
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTCLOCK
      USE TRACER_MOD,   ONLY : STT
      USE CMN_SIZE_MOD
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN) :: STEPLEN ! Step length (min)
!
! !REVISION HISTORY: 
!  16 Oct 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,J,L,LCLOCK
      REAL*8, PARAMETER :: INCRATE = 0.5d-9 ! vv/day increase
      
      !=================================================================
      ! SET_CLOCK_TRAC begins here!
      !=================================================================

      IF (CLOCKMR.lt.TINY(1d0)) THEN
         ! Reset clock tracer everywhere
         LCLOCK = LLPAR
      ELSE
         ! Only set bottom 5 levels
         LCLOCK = 5
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO J=1,JJPAR
      DO I=1,IIPAR
      DO L=1,LCLOCK
         STT(I,J,L,IDTCLOCK) = CLOCKMR
      ENDDO ! L
      ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO

      ! Increment clock mixing ratio for next step
      CLOCKMR = CLOCKMR + (INCRATE*STEPLEN/(60d0*24d0))

      END SUBROUTINE SET_CLOCK_TRAC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_h2o_trac
!
! !DESCRIPTION: Subroutine SET\_H2O\_TRAC sets the H2O tracer throughout 
!  the selected domain (either troposphere only or the full grid).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_H2O_TRAC ( SETSTRAT )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTH2O
      USE TRACER_MOD,   ONLY : STT
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE DAO_MOD,      ONLY : SPHU, AD, T, RH
      USE CHEMGRID_MOD, ONLY : ITS_IN_THE_TROP
      USE LOGICAL_MOD,  ONLY : LACTIVEH2O
      USE CMN_SIZE_MOD
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)          :: SETSTRAT   ! Set strat H2O?
!
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,J,L

      REAL*8            :: ESAT, SHMB, PRES, TEMP
      REAL*8, PARAMETER :: A =  23.5518d0
      REAL*8, PARAMETER :: B =  2937.4d0
      REAL*8, PARAMETER :: C = -4.9283d0

      LOGICAL           :: READ_SPHU

      !=================================================================
      ! SET_H2O_TRAC begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PRES, TEMP, ESAT, SHMB, READ_SPHU )
      DO J=1,JJPAR
      DO I=1,IIPAR
      DO L=1,LLPAR
         READ_SPHU = (ITS_IN_THE_TROP(I,J,L).or.SETSTRAT.or.
     &                             (.not.LACTIVEH2O))
         IF (READ_SPHU) THEN
            STT(I,J,L,IDTH2O) = SPHU(I,J,L) * 1.d-3 * AD(I,J,L)
         ELSE
            ! Set specific humidity to transported H2O
            SPHU(I,J,L) = STT(I,J,L,IDTH2O) * 1.d3 / AD(I,J,L)

            ! Pressure at midpoint of box (I,J,L)
            PRES = GET_PCENTER(I,J,L)

            ! Temperature at grid box (I,J,L)
            TEMP = T(I,J,L)

            ! Saturation water vapor pressure in mbar 
            ! (from NASA GTE PEM-Tropics handbook)
            ESAT = ( 10d0**( A - ( B / TEMP ) ) ) * ( TEMP**C )
               
            ! Specific humidity in mb
            SHMB = SPHU(I,J,L) * 1.6072d-3 * PRES
               
            ! Relative humidity as a percentage
            RH(I,J,L) = ( SHMB / ESAT ) * 100d0 
         ENDIF
      ENDDO ! L
      ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO

      END SUBROUTINE SET_H2O_TRAC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ucx_h2so4phot
!
! !DESCRIPTION: Subroutine UCX\_H2SO4PHOT propagates the calculated H2SO4
!  photolysis (J) rate at the top of the chemistry grid through to the top
!  of the transport grid, approximating H2SO4 photolysis in the mesosphere.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE UCX_H2SO4PHOT( SUNCOS )
!
! !USES:
!
      USE CHEMGRID_MOD,    ONLY : GET_CHEMGRID_LEVEL
      USE COMODE_LOOP_MOD, ONLY : IRM, NKSO4PHOT, NCSGAS
      USE COMODE_LOOP_MOD, ONLY : DEFPRAT, NAMEGAS, NRATES
      USE TIME_MOD,        ONLY : GET_TS_CHEM
      USE TRACER_MOD,      ONLY : STT,TRACER_MW_G
      USE TRACERID_MOD,    ONLY : IDTSO2,IDTSO4
      USE FAST_JX_MOD,     ONLY : FJXFUNC
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN)           :: SUNCOS(MAXIJ) ! Cosine of SZA
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  17 Aug 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: I, J, L, IJWINDOW
      REAL*8                :: GMU,SO4_IN,PHOTDELTA,DTCHEM
      REAL*8                :: SO4_DELTA
      LOGICAL               :: DAYCOLUMN
      INTEGER               :: LMINPHOT
      INTEGER,SAVE          :: LOCAL_NKSO4PHOT
      CHARACTER*8,SAVE      :: LOCAL_NAME
      INTEGER,SAVE          :: LOCAL_IPHOT
      INTEGER,SAVE          :: IFNC,IBRCH
      REAL*8,SAVE           :: RELWT
      LOGICAL,SAVE          :: FIRST=.TRUE.
      INTEGER               :: ICS

      !=================================================================
      ! UCX_H2SO4PHOT begins here!
      !=================================================================

      IF (FIRST) THEN
         FIRST = .FALSE.
         DO ICS = 1,NCSGAS
            IF (NKSO4PHOT(ICS).ne.0) THEN
               LOCAL_NKSO4PHOT = NKSO4PHOT(ICS)
               LOCAL_NAME = NAMEGAS(IRM(1,LOCAL_NKSO4PHOT,ICS))
               LOCAL_IPHOT = LOCAL_NKSO4PHOT - NRATES(ICS)
               IFNC = DEFPRAT(LOCAL_NKSO4PHOT,ICS) + 0.01d0
               IBRCH = 10.d0*(DEFPRAT(LOCAL_NKSO4PHOT,ICS)-IFNC) + 0.5d0
               RELWT = TRACER_MW_G(IDTSO2)/TRACER_MW_G(IDTSO4)
            ENDIF
         ENDDO 
      ENDIF

      ! Allow for the possibility of variable timestep
      DTCHEM = GET_TS_CHEM()*60.d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, IJWINDOW, DAYCOLUMN, GMU )
!$OMP+PRIVATE( LMINPHOT, PHOTDELTA, SO4_IN, SO4_DELTA )
      DO I=1,IIPAR
      DO J=1,JJPAR
      IJWINDOW = (J-1)*IIPAR + I
      GMU = SUNCOS(IJWINDOW)
      DAYCOLUMN = (GMU.gt.0d0)
      IF (DAYCOLUMN) THEN
         LMINPHOT = GET_CHEMGRID_LEVEL(I,J)
         ! Retrieve photolysis rate as a fraction of gaseous SO4
         PHOTDELTA = FJXFUNC(I,J,LMINPHOT,LOCAL_IPHOT,IBRCH,LOCAL_NAME)
     &                                                         * DTCHEM
         PHOTDELTA = MIN(1.d0,PHOTDELTA)
         DO L=LMINPHOT+1,LLPAR
            ! Apply photolysis to SO4
            ! First retrieve gaseous fraction
            SO4_IN = STT(I,J,L,IDTSO4)*SO4_PHOTFRAC(I,J,L)
            SO4_DELTA = PHOTDELTA*SO4_IN
            ! Remove from SO4
            STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4) - SO4_DELTA
            ! Add to SO2. Note change in molar mass
            STT(I,J,L,IDTSO2) = STT(I,J,L,IDTSO2) + (SO4_DELTA*RELWT)
         ENDDO
      ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE UCX_H2SO4PHOT
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
!!  line of description. Routine length is 51 lines.
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
! !IROUTINE: init_ucx
!
! !DESCRIPTION: Subroutine INIT\_UCX initializes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_UCX( am_I_Root, YYYYMMDD, HHMMSS )
!
! !USES:
!
      USE TRACERID_MOD
      USE LOGICAL_MOD
      USE TIME_MOD,     ONLY : GET_YEAR
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
      LOGICAL :: IT_EXISTS, USE2DDATA
      CHARACTER(LEN=255) :: DBGMSG, GRIDSPEC

      !=================================================================
      ! INIT_UCX begins here!
      !=================================================================

      WRITE( 6,'(a)') REPEAT( '=', 79 )
      WRITE( 6,'(a)') 'U N I F I E D   C H E M I S T R Y'
      WRITE( 6,'(a)') 'Routines written by SEBASTIAN D. EASTHAM'
      WRITE( 6,'(a)') REPEAT( '=', 79 )

      USE2DDATA = (LSETOCS.or.LSETCFC.or.LSETCL.or.LSETBR.or.LSETH2SO4)
      USE2DDATA = (USE2DDATA.or.STRAT2DCH4.or.LSETBRSTRAT)
      USE2DDATA = (USE2DDATA.or.LSETNOYSTRAT)
      USE2DDATA = (USE2DDATA.or.LSETN2O.or.LBASICEMIS)

#if !defined(GRID4x5) && !defined(GRID2x25)
      IF (LUCX.and.USE2DDATA) THEN
       DBGMSG = 'Zonal means preprocessed only for 2x2.5 and 4x5 grid.'
       CALL ERROR_STOP( DBGMSG, 'INIT_UCX (UCX_mod.f)!' )
      ELSE
       AVG_FILE_ROOT = ''
      ENDIF
#elif defined(GRID2x25)
      GRIDSPEC = 'Grid2x25/InitCFC_'
#elif defined(GRID4x5)
      GRIDSPEC = 'Grid4x5/InitCFC_'
#endif

      !-----SDE DEBUG
      ! Check - and write - the `record length' of AVG_FILE_ROOT
      WRITE(6,*) TRIM(TRAC_FILE_ROOT)
      WRITE(6,*) '/DiurnalAvg/'
      WRITE(6,*) TRIM(GRIDSPEC)
      WRITE(6,*) TRIM(TRAC_FILE_ROOT), '/DiurnalAvg/', TRIM(GRIDSPEC)
      !-----SDE DEBUG

      ! Determine folder paths from root folder
      WRITE(     AVG_FILE_ROOT,'(a,a,a)') TRIM(TRAC_FILE_ROOT),
     &   '/DiurnalAvg/', TRIM(GRIDSPEC)
      WRITE(    NOON_FILE_ROOT,'(a,a,a)') TRIM(TRAC_FILE_ROOT), 
     &   '/NoonTime/', TRIM(GRIDSPEC)
      WRITE(MONTREAL_FILE_ROOT,'(a,a)') TRIM(TRAC_FILE_ROOT), 
     &   '/Montreal/'
  
      IF ( LPRT ) THEN
         WRITE(DBGMSG,'(a,a)') '### UCX: Reading O1D/O3P from ',
     &     TRIM(NOON_FILE_ROOT)
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,a)') '### UCX: Initializing CFCs from ',
     &     TRIM(AVG_FILE_ROOT)
         CALL DEBUG_MSG( DBGMSG )
         WRITE(DBGMSG,'(a,a)') '### UCX: Projecting CFCs from ',
     &     TRIM(MONTREAL_FILE_ROOT)
         CALL DEBUG_MSG( DBGMSG )
      ENDIF
      
      MONTREAL_FILE_ROOT = TRIM(TRAC_FILE_ROOT)
     &      // '/Montreal/'

      ! Allocate array of input pressure levels
      ALLOCATE( AER_PLEVS( AER_NLEVS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AER_PLEVS' )

      ! Set input pressure levels (hPa)
      AER_PLEVS = (/ 0.2200d+00, 0.2600d+00, 0.3100d+00, 0.3600d+00,
     &               0.4300d+00, 0.5100d+00, 0.6000d+00, 0.7100d+00,
     &               0.8400d+00, 0.9900d+00, 1.1700d+00, 1.3800d+00,
     &               1.6300d+00, 1.9300d+00, 2.2800d+00, 2.6900d+00,
     &               3.1800d+00, 3.7600d+00, 4.4400d+00, 5.2500d+00,
     &               6.2000d+00, 7.3200d+00, 8.6500d+00, 1.0220d+01,
     &               1.2070d+01, 1.4260d+01, 1.6850d+01, 1.9910d+01,
     &               2.3520d+01, 2.7780d+01, 3.2820d+01, 3.8770d+01,
     &               4.5810d+01, 5.4110d+01, 6.3930d+01, 7.5220d+01,
     &               8.9220d+01, 1.0540d+02, 1.2451d+02, 1.4710d+02,
     &               1.7377d+02, 2.0529d+02, 2.4252d+02, 2.8650d+02,
     &               3.3847d+02, 3.9985d+02, 4.7237d+02, 5.5804d+02,
     &               6.5924d+02, 7.7880d+02, 9.2004d+02 /)
 
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

      ! Set clock tracer mixing ratio to zero
      CLOCKMR = 0d0

      ! Set "current" month to impossible value
      SFC_MONTH = 0
      OFFSETYEAR = 0
      IF (.not.LFUTURECFC) THEN
         CFCYEAR = GET_YEAR()
      ENDIF

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
         TRAC_EMIT(N_EMIT) = 'CCl4'

         ! CH3Cl
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3Cl
         TRAC_EMIT(N_EMIT) = 'CH3Cl'

         ! CH3CCl3
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCH3CCl3
         TRAC_EMIT(N_EMIT) = 'CH3CCl3'
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
         TRAC_EMIT(N_EMIT) = 'CH3Br'
      ENDIF

      IF (IDTH2.gt.0) THEN
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH2
         TRAC_EMIT(N_EMIT) = 'H2'
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

      ! Initialize NOx coefficient arrays
      ALLOCATE( NOX_O( IIPAR, JJPAR, LLPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOX_O' )
      NOX_O = 0d0

      ALLOCATE( NOX_J( IIPAR, JJPAR, LLPAR, 4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOX_J' )
      NOX_J = 0d0

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
      !          7 - H2O
      ALLOCATE( AERFRAC( IIPAR, JJPAR, LLPAR, 7 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AERFRAC' )
      AERFRAC = 0d0

      ALLOCATE( AERFRACIND( 7 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AERFRACIND' )
      AERFRACIND(1) = IDTSO4
      AERFRACIND(2) = IDTHNO3
      AERFRACIND(3) = IDTHCl
      AERFRACIND(4) = IDTHOCl
      AERFRACIND(5) = IDTHBr
      AERFRACIND(6) = IDTHOBr
      AERFRACIND(7) = IDTH2O
    
      ! Sticking coefficients for PSC reactions on SLA 
      ALLOCATE( KHETI_SLA( IIPAR, JJPAR, LLPAR, 11 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'KHETI_SLA' )
      KHETI_SLA = 0d0
     
      ! Read in PSC state
      CALL READ_PSC_FILE( am_I_Root, YYYYMMDD,  
     &                    HHMMSS,    IT_EXISTS )

      IF ((am_I_Root).and.(.not.IT_EXISTS)) THEN
         WRITE(6,*) 
     &  '    - UCX: PSC restart not found, initialize PSC-free'
      ENDIF

      ! H2SO4 photolysis rate at the top of the chemgrid
      ALLOCATE( SO4_TOPPHOT( IIPAR,JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_TOPPHOT' )
      SO4_TOPPHOT = 0.d0

      END SUBROUTINE INIT_UCX
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_ucx
!
! !DESCRIPTION: Subroutine CLEANUP\_UCX deallocates module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_UCX
!
! !REVISION HISTORY: 
!  04 Apr 2013 - S. D. Eastham - Initial version
!  13 Apr 2013 - S. D. Eastham - Added PSC arrays
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_UCX begins here!
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
      IF ( ALLOCATED( AER_MR     ) ) DEALLOCATE( AER_MR     )
      IF ( ALLOCATED( AER_PLEVS  ) ) DEALLOCATE( AER_PLEVS  )
      IF ( ALLOCATED( NOX_O      ) ) DEALLOCATE( NOX_O      )
      IF ( ALLOCATED( NOX_J      ) ) DEALLOCATE( NOX_J      )
      IF ( ALLOCATED( SO4_TOPPHOT) ) DEALLOCATE( SO4_TOPPHOT)

      END SUBROUTINE CLEANUP_UCX
!
!EOC
      END MODULE UCX_MOD
