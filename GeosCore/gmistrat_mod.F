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
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: SET_INITIAL_MIXRATIOS
      PUBLIC  :: EMISS_BASIC
      PUBLIC  :: SETTLE_STRAT_AER
      PUBLIC  :: CALC_H2SO4_GAS
      PUBLIC  :: SO4_PHOTFRAC
      PUBLIC  :: CALC_PSC
      PUBLIC  :: INIT_GMISTRAT
      PUBLIC  :: CLEANUP_GMISTRAT
!
! PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: APPLY_2DTRAC
      PRIVATE :: READ_SFC
      PRIVATE :: GET_MONTREAL
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
      ! N_PSC           : Number of PSC types
      ! I_STS           : Index of STS PSCs
      ! I_PRT           : Index of particulate PSCs
      ! AVOGADRO        : Avogadro's number (#/mole)
      !
      !=================================================================

      INTEGER, PARAMETER            :: AER_NLEVS=51
      REAL*8, PARAMETER             :: T_STS=240.0d0
      INTEGER, PARAMETER            :: N_PSC=2
      INTEGER, PARAMETER            :: I_STS=1
      INTEGER, PARAMETER            :: I_PRT=2
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
      !
      ! Arrays
      !
      ! AER_MR             : AER 2D mixing ratios (v/v)
      ! AER_PLEVS          : Pressure levels of AER data (hPa)
      ! IDT_EMIT           : IDs of emitted tracers
      ! GRID_EMIT          : Surface CFC mixing ratios
      ! SFC_MONTH          : Month of last read-in
      ! H2SO4_GASFRAC      : Fraction of SO4 considered gaseous
      ! RAD_PSC            : PSC aerosol radius (m)
      ! KG_PSC             : PSC mass (kg/box)
      ! SAD_PSC            : PSC surface area density (m2/m3)
      ! OD_PSC             : PSC optical depth (unitless)
      ! NDENS_PSC          : PSC aerosol number density (#/m3)
      ! RHO_PSC            : PSC aerosol mass density (kg/m3 aerosol)
      ! STATE_PSC          : Gridbox PSC type (see Kirner et al)
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

      ! Arrays
      REAL*8,DIMENSION(:,:),ALLOCATABLE     :: AER_MR
      REAL*8,DIMENSION(:),ALLOCATABLE       :: AER_PLEVS
      INTEGER,DIMENSION(MAX_EMIT)           :: IDT_EMIT
      REAL*8,DIMENSION(:,:),ALLOCATABLE     :: GRID_EMIT
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE   :: H2SO4_GASFRAC
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: RAD_PSC
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: KG_PSC
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: SAD_PSC
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: OD_PSC
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: NDENS_PSC
      REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE :: RHO_PSC
      INTEGER,DIMENSION(:,:,:),ALLOCATABLE  :: STATE_PSC
 
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
!     (25) BrCl                = BrCl                 < (LSETCL || LSETBR)
!     (26) Br2                 = Br2                  < LSETBR
!     (27) Br                  = Br                   < LSETBR
!     (28) BrO                 = BrO                  < LSETBR
!     (29) HOBr                = HOBr                 < LSETBR
!     (30) HBr                 = HBr                  < LSETBR
!     (31) BrNO2               = BrNO2                < LSETBR
!     (32) BrNO3               = BrNO3                < LSETBR
!     (33) CHBr3               = CHBr3                < LSETBR
!     (34) CH2Br2              = CH2Br2               < LSETBR
!     (35) CH3Br               = CH3Br                < LSETBR
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
      USE LOGICAL_MOD,          ONLY : LSETH2SO4,LGMISTRAT
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
      USE2DDATA = (USE2DDATA.or.STRAT2DCH4)
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
         IF (LSETBR) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|Br2      |Br       |BrO      |'
               WRITE(6,*) '|HOBr     |HBr      |BrNO2    |'
               WRITE(6,*) '|BrONO2   |CHBr3    |CH2Br2   |'
               WRITE(6,*) '|CH3Br    |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTBr2
            TRAC_SET  = 'BR2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         
            TRAC_IDX = IDTBr
            TRAC_SET  = 'BR'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         
            TRAC_IDX = IDTBrO
            TRAC_SET  = 'BRO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         
            TRAC_IDX = IDTHOBr
            TRAC_SET  = 'HOBR'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         
            TRAC_IDX = IDTHBr
            TRAC_SET  = 'HBR'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTBrNO2
            TRAC_SET  = 'BRNO2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTBrNO3
            TRAC_SET  = 'BRNO3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCHBr3 
            TRAC_SET  = 'CHBR3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCH2Br2
            TRAC_SET  = 'CH2BR2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCH3Br
            TRAC_SET  = 'CH3BR'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         ENDIF

         ! Br-Cl (special case)
         IF (LSETCL.or.LSETBR) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|BrCl     |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTBrCl
            TRAC_SET  = 'BRCL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
         ENDIF

         ! Chlorine species (except CFCS, HCFCs and Halons)
         IF (LSETCL) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|Cl2      |ClOx     |Cl2O2    |'
               WRITE(6,*) '|ClNO2    |ClONO2   |HCl      |'
               WRITE(6,*) '|HOCl     |OClO     |CH3Cl    |'
               WRITE(6,*) '|CH3CCl3  |CCl4     |ClOO     |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTClOx
            TRAC_SET  = 'CL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
            TRAC_SET  = 'CLO'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
     
            TRAC_IDX = IDTCl2
            TRAC_SET  = 'CL2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCl2O2
            TRAC_SET  = 'CL2O2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTClNO2
            TRAC_SET  = 'CLNO2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)

            TRAC_IDX = IDTClNO3
            TRAC_SET  = 'CLNO3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)

            TRAC_IDX = IDTHCl
            TRAC_SET  = 'HCL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTHOCl
            TRAC_SET  = 'HOCL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTOClO
            TRAC_SET  = 'OCLO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)

            TRAC_IDX = IDTClOO
            TRAC_SET  = 'CLOO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)

            ! Handle chloroalkanes      
            TRAC_IDX = IDTCH3Cl
            TRAC_SET  = 'CH3CL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCH3CCl3
            TRAC_SET  = 'CH3CCL3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCCl4
            TRAC_SET  = 'CCL4'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         ENDIf
         
         ! CFCs, HCFCs and Halons
         IF (LSETCFC) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|H-1211   |H-1301   |H-2402   |'
               WRITE(6,*) '|CFC-11   |CFC-12   |CFC-X    |'
               WRITE(6,*) '|HCFC-X   |H-1202   |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            ! H1202, 1211, 1301 and 2402
            TRAC_IDX = IDTH1202 
            TRAC_SET  = 'CBR2F2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTH1211  
            TRAC_SET  = 'CBRCLF2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTH1301 
            TRAC_SET  = 'CBRF3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTH2402  
            TRAC_SET  = 'C2BR2F4'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            ! CFC  11, 12, and [113/114/115]
            TRAC_IDX = IDTCFC11
            TRAC_SET  = 'CFCL3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)   
    
            TRAC_IDX = IDTCFC12 
            TRAC_SET  = 'CF2CL2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
   
            TRAC_IDX = IDTCFC
            TRAC_SET  = 'C2CL3F3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.) 
            TRAC_SET  = 'C2CL2F4'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
            TRAC_SET  = 'C2CLF5'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
     
            ! HCFC 22, 141b and 142b
            TRAC_IDX = IDTHCFC
            TRAC_SET  = 'CHCLF2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
            TRAC_SET  = 'C2H3FCL2'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
            TRAC_SET  = 'C2H3F2CL'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
         ENDIF

         IF (LSETOCS) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|OCS      |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTOCS
            TRAC_SET  = 'OCS'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         ENDIF

         IF (LSETH2SO4) THEN
            IF ( am_I_Root ) THEN
               WRITE(6,*) '|H2SO4    |         |         |'
               WRITE(6,*) '|---------|---------|---------|'
            ENDIF
            TRAC_IDX = IDTSO4
            ! Actually H2SO4 aerosol
            TRAC_SET  = 'SO4'
            CALL APPLY_2DTRAC(.FALSE.,.TRUE.)
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
         CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
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
            CALL APPLY_2DTRAC(.FALSE.,.TRUE.)
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
      SUBROUTINE APPLY_2DTRAC( TRAC_ADD, STRAT_ONLY )
!
! !USES:
!
      USE FILE_MOD,      ONLY : IOERROR
      USE PRESSURE_MOD,  ONLY : GET_PCENTER
      USE TRACER_MOD,    ONLY : STT, TCVV
      USE CHEMGRID_MOD,  ONLY : ITS_IN_THE_STRATMESO
      USE DAO_MOD,       ONLY : AD
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)              :: TRAC_ADD     ! Add (instead of overwrite) 
      LOGICAL, INTENT(IN)              :: STRAT_ONLY   ! Only write to stratosphere 
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

      !=================================================================
      ! APPLY_2DTRAC begins here!
      !=================================================================

      TRAC_FILE = TRIM(TRAC_FILE_ROOT)//TRIM(TRAC_SET )
      TRAC_FILE = TRIM(TRAC_FILE)//'_'//TRIM(MONTH_STR)//'.dat'

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
 
      ! Scan through STT, element by element
      DO J = 1,JJPAR
         DO I = 1,IIPAR
            ! Vertcount is the layer count for the input, where layer 1
            ! is at the *top* of the atmosphere
            VERTCOUNT = 51
            EXTRAP = .TRUE.
            DO L = 1, LLPAR
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

      ! Allocated module variables, if necessary
      IF (.not.ALLOCATED(GRID_EMIT)) THEN
         ! Set "current" month to impossible value
         SFC_MONTH = 0

         CALL INIT_GMISTRAT
      ENDIF

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
            CASE ('GENERIC_CFC')
               ! Add 3 CFCs
               LOCAL_NAME = TRIM('CFC113')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('CFC114')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('CFC115')
               CALL GET_MONTREAL(LOCAL_NAME  ,TARG_LINE,N)
            CASE ('GENERIC_HCFC')
               ! Add 3 HCFCs
               LOCAL_NAME = TRIM('HCFC22')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('HCFC141b')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
               LOCAL_NAME = TRIM('HCFC142b')
               CALL GET_MONTREAL(LOCAL_NAME,TARG_LINE,N)
            CASE ('H2402_YEARLY')
               ! Limited data available
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
!  of stratospheric sulfates, nitrates, black carbon and organic carbon on a
!  total mass basis. It is copied largely from GRAV\_SETTLING in
!  sulfate\_mod.F, using a mean radius of 0.1 micrometers for STS (which
!  all stratospheric sulfate is assumed to be) and the relevant aerosol radii
!  calculated in aerosol\_mod.F for the rest. All of this is ignored if APM
!  is active.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETTLE_STRAT_AER
!
! !USES:
!
      USE DAO_MOD,       ONLY : T, BXHEIGHT
      USE PRESSURE_MOD,  ONLY : GET_PCENTER
      USE TRACER_MOD,    ONLY : STT, XNUMOL
      USE TRACERID_MOD,  ONLY : IDTSO4,IDTBCPI
      USE TIME_MOD,      ONLY : GET_ELAPSED_SEC, GET_TS_CHEM
      USE COMODE_MOD,    ONLY : JLOP,WERADIUS,ERADIUS
      USE LOGICAL_MOD,   ONLY : LGMISTRAT
      USE CHEMGRID_MOD,  ONLY : ITS_IN_THE_STRATMESO
      USE CHEMGRID_MOD,  ONLY : ITS_IN_THE_TROP

      USE CMN_SIZE_MOD        ! Size parameters
      USE CMN_DIAG_MOD        ! ND44
      USE CMN_GCTM_MOD        ! g0
      USE ERROR_MOD,   ONLY : IT_IS_NAN,ERROR_STOP     ! Test for NaN
!
! !REVISION HISTORY: 
!  11 Apr 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Common to all particles
      INTEGER                :: I,      J,     L,        DTCHEM
      REAL*8                 :: DELZ,   DELZ1
      REAL*8                 :: P,      TEMP,  DP,       PDP     
      REAL*8                 :: CONST,  SLIP,  VISC,     FAC1
      REAL*8                 :: FAC2,   FLUX,  AREA_CM2, RHB
      REAL*8                 :: RCM
      REAL*8                 :: TOT1,   TOT2
      INTEGER                :: KLOOP
      LOGICAL                :: STSBOX, NATBOX
     
      ! Specific to each class
      REAL*8                 :: RWET(2),CONST_V(2)
      REAL*8                 :: RHO(2),RATIO_R(2),REFF(2)
      REAL*8                 :: VTS(LLPAR,2),INMASS(4)
      REAL*8                 :: OUTMASS(4)
 
      ! Used for old Seinfeld & Pandis slip factor calc
      REAL*8                 :: sp_Lambda, sp_Num

      ! H2SO4 aqueous fraction
      INTEGER                :: LOCALSO4AQ
      
      ! Parameters
      REAL*8, PARAMETER      :: SALTDEN = 2200.0d0 ! density (kg/m3)
      REAL*8, PARAMETER      :: H2SO4DEN = 1840.0d0
      REAL*8, PARAMETER      :: BCDEN = 1000.d0 ! density (kg/m3)
      REAL*8, PARAMETER      :: R_STRAT_STS = 2.0d-7 ! Mean STS rad. (m)

      ! Indexing
      INTEGER, PARAMETER     :: IBC = 1
      INTEGER, PARAMETER     :: ISO4= 2
      INTEGER, PARAMETER     :: NSETTLE = 2
      INTEGER                :: IAERO
      LOGICAL                :: RUNCALC

      !=================================================================
      ! SETTLE_STRAT_AER begins here!
      !=================================================================

      ! Return if it's the start of the run
      IF ( GET_ELAPSED_SEC() == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,     L,    VTS,  P,        TEMP, RHB,  RWET ) 
!$OMP+PRIVATE( RATIO_R, RHO,   DP,   PDP,  CONST,    SLIP, VISC )
!$OMP+PRIVATE( DELZ,    DELZ1, TOT1, TOT2, AREA_CM2, FLUX, IAERO      )
!$OMP+PRIVATE( STSBOX, NATBOX, LOCALSO4AQ, SP_NUM, SP_LAMBDA, RUNCALC )
!$OMP+PRIVATE( CONST_V )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR       

         ! Initialize 
         DO L = 1, LLPAR
            VTS(L,:) = 0d0
         ENDDO

         ! Loop over levels
         DO L = 1, LLPAR

            ! Temperature [K]
            TEMP    = T(I,J,L)
 
            RUNCALC = ITS_IN_THE_STRATMESO(I,J,L)
            LOCALSO4AQ = 1.d0 - H2SO4_GASFRAC(I,J,L)

            IF (RUNCALC) THEN
               ! Need to translate for radii
               KLOOP = JLOP(I,J,L)

               !-----SDE DEBUG
               ! Replace these with actual values!
               ! Can calculate density from mass/vol
               RWET(IBC) = WERADIUS(KLOOP,2)*1.d-2
               RHO(IBC) = BCDEN
               !-----SDE DEBUG

               IF (LOCALSO4AQ.gt.TINY(0d0)) THEN
                  ! Sulfate now takes the form of STS
                  RWET(ISO4) = R_STRAT_STS ! Use simple mean radius for STS for now
                  RHO(ISO4) = H2SO4DEN ! kg/m3
                  NATBOX = (TEMP.le.T_STS)
               ELSE
                  ! Only sediment if aerosol form
                  RWET(ISO4) = -1d0
                  RHO(ISO4) = -1d0
               ENDIF
            ENDIF

            IF (.not.RUNCALC) THEN
               VTS(L,:) = 0d0
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
               ! Correct for partial gas content
               VTS(L,ISO4) = VTS(L,ISO4) * LOCALSO4AQ
            ENDIF ! RUNCALC
         ENDDO

         ! Handle model top condition
         L = LLPAR
         DELZ = BXHEIGHT(I,J,L)
         
         DO IAERO=1,2
            CONST_V(IAERO) = 1.d0/(1.d0+DTCHEM*VTS(L,IAERO)/DELZ)
         ENDDO

         STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4) * CONST_V(ISO4)
         STT(I,J,L,IDTBCPI)= STT(I,J,L,IDTBCPI)* CONST_V(IBC) 

         DO L = LLPAR-1,1,-1
            IF (ITS_IN_THE_TROP(I,J,L+1)) CYCLE
            DELZ  = BXHEIGHT(I,J,L)
            DELZ1 = BXHEIGHT(I,J,L+1)
            STT(I,J,L,IDTSO4) = 1.d0/(1.d0+DTCHEM*VTS(L,ISO4)/DELZ)
     &              * (STT(I,J,L,IDTSO4)+DTCHEM*VTS(L+1,ISO4)/DELZ1
     &              * STT(I,J,L+1,IDTSO4))
            STT(I,J,L,IDTBCPI) = 1.d0/(1.d0+DTCHEM*VTS(L,IBC)/DELZ)
     &              * (STT(I,J,L,IDTBCPI)+DTCHEM*VTS(L+1,IBC)/DELZ1
     &              * STT(I,J,L+1,IDTBCPI)) 
         ENDDO

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
         IF (.not.ALLOCATED(H2SO4_GASFRAC)) THEN 
            CALL INIT_GMISTRAT
         ENDIF
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PCENTER, TCENTER, INVAIR ) 
!$OMP+PRIVATE( H2SO4SUM, GF_DIFF )
!$OMP+PRIVATE( GF_PP, GF_INVT, GF_CFACTOR, GF_LOGPSULFATE, GF_PVAP )
!$OMP+SCHEDULE( DYNAMIC )

      DO J=1, JJPAR
      DO I=1, IIPAR
      DO L=1, LLPAR
         ! Only interested in low-pressure boxes
         PCENTER = GET_PCENTER(I,J,L)
         TCENTER = T(I,J,L)
         INVAIR  = AIRMW / AD(I,J,L)
         IF (PCENTER.ge.1.d2) THEN
            H2SO4_GASFRAC(I,J,L) = 0d0
            CYCLE
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
               H2SO4_GASFRAC(I,J,L) = 0.d0
            ELSEIF (GF_DIFF .lt. GF_RANGE) THEN
               H2SO4_GASFRAC(I,J,L) = GF_DIFF/GF_RANGE
            ELSE
               H2SO4_GASFRAC(I,J,L) = 1.d0
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

!$OMP END PARALLEL DO
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

      IF (.not.ALLOCATED(H2SO4_GASFRAC)) THEN 
         CALL CALC_H2SO4_GAS
      ENDIF

      SO4_PHOTFRAC = H2SO4_GASFRAC(I,J,L)

      END FUNCTION SO4_PHOTFRAC
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_psc
!
! !DESCRIPTION: Subroutine CALC\_PSC calculates PSC properties throughout the
!  chemistry grid using the thermodynamic parameterization described in
!  Kirner et al. (`Simulation of polar stratospheric clouds in the
!  chemistry-climate-model EMAC via the submodel PSC', Geosci. Mod. Dev.,
!  4, 169-182, 2011). This routine calculates solid PSC properties directly
!  before handing over to CALC\_STS for liquid aerosol properties.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_PSC
!
! !USES:
!
      USE DAO_MOD,          ONLY : T, AD, CLDF, QI, SPHU
      USE PRESSURE_MOD,     ONLY : GET_PCENTER
      USE LOGICAL_MOD,      ONLY : LPRT
      USE ERROR_MOD,        ONLY : DEBUG_MSG
      USE GRID_MOD,         ONLY : GET_YEDGE
      USE TRACER_MOD,       ONLY : STT, TRACER_MW_G
      USE TRACERID_MOD,     ONLY : IDTHNO3
      USE CHEMGRID_MOD,     ONLY : ITS_IN_THE_STRATMESO

      USE CMN_GCTM_MOD,     ONLY : AIRMW, Rd
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
      REAL*8,PARAMETER      :: PSC_PMIN   =   5.0d3

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

      ! Local conditions
      REAL*8                :: TCENTER, PCENTER, PCENTER_PA, DENAIR
      REAL*8                :: INVAIR, TINV

      ! Gridbox mixing ratios and partial pressures
      REAL*8                :: HNO3SUM, H2OSUM
      REAL*8                :: HNO3PP,  H2OPP
      REAL*8                :: PSATHNO3, PSATH2O

      ! Gridbox phase data
      REAL*8                :: HNO3_BOX_NAT, H2O_BOX_ICE

      ! Grid box location
      REAL*8                :: BOX_LAT_S, BOX_LAT_N
      LOGICAL               :: IS_VALID

      ! Loop variables
      INTEGER               :: I, J, L

      ! Loop limits for column
      INTEGER               :: MINLEV, MAXLEV

      !=================================================================
      ! CALC_PSC begins here!
      !=================================================================

      !-----SDE DEBUG
      ! Function stub - fill out later (SDE 04/13/13)
      WRITE(6,*) 'SDE DEBUG - FUNCTION NOT READY!'
      RETURN
      !-----SDE DEBUG

      ! Allocate if not ready
      IF (.not.ALLOCATED(RAD_PSC)) THEN
         CALL INIT_GMISTRAT
      ENDIF

      IF (LPRT) THEN
         CALL DEBUG_MSG( '### GMISTRAT: start CALC_PSC' )
      ENDIF

      ! First step is to determine what PSCs can form, if any
      ! Loop over latitude boxes first
      DO J=1,JJPAR
         BOX_LAT_S = GET_YEDGE(1,J,1)
         BOX_LAT_N = GET_YEDGE(1,J+1,1)
         IS_VALID = ((BOX_LAT_S.lt.PSC_MINLAT).or.
     &               (BOX_LAT_N.gt.PSC_MAXLAT))
         IF (.not.IS_VALID) THEN
            STATE_PSC(:,J,:) = 0
         ELSE
            ! Loop over longitudinal boxes
            DO I=1,IIPAR
               ! Set impossible loop limits
               MINLEV=LLPAR+1
               MAXLEV=0
               LLOOP: DO L=1,LLPAR

                  ! Get local conditions
                  PCENTER = GET_PCENTER(I,J,L) 
                  PCENTER_PA = PCENTER*1.d2
                  TINV = 1.d0/TCENTER
                  TCENTER = T(I,J,L)

                  ! Apply other limits from Kirner et al.
                  IS_VALID = ((.not.ITS_IN_THE_STRATMESO(I,J,L)).or.
     &                                  (PCENTER_PA.lt.PSC_PMIN).or.
     &                                  (PCENTER_PA.gt.PSC_PMAX))

                  IF (.not.IS_VALID) THEN
                     ! No PSCs
                     STATE_PSC(I,J,L) = 0
                     CYCLE LLOOP
                  ENDIF

                  IF (TCENTER.gt.T_MAX) THEN
                     ! STS only
                     STATE_PSC(I,J,L) = 1
                     CYCLE LLOOP
                  ENDIF
 
                  ! Calculate local air density
                  DENAIR = AVOGADRO * PCENTER_PA / (TCENTER*Rd)
                  INVAIR = AIRMW / AD(I,J,L)

                  !-----SDE DEBUG
                  ! Get available HNO3
                  ! Ignoring NIT for now
                  HNO3SUM = STT(I,J,L,IDTHNO3)*INVAIR/
     &                       TRACER_MW_G(IDTHNO3)
                  !-----SDE DEBUG

                  !-----SDE DEBUG
                  ! Revert to static H2O for now
                  H2OSUM = SPHU(I,J,L)*1.d-3*AD(I,J,L)
                  !-----SDE DEBUG

                  ! Calculate partial pressures (Pa)
                  HNO3PP = PCENTER_PA * HNO3SUM
                  H2OPP  = PCENTER_PA * H2OSUM

                  ! Calculate saturation pressures
                  PSATHNO3 = 10.d0**(((TSATHNO3_A+(TCENTER*TSATHNO3_B))*
     &            LOG10(H2OPP*760.0d0/101325.0d0)) + (TSATHNO3_C + 
     &            (TSATHNO3_D*TINV) + (TSATHNO3_E*TCENTER)))*
     &            101325.0d0/760.0d0

                  PSATH2O = 10.0d0**((TSATH2O_A*TINV)+TSATH2O_B)
                  
                  ! Use local ice mass ratio from GEOS-5 data
                  H2O_BOX_ICE = QI(I,J,L)*CLDF(L,I,J)*AIRMW/18.0d0

                  ! Ice exists/possible?
                  IF ((H2OPP.gt.PSATH2O).or.
     &                   (H2O_BOX_ICE.gt.TINY(0d0))) THEN
                     STATE_PSC(I,J,L) = 3
                  ELSE
                     ! If ice not possible could still have NAT
                     IF (STATE_PSC(I,J,L).gt.1) THEN
                        ! Had PSCs here last time
                        IF (((HNO3PP.gt.PSATHNO3).and.
     &                       (PSATHNO3.gt.PSATH2O)).or.
     &                       (PCENTER.lt.P_MAXNAT)) THEN
                           STATE_PSC(I,J,L) = 2
                        ELSE
                           ! STS only
                           STATE_PSC(I,J,L) = 1
                        ENDIF 
                     ELSE
                        ! STS only
                        STATE_PSC(I,J,L) = 1
                     ENDIF
                  ENDIF
               ENDDO LLOOP ! L
                  ! Check for NAT formation
            ENDDO ! I
         ENDIF ! IS_POLAR
      ENDDO ! J

      END SUBROUTINE CALC_PSC
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
      SUBROUTINE INIT_GMISTRAT
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

      !=================================================================
      ! INIT_GMISTRAT begins here!
      !=================================================================

      WRITE( 6,'(a)') REPEAT( '=', 79 )
      WRITE( 6,'(a)') 'U N I F I E D   C H E M I S T R Y'
      WRITE( 6,'(a)') 'Routines written by SEBASTIAN D. EASTHAM'
      WRITE( 6,'(a)') 'Last Modification Date: 05/09/13'
      WRITE( 6,'(a)') REPEAT( '=', 79 )

      ! Initialize H2SO4 gas fraction
      ALLOCATE( H2SO4_GASFRAC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'H2SO4_GAS_FRAC' )
      H2SO4_GASFRAC = 0d0

      ! Determine which tracers are to be emitted
      N_EMIT = 0
      IDT_EMIT = 0

      IF (LCFCEMIS) THEN
         ! R-113/114/115
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFC
         TRAC_EMIT(N_EMIT) = 'GENERIC_CFC'

         ! HCFC-22/141b/142b
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTHCFC
         TRAC_EMIT(N_EMIT) = 'GENERIC_HCFC'

         ! R-11
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFC11 
         TRAC_EMIT(N_EMIT) = 'CFC11'

         ! R-12
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCFC12 
         TRAC_EMIT(N_EMIT) = 'CFC12'

         ! H-1202
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH1202 
         TRAC_EMIT(N_EMIT) = 'CBR2F2'

         ! H-1211
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH1211  
         TRAC_EMIT(N_EMIT) = 'H1211'

         ! H-1301
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH1301 
         TRAC_EMIT(N_EMIT) = 'H1301'

         ! H-2402
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTH2402  
         TRAC_EMIT(N_EMIT) = 'H2402_YEARLY'
      ENDIF

      IF (LCLEMIS) THEN
         ! CCl4
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCCl4
         TRAC_EMIT(N_EMIT) = 'CCL4'

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
      ALLOCATE( RAD_PSC( IIPAR, JJPAR, LLPAR, N_PSC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RAD_PSC' )
      RAD_PSC = 0d0
      
      ALLOCATE( KG_PSC( IIPAR, JJPAR, LLPAR, N_PSC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'KG_PSC' )
      KG_PSC = 0d0
      
      ALLOCATE( SAD_PSC( IIPAR, JJPAR, LLPAR, N_PSC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SAD_PSC' )
      SAD_PSC = 0d0
      
      ALLOCATE( OD_PSC( IIPAR, JJPAR, LLPAR, N_PSC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OD_PSC' )
      OD_PSC = 0d0
      
      ALLOCATE( NDENS_PSC( IIPAR, JJPAR, LLPAR, N_PSC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NDENS_PSC' )
      NDENS_PSC = 0d0
      
      ALLOCATE( RHO_PSC( IIPAR, JJPAR, LLPAR, N_PSC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RHO_PSC' )
      RHO_PSC = 0d0
      
      ALLOCATE( STATE_PSC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'STATE_PSC' )
      STATE_PSC = 0d0
      
      !-----SDE DEBUG
      WRITE(6,'(a)') 'WARNING: Need to add subroutine to read in PSCs!'
      !-----SDE DEBUG

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
      IF ( ALLOCATED( RAD_PSC    ) ) DEALLOCATE( RAD_PSC    )
      IF ( ALLOCATED( KG_PSC     ) ) DEALLOCATE( KG_PSC     )
      IF ( ALLOCATED( OD_PSC     ) ) DEALLOCATE( OD_PSC     )
      IF ( ALLOCATED( STATE_PSC  ) ) DEALLOCATE( STATE_PSC  )
      IF ( ALLOCATED( NDENS_PSC  ) ) DEALLOCATE( NDENS_PSC  )

      END SUBROUTINE CLEANUP_GMISTRAT
!
!EOC
      END MODULE GMISTRAT_MOD
