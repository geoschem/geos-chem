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
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: SET_INITIAL_MIXRATIOS
      PUBLIC  :: EMISS_BASIC
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
      !
      !=================================================================

      INTEGER, PARAMETER            :: AER_NLEVS=51

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
      REAL*8,DIMENSION(:,:),ALLOCATABLE :: AER_MR
      REAL*8,DIMENSION(:),ALLOCATABLE   :: AER_PLEVS
      INTEGER,DIMENSION(MAX_EMIT)       :: IDT_EMIT
      REAL*8,DIMENSION(:,:),ALLOCATABLE :: GRID_EMIT
 
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
!     (7 ) CFCl3               = CFC-11               < LSETCFC
!     (8 ) CF2Cl2              = CFC-12               < LSETCFC
!     (9 ) CF2Br2              = Halon 1202           < LSETCFC
!     (10) CF2ClBr             = Halon 1211           < LSETCFC
!     (11) CBrF3               = Halon 1301           < LSETCFC
!     (12) C2Br2F4             = Halon 2402           < LSETCFC
!     (13) Cl2                 = Cl2                  < LSETCL
!     (14) ClOx                = Cl + ClO             < LSETCL
!     (15) CCl4                = CCl4                 < LSETCL
!     (16) CH3Cl               = CH3Cl                < LSETCL
!     (17) CH3CCl3             = CH3CCl3              < LSETCL
!     (18) HCl                 = HCl                  < LSETCL
!     (19) HOCl                = HOCl                 < LSETCL
!     (20) Cl2O2               = Cl2O2                < LSETCL
!     (21) ClNO2               = ClNO2                < LSETCL
!     (22) ClONO2 (ClNO3)      = ClNO3                < LSETCL
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
            TRAC_IDX = IDTCF2Br2
            TRAC_SET  = 'CBR2F2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCBrClF2
            TRAC_SET  = 'CBRCLF2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCBrF3
            TRAC_SET  = 'CBRF3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTC2Br2F4
            TRAC_SET  = 'C2BR2F4'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            ! CFC  11, 12, and [113/114/115]
            TRAC_IDX = IDTCFCl3
            TRAC_SET  = 'CFCL3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)   
    
            TRAC_IDX = IDTCF2Cl2
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
         IDT_EMIT( N_EMIT) = IDTCFCl3
         TRAC_EMIT(N_EMIT) = 'CFC11'

         ! R-12
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCF2Cl2
         TRAC_EMIT(N_EMIT) = 'CFC12'

         ! H-1202
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCF2Br2
         TRAC_EMIT(N_EMIT) = 'CBR2F2'

         ! H-1211
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCBrClF2
         TRAC_EMIT(N_EMIT) = 'H1211'

         ! H-1301
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCBrF3
         TRAC_EMIT(N_EMIT) = 'H1301'

         ! H-2402
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTC2Br2F4
         TRAC_EMIT(N_EMIT) = 'H2402_YEARLY'
      ENDIF

      IF (LCLEMIS) THEN
         ! CCl4
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCCl4
         TRAC_EMIT(N_EMIT) = 'CCL4'

         ! CH3CCl3
         N_EMIT = N_EMIT + 1
         IDT_EMIT( N_EMIT) = IDTCCl4
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
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_GMISTRAT begins here!
      !=================================================================

      IF ( ALLOCATED( GRID_EMIT  ) ) DEALLOCATE( GRID_EMIT  )

      END SUBROUTINE CLEANUP_GMISTRAT
!
!EOC
      END MODULE GMISTRAT_MOD
