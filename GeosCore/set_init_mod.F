!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: set_init_mod
!     
! !DESCRIPTION: Module SET\_INIT\_MOD contains routines and variables which
!  can be used to set initial mixing ratios for specified tracers.
!\\   
!\\   
! !INTERFACE: 
!
      MODULE SET_INIT_MOD
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
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: SET_INITIAL_MIXRATIOS
!
! PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: APPLY_2DTRAC 
!
! !REMARKS:
!  Global concentrations of each of the following tracers are read from
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
!     (9 ) CF2ClBr             = Halon 1211           < LSETCFC
!     (10) CBrF3               = Halon 1301           < LSETCFC
!     (11) C2Br2F4             = Halon 2402           < LSETCFC
!     (12) Cl2                 = Cl2                  < LSETCL
!     (13) ClOx                = Cl + ClO             < LSETCL
!     (14) CCl4                = CCl4                 < LSETCL
!     (15) CH3Cl               = CH3Cl                < LSETCL
!     (16) CH3CCl3             = CH3CCl3              < LSETCL
!     (17) HCl                 = HCl                  < LSETCL
!     (18) HOCl                = HOCl                 < LSETCL
!     (19) Cl2O2               = Cl2O2                < LSETCL
!     (20) ClNO2               = ClNO2                < LSETCL
!     (21) ClONO2 (ClNO3)      = ClNO3                < LSETCL
!     (22) OClO                = OClO                 < LSETCL
!     (23) Br2                 = Br2                  < LSETBR
!     (24) Br                  = Br                   < LSETBR
!     (25) BrO                 = BrO                  < LSETBR
!     (26) HOBr                = HOBr                 < LSETBR
!     (27) HBr                 = HBr                  < LSETBR
!     (28) BrNO2               = BrNO2                < LSETBR
!     (29) BrNO3               = BrNO3                < LSETBR
!     (30) CHBr3               = CHBr3                < LSETBR
!     (31) CH2Br2              = CH2Br2               < LSETBR
!     (32) CH3Br               = CH3Br                < LSETBR
!                                                                             .
!
! !REVISION HISTORY: 
!  26 Mar 2013 - S. D. Eastham - Initial version
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
      !
      ! Arrays
      !
      ! AER_MR             : AER 2D mixing ratios (v/v)
      ! AER_PLEVS          : Pressure levels of AER data (hPa)
      !
      ! Strings
      !
      ! TRAC_FILE          : Input filename
      ! MONTH_STR          : Current month
      ! TRAC_NAME          : Tracer name to be overwritten
      !
      ! Logicals
      ! 
      ! TRAC_ADD            : Add to (rather than overwrite) tracer
      ! 
      !=================================================================

      ! Scalars
      INTEGER                           :: TRAC_IDX

      ! Arrays
      REAL*8,DIMENSION(:,:),ALLOCATABLE :: AER_MR
      REAL*8,DIMENSION(:),ALLOCATABLE   :: AER_PLEVS
 
      ! Strings
      CHARACTER(LEN=255)                :: TRAC_FILE
      CHARACTER(LEN=2)                  :: MONTH_STR
      CHARACTER(LEN=255)                :: TRAC_NAME

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
! !DESCRIPTION: Subroutine SET\_INITIAL\_MIXRATIOS is the public interface 
!  of the module, fielding the call to initialize the 2D mixing ratios
!  of whichever tracers have been designated.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_INITIAL_MIXRATIOS
!
! !USES:
!
      USE LOGICAL_MOD,          ONLY : LSETCH4,  LSETOCS,  LSETCFC
      USE LOGICAL_MOD,          ONLY : LSETCL,   LSETBR,   LSETN2O
      USE LOGICAL_MOD,          ONLY : LSETH2SO4,LGMISTRAT
      USE TRACERID_MOD
      USE TRACER_MOD,           ONLY : N_TRACERS, STT
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
      INTEGER            :: AS, N, J
      REAL*8             :: PLAST, PCURR, PNEXT
      REAL*8             :: YLAT
      CHARACTER(LEN=255) :: MSG
      LOGICAL            :: USE2DDATA
      REAL*8             :: C3090S, C0030S, C0030N, C3090N
      INTEGER            :: CH4_YEAR

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

      ! Retrieval from GC data fields not yet implemented
      WRITE(6,*) '<<<<2D data read still necessary for N2O>>>>'
      USE2DDATA = (USE2DDATA.or.LSETN2O)

#if !defined(GRID4x5) && !defined(GRID2x25)
      IF (LGMISTRAT.and.USE2DDATA) THEN
       MSG = 'Zonal means preprocessed only for 2x2.5 and 4x5 grid.'
       CALL ERROR_STOP( MSG, 'SET_INITIAL_MIXRATIOS (set_init_mod.f)!' )
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

         WRITE(6,*) 'Setting initial values for:'
            WRITE(6,*) '|---------|---------|---------|'

         ! Bromine species (except Halons)
         IF (LSETBR) THEN
            WRITE(6,*) '|Br2      |Br       |BrO      |'
            WRITE(6,*) '|HOBr     |HBr      |BrNO2    |'
            WRITE(6,*) '|BrONO2   |CHBr3    |CH2Br2   |'
            WRITE(6,*) '|CH3Br    |         |         |'
            WRITE(6,*) '|---------|---------|---------|'
            TRAC_IDX = IDTBr2
            TRAC_NAME = 'BR2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         
            TRAC_IDX = IDTBr
            TRAC_NAME = 'BR'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         
            TRAC_IDX = IDTBrO
            TRAC_NAME = 'BRO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         
            TRAC_IDX = IDTHOBr
            TRAC_NAME = 'HOBR'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         
            TRAC_IDX = IDTHBr
            TRAC_NAME = 'HBR'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTBrNO2
            TRAC_NAME = 'BRNO2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTBrNO3
            TRAC_NAME = 'BRNO3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCHBr3 
            TRAC_NAME = 'CHBR3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCH2Br2
            TRAC_NAME = 'CH2BR2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCH3Br
            TRAC_NAME = 'CH3BR'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         ENDIF

         ! Chlorine species (except CFCS, HCFCs and Halons)
         IF (LSETCL) THEN
            WRITE(6,*) '|Cl2      |ClOx     |Cl2O2    |'
            WRITE(6,*) '|ClNO2    |ClONO2   |HCl      |'
            WRITE(6,*) '|HOCl     |OClO     |CH3Cl    |'
            WRITE(6,*) '|CH3CCl3  |CCl4     |         |'
            WRITE(6,*) '|---------|---------|---------|'
            TRAC_IDX = IDTClOx
            TRAC_NAME = 'CL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
            TRAC_NAME = 'CLO'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
     
            TRAC_IDX = IDTCl2
            TRAC_NAME = 'CL2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCl2O2
            TRAC_NAME = 'CL2O2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTClNO2
            TRAC_NAME = 'CLNO2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)

            TRAC_IDX = IDTClNO3
            TRAC_NAME = 'CLNO3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)

            TRAC_IDX = IDTHCl
            TRAC_NAME = 'HCL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTHOCl
            TRAC_NAME = 'HOCL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTOClO
            TRAC_NAME = 'OCLO'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)

            ! Handle chloroalkanes      
            TRAC_IDX = IDTCH3Cl
            TRAC_NAME = 'CH3CL'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCH3CCl3
            TRAC_NAME = 'CH3CCL3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCCl4
            TRAC_NAME = 'CCL4'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         ENDIf
         
         ! CFCs, HCFCs and Halons
         IF (LSETCFC) THEN
            WRITE(6,*) '|H-1201   |H-1311   |H-2402   |'
            WRITE(6,*) '|CFC-11   |CFC-12   |CFC-X    |'
            WRITE(6,*) '|HCFC-X   |         |         |'
            WRITE(6,*) '|---------|---------|---------|'
            ! H1201, 1311 and 2402
            TRAC_IDX = IDTCBrClF2
            TRAC_NAME = 'CBRCLF2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTCBrF3
            TRAC_NAME = 'CBRF3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            TRAC_IDX = IDTC2Br2F4
            TRAC_NAME = 'C2BR2F4'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
     
            ! CFC  11, 12, and [113/114/115]
            TRAC_IDX = IDTCFCl3
            TRAC_NAME = 'CFCL3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)   
    
            TRAC_IDX = IDTCF2Cl2
            TRAC_NAME = 'CF2CL2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
   
            TRAC_IDX = IDTCFC
            TRAC_NAME = 'C2CL3F3'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.) 
            TRAC_NAME = 'C2CL2F4'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
            TRAC_NAME = 'C2CLF5'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
     
            ! HCFC 22, 141b and 142b
            TRAC_IDX = IDTHCFC
            TRAC_NAME = 'CHCLF2'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
            TRAC_NAME = 'C2H3FCL2'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
            TRAC_NAME = 'C2H3F2CL'
            CALL APPLY_2DTRAC(.TRUE.,.FALSE.)
         ENDIF

         IF (LSETOCS) THEN
            WRITE(6,*) '|OCS      |         |         |'
            WRITE(6,*) '|---------|---------|---------|'
            TRAC_IDX = IDTOCS
            TRAC_NAME = 'OCS'
            CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
         ENDIF

         IF (LSETH2SO4) THEN
            WRITE(6,*) '|H2SO4    |         |         |'
            WRITE(6,*) '|---------|---------|---------|'
            TRAC_IDX = IDTSO4
            ! Actually H2SO4 aerosol
            TRAC_NAME = 'SO4'
            CALL APPLY_2DTRAC(.FALSE.,.TRUE.)
         ENDIF
      ENDIF ! USE2DDATA

      ! N2O and CH4 are treated differently
      ! Take these from the existing GEOS-Chem 3D input data,
      ! rather than from 2D data

      IF (LSETN2O) THEN
         WRITE(6,*) '|N2O      |         |         |'
         WRITE(6,*) '|---------|---------|---------|'

         ! SDE DEBUG
         WRITE(6,*) 'SDE DEBUG: N2O STILL BEING SET IN 2D!'

         IF ( LPRT ) CALL DEBUG_MSG( '### SET_INIT: Forcing 2D N2O' )

         TRAC_IDX = IDTN2O
         TRAC_NAME = 'N2O'
         CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
      ENDIF

      IF (LSETCH4) THEN
         WRITE(6,*) '|CH4      |         |         |'
         WRITE(6,*) '|---------|---------|---------|'

         ! Loop over boxes per grid block
         DO J=1,JJPAR

            ! Grid-box latitude index
            YLAT  = GET_YMID( 1,J,1 )
            
            ! Pick the CH4 concentration [ppbv] for the proper lat bin
            ! CH4 values are read in "chemdr.f"
            IF ( LFUTURE ) THEN
               CH4_YEAR = GET_FUTURE_YEAR()
            ELSE
               CH4_YEAR = GET_YEAR()
            ENDIF

            ! Get CH4 [ppbv] in 4 latitude bins for each year
            CALL GET_GLOBAL_CH4( CH4_YEAR, .TRUE., C3090S,
     &                           C0030S,   C0030N, C3090N, .FALSE. )

            IF ( YLAT < -30d0 ) THEN
               STT(:,J,:,IDTCH4) = C3090S*1.d-9
            ELSE IF ( YLAT >= -30d0 .and. YLAT < 0d0  ) THEN
               STT(:,J,:,IDTCH4) = C0030S*1.d-9
            ELSE IF ( YLAT >=   0d0 .and. YLAT < 30d0 ) THEN
               STT(:,J,:,IDTCH4) = C0030N*1.d-9
            ELSE
               STT(:,J,:,IDTCH4) = C3090N*1.d-9
            ENDIF

!            ! Convert from [ppbv CH4] to [molec CH4/cm3]
!            CBLK(KLOOP,ICH4) = CBLK(KLOOP,ICH4) *1d-9 * AIRDENS(JLOOP)
         ENDDO

         ! Code to set using AER 2D data
!         TRAC_IDX = IDTCH4
!         TRAC_NAME = 'CH4'
!         CALL APPLY_2DTRAC(.FALSE.,.FALSE.)
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
      USE TRACER_MOD,    ONLY : STT ! v/v
      USE CMN_SIZE_MOD,  ONLY : IIPAR, JJPAR, LLPAR
      USE CHEMGRID_MOD,  ONLY : ITS_IN_THE_STRAT
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
      REAL*8  :: TRACMR
      INTEGER :: VERTCOUNT
      LOGICAL :: FOUNDLEV,EXTRAP
      REAL*8  :: TRACMULT
      LOGICAL :: SET_CELL
      INTEGER :: IU_FILE

      !=================================================================
      ! APPLY_2DTRAC begins here!
      !=================================================================

      TRAC_FILE = TRIM(TRAC_FILE_ROOT)//TRIM(TRAC_NAME)
      TRAC_FILE = TRIM(TRAC_FILE)//'_'//TRIM(MONTH_STR)//'.dat'

      ! Get a free LUN
      IU_FILE = findFreeLUN()

      IOS = 1
      OPEN( IU_FILE, FILE=TRIM( TRAC_FILE ),STATUS='OLD',IOSTAT=IOS) 
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) 'SET_INIT_MOD: Could not read ', TRIM( TRAC_FILE )
         CALL IOERROR( IOS, IU_FILE,'SET_INIT_MOD:APPLY_2DTRAC')
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
                  ! Only want to set strat cells
                  SET_CELL = ITS_IN_THE_STRAT(I,J,L)
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
               ! Send to STT (currently v/v)
               IF (SET_CELL) THEN
                  IF (TRAC_ADD) THEN
                     STT(I,J,L,TRAC_IDX) = STT(I,J,L,TRAC_IDX)+TRACMR
                  ELSE
                     STT(I,J,L,TRAC_IDX) = TRACMR
                  ENDIF
                ENDIF
            ENDDO
         ENDDO
      ENDDO
 
      ! Return to calling program
      END SUBROUTINE APPLY_2DTRAC
!EOC

      END MODULE SET_INIT_MOD
