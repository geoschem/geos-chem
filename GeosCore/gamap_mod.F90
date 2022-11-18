#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gamap_mod.F90
!
! !DESCRIPTION: Module GAMAP\_MOD contains routines to create GAMAP
!  "tracerinfo.dat" and "diaginfo.dat" files which are customized to each
!   particular GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
MODULE GAMAP_MOD
!
! !USES:
!
  USE CMN_DIAG_MOD        ! Diagnostic parameters
  USE inquireMod,    ONLY : findFreeLUN
  USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: DO_GAMAP
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: CREATE_DINFO
  PRIVATE :: CREATE_TINFO
  PRIVATE :: WRITE_TINFO
  PRIVATE :: WRITE_SEPARATOR
  PRIVATE :: INIT_DIAGINFO
  PRIVATE :: INIT_TRACERINFO
  PRIVATE :: INIT_GAMAP
  PRIVATE :: CLEANUP_GAMAP
!
! !REMARKS:
!  For more information, please see the GAMAP Online Users' Manual:
!   http://acmg.seas.harvard.edu/gamap/doc/index.html
!
! !REVISION HISTORY:
!  03 May 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! For "diaginfo.dat"
  INTEGER,           PARAMETER   :: MAXCAT    = 160
  INTEGER,           PARAMETER   :: SPACING   = 1000
  INTEGER                        :: NCATS
  INTEGER,           ALLOCATABLE :: OFFSET(:)
  CHARACTER(LEN=40), ALLOCATABLE :: CATEGORY(:)
  CHARACTER(LEN=40), ALLOCATABLE :: DESCRIPT(:)
  CHARACTER(LEN=255)             :: DFILE

  ! For "tracerinfo.dat"
  INTEGER,           PARAMETER   :: MAXDIAG   = MAX_DIAG
  INTEGER,           PARAMETER   :: MAXTRACER = 2 * MAX_TRACER
  INTEGER,           ALLOCATABLE :: NTRAC(:)
  INTEGER,           ALLOCATABLE :: INDEX(:,:)
  REAL*4,            ALLOCATABLE :: MWT(:,:)
  REAL*4,            ALLOCATABLE :: SCALE(:,:)
  CHARACTER(LEN=40), ALLOCATABLE :: NAME(:,:)
  CHARACTER(LEN=40), ALLOCATABLE :: FNAME(:,:)
  CHARACTER(LEN=40), ALLOCATABLE :: UNIT(:,:)
  CHARACTER(LEN=255)             :: TFILE

  ! Other variables
  CHARACTER(LEN=16)              :: STAMP
  CHARACTER(LEN=40)              :: SIM_NAME

  ! Species ID flags
  INTEGER :: id_NK1

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_gamap
!
! !DESCRIPTION: Subroutine DO\_GAMAP is the driver program for creating
!  the customized GAMAP files "diaginfo.dat" and "tracerinfo.dat".
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_GAMAP( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE TIME_MOD,       ONLY : SYSTEM_TIMESTAMP
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  03 May 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GC_SUCCESS
    
    ! Exit if this is a dry-run
    IF ( Input_Opt%DryRun ) RETURN

    ! Allocate and initialize variables
    CALL INIT_GAMAP( Input_Opt, State_Chm, RC )

    ! Create a timestamp with the system date & time
    STAMP    = SYSTEM_TIMESTAMP()

    ! Get simulation name
    SIM_NAME = Input_Opt%SimulationName

    ! Write "diaginfo.dat" file
    CALL CREATE_DINFO()

    ! Write "tracerinfo.dat" file
    CALL CREATE_TINFO( Input_Opt, State_Chm, RC )

    ! Deallocate variables
    CALL CLEANUP_GAMAP

  END SUBROUTINE DO_GAMAP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: create_dinfo
!
! !DESCRIPTION: Subroutine CREATE\_DINFO writes information about diagnostic
!  categories to a customized "diaginfo.dat" file. (bmy, 5/3/05)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CREATE_DINFO()
!
! !USES:
!
    USE FILE_MOD, ONLY : IOERROR
!
! !REVISION HISTORY:
!  03 May 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: IOS, N
    INTEGER :: IU_FILE

    !=================================================================
    ! CREATE_DINFO begins here!
    !=================================================================

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open "diaginfo.dat" file for output
    OPEN( IU_FILE, FILE=TRIM( DFILE ), STATUS='UNKNOWN', IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'create_dinfo:1' )

    ! Write file header
    WRITE( IU_FILE, '(a)' ) '#' // REPEAT( '=', 78 )
    WRITE( IU_FILE,  100  ) STAMP
    WRITE( IU_FILE,  105  ) TRIM( SIM_NAME )
    WRITE( IU_FILE,  110  )
    WRITE( IU_FILE,  115  )
    WRITE( IU_FILE, '(a)' ) '# File Format:'
    WRITE( IU_FILE, '(a)' ) '# ' // REPEAT( '-', 77 )
    WRITE( IU_FILE,  120  )
    WRITE( IU_FILE,  125  )
    WRITE( IU_FILE,  130  )
    WRITE( IU_FILE,  135  )
    WRITE( IU_FILE,  125  )
    WRITE( IU_FILE, '(a)' ) '#'
    WRITE( IU_FILE,  140  ) SPACING
    WRITE( IU_FILE, '(a)' ) '#' // REPEAT( '=', 78 )

    ! Loop over categories
    DO N = 1, NCATS

       ! Write one line to "diaginfo.dat" file
       WRITE( IU_FILE, 145, IOSTAT=IOS ) &
            OFFSET(N), ADJUSTL( CATEGORY(N) ), ADJUSTL( DESCRIPT(N) )

       ! Error check
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'create_dinfo:1' )
    ENDDO

    ! FORMAT strings
100 FORMAT( '# diaginfo.dat: Created by GEOS-CHEM at ', a,     /,'#' )
105 FORMAT( '# ****** CUSTOMIZED FOR ', a, ' SIMULATION *****',/,'#' )
110 FORMAT( '# This file contains category names and the offsets',   &
            ' which they are stored'                                 )
115 FORMAT( '# in file "tracerinfo.dat".  This file is read into',   &
            ' GAMAP by routine', /, '# "ctm_diaginfo.pro".',   /,'#' )
120 FORMAT( '# OFFSET    (I8 )  Constant to add to tracer numbers',  &
            ' in order to distinguish', /, '#', 18x, 'for the given',&
            ' diagnostic category, as stored in file', /, '#', 18x,  &
            '"tracerinfo.dat".  OFFSET may be up to 8 digits long.'  )
125 FORMAT( '#  --       (1X )  1-character spacer'                  )
130 FORMAT( '# CATEGORY  (A40)  Category name for CTM diagnostics.', &
            ' NOTE: The category name', /, '#', 18x,                 &
            'can be up to 40 chars long, but historically the',      &
            ' GEOS-CHEM', /,'#', 18x, 'and GISS models have used an',&
            ' 8-character category name.'                            )
135 FORMAT( '# COMMENT   (A  )  Descriptive comment string',   /,'#' )
140 FORMAT('##### SPACING BETWEEN DIAGNOSTIC CATEGORY OFFSETS = ',i8 )
145 FORMAT( i8, 1x, a40, 1x, a )

    ! Close file
    CLOSE( IU_FILE )

  END SUBROUTINE CREATE_DINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: create_tinfo
!
! !DESCRIPTION: Subroutine CREATE\_TINFO writes information about tracers to
!  a customized tracerinfo.dat" file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CREATE_TINFO( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE FILE_MOD,           ONLY : IOERROR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Apr 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: D, IOS, N, T
    INTEGER           :: IU_FILE
    REAL*4            :: SCALE_NEW
    CHARACTER(LEN=2)  :: C
    CHARACTER(LEN=40) :: UNIT_NEW, NAME_NEW, FNAME_NEW

    !=================================================================
    ! CREATE_TINFO begins here!
    !=================================================================

    ! Assume success
    RC      =  GC_SUCCESS

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open "tracerinfo.dat" file for output
    OPEN( IU_FILE, FILE=TRIM( TFILE ), STATUS='UNKNOWN', IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'create_tinfo:1' )

    ! Write file header
    WRITE( IU_FILE, '(a)' ) '#' // REPEAT( '=', 78 )
    WRITE( IU_FILE,  100  ) STAMP
    WRITE( IU_FILE,  105  ) TRIM( SIM_NAME )
    WRITE( IU_FILE,  110  )
    WRITE( IU_FILE,  115  )
    WRITE( IU_FILE, '(a)' ) '# File Format:'
    WRITE( IU_FILE, '(a)' ) '# ' // REPEAT( '-', 77 )
    WRITE( IU_FILE,  120  )
    WRITE( IU_FILE,  125  )
    WRITE( IU_FILE,  130  )
    WRITE( IU_FILE,  135  )
    WRITE( IU_FILE,  145  )
    WRITE( IU_FILE,  150  )
    WRITE( IU_FILE,  125  )
    WRITE( IU_FILE,  155  )

    ! FORMAT strings
100 FORMAT( '# tracerinfo.dat: Created by GEOS-CHEM at ', a,   /,'#' )
105 FORMAT( '# ****** CUSTOMIZED FOR ', a, ' SIMULATION *****',/,'#' )
110 FORMAT( '# This file contains name weight and index',            &
            ' information about GEOS-CHEM'                           )
115 FORMAT( '# tracers.  It is read by routine ',                    &
            '"ctm_tracerinfo.pro" of the GAMAP package.',      /,'#' )
120 FORMAT( '# NAME     (A8   )  Tracer name (up to 8 chars)'        )
125 FORMAT( '#  --      (1X   )  1-character spacer'                 )
130 FORMAT( '# FULLNAME (A30  )  Full tracer name (up to 30 chars)'  )
135 FORMAT( '# MOLWT    (E10.0)  Molecular weight (kg/mole)'         )
145 FORMAT( '# TRACER   (I9   )  Tracer number (up to 9 digits)'     )
150 FORMAT( '# SCALE    (E10.3)  Standard scale factor to convert',  &
            ' to unit given below'                                   )
155 FORMAT( '# UNIT     (A40  )  Unit string',                 /,'#' )

    !-------------------------------------
    ! 0: Tracers [ppbv]
    !-------------------------------------

    ! Write separator line
    CALL WRITE_SEPARATOR( IU_FILE, 0 )

    ! Loop over tracers
    DO T = 1, NTRAC(45)

       ! GAMAP tracer number
       N = ( SPACING * 0 ) + T

       ! Write tracers [ppbv] to "tracerinfo.dat" file
       CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), &
                         SCALE(T,45), UNIT(T,45), N )
    ENDDO

    !-------------------------------------
    ! SPACING*1: Tracers [molec/cm2/s]
    !-------------------------------------

    ! Write separator line
    CALL WRITE_SEPARATOR( IU_FILE, 100 )

    ! Loop over tracers
    DO T = 1, NTRAC(45)

       ! GAMAP tracer number
       N = ( SPACING * 1 ) + T

       ! New scale
       SCALE_NEW = 1.0e0

       ! New unit
       IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
          UNIT_NEW = 'atoms C/cm2/s'
       ELSE
          UNIT_NEW = 'molec/cm2/s'
       ENDIF

       ! Write tracers [molec/cm2/s] to "tracerinfo.dat"
       CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), &
                         SCALE_NEW,   UNIT_NEW, N )
    ENDDO

    !-------------------------------------
    ! SPACING*2: Tracers [molec/cm2]
    !-------------------------------------

    ! Write separator line
    CALL WRITE_SEPARATOR( IU_FILE, 200 )

    ! Loop over tracers
    DO T = 1, NTRAC(45)

       ! GAMAP tracer number
       N = ( SPACING * 2 ) + T

       ! New scale
       SCALE_NEW = 1.0e0

       ! New unit
       IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
          UNIT_NEW = 'atoms C/cm2'
       ELSE
          UNIT_NEW = 'molec/cm2'
       ENDIF

       ! Write tracers [molec/cm2] to "tracerinfo.dat"
       CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), &
                         SCALE_NEW,   UNIT_NEW, N )
    ENDDO

    !-------------------------------------
    ! SPACING*3: Tracers [kg/s]
    !-------------------------------------

    ! Write separator line
    CALL WRITE_SEPARATOR( IU_FILE, 300 )

    ! Loop over tracers
    DO T = 1, NTRAC(45)

       ! GAMAP tracer number
       N = ( SPACING * 3 ) + T

       ! New scale
       SCALE_NEW = 1.0e0

       ! New unit
       IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
          UNIT_NEW = 'kg C/s'
       ELSE
          UNIT_NEW = 'kg/s'
       ENDIF

       ! Write tracers [kg/s] to "tracerinfo.dat"
       CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), &
                         SCALE_NEW,   UNIT_NEW, N )
    ENDDO

    !-------------------------------------
    ! SPACING*4: Tracers [kg]
    !-------------------------------------

    ! Write separator line
    CALL WRITE_SEPARATOR( IU_FILE, 400 )

    ! Loop over tracers
    DO T = 1, NTRAC(45)

       ! GAMAP tracer number
       N = ( SPACING * 4 ) + T

       ! New scale
       SCALE_NEW = 1.0e0

       ! New unit
       IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
          UNIT_NEW = 'kg C'
       ELSE
          UNIT_NEW = 'kg'
       ENDIF

       ! Write tracers [kg] to "tracerinfo.dat"
       CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), &
                         SCALE_NEW,   UNIT_NEW, N )
    ENDDO

    !-------------------------------------
    ! All other diagnostics
    !-------------------------------------
    DO D = 1, MAXDIAG

       ! If tracers are defined then...
       IF ( NTRAC(D) > 0 ) THEN

          ! Skip ND45, we already wrote tracers above
          IF ( D == 45 ) CYCLE

          ! Write separator
          CALL WRITE_SEPARATOR( IU_FILE, D )

          ! Write tracers to file
          DO T = 1, NTRAC(D)
             CALL WRITE_TINFO( IU_FILE, NAME(T,D), FNAME(T,D), MWT(T,D), &
                               SCALE(T,D), UNIT(T,D), INDEX(T,D) )
          ENDDO
       ENDIF
    ENDDO

  END SUBROUTINE CREATE_TINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_tinfo
!
! !DESCRIPTION: Subroutine WRITE\_TINFO writes one line to the customized
!  "tracerinfo.dat" file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WRITE_TINFO( IU_FILE, NAME, FNAME, MWT, SCALE, UNIT, N )
!
! !USES:
!
    USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: IU_FILE ! Logical unit number
    CHARACTER(LEN=*), INTENT(IN) :: NAME    ! GAMAP short tracer name
    CHARACTER(LEN=*), INTENT(IN) :: FNAME   ! GAMAP long tracer name
    REAL*4,           INTENT(IN) :: MWT     ! Molecular weight [kg/mole]
    INTEGER,          INTENT(IN) :: N       ! Tracer number
    REAL*4,           INTENT(IN) :: SCALE   ! GAMAP scale factor
    CHARACTER(LEN=*), INTENT(IN) :: UNIT    ! Unit string
!
! !REVISION HISTORY:
!  03 May 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: IOS

    !=================================================================
    ! WRITE_TINFO begins here!
    !=================================================================

    ! Write one line to "tracerinfo.dat" file
    WRITE( IU_FILE, 100, IOSTAT=IOS ) &
         ADJUSTL( NAME ), ADJUSTL( FNAME ), MWT, 1, N, SCALE, TRIM( UNIT )

    ! Error check
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_tinfo:1' )

    ! FORMAT string
100 FORMAT( a8, 1x, a30, es10.3, i3, i9, es10.3, 1x, a )

  END SUBROUTINE WRITE_TINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_separator

!
! !DESCRIPTION: Subroutine WRITE\_SEPARATOR writes a separator block to the
!  customized "tracerinfo.dat" file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WRITE_SEPARATOR( IU_FILE, DIAG )
!
! !USES:
!
    USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: IU_FILE
    INTEGER, INTENT(IN) :: DIAG   ! GEOS-Chem diagnostic number
!
! !REVISION HISTORY:
!  03 May 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: IOS
    CHARACTER(LEN=79) :: SEPARATOR

    !=================================================================
    ! WRITE_SEPARATOR begins here!
    !=================================================================

    ! Create separator string
    SEPARATOR = '#' // REPEAT( '=', 78 )

    ! Write separator string
    WRITE( IU_FILE, '(a)', IOSTAT=IOS ) SEPARATOR
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_separator:1' )

    ! Write the appropriate message
    SELECT CASE( DIAG )
    CASE( 0 )
       WRITE( IU_FILE, 100, IOSTAT=IOS )
    CASE( 100 )
       WRITE( IU_FILE, 110, IOSTAT=IOS )
    CASE( 200 )
       WRITE( IU_FILE, 120, IOSTAT=IOS )
    CASE( 300 )
       WRITE( IU_FILE, 130, IOSTAT=IOS )
    CASE( 400 )
       WRITE( IU_FILE, 140, IOSTAT=IOS )
    CASE( 500 )
       WRITE( IU_FILE, 150, IOSTAT=IOS )
    CASE( 600 )
       WRITE( IU_FILE, 160, IOSTAT=IOS )
    CASE( 700 )          !(win, 7/9/09)
       WRITE( IU_FILE, 190, IOSTAT=IOS ) !sfarina - this can stay
    CASE( 999 )
       WRITE( IU_FILE, 170, IOSTAT=IOS )
    CASE DEFAULT
       WRITE( IU_FILE, 180, IOSTAT=IOS ) DIAG
    END SELECT

    ! Error check
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_separator:2' )

    ! Write separator string
    WRITE( IU_FILE, '(a)', IOSTAT=IOS ) SEPARATOR
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_separator:3' )

    ! FORMAT strings
100 FORMAT( '# GEOS-CHEM tracers [ppbv]'             )
110 FORMAT( '# GEOS-CHEM tracers [molec/cm2/s]'      )
120 FORMAT( '# GEOS-CHEM tracers [molec/cm2]'        )
130 FORMAT( '# GEOS-CHEM tracers [kg/s]'             )
140 FORMAT( '# GEOS-CHEM tracers [kg]'               )
150 FORMAT( '# GEOS-CHEM tracers [ug/m3]'            )
160 FORMAT( '# SOA GPROD & APROD [kg/kg]'            )
170 FORMAT( '# Quantities for Soil NOx restart file' )
180 FORMAT( '# ND', i2.2, ' diagnostic quantities'   )
190 FORMAT( '# TOMAS aerosol rate [cm-3 s-1]'        )   !(win, 7/9/09)

  END SUBROUTINE WRITE_SEPARATOR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diaginfo
!
! !DESCRIPTION: Subroutine INIT\_DIAGINFO initializes the CATEGORY, DESCRIPT,
!  and OFFSET variables, which are used to define the "diaginfo.dat" file
!  for GAMAP.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_DIAGINFO()
!
! !REVISION HISTORY:
!  17 Oct 1996 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    !=================================================================
    ! INIT_DIAGINFO begins here!
    !=================================================================

    N           = 1
    CATEGORY(N) = 'IJ-AVG-$'
    DESCRIPT(N) = 'Tracer concentration'
    OFFSET(N)   = SPACING * 0

    N           = N + 1
    CATEGORY(N) = 'IJ-24H-$'
    DESCRIPT(N) = '24-hr avg tracer conc.'
    OFFSET(N)   = SPACING * 0

    N           = N + 1
    CATEGORY(N) = 'INST-MAP'
    DESCRIPT(N) = 'Instantaneous tracer'
    OFFSET(N)   = SPACING * 0

    N           = N + 1
    CATEGORY(N) = 'DUSTSRCE'
    DESCRIPT(N) = 'Dust emission'
    OFFSET(N)   = SPACING * 4

#ifdef TOMAS
    N           = N + 1               !(win, 7/9/09)
    CATEGORY(N) = 'NK-EMISS'
    DESCRIPT(N) = 'Size-res aerosol number emissions'
    OFFSET(N)   = SPACING * 4

    N           = N + 1               !(win, 7/9/09)
    CATEGORY(N) = 'SF-EMISS'
    DESCRIPT(N) = 'Size-res sulfate mass emissions'
    OFFSET(N)   = SPACING * 4

    N           = N + 1               !(win, 7/9/09)
    CATEGORY(N) = 'SS-EMISS'
    DESCRIPT(N) = 'Size-res sea-salt mass emissions'
    OFFSET(N)   = SPACING * 4

    N           = N + 1               !(win, 7/9/09)
    CATEGORY(N) = 'ECIL-SRC'
    DESCRIPT(N) = 'Size-res H-phillic EC mass emissions'
    OFFSET(N)   = SPACING * 4

    N           = N + 1               !(win, 7/9/09)
    CATEGORY(N) = 'ECOB-SRC'
    DESCRIPT(N) = 'Size-res H-phobic EC mass emissions'
    OFFSET(N)   = SPACING * 4

    N           = N + 1               !(win, 7/9/09)
    CATEGORY(N) = 'OCIL-SRC'
    DESCRIPT(N) = 'Size-res H-phillic OC mass emissions'
    OFFSET(N)   = SPACING * 4

    N           = N + 1               !(win, 7/9/09)
    CATEGORY(N) = 'OCOB-SRC'
    DESCRIPT(N) = 'Size-res H-phobic OC mass emissions'
    OFFSET(N)   = SPACING * 4

    N           = N + 1               !(win, 7/9/09)
    CATEGORY(N) = 'DUST-SRC'
    DESCRIPT(N) = 'Size-res dust mass emissions'
    OFFSET(N)   = SPACING * 4

    N           = N + 1
    CATEGORY(N) = 'TMS-COND'
    DESCRIPT(N) = 'Condensation rate'
    OFFSET(N)   = SPACING * 6

    N           = N + 1
    CATEGORY(N) = 'TMS-COAG'
    DESCRIPT(N) = 'Condensation rate'
    OFFSET(N)   = SPACING * 6

    N           = N + 1
    CATEGORY(N) = 'TMS-NUCL'
    DESCRIPT(N) = 'Condensation rate'
    OFFSET(N)   = SPACING * 6

    N           = N + 1
    CATEGORY(N) = 'TMS-AQOX'
    DESCRIPT(N) = 'Condensation rate'
    OFFSET(N)   = SPACING * 6

    N           = N + 1
    CATEGORY(N) = 'AERO-FIX'
    DESCRIPT(N) = 'Condensation rate'
    OFFSET(N)   = SPACING * 6

    N           = N + 1
    CATEGORY(N) = 'TMS-SOA'
    DESCRIPT(N) = 'Condensation rate'
    OFFSET(N)   = SPACING * 6

    N           = N + 1           !(win, 7/9/09)
    CATEGORY(N) = 'TOMAS-3D'
    DESCRIPT(N) = '3-D aerosol rate of a selected microphysics'
    OFFSET(N)   = SPACING * 7
#endif

    N           = N + 1
    CATEGORY(N) = 'PEDGE-$'
    DESCRIPT(N) = 'Pressure at level edges'
    OFFSET(N)   = SPACING * 10 

    N           = N + 1
    CATEGORY(N) = 'DAO-FLDS'
    DESCRIPT(N) = 'GMAO 2-D met fields'
    OFFSET(N)   = SPACING * 11

    N           = N + 1
    CATEGORY(N) = 'DAO-3D-$'
    DESCRIPT(N) = 'GMAO 3-D met fields'
    OFFSET(N)   = SPACING * 12

    N           = N + 1
    CATEGORY(N) = 'OD-MAP-$'
    DESCRIPT(N) = 'Optical Depths'
    OFFSET(N)   = SPACING * 14

    N           = N + 1
    CATEGORY(N) = 'CHEM-L=$'
    DESCRIPT(N) = 'Chemical Prod/Loss'
    OFFSET(N)   = SPACING * 16

    N           = N + 1
    CATEGORY(N) = 'TIME-SER'
    DESCRIPT(N) = 'Timeseries quantities'
    OFFSET(N)   = SPACING * 19

    N           = N + 1
    CATEGORY(N) = 'BXHGHT-$'
    DESCRIPT(N) = 'Boxheight, airmass, etc'
    OFFSET(N)   = SPACING * 24

    N           = N + 1
    CATEGORY(N) = 'PBLDEPTH'
    DESCRIPT(N) = 'Afternoon PBL height'
    OFFSET(N)   = SPACING * 27

    N           = N + 1
    CATEGORY(N) = 'HG-SRCE'
    DESCRIPT(N) = 'Hg emissions'
    OFFSET(N)   = SPACING * 34

    N           = N + 1
    CATEGORY(N) = 'HG0-ANTH'
    DESCRIPT(N) = 'Hg(0) Anthro Emissions'
    OFFSET(N)   = SPACING * 61

    N           = N + 1
    CATEGORY(N) = 'HG0-AQUA'
    DESCRIPT(N) = 'Hg(0) Ocean Mass'
    OFFSET(N)   = SPACING * 62

    N           = N + 1
    CATEGORY(N) = 'HG0-RECY'
    DESCRIPT(N) = 'Hg(0) Land Prompt Recycling'
    OFFSET(N)   = SPACING * 63

    N           = N + 1
    CATEGORY(N) = 'HGNET-OC'
    DESCRIPT(N) = 'Hg(0) Net Ocean Emissions'
    OFFSET(N)   = SPACING * 64

    N           = N + 1
    CATEGORY(N) = 'HG0-GEOG'
    DESCRIPT(N) = 'Hg(0) Geogenic Emissions'
    OFFSET(N)   = SPACING * 65

    N           = N + 1
    CATEGORY(N) = 'HG2-ANTH'
    DESCRIPT(N) = 'Hg(II) Anthro Emissions'
    OFFSET(N)   = SPACING * 66

    N           = N + 1
    CATEGORY(N) = 'HG2-AQUA'
    DESCRIPT(N) = 'Hg(II) Ocean Mass'
    OFFSET(N)   = SPACING * 67

    N           = N + 1
    CATEGORY(N) = 'HG2-SINK'
    DESCRIPT(N) = 'Hg(II) Sinking'
    OFFSET(N)   = SPACING * 68

    N           = N + 1
    CATEGORY(N) = 'HGP-ANTH'
    DESCRIPT(N) = 'Hg(P) Anthro Emissions'
    OFFSET(N)   = SPACING * 69

    N           = N + 1
    CATEGORY(N) = 'HGT-AQUA'
    DESCRIPT(N) = 'Hg Total Ocean Mass'
    OFFSET(N)   = SPACING * 70

    N           = N + 1
    CATEGORY(N) = 'HGP-AQUA'
    DESCRIPT(N) = 'Hg(P) Ocean Mass'
    OFFSET(N)   = SPACING * 71

    N           = N + 1
    CATEGORY(N) = 'HGP-CONV'
    DESCRIPT(N) = 'Hg(II) Aqeous Conversion to Particle'
    OFFSET(N)   = SPACING * 72

    N           = N + 1
    CATEGORY(N) = 'HG0-BURN'
    DESCRIPT(N) = 'Hg(0) Biomass Burning Emissions'
    OFFSET(N)   = SPACING * 73

    N           = N + 1
    CATEGORY(N) = 'HG0-VEGT'
    DESCRIPT(N) = 'Hg(0) Veg Transp. Emissions'
    OFFSET(N)   = SPACING * 74

    N           = N + 1
    CATEGORY(N) = 'HG0-SOIL'
    DESCRIPT(N) = 'Hg(0) Soil Emissions'
    OFFSET(N)   = SPACING * 75

    N           = N + 1
    CATEGORY(N) = 'HG0-FXUP'
    DESCRIPT(N) = 'Hg(0) Gross Ocean Up Flux'
    OFFSET(N)   = SPACING * 76

    N           = N + 1
    CATEGORY(N) = 'HG0-FXDN'
    DESCRIPT(N) = 'Hg(0) Gross Ocean Down Flux'
    OFFSET(N)   = SPACING * 77

    N           = N + 1
    CATEGORY(N) = 'HG0-SNOW'
    DESCRIPT(N) = 'Hg(0) Snow Emissions'
    OFFSET(N)   = SPACING * 78

    N           = N + 1
    CATEGORY(N) = 'HG-NETOX'
    DESCRIPT(N) = 'Hg(0) to Hg(II) Net Oxidation'
    OFFSET(N)   = SPACING * 79

    N           = N + 1
    CATEGORY(N) = 'HG2-OXOH'
    DESCRIPT(N) = 'Hg(0) Oxidation to Hg(II) by OH'
    OFFSET(N)   = SPACING * 80

    N           = N + 1
    CATEGORY(N) = 'HG2-OXO3'
    DESCRIPT(N) = 'Hg(0) Oxidation to Hg(II) by O3'
    OFFSET(N)   = SPACING * 81

    N           = N + 1
    CATEGORY(N) = 'HG2-SALT'
    DESCRIPT(N) = 'Hg(II) Loss by Sea Salt'
    OFFSET(N)   = SPACING * 82

    N           = N + 1
    CATEGORY(N) = 'HG2-SSRX'
    DESCRIPT(N) = 'Hg(0) Sea Salt Reaction Rate'
    OFFSET(N)   = SPACING * 83

    N           = N + 1
    CATEGORY(N) = 'HG2-OXBR'
    DESCRIPT(N) = 'Hg(0) Oxidation to Hg(II) by Br'
    OFFSET(N)   = SPACING * 84

    N           = N + 1
    CATEGORY(N) = 'PL-HG2-$'
    DESCRIPT(N) = 'Prod / loss of Hg2'
    OFFSET(N)   = SPACING * 35

    N           = N + 1
    CATEGORY(N) = 'DRYD-FLX'
    DESCRIPT(N) = 'Drydep fluxes'
    OFFSET(N)   = SPACING * 36

    N           = N + 1
    CATEGORY(N) = 'DRYD-VEL'
    DESCRIPT(N) = 'Drydep velocities'
    OFFSET(N)   = SPACING * 37

    N           = N + 1
    CATEGORY(N) = 'CO2-SRCE'
    DESCRIPT(N) = 'CO2 fluxes'
    OFFSET(N)   = SPACING * 40

    N           = N + 1
    CATEGORY(N) = 'OCEAN-HG'
    DESCRIPT(N) = 'Oceanic Hg emissions'
    OFFSET(N)   = SPACING * 41

    N           = N + 1
    CATEGORY(N) = 'LFLASH-$'
    DESCRIPT(N) = 'Lightning flash rates'
    OFFSET(N)   = SPACING * 42

    N           = N + 1
    CATEGORY(N) = 'IJ-SOA-$'
    DESCRIPT(N) = 'SOA concentrations'
    OFFSET(N)   = SPACING * 43

    ! For mercury simulation only so we can use same spacing. (ccc, 5/21/10)
    N           = N + 1
    CATEGORY(N) = 'SNOW-HG'
    DESCRIPT(N) = 'Reducible Hg mass in ocean snow and ice'
    OFFSET(N)   = SPACING * 49

    ! Non-reducible snow Hg reservoir
    N           = N + 1
    CATEGORY(N) = 'SNOW-HGN'
    DESCRIPT(N) = 'Non-reducible Hg mass in ocean snow and ice'
    OFFSET(N)   = SPACING * 50

    N           = N + 1
    CATEGORY(N) = 'SNHG-LN'
    DESCRIPT(N) = 'Reducible Hg mass in land snow and ice'
    OFFSET(N)   = SPACING * 51

    N           = N + 1
    CATEGORY(N) = 'SNHGN-LN'
    DESCRIPT(N) = 'Non-reducible Hg mass in snow and ice'
    OFFSET(N)   = SPACING * 52

    !(FP 06/19/2009)
    N           = N + 1
    CATEGORY(N) = 'THETA-$'
    DESCRIPT(N) = 'Potential temperature'
    OFFSET(N)   = SPACING * 57

    N           = N + 1
    CATEGORY(N) = 'SHIP-$$$'
    DESCRIPT(N) = 'Ship diagnostics'
    OFFSET(N)   = SPACING * 59

    N           = N + 1
    CATEGORY(N) = 'RST-SOIL'
    DESCRIPT(N) = 'Soil NOx restart file'
    OFFSET(N)   = SPACING * 60

    ! Added for POPs (clf, 3/2/2011)
    N           = N + 1
    CATEGORY(N) = 'PG-SRCE'
    DESCRIPT(N) = 'POPs Emissions'
    OFFSET(N)   = SPACING * 52

    N           = N + 1
    CATEGORY(N) = 'PG-PP-$'
    DESCRIPT(N) = 'POP partitioning/oxidation'
    OFFSET(N)   = SPACING * 53

    N           = N + 1
    CATEGORY(N) = 'RADMAP-$'
    DESCRIPT(N) = 'Radiative diagnostics'
    OFFSET(N)   = SPACING * 72

    ! Number of categories
    NCATS = N

  END SUBROUTINE INIT_DIAGINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tracerinfo
!
! !DESCRIPTION: Subroutine INIT\_TRACERINFO initializes the NAME, FNAME,
!  MWT, INDEX, UNIT arrays which are used to define the
!  "tracerinfo.dat" file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_TRACERINFO( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE CMN_FJX_MOD,        ONLY : W_
    USE CMN_SIZE_MOD,       ONLY : NRHAER, NDUST, NSTRATAER
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
#ifdef TOMAS
    USE TOMAS_MOD,          ONLY : IBINS, ICOMP,   IDIAG  !(win, 7/14/09)
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  17 Oct 1996 - R. Yantosca & P. Le Sager - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: N, NN, NYMDb, NHMSb, T, T0, T1
    INTEGER                 :: nAdvect, nDryDep, nSpecies
    LOGICAL                 :: DO_TIMESERIES
    CHARACTER(LEN=40)       :: NAME5, NAME6
    REAL(fp),      SAVE     :: POP_XMW
    REAL(fp)                :: DUM
    INTEGER                 :: N_Hg_CATS
    INTEGER                 :: SPC_INDEX

    ! Objects
    TYPE(Species), POINTER  :: SpcInfo

    !=================================================================
    ! INIT_TRACERINFO begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    SpcInfo => NULL()

    ! Initialization
    DUM = 1.0

    ! Set a flag if any timeseries diagnostics are turned on
    DO_TIMESERIES = ( Input_Opt%DO_ND51 .or. Input_Opt%DO_ND51b )

    ! Number of total species
    nSpecies  = State_Chm%nSpecies

    ! Number of advected species
    nAdvect   = State_Chm%nAdvect

    ! Number of drydep species
    nDryDep   = State_Chm%nDryDep

    !----------------------------------
    ! General tracer information
    !----------------------------------

    ! NOTE: Write metadata for all GEOS-Chem species to tracerinfo.dat
    ! so that GAMAP will be able to read the new netCDF restart file.
    ! (bmy, 1/23/17)
    DO T = 1, nSpecies

       ! Set pointer to this species information
       ! NOTE: requires 1:1 tracer index to species index mapping
       SpcInfo => State_Chm%SpcData(T)%Info

       ! Store quantities for each tracer
       NAME (T,45) = SpcInfo%Name
       FNAME(T,45) = TRIM( NAME(T,45) ) // ' tracer'
       SCALE(T,45) = 1.0e+9
       INDEX(T,45) = T   !changed from N (phs, 3/19/03)
       MWT  (T,45) = SpcInfo%MW_g * 1.e-3_fp
       UNIT(T,45) = 'ppbv'

       ! Special handling for Rn-Pb-Be simulation (bmy, 5/11/05)
       IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN
          SELECT CASE( T )
          CASE( 1 )
             UNIT (T,45) = 'mBq/SCM'
             SCALE(T,45) = 5.6397E+22
          CASE( 2, 3 )
             UNIT (T,45) = 'mBq/SCM'
             SCALE(T,45) = 2.6141E+19
          CASE( 4, 5, 6, 7 )
             UNIT (T,45) = 'mBq/SCM'
             SCALE(T,45) = 4.0513E+21
          CASE DEFAULT
             ! Nothing
          END SELECT
       ENDIF

       ! For mercury, print as pptv (bmy, 1/24/06)
       IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
          UNIT(T,45)  = 'pptv'
          SCALE(T,45) = 1.0e+12
       ENDIF

    ENDDO

    ! Number of ND45 tracers, now all species (bmy, 1/23/17)
    NTRAC(45) = nSpecies


    IF ( DO_TIMESERIES ) THEN

       !-------------------------------------      
       ! Optical depths
       !-------------------------------------

       ! Number of tracers
       NTRAC(21) = 8 + (NRHAER+NDUST)*3 + (NRHAER*2) + (NSTRATAER*3)

       ! Loop over tracers
       DO T = 1, NTRAC(21)

          ! Define quantities
          UNIT (T,21) = 'unitless'
          INDEX(T,21) = T + ( SPACING * 14 )
          MWT  (T,21) = 0e0
          SCALE(T,21) = 1e0

          ! Get name long-name (and sometimes, unit)
          SELECT CASE( T )
          CASE( 1  )
             NAME (T,21) = 'OPTD'
             FNAME(T,21) = 'Cloud optical depth'
          CASE( 2  )
             ! GEOS-3, GEOS-4: CLDTOT
             NAME (T,21) = 'CLDTOT'
             FNAME(T,21) = '3-D cloud frc'
          CASE( 3  )
             NAME (T,21) = 'OBSOLETE1'
             FNAME(T,21) = 'This diagnostic is obsolete'
          CASE( 4  )
             NAME (T,21) = 'OPD'
             FNAME(T,21) = 'Mineral dust opt depth'
          CASE( 5  )
             NAME (T,21) = 'SD'
             FNAME(T,21) = 'Mineral dust surface area'
             UNIT (T,21) = 'cm2/cm3'
          CASE( 6  )
             NAME (T,21) = 'OPSO4'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Sulfate opt depth (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 7  )
             NAME (T,21) = 'HGSO4'
             FNAME(T,21) = 'Hygr growth of SO4'
          CASE( 8  )
             NAME (T,21) = 'SSO4'
             FNAME(T,21) = 'Sulfate surface area'
             UNIT (T,21) = 'cm2/cm3'
          CASE( 9  )
             NAME (T,21) = 'OPBC'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Black carbon opt depth (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 10 )
             NAME (T,21) = 'HGBC'
             FNAME(T,21) = 'Hygr growth of BC'
          CASE( 11 )
             NAME (T,21) = 'SBC'
             FNAME(T,21) = 'BC surface area'
             UNIT (T,21) = 'cm2/cm3'
          CASE( 12 )
             NAME (T,21) = 'OPOC'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Org carbon opt depth (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 13 )
             NAME (T,21) = 'HGOC'
             FNAME(T,21) = 'Hygr growth of OC'
          CASE( 14 )
             NAME (T,21) = 'SOC'
             FNAME(T,21) = 'OC surface area'
             UNIT (T,21) = 'cm2/cm3'
          CASE( 15 )
             NAME (T,21) = 'OPSSa'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Seasalt (accum) opt depth (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 16 ) 
             NAME (T,21) = 'HGSSa'
             FNAME(T,21) = 'Hygr growth of seasalt (accum)'
          CASE( 17 )
             NAME (T,21) = 'SSSa'
             FNAME(T,21) = 'Seasalt (accum) surface area'
             UNIT (T,21) = 'cm2/cm3'
          CASE( 18 )
             NAME (T,21) = 'OPSSc'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Seasalt (coarse) opt depth (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 19 ) 
             NAME (T,21) = 'HGSSc'
             FNAME(T,21) = 'Hygr growth of seasalt (coarse)'
          CASE( 20 )
             NAME (T,21) = 'SSSc'
             FNAME(T,21) = 'Seasalt (coarse) surface area'
             UNIT (T,21) = 'cm2/cm3'
          CASE( 21 )
             NAME (T,21) = 'OPD1'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Dust bin 1 AOD (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 22 )
             NAME (T,21) = 'OPD2'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Dust bin 2 AOD (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 23 )
             NAME (T,21) = 'OPD3'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Dust bin 3 AOD (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 24 )
             NAME (T,21) = 'OPD4'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Dust bin 4 AOD (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 25 )
             NAME (T,21) = 'OPD5'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Dust bin 5 AOD (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 26 )
             NAME (T,21) = 'OPD6'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Dust bin 6 AOD (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 27 )
             NAME (T,21) = 'OPD7'//Input_Opt%STRWVSELECT(1)
             FNAME(T,21) = 'Dust bin 7 AOD (' &
                           // Input_Opt%STRWVSELECT(1) //'nm)'
          CASE( 28 )
             NAME (T,21) = 'OPSO4'//Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'Seasalt (accum) opt depth (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 29 )
             NAME (T,21) = 'OPBC' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'BC opt depth (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 30 )
             NAME (T,21) = 'OPOC' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'OC opt depth (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 31 )
             NAME (T,21) = 'OPSSa' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'Seasalt (accum) opt depth (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 32 )
             NAME (T,21) = 'OPSSc' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'seasalt (coarse) opt depth (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 33 )
             NAME (T,21) = 'OPD1' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'Dust bin 1 AOD (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 34 )
             NAME (T,21) = 'OPD2' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'Dust bin 2 AOD (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 35 )
             NAME (T,21) = 'OPD3' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'Dust bin 3 AOD (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 36 )
             NAME (T,21) = 'OPD4' // Input_Opt%STRWVSELECT(2) 
             FNAME(T,21) = 'Dust bin 4 AOD (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 37 )
             NAME (T,21) = 'OPD5' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'Dust bin 5 AOD (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 38 )
             NAME (T,21) = 'OPD6' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'Dust bin 6 AOD (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 39 )
             NAME (T,21) = 'OPD7' // Input_Opt%STRWVSELECT(2)
             FNAME(T,21) = 'Dust bin 7 AOD (' &
                            // Input_Opt%STRWVSELECT(2) // 'nm)'
          CASE( 40 )
             NAME (T,21) = 'OPSO4' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'Seasalt (accum) opt depth (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 41 )
             NAME (T,21) = 'OPBC' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'BC opt depth (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 42 )
             NAME (T,21) = 'OPOC' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'OC opt depth (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 43 )
             NAME (T,21) = 'OPSSa' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'Seasalt (accum) opt depth (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 44 )
             NAME (T,21) = 'OPSSc' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'seasalt (coarse) opt depth (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 45 )
             NAME (T,21) = 'OPD1' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'Dust bin 1 AOD (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 46 )
             NAME (T,21) = 'OPD2' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'Dust bin 2 AOD (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 47 )
             NAME (T,21) = 'OPD3' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'Dust bin 3 AOD (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 48 )
             NAME (T,21) = 'OPD4' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'Dust bin 4 AOD (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 49 )
             NAME (T,21) = 'OPD5' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'Dust bin 5 AOD (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 50 )
             NAME (T,21) = 'OPD6' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'Dust bin 6 AOD (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'
          CASE( 51 )
             NAME (T,21) = 'OPD7' // Input_Opt%STRWVSELECT(3)
             FNAME(T,21) = 'Dust bin 7 AOD (' &
                            // Input_Opt%STRWVSELECT(3) // 'nm)'

          ! For the UCX simualtion
          CASE( 52 )
             NAME (T,21) = 'ODSLA'
             FNAME(T,21) = 'SLA AOD (600 nm)'
          CASE( 53 )
             NAME (T,21) = 'SASLA'
             FNAME(T,21) = 'SLA surface area'
             UNIT (T,21) = 'cm2/cm3'
          CASE( 54 )
             NAME (T,21) = 'NDSLA'
             FNAME(T,21) = 'SLA number density'
             UNIT (T,21) = 'num/cm3'
          CASE( 55 )
             NAME (T,21) = 'ODSPA'
             FNAME(T,21) = 'PSC type 1a/2 AOD (600 nm)'
          CASE( 56 )
             NAME (T,21) = 'SASPA'
             FNAME(T,21) = 'PSC type 1a/2 surface area'
             UNIT (T,21) = 'cm2/cm3'
          CASE( 57 )
             NAME (T,21) = 'NDSPA'
             FNAME(T,21) = 'SPA number density'
             UNIT (T,21) = 'num/cm3'
          CASE( 58 )
             NAME (T,21) = 'ISOPAOD'
             FNAME(T,21) = 'Isoprene AOD (550 nm)'
          CASE( 59 )
             NAME (T,21) = 'AQAVOL'
             FNAME(T,21) = 'Aq. aerosol vol (cm3/cm3)'
          CASE( 60 )
             NAME (T,21) = 'OBSOLETE2'
             FNAME(T,21) = 'This diagnostic is obsolete'
          CASE DEFAULT
             ! Nothing
          END SELECT

       ENDDO

       !-------------------------------------      
       ! Surface pressure (ND31)
       !-------------------------------------   
       T           = 1
       NTRAC(31)   = T
       NAME (T,31) = 'PSURF'
       FNAME(T,31) = 'Surface pressure'
       UNIT (T,31) = 'hPa'
       INDEX(T,31) = T + ( SPACING * 10 )
       MWT  (T,31) = 0e0
       SCALE(T,31) = 1e0

       !-------------------------------------      
       ! Afternoon-average boundary 
       ! layer heights (ND41) + timeseries
       !-------------------------------------      

       ! Number of tracers
       NTRAC(41) = 2

       ! Loop over tracers
       DO T = 1, NTRAC(41)

          ! Get name and unit for each met field
          SELECT CASE( T )
          CASE( 1 )
             NAME(T,41) = 'PBL-M'
             UNIT(T,41) = 'm'
          CASE( 2 )
             NAME(T,41) = 'PBL-L'
             UNIT(T,41) = 'levels'
          CASE DEFAULT
             ! Nothing
          END SELECT

          ! Define the rest of the quantities
          FNAME(T,41) = 'PBL depth'
          INDEX(T,41) = T + ( SPACING * 27 )
          MWT  (T,41) = 0e0
          SCALE(T,41) = 1e0

       ENDDO

       !-------------------------------------      
       ! Chemically-produced OH, etc (ND43)
       !-------------------------------------      

       ! Number of tracers
       NTRAC(43) = 4

       ! Loop over tracers
       DO T = 1, NTRAC(43)

          ! Get name and unit for each met field
          SELECT CASE( T )
          CASE( 1 )
             NAME (T,43) = 'OH'
             UNIT (T,43) = 'molec/cm3'
             MWT  (T,43) = 17e-3
             INDEX(T,43) = 1 + ( SPACING * 16 )
          CASE( 2 )
             NAME(T,43)  = 'HO2'
             UNIT(T,43)  = 'v/v'
             MWT (T,43)  = 33e-3
             INDEX(T,43) = 3 + ( SPACING * 16 )
          CASE( 3 )
             NAME(T,43)  = 'O1D'
             UNIT(T,43)  = 'molec/cm3'
             MWT (T,43)  = 16e-3
             INDEX(T,43) = 4 + ( SPACING * 16 )
          CASE( 4 )
             NAME(T,43)  = 'O'
             UNIT(T,43)  = 'molec/cm3'
             MWT (T,43)  = 16e-3
             INDEX(T,43) = 5 + ( SPACING * 16 )   
          CASE DEFAULT
             ! Nothing
          END SELECT

          ! Define the rest of the quantities
          FNAME(T,43) = 'Chemically produced ' // TRIM( NAME(T,43) )
          SCALE(T,43) = 1e0
       ENDDO

    ENDIF ! DO_TIMESERIES

#ifdef TOMAS
    !-------------------------------------
    ! Dry deposition fluxes     (ND44)
    ! Dry deposition velocities (ND44)
    !-------------------------------------
    IF ( ND44 > 0 ) THEN

       ! Loop over # of dry depositing species
       DO N = 1, nDryDep + (ICOMP-IDIAG)* IBINS

          IF ( N <= nDryDep ) THEN

             !-----------------------------------------------
             ! Drydep flux
             !-----------------------------------------------
             T             = State_Chm%Map_DryDep(N)
             SpcInfo       => State_Chm%SpcData(T)%Info
             NAME (N,44)   = TRIM( SpcInfo%Name ) // 'df'
             FNAME(N,44)   = TRIM( SpcInfo%Name ) // ' drydep flux'
             MWT  (N,44)   = SpcInfo%MW_g * 1.e-3_fp
             SCALE(N,44)   = 1.0e0
             INDEX(N,44)   = T + ( SPACING * 36 )
             NTRAC(44)     = NTRAC(44) + 1

             ! For the Rn simulation, unit is kg/s (bmy, 2/22/08)
             IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN
                UNIT(N,44) = 'kg/s'
             ELSE
                UNIT(N,44) = 'molec/cm2/s'
             ENDIF

             !-----------------------------------------------
             ! Drydep velocity
             !-----------------------------------------------
             NN            = N + nDryDep  +(ICOMP-IDIAG)*IBINS
             NAME (NN,44)  = TRIM( SpcInfo%Name ) // 'dv'
             FNAME(NN,44)  = TRIM( SpcInfo%Name ) // ' drydep velocity'
             MWT  (NN,44)  = SpcInfo%MW_g * 1.e-3_fp
             UNIT (NN,44)  = 'cm/s'
             SCALE(NN,44)  = 1.0e0
             INDEX(NN,44)  = T + ( SPACING * 37 )
             NTRAC(44)     = NTRAC(44) + 1

          ELSE

             !-----------------------------------------------
             ! Drydep flux: TOMAS aerosol tracers
             !-----------------------------------------------
             T = id_NK1 + IBINS - 1 + ( N - nDryDep )
             SpcInfo       => State_Chm%SpcData(T)%Info

             ! Drydep flux (extended deposited species)
             NAME (N,44)  = TRIM( SpcInfo%Name ) // 'df'
             FNAME(N,44)  = TRIM( SpcInfo%Name ) // ' drydep flux'
             UNIT (N,44)  = 'molec/cm2/s'
             MWT  (N,44)  = SpcInfo%MW_g * 1.e-3_fp
             SCALE(N,44)  = 1.0e0
             INDEX(N,44)  = T + ( SPACING * 36 )
             NTRAC(44)    = NTRAC(44) + 1

             !-----------------------------------------------
             ! Drydep velocity all other G-C tracers
             !-----------------------------------------------
             NN            = N + nDryDep + (ICOMP-IDIAG)* IBINS
             NAME (NN,44)  = TRIM( SpcInfo%Name ) // 'dv'
             FNAME(NN,44)  = TRIM( SpcInfo%Name ) // ' drydep velocity'
             MWT  (NN,44)  = SpcInfo%MW_g * 1.e-3_fp
             UNIT (NN,44)  = 'cm/s'
             SCALE(NN,44)  = 1.0e0
             INDEX(NN,44)  = T + ( SPACING * 37 )
             NTRAC(44)     = NTRAC(44) + 1

          ENDIF

       ENDDO

    ENDIF
#endif

    !-------------------------------------
    ! Timeseries diagnostics
    !-------------------------------------
    IF ( DO_TIMESERIES ) THEN

       ! Number of tracers
       ! Increased from 26 to 32 (mpb,2009)
       NTRAC(48) = 17

       ! Loop over tracers
       DO T = 1, NTRAC(48)

          ! Define quantities
          INDEX(T,48) = T + ( SPACING * 19 )
          MWT  (T,48) = 0e0
          SCALE(T,48) = 1e0

          ! Get name, long-name (and others where necessary)
          SELECT CASE( T )
          CASE( 1  )
             NAME (T,48) = 'OH'
             FNAME(T,48) = 'OH time series'
             UNIT (T,48) = 'molec/cm3'
          CASE( 2  )
             NAME (T,48) = 'NOy'
             FNAME(T,48) = 'NOy time series'
             UNIT (T,48) = 'ppbv'
             SCALE(T,48) = 1.0e+9
          CASE( 3  )
             NAME (T,48) = 'RH'
             FNAME(T,48) = 'Relative humidity'
             UNIT (T,48) = '%'
          CASE( 4  )
             NAME (T,48) = 'CF'
             FNAME(T,48) = 'Cloud fraction'
             UNIT (T,48) = 'unitless'
          CASE( 5  )
             NAME (T,48) = 'ODCOL'
             FNAME(T,48) = 'Cloud optical depth'
             UNIT (T,48) = 'unitless'
          CASE( 6  )
             NAME (T,48) = 'CThgt'
             FNAME(T,48) = 'Cloud top height'
             UNIT (T,48) = 'hPa'
          CASE( 7  )
             NAME (T,48) = 'AIRDEN'
             FNAME(T,48) = 'Air density'
             UNIT (T,48) = 'molec/cm3'
          CASE( 8  )
             NAME (T,48) = 'TOTSALT'
             FNAME(T,48) = 'Total seasalt tracer'
             UNIT (T,48) = 'ppbv'
             MWT  (T,48) = 36e-3
             SCALE(T,48) = 1.0e+9
          CASE( 9  )
             NAME (T,48) = 'D_LAI'
             FNAME(T,48) = 'Daily LAI'
             UNIT (T,48) = 'm2/m2'
          CASE( 10 )
             NAME (T,48) = 'NO_SOIL'
             FNAME(T,48) = 'Soil NO'
             UNIT (T,48) = 'molec/cm2/s'
          CASE( 11 )
             NAME (T,48) = 'NO_FERT'
             FNAME(T,48) = 'Fertilizer NO'
             UNIT (T,48) = 'molec/cm2/s'
          CASE( 12 )
             NAME (T,48) = 'P_DRY'
             FNAME(T,48) = 'Dry period'
             UNIT (T,48) = 'hours'
          CASE( 13 )
             NAME (T,48) = 'PFACTOR'
             FNAME(T,48) = 'Pulsing factor'
             UNIT (T,48) = 'unitless'
          CASE( 14 )
             NAME (T,48) = 'GWET_PRV'
             FNAME(T,48) = 'Soil moisture'
             UNIT (T,48) = 'unitless'
          CASE( 15 )
             NAME (T,48) = 'SZA'
             FNAME(T,48) = 'Cosine Solar Zenith Angle'
             UNIT (T,48) = 'unitless'
          CASE( 16 )
             NAME (T,48) = 'GLYX'
             FNAME(T,48) = 'Glyoxal time series'
             UNIT (T,48) = 'molec/cm3'
          CASE( 17 )
             NAME (T,48) = 'AERMASS'
             FNAME(T,48) = 'Aerosol mass time series'
             UNIT (T,48) = 'ug/m3'
          CASE DEFAULT
             ! Nothing
          END SELECT

       ENDDO
    ENDIF

#ifdef TOMAS
    !-------------------------------------
    ! TOMAS 3-D rate (ND60)
    !-------------------------------------
    IF ( ND60 > 0 ) THEN

       ! # of tracers
       NTRAC(60) = TMAX(60)

       ! Loop over tracers
       DO T = 1, NTRAC(60)

          ! Set pointer to this species information
          ! NOTE: requires 1:1 tracer index to species index mapping
          SpcInfo => State_Chm%SpcData(T)%Info

          IF ( ( T .ge. id_NK1       ) .and. &
               ( T .lt. id_NK1+IBINS ) ) THEN
             UNIT(T,60) = 'no.'
          ELSE
             UNIT(T,60) = 'kg'
          ENDIF

          ! Define the rest of the quantities
          NAME (T,60) = SpcInfo%Name
          FNAME(T,60) = 'Rates from ND60 TOMAS diagnostic'
          INDEX(T,60) = T + ( SPACING * 6 )
          MWT  (T,60) = 1e0
          SCALE(T,60) = 1e0
       ENDDO

    ENDIF

    !-------------------------------------
    ! TOMAS 3-D rate (ND61)                !(win, 7/9/09)
    !-------------------------------------
    IF ( ND61 > 0 ) THEN
       NTRAC(61) = 2

       ! Loop over tracers
       DO T = 1, NTRAC(61)

          ! Get name and long name for each tracer
          SELECT CASE( T )
          CASE( 1 )
             NAME (T,61) = 'NUC10NM'
             FNAME(T,61) = 'Particle Dp>10nm formation rate'
          CASE( 2 )
             NAME (T,61) = 'NUC_1NM'
             FNAME(T,61) = 'Nucleation rate (cluster size)'
          CASE DEFAULT
             !Nothing
          END SELECT

          ! Define the rest of the quantities
          UNIT (T,61) = 'cm-3s-1'
          INDEX(T,61) = T + ( SPACING * 7 )
          MWT  (T,61) = 1e0
          SCALE(T,61) = 1e0
       ENDDO
    ENDIF

    !-------------------------------------
    ! Family prod & loss (ND65)
    !-------------------------------------
    IF ( Input_Opt%DO_SAVE_PL ) THEN

       ! Number of P/L families
       NTRAC(65) = Input_Opt%NFAM

       ! Loop over each P/L family
       DO T = 1, NTRAC(65)
          NAME (T,65) = Input_Opt%FAM_NAME( T )
          FNAME(T,65) = TRIM( NAME(T,65) ) // ' P/L family'
          INDEX(T,65) = T + ( SPACING * 17 )
          MWT  (T,65) = 0e0
          SCALE(T,65) = 1e0

          ! Unit for Tag O3 is kg/s, otherwise molec/cm3/s
          IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
             UNIT(T,65) = 'kg/s'
          ELSE
             UNIT(T,65) = 'molec/cm3/s'
          ENDIF

          ! For sulfur tracking species, use 'kg S' unit and
          ! give a descriptive FNAME. (win, 8/9/09)
          IF ( NAME(T,65) == 'PSO21' ) THEN
             FNAME(T,65) = 'P(SO2) from DMS+OH (Online)'
             UNIT(T,65)  = 'kg S'
             MWT(T,65)   = 64e-3
          ELSE IF ( NAME(T,65) == 'PSO22' ) THEN
             FNAME(T,65) = 'P(SO2) from DMS+NO3 (Online)'
             UNIT(T,65)  = 'kg S'
             MWT(T,65)   = 64e-3
          ELSE IF ( NAME(T,65) == 'PSO4' ) THEN
             FNAME(T,65) = 'P(H2SO4) from SO2+OH (Online)'
             UNIT(T,65)  = 'kg S'
          ELSE IF ( NAME(T,65) == 'PMSA' ) THEN
             FNAME(T,65) = 'P(MSA) from DMS (Online)'
             UNIT(T,65) = 'kg S'
          ENDIF

       ENDDO

    ENDIF
#endif

    IF ( DO_TIMESERIES ) THEN

       !-------------------------------------      
       ! 3-D GMAO met fields (ND66)
       ! also for timeseries diagnostics
       !-------------------------------------      

       ! Number of tracers
       NTRAC(66) = 7

       ! Loop over tracers
       DO T = 1, NTRAC(66)

          ! Get name and unit for each met field
          SELECT CASE( T )
          CASE( 1  )
             NAME(T,66) = 'UWND'
             UNIT(T,66) = 'm/s'
          CASE( 2  )
             NAME(T,66) = 'VWND'
             UNIT(T,66) = 'm/s'
          CASE( 3  )
             NAME(T,66) = 'TMPU'
             UNIT(T,66) = 'K'
          CASE( 4 )
             NAME(T,66) = 'SPHU'
             UNIT(T,66) = 'g/kg'
          CASE( 5 )
             NAME(T,66) = 'CMFMC'
             UNIT(T,66) = 'kg/m2/s'
          CASE( 6 )
             NAME(T,66) = 'DTRAIN'
             UNIT(T,66) = 'kg/m2/s'
          CASE( 7 )
             NAME(T,66) = 'OMEGA'
             UNIT(T,66) = 'Pa/s'
          CASE DEFAULT
             ! Nothing
          END SELECT

          ! Define the rest of the quantities
          FNAME(T,66) = 'GMAO ' // TRIM( NAME(T,66) ) // ' field'
          INDEX(T,66) = T + ( SPACING * 12 )
          MWT  (T,66) = 0e0
          SCALE(T,66) = 1e0
       ENDDO

       !-------------------------------------      
       ! 2-D GMAO met fields (ND67)
       !-------------------------------------      

       ! Number of tracers
       NTRAC(67) = 23 ! (Lin, 03/31/09)

       ! Loop over tracers
       DO T = 1, NTRAC(67)

          ! Get name and unit for each met field
          SELECT CASE( T )
          CASE( 1  )
             NAME(T,67) = 'HFLUX'
             UNIT(T,67) = 'W/m2'
          CASE( 2  )
             NAME(T,67) = 'RADSWG'
             UNIT(T,67) = 'W/m2'
          CASE( 3  )
             NAME(T,67) = 'PREACC'
             UNIT(T,67) = 'mm/day'
          CASE( 4  )
             NAME(T,67) = 'PRECON'
             UNIT(T,67) = 'mm/day'
          CASE( 5  )
             NAME(T,67) = 'TS'
             UNIT(T,67) = 'K'
          CASE( 6  )
             NAME(T,67) = 'RADSWT'
             UNIT(T,67) = 'W/m2'
          CASE( 7  )
             NAME(T,67) = 'USTAR'
             UNIT(T,67) = 'm/s'
          CASE( 8  )
             NAME(T,67) = 'Z0'
             UNIT(T,67) = 'm'
          CASE( 9  )
             NAME(T,67) = 'PBL'
             UNIT(T,67) = 'm'
          CASE( 10 )
             NAME(T,67) = 'CLDFRC'
             UNIT(T,67) = 'unitless'
          CASE( 11 )
             NAME(T,67) = 'U10M'
             UNIT(T,67) = 'm/s'
          CASE( 12 )
             NAME(T,67) = 'V10M'
             UNIT(T,67) = 'm/s'
          CASE( 13 )
             NAME(T,67) = 'PS_PBL'
             UNIT(T,67) = 'hPa'
          CASE( 14 )
             NAME(T,67) = 'ALBD'
             UNIT(T,67) = 'unitless'
          CASE( 15 )
             NAME(T,67) = 'PHIS'
             UNIT(T,67) = 'm'
          CASE( 16 )
             NAME(T,67) = 'CLDTOP'
             UNIT(T,67) = 'level'
          CASE( 17 )
             NAME(T,67) = 'TROPPRAW'
             UNIT(T,67) = 'hPa'
          CASE( 18 )
             NAME(T,67) = 'SLP'
             UNIT(T,67) = 'hPa'
          CASE( 19 )
             NAME(T,67) = 'TSKIN'
             UNIT(T,67) = 'K'
          CASE( 20 )
             NAME(T,67) = 'PARDF'
             UNIT(T,67) = 'W/m2'
          CASE( 21 )
             NAME(T,67) = 'PARDR'
             UNIT(T,67) = 'W/m2'
          CASE( 22 )
             NAME(T,67) = 'GWET'
             UNIT(T,67) = 'unitless'
          CASE( 23 )
             ! Add EFLUX (Lin, 05/16/08)
             NAME(T,67) = 'EFLUX'
             UNIT(T,67) = 'W/m2'
          CASE DEFAULT
             ! Nothing
          END SELECT

          ! Define the rest of the quantities
          FNAME(T,67) = 'GMAO ' // TRIM( NAME(T,67) ) // ' field'
          INDEX(T,67) = T + ( SPACING * 11 )
          MWT  (T,67) = 0e0
          SCALE(T,67) = 1e0
       ENDDO

       !-------------------------------------      
       ! Grid box heights and related 
       ! quantities (ND68) + timeseries
       !-------------------------------------      

       ! Number of tracers
       NTRAC(68) = 8

       ! Loop over tracers
       DO T = 1, NTRAC(68)

          ! Get name and unit for each met field
          SELECT CASE( T )
          CASE( 1 )
             NAME (T,68) = 'BXHEIGHT'
             FNAME(T,68) = 'Grid box height'
             UNIT (T,68) = 'm'
          CASE( 2 )
             NAME (T,68) = 'AD'
             FNAME(T,68) = 'Air mass in grid box'
             UNIT (T,68) = 'kg'
          CASE( 3 )
             NAME (T,68) = 'AVGW'
             FNAME(T,68) = 'Mixing ratio of H2O vapor'
             UNIT (T,68) = 'v/v'
          CASE( 4 )
             NAME (T,68) = 'AIRNUMDEN'
             FNAME(T,68) = 'Dry air number density'
             UNIT (T,68) = 'molec air/cm3'
          CASE( 5 )
             NAME (T,68) = 'T'
             FNAME(T,68) = 'Temperature'
             UNIT (T,68) = 'K'
          CASE( 6 )
             NAME (T,68) = 'PMID'
             FNAME(T,68) = 'Pressure at average pressure level'
             UNIT (T,68) = 'hPa'
          CASE( 7 )
             NAME (T,68) = 'PEDGE'
             FNAME(T,68) = 'Pressure at grid box lower edge'
             UNIT (T,68) = 'hPa'
          CASE( 8 )
             NAME (T,68) = 'RH'
             FNAME(T,68) = 'Relative humidity'
             UNIT (T,68) = '%'
          CASE DEFAULT
             ! Nothing
          END SELECT

          ! Define the rest of the quantities
          INDEX(T,68) = T + ( SPACING * 24 )
          MWT  (T,68) = 0e0
          SCALE(T,68) = 1e0
       ENDDO

    ENDIF ! DO_TIMERSERIES

    ! Nullify pointer
    SpcInfo  => NULL()

  END SUBROUTINE INIT_TRACERINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gamap
!
! !DESCRIPTION: Subroutine INIT\_GAMAP allocates and initializes all
!  module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_GAMAP( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,     ONLY : ALLOC_ERR
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Chm_Mod, ONLY : Ind_
    USE TIME_MOD,      ONLY : EXPAND_DATE, GET_NHMSb, GET_NYMDb
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  22 Apr 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: AS, NYMDb, NHMSb
    INTEGER                :: MAXTRACER_HG, N

    !=================================================================
    ! INIT_GAMAP begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Save from arguments to module variables
    DFILE = Input_Opt%GAMAP_DIAGINFO
    TFILE = Input_Opt%GAMAP_TRACERINFO

    ! Get starting date & time
    NYMDb = GET_NYMDb()
    NHMSb = GET_NHMSb()

    ! Replace any date/time tokens in the file names
    CALL EXPAND_DATE( DFILE, NYMDb, NHMSb )
    CALL EXPAND_DATE( TFILE, NYMDb, NHMSb )

    !=================================================================
    ! Define species ID flags
    !=================================================================
    id_NK1   = Ind_('NK1'  )

    !=================================================================
    ! Allocate module arrays
    !=================================================================

    ALLOCATE( OFFSET( MAXCAT ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'OFFSET' )
    OFFSET = 0

    ALLOCATE( CATEGORY( MAXCAT ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CATEGORY' )
    CATEGORY = ''

    ALLOCATE( DESCRIPT( MAXCAT ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'DESCRIPT' )
    DESCRIPT = ''

    ALLOCATE( NTRAC( MAXDIAG ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'NTRAC' )
    NTRAC = 0

    IF ( Input_Opt%ITS_A_MERCURY_SIM .and. Input_Opt%LSPLIT) THEN

       MAXTRACER_HG = 480

       ALLOCATE( INDEX( MAXTRACER_HG, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'INDEX' )
       INDEX = 0

       ALLOCATE( MWT( MAXTRACER_HG, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'MWT' )
       MWT = 0.0

       ALLOCATE( SCALE( MAXTRACER_HG, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'SCALE' )
       SCALE = 0.0

       ALLOCATE( NAME( MAXTRACER_HG, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'NAME' )
       NAME = ''

       ALLOCATE( FNAME( MAXTRACER_HG, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'FNAME' )
       FNAME = ''

       ALLOCATE( UNIT( MAXTRACER_HG, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'UNIT' )
       UNIT = ''

    ELSE

       ALLOCATE( INDEX( MAXTRACER, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'INDEX' )
       INDEX = 0

       ALLOCATE( MWT( MAXTRACER, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'MWT' )
       MWT = 0.0

       ALLOCATE( SCALE( MAXTRACER, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'SCALE' )
       SCALE = 0.0

       ALLOCATE( NAME( MAXTRACER, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'NAME' )
       NAME = ''

       ALLOCATE( FNAME( MAXTRACER, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'FNAME' )
       FNAME = ''

       ALLOCATE( UNIT( MAXTRACER, MAXDIAG ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'UNIT' )
       UNIT = ''

    ENDIF

    !=================================================================
    ! Initialize arrays for "diaginfo.dat" & "tracerinfo.dat" files
    !=================================================================

    ! Initialize arrays for "diaginfo.dat"
    CALL INIT_DIAGINFO()

    ! Initialize arrays for "tracerinfo.dat"
    CALL INIT_TRACERINFO( Input_Opt, State_Chm, RC )

  END SUBROUTINE INIT_GAMAP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_gamap
!
! !DESCRIPTION: Subroutine CLEANUP\_GAMAP deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_GAMAP
!
! !REVISION HISTORY:
!  25 Apr 2005 - R. Yantosca - Initial version
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( CATEGORY ) ) DEALLOCATE( CATEGORY )
    IF ( ALLOCATED( DESCRIPT ) ) DEALLOCATE( DESCRIPT )
    IF ( ALLOCATED( FNAME    ) ) DEALLOCATE( FNAME    )
    IF ( ALLOCATED( INDEX    ) ) DEALLOCATE( INDEX    )
    IF ( ALLOCATED( MWT      ) ) DEALLOCATE( MWT      )
    IF ( ALLOCATED( NAME     ) ) DEALLOCATE( NAME     )
    IF ( ALLOCATED( NTRAC    ) ) DEALLOCATE( NTRAC    )
    IF ( ALLOCATED( OFFSET   ) ) DEALLOCATE( OFFSET   )
    IF ( ALLOCATED( SCALE    ) ) DEALLOCATE( SCALE    )
    IF ( ALLOCATED( UNIT     ) ) DEALLOCATE( UNIT     )

  END SUBROUTINE CLEANUP_GAMAP
!EOC
END MODULE GAMAP_MOD
#endif
