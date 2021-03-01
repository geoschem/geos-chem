!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: input_gtmm_mod
!
! !DESCRIPTION: Module INPUT\_GTMM\_MOD contains subroutines to read 
!  input file for GTMM.
!\\
!\\
! !INTERFACE:
!
MODULE INPUT_GTMM_MOD
!
! !USES:
!
  USE defineConstants
  USE defineArrays

  IMPLICIT NONE
  PRIVATE
!
! !PRIVATE DATA MEMBERS:
! 
  INTEGER, PARAMETER :: FIRSTCOL = 38
  INTEGER, PARAMETER :: HEADER   = 8
  CHARACTER(LEN=255) :: FILENAME = 'input.gtmm'
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: READ_GTMM_INPUT_FILE
!
! !REVISION HISTORY:
!  17 Dec 2009 - C. Carouge   - Initial version
!EOP
!-----------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gtmm_input_file 
!
! !DESCRIPTION: Subroutine READ\_GTMM\_INPUT\_FILE reads the input.gtmm file. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_GTMM_INPUT_FILE
!
! !REVISION HISTORY:
!  17 Dec 2009 - C. Carouge  -- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    LOGICAL            :: EOF
    INTEGER            :: IOS, I, IU_GTMM
    CHARACTER(LEN=255) :: LINE

    IU_GTMM=20
    ! Open file
    OPEN( IU_GTMM, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       print*, 'Problem opening the input file for GTMM'
       stop
    ENDIF

    !Read header
    DO i=1,HEADER
       READ(IU_GTMM, '(a)') LINE
    ENDDO

    !Read data
    !Read a line from the file, exit if EOF
    READ(IU_GTMM,'(a)') LINE
       
    !Keep only the data part of the line
    LINE=LINE(FIRSTCOL:)

    READ(LINE,*) NPPequilibriumYear
       
    !Read a line from the file, exit if EOF
    READ(IU_GTMM,'(a)') LINE
       
    !Keep only the data part of the line
    LINE=LINE(FIRSTCOL:)

    READ(LINE,*) HgPoolsequilibriumYear

    !Initialize indYear too.
    indYear = HgPoolsequilibriumYear
       
    !Read a line from the file, exit if EOF
    READ(IU_GTMM,'(a)') LINE
       
    !Keep only the data part of the line
    LINE=LINE(FIRSTCOL:)

    READ(LINE,*) preindYear
    
    !Read a line from the file, exit if EOF
    READ(IU_GTMM,'(a)') LINE
       
    !Keep only the data part of the line
    LINE=LINE(FIRSTCOL:)

    READ(LINE,*) LRESTART
    
    !Read a line from the file, exit if EOF
    READ(IU_GTMM,'(a)') LINE
       
    !Keep only the data part of the line
    LINE=LINE(FIRSTCOL:)

    READ(LINE,*) restartfile
    
    CLOSE(IU_GTMM)
  END SUBROUTINE READ_GTMM_INPUT_FILE
!EOC
END MODULE INPUT_GTMM_MOD
