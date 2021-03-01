#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: bpch2_mod.F90
!
! !DESCRIPTION: Module BPCH2\_MOD contains the routines used to read data
!  from and write data to binary punch (BPCH) file format (v2.0).
!\\
!\\
! !INTERFACE:
!
MODULE BPCH2_MOD
!
! !USES:
!
  USE inquireMod, ONLY : findFreeLUN
  USE inquireMod, ONLY : I_Am_UnOPENed

  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: OPEN_BPCH2_FOR_READ
  PUBLIC  :: OPEN_BPCH2_FOR_WRITE
  PUBLIC  :: BPCH2
  PUBLIC  :: BPCH2_HDR
  PUBLIC  :: GET_MODELNAME
  PUBLIC  :: GET_HALFPOLAR
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE WILL BE A STUB UNLESS GEOS-Chem IS COMPILED    %%%
!  %%%  WITH THE BPCH_DIAG=y OPTION. (bmy, 10/4/19)                      %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Open_Bpch2_For_Read
!
! !DESCRIPTION: Subroutine OPEN\_BPCH2\_FOR\_READ opens a binary punch file
!  (version 2.0 format) for reading only.  Also reads FTI and TITLE strings.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OPEN_BPCH2_FOR_READ( IUNIT, FILENAME, TITLE )
!
! !USES:
!
    USE ERROR_MOD, ONLY : ERROR_STOP
    USE FILE_MOD,  ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN)            :: IUNIT     ! LUN for file I/O
    CHARACTER(LEN=*),  INTENT(IN)            :: FILENAME  ! File to open
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=80), INTENT(OUT), OPTIONAL :: TITLE     ! File title string
!
! !REMARKS:
!  ###########################################################################
!  ##### BINARY PUNCH INPUT IS BEING PHASED OUT.  MOST INPUT IS NOW READ #####
!  ##### FROM COARDS-COMPLIANT netCDF FILES VIA HEMCO (bmy, 4/1/15)      #####
!  ###########################################################################
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                                  :: IOS
    CHARACTER(LEN=40)                        :: FTI
    CHARACTER(LEN=80)                        :: TMP_TITLE

    !=================================================================
    ! OPEN_BPCH2_FOR_READ begins here!
    !=================================================================

    ! Open file for input -- readonly
    OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='OLD', &
          IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

    ! Error check
    IF ( IOS /= 0 ) THEN
       WRITE(6,*)'Error opening filename=',trim(filename)
       CALL FLUSH(6)
       CALL IOERROR( IOS, IUNIT, 'open_bpch2_for_read:1')
    ENDIF

    ! Read file type identifier
    READ( IUNIT, IOSTAT=IOS ) FTI

    ! Error check
    IF ( IOS /= 0 ) THEN
       WRITE(6,*)'Error reading FTI for filename=',trim(filename)
       CALL FLUSH(6)
       CALL IOERROR( IOS, IUNIT, 'open_bpch2_for_read:2' )
    ENDIF

    ! Stop if this is not a binary punch file
    IF ( TRIM( FTI ) /= 'CTM bin 02' ) THEN
       WRITE(6,*)'Error filename=',trim(filename)
       CALL FLUSH(6)
       CALL ERROR_STOP( 'Invalid file format!', &
                         'OPEN_BPCH2_FOR_READ (bpch2_mod.F90)')
    ENDIF

    ! Read top title
    READ( IUNIT, IOSTAT=IOS ) TMP_TITLE

    ! Error check
    IF ( IOS /= 0 ) THEN
       WRITE(6,*)'Error reading filename=',trim(filename)
       CALL FLUSH(6)
       CALL IOERROR( IOS, IUNIT, 'open_bpch2_for_read:3' )
    ENDIF

    ! Copy value of TMP_TITLE to TITLE for return
    IF ( PRESENT( TITLE ) ) TITLE = TMP_TITLE

  END SUBROUTINE OPEN_BPCH2_FOR_READ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Open_Bpch2_For_Write
!
! !DESCRIPTION: Subroutine OPEN\_BPCH2\_FOR\_WRITE opens a binary punch file
!  (version 2.0) for writing.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OPEN_BPCH2_FOR_WRITE( IUNIT, FILENAME, TITLE )
!
! !USES:
!
    USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN)            :: IUNIT     ! LUN for file I/O
    CHARACTER(LEN=*),  INTENT(IN)            :: FILENAME  ! Name of file
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=80), INTENT(OUT), OPTIONAL :: TITLE     ! File title string
!
! !REMARKS:
!  ###########################################################################
!  ##### BINARY PUNCH INPUT IS BEING PHASED OUT.  MOST INPUT IS NOW READ #####
!  ##### FROM COARDS-COMPLIANT netCDF FILES VIA HEMCO (bmy, 4/1/15)      #####
!  ###########################################################################
!
! !REVISION HISTORY:
!  30 Jul 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: IOS
    CHARACTER(LEN=80) :: TMP_TITLE

    !=================================================================
    ! OPEN_BPCH2_FOR_WRITE begins here!
    !=================================================================

    ! If TITLE is not passed, create a default title string
    IF ( PRESENT( TITLE ) ) THEN
       TMP_TITLE = TITLE
    ELSE
       TMP_TITLE = 'GEOS-CHEM binary punch file v. 2.0'
    ENDIF

    ! Open file for output
    OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN', &
          IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

    ! Error check
    IF ( IOS /= 0 ) THEN
       WRITE(6,*) ' '
       WRITE(6,*) "CANNOT WRITE : " // FILENAME
       CALL IOERROR( IOS, IUNIT,'open_bpch2_for_write:1')
    ENDIF

    ! Write the top-of-file title to disk
    CALL BPCH2_HDR( IUNIT, TMP_TITLE )

  END SUBROUTINE OPEN_BPCH2_FOR_WRITE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Bpch2_hdr
!
! !DESCRIPTION: Subroutine BPCH2\_HDR writes a header at the top of the binary
!  punch file, version 2.0.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE BPCH2_HDR ( IUNIT, TITLE )
!
! !USES:
!
    USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN) :: IUNIT   ! LUN for file I/O
    CHARACTER(LEN=80), INTENT(IN) :: TITLE   ! Top-of-file title string
!
! !REMARKS:
!  ###########################################################################
!  ##### BINARY PUNCH INPUT IS BEING PHASED OUT.  MOST INPUT IS NOW READ #####
!  ##### FROM COARDS-COMPLIANT netCDF FILES VIA HEMCO (bmy, 4/1/15)      #####
!  ###########################################################################
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                       :: IOS
    CHARACTER(LEN=40)             :: FTI = 'CTM bin 02'

    !=================================================================
    ! BPCH2_HDR begins here!
    !
    ! Write header information to binary punch file
    ! Also be sure to trap I/O Error conditions
    !=================================================================

    WRITE ( IUNIT, IOSTAT=IOS ) FTI
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2_hdr:1' )

    WRITE ( IUNIT, IOSTAT=IOS ) TITLE
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2_hdr:2' )

  END SUBROUTINE BPCH2_HDR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Bpch2
!
! !DESCRIPTION: Subroutine BPCH2 writes binary punch file (version 2.0) to
!  disk.  Information about the model grid is also stored with each data block.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE BPCH2( IUNIT,     MODELNAME, LONRES,   LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY, NTRACER,  &
                    UNIT,      TAU0,      TAU1,     RESERVED, &
                    NI,        NJ,        NL,       IFIRST,   &
                    JFIRST,    LFIRST,    ARRAY )
!
! !USES:
!
    USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN) :: IUNIT            ! LUN for file I/O
    CHARACTER(LEN=20), INTENT(IN) :: MODELNAME        ! Met field type
    REAL*4,            INTENT(IN) :: LONRES           ! Lon resolution [deg]
    REAL*4,            INTENT(IN) :: LATRES           ! Lat resolution [deg]
    INTEGER,           INTENT(IN) :: HALFPOLAR        ! 1/2-size polar boxes?
    INTEGER,           INTENT(IN) :: CENTER180        ! 1st box center -180?
    CHARACTER(LEN=40), INTENT(IN) :: CATEGORY         ! Diag. category name
    INTEGER,           INTENT(IN) :: NTRACER          ! Tracer index #
    CHARACTER(LEN=40), INTENT(IN) :: UNIT             ! Unit string
    REAL(f8),          INTENT(IN) :: TAU0             ! TAU values @ start &
    REAL(f8),          INTENT(IN) :: TAU1             !  end of diag interval
    CHARACTER(LEN=40), INTENT(IN) :: RESERVED         ! Extra string
    INTEGER,           INTENT(IN) :: NI, NJ, NL       ! Dimensions of ARRAY
    INTEGER,           INTENT(IN) :: IFIRST           ! (I,J,L) indices of
    INTEGER,           INTENT(IN) :: JFIRST           !  the first grid box
    INTEGER,           INTENT(IN) :: LFIRST           !  in Fortran notation
    REAL*4,            INTENT(IN) :: ARRAY(NI,NJ,NL)  ! Data array
!
! !REMARKS:
!  ############################################################################
!  ##### BINARY PUNCH OUTPUT IS BEING PHASED OUT. (bmy, 4/1/15)           #####
!  ############################################################################
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                       :: I, J, L, NSKIP, IOS
!
! !DEFINED PARAMETERS:
!
    INTEGER, PARAMETER            :: BYTES_PER_NUMBER = 4
    INTEGER, PARAMETER            :: END_OF_RECORD    = 8

    !=================================================================
    ! BPCH2 begins here!!
    !
    ! Compute the number of bytes to skip between the end of one
    ! data block and the beginning of the next data header line
    !=================================================================
    NSKIP = ( BYTES_PER_NUMBER * ( NI * NJ * NL ) ) + END_OF_RECORD

    !=================================================================
    ! Write data block to binary punch file
    ! Check for I/O errors
    !=================================================================
    WRITE( IUNIT, IOSTAT=IOS ) &
         MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

    IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:1' )

    WRITE( IUNIT, IOSTAT = IOS ) &
         CATEGORY, NTRACER,  UNIT, TAU0,   TAU1,   RESERVED, &
         NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST, NSKIP

    IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:2' )

    WRITE( IUNIT, IOSTAT=IOS ) &
         ( ( ( ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

    IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:3' )

  END SUBROUTINE BPCH2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Modelname
!
! !DESCRIPTION: Function GET\_MODELNAME returns the proper value of MODELNAME
!  for current met field type.  MODELNAME is written to the binary punch file
!  and is also used by the GAMAP package.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_MODELNAME( Input_Opt, State_Grid ) &
       RESULT( MODELNAME )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State
!
! !RETURN VALUE:
!
    CHARACTER(LEN=20) :: MODELNAME   ! Model name for the current met field
!
! !REMARKS:
!  We now read many data files via HEMCO, so we don't have much of a need
!  of constructing file names w/in the code.  This routine is now pretty
!  much obsolete and is slated for eventual removal.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( EXTERNAL_FORCING ) || defined ( EXTERNAL_GRID )
    MODELNAME = 'EXTERNAL'
#else
    IF ( State_Grid%NZ == 47 ) THEN
       MODELNAME = TRIM(Input_Opt%MetField) // '_47L'
    ELSE
       MODELNAME = TRIM(Input_Opt%MetField)
    ENDIF
#endif

  END FUNCTION GET_MODELNAME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Halfpolar
!
! !DESCRIPTION: Function GET\_HALFPOLAR returns 1 if the current grid has
!  half-sized polar boxes (e.g. GEOS) or zero otherwise (e.g. GCAP).
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_HALFPOLAR() RESULT( HALFPOLAR )
!
! !RETURN VALUE:
!
    INTEGER :: HALFPOLAR  ! =1 if we have half-sized polar boxes, =0 if not
!
! !REVISION HISTORY:
!  28 Jun 2005 - S. Wu & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! All GEOS grids have half-sized polar boxes
    HALFPOLAR = 1

  END FUNCTION GET_HALFPOLAR
!EOC
END MODULE BPCH2_MOD
#endif
