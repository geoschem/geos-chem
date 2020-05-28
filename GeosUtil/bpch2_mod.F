!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: bpch2_mod.F
!
! !DESCRIPTION: Module BPCH2\_MOD contains the routines used to read data
!  from and write data to binary punch (BPCH) file format (v. 2.0).
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
      PUBLIC  :: GET_TAU0

      INTERFACE GET_TAU0
         MODULE PROCEDURE GET_TAU0_6A
      END INTERFACE
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GET_TAU0_6A
!
! !REMARKS:
!  ###########################################################################
!  ##### BINARY PUNCH INPUT IS BEING PHASED OUT.  MOST INPUT IS NOW READ #####
!  ##### FROM COARDS-COMPLIANT netCDF FILES VIA HEMCO (bmy, 4/1/15)      #####
!  ###########################################################################
!
! !REVISION HISTORY:
!  (1 ) Added routine GET_TAU0 (bmy, 7/20/00)
!  (2 ) Added years 1985-2001 for routine GET_TAU0 (bmy, 8/1/00)
!  (3 ) Use IOS /= 0 criterion to also check for EOF (bmy, 9/12/00)
!  (4 ) Removed obsolete code in "read_bpch2.f" (bmy, 12/18/00)
!  (5 ) Correct error for 1991 TAU values in GET_TAU0 (bnd, bmy, 1/4/01)
!  (6 ) BPCH2_MOD is now independent of any GEOS-CHEM size parameters.
!        (bmy, 4/18/01)
!  (7 ) Now have 2 versions of "GET_TAU0" overloaded by an interface.  The
!        original version takes 2 arguments (MONTH, YEAR).  The new version
!        takes 3 arguments (MONTH, DAY, YEAR). (bmy, 8/22/01)
!  (8 ) Updated comments (bmy, 9/4/01)
!  (9 ) Renamed GET_TAU0_3A to GET_TAU0_6A, and updated the GET_TAU0
!        interface.  Also updated comments (bmy, 9/26/01)
!  (10) Now use special model name for GEOS-3 w/ 30 layers (bmy, 10/9/01)
!  (11) Minor bug fix in GET_TAU0_2A.  Also deleted obsolete code from 9/01.
!        (bmy, 11/15/01)
!  (12) Moved routines JULDAY, MINT, CALDATE to "julian_mod.f".  Now
!        references routine JULDAY from "julday_mod.f".  Also added code
!        for GEOS-4/fvDAS model type. (bmy, 11/20/01)
!  (23) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Also add MODULE INTERFACES section,
!        since we have an interface here. (bmy, 5/28/02)
!  (24) Added OPEN_BPCH2_FOR_READ and OPEN_BPCH2_FOR_WRITE.  Also now
!        reference IU_FILE and IOERROR from "file_mod.f". (bmy, 7/30/02)
!  (25) Now references "error_mod.f".  Also obsoleted routine GET_TAU0_2A.
!        (bmy, 10/15/02)
!  (26) Made modification in READ_BPCH2 for 1x1 nested grids (bmy, 3/11/03)
!  (27) Modifications for GEOS-4, 30-layer grid (bmy, 11/3/03)
!  (28) Added cpp switches for GEOS-4 1x125 grid (bmy, 12/1/04)
!  (29) Modified for GCAP and GEOS-5 met fields.  Added function
!        GET_HALFPOLAR. (bmy, 6/28/05)
!  (30) Added GET_NAME_EXT_2D to get filename extension for files which do
!        not contain any vertical information (bmy, 8/16/05)
!  (31) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (32) Renamed GRID30LEV to GRIDREDUCED.  Also increase TEMPARRAY in
!        READ_BPCH2 for GEOS-5 vertical levels. (bmy, 2/16/07)
!  (33) Modifications for GEOS-5 nested grids (bmy, 11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Added modifications for MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  19 Jul 2012 - R. Yantosca - Bug fix in GET_NAME_EXT_2D
!  03 Aug 2012 - R. Yantosca - Reference file LUN routines from inquireMod.F90
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  11 Aug 2015 - R. Yantosca - Add modifications for MERRA2 data
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
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
!  (1 ) Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  06 Aug 2012 - R. Yantosca - Do not call findFreeLun() in this subroutine
!                              but instead in the calling routine and pass
!                              the file unit as an argument.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TPBC ) || defined( BPCH_TIMESER )
      INTEGER                                  :: IOS
      CHARACTER(LEN=40)                        :: FTI
      CHARACTER(LEN=80)                        :: TMP_TITLE

      !=================================================================
      ! OPEN_BPCH2_FOR_READ begins here!
      !=================================================================

      ! Open file for input -- readonly
      OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='OLD',
     &      IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

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
         CALL ERROR_STOP( 'Invalid file format!',
     &                    'OPEN_BPCH2_FOR_READ (bpch2_mod.f)')
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

#endif
      ! Return to calling program
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
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  06 Aug 2012 - R. Yantosca - Do not call findFreeLun() in this subroutine
!                              but instead in the calling routine and pass
!                              the file unit as an argument.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TPBC ) || defined( BPCH_TIMESER )

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
      OPEN( IUNIT,      FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &      IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

      ! Error check
      IF ( IOS /= 0 ) THEN
         WRITE(6,*) ' '
         WRITE(6,*) "CANNOT WRITE : " // FILENAME
         CALL IOERROR( IOS, IUNIT,'open_bpch2_for_write:1')
      ENDIF


      ! Write the top-of-file title to disk
      CALL BPCH2_HDR( IUNIT, TMP_TITLE )

#endif
      ! Return to calling program
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
!  (1 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (2 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!  (3 ) Now reference IOERROR from "file_mod.f". (bmy, 6/26/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TPBC ) || defined( BPCH_TIMESER )

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

#endif

      ! Return to calling program
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
      SUBROUTINE BPCH2( IUNIT,     MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NTRACER,
     &                  UNIT,      TAU0,      TAU1,     RESERVED,
     &                  NI,        NJ,        NL,       IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
      ! Arguments
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
!  (1 ) Added indices to IOERROR calls (e.g. "bpch2:1", "bpch2:2", etc.)
!        (bmy, 10/4/99)
!  (2 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (3 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!  (4 ) Now reference IOERROR from "file_mod.f". (bmy, 6/26/02)
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!EOP
!------------------------------------------------------------------------------
!BOC
#if defined( BPCH_DIAG ) || defined( BPCH_TPBC ) || defined( BPCH_TIMESER )
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
      WRITE( IUNIT, IOSTAT=IOS )
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:1' )

      WRITE( IUNIT, IOSTAT = IOS )
     &     CATEGORY, NTRACER,  UNIT, TAU0,   TAU1,   RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:2' )

      WRITE( IUNIT, IOSTAT=IOS )
     &     ( ( ( ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:3' )

      !=================================================================
      ! Return to calling program
      !=================================================================
#endif
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
      FUNCTION GET_MODELNAME( Input_Opt, State_Grid )
     &    RESULT( MODELNAME )
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
!  (1 ) Now use special model name for GEOS-3 w/ 30 layers (bmy, 10/9/01)
!  (2 ) Added modelname for GEOS-4/fvDAS model type (bmy, 11/20/01)
!  (3 ) Added "GEOS4_30L" for reduced GEOS-4 grid.  Also now use C-preprocessor
!        switch "GRID30LEV" instead of IF statements. (bmy, 11/3/03)
!  (4 ) Updated for GCAP and GEOS-5 met fields.  Rearranged coding for
!        simplicity. (swu, bmy, 5/24/05)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Rename GRID30LEV to GRIDREDUCED (bmy, 2/7/07)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  13 Aug 2010 - R. Yantosca - Added MERRA model names
!  01 Feb 2012 - R. Yantosca - Added GEOS-5.7.x model names
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  11 Dec 2012 - R. Yantosca - Bug fix: Need to specify both EXTERNAL_GRID and
!                              EXTERNAL_FORCING Cpp switches
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  11 Aug 2015 - R. Yantosca - Add support for MERRA2 data
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

      ! Return to calling program
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
!  20 Nov 2009 - R. Yantosca         - Added ProTeX header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! All GEOS grids have half-sized polar boxes
      HALFPOLAR = 1

      ! Return to calling program
      END FUNCTION GET_HALFPOLAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Tau0_6a
!
! !DESCRIPTION: Function GET\_TAU0\_6A returns the corresponding TAU0 value
!  for the first day of a given MONTH of a given YEAR.  This is necessary to
!  index monthly mean binary punch files, which are used as input to GEOS-Chem.
!\\
!\\
!  This function takes 3 mandatory arguments (MONTH, DAY, YEAR) and 3
!  optional arguments (HOUR, MIN, SEC).  It is intended to replace the current
!  2-argument version of GET\_TAU0.  The advantage being that GET\_TAU0\_6A
!  can compute a TAU0 for any date and time in the GEOS-Chem epoch, rather
!  than just the first day of each month.  Overload this w/ an interface so
!  that the user can also choose the version of GET\_TAU0 w/ 2 arguments
!  (MONTH, YEAR), which is the prior version.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TAU0_6A( MONTH, DAY, YEAR,
     &                      HOUR,  MIN, SEC  ) RESULT( THIS_TAU0 )
!
! !USES:
!
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE JULDAY_MOD, ONLY : JULDAY
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)           :: MONTH
      INTEGER, INTENT(IN)           :: DAY
      INTEGER, INTENT(IN)           :: YEAR
      INTEGER, INTENT(IN), OPTIONAL :: HOUR
      INTEGER, INTENT(IN), OPTIONAL :: MIN
      INTEGER, INTENT(IN), OPTIONAL :: SEC
!
! !RETURN VALUE:
!
      REAL(f8)                        :: THIS_TAU0   ! TAU0 timestamp
!
! !REMARKS:
!  TAU0 is hours elapsed since 00:00 GMT on 01 Jan 1985.
!
! !REVISION HISTORY:
!  (1 ) 1985 is the first year of the GEOS epoch.
!  (2 ) Add TAU0 values for years 1985-2001 (bmy, 8/1/00)
!  (3 ) Correct error for 1991 TAU values.  Also added 2002 and 2003.
!        (bnd, bmy, 1/4/01)
!  (4 ) Updated comments  (bmy, 9/26/01)
!  (5 ) Now references JULDAY from "julday_mod.f" (bmy, 11/20/01)
!  (6 ) Now references ERROR_STOP from "error_mod.f"  (bmy, 10/15/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: TMP_HOUR, TMP_MIN, TMP_SEC
      REAL(f8) :: DAYS

      !=================================================================
      ! GET_TAU0_6A begins here!
      !=================================================================

      ! Error checking
      IF ( MONTH < 1 .or. MONTH > 12 ) THEN
         CALL ERROR_STOP ( 'Invalid MONTH selection!', 'GET_TAU0' )
      ENDIF

      ! Error checking
      IF ( DAY < 1 .or. DAY > 31 ) THEN
         CALL ERROR_STOP ( 'Invalid DAY selection!', 'GET_TAU0' )
      ENDIF

      ! If HOUR isn't passed, default to 0
      IF ( PRESENT( HOUR ) ) THEN
         TMP_HOUR = HOUR
      ELSE
         TMP_HOUR = 0
      ENDIF

      ! If MIN isn't passed, default to 0
      IF ( PRESENT( MIN ) ) THEN
         TMP_MIN = MIN
      ELSE
         TMP_MIN = 0
      ENDIF

      ! If SEC isn't passed, default to 0
      IF ( PRESENT( SEC ) ) THEN
         TMP_SEC = SEC
      ELSE
         TMP_SEC = 0
      ENDIF

      ! Number of days since midnight on 1/1/1985
      THIS_TAU0 = JULDAY( YEAR, MONTH, DBLE( DAY ) ) - 2446066.5e+0_f8

      ! Multiply by 24 to get hours since 1/1/1985
      ! Also add in the hours elapsed since midnight on this date
      THIS_TAU0 = ( THIS_TAU0 * 24e+0_f8 ) + ( TMP_HOUR             ) +
     &            ( TMP_MIN   / 60e+0_f8 ) + ( TMP_SEC / 3600e+0_f8 )

      ! Return to calling program
      END FUNCTION GET_TAU0_6A
!EOC

      ! End of module
      END MODULE BPCH2_MOD
