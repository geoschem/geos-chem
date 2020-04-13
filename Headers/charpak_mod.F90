!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: charpak_mod.F90
!
! !DESCRIPTION: Module CHARPAK\_MOD contains routines from the CHARPAK
!  string and character manipulation package used by GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE Charpak_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CleanText
  PUBLIC  :: CntMat
  PUBLIC  :: CopyTxt
  PUBLIC  :: CStrip
  PUBLIC  :: IsDigit
  PUBLIC  :: ReadOneLine
  PUBLIC  :: StrRepl
  PUBLIC  :: StrSplit
  PUBLIC  :: StrSqueeze
  PUBLIC  :: To_UpperCase
  PUBLIC  :: TranLc
  PUBLIC  :: TranUc
  PUBLIC  :: Txtext
  PUBLIC  :: WordWrapPrint
!
! !PRIVATE MEMBER FUNCTIONS
!
!
! !REMARKS:
!  CHARPAK routines by Robert D. Stewart, 1992.  Subsequent modifications
!  made for GEOS-CHEM by Bob Yantosca (1998, 2002, 2004).
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS
!
  ! Maximum string length
  INTEGER, PARAMETER, PUBLIC :: MAXSTRLEN = 500

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CntMat
!
! !DESCRIPTION: Counts the number of characters in str1 that match
!  a character in str2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CntMat( Str1, Str2, Imat, Locations )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) ::  Str1             ! Text to scan
    CHARACTER(LEN=*), INTENT(IN) ::  Str2             ! Character to match
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: imat             ! Number of matches
    INTEGER,          OPTIONAL    :: Locations(255)   ! Positions of matches
!
! !REVISION HISTORY:
!     DATE:   JAN. 6, 1995
!     AUTHOR: R.D. STEWART
!     COMMENTS: Revised slightly (2-5-1996) so that trailing
!               blanks in str1 are ignored.  Revised again
!               on 3-6-1996.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: L1, L2, i, j
    LOGICAL :: again

    ! Arrays
    INTEGER :: TmpLocations(255)

    ! Initialize
    TmpLocations = 0
    L1           = MAX(1,LEN_TRIM(str1))
    L2           = LEN(str2)
    imat         = 0

    DO i=1,L1
       again = .true.
       j = 1
       DO WHILE (again)
          IF (str2(j:j).EQ.str1(i:i)) THEN
             imat               = imat+1
             TmpLocations(imat) = i
             again              = .false.
          ELSEIF (j.LT.L2) THEN
             j=j+1
          ELSE
             again = .false.
          ENDIF
       ENDDO
    ENDDO

    ! Return positions where matches occured (OPTIONAL)
    IF ( PRESENT( Locations ) ) Locations = TmpLocations

  END SUBROUTINE CntMat
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CopyTxt
!
! !DESCRIPTION: Write all of the characters in str1 into variable
!               str2 beginning at column, col.  If the length of str1
!               + col is longer than the number of characters str2
!               can store, some characters will not be transfered to
!               str2.  Any characters already existing in str2 will
!               will be overwritten.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CopyTxt( col, str1, str2 )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)    :: col
    CHARACTER(LEN=*), INTENT(IN)    :: str1
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(INOUT) :: str2
!
! !REVISION HISTORY:
!     DATE:   DEC. 24, 1993
!     AUTHOR: R.D. STEWART
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: ilt1,i1,i,j,ic

    i1 = LEN(str2)
    IF (i1.GT.0) THEN
       ilt1 = LEN(str1)
       IF (ilt1.GT.0) THEN
          ic = MAX0(col,1)
          i = 1
          j = ic
          DO WHILE ((i.LE.ilt1).and.(j.LE.i1))
             str2(j:j) = str1(i:i)
             i = i + 1
             j = ic + (i-1)
          ENDDO
       ENDIF
    ENDIF

  END SUBROUTINE CopyTxt
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cstrip
!
! !DESCRIPTION: Strip blanks and null characters for the variable TEXT.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CStrip( text, KeepSpaces )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          OPTIONAL      :: KeepSpaces ! If =T, then keep spaces
                                                  !  but skip all other
                                                  !  non-printing chars
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(INOUT) :: TEXT       ! Text to be modified
!
! !REMARKS:
!  The original "text" is destroyed upon exit.
!
! !REVISION HISTORY:
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER          :: ilen, iasc, icnt, i, Start
    CHARACTER(LEN=1) :: ch

    ! Default: Skip space characters
    Start = 32

    ! If KEEPSPACES=T then skip all non-printing characters,
    ! but keep space characters. (bmy, 1/30/18)
    IF ( PRESENT( KeepSpaces ) ) THEN
       IF ( KeepSpaces ) Start = 31
    ENDIF

    ilen = LEN(text)
    IF (ilen.GT.1) THEN
       icnt = 1
       DO i=1,ilen
          iasc = ICHAR(text(i:i))

          ! Keep characters between these limits
          IF ( ( iasc > Start ).AND. (iasc < 255 ) ) THEN
             ch = text(i:i)
             text(icnt:icnt) = ch
             icnt = icnt + 1
          ENDIF
       ENDDO
       ! Fill remainder of text with blanks
       DO i=icnt,ilen
          text(i:i) = ' '
       ENDDO
    ENDIF

  END SUBROUTINE CStrip
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsDigit
!
! !DESCRIPTION: Returned as true if ch is a numeric character (i.e., one of
!  the numbers from 0 to 9).
!\\
!\\
! !INTERFACE:
!
  FUNCTION IsDigit( ch ) RESULT( lnum )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=1), INTENT(IN) :: ch
!
! !RETURN VALUE:
!
    LOGICAL                      :: lnum
!
! !REMARKS:
!  NOTE: Changed name from ISNUM to ISDIGIT (bmy, 7/15/04)
!
! !REVISION HISTORY:
!     DATE:   NOV. 11, 1993
!     AUTHOR: R.D. STEWART
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER iasc

    iasc = ICHAR(ch)
    lnum = .FALSE.
    IF ((iasc.GE.48).AND.(iasc.LE.57)) THEN
       lnum = .TRUE.
    ENDIF

  END FUNCTION IsDigit
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: StrRepl
!
! !DESCRIPTION: Subroutine StrRepl replaces all instances of PATTERN within
!  a string STR with replacement text REPLTXT.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE StrRepl( Str, Pattern, ReplTxt )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: Pattern   ! Pattern to search for
    CHARACTER(LEN=*), INTENT(IN)    :: ReplTxt   ! Text to replace
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(INOUT) :: Str       ! String to be manipulated
!
! !REMARKS:
!  PATTERN and REPLTXT can now have a different number of characters.
!
! !REVISION HISTORY:
!  25 Jun 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Local variables
    INTEGER :: I1, I2

    !=================================================================
    ! StrRepl begins here!
    !=================================================================
    DO

       ! I1 is the first character that matches the search pattern;
       ! it must be 1 or larger.  Otherwise exit the routine.
       I1 = INDEX( Str, Pattern )
       IF ( I1 < 1 ) RETURN

       ! Replace the text.  I2 is the starting position of the
       ! string following the point of text replacement.
       I2 = I1 + LEN( Pattern )
       Str = Str(1:I1-1) // ReplTxt // Str(I2:)

    ENDDO

  END SUBROUTINE StrRepl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: StrSplit
!
! !DESCRIPTION: Subroutine STRSPLIT returns substrings in a string, separated
!  by a separator character (similar to IDL's StrSplit function).  This is
!  mainly a convenience wrapper for CHARPAK routine TxtExt.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE StrSplit( Str, Sep, Result, N_SubStrs )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: STR           ! String to be searched
    CHARACTER(LEN=1), INTENT(IN)  :: SEP           ! Separator character
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(OUT) :: Result(255)   ! Returned substrings
    INTEGER,          OPTIONAL    :: N_SubStrs     ! # of substrings
!
! !REVISION HISTORY:
!  11 Jul 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, IFLAG, COL
    CHARACTER(LEN=2047) :: WORD

    !=======================================================================
    ! STRSPLIT begins here!
    !=======================================================================

    ! Initialize
    I         = 0
    COL       = 1
    IFLAG     = 0
    RESULT(:) = ''

    ! Loop until all matches found, or end of string
    DO WHILE ( IFLAG == 0 )

       ! Look for strings beteeen separator string
       CALL TXTEXT ( SEP, TRIM( STR ), COL, WORD, IFLAG )

       ! Store substrings in RESULT array
       I         = I + 1
       RESULT(I) = TRIM( WORD )

    ENDDO

    ! Optional argument: return # of substrings found
    IF ( PRESENT( N_SUBSTRS ) ) N_SUBSTRS = I

  END SUBROUTINE StrSplit
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: StrSqueeze
!
! !DESCRIPTION: Subroutine STRSQUEEZE strips white space from both ends of a
!  string.  White space in the middle of the string (i.e. between characters)
!  will be preserved as-is.  Somewhat similar (though not exactly) to IDL's
!  STRCOMPRESS function.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE StrSqueeze( Str )
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(INOUT) :: Str   ! String to be squeezed
!
! !REVISION HISTORY:
!  11 Jul 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! STRSQUEEZE begins here!
    !=================================================================
    Str = ADJUSTR( TRIM( Str ) )
    Str = ADJUSTL( TRIM( Str ) )

  END SUBROUTINE StrSqueeze
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TranLc
!
! !DESCRIPTION: Tranlate a character variable to all lowercase letters.
!               Non-alphabetic characters are not affected.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TranLc( text )
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*) :: text
!
! !REMARKS:
!  The original "text" is destroyed.
!
! !REVISION HISTORY:
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: iasc,i,ilen

    ilen = LEN(text)
    DO I=1,ilen
       iasc = ICHAR(text(i:i))
       IF ((iasc.GT.64).AND.(iasc.LT.91)) THEN
          text(i:i) = CHAR(iasc+32)
       ENDIF
    ENDDO

  END SUBROUTINE TRANLC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TranUc
!
! !DESCRIPTION: Tranlate a character variable to all upper case letters.
!               Non-alphabetic characters are not affected.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TranUc( text )
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*) :: text
!
! !REMARKS:
!  The original "text" is destroyed.
!
! !REVISION HISTORY:
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: iasc,i,ilen

    ilen = LEN(text)
    DO i=1,ilen
       iasc = ICHAR(text(i:i))
       IF ((iasc.GT.96).AND.(iasc.LT.123)) THEN
          text(i:i) = CHAR(iasc-32)
       ENDIF
    ENDDO

  END SUBROUTINE TRANUC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TxtExt
!
! !DESCRIPTION: TxtExt extracts a sequence of characters from
!               text and transfers them to word.  The extraction
!               procedure uses a set of character "delimiters"
!               to denote the desired sequence of characters.
!               For example if ch=' ', the first character sequence
!               bracketed by blank spaces will be returned in word.
!               The extraction procedure begins in column, col,
!               of TEXT.  If text(col:col) = ch (any character in
!               the character string), the text is returned beginning
!               with col+1 in text (i.e., the first match with ch
!               is ignored).
!\\
!\\
!               After completing the extraction, col is incremented to
!               the location of the first character following the
!               end of the extracted text.
!\\
!\\
!               A status flag is also returned with the following
!               meaning(s)
!\\
!\\
!               IF iflg = -1, found a text block, but no more characters
!                             are available in TEXT
!                  iflg = 0,  task completed sucessfully (normal term)
!                  iflg = 1,  ran out of text before finding a block of
!                             text.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TxtExt(ch,text,col,word,iflg)
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: ch,text
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: col
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(OUT)   :: word
    INTEGER                         :: iflg
!
! !REMARKS:
!  TxtExt is short for Text Extraction.  This routine provides a set of
!  powerful line-by-line text search and extraction capabilities in
!  standard FORTRAN.
!
! !REVISION HISTORY:
!      AUTHOR: Robert D. Stewart
!        DATE: Jan. 1st, 1995
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: Tmax,T1,T2,imat
    LOGICAL :: again,prev

    ! Length of text
    Tmax = LEN(text)

    ! Fill Word with blanks
    WORD = REPEAT( ' ', LEN( WORD ) )

    IF (col.GT.Tmax) THEN
       ! Text does not contain any characters past Tmax.
       ! Reset col to one and return flag = {error condition}
       iflg = 1
       col = 1
    ELSEIF (col.EQ.Tmax) THEN
       ! End of TEXT reached
       CALL CntMat(ch,text(Tmax:Tmax),imat)
       IF (imat.EQ.0) THEN
          ! Copy character into Word and set col=1
          CALL CopyTxt(1,Text(Tmax:Tmax),Word)
          col = 1
          iflg = -1
       ELSE
          ! Same error condition as if col.GT.Tmax
          iflg = 1
       ENDIF
    ELSE
       ! Make sure column is not less than 1
       IF (col.LT.1) col=1
       CALL CntMat(ch,text(col:col),imat)
       IF (imat.GT.0) THEN
          prev=.true.
       ELSE
          prev=.false.
       ENDIF
       T1=col
       T2 = T1

       again = .true.
       DO WHILE (again)
          ! Check for a match with a character in ch
          CALL CntMat(ch,text(T2:T2),imat)
          IF (imat.GT.0) THEN
             ! Current character in TEXT matches one (or more) of the
             ! characters in ch.
             IF (prev) THEN
                IF (T2.LT.Tmax) THEN
                   ! Keep searching for a block of text
                   T2=T2+1
                   T1=T2
                ELSE
                   ! Did not find any text blocks before running
                   ! out of characters in TEXT.
                   again=.false.
                   iflg=1
                ENDIF
             ELSE
                 ! Previous character did not match ch, so terminate.
                 ! NOTE: This is "NORMAL" termination of the loop
                again=.false.
                T2=T2-1
                iflg = 0
             ENDIF
          ELSEIF (T2.LT.Tmax) THEN
             ! Add a letter to the current block of text
             prev = .false.
             T2=T2+1
          ELSE
             ! Reached the end of the characters in TEXT before reaching
             ! another delimiting character.  A text block was identified
             ! however.
             again=.false.
             iflg=-1
          ENDIF
       ENDDO

       IF (iflg.EQ.0) THEN
          ! Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = T2+1
       ELSE
          ! Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = 1
       ENDIF
    ENDIF

  END SUBROUTINE TxtExt
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: To_Uppercase
!
! !DESCRIPTION: Converts a string to uppercase, so that we can reliably
!  do string matching.
!\\
!\\
! !INTERFACE:
!
  FUNCTION To_UpperCase( Text ) RESULT( UpCaseText )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: Text         ! Input test
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255)           :: UpCaseText   ! Output text, uppercase
!
! !REMARKS:
!  Code originally from routine TRANUC (Author: R. D. Stewart, 19 May 1992)
!
! !REVISION HISTORY:
!  26 Jun 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: C, Ascii

    !=======================================================================
    ! Convert to uppercase
    !=======================================================================

    ! Initialize
    UpCaseText = Text

    ! Loop over all characters
    DO C = 1, LEN_TRIM( UpCaseText )

       ! Get the ASCII code for each character
       Ascii = ICHAR( UpCaseText(C:C) )

       ! If lowercase, convert to uppercase
       IF ( Ascii > 96 .and. Ascii < 123 ) THEN
          UpCaseText(C:C) = CHAR( Ascii - 32 )
       ENDIF
    ENDDO

  END FUNCTION To_UpperCase
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadOneLine
!
! !DESCRIPTION: Subroutine READ\_ONE\_LINE reads a line from the input file.
!  If the global variable VERBOSE is set, the line will be printed to stdout.
!  READ\_ONE\_LINE can trap an unexpected EOF if LOCATION is passed.
!  Otherwise, it will pass a logical flag back to the calling routine,
!  where the error trapping will be done.
!\\
!\\
! !INTERFACE:
!
  FUNCTION ReadOneLine( fId, EndOfFile, IoStatus, Squeeze ) RESULT( Line )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)      :: fId        ! File unit number
    LOGICAL, OPTIONAL        :: Squeeze    ! Call Strsqueeze?
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT)     :: EndOfFile  ! Denotes EOF condition
    INTEGER, INTENT(OUT)     :: IoStatus   ! I/O status code
!
! !RETURN VALUE:
!
    CHARACTER(LEN=MAXSTRLEN) :: Line       ! Single line from the input file
!
! !REMARKS:
!  Mostly used by routines in the History/ folder.
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version, based on GEOS-Chem
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Initialize
    !=================================================================
    EndOfFile = .FALSE.
    IoStatus  = 0
    Line      = ''

    !=================================================================
    ! Read data from the file
    !=================================================================

    ! Read a line from the file
    READ( fId, '(a)', IOSTAT=IoStatus ) Line

    ! IO Status < 0: EOF condition
    IF ( IoStatus < 0 ) THEN
       EndOfFile = .TRUE.
       RETURN
    ENDIF

    ! If desired, call StrSqueeze to strip leading and trailing blanks
    IF ( PRESENT( Squeeze ) ) THEN
       IF ( Squeeze ) THEN
          CALL StrSqueeze( Line )
       ENDIF
    ENDIF

  END FUNCTION ReadOneLine
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CleanText
!
! !DESCRIPTION: Strips commas, apostrophes, spaces, and tabs from a string.
!\\
!\\
! !INTERFACE:
!
  FUNCTION CleanText( Str ) RESULT( CleanStr )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: Str        ! Original string
!
! !RETURN VALUE
!
    CHARACTER(LEN=255)           :: CleanStr   ! Cleaned-up string
!
! !REMARKS:
!  Mostly used by routines in the History/ folder.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialize
    CleanStr = Str

    ! Strip out non-printing characters (e.g. tabs)
    CALL CStrip    ( CleanStr           )

    ! Remove commas and quotes
    CALL StrRepl   ( CleanStr, ",", " " )
    CALL StrRepl   ( CleanStr, "'", " " )

    ! Remove leading and trailing spaces
    CALL StrSqueeze( CleanStr           )

  END FUNCTION CleanText
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: WordWrapPrint
!
! !DESCRIPTION: Prints a text string wrapped to a specified line width.
!  Useful for displaying error and warning messages.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WordWrapPrint( Text, LineWidth, Delimiter )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: Text        ! Text to print
    INTEGER,          INTENT(IN) :: LineWidth   ! Width (characters) of lines
    CHARACTER(LEN=1), OPTIONAL   :: Delimiter   ! Delimiter between words
!
! !REMARKS:
!  The default DELIMITER is the space (" ") character.
!
! !REVISION HISTORY:
!  20 Dec 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER          :: C, S, B, Matches, Length

    ! Arrays
    INTEGER          :: BreakPts(100)
    INTEGER          :: SpaceLoc(500)

    ! Strings
    CHARACTER(LEN=1) :: Delim

    !=======================================================================
    ! WordWrapPrint begins here!
    !=======================================================================

    ! SpaceLoc is the array of where delimiters (usually the " "
    ! character) occur in the text, and S is its index.
    S           = 1
    SpaceLoc    = 0

    ! BreakPts is the array of where line breaks occur
    ! and B is its index.
    BreakPts    = 0
    B           = 1
    BreakPts(B) = 1

    ! Delimiter for separating words (will be the space character by default)
    IF ( PRESENT( Delimiter ) ) THEN
       Delim = Delimiter
    ELSE
       Delim = ' '
    ENDIF

    ! Find the Location of spaces in the text
    CALL CntMat( Text, ' ', Matches, SpaceLoc )

    ! Loop through the number of matches
    DO

       ! Move to the next delimiter location
       S = S + 1

       ! Compute the length of the line
       Length = SpaceLoc(S) - BreakPts(B)

       ! If the length of this segment is greater than the requested
       ! line length, store the position of this line break
       IF ( Length > LineWidth ) THEN
          B           = B             + 1
          BreakPts(B) = SpaceLoc(S-1) + 1
       ENDIF

       ! If we have exceeded the number of delimiters in the text, then set
       ! the last breakpoint at the end of the text and exit the loop.
       IF ( S > Matches ) THEN
          B           = B + 1
          BreakPts(B) = LEN_TRIM( Text ) + 1
          EXIT
       ENDIF

    ENDDO

    ! Print each line
    DO C = 1, B-1
       WRITE( 6, '(a)' ) Text( BreakPts(C):BreakPts(C+1)-1 )
    ENDDO

  END SUBROUTINE WordWrapPrint
!EOC
END MODULE CharPak_Mod
