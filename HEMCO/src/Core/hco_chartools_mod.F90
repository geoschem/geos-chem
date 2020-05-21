!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!!BOP
!
! !MODULE: hco_chartools_mod.F90
!
! !DESCRIPTION: Module HCO\_CHARTOOLS\_MOD contains a collection of
! helper routines to handle character strings and parse tokens. It also
! contains definitions of special characters such as space, tab and
! comment.
! \\
! !INTERFACE:
!
MODULE HCO_CharTools_Mod
!
! !USES:
!
  USE HCO_Error_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCO_CharSplit
  PUBLIC :: HCO_CharMatch
  PUBLIC :: HCO_CharParse
  PUBLIC :: HCO_GetBase
  PUBLIC :: IsInWord
  PUBLIC :: NextCharPos
  PUBLIC :: GetNextLine
  PUBLIC :: HCO_READLINE
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller   - Initialization
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jul 2014 - C. Keller   - Characters added
!  17 Oct 2014 - C. Keller   - Added parser and tokens ROOT, MET, and RES.
!  12 Aug 2015 - R. Yantosca - Add new value of DEF_MET for MERRA2
!  12 Aug 2015 - R. Yantosca - Add new value of DEF_RES for 0.5 x 0.625 grids
!  20 Sep 2015 - C. Keller   - Moved tokens to hco_extlist_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE HCO_CharSplit
     MODULE PROCEDURE HCO_CharSplit_R8
     MODULE PROCEDURE HCO_CharSplit_R4
     MODULE PROCEDURE HCO_CharSplit_INT
  END INTERFACE
!
! !MODULE PARAMETER:
!
  ! Fixed characters
  CHARACTER(LEN=1), PARAMETER, PUBLIC :: HCO_SPC     = ' '
  CHARACTER(LEN=1), PARAMETER, PUBLIC :: HCO_TAB     = ACHAR(9)
  CHARACTER(LEN=1), PARAMETER, PUBLIC :: HCO_CMT     = '#'
!
! !PRIVATE MEMBER FUNCTIONS:
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CharSplit_R8
!
! !DESCRIPTION: Subroutine HCO\_CharSplit\_R8 splits the passed character
! string into N real8 values, using character SEP as separator. Wildcard
! values (WC) are set to -999.
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CharSplit_R8( CharStr, SEP, WC, Reals, N, RC )
!
! !USES:
!
    USE CharPak_Mod,  ONLY : StrSplit
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: CharStr   ! Character string
    CHARACTER(LEN=1), INTENT(IN   ) :: SEP       ! Separator
    CHARACTER(LEN=1), INTENT(IN   ) :: WC        ! Wildcard character
!
! !OUTPUT PARAMETERS:
!
    REAL(dp),         INTENT(  OUT) :: Reals(:)  ! Output values
    INTEGER,          INTENT(  OUT) :: N         ! # of valid values
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC         ! Return code
!
! !REVISION HISTORY:
!  18 Sep 2013 - C. Keller - Initial version (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    INTEGER            :: I
    CHARACTER(LEN=255) :: SUBSTR(255)
    CHARACTER(LEN=255) :: LOC

    !=================================================================
    ! HCO_CharSplit_R8 begins here!
    !=================================================================

    ! Enter
    LOC = 'HCO_CharSplit_R8 (HCO_CHARTOOLS_MOD.F90)'

    ! Init
    Reals(:) = -999_dp

    ! Extract strings to be translated into integers
    !CALL STRSPLIT( CharStr, TRIM(SEP), SUBSTR, N )
    CALL STRSPLIT( CharStr, SEP, SUBSTR, N )
    IF ( N > SIZE(Reals,1) ) THEN
       WRITE(*,*) 'Too many substrings - error in ', TRIM(LOC)
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Return here if no entry found
    IF ( N == 0 ) RETURN

    ! Pass all extracted strings to integer vector. Replace wildcard
    ! character with -999!
    DO I = 1, N
       IF ( TRIM(SUBSTR(I)) == TRIM(WC) ) THEN
          Reals(I) = -999_dp
       ELSEIF ( TRIM(SUBSTR(I)) == '-' ) THEN
          Reals(I) = -999_dp
       ELSE
          READ( SUBSTR(I), * ) Reals(I)
       ENDIF
    ENDDO

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_CharSplit_R8
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CharSplit_R4
!
! !DESCRIPTION: Subroutine HCO\_CharSplit\_R4 splits the passed character
! string into N real4 values, using character SEP as separator. Wildcard
! values (WC) are set to -999.
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CharSplit_R4( CharStr, SEP, WC, Reals, N, RC )
!
! !USES:
!
    USE CharPak_Mod,  ONLY : StrSplit
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: CharStr    ! Character string
    CHARACTER(LEN=1), INTENT(IN   ) :: SEP        ! Separator
    CHARACTER(LEN=1), INTENT(IN   ) :: WC         ! Wildcard character
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(  OUT) :: Reals(:)   ! Output values
    INTEGER,          INTENT(  OUT) :: N          ! # of valid values
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC         ! Return code
!
! !REVISION HISTORY:
!  18 Sep 2013 - C. Keller - Initial version (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I
    CHARACTER(LEN=255)    :: SUBSTR(255)
    CHARACTER(LEN=255)    :: LOC

    !=================================================================
    ! HCO_CharSplit_R4 begins here!
    !=================================================================

    ! Enter
    LOC = 'HCO_CharSplit_R4 (HCO_CHARTOOLS_MOD.F90)'

    ! Init
    Reals(:) = -999_sp

    ! Extract strings to be translated into integers
    !CALL STRSPLIT( CharStr, TRIM(SEP), SUBSTR, N )
    CALL STRSPLIT( CharStr, SEP, SUBSTR, N )
    IF ( N > SIZE(Reals,1) ) THEN
       WRITE(*,*) 'Too many substrings - error in ', TRIM(LOC)
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Return here if no entry found
    IF ( N == 0 ) RETURN

    ! Pass all extracted strings to integer vector. Replace wildcard
    ! character with -999!
    DO I = 1, N
       IF ( TRIM(SUBSTR(I)) == TRIM(WC) ) THEN
          Reals(I) = -999_sp
       ELSEIF ( TRIM(SUBSTR(I)) == '-' ) THEN
          Reals(I) = -999_sp
       ELSE
          READ( SUBSTR(I), * ) Reals(I)
       ENDIF
    ENDDO

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_CharSplit_R4
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CharSplit_Int
!
! !DESCRIPTION: Subroutine HCO\_CharSplit\_Int splits the passed character
! string into N integers, using character SEP as separator. Wildcard
! values (WC) are set to -999.
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CharSplit_INT( CharStr, SEP, WC, Ints, N, RC )
!
! !USES:
!
    USE CharPak_Mod, ONLY : StrSplit
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: CharStr    ! Character string
    CHARACTER(LEN=1), INTENT(IN   ) :: SEP        ! Separator
    CHARACTER(LEN=1), INTENT(IN   ) :: WC         ! Wildcard character
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT) :: Ints(:)    ! Output values
    INTEGER,          INTENT(  OUT) :: N          ! # of valid values
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC         ! Return code
!
! !REVISION HISTORY:
!  18 Sep 2013 - C. Keller - Initial version (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I
    CHARACTER(LEN=255)    :: SUBSTR(255)
    CHARACTER(LEN=255)    :: LOC

    !=================================================================
    ! HCO_CharSplit_INT begins here!
    !=================================================================

    ! Enter
    LOC = 'HCO_CharSplit_Int (HCO_CHARTOOLS_MOD.F90)'

    ! Init
    Ints(:) = -999

    ! If input string is wildcard or otherwise empty, return here.
    IF ( TRIM(CharStr) == TRIM(WC) .OR. &
         TRIM(CharStr) == '-'            ) THEN
       N = 0
       RETURN
    ENDIF

    ! Extract strings to be translated into integers
    !CALL STRSPLIT( CharStr, TRIM(SEP), SUBSTR, N )
    CALL STRSPLIT( CharStr, SEP, SUBSTR, N )
    IF ( N > SIZE(Ints,1) ) THEN
       WRITE(*,*) 'Too many substrings - error in ', TRIM(LOC)
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Return here if no entry found
    IF ( N == 0 ) RETURN

    ! Pass all extracted strings to integer vector.
    DO I = 1, N
       READ( SUBSTR(I), * ) Ints(I)
    ENDDO

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_CharSplit_INT
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CharMatch
!
! !DESCRIPTION: Subroutine HCO\_CharMatch returns the index of each
! vector element of vec1 in vec2. nnmatch denotes the number of
! vec1 elements which have a matching counterpart in vec2.
! For example, if vec1 is (/ 'NO', 'CO', 'ALK4', 'HBr' /), and
! vec2 is (/ 'CO', 'NO', 'CH3Br' /), then matchidx becomes
! (/ 2, 1, -1, -1 /) and nnmatch is 2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CharMatch( vec1, n1, vec2, n2, matchidx, nnmatch )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: vec1(n1)     ! char. vector 1
    INTEGER,          INTENT(IN   ) :: n1           ! len of vec1
    CHARACTER(LEN=*), INTENT(IN   ) :: vec2(n2)     ! char. vector 2
    INTEGER,          INTENT(IN   ) :: n2           ! len of vec2
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT) :: matchidx(n1) ! index of vec2 in vec1
    INTEGER,          INTENT(  OUT) :: nnmatch      ! # of matches
!
! !REVISION HISTORY:
!  18 Sep 2013 - C. Keller - Initial version (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J

    !=================================================================
    ! HCO_CharMatch begins here!
    !=================================================================

    ! Init
    nnmatch = 0

    ! Do for every element in vec1
    DO I = 1, n1

       ! Default = no match
       matchidx(I) = -1

       DO J = 1, n2
          IF ( TRIM(vec1(I)) == TRIM(vec2(J)) ) THEN
             matchidx(I) = J
             nnmatch     = nnmatch + 1
             EXIT
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE HCO_CharMatch
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CharParse
!
! !DESCRIPTION: Routine HCO\_CharParse parses the provided character string
! by searching for tokens such as \$ROOT, \$YYYY, etc., within the string and
! replacing those values by the intendend characters.
!\\
!\\
! The following list shows the 'default' HEMCO tokens. These are available
! in any HEMCO simulation. Tokens \$ROOT, \$MET, and \$RES are internally
! stored as a HEMCO option in module hco\_extlist\_mod.F90 (see subroutine
! HCO\_SetDefaultToken).
! \begin{itemize}
! \item \$ROOT: will be replaced by the root path specified in the settings
! section of the configuration file.
! \item \$MET: will be replaced by the met-field token.
! \item \$RES: will be replaced by the resolution token.
! \item \$YYYY: will be replaced by the (4-digit) year according to the
! source time settings set in the configuration file.
! \item \$MM: will be replaced by the (2-digit) month according to the
! source time settings set in the configuration file.
! \item \$DD: will be replaced by the (2-digit) day according to the
! source time settings set in the configuration file.
! \item \$HH: will be replaced by the (2-digit) hour according to the
! source time settings set in the configuration file.
! \item \$MN: will be replaced by the (2-digit) minute.
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE HCO_CharParse ( HcoConfig, str, yyyy, mm, dd, hh, mn, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod, ONLY  : HCO_GetOpt, HCO_Root
    USE HCO_Types_Mod,   ONLY  : ConfigObj
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER        :: HcoConfig
    INTEGER,          INTENT(IN   )  :: yyyy  ! replace $YYYY with this value
    INTEGER,          INTENT(IN   )  :: mm    ! replace $MM with this value
    INTEGER,          INTENT(IN   )  :: dd    ! replace $DD with this value
    INTEGER,          INTENT(IN   )  :: hh    ! replace $HH with this value
    INTEGER,          INTENT(IN   )  :: mn    ! replace $MN with this value
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT)  :: str   ! string to be parsed
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! return code
!
! !REVISION HISTORY:
!  01 Oct 2014 - C. Keller - Initial version
!  20 Sep 2015 - C. Keller - Tokens can now be any option setting set in the
!                            HEMCO configuration file.
!  07 Jul 2017 - C. Keller - Extended list of token delimiters.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = 'HCO_CharParse (HCO_CharTools_Mod.F90)'
    CHARACTER(LEN=255)  :: TOKEN
    CHARACTER(LEN=2047) :: TMPSTR, BEFORE, AFTER
    INTEGER             :: I, LN, IDX, OFF
    CHARACTER(LEN=4)    :: str4
    CHARACTER(LEN=2)    :: str2
    CHARACTER(LEN=1)    :: SEP

    !=================================================================
    ! HCO_CharParse begins here
    !=================================================================

    ! Get characters
    SEP = HCO_GetOpt(HcoConfig%ExtList,'Separator')

    ! Check for year token
    !-------------------------------------------------------------------
    DO
       IDX = INDEX( str, '$YYYY' )
       IF ( IDX <= 0 ) EXIT
       LN = LEN(str)
       IF ( IDX > 1 ) THEN
          BEFORE = str(1:(IDX-1))
       ELSE
          BEFORE = ''
       ENDIF
       OFF   = 5
       AFTER = str((IDX+OFF):LN)

       WRITE(str4,'(i4.4)') yyyy

       ! Updated string
       str = TRIM(BEFORE) // TRIM(str4) // TRIM(AFTER)
    ENDDO

    ! Check for month token
    !-------------------------------------------------------------------
    DO
       IDX = INDEX( str, '$MM' )
       IF ( IDX <= 0 ) EXIT
       LN = LEN(str)
       IF ( IDX > 1 ) THEN
          BEFORE = str(1:(IDX-1))
       ELSE
          BEFORE = ''
       ENDIF
       OFF   = 3
       AFTER = str((IDX+OFF):LN)

       WRITE(str2,'(i2.2)') mm

       ! Updated string
       str = TRIM(BEFORE) // TRIM(str2) // TRIM(AFTER)
    ENDDO

    ! Check for day token
    !-------------------------------------------------------------------
    DO
       IDX = INDEX( str, '$DD' )
       IF ( IDX <= 0 ) EXIT
       LN = LEN(str)
       IF ( IDX > 1 ) THEN
          BEFORE = str(1:(IDX-1))
       ELSE
          BEFORE = ''
       ENDIF
       OFF   = 3
       AFTER = str((IDX+OFF):LN)

       WRITE(str2,'(i2.2)') dd

       ! Updated string
       str = TRIM(BEFORE) // TRIM(str2) // TRIM(AFTER)
    ENDDO

    ! Check for hour token
    !-------------------------------------------------------------------
    DO
       IDX = INDEX( str, '$HH' )
       IF ( IDX <= 0 ) EXIT
       LN = LEN(str)
       IF ( IDX > 1 ) THEN
          BEFORE = str(1:(IDX-1))
       ELSE
          BEFORE = ''
       ENDIF
       OFF   = 3
       AFTER = str((IDX+OFF):LN)

       WRITE(str2,'(i2.2)') hh

       ! Updated string
       str = TRIM(BEFORE) // TRIM(str2) // TRIM(AFTER)
    ENDDO

    ! Check for minute token
    !-------------------------------------------------------------------
    DO
       IDX = INDEX( str, '$MN' )
       IF ( IDX <= 0 ) EXIT
       LN = LEN(str)
       IF ( IDX > 1 ) THEN
          BEFORE = str(1:(IDX-1))
       ELSE
          BEFORE = ''
       ENDIF
       OFF   = 3
       AFTER = str((IDX+OFF):LN)

       WRITE(str2,'(i2.2)') mn

       ! Updated string
       str = TRIM(BEFORE) // TRIM(str2) // TRIM(AFTER)
    ENDDO

    ! Check for root token
    !-------------------------------------------------------------------
    IDX = INDEX( str, '$ROOT' )
    IF ( IDX > 0 ) THEN
       LN = LEN(str)
       IF ( IDX > 1 ) THEN
          BEFORE = str(1:(IDX-1))
       ELSE
          BEFORE = ''
       ENDIF
       OFF   = 5
       AFTER = str((IDX+OFF):LN)

       ! Updated string
       str = TRIM(BEFORE) // TRIM(HCO_ROOT(HcoConfig)) // TRIM(AFTER)
    ENDIF

    ! Check for any other token
    !-------------------------------------------------------------------
    DO
       IDX = INDEX( str, '$' )
       IF ( IDX <= 0 ) EXIT
       LN = LEN(TRIM(str))

       ! Determine token name:
       ! Find end of token by starting at the first character after the
       ! token and then advance in string until a 'token end character'
       ! is encountered.
       DO I = IDX+1,LN

          ! Special case that end of string is encountered:
          IF ( I == LN ) THEN
             TOKEN = str( (IDX+1) : I )
             OFF   = I
             EXIT
          ENDIF

          ! Scan for token cap:
          IF ( str(I:I) == ' '  .OR. &
               str(I:I) == '.'  .OR. &
               str(I:I) == ':'  .OR. &
               str(I:I) == '$'  .OR. &
               str(I:I) == '%'  .OR. &
               str(I:I) == '+'  .OR. &
               str(I:I) == '*'  .OR. &
               str(I:I) == '/'  .OR. &
               str(I:I) == '^'  .OR. &
               str(I:I) == '_'  .OR. &
               str(I:I) == '*'  .OR. &
               str(I:I) == '/'  .OR. &
               str(I:I) == '^'  .OR. &
               str(I:I) == '-'  .OR. &
               str(I:I) == 'x'  .OR. &
               str(I:I) == '('  .OR. &
               str(I:I) == ')'  .OR. &
               str(I:I) == '0'  .OR. &
               str(I:I) == '1'  .OR. &
               str(I:I) == '2'  .OR. &
               str(I:I) == '3'  .OR. &
               str(I:I) == '4'  .OR. &
               str(I:I) == '5'  .OR. &
               str(I:I) == '6'  .OR. &
               str(I:I) == '7'  .OR. &
               str(I:I) == '8'  .OR. &
               str(I:I) == '9'  .OR. &
               str(I:I) == SEP        ) THEN

             TOKEN = str( (IDX+1) : (I-1) )
             OFF   = I
             EXIT
          ENDIF
       ENDDO

       IF ( IDX > 1 ) THEN
          BEFORE = str(1:(IDX-1))
       ELSE
          BEFORE = ''
       ENDIF
       IF ( OFF >= LN ) THEN
          AFTER = ''
       ELSE
          AFTER = str(OFF:LN)
       ENDIF

       ! Update string
       str = TRIM(BEFORE) // &
             TRIM(HCO_GetOpt(HcoConfig%ExtList,TOKEN)) // &
             TRIM(AFTER)

    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_CharParse
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GetBase
!
! !DESCRIPTION: Routine HCO\_GetBase returns the base location of the given
! file. This is the entire file path up to the last forward slash, e.g. for
! file '/home/dir/Config.rc', the base is '/home/dir/'
!
! !INTERFACE:
!
  SUBROUTINE HCO_GetBase ( str, base, RC )
!
! !USES:
!
    USE CharPak_Mod, ONLY : StrSplit
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )  :: str   ! string to be checked
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT)  :: base  ! base
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! return code
!
! !REVISION HISTORY:
!  16 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, N
    CHARACTER(LEN=255)  :: SUBSTR(255)

    !=================================================================
    ! HCO_GetBase begins here
    !=================================================================

    CALL STRSPLIT( str, '/', SUBSTR, N )
    IF ( N <= 1 ) THEN
       base = '.'
    ELSE
       base = '/' // TRIM(SUBSTR(1))
       IF ( N > 2 ) THEN
          DO I = 2,(N-1)
             base = TRIM(base) // '/' // TRIM(SUBSTR(I))
          ENDDO
       ENDIF
    ENDIF


    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_GetBase
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsInWord
!
! !DESCRIPTION: Function IsInWord checks if the word InString contains the
! sequence of SearchString.
!\\
! !INTERFACE:
!
  FUNCTION IsInWord( InString, SearchString ) RESULT ( Cnt )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: InString
    CHARACTER(LEN=*), INTENT(IN   ) :: SearchString
!
! !RETURN VALUE:
!
    LOGICAL                         :: Cnt
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------

    Cnt = INDEX( TRIM(InString), TRIM(SearchString) ) > 0

  END FUNCTION IsInWord
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NextCharPos
!
! !DESCRIPTION: Function NextCharPos returns the position of the next
! occurrence of character CHR in word WORD, starting from position START.
! Returns -1 if the word does not contain CHR at all (after position START).
!\\
!\\
! !INTERFACE:
!
  FUNCTION NextCharPos ( WORD, CHR, START ) RESULT ( POS )
!
! !USES:
!
!
! !INPUT ARGUMENTS:
!
    CHARACTER(LEN=*),           INTENT(IN)  :: WORD
    CHARACTER(LEN=1),           INTENT(IN)  :: CHR
    INTEGER,          OPTIONAL, INTENT(IN)  :: START
!
! !RETURN ARGUMENT:
!
    INTEGER                                 :: POS
!
! !REVISION HISTORY:
!  09 Jul 2014 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER                         :: LNG, N, BEG

    !=================================================================
    ! NextCharPos begins here
    !=================================================================

    ! Initialize
    POS = -1

    ! Get first index
    IF ( PRESENT(START) ) THEN
       BEG = START
    ELSE
       BEG = 1
    ENDIF

    ! Lenght of word
    LNG = LEN(TRIM(WORD))

    ! Error traps
    IF ( BEG > LNG ) RETURN

    ! Search for occurrence of CHR
    DO N = BEG, LNG
       IF ( WORD(N:N) == CHR ) THEN
          POS = N
          EXIT
       ENDIF
    ENDDO

  END FUNCTION NextCharPos
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetNextLine
!
! !DESCRIPTION: Subroutine GetNextLine returns the next line.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetNextLine( LUN, LINE, EOF, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: LUN         ! Stream to read from
!
! !OUTPUT PARAMETERS
!
    CHARACTER(LEN=*), INTENT(  OUT) :: LINE        ! Next (valid) line in stream
!
! !INPUT/OUTPUT PARAMETERS
!
    LOGICAL,          INTENT(INOUT) :: EOF         ! End of file encountered?
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  10 Apr 2015 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER             :: IOS
    CHARACTER(LEN=4095) :: DUM

    !=================================================================
    ! GetNextLine begins here
    !=================================================================

    ! Init
    RC = HCO_SUCCESS

    ! Repeat until valid line is encountered
    DO
       CALL HCO_ReadLine( LUN, DUM, EOF, RC )
       IF ( EOF .OR. RC /= HCO_SUCCESS ) RETURN

       ! Skip if empty or commented line
       IF ( TRIM(DUM) == ''      ) CYCLE
       IF ( DUM(1:1)  == HCO_CMT ) CYCLE

       ! Make sure that character string DUM is not longer than LINE
       IF ( LEN_TRIM(DUM) > LEN(LINE) ) THEN
          WRITE( 6, '(a)' ) REPEAT( '=', 79 )
          WRITE( 6, * ) ' Line is too long - cannot copy into output argument '
          WRITE( 6, * ) TRIM(DUM)
          WRITE( 6, * ) ' '
          WRITE( 6, * ) ' To fix this, increase length of argument `LINE` in '
          WRITE( 6, * ) ' the subprogram which is calling '
          WRITE( 6, * ) ' HCO_ReadLine (hco_chartools_mod.F90)'
          RC = HCO_FAIL
          RETURN
          WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       ELSE
          ! If we get here, exit loop
          LINE = DUM
          EXIT
       ENDIF

    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetNextLine
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ReadLine
!
! !DESCRIPTION: Subroutine HCO\_Line reads a line from the provided stream.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ReadLine( LUN, LINE, EOF, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: LUN      ! Stream LUN
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(INOUT) :: LINE     ! Line
    LOGICAL,          INTENT(INOUT) :: EOF      ! End of file?
    INTEGER,          INTENT(INOUT) :: RC       ! Return code
!
! !REVISION HISTORY:
!  18 Sep 2013 - C. Keller - Initial version (adapted from B. Yantosca's code)
!  15 Jul 2014 - R. Yantosca - Remove dependency on routine IOERROR
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: IOS
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=4095) :: DUM

    !=================================================================
    ! HCO_ReadLine begins here!
    !=================================================================

    ! Initialize
    EOF = .FALSE.
    RC  = HCO_SUCCESS

    ! Read a line from the file
    READ( LUN, '(a)', IOSTAT=IOS ) DUM

    ! IO Status < 0: EOF condition
    IF ( IOS < 0 ) THEN
       EOF = .TRUE.
       RETURN
    ENDIF

    ! IO Status > 0: true I/O error condition
    IF ( IOS > 0 ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, 100   ) IOS
100    FORMAT( 'ERROR ', i5, ' in HCO_Readline (hco_chartools_mod.F90)' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Make sure that character string DUM is not longer than LINE
    IF ( LEN(TRIM(DUM)) > LEN(LINE) ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, * ) ' Line is too long - cannot read line ', TRIM(DUM)
       WRITE( 6, * ) ' '
       WRITE( 6, * ) ' To fix this, increase length of argument `LINE` in '
       WRITE( 6, * ) ' HCO_ReadLine (hco_chartools_mod.F90)'
       RC = HCO_FAIL
       RETURN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ELSE
       LINE = DUM(1:LEN(LINE))
    ENDIF

  END SUBROUTINE HCO_ReadLine
!EOC
END MODULE HCO_CharTools_Mod
