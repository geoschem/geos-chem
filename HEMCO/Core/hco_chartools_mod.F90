!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!!BOP
!
! !MODULE: hco_chartools_mod.F90
!
! !DESCRIPTION: Module HCO\_CHARTOOLS\_MOD contains a collection of 
! helper routines to handle character strings and parse tokens. It also 
! contains definitions of special characters, including space, tab, wildcard, 
! separator, colon, and comment. Some of these characters (wildcard, separator, 
! colon, as well as tokens ROOT, MET, RES) may be manually set in the settings 
! section of the HEMCO configuration file.
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
  PUBLIC :: IsInWord
  PUBLIC :: NextCharPos
  PUBLIC :: HCO_WCD
  PUBLIC :: HCO_SPC
  PUBLIC :: HCO_SEP
  PUBLIC :: HCO_COL
  PUBLIC :: HCO_CMT
  PUBLIC :: HCO_TAB
  PUBLIC :: HCO_ROOTTOKEN
  PUBLIC :: HCO_METTOKEN
  PUBLIC :: HCO_RESTOKEN
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller   - Initialization
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jul 2014 - C. Keller   - Characters added
!  17 Oct 2014 - C. Keller   - Added parser and tokens ROOT, MET, and RES.
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
  CHARACTER(LEN=1), PARAMETER :: DEF_SPACE     = ' '
  CHARACTER(LEN=1), PARAMETER :: DEF_TAB       = ACHAR(9)
  CHARACTER(LEN=1), PARAMETER :: DEF_COMMENT   = '#'

  ! Default values for characters that can be changed 
  ! through the configuration file
  CHARACTER(LEN=1), PARAMETER :: DEF_COLON     = ':'
  CHARACTER(LEN=1), PARAMETER :: DEF_SEPARATOR = '/'
  CHARACTER(LEN=1), PARAMETER :: DEF_WILDCARD  = '*'

  !---------------------------------------------------------------------------
  ! Default tokens
  ! HEMCO has three tokens that can be specified in the HEMCO configuration
  ! file: ROOT (root directory), MET (met fields), and RES (horizontal 
  ! resolution). These tokens can be used in file names to be dynamically
  ! replaced, e.g. file.$MET.$RES.nc becomes file.geos5.4x5.nc if MET is set
  ! to 'geos5' and RES to '4x5'. Tokens are also allowed for dates ($YYYY,
  ! $MM, $DD, $HH, see routine HCO_CharParse).
  ! The default tokens below will be used if by default, i.e. if the
  ! corresponding token is not specified in the HEMCO configuration file.
  !---------------------------------------------------------------------------

  ! Default root directory
  CHARACTER(LEN=1023), PARAMETER :: DEF_ROOT = '/path/do/dir'

  ! Default met field token
#if defined( GEOS_FP )
  CHARACTER(LEN=15),   PARAMETER :: DEF_MET = 'geosfp'
#elif defined( GEOS_5 )
  CHARACTER(LEN=15),   PARAMETER :: DEF_MET = 'geos5'
#elif defined( GEOS_4 )
  CHARACTER(LEN=15),   PARAMETER :: DEF_MET = 'geos4'
#elif defined( MERRA )
  CHARACTER(LEN=15),   PARAMETER :: DEF_MET = 'merra'
#elif defined( GCAP )
  CHARACTER(LEN=15),   PARAMETER :: DEF_MET = 'gcap'
#else
  CHARACTER(LEN=15),   PARAMETER :: DEF_MET = 'unknown_model'
#endif 

  ! Default resolution token
#if defined( GRID4x5 )
  CHARACTER(LEN=15),   PARAMETER :: DEF_RES = '4x5'
#elif defined( GRID2x25 )
  CHARACTER(LEN=15),   PARAMETER :: DEF_RES = '2x25'
#elif defined( GRID1x125 )
  CHARACTER(LEN=15),   PARAMETER :: DEF_RES = '1x125'
#elif defined( GRID05x0666 )
  CHARACTER(LEN=15),   PARAMETER :: DEF_RES = '05x0666'
#else
  CHARACTER(LEN=15),   PARAMETER :: DEF_RES = 'unknown_res'
#endif
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
    CALL STRSPLIT( CharStr, TRIM(SEP), SUBSTR, N )
    IF ( N > SIZE(Reals,1) ) THEN
       CALL HCO_ERROR( 'Too many substrings!', RC, THISLOC=LOC )
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
    CALL STRSPLIT( CharStr, TRIM(SEP), SUBSTR, N )
    IF ( N > SIZE(Reals,1) ) THEN
       CALL HCO_ERROR( 'Too many substrings!', RC, THISLOC=LOC )
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

    ! Extract strings to be translated into integers 
    CALL STRSPLIT( CharStr, TRIM(SEP), SUBSTR, N )
    IF ( N > SIZE(Ints,1) ) THEN
       CALL HCO_ERROR( 'Too many substrings!', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Return here if no entry found
    IF ( N == 0 ) RETURN 

    ! Pass all extracted strings to integer vector. Replace wildcard
    ! character with -999!
    DO I = 1, N
       IF ( TRIM(SUBSTR(I)) == TRIM(WC) ) THEN
          Ints(I) = -999
       ELSEIF ( TRIM(SUBSTR(I)) == '-' ) THEN
          Ints(I) = -999
       ELSE
          READ( SUBSTR(I), * ) Ints(I) 
       ENDIF
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
! by searching for tokens such as $ROOT, $YYYY, etc., within the string and 
! replacing those values by the intendend characters. 
!\\
!\\
! At the moment, the following tokens are searched and replaced:
!\begin{itemize}
!\item \$ROOT: will be replaced by the root path specified in the settings
! section of the configuration file.
!\item \$YYYY: will be replaced by the (4-digit) year according to the 
! source time settings set in the configuration file.
!\item \$MM: will be replaced by the (2-digit) month according to the
! source time settings set in the configuration file.
!\item \$DD: will be replaced by the (2-digit) day according to the
! source time settings set in the configuration file.
!\item \$HH: will be replaced by the (2-digit) hour according to the
! source time settings set in the configuration file.
!\end{itemize}
! !INTERFACE:
!
  SUBROUTINE HCO_CharParse ( str, yyyy, mm, dd, hh, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   )  :: yyyy  ! replace $YYYY with this value 
    INTEGER,          INTENT(IN   )  :: mm    ! replace $MM with this value 
    INTEGER,          INTENT(IN   )  :: dd    ! replace $DD with this value
    INTEGER,          INTENT(IN   )  :: hh    ! replace $HH with this value
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
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = 'HCO_CharParse (HCO_CharTools_Mod.F90)'
    CHARACTER(LEN=2047) :: TMPSTR, BEFORE, AFTER
    INTEGER             :: LN, IDX, OFF
    CHARACTER(LEN=4)    :: str4
    CHARACTER(LEN=2)    :: str2

    !=================================================================
    ! HCO_CharParse begins here
    !=================================================================

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
       str = TRIM(BEFORE) // TRIM(HCO_ROOTTOKEN()) // TRIM(AFTER)
    ENDIF

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

    ! Check for met. model token 
    !-------------------------------------------------------------------
    DO
       IDX = INDEX( str, '$MET' )
       IF ( IDX <= 0 ) EXIT 
       LN = LEN(str)
       IF ( IDX > 1 ) THEN
          BEFORE = str(1:(IDX-1))
       ELSE
          BEFORE = ''
       ENDIF
       OFF   = 4
       AFTER = str((IDX+OFF):LN)

       ! Updated string
       str = TRIM(BEFORE) // TRIM(HCO_METTOKEN()) // TRIM(AFTER)
    ENDDO

    ! Check for met. resolution token 
    !-------------------------------------------------------------------
    DO
       IDX = INDEX( str, '$RES' )
       IF ( IDX <= 0 ) EXIT 
       LN = LEN(str)
       IF ( IDX > 1 ) THEN
          BEFORE = str(1:(IDX-1))
       ELSE
          BEFORE = ''
       ENDIF
       OFF   = 4
       AFTER = str((IDX+OFF):LN)

       ! Updated string
       str = TRIM(BEFORE) // TRIM(HCO_RESTOKEN()) // TRIM(AFTER)
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
! !ROUTINE: HCO_WCD
!
! !DESCRIPTION: Function HCO\_WCD returns the HEMCO WILDCARD character. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_WCD() RESULT( WILDCARD )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !ARGUMENTS:
!
    CHARACTER(LEN=1) :: WILDCARD 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
    LOGICAL,          SAVE      :: FIRST = .TRUE.
    CHARACTER(LEN=1), SAVE      :: WCD
    LOGICAL                     :: FOUND
    INTEGER                     :: myRC

    !======================================================================
    ! HCO_WCD begins here 
    !======================================================================

    ! On first call, check if WildCard character has been set in settings.
    ! Use default value otherwise.
    IF ( FIRST ) THEN 
       CALL GetExtOpt( 0, 'Wildcard', OptValChar=WCD, Found=FOUND, RC=myRC )
       IF ( .NOT. FOUND .OR. myRC /= HCO_SUCCESS ) WCD = DEF_WILDCARD 
       FIRST = .FALSE.
    ENDIF

    ! Return
    WILDCARD = WCD

  END FUNCTION HCO_WCD
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_SPC
!
! !DESCRIPTION: Function HCO\_SPC returns the HEMCO space character. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_SPC() RESULT( SPACE )
!
! !ARGUMENTS:
!
    CHARACTER(LEN=1) :: SPACE 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_SPC begins here 
    !======================================================================

    ! Return
    SPACE = DEF_SPACE 

  END FUNCTION HCO_SPC
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_SEP
!
! !DESCRIPTION: Function HCO\_SEP returns the HEMCO SEPARATOR character. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_SEP() RESULT( SEPARATOR )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !ARGUMENTS:
!
    CHARACTER(LEN=1) :: SEPARATOR 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
    LOGICAL,          SAVE      :: FIRST = .TRUE.
    CHARACTER(LEN=1), SAVE      :: SEP
    LOGICAL                     :: FOUND
    INTEGER                     :: myRC

    !======================================================================
    ! HCO_SEP begins here 
    !======================================================================

    ! On first call, check if Separator character has been set in settings.
    ! Use default value otherwise.
    IF ( FIRST ) THEN 
       CALL GetExtOpt( 0, 'Separator', OptValChar=SEP, Found=FOUND, RC=myRC )
       IF ( .NOT. FOUND .OR. myRC /= HCO_SUCCESS ) SEP = DEF_SEPARATOR
       FIRST = .FALSE.
    ENDIF

    ! Return wildcard character
    SEPARATOR = SEP

  END FUNCTION HCO_SEP
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_COL
!
! !DESCRIPTION: Function HCO\_COL returns the HEMCO COLON character. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_COL() RESULT( COLON )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !ARGUMENTS:
!
    CHARACTER(LEN=1) :: COLON 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
    LOGICAL,          SAVE      :: FIRST = .TRUE.
    CHARACTER(LEN=1), SAVE      :: COL
    LOGICAL                     :: FOUND
    INTEGER                     :: myRC

    !======================================================================
    ! HCO_COL begins here 
    !======================================================================

    ! On first call, check if Colon character has been set in settings.
    ! Use default value otherwise.
    IF ( FIRST ) THEN 
       CALL GetExtOpt( 0, 'Colon', OptValChar=COL, Found=FOUND, RC=myRC )
       IF ( .NOT. FOUND .OR. myRC /= HCO_SUCCESS ) COL = DEF_COLON
       FIRST = .FALSE.
    ENDIF

    ! Return 
    COLON = COL

  END FUNCTION HCO_COL
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_ROOTTOKEN
!
! !DESCRIPTION: Function HCO\_ROOTTOKEN returns the HEMCO root character
! as specified in the HEMCO configuration file settings (ROOT: /set/path). 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_ROOTTOKEN() RESULT( ROOTOUT )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !ARGUMENTS:
!
    CHARACTER(LEN=2047) :: ROOTOUT
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
    LOGICAL,             SAVE   :: FIRST = .TRUE.
    CHARACTER(LEN=1023), SAVE   :: ROOT 
    LOGICAL                     :: FOUND
    INTEGER                     :: myRC

    !======================================================================
    ! HCO_ROOTTOKEN begins here 
    !======================================================================

    ! On first call, check if Colon character has been set in settings.
    ! Use default value otherwise.
    IF ( FIRST ) THEN 
       CALL GetExtOpt( 0, 'ROOT', OptValChar=ROOT, Found=FOUND, RC=myRC )
       IF ( .NOT. FOUND .OR. myRC /= HCO_SUCCESS ) ROOT = DEF_ROOT
       FIRST = .FALSE.
    ENDIF

    ! Return 
    ROOTOUT = ROOT

  END FUNCTION HCO_ROOTTOKEN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_METTOKEN
!
! !DESCRIPTION: Function HCO\_METTOKEN returns the HEMCO met field character
! (e.g. GEOS_FP) specified in the HEMCO configuration file settings (e.g. 
! MET: GEOS_FP). If not set in the HEMCO config. file, a default value is 
! taken based on the compiler switches.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_METTOKEN() RESULT( METOUT )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !ARGUMENTS:
!
    CHARACTER(LEN=15) :: METOUT
!
! !REVISION HISTORY:
!  17 Oct 2014 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
    LOGICAL,           SAVE   :: FIRST = .TRUE.
    CHARACTER(LEN=15), SAVE   :: MET
    LOGICAL                   :: FOUND
    INTEGER                   :: myRC

    !======================================================================
    ! HCO_METTOKEN begins here 
    !======================================================================

    ! On first call, check if Colon character has been set in settings.
    ! Use default value otherwise.
    IF ( FIRST ) THEN 
       CALL GetExtOpt( 0, 'MET', OptValChar=MET, Found=FOUND, RC=myRC )
       IF ( .NOT. FOUND .OR. myRC /= HCO_SUCCESS ) MET = DEF_MET
       FIRST = .FALSE.
    ENDIF

    ! Return 
    METOUT = MET

  END FUNCTION HCO_METTOKEN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_RESTOKEN
!
! !DESCRIPTION: Function HCO\_RESTOKEN returns the HEMCO resolution character
! (e.g. 2x25) specified in the HEMCO configuration file settings (e.g. 
! RES: 2x25). If not set in the HEMCO config. file, a default value is 
! taken based on the compiler switches.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_RESTOKEN() RESULT( RESOUT )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !ARGUMENTS:
!
    CHARACTER(LEN=15) :: RESOUT
!
! !REVISION HISTORY:
!  17 Oct 2014 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
    LOGICAL,           SAVE   :: FIRST = .TRUE.
    CHARACTER(LEN=15), SAVE   :: RES
    LOGICAL                   :: FOUND
    INTEGER                   :: myRC

    !======================================================================
    ! HCO_RESTOKEN begins here 
    !======================================================================

    ! On first call, check if Colon character has been set in settings.
    ! Use default value otherwise.
    IF ( FIRST ) THEN 
       CALL GetExtOpt( 0, 'RES', OptValChar=RES, Found=FOUND, RC=myRC )
       IF ( .NOT. FOUND .OR. myRC /= HCO_SUCCESS ) RES = DEF_RES
       FIRST = .FALSE.
    ENDIF

    ! Return 
    RESOUT = RES

  END FUNCTION HCO_RESTOKEN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_CMT
!
! !DESCRIPTION: Function HCO\_CMT returns the HEMCO COMMENT character. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_CMT() RESULT( COMMENT )
!
! !ARGUMENTS:
!
    CHARACTER(LEN=1) :: COMMENT 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_CMT begins here 
    !======================================================================

    ! Return wildcard character
    COMMENT = DEF_COMMENT 

  END FUNCTION HCO_CMT
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_TAB
!
! !DESCRIPTION: Function HCO\_TAB returns the HEMCO TAB character. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_TAB() RESULT( TAB )
!
! !ARGUMENTS:
!
    CHARACTER(LEN=1) :: TAB
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_TAB begins here 
    !======================================================================

    ! Return
    TAB = DEF_TAB

  END FUNCTION HCO_TAB
!EOC
END MODULE HCO_CharTools_Mod
