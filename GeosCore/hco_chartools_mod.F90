!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!
! !MODULE: hco_chartools_mod 
!
! !DESCRIPTION: Module HCO\_CHARTOOLS\_MOD contains a collection of 
! helper routines to handle character strings. It also contains definitions
! of special characters, including space, tab, wildcard, separator, colon, 
! and comment. Some of these characters (wildcard, separator, colon) may be
! manually set in the settings section of the HEMCO configuration file.
! \\
! !INTERFACE: 
!
      MODULE HCO_CHARTOOLS_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: HCO_CharSplit
      PUBLIC :: HCO_CharMatch
      PUBLIC :: IsInWord
      PUBLIC :: NextCharPos
      PUBLIC :: HCO_WCD
      PUBLIC :: HCO_SPC
      PUBLIC :: HCO_SEP
      PUBLIC :: HCO_COL
      PUBLIC :: HCO_CMT
      PUBLIC :: HCO_TAB
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller - Initialization
!
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
!
! !PRIVATE MEMBER FUNCTIONS:
!
      CONTAINS
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CharSplit_R8 
!
! !DESCRIPTION: Subroutine HCO_CharSplit\_R8 splits the passed character
! string into N real8 values, using character SEP as separator. Wildcard
! values (WC) are set to -999. 
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_CharSplit_R8 ( CharStr, SEP, WC, Reals, N, RC ) 
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : STRSPLIT
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN   )       :: CharStr ! Character string
      CHARACTER(LEN=1), INTENT(IN   )       :: SEP     ! Separator
      CHARACTER(LEN=1), INTENT(IN   )       :: WC      ! Wildcard character
!
! !OUTPUT PARAMETERS:
!
      REAL*8,           INTENT(  OUT)       :: Reals(:) ! Output values
      INTEGER,          INTENT(  OUT)       :: N        ! # of valid values
      INTEGER,          INTENT(INOUT)       :: RC       ! Return code
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CharSplit_R4 
!
! !DESCRIPTION: Subroutine HCO_CharSplit\_R4 splits the passed character
! string into N real4 values, using character SEP as separator. Wildcard
! values (WC) are set to -999. 
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_CharSplit_R4 ( CharStr, SEP, WC, Reals, N, RC ) 
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : STRSPLIT
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN   )       :: CharStr 
      CHARACTER(LEN=1), INTENT(IN   )       :: SEP 
      CHARACTER(LEN=1), INTENT(IN   )       :: WC 
!
! !OUTPUT PARAMETERS:
!
      REAL*4,           INTENT(  OUT)       :: Reals(:)
      INTEGER,          INTENT(  OUT)       :: N
      INTEGER,          INTENT(INOUT)       :: RC
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CharSplit_INT 
!
! !DESCRIPTION: Subroutine HCO_CharSplit\_R8 splits the passed character
! string into N integers, using character SEP as separator. Wildcard
! values (WC) are set to -999. 
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_CharSplit_INT ( CharStr, SEP, WC, Ints, N, RC ) 
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : STRSPLIT
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN   )       :: CharStr 
      CHARACTER(LEN=1), INTENT(IN   )       :: SEP 
      CHARACTER(LEN=1), INTENT(IN   )       :: WC 
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(  OUT)       :: Ints(:)
      INTEGER,          INTENT(  OUT)       :: N
      INTEGER,          INTENT(INOUT)       :: RC
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CharMatch
!
! !DESCRIPTION: Subroutine HCO_CharMatch returns the index of each 
! vector element of vec1 in vec2. nnmatch denotes the number of
! vec1 elements which have a matching counterpart in vec2.
! For example, if vec1 is (/ 'NO', 'CO', 'ALK4', 'HBr' /), and 
! vec2 is (/ 'CO', 'NO', 'CH3Br' /), then matchidx becomes
! (/ 2, 1, -1, -1 /) and nnmatch is 2. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_CharMatch (vec1, n1, vec2, n2, matchidx, nnmatch )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN   )       :: vec1(n1)     ! char. vector 1 
      INTEGER,          INTENT(IN   )       :: n1           ! len of vec1
      CHARACTER(LEN=*), INTENT(IN   )       :: vec2(n2)     ! char. vector 2
      INTEGER,          INTENT(IN   )       :: n2           ! len of vec2
      INTEGER,          INTENT(  OUT)       :: matchidx(n1) ! index of vec2 in vec1
      INTEGER,          INTENT(  OUT)       :: nnmatch      ! # of matches
! 
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller - Initial version (update) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: I, J 

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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      FUNCTION IsInWord ( InString, SearchString ) RESULT ( Cnt ) 
!
! !USES:
!
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*), INTENT(IN   )    :: InString
      CHARACTER(LEN=*), INTENT(IN   )    :: SearchString
      LOGICAL                            :: Cnt
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------

      Cnt = INDEX( TRIM(InString), TRIM(SearchString) ) > 0

      END FUNCTION IsInWord
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      END MODULE HCO_CHARTOOLS_MOD
!EOM
