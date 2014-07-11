!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_tools_mod.F90
!
! !DESCRIPTION: Module HCO\_Tools\_Mod contains a collection of 
! miscellaneous helper routines used by the HEMCO modules.
! \\
! !INTERFACE: 
!
MODULE HCO_Tools_Mod
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
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES: 
!
  INTERFACE HCO_CharSplit
     MODULE PROCEDURE HCO_CharSplit_R8
     MODULE PROCEDURE HCO_CharSplit_R4
     MODULE PROCEDURE HCO_CharSplit_Int
  END INTERFACE HCO_CharSplit

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
    USE CharPak_Mod, ONLY : StrSplit
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN   ) :: CharStr  ! Character string
    CHARACTER(LEN=1), INTENT(IN   ) :: SEP      ! Separator
    CHARACTER(LEN=1), INTENT(IN   ) :: WC       ! Wildcard character
!
! !OUTPUT PARAMETERS:
!
    REAL(dp),         INTENT(  OUT) :: Reals(:) ! Output values
    INTEGER,          INTENT(  OUT) :: N        ! # of valid values
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC       ! Return code
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
    LOC = 'HCO_CharSplit_R8 (hco_tools_mod.F90)'

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
    USE CharPak_Mod, ONLY: StrSplit
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN   ) :: CharStr 
    CHARACTER(LEN=1), INTENT(IN   ) :: SEP 
    CHARACTER(LEN=1), INTENT(IN   ) :: WC 
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(  OUT) :: Reals(:)
    INTEGER,          INTENT(  OUT) :: N
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC
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
    LOC = 'HCO_CharSplit_R4 (hco_tools_mod.F90)'

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
! !DESCRIPTION: Subroutine HCO\_CharSplit\_INT splits the passed character
! string into N integers, using character SEP as separator. Wildcard
! values (WC) are set to -999. 
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CharSplit_Int( CharStr, SEP, WC, Ints, N, RC ) 
!
! !USES:
!
    USE CharPak_Mod, ONLY: StrSplit
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN   ) :: CharStr 
    CHARACTER(LEN=1), INTENT(IN   ) :: SEP 
    CHARACTER(LEN=1), INTENT(IN   ) :: WC 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT) :: Ints(:)
    INTEGER,          INTENT(  OUT) :: N
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC
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
    LOC = 'HCO_CharSplit_Int (hco_tools_mod.F90)'

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
  FUNCTION IsInWord ( InString, SearchString ) RESULT ( Cnt ) 
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )    :: InString
    CHARACTER(LEN=*), INTENT(IN   )    :: SearchString
!
! !RETURN VALUE
!
    LOGICAL                            :: Cnt
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------

    Cnt = INDEX( TRIM(InString), TRIM(SearchString) ) > 0

  END FUNCTION IsInWord
!EOC
END MODULE HCO_Tools_Mod
