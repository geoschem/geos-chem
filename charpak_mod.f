! $Id: charpak_mod.f,v 1.3 2004/12/02 21:48:33 bmy Exp $
      MODULE CHARPAK_MOD
!
!******************************************************************************
!  Module CHARPAK_MOD contains routines from the CHARPAK string and character
!  manipulation package used by GEOS-CHEM (bmy, 10/15/01, 7/20/04)
!
!  CHARPAK routines by Robert D. Stewart, 1992.  Subsequent modifications 
!  made for GEOS-CHEM by Bob Yantosca (1998, 2002, 2004).
!
!  Module Routines:
!  ============================================================================
!  (1 ) CNTMAT     : Counts # of chars in STR1 that match a char in STR2
!  (2 ) COPYTXT    : Writes chars from STR1 into STR2
!  (3 ) CSTRIP     : Strip blanks and null characters from a string 
!  (4 ) ISDIGIT    : Returns TRUE if a character is a numeric digit           
!  (5 ) STRREPL    : Replaces characters w/in a string with replacement text
!  (6 ) STRSPLIT   : Convenience wrapper for TXTEXT
!  (7 ) STRSQUEEZE : Squeezes text by removing white space from both ends
!  (8 ) TRANLC     : Translates text to LOWERCASE
!  (9 ) TRANUC     : Translates text to UPPERCASE
!  (10) TXT2INUM   : Converts a string of characters into an integer number
!  (11) TXTEXT     : Extracts a sequence of characters from a string
!
!  GEOS-CHEM modules referenced by charpak_mod.f
!  ============================================================================
!  none
!
!  NOTES:
!  (1 ) Moved "cntmat.f", "copytxt.f", "cstrip.f", "fillstr.f", "txt2inum.f",
!        "txtext.f", into this F90 module for easier bookkeeping 
!        (bmy, 10/15/01)
!  (2 ) Moved "tranuc.f" into this F90 module (bmy, 11/15/01)
!  (3 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (4 ) Wrote a new file "strrepl.f", which replaces a character pattern
!        within a string with replacement text.  Moved "tranlc.f" into
!        this module.  Replaced calls to function LENTRIM with F90 
!        intrinsic function LEN_TRIM.  Removed function FILLSTR and
!        replaced it w/ F90 intrinsic REPEAT. (bmy, 6/25/02)
!  (5 ) Added routine STRSPLIT as a wrapper for TXTEXT.  Also added
!        routines STRREPL and STRSQUEEZE. (bmy, 7/30/02)
!  (6 ) Added function ISDIGIT.  Also replace LEN_TRIM with LEN in routine
!        STRREPL, to allow us to replace tabs w/ spaces. (bmy, 7/20/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CntMat(str1,str2,imat)
C
C     Count the number of characters in str1 that match
C     a character in str2.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C          LENTRIM                CharPak
C
C     DATE:   JAN. 6, 1995
C     AUTHOR: R.D. STEWART
C     COMMENTS: Revised slightly (2-5-1996) so that trailing
C               blanks in str1 are ignored.  Revised again
C               on 3-6-1996.
C
      CHARACTER*(*) str1,str2
      INTEGER imat
      INTEGER L1,L2,i,j
      LOGICAL again

      L1 = MAX(1,LEN_TRIM(str1))
      L2 = LEN(str2)
      imat = 0
      DO i=1,L1
        again = .true.
        j = 1
        DO WHILE (again)
          IF (str2(j:j).EQ.str1(i:i)) THEN
            imat = imat+1
            again = .false.
          ELSEIF (j.LT.L2) THEN
            j=j+1
          ELSE
            again = .false.
          ENDIF
        ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CntMat

!------------------------------------------------------------------------------

      SUBROUTINE CopyTxt(col,str1,str2)
C
c     PURPOSE: Write all of the characters in str1 into variable
C              str2 beginning at column, col.  If the length of str1
C              + col is longer than the number of characters str2
C              can store, some characters will not be transfered to
C              str2.  Any characters already existing in str2 will
C              will be overwritten.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
C     DATE:   DEC. 24, 1993
C     AUTHOR: R.D. STEWART
C
      CHARACTER*(*) str2,str1
      INTEGER col,ilt1,i1,i,j,ic

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

      ! Return to calling program
      END SUBROUTINE CopyTxt

!------------------------------------------------------------------------------

      SUBROUTINE CSTRIP(text)
C
C     PURPOSE: Strip blanks and null characters for the variable TEXT.
C
C     COMMENTS: The original "text" is destroyed upon exit.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
C      AUTHOR: Robert D. Stewart
C        DATE: May 19, 1992
C
      CHARACTER*(*) TEXT
      INTEGER ilen,iasc,icnt,i
      CHARACTER*1 ch

      ilen = LEN(text)
      IF (ilen.GT.1) THEN
        icnt = 1
        DO i=1,ilen
          iasc = ICHAR(text(i:i))
          IF ((iasc.GT.32).AND.(iasc.LT.255)) THEN
C           Keep character
            ch = text(i:i)
            text(icnt:icnt) = ch
            icnt = icnt + 1
          ENDIF
        ENDDO
C       Fill remainder of text with blanks
        DO i=icnt,ilen
          text(i:i) = ' '
        ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CSTRIP

!------------------------------------------------------------------------------

      FUNCTION ISDIGIT( ch ) RESULT( LNUM )
C
C     Returned as true if ch is a numeric character (i.e., one of
C     the numbers from 0 to 9).
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
C     DATE:   NOV. 11, 1993
C     AUTHOR: R.D. STEWART
C
C     NOTE: Changed name from ISNUM to ISDIGIT (bmy, 7/15/04)
C  
      CHARACTER*1 ch
      INTEGER iasc
      LOGICAL lnum

      iasc = ICHAR(ch)
      lnum = .FALSE.
      IF ((iasc.GE.48).AND.(iasc.LE.57)) THEN
        lnum = .TRUE.
      ENDIF

      ! Return to calling program
      END FUNCTION ISDIGIT

!------------------------------------------------------------------------------

      SUBROUTINE StrRepl( STR, PATTERN, REPLTXT )

      !=================================================================
      ! Subroutine STRREPL replaces all instances of PATTERN within
      ! a string STR with replacement text REPLTXT. 
      ! (bmy, 6/25/02, 7/20/04)
      !
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) STR     : String to be searched
      ! (2 ) PATTERN : Pattern of characters to replace w/in STR
      ! (3 ) REPLTXT : Replacement text for PATTERN
      !
      ! Arguments as Output:
      ! ----------------------------------------------------------------
      ! (1 ) STR     : String with new replacement text 
      !
      ! NOTES
      ! (1 ) REPLTXT must have the same # of characters as PATTERN.
      ! (2 ) Replace LEN_TRIM with LEN (bmy, 7/20/04)
      !=================================================================

      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: STR
      CHARACTER(LEN=*), INTENT(IN)    :: PATTERN, REPLTXT
      
      ! Local variables
      INTEGER                         :: I1, I2

      !=================================================================
      ! STRREPL begins here!
      !=================================================================

      ! Error check: make sure PATTERN and REPLTXT have the same # of chars
      IF ( LEN( PATTERN ) /= LEN( REPLTXT ) ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 
     &    'STRREPL: PATTERN and REPLTXT must have same # of characters!'
         WRITE( 6, '(a)' ) 'STOP in STRREPL (charpak_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      ! Loop over all instances of PATTERN in STR
      DO 

         ! I1 is the starting location of PATTERN w/in STR  
         I1 = INDEX( STR, PATTERN )

         ! If pattern is not found, then return to calling program
         IF ( I1 < 1 ) RETURN

         ! I2 is the ending location of PATTERN w/in STR
         I2 = I1 + LEN_TRIM( PATTERN ) - 1
      
         ! Replace text
         STR(I1:I2) = REPLTXT

      ENDDO
         
      ! Return to calling program
      END SUBROUTINE StrRepl

!------------------------------------------------------------------------------

      SUBROUTINE StrSplit( STR, SEP, RESULT, N_SUBSTRS )
!
!******************************************************************************
!  Subroutine STRSPLIT returns substrings in a string, separated by a 
!  separator character (similar to IDL's StrSplit function).  This is mainly
!  a convenience wrapper for CHARPAK routine TxtExt. (bmy, 7/11/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) STR       (CHARACTER*(*)) : String to be searched (variable length)  
!  (2 ) SEP       (CHARACTER*1  ) : Separator character
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) RESULT    (CHARACTER*255) : Array containing substrings (255 elements)
!  (4 ) N_SUBSTRS (INTEGER      ) : Number of substrings returned (optional)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)            :: STR
      CHARACTER(LEN=1), INTENT(IN)            :: SEP
      CHARACTER(LEN=*), INTENT(OUT)           :: RESULT(255)
      INTEGER,          INTENT(OUT), OPTIONAL :: N_SUBSTRS

      ! Local variables
      INTEGER                                 :: I, IFLAG, COL
      CHARACTER (LEN=255)                     :: WORD

      !=================================================================
      ! STRSPLIT begins here!
      !=================================================================

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

      ! Return to calling program
      END SUBROUTINE StrSplit

!------------------------------------------------------------------------------

      SUBROUTINE StrSqueeze( STR )
!
!******************************************************************************
!  Subroutine STRSQUEEZE strips white space from both ends of a string.  
!  White space in the middle of the string (i.e. between characters) will
!  be preserved as-is.  Somewhat similar (though not exactly) to IDL's 
!  STRCOMPRESS function. (bmy, 7/11/02)
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) STR (CHAR*(*)) : String to be squeezed (will be overwritten in place!)
!
!  NOTES:
!******************************************************************************
!      
      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: STR

      !=================================================================
      ! STRSQUEEZE begins here!
      !=================================================================
      STR = ADJUSTR( TRIM( STR ) )
      STR = ADJUSTL( TRIM( STR ) )

      ! Return to calling program
      END SUBROUTINE StrSqueeze
      
!------------------------------------------------------------------------------

      SUBROUTINE TRANLC(text)
C
C     PURPOSE: Tranlate a character variable to all lowercase letters.
C              Non-alphabetic characters are not affected.
C
C    COMMENTS: The original "text" is destroyed.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
C      AUTHOR: Robert D. Stewart
C        DATE: May 19, 1992
C
      CHARACTER*(*) text
      INTEGER iasc,i,ilen

      ilen = LEN(text)
      DO I=1,ilen
        iasc = ICHAR(text(i:i))
        IF ((iasc.GT.64).AND.(iasc.LT.91)) THEN
          text(i:i) = CHAR(iasc+32)
        ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE TRANLC

!------------------------------------------------------------------------------

      SUBROUTINE TRANUC(text)
C
C     PURPOSE: Tranlate a character variable to all upper case letters.
C              Non-alphabetic characters are not affected.
C
C    COMMENTS: The original "text" is destroyed.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
C      AUTHOR: Robert D. Stewart
C        DATE: May 19, 1992
C
      CHARACTER*(*) text
      INTEGER iasc,i,ilen

      ilen = LEN(text)
      DO i=1,ilen
        iasc = ICHAR(text(i:i))
        IF ((iasc.GT.96).AND.(iasc.LT.123)) THEN
          text(i:i) = CHAR(iasc-32)
        ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE TRANUC

!------------------------------------------------------------------------------

      SUBROUTINE Txt2Inum(fmat,txt,Inum,iflg)
C
C     <Txt2Inum> attempts to convert the string of characters
C     in txt into a integer number.  fmat is the
C     VALID format specifier to use in the internal read
C     statement.  iflg is returned as a status flag indicating
C     the success or failure of the operation.  iflg <=0 if the
C     operation was successful, and > 0 if it failed.
C
C     COMMENTS: Generally, the Fxx.0 format should be used to convert
C               string of characters to a number.
C
C      AUTHOR: Robert D. Stewart
C        DATE: DEC 24, 1992
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        N/A
C
      CHARACTER*(*) txt,fmat
      INTEGER inum
      INTEGER iflg

      READ(txt,fmt=fmat,iostat=iflg) inum

      ! Return to calling program
      END SUBROUTINE Txt2Inum

!------------------------------------------------------------------------------

      SUBROUTINE TxtExt(ch,text,col,word,iflg)
C
C     PURPOSE: TxtExt extracts a sequence of characters from
C              text and transfers them to word.  The extraction
C              procedure uses a set of character "delimiters"
C              to denote the desired sequence of characters.
C              For example if ch=' ', the first character sequence
C              bracketed by blank spaces will be returned in word.
C              The extraction procedure begins in column, col,
C              of TEXT.  If text(col:col) = ch (any character in
C              the character string), the text is returned beginning
C              with col+1 in text (i.e., the first match with ch
C              is ignored).
C
C              After completing the extraction, col is incremented to
C              the location of the first character following the
C              end of the extracted text.
C
C              A status flag is also returned with the following
C              meaning(s)
C
C                 IF iflg = -1, found a text block, but no more characters
C                               are available in TEXT
C                    iflg = 0, task completed sucessfully (normal term)
C                    iflg = 1, ran out of text before finding a block of
C                              text.
C
C       COMMENTS: TxtExt is short for Text Extraction.  This routine
C                 provides a set of powerful line-by-line
C                 text search and extraction capabilities in
C                 standard FORTRAN.
C
C     CODE DEPENDENCIES:
C      Routine Name                  File
C        CntMat                    CHARPAK.FOR
C        TxtExt                    CHARPAK.FOR
C        FillStr                   CHARPAK.FOR
C        CopyTxt                   CHARPAK.FOR
C
C        other routines are indirectly called.
C      AUTHOR: Robert D. Stewart
C        DATE: Jan. 1st, 1995
C
C      REVISIONS: FEB 22, 1996.  Slight bug fix (introduced by a
C        (recent = FLIB 1.04) change in the CntMat routine)
C        so that TxtExt correctlyhandles groups of characters
C        delimited by blanks).
C
C      MODIFICATIONS by Bob Yantosca (6/25/02)
C        (1) Replace call to FILLSTR with F90 intrinsic REPEAT
C
      CHARACTER*(*) ch,text,word
      INTEGER col,iflg
      INTEGER Tmax,T1,T2,imat
      LOGICAL again,prev

C     Length of text
      Tmax = LEN(text)

C     Fill Word with blanks
      WORD = REPEAT( ' ', LEN( WORD ) )
      
      IF (col.GT.Tmax) THEN
C       Text does not contain any characters past Tmax.
C       Reset col to one and return flag = {error condition}
        iflg = 1
        col = 1
      ELSEIF (col.EQ.Tmax) THEN
C       End of TEXT reached
        CALL CntMat(ch,text(Tmax:Tmax),imat)
        IF (imat.EQ.0) THEN
C         Copy character into Word and set col=1
          CALL CopyTxt(1,Text(Tmax:Tmax),Word)
          col = 1
          iflg = -1
        ELSE
C         Same error condition as if col.GT.Tmax
          iflg = 1
        ENDIF
      ELSE
C       Make sure column is not less than 1
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
C         Check for a match with a character in ch
          CALL CntMat(ch,text(T2:T2),imat)
          IF (imat.GT.0) THEN
C           Current character in TEXT matches one (or more) of the
C           characters in ch.
            IF (prev) THEN
              IF (T2.LT.Tmax) THEN
C               Keep searching for a block of text
                T2=T2+1
                T1=T2
              ELSE
C               Did not find any text blocks before running
C               out of characters in TEXT.
                again=.false.
                iflg=1
              ENDIF
            ELSE
C             Previous character did not match ch, so terminate.
C             NOTE: This is "NORMAL" termination of the loop
              again=.false.
              T2=T2-1
              iflg = 0
            ENDIF
          ELSEIF (T2.LT.Tmax) THEN
C           Add a letter to the current block of text
            prev = .false.
            T2=T2+1
          ELSE
C           Reached the end of the characters in TEXT before reaching
C           another delimiting character.  A text block was identified
C           however.
            again=.false.
            iflg=-1
          ENDIF
        ENDDO

        IF (iflg.EQ.0) THEN
C         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = T2+1
        ELSE
C         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = 1
        ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE TxtExt

!------------------------------------------------------------------------------

      END MODULE CHARPAK_MOD
