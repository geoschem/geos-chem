!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: charpak_mod
!
! !DESCRIPTION: Module CHARPAK\_MOD contains routines from the CHARPAK 
!  string and character manipulation package used by GEOS-Chem.
!\\   
!\\
! !INTERFACE: 
!
      MODULE CHARPAK_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CNTMAT     
      PUBLIC  :: COPYTXT    
      PUBLIC  :: CSTRIP     
      PUBLIC  :: ISDIGIT    
      PUBLIC  :: STRREPL    
      PUBLIC  :: STRSPLIT   
      PUBLIC  :: STRSQUEEZE 
      PUBLIC  :: TRANLC     
      PUBLIC  :: TRANUC     
      PUBLIC  :: TXT2INUM   
      PUBLIC  :: TXTEXT     
!
! !REMARKS:
!  CHARPAK routines by Robert D. Stewart, 1992.  Subsequent modifications 
!  made for GEOS-CHEM by Bob Yantosca (1998, 2002, 2004).
!
! !REVISION HISTORY:
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
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------

      SUBROUTINE CntMat(str1,str2,imat)
!
!     Count the number of characters in str1 that match
!     a character in str2.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!          LENTRIM                CharPak
!
!     DATE:   JAN. 6, 1995
!     AUTHOR: R.D. STEWART
!     COMMENTS: Revised slightly (2-5-1996) so that trailing
!               blanks in str1 are ignored.  Revised again
!               on 3-6-1996.
!
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
!
c     PURPOSE: Write all of the characters in str1 into variable
!              str2 beginning at column, col.  If the length of str1
!              + col is longer than the number of characters str2
!              can store, some characters will not be transfered to
!              str2.  Any characters already existing in str2 will
!              will be overwritten.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   DEC. 24, 1993
!     AUTHOR: R.D. STEWART
!
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
!
!     PURPOSE: Strip blanks and null characters for the variable TEXT.
!
!     COMMENTS: The original "text" is destroyed upon exit.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!
      CHARACTER*(*) TEXT
      INTEGER ilen,iasc,icnt,i
      CHARACTER*1 ch

      ilen = LEN(text)
      IF (ilen.GT.1) THEN
        icnt = 1
        DO i=1,ilen
          iasc = ICHAR(text(i:i))
          IF ((iasc.GT.32).AND.(iasc.LT.255)) THEN
!           Keep character
            ch = text(i:i)
            text(icnt:icnt) = ch
            icnt = icnt + 1
          ENDIF
        ENDDO
!       Fill remainder of text with blanks
        DO i=icnt,ilen
          text(i:i) = ' '
        ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CSTRIP

!------------------------------------------------------------------------------

      FUNCTION ISDIGIT( ch ) RESULT( LNUM )
!
!     Returned as true if ch is a numeric character (i.e., one of
!     the numbers from 0 to 9).
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!     DATE:   NOV. 11, 1993
!     AUTHOR: R.D. STEWART
!
!     NOTE: Changed name from ISNUM to ISDIGIT (bmy, 7/15/04)
!  
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
!
!     PURPOSE: Tranlate a character variable to all lowercase letters.
!              Non-alphabetic characters are not affected.
!
!    COMMENTS: The original "text" is destroyed.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!
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
!
!     PURPOSE: Tranlate a character variable to all upper case letters.
!              Non-alphabetic characters are not affected.
!
!    COMMENTS: The original "text" is destroyed.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
!      AUTHOR: Robert D. Stewart
!        DATE: May 19, 1992
!
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
!
!     <Txt2Inum> attempts to convert the string of characters
!     in txt into a integer number.  fmat is the
!     VALID format specifier to use in the internal read
!     statement.  iflg is returned as a status flag indicating
!     the success or failure of the operation.  iflg <=0 if the
!     operation was successful, and > 0 if it failed.
!
!     COMMENTS: Generally, the Fxx.0 format should be used to convert
!               string of characters to a number.
!
!      AUTHOR: Robert D. Stewart
!        DATE: DEC 24, 1992
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        N/A
!
      CHARACTER*(*) txt,fmat
      INTEGER inum
      INTEGER iflg

      READ(txt,fmt=fmat,iostat=iflg) inum

      ! Return to calling program
      END SUBROUTINE Txt2Inum

!------------------------------------------------------------------------------

      SUBROUTINE TxtExt(ch,text,col,word,iflg)
!
!     PURPOSE: TxtExt extracts a sequence of characters from
!              text and transfers them to word.  The extraction
!              procedure uses a set of character "delimiters"
!              to denote the desired sequence of characters.
!              For example if ch=' ', the first character sequence
!              bracketed by blank spaces will be returned in word.
!              The extraction procedure begins in column, col,
!              of TEXT.  If text(col:col) = ch (any character in
!              the character string), the text is returned beginning
!              with col+1 in text (i.e., the first match with ch
!              is ignored).
!
!              After completing the extraction, col is incremented to
!              the location of the first character following the
!              end of the extracted text.
!
!              A status flag is also returned with the following
!              meaning(s)
!
!                 IF iflg = -1, found a text block, but no more characters
!                               are available in TEXT
!                    iflg = 0, task completed sucessfully (normal term)
!                    iflg = 1, ran out of text before finding a block of
!                              text.
!
!       COMMENTS: TxtExt is short for Text Extraction.  This routine
!                 provides a set of powerful line-by-line
!                 text search and extraction capabilities in
!                 standard FORTRAN.
!
!     CODE DEPENDENCIES:
!      Routine Name                  File
!        CntMat                    CHARPAK.FOR
!        TxtExt                    CHARPAK.FOR
!        FillStr                   CHARPAK.FOR
!        CopyTxt                   CHARPAK.FOR
!
!        other routines are indirectly called.
!      AUTHOR: Robert D. Stewart
!        DATE: Jan. 1st, 1995
!
!      REVISIONS: FEB 22, 1996.  Slight bug fix (introduced by a
!        (recent = FLIB 1.04) change in the CntMat routine)
!        so that TxtExt correctlyhandles groups of characters
!        delimited by blanks).
!
!      MODIFICATIONS by Bob Yantosca (6/25/02)
!        (1) Replace call to FILLSTR with F90 intrinsic REPEAT
!
      CHARACTER*(*) ch,text,word
      INTEGER col,iflg
      INTEGER Tmax,T1,T2,imat
      LOGICAL again,prev

!     Length of text
      Tmax = LEN(text)

!     Fill Word with blanks
      WORD = REPEAT( ' ', LEN( WORD ) )
      
      IF (col.GT.Tmax) THEN
!       Text does not contain any characters past Tmax.
!       Reset col to one and return flag = {error condition}
        iflg = 1
        col = 1
      ELSEIF (col.EQ.Tmax) THEN
!       End of TEXT reached
        CALL CntMat(ch,text(Tmax:Tmax),imat)
        IF (imat.EQ.0) THEN
!         Copy character into Word and set col=1
          CALL CopyTxt(1,Text(Tmax:Tmax),Word)
          col = 1
          iflg = -1
        ELSE
!         Same error condition as if col.GT.Tmax
          iflg = 1
        ENDIF
      ELSE
!       Make sure column is not less than 1
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
!         Check for a match with a character in ch
          CALL CntMat(ch,text(T2:T2),imat)
          IF (imat.GT.0) THEN
!           Current character in TEXT matches one (or more) of the
!           characters in ch.
            IF (prev) THEN
              IF (T2.LT.Tmax) THEN
!               Keep searching for a block of text
                T2=T2+1
                T1=T2
              ELSE
!               Did not find any text blocks before running
!               out of characters in TEXT.
                again=.false.
                iflg=1
              ENDIF
            ELSE
!             Previous character did not match ch, so terminate.
!             NOTE: This is "NORMAL" termination of the loop
              again=.false.
              T2=T2-1
              iflg = 0
            ENDIF
          ELSEIF (T2.LT.Tmax) THEN
!           Add a letter to the current block of text
            prev = .false.
            T2=T2+1
          ELSE
!           Reached the end of the characters in TEXT before reaching
!           another delimiting character.  A text block was identified
!           however.
            again=.false.
            iflg=-1
          ENDIF
        ENDDO

        IF (iflg.EQ.0) THEN
!         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = T2+1
        ELSE
!         Copy characters into WORD and set col for return
          CALL CopyTxt(1,Text(T1:T2),Word)
          col = 1
        ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE TxtExt
!EOC
      END MODULE CHARPAK_MOD
