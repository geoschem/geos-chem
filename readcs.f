! $Id: readcs.f,v 1.2 2004/03/24 20:52:31 bmy Exp $
      SUBROUTINE READCS(LDEBUG)
      
      ! References to F90 modules (bmy, 10/15/02)
      USE ERROR_MOD, ONLY : ERROR_STOP

      IMPLICIT NONE

C**********************************************************************
C                                                                     *
C  HARVARD TROPOSPHERIC CHEMISTRY MODULE FOR 3-D APPLICATIONS         *
C  by Larry Horowitz, Jinyou Liang, Gerry Gardner, Prof. Daniel Jacob *
C  of HARVARD UNIVERSITY    (Release V2.0)                            *
C**********************************************************************
C
C**************************************************************
C                                                              
C This block reads in 8-column data that jcx just coded in. The
C data is saved in XSECT[IWL, ISPEC, ITEMP, IBRCH].
C
C**************************************************************

#     include "CMN_SIZE"
#     include "comsol.h"
#     include "CMN_SETUP" ! (bmy, 4/3/02)
      LOGICAL LDEBUG

      CHARACTER*8 CHARSPEC, TITLE, BLANK
      INTEGER NTEMP(MXSPE,MXBRCH), NBRCH (MXSPE)
      REAL*8 XXSECT(MXWL)
      
      INTEGER IX,IWL,J,MCS,ISPEC,LCS,NWL,NSPEC,MTEMP,I,IBRCH,ITEMP
      
      REAL*8 XA
      
      DATA BLANK/'        '/
      DATA LCS/40/
      DATA NWL/MXWL/
      DATA MCS/41/
      CHARACTER(LEN=255) :: FILENAME

      ! Prior to 4/3/02:
      !OPEN(LCS, FILE='8col.dat',STATUS='old')

      ! Define file name
      FILENAME = TRIM( DATA_DIR ) // 'slowj_200203/8col.dat' 

      ! Now read file directly from DATA_DIR/slowj_200203 (bmy, 4/3/02)
      OPEN( LCS, FILE=TRIM( FILENAME ), STATUS='old' )

      ISPEC=0
      NSPEC=0
      IBRCH=0
      ITEMP=0
C
C Zero out ARRAYS
C
      DO I =1, MXSPE
         NBRCH (I) = 0
         DO J =1, MXBRCH
            NTEMP(MXSPE,MXBRCH) = 0
         ENDDO
      ENDDO
C
C Begin a reading loop
C
 100  CONTINUE
C
C Check if the # of wavelength is as specified
C
      IF(LDEBUG) REWIND MCS

      READ(LCS, 1001, END=250) CHARSPEC, IBRCH, MTEMP
      IF(LDEBUG) WRITE (MCS, 1001) CHARSPEC, IBRCH, MTEMP
      DO 120 I=1,8
         IF('+' .EQ. CHARSPEC(I:I)) THEN
            DO 110 J=I, 8
               CHARSPEC(J:J)=' '
  110       CONTINUE
            GOTO 130
         ENDIF
 120  CONTINUE
 130  CONTINUE

      IF(LDEBUG) WRITE(MCS, "('LINE 20')")
      DO 160 I=1, NSPEC
         IF (CNAME(I).EQ.CHARSPEC) THEN
C It's an old species
            ISPEC= I
            IF(IBRCH.EQ.NBRCH(I)) THEN
C For new temperature
               IX=1
               DO 140 J=1, NTEMP(I,IBRCH)
                  IF (ABS(MTEMP-TARRAY(I,IBRCH,J)).LT.1.D-20) THEN
                     IX=0
                     ITEMP=J
                     GOTO 150
                  ENDIF
 140           CONTINUE
               IF(LDEBUG) WRITE(MCS, "('LINE 35')")
               NTEMP(ISPEC, IBRCH)= NTEMP(I,IBRCH)+ IX
               ITEMP= NTEMP(I,IBRCH)
 150           CONTINUE
               TARRAY(I,IBRCH,ITEMP)= REAL(MTEMP)
            ELSE
C For new branch
C The following lines suppose the branchs arrange in increasing order!
               NBRCH(I)= IBRCH
               IF(LDEBUG) WRITE(MCS, '("LINE 37")')
C New temperature and new branch
               ITEMP=1
               NTEMP(I, IBRCH)=ITEMP
               TARRAY(I,IBRCH,ITEMP)= REAL(MTEMP)
               IF(LDEBUG) WRITE(MCS, "('LINE 38')")
            ENDIF
*** END FOR old species manipulation
            GOTO 170
         ENDIF
 160  CONTINUE
      NSPEC=NSPEC + 1
      ISPEC= NSPEC
      ITEMP= 1
      CNAME(ISPEC)= CHARSPEC
C Ibrch is supposed 1 here, otherwise need revise
      NBRCH(ISPEC)= IBRCH     
      NTEMP(ISPEC,IBRCH)= 1
      TARRAY(ISPEC,IBRCH,ITEMP)= REAL(MTEMP)
 170  CONTINUE
C Read in wavelength, cross section, quantum yield.
      READ(LCS, 1003, ERR=2000) XXSECT
      IF(LDEBUG) WRITE(MCS, 1003) XXSECT 
C Next reaction begins if wavelength is less than 1.e-20.
      DO 180 IWL=1, NWL
         XSECT(IWL, ISPEC, ITEMP, IBRCH)=XXSECT(IWL)
 180  CONTINUE      
      READ (LCS, 1002, END=250) TITLE
C Continue reading.
      IF(TITLE.NE.'end     ' .AND. TITLE.NE.'END     ') GOTO 100
**** end reading
**** Rearrange TARRAY/XSECT to ensure temp. is in increasing order.
      DO 230 ISPEC=1, NSPEC
         DO 220 IBRCH=1, NBRCH(ISPEC)
            DO 210 ITEMP=1, NTEMP(ISPEC, IBRCH)-1
               DO 200 I= ITEMP+1, NTEMP(ISPEC, IBRCH)
**** POSSIBLE ERROR
                  IF(TARRAY(ISPEC,IBRCH,I).LT.1.D-20) GOTO 240
                  IF(TARRAY(ISPEC,IBRCH,I).LT.
     x                 TARRAY(ISPEC,IBRCH,ITEMP)) THEN
                     
                     XA=TARRAY(ISPEC,IBRCH,ITEMP)
                     TARRAY(ISPEC,IBRCH,ITEMP)=TARRAY(ISPEC,IBRCH,I)
                     TARRAY(ISPEC,IBRCH,I)= XA 
                     DO 190 J=1, NWL
                        XA= XSECT(J,ISPEC,ITEMP,IBRCH)
                        XSECT(J,ISPEC,ITEMP,IBRCH)=
     x                       XSECT(J,ISPEC,I,IBRCH)
                        XSECT(J,ISPEC,I,IBRCH)= XA
 190                 CONTINUE
                  ENDIF
 200           CONTINUE
 210        CONTINUE
 220     CONTINUE
 230  CONTINUE
      GOTO 250
 240  CONTINUE
***** ERROR in tarray
C     WRITE(MCS,1004 ) (ispec,ibrch,(TARRAY(ISPEC,IBRCH,J), J=1,I))
      CALL ERROR_STOP( 'STOP_300', 'readcs.f' )
 250  CONTINUE
      CLOSE (LCS)
***  for debug only.
      IF (LDEBUG) THEN
         WRITE(*,*) '8col.dat fine'
         WRITE(*,1010) NWL, ISPEC, IBRCH, ITEMP
      ENDIF
      RETURN
 1001 FORMAT(A8, 43X, I1, 5X, I3, 5x, a8)
 1002 FORMAT(A8)
 1003 FORMAT(8E10.3)
 1004 FORMAT(' ntemp is incorrect: is/ib/tarray= ',2I4,5E10.1)
 1010 FORMAT(1X,'LAST input of NWL,IS, IB, IT are: '/1x,4I8)
 2000 CONTINUE
C Err in format or code
      WRITE (MCS, 1001) CHARSPEC, IBRCH, MTEMP
      WRITE(*,'("XSECT FORMAT ERR")')
      CALL ERROR_STOP( 'STOP 1001', 'readcs.f' )

      ! Return to calling program
      END SUBROUTINE READCS
