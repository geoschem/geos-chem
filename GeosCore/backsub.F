! $Id: backsub.f,v 1.1 2009/09/16 14:06:41 bmy Exp $
      SUBROUTINE BACKSUB
!
!******************************************************************************
!  Subroutine BACKSUB does the back-substitution on the decomposed matrix.
!  (M. Jacobson 1997; bdf, bmy, 4/1/03, 7/9/03)
!
!  NOTES:
!  (1 ) Comment out counter variable NUM_BACKSUB, you can get the same info
!        w/ a profiling run. (bmy, 7/9/03)
!******************************************************************************
!
      USE CMN_SIZE_MOD  ! Size parameters
      USE COMODE_LOOP_MOD  ! SMVGEAR II arrays

      IMPLICIT NONE
#     include "define.h"
#     include "define.h"

      ! Local variables
      INTEGER IJ,I,KZT,KL5,KH5,KL4,KH4,KL3,KH3,KL2,KH2,KL1,KH1,KC
      INTEGER J0,IJ0,IJ1,IJ2,IJ3,IJ4,J1,J2,J3,J4,K,MZT,ML5,MH5,ML4,MH4
      INTEGER ML3,MH3,ML2,MH2,ML1,MH1,MC
C
C *********************************************************************
C ************        WRITTEN BY MARK JACOBSON (1993)      ************
C ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
C ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
C ***                         (650) 723-6836                        *** 
C *********************************************************************
C
C     BBBBBBB     A     CCCCCCC  K     K  SSSSSSS  U      U  BBBBBBB
C     B     B    A A    C        K   K    S        U      U  B     B
C     BBBBBBB   A   A   C        K K      SSSSSSS  U      U  BBBBBBB
C     B     B  AAAAAAA  C        K   K          S  U      U  B     B
C     BBBBBBB A       A CCCCCCC  K     K  SSSSSSS  UUUUUUUU  BBBBBBB
C
C *********************************************************************
C ******* PERFORM BACK-SUBSTITUTIONS ON THE DECOMPOSED MATRIX   ******* 
C *********************************************************************

C *********************************************************************
C * THIS SUBROUTINE SOLVES THE LINEAR SET OF EQUATIONS Ax = B FOR x,  *
C * THE CORRECTION VECTOR, WHERE "A" IS THE L-U DECOMPOSTION OF THE   *
C * ORIGINAL MATRIX,                                                  * 
C *                                                                   *
C *                       P = I - H x Bo x J,                         * 
C *                                                                   *
C * I = IDENTITY MATRIX, H = TIME-STEP, Bo = A COEFFICIENT THAT       * 
C * DEPENDS ON THE ORDER OF THE INTEGRATION METHOD, AND J IS THE      * 
C * MATRIX OF PARTIAL DERIVATIVES. B IS SENT FROM SMVGEAR AS A        * 
C * CORRECTED VALUE OF THE FIRST DERIVATIVES OF THE ORDINARY DIFFER-  * 
C * ENTIAL EQUATIONS. SUBROUTINE DECOMP.F SOLVED FOR "A", THE         * 
C * DECOMPOSED MATRIX. SEE PRESS ET AL. (1992) NUMERICAL RECIPES.     *
C * CAMBRIDGE UNIVERSITY PRESS, FOR A BETTER DESCRIPTION OF THE BACK- *
C * SUBSTITUTION PROCESS.                                             *
C *                                                                   *
C * THIS BACK-SUBSTITUTION PROCESS USES SPARSE-MATRIX TECHNIQUES,     *
C * VECTORIZES AROUND THE GRID-CELL DIMENSION, AND USES NO PARTIAL    *
C * PIVOTING. TESTS BY SHERMAN & HINDMARSH (1980) LAWRENCE LIVERMORE  * 
C * REP. UCRL-84102 AND BY US HAVE CONFIRMED THAT THE REMOVAL OF      *
C * PARTIAL PIVOTING HAS LITTLE EFFECT ON RESULTS.                    *
C *                                                                   *
C * HOW TO CALL SUBROUTINE:                                           *
C * ----------------------                                            *
C *  CALL BACKSUB.F FROM SMVGEAR.F WITH                               * 
C *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
C *     NCSP = NCS        FOR DAYTIME   GAS CHEM                      *  
C *     NCSP = NCS   +ICS FOR NIGHTTIME GAS CHEM                      *  
C *********************************************************************
C
C *********************************************************************
C *                        BACKSUB LOOP # 1                           *
C * FIRST, ADJUST RIGHT SIDE OF Ax = B USING LOWER TRIANGULAR MATRIX  * 
C *********************************************************************
C SUM 1,2,3,4, OR 5 TERMS AT A TIME TO IMPROVE VECTORIZATION.
C
C KTLOOP   = NUMBER OF GRID-CELLS IN A GRID-BLOCK
C ISCHAN   = ORDER OF MATRIX
C CC2      = ARRAY HOLDING VALUES OF DECOMPOSED MATRIX. 
C GLOSS    = ARRAY INITIALLY HOLDING RIGHT SIDE OF EQUATION. THESE
C            VALUES ARE CONVERTED TO THE SOLUTION DURING BACK-SUBSTITUTION.
C KZEROA,..= ARRAYS IDENTIFYING TERMS IN GLOSS ARRAY 
C
      IJ            = 1 
      DO 310 KZT    = KZTLO(NCSP), KZTHI(NCSP)
         I            = IKZTOT(KZT)
         KL5          = KBL5(  KZT)
         KH5          = KBH5(  KZT)
         KL4          = KBL4(  KZT)
         KH4          = KBH4(  KZT)
         KL3          = KBL3(  KZT)
         KH3          = KBH3(  KZT)
         KL2          = KBL2(  KZT)
         KH2          = KBH2(  KZT)
         KL1          = KBL1(  KZT)
         KH1          = KBH1(  KZT)
C
C ***********************  SUM 5 TERMS AT A TIME  ********************* 
C
       DO 105 KC    = KL5, KH5
          IJ0         = IJ
          IJ1         = IJ + 1
          IJ2         = IJ + 2 
          IJ3         = IJ + 3 
          IJ4         = IJ + 4 
          IJ          = IJ + 5      
          J0          = KZEROA(KC)
          J1          = KZEROB(KC) 
          J2          = KZEROC(KC)
          J3          = KZEROD(KC)
          J4          = KZEROE(KC)
          DO 100 K    = 1, KTLOOP
             GLOSS(K,I) = GLOSS(K,I)
     1                  - CC2(K,IJ0) * GLOSS(K,J0)
     2                  - CC2(K,IJ1) * GLOSS(K,J1)
     3                  - CC2(K,IJ2) * GLOSS(K,J2)
     4                  - CC2(K,IJ3) * GLOSS(K,J3)
     5                  - CC2(K,IJ4) * GLOSS(K,J4)
 100      CONTINUE
 105   CONTINUE
C
C ***********************  SUM 4 TERMS AT A TIME  ********************* 
C     
       DO 155 KC    = KL4, KH4 
          IJ0         = IJ
          IJ1         = IJ + 1
          IJ2         = IJ + 2 
          IJ3         = IJ + 3 
          IJ          = IJ + 4      
          J0          = KZEROA(KC)
          J1          = KZEROB(KC) 
          J2          = KZEROC(KC)
          J3          = KZEROD(KC)
          DO 150 K    = 1, KTLOOP
             GLOSS(K,I) = GLOSS(K,I)
     1                  - CC2(K,IJ0) * GLOSS(K,J0)
     2                  - CC2(K,IJ1) * GLOSS(K,J1)
     3                  - CC2(K,IJ2) * GLOSS(K,J2)
     4                  - CC2(K,IJ3) * GLOSS(K,J3)
 150      CONTINUE
 155   CONTINUE
C
C ***********************  SUM 3 TERMS AT A TIME  ********************* 
C
       DO 205 KC    = KL3, KH3 
          IJ0         = IJ
          IJ1         = IJ + 1
          IJ2         = IJ + 2 
          IJ          = IJ + 3       
          J0          = KZEROA(KC)
          J1          = KZEROB(KC) 
          J2          = KZEROC(KC)
          DO 200 K    = 1, KTLOOP
             GLOSS(K,I) = GLOSS(K,I)
     1                  - CC2(K,IJ0) * GLOSS(K,J0)
     2                  - CC2(K,IJ1) * GLOSS(K,J1)
     3                  - CC2(K,IJ2) * GLOSS(K,J2)
 200      CONTINUE
 205   CONTINUE
C
C ***********************  SUM 2 TERMS AT A TIME  ********************* 
C
       DO 255 KC    = KL2, KH2 
          IJ0         = IJ
          IJ1         = IJ + 1
          IJ          = IJ + 2        
          J0          = KZEROA(KC)
          J1          = KZEROB(KC) 
          DO 250 K    = 1, KTLOOP
             GLOSS(K,I) = GLOSS(K,I)
     1                  - CC2(K,IJ0) * GLOSS(K,J0)
     2                  - CC2(K,IJ1) * GLOSS(K,J1)
 250      CONTINUE
 255   CONTINUE
C
C ***********************  SUM 1 TERM AT A TIME  ********************** 
C
       DO 305 KC    = KL1, KH1 
          IJ0         = IJ
          IJ          = IJ + 1         
          J0          = KZEROA(KC)
          DO 300 K    = 1, KTLOOP
             GLOSS(K,I) = GLOSS(K,I)
     1                  - CC2(K,IJ0) * GLOSS(K,J0)
 300      CONTINUE
 305   CONTINUE
 310  CONTINUE
C
C *********************************************************************
C *                         BACKSUB LOOP # 2                          * 
C * BACKSUBSTITE WITH UPPER TRIANGULAR MATRIX TO FIND SOLUTION        *
C *********************************************************************
C AGAIN, SUM UP SEVERAL TERMS AT A TIME TO IMPROVE VECTORIZATION.
C VDIAG  = DIAGONAL TERM FROM L-U DECOMPOSTION.
C GLOSS  = SOLUTION ON OUTPUT
C
      DO 710 I       = ISCHAN, 1, -1
         MZT           = IMZTOT(I,NCSP)  
         IF (MZT.GT.0) THEN
            ML5          = MBL5(  MZT)
            MH5          = MBH5(  MZT)
            ML4          = MBL4(  MZT)
            MH4          = MBH4(  MZT)
            ML3          = MBL3(  MZT)
            MH3          = MBH3(  MZT)
            ML2          = MBL2(  MZT)
            MH2          = MBH2(  MZT)
            ML1          = MBL1(  MZT)
            MH1          = MBH1(  MZT)
C
C ***********************  SUM 5 TERMS AT A TIME  ********************* 
C
            DO 405 MC    = ML5, MH5
               IJ0         = IJ
               IJ1         = IJ + 1
               IJ2         = IJ + 2 
               IJ3         = IJ + 3 
               IJ4         = IJ + 4 
               IJ          = IJ + 5      
               J0          = MZEROA(MC)
               J1          = MZEROB(MC) 
               J2          = MZEROC(MC)
               J3          = MZEROD(MC)
               J4          = MZEROE(MC)
               DO 400 K    = 1, KTLOOP
                  GLOSS(K,I) = GLOSS(K,I)
     1                       - CC2(K,IJ0) * GLOSS(K,J0)
     2                       - CC2(K,IJ1) * GLOSS(K,J1)
     3                       - CC2(K,IJ2) * GLOSS(K,J2)
     4                       - CC2(K,IJ3) * GLOSS(K,J3)
     5                       - CC2(K,IJ4) * GLOSS(K,J4)
 400           CONTINUE
 405        CONTINUE
C
C ***********************  SUM 4 TERMS AT A TIME  ********************* 
C
            DO 455 MC    = ML4, MH4 
               IJ0         = IJ
               IJ1         = IJ + 1
               IJ2         = IJ + 2 
               IJ3         = IJ + 3 
               IJ          = IJ + 4      
               J0          = MZEROA(MC)
               J1          = MZEROB(MC) 
               J2          = MZEROC(MC)
               J3          = MZEROD(MC)
               DO 450 K    = 1, KTLOOP
                  GLOSS(K,I) = GLOSS(K,I)
     1                       - CC2(K,IJ0) * GLOSS(K,J0)
     2                       - CC2(K,IJ1) * GLOSS(K,J1)
     3                       - CC2(K,IJ2) * GLOSS(K,J2)
     4                       - CC2(K,IJ3) * GLOSS(K,J3)
 450           CONTINUE
 455        CONTINUE
C
C ***********************  SUM 3 TERMS AT A TIME  ********************* 
C
            DO 505 MC    = ML3, MH3 
               IJ0         = IJ
               IJ1         = IJ + 1
               IJ2         = IJ + 2 
               IJ          = IJ + 3      
               J0          = MZEROA(MC)
               J1          = MZEROB(MC) 
               J2          = MZEROC(MC)
               DO 500 K    = 1, KTLOOP
                  GLOSS(K,I) = GLOSS(K,I)
     1                       - CC2(K,IJ0) * GLOSS(K,J0)
     2                       - CC2(K,IJ1) * GLOSS(K,J1)
     3                       - CC2(K,IJ2) * GLOSS(K,J2)
 500           CONTINUE
 505        CONTINUE
C     
C ***********************  SUM 2 TERMS AT A TIME  ********************* 
C
            DO 555 MC    = ML2, MH2 
               IJ0         = IJ
               IJ1         = IJ + 1
               IJ          = IJ + 2      
               J0          = MZEROA(MC)
               J1          = MZEROB(MC) 
               DO 550 K    = 1, KTLOOP
                  GLOSS(K,I) = GLOSS(K,I)
     1                       - CC2(K,IJ0) * GLOSS(K,J0)
     2                       - CC2(K,IJ1) * GLOSS(K,J1)
 550           CONTINUE
 555        CONTINUE
C     
C ***********************  SUM 1 TERM AT A TIME  ********************** 
C
            DO 605 MC    = ML1, MH1 
               IJ0         = IJ
               IJ          = IJ + 1       
               J0          = MZEROA(MC)
               DO 600 K    = 1, KTLOOP
                  GLOSS(K,I) = GLOSS(K,I)
     1                       - CC2(K,IJ0) * GLOSS(K,J0)
 600           CONTINUE
 605        CONTINUE
         ENDIF
C      ENDIF MZT.GT.0
C
C ***************  ADJUST GLOSS WITH DIAGONAL ELEMENT  **************** 
C
         DO 700 K    = 1, KTLOOP
            GLOSS(K,I) = GLOSS(K,I) * VDIAG(K,I)
 700     CONTINUE
 710  CONTINUE
C     
C *********************************************************************
C ******************** END OF SUBROUTINE BACKSUB **********************
C *********************************************************************
C
      RETURN
      END SUBROUTINE BACKSUB
