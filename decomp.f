! $Id: decomp.f,v 1.1 2003/06/30 20:26:09 bmy Exp $
      SUBROUTINE DECOMP
!
!******************************************************************************
!  Subroutine DECOMP decomposes the sparse matrix for the SMVGEAR II solver.
!  (M. Jacobson, 1997; bdf, bmy, 4/18/03)
!
!  NOTES:
!  (1 ) Now use & as F90 continuation character.  Now also force double
!        precision with the "D" exponent. (bmy, 4/18/03)
!******************************************************************************
!
      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "comode.h"
C
C *********************************************************************
C ************        WRITTEN BY MARK JACOBSON (1993)      ************
C ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
C ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
C ***                         (650) 723-6836                        *** 
C *********************************************************************
C
C         DDDDDDD  EEEEEEE  CCCCCCC  OOOOOOO  M     M  PPPPPPP 
C         D     D  E        C        O     O  MM   MM  P     P 
C         D     D  EEEEEEE  C        O     O  M M M M  PPPPPPP 
C         D     D  E        C        O     O  M  M  M  P  
C         DDDDDDD  EEEEEEE  CCCCCCC  OOOOOOO  M     M  P
C 
C *********************************************************************
C ************** DECOMPOSE THE SPARSE MATRIX **************************
C *********************************************************************
C
C *********************************************************************
C * THIS SUBROUTINE DECOMPOSES THE MATRIX "P" INTO THE MATRIX "A" IN  *
C * ORDER TO SOLVE THE LINEAR SET OF EQUATIONS Ax = B FOR x, WHICH IS *
C * A CORRECTION VECTOR. Ax = B IS SOLVED IN SUBROUTINE BACKSUB.F     *
C * ABOVE, THE ORIGINAL MATRIX "P" IS                                 *
C *                                                                   *
C *                       P = I - H x Bo x J,                         * 
C *                                                                   *
C * WHERE I = IDENTITY MATRIX, H = TIME-STEP, Bo = A COEFFICIENT THAT * 
C * DEPENDS ON THE ORDER OF THE INTEGRATION METHOD, AND J IS THE      * 
C * MATRIX OF PARTIAL DERIVATIVES. SEE PRESS ET AL. (1992) NUMERICAL  *
C * RECIPES CAMBRIDGE UNIVERSITY PRESS, FOR A BETTER DESCRIPTION OF   *
C * THE L-U DECOMPOSTION PROCESS                                      *
C *                                                                   *
C * THIS L-U DECOMPOSTION PROCESS USES SPARSE-MATRIX TECHNIQUES,      *
C * VECTORIZES AROUND THE GRID-CELL DIMENSION, AND USES NO PARTIAL    *
C * PIVOTING. TESTS BY SHERMAN & HINDMARSH (1980) LAWRENCE LIVERMORE  * 
C * REP. UCRL-84102 AND BY US HAVE CONFIRMED THAT THE REMOVAL OF      *
C * PARTIAL PIVOTING HAS LITTLE EFFECT ON RESULTS.                    *
C *                                                                   *
C * HOW TO CALL SUBROUTINE:                                           *
C * ----------------------                                            *
C *  CALL DECOMP.F FROM SMVGEAR.F WITH                                * 
C *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
C *     NCSP = NCS        FOR DAYTIME   GAS CHEM                      *  
C *     NCSP = NCS   +ICS FOR NIGHTTIME GAS CHEM                      *  
C *********************************************************************
C
C KTLOOP   = NUMBER OF GRID-CELLS IN A GRID-BLOCK
C ISCHAN   = ORIGINAL ORDER OF MATRIX  
C CC2      = ARRAY OF IARRAY UNITS HOLDING VALUES OF EACH MATRIX
C            POSITION ACTUALLY USED. ORIGINALLY,
C            CC2 = P = I - DELT * ASET(NQQ,1) * PARTIAL DERIVATIVES.   
C            HOWEVER, CC2 IS DECOMPOSED HERE
C
C *********************************************************************
C ***                  FIRST LOOP OF L-U DECOMPOSTION               *** 
C *********************************************************************
C SUM 1,2,3,4, OR 5 TERMS AT A TIME TO IMPROVE VECTORIZATION 
C
      INTEGER J,IJT,IJ,IL5,IH5,IL4,IH4,IL3,IH3,IL2,IH2,IL1,IH1
      INTEGER IC,IK0,IK1,IK2,IK3,IK4,KJ0,KJ1,KJ2,KJ3,KJ4,K,IAR
      INTEGER JL,JH,JC,IJA

      ! bdf timing calculations
      NUM_DECOMP = NUM_DECOMP + 1

      DO 510 J       = 1, ISCHAN
         DO 310 IJT    = IJTLO(J,NCSP), IJTHI(J,NCSP)
            IJ           = IJVAL(IJT)
            IL5          = IDL5( IJT)
            IH5          = IDH5( IJT)
            IL4          = IDL4( IJT)  
            IH4          = IDH4( IJT)
            IL3          = IDL3( IJT)  
            IH3          = IDH3( IJT)
            IL2          = IDL2( IJT)  
            IH2          = IDH2( IJT)
            IL1          = IDL1( IJT)  
            IH1          = IDH1( IJT)

C ********************* SUM 5 TERMS AT A TIME *************************
C
            DO 105 IC    = IL5, IH5
               IK0         = IKDECA(IC)
               IK1         = IKDECB(IC)
               IK2         = IKDECC(IC)
               IK3         = IKDECD(IC)
               IK4         = IKDECE(IC)
               KJ0         = KJDECA(IC)
               KJ1         = KJDECB(IC)
               KJ2         = KJDECC(IC)
               KJ3         = KJDECD(IC)
               KJ4         = KJDECE(IC)
               DO 100 K    = 1, KTLOOP 
                  CC2(K,IJ)  = CC2(K,IJ) - CC2(K,IK0) * CC2(K,KJ0) 
     &                       - CC2(K,IK1) * CC2(K,KJ1) 
     &                       - CC2(K,IK2) * CC2(K,KJ2) 
     &                       - CC2(K,IK3) * CC2(K,KJ3) 
     &                       - CC2(K,IK4) * CC2(K,KJ4) 
 100           CONTINUE 
 105        CONTINUE 
C
C ********************* SUM 4 TERMS AT A TIME *************************
C
            DO 155 IC    = IL4, IH4
               IK0         = IKDECA(IC)
               IK1         = IKDECB(IC)
               IK2         = IKDECC(IC)
               IK3         = IKDECD(IC)
               KJ0         = KJDECA(IC)
               KJ1         = KJDECB(IC)
               KJ2         = KJDECC(IC)
               KJ3         = KJDECD(IC)
               DO 150 K    = 1, KTLOOP 
                  CC2(K,IJ)  = CC2(K,IJ) - CC2(K,IK0) * CC2(K,KJ0) 
     &                       - CC2(K,IK1) * CC2(K,KJ1) 
     &                       - CC2(K,IK2) * CC2(K,KJ2) 
     &                       - CC2(K,IK3) * CC2(K,KJ3) 
 150           CONTINUE 
 155        CONTINUE 
C
C ********************* SUM 3 TERMS AT A TIME *************************
C
            DO 205 IC    = IL3, IH3
               IK0         = IKDECA(IC)
               IK1         = IKDECB(IC)
               IK2         = IKDECC(IC)
               KJ0         = KJDECA(IC)
               KJ1         = KJDECB(IC)
               KJ2         = KJDECC(IC)
               DO 200 K    = 1, KTLOOP 
                  CC2(K,IJ)  = CC2(K,IJ) - CC2(K,IK0) * CC2(K,KJ0) 
     &                       - CC2(K,IK1) * CC2(K,KJ1) 
     &                       - CC2(K,IK2) * CC2(K,KJ2) 
 200           CONTINUE
 205        CONTINUE
C
C ********************* SUM 2 TERMS AT A TIME *************************
C
            DO 255 IC    = IL2, IH2
               IK0         = IKDECA(IC)
               IK1         = IKDECB(IC)
               KJ0         = KJDECA(IC)
               KJ1         = KJDECB(IC)
               DO 250 K    = 1, KTLOOP 
                  CC2(K,IJ)  = CC2(K,IJ) - CC2(K,IK0) * CC2(K,KJ0) 
     &                       - CC2(K,IK1) * CC2(K,KJ1) 
 250           CONTINUE 
 255        CONTINUE 
C     
C ********************* SUM 1 TERM  AT A TIME *************************
C
            DO 305 IC    = IL1, IH1
               IK0         = IKDECA(IC)
               KJ0         = KJDECA(IC)
               DO 300 K    = 1, KTLOOP 
                  CC2(K,IJ)  = CC2(K,IJ) - CC2(K,IK0) * CC2(K,KJ0) 
 300           CONTINUE
 305        CONTINUE
C     
 310   CONTINUE
C     
C *********************************************************************
C *    VDIAG = 1 / CURRENT DIAGONAL TERM OF THE DECOMPOSED MATRIX     *
C *********************************************************************
C
       IAR          = JARRDIAG(J,NCSP)
       DO 400 K     = 1, KTLOOP 
          VDIAG(K,J)  = 1.0d0 / CC2(K,IAR) 
 400   CONTINUE
C
C *********************************************************************
C ***               SECOND LOOP OF DECOMPOSTION                     *** 
C *********************************************************************
C JZEROA  = IDENTIFIES THE ARRAY POSITION OF EACH JLOZ1..JHIZ1 TERM
C
       JL           = JLOZ1(J,NCSP)
       JH           = JHIZ1(J,NCSP)
       DO 505 JC    = JL, JH
          IJA         = JZEROA(JC)
          DO 500 K    = 1, KTLOOP 
             CC2(K,IJA) = CC2(K,IJA) * VDIAG(K,J)  
 500      CONTINUE
 505   CONTINUE
C
 510  CONTINUE
C
C *********************************************************************
C ********************* END OF SUBROUTINE DECOMP  ********************* 
C *********************************************************************
C
      RETURN
      END SUBROUTINE DECOMP
