! $Id: pderiv.f,v 1.1 2009/09/16 14:06:16 bmy Exp $
      SUBROUTINE PDERIV
!
!******************************************************************************
!  Subroutine PDERIV places the partial differential equations into a matrix
!  for SMVGEAR II. (M. Jacobson, 1997; bdf, bmy, 4/18/03)
!
!  NOTES:
!  (1 ) Now force double-precision w/ "D" exponents (bmy, 4/18/03)
!******************************************************************************
!
#     include "define.h"
      USE CMN_SIZE_MOD  ! Size parameters
      USE COMODE_LOOP_MOD  ! SMVGEAR II arrays

      IMPLICIT NONE

C
C *********************************************************************
C ************        WRITTEN BY MARK JACOBSON (1993)      ************
C ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
C ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
C ***                         (650) 723-6836                        *** 
C *********************************************************************
C
C       PPPPPPP  DDDDDDD  EEEEEEE  RRRRRRR  IIIIIII  V       V 
C       P     P  D     D  E        R     R     I      V     V
C       PPPPPPP  D     D  EEEEEEE  RRRRRRR     I       V   V 
C       P        D     D  E        R  R        I        V V  
C       P        DDDDDDD  EEEEEEE  R    R   IIIIIII      V  
C
C *********************************************************************
C * THIS SUBROUTINE PUTS THE PARTIAL DERIVATIVES OF EACH ORDINARY     *  
C * DIFFERENTIAL EQUATION INTO A MATRIX. THE FORM OF THE MATRIX       *
C * EQUATION IS                                                       *
C *                       P = I - H x Bo x J                          *
C *                                                                   *
C * WHERE I = IDENTITY MATRIX, H = TIME-STEP, Bo = COEFFICIENT        *
C * CORRESPONDING TO THE ORDER OF THE METHOD, AND J IS THE JACOBIAN   *
C * MATRIX OF PARTIAL DERIVATIVES.                                    * 
C *                                                                   *
C * HOW TO CALL SUBROUTINE:                                           *
C * ----------------------                                            *
C *  CALL PDERIV.F FROM SMVGEAR.F WITH                                * 
C *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
C *     NCSP = NCS        FOR DAYTIME   GAS CHEM                      *  
C *     NCSP = NCS   +ICS FOR NIGHTTIME GAS CHEM                      *  
C *********************************************************************
C
C *********************************************************************
C *                        INITIALIZE MATRIX                          * 
C *********************************************************************
C CC2      = ARRAY OF IARRAY UNITS HOLDING VALUES OF EACH MAXTRIX
C            POSITION ACTUALLY USED.
C            CC2 = P = I - DELT * ASET(NQQ,1) * PARTIAL DERIVATIVES.   
C URATE    = TERM OF JACOBIAN (J) = PARTIAL DERIVATIVE 
C IARRAY   = TOTAL NUMBER OF MATRIX POSITIONS FILLED AFTER MAT. PROCESSES
C IRMA,B,C = SPECIES # OF EACH REACTANT 
C ISCHAN   = ORIGINAL ORDER OF MATRIX  
C KTLOOP   = NUMBER OF GRID-CELLS IN A GRID-BLOCK
C NONDIAG1 = 1 + # OF FINAL MATRIX POSITIONS, EXCLUDING DIAGONAL TERMS,
C            FILLED AFTER ALL MATRIX PROCESSES.
C NPDERIV  = COUNTER OF NUMBER OF TIMES THIS ROUTINE IS CALLED
C R1DELT   = -ASET(NQQ,1) * TIME STEP = -COEFFICIENT OF METHOD * DT 
C RRATE    = REACTION RATE COEFFICIENT
C
C EXAMPLE OF HOW PARTIAL DERIVATIVES ARE PLACED IN AN ARRAY:
C ---------------------------------------------------------
C
C SPECIES:          A,   B,   C
C CONCENTRATIONS:  [A], [B], [C]       
C
C REACTIONS:    1) A           --> B      J 
C               2) A  + B      --> C      K1 
C               3  A  + B + C  --> D      K2  
C
C FIRST         d[A] / dt  =  -J[A] - K1[A][B] - K2[A][B][C]
C DERIVATIVES:  d[B] / dt  =  +J[A] - K1[A][B] - K2[A][B][C]
C               d[C] / dt  =        + K1[A][B] - K2[A][B][C]
C               d[D] / dt  =                   + K2[A][B][C]
C
C PREDICTOR MATRIX (P) = I - h * b * J:
C  J = JACOBIAN MATRIX OF PARTIAL DERIVATES
C  I = IDENTITY MATRIX
C  h = TIME-STEP
C  b = COEFFICIENT OF METHOD
C  R = h * b = -R1DELT  
C
C                A                 B                 C                 D 
C     ___________________________________________________________________
C   | 
C A | 1-R(-J-K1[B]-K2[B][C])  -R(-K1[A]-K2[A][C])   -R(-K2[A][B])       0    
C   | 
C B |  -R(+J-K1[B]-K2[B][C]) 1-R(-K1[A]-K2[A][C])   -R(-K2[A][B])       0  
C   | 
C C |  -R(  +K1[B]-K2[B][C])  -R(+K1[A]-K2[A][C])  1-R(-K2[A][B])       0 
C   |
C D |  -R(        +K2[B][C])  -R(      +K2[A][C])   -R(+K2[A][B])       1 
C    
C
C *********************************************************************
C *********         CALCULATE PARTIAL DERIVATIVES            **********
C *********    AND SUM UP PARTIAL DERIVATIVE LOSS TERMS      ********** 
C *********************************************************************
C
      INTEGER IARRY,NONDIAG,NONDIAG1,NPDL,NPDH,NKN,JA,JB,JC,K,IAR,N
      INTEGER IAL
      REAL*8 FRACR1

      NPDERIV             = NPDERIV + 1
      IARRY               = IARRAY(NCSP) 
      NONDIAG             = IARRY - ISCHAN  
      NONDIAG1            = NONDIAG + 1
      NFDH1               = NFDH2 + IONER(NCSP) 
      NPDL                = NPDLO(NCSP)
      NPDH                = NPDHI(NCSP)
C
C *********************************************************************
C *    PARTIAL DERIVATIVES FOR RATES WITH THREE ACTIVE LOSS TERMS     *
C *********************************************************************
C
      DO 105 NKN       = 1, NFDH3
         JA              = IRMA(NKN)
         JB              = IRMB(NKN)
         JC              = IRMC(NKN)
         DO 100 K        = 1, KTLOOP
            URATE(K,NKN,1) = RRATE(K,NKN) * CNEW(K,JB) * CNEW(K,JC) 
            URATE(K,NKN,2) = RRATE(K,NKN) * CNEW(K,JA) * CNEW(K,JC)
            URATE(K,NKN,3) = RRATE(K,NKN) * CNEW(K,JA) * CNEW(K,JB)
 100     CONTINUE 
 105  CONTINUE 
C
C *********************************************************************
C *    PARTIAL DERIVATIVES FOR RATES WITH TWO ACTIVE LOSS TERMS       *
C *********************************************************************
C
      DO 155 NKN       = NFDL2, NFDH2
         JA              = IRMA(NKN)
         JB              = IRMB(NKN)
         DO 150 K        = 1, KTLOOP
            URATE(K,NKN,1) = RRATE(K,NKN) * CNEW(K,JB) 
            URATE(K,NKN,2) = RRATE(K,NKN) * CNEW(K,JA) 
 150     CONTINUE 
 155  CONTINUE 
C
C *********************************************************************
C *    PARTIAL DERIVATIVES FOR RATES WITH ONE ACTIVE LOSS TERM        *
C *********************************************************************
C
      DO 205 NKN       = NFDL1, NFDH1 
         DO 200 K        = 1, KTLOOP
            URATE(K,NKN,1) = RRATE(K,NKN) 
 200     CONTINUE
 205  CONTINUE
C
C *********************************************************************
C * PUT PARTIAL DERIVATIVES PRODUCTION AND LOSS TERMS IN MATRIX ARRAY * 
C *********************************************************************
C FRACPL = -1. FOR ALL REACTANTS 
C        = +1. OR +FRACTION FOR ALL PRODUCTS
C   
      DO 255 IAR     = 1, NONDIAG
         DO 250 K      = 1, KTLOOP
            CC2(K,IAR)   = 0.d0
 250     CONTINUE
 255  CONTINUE
C
      DO 305 IAR     = NONDIAG1, IARRY
         DO 300 K      = 1, KTLOOP
            CC2(K,IAR)   = 1.d0
 300     CONTINUE
 305  CONTINUE
C
      DO 405 N       = NPDL, NPDH
         NKN           = NKPDTERM(N) 
         IAR           = IPOSPD(  N)
         IAL           = IIALPD(  N)
         FRACR1        = FRACPL(  N) * R1DELT  
         DO 400 K      = 1, KTLOOP 
            CC2(K,IAR)   = CC2(K,IAR) + FRACR1 * URATE(K,NKN,IAL)   
 400     CONTINUE 
 405  CONTINUE 
C
C *********************************************************************
C ********************* END OF SUBROUTINE PDERIV **********************
C *********************************************************************
C
      RETURN 
      END SUBROUTINE PDERIV
