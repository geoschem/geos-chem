! $Id: pderiv.f,v 1.1 2009/09/16 14:06:16 bmy Exp $
      SUBROUTINE PDERIV
!
!******************************************************************************
!  Subroutine PDERIV places the partial differential equations into a matrix
!  for SMVGEAR II. (M. Jacobson, 1997; bdf, bmy, 4/18/03)
!
!  NOTES:
!  (1 ) Now force double-precision w/ "D" exponents (bmy, 4/18/03)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE CMN_SIZE_MOD  ! Size parameters
      USE COMODE_LOOP_MOD  ! SMVGEAR II arrays

      IMPLICIT NONE

!
! *********************************************************************
! ************        WRITTEN BY MARK JACOBSON (1993)      ************
! ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
! ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
!       PPPPPPP  DDDDDDD  EEEEEEE  RRRRRRR  IIIIIII  V       V 
!       P     P  D     D  E        R     R     I      V     V
!       PPPPPPP  D     D  EEEEEEE  RRRRRRR     I       V   V 
!       P        D     D  E        R  R        I        V V  
!       P        DDDDDDD  EEEEEEE  R    R   IIIIIII      V  
!
! *********************************************************************
! * THIS SUBROUTINE PUTS THE PARTIAL DERIVATIVES OF EACH ORDINARY     *  
! * DIFFERENTIAL EQUATION INTO A MATRIX. THE FORM OF THE MATRIX       *
! * EQUATION IS                                                       *
! *                       P = I - H x Bo x J                          *
! *                                                                   *
! * WHERE I = IDENTITY MATRIX, H = TIME-STEP, Bo = COEFFICIENT        *
! * CORRESPONDING TO THE ORDER OF THE METHOD, AND J IS THE JACOBIAN   *
! * MATRIX OF PARTIAL DERIVATIVES.                                    * 
! *                                                                   *
! * HOW TO CALL SUBROUTINE:                                           *
! * ----------------------                                            *
! *  CALL PDERIV.F FROM SMVGEAR.F WITH                                * 
! *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
! *     NCSP = NCS        FOR DAYTIME   GAS CHEM                      *  
! *     NCSP = NCS   +ICS FOR NIGHTTIME GAS CHEM                      *  
! *********************************************************************
!
! *********************************************************************
! *                        INITIALIZE MATRIX                          * 
! *********************************************************************
! CC2      = ARRAY OF IARRAY UNITS HOLDING VALUES OF EACH MAXTRIX
!            POSITION ACTUALLY USED.
!            CC2 = P = I - DELT * ASET(NQQ,1) * PARTIAL DERIVATIVES.   
! URATE    = TERM OF JACOBIAN (J) = PARTIAL DERIVATIVE 
! IARRAY   = TOTAL NUMBER OF MATRIX POSITIONS FILLED AFTER MAT. PROCESSES
! IRMA,B,C = SPECIES # OF EACH REACTANT 
! ISCHAN   = ORIGINAL ORDER OF MATRIX  
! KTLOOP   = NUMBER OF GRID-CELLS IN A GRID-BLOCK
! NONDIAG1 = 1 + # OF FINAL MATRIX POSITIONS, EXCLUDING DIAGONAL TERMS,
!            FILLED AFTER ALL MATRIX PROCESSES.
! NPDERIV  = COUNTER OF NUMBER OF TIMES THIS ROUTINE IS CALLED
! R1DELT   = -ASET(NQQ,1) * TIME STEP = -COEFFICIENT OF METHOD * DT 
! RRATE    = REACTION RATE COEFFICIENT
!
! EXAMPLE OF HOW PARTIAL DERIVATIVES ARE PLACED IN AN ARRAY:
! ---------------------------------------------------------
!
! SPECIES:          A,   B,   C
! CONCENTRATIONS:  [A], [B], [C]       
!
! REACTIONS:    1) A           --> B      J 
!               2) A  + B      --> C      K1 
!               3  A  + B + C  --> D      K2  
!
! FIRST         d[A] / dt  =  -J[A] - K1[A][B] - K2[A][B][C]
! DERIVATIVES:  d[B] / dt  =  +J[A] - K1[A][B] - K2[A][B][C]
!               d[C] / dt  =        + K1[A][B] - K2[A][B][C]
!               d[D] / dt  =                   + K2[A][B][C]
!
! PREDICTOR MATRIX (P) = I - h * b * J:
!  J = JACOBIAN MATRIX OF PARTIAL DERIVATES
!  I = IDENTITY MATRIX
!  h = TIME-STEP
!  b = COEFFICIENT OF METHOD
!  R = h * b = -R1DELT  
!
!                A                 B                 C                 D 
!     ___________________________________________________________________
!   | 
! A | 1-R(-J-K1[B]-K2[B][C])  -R(-K1[A]-K2[A][C])   -R(-K2[A][B])       0    
!   | 
! B |  -R(+J-K1[B]-K2[B][C]) 1-R(-K1[A]-K2[A][C])   -R(-K2[A][B])       0  
!   | 
! C |  -R(  +K1[B]-K2[B][C])  -R(+K1[A]-K2[A][C])  1-R(-K2[A][B])       0 
!   |
! D |  -R(        +K2[B][C])  -R(      +K2[A][C])   -R(+K2[A][B])       1 
!    
!
! *********************************************************************
! *********         CALCULATE PARTIAL DERIVATIVES            **********
! *********    AND SUM UP PARTIAL DERIVATIVE LOSS TERMS      ********** 
! *********************************************************************
!
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
!
! *********************************************************************
! *    PARTIAL DERIVATIVES FOR RATES WITH THREE ACTIVE LOSS TERMS     *
! *********************************************************************
!
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
!
! *********************************************************************
! *    PARTIAL DERIVATIVES FOR RATES WITH TWO ACTIVE LOSS TERMS       *
! *********************************************************************
!
      DO 155 NKN       = NFDL2, NFDH2
         JA              = IRMA(NKN)
         JB              = IRMB(NKN)
         DO 150 K        = 1, KTLOOP
            URATE(K,NKN,1) = RRATE(K,NKN) * CNEW(K,JB) 
            URATE(K,NKN,2) = RRATE(K,NKN) * CNEW(K,JA) 
 150     CONTINUE 
 155  CONTINUE 
!
! *********************************************************************
! *    PARTIAL DERIVATIVES FOR RATES WITH ONE ACTIVE LOSS TERM        *
! *********************************************************************
!
      DO 205 NKN       = NFDL1, NFDH1 
         DO 200 K        = 1, KTLOOP
            URATE(K,NKN,1) = RRATE(K,NKN) 
 200     CONTINUE
 205  CONTINUE
!
! *********************************************************************
! * PUT PARTIAL DERIVATIVES PRODUCTION AND LOSS TERMS IN MATRIX ARRAY * 
! *********************************************************************
! FRACPL = -1. FOR ALL REACTANTS 
!        = +1. OR +FRACTION FOR ALL PRODUCTS
!   
      DO 255 IAR     = 1, NONDIAG
         DO 250 K      = 1, KTLOOP
            CC2(K,IAR)   = 0.d0
 250     CONTINUE
 255  CONTINUE
!
      DO 305 IAR     = NONDIAG1, IARRY
         DO 300 K      = 1, KTLOOP
            CC2(K,IAR)   = 1.d0
 300     CONTINUE
 305  CONTINUE
!
      DO 405 N       = NPDL, NPDH
         NKN           = NKPDTERM(N) 
         IAR           = IPOSPD(  N)
         IAL           = IIALPD(  N)
         FRACR1        = FRACPL(  N) * R1DELT  
         DO 400 K      = 1, KTLOOP 
            CC2(K,IAR)   = CC2(K,IAR) + FRACR1 * URATE(K,NKN,IAL)   
 400     CONTINUE 
 405  CONTINUE 
!
! *********************************************************************
! ********************* END OF SUBROUTINE PDERIV **********************
! *********************************************************************
!
      RETURN 
      END SUBROUTINE PDERIV
