! $Id: subfun.f,v 1.1 2009/09/16 14:06:04 bmy Exp $
      SUBROUTINE SUBFUN
!
!******************************************************************************
!  Subroutine SUBFUN evaluates the first derivative of each ODE for SMVGEAR II.
!  (M. Jacobson, 1997; bdf, bmy, 4/1/03)
!
!  NOTES:
!  (1 ) Now force double-precision with the "D" exponent (bmy, 4/18/03)
!******************************************************************************
!
      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! SMVGEAR II arrays
C
C *********************************************************************
C ************        WRITTEN BY MARK JACOBSON (1993)      ************
C ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
C ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
C ***                         (650) 723-6836                        *** 
C *********************************************************************
C
C      SSSSSSS  U     U  BBBBBBB  FFFFFFF  U     U  N     N  
C      S        U     U  B     B  F        U     U  NN    N 
C      SSSSSSS  U     U  BBBBBBB  FFF      U     U  N  N  N 
C            S  U     U  B     B  F        U     U  N    NN
C      SSSSSSS  UUUUUUU  BBBBBBB  F        UUUUUUU  N     N 
C
C *********************************************************************
C *  THIS SUBROUTINE EVALUATES THE FIRST DERIVATIVE OF EACH ORDINARY  *  
C *                  DIFFERENTIAL EQUATION (ODE)                      * 
C *                                                                   *
C * HOW TO CALL SUBROUTINE:                                           *
C * ----------------------                                            *
C *  CALL SUBFUN.F FROM SMVGEAR.F WITH                                * 
C *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
C *     NCSP = NCS        FOR DAYTIME   GAS CHEM                      *  
C *     NCSP = NCS   +ICS FOR NIGHTTIME GAS CHEM                      *  
C *********************************************************************
C
C EXAMPLE
C -------
C
C SPECIES:         A,   B,   C
C CONCENTRATIONS: [A], [B], [C]
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
C *********************************************************************
C
C CONCMULT  = PRODUCT OF CONCENTRATIONS IN A RATE. IF TWO  
C             CONSECUTIVE REACTIONS HAVE THE SAME SPECIES REACTING
C             (EG A + B --> C AND A + B --> D + E) THEN USE THE 
C             SAME VALUE OF CONCMULT FOR BOTH REACTIONS.
C CNEW      = INIT (AND FINAL) SPECIES CONC (# CM-3-AIR OR MOLES L-1-H2O)
C GLOSS     = FIRST DERIVATIVE = SUM OF PROD. MINUS LOSS RATES FOR A SPECIES
C IRMA,B,C  = LOCATES REORDERED ACTIVE SPECIES NUMBERS  
C ISCHAN    = NUMBER OF ODES.
C LOSSRA..  = REAORDERED REACTION RATE NUMBERS FOR EACH LOSS (AND PROD) TERM
C KTLOOP    = NUMBER OF GRID-CELLS IN A GRID-BLOCK
C NSUBFUN   = COUNTS THE NUMBER OF TIMES THIS ROUTINE IS CALLED
C RRATE     = FORWARD RATE COEFFICIENT 
C           = S-1                                 FOR RATES WITH 1 REACTANT
C           = L-H2O MOLE-1 S-1  OR CM**3 #-1 S-1  FOR RATES WITH 2 REACTANTS 
C           = L**2-H2O M-2 S-1  OR CM**6 #-2 S-1  FOR RATES WITH 3 REACTANTS 
C TRATE     = REACTION RATE  MOLES L-1 -H2O S-1 OR # CM-3 S-1 
C 
C *********************************************************************
C *                      SET RATES OF REACTION                        *
C *********************************************************************
C
C
      ! Local variables
      INTEGER NKN,JA,JB,JC,NH,K,NK2,NH2,JSPC,NPL,NL5,NH5,NL4,NH4,NL3,NH3
      INTEGER NL2,NL1,NH1,NC,NK0,NK1,NK3,NK4,N
      INTEGER NK,I,JNEW,KLOOP
      REAL*8 CONCMULT,FRACN

      NSUBFUN        = NSUBFUN + 1
      NFDH1          = NFDH2 + IONER(NCSP) 
C
C *********************************************************************
C *     FIRST DERIVATIVES FOR RATES WITH THREE ACTIVE LOSS TERMS      *
C *********************************************************************
C

      DO 102 NKN     = 1, NFDH3  
       JA            = IRMA(NKN)
       JB            = IRMB(NKN)
       JC            = IRMC(NKN)
       NH            = NKN + NALLR
       DO 100 K      = 1, KTLOOP
        TRATE(K,NKN) = RRATE(K,NKN)*CNEW(K,JA)*CNEW(K,JB)*CNEW(K,JC) 
        TRATE(K,NH) = -TRATE(K,NKN) 
 100   CONTINUE
 102  CONTINUE

C
C *********************************************************************
C *     FIRST DERIVATIVES FOR RATES WITH TWO ACTIVE LOSS TERMS        *
C *********************************************************************
C

      DO 152 NKN     = NFDL2, NFDREP
       JA            = IRMA(NKN)
       JB            = IRMB(NKN)
       NH            = NKN + NALLR
       DO 150 K      = 1, KTLOOP
        TRATE(K,NKN) = RRATE(K,NKN) * CNEW(K,JA) * CNEW(K,JB) 
        TRATE(K,NH) = -TRATE(K,NKN)
 150   CONTINUE 
 152  CONTINUE 

C
C *********************************************************************
C *     FIRST DERIVATIVES FOR RATES WITH TWO ACTIVE LOSS TERMS AND    *
C *     WHERE THE SUBSEQUENT REACTION HAS THE SAME REACTANTS BUT A    *
C *     DIFFERENT RATE.                                               *
C *********************************************************************
C
      DO 202 NKN     = NFDREP1, NFDH2, 2
       JA            = IRMA(NKN)
       JB            = IRMB(NKN)
       NK2           = NKN + 1
       NH            = NKN + NALLR
       NH2           = NK2 + NALLR
       DO 200 K      = 1, KTLOOP
        CONCMULT     = CNEW(K,JA)   * CNEW(K,JB) 
        TRATE(K,NKN) = RRATE(K,NKN) * CONCMULT
        TRATE(K,NK2) = RRATE(K,NK2) * CONCMULT
        TRATE(K,NH)  = -TRATE(K,NKN) 
        TRATE(K,NH2) = -TRATE(K,NK2) 
 200   CONTINUE 
 202  CONTINUE 

C
C *********************************************************************
C *     FIRST DERIVATIVES FOR RATES WITH ONE ACTIVE LOSS TERM         *
C *********************************************************************
C
      DO 252 NKN     = NFDL1, NFDH1
       JA            = IRMA(NKN)
       NH            = NKN + NALLR
       DO 250 K      = 1, KTLOOP
        TRATE(K,NKN) = RRATE(K,NKN) * CNEW(K,JA) 
        TRATE(K,NH) = -TRATE(K,NKN) 
 250   CONTINUE 
 252  CONTINUE 

C
C *********************************************************************
C *                  INITIALIZE FIRST DERIVATIVE = 0                  *
C *********************************************************************
C

      DO 302 JSPC      = 1, ISCHAN
       DO 300 K        = 1, KTLOOP
        GLOSS(K,JSPC)  = 0.d0
 300   CONTINUE
 302  CONTINUE

C
C *********************************************************************
C * SUM NET (NOT REPRODUCED) KINETIC AND PHOTO GAINS AND LOSSES FOR   *
C * EACH SPECIES.                                                     *
C *********************************************************************
C SUM 1,2,3,4, OR 5 TERMS AT A TIME TO IMPROVE VECTORIZATION.
C
      DO 554 NPL       = NPLLO(NCSP), NPLHI(NCSP)
       JSPC            = JSPNPL(NPL)
       NL5             = NPL5(  NPL)
       NH5             = NPH5(  NPL)
       NL4             = NPL4(  NPL)
       NH4             = NPH4(  NPL)
       NL3             = NPL3(  NPL)
       NH3             = NPH3(  NPL)
       NL2             = NPL2(  NPL)
       NH2             = NPH2(  NPL)
       NL1             = NPL1(  NPL)
       NH1             = NPH1(  NPL)
C
C ***********************  SUM 5 TERMS AT A TIME  ********************* 
C
       DO 352 NC       = NL5, NH5
        NK0            = LOSSRA(NC)
        NK1            = LOSSRB(NC)
        NK2            = LOSSRC(NC)
        NK3            = LOSSRD(NC)
        NK4            = LOSSRE(NC)
        DO 350 K       = 1, KTLOOP
         GLOSS(K,JSPC) = GLOSS(K,JSPC) - TRATE(K,NK0)       
     1                 - TRATE(K,NK1)  - TRATE(K,NK2)
     2                 - TRATE(K,NK3)  - TRATE(K,NK4)  
 350    CONTINUE
 352   CONTINUE

C
C ***********************  SUM 4 TERMS AT A TIME  ********************* 
C
       DO 402 NC       = NL4, NH4 
        NK0            = LOSSRA(NC)
        NK1            = LOSSRB(NC)
        NK2            = LOSSRC(NC)
        NK3            = LOSSRD(NC)
        DO 400 K       = 1, KTLOOP
         GLOSS(K,JSPC) = GLOSS(K,JSPC) - TRATE(K,NK0)       
     1                 - TRATE(K,NK1)  - TRATE(K,NK2)
     2                 - TRATE(K,NK3)  
 400    CONTINUE
 402   CONTINUE
C
C ***********************  SUM 3 TERMS AT A TIME  ********************* 
C
       DO 452 NC       = NL3, NH3  
        NK0            = LOSSRA(NC)
        NK1            = LOSSRB(NC)
        NK2            = LOSSRC(NC)
        DO 450 K       = 1, KTLOOP
         GLOSS(K,JSPC) = GLOSS(K,JSPC) - TRATE(K,NK0)       
     1                 - TRATE(K,NK1)  - TRATE(K,NK2)
 450    CONTINUE
 452   CONTINUE
C
C ***********************  SUM 2 TERMS AT A TIME  ********************* 
C
       DO 502 NC       = NL2, NH2   
        NK0            = LOSSRA(NC)
        NK1            = LOSSRB(NC)
        DO 500 K       = 1, KTLOOP
         GLOSS(K,JSPC) = GLOSS(K,JSPC) - TRATE(K,NK0) 
     1                 - TRATE(K,NK1)      
 500    CONTINUE
 502   CONTINUE
C
C ***********************  SUM 1 TERM AT A TIME  ********************** 
C
       DO 552 NC       = NL1, NH1    
        NK0            = LOSSRA(NC)
        DO 550 K       = 1, KTLOOP
         GLOSS(K,JSPC) = GLOSS(K,JSPC) - TRATE(K,NK0)       
 550    CONTINUE
 552   CONTINUE
 554  CONTINUE
C
C *********************************************************************
C *  SUM PRODUCTION TERM FOR REACTIONS WHERE PRODUCTS FRACTIONATED    *
C *********************************************************************
C
      DO 802 N         = NFRLO(NCSP), NFRHI(NCSP)
       JSPC            = JSPCNFR(N)
       NKN             = NKNFR(  N)
       FRACN           = FRACNFR(N)
       DO 800 K        = 1, KTLOOP
        GLOSS(K,JSPC)  = GLOSS(K,JSPC) + FRACN * TRATE(K,NKN)       
 800   CONTINUE 
 802  CONTINUE 

C
C *********************************************************************
C **********************  END OF SUBROUTINE SUBFUN  *******************
C *********************************************************************
C
      RETURN
      END SUBROUTINE SUBFUN
