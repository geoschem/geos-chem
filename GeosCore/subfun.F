! $Id: subfun.f,v 1.1 2009/09/16 14:06:04 bmy Exp $
      SUBROUTINE SUBFUN
!
!******************************************************************************
!  Subroutine SUBFUN evaluates the first derivative of each ODE for SMVGEAR II.
!  (M. Jacobson, 1997; bdf, bmy, 4/1/03)
!
!  NOTES:
!  (1 ) Now force double-precision with the "D" exponent (bmy, 4/18/03)
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
!      SSSSSSS  U     U  BBBBBBB  FFFFFFF  U     U  N     N  
!      S        U     U  B     B  F        U     U  NN    N 
!      SSSSSSS  U     U  BBBBBBB  FFF      U     U  N  N  N 
!            S  U     U  B     B  F        U     U  N    NN
!      SSSSSSS  UUUUUUU  BBBBBBB  F        UUUUUUU  N     N 
!
! *********************************************************************
! *  THIS SUBROUTINE EVALUATES THE FIRST DERIVATIVE OF EACH ORDINARY  *  
! *                  DIFFERENTIAL EQUATION (ODE)                      * 
! *                                                                   *
! * HOW TO CALL SUBROUTINE:                                           *
! * ----------------------                                            *
! *  CALL SUBFUN.F FROM SMVGEAR.F WITH                                * 
! *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
! *     NCSP = NCS        FOR DAYTIME   GAS CHEM                      *  
! *     NCSP = NCS   +ICS FOR NIGHTTIME GAS CHEM                      *  
! *********************************************************************
!
! EXAMPLE
! -------
!
! SPECIES:         A,   B,   C
! CONCENTRATIONS: [A], [B], [C]
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
! *********************************************************************
!
! CONCMULT  = PRODUCT OF CONCENTRATIONS IN A RATE. IF TWO  
!             CONSECUTIVE REACTIONS HAVE THE SAME SPECIES REACTING
!             (EG A + B --> C AND A + B --> D + E) THEN USE THE 
!             SAME VALUE OF CONCMULT FOR BOTH REACTIONS.
! CNEW      = INIT (AND FINAL) SPECIES CONC (# CM-3-AIR OR MOLES L-1-H2O)
! GLOSS     = FIRST DERIVATIVE = SUM OF PROD. MINUS LOSS RATES FOR A SPECIES
! IRMA,B,C  = LOCATES REORDERED ACTIVE SPECIES NUMBERS  
! ISCHAN    = NUMBER OF ODES.
! LOSSRA..  = REAORDERED REACTION RATE NUMBERS FOR EACH LOSS (AND PROD) TERM
! KTLOOP    = NUMBER OF GRID-CELLS IN A GRID-BLOCK
! NSUBFUN   = COUNTS THE NUMBER OF TIMES THIS ROUTINE IS CALLED
! RRATE     = FORWARD RATE COEFFICIENT 
!           = S-1                                 FOR RATES WITH 1 REACTANT
!           = L-H2O MOLE-1 S-1  OR CM**3 #-1 S-1  FOR RATES WITH 2 REACTANTS 
!           = L**2-H2O M-2 S-1  OR CM**6 #-2 S-1  FOR RATES WITH 3 REACTANTS 
! TRATE     = REACTION RATE  MOLES L-1 -H2O S-1 OR # CM-3 S-1 
! 
! *********************************************************************
! *                      SET RATES OF REACTION                        *
! *********************************************************************
!
!
      ! Local variables
      INTEGER NKN,JA,JB,JC,NH,K,NK2,NH2,JSPC,NPL,NL5,NH5,NL4,NH4,NL3,NH3
      INTEGER NL2,NL1,NH1,NC,NK0,NK1,NK3,NK4,N
      INTEGER NK,I,JNEW,KLOOP
      REAL*8 CONCMULT,FRACN

      NSUBFUN        = NSUBFUN + 1
      NFDH1          = NFDH2 + IONER(NCSP) 
!
! *********************************************************************
! *     FIRST DERIVATIVES FOR RATES WITH THREE ACTIVE LOSS TERMS      *
! *********************************************************************
!

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

!
! *********************************************************************
! *     FIRST DERIVATIVES FOR RATES WITH TWO ACTIVE LOSS TERMS        *
! *********************************************************************
!

      DO 152 NKN     = NFDL2, NFDREP
       JA            = IRMA(NKN)
       JB            = IRMB(NKN)
       NH            = NKN + NALLR
       DO 150 K      = 1, KTLOOP
        TRATE(K,NKN) = RRATE(K,NKN) * CNEW(K,JA) * CNEW(K,JB) 
        TRATE(K,NH) = -TRATE(K,NKN)
 150   CONTINUE 
 152  CONTINUE 

!
! *********************************************************************
! *     FIRST DERIVATIVES FOR RATES WITH TWO ACTIVE LOSS TERMS AND    *
! *     WHERE THE SUBSEQUENT REACTION HAS THE SAME REACTANTS BUT A    *
! *     DIFFERENT RATE.                                               *
! *********************************************************************
!
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

!
! *********************************************************************
! *     FIRST DERIVATIVES FOR RATES WITH ONE ACTIVE LOSS TERM         *
! *********************************************************************
!
      DO 252 NKN     = NFDL1, NFDH1
       JA            = IRMA(NKN)
       NH            = NKN + NALLR
       DO 250 K      = 1, KTLOOP
        TRATE(K,NKN) = RRATE(K,NKN) * CNEW(K,JA) 
        TRATE(K,NH) = -TRATE(K,NKN) 
 250   CONTINUE 
 252  CONTINUE 

!
! *********************************************************************
! *                  INITIALIZE FIRST DERIVATIVE = 0                  *
! *********************************************************************
!

      DO 302 JSPC      = 1, ISCHAN
       DO 300 K        = 1, KTLOOP
        GLOSS(K,JSPC)  = 0.d0
 300   CONTINUE
 302  CONTINUE

!
! *********************************************************************
! * SUM NET (NOT REPRODUCED) KINETIC AND PHOTO GAINS AND LOSSES FOR   *
! * EACH SPECIES.                                                     *
! *********************************************************************
! SUM 1,2,3,4, OR 5 TERMS AT A TIME TO IMPROVE VECTORIZATION.
!
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
!
! ***********************  SUM 5 TERMS AT A TIME  ********************* 
!
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

!
! ***********************  SUM 4 TERMS AT A TIME  ********************* 
!
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
!
! ***********************  SUM 3 TERMS AT A TIME  ********************* 
!
       DO 452 NC       = NL3, NH3  
        NK0            = LOSSRA(NC)
        NK1            = LOSSRB(NC)
        NK2            = LOSSRC(NC)
        DO 450 K       = 1, KTLOOP
         GLOSS(K,JSPC) = GLOSS(K,JSPC) - TRATE(K,NK0)       
     1                 - TRATE(K,NK1)  - TRATE(K,NK2)
 450    CONTINUE
 452   CONTINUE
!
! ***********************  SUM 2 TERMS AT A TIME  ********************* 
!
       DO 502 NC       = NL2, NH2   
        NK0            = LOSSRA(NC)
        NK1            = LOSSRB(NC)
        DO 500 K       = 1, KTLOOP
         GLOSS(K,JSPC) = GLOSS(K,JSPC) - TRATE(K,NK0) 
     1                 - TRATE(K,NK1)      
 500    CONTINUE
 502   CONTINUE
!
! ***********************  SUM 1 TERM AT A TIME  ********************** 
!
       DO 552 NC       = NL1, NH1    
        NK0            = LOSSRA(NC)
        DO 550 K       = 1, KTLOOP
         GLOSS(K,JSPC) = GLOSS(K,JSPC) - TRATE(K,NK0)       
 550    CONTINUE
 552   CONTINUE
 554  CONTINUE
!
! *********************************************************************
! *  SUM PRODUCTION TERM FOR REACTIONS WHERE PRODUCTS FRACTIONATED    *
! *********************************************************************
!
      DO 802 N         = NFRLO(NCSP), NFRHI(NCSP)
       JSPC            = JSPCNFR(N)
       NKN             = NKNFR(  N)
       FRACN           = FRACNFR(N)
       DO 800 K        = 1, KTLOOP
        GLOSS(K,JSPC)  = GLOSS(K,JSPC) + FRACN * TRATE(K,NKN)       
 800   CONTINUE 
 802  CONTINUE 

!
! *********************************************************************
! **********************  END OF SUBROUTINE SUBFUN  *******************
! *********************************************************************
!
      RETURN
      END SUBROUTINE SUBFUN
