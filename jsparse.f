! $Id: jsparse.f,v 1.2 2003/07/08 15:32:25 bmy Exp $
      SUBROUTINE JSPARSE
!
!******************************************************************************
!  Subroutine JSPARSE sets up the sparse-matrix arrays for SMVGEAR II.
!  (M. Jacobson 1993; bdf, bmy, 4/18/03)
!
!  NOTES:
!  (1 ) For GEOS-CHEM we had to remove T3 from "comode.h" and to declare it
!        allocatable in "comode_mod.f".  This allows us to only allocate it
!        if we are doing a fullchem run.  Write list of repeat reactants to 
!        and change in moles to "smv2.log".  Now call GEOS_CHEM_STOP to
!        deallocate all arrays and stop the run safely.  Now force double
!        precision with "D" exponents. (bmy, 4/18/03)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD, ONLY : T3
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP

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
C        JJ  SSSSSSS  PPPPPPP     A      RRRRRRR  SSSSSSS  EEEEEEE
C         J  S        P     P    A A     R     R  S        E 
C         J  SSSSSSS  PPPPPPP   A   A    RRRRRRR  SSSSSSS  EEEEEEE
C   J     J        S  P        AAAAAAA   R  R           S  E
C   JJJJJJJ  SSSSSSS  P       A       A  R    R   SSSSSSS  EEEEEEE  
C
C *********************************************************************
C * THIS ROUTINE SETS UP SPARSE-MATRIX AND OTHER ARRAYS FOR SMVGEAR   *
C * (SPARSE-MATRIX VECTORIZED GEAR-CODE. IT SETS ARRAYS FOR GAS-      *
C * -PHASE, AQUEOUS-PHASE, AND ANY OTHER TYPE OF CHEMISTRY. IT ALSO   *
C * SETS ARRAYS FOR BOTH DAY AND NIGHT CHEMISTRY OF EACH TYPE.        *
C *                                                                   *
C * HOW TO CALL SUBROUTINE:                                           *
C * ----------------------                                            *
C *  CALL JSPARSE.F FROM READCHEM.F WITH                              * 
C *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
C *********************************************************************
C
C *********************************************************************
C ******* SETS UP ARRAYS FOR GAS- AND AQUEOUS-PHASE CHEMISTRY  ******** 
C * INCLUDES ARRAYS FOR CALCULATING FIRST DERIVATIVES, PARTIAL DERIV- *
C * ATIVES, MATRIX DECOMPOSTION, AND MATRIX BACK-SUBSTITUTION. FIRST, *
C * JSPARSE RE-ORDERS THE ORDINARY DIFFERENTIAL EQUATIONS TO MAXIMIZE *
C * THE NUMBER OF ZEROS IN THE MATRIX OF PARTIAL DERIVATIVES. IT      *
C * LATER SETS ARRAYS TO ELIMINATE ALL CALCULATIONS INVOLVING A ZERO. * 
C *********************************************************************
* 
C NTSPEC    = TOTAL NUMBER OF ACTIVE + INACTIVE SPECIES.
C NSPEC     = TOTAL NUMBER OF ACTIVE SPECIES.
C NMREAC    = 3 = MAXIMUM NUMBER OF ACTIVE REACTANTS IN A REACTION 
C NALLREAC  = 4 = TOTAL REACTANT POSITIONS IN A REACTION 
C NMPROD    = 5 = MAXIMUN NUMBER OF ACTIVE PRODUCTS IN A REACTION 
C NPRODLO   = NALLREAC  + 1 = LOWEST PRODUCT POSITION NUMBER. 
C NPRODHI   = NALLREAC + NMPROD = HIGHEST PRODUCT POSITION NUMBER. 
C
C *********************************************************************
C * DETERMINE HOW MANY PARTIAL DERIV TERMS ARE NEEDED FOR EACH SPECIES*
C *********************************************************************
C IFREPRO   = 1 THEN SPECIES IS LOST AND REPRODUCED IN REACTION NK 
C IRM       = SPECIES # OF EACH REACT OR PRODUCT IN EACH NK REACTION
C ISAPORL   = COUNTS PARTIAL DERIVATIVE TERMS FOR EACH SPECIES
C FKOEF     = 1, 2, FRACTION, OR MORE = # OF A GIVEN REACTANT OR PRODUCTS
C             E.G. REACTION      A + B  --> 2C + 0.34D  + D 
C             VALUE OF FKOEF     1   1      2    0.34     1     
C NCS       = 1..NCSGAS FOR GAS CHEMISTRY                            
C NCSP      = NCS        FOR DAYTIME   GAS CHEM 
C           = NCS   +ICS FOR NIGHTTIME GAS CHEM            
C NK        = REACTION # OF EACH REACTION 
C NRATES    = NUMBER OF KINETIC (NON-PHOTO) RATE COEFFICIENTS
C NTRATES   = NUMBER OF KINETIC PLUS PHOTO  RATE COEFFICIENTS
C NALLRAT   = NUMBER OF KINETIC PLUS PHOTO REACTION RATES  
C

      INTEGER NREPT,I,J,NAR,NK,K,IREACT,L,IPO,NOCHANG,JOLD,JNEW
      INTEGER MINVALU,IMINOLD,IMINNEW,INEW,IOLD,NKLAST,IAL,IRE
      INTEGER NMO,NOL,ISDIFF,IB,JSPCL,ISPC1,ISPC2,ISPC3,IAP,IPROD
      INTEGER IPR,LFRAC,NGN,KPRODS,KDIF,NPL,IC,NK1,NTWO,ICB,ICD
      INTEGER NKN,IGR,ISP,NSP,NGR,NGTSUM,NLTSUM,NGSUM,NLSUM,NGFSUM
      INTEGER N,JGAS,NA,IHIREAC,JAL,JRE,JPR
      INTEGER KNUMPORL,NCCOUNT,NREMAIN,NFIVE,NFOUR,NTHREE,NONE,MC
      INTEGER IR,JR,IAR,JP,JSPC

      REAL*8 RFRAC,ALFRAC,DIFF,TNUMGNA,TNUMGN
      REAL*8 TNUMLS,SUMGN,TSUMGNA,TNUMLSA

      INTEGER, SAVE :: NPLTOT,NPLFUN,NFRCOUN,NPDCOUN

      NCSP                        = NCS + ICS
      NREPT                       = 0
C
       DO 30 I                     = 1, MXGSAER 
        ISAPORL( I)                = 0
 30    CONTINUE 
C
       DO 33 I                     = 1, MAXGL 
        NEWNK(I)                   = 0
 33    CONTINUE
C
       DO 42 I                     = 1, MXGSAER
        DO 41 J                    = 1, MXGSAER
         ISPARDER(I,J)             = 0
 41     CONTINUE
 42    CONTINUE
C
       DO 100 NAR                 = 1, NALLRAT(NCS)
        NK                        = NCEQUAT(NAR,NCS) 
        IF (NK.LE.NRATES(NCS))      NALLRAT(NCSP) = NAR
        DO 60 K                   = 1, NMREAC  
         IREACT                   = IRM(K,NK,NCS)
         IF (IREACT.GT.0.AND.IREACT.LE.NSPEC(NCS)) THEN
          DO 50 L                 = 1, NPRODHI  
           IPO                    = IRM(L,NK,NCS)
           IF ((L.LE.NMREAC.OR.L.GE.NPRODLO).AND.IPO.GT.0.AND.
     1          IPO.LE.NSPEC(NCS)) ISPARDER(IPO,IREACT) = 1 
 50       CONTINUE 
         ENDIF 
 60     CONTINUE 
 100   CONTINUE
C      CONTINUE NAR = 1, NALLRAT
C
       DO 72 IREACT                = 1, NTSPEC(NCS)
        DO 70 IPO                  = 1, NTSPEC(NCS)
         IF (ISPARDER(IPO,IREACT).EQ.1) ISAPORL(IPO)=ISAPORL(IPO)+1
 70     CONTINUE 
 72    CONTINUE 
C
C *********************************************************************
C *  RE-ARRAGE SPECIES ARRAY SO THAT ALL SPECIES WITH AT LEAST ONE    * 
C *  PARTIAL DERIVATIVE TERM APPEAR FIRST, AND THOSE WITH ZERO        *
C *  APPEAR LAST.                                                     * 
C *********************************************************************
C ISCHANG = NUMBER OF ORIGINAL NSPEC SPECIES WITH AT LEAST ONE PD TERM. 
C INEWOLD = ORIGINAL SPECIES NUMBER OF EACH NEW JNEW SPECIES 
C MAPPL   = NEW SPECIES NUMBER FOR CHEMISTRY OF EACH ORIGINAL JOLD SPECIES 
C
       NOCHANG                 = NSPEC(NCS) 
       DO 110 JOLD             = 1, NTSPEC(NCS)
        IF (JOLD.GT.NSPEC(NCS)) THEN 
         MAPPL(JOLD,NCS)       = JOLD
         INEWOLD(JOLD,NCS)     = JOLD   
        ELSEIF (ISAPORL(JOLD).GT.0) THEN
         ISCHANG(NCS)          = ISCHANG(NCS) + 1
         JNEW                  = ISCHANG(NCS) 
         INEWOLD(JNEW,NCS)     = JOLD   
         MAPPL(JOLD,NCS)       = JNEW  
        ELSE
         INEWOLD(NOCHANG,NCS)  = JOLD   
         MAPPL(JOLD,NCS)       = NOCHANG  
         NOCHANG               = NOCHANG - 1
        ENDIF
 110   CONTINUE
C
C *********************************************************************
C *  RE-ARRAGE SPECIES IN ISCHANG ARRAY SO THAT SPECIES WITH THE      *
C *  FEWEST PARTIAL DERIVATIVE TERMS COMBINED ARE PLACED FIRST,       *
C *  AND THOSE WITH THE MOST APPEAR LAST. HOWEVER, SPECIES WITH ZERO  *
C *  PARTIAL DERIVATIVE TERMS STILL APPEAR AFTER ALL ISCHANG SPECIES  *
C *********************************************************************
C 
       DO 117 JNEW             = 1, ISCHANG(NCS)
        JOLD                   = INEWOLD(JNEW,NCS)
        MINVALU                = ISAPORL(JOLD)
        IMINOLD                = JOLD 
        IMINNEW                = JNEW
        DO 115 INEW            = JNEW+1, ISCHANG(NCS)
         IOLD                  = INEWOLD(INEW,NCS)
         IF (ISAPORL(IOLD).LT.MINVALU) THEN
          MINVALU              = ISAPORL(IOLD)
          IMINOLD              = IOLD  
          IMINNEW              = INEW
         ENDIF
 115    CONTINUE
        INEWOLD(IMINNEW,NCS)   = JOLD  
        INEWOLD(JNEW,NCS)      = IMINOLD  
        MAPPL(JOLD,NCS)        = IMINNEW    
        MAPPL(IMINOLD,NCS)     = JNEW    
 117   CONTINUE
C
C *********************************************************************
C *                    COUNT GROSS AND NET LOSS                       * 
C *********************************************************************
C IONER    = NUMBER OF REACTIONS WITH ONE ACTIVE REACTANT
C ITWOR    = NUMBER OF REACTIONS WITH TWO ACTIVE REACTANTS 
C ITHRR    = NUMBER OF REACTIONS WITH THREE ACTIVE REACTANTS 
C NKONER   = REACTION NUMBER OF EACH IONER REACTION 
C NKTWOR   = REACTION NUMBER OF EACH ITWOR REACTION 
C NKTHRR   = REACTION NUMBER OF EACH ITHRR REACTION 
C NUMLOST  = EVERY OCCURENCE OF A LOSS (ACTIVE & INACTIVE SPEC) 
C NUMLOSS  = EVERY NET OCCURENCE OF A LOSS WHERE THE SPECIES IS NOT 
C            REPRODUCED IN THE SAME REACTION. (ACTIVE & INACTIVE SPECIES)
C JLOSS    = REACTION NUMBER OF EACH NET LOSS OCCURRENCE
C IRM2     = IDENTIFIES EACH NEW ACTIVE SPECIES NUMBER IN EACH REACTION
C NUMKIAL  = NUMBER OF REACTIONS WITH EITHER 1, 2, OR 3 ACTIVE REACTANTS  
C NKSDT    = REACTION NUMBER OF EACH NUMKIAL REACTION 
C NRUSE    = 1,2,3 IF REACTION HAS 1, 2, OR 3 ACTIVE REACTANTS, RESPECTIVELY.
C NRREP    = 0 FOR EACH OF TWO REACTIONS WHERE THE REACTANTS ARE IDENTICAL.
C            IF MORE THAN TWO REACTIONS HAVE THE SAME REACTANTS, NRREP = 0
C            FOR THE FIRST TWO REACTIONS ONLY.
C          = 1,2,3 IF REACTION HAS 1, 2, OR 3 REACTANTS, RESPECTIVELY.
C NMOTH    = # OF OCCURRENCES WHERE INACTIVE SPEC APPEARS IN RATE EQUATION
C            EXCLUDES THIRD BODIES IN ARRAY NM3BOD (E.G., O2, N2, M, WHEN
C            THESE SPECIES DO NOT LOSE CONCENTRATION IN THE REACTION)
C NREACOTH = REACTION NUMBER OF EACH NMOTH OCCURRENCE
C LGASBINO = OLD SPECIES NUMBER OF EACH INACTIVE SPECIES  
C
       NOLOSP(NCSP)              = 0
       NKLAST                    = 0
C
       DO 230 NAR                = 1, NALLRAT(NCS)
        NK                       = NCEQUAT(NAR,NCS)
C
C *********************************************************************
C ***  DETERMINE OCCURRENCES OF INACTIVE SPECIES IN RATE EQUATIONS  ***
C *             SET ARRAY TO IDENTIFY ACTIVE LOSS SPECIES             *
C *********************************************************************
C
        IAL                      = 0
C
        DO 157 JSPC              = 1, MXGSAER
         APORL(JSPC)             = 0.d0
 157    CONTINUE
C
        DO 158 J                 = 1, NMREAC
         IREACT                  = IRM(J,NK,NCS)
         IF (IREACT.GT.0) THEN
          IRE                    = MAPPL(IREACT,NCS)
C
          APORL(IRE)             = APORL(IRE) - 1.d0
          NUMLOST(IRE,NCS)       = NUMLOST(IRE,NCS) + 1
C
          IF (IRE.LE.NSPEC(NCS)) THEN
C
           IAL                   = IAL + 1
           IRM2(IAL,NK,NCS)      = IRE
C
          ELSEIF (IRE.GT.NSPEC(NCS)) THEN
C
           IF (NK.LE.NRATES(NCS)) THEN
            NMOTH(NCS)           = NMOTH(NCS) + 1
            NMO                  = NMOTH(NCS)
            NREACOTH(NMO,NCS)    = NK
            LGASBINO(NMO,NCS)    = IREACT
           ELSE
            NOLOSP(NCS)          = NOLOSP(NCS) + 1
            NOL                  = NOLOSP(NCS)
            NKNLOSP(NOL,NCS)     = NK
            LOSINACP(NOL,NCS)    = IREACT
           ENDIF
C
          ENDIF
         ENDIF
C
 158    CONTINUE
C
C *********************************************************************
C *  SET ARRAYS TO IDENTIFY REACTIONS WITH AT LEAST ONE ACTIVE LOSS   *
C *********************************************************************
C
        IF (IAL.GT.0) THEN
         NRUSE(NK,NCS)      = IAL 
         NRREP(NK,NCS)      = IAL
C
         IF (IAL.EQ.1) THEN
          IONER(NCS)             = IONER(NCS) + 1
          NKONER(IONER(NCS),NCS) = NK 
         ELSEIF (IAL.EQ.2) THEN
          ITWOR(NCS)             = ITWOR(NCS) + 1
          NKTWOR(ITWOR(NCS),NCS) = NK 
         ELSEIF (IAL.EQ.3) THEN 
          ITHRR(NCS)             = ITHRR(NCS) + 1
          NKTHRR(ITHRR(NCS),NCS) = NK 
         ENDIF
C
C *********************************************************************
C * COMPARE TWO CONSECUTIVE REACTIONS. IF THE SPECIES (BUT NOT RATES) * 
C * ARE THE SAME, THEN SAVE MULTIPLICATIONS IN SUBFUN.F               *
C *********************************************************************
C
         IF (NKLAST.GT.0) THEN 
          IF (NRUSE(NKLAST,NCS).EQ.IAL) THEN  
           ISDIFF           = 0 
           DO 150 IB        = 1, IAL  
            JSPCL           = IRM2(IB,NKLAST,NCS) 
            JSPC            = IRM2(IB,NK    ,NCS) 
            IF (JSPCL.NE.JSPC) ISDIFF = 1 
 150       CONTINUE 
           IF (ISDIFF.EQ.0.AND.NRREP(NKLAST,NCS).NE.0) THEN 
            NRREP(NK,NCS)     = 0
            NRREP(NKLAST,NCS) = 0
            NREPT           = NREPT + 1
            ISPC1           = IRM2(1,NK,NCS) 
            ISPC2           = IRM2(2,NK,NCS) 
            ISPC3           = IRM2(3,NK,NCS)
            IF (ISPC1.GT.0) ISPC1 = INEWOLD(ISPC1,NCS)
            IF (ISPC2.GT.0) ISPC2 = INEWOLD(ISPC2,NCS)
            IF (ISPC3.GT.0) ISPC3 = INEWOLD(ISPC3,NCS)
            WRITE(IO93,155) NREPT, NK,NAMENCS(ISPC1,NCS), 
     1                   NAMENCS(ISPC2,NCS), NAMENCS(ISPC3,NCS)
 155        FORMAT('REPEAT REACTANTS: ',I5,I5,3(1X,A14))
           ENDIF 
          ENDIF 
         ENDIF 
C
C *********************************************************************
C *   DETERMINE THE NUMBER OF REACTIONS WITH ZERO ACTIVE LOSS TERMS   * 
C *********************************************************************
C NOLOSRAT = NUMBER OF ACTIVE REACTIONS WITH NO LOSS TERMS 
C NOLOSRN  = REACTION NUMBER OF EACH REACTION WITH NO LOSS TERMS
C

        ELSEIF (IAL.EQ.0) THEN
         NOLOSRAT(NCS)         = NOLOSRAT(NCS) + 1
         NOL                   = NOLOSRAT(NCS)
         NOLOSRN(NOL,NCS)      = NK
        ENDIF
C       ENDIF IAL.GT.0 
C
C *********************************************************************
C * COUNT GROSS AND NET PRODUCTION AND SET A PARTIAL DERIVATIVE ARRAY * 
C *********************************************************************
C NUMGAINT = EVERY OCCURENCE OF A PRODUCTION (ACTIVE & INACTIVE SPEC) 
C NUMGAIN  = EVERY NET OCCURENCE OF A PRODUCTION WHERE THE SPECIES IS 
C            NOT LOST IN THE SAME REACTION. (ACTIVE & INACTIVE SPEC)
C IAPROD   = NUMBER OF ACTIVE PRODUCTS IN EACH NK REACTION. USED
C            TO CALCULATE PARTIAL DERIVATIVES IN PDERIV.F. 
C IRM2     = NEW SPECIES # OF EACH ACTIVE PRODUCT IN EACH NK REACTION
C
        IAP                       = NPRODLO - 1
        DO 210 K                  = NPRODLO, NPRODHI  
         IPROD                    = IRM(K,NK,NCS)
         IF (IPROD.GT.0) THEN
          IPR                     = MAPPL(IPROD,NCS)
          RFRAC                   = FKOEF(K,NK,NCS)
          LFRAC                   = INT(RFRAC + SMAL1) 
          ALFRAC                  = FLOAT(LFRAC)
          DIFF                    = ABS(RFRAC-ALFRAC)
C
C ******************** PRODUCTION TERM IS A FRACTION ******************
C
          IF (DIFF.GT.SMAL1) THEN 
           IF (IPR.LE.NSPEC(NCS)) THEN 
            NGNFRAC(NCS)          = NGNFRAC(NCS) + 1 
            NGN                   = NGNFRAC(NCS) 
            IGNFRAC( NGN,NCS)     = IPR 
            NKGNFRAC(NGN,NCS)     = NK 
            FRACP(   NGN,NCS)     = RFRAC  
           ENDIF 
           KPRODS                 = 1
           NUMGFRT( IPR,NCS)      = NUMGFRT( IPR,NCS) + 1
           FRACGAIN(IPR,NCS)      = FRACGAIN(IPR,NCS) + RFRAC 
C
C ******************* PRODUCTION TERM IS NON-FRACTION *****************
C
          ELSE
           APORL(IPR)             = APORL(IPR) + RFRAC
           KPRODS                 = LFRAC
           NUMGAINT(IPR,NCS)      = NUMGAINT(IPR,NCS) + LFRAC
           FKOEF(K,NK,NCS)        = 1.d0
          ENDIF  
C
C ******************* IDENTIFY ALL PRODUCTION TERMS *******************
C
          IF (IPR.LE.NSPEC(NCS)) THEN
           DO 170 L               = 1, KPRODS
            IAP                   = IAP + 1
            IAPROD(NK,NCS)        = IAP
            IRM2(IAP,NK,NCS)      = IPR
            FK2( IAP,NK,NCS)      = FKOEF(K,NK,NCS)
 170       CONTINUE  
          ENDIF
C
         ENDIF
C
 210    CONTINUE
C
C *********************************************************************
C *  FIND NET PROD AND LOSS TERMS FOR ALL BUT FRACTIONATED PRODUCTS   * 
C *********************************************************************
C
         DO 220 JSPC              = 1, NTSPEC(NCS)
          IF (ABS(APORL(JSPC)).LT.SMAL1) THEN
           KDIF                   = 0 
C
          ELSEIF (APORL(JSPC).GT.0.) THEN 
           KDIF                   = INT(APORL(JSPC) + 0.00001)
           DO 190 L               = 1, KDIF 
            NUMGAIN(JSPC,NCS)     = NUMGAIN(JSPC,NCS) + 1
            NUMPORL(JSPC,NCS)     = NUMPORL(JSPC,NCS) + 1
            NPL                   = NUMPORL(JSPC,NCS)
            JPORL(JSPC,NPL,NCS)   = NK + NTRATES(NCS)  
 190       CONTINUE  
          ELSE 
           KDIF                   = -INT(APORL(JSPC) - 0.00001)
           DO 140 L               = 1, KDIF  
            NUMLOSS(JSPC,NCS)     = NUMLOSS(JSPC,NCS) + 1
            NUMPORL(JSPC,NCS)     = NUMPORL(JSPC,NCS) + 1
            NPL                   = NUMPORL(JSPC,NCS)
            JPORL(JSPC,NPL,NCS)   = NK 
 140       CONTINUE
          ENDIF 
C
          IF (NK.LE.NRATES(NCS)) THEN 
           NUMLOSS(JSPC,NCSP)     = NUMLOSS(JSPC,NCS)
           NUMGAIN(JSPC,NCSP)     = NUMGAIN(JSPC,NCS)
           NUMPORL(JSPC,NCSP)     = NUMPORL(JSPC,NCS)
          ENDIF 
C
 220     CONTINUE 
C
         IF (NK.LE.NRATES(NCS)) THEN 
          NOLOSRAT(NCSP)          = NOLOSRAT(NCS) 
          NGNFRAC( NCSP)          = NGNFRAC( NCS)
          IONER(   NCSP)          = IONER(   NCS)
         ENDIF 
C
         NKLAST                   = NK 
C
 230   CONTINUE
C      CONTINUE N = 1, NTRATES
C
C *********************************************************************
C * SET ARRAY FOR REORDERING RATES FROM 3..2..1..0 BODY REACTIONS     *
C *********************************************************************
C INOREP   = LAST REORDERED REACTION NUMBER PRIOR TO SETS OF TWO
C            REACTIONS WITH TWO REACTANTS  
C NOLDFNEW = OLD REACTION RATE # CORRESP. TO EACH REORDERED REACTION
C NEWFOLD  = NEW REACTION RATE # CORRESP. TO EACH ORIGINAL RATE NUMBER
C
       IC                 = 0
       DO 235 I           = 1, ITHRR(NCS)   
        IC                = IC + 1
        NK                = NKTHRR(I,NCS)
        NK1               = NK + NTRATES(NCS) 
        NOLDFNEW(IC, NCS) = NK
        NEWFOLD( NK, NCS) = IC
        NEWFOLD( NK1,NCS) = IC + NALLRAT(NCS) 
 235   CONTINUE 
C
       NTWO                = ITHRR(NCS) + ITWOR(NCS) 
       ICB                 = NTWO + 1 
       DO 237 I            = 1, ITWOR(NCS)   
        NK                 = NKTWOR(I,NCS)
        NK1                = NK + NTRATES(NCS) 
        IF (NRREP(NK,NCS).GT.0) THEN  
         IC                = IC + 1
         ICD               = IC
        ELSE 
         ICB               = ICB - 1
         ICD               = ICB
        ENDIF 
        NOLDFNEW(ICD, NCS) = NK
        NEWFOLD( NK,  NCS) = ICD 
        NEWFOLD( NK1, NCS) = ICD + NALLRAT(NCS) 
 237   CONTINUE 
C
       INOREP(NCS)         = IC 
       IC                  = NTWO 
       DO 239 I            = 1, IONER(NCS)   
        IC                 = IC + 1
        NK                 = NKONER(I,NCS)
        NK1                = NK + NTRATES(NCS) 
        NOLDFNEW(IC, NCS)  = NK
        NEWFOLD( NK, NCS)  = IC
        NEWFOLD( NK1,NCS)  = IC + NALLRAT(NCS) 
 239   CONTINUE 
C
       DO 241 I            = 1, NOLOSRAT(NCS)
        IC                 = IC + 1
        NK                 = NOLOSRN(I,NCS)
        NK1                = NK + NTRATES(NCS)
        NOLDFNEW(IC, NCS)  = NK
        NEWFOLD( NK, NCS)  = IC
        NEWFOLD( NK1,NCS)  = IC + NALLRAT(NCS)
 241   CONTINUE
C
       IF (IC.NE.NALLRAT(NCS)) THEN
        WRITE(6,245) IC, NALLRAT(NCS)
        CALL GEOS_CHEM_STOP
       ENDIF
C
C *********************************************************************
C                SET A SLIGHTLY MORE EFFICIENT PHOTO ARRAY 
C *********************************************************************
C
       DO 243 J          = 1, JPHOTRAT(NCS)
        NK               = NKPHOTRAT(J,NCS)
        NKN              = NEWFOLD(NK,NCS)
        NKNPHOTRT(J,NCS) = NKN
 243   CONTINUE
C
 245   FORMAT('JSPARSE: IC NE NALLRAT =',2(I5))
C
C *********************************************************************
C ****** DETERMINE NUMBER OF SPECIES WITH GROSS/NET LOSSES/GAINS ******
C *********************************************************************
C NSPCSOLV = # OF ACTIVE SPECIES WITH AT LEAST ONE GROSS LOSS
C ISOLVSPC = SPECIES NUMBER OF EACH NSPCSOLV SPECIES
C ISGAINR  = # OF ACTIVE SPECIES WITH AT LEAST ONE NET CHEM GAIN 
C IGAINR   = SPECIES NUMBER OF EACH ISGAINR SPECIES
C ISGAINE  = # OF ACTIVE SPECIES WITH AT LEAST 1 NET CHEM GAIN 
C IGAINR   = SPECIES NUMBER OF EACH ISGAINR SPECIES
C NOGAINE  = # OF ACTIVE SPECIES WITH ZERO NET CHEM OR GAINS 
C NGAINE   = SPECIES NUMBER OF EACH NOGAINE SPECIES
C ISPORL   = # OF ACTIVE SPECIES WITH AT LEAST ONE NET PRODUCTION
C            OR LOSS TERM FOR SMVGEAR.
C IPORL    = SPECIES NUMBER OF EACH ISPORL SPECIES
C
       DO 300 JOLD             = 1, NSPEC(NCS) 
        JNEW                   = MAPPL(JOLD,NCS)
C
        IF (NUMGAIN(JNEW,NCS).GT.0) THEN
         ISGAINR(NCS)          = ISGAINR(NCS) + 1
         IGR                   = ISGAINR(NCS)
         IGAINR(IGR,NCS)       = JNEW  
        ENDIF 
C
        IF (NUMPORL(JNEW,NCS).GT.0) THEN
         ISPORL(NCS)           = ISPORL(NCS) + 1 
         ISP                   = ISPORL(NCS)
         IPORL(ISP,NCS)        = JNEW 
        ENDIF
C
        IF (NUMLOST(JNEW,NCS).GT.0) THEN
         NSPCSOLV(NCS)         = NSPCSOLV(NCS) + 1
         NSP                   = NSPCSOLV(NCS)
         ISOLVSPC(NSP,NCS)     = JNEW  
        ENDIF
C
        IF (NUMGAIN(JNEW,NCS).GT.0.OR.FRACGAIN(JNEW,NCS).GT.0) THEN
         ISGAINE(NCS)         = ISGAINE(NCS) + 1
         IGR                  = ISGAINE(NCS)
         IGAINE(IGR,NCS)      = JNEW  
        ELSEIF (NUMLOSS(JNEW,NCS).GT.0) THEN 
         NOGAINE(NCS)         = NOGAINE(NCS) + 1
         NGR                  = NOGAINE(NCS)
         NGAINE(NGR,NCS)      = JNEW    
        ENDIF
C
 300   CONTINUE
C
C *********************************************************************
C ********  CHECK DIMENSIONS RESULTING FROM GAINS AND LOSSES  *********
C *********************************************************************
C
       NGTSUM   = 0
       NLTSUM   = 0
       NGSUM    = 0
       NLSUM    = 0
       NGFSUM   = 0
       DO 260 K = 1, NTSPEC(NCS)
        J       = INEWOLD(K,NCS) 
        NGTSUM  = NGTSUM + NUMGAINT(K,NCS) 
        NLTSUM  = NLTSUM + NUMLOST( K,NCS) 
        NGSUM   = NGSUM  + NUMGAIN( K,NCS) 
        NLSUM   = NLSUM  + NUMLOSS( K,NCS) 
        NGFSUM  = NGFSUM + NUMGFRT( K,NCS) 
        IF (NUMGAINT(K,NCS)   .GT.   MAXGL .OR.
     1      NUMLOST( K,NCS)   .GT.   MAXGL) THEN 
         WRITE(6,280) NAMENCS(J,NCS), NUMGAINT(K,NCS), NUMLOST(K,NCS)
         CALL GEOS_CHEM_STOP
        ENDIF
 260   CONTINUE
C
       IF (IOREAC.EQ.1) THEN
        WRITE(IO93,*)
        WRITE(IO93,240)
        DO 270 K = 1, NTSPEC(NCS) 
         J       = INEWOLD(K,NCS) 
         WRITE(IO93,250)NAMENCS( J,NCS),NUMGAINT(K,NCS),NUMGAIN( K,NCS),
     1                NUMLOST( K,NCS),NUMLOSS( K,NCS),NUMGAINT(K,NCS) 
     2               -NUMLOST( K,NCS)-NUMGAIN( K,NCS)+NUMLOSS( K,NCS),
     3                FRACGAIN(K,NCS),NUMGFRT( K,NCS)
 270    CONTINUE
        WRITE(IO93,250) 'OVERALL       ',NGTSUM, NGSUM, NLTSUM, NLSUM,
     1               NGTSUM - NLTSUM - NGSUM + NLSUM, 0., NGFSUM
       ENDIF
C
       IF (NMOTH(  NCS).GT.MAXGL2.OR.NOLOSP(NCS).GT.MAXGL3.OR.
     1     NGNFRAC(NCS).GT.MAXGL) THEN
        WRITE(6,275) MAXGL2, NMOTH(  NCS), MAXGL3, NOLOSP(NCS),
     1               MAXGL,  NGNFRAC(NCS)
        CALL GEOS_CHEM_STOP
       ENDIF
C
C *********************************************************************
C *       CHECK WHETHER CHEMICAL SYSTEM IS ATOM-CONSERVATIVE          *
C *********************************************************************
C JMBCOMP = SPECIES NUMBER FOR EACH SPECIES IN A MASS BAL. GROUP
C MBCOMP  = COUNTS THE NUMBER OF MASS BALANCE SPECIES IN EACH M.B GROUP
C NMASBAL = NUMBER OF MASS BALANCE GROUPS (E.G. S, N, C ARE GROUPS)
C WTMB(1) = NUMBER OF ATOMS OF A GIVEN MASS BALANCE SPECIES PER MOLECULE 
C
       WRITE(IO93,360) CHEMTYP(NCS)
C
       IF (NCS.LE.NCSGAS) THEN
C 
C ----------------------------   GAS-PHASE   -------------------------- 
C
        DO 385 N     = 1, NMASBAL 
         IF (MBCOMP(N,MB1).GT.0) THEN 
          TNUMGN     = 0
          TNUMLS     = 0
          WRITE(IO93,325) NAMEMB(N)
          DO 380 J   = 1, MBCOMP(N,MB1) 
           JGAS      = JMBCOMP(N,J,MB1)
           JNEW      = MAPPL(JGAS,NCS)
           SUMGN     = NUMGAIN(JNEW,NCS) + FRACGAIN(JNEW,NCS) 
           TNUMGNA   = SUMGN             * WTMB(N,JGAS,MB1)  
           TNUMLSA   = NUMLOSS(JNEW,NCS) * WTMB(N,JGAS,MB1) 
           TNUMGN    = TNUMGN + TNUMGNA
           TNUMLS    = TNUMLS + TNUMLSA
           WRITE(IO93,320) NAMEGAS(JGAS), TNUMGNA, TNUMLSA, 0 
 380      CONTINUE
          WRITE(IO93,370) TNUMGN, TNUMLS, TNUMGN - TNUMLS 
         ENDIF 
 385    CONTINUE
       ENDIF
C
       WRITE(IO93,375) NALLRAT(NCSP), NALLRAT(NCS) - NALLRAT(NCSP),
     1                 NALLRAT(NCS) 
C
 360   FORMAT(/'CHANGE IN MOLES DUE TO ',A14,' CHEMISTRY')
 325   FORMAT('MASS BALANCE GROUP              = ',A14)
 320   FORMAT('GAINS/LOSSES FOR ',A14,' = ',2(F8.3),I5)     
 370   FORMAT('TOTAL GAINS - LOSSES            = ',3(F8.3)) 
 375   FORMAT(/'# KINETIC REACTIONS: ',I5,' PHOTORATES: ',I5,
     1        ' TOTAL: ',I5) 
 240   FORMAT('SPEC           NUMGT  NUMG  NUMLT NUML   NGT-NLT-', 
     1        'NG+NL FRACGN NUMGFT') 
 250   FORMAT(A14,4(2X,I4),7X,I4,3X,F8.3,I5)
 280   FORMAT('GEARSET: SPEC ',A6,' DIMENS EXCEEDED. EITHER NUMGAINT ', 
     1        'NUMLOSS,NUMGAIN, OR NUMLOST > MAXGL ',
     2        4(I3,1X)) 
 275   FORMAT('JSPARSE: ONE OF THE DIMENSIONS BELOW IS TOO SMALL:',/,
     1        'DIMENSION: MAXGL2   =  ',I4,' VARIABLE: NMOTH    = ',I4/  
     2        'DIMENSION: MAXGL3   =  ',I4,' VARIABLE: NOLOSP   = ',I4/
     3        'DIMENSION: MAXGL    =  ',I4,' VARIABLE: NGNFRAC  = ',I4)  
C
C *********************************************************************
C *********************************************************************
C **        SET ARRAYS TO TAKE ADVANTAGE OF SPARSE MATRICES          ** 
C *********************************************************************
C *********************************************************************
C
C IFSUN  = 1 THEN DAY-CHEMISTRY;  = 2 THEN NIGHT CHEMISTRY
C NCSP   = NCS       FOR DAYTIME   TROP-GAS, STRAT-GAS CHEM  
C NCSP   = NCS + ICP FOR NIGHTTIME TROP-GAS, STRAT-GAS CHEM  
C
C LZERO    = 1 IF AN ARRAY SPOT IS FILLED WITH A NON-ZERO VALUE. LZERO
C            IS UPDATED AS WE SIMULATE THE ORDER OF CALCULATIONS DURING
C            A PRACTICE L-U DECOMPOSITION
C MXGSAER  = LARGER OF IGAS, IAERTY
C
C
      IF (IFNONE.EQ.0) THEN
       IFNONE                 = 1
       NPLFUN                 = 0 
       NFRCOUN                = 0 
       NPDCOUN                = 0 
       NPLTOT                 = 0 
      ENDIF
C
      DO 700 IFSUN            = 1, 2 
       NCSP                   = (IFSUN - 1) * ICS + NCS
C
       DO 517 I               = 1, MXGSAER
        DO 515 J              = 1, MXGSAER
         LZERO(J,I)           = 0
 515    CONTINUE
        LZERO(I,I)            = 1
 517   CONTINUE 
C
       DO 504 NA              = 1, NALLRAT(NCSP)
        NK                    = NCEQUAT(NA,NCS)
        IHIREAC               = NRUSE(  NK,NCS)
        DO 502 IAL            = 1, IHIREAC
         IRE                  = IRM2(IAL,NK,NCS)
         DO 490 JAL           = 1, IHIREAC 
          JRE                 = IRM2(JAL,NK,NCS)
          LZERO(JRE,IRE)      = 1
 490     CONTINUE
         DO 500 IAP           = NPRODLO, IAPROD(NK,NCS)
          JPR                 = IRM2(IAP,NK,NCS)
          LZERO(JPR,IRE)      = 1
 500     CONTINUE
 502    CONTINUE
 504   CONTINUE
C
C *********************************************************************
C *   SET DECOMPOSITION AND BACK-SUBSTITUTION SPARSE-MATRIX ARRAYS    *
C *********************************************************************
C
       CALL KSPARSE

C
C *********************************************************************
C *    SET ARRAYS TO IMPROVE EFFICIENCY OF FIRST-DERIVATIVE CALCS     * 
C *********************************************************************
C *********************************************************************
C **   SET ARRAYS FOR KINETIC AND PHOTO PRODUCTION AND LOSS RATES    **
C *********************************************************************
C
       NPLLO(NCSP)         = NPLTOT + 1
       DO 670 I            = 1, ISPORL(NCS)
        JSPC               = IPORL(I,NCS)
        KNUMPORL           = NUMPORL(JSPC,NCSP) 
        NCCOUNT            = 0 
        NPLTOT             = NPLTOT + 1
        NREMAIN            = KNUMPORL
        NFIVE              = (NREMAIN + 0.0001) / 5 
        NREMAIN            =  NREMAIN - NFIVE   * 5 
        NFOUR              = (NREMAIN + 0.0001) / 4 
        NREMAIN            =  NREMAIN - NFOUR   * 4
        NTHREE             = (NREMAIN + 0.0001) / 3  
        NREMAIN            =  NREMAIN - NTHREE  * 3 
        NTWO               = (NREMAIN + 0.0001) / 2   
        NREMAIN            =  NREMAIN - NTWO    * 2  
        NONE               = (NREMAIN + 0.0001)  
        NREMAIN            =  NREMAIN - NONE
C
        JSPNPL(NPLTOT)     = JSPC 
        NPL5(  NPLTOT)     = NPLFUN       + 1
        NPH5(  NPLTOT)     = NPLFUN       + NFIVE  
        NPL4(  NPLTOT)     = NPH5(NPLTOT) + 1
        NPH4(  NPLTOT)     = NPH5(NPLTOT) + NFOUR   
        NPL3(  NPLTOT)     = NPH4(NPLTOT) + 1
        NPH3(  NPLTOT)     = NPH4(NPLTOT) + NTHREE
        NPL2(  NPLTOT)     = NPH3(NPLTOT) + 1
        NPH2(  NPLTOT)     = NPH3(NPLTOT) + NTWO
        NPL1(  NPLTOT)     = NPH2(NPLTOT) + 1
        NPH1(  NPLTOT)     = NPH2(NPLTOT) + NONE
        NPLFUN             = NPH1(NPLTOT)
C
        DO 649 N           = 1, KNUMPORL 
         NK                = JPORL(JSPC,N,NCS) 
         NEWNK(N)          = NEWFOLD(NK,NCS)  
 649    CONTINUE 
C
        DO 651 MC          = NPL5(NPLTOT), NPH5(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         LOSSRB(MC)        = NEWNK(NCCOUNT+2) 
         LOSSRC(MC)        = NEWNK(NCCOUNT+3) 
         LOSSRD(MC)        = NEWNK(NCCOUNT+4) 
         LOSSRE(MC)        = NEWNK(NCCOUNT+5) 
         NCCOUNT           = NCCOUNT + 5
 651    CONTINUE 
C
        DO 652 MC          = NPL4(NPLTOT), NPH4(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         LOSSRB(MC)        = NEWNK(NCCOUNT+2) 
         LOSSRC(MC)        = NEWNK(NCCOUNT+3) 
         LOSSRD(MC)        = NEWNK(NCCOUNT+4) 
         NCCOUNT           = NCCOUNT + 4  
 652    CONTINUE 
C
        DO 653 MC          = NPL3(NPLTOT), NPH3(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         LOSSRB(MC)        = NEWNK(NCCOUNT+2) 
         LOSSRC(MC)        = NEWNK(NCCOUNT+3) 
         NCCOUNT           = NCCOUNT + 3   
 653    CONTINUE 
C
        DO 654 MC          = NPL2(NPLTOT), NPH2(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         LOSSRB(MC)        = NEWNK(NCCOUNT+2) 
         NCCOUNT           = NCCOUNT + 2    
 654    CONTINUE 
C
        DO 656 MC          = NPL1(NPLTOT), NPH1(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         NCCOUNT           = NCCOUNT + 1     
 656    CONTINUE 
C
 670   CONTINUE 
       NPLHI(NCSP)         = NPLTOT
C
C *********************************************************************
C *              SET ARRAY FOR FRACTIONATED PRODUCTS                  *  
C *********************************************************************
C
       NFRLO(NCSP)          = NFRCOUN + 1 
       DO 695 I             = 1, NGNFRAC(NCSP)
        JSPC                = IGNFRAC(I,NCS)
        NFRCOUN             = NFRCOUN + 1 
        JSPCNFR(NFRCOUN)    = JSPC 
        NK                  = NKGNFRAC(I,NCS)  
        NKNFR(  NFRCOUN)    = NEWFOLD(NK,NCS)
        FRACNFR(NFRCOUN)    = FRACP(I,NCS)
 695   CONTINUE 
       NFRHI(NCSP)          = NFRCOUN
C
C *********************************************************************
C * SET ARRAYS TO IMPROVE EFFICIENCY OF PARTIAL DERIVATIVE CALCS      * 
C *********************************************************************
C
       NPDLO(NCSP)           = NPDCOUN + 1
C
       DO 974 NA             = 1, NALLRAT(NCSP) 
        NK                   = NCEQUAT(NA,NCS) 
        IHIREAC              = NRUSE(  NK,NCS) 
C
        DO 972 IAL           = 1, IHIREAC
         IR                  = IRM2(IAL,NK,NCS)
         DO 960 JAL          = 1, IHIREAC 
          JR                 = IRM2(JAL,NK,NCS)
          IAR                = JARRAYPT(JR,IR)
          NPDCOUN            = NPDCOUN + 1 
          NKPDTERM(NPDCOUN)  = NEWFOLD(NK,NCS)  
          IPOSPD(  NPDCOUN)  = IAR 
          IIALPD(  NPDCOUN)  = IAL  
          FRACPL(  NPDCOUN)  = -1.
 960     CONTINUE
C
         DO 970 IAP          = NPRODLO, IAPROD(NK,NCS)
          JP                 = IRM2(IAP,NK,NCS)
          IAR                = JARRAYPT(JP,IR)
          NPDCOUN            = NPDCOUN + 1 
          NKPDTERM(NPDCOUN)  = NEWFOLD(NK,NCS)  
          IPOSPD(  NPDCOUN)  = IAR 
          IIALPD(  NPDCOUN)  = IAL  
          FRACPL(  NPDCOUN)  = FK2(IAP,NK,NCS)  
 970     CONTINUE
 972    CONTINUE
 974   CONTINUE
C
       NPDHI(NCSP)          = NPDCOUN
C
C *********************************************************************
C **        CHECK DIMENSIONS AND PRINT OUT ARRAY SAVINGS             ** 
C *********************************************************************
C
       IF (NPLTOT   .GT. MXCOUNT4  .OR. NPLFUN   .GT. MXCOUNT4 .OR.
     3     NFRCOUN  .GT. MXCOUNT4 .OR.  NPDCOUN  .GT. MXCOUNT2) THEN
        WRITE(6,645) MXCOUNT4, NPLTOT,    MXCOUNT4, NPLFUN,
     2               MXCOUNT4, NFRCOUN,   MXCOUNT2, NPDCOUN
        CALL GEOS_CHEM_STOP
       ENDIF
C
 700  CONTINUE
C     CONTINUE IFSUN = 1, 2
C
 645  FORMAT('ONE OF THE DIMENSIONS BELOW IS TOO SMALL:',/,
     1       'DIMENSION: MXCOUNT4 =  ',I5,' VARIABLE: NPLTOT   = ',I5,/,
     2       'DIMENSION: MXCOUNT4 =  ',I5,' VARIABLE: NPLFUN   = ',I5,/,
     3       'DIMENSION: MXCOUNT4 =  ',I5,' VARIABLE: NFRCOUN  = ',I5,/,
     4       'DIMENSION: MXCOUNT2 =  ',I5,' VARIABLE: NPDCOUN  = ',I5)
C
C *********************************************************************
C ********************** END OF SUBROUTINE JSPARSE ********************
C *********************************************************************
C
      RETURN                                                             
      END SUBROUTINE JSPARSE
