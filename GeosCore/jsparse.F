! $Id: jsparse.f,v 1.1 2009/09/16 14:06:23 bmy Exp $
      SUBROUTINE JSPARSE( am_I_Root )
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
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument even when
!                              called from the traditional driver "main.F"
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD, ONLY : T3
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP

      USE CMN_SIZE_MOD  ! Size parameters
      USE COMODE_LOOP_MOD  ! SMVGEAR II arrays

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! *********************************************************************
! ************        WRITTEN BY MARK JACOBSON (1993)      ************
! ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
! ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!                                                                         
!        JJ  SSSSSSS  PPPPPPP     A      RRRRRRR  SSSSSSS  EEEEEEE
!         J  S        P     P    A A     R     R  S        E 
!         J  SSSSSSS  PPPPPPP   A   A    RRRRRRR  SSSSSSS  EEEEEEE
!   J     J        S  P        AAAAAAA   R  R           S  E
!   JJJJJJJ  SSSSSSS  P       A       A  R    R   SSSSSSS  EEEEEEE  
!
! *********************************************************************
! * THIS ROUTINE SETS UP SPARSE-MATRIX AND OTHER ARRAYS FOR SMVGEAR   *
! * (SPARSE-MATRIX VECTORIZED GEAR-CODE. IT SETS ARRAYS FOR GAS-      *
! * -PHASE, AQUEOUS-PHASE, AND ANY OTHER TYPE OF CHEMISTRY. IT ALSO   *
! * SETS ARRAYS FOR BOTH DAY AND NIGHT CHEMISTRY OF EACH TYPE.        *
! *                                                                   *
! * HOW TO CALL SUBROUTINE:                                           *
! * ----------------------                                            *
! *  CALL JSPARSE.F FROM READCHEM.F WITH                              * 
! *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
! *********************************************************************
!
! *********************************************************************
! ******* SETS UP ARRAYS FOR GAS- AND AQUEOUS-PHASE CHEMISTRY  ******** 
! * INCLUDES ARRAYS FOR CALCULATING FIRST DERIVATIVES, PARTIAL DERIV- *
! * ATIVES, MATRIX DECOMPOSTION, AND MATRIX BACK-SUBSTITUTION. FIRST, *
! * JSPARSE RE-ORDERS THE ORDINARY DIFFERENTIAL EQUATIONS TO MAXIMIZE *
! * THE NUMBER OF ZEROS IN THE MATRIX OF PARTIAL DERIVATIVES. IT      *
! * LATER SETS ARRAYS TO ELIMINATE ALL CALCULATIONS INVOLVING A ZERO. * 
! *********************************************************************
* 
! NTSPEC    = TOTAL NUMBER OF ACTIVE + INACTIVE SPECIES.
! NSPEC     = TOTAL NUMBER OF ACTIVE SPECIES.
! NMREAC    = 3 = MAXIMUM NUMBER OF ACTIVE REACTANTS IN A REACTION 
! NALLREAC  = 4 = TOTAL REACTANT POSITIONS IN A REACTION 
! NMPROD    = 5 = MAXIMUN NUMBER OF ACTIVE PRODUCTS IN A REACTION 
! NPRODLO   = NALLREAC  + 1 = LOWEST PRODUCT POSITION NUMBER. 
! NPRODHI   = NALLREAC + NMPROD = HIGHEST PRODUCT POSITION NUMBER. 
!
! *********************************************************************
! * DETERMINE HOW MANY PARTIAL DERIV TERMS ARE NEEDED FOR EACH SPECIES*
! *********************************************************************
! IFREPRO   = 1 THEN SPECIES IS LOST AND REPRODUCED IN REACTION NK 
! IRM       = SPECIES # OF EACH REACT OR PRODUCT IN EACH NK REACTION
! ISAPORL   = COUNTS PARTIAL DERIVATIVE TERMS FOR EACH SPECIES
! FKOEF     = 1, 2, FRACTION, OR MORE = # OF A GIVEN REACTANT OR PRODUCTS
!             E.G. REACTION      A + B  --> 2C + 0.34D  + D 
!             VALUE OF FKOEF     1   1      2    0.34     1     
! NCS       = 1..NCSGAS FOR GAS CHEMISTRY                            
! NCSP      = NCS        FOR DAYTIME   GAS CHEM 
!           = NCS   +ICS FOR NIGHTTIME GAS CHEM            
! NK        = REACTION # OF EACH REACTION 
! NRATES    = NUMBER OF KINETIC (NON-PHOTO) RATE COEFFICIENTS
! NTRATES   = NUMBER OF KINETIC PLUS PHOTO  RATE COEFFICIENTS
! NALLRAT   = NUMBER OF KINETIC PLUS PHOTO REACTION RATES  
!

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
!
       DO 30 I                     = 1, MXGSAER 
        ISAPORL( I)                = 0
 30    CONTINUE 
!
       DO 33 I                     = 1, MAXGL 
        NEWNK(I)                   = 0
 33    CONTINUE
!
       DO 42 I                     = 1, MXGSAER
        DO 41 J                    = 1, MXGSAER
         ISPARDER(I,J)             = 0
 41     CONTINUE
 42    CONTINUE
!
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
!      CONTINUE NAR = 1, NALLRAT
!
       DO 72 IREACT                = 1, NTSPEC(NCS)
        DO 70 IPO                  = 1, NTSPEC(NCS)
         IF (ISPARDER(IPO,IREACT).EQ.1) ISAPORL(IPO)=ISAPORL(IPO)+1
 70     CONTINUE 
 72    CONTINUE 
!
! *********************************************************************
! *  RE-ARRAGE SPECIES ARRAY SO THAT ALL SPECIES WITH AT LEAST ONE    * 
! *  PARTIAL DERIVATIVE TERM APPEAR FIRST, AND THOSE WITH ZERO        *
! *  APPEAR LAST.                                                     * 
! *********************************************************************
! ISCHANG = NUMBER OF ORIGINAL NSPEC SPECIES WITH AT LEAST ONE PD TERM. 
! INEWOLD = ORIGINAL SPECIES NUMBER OF EACH NEW JNEW SPECIES 
! MAPPL   = NEW SPECIES NUMBER FOR CHEMISTRY OF EACH ORIGINAL JOLD SPECIES 
!
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
!
! *********************************************************************
! *  RE-ARRAGE SPECIES IN ISCHANG ARRAY SO THAT SPECIES WITH THE      *
! *  FEWEST PARTIAL DERIVATIVE TERMS COMBINED ARE PLACED FIRST,       *
! *  AND THOSE WITH THE MOST APPEAR LAST. HOWEVER, SPECIES WITH ZERO  *
! *  PARTIAL DERIVATIVE TERMS STILL APPEAR AFTER ALL ISCHANG SPECIES  *
! *********************************************************************
! 
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
!
! *********************************************************************
! *                    COUNT GROSS AND NET LOSS                       * 
! *********************************************************************
! IONER    = NUMBER OF REACTIONS WITH ONE ACTIVE REACTANT
! ITWOR    = NUMBER OF REACTIONS WITH TWO ACTIVE REACTANTS 
! ITHRR    = NUMBER OF REACTIONS WITH THREE ACTIVE REACTANTS 
! NKONER   = REACTION NUMBER OF EACH IONER REACTION 
! NKTWOR   = REACTION NUMBER OF EACH ITWOR REACTION 
! NKTHRR   = REACTION NUMBER OF EACH ITHRR REACTION 
! NUMLOST  = EVERY OCCURENCE OF A LOSS (ACTIVE & INACTIVE SPEC) 
! NUMLOSS  = EVERY NET OCCURENCE OF A LOSS WHERE THE SPECIES IS NOT 
!            REPRODUCED IN THE SAME REACTION. (ACTIVE & INACTIVE SPECIES)
! JLOSS    = REACTION NUMBER OF EACH NET LOSS OCCURRENCE
! IRM2     = IDENTIFIES EACH NEW ACTIVE SPECIES NUMBER IN EACH REACTION
! NUMKIAL  = NUMBER OF REACTIONS WITH EITHER 1, 2, OR 3 ACTIVE REACTANTS  
! NKSDT    = REACTION NUMBER OF EACH NUMKIAL REACTION 
! NRUSE    = 1,2,3 IF REACTION HAS 1, 2, OR 3 ACTIVE REACTANTS, RESPECTIVELY.
! NRREP    = 0 FOR EACH OF TWO REACTIONS WHERE THE REACTANTS ARE IDENTICAL.
!            IF MORE THAN TWO REACTIONS HAVE THE SAME REACTANTS, NRREP = 0
!            FOR THE FIRST TWO REACTIONS ONLY.
!          = 1,2,3 IF REACTION HAS 1, 2, OR 3 REACTANTS, RESPECTIVELY.
! NMOTH    = # OF OCCURRENCES WHERE INACTIVE SPEC APPEARS IN RATE EQUATION
!            EXCLUDES THIRD BODIES IN ARRAY NM3BOD (E.G., O2, N2, M, WHEN
!            THESE SPECIES DO NOT LOSE CONCENTRATION IN THE REACTION)
! NREACOTH = REACTION NUMBER OF EACH NMOTH OCCURRENCE
! LGASBINO = OLD SPECIES NUMBER OF EACH INACTIVE SPECIES  
!
       NOLOSP(NCSP)              = 0
       NKLAST                    = 0
!
       DO 230 NAR                = 1, NALLRAT(NCS)
        NK                       = NCEQUAT(NAR,NCS)
!
! *********************************************************************
! ***  DETERMINE OCCURRENCES OF INACTIVE SPECIES IN RATE EQUATIONS  ***
! *             SET ARRAY TO IDENTIFY ACTIVE LOSS SPECIES             *
! *********************************************************************
!
        IAL                      = 0
!
        DO 157 JSPC              = 1, MXGSAER
         APORL(JSPC)             = 0.d0
 157    CONTINUE
!
        DO 158 J                 = 1, NMREAC
         IREACT                  = IRM(J,NK,NCS)
         IF (IREACT.GT.0) THEN
          IRE                    = MAPPL(IREACT,NCS)
!
          APORL(IRE)             = APORL(IRE) - 1.d0
          NUMLOST(IRE,NCS)       = NUMLOST(IRE,NCS) + 1
!
          IF (IRE.LE.NSPEC(NCS)) THEN
!
           IAL                   = IAL + 1
           IRM2(IAL,NK,NCS)      = IRE
!
          ELSEIF (IRE.GT.NSPEC(NCS)) THEN
!
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
!
          ENDIF
         ENDIF
!
 158    CONTINUE
!
! *********************************************************************
! *  SET ARRAYS TO IDENTIFY REACTIONS WITH AT LEAST ONE ACTIVE LOSS   *
! *********************************************************************
!
        IF (IAL.GT.0) THEN
         NRUSE(NK,NCS)      = IAL 
         NRREP(NK,NCS)      = IAL
!
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
!
! *********************************************************************
! * COMPARE TWO CONSECUTIVE REACTIONS. IF THE SPECIES (BUT NOT RATES) * 
! * ARE THE SAME, THEN SAVE MULTIPLICATIONS IN SUBFUN.F               *
! *********************************************************************
!
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
            IF ( am_I_Root ) THEN
               WRITE(IO93,155) NREPT, NK,NAMENCS(ISPC1,NCS), 
     &            NAMENCS(ISPC2,NCS), NAMENCS(ISPC3,NCS)
            ENDIF
 155        FORMAT('REPEAT REACTANTS: ',I5,I5,3(1X,A14))
           ENDIF 
          ENDIF 
         ENDIF 
!
! *********************************************************************
! *   DETERMINE THE NUMBER OF REACTIONS WITH ZERO ACTIVE LOSS TERMS   * 
! *********************************************************************
! NOLOSRAT = NUMBER OF ACTIVE REACTIONS WITH NO LOSS TERMS 
! NOLOSRN  = REACTION NUMBER OF EACH REACTION WITH NO LOSS TERMS
!

        ELSEIF (IAL.EQ.0) THEN
         NOLOSRAT(NCS)         = NOLOSRAT(NCS) + 1
         NOL                   = NOLOSRAT(NCS)
         NOLOSRN(NOL,NCS)      = NK
        ENDIF
!       ENDIF IAL.GT.0 
!
! *********************************************************************
! * COUNT GROSS AND NET PRODUCTION AND SET A PARTIAL DERIVATIVE ARRAY * 
! *********************************************************************
! NUMGAINT = EVERY OCCURENCE OF A PRODUCTION (ACTIVE & INACTIVE SPEC) 
! NUMGAIN  = EVERY NET OCCURENCE OF A PRODUCTION WHERE THE SPECIES IS 
!            NOT LOST IN THE SAME REACTION. (ACTIVE & INACTIVE SPEC)
! IAPROD   = NUMBER OF ACTIVE PRODUCTS IN EACH NK REACTION. USED
!            TO CALCULATE PARTIAL DERIVATIVES IN PDERIV.F. 
! IRM2     = NEW SPECIES # OF EACH ACTIVE PRODUCT IN EACH NK REACTION
!
        IAP                       = NPRODLO - 1
        DO 210 K                  = NPRODLO, NPRODHI  
         IPROD                    = IRM(K,NK,NCS)
         IF (IPROD.GT.0) THEN
          IPR                     = MAPPL(IPROD,NCS)
          RFRAC                   = FKOEF(K,NK,NCS)
          LFRAC                   = INT(RFRAC + SMAL1) 
          ALFRAC                  = FLOAT(LFRAC)
          DIFF                    = ABS(RFRAC-ALFRAC)
!
! ******************** PRODUCTION TERM IS A FRACTION ******************
!
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
!
! ******************* PRODUCTION TERM IS NON-FRACTION *****************
!
          ELSE
           APORL(IPR)             = APORL(IPR) + RFRAC
           KPRODS                 = LFRAC
           NUMGAINT(IPR,NCS)      = NUMGAINT(IPR,NCS) + LFRAC
           FKOEF(K,NK,NCS)        = 1.d0
          ENDIF  
!
! ******************* IDENTIFY ALL PRODUCTION TERMS *******************
!
          IF (IPR.LE.NSPEC(NCS)) THEN
           DO 170 L               = 1, KPRODS
            IAP                   = IAP + 1
            IAPROD(NK,NCS)        = IAP
            IRM2(IAP,NK,NCS)      = IPR
            FK2( IAP,NK,NCS)      = FKOEF(K,NK,NCS)
 170       CONTINUE  
          ENDIF
!
         ENDIF
!
 210    CONTINUE
!
! *********************************************************************
! *  FIND NET PROD AND LOSS TERMS FOR ALL BUT FRACTIONATED PRODUCTS   * 
! *********************************************************************
!
         DO 220 JSPC              = 1, NTSPEC(NCS)
          IF (ABS(APORL(JSPC)).LT.SMAL1) THEN
           KDIF                   = 0 
!
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
!
          IF (NK.LE.NRATES(NCS)) THEN 
           NUMLOSS(JSPC,NCSP)     = NUMLOSS(JSPC,NCS)
           NUMGAIN(JSPC,NCSP)     = NUMGAIN(JSPC,NCS)
           NUMPORL(JSPC,NCSP)     = NUMPORL(JSPC,NCS)
          ENDIF 
!
 220     CONTINUE 
!
         IF (NK.LE.NRATES(NCS)) THEN 
          NOLOSRAT(NCSP)          = NOLOSRAT(NCS) 
          NGNFRAC( NCSP)          = NGNFRAC( NCS)
          IONER(   NCSP)          = IONER(   NCS)
         ENDIF 
!
         NKLAST                   = NK 
!
 230   CONTINUE
!      CONTINUE N = 1, NTRATES
!
! *********************************************************************
! * SET ARRAY FOR REORDERING RATES FROM 3..2..1..0 BODY REACTIONS     *
! *********************************************************************
! INOREP   = LAST REORDERED REACTION NUMBER PRIOR TO SETS OF TWO
!            REACTIONS WITH TWO REACTANTS  
! NOLDFNEW = OLD REACTION RATE # CORRESP. TO EACH REORDERED REACTION
! NEWFOLD  = NEW REACTION RATE # CORRESP. TO EACH ORIGINAL RATE NUMBER
!
       IC                 = 0
       DO 235 I           = 1, ITHRR(NCS)   
        IC                = IC + 1
        NK                = NKTHRR(I,NCS)
        NK1               = NK + NTRATES(NCS) 
        NOLDFNEW(IC, NCS) = NK
        NEWFOLD( NK, NCS) = IC
        NEWFOLD( NK1,NCS) = IC + NALLRAT(NCS) 
 235   CONTINUE 
!
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
!
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
!
       DO 241 I            = 1, NOLOSRAT(NCS)
        IC                 = IC + 1
        NK                 = NOLOSRN(I,NCS)
        NK1                = NK + NTRATES(NCS)
        NOLDFNEW(IC, NCS)  = NK
        NEWFOLD( NK, NCS)  = IC
        NEWFOLD( NK1,NCS)  = IC + NALLRAT(NCS)
 241   CONTINUE
!
       IF (IC.NE.NALLRAT(NCS)) THEN
        WRITE(6,245) IC, NALLRAT(NCS)
        CALL GEOS_CHEM_STOP
       ENDIF
!
! *********************************************************************
!                SET A SLIGHTLY MORE EFFICIENT PHOTO ARRAY 
! *********************************************************************
!
       DO 243 J          = 1, JPHOTRAT(NCS)
         NKNPHOTRT(J,NCS) = NKPHOTRAT(J,NCS)
 243   CONTINUE
!
 245   FORMAT('JSPARSE: IC NE NALLRAT =',2(I5))
!
! *********************************************************************
! ****** DETERMINE NUMBER OF SPECIES WITH GROSS/NET LOSSES/GAINS ******
! *********************************************************************
! NSPCSOLV = # OF ACTIVE SPECIES WITH AT LEAST ONE GROSS LOSS
! ISOLVSPC = SPECIES NUMBER OF EACH NSPCSOLV SPECIES
! ISGAINR  = # OF ACTIVE SPECIES WITH AT LEAST ONE NET CHEM GAIN 
! IGAINR   = SPECIES NUMBER OF EACH ISGAINR SPECIES
! ISGAINE  = # OF ACTIVE SPECIES WITH AT LEAST 1 NET CHEM GAIN 
! IGAINR   = SPECIES NUMBER OF EACH ISGAINR SPECIES
! NOGAINE  = # OF ACTIVE SPECIES WITH ZERO NET CHEM OR GAINS 
! NGAINE   = SPECIES NUMBER OF EACH NOGAINE SPECIES
! ISPORL   = # OF ACTIVE SPECIES WITH AT LEAST ONE NET PRODUCTION
!            OR LOSS TERM FOR SMVGEAR.
! IPORL    = SPECIES NUMBER OF EACH ISPORL SPECIES
!
       DO 300 JOLD             = 1, NSPEC(NCS) 
        JNEW                   = MAPPL(JOLD,NCS)
!
        IF (NUMGAIN(JNEW,NCS).GT.0) THEN
         ISGAINR(NCS)          = ISGAINR(NCS) + 1
         IGR                   = ISGAINR(NCS)
         IGAINR(IGR,NCS)       = JNEW  
        ENDIF 
!
        IF (NUMPORL(JNEW,NCS).GT.0) THEN
         ISPORL(NCS)           = ISPORL(NCS) + 1 
         ISP                   = ISPORL(NCS)
         IPORL(ISP,NCS)        = JNEW 
        ENDIF
!
        IF (NUMLOST(JNEW,NCS).GT.0) THEN
         NSPCSOLV(NCS)         = NSPCSOLV(NCS) + 1
         NSP                   = NSPCSOLV(NCS)
         ISOLVSPC(NSP,NCS)     = JNEW  
        ENDIF
!
        IF (NUMGAIN(JNEW,NCS).GT.0.OR.FRACGAIN(JNEW,NCS).GT.0) THEN
         ISGAINE(NCS)         = ISGAINE(NCS) + 1
         IGR                  = ISGAINE(NCS)
         IGAINE(IGR,NCS)      = JNEW  
        ELSEIF (NUMLOSS(JNEW,NCS).GT.0) THEN 
         NOGAINE(NCS)         = NOGAINE(NCS) + 1
         NGR                  = NOGAINE(NCS)
         NGAINE(NGR,NCS)      = JNEW    
        ENDIF
!
 300   CONTINUE
!
! *********************************************************************
! ********  CHECK DIMENSIONS RESULTING FROM GAINS AND LOSSES  *********
! *********************************************************************
!
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
!
       IF (IOREAC.EQ.1) THEN
        IF(am_I_Root) WRITE(IO93,*)
        IF(am_I_Root) WRITE(IO93,240)
        DO 270 K = 1, NTSPEC(NCS) 
         J       = INEWOLD(K,NCS) 
         IF ( am_I_Root ) THEN
            WRITE(IO93,250)
     &           NAMENCS( J,NCS),NUMGAINT(K,NCS),NUMGAIN( K,NCS),
     &           NUMLOST( K,NCS),NUMLOSS( K,NCS),NUMGAINT(K,NCS) 
     &           -NUMLOST( K,NCS)-NUMGAIN( K,NCS)+NUMLOSS( K,NCS),
     &           FRACGAIN(K,NCS),NUMGFRT( K,NCS)
         ENDIF
 270    CONTINUE
        IF ( am_I_Root ) THEN
           WRITE(IO93,250) 'OVERALL       ',
     &          NGTSUM, NGSUM, NLTSUM, NLSUM,
     &          NGTSUM - NLTSUM - NGSUM + NLSUM, 0., NGFSUM
        ENDIF
       ENDIF
!
       IF (NMOTH(  NCS).GT.MAXGL2.OR.NOLOSP(NCS).GT.MAXGL3.OR.
     1     NGNFRAC(NCS).GT.MAXGL) THEN
        WRITE(6,275) MAXGL2, NMOTH(  NCS), MAXGL3, NOLOSP(NCS),
     1               MAXGL,  NGNFRAC(NCS)
        CALL GEOS_CHEM_STOP
       ENDIF
!
! *********************************************************************
! *       CHECK WHETHER CHEMICAL SYSTEM IS ATOM-CONSERVATIVE          *
! *********************************************************************
! JMBCOMP = SPECIES NUMBER FOR EACH SPECIES IN A MASS BAL. GROUP
! MBCOMP  = COUNTS THE NUMBER OF MASS BALANCE SPECIES IN EACH M.B GROUP
! NMASBAL = NUMBER OF MASS BALANCE GROUPS (E.G. S, N, C ARE GROUPS)
! WTMB(1) = NUMBER OF ATOMS OF A GIVEN MASS BALANCE SPECIES PER MOLECULE 
!
       IF(am_I_Root) WRITE(IO93,360) CHEMTYP(NCS)
!
       IF (NCS.LE.NCSGAS) THEN
! 
! ----------------------------   GAS-PHASE   -------------------------- 
!
        DO 385 N     = 1, NMASBAL 
         IF (MBCOMP(N,MB1).GT.0) THEN 
          TNUMGN     = 0
          TNUMLS     = 0
          IF(am_I_Root) WRITE(IO93,325) NAMEMB(N)
          DO 380 J   = 1, MBCOMP(N,MB1) 
           JGAS      = JMBCOMP(N,J,MB1)
           JNEW      = MAPPL(JGAS,NCS)
           SUMGN     = NUMGAIN(JNEW,NCS) + FRACGAIN(JNEW,NCS) 
           TNUMGNA   = SUMGN             * WTMB(N,JGAS,MB1)  
           TNUMLSA   = NUMLOSS(JNEW,NCS) * WTMB(N,JGAS,MB1) 
           TNUMGN    = TNUMGN + TNUMGNA
           TNUMLS    = TNUMLS + TNUMLSA
           IF ( am_I_Root ) THEN
              WRITE(IO93,320) NAMEGAS(JGAS), TNUMGNA, TNUMLSA, 0 
           ENDIF
 380      CONTINUE
          IF ( am_I_Root ) THEN
             WRITE(IO93,370) TNUMGN, TNUMLS, TNUMGN - TNUMLS 
          ENDIF
         ENDIF 
 385    CONTINUE
       ENDIF
!
       IF(am_I_Root) WRITE(IO93,375) 
     &                 NALLRAT(NCSP), NALLRAT(NCS) - NALLRAT(NCSP),
     &                 NALLRAT(NCS) 
!
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
!
! *********************************************************************
! *********************************************************************
! **        SET ARRAYS TO TAKE ADVANTAGE OF SPARSE MATRICES          ** 
! *********************************************************************
! *********************************************************************
!
! IFSUN  = 1 THEN DAY-CHEMISTRY;  = 2 THEN NIGHT CHEMISTRY
! NCSP   = NCS       FOR DAYTIME   TROP-GAS, STRAT-GAS CHEM  
! NCSP   = NCS + ICP FOR NIGHTTIME TROP-GAS, STRAT-GAS CHEM  
!
! LZERO    = 1 IF AN ARRAY SPOT IS FILLED WITH A NON-ZERO VALUE. LZERO
!            IS UPDATED AS WE SIMULATE THE ORDER OF CALCULATIONS DURING
!            A PRACTICE L-U DECOMPOSITION
! MXGSAER  = LARGER OF IGAS, IAERTY
!
!
      IF (IFNONE.EQ.0) THEN
       IFNONE                 = 1
       NPLFUN                 = 0 
       NFRCOUN                = 0 
       NPDCOUN                = 0 
       NPLTOT                 = 0 
      ENDIF
!
      DO 700 IFSUN            = 1, 2 
       NCSP                   = (IFSUN - 1) * ICS + NCS
!
       DO 517 I               = 1, MXGSAER
        DO 515 J              = 1, MXGSAER
         LZERO(J,I)           = 0
 515    CONTINUE
        LZERO(I,I)            = 1
 517   CONTINUE 
!
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
!
! *********************************************************************
! *   SET DECOMPOSITION AND BACK-SUBSTITUTION SPARSE-MATRIX ARRAYS    *
! *********************************************************************
!
       CALL KSPARSE(am_I_Root)

!
! *********************************************************************
! *    SET ARRAYS TO IMPROVE EFFICIENCY OF FIRST-DERIVATIVE CALCS     * 
! *********************************************************************
! *********************************************************************
! **   SET ARRAYS FOR KINETIC AND PHOTO PRODUCTION AND LOSS RATES    **
! *********************************************************************
!
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
!
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
!
        DO 649 N           = 1, KNUMPORL 
         NK                = JPORL(JSPC,N,NCS) 
         NEWNK(N)          = NEWFOLD(NK,NCS)  
 649    CONTINUE 
!
        DO 651 MC          = NPL5(NPLTOT), NPH5(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         LOSSRB(MC)        = NEWNK(NCCOUNT+2) 
         LOSSRC(MC)        = NEWNK(NCCOUNT+3) 
         LOSSRD(MC)        = NEWNK(NCCOUNT+4) 
         LOSSRE(MC)        = NEWNK(NCCOUNT+5) 
         NCCOUNT           = NCCOUNT + 5
 651    CONTINUE 
!
        DO 652 MC          = NPL4(NPLTOT), NPH4(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         LOSSRB(MC)        = NEWNK(NCCOUNT+2) 
         LOSSRC(MC)        = NEWNK(NCCOUNT+3) 
         LOSSRD(MC)        = NEWNK(NCCOUNT+4) 
         NCCOUNT           = NCCOUNT + 4  
 652    CONTINUE 
!
        DO 653 MC          = NPL3(NPLTOT), NPH3(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         LOSSRB(MC)        = NEWNK(NCCOUNT+2) 
         LOSSRC(MC)        = NEWNK(NCCOUNT+3) 
         NCCOUNT           = NCCOUNT + 3   
 653    CONTINUE 
!
        DO 654 MC          = NPL2(NPLTOT), NPH2(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         LOSSRB(MC)        = NEWNK(NCCOUNT+2) 
         NCCOUNT           = NCCOUNT + 2    
 654    CONTINUE 
!
        DO 656 MC          = NPL1(NPLTOT), NPH1(NPLTOT)
         LOSSRA(MC)        = NEWNK(NCCOUNT+1) 
         NCCOUNT           = NCCOUNT + 1     
 656    CONTINUE 
!
 670   CONTINUE 
       NPLHI(NCSP)         = NPLTOT
!
! *********************************************************************
! *              SET ARRAY FOR FRACTIONATED PRODUCTS                  *  
! *********************************************************************
!
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
!
! *********************************************************************
! * SET ARRAYS TO IMPROVE EFFICIENCY OF PARTIAL DERIVATIVE CALCS      * 
! *********************************************************************
!
       NPDLO(NCSP)           = NPDCOUN + 1
!
       DO 974 NA             = 1, NALLRAT(NCSP) 
        NK                   = NCEQUAT(NA,NCS) 
        IHIREAC              = NRUSE(  NK,NCS) 
!
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
!
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
!
       NPDHI(NCSP)          = NPDCOUN
!
! *********************************************************************
! **        CHECK DIMENSIONS AND PRINT OUT ARRAY SAVINGS             ** 
! *********************************************************************
!
       IF (NPLTOT   .GT. MXCOUNT4  .OR. NPLFUN   .GT. MXCOUNT4 .OR.
     3     NFRCOUN  .GT. MXCOUNT4 .OR.  NPDCOUN  .GT. MXCOUNT2) THEN
        WRITE(6,645) MXCOUNT4, NPLTOT,    MXCOUNT4, NPLFUN,
     2               MXCOUNT4, NFRCOUN,   MXCOUNT2, NPDCOUN
        CALL GEOS_CHEM_STOP
       ENDIF
!
 700  CONTINUE
!     CONTINUE IFSUN = 1, 2
!
 645  FORMAT('ONE OF THE DIMENSIONS BELOW IS TOO SMALL:',/,
     1       'DIMENSION: MXCOUNT4 =  ',I5,' VARIABLE: NPLTOT   = ',I5,/,
     2       'DIMENSION: MXCOUNT4 =  ',I5,' VARIABLE: NPLFUN   = ',I5,/,
     3       'DIMENSION: MXCOUNT4 =  ',I5,' VARIABLE: NFRCOUN  = ',I5,/,
     4       'DIMENSION: MXCOUNT2 =  ',I5,' VARIABLE: NPDCOUN  = ',I5)
!
! *********************************************************************
! ********************** END OF SUBROUTINE JSPARSE ********************
! *********************************************************************
!
      RETURN                                                             
      END SUBROUTINE JSPARSE
