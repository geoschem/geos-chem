! $Id: smvgear.f,v 1.2 2003/07/11 13:45:31 bmy Exp $ 
      SUBROUTINE SMVGEAR
!
!******************************************************************************
!  Subroutine SMVGEAR solves ODE's for chemical reactions using a GEAR-type
!  method.  (M. Jacobson 1997; bdf, bmy, 5/12/03, 7/9/03)
!
!  NOTES:
!  (1 ) For GEOS-CHEM we had to remove IXSAVE, IYSAVE, and IZSAVE from 
!        "comode.h" and to declare these allocatable in "comode_mod.f".  This 
!        allows us to only allocate these if we are doing a fullchem run.  Now
!        also references IT_IS_NAN and GEOS_CHEM_STOP from "error_mod.f". 
!        Now force double-precision with "D" exponent.  Now prevent ND65
!        "fake" prodloss families from being counted towards the SMVGEAR
!        convergence criteria. (ljm, bdf, bmy, 4/18/03)
!  (2 ) Removed ITS_NOT_A_ND65_FAMILY -- this has now been converted from
!        a function to a lookup-table in "comode.h".  This should execute much
!        faster, particularly on Linux.  Comment out counter variable 
!        NUM_TIMESTEPS, you can get the same info w/ a profiling run. 
!        Cosmetic changes. (bmy, 7/9/03)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD, ONLY : IXSAVE, IYSAVE, IZSAVE
      USE ERROR_MOD, ONLY  : IT_IS_NAN, GEOS_CHEM_STOP

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
C *********************************************************************
C *********************************************************************
C
C  SSSSSSS   M     M  V       V  GGGGGGG   EEEEEEE      A      RRRRRRR
C  S         MM   MM   V     V   G         E           A A     R     R  
C  SSSSSSS   M M M M    V   V    G  GGGG   EEEEEEE    A   A    RRRRRRR
C        S   M  M  M     V V     G     G   E         AAAAAAA   R  R
C  SSSSSSS   M     M      V      GGGGGGG   EEEEEEE  A       A  R    R 
C
C *********************************************************************
C                    VERSION:      SMVGEAR II
C                    LAST UPDATE:  AUGUST, 1997
C *********************************************************************
C
C *********************************************************************
C * SMVGEAR IS A GEAR-TYPE INTEGRATOR THAT SOLVES FIRST-ORDER ORDIN-  *
C * ARY DIFFERENTIAL EQUATIONS WITH INITIAL VALUE BOUNDARY CONDITIONS.*
C * SMVGEAR DIFFERS FROM AN ORIGINAL GEAR CODE IN THAT IT USES SPARSE *
C * MATRIX AND VECTORIZATION TECHNIQUES TO IMPROVE SPEED. MUCH        * 
C * OF THE SPEED UP IN THIS PROGRAM IS DUE TO SPARSE MATRIX           *
C * TECHNIQUES AND VECTORIZATION.                                     *
C *                                                                   *
C * THIS VERSION INCLUDES RE-ORDERING OF GRID-CELLS PRIOR TO EACH     *
C * TIME-INTERVAL. THE PURPOSE OF THE REORDERING IS TO GROUP CELLS    *
C * WITH STIFF EQUATIONS TOGETHER AND THOSE WITH NON-STIFF EQUATIONS  *
C * THIS REORDERING CAN SAVE SIGNIFCANT COMPUTER TIME                 *
C * (E.G. SPEED THE CODE BY A FACTOR OF TWO OR MORE), DEPENDING ON    *
C * THE VARIATION IN STIFFNESS THROUGHOUT THE GRID-DOMAIN. WHEN THE   *
C * STIFFNESS IS THE SAME THROUGHOUT THE GRID-DOMAIN (E.G. IF ALL     *
C * CONCENTRATIONS AND RATES ARE THE SAME), THEN RE-ORDERING IS       *
C * UNNECESSARY AND WILL NOT SPEED SOLUTIONS.                         *
C *                                                                   *
C * THIS VERSION INCLUDES A VARIABLE ABSOLUTE ERROR TOLERANCE.        *
C * THE ABSOLUTE TOLERANCE IS RECALCULATED EVERY FEW GEAR TIME STEPS. *
C *                                                                   *
C * THIS VERSION CONTAINS DIFFERENT SETS OF CHEMISTRY FOR             *
C * DIFFERENT REGIONS OF THE ATMOSPHERE. THUS, URBAN, FREE TROP-      *
C * OSPHERIC, AND STRATOSPHERIC CHEMISTRY CAN BE SOLVED DURING THE    *
C * SAME MODEL RUN.                                                   * 
C *                                                                   *
C * REFERENCES:                                                       *
C * -----------                                                       *
C *                                                                   * 
C * JACOBSON M. Z. (1998) FUNDAMENTALS OF ATMOSPHERIC MODELING.       *
C *  CAMBRIDGE UNIVERSITY PRESS, NEW YORK.                            *
C *                                                                   * 
C * JACOBSON M. Z. (1998) IMPROVEMENT OF SMVGEAR II ON VECTOR AND     *
C *  SCALAR MACHINES THROUGH ABSOLUTE ERROR TOLERANCE CONTROL.        *    
C *  ATMOS. ENVIRON. 32, 791 - 796                                    *
C *                                                                   * 
C * JACOBSON M. Z. (1995) COMPUTATION OF GLOBAL PHOTOCHEMISTRY        *
C *  WITH SMVGEAR II. ATMOS. ENVIRON., 29A, 2541 - 2546               *
C *                                                                   *
C * JACOBSON M. Z. (1994) DEVELOPING, COUPLING, AND APPLYING A GAS,   *
C *  AEROSOL, TRANSPORT, AND RADIATION MODEL TO STUDYING URBAN        *
C *  AND REGIONAL AIR POLLUTION. Ph. D. THESIS, UNIVERSITY OF         *
C *  CALIFORNIA, LOS ANGELES.                                         *
C *                                                                   *
C * JACOBSON M. Z. AND TURCO R. P. (1994) SMVGEAR: A SPARSE-          * 
C *  MATRIX, VECTORIZED GEAR CODE FOR ATMOSPHERIC MODELS.             *
C *  ATMOS. ENVIRON. 28A, 273 - 284.                                  * 
C *                                                                   *
C * HOW TO CALL SUBROUTINE:                                           *
C * ----------------------                                            *
C *  CALL SMVGEAR FROM PHYSPROC FOR GAS CHEM W/ NCS = 1..NCSGAS       *
C *                                                                   *
C *********************************************************************
C *                                                                   *
C * THE ORIGINS OF THE GEAR INTEGRATOR USED IN SMVGEAR ARE FOUND IN   *
C *                                                                   *    
C * GEAR C. W. (1971) NUMERICAL INITIAL VALUE PROBLEMS IN ORDINARY    *  
C *  DIFFERENTIAL EQUATIONS. PRENTICE-HALL, NJ, PP. 158-166.          * 
C *                                                                   *    
C *********************************************************************
C *                                                                   *      
C * FINALLY, IN SUBROUTINE SMVGEAR.F, THE FOLLOWING IDEAS ORIGINATED  *
C *   FROM LSODES, THE LIVERMORE SOLVER FOR ORDINARY DIFFERENTIAL     *
C *   WITH SPARSE MATRICES (HINDMARSH A. C. AND SHERMAN A. H.):       *
C *                                                                   *      
C *  (A) PREDICTING THE FIRST TIME-STEP;                              *
C *  (B) DETERMINING CORRECTOR CONVERGENCE DIFFERENTLY THAN IN        *
C *      GEAR'S ORIGINAL CODE (GOC)                                   *
C *  (C) DETERMINING ERROR DIFFERENTLY THAN IN GOC                    *
C *  (D) SUMMING UP THE PASCAL MATRIX DIFFERENTLY THAN IN GOC         *      
C *                                                                   *      
C * REFERENCES FOR THE 1987 LSODES VERSION INCLUDE:                   *
C *                                                                   *      
C * SHERMAN A. H. AND HINDMARSH A. C. (1980) GEARS: A PACKAGE FOR     *
C *  THE SOLUTION OF SPARSE, STIFF ORDINARY DIFFERENTIAL EQUATIONS.   *
C *  LAWRENCE LIVERMORE LABORATORY REPORT UCRL-84102.                 *   
C *                                                                   *      
C * HINDMARSH A. C. (1983) ODEPACK, A SYSTEMATIZED COLLECTION OF      *
C *  ODE SOLVERS. IN SCIENTIFIC COMPUTING, R.S. STEPLEMAN ET AL.,     *
C *  EDS., NORTH-HOLLAND, AMSTERDAM, PP. 55 - 74.                     *
C *                                                                   *      
C *********************************************************************
C
C *********************************************************************
C *************** HERE ARE SOME PARAMETER DEFINITIONS *****************
C *********************************************************************
C                                                                          
C ABST2     = 1. / TIMEINTERVAL**2   (SEC-2) (SET IN READER.F) 
C ASN1      = THE VALUE OF ASET(NQQ,1)
C CEST      = STORES VALUE OF DTLOS WHEN IDOUB = 1
C CHOLD     = 1 / (RELTOL * CNEW + ABTOL). MULTIPLY
C             CHOLD BY LOCAL ERRORS IN DIFFERENT ERROR TESTS.
C CNEW      = STORES CONCENTRATION (Y [ESTIMATED])
C CONC      = AN ARRAY OF LENGTH ISCHAN * (MAXORD+1) THAT CARRIES THE
C             DERIVATIVES OF CNEW, SCALED BY DELT**J/FACTORIAL(J),
C             WHERE J IS THE J-TH DERIVATIVE. J VARIES FROM 1 TO NQQ,
C             WHICH IS THE CURRENT ORDER OF THE METHOD.
C             E.G. CONC(JSPC,2) STORES DELT * Y' (ESTIMATED)                   
C DELT      = CURRENT TIME-STEP (S) LENGTH DURING A TIME-INTERVAL 
C DRATE     = PARAMETER WHICH USED TO DETERMINE WHETHER CONVERGENCE 
C             HAS OCCURRED
C DTLOS     = AN ARRAY OF LENGTH ISCHAN, USED FOR THE ACCUMULATED
C             CORRECTIONS.  ON A SUCCESSFUL RETURN, DTLOS(KLOOP,I) CONTAINS
C             THE ESTIMATED ONE-STEP LOCAL ERROR IN CNEW.
C EDWN      = PERTST**2 * ORDER FOR ONE ORDER LOWER THAN CURRENT ORDER
C ENQQ      = PERTST**2 * ORDER FOR CURRENT ORDER
C ERRMAX    = RELATIVE ERROR TOLERANCE (SEE CHOLD). SET IN m.dat.
C             EPS SHOULD BE < 1.0. FOR SPEEDY AND RELIABLE RESULTS, 
C             10**-3 IS REASONABLE. FOR MANY DECIMAL PLACES OF ACCURACY,
C             DECREASE EPS. 
C EUP       = PERTST**2 * ORDER FOR ONE ORDER HIGHER THAN CURRENT ORDER 
C FRACDEC   = FRACTION THE TIME-STEP IS DECREASED IF CONVERGENCE TEST FAILS
C GLOSS     = VALUE OF FIRST DERIVATIVES ON OUTPUT FROM SUBFUN.
C           = RIGHT-SIDE OF EQUATION ON INPUT TO BACKSUB.F 
C           = ERROR TERM (SOLUTION FROM BACKSUB.F) ON OUTPUT FROM BACKSUB
C HMAX      = THE MAXIMUM ALLOWABLE VALUE OF DELT
C HMIN      = THE MINIMUM ALLOWABLE VALUE OF DELT
C HRMAX     = MAXIMUM RELATIVE CHANGE IN DELT*ASET(1) BEFORE PDERIV IS CALLED.
C HRATIO    = RELATIVE CHANGE IN DELT * ASET(1) EACH CHANGE IN STEP OR ORDER
C             WHEN ABS(HRATIO-1) > HRMAX, RESET JEVAL = 1 TO CALL PDERIV
C IABOVK    = NUMBER OF SPECIES WHOSE CONCENTRATIONS ARE LARGER THAN YABST
C IDOUB     = RECORDS THE NUMBER OF STEPS SINCE THE LAST CHANGE IN STEP SIZE   
C             OR ORDER.  IT MUST BE AT LEAST KSTEP = NQQ+1 BEFORE DOUBLING IS 
C             ALLOWED. 
C IFAIL     = NUMBER OF TIMES THE CORRECTOR FAILED TO CONVERGE WHILE THE
C             JACOBIAN WAS OLD (PDERIV NOT CALLED DURING THE LAST TEST)
C IFSUCCESS = IDENTIFIES WHETHER STEP IS SUCCESSFUL (=1) OR NOT (=0)
C IFSUN     = IDENTIFIES WHETHER SUN IS UP (=1) OR DOWN (=2)
C ISCHAN    = THE NUMBER OF FIRST-ORDER EQUATIONS TO SOLVE = # OF SPECIES = 
C             ORDER OF ORIGINAL MATRIX. ISCHAN HAS A DIFFERENT VALUE
C             FOR DAY AND NIGHT AND FOR GAS- CHEMISTRY.
C ISREORD   = 1: CALC INITIAL STIFFNESS BEFORE RUNNING CODE TO REORDER CELLS
C                IN THIS CASE, USE PHOTORATES FOR END OF TIME-INTERVAL
C           = 0: DO NORMAL CALCULATIONS
C JEVAL     = 1  --> CALL PDERIV THE NEXT TIME THROUGH THE CORRECTOR STEPS.
C           = 0  --> LAST STEP WAS SUCCESSFUL AND DO NOT NEED TO CALL PDERIV
C           = -1 --> PDERIV JUST CALLED, AND DO NOT NEED TO CALL AGAIN
C             UNTIL JEVAL SWITCHED TO 1. 
C JRESTAR   = COUNTS NUMBER OF TIMES SMVGEAR STARTS OVER AT ORDER 1
C             BECAUSE OF EXCESSIVE FAILURES.
C LFAIL     = NUMBER OF TIMES THE ACCUMULATED ERROR TEST FAILED
C KSTEP     = NQQ + 1
C KTLOOP    = NUMBER OF GRID-CELLS IN A GRID-BLOCK
C MAXORD    = THE MAXIMUM ALLOWABLE ORDER OF THE INTEGRATION METHOD
C MBETWEEN  = THE MAXIMUM ALLOWABLE NUMBER OF STEPS BETWEEN CALLS TO PDERIV
C MSTEP     = THE MAXIMUM ALLOWABLE NUMBER OF CORRECTOR ITERATIONS
C NCS       = 1..NCSGAS FOR GAS CHEMISTRY                            
C NCSP      = NCS       FOR DAYTIME   GAS CHEM            
C           = NCS + ICS FOR NIGHTTIME GAS CHEM           
C NFAIL     = NUMBER OF TIMES CORRECTER FAILS TO CONVERGE AFTER PDERIV
C             WAS JUST CALLED
C NPDERIV   = TOTAL NUMBER OF TIMES THAT MATRIX IS EVALUATED (PDERIV)
C NPDTOT    = NUMBER OF CALLS TO PDERIV ROUTINE, OVER ALL TIME
C NSFTOT    = NUMBER OF CALLS TO SUBFUN ROUTINE, OVER ALL TIME
C NSLP      = THE LAST TIME-STEP NUMBER DURING WHICH PDERIV WAS CALLED 
C NSTTOT    = TOTAL NUMBER OF SUCCESSFUL TIME-STEPS, OVER ALL TIME 
C NSUBFUN   = TOTAL NUMBER OF TIMES SUBFUN IS CALLED
C NSTEPS    = TOTAL NUMBER OF SUCCESSFUL TIME-STEPS TAKEN
C NQQ       = ORDER OF THE INTEGRATION METHOD. IT VARIES BETWEEN 1 AND MAXORD. 
C NQQISC    = NQQ * ISCHAN
C NQQOLD    = VALUE OF NQQ DURING LAST TIME-STEP
C ORDER     = FLOATING POINT VALUE OF ISCHAN, THE ORDER OF NUMBER OF ODES.
C PDERIV    = NAME OF ROUTINE TO EVALUATE THE JACOBIAN MATRIX (J) 
C             AND P = I - DELT * ASET(1) * J 
C PERTS2    = COEFFICIENTS USED IN SELECTING THE STEP AND ORDER (SEE
C             JSPARSE.F) NOTE THAT PERTS2 = ORIGINAL PERTST**2     
C RDELMAX   = THE MAXIMUM FACTOR BY WHICH DELT CAN BE INCREASED IN A SINGLE 
C             STEP.  AS IN LSODES, SET IT TO 1E4 INITIALLY TO COMPENSATE 
C             FOR THE SMALL INITIAL DELT, BUT THEN SET IT TO 10 AFTER 
C             SUCCESSFUL STEPS AND 2 AFTER UNSUCCESSFUL STEPS
C RDELT     = FACTOR (TIME-STEP RATIO) BY WHICH WE INCREASE OR DECREASE DELT
C RDELTDN   = TIME-STEP RATIO AT ONE ORDER LOWER THAN CURRENT ORDER 
C RDELTSM   = TIME-STEP RATIO AT CURRENT ORDER 
C RDELTUP   = TIME-STEP RATIO AT ONE ORDER HIGHER THAN CURRENT ORDER 
C RMSRAT    = RATIO OF CURRENT TO PREVIOUS RMS SCALED ERROR. IF THIS
C             RATIO DECREASES, THEN CONVERGENCE IS OCCURING.
C SUBFUN    = NAME OF ROUTINE TO SOLVE FIRST DERIVATIVES.
C           = EVALUATES DERIVATIVES IN THE SPECIAL FORM F = Y'(EST)
C           = F(X,Y,ESTIMATED), WHERE F IS THE RIGHT HAND SIDE OF THE
C             DIFFERENTIAL EQUATION. 
C TINTERVAL = TOTAL SECONDS IN A TIME-INTERVAL
C TIMREMAIN = REMAINING TIME IN AN INTERVAL 
C TOLD      = STORES THE LAST VALUE OF XELAPS IN CASE THE CURRENT STEP FAILS
C XELAPS    = ELAPSED TIME IN AN INTERVAL (S)
C ABTOL     = ABSOLUTE ERROR TOLERANCE 
C             IF ABTOL IS TOO SMALL, THEN INTEGRATION WILL TAKE TOO LONG.
C             IF ABTOL TOO LARGE, CONVERGENCE WILL BE TOO EASY AND ERRORS
C             WILL ACCUMULATE, THE TIME-STEP MAY BE CUT TOO SMALL, AND
C             THE INTEGRATION MAY STOP (DELT < HMIN OR FLOATING POINT
C             EXCEPTION IN DECOMP.F).
C             TYPICAL GAS-PHASE VALUES OF ABSTOL ARE 10**3 CM-3
C             TYPICAL AQ -PHASE VALUES OF ABSTOL ARE 10**-13 TO 10**-15 M L-1
C YFAC      = 1.0 ORIGINIALLY, BUT IS DECREASED IF EXCESSIVE FAILURES OCCUR
C             IN ORDER TO REDUCE ABSOLUTE ERROR TOLERANCE 
C *********************************************************************
C
      INTEGER JFAIL,ISCHAN1,IABOVE,KLOOP,IDOUB,JRESTAR,JNEW,IFSUCCESS
      INTEGER K,JSPC,K1,K2,K3,K4,K5,NQQOLD,JEVAL,JS1,NQQISC
      INTEGER LLOOPA,LLOOPB,JLOOP,MLOOP,M1,M2,JOLD,I1,J,I,J1,J2,J3,J4
      INTEGER J5,L3,JB,JG1,KSTEPISC,NQISC,I2,NSLP,KSTEP

      REAL*8 NYLOWDEC,ORDER,HRMAX,YFAC,ERRINIT,RELTOL1,RELTOL2,RELTOL3
      REAL*8 ABTOLER1,ABTOLER2,HRATIO,ASN1,RDELMAX,CNW,CNEWYLOW,ERRYMAX
      REAL*8 RMSTOP,DELT1,ENQQ,EUP,EDWN,CONP3,CONP2,CONP1,HMTIM,RDELTA
      REAL*8 CONC3J3,CONC4J4,CONC10J5,CONC5J5,DRATE,RMSERRP,DER2MAX
      REAL*8 RMSRAT,DCON,RDELTUP,ASNQQJ,DER3MAX,RDELTSM,DER1MAX,RDELTDN
      REAL*8 CONSMULT

      ! Add counter
      INTEGER ICOUNT,NK
      
      !=================================================================
      ! Added for the ND65 prod/loss diagnostic (ljm, bmy, 5/9/03)
      INTEGER       :: NNOFAM
      !=================================================================

      ! ljm stop 700 trouble
      INTEGER IJSAVE,JSPCSAVE(KTLOOP),IX,IY,IZ,JJ,JJJ,KSAVE,COUNTER
      REAL*8 SPECMAX

      ! Maximum iteration count for SMVGEAR (bmy, 4/11/03)
      INTEGER, PARAMETER :: MAX_ITERATIONS = 9999

      !=================================================================
      ! SMVGEAR begins here!
      !=================================================================
      COUNTER   = 0
      ICOUNT    = 0
      NSUBFUN   = 0
      NPDERIV   = 0
      NSTEPS    = 0
      IFAIL     = 0
      JFAIL     = 0
      LFAIL     = 0
      NFAIL     = 0
      NYLOWDEC  = 0
      TINTERVAL = TIMEINTV(NCS)
      ISCHAN    = ISCHANG( NCS)
      ISCHAN1   = ISCHAN - 1

      !=================================================================
      ! Added for the ND65 prod/loss diagnostic, in order to prevent
      ! ND65 prod/loss families from being counted towards the 
      ! SMVGEAR convergence criteria. (ljm, bmy, 5/9/03)
      NNOFAM    = ISCHAN - NFAMILIES
      ORDER     = REAL( NNOFAM )
      !=================================================================
C
      IABOVE    = ORDER * 0.4d0
C
      DO 115 KLOOP   = 1, KTLOOP
       IABOVK(KLOOP) = IABOVE
 115  CONTINUE
C
      HRMAX     = 0.3d0
      HMAX      = HMAXUSE( NCSP)
      YFAC      = 1.0d0
      ERRINIT   = MIN(ERRMAX(NCS),1.0D-03)
C
C *********************************************************************
C            START TIME INTERVAL OR RE-ENTER AFTER TOTAL FAILURE
C *********************************************************************
C
 120  IDOUB     = 2
      NSLP      = MBETWEEN
      JRESTAR   = 0 
      DELT      = 0.d0
      XELAPS    = 0.d0
      XELAPLAST = -1.d0 
      TOLD      = 0.d0
      TIMREMAIN = TINTERVAL
      RELTOL1   = YFAC   / ERRINIT
      RELTOL2   = YFAC   / ERRMAX(NCS)
      RELTOL3   = 1.d0   / ERRMAX(NCS)
      ABTOLER1  = ABTOL(6,NCS) * RELTOL1
      ABTOLER2  = ABTOL(6,NCS) * RELTOL2
C
C *********************************************************************
C                  INITIALIZE CONCENTRATION ARRAY 
C *********************************************************************
C CORIG = ORIGINAL CONCENTRATIONS, WHICH DO NOT CHANGE IN SMVGEAR
C CNEW  = FINAL CONCENTRATIONS, CALCULATED IN SMVGEAR
C

      DO 129 JNEW         = 1, ISCHAN
       DO 127 KLOOP       = 1, KTLOOP
        CNEW( KLOOP,JNEW) = CORIG(KLOOP,JNEW)
 127   CONTINUE
 129  CONTINUE

C
C *********************************************************************
C  RE-ENTER HERE IF TOTAL FAILURE OR IF RESTARTING WITH NEW CELL BLOCK 
C *********************************************************************
C
 140  HRATIO    = 0.d0
      ASN1      = 1.d0
      IFSUCCESS = 1
      RDELMAX   = 1.0d+04

C *********************************************************************
C                         INITIALIZE PHOTRATES 
C *********************************************************************
C
      ! Called for photorates with no active loss terms (bdf, 4/18/03)
      IF (IFSUN.EQ.1) CALL UPDATE
C
C *********************************************************************
C               INITIALIZE FIRST DERIVATIVE FOR CHEMISTRY 
C *********************************************************************
C
      CALL SUBFUN
C
C *********************************************************************
C                DETERMINE INITIAL ABSOLUTE ERROR TOLERANCE 
C *********************************************************************
C IABOVK  = NUMBER OF SPECIES WHOSE CONCENTRATIONS ARE LARGER THAN YABST
C ISREORD = 1: CALC INITIAL STIFFNESS BEFORE RUNNING CODE TO REORDER CELLS
C              IN THIS CASE, USE PHOTORATES FOR END OF TIME-INTERVAL
C         = 2: DO NORMAL CALCULATIONS
C KGRP    = COUNTS NUMBER OF CONCENTRATIONS ABOVE ABTOL(I), I = 1..  
C YABST   = ABSOLUTE ERROR TOLERANCE (MOLEC. CM-3 FOR GASES) 
C ABTOL   = PRE-DEFINED ABSOLUTE ERROR TOLERANCES 
C
      DO 142 KLOOP    = 1, KTLOOP
       ERRHOLD(KLOOP) = 0.d0
 142  CONTINUE 
C
      IF (ISREORD.NE.1) THEN
C
       DO 134 K              = 1, 5        
        DO 132 KLOOP         = 1, KTLOOP
         KGRP(KLOOP,K)       = 0
 132    CONTINUE
 134   CONTINUE
C
       DO 136 JSPC           = 1, ISCHAN
        !==============================================================
        ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
        ! prod/loss species from being counted towards the convergence
        ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)         
        IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN
           DO 135 KLOOP         = 1, KTLOOP
              CNW                 = CNEW(KLOOP,JSPC)
              IF (CNW.GT.ABTOL(1,NCS)) THEN
                 KGRP(KLOOP,1)      = KGRP(KLOOP,1) + 1
              ELSEIF (CNW.GT.ABTOL(2,NCS)) THEN
                 KGRP(KLOOP,2)      = KGRP(KLOOP,2) + 1
              ELSEIF (CNW.GT.ABTOL(3,NCS)) THEN
                 KGRP(KLOOP,3)      = KGRP(KLOOP,3) + 1
              ELSEIF (CNW.GT.ABTOL(4,NCS)) THEN
                 KGRP(KLOOP,4)      = KGRP(KLOOP,4) + 1
              ELSEIF (CNW.GT.ABTOL(5,NCS)) THEN
                 KGRP(KLOOP,5)      = KGRP(KLOOP,5) + 1
              ENDIF
 135       CONTINUE
        ENDIF
        !==============================================================
 136   CONTINUE
C

       DO 137 KLOOP         = 1, KTLOOP
        K1                  = KGRP(KLOOP,1)
        K2                  = KGRP(KLOOP,2) + K1
        K3                  = KGRP(KLOOP,3) + K2
        K4                  = KGRP(KLOOP,4) + K3
        K5                  = KGRP(KLOOP,5) + K4
        IF (K1.GT.IABOVK(KLOOP)) THEN
         YABST(KLOOP)       = ABTOL(1,NCS)
        ELSEIF (K2.GT.IABOVK(KLOOP)) THEN
         YABST(KLOOP)       = ABTOL(2,NCS)
        ELSEIF (K3.GT.IABOVK(KLOOP)) THEN
         YABST(KLOOP)       = ABTOL(3,NCS)
        ELSEIF (K4.GT.IABOVK(KLOOP)) THEN
         YABST(KLOOP)       = ABTOL(4,NCS)
        ELSEIF (K5.GT.IABOVK(KLOOP)) THEN
         YABST(KLOOP)       = ABTOL(5,NCS)
        ELSE
         YABST(KLOOP)       = ABTOL(6,NCS)
        ENDIF
 137   CONTINUE
C
       DO 139 JSPC      = 1, ISCHAN
         !=============================================================
         ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
         ! prod/loss species from being counted towards the convergence
         ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)
         IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN 
            DO 138 KLOOP    = 1, KTLOOP
               CNEWYLOW       = CNEW(KLOOP,JSPC) + YABST(KLOOP) *RELTOL1
               ERRYMAX        = GLOSS(KLOOP,JSPC) / CNEWYLOW
               ERRHOLD(KLOOP) = ERRHOLD(KLOOP) + ERRYMAX * ERRYMAX
 138        CONTINUE
         ENDIF
         !=============================================================
 139   CONTINUE

C
      ELSE
C
C *********************************************************************
C          USE LOWEST ABSOLUTE ERROR TOLERANCE WHEN REORDERING 
C          IF REORDERING, SET ERRMX2 THEN RETURN TO PHYSPROC.F 
C *********************************************************************
C ABTOLER1 = YFAC * ABTOL(6,NCS) / MIN(ERRMAX,1.0E-03) 
C 
       DO 144 JSPC      = 1, ISCHAN  
         !=============================================================
         ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
         ! prod/loss species from being counted towards the convergence
         ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)
         IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN
            DO 143 KLOOP    = 1, KTLOOP 
               ERRYMAX        = GLOSS(KLOOP,JSPC)/
     &                          (CNEW(KLOOP,JSPC)+ABTOLER1)
               ERRHOLD(KLOOP) = ERRHOLD(KLOOP) + ERRYMAX * ERRYMAX
 143        CONTINUE
         ENDIF
         !==============================================================
 144   CONTINUE
C
       IF (ISREORD.EQ.1) THEN 
        DO 150 KLOOP           = 1, KTLOOP 
         ERRMX2(JLOOPLO+KLOOP) = ERRHOLD(KLOOP)
 150    CONTINUE
C
        RETURN
       ENDIF
      ENDIF

C
C *********************************************************************
C               CALCULATE INITIAL TIME STEP SIZE (S) 
C *********************************************************************
C SQRT(ERRHOLD / [ERRINIT * ORDER]) = RMSNORM OF ERROR SCALED TO ERRINIT 
C                                * CNEW + ABTOL/RELTOL
C
      RMSTOP         = 0.d0
C
      DO 151 KLOOP   = 1, KTLOOP
       IF (ERRHOLD(KLOOP).GT.RMSTOP) RMSTOP = ERRHOLD(KLOOP)   
 151  CONTINUE
C
      DELT1        = SQRT(ERRINIT / (ABST2(NCS) + RMSTOP / ORDER)) 
      DELT         = MAX(MIN(DELT1,TIMREMAIN,HMAX),HMIN) 
C
C *********************************************************************
C                      SET INITIAL ORDER TO 1
C *********************************************************************
C
      NQQOLD       = 0
      NQQ          = 1 
      JEVAL        = 1
      RDELT        = 1.0d0
C
C *********************************************************************
C *   STORE INITIAL CONCENTRATION AND FIRST DERIVATIVES x TIME-STEP   * 
C *********************************************************************
C
      DO 155 JSPC        = 1, ISCHAN 
       JS1               = ISCHAN + JSPC
       DO 154 KLOOP      = 1, KTLOOP
        CONC(KLOOP,JSPC) = CNEW(KLOOP,JSPC)
        CONC(KLOOP,JS1)  = DELT * GLOSS(KLOOP,JSPC) 
 154   CONTINUE
 155  CONTINUE

C
C *********************************************************************
C ** UPDATE COEFFICIENTS OF THE ORDER. NQQ IS THE ORDER. ASET AND    **
C ** PERTS2 ARE DEFINED IN SUBROUTINE KSPARSE. NOTE THAT PERTS2      **
C ** IS THE ORIGINAL PERTST**2                                       **
C *********************************************************************
C
 170  IF (NQQ.NE.NQQOLD) THEN
       NQQOLD                = NQQ 
       KSTEP                 = NQQ + 1
       HRATIO                = HRATIO * ASET(NQQ,1) / ASN1
       ASN1                  = ASET(NQQ,1)  
       ENQQ                  = PERTS2(NQQ,1) * ORDER
       EUP                   = PERTS2(NQQ,2) * ORDER 
       EDWN                  = PERTS2(NQQ,3) * ORDER  
       CONP3                 = 1.4d0 / ( EUP**ENQQ3(NQQ))
       CONP2                 = 1.2d0 / (ENQQ**ENQQ2(NQQ))
       CONP1                 = 1.3d0 / (EDWN**ENQQ1(NQQ))
       NQQISC                = NQQ * ISCHAN 
      ENDIF
      counter=counter+1
C
C *********************************************************************
C   LIMIT SIZE OF RDELT, THEN RECALCULATE NEW TIME STEP AND UPDATE
C HRATIO. USE HRATIO TO DETERMINE WHETHER PDERIV SHOULD BE CALLED AGAIN 
C *********************************************************************
C
      HMTIM         = MIN(HMAX,TIMREMAIN) 
      RDELT         = MIN(RDELT,RDELMAX,HMTIM/DELT)
      DELT          = DELT   * RDELT 
      HRATIO        = HRATIO * RDELT 
      XELAPS        = XELAPS + DELT  

C
      IF (ABS(HRATIO-1.0).GT.HRMAX.OR.NSTEPS.GE.NSLP) JEVAL = 1
C
C *********************************************************************
C      IF TIME STEP < HMIN, TIGHTEN ABSOLOUTE ERROR TOLERANCE AND 
C          RESTART INTEGRATION AT BEGINNING OF TIME INTERVAL
C *********************************************************************
C
      IF (DELT.LT.HMIN) THEN
       WRITE(6,233)DELT,KBLK,KTLOOP,NCS,TIME,TIMREMAIN,YFAC,ERRMAX(NCS)
       NYLOWDEC  = NYLOWDEC + 1
       YFAC      = YFAC  * 0.01d0
C
       IF (NYLOWDEC.EQ.10) THEN
        LLOOPA      = 1
        LLOOPB      = KTLOOP
        WRITE(6,234)
C
        DO 177 KLOOP = 1, KTLOOP
         JLOOP       = JREORDER(JLOOPLO+KLOOP)
         K           = (JLOOP - 1) / NLOOP + 1
         MLOOP       = JLOOP - (K - 1) * NLOOP
         M1          = (MLOOP - 1) / NLONG + 1
         M2          = MLOOP - (M1 - 1) * NLONG
         WRITE(6,685) M1, M2, K, ERRHOLD(KLOOP)
 177    CONTINUE
C
        DO 178 JNEW = 1, ISCHAN
         JOLD       = INEWOLD(JNEW,NCS)
         WRITE(6,690) JNEW, NCS, NAMENCS(JOLD,NCS),CORIG(LLOOPA,JNEW),
     1                CORIG(LLOOPB,JNEW)
 178    CONTINUE

        ! Stop run w/ error msg
        CALL GEOS_CHEM_STOP
       ENDIF
C
       GOTO 120
      ENDIF

C
C *********************************************************************
C * IF THE DELT IS DIFFERENT THAN DURING THE LAST STEP (IF RDELT NE   *
C * 1), THEN SCALE THE DERIVATIVES                                    *
C *********************************************************************
C
      IF (RDELT.NE.1.0) THEN
       RDELTA            = 1.0d0
       I1                = 1
       DO 184 J          = 2, KSTEP 
        RDELTA           = RDELTA * RDELT 
        I1               = I1 + ISCHAN 
        DO 182 I         = I1, I1 + ISCHAN1
         DO 180 KLOOP    = 1, KTLOOP
          CONC(KLOOP,I)  = CONC(KLOOP,I) * RDELTA  
 180     CONTINUE
 182    CONTINUE
 184   CONTINUE
      ENDIF
C
C *********************************************************************
C * UPDATE PHOTO RATES BECAUSE THE TIME CHANGED.                      *
C * NOTE THAT A TIME CHANGE COULD CORRESPOND TO EITHER A SUCCESSFUL   *
C * OR FAILED STEP                                                    * 
C *********************************************************************
C
      ! Called for photorates with no active loss terms (bdf, 4/18/03)
      IF (IFSUN.EQ.1.AND.XELAPS.NE.XELAPLAST) CALL UPDATE
C
C *********************************************************************
C * IF THE LAST STEP WAS SUCCESSFUL, RESET RDELMAX = 10 AND UPDATE    *
C * THE CHOLD ARRAY WITH CURRENT VALUES OF CNEW.                      * 
C *********************************************************************
C
      IF (IFSUCCESS.EQ.1) THEN
       RDELMAX             = 10.d0
C
C *********************************************************************
C                DETERMINE NEW ABSOLUTE ERROR TOLERANCE 
C *********************************************************************
C KGRP    = COUNTS NUMBER OF CONCENTRATIONS ABOVE ABTOL(I), I = 1..  
C YABST   = ABSOLUTE ERROR TOLERANCE (MOLEC. CM-3 FOR GASES) 
C ABTOL   = PRE-DEFINED ABSOLUTE ERROR TOLERANCES 
C
       IF (MOD(NSTEPS,3).EQ.2) THEN
        DO 203 K              = 1, 5        
         DO 201 KLOOP         = 1, KTLOOP
          KGRP(KLOOP,K)       = 0
 201     CONTINUE
 203    CONTINUE
C
        DO 207 JSPC           = 1, ISCHAN
         !==============================================================
         ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
         ! prod/loss species from being counted towards the convergence
         ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)
         IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN 
            DO 205 KLOOP         = 1, KTLOOP
               CNW                 = CNEW(KLOOP,JSPC)
               IF (CNW.GT.ABTOL(1,NCS)) THEN
                  KGRP(KLOOP,1)      = KGRP(KLOOP,1) + 1
               ELSEIF (CNW.GT.ABTOL(2,NCS)) THEN
                  KGRP(KLOOP,2)      = KGRP(KLOOP,2) + 1
               ELSEIF (CNW.GT.ABTOL(3,NCS)) THEN
                  KGRP(KLOOP,3)      = KGRP(KLOOP,3) + 1
               ELSEIF (CNW.GT.ABTOL(4,NCS)) THEN
                  KGRP(KLOOP,4)      = KGRP(KLOOP,4) + 1
               ELSEIF (CNW.GT.ABTOL(5,NCS)) THEN
                  KGRP(KLOOP,5)      = KGRP(KLOOP,5) + 1
               ENDIF
 205        CONTINUE
         ENDIF
         !==============================================================
 207    CONTINUE
C
        DO 209 KLOOP         = 1, KTLOOP
         K1                  = KGRP(KLOOP,1)
         K2                  = KGRP(KLOOP,2) + K1
         K3                  = KGRP(KLOOP,3) + K2 
         K4                  = KGRP(KLOOP,4) + K3 
         K5                  = KGRP(KLOOP,5) + K4 
         IF (K1.GT.IABOVK(KLOOP)) THEN
          YABST(KLOOP)       = ABTOL(1,NCS)
         ELSEIF (K2.GT.IABOVK(KLOOP)) THEN
          YABST(KLOOP)       = ABTOL(2,NCS)
         ELSEIF (K3.GT.IABOVK(KLOOP)) THEN
          YABST(KLOOP)       = ABTOL(3,NCS)
         ELSEIF (K4.GT.IABOVK(KLOOP)) THEN
          YABST(KLOOP)       = ABTOL(4,NCS)
         ELSEIF (K5.GT.IABOVK(KLOOP)) THEN
          YABST(KLOOP)       = ABTOL(5,NCS)
         ELSE
          YABST(KLOOP)       = ABTOL(6,NCS)
         ENDIF
 209    CONTINUE
       ENDIF
C
       DO 213 JSPC         = 1, ISCHAN
        !===============================================================
        ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
        ! prod/loss species from being counted towards the convergence
        ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)
        IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN
           DO 211 KLOOP         = 1, KTLOOP
             CHOLD(KLOOP,JSPC)  = RELTOL3 / (MAX(CNEW(KLOOP,JSPC),0.D0)
     1                          + YABST(KLOOP) * RELTOL2)
 211      CONTINUE 
        ENDIF
        !===============================================================
 213   CONTINUE 
C
      ENDIF 

C     ENDIF IFSUCCESS.EQ.1
C
C *********************************************************************
C * COMPUTE THE PREDICTED CONCENTRATION AND DERIVATIVES BY MULTIPLY-  *
C * ING PREVIOUS VALUES BY THE PASCAL TRIANGLE MATRIX.                * 
C *********************************************************************
C THIS SET OF OPERATIONS IS EQUIVALENT TO THE REVERSE OF LOOP 419.
C THE EXPANSION OF THE PASCAL TRIANGLE MATRIX WAS CALCULATED BY B. SCHWARTZ.
C THE FIRST DERIVATIVE MULTIPLIED BY THE TIME STEP IS THE SUM 
C OF TERMS ADDED TO CONC(KLOOP,I)
C
      IF (NQQ.EQ.1) THEN
       DO 236 I          = 1, ISCHAN
        J                = I + ISCHAN
        DO 235 KLOOP     = 1, KTLOOP
         CONC(KLOOP,I)   = CONC(KLOOP,I) + CONC(KLOOP,J)
 235    CONTINUE
 236   CONTINUE
C 
      ELSEIF (NQQ.EQ.2) THEN
C
       DO 238 I         = 1, ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN
        DO 237 KLOOP    = 1, KTLOOP
         CONC(KLOOP, I) = CONC(KLOOP, I) + CONC(KLOOP,J1)
     1                  +                  CONC(KLOOP,J2)
         CONC(KLOOP,J1) = CONC(KLOOP,J1) + CONC(KLOOP,J2) * 2.d0
 237    CONTINUE
 238   CONTINUE
C
      ELSEIF (NQQ.EQ.3) THEN
C
       DO 240 I         = 1,   ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN
        J3              = J2 + ISCHAN
        DO 239 KLOOP    = 1, KTLOOP
         CONC3J3        = CONC(KLOOP,J3) * 3.d0
         CONC(KLOOP, I) = CONC(KLOOP, I) + CONC(KLOOP,J1) 
     1                  + CONC(KLOOP,J2) + CONC(KLOOP,J3)
         CONC(KLOOP,J1) = CONC(KLOOP,J1) + CONC(KLOOP,J2)*2.d0 + CONC3J3 
         CONC(KLOOP,J2) = CONC(KLOOP,J2) + CONC3J3
 239    CONTINUE
 240   CONTINUE
C
      ELSEIF (NQQ.EQ.4) THEN
C
       DO 242 I         = 1, ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN
        J3              = J2 + ISCHAN
        J4              = J3 + ISCHAN
        DO 241 KLOOP    = 1, KTLOOP
         CONC3J3        = CONC(KLOOP,J3) * 3.d0 
         CONC4J4        = CONC(KLOOP,J4) * 4.d0
         CONC(KLOOP, I) = CONC(KLOOP, I) + CONC(KLOOP,J1) 
     1                  + CONC(KLOOP,J2) + CONC(KLOOP,J3) 
     2                  + CONC(KLOOP,J4)
         CONC(KLOOP,J1) = CONC(KLOOP,J1) + CONC(KLOOP,J2)*2.d0 + CONC3J3 
     1                                   + CONC4J4
         CONC(KLOOP,J2) = CONC(KLOOP,J2) + CONC3J3 + CONC(KLOOP,J4)*6.d0 
         CONC(KLOOP,J3) = CONC(KLOOP,J3) + CONC4J4
 241    CONTINUE 
 242   CONTINUE
C
      ELSEIF (NQQ.EQ.5) THEN
C
       DO 244 I         = 1, ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN
        J3              = J2 + ISCHAN
        J4              = J3 + ISCHAN
        J5              = J4 + ISCHAN
        DO 243 KLOOP    = 1, KTLOOP
         CONC3J3        = CONC(KLOOP,J3) * 3.d0 
         CONC4J4        = CONC(KLOOP,J4) * 4.d0 
         CONC5J5        = CONC(KLOOP,J5) * 5.d0 
         CONC10J5       = CONC5J5 + CONC5J5 
         CONC(KLOOP, I) = CONC(KLOOP, I) + CONC(KLOOP,J1) 
     1                  + CONC(KLOOP,J2) + CONC(KLOOP,J3) 
     2                  + CONC(KLOOP,J4) + CONC(KLOOP,J5)
         CONC(KLOOP,J1) = CONC(KLOOP,J1) + CONC(KLOOP,J2)*2.d0 + CONC3J3 
     1                                   + CONC4J4 + CONC5J5
         CONC(KLOOP,J2) = CONC(KLOOP,J2) + CONC3J3 + CONC(KLOOP,J4)*6.d0 
     1                                   + CONC10J5
         CONC(KLOOP,J3) = CONC(KLOOP,J3) + CONC4J4 + CONC10J5
         CONC(KLOOP,J4) = CONC(KLOOP,J4) + CONC5J5
 243    CONTINUE
 244   CONTINUE
      ENDIF

C
C *********************************************************************
C ************************** CORRECTION LOOP **************************
C * TAKE UP TO 3 CORRECTOR ITERATIONS. TEST CONVERGENCE BY REQUIRING  *
C * THAT CHANGES BE LESS THAN THE RMS NORM WEIGHTED BY CHOLD.         * 
C * ACCUMULATE THE CORRECTION IN THE ARRAY DTLOS(). IT EQUALS THE     *
C * THE J-TH DERIVATIVE OF CONC() MULTIPLIED BY DELT**KSTEP /         *
C * (FACTORIAL(KSTEP-1)*ASET(KSTEP)); THUS, IT IS PROPORTIONAL TO THE *
C * ACTUAL ERRORS TO THE LOWEST POWER OF DELT PRESENT (DELT**KSTEP)   *
C *********************************************************************
C
 220  L3                  = 0
      DO 232 JSPC         = 1, ISCHAN
       DO 230 KLOOP       = 1, KTLOOP
        CNEW(KLOOP,JSPC)  = CONC(KLOOP,JSPC)
        DTLOS(KLOOP,JSPC) = 0.d0
 230   CONTINUE
 232  CONTINUE

C
C *********************************************************************
C * IF JEVAL = 1, RE-EVALUATE PREDICTOR MATRIX P = I - H * ASET(1) *J *
C * BEFORE STARTING THE CORRECTOR ITERATION. AFTER CALLING PDERIV,    * 
C * SET JEVAL = -1 TO PREVENT RECALLING PDERIV UNLESS NECESSARY LATER.*
C * CALL DECOMP TO DECOMPOSE THE MATRIX                               *
C *********************************************************************
C
      IF (JEVAL.EQ.1) THEN
       R1DELT   = -ASN1 * DELT
C
       CALL PDERIV
C     
       CALL DECOMP 
       JEVAL    = -1 
       HRATIO   = 1.0d0
       NSLP     = NSTEPS + MBETWEEN
       DRATE    = 0.7d0
      ENDIF

C
C *********************************************************************
C *   EVALUATE THE FIRST DERIVATIVE USING CORRECTED VALUES OF CNEW    * 
C *********************************************************************
C     
 270  CALL SUBFUN
C
C *********************************************************************
C * IN THE CASE OF THE CHORD METHOD, COMPUTE ERROR (GLOSS) FROM THE   * 
C * CORRECTED CALCULATION OF THE FIRST DERIVATIVE                     * 
C *********************************************************************
C
      DO 362 JSPC         = 1,  ISCHAN
       J                  = JSPC + ISCHAN
       DO 360 KLOOP       = 1, KTLOOP
        GLOSS(KLOOP,JSPC) = DELT * GLOSS(KLOOP,JSPC)  
     1                    - (CONC(KLOOP,J) + DTLOS(KLOOP,JSPC))
 360   CONTINUE
 362  CONTINUE
C
C *********************************************************************
C * SOLVE THE LINEAR SYSTEM OF EQUATIONS WITH THE CORRECTOR ERROR.    *
C * BACKSUB.F SOLVES BACKSUBSTITUTION OVER MATRIX OF PARTIAL DERIVS.  * 
C *********************************************************************
C
      CALL BACKSUB
C
C *********************************************************************
C * SUM-UP THE ACCUMULATED ERROR, CORRECT THE CONCENTRATION WITH THE  * 
C * ERROR, AND BEGIN TO CALCULATE THE RMSNORM OF THE ERROR RELATIVE   *  
C * TO CHOLD.                                                         *
C *********************************************************************
C
      DO 365 KLOOP     = 1, KTLOOP
       DELY(KLOOP)     = 0.d0
 365  CONTINUE
C
      specmax = 0.0d0
      IF (ASN1.EQ.1.0) THEN 
       DO 368 I         = 1, ISCHAN
           DO 366 KLOOP    = 1, KTLOOP
              DTLOS(KLOOP,I) = DTLOS(KLOOP,I)  + GLOSS(KLOOP,I)
              CNEW(KLOOP,I)  = CONC(KLOOP,I)   + DTLOS(KLOOP,I)
              !=========================================================
              ! Added for the ND65 prod/loss diagnostic.  This prevents 
              ! ND65 prod/loss species from being counted towards the 
              ! convergence criteria for SMVGEAR (ljm, bmy, 5/9/03)
              IF ( ITS_NOT_A_ND65_FAMILY(I) ) THEN
                 ERRYMAX        = GLOSS(KLOOP,I)  * CHOLD(KLOOP,I)
                 DELY(KLOOP)    = DELY(KLOOP)     + ERRYMAX * ERRYMAX
              ENDIF
              !=========================================================
 366       CONTINUE
 368   CONTINUE
      ELSE
       DO 372 I         = 1, ISCHAN
           DO 370 KLOOP    = 1, KTLOOP
              DTLOS(KLOOP,I) = DTLOS(KLOOP,I)  + GLOSS(KLOOP,I)
              CNEW(KLOOP,I)  = CONC(KLOOP,I)   + ASN1  *  DTLOS(KLOOP,I)
              !=========================================================
              ! Added for the ND65 prod/loss diagnostic.  This prevents 
              ! ND65 prod/loss species from being counted towards the 
              ! convergence criteria for SMVGEAR (ljm, bmy, 5/9/03)
              IF ( ITS_NOT_A_ND65_FAMILY(I) ) THEN
                 ERRYMAX        = GLOSS(KLOOP,I)  * CHOLD(KLOOP,I)
                 DELY(KLOOP)    = DELY(KLOOP)     + ERRYMAX * ERRYMAX
              ENDIF
              !=========================================================
 370       CONTINUE
 372   CONTINUE
      ENDIF

C
C *********************************************************************
C * SET THE PREVIOUS RMS ERROR AND CALCULATE THE NEW RMS ERROR.       *
C * IF DCON < 1, THEN SUFFICIENT CONVERGENCE HAS OCCURRED. OTHERWISE, *
C * IF THE RATIO OF THE CURRENT TO PREVIOUS RMSERR IS DECREASING,     *
C * ITERATE MORE. IF IT IS NOT, THEN THE CONVERGENCE TEST FAILED      *
C *********************************************************************
C
      RMSERRP           = RMSERR
      DER2MAX           = 0.d0
C
      ksave=0d0
      DO 427 KLOOP      = 1, KTLOOP
       IF (DELY(KLOOP).GT.DER2MAX) DER2MAX = DELY(KLOOP)   
       ijsave=jlooplo+kloop
       ksave=kloop
 427  CONTINUE
C
      RMSERR             = SQRT(DER2MAX / ORDER)
C
      L3                 = L3 + 1
C
      IF (L3.GT.1) THEN
       RMSRAT            = RMSERR / RMSERRP
       DRATE             = MAX(0.2d0 * DRATE, RMSRAT)
      ENDIF
C
      DCON               = RMSERR * MIN(CONPST(NQQ),CONP15(NQQ)*DRATE)
C
C *********************************************************************
C       IF CONVERGENCE OCCURS, GO ON TO CHECK ACCUMULATED ERROR
C *********************************************************************
C
      IF (DCON .LE. 1.0) THEN
       GOTO 390 
C
C *********************************************************************
C  IF NONCONVERGENCE AFTER ONE STEP, RE-EVALUATE FIRST DERIVATIVE WITH
C                           NEW VALUES OF CNEW
C *********************************************************************
C
C     ELSEIF (L3.LT.MSTEP.AND.(L3.EQ.1.OR.RMSRAT.LE.0.9)) THEN
      ELSEIF (L3.EQ.1) THEN
       GOTO 270
C
C *********************************************************************
C *             THE CORRECTOR ITERATION FAILED TO CONVERGE            *
C * IF THE JACOBIAN MATRIX IS MORE THAN ONE STEP OLD, UPDATE THE      *
C * JACOBIAN AND TRY CONVERGENCE AGAIN. IF THE JACOBIAN IS CURRENT,   * 
C * THEN REDUCE THE TIME-STEP, RE-SET THE ACCUMULATED DERIVATIVES TO  *
C * THEIR VALUES BEFORE THE FAILED STEP, AND RETRY WITH THE SMALLER   *
C * STEP.                                                             *
C *********************************************************************
C

      ELSEIF (JEVAL .EQ. 0) THEN
       IFAIL           = IFAIL + 1
       JEVAL           = 1
       GOTO 220
      ENDIF
C
      NFAIL            = NFAIL + 1
      RDELMAX          = 2.0d0
      JEVAL            = 1
      IFSUCCESS        = 0
      XELAPS           = TOLD
      RDELT            = FRACDEC
C
C *********************************************************************
C               SUBTRACT OFF DERIVATIVES PREVIOUSLY ADDED 
C *********************************************************************
C THIS SET OF OPERATIONS IS EQUIVALENT TO LOOP 419.
C
      IF (NQQ.EQ.1) THEN
       DO 376 I        = 1, ISCHAN
        J              = I + ISCHAN
        DO 375 KLOOP   = 1, KTLOOP
         CONC(KLOOP,I) = CONC(KLOOP,I) - CONC(KLOOP,J)
 375    CONTINUE
 376   CONTINUE
C
      ELSEIF (NQQ.EQ.2) THEN
C
       DO 378 I         = 1, ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN 
        DO 377 KLOOP    = 1, KTLOOP
         CONC(KLOOP, I) = CONC(KLOOP, I) - CONC(KLOOP,J1)-CONC(KLOOP,J2) 
         CONC(KLOOP,J1) = CONC(KLOOP,J1) - CONC(KLOOP,J2) * 2.d0 
 377    CONTINUE
 378   CONTINUE
C
      ELSEIF (NQQ.EQ.3) THEN
       DO 380 I         = 1,   ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN
        J3              = J2 + ISCHAN
        DO 379 KLOOP    = 1, KTLOOP
         CONC3J3        = CONC(KLOOP,J3) * 3.d0 
         CONC(KLOOP, I) = CONC(KLOOP, I) - CONC(KLOOP,J1) 
     1                  - CONC(KLOOP,J2) - CONC(KLOOP,J3)
         CONC(KLOOP,J1) = CONC(KLOOP,J1) - CONC(KLOOP,J2)*2.d0 - CONC3J3 
         CONC(KLOOP,J2) = CONC(KLOOP,J2) - CONC3J3
 379    CONTINUE
 380   CONTINUE
C
      ELSEIF (NQQ.EQ.4) THEN
C
       DO 382 I         = 1, ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN
        J3              = J2 + ISCHAN
        J4              = J3 + ISCHAN
        DO 381 KLOOP    = 1, KTLOOP
         CONC3J3        = CONC(KLOOP,J3) * 3.d0 
         CONC4J4        = CONC(KLOOP,J4) * 4.d0 
         CONC(KLOOP, I) = CONC(KLOOP, I) - CONC(KLOOP,J1) 
     1                  - CONC(KLOOP,J2) - CONC(KLOOP,J3)
     2                  - CONC(KLOOP,J4)
         CONC(KLOOP,J1) = CONC(KLOOP,J1) - CONC(KLOOP,J2)*2.d0 - CONC3J3 
     1                                   - CONC4J4
         CONC(KLOOP,J2) = CONC(KLOOP,J2) - CONC3J3 - CONC(KLOOP,J4)*6.d0 
         CONC(KLOOP,J3) = CONC(KLOOP,J3) - CONC4J4
 381    CONTINUE 
 382   CONTINUE
C
      ELSEIF (NQQ.EQ.5) THEN
C
       DO 384 I         = 1, ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN
        J3              = J2 + ISCHAN
        J4              = J3 + ISCHAN
        J5              = J4 + ISCHAN
        DO 383 KLOOP    = 1, KTLOOP
         CONC3J3        = CONC(KLOOP,J3) * 3.d0 
         CONC4J4        = CONC(KLOOP,J4) * 4.d0
         CONC5J5        = CONC(KLOOP,J5) * 5.d0
         CONC10J5       = CONC5J5 + CONC5J5 
         CONC(KLOOP, I) = CONC(KLOOP, I) - CONC(KLOOP,J1) - 
     1                    CONC(KLOOP,J2) - CONC(KLOOP,J3) - 
     2                    CONC(KLOOP,J4) - CONC(KLOOP,J5)
         CONC(KLOOP,J1) = CONC(KLOOP,J1) - CONC(KLOOP,J2)*2.d0 - CONC3J3 
     1                                   - CONC4J4 - CONC5J5
         CONC(KLOOP,J2) = CONC(KLOOP,J2) - CONC3J3 - CONC(KLOOP,J4)*6.d0 
     2                                   - CONC10J5
         CONC(KLOOP,J3) = CONC(KLOOP,J3) - CONC4J4 - CONC10J5
         CONC(KLOOP,J4) = CONC(KLOOP,J4) - CONC5J5
 383    CONTINUE
 384   CONTINUE
      ENDIF
C
      GOTO 170
C
C *********************************************************************
C *               THE CORRECTOR ITERATION CONVERGED                   *  
C * SET JEVAL = 0 SO THAT IT DOES NOT NEED TO BE CALLED THE NEXT STEP *
C * IF ALL ELSE GOES WELL. NEXT, TEST THE ACCUMULATED ERROR FROM THE  * 
C *                    CONVERGENCE PROCESS, ABOVE                     * 
C *********************************************************************
C
 390  JEVAL          = 0
C
      IF (L3.GT.1) THEN
       DO 395 KLOOP  = 1, KTLOOP
        DELY(KLOOP)  = 0.d0
 395   CONTINUE
C
       DO 402 JSPC   = 1, ISCHAN     
        !===============================================================
        ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
        ! prod/loss species from being counted towards the convergence
        ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)
        IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN
           DO 400 KLOOP = 1, KTLOOP 
              ERRYMAX     = DTLOS(KLOOP,JSPC) * CHOLD(KLOOP,JSPC)
              DELY(KLOOP) = DELY(KLOOP) + ERRYMAX * ERRYMAX
 400       CONTINUE
        ENDIF
        !===============================================================
 402   CONTINUE
C
       DER2MAX       = 0.d0
C
       DO 405 KLOOP  = 1, KTLOOP
        IF (DELY(KLOOP).GT.DER2MAX) DER2MAX = DELY(KLOOP)   
 405   CONTINUE 
      ENDIF
C
C *********************************************************************
C *                 THE ACCUMULATED ERROR TEST FAILED                 *  
C * IN ALL CASES, RE-SET THE DERIVATIVES TO THEIR VALUES BEFORE THE   *
C * LAST TIME-STEP. NEXT                                              *    
C *                                                                   * 
C * (A)   RE-ESTIMATE A TIME-STEP AT THE SAME OR ONE LOWER ORDER AND  *   
C *        RETRY THE STEP.                                            *  
C * (B)   IF THE FIRST ATTEMPTS FAIL, RETRY THE STEP AT FRACDEC x     * 
C *        THE PRIOR STEP                                             *  
C * (C)   IF THIS FAILS, RE-SET THE ORDER TO 1 AND GO BACK TO THE     * 
C *        BEGINNING, AT ORDER = 1, BECAUSE ERRORS OF THE WRONG ORDER *
C *        HAVE ACCUMULATED                                           *
C *********************************************************************
C 
      IF (DER2MAX.GT.ENQQ) THEN
       XELAPS           = TOLD
       LFAIL            = LFAIL + 1
       JFAIL            = JFAIL  + 1
C
       I1               = NQQISC + 1
       DO 419 JB        = 1, NQQ
        I1              = I1 - ISCHAN
        DO 417 I        = I1, NQQISC
         J              = I + ISCHAN
         DO 415 KLOOP   = 1, KTLOOP
          CONC(KLOOP,I) = CONC(KLOOP,I) - CONC(KLOOP,J)
 415     CONTINUE
 417    CONTINUE
 419   CONTINUE
C
       RDELMAX          = 2.0d0
       IF (JFAIL.LE.6) THEN
        IFSUCCESS       = 0 
        RDELTUP         = 0.0d0
        GOTO 540
       ELSEIF (JFAIL.LE.20) THEN
        IFSUCCESS       = 0 
        RDELT           = FRACDEC
        GOTO 170 
       ELSE
        DELT            = DELT * 0.1d0
        RDELT           = 1.
        JFAIL           = 0
        JRESTAR         = JRESTAR + 1
        IDOUB           = 5
C
        DO 432 JSPC        = 1, ISCHAN
         DO 430 KLOOP      = 1, KTLOOP
          CNEW(KLOOP,JSPC) = CONC(KLOOP,JSPC)
 430     CONTINUE
 432    CONTINUE
C
        WRITE(6,670) DELT, XELAPS
        !print*, 'kblk = ', kblk      !gcc
        IF (JRESTAR.EQ.100) THEN
         WRITE(6,680)
         CALL GEOS_CHEM_STOP
        ENDIF
C
        GOTO 140
       ENDIF
C
      ELSE
C
C *********************************************************************
C *             ALL SUCCESSFUL STEPS COME THROUGH HERE                * 
C *                                                                   * 
C * AFTER A SUCCESSFUL STEP, UPDATE THE CONCENTRATION AND ALL DERIV-  *
C * ATIVES, RESET TOLD, SET IFSUCCESS = 1, INCREMENT NSTEPS, AND      *
C * RESET JFAIL = 0.                                                  *
C *********************************************************************
C
       !-----------------------------------------------------------------
       ! Prior to 7/9/03:
       ! Comment out counter variable NUM_TIMESTEPS, you can get the same 
       ! info w/ a profiling run. (bmy, 7/9/03)         
       !! timing calculations (bdf, 4/18/03)
       !num_timesteps=num_timesteps+1
       !-----------------------------------------------------------------

       JFAIL            = 0
       IFSUCCESS        = 1
       NSTEPS           = NSTEPS + 1
       TOLD             = XELAPS
C
C *********************************************************************
C
       I1               = 1 
       DO 474 J         = 2, KSTEP
        I1              = I1 + ISCHAN
        ASNQQJ          = ASET(NQQ,J) 
        DO 472 JSPC     = 1, ISCHAN
         I              = JSPC + I1 - 1 
         DO 470 KLOOP   = 1, KTLOOP
          CONC(KLOOP,I) = CONC(KLOOP,I) + ASNQQJ * DTLOS(KLOOP,JSPC)
 470     CONTINUE
 472    CONTINUE
 474   CONTINUE
C
       IF (ASN1.EQ.1.0) THEN 
        DO 473 JSPC         = 1, ISCHAN
         DO 471 KLOOP       = 1, KTLOOP
          CONC( KLOOP,JSPC) = CONC( KLOOP,JSPC) + DTLOS(KLOOP,JSPC)
 471     CONTINUE 
 473    CONTINUE 
       ELSE
        DO 477 JSPC         = 1, ISCHAN
         DO 475 KLOOP       = 1, KTLOOP
          CONC( KLOOP,JSPC) = CONC( KLOOP,JSPC) + ASN1*DTLOS(KLOOP,JSPC)
 475     CONTINUE 
 477    CONTINUE 
       ENDIF 
C
C *********************************************************************
C          EXIT SMVGEAR IF A TIME INTERVAL HAS BEEN COMPLETED
C *********************************************************************
C
       TIMREMAIN        = TINTERVAL - XELAPS
       IF (TIMREMAIN.LE.1.0d-06) GOTO 650

       ! Increment counter of internal timesteps
       ICOUNT = ICOUNT + 1

       ! STOP 700 error -- nonconvergence after many tries
       IF ( DELT < HMIN .or. ICOUNT > MAX_ITERATIONS ) THEN

          WRITE( 6, '(/,a)' ) 'SMVGEAR ERROR -- nonconvergence!'
          WRITE( 6, '(  a)' ) '---------------------------------'

          ! Write DELT and HMIN
          WRITE( 6, 231 ) DELT, HMIN
 231      FORMAT( 'DELT = ', ES10.3, ' HMIN = ', ES10.3 )

          ! List all kinetic and photo reactions
          WRITE( 6, '(/,a)' ) 'Kinetic and Photolysis Reactions:'
          WRITE( 6, '(  a)' ) '---------------------------------'

          DO NK = 1, NALLRAT(NCS)
             WRITE( 6, 248 ) NK, RRATE(KSAVE,NK),TRATE(KSAVE,NK)
 248         FORMAT( 'Rxn #:', i5, ' RRATE = ', es13.6 ,
     x            '  TRATE = ', es13.6 )
          ENDDO

          ! List various SMVGEAR parameters
          WRITE( 6, * ) 'RDELT       = ', RDELT
          WRITE( 6, * ) 'TIMREMAIN   = ', TIMREMAIN
          WRITE( 6, * ) 'HMAX        = ', HMAX
          WRITE( 6, * ) 'RDELMAX     = ', RDELMAX
          WRITE( 6, * ) 'HRATIO      = ', HRATIO

          ! Write offending grid box
          IX = IXSAVE(IJSAVE)
          IY = IYSAVE(IJSAVE)
          IZ = IZSAVE(IJSAVE)         
          WRITE( 6, * ) 'TROUBLE BOX = ', IX, IY, IZ

          ! Nonconvergence after too many iterations
          IF ( ICOUNT > MAX_ITERATIONS ) THEN
             WRITE( 6, * ) 'ICOUNT      = ', ICOUNT
             WRITE( 6, * ) 'Too many iterations!'
          ENDIF
          
          ! Stop w/ error msg 
          WRITE( 6, '(/,a)' ) 'STOP 700 in smvgear.f'
          CALL GEOS_CHEM_STOP
       ENDIF
C
C *********************************************************************
C * IDOUB COUNTS THE NUMBER OF SUCCESSFUL STEPS BEFORE RE-TESTING THE *
C * STEP-SIZE AND ORDER                                               *
C *                                                                   * 
C * IF IDOUB > 1, DECREASE IDOUB AND GO ON TO THE NEXT TIME-STEP WITH *
C *               THE CURRENT STEP-SIZE AND ORDER.                    * 
C * IF IDOUB = 1, STORE THE VALUE OF THE ERROR (DTLOS) FOR THE TIME-  *
C *               STEP PREDICTION, WHICH WILL OCCUR WHEN IDOUB = 0,   *
C *               BUT GO ON TO THE NEXT STEP WITH THE CURRENT STEP-   *
C *               SIZE AND ORDER.                                     * 
C * IF IDOUB = 0, TEST THE TIME-STEP AND ORDER FOR A CHANGE.          *  
C *********************************************************************
C
       IF (IDOUB.GT.1) THEN
        IDOUB               = IDOUB - 1
        IF (IDOUB.EQ.1) THEN
         DO 527 JSPC        = 1, ISCHAN, 2 
          JG1               = JSPC + 1
          DO 525 KLOOP      = 1, KTLOOP
           CEST(KLOOP,JSPC) = DTLOS(KLOOP,JSPC)
           CEST(KLOOP,JG1)  = DTLOS(KLOOP,JG1)
 525      CONTINUE
 527     CONTINUE
        ENDIF
        RDELT               = 1.0d0
        GOTO 170
       ENDIF
C
      ENDIF 
C     ENDIF DER2MAX.GT.ENQQ
C
C *********************************************************************
C *         TEST WHETHER TO CHANGE THE STEP-SIZE AND ORDER            * 
C * DETERMINE THE TIME-STEP AT (A) ONE ORDER LOWER THAN, (B) THE SAME *
C * ORDER AS, AND (C) ONE ORDER HIGHER THAN THE CURRENT ORDER. IN THE *
C * CASE OF MULTIPLE GRID-CELLS IN A GRID-BLOCK, FIND THE MINIMUM     *
C * STEP-SIZE AMONG ALL THE CELLS FOR EACH OF THE ORDERS. THEN, IN    *
C * ALL CASES, CHOOSE THE LONGEST TIME-STEP AMONG THE THREE STEPS     *
C * PAIRED WITH ORDERS, AND CHOOSE THE ORDER ALLOWING THIS LONGEST    *
C * STEP.                                                             *
C *********************************************************************
C
C *********************************************************************
C * ESTIMATE THE TIME-STEP RATIO (RDELTUP) AT ONE ORDER HIGHER THAN   *
C * THE CURRENT ORDER. IF NQQ >= MAXORD, THEN WE DO NOT ALLOW THE     * 
C *                        ORDER TO INCREASE.                         *
C *********************************************************************
C
      IF (NQQ.LT.MAXORD) THEN  
       DO 542 KLOOP  = 1, KTLOOP
        DELY(KLOOP)  = 0.d0
 542   CONTINUE
C
       DO 545 JSPC   = 1, ISCHAN     
        !===============================================================
        ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
        ! prod/loss species from being counted towards the convergence
        ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)
        IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN
           DO 544 KLOOP = 1, KTLOOP 
              ERRYMAX     = (DTLOS(KLOOP,JSPC) - CEST(KLOOP,JSPC)) *
     1                       CHOLD(KLOOP,JSPC)
              DELY(KLOOP) = DELY(KLOOP) + ERRYMAX * ERRYMAX
              if (errymax .gt. specmax) then
                 errymax = specmax
                 jspcsave(kloop) = i
              endif
 544       CONTINUE
        ENDIF
        !===============================================================
 545   CONTINUE
C
       DER3MAX       = 0.d0
C
       DO 546 KLOOP  = 1, KTLOOP
        IF (DELY(KLOOP).GT.DER3MAX) DER3MAX = DELY(KLOOP)
 546   CONTINUE
C
       RDELTUP       = 1.0d0 / (CONP3*DER3MAX**ENQQ3(NQQ)+1.4d-6)
      ELSE
       RDELTUP       = 0.0d0
      ENDIF
C
C *********************************************************************
C *    ESTIMATE THE TIME-STEP RATIO (RDELTSM) AT THE CURRENT ORDER    *
C *      WE CALCULATED DER2MAX DURING THE ERROR TESTS EARLIER         *  
C *********************************************************************
C
 540  RDELTSM        = 1.0d0 / (CONP2*DER2MAX**ENQQ2(NQQ)+1.2d-6)
C
C *********************************************************************
C * ESTIMATE THE TIME-STEP RATIO (RDELTDN) AT ONE ORDER LOWER THAN    *
C * THE CURRENT ORDER. IF NQQ = 1, THEN WE CANNOT TEST A LOWER ORDER. * 
C *********************************************************************
C
      IF (NQQ.GT.1) THEN
       DO 552 KLOOP  = 1, KTLOOP
        DELY(KLOOP)  = 0.d0
 552   CONTINUE
C
       KSTEPISC      = (KSTEP - 1) * ISCHAN
       DO 555 JSPC   = 1, ISCHAN     
        !===============================================================
        ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
        ! prod/loss species from being counted towards the convergence
        ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)
        IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN
           I            = JSPC + KSTEPISC
           DO 554 KLOOP = 1, KTLOOP 
              ERRYMAX     = CONC(KLOOP,I) * CHOLD(KLOOP,JSPC)
              DELY(KLOOP) = DELY(KLOOP)   + ERRYMAX * ERRYMAX
 554       CONTINUE
        ENDIF
        !===============================================================
 555   CONTINUE
C
       DER1MAX       = 0.d0
C
       DO 556 KLOOP = 1, KTLOOP
        IF (DELY(KLOOP).GT.DER1MAX) DER1MAX = DELY(KLOOP)   
 556   CONTINUE
        
       RDELTDN       = 1.0d0 / (CONP1*DER1MAX**ENQQ1(NQQ)+1.3d-6)
C
      ELSE
       RDELTDN       = 0.d0
      ENDIF
C
C *********************************************************************
C * FIND THE LARGEST OF THE PREDICTED TIME-STEPS RATIOS OF EACH ORDER * 
C *********************************************************************
C
      RDELT        = MAX(RDELTUP,RDELTSM,RDELTDN)
C
C *********************************************************************
C * IF THE LAST STEP WAS SUCCESSFUL AND RDELT IS SMALL, KEEP THE      *
C * CURRENT STEP AND ORDER AND ALLOW THREE SUCCESSFUL STEPS BEFORE    *
C * RE-CHECKING THE TIME-STEP AND ORDER.                              *
C *********************************************************************
C
      IF (RDELT.LT.1.1.AND.IFSUCCESS.EQ.1) THEN
       IDOUB         = 3
       GOTO 170
C
C *********************************************************************
C * IF THE MAXIMUM TIME-STEP RATIO IS THAT OF ONE ORDER LOWER THAN    *  
C * THE CURRENT ORDER, DECREASE THE ORDER. DO NOT MINIMIZE RDELT      *
C * TO =< 1 WHEN IFSUCCESS = 0 SINCE THIS IS LESS EFFICIENT.          *
C *********************************************************************
C
      ELSEIF (RDELT.EQ.RDELTDN) THEN
       NQQ              = NQQ - 1
C
C *********************************************************************
C * IF THE MAXIMUM TIME-STEP RATIO IS THAT OF ONE ORDER HIGHER THAN   *  
C * THE CURRENT ORDER, INCREASE THE ORDER AND ADD A DERIVATIVE TERM   *
C * FOR THE HIGHER ORDER.                                             *  
C *********************************************************************
C
      ELSEIF (RDELT.EQ.RDELTUP) THEN
       CONSMULT         = ASET(NQQ,KSTEP) / FLOAT(KSTEP)
       NQQ              = KSTEP
       NQISC            = NQQ * ISCHAN
       DO 602 JSPC      = 1, ISCHAN, 2     
        JG1             = JSPC + 1
        I1              = JSPC + NQISC 
        I2              = JG1  + NQISC 
        DO 600 KLOOP    = 1, KTLOOP 
         CONC(KLOOP,I1) = DTLOS(KLOOP,JSPC) * CONSMULT
         CONC(KLOOP,I2) = DTLOS(KLOOP,JG1)  * CONSMULT
 600    CONTINUE
 602   CONTINUE
      ENDIF
C
C *********************************************************************
C * IF THE LAST TWO STEPS HAVE FAILED, RE-SET IDOUB TO THE CURRENT    *
C * ORDER + 1. DO NOT MINIMIZE RDELT IF JFAIL.GE.2 SINCE TESTS SHOW   *
C * THAT THIS MERELY LEADS TO ADDITIONAL COMPUTATIONS.                *
C *********************************************************************
C
      IDOUB             = NQQ + 1
C
      GOTO 170  
C
C *********************************************************************
C *                      UPDATE COUNTERS                              *
C *********************************************************************
C
 650  NSFTOT             = NSFTOT   + NSUBFUN 
      NPDTOT             = NPDTOT   + NPDERIV
      NSTTOT             = NSTTOT   + NSTEPS 
      IFAILTOT           = IFAILTOT + IFAIL
      NFAILTOT           = NFAILTOT + NFAIL
      LFAILTOT           = LFAILTOT + LFAIL

C
C *********************************************************************
C *        SET FINAL CONCENTRATION FOR RUN AND UPDATE COUNTERS        * 
C *********************************************************************
C
      DO JSPC        = 1, ISCHAN
      DO KLOOP       = 1, KTLOOP

         ! Stop with an error message if NaN's are encountered
         ! (bmy, pip, 4/27/00)
         IF ( IT_IS_NAN( CNEW(KLOOP,JSPC) ) ) THEN
            DO NK = 1, NALLRAT(NCS)
               WRITE( 6, 249 ) NK, RRATE(KSAVE,NK),TRATE(KSAVE,NK)
 249           FORMAT( 'Rxn #:', i5, ' RRATE = ', es13.6 ,
     x              '  TRATE = ', es13.6 )
            ENDDO
            write(6,*) 'sum of rrate = ',sum(rrate)
            PRINT*, 'SMVGEAR: CNEW is NaN!'
            PRINT*, 'Species index : ', JSPC
            PRINT*, 'Grid Box      : ', IXSAVE(KLOOP+JLOOPLO),
     &           IYSAVE(KLOOP+JLOOPLO), IZSAVE(KLOOP+JLOOPLO)
            PRINT*, 'STOP in smvgear.f!'

            ! Stop the run and deallocate all arrays 
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Reset negatives to a very small positive number
         CNEW(KLOOP,JSPC) = MAX(CNEW(KLOOP,JSPC),SMAL2)
      ENDDO
      ENDDO
C
C
C *********************************************************************
C                             FORMATS
C *********************************************************************
C
 233  FORMAT('SMVGEAR: DELT= ',1PE8.2,' TOO LOW DEC YFAC. KBLK, ', 
     1       'KTLOOP, NCS, TIME, TIMREMAIN, YFAC, ', 
     2       'EPS = ',/3(1X,I4),2X,4(1PE9.3,1X))
 234  FORMAT('SMVGEAR: TOO MANY DECREASES OF YFAC ')
 670  FORMAT('DELT DEC TO =',E13.5,'; TIME ',E13.5,' BECAUSE ',
     1       'EXCESSIVE ERRORS')
 680  FORMAT('SMVGEAR: STOP BECAUSE OF EXCESSIVE ERRORS.')
 685  FORMAT('M1,M2,K,ERR = ',3(I4),2X,1PE10.4)
 690  FORMAT('CONC WHEN STOP = ',2(I4,1X),A14,2(1X,1PE10.2))
C
C *********************************************************************
C ***************     END OF SUBROUTINE SMVGEAR     *******************
C *********************************************************************
C
      RETURN

!------------------------------------------------------------------------------
! Prior to 7/8/03:
! Now use a lookup table to test for ND65 families (bmy, 7/8/03)
!      !=================================================================
!      ! INTERNAL FUNCTIONS -- can "see" all variable 
!      ! declarations that are made w/in routine SMVGEAR
!      !=================================================================
!      CONTAINS
!
!!-----------------------------------------------------------------------------
!
!      !FUNCTION ITS_NOT_A_ND65_FAMILY( JSPC ) RESULT( LCOUNT )
!      FUNCTION ND65_FAMILY_TEST( JSPC ) RESULT( LCOUNT )
!!
!!******************************************************************************
!!  Function ITS_NOT_A_ND65_FAMILY returns TRUE if a species is not a ND65
!!  prod.loss family, or FALSE otherwise.  This is needed to keep the "fake"
!!  ND65 prod/loss family names from being included in the convergence
!!  criteria for the SMVGEAR II solver. (ljm, bmy, 5/12/03)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) JSPC (INTEGER) : SMVGEAR species index
!!
!!  NOTES:
!!******************************************************************************
!!
!      ! Arguments
!      INTEGER, INTENT(IN) :: JSPC
!
!      ! Local variables
!      INTEGER             :: N
!
!      ! Function value
!      LOGICAL             :: LCOUNT
!
!      !=================================================================
!      ! ITS_NOT_A_ND65_FAMILY begins here!
!      !=================================================================
!      LCOUNT = .TRUE.
!
!      ! Test if this species is a ND65 prod/loss family name
!      ! For now just assume urban chemistry
!      DO N = 1, NFAMILIES
!         IF ( JSPC == MAPPL(IFAM(N),NCSURBAN) ) THEN
!            LCOUNT = .FALSE.
!            EXIT
!         ENDIF
!      ENDDO  
!
!      ! Return to SMVGEAR
!      !END FUNCTION ITS_NOT_A_ND65_FAMILY
!      END FUNCTION ND65_FAMILY_TEST
!
!-----------------------------------------------------------------------------

      END SUBROUTINE SMVGEAR
