! $Id: smvgear.f,v 1.1 2009/09/16 14:06:08 bmy Exp $
      SUBROUTINE SMVGEAR( am_I_root, Input_Opt, RC )
!
!******************************************************************************
!  Subroutine SMVGEAR solves ODE's for chemical reactions using a GEAR-type
!  method.  (M. Jacobson 1997; bdf, bmy, 5/12/03, 11/1/07)
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
!  (3 ) Declare NSTEPS, KLOOP, IABOVK as local variables since they are only 
!        ever used w/in "smvgear.f"  Also reference ERRMX2 from "comode_mod.f" 
!  (4 ) Now declare DELY, ERRHOLD, YABST, as local variables, since these are
!        used only w/in this routine and nowhere else --- also removed these
!        from /DKBLOOP/ and /DKBLOOP5/ in "comode.h". (bmy, 7/28/03)
!  (5 ) Increase max allowable iteration count to 99999 (mje, bmy, 9/15/03)
!  (6 ) Added trap for negative CNEW values if iteration passes both local and
!        global error tests.  If negative values are found, go back to CNEW 
!        values at the end of the previous successful march, and try again 
!        with smaller timestep.  Now stop at the end of the code if negative 
!        values are encountered in CNEW.  Do not reset negative CNEW values
!        to zero.  (tmf, 11/1/07)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,         ONLY : ERRMX2, IXSAVE, IYSAVE, IZSAVE
      USE COMODE_LOOP_MOD
      USE ERROR_MOD,          ONLY  : IT_IS_NAN, GEOS_CHEM_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! *********************************************************************
! ************        WRITTEN BY MARK JACOBSON (1993)      ************
! ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
! ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
! *********************************************************************
! *********************************************************************
!
!  SSSSSSS   M     M  V       V  GGGGGGG   EEEEEEE      A      RRRRRRR
!  S         MM   MM   V     V   G         E           A A     R     R  
!  SSSSSSS   M M M M    V   V    G  GGGG   EEEEEEE    A   A    RRRRRRR
!        S   M  M  M     V V     G     G   E         AAAAAAA   R  R
!  SSSSSSS   M     M      V      GGGGGGG   EEEEEEE  A       A  R    R 
!
! *********************************************************************
!                    VERSION:      SMVGEAR II
!                    LAST UPDATE:  AUGUST, 1997
! *********************************************************************
!
! *********************************************************************
! * SMVGEAR IS A GEAR-TYPE INTEGRATOR THAT SOLVES FIRST-ORDER ORDIN-  *
! * ARY DIFFERENTIAL EQUATIONS WITH INITIAL VALUE BOUNDARY CONDITIONS.*
! * SMVGEAR DIFFERS FROM AN ORIGINAL GEAR CODE IN THAT IT USES SPARSE *
! * MATRIX AND VECTORIZATION TECHNIQUES TO IMPROVE SPEED. MUCH        * 
! * OF THE SPEED UP IN THIS PROGRAM IS DUE TO SPARSE MATRIX           *
! * TECHNIQUES AND VECTORIZATION.                                     *
! *                                                                   *
! * THIS VERSION INCLUDES RE-ORDERING OF GRID-CELLS PRIOR TO EACH     *
! * TIME-INTERVAL. THE PURPOSE OF THE REORDERING IS TO GROUP CELLS    *
! * WITH STIFF EQUATIONS TOGETHER AND THOSE WITH NON-STIFF EQUATIONS  *
! * THIS REORDERING CAN SAVE SIGNIFCANT COMPUTER TIME                 *
! * (E.G. SPEED THE CODE BY A FACTOR OF TWO OR MORE), DEPENDING ON    *
! * THE VARIATION IN STIFFNESS THROUGHOUT THE GRID-DOMAIN. WHEN THE   *
! * STIFFNESS IS THE SAME THROUGHOUT THE GRID-DOMAIN (E.G. IF ALL     *
! * CONCENTRATIONS AND RATES ARE THE SAME), THEN RE-ORDERING IS       *
! * UNNECESSARY AND WILL NOT SPEED SOLUTIONS.                         *
! *                                                                   *
! * THIS VERSION INCLUDES A VARIABLE ABSOLUTE ERROR TOLERANCE.        *
! * THE ABSOLUTE TOLERANCE IS RECALCULATED EVERY FEW GEAR TIME STEPS. *
! *                                                                   *
! * THIS VERSION CONTAINS DIFFERENT SETS OF CHEMISTRY FOR             *
! * DIFFERENT REGIONS OF THE ATMOSPHERE. THUS, URBAN, FREE TROP-      *
! * OSPHERIC, AND STRATOSPHERIC CHEMISTRY CAN BE SOLVED DURING THE    *
! * SAME MODEL RUN.                                                   * 
! *                                                                   *
! * REFERENCES:                                                       *
! * -----------                                                       *
! *                                                                   * 
! * JACOBSON M. Z. (1998) FUNDAMENTALS OF ATMOSPHERIC MODELING.       *
! *  CAMBRIDGE UNIVERSITY PRESS, NEW YORK.                            *
! *                                                                   * 
! * JACOBSON M. Z. (1998) IMPROVEMENT OF SMVGEAR II ON VECTOR AND     *
! *  SCALAR MACHINES THROUGH ABSOLUTE ERROR TOLERANCE CONTROL.        *    
! *  ATMOS. ENVIRON. 32, 791 - 796                                    *
! *                                                                   * 
! * JACOBSON M. Z. (1995) COMPUTATION OF GLOBAL PHOTOCHEMISTRY        *
! *  WITH SMVGEAR II. ATMOS. ENVIRON., 29A, 2541 - 2546               *
! *                                                                   *
! * JACOBSON M. Z. (1994) DEVELOPING, COUPLING, AND APPLYING A GAS,   *
! *  AEROSOL, TRANSPORT, AND RADIATION MODEL TO STUDYING URBAN        *
! *  AND REGIONAL AIR POLLUTION. Ph. D. THESIS, UNIVERSITY OF         *
! *  CALIFORNIA, LOS ANGELES.                                         *
! *                                                                   *
! * JACOBSON M. Z. AND TURCO R. P. (1994) SMVGEAR: A SPARSE-          * 
! *  MATRIX, VECTORIZED GEAR CODE FOR ATMOSPHERIC MODELS.             *
! *  ATMOS. ENVIRON. 28A, 273 - 284.                                  * 
! *                                                                   *
! * HOW TO CALL SUBROUTINE:                                           *
! * ----------------------                                            *
! *  CALL SMVGEAR FROM PHYSPROC FOR GAS CHEM W/ NCS = 1..NCSGAS       *
! *                                                                   *
! *********************************************************************
! *                                                                   *
! * THE ORIGINS OF THE GEAR INTEGRATOR USED IN SMVGEAR ARE FOUND IN   *
! *                                                                   *    
! * GEAR C. W. (1971) NUMERICAL INITIAL VALUE PROBLEMS IN ORDINARY    *  
! *  DIFFERENTIAL EQUATIONS. PRENTICE-HALL, NJ, PP. 158-166.          * 
! *                                                                   *    
! *********************************************************************
! *                                                                   *      
! * FINALLY, IN SUBROUTINE SMVGEAR.F, THE FOLLOWING IDEAS ORIGINATED  *
! *   FROM LSODES, THE LIVERMORE SOLVER FOR ORDINARY DIFFERENTIAL     *
! *   WITH SPARSE MATRICES (HINDMARSH A. C. AND SHERMAN A. H.):       *
! *                                                                   *      
! *  (A) PREDICTING THE FIRST TIME-STEP;                              *
! *  (B) DETERMINING CORRECTOR CONVERGENCE DIFFERENTLY THAN IN        *
! *      GEAR'S ORIGINAL CODE (GOC)                                   *
! *  (C) DETERMINING ERROR DIFFERENTLY THAN IN GOC                    *
! *  (D) SUMMING UP THE PASCAL MATRIX DIFFERENTLY THAN IN GOC         *      
! *                                                                   *      
! * REFERENCES FOR THE 1987 LSODES VERSION INCLUDE:                   *
! *                                                                   *      
! * SHERMAN A. H. AND HINDMARSH A. C. (1980) GEARS: A PACKAGE FOR     *
! *  THE SOLUTION OF SPARSE, STIFF ORDINARY DIFFERENTIAL EQUATIONS.   *
! *  LAWRENCE LIVERMORE LABORATORY REPORT UCRL-84102.                 *   
! *                                                                   *      
! * HINDMARSH A. C. (1983) ODEPACK, A SYSTEMATIZED COLLECTION OF      *
! *  ODE SOLVERS. IN SCIENTIFIC COMPUTING, R.S. STEPLEMAN ET AL.,     *
! *  EDS., NORTH-HOLLAND, AMSTERDAM, PP. 55 - 74.                     *
! *                                                                   *      
! *********************************************************************
!
! *********************************************************************
! *************** HERE ARE SOME PARAMETER DEFINITIONS *****************
! *********************************************************************
!                                                                          
! ABST2     = 1. / TIMEINTERVAL**2   (SEC-2) (SET IN READER.F) 
! ASN1      = THE VALUE OF ASET(NQQ,1)
! CEST      = STORES VALUE OF DTLOS WHEN IDOUB = 1
! CHOLD     = 1 / (RELTOL * CNEW + ABTOL). MULTIPLY
!             CHOLD BY LOCAL ERRORS IN DIFFERENT ERROR TESTS.
! CNEW      = STORES CONCENTRATION (Y [ESTIMATED])
! CONC      = AN ARRAY OF LENGTH ISCHAN * (MAXORD+1) THAT CARRIES THE
!             DERIVATIVES OF CNEW, SCALED BY DELT**J/FACTORIAL(J),
!             WHERE J IS THE J-TH DERIVATIVE. J VARIES FROM 1 TO NQQ,
!             WHICH IS THE CURRENT ORDER OF THE METHOD.
!             E.G. CONC(JSPC,2) STORES DELT * Y' (ESTIMATED)                   
! DELT      = CURRENT TIME-STEP (S) LENGTH DURING A TIME-INTERVAL 
! DRATE     = PARAMETER WHICH USED TO DETERMINE WHETHER CONVERGENCE 
!             HAS OCCURRED
! DTLOS     = AN ARRAY OF LENGTH ISCHAN, USED FOR THE ACCUMULATED
!             CORRECTIONS.  ON A SUCCESSFUL RETURN, DTLOS(KLOOP,I) CONTAINS
!             THE ESTIMATED ONE-STEP LOCAL ERROR IN CNEW.
! EDWN      = PERTST**2 * ORDER FOR ONE ORDER LOWER THAN CURRENT ORDER
! ENQQ      = PERTST**2 * ORDER FOR CURRENT ORDER
! ERRMAX    = RELATIVE ERROR TOLERANCE (SEE CHOLD). SET IN m.dat.
!             EPS SHOULD BE < 1.0. FOR SPEEDY AND RELIABLE RESULTS, 
!             10**-3 IS REASONABLE. FOR MANY DECIMAL PLACES OF ACCURACY,
!             DECREASE EPS. 
! EUP       = PERTST**2 * ORDER FOR ONE ORDER HIGHER THAN CURRENT ORDER 
! FRACDEC   = FRACTION THE TIME-STEP IS DECREASED IF CONVERGENCE TEST FAILS
! GLOSS     = VALUE OF FIRST DERIVATIVES ON OUTPUT FROM SUBFUN.
!           = RIGHT-SIDE OF EQUATION ON INPUT TO BACKSUB.F 
!           = ERROR TERM (SOLUTION FROM BACKSUB.F) ON OUTPUT FROM BACKSUB
! HMAX      = THE MAXIMUM ALLOWABLE VALUE OF DELT
! HMIN      = THE MINIMUM ALLOWABLE VALUE OF DELT
! HRMAX     = MAXIMUM RELATIVE CHANGE IN DELT*ASET(1) BEFORE PDERIV IS CALLED.
! HRATIO    = RELATIVE CHANGE IN DELT * ASET(1) EACH CHANGE IN STEP OR ORDER
!             WHEN ABS(HRATIO-1) > HRMAX, RESET JEVAL = 1 TO CALL PDERIV
! IABOVK    = NUMBER OF SPECIES WHOSE CONCENTRATIONS ARE LARGER THAN YABST
! IDOUB     = RECORDS THE NUMBER OF STEPS SINCE THE LAST CHANGE IN STEP SIZE   
!             OR ORDER.  IT MUST BE AT LEAST KSTEP = NQQ+1 BEFORE DOUBLING IS 
!             ALLOWED. 
! IFAIL     = NUMBER OF TIMES THE CORRECTOR FAILED TO CONVERGE WHILE THE
!             JACOBIAN WAS OLD (PDERIV NOT CALLED DURING THE LAST TEST)
! IFSUCCESS = IDENTIFIES WHETHER STEP IS SUCCESSFUL (=1) OR NOT (=0)
! IFSUN     = IDENTIFIES WHETHER SUN IS UP (=1) OR DOWN (=2)
! ISCHAN    = THE NUMBER OF FIRST-ORDER EQUATIONS TO SOLVE = # OF SPECIES = 
!             ORDER OF ORIGINAL MATRIX. ISCHAN HAS A DIFFERENT VALUE
!             FOR DAY AND NIGHT AND FOR GAS- CHEMISTRY.
! ISREORD   = 1: CALC INITIAL STIFFNESS BEFORE RUNNING CODE TO REORDER CELLS
!                IN THIS CASE, USE PHOTORATES FOR END OF TIME-INTERVAL
!           = 0: DO NORMAL CALCULATIONS
! JEVAL     = 1  --> CALL PDERIV THE NEXT TIME THROUGH THE CORRECTOR STEPS.
!           = 0  --> LAST STEP WAS SUCCESSFUL AND DO NOT NEED TO CALL PDERIV
!           = -1 --> PDERIV JUST CALLED, AND DO NOT NEED TO CALL AGAIN
!             UNTIL JEVAL SWITCHED TO 1. 
! JRESTAR   = COUNTS NUMBER OF TIMES SMVGEAR STARTS OVER AT ORDER 1
!             BECAUSE OF EXCESSIVE FAILURES.
! LFAIL     = NUMBER OF TIMES THE ACCUMULATED ERROR TEST FAILED
! KSTEP     = NQQ + 1
! KTLOOP    = NUMBER OF GRID-CELLS IN A GRID-BLOCK
! MAXORD    = THE MAXIMUM ALLOWABLE ORDER OF THE INTEGRATION METHOD
! MBETWEEN  = THE MAXIMUM ALLOWABLE NUMBER OF STEPS BETWEEN CALLS TO PDERIV
! MSTEP     = THE MAXIMUM ALLOWABLE NUMBER OF CORRECTOR ITERATIONS
! NCS       = 1..NCSGAS FOR GAS CHEMISTRY                            
! NCSP      = NCS       FOR DAYTIME   GAS CHEM            
!           = NCS + ICS FOR NIGHTTIME GAS CHEM           
! NFAIL     = NUMBER OF TIMES CORRECTER FAILS TO CONVERGE AFTER PDERIV
!             WAS JUST CALLED
! NPDERIV   = TOTAL NUMBER OF TIMES THAT MATRIX IS EVALUATED (PDERIV)
! NPDTOT    = NUMBER OF CALLS TO PDERIV ROUTINE, OVER ALL TIME
! NSFTOT    = NUMBER OF CALLS TO SUBFUN ROUTINE, OVER ALL TIME
! NSLP      = THE LAST TIME-STEP NUMBER DURING WHICH PDERIV WAS CALLED 
! NSTTOT    = TOTAL NUMBER OF SUCCESSFUL TIME-STEPS, OVER ALL TIME 
! NSUBFUN   = TOTAL NUMBER OF TIMES SUBFUN IS CALLED
! NSTEPS    = TOTAL NUMBER OF SUCCESSFUL TIME-STEPS TAKEN
! NQQ       = ORDER OF THE INTEGRATION METHOD. IT VARIES BETWEEN 1 AND MAXORD. 
! NQQISC    = NQQ * ISCHAN
! NQQOLD    = VALUE OF NQQ DURING LAST TIME-STEP
! ORDER     = FLOATING POINT VALUE OF ISCHAN, THE ORDER OF NUMBER OF ODES.
! PDERIV    = NAME OF ROUTINE TO EVALUATE THE JACOBIAN MATRIX (J) 
!             AND P = I - DELT * ASET(1) * J 
! PERTS2    = COEFFICIENTS USED IN SELECTING THE STEP AND ORDER (SEE
!             JSPARSE.F) NOTE THAT PERTS2 = ORIGINAL PERTST**2     
! RDELMAX   = THE MAXIMUM FACTOR BY WHICH DELT CAN BE INCREASED IN A SINGLE 
!             STEP.  AS IN LSODES, SET IT TO 1E4 INITIALLY TO COMPENSATE 
!             FOR THE SMALL INITIAL DELT, BUT THEN SET IT TO 10 AFTER 
!             SUCCESSFUL STEPS AND 2 AFTER UNSUCCESSFUL STEPS
! RDELT     = FACTOR (TIME-STEP RATIO) BY WHICH WE INCREASE OR DECREASE DELT
! RDELTDN   = TIME-STEP RATIO AT ONE ORDER LOWER THAN CURRENT ORDER 
! RDELTSM   = TIME-STEP RATIO AT CURRENT ORDER 
! RDELTUP   = TIME-STEP RATIO AT ONE ORDER HIGHER THAN CURRENT ORDER 
! RMSRAT    = RATIO OF CURRENT TO PREVIOUS RMS SCALED ERROR. IF THIS
!             RATIO DECREASES, THEN CONVERGENCE IS OCCURING.
! SUBFUN    = NAME OF ROUTINE TO SOLVE FIRST DERIVATIVES.
!           = EVALUATES DERIVATIVES IN THE SPECIAL FORM F = Y'(EST)
!           = F(X,Y,ESTIMATED), WHERE F IS THE RIGHT HAND SIDE OF THE
!             DIFFERENTIAL EQUATION. 
! TINTERVAL = TOTAL SECONDS IN A TIME-INTERVAL
! TIMREMAIN = REMAINING TIME IN AN INTERVAL 
! TOLD      = STORES THE LAST VALUE OF XELAPS IN CASE THE CURRENT STEP FAILS
! XELAPS    = ELAPSED TIME IN AN INTERVAL (S)
! ABTOL     = ABSOLUTE ERROR TOLERANCE 
!             IF ABTOL IS TOO SMALL, THEN INTEGRATION WILL TAKE TOO LONG.
!             IF ABTOL TOO LARGE, CONVERGENCE WILL BE TOO EASY AND ERRORS
!             WILL ACCUMULATE, THE TIME-STEP MAY BE CUT TOO SMALL, AND
!             THE INTEGRATION MAY STOP (DELT < HMIN OR FLOATING POINT
!             EXCEPTION IN DECOMP.F).
!             TYPICAL GAS-PHASE VALUES OF ABSTOL ARE 10**3 CM-3
!             TYPICAL AQ -PHASE VALUES OF ABSTOL ARE 10**-13 TO 10**-15 M L-1
! YFAC      = 1.0 ORIGINIALLY, BUT IS DECREASED IF EXCESSIVE FAILURES OCCUR
!             IN ORDER TO REDUCE ABSOLUTE ERROR TOLERANCE 
! *********************************************************************
!
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
  
      !====================================
      ! Additional variable declarations
      !====================================

      ! Add counter
      INTEGER            :: ICOUNT, NK
      
      !-----------------------------------------------------------------
      ! Added for the ND65 prod/loss diagnostic (ljm, bmy, 5/9/03)
      INTEGER            :: NNOFAM
      !-----------------------------------------------------------------

      ! ljm stop 700 trouble
      INTEGER            :: IJSAVE, JSPCSAVE(KTLOOP)
      INTEGER            :: IX, IY, IZ, JJ, JJJ, KSAVE, COUNTER
      REAL*8             :: SPECMAX

      ! Maximum iteration count for SMVGEAR (bmy, 4/11/03)
      INTEGER, PARAMETER :: MAX_ITERATIONS = 99999
      
      ! Variables from "comode.h" which are only ever used in "smvgear.f"
      ! Remove them from "comode.h" and the THREADPRIVATE declarations
      ! (bmy, 7/28/03)
      INTEGER            :: NSTEPS
      INTEGER            :: KGRP(KBLOOP,5), IABOVK(KBLOOP)   
      REAL*8             :: DELY(KBLOOP),   ERRHOLD(KBLOOP)
      REAL*8             :: YABST(KBLOOP)

      !-----------------------------------------------------------------------
      ! %%%%% MODIFICATION TO PREVENT NEGATIVE CNEW (tmf, 11/1/07) %%%%%
      ! INEG = flag for identifying negative values in CNEW
      ! If IDTFORCE==1, then use DTFORCE as timestep instead of DELT
      ! CPREVM stores CNEW at the end of the previous successful march
      INTEGER            :: INEG, IDTFORCE 
      REAL*8             :: DTFORCE  
      REAL*8             :: CPREVM( KBLOOP, MXGSAER )
      !-----------------------------------------------------------------------

      !=================================================================
      ! SMVGEAR begins here!
      !=================================================================
      
      ! Assume success
      RC        = GIGC_SUCCESS

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

      !-----------------------------------------------------------------------
      ! %%%%% MODIFICATION TO PREVENT NEGATIVE CNEW (tmf, 11/1/07) %%%%%
      ! Initialize 
      IDTFORCE = 0
      DTFORCE  = 0.d0
      !-----------------------------------------------------------------------
      ! Added for the ND65 prod/loss diagnostic, in order to prevent
      ! ND65 prod/loss families from being counted towards the 
      ! SMVGEAR convergence criteria. (ljm, bmy, 5/9/03)
      NNOFAM    = ISCHAN - NFAMILIES
      ORDER     = REAL( NNOFAM )
      !-----------------------------------------------------------------------
!
      IABOVE    = ORDER * 0.4d0
!
      DO 115 KLOOP   = 1, KTLOOP
       IABOVK(KLOOP) = IABOVE
 115  CONTINUE
!
      HRMAX     = 0.3d0
      HMAX      = HMAXUSE( NCSP)
      YFAC      = 1.0d0
      ERRINIT   = MIN(ERRMAX(NCS),1.0D-03)
!
! *********************************************************************
!            START TIME INTERVAL OR RE-ENTER AFTER TOTAL FAILURE
! *********************************************************************
!
      ! EXPLANATORY NOTE: Internal timestep loop begins here (tmf, 11/1/07)
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
!
! *********************************************************************
!                  INITIALIZE CONCENTRATION ARRAY 
! *********************************************************************
! CORIG  = ORIGINAL CONCENTRATIONS, WHICH DO NOT CHANGE IN SMVGEAR
! CNEW   = FINAL CONCENTRATIONS, CALCULATED IN SMVGEAR
! CPREVM = VALUE OF CNEW AT END OF LAST SUCCESSFUL MARCH
!

      DO 129 JNEW         = 1, ISCHAN
       DO 127 KLOOP       = 1, KTLOOP
        CNEW( KLOOP,JNEW) = CORIG(KLOOP,JNEW)
        !---------------------------------------------------------------------
        ! %%%%% MODIFICATION TO PREVENT NEGATIVE CNEW (tmf, 11/1/07) %%%%%
        ! Initialize CPREVM with CORIG before 1st internal march
        CPREVM( KLOOP,JNEW) = CORIG(KLOOP,JNEW)
        !---------------------------------------------------------------------
 127   CONTINUE
 129  CONTINUE
      
!
! *********************************************************************
!  RE-ENTER HERE IF TOTAL FAILURE OR IF RESTARTING WITH NEW CELL BLOCK 
! *********************************************************************
!
 140  HRATIO    = 0.d0
      ASN1      = 1.d0
      IFSUCCESS = 1
      RDELMAX   = 1.0d+04

! *********************************************************************
!                         INITIALIZE PHOTRATES 
! *********************************************************************
!
      ! Called for photorates with no active loss terms (bdf, 4/18/03)
      IF (IFSUN.EQ.1) CALL UPDATE
!
! *********************************************************************
!               INITIALIZE FIRST DERIVATIVE FOR CHEMISTRY 
! *********************************************************************
!
      CALL SUBFUN


!
! *********************************************************************
!                DETERMINE INITIAL ABSOLUTE ERROR TOLERANCE 
! *********************************************************************
! IABOVK  = NUMBER OF SPECIES WHOSE CONCENTRATIONS ARE LARGER THAN YABST
! ISREORD = 1: CALC INITIAL STIFFNESS BEFORE RUNNING CODE TO REORDER CELLS
!              IN THIS CASE, USE PHOTORATES FOR END OF TIME-INTERVAL
!         = 2: DO NORMAL CALCULATIONS
! KGRP    = COUNTS NUMBER OF CONCENTRATIONS ABOVE ABTOL(I), I = 1..  
! YABST   = ABSOLUTE ERROR TOLERANCE (MOLEC. CM-3 FOR GASES) 
! ABTOL   = PRE-DEFINED ABSOLUTE ERROR TOLERANCES 
!
      DO 142 KLOOP    = 1, KTLOOP
       ERRHOLD(KLOOP) = 0.d0
 142  CONTINUE 
!
      ! EXPLANATORY NOTE: 
      ! We don't reorder: IFREORD=0 in mglob.dat, (tmf, 11/1/07)
      IF (ISREORD.NE.1) THEN
!
       DO 134 K              = 1, 5        
        DO 132 KLOOP         = 1, KTLOOP
         KGRP(KLOOP,K)       = 0
 132    CONTINUE
 134   CONTINUE
!
       DO 136 JSPC           = 1, ISCHAN
        !---------------------------------------------------------------------
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
        !---------------------------------------------------------------------
 136   CONTINUE

!

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
!
       DO 139 JSPC      = 1, ISCHAN
         !--------------------------------------------------------------------
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
         !--------------------------------------------------------------------
 139   CONTINUE

!
      ELSE

! EXPLANATORY NOTE: This is not used (tmf, 11/1/07)
!
! *********************************************************************
!          USE LOWEST ABSOLUTE ERROR TOLERANCE WHEN REORDERING 
!          IF REORDERING, SET ERRMX2 THEN RETURN TO PHYSPROC.F 
! *********************************************************************
! ABTOLER1 = YFAC * ABTOL(6,NCS) / MIN(ERRMAX,1.0E-03) 
! 
       DO 144 JSPC      = 1, ISCHAN  
         !--------------------------------------------------------------------
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
         !--------------------------------------------------------------------
 144   CONTINUE
!
       IF (ISREORD.EQ.1) THEN 
        DO 150 KLOOP           = 1, KTLOOP 
         ERRMX2(JLOOPLO+KLOOP) = ERRHOLD(KLOOP)
 150    CONTINUE
!
        RETURN
       ENDIF
      ENDIF

! EXPLANATORY NOTE (tmf, 11/1/07)
! This is the end of the ISREORDER condition.

!
! *********************************************************************
!               CALCULATE INITIAL TIME STEP SIZE (S) 
! *********************************************************************
! SQRT(ERRHOLD / [ERRINIT * ORDER]) = RMSNORM OF ERROR SCALED TO ERRINIT 
!                                * CNEW + ABTOL/RELTOL
!
      RMSTOP         = 0.d0
!
      DO 151 KLOOP   = 1, KTLOOP
       IF (ERRHOLD(KLOOP).GT.RMSTOP) RMSTOP = ERRHOLD(KLOOP)   
 151  CONTINUE
!
      DELT1    = SQRT(ERRINIT / (ABST2(NCS) + RMSTOP / ORDER)) 

      !-----------------------------------------------------------------------
      ! %%%%% MODIFICATION TO PREVENT NEGATIVE CNEW (tmf, 11/1/07) %%%%%
      !
      ! If IDTFORCE==0 then compute DELT w/ the original method
      ! If IDTFORCE==1 then manually set DELT to DTFORCE 
      IF ( IDTFORCE == 0 ) THEN 
         DELT  = MAX( MIN( DELT1, TIMREMAIN, HMAX ), HMIN ) 
      ELSE
         DELT  = DTFORCE   
      ENDIF

      ! Reset IDTFORCE for next internal timestep march
      IDTFORCE = 0
      !-----------------------------------------------------------------------
!
! *********************************************************************
!                      SET INITIAL ORDER TO 1
! *********************************************************************
!
      NQQOLD       = 0
      NQQ          = 1 
      JEVAL        = 1
      RDELT        = 1.0d0

!
! *********************************************************************
! *   STORE INITIAL CONCENTRATION AND FIRST DERIVATIVES x TIME-STEP   * 
! *********************************************************************
!
      DO 155 JSPC        = 1, ISCHAN 
       JS1               = ISCHAN + JSPC
       DO 154 KLOOP      = 1, KTLOOP
        CONC(KLOOP,JSPC) = CNEW(KLOOP,JSPC)
        CONC(KLOOP,JS1)  = DELT * GLOSS(KLOOP,JSPC) 
 154   CONTINUE
 155  CONTINUE

!
! *********************************************************************
! ** UPDATE COEFFICIENTS OF THE ORDER. NQQ IS THE ORDER. ASET AND    **
! ** PERTS2 ARE DEFINED IN SUBROUTINE KSPARSE. NOTE THAT PERTS2      **
! ** IS THE ORIGINAL PERTST**2                                       **
! *********************************************************************
!

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

!
! *********************************************************************
!   LIMIT SIZE OF RDELT, THEN RECALCULATE NEW TIME STEP AND UPDATE
! HRATIO. USE HRATIO TO DETERMINE WHETHER PDERIV SHOULD BE CALLED AGAIN 
! *********************************************************************
!
      HMTIM         = MIN(HMAX,TIMREMAIN) 
      RDELT         = MIN(RDELT,RDELMAX,HMTIM/DELT)
      DELT          = DELT   * RDELT 
      HRATIO        = HRATIO * RDELT 
      XELAPS        = XELAPS + DELT  


!
      IF (ABS(HRATIO-1.0).GT.HRMAX.OR.NSTEPS.GE.NSLP) JEVAL = 1
!
! *********************************************************************
!      IF TIME STEP < HMIN, TIGHTEN ABSOLOUTE ERROR TOLERANCE AND 
!          RESTART INTEGRATION AT BEGINNING OF TIME INTERVAL
! *********************************************************************
!
      IF (DELT.LT.HMIN) THEN
! Add CPU # to output (bmy, 3/8/13)
!       WRITE(6,233)DELT,KBLK,KTLOOP,NCS,TIME,TIMREMAIN,YFAC,ERRMAX(NCS)
       WRITE(6,233) DELT,Input_Opt%myCpu,KBLK,KTLOOP,NCS,TIME,
     &             TIMREMAIN,YFAC,ERRMAX(NCS)
       NYLOWDEC  = NYLOWDEC + 1
       YFAC      = YFAC  * 0.01d0
!
       IF (NYLOWDEC.EQ.10) THEN
        LLOOPA      = 1
        LLOOPB      = KTLOOP
        WRITE(6,234)
!
        DO 177 KLOOP = 1, KTLOOP
         JLOOP       = JREORDER(JLOOPLO+KLOOP)
         K           = (JLOOP - 1) / NLOOP + 1
         MLOOP       = JLOOP - (K - 1) * NLOOP
         M1          = (MLOOP - 1) / NLONG + 1
         M2          = MLOOP - (M1 - 1) * NLONG
! Add CPU # to output (bmy, 3/8/13)
!         WRITE(6,685) M1, M2, K, ERRHOLD(KLOOP)
         WRITE(6,685) Input_Opt%myCpu, M1, M2, K, ERRHOLD(KLOOP)
 177    CONTINUE
!
        DO 178 JNEW = 1, ISCHAN
         JOLD       = INEWOLD(JNEW,NCS)
! Add CPU # to output (bmy, 3/8/13)
!         WRITE(6,690) JNEW, NCS, NAMENCS(JOLD,NCS),CORIG(LLOOPA,JNEW),
!     1                CORIG(LLOOPB,JNEW)
         WRITE(6,690) Input_Opt%myCpu, JNEW, NCS, 
     &                NAMENCS(JOLD,NCS),CORIG(LLOOPA,JNEW),
     &                CORIG(LLOOPB,JNEW)
 178    CONTINUE

        ! Stop run w/ error msg
        CALL GEOS_CHEM_STOP
       ENDIF
!
       GOTO 120
      ENDIF

!
! *********************************************************************
! * IF THE DELT IS DIFFERENT THAN DURING THE LAST STEP (IF RDELT NE   *
! * 1), THEN SCALE THE DERIVATIVES                                    *
! *********************************************************************
!
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
!
! *********************************************************************
! * UPDATE PHOTO RATES BECAUSE THE TIME CHANGED.                      *
! * NOTE THAT A TIME CHANGE COULD CORRESPOND TO EITHER A SUCCESSFUL   *
! * OR FAILED STEP                                                    * 
! *********************************************************************
!
      ! Called for photorates with no active loss terms (bdf, 4/18/03)
      IF (IFSUN.EQ.1.AND.XELAPS.NE.XELAPLAST) CALL UPDATE
!
! *********************************************************************
! * IF THE LAST STEP WAS SUCCESSFUL, RESET RDELMAX = 10 AND UPDATE    *
! * THE CHOLD ARRAY WITH CURRENT VALUES OF CNEW.                      * 
! *********************************************************************
!
      IF (IFSUCCESS.EQ.1) THEN
       RDELMAX             = 10.d0
!
! *********************************************************************
!                DETERMINE NEW ABSOLUTE ERROR TOLERANCE 
! *********************************************************************
! KGRP    = COUNTS NUMBER OF CONCENTRATIONS ABOVE ABTOL(I), I = 1..  
! YABST   = ABSOLUTE ERROR TOLERANCE (MOLEC. CM-3 FOR GASES) 
! ABTOL   = PRE-DEFINED ABSOLUTE ERROR TOLERANCES 
!
! EXPLANATORY NOTES (tmf, 11/1/07)
! (1) If the previous step is successful, then redetermine
!     the absolute error tolerance every on 3rd step. 
! *********************************************************************
!
       IF (MOD(NSTEPS,3).EQ.2) THEN
        DO 203 K              = 1, 5        
         DO 201 KLOOP         = 1, KTLOOP
          KGRP(KLOOP,K)       = 0
 201     CONTINUE
 203    CONTINUE
!
        DO 207 JSPC           = 1, ISCHAN
         !--------------------------------------------------------------------
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
         !--------------------------------------------------------------------
 207    CONTINUE
!
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
!
! EXPLANATORY NOTES: (tmf, 11/1/07)
! (1) What is CHOLD?
!     CHOLD = 1./( (Relative tolerance) * N + (Absolute tolerance * Yfac) )
       DO 213 JSPC         = 1, ISCHAN
        !---------------------------------------------------------------------
        ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
        ! prod/loss species from being counted towards the convergence
        ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)
        IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN
           DO 211 KLOOP         = 1, KTLOOP
             CHOLD(KLOOP,JSPC)  = RELTOL3 / (MAX(CNEW(KLOOP,JSPC),0.D0)
     1                          + YABST(KLOOP) * RELTOL2)
 211      CONTINUE 
        ENDIF
        !---------------------------------------------------------------------
 213   CONTINUE 
!
      ENDIF 

!     ENDIF IFSUCCESS.EQ.1
!
! *********************************************************************
! * COMPUTE THE PREDICTED CONCENTRATION AND DERIVATIVES BY MULTIPLY-  *
! * ING PREVIOUS VALUES BY THE PASCAL TRIANGLE MATRIX.                * 
! *********************************************************************
! THIS SET OF OPERATIONS IS EQUIVALENT TO THE REVERSE OF LOOP 419.
! THE EXPANSION OF THE PASCAL TRIANGLE MATRIX WAS CALCULATED BY B. SCHWARTZ.
! THE FIRST DERIVATIVE MULTIPLIED BY THE TIME STEP IS THE SUM 
! OF TERMS ADDED TO CONC(KLOOP,I)
!
      IF (NQQ.EQ.1) THEN
       DO 236 I          = 1, ISCHAN
        J                = I + ISCHAN
        DO 235 KLOOP     = 1, KTLOOP
         CONC(KLOOP,I)   = CONC(KLOOP,I) + CONC(KLOOP,J)
 235    CONTINUE
 236   CONTINUE
! 
      ELSEIF (NQQ.EQ.2) THEN
!
       DO 238 I         = 1, ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN
        DO 237 KLOOP    = 1, KTLOOP
         CONC(KLOOP, I) = CONC(KLOOP, I) + CONC(KLOOP,J1)
     1                  +                  CONC(KLOOP,J2)
         CONC(KLOOP,J1) = CONC(KLOOP,J1) + CONC(KLOOP,J2) * 2.d0
 237    CONTINUE
 238   CONTINUE
!
      ELSEIF (NQQ.EQ.3) THEN
!
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
!
      ELSEIF (NQQ.EQ.4) THEN
!
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
!
      ELSEIF (NQQ.EQ.5) THEN
!
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

!
! *********************************************************************
! ************************** CORRECTION LOOP **************************
! * TAKE UP TO 3 CORRECTOR ITERATIONS. TEST CONVERGENCE BY REQUIRING  *
! * THAT CHANGES BE LESS THAN THE RMS NORM WEIGHTED BY CHOLD.         * 
! * ACCUMULATE THE CORRECTION IN THE ARRAY DTLOS(). IT EQUALS THE     *
! * THE J-TH DERIVATIVE OF CONC() MULTIPLIED BY DELT**KSTEP /         *
! * (FACTORIAL(KSTEP-1)*ASET(KSTEP)); THUS, IT IS PROPORTIONAL TO THE *
! * ACTUAL ERRORS TO THE LOWEST POWER OF DELT PRESENT (DELT**KSTEP)   *
! *********************************************************************
!
 220  L3                  = 0
      DO 232 JSPC         = 1, ISCHAN
       DO 230 KLOOP       = 1, KTLOOP
        CNEW(KLOOP,JSPC)  = CONC(KLOOP,JSPC)
        DTLOS(KLOOP,JSPC) = 0.d0
 230   CONTINUE
 232  CONTINUE

!
! *********************************************************************
! * IF JEVAL = 1, RE-EVALUATE PREDICTOR MATRIX P = I - H * ASET(1) *J *
! * BEFORE STARTING THE CORRECTOR ITERATION. AFTER CALLING PDERIV,    * 
! * SET JEVAL = -1 TO PREVENT RECALLING PDERIV UNLESS NECESSARY LATER.*
! * CALL DECOMP TO DECOMPOSE THE MATRIX                               *
! *********************************************************************
!
      IF (JEVAL.EQ.1) THEN
       R1DELT   = -ASN1 * DELT
!
       CALL PDERIV
!     
       CALL DECOMP 
       JEVAL    = -1 
       HRATIO   = 1.0d0
       NSLP     = NSTEPS + MBETWEEN
       DRATE    = 0.7d0

      ENDIF

!
! *********************************************************************
! *   EVALUATE THE FIRST DERIVATIVE USING CORRECTED VALUES OF CNEW    * 
! *********************************************************************
!     
 270  CALL SUBFUN

!
! *********************************************************************
! * IN THE CASE OF THE CHORD METHOD, COMPUTE ERROR (GLOSS) FROM THE   * 
! * CORRECTED CALCULATION OF THE FIRST DERIVATIVE                     * 
! *                                                                   *
! * EXPLANATORY NOTES (tmf, 11/1/07)                                  *
! * (1) GLOSS now changes from 1st derivative to error here           *
! * (2) GLOSS is now the B matrix in Px = B                           *
! *      (Equation 3 in Jacobson & Turco 1994)                        *
! * (3) DTLOS is the accumulation of deltaN (x array in Px = B)       *
! *********************************************************************
!
      DO 362 JSPC         = 1,  ISCHAN
       J                  = JSPC + ISCHAN
       DO 360 KLOOP       = 1, KTLOOP
        GLOSS(KLOOP,JSPC) = DELT * GLOSS(KLOOP,JSPC)  
     1                    - (CONC(KLOOP,J) + DTLOS(KLOOP,JSPC))
 360   CONTINUE
 362  CONTINUE


!
! *********************************************************************
! * SOLVE THE LINEAR SYSTEM OF EQUATIONS WITH THE CORRECTOR ERROR.    *
! * BACKSUB.F SOLVES BACKSUBSTITUTION OVER MATRIX OF PARTIAL DERIVS.  * 
! *                                                                   *
! * EXPLANATORY NOTES (tmf, 11/1/07)                                  *
! * (1) BACKSUB solves Px = B.                                        *
! * (2) GLOSS is now the solution B array.                            *
! *********************************************************************
!
      CALL BACKSUB
!
! *********************************************************************
! * SUM-UP THE ACCUMULATED ERROR, CORRECT THE CONCENTRATION WITH THE  * 
! * ERROR, AND BEGIN TO CALCULATE THE RMSNORM OF THE ERROR RELATIVE   *  
! * TO CHOLD.                                                         *
! *                                                                   *
! * EXPLANATORY NOTES (tmf, 11/1/07)                                  *
! * In loops 366..368 and 370..372 we do the following:               *
! * (1) Calculate local error                                         * 
! * (2) Update DTLOS = accumulation of delta N                        *
! * (3) Update CNEW = concentration array                             *
! *********************************************************************
!
      DO 365 KLOOP     = 1, KTLOOP
       DELY(KLOOP)     = 0.d0
 365  CONTINUE
!
      specmax = 0.0d0
      IF (ASN1.EQ.1.0) THEN 
       DO 368 I         = 1, ISCHAN
           DO 366 KLOOP    = 1, KTLOOP
              DTLOS(KLOOP,I) = DTLOS(KLOOP,I)  + GLOSS(KLOOP,I)
              CNEW(KLOOP,I)  = CONC(KLOOP,I)   + DTLOS(KLOOP,I)
              !---------------------------------------------------------------
              ! Added for the ND65 prod/loss diagnostic.  This prevents 
              ! ND65 prod/loss species from being counted towards the 
              ! convergence criteria for SMVGEAR (ljm, bmy, 5/9/03)
              IF ( ITS_NOT_A_ND65_FAMILY(I) ) THEN
                 ERRYMAX        = GLOSS(KLOOP,I)  * CHOLD(KLOOP,I)
                 DELY(KLOOP)    = DELY(KLOOP)     + ERRYMAX * ERRYMAX
              ENDIF
              !---------------------------------------------------------------
 366       CONTINUE
 368   CONTINUE
      ELSE
       DO 372 I         = 1, ISCHAN
           DO 370 KLOOP    = 1, KTLOOP
              DTLOS(KLOOP,I) = DTLOS(KLOOP,I)  + GLOSS(KLOOP,I)
              CNEW(KLOOP,I)  = CONC(KLOOP,I)   + ASN1  *  DTLOS(KLOOP,I)
              !---------------------------------------------------------------
              ! Added for the ND65 prod/loss diagnostic.  This prevents 
              ! ND65 prod/loss species from being counted towards the 
              ! convergence criteria for SMVGEAR (ljm, bmy, 5/9/03)
              IF ( ITS_NOT_A_ND65_FAMILY(I) ) THEN
                 ERRYMAX        = GLOSS(KLOOP,I)  * CHOLD(KLOOP,I)
                 DELY(KLOOP)    = DELY(KLOOP)     + ERRYMAX * ERRYMAX
              ENDIF
              !---------------------------------------------------------------
 370       CONTINUE
 372   CONTINUE
      ENDIF


!
! *********************************************************************
! * SET THE PREVIOUS RMS ERROR AND CALCULATE THE NEW RMS ERROR.       *
! * IF DCON < 1, THEN SUFFICIENT CONVERGENCE HAS OCCURRED. OTHERWISE, *
! * IF THE RATIO OF THE CURRENT TO PREVIOUS RMSERR IS DECREASING,     *
! * ITERATE MORE. IF IT IS NOT, THEN THE CONVERGENCE TEST FAILED      *
! *                                                                   *
! * EXPLANATORY NOTES (tmf, 11/1/07)                                  *
! * (1) ORDER is the # of species, excluding ND65 families            *
! * (2) RMSRAT = (local error this march) / (local error last march)  *
! *********************************************************************
!
      RMSERRP           = RMSERR
      DER2MAX           = 0.d0
!
      ksave=0d0
      DO 427 KLOOP      = 1, KTLOOP
       !-----------------------------------------------------------------------
       ! %%%%% MODIFICATION TO PREVENT NEGATIVE CNEW (tmf, 11/1/07) %%%%%
       ! Put STOP 700 debug variables into IF statement 
       !IF (DELY(KLOOP).GT.DER2MAX) DER2MAX = DELY(KLOOP) 
       !-----------------------------------------------------------------------
       IF (DELY(KLOOP).GT.DER2MAX) THEN
          DER2MAX = DELY(KLOOP)   
          ijsave=jlooplo+kloop
          ksave=kloop
       ENDIF
 427  CONTINUE
!
      RMSERR             = SQRT(DER2MAX / ORDER)
!
      L3                 = L3 + 1
!
      IF (L3.GT.1) THEN
       RMSRAT            = RMSERR / RMSERRP    
       DRATE             = MAX(0.2d0 * DRATE, RMSRAT)
      ENDIF
!
      DCON               = RMSERR * MIN(CONPST(NQQ),CONP15(NQQ)*DRATE)
!
! *********************************************************************
!       IF CONVERGENCE OCCURS, GO ON TO CHECK ACCUMULATED ERROR
! *********************************************************************
!


      ! EXPLANATORY NOTE (tmf, 11/1/07)
      ! This is where we check for local error convergence
      IF (DCON .LE. 1.0) THEN

       ! EXPLANATORY NOTE (tmf, 11/1/07)
       ! Go on to check global error
       GOTO 390 
!
! *********************************************************************
!  IF NONCONVERGENCE AFTER ONE STEP, RE-EVALUATE FIRST DERIVATIVE WITH
!                           NEW VALUES OF CNEW
! *********************************************************************
!
!     ELSEIF (L3.LT.MSTEP.AND.(L3.EQ.1.OR.RMSRAT.LE.0.9)) THEN
      ELSEIF (L3.EQ.1) THEN

       GOTO 270
!
! *********************************************************************
! *             THE CORRECTOR ITERATION FAILED TO CONVERGE            *
! * IF THE JACOBIAN MATRIX IS MORE THAN ONE STEP OLD, UPDATE THE      *
! * JACOBIAN AND TRY CONVERGENCE AGAIN. IF THE JACOBIAN IS CURRENT,   * 
! * THEN REDUCE THE TIME-STEP, RE-SET THE ACCUMULATED DERIVATIVES TO  *
! * THEIR VALUES BEFORE THE FAILED STEP, AND RETRY WITH THE SMALLER   *
! * STEP.                                                             *
! *********************************************************************
!

      ELSEIF (JEVAL .EQ. 0) THEN

       ! EXPLANATORY NOTE (tmf, 11/1/07)
       ! If Jacobian is old, reevaluate Jacobian with 
       !   CNEW = CNEW_old + dt*GLOSS(CNEW_old), 
       ! and reset L3 = 0.
       IFAIL           = IFAIL + 1
       JEVAL           = 1

       GOTO 220
      ENDIF

!
 399  NFAIL            = NFAIL + 1
      RDELMAX          = 2.0d0
      JEVAL            = 1
      IFSUCCESS        = 0
      XELAPS           = TOLD
      RDELT            = FRACDEC
!
! *********************************************************************
!               SUBTRACT OFF DERIVATIVES PREVIOUSLY ADDED 
! *********************************************************************
! THIS SET OF OPERATIONS IS EQUIVALENT TO LOOP 419.
!
      ! EXPLANATORY NOTE (tmf, 11/1/07)
      ! The local error is still not converging
      IF (NQQ.EQ.1) THEN
       DO 376 I        = 1, ISCHAN
        J              = I + ISCHAN
        DO 375 KLOOP   = 1, KTLOOP
         CONC(KLOOP,I) = CONC(KLOOP,I) - CONC(KLOOP,J)
 375    CONTINUE
 376   CONTINUE
!
      ELSEIF (NQQ.EQ.2) THEN
!
       DO 378 I         = 1, ISCHAN
        J1              = I  + ISCHAN
        J2              = J1 + ISCHAN 
        DO 377 KLOOP    = 1, KTLOOP
         CONC(KLOOP, I) = CONC(KLOOP, I) - CONC(KLOOP,J1)-CONC(KLOOP,J2) 
         CONC(KLOOP,J1) = CONC(KLOOP,J1) - CONC(KLOOP,J2) * 2.d0 
 377    CONTINUE
 378   CONTINUE
!
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
!
      ELSEIF (NQQ.EQ.4) THEN
!
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
!
      ELSEIF (NQQ.EQ.5) THEN
!
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

      GOTO 170
!
! *********************************************************************
! *               THE CORRECTOR ITERATION CONVERGED                   *  
! * SET JEVAL = 0 SO THAT IT DOES NOT NEED TO BE CALLED THE NEXT STEP *
! * IF ALL ELSE GOES WELL. NEXT, TEST THE ACCUMULATED ERROR FROM THE  * 
! *                    CONVERGENCE PROCESS, ABOVE                     * 
! *********************************************************************
!
 390  JEVAL          = 0
!
      IF (L3.GT.1) THEN
       DO 395 KLOOP  = 1, KTLOOP
        DELY(KLOOP)  = 0.d0
 395   CONTINUE
!
       DO 402 JSPC   = 1, ISCHAN     
        !---------------------------------------------------------------------
        ! Added for the ND65 prod/loss diagnostic.  This prevents ND65
        ! prod/loss species from being counted towards the convergence
        ! criteria for the SMVGEAR solver (ljm, bmy, 5/9/03)
        IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN
           DO 400 KLOOP = 1, KTLOOP 
              ERRYMAX     = DTLOS(KLOOP,JSPC) * CHOLD(KLOOP,JSPC)
              DELY(KLOOP) = DELY(KLOOP) + ERRYMAX * ERRYMAX
 400       CONTINUE
        ENDIF
        !---------------------------------------------------------------------
 402   CONTINUE
!
       DER2MAX       = 0.d0
!
       DO 405 KLOOP  = 1, KTLOOP
        IF (DELY(KLOOP).GT.DER2MAX) DER2MAX = DELY(KLOOP)   
 405   CONTINUE 
      ENDIF
!
! *********************************************************************
! *                 THE ACCUMULATED ERROR TEST FAILED                 *  
! * IN ALL CASES, RE-SET THE DERIVATIVES TO THEIR VALUES BEFORE THE   *
! * LAST TIME-STEP. NEXT                                              *    
! *                                                                   * 
! * (A)   RE-ESTIMATE A TIME-STEP AT THE SAME OR ONE LOWER ORDER AND  *   
! *        RETRY THE STEP.                                            *  
! * (B)   IF THE FIRST ATTEMPTS FAIL, RETRY THE STEP AT FRACDEC x     * 
! *        THE PRIOR STEP                                             *  
! * (C)   IF THIS FAILS, RE-SET THE ORDER TO 1 AND GO BACK TO THE     * 
! *        BEGINNING, AT ORDER = 1, BECAUSE ERRORS OF THE WRONG ORDER *
! *        HAVE ACCUMULATED                                           *
! *********************************************************************
! 
      ! EXPLANATORY NOTE (tmf, 11/1/07)
      ! Check if global error is tolerable
      IF (DER2MAX.GT.ENQQ) THEN

       XELAPS           = TOLD
       LFAIL            = LFAIL + 1
       JFAIL            = JFAIL  + 1
!
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
!
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
!
        DO 432 JSPC        = 1, ISCHAN
         DO 430 KLOOP      = 1, KTLOOP
          CNEW(KLOOP,JSPC) = CONC(KLOOP,JSPC)
 430     CONTINUE
 432    CONTINUE
!
        WRITE(6,670) DELT, XELAPS
        !print*, 'kblk = ', kblk      !gcc
        IF (JRESTAR.EQ.100) THEN
         WRITE(6,680)
         CALL GEOS_CHEM_STOP
        ENDIF
!
        GOTO 140
       ENDIF    ! end of JFAIL condition
!
      ELSE

!
! *********************************************************************
! *             ALL SUCCESSFUL STEPS COME THROUGH HERE                * 
! *                                                                   * 
! * AFTER A SUCCESSFUL STEP, UPDATE THE CONCENTRATION AND ALL DERIV-  *
! * ATIVES, RESET TOLD, SET IFSUCCESS = 1, INCREMENT NSTEPS, AND      *
! * RESET JFAIL = 0.                                                  *
! *********************************************************************
!
!-----------------------------------------------------------------------------
! %%%%% MODIFICATION TO PREVENT NEGATIVE CNEW (tmf, 11/1/07) %%%%%
!
! Even if both local and global error tests are passed, if any of the values 
! in the CNEW array are less than zero, start over with a smaller time step, 
! and re-evaluate the Jacobian (i.e. GOTO 140).
!
      ! Initialize counter
 820  INEG = 0

      ! Count the # of negatives in CNEW
      DO JSPC  = 1, ISCHAN
      DO KLOOP = 1, KTLOOP
         IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN
            IF ( CNEW(KLOOP,JSPC) < 0.d0 ) INEG = INEG + 1
         ENDIF
      ENDDO
      ENDDO
       
      ! If there are negatives in CNEW then ...
      IF ( INEG > 0 ) THEN

         ! ... Set CNEW to CPREVM.  CPREVM stores the values that were in
         ! in the CNEW array at the end of the previous successful march.
         DO JSPC  = 1, ISCHAN
         DO KLOOP = 1, KTLOOP
            CNEW(KLOOP,JSPC) = CPREVM(KLOOP,JSPC)
         ENDDO
         ENDDO

         ! ... Then reset the various time & tolerance variables 
         XELAPS          = TOLD
         JEVAL           = 1
         IFSUCCESS       = 0
         DTFORCE         = DELT * FRACDEC
         IDTFORCE        = 1
         RDELT           = 1.d0
         JFAIL           = 0
         JRESTAR         = JRESTAR + 1
         IDOUB           = 5
         L3 = 0

         ! ... Then try another march
         GOTO 140

      ENDIF
!------------------------------------------------------------------------------

       JFAIL            = 0
       IFSUCCESS        = 1
       NSTEPS           = NSTEPS + 1
       TOLD             = XELAPS
!
! *********************************************************************
!
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
!
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
!
! *********************************************************************
!          EXIT SMVGEAR IF A TIME INTERVAL HAS BEEN COMPLETED
! *********************************************************************
!
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
!
! *********************************************************************
! * IDOUB COUNTS THE NUMBER OF SUCCESSFUL STEPS BEFORE RE-TESTING THE *
! * STEP-SIZE AND ORDER                                               *
! *                                                                   * 
! * IF IDOUB > 1, DECREASE IDOUB AND GO ON TO THE NEXT TIME-STEP WITH *
! *               THE CURRENT STEP-SIZE AND ORDER.                    * 
! * IF IDOUB = 1, STORE THE VALUE OF THE ERROR (DTLOS) FOR THE TIME-  *
! *               STEP PREDICTION, WHICH WILL OCCUR WHEN IDOUB = 0,   *
! *               BUT GO ON TO THE NEXT STEP WITH THE CURRENT STEP-   *
! *               SIZE AND ORDER.                                     * 
! * IF IDOUB = 0, TEST THE TIME-STEP AND ORDER FOR A CHANGE.          *  
! *********************************************************************
!
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

        !---------------------------------------------------------------------
        ! %%%%% MODIFICATION TO PREVENT NEGATIVE CNEW (tmf, 11/1/07) %%%%%
        !
        ! Store successful CNEW in CPREVM.  CPREVM will be used in a future
        ! march to re-initialize the CNEW array if negatives are found.
        !
        DO JSPC  = 1, ISCHAN 
        DO KLOOP = 1, KTLOOP
           CPREVM(KLOOP,JSPC) = CNEW(KLOOP,JSPC)
        ENDDO
        ENDDO
        !---------------------------------------------------------------------

        GOTO 170
       ENDIF
!
      ENDIF 
!     ENDIF DER2MAX.GT.ENQQ
!
! *********************************************************************
! *         TEST WHETHER TO CHANGE THE STEP-SIZE AND ORDER            * 
! * DETERMINE THE TIME-STEP AT (A) ONE ORDER LOWER THAN, (B) THE SAME *
! * ORDER AS, AND (C) ONE ORDER HIGHER THAN THE CURRENT ORDER. IN THE *
! * CASE OF MULTIPLE GRID-CELLS IN A GRID-BLOCK, FIND THE MINIMUM     *
! * STEP-SIZE AMONG ALL THE CELLS FOR EACH OF THE ORDERS. THEN, IN    *
! * ALL CASES, CHOOSE THE LONGEST TIME-STEP AMONG THE THREE STEPS     *
! * PAIRED WITH ORDERS, AND CHOOSE THE ORDER ALLOWING THIS LONGEST    *
! * STEP.                                                             *
! *********************************************************************
!
! *********************************************************************
! * ESTIMATE THE TIME-STEP RATIO (RDELTUP) AT ONE ORDER HIGHER THAN   *
! * THE CURRENT ORDER. IF NQQ >= MAXORD, THEN WE DO NOT ALLOW THE     * 
! *                        ORDER TO INCREASE.                         *
! *********************************************************************
!
      IF (NQQ.LT.MAXORD) THEN  
       DO 542 KLOOP  = 1, KTLOOP
        DELY(KLOOP)  = 0.d0
 542   CONTINUE
!
       DO 545 JSPC   = 1, ISCHAN     
        !---------------------------------------------------------------------
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
        !---------------------------------------------------------------------
 545   CONTINUE
!
       DER3MAX       = 0.d0
!
       DO 546 KLOOP  = 1, KTLOOP
        IF (DELY(KLOOP).GT.DER3MAX) DER3MAX = DELY(KLOOP)
 546   CONTINUE
!
       RDELTUP       = 1.0d0 / (CONP3*DER3MAX**ENQQ3(NQQ)+1.4d-6)
      ELSE
       RDELTUP       = 0.0d0
      ENDIF
!
! *********************************************************************
! *    ESTIMATE THE TIME-STEP RATIO (RDELTSM) AT THE CURRENT ORDER    *
! *      WE CALCULATED DER2MAX DURING THE ERROR TESTS EARLIER         *  
! *********************************************************************
!
 540  RDELTSM        = 1.0d0 / (CONP2*DER2MAX**ENQQ2(NQQ)+1.2d-6)
!
! *********************************************************************
! * ESTIMATE THE TIME-STEP RATIO (RDELTDN) AT ONE ORDER LOWER THAN    *
! * THE CURRENT ORDER. IF NQQ = 1, THEN WE CANNOT TEST A LOWER ORDER. * 
! *********************************************************************
!
      IF (NQQ.GT.1) THEN
       DO 552 KLOOP  = 1, KTLOOP
        DELY(KLOOP)  = 0.d0
 552   CONTINUE
!
       KSTEPISC      = (KSTEP - 1) * ISCHAN
       DO 555 JSPC   = 1, ISCHAN     
        !---------------------------------------------------------------------
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
        !---------------------------------------------------------------------
 555   CONTINUE
!
       DER1MAX       = 0.d0
!
       DO 556 KLOOP = 1, KTLOOP
        IF (DELY(KLOOP).GT.DER1MAX) DER1MAX = DELY(KLOOP)   
 556   CONTINUE
        
       RDELTDN       = 1.0d0 / (CONP1*DER1MAX**ENQQ1(NQQ)+1.3d-6)
!
      ELSE
       RDELTDN       = 0.d0
      ENDIF
!
! *********************************************************************
! * FIND THE LARGEST OF THE PREDICTED TIME-STEPS RATIOS OF EACH ORDER * 
! *********************************************************************
!
      RDELT        = MAX(RDELTUP,RDELTSM,RDELTDN)
!
! *********************************************************************
! * IF THE LAST STEP WAS SUCCESSFUL AND RDELT IS SMALL, KEEP THE      *
! * CURRENT STEP AND ORDER AND ALLOW THREE SUCCESSFUL STEPS BEFORE    *
! * RE-CHECKING THE TIME-STEP AND ORDER.                              *
! *********************************************************************
!
      IF (RDELT.LT.1.1.AND.IFSUCCESS.EQ.1) THEN
       IDOUB         = 3
       GOTO 170
!
! *********************************************************************
! * IF THE MAXIMUM TIME-STEP RATIO IS THAT OF ONE ORDER LOWER THAN    *  
! * THE CURRENT ORDER, DECREASE THE ORDER. DO NOT MINIMIZE RDELT      *
! * TO =< 1 WHEN IFSUCCESS = 0 SINCE THIS IS LESS EFFICIENT.          *
! *********************************************************************
!
      ELSEIF (RDELT.EQ.RDELTDN) THEN
       NQQ              = NQQ - 1

!
! *********************************************************************
! * IF THE MAXIMUM TIME-STEP RATIO IS THAT OF ONE ORDER HIGHER THAN   *  
! * THE CURRENT ORDER, INCREASE THE ORDER AND ADD A DERIVATIVE TERM   *
! * FOR THE HIGHER ORDER.                                             *  
! *********************************************************************
!
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
!
! *********************************************************************
! * IF THE LAST TWO STEPS HAVE FAILED, RE-SET IDOUB TO THE CURRENT    *
! * ORDER + 1. DO NOT MINIMIZE RDELT IF JFAIL.GE.2 SINCE TESTS SHOW   *
! * THAT THIS MERELY LEADS TO ADDITIONAL COMPUTATIONS.                *
! *********************************************************************
!
      IDOUB             = NQQ + 1
!

      GOTO 170  
!
! *********************************************************************
! *                      UPDATE COUNTERS                              *
! *********************************************************************
!
 650  NSFTOT             = NSFTOT   + NSUBFUN 
      NPDTOT             = NPDTOT   + NPDERIV
      NSTTOT             = NSTTOT   + NSTEPS 
      IFAILTOT           = IFAILTOT + IFAIL
      NFAILTOT           = NFAILTOT + NFAIL
      LFAILTOT           = LFAILTOT + LFAIL

!
! *********************************************************************
! *        SET FINAL CONCENTRATION FOR RUN AND UPDATE COUNTERS        * 
! *********************************************************************
!
      DO JSPC        = 1, ISCHAN
      DO KLOOP       = 1, KTLOOP

         ! Stop with an error message if NaN's are encountered
         ! (bmy, pip, 4/27/00)
         IF ( IT_IS_NAN( CNEW(KLOOP,JSPC) ) ) THEN
            DO NK = 1, NALLRAT(NCS)
               WRITE( 6, 249 ) NK, RRATE(KLOOP,NK),TRATE(KLOOP,NK),KLOOP
 249           FORMAT( 'Rxn #:', i5, ' RRATE = ', es13.6 ,
     x              '  TRATE = ', es13.6, 'KLOOP = ', I6 )
            ENDDO
            write(6,*) 'sum of rrate = ',sum(rrate)
            PRINT*, 'SMVGEAR: CNEW is NaN!'
            PRINT*, 'Species index : ', JSPC
            PRINT*, 'Species NAME : ', NAMEGAS(JSPC)

            PRINT*, 'Grid Box      : ', IXSAVE(KLOOP+JLOOPLO),
     &           IYSAVE(KLOOP+JLOOPLO), IZSAVE(KLOOP+JLOOPLO)
            PRINT*, 'STOP AT END OF smvgear.f!'

            ! Stop the run and deallocate all arrays 
            CALL GEOS_CHEM_STOP
         ENDIF

!------------------------------------------------------------------------------
! %%%%% MODIFICATION TO PREVENT NEGATIVE CNEW (tmf, 11/1/07) %%%%%
!
! The previous code would reset negative CNEW to a small positive #.  We don't
! want to do that anymore.  If any negative CNEW values exist at the end of
! an internal march, we reset the timestep and re-evaluate the Jacobian.
! (see at CONTINUE statement 820).
!
! However, if any negative values in CNEW still exist at this point, we will
! stop the run and print out debug information.


        IF ( ITS_NOT_A_ND65_FAMILY(JSPC) ) THEN

         IF ( CNEW(KLOOP,JSPC) .LT. 0d0 ) THEN
            DO NK = 1, NALLRAT(NCS)
               WRITE( 6, 249 ) NK, RRATE(KSAVE,NK),TRATE(KSAVE,NK)
 250           FORMAT( 'Rxn #:', i5, ' RRATE = ', es13.6 ,
     x              '  TRATE = ', es13.6 )
            ENDDO
            write(6,*) 'sum of rrate = ',sum(rrate)
            PRINT*, 'SMVGEAR: CNEW is negative!'
            PRINT*, 'Species index : ', JSPC, 'CNEW =', CNEW(KLOOP,JSPC)
            PRINT*, 'Grid Box      : ', IXSAVE(KLOOP+JLOOPLO),
     &           IYSAVE(KLOOP+JLOOPLO), IZSAVE(KLOOP+JLOOPLO)
            PRINT*, 'STOP in smvgear.f!'

            ! Stop the run and deallocate all arrays 
            CALL GEOS_CHEM_STOP
         ENDIF

        ENDIF

          ! Comment this line out, we don't want to reset CNEW anymore
          !! Reset negatives to a very small positive number
          !CNEW(KLOOP,JSPC) = MAX(CNEW(KLOOP,JSPC),SMAL2)

!------------------------------------------------------------------------------

      ENDDO
      ENDDO


!
!
! *********************************************************************
!                             FORMATS
! *********************************************************************
!
      ! Also add CPU # to formats 233, 685, 690 (bmy, 3/8/13)
 233  FORMAT('SMVGEAR: DELT= ',1PE8.2,' TOO LOW DEC YFAC. CPU, KBLK, ', 
     1       'KTLOOP, NCS, TIME, TIMREMAIN, YFAC, ', 
     2       'EPS = ',/4(1X,I4),2X,4(1PE9.3,1X))
 234  FORMAT('SMVGEAR: TOO MANY DECREASES OF YFAC ')
 670  FORMAT('DELT DEC TO =',E13.5,'; TIME ',E13.5,' BECAUSE ',
     1       'EXCESSIVE ERRORS')
 680  FORMAT('SMVGEAR: STOP BECAUSE OF EXCESSIVE ERRORS.')
 685  FORMAT('CPU,M1,M2,K,ERR = ',4(I4),2X,1PE10.4)
 690  FORMAT('CPU,CONC WHEN STOP = ',3(I4,1X),A14,2(1X,1PE10.2))
!
! *********************************************************************
! ***************     END OF SUBROUTINE SMVGEAR     *******************
! *********************************************************************
!
      RETURN
      END SUBROUTINE SMVGEAR
