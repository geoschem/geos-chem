!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: physproc
!
! !DESCRIPTION: Subroutine PHYSPROC is the driver for SMVGEAR II chemistry.  
!  It calls both CALCRATE to compute the rxn rates and the SMVGEAR solver 
!  routine.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PHYSPROC( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,         ONLY : ABSHUM
      USE COMODE_MOD,         ONLY : AIRDENS
      USE COMODE_MOD,         ONLY : CSPEC
      USE COMODE_MOD,         ONLY : CSUMA
      USE COMODE_MOD,         ONLY : CSUMC
      USE COMODE_MOD,         ONLY : ERRMX2
      USE COMODE_MOD,         ONLY : IXSAVE
      USE COMODE_MOD,         ONLY : IYSAVE
      USE COMODE_MOD,         ONLY : T3
      USE COMODE_LOOP_MOD
      USE CHEMISTRY_MOD,      ONLY : GCKPP_DRIVER
      USE GCKPP_GLOBAL,       ONLY : NTT
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : TIMESTAMP_STRING

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! *********************************************************************
! ************        WRITTEN BY MARK JACOBSON (1993)      ************
! ***             (C) COPYRIGHT, 1993 BY MARK Z. JACOBSON           *** 
! ***       U.S. COPYRIGHT OFFICE REGISTRATION NO. TXu 670-279      *** 
! ***                         (650) 723-6836                        *** 
! *********************************************************************
!
! PPPPPPP  H     H  Y     Y  SSSSSSS  PPPPPPP RRRRRRR  OOOOOOO  CCCCCCC
! P     P  H     H   Y   Y   S        P     P R     R  O     O  C 
! PPPPPPP  HHHHHHH     Y     SSSSSSS  PPPPPPP RRRRRRR  O     O  C 
! P        H     H     Y           S  P       R  R     O     O  C 
! P        H     H     Y     SSSSSSS  P       P    P   OOOOOOO  CCCCCCC
!
! *********************************************************************
! * THIS SUBROUTINE CALLS CALCRATE.F AND SMVGEAR.F. TO SOLVE GAS-     *
! * PHASE CHEMICAL EQUATIONS. THE ROUTINE DIVIDES THE GRID DOMAIN     *
! * INTO GRID BLOCKS, AND THE CODE VECTORIZES AROUND THE NUMBER OF    *
! * GRID CELLS IN EACH BLOCK.                                         *
! *                                                                   *
! *********************************************************************
!                                                                             .
! *********************************************************************
! ****************         UPDATE 24-HOUR CLOCK        ****************
! *********************************************************************
! CHEMINTV = TIME INTERVAL FOR CHEMISTRY
! IRCHEM   = COUNTS # CHEMINTV TIME-INTERVALS
! 
! !REVISION HISTORY: 
!  03 Jan 1993 - M. Jacobson - Initial version
!  (1 ) For GEOS-CHEM we had to remove ABSHUM, AIRDENS, CSPEC, IXSAVE, IYSAVE,
!        and T3 from "comode.h" and to declare these allocatable in F90 module
!        "comode_mod.f".  This allows us to only allocate these if we are 
!        doing a fullchem run.  Now references TIMESTAMP_STRING from 
!        "time_mod.f".  Now pass SUNCOSB via the arg list.  Now force double
!        precision with the "D" exponent. (bmy, 4/18/03)
!  (2 ) Comment out section that computes photorates from original SMVGEAR II
!        file "photrate.dat"...this is not needed.  Remove TFROMID, it's not
!        used anywhere else.  Remove references to LASTCHEM, this is mpt 
!        initialized anywhere.  Now reference CSUMA, CSUMC, ERRMX2 from
!        "comode_mod.f". (bmy, 7/30/03)
!  (3 ) LINUX has a problem putting a function call w/in a WRITE statement.  
!        Now save output from TIMESTAMP_STRING to STAMP and print that.
!        (bmy, 9/29/03)
!  (4 ) Fixed case of small KULOOP (phs, 10/5/07)
!  (5 ) Now only get the rx rates if not using SMVGEAR (phs,ks,dhk, 09/15/09)
!  (6 ) Now calls KPP after calculating the reaction rates to save memory.
!       (ccc, 12/9/09)
!  (7 ) Remove obsolete print statements & formats (bmy, 12/18/09)
!  (8 ) Now call GCKPP_DRIVER with NSPEC(1), which is the # of active species 
!        for urban chemistry. (ccc, 1/20/10)! 
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!  27 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOSmid
!  01 Mar 2013 - R. Yantosca - Now set RC to GIGC_SUCCESS for default
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER IDLAST,IMLAST,NMINADD,NHRADD,NDAYADD,NMONADD,NYEARAD
      INTEGER MLOOP,LOREORD,ILNCS,IHNCS,JLOOP,NBLOCKUSE,IAVBLOK,IAVGSIZE
      INTEGER JLOOPN,JOLD,JNEW,KLOOP,NSUNRISE,NSUNSET,JLOOPC,NNORISE
      INTEGER JLLAST,IT,IRADD,LVAL,IRVAL,IRADD1,JREORD,IPAR,JPAR,JPAR1
      INTEGER NSUMBLOCK,NCELLROW,NBLOCKROW,ICG,I,NGCOUNT,NGHI,IAVG
      INTEGER IREMAIN,IUSESIZE,NREBLOCK,L
      INTEGER IX,IY,IJWINDOW,KBLK2,NK

      INTEGER COUNTER,JGAS
      REAL*8 S1CON,S2CON,ARGS,CONSTQ,SNOON,CONTEMP,DIFCONC,PLODYN
      REAL*8 PR3,RIS,RST,TBEGIN,TFINISH,VALLOW,CLO1,CLO2,SUMFRACS
      REAL*8 SUMRMS,SUMHI,SUMRMSH,CMOD,CGOOD,FRACDIF,FRACABS,AVGERR
      REAL*8 RMSCUR,AVGHI,RMSCURH,FSTEPT,FITS,TSTEPIT,PHIDYN
      REAL*8 GMU

      INTEGER :: AS
      ! For LINUX fix (bmy, 9/29/03)
      CHARACTER(LEN=16) :: STAMP

      REAL :: START,FINISH

      LOGICAL :: LKPP

      !=================================================================
      ! PHYSPROC begins here!
      !=================================================================

      ! Assume success
      RC      = GIGC_SUCCESS

      ! Copy values from Input_Opt to local variables
      LKPP    = Input_Opt%LKPP

      ! Initialize
      IRCHEM  = IRCHEM + 1
      TIME    = TIME        + CHEMINTV
      TSPMIDC = MOD(TSPMIDC + CHEMINTV,SCDAY) 
      COUNTER = 0

      ! Return if we have turned off SMVGEAR
      IF (IFSOLVE.EQ.0) RETURN

      ! Echo timestamp
      STAMP = TIMESTAMP_STRING()
      IF ( am_I_Root ) WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - PHYSPROC: Trop chemistry at ', a )
!                                                                       
! *********************************************************************
! *********************************************************************
! NCS       = 1..NCSGAS --> DO GAS CHEMISTRY
! LOREORD   = 1 IF REORDERING; = 2 IF NO REORDERING
!


      IF (IFREORD.EQ.1.AND.NTLOOP.GT.1) THEN
       LOREORD   = 1
      ELSE 
       LOREORD   = 2
      ENDIF 
!
      ILNCS      = 1
      IHNCS      = NCSGAS 
!
! *********************************************************************
! *         REORDER CELLS AND BLOCKS THEN SOLVE CHEMICAL ODES         *
! *********************************************************************
! ISREORD   = 1: THEN REORDER GRID CELLS AND GRID BLOCKS FOR CHEMISTY;
!           = 2: SOLVE CHEMISTRY 
! JREORDER  = GIVES ORIGINAL  GRID-CELL FROM RE-ORDERED GRID-CELL
! LREORDER  = JREORDER
! NBLOCKUSE = # OF ORIGINAL BLOCKS (ISREORD EQ LOREORD) OR
!             # OF BLOCKS AFTER REORDERING (ISREORD NE LOREORD)
! NCS       = 1..NCSGAS FOR GAS CHEMISTRY
! NCSP      = NCS       FOR DAYTIME   GAS CHEM
!           = NCS + ICS FOR NIGHTTIME GAS CHEM
!
      DO 860 NCS            = ILNCS, IHNCS
!
       DO 855 ISREORD       = LOREORD, 2
!
        !write(6,*) 'value of isreord= ',isreord
        IF (ISREORD.EQ.LOREORD) THEN
!
! *********************************************************************
!                   DETERMINE BLOCK SIZES FOR CHEMISTRY
! *********************************************************************
!                CHEMISTRY IN ONE REGION OF THE ATMOSPHERE
! *********************************************************************
! IGLOBCHEM = -2 --> SOLVE ALL GAS CHEMISTRY WITH COMBINATION OF U/R/S SETS
!           = -1 --> SOLVE ALL GAS CHEMISTRY WITH COMBINATION OF R/S SETS
!           = 0  --> SOLVE ALL GAS CHEMISTRY WITH EITHER U, R, OR S SETS
!           = 1  --> SOLVE EACH REGION SEPARATELY WITH U, R, OR S SET
!
          IF (IGLOBCHEM.LE.0) THEN
            !NTLOOPUSE           = NTLOOPNCS(NCS)

           ! updated ntloop calc in ruralbox.f
            NTLOOPUSE          = NTLOOP  
           DO 320 JLOOP        = 1, NTLOOPUSE
 320        JREORDER(JLOOP)    = JLOOP
!
          ELSE

!
! *********************************************************************
!        GLOBAL CHEMISTRY - ASSUME THREE REGIONS OF THE ATMOSPHERE
!                   URBAN, TROPOSPHERIC, STRATOSPHERIC
! *********************************************************************
! NCS     = 1..NCSGAS FOR GAS CHEMISTRY
! PRESS3  = MODEL VERTICAL LAYER CENTER PRESSURE (MB)
! PLOURB  = PRES (MB), BELOW WHICH URBAN, URBAN/TROP, OR ALL CHEMISTRY OCCURS
! PLOTROP = PRES (MB), BELOW WHICH TROP,  URBAN/TROP, OR ALL CHEMISTRY OCCURS
!         =            ABOVE WHICH STRAT              OR ALL CHEMISTRY OCCURS
!
             IF (NCS.EQ.NCSURBAN) THEN
                NTLOOPUSE          = NTLOOPNCS(NCS)
                DO JLOOP        = 1, NTLOOPUSE
                   JREORDER(JLOOP)    = NCSLOOP(JLOOP,NCS)
                ENDDO
             ELSEIF (NCS.EQ.NCSTROP) THEN
                NTLOOPUSE          = NTLOOPNCS(NCS)
                DO JLOOP        = 1, NTLOOPUSE
                   JREORDER(JLOOP)    = NCSLOOP(JLOOP,NCS)
                ENDDO
             ELSEIF (NCS.EQ.NCSSTRAT) THEN
                NTLOOPUSE          = NTLOOPNCS(NCS)
                DO JLOOP        = 1, NTLOOPUSE
                   JREORDER(JLOOP)    = NCSLOOP(JLOOP,NCS)
                ENDDO
             ENDIF
!
          ENDIF
!         ENDIF IGLOBCHEM.EQ.0
!
! *********************************************************************
!             DETERMINE ORIGINAL NUMBER OF GRID BLOCKS
! *********************************************************************
! NBLOCKUSE = ORIGINAL NUMBER OF GRID BLOCKS FOR PREDICTING STIFFNESS
! IUSESIZE  = # OF GRID CELLS IN EACH GRID BLOCK
! NBLOCKUSE = HERE, TOTAL NUMBER OF GRID CELLS FOR CHEMISTRY CALCULATIONS
! JLOWVAR   = LOWEST GRID CELL NUMBER - 1 IN EACH GRID BLOCK
!
          ! Comment out write statements for now (bmy, 4/1/03)
          !write(6,*) 'in physproc, iglobchem= ',iglobchem
          !write(6,*) 'val of ntloopuse= ',ntloopuse

          NBLOCKUSE          = 1 + NTLOOPUSE / (KULOOP    + 0.0001d0)
          IAVBLOK            = 1 + NTLOOPUSE / (NBLOCKUSE + 0.0001d0)
          IAVGSIZE           = MIN0(IAVBLOK,KULOOP)
!
          JLOOPLO            = 0
          IREMAIN            = NTLOOPUSE
!
          DO 200 KBLK        = 1, NBLOCKUSE
           IUSESIZE          = MIN(IAVGSIZE,MAX(IREMAIN,0))
           JLOWVAR(KBLK)     = JLOOPLO
           KTLPVAR(KBLK)     = IUSESIZE
           IREMAIN           = IREMAIN - IUSESIZE
           JLOOPLO           = JLOOPLO + IUSESIZE
 200      CONTINUE

          ! Added fix for small (1 to 3) KULOOP (10/5/07, phs)
          IF (IREMAIN /= 0) THEN
             DO WHILE ( IREMAIN /= 0 )
                NBLOCKUSE          = NBLOCKUSE + 1
                IUSESIZE           = MIN(IAVGSIZE,MAX(IREMAIN,0))
                JLOWVAR(NBLOCKUSE) = JLOOPLO
                KTLPVAR(NBLOCKUSE) = IUSESIZE
                IREMAIN            = IREMAIN - IUSESIZE
                JLOOPLO            = JLOOPLO + IUSESIZE
             END DO
          ENDIF

!
! *********************************************************************
!                  NUMBER OF GRID BLOCKS AFTER REORDERING 
! *********************************************************************
!
         ELSE 
          NBLOCKUSE          = NREBLOCK
         ENDIF
!        ENDIF ISREORD.EQ.LOREORD
!
! *********************************************************************
!                          SET LREORDER ARRAY
! *********************************************************************
! LREORDER = GIVES ORIGINAL GRID CELL FROM RE-ORDERED CELL
!

         DO 340 JLOOPN      = 1, NTLOOPUSE
 340      LREORDER(JLOOPN)  = JREORDER(JLOOPN)
!
! *********************************************************************
!                   START GRID BLOCK LOOP
! *********************************************************************
!
!--- Moved from chemdr.f (ccc, 12/9/09)
          IF ( LKPP) NTT = NTTLOOP

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JLOOP,KLOOP,KBLK2,JNEW,JOLD)
!$OMP+SCHEDULE( DYNAMIC )

         DO 640 KBLK2       = 1, NBLOCKUSE
          KBLK              = KBLK2
          JLOOPLO           = JLOWVAR(KBLK)
          KTLOOP            = KTLPVAR(KBLK)
!

          IF (KTLOOP.EQ.0) GOTO 640

! *********************************************************************
! *  PLACE LARGE DOMAIN GAS ARRAY (# CM-3) INTO SMALLER BLOCK ARRAY   *
! *********************************************************************
! CINIT  = INITIAL CONC (MOLEC. CM-3), USED TO CALCULATE RATES IN CALCRATE
! CORIG  = INITIAL CONC (MOLEC. CM-3), USED TO RESTART SMVGEAR IF A FAILURE
! NTSPEC = NUMBER OF ACTIVE PLUS INACTIVE GASES.
! MAPPL  = MAPS ORIGINAL SPECIES NUMBERS TO SPECIES NUMBERS
!          RE-ORDERED FOR CHEMISTRY.
!
             DO 572 JOLD         = 1, NTSPEC(NCS)
                JNEW               = MAPPL(JOLD,NCS)
                DO 570 KLOOP       = 1, KTLOOP
                   JLOOP             = JREORDER(JLOOPLO+KLOOP)
                   CBLK( KLOOP,JOLD) = CSPEC(JLOOP,JOLD)
                   CORIG(KLOOP,JNEW) = CSPEC(JLOOP,JOLD)
 570            CONTINUE
 572         CONTINUE
!
! *********************************************************************
! *              CALCULATE REACTION RATE COEFFICIENTS                 *
! *********************************************************************
!
          CALL CALCRATE( am_I_Root, Input_Opt, 
     &                   State_Met, State_Chm, RC )
!
! *********************************************************************
! *                        SOLVE CHEMICAL ODES                        *
! *********************************************************************
!
!--- Move call to KPP here from chemdr.f to save memory space
!    (ccc, 12/9/09)
!    NSPEC(1) is the # of active species for urban chemistry.
!    (ccc, 01/20/10)
          call CPU_TIME(start)
          IF ( LKPP ) THEN
             CALL GCKPP_DRIVER(KTLOOP, JLOOPLO, RRATE_FOR_KPP, NSPEC(1))

          ELSE
             CALL SMVGEAR( am_I_Root, Input_Opt, RC )
          call CPU_TIME(finish)

!
! *********************************************************************
! * REPLACE BLOCK CONCENTRATIONS (# CM-3) INTO DOMAIN CONCENTRATIONS  *
! *********************************************************************
! ISREORD = 1: CALL CALCRATE TO FIND STIFFNESS OF EACH GRID-CELL     
! ISREORD = 2: SET CHEMISTRY RATES AND SOLVE EQUATIONS                
! GQSCHEM   = CHANGE IN (# OF MOLECULES) OVER THE ENTIRE GRID OF EACH 
!             SPECIES AS A RESULT OF CHEMISTRY.
! GRIDVH   = GRID CELL VOLUME (CM3) 
! CNEW     = # CM-3
! C        = # CM-3 
!
             IF (ISREORD.EQ.2) THEN
                DO 620 JNEW         = 1, ISCHANG(NCS)
                   JOLD               = INEWOLD(JNEW,NCS)
                   DO 620 KLOOP       = 1, KTLOOP
                      JLOOP             = JREORDER(JLOOPLO+KLOOP)
                      CSPEC(JLOOP,JOLD) = MAX(CNEW(KLOOP,JNEW),SMAL2)
 620            CONTINUE
             ENDIF
          ENDIF
!
          
 640      CONTINUE
!$OMP END PARALLEL DO

!------------------------------------------------------------------------------
! Prior to 4/28/10:
! Comment out code for reordering grid boxes.  We don't use this as it is
! not computationally efficient for our situation. (bmy, 4/28/10)
!C        CONTINUE KBLK = 1, NBLOCKUSE
!C
!C *********************************************************************
!C *           REORDER GRID-CELLS FROM LEAST TO MOST STIFF             *
!C *********************************************************************
!C AT SUNRISE/SET REORDER CELLS SO THOSE WITH SIMILAR SUNRISE GROUP TOGETHER
!C OTHERWISE,  REORDER CELLS SO THOSE WITH SIMILAR STIFFNESS GROUP TOGETHER
!C JREORDER  = GIVES ORIGINAL  GRID-CELL FROM RE-ORDERED GRID-CELL
!C LREORDER  = GIVES ORIGINAL GRID CELL FROM RE-ORDERED CELL, EXCEPT, 
!C             WHEN CELL IS A VIRTUAL BOUNDARY CELL, THEN LREORDER 
!C             GIVES ORIGINAL EDGE CELL FROM RE-ORDERED V. B. CELL 
!C JLOOPC    = IDENTIFIES AN EDGE CELL FOR EACH VIRTUAL BOUNDARY CELL; 
!C             OTHERWISE, IDENTIFIES THE JLOOP CELL 
!C
!
!         !==============================================================
!         ! New additions to reordering based on Loretta's implementation
!         ! (bdf, 4/1/03)
!         !==============================================================
!         IF (ISREORD.EQ.1) THEN
!            NSUNRISE = 0
!            NSUNSET  = 0
!
!            DO 660 JLOOP = 1, NTLOOPUSE
!               JLOOPC    = LREORDER(JLOOP)
!               IX        = IXSAVE(JLOOPC)
!               IY        = IYSAVE(JLOOPC)
!               IJWINDOW  = (IY-1)*IIPAR + IX
!
!               IF( State_Met%SUNCOSmid(IX,IY) .GT. -.25 .AND.
!     &             State_Met%SUNCOSmid(IX,IY) .LT. .25 ) THEN
!                  ITWO(JLOOP)  = 1
!                  NSUNRISE     = NSUNRISE + 1
!                  CSUMA(JLOOP) = SUNCOS(IJWINDOW) - 
!     &                           ABS( SUNCOSB(IJWINDOW) )
!               ELSE
!                  ITWO(JLOOP)  = 0
!                  CSUMA(JLOOP) = ERRMX2(JLOOP)
!               ENDIF
! 660        CONTINUE
!
!            NNORISE = NTLOOPUSE - NSUNRISE - NSUNSET
!
!            DO 670 JLOOP = 1, NTLOOPUSE
!               LREORDER(JLOOP) = JREORDER(JLOOP)
!               CSUMC(   JLOOP) = CSUMA(   JLOOP)
! 670        CONTINUE
!C
!C *********************************************************************
!C    REORDER GRID-CELLS SO ALL CELLS WHERE SUNSET OCCURS ARE AT END
!C         LREORDER AND CSUMC ARE USED HERE ONLY TO STORE VALUES
!C        OF JREORDER AND CSUMA TEMPORARILY AND ARE USED ELSEWHERE
!C                       FOR A DIFFERENT PURPOSE. 
!C *********************************************************************
!C
!          JLLAST             = NTLOOPUSE 
!          DO 700 JLOOP       = 1, NTLOOPUSE
!           IF (ITWO(JLOOP).EQ.2) THEN
!            JREORDER(JLLAST) = LREORDER(JLOOP) 
!            CSUMA(   JLLAST) = CSUMC(   JLOOP) 
!            JLLAST           = JLLAST - 1 
!           ENDIF
! 700      CONTINUE
!C
!C *********************************************************************
!C    NOW REORDER GRID-CELLS SO ALL CELLS WHERE SUNRISE OCCURS ARE 
!C           IMMEDIATELY BEFORE CELLS WHERE SUNSET OCCURS 
!C *********************************************************************
!C
!          DO 705 JLOOP       = 1, NTLOOPUSE
!           IF (ITWO(JLOOP).EQ.1) THEN
!            JREORDER(JLLAST) = LREORDER(JLOOP) 
!            CSUMA(   JLLAST) = CSUMC(   JLOOP) 
!            JLLAST           = JLLAST - 1 
!           ENDIF
! 705      CONTINUE
!C
!C *********************************************************************
!C FINALLY, PLACE ALL OTHER GRID CELLS BEFORE SUNRISE AND SUNSET CELLS.
!C              JLLAST WILL EQUAL ZERO AFTER THIS LOOP 
!C *********************************************************************
!C
!          DO 710 JLOOP       = 1, NTLOOPUSE
!           IF (ITWO(JLOOP).EQ.0) THEN
!            JREORDER(JLLAST) = LREORDER(JLOOP) 
!            CSUMA(   JLLAST) = CSUMC(   JLOOP) 
!            JLLAST           = JLLAST - 1 
!           ENDIF
! 710      CONTINUE 
!C
!C *********************************************************************
!C REORDER GRID-CELLS IN THREE STEPS:
!C   1) WHERE NO SUNRISE/SET, FROM LEAST TO MOST STIFF 
!C      (SMALLER ERRMX2 (CSUMA) -->LESS STIFF) 
!C      CSUMA = ERRMX2
!C   2) WHERE SUNRISE OCCURS, FROM TIME OF SUNRISE 
!C      CSUMA = TIME OF SUNRISE (IN SECONDS PAST MIDNIGHT)
!C   3) WHERE SUNSET  OCCURS, FROM TIME OF SUNSET  
!C      CSUMA = TIME OF SUNSET (IN SECONDS PAST MIDNIGHT)
!C
!C SORT USING HEAPSORT ROUTINE (NUMERICAL RECIPES), AN N(logb2)N PROCESS 
!C THIS REORDERING SCHEME IS VERY FAST, ALTHOUGH COMPLICATED.
!C ERRMX2 FROM SMVGEAR: DENOTES STIFFNESS (LARGER VALUE --> MORE STIFF).
!C *********************************************************************
!C
!          DO 760 IT             = 1, 3 
!           IF (IT.EQ.1) THEN
!            IRADD               = 0.d0
!            LVAL                = IRADD + NNORISE  * 0.5d0 + 1
!            IRVAL               = IRADD + NNORISE
!           ELSEIF (IT.EQ.2) THEN
!            IRADD               = NNORISE 
!            LVAL                = IRADD + NSUNRISE * 0.5d0 + 1
!            IRVAL               = IRADD + NSUNRISE
!           ELSEIF (IT.EQ.3) THEN
!            IRADD               = NNORISE + NSUNRISE  
!            LVAL                = IRADD + NSUNSET  * 0.5d0 + 1
!            IRVAL               = IRADD + NSUNSET 
!           ENDIF 
!C
!           IRADD1               = IRADD + 1
!C
!           IF (IRVAL.GT.IRADD1) THEN
!C
! 800        IF (LVAL.GT.IRADD1) THEN
!             LVAL                = LVAL - 1
!             VALLOW              = CSUMA(   LVAL)
!             JREORD              = JREORDER(LVAL)        
!            ELSE 
!             VALLOW              = CSUMA(   IRVAL)
!             JREORD              = JREORDER(IRVAL)        
!             CSUMA(   IRVAL)     = CSUMA(   IRADD1)
!             JREORDER(IRVAL)     = JREORDER(IRADD1) 
!             IRVAL               = IRVAL - 1
!             IF (IRVAL.EQ.IRADD1) THEN
!              CSUMA(    IRADD1)  = VALLOW
!              JREORDER( IRADD1)  = JREORD
!              GOTO 760 
!             ENDIF
!            ENDIF
!            IPAR                 = LVAL
!            JPAR                 = LVAL + LVAL - IRADD 
!C
! 820        IF (JPAR.LE.IRVAL) THEN
!             IF (JPAR.LT.IRVAL) THEN
!              JPAR1              = JPAR + 1
!              IF (CSUMA(JPAR).LT.CSUMA(JPAR1)) JPAR = JPAR1
!             ENDIF
!             IF (VALLOW.LT.CSUMA(JPAR)) THEN
!              CSUMA(   IPAR)     = CSUMA(   JPAR)
!              JREORDER(IPAR)     = JREORDER(JPAR)
!              IPAR               = JPAR
!              JPAR               = JPAR + JPAR - IRADD 
!              GOTO 820 
!             ENDIF
!            ENDIF
!C
!            CSUMA(   IPAR)       = VALLOW
!            JREORDER(IPAR)       = JREORD 
!            GOTO 800  
!C
!           ENDIF
!C          ENDIF IRVAL.GT.0
! 760      CONTINUE 
!C
!C *********************************************************************
!C  DETERMINE HOW MANY BLOCKS ARE NEEDED IN EACH REORDER GROUP (SUNRISE,
!C                 SUNSET, STIFFNESS) AFTER REORDERING
!C *********************************************************************
!C NBLOCKROW = # BLOCKS OF EACH REORDER GROUP (STIFFNESS, SUNRISE, SUNSET)
!C IUSESIZE  = # OF GRID CELLS IN EACH GRID BLOCK
!C NCELLROW  = # OF GRID CELLS IN EACH REORDER GROUP
!C NNORISE   = # OF STIFFNESS (NON-SUNRISE, NON-SUNSET) CELLS
!C NSUNRISE  = # OF SUNRISE CELLS
!C NSUNSET   = # OF SUNSET  CELLS
!C NREBLOCK  = COUNTS NUMBER OF NEW BLOCKS
!C
!
!          NSUMBLOCK        = 0
!          NREBLOCK         = 0
!          JLOOPLO          = 0 
!C
!          !write(6,*) 'norise,sunrise,sunset=',nnorise,nsunrise,nsunset
!          DO 770 IT        = 1, 3
!           IF (IT.EQ.1) THEN
!            NCELLROW       = NNORISE
!            NBLOCKROW      = 1 + NCELLROW / (KULOOP    + 0.0001d0)
!           ELSEIF (IT.EQ.2) THEN
!            NCELLROW       = NSUNRISE
!            !NBLOCKROW      = 1 + NCELLROW * 3./ (KULOOP  + 0.0001d0)
!           NBLOCKROW      = 1 + NCELLROW / (KULOOP  + 0.0001d0)
!           ELSEIF (IT.EQ.3) THEN
!            NCELLROW       = NSUNSET
!            !NBLOCKROW      = 1 + NCELLROW * 3./ (KULOOP  + 0.0001d0)
!           NBLOCKROW      = 1 + NCELLROW / (KULOOP  + 0.0001d0)
!           ENDIF
!C
!           NSUMBLOCK       = NSUMBLOCK + NBLOCKROW
!C
!           IF (NSUMBLOCK.GT.MXBLOCK) THEN
!!            write(6,*) 'val of mxblock= ',mxblock
!!            WRITE(6,*)'PHYSPROC: NSUMBLOCK>MXBLOCK. INCREASE MXBLOCK ',
!!     1                 NSUMBLOCK, NNORISE, NSUNRISE, NSUNSET, KULOOP 
!            STOP
!           ENDIF
!C
!           IF (NCELLROW.EQ.0) THEN
!            NBLOCKROW      = 0
!           ELSE
!            IAVBLOK            = 1 + NCELLROW / (NBLOCKROW + 0.0001d0)
!            IAVGSIZE           = MIN(IAVBLOK,KULOOP)
!            IREMAIN            = NCELLROW
!C
!            !write(6,*) 'it,nblockrow,iavesize= ',it,nblockrow,iavgsize
!            DO 765 KBLK        = 1, NBLOCKROW
!             NREBLOCK          = NREBLOCK + 1
!             IUSESIZE          = MIN(IAVGSIZE,MAX(IREMAIN,0))
!             JLOWVAR(NREBLOCK) = JLOOPLO
!             KTLPVAR(NREBLOCK) = IUSESIZE
!             IREMAIN           = IREMAIN - IUSESIZE
!             JLOOPLO           = JLOOPLO + IUSESIZE
! 765        CONTINUE
!           ENDIF
! 770      CONTINUE
!C770      CONTINUE IT = 1, 3
!C
!        ENDIF
!C       ENDIF ISREORD.EQ.1 
!C
!------------------------------------------------------------------------------
 855   CONTINUE
 860  CONTINUE
!
! *********************************************************************
! ********************* END OF SUBROUTINE PHYSPROC.F ******************
! *********************************************************************
!
      RETURN
      END SUBROUTINE PHYSPROC
