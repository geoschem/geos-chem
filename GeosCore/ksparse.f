! $Id: ksparse.f,v 1.1 2009/09/16 14:06:22 bmy Exp $
      SUBROUTINE KSPARSE
!
!******************************************************************************
!  Subroutine KSPARSE sets up the sparse-matrix arrays, and also arrays for
!  day & night chemistry for SMVGEAR II. (M. Jacobson 1997; bdf, bmy, 4/18/03)
!
!  NOTES:
!  (1 ) Now direct some output to "smv2.log" file.  Now call GEOS_CHEM_STOP
!        to deallocate all arrays and stop the run safely.  Now also force
!        double-precision with "D" exponents. (bmy, 4/18/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

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
C   K    K  SSSSSSS  PPPPPPP     A      RRRRRRR  SSSSSSS  EEEEEEE
C   K  K    S        P     P    A A     R     R  S        E 
C   KK      SSSSSSS  PPPPPPP   A   A    RRRRRRR  SSSSSSS  EEEEEEE
C   K  K          S  P        AAAAAAA   R  R           S  E
C   K    K  SSSSSSS  P       A       A  R    R   SSSSSSS  EEEEEEE  
C
C *********************************************************************
C * THIS ROUTINE SETS UP SPARSE-MATRIX AND OTHER ARRAYS. IT ALSO      *
C * SETS ARRAYS FOR GAS-PHASE, AQUEOUS-PHASE, OR ANY OTHER TYPE       *
C * OF CHEMISTRY. FINALLY, IT SETS ARRAYS FOR BOTH DAY AND NIGHT      *
C * CHEMISTRY OF EACH TYPE.                                           *
C *                                                                   *
C * HOW TO CALL SUBROUTINE:                                           *
C * ----------------------                                            *
C *  CALL KSPARSE.F FROM JSPARSE.F WITH                               * 
C *     NCS  = 1..NCSGAS FOR GAS CHEMISTRY                            *
C *********************************************************************
C
C *********************************************************************
C * SETS UP ARRAYS FOR DECOMPOSITION / BACK-SUBSTITUTION OF SPARSE    *
C * MATRICES BY REMOVING ALL CALCULATIONS INVOLVING A ZERO.           *
C *********************************************************************
C
C *********************************************************************
C *********************************************************************
C **        SET ARRAYS TO TAKE ADVANTAGE OF SPARSE MATRICES          ** 
C *********************************************************************
C *********************************************************************
C
C IFSUN    = 1 THEN DAY-CHEMISTRY;  = 2 THEN NIGHT CHEMISTRY
C NCSP     = NCS        FOR DAYTIME   GAS CHEM  
C NCSP     = NCS   +ICS FOR NIGHTTIME GAS CHEM  
C
C KOUNT0A  = # INITIAL MATRIX SPOTS FILLED W/O  SPARSE-MATRIX REDUCTIONS 
C KOUNT0   = # INITIAL MATRIX SPOTS FILLED WITH SPARSE-MATRIX REDUCTIONS 
C KNTARRAY = # FINAL MATRIX SPOTS FILLED W/O  SPARSE-MATRIX REDUCTIONS 
C IARRAY2  = # FINAL MATRIX SPOTS FILLED WITH SPARSE-MATRIX REDUCTIONS 
C ICNTA    = # OPERATIONS IN DECOMP LOOP 1 W/O  SPARSE-MATRIX REDUCTIONS 
C ICNTB    = # OPERATIONS IN DECOMP LOOP 1 WITH SPARSE-MATRIX REDUCTIONS 
C JCNTA    = # OPERATIONS IN DECOMP LOOP 2 W/O  SPARSE-MATRIX REDUCTIONS 
C JCNTB    = # OPERATIONS IN DECOMP LOOP 2 WITH SPARSE-MATRIX REDUCTIONS 
C KCNTA    = # OPERATIONS IN BACK-SUP LOOP 1 W/O  SPARSE-MATRIX REDUCTIONS 
C KCNTB    = # OPERATIONS IN BACK-SUB LOOP 1 WITH SPARSE-MATRIX REDUCTIONS 
C MCNTA    = # OPERATIONS IN BACK-SUP LOOP 2 W/O  SPARSE-MATRIX REDUCTIONS 
C MCNTB    = # OPERATIONS IN BACK-SUB LOOP 2 WITH SPARSE-MATRIX REDUCTIONS 
C
C LZERO    = 1 IF AN ARRAY SPOT IS FILLED WITH A NON-ZERO VALUE. LZERO
C            IS UPDATED AS WE SIMULATE THE ORDER OF CALCULATIONS DURING
C            A PRACTICE L-U DECOMPOSITION
C

      INTEGER KOUNT0A,KOUNT0,ICNTA,ICNTB
      INTEGER KCNTA,KCNTB,MCNTA,MCNTV,IARRAY2,J,K,J1,I,I1,I2,KNTARRAY
      INTEGER IZIL,NREMAIN,NFIVE,NFOUR,NTHREE,NTWO,NONE,IC,KA,KB,KC,KD
      ! Bug fix (gcc)
      !INTEGER IA,KZIL,MC,JCNTA,JCNTB,MCNTA,MCNTB,KE,MZIL
      INTEGER IA,KZIL,MC,JCNTA,JCNTB,MCNTB,KE,MZIL

      INTEGER, SAVE :: MCNT,KCNT,ICNT,JCNT,MZTOT,IJTOT,KZTOT,IDECOMP
      INTEGER, SAVE :: MCCOUNT,ICCOUNT,JCCOUNT,KCCOUNT,KBSUB,MBSUB

       IF (IFNEVER.EQ.0) THEN
        IFNEVER               = 1
        ICNT                  = 0 
        JCNT                  = 0 
        KCNT                  = 0 
        MCNT                  = 0 
        ICCOUNT               = 0
        JCCOUNT               = 0
        KCCOUNT               = 0 
        MCCOUNT               = 0 
        IDECOMP               = 0
        KBSUB                 = 0
        MBSUB                 = 0
        IJTOT                 = 0 
        KZTOT                 = 0 
        MZTOT                 = 0 
       ENDIF
C
       KOUNT0A                = 0
       KOUNT0                 = 0
       ICNTA                  = 0
       ICNTB                  = 0
       JCNTA                  = 0
       JCNTB                  = 0
       KCNTA                  = 0
       KCNTB                  = 0
       MCNTA                  = 0
       MCNTB                  = 0
       IARRAY2                = 0 
C
       DO 522 J               = 1, ISCHANG(NCS)
        DO 520 K              = 1, ISCHANG(NCS)
         KOUNT0A              = KOUNT0A + 1
         IF (LZERO(K,J).EQ.1) KOUNT0 = KOUNT0 + 1
         JARRAYPT(K,J)        = 0
 520    CONTINUE
 522   CONTINUE
C
C *********************************************************************
C **                ARRAYS FOR DECOMPOSITION (LUDCMP)                ** 
C *********************************************************************
C IZILCH = # OF CALCULATIONS WITH NON-ZERO VALUES DURING MATRIX DECOMP
C IZERO  = EACH OCCURRENCE OF EACH IZILCH CALCULATION
C
       DO 562 J                  = 1, ISCHANG(NCS)
        JZILCH(J)                = 0 
        J1                       = J - 1
C
C ------------------- FIRST LOOP OF DECOMPOSTION ----------------------
C
        DO 542 I            = 2, ISCHANG(NCS) 
         IZILCH(J,I)        = 0 
         I1                 = J1 
         IF (I.LE.J1) I1    = I - 1
         DO 540 K           = 1, I1
          ICNTA             = ICNTA + 1
          IF (LZERO(I,K).EQ.1.AND.LZERO(K,J).EQ.1) THEN
           IZILCH(J,I)      = IZILCH(J,I) + 1
           ICNT             = ICNT  + 1
           ICNTB            = ICNTB + 1
           IZEROK(ICNT)     = K   
           LZERO(I,J)       = 1 
          ENDIF
 540     CONTINUE
 542    CONTINUE
C
C ------------------- SECOND LOOP OF DECOMPOSTION ---------------------
C
C JZILCH  = # OF CALCULATIONS WITH NON-ZERO VALUES TO FILL LOWER
C           PART OF DECOMPOSED MATRIX   
C
        DO 560 I            = J+1, ISCHANG(NCS) 
         JCNTA              = JCNTA + 1
         IF (LZERO(I,J).EQ.1) THEN
          JZILCH(J)         = JZILCH(J) + 1
          JCNT              = JCNT  + 1
          JCNTB             = JCNTB + 1
          JZERO(JCNT)       = I  
         ENDIF 
 560    CONTINUE 
 562   CONTINUE 
C
C *********************************************************************
C **              ARRAYS FOR BACK-SUBSTITUTION (LUBKSB)              ** 
C *********************************************************************
C JZILCH AND KZILCH HAVE SAME NUMBER OF TOTAL ELEMENTS
C BOTH CONTAIN NON-ZEROS IN LOWER TRIANGLULAR MATRIX 
C
C
C ------------------ FIRST LOOP OF BACK-SUBSTITUTION ------------------
C
       DO 572 I             = 2, ISCHANG(NCS)
        KZILCH(I)           = 0 
        I1                  = I - 1
        DO 570 J            = 1, I1    
         KCNTA              = KCNTA + 1
         IF (LZERO(I,J).EQ.1) THEN 
          KZILCH(I)         = KZILCH(I) + 1
          KCNTB             = KCNTB     + 1
          KCNT              = KCNT      + 1 
          IARRAY2           = IARRAY2   + 1
          KZERO(KCNT)       = J
          JARRAYPT(I,J)     = IARRAY2 
         ENDIF
 570    CONTINUE 
 572   CONTINUE 
C
C ----------------- SECOND LOOP OF BACK-SUBSTITUTION ------------------
C
C MZILCH CONTAINS NON-ZEROS FOR UPPER TRIANGULAR MATRIX, WHERE BACK-
C SUBSTITUTION OCCURS. 
C
       DO 577 I             = ISCHANG(NCS), 1, -1
        MZILCH(I)           = 0 
        I2                  = I + 1
        DO 575 J            = I+1, ISCHANG(NCS)
         MCNTA              = MCNTA + 1
         IF (LZERO(I,J).EQ.1) THEN 
          MZILCH(I)         = MZILCH(I) + 1
          MCNTB             = MCNTB     + 1
          MCNT              = MCNT      + 1
          IARRAY2           = IARRAY2   + 1
          MZERO(MCNT)       = J
          JARRAYPT(I,J)     = IARRAY2 
         ENDIF
 575    CONTINUE
 577   CONTINUE
C
C *********************************************************************
C * FILL JARRAYPT WITH REMAINING ARRAY POINTS (ALONG DIAGONAL)        *
C *********************************************************************
C
       DO 580 I             = 1, ISCHANG(NCS) 
        IARRAY2             = IARRAY2 + 1
        JARRAYPT(I,I)       = IARRAY2 
 580   CONTINUE
C 
       IARRAY(NCSP)         = IARRAY2 
       KNTARRAY             = KCNTA + MCNTA + ISCHANG(NCS)
C
C *********************************************************************
C *** CHANGE IZERO AND JZERO ARRAYS SO THEIR VALUES POINT TO NEW    ***
C ***              ARRAY POSITIONS DEFINED IN JARRAYPT              *** 
C *********************************************************************
C
C JARRAYPT = IDENTIFIES THE ONE-DIMENSIONAL ARRAY POINT FOR EACH TWO-
C            DIMENSIONAL POINT I,J
C IARRAY   = THE LENGTH OF THE ONE-DIMENSIONAL ARRAY HOLDING ALL
C            SPARSE MATRIX POINTS = SPARSE-MATRIX DIMENSION
C IZER2    = USED TO IDENTIFY THE 1-D ARRAY POINT FOR EACH K,J VALUE
C            FOUND IN THE FIRST MAJOR LOOP OF MATRIX DECOMPOSITION
C IZERO    = USED TO FIND THE 1-D ARRAY POINT FOR EACH I,K VALUE
C            FOUND IN THE SAME LOOP.
C
       DO 595 J                  = 1, ISCHANG(NCS)
C
C ------------------- FIRST LOOP OF DECOMPOSTION ----------------------
C
        IJTLO(J,NCSP)             = IJTOT + 1
        DO 605 I                  = 2, ISCHANG(NCS)
         IZIL                     = IZILCH(J,I) 
         IF (IZIL.GT.0) THEN 
          IJTOT                   = IJTOT + 1
          NREMAIN                 = IZIL
          NFIVE                   = (NREMAIN + 0.0001d0) / 5 
          NREMAIN                 =  NREMAIN - NFIVE   * 5 
          NFOUR                   = (NREMAIN + 0.0001d0) / 4 
          NREMAIN                 =  NREMAIN - NFOUR   * 4
          NTHREE                  = (NREMAIN + 0.0001d0) / 3  
          NREMAIN                 =  NREMAIN - NTHREE  * 3 
          NTWO                    = (NREMAIN + 0.0001d0) / 2   
          NREMAIN                 =  NREMAIN - NTWO    * 2  
          NONE                    = (NREMAIN + 0.0001d0)  
          NREMAIN                 =  NREMAIN - NONE
C

          IJVAL(IJTOT)            = JARRAYPT(I,J) 
          IDL5( IJTOT)            = IDECOMP     + 1
          IDH5( IJTOT)            = IDECOMP     + NFIVE  
          IDL4( IJTOT)            = IDH5(IJTOT) + 1
          IDH4( IJTOT)            = IDH5(IJTOT) + NFOUR   
          IDL3( IJTOT)            = IDH4(IJTOT) + 1
          IDH3( IJTOT)            = IDH4(IJTOT) + NTHREE
          IDL2( IJTOT)            = IDH3(IJTOT) + 1
          IDH2( IJTOT)            = IDH3(IJTOT) + NTWO
          IDL1( IJTOT)            = IDH2(IJTOT) + 1
          IDH1( IJTOT)            = IDH2(IJTOT) + NONE
          IDECOMP                 = IDH1(IJTOT)
C
          DO 601 IC          = IDL5(IJTOT), IDH5(IJTOT)
           KA                = IZEROK(ICCOUNT+1) 
           KB                = IZEROK(ICCOUNT+2) 
           KC                = IZEROK(ICCOUNT+3) 
           KD                = IZEROK(ICCOUNT+4) 
           KE                = IZEROK(ICCOUNT+5) 
           ICCOUNT           = ICCOUNT + 5
           IKDECA(IC)        = JARRAYPT(I,KA)
           IKDECB(IC)        = JARRAYPT(I,KB)
           IKDECC(IC)        = JARRAYPT(I,KC)
           IKDECD(IC)        = JARRAYPT(I,KD)
           IKDECE(IC)        = JARRAYPT(I,KE)
           KJDECA(IC)        = JARRAYPT(KA,J)
           KJDECB(IC)        = JARRAYPT(KB,J)
           KJDECC(IC)        = JARRAYPT(KC,J)
           KJDECD(IC)        = JARRAYPT(KD,J)
           KJDECE(IC)        = JARRAYPT(KE,J)
 601      CONTINUE 
C
          DO 602 IC          = IDH5(IJTOT) + 1, IDH4(IJTOT)
           KA                = IZEROK(ICCOUNT+1) 
           KB                = IZEROK(ICCOUNT+2) 
           KC                = IZEROK(ICCOUNT+3) 
           KD                = IZEROK(ICCOUNT+4) 
           ICCOUNT           = ICCOUNT + 4 
           IKDECA(IC)        = JARRAYPT(I,KA)
           IKDECB(IC)        = JARRAYPT(I,KB)
           IKDECC(IC)        = JARRAYPT(I,KC)
           IKDECD(IC)        = JARRAYPT(I,KD)
           KJDECA(IC)        = JARRAYPT(KA,J)
           KJDECB(IC)        = JARRAYPT(KB,J)
           KJDECC(IC)        = JARRAYPT(KC,J)
           KJDECD(IC)        = JARRAYPT(KD,J)
 602      CONTINUE 
C
          DO 603 IC          = IDH4(IJTOT) + 1, IDH3(IJTOT)
           KA                = IZEROK(ICCOUNT+1) 
           KB                = IZEROK(ICCOUNT+2) 
           KC                = IZEROK(ICCOUNT+3) 
           ICCOUNT           = ICCOUNT + 3 
           IKDECA(IC)        = JARRAYPT(I,KA)
           IKDECB(IC)        = JARRAYPT(I,KB)
           IKDECC(IC)        = JARRAYPT(I,KC)
           KJDECA(IC)        = JARRAYPT(KA,J)
           KJDECB(IC)        = JARRAYPT(KB,J)
           KJDECC(IC)        = JARRAYPT(KC,J)
 603      CONTINUE 
C
          DO 604 IC          = IDH3(IJTOT) + 1, IDH2(IJTOT)
           KA                = IZEROK(ICCOUNT+1) 
           KB                = IZEROK(ICCOUNT+2) 
           ICCOUNT           = ICCOUNT + 2 
           IKDECA(IC)        = JARRAYPT(I,KA)
           IKDECB(IC)        = JARRAYPT(I,KB)
           KJDECA(IC)        = JARRAYPT(KA,J)
           KJDECB(IC)        = JARRAYPT(KB,J)
 604      CONTINUE 
C
          DO 606 IC          = IDH2(IJTOT) + 1, IDH1(IJTOT)
           KA                = IZEROK(ICCOUNT+1) 
           ICCOUNT           = ICCOUNT + 1
           IKDECA(IC)        = JARRAYPT(I,KA)
           KJDECA(IC)        = JARRAYPT(KA,J)
 606      CONTINUE

         ENDIF
 605    CONTINUE
C

        IJTHI(J,NCSP)           = IJTOT
C
C ------------------ DIAGONAL TERM OF DECOMPOSTION --------------------
C
        JARRDIAG(J,NCSP)        = JARRAYPT(J,J)
C
C ------------------- SECOND LOOP OF DECOMPOSTION ---------------------
C
        JLOZ1(J,NCSP)           = JCCOUNT + 1
        DO 635 I                = 1, JZILCH(J)
         JCCOUNT                = JCCOUNT + 1
         IA                     = JZERO(JCCOUNT)
         JZEROA(JCCOUNT)        = JARRAYPT(IA,J)
 635    CONTINUE
        JHIZ1(J,NCSP)           = JCCOUNT 
C
 595   CONTINUE

C
C *********************************************************************
C ** CREATE MORE BACK-SUBSTITUTION ARRAYS TO INCREASE EFFICIENCY     ** 
C *********************************************************************
C
C ------------------ FIRST LOOP OF BACK-SUBSTITUTION ------------------
C
       KZTLO(NCSP)          = KZTOT + 1 
       DO 620 I             = 2, ISCHANG(NCS)
        KZIL                = KZILCH(I) 
        IF (KZIL.GT.0) THEN 
         KZTOT              = KZTOT + 1
         NREMAIN            = KZIL
         NFIVE              = (NREMAIN + 0.0001d0) / 5 
         NREMAIN            =  NREMAIN - NFIVE   * 5 
         NFOUR              = (NREMAIN + 0.0001d0) / 4 
         NREMAIN            =  NREMAIN - NFOUR   * 4
         NTHREE             = (NREMAIN + 0.0001d0) / 3  
         NREMAIN            =  NREMAIN - NTHREE  * 3 
         NTWO               = (NREMAIN + 0.0001d0) / 2   
         NREMAIN            =  NREMAIN - NTWO    * 2  
         NONE               = (NREMAIN + 0.0001d0)  
         NREMAIN            =  NREMAIN - NONE
C

         IKZTOT(KZTOT)      = I 
         KBL5( KZTOT)       = KBSUB       + 1
         KBH5( KZTOT)       = KBSUB       + NFIVE  
         KBL4( KZTOT)       = KBH5(KZTOT) + 1
         KBH4( KZTOT)       = KBH5(KZTOT) + NFOUR   
         KBL3( KZTOT)       = KBH4(KZTOT) + 1
         KBH3( KZTOT)       = KBH4(KZTOT) + NTHREE
         KBL2( KZTOT)       = KBH3(KZTOT) + 1
         KBH2( KZTOT)       = KBH3(KZTOT) + NTWO
         KBL1( KZTOT)       = KBH2(KZTOT) + 1
         KBH1( KZTOT)       = KBH2(KZTOT) + NONE
         KBSUB              = KBH1(KZTOT)
C
         DO 611 KC          = KBL5(KZTOT), KBH5(KZTOT)
          KZEROA(KC)        = KZERO(KCCOUNT+1) 
          KZEROB(KC)        = KZERO(KCCOUNT+2) 
          KZEROC(KC)        = KZERO(KCCOUNT+3) 
          KZEROD(KC)        = KZERO(KCCOUNT+4) 
          KZEROE(KC)        = KZERO(KCCOUNT+5) 
          KCCOUNT           = KCCOUNT + 5
 611     CONTINUE 
C
         DO 612 KC          = KBL4(KZTOT), KBH4(KZTOT)
          KZEROA(KC)        = KZERO(KCCOUNT+1) 
          KZEROB(KC)        = KZERO(KCCOUNT+2) 
          KZEROC(KC)        = KZERO(KCCOUNT+3) 
          KZEROD(KC)        = KZERO(KCCOUNT+4) 
          KCCOUNT           = KCCOUNT + 4 
 612     CONTINUE 
C
         DO 613 KC          = KBL3(KZTOT), KBH3(KZTOT)
          KZEROA(KC)        = KZERO(KCCOUNT+1) 
          KZEROB(KC)        = KZERO(KCCOUNT+2) 
          KZEROC(KC)        = KZERO(KCCOUNT+3) 
          KCCOUNT           = KCCOUNT + 3  
 613     CONTINUE 
C

         DO 614 KC          = KBL2(KZTOT), KBH2(KZTOT)
          KZEROA(KC)        = KZERO(KCCOUNT+1) 
          KZEROB(KC)        = KZERO(KCCOUNT+2) 
          KCCOUNT           = KCCOUNT + 2  
 614     CONTINUE 
C
         DO 615 KC          = KBL1(KZTOT), KBH1(KZTOT)
          KZEROA(KC)        = KZERO(KCCOUNT+1) 
          KCCOUNT           = KCCOUNT + 1  
 615     CONTINUE 
        ENDIF 
 620   CONTINUE 
       KZTHI(NCSP)          = KZTOT 
C
C ----------------- SECOND LOOP OF BACK-SUBSTITUTION ------------------
C

       DO 640 I             = ISCHANG(NCS), 1, -1
        MZIL                = MZILCH(I) 
        IF (MZIL.GT.0) THEN 
         MZTOT              = MZTOT + 1
         NREMAIN            = MZIL
         NFIVE              = (NREMAIN + 0.0001d0) / 5 
         NREMAIN            =  NREMAIN - NFIVE   * 5 
         NFOUR              = (NREMAIN + 0.0001d0) / 4 
         NREMAIN            =  NREMAIN - NFOUR   * 4
         NTHREE             = (NREMAIN + 0.0001d0) / 3  
         NREMAIN            =  NREMAIN - NTHREE  * 3 
         NTWO               = (NREMAIN + 0.0001d0) / 2   
         NREMAIN            =  NREMAIN - NTWO    * 2  
         NONE               = (NREMAIN + 0.0001d0)  
         NREMAIN            =  NREMAIN - NONE
C
         IMZTOT(I,NCSP)     = MZTOT 
         MBL5(  MZTOT)      = MBSUB       + 1
         MBH5(  MZTOT)      = MBSUB       + NFIVE  
         MBL4(  MZTOT)      = MBH5(MZTOT) + 1
         MBH4(  MZTOT)      = MBH5(MZTOT) + NFOUR   
         MBL3(  MZTOT)      = MBH4(MZTOT) + 1
         MBH3(  MZTOT)      = MBH4(MZTOT) + NTHREE
         MBL2(  MZTOT)      = MBH3(MZTOT) + 1
         MBH2(  MZTOT)      = MBH3(MZTOT) + NTWO
         MBL1(  MZTOT)      = MBH2(MZTOT) + 1
         MBH1(  MZTOT)      = MBH2(MZTOT) + NONE
         MBSUB              = MBH1(MZTOT)
C
         DO 631 MC          = MBL5(MZTOT), MBH5(MZTOT)
          MZEROA(MC)        = MZERO(MCCOUNT+1) 
          MZEROB(MC)        = MZERO(MCCOUNT+2) 
          MZEROC(MC)        = MZERO(MCCOUNT+3) 
          MZEROD(MC)        = MZERO(MCCOUNT+4) 
          MZEROE(MC)        = MZERO(MCCOUNT+5) 
          MCCOUNT           = MCCOUNT + 5
 631     CONTINUE 
C
         DO 632 MC          = MBL4(MZTOT), MBH4(MZTOT)
          MZEROA(MC)        = MZERO(MCCOUNT+1) 
          MZEROB(MC)        = MZERO(MCCOUNT+2) 
          MZEROC(MC)        = MZERO(MCCOUNT+3) 
          MZEROD(MC)        = MZERO(MCCOUNT+4) 
          MCCOUNT           = MCCOUNT + 4 
 632     CONTINUE 
C
         DO 633 MC          = MBL3(MZTOT), MBH3(MZTOT)
          MZEROA(MC)        = MZERO(MCCOUNT+1) 
          MZEROB(MC)        = MZERO(MCCOUNT+2) 
          MZEROC(MC)        = MZERO(MCCOUNT+3) 
          MCCOUNT           = MCCOUNT + 3  
 633     CONTINUE 
C
         DO 634 MC          = MBL2(MZTOT), MBH2(MZTOT)
          MZEROA(MC)        = MZERO(MCCOUNT+1) 
          MZEROB(MC)        = MZERO(MCCOUNT+2) 
          MCCOUNT           = MCCOUNT + 2  
 634     CONTINUE 
C
         DO 636 MC          = MBL1(MZTOT), MBH1(MZTOT)
          MZEROA(MC)        = MZERO(MCCOUNT+1) 
          MCCOUNT           = MCCOUNT + 1  
 636     CONTINUE 
        ENDIF 
 640   CONTINUE 
C
C *********************************************************************
C **        CHECK DIMENSIONS AND PRINT OUT ARRAY SAVINGS             ** 
C *********************************************************************
C
       IF (ICNT    .GT. MXCOUNT2  .OR. JCNT     .GT. MXCOUNT3   .OR. 
     1     KCNT    .GT. MXCOUNT3  .OR. MCNT     .GT. MXCOUNT3   .OR. 
     2     ICCOUNT .GT. MXCOUNT2  .OR. JCCOUNT  .GT. MXCOUNT3   .OR. 
     3     KCCOUNT .GT. MXCOUNT3  .OR. MCCOUNT  .GT. MXCOUNT3   .OR. 
     4     IJTOT   .GT. MXCOUNT3  .OR. IDECOMP  .GT. MXCOUNT3   .OR.  
     5     KZTOT   .GT. MXCOUNT4  .OR. KBSUB    .GT. MXCOUNT4   .OR.  
     6     MZTOT   .GT. MXCOUNT4  .OR. MBSUB    .GT. MXCOUNT4   .OR.  
     7     IARRAY2 .GT. MXARRAY) THEN
C
        WRITE(6,705)
     1  MXCOUNT2, ICNT,        MXCOUNT3, JCNT,
     2  MXCOUNT3, KCNT,        MXCOUNT3, MCNT,
     3  MXCOUNT2, ICCOUNT,     MXCOUNT3, JCCOUNT,
     4  MXCOUNT3, KCCOUNT,     MXCOUNT3, MCCOUNT,
     5  MXCOUNT3, IJTOT,       MXCOUNT3, IDECOMP,
     6  MXCOUNT4, KZTOT,       MXCOUNT4, KBSUB,
     7  MXCOUNT4, MZTOT,       MXCOUNT4, MBSUB,
     8  MXARRAY,  IARRAY2     
        CALL GEOS_CHEM_STOP
       ENDIF
C
 705   FORMAT('KSPARSE: ONE OF THE DIMENSIONS BELOW IS TOO SMALL:',/,
     1        'DIMENSION: MXCOUNT2 = ',I5,' VARIABLE: ICNT     = ',I5,/,  
     2        'DIMENSION: MXCOUNT3 = ',I5,' VARIABLE: JCNT     = ',I5,/,  
     3        'DIMENSION: MXCOUNT3 = ',I5,' VARIABLE: KCNT     = ',I5,/,  
     4        'DIMENSION: MXCOUNT3 = ',I5,' VARIABLE: MCNT     = ',I5,/,  
     5        'DIMENSION: MXCOUNT2 = ',I5,' VARIABLE: ICCOUNT  = ',I5,/,  
     6        'DIMENSION: MXCOUNT3 = ',I5,' VARIABLE: JCCOUNT  = ',I5,/,  
     7        'DIMENSION: MXCOUNT3 = ',I5,' VARIABLE: KCCOUNT  = ',I5,/,  
     8        'DIMENSION: MXCOUNT3 = ',I5,' VARIABLE: MCCOUNT  = ',I5,/,  
     9        'DIMENSION: MXCOUNT3 = ',I5,' VARIABLE: IJTOT    = ',I5,/, 
     1        'DIMENSION: MXCOUNT3 = ',I5,' VARIABLE: IDECOMP  = ',I5,/, 
     2        'DIMENSION: MXCOUNT4 = ',I5,' VARIABLE: KZTOT    = ',I5,/, 
     3        'DIMENSION: MXCOUNT4 = ',I5,' VARIABLE: KBSUB    = ',I5,/, 
     4        'DIMENSION: MXCOUNT4 = ',I5,' VARIABLE: MZTOT    = ',I5,/, 
     5        'DIMENSION: MXCOUNT4 = ',I5,' VARIABLE: MBSUB    = ',I5,/, 
     6        'DIMENSION: MXARRAY  = ',I5,' VARIABLE: IARRAY2  = ',I5)
C
        WRITE(IO93,655)NCSP,KOUNT0A,KOUNT0,KNTARRAY,IARRAY2,ICNTA,ICNTB,
     1               JCNTA,JCNTB,KCNTA,KCNTB,MCNTA,MCNTB 
C
 655    FORMAT(/'PARAM    POSS MATRIX POINTS -- NONZEROS -- NCSP=',I4/
     1          'INITMAT  ',4X,I8,9X,I8/                               
     2          'FINMAT   ',4X,I8,9X,I8/   
     3          'DECOMP1  ',4X,I8,9X,I8/ 
     4          'DECOMP2  ',4X,I8,9X,I8/ 
     5          'BACKSB1  ',4X,I8,9X,I8/ 
     6          'BACKSB2  ',4X,I8,9X,I8/)  
C
C *********************************************************************
C *           SET COEFFICIENTS OF THE INTEGRATION METHOD              * 
C *********************************************************************
C
C PARAMETERS USED IN SMVGEAR
C --------------------------
C PERTST   = COEFFICIENTS USED TO SELECT THE STEP-SIZE AND ORDER. THUS,
C            ONLY ABOUT ONE-PERCENT ACCURACY NEEDED. SEE GEAR(1971)
C            OR HINDMARSH '73 UCID-30059.  
C ASET     = PARAMETERS FOR DETERMINING THE ORDER OF THE INTEGRATION METHOD
C            AND FOR CALCULATION THE MATRIX P.
C MSTEP    = MAXIMUM NUMBER OF CORRECTOR ITERATIONS ALLOWED
C HMIN     = MINIMUM TIME-STEP ALLOWED (SEC)
C MAXORD   = MAXIMUM ORDER OF THE METHOD USED
C MBETWEEN = MAXIMUM NUMBER OF STEPS BETWEEN CALLS TO PDERIV
C NQQ      = ORDER OF THE INTEGRATION METHOD
C
      IF (IFDID.EQ.0) THEN 
       IFDID = 1 
C
       ! Now force double-precision with "D" exponents (bmy, 4/18/03)
       DATA PERTST / 
     1  2.0d0, 4.5d0, 7.333d0, 10.42d0,   13.7d0,    17.15d0,     1.0d0,
     3  3.0d0, 6.0d0, 9.167d0, 12.5d0,    15.98d0,    1.0d0,      1.0d0,
     5  1.0d0, 1.0d0, 0.5d0,    0.1667d0,  0.04133d0, 0.008267d0, 1.0d0/ 

C
C ADAMS-MOULTON COEFFICIENTS
C
C    2       2.0, 12.0, 24.0,   37.89,   53.33,    70.08,    87.97,
C    4      12.0, 24.0, 37.89,  53.33,   70.08,    87.97,     1.0,
C    6       1.0,  1.0,  2.0,    1.0,     0.3157,   0.07407,  0.0139 / 
C
       MSTEP          = 3
       HMIN           = 1.0d-15 
       MAXORD         = 5
       MBETWEEN       = 50
C
       DO 800 NQQ     = 1, 7
        ENQQ1(NQQ)    = 0.5d0 / FLOAT(NQQ    )
        ENQQ2(NQQ)    = 0.5d0 / FLOAT(NQQ + 1) 
        ENQQ3(NQQ)    = 0.5d0 / FLOAT(NQQ + 2)
        CONPST(NQQ)   = 1.0d0 / (PERTST(NQQ,1) * ENQQ3(NQQ)) 
        CONP15(NQQ)   = 1.5d0 * CONPST(NQQ)
        PERTS2(NQQ,1) = PERTST(NQQ,1) * PERTST(NQQ,1)
        PERTS2(NQQ,2) = PERTST(NQQ,2) * PERTST(NQQ,2)
        PERTS2(NQQ,3) = PERTST(NQQ,3) * PERTST(NQQ,3)
 800   CONTINUE
C
       DO 830 I2   = 1, 6 
        ASET(I2,2) = 1.0d0
        ASET(I2,8) = 0.d0
 830   CONTINUE
C
       ASET(1,1)   = 1.0d0                                                    
C
       ASET(2,1)   = 2.0d0    /    3.0d0 
       ASET(2,3)   = 1.0d0    /    3.0d0  
C
       ASET(3,1)   = 6.0d0    /   11.0d0
       ASET(3,3)   = 6.0d0    /   11.0d0 
       ASET(3,4)   = 1.0d0    /   11.0d0                                     
C
       ASET(4,1)   = 12.0d0   /   25.0d0                                      
       ASET(4,3)   =   .70d0    
       ASET(4,4)   =   .20d0   
       ASET(4,5)   =   .020d0  
C
       ASET(5,1)   =   60.0d0 /  137.0d0
       ASET(5,3)   =  225.0d0 /  274.0d0  
       ASET(5,4)   =   85.0d0 /  274.0d0  
       ASET(5,5)   =   15.0d0 /  274.0d0  
       ASET(5,6)   =    1.0d0 /  274.0d0                                      
C
       ASET(6,1)   =  180.0d0 /  441.0d0 
       ASET(6,3)   =  406.0d0 /  441.0d0   
       ASET(6,4)   =  735.0d0 / 1764.0d0  
       ASET(6,5)   =  175.0d0 / 1764.0d0 
       ASET(6,6)   =   21.0d0 / 1764.0d0  
       ASET(6,7)   =    1.0d0 / 1764.0d0  
C
      ENDIF
C     ENDIF IFDID.EQ.0
C
C *********************************************************************
C ********************** END OF SUBROUTINE KSPARSE ********************
C *********************************************************************
C
      RETURN                                                             
      END SUBROUTINE KSPARSE
