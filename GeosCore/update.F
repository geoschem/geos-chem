! $Id: update.f,v 1.1 2009/09/16 14:05:58 bmy Exp $
      SUBROUTINE UPDATE
!
!******************************************************************************
!  Subroutine UPDATE updates rxn rates for each timestep for SMVGEAR II. 
!  (M. Jacobson, 1997, bdf, bmy, 4/18/03)
!
!  NOTES:
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
! ***                        (650) 723-6836                         *** 
! *********************************************************************
!
!        U     U  PPPPPPP  DDDDDD      A    TTTTTTT  EEEEEEE
!        U     U  P     P  D     D    A A      T     E  
!        U     U  PPPPPPP  D     D   A   A     T     EEEEEEE 
!        U     U  P        D     D  AAAAAAA    T     E  
!        UUUUUUU  P        DDDDDD  A       A   T     EEEEEEE
!
! *********************************************************************
! * THIS SUBROUTINE UPDATES PHOTORATES AND ARB EMISSIONS RATES FOR    *
! * EACH TIME-STEP. PHOTORATES ARE INCLUDED IN FIRST AND PARTIAL      *
! * DERIVATIVE EQUATIONS WHILE EMISSIONS RATES ARE INCLUDED IN FIRST  *
! * DERIVATE EQUATIONS ONLY. SINCE THE EMISSIONS RATES ARE CONSTANT   * 
! * FOR A GIVEN TIME STEP AND LOCATION (ALTHOUGH THEY CHANGE EACH     *
! * TIME STEP AND LOCATION, THEY ARE PUT INTO THE FIRST DERIVATIVE    *
! * TERM OF SUBFUN.F ONLY (NOT INTO PARTIAL DERIVATIVE TERMS. EVERY   *
! * INTEGRATION TIME-STEP, EMISSIONS ARE RECALCULATED.                * 
! *********************************************************************
!
! *********************************************************************
! * UPDATE PHOTO-RATES AND OTHER PARMETERS BECAUSE THE TIME CHANGED.  *
! * NOTE THAT A TIME CHANGE COULD CORRESPOND TO EITHER A SUCCESSFUL   *
! * OR FAILED STEP                                                    * 
! *********************************************************************
! RRATE    = PRATE1 + XELAPS * (PRATE - PRATE1) 
! XELAPS   = ELAPSED TIME DURING INTERVAL  
! IFPRAT   = 1: USE SCALED PHOTORATES FROM photrate.dat (ITESTGEAR.EQ.0) 
!          = 0: USE PHOTORATES FROM globchem.dat (ITESTGEAR > 0)  
!
! *********************************************************************
! **************          UPDATE PHOTORATES             *************** 
! ****************** INTERPOLATE BETWEEN TWO VALUES *******************
! *********************************************************************
!
      ! Local variables
      INTEGER J,NKN,KLOOP,I,NK,NH,ISPC1,ISPC2,ISPC3

      REAL*8 TOFDAY,HOURANG,SINFUNC
!
! *********************************************************************
! *      SET RATES WHERE PHOTOREACTION HAS NO ACTIVE LOSS TERM        *
! *********************************************************************
! JOLD = MAPPL(JOLD) FOR INACTIVE SPECIES 
!
      DO 80 I            = 1, NOLOSP(NCSP) 
         NK                = NKNLOSP(I,NCS)
         NKN               = NEWFOLD(NK,NCS) 
         NH                = NKN + NALLR 
         DO 79 KLOOP       = 1, KTLOOP
            TRATE(KLOOP,NKN) =  RRATE(KLOOP,NKN)
            TRATE(KLOOP,NH)  = -RRATE(KLOOP,NKN)
 79      CONTINUE
 80   CONTINUE
!
! *********************************************************************
! *              PRINT OUT CHEMICAL RATES AND STOP                    *
! *********************************************************************
!
      IF (IPRATES.EQ.1) THEN
      if ( jlooplo == 744 ) then
       DO 90 I           = 1, NALLRAT(NCS)
        NK               = NCEQUAT(I,NCS)
        NKN              = NEWFOLD(NK,NCS) 
        ISPC1            = IRM(1,NK,NCS) 
        ISPC2            = IRM(2,NK,NCS) 
        ISPC3            = IRM(3,NK,NCS)
        IF (ISPC3.LT.0)          ISPC3 = 0 
        IF (ISPC1.GT.NSPEC(NCS)) ISPC1 = 0
        IF (ISPC2.GT.NSPEC(NCS)) ISPC2 = 0
        IF (ISPC3.GT.NSPEC(NCS)) ISPC3 = 0
        WRITE(6,95)I,NK,NKN,NAMENCS(ISPC1,NCS), NAMENCS(ISPC2,NCS), 
     1             NAMENCS(ISPC3,NCS), RRATE(1,NKN)
 90    CONTINUE
       STOP
       endif
      ENDIF
 95   FORMAT(I3,1X,I3,1X,I3,1X,3A15,1X,1PE13.6)
!
! *********************************************************************
! ******************** END OF SUBROUTINE UPDATE.F *********************
! *********************************************************************
!
      RETURN
      END SUBROUTINE UPDATE

