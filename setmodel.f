! $Id: setmodel.f,v 1.1 2003/06/30 20:26:06 bmy Exp $
      SUBROUTINE SETMODEL
!
!******************************************************************************
!  Subroutine SETMODEL computes the number of grid blocks that are needed.
!  (M. Jacobson 1997; bdf, bmy, 4/18/03)
!
!  NOTES:
!  (1 ) Originally, this routine also computed other meteorological parameters
!        such as horizontal & vertical coordinates, sun angles, etc.  These
!        are now computed elsewhere in GEOS-CHEM so this code has now been
!        removed.  The only code left is the code which determines the number
!        of grid blocks used for the parallelization.  Now force double-
!        precision with the "D" exponent. (bdf, bmy, 4/18/03)
!******************************************************************************
!
      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! SMVGEAR II arrays
C
C *********************************************************************
C ************        WRITTEN BY MARK JACOBSON (1993-4)    ************
C ***            (C) COPYRIGHT, 1993-4 BY MARK Z. JACOBSON          *** 
C ***               EXCEPT FOR DENOTED EXCERPTED PORTIONS           *** 
C ***                         (650) 650-6836                        *** 
C *********************************************************************
C
C SSSSSSS  EEEEEEE  TTTTTTT  M     M  OOOOOOO DDDDDD  EEEEEEE  L
C S        E           T     M M M M  O     O D     D E        L
C SSSSSSS  EEEEEEE     T     M  M  M  O     O D     D EEEEEEE  L
C       S  E           T     M     M  O     O D     D E        L 
C SSSSSSS  EEEEEEE     T     M     M  OOOOOOO DDDDDD  EEEEEEE  LLLLLLL
C
C *********************************************************************
C *       THIS SUBROUTINE INITIALIZES METEOROLOGICAL PARAMETERS       *
C *********************************************************************
C
      ! Local variables
      INTEGER :: IAVBLOK, IAVGSIZE, IREMAIN, JADD
C
C *********************************************************************
C *         DETERMINE HOW MANY BLOCKS OF GRID POINTS ARE NEEDED       *
C *********************************************************************
C
      KULOOP             = MIN(KULOOP,KBLOOP,NTLOOP) 
      NBLOCKS            = 1 + NTTLOOP / (KULOOP  + 0.0001d0)
      IAVBLOK            = 1 + NTTLOOP / (NBLOCKS + 0.0001d0)
      IAVGSIZE           = MIN0(IAVBLOK,KULOOP)
      JLOOPLO            = 0 
      IREMAIN            = NTTLOOP
C
      DO 260 KBLK        = 1, NBLOCKS 
         JADD            = MIN0(IAVGSIZE,IREMAIN)
         JLOFIXED(KBLK)  = JLOOPLO
         JHIFIXED(KBLK)  = JADD
         IREMAIN         = IREMAIN - JADD
         JLOOPLO         = JLOOPLO + JADD
 260  CONTINUE
C
C MAKE SURE MXBLOCK IS SUFFICIENTLY LARGE SINCE NBLOCKUSE CHANGES IN
C PHYSPROC.F
C
      IF (NBLOCKS+15.GT.MXBLOCK) THEN
         WRITE(6,*)'READER: NBLOCKS+15>MXBLOCKS ',NBLOCKS+15, MXBLOCK 
         STOP
      ENDIF
C
C *********************************************************************
C ******************** END OF SUBROUTINE SETMODEL.F *******************
C *********************************************************************
C
      RETURN
      END SUBROUTINE SETMODEL
