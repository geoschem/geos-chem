! $Id: setmodel.f,v 1.1 2009/09/16 14:06:09 bmy Exp $
      SUBROUTINE SETMODEL( am_I_Root )
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
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE CMN_SIZE_MOD  ! Size parameters
      USE COMODE_LOOP_MOD  ! SMVGEAR II arrays

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! *********************************************************************
! ************        WRITTEN BY MARK JACOBSON (1993-4)    ************
! ***            (C) COPYRIGHT, 1993-4 BY MARK Z. JACOBSON          *** 
! ***               EXCEPT FOR DENOTED EXCERPTED PORTIONS           *** 
! ***                         (650) 650-6836                        *** 
! *********************************************************************
!
! SSSSSSS  EEEEEEE  TTTTTTT  M     M  OOOOOOO DDDDDD  EEEEEEE  L
! S        E           T     M M M M  O     O D     D E        L
! SSSSSSS  EEEEEEE     T     M  M  M  O     O D     D EEEEEEE  L
!       S  E           T     M     M  O     O D     D E        L 
! SSSSSSS  EEEEEEE     T     M     M  OOOOOOO DDDDDD  EEEEEEE  LLLLLLL
!
! *********************************************************************
! *       THIS SUBROUTINE INITIALIZES METEOROLOGICAL PARAMETERS       *
! *********************************************************************
!
      ! Local variables
      INTEGER :: IAVBLOK, IAVGSIZE, IREMAIN, JADD
!
! *********************************************************************
! *         DETERMINE HOW MANY BLOCKS OF GRID POINTS ARE NEEDED       *
! *********************************************************************
!
      KULOOP             = MIN(KULOOP,KBLOOP,NTLOOP) 
      NBLOCKS            = 1 + NTTLOOP / (KULOOP  + 0.0001d0)
      IAVBLOK            = 1 + NTTLOOP / (NBLOCKS + 0.0001d0)
      IAVGSIZE           = MIN0(IAVBLOK,KULOOP)
      JLOOPLO            = 0 
      IREMAIN            = NTTLOOP
!
      DO 260 KBLK        = 1, NBLOCKS 
         JADD            = MIN0(IAVGSIZE,IREMAIN)
         JLOFIXED(KBLK)  = JLOOPLO
         JHIFIXED(KBLK)  = JADD
         IREMAIN         = IREMAIN - JADD
         JLOOPLO         = JLOOPLO + JADD
 260  CONTINUE
!
! MAKE SURE MXBLOCK IS SUFFICIENTLY LARGE SINCE NBLOCKUSE CHANGES IN
! PHYSPROC.F
!
      IF (NBLOCKS+15.GT.MXBLOCK) THEN
         IF ( am_I_Root ) THEN
            WRITE(6,*)'READER: NBLOCKS+15>MXBLOCKS ',
     &           NBLOCKS+15, MXBLOCK 
         ENDIF
         STOP
      ENDIF
!
! *********************************************************************
! ******************** END OF SUBROUTINE SETMODEL.F *******************
! *********************************************************************
!
      RETURN
      END SUBROUTINE SETMODEL
