! $Id: readcshcn.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE READCSHCN( XSECT_HCN )

      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

      IMPLICIT NONE

C**********************************************************************
C                                                                     *
C  HARVARD TROPOSPHERIC CHEMISTRY MODULE FOR 3-D APPLICATIONS         *
C  by Larry Horowitz, Jinyou Liang, Gerry Gardner, Prof. Daniel Jacob *
C  of HARVARD UNIVERSITY    (Release V2.0)                            *
C**********************************************************************
C
C**************************************************************
C                                                              
C This block reads in 8-column data that jcx just coded in. The
C data is saved in XSECT[IWL, ISPEC, ITEMP, IBRCH].
C
C**************************************************************

#     include "CMN_SIZE"
#     include "comsol.h"

      ! Arguments
      REAL*8, INTENT(OUT) :: XSECT_HCN(MXWL)

      ! Local variables
      CHARACTER*8         :: CHARSPEC 

      
      INTEGER K, ISPEC, MTEMP,I,IBRCH,ITEMP
      
      !=================================================================
      ! READCSHCN begins here!
      !=================================================================
      OPEN( 40, FILE='8col.dat.hcn',STATUS='old', ERR=700 )

      ! Begin a reading loop
 100  CONTINUE

      ! Check if the # of wavelength is as specified
      READ( 40, 1001, END=800 ) CHARSPEC, IBRCH, MTEMP
      READ( 40, 1003, END=800 ) ( XSECT_HCN(K), K=1, MXWL )

      CLOSE( 40 )

      ! Exit
      RETURN
 
      ! Trap file open errors
 700  CONTINUE
      CALL ERROR_STOP( 'Open error in "8col.dat.hcn"', 'readcshcn.f' )

      ! Trap file read errors
 800  CONTINUE
      CALL ERROR_STOP( 'Read error in "8col.dat.hcn"', 'readcshcn.f' )

      ! FORMATS
 1001 FORMAT(A8, 43X, I1, 5X, I3, 5x, a8)
 1003 FORMAT(8E10.3)

      ! Return to calling program
      END SUBROUTINE READCSHCN
