! $Id: xltmmp.f,v 1.1 2009/11/30 20:11:57 ccarouge Exp $
      FUNCTION XLTMMP( I, J, IJLOOP ) RESULT( VALUE )
!
!******************************************************************************
!  The new XLTMMP passes the value of the DAO meterological field 
!  TS(IIPAR,JJPAR) back to the calling subroutine.  This preserves the 
!  functionality of the H/G/I CTM function XLTMMP. (bmy, 1/30/98, 8/4/05)
!
!  NOTES
!  (1 ) XLTMMP is written in Fixed-Form Fortran 90.
!  (2 ) I, J are the long/lat indices of the grid box.  IJLOOP is passed
!        in order to maintain compatibility with the H/G/I subroutines,
!        but is not used. 
!  (3 ) TS is passed to XLTMMP via the "CMN_TS" include file.
!  (4 ) Use C-preprocessor #include statement to include CMN_SIZE, which 
!        has IIPAR, JJPAR, LLPAR, IGLOB, JGLOB, LGLOB.
!  (4 ) Now reference TS from "dao_mod.f" instead of from common block
!        header file "CMN_TS". (bmy, 6/23/00) 
!  (5 ) Eliminated obsolete code from 6/23/00 (bmy, 8/31/00)
!  (6 ) Now declare XLTMMP as REAL*8 w/in program body.  Also updated 
!        comments. (bmy, 9/26/01)
!  (7 ) Remove obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (8 ) IJLOOP is now not declared optional...this facilitates compiling with
!        -C on Altix (psk, bmy, 7/20/04)
!  (9 ) Now make IJLOOP an optional argument; it's only kept for backwards
!        compatibility w/ older code (bmy, 8/4/05)
!******************************************************************************
! 
      ! References to F90 modules
      USE DAO_MOD, ONLY : TS
      
      IMPLICIT NONE

#     include "CMN_SIZE"

      ! Arguments
      INTEGER, INTENT(IN)           :: I, J
      INTEGER, INTENT(IN), OPTIONAL :: IJLOOP

      ! Function value
      REAL*8                        :: VALUE

      !=================================================================
      ! XLTMMP begins here!      
      !=================================================================
      VALUE = TS(I,J)

      ! Return to calling program
      END FUNCTION XLTMMP
