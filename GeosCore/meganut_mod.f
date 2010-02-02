! $Id: meganut_mod.f,v 1.3 2010/02/02 16:57:52 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: MEGANUT\_MOD
!
! !DESCRIPTION: Module MEGANUT\_MOD contains functions used by MEGAN.
!\\
!\\
! !INTERFACE
!
      MODULE MEGANUT_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: XLTMMP
      PUBLIC :: XLPARDF
      PUBLIC :: XLPARDR
!
! !REVISION HISTORY
!  20 Nov 2009 - C. Carouge - Create the module with xltmmp, xlpardf and 
!                             xlpardr functions.
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: XLTMMP
!
! !DESCRIPTION: Function XLTMMP passes the value of the DAO meterological 
!  field TS(IIPAR,JJPAR) back to the calling subroutine.  This preserves the 
!  functionality of the H/G/I CTM function XLTMMP.  XLTMMP is written in 
!  Fixed-Form Fortran 90.  I, J are the long/lat indices of the grid box.  
!  IJLOOP is passed in order to maintain compatibility with the H/G/I 
!  subroutines, but is not used. 
!\\
!\\
! !INTERFACE:
      FUNCTION XLTMMP( I, J, IJLOOP ) RESULT( VALUE )
!
! !USES:
!
      USE DAO_MOD, ONLY : TS
      
#     include "CMN_SIZE"
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)           :: I, J
      INTEGER, INTENT(IN), OPTIONAL :: IJLOOP
!
! !REVISION HISTORY:
!
!                              Use C-preprocessor #include statement to 
!                              include CMN_SIZE, which has IIPAR, JJPAR, 
!                              LLPAR, IGLOB, JGLOB, LGLOB.
!  23 Jun 2000 - R. Yantosca - Now reference TS from "dao_mod.f" instead of 
!                              from common block header file "CMN_TS". 
!  31 Aug 2000 - R. Yantosca - Eliminated obsolete code from 6/23/00
!  26 Sep 2001 - R. Yantosca - Now declare XLTMMP as REAL*8 w/in program body.
!                              Also updated comments.
!  24 Oct 2001 - R. Yantosca - Remove obsolete commented out code from 9/01
!  20 Jul 2004 - R. Yantosca - IJLOOP is now not declared optional...this 
!                              facilitates compiling with -C on Altix
!  04 Aug 2005 - R. Yantosca - Now make IJLOOP an optional argument; it's only 
!                              kept for backwards compatibility w/ older code
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES
!
      ! Function value
      REAL*8                        :: VALUE

      !=================================================================
      ! XLTMMP begins here!      
      !=================================================================
      VALUE = TS(I,J)

      ! Return to calling program
      END FUNCTION XLTMMP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: XLPARDR
!
! !DESCRIPTION: Function XLPARDR passes the value of the DAO meterological 
!  field PARDR(IIPAR,JJPAR) back to the calling subroutine.  This preserves 
!  the functionality of the H/G/I CTM function PARDR.  I, J are the long/lat 
!  indices of the grid box. IJLOOP is passed in order to maintain compatibility
!  with the H/G/I subroutines, but is not used.
!\\
!\\
!!INTERFACE
!
      FUNCTION XLPARDR( I, J, IJLOOP ) RESULT( VALUE )
!
! !USES
!
      USE DAO_MOD, ONLY : PARDR 
      
#     include "CMN_SIZE"
!
! !INPUT PARAMETERS 
!
      INTEGER, INTENT(IN)           :: I, J
      INTEGER, INTENT(IN), OPTIONAL :: IJLOOP
!
! !REVISION HISTORY
!
! 20 Nov 2009 - M. Barkley - Original version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      REAL*8             :: VALUE

      !=================================================================
      ! XLTMMP begins here!      
      !=================================================================
      VALUE = PARDR(I,J)
 
      ! Return to calling program
      END FUNCTION XLPARDR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: XLPARDF
!
! !DESCRIPTION: Function XLPARDF passes the value of the DAO meterological 
!  field PARDF(IIPAR,JJPAR) back to the calling subroutine.  This preserves 
!  the functionality of the H/G/I CTM function PARDF.  I, J are the long/lat 
!  indices of the grid box. IJLOOP is passed in order to maintain compatibility
!  with the H/G/I subroutines, but is not used.
!\\
!\\
!!INTERFACE
!
      FUNCTION XLPARDF( I, J, IJLOOP ) RESULT( VALUE )
!
! !USES
!
      USE DAO_MOD, ONLY : PARDF 
      
#     include "CMN_SIZE"
!
! !INPUT PARAMETERS 
!
      INTEGER, INTENT(IN)           :: I, J
      INTEGER, INTENT(IN), OPTIONAL :: IJLOOP
!
! !REVISION HISTORY
!  20 Nov 2009 - M. Barkley - Original version
!!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Function value
      REAL*8             :: VALUE

      !=================================================================
      ! XLPARDF begins here!      
      !=================================================================
      VALUE = PARDF(I,J)
 
      ! Return to calling program
      END FUNCTION XLPARDF
!EOC
      END MODULE MEGANUT_MOD
