!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_mod.F
!
! !DESCRIPTION: Module CMN is the remnant of header file "CMN", which once
!  held many global variables, but now is reduced to only a couple.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR

      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      ! LPAUSE  stores the annual mean tropopause (bmy, 12/6/99)
      INTEGER :: LPAUSE(IIPAR,JJPAR)

      ! IFLX stores the flags for ND27 diagnostic (bmy, 12/6/99)
      ! This is mostly obsolete by now (bmy, 8/24/11)
      INTEGER :: IFLX(IIPAR,JJPAR)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%% NOTE: THIS MODULE IS DEPRECATED.  IT WILL BE REMOVED IN THE FUTURE  %%% 
!  %%% COMBINE THESE ARRAYS INTO ANOTHER MODULE FILE.  FOR NOW, KEEP IT.   %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!EOP
!-----------------------------------------------------------------------------
!BOC
      END MODULE CMN_MOD
!EOC
