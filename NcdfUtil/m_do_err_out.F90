! $Id: m_do_err_out.F90,v 1.1 2009/08/04 14:52:04 bmy Exp $
!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_Do_Err_Outc
!
! !INTERFACE:
!
   module m_Do_Err_Out
!
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
   public  Do_Err_Out
!
! !DESCRIPTION: Provides a routine to print an error message and exit the code.
!\\
!\\
! !AUTHOR: 
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Err_Out
!
! !INTERFACE:
!
      subroutine Do_Err_Out  &
        (err_msg, err_do_stop, err_num_ints, err_int1, err_int2,  &
         err_num_reals, err_real1, err_real2)
!
      implicit none
!
! !INPUT PARAMETERS:
!!     err_msg       : error message to be printed out
!!     err_do_stop   : do stop on error?
!!     err_num_ints  : number of integers to be printed out (0, 1, or 2)
!!     err_int1      : integer 1 to print out
!!     err_int2      : integer 2 to print out
!!     err_num_reals : number of reals to be printed out (0, 1, or 2)
!!     err_real1     : real 1 to print out
!!     err_real2     : real 2 to print out
      character (len=*), intent(in) :: err_msg
      logical          , intent(in) :: err_do_stop
      integer          , intent(in) :: err_num_ints
      integer          , intent(in) :: err_int1
      integer          , intent(in) :: err_int2
      integer          , intent(in) :: err_num_reals
      real*8           , intent(in) :: err_real1
      real*8           , intent(in) :: err_real2
!
! !DESCRIPTION: Outputs error messages, and exits if requested.
!\\
!\\
! !AUTHOR: 
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      Write (6,*)
      Write (6,*) &
        '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

      Write (6,*) '!! ' // Trim (err_msg)

      if (err_num_ints == 1) then
         Write (6,*) '!! ', err_int1
      else if (err_num_ints == 2) then
         Write (6,*) '!! ', err_int1, err_int2
      end if

      if (err_num_reals == 1) then
         Write (6,*) '!! ', err_real1
      else if (err_num_reals == 2) then
         Write (6,*) '!! ', err_real1, err_real2
      end if

      Write (6,*) &
        '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      Write (6,*)

      if (err_do_stop) then
        stop "Code stopped from Do_Err_Out."
      end if

      return

      end subroutine Do_Err_Out
!EOC
!------------------------------------------------------------------------
end module m_Do_Err_Out
