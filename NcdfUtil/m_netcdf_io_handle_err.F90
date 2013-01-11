! $Id: m_netcdf_io_handle_err.F90,v 1.1 2009/08/04 14:52:04 bmy Exp $
!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_handle_err
!
! !INTERFACE:
!
      module m_netcdf_io_handle_err
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  Nchandle_Err
!
! !DESCRIPTION: Provides a routine to handle error messages.
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
! !IROUTINE: Nchandle_Err
!
! !INTERFACE:
!
      subroutine Nchandle_Err (ierr)
!
! !USES:
!
      use m_do_err_out
!
      implicit none
!
      include "netcdf.inc"
!
! !INPUT PARAMETERS:
!     ierr : netCDF error number
      integer, intent (in)   :: ierr
!
! !DESCRIPTION: Handles netCDF errors. Prints out a message and then exit. 
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
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
!
      err_msg = 'In Nchandle_Err:  ' // Nf_Strerror (ierr)

      call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

      return

      end subroutine Nchandle_Err
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_handle_err

