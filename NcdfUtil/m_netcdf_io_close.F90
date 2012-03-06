! $Id: m_netcdf_io_close.F90,v 1.1 2009/08/04 14:52:04 bmy Exp $
!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_close
!
! !INTERFACE:
!
      module m_netcdf_io_close
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  Nccl
      public  Nccl_Noerr
!
! !DESCRIPTION: Routines to close a netCDF file.
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
! !IROUTINE: Nccl
!
! !INTERFACE:
!
      subroutine Nccl (ncid)
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
!!    ncid : netCDF file id
      integer, intent (in)   :: ncid
!
! !DESCRIPTION: Closes a netCDF file with file id ncid.
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
      integer             :: ierr
!
      ierr = Nf_Close (ncid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Nccl:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Nccl
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nccl_Noerr
!
! !INTERFACE:
!
      subroutine Nccl_Noerr (ncid)
!
      implicit none
!
      include "netcdf.inc"
!
! !INPUT PARAMETERS:
!!    ncid : netCDF file id
      integer, intent (in)   :: ncid
!
! !DESCRIPTION: Closes a netCDF file (with file id ncid) if it is open and
!  suppresses Ncclos error messages/exit if it is not.
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
     integer             :: ierr
!
      ierr = Nf_Close (ncid)

      return

      end subroutine Nccl_Noerr
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_close

