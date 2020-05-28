!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_close.F90
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
!  10 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  10 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
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
!!  ncid : netCDF file id
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
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
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
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
!!  ncid : netCDF file id
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
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    integer             :: ierr
!
    ierr = Nf_Close (ncid)

    return

  end subroutine Nccl_Noerr
!EOC
end module m_netcdf_io_close

