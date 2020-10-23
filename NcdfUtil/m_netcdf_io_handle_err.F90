!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_handle_err.F90
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
! !REMARKS:
!  This file is based on code from NASA/GSFC, SIVO, Code 610.3
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!-----------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
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
!   ierr : netCDF error number
    integer, intent (in)   :: ierr
!
! !DESCRIPTION: Handles netCDF errors. Prints out a message and then exit.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
!
    err_msg = 'In Nchandle_Err:  ' // Nf_Strerror (ierr)

    call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

    return

  end subroutine Nchandle_Err
!EOC
end module m_netcdf_io_handle_err

