!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_get_dimlen
!
! !INTERFACE:
!
module m_netcdf_io_get_dimlen
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public  Ncget_Dimlen
  public  Ncget_Unlim_Dimlen
!
! !DESCRIPTION: Provides routines to obtain the length of a given dimension.
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
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncget_Dimlen
!
! !INTERFACE:
!
  subroutine Ncget_Dimlen (ncid, dim_name, dim_len )
!
! !USES:
!
    use m_do_err_out
!
    implicit none
!
    include 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  dim_name : netCDF dimension name
!!  ncid     : netCDF file id
    character (len=*), intent(in) :: dim_name
    integer,           intent(in) :: ncid
!
! !OUTPUT PARAMETERS:
!!  dim_len: netCDF dimension length
    integer,           intent(out)   :: dim_len
!
! !DESCRIPTION: Returns the length of a given netCDF dimension.
!               If err\_stop is set to FALSE, -1 is returned if
!               the given dimension cannot be found. Otherwise,
!               an error is prompted and the program stops.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: dimid
    integer             :: ierr

    ierr = Nf_Inq_Dimid  (ncid, dim_name, dimid)

    if (ierr /= NF_NOERR ) then
       err_msg = 'In Ncget_Dimlen #1:  ' // Trim (dim_name) // &
                 ', ' // Nf_Strerror (ierr)
       call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
    end if

    ierr = Nf_Inq_Dimlen (ncid, dimid, dim_len)

    if (ierr /= NF_NOERR) then
       err_msg = 'In Ncget_Dimlen #2:  ' // Nf_Strerror (ierr)
       call Do_Err_Out (err_msg, .true., 2, ncid, dimid, 0, 0.0d0, 0.0d0)
    end if

    return
  end subroutine Ncget_Dimlen
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncget_Unlim_Dimlen
!
! !INTERFACE:
!
  subroutine Ncget_Unlim_Dimlen (ncid, udim_len)
!
! !USES:
!
    use m_do_err_out
!
    implicit none
!
    include 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid     : netCDF file id
    integer,           intent(in) :: ncid
!
! !OUTPUT PARAMETERS:
!!  udim_len : netCDF unlimited dimension length
    integer,           intent(out) :: udim_len
!
! !DESCRIPTION: Returns the length of the unlimited netCDF dimension.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: ierr
    integer             :: udimid
!
    ierr = Nf_Inq_Unlimdim (ncid, udimid)

    if (ierr /= NF_NOERR) then
       err_msg = 'In Ncget_Unlim_Dimlen #1:  ' // Nf_Strerror (ierr)
       call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
    end if

    ierr = Nf_Inq_Dimlen (ncid, udimid, udim_len)

    if (ierr /= NF_NOERR) then
       err_msg = 'In Ncget_Unlim_Dimlen #2:  ' // Nf_Strerror (ierr)
       call Do_Err_Out (err_msg, .true., 2, ncid, udimid, 0, 0.0d0, 0.0d0)
    end if

    return

  end subroutine Ncget_Unlim_Dimlen
!EOC
end module m_netcdf_io_get_dimlen
