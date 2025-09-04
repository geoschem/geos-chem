!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_get_dimlen
!
! !INTERFACE:
!
MODULE m_netcdf_io_get_dimlen
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Ncget_Dimlen
  PUBLIC :: Ncget_Unlim_Dimlen
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
  SUBROUTINE Ncget_Dimlen(ncid, dim_name, dim_len)
!
! !USES:
!
    use netCDF
    use m_do_err_out
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
    CHARACTER(len=512) :: err_msg
    INTEGER            :: dimid
    INTEGER            :: ierr

    ierr = NF90_Inq_Dimid(ncid, dim_name, dimid)

    IF (ierr /= NF90_NOERR ) THEN
       err_msg = 'In Ncget_Dimlen #1:  ' // TRIM(dim_name) // &
                 ', ' // NF90_Strerror (ierr)
       CALL Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
    ENDIF

    ierr = NF90_Inquire_Dimension(ncid, dimid, len=dim_len)

    IF (ierr /= NF90_NOERR) THEN
       err_msg = 'In Ncget_Dimlen #2:  ' // NF90_Strerror (ierr)
       CALL Do_Err_Out (err_msg, .true., 2, ncid, dimid, 0, 0.0d0, 0.0d0)
    ENDIF

  END SUBROUTINE Ncget_Dimlen
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
  SUBROUTINE Ncget_Unlim_Dimlen (ncid, udim_len)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
!
! !INPUT PARAMETERS:
!!  ncid     : netCDF file id
    INTEGER, INTENT(IN)  :: ncid
!
! !OUTPUT PARAMETERS:
!!  udim_len : netCDF unlimited dimension length
    INTEGER, INTENT(OUT) :: udim_len
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
    CHARACTER(len=512) :: err_msg
    INTEGER            :: ierr, udim_id

    udim_len = -1
    ierr = NF90_Inquire(ncid, unlimitedDimId=udim_id)
    IF ( ierr /= NF90_NOERR ) THEN
       ierr = NF90_Inquire_Dimension( ncid, udim_id, len=udim_len )
    ENDIF

  END SUBROUTINE Ncget_Unlim_Dimlen
!EOC
END MODULE m_netcdf_io_get_dimlen
