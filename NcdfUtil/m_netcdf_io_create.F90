!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_create.F90
!
! !INTERFACE:
!
module m_netcdf_io_create
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public  Nccr_Wr
  public  Ncdo_Sync
!
! !DESCRIPTION: Routines for creating and syncronizing netCDF files.
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
! !IROUTINE: Nccr_Wr
!
! !INTERFACE:
!
  subroutine Nccr_Wr (ncid, filname, WRITE_NC4)
!
! !USES:
!
    use netCDF
    use m_do_err_out
!
! !INPUT PARAMETERS:
!   ncid    : opened netCDF file id
!   filname : name of netCDF file to open for writing
    integer          , intent(INOUT) :: ncid
    character (len=*), intent(IN)    :: filname
    LOGICAL, OPTIONAL, INTENT(IN)    :: WRITE_NC4
!
! !DESCRIPTION: Creates a netCDF file for writing and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REMARKS:
!  If the netCDF4 library is used, then the NF90_CLOBBER flag will write
!  a classic (i.e. netCDF3) file.  Use OR(NF_NETCDF4,NF_CLASSIC_MODEL) to
!  create netCDF-4 file that supports compression and uses "classic"
!  netcdf data model (no groups, no user-defined types)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=128) :: err_msg
    integer             :: ierr
    INTEGER             :: mode
    LOGICAL             :: TMP_NC4
!
    ! Save the value of the optional WRITE_NC4 variable in
    ! a local shadow variable (bmy, 11/7/11)
    IF ( PRESENT( WRITE_NC4 ) ) THEN
       TMP_NC4 = WRITE_NC4
    ELSE
       TMP_NC4 = .FALSE.
    ENDIF

    IF ( TMP_NC4 ) THEN
#ifdef NC_HAS_COMPRESSION )
       mode = IOR( NF90_NETCDF4, NF90_CLASSIC_MODEL )       ! netCDF4 file
       ierr = NF90_Create(filname, mode, ncid)              !  w/ compression
#else
       ierr = NF90_Create(filname, NF90_64BIT_OFFSET, ncid) ! netCDF4 file
                                                            !  w/o compression
#endif
    ELSE
       ierr = NF90_Create(filname, NF90_CLOBBER, ncid)      ! netCDF3 file
    ENDIF

    if (ierr /= NF90_NOERR) then
       err_msg = 'In Nccr_Wr, cannot create:  ' // Trim (filname)
       call Do_Err_Out (err_msg, .true., 0, 0, 0, 0 , 0.0d0, 0.0d0)
    end if

    return

  end subroutine Nccr_Wr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncdo_Sync
!
! !INTERFACE:
!
  subroutine Ncdo_Sync (ncid)
!
! !USES:
!
    use netCDF
    use m_do_err_out
!
! !INPUT PARAMETERS:
!!  ncid : netCDF file id
    integer, intent(in)   :: ncid
!
! !DESCRIPTION: Synchronizes a netCDF file.
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
    character (len=128) :: err_msg
    integer             :: ierr
!
    ierr = NF90_Sync (ncid)

    if (ierr /= NF90_NOERR) then
       err_msg = 'In Ncdo_Sync:  ' // NF90_Strerror (ierr)
       call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
    end if

  end subroutine Ncdo_Sync
!EOC
end module m_netcdf_io_create
