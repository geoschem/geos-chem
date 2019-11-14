!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
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
! !REVISION HISTORY:
!  07 Nov 2011 - R. Yantosca - Also give the option to create a netCDF4 file
!  10 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  10 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
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
! !IROUTINE: Nccr_Wr
!
! !INTERFACE:
!
  subroutine Nccr_Wr (ncid, filname, WRITE_NC4)
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
!   ncid    : opened netCDF file id
!   filname : name of netCDF file to open for writing
    integer          , intent(in)   :: ncid
    character (len=*), intent(in)   :: filname
    LOGICAL, OPTIONAL, INTENT(IN)   :: WRITE_NC4
!
! !DESCRIPTION: Creates a netCDF file for writing and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REMARKS:
!  If the netCDF4 library is used, then the NF_CLOBBER flag will write
!  a classic (i.e. netCDF3) file.  Use OR(NF_NETCDF4,NF_CLASSIC_MODEL) to
!  create netCDF-4 file that supports compression and uses "classic" netcdf data model
!  (no groups, no user-defined types)
!
! !REVISION HISTORY:
!  Initial code.
!  07 Nov 2011 - R. Yantosca - Also give the option to create a netCDF4 file
!                              by passing the optional WRITE_NC4 argument
!  17 Feb 2017 - C. Holmes   - Use netCDF-4 classic model for netCDF-4 files
!  01 Mar 2017 - R. Yantosca - Add an #ifdef to enable netCDF4 compression
!                              only if the library has nf_def_var_deflate
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
#if defined( NC_HAS_COMPRESSION )
       mode = IOR( NF_NETCDF4, NF_CLASSIC_MODEL )         ! netCDF4 file
       ierr = Nf_Create (filname, mode, ncid)             !  w/ compression
#else
       ierr = Nf_Create (filname, NF_64BIT_OFFSET, ncid)  ! netCDF4 file
                                                          !  w/o compression
#endif
    ELSE
       ierr = Nf_Create (filname, NF_CLOBBER, ncid)       ! netCDF3 file
    ENDIF

    if (ierr /= NF_NOERR) then
       err_msg = 'In Nccr_Wr, cannot create:  ' // Trim (filname)
       call Do_Err_Out (err_msg, .true., 0, 0, 0, 0 , 0.0d0, 0.0d0)
    end if

    return

  end subroutine Nccr_Wr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
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
    use m_do_err_out
!
    implicit none
!
    include "netcdf.inc"
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
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=128) :: err_msg
    integer             :: ierr
!
    ierr = Nf_Sync (ncid)

    if (ierr /= NF_NOERR) then
       err_msg = 'In Ncdo_Sync:  ' // Nf_Strerror (ierr)
       call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
    end if

    return

  end subroutine Ncdo_Sync
!EOC
end module m_netcdf_io_create
