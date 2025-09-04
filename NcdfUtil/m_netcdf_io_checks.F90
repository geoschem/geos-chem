!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_checks.F90
!
! !INTERFACE:
!
module m_netcdf_io_checks
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public  Ncdoes_Udim_Exist
  public  Ncdoes_Var_Exist
  public  Ncdoes_Attr_Exist
!
! !DESCRIPTION: Routines to check if a netCDF file contains a specified
!  variable.
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
!-------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Udim_Exist
!
! !INTERFACE:
!
  function Ncdoes_Udim_Exist (ncid)
!
    use netCDF
!
! !INPUT PARAMETERS:
!!  ncid : netCDF file id to check
    integer, intent (in)   :: ncid
!
! !DESCRIPTION: Checks a given netCDF file to see if it contains an
!  unlimited dimension.
!\\
!\\
! !RETURN VALUE:
    logical :: Ncdoes_Udim_Exist
!
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
    integer :: ierr, udim_id

    Ncdoes_Udim_Exist = .false.
    ierr = NF90_Inquire(ncid, unlimitedDimId=udim_id)
    IF ( ierr /= NF90_NOERR ) Ncdoes_Udim_Exist = .true.

  end function Ncdoes_Udim_Exist
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Var_Exist
!
! !INTERFACE:
!
  function Ncdoes_Var_Exist (ncid, varname)
!
    use netCDF
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id       to check
!!  varname : netCDF variable name to check
    integer,           intent (in)   :: ncid
    character (len=*), intent (in)   :: varname
!
! !DESCRIPTION: Checks a given netCDF file to see if a given netCDF variable
!  exists in it.
!\\
!\\
! !RETURN VALUE:
    logical :: Ncdoes_Var_Exist
!
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
    integer :: ierr
    integer :: varid
!
    ierr = NF90_Inq_Varid(ncid, varname, varid)
    Ncdoes_Var_Exist = .false.
    if (ierr == NF90_NOERR) Ncdoes_Var_Exist = .true.

  end function Ncdoes_Var_Exist
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Attr_Exist
!
! !INTERFACE:
!
  function Ncdoes_Attr_Exist(ncid, varname, attname, attType)
!
    use netCDF
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id       to check
!!  varname : netCDF variable name to check
!!  attname : netCDF attribute name to check
    integer,           intent (in)   :: ncid
    character (len=*), intent (in)   :: varname
    character (len=*), intent (in)   :: attname
!
! !OUTPUT PARAMETERS:
!
!! attType  : Attribute type.  This value is will be set to one of the
!! following: NF_BYTE, NF_CHAR, NF_SHORT, NF_INT, NF_FLOAT, or NF_DOUBLE.
    INTEGER,           INTENT(OUT)   :: attType
!
! !DESCRIPTION: Checks a given netCDF file to see if a given netCDF variable
!  exists in it.
!\\
!\\
! !RETURN VALUE:
    logical :: Ncdoes_Attr_Exist
!
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
    INTEGER :: ierr, varId, attLen, attNum

    ! Init
    Ncdoes_Attr_Exist = .false.
    attType           = -1

    ! First check the variable
    ierr = NF90_Inq_Varid (ncid, varname, varid)

    ! Check the attribute if variable was found
    IF ( ierr == NF90_NOERR ) THEN
       ierr = NF90_Inquire_Attribute( ncId,    varId,  attName,  &
                                      attType, attLen, attNum   )
       IF ( ierr == NF90_NOERR ) THEN
          NcDoes_Attr_Exist = .TRUE.
       ENDIF
    ENDIF

    return

  end function Ncdoes_Attr_Exist
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Dim_Exist
!
! !INTERFACE:
!
  function Ncdoes_Dim_Exist (ncid, dimname )
!
    use netCDF
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id        to check
!!  dimname : netCDF dimenison name to check
    integer,           intent (in)   :: ncid
    character (len=*), intent (in)   :: dimname
!
! !DESCRIPTION: Checks a given netCDF file to see if a given netCDF variable
!  exists in it.
!\\
!\\
! !RETURN VALUE:
    logical :: Ncdoes_Dim_Exist
!
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    integer :: ierr
    integer :: dimid

    ! First check the variable
    ierr = NF90_Inq_Dimid(ncid, dimname, dimid)

    ! Check the attribute if variable was found
    Ncdoes_Dim_Exist = .false.
    if (ierr == NF90_NOERR) Ncdoes_Dim_Exist = .true.

    return

  end function Ncdoes_Dim_Exist
!EOC
end module m_netcdf_io_checks
