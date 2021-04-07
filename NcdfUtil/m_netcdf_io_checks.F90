!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
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
! !REVISION HISTORY:
!  10 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  10 Jul 2014 - R. Yantosca - Cosmetic changes to ProTeX headers
!EOP
!-------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Udim_Exist
!
! !INTERFACE:
!
  function Ncdoes_Udim_Exist (ncid)
!
    implicit none
!
    include "netcdf.inc"
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
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    integer :: ierr
    integer :: udimid
!
    ierr = Nf_Inq_Unlimdim (ncid, udimid)

    if (ierr == NF_NOERR) then
       Ncdoes_Udim_Exist = .true.
    else
       Ncdoes_Udim_Exist = .false.
    end if

    return

  end function Ncdoes_Udim_Exist
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Var_Exist
!
! !INTERFACE:
!
  function Ncdoes_Var_Exist (ncid, varname)
!
    implicit none
!
    include "netcdf.inc"
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
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    integer :: ierr
    integer :: varid
!
    ierr = Nf_Inq_Varid (ncid, varname, varid)

    if (ierr == NF_NOERR) then
       Ncdoes_Var_Exist = .true.
    else
       Ncdoes_Var_Exist = .false.
    end if

    return

  end function Ncdoes_Var_Exist
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Attr_Exist
!
! !INTERFACE:
!
  function Ncdoes_Attr_Exist (ncid, varname, attname, attType)
!
    implicit none
!
    include "netcdf.inc"
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
!  Initial code.
!  03 Oct 2014 - C.Keller - Now check for int, real and character attributes
!  20 Feb 2015 - R. Yantosca - Now use NF_ATT_INQ function, it's more robust
!  20 Feb 2015 - R. Yantosca - Now return attribute type to calling routine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    integer :: ierr
    integer :: varid
    INTEGER :: attLen

    ! Init
    Ncdoes_Attr_Exist = .false.
    attType           = -1

    ! First check the variable
    ierr = Nf_Inq_Varid (ncid, varname, varid)

    ! Check the attribute if variable was found
    IF ( ierr == NF_NOERR ) THEN
       ierr = Nf_Inq_Att( ncId, varId, attName, attType, attLen )
       IF ( ierr == NF_NOERR ) THEN
          NcDoes_Attr_Exist = .TRUE.
       ENDIF
    ENDIF

    return

  end function Ncdoes_Attr_Exist
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Dim_Exist
!
! !INTERFACE:
!
  function Ncdoes_Dim_Exist (ncid, dimname )
!
    implicit none
!
    include "netcdf.inc"
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
!  Initial code.
!
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    integer :: ierr
    integer :: dimid

    ! First check the variable
    ierr = Nf_Inq_Dimid (ncid, dimname, dimid)

    ! Check the attribute if variable was found
    if (ierr == NF_NOERR) then
       Ncdoes_Dim_Exist = .true.
    else
       Ncdoes_Dim_Exist = .false.
    end if

    return

  end function Ncdoes_Dim_Exist
!EOC
end module m_netcdf_io_checks
