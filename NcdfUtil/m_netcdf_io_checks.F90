! $Id: m_netcdf_io_checks.F90,v 1.1 2009/08/04 14:52:04 bmy Exp $
!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_checks
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
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
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
!!    ncid : netCDF file id to check
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
!-------------------------------------------------------------------------
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
!-------------------------------------------------------------------------
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
!!    ncid    : netCDF file id       to check
!!    varname : netCDF variable name to check
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
!-------------------------------------------------------------------------
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
!-------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Attr_Exist
!
! !INTERFACE:
!
      function Ncdoes_Attr_Exist (ncid, varname, attname)
!
      implicit none
!
      include "netcdf.inc"
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id       to check
!!    varname : netCDF variable name to check
!!    attname : netCDF attribute name to check
      integer,           intent (in)   :: ncid
      character (len=*), intent (in)   :: varname
      character (len=*), intent (in)   :: attname
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
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: ierr
      integer :: varid
      integer :: tmpout

      ! Init
      Ncdoes_Attr_Exist = .false.

      ! First check the variable
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      ! Check the attribute if variable was found
      if (ierr == NF_NOERR) then

         ierr = Nf_Get_Att_Int( ncid, varid, attname, tmpout )

         if ( ierr == NF_NOERR ) then
            Ncdoes_Attr_Exist = .true.
         end if

      end if

      return

      end function Ncdoes_Attr_Exist
!EOC
!------------------------------------------------------------------------
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
!!    ncid    : netCDF file id        to check
!!    dimname : netCDF dimenison name to check
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
!-------------------------------------------------------------------------
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
!------------------------------------------------------------------------
end module m_netcdf_io_checks
