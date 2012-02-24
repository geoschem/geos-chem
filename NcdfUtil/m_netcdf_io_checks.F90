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
!------------------------------------------------------------------------
end module m_netcdf_io_checks

