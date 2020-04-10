! $Id: m_netcdf_io_write.F90,v 1.1 2009/08/04 14:52:05 bmy Exp $
!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_write
!
! !INTERFACE:
!
      module m_netcdf_io_write
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      ! Public interface
      PUBLIC :: NcWr

      ! Private methods overloaded by public interface
      ! (see below for info about these routines & the arguments they take)
      INTERFACE NcWr
         MODULE PROCEDURE Ncwr_Scal_R4
         MODULE PROCEDURE Ncwr_Scal_R8
         MODULE PROCEDURE Ncwr_Scal_Int
         MODULE PROCEDURE Ncwr_1d_R8
         MODULE PROCEDURE Ncwr_1d_R4
         MODULE PROCEDURE Ncwr_1d_Int
         MODULE PROCEDURE Ncwr_1d_Char
         MODULE PROCEDURE Ncwr_2d_R8
         MODULE PROCEDURE Ncwr_2d_R4
         MODULE PROCEDURE Ncwr_2d_Int
         MODULE PROCEDURE Ncwr_2d_Char
         MODULE PROCEDURE Ncwr_3d_R8
         MODULE PROCEDURE Ncwr_3d_R4
         MODULE PROCEDURE Ncwr_3d_Int
         MODULE PROCEDURE Ncwr_4d_R8
         MODULE PROCEDURE Ncwr_4d_R4
         MODULE PROCEDURE Ncwr_4d_Int
         MODULE PROCEDURE Ncwr_5d_R8
         MODULE PROCEDURE Ncwr_5d_R4
         MODULE PROCEDURE Ncwr_6d_R8
         MODULE PROCEDURE Ncwr_6d_R4

      END INTERFACE
!
! !DESCRIPTION: Routines for writing variables in a netCDF file.
!\\
!\\
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Add REAL*8 and REAL*4 versions of all
!                              NCWR_* routines.
!  20 Dec 2011 - R. Yantosca - Added Ncwr_4d_Int
!  20 Dec 2011 - R. Yantosca - Make process more efficient by not casting
!                              to temporary variables before file write
!EOP
!-------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_Scal_R4
!
! !INTERFACE:
!
      subroutine Ncwr_Scal_R4(varwr_scal, ncid, varname)
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
!!     ncid       : netCDF file id to write variable to
!!     varname    : netCDF variable name
!!     varwr_scal : variable to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      real*4           , intent(in)   :: varwr_scal

!
! !DESCRIPTION: Writes out a netCDF real scalar variable.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!  25 Aug 2017 - R. Yantosca - Renamed to NcWr_Scal_R4 and takes a
!                              REAL*4 argument
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal_R4 #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Var_Real (ncid, varid, varwr_scal)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal+R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_Scal_R4
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_Scal_R8
!
! !INTERFACE:
!
      subroutine Ncwr_Scal_R8 (varwr_scal, ncid, varname)
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
!!     ncid       : netCDF file id to write variable to
!!     varname    : netCDF variable name
!!     varwr_scal : variable to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      real*8           , intent(in)   :: varwr_scal
!
! !DESCRIPTION: Writes out a netCDF real scalar variable.
!\\
!\\
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
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal_R8 #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Var_Double(ncid, varid, varwr_scal)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_Scal_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_Scal_Int
!
! !INTERFACE:
!
      subroutine Ncwr_Scal_Int (varwr_scali, ncid, varname)
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
!!    ncid       : netCDF file id to write variable to
!!    varname    : netCDF variable name
!!    varwr_scali : integer variable to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: varwr_scali
!
! !DESCRIPTION: Writes out a netCDF integer scalar variable.
!\\
!\\
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
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Var_Int (ncid, varid, varwr_scali)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal_Int #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_Scal_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_1d_R8
!
! !INTERFACE:
!
      subroutine Ncwr_1d_R8 (varwr_1d, ncid, varname, strt1d, cnt1d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt1d   : vector specifying the index in varwr_1d where
!!               the first of the data values will be written
!!    cnt1d    : varwr_1d dimension
!!    varwr_1d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
      real*8           , intent(in)   :: varwr_1d(cnt1d(1))
!
! !DESCRIPTION: Writes out a 1D netCDF real array and does some error
!  checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to Ncrd_1d_R8.  REAL*8 version.
!  20 Dec 2011 - R. Yantosca - Now write varwr_1d directly to file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_PUT_VARA_DOUBLE
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_R8 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Double (ncid, varid, strt1d, cnt1d, varwr_1d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_1d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_1d_R4
!
! !INTERFACE:
!
      subroutine Ncwr_1d_R4 (varwr_1d, ncid, varname, strt1d, cnt1d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt1d   : vector specifying the index in varwr_1d where
!!               the first of the data values will be written
!!    cnt1d    : varwr_1d dimension
!!    varwr_1d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
      real*4           , intent(in)   :: varwr_1d(cnt1d(1))
!
! !DESCRIPTION: Writes out a 1D netCDF real array and does some error
!  checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to Ncwr_1d_R4.  REAL*4 version.
!  20 Dec 2011 - R. Yantosca - Now write varwr_1d directly to file
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_R4 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Real (ncid, varid, strt1d, cnt1d, varwr_1d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_1d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_1d_Int
!
! !INTERFACE:
!
      subroutine Ncwr_1d_Int (varwr_1di, ncid, varname, strt1d, cnt1d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt1d   : vector specifying the index in varwr_1di where
!!               the first of the data values will be written
!!    cnt1d    : varwr_1di dimension
!!    varwr_1di : intger array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
      integer          , intent(in)   :: varwr_1di(cnt1d(1))
!
! !DESCRIPTION: Writes out a 1D netCDF integer array and does some error
! checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_Int #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Int (ncid, varid, strt1d, cnt1d, varwr_1di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_Int #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_1d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_2d_R8
!
! !INTERFACE:
!
      subroutine Ncwr_2d_R8 (varwr_2d, ncid, varname, strt2d, cnt2d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt2d   : vector specifying the index in varwr_2d where
!!               the first of the data values will be written
!!    cnt2d    : varwr_2d dimensions
!!    varwr_2d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
      real*8           , intent(in)   :: varwr_2d(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION: Writes out a 2D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to Ncwr_2d_R8.  REAL*8 version.
!  20 Dec 2011 - R. Yantosca - Now write varwr_2d directly to file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_PUT_VARA_DOUBLE
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_R8 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Double (ncid, varid, strt2d, cnt2d, varwr_2d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_2d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_2d_R4
!
! !INTERFACE:
!
      subroutine Ncwr_2d_R4 (varwr_2d, ncid, varname, strt2d, cnt2d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt2d   : vector specifying the index in varwr_2d where
!!               the first of the data values will be written
!!    cnt2d    : varwr_2d dimensions
!!    varwr_2d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
      real*4           , intent(in)   :: varwr_2d(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION: Writes out a 2D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to Ncwr_2d_R4.  REAL*4 version.
!  20 Dec 2011 - R. Yantosca - Now write varwr_2d directly to file
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_R4 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Real (ncid, varid, strt2d, cnt2d, varwr_2d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_2d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_2d_Int
!
! !INTERFACE:
!
      subroutine Ncwr_2d_Int (varwr_2di, ncid, varname, strt2d, cnt2d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt2d   : vector specifying the index in varwr_2di where
!!               the first of the data values will be written
!!    cnt2d    : varwr_2di dimensions
!!    varwr_2di : intger array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
      integer          , intent(in)   :: varwr_2di(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION: Writes out a 2D netCDF integer array and does some error
!   checking.
!\\
!\\
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
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_Int #1:  ' // Trim (varname) //  &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Int (ncid, varid, strt2d, cnt2d, varwr_2di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_Int #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_2d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_3d_R8
!
! !INTERFACE:
!
      subroutine Ncwr_3d_R8 (varwr_3d, ncid, varname, strt3d, cnt3d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt3d   : vector specifying the index in varwr_3d where
!!               the first of the data values will be written
!!    cnt3d    : varwr_3d dimensions
!!    varwr_3d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
      real*8           , intent(in)   :: varwr_3d(cnt3d(1), cnt3d(2), cnt3d(3))
!
! !DESCRIPTION: Writes out a 3D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to NcWr_3d_R8.  REAL*8 version
!  20 Dec 2011 - R. Yantosca - Now write varwr_3d directly to file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_PUT_VARA_DOUBLE
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d_R8 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Double (ncid, varid, strt3d, cnt3d, varwr_3d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_3d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_3d_R4
!
! !INTERFACE:
!
      subroutine Ncwr_3d_R4 (varwr_3d, ncid, varname, strt3d, cnt3d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt3d   : vector specifying the index in varwr_3d where
!!               the first of the data values will be written
!!    cnt3d    : varwr_3d dimensions
!!    varwr_3d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
      real*4           , intent(in)   :: varwr_3d(cnt3d(1), cnt3d(2), cnt3d(3))
!
! !DESCRIPTION: Writes out a 3D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to NcWr_3d_R4.  REAL*4 version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d_R4 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Real (ncid, varid, strt3d, cnt3d, varwr_3d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_3d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_3d_Int
!
! !INTERFACE:
!
      subroutine Ncwr_3d_Int (varwr_3di, ncid, varname, strt3d, cnt3d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt3d   : vector specifying the index in varwr_3di where
!!               the first of the data values will be written
!!    cnt3d    : varwr_3di dimensions
!!    varwr_3di : intger array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
      integer          , intent(in)   :: varwr_3di(cnt3d(1), cnt3d(2), cnt3d(3))
!
! !DESCRIPTION: Writes out a 3D netCDF integer array and does some error
!  checking.
!\\
!\\
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
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr = Nf_Put_Vara_Int (ncid, varid, strt3d, cnt3d, varwr_3di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d_Int #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_3d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_4d_R8
!
! !INTERFACE:
!
      subroutine Ncwr_4d_R8 (varwr_4d, ncid, varname, strt4d, cnt4d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt4d   : vector specifying the index in varwr_4d where
!!               the first of the data values will be written
!!    cnt4d    : varwr_4d dimensions
!!    varwr_4d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt4d(4)
      integer          , intent(in)   :: cnt4d (4)
      real*8           , intent(in)   :: varwr_4d(cnt4d(1), cnt4d(2), &
                                                  cnt4d(3), cnt4d(4))
!
! !DESCRIPTION: Writes out a 4D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to NcWr_3d_R8.  REAL*8 version
!  20 Dec 2011 - R. Yantosca - Now write varwr_4d directly to file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_PUT_VARA_DOUBLE
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_4d_R8 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr = Nf_Put_Vara_Double (ncid, varid, strt4d, cnt4d, varwr_4d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_4d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_4d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_4d_R4
!
! !INTERFACE:
!
      subroutine Ncwr_4d_R4 (varwr_4d, ncid, varname, strt4d, cnt4d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt4d   : vector specifying the index in varwr_4d where
!!               the first of the data values will be written
!!    cnt4d    : varwr_4d dimensions
!!    varwr_4d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt4d(4)
      integer          , intent(in)   :: cnt4d (4)
      real*4           , intent(in)   :: varwr_4d(cnt4d(1), cnt4d(2), &
                                                  cnt4d(3), cnt4d(4))
!
! !DESCRIPTION: Writes out a 4D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to NcWr_3d_R8.  REAL*8 version
!  20 Dec 2011 - R. Yantosca - Now write varwr_4d directly to file
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_4d_R4 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr = Nf_Put_Vara_Real (ncid, varid, strt4d, cnt4d, varwr_4d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_4d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_4d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_3d_Int
!
! !INTERFACE:
!
      subroutine Ncwr_4d_Int (varwr_4di, ncid, varname, strt4d, cnt4d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt3d   : vector specifying the index in varwr_3di where
!!               the first of the data values will be written
!!    cnt3d    : varwr_3di dimensions
!!    varwr_3di : intger array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt4d(4)
      integer          , intent(in)   :: cnt4d (4)
      integer          , intent(in)   :: varwr_4di(cnt4d(1), cnt4d(2), &
                                                   cnt4d(3), cnt4d(4))
!
! !DESCRIPTION: Writes out a 3D netCDF integer array and does some error
!  checking.
!\\
!\\
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
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_4d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr = Nf_Put_Vara_Int (ncid, varid, strt4d, cnt4d, varwr_4di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_4d_Int #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_4d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_5d_R8
!
! !INTERFACE:
!
      subroutine Ncwr_5d_R8 (varwr_5d, ncid, varname, strt5d, cnt5d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt5d   : vector specifying the index in varwr_5d where
!!               the first of the data values will be written
!!    cnt5d    : varwr_5d dimensions
!!    varwr_5d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt5d(5)
      integer          , intent(in)   :: cnt5d (5)
      real*8           , intent(in)   :: varwr_5d(cnt5d(1), cnt5d(2), &
                                                  cnt5d(3), cnt5d(4), &
                                                  cnt5d(5))
!
! !DESCRIPTION: Writes out a 5D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to NcWr_5d_R8.  REAL*8 version
!  20 Dec 2011 - R. Yantosca - Now write varwr_5d directly to file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_PUT_VARA_DOUBLE
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_5d_R8 #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Double (ncid, varid, strt5d, cnt5d, varwr_5d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_5d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_5d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_5d_R4
!
! !INTERFACE:
!
      subroutine Ncwr_5d_R4 (varwr_5d, ncid, varname, strt5d, cnt5d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt5d   : vector specifying the index in varwr_5d where
!!               the first of the data values will be written
!!    cnt5d    : varwr_5d dimensions
!!    varwr_5d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt5d(5)
      integer          , intent(in)   :: cnt5d (5)
      real*4           , intent(in)   :: varwr_5d(cnt5d(1), cnt5d(2), &
                                                  cnt5d(3), cnt5d(4), &
                                                  cnt5d(5))
!
! !DESCRIPTION: Writes out a 5D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to NcWr_5d_R4.  REAL*4 version
!  20 Dec 2011 - R. Yantosca - Now write var5d_tmp directly to file

!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_5d_R4 #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Real (ncid, varid, strt5d, cnt5d, varwr_5d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_5d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_5d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_6d_R8
!
! !INTERFACE:
!
      subroutine Ncwr_6d_R8 (varwr_6d, ncid, varname, strt6d, cnt6d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt6d   : vector specifying the index in varwr_6d where
!!               the first of the data values will be written
!!    cnt6d    : varwr_6d dimensions
!!    varwr_6d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt6d(6)
      integer          , intent(in)   :: cnt6d (6)
      real*8           , intent(in)   :: varwr_6d(cnt6d(1), cnt6d(2), &
                                                  cnt6d(3), cnt6d(4), &
                                                  cnt6d(5), cnt6d(6))
!
! !DESCRIPTION: Writes out a 6D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to NcWr_6d_R8.  REAL*8 version.
!  20 Dec 2011 - R. Yantosca - Now write varwr_6d directly to file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_PUT_VARA_DOUBLE
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_6d_R8 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Double (ncid, varid, strt6d, cnt6d, varwr_6d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_6d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_6d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_6d_R4
!
! !INTERFACE:
!
      subroutine Ncwr_6d_R4 (varwr_6d, ncid, varname, strt6d, cnt6d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt6d   : vector specifying the index in varwr_6d where
!!               the first of the data values will be written
!!    cnt6d    : varwr_6d dimensions
!!    varwr_6d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt6d(6)
      integer          , intent(in)   :: cnt6d (6)
      real*4           , intent(in)   :: varwr_6d(cnt6d(1), cnt6d(2), &
                                                  cnt6d(3), cnt6d(4), &
                                                  cnt6d(5), cnt6d(6))
!
! !DESCRIPTION: Writes out a 6D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Renamed to NcWr_6d_R4.  REAL*4 version.
!  20 Dec 2011 - R. Yantosca - Now write varwr_6d directly to file
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_6d_R4 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Real (ncid, varid, strt6d, cnt6d, varwr_6d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_6d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_6d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_1d_Char
!
! !INTERFACE:
!
      subroutine Ncwr_1d_Char (varwr_1dc, ncid, varname, strt1d, cnt1d)
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
!!    ncid     : netCDF file id to write array output data to
!!    varname  : netCDF variable name for array
!!    strt1d   : vector specifying the index in varwr_1dc where
!!               the first of the data values will be written
!!    cnt1d    : varwr_1dc dimension
!!    varwr_1dc : intger array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
      character (len=1), intent(in)   :: varwr_1dc(cnt1d(1))
!
! !DESCRIPTION: Writes out a 1D netCDF character array and does some error
!   checking.
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
!BOC
!
! !LOCAL VARIABLES:
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_Char #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Text (ncid, varid, strt1d, cnt1d, varwr_1dc)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_Char #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_1d_Char
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_2d_Char
!
! !INTERFACE:
!
      subroutine Ncwr_2d_Char (char_2d, ncid, tvarname, strt2d, cnt2d)
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
!!    ncid     : netCDF file id to write text to
!!    tvarname : netCDF variable name for text
!!    strt2d   : vector specifying the index in char_2d where
!!               the first of the data values will be written
!!    cnt2d    : char_2d dimensions
!!    char_2d  : text to write
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: tvarname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
      character (len=1), intent(in)   :: char_2d(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION: Writes out a 2D netCDF character array and does some error
!  checking.
!\\
!\\
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
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: tvarid
!
      ierr = Nf_Inq_Varid (ncid, tvarname, tvarid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_Char #1:  ' // Trim (tvarname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr = Nf_Put_Vara_Text (ncid, tvarid, strt2d, cnt2d, char_2d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_Char #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, tvarid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncwr_2d_Char
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_write

