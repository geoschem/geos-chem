! $Id: m_netcdf_io_read.F90,v 1.1 2009/08/04 14:52:05 bmy Exp $
!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_read
!
! !INTERFACE:
!
      MODULE m_netcdf_io_read
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      ! Public interface
      PUBLIC :: NcRd

      ! Private methods overloaded by public interface
      ! (see below for info about these routines & the arguments they take)
      INTERFACE NcRd
         MODULE PROCEDURE Ncrd_Scal
         MODULE PROCEDURE Ncrd_Scal_Int
         MODULE PROCEDURE Ncrd_1d_R8
         MODULE PROCEDURE Ncrd_1d_R4
         MODULE PROCEDURE Ncrd_1d_Int
         MODULE PROCEDURE Ncrd_1d_Char
         MODULE PROCEDURE Ncrd_2d_R8
         MODULE PROCEDURE Ncrd_2d_R4
         MODULE PROCEDURE Ncrd_2d_Int
         MODULE PROCEDURE Ncrd_2d_Char
         MODULE PROCEDURE Ncrd_3d_R8
         MODULE PROCEDURE Ncrd_3d_R4
         MODULE PROCEDURE Ncrd_3d_Int
         MODULE PROCEDURE Ncrd_4d_R8
         MODULE PROCEDURE Ncrd_4d_R4
         MODULE PROCEDURE Ncrd_4d_Int
         MODULE PROCEDURE Ncrd_5d_R8
         MODULE PROCEDURE Ncrd_5d_R4
         MODULE PROCEDURE Ncrd_6d_R8
         MODULE PROCEDURE Ncrd_6d_R4
         MODULE PROCEDURE Ncrd_7d_R8
         MODULE PROCEDURE Ncrd_7d_R4
      END INTERFACE
!
! !DESCRIPTION: Routines for reading variables in a netCDF file.
!\\
!\\
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!  03 Jul 2008 - R. Yantosca - Now overload all module methods with a
!                              single public interface.
!  26 Oct 2011 - R. Yantosca - Add REAL*8 and REAL*4 versions of all
!                              NCRD_* routines.
!  20 Dec 2011 - R. Yantosca - Added Ncwr_4d_Int
!  20 Dec 2011 - R. Yantosca - Make process more efficient by not casting
!                              to temporary variables after file read
!  04 Feb 2015 - C. Keller   - Added 7d reading routines.
!EOP
!-------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_Scal
!
! !INTERFACE:
!
      subroutine Ncrd_Scal (varrd_scal, ncid, varname)
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
!!    ncid       : netCDF file id to read variable from
!!    varname    : netCDF variable name
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
!
! !OUTPUT PARAMETERS:
!!    varrd_scal : variable to fill
      real*8           , intent(out)  :: varrd_scal
!
! !DESCRIPTION: Reads in a netCDF scalar variable.
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
      real*4              :: varrd_scal_tmp
!
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_Scal #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Var_Real   (ncid, varid, varrd_scal_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_Scal #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      varrd_scal = varrd_scal_tmp

      return

      end subroutine Ncrd_Scal
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_Scal_Int
!
! !INTERFACE:
!
      subroutine Ncrd_Scal_Int (varrd_scali, ncid, varname)
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
!!    ncid       : netCDF file id to read variable from
!!    varname    : netCDF variable name
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
!
! !OUTPUT PARAMETERS:
!!    varrd_scali : integer variable to fill
      integer          , intent(out)  :: varrd_scali
!
! !DESCRIPTION: Reads in a netCDF integer scalar variable.
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
        err_msg = 'In Ncrd_Scal_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Var_Int (ncid, varid, varrd_scali)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_Scal_Int #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncrd_Scal_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_1d_R8
!
! !INTERFACE:
!
      subroutine Ncrd_1d_R8 (varrd_1d, ncid, varname, strt1d, cnt1d,   &
                             err_stop, stat)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt1d   : vector specifying the index in varrd_1d where
!!               the first of the data values will be read
!!    cnt1d    : varrd_1d dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
      logical, optional, intent(in)   :: err_stop
!
! !OUTPUT PARAMETERS:
!!    varrd_1d : array to fill
      real*8           , intent(out)  :: varrd_1d(cnt1d(1))
      integer, optional, intent(out)  :: stat
!
! !DESCRIPTION: Reads in a 1D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_1d_R8.  REAL*8 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_1d directly from file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_GET_VARA_DOUBLE
!  24 Jan 2013 - C. Keller   - Added optional input arguments err_stop
!                              and stat
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
      logical             :: dostop

      ! set dostop flag
      if ( present ( err_stop ) ) then
        dostop = err_stop
      else
        dostop = .true.
      endif

      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        if ( dostop ) then
          err_msg = 'In Ncrd_1d_R8 #1:  ' // Trim (varname) // &
                     ', ' // Nf_Strerror (ierr)
          call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
        else
          varrd_1d(:) = -999d0
          if ( present ( stat ) ) stat = 1
          return
        end if
      end if

      ierr =  Nf_Get_Vara_Double (ncid, varid, strt1d, cnt1d, varrd_1d)

      if (ierr /= NF_NOERR) then
        if ( dostop ) then
          err_msg = 'In Ncrd_1d_R8 #2:  ' // Nf_Strerror (ierr)
          call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
        else
          varrd_1d(:) = -999d0
          if ( present ( stat ) ) stat = 2
          return
        endif
      end if

      ! set stat to 0 (= success)
      if ( present ( stat ) ) stat = 0

      end subroutine Ncrd_1d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_1d_R4
!
! !INTERFACE:
!
      subroutine Ncrd_1d_R4 (varrd_1d, ncid, varname, strt1d, cnt1d, &
                             err_stop, stat)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt1d   : vector specifying the index in varrd_1d where
!!               the first of the data values will be read
!!    cnt1d    : varrd_1d dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
      logical, optional, intent(in)   :: err_stop
!
! !OUTPUT PARAMETERS:
!!    varrd_1d : array to fill
      real*4           , intent(out)  :: varrd_1d(cnt1d(1))
      integer, optional, intent(out)  :: stat
!
! !DESCRIPTION: Reads in a 1D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_1d_R4.  REAL*4 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_1d directly from file
!  24 Jan 2013 - C. Keller   - Added optional input arguments err_stop
!                              and stat
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
      logical             :: dostop

      ! set dostop flag
      if ( present ( err_stop ) ) then
        dostop = err_stop
      else
        dostop = .true.
      endif

      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        if ( dostop ) then
           err_msg = 'In Ncrd_1d_R4 #1:  ' // Trim (varname) // &
                      ', ' // Nf_Strerror (ierr)
           call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
        else
          varrd_1d(:) = -999.0
          if ( present ( stat ) ) stat = 1
          return
        end if
      end if

      ierr =  Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, varrd_1d)

      if (ierr /= NF_NOERR) then
        if ( dostop ) then
          err_msg = 'In Ncrd_1d_R4 #2:  ' // Nf_Strerror (ierr)
          call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
        else
          varrd_1d(:) = -999.0
          if ( present ( stat ) ) stat = 2
          return
        endif
      end if

      ! set stat to 0 (= success)
      if ( present ( stat ) ) stat = 0
      return

      end subroutine Ncrd_1d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_1d_Int
!
! !INTERFACE:
!
      subroutine Ncrd_1d_Int (varrd_1di, ncid, varname, strt1d, cnt1d, &
                              err_stop, stat)
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
!
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt1d   : vector specifying the index in varrd_1di where
!!               the first of the data values will be read
!!    cnt1d    : varrd_1di dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
      logical, optional, intent(in)   :: err_stop
!
! !OUTPUT PARAMETERS:
!!    varrd_1di : intger array to fill
      integer          , intent(out)  :: varrd_1di(cnt1d(1))
      integer, optional, intent(out)  :: stat
!
! !DESCRIPTION: Reads in a 1D netCDF integer array and does some error
!  checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      character (len=512) :: err_msg
      integer             :: ierr
      integer             :: varid
      logical             :: dostop

      ! set dostop flag
      if ( present ( err_stop ) ) then
        dostop = err_stop
      else
        dostop = .true.
      endif

      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        if ( dostop ) then
          err_msg = 'In Ncrd_1d_Int #1:  ' // Trim (varname) // &
                    ', ' // Nf_Strerror (ierr)
          call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
        else
          varrd_1di(:) = -999
          if ( present ( stat ) ) stat = 1
          return
        end if
      end if

      ierr = Nf_Get_Vara_Int (ncid, varid, strt1d, cnt1d, varrd_1di)

      if (ierr /= NF_NOERR) then
        if ( dostop ) then
          err_msg = 'In Ncrd_1d_Int #2:  ' // Nf_Strerror (ierr)
          call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
        else
          varrd_1di(:) = -999
          if ( present ( stat ) ) stat = 2
          return
        endif
      end if

      ! set stat to 0 (= success)
      if ( present ( stat ) ) stat = 0

      return

      end subroutine Ncrd_1d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_2d_R8
!
! !INTERFACE:
!
      subroutine Ncrd_2d_R8 (varrd_2d, ncid, varname, strt2d, cnt2d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt2d   : vector specifying the index in varrd_2d where
!!               the first of the data values will be read
!!    cnt2d    : varrd_2d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
!
! !OUTPUT PARAMETERS:
!!    varrd_2d : array to fill
      real*8           , intent(out)  :: varrd_2d(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION: Reads in a 2D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_2d_R8.  REAL*8 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_2d directly from file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_GET_VARA_DOUBLE
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
        err_msg = 'In Ncrd_2d_R8 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Double (ncid, varid, strt2d, cnt2d, varrd_2d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_2d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_2d_R4
!
! !INTERFACE:
!
      subroutine Ncrd_2d_R4 (varrd_2d, ncid, varname, strt2d, cnt2d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt2d   : vector specifying the index in varrd_2d where
!!               the first of the data values will be read
!!    cnt2d    : varrd_2d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
!
! !OUTPUT PARAMETERS:
!!    varrd_2d : array to fill
      real*4           , intent(out)  :: varrd_2d(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION: Reads in a 2D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_2d_R4.  REAL*4 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_2d directly from file
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
        err_msg = 'In Ncrd_2d_R4 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Real (ncid, varid, strt2d, cnt2d, varrd_2d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_2d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_2d_Int
!
! !INTERFACE:
!
      subroutine Ncrd_2d_Int (varrd_2di, ncid, varname, strt2d, cnt2d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt2d   : vector specifying the index in varrd_2d where
!!               the first of the data values will be read
!!    cnt2d    : varrd_2di dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
!
! !OUTPUT PARAMETERS:
!!    varrd_2di : intger array to fill
      integer          , intent(out)  :: varrd_2di(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION: Reads in a 2D netCDF integer array and does some error
!  checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
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
        err_msg = 'In Ncrd_2d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Int (ncid, varid, strt2d, cnt2d, varrd_2di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d_Int #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_2d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_3d_R8
!
! !INTERFACE:
!
      subroutine Ncrd_3d_R8 (varrd_3d, ncid, varname, strt3d, cnt3d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt3d   : vector specifying the index in varrd_3d where
!!               the first of the data values will be read
!!    cnt3d    : varrd_3d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
!
! !OUTPUT PARAMETERS:
!!    varrd_3d : array to fill
      real*8           , intent(out)  :: varrd_3d(cnt3d(1), cnt3d(2), &
                                                  cnt3d(3))
!
! !DESCRIPTION: Reads in a 3D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_3d_R8.  REAL*8 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_3d directly from file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_GET_VARA_DOUBLE
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
        err_msg = 'In Ncrd_3d_R8 #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Double (ncid, varid, strt3d, cnt3d, varrd_3d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_3d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_3d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_3d_R4
!
! !INTERFACE:
!
      subroutine Ncrd_3d_R4 (varrd_3d, ncid, varname, strt3d, cnt3d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt3d   : vector specifying the index in varrd_3d where
!!               the first of the data values will be read
!!    cnt3d    : varrd_3d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
!
! !OUTPUT PARAMETERS:
!!    varrd_3d : array to fill
      real*4           , intent(out)  :: varrd_3d(cnt3d(1), cnt3d(2), &
                                                  cnt3d(3))
!
! !DESCRIPTION: Reads in a 3D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_3d_R4.  REAL*4 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_3d directly from file
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
        err_msg = 'In Ncrd_3d_R4 #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, varrd_3d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_3d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_3d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_3d_Int
!
! !INTERFACE:
!
      subroutine Ncrd_3d_Int (varrd_3di, ncid, varname, strt3d, cnt3d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt3d   : vector specifying the index in varrd_3d where
!!               the first of the data values will be read
!!    cnt3d    : varrd_3di dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
!
! !OUTPUT PARAMETERS:
!!    varrd_3di : intger array to fill
      integer          , intent(out)  :: varrd_3di(cnt3d(1), cnt3d(2), &
                                                   cnt3d(3))
!
! !DESCRIPTION: Reads in a 3D netCDF integer array and does some error
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
        err_msg = 'In Ncrd_3d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Int (ncid, varid, strt3d, cnt3d, varrd_3di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_3d_Int #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_3d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_4d_R8
!
! !INTERFACE:
!
      subroutine Ncrd_4d_R8 (varrd_4d, ncid, varname, strt4d, cnt4d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt4d   : vector specifying the index in varrd_4d where
!!               the first of the data values will be read
!!    cnt4d    : varrd_4d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt4d(4)
      integer          , intent(in)   :: cnt4d (4)
!
! !OUTPUT PARAMETERS:
!!    varrd_4d : array to fill
      real*8           , intent(out)  :: varrd_4d(cnt4d(1), cnt4d(2), &
                                                  cnt4d(3), cnt4d(4))
!
! !DESCRIPTION: Reads in a 4D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_4d_R8.  REAL*8 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_4d directly from file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_GET_VARA_DOUBLE
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
        err_msg = 'In Ncrd_4d_R8 #1:  ' // Trim (varname) // &
                    ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr =  Nf_Get_Vara_Double (ncid, varid, strt4d, cnt4d, varrd_4d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_4d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_4d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_4d_R4
!
! !INTERFACE:
!
      subroutine Ncrd_4d_R4 (varrd_4d, ncid, varname, strt4d, cnt4d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt4d   : vector specifying the index in varrd_4d where
!!               the first of the data values will be read
!!    cnt4d    : varrd_4d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt4d(4)
      integer          , intent(in)   :: cnt4d (4)
!
! !OUTPUT PARAMETERS:
!!    varrd_4d : array to fill
      real*4           , intent(out)  :: varrd_4d(cnt4d(1), cnt4d(2), &
                                                  cnt4d(3), cnt4d(4))
!
! !DESCRIPTION: Reads in a 4D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_4d_R4.  REAL*4 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_4d directly from file
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
        err_msg = 'In Ncrd_4d_R4 #1:  ' // Trim (varname) // &
                    ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  Nf_Get_Vara_Real (ncid, varid, strt4d, cnt4d, varrd_4d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_4d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_4d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_4d_Int
!
! !INTERFACE:
!
      subroutine Ncrd_4d_Int (varrd_4di, ncid, varname, strt4d, cnt4d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt3d   : vector specifying the index in varrd_3d where
!!               the first of the data values will be read
!!    cnt3d    : varrd_3di dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt4d(4)
      integer          , intent(in)   :: cnt4d (4)
!
! !OUTPUT PARAMETERS:
!!    varrd_3di : intger array to fill
      integer          , intent(out)  :: varrd_4di(cnt4d(1), cnt4d(2), &
                                                   cnt4d(3), cnt4d(4))
!
! !DESCRIPTION: Reads in a 3D netCDF integer array and does some error
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
        err_msg = 'In Ncrd_3d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Int (ncid, varid, strt4d, cnt4d, varrd_4di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_3d_Int #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_4d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_5d_R8
!
! !INTERFACE:
!
      subroutine Ncrd_5d_R8 (varrd_5d, ncid, varname, strt5d, cnt5d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt5d   : vector specifying the index in varrd_5d where
!!               the first of the data values will be read
!!    cnt5d    : varrd_5d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt5d(5)
      integer          , intent(in)   :: cnt5d (5)
!
! !OUTPUT PARAMETERS:
!!    varrd_5d : array to fill
      real*8         , intent(out)  :: varrd_5d(cnt5d(1), cnt5d(2), &
                                                cnt5d(3), cnt5d(4), &
                                                cnt5d(5))
!
! !DESCRIPTION: Reads in a 5D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_45_R8.  REAL*8 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_5d directly from file
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_GET_VARA_DOUBLE
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
        err_msg = 'In Ncrd_5d_R8 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Double (ncid, varid, strt5d, cnt5d, varrd_5d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_5d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_5d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_5d_R4
!
! !INTERFACE:
!
      subroutine Ncrd_5d_R4 (varrd_5d, ncid, varname, strt5d, cnt5d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt5d   : vector specifying the index in varrd_5d where
!!               the first of the data values will be read
!!    cnt5d    : varrd_5d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt5d(5)
      integer          , intent(in)   :: cnt5d (5)
!
! !OUTPUT PARAMETERS:
!!    varrd_5d : array to fill
      real*4          , intent(out)  :: varrd_5d(cnt5d(1), cnt5d(2), &
                                                cnt5d(3), cnt5d(4), &
                                                cnt5d(5))
!
! !DESCRIPTION: Reads in a 5D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_45_R4.  REAL*4 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_5d directly from file
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
        err_msg = 'In Ncrd_5d_R4 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Real (ncid, varid, strt5d, cnt5d, varrd_5d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_5d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_5d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_6d_R8
!
! !INTERFACE:
!
      subroutine Ncrd_6d_R8 (varrd_6d, ncid, varname, strt6d, cnt6d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt5d   : vector specifying the index in varrd_5d where
!!               the first of the data values will be read
!!    cnt5d    : varrd_5d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt6d(6)
      integer          , intent(in)   :: cnt6d (6)
!
! !OUTPUT PARAMETERS:
!!    varrd_5d : array to fill
      real*8         , intent(out)  :: varrd_6d(cnt6d(1), cnt6d(2), &
                                                cnt6d(3), cnt6d(4), &
                                                cnt6d(5), cnt6d(6))
!
! !DESCRIPTION: Reads in a 5D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Initial version
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_GET_VARA_DOUBLE
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
        err_msg = 'In Ncrd_6d_R8 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Double (ncid, varid, strt6d, cnt6d, varrd_6d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_6d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_6d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_6d_R4
!
! !INTERFACE:
!
      subroutine Ncrd_6d_R4 (varrd_6d, ncid, varname, strt6d, cnt6d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt5d   : vector specifying the index in varrd_5d where
!!               the first of the data values will be read
!!    cnt5d    : varrd_5d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt6d(6)
      integer          , intent(in)   :: cnt6d (6)
!
! !OUTPUT PARAMETERS:
!!    varrd_5d : array to fill
      real*4          , intent(out)  :: varrd_6d(cnt6d(1), cnt6d(2), &
                                                 cnt6d(3), cnt6d(4), &
                                                 cnt6d(5), cnt6d(6))
!
! !DESCRIPTION: Reads in a 5D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  26 Oct 2011 - R. Yantosca - Renamed to Ncrd_45_R4.  REAL*4 version.
!  20 Dec 2011 - R. Yantosca - Now read varrd_5d directly from file
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
        err_msg = 'In Ncrd_6d_R4 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Real   (ncid, varid, strt6d, cnt6d, varrd_6d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_6d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_6d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_7d_R8
!
! !INTERFACE:
!
      subroutine Ncrd_7d_R8 (varrd_7d, ncid, varname, strt7d, cnt7d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt7d   : vector specifying the index in varrd_7d where
!!               the first of the data values will be read
!!    cnt7d    : varrd_7d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt7d(7)
      integer          , intent(in)   :: cnt7d (7)
!
! !OUTPUT PARAMETERS:
!!    varrd_5d : array to fill
      real*8         , intent(out)  :: varrd_7d(cnt7d(1), cnt7d(2), &
                                                cnt7d(3), cnt7d(4), &
                                                cnt7d(5), cnt7d(6), &
                                                cnt7d(7))
!
! !DESCRIPTION: Reads in a 7D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Initial version
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_GET_VARA_DOUBLE
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
        err_msg = 'In Ncrd_7d_R8 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Double (ncid, varid, strt7d, cnt7d, varrd_7d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_7d_R8 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_7d_R8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_7d_R4
!
! !INTERFACE:
!
      subroutine Ncrd_7d_R4 (varrd_7d, ncid, varname, strt7d, cnt7d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt7d   : vector specifying the index in varrd_7d where
!!               the first of the data values will be read
!!    cnt7d    : varrd_7d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt7d(7)
      integer          , intent(in)   :: cnt7d (7)
!
! !OUTPUT PARAMETERS:
!!    varrd_7d : array to fill
      real*4          , intent(out)  :: varrd_7d(cnt7d(1), cnt7d(2), &
                                                 cnt7d(3), cnt7d(4), &
                                                 cnt7d(5), cnt7d(6), &
                                                 cnt7d(7))
!
! !DESCRIPTION: Reads in a 7D netCDF real array and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  20 Dec 2011 - R. Yantosca - Initial version
!  20 Dec 2011 - R. Yantosca - Now use netCDF function NF_GET_VARA_DOUBLE
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
        err_msg = 'In Ncrd_7d_R4 #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Real (ncid, varid, strt7d, cnt7d, varrd_7d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_7d_R4 #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_7d_R4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_1d_Char
!
! !INTERFACE:
!
      subroutine Ncrd_1d_Char (varrd_1dc, ncid, varname, strt1d, cnt1d)
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
!
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt1d   : vector specifying the index in varrd_1dc where
!!               the first of the data values will be read
!!    cnt1d    : varrd_1dc dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
!
! !OUTPUT PARAMETERS:
!!    varrd_1dc : intger array to fill
      character (len=1), intent(out)  :: varrd_1dc(cnt1d(1))
!
! !DESCRIPTION: Reads in a 1D netCDF character array and does some error
!  checking.
!\\
!\\ !AUTHOR:
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
        err_msg = 'In Ncrd_1d_Char #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Text (ncid, varid, strt1d, cnt1d, varrd_1dc)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_1d_Char #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_1d_Char
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_2d_Char
!
! !INTERFACE:
!
      subroutine Ncrd_2d_Char (varrd_2dc, ncid, varname, strt2d, cnt2d)
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
!!    ncid     : netCDF file id to read array input data from
!!    varname  : netCDF variable name for array
!!    strt2d   : vector specifying the index in varrd_2dc where
!!               the first of the data values will be read
!!    cnt2d    : varrd_2dc dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
!
! !OUTPUT PARAMETERS:
!!    varrd_2dc : charcter array to fill
      character        , intent(out)  :: varrd_2dc(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION: Reads in a 2D netCDF character array and does some error
!  checking.
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
        err_msg = 'In Ncrd_2d_Char #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Text (ncid, varid, strt2d, cnt2d, varrd_2dc)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d_Char #2:  ' // Nf_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      end subroutine Ncrd_2d_Char
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_read

