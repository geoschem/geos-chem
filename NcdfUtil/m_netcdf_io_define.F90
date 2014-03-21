! $Id: m_netcdf_io_define.F90,v 1.1 2009/08/04 14:52:04 bmy Exp $
!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_define
!
! !INTERFACE:
!
      module m_netcdf_io_define
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  NcDef_dimension
      public  NcDef_variable
      public  NcSetFill
      public  NcEnd_def

      PUBLIC :: NcDef_glob_attributes
      INTERFACE NcDef_glob_attributes
         MODULE PROCEDURE NcDef_glob_attributes_c
         MODULE PROCEDURE NcDef_glob_attributes_i
         MODULE PROCEDURE NcDef_glob_attributes_r4
         MODULE PROCEDURE NcDef_glob_attributes_r8
         MODULE PROCEDURE NcDef_glob_attributes_i_arr
         MODULE PROCEDURE NcDef_glob_attributes_r4_arr
         MODULE PROCEDURE NcDef_glob_attributes_r8_arr
      END INTERFACE

      PUBLIC :: NcDef_var_attributes
      INTERFACE NcDef_var_attributes
         MODULE PROCEDURE NcDef_var_attributes_c
         MODULE PROCEDURE NcDef_var_attributes_i
         MODULE PROCEDURE NcDef_var_attributes_r4
         MODULE PROCEDURE NcDef_var_attributes_r8
         MODULE PROCEDURE NcDef_var_attributes_i_arr
         MODULE PROCEDURE NcDef_var_attributes_r4_arr
         MODULE PROCEDURE NcDef_var_attributes_r8_arr
      END INTERFACE
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: NcDef_glob_attributes_c
      PRIVATE :: NcDef_glob_attributes_i
      PRIVATE :: NcDef_glob_attributes_r4
      PRIVATE :: NcDef_glob_attributes_r8
      PRIVATE :: NcDef_glob_attributes_i_arr
      PRIVATE :: NcDef_glob_attributes_r4_arr
      PRIVATE :: NcDef_glob_attributes_r8_arr
      PRIVATE :: NcDef_var_attributes_c
      PRIVATE :: NcDef_var_attributes_i
      PRIVATE :: NcDef_var_attributes_r4
      PRIVATE :: NcDef_var_attributes_r8
      PRIVATE :: NcDef_var_attributes_i_arr
      PRIVATE :: NcDef_var_attributes_r4_arr
      PRIVATE :: NcDef_var_attributes_r8_arr
!
! !DESCRIPTION: Provides netCDF utility routines to define dimensions, 
!  variables and attributes.
!\\
!\\
! !AUTHOR: 
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!  26 Sep 2013 - R. Yantosca - Add routines to save attributes of different
!                              numerical types
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_dimension
!
! !INTERFACE:
!
      subroutine NcDef_dimension(ncid,name,len,id)
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
!!    ncid  : netCDF file id
!!    name  : dimension name
!!    len   : dimension number
      character (len=*), intent(in) :: name
      integer,           intent(in) :: ncid, len
!
! !OUTPUT PARAMETERS:
!!    id    : dimension id
      integer,           intent(out) :: id
!
! !DESCRIPTION: Defines dimension.
!\\
!\\
! !AUTHOR: 
!  Jules Kouatchou and Maharaj Bhat
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer :: ierr
!
      ierr = Nf_Def_Dim (ncid, name, len, id)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Def_Dim: can not define dimension : '// Trim (name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return
      end subroutine NcDef_dimension
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_variable
!
! !INTERFACE:
!
      subroutine NcDef_variable(ncid,name,type,ndims,dims,var_id)
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
!
!!    ncid   : netCDF file id
!!    name   : name of the variable
!!    type   : type of the variable 
!!             (NF_FLOAT, NF_CHAR, NF_INT, NF_DOUBLE, NF_BYTE, NF_SHORT)
!!    ndims  : number of dimensions of the variable
!!    dims   : netCDF dimension id of the variable
!!    varid  : netCDF varid id

      character (len=*), intent(in) :: name
      integer,           intent(in) :: ncid, ndims, var_id
      integer,           intent(in) :: dims(ndims)
      integer,           intent(in) :: type
!
! !DESCRIPTION: Defines a netCDF variable.
!\\
!\\
! !AUTHOR: 
!  Jules Kouatchou and Maharaj Bhat
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer ::  ierr
!
      ierr = Nf_Def_Var (ncid, name, type, ndims, dims, var_id)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Def_Var: can not define variable : '// Trim (name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_variable
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes
!
! !INTERFACE:
!
      subroutine NcDef_var_attributes_c(ncid,var_id,att_name,att_val)
!
! !USES:
!
      use m_do_err_out
!
      implicit none
      include 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    var_id  : netCDF variable id
!!    att_name: attribute name
!!    att_val : attribute value
      character (len=*), intent(in) :: att_name, att_val
      integer,           intent(in) :: ncid, var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: CHARACTER.
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  Initial code.
!  26 Sep 2013 - R. Yantosca - Renamed to NcDef_var_attributes_c and made
!                              into a PRIVATE array so we can overload it
!EOP
!=====-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             ::  mylen, ierr
!
      mylen = len(att_val)
      ierr = Nf_Put_Att_Text (ncid, var_id, att_name, mylen, att_val)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_var_attributes_c: can not define attribute : ' // &
                   Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_var_attributes_c
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_i
!
! !INTERFACE:
!
      subroutine NcDef_var_attributes_i(ncid,var_id,att_name,att_val)
!
! !USES:
!
      use m_do_err_out
!
      implicit none
      include 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    var_id  : netCDF variable id
!!    att_name: attribute name
!!    att_val : attribute value
      INTEGER,           INTENT(IN) :: att_val
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid, var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: INTEGER.
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = 1
      ierr  = Nf_Put_Att_Real( ncid,   var_id, att_name, &
                               NF_INT, mylen,  att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_var_attributes_i: can not define attribute : ' // &
                    Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_var_attributes_i
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r4
!
! !INTERFACE:
!
      subroutine NcDef_var_attributes_r4(ncid,var_id,att_name,att_val)
!
! !USES:
!
      use m_do_err_out
!
      implicit none
      include 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    var_id  : netCDF variable id
!!    att_name: attribute name
!!    att_val : attribute value
      REAL*4,            INTENT(IN) :: att_val
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid, var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: REAL*4.
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = 1
      ierr  = Nf_Put_Att_Real( ncid,     var_id, att_name, &
                               NF_FLOAT, mylen,  att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_var_attributes_r4: can not define attribute : ' // &
                    Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_var_attributes_r4
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r8
!
! !INTERFACE:
!
      subroutine NcDef_var_attributes_r8(ncid,var_id,att_name,att_val)
!
! !USES:
!
      use m_do_err_out
!
      implicit none
      include 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    var_id  : netCDF variable id
!!    att_name: attribute name
!!    att_val : attribute value
      REAL*8,            INTENT(IN) :: att_val
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid, var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: REAL*4.
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  20 Sep 2013 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             ::  mylen, ierr
!
      mylen = 1
      ierr  = Nf_Put_Att_Double( ncid,      var_id, att_name, &
                                 NF_DOUBLE, mylen,  att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_var_attributes_r8: can not define attribute : ' // &
                    Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_var_attributes_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_i_arr
!
! !INTERFACE:
!
      subroutine NcDef_var_attributes_i_arr(ncid,var_id,att_name,att_val)
!
! !USES:
!
      use m_do_err_out
!
      implicit none
      include 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    var_id  : netCDF variable id
!!    att_name: attribute name
!!    att_val : attribute value
      INTEGER,           INTENT(IN) :: att_val(:)
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid, var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: INTEGER vector.
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = SIZE( att_val )
      ierr  = Nf_Put_Att_Real( ncid,   var_id, att_name, &
                               NF_INT, mylen,  att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_var_attributes_i_arr: can not define attribute : ' &
                    // Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_var_attributes_i_arr
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r4_arr
!
! !INTERFACE:
!
      subroutine NcDef_var_attributes_r4_arr(ncid,var_id,att_name,att_val)
!
! !USES:
!
      use m_do_err_out
!
      implicit none
      include 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    var_id  : netCDF variable id
!!    att_name: attribute name
!!    att_val : attribute value
      REAL*4,            INTENT(IN) :: att_val(:)
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid, var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: REAL*4 vector
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = SIZE( att_val )
      ierr  = Nf_Put_Att_Real( ncid,     var_id, att_name, &
                               NF_FLOAT, mylen,  att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_var_attributes_r4_arr: can not define attribute : ' &
                    // Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_var_attributes_r4_arr
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r8_arr
!
! !INTERFACE:
!
      subroutine NcDef_var_attributes_r8_arr(ncid,var_id,att_name,att_val)
!
! !USES:
!
      use m_do_err_out
!
      implicit none
      include 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    var_id  : netCDF variable id
!!    att_name: attribute name
!!    att_val : attribute value
      REAL*8,            INTENT(IN) :: att_val(:)
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid, var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: REAL*8 vector
!\\
!\\
! !AUTHOR: 
!  Jules Kouatchou and Maharaj Bhat
!
! !REVISION HISTORY:
!  20 Sep 2013 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             ::  mylen, ierr
!
      mylen = size( att_val )
      ierr  = Nf_Put_Att_Double( ncid,      var_id, att_name, &
                                 NF_DOUBLE, mylen,  att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_var_attributes_r4_arr: can not define attribute : '&
                     // Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_var_attributes_r8_arr
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_c
!
! !INTERFACE:
!
      subroutine NcDef_glob_attributes_c(ncid,att_name,att_val)
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
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
      character (len=*), intent(in) :: att_val
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid
!
! !DESCRIPTION: Defines global attributes of type: CHARACTER
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca( based on code by Jules Kouatchou)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             ::  mylen, ierr
!
      mylen = len(att_val)
      ierr = Nf_Put_Att_Text (ncid, NF_GLOBAL, att_name, mylen, att_val)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_glob_attributes_c: can not define attribute : ' // &
                   Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_glob_attributes_c
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_i
!
! !INTERFACE:
!
      subroutine NcDef_glob_attributes_i(ncid,att_name,att_val)
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
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
      INTEGER,           intent(in) :: att_val
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid
!
! !DESCRIPTION: Defines global attributes of type: INTEGER
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca( based on code by Jules Kouatchou)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = 1
      ierr  = Nf_Put_Att_Int( ncid,   NF_GLOBAL, att_name, &
                              NF_INT, mylen,     att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_glob_attributes_i: can not define attribute : ' // &
                   Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_glob_attributes_i
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r4
!
! !INTERFACE:
!
      subroutine NcDef_glob_attributes_r4(ncid,att_name,att_val)
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
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
      REAL*4,            intent(in) :: att_val
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid
!
! !DESCRIPTION: Defines global attributes of type: REAL*4
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca( based on code by Jules Kouatchou)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = 1
      ierr  = Nf_Put_Att_Real( ncid,     NF_GLOBAL, att_name, &
                               NF_FLOAT, mylen,     att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_glob_attributes_r4: can not define attribute : ' // &
                   Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_glob_attributes_r4
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r8
!
! !INTERFACE:
!
      subroutine NcDef_glob_attributes_r8(ncid,att_name,att_val)
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
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
      REAL*8,            intent(in) :: att_val
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid
!
! !DESCRIPTION: Defines global attributes of type: REAL*4
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca( based on code by Jules Kouatchou)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = 1
      ierr  = Nf_Put_Att_Double( ncid,     NF_GLOBAL, att_name, &
                                 NF_FLOAT, mylen,     att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_glob_attributes_r8: can not define attribute : ' // &
                   Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_glob_attributes_r8
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_i_arr
!
! !INTERFACE:
!
      subroutine NcDef_glob_attributes_i_arr(ncid,att_name,att_val)
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
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
      INTEGER,           intent(in) :: att_val(:)
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid
! 
! !DESCRIPTION: Defines global attributes of type: INTEGER vector
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca( based on code by Jules Kouatchou)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = SIZE( att_val )
      ierr  = Nf_Put_Att_Int( ncid,   NF_GLOBAL, att_name, &
                              NF_INT, mylen,     att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_glob_attributes_i_arr: can not define attribute : ' &
                    // Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_glob_attributes_i_arr
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r4_arr
!
! !INTERFACE:
!
      subroutine NcDef_glob_attributes_r4_arr(ncid,att_name,att_val)
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
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
      REAL*4,            intent(in) :: att_val(:)
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid
!
! !DESCRIPTION: Defines global attributes of type: REAL*4 vector
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca( based on code by Jules Kouatchou)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = SIZE( att_val )
      ierr  = Nf_Put_Att_Real( ncid,     NF_GLOBAL, att_name, &
                               NF_FLOAT, mylen,     att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_glob_attributes_r4_arr: can not define attribute : ' &
                   // Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_glob_attributes_r4_arr
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r8_arr
!
! !INTERFACE:
!
      subroutine NcDef_glob_attributes_r8_arr(ncid,att_name,att_val)
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
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
      REAL*8,            intent(in) :: att_val(:)
      character (len=*), intent(in) :: att_name
      integer,           intent(in) :: ncid
!
! !DESCRIPTION: Defines global attributes of type: REAL*8 vector
!\\
!\\
! !AUTHOR: 
!  Bob Yantosca( based on code by Jules Kouatchou)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=128) :: err_msg
      integer             :: mylen, ierr
!
      mylen = SIZE( att_val )
      ierr  = Nf_Put_Att_Double( ncid,     NF_GLOBAL, att_name, &
                                 NF_FLOAT, mylen,     att_val )

      if (ierr.ne.NF_NOERR) then
         err_msg = 'NcDef_glob_attributes_r8_arr: can not define attribute : ' &
                    // Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_glob_attributes_r8_arr
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcSetFill
!
! !INTERFACE:
!
      subroutine NcSetFill(ncid,ifill,omode)
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
!
      integer,           intent(in) :: ncid, ifill,omode
!
! !DESCRIPTION: Sets fill method.
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
      character (len=128) :: err_msg
      integer             ::  mylen, ierr
!
      ierr = Nf_Set_Fill (ncid, NF_NOFILL, omode)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Put_Att_Text: Error in omode  '
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcSetFill
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcEnd_def
!
! !INTERFACE:
!
      subroutine NcEnd_def(ncid)
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
!
      integer,           intent(in) :: ncid
!
! !DESCRIPTION: Ends definitions of variables and their attributes.
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
      character (len=128) :: err_msg
      integer             ::  ierr
!
      ierr = Nf_Enddef (ncid)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Put_Att_Text: Error in closing global attribute'
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcEnd_def
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_define
