!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_define.F90
!
! !INTERFACE:
!
MODULE m_netcdf_io_define
!
! !USES:
!
  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: NcDef_Dimension
  PUBLIC :: NcDef_Variable
  PUBLIC :: NcSetFill
  PUBLIC :: NcEnd_Def
  PUBLIC :: NcBegin_Def

  PUBLIC :: NcDef_glob_attributes
  INTERFACE NcDef_glob_attributes
     MODULE PROCEDURE NcDef_glob_attributes_c
     MODULE PROCEDURE NcDef_glob_attributes_i
     MODULE PROCEDURE NcDef_glob_attributes_r4
     MODULE PROCEDURE NcDef_glob_attributes_r8
     MODULE PROCEDURE NcDef_glob_attributes_i_arr
     MODULE PROCEDURE NcDef_glob_attributes_r4_arr
     MODULE PROCEDURE NcDef_glob_attributes_r8_arr
  END INTERFACE NcDef_glob_attributes

  PUBLIC :: NcDef_var_attributes
  INTERFACE NcDef_var_attributes
     MODULE PROCEDURE NcDef_var_attributes_c
     MODULE PROCEDURE NcDef_var_attributes_i
     MODULE PROCEDURE NcDef_var_attributes_r4
     MODULE PROCEDURE NcDef_var_attributes_r8
     MODULE PROCEDURE NcDef_var_attributes_i_arr
     MODULE PROCEDURE NcDef_var_attributes_r4_arr
     MODULE PROCEDURE NcDef_var_attributes_r8_arr
  END INTERFACE NcDef_var_attributes
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
!  14 May 2014 - R. Yantosca - Add function NcBegin_Def to reopen define mode
!  14 May 2014 - R. Yantosca - Now use F90 free formatting
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
! !IROUTINE: NcDef_dimension
!
! !INTERFACE:
!
  SUBROUTINE NcDef_dimension(ncid,name,len,id,unlimited)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid  : netCDF file id
!!  name  : dimension name
!!  len   : dimension number
    CHARACTER (LEN=*), INTENT(IN)  :: name
    INTEGER,           INTENT(IN)  :: ncid, len
    LOGICAL, OPTIONAL, INTENT(IN)  :: unlimited
!
! !OUTPUT PARAMETERS:
!!  id    : dimension id
    INTEGER,           INTENT(OUT) :: id

    INTEGER  :: len0
!
! !DESCRIPTION: Defines dimension.
!\\
!\\
! !AUTHOR:
!  Jules Kouatchou and Maharaj Bhat
!
! !REVISION HISTORY:
!  Initial code.
!  18 May 2018 - C. Holmes - Add support for unlimited dimensions
!  25 Jun 2018 - R. Yantosca - Fixed typo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (len=512) :: err_msg
    INTEGER :: ierr

    ! If unlimited variable is present and true,
    ! then make this dimension unlimited
    len0 = len
    if (present(unlimited)) then
       if (unlimited) then
          len0 = NF_UNLIMITED
       endif
    endif

    ierr = Nf_Def_Dim (ncid, name, len0, id)

    IF (ierr.ne.NF_NOERR) then
       err_msg = 'Nf_Def_Dim: can not define dimension : '// Trim (name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_dimension
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_variable
!
! !INTERFACE:
!
  SUBROUTINE NcDef_variable(ncid,name,type,ndims,dims,var_id,compress)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!
!!  ncid   : netCDF file id
!!  name   : name of the variable
!!  type   : type of the variable
!!           (NF_FLOAT, NF_CHAR, NF_INT, NF_DOUBLE, NF_BYTE, NF_SHORT)
!!  ndims  : number of dimensions of the variable
!!  dims   : netCDF dimension id of the variable
    CHARACTER (LEN=*), INTENT(IN)  :: name
    INTEGER,           INTENT(IN)  :: ncid, ndims
    INTEGER,           INTENT(IN)  :: dims(ndims)
    INTEGER,           INTENT(IN)  :: type
    LOGICAL, OPTIONAL, INTENT(IN)  :: compress
!
! !OUTPUT PARAMETERS:
!
!!  varid  : netCDF variable id returned by NF_DEF_VAR
    INTEGER,           INTENT(OUT) :: var_id
!
! !DESCRIPTION: Defines a netCDF variable.
!\\
!\\
! !AUTHOR:
!  Jules Kouatchou and Maharaj Bhat
!
! !REVISION HISTORY:
!  Initial code.
!  17 Feb 2017 - C. Holmes   - Enable netCDF-4 compression
!  01 Mar 2017 - R. Yantosca - Add an #ifdef to enable netCDF4 compression
!                              only if the library has nf_def_var_deflate
!  10 May 2017 - R. Yantosca - Bug fix: var_id needs to be INTENT(OUT),
!                              because it's returned from NF_DEF_VAR
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer ::  ierr
    logical ::  doStop
    ! Compression settings
    ! choose deflate_level=1 for fast, minimal compression.
    ! Informal testing suggests minimal benefit from higher compression level
    integer, parameter :: shuffle=1, deflate=1, deflate_level=1
!
    ierr = Nf_Def_Var (ncid, name, type, ndims, dims, var_id)

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'Nf_Def_Var: can not define variable : '// Trim (name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

#if defined( NC_HAS_COMPRESSION )
    !=====================================================================
    ! If the optional "compress" variable is used and set to TRUE,
    ! then enable variable compression (cdh, 0/17/17)
    !
    ! NOTE: We need to block this out with an #ifdef because some
    ! netCDF installations might lack the nf_def_var_deflate function
    ! which would cause a compile-time error. (bmy, 3/1/17)
    !=====================================================================
    if (present(Compress)) then

       if (Compress) then

          ! Set compression
          ierr = nf_def_var_deflate( ncid, var_id,  shuffle,       &
                                           deflate, deflate_level )

          ! Check for errors.
          ! No message will be generated if the error is simply that the
          ! file is not netCDF-4
          ! (i.e. netCDF-3 don't support compression)
          IF ( (ierr.ne.NF_NOERR) .and. (ierr.ne.NF_ENOTNC4)) THEN

             ! Errors enabling compression will not halt the program
             doStop = .False.

             ! Print error
             err_msg = 'Nf_Def_Var_Deflate: can not compress variable : '// Trim (name)
             CALL Do_Err_Out (err_msg, doStop, 0, 0, 0, 0, 0.0d0, 0.0d0)
          END IF

       endif
    endif
#endif

  END SUBROUTINE NcDef_variable
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_c(ncid,var_id,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT none
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  var_id  : netCDF variable id
!!  att_name: attribute name
!!  att_val : attribute value
    CHARACTER (LEN=*), INTENT(IN) :: att_name, att_val
    INTEGER,           INTENT(IN) :: ncid, var_id
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
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: mylen, ierr
!
    mylen = LEN(att_val)
    ierr = Nf_Put_Att_Text (ncid, var_id, att_name, mylen, att_val)

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_var_attributes_c: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_c
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_i
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_i(ncid,var_id,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  var_id  : netCDF variable id
!!  att_name: attribute name
!!  att_val : attribute value
    INTEGER,           INTENT(IN) :: att_val
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid, var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: INTEGER.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!  12 Jun 2017 - R. Yantosca - Bug fix, should call NF_PUT_ATT_INT
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: mylen, ierr
!
    mylen = 1
    ierr  = Nf_Put_Att_Int( ncid,   var_id, att_name, &
                            NF_INT, mylen,  att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_var_attributes_i: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_i
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r4
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_r4(ncid,var_id,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  var_id  : netCDF variable id
!!  att_name: attribute name
!!  att_val : attribute value
    REAL*4,            INTENT(IN) :: att_val
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid, var_id
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
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: mylen, ierr
!
    mylen = 1
    ierr  = Nf_Put_Att_Real( ncid,     var_id, att_name, &
                             NF_FLOAT, mylen,  att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_var_attributes_r4: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_r4
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r8
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_r8(ncid,var_id,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT none
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  var_id  : netCDF variable id
!!  att_name: attribute name
!!  att_val : attribute value
    REAL*8,            INTENT(IN) :: att_val
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid, var_id
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
    CHARACTER (LEN=512) :: err_msg
    INTEGER             ::  mylen, ierr
!
    mylen = 1
    ierr  = Nf_Put_Att_Double( ncid,      var_id, att_name, &
                               NF_DOUBLE, mylen,  att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_var_attributes_r8: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_r8
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_i_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_i_arr(ncid,var_id,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT none
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  var_id  : netCDF variable id
!!  att_name: attribute name
!!  att_val : attribute value
    INTEGER,           INTENT(IN) :: att_val(:)
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid, var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: INTEGER vector.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!  12 Jun 2017 - R. Yantosca - Bug fix, should call NF_PUT_ATT_INT
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: mylen, ierr
!
    mylen = SIZE( att_val )
    ierr  = Nf_Put_Att_Int( ncid,   var_id, att_name, &
                            NF_INT, mylen,  att_val )

    iF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_var_attributes_i_arr: can not define attribute : ' &
            // TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_i_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r4_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_r4_arr(ncid,var_id,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT none
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  var_id  : netCDF variable id
!!  att_name: attribute name
!!  att_val : attribute value
    REAL*4,            INTENT(IN) :: att_val(:)
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid, var_id
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
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: mylen, ierr
!
    mylen = SIZE( att_val )
    ierr  = Nf_Put_Att_Real( ncid,     var_id, att_name, &
                             NF_FLOAT, mylen,  att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_var_attributes_r4_arr: can not define attribute : ' &
                    // TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_r4_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r8_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_r8_arr(ncid,var_id,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    var_id  : netCDF variable id
!!    att_name: attribute name
!!    att_val : attribute value
    REAL*8,            INTENT(IN) :: att_val(:)
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid, var_id
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
    character (len=512) :: err_msg
    integer             ::  mylen, ierr
!
    mylen = size( att_val )
    ierr  = Nf_Put_Att_Double( ncid,      var_id, att_name, &
                               NF_DOUBLE, mylen,  att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_var_attributes_r4_arr: can not define attribute : '&
                     // Trim (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_r8_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_c
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_c(ncid,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  att_name: attribute name
!!  att_val : attribute value
!
    CHARACTER (LEN=*), INTENT(IN) :: att_val
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid
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
    CHARACTER (LEN=512) :: err_msg
    INTEGER             ::  mylen, ierr
!
    mylen = len(att_val)
    ierr = Nf_Put_Att_Text (ncid, NF_GLOBAL, att_name, mylen, att_val)

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_c: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_c
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_i
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_i(ncid,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT none
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
    INTEGER,           INTENT(IN) :: att_val
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: NCID
!
! !DESCRIPTION: Defines global attributes of type: INTEGER
!\\
!\\
! !AUTHOR:
!  Bob Yantosca( based on code by Jules Kouatchou)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!  12 Jun 2017 - R. Yantosca - Bug fix, should call NF_PUT_ATT_INT
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: mylen, ierr
!
    mylen = 1
    ierr  = Nf_Put_Att_Int( ncid,   NF_GLOBAL, att_name, &
                            NF_INT, mylen,     att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_i: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_i
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r4
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_r4(ncid,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  att_name: attribute name
!!  att_val : attribute value
!
    REAL*4,            INTENT(IN) :: att_val
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid
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
    character (len=512) :: err_msg
    integer             :: mylen, ierr
!
    mylen = 1
    ierr  = Nf_Put_Att_Real( ncid,     NF_GLOBAL, att_name, &
                             NF_FLOAT, mylen,     att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_r4: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_r4
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r8
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_r8(ncid,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
    REAL*8,            INTENT(IN) :: att_val
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid
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
    character (len=512) :: err_msg
    integer             :: mylen, ierr
!
    mylen = 1
    ierr  = Nf_Put_Att_Double( ncid,     NF_GLOBAL, att_name, &
                               NF_FLOAT, mylen,     att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_r8: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_r8
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_i_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_i_arr(ncid,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
    INTEGER,           INTENT(IN) :: att_val(:)
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid
!
! !DESCRIPTION: Defines global attributes of type: INTEGER vector
!\\
!\\
! !AUTHOR:
!  Bob Yantosca( based on code by Jules Kouatchou)
!
! !REVISION HISTORY:
!  26 Sep 2013 - R. Yantosca - Initial version
!  12 Jun 2017 - R. Yantosca - Bug fix: Should call NF_PUT_ATT_INT
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: mylen, ierr
!
    mylen = SIZE( att_val )
    ierr  = Nf_Put_Att_Int( ncid,   NF_GLOBAL, att_name, &
                            NF_INT, mylen,     att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_i_arr: can not define attribute : ' &
            // Trim (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_i_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r4_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_r4_arr(ncid,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT none
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  att_name: attribute name
!!  att_val : attribute value
!
    REAL*4,            INTENT(IN) :: att_val(:)
    CHARACTER (LEN=*), INTENT(IN) :: att_name
    INTEGER,           INTENT(IN) :: ncid
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
    character (len=512) :: err_msg
    integer             :: mylen, ierr
!
    mylen = SIZE( att_val )
    ierr  = Nf_Put_Att_Real( ncid,     NF_GLOBAL, att_name, &
                             NF_FLOAT, mylen,     att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_r4_arr: can not define attribute : ' &
              // TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_r4_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r8_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_r8_arr(ncid,att_name,att_val)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT none
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  att_name: attribute name
!!  att_val : attribute value
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
    character (len=512) :: err_msg
    integer             :: mylen, ierr
!
    mylen = SIZE( att_val )
    ierr  = Nf_Put_Att_Double( ncid,     NF_GLOBAL, att_name, &
                               NF_FLOAT, mylen,     att_val )

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_r8_arr: can not define attribute : ' &
            // TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_r8_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcSetFill
!
! !INTERFACE:
!
  SUBROUTINE NcSetFill(ncid,ifill,omode)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(in) :: ncid, ifill,omode
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
    character (len=512) :: err_msg
    integer             ::  mylen, ierr
!
    ierr = Nf_Set_Fill (ncid, NF_NOFILL, omode)

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'Nf_Set_FIll: Error in omode  '
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcSetFill
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcEnd_Def
!
! !INTERFACE:
!
  SUBROUTINE NcEnd_Def(ncid)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT NONE
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: ncid
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
    CHARACTER (LEN=512) :: err_msg
    INTEGER             ::  ierr
!
    ierr = Nf_Enddef (ncid)

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'Nf_EndDef: Error in closing netCDF define mode!'
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcEnd_def
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcBegin_Def
!
! !INTERFACE:
!
  SUBROUTINE NcBegin_Def(ncid)
!
! !USES:
!
    USE m_do_err_out
!
    IMPLICIT none
!
    INCLUDE 'netcdf.inc'
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: ncid
!
! !DESCRIPTION: Opens (or re-opens) netCDF define mode, where variables
!  and attributes can be defined.
!\\
!\\
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  14 May 2014 - R. Yantosca - Initial version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: ierr
!
    ierr = Nf_Redef (ncid)

    IF (ierr.ne.NF_NOERR) THEN
       err_msg = 'Nf_ReDef: Error in opening netCDF define mode!'
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcBegin_Def
!EOC
END MODULE m_netcdf_io_define
