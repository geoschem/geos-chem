!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  PRIVATE
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
!  See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: NcDef_dimension
!
! !INTERFACE:
!
  SUBROUTINE NcDef_dimension(ncid,name,len,id,unlimited)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
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
          len0 = NF90_UNLIMITED
       endif
    endif

    ierr = NF90_Def_Dim(ncid, name, len0, id)

    IF (ierr.ne.NF90_NOERR) then
       err_msg = 'NF90_Def_Dim: can not define dimension : '// Trim (name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_dimension
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_variable
!
! !INTERFACE:
!
  SUBROUTINE NcDef_variable(ncid, name, xtype, ndims, dims, var_id, compress)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
!
! !INPUT PARAMETERS:
!
!!  ncid   : netCDF file id
!!  name   : name of the variable
!!  type   : type of the variable
!!           (NF90_FLOAT,  NF90_CHAR, NF90_INT,
!!            NF90_DOUBLE, NF90_BYTE, NF90_SHORT)
!!  ndims  : number of dimensions of the variable
!!  dims   : netCDF dimension id of the variable
    CHARACTER (LEN=*), INTENT(IN)  :: name
    INTEGER,           INTENT(IN)  :: ncid, ndims
    INTEGER,           INTENT(IN)  :: dims(ndims)
    INTEGER,           INTENT(IN)  :: xtype
    LOGICAL, OPTIONAL, INTENT(IN)  :: compress
!
! !OUTPUT PARAMETERS:
!
!!  varid  : netCDF variable id returned by NF90_DEF_VAR
    INTEGER,           INTENT(OUT) :: var_id
!
! !DESCRIPTION: Defines a netCDF variable.
!\\
!\\
! !AUTHOR:
!  Jules Kouatchou and Maharaj Bhat
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character(len=512) :: err_msg
    integer            :: ierr
    logical            :: doStop

#ifdef NC_HAS_COMPRESSION
    !=====================================================================
    ! Create a compressed (deflated) netCDF variable
    !
    ! NOTE: We need to block this out with an #ifdef because some
    ! netCDF installations might lack the deflation capability,
    ! which would cause a compile-time error. (bmy, 3/1/17)
    !========================================================================
    IF ( PRESENT( Compress ) ) then

       !------------------------------------------------------------------
       ! If COMPRESS is passed as an optional argument, and is TRUE,
       ! then define the variable with deflate_level=1.  Higher values
       ! of deflate_level yield minimal additiional benefit.
       !
       ! ALSO NOTE: Newer versions of netCDF balk when you try to compress
       ! a scalar variable.  This generates an annoying warning message.
       ! To avoid this, only compress array variables. (bmy, 11/30/20)
       !-------------------------------------------------------------------
       IF ( Compress .and. ndims > 0 ) THEN

          ! Create deflated variable
          ierr = NF90_Def_Var( ncid, name, xtype, dims, var_id,              &
                               shuffle=.TRUE., deflate_level=1 )

          ! Check for errors.
          ! No message will be generated if the error is simply that the
          ! file is not netCDF-4 (as netCDF-3 doesn't support compression)
          IF ( (ierr.ne.NF90_NOERR) .and. (ierr.ne.NF90_ENOTNC4)) THEN

             ! Errors enabling compression will not halt the program
             doStop = .False.

             ! Print error
             err_msg = 'NF90_Def_Var: can not create compressed variable : '//&
                        Trim(name)
             CALL Do_Err_Out (err_msg, doStop, 0, 0, 0, 0, 0.0d0, 0.0d0)
          END IF

          ! Return successfully
          RETURN
       ENDIF
    ENDIF
#endif

    !========================================================================
    ! Create an uncompressed netCDF variable if:
    ! (1) COMPRESS is not passed as an optional argument
    ! (2) COMPRESS is passed as an optional argument but is FALSE
    ! (3) The variable is a scalar (ndims == 0)
    !========================================================================
    ierr = NF90_Def_Var( ncid, name, xtype, dims, var_id )
    IF ( ierr /= NF90_NOERR ) THEN
       err_msg = 'NF90_Def_Var_Deflate: can not create variable : '// &
            Trim (name)
       CALL Do_Err_Out (err_msg, doStop, 0, 0, 0, 0, 0.0d0, 0.0d0)
    ENDIF

  END SUBROUTINE NcDef_variable
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_c(ncid, var_id, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
!
! !INPUT PARAMETERS:
!!  ncid    : netCDF file id
!!  var_id  : netCDF variable id
!!  att_name: attribute name
!!  att_val : attribute value
    CHARACTER (LEN=*), INTENT(IN) :: att_name, att_val
    INTEGER,           INTENT(IN) :: ncid,     var_id
!
! !DESCRIPTION: Defines a netCDF variable attribute of type: CHARACTER.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: ierr
!
    ierr = NF90_Put_Att(ncid, var_id, att_name, att_val)

    IF (ierr /= NF90_NOERR) THEN
       err_msg = 'NcDef_var_attributes_c: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_c
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_i
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_i(ncid, var_id, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: mylen, ierr
!
    ierr  = NF90_Put_Att( ncid, var_id, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_var_attributes_i: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_i
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r4
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_r4(ncid, var_id, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: ierr
!
    ierr  = NF90_Put_Att( ncid, var_id, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_var_attributes_r4: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_r4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r8
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_r8(ncid, var_id, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: ierr
!
    ierr  = NF90_Put_Att( ncid, var_id, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_var_attributes_r8: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_r8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_i_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_i_arr(ncid, var_id, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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

!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: ierr
!
    ierr  = NF90_Put_Att( ncid, var_id, att_name, att_val )

    iF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_var_attributes_i_arr: can not define attribute : ' &
            // TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_i_arr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r4_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_r4_arr(ncid, var_id, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: ierr
!
    ierr  = NF90_Put_Att( ncid, var_id, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_var_attributes_r4_arr: can not define attribute : ' &
                    // TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_r4_arr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes_r8_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_var_attributes_r8_arr(ncid, var_id, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: ierr
!
    ierr  = NF90_Put_Att( ncid, var_id, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_var_attributes_r4_arr: can not define attribute : '&
                     // Trim (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_var_attributes_r8_arr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_c
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_c(ncid, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: ierr
!
    ierr = NF90_Put_Att(ncid, NF90_GLOBAL, att_name, att_val)

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_c: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_c
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_i
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_i(ncid, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: ierr
!
    ierr  = NF90_Put_Att( ncid, NF90_GLOBAL, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_i: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_i
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r4
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_r4(ncid, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: ierr
!
    ierr  = NF90_Put_Att( ncid, NF90_GLOBAL, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_r4: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_r4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r8
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_r8(ncid, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: ierr
!
    ierr  = NF90_Put_Att( ncid, NF90_GLOBAL, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_r8: can not define attribute : ' // &
            TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_r8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_i_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_i_arr(ncid, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: ierr
!
    ierr  = NF90_Put_Att( ncid, NF90_GLOBAL, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_i_arr: can not define attribute : ' &
            // Trim (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_i_arr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: ierr
!
    ierr  = NF90_Put_Att( ncid, NF90_GLOBAL, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_r4_arr: can not define attribute : ' &
              // TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_r4_arr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes_r8_arr
!
! !INTERFACE:
!
  SUBROUTINE NcDef_glob_attributes_r8_arr(ncid, att_name, att_val)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: ierr
!
    ierr  = NF90_Put_Att( ncid, NF90_GLOBAL, att_name, att_val )

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NcDef_glob_attributes_r8_arr: can not define attribute : ' &
            // TRIM (att_name)
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcDef_glob_attributes_r8_arr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcSetFill
!
! !INTERFACE:
!
  SUBROUTINE NcSetFill(ncid, ifill, omode)
!
! !USES:
!
    USE netCDF
    USE m_do_err_out
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN   ) :: ncid, ifill
    INTEGER, INTENT(INOUT) :: omode
!
! !DESCRIPTION: Sets fill method.
!\\
!\\
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: ierr
!
    ierr = NF90_Set_Fill(ncid, ifill, omode)

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NF90_Set_FIll: Error in omode  '
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcSetFill
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    CHARACTER (LEN=512) :: err_msg
    INTEGER             :: ierr
!
    ierr = NF90_Enddef(ncid)

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NF90_EndDef: Error in closing netCDF define mode!'
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcEnd_def
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    USE netCDF
    USE m_do_err_out
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=512) :: err_msg
    integer             :: ierr
!
    ierr = NF90_Redef (ncid)

    IF (ierr.ne.NF90_NOERR) THEN
       err_msg = 'NF90_ReDef: Error in opening netCDF define mode!'
       CALL Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
    END IF

  END SUBROUTINE NcBegin_Def
!EOC
END MODULE m_netcdf_io_define
