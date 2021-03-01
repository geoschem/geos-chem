!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_readattr.F90
!
! !INTERFACE:
!
MODULE m_netcdf_io_readattr
!
! !USES:

  USE m_do_err_out

  IMPLICIT NONE
  PRIVATE

  INCLUDE "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: NcGet_Var_Attributes
  INTERFACE NcGet_Var_Attributes
     MODULE PROCEDURE NcGet_Var_Attr_C
     MODULE PROCEDURE NcGet_Var_Attr_I4
     MODULE PROCEDURE NcGet_Var_Attr_R4
     MODULE PROCEDURE NcGet_Var_Attr_R8
     MODULE PROCEDURE NcGet_Var_Attr_I4_arr
     MODULE PROCEDURE NcGet_Var_Attr_R4_arr
     MODULE PROCEDURE NcGet_Var_Attr_R8_arr
  END INTERFACE

  PUBLIC :: NcGet_Glob_Attributes
  INTERFACE NcGet_Glob_Attributes
     MODULE PROCEDURE NcGet_Glob_Attr_C
     MODULE PROCEDURE NcGet_Glob_Attr_I4
     MODULE PROCEDURE NcGet_Glob_Attr_R4
     MODULE PROCEDURE NcGet_Glob_Attr_R8
     MODULE PROCEDURE NcGet_Glob_Attr_I4_arr
     MODULE PROCEDURE NcGet_Glob_Attr_R4_arr
     MODULE PROCEDURE NcGet_Glob_Attr_R8_arr
  END INTERFACE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcGet_Var_Attr_C
  PRIVATE :: NcGet_Var_Attr_I4
  PRIVATE :: NcGet_Var_Attr_R4
  PRIVATE :: NcGet_Var_Attr_R8
  PRIVATE :: NcGet_Var_Attr_I4_arr
  PRIVATE :: NcGet_Var_Attr_R4_arr
  PRIVATE :: NcGet_Var_Attr_R8_arr
  PRIVATE :: NcGet_Glob_Attr_C
  PRIVATE :: NcGet_Glob_Attr_I4
  PRIVATE :: NcGet_Glob_Attr_R4
  PRIVATE :: NcGet_Glob_Attr_R8
  PRIVATE :: NcGet_Glob_Attr_I4_arr
  PRIVATE :: NcGet_Glob_Attr_R4_arr
  PRIVATE :: NcGet_Glob_Attr_R8_arr
!
! !DESCRIPTION: Provides netCDF utility routines to read both netCDF
!  variable attributes and global attributes.  Individual routines for
!  reading attributes of different types are overloaded with F90
!  interfaces.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  30 Apr 2012 - R. Yantosca - Modified for compatibility with netCDF-3
!  30 Apr 2012 - R. Yantosca - Added comments
!  26 Sep 2013 - R. Yantosca - Add routines for reading vector attributes
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
! !IROUTINE: NcGet_Var_Attr_C
!
! !DESCRIPTION: Returns a variable attribute of type CHARACTER.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Var_Attr_C( fid, varName, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId        ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: varName    ! netCDF variable name
    CHARACTER(LEN=*), INTENT(IN)  :: attName    ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(OUT) :: attValue   ! Attribute value
!
! !DESCRIPTION: Reads a variable attribute (CHARACTER type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_TEXT,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg
    INTEGER            :: status, vId

    ! Zero return value
    attValue = ''

    ! Check if VARNAME is a valid variable
    status = Nf_Inq_Varid ( fId, varName, vId )

    ! Exit w/ error message if VARNAME is not valid
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_C: ' // TRIM( varName )        // &
                 ', '                   // Nf_Strerror( status )
       CALL Do_Err_Out ( errMsg, .TRUE., 1, fId, 0, 0, 0.0d0, 0.0d0)
    ENDIF

    !  Get the attribute
    status = Nf_Get_Att_Text( fId, vId, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_C: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Var_Attr_C
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Var_Attr_I4
!
! !DESCRIPTION: Returns a variable attribute of type INTEGER*4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Var_Attr_I4( fid, varName, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId        ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: varName    ! netCDF variable name
    CHARACTER(LEN=*), INTENT(IN)  :: attName    ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: attValue   ! Attribute value
!
! !DESCRIPTION: Reads a variable attribute (INTEGER type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_INT,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg
    INTEGER            :: status, vId

    ! Zero return value
    attValue = 0

    ! Check if VARNAME is a valid variable
    status = Nf_Inq_Varid ( fId, varName, vId )

    ! Exit w/ error message if VARNAME is not valid
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_I4: ' // TRIM( varName )        // &
                 ', '                   // Nf_Strerror( status )
       CALL Do_Err_Out ( errMsg, .TRUE., 1, fId, 0, 0, 0.0d0, 0.0d0)
    ENDIF

    ! Get the attribute
    status = Nf_Get_Att_Int( fId, vId, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_I4: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Var_Attr_I4
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Var_Attr_R4
!
! !DESCRIPTION: Returns a variable attribute of type REAL*4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Var_Attr_R4( fid, varName, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId        ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: varName    ! netCDF variable name
    CHARACTER(LEN=*), INTENT(IN)  :: attName    ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    REAL*4,           INTENT(OUT) :: attValue   ! Attribute value
!
! !DESCRIPTION: Reads a variable attribute (REAL*4 type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_REAL,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg
    INTEGER            :: status, vId

    ! Zero return value
    attValue = 0e0

    ! Check if VARNAME is a valid variable
    status = Nf_Inq_Varid ( fId, varName, vId )

    ! Exit w/ error message if VARNAME is not valid
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_R4: ' // TRIM( varName )        // &
                 ', '                   // Nf_Strerror( status )
       CALL Do_Err_Out ( errMsg, .TRUE., 1, fId, 0, 0, 0.0d0, 0.0d0)
    ENDIF

    ! Get the attribute
    status = Nf_Get_Att_Real( fId, vId, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_R4: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Var_Attr_R4
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Var_Attr_R8
!
! !DESCRIPTION: Returns a variable attribute of type REAL*8.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Var_Attr_R8( fid, varName, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId        ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: varName    ! netCDF variable name
    CHARACTER(LEN=*), INTENT(IN)  :: attName    ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    REAL*8,           INTENT(OUT) :: attValue   ! Attribute value
!
! !DESCRIPTION: Reads a variable attribute (REAL*4 type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use internal function NF_GET_ATT_DOUBLE,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg
    INTEGER            :: status, vId

    ! Zero return value
    attValue = 0d0

    ! Check if VARNAME is a valid variable
    status = Nf_Inq_Varid ( fId, varName, vId )

    ! Exit w/ error message if VARNAME is not valid
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_R8: ' // TRIM( varName )        // &
                 ', '                   // Nf_Strerror( status )
       CALL Do_Err_Out ( errMsg, .TRUE., 1, fId, 0, 0, 0.0d0, 0.0d0)
    ENDIF

    ! Get the attribute
    status = Nf_Get_Att_Double( fId, vId, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_R8: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Var_Attr_R8
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Var_Attr_I4_arr
!
! !DESCRIPTION: Returns a vector variable attribute of type INTEGER*4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Var_Attr_I4_arr( fid, varName, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId          ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: varName      ! netCDF variable name
    CHARACTER(LEN=*), INTENT(IN)  :: attName      ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: attValue(:)  ! Attribute value
!
! !DESCRIPTION: Reads a variable attribute (INTEGER type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_INT,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg
    INTEGER            :: status, vId

    ! Zero return value
    attValue = 0

    ! Check if VARNAME is a valid variable
    status = Nf_Inq_Varid ( fId, varName, vId )

    ! Exit w/ error message if VARNAME is not valid
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_I4_arr: ' // TRIM( varName )        // &
                 ', '                        // Nf_Strerror( status )
       CALL Do_Err_Out ( errMsg, .TRUE., 1, fId, 0, 0, 0.0d0, 0.0d0)
    ENDIF

    ! Get the attribute
    status = Nf_Get_Att_Int( fId, vId, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_I4_arr: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Var_Attr_I4_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Var_Attr_R4_arr
!
! !DESCRIPTION: Returns a vector variable attribute of type REAL*4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Var_Attr_R4_arr( fid, varName, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId          ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: varName      ! netCDF variable name
    CHARACTER(LEN=*), INTENT(IN)  :: attName      ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    REAL*4,           INTENT(OUT) :: attValue(:)  ! Attribute value
!
! !DESCRIPTION: Reads a variable attribute (REAL*4 type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_REAL,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg
    INTEGER            :: status, vId

    ! Zero return value
    attValue = 0e0

    ! Check if VARNAME is a valid variable
    status = Nf_Inq_Varid ( fId, varName, vId )

    ! Exit w/ error message if VARNAME is not valid
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_R4_arr: ' // TRIM( varName )        // &
                 ', '                        // Nf_Strerror( status )
       CALL Do_Err_Out ( errMsg, .TRUE., 1, fId, 0, 0, 0.0d0, 0.0d0)
    ENDIF

    ! Get the attribute
    status = Nf_Get_Att_Real( fId, vId, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_R4_arr: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Var_Attr_R4_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Var_Attr_R8_arr
!
! !DESCRIPTION: Returns a vector variable attribute of type REAL*8.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Var_Attr_R8_arr( fid, varName, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId          ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: varName      ! netCDF variable name
    CHARACTER(LEN=*), INTENT(IN)  :: attName      ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    REAL*8,           INTENT(OUT) :: attValue(:)  ! Attribute value
!
! !DESCRIPTION: Reads a variable attribute (REAL*4 type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use internal function NF_GET_ATT_DOUBLE,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg
    INTEGER            :: status, vId

    ! Zero return value
    attValue = 0d0

    ! Check if VARNAME is a valid variable
    status = Nf_Inq_Varid ( fId, varName, vId )

    ! Exit w/ error message if VARNAME is not valid
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_R8_arr: ' // TRIM( varName )        // &
                 ', '                        // Nf_Strerror( status )
       CALL Do_Err_Out ( errMsg, .TRUE., 1, fId, 0, 0, 0.0d0, 0.0d0)
    ENDIF

    ! Get the attribute
    status = Nf_Get_Att_Double( fId, vId, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Var_Attr_R8_arr: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Var_Attr_R8_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Glob_Attr_C
!
! !DESCRIPTION: Returns a variable attribute of type CHARACTER.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Glob_Attr_C( fid, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId        ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: attName    ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(OUT) :: attValue   ! Attribute value
!
! !DESCRIPTION: Reads a global attribute (CHARACTER type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_TEXT,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg, varName
    INTEGER            :: status

    ! Zero return value
    attValue = ''

    ! Get the attribute
    status = Nf_Get_Att_Text( fId, NF_GLOBAL, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Glob_Attr_C: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Glob_Attr_C
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Glob_Attr_I4
!
! !DESCRIPTION: Returns a variable attribute of type INTEGER*4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Glob_Attr_I4( fid, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId        ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: attName    ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: attValue   ! Attribute value
!
! !DESCRIPTION: Reads a global attribute (INTEGER type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_INT,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg, varName
    INTEGER            :: status

    ! Zero return value
    attValue = 0

    ! Get the attribute
    status = Nf_Get_Att_Int( fId, NF_GLOBAL, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Glob_Attr_I4: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Glob_Attr_I4
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Glob_Attr_R4
!
! !DESCRIPTION: Returns a variable attribute of type REAL*4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Glob_Attr_R4( fid, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId        ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: attName    ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    REAL*4,           INTENT(OUT) :: attValue   ! Attribute value
!
! !DESCRIPTION: Reads a global attribute (REAL*4 type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_REAL,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg, varName
    INTEGER            :: status

    ! Zero return value
    attValue = 0e0

    ! Get the attribute
    status = Nf_Get_Att_Real( fId, NF_GLOBAL, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Glob_Attr_R4: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Glob_Attr_R4
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Glob_Attr_R8
!
! !DESCRIPTION: Returns a variable attribute of type REAL*8.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Glob_Attr_R8( fid, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId        ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: attName    ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    REAL*8,           INTENT(OUT) :: attValue   ! Attribute value
!
! !DESCRIPTION: Reads a global attribute (REAL*8 type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_DOUBLE,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg, varName
    INTEGER            :: status

    ! Zero return value
    attValue = 0d0

    ! Get the attribute
    status = Nf_Get_Att_Double( fId, NF_GLOBAL, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Glob_Attr_R8: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Glob_Attr_R8
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Glob_Attr_I4_arr
!
! !DESCRIPTION: Returns a variable attribute of type INTEGER*4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Glob_Attr_I4_arr( fid, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId          ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: attName      ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: attValue(:)  ! Attribute value
!
! !DESCRIPTION: Reads a global attribute (INTEGER type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_INT,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg, varName
    INTEGER            :: status

    ! Zero return value
    attValue = 0

    ! Get the attribute
    status = Nf_Get_Att_Int( fId, NF_GLOBAL, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Glob_Attr_I4_arr: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Glob_Attr_I4_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Glob_Attr_R4_arr
!
! !DESCRIPTION: Returns a variable attribute of type REAL*4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Glob_Attr_R4_arr( fid, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId          ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: attName      ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    REAL*4,           INTENT(OUT) :: attValue(:)  ! Attribute value
!
! !DESCRIPTION: Reads a global attribute (REAL*4 type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_REAL,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg, varName
    INTEGER            :: status

    ! Zero return value
    attValue = 0e0

    ! Get the attribute
    status = Nf_Get_Att_Real( fId, NF_GLOBAL, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Glob_Attr_R4_arr: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Glob_Attr_R4_arr
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcGet_Glob_Attr_R8
!
! !DESCRIPTION: Returns a variable attribute of type REAL*8.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcGet_Glob_Attr_R8_arr( fid, attName, attValue )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId          ! netCDF file ID
    CHARACTER(LEN=*), INTENT(IN)  :: attName      ! Name of variable attribute
!
! !OUTPUT PARAMETERS:
!
    REAL*8,           INTENT(OUT) :: attValue(:)  ! Attribute value
!
! !DESCRIPTION: Reads a global attribute (REAL*8 type) from a netCDF file.
!\\
!\\
! !AUTHOR:
!  Bob Yantosca (based on code by Jules Kouatchou and Maharaj Bhat)
!
! !REVISION HISTORY:
!  25 Jan 2012 - R. Yantosca - Initial version
!  31 Jan 2012 - R. Yantosca - Zero attValue before reading attributes
!  30 Apr 2012 - R. Yantosca - Use netCDF library function NF_GET_ATT_DOUBLE,
!                              which is compatible w/ netCDF3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: errMsg, varName
    INTEGER            :: status

    ! Zero return value
    attValue = 0d0

    ! Get the attribute
    status = Nf_Get_Att_Double( fId, NF_GLOBAL, attName, attValue )

    ! Exit w/ error message if unsuccessful
    IF ( status /= NF_NOERR ) THEN
       errMsg = 'In NcGet_Glob_Attr_R8_arr: cannot read attribute : ' // &
                 TRIM( attName )
       CALL Do_Err_Out( errMsg, .TRUE., 0, 0, 0, 0, 0.0d0, 0.0d0 )
    endif

  END SUBROUTINE NcGet_Glob_Attr_R8_arr
!EOC
END MODULE m_netcdf_io_readattr
