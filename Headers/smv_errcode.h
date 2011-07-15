!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: smv_errcode.h
!
! !DESCRIPTION: This include file contains the various success or failure
!  parameters for the GEOS-Chem column chemistry code.
!\\
!\\
! !DEFINED PARAMETERS: 
!
      ! Return w/ success
      INTEGER, PARAMETER :: SMV_SUCCESS          =  0

      ! Return w/ failure
      INTEGER, PARAMETER :: SMV_FAILURE          = -1

      !----------------------------------------------------------------------
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% NOTE: The following are deprecated and will be removed later %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !
      !! Return codes for emissions readers & internal state routines
      !INTEGER, PARAMETER :: SMV_FAIL_EM_RD_INIT  = -1000
      !INTEGER, PARAMETER :: SMV_FAIL_EM_RD_RUN   = -1001
      !INTEGER, PARAMETER :: SMV_FAIL_EM_RD_FINAL = -1002
      !INTEGER, PARAMETER :: SMV_FAIL_EM_INIT     = -1003
      !INTEGER, PARAMETER :: SMV_FAIL_EM_RUN      = -1004
      !INTEGER, PARAMETER :: SMV_FAIL_EM_FINAL    = -1005
      !
      !! Return codes for emissions
      INTEGER, PARAMETER :: SMV_FAIL_GINOUX      = -1502 
      INTEGER, PARAMETER :: SMV_FAIL_EMDUSTBOX   = -1503
      !INTEGER, PARAMETER :: SMV_FAIL_CANOPYNOX   = -1900
      !INTEGER, PARAMETER :: SMV_FAIL_SOILNOX     = -1902
      !----------------------------------------------------------------------
!
! !REVISION HISTORY: 
!  20 Mar 2009 - R. Yantosca - Initial version
!  15 Jul 2009 - R. Yantosca - Updated w/ error codes for drydep,
!                              wetdep, and PBL mixing routines
!  03 Nov 2009 - R. Yantosca - Added error codes for column & interface
!  14 Dec 2009 - R. Yantosca - Added error code for unit conversion
!  01 Feb 2010 - R. Yantosca - Added error code for ISORROPIA ATE code
!  06 May 2010 - R. Yantosca - Deleted redundant error codes
!  03 Jun 2010 - R. Yantosca - Deleted error codes for SCHEM routines
!EOP
!------------------------------------------------------------------------------
