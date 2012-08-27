!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: smv_errcode_mod.F
!
! !DESCRIPTION: This include file contains the various success or failure
!  parameters for the GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE: 
!
      MODULE SMV_ERRCODE_MOD
!
! !USES:
!
      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS: 
!
      ! Return w/ success
      INTEGER, PARAMETER :: SMV_SUCCESS =  0

      ! Return w/ failure
      INTEGER, PARAMETER :: SMV_FAILURE = -1

!------------------------------------------------------------------------------
! Added from BMY's column module code
      INTEGER, PARAMETER :: SMV_FAIL_INTERFACE   = -2

      ! Return codes for emissions readers & internal state routines
      INTEGER, PARAMETER :: SMV_FAIL_EM_RD_INIT  = -1000
      INTEGER, PARAMETER :: SMV_FAIL_EM_RD_RUN   = -1001
      INTEGER, PARAMETER :: SMV_FAIL_EM_RD_FINAL = -1002
      INTEGER, PARAMETER :: SMV_FAIL_EM_INIT     = -1003
      INTEGER, PARAMETER :: SMV_FAIL_EM_RUN      = -1004
      INTEGER, PARAMETER :: SMV_FAIL_EM_FINAL    = -1005

      ! Return codes for initialization routines
      INTEGER, PARAMETER :: SMV_FAIL_INIT_COL    = -1100
      INTEGER, PARAMETER :: SMV_FAIL_INIT_SCHEM  = -1101
      INTEGER, PARAMETER :: SMV_FAIL_SCHEM_READ  = -1102
      INTEGER, PARAMETER :: SMV_FAIL_READER      = -1103
      INTEGER, PARAMETER :: SMV_FAIL_READCHEM    = -1104
      INTEGER, PARAMETER :: SMV_FAIL_JSPARSE     = -1105
      INTEGER, PARAMETER :: SMV_FAIL_KSPARSE     = -1106
      INTEGER, PARAMETER :: SMV_FAIL_ILISOPOH    = -1107
      INTEGER, PARAMETER :: SMV_FAIL_SETTRACE    = -1108
      INTEGER, PARAMETER :: SMV_FAIL_SETEMDEP    = -1109
      INTEGER, PARAMETER :: SMV_FAIL_COLUMN      = -1150
      
      ! Return codes for chemistry setup routines
      INTEGER, PARAMETER :: SMV_FAIL_GASCONC     = -1200
      INTEGER, PARAMETER :: SMV_FAIL_PARTITION   = -1201
      INTEGER, PARAMETER :: SMV_FAIL_RDUST       = -1202
      INTEGER, PARAMETER :: SMV_FAIL_RDAER       = -1203
      INTEGER, PARAMETER :: SMV_FAIL_SETEMIS     = -1204

      ! Return codes for FAST-J routines
      INTEGER, PARAMETER :: SMV_FAIL_INPHOT      = -1300
      INTEGER, PARAMETER :: SMV_FAIL_FAST_J      = -1301
      INTEGER, PARAMETER :: SMV_FAIL_PHOTOJ      = -1302
      INTEGER, PARAMETER :: SMV_FAIL_SET_PROF    = -1303
      INTEGER, PARAMETER :: SMV_FAIL_JVALUE      = -1304
      INTEGER, PARAMETER :: SMV_FAIL_JRATET      = -1305 
      INTEGER, PARAMETER :: SMV_FAIL_OPMIE       = -1306

      ! Return codes for chemistry routines
      INTEGER, PARAMETER :: SMV_FAIL_PHYSPROC    = -1400
      INTEGER, PARAMETER :: SMV_FAIL_CALCRATE    = -1401
      INTEGER, PARAMETER :: SMV_FAIL_SMVGEAR     = -1402
      INTEGER, PARAMETER :: SMV_FAIL_SUBFUN      = -1403
      INTEGER, PARAMETER :: SMV_FAIL_PDERIV      = -1404
      INTEGER, PARAMETER :: SMV_FAIL_DECOMP      = -1405
      INTEGER, PARAMETER :: SMV_FAIL_BACKSUB     = -1406
      INTEGER, PARAMETER :: SMV_FAIL_LUMP        = -1407
      INTEGER, PARAMETER :: SMV_FAIL_SCHEM       = -1420
      INTEGER, PARAMETER :: SMV_FAIL_CHEMDR      = -1421

      ! Return codes for aerosol chemistry
      INTEGER, PARAMETER :: SMV_FAIL_RPMARES     = -1500
      INTEGER, PARAMETER :: SMV_FAIL_CHEMSEASALT = -1500
      INTEGER, PARAMETER :: SMV_FAIL_DUST        = -1501
      INTEGER, PARAMETER :: SMV_FAIL_GINOUX      = -1502 
      INTEGER, PARAMETER :: SMV_FAIL_EMDUSTBOX   = -1503

      ! Return codes for deposition
      INTEGER, PARAMETER :: SMV_FAIL_DRYDEP      = -1600
      INTEGER, PARAMETER :: SMV_FAIL_DEPVEL      = -1602
      INTEGER, PARAMETER :: SMV_FAIL_DRYFLX      = -1603

      ! Return codes for PBL mixing
      INTEGER, PARAMETER :: SMV_FAIL_COMPUTE_PBL = -1700
      INTEGER, PARAMETER :: SMV_FAIL_DO_PBL_MIX  = -1701

      ! Return codes for wet deposition
      INTEGER, PARAMETER :: SMV_FAIL_INIT_WD     = -1800
      INTEGER, PARAMETER :: SMV_FAIL_INIT_WDF    = -1801
      INTEGER, PARAMETER :: SMV_FAIL_MAKE_QQ     = -1802
      INTEGER, PARAMETER :: SMV_FAIL_WETDEP      = -1803

      ! Return codes for soil NOx emissions
      INTEGER, PARAMETER :: SMV_FAIL_CANOPYNOX   = -1900
      INTEGER, PARAMETER :: SMV_FAIL_SOILNOX     = -1902

!------------------------------------------------------------------------------

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
!  03 Aug 2011 - M. Long     - Converted from Header file to Module
!EOP
!------------------------------------------------------------------------------
!BOC
	END MODULE SMV_ERRCODE_MOD
!EOC
