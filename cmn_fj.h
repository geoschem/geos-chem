! $Id: cmn_fj.h,v 1.4 2009/05/06 14:14:47 ccarouge Exp $
!
!******************************************************************************
!  CMN_FJ.H -- Header file containing parameters and common
!  blocks used to interface between Harvard chemistry and UC-Irvine 
!  Fast-J photolysis programs.
!  
!  Based on code from Oliver Wild (9 Jul 1999)
!
!  NOTES:
!  (1 ) Uses Fortran 90 declarations for parameters and variables
!  (2 ) Pass CTM size parameters and preprocessor switches via CMN_SIZE.
!  (3 ) Update JPMAX for new chemistry mechanism (amf, bmy, 4/20/00)
!  (4 ) Return JPMAX to original setting (bmy, 9/25/00)
!  (5 ) Return JPMAX to 55 for peroxy recycling (again) (bmy, 12/20/00)
!  (6 ) Now need to use the window parameters IIPAR,JJPAR,LLPAR (bmy, 9/25/01)
!  (7 ) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!        compilation. (bmy, 6/25/02)
!  (8 ) Replaced ESIG array with ETAA and ETAB arrays for the hybrid
!        pressure formulation.  Also deleted PREST, since we don't need that
!        anymore. (bmy, 8/23/02)
!  (9 ) Removed ETAA and ETAB arrays.  We now compute PJ directly from the 
!        GET_PEDGE routine. (bmy, 10/30/07)
!  (10) Increase photolysis rxns JPMAX = 79 (tmf, 1/7/09)
!
!          - Bob Yantosca [bmy@io.as.harvard.edu], 30 Oct 2007
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Global array sizes in longitude, latitude, altitude
      INTEGER, PARAMETER :: IPAR = IIPAR
      INTEGER, PARAMETER :: JPAR = JJPAR
      INTEGER, PARAMETER :: LPAR = LLPAR

      ! max # of photolysis rxns = 4 + IPHOT (see comode.h)
      INTEGER, PARAMETER :: JPMAX = 79

      ! Variables for number of layers and number of photolysis rxns
      INTEGER            :: JPNL, JPPJ       
      COMMON /FJ_INTEG/     JPNL, JPPJ

      ! Branches for photolysis species
      INTEGER            :: BRANCH        
      COMMON /FJ_BRANCH/    BRANCH(JPMAX)

      ! Names of photolysis species
      CHARACTER (LEN=4)  :: RNAMES
      COMMON /FJ_NAME/      RNAMES(JPMAX)

      ! Mapping array from Harvard species names to UCI species names
      INTEGER            :: RINDEX                
      COMMON /FJ_INDX/      RINDEX(JPMAX)

      ! Output J-values
      REAL*8             :: ZPJ        
      COMMON /FJ_VALUE/     ZPJ(LPAR,JPMAX,IPAR,JPAR)
