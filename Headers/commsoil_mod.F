!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: commsoil_mod
!
! !DESCRIPTION: Module commsoil\_mod contains global variables for the 
!  soil NOx emissions routines.
!\\
!\\
! !INTERFACE: 
!
      MODULE COMMSOIL_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, MAXIJ, NTYPE, NVEGTYPE
      USE CMN_SIZE_MOD, ONLY : IGLOB, JGLOB

      IMPLICIT NONE
      PUBLIC

#     include "define.h"
!
! !DEFINED PARAMETERS: 
!      
      ! The defined soil types (Olson soil types)
      INTEGER, PARAMETER :: NSOIL = 11

      ! Number of MODIS/Koppen biome types
      INTEGER, PARAMETER :: NSOILB = 24

#if   defined( GRID4x5  )

      ! There are 1118 land boxes for the 4 x 5 GLOBAL GRID
      INTEGER, PARAMETER :: NLAND = 1118

#elif defined( GRID2x25 )

      ! There are 3920 land boxes for the 2 x 2.5 GLOBAL GRID
      INTEGER, PARAMETER :: NLAND = 3920 

#elif defined( GRID1x125 )

      !%%% NOTE: still to be determined
      INTEGER, PARAMETER :: NLAND = 9999 

#elif defined( GRID1x1 ) && defined( NESTED_CH )

      ! There are 2861 land points for the 1x1 CHINA nested grid
      INTEGER, PARAMETER :: NLAND = 2861

#elif defined( GRID1x1 ) && defined( NESTED_NA )

      ! There are 2118  land points for the 1x1 N. AMERICA nested grid
      INTEGER, PARAMETER :: NLAND = 2118

#elif defined( GRID1x1 )

      ! There are 17174 land points for the 1x1 GLOBAL grid
      INTEGER, PARAMETER :: NLAND=17174

#elif defined( GRID05x0666 ) && defined( NESTED_CH )

      ! There are 8261 land points for the 0.5 x 0.666 CHINA nested grid
      INTEGER, PARAMETER :: NLAND = 8261
      
#elif defined( GRID05x0666 ) && defined( NESTED_NA )

      !%%% NOTE: still to be determined
      INTEGER, PARAMETER :: NLAND = 8568   

#elif defined( GRID05x0666 ) && defined( NESTED_EU )

      !%%% NOTE: still to be determined
      INTEGER, PARAMETER :: NLAND = 5536

#elif defined( GRID025x03125 ) && defined( NESTED_CH )

      !%%% NOTE: still to be determined, use fudge factor for now
      INTEGER, PARAMETER :: NLAND = 5536

#elif defined( GRID025x03125 ) && defined( SEAC4RS )

      !%%% NOTE: still to be determined, use fudge factor for now
      INTEGER, PARAMETER :: NLAND = 5536

#elif defined( EXTERNAL_GRID ) && defined( DEVEL )
      !%%% NOTE: still to be determined!!!
      INTEGER, PARAMETER :: NLAND = 5536

#endif
!
! !PUBLIC DATA MEMBERS:
!
      ! Land types, Yienger & Levy (Used by the drydep code):
      ! water/desert/ice//Trop. Rain. Forst.//conifers//dry deciduous//
      ! other deciduous//woodland//grassland//agriculture (other than rice)
      ! rice paddies//wetland/tundra

      ! i,j of the grid
      INTEGER :: INDEXSOIL(2,NLAND) 

      ! New Soil NOx model arrays
      REAL*8 :: SOILFERT(IGLOB,JGLOB,366)
      REAL*4 :: CLIM(IGLOB,JGLOB,2)      !Fraction of arid (layer 1) and non-arid (layer 2) land
      REAL*4 :: LAND2(IGLOB,JGLOB,24)    ! MODIS LANDTYPE
      REAL*4 :: DRYPERIOD(IGLOB,JGLOB)   ! Dry period length
      REAL*4 :: PFACTOR(IGLOB,JGLOB)     ! Pulse Factor
      REAL*4 :: GWET_PREV(IGLOB,JGLOB)   ! Pulse Factor
      REAL*8 :: SOILEXC(NSOILB)

      ! Stores output
#if defined( DEVEL )
      REAL*8, ALLOCATABLE  :: SOILNOX(:,:)
#else
      REAL*8  :: SOILNOX(IIPAR,JJPAR) 
#endif

      !Soil NOx Diag49
      REAL*8 INST_SOIL(IGLOB,JGLOB) !instantaneous soil nox
      REAL*8 INST_FERT(IGLOB,JGLOB) !instantaneous fert nox

      ! Olson->soil type,nvegtype in commbio.h
      INTEGER :: NCONSOIL(NVEGTYPE) 

      ! Track NOx within canopy dry dep.
#if defined( DEVEL )
      REAL*8, ALLOCATABLE  :: CANOPYNOX(:,:)
#else
      REAL*8  :: CANOPYNOX(MAXIJ,NSOILB) 
#endif

      !Soil NOx deposited N arrays
      REAL*8 DRY_NO2(IGLOB,JGLOB)
      REAL*8 DRY_PAN(IGLOB,JGLOB)
      REAL*8 DRY_HNO3(IGLOB,JGLOB)
      REAL*8 DRY_NH3(IGLOB,JGLOB)
      REAL*8 DRY_NH4(IGLOB,JGLOB)
      REAL*8 DRY_NIT(IGLOB,JGLOB)
      REAL*8 DEP_RESERVOIR(IGLOB,JGLOB)
      REAL*8 WET_HNO3(IGLOB,JGLOB)
      REAL*8 WET_NH3(IGLOB,JGLOB)
      REAL*8 WET_NH4(IGLOB,JGLOB)
      REAL*8 WET_NIT(IGLOB,JGLOB)

      !MODIS/Koppen resistance values
      INTEGER  :: SNIMODIS(NSOILB)
      INTEGER  :: SNIRI  (NSOILB)       
      INTEGER  :: SNIRLU (NSOILB)      
      INTEGER  :: SNIRAC (NSOILB)      
      INTEGER  :: SNIRGSS(NSOILB)     
      INTEGER  :: SNIRGSO(NSOILB)     
      INTEGER  :: SNIRCLS(NSOILB)     
      INTEGER  :: SNIRCLO(NSOILB) 
      INTEGER  :: SNIVSMAX(NSOILB) 

!
! !REMARKS:
!  Original code from:
!    HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP                       
!    MODULE FOR SOIL NOx EMISSIONS                                      
!    by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob               
!    (Release V2.1) 
!                                                                             .
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS DEPRECATED!  IT WILL BE REMOVED WHEN WE  %%%
!  %%%  REPLACE THE SOIL NOX EMISSIONS ALGORITHM (IN THE NEAR FUTURE  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! !
! !REVISION HISTORY:
!  (1 ) Be sure to force double precision with the DBLE function            
!        and the "D" exponent, wherever necessary (bmy, 10/6/99)             
!  (2 ) Changed RCS ID tag comment character from "C" to "!" to allow 
!        freeform compilation.  Also added & continuation characters in 
!        column 73 to allow header files to be included in F90 freeform 
!        files. Updated comments, cosmetic changes. (bmy, 6/25/02)
!  (3 ) Now use cpp switches to define 1x1 parameters.  Also added
!        space in the #ifdef block for the 1x125 grid (bmy, 12/1/04)
!  (4 ) Bug fix: 2681 should be 2861 in NLAND (bmy, 9/22/06)
!  (5 ) Set # of land boxes for GEOS-5 nested grids (yxw, dan, bmy, 11/6/08)
!  (6 ) Set # of land boxes for GEOS-5 EUROPE nested grid (amv, 10/19/09)
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!  30 Aug 2012 - J.D. Maasakkers - Removed all obsolete old soil NOx code data
!EOP
!------------------------------------------------------------------------------
!BOC
#if defined( DEVEL )
      CONTAINS

      SUBROUTINE SET_COMMSOIL_MOD

      IMPLICIT NONE

      INTEGER AS

      ALLOCATE(
     &     SOILNOX(IIPAR,JJPAR),
     &     CANOPYNOX(MAXIJ,NTYPE), 
     &     STAT=AS)

      END SUBROUTINE SET_COMMSOIL_MOD
#endif
      END MODULE COMMSOIL_MOD
!EOC
