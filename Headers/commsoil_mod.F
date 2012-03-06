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

      IMPLICIT NONE
      PUBLIC

#     include "define.h"
!
! !DEFINED PARAMETERS: 
!      
      ! The defined soil types
      INTEGER, PARAMETER :: NSOIL = 11

      ! Number of soil pulsing types
      INTEGER, PARAMETER :: NPULSE = 3  

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

#endif
!
! !PUBLIC DATA MEMBERS:
!
      ! Land types:
      ! water/desert/ice//Trop. Rain. Forst.//conifers//dry deciduous//
      ! other deciduous//woodland//grassland//agriculture (other than rice)
      ! rice paddies//wetland/tundra

      ! i,j of the grid
      INTEGER :: INDEXSOIL(2,NLAND) 

      ! Tracking of wet/dry & three types of pulsing (Y&L, 94)
      REAL*8  :: SOILPULS(NPULSE+1,NLAND) 

      ! Two month observed precip
      REAL*8  :: SOILPREP(2,NLAND)  

      ! Fertilizers
      REAL*8  :: SOILFERT(NLAND)    

      ! Pulsing factors
      REAL*8  :: PULSFACT(NPULSE)

      ! Pulsing decay per timestep
      REAL*8  :: PULSDECAY(NPULSE)

      ! Stores output
#if defined( DEVEL )
      REAL*8, ALLOCATABLE  :: SOILNOX(:,:)
#else
      REAL*8  :: SOILNOX(IIPAR,JJPAR) 
#endif

      ! Olson->soil type,nvegtype in commbio.h
      INTEGER :: NCONSOIL(NVEGTYPE) 

      ! Track NOx within canopy dry dep.
#if defined( DEVEL )
      REAL*8, ALLOCATABLE  :: CANOPYNOX(:,:)
#else
      REAL*8  :: CANOPYNOX(MAXIJ,NTYPE) 
#endif

      ! Canopy wind extinction coeff from Y&L
      REAL*8  :: SOILTA(NSOIL)
      REAL*8  :: SOILTB(NSOIL)
      REAL*8  :: SOILAW(NSOIL)
      REAL*8  :: SOILAD(NSOIL)
      REAL*8  :: SOILEXC(NSOIL)     

      ! The correct sequence of PULSFACT is 5, 10, 15 
      DATA PULSFACT  / 5.D0,    10.D0,   15.D0   /

      ! PULSDECAY now contains the correct decay factors from Yienger & Levy
      DATA PULSDECAY / 0.805D0, 0.384D0, 0.208D0 / 

      ! SOILTA = Coefficient used to convert from surface temperture to  
      !          soil temperature     
      DATA SOILTA /0.D0,   0.84D0, 0.84D0, 0.84D0, 0.84D0,              
     &             0.66D0, 0.66D0, 1.03D0, 1.03D0, 0.92D0,              
     &             0.66D0/

      ! SOILTB = Coefficient used to convert from surface temperture to  
      !          soil temperature   
      DATA SOILTB /0.D0,   3.6D0,  3.6D0,  3.6D0,  3.6D0,               
     &             8.8D0,  8.8D0,  2.9D0,  2.9D0,  4.4D0,               
     &             8.8D0/
   
      ! SOILAW = Wet biome coefficient   
      DATA SOILAW /0.D0,   2.6D0,  0.03D0, 0.06D0, 0.03D0,              
     &             0.17D0, 0.36D0, 0.36D0, 0.36D0, 0.003D0,             
     &             0.05D0/

      ! SOILAD = Dry biome coefficient  
      DATA SOILAD /0.D0,   8.6D0,  0.22D0, 0.40D0, 0.22D0,              
     &             1.44D0, 2.65D0, 2.65D0, 2.65D0, 0.003D0,             
     &             0.37D0/

      ! SOILEXC = Canopy wind extinction coeff.  
      DATA SOILEXC /0.1D0, 4.D0,   4.D0,   4.D0,   4.D0,                
     &              2.D0,  1.D0,   2.D0,   2.D0,   0.5D0,               
     &              0.1D0/
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
