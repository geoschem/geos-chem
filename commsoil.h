! $Id: commsoil.h,v 1.4 2006/10/17 17:51:10 bmy Exp $
!
!**********************************************************************
!                                                                     *
!  HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP                       *
!  MODULE FOR SOIL NOx EMISSIONS                                      *
!  by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob               *
!  (Release V2.1)                                                     *
!                                                                     *
!  Contact person: Bob Yantosca (bmy@io.harvard.edu)                  *
!                                                                     *
!**********************************************************************
! NOTES:
! (1 ) Be sure to force double precision with the DBLE function            
!       and the "D" exponent, wherever necessary (bmy, 10/6/99)             
! (2 ) Changed RCS ID tag comment character from "C" to "!" to allow 
!       freeform compilation.  Also added & continuation characters in 
!       column 73 to allow header files to be included in F90 freeform 
!       files. Updated comments, cosmetic changes. (bmy, 6/25/02)
! (3 ) Now use cpp switches to define 1x1 parameters.  Also added
!       space in the #ifdef block for the 1x125 grid (bmy, 12/1/04)
! (4 ) Bug fix: 2681 should be 2861 in NLAND (bmy, 9/22/06)
!**********************************************************************
!
! header file for soil NOx emissions
      
      INTEGER NLAND, NPULSE, NSOIL

! GEOS-CHEM 4 x 5   grid has  1118 grid boxes on land (bmy, 1/23/98)
! GEOS-CHEM 2 x 2.5 grid has  3920 grid boxes on land (bmy, 1/23/98)
! GEOS-CHEM 1 x 1   grid has 17174 grid boxes on land (bmy, 8/7/00)
! GEOS-CHEM 1 x 1   nested China Grid has 2681 land boxes (bmy, 3/11/03)
#if   defined( GRID4x5  )
      PARAMETER( NLAND=1118, NPULSE=3 ) 

#elif defined( GRID2x25 )
      PARAMETER( NLAND=3920, NPULSE=3 ) 

#elif defined( GRID1x125 )
      ! NOTE: NEED TO DEFINE THESE!!!
      PARAMETER( NLAND=????, NPULSE=3 )       

#elif defined( GRID1x1 )

      ! There are 2861  land points for the CHINA nested grid
      ! There are 2118  land points for the N. AMERICA nested grid
      ! There are 17174 land points for the global grid
#if   defined( NESTED_CH )
      !-------------------------------------
      ! Prior to 9/22/06:
      !PARAMETER( NLAND=2681, NPULSE=3 ) 
      !-------------------------------------
      PARAMETER( NLAND=2861, NPULSE=3 ) 
#elif defined( NESTED_NA )
      PARAMETER( NLAND=2118, NPULSE=3 )
#else   
      PARAMETER( NLAND=17174, NPULSE=3 )
#endif

#endif

      PARAMETER (NSOIL=11)      !the defined soil types
! water/desert/ice//Trop. Rain. Forst.//conifers//dry deciduous//
! other deciduous//woodland//grassland//agriculture (other than rice)
! rice paddies//wetland/tundra
      INTEGER INDEXSOIL(2,NLAND) !i,j of the grid
      REAL*8 SOILPULS(NPULSE+1,NLAND) 
	 !tracking of wet/dry & three types of pulsing (Y&L, 94)
      REAL*8 SOILPREP(2,NLAND)  !two month observed precip
      REAL*8 SOILFERT(NLAND)    !ferterlizers
      REAL*8 PULSFACT(NPULSE)   !pulsing factors
      REAL*8 PULSDECAY(NPULSE)  !pulsing decay per timestep
      REAL*8 SOILNOX(IGLOB,JGLOB) !stores output

      INTEGER NCONSOIL(NVEGTYPE) !olson->soil type,nvegtype in commbio.h
      REAL*8 CANOPYNOX(MAXIJ,NTYPE) !track NOx within canopy dry dep.
      REAL*8 SOILTA(NSOIL),SOILTB(NSOIL),SOILAW(NSOIL),SOILAD(NSOIL)
      REAL*8 SOILEXC(NSOIL)     !canopy wind extinction coeff.

      COMMON /SOIL/ SOILNOX,  INDEXSOIL, NCONSOIL, SOILPULS,            &
     &              SOILPREP, SOILFERT,  CANOPYNOX

      ! The correct sequence of PULSFACT is 5, 10, 15 
      DATA PULSFACT  / 5.D0,    10.D0,   15.D0   /

      ! PULSDECAY now contains the correct decay factors from Yienger & Levy
      DATA PULSDECAY / 0.805D0, 0.384D0, 0.208D0 / 

      ! SOILTA = Coefficient used to convert from surface temperture to  
      !          soil temperature     
      DATA SOILTA /0.D0,   0.84D0, 0.84D0, 0.84D0, 0.84D0,              &
     &             0.66D0, 0.66D0, 1.03D0, 1.03D0, 0.92D0,              &
     &             0.66D0/

      ! SOILTB = Coefficient used to convert from surface temperture to  
      !          soil temperature   
      DATA SOILTB /0.D0,   3.6D0,  3.6D0,  3.6D0,  3.6D0,               &
     &             8.8D0,  8.8D0,  2.9D0,  2.9D0,  4.4D0,               &
     &             8.8D0/
   
      ! SOILAW = Wet biome coefficient   
      DATA SOILAW /0.D0,   2.6D0,  0.03D0, 0.06D0, 0.03D0,              &
     &             0.17D0, 0.36D0, 0.36D0, 0.36D0, 0.003D0,             &
     &             0.05D0/

      ! SOILAD = Dry biome coefficient  
      DATA SOILAD /0.D0,   8.6D0,  0.22D0, 0.40D0, 0.22D0,              &
     &             1.44D0, 2.65D0, 2.65D0, 2.65D0, 0.003D0,             &
     &             0.37D0/

      ! SOILEXC = Canopy wind extinction coeff.  
      DATA SOILEXC /0.1D0, 4.D0,   4.D0,   4.D0,   4.D0,                &
     &              2.D0,  1.D0,   2.D0,   2.D0,   0.5D0,               &
     &              0.1D0/

      INTEGER ISOILDIAG, JSOILDIAG
      DATA ISOILDIAG, JSOILDIAG /18,33/ !68% ag
