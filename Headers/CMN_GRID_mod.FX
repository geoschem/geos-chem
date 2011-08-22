!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_GRID_mod.F
!
! !DESCRIPTION: CMN\_SIZE contains size parameters for GEOS-Chem arrays.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_GRID_MOD
! !USES:
!
      ! C Preprocessor #define statements for conditional compilation

#     include "define.h"
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! GRID PARAMETERS
      !
      ! IGLOB      = global longitude dimension
      ! JGLOB      = global latitude dimension
      ! LGLOB      = max number of sigma levels 
      ! IIPAR      = window longitude dimension
      ! JJPAR      = window latitude dimension
      ! LLPAR      = window vertical dimension
      ! LLTROP     = maximum number of tropospheric levels for variable
      !              tropopause
      ! LLTROP_FIX = number of tropospheric levels for offline simulations
      ! PTOP       = model top pressure (mb)
      !
      ! Most of the time, GEOS-CHEM is used for global simulations.
      ! In this case, then IIPAR=IGLOB, JJPAR=JGLOB, LLPAR=LGLOB.
      !
      ! For nested grids, then IIPAR<IGLOB, JJPAR<JGLOB, LLPAR<LGLOB.
      !=================================================================
#if   defined( GCAP ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! GCAP: 4 x 5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 72
      INTEGER, PARAMETER :: JGLOB      = 45
      INTEGER, PARAMETER :: LGLOB      = 23
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
      INTEGER, PARAMETER :: LLPAR      = LGLOB
      INTEGER, PARAMETER :: LLTROP     = 12   
      INTEGER, PARAMETER :: LLTROP_FIX = LLTROP     
      REAL*8,  PARAMETER :: PTOP       = 0.002d0

#elif defined( GEOS_3 ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! GEOS-3: 4 x 5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 72
      INTEGER, PARAMETER :: JGLOB      = 46
      INTEGER, PARAMETER :: LGLOB      = 48            
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 30        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP     = 24        !  -- 30 levels
#else
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP     = 25        !  -- 48 levels
#endif
      INTEGER, PARAMETER :: LLTROP_FIX = 20 
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_3 ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! GEOS-3: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 144
      INTEGER, PARAMETER :: JGLOB      = 91
      INTEGER, PARAMETER :: LGLOB      = 48            
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 30        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP     = 24        !  -- 30 levels
#else
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP     = 25        !  -- 48 levels
#endif 
      INTEGER, PARAMETER :: LLTROP_FIX = 20
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_3 ) && defined( GRID1x1 )

      !-----------------------------------------------------------------
      ! GEOS-3: 1 x 1 
      !-----------------------------------------------------------------
#if   defined( NESTED_CH )
      INTEGER, PARAMETER :: IGLOB      = 81        ! CHINA NESTED 1x1
      INTEGER, PARAMETER :: JGLOB      = 67
      INTEGER, PARAMETER :: LGLOB      = 48           
#elif defined( NESTED_NA )
      INTEGER, PARAMETER :: IGLOB      = 101       ! N. AMERICA NESTED 1x1
      INTEGER, PARAMETER :: JGLOB      = 51    
      INTEGER, PARAMETER :: LGLOB      = 48    
#else
      INTEGER, PARAMETER :: IGLOB      = 360       ! GLOBAL GRID 1x1
      INTEGER, PARAMETER :: JGLOB      = 181
      INTEGER, PARAMETER :: LGLOB      = 48                 

#endif
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 30        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP     = 24        !  -- 30 levels  
#else
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP     = 25        !  -- 48 levels
#endif
      INTEGER, PARAMETER :: LLTROP_FIX = 20
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_4 ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! GEOS-4: 4 x 5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 72
      INTEGER, PARAMETER :: JGLOB      = 46
      INTEGER, PARAMETER :: LGLOB      = 55           
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED ) 
      INTEGER, PARAMETER :: LLPAR      = 30        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP     = 22        !  -- 30 levels
#else
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP     = 23        !  -- 55 levels
#endif
      INTEGER, PARAMETER :: LLTROP_FIX = 17
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_4 ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! GEOS-4: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 144
      INTEGER, PARAMETER :: JGLOB      = 91
      INTEGER, PARAMETER :: LGLOB      = 55
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 30        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP     = 22        !  -- 30 levels
#else 
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP     = 23        !  -- 55 levels
#endif
      INTEGER, PARAMETER :: LLTROP_FIX = 17
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_4 ) && defined( GRID1x125 )

      !-----------------------------------------------------------------
      ! GEOS-4: 1 x 1.2.5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 288
      INTEGER, PARAMETER :: JGLOB      = 181
      INTEGER, PARAMETER :: LGLOB      = 55
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 30        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP     = 22        !  -- 30 levels
#else 
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP     = 23        !  -- 55 levels
#endif
      INTEGER, PARAMETER :: LLTROP_FIX = 17
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_5 ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! GEOS-5: 4 x 5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 72
      INTEGER, PARAMETER :: JGLOB      = 46
      INTEGER, PARAMETER :: LGLOB      = 72         
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED ) 
      INTEGER, PARAMETER :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_5 ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! GEOS-5: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 144
      INTEGER, PARAMETER :: JGLOB      = 91
      INTEGER, PARAMETER :: LGLOB      = 72
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else 
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_5 ) && defined( GRID1x125 )

      !-----------------------------------------------------------------
      ! GEOS-5: 1 x 1.25
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 288
      INTEGER, PARAMETER :: JGLOB      = 181
      INTEGER, PARAMETER :: LGLOB      = 72
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38
#else 
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels    
      INTEGER, PARAMETER :: LLTROP     = 40
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_5 ) && defined( GRID05x0666 )

      !-----------------------------------------------------------------
      ! GEOS-5: 0.5 x 0.666
      !-----------------------------------------------------------------
#if   defined( NESTED_CH )
      INTEGER, PARAMETER :: IGLOB      = 121       ! NESTED CHINA   0.5x0.666
      INTEGER, PARAMETER :: JGLOB      = 133
      INTEGER, PARAMETER :: LGLOB      = 72
#elif defined( NESTED_NA )
      INTEGER, PARAMETER :: IGLOB      = 151       ! NESTED N.AMER. 0.5x0.666
      INTEGER, PARAMETER :: JGLOB      = 121
      INTEGER, PARAMETER :: LGLOB      = 72
#elif defined( NESTED_EU )
      INTEGER, PARAMETER :: IGLOB      = 121       ! NESTED EUROPE  0.5x0.666
      INTEGER, PARAMETER :: JGLOB      = 81
      INTEGER, PARAMETER :: LGLOB      = 72
#endif
      INTEGER, PARAMETER :: IIPAR      = IGLOB  
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38
#else
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( MERRA ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! MERRA: 4 x 5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 72
      INTEGER, PARAMETER :: JGLOB      = 46
      INTEGER, PARAMETER :: LGLOB      = 72         
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED ) 
      INTEGER, PARAMETER :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels 
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( MERRA ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! MERRA: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGLOB      = 144
      INTEGER, PARAMETER :: JGLOB      = 91
      INTEGER, PARAMETER :: LGLOB      = 72
      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else 
      INTEGER, PARAMETER :: LLPAR      = LGLOB     ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#endif

      !-----------------------------------------------------------------
      ! For GEOS 1x1 files
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: I1x1       = 360
      INTEGER, PARAMETER :: J1x1       = 181

      !-----------------------------------------------------------------
      ! For GEOS 05x0666 files
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: I05x0666   = 540 
      INTEGER, PARAMETER :: J05x0666   = 361

      !=================================================================
      ! DISIZE = size (in degrees) of a longitude grid box
      ! DJSIZE = size (in degrees) of a latitude  grid box
      !=================================================================
#if   defined( GRID4x5 ) 
      REAL*8, PARAMETER :: DISIZE = 5.0d0
      REAL*8, PARAMETER :: DJSIZE = 4.0d0
#elif defined( GRID2x25 )
      REAL*8, PARAMETER :: DISIZE = 2.5d0 
      REAL*8, PARAMETER :: DJSIZE = 2.0d0
#elif defined( GRID1x125 )
      REAL*8, PARAMETER :: DISIZE = 1.25d0 
      REAL*8, PARAMETER :: DJSIZE = 1.0d0
#elif defined( GRID1x1 )
      REAL*8, PARAMETER :: DISIZE = 1.0d0 
      REAL*8, PARAMETER :: DJSIZE = 1.0d0
#elif defined( GRID05x0666 )
      REAL*8, PARAMETER :: DISIZE = 2d0/3d0
      REAL*8, PARAMETER :: DJSIZE = 0.5d0
#endif

      END MODULE CMN_GRID_MOD
