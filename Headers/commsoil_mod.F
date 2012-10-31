!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: commsoil_mod
!
! !DESCRIPTION: Module COMMSOIL\_MOD contains global variables for the 
!  soil NOx emissions routines.  This has been updated to the new Soil NOx
!  algorithm (2012).
!\\
!\\
! !INTERFACE: 
!
      MODULE COMMSOIL_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, MAXIJ

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
!
! !PUBLIC DATA MEMBERS:
!
      !=================================================================
      ! The following arrays depend on longitude & latitude
      !=================================================================

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! All arrays that depend on (IIPAR,JJPAR) need to be made
      ! ALLOCATABLE.  The sizes of IIPAR, JJPAR will be determined by
      ! the ESMF interface so we cannot assume that they are fixed
      ! parameters. (bmy, 10/30/12)
      !-----------------------------------------------------------------

      ! Soil NOx emissions [molec/cm2/s]
      REAL*8,  ALLOCATABLE :: SOILNOX      (:,:  )

      ! Soil fertilizer 
      REAL*8,  ALLOCATABLE :: SOILFERT     (:,:,:)

      ! Fraction of arid (layer 1) and non-arid (layer 2) land
      REAL*4,  ALLOCATABLE :: CLIM         (:,:,:)
                                         
      ! MODIS landtype
      REAL*4,  ALLOCATABLE :: LAND2        (:,:,:)

      ! Dry period length
      REAL*4,  ALLOCATABLE :: DRYPERIOD    (:,:  )

      ! Pulse factors
      REAL*4,  ALLOCATABLE :: PFACTOR      (:,:  )
      REAL*4,  ALLOCATABLE :: GWET_PREV    (:,:  )

      ! Instantaneous soil NOx and fertilizer
      REAL*8,  ALLOCATABLE :: INST_SOIL    (:,:  )
      REAL*8,  ALLOCATABLE :: INST_FERT    (:,:  )

      ! NOx in the canopy, used in dry deposition
      REAL*8,  ALLOCATABLE :: CANOPYNOX    (:,:  )

      ! Soil NOx deposited N arrays
      REAL*8,  ALLOCATABLE :: DRY_NO2      (:,:  )
      REAL*8,  ALLOCATABLE :: DRY_PAN      (:,:  )
      REAL*8,  ALLOCATABLE :: DRY_HNO3     (:,:  )
      REAL*8,  ALLOCATABLE :: DRY_NH3      (:,:  )
      REAL*8,  ALLOCATABLE :: DRY_NH4      (:,:  )
      REAL*8,  ALLOCATABLE :: DRY_NIT      (:,:  )
      REAL*8,  ALLOCATABLE :: DEP_RESERVOIR(:,:  )
      REAL*8,  ALLOCATABLE :: WET_HNO3     (:,:  ) 
      REAL*8,  ALLOCATABLE :: WET_NH3      (:,:  )
      REAL*8,  ALLOCATABLE :: WET_NH4      (:,:  )
      REAL*8,  ALLOCATABLE :: WET_NIT      (:,:  )

#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice in the std GEOS-Chem is to declare arrays
      ! with IIPAR, JJPAR as fixed parameters. (bmy, 10/30/12)
      !-----------------------------------------------------------------

      ! SOIL NOx emissions [molec/cm2/s]
      REAL*8               :: SOILNOX      (IIPAR,JJPAR    ) 

      ! New Soil NOx model arrays
      REAL*8               :: SOILFERT     (IIPAR,JJPAR,366)

      ! Fraction of arid (layer 1) and non-arid (layer 2) land
      REAL*4               :: CLIM         (IIPAR,JJPAR,2  )      

      ! MODIS landtype
      REAL*4               :: LAND2        (IIPAR,JJPAR,24 )    

      ! Dry period length
      REAL*4               :: DRYPERIOD    (IIPAR,JJPAR    )   

      ! Pulsing factors
      REAL*4               :: PFACTOR      (IIPAR,JJPAR    )   
      REAL*4               :: GWET_PREV    (IIPAR,JJPAR    )   

      ! Instantaneous soil NOx and fertilizer
      REAL*8               :: INST_SOIL    (IIPAR,JJPAR    ) 
      REAL*8               :: INST_FERT    (IIPAR,JJPAR    ) 

      ! NOx in the canopy, used in dry deposition
      REAL*8               :: CANOPYNOX    (MAXIJ,NSOILB   ) 

      ! Soil NOx deposited N arrays
      REAL*8               :: DRY_NO2      (IIPAR,JJPAR    )  
      REAL*8               :: DRY_PAN      (IIPAR,JJPAR    )  
      REAL*8               :: DRY_HNO3     (IIPAR,JJPAR    )
      REAL*8               :: DRY_NH3      (IIPAR,JJPAR    )
      REAL*8               :: DRY_NH4      (IIPAR,JJPAR    )
      REAL*8               :: DRY_NIT      (IIPAR,JJPAR    )
      REAL*8               :: DEP_RESERVOIR(IIPAR,JJPAR    )
      REAL*8               :: WET_HNO3     (IIPAR,JJPAR    )
      REAL*8               :: WET_NH3      (IIPAR,JJPAR    )
      REAL*8               :: WET_NH4      (IIPAR,JJPAR    )
      REAL*8               :: WET_NIT      (IIPAR,JJPAR    )

#endif

      !=================================================================
      ! The following arrays do not depend on longitude & latitude
      !=================================================================

      REAL*8               :: SOILEXC      (NSOILB)

      ! MODIS/Koppen resistance values
      INTEGER              :: SNIMODIS     (NSOILB)
      INTEGER              :: SNIRI        (NSOILB)       
      INTEGER              :: SNIRLU       (NSOILB)      
      INTEGER              :: SNIRAC       (NSOILB)      
      INTEGER              :: SNIRGSS      (NSOILB)     
      INTEGER              :: SNIRGSO      (NSOILB)     
      INTEGER              :: SNIRCLS      (NSOILB)     
      INTEGER              :: SNIRCLO      (NSOILB) 
      INTEGER              :: SNIVSMAX     (NSOILB) 
!
! !REMARKS:
!  Updated to new Soil NOx algorithm (2012).  See:
!  http://wiki.seas.harvard.edu/geos-chem/index.php/Soil_NOx_Emissions
!
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
!  30 Oct 2012 - R. Yantosca - Removed obsolete NLAND parameter, that cannot
!                              be used with the Grid-Independent GEOS-Chem
!  30 Oct 2012 - R. Yantosca - Now make all arrays that depend on lon &
!                              lat into ALLOCATABLE arrays (for GIGC code)
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_commsoil_mod
!
! !DESCRIPTION: Routine SET_COMMSOIL_MOD allocates all module arrays
!  with the longitude and latitude values IIPAR and JJPAR.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_COMMSOIL_MOD( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!  This is used for the Grid-Independent GEOS-Chem.  We cannot assume that
!  IIPAR and JJPAR will be fixed parameters, since these would be determined
!  from the interface to the external GCM.
!                                                                             .
!  May need to add better error checking 
! 
! !REVISION HISTORY: 
!  30 Oct 2012 - R. Yantosca - Now allocate all arrays depending on lon & lat
!  30 Oct 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )

      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! For the GIGC, we need to make these arrays ALLOCATABLE, becuase
      ! the ESMF interface to the GEOS-5 GCM will determine the values
      ! of IIPAR and JJPAR.  We cannot use fixed IIPAR and JJPAR
      ! parameters as is done for the traditional GEOS-Chem.
      !-----------------------------------------------------------------

      ALLOCATE( CANOPYNOX    ( IIPAR*JJPAR, NSOILB ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( SOILNOX      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
                                                    
      ALLOCATE( SOILFERT     ( IIPAR, JJPAR, 366   ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( CLIM         ( IIPAR, JJPAR, 2     ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( LAND2        ( IIPAR, JJPAR, 24    ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( DRYPERIOD    ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( PFACTOR      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( GWET_PREV    ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( INST_SOIL    ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( INST_FERT    ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( DRY_NO2      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( DRY_PAN      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( DRY_HNO3     ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( DRY_NH3      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( DRY_NH4      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
                                                   
      ALLOCATE( DRY_NIT      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( DEP_RESERVOIR( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( WET_HNO3     ( IIPAR, JJPAR        ), STAT=RC ) 
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( WET_NH3      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( WET_NH4      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( WET_NIT      ( IIPAR, JJPAR        ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

#endif

      ! Zero arrays
      SOILNOX       = 0d0
      SOILFERT      = 0d0
      CLIM          = 0e0
      LAND2         = 0e0
      DRYPERIOD     = 0e0
      PFACTOR       = 0e0
      GWET_PREV     = 0e0
      INST_SOIL     = 0d0
      INST_FERT     = 0d0
      CANOPYNOX     = 0d0
      DRY_NO2       = 0d0
      DRY_PAN       = 0d0
      DRY_HNO3      = 0d0
      DRY_NH3       = 0d0
      DRY_NH4       = 0d0
      DRY_NIT       = 0d0
      DEP_RESERVOIR = 0d0
      WET_HNO3      = 0d0 
      WET_NH3       = 0d0
      WET_NH4       = 0d0
      WET_NIT       = 0d0

      END SUBROUTINE SET_COMMSOIL_MOD
!EOC
      END MODULE COMMSOIL_MOD
