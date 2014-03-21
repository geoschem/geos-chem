!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      PRIVATE
!
! !DEFINED PARAMETERS: 
!      
      ! The defined soil types (Olson soil types)
      INTEGER, PUBLIC, PARAMETER :: NSOIL = 11

      ! Number of MODIS/Koppen biome types
      INTEGER, PUBLIC, PARAMETER :: NSOILB = 24
!
! !PUBLIC DATA MEMBERS:
!
      !=================================================================
      ! The following arrays depend on longitude & latitude
      !=================================================================

      ! Soil NOx emissions [molec/cm2/s]
      REAL*8,  PUBLIC, ALLOCATABLE :: SOILNOX      (:,:  )

      ! Soil fertilizer 
      REAL*8,  PUBLIC, ALLOCATABLE :: SOILFERT     (:,:,:)

      ! Fraction of arid (layer 1) and non-arid (layer 2) land
      REAL*4,  PUBLIC, ALLOCATABLE :: CLIM         (:,:,:)
                                         
      ! MODIS landtype
      REAL*4,  PUBLIC, ALLOCATABLE :: LAND2        (:,:,:)

      ! Dry period length
      REAL*4,  PUBLIC, ALLOCATABLE :: DRYPERIOD    (:,:  )

      ! Pulse factors
      REAL*4,  PUBLIC, ALLOCATABLE :: PFACTOR      (:,:  )
      REAL*4,  PUBLIC, ALLOCATABLE :: GWET_PREV    (:,:  )

      ! Instantaneous soil NOx and fertilizer
      REAL*8,  PUBLIC, ALLOCATABLE :: INST_SOIL    (:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: INST_FERT    (:,:  )

      ! NOx in the canopy, used in dry deposition
      REAL*8,  PUBLIC, ALLOCATABLE :: CANOPYNOX    (:,:  )

      ! Soil NOx deposited N arrays
      REAL*8,  PUBLIC, ALLOCATABLE :: DRY_NO2      (:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: DRY_PAN      (:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: DRY_HNO3     (:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: DRY_NH3      (:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: DRY_NH4      (:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: DRY_NIT      (:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: DEP_RESERVOIR(:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: WET_HNO3     (:,:  ) 
      REAL*8,  PUBLIC, ALLOCATABLE :: WET_NH3      (:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: WET_NH4      (:,:  )
      REAL*8,  PUBLIC, ALLOCATABLE :: WET_NIT      (:,:  )

      !=================================================================
      ! The following arrays do not depend on longitude & latitude
      !=================================================================

      REAL*8,  PUBLIC              :: SOILEXC      (NSOILB)

      ! MODIS/Koppen resistance values
      INTEGER, PUBLIC              :: SNIMODIS     (NSOILB)
      INTEGER, PUBLIC              :: SNIRI        (NSOILB)       
      INTEGER, PUBLIC              :: SNIRLU       (NSOILB)      
      INTEGER, PUBLIC              :: SNIRAC       (NSOILB)      
      INTEGER, PUBLIC              :: SNIRGSS      (NSOILB)     
      INTEGER, PUBLIC              :: SNIRGSO      (NSOILB)     
      INTEGER, PUBLIC              :: SNIRCLS      (NSOILB)     
      INTEGER, PUBLIC              :: SNIRCLO      (NSOILB) 
      INTEGER, PUBLIC              :: SNIVSMAX     (NSOILB) 
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: Init_COMMSOIL
      PUBLIC :: Cleanup_COMMSOIL
!
! !REMARKS:
!  Updated to new Soil NOx algorithm (2012).  See:
!  http://wiki.seas.harvard.edu/geos-chem/index.php/Soil_NOx_Emissions
!                                                                             .
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%% NOTE: THIS MODULE IS DEPRECATED.  WE MAY END UP PLACING THE MODULE %%%
!  %%% VARIABLES INTO A DERIVED TYPE FOR THE NEW EMISSIONS COMPONENT IN   %%%
!  %%% THE NEAR FUTURE. (bmy, 11/19/12)                                   %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
! !IROUTINE: init_commsoil
!
! !DESCRIPTION: Routine INIT\_COMMSOIL allocates all module arrays
!  with the longitude and latitude values IIPAR and JJPAR.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_COMMSOIL( am_I_Root, RC )
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

      END SUBROUTINE Init_COMMSOIL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_commsoil
!
! !DESCRIPTION: Subroutine CLEANUP\_COMMSOIL deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_COMMSOIL( am_I_Root, RC )
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
      INTEGER, INTENT(OUT) :: RC          ! Return code
! 
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS
      
      ! Deallocate arrays

      IF( ALLOCATED( CANOPYNOX     ) ) DEALLOCATE( CANOPYNOX     )
      IF( ALLOCATED( SOILNOX       ) ) DEALLOCATE( SOILNOX       )
      IF( ALLOCATED( SOILFERT      ) ) DEALLOCATE( SOILFERT      )
      IF( ALLOCATED( CLIM          ) ) DEALLOCATE( CLIM          )
      IF( ALLOCATED( LAND2         ) ) DEALLOCATE( LAND2         )
      IF( ALLOCATED( DRYPERIOD     ) ) DEALLOCATE( DRYPERIOD     )
      IF( ALLOCATED( PFACTOR       ) ) DEALLOCATE( PFACTOR       )
      IF( ALLOCATED( GWET_PREV     ) ) DEALLOCATE( GWET_PREV     )
      IF( ALLOCATED( INST_SOIL     ) ) DEALLOCATE( INST_SOIL     )
      IF( ALLOCATED( INST_FERT     ) ) DEALLOCATE( INST_FERT     )
      IF( ALLOCATED( DRY_NO2       ) ) DEALLOCATE( DRY_NO2       )
      IF( ALLOCATED( DRY_PAN       ) ) DEALLOCATE( DRY_PAN       )
      IF( ALLOCATED( DRY_HNO3      ) ) DEALLOCATE( DRY_HNO3      )
      IF( ALLOCATED( DRY_NH3       ) ) DEALLOCATE( DRY_NH3       )
      IF( ALLOCATED( DRY_NH4       ) ) DEALLOCATE( DRY_NH4       )
      IF( ALLOCATED( DRY_NIT       ) ) DEALLOCATE( DRY_NIT       )
      IF( ALLOCATED( DEP_RESERVOIR ) ) DEALLOCATE( DEP_RESERVOIR )
      IF( ALLOCATED( WET_HNO3      ) ) DEALLOCATE( WET_HNO3      ) 
      IF( ALLOCATED( WET_NH3       ) ) DEALLOCATE( WET_NH3       )
      IF( ALLOCATED( WET_NH4       ) ) DEALLOCATE( WET_NH4       )
      IF( ALLOCATED( WET_NIT       ) ) DEALLOCATE( WET_NIT       )

      END SUBROUTINE Cleanup_COMMSOIL
!EOC
      END MODULE COMMSOIL_MOD
