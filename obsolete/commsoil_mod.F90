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
  IMPLICIT NONE
  PRIVATE
!
! !DEFINED PARAMETERS: 
!      
  ! The defined soil types (Olson soil types)
  INTEGER, PUBLIC, PARAMETER :: NSOIL = 11

  ! Number of MODIS/Koppen biome types
  INTEGER, PUBLIC, PARAMETER :: NSOILB = 24

  ! MODIS/Koppen resistance values
  INTEGER, PUBLIC, PARAMETER :: SNIMODIS(NSOILB) =              (/ &
       1,    2,    3,    4,    5,    6,    7,    8,    9,   10,    &
      11,   12,   13,   14,   15,   16,   17,   18,   19,   20,    &
      21,   22,   23,   24                                       /)

  INTEGER, PUBLIC, PARAMETER :: SNIRI   (NSOILB) =              (/ &
    9999,  200, 9999, 9999, 9999, 9999,  200,  200,  200,  200,    &
     200,  200,  200,  200,  200,  200,  200,  400,  400,  200,    &
     200,  200, 9999,  200                                       /)

  INTEGER, PUBLIC, PARAMETER :: SNIRLU  (NSOILB) =              (/ &
    9999, 9000, 9999, 9999, 9999, 9999, 9000, 9000, 9000, 9000,    &
    9000, 9000, 9000, 9000, 9000, 1000, 9000, 9000, 9000, 9000,    &
    1000, 9000, 9999, 9000                                       /)

  INTEGER, PUBLIC, PARAMETER :: SNIRAC  (NSOILB) =              (/ &
       0,  300,    0,    0,    0,    0,  100,  100,  100,  100,    &
     100,  100,  100,  100, 2000, 2000, 2000, 2000, 2000, 2000,    &
    2000,  200,  100,  200                                      /)

  INTEGER, PUBLIC, PARAMETER :: SNIRGSS (NSOILB) =              (/ &
       0,    0,  100, 1000,  100, 1000,  350,  350,  350,  350,    &
     350,  350,  350,  350,  500,  200,  500,  500,  500,  500,    &
     200,  150,  400,  150                                      /)

  INTEGER, PUBLIC, PARAMETER :: SNIRGSO (NSOILB) =              (/ &
    2000, 1000, 3500,  400, 3500,  400,  200,  200,  200,  200,    &
     200,  200,  200,  200,  200,  200,  200,  200,  200,  200,    &
     200,  150,  300,  150                                       /)

  INTEGER, PUBLIC, PARAMETER :: SNIRCLS (NSOILB) =              (/ &
    9999, 2500, 9999, 9999, 9999, 9999, 2000, 2000, 2000, 2000,    &
    2000, 2000, 2000, 2000, 2000, 9999, 2000, 2000, 2000, 2000,    &
    9999, 2000, 9999, 2000                                       /)

  INTEGER, PUBLIC, PARAMETER :: SNIRCLO (NSOILB) =              (/ &
    9999, 1000, 1000, 9999, 1000, 9999, 1000, 1000, 1000, 1000,    & 
    1000, 1000, 1000, 1000, 1000, 9999, 1000, 1000, 1000, 1000,    &
    9999, 1000, 9999, 1000                                       /)

  INTEGER, PUBLIC, PARAMETER :: SNIVSMAX(NSOILB) =              (/ &
      10,  100,  100,   10,  100,   10,  100,  100,  100,  100,    &
     100,  100,  100,  100,  100,  100,  100,  100,  100,  100,    &
     100,  100,  100,  100                                       /)
!
! !PUBLIC DATA MEMBERS:
!
  !========================================================================
  ! The following arrays depend on longitude & latitude
  !========================================================================

  INTEGER, PUBLIC                      :: Nx  ! # of lons (x-dimension) in grid
  INTEGER, PUBLIC                      :: Ny  ! # of lats (y-dimension) in grid

  ! Soil NOx emissions [molec/cm2/s]
  REAL*8,  PUBLIC, ALLOCATABLE         :: SOILNOX      (:,:  )

  ! Soil fertilizer 
  REAL*8,  PUBLIC, ALLOCATABLE         :: SOILFERT     (:,:,:)

  ! Fraction of arid (layer 1) and non-arid (layer 2) land
  REAL*4,  PUBLIC, ALLOCATABLE         :: CLIM         (:,:,:)
                                         
  ! MODIS landtype
  REAL*4,  PUBLIC, ALLOCATABLE         :: LAND2        (:,:,:)

  ! Dry period length
  REAL*4,  PUBLIC, ALLOCATABLE         :: DRYPERIOD    (:,:  )

 ! Pulse factors
  REAL*4,  PUBLIC, ALLOCATABLE         :: PFACTOR      (:,:  )
  REAL*4,  PUBLIC, ALLOCATABLE         :: GWET_PREV    (:,:  )

  ! Instantaneous soil NOx and fertilizer
  REAL*8,  PUBLIC, ALLOCATABLE         :: INST_SOIL    (:,:  )
  REAL*8,  PUBLIC, ALLOCATABLE         :: INST_FERT    (:,:  )

  ! NOx in the canopy, used in dry deposition
  REAL*8,  PUBLIC, ALLOCATABLE         :: CANOPYNOX    (:,:  )

  ! Soil NOx deposited N arrays
  REAL*8,  PUBLIC, ALLOCATABLE, TARGET :: DEP_RESERVOIR(:,:  )

  ! Sum N on the fly (ckeller, 14/04/02)
  REAL*8,  PUBLIC, ALLOCATABLE, TARGET :: DRY_TOTN     (:,:  )
  REAL*8,  PUBLIC, ALLOCATABLE, TARGET :: WET_TOTN     (:,:  )
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
!  %%% NOTE: THIS MODULE IS DEPRECATED.  MANY OF THE FIELDS NOW HAVE BEEN %%%
!  %%% MADE OBSOLETE BY THE HEMCO EMISSIONS COMPONENT.  WE NEED TO KEEP   %%%
!  %%% THIS HERE FOR get_ndep_mod.F.  WE MAY BE ABLE TO FOLD VARIABLES    %%%
!  %%% INTO get_ndep_mod.F LATER ON. (bmy, 7/25/14)                       %%%
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
!  25 Jun 2014 - R. Yantosca - Now declare MODIS/Koppen resistances as
!                              PARAMETERS.  This facilitates I/O w/ ESMF.
!  25 Jun 2014 - R. Yantosca - Renamed to commsoil_mod.F90; Use F90 free-format
!  25 Jul 2014 - R. Yantosca - Remove variables made obsolete by HEMCO
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
! !IROUTINE: Init_commsoil
!
! !DESCRIPTION: Routine INIT_COMMSOIL allocates all module arrays
!  with the longitude and latitude values IIPAR and JJPAR.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_COMMSOIL( am_I_Root, arg_Nx, arg_Ny, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    INTEGER, INTENT(IN)  :: arg_Nx      ! # of lons (x-dimension) in grid
    INTEGER, INTENT(IN)  :: arg_Ny      ! # of lats (y-dimension) in grid
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
!  25 Jun 2014 - R. Yantosca - Now use F90 free-format indentation
!  25 Jul 2014 - R. Yantosca - Remove variables made obsolete by HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Store local copies of size dimensions
    Nx = arg_Nx
    Ny = arg_Ny
 
    ! Allocate arrays
    ALLOCATE( CANOPYNOX    ( Nx*Ny, NSOILB ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( SOILNOX      ( Nx, Ny        ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
                                                    
    ALLOCATE( SOILFERT     ( Nx, Ny, 366   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( CLIM         ( Nx, Ny, 2     ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( LAND2        ( Nx, Ny, 24    ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( DRYPERIOD    ( Nx, Ny        ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( PFACTOR      ( Nx, Ny        ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( GWET_PREV    ( Nx, Ny        ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( INST_SOIL    ( Nx, Ny        ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( INST_FERT    ( Nx, Ny        ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( DEP_RESERVOIR( Nx, Ny        ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ! ckeller (14/04/02)
    ALLOCATE( WET_TOTN     ( Nx, Ny        ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( DRY_TOTN     ( Nx, Ny        ), STAT=RC )
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
    DEP_RESERVOIR = 0d0
    WET_TOTN      = 0d0
    DRY_TOTN      = 0d0

  END SUBROUTINE Init_COMMSOIL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!  25 Jun 2014 - R. Yantosca - Now use F90 free-format indentation
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC = GIGC_SUCCESS
      
     ! Deallocate arrays
    IF ( ALLOCATED( CANOPYNOX     ) ) DEALLOCATE( CANOPYNOX     )
    IF ( ALLOCATED( SOILNOX       ) ) DEALLOCATE( SOILNOX       )
    IF ( ALLOCATED( SOILFERT      ) ) DEALLOCATE( SOILFERT      )
    IF ( ALLOCATED( CLIM          ) ) DEALLOCATE( CLIM          )
    IF ( ALLOCATED( LAND2         ) ) DEALLOCATE( LAND2         )
    IF ( ALLOCATED( DRYPERIOD     ) ) DEALLOCATE( DRYPERIOD     )
    IF ( ALLOCATED( PFACTOR       ) ) DEALLOCATE( PFACTOR       )
    IF ( ALLOCATED( GWET_PREV     ) ) DEALLOCATE( GWET_PREV     )
    IF ( ALLOCATED( INST_SOIL     ) ) DEALLOCATE( INST_SOIL     )
    IF ( ALLOCATED( INST_FERT     ) ) DEALLOCATE( INST_FERT     )
    IF ( ALLOCATED( DEP_RESERVOIR ) ) DEALLOCATE( DEP_RESERVOIR )
    IF ( ALLOCATED( DRY_TOTN      ) ) DEALLOCATE( DRY_TOTN      )
    IF ( ALLOCATED( WET_TOTN      ) ) DEALLOCATE( WET_TOTN      )

  END SUBROUTINE Cleanup_COMMSOIL
!EOC
END MODULE COMMSOIL_MOD
