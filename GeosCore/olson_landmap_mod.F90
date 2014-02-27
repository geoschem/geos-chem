!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: olson_landmap_mod
!
! !DESCRIPTION: Module OLSON\_LANDMAP\_MOD reads the Olson land map and
!  computes the IREG, ILAND, and IUSE arrays.  This module was written to
!  facilitate Grid-Independent GEOS-Chem development while still keeping 
!  backwards compatibility with existing legacy code.  It replaces the old 
!  routine rdland.F.
!\\
!\\
! !INTERFACE: 
!
MODULE Olson_LandMap_Mod
!
! !USES:
!
  USE CMN_GCTM_MOD                      ! Physical constants
  USE CMN_SIZE_MOD                      ! Size parameters
  USE DIRECTORY_MOD                     ! Disk directory paths   
  USE ERROR_MOD                         ! Error checking routines
  USE GRID_MOD                          ! Horizontal grid definition
  USE LOGICAL_MOD                       ! Logical switches
  USE MAPPING_MOD                       ! Mapping weights & areas

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_Olson_Landmap
  PUBLIC  :: Compute_Olson_Landmap
  PUBLIC  :: Cleanup_Olson_LandMap
!
! !REMARKS:
!  The Olson land types are as follows:
!  ============================================================================
!   0 Water              25 Deciduous           50 Desert
!   1 Urban              26 Deciduous           51 Desert
!   2 Shrub              27 Conifer             52 Steppe
!   3 ---                28 Dwarf forest        53 Tundra
!   4 ---                29 Trop. broadleaf     54 rainforest
!   5 ---                30 Agricultural        55 mixed wood/open
!   6 Trop. evergreen    31 Agricultural        56 mixed wood/open
!   7 ---                32 Dec. woodland       57 mixed wood/open
!   8 Desert             33 Trop. rainforest    58 mixed wood/open
!   9 ---                34 ---                 59 mixed wood/open
!  10 ---                35 ---                 60 conifers
!  11 ---                36 Rice paddies        61 conifers
!  12 ---                37 agric               62 conifers
!  13 ---                38 agric               63 Wooded tundra
!  14 ---                39 agric.              64 Moor
!  15 ---                40 shrub/grass         65 coastal
!  16 Scrub              41 shrub/grass         66 coastal
!  17 Ice                42 shrub/grass         67 coastal
!  18 ---                43 shrub/grass         68 coastal
!  19 ---                44 shrub/grass         69 desert
!  20 Conifer            45 wetland             70 ice
!  21 Conifer            46 scrub               71 salt flats
!  22 Conifer            47 scrub               72 wetland
!  23 Conifer/Deciduous  48 scrub               73 water
!  24 Deciduous/Conifer  49 scrub
!                                                                             .
!                                                                             .
!  Arrays computed by olson_landmap_mod.F90
!  ============================================================================
!  (1) IREG   (in CMN_DEP_mod.F): # of Olson land types per G-C grid box 
!  (2) ILAND  (in CMN_DEP_mod.F): List of all Olson land types in G-C grid box
!  (3) IUSE   (in CMN_DEP_mod.F): Coverage of each Olson type in G-C grid box
!  (4) IJREG  (in CMN_VEL_mod.F): %%%%% OBSOLETE: NOW REPLACED BY IREG  %%%%%
!  (5) IJLAND (in CMN_VEL_mod.F): %%%%% OBSOLETE: NOW REPLACED BY ILAND %%%%%
!  (6) IJUSE  (in CMN_VEL_mod.F): %%%%% OBSOLETE: NOW REPLACED BY IUSE  %%%%%
!  (7) FRCLND (in CMN_DEP_mod.F): Fraction of G-C grid box that is not water
!                                                                             .
!  NOTES: 
!  (1) IREG,  ILAND,  IUSE  are used by the soil NOx emissions routines
!  (2) IJREG, IJLAND, IJUSE are used by the drydep routines (legacy code)
!  (3) FRCLND               is  used by various GEOS-Chem routines
!                                                                             .
!                                                                             .
!  BUG IN THE OLD "rdland.F" FOR 2 X 2.5 DEGREE RESOLUTION
!  ============================================================================
!  This module ("olson_landmap_mod.F") replaces the old routine "rdland.F", 
!  which previously read in the Olson landtype data from the ASCII format
!  file named "vegtype.global".  There used to be a different "vegtype.global"
!  file for each different horizontal grid resolution.
!                                                                             .
!  The "vegtype.global" stored the following quantities, such that values
!  for a single grid box were saved on a single line:
!                                                                             .
!    I, J, IREG(I,J), ILAND(I,J,K), IUSE(I,J,K)  (where K=1,IREG(I,J))
!                                                                             .
!  Routine "rdland.F" reads these quantities from "vegtype.global" assuming 
!  there were 20 integer characters on a single line (i.e. using Fortran
!  FORMAT '(20i4)').   However, ~ 12 lines of the 2 x 2.5 "vegtype.global"
!  file contained more than 20 integer values.  This caused "rdland.F", 
!  to read in the values from these lines improperly, which in turn caused
!  the IREG, ILAND, IUSE, IJREG, IJLAND, IJUSE, and FRCLND arrays to be
!  improperly initialized for the grid boxes corresponding to these
!  lines in the "vegtype.global" file.
!                                                                             .
!  Bob Yantosca has validated that "olson_landmap_mod.F" returns results
!  100% identical to the "vegtype.global" file.  Therefore, if you want
!  to compare the output of model simulations using "olson_landmap_mod.F" 
!  the output of simulations using "rdland.F", you will see a slight 
!  difference in the MCL lifetime and tracer concentrations.
!                                                                             .
!  If you need to run a GEOS-Chem simulation with an older version of the
!  code using "rdland.F", then this bug may be corrected by changing the
!  line of code:
!                                                                             .
!      101  FORMAT(20I4)
!                                                                             .
!  to:
!                                                                             .
!     #if   defined( GRID2x25 )
!      101  FORMAT(25I4)
!     #else
!      100  FORMAT(20I4)
!     #endif
!                                                                             .
!  This is more or less a moot point, as "olson_landmap_mod.F" will be
!  installed into GEOS-Chem v9-01-03 and higher versions.
!                                                                             .
!                                                                             .
!  NOTE FOR 0.5 x 0.666 grids
!  ============================================================================
!  As of 21 Mar 2012, the IUSE values computed by "olson_landmap_mod.F90"
!  may slightly differ from those specified in the "vegtype.global" files
!  for 0.5 x 0.666 nested grids.  We attribute this to roundoff error caused
!  by the the longitude spacing being an irrational number (0.6666666...). 
!  We are still investigating.
!
! !REVISION HISTORY:
!  13 Mar 2012 - R. Yantosca - Initial version
!  19 Mar 2012 - R. Yantosca - Minor last-minute bug fixes
!  21 Mar 2012 - R. Yantosca - Now use REAL*4 for computations
!  22 Mar 2012 - R. Yantosca - Now read surface area from the file
!  22 Mar 2012 - R. Yantosca - Now make lon, lat, OLSON, A_CM2 allocatable
!  22 Mar 2012 - R. Yantosca - Now define I_OLSON, J_OLSON, N_OLSON, D_LON,
!                              and D_LAT in routine Init_Olson_LandMap
!  27 Mar 2012 - R. Yantosca - Now reference USE_OLSON_2001 from logical_mod.F
!  02 Apr 2012 - R. Yantosca - Now reference mapping_mod.F90
!  02 Apr 2012 - R. Yantosca - Moved routine GET_MAP_WT to mapping_mod.F90
!  02 Apr 2012 - R. Yantosca - Now Save mapping info for later use
!  09 Apr 2012 - R. Yantosca - Removed IJREG, IJUSE, IJLAND; these are now
!                              replaced by IREG, IUSE, ILAND arrays
!  09 Apr 2012 - R. Yantosca - Removed reference to CMN_VEL_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Scalars
  INTEGER              :: I_OLSON       ! # of lons (0.5 x 0.5)
  INTEGER              :: J_OLSON       ! # of lats (0.5 x 0.5)
  INTEGER              :: N_OLSON       ! Number of Olson land types 
  REAL*8               :: D_LON         ! Delta longitude, Olson grid [degrees]
  REAL*8               :: D_LAT         ! Delta latitude,  Olson grid [degrees]

  ! Arrays
  REAL*4,  ALLOCATABLE :: lon  (:    )  ! Lon centers, Olson grid [degrees]
  REAL*4,  ALLOCATABLE :: lat  (  :  )  ! Lat centers, Olson grid [degrees]
  INTEGER, ALLOCATABLE :: OLSON(:,:,:)  ! Olson land types
  REAL*4,  ALLOCATABLE :: A_CM2(:,:,:)  ! Surface areas [cm2]

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_olson_landmap
!
! !DESCRIPTION: Subroutine COMPUTE\_OLSON\_LANDMAP computes the GEOS-Chem
!  arrays IREG, ILAND, IUSE (and corresponding 1-D arrays IJREG, IJLAND, 
!  IJUSE) on-the-fly from the Olson Land map file.  This routine, which is
!  intended to facilitate the Grid-Independent GEOS-Chem, replaces
!  the old rdland.F, which read from pre-computed "vegtype.global" files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Olson_LandMap( am_I_Root, mapping, State_Met )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)    :: am_I_Root    ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MapWeight), POINTER       :: mapping(:,:) ! "fine" -> "coarse" mapping
    TYPE(MetState),  INTENT(INOUT) :: State_Met    ! Meteorology State object
!
! !REMARKS:
!  This routine supplies arrays that are required for legacy code routines:
!  (1) IREG,  ILAND,  IUSE are used by the Soil NOx routines
!  (2) IJREG, IJLAND, IJUSE are used by the dry deposition routines
! 
! !REVISION HISTORY: 
!  13 Mar 2012 - R. Yantosca - Initial version
!  19 Mar 2012 - R. Yantosca - Reorder ILAND, IUSE, IJLAND, IJUSE to be
!                              consistent w/ the leaf area indices
!  19 Mar 2012 - R. Yantosca - Compute the FRCLND array (from CMN_DEP_mod.F)
!  21 Mar 2012 - R. Yantosca - Now use REAL*4 for computation, to reduce
!                              roundoff errors at high-resolution
!  22 Mar 2012 - R. Yantosca - Now get surface area directly from variable
!                              A_CM2 (read from disk) instead of computing it
!  02 Apr 2012 - R. Yantosca - Now pass MAP (mapping weight object) via the
!                              arg list, to save the mapping info for later
!  09 Apr 2012 - R. Yantosca - Remove IJLOOP variable
!  09 Apr 2012 - R. Yantosca - Now do not compute IJREG, IJLAND, IJUSE; these
!                              are replaced by IREG, ILAND, IUSE arrays
!  17 Apr 2012 - R. Yantosca - Rename "map" object to "mapping" to avoid name
!                              confusion with an F90 intrinsic function
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  29 Nov 2012 - R. Yantosca - Added am_I_Root argument
!  12 Dec 2012 - R. Yantosca - Now get IREG, ILAND, IUSE from State_Met
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL :: isGlobal
    INTEGER :: I,         J,         II,       III
    INTEGER :: JJ,        T,         N,        type
    INTEGER :: uniqOlson, sumIuse,   C,        IG
    REAL*4  :: xedge_w,   xedge_e,   yedge_s,  yedge_n
    REAL*4  :: xedgeC_w,  xedgeC_e,  yedgeC_s, yedgeC_n
    REAL*4  :: dxdy,      dxdy4,     mapWt,    area
    REAL*4  :: sumArea
    
    ! Generic arrays
    INTEGER :: maxIuse(1)
    
    ! Arrays on the Olson land map NATIVE GRID
    INTEGER :: indLon  (I_OLSON                     ) ! Index array for lons
    INTEGER :: shiftLon(I_OLSON                     ) ! Shifted indLon array
    REAL*4  :: lonedge (I_OLSON+1                   ) ! Lon edges   [degrees]
    REAL*4  :: latedge (          J_OLSON+1         ) ! Lat edges   [degrees]
    
    ! Arrays on the GEOS-CHEM GRID                 
    INTEGER :: ctOlson (IIPAR,    JJPAR, 0:N_OLSON-1) ! Count of land types/box
    REAL*4  :: frOlson (IIPAR,    JJPAR, 0:N_OLSON-1) ! Frac of land types/box
    INTEGER :: ordOlson(IIPAR,    JJPAR, 0:N_OLSON-1) ! Order of land types

    ! Pointers
    INTEGER, POINTER :: IREG(:,:)
    INTEGER, POINTER :: ILAND(:,:,:)
    INTEGER, POINTER :: IUSE(:,:,:)
    REAL*8,  POINTER :: FRCLND(:,:)

    !======================================================================
    ! NATIVE GRID parameters (i.e. 0.5 x 0.5 "GENERIC")
    !======================================================================

    ! Be lazy, construct lon edges from lon centers
    DO I = 1, I_OLSON
       lonedge(I)      = DBLE( lon(I) ) - ( D_LON * 0.5d0 )
       indLon(I)       = I
    ENDDO
    lonedge(I_OLSON+1) = lonedge(I_OLSON) + D_LON
    
    ! Be lazy, construct lat edges from lat centers
    DO J = 1, J_OLSON
       latedge(J)      = DBLE( lat(J) ) - ( D_LAT * 0.5d0 )
    ENDDO
    latedge(J_OLSON+1) = latedge(J_OLSON) + D_LAT
    
    ! Shift longitudes by 2 degrees to the west for date-line handling
    shiftLon           = CSHIFT( indLon, -20 )

    !======================================================================
    ! Initialize variables outside of the main loop 
    !======================================================================

    ! Initialize pointers
    IREG    => State_Met%IREG
    ILAND   => State_Met%ILAND
    IUSE    => State_Met%IUSE
    FRCLND  => State_Met%FRCLND

    IREG     = 0
    ILAND    = 0
    IUSE     = 0
    FRCLND   = 1000e0
    ctOlson  = 0
    frOlson  = 0e0
    ordOlson = -999
    isGlobal = ( .not. ITS_A_NESTED_GRID() )

    !======================================================================
    ! Loop over all GEOS-CHEM GRID BOXES and initialize variables
    !======================================================================
    !$OMP PARALLEL DO                                                  &
    !$OMP DEFAULT( SHARED )                                            &
    !$OMP PRIVATE( I,        J,         xedgeC_w, yedgeC_s, xedgeC_e ) &
    !$OMP PRIVATE( yedgeC_n, dxdy4,     sumArea,  JJ,       III      ) &
    !$OMP PRIVATE( dxdy,     mapWt,     II,       xedge_w,  yedge_s  ) &
    !$OMP PRIVATE( xedge_e,  yedge_n,   area,     type,     maxIuse  ) &
    !$OMP PRIVATE( sumIUse,  uniqOlson, C,        IG                 )
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Global lon index (needed for when running in ESMF)
       IG = I + I_LO - 1

       ! Edges of this GEOS-CHEM GRID box
       xedgeC_w  = GET_XEDGE( I,   J,   1 )          ! W edge
       yedgeC_s  = GET_YEDGE( I,   J,   1 )          ! S edge
       xedgeC_e  = GET_XEDGE( I+1, J,   1 )          ! E edge
       yedgeC_n  = GET_YEDGE( I,   J+1, 1 )          ! N edge
       
       ! "Area" of the GEOS-CHEM GRID box in degrees (DLON * DLAT)
       dxdy4     = ( xedgeC_e - xedgeC_w ) * ( yedgeC_n - yedgeC_s )
     
       ! Zero the summing array
       sumArea   = 0e0

       ! Reset counter of olson land types found per box
       uniqOlson = 0e0

       ! Counter for mapping object
       C         = 0

       !===================================================================
       ! Find each 0.5 x 0.5 NATIVE GRID BOX that fits into the GEOS-CHEM
       ! GRID BOX.  Keep track of the land types and coverage fractions.
       !===================================================================
       DO JJ  = 1, J_OLSON
       DO III = 1, I_OLSON

          ! Initialize
          dxdy       = 0e0
          mapWt      = 0e0
       
          ! Find the NATIVE GRID longitude index for use below.  Account for 
          ! the first GEOS-CHEM GRID box, which straddles the date line.
          IF ( isGlobal .and.  IG == 1 ) THEN
             II      = shiftLon(III)
          ELSE
             II      = indLon(III)
          ENDIF
          
          ! Edges of this NATIVE GRID box
          xedge_w    = lonedge(II  )                ! W edge
          yedge_s    = latedge(JJ  )                ! S edge
          xedge_e    = lonedge(II+1)                ! E edge
          yedge_n    = latedge(JJ+1)                ! N edge

          ! Because the first GEOS-CHEM GRID BOX straddles the date line,
          ! we have to adjust the W and E edges of the NATIVE GRID BOX to
          ! be in monotonically increasing order.  This will prevent
          ! erronous results from being returned by GET_MAP_WT below.
          IF ( isGlobal .and. IG == 1 .and. II >= shiftLon(1) )  THEN
             xedge_w = xedge_w - 360e0
             xedge_e = xedge_e - 360e0
          ENDIF
         
          ! "Area" of the GEOS-CHEM GRID BOX in degrees (DLON * DLAT)
          dxdy       = ( xedge_e - xedge_w ) * ( yedge_n - yedge_s )

          ! Get the mapping weight (i.e. The fraction of the NATIVE 
          ! GRID BOX that lies w/in the GEOS-CHEM GRID BOX)
          CALL GET_MAP_WT( xedge_w, xedge_e, xedgeC_w, xedgeC_e,  &
                           yedge_s, yedge_n, yedgeC_s, yedgeC_n,  &
                           mapWt                                 )

          ! Skip unless part (or all) of the NATIVE GRID BOX
          ! actually fits into the GEOS-CHEM GRID BOX
          IF ( mapWt <= 0e0 .or. mapWt > 1e0 ) CYCLE

          ! Area of the NATIVE GRID BOX that lies w/in the GEOS-CHEM GRID BOX
          area              = A_CM2(II,JJ,1) * mapWt
           
          ! Keep a total of the area
          sumArea           = sumArea + area

          ! Olson land map type on the NATIVE GRID
          type              = OLSON(II,JJ,1)
           
          ! Increment count of Olson types
          ctOlson(I,J,type) = ctOlson(I,J,type) + 1

          ! Add area covered by this olson type
          frOlson(I,J,type) = frOlson(I,J,type) + area

          ! Preserve ordering for backwards-compatibility w/ LAI data
          IF ( ordOlson(I,J,type) < 0 ) THEN 

             ! Counter of land types we have encountered for the first time
             uniqOlson          = uniqOlson + 1

             ! Record the order in which this land type was first encountered
             ordOlson(I,J,type) = uniqOlson

          ENDIF

          ! Save mapping information for later use in modis_lai_mod.F90
          ! in order to prepare the State_Met%XLAI array for use with the 
          ! legacy dry-deposition and soil NOx emissions codes.
          C                     = C + 1
          mapping(I,J)%count    = C
          mapping(I,J)%II(C)    = II
          mapping(I,J)%JJ(C)    = JJ
          mapping(I,J)%olson(C) = type
          mapping(I,J)%area(C)  = area
          mapping(I,J)%sumarea  = sumarea

       ENDDO
    ENDDO

       !===================================================================
       ! Construct GEOS-Chem type output arrays from the binning that we 
       ! just have completed.  Preserve the ordering from "vegtype.global"
       ! for backwards compatibility w/ existing code.
       !===================================================================
     
       ! Land type index for ILAND & IUSE
       maxIUse = 0

       ! Loop over all land types
       DO T = 0, N_OLSON-1

          ! Save the ordering of Olson land types for later use 
          ! by routines in the module modis_lai_mod.F90
          mapping(I,J)%ordOlson(T) = ordOlson(I,J,T)

          ! Normalize the land type coverage 
          frOlson(I,J,T)                =  &
               INT( ( ( frOlson(I,J,T) / sumArea ) * 1e3 ) + 0.5e0 )
 
          ! If land type T is represented in this box ...
          IF ( ctOlson(I,J,T) > 0 .and. ordOlson(I,J,T) > 0 ) THEN 
 
             ! Increment the count of Olson types in the box 
             IREG(I,J)                  = IREG(I,J) + 1
             
             ! Save land type into ILAND
             ILAND(I,J,ordOlson(I,J,T)) = T
             
             ! Save the fraction (in mils) of this land type
             IUSE(I,J,ordOlson(I,J,T))  = frOlson(I,J,T)

          ENDIF
       ENDDO

       ! Land type with the largest coverage in the GEOS-CHEM GRID BOX
       maxIuse = MAXLOC( IUSE( I, J, 1:IREG(I,J) ) )

       ! Sum of all land types in the GEOS-CHEM GRID BOX (should be 1000)
       sumIUse = SUM   ( IUSE( I, J, 1:IREG(I,J) ) )

       ! Make sure everything adds up to 1000.  If not, then adjust
       ! the land type w/ the largest coverage accordingly.
       ! This follows the algorithm from "regridh_lai.pro".
       IF ( sumIUse /= 1000 ) THEN
          IUSE(I,J,maxIUse) = IUSE(I,J,maxIUse) &
                            + ( 1000 - sumIUse )
       ENDIF
      
       ! Loop over land types in the GEOS-CHEM GRID BOX
       DO T = 1, IREG(I,J)

          ! If the current Olson land type is water (type 0),
          ! subtract the coverage fraction (IUSE) from FRCLND.
          IF ( ILAND(I,J,T) == 0 ) THEN
             FRCLND(I,J) = FRCLND(I,J)  - IUSE(I,J,T)
          ENDIF
       ENDDO

       ! Normalize FRCLND into the range of 0-1
       ! NOTE: Use REAL*4 for backwards compatibility w/ existing code!
       FRCLND(I,J) = FRCLND(I,J) / 1000e0

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  
!### Save code here for debugging
!###    do j = 1, jjpar
!###    do i = 1, iipar
!###       write( 800, '(2i5, f13.6)' ) i, j, frclnd(i,j)
!###       
!###       write( 810, '(25i4)'       ) i, j, ireg(i,j),                 &
!###                                    ( iland(i,j,t), t=1,ireg(i,j) ), &
!###                                    ( iuse (i,j,t), t=1,ireg(i,j) )
!###    enddo
!###    enddo
!###   
!###    ! ### DEBUG OUTPUT
!###    C = map(23,34)%count
!###    print*, '### count   : ', C
!###    print*, '### II      : ', map(23,34)%II(1:C)
!###    print*, '### JJ      : ', map(23,34)%JJ(1:C)
!###    print*, '### area    : ', map(23,34)%area(1:C)
!###    print*, '### sumarea : ', map(23,34)%sumarea

  END SUBROUTINE Compute_Olson_LandMap
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_olson_landmap
!
! !DESCRIPTION: Subroutine INIT\_OLSON\_LANDMAP reads Olson land map 
! information from disk (in netCDF format).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Olson_LandMap( am_I_Root, DATA_DIR_1x1 )
!
! !USES:
!
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close
    
    IMPLICIT NONE
    
#   include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN) :: am_I_Root
    CHARACTER(LEN=255), INTENT(IN) :: DATA_DIR_1x1
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!
! !REVISION HISTORY:
!  13 Mar 2012 - R. Yantosca - Initial version
!  22 Mar 2012 - R. Yantosca - Also read in surface areas [m2] from file
!  27 Mar 2012 - R. Yantosca - Now read the "units" attribute of each variable
!  27 Mar 2012 - R. Yantosca - Now echo file I/O status info to stdout
!  27 Mar 2012 - R. Yantosca - Now can read Olson 1992 or Olson 2001 land map
!  29 Nov 2012 - R. Yantosca - Add am_I_Root to the argument list
!  26 Feb 2013 - M. Long     - Now pass DATA_DIR_1x1 via the argument list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !======================================================================
    ! Variable declarations
    !======================================================================
    
    ! Scalars
    INTEGER            :: I, J               ! Loop indices
    INTEGER            :: fId                ! netCDF file ID
    INTEGER            :: as                 ! Allocation status
    
    ! Character strings
    CHARACTER(LEN=255) :: nc_dir             ! netCDF directory name
    CHARACTER(LEN=255) :: nc_file            ! netCDF file name
    CHARACTER(LEN=255) :: nc_path            ! netCDF path name
    CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
    CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
    CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
     
    ! Arrays for netCDF start and count values
    INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
    INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays 
     
    !======================================================================
    ! Initialize variables
    !======================================================================
    IF ( USE_OLSON_2001 ) THEN

       !--------------------------------
       ! Settings for Olson 2001 grid
       !--------------------------------
       I_OLSON = 1440                                     ! # lons (0.25x0.25)
       J_OLSON = 720                                      ! # lats (0.25x0.25)
       N_OLSON = 74                                       ! # of land types
       D_LON   = 0.25d0                                   ! Delta lon [degrees]
       D_LAT   = 0.25d0                                   ! Delta lat [degrees]
       nc_file = 'Olson_2001_Land_Map.025x025.generic.nc' ! Input file name

    ELSE

       !--------------------------------
       ! Settings for Olson 1992 grid
       !--------------------------------
       I_OLSON = 720                                      ! # lons (0.5x0.5)
       J_OLSON = 360                                      ! # lats (0.5x0.5)
       N_OLSON = 74                                       ! # of land types
       D_LON   = 0.5d0                                    ! Delta lon [degrees]
       D_LAT   = 0.5d0                                    ! Delta lat [degrees]
       nc_file = 'Olson_1992_Land_Map.05x05.generic.nc'   ! Input file name

    ENDIF

    ! Allocate arrays
    ALLOCATE( lon  ( I_OLSON             ), STAT=as )  
    ALLOCATE( lat  ( J_OLSON             ), STAT=as )
    ALLOCATE( OLSON( I_OLSON, J_OLSON, 1 ), STAT=as ) 
    ALLOCATE( A_CM2( I_OLSON, J_OLSON, 1 ), STAT=as )

    !======================================================================
    ! Open and read data from the netCDF file
    !======================================================================

    ! Construct file path from directory & file name
    nc_dir  = TRIM( DATA_DIR_1x1 ) // 'Olson_Land_Map_201203/'
    nc_path = TRIM( nc_dir ) // TRIM( nc_file )

    ! Open file for read
    CALL Ncop_Rd( fId, TRIM(nc_path) )
     
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 100 ) REPEAT( '%', 79 )
       WRITE( 6, 110 ) TRIM(nc_file)
       WRITE( 6, 120 ) TRIM(nc_dir)
    ENDIF

    !----------------------------------------
    ! VARIABLE: lon
    !----------------------------------------
     
    ! Variable name
    v_name = "lon"
    
    ! Read lon from file
    st1d   = (/ 1       /)
    ct1d   = (/ I_OLSON /)
    CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )
 
    ! Read the lon:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)     
    ENDIF

    !----------------------------------------
    ! VARIABLE: lat
    !----------------------------------------
    
    ! Variable name
    v_name = "lat"
    
    ! Read lat from file
    st1d   = (/ 1       /)
    ct1d   = (/ J_OLSON /)
    CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )
     
    ! Read the lat:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val) 
    ENDIF

    !----------------------------------------
    ! VARIABLE: OLSON
    !----------------------------------------
    
    ! Variable name
    v_name = "OLSON"
    
    ! Read OLSON from file
    st3d   = (/ 1,       1,       1 /)
    ct3d   = (/ I_OLSON, J_OLSON, 1 /)
    CALL NcRd( OLSON, fId, TRIM(v_name), st3d, ct3d )

    ! Read the OLSON:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val) 
    ENDIF

    !----------------------------------------
    ! VARIABLE: DXYP 
    ! Convert from m2 to cm2; store as A_CM2 
    !----------------------------------------
    
    ! Variable name
    v_name = "DXYP"
    
    ! Read OLSON from file
    st3d   = (/ 1,       1,       1 /)
    ct3d   = (/ I_OLSON, J_OLSON, 1 /)
    CALL NcRd( A_CM2, fId, TRIM(v_name), st3d, ct3d )
    
    ! Read the DXYP:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val) 
    ENDIF

    ! Convert from [m2] to [cm2]
    A_CM2  = A_CM2 * 1e4

    !=================================================================
    ! Cleanup and quit
    !=================================================================
    
    ! Close netCDF file
    CALL NcCl( fId )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 140 )
       WRITE( 6, 100 ) REPEAT( '%', 79 )
    ENDIF

    ! FORMAT statements
100 FORMAT( a                                              )
110 FORMAT( '%% Opening file  : ',         a               )
120 FORMAT( '%%  in directory : ',         a, / , '%%'     )
130 FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
140 FORMAT( '%% Successfully closed file!'                 )

  END SUBROUTINE Init_Olson_LandMap
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_olson_landmap
!
! !DESCRIPTION: Subroutine CLEANUP\_OLSON\_LANDMAP deallocates all allocated
!  global module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Olson_LandMap( am_I_Root )
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN) :: am_I_Root   ! Are we on the root CPU?
!
! !REVISION HISTORY:'
!  22 Mar 2012 - R. Yantosca - Initial version
!  29 Nov 2012 - R. Yantosca - Add am_I_Root as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( lon   ) ) DEALLOCATE( lon   )
    IF ( ALLOCATED( lat   ) ) DEALLOCATE( lat   )
    IF ( ALLOCATED( OLSON ) ) DEALLOCATE( OLSON )
    IF ( ALLOCATED( A_CM2 ) ) DEALLOCATE( A_CM2 )

  END SUBROUTINE Cleanup_Olson_LandMap
!EOC
END MODULE Olson_LandMap_Mod
