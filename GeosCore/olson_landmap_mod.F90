!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: olson_landmap_mod.F90
!
! !DESCRIPTION: Module OLSON\_LANDMAP\_MOD reads the Olson land map and
!  computes the IREG, ILAND, IUSE, and FRCLND State\_Met arrays. 
!\\
!\\
! !INTERFACE: 
!
MODULE Olson_LandMap_Mod
!
! !USES:
!
  USE CMN_SIZE_MOD                      ! Size parameters
  USE ERROR_MOD                         ! Error checking routines
  USE GC_GRID_MOD                       ! Horizontal grid definition
  USE MAPPING_MOD                       ! Mapping weights & areas
  USE PhysConstants                     ! Physical constants
  USE PRECISION_MOD                     ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
!#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
  PUBLIC  :: Compute_Olson_Landmap_GCHP
!#else
  PUBLIC  :: Init_Olson_Landmap
  PUBLIC  :: Compute_Olson_Landmap
  PUBLIC  :: Cleanup_Olson_LandMap
!#endif
!
! !REMARKS:
!  Eloise Marais and the GEOS-Chem Support Team updated the Olson 2001 
!  landcover dataset and corresponding GEOS-Chem modules in 2012. The Olson
!  2001 landmap superceded the Olson 1992 landmap starting following v9-01-03.
!  The option to use Olson 1992 has been removed from v11-02 and later. 
!  The following text is taken from the data processing README:
!
!  "The Olson 2001 landcover map is at a native resolution of 1km x 1km. I've 
!  identified the dominant vegetation types in each 0.25x0.25 degree gridbox 
!  and use this as input to GEOS-Chem. 
!  
!  The Olson 2001 landcover map also has 96 vegetation types compared with 74 
!  for Olson 1992. For the most part vegetation types 75-96 are either not 
!  dominant vegetation types at 0.25x0.25 degrees or they are crop types that 
!  I lump with other similar vegetation types (either crops or mixed 
!  forest/field vegetation) so that the Olson 2001 landcover dataset at 
!  0.25x0.25 degrees has 74 vegetation types.
!  
!  There are also new vegetation types that are defined in the Olson 2001 
!  dataset from 1-74 that were previously listed as "not used" in 
!  drydep.table. These are assigned appropriate deposition ID # and z0 values. 
!  
!  Vegetation types that are listed as "not used" in the updated drydep.table 
!  dataset are those that are not dominant at 0.25x0.25, but may be present 
!  in the 1kmx1km dataset."
!  
!   The following table shows the the translation between the Olson 2001 and 
!   Olson 1992 land maps:
!
!   Olson 2001				        Olson 1992 	# in
!   LC# Description    		                Equivalent      Dry deposition
!   ==========================================================================
!   1	Urban					1		2
!   2	Low Sparse Grassland			2		3
!   3	Coniferous Forest			3		4
!   4	Deciduous Conifer Forest		4		5
!   5	Deciduous Broadleaf Forest		5		6
!   6	Evergreen Broadleaf Forests		6		7
!   7	Tall Grasses and Shrubs			7		8
!   8	Bare Desert				8		9
!   9	Upland Tundra				9		10
!   10	Irrigated Grassland			10		11
!   11	Semi Desert				11		12
!   12	Glacier Ice				12		13
!   13	Wooded Wet Swamp			13		14
!   14	Inland Water				0		1
!   15	Sea Water				0		1
!   16	Shrub Evergreen				16		17
!   17	Shrub Deciduous				18		19
!   18	Mixed Forest and Field			none present	
!   19	Evergreen Forest and Fields		19		20
!   20	Cool Rain Forest			20		21
!   21	Conifer Boreal Forest			21		22
!   22	Cool Conifer Forest			22		23
!   23	Cool Mixed Forest			23		24
!   24	Mixed Forest				24		25
!   25	Cool Broadleaf Forest			25		26
!   26	Deciduous Broadleaf Forest		26		27
!   27	Conifer Forest				27		28
!   28	Montane Tropical Forests		28		29
!   29	Seasonal Tropical Forest		29		30
!   30	Cool Crops and Towns			30		31
!   31	Crops and Town				31		32
!   32	Dry Tropical Woods			32		33
!   33	Tropical Rainforest			33		34
!   34	Tropical Degraded Forest		34		35
!   35	Corn and Beans Cropland			35		36
!   36	Rice Paddy and Field			36		37
!   37	Hot Irrigated Cropland			37		38
!   38	Cool Irrigated Cropland			38		39
!   39	Cold Irrigated Cropland			none present	
!   40	Cool Grasses and Shrubs			40		41
!   41	Hot and Mild Grasses and Shrubs		41		42
!   42	Cold Grassland				42		43
!   43	Savanna (Woods)				43		44
!   44	Mire, Bog, Fen				44		45
!   45	Marsh Wetland				45		46
!   46	Mediterranean Scrub			46		47
!   47	Dry Woody Scrub				47		48
!   48	Dry Evergreen Woods			none present	
!   49	Volcanic Rock				none present	
!   50	Sand Desert				none present	
!   51	Semi Desert Shrubs			51		52
!   52	Semi Desert Sage			52		53
!   53	Barren Tundra				53		54
!   54	Cool Southern Hemisphere Mixed Forests	54		55
!   55	Cool Fields and Woods			55		56
!   56	Forest and Field			56		57
!   57	Cool Forest and Field			57		58
!   58	Fields and Woody Savanna		58		59
!   59	Succulent and Thorn Scrub		59		60
!   60	Small Leaf Mixed Woods			60		61
!   61	Deciduous and Mixed Boreal Forest	61		62
!   62	Narrow Conifers				62		63
!   63	Wooded Tundra				63		64
!   64	Heath Scrub				64		65
!   65	Coastal Wetland, NW			none present	
!   66	Coastal Wetland, NE			none present	
!   67	Coastal Wetland, SE			none present	
!   68	Coastal Wetland, SW			none present	
!   69	Polar and Alpine Desert			69		70
!   70	Glacier Rock				none present	
!   71	Salt Playas				none present	
!   72	Mangrove				72		73
!   73	Water and Island Fringe			none present	
!   74	Land, Water, and Shore (see Note 1)	none present	
!   75	Land and Water, Rivers (see Note 1)	none present	
!   76	Crop and Water Mixtures			36		37
!   77	Southern Hemisphere Conifers		none present	
!   78	Southern Hemisphere Mixed Forest	32		33
!   79	Wet Sclerophylic Forest			26		27
!   80	Coastline Fringe			none present	
!   81	Beaches and Dunes			none present	
!   82	Sparse Dunes and Ridges			none present	
!   83	Bare Coastal Dunes			none present	
!   84	Residual Dunes and Beaches		none present	
!   85	Compound Coastlines			none present	
!   86	Rocky Cliffs and Slopes			none present	
!   87	Sandy Grassland and Shrubs		none present	
!   88	Bamboo					none present	
!   89	Moist Eucalyptus			26		27
!   90	Rain Green Tropical Forest		33		34
!   91	Woody Savanna				43		44
!   92	Broadleaf Crops				29		30
!   93	Grass Crops				41		42
!   94	Crops, Grass, Shrubs			41		42
!   95	Evergreen Tree Crop			33		34
!   96	Deciduous Tree Crop			33		34
!                                                                             
!  Arrays computed by olson_landmap_mod.F90
!  ============================================================================
!  (1) IREG   (in CMN_DEP_mod.F): # of Olson land types per GC grid box 
!  (2) ILAND  (in CMN_DEP_mod.F): List of all Olson land types in GC grid box
!  (3) IUSE   (in CMN_DEP_mod.F): Coverage of each Olson type in GC grid box
!  (4) FRCLND (in CMN_DEP_mod.F): Fraction of G-C grid box that is not water
!                           
!  The variables are defined as follows:
!      State_Met%IREG(I,J)    : # of land types in horizontal grid cell (I,J)
!      State_Met%ILAND(I,J,T) : Land type ID for land types T=1,IREG(I,J)
!      State_Met%IUSE(I,J,T)  : Fraction area (per mil) occupied by land types
!                               T=1,IREG(I,J) 
!      State_Met%FRCLND(I,J)  : Fraction area occupied by land for cell (I,J)
!                                                  
!  NOTES: 
!  (1) IREG, ILAND, and IUSE are used by the soil NOx emissions routines
!  (2) FRCLND is used by various GEOS-Chem routines
!                                                                             
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
!  20 Mar 2014 - R. Yantosca - Speed up Olson computation by skipping boxes
!  24 Jun 2014 - R. Yantosca - Remove references to logical_mod.F
!  17 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  18 Oct 2016 - E. Lundgren - Add GCHP routine for computing landmap variables
!  02 Nov 2016 - E. Lundgren - Remove N_OLSON since same as global NSURFTYPE
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  13 Sep 2017 - M. Sulprizio- Remove Input_Opt%USE_OLSON_2001. Olson 2001 is
!                              now the default.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Scalars
  INTEGER              :: I_OLSON       ! # of lons (0.5 x 0.5)
  INTEGER              :: J_OLSON       ! # of lats (0.5 x 0.5)
  REAL(fp)             :: D_LON         ! Delta longitude, Olson grid [degrees]
  REAL(fp)             :: D_LAT         ! Delta latitude,  Olson grid [degrees]

  ! Arrays
  REAL*4,  ALLOCATABLE :: lon  (:    )  ! Lon centers, Olson grid [degrees]
  REAL*4,  ALLOCATABLE :: lat  (  :  )  ! Lat centers, Olson grid [degrees]
  INTEGER, ALLOCATABLE :: OLSON(:,:,:)  ! Olson land types
  REAL*4,  ALLOCATABLE :: A_CM2(:,:,:)  ! Surface areas [cm2]

CONTAINS
!EOC
!#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_olson_landmap_gchp
!
! !DESCRIPTION: Subroutine COMPUTE\_OLSON\_LANDMAP\_GCHP computes the 
!  GEOS-Chem State\_Met variables that are dependent on the Olson Landmap, 
!  specifically IREG, ILAND, IUSE, and FRCLND.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Olson_Landmap_GCHP( am_I_Root, State_Met, RC )
!
! !USES:
!
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)    :: am_I_Root    ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),  INTENT(INOUT) :: State_Met    ! Meteorology State object
    INTEGER,         INTENT(INOUT) :: RC
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  27 Sep 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER        :: I, J, T
    INTEGER        :: typeCounter, maxFracInd(1), sumIUSE

    !======================================================================
    ! Initialize
    !======================================================================

    
    ! Loop over all grid cells to set State_Met variables
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Initialize fraction land for this grid cell
       State_Met%FRCLND(I,J) = 1.e+0_fp ! Initialized as all land

       ! Initialize local variables
       typeCounter = 0       ! Tally of number of types found in the cell
       maxFracInd  = 0       ! type index with greatest coverage
       sumIUSE     = 0       ! total coverage across all types [mil]

       ! Loop over all landmap types to set IREG, ILAND, and IUSE
       DO T = 1, NSURFTYPE

          ! If this type has non-zero coverage in this grid box, update vars
          IF ( State_Met%LandTypeFrac(I,J,T) > 0.e+0_fp ) THEN

             ! Increment number of types in this cell
             typeCounter = typeCounter + 1

             ! Set IREG to number of types
             State_Met%IREG(I,J) = typeCounter

             ! Store type index in ILAND array for this grid cell.
             ! Use 0-based index for compatibility with legacy drydep code.
             State_Met%ILAND(I,J,typeCounter) = T-1
             
             ! Store fractional coverage in IUSE array for this grid cell.
             ! Units are [mil] for compatibility with legacy drydep code.
             State_Met%IUSE(I,J,typeCounter) = State_Met%LandTypeFrac(I,J,T) &
                                               * 1000

             ! If this type is water, set fraction land
             IF ( T .eq. 1 ) THEN
                State_Met%FRCLND(I,J) = 1.e+0_fp                          &
                                        - State_Met%LandTypeFrac(I,J,T)
             ENDIF

          ENDIF
       ENDDO

       ! Get IUSE type index with maximum coverage [mil]
       maxFracInd  = MAXLOC(State_Met%IUSE(I,J,1:State_Met%IREG(I,J)))

       ! Force IUSE to sum to 1000 by updating max value if necessary
       sumIUSE =  SUM(State_Met%IUSE(I,J,1:State_Met%IREG(I,J)))
       IF ( sumIUSE /= 1000 ) THEN
          State_Met%IUSE(I,J,maxFracInd) = State_Met%IUSE(I,J,maxFracInd) &
                                           + ( 1000 - sumIUSE )

       ENDIF

    ENDDO
    ENDDO

  END SUBROUTINE Compute_Olson_Landmap_GCHP
!EOC
!#else
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
  SUBROUTINE Init_Olson_LandMap( am_I_Root, Input_Opt, RC )
!
! !USES:
!

    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close
    
    IMPLICIT NONE
    
#   include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
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
!  24 Jun 2014 - R. Yantosca - Now accept Input_Opt, RC via the arg list
!  05 Mar 2015 - R. Yantosca - Now read data w/r/t ExtData/CHEM_INPUTS
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

    ! Assume success 
    RC = GC_SUCCESS

    I_OLSON = 1440                                     ! # lons (0.25x0.25)
    J_OLSON = 720                                      ! # lats (0.25x0.25)
    D_LON   = 0.25e+0_fp                               ! Delta lon [degrees]
    D_LAT   = 0.25e+0_fp                               ! Delta lat [degrees]
    nc_file = 'Olson_2001_Land_Map.025x025.generic.nc' ! Input file name

    ! Allocate arrays
    ALLOCATE( lon  ( I_OLSON             ), STAT=as )  
    ALLOCATE( lat  ( J_OLSON             ), STAT=as )
    ALLOCATE( OLSON( I_OLSON, J_OLSON, 1 ), STAT=as ) 
    ALLOCATE( A_CM2( I_OLSON, J_OLSON, 1 ), STAT=as )

    !======================================================================
    ! Open and read data from the netCDF file
    !======================================================================

    ! Construct file path from directory & file name
    nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // 'Olson_Land_Map_201203/'
    nc_path = TRIM( nc_dir )                    // TRIM( nc_file )

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
    USE State_Met_Mod,      ONLY : MetState
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
!  20 Mar 2014 - R. Yantosca - Add shunts in lat & lon to reduce wall time
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
    INTEGER :: ctOlson (IIPAR, JJPAR, 0:NSURFTYPE-1) ! Count of land types/box
    REAL*4  :: frOlson (IIPAR, JJPAR, 0:NSURFTYPE-1) ! Frac of land types/box
    INTEGER :: ordOlson(IIPAR, JJPAR, 0:NSURFTYPE-1) ! Order of land types

    ! Pointers
    INTEGER,  POINTER :: IREG(:,:)
    INTEGER,  POINTER :: ILAND(:,:,:)
    INTEGER,  POINTER :: IUSE(:,:,:)
    REAL(fp), POINTER :: FRCLND(:,:)
!
! !DEFINED PARAMETERS:
!
! The following parameters are used to skip over Olson NATIVE GRID boxes
! that are too far away from the GEOS-CHEM GRID BOX.  This can speed up
! the Olson computation by a factor of 100 or more!
!
#if defined( GRID05x0625 ) || defined( GRID025x03125 )
    REAL(fp), PARAMETER :: latThresh = 1e+0_fp   ! Lat threshold, nested grid
    REAL(fp), PARAMETER :: lonThresh = 1e+0_fp   ! Lon threshold, nested grid
#else
    REAL(fp), PARAMETER :: latThresh = 5e+0_fp   ! Lat threshold, global
    REAL(fp), PARAMETER :: lonThresh = 6e+0_fp   ! Lon threshold, global
#endif

    !======================================================================
    ! NATIVE GRID parameters (i.e. 0.5 x 0.5 "GENERIC")
    !======================================================================

    ! Be lazy, construct lon edges from lon centers
    DO I = 1, I_OLSON
       lonedge(I)      = DBLE( lon(I) ) - ( D_LON * 0.5e+0_fp )
       indLon(I)       = I
    ENDDO
    lonedge(I_OLSON+1) = lonedge(I_OLSON) + D_LON
    
    ! Be lazy, construct lat edges from lat centers
    DO J = 1, J_OLSON
       latedge(J)      = DBLE( lat(J) ) - ( D_LAT * 0.5e+0_fp )
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
    !$OMP PRIVATE( sumIUse,  uniqOlson, C,        IG                 ) &
    !$OMP SCHEDULE( DYNAMIC )
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
       ! Find each NATIVE GRID BOX that fits into the GEOS-CHEM GRID BOX.  
       ! Keep track of the land types and coverage fractions.
       !===================================================================

       ! Loop over latitudes on the NATIVE GRID
       DO JJ  = 1, J_OLSON

          ! Latitude edges of this NATIVE GRID box
          yedge_s    = latedge(JJ  )                ! S edge
          yedge_n    = latedge(JJ+1)                ! N edge

          !%%%%%% LATITUDE SHUNT TO REDUCE WALL TIME (bmy, 3/20/14) %%%%%%%%%%
          !%%%
          !%%% Skip further computations unless we are within LATTHRESH 
          !%%% degrees of the western edge of box (I,J).  This prevents 
          !%%% excess computations and subroutine calls.
          !%%%
          IF ( ABS( yedge_s - yedgeC_s ) > latThresh ) CYCLE
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! Loop over longitudes on the NATIVE GRID
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
             xedge_e    = lonedge(II+1)                ! E edge

             ! Because the first GEOS-CHEM GRID BOX straddles the date line,
             ! we have to adjust the W and E edges of the NATIVE GRID BOX to
             ! be in monotonically increasing order.  This will prevent
             ! erronous results from being returned by GET_MAP_WT below.
             IF ( isGlobal .and. IG == 1 .and. II >= shiftLon(1) )  THEN
                xedge_w = xedge_w - 360e0
                xedge_e = xedge_e - 360e0
             ENDIF
         
             !%%%%%% LONGITUDE SHUNT TO REDUCE WALL TIME (bmy, 3/20/14) %%%%%%
             !%%%
             !%%% Skip further computations unless we are within LONTHRESH 
             !%%% degrees of the western edge of box (I,J).  This prevents 
             !%%% excess computations and subroutine calls.
             !%%%
             IF ( ABS( xedge_w - xedgeC_w ) > lonThresh ) CYCLE
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             ! "Area" of the NATIVE GRID BOX in degrees (DLON * DLAT)
             dxdy       = ( xedge_e - xedge_w ) * ( yedge_n - yedge_s )

             ! Get the mapping weight (i.e. The fraction of the NATIVE 
             ! GRID BOX that lies w/in the GEOS-CHEM GRID BOX)
             CALL GET_MAP_WT( xedge_w, xedge_e, xedgeC_w, xedgeC_e,  &
                              yedge_s, yedge_n, yedgeC_s, yedgeC_n,  &
                              mapWt                                 )

             ! Skip unless part (or all) of the NATIVE GRID BOX
             ! actually fits into the GEOS-CHEM GRID BOX
             IF ( mapWt <= 0e0 .or. mapWt > 1e0 ) CYCLE

             ! Area of the NATIVE GRID BOX that lies
             ! within the GEOS-CHEM GRID BOX
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
       DO T = 0, NSURFTYPE-1

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
!###    C = mapping(23,34)%count
!###    print*, '### count   : ', C
!###    print*, '### II      : ', mapping(23,34)%II(1:C)
!###    print*, '### JJ      : ', mapping(23,34)%JJ(1:C)
!###    print*, '### area    : ', mapping(23,34)%area(1:C)
!###    print*, '### sumarea : ', mapping(23,34)%sumarea

  END SUBROUTINE Compute_Olson_LandMap
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
!#endif
END MODULE Olson_LandMap_Mod
