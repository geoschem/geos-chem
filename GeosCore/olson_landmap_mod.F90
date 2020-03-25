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
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_LandTypeFrac
  PUBLIC  :: Compute_Olson_Landmap
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
!  The following table shows the the translation between the Olson 2001 and
!  Olson 1992 land maps:
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
!  (1) State_Met%IREG(I,J)    : # of land types in horizontal grid cell (I,J)
!  (2) State_Met%ILAND(I,J,T) : Land type ID for land types T=1,IREG(I,J)
!  (3) State_Met%IUSE(I,J,T)  : Fraction area (per mil) occupied by land types
!                               T=1,IREG(I,J)
!  (4) State_Met%FRCLND(I,J)  : Fraction area occupied by land for cell (I,J)
!
! !REVISION HISTORY:
!  13 Mar 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: compute_olson_landmap
!
! !DESCRIPTION: Subroutine COMPUTE\_OLSON\_LANDMAP computes the
!  GEOS-Chem State\_Met variables that are dependent on the Olson Landmap,
!  specifically IREG, ILAND, IUSE, and FRCLND.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Olson_Landmap( Input_Opt, State_Grid, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,   ONLY : NSURFTYPE
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),  INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState),  INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),  INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  27 Sep 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER   :: I,           J,             T
    INTEGER   :: typeCounter, maxFracInd(1), sumIUSE

    ! Arrays

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Compute_Olson_Landmap begins here!
    !=======================================================================

    ! Initialize
    RC     = GC_SUCCESS
    ErrMsg = ''
    ThisLoc = &
     '-> at Compute_Olson_Landmap (in module GeosCore/olson_landmap_mod.F90)'

    !-----------------------------------------------------------------------
    ! Loop over all grid cells to set State_Met variables
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize fraction land for this grid cell
       State_Met%FRCLND(I,J) = 1.e+0_fp ! Initialized as all land

       ! Initialize local variables
       typeCounter = 0       ! Tally of number of types found in the cell
       maxFracInd  = 0       ! type index with greatest coverage
       sumIUSE     = 0       ! total coverage across all types [mil]

       !--------------------------------------------------------------------
       ! Loop over all landmap types to set IREG, ILAND, and IUSE
       !--------------------------------------------------------------------
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
                State_Met%FRCLND(I,J) = 1.e+0_fp                             &
                                        - State_Met%LandTypeFrac(I,J,T)
             ENDIF

          ENDIF
       ENDDO

       !--------------------------------------------------------------------
       ! Make sure that State_Met%IUSE sums up to 1000 (per mil)
       !--------------------------------------------------------------------

       ! Get IUSE type index with maximum coverage [mil]
       ! (NOTE: MAXLOC returns a vector with 1 element)
       maxFracInd = MAXLOC( State_Met%IUSE(I,J,1:State_Met%IREG(I,J)) )

       ! Make sure we find the index of IUSE with maximum coverage
       IF ( maxFracInd(1) > 0 ) THEN

          ! Force IUSE to sum to 1000 by updating max value if necessary
          sumIUSE =  SUM(State_Met%IUSE(I,J,1:State_Met%IREG(I,J)))
          IF ( sumIUSE /= 1000 ) THEN
             State_Met%IUSE(I,J,maxFracInd(1)) =                             &
             State_Met%IUSE(I,J,maxFracInd(1)) + ( 1000 - sumIUSE )
          ENDIF

       ELSE

          ! If we could not find the index IUSE with maximum coverage,
          ! then this indicates a potential problem with the regridding.
          ! Throw an error and exit the routine here.
          WRITE( ErrMsg, 100 ) I, J
 100      FORMAT( 'Error: State_Met%IUSE is zero at grid box (',             &
                  i6, ', ', i6,                                              &
                  '!  This indicates a potential regridding problem! ' )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN

       ENDIF

    ENDDO
    ENDDO

  END SUBROUTINE Compute_Olson_Landmap
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landtypefrac
!
! !DESCRIPTION: Attaches pointers from the MODIS XLAI data read in by
!  HEMCO to the LandTypeFrac field of State\_Met.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_LandTypeFrac( Input_Opt, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,      ONLY : NSURFTYPE
    USE ErrCode_Mod
    USE Hco_Interface_Mod, ONLY : HcoState
    USE Hco_EmisList_Mod,  ONLY : Hco_GetPtr
    USE Input_Opt_Mod,     ONLY : OptInput
    USE State_Met_Mod,     ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This follows the same methodology for GCHP, except that GCHP obtains the
!  land type fractions via MAPL, and here we obtain them via HEMCO.
!
! !REVISION HISTORY:
!  13 Feb 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: T

    ! Strings
    CHARACTER(LEN=10)  :: Name
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Pointers
    REAL(f4), POINTER  :: Ptr2D(:,:)

    !=======================================================================
    ! Init_LandTypeFrac begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
      ' -> at Init_LandTypeFrac (in module "GeosCore/olson_landmap_mod.F90'

    ! Free pointer
    Ptr2D => NULL()

    ! Loop over the number of Olson land types
    DO T = 1, NSURFTYPE

       ! Get the HEMCO pointer to each Olson landtype mask
       ! (variable names are LANDTYPE00, LANDTYPE01 .. LANDTYPE72)
       WRITE( Name, 100 ) T-1
 100   FORMAT( 'LANDTYPE', i2.2 )
       CALL HCO_GetPtr( HcoState, Name, Ptr2D, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not get pointer to HEMCO field: ' // TRIM( Name )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Met%LandTypeFrac
       State_Met%LandTypeFrac(:,:,T) = Ptr2D

       ! Free pointer
       Ptr2D => NULL()

    ENDDO

  END SUBROUTINE Init_LandTypeFrac
!EOC
END MODULE Olson_LandMap_Mod
