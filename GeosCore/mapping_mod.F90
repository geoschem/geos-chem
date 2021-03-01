!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mapping_mod.F90
!
! !DESCRIPTION: Module MAPPING\_MOD contains a derived-type object to compute
!  and save the mapping weight (i.e. fraction of each "fine" grid box that
!  fits into the "coarse" grid box") and areal mapping (i.e. the area of each
!  "fine" grid box contained within a "coarse" grid box).
!\\
!\\
! !INTERFACE:
!
MODULE Mapping_Mod
!
! !USES:
!
  USE ERROR_MOD                       ! Error handling routines
  USE PRECISION_MOD                   ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
  PUBLIC :: MapWeight
  TYPE MapWeight
     INTEGER          :: count        ! # of "fine" boxes per "coarse" box
     INTEGER, POINTER :: II(:)        ! Longitude indices,  "fine"   grid
     INTEGER, POINTER :: JJ(:)        ! Latitude  indices,  "fine"   grid
     INTEGER, POINTER :: olson(:)     ! Olson land type,    "fine"   grid
     INTEGER, POINTER :: ordOlson(:)  ! Ordering of Olson land types
     REAL*4,  POINTER :: area(:)      ! Surface areas,      "fine"   grid
     REAL*4           :: sumarea      ! Total surface area, "coarse" grid
  END TYPE MapWeight

!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_Mapping
  PUBLIC  :: Get_Map_Wt
  PUBLIC  :: Cleanup_Mapping
!
! !REMARKS:
!  The mapping weights and areal mapping are initialized when the Olson
!  land map is read from disk (in olson_landmap_mod.F90).  They are used
!  again when the MODIS leaf area index data is prepared for input into
!  GEOS-Chem's (legacy) dry deposition module.
!                                                                             .
!  Also, we do not define the mapping weight object within this module.
!  This allows you to create more than one mapping weight object for
!  different native grids (e.g. 0.5 x 0.5 and 0.25 x 0.25, etc.)
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
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
! !IROUTINE: init_mapping
!
! !DESCRIPTION: Subroutine INIT\_MAPPING allocates and initializes a
!  derived-type object containing grid mapping information.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Mapping( Input_Opt, I_FINE, J_FINE, I_COARSE,  J_COARSE, &
                           mapping,   RC  )
!
! !USES:
!
    USE CMN_SIZE_Mod,       ONLY : NSURFTYPE
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt ! Input Options object
    INTEGER,        INTENT(IN) :: I_FINE    ! # of lons on the "fine" grid
    INTEGER,        INTENT(IN) :: J_FINE    ! # of lats on the,"fine" grid
    INTEGER,        INTENT(IN) :: I_COARSE  ! # of lons on the "coarse" grid
    INTEGER,        INTENT(IN) :: J_COARSE  ! # of lats on the "coarse" grid
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MapWeight), POINTER, INTENT(INOUT) :: mapping(:,:) !"fine" -> "coarse"
!
! !OUTPUT PARAMETERS:
!
      INTEGER,       INTENT(OUT) :: RC      ! Success or failure?
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: I,   J,   FINE_PER_COARSE, ADD, as
    INTEGER :: as1, as2, as3,             as4, as5

    ! Arrays
    INTEGER :: err(I_COARSE,J_COARSE)

    !======================================================================
    ! INIT_MAPPING begins here!
    !======================================================================

    ! Initialize
    err = 0
    RC = GC_SUCCESS

    ! Define a number of extra boxes to add to FINE_PER_COARSE
    ! in order to prevent out-of-bounds errors
    ADD = 10

    ! Number of "fine" grid boxes that fit into the "coarse" grid box
    ! (roughly, use the +ADD to make sure it is big enough)
    FINE_PER_COARSE = ( DBLE( I_FINE ) / DBLE( I_COARSE ) ) &
                    * ( DBLE( J_FINE ) / DBLE( J_COARSE ) ) + ADD

    ! For saving mapping weights
    IF ( .not. ASSOCIATED( mapping ) ) THEN

       ! Allocate the mapping weight object
       ALLOCATE( mapping( I_COARSE, J_COARSE ), STAT=as )
       IF ( as /= 0 ) CALL ALLOC_ERR( 'map' )

       ! Populate the mapping weight object
       !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J )
       DO J = 1, J_COARSE
       DO I = 1, I_COARSE

          ! Allocate sub-fields of MAPPING object
          ALLOCATE( mapping(I,J)%ii      ( FINE_PER_COARSE ), STAT=as1 )
          ALLOCATE( mapping(I,J)%jj      ( FINE_PER_COARSE ), STAT=as2 )
          ALLOCATE( mapping(I,J)%olson   ( FINE_PER_COARSE ), STAT=as3 )
          ALLOCATE( mapping(I,J)%ordOlson( 0:NSURFTYPE-1   ), STAT=as4 )
          ALLOCATE( mapping(I,J)%area    ( FINE_PER_COARSE ), STAT=as5 )

          ! Check for allocation error in such a way
          ! as to not interfere w/ the parallel loop
          IF ( as1 + as2 + as3 + as4 + as5 /= 0 ) THEN
             err(I,J) = 1
             EXIT
          ENDIF

          ! Initialize sub-fields
          IF ( err(I,J) == 0 ) THEN
             mapping(I,J)%count    = 0
             mapping(I,J)%ii       = 0
             mapping(I,J)%jj       = 0
             mapping(I,J)%olson    = 0
             mapping(I,J)%ordOlson = 0
             mapping(I,J)%area     = 0e0
             mapping(I,J)%sumarea  = 0e0
          ENDIF
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! Stop if any sub-field of MAPPING was not allocated properly
    IF ( SUM( err ) > 0 ) THEN
       CALL ERROR_STOP( 'Error allocating sub-fields of MAPPING object!', &
                        'Init_Mapping (mapping_mod.F90)' )
    ENDIF

    ! Stop w/ error if the MAPPING object is not dimensioned properly
    IF ( SIZE( mapping, 1 ) /= I_COARSE  .or. &
         SIZE( mapping, 2 ) /= J_COARSE ) THEN
       CALL ERROR_STOP( 'MAPPING object has incorrect dimensions!', &
                        'Init_Mapping (mapping_mod.F90)' )
    ENDIF

  END SUBROUTINE Init_Mapping
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_map_wt
!
! !DESCRIPTION: Subroutine GET\_MAP\_Wt returns the "mapping weight", that
!  is, the fraction that each "fine" grid box fits into each "coarse" grid
!  box.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Map_Wt( xedge_w, xedge_e, xedgeC_w, xedgeC_e,   &
                         yedge_s, yedge_n, yedgeC_s, yedgeC_n,   &
                         mapWt                                  )
!
! !INPUT PARAMETERS:
!
    REAL*4, INTENT(IN)  :: xedge_w,  xedge_e    ! Lon edges, fine grid
    REAL*4, INTENT(IN)  :: xedgeC_w, xedgeC_e   ! Lon edges, coarse grid
    REAL*4, INTENT(IN)  :: yedge_s,  yedge_n    ! Lat edges, fine grid
    REAL*4, INTENT(IN)  :: yedgeC_s, yedgeC_n   ! Lat edges, coarse grid
    REAL*4, INTENT(OUT) :: mapWt                ! Mapping weight
!
! !REMARKS:
!  Follows the algorithm from GAMAP routine ctm_getweight.pro
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: ox1, ox2, nx1, nx2, ov1, ov2, xOverLap
    REAL(fp) :: oy1, oy2, ny1, ny2,           yOverLap

    !======================================================================
    ! Get overlap in longitude
    !======================================================================

    ! OX1, OX2 are the lon edges of the "fine" grid box
    ox1 = xedge_w
    ox2 = xedge_e

    ! NX1, NX2 are the lon edges of the coarse grid box
    nx1 = xedgeC_w
    nx2 = xedgeC_e

    ! Deal with over-the-dateline cases (phs, 9/26/07)
    ! That fixes a problem when going from GEOS-5
    ! 0.66667 x 0.5 to GENERIC 1 x 1.
    ! Maybe it fixes also the kludges below ?? ## need checking
    if ( ox2 .lt. nx1 ) then
       ox1 = ox1 + 360e0
       ox2 = ox2 + 360e0
    endif

    if ( ox1 .gt. nx2 ) then
       ox1 = ox1 - 360e0
       ox2 = ox2 - 360e0
    endif

    ! convert to equivalent longitudes where necessary
    if ( nx1 < -90. .AND. ox1 > 0. ) nx1 = nx1 + 360e0
    if ( nx2 < -90. .AND. ox2 > 0. ) nx2 = nx2 + 360e0

    ! OV1 is the greater of OX1 and NX1
    ! OV2 is the lesser of OX2 and NX2
    ov1 = MAX( ox1, nx1 )
    ov2 = MIN( nx2, ox2 )

    ! XOVERLAP is the fraction of the old (fine) grid box that
    ! occupies the new (coarse) grid box in the longitude
    xOverLap = ( ov2 - ov1 ) / ( ox2 - ox1 )

    ! If XOVERLAP is not in the range of 0-1, then it means that the "fine"
    ! grid box lies completely outside the "coarse" grid box (in longitude).
    ! Set to zero to avoid erroneous results in the calling routine.
    if ( xOverLap < 0e0 .or. xOverLap > 1e0 ) xOverlap = 0e0

    !======================================================================
    ! Get overlap in latitude
    !======================================================================

    ! OY1 and OY2 are lat edges of the "fine" grid
    oy1 = yedge_s
    oy2 = yedge_n

    ! NY1 and NY2 are consecutive Y-edges for the coarse
    ny1 = yedgeC_s
    ny2 = yedgeC_n

    ! OV1 is the greater of OY1 and NY1
    ! OV2 is the lesser of OY2 and NY2
    ov1 = MAX( oy1, ny1 )
    ov2 = MIN( ny2, oy2 )

    ! YOVERLAP is the fraction of the old (fine) grid box that
    ! occupies the new (coarse) grid box in latitude
    yoverlap = ( ov2 - ov1 ) / ( oy2 - oy1 )

    ! If YOVERLAP is not in the range of 0-1, then it means that the "fine"
    ! grid box lies completely outside the "coarse" grid box (in latitude).
    ! Set to zero to avoid erroneous results in the calling routine.
    if ( yOverLap < 0e0 .or. yOverLap > 1e0 ) yOverlap = 0e0

    ! Resultant mapping weight
    mapWt = xOverLap * yOverLap

  END SUBROUTINE Get_Map_Wt
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_mapping
!
! !DESCRIPTION: Subroutine CLEANUP\_MAPPING deallocates memory from a
!  derived-type object containing mapping information.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Mapping( mapping )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MapWeight), POINTER, INTENT(INOUT) :: mapping(:,:)
!
! !REVISION HISTORY:'
!  03 Mar 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER :: I, J

    ! Test if MAP has been allocated memory
    IF ( ASSOCIATED( mapping ) ) THEN

       ! First deallocate the pointer fields of the MAP object
       !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J )
       DO J = 1, SIZE( mapping, 2 )
       DO I = 1, SIZE( mapping, 1 )
          DEALLOCATE( mapping(I,J)%ii       )
          DEALLOCATE( mapping(I,J)%jj       )
          DEALLOCATE( mapping(I,J)%olson    )
          DEALLOCATE( mapping(I,J)%ordOlson )
          DEALLOCATE( mapping(I,J)%area     )
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Then deallocate the MAP object itself
       DEALLOCATE( mapping )
    ENDIF

  END SUBROUTINE Cleanup_Mapping
!EOC
END MODULE Mapping_Mod
