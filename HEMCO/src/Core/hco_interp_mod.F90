!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_interp_mod.F90
!
! !DESCRIPTION: Module HCO\_INTERP\_MOD contains routines to interpolate
! input data onto the HEMCO grid. This module contains routine for
! horizontal regridding between regular grids (MAP\_A2A), as well as
! vertical interpolation amongst GEOS model levels (full <--> reduced).
!\\
!\\
! Regridding is supported for concentration quantities (default) and
! index-based values. For the latter, the values in the regridded grid
! boxes correspond to the value of the original grid that contrbutes most
! to the given box.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_Interp_Mod
!
! !USES:
!
  USE HCO_Types_Mod
  USE HCO_Error_Mod
  USE HCO_State_Mod,       ONLY : Hco_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: ModelLev_Check
  PUBLIC  :: ModelLev_Interpolate
  PUBLIC  :: REGRID_MAPA2A
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PRIVATE :: GEOS5_TO_GEOS4_LOWLEV
  PRIVATE :: COLLAPSE
  PRIVATE :: INFLATE
!
! !REVISION HISTORY:
!  30 Dec 2014 - C. Keller - Initialization
!  03 Feb 2015 - C. Keller   - Added REGRID_MAPA2A (from hcoio_dataread_mod.F90).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE VARIABLES:
!
  ! AP parameter of native GEOS-5 grid. Needed to remap GEOS-5 data from native
  ! onto the reduced vertical grid.
  REAL(hp), TARGET :: G5_EDGE_NATIVE(73) = (/                          &
              0.000000e+00_hp, 4.804826e-02_hp, 6.593752e+00_hp, 1.313480e+01_hp, &
              1.961311e+01_hp, 2.609201e+01_hp, 3.257081e+01_hp, 3.898201e+01_hp, &
              4.533901e+01_hp, 5.169611e+01_hp, 5.805321e+01_hp, 6.436264e+01_hp, &
              7.062198e+01_hp, 7.883422e+01_hp, 8.909992e+01_hp, 9.936521e+01_hp, &
              1.091817e+02_hp, 1.189586e+02_hp, 1.286959e+02_hp, 1.429100e+02_hp, &
              1.562600e+02_hp, 1.696090e+02_hp, 1.816190e+02_hp, 1.930970e+02_hp, &
              2.032590e+02_hp, 2.121500e+02_hp, 2.187760e+02_hp, 2.238980e+02_hp, &
              2.243630e+02_hp, 2.168650e+02_hp, 2.011920e+02_hp, 1.769300e+02_hp, &
              1.503930e+02_hp, 1.278370e+02_hp, 1.086630e+02_hp, 9.236572e+01_hp, &
              7.851231e+01_hp, 6.660341e+01_hp, 5.638791e+01_hp, 4.764391e+01_hp, &
              4.017541e+01_hp, 3.381001e+01_hp, 2.836781e+01_hp, 2.373041e+01_hp, &
              1.979160e+01_hp, 1.645710e+01_hp, 1.364340e+01_hp, 1.127690e+01_hp, &
              9.292942e+00_hp, 7.619842e+00_hp, 6.216801e+00_hp, 5.046801e+00_hp, &
              4.076571e+00_hp, 3.276431e+00_hp, 2.620211e+00_hp, 2.084970e+00_hp, &
              1.650790e+00_hp, 1.300510e+00_hp, 1.019440e+00_hp, 7.951341e-01_hp, &
              6.167791e-01_hp, 4.758061e-01_hp, 3.650411e-01_hp, 2.785261e-01_hp, &
              2.113490e-01_hp, 1.594950e-01_hp, 1.197030e-01_hp, 8.934502e-02_hp, &
              6.600001e-02_hp, 4.758501e-02_hp, 3.270000e-02_hp, 2.000000e-02_hp, &
              1.000000e-02_hp /)

  ! AP parameter of native GEOS-4 grid. Needed to remap GEOS-4 data from native
  ! onto the reduced vertical grid.
  REAL(hp), TARGET :: G4_EDGE_NATIVE(56) = (/       &
                    0.000000_hp,   0.000000_hp,  12.704939_hp, &
                   35.465965_hp,  66.098427_hp, 101.671654_hp, &
                  138.744400_hp, 173.403183_hp, 198.737839_hp, &
                  215.417526_hp, 223.884689_hp, 224.362869_hp, &
                  216.864929_hp, 201.192093_hp, 176.929993_hp, &
                  150.393005_hp, 127.837006_hp, 108.663429_hp, &
                   92.365662_hp,  78.512299_hp,  66.603378_hp, &
                   56.387939_hp,  47.643932_hp,  40.175419_hp, &
                   33.809956_hp,  28.367815_hp,  23.730362_hp, &
                   19.791553_hp,  16.457071_hp,  13.643393_hp, &
                   11.276889_hp,   9.292943_hp,   7.619839_hp, &
                    6.216800_hp,   5.046805_hp,   4.076567_hp, &
                    3.276433_hp,   2.620212_hp,   2.084972_hp, &
                    1.650792_hp,   1.300508_hp,   1.019442_hp, &
                    0.795134_hp,   0.616779_hp,   0.475806_hp, &
                    0.365041_hp,   0.278526_hp,   0.211349_hp, &
                    0.159495_hp,   0.119703_hp,   0.089345_hp, &
                    0.066000_hp,   0.047585_hp,   0.032700_hp, &
                    0.020000_hp,   0.010000_hp /)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Regrid_MAPA2A
!
! !DESCRIPTION: Subroutine Regrid\_MAPA2A regrids input array NcArr onto
! the simulation grid and stores the data in list container Lct. Horizontal
! regridding is performed using MAP\_A2A algorithm. Vertical interpolation
! between GEOS levels (full vs. reduced, GEOS-5 vs. GEOS-4), is also
! supported.
!\\
!\\
! This routine can remap concentrations and index-based quantities.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE REGRID_MAPA2A( HcoState, NcArr, LonE, LatE, Lct, RC )
!
! !USES:
!
    USE REGRID_A2A_Mod,     ONLY : MAP_A2A
    USE HCO_FileData_Mod,   ONLY : FileData_ArrCheck
    USE HCO_UNIT_MOD,       ONLY : HCO_IsIndexData
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState          ! HEMCO state object
    REAL(sp),         POINTER        :: NcArr(:,:,:,:)    ! 4D input data
    REAL(hp),         POINTER        :: LonE(:)           ! Input grid longitude edges
    REAL(hp),         POINTER        :: LatE(:)           ! Input grid latitude edges
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER        :: Lct               ! HEMCO list container
    INTEGER,          INTENT(INOUT)  :: RC                ! Success or failure?
!
! !REVISION HISTORY:
!  03 Feb 2015 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: nLonEdge, nLatEdge
    INTEGER                 :: NX, NY, NZ, NLEV, NTIME, NCELLS
    INTEGER                 :: I, J, L, T, AS, I2
    INTEGER                 :: nIndex
    REAL(sp), ALLOCATABLE   :: LonEdgeI(:)
    REAL(sp), ALLOCATABLE   :: LatEdgeI(:)
    REAL(sp)                :: LonEdgeO(HcoState%NX+1)
    REAL(sp)                :: LatEdgeO(HcoState%NY+1)

    REAL(sp), POINTER       :: ORIG_2D(:,:)
    REAL(sp), POINTER       :: REGR_2D(:,:)
    REAL(sp), POINTER       :: REGR_4D(:,:,:,:)

    REAL(sp), ALLOCATABLE, TARGET :: FRACS(:,:,:,:)
    REAL(hp), ALLOCATABLE         :: REGFRACS(:,:,:,:)
    REAL(hp), ALLOCATABLE         :: MAXFRACS(:,:,:,:)
    REAL(hp), ALLOCATABLE         :: INDECES(:,:,:,:)
    REAL(hp), ALLOCATABLE         :: UNIQVALS(:)
    REAL(hp)                      :: IVAL
    LOGICAL                       :: IsIndex

    LOGICAL                 :: VERB
    CHARACTER(LEN=255)      :: MSG
    CHARACTER(LEN=255)      :: LOC = 'ModelLev_Interpolate (hco_interp_mod.F90)'

    !=================================================================
    ! REGRID_MAPA2A begins here
    !=================================================================

    ! Init
    ORIG_2D => NULL()
    REGR_2D => NULL()
    REGR_4D => NULL()

    ! Check for verbose mode
    verb = HCO_IsVerb(HcoState%Config%Err,  3 )

    ! get longitude / latitude sizes
    nLonEdge = SIZE(LonE,1)
    nLatEdge = SIZE(LatE,1)

    ! Write input grid edges to shadow variables so that map_a2a accepts them
    ! as argument.
    ! Also, for map_a2a, latitudes have to be sines...
    ALLOCATE(LonEdgeI(nlonEdge), LatEdgeI(nlatEdge), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'alloc error LonEdgeI/LatEdgeI', RC, THISLOC=LOC )
       RETURN
    ENDIF
    LonEdgeI(:) = LonE
    LatEdgeI(:) = SIN( LatE * HcoState%Phys%PI_180 )

    ! Get output grid edges from HEMCO state
    LonEdgeO(:) = HcoState%Grid%XEDGE%Val(:,1)
    LatEdgeO(:) = HcoState%Grid%YSIN%Val(1,:)

    ! Get input array sizes
    NX     = size(ncArr,1)
    NY     = size(ncArr,2)
    NLEV   = size(ncArr,3)
    NTIME  = size(ncArr,4)
    NCELLS = NX * NY * NLEV * NTIME

    ! Are these index-based data? If so, need to remap the fraction (1 or 0)
    ! of every value independently. For every grid box, the value with the
    ! highest overlap (closest to 1) is taken.
    IsIndex = HCO_IsIndexData(Lct%Dct%Dta%OrigUnit)

    IF ( IsIndex ) THEN

       ! Allocate working arrays:
       ! - FRACS contains the fractions on the original grid. These are
       !   binary (1 or 0).
       ! - MAXFRACS stores the highest used fraction for each output grid
       !   box. Will be updated continously.
       ! - INDECES is the output array holding the index-based remapped
       !   values. Will be updated continuously.
       ! - UNIQVALS is a vector holding all unique values of the input
       !   array (NINDEX is the number of unique values).
       !
       ! ckeller, 9/24/15: Extend vertical axis of MAXFRACS, REGFRACS, and
       ! INDECES to HcoState%NZ+1 for fields that are on edges instead of
       ! mid-points.
       ALLOCATE( FRACS(NX,NY,NLEV,NTIME), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( 'alloc error FRACS', RC, THISLOC=LOC )
          RETURN
       ENDIF
       ALLOCATE( MAXFRACS(HcoState%NX,HcoState%NY,HcoState%NZ+1,NTIME), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( 'alloc error MAXFRACS', RC, THISLOC=LOC )
          RETURN
       ENDIF
       ALLOCATE( REGFRACS(HcoState%NX,HcoState%NY,HcoState%NZ+1,NTIME), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( 'alloc error INDECES', RC, THISLOC=LOC )
          RETURN
       ENDIF
       ALLOCATE( INDECES(HcoState%NX,HcoState%NY,HcoState%NZ+1,NTIME), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( 'alloc error INDECES', RC, THISLOC=LOC )
          RETURN
       ENDIF
       ALLOCATE( UNIQVALS(NCELLS), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( 'alloc error INDECES', RC, THISLOC=LOC )
          RETURN
       ENDIF
       FRACS    = 0.0_sp
       REGFRACS = 0.0_hp
       MAXFRACS = 0.0_hp
       INDECES  = 0.0_hp
       UNIQVALS = 0.0_hp

       ! Get unique values. Loop over all input data values and add
       ! them to UNIQVALS vector if UNIQVALS doesn't hold that same value
       ! yet.
       NINDEX = 0
       DO T = 1, NTIME
       DO L = 1, NLEV
       DO J = 1, NY
       DO I = 1, NX

          ! Current value
          IVAL = NcArr(I,J,L,T)

          ! Check if value already exists in UNIQVALS
          IF ( NINDEX > 0 ) THEN
             IF ( ANY(UNIQVALS(1:NINDEX) == IVAL) ) CYCLE
          ENDIF

          ! Add to UNIQVALS
          NINDEX = NINDEX + 1
          UNIQVALS(NINDEX) = IVAL
       ENDDO
       ENDDO
       ENDDO
       ENDDO

       ! Verbose mode
       IF ( verb ) THEN
          MSG = 'Do index based regridding for field ' // TRIM(Lct%Dct%cName)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,*) '   - Number of indeces: ', NINDEX
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

    ELSE
       NINDEX = 1
    ENDIF

    ! Define array to put horizontally regridded data onto. If this
    ! is 3D data, we first regrid all vertical levels horizontally
    ! and then pass these data to the list container. In this second
    ! step, levels may be deflated/collapsed.

    ! 2D data is directly passed to the data container
    IF ( Lct%Dct%Dta%SpaceDim <= 2 ) THEN
       CALL FileData_ArrCheck( HcoState%Config, Lct%Dct%Dta, &
                               HcoState%NX, HcoState%NY, NTIME, RC )
       IF ( RC /= 0 ) RETURN
    ENDIF

    ! 3D data and index data is first written into a temporary array,
    ! REGR_4D.
    IF ( Lct%Dct%Dta%SpaceDim == 3 .OR. IsIndex ) THEN
       ALLOCATE( REGR_4D(HcoState%NX,HcoState%NY,NLEV,NTIME), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( 'alloc error REGR_4D', RC, THISLOC=LOC )
          RETURN
       ENDIF
       REGR_4D = 0.0_hp
    ENDIF

    ! Do regridding for every index value. If it's not index data, this loop
    ! is executed only once (NINDEX=1).
    DO I = 1, NINDEX

       ! For index based data, create fractions array for the given index.
       IF ( IsIndex ) THEN
          IVAL = UNIQVALS(I)
          WHERE( ncArr == IVAL )
             FRACS = 1.0_sp
          ELSEWHERE
             FRACS = 0.0_sp
          END WHERE
       ENDIF

       ! Regrid horizontally
       DO T = 1, NTIME
       DO L = 1, NLEV

          ! Point to 2D slices to be regridded:
          ! - Original 2D array
          IF ( IsIndex ) THEN
            ORIG_2D => FRACS(:,:,L,T)
          ELSE
            ORIG_2D => ncArr(:,:,L,T)
          ENDIF

          ! - Regridded 2D array
          IF ( Lct%Dct%Dta%SpaceDim <= 2 .AND. .NOT. IsIndex ) THEN
             REGR_2D => Lct%Dct%Dta%V2(T)%Val(:,:)
          ELSE
             REGR_2D => REGR_4D(:,:,L,T)
          ENDIF

          ! Do the regridding
          CALL MAP_A2A( NX,      NY, LonEdgeI,    LatEdgeI, ORIG_2D,  &
                        HcoState%NX, HcoState%NY, LonEdgeO, LatEdgeO, &
                        REGR_2D, 0, 0, HCO_MISSVAL )
          ORIG_2D => NULL()
          REGR_2D => NULL()

       ENDDO !L
       ENDDO !T

       ! Eventually inflate/collapse levels onto simulation levels.
       IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
          CALL ModelLev_Interpolate( HcoState, REGR_4D, Lct, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! For index based data, map fractions back to corresponding value.
       ! Array INDECES holds the index-based remapped values. Set INDECES
       ! to current index value in every grid box where the regridded
       ! fraction of this index is higher than any previous fraction
       ! (array MAXFRACS stores the highest used fraction in each grid box).
       IF ( IsIndex ) THEN

          ! Reset
          REGFRACS = 0.0_hp

          ! 3D data written to Lct needs to be mapped back onto REGR_4D.
          IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
             DO T = 1, NTIME
                NZ = SIZE(Lct%Dct%Dta%V3(T)%Val,3)
                REGFRACS(:,:,1:NZ,T) = Lct%Dct%Dta%V3(T)%Val(:,:,:)
             ENDDO
          ELSE
             REGFRACS(:,:,1:NLEV,:) = REGR_4D(:,:,:,:)
          ENDIF

          ! REGR_4D are the remapped fractions.
          DO T  = 1, NTIME
          DO L  = 1, NLEV
          DO J  = 1, HcoState%NY
          DO I2 = 1, HcoState%NX
             IF ( REGFRACS(I2,J,L,T) > MAXFRACS(I2,J,L,T) ) THEN
                MAXFRACS(I2,J,L,T) = REGR_4D(I2,J,L,T)
                INDECES (I2,J,L,T) = IVAL
             ENDIf
          ENDDO
          ENDDO
          ENDDO
          ENDDO

!------------------------------------------------------------------------------
! Prior to 9/29/16:
!          ! This code is preblematic in Gfortran.  Replace it with the
!          ! explicit DO loops above.  Leave this here for reference.
!          ! (sde, bmy, 9/21/16)
!          WHERE ( REGFRACS > MAXFRACS )
!             MAXFRACS = REGR_4D
!             INDECES  = IVAL
!          END WHERE
!------------------------------------------------------------------------------
       ENDIF

    ENDDO !I

    ! For index values, pass index data to data container.
    IF ( IsIndex ) THEN
       IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
          DO T = 1, NTIME
             NZ = SIZE(Lct%Dct%Dta%V3(T)%Val,3)
             Lct%Dct%Dta%V3(T)%Val(:,:,:) = INDECES(:,:,1:NZ,T)
          ENDDO
       ELSE
          DO T = 1, NTIME
             Lct%Dct%Dta%V2(T)%Val(:,:)   = INDECES(:,:,1,T)
          ENDDO
       ENDIF
    ENDIF

    ! Cleanup
    DEALLOCATE(LonEdgeI, LatEdgeI)
    IF ( ASSOCIATED( REGR_4D  ) ) DEALLOCATE( REGR_4D  )
    IF ( ALLOCATED ( FRACS    ) ) DEALLOCATE( FRACS    )
    IF ( ALLOCATED ( REGFRACS ) ) DEALLOCATE( REGFRACS )
    IF ( ALLOCATED ( MAXFRACS ) ) DEALLOCATE( MAXFRACS )
    IF ( ALLOCATED ( INDECES  ) ) DEALLOCATE( INDECES  )
    IF ( ALLOCATED ( UNIQVALS ) ) DEALLOCATE( UNIQVALS )

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE REGRID_MAPA2A
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ModelLev_Check
!
! !DESCRIPTION: Subroutine ModelLev\_Check checks if the passed number of
! vertical levels indicates that these are model levels or not.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ModelLev_Check( HcoState, nLev, IsModelLev, RC )
!
! !USES:
!
    USE HCO_FileData_Mod,   ONLY : FileData_ArrCheck
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState          ! HEMCO state object
    INTEGER,          INTENT(IN   )  :: nlev              ! number of levels
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(INOUT)  :: IsModelLev        ! Are these model levels?
    INTEGER,          INTENT(INOUT)  :: RC                ! Success or failure?
!
! !REVISION HISTORY:
!  29 Sep 2015 - C. Keller   - Initial version
!  22 May 2017 - R. Yantosca - Bug fix: Add MERRA2 to the #elif statement

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: nz

    !=================================================================
    ! ModelLev_Check begins here
    !=================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! If IsModelLev is already TRUE, nothing to do
    IF ( IsModelLev ) RETURN

    ! Shadow number of vertical levels on grid
    nz = HcoState%NZ

    ! Assume model levels if input data levels correspond to # of grid
    ! levels or levels + 1 (edges)
    IF ( nlev == nz .OR. nlev == nz + 1 ) THEN
       IsModelLev = .TRUE.
       RETURN
    ENDIF

    ! Other supported levels that depend on compiler flags
    ! Full grid
    IF ( nz == 72 ) THEN
       IF ( nlev <= 73 ) THEN
          IsModelLev = .TRUE.
       ENDIF

    ! Reduced grid
    ELSEIF ( nz == 47 ) THEN
       IF ( nlev == 72 .OR. &
            nlev == 73 .OR. &
            nlev <= 47       ) THEN
          IsModelLev = .TRUE.
       ENDIF
    ENDIF

  END SUBROUTINE ModelLev_Check
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ModelLev_Interpolate
!
! !DESCRIPTION: Subroutine ModelLev\_Interpolate puts 3D data from an
! arbitrary number of model levels onto the vertical levels of the simulation
! grid. Since the input data is already on model levels, this is only to
! inflate/collapse fields between native/reduced vertical levels, e.g. from
! 72 native GEOS-5 levels onto the reduced 47 levels. The vertical
! interpolation scheme depends on compiler switches. If none of the compiler
! switches listed below is used, no vertical interpolation is performed,
! e.g. the vertical levels of the input grid are retained.
!\\
!\\
! The input data (REGR\_4D) is expected to be already regridded horizontally.
! The 4th dimension of REGR\_4D denotes time.
!\\
!\\
! The 3rd dimension of REGR\_3D holds the vertical levels. It is assumed that
! these are model levels, starting at the surface (level 1). If the input
! data holds 72 input levels, this is interpreted as native data and will
! be collapsed onto the reduced grid. If the input data holds X <=47 levels,
! these levels are interpreted as levels 1-X of the reduced grid. In other
! words, input data with 33 levels will be interpreted as 33 levels on the
! reduced grid, and the data is accordingly mapped onto the simulation grid.
! If data becomes inflated or collapsed, the output data will always extent
! over all vertical levels of the simulation grid. If necessary, the unused
! upper levels will be filled with zeros. If no data interpolation is needed,
! the vertical extent of the output data is limited to the number of used
! levels. For instance, if the input data has 5 vertical levels, the output
! array will only extent over those 5 (bottom) levels.
!\\
!\\
! Currently, this routine can remap the following combinations:
!\begin{itemize}
! \item Native  GEOS-5 onto reduced GEOS-5 (72 --> 47 levels)
! \item Reduced GEOS-5 onto native  GEOS-5 (47 --> 72 levels)
! \item Native  GEOS-4 onto reduced GEOS-4 (55 --> 30 levels)
! \item Reduced GEOS-4 onto native  GEOS-4 (30 --> 55 levels)
! \item Native  GEOS-5 onto native  GEOS-4 (72 --> 55 levels)
! \item Reduced GEOS-5 onto native  GEOS-4 (47 --> 55 levels)
! \item Native  GEOS-5 onto reduced GEOS-4 (72 --> 30 levels)
! \item Reduced GEOS-5 onto reduced GEOS-4 (47 --> 30 levels)
!\end{itemize}
! Interpolation from GEOS-5 onto GEOS-4 levels is currently not supported.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ModelLev_Interpolate( HcoState, REGR_4D, Lct, RC )
!
! !USES:
!
    USE HCO_FileData_Mod,   ONLY : FileData_ArrCheck
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState          ! HEMCO state object
    REAL(sp),         POINTER        :: REGR_4D(:,:,:,:)  ! 4D input data
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER        :: Lct               ! HEMCO list container
    INTEGER,          INTENT(INOUT)  :: RC                ! Success or failure?
!
! !REVISION HISTORY:
!  30 Dec 2014 - C. Keller   - Initial version
!  24 Feb 2015 - R. Yantosca - Now exit if vertical interpolation isn't needed
!  12 Aug 2015 - R. Yantosca - Vertically remap MERRA2 as we do for GEOS-FP
!  24 Sep 2015 - C. Keller   - Added interpolation on edges.
!  06 Dec 2015 - C. Keller   - Pass # of GEOS-5 levels to be mapped onto GEOS-4
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: nx, ny, nz, nt
    INTEGER                 :: minlev, nlev, nout
    INTEGER                 :: L, T, NL
    INTEGER                 :: OS
    INTEGER                 :: G5T4
    LOGICAL                 :: verb, infl, clps
    LOGICAL                 :: DONE
    CHARACTER(LEN=255)      :: MSG

    !=================================================================
    ! ModelLev_Interpolate begins here
    !=================================================================

    ! Enter
    CALL HCO_ENTER (HcoState%Config%Err,&
                   'ModelLev_Interpolate (hco_interp_mod.F90)' , RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Check for verbose mode
    verb = HCO_IsVerb(HcoState%Config%Err,  3 )
    IF ( verb ) THEN
       MSG = 'Vertically interpolate model levels: '//TRIM(Lct%Dct%cName)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Get HEMCO grid dimensions
    nx = HcoState%NX
    ny = HcoState%NY
    nz = HcoState%NZ

    ! Variable G5T4 is the # of GEOS-5 levels that need to be mapped
    ! onto GEOS-4 levels.
    G5T4 = 0

    ! Input data must be on horizontal HEMCO grid
    IF ( SIZE(REGR_4D,1) /= nx ) THEN
       WRITE(MSG,*) 'x dimension mismatch ', TRIM(Lct%Dct%cName), &
          ': ', nx, SIZE(REGR_4D,1)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF
    IF ( SIZE(REGR_4D,2) /= ny ) THEN
       WRITE(MSG,*) 'y dimension mismatch ', TRIM(Lct%Dct%cName), &
          ': ', ny, SIZE(REGR_4D,2)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

    ! Get vertical and time dimension of input data
    nlev = SIZE(REGR_4D,3)
    nt   = SIZE(REGR_4D,4)

    ! Vertical interpolation done?
    DONE = .FALSE.

    !===================================================================
    ! If no vertical interpolation is needed, then (1) save the 4D
    ! input data array to to the HEMCO list container object and
    ! (2) exit this subroutine.
    !===================================================================
    IF ( ( nlev == nz ) .OR. ( nlev == nz+1 ) ) THEN

       CALL FileData_ArrCheck( HcoState%Config, Lct%Dct%Dta, nx, ny, nlev, nt, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       DO T = 1, nt
          Lct%Dct%Dta%V3(T)%Val(:,:,:) = REGR_4D(:,:,:,T)
       ENDDO

       ! Verbose
       IF ( HCO_IsVerb(HcoState%Config%Err, 3) ) THEN
          MSG = '# of input levels = # of output levels - passed as is.'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       ! Done!
       DONE = .TRUE.
    ENDIF

    !===================================================================
    ! Do vertical regridding:
    !===================================================================
    IF ( .NOT. DONE ) THEN

       !----------------------------------------------------------------
       ! Native levels
       !----------------------------------------------------------------
       IF ( nz == 72 ) THEN

          ! Determine number of output levels. If the input data has
          ! 47 or less levels, it is assumed to represent reduced
          ! GEOS-5 levels and data is mapped accordingly. If input data
          ! has more than 47 levels, it cannot be on the reduced grid
          ! and mapping is done 1:1
          IF ( nlev > 36 .AND. nlev <= 48 ) THEN
             IF ( nlev == 48 ) THEN
                nz   = nz + 1
                nout = nz
                NL   = 37
             ELSE
                nout = nz
                NL   = 36
             ENDIF
          ELSE
             nout = nlev
             NL   = nout
          ENDIF

          ! Make sure output array is allocated
          CALL FileData_ArrCheck( HcoState%Config, Lct%Dct%Dta, nx, ny, nout, nt, RC )

          ! Do for every time slice
          DO T = 1, nt

             ! Levels that are passed level-by-level.
             DO L = 1, NL
                Lct%Dct%Dta%V3(T)%Val(:,:,L) = REGR_4D(:,:,L,T)
             ENDDO !L

             ! If needed, inflate from reduced GEOS-5 grid onto native GEOS-5
             IF ( ( NL == 36 .AND. nz == 72 ) .OR. &
                  ( NL == 37 .AND. nz == 73 )       ) THEN
                ! Distribute over 2 levels (e.g. level 38 into 39-40):
                CALL INFLATE( Lct, REGR_4D, NL+1 , NL+1, 2, T )
                CALL INFLATE( Lct, REGR_4D, NL+2 , NL+3, 2, T )
                CALL INFLATE( Lct, REGR_4D, NL+3 , NL+5, 2, T )
                CALL INFLATE( Lct, REGR_4D, NL+4 , NL+7, 2, T )
                ! Distribute over 4 levels:
                CALL INFLATE( Lct, REGR_4D, NL+5 , NL+9, 4, T )
                CALL INFLATE( Lct, REGR_4D, NL+6 , NL+13, 4, T )
                CALL INFLATE( Lct, REGR_4D, NL+7 , NL+17, 4, T )
                CALL INFLATE( Lct, REGR_4D, NL+8 , NL+21, 4, T )
                CALL INFLATE( Lct, REGR_4D, NL+9 , NL+25, 4, T )
                CALL INFLATE( Lct, REGR_4D, NL+10, NL+29, 4, T )
                CALL INFLATE( Lct, REGR_4D, NL+11, NL+33, 4, T )
             ENDIF

          ENDDO ! T

          ! Verbose
          IF ( HCO_IsVerb(HcoState%Config%Err, 3) ) THEN
             WRITE(MSG,*) 'Mapped ', nlev, ' levels onto native GEOS-5 levels.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ! Done!
          DONE = .TRUE.

       !----------------------------------------------------------------
       ! Reduced levels
       !----------------------------------------------------------------
       ELSEIF ( nz == 47 ) THEN

          ! Determine number of output levels. If input data is on the
          ! native grid, we collapse them onto the reduced GEOS-5 grid.
          ! In all other cases, we assume the input data is already on
          ! the reduced levels and mappings occurs 1:1.
          IF ( nlev == 72 ) THEN
             nout = nz
             NL   = 36
          ELSEIF ( nlev == 73 ) THEN
             nz   = nz + 1
             nout = nz
             NL   = 37
          ELSEIF ( nlev > 47 ) THEN
             MSG = 'Can only remap from native onto reduced GEOS-5 if '// &
                   'input data has exactly 72 or 73 levels: '//TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ELSE
             nout = nlev
             NL   = nout
          ENDIF

          ! Make sure output array is allocated
          CALL FileData_ArrCheck( HcoState%Config, Lct%Dct%Dta, nx, ny, nout, nt, RC )

          ! Do for every time slice
          DO T = 1, nt

             ! Levels that are passed level-by-level.
             DO L = 1, NL
                Lct%Dct%Dta%V3(T)%Val(:,:,L) = REGR_4D(:,:,L,T)
             ENDDO !L

             ! If needed, collapse from native GEOS-5 onto reduced GEOS-5
             IF ( nlev == 72 .OR. nlev == 73 ) THEN

                ! Add one level offset if these are edges
                IF ( nlev == 73 ) THEN
                   OS = 1
                ELSE
                   OS = 0
                ENDIF

                ! Collapse two levels (e.g. levels 39-40 into level 38):
                CALL COLLAPSE( Lct, REGR_4D, 37+OS, 37+OS, 2, T, 5 )
                CALL COLLAPSE( Lct, REGR_4D, 38+OS, 39+OS, 2, T, 5 )
                CALL COLLAPSE( Lct, REGR_4D, 39+OS, 41+OS, 2, T, 5 )
                CALL COLLAPSE( Lct, REGR_4D, 40+OS, 43+OS, 2, T, 5 )
                ! Collapse four levels:
                CALL COLLAPSE( Lct, REGR_4D, 41+OS, 45+OS, 4, T, 5 )
                CALL COLLAPSE( Lct, REGR_4D, 42+OS, 49+OS, 4, T, 5 )
                CALL COLLAPSE( Lct, REGR_4D, 43+OS, 53+OS, 4, T, 5 )
                CALL COLLAPSE( Lct, REGR_4D, 44+OS, 57+OS, 4, T, 5 )
                CALL COLLAPSE( Lct, REGR_4D, 45+OS, 61+OS, 4, T, 5 )
                CALL COLLAPSE( Lct, REGR_4D, 46+OS, 65+OS, 4, T, 5 )
                CALL COLLAPSE( Lct, REGR_4D, 47+OS, 69+OS, 4, T, 5 )

             ENDIF
          ENDDO ! T

          ! Verbose
          IF ( HCO_IsVerb(HcoState%Config%Err, 3) ) THEN
             WRITE(MSG,*) 'Mapped ', nlev, ' levels onto reduced GEOS-5 levels.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ! Done!
          DONE = .TRUE.
       ENDIF

    ENDIF ! Vertical regridding required

    !===================================================================
    ! For all other cases, do not do any vertical regridding
    !===================================================================
    IF ( .NOT. DONE ) THEN
       CALL FileData_ArrCheck( HcoState%Config, Lct%Dct%Dta, nx, ny, nlev, nt, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       DO T = 1, nt
          Lct%Dct%Dta%V3(T)%Val(:,:,:) = REGR_4D(:,:,:,T)
       ENDDO

       ! Verbose
       IF ( HCO_IsVerb(HcoState%Config%Err, 3) ) THEN
          WRITE(MSG,*) 'Could not find vertical interpolation key - ', &
                       'kept the original ', nlev, ' levels.'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       ! Done!
       DONE = .TRUE.
    ENDIF

    !===================================================================
    ! Error check / verbose mode
    !===================================================================
    IF ( DONE ) THEN
      IF ( HCO_IsVerb(HcoState%Config%Err, 2) ) THEN
          WRITE(MSG,*) 'Did vertical regridding for ',TRIM(Lct%Dct%cName),':'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,*) 'Number of original levels: ', nlev
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,*) 'Number of output levels: ', SIZE(Lct%Dct%Dta%V3(1)%Val,3)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
    ELSE
       WRITE(MSG,*) 'Vertical regridding failed: ',TRIM(Lct%Dct%cName)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

    ! Return w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE ModelLev_Interpolate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS5_TO_GEOS4_LOWLEV
!
! !DESCRIPTION: Helper routine to map the lowest 28 GEOS-5 levels onto the
! lowest 11 GEOS-4 levels. The individual level weights were calculated
! offline and are hard-coded here.
! These are the edge pressure values on the lowest 28 GEOS-5 levels:
! 1013.25, 998.05, 982.76, 967.47, 952.19, 936.91
!  921.62, 906.34, 891.05, 875.77, 860.49, 845.21,
!  829.92, 809.55, 784.08, 758.62, 733.15, 707.69,
!  682.23, 644.05, 605.87, 567.70, 529.54, 491.40,
!  453.26, 415.15, 377.07, 339.00, 288.92
!
! And these are the edge pressure values on the lowest 12 GEOS-4 levels:
! 1013.25, 998.16, 968.49, 914.79, 841.15, 752.89,
!  655.96, 556.85, 472.64, 401.14, 340.43, 288.92
!
! The value at every given GEOS-4 level is determined from the GEOS-5 values
! by multiplying the (GEOS-5) input data by the normalized level weights. For
! instance, the first GEOS-5 level is the only level contributing to the 1st
! GEOS-4 level. For the 2nd GEOS-4 level, contributions from GEOS-5 levels
! 1-3 are used. Of GEOS-5 level 1, only 0.7% lies in level 2 of GEOS-4 (99.3%
! is in GEOS-4 level 1), whereas 100% of GEOS-5 level 2 and 93.3% of GEOS-5
! level 3 contribute to GEOS-4 level 2. The corresponding normalized weights
! become 0.00378,0.515, and 0.481, respectively.
!\\
!\\
! The weights don't always add up to exactly 1.00 due to rounding errors.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS5_TO_GEOS4_LOWLEV( HcoState, Lct, REGR_4D, NZ, T, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState          ! HEMCO state object
    REAL(sp),         POINTER        :: REGR_4D(:,:,:,:)  ! 4D input data
    INTEGER,          INTENT(IN)     :: T                 ! Time index
    INTEGER,          INTENT(IN)     :: NZ                ! # of vertical levels to remap. Must be 28 or 29
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER        :: Lct               ! HEMCO list container
    INTEGER,          INTENT(INOUT)  :: RC                ! Return code
!
! !REVISION HISTORY:
!  07 Jan 2015 - C. Keller   - Initial version.
!  24 Sep 2015 - C. Keller   - Added option to interpolate edges.
!  06 Dec 2015 - C. Keller   - Added input argument NZ
!EOP
!------------------------------------------------------------------------------
!BOC
    REAL(hp)           :: WGHT
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'GEOS5_TO_GEOS4_LOWLEV (hco_interp_mod.F90)'

    !=================================================================
    ! GEOS5_TO_GEOS4_LOWLEV begins here
    !=================================================================

    ! Check number of levels to be used
    IF ( NZ /= 28 .AND. NZ /= 29 ) THEN
       MSG = 'Cannot map GEOS-5 onto GEOS-4 data, number of levels must be 28 or 29: '//TRIM(Lct%Dct%cName)
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Error check: make sure array REGR_4D has at least NZ levels
    IF ( SIZE(REGR_4D,3) < NZ ) THEN
       WRITE(MSG,*) 'Cannot map GEOS-5 onto GEOS-4 data, original data has not enough levels: ', &
          TRIM(Lct%Dct%cName), ' --> ', SIZE(REGR_4D,3), ' smaller than ', NZ
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Map 28 GEOS-5 levels onto 11 GEOS-4 levels (grid midpoints):
    IF ( NZ == 28 ) THEN

       ! Reset
       Lct%Dct%Dta%V3(T)%Val(:,:,1:11) = 0.0_sp

       ! Level 1:
       Lct%Dct%Dta%V3(T)%Val(:,:, 1) = REGR_4D(:,:,1,T)

       ! Level 2:
       Lct%Dct%Dta%V3(T)%Val(:,:, 2) = 3.78e-3_sp * REGR_4D(:,:, 1,T) &
                                     + 0.515_sp   * REGR_4D(:,:, 2,T) &
                                     + 0.481_sp   * REGR_4D(:,:, 3,T)

       ! Level 3:
       Lct%Dct%Dta%V3(T)%Val(:,:, 3) = 1.88e-2_sp * REGR_4D(:,:, 3,T) &
                                     + 0.285_sp   * REGR_4D(:,:, 4,T) &
                                     + 0.285_sp   * REGR_4D(:,:, 5,T) &
                                     + 0.285_sp   * REGR_4D(:,:, 6,T) &
                                     + 0.127_sp   * REGR_4D(:,:, 7,T)

       ! Level 4:
       Lct%Dct%Dta%V3(T)%Val(:,:, 4) = 0.115_sp   * REGR_4D(:,:, 7,T) &
                                     + 0.208_sp   * REGR_4D(:,:, 8,T) &
                                     + 0.208_sp   * REGR_4D(:,:, 9,T) &
                                     + 0.208_sp   * REGR_4D(:,:,10,T) &
                                     + 0.208_sp   * REGR_4D(:,:,11,T) &
                                     + 5.51e-2_sp * REGR_4D(:,:,12,T)

       ! Level 5:
       Lct%Dct%Dta%V3(T)%Val(:,:, 5) = 0.189_sp   * REGR_4D(:,:,12,T) &
                                     + 0.253_sp   * REGR_4D(:,:,13,T) &
                                     + 0.253_sp   * REGR_4D(:,:,14,T) &
                                     + 0.253_sp   * REGR_4D(:,:,15,T) &
                                     + 5.68e-2_sp * REGR_4D(:,:,16,T)

       ! Level 6:
       Lct%Dct%Dta%V3(T)%Val(:,:, 6) = 0.224_sp   * REGR_4D(:,:,16,T) &
                                     + 0.289_sp   * REGR_4D(:,:,17,T) &
                                     + 0.289_sp   * REGR_4D(:,:,18,T) &
                                     + 0.199_sp   * REGR_4D(:,:,19,T)

       ! Level 7:
       Lct%Dct%Dta%V3(T)%Val(:,:, 7) = 0.120_sp   * REGR_4D(:,:,19,T) &
                                     + 0.385_sp   * REGR_4D(:,:,20,T) &
                                     + 0.385_sp   * REGR_4D(:,:,21,T) &
                                     + 0.110_sp   * REGR_4D(:,:,22,T)

       ! Level 8:
       Lct%Dct%Dta%V3(T)%Val(:,:, 8) = 0.324_sp   * REGR_4D(:,:,22,T) &
                                     + 0.453_sp   * REGR_4D(:,:,23,T) &
                                     + 0.223_sp   * REGR_4D(:,:,24,T)

       ! Level 9:
       Lct%Dct%Dta%V3(T)%Val(:,:, 9) = 0.271_sp   * REGR_4D(:,:,24,T) &
                                     + 0.533_sp   * REGR_4D(:,:,25,T) &
                                     + 0.196_sp   * REGR_4D(:,:,26,T)

       ! Level 10:
       Lct%Dct%Dta%V3(T)%Val(:,:,10) = 0.396_sp   * REGR_4D(:,:,26,T) &
                                     + 0.604_sp   * REGR_4D(:,:,27,T)

       ! Level 11:
       Lct%Dct%Dta%V3(T)%Val(:,:,11) = 3.63e-2_sp * REGR_4D(:,:,27,T) &
                                     + 0.964_sp   * REGR_4D(:,:,28,T)

    ! Map 29 GEOS-5 levels onto 12 GEOS-4 levels (grid edges):
    ELSEIF ( NZ == 29 ) THEN

       ! Reset
       Lct%Dct%Dta%V3(T)%Val(:,:,1:12) = 0.0_sp

       ! Level 1
       Lct%Dct%Dta%V3(T)%Val(:,:, 1) = REGR_4D(:,:,1,T)

       ! Level 2:
       Lct%Dct%Dta%V3(T)%Val(:,:, 2) = 5.01e-3_sp * REGR_4D(:,:, 1,T) &
                                     + 0.680_sp   * REGR_4D(:,:, 2,T) &
                                     + 0.314_sp   * REGR_4D(:,:, 3,T)

       ! Level 3:
       Lct%Dct%Dta%V3(T)%Val(:,:, 3) = 0.197_sp   * REGR_4D(:,:, 3,T) &
                                     + 0.366_sp   * REGR_4D(:,:, 4,T) &
                                     + 0.366_sp   * REGR_4D(:,:, 5,T) &
                                     + 6.98e-2_sp * REGR_4D(:,:, 6,T)

       ! Level 4:
       Lct%Dct%Dta%V3(T)%Val(:,:, 4) = 0.194_sp   * REGR_4D(:,:, 6,T) &
                                     + 0.240_sp   * REGR_4D(:,:, 7,T) &
                                     + 0.240_sp   * REGR_4D(:,:, 8,T) &
                                     + 0.240_sp   * REGR_4D(:,:, 9,T) &
                                     + 8.55e-2_sp * REGR_4D(:,:,10,T)

       ! Level 5:
       Lct%Dct%Dta%V3(T)%Val(:,:, 5) = 0.139_sp   * REGR_4D(:,:,10,T) &
                                     + 0.216_sp   * REGR_4D(:,:,11,T) &
                                     + 0.216_sp   * REGR_4D(:,:,12,T) &
                                     + 0.216_sp   * REGR_4D(:,:,13,T) &
                                     + 0.214_sp   * REGR_4D(:,:,14,T)

       ! Level 6:
       Lct%Dct%Dta%V3(T)%Val(:,:, 6) = 2.20e-2_sp * REGR_4D(:,:,14,T) &
                                     + 0.275_sp   * REGR_4D(:,:,15,T) &
                                     + 0.275_sp   * REGR_4D(:,:,16,T) &
                                     + 0.275_sp   * REGR_4D(:,:,17,T) &
                                     + 0.173_sp   * REGR_4D(:,:,18,T)

       ! Level 7:
       Lct%Dct%Dta%V3(T)%Val(:,:, 7) = 0.130_sp   * REGR_4D(:,:,18,T) &
                                     + 0.345_sp   * REGR_4D(:,:,19,T) &
                                     + 0.345_sp   * REGR_4D(:,:,20,T) &
                                     + 0.170_sp   * REGR_4D(:,:,21,T)

       ! Level 8:
       Lct%Dct%Dta%V3(T)%Val(:,:, 8) = 0.214_sp   * REGR_4D(:,:,21,T) &
                                     + 0.416_sp   * REGR_4D(:,:,22,T) &
                                     + 0.370_sp   * REGR_4D(:,:,23,T)

       ! Level 9:
       Lct%Dct%Dta%V3(T)%Val(:,:, 9) = 5.49e-2_sp * REGR_4D(:,:,23,T) &
                                     + 0.490_sp   * REGR_4D(:,:,24,T) &
                                     + 0.455_sp   * REGR_4D(:,:,25,T)

       ! Level 10:
       Lct%Dct%Dta%V3(T)%Val(:,:,10) = 4.06e-2_sp * REGR_4D(:,:,25,T) &
                                     + 0.576_sp   * REGR_4D(:,:,26,T) &
                                     + 0.383_sp   * REGR_4D(:,:,27,T)

       ! Level 11:
       Lct%Dct%Dta%V3(T)%Val(:,:,11) = 0.254_sp   * REGR_4D(:,:,27,T) &
                                     + 0.746_sp   * REGR_4D(:,:,28,T)

       ! Level 12:
       Lct%Dct%Dta%V3(T)%Val(:,:,12) = 1.60e-2_sp * REGR_4D(:,:,28,T) &
                                     + 0.984_sp   * REGR_4D(:,:,29,T)

    ENDIF

    ! Return with success
    RC = HCO_SUCCESS

  END SUBROUTINE GEOS5_TO_GEOS4_LOWLEV
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: COLLAPSE
!
! !DESCRIPTION: Helper routine to collapse input levels onto the output grid.
! The input data is weighted by the grid box thicknesses defined on top of
! this module. The input parameter T determines the time slice to be considered,
! and MET denotes the met field type of the input data (4 = GEOS-4 levels, GEOS-5
! otherwise).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE COLLAPSE ( Lct, REGR_4D, OutLev, InLev1, NLEV, T, MET )
!
! !INPUT PARAMETERS:
!
    REAL(sp),         POINTER        :: REGR_4D(:,:,:,:)  ! 4D input data
    INTEGER,          INTENT(IN)     :: OutLev
    INTEGER,          INTENT(IN)     :: InLev1
    INTEGER,          INTENT(IN)     :: NLEV
    INTEGER,          INTENT(IN)     :: T
    INTEGER,          INTENT(IN)     :: MET               ! 4=GEOS-4, else GEOS-5
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER        :: Lct               ! HEMCO list container
!
! !REVISION HISTORY:
!  30 Dec 2014 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER               :: I, NZ, ILEV, TOPLEV
    REAL(hp)              :: THICK
    REAL(hp), POINTER     :: EDG(:)
    REAL(hp), ALLOCATABLE :: WGT(:)

    !=================================================================
    ! COLLAPSE begins here
    !=================================================================

    ! Init
    EDG => NULL()

    ! Reset
    Lct%Dct%Dta%V3(T)%Val(:,:,OutLev) = 0.0_hp

    ! Don't do anything if there are not enough levels in REGR_4D
    NZ = SIZE(REGR_4D,3)
    IF ( NZ < InLev1 ) RETURN

    ! Get maximum level to be used for pressure thickness calculations.
    TOPLEV = InLev1 + ( NLEV-1 )

    ! Get pointer to grid edges on the native input grid
    IF ( Met == 4 ) THEN
       EDG => G4_EDGE_NATIVE(InLev1:TOPLEV)
    ELSE
       EDG => G5_EDGE_NATIVE(InLev1:TOPLEV)
    ENDIF

    ! Thickness of output level
    THICK = EDG(1) - EDG(NLEV)

    ! Get level weights
    ALLOCATE(WGT(NLEV))
    WGT = 0.0
    DO I = 1, NLEV-1
       WGT(I) = ( EDG(I) - EDG(I+1) ) / THICK
    ENDDO

    ! Pass levels to output data, one after each other
    Lct%Dct%Dta%V3(T)%Val(:,:,OutLev) = REGR_4D(:,:,InLev1,T) * WGT(1)
    DO I = 1, NLEV-1
       ILEV = InLev1 + I
       IF ( NZ < ILEV ) EXIT
       Lct%Dct%Dta%V3(T)%Val(:,:,OutLev) = Lct%Dct%Dta%V3(T)%Val(:,:,OutLev) &
                                         + ( REGR_4D(:,:,ILEV,T) * WGT(I+1) )
    ENDDO

    ! Cleanup
    DEALLOCATE(WGT)
    EDG => NULL()

  END SUBROUTINE COLLAPSE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INFLATE
!
! !DESCRIPTION: Helper routine to inflate input levels onto the output grid.
! The values on the input data are evenly distributed amongst all output
! levels.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INFLATE ( Lct, REGR_4D, InLev, OutLev1, NLEV, T )
!
! !INPUT PARAMETERS:
!
    REAL(sp),         POINTER        :: REGR_4D(:,:,:,:)  ! 4D input data
    INTEGER,          INTENT(IN)     :: InLev
    INTEGER,          INTENT(IN)     :: OutLev1
    INTEGER,          INTENT(IN)     :: NLEV
    INTEGER,          INTENT(IN)     :: T
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER        :: Lct               ! HEMCO list container
!
! !REVISION HISTORY:
!  30 Dec 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: I, NZ, ILEV

    !=================================================================
    ! INFLATE begins here
    !=================================================================

    ! Get input data array
    NZ = SIZE(REGR_4D,3)

    ! Do for every output level
    DO I = 1, NLEV

       ! Current output level
       ILEV = OutLev1 + I - 1

       ! If input level is beyond vert. extent of input data, set output
       ! data to zero.
       IF ( InLev > NZ ) THEN
          Lct%Dct%Dta%V3(T)%Val(:,:,ILEV) = 0.0_hp

       ! Otherwise, evenly distribute input data
       ELSE
          Lct%Dct%Dta%V3(T)%Val(:,:,ILEV) = REGR_4D(:,:,InLev,T)
       ENDIF
    ENDDO

  END SUBROUTINE INFLATE
!EOC
END MODULE HCO_Interp_Mod
