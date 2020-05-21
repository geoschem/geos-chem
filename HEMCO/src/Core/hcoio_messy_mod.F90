!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_messy_mod.F90
!
! !DESCRIPTION: Module HCOIO\_MESSY\_MOD interfaces HEMCO with the regridding
! tool NCREGRID of the Modular Earth Submodel System (MESSy). This regridding
! scheme is used for vertical regridding and/or for index data, i.e. data with
! discrete values (e.g. land type integers). This code currently only works for
! rectilinear (regular lon-lat) grids but can be extended to support curvilinear
! grids.
!\\
! REFERENCES:
! \begin{itemize}
! \item Joeckel, P. Technical note: Recursive rediscretisation of geo-
! scientific data in the Modular Earth Submodel System (MESSy), ACP, 6,
! 3557--3562, 2006.
! \end{itemize}
! !INTERFACE:
!
MODULE HCOIO_MESSY_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_TYPES_MOD,        ONLY : ListCont
  USE HCO_STATE_MOD,        ONLY : Hco_State
  USE MESSY_NCREGRID_BASE,  ONLY : NARRAY, AXIS

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCO_MESSY_REGRID
!
! !PRIVATE MEMBER FUNCTIONS:
!
  ! Set this value to TRUE if you want to reduce the output array
  ! to the minimum required number of vertical levels.
  LOGICAL, PARAMETER            :: ReduceVert = .FALSE.
!
! !MODULE INTERFACES:
!
! !REVISION HISTORY:
!  24 Jun 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:

  CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Hco_Messy_Regrid
!
! !DESCRIPTION: This is the wrapper routine to regrid a 4D input array
! NcArr (x,y,z,t) onto the HEMCO emissions grid (defined in HcoState)
! using the regridding tool NCREGRID. LonEdge, LatEdge and LevEdge are
! the grid point edges of the input grid. The data is written into list
! container Lct.
!\\
!\\
! If the input grid is 2D (horizontal only), LevEdge must not be specified
! (null pointer) and the data is regridded in the horizontal only. If the
! input grid has only one vertical level, it is assumed that this is the
! surface level and the output data is 3D but with only one vertical level.
!\\
!\\
! For input data with more than one vertical level, the data is mapped
! onto the entire 3D grid. The module parameter ReduceVert can be used to
! cap the output data at the lowest possible level. For example, if the
! input grid only covers three surface levels with a minimum sigma value
! of 0.75, vertical regridding is performed within this sigma range (1-0.75)
! and the output grid is reduced accordingly. This option is not used in the
! standard HEMCO setup because problems can arise if the data array of a
! given container suddently changes its size (i.e. when updated data covers
! more/less vertical levels than the data beforehand).
!\\
!\\
! The input argument IsModelLev denotes whether or not the vertical
! coordinates of the input data are on model levels. If set to yes and LevEdge
! is not provided (i.e. a nullified pointer), the MESSy regridding routines are
! only used for the horizontal remapping and subroutine ModelLev\_Interpolate
! (module hco\_interp\_mod.F90) is used for the vertical remapping.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_MESSY_REGRID ( HcoState,  NcArr,               &
                                LonEdge,   LatEdge,    LevEdge, &
                                Lct,       IsModelLev, RC        )
!
! !USES:
!
  USE HCO_FILEDATA_MOD,     ONLY : FileData_ArrCheck
  USE HCO_UNIT_MOD,         ONLY : HCO_IsIndexData
  USE HCO_INTERP_MOD,       ONLY : ModelLev_Interpolate
  USE MESSY_NCREGRID_BASE,  ONLY : RG_INT, RG_IDX
  USE MESSY_NCREGRID_BASE,  ONLY : NREGRID
  USE MESSY_NCREGRID_BASE,  ONLY : INIT_NARRAY
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState        ! HEMCO obj.
    REAL(sp),         POINTER        :: ncArr(:,:,:,:)  ! Input array(x,y,z,t)
    REAL(hp),         POINTER        :: LonEdge(:)      ! lon edges
    REAL(hp),         POINTER        :: LatEdge(:)      ! lat edges
    REAL(hp),         POINTER        :: LevEdge(:,:,:)  ! sigma level edges
    TYPE(ListCont),   POINTER        :: Lct             ! Target list container
    LOGICAL,          INTENT(IN   )  :: IsModelLev      ! Are these model levels?
    INTEGER,          INTENT(INOUT)  :: RC              ! Return code
!
! !REVISION HISTORY:
!  27 Jun 2014 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ROUTINE ARGUMENTS:
!
    TYPE(narray), POINTER         :: narr_src(:)
    TYPE(narray), POINTER         :: narr_dst(:)
    TYPE(axis),   POINTER         :: axis_src(:)
    TYPE(axis),   POINTER         :: axis_dst(:)
    INTEGER,      POINTER         :: rg_type(:)
    INTEGER,      POINTER         :: rcnt   (:)
    REAL(dp),     POINTER         :: sovl   (:)
    REAL(dp),     POINTER         :: dovl   (:)
    REAL(hp),     POINTER         :: lon    (:)
    REAL(hp),     POINTER         :: lat    (:)
    REAL(hp),     POINTER         :: sigma  (:,:,:)
    REAL(sp),     POINTER         :: ArrIn  (:,:,:,:)
    REAL(sp),     POINTER         :: ArrOut (:,:,:,:)
    REAL(hp), ALLOCATABLE, TARGET :: sigout (:,:,:)
    REAL(hp)                      :: sigMin
    INTEGER                       :: NZIN, NZOUT, NTIME
    INTEGER                       :: NXIN, NYIN
    INTEGER                       :: I, L, AS
    INTEGER                       :: NCALLS
    CHARACTER(LEN=255)            :: MSG, LOC
    LOGICAL                       :: SameGrid, verb

    !=================================================================
    ! HCO_MESSY_REGRID begins here
    !=================================================================

    ! For error handling
    LOC = 'HCO_MESSY_REGRID (HCOI_MESSY_MOD.F90)'
    CALL HCO_ENTER ( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Init
    narr_src => NULL()
    narr_dst => NULL()
    axis_src => NULL()
    axis_dst => NULL()
    rg_type  => NULL()
    rcnt     => NULL()
    sovl     => NULL()
    dovl     => NULL()
    lon      => NULL()
    lat      => NULL()
    sigma    => NULL()
    ArrIn    => NULL()
    ArrOut   => NULL()

    ! verbose?
    verb = HCO_IsVerb(HcoState%Config%Err,3)

    ! Horizontal dimension of input data
    NXIN = SIZE(NcArr,1)
    NYIN = SIZE(NcArr,2)

    ! Number of vertical levels of input data
    NZIN = SIZE(NcArr,3)

    ! Number of time slices. All time slices will be regridded
    ! simultaneously.
    NTIME = SIZE(NcArr,4)

    ! Error check: data must be 2D or 3D.
    IF ( Lct%Dct%Dta%SpaceDim /= 2 .AND. &
         Lct%Dct%Dta%SpaceDim /= 3         ) THEN
       MSG = 'Can only regrid 2D or 3D data: ' // TRIM(Lct%Dct%cName)
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Shortcut if input field is already on output grid: directly
    ! pass data to Lct.
    ! NOTE: if the number of input levels matches the number of output
    ! levels (and the horizontal dimensions agree as well), it is
    ! assumed that they are on the same grid!
    !-----------------------------------------------------------------

    ! Input grid = output grid?
    SameGrid = .FALSE.

    ! Horizontal dimensions have to match
    IF ( (NXIN == HcoState%NX) .AND. (NYIN == HcoState%NY) ) THEN

       ! Vertical dimensions have to match or be 1
       IF ( NZIN == 1 .OR. NZIN == HcoState%NZ ) THEN

          ! Assume same grid
          SameGrid = .TRUE.

          ! Check for same boundaries. Otherwise falsify SameGrid
          ! Lon ...
          IF ( MINVAL(LonEdge) /= MINVAL(HcoState%Grid%XEDGE%Val) .OR. &
               MAXVAL(LonEdge) /= MAXVAL(HcoState%Grid%XEDGE%Val) ) THEN
             SameGrid = .FALSE.
          ENDIF
          ! ... Lat ...
          IF ( MINVAL(LatEdge) /= MINVAL(HcoState%Grid%YEDGE%Val) .OR. &
               MAXVAL(LatEdge) /= MAXVAL(HcoState%Grid%YEDGE%Val) ) THEN
             SameGrid = .FALSE.
          ENDIF
          ! ... Lev
          IF ( NZIN > 1 ) THEN
             ! TODO: Eventually need to add level check here.
             ! For now, assume that input levels are equal to
             ! output level if they have the same dimension!
          ENDIF
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Define number of vertical levels on output grid
    !-----------------------------------------------------------------

    ! Error check. If we are not on the same grid, LevEdge must be
    ! provided or IsModelLev must be set to true for 3D data.
    IF ( NZIN > 1 .AND. .NOT. SameGrid ) THEN
       IF ( .NOT. ASSOCIATED(LevEdge) .AND. .NOT. IsModelLev ) THEN
          MSG = 'Cannot regrid '//TRIM(Lct%Dct%cName)//'. Either level '//&
                'edges must be provided or data must be on model levels.'
          CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF
    ENDIF

    NZOUT = NZIN
    IF ( .NOT. SameGrid .AND. ASSOCIATED(LevEdge) ) THEN

       ! Calculate sigma level for each grid point on the output grid.
       ! This is the pressure at location i,j,l normalized by surface
       ! pressure @ i,j: sigma(i,j,l) = p(i,j,l) / ps(i,j)
       ALLOCATE(sigout(HcoState%NX,HcoState%NY,HcoState%NZ+1),STAT=AS)
       IF ( AS/= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate sigout', RC )
          RETURN
       ENDIF
       DO l = 1, HcoState%NZ+1
          sigout(:,:,l) = HcoState%Grid%PEDGE%Val(:,:,l) &
                        / HcoState%Grid%PEDGE%Val(:,:,1)
       ENDDO

       ! Now find first level on output grid where all sigma levels are
       ! lower than the lowest sigma value on the input grid. This is
       ! the highest level that needs to be considered for regridding.
       IF ( ReduceVert ) THEN
          sigMin = MINVAL(LevEdge)
          DO l = 1, HcoState%NZ+1
             IF ( MINVAL(sigout(:,:,l)) < sigMin ) EXIT
          ENDDO

          ! The output grid is at grid center, so use l-1. Must be at least one.
          NZOUT = max(1,l-1)

       ! Use full vertical grid if vertical levels shall not be restricted
       ! to range of input data (default).
       ELSE
          NZOUT = HcoState%NZ
       ENDIF
    ENDIF

    ! verbose mode
    IF ( verb ) THEN
       MSG = 'Do MESSy regridding: ' // TRIM(Lct%Dct%cName)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - SameGrid     ? ', SameGrid
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Model levels ? ', IsModelLev
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    !-----------------------------------------------------------------
    ! Make sure output array is defined & allocated
    !-----------------------------------------------------------------
    IF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN
       CALL FileData_ArrCheck( HcoState%Config, &
                               Lct%Dct%Dta, HcoState%NX, HcoState%NY, &
                               NTIME,       RC                         )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ELSEIF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
       CALL FileData_ArrCheck( HcoState%Config, &
                               Lct%Dct%Dta, HcoState%NX, HcoState%NY, &
                               NZOUT,       NTIME,       RC            )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Do straight-forward mapping if input grid = output grid
    !-----------------------------------------------------------------
    IF ( SameGrid ) THEN
       MSG = 'Input grid seems to match output grid. ' // &
             'No regridding is performed: ' // TRIM(Lct%Dct%cName)
       CALL HCO_WARNING( HcoState%Config%Err, MSG, RC )

       ! For every time slice...
       DO I = 1, NTIME
       DO L = 1, NZOUT

          ! 3D data
          IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
             Lct%Dct%Dta%V3(I)%Val(:,:,L) = NcArr(:,:,L,I)

          ! 2D data
          ELSE
             Lct%Dct%Dta%V2(I)%Val(:,:) = NcArr(:,:,L,I)
          ENDIF
       ENDDO
       ENDDO

       ! All done!
       CALL HCO_LEAVE ( HcoState%Config%Err, RC )
       RETURN
    ENDIF

    !=================================================================
    ! MESSy regridding follows below
    !=================================================================

    !-----------------------------------------------------------------
    ! Source grid description.
    ! This creates a MESSy axis object for the source grid.
    !-----------------------------------------------------------------
    lon   => LonEdge
    lat   => LatEdge
    sigma => LevEdge

    CALL AXIS_CREATE( HcoState, lon, lat, sigma, axis_src, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Free pointer
    lon   => NULL()
    lat   => NULL()
    sigma => NULL()

    !-----------------------------------------------------------------
    ! Destination grid description.
    ! This creates a MESSy axis object for the target (=HEMCO) grid.
    !-----------------------------------------------------------------

    ! Get horizontal grid directly from HEMCO state
    lon   => HcoState%Grid%XEDGE%Val(:,1)
    lat   => HcoState%Grid%YEDGE%Val(1,:)
    IF( ASSOCIATED(LevEdge) ) sigma => sigout(:,:,1:NZOUT+1)

    CALL AXIS_CREATE( HcoState, lon, lat, sigma, axis_dst, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Free pointer
    lon   => NULL()
    lat   => NULL()
    sigma => NULL()
    IF ( ALLOCATED(sigout) ) DEALLOCATE(sigout)

    !-----------------------------------------------------------------
    ! Set all other regridding parameter
    !-----------------------------------------------------------------

    ! rg_type denotes the regridding type for each array (i.e. time
    ! slice). Set to 'intensive quantity' for all concentrations (incl.
    ! unitless) data. Set to 'index distribution' for data marked as
    ! index data in the configuration file. This will remap discrete
    ! values without interpolation, i.e. each grid box on the new
    ! grid holds the value with most overlap in the original grid.
    ALLOCATE(rg_type(NTIME))
    IF ( HCO_IsIndexData(Lct%Dct%Dta%OrigUnit) ) THEN
       rg_type(:) = RG_IDX
       IF ( verb ) THEN
          MSG = ' - Remap as index data.'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
    ELSE
       rg_type(:) = RG_INT
       IF ( verb ) THEN
          MSG = ' - Remap as concentration data.'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Number of times the regridding need to be performed. If input
    ! data is on model levels, the data is only horizontally regridded,
    ! e.g. the regridding routine is called for every horizontal level
    ! separately. The vertical interpolation is done afterwards using
    ! routine ModelLev_Interpolate.
    !-----------------------------------------------------------------
    NCALLS = 1
    IF ( IsModelLev .AND. .NOT. ASSOCIATED(LevEdge) .AND. NZIN > 1 ) THEN
       NCALLS = NZIN
       ALLOCATE(ArrOut(HcoState%NX,HcoState%NY,NZIN,NTIME),STAT=AS)
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate ArrOut', RC )
          RETURN
       ENDIF
       ArrOut = 0.0_sp
    ENDIF

    ! Do for all level batches ...
    DO I = 1, NCALLS

       ! ArrIn is the input array to be used
       IF ( NCALLS /= 1 ) THEN
          ArrIn => NcArr(:,:,I:I,:)
       ELSE
          ArrIn => NcArr
       ENDIF

       !-----------------------------------------------------------------
       ! Map input array onto MESSy array. Different time slices are
       ! stored as individual vector elements of narr_src.
       !-----------------------------------------------------------------
       CALL HCO2MESSY( HcoState, ArrIn, narr_src, axis_src, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       !-----------------------------------------------------------------
       ! Do the regridding
       !-----------------------------------------------------------------
       CALL NREGRID(s=narr_src,      sax=axis_src, dax=axis_dst, d=narr_dst, &
                    rg_type=rg_type, sovl=sovl,    dovl=dovl,    rcnt=rcnt    )

       !-----------------------------------------------------------------
       ! Map the destination array narr_dst onto the data vector in the
       ! HEMCO list container or onto the temporary array ArrOut.
       !-----------------------------------------------------------------
       CALL MESSY2HCO( HcoState, narr_dst, Lct, I, ArrOut, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Cleanup
       ArrIn => NULL()

    ENDDO !NCALLS

    !-----------------------------------------------------------------
    ! If these are model levels, do vertical interpolation now
    !-----------------------------------------------------------------
    IF ( ASSOCIATED(ArrOut) ) THEN
       CALL ModelLev_Interpolate( HcoState, ArrOut, Lct, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       DEALLOCATE(ArrOut)
    ENDIF

    !-----------------------------------------------------------------
    ! Cleanup
    !-----------------------------------------------------------------
    DEALLOCATE( sovl, dovl, rcnt, STAT=AS)
    IF(AS/=0) THEN
       CALL HCO_ERROR(HcoState%Config%Err,'DEALLOCATION ERROR 1', RC )
       RETURN
    ENDIF
    NULLIFY(sovl, dovl, rcnt)

    DO I=1, SIZE(narr_dst)
       CALL INIT_NARRAY(narr_dst(I))
    ENDDO
    DEALLOCATE(narr_dst, STAT=AS)
    IF(AS/=0) THEN
       CALL HCO_ERROR(HcoState%Config%Err,'DEALLOCATION ERROR 3', RC )
       RETURN
    ENDIF
    NULLIFY(narr_dst)

    DO I=1, SIZE(narr_src)
       CALL INIT_NARRAY(narr_src(I))
    ENDDO
    DEALLOCATE(narr_src, STAT=AS)
    IF(AS/=0) THEN
       CALL HCO_ERROR(HcoState%Config%Err,'DEALLOCATION ERROR 2', RC )
       RETURN
    ENDIF
    NULLIFY(narr_src)

    DEALLOCATE( rg_type )
    rg_type => NULL()

    CALL AXIS_DELETE( axis_src, axis_dst, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE HCO_MESSY_REGRID
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Axis_Create
!
! !DESCRIPTION: Subroutine AXIS\_CREATE creates a MESSy axis type
! from the grid defined by mid points Lon, Lat, Lev. Lev must be in
! (unitless) sigma coordinates: sigma(i,j,l) = p(i,j,l) / p\_surface(i,j)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AXIS_CREATE( HcoState, lon, lat, lev, ax, RC )
!
! !USES:
!
  USE MESSY_NCREGRID_BASE,  ONLY : INIT_AXIS
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState
    TYPE(ListCont),   POINTER                 :: Lct
    REAL(hp),         POINTER                 :: Lon(:)
    REAL(hp),         POINTER                 :: Lat(:)
    REAL(hp),         POINTER                 :: Lev(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(axis),       POINTER                 :: ax(:)
    INTEGER,          INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  22 Jun 2014 - C. Keller - Initial version (from messy_ncregrid_geohyb.f90)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ROUTINE ARGUMENTS:
!
    INTEGER                :: fID, I, J, N, status
    INTEGER                :: LOW, UPP
    INTEGER                :: XLON, YLAT, XLEV, YLEV, ZLEV
    INTEGER                :: ndp, nlev, cnt
    INTEGER                :: vtype
    INTEGER                :: ndep_lon, ndep_lat
    CHARACTER(LEN=255)     :: MSG, LOC

    !=================================================================
    ! AXIS_CREATE begins here
    !=================================================================

    ! For error handling
    LOC = 'AXIS_CREATE (HCOI_MESSY_MOD.F90)'
    CALL HCO_ENTER ( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! Pass horizontal grid dimensions to local variables. For now,
    ! assume all grids to be regular, i.e. no curvilinear grids.
    ! In this case, lon and lat dimensions are 1D-vectors.
    ! ----------------------------------------------------------------

    ! ndep_lon and ndep_lat are the axis number of lon and lat. Needed
    ! if we have to set dependencies for the vertical axis.
    ndep_lon = 0
    ndep_lat = 0

    ! Count axes and get grid dimensions
    ! Last axis is always 'free' dimension. This is required for the
    ! regridding to work properly.
    N = 1
    IF ( ASSOCIATED(lon) ) N = N + 1
    IF ( ASSOCIATED(lat) ) N = N + 1
    IF ( ASSOCIATED(lev) ) N = N + 1

    ! ALLOCATE AXIS
    ALLOCATE(ax(N), STAT=status)
    IF ( status /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate axis', RC )
       RETURN
    ENDIF
    DO I=1, N
       CALL INIT_AXIS(ax(I))
    ENDDO

    ! N is the current axis
    N = 0

    ! ----------------------------------------------------------------
    ! Assign longitude: this is always the first dimension
    ! ----------------------------------------------------------------
    IF ( ASSOCIATED(lon) ) THEN
       N        = N + 1
       ax(N)%lm = .true.     ! LONGITUDE IS MODULO AXIS

       ! Axis dimension
       XLON = SIZE(lon,1)

       ! FOR NOW, ASSUME NO DEPENDENCIES. NEED TO EDIT HERE
       ! IF WE WANT TO USE CURVILINEAR GRIDS
       ax(N)%ndp    = 1          ! LONGITUDE IS ...
       ALLOCATE(ax(N)%dep(1), STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate lon dependencies', RC )
          RETURN
       ENDIF
       ax(N)%dep(1) = N          ! ... INDEPENDENT
       ndep_lon     = N

       ax(N)%dat%n = 1          ! 1 dimension
       ALLOCATE(ax(N)%dat%dim(ax(N)%dat%n), STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate lon dimensions', RC )
          RETURN
       ENDIF
       ax(N)%dat%dim(:) = XLON

       ALLOCATE(ax(N)%dat%vd(XLON),STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate lon axis', RC )
          RETURN
       ENDIF
       ax(N)%dat%vd(:) = lon

    ENDIF !lon

    ! ----------------------------------------------------------------
    ! Assign latitude: this is always the second dimension
    ! ----------------------------------------------------------------
    IF ( ASSOCIATED(lat) ) THEN
       N        = N + 1
       ax(N)%lm = .false.    ! LATITUDE IS NON-MODULO AXIS

       ! Axis dimension
       YLAT = SIZE(lat,1)

       ! FOR NOW, ASSUME NO DEPENDENCIES. NEED TO EDIT HERE
       ! IF WE WANT TO USE CURVILINEAR GRIDS
       ax(N)%ndp    = 1          ! LATITUDE IS ...
       ALLOCATE(ax(N)%dep(1), STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate lat dependencies', RC )
          RETURN
       ENDIF
       ax(N)%dep(1) = N          ! ... INDEPENDENT
       ndep_lat     = N

       ax(N)%dat%n = 1          ! 1 dimension
       ALLOCATE(ax(N)%dat%dim(ax(N)%dat%n), STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate lat dimensions', RC )
          RETURN
       ENDIF
       ax(N)%dat%dim(:) = YLAT

       ALLOCATE(ax(N)%dat%vd(YLAT),STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate lat axis', RC )
          RETURN
       ENDIF
       ax(N)%dat%vd(:) = lat

       ! TAKE INTO ACCOUNT SPHERICAL GEOMETRY ...
       ax(N)%dat%vd = COS( ( (ax(N)%dat%vd - 90.0_dp) / 180.0_dp ) * &
                           HcoState%Phys%PI )

    ENDIF !lat

    ! ----------------------------------------------------------------
    ! Assign vertical levels (if defined): this is the 3rd dimension.
    ! The vertical axis is assumed to be in sigma coordinates.
    ! ----------------------------------------------------------------
    IF ( ASSOCIATED(lev) ) THEN

       ! -------------------------------------------------------------
       ! Initialize vertical axis. Set dependencies on other axis.
       ! Define dimensions in the following order: lev, lat, lon.
       ! (first dimension has to be the 'own' dimension!).
       ! -------------------------------------------------------------
       N        = N + 1
       ax(N)%lm = .false.     ! VERTICAL AXIS IS NON-MODULO AXIS

       ! Axis dimension
       XLEV = SIZE(lev,1)
       YLEV = SIZE(lev,2)
       ZLEV = SIZE(lev,3)

       ! Sanity check: if XLEV and YLEV are > 1, they must correspond
       ! to the lon/lat axis defined above
       IF ( XLEV > 1 .AND. XLEV /= (XLON-1) ) THEN
          CALL HCO_ERROR ( 'level lon has wrong dimension', RC )
          RETURN
       ENDIF
       IF ( YLEV > 1 .AND. YLEV /= (YLAT-1) ) THEN
          CALL HCO_ERROR ( 'level lat has wrong dimension', RC )
          RETURN
       ENDIF

       ! Set dependencies. First dimension must be vertical axis!
       ndp = 1
       IF ( YLEV > 1 ) ndp = ndp + 1
       IF ( XLEV > 1 ) ndp = ndp + 1

       ax(N)%ndp = ndp
       ALLOCATE(ax(N)%dep(ax(N)%ndp), STAT=status)
       IF ( status /= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate lev dependencies', RC )
          RETURN
       ENDIF

       ! The variable dat%n holds the number of axis that level depends
       ! upon, and dat%dim are the corresponding axis dimensions (lengths).
       ax(N)%dat%n = ndp
       ALLOCATE(ax(N)%dat%dim(ax(N)%dat%n), STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate lat dimensions', RC )
          RETURN
       ENDIF

       ! Set axis indices and lengths:
       ! - First dimension is vertical axis
       cnt                = 1
       ax(N)%dep(cnt)     = N
       ax(N)%dat%dim(cnt) = ZLEV
       nlev               = ZLEV  ! number of grid cells
       ! - Next dimension is latitude
       IF ( YLEV > 1 ) THEN
          cnt                = cnt + 1
          ax(N)%dep(cnt)     = ndep_lat
          ax(N)%dat%dim(cnt) = YLEV
          nlev               = nlev * YLEV
       ENDIF
       ! - Next dimension is longitude
       IF ( XLEV > 1 ) THEN
          cnt                = cnt + 1
          ax(N)%dep(cnt)     = ndep_lon
          ax(N)%dat%dim(cnt) = XLEV
          nlev               = nlev * XLEV
       ENDIF

       ALLOCATE(ax(N)%dat%vd(nlev),STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate lat axis', RC )
          RETURN
       ENDIF

       ! -------------------------------------------------------------
       ! Pass vertical sigma coordinates to vector
       ! -------------------------------------------------------------
       LOW = 1
       UPP = 1
       DO J = 1, YLEV !lat
       DO I = 1, XLEV !lon
          UPP                   = LOW + ZLEV - 1 ! Upper index
          ax(N)%dat%vd(LOW:UPP) = lev(I,J,:)     ! Pass to vector
          LOW                   = UPP + 1        ! Next lower index
       ENDDO
       ENDDO

    END IF !lev

    CALL HCO_LEAVE ( HcoState%config%Err, RC )

  END SUBROUTINE AXIS_CREATE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Axis_Delete
!
! !DESCRIPTION: Subroutine AXIS\_DELETE deletes the specified MESSy
! axis.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AXIS_DELETE ( ax1, ax2, RC )
!
! !USES:
!
  USE MESSY_NCREGRID_BASE,  ONLY : INIT_AXIS
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(axis),       POINTER            :: ax1(:)
    TYPE(axis),       POINTER            :: ax2(:)
    INTEGER,          INTENT(INOUT)      :: RC
!
! !REVISION HISTORY:
!  28 Aug 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ROUTINE ARGUMENTS:
!
    INTEGER     :: I, status

    !=================================================================
    ! AXIS_DELETE begins here
    !=================================================================

    DO I=1, SIZE(ax1)
       CALL INIT_AXIS(ax1(I))
       CALL INIT_AXIS(ax2(I))
    ENDDO
    DEALLOCATE(ax1, ax2, STAT=status)
    IF(status/=0) THEN
       CALL HCO_ERROR('AXIS DEALLOCATION ERROR', RC )
       RETURN
    ENDIF
    NULLIFY(ax1)
    NULLIFY(ax2)

  END SUBROUTINE AXIS_DELETE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Hco2Messy
!
! !DESCRIPTION: Subroutine HCO2MESSY converts a HEMCO data array into a
! messy array-structure.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO2MESSY( HcoState, InArr, narr, ax, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState
    REAL(sp),         POINTER        :: InArr(:,:,:,:)
    TYPE(narray),     POINTER        :: narr(:)
    TYPE(axis),       POINTER        :: ax(:)
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  27 Jun 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ROUTINE ARGUMENTS:
!
    INTEGER            :: NCELLS, NX, NY, NZ, NT
    INTEGER            :: TMP
    INTEGER            :: status
    INTEGER            :: J, L, T, LOW, UPP
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! HCO2MESSY begins here
    !=================================================================

    ! For error handling
    LOC = 'HCO2MESSY (HCOI_MESSY_MOD.F90)'

    ! ----------------------------------------------------------------
    ! Number of grid cells
    ! ----------------------------------------------------------------
    NX = SIZE(InArr,1)
    NY = SIZE(InArr,2)
    NZ = SIZE(InArr,3)
    NT = SIZE(InArr,4)
    NCELLS = NX * NY * NZ

    ! create
    ALLOCATE(narr(NT),STAT=status)
    IF(status/=0) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'narr allocation error', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Set MESSy vector dimensions. Copy from source axis object.
    ! Array dimensions are always one less than axis interface dimensions!
    ! For last ('free') dimension, set dimension length to 1.
    DO T = 1, NT
       narr(T)%n = size(ax)
       ALLOCATE(narr(T)%dim(narr(T)%n),STAT=status)
       IF(status/=0) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate array dims', RC, THISLOC=LOC )
          RETURN
       ENDIF

       DO J = 1, narr(T)%n
          IF ( ax(J)%dat%n > 0 ) THEN
             tmp = ax(J)%dat%dim(1)-1
          ELSE
             tmp = 1
          ENDIF
          narr(T)%dim(J) = tmp
       ENDDO

       ALLOCATE(narr(T)%vd(NCELLS),STAT=status)
       IF(status/=0) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate array', RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDDO !T

    ! ----------------------------------------------------------------
    ! Pass HEMCO data array in slices to MESSy data vector.
    ! The MESSy data vector is an 'unfolded' HEMCO data array, with
    ! the longitude dimension changing the fastest, then followed by
    ! latitude, and level. Hence, we pass the HEMCO array in slices
    ! along the longitude axis to the MESSy vector, starting with the
    ! slice at lat/lev position 1/1, followed by 2/1, 3/1, etc.
    ! For now, the MESSy vector is always of type double. In future,
    ! we may want to set it to the same type as the HEMCO array and
    ! use pointers to the HEMCO data slices!
    ! ----------------------------------------------------------------

    ! Vector index counter
    LOW = 1
    UPP = 1

    ! Loop over all higher dimensions. Don't loop over lon because we pass all
    ! lon-values in slices.
    DO L = 1, NZ  ! NZ is 1 for 2D arrays
    DO J = 1, NY
       UPP = LOW + NX - 1   ! Upper index of slice to be filled
       DO T = 1, NT
          narr(T)%vd(LOW:UPP) = InArr(:,J,L,T)   ! Pass this slice to vector
       ENDDO
       LOW = UPP + 1        ! Next slice begins right after this one
    ENDDO !J
    ENDDO !L

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO2MESSY
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Messy2Hco
!
! !DESCRIPTION: Subroutine MESSY2HCO converts a MESSy array structure
! into a HEMCO data array. This is basically the reverse function of
! HCO2MESSY.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MESSY2HCO( HcoState, narr, Lct, LEV, Ptr4D, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                :: HcoState
    TYPE(narray),     POINTER                :: narr(:)
    TYPE(ListCont),   POINTER                :: Lct
    INTEGER,          INTENT(IN   )          :: LEV
    REAL(sp),         POINTER                :: Ptr4D(:,:,:,:)
    INTEGER,          INTENT(INOUT)          :: RC
!
! !REVISION HISTORY:
!  27 Jun 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ROUTINE ARGUMENTS:
!
    INTEGER            :: J, L, T
    INTEGER            :: LOW, UPP
    INTEGER            :: NX, NY, NZ, NT, Z1, Z2
    CHARACTER(LEN=255) :: MSG

    !=================================================================
    ! MESSY2HCO begins here
    !=================================================================

    ! Grid dimensions
    NX = HcoState%NX
    NY = HcoState%NY
    NT = SIZE(narr)

    ! Vertical levels to be filled on output array. Only level LEV
    ! is filled if NC is different than one!
    IF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN
       Z1 = 1
       Z2 = 1
    ELSEIF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
       IF ( ASSOCIATED(Ptr4D) ) THEN
          Z1 = LEV
          Z2 = LEV
       ELSE
          Z1 = 1
          Z2 = SIZE(Lct%Dct%Dta%V3(1)%Val,3)
       ENDIF
    ENDIF

    ! Vector index counter
    LOW = 1
    UPP = 1

    ! If temporary array Ptr4D is provided, make sure that its
    ! dimensions are correct
    IF ( ASSOCIATED(Ptr4D) ) THEN
       NZ = Z2 - Z1 + 1
       IF ( ( SIZE(Ptr4D,1) /= NX ) .OR. &
            ( SIZE(Ptr4D,2) /= NY ) .OR. &
            ( SIZE(Ptr4D,3) /= NZ ) .OR. &
            ( SIZE(Ptr4D,4) /= NT )       ) THEN
          WRITE(MSG,*) 'Temporary pointer has wrong dimensions: ', &
                       TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, &
                           THISLOC='MESSY2HCO (hcoio_messy_mod.F90)' )
          RETURN
       ENDIF
    ENDIF

    ! Loop over all higher dimensions
    DO L = Z1, Z2
    DO J = 1, NY

       ! Upper index of slice to be filled
       UPP = LOW + NX - 1

       ! Do for every time slice
       DO T = 1, NT

          ! If provided, write into temporary array
          IF ( ASSOCIATED(Ptr4D) ) THEN
             Ptr4D(:,J,L,T) = narr(T)%vd(LOW:UPP)

          ! If temporary array Ptr4D is not provided, write directly
          ! into list container array
          ELSE

             ! 2D
             IF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN
                Lct%Dct%Dta%V2(T)%Val(:,J) = narr(T)%vd(LOW:UPP)
             ! 3D
             ELSE
                Lct%Dct%Dta%V3(T)%Val(:,J,L) = narr(T)%vd(LOW:UPP)
             ENDIF

          ENDIF
       ENDDO !NT

       ! Next slice begins right after this one
       LOW = UPP + 1
    ENDDO !J
    ENDDO !L

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE MESSY2HCO
!EOC
END MODULE HCOIO_MESSY_MOD
