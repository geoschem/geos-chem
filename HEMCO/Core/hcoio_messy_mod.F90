!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_messy_mod.F90 
!
! !DESCRIPTION: Module HCOIO\_MESSY\_MOD interfaces HEMCO with the regridding
! tool NCREGRID of the Modular Earth Submodel System (MESSy).
! Note that for now, this code only works for rectilinear (lon-lat)
! grids.
! TODO: insert vertical regridding capabilities. 
!\\
! !REFERENCES: 
! \begin{itemize}
! \item Joeckel, P. Technical note: Recursive rediscretisation of geo-
! scientific data in the Modular Earth Submodel System (MESSy), ACP, 6,
! 3557--3562, 2006.
! \end{itemize}
!\\
! !INTERFACE: 
!
MODULE HCOIO_MESSY_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD,        ONLY : Hco_State
  USE HCO_DATACONT_MOD,     ONLY : ListCont
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
! !IROUTINE: HCO_MESSY_REGRID
!
! !DESCRIPTION: This is the wrapper routine to regrid a 4D input array
! NcArr (x,y,z,t) onto the HEMCO emissions grid (defined in HcoState)
! using the regridding tool NCREGRID. LonEdge, LatEdge and LevEdge are 
! the grid point edges of the input grid. The data is written into list 
! container Lct.
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_MESSY_REGRID ( am_I_Root, HcoState, NcArr,   &
                                LonEdge,   LatEdge,  LevEdge, &
                                Lct,       RC                  ) 
!
! !USES:
!
  USE HCO_FILEDATA_MOD,     ONLY : FileData_ArrCheck
  USE HCO_UNIT_MOD,         ONLY : HCO_IsIndexData
  USE MESSY_NCREGRID_BASE,  ONLY : RG_INT, RG_IDX
  USE MESSY_NCREGRID_BASE,  ONLY : NREGRID
  USE MESSY_NCREGRID_BASE,  ONLY : INIT_NARRAY 
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root       ! Root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState        ! HEMCO obj.
    REAL(sp),         POINTER        :: ncArr(:,:,:,:)  ! Input array(x,y,z,t)
    REAL(hp),         POINTER        :: LonEdge(:)      ! lon edges
    REAL(hp),         POINTER        :: LatEdge(:)      ! lat edges
    REAL(hp),         POINTER        :: LevEdge(:,:,:)  ! sigma pressure edges
    TYPE(ListCont),   POINTER        :: Lct             ! Target list container
    INTEGER,          INTENT(INOUT)  :: RC              ! Return code
!
! !REVISION HISTORY:
!  27 Jun 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
    TYPE(narray), POINTER   :: narr_src(:)    => NULL()
    TYPE(narray), POINTER   :: narr_dst(:)    => NULL()
    TYPE(axis),   POINTER   :: axis_src(:)    => NULL()
    TYPE(axis),   POINTER   :: axis_dst(:)    => NULL()
    INTEGER,      POINTER   :: rg_type(:)     => NULL()
    INTEGER,      POINTER   :: rcnt   (:)     => NULL()
    REAL(dp),     POINTER   :: sovl   (:)     => NULL()
    REAL(dp),     POINTER   :: dovl   (:)     => NULL()
    REAL(hp),     POINTER   :: lon    (:)     => NULL()
    REAL(hp),     POINTER   :: lat    (:)     => NULL()
    REAL(hp),     POINTER   :: lev    (:,:,:) => NULL()
    INTEGER                 :: NLEV, NTIME
    INTEGER                 :: I, L, status
    CHARACTER(LEN=255)      :: MSG, LOC
    LOGICAL                 :: SameGrid

    ! testing only (avoid vertical regridding)
    real(sp), pointer :: tmparr(:,:,:,:) => NULL()

    !=================================================================
    ! HCO_MESSY_REGRID begins here
    !=================================================================

    ! For error handling
    LOC = 'HCO_MESSY_REGRID (HCOI_MESSY_MOD.F90)'
    CALL HCO_ENTER ( LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Number of vertical levels of input data
    NLEV = SIZE(NcArr,3)

    ! Number of time slices. All time slices will be regridded 
    ! simultaneously.
    NTIME = SIZE(NcArr,4)
    
    ! Error check: data must be 2D or 3D.
    IF ( Lct%Dct%Dta%SpaceDim /= 2 .AND. &
         Lct%Dct%Dta%SpaceDim /= 3         ) THEN
       MSG = 'Can only regrid 2D or 3D data: ' // TRIM(Lct%Dct%cName)
       CALL HCO_ERROR ( MSG, RC )
       RETURN 
    ENDIF

    !-----------------------------------------------------------------
    ! Make sure output array is defined & allocated
    !-----------------------------------------------------------------
    IF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN
       CALL FileData_ArrCheck( Lct%Dct%Dta, HcoState%NX, HcoState%NY, &
                               NTIME,       RC                         )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ELSEIF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
       CALL FileData_ArrCheck( Lct%Dct%Dta, HcoState%NX, HcoState%NY, &
                               HcoState%NZ, NTIME,       RC            ) 
       IF ( RC /= HCO_SUCCESS ) RETURN   
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
    IF ( SIZE(NcArr,1) == HcoState%NX .AND. &
         SIZE(NcArr,2) == HcoState%NY        ) THEN

       ! Vertical dimensions have to match or be 1
       IF ( NLEV == 1 .OR. NLEV == HcoState%NZ ) THEN

          ! Assume same grid
          SameGrid = .TRUE.

          ! Check for same boundaries. Otherwise (re-) falsify SameGrid
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
          IF ( NLEV > 1 ) THEN
             ! TODO: Eventually need to add level check here. 
             ! For now, assume that input levels are equal to 
             ! output level if they have the same dimension!
          ENDIF
       ENDIF
    ENDIF

    ! If grids match, pass data to list container.
    IF ( SameGrid ) THEN
       MSG = 'Input grid seems to match output grid. ' // &
             'No regridding is performed: ' // TRIM(Lct%Dct%cName)
       CALL HCO_WARNING( MSG, RC )

       ! For every time slice...
       DO I = 1, NTIME
       DO L = 1, NLEV

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
       CALL HCO_LEAVE ( RC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Source grid description.
    ! This creates a MESSy axis object for the source grid.
    !-----------------------------------------------------------------
    lon => LonEdge
    lat => LatEdge
    ! vertical regridding not supported for now!
!    lev => LevEdge

    CALL AXIS_CREATE( am_I_Root, HcoState, lon, lat, lev, axis_src, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Free pointer
    lon => NULL()
    lat => NULL()
    lev => NULL()

    !-----------------------------------------------------------------
    ! Destination grid description.
    ! This creates a MESSy axis object for the target (=HEMCO) grid.
    !-----------------------------------------------------------------

    lon => HcoState%Grid%XEDGE%Val(:,1)
    lat => HcoState%Grid%YEDGE%Val(1,:)
    ! vertical regridding not supported for now!
!    IF ( ASSOCIATED(LevEdge) ) lev => HcoState%Grid%ZSIGMA%Val(:,:,:)

    CALL AXIS_CREATE( am_I_Root, HcoState, lon, lat, lev, axis_dst, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Free pointer
    lon => NULL()
    lat => NULL()
    lev => NULL()

    !-----------------------------------------------------------------
    ! Set all other regridding parameter
    !-----------------------------------------------------------------

    ! rg_type denotes the regridding type for each array (i.e. time 
    ! slice). Set to 'intensive quantity' for all concentrations (incl.
    ! unitless) data. Set this to 'index distribution' for data marked
    ! as index data in the configuration file.
    ALLOCATE(rg_type(NTIME))
    IF ( HCO_IsIndexData(Lct%Dct%Dta%OrigUnit) ) THEN
       rg_type(:) = RG_IDX
    ELSE
       rg_type(:) = RG_INT 
    ENDIF

    ! temporary level loop to avoid vertical regridding...
    ! TODO: remove
    if ( associated(levedge) ) then
       nlev = size(levedge,3) - 1
    else
       nlev = 1
    endif
    do l = 1,nlev
       tmpArr => ncArr(:,:,l:l,:)

    !-----------------------------------------------------------------
    ! Map input array onto MESSy array. Different time slices are
    ! stored as individual vector elements of narr_src.
    !-----------------------------------------------------------------
!    CALL HCO2MESSY( am_I_Root, NcArr, narr_src, axis_src, RC )
    CALL HCO2MESSY( am_I_Root, tmpArr, narr_src, axis_src, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Do the regridding
    !-----------------------------------------------------------------
    CALL NREGRID(s=narr_src,      sax=axis_src, dax=axis_dst, d=narr_dst, &
                 rg_type=rg_type, sovl=sovl,    dovl=dovl,    rcnt=rcnt    )

    !-----------------------------------------------------------------
    ! Map the destination array narr_dst onto the data vector in the 
    ! HEMCO list container.
    !-----------------------------------------------------------------
    ! TODO: remove level index l
    CALL MESSY2HCO( am_I_Root, HcoState, narr_dst, Lct, l, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN      
      
    !-----------------------------------------------------------------
    ! Cleanup
    !-----------------------------------------------------------------
    DEALLOCATE( sovl, dovl, rcnt, STAT=status)
    IF(status/=0) THEN
       CALL HCO_ERROR('DEALLOCATION ERROR 1', RC )
       RETURN
    ENDIF
    NULLIFY(sovl, dovl, rcnt)
    
    DO I=1, SIZE(narr_dst)
       CALL INIT_NARRAY(narr_dst(I))
    ENDDO
    DEALLOCATE(narr_dst, STAT=status)
    IF(status/=0) THEN
       CALL HCO_ERROR('DEALLOCATION ERROR 3', RC )
       RETURN
    ENDIF
    NULLIFY(narr_dst)
    
    DO I=1, SIZE(narr_src)
       CALL INIT_NARRAY(narr_src(I))
    ENDDO    
    DEALLOCATE(narr_src, STAT=status)
    IF(status/=0) THEN
       CALL HCO_ERROR('DEALLOCATION ERROR 2', RC )
       RETURN
    ENDIF
    NULLIFY(narr_src)

    ! end temporary loop
    enddo
    
    DEALLOCATE( rg_type )
    rg_type => NULL()

    CALL AXIS_DELETE( axis_src, axis_dst, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCO_MESSY_REGRID
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AXIS_CREATE
!
! !DESCRIPTION: Subroutine AXIS\_CREATE creates a MESSy axis type 
! from the grid defined by mid points Lon, Lat, Lev.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AXIS_CREATE ( am_I_Root, HcoState, lon, lat, lev, ax, RC ) 
!
! !USES:
!
  USE MESSY_NCREGRID_BASE,  ONLY : INIT_AXIS 
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root
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
    CALL HCO_ENTER ( LOC, RC )
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
       CALL HCO_ERROR ( 'Cannot allocate axis', RC )
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
          CALL HCO_ERROR ( 'Cannot allocate lon dependencies', RC )
          RETURN
       ENDIF
       ax(N)%dep(1) = N          ! ... INDEPENDENT
       ndep_lon     = N
       
       ax(N)%dat%n = 1          ! 1 dimension
       ALLOCATE(ax(N)%dat%dim(ax(N)%dat%n), STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate lon dimensions', RC )
          RETURN
       ENDIF
       ax(N)%dat%dim(:) = XLON
    
       ALLOCATE(ax(N)%dat%vd(XLON),STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate lon axis', RC )
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
          CALL HCO_ERROR ( 'Cannot allocate lat dependencies', RC )
          RETURN
       ENDIF
       ax(N)%dep(1) = N          ! ... INDEPENDENT
       ndep_lat     = N
      
       ax(N)%dat%n = 1          ! 1 dimension
       ALLOCATE(ax(N)%dat%dim(ax(N)%dat%n), STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate lat dimensions', RC )
          RETURN
       ENDIF
       ax(N)%dat%dim(:) = YLAT

       ALLOCATE(ax(N)%dat%vd(YLAT),STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate lat axis', RC )
          RETURN
       ENDIF
       ax(N)%dat%vd(:) = lat 
 
       ! TAKE INTO ACCOUNT SPHERICAL GEOMETRY ...
       ax(N)%dat%vd = COS( ( (ax(N)%dat%vd - 90.0_dp) / 180.0_dp ) * &
                           HcoState%Phys%PI )

    ENDIF !lat
 
    ! ----------------------------------------------------------------
    ! Assign vertical levels (if defined): this is the 3rd dimension.
    ! The vertical axis is always in hybrid sigma coordinates.
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
       IF ( XLEV > 1 .AND. XLEV /= XLON ) THEN
          CALL HCO_ERROR ( 'level lon has wrong dimension', RC )
          RETURN
       ENDIF
       IF ( YLEV > 1 .AND. YLEV /= YLAT ) THEN
          CALL HCO_ERROR ( 'level lat has wrong dimension', RC )
          RETURN
       ENDIF

       ! Set dependencies. First dimension must be vertical axis!
       ndp = 1
       IF ( YLEV > 1 ) ndp = ndp + 1
       IF ( XLEV > 1 ) ndp = ndp + 1

       ax(N)%ndp = ndp 
       ALLOCATE(ax(N)%dep(ax(N)%ndp), STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate lev dependencies', RC )
          RETURN
       ENDIF
    
       ! The variable dat%n holds the number of axis level depends upon, 
       ! and dat%dim are the corresponding axis dimensions (lengths).
       ax(N)%dat%n = ndp
       ALLOCATE(ax(N)%dat%dim(ax(N)%dat%n), STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate lat dimensions', RC )
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
          ax(N)%dat%dim(cnt) = YLAT
          nlev               = nlev * YLAT
       ENDIF
       ! - Next dimension is longitude
       IF ( XLEV > 1 ) THEN
          cnt                = cnt + 1
          ax(N)%dep(cnt)     = ndep_lon
          ax(N)%dat%dim(cnt) = XLON
          nlev               = nlev * XLON
       ENDIF

       ALLOCATE(ax(N)%dat%vd(nlev),STAT=status)
       IF ( status/= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate lat axis', RC )
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

    CALL HCO_LEAVE ( RC )

  END SUBROUTINE AXIS_CREATE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AXIS_DELETE
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
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO2MESSY
!
! !DESCRIPTION: Subroutine HCO2MESSY converts a HEMCO data array into a 
! messy array-structure.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO2MESSY( am_I_Root, InArr, narr, ax, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
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
       CALL HCO_ERROR( 'narr allocation error', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Set MESSy vector dimensions. Copy from source axis object.
    ! Array dimensions are always one less than axis interface dimensions!
    ! For last ('free') dimension, set dimension length to 1.
    DO T = 1, NT
       narr(T)%n = size(ax)
       ALLOCATE(narr(T)%dim(narr(T)%n),STAT=status)
       IF(status/=0) THEN
          CALL HCO_ERROR( 'Cannot allocate array dims', RC, THISLOC=LOC )
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
          CALL HCO_ERROR( 'Cannot allocate array', RC, THISLOC=LOC )
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
! !IROUTINE: MESSY2HCO
!
! !DESCRIPTION: Subroutine MESSY2HCO converts a MESSy array structure 
! into a HEMCO data array. This is basically the reverse function of
! HCO2MESSY.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MESSY2HCO( am_I_Root, HcoState, narr, Lct, lev, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
    TYPE(HCO_State),  POINTER        :: HcoState
    TYPE(narray),     POINTER        :: narr(:)
    TYPE(ListCont),   POINTER        :: Lct 
    INTEGER,          INTENT(IN   )  :: lev ! temporary level index
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
    INTEGER        :: J, L, T
    INTEGER        :: LOW, UPP
    INTEGER        :: NX, NY, NZ, NT

    !=================================================================
    ! MESSY2HCO begins here
    !=================================================================

    ! Grid dimensions
    NX = HcoState%NX
    NY = HcoState%NY 
    IF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN
       NZ = 1
    ELSEIF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
       NZ = HcoState%NZ
    ENDIF
    NT = SIZE(narr)

    ! Vector index counter
    LOW = 1
    UPP = 1

    ! Loop over all higher dimensions
    ! TODO: re-establish loop over vertical dimension
    DO L = 1, 1
!    DO L = 1, NZ  ! NZ is 1 for 2D arrays
    DO J = 1, NY
    
       ! Upper index of slice to be filled
       UPP = LOW + NX - 1
      
       ! Do for every time slice
       DO T = 1, NT
          ! 2D
          IF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN
             Lct%Dct%Dta%V2(T)%Val(:,J) = narr(T)%vd(LOW:UPP)
          ! 3D
          ELSE   
             Lct%Dct%Dta%V3(T)%Val(:,J,LEV) = narr(T)%vd(LOW:UPP)
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
