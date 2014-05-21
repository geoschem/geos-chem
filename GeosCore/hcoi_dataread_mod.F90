!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_dataread_mod.F90 
!
! !DESCRIPTION: Module HCO\_DATAREAD\_MOD controls data processing (file
! reading, unit conversion, regridding) for HEMCO.
!\\
!\\
! !INTERFACE: 
!
MODULE HCOI_DATAREAD_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD,       ONLY : Hco_State
  USE HCO_DATACONT_MOD,    ONLY : ListCont

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOI_DATAREAD
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !MODULE INTERFACES:
!
! !REVISION HISTORY:
!  22 Aug 2013 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
#if defined(ESMF_)
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_DATAREAD 
!
! !DESCRIPTION: Interface routine between ESMF and HEMCO to obtain
! the data array for a given HEMCO data container. 
!
! NOTE/TODO: For now, all arrays are copied into the HEMCO data array. 
! We may directly point to the ESMF arrays in future. In this case, we
! may have to force the target ID to be equal to the container ID (the
! target ID is set in hco\_config\_mod).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOI_DATAREAD ( am_I_Root, HcoState, Lct, RC ) 
!
! !USES:
!
      USE HCO_FILEDATA_MOD, ONLY : FileData_ArrCheck2D
      USE HCO_FILEDATA_MOD, ONLY : FileData_ArrCheck3D
      USE ESMF_MOD
      USE MAPL_MOD
# include "MAPL_Generic.h"
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State),  POINTER        :: HcoState
      TYPE(ListCont),   POINTER        :: Lct 
      INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  28 Aug 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      INTEGER                    :: II, JJ, LL, TT
      INTEGER                    :: I, J, L, T
      INTEGER                    :: STAT
      REAL,             POINTER  :: Ptr3D(:,:,:)   => NULL() 
      REAL,             POINTER  :: Ptr2D(:,:)     => NULL() 
      TYPE(ESMF_State), POINTER  :: IMPORT         => NULL()
      CHARACTER(LEN=255)         :: LOC

      !=================================================================
      ! HCOI_DATAREAD begins here
      !=================================================================

      ! For error handling
      LOC = 'HCOI_DATAREAD (HCOI_DATAREAD_MOD.F90)'
      CALL HCO_ENTER ( LOC, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Point to ESMF IMPORT object
      IMPORT => HcoState%IMPORT

      !-----------------------------------------------------------------
      ! Read 3D data from ESMF 
      !-----------------------------------------------------------------
      IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN

         ! Get data
         CALL MAPL_GetPointer ( IMPORT, Ptr3D, &
                                TRIM(Lct%Dct%cName), RC=STAT )

         ! Check for MAPL error
         IF( MAPL_VRFY(STAT,LOC,1) ) THEN
            CALL HCO_ERROR ( 'Cannot get xyz pointer', RC ) 
            RETURN
         ENDIF

         ! Get array dimensions 
         II = SIZE(Ptr3D,1)
         JJ = SIZE(Ptr3D,2) 
         LL = SIZE(Ptr3D,3)
         TT = 1 

         ! Allocate HEMCO array if not yet defined
         CALL FileData_ArrCheck3D( Lct%Dct%Dta, II, JJ, LL, TT, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! Copy data and cast to real*8
         Lct%Dct%Dta%V3(1)%Val(:,:,:) = Ptr3D(:,:,:)

      !-----------------------------------------------------------------
      ! Read 2D data from ESMF 
      !-----------------------------------------------------------------
      ELSEIF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN

         ! Get data
         CALL MAPL_GetPointer ( IMPORT, Ptr2D, &
                                TRIM(Lct%Dct%cName), RC=STAT )


         ! Check for MAPL error 
         IF( MAPL_VRFY(STAT,LOC,2) ) THEN
            CALL HCO_ERROR ( 'Cannot get xy pointer', RC ) 
            RETURN
         ENDIF

         ! Get array dimensions 
         II = SIZE(Ptr2D,1)
         JJ = SIZE(Ptr2D,2) 
         LL = 1 
         TT = 1 

         ! Allocate HEMCO array if not yet defined
         CALL FileData_ArrCheck2D( Lct%Dct%Dta, II, JJ, TT, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! Copy to HEMCO container
         Lct%Dct%Dta%V2(1)%Val(:,:) = Ptr2D(:,:)

      ENDIF  
 
      !-----------------------------------------------------------------
      ! Cleanup and leave 
      !-----------------------------------------------------------------
      Ptr3D  => NULL()
      Ptr2D  => NULL()
      IMPORT => NULL()   

      ! Return w/ success
      CALL HCO_LEAVE ( RC )

      END SUBROUTINE HCOI_DATAREAD
!EOC
#else
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_DATAREAD 
!
! !DESCRIPTION: Reads a netCDF file and returns the regridded array in proper
! units. This is an interface to the GEOS-Chem reading and regridding 
! routines. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOI_DATAREAD ( am_I_Root, HcoState, Lct, RC ) 
!
! !USES:
!
      USE NCDF_MOD,           ONLY : NC_OPEN,        NC_CLOSE
      USE NCDF_MOD,           ONLY : NC_READ_VAR,    NC_READ_ARR
      USE NCDF_MOD,           ONLY : NC_READ_GRID
      USE HCO_UNIT_MOD,       ONLY : HCO_UNIT_CHANGE, HCO_UNIT_SCALCHECK
      USE REGRID_A2A_MOD,     ONLY : MAP_A2A
      USE HCO_FILEDATA_MOD,   ONLY : FileData_ArrCheck2D
      USE HCO_FILEDATA_MOD,   ONLY : FileData_ArrCheck3D
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State),  POINTER        :: HcoState
      TYPE(ListCont),   POINTER        :: Lct
      INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      CHARACTER(LEN=255)    :: lonUnit, latUnit, levUnit, ncUnit
      CHARACTER(LEN=255)    :: MSG 
      INTEGER               :: I, L, T
      INTEGER               :: NX, NY
      INTEGER               :: NCRC, Flag
      INTEGER               :: ncLun
      INTEGER               :: nLon,   nLat,  nLev, nTime
      INTEGER               :: lev1,   lev2 
      INTEGER               :: tidx1,  tidx2, ncYr, ncMt
      INTEGER               :: HcoID
      REAL(sp), POINTER     :: ncArr(:,:,:,:) => NULL()
      REAL(sp), POINTER     :: ORIG_2D(:,:)   => NULL()
      REAL(hp), POINTER     :: REGR_2D(:,:)   => NULL()
      REAL(sp), ALLOCATABLE :: XEDGE_IN (:)
      REAL(sp), ALLOCATABLE :: YSIN_IN  (:)
      REAL(sp)              :: XEDGE_OUT(HcoState%NX+1) 
      REAL(sp)              :: YSIN_OUT (HcoState%NY+1)
      LOGICAL               :: verb, forcescal, IsPerArea

      !=================================================================
      ! HCOI_DATAREAD begins here
      !=================================================================

      ! Enter
      CALL HCO_ENTER ('HCOI_DATAREAD (HCOI_DATAREAD_MOD.F90)' , RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      
      ! Check for verbose mode and if scale factors have to be in scale
      ! factor units.
      verb      = HCO_VERBOSE_CHECK()
      forcescal = HCO_FORCESCAL_CHECK()
 
      ! Copy horizontal grid dimensions from HEMCO state object
      NX = HcoState%NX
      NY = HcoState%NY

      ! Verbose mode
      IF ( verb ) THEN
         Write(MSG,*) 'Reading file ', TRIM(Lct%Dct%Dta%ncFile)
         CALL HCO_MSG(MSG,SEP1='-')
      ENDIF 

      ! ----------------------------------------------------------------
      ! Open netCDF
      ! ----------------------------------------------------------------
      CALL NC_OPEN ( TRIM(Lct%Dct%Dta%ncFile), ncLun )

      ! ----------------------------------------------------------------
      ! Extract time slice information
      ! This determines the lower and upper time slice index (tidx1 
      ! and tidx2) to be read based upon the time slice information 
      ! extracted from the file and the time stamp settings set in the
      ! HEMCO configuration file.
      ! ----------------------------------------------------------------
      CALL GET_TIMEIDX ( am_I_Root, HcoState, Lct, &
                         ncLun,     tidx1,    tidx2,   &
                         ncYr,      ncMt,     RC        )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! ----------------------------------------------------------------
      ! Read data 
      ! ----------------------------------------------------------------

      ! Extract grid dimension
      CALL NC_READ_VAR ( ncLun, 'lat', nLat, LatUnit, RC=NCRC )
      IF ( NCRC /= 0 .OR. nLat == 0 ) THEN
         CALL HCO_ERROR( 'NC_READ_LAT', RC )
         RETURN 
      ENDIF

      CALL NC_READ_VAR ( ncLun, 'lon', nLon, LonUnit, RC=NCRC )
      IF ( NCRC /= 0 .OR. nLon ) THEN
         CALL HCO_ERROR( 'NC_READ_LON', RC )
         RETURN 
      ENDIF

      CALL NC_READ_VAR ( ncLun, 'lev', nLev, LevUnit, RC=NCRC )
      IF ( NCRC /= 0 ) THEN
         CALL HCO_ERROR( 'NC_READ_LEV', RC )
         RETURN 
      ENDIF
      IF ( nLev == 0 ) THEN
         CALL NC_READ_VAR ( ncLun, 'height', nLev, LevUnit, RC=NCRC )
         IF ( NCRC /= 0 ) THEN
            CALL HCO_ERROR( 'NC_READ_LEV', RC )
            RETURN 
         ENDIF
      ENDIF 

      ! Sanity check of lon/lat units
      IF ( INDEX( LonUnit, 'degrees_east' ) == 0 ) THEN
         MSG = 'illegal longitude unit in ' // &
               TRIM(Lct%Dct%Dta%ncFile)
         CALL HCO_ERROR ( MSG, RC )
         RETURN
      ENDIF
      IF ( INDEX( LatUnit, 'degrees_north' ) == 0 ) THEN
         MSG = 'illegal latitude unit in ' // &
               TRIM(Lct%Dct%Dta%ncFile)
         CALL HCO_ERROR ( MSG, RC )
         RETURN
      ENDIF

      ! Sanity check of vertical dimensions
      ! ==> numbers of vertical levels must not exceed HEMCO state
      ! levels (only horizontal regridding supported so far!)
      ! Also check if dimensionality agrees with settings in input
      ! file!
      IF ( nLev > HcoState%NZ ) THEN
         MSG = 'Too many vert. levels in ' // TRIM(Lct%Dct%Dta%ncFile)
         CALL HCO_ERROR ( MSG, RC )
         RETURN 
      ENDIF

      ! Set level indeces to be read
      ! NOTE: for now, always read all levels. Edit here to read
      ! only particular levels.
      ! Also do sanity check whether or not vertical dimension does 
      ! agree with space dimension specified in configuration file
      IF ( nLev > 0 ) THEN
         lev1 = 1
         lev2 = nlev
         IF ( Lct%Dct%Dta%SpaceDim <= 2 ) THEN
            MSG = 'Found vertical axis but dimension is set to xy - ' // &
                  'will only use first level: ' // TRIM(Lct%Dct%Dta%ncFile)
            CALL HCO_WARNING ( MSG, RC )
            Lct%Dct%Dta%SpaceDim = 3
            lev2 = 1
         ENDIF
      ELSE
         lev1 = 0
         lev2 = 0
         IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
            MSG = 'This is not 3D data: ' // TRIM(Lct%Dct%Dta%ncFile)
            CALL HCO_WARNING ( MSG, RC )
            Lct%Dct%Dta%SpaceDim = 2
         ENDIF
      ENDIF

      ! Read array
      CALL NC_READ_ARR ( fID     = ncLun,              &
                         ncVar   = Lct%Dct%Dta%ncPara, &
                         lon1    = 1,                  &
                         lon2    = nlon,               &
                         lat1    = 1,                  &
                         lat2    = nlat,               &
                         lev1    = lev1,               &
                         lev2    = lev2,               &
                         time1   = tidx1,              &
                         time2   = tidx2,              &
                         ncArr   = ncArr,              &
                         varUnit = ncUnit,             &
                         RC      = NCRC                 )
      IF ( NCRC /= 0 ) THEN
         CALL HCO_ERROR( 'NC_READ_ARRAY', RC )
         RETURN 
      ENDIF

      ! ----------------------------------------------------------------
      ! Close netCDF
      ! ----------------------------------------------------------------
      CALL NC_CLOSE ( ncLun )

      !-----------------------------------------------------------------
      ! Convert to HEMCO units 
      !-----------------------------------------------------------------

      ! Convert to HEMCO units. This is kg/m2/s for fluxes and kg/m3 
      ! for concentrations. Ignore this if no species ID defined,
      ! i.e. for scale factors.
      HcoID = Lct%Dct%HcoID 
      IF ( HcoID > 0 ) THEN

         CALL HCO_UNIT_CHANGE(                                       &
                 Array         = ncArr,                              &
                 Units         = ncUnit,                             &
                 MW_IN         = HcoState%Spc(HcoID)%MW_g,           & 
                 MW_OUT        = HcoState%Spc(HcoID)%EmMW_g,         & 
                 MOLEC_RATIO   = HcoState%Spc(HcoID)%MolecRatio,     & 
                 YYYY          = ncYr,                               &
                 MM            = ncMt,                               &
                 IsPerArea     = IsPerArea,                          &
                 RC            = NCRC                                 )

         IF ( NCRC /= 0 ) THEN
             CALL HCO_ERROR('CHANGE_UNITS', RC )
             RETURN 
         ENDIF

         IF ( .NOT. IsPerArea ) THEN
            CALL NORMALIZE_AREA( HcoState, ncArr, Lct%Dct%Dta%ncFile, RC )
            IF ( RC /= HCO_SUCCESS ) RETURN 
         ENDIF
     
      ! For scale factors and masks, check if units are indeed unitless. 
      ELSE
         Flag = HCO_UNIT_SCALCHECK( ncUnit )
         MSG = 'Scale factor does not appear to be unitless: ' // & 
               TRIM(Lct%Dct%Dta%ncFile)
         IF ( Flag == 1 ) THEN
            CALL HCO_WARNING( MSG, RC )
         ELSEIF ( Flag == 2 ) THEN
            IF ( ForceScal ) THEN
               MSG = 'Illegal scale factor unit: ' // TRIM(Lct%Dct%Dta%ncFile)
               CALL HCO_ERROR( MSG, RC )
               RETURN
            ELSE
               CALL HCO_WARNING( MSG, RC )
            ENDIF
         ENDIF

      ENDIF

      !-----------------------------------------------------------------
      ! Regrid onto emissions grid 
      !-----------------------------------------------------------------

      ! Get input grid edges from netCDF file 
      ! Note: use error syntax from ncdf_mod!
      ALLOCATE( XEDGE_IN(NLON+1) )
      ALLOCATE( YSIN_IN (NLAT+1) )
      XEDGE_IN(:) = 0.0_sp
      YSIN_IN (:) = 0.0_sp
      CALL NC_READ_GRID( nLon,     nLat,    Lct%Dct%Dta%ncFile, &
                         XEDGE_IN, YSIN_IN, NCRC                     )
      IF ( NCRC /= 0 ) THEN
         CALL HCO_ERROR ( 'NC_READ_GRID', RC )
         RETURN 
      ENDIF

      ! Get output grid edges from HEMCO state
      XEDGE_OUT(:) = HcoState%Grid%XEDGE(:,1)
      YSIN_OUT (:) = HcoState%Grid%YSIN (1,:) 
  
      ! Reset nlev and ntime to effective array sizes
      nlev  = size(ncArr,3)
      ntime = size(ncArr,4)

      ! Allocate output array if not yet defined
      IF ( Lct%Dct%Dta%SpaceDim <= 2 ) THEN
         CALL FileData_ArrCheck2D( Lct%Dct%Dta, nx, ny, ntime, RC ) 
         IF ( RC /= 0 ) RETURN
      ELSE
         CALL FileData_ArrCheck3D( Lct%Dct%Dta, nx, ny, nlev, ntime, RC ) 
         IF ( RC /= 0 ) RETURN
      ENDIF

      ! Do regridding
      DO T = 1, NTIME
      DO L = 1, NLEV 

         ! Point to 2D slices to be regridded
         ORIG_2D => ncArr(:,:,L,T)
         IF ( Lct%Dct%Dta%SpaceDim <= 2 ) THEN
            REGR_2D => Lct%Dct%Dta%V2(T)%Val(:,:)
         ELSE
            REGR_2D => Lct%Dct%Dta%V3(T)%Val(:,:,L)
         ENDIF

         ! Do the regridding
         CALL MAP_A2A( NLON,  NLAT,  XEDGE_IN,  YSIN_IN,  ORIG_2D, &
                       NX, NY, XEDGE_OUT, YSIN_OUT, REGR_2D, 0, 0 )

         ORIG_2D => NULL()
         REGR_2D => NULL()

      ENDDO !L
      ENDDO !T

      !-----------------------------------------------------------------
      ! Cleanup and leave 
      !-----------------------------------------------------------------
      IF ( ASSOCIATED(ncArr    ) ) DEALLOCATE ( ncArr    )
      IF ( ALLOCATED (XEDGE_IN ) ) DEALLOCATE ( XEDGE_IN )
      IF ( ALLOCATED (YSIN_IN  ) ) DEALLOCATE ( YSIN_IN  )

      ! Return w/ success
      CALL HCO_LEAVE ( RC ) 

      END SUBROUTINE HCOI_DATAREAD
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_TIMEIDX
!
! !DESCRIPTION: Returns the lower and upper time slice index (tidx1
! and tidx2, respectively) to be read. These values are determined 
! based upon the time slice information extracted from the netCDF file, 
! the time stamp settings set in the config. file, and the current 
! simulation date.
! Also return the time slice year and month, as these values may be
! used for unit conversion! 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_TIMEIDX ( am_I_Root, HcoState, Lct,     &
                               ncLun,     tidx1,    tidx2,   &
                               ncYr,      ncMt,     RC        )
!
! !USES:
!
      USE NCDF_MOD,           ONLY : NC_READ_TIME_YYYYMMDDhh
      USE HCO_TIDX_MOD,       ONLY : HCO_GetPrefTimeAttr
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root ! Root CPU?
      TYPE(HCO_State),  POINTER        :: HcoState  ! HcoState object
      TYPE(ListCont),   POINTER        :: Lct   ! List container
      INTEGER,          INTENT(IN   )  :: ncLun     ! open ncLun
      INTEGER,          INTENT(  OUT)  :: tidx1     ! lower time idx
      INTEGER,          INTENT(  OUT)  :: tidx2     ! upper time idx
      INTEGER,          INTENT(  OUT)  :: ncYr      ! time slice year
      INTEGER,          INTENT(  OUT)  :: ncMt      ! time slice month
      INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      CHARACTER(LEN=255)    :: MSG
      CHARACTER(LEN=1023)   :: MSG_LONG
      INTEGER               :: nTime,  T, NCRC 
      INTEGER               :: prefYr, prefMt, prefDy, prefHr
      INTEGER               :: refYear
      INTEGER               :: prefYMDh
      INTEGER, POINTER      :: availYMDh(:) => NULL() 
      LOGICAL               :: verb

      !=================================================================
      ! GET_TIMEIDX begins here
      !=================================================================

      ! Init 
      CALL HCO_ENTER ('GET_TIMEIDX (HCO_DATAREAD_MOD.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      verb = HCO_VERBOSE_CHECK() 
 
      ! ---------------------------------------------------------------- 
      ! Extract netCDF time slices (YYYYMMDDhh) 
      ! ----------------------------------------------------------------
      CALL NC_READ_TIME_YYYYMMDDhh ( ncLun, nTime,    availYMDH, &
                                     refYear=refYear, RC=NCRC     )     
      IF ( NCRC /= 0 ) THEN
         CALL HCO_ERROR( 'NC_READ_TIME_YYYYMMDDhh', RC )
         RETURN 
      ENDIF

      ! Return warning if reference year prior to 1801: it seems like
      ! the time slices may be off by one day!
      IF ( refYear <= 1900 ) THEN
         msg = 'ncdf reference year is prior to 1901 - ' // &
               'time stamps may be wrong!'
         CALL HCO_WARNING ( MSG, RC )
      ENDIF

      ! verbose mode 
      IF ( verb ) THEN
         write(MSG,'(A30,I12)') '# time slices read: ', nTime
         CALL HCO_MSG(MSG)
         IF ( nTime > 0 ) THEN
            write(MSG,'(A30,I12,I12)') '# time slice range: ', &
                                       availYMDH(1), availYMDH(nTime) 
            CALL HCO_MSG(MSG)
         ENDIF
      ENDIF

      ! ---------------------------------------------------------------- 
      ! Select time slices to read
      ! ---------------------------------------------------------------- 

      ! ================================================================
      ! Case 1: Only one time slice available. 
      ! ================================================================
      IF ( nTime == 1 ) THEN
         tidx1 = 1
         tidx2 = 1

      ! ================================================================
      ! Case 2: More than one time slice available. Determine lower 
      ! and upper time slice index from file & HEMCO settings. 
      ! ================================================================
      ELSEIF ( nTime > 1 ) THEN

         ! Init
         tidx1 = -1
         tidx2 = -1 

         ! ------------------------------------------------------------- 
         ! Get preferred time stamp to read based upon the specs set
         ! in the config. file. 
         ! This can return value -1 for prefHr, indicating that all 
         ! hourly slices shall be read. Will always return a valid 
         ! attribute for Yr, Mt, and Dy.
         ! ------------------------------------------------------------- 
         CALL HCO_GetPrefTimeAttr ( Lct,    prefYr, prefMt, &
                                    prefDy, prefHr, RC       )
         IF ( RC /= HCO_SUCCESS ) RETURN
         prefYMDh = prefYr*1000000 + prefMt*10000 + &
                    prefDy*100     + max(prefHr,0)

         ! verbose mode
         IF ( verb ) THEN
            write(MSG,'(A30,I12)') 'preferred datetime: ', prefYMDh
            CALL HCO_MSG(MSG)
         ENDIF

         ! ------------------------------------------------------------- 
         ! Check if preferred datetime is within the range of available
         ! time slices. In this case, set tidx1 to the index of the 
         ! closest time slice that is not in the future. 
         ! ------------------------------------------------------------- 
         CALL Check_availYMDh ( nTime, availYMDh, prefYMDh, tidx1 )

         ! ------------------------------------------------------------- 
         ! If tidx1 couldn't be set in the call above, re-adjust 
         ! preferred year to the closest available year in the 
         ! time slices. Then repeat the check. 
         ! ------------------------------------------------------------- 
         IF ( tidx1 < 0 ) THEN
            CALL prefYMDh_adjustYear ( nTime, availYMDh, prefYMDh )

            ! verbose mode 
            IF ( verb ) THEN
               write(MSG,'(A30,I12)') 'adjusted preferred datetime: ', prefYMDh
               CALL HCO_MSG(MSG)
            ENDIF

            CALL Check_availYMDh ( nTime, availYMDh, prefYMDh, tidx1 )
         ENDIF

         ! ------------------------------------------------------------- 
         ! If tidx1 still isn't defined yet, i.e. prefYMDh is still 
         ! outside the range of availYMDh, set tidx1 to the closest
         ! available date. This must be 1 or nTime!
         ! ------------------------------------------------------------- 
         IF ( tidx1 < 0 ) THEN
            IF ( prefYMDh < availYMDh(1) ) THEN
               tidx1 = 1
            ELSE
               tidx1 = nTime
            ENDIF
         ENDIF
 
         ! verbose mode 
         IF ( verb ) THEN
            WRITE(MSG,'(A30,I12)') 'selected tidx1: ', tidx1
            CALL HCO_MSG(MSG)
         ENDIF

         ! ------------------------------------------------------------- 
         ! Now need to set upper time slice index tidx2. This index
         ! is only different from tidx1 if multiple hourly slices are
         ! read (--> prefHr = -1). In this case, check if there are 
         ! multiple time slices for the selected date (y/m/d). 
         ! ------------------------------------------------------------- 
         IF ( prefHr < 0 ) THEN
            CALL SET_TIDX2 ( nTime, availYMDH, tidx1, tidx2 ) 

            ! verbose mode 
            IF ( verb ) THEN
               WRITE(MSG,'(A30,I12)') 'selected tidx2: ', tidx1
               CALL HCO_MSG(MSG)
            ENDIF
         ELSE
            tidx2 = tidx1
         ENDIF 

      ! ================================================================
      ! Case 3: No time slice available. Set both indeces to zero. 
      ! ================================================================
      ELSE
         tidx1 = 0
         tidx2 = 0 
      ENDIF

      !-----------------------------------------------------------------
      ! If multiple time slices are read, extract time interval between
      ! time slices in memory (in hours). This is to make sure that the
      ! cycling between the slices will be done at the correct rate 
      ! (e.g. every hour, every 3 hours, ...).
      !-----------------------------------------------------------------
      IF ( tidx2 > tidx1 ) THEN
         Lct%Dct%Dta%DeltaT = availYMDh(tidx1+1) - availYMDh(tidx1)
      ELSE
         Lct%Dct%Dta%DeltaT = 0
      ENDIF

      ! verbose 
      IF ( verb ) THEN
         if ( nTime > 1 ) THEN
            write(MSG,'(A30,I12)') 'corresponding datetime 1: ', availYMDh(tidx1)
            CALL HCO_MSG(MSG)
            write(MSG,'(A30,I12)') 'corresponding datetime 2: ', availYMDh(tidx2)
            CALL HCO_MSG(MSG)
         endif
         write(MSG,'(A30,I12)') 'assigned delta t [h]: ', Lct%Dct%Dta%DeltaT 
         CALL HCO_MSG(MSG)
      ENDIF

      ! ----------------------------------------------------------------
      ! TODO: set time brackets 
      ! --> In future, we may want to set time brackets denoting the 
      ! previous and next time slice available in the netCDF file. This
      ! may become useful for temporal interpolations and more efficient
      ! data update calls (only update if new time slice is available). 
      ! ----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Prepare output, cleanup and leave 
      !-----------------------------------------------------------------

      ! ncYr and ncMt are the year and month fo the time slice to be
      ! used. These values may be required to convert units to 'per
      ! seconds'.
      IF ( tidx1 > 0 ) THEN
         ncYr = FLOOR( MOD(availYMDh(tidx1),10000000000) / 1.0d6 )
         ncMt = FLOOR( MOD(availYMDh(tidx1),1000000)     / 1.0d4 )
      ELSE
         ncYr = 0
         ncMt = 0
      ENDIF

      IF ( ASSOCIATED(availYMDh) ) DEALLOCATE(availYMDh)

      ! Return w/ success
      CALL HCO_LEAVE ( RC ) 

      END SUBROUTINE GET_TIMEIDX 
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_availYMDh  
!
! !DESCRIPTION: Checks if prefYMDh is within the range of availYMDh
! and returns the location of the closest vector element that is in
! the past as tidx1. tidx1 is set to -1 otherwise. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Check_availYMDh ( N, availYMDh, prefYMDh, tidx1 )
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER, INTENT(IN)   :: N
      INTEGER, INTENT(IN)   :: availYMDh(N)
      INTEGER, INTENT(IN)   :: prefYMDh
      INTEGER, INTENT(OUT)  :: tidx1
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      INTEGER               :: I

      !=================================================================
      ! Check_availYMDh begins here
      !=================================================================

      ! Init
      tidx1 = -1
 
      ! Return if preferred datetime not within the vector range
      IF ( prefYMDh < availYMDh(1) .OR. prefYMDh > availYMDh(N) ) RETURN

      ! get closest index that is not in the future
      DO I = 1, N
         IF ( availYMDh(I) == prefYMDh ) THEN
            tidx1 = I
            EXIT
         ENDIF

         IF ( availYMDh(I+1) > prefYMDh ) THEN
            tidx1 = I
            EXIT
         ENDIF
      ENDDO

      END SUBROUTINE Check_availYMDh 
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: PrefYMDh_AdjustYear 
!
! !DESCRIPTION: Adjusts the year in PrefYMDh to the closest available
! year in availYMDh 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE prefYMDh_adjustYear ( N, availYMDh, prefYMDh ) 
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER, INTENT(IN)     :: N
      INTEGER, INTENT(IN)     :: availYMDh(N)
      INTEGER, INTENT(INOUT)  :: prefYMDh

!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      INTEGER :: oldYear, prefMDh, newYear

      CONTINUE

      !=================================================================
      ! prefYMDH_adjustYear begins here! 
      !=================================================================

      ! Extract old year as well as oldMDh
      oldYear = FLOOR( MOD(prefYMDh,10000000000) / 1.0d6 )
      prefMDh = prefYMDh - (oldYear*1000000)

      ! Get new year: this is just the closest available year
      IF ( prefYMDh < availYMDh(1) ) THEN
         newYear = FLOOR( MOD(availYMDh(1),10000000000) / 1.0d6 )
      ELSE
         newYear = FLOOR( MOD(availYMDh(N),10000000000) / 1.0d6 )
      ENDIF

      ! Update prefYMDh
      prefYMDh = newYear*1000000 + prefMDh 

      END SUBROUTINE prefYMDh_adjustYear
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SET_TIDX2 
!
! !DESCRIPTION: sets the upper time slice index by selecting the range
! of all elements in availYMDh with the same date (year,month,day) as
! availYMDh(tidx1). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_TIDX2 ( N, availYMDh, tidx1, tidx2 ) 
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER, INTENT(IN)     :: N
      INTEGER, INTENT(IN)     :: availYMDh(N)
      INTEGER, INTENT(IN)     :: tidx1 
      INTEGER, INTENT(  OUT)  :: tidx2 
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      INTEGER :: YMD, I, IYMD

      !=================================================================
      ! SET_TIDX2 begins here! 
      !=================================================================

      ! Init
      tidx2 = tidx1

      ! Sanity check
      IF ( tidx1 == N ) RETURN

      ! Get wanted YMD
      YMD = floor(availYMDh(tidx1) / 1d2)

      ! See how many more tile slices with the same YMD exist from index
      ! tidx1 onwards.
      DO I = tidx1, N
         iYMD = floor(availYMDh(I) / 1d2)
         IF ( iYMD == YMD ) THEN
            tidx2 = I
         ELSEIF ( iYMD > YMD ) THEN
            EXIT
         ENDIF
      ENDDO

      END SUBROUTINE SET_TIDX2 
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NORMALIZE_AREA 
!
! !DESCRIPTION: Subroutine NORMALIZE\_AREA normalizes the given array
! by the surface area calculated from the given netCDF file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NORMALIZE_AREA ( HcoState, Array, ncFile, RC ) 
!
! !USES:
!
      USE NCDF_MOD,                 ONLY : NC_READ_GRID
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(HCO_State),  POINTER         :: HcoState
      REAL(sp),         POINTER         :: Array(:,:,:,:) ! Data
      CHARACTER(LEN=*), INTENT(IN   )   :: ncFile         ! netCDF file
      INTEGER,          INTENT(INOUT)   :: RC             ! Return code
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      REAL(hp)              :: DLAT, AREA
      INTEGER               :: J, NLON,  NLAT, NCRC
      REAL(sp), ALLOCATABLE :: XEDGE(:), YSIN(:)
      CHARACTER(LEN=255)    :: MSG, LOC

      !=================================================================
      ! NORNALIZE_AREA begins here! 
      !=================================================================

      ! Initialize
      LOC = 'NORMALIZE_AREA ( HCOI_DATAREAD_MOD.F90 )'

      ! Get array grid dimensions
      NLON = SIZE(ARRAY,1)
      NLAT = SIZE(ARRAY,2)

      ! Allocate grid edges
      ALLOCATE ( XEDGE(NLON+1) )
      ALLOCATE ( YSIN (NLAT+1) )
      XEDGE = 0_sp
      YSIN  = 0_sp

      ! Read input grid specifications
      CALL NC_READ_GRID( NLON, NLAT, TRIM(NCFILE), XEDGE, YSIN, RC)
      IF ( RC /= 0 ) THEN
         WRITE(MSG,*) 'Cannot read grid of file ' // TRIM(NCFILE) 
         CALL HCO_ERROR ( MSG, RC, THISLOC=LOC ) 
         RETURN
      ENDIF
 
      ! Loop over all latitudes
      DO J = 1, NLAT
         ! get grid box area in m2 for grid box with lower and upper latitude llat/ulat:
         ! Area = 2 * PI * Re^2 * DLAT / NLON, where DLAT = abs( sin(ulat) - sin(llat) ) 
         DLAT = ABS( YSIN(J+1) - YSIN(J) )
         AREA = ( 2_hp * HcoState%Phys%PI * DLAT * HcoState%Phys%Re**2 ) / REAL(NLON,hp)

         ! convert array data to m-2
         ARRAY(:,J,:,:) = ARRAY(:,J,:,:) / AREA 
      ENDDO

      ! Prompt a warning
      WRITE(MSG,*) 'No area unit found in ' // TRIM(NCFILE) // ' - convert to m-2!'
      CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )

      ! Deallocate arrays
      DEALLOCATE ( XEDGE, YSIN )

      END SUBROUTINE NORMALIZE_AREA 
!EOC
#endif
END MODULE HCOI_DATAREAD_MOD
