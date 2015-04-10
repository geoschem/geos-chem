!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_dataread_mod.F90 
!
! !DESCRIPTION: Module HCOIO\_DataRead\_Mod controls data processing 
! (file reading, unit conversion, regridding) for HEMCO.
!\\
!\\
! Currently, HEMCO can read data from the following data sources:
! \begin{itemize}
! \item Gridded data from netCDF file. The netCDF file should adhere to
!   the COARDS conventions and can hold data on resolutions different 
!   than then simulation grid. Regridding is performed as part of the 
!   data reading. Routine HCOIO\_DataRead is the driver routine to read
!   data from netCDF file. In an ESMF environment, this routine simply
!   calls down to the MAPL/ESMF - generic I/O routines. In a non-ESMF
!   environment, the HEMCO generic reading and remapping algorithms are
!   used. Those support vertical regridding, unit conversion, and more
!   (see below).
! \item Scalar data directly specified in the HEMCO configuration file.
!   If multiple values - separated by the separator sign (/) - are 
!   provided, they are interpreted as temporally changing values:
!   7 values = Sun, Mon, ..., Sat; 12 values = Jan, Feb, ..., Dec;
!   24 values = 0am, 1am, ..., 23pm (local time!). For masks, exactly
!   four values must be provided, interpreted as lower left and upper
!   right mask box corners (lon1/lat1/lon2/lat2).
! \item Country-specific data specified in a separate ASCII file. This
!   file must end with the suffix '.txt' and hold the country specific
!   values listed by country ID. The IDs must correspond to the IDs of
!   a corresponding (netCDF) mask file. The container name of this mask 
!   file must be given in the first line of the file, and must be listed 
!   HEMCO configuration file. ID 0 is reserved for the default values,
!   applied to all countries with no specific values listed. The .txt 
!   file must be structured as follows:
!
!   \# Country mask file name
!   CountryMask
!
!   \# CountryName CountryID CountryValues
!   DEFAULT   0 1.0/1.0/1.0/1.0/1.0/1.0/1/0
!   USA     840 0.8/0.9/1.0/1.1/1.2/1.1/0.9
!
!   The CountryValues are interpreted the same way as scalar values, 
!   except that they are applied to all grid boxes with the given 
!   country ID. 
! \end{itemize}
!\\
!\\
! Outside of an ESMF environment, the GEOS-Chem netCDF reading utilities
! are used to read netCDF data from disk. The selection of the time slice
! to be read depends on the current simulation time and the datetime 
! settings assigned to a given data container (set in the configuration
! file). These settings include:
! \begin{itemize}
! \item datetime range (srcTime), given as YYYY/MM/DD/hh. These can be 
!  given as fixed date (e.g. 2010/1/1/0), ranges (e.g. 1990-2010/1/1/0 
!  or 2000-2100/1-12/0/0-23), or using tokens (e.g. $YYYY/$MM/1/0).
!  Data is automatically updated if a 'dynamic' time attribute is given.
!  For instance, for attribute $YYYY/1/1/0 the file will be updated 
!  every year, attribute $YYYY/$MM/$DD/0 every day, etc. The date time
!  tokens are replaced with the current simulation datetime. If a range
!  is provided, only time stamps within the given range are being used.
!  If there is no exact match between the preferred datetime (determined
!  from srcTime) and the time slices in the netCDF file, the cycle flag
!  determines what time slice index is selected. 
! \item Cycling behavior. This determines what to do if there is no 
!  exact match between preferred datetime and available datetimes of a
!  file. The options are cycling (C, default), range (R), exact (E), 
!  and interpolation (I). If cycling is used, data is recycled if the
!  simulation time is outside of the available data time interval. 
!  If cycling is set to range, a data container is ignored if the
!  *simulation* time is outside the given range. For example, if the
!  range is set to 2010-2015/1-12/1/0, this data container is used 
!  for simulation dates between 2010 and 2015. If the actual netCDF
!  file data is outside that range, the closest available time slice
!  is selected using the algorithm described below.
!  If cycling is set to exact, HEMCO returns w/ an error if no time
!  slice can be found in the netCDF file that exactly matches the
!  preferred time slices.
!  Finally, if interpolation is selected, data will be interpolated
!  between two time slices if the current preferred datetime cannot
!  be found in the input data. 
!  If the preferred datetime does not match with any of the ncdf 
!  datetimes, the following method is used to select the time slice 
!  to be used:
!  If the preferred datetime is within the range of the available
!  dates, the closest available time stamp in the past is used in 
!  most cases. For example, assume a file contains January data 
!  between 2005 and 2010, and a simulation starts on July 2007. In 
!  this case, the data from Jan 2007 will be used and replaced with 
!  Jan 2008 as soon as the simulation date changes to 2008. 
!  If the datetimes of the netCDF file contain discontinuities (e.g.
!  don't have the same time interval between all time stamps), an 
!  attempt is made to maintain the highest cycling frequency. For 
!  instance, if a file contains monthly data for years 2005 and 2020 
!  and the srcTime attribute is set to $YYYY/1-12/1/0. For July 2008, 
!  this will use the data from July 2005, and not December 2005 (which 
!  would be the closest date). Data is updated to Aug 2005 once the 
!  simulation time changes to Aug 2008.
!  The same principles are applied if the current datetime is outside 
!  the data range of the ncdf file, but in this case it is also 
!  possible to use datetimes in the future.
!  If the interpolation option is enabled, the second time slice is
!  selected in the same manner as for slice 1, but for future dates.
!  If there is no additional time step in a file, e.g. the selected
!  time slice is the only/last one of a file, a check is made to see
!  if there is another file with the same filename but other datetime
!  tokens. If this is the case, interpolation will be performed between
!  these two files. Note that at the moment, the first time slice of
!  the second file is always used!!
!  For example, assume data for years 2005 and 2010 is stored in files 
!  file2005.nc and file2010.nc, respectively. The following entries 
!  in the HEMCO configuration file will enable annual interpolation
!  between these two files:
!  0 TEST file$YYYY.nc VAL 2005-2010/1/1/0 I ...
! \end{itemize} 
! !INTERFACE: 
!
MODULE HCOIO_DataRead_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_CharTools_Mod
  USE HCO_State_Mod,       ONLY : Hco_State
  USE HCO_DataCont_Mod,    ONLY : ListCont

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOIO_DataRead
  PUBLIC  :: HCOIO_ReadOther
!
! !REVISION HISTORY:
!  22 Aug 2013 - C. Keller   - Initial version
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
#if defined(ESMF_)
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_DataRead (ESMF/MAPL version)
!
! !DESCRIPTION: Interface routine between ESMF and HEMCO to obtain
! the data array for a given HEMCO data container. The data is obtained
! through the ExtData interface. The HEMCO source file attribute is taken
! to identify the ExtData pointer name.
!\\
!\\
! NOTE/TODO: For now, all arrays are copied into the HEMCO data array. 
! We may directly point to the ESMF arrays in future. 
!\\
!\\
! !INTERFACE:
  !
  SUBROUTINE HCOIO_DataRead( am_I_Root, HcoState, Lct, CloseFile, LUN, RC ) 
!
! !USES:
!
    USE ESMF
    USE MAPL_mod
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrInit

# include "MAPL_Generic.h"
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
    TYPE(HCO_State),  POINTER        :: HcoState
    TYPE(ListCont),   POINTER        :: Lct 
    LOGICAL,          INTENT(IN   )  :: CloseFile  ! Close file after reading?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: LUN
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  28 Aug 2013 - C. Keller   - Initial version
!  27 Aug 2014 - R. Yantosca - Err msg now displays hcoio_dataread_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER                    :: II, JJ, LL, TT
    INTEGER                    :: I, J, L, T
    INTEGER                    :: STAT
    REAL,             POINTER  :: Ptr3D(:,:,:)   => NULL() 
    REAL,             POINTER  :: Ptr2D(:,:)     => NULL() 
    TYPE(ESMF_State), POINTER  :: IMPORT         => NULL()
    CHARACTER(LEN=255)         :: MSG
    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCOIO_DATAREAD (hcoi_dataread_mod.F90)'
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam

    !=================================================================
    ! HCOIO_DATAREAD begins here
    !=================================================================

    ! For error handling
    Iam = LOC
    CALL HCO_ENTER ( LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Initialize output values
    LUN = -1

    ! Point to ESMF IMPORT object
    IMPORT => HcoState%IMPORT
    ASSERT_(ASSOCIATED(IMPORT))

    ! Verbose?
    IF ( HCO_IsVerb(2) ) THEN
       MSG = 'Reading from ExtData: ' // TRIM(Lct%Dct%Dta%ncFile)
       CALL HCO_MSG(MSG)
    ENDIF

    !-----------------------------------------------------------------
    ! Read 3D data from ESMF 
    !-----------------------------------------------------------------
    IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN

       ! Get data
       CALL MAPL_GetPointer( IMPORT, Ptr3D, &
                             TRIM(Lct%Dct%Dta%ncFile), RC=STAT )

       ! Check for MAPL error
       IF( STAT /= ESMF_SUCCESS ) THEN 
          MSG = 'Cannot get xyz pointer: ' // TRIM(Lct%Dct%Dta%ncFile)
          CALL HCO_ERROR ( MSG, RC ) 
          RETURN
       ENDIF

       ! Get array dimensions 
       II = SIZE(Ptr3D,1)
       JJ = SIZE(Ptr3D,2) 
       LL = SIZE(Ptr3D,3)
       TT = 1 

       ! Define HEMCO array if not yet defined.
       IF ( .NOT. ASSOCIATED(Lct%Dct%Dta%V3) ) THEN

          ! Use pointer if types match
          CALL FileData_ArrInit( Lct%Dct%Dta, TT, 0, 0, 0, RC )
          !CALL FileData_ArrInit( Lct%Dct%Dta, TT, II, JJ, LL, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Pointer to data. HEMCO expects data to have surface level at
       ! index 1 ('up').
       Lct%Dct%Dta%V3(1)%Val => Ptr3D(:,:,LL:1:-1)
       !Lct%Dct%Dta%V3(1)%Val(:,:,:) = Ptr3D(:,:,LL:1:-1)

    !-----------------------------------------------------------------
    ! Read 2D data from ESMF 
    !-----------------------------------------------------------------
    ELSEIF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN

       ! Get data
       CALL MAPL_GetPointer( IMPORT, Ptr2D, &
                             TRIM(Lct%Dct%Dta%ncFile), RC=STAT )

       ! Check for MAPL error 
       IF( STAT /= ESMF_SUCCESS ) THEN 
          MSG = 'Cannot get xy pointer: ' // TRIM(Lct%Dct%Dta%ncFile)
          CALL HCO_ERROR ( MSG, RC ) 
          RETURN
       ENDIF

       ! Get array dimensions 
       II = SIZE(Ptr2D,1)
       JJ = SIZE(Ptr2D,2) 
       LL = 1 
       TT = 1 

       ! Define HEMCO array pointer if not yet defined
       IF ( .NOT. ASSOCIATED(Lct%Dct%Dta%V2) ) THEN
          CALL FileData_ArrInit( Lct%Dct%Dta, TT, 0, 0, RC )
          !CALL FileData_ArrInit( Lct%Dct%Dta, TT, II, JJ, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Pointer to data
       Lct%Dct%Dta%V2(1)%Val => Ptr2D
       !Lct%Dct%Dta%V2(1)%Val = Ptr2D

    ENDIF
 
    !-----------------------------------------------------------------
    ! Cleanup and leave 
    !-----------------------------------------------------------------
    Ptr3D  => NULL()
    Ptr2D  => NULL()
    IMPORT => NULL()   

    ! Return w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOIO_DataRead
!EOC
#else
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_DataRead
!
! !DESCRIPTION: Reads a netCDF file and returns the regridded array in proper
! units. This routine uses the HEMCO generic data reading and regridding
! routines.
!\\
!\\
! Two different regridding algorithm are used: NCREGRID for 3D data with
! vertical regridding, and map\_a2a for all other data. map\_a2a also
! supports index-based remapping, while this feature is currently not
! possible in combination with NCREGRID.
!\\
!\\
! 3D data is vertically regridded onto the simulation grid on the sigma 
! interface levels. In order to calculate these levels correctly, the netCDF 
! vertical coordinate description must adhere to the CF - conventions. See 
! routine NC\_Get\_Sigma\_Levels in Ncdf\_Mod for more details.
!\\
!\\
! A simpler vertical interpolation scheme is used if (a) the number of 
! vertical levels of the input data corresponds to the number of levels
! on the simulation grid (direct mapping, no remapping), (b) the vertical
! level variable name (long\_name) contains the word "GEOS-Chem level". In
! the latter case, the vertical levels of the input data is interpreted as
! GEOS vertical levels and mapped onto the simulation grid using routine
! ModelLev\_Interpolate. 
!\\
!\\
! Argument CloseFile can be set to false to avoid closing the file. 
! Argument LUN can be used to read data from a previously opened stream. If
! the input value of LUN is greater than zero, the source file associated
! with the passed list container Lct is not being opened but the data is 
! read from stream LUN. The returned LUN value is equal to the LUN of the 
! file just being used if CloseFile is set to .FALSE., and to -1 otherwise. 
!\\
!\\  
! !INTERFACE:
!
  SUBROUTINE HCOIO_DataRead( am_I_Root, HcoState, Lct, CloseFile, LUN, RC ) 
!
! !USES:
!
    USE Ncdf_Mod,           ONLY : NC_Open
    USE Ncdf_Mod,           ONLY : NC_Close
    USE Ncdf_Mod,           ONLY : NC_Read_Var
    USE Ncdf_Mod,           ONLY : NC_Read_Arr
    USE Ncdf_Mod,           ONLY : NC_Get_Grid_Edges
    USE Ncdf_Mod,           ONLY : NC_Get_Sigma_Levels
    USE Ncdf_Mod,           ONLY : NC_ISMODELLEVEL
    USE HCO_Unit_Mod,       ONLY : HCO_Unit_Change
    USE HCO_Unit_Mod,       ONLY : HCO_Unit_ScalCheck
    USE HCO_Unit_Mod,       ONLY : HCO_IsUnitless
    USE HCO_Unit_Mod,       ONLY : HCO_IsIndexData
    USE HCO_Unit_Mod,       ONLY : HCO_UnitTolerance
    USE HCO_GeoTools_Mod,   ONLY : HCO_ValidateLon
    USE HCO_FileData_Mod,   ONLY : FileData_ArrCheck
    USE HCO_FileData_Mod,   ONLY : FileData_Cleanup
    USE HCOIO_MESSY_MOD,    ONLY : HCO_MESSY_REGRID
    USE HCO_INTERP_MOD,     ONLY : REGRID_MAPA2A 
    USE HCO_INTERP_MOD,     ONLY : ModelLev_Interpolate 
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_Get
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object
    TYPE(ListCont),   POINTER        :: Lct        ! HEMCO list container
    LOGICAL,          INTENT(IN   )  :: CloseFile  ! Close file after reading?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: LUN        ! LUN of file.
    INTEGER,          INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller   - Initial version
!  27 Aug 2014 - R. Yantosca - Err msg now displays hcoio_dataread_mod.F90
!  01 Oct 2014 - C. Keller   - Added file name parser
!  03 Oct 2014 - C. Keller   - Added vertical regridding capability
!  12 Dec 2014 - C. Keller   - Don't do vertical regridding if data is already
!                              on GEOS-Chem levels. 
!  31 Dec 2014 - C. Keller   - Now call ModelLev_Interpolate for model remapping
!                              of model levels.
!  15 Jan 2015 - C. Keller   - Now allow model level interpolation in 
!                              combination with MESSy (horizontal) regridding.
!  03 Feb 2015 - C. Keller   - Moved map_a2a regridding to hco_interp_mod.F90.
!  24 Mar 2015 - C. Keller   - Added arguments LUN and CloseFile.
!  27 Mar 2015 - R. Yantosca - Now use a FORMAT statement when printing the
!                              filename to the Unix stdout.
!  08 Apr 2015 - R. Yantosca - Bug fix: set KeepSpec=.TRUE. if there is no
!                              species in the container.  This prevents
!                              diffs in output in sp vs mp runs.
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)            :: thisUnit, LevUnit, LevName
    CHARACTER(LEN=255)            :: MSG 
    CHARACTER(LEN=1023)           :: srcFile, srcFile2
    INTEGER                       :: NX, NY
    INTEGER                       :: NCRC, Flag, AS
    INTEGER                       :: ncLun, ncLun2
    INTEGER                       :: nlon,   nlat,  nlev, nTime
    INTEGER                       :: lev1,   lev2,  dir 
    INTEGER                       :: tidx1,  tidx2,  ncYr,  ncMt
    INTEGER                       :: tidx1b, tidx2b, ncYr2, ncMt2
    INTEGER                       :: HcoID
    INTEGER                       :: nlatEdge, nlonEdge
    REAL(hp)                      :: MW_g, EmMW_g, MolecRatio
    REAL(sp)                      :: wgt1,   wgt2
    REAL(sp), POINTER             :: ncArr(:,:,:,:)   => NULL()
    REAL(sp), POINTER             :: ncArr2(:,:,:,:)  => NULL()
    REAL(hp), POINTER             :: SigEdge(:,:,:)   => NULL()
    REAL(hp), POINTER             :: SigLev (:,:,:)   => NULL()
    REAL(hp), POINTER             :: LonMid   (:)     => NULL()
    REAL(hp), POINTER             :: LatMid   (:)     => NULL()
    REAL(hp), POINTER             :: LevMid   (:)     => NULL()
    REAL(hp), POINTER             :: LonEdge  (:)     => NULL()
    REAL(hp), POINTER             :: LatEdge  (:)     => NULL()
    REAL(hp)                      :: UnitFactor
    LOGICAL                       :: KeepSpec
    LOGICAL                       :: FOUND
    LOGICAL                       :: IsModelLevel
    INTEGER                       :: UnitTolerance
    INTEGER                       :: AreaFlag, TimeFlag 
    INTEGER                       :: YMDha, YMDhb, YMDh1 
    INTEGER                       :: oYMDh1, oYMDh2

    ! Use MESSy regridding routines?
    LOGICAL                       :: UseMESSy

    !=================================================================
    ! HCOIO_DATAREAD begins here
    !=================================================================

    ! Enter
    CALL HCO_ENTER ('HCOIO_DATAREAD (hcoio_dataread_mod.F90)' , RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get unit tolerance set in configuration file
    UnitTolerance = HCO_UnitTolerance()
 
    ! For convenience, copy horizontal grid dimensions from HEMCO 
    ! state object
    NX = HcoState%NX
    NY = HcoState%NY

    ! ----------------------------------------------------------------
    ! Parse source file name. This will replace all tokens ($ROOT, 
    ! ($YYYY), etc., with valid values.
    ! ----------------------------------------------------------------
    CALL SrcFile_Parse ( am_I_Root, HcoState, Lct, srcFile, FOUND, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If file not found, return w/ error. No error if cycling attribute is 
    ! select to range. In that case, just make sure that array is empty.
    IF ( .NOT. FOUND ) THEN 
       IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE ) THEN
          CALL FileData_Cleanup( Lct%Dct%Dta, DeepClean=.FALSE. )
          MSG = 'No valid file found for current simulation time - data '// &
                'will be ignored - ' // TRIM(Lct%Dct%cName) 
          CALL HCO_WARNING ( MSG, RC, WARNLEV=1 )
          CALL HCO_LEAVE ( RC ) 
          RETURN
       ELSE
          MSG = 'File does not exist: ' // TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! Open netCDF
    ! ----------------------------------------------------------------
    IF ( LUN > 0 ) THEN
       ncLun = LUN

       ! Verbose mode
       IF ( HCO_IsVerb(2) ) THEN
          WRITE(MSG,*) '- Reading from existing stream: ', TRIM(srcFile)
          CALL HCO_MSG(MSG)
       ENDIF

    ELSE
       CALL NC_OPEN ( TRIM(srcFile), ncLun )

       ! Verbose mode
       IF ( HCO_IsVerb(1) ) THEN
          WRITE(MSG,*) '- Opening file: ', TRIM(srcFile)
          CALL HCO_MSG(MSG)
       ENDIF

       ! Also write to standard output
       WRITE( 6, 100 ) TRIM( srcFile )
 100   FORMAT( 'HEMCO: Opening ', a )

    ENDIF

    ! ----------------------------------------------------------------
    ! Extract time slice information
    ! This determines the lower and upper time slice index (tidx1 
    ! and tidx2) to be read based upon the time slice information 
    ! extracted from the file and the time stamp settings set in the
    ! HEMCO configuration file. Multiple time slices are only selected
    ! for weekdaily data or for 'autodetected' hourly data (using the
    ! wildcard character in the configuration file time attribute) or
    ! if data shall be interpolated between two (consecutive) time 
    ! slices. The weights to be assigned to those two time slices is
    ! also calculated in GET_TIMEIDX and returned as variables wgt1 
    ! and wgt2, respectively.
    ! ----------------------------------------------------------------
    CALL GET_TIMEIDX ( am_I_Root, HcoState, Lct,     &
                       ncLun,     tidx1,    tidx2,   &
                       wgt1,      wgt2,     oYMDh1,  &
                       YMDha,     YMDh1,    RC        )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Check for negative tidx1. tidx1 can still be negative if: 
    ! (a) CycleFlag is set to range and the current simulation 
    ! time is outside of the data time range. In this case, we 
    ! prompt a warning and make sure that there is no data 
    ! associated with this FileData container.
    ! (b) CycleFlag is set to exact and none of the data time 
    ! stamps matches the current simulation time exactly. Return 
    ! with error!
    !-----------------------------------------------------------------
    IF ( tidx1 < 0 ) THEN
       IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT ) THEN
          MSG = 'Exact time not found in ' // TRIM(srcFile) 
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ELSEIF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_CYCLE ) THEN
          MSG = 'Invalid time index: ' // TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ELSEIF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE ) THEN
          CALL FileData_Cleanup( Lct%Dct%Dta, DeepClean=.FALSE.)
          MSG = 'Simulation time is outside of time range provided for '//&
               TRIM(Lct%Dct%cName) // ' - data is ignored!'
          CALL HCO_WARNING ( MSG, RC, WARNLEV=1 )
          IF ( CloseFile ) THEN
             CALL NC_CLOSE ( ncLun )
             LUN = -1
          ELSE
             LUN = ncLUN
          ENDIF
          CALL HCO_LEAVE ( RC ) 
          RETURN
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! Read grid 
    ! ----------------------------------------------------------------

    ! Extract longitude midpoints
    CALL NC_READ_VAR ( ncLun, 'lon', nlon, thisUnit, LonMid, NCRC )
    IF ( NCRC /= 0 .OR. nlon == 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_LON', RC )
       RETURN 
    ENDIF
    IF ( INDEX( thisUnit, 'degrees_east' ) == 0 ) THEN
       MSG = 'illegal longitude unit in ' // &
            TRIM(srcFile)
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF
    ! Make sure longitude is steadily increasing.
    CALL HCO_ValidateLon( nlon, LonMid, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    
    ! Extract latitude midpoints
    CALL NC_READ_VAR ( ncLun, 'lat', nlat, thisUnit, LatMid, NCRC )
    IF ( NCRC /= 0 .OR. nlat == 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_LAT', RC )
       RETURN 
    ENDIF
    IF ( INDEX( thisUnit, 'degrees_north' ) == 0 ) THEN
       MSG = 'illegal latitude unit in ' // TRIM(srcFile)
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    ! Get level index if we are dealing with 3D data
    IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN

       ! Try to extract level midpoints
       LevName = 'lev'
       CALL NC_READ_VAR ( ncLun, LevName, nlev, LevUnit, LevMid, NCRC )
       IF ( NCRC /= 0 ) THEN
          CALL HCO_ERROR( 'NC_READ_LEV', RC )
          RETURN 
       ENDIF
       IF ( nlev == 0 ) THEN
          LevName = 'height'
          CALL NC_READ_VAR ( ncLun, LevName, nlev, LevUnit, LevMid, NCRC )
          IF ( NCRC /= 0 ) THEN
             CALL HCO_ERROR( 'NC_READ_LEV', RC )
             RETURN 
          ENDIF
       ENDIF

       ! Error check
       IF ( nlev == 0 ) THEN
          MSG = 'Source data of '//TRIM(Lct%Dct%cName)//' is not 3D: '//&
                TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Are these model levels? This will only return true if the long
       ! name of the level variable contains "GEOS-Chem level".
       IsModelLevel = NC_ISMODELLEVEL( ncLun, LevName )

       ! Set level indeces to be read
       ! NOTE: for now, always read all existing levels. Edit here to
       ! read only particular levels.
       lev1 = 1
       lev2 = nlev

       ! If # of levels are exactly # of simulation levels, assume that 
       ! they are on model levels. 
       ! This should probably be removed eventually, as it's better to 
       ! explicitly state model levels via the level long name 
       ! "GEOS-Chem level" (see above)!
       ! (ckeller, 12/12/14).
       IF ( nlev == HcoState%NZ ) IsModelLevel = .TRUE.

    ! For 2D data, set lev1 and lev2 to zero. This will ignore
    ! the level dimension in the netCDF reading call that follows.
    ELSE 
       nlev        = 0 
       lev1        = 0
       lev2        = 0
       IsModelLevel = .FALSE.
    ENDIF

    ! ----------------------------------------------------------------
    ! Read data
    ! ----------------------------------------------------------------

    ! Verbose mode
    IF ( HCO_IsVerb(2) ) THEN
       WRITE(MSG,*) 'Reading variable ', TRIM(Lct%Dct%Dta%ncPara)
       CALL HCO_MSG(MSG)
    ENDIF

    CALL NC_READ_ARR( fID     = ncLun,              &
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
                      varUnit = thisUnit,           &
                      wgt1    = wgt1,               &
                      wgt2    = wgt2,               &
                      MissVal = HCO_MISSVAL,        &
                      RC      = NCRC                 )

    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_ARRAY', RC )
       RETURN 
    ENDIF

    ! Check for missing values: set base emissions and masks to 0, and
    ! scale factors to 1. This will make sure that these entries will
    ! be ignored.
    CALL CheckMissVal ( Lct, ncArr )

    !-----------------------------------------------------------------
    ! Eventually do interpolation between files. This is a pretty 
    ! crude implementation for data interpolation between different 
    ! files. It is only applied to data that is marked as interpolated
    ! data and if no appropriate interpolation date could be found in
    ! the first file. This will only be the case if the preferred date-
    ! time is outside the file range.
    ! be found in the we check here if there exist 
    ! another file with the same data tokens, etc. but for a future 
    ! date. If this is the case 
    !-----------------------------------------------------------------
    IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_INTER .AND. wgt1 < 0.0_sp ) THEN

       ! Check if there exists another file for a future date 
       CALL SrcFile_Parse ( am_I_Root, HcoState, Lct, srcFile2, &
                            FOUND, RC, FUTURE=.TRUE. )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! If found, read data. Assume that all meta-data is the same.
       IF ( FOUND ) THEN

          ! Open file
          CALL NC_OPEN ( TRIM(srcFile2), ncLun2 )

          ! Define time stamp to be read. Use this call only
          ! to get the datetime of the first time slice (YMDh1).
          ! All other values will be ignored and reset below.
          CALL GET_TIMEIDX ( am_I_Root, HcoState, Lct,    &
                             ncLun2,    tidx1,    tidx2,  &
                             wgt1,      wgt2,     oYMDh2, & 
                             YMDhb,     YMDh1,    RC       )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Always read first time slice
          tidx1 = 1
          tidx2 = 1
          wgt1  = -1.0_sp
          wgt2  = -1.0_sp

          ! Read data and write into array ncArr2 
          CALL NC_READ_ARR( fID     = ncLun,              &
                            ncVar   = Lct%Dct%Dta%ncPara, &
                            lon1    = 1,                  &
                            lon2    = nlon,               &
                            lat1    = 1,                  &
                            lat2    = nlat,               &
                            lev1    = lev1,               &
                            lev2    = lev2,               &
                            time1   = tidx1,              &
                            time2   = tidx2,              &
                            ncArr   = ncArr2,             &
                            varUnit = thisUnit,           &
                            wgt1    = wgt1,               &
                            wgt2    = wgt2,               &
                            MissVal = HCO_MISSVAL,        &
                            RC      = NCRC                 )
          IF ( NCRC /= 0 ) THEN
             CALL HCO_ERROR( 'NC_READ_ARRAY (2)', RC )
             RETURN 
          ENDIF

          ! Eventually fissing values
          CALL CheckMissVal ( Lct, ncArr2 )

          ! Calculate weights to be applied to ncArr2 and ncArr1. These
          ! weights are calculated based on the originally preferred 
          ! datetime oYMDh1 and the selected datetime of file 1 (YMDha)
          ! and file 2 (YMDh1)
          CALL GetWeights ( YMDha, YMDh1, oYMDh1, wgt1, wgt2 ) 

          ! Apply weights
          ncArr = (wgt1 * ncArr) + (wgt2 * ncArr2)

          ! Verbose
          IF ( HCO_IsVerb(2) ) THEN
             MSG = 'Interpolated data between two files:'
             CALL HCO_MSG(MSG)
             MSG = '- File 1: ' // TRIM(srcFile)
             CALL HCO_MSG(MSG)
             WRITE(MSG,*) '   Time stamp used: ', YMDha
             CALL HCO_MSG(MSG)
             WRITE(MSG,*) '   Applied weight: ', wgt1
             CALL HCO_MSG(MSG)
             MSG = '- File 2: ' // TRIM(srcFile2)
             CALL HCO_MSG(MSG)
             WRITE(MSG,*) '   Time stamp used: ', YMDh1
             CALL HCO_MSG(MSG)
             WRITE(MSG,*) '   Applied weight: ', wgt2
             CALL HCO_MSG(MSG)
          ENDIF

          ! Cleanup
          IF ( ASSOCIATED(ncArr2) ) DEALLOCATE(ncArr2) 

          ! Close file
          CALL NC_CLOSE ( ncLun2 )
       ENDIF !FOUND
    ENDIF

    !-----------------------------------------------------------------
    ! Convert to HEMCO units 
    ! HEMCO data are all in kg/m2/s for fluxes and kg/m3 for 
    ! concentrations. Unit conversion is performed based on the
    ! unit on the input file and the srcUnit attribute given in the
    ! configuration file. By default, HEMCO will attempt to convert
    ! the units found in the input file to the standard quantities
    ! for mass (kg), area (m2 or m3), and time (s). For instance,
    ! g/cm2/hr will be converted to kg/m2/s. The exceptions to this
    ! rule are:
    ! 1. If srcUnit is set to '1', the input data are expected to
    !    be unitless. If the units string on the input file is none 
    !    of the units recognized by HEMCO as unitless, an error is
    !    returned if the unit tolerance setting is set to zero, or 
    !    a warning is prompted if unit tolerance is greater than zero. 
    ! 2. If srcUnit is set to 'count', no unit conversion is performed
    !    and data will be treated as 'index' data, e.g. regridding will
    !    preserve the absolute values.
    !
    ! Special attention needs to be paid to species that are emitted
    ! in quantities other than species molecules, e.g. molecules 
    ! carbon. For these species, the species MW differs from the 
    ! 'emitted MW', and the molecular ratio determines how many 
    ! molecules are being emitted per molecule species. By default, 
    ! HEMCO will attempt to convert all input data to kg emitted 
    ! species. If a species is emitted as kgC/m2/s and the input data 
    ! is in kg/m2/s, the mass will be adjusted based on the molecular 
    ! ratio and the ratio of emitted MW / species MW. Only input data
    ! that is already in kgC/m2/s will not be converted!
    ! This behavior can be avoided by explicitly setting the srcUnit
    ! to the same value as the input data unit. In this case, HEMCO 
    ! will not convert between species MW and emitted MW. 
    ! This is useful for cases where the input data does not
    ! contain data of the actual species, e.g. if VOC emissions are 
    ! calculated from scaled CO emissions. The scale factors then
    ! must include the conversion from CO to the VOC of interest! 
    !-----------------------------------------------------------------

    ! If OrigUnit is set to wildcard character: use unit from source file
    IF ( TRIM(Lct%Dct%Dta%OrigUnit) == HCO_WCD() ) THEN
       Lct%Dct%Dta%OrigUnit = TRIM(thisUnit)
    ENDIF

    ! If OrigUnit is set to '1' or to 'count', perform no unit 
    ! conversion.
    IF ( HCO_IsUnitLess(Lct%Dct%Dta%OrigUnit)  .OR. &
         HCO_IsIndexData(Lct%Dct%Dta%OrigUnit)       ) THEN

       ! Check if file unit is also unitless. This will return 0 for
       ! unitless, 1 for HEMCO emission unit, 2 for HEMCO conc. unit, 
       ! -1 otherwise.
       Flag = HCO_UNIT_SCALCHECK( thisUnit )
      
       ! Return with error if: (1) thisUnit is recognized as HEMCO unit and 
       ! unit tolerance is set to zero; (2) thisUnit is neither unitless nor
       ! a HEMCO unit and unit tolerance is set to zero or one.
       ! The unit tolerance is defined in the configuration file.
       IF ( Flag /= 0 .AND. UnitTolerance == 0 ) THEN
          MSG = 'Illegal unit: ' // TRIM(thisUnit) // '. File: ' // &
                TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Prompt a warning if thisUnit is not recognized as unitless.
       IF ( Flag /= 0 ) THEN 
          MSG = 'Data is treated as unitless, but file attribute suggests ' // &
                'it is not: ' // TRIM(thisUnit) // '. File: ' // TRIM(srcFile)
          CALL HCO_WARNING( MSG, RC, WARNLEV=1 )
       ENDIF

       ! Verbose mode
       IF ( HCO_IsVerb(2) ) THEN
          WRITE(MSG,*) 'Based on srcUnit attribute, no unit conversion is ', &
                       'performed: ', TRIM(Lct%Dct%Dta%OrigUnit)
          CALL HCO_MSG(MSG)
       ENDIF

    ! Convert to HEMCO units in all other cases. 
    ELSE

       ! For zero unit tolerance, make sure that thisUnit matches 
       ! with unit set in configuration file. For higher unit
       ! tolerances, prompt a level 3 warning. 
       IF ( TRIM(Lct%Dct%Dta%OrigUnit) /= TRIM(thisUnit) ) THEN
          MSG = 'File units do not match: ' // TRIM(thisUnit) // &
                ' vs. ' // TRIM(Lct%Dct%Dta%OrigUnit)    // &
                '. File: ' // TRIM(srcFile)

          IF ( UnitTolerance == 0 ) THEN
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ELSE
             CALL HCO_WARNING( MSG, RC, WARNLEV=3 )
          ENDIF
       ENDIF 

       ! Mirror species properties needed for unit conversion
       HcoID = Lct%Dct%HcoID
       IF ( HcoID > 0 ) THEN

          ! Emitted species molecular weight
          EmMW_g     = HcoState%Spc(HcoID)%EmMW_g
          MW_g       = HcoState%Spc(HcoID)%MW_g
          MolecRatio = HcoState%Spc(HcoID)%MolecRatio

          ! Species molecular weight and molecular ratio to be
          ! applied. Set to 1.0 if source unit matches input units.
          IF ( TRIM(Lct%Dct%Dta%OrigUnit) == TRIM(thisUnit) ) THEN
             !MW_g       = EmMW_g
             !MolecRatio = 1.0_hp
             KeepSpec    = .TRUE.
          ELSE
             !MW_g       = HcoState%Spc(HcoID)%MW_g
             !MolecRatio = HcoState%Spc(HcoID)%MolecRatio
             KeepSpec    = .FALSE.
          ENDIF

       ! If there is no species associated with this container, 
       ! it won't be possible to do unit conversion of mass. 
       ! This will cause an error if the input data is not in 
       ! units of kg already!
       ELSE
          KeepSpec   = .TRUE.
          MW_g       = -999.0_hp
          EmMW_g     = -999.0_hp
          MolecRatio = -999.0_hp
       ENDIF

       ! Now convert to HEMCO units. This attempts to convert mass, 
       ! area/volume and time to HEMCO standards (kg, m2/m3, s).
       ncYr  = FLOOR( MOD(oYMDh1,10000000000) / 1.0d6 )
       ncMt  = FLOOR( MOD(oYMDh1,1000000)     / 1.0d4 )
       IF ( ncYr == 0 ) CALL HcoClock_Get( cYYYY = ncYr, RC=RC ) 
       IF ( ncMt == 0 ) CALL HcoClock_Get( cMM   = ncMt, RC=RC ) 

       ! Verbose mode
       IF ( HCO_IsVerb(3) ) THEN
          WRITE(MSG,*) 'Unit conversion settings: ' 
          CALL HCO_MSG(MSG)
          WRITE(MSG,*) '- Species MW         : ', MW_g
          CALL HCO_MSG(MSG)
          WRITE(MSG,*) '- emitted compound MW: ', EmMW_g
          CALL HCO_MSG(MSG)
          WRITE(MSG,*) '- molecular ratio    : ', MolecRatio 
          CALL HCO_MSG(MSG)
          WRITE(MSG,*) '- keep input species : ', KeepSpec 
          CALL HCO_MSG(MSG)
          WRITE(MSG,*) '- Year, month        : ', ncYr, ncMt 
          CALL HCO_MSG(MSG)
       ENDIF

       CALL HCO_UNIT_CHANGE(                &
            Array         = ncArr,          &
            Units         = thisUnit,       &
            MW_IN         = MW_g,           & 
            MW_OUT        = EmMW_g,         & 
            MOLEC_RATIO   = MolecRatio,     & 
            KeepSpec      = KeepSpec,       & 
            YYYY          = ncYr,           &
            MM            = ncMt,           &
            AreaFlag      = AreaFlag,       &
            TimeFlag      = TimeFlag,       &
            FACT          = UnitFactor,     &
            RC            = RC               )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot convert units for ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG , RC )
          RETURN 
       ENDIF

       ! Verbose mode
       IF ( UnitFactor /= 1.0_hp ) THEN
          IF ( HCO_IsVerb(1) ) THEN
             WRITE(MSG,*) 'Data was in units of ', TRIM(thisUnit), &
                          ' - converted to HEMCO units by applying ', &
                          'scale factor ', UnitFactor
             CALL HCO_MSG(MSG)
          ENDIF
       ELSE
          IF ( HCO_IsVerb(2) ) THEN
             WRITE(MSG,*) 'Data was in units of ', TRIM(thisUnit), &
                          ' - unit conversion factor is ', UnitFactor 
             CALL HCO_MSG(MSG)
          ENDIF
       ENDIF

       ! Check for valid unit combinations, i.e. emissions must be kg/m2/s, 
       ! concentrations kg/m3. Eventually multiply by emission time step
       ! or divide by area to obtain those values.
      
       ! Concentration data
       IF ( AreaFlag == 3 .AND. TimeFlag == 0 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.
  
       ! If concentration data is per second (kg/m3/s), multiply by emission 
       ! time step to get concentration (kg/m3).
       ELSEIF ( AreaFlag == 3 .AND. TimeFlag == 1 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.

          ncArr = ncArr * HcoState%TS_EMIS
          MSG = 'Data converted from kg/m3/s to kg/m3: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(thisUnit)
          CALL HCO_WARNING( MSG, RC, WARNLEV=1 )
 
       ! Unitless data
       ELSEIF ( AreaFlag == -1 .AND. TimeFlag == -1 ) THEN
          ! nothing do to

       ! Emission data
       ELSEIF ( AreaFlag == 2 .AND. TimeFlag == 1 ) THEN
          ! nothing do to
 
       ! Emission data that is not per time (kg/m2): convert to kg/m2/s
       ELSEIF ( AreaFlag == 2 .AND. TimeFlag == 0 ) THEN
          ncArr = ncArr / HcoState%TS_EMIS
          MSG = 'Data converted from kg/m2 to kg/m2/s: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(thisUnit)
          CALL HCO_WARNING( MSG, RC, WARNLEV=1 )

       ! Emission data that is not per area (i.e. kg/s) needs to be converted
       ! to per area manually.
       ELSEIF ( AreaFlag == 0 .AND. TimeFlag == 1 ) THEN

          ! Get lat edges: those are read from file if possible, otherwise
          ! calculated from the lat midpoints.
          ! ==> Sine of lat is needed. Do conversion right here.
          CALL NC_GET_GRID_EDGES ( ncLun, 2, LatMid,   nlat, &
                                   LatEdge,  nlatEdge, NCRC   )
          IF ( NCRC /= 0 ) THEN
             MSG = 'Cannot read lat edge of ' // TRIM(srcFile)
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ENDIF 

          ! Now normalize data by area calculated from lat edges.
          CALL NORMALIZE_AREA( HcoState, ncArr,   nlon, &
                               LatEdge,  srcFile, RC     )
          IF ( RC /= HCO_SUCCESS ) RETURN

       ! All other combinations are invalid
       ELSE
          MSG = 'Unit must be unitless, emission or concentration: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(thisUnit)
          CALL HCO_ERROR ( MSG, RC )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Get horizontal grid edges 
    !-----------------------------------------------------------------

    ! Get longitude edges and make sure they are steadily increasing.
    CALL NC_GET_GRID_EDGES ( ncLun, 1, LonMid,   nlon, &
                             LonEdge,  nlonEdge, NCRC   )
    IF ( NCRC /= 0 ) THEN
       MSG = 'Cannot read lon edge of ' // TRIM(srcFile)
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF 
    CALL HCO_ValidateLon( nlonEdge, LonEdge, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get latitude edges (only if they have not been read yet
    ! for unit conversion)
    IF ( .NOT. ASSOCIATED( LatEdge ) ) THEN
       CALL NC_GET_GRID_EDGES ( ncLun, 2, LatMid,   nlat, &
                                LatEdge,  nlatEdge, NCRC   )
       IF ( NCRC /= 0 ) THEN
          MSG = 'Cannot read lat edge of ' // TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Determine regridding algorithm to be applied: use NCREGRID from
    ! MESSy only if we need to regrid vertical levels. For all other 
    ! fields, use the much faster map_a2a.
    ! Perform no vertical regridding if the vertical levels are model
    ! levels. Model levels are assumed to start at the surface, i.e.
    ! the first input level must correspond to the surface level. The 
    ! total number of vertical levels must not match the number of 
    ! vertical levels on the simulation grid. Data is not extrapolated
    ! beyond the existing levels.
    ! Vertical regridding based on NCREGRID will always map the input 
    ! data onto the entire simulation grid (no extrapolation beyond
    ! the vertical input coordinates).
    ! Index-based remapping can currently not be done with the MESSy
    ! routines, i.e. it is not possible to vertically regrid index- 
    ! based data.
    !-----------------------------------------------------------------

    UseMESSy = .FALSE.
    IF ( nlev > 1 .AND. .NOT. IsModelLevel ) THEN 
       UseMESSy = .TRUE.
    ENDIF
    IF ( HCO_IsIndexData(Lct%Dct%Dta%OrigUnit) .AND. UseMESSy ) THEN
       MSG = 'Cannot do MESSy regridding for index data: ' // &
             TRIM(srcFile)
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Use MESSy regridding
    !-----------------------------------------------------------------
    IF ( UseMESSy ) THEN
       IF ( HCO_IsVerb(2) ) THEN
          WRITE(MSG,*) '  ==> Use MESSy regridding (NCREGRID)'
          CALL HCO_MSG(MSG)
       ENDIF

       ! If we do MESSy regridding, we can only do one time step 
       ! at a time at the moment!
       IF ( tidx1 /= tidx2 ) THEN
          MSG = 'Cannot do MESSy regridding for more than one time step; ' &
                // TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       !--------------------------------------------------------------
       ! Eventually get sigma levels
       ! Vertical regridding is performed on sigma interface levels:
       ! sigma(i,j,l) = p(i,j,l) / ps(i,j)
       ! NC_Get_Sigma_Levels attempts to create the sigma levels from
       ! the content of the netCDF file. 
       ! For now, it is assumed that all input data is on vertical 
       ! mid-point levels, and the interface values are calculated 
       ! by linear interpolation of the mid-point values in a second
       ! step.
       ! For model levels, the sigma levels don't need to be known
       ! as vertical interpolation will be done based on subroutine
       ! ModelLev_Interpolate (within HCO_MESSY_REGRID).
       !--------------------------------------------------------------
       IF ( nlev > 1 .AND. .NOT. IsModelLevel ) THEN

          ! Get sigma levels
          CALL NC_Get_Sigma_Levels ( fID     = ncLun,   &
                                     ncFile  = srcFile, &
                                     levName = LevName, &
                                     lon1    = 1,       &
                                     lon2    = nlon,    &
                                     lat1    = 1,       &
                                     lat2    = nlat,    &
                                     lev1    = 1,       &
                                     lev2    = nlev,    &
                                     time    = tidx1,   &
                                     SigLev  = SigLev,  &
                                     Dir     = dir,     &
                                     RC      = NCRC      )
          IF ( NCRC /= 0 ) THEN
             CALL HCO_ERROR( 'Cannot read sigma levels of '//TRIM(srcFile), RC )
             RETURN
          ENDIF

          ! Interpolate onto edges
          CALL SigmaMidToEdges ( am_I_Root, SigLev, SigEdge, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
         
          ! Sigma levels are not needed anymore
          IF ( ASSOCIATED(SigLev) ) DEALLOCATE(SigLev)

          !-----------------------------------------------------------
          ! Flip vertical axis if positive axis is 'down', i.e. level
          ! index 1 is the top of the atmosphere
          !-----------------------------------------------------------
          IF ( dir == -1 ) THEN
             SigEdge(:,:,:  ) = SigEdge(:,:,nlev+1:1:-1  )
             NcArr  (:,:,:,:) = NcArr  (:,:,nlev  :1:-1,:)
          ENDIF

       ENDIF ! nlev>1

       ! Now do the regridding
       CALL HCO_MESSY_REGRID ( am_I_Root, HcoState,     NcArr,   &
                               LonEdge,   LatEdge,      SigEdge, &
                               Lct,       IsModelLevel, RC        )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Cleanup
       IF ( ASSOCIATED(SigEdge) ) DEALLOCATE(SigEdge)

    !-----------------------------------------------------------------
    ! Use map_a2a regridding
    !-----------------------------------------------------------------
    ELSE
       IF ( HCO_IsVerb(2) ) THEN
          WRITE(MSG,*) '  ==> Use map_a2a regridding'
          CALL HCO_MSG(MSG)
       ENDIF

       CALL REGRID_MAPA2A ( am_I_Root, HcoState, NcArr, &
                            LonEdge,   LatEdge,  Lct,   RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 

    ENDIF

    ! ----------------------------------------------------------------
    ! Close netCDF
    ! ----------------------------------------------------------------
    IF ( CloseFile ) THEN
       CALL NC_CLOSE ( ncLun )
       LUN = -1
    ELSE
       LUN = ncLun
    ENDIF      

    !-----------------------------------------------------------------
    ! Cleanup and leave 
    !-----------------------------------------------------------------
    IF ( ASSOCIATED ( ncArr   ) ) DEALLOCATE ( ncArr   )
    IF ( ASSOCIATED ( LonMid  ) ) DEALLOCATE ( LonMid  )
    IF ( ASSOCIATED ( LatMid  ) ) DEALLOCATE ( LatMid  )
    IF ( ASSOCIATED ( LevMid  ) ) DEALLOCATE ( LevMid  )
    IF ( ASSOCIATED ( LonEdge ) ) DEALLOCATE ( LonEdge )
    IF ( ASSOCIATED ( LatEdge ) ) DEALLOCATE ( LatEdge )

    ! Return w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE HCOIO_DataRead
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_TimeIdx
!
! !DESCRIPTION: Returns the lower and upper time slice index (tidx1
! and tidx2, respectively) to be read. These values are determined 
! based upon the time slice information extracted from the netCDF file, 
! the time stamp settings set in the config. file, and the current 
! simulation date.
!\\
!\\
! Return arguments wgt1 and wgt2 denote the weights to be given to
! the two time slices. This is only of relevance for data that shall
! be interpolated between two (not necessarily consecutive) time slices. 
! In all other cases, the returned weights are negative and will be 
! ignored.
!\\
!\\
! Also returns the time slice year and month, as these values may be
! used for unit conversion. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_TIMEIDX( am_I_Root, HcoState, Lct,     &
                          ncLun,     tidx1,    tidx2,   &
                          wgt1,      wgt2,     oYMDh,   &
                          YMDh,      YMDh1,    RC        )
!
! !USES:
!
    USE Ncdf_Mod,      ONLY : NC_Read_Time_YYYYMMDDhh
    USE HCO_tIdx_Mod,  ONLY : HCO_GetPrefTimeAttr
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root ! Root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState  ! HcoState object
    TYPE(ListCont),   POINTER        :: Lct       ! List container
    INTEGER,          INTENT(IN   )  :: ncLun     ! open ncLun
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)  :: tidx1     ! lower time idx
    INTEGER,          INTENT(  OUT)  :: tidx2     ! upper time idx
    REAL(sp),         INTENT(  OUT)  :: wgt1      ! weight to tidx1
    REAL(sp),         INTENT(  OUT)  :: wgt2      ! weight to tidx2
    INTEGER,          INTENT(  OUT)  :: oYMDh     ! preferred time slice 
    INTEGER,          INTENT(  OUT)  :: YMDh      ! selected time slice 
    INTEGER,          INTENT(  OUT)  :: YMDh1     ! 1st time slice in file 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!  27 Feb 2015 - C. Keller - Added weigths
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
    CHARACTER(LEN=255)    :: MSG
    CHARACTER(LEN=1023)   :: MSG_LONG
    INTEGER               :: tidx1a
    INTEGER               :: nTime,  T, CNT, NCRC 
    INTEGER               :: prefYr, prefMt, prefDy, prefHr
    INTEGER               :: refYear
    INTEGER               :: origYMDh, prefYMDh
    INTEGER, POINTER      :: availYMDh(:) => NULL() 
    LOGICAL               :: ExitSearch 
    LOGICAL               :: verb

    !=================================================================
    ! GET_TIMEIDX begins here
    !=================================================================

    ! Init 
    CALL HCO_ENTER ('GET_TIMEIDX (hco_dataread_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    verb = HCO_IsVerb(3)

    ! Initialize
    wgt1  = -1.0_sp
    wgt2  = -1.0_sp
    oYMDh = 0
    YMDh  = 0
    YMDh1 = 0
 
    ! ---------------------------------------------------------------- 
    ! Extract netCDF time slices (YYYYMMDDhh) 
    ! ----------------------------------------------------------------
    CALL NC_READ_TIME_YYYYMMDDhh ( ncLun, nTime,    availYMDH, &
                                   refYear=refYear, RC=NCRC     )     
    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_TIME_YYYYMMDDhh', RC )
       RETURN 
    ENDIF

    ! Return warning if netCDF reference year prior to 1901: it seems 
    ! like there are some problems with that and the time slices can be 
    ! off by one day!
    IF ( (refYear <= 1900) .AND. (nTime > 0) ) THEN
       MSG = 'ncdf reference year is prior to 1901 - ' // &
            'time stamps may be wrong!'
       CALL HCO_WARNING ( MSG, RC, WARNLEV=1 )
    ENDIF

    ! verbose mode 
    IF ( verb ) THEN
       write(MSG,'(A30,I12)') '# time slices found: ', nTime
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

    ! ---------------------------------------------------------------- 
    ! Get preferred time stamp to read based upon the specs set in the
    ! config. file. 
    ! This can return value -1 for prefHr, indicating that all  
    ! corresponding time slices shall be read.
    ! This call will return -1 for all date attributes if the 
    ! simulation date is outside of the data range given in the 
    ! configuration file.
    ! ---------------------------------------------------------------- 
    CALL HCO_GetPrefTimeAttr ( Lct, prefYr, prefMt, prefDy, prefHr, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Check if we are outside of provided range
    IF ( prefYr < 0 .OR. prefMt < 0 .OR. prefDy < 0 ) THEN
     
       ! This should only happen for 'range' data 
       IF ( Lct%Dct%Dta%CycleFlag /= HCO_CFLAG_RANGE ) THEN
          MSG = 'Cannot get preferred datetime for ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! If this part of the code gets executed, the data associated 
       ! with this container shall not be used at the current date.
       ! To do so, set the time indeces to -1 and leave right here.
       tidx1 = -1
       tidx2 = -1
 
       ! Leave w/ success
       CALL HCO_LEAVE( RC )
       RETURN 
    ENDIF

    ! origYMDh is the preferred datetime. Store into shadow variable
    ! prefYMDh. prefYMDh may be adjusted if origYMDh is outside of the
    ! netCDF datetime range.
    origYMDh = prefYr*1000000 + prefMt*10000 + &
               prefDy*100 + max(prefHr,0)
    prefYMDh = origYMDh

    ! verbose mode
    IF ( verb ) THEN
       write(MSG,'(A30,I12)') 'preferred datetime: ', prefYMDh
       CALL HCO_MSG(MSG)
    ENDIF

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
       tidx1   = -1
       tidx2   = -1 

       ! ------------------------------------------------------------- 
       ! Check if preferred datetime prefYMDh is within the range
       ! available time slices, e.g. it falls within the interval
       ! of availYMDh. In this case, set tidx1 to the index of the 
       ! closest time slice that is not in the future. 
       ! ------------------------------------------------------------- 
       CALL Check_AvailYMDh ( Lct, nTime, availYMDh, prefYMDh, tidx1a )

       ! ------------------------------------------------------------- 
       ! Check if we need to continue search. Even if the call above
       ! returned a time slice, it may be possible to continue looking
       ! for a better suited time stamp. This is only the case if
       ! there are discontinuities in the time stamps, e.g. if a file
       ! contains monthly data for 2005 and 2020. In that case, the
       ! call above would return the index for Dec 2005 for any 
       ! simulation date between 2005 and 2010 (e.g. July 2010),
       ! whereas it makes more sense to use July 2005 (and eventually
       ! interpolate between the July 2005 and July 2020 data).
       ! The IsClosest command checks if there are any netCDF time
       ! stamps (prior to the selected one) that are closer to each
       ! other than the difference between the preferred time stamp
       ! prefYMDh and the currently selected time stamp 
       ! availYMDh(tidx1a). In that case, it continues the search by
       ! updating prefYMDh so that it falls within the range of the
       ! 'high-frequency' interval.
       ! ------------------------------------------------------------- 
       ExitSearch = .FALSE.
       IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT ) THEN
          ExitSearch = .TRUE.
       ELSE IF ( tidx1a > 0 ) THEN 
          ExitSearch = IsClosest( prefYMDh, availYMDh, nTime, tidx1a )
       ENDIF 

       ! Write to tidx1 if this is the best match. 
       IF ( ExitSearch ) THEN
          tidx1 = tidx1a

       ! ------------------------------------------------------------- 
       ! If search shall be continued, adjust preferred year, then 
       ! month, then day to the closest available year (month, day) 
       ! in the time slices, and check if this is a better match.
       ! ------------------------------------------------------------- 
       ELSE
         
          ! Adjust year, month, and day (in this order).
          CNT  = 0
          DO 
             CNT = CNT + 1
             IF ( ExitSearch .OR. CNT > 3 ) EXIT

             ! Adjust prefYMDh at the given level (1=Y, 2=M, 3=D)
             CALL prefYMDh_Adjust ( nTime, availYMDh, prefYMDh, CNT, tidx1a )

             ! verbose mode 
             IF ( verb ) THEN
                write(MSG,'(A30,I12)') 'adjusted preferred datetime: ', prefYMDh
                CALL HCO_MSG(MSG)
             ENDIF
      
             ! check for time stamp with updated date/time
             CALL Check_AvailYMDh ( Lct, nTime, availYMDh, prefYMDh, tidx1a )
 
             ! Can we leave now?
             ExitSearch = IsClosest( prefYMDh, availYMDh, nTime, tidx1a )
             IF ( ExitSearch ) tidx1 = tidx1a 
 
          ENDDO
       ENDIF   

       ! ------------------------------------------------------------- 
       ! If tidx1 still isn't defined, i.e. prefYMDh is still 
       ! outside the range of availYMDh, set tidx1 to the closest
       ! available date. This must be 1 or nTime! 
       ! ------------------------------------------------------------- 
       IF ( .NOT. ExitSearch ) THEN 
          IF ( prefYMDh < availYMDh(1) ) THEN
             tidx1 = 1
          ELSE
             tidx1 = nTime
          ENDIF
       ENDIF
 
       ! ------------------------------------------------------------- 
       ! If we are dealing with weekday data, pick the slice to be
       ! used based on the current day of week. 
       ! The ncDys flag has been set in subroutine HCO_ExtractTime
       ! (hco_tidx_mod.F90) based upon the time attributes set in the
       ! configuration file. It can have the following values:
       ! >0  : specific days are given.
       ! -1  : wildcard (autodetect)
       ! -10 : WD (weekday). 
       ! -999: determine from current simulation day.
       ! For specific days or if determined from the current datetime
       ! (flags >0 or -999), the weekday is not taken into account.
       ! If auto-detection is enabled, days are treated as weekday if
       ! (and only if) there are exactly 7 time slices. Otherwise, they
       ! are interpreted as 'regular' day data. 
       ! If flag is set to -10, e.g. time attribute is 'WD', the current 
       ! time index is assumed to hold Sunday data, with the following 
       ! six slices being Mon, Tue, ..., Sat. For weekdaily data, all 
       ! seven time slices will be read into memory so that at any given
       ! time, the local weekday can be taken (weekdaily data is always 
       ! assumed to be in local time).
       ! ------------------------------------------------------------- 

       ! Day flag is -1: wildcard
       IF ( Lct%Dct%Dta%ncDys(1) == -1 .AND. nTime == 7 ) THEN
          tidx1 = 1
          tidx2 = nTime 

          ! Make sure data is treated in local time
          Lct%Dct%Dta%IsLocTime = .TRUE.

       ! Day flag is -10: WD
       ELSEIF ( Lct%Dct%Dta%ncDys(1) == -10 ) THEN

          ! There must be at least 7 time slices 
          IF ( nTime < 7 ) THEN
             MSG = 'Data must have exactly 7 time slices '// &
                   'if you set day attribute to WD: '//TRIM(Lct%Dct%cName)
             CALL HCO_ERROR ( MSG, RC )
             RETURN
          ENDIF
 
          ! If there are exactly seven time slices, interpret them as
          ! the seven weekdays. 
          IF ( nTime == 7 ) THEN
             tidx1 = 1
             tidx2 = 7
             
          ! If there are more than 7 time slices, interpret the current
          ! selected index as sunday of the current time frame (e.g. sunday
          ! data of current month), and select the time slice index
          ! accordingly. This requires that there are at least 6 more time
          ! slices following the current one. 
          ELSE
             IF ( tidx1 < 0 ) THEN
                WRITE(MSG,*) 'Cannot get weekday slices for: ', &
                   TRIM(Lct%Dct%cName), '. Cannot find first time slice.' 
                CALL HCO_ERROR ( MSG, RC )
                RETURN
             ENDIF

             IF ( (tidx1+6) > nTime ) THEN
                WRITE(MSG,*) 'Cannot get weekday for: ',TRIM(Lct%Dct%cName), &
                   '. There are less than 6 additional time slices after ',  &
                   'selected start date ', availYMDh(tidx1)
                CALL HCO_ERROR ( MSG, RC )
                RETURN
             ENDIF
             tidx2 = tidx1 + 6
          ENDIF

          ! Make sure data is treated in local time
          Lct%Dct%Dta%IsLocTime = .TRUE.

       ENDIF

       ! ------------------------------------------------------------- 
       ! Now need to set upper time slice index tidx2. This index
       ! is only different from tidx1 if:
       ! (1) We interpolate between two time slices, i.e. TimeCycle
       !     attribute is set to 'I'. In this case, we simply pick 
       !     the next higher time slice index and calculate the 
       !     weights for time1 and time2 based on the current time.
       ! (2) Multiple hourly slices are read (--> prefHr = -1 or -10, 
       !     e.g. hour attribute in config. file was set to wildcard 
       !     character or data is in local hours). In this case, 
       !     check if there are multiple time slices for the selected 
       !     date (y/m/d).
       ! tidx2 has already been set to proper value above if it's
       ! weekday data.
       ! -------------------------------------------------------------
       IF ( tidx2 < 0 ) THEN

          ! Interpolate between dates
          IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_INTER ) THEN
        
             CALL GetIndex2Interp( am_I_Root, Lct,      nTime,    &
                                   availYMDh, prefYMDh, origYMDh, &
                                   tidx1,     tidx2,    wgt1,     &
                                   wgt2,      RC                   )
             IF ( RC /= HCO_SUCCESS ) RETURN

          ! Check for multiple hourly data
          ELSEIF ( tidx1 > 0 .AND. prefHr < 0 ) THEN
             CALL SET_TIDX2 ( nTime, availYMDH, tidx1, tidx2 )    

             ! Denote as local time if necessary
             IF ( Lct%Dct%Dta%ncHrs(1) == -10 ) THEN
                Lct%Dct%Dta%IsLocTime = .TRUE.
             ENDIF
          ELSE
             tidx2 = tidx1
          ENDIF
       ENDIF   

    ! ================================================================
    ! Case 3: No time slice available. Set both indeces to zero. Data
    ! with no time stamp must have CycleFlag 'Cycling'.
    ! ================================================================
    ELSE
       IF ( Lct%Dct%Dta%CycleFlag /= HCO_CFLAG_CYCLE ) THEN
          MSG = 'Field has no time/date variable - cycle flag must' // &
                'be set to `C` in the HEMCO configuration file:'    // &
                TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC )
          RETURN
       ENDIF

       tidx1 = 0
       tidx2 = 0 
    ENDIF

    !-----------------------------------------------------------------
    ! Sanity check: if CycleFlag is set to 'Exact', the file time stamp
    ! must exactly match the current time.
    !-----------------------------------------------------------------
    IF ( (Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT) .AND. (tidx1 > 0) ) THEN
       IF ( availYMDh(tidx1) /= prefYMDh ) THEN
          tidx1 = -1
          tidx2 = -1
       ENDIF
    ENDIF
    
    !-----------------------------------------------------------------
    ! If multiple time slices are read, extract time interval between
    ! time slices in memory (in hours). This is to make sure that the
    ! cycling between the slices will be done at the correct rate 
    ! (e.g. every hour, every 3 hours, ...).
    !-----------------------------------------------------------------
    IF ( (tidx2>tidx1) .AND. (Lct%Dct%Dta%CycleFlag/=HCO_CFLAG_INTER) ) THEN
       Lct%Dct%Dta%DeltaT = YMDh2hrs( availYMDh(tidx1+1) - availYMDh(tidx1) )
    ELSE
       Lct%Dct%Dta%DeltaT = 0
    ENDIF

    ! verbose mode 
    IF ( verb ) THEN
       WRITE(MSG,'(A30,I12)') 'selected tidx1: ', tidx1
       CALL HCO_MSG(MSG)
       IF ( tidx1 > 0 ) THEN
          WRITE(MSG,'(A30,I12)') 'corresponding datetime 1: ', availYMDh(tidx1)
          CALL HCO_MSG(MSG)
          IF ( wgt1 >= 0.0_sp ) THEN
             WRITE(MSG,*) 'weight1: ', wgt1
             CALL HCO_MSG(MSG)
          ENDIF
       ENDIF

       IF ( (tidx2 /= tidx1) ) THEN
          WRITE(MSG,'(A30,I12)') 'selected tidx2: ', tidx2
          CALL HCO_MSG(MSG)
          WRITE(MSG,'(A30,I12)') 'corresponding datetime 2: ', availYMDh(tidx2)
          CALL HCO_MSG(MSG)
          IF ( wgt1 >= 0.0_sp ) THEN
             WRITE(MSG,*) 'weight2: ', wgt2
             CALL HCO_MSG(MSG)
          ENDIF
       ENDIF

       WRITE(MSG,'(A30,I12)') 'assigned delta t [h]: ', Lct%Dct%Dta%DeltaT 
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) 'local time? ', Lct%Dct%Dta%IsLocTime
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
       YMDh  = availYMDh(tidx1)
       YMDh1 = availYMDh(1)
       oYMDh = origYMDh
    ENDIF

    IF ( ASSOCIATED(availYMDh) ) DEALLOCATE(availYMDh)

    ! Return w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE GET_TIMEIDX
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_AvailYMDh  
!
! !DESCRIPTION: Checks if prefYMDh is within the range of availYMDh
! and returns the location of the closest vector element that is in
! the past (--> tidx1). tidx1 is set to -1 otherwise. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_AvailYMDh( Lct, N, availYMDh, prefYMDh, tidx1 )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER      :: Lct 
    INTEGER,          INTENT(IN)   :: N
    INTEGER,          INTENT(IN)   :: availYMDh(N)
    INTEGER,          INTENT(IN)   :: prefYMDh
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)  :: tidx1
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER :: I

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

       ! Check if next time slice is in the future, in which case the
       ! current slice is selected. Don't do this for a CycleFlag of
       ! 3 (==> exact match).
!       IF ( availYMDh(I+1) > prefYMDh ) THEN 
       IF ( (availYMDh(I+1)        >  prefYMDh       ) .AND. &
            (Lct%Dct%Dta%CycleFlag /= HCO_CFLAG_EXACT) ) THEN
          tidx1 = I
          EXIT
       ENDIF
    ENDDO

  END SUBROUTINE Check_AvailYMDh
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: prefYMDh_Adjust
!
! !DESCRIPTION: Adjusts prefYMDh to the closest available time attribute. Can
! be adjusted for year (level=1), month (level=2), or day (level=3).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE prefYMDh_Adjust( N, availYMDh, prefYMDh, level, tidx1 ) 
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)     :: N
    INTEGER, INTENT(IN)     :: availYMDh(N)
    INTEGER, INTENT(IN)     :: level
    INTEGER, INTENT(IN)     :: tidx1
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)  :: prefYMDh
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!  17 Jul 2014 - C. Keller - Now allow to adjust year, month, or day. 
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    
    INTEGER          :: I, IMIN, IMAX
    INTEGER          :: origYr,  origMt,  origDy, origHr
    INTEGER          :: refAttr, tmpAttr, newAttr
    INTEGER          :: iDiff,   minDiff
    INTEGER(8)       :: modVal
    REAL(dp)         :: div

    !=================================================================
    ! prefYMDh_Adjust begins here! 
    !=================================================================

    ! Get original Yr, Mt, Dy and Hr
    origYr = FLOOR( MOD(prefYMDh, 10000000000) / 1.0d6 )
    origMt = FLOOR( MOD(prefYMDh, 1000000    ) / 1.0d4 )
    origDy = FLOOR( MOD(prefYMDh, 10000      ) / 1.0d2 )
    origHr = FLOOR( MOD(prefYMDh, 100        ) / 1.0d0 )

    ! Extract new attribute from availYMDh and insert into prefYMDh. Pick
    ! closest available value.
    SELECT CASE ( level ) 
       ! --- Year
       CASE ( 1 )
          modVal  = 10000000000
          div     = 1.0d6
          refAttr = origYr

       ! --- Month
       CASE ( 2 )
          modVal  = 1000000
          div     = 1.0d4
          refAttr = origMt

       ! --- Day 
       CASE ( 3 )
          modVal  = 10000
          div     = 1.0d2
          refAttr = origMt

       CASE DEFAULT
          RETURN
    END SELECT

    ! Maximum loop number:
    ! If tidx1 is already set, only search values in the past.
    IF ( tidx1 > 0 ) THEN
       IMIN = 1
       IMAX = tidx1

    ! If tidx1 is not yet set, prefYMDh must be outside the range of availYMDh.
    ! Pick only the closest available time stamp.
    ELSE
       IF ( prefYMDh > availYMDh(1) ) THEN
          IMIN = N
          IMAX = N
       ELSE
          IMIN = 1
          IMAX = 1
       ENDIF
    ENDIF

    ! Select current minimum value
    minDiff = 10000000000000000
    newAttr = -1
    DO I = IMIN, IMAX 
       tmpAttr = FLOOR( MOD(availYMDh(I),modVal) / div )
       iDiff   = ABS( tmpAttr - refAttr )
       IF ( iDiff < minDiff ) THEN
          newAttr = tmpAttr
          minDiff = iDiff
       ENDIF
    ENDDO

    ! Just reuse current value if no better value could be found
    IF ( newAttr < 0 ) THEN
       newAttr = refAttr
    ENDIF

    ! Update variable
    ! --- Year
    IF ( level == 1 ) THEN
       prefYMDh = newAttr * 1000000 + origMt * 10000 + origDy * 100 + origHr

    ! --- Month 
    ELSEIF ( level == 2 ) THEN
       prefYMDh = origYr * 1000000 + newAttr * 10000 + origDy * 100 + origHr

    ! --- Day
    ELSEIF ( level == 3 ) THEN
       prefYMDh = origYr * 1000000 + origMt * 10000 + newAttr * 100 + origHr
    ENDIF

  END SUBROUTINE prefYMDh_Adjust
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_tIdx2 
!
! !DESCRIPTION: sets the upper time slice index by selecting the range
! of all elements in availYMDh with the same date (year,month,day) as
! availYMDh(tidx1). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_tIdx2( N, availYMDh, tidx1, tidx2 ) 
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: N
    INTEGER, INTENT(IN)  :: availYMDh(N)
    INTEGER, INTENT(IN)  :: tidx1 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: tidx2 
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
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

  END SUBROUTINE Set_tIdx2
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsClosest 
!
! !DESCRIPTION: function IsClosest returns true if the selected time index
! is the 'closest' one. It is defined as being closest if: 
! (a) the currently selected index exactly matches the preferred one.
! (b) the time gap between the preferred time stamp and the currently selected 
! index is at least as small as any other gap of consecutive prior time stamps.
!\\
!\\
! !INTERFACE:
!
  FUNCTION IsClosest ( prefYMDh, availYMDh, nTime, ctidx1 ) RESULT ( Closest )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: prefYMDh 
    INTEGER, INTENT(IN)  :: availYMDh(nTime)
    INTEGER, INTENT(IN)  :: nTime
    INTEGER, INTENT(IN)  :: ctidx1
!
! !OUTPUT PARAMETERS:
!
    LOGICAL              :: Closest
!
! !REVISION HISTORY:
!  03 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER :: N
    INTEGER :: diff, idiff

    !=================================================================
    ! IsClosest begins here! 
    !=================================================================

    ! Init
    Closest = .TRUE.

    ! It's not closest if index is not defined
    IF ( ctidx1 <= 0 ) THEN
       Closest = .FALSE.
       RETURN
    ENDIF

    ! It's closest if it is the first index
    IF ( ctidx1 == 1 ) RETURN

    ! It's closest if it matches date exactly
    IF ( availYMDh(ctidx1) == prefYMDh ) RETURN

    ! It's closest if current select one is in the future
    IF ( availYMDh(ctidx1) > prefYMDh ) RETURN

    ! Check if any of the time stamps in the past have closer intervals
    ! than the current select time stamp to it's previous one
    diff = prefYMDh - availYMDh(ctidx1)
    DO N = 2, ctidx1
       idiff = availYMDh(N) - availYMDh(N-1)
       IF ( idiff < diff ) THEN
          Closest = .FALSE.
          RETURN
       ENDIF
    ENDDO

  END FUNCTION IsClosest 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetIndex2Interp 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetIndex2Interp ( am_I_Root, Lct, nTime, availYMDh, &
                               prefYMDh,  origYMDh,   tidx1,     &
                               tidx2,     wgt1,       wgt2,  RC   ) 
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
    TYPE(ListCont),   POINTER       :: Lct
    INTEGER,          INTENT(IN)    :: nTime
    INTEGER,          INTENT(IN)    :: availYMDh(nTime)
    INTEGER,          INTENT(IN)    :: prefYMDh
    INTEGER,          INTENT(IN)    :: origYMDh
    INTEGER,          INTENT(IN)    :: tidx1
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: tidx2
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(INOUT) :: wgt1
    REAL(sp),         INTENT(INOUT) :: wgt2
    INTEGER,          INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  02 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER             :: I, tmpYMDh
    LOGICAL             :: verb
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = 'GetIndex2Interp (hcoio_dataread_mod.F90)'

    !=================================================================
    ! GetIndex2Interp begins here
    !=================================================================

    ! Verbose mode?
    verb = HCO_IsVerb(3) 

    ! If the originally wanted datetime was below the available data
    ! range, set all weights to the first index. 
    IF ( origYMDh <= availYMDh(1) ) THEN
       tidx2 = tidx1 
       wgt1  = 1.0_sp
       wgt2  = 0.0_sp

    ! If the originally wnated datetime is beyond the available data
    ! range, set tidx2 to tidx1 but leave weights in their original 
    ! values (-1.0). The reason is that we will attempt to interpolate
    ! between a second file, which is only done if the weights are 
    ! negative. 
    ELSEIF ( origYMDh >= availYMDh(nTime) ) THEN 
       tidx2 = tidx1 

    ! No interpolation needed if there is a time slices that exactly 
    ! matches the (originally) preferred datetime.
    ELSEIF( origYMDh == availYMDh(tidx1) ) THEN
       tidx2 = tidx1 
       wgt1  = 1.0_sp
       wgt2  = 0.0_sp

    ! If we are inside the data range but none of the time slices 
    ! matches the preferred datetime, get the second time slices that
    ! shall be used for data interpolation. This not necessarily needs
    ! to be the consecutive time slice. For instance, imagine a data
    ! set that contains montlhly data for years 2005 and 2010. For
    ! Feb 2007, we would want to interpolate between Feb 2005 and Feb 
    ! 2010 data. The index tidx1 already points to Feb 2005, but the
    ! upper index tidx2 needs to be set accordingly.
    ELSE

       ! Init
       tidx2 = -1

       ! Search for a time slice in the future that has the same 
       ! month/day/hour as currently selected time slice.
       tmpYMDh = availYMDh(tidx1)
       DO 
          ! Increase by one year
          tmpYMDh = tmpYMDh + 1000000
 
          ! Exit if we are beyond available dates
          IF ( tmpYMDh > availYMDh(nTime) ) EXIT
 
          ! Check if there is a time slice with that date
          DO I = tidx1,nTime
             write(*,*) 'comparing against ', availYMDh(I)
             IF ( tmpYMDh == availYMDh(I) ) THEN
                write(*,*) 'match!!'
                tidx2 = I
                EXIT
             ENDIF
          ENDDO
          IF ( tidx2 > 0 ) EXIT
       ENDDO 

       ! Repeat above but now only modify month. 
       IF ( tidx2 < 0 ) THEN
          tmpYMDh = availYMDh(tidx1)
          DO 
             ! Increase by one month
             tmpYMDh = tmpYMDh + 10000
           
             ! Exit if we are beyond available dates
             IF ( tmpYMDh > availYMDh(nTime) ) EXIT
    
             ! Check if there is a time slice with that date
             DO I = tidx1,nTime
                IF ( tmpYMDh == availYMDh(I) ) THEN
                   tidx2 = I
                   EXIT
                ENDIF
             ENDDO
             IF ( tidx2 > 0 ) EXIT
          ENDDO 
       ENDIF

       ! Repeat above but now only modify day 
       IF ( tidx2 < 0 ) THEN
          tmpYMDh = availYMDh(tidx1)
          DO 
             ! Increase by one day
             tmpYMDh = tmpYMDh + 100
           
             ! Exit if we are beyond available dates
             IF ( tmpYMDh > availYMDh(nTime) ) EXIT
    
             ! Check if there is a time slice with that date
             DO I = tidx1,nTime
                IF ( tmpYMDh == availYMDh(I) ) THEN
                   tidx2 = I
                   EXIT
                ENDIF
             ENDDO
             IF ( tidx2 > 0 ) EXIT
          ENDDO 
       ENDIF          

       ! If all of those tests failed, simply get the next time
       ! slice. 
       IF ( tidx2 < 0 ) THEN
          tidx2 = tidx1 + 1

          ! Prompt warning
          WRITE(MSG,*) 'Having problems in finding the next time slice ', &
                'to interpolate from, just take the next available ',     &
                'slice. Interpolation will be performed from ',           &
                availYMDh(tidx1), ' to ', availYMDh(tidx2), '. Data ',    &
                'container: ', TRIM(Lct%Dct%cName)
          CALL HCO_WARNING(MSG, RC, WARNLEV=1, THISLOC=LOC)
       ENDIF
       
       ! Calculate weights wgt1 and wgt2 to be given to slice 1 and 
       ! slice2, respectively.
       CALL GetWeights ( availYMDh(tidx1), availYMDh(tidx2), origYMDh, wgt1, wgt2 ) 

    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetIndex2Interp 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetWeights 
!
! !DESCRIPTION: Helper function to get the interpolation weights between
! two datetime intervals (int1, int2) and for a given time cur. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetWeights ( int1, int2, cur, wgt1, wgt2 )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   )   :: int1, int2, cur 
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(  OUT)   :: wgt1, wgt2 
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER               :: diff1, diff2 

    !=================================================================
    ! GetWeights begins here! 
    !=================================================================

    ! Check if outside of range
    IF ( cur <= int1 ) THEN
       wgt1 = 1.0_sp
    ELSEIF ( cur >= int2 ) THEN
       wgt1 = 0.0_sp
    ELSE
       diff1 = int2 - cur 
       diff2 = int2 - int1 
       wgt1  = REAL(diff1,kind=sp) / REAL(diff2,kind=sp)
    ENDIF

    ! second weight is just complement of wgt1
    wgt2  = 1.0_sp - wgt1 

  END SUBROUTINE GetWeights 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: YMDh2hrs
!
! !DESCRIPTION: returns the hours of element YMDh. For simplicity, 30 days are
! assigned to every month. At the moment, this routine is only called to
! determine the time interval between two emission time slices (DeltaT) and 
! this approximation is good enough.
!\\
!\\
! !INTERFACE:
!
  FUNCTION YMDh2hrs ( YMDh ) RESULT ( hrs ) 
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: YMDh
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER              :: hrs
!
! !REVISION HISTORY:
!  26 Jan 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! YMDh2hrs begins here! 
    !=================================================================

    hrs = FLOOR( MOD(YMDh, 10000000000) / 1.0d6 ) * 8760 + &
          FLOOR( MOD(YMDh, 1000000    ) / 1.0d4 ) * 720  + &
          FLOOR( MOD(YMDh, 10000      ) / 1.0d2 ) * 24   + &
          FLOOR( MOD(YMDh, 100        ) / 1.0d0 )

  END FUNCTION YMDh2hrs 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Normalize_Area 
!
! !DESCRIPTION: Subroutine Normalize\_Area normalizes the given array
! by the surface area calculated from the given netCDF file. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Normalize_Area( HcoState, Array, nlon, LatEdge, FN, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER         :: HcoState           ! HEMCO state object
    INTEGER,          INTENT(IN   )   :: nlon               ! # of lon midpoints
    REAL(hp),         POINTER         :: LatEdge(:)         ! lat edges 
    CHARACTER(LEN=*), INTENT(IN   )   :: FN                 ! filename
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),         POINTER         :: Array(:,:,:,:) ! Data
    INTEGER,          INTENT(INOUT)   :: RC             ! Return code
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    REAL(hp)              :: DLAT, AREA
    REAL(dp)              :: PI_180
    INTEGER               :: NLAT, J
    CHARACTER(LEN=255)    :: MSG, LOC

    !=================================================================
    ! NORNALIZE_AREA begins here! 
    !=================================================================

    ! Initialize
    LOC    = 'NORMALIZE_AREA (hcoio_dataread_mod.F90 )'
    PI_180 = HcoState%Phys%PI / 180.0_dp

    ! Check array size
    NLAT = SIZE(LatEdge,1) - 1
    
    IF ( SIZE(Array,1) /= nlon ) THEN
       MSG = 'Array size does not agree with nlon: ' // TRIM(FN)
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( SIZE(Array,2) /= NLAT ) THEN
       MSG = 'Array size does not agree with nlat: ' // TRIM(FN)
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Loop over all latitudes
    DO J = 1, NLAT
       ! get grid box area in m2 for grid box with lower and upper latitude llat/ulat:
       ! Area = 2 * PI * Re^2 * DLAT / nlon, where DLAT = abs( sin(ulat) - sin(llat) ) 
       DLAT = ABS( SIN(LatEdge(J+1)*PI_180) - SIN(LatEdge(J)*PI_180) )
       AREA = ( 2_hp * HcoState%Phys%PI * DLAT * HcoState%Phys%Re**2 ) / REAL(nlon,hp)

       ! convert array data to m-2
       ARRAY(:,J,:,:) = ARRAY(:,J,:,:) / AREA 
    ENDDO

    ! Prompt a warning
    WRITE(MSG,*) 'No area unit found in ' // TRIM(FN) // ' - convert to m-2!'
    CALL HCO_WARNING ( MSG, RC, WARNLEV=1, THISLOC=LOC )

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Normalize_Area
!EOC
#endif
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_ReadOther
!
! !DESCRIPTION: Subroutine HCOIO\_ReadOther is a wrapper routine to
! read data from sources other than netCDF.
!\\
!\\
! If a file name is given (ending with '.txt'), the data are assumed
! to hold country-specific values (e.g. diurnal scale factors). In all
! other cases, the data is directly read from the configuration file
! (scalars).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_ReadOther ( am_I_Root, HcoState, Lct, RC ) 
!
! !USES:
!
!
! !INPUT PARAMTERS:
!
    LOGICAL,         INTENT(IN   )    :: am_I_Root
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC 
!
! !REVISION HISTORY:
!  22 Dec 2014 - C. Keller: Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: MSG

    !======================================================================
    ! HCOIO_ReadOther begins here
    !======================================================================
  
    ! Error check: data must be in local time 
    IF ( .NOT. Lct%Dct%Dta%IsLocTime ) THEN
       MSG = 'Cannot read data from file that is not in local time: ' // &
             TRIM(Lct%Dct%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC='HCOIO_ReadOther (hcoio_dataread_mod.F90)' )
       RETURN
    ENDIF

    ! Read an ASCII file as country values
    IF ( INDEX( TRIM(Lct%Dct%Dta%ncFile), '.txt' ) > 0 ) THEN
       CALL HCOIO_ReadCountryValues ( am_I_Root, HcoState, Lct, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ! Directly read from configuration file otherwise
    ELSE
       CALL HCOIO_ReadFromConfig ( am_I_Root, HcoState, Lct, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_ReadOther
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_ReadCountryValues
!
! !DESCRIPTION: Subroutine HCOIO\_ReadCountryValues
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_ReadCountryValues ( am_I_Root, HcoState, Lct, RC ) 
!
! !USES:
!
    USE inquireMod,         ONLY : findFreeLUN
    USE HCO_CHARTOOLS_MOD,  ONLY : HCO_CMT, HCO_SPC, NextCharPos
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE HCO_FileData_Mod,   ONLY : FileData_ArrCheck
!
! !INPUT PARAMTERS:
!
    LOGICAL,         INTENT(IN   )    :: am_I_Root
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC 
!
! !REVISION HISTORY:
!  22 Dec 2014 - C. Keller: Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: IUFILE, IOS
    INTEGER               :: ID1, ID2, I, NT, CID, NLINE
    REAL(sp), POINTER     :: CNTR(:,:) => NULL()
    INTEGER,  ALLOCATABLE :: CIDS(:,:)
    REAL(hp), POINTER     :: Vals(:) => NULL()
    LOGICAL               :: Verb
    CHARACTER(LEN=2047)   :: LINE
    CHARACTER(LEN=255)    :: MSG, DUM, CNT
    CHARACTER(LEN=255)    :: LOC = 'HCOIO_ReadCountryValues (hcoio_dataread_mod.F90)'

    !======================================================================
    ! HCOIO_ReadCountryValues begins here
    !======================================================================
   
    ! verbose mode? 
    Verb = HCO_IsVerb(2) 
   
    ! Verbose
    IF ( Verb ) THEN
       MSG = 'Use country-specific values for ' // TRIM(Lct%Dct%cName)
       CALL HCO_MSG(MSG)
       MSG = '- Source file: ' // TRIM(Lct%Dct%Dta%ncFile)
       CALL HCO_MSG(MSG)
    ENDIF

    ! Open file
    IUFILE = FindFreeLun()
    OPEN ( IUFILE, FILE=TRIM( Lct%Dct%Dta%ncFile ), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       MSG = 'Cannot open ' // TRIM(Lct%Dct%Dta%ncFile)
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN 
    ENDIF 

    ! Repeat for every line
    NLINE = 0
    DO

       ! Read line
       READ( IUFILE, '(a)', IOSTAT=IOS ) LINE
    
       ! End of file?
       IF ( IOS < 0 ) EXIT

       ! Error?
       IF ( IOS > 0 ) THEN
          MSG = 'Error reading ' // TRIM(Lct%Dct%Dta%ncFile)
          MSG = TRIM(MSG) // ' - last valid line: ' // TRIM(LINE)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Skip commented lines and/or empty lines
       IF ( TRIM(LINE) == '' ) CYCLE
       IF ( LINE(1:1) == HCO_CMT() ) CYCLE

       ! First (valid) line holds the name of the mask container
       IF ( NLINE == 0 ) THEN

          ! Get pointer to mask. Convert to integer
          CALL HCO_GetPtr ( am_I_Root, TRIM(LINE), CNTR, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          ALLOCATE( CIDS(HcoState%NX, HcoState%NY), STAT=IOS )
          IF ( IOS /= 0 ) THEN
             CALL HCO_ERROR( 'Cannot allocate CIDS', RC, THISLOC=LOC )
             RETURN
          ENDIF
          CIDS = NINT(CNTR)

          ! Verbose
          IF ( HCO_IsVerb(3) ) THEN
             MSG = '- Use ID mask ' // TRIM(LINE)
             CALL HCO_MSG(MSG)
          ENDIF

          ! Go to next line
          NLINE = NLINE + 1
          CYCLE
       ENDIF

       ! Get first space character to skip country name.
       ! We assume here that a country name is given right at the
       ! beginning of the line, e.g. 'USA 744 1.05/1.02/...'
       ID1 = NextCharPos( LINE, HCO_SPC() )
       CNT = LINE(1:ID1)

       ! Get country ID
       DO I = ID1, LEN(LINE)
          IF ( LINE(I:I) /= HCO_SPC() ) EXIT
       ENDDO
       ID1 = I
       ID2 = NextCharPos( LINE, HCO_SPC(), START=ID1 )

       IF ( ID2 >= LEN(LINE) .OR. ID2 < 0 ) THEN
          MSG = 'Cannot extract country ID from: ' // TRIM(LINE)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       DUM = LINE(ID1:ID2)
       READ( DUM, * ) CID

       ! Extract data values
       ID1  = ID2+1
       ID2  = LEN(LINE)
       LINE = LINE(ID1:ID2)
       CALL GetDataVals( am_I_Root, HcoState, Lct, LINE, Vals, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Check data / array dimensions
       NT = SIZE(Vals,1)
       CALL FileData_ArrCheck( Lct%Dct%Dta, HcoState%NX, HcoState%NY, NT, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Pass to data array. If the country ID is larger than zero, fill 
       ! only those grid boxes. Otherwise, fill all grid boxes that have 
       ! not yet been filled. 
       DO I = 1, NT
          IF ( CID == 0 ) THEN
             WHERE ( Lct%Dct%Dta%V2(I)%Val <= 0.0_sp )
                Lct%Dct%Dta%V2(I)%Val = Vals(I)
             ENDWHERE
          ELSE
             WHERE ( CIDS == CID )
                Lct%Dct%Dta%V2(I)%Val = Vals(I)
             ENDWHERE
          ENDIF
       ENDDO

       ! Verbose
       IF ( HCO_IsVerb(3) ) THEN
          WRITE(MSG,*) '- Obtained values for ',TRIM(CNT),' ==> ID:', CID
          CALL HCO_MSG(MSG)
       ENDIF
 
       ! Cleanup
       IF ( ASSOCIATED(Vals) ) DEALLOCATE( Vals )
       Vals => NULL()

       ! Update # of read lines
       NLINE = NLINE + 1
    ENDDO

    ! Close file
    CLOSE ( IUFILE )

    ! Data is 2D
    Lct%Dct%Dta%SpaceDim  = 2

    ! Make sure data is in local time
    IF ( .NOT. Lct%Dct%Dta%IsLocTime ) THEN
       Lct%Dct%Dta%IsLocTime = .TRUE.
       MSG = 'Data assigned to mask regions will be treated in local time: '//&
              TRIM(Lct%Dct%cName)
       CALL HCO_WARNING( MSG, RC, WARNLEV=2, THISLOC=LOC )
    ENDIF

    ! Cleanup
    Cntr => NULL()
    IF ( ALLOCATED(CIDS) ) DEALLOCATE ( CIDS )

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_ReadCountryValues
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_ReadFromConfig 
!
! !DESCRIPTION: Subroutine HCOIO\_ReadFromConfig reads data directly from 
! the configuration file (instead of reading it from a netCDF file).
! These data is always assumed to be spatially uniform, but it is possible
! to specify multiple time slices by separating the individual time slice
! values by the HEMCO separator sign ('/' by default). The time dimension
! of these data is either determined from the srcTime attribute or estimated
! from the number of time slices provided. For example, if no srcTime is 
! specified and 24 time slices are provided, data is assumed to represent
! hourly data. Similarly, data is assumed to represent weekdaily or monthly
! data for 7 or 12 time slices, respectively.
!\\
!\\
! If the srcTime attribute is defined, the time slices are determined from
! this attribute. Only one time dimension (year, month, day, or hour) can
! be defined for scalar fields!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_ReadFromConfig ( am_I_Root, HcoState, Lct, RC ) 
!
! !USES:
!
    USE HCO_FILEDATA_MOD,   ONLY : FileData_ArrCheck
!
! !INPUT PARAMTERS:
!
    LOGICAL,         INTENT(IN   )    :: am_I_Root
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC 
!
! !REVISION HISTORY:
!  24 Jul 2014 - C. Keller: Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, NT
    REAL(hp), POINTER  :: Vals(:) => NULL()
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'HCOIO_ReadFromConfig (hcoio_dataread_mod.F90)'

    !======================================================================
    ! HCOIO_ReadFromConfig begins here
    !======================================================================
   
    ! Verbose
    IF ( HCO_IsVerb(2) ) THEN
       WRITE(MSG, *) 'Read from config file: ', TRIM(Lct%Dct%cName)
       CALL HCO_MSG(MSG)
    ENDIF

    !-------------------------------------------------------------------
    ! Get data values for this time step.
    !-------------------------------------------------------------------
    CALL GetDataVals ( am_I_Root, HcoState, Lct, Lct%Dct%Dta%ncFile, Vals, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
 
    !-------------------------------------------------------------------
    ! Copy data into array.
    !-------------------------------------------------------------------

    ! Number of values
    NT = SIZE(Vals,1)

    ! For masks, interpret data as mask corners (lon1/lat1/lon2/lat2) 
    ! with no time dimension 
    IF ( Lct%Dct%DctType == HCO_DCTTYPE_MASK ) THEN

       ! Make sure data is allocated
       CALL FileData_ArrCheck( Lct%Dct%Dta, HcoState%NX, HcoState%NY, 1, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Fill array: 1.0 within grid box, 0.0 outside.
       CALL FillMaskBox ( am_I_Root, HcoState, Lct, Vals, RC )

       ! Data is 2D
       Lct%Dct%Dta%SpaceDim = 2

    ! For base emissions and scale factors, interpret data as scalar 
    ! values with a time dimension.
    ELSE

       CALL FileData_ArrCheck( Lct%Dct%Dta, 1, 1, NT, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       DO I = 1, NT
          Lct%Dct%Dta%V2(I)%Val(1,1) = Vals(I)
       ENDDO

       ! Data is 1D
       Lct%Dct%Dta%SpaceDim  = 1

       ! Make sure data is in local time
       IF ( .NOT. Lct%Dct%Dta%IsLocTime ) THEN
          Lct%Dct%Dta%IsLocTime = .TRUE.
          MSG = 'Scale factors read from file are treated as local time: '// &
                 TRIM(Lct%Dct%cName)
          CALL HCO_WARNING( MSG, RC, WARNLEV=2, THISLOC=LOC )
       ENDIF

    ENDIF

    ! Cleanup
    IF ( ASSOCIATED(Vals) ) DEALLOCATE(Vals)

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_ReadFromConfig 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetDataVals 
!
! !DESCRIPTION: Subroutine GetDataVals extracts the data values from ValStr
! and writes them into vector Vals. ValStr is typically a character string
! read from an external ASCII file or directly from the HEMCO configuration
! file. Depending on the time specifications provided in the configuration
! file, Vals will be filled with only a subset of the values of ValStr.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetDataVals ( am_I_Root, HcoState, Lct, ValStr, Vals, RC ) 
!
! !USES:
!
    USE HCO_CHARTOOLS_MOD,  ONLY : HCO_CharSplit
    USE HCO_CHARTOOLS_MOD,  ONLY : HCO_WCD, HCO_SEP
    USE HCO_UNIT_MOD,       ONLY : HCO_Unit_Change
    USE HCO_tIdx_Mod,       ONLY : HCO_GetPrefTimeAttr
!
! !INPUT PARAMTERS:
!
    LOGICAL,          INTENT(IN   )   :: am_I_Root
    TYPE(HCO_State),  POINTER         :: HcoState    ! HEMCO state
    CHARACTER(LEN=*), INTENT(IN   )   :: ValStr
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC 
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         POINTER         :: Vals(:)
!
! !REVISION HISTORY:
!  22 Dec 2014 - C. Keller: Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)           :: MW_g,  EmMW_g, MolecRatio
    INTEGER            :: HcoID
    INTEGER            :: I, N, NUSE, AS
    INTEGER            :: IDX1, IDX2
    INTEGER            :: AreaFlag, TimeFlag, Check
    INTEGER            :: prefYr, prefMt, prefDy, prefHr
    REAL(hp)           :: UnitFactor 
    REAL(hp)           :: FileVals(100)
    REAL(hp), POINTER  :: FileArr(:,:,:,:) => NULL()
    LOGICAL            :: IsPerArea
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'GetDataVals (hcoio_dataread_mod.F90)'

    !======================================================================
    ! GetDataVals begins here
    !======================================================================
   
    ! Shadow molecular weights and molec. ratio (needed for
    ! unit conversion during file read)
    HcoID = Lct%Dct%HcoID
    IF ( HcoID > 0 ) THEN
       MW_g       = HcoState%Spc(HcoID)%MW_g
       EmMW_g     = HcoState%Spc(HcoID)%EmMW_g
       MolecRatio = HcoState%Spc(HcoID)%MolecRatio
    ELSE
       MW_g       = -999.0_hp 
       EmMW_g     = -999.0_hp 
       MolecRatio = -999.0_hp
    ENDIF

    ! Read data into array
    CALL HCO_CharSplit ( ValStr, HCO_SEP(), HCO_WCD(), FileVals, N, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ error if no scale factor defined
    IF ( N == 0 ) THEN
       MSG = 'Cannot read data: ' // TRIM(Lct%Dct%cName) // &
             ': ' // TRIM(ValStr)
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC)
       RETURN 
    ENDIF

    ! ---------------------------------------------------------------- 
    ! For masks, assume that values represent the corners of the mask
    ! box, e.g. there must be four values. Masks are time-independent
    ! and unitless
    ! ---------------------------------------------------------------- 
    IF ( Lct%Dct%DctType == HCO_DCTTYPE_MASK ) THEN

       ! There must be exactly four values
       IF ( N /= 4 ) THEN
          MSG = 'Mask values are not lon1/lat1/lon2/lat2: ' // &
                TRIM(ValStr) // ' --> ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Pass to FileArr array (will be used below)
       NUSE = 4
       ALLOCATE( FileArr(1,1,1,NUSE), STAT=AS )
       IF ( AS /= 0 ) THEN
          MSG = 'Cannot allocate FileArr'
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       FileArr(1,1,1,:) = FileVals(1:NUSE) 

    ! ----------------------------------------------------------------
    ! For non-masks, the data is interpreted as uniform values with
    ! a time dimension. Need to select the time slices to be used at
    ! this time (depending on the provided time attributes), as well
    ! as to ensure that values are in the correct units.
    ! Use all time slices unless a time interval is provided in
    ! attribute srcTime of the configuration file.
    ! ---------------------------------------------------------------- 
    ELSE

       ! Get the preferred times, i.e. the preferred year, month, day, 
       ! or hour (as specified in the configuration file).
       CALL HCO_GetPrefTimeAttr ( Lct, prefYr, prefMt, prefDy, prefHr, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
   
       ! Currently, data read directly from the configuration file can only
       ! represent one time dimension, i.e. it can only be yearly, monthly,
       ! daily (or hourly data, but this is read all at the same time). 
   
       ! Annual data 
       IF ( Lct%Dct%Dta%ncYrs(1) /= Lct%Dct%Dta%ncYrs(2) ) THEN
          ! Error check
          IF ( Lct%Dct%Dta%ncMts(1) /= Lct%Dct%Dta%ncMts(2) .OR. & 
               Lct%Dct%Dta%ncDys(1) /= Lct%Dct%Dta%ncDys(2) .OR. & 
               Lct%Dct%Dta%ncHrs(1) /= Lct%Dct%Dta%ncHrs(2)       ) THEN
             MSG = 'Data must not have more than one time dimension: ' // &
                    TRIM(Lct%Dct%cName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
   
          CALL GetSliceIdx ( Lct, 1, prefYr, IDX1, RC ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
          IDX2 = IDX1
          NUSE = 1
   
       ! Monthly data
       ELSEIF ( Lct%Dct%Dta%ncMts(1) /= Lct%Dct%Dta%ncMts(2) ) THEN
          ! Error check
          IF ( Lct%Dct%Dta%ncDys(1) /= Lct%Dct%Dta%ncDys(2) .OR. & 
               Lct%Dct%Dta%ncHrs(1) /= Lct%Dct%Dta%ncHrs(2)       ) THEN
             MSG = 'Data must only have one time dimension: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
   
          CALL GetSliceIdx ( Lct, 2, prefMt, IDX1, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IDX2 = IDX1
          NUSE = 1
   
       ! Daily data
       ELSEIF ( Lct%Dct%Dta%ncDys(1) /= Lct%Dct%Dta%ncDys(2) ) THEN
          ! Error check
          IF ( Lct%Dct%Dta%ncHrs(1) /= Lct%Dct%Dta%ncHrs(2) ) THEN
             MSG = 'Data must only have one time dimension: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
   
          CALL GetSliceIdx ( Lct, 3, prefDy, IDX1, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IDX2 = IDX1
          NUSE = 1
   
       ! All other cases (incl. hourly data): read all time slices).
       ELSE
          IDX1 = 1
          IDX2 = N
          NUSE = N
       ENDIF
   
       ! ---------------------------------------------------------------- 
       ! Read selected time slice(s) into data array
       ! ----------------------------------------------------------------
       IF ( IDX2 > N ) THEN
          WRITE(MSG,*) 'Index ', IDX2, ' is larger than number of ', &
                       'values found: ', TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ALLOCATE( FileArr(1,1,1,NUSE), STAT=AS )
       IF ( AS /= 0 ) THEN
          MSG = 'Cannot allocate FileArr'
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    
       ! IDX1 becomes -1 for data that is outside of the valid range
       ! (and no time cycling enabled). In this case, make sure that
       ! scale factor is set to zero.
       IF ( IDX1 < 0 ) THEN
          IF ( Lct%Dct%DctType == HCO_DCTTYPE_BASE ) THEN
             FileArr(1,1,1,:) = 0.0_hp
             MSG = 'Base field outside of range - set to zero: ' // &
                   TRIM(Lct%Dct%cName)
             CALL HCO_WARNING ( MSG, RC, WARNLEV=1, THISLOC=LOC )
          ELSE
             FileArr(1,1,1,:) = 1.0_hp
             MSG = 'Scale factor outside of range - set to one: ' // &
                   TRIM(Lct%Dct%cName)
             CALL HCO_WARNING ( MSG, RC, WARNLEV=1, THISLOC=LOC )
          ENDIF
       ELSE
          FileArr(1,1,1,:) = FileVals(IDX1:IDX2)
       ENDIF
   
       ! ---------------------------------------------------------------- 
       ! Convert data to HEMCO units 
       ! ---------------------------------------------------------------- 
       CALL HCO_UNIT_CHANGE( Array         = FileArr,                    &
                             Units         = TRIM(Lct%Dct%Dta%OrigUnit), &
                             MW_IN         = MW_g,                       &
                             MW_OUT        = EmMW_g,                     &
                             MOLEC_RATIO   = MolecRatio,                 &
                             YYYY          = -999,                       &
                             MM            = -999,                       &
                             AreaFlag      = AreaFlag,                   &
                             TimeFlag      = TimeFlag,                   &
                             FACT          = UnitFactor,                 &
                             RC            = RC                           )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! testing only
       IF ( UnitFactor /= 1.0_hp ) THEN
             WRITE(MSG,*) 'Data was in units of ', TRIM(Lct%Dct%Dta%OrigUnit), &
                          ' - converted to HEMCO units by applying ', &
                          'scale factor ', UnitFactor
             write(*,*) TRIM(MSG) 
       ENDIF
 
       ! Verbose mode
       IF ( UnitFactor /= 1.0_hp ) THEN
          IF ( HCO_IsVerb(1) ) THEN
             WRITE(MSG,*) 'Data was in units of ', TRIM(Lct%Dct%Dta%OrigUnit), &
                          ' - converted to HEMCO units by applying ', &
                          'scale factor ', UnitFactor
             CALL HCO_MSG(MSG)
          ENDIF
       ELSE
          IF ( HCO_IsVerb(2) ) THEN
             WRITE(MSG,*) 'Data was in units of ', TRIM(Lct%Dct%Dta%OrigUnit), &
                          ' - unit conversion factor is ', UnitFactor 
             CALL HCO_MSG(MSG)
          ENDIF
       ENDIF
 
       ! Data must be ... 
       ! ... concentration ...
       IF ( AreaFlag == 3 .AND. TimeFlag == 0 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.
   
       ELSEIF ( AreaFlag == 3 .AND. TimeFlag == 1 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.
          FileArr = FileArr * HcoState%TS_EMIS
          MSG = 'Data converted from kg/m3/s to kg/m3: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(Lct%Dct%Dta%OrigUnit)
          CALL HCO_WARNING ( MSG, RC, WARNLEV=1, THISLOC=LOC )
   
       ! ... emissions or unitless ...
       ELSEIF ( (AreaFlag == -1 .AND. TimeFlag == -1) .OR. &
                (AreaFlag ==  2 .AND. TimeFlag ==  1)       ) THEN
          Lct%Dct%Dta%IsConc = .FALSE.
   
       ! ... invalid otherwise:
       ELSE
          MSG = 'Unit must be unitless, emission or concentration: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(Lct%Dct%Dta%OrigUnit)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
   
       ! Auto-detect delta t [in hours] between time slices.
       ! Scale factors can be:
       ! length 1 : constant
       ! length 7 : weekday factors: Sun, Mon, ..., Sat
       ! length 12: monthly factors: Jan, Feb, ..., Dec
       ! length 24: hourly  factors: 12am, 1am, ... 11pm
       IF ( NUSE == 1 ) THEN
          Lct%Dct%Dta%DeltaT = 0
       ELSEIF ( NUSE == 7 ) THEN
          Lct%Dct%Dta%DeltaT = 24
       ELSEIF ( NUSE == 12 ) THEN
          Lct%Dct%Dta%DeltaT = 720 
       ELSEIF ( NUSE == 24 ) THEN
          Lct%Dct%Dta%DeltaT = 1 
       ELSE
          MSG = 'Factor must be of length 1, 7, 12, or 24!' // &
                 TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC)
          RETURN 
       ENDIF
   
    ENDIF ! Masks vs. non-masks
   
    ! Copy data into output array.
    IF ( ASSOCIATED(Vals) ) DEALLOCATE( Vals )
    ALLOCATE( Vals(NUSE), STAT=AS )
    IF ( AS /= 0 ) THEN
       MSG = 'Cannot allocate Vals'
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    Vals(:) = FileArr(1,1,1,:)

    ! Cleanup
    IF ( ASSOCIATED(FileArr) ) DEALLOCATE(FileArr)
    FileArr => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetDataVals 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FillMaskBox
!
! !DESCRIPTION: Subroutine FillMaskBox fills the data array of the passed list 
! container Lct according to the mask region provided in Vals. Vals contains
! the mask region of interest, denoted by the lower left and upper right grid
! box corners: lon1, lat1, lon2, lat2. The data array of Lct is filled such 
! that all grid boxes are set to 1 whose mid-point is inside of the given box 
! range.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FillMaskBox ( am_I_Root, HcoState, Lct, Vals, RC )
!
! !USES:
!
!
! !INPUT PARAMTERS:
!
    LOGICAL,          INTENT(IN   )   :: am_I_Root
    TYPE(HCO_State),  POINTER         :: HcoState    ! HEMCO state
    REAL(hp)        , POINTER         :: Vals(:)
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER         :: Lct
    INTEGER,          INTENT(INOUT)   :: RC 
!
! !REVISION HISTORY:
!  29 Dec 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J
    INTEGER            :: LON1, LON2, LAT1, LAT2
    CHARACTER(LEN=255) :: LOC = 'FillMaskBox (HCOIO_DataRead_Mod.F90)'

    !=================================================================
    ! FillMaskBox begins here! 
    !=================================================================

    ! Extract lon1, lon2, lat1, lat2
    LON1 = VALS(1)
    LAT1 = VALS(2)
    LON2 = VALS(3)
    LAT2 = VALS(4)

    ! Check for every grid box if mid point is within mask region. 
    ! Set to 1.0 if this is the case.
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( I, J )    &
!$OMP SCHEDULE( DYNAMIC )
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX
    
       IF ( HcoState%Grid%XMID%Val(I,J) >= LON1 .AND. &
            HcoState%Grid%XMID%Val(I,J) <= LON2 .AND. &
            HcoState%Grid%YMID%Val(I,J) >= LAT1 .AND. &
            HcoState%Grid%YMID%Val(I,J) <= LAT2        ) THEN

          Lct%Dct%Dta%V2(1)%Val(I,J) = 1.0_sp
       ENDIF 

    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE FillMaskBox
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetSliceIdx 
!
! !DESCRIPTION: gets the time slice index to be used for data directly
! read from the HEMCO configuration file. prefDt denotes the preferred
! time attribute (year, month, or day). DtType is used to identify the 
! time attribute type (1=year, 2=month, 3=day). The time slice index will 
! be selected based upon those two variables. IDX is the selected time 
! slice index. It will be set to -1 if the current simulation date
! is outside of the specified time range and the time cycle attribute is 
! not enabled for this field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetSliceIdx ( Lct, DtType, prefDt, IDX, RC ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER                 :: Lct
    INTEGER,          INTENT(IN   )           :: DtType
    INTEGER,          INTENT(IN   )           :: prefDt
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: IDX
    INTEGER,          INTENT(INOUT)           :: RC 
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER            :: lowDt, uppDt
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'GetSliceIdx (HCOIO_DataRead_Mod.F90)'

    !=================================================================
    ! GetSliceIdx begins here! 
    !=================================================================

    ! Init
    RC = HCO_SUCCESS

    ! Get upper and lower time range
    IF ( DtType == 1 ) THEN
       lowDt = Lct%Dct%Dta%ncYrs(1)
       uppDt = Lct%Dct%Dta%ncYrs(2)
    ELSEIF ( DtType == 2 ) THEN
       lowDt = Lct%Dct%Dta%ncMts(1)
       uppDt = Lct%Dct%Dta%ncMts(2)
    ELSEIF ( DtType == 3 ) THEN
       lowDt = Lct%Dct%Dta%ncDys(1)
       uppDt = Lct%Dct%Dta%ncDys(2)
    ELSE
       WRITE(MSG,*) "DtType must be one of 1, 2, 3: ", DtType
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN 
    ENDIF

    ! Check for cycle flags:

    ! Data cycle set to range or exact date: in these cases, the 
    ! the preferred date will be equal to the current date, so 
    ! check if the preferred date is indeed within the available 
    ! range (lowDt, uppDt).
    ! For data only to be used within the specified range, set 
    ! index to -1. This will force the scale factors to be set to
    ! zero!
    IF ( prefDt < lowDt .OR. prefDt > uppDt ) THEN
       IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT ) THEN ! Exact match
          MSG = 'Data is not on exact date: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN 
       ELSEIF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE ) THEN ! w/in range
          IDX = -1
          RETURN
       ELSE
          ! this here should never happen, since for a cycle flag of 1,
          ! the preferred date should always be restricted to the range
          ! of available time stamps.
          MSG = 'preferred date is outside of range: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN 
       ENDIF
    ENDIF
    
    ! If the code makes it to here, prefDt is within the available data range
    ! and we simply get the wanted index from the current index and the lowest
    ! available index. 
    IDX = prefDt - lowDt + 1

  END SUBROUTINE GetSliceIdx
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SrcFile_Parse
!
! !DESCRIPTION: Routine SrcFile\_Parse parses the source file name ('ncFile')
! of the provided list container Lct. In particular, it searches for tokens 
! such as $ROOT, $YYYY, etc., within the file name and replaces those values 
! with the intendend characters. The parsed file name is returned in string
! srcFile, while the original file name is retained in Lct.
!\\
!\\
! It now also checks if the file exists. If the file does not exist and the
! file name contains date tokens, it tries to adjust the file name to the
! closest available date in the past. The optional flag FUTURE can be used
! to denote that the next available file in the future shall be selected,
! even if there is a file that exactly matches the preferred date time. This
! is useful for interpolation between fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SrcFile_Parse ( am_I_Root, HcoState, Lct, srcFile, FOUND, RC, FUTURE )
!
! !USES:
!
    USE HCO_TIDX_MOD,         ONLY : HCO_GetPrefTimeAttr
    USE HCO_TIDX_MOD,         ONLY : tIDx_IsInRange 
    USE HCO_CHARTOOLS_MOD,    ONLY : HCO_CharParse
    USE HCO_CLOCK_MOD,        ONLY : HcoClock_Get
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root  ! Root CPU?
    TYPE(HCO_State),  POINTER                 :: HcoState   ! HEMCO state object
    TYPE(ListCont),   POINTER                 :: Lct        ! HEMCO list container
    LOGICAL,          INTENT(IN   ), OPTIONAL :: FUTURE     ! If needed, update
                                                            ! date tokens to future 
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT)           :: srcFile    ! output string
    LOGICAL,          INTENT(  OUT)           :: FOUND      ! Does file exist?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC         ! return code
!
! !REVISION HISTORY:
!  01 Oct 2014 - C. Keller - Initial version
!  23 Feb 2015 - C. Keller - Now check for negative return values in
!                            HCO_GetPrefTimeAttr 
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER :: INC,     CNT,    TYPCNT, TYP,   NEWTYP
    INTEGER :: prefYr,  prefMt, prefDy, prefHr
    INTEGER :: origYr,  origMt, origDy, origHr
    LOGICAL :: hasFile, hasYr,  hasMt,  hasDy, hasHr
    LOGICAL :: nextTyp

    ! maximum # of iterations for file search
    INTEGER, PARAMETER :: MAXIT = 10000

    !=================================================================
    ! SrcFile_Parse
    !=================================================================

    ! Initialize to input string
    srcFile = Lct%Dct%Dta%ncFile

    ! Get preferred dates (to be passed to parser
    CALL HCO_GetPrefTimeAttr ( Lct, prefYr, prefMt, prefDy, prefHr, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Make sure dates are not negative 
    IF ( prefYr <= 0 ) THEN
       CALL HcoClock_Get( cYYYY = prefYr, RC = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF
    IF ( prefMt <= 0 ) THEN
       CALL HcoClock_Get( cMM   = prefMt, RC = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF
    IF ( prefDy <= 0 ) THEN
       CALL HcoClock_Get( cDD   = prefDy, RC = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF
    IF ( prefHr <  0 ) THEN
       CALL HcoClock_Get( cH    = prefHr, RC = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF 

    ! Call the parser
    CALL HCO_CharParse ( srcFile, prefYr, prefMt, prefDy, prefHr, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Check if file exists
    INQUIRE( FILE=TRIM(srcFile), EXIST=HasFile )

    ! If we are looking for a future file, force HasFile to be false.
    IF ( PRESENT(FUTURE) ) THEN
       IF ( FUTURE ) HasFile = .FALSE.
    ENDIF

    ! If file does not exist, check if we can adjust prefYr, prefMt, etc.
    IF ( .NOT. HasFile .AND. Lct%Dct%DctType /= HCO_CFLAG_EXACT ) THEN

       ! Check if any token exist
       HasYr = ( INDEX(TRIM(Lct%Dct%Dta%ncFile),'YYYY') > 0 )
       HasMt = ( INDEX(TRIM(Lct%Dct%Dta%ncFile),'MM'  ) > 0 )
       HasDy = ( INDEX(TRIM(Lct%Dct%Dta%ncFile),'DD'  ) > 0 )
       HasHr = ( INDEX(TRIM(Lct%Dct%Dta%ncFile),'HH'  ) > 0 )

       ! Search for file
       IF ( HasYr .OR. HasMt .OR. HasDy .OR. HasHr ) THEN

          ! Date increments
          INC = -1
          IF ( PRESENT(FUTURE) ) THEN
             IF ( FUTURE ) INC = 1
          ENDIF

          ! Initialize counters
          CNT = 0
          
          ! Type is the update type (see below)
          TYP = 0

          ! Mirror preferred variables
          origYr = prefYr
          origMt = prefMt
          origDy = prefDy
          origHr = prefHr

          ! Do until file is found or counter exceeds threshold
          DO WHILE ( .NOT. HasFile )
            
             ! Inrease counter
             CNT = CNT + 1
             IF ( CNT > MAXIT ) EXIT
 
             ! Increase update type if needed:
             nextTyp = .FALSE.
 
             ! Type 0: Initialization
             IF ( TYP == 0 ) THEN
                nextTyp = .TRUE.
             ! Type 1: update hour only
             ELSEIF ( TYP == 1 .AND. TYPCNT > 24 ) THEN
                nextTyp = .TRUE.
             ! Type 2: update day only
             ELSEIF ( TYP == 2 .AND. TYPCNT > 31 ) THEN
                nextTyp = .TRUE.
             ! Type 3: update month only
             ELSEIF ( TYP == 3 .AND. TYPCNT > 12 ) THEN
                nextTyp = .TRUE.
             ! Type 4: update year only
             ELSEIF ( TYP == 4 .AND. TYPCNT > 300 ) THEN
                nextTyp = .TRUE.
             ! Type 5: update hour and day 
             ELSEIF ( TYP == 5 .AND. TYPCNT > 744 ) THEN
                nextTyp = .TRUE.
             ! Type 6: update day and month 
             ELSEIF ( TYP == 6 .AND. TYPCNT > 372 ) THEN
                nextTyp = .TRUE.
             ! Type 7: update month and year
             ELSEIF ( TYP == 7 .AND. TYPCNT > 3600 ) THEN
                EXIT
             ENDIF

             ! Get next type
             IF ( nextTyp ) THEN
                NEWTYP = -1
                IF     ( hasHr .AND. TYP < 1 ) THEN
                   NEWTYP = 1
                ELSEIF ( hasDy .AND. TYP < 2 ) THEN
                   NEWTYP = 2
                ELSEIF ( hasMt .AND. TYP < 3 ) THEN
                   NEWTYP = 3
                ELSEIF ( hasYr .AND. TYP < 4 ) THEN
                   NEWTYP = 4
                ELSEIF ( hasDy .AND. TYP < 2 ) THEN
                   NEWTYP = 5
                ELSEIF ( hasDy .AND. TYP < 2 ) THEN
                   NEWTYP = 6
                ELSEIF ( hasDy .AND. TYP < 2 ) THEN
                   NEWTYP = 7
                ENDIF
    
                ! Exit if no other type found
                IF ( NEWTYP < 0 ) EXIT
 
                ! This is the new type, reset type counter
                TYP    = NEWTYP
                TYPCNT = 0

                ! Make sure we reset all values 
                prefYr = origYr
                prefMt = origMt
                prefDy = origDy
                prefHr = origHr

             ENDIF

             ! Update preferred datetimes
             SELECT CASE ( TYP ) 
                ! Adjust hour only
                CASE ( 1 ) 
                   prefHr = prefHr + INC                      
                ! Adjust day only
                CASE ( 2 )
                   prefDy = prefDy + INC                      
                ! Adjust month only
                CASE ( 3 )
                   prefMt = prefMt + INC                      
                ! Adjust year only
                CASE ( 4 )
                   prefYr = prefYr + INC                      
                ! Adjust hour and day 
                CASE ( 5 )
                   prefHr = prefHr + INC
                   IF ( MOD(TYPCNT,24) == 0 ) prefDy = prefDy + INC
                ! Adjust day and month 
                CASE ( 6 )
                   prefDy = prefDy + INC                      
                   IF ( MOD(TYPCNT,31) == 0 ) prefMt = prefMt + INC
                ! Adjust month and year
                CASE ( 7 )
                   prefMt = prefMt + INC                      
                   IF ( MOD(TYPCNT,12) == 0 ) prefYr = prefYr + INC
                CASE DEFAULT
                   EXIT
             END SELECT

             ! Check if we need to adjust a year/month/day/hour
             IF ( prefHr < 0 ) THEN
                prefHr = 23 
                prefDy = prefDy - 1
             ENDIF
             IF ( prefHr > 23 ) THEN
                prefHr = 0 
                prefDy = prefDy + 1
             ENDIF
             IF ( prefDy < 1  ) THEN
                prefDy = 31
                prefMt = prefMt - 1
             ENDIF
             IF ( prefDy > 31 ) THEN
                prefDy = 1 
                prefMt = prefMt + 1
             ENDIF
             IF ( prefMt < 1  ) THEN
                prefMt = 12
                prefYr = prefYr - 1
             ENDIF
             IF ( prefMt > 12 ) THEN
                prefMt = 1 
                prefYr = prefYr + 1
             ENDIF
 
             ! Mirror original file          
             srcFile = Lct%Dct%Dta%ncFile

             ! Call the parser with adjusted values
             CALL HCO_CharParse ( srcFile, prefYr, prefMt, prefDy, prefHr, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! Check if this file exists
             INQUIRE( FILE=TRIM(srcFile), EXIST=HasFile )

             ! Update counter
             TYPCNT = TYPCNT + 1
          ENDDO
       ENDIF
    ENDIF

    ! Additional check for data with a given range: make sure that the selected
    ! field is not outside of the given range
    IF ( HasFile .AND. ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE ) ) THEN
       HasFile = TIDX_IsInRange ( Lct, prefYr, prefMt, prefDy, prefHr ) 
    ENDIF

    ! Return variable
    FOUND = HasFile

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE SrcFile_Parse
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SigmaMidToEdges
!
! !DESCRIPTION: Helper routine to interpolate sigma mid point values to edges.
! A simple linear interpolation is performed.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SigmaMidToEdges ( am_I_Root, SigMid, SigEdge, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root       ! Root CPU?
    REAL(hp),         POINTER                 :: SigMid(:,:,:)   ! sigma levels
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         POINTER                 :: SigEdge(:,:,:)  ! sigma edges
    INTEGER,          INTENT(  OUT)           :: RC              ! return code
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER            :: L, AS
    INTEGER            :: nx, ny, nz
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'SigmaMidToEdges (HCOIO_DataRead_Mod.F90)'

    !=================================================================
    ! SigmaMidToEdges begins here! 
    !=================================================================

    ! Allocate space as required
    nx = SIZE(SigMid,1)
    ny = SIZE(SigMid,2)
    nz = SIZE(SigMid,3)
    IF ( ASSOCIATED(SigEdge) ) DEALLOCATE(SigEdge)
    ALLOCATE(SigEdge(nx,ny,nz+1),STAT=AS)
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( 'Allocate SigEdge', RC, THISLOC=LOC )
       RETURN
    ENDIF
    SigEdge = 0.0_hp

    ! Calculate sigma edges by linear interpolation (symmetric mid-points)
    DO L = 1, nz-1
       SigEdge(:,:,L+1) = ( SigMid(:,:,L) + SigMid(:,:,L+1) ) / 2.0_hp
    ENDDO

    ! Get outermost values:
    SigEdge(:,:,1   ) = SigMid(:,:,1 ) - ( SigEdge(:,:,2) - SigMid(:,:,1)   )
    SigEdge(:,:,nz+1) = SigMid(:,:,nz) + ( SigMid(:,:,nz) - SigEdge(:,:,nz) )

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE SigmaMidToEdges
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckMissVal 
!
! !DESCRIPTION: Checks for missing values in the passed array. Missing values
! of base emissions and masks are set to 0, missing values of scale factors
! are set to 1. 
!\\
! !INTERFACE:
!
  SUBROUTINE CheckMissVal ( Lct, Arr )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER                 :: Lct
    REAL(sp),         POINTER                 :: Arr(:,:,:,:)   
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    !=================================================================
    ! CheckMissVal begins here! 
    !=================================================================

    ! Error trap
    IF ( .NOT. ASSOCIATED(Arr) ) RETURN
 
    IF ( ANY(Arr == HCO_MISSVAL) ) THEN
       ! Base emissions
       IF ( Lct%Dct%DctType == HCO_DCTTYPE_BASE ) THEN
          WHERE(Arr == HCO_MISSVAL) Arr = 0.0_sp
       ! Scale factor
       ELSEIF ( Lct%Dct%DctType == HCO_DCTTYPE_SCAL ) THEN
          WHERE(Arr == HCO_MISSVAL) Arr = 1.0_sp
       ! Mask
       ELSEIF ( Lct%Dct%DctType == HCO_DCTTYPE_MASK ) THEN
          WHERE(Arr == HCO_MISSVAL) Arr = 0.0_sp
       ENDIF
    ENDIF

  END SUBROUTINE CheckMissVal 
!EOC
END MODULE HCOIO_DataRead_Mod
