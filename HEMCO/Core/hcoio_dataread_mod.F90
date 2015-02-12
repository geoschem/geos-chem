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
  SUBROUTINE HCOIO_DataRead( am_I_Root, HcoState, Lct, RC ) 
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
!
! !INPUT/OUTPUT PARAMETERS:
!
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
    LOGICAL                    :: verb
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

    ! Point to ESMF IMPORT object
    IMPORT => HcoState%IMPORT
    ASSERT_(ASSOCIATED(IMPORT))

    ! Verbose?
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root
    IF ( verb ) THEN
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
! !INTERFACE:
!
  SUBROUTINE HCOIO_DataRead( am_I_Root, HcoState, Lct, RC ) 
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

    USE Regrid_A2A_Mod,     ONLY : MAP_A2A
    USE HCO_INTERP_MOD,     ONLY : ModelLev_Interpolate 
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object
    TYPE(ListCont),   POINTER        :: Lct        ! HEMCO list container
!
! !INPUT/OUTPUT PARAMETERS:
!
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
!  15 Jan 2015 - C. Keller   - Now allow model level interpolation in combination
!                              with MESSy (horizontal) regridding.
!  03 Feb 2015 - C. Keller   - Moved map_a2a regridding to hco_interp_mod.F90.
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)            :: thisUnit, LevUnit, LevName
    CHARACTER(LEN=255)            :: MSG 
    CHARACTER(LEN=1023)           :: srcFile
    INTEGER                       :: NX, NY
    INTEGER                       :: NCRC, Flag, AS
    INTEGER                       :: ncLun
    INTEGER                       :: nlon,   nlat,  nlev, nTime
    INTEGER                       :: lev1,   lev2,  dir 
    INTEGER                       :: tidx1,  tidx2, ncYr, ncMt
    INTEGER                       :: HcoID
    INTEGER                       :: nlatEdge, nlonEdge
    REAL(sp), POINTER             :: ncArr(:,:,:,:)   => NULL()
    REAL(hp), POINTER             :: SigEdge(:,:,:)   => NULL()
    REAL(hp), POINTER             :: SigLev (:,:,:)   => NULL()
    REAL(hp), POINTER             :: LonMid   (:)     => NULL()
    REAL(hp), POINTER             :: LatMid   (:)     => NULL()
    REAL(hp), POINTER             :: LevMid   (:)     => NULL()
    REAL(hp), POINTER             :: LonEdge  (:)     => NULL()
    REAL(hp), POINTER             :: LatEdge  (:)     => NULL()
    LOGICAL                       :: verb
    LOGICAL                       :: IsModelLevel
    REAL(hp)                      :: MW_g, EmMW_g, MolecRatio
    INTEGER                       :: UnitTolerance
    INTEGER                       :: AreaFlag, TimeFlag 

    ! Use MESSy regridding routines?
    LOGICAL                       :: UseMESSy


    !=================================================================
    ! HCOIO_DATAREAD begins here
    !=================================================================

    ! Enter
    CALL HCO_ENTER ('HCOIO_DATAREAD (hcoio_dataread_mod.F90)' , RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    
    ! Check for verbose mode
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root

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
    CALL SrcFile_Parse ( am_I_Root, HcoState, Lct, srcFile, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    IF ( verb ) THEN
       Write(MSG,*) '- Reading file ', TRIM(srcFile)
       CALL HCO_MSG(MSG)
    ENDIF

    ! ----------------------------------------------------------------
    ! Open netCDF
    ! ----------------------------------------------------------------
    CALL NC_OPEN ( TRIM(srcFile), ncLun )

    ! ----------------------------------------------------------------
    ! Extract time slice information
    ! This determines the lower and upper time slice index (tidx1 
    ! and tidx2) to be read based upon the time slice information 
    ! extracted from the file and the time stamp settings set in the
    ! HEMCO configuration file. Multiple time slices are only selected
    ! for weekdaily data or for 'autodetected' hourly data (using the
    ! wildcard character in the configuration file time attribute).
    ! ----------------------------------------------------------------
    CALL GET_TIMEIDX ( am_I_Root, HcoState, Lct,     &
                       ncLun,     tidx1,    tidx2,   &
                       ncYr,      ncMt,     RC        )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Check for negative tidx1. tidx1 can still be negative if: 
    ! (a) CycleFlag is set to 2 and the current simulation time is 
    ! outside of the data time range. In this case, we prompt a 
    ! warning and make sure that there is no data associated with
    ! this FileData container.
    ! (b) CycleFlag is set to 3 and none of the data time stamps 
    ! matches the current simulation time exactly. Return with 
    ! error!
    !-----------------------------------------------------------------
    IF ( tidx1 < 0 ) THEN
       IF ( Lct%Dct%Dta%CycleFlag == 3 ) THEN
          MSG = 'Exact time not found in ' // TRIM(srcFile) 
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ELSEIF ( Lct%Dct%Dta%CycleFlag == 1 ) THEN
          MSG = 'Invalid time index: ' // TRIM(srcFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ELSEIF ( Lct%Dct%Dta%CycleFlag == 2 ) THEN
          CALL FileData_Cleanup( Lct%Dct%Dta, DeepClean=.FALSE.)
          MSG = 'Simulation time is outside of time range provided for '//&
               TRIM(Lct%Dct%cName) // ' - data is ignored!'
          CALL HCO_WARNING ( MSG, RC )
          CALL NC_CLOSE ( ncLun ) 
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
                      RC      = NCRC                 )

    IF ( NCRC /= 0 ) THEN
       CALL HCO_ERROR( 'NC_READ_ARRAY', RC )
       RETURN 
    ENDIF

    !-----------------------------------------------------------------
    ! Convert to HEMCO units 
    !-----------------------------------------------------------------

    ! Convert to HEMCO units. This is kg/m2/s for fluxes and kg/m3 
    ! for concentrations.
    ! The srcUnit attribute of the configuration file determines to
    ! which fields unit conversion is applied:

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
          MSG = 'Data does not appear to be unitless: ' // &
                TRIM(thisUnit) // '. File: ' // TRIM(srcFile)
          CALL HCO_WARNING( MSG, RC )
       ENDIF

    ! Convert to HEMCO units in all other cases. 
    ELSE

       ! For zero unit tolerance, make sure that thisUnit matches 
       ! with unit set in configuration file!
       ! Otherwise, prompt at least a warning.
       IF ( TRIM(Lct%Dct%Dta%OrigUnit) /= TRIM(thisUnit) ) THEN
          MSG = 'File units do not match: ' // TRIM(thisUnit) // &
                ' vs. ' // TRIM(Lct%Dct%Dta%OrigUnit)    // &
                '. File: ' // TRIM(srcFile)

          IF ( UnitTolerance == 0 ) THEN
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ELSE
             CALL HCO_WARNING( MSG, RC )
          ENDIF
       ENDIF 

       ! Mirror species properties needed for unit conversion
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

       ! Now convert to HEMCO units. This only attempts to convert
       ! mass, area/volume and time to HEMCO standards (kg, m2/m3, s).
       CALL HCO_UNIT_CHANGE(                &
            Array         = ncArr,          &
            Units         = thisUnit,       &
            MW_IN         = MW_g,           & 
            MW_OUT        = EmMW_g,         & 
            MOLEC_RATIO   = MolecRatio,     & 
            YYYY          = ncYr,           &
            MM            = ncMt,           &
            AreaFlag      = AreaFlag,       &
            TimeFlag      = TimeFlag,       &
            RC            = RC               )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot convert units for ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG , RC )
          RETURN 
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
          CALL HCO_WARNING( MSG, RC )
 
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
          CALL HCO_WARNING( MSG, RC )

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
       IF ( verb ) THEN
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
       IF ( verb ) THEN
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
    CALL NC_CLOSE ( ncLun )
      
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
! Also returns the time slice year and month, as these values may be
! used for unit conversion! 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_TimeIdx( am_I_Root, HcoState, Lct,     &
                          ncLun,     tidx1,    tidx2,   &
                          ncYr,      ncMt,     RC        )
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
    INTEGER,          INTENT(  OUT)  :: ncYr      ! time slice year
    INTEGER,          INTENT(  OUT)  :: ncMt      ! time slice month
!
! !INPUT/OUTPUT PARAMETERS:
!
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
    INTEGER               :: nTime,  T, CNT, NCRC 
    INTEGER               :: prefYr, prefMt, prefDy, prefHr
    INTEGER               :: refYear
    INTEGER               :: prefYMDh
    INTEGER, POINTER      :: availYMDh(:) => NULL() 
    LOGICAL               :: verb

    !=================================================================
    ! GET_TIMEIDX begins here
    !=================================================================

    ! Init 
    CALL HCO_ENTER ('GET_TIMEIDX (hco_dataread_mod.F90)', RC )
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

    ! ---------------------------------------------------------------- 
    ! Get preferred time stamp to read based upon the specs set
    ! in the config. file. 
    ! This can return value -1 for prefHr, indicating that all  
    ! corresponding time slices shall be read.
    ! This call will return -1 for all dates if the simulation date is
    ! outside of the data range given in the configuration file.
    ! ---------------------------------------------------------------- 
    CALL HCO_GetPrefTimeAttr ( Lct, prefYr, prefMt, prefDy, prefHr, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Check if we are outside of provided range
    IF ( prefYr < 0 .OR. prefMt < 0 .OR. prefDy < 0 ) THEN
     
       ! This should only happen for 'range' data (cycle flag is 2). 
       IF ( Lct%Dct%Dta%CycleFlag /= 2 ) THEN
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

    ! prefYMDh is the preferred datetime
    prefYMDh = prefYr*1000000 + prefMt*10000 + &
               prefDy*100 + max(prefHr,0)

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
       tidx1 = -1
       tidx2 = -1 

       ! ------------------------------------------------------------- 
       ! Check if preferred datetime is within the range of available
       ! time slices. In this case, set tidx1 to the index of the 
       ! closest time slice that is not in the future. If CycleFlag
       ! is set to 3 (= exact match), tidx1 is only set if the file
       ! time stamp exactly matches with prefYMDh!
       ! ------------------------------------------------------------- 
       CALL Check_availYMDh ( Lct, nTime, availYMDh, prefYMDh, tidx1 )

       ! ------------------------------------------------------------- 
       ! If tidx1 couldn't be set in the call above, re-adjust 
       ! preferred year to the closest available year in the 
       ! time slices. Then repeat the check. Don't do this for exact
       ! dates. 
       ! ------------------------------------------------------------- 
       IF ( Lct%Dct%Dta%CycleFlag /= 3 ) THEN
         
          ! Adjust year, month, and day (in this order).
          CNT  = 0
          DO 
             CNT = CNT + 1
             IF ( tidx1 > 0 .OR. CNT > 3 ) EXIT

             ! Adjust prefYMDh at the given level (1=Y, 2=M, 3=D)
             CALL prefYMDh_Adjust ( nTime, availYMDh, prefYMDh, CNT )

             ! verbose mode 
             IF ( verb ) THEN
                write(MSG,'(A30,I12)') 'adjusted preferred datetime: ', prefYMDh
                CALL HCO_MSG(MSG)
             ENDIF
      
             CALL Check_availYMDh ( Lct, nTime, availYMDh, prefYMDh, tidx1 )
   
          ENDDO
       ENDIF
   
       ! ------------------------------------------------------------- 
       ! If tidx1 still isn't defined, i.e. prefYMDh is still 
       ! outside the range of availYMDh, set tidx1 to the closest
       ! available date. This must be 1 or nTime! 
       ! ------------------------------------------------------------- 
       IF ( tidx1 < 0 .AND. Lct%Dct%Dta%CycleFlag /= 3 ) THEN
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
             tidx2 = tidx1 + 6
             
          ! If there are more than 7 time slices, interpret the current
          ! selected index as sunday of the current time frame (e.g. sunday
          ! data of current month), and select the time slice index
          ! accordingly. This requires that there are at least 6 more time
          ! slices following the current one. 
          ELSE
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

       ! verbose mode 
       IF ( verb ) THEN
          WRITE(MSG,'(A30,I12)') 'selected tidx1: ', tidx1
          CALL HCO_MSG(MSG)
       ENDIF

       ! ------------------------------------------------------------- 
       ! Now need to set upper time slice index tidx2. This index
       ! is only different from tidx1 if multiple hourly slices are
       ! read (--> prefHr = -1 or -10, e.g. hour attribute in config. 
       ! file was set to wildcard character or data is in local hours). 
       ! In this case, check if there are multiple time slices for the 
       ! selected date (y/m/d).
       ! tidx2 has already been set to proper value above if it's
       ! weekday data.
       ! -------------------------------------------------------------
       IF ( tidx2 < 0 ) THEN

          ! Check for multiple hourly data
          IF ( tidx1 > 0 .AND. prefHr < 0 ) THEN
             CALL SET_TIDX2 ( nTime, availYMDH, tidx1, tidx2 )    

             ! Denote as local time if necessary
             IF ( Lct%Dct%Dta%ncHrs(1) == -10 ) THEN
                Lct%Dct%Dta%IsLocTime = .TRUE.
             ENDIF
          ELSE
             tidx2 = tidx1
          ENDIF
       ENDIF   

       ! verbose mode 
       IF ( (tidx2 /= tidx1) .AND. verb ) THEN
          WRITE(MSG,'(A30,I12)') 'selected tidx2: ', tidx2
          CALL HCO_MSG(MSG)
       ENDIF

    ! ================================================================
    ! Case 3: No time slice available. Set both indeces to zero. 
    ! ================================================================
    ELSE
       tidx1 = 0
       tidx2 = 0 
    ENDIF

    !-----------------------------------------------------------------
    ! Sanity check: if CycleFlag is set to 3, the file time stamp
    ! must exactly match the current time.
    !-----------------------------------------------------------------
    IF ( Lct%Dct%Dta%CycleFlag == 3 .AND. tidx1 > 0 ) THEN
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
    IF ( tidx2 > tidx1 ) THEN
       Lct%Dct%Dta%DeltaT = YMDh2hrs( availYMDh(tidx1+1) - availYMDh(tidx1) )
    ELSE
       Lct%Dct%Dta%DeltaT = 0
    ENDIF

    ! verbose 
    IF ( verb .and. tidx1 > 0 ) THEN
       write(MSG,'(A30,I12)') 'corresponding datetime 1: ', availYMDh(tidx1)
       CALL HCO_MSG(MSG)
       if ( tidx2 > tidx1 ) THEN
          write(MSG,'(A30,I12)') 'corresponding datetime 2: ', availYMDh(tidx2)
          CALL HCO_MSG(MSG)
       endif
       write(MSG,'(A30,I12)') 'assigned delta t [h]: ', Lct%Dct%Dta%DeltaT 
       CALL HCO_MSG(MSG)
       write(MSG,*) 'local time? ', Lct%Dct%Dta%IsLocTime
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

  END SUBROUTINE Get_TimeIdx
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
       IF ( availYMDh(I+1) > prefYMDh .AND. Lct%Dct%Dta%CycleFlag < 3 ) THEN
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
  SUBROUTINE prefYMDh_Adjust( N, availYMDh, prefYMDh, level ) 
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)     :: N
    INTEGER, INTENT(IN)     :: availYMDh(N)
    INTEGER, INTENT(IN)     :: level
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
    INTEGER  :: IDX, origYr, origMt, origDy, origHr, newAttr

    !=================================================================
    ! prefYMDh_Adjust begins here! 
    !=================================================================

    ! Are we taking the first or the last element of the available
    ! time slice?
    IF ( prefYMDh < availYMDh(1) ) THEN
       IDX = 1 
    ELSE
       IDX = N
    ENDIF

    ! Get original Yr, Mt, Dy and Hr
    origYr = FLOOR( MOD(prefYMDh, 10000000000) / 1.0d6 )
    origMt = FLOOR( MOD(prefYMDh, 1000000    ) / 1.0d4 )
    origDy = FLOOR( MOD(prefYMDh, 10000      ) / 1.0d2 )
    origHr = FLOOR( MOD(prefYMDh, 100        ) / 1.0d0 )

    ! Extract new attribute from availYMDh and insert into prefYMDh
    ! --- Year
    IF ( level == 1 ) THEN
       newAttr  = FLOOR( MOD(availYMDh(IDX),10000000000) / 1.0d6 )
       prefYMDh = newAttr * 1000000 + origMt * 10000 + origDy * 100 + origHr

    ! --- Month 
    ELSEIF ( level == 2 ) THEN
       newAttr  = FLOOR( MOD(availYMDh(IDX),1000000) / 1.0d4 )
       prefYMDh = origYr * 1000000 + newAttr * 10000 + origDy * 100 + origHr

    ! --- Day
    ELSEIF ( level == 3 ) THEN
       newAttr  = FLOOR( MOD(availYMDh(IDX),10000) / 1.0d2 )
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
    CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )

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
    Verb = HCO_VERBOSE_CHECK() .and. am_I_Root
   
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
          IF ( Verb ) THEN
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
       IF ( verb ) THEN
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
       CALL HCO_WARNING( MSG, RC, THISLOC=LOC )
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
    USE HCO_CHARTOOLS_MOD,  ONLY : HCO_CharSplit
    USE HCO_CHARTOOLS_MOD,  ONLY : HCO_WCD, HCO_SEP
    USE HCO_UNIT_MOD,       ONLY : HCO_Unit_Change
    USE HCO_tIdx_Mod,       ONLY : HCO_GetPrefTimeAttr
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
    INTEGER            :: I, I1, I2, J1, J2, NT
    REAL(hp), POINTER  :: Vals(:) => NULL()
    LOGICAL            :: Verb
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'HCOIO_ReadFromConfig (hcoio_dataread_mod.F90)'

    !======================================================================
    ! HCOIO_ReadFromConfig begins here
    !======================================================================
   
    ! verbose mode? 
    Verb = HCO_VERBOSE_CHECK() .and. am_I_Root
   
    ! Verbose
    IF ( Verb ) THEN
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
    IF ( Lct%Dct%DctType == 3 ) THEN

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
          CALL HCO_WARNING( MSG, RC, THISLOC=LOC )
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
    REAL(hp)           :: FileVals(100)
    REAL(hp), POINTER  :: FileArr(:,:,:,:) => NULL()
    LOGICAL            :: Verb, IsPerArea
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'GetDataVals (hcoio_dataread_mod.F90)'

    !======================================================================
    ! GetDataVals begins here
    !======================================================================
   
    ! verbose mode? 
    Verb = HCO_VERBOSE_CHECK() .and. am_I_Root
   
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
    IF ( Lct%Dct%DctType == 3 ) THEN

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
          FileArr(1,1,1,:) = 0.0_hp
          MSG = 'Scale factor outside of range - set to zero: ' // &
                TRIM(Lct%Dct%cName)
          CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
       ELSE
          FileArr(1,1,1,:) = FileVals(IDX1:IDX2)
       ENDIF
   
       ! ---------------------------------------------------------------- 
       ! Convert data to HEMCO units 
       ! ---------------------------------------------------------------- 
       CALL HCO_Unit_Change( Array       = FileArr,                    &
                             Units       = TRIM(Lct%Dct%Dta%OrigUnit), &
                             MW_IN       = MW_g,                       &
                             MW_OUT      = EmMW_g,                     &
                             MOLEC_RATIO = MolecRatio,                 &
                             YYYY        = -999,                       &
                             MM          = -999,                       &
                             AreaFlag    = AreaFlag,                   &
                             TimeFlag    = TimeFlag,                   &
                             RC          = RC                           )
       IF ( RC /= HCO_SUCCESS ) RETURN
   
       ! Data must be ... 
       ! ... concentration ...
       IF ( AreaFlag == 3 .AND. TimeFlag == 0 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.
   
       ELSEIF ( AreaFlag == 3 .AND. TimeFlag == 1 ) THEN
          Lct%Dct%Dta%IsConc = .TRUE.
          FileArr = FileArr * HcoState%TS_EMIS
          MSG = 'Data converted from kg/m3/s to kg/m3: ' // &
                TRIM(Lct%Dct%cName) // ': ' // TRIM(Lct%Dct%Dta%OrigUnit)
          CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
   
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

    ! Data cycle set to 2 (within range) or 3 (exact date): in these 
    ! cases, the preferred date will be equal to the current date, so 
    ! check if the preferred date is indeed within the available range 
    ! (lowDt, uppDt).
    ! For data only to be used within the specified range, set index 
    ! to -1. This will force the scale factors to be set to zero!
    IF ( prefDt < lowDt .OR. prefDt > uppDt ) THEN
       IF ( Lct%Dct%Dta%CycleFlag == 3 ) THEN ! Exact match
          MSG = 'Data is not on exact date: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN 
       ELSEIF ( Lct%Dct%Dta%CycleFlag == 2 ) THEN ! w/in range
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
! !INTERFACE:
!
  SUBROUTINE SrcFile_Parse ( am_I_Root, HcoState, Lct, srcFile, RC )
!
! !USES:
!
    USE HCO_TIDX_MOD,         ONLY : HCO_GetPrefTimeAttr
    USE HCO_CHARTOOLS_MOD,    ONLY : HCO_CharParse
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object
    TYPE(ListCont),   POINTER        :: Lct        ! HEMCO list container
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT)  :: srcFile    ! output string
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! return code
!
! !REVISION HISTORY:
!  01 Oct 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER :: prefYr, prefMt, prefDy, prefHr

    !=================================================================
    ! SrcFile_Parse
    !=================================================================

    ! Initialize to input string
    srcFile = Lct%Dct%Dta%ncFile

    ! Get preferred dates (to be passed to parser
    CALL HCO_GetPrefTimeAttr ( Lct, prefYr, prefMt, prefDy, prefHr, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Call the parser
    CALL HCO_CharParse ( srcFile, prefYr, prefMt, prefDy, prefHr, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

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
END MODULE HCOIO_DataRead_Mod
