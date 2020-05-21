!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_dataread_mod.F90
!
! !DESCRIPTION: Module HCOIO\_DataRead\_Mod controls data processing
! for HEMCO.
! Depending on the model environment (standard, ESMF, etc.), it invokes
! the corresponding routines to read the data from file, convert units
! as required, and interpolate the data onto the model grid.
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
!
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
  USE HCO_Types_Mod
  USE HCO_Error_Mod
  USE HCO_CharTools_Mod
  USE HCO_State_Mod,       ONLY : Hco_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOIO_DataRead
!
! !REVISION HISTORY:
!  22 Aug 2013 - C. Keller   - Initial version
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  22 Feb 2016 - C. Keller   - Environment specific routines are now
!                              in respective modules.
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_DataRead
!
! !DESCRIPTION: Routine HCOIO\_DataRead invokes the appropriate data reading
! routines for the given model environment.
!\\
!\\
! !INTERFACE:
  !
  SUBROUTINE HCOIO_DataRead( HcoState, Lct, RC )
!
! !USES:
!
#if defined(ESMF_)
    USE HCOIO_READ_ESMF_MOD,  ONLY : HCOIO_READ_ESMF
#else
    USE HCOIO_READ_STD_MOD,   ONLY : HCOIO_READ_STD
#endif
!
! !INPUT PARAMETERS:
!
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
!  22 Feb 2016 - C. Keller   - Now calls down to model-specific routines.
!  24 Mar 2016 - C. Keller   - Removed LUN and CloseFile. Not needed any more.
!  05 Oct 2018 - R. Yantosca - Read files with cycle flag "E" just once.
!                              Also improve the error trapping.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=512) :: MSG

    ! SAVEd scalars
    LOGICAL, SAVE      :: doPrintWarning = .TRUE.

    !=======================================================================
    ! HCOIO_DATAREAD begins here
    !=======================================================================

    ! Assume success until proven otherwise
    RC = HCO_SUCCESS

    ! Denote we are entering this routine
    CALL HCO_ENTER( HcoState%Config%Err,                                     &
                    'HCOIO_DataRead (hcoio_dataread_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

#if defined( ESMF_ ) || defined( MODEL_GEOS )

    !-----------------------------------------------------------------------
    ! %%%%% ESMF environment (e.g. GCHP, GEOS-DAS) %%%%%%
    !-----------------------------------------------------------------------

    ! Call ESMF I/O routine.
    CALL HCOIO_READ_ESMF( HcoState, Lct, RC )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Error encountered in routine HCOIO_Read_ESMF!'
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

#else

    !-----------------------------------------------------------------------
    ! %%%%% Standard environment (e.g. GEOS-Chem "Classic") %%%%%
    !-----------------------------------------------------------------------

    ! If the file has cycle flag "E" (e.g. it's a restart file), then we will
    ! read it only once and then never again.  If the file has already been
    ! read on a previous call, then don't call HCOIO_READ_STD. (bmy, 10/4/18)
    IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT .and.                      &
         Lct%Dct%Dta%UpdtFlag  == HCO_UFLAG_ONCE  .and.                      &
         Lct%Dct%Dta%isTouched                          ) THEN

       ! Print a warning message only once
       IF ( doPrintWarning ) THEN
          doPrintWarning = .FALSE.
          MSG = 'No further attempts will be made to read file: ' //         &
                TRIM( Lct%Dct%Dta%NcFile )
          CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, WARNLEV=1 )
       ENDIF

       ! Return without reading
       CALL HCO_LEAVE( HcoState%Config%Err, RC )
       RETURN
    ENDIF

    ! Call standard I/O routines
    CALL HCOIO_READ_STD( HcoState, Lct, RC )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Error encountered in routine HCOIO_Read_Std!'
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

#endif

    ! Denote we are leaving this routine
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE HCOIO_DataRead
!EOC
END MODULE HCOIO_DataRead_Mod
