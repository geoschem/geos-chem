!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: modis_lai_mod
!
! !DESCRIPTION: Module MODIS\_LAI\_MOD reads the MODIS LAI data at native
!  resolution (either 0.25 x 0.25 or 0.5 x 0.5, in netCDF format) and rebins 
!  it to the proper GEOS-Chem LAI arrays.  This module eliminates the need 
!  for the following GEOS-Chem modules, routines, and data files:
!
! \begin{itemize}
! \item lai\_mod.F
! \item readlai.F
! \item rdlai.F
! \item findmon.F
! \item The \texttt{lai*.global} input files
! \item CMN\_VEL\_mod.F
! \end{itemize}
!
! !INTERFACE: 
!
MODULE Modis_Lai_Mod
!
! !USES:
!
  USE CMN_SIZE_Mod                                ! Size parameters
  USE Directory_Mod                               ! Disk directory paths   
  USE Error_Mod                                   ! Error checking routines
  USE Logical_Mod                                 ! Logical switches
  USE Mapping_Mod                                 ! Mapping weights & areas
  USE Time_Mod                                    ! EXPAND_DATE

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
  INTEGER, PUBLIC              :: DAYS_BTW_MON    ! Days btw LAI midmonths
  REAL*8,  PUBLIC, ALLOCATABLE :: GC_LAI   (:,:)  ! Daily        LAI, G-C grid
  REAL*8,  PUBLIC, ALLOCATABLE :: GC_LAI_PM(:,:)  ! Prev month's LAI, G-C grid
  REAL*8,  PUBLIC, ALLOCATABLE :: GC_LAI_CM(:,:)  ! Curr month's LAI, G-C grid
  REAL*8,  PUBLIC, ALLOCATABLE :: GC_LAI_NM(:,:)  ! Next month's LAI, G-C grid
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Compute_Modis_Lai
  PUBLIC  :: Read_Modis_Lai
  PUBLIC  :: Find_Lai_Month
  PUBLIC  :: Init_Modis_Lai
  PUBLIC  :: Cleanup_Modis_Lai
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: RoundOff
!
! !REMARKS:
!  Functionality of this module:
!  ===========================================================================
!  If you are using the Olson 1992 land map, then this module will pick the
!  MODIS LAI data at 0.5 x 0.5 native resolution.  This is because the legacy
!  code assumed a direct correspondence between the Olson 1992 land map and
!  the MODIS LAI data.  Similarly, if you are using the Olson 2001 land map,
!  then this module will pick the MODIS LAI data at 0.25 x 0.25 resolution.
!                                                                             .
!  Follows the same algorithm as in the IDL codes used to regrid MODIS LAI
!  data (regridmodis_lai_v5.pro; contact GEOS-Chem Support team).
!                                                                             .
!                                                                             .
!  Historical background of how LAI data have been used in GEOS-Chem:
!  ===========================================================================
!  Note that GEOS-Chem (as of April 2012) uses LAI data from two separate
!  sources.  The dry deposition and soil NOx modules rely on the data from 
!  "lai*.global" ASCII files.  These files (which are pre-processed offline 
!  by IDL codes) are generated for each specific GEOS-Chem grid configuration
!  (e.g. 4x5, 2x25, 0.5x0.666 nested grids).  These files are read from disk 
!  by routine RDLAI, which saves the LAI data into the XLAI and XYLAI arrays.
!  XLAI and XYLAI store the leaf area index as a function of Olson land type 
!  (cf Olson 1992 land map).
!                                                                             .
!  However, the MEGAN biogenic emissions code relies on LAI data stored at 
!  1x1 resolution stored in bpch format.  These binary files are read by 
!  routine RDISOLAI (and other underlying routines in lai_mod.F), and are
!  regridded on-the-fly to the current GEOS-Chem grid resolution.
!                                                                             .
!  Therefore, these two sources of LAI data present an inconsistency that 
!  should be resolved.  Also, for the Grid-Indpendent GEOS-Chem project, 
!  we must move away from ASCII files (which prevent interfacing with 
!  external GCMs).  We also cannot assume any particular horizontal grid, 
!  since that is now to be specified at the start of the simulation.
!                                                                             .
!  Also, to facilitate simulations at ultra-fine horizontal resolution, we 
!  will eventually adopt the Olson 2001 land map, which has a native 
!  resolution of 0.25 x 0.25 degrees, and likewise use an updated version 
!  of the MODIS LAI data at 0.25 x 0.25 resolution.
!                                                                             .
!  To resolve these issues, we have created a new module (modis_lai_mod.F90)
!  which reads from the MODIS LAI data in netCDF format at the native 
!  resolution and then regrids the LAI data to GEOS-Chem resolution on-the-
!  fly.  The XLAI array is populated for backwards compatibility with the 
!  existing legacy codes.  The LAI arrays used for MEGAN (ISOLAI, PMISOLAI, 
!  MISOLAI, and NMISOLAI) are now replaced by arrays GC_LAI, GC_LAI_PM, 
!  GC_LAI_CM, and GC_LAI_NM) from modis_lai_mod.F.
!                                                                             .
!  We have validated that the new scheme generates identical XLAI arrays 
!  w/r/t the old scheme.  The arrays GC_LAI etc. differ from the ISOLAI etc. 
!  arrays slightly (but generally agree to within 0.001).  This is due to 
!  the fact that the ISOLAI arrays were regridded from 1 x 1 native
!  resolution, but now we are regridding from much finer resolution 
!  (either 0.5 x 0.5 or 0.25 x 0.25).
!                                                                             .
!  NOTES:
!  (1) At the present time, we have removed all references to the obsolete 
!      XYLAI array and its parent module CMN_VEL_mod.F.
!  (2) At the present time, we have not yet disabled the RDISOLAI function.  
!      We will do so in the future, and will validate this with a separate 
!      benchmark.
!  (3) As of December 2012, XLAI and XLAI2 have been moved out of obsolete
!      module Headers/CMN_DEP_mod.F and are now carried as part of the 
!      Meteorology State object (State_Met).  This modification was made
!      to facilitate the Grid-Independent GEOS-Chem (GIGC) project.
!      
!      -- Bob Yantosca (geos-chem-support@as.harvard.edu), 13 Dec 2012
!                                                                             .
!                                                                             .
!  LAI arrays and where they are (or will be) used in GEOS-Chem:
!  ===========================================================================
!  (1) State_Met%XLAI  --> Used in dry deposition routine DEPVEL
!  (2) State_Met%XLAI2 --> Used to compute XLAI
!  (3) XYLAI           --> %%% OBSOLETE: REMOVED, NOW REPLACED BY XLAI %%%
!  (4) GC_LAI          --> Intended replacement for ISOLAI   (from lai_mod.F)
!  (5) GC_LAI_PM       --> Intended replacement for PMISOLAI (from lai_mod.F)
!  (6) GC_LAI_CM       --> Intended replacement for MISOLAI  (from lai_mod.F)
!  (7) GC_LAI_NM       --> Intended replacement for NMISOLAI (from lai_mod.F)
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
!  05 Apr 2012 - R. Yantosca - Added descriptive comments
!  09 Apr 2012 - R. Yantosca - Fixed error in ROUNDOFF function that caused
!                              numbers to be rounded up incorrectly.
!  09 Apr 2012 - R. Yantosca - Changed variables to REAL*8
!  09 Apr 2012 - R. Yantosca - Now set MODIS_START and MODIS_END depending
!                              on which version of MODIS LAI we are using
!  13 Dec 2012 - R. Yantosca - Remove reference to obsolete CMN_DEP_mod.F;
!                              XLAI, XLAI2 now are carried in State_Met
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Scalars
  INTEGER              :: I_MODIS             ! # of longitudes, MODIS grid
  INTEGER              :: J_MODIS             ! # of latitudes,  MODIS grid
  INTEGER              :: MODIS_START         ! First year of MODIS data  
  INTEGER              :: MODIS_END           ! Last  year of MODIS data
                                              
  ! Arrays                                    
  REAL*4,  ALLOCATABLE :: MODIS_LAI   (:,:)   ! Daily LAI on the MODIS grid
  REAL*4,  ALLOCATABLE :: MODIS_LAI_PM(:,:)   ! MODIS LAI for previous month 
  REAL*4,  ALLOCATABLE :: MODIS_LAI_CM(:,:)   ! MODIS LAI for current month 
  REAL*4,  ALLOCATABLE :: MODIS_LAI_NM(:,:)   ! MODIS LAI for next month 

  ! specify midmonth day for year 2000
  INTEGER, PARAMETER   :: startDay(13) = (/  15,  45,  74, 105,      &
                                            135, 166, 196, 227,      &
                                            258, 288, 319, 349, 380/)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_modis_lai
!
! !DESCRIPTION: Subroutine COMPUTE\_MODIS\_LAI computes the daily MODIS leaf
!  area indices for GEOS-Chem directly from the native grid resolution 
!  (0.25 x 0.25 or 0.5 x 0.5).  The XLAI array (used in the legacy soil NOx 
!  and dry deposition routines) are populated accordingly.  The XYLAI array
!  is now obsolete and has been replaced by XLAI.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Modis_Lai( am_I_Root, State_Met,    doy, mm,  &
                                mapping,   wasModisRead, RC       )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)  :: am_I_Root     ! Are we on the root CPU?
    TYPE(MetState),  INTENT(IN)  :: State_Met     ! Meteorology State object
    INTEGER,         INTENT(IN)  :: doy           ! Day of year
    INTEGER,         INTENT(IN)  :: mm            ! Month for LAI data
    TYPE(MapWeight), POINTER     :: mapping(:,:)  ! "fine" -> "coarse" grid map
    LOGICAL,         INTENT(IN)  :: wasModisRead  ! Was LAI data just read in?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT) :: RC            ! Success or failure?
!
! !REMARKS:
!  Uses same algorithm as RDISOLAI in the existing lai_mod.F.
!
! !REVISION HISTORY: 
!  03 Apr 2012 - R. Yantosca - Initial version
!  05 Apr 2012 - R. Yantosca - Renamed arg "doMonthly" to "wasModisRead"
!  09 Apr 2012 - R. Yantosca - Changed variables to REAL*8
!  09 Apr 2012 - R. Yantosca - Now follows same algorithm as rdlai.F for
!                              populating XLAI array
!  09 Apr 2012 - R. Yantosca - Remove refs to CMN_VEL_mod.F and XYLAI array;
!                              these are now obsolete
!  17 Apr 2012 - R. Yantosca - Now rename "map" object to "mapping" to avoid
!                              name confusion w/ an F90 intrinsic function
!  13 Dec 2012 - R. Yantosca - Add am_I_Root, State_Met, RC arguments
!  13 Dec 2012 - R. Yantosca - XLAI, XLAI2 are now carried in State_Met
!                              instead of in obsolete Headers/CMN_DEP_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE :: FIRST = .TRUE.

    ! Scalars
    INTEGER :: I,     J,    IMUL,    ITD,  IJLOOP
    INTEGER :: C,     II,   JJ,      type, K
    REAL*8  :: mapWt, area, sumArea, DMON, DITD, DIMUL

    ! Arrays
    REAL*8  :: tempArea (0:NVEGTYPE-1)
    REAL*8  :: tempLai  (0:NVEGTYPE-1)
    REAL*8  :: tempLaiCm(0:NVEGTYPE-1)
    REAL*8  :: tempLaiNm(0:NVEGTYPE-1)

    !======================================================================
    ! Interpolate the LAI data on the MODIS grid to current day
    ! Use same algorithm as in routines RDISOLAI (in lai_mod.F)
    !======================================================================
    
    ! Assume success
    RC                = GIGC_SUCCESS

    ! IMUL is days since midmonth
    ! ITD  is days between midmonths
    IF ( doy < startDay(1) ) THEN
       IMUL           = 365 + doy - startDay(12) 
       ITD            = 31
    ELSE
       IMUL           = doy            - startDay(mm)
       ITD            = startDay(mm+1) - startDay(mm)
    ENDIF

    ! Archive the days between midmonths in the LAI data
    DAYS_BTW_MON      = ITD

    ! Cast ITD, IMUL to REAL*8
    DITD              = DBLE( ITD  )
    DIMUL             = DBLE( IMUL )

    ! Fraction of the LAI month that we are in
    DMON              = REAL( IMUL ) / REAL( ITD ) 
       
    ! Interpolate to daily LAI value, on the MODIS grid
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) & 
    !$OMP PRIVATE( I, J   )
    DO J = 1, J_MODIS
    DO I = 1, I_MODIS
       MODIS_LAI(I,J) = MODIS_LAI_CM(I,J)  &
                      + ( ( MODIS_LAI_NM(I,J) - MODIS_LAI_CM(I,J) ) * DMON )
    ENDDO
    ENDDO 
    !$OMP END PARALLEL DO

    !======================================================================
    ! Bin data from the "fine" MODIS grid to the "coarse" GEOS-Chem grid.
    ! Populate arrays for backwards-compatibility w/ existing routines
    ! Use same algorithm as in routine rdlai.F
    !======================================================================
    !$OMP PARALLEL DO                                                 &
    !$OMP DEFAULT( SHARED                                           ) &
    !$OMP PRIVATE( I,         J,       tempArea, tempLai, tempLaiCm ) &
    !$OMP PRIVATE( tempLaiNm, sumArea, IJLOOP,   C,       II        ) & 
    !$OMP PRIVATE( JJ,        type,    area,     K                  )      
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Initialize
       tempArea             = 0d0
       tempLai              = 0d0
       tempLaiCm            = 0d0
       tempLaiNm            = 0d0
       sumArea              = mapping(I,J)%sumarea
       IJLOOP               = ( (J-1) * IIPAR ) + I
       GC_LAI(I,J)          = 0d0

       ! If a new month of MODIS LAI data was just read from disk,
       ! then also initialize the appropriate data arrays here.
       IF ( wasModisRead ) THEN
          GC_LAI_PM(I,J)    = 0d0
          GC_LAI_CM(I,J)    = 0d0
          GC_LAI_NM(I,J)    = 0d0
       ENDIF

       !-------------------------------------------------------------------
       ! Sum up the leaf area indices from all of the the "fine" grid 
       ! boxes (II,JJ) that are located within "coarse" grid box (I,J)
       !-------------------------------------------------------------------
       DO C = 1, mapping(I,J)%count

          ! Extract fields from MAP object
          II                = mapping(I,J)%II(C)
          JJ                = mapping(I,J)%JJ(C)
          type              = mapping(I,J)%olson(C)
          area              = mapping(I,J)%area(C)

          ! Sum of areas corresponding to each Olson
          ! for "coarse" GEOS-Chem grid box (I,J)
          tempArea(type)    = tempArea(type)  + area 

          ! Compute the total leaf area in "coarse" GEOS-Chem 
          ! grid box (I,J) corresponding to each Olson land type
          tempLaiCm(type)   = tempLaiCm(type) + ( MODIS_LAI_CM(II,JJ) * area )
          tempLaiNm(type)   = tempLaiNm(type) + ( MODIS_LAI_NM(II,JJ) * area )

          ! Compute the total leaf area in "coarse" GEOS-Chem
          ! grid box (I,J), irrespective of Olson land type
          GC_LAI(I,J)       = GC_LAI(I,J)     + ( MODIS_LAI(II,JJ)    * area )

          ! If a new month of MODIS LAI data was just read from disk,
          ! then also compute the corresponding total leaf areas in the 
          ! "coarse" GEOS-Chem grid box (I,J).
          IF ( wasModisRead ) THEN
             GC_LAI_PM(I,J) = GC_LAI_PM(I,J)  + ( MODIS_LAI_PM(II,JJ) * area )
             GC_LAI_CM(I,J) = GC_LAI_CM(I,J)  + ( MODIS_LAI_CM(II,JJ) * area )
             GC_LAI_NM(I,J) = GC_LAI_NM(I,J)  + ( MODIS_LAI_NM(II,JJ) * area )
          ENDIF
       ENDDO

       !-------------------------------------------------------------------
       ! Compute the resultant (i.e. for all land types) daily-interpolated 
       ! LAI for the "coarse" GEOS-Chem grid box (I,J).  DAILY_LAI is
       ! a replacement for the ISOLAI array from "lai_mod.F".
       !-------------------------------------------------------------------

       ! Convert leaf area [cm2 leaf] to LAI [cm2 leaf/cm2 grid box]
       ! grid box (I,J), irrespective of Olson land type
       GC_LAI(I,J)       = GC_LAI(I,J)    / sumArea

       ! If a new month of MODIS LAI data was just read from disk,
       ! then also convert the appropriate arrays to leaf area index
       ! [cm2 leaf/cm2 grid box].
       IF ( wasModisRead ) THEN 
          GC_LAI_PM(I,J) = GC_LAI_PM(I,J) / sumArea
          GC_LAI_CM(I,J) = GC_LAI_CM(I,J) / sumArea
          GC_LAI_NM(I,J) = GC_LAI_NM(I,J) / sumArea 
       ENDIF

       !-------------------------------------------------------------------
       ! Compute the LAI for each Olson land type at GEOS-Chem grid box
       ! (I,J).  These will be used to populate the XLAI & XYLAI arrays.
       !-------------------------------------------------------------------
       DO C = 0, NVEGTYPE-1
          
          ! Skip land types that are not in "coarse" grid box (I,J)
          IF ( tempArea(C) > 0d0 ) THEN
          
             ! Ordering for ILAND, IUSE, XLAI, XYLAI etc arrays
             K = mapping(I,J)%ordOlson(C)

             ! Convert leaf area [cm2 leaf] to LAI [cm2 leaf/cm2 grid box]
             tempLaiCm(C) = tempLaiCm(C) / tempArea(C)
             tempLaiNm(C) = tempLaiNm(C) / tempArea(C)

             ! Round off to 1 digit of precision to mimic the fact that
             ! the LAI in the lai*.global files only had one decimal point
             tempLaiCm(C) = RoundOff( tempLaiCm(C), 1 )
             tempLaiNm(C) = RoundOff( tempLaiNm(C), 1 )

             ! This IF statement mimics the algorithm in the obsolete
             ! routine rdlai.F.  We need to keep the same algorithm
             ! for backwards compatibility.
             IF ( FIRST ) THEN 

                !----------------------------------------------------------
                ! %%%%% START OF SIMULATION, FIRST LAI DATA READ %%%%%
                !
                ! Follow original algorithm in the old rdland.F
                ! (1) XLAI  gets read in as the current month of LAI data
                ! (2) XLAI2 gets read in as the next    month of LAI data
                ! (3) XLAI2 is recomputed as the Delta-LAI this month
                !     ([ next month - this month ] / # of days in month)
                ! (4) XLAI is incremented by the amount 
                !      ( Delta-LAI ) * # of days since start of month
                ! 
                ! NOTE: As of Dec 2012, XLAI and XLAI2 are now part of 
                ! the Meteorology State object (bmy, 12/13/12)
                !----------------------------------------------------------
                State_Met%XLAI2(I,J,K) = ( tempLaiNm(C) - tempLaiCm(C)    ) &
                                       / ( DITD                           )
                State_Met%XLAI (I,J,K) = ( tempLaiCm(C)                   ) & 
                                       + ( State_Met%XLAI2(I,J,K) * DIMUL )

             ELSE

                IF ( wasModisRead ) THEN

                   !----------------------------------------------------------
                   ! %%%% SUBSEQUENT DATA READ @ START OF NEW LAI MONTH %%%%
                   !
                   ! Follow original algorithm in the old rdland.F
                   ! (1) XLAI  gets read in as the current month of LAI data
                   ! (2) XLAI2 is computed as the Delta-LAI this month
                   !     (i.e. [next month - this month ] / # of days
                   !----------------------------------------------------------
                   State_Met%XLAI2(I,J,K) = ( tempLaiNm(C) - tempLaiCm(C) ) &
                                          / ( DITD                        )
                   State_Met%XLAI (I,J,K) = ( tempLaiCm(C)                )
 
                ELSE
                
                   !----------------------------------------------------------
                   ! %%%%% ALL OTHER TIMES OF THE MONTH (NO DATA READS %%%%%
                   !
                   ! Follow original algorithm in the old rdland.F
                   ! (1) Increment LAI by the this month's Delta-LAI
                   !----------------------------------------------------------
                   State_Met%XLAI(I,J,K)  = State_Met%XLAI (I,J,K) &
                                          + State_Met%XLAI2(I,J,K)

                ENDIF
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Save
    FIRST = .FALSE.

  END SUBROUTINE Compute_Modis_Lai
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_modis_lai
!
! !DESCRIPTION: Subroutine READ\_MODIS\_LAI reads the MODIS LAI from disk
!  (in netCDF format) for the current month, and for next month.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Read_Modis_Lai( yyyy, mm, wasModisRead )
!
! !USES:
!
    USE m_netcdf_io_open                       ! netCDF file open
    USE m_netcdf_io_read                       ! netCDF read
    USE m_netcdf_io_readattr                   ! netCDF attribute reads
    USE m_netcdf_io_close                      ! netCDF file close
    
#   include "netcdf.inc"                       ! netCDF settings & parameters
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: yyyy               ! Year for LAI data
    INTEGER, INTENT(IN)  :: mm                 ! Month for LAI data
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT) :: wasModisRead       ! Was LAI data just read in?
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
!  05 Apr 2012 - R. Yantosca - Renamed arg "doMonthly" to "wasModisRead"
!  05 Jun 2013 - R. Yantosca - Bug fix, use "mm" for current month index
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: fId                ! netCDF file ID
    INTEGER              :: Pyyyy              ! Previous LAI year
    INTEGER              :: Nyyyy              ! Next     LAI year
    INTEGER              :: Pmm                ! Previous LAI month
    INTEGER              :: Nmm                ! Next     LAI month
    INTEGER              :: yyyymmdd           ! Date variable
                         
    ! Character strings  
    CHARACTER(LEN=255)   :: nc_dir             ! netCDF directory name
    CHARACTER(LEN=255)   :: nc_file            ! netCDF file name
    CHARACTER(LEN=255)   :: nc_path            ! netCDF path name
    CHARACTER(LEN=255)   :: nc_tmpl            ! netCDF file name template
    CHARACTER(LEN=255)   :: v_name             ! netCDF variable name 
    CHARACTER(LEN=255)   :: a_name             ! netCDF attribute name
    CHARACTER(LEN=255)   :: a_val              ! netCDF attribute value
                         
    ! Arrays for netCDF  start and count values
    INTEGER              :: st1d(1), ct1d(1)   ! For 1D arrays    
    INTEGER              :: st3d(3), ct3d(3)   ! For 3D arrays 

    ! SAVEd variables
    INTEGER, SAVE        :: mmLast = -1

    !======================================================================
    ! Test if it is time to read data
    !======================================================================

    ! If we enter a new LAI month, then read MODIS LAI from disk
    ! Otherwise, just exit, since it is not time to read data yet
    IF ( mm /= mmLast ) THEN
       wasModisRead = .TRUE.
       mmLast       = mm
    ELSE
       wasModisRead = .FALSE.
       RETURN
    ENDIF

    ! Save for next iteration

    !======================================================================
    ! Initialize variables
    !======================================================================

    ! Filename template 
    IF ( USE_OLSON_2001 ) THEN
       nc_tmpl = 'For_Olson_2001/MODIS.LAIv.V5.generic.025x025.YYYY.nc'
    ELSE
       nc_tmpl = 'For_Olson_1992/MODIS.LAIv.V5.generic.05x05.YYYY.nc'   
    ENDIF

    ! Construct file path from directory & file name
    nc_dir  = TRIM( DATA_DIR_1x1 ) // 'MODIS_LAI_201204/'

    !======================================================================
    ! Read current month's LAI
    !======================================================================

    ! Test if yyyy is w/in the valid range of MODIS data
    IF ( yyyy >= MODIS_START .and. yyyy <= MODIS_END ) THEN

       ! Here, yyyy lies w/in the MODIS data timespan
       nc_file  = nc_tmpl
       yyyymmdd = yyyy*10000 + mm*100 + 01
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )

    ELSE

       ! Here, yyyy lies outside the MODIS data timespan,
       ! so we have to read data from a different year.
       IF ( USE_OLSON_2001 ) THEN
          IF ( yyyy > MODIS_END ) THEN                    !%%% OLSON 2001 %%%
             yyyymmdd = MODIS_END*10000   + mm*100 + 01   ! Use final year
          ELSE IF ( yyyy < MODIS_START ) THEN             !
             yyyymmdd = MODIS_START*10000 + mm*100 + 01   ! Use 1st year
          ENDIF
       ELSE                                               !%%% OLSON 1992 %%%
          yyyymmdd = 19850001 + mm*100                    ! Use climatology
       ENDIF

       ! Expand date tokens in filename
       nc_file  = nc_tmpl
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )

    ENDIF

    ! Open file for read
    nc_path = TRIM( nc_dir ) // TRIM( nc_file )
    CALL Ncop_Rd( fId, TRIM(nc_path) )
     
    ! Echo info to stdout
    WRITE( 6, 100 ) REPEAT( '%', 79 )
    WRITE( 6, 110 ) TRIM(nc_file)
    WRITE( 6, 120 ) TRIM(nc_dir)

    ! Variable name
    v_name = "MODIS"
    
    ! Read OLSON from file
    st3d   = (/ 1,       1,       mm /)
    ct3d   = (/ I_MODIS, J_MODIS, 1  /)
    CALL NcRd( MODIS_LAI_CM, fId, TRIM(v_name), st3d, ct3d )

    ! Read the OLSON:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val), mm

    ! Close netCDF file
    CALL NcCl( fId )
    
    ! Echo info to stdout
    WRITE( 6, 140 )

    !======================================================================
    ! Read previous month's LAI
    !======================================================================

    ! Previous LAI month
    Pmm = mm - 1

    ! Year corresponding to previous LAI month (readjust for January)
    IF ( Pmm == 0 ) THEN
       Pmm   = 12
       Pyyyy = yyyy - 1
    ELSE
       Pyyyy = yyyy
    ENDIF

    ! Test if Pyyy is w/in the valid range of MODIS data
    IF ( Pyyyy >= MODIS_START .and. Pyyyy <= MODIS_END ) THEN

       ! Here, Pyyyy lies w/in the MODIS data timespan
       nc_file  = nc_tmpl
       yyyymmdd = Pyyyy*10000 + Pmm*100 + 01
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )

    ELSE

       ! Here, yyyy lies outside the MODIS data timespan,
       ! so we have to read data from a different year.
       IF ( USE_OLSON_2001 ) THEN
          IF ( Pyyyy > MODIS_END ) THEN                   !%%% OLSON 2001 %%%
             yyyymmdd = MODIS_END*10000   + Pmm*100 + 01  ! Use final year
          ELSE IF ( Pyyyy < MODIS_START ) THEN            !
             yyyymmdd = MODIS_START*10000 + Pmm*100 + 01  ! Use 1st year
          ENDIF
       ELSE                                               !%%% OLSON 1992 %%%
          yyyymmdd = 19850001 + Pmm*100                   ! Use climatology
       ENDIF

       ! Expand date tokens in filename
       nc_file  = nc_tmpl
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )
          
    ENDIF

    ! Open file for read
    nc_path = TRIM( nc_dir ) // TRIM( nc_file )
    CALL Ncop_Rd( fId, TRIM(nc_path) )
     
    ! Echo info to stdout
    WRITE( 6, 100 ) REPEAT( '%', 79 )
    WRITE( 6, 110 ) TRIM(nc_file)
    WRITE( 6, 120 ) TRIM(nc_dir)

    ! Variable name
    v_name = "MODIS"
    
    ! Read OLSON from file
    st3d   = (/ 1,       1,       Pmm /)
    ct3d   = (/ I_MODIS, J_MODIS, 1   /)
    CALL NcRd( MODIS_LAI_PM, fId, TRIM(v_name), st3d, ct3d )

    ! Read the OLSON:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val), Pmm

    ! Close netCDF file
    CALL NcCl( fId )
    
    ! Echo info to stdout
    WRITE( 6, 140 )

    !======================================================================
    ! Read next month's LAI
    !======================================================================

    ! Next LAI month
    Nmm = mm + 1

    ! Year corresponding to next LAI month (readjust for December)
    IF ( Nmm == 13 ) THEN
       Nmm   = 1
       Nyyyy = yyyy + 1
    ELSE
       Nyyyy = yyyy
    ENDIF

    ! Test if Nyyy is w/in the valid range of MODIS data
    IF ( Nyyyy >= MODIS_START .and. Nyyyy <= MODIS_END ) THEN

       ! Here, Nyyyy lies w/in the MODIS data timespan
       nc_file  = nc_tmpl
       yyyymmdd = Nyyyy*10000 + Nmm*100 + 01
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )

    ELSE

       ! Here, yyyy lies outside the MODIS data timespan,
       ! so we have to read data from a different year.
       IF ( USE_OLSON_2001 ) THEN
          IF ( Nyyyy > MODIS_END ) THEN                   !%%% OLSON 2001 %%%
             yyyymmdd = MODIS_END*10000   + Nmm*100 + 01  ! Use final year
          ELSE IF ( Nyyyy < MODIS_START ) THEN            !
             yyyymmdd = MODIS_START*10000 + Nmm*100 + 01  ! Use 1st year
          ENDIF
       ELSE                                               !%%% OLSON 1992 %%%
          yyyymmdd = 19850001 + Nmm*100                   ! Use climatology
       ENDIF

       ! Expand date tokens in filename
       nc_file  = nc_tmpl
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )
          
    ENDIF

    ! Open file for read
    nc_path = TRIM( nc_dir ) // TRIM( nc_file )
    CALL Ncop_Rd( fId, TRIM(nc_path) )
     
    ! Echo info to stdout
    WRITE( 6, 100 ) REPEAT( '%', 79 )
    WRITE( 6, 110 ) TRIM(nc_file)
    WRITE( 6, 120 ) TRIM(nc_dir)

    ! Variable name
    v_name = "MODIS"
    
    ! Read OLSON from file
    st3d   = (/ 1,       1,       Nmm /)
    ct3d   = (/ I_MODIS, J_MODIS, 1   /)
    CALL NcRd( MODIS_LAI_NM, fId, TRIM(v_name), st3d, ct3d )

    ! Read the OLSON:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val), Nmm

    ! Close netCDF file
    CALL NcCl( fId )
    
    ! Echo info to stdout
    WRITE( 6, 140 )
    WRITE( 6, 100 ) REPEAT( '%', 79 )

    ! FORMAT statements
100 FORMAT( a                                                                 )
110 FORMAT( '%% Opening file  : ',         a                                  )
120 FORMAT( '%%  in directory : ',         a, / , '%%'                        )
130 FORMAT( '%% Successfully read ',       a, ' [', a, '] for month = ', i2.2 )
140 FORMAT( '%% Successfully closed file!'                                    )

  END SUBROUTINE Read_Modis_Lai
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: find_lai_month
!
! !DESCRIPTION: Function FIND\_LAI\_MONTH returns the corresponding LAI 
!  month and year for the current calendar date.  Note that the LAI data
!  starts at mid-month.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Find_Lai_Month( doy, month, year, mm, yyyy )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: doy       ! Current day of year
    INTEGER, INTENT(IN)  :: month     ! Current month
    INTEGER, INTENT(IN)  :: year      ! Current year
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: mm        ! Output month for LAI data
    INTEGER, INTENT(OUT) :: yyyy      ! Output year  for LAI data
!
! !REVISION HISTORY: 
!  05 Jan 1994 - Y. H. Wang, G.M. Gardner, D. Jacob - Initial version
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!  (2 ) Add the current simulation year as input & the current LAI as output.
!       This is necessary for reading in MODIS LAI (mpb,2009).
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  03 Apr 2012 - R. Yantosca - Renamed to FIND_LAI_MONTH; made PUBLIC
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( doy < startDay(1) ) THEN

       ! LAI month & year in December of the preceding calendar year
       mm   = 12
       yyyy = year - 1       

    ELSE IF ( doy < startDay(month) ) THEN

       ! LAI month is the preceding month of this year
       mm   = month - 1
       yyyy = year           

    ELSE

       ! LAI month is the the current month of this year
       mm   = month
       yyyy = year           

    ENDIF

  END SUBROUTINE Find_Lai_Month
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RoundOff
!
! !DESCRIPTION: Rounds a number X to N decimal places of precision.
!\\
!\\
! !INTERFACE:
!
  FUNCTION RoundOff( X, N ) RESULT( Y )
!
! !INPUT PARAMETERS:
! 
    REAL*8,  INTENT(IN) :: X   ! Number to be rounded
    INTEGER, INTENT(IN) :: N   ! Number of decimal places to keep
!
! !RETURN VALUE:
!
    REAL*8              :: Y   ! Number rounded to N decimal places
!
! !REMARKS:
!  The algorithm to round X to N decimal places is as follows:
!  (1) Multiply X by 10**(N+1)
!  (2) If X < 0, then add -5 to X; otherwise add 5 to X
!  (3) Take the integer part of X
!  (4) Divide X by 10**(N+1)
!  (5) Truncate X to N decimal places: INT( X * 10**N ) / 10**N
!                                                                             .
!  Rounding algorithm from: Hultquist, P.F, "Numerical Methods for Engineers 
!   and Computer Scientists", Benjamin/Cummings, Menlo Park CA, 1988, p. 20.
!                                                                             .
!  Truncation algorithm from: http://en.wikipedia.org/wiki/Truncation
!                                                                             .
!  The two algorithms have been merged together for efficiency.
!
! !REVISION HISTORY:
!  06 Apr 2012 - R. Yantosca - Initial version
!  09 Apr 2012 - R. Yantosca - Changed all variables & arguments to REAL*8
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    REAL*8 :: TEN_TO_THE_N                   ! Term for 10**N
    REAL*8 :: TEN_TO_THE_Np1                 ! Term for 10**(N+1)

    ! Pre-compute exponential terms
    TEN_TO_THE_N   = 10d0**N
    TEN_TO_THE_Np1 = 10d0**(N+1)
    
    ! Steps (1) through (4) above
    Y = INT( ( X * TEN_TO_THE_Np1 ) + SIGN( 5d0, X ) ) / TEN_TO_THE_Np1

    ! Step (5) above
    Y = INT( Y * TEN_TO_THE_N ) / TEN_TO_THE_N
  
  END FUNCTION RoundOff
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_modis
!
! !DESCRIPTION: Subroutine INIT\_MODIS\_LAI initializes and allocates
!  all module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Modis_Lai()
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: as

    !======================================================================
    ! Allocate arrays on the "coarse" GEOS-Chem grid
    !======================================================================
    ALLOCATE( GC_LAI( IIPAR, JJPAR ), STAT=as ) 
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'GC_LAI' )
    GC_LAI = 0d0

    ALLOCATE( GC_LAI_PM( IIPAR, JJPAR ), STAT=as ) 
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'GC_LAI_PM' )
    GC_LAI_PM = 0d0

    ALLOCATE( GC_LAI_CM( IIPAR, JJPAR ), STAT=as ) 
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'GC_LAI_CM' )
    GC_LAI_CM = 0d0

    ALLOCATE( GC_LAI_NM( IIPAR, JJPAR ), STAT=as ) 
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'GC_LAI_NM' )
    GC_LAI_NM = 0d0

    !======================================================================
    ! Allocate arrays on the "fine" MODIS grid grid
    !======================================================================
    IF ( USE_OLSON_2001 ) THEN
       I_MODIS     = 1440             ! For Olson 2001, use MODIS LAI
       J_MODIS     = 720              ! on the 0.25 x 0.25 native grid
       MODIS_START = 2005             ! First year of MODIS data  
       MODIS_END   = 2009             ! Last  year of MODIS data
    ELSE
       I_MODIS     = 720              ! For Olson 1992, use MODIS LAI
       J_MODIS     = 360              ! on the 0.5 x 0.5 native grid
       MODIS_START = 2000             ! First year of MODIS data  
       MODIS_END   = 2008             ! Last  year of MODIS data
    ENDIF

    ALLOCATE( MODIS_LAI( I_MODIS, J_MODIS ), STAT=as ) 
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'MODIS_LAI' )
    MODIS_LAI = 0d0

    ALLOCATE( MODIS_LAI_PM( I_MODIS, J_MODIS ), STAT=as ) 
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'MODIS_LAI_PM' )
    MODIS_LAI_PM = 0d0

    ALLOCATE( MODIS_LAI_CM( I_MODIS, J_MODIS ), STAT=as ) 
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'MODIS_LAI_CM' )
    MODIS_LAI_CM = 0d0

    ALLOCATE( MODIS_LAI_NM( I_MODIS, J_MODIS ), STAT=as ) 
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'MODIS_LAI_NM' )
    MODIS_LAI_NM = 0d0

  END SUBROUTINE Init_Modis_Lai
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_modis_lai
!
! !DESCRIPTION: Subroutine CLEANUP\_MODIS\_LAI deallocates all
!  previously-allocated module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Modis_Lai
!
! !REVISION HISTORY:'
!  03 Apr 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( GC_LAI       ) ) DEALLOCATE( GC_LAI       )
    IF ( ALLOCATED( GC_LAI_PM    ) ) DEALLOCATE( GC_LAI_PM    )
    IF ( ALLOCATED( GC_LAI_CM    ) ) DEALLOCATE( GC_LAI_CM    )
    IF ( ALLOCATED( GC_LAI_NM    ) ) DEALLOCATE( GC_LAI_NM    )
    IF ( ALLOCATED( MODIS_LAI    ) ) DEALLOCATE( MODIS_LAI    )
    IF ( ALLOCATED( MODIS_LAI_PM ) ) DEALLOCATE( MODIS_LAI_PM )
    IF ( ALLOCATED( MODIS_LAI_CM ) ) DEALLOCATE( MODIS_LAI_CM )
    IF ( ALLOCATED( MODIS_LAI_NM ) ) DEALLOCATE( MODIS_LAI_NM )

  END SUBROUTINE Cleanup_Modis_Lai
!EOC
END MODULE Modis_Lai_Mod
