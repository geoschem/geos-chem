!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: modis_lai_mod
!
! !DESCRIPTION: Module MODIS\_LAI\_MOD reads the MODIS LAI data at native
!  resolution (either 0.25 x 0.25 or 0.5 x 0.5) and rebins it to the proper
!  GEOS-Chem leaf area index arrays.  This module eliminates the need for the
!  following GEOS-Chem modules and routines:
!
! \begin{itemize}
! \item lai_mod.F
! \item readlai.F
! \item rdlai.F
! \item The \texttt{lai*.global} input files
! \end{itemize}
!\\
!\\
! !INTERFACE: 
!
MODULE Modis_Lai_Mod
!
! !USES:
!
  USE CMN_SIZE_Mod                                ! Size parameters
  USE CMN_DEP_Mod                                 ! IREG, ILAND, IUSE, FRCLND
  USE CMN_VEL_Mod                                 ! IJREG, IJLAND, IJUSE
  USE Directory_Mod                               ! Disk directory paths   
  USE Error_Mod                                   ! Error checking routines
  USE Logical_Mod                                 ! Logical switches
  USE Mapping_Mod                                 ! Mapping weights
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
!
! !REMARKS:
!  If you are using the Olson 1992 land map, then this module will pick the
!  MODIS LAI data at 0.5 x 0.5 native resolution.  This is because the legacy
!  code assumed a direct correspondence between the Olson 1992 land map and
!  the MODIS LAI data.  Similarly, if you are using the Olson 2001 land map,
!  then thismodule will pick the MODIS LAI data at 0.25 x 0.25 resolution.
!                                                                             .
!  Follows the same algorithm as in the IDL codes used to regrid MODIS LAI
!  data (regridmodis_lai_v5.pro; contact GEOS-Chem Support team).
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Scalars
  INTEGER              :: I_MODIS             ! # of longitudes, MODIS grid
  INTEGER              :: J_MODIS             ! # of latitudes,  MODIS grid
                                              
  ! Arrays                                    
  REAL*8,  ALLOCATABLE :: MODIS_LAI   (:,:)   ! Daily LAI on the MODIS grid
  REAL*8,  ALLOCATABLE :: MODIS_LAI_PM(:,:)   ! MODIS LAI for previous month 
  REAL*8,  ALLOCATABLE :: MODIS_LAI_CM(:,:)   ! MODIS LAI for current month 
  REAL*8,  ALLOCATABLE :: MODIS_LAI_NM(:,:)   ! MODIS LAI for next month 
!                                             
! !DEFINED PARAMETERS:                        
!                                             
  INTEGER, PARAMETER   :: MODIS_START = 2000  ! First year of MODIS data  
  INTEGER, PARAMETER   :: MODIS_END   = 2008  ! Last  year of MODIS data

  ! specify midmonth day for year 2000
  INTEGER, PARAMETER   :: startDay(13) = (/  15,  45,  74, 105,      &
                                            135, 166, 196, 227,      &
                                            258, 288, 319, 349, 380/)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_modis_lai
!
! !DESCRIPTION: Subroutine COMPUTE\_MODIS\_LAI computes the daily MODIS leaf
!  area indices for GEOS-Chem directly from the native grid resolution 
!  (0.25 x 0.25 or 0.5 x 0.5).  Arrays XLAI and XYLAI (used in the legacy
!  soil NOx and dry deposition routines) are populated accordingly.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Modis_Lai( doy, mm, map, doMonthly )
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN) :: doy         ! Day of year
    INTEGER,         INTENT(IN) :: mm          ! Month for LAI data
    TYPE(MapWeight), POINTER    :: map(:,:)    ! "fine" -> "coarse" grid map
    LOGICAL,         INTENT(IN) :: doMonthly   ! Interpolate monthly files? 
!
! !REMARKS:
!  Uses same algorithm as RDISOLAI in the existing lai_mod.F.
!
! !REVISION HISTORY: 
!  03 Apr 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: I,      J,       IMUL,   ITD
    INTEGER :: IJLOOP
    INTEGER :: C,      II,      JJ,     type
    REAL*4  :: mapWt,  area,    sumArea
    REAL*8  :: FRAC,   leafArea
    
    ! Arrays
    REAL*8  :: tempArea(0:NVEGTYPE-1)
    REAL*8  :: tempLai (0:NVEGTYPE-1)

    !======================================================================
    ! Interpolate the LAI data on the MODIS grid to current day
    !======================================================================

    ! IMUL is days since midmonth
    ! ITD  is days betw een midmonths
    IF ( doy < startDay(1) ) THEN
       IMUL           = 365 + doy - startDay(12) 
       ITD            = 31
    ELSE
       IMUL           = doy            - startDay(mm)
       ITD            = startDay(mm+1) - startDay(mm)
    ENDIF

    ! Archive the days between midmonths in the LAI data
    DAYS_BTW_MON      = ITD

    ! Fraction of the LAI month that we are in
    FRAC              = DBLE( IMUL ) / DBLE( ITD ) 
       
    ! Interpolate to daily LAI value, on the MODIS grid
    !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J )
    DO J = 1, J_MODIS
    DO I = 1, I_MODIS
       MODIS_LAI(I,J) = MODIS_LAI_CM(I,J)  &
                      + ( FRAC * ( MODIS_LAI_NM(I,J) - MODIS_LAI_CM(I,J) ) )
    ENDDO
    ENDDO 
    !$OMP END PARALLEL DO

    !======================================================================
    ! Bin data from the "fine" MODIS grid to the "coarse" GEOS-Chem grid.
    ! Populate arrays for backwards-compatibility w/ existing codes.
    !======================================================================
    !$OMP PARALLEL DO                                             &
    !$OMP DEFAULT( SHARED                                       ) &
    !$OMP PRIVATE( I, J,  tempArea, tempLai, sumArea, IJLOOP   ) &
    !$OMP PRIVATE( C, II, JJ,       type,    area,     leafArea ) 
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Initialize
       tempArea             = 0d0
       tempLai              = 0d0
       sumArea              = map(I,J)%sumarea
       IJLOOP               = ( (J-1) * IIPAR ) + I
       GC_LAI   (I,J)       = 0d0

       IF ( doMonthly ) THEN
          GC_LAI_PM(I,J)    = 0d0
          GC_LAI_CM(I,J)    = 0d0
          GC_LAI_NM(I,J)    = 0d0
       ENDIF

       !-------------------------------------------------------------------
       ! Sum up the leaf area indices from all of the the "fine" grid 
       ! boxes (II,JJ) that are located within "coarse" grid box (I,J)
       !-------------------------------------------------------------------
       DO C = 1, map(I,J)%count

          ! Extract fields from MAP object
          II                = map(I,J)%II(C)
          JJ                = map(I,J)%JJ(C)
          type              = map(I,J)%olson(C)
          area              = map(I,J)%area(C)

          ! Sum of areas corresponding to each Olson
          ! for "coarse" GEOS-Chem grid box (I,J)
          tempArea(type)    = tempArea(type) + area 

          ! Compute the total leaf area in "coarse" GEOS-Chem 
          ! grid box (I,J) corresponding to each Olson land type
          tempLai(type)     = tempLai(type)  + ( MODIS_LAI(II,JJ)    * area )

          ! Compute the total leaf area in "coarse" GEOS-Chem
          ! grid box (I,J), irrespective of Olson land type
          GC_LAI(I,J)       = GC_LAI(I,J)    + ( MODIS_LAI(II,JJ)    * area )

          IF ( doMonthly ) THEN
             GC_LAI_PM(I,J) = GC_LAI_PM(I,J) + ( MODIS_LAI_PM(II,JJ) * area )
             GC_LAI_CM(I,J) = GC_LAI_CM(I,J) + ( MODIS_LAI_CM(II,JJ) * area )
             GC_LAI_NM(I,J) = GC_LAI_NM(I,J) + ( MODIS_LAI_NM(II,JJ) * area )
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

       IF ( doMonthly ) THEN 
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
          
             ! Convert leaf area [cm2 leaf] to LAI [cm2 leaf/cm2 grid box]
             tempLai(C) = tempLai(C) / tempArea(C)

             ! Save into GEOS-Chem arrays for backwards compatibility
             XLAI ( I, J,   map(I,J)%ordOlson(C) ) = tempLai(C)
             XYLAI( IJLOOP, map(I,J)%ordOlson(C) ) = tempLai(C)

          ENDIF
       ENDDO

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Compute_Modis_Lai
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Read_Modis_Lai( yyyy, mm, doMonthly )
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
    LOGICAL, INTENT(OUT) :: doMonthly
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
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

    ! If we enter a new LAI month, then read the data below
    ! Otherwise, just exit 
    IF ( mm /= mmLast ) THEN
       doMonthly = .TRUE.
       mmLast    = mm
    ELSE
       doMonthly = .FALSE.
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

    ! Pick either the proper year of MODIS data (2000-2008),
    ! or if outside that range, the climatology (1985)
    IF ( yyyy >= MODIS_START .and. yyyy <= MODIS_END ) THEN
       nc_file  = nc_tmpl
       yyyymmdd = yyyy*10000 + mm*100 + 01
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )
    ELSE
       nc_file  = nc_tmpl
       yyyymmdd = 19850001 + mm*100
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

    ! Previous month
    Pmm = mm - 1

    ! Previous year (readjust accordingly if this month is January)
    IF ( Pmm == 0 ) THEN
       Pmm   = 12
       Pyyyy = yyyy - 1
    ELSE
       Pyyyy = yyyy
    ENDIF

    ! Pick either the proper year of MODIS data (2000-2008),
    ! or if outside that range, the climatology (1985)
    IF ( Pyyyy >= MODIS_START .and. Pyyyy <= MODIS_END ) THEN
       nc_file  = nc_tmpl
       yyyymmdd = Pyyyy*10000 + Pmm*100 + 01
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )
    ELSE
       nc_file  = nc_tmpl
       yyyymmdd = 19850001 + Pmm*100
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

    ! Previous month
    Nmm = mm + 1

    ! Previous year (readjust accordingly if this month is December)
    IF ( Nmm == 13 ) THEN
       Nmm   = 1
       Nyyyy = yyyy + 1
    ELSE
       Nyyyy = yyyy
    ENDIF

    ! Pick either the proper year of MODIS data (2000-2008),
    ! or if outside that range, the climatology (1985)
    IF ( Nyyyy >= MODIS_START .and. Nyyyy <= MODIS_END ) THEN
       nc_file  = nc_tmpl
       yyyymmdd = Nyyyy*10000 + Nmm*100 + 01
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )
    ELSE
       nc_file  = nc_tmpl
       yyyymmdd = 19850001 + Nmm*100
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!  03 Apr 2012 - R. Yantosca - Copied to this module as a PRIVATE routine
!  03 Apr 2012 - R. Yantosca - Renamed to FIND_LAI_MONTH
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_modis
!
! !DESCRIPTION: Subroutine INIT\_OLSON\_LANDMAP reads Olson land map 
! information from disk (in netCDF format).
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
       I_MODIS = 1440             ! For Olson 2001, use MODIS LAI
       J_MODIS = 720              ! on the 0.25 x 0.25 native grid
    ELSE
       I_MODIS = 720              ! For Olson 1992, use MODIS LAI
       J_MODIS = 360              ! on the 0.5 x 0.5 native grid
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
    IF ( ALLOCATED( MODIS_LAI_CM ) ) DEALLOCATE( MODIS_LAI_CM )
    IF ( ALLOCATED( MODIS_LAI_NM ) ) DEALLOCATE( MODIS_LAI_NM )

  END SUBROUTINE Cleanup_Modis_Lai
!EOC
END MODULE Modis_Lai_Mod
