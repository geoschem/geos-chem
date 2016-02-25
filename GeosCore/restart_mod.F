!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: restart_mod
!
! !DESCRIPTION: Module RESTART\_MOD contains variables and routines which 
!  are used to read and write restart files for GEOS-Chem advected tracers
!  in units of [v/v] mixing ratio, and chemical species in concentration 
!  units of [molec/cm3]. [v/v] represents moles tracer per moles dry air.
!\\
!\\
! !INTERFACE:
!
      MODULE RESTART_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD
      USE GIGC_ErrCode_Mod
      USE inquireMod, ONLY : findFreeLUN
      USE Precision_Mod

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS
!
      PUBLIC  :: SET_RESTART           ! called from input_mod
      PUBLIC  :: READ_GC_RESTART
      PUBLIC  :: WRITE_GC_RESTART
      PUBLIC  :: READ_CSPEC_FILE       ! For binary-punch read
      PUBLIC  :: INIT_GC_RESTART_NC    ! For NetCDF write
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CHECK_DIMENSIONS
      PRIVATE :: CONVERT_TRACER_TO_VV
      PRIVATE :: READ_GC_RESTART_BPCH  ! For binary-punch read
      PRIVATE :: COPY_STT              ! For binary-punch read
      PRIVATE :: CHECK_DATA_BLOCKS     ! For binary-punch read
      PRIVATE :: WRITE_GC_RESTART_BPCH ! For binary-punch write
      PRIVATE :: MAKE_CSPEC_FILE       ! For binary-punch write
      PRIVATE :: READ_GC_RESTART_NC    ! For NetCDF read
      PRIVATE :: WRITE_GC_RESTART_NC   ! For NetCDF write
!
! !DEFINED PARAMETERS:
!
      ! Logical indicator of whether species data is read from NetCDF file 
      ! (true) or if default values (previously from globchem.dat) are 
      ! used (false). This parameter is only used if using NetCDF files.
      INTEGER, SAVE, PUBLIC :: SPC_IN_NC_RST  

      ! Full path name of the advected tracer restart file (INPUT)
      CHARACTER(LEN=255)    :: INPUT_RESTART_FILE  

      ! Full path name of (w/ replaceable tokens) of the 
      ! advected tracer restart file (OUTPUT)
      CHARACTER(LEN=255)    :: OUTPUT_RESTART_FILE 
!
! !REVISION HISTORY:
!  25 Jun 2002 - R. Yantosca - Initial version
!  (1 ) Moved routines "make_restart_file.f"" and "read_restart_file.f" into
!        this module.  Also now internal routines to "read_restart_file.f"
!        are now a part of this module.  Now reference "file_mod.f" to get
!        file unit numbers and error checking routines. (bmy, 6/25/02)
!  (2 ) Now reference AD from "dao_mod.f".  Now reference "error_mod.f".
!        Also added minor bug fix for ALPHA platform. (bmy, 10/15/02)
!  (3 ) Now references "grid_mod.f" and the new "time_mod.f" (bmy, 2/11/03)
!  (4 ) Added error-check and cosmetic changes (bmy, 4/29/03)
!  (5 ) Removed call to COPY_STT_FOR_OX, it's obsolete (bmy, 8/18/03)
!  (6 ) Add fancy output (bmy, 4/26/04)
!  (7 ) Added routine SET_RESTART.  Now reference "logical_mod.f" and
!        "tracer_mod.f" (bmy, 7/20/04)
!  (8 ) Removed obsolete routines TRUE_TRACER_INDEX and COPY_DATA_FOR_CO_OH
!        (bmy, 6/28/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now pass TAU via the arg list in MAKE_RESTART_FILE (bmy, 12/15/05)
!  (11) Add MAKE_CSPEC_FILE and READ_CSPEC_FILE routines to save and read
!        CSPEC_FULL restart files (dkh, 02/12/09)
!  11 Jul 2011 - R. Yantosca - Corrected mis-indexing problem w/ the 
!                              CSPEC restart file
!  21 Jul 2011 - M. Long     - Now include F77_CMN_SIZE instead of CMN_SIZE
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  19 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  09 Feb 2016 - E. Lundgren - Add NetCDF restart file routines
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
! !IROUTINE: init_gc_restart_nc
!
! !DESCRIPTION: Subroutine INIT\_GC\_RESTART\_NC initializes the GEOS-Chem
!  restart file collection and populates it with containers. This collection
!  includes restart data that traditionally would be found in three different 
!  binary-punch format files: tracers, species, and ocean mercury.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GC_RESTART_NC( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE COMODE_LOOP_MOD,    ONLY : NTSPEC, IFAM, NAMEGAS
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : AREA_M2
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : ID_Hg0, ID_Hg2, ID_HgP
      USE TRACERID_MOD,       ONLY : ID_Hg_tot, N_Hg_CATS

      USE HCO_Error_Mod
      USE HCO_Diagn_Mod
!
! !INPUT PARAMETERS: 
!
      LOGICAL,          INTENT(IN)    :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput),   INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  09 Feb 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local Variables  
      INTEGER            :: GC_RST_Collection, N, M, cID
      INTEGER            :: DeltaYMD, DeltaHMS 
      REAL(sp)           :: TS
      REAL(fp), POINTER  :: GridAreas_ptr(:,:) => NULL()
      CHARACTER(LEN=60)  :: ContainerName
      CHARACTER(LEN=60)  :: Prefix, Suffix, TracerIDStr
      CHARACTER(LEN=255) :: MSG, LOC = 
     &                      'INIT_GC_RESTART_NC (restart_mod.F90)'

      !=================================================================
      ! INIT_GC_RESTART_NC begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Define collection variables
      GridAreas_ptr    => AREA_M2(:,:,1)
      TS     =  GET_TS_CHEM() * 60.0_sp

      !------------------------------------------------------------------
      ! Create restart collection for GEOS-Chem.  This will keep the
      ! GEOS-Chem restart separate from other output.
      !------------------------------------------------------------------
#if defined ( ESMF_ )
      deltaYMD = 0 
      deltaHMS = 1
#else
      deltaYMD = 99999999
      deltaHMS = 999999
#endif     
      ! Need to adjust this for restart file
      CALL DiagnCollection_Create( am_I_Root,                      
     &                              NX           = IIPAR,           
     &                              NY           = JJPAR,           
     &                              NZ           = LLPAR,           
     &                              TS           = TS,              
     &                              AM2          = GridAreas_ptr,       
     &                              COL          = GC_RST_Collection,    
     &                              PREFIX       = 'GEOSChem_restart',     
     &                              deltaYMD     = deltaYMD,       
     &                              deltaHMS     = deltaHMS,       
     &                              OutTimeStamp = HcoDiagnEnd,     
     &                              RC           = RC                )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP( 'Error creating GEOS-Chem restart collection',
     &                    LOC ) 
      ENDIF
      
      ! Cleanup
      GridAreas_ptr => NULL()
      
      ! Save GEOS-Chem restart collection ID in Input_Opt
      Input_Opt%GC_RST_COLLECTION = GC_RST_Collection

      !-----------------------------------------------------------------
      ! Create restart containers for tracers
      !-----------------------------------------------------------------
      
      ! Initialize container id
      cID = 0

      ! Loop over tracers
      DO N = 1, Input_Opt%N_TRACERS
         
         ! Define container name and id
         ContainerName = 'TRC_' // TRIM( Input_Opt%TRACER_NAME(N) )
         cID = cID + 1

         ! Create container
         CALL Diagn_Create( am_I_Root,                     
     &                      Col       = GC_RST_Collection,         
     &                      cId       = cId,                 
     &                      cName     = TRIM( ContainerName ), 
     &                      SpaceDim  =  3,                
     &                      OutUnit   = 'mol mol-1',           
     &                      OutOper   = 'Instantaneous', 
     &                      RC        = RC )
         
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot create restart tracer container: ' 
     &            // TRIM(ContainerName)
            CALL ERROR_STOP( MSG, LOC ) 
         ENDIF  
      ENDDO

      !-----------------------------------------------------------------
      ! Create restart containers for species
      !-----------------------------------------------------------------
      ! Only save species if chemistry simulation
      IF ( Input_Opt%LCHEM .AND. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN 

         ! Loop over the total # of species.  This also includes the "fake"
         ! prod/loss family species for the ND65 diagnostic.
         DO N = 1 , NTSPEC(1)
         
            ! Ignore ND65 handling pending further netcdf diagnostics work.
            ! If the species is included in the IFAM array, which denotes 
            ! the index # of species in the CSPEC array that are "fake" ND65 
            ! prod/loss families, then skip. 
            !IF ( ANY( IFAM == N ) ) THEN
            !   CYCLE
            !ENDIF
         
            ! Define container name and ID
            ContainerName = 'SPC_' // TRIM( NAMEGAS(N) )
            cID = cID + 1
         
            ! Create container
            CALL Diagn_Create( am_I_Root,                     
     &                         Col       = GC_RST_Collection,         
     &                         cId       = cID,                 
     &                         cName     = TRIM( ContainerName ), 
     &                         SpaceDim  =  3,                
     &                         OutUnit   = 'mol/mol',           
     &                         OutOper   = 'Instantaneous', 
     &                         RC        = RC )
            
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot create restart species container: ' 
     &               // TRIM(ContainerName)
               CALL ERROR_STOP( MSG, LOC ) 
            ENDIF  
         ENDDO
      ENDIF

      !-----------------------------------------------------------------
      ! If Hg simulation, create restart containers for mercury
      !-----------------------------------------------------------------
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         !--------------------------------------------------------------
         ! Total Hg(0), Hg(II), and Hg(P) in ocean
         !--------------------------------------------------------------
         DO M = 1, 3
         
            SELECT CASE( M )
               CASE ( 1 )
                  Prefix = '_Hg0_' 
               CASE ( 2 )
                  Prefix = '_Hg2_'
               CASE ( 3 )
                  Prefix = '_HgP_' 
            END SELECT
         
            ! Define container name and ID
            ContainerName = 'OCEAN_HG' // TRIM( Prefix ) // 'Total' 
            cID = cID + 1
            
            ! Create container
            CALL Diagn_Create( am_I_Root,                     
     &                         Col       = GC_RST_Collection,         
     &                         cId       = cID,                 
     &                         cName     = TRIM( ContainerName ), 
     &                         SpaceDim  =  2,                
     &                         OutUnit   = 'kg',           
     &                         OutOper   = 'Instantaneous', 
     &                         RC        = RC )
            
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot create restart species container: ' 
     &               // TRIM(ContainerName)
               CALL ERROR_STOP( MSG, LOC ) 
            ENDIF  
         ENDDO

         !-----------------------------------------------------------
         ! Hg snowpack on land and ocean
         !-----------------------------------------------------------
         DO M = 1, 4
            SELECT CASE( M )
               CASE ( 1 )
                  Prefix = 'SNOW_HG_OCN_Hg0_'  ! Reducible on ocean
               CASE ( 2 )
                  Prefix = 'SNOW_HGN_OCN_Hg0_' ! Non-reducible on ocean
               CASE ( 3 )
                  Prefix = 'SNOW_HG_LND_Hg0_'  ! Reducible on land
               CASE ( 4 )
                  Prefix = 'SNOW_HGN_LND_Hg0_' ! Non-reducible on land
            END SELECT
         
            ! Loop over all categories to include tagged Hg
            DO N = 1, N_Hg_CATS
            
               ! Define container name and id
               ! (temperarily append tracer # until name array exists)
               write( TracerIDStr, '(I5)' ) N
               ContainerName = TRIM( Prefix ) // TRIM( TracerIDStr )
               cID = cId + 1
            
               ! Create container
               CALL Diagn_Create( am_I_Root,                     
     &                            Col       = GC_RST_Collection,         
     &                            cId       = cID,                 
     &                            cName     = TRIM( ContainerName ), 
     &                            SpaceDim  =  2,                
     &                            OutUnit   = 'kg',           
     &                            OutOper   = 'Instantaneous', 
     &                            RC        = RC )
               
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot create restart species container: ' 
     &                  // TRIM(ContainerName)
                  CALL ERROR_STOP( MSG, LOC ) 
               ENDIF  
            ENDDO
         ENDDO

         !--------------------------------------------------------------
         ! Additional restart containers for tagged Hg simulation 
         !--------------------------------------------------------------
         IF ( Input_Opt%LSPLIT ) THEN

            !-----------------------------------------------------------
            ! Tagged Hg(0), Hg(II), and Hg(P) in ocean
            !-----------------------------------------------------------
            DO M = 1, 3
               SELECT CASE( M )
                  CASE ( 1 )
                     Prefix = '_Hg0_' 
                  CASE ( 2 )
                     Prefix = '_Hg2_'
                  CASE ( 3 )
                     Prefix = '_HgP_' 
               END SELECT
            
               DO N = 2, N_Hg_CATS
               
                  ! Define container name and ID
                  ! (temperarily append tracer # until name array exists)
                  write( TracerIDStr, '(I5)' ) N
                  Suffix = TRIM( TracerIDStr )
                  ContainerName = 'OCEAN_HG' // TRIM( Prefix ) 
     &                            // TRIM( Suffix )
                  cID = cId + 1
               
                  ! Create container
                  CALL Diagn_Create( am_I_Root,                     
     &                               Col       = GC_RST_Collection,         
     &                               cId       = cID,                 
     &                               cName     = TRIM( ContainerName ), 
     &                               SpaceDim  =  2,                
     &                               OutUnit   = 'kg',           
     &                               OutOper   = 'Instantaneous', 
     &                               RC        = RC )
                  
                  IF ( RC /= HCO_SUCCESS ) THEN
                     MSG = 'Cannot create restart species container: ' 
     &                     // TRIM(ContainerName)
                     CALL ERROR_STOP( MSG, LOC ) 
                  ENDIF  
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      END SUBROUTINE INIT_GC_RESTART_NC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_restart
!
! !DESCRIPTION: Subroutine SET\_RESTART initializes the variables 
!  INPUT\_RESTART\_FILE and OUTPUT\_RESTART\_FILE with the values read from 
!  the \texttt{input.geos} file.  These specify the names of the input and
!  output restart files for GEOS-Chem advected tracers.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_RESTART( INFILE, OUTFILE )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=255) :: INFILE    ! Advected tracer input  restart file
      CHARACTER(LEN=255) :: OUTFILE   ! Advected tracer output restart file 
! 
! !REVISION HISTORY: 
!  09 Jul 2004 - R. Yantosca - Initial version
!  11 Jul 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! SET_RESTART begins here
      !=================================================================
      INPUT_RESTART_FILE  = INFILE
      OUTPUT_RESTART_FILE = OUTFILE
     
      END SUBROUTINE SET_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gc_restart
!
! !DESCRIPTION: Subroutine READ\_GC\_RESTART is wrapper subroutine
!  to read restart files in either binary-punch (bpch) or NetCDF format.
!  If reading bpch, indicated by setting compile option BPCH=y, then only 
!  tracer concentrations are read in by this routine. If reading NetCDF,
!  indicated by setting compile option NETCDF=y, then both tracers and 
!  species are read in using this routine. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_GC_RESTART( am_I_Root, YYYYMMDD,  HHMMSS,
     &                            Input_Opt, State_Met, State_Chm, 
     &                            RC ) 
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE OCEAN_MERCURY_MOD,  ONLY : READ_OCEAN_Hg_RESTART
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root  ! Are we on the root CPU?
      INTEGER,        INTENT(IN)    :: YYYYMMDD   ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)    :: HHMMSS     ! hh:mm:ss   GMT time
      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC         ! Failure or success
! 
! !REVISION HISTORY: 
!
!  09 Feb 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! READ_GC_RESTART begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

#if defined( BPCH_RST_IN )
      ! Read binary-punch restart files
      CALL READ_GC_RESTART_BPCH( YYYYMMDD,  HHMMSS,   Input_Opt, 
     &                           State_Met, State_Chm )

      ! Read ocean Hg initial conditions (if needed)
      IF ( Input_Opt%ITS_A_MERCURY_SIM .AND.
     &     Input_Opt%LDYNOCEAN ) THEN
         CALL READ_OCEAN_Hg_RESTART( YYYYMMDD, HHMMSS, Input_Opt )
      ENDIF
#else
      ! Read NetCDF restart file
      CALL READ_GC_RESTART_NC( am_I_Root, YYYYMMDD, HHMMSS, 
     &                         Input_Opt, State_Met, State_Chm, RC )
#endif

      !### Debug
      IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN
         CALL DEBUG_MSG('### READ_GC_RESTART: read file')
      ENDIF

      END SUBROUTINE READ_GC_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gc_restart_nc
!
! !DESCRIPTION: Subroutine READ\_GC\_RESTART\_NC initializes GEOS-Chem 
!  concentrations of advected tracers and species (if applicable) from 
!  a restart file in NetCDF format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_GC_RESTART_NC( am_I_Root,  YYYYMMDD,  HHMMSS,
     &                               Input_Opt, State_Met, State_Chm, 
     &                               RC ) 
!
! !USES:
!     
      USE COMODE_LOOP_MOD,    ONLY : NTSPEC, NAMEGAS
      USE m_netcdf_io_open                         ! netCDF file open
      USE m_netcdf_io_read                         ! netCDF read
      USE m_netcdf_io_readattr                     ! netCDF attribute reads
      USE m_netcdf_io_get_dimlen                   ! netCDF get dimensions
      USE m_netcdf_io_close                        ! netCDF file close
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PHYSCONSTANTS,      ONLY : BOLTZ
      USE TIME_MOD,           ONLY : EXPAND_DATE 
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root  ! Are we on the root CPU?
      INTEGER,        INTENT(IN)    :: YYYYMMDD   ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)    :: HHMMSS     ! hh:mm:ss   GMT time
      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
!
      include "netcdf.inc"
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS: 
!
      INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
!
! !REVISION HISTORY: 
!
!  09 Feb 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, L, N         ! lon, lat, lev, var index
      CHARACTER(LEN=10)    :: SPC_NAME           ! species name 
      CHARACTER(LEN=10)    :: UNIT
      REAL(fp)             :: Background_VV      ! background species conc
      REAL(fp)             :: CONV_FACTOR        ! convert mol/mol -> molec/cm3

      ! NetCDF file info
      INTEGER              :: fId                ! netCDF file ID
      INTEGER              :: NX, NY, NZ, NT     ! netCDF file dimensions
      CHARACTER(LEN=255)   :: nc_file            ! filename

      ! NetCDF variable info
      CHARACTER(LEN=255)   :: v_name             ! variable name 
      CHARACTER(LEN=255)   :: v_attr             ! attribute name
      CHARACTER(LEN=255)   :: v_unit             ! unit value
      INTEGER              :: v_id               ! variable id
      INTEGER              :: ierr               ! 0 means variable in file
      INTEGER              :: st4d(4), ct4d(4)   ! start and count values
                                                 ! for 4D variable array
      ! Temporary arrays and pointers
      REAL*4, TARGET  :: Temp3D( IIPAR, JJPAR, LLPAR )
      REAL*4, POINTER :: Ptr3D(:,:,:) => NULL()

      ! Other
      CHARACTER(LEN=255)   :: LOC

      !=================================================================
      ! READ_GC_RESTART_NC begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS
      
      ! Name of this routine
      LOC = 'READ_GC_RESTART_NC (restart_mod.F)'

      ! Copy input filename to a local variable
      nc_file = TRIM( INPUT_RESTART_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

      !=================================================================
      ! Open the GEOS-Chem restart file and prepare variables for read
      !=================================================================

      ! Open file for read
      CALL Ncop_Rd( fID, TRIM( nc_file ) )

      ! Read the dimensions from the netCDF file
      CALL NcGet_DimLen( fID, 'lon',   NX )
      CALL NcGet_DimLen( fID, 'lat',   NY )
      CALL NcGet_DimLen( fID, 'lev',   NZ )
      CALL NcGet_DimLen( fID, 'time',  NT )

      ! Make sure the dimensions of the file are valid
      CALL Check_Dimensions( lon=NX, lat=NY, lev=NZ, time=NT,  
     &                       time_expected=1, LOC=LOC )

      ! Set 4D start and count indices for dims: lon, lat, lev, time
      st4d   = (/ 1, 1, 1, 1 /)
      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)

      ! Write read message to log
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T '
     &                    // '( NetCDF )'
      WRITE( 6, 100 ) TRIM( nc_file )
 100  FORMAT( 'READ_RESTART_FILE: Reading ', a )

      !=================================================================
      ! Read tracer data, convert units, and store in State_Chm
      !=================================================================

      ! Print header for min/max concentration to log
      WRITE( 6, 110 )
 110  FORMAT( /,'Min and Max of each tracer in file [mol mol-1]:' )

      ! Loop over tracers
      DO N = 1, Input_Opt%N_TRACERS
 
         ! Initialize local array for this tracer and point to it
         Temp3D = 0.0
         Ptr3D => Temp3D

         ! Define variable name
         v_name = 'TRC_' // TRIM( Input_Opt%TRACER_NAME(N) )
         
         ! Read variable from file and store in local array
         CALL NcRd( Ptr3D, fId, TRIM(v_name), st4d, ct4d )

         ! Read the units attribute
         v_attr = "units"
         CALL NcGet_Var_Attributes( fId, TRIM(v_name), TRIM(v_attr), 
     &                              v_unit )

         ! Convert units to v/v if not already
         CALL CONVERT_TRACER_TO_VV( N, Ptr3D, v_unit )

         ! Print the min & max of each tracer as it is read from the file
         WRITE( 6, 120 ) N, TRIM( Input_Opt%TRACER_NAME(N) ),
     &                   MINVAL( Ptr3D ), MAXVAL( Ptr3D )
 120     FORMAT( 'Tracer ', i3, ', ', a8, ': Min = ', es15.9,
     &           '  Max = ',es15.9)
         
         ! Convert tracer concentration from [v/v dry] to [kg/kg dry] and 
         ! store in State_Chm%Tracers
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            State_Chm%Tracers(I,J,L,N) = Ptr3D(I,J,L) / 
     &                                   Input_Opt%TCVV(N) 
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      
         ! Set State_Chm units
         State_Chm%Trac_Units = 'kg/kg dry'

         ! Nullify pointer
         Ptr3D => NULL()
         
      ENDDO

      !=================================================================
      ! Print mean atmospheric kg/kg for each tracer
      !=================================================================

      WRITE( 6, 130 )
 130  FORMAT( /, 'Mean mass mixing ratio per tracer [kg/kg]: ' ) 
      
      DO N = 1, Input_Opt%N_TRACERS
      
         ! For carbon, be sure to use correct unit string
         IF ( INT( Input_Opt%TRACER_MW_G(N) + 0.5 ) == 12 ) THEN
            UNIT = 'kg C/kg'
         ELSE
            UNIT = 'kg/kg  '
         ENDIF
      
         ! Print totals
         WRITE( 6, 150 ) N, TRIM( Input_Opt%TRACER_NAME(N) ), 
     &                   SUM( State_Chm%TRACERS(:,:,:,N) ) / 
     &                   ( IIPAR * JJPAR * LLPAR), ADJUSTL( UNIT )
 150     FORMAT( 'Tracer ', i3, ', ', a8, ':  ', es15.9, 1x, a7)
      ENDDO
            
      !=================================================================
      ! Read species concentrations and store in array
      !=================================================================

      ! Only read species if chemistry simulation
      IF ( Input_Opt%LCHEM .AND. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN 
         
         ! Check if species data is in file by looking for first species
         v_name = 'SPC_' // TRIM( NAMEGAS(1) )
         ierr = Nf_Inq_Varid (fId, v_name, v_id)

         ! Set species in restart file logical flag and write to log
         IF ( ierr /= NF_NOERR ) THEN
            SPC_IN_NC_RST = .FALSE.
            WRITE( 6, 160 )
 160        FORMAT( /,'No species data found in restart file.',/,
     &                'Using default background of each species',
     &                ' [mol mol-1]:' )
         ELSE
            SPC_IN_NC_RST = .TRUE.
            WRITE( 6, 170 )
 170        FORMAT( /,'Species data found in restart file,',/,
     &                'Min and Max of each species after conversion',
     &                ' to [molec/cm3]: ')
         ENDIF

         ! Initialize species to bad number
         State_Chm%Species = -999.e+0_fp

         ! Loop over species
         DO N = 1, NTSPEC(1)
         
            ! If the species is included in the IFAM array, which denotes 
            ! the index # of species in the CSPEC array that are "fake" ND65 
            ! prod/loss families, then skip. 
            !IF ( ANY( IFAM == N ) ) THEN
            !   CYCLE
            !ENDIF
              
            ! Define species name
            SPC_NAME = TRIM( NAMEGAS(N) )

            ! Read species data from file if it is there, otherwise set
            ! default values 
            IF ( SPC_IN_NC_RST ) THEN

               ! Initialize temporary array for this tracer and point to it
               ! use 1e-30 as min (hotp 2/25/09)
               ! NOTE: the 1.0e-30 minimum is for [molec/cm3] in legacy
               Temp3D = 1.0e-30_fp 
               Ptr3D => Temp3D

               ! Define variable name
               v_name = 'SPC_' // TRIM( SPC_NAME )
               
               ! Read variable from file and store in local array
               CALL NcRd( Ptr3D, fId, TRIM(v_name), st4d, ct4d )

               ! Convert [mol/mol] to [molec/cm3]
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, CONV_FACTOR )
               ! Loop over all chemistry grid boxes
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
               
                  ! Set box-dependent unit conversion factor 
                  CONV_FACTOR = 
     &                    State_Met%PMID_DRY(I,J,L) * 1000e+0_fp /
     &                    ( State_Met%T(I,J,L) * BOLTZ * 1e+7_fp )
                                       
                  ! Convert units from [mol/mol] to [molec/cm3/box]
                  State_Chm%Species(I,J,L,N) = 
     &                    Ptr3D(I,J,L) * CONV_FACTOR
               
                  ! Make a small number if concentration is very low
                  State_Chm%Species(I,J,L,N) = 
     &                 MAX( State_Chm%Species(I,J,L,N), 1.0e-30_fp )
               
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

               ! Print the min and max of this species in [molec/cm3]
               WRITE( 6, 180 ) N, TRIM( SPC_NAME ),
     &                         MINVAL( State_Chm%Species(:,:,:,N) ),
     &                         MAXVAL( State_Chm%Species(:,:,:,N) )
 180           FORMAT('Species ', i3, ', ', a9, ': MIN = ', 
     &                 es15.9, ' MAX = ', es15.9)

               ! Nullify pointer
               Ptr3D => NULL()

            ELSE

               ! If species data is not in the restart file, then assign
               ! default background values in units of [mol/mol].
               ! These values were traditionally stored in globchem.dat
               ! and assigned along with a unit conversion to [molec/cm3]
               ! in gasconc called from chemdr. Now chemdr converts these
               ! initial [mol/mol] values to [molec/cm3]. 
               ! Eventually store these values in the species database.
               SELECT CASE ( TRIM( SPC_NAME ) )
                  CASE ('CH2O', 'H2O2', 'HNO2', 'HNO3', 'HNO4') 
                     Background_VV = 4.0e-15_fp
                  CASE ('HO2', 'MO2', 'MP', 'N2O5', 'NO3', 'OH')
                     Background_VV = 4.0e-15_fp
                  CASE ('CH4')
                     Background_VV = 1.7e-06_fp
                  CASE ('CO')
                     Background_VV = 1.0e-07_fp
                  CASE ('CO2')
                     Background_VV = 3.55e-04_fp
                  CASE ('H2')
                     Background_VV = 5.0e-07_fp
                  CASE ('H2O')
                     Background_VV = 1.839e-02_fp
                  CASE ('M')
                     Background_VV = 0.0e+0_fp
                  CASE ('N2')
                     Background_VV = 7.808e-01_fp
                  CASE ('NO', 'NO2')
                     Background_VV = 4.0e-13_fp
                  CASE ('N2O')
                     Background_VV = 3.0e-07_fp
                  CASE ('EMISSION')
                     Background_VV = 1.0e+0_fp
                  CASE ('H', 'O', 'N')
                     Background_VV = 4.0e-20_fp
                  CASE ('O1D')
                     IF ( Input_Opt%LUCX ) THEN
                        Background_VV = 1.0e-15_fp ! benchmark and ucx sims
                     ELSE
                        Background_VV = 4.0e-22_fp
                     ENDIF
                  CASE ('O2')
                     Background_VV = 2.095e-01_fp
                  CASE ('O3')
                     Background_VV = 2.0e-08_fp
                  CASE ('OCS')
                     Background_VV = 9.0e-15_fp
                  CASE DEFAULT
                     Background_VV = 1.0e-20_fp
               END SELECT
               
               ! Set 'DRY' species to all zeroes
               IF ( SPC_NAME(1:3) == 'DRY') THEN
                  Background_VV = 0.0e+0_fp
               ENDIF

               ! Put this array into State_Chm%Species [mol/mol], but
               ! only up to LLCHEM levels
               State_Chm%Species(:,:,1:LLCHEM,N) = Background_VV
               
               ! Print the background value of this species
               WRITE( 6, 190 ) N, TRIM( SPC_NAME ),
     &                         Background_VV
 190           FORMAT('Species ', i3, ', ', a9, ': Background = ', 
     &                 es15.9)

               ! NOTE: keep State_Chm%Species background values as
               ! mol/mol here to preserve location of unit conversion
               ! used for binary-punch restart files. This is temporary
               ! to minimize impact of the introduction of NetCDF restart
               ! files and to minimize differences with FlexChem.

            ENDIF
         ENDDO
      ENDIF

      !=================================================================
      ! Clean up
      !=================================================================

      ! DEBUGGING
      WRITE( 6, 200 ) MINVAL(State_Met%T), MAXVAL(State_Met%T)
 200  FORMAT('State_Chm%T : MIN =', es15.9, ' MAX = ', es15.9 )
      WRITE( 6, 210 ) MINVAL(State_Met%PMID_DRY), 
     &                MAXVAL(State_Met%PMID_DRY)
 210  FORMAT('State_Chm%PMID_DRY : MIN =', es15.9,
     &       ' MAX = ', es15.9 )
      ! END DEBUGGING

      ! Mark end of section in log
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Close netCDF file
      CALL NcCl( fId )
      
      !### Debug
      IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN
         CALL DEBUG_MSG('### READ_GC_RESTART_NC: read file')
      ENDIF

      END SUBROUTINE READ_GC_RESTART_NC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gc_restart_bpch
!
! !DESCRIPTION: Subroutine READ\_GC\_RESTART\_FILE\_BPCH initializes GEOS-Chem 
!  advected tracer concentrations from a restart file in binary punch file 
!  format
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_GC_RESTART_BPCH( YYYYMMDD,  HHMMSS,
     &                                 Input_Opt, State_Met, 
     &                                 State_Chm ) 
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_READ
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE UNITCONV_MOD
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: YYYYMMDD   ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)    :: HHMMSS     ! hh:mm:ss   GMT time
      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
! 
! !REVISION HISTORY: 
!  27 May 1999 - R. Yantosca - Initial version
!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
!        Also reorganize some print statements  (bmy, 10/25/99)
!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!  (5 ) Broke up sections of code into internal subroutines.  Also updated
!        comments & cleaned up a few things. (bmy, 7/17/00)
!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
!  (10) Added updates from amf for tagged Ox run.  Also updated comments
!        and made some cosmetic changes (bmy, 7/3/01)
!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
!        accordingly. (amf, bmy, 9/6/01)
!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
!  (13) Updated comments (bmy, 1/25/02)
!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
!  (17) Add fancy output string (bmy, 4/26/04)
!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
!        and "tracer_mod.f" (bmy, 7/20/04)
!  (19) Remove code for obsolete CO-OH simulation.  Also remove references
!        to CMN_DIAG and TRCOFFSET.   Change tracer name format string to 
!        A10. (bmy, 6/24/05)
!  (20) Updated comments (bmy, 12/16/05)
!  11 Jul 2011 - R. Yantosca - Added ProTeX headers
!  21 Jul 2011 - M. Long     - Now include F77_CMN_SIZE instead of CMN_SIZE
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  25 Mar 2013 - R. Yantosca - Now accept Input_Opt, State_Chm arguments
!  16 Dec 2014 - M. Yannetti - Changed ZTAU0 and ZTAU1 to REAL*8
!  17 Feb 2015 - E. Lundgren - Change main units from kg to kg/kg.
!                              Write mean kg/kg per tracer to log.
!  09 Feb 2016 - E. Lundgren - Change name from READ_RESTART_FILE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL             :: LSPLIT
      INTEGER             :: I, IOS, IU_RST, J, L, N
      INTEGER             :: NCOUNT(NNPAR) 
      REAL*4              :: TRACER(IIPAR,JJPAR,LLPAR)
      REAL(fp)            :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,        NJ,       NL
      INTEGER             :: IFIRST,    JFIRST,   LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED
      CHARACTER(LEN=255)  :: LOC = 'READ_GC_RESTART_BPCH'
     &                             // ' (restart_mod.F)'

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)

      !=================================================================
      ! READ_GC_RESTART_BPCH begins here!
      !=================================================================

      ! Get fields from Input_Opt
      LSPLIT    = Input_Opt%LSPLIT

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0

      ! Initialize pointer for Chemistry State object [kg/kg] (ewl, 2/17/15)
      STT => State_Chm%Tracers

      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================

      ! Find a free file LUN
      IU_RST   = findFreeLUN()
      
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_RESTART_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T ( BPCH )'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_RESTART_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      ! Echo more output
      WRITE( 6, 110 )
 110  FORMAT( /, 'Min and Max of each tracer in file [v/v]:' )

      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      DO 
         READ( IU_RST, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:4' )

         READ( IU_RST, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:5')

         READ( IU_RST, IOSTAT=IOS ) 
     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_restart_file:6')

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================
  
         ! Only process concentration data (i.e. mixing ratio)
         IF ( CATEGORY(1:8) == 'IJ-AVG-$' ) THEN 

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            CALL CHECK_DIMENSIONS( lon=NI, lat=NJ, lev=NL, LOC=LOC )

            ! Convert TRACER from its native units (eg. [ppbv]) to 
            ! [v/v] mixing ratio, and write min/max to file
            CALL CONVERT_TRACER_TO_VV( NTRACER, TRACER, UNIT )

            ! Print the min & max of each tracer
            WRITE( 6, 120 ) NTRACER,  
     &                      TRIM( Input_Opt%TRACER_NAME(NTRACER) ),
     &                      MINVAL( TRACER ), MAXVAL( TRACER )
 120        FORMAT( 'Tracer ', i3, ', ', a8, ': Min = ', es15.9, 
     &              '  Max = ',  es15.9 )

            ! Convert TRACER from [v/v] to [kg/kg] and copy into STT array
            CALL COPY_STT( NTRACER,   TRACER,   NCOUNT,
     &                     Input_Opt, State_Met, State_Chm )

         ENDIF

      ENDDO

      !=================================================================
      ! Examine data blocks and return
      !=================================================================

      ! Check for missing or duplicate data blocks
      CALL CHECK_DATA_BLOCKS( Input_Opt%N_TRACERS, NCOUNT )

      ! Close file
      CLOSE( IU_RST )      

      ! Print mean atmospheric kg/kg for each tracer
      WRITE( 6, 130 )
 130  FORMAT( /, 'Mean mass mixing ratio of each tracer [kg/kg]: ' ) 

      DO N = 1, Input_Opt%N_TRACERS

         ! For carbon, be sure to use correct unit string
         IF ( INT( Input_Opt%TRACER_MW_G(N) + 0.5 ) == 12 ) THEN
            UNIT = 'kg C/kg'
         ELSE
            UNIT = 'kg/kg  '
         ENDIF

         ! Print totals
         WRITE( 6, 150 ) N, TRIM( Input_Opt%TRACER_NAME(N) ), 
     &                   SUM( STT(:,:,:,N) ) / (NI*NJ*NL), 
     &                   ADJUSTL( UNIT )
 150     FORMAT( 'Tracer ', i3, ', ', a8, ':  ', es15.9, 1x, a7)
      ENDDO

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE READ_GC_RESTART_BPCH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_gc_restart
!
! !DESCRIPTION: Subroutine WRITE\_GC\_RESTART creates restart files  
!  for GEOS-Chem advected tracers [units: vol tracer /vol dry air] by
!  calling subroutines that write to either binary punch or NetCDF format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_GC_RESTART( am_I_Root, Input_Opt, 
     &                             YYYYMMDD,  HHMMSS, TAU,       
     &                             State_Chm, State_Met, RC )
! 
! !USES:
!
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE OCEAN_MERCURY_MOD,  ONLY : MAKE_OCEAN_Hg_RESTART
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)  :: YYYYMMDD    ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)  :: HHMMSS      ! hh:mm:ss GMT time
      REAL(f8),       INTENT(IN)  :: TAU         ! TAU value (hrs from 1/1/1985)
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC        ! Success or failure?
! 
! !REVISION HISTORY: 
!  09 Feb 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local Variables      
      INTEGER           :: I, J, L, N
      REAL(f4), TARGET  :: TracerVV
     &                     (IIPAR,JJPAR,LLPAR,Input_Opt%N_TRACERS)
      REAL(f4), POINTER :: TracerVV_ptr(:,:,:,:)


      !=================================================================
      ! WRITE_GC_RESTART begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Calculate tracer concentration in v/v dry air
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, Input_Opt%N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

          ! Convert from mass mixing ratio [kg/kg dry air] to 
          ! molar mixing ratio [v/v dry air]
          TracerVV(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) * 
     &                         Input_Opt%TCVV(N)
         
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Set pointer
      TracerVV_ptr => TracerVV

      !=================================================================
      ! Call routines to write binary-punch or NetCDF restart files 
      !=================================================================
#if defined( BPCH_RST_OUT )                  
      ! Write GEOS-Chem restart file in bpch format
      CALL WRITE_GC_RESTART_BPCH( am_I_Root, Input_Opt, 
     &                             YYYYMMDD, HHMMSS, 
     &                             TAU,      TracerVV_ptr, RC )

      ! Make ocean mercury restart file in bpch format
      IF ( Input_Opt%ITS_A_MERCURY_SIM .and. Input_Opt%LDYNOCEAN ) THEN
         CALL MAKE_OCEAN_Hg_RESTART( YYYYMMDD, HHMMSS, TAU, Input_Opt )
      ENDIF

      ! Save species concentrations (CSPEC_FULL) (bpch) (dkh, 02/12/09)
      IF ( Input_Opt%LCHEM .and. Input_Opt%LSVCSPEC ) THEN 
         CALL MAKE_CSPEC_FILE( am_I_Root, Input_Opt,  YYYYMMDD, 
     &                         HHMMSS,    State_Chm,  RC    )
      ENDIF 
#else
      ! Write GEOS-Chem restart file in NetCDF format
      CALL WRITE_GC_RESTART_NC( am_I_Root, YYYYMMDD, HHMMSS, Input_Opt, 
     &                          State_Chm, State_Met, TracerVV_ptr, RC )
#endif

      ! Free pointer
      NULLIFY( TracerVV_ptr )

      !### Debug
      IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN 
         CALL DEBUG_MSG('### WRITE_GC_RESTART: read file')
      ENDIF

      END SUBROUTINE WRITE_GC_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_gc_restart_nc
!
! !DESCRIPTION: Subroutine WRITE\_GC\_RESTART\_NC creates restart files for 
!  GEOS-Chem advected tracers [units: vol tracer /vol dry air] in NetCDF
!  format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_GC_RESTART_NC( am_I_Root, YYYYMMDD,  HHMMSS,
     &                                Input_Opt, State_Chm, State_Met,
     &                                TracerVV,  RC )
!
! !USES:
!
      USE COMODE_LOOP_MOD,    ONLY : NTSPEC, NAMEGAS
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PHYSCONSTANTS,      ONLY : BOLTZ
      USE TIME_MOD,           ONLY : EXPAND_DATE

      ! HEMCO
      USE HCO_STATE_MOD,      ONLY : HCO_STATE
      USE HCO_Error_Mod
      USE HCO_Diagn_Mod,      ONLY : Diagn_Update
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState
      USE HCOIO_Diagn_Mod,    ONLY : HCOIO_Diagn_WriteOut
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      INTEGER,        INTENT(IN)  :: YYYYMMDD    ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)  :: HHMMSS      ! hh:mm:ss GMT time
      REAL(f4),       INTENT(IN)  :: TracerVV(:,:,:,:) ! Tracer concentrations
                                                       ! [v/v dry air]!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC        ! Success or failure?
! 
! !REVISION HISTORY: 
!  09 Feb 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local Variables      
      INTEGER              :: N, I, J, L, HCRC
      INTEGER              :: GC_RST_Collection
      REAL(fp)             :: CONV_FACTOR        ! convert molec/cm3 -> mol/mol
      REAL*4               :: SMALL_NUM
      REAL*4, TARGET       :: Temp3D(IIPAR, JJPAR, LLPAR)
      REAL*4, POINTER      :: Ptr3D(:,:,:) 
      CHARACTER(LEN=60)    :: ContainerName
      CHARACTER(LEN=255)   :: MSG, LOC = 
     &                       'WRITE_GC_NC_RESTART (restart_mod.F90)'
      CHARACTER(LEN=255)   :: FILENAME


      ! Pointer to HEMCO state (is this needed?)
      TYPE(HCO_STATE), POINTER :: HcoState => NULL()

      !=================================================================
      ! WRITE_GC_RESTART_NC begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Set the GEOS-Chem restart collection ID
      GC_RST_Collection = Input_Opt%GC_RST_COLLECTION

      ! Set minimum value of species concentration [molec/cm3/box]
      SMALL_NUM = 1.e-30_sp

      !=================================================================
      ! Construct output restart filename
      !=================================================================

      ! Copy the output restart file name into a local variable
      FILENAME = TRIM( OUTPUT_RESTART_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Write to log
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   O U T P U T ' 
     &                   // '( NetCDF )'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'WRITE_GC_RESTART_NC: Writing ', a )
      
      !-----------------------------------------------------------------
      ! Populate restart tracer containers before writing tracer data
      !-----------------------------------------------------------------
      ! Print min and max to log for testing
      WRITE( 6, 110 )
 110  FORMAT( /, 'Min and Max of each tracer in file [mol/mol]:' )

      DO N = 1, Input_Opt%N_TRACERS
      
         ! Define container name
         ContainerName = 'TRC_' // TRIM( Input_Opt%TRACER_NAME(N) )

         ! Point to this tracer concentration array
         Temp3D = TracerVV(:,:,:,N)
         Ptr3D => Temp3D
 
         ! Populate with data
         CALL Diagn_Update( am_I_Root,                   
     &                      cName   = TRIM( ContainerName ),
     &                      COL     = GC_RST_Collection,
     &                      Array3D = Ptr3D,             
     &                      RC      = HCRC )
      
         ! Print the min & max of each tracer written to file
         WRITE( 6, 120 ) N, TRIM( Input_Opt%TRACER_NAME(N) ),
     &                   MINVAL( Temp3D ), MAXVAL( Temp3D )
 120     FORMAT( 'Tracer ', i3, ', ', a8, ': Min = ', es15.9, 
     &           '  Max = ',es15.9)

         ! Free the pointer
         Ptr3D => NULL()
      
         IF ( HCRC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot update restart tracer container: ' 
     &            // TRIM( ContainerName )
            CALL ERROR_STOP( MSG, LOC ) 
         ENDIF

      ENDDO

      !-----------------------------------------------------------------
      ! Populate restart species containers before writing species data
      !-----------------------------------------------------------------
      ! Only save species if doing full chemistry simulation
      IF ( Input_Opt%LCHEM .AND. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

         ! Print min/max species conc to log in [molec/cm3]
         ! to compare with legacy bpch cspec restart file for testing
         WRITE( 6, 130 )
 130     FORMAT( /, 'Min and Max of each species before unit ',
     &              'conversion [molec/cm3]:' )
         DO N = 1 , NTSPEC(1)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               IF ( State_Chm%Species(I,J,L,N) >= SMALL_NUM ) THEN
                  Temp3D(I,J,L) = State_Chm%Species(I,J,L,N)
               ELSE
                  Temp3D(I,J,L) = SMALL_NUM  
               ENDIF
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            WRITE( 6, 140 ) N, TRIM( NAMEGAS(N) ),
     &                      MINVAL( Temp3D ), 
     &                      MAXVAL( Temp3D )
 140        FORMAT( 'Species ', i3, ', ', a9, ': Min = ', es15.9, 
     &              '  Max = ', es15.9 )

         ENDDO

         ! DEBUGGING
         WRITE( 6, 200 ) MINVAL(State_Met%T), MAXVAL(State_Met%T)
 200     FORMAT('State_Chm%T : MIN =', es15.9, ' MAX = ', es15.9 )
         WRITE( 6, 210 ) MINVAL(State_Met%PMID_DRY), 
     &                   MAXVAL(State_Met%PMID_DRY)
 210     FORMAT('State_Chm%PMID_DRY : MIN =', es15.9, 
     &          ' MAX = ', es15.9 )
         !END DEBUGGING

         ! Convert species values locally to [mol/mol] for output
         ! in NetCDF. Print min and max values to log.
         WRITE( 6, 150 )
 150     FORMAT( /, 'Min and Max of each species written to file',
     &              ' [mol/mol]:' )
         
         DO N = 1 , NTSPEC(1)
         
            ! Remove check for ND65 prod/loss families pending work on ND65 
            ! diagnostics and incorporation of FlexChem (ewl, 2/4/16)
            !! If the species is included in the IFAM array, which denotes 
            !! the index # of species in the CSPEC array that are "fake" ND65 
            !! prod/loss families, then skip. 
            !IF ( ANY( IFAM == N ) ) THEN
            !   CYCLE
            !ENDIF
         
            ! Clear the species concentration array to bad values
            Temp3D(:,:,:) = -999.e-0_sp
         
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, CONV_FACTOR )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! NOTE: species is output to the binary-punch "spec_rst"
               ! restart file as molec/cm3 but is output to the NetCDF
               ! restart file as mol/mol. This results in small differences
               ! in the initial molec/cm3 species values at the beginning
               ! of the next run that uses the restart file. The difference
               ! is due to the pressure-dependency of the unit conversion 
               ! between molec/cm3 and mol/mol. 
               !
               ! While the met-fields are the same for the end time of 
               ! this run and the start time of the next run, the pressures 
               ! are not the same. This run calculates the dry pressure 
               ! using the "floating pressure" which is influenced by the 
               ! pressure fixer in transport and specific humidity which 
               ! is variable in GEOS-Chem.
               ! 
               ! The next run will also use the "floating pressure" but 
               ! that pressure will be slightly different. The unit 
               ! conversion will happen in chemistry which is at a 
               ! different point in the timestep as this (after transport)
               ! and there will have been fewer iterations for the 
               ! pressure to be "fixed".

               ! Set box-dependent unit conversion factor to convert
               ! from [molec/cm3] to [mol/mol]
               CONV_FACTOR = State_Met%T(I,J,L) * BOLTZ * 1e+7_fp /
     &                       ( State_Met%PMID_DRY(I,J,L) * 1000e+0_fp )

               ! Store species conc in [mol/mol] in local array
               IF ( State_Chm%Species(I,J,L,N) > SMALL_NUM ) THEN
                  Temp3D(I,J,L) = State_Chm%Species(I,J,L,N) *
     &                            CONV_FACTOR
               ELSE
                  Temp3D(I,J,L) = 1.0e-30_fp * CONV_FACTOR  
               ENDIF

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Point to the array
            Ptr3D => Temp3D
            
            ! Define container name
            ContainerName = 'SPC_' // TRIM( NAMEGAS(N) )
            
            ! Populate container
            CALL Diagn_Update( am_I_Root,                   
     &                         cName   = TRIM( ContainerName ),
     &                         COL     = GC_RST_Collection,
     &                         Array3D = Ptr3D,             
     &                         RC      = HCRC )
            
            ! Free the pointer
            Ptr3D => NULL()

            ! Print the min & max of each species written to file
            WRITE( 6, 160 ) N, TRIM( NAMEGAS(N) ),
     &                      MINVAL( Temp3D ), MAXVAL( Temp3D )
 160        FORMAT( 'Species ', i3, ', ', a9, ': Min = ', es15.9, 
     &              '  Max = ', es15.9 )
                        
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot update restart species container: ' 
     &               // TRIM(ContainerName)
               CALL ERROR_STOP( MSG, LOC ) 
            ENDIF  
         ENDDO
      ENDIF

!      !-----------------------------------------------------------------
!      ! If Hg simulation, create restart containers for mercury
!      !-----------------------------------------------------------------
!      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
!
!         !--------------------------------------------------------------
!         ! Total Hg(0), Hg(II), and Hg(P) in ocean
!         !--------------------------------------------------------------
!         DO M = 1, 3
!         
!            SELECT CASE( M )
!               CASE ( 1 )
!                  Prefix = '_Hg0_' 
!               CASE ( 2 )
!                  Prefix = '_Hg2_'
!               CASE ( 3 )
!                  Prefix = '_HgP_' 
!            END SELECT
!         
!            ! Define container name and ID
!            ContainerName = 'OCEAN_HG' // TRIM( Prefix ) // 'Total' 
!            cID = cId + 1
!            
!            ! Create container

!            CALL Diagn_Create( am_I_Root,                     
!     &                         Col       = GC_RST_Collection,         
!     &                         cId       = cID,                 
!     &                         cName     = TRIM( ContainerName ), 
!     &                         SpaceDim  =  2,                
!     &                         OutUnit   = 'kg',           
!     &                         OutOper   = 'Instantaneous', 
!     &                         RC        = RC )
!            
!            IF ( RC /= HCO_SUCCESS ) THEN
!               MSG = 'Cannot create restart species container: ' 
!     &               // TRIM(ContainerName)
!               CALL ERROR_STOP( MSG, LOC ) 
!            ENDIF  
!         ENDDO
!
!         !-----------------------------------------------------------
!         ! Hg snowpack on land and ocean
!         !-----------------------------------------------------------
!         DO M = 1, 4
!            SELECT CASE( M )
!               CASE ( 1 )
!                  Prefix = 'SNOW_HG_OCN_Hg0_'  ! Reducible on ocean
!               CASE ( 2 )
!                  Prefix = 'SNOW_HGN_OCN_Hg0_' ! Non-reducible on ocean
!               CASE ( 3 )
!                  Prefix = 'SNOW_HG_LND_Hg0_'  ! Reducible on land
!               CASE ( 4 )
!                  Prefix = 'SNOW_HGN_LND_Hg0_' ! Non-reducible on land
!            END SELECT
!         
!            ! Loop over all categories to include tagged Hg
!            DO N = 1, N_Hg_CATS
!            
!               ! Define container name and id
!               ! (temperarily append tracer # until name array exists)
!               write( TracerIDStr, '(I5)' ) N
!               ContainerName = TRIM( Prefix ) // TRIM( TracerIDStr )
!               cID = cId + 1
!            
!               ! Create container
!               CALL Diagn_Create( am_I_Root,                     
!     &                            Col       = GC_RST_Collection,         
!     &                            cId       = cID,                 
!     &                            cName     = TRIM( ContainerName ), 
!     &                            SpaceDim  =  2,                
!     &                            OutUnit   = 'kg',           
!     &                            OutOper   = 'Instantaneous', 
!     &                            RC        = RC )
!               
!               IF ( RC /= HCO_SUCCESS ) THEN
!                  MSG = 'Cannot create restart species container: ' 
!     &                  // TRIM(ContainerName)
!                  CALL ERROR_STOP( MSG, LOC ) 
!               ENDIF  
!            ENDDO
!         ENDDO
!
!         !--------------------------------------------------------------
!         ! Additional restart containers for tagged Hg simulation 
!         !--------------------------------------------------------------
!         IF ( Input_Opt%LSPLIT ) THEN
!
!            !-----------------------------------------------------------
!            ! Tagged Hg(0), Hg(II), and Hg(P) in ocean
!            !-----------------------------------------------------------
!            DO M = 1, 3
!               SELECT CASE( M )
!                  CASE ( 1 )
!                     Prefix = '_Hg0_' 
!                  CASE ( 2 )
!                     Prefix = '_Hg2_'
!                  CASE ( 3 )
!                     Prefix = '_HgP_' 
!               END SELECT
!            
!               DO N = 2, N_Hg_CATS
!               
!                  ! Define container name and ID
!                  ! (temperarily append tracer # until name array exists)
!                  write( TracerIDStr, '(I5)' ) N
!                  Suffix = TRIM( TracerIDStr )
!                  ContainerName = 'OCEAN_HG' // TRIM( Prefix ) 
!     &                            // TRIM( Suffix )
!                  cID = cId + 1
!               
!                  ! Create container
!                  CALL Diagn_Create( am_I_Root,                     
!     &                               Col       = GC_RST_Collection,         
!     &                               cId       = cID,                 
!     &                               cName     = TRIM( ContainerName ), 
!     &                               SpaceDim  =  2,                
!     &                               OutUnit   = 'kg',           
!     &                               OutOper   = 'Instantaneous', 
!     &                               RC        = RC )
!                  
!                  IF ( RC /= HCO_SUCCESS ) THEN
!                     MSG = 'Cannot create restart species container: ' 
!     &                     // TRIM(ContainerName)
!                     CALL ERROR_STOP( MSG, LOC ) 
!                  ENDIF  
!               ENDDO
!            ENDDO
!         ENDIF
!      ENDIF

      !-----------------------------------------------------------------
      ! Write out NetCDF using HEMCO 
      !-----------------------------------------------------------------
      ! Get pointer to HEMCO state object.
      CALL GetHcoState( HcoState )
      IF ( .NOT. ASSOCIATED(HcoState) ) THEN
         CALL ERROR_STOP( 'Cannot get HEMCO state object', LOC )
      ENDIF

      CALL HCOIO_Diagn_WriteOut( am_I_Root, 
     &                           HcoState,                                 
     &                           ForceWrite  = .TRUE.,                    
     &                           UsePrevTime = .FALSE.,                    
     &                           COL         = GC_RST_Collection,  
     &                           RC          = RC                )

      !### Debug
      IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN 
         CALL DEBUG_MSG('### WRITE_GC_RESTART_NC: wrote file')
      ENDIF

      END SUBROUTINE WRITE_GC_RESTART_NC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_gc_restart_bpch
!
! !DESCRIPTION: Subroutine WRITE\_GC\_RESTART\_BPCH creates restart files 
!  for GEOS-Chem advected tracers [units: vol tracer /vol dry air] in binary 
!  punch (bpch) format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_GC_RESTART_BPCH( am_I_Root, Input_Opt, 
     &                                  YYYYMMDD,  HHMMSS,   
     &                                  TAU,       TracerVV, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : BPCH2,         GET_MODELNAME
      USE BPCH2_MOD,          ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_XOFFSET,   GET_YOFFSET
      USE TIME_MOD,           ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)  :: YYYYMMDD    ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)  :: HHMMSS      ! hh:mm:ss GMT time
      REAL(f8),       INTENT(IN)  :: TAU         ! TAU value (hrs from 1/1/1985)
      REAL(f4),       INTENT(IN)  :: TracerVV(:,:,:,:) ! Tracer concentrations
                                                       ! [v/v dry air]
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC        ! Success or failure?
! 
! !REVISION HISTORY: 
!  27 May 1999 - R. Yantosca - Initial version
!  (1 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!  (2 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
!        (bmy, 6/22/00)
!  (3 ) Now do not write more than NTRACE data blocks to disk.  
!        Also updated comments. (bmy, 7/17/00)
!  (4 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (5 ) Added to "restart_mod.f".  Also now save the entire grid to the
!        restart file. (bmy, 6/24/02)
!  (6 ) Bug fix: Remove duplicate definition of MM.  This causes compile-time
!        problems on the ALPHA platform. (gcc, bmy, 11/6/02)
!  (7 ) Now references functions GET_OFFSET, GET_YOFFSET from "grid_mod.f".
!        Now references function GET_TAU from "time_mod.f".  Now added a call 
!        to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
!  (8 ) Cosmetic changes (bmy, 4/29/03)
!  (9 ) Now reference STT, N_TRACERS, TCVV from "tracer_mod.f".  Also now
!        remove hardwired output restart filename.   Now references LPRT
!        from "logical_mod.f". (bmy, 7/20/04)
!  (10) Remove references to CMN_DIAG and TRCOFFSET.  Now call GET_HALFPOLAR 
!        from "bpch2_mod.f" to get the HALFPOLAR flag value for GEOS or GCAP 
!        grids. (bmy, 6/28/05)
!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (12) Add TAU to the argument list (bmy, 12/16/05)
!  11 Jul 2011 - R. Yantosca - Added ProTeX headers
!  21 Jul 2011 - M. Long     - Now include F77_CMN_SIZE instead of CMN_SIZE
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  17 Feb 2015 - E. Lundgren - Change program main units from kg to kg/kg
!  02 Feb 2016 - E. Lundgren - Change name from MAKE_RESTART_FILE
!  09 Feb 2016 - E. Lundgren - Pass TracerVV ptr, remove State_Met and
!                              State_Chm from arg list, and print min/max
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local Variables      
      INTEGER              :: I0, IOS, J0, N
      INTEGER              :: YYYY, MM, DD, HH, SS
      REAL(f4)             :: ThisTracerVV(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER              :: IU_RST
      REAL(f4)             :: LONRES, LATRES
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180 = 1
      
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 

      !=================================================================
      ! WRITE_GC_RESTART_BPCH begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Restart File: ' // 
     &           'Instantaneous Tracer Concentrations (v/v)'
      UNIT     = 'v/v'
      CATEGORY = 'IJ-AVG-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Call GET_HALFPOLAR to return the proper value
      ! for either GCAP or GEOS grids (bmy, 6/28/05)
      HALFPOLAR = GET_HALFPOLAR()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the restart file for output -- binary punch format
      !=================================================================

      ! Find a free file LUN
      IU_RST   = findFreeLUN()

      ! Copy the output restart file name into a local variable
      FILENAME = TRIM( OUTPUT_RESTART_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      WRITE( 6, 100 ) TRIM( FILENAME ), IU_RST
 100  FORMAT( '     - MAKE_RESTART_FILE: Writing ', a, ' on unit ', i4 )

      ! Open restart file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   O U T P U T ' 
     &                   // '( BPCH )'
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( 'WRITE_GC_RESTART_BPCH: Writing ', a )
      
      ! Echo more output
      WRITE( 6, 120 )
 120  FORMAT( /, 'Min and Max of each tracer [v/v]:' )

      !=================================================================
      ! Write each tracer to the restart file
      !=================================================================

      ! Loop over all tracers
      DO N = 1, Input_Opt%N_TRACERS

         ! Store data for this tracer only
         ThisTracerVV = TracerVV(:,:,:,N)

         ! Write tracer data to file
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,      LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY,    N,
     &               UNIT,      TAU,       TAU,         RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,       I0+1,            
     &               J0+1,      1,         ThisTracerVV )

         ! Print the min & max of each tracer written to file
         WRITE( 6, 130 ) N, TRIM( Input_Opt%TRACER_NAME(N) ),
     &                   MINVAL( ThisTracerVV ), MAXVAL( ThisTracerVV )
 130     FORMAT( 'Tracer ', i3, ', ', a8, ': Min = ', es15.9,
     &           '  Max = ',es15.9)

      ENDDO  

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN
         CALL DEBUG_MSG( '### WRITE_GC_RESTART_BPCH: wrote file' )
      ENDIF

      END SUBROUTINE WRITE_GC_RESTART_BPCH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_tracer_to_vv
!
! !DESCRIPTION: Subroutine CONVERT\_TRACER\_TO\_VV converts the TRACER array 
!  from its natural units (e.g. ppbv, ppmv) as read from the restart file to 
!  v/v mixing ratio. It also writes tracer min and max mixing ratio to file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CONVERT_TRACER_TO_VV( NTRACER, TRACER, UNIT )
!
! !USES:
!
      USE CHARPAK_MOD, ONLY : TRANUC
!
! !INPUT PARAMETERS: 
!
      ! Tracer number and units 
      INTEGER,          INTENT(IN) :: NTRACER  
      CHARACTER(LEN=*), INTENT(IN) :: UNIT
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Array containing tracer concentrations
      REAL*4,        INTENT(INOUT) :: TRACER(IIPAR,JJPAR,LLPAR)  
! 
! !REVISION HISTORY: 
!  (1 ) Added to "restart_mod.f".  Can now also convert from ppm or ppmv
!        to v/v mixing ratio. (bmy, 6/25/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (3 ) Remove obsolete reference to CMN (bmy, 6/24/05)
!  11 Jul 2011 - R. Yantosca - Added ProTeX headers
!  21 Jul 2011 - M. Long     - Now include F77_CMN_SIZE instead of CMN_SIZE
!  09 Feb 2016 - E. Lundgren - Move min/max log print to calling routines
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CONVERT_TRACER_TO_VV begins here!
      !=================================================================

      ! Convert UNIT to uppercase
      CALL TRANUC( UNIT )
      
      ! Convert from the current unit to v/v
      SELECT CASE ( TRIM( UNIT ) )

         CASE ( '', 'V/V', 'MOL MOL-1', 'MOL/MOL' )
            ! Do nothing, TRACER is already in v/v

         CASE ( 'PPM', 'PPMV', 'PPMC' ) 
            TRACER = TRACER * 1d-6

         CASE ( 'PPB', 'PPBV', 'PPBC' ) 
            TRACER = TRACER * 1d-9

         CASE ( 'PPT', 'PPTV', 'PPTC' )
            TRACER = TRACER * 1d-12

         CASE DEFAULT
            WRITE( 6, '(a)' ) 'Incompatible units in punch file!'
            WRITE( 6, '(a)' ) 'STOP in CONVERT_TRACER_TO_VV'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP

      END SELECT

      END SUBROUTINE CONVERT_TRACER_TO_VV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dimensions
!
! !DESCRIPTION: Subroutine CHECK\_DIMENSIONS checks to see if dimensions 
!  read from the netCDF file match the defined GEOS-Chem dimensions.  If not, 
!  then it will stop the GEOS-Chem simulation with an error message.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Check_Dimensions( lon, lat, lev, time, 
     &                             time_expected, LOC )
!
! !INPUT PARAMETERS:
!
      INTEGER,          OPTIONAL, INTENT(IN)  :: lon           ! Lon dimension
      INTEGER,          OPTIONAL, INTENT(IN)  :: lat           ! Lat dimension
      INTEGER,          OPTIONAL, INTENT(IN)  :: lev           ! Alt dimension
      INTEGER,          OPTIONAL, INTENT(IN)  :: time          ! Time dimension
      INTEGER,          OPTIONAL, INTENT(IN)  :: time_expected ! Expected # of 
                                                               !  time slots
      CHARACTER(LEN=*), INTENT(IN)            :: LOC           ! Caller location
! 
! !REMARKS:
!  Call this routine with keyword arguments, e.g
!     CALL CHECK_DIMENSION( lon=X,  lat=Y,           lev=Z,         &
!                           time=T, time_expected=8, caller=caller )
!
! !REVISION HISTORY:
! 
! !REVISION HISTORY: 
!  (1 ) Added to "restart_mod.f".  Now no longer allow initialization with 
!        less than a globally-sized data block. (bmy, 6/25/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  11 Jul 2011 - R. Yantosca - Added ProTeX headers
!  21 Jul 2011 - M. Long     - Now include F77_CMN_SIZE instead of CMN_SIZE
!  09 Feb 2016 - E. Lundgren - Add args and all args optional except location
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Error message string
      CHARACTER(LEN=255) :: errMsg                  
      
      !=================================================================
      ! CHECK_DIMENSIONS begins here!
      !=================================================================

      ! Error check latitude dimension 
      IF ( PRESENT( lon ) ) THEN
         IF ( lon /= IIPAR ) THEN
            errMsg = 'Incorrect lon dimension in restart file'
            CALL ERROR_STOP( errMsg, LOC )
         ENDIF
      ENDIF


      ! Error check longitude dimension 
      IF ( PRESENT( lat ) ) THEN
         IF ( lat /= JJPAR ) THEN
            errMsg = 'Incorrect lat dimension in restart file'
            CALL ERROR_STOP( errMsg, LOC )
         ENDIF
      ENDIF
      
      
      ! Error check level dimension 
      IF ( PRESENT( lev ) ) THEN
         IF ( lev /= LLPAR ) THEN
            errMsg = 'Incorrect lev dimension in restart file'
            CALL ERROR_STOP( errMsg, LOC )
         ENDIF
      ENDIF

      ! Error check time dimension 
      IF ( PRESENT( time ) .and. PRESENT( time_expected ) ) THEN
         IF ( time /= time_expected ) THEN
            errMsg = 'Incorrect time dimension in restart file'
            CALL ERROR_STOP( errMsg, LOC )
         ENDIF
      ENDIF

      END SUBROUTINE Check_Dimensions
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: copy_stt
! !DESCRIPTION: Subroutine COPY\_STT converts tracer concentrations from [v/v] 
!  to [kg/kg] and then copies the results into the STT tracer array. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COPY_STT( NTRACER,   TRACER,    NCOUNT,
     &                     Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: NTRACER                   ! Tracer #
      REAL*4,         INTENT(IN)    :: TRACER(IIPAR,JJPAR,LLPAR) ! Tracers from 
                                                                 !rst file [v/v]
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,        INTENT(INOUT) :: NCOUNT(NNPAR)       ! # of data blocks
!                                                          !  for each tracer
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY: 
!  (1 ) Added to "restart_mod.f".  Also added parallel loops. (bmy, 6/25/02)
!  (2 ) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
!  (3 ) Now exit if N is out of range (bmy, 4/29/03)
!  (4 ) Now references N_TRACERS, STT & TCVV from "tracer_mod.f" (bmy, 7/20/04)
!  (5 ) Remove call to TRUE_TRACER_INDEX (bmy, 6/24/05)
!  11 Jul 2011 - R. Yantosca - Added ProTeX headers
!  21 Jul 2011 - M. Long     - Now include F77_CMN_SIZE instead of CMN_SIZE
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now pass State_Chm, Input_Opt as arguments
!  17 Feb 2015 - E. Lundgren - Change units from kg to kg/kg
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L, N

      !=================================================================
      ! COPY_STT begins here!
      !=================================================================

      ! Tracer number
      N = NTRACER

      ! Exit if N is out of range
      IF ( N < 1 .or. N > Input_Opt%N_TRACERS ) RETURN

      ! Convert tracer concentration from [v/v dry] to [kg/kg dry] and 
      ! store in State_Chm%Tracers
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         State_Chm%Tracers(I,J,L,N) = TRACER(I,J,L) / Input_Opt%TCVV(N) 

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Set State_Chm units
      State_Chm%Trac_Units = 'kg/kg dry'

      ! Increment the # of records found for tracer N
      NCOUNT(N) = NCOUNT(N) + 1

      END SUBROUTINE COPY_STT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_data_blocks
!
! !DESCRIPTION: Subroutine CHECK\_DATA\_BLOCKS checks to see if we have 
!  multiple or missing data blocks for a given tracer. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_DATA_BLOCKS( NTRACE, NCOUNT )
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: NTRACE          ! # of advected tracers 
      INTEGER, INTENT(IN) :: NCOUNT(NNPAR)   ! Data blocks found per tracer
! 
! !REVISION HISTORY: 
!  25 Jun 2002 - R. Yantosca - Initial version
!  (1 ) Added to "restart_mod.f".  Also now use F90 intrinsic REPEAT to
!        write a long line of "="'s to the screen. (bmy, 6/25/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  11 Jul 2011 - R. Yantosca - Added ProTeX headers
!  21 Jul 2011 - M. Long     - Now include F77_CMN_SIZE instead of CMN_SIZE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N

      !=================================================================
      ! CHECK_DATA_BLOCKS begins here! 
      !=================================================================

      ! Loop over all tracers
      DO N = 1, NTRACE

         ! Stop if a tracer has more than one data block 
         IF ( NCOUNT(N) > 1 ) THEN 
            WRITE( 6, 100 ) N
            WRITE( 6, 120 ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP
         ENDIF
         
         ! Stop if a tracer has no data blocks 
         IF ( NCOUNT(N) == 0 ) THEN
            WRITE( 6, 110 ) N
            WRITE( 6, 120 ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP
         ENDIF
      ENDDO

      ! FORMAT statements
 100  FORMAT( 'More than one record found for tracer : ', i4 )
 110  FORMAT( 'No records found for tracer : ',           i4 ) 
 120  FORMAT( 'STOP in CHECK_DATA_BLOCKS (restart_mod.f)'    )

      END SUBROUTINE CHECK_DATA_BLOCKS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_cspec_file
!
! !DESCRIPTION: Subroutine READ\_CSPEC\_FILE initializes GEOS-Chem chemical
!  species concentrations from a checkpoint file (binary punch file format).
!  The chemical species are saved into the CSPEC array, which is used by
!  the SMVGEAR and KPP solvers.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CSPEC_FILE( am_I_Root, Input_Opt, YYYYMMDD, 
     &                            HHMMSS,    IT_EXISTS, State_Chm, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_READ
      USE COMODE_LOOP_MOD
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE COMODE_MOD,         ONLY : JLOP
      USE FILE_MOD,           ONLY : IOERROR
      USE FILE_MOD,           ONLY : FILE_EXISTS
      USE TIME_MOD,           ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)    :: YYYYMMDD    ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)    :: HHMMSS      ! hh:mm:ss   GMT time
!     
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!     
! !OUTPUT PARAMETERS: 
!
      LOGICAL,        INTENT(OUT)   :: IT_EXISTS   ! Does CSPEC restart exist?
      INTEGER,        INTENT(OUT)   :: RC 
! 
! !REVISION HISTORY: 
!  30 Aug 2004 - D. Henze    - Initial version, based on READ_RESTART  
!  11 Jul 2011 - R. Yantosca - Now skip over ND65 families when reading 
!                              species from the restart.cspec.YYYYMMDDhh file
!  11 Jul 2011 - R. Yantosca - Added ProTeX headers
!  12 Jul 2011 - R. Yantosca - Now read species name from RESERVED slot
!                              of bpch file as an extra error check
!  21 Jul 2011 - M. Long     - Now include F77_CMN_SIZE instead of CMN_SIZE
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  04 Apr 2013 - R. Yantosca - Now declare State_Chm INTENT(INOUT)
!  22 Aug 2013 - R. Yantosca - Add Input_Opt, RC arguments
!  22 Aug 2013 - R. Yantosca - Now get species restart filename from Input_Opt
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete CMN_mod.F
!  16 Dec 2014 - M. Yannetti - Changed ZTAU0 and ZTAU1 to REAL*8
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, IOS, IU_RST, J, L, N, JLOOP, NN, NTL
      INTEGER             :: NCOUNT(NNPAR) 
      REAL*4		  :: TMP(ILONG,ILAT,IPVERT)

      REAL(fp)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: MSG

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1  ! has to be REAL*8
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 1/24/12)
      REAL(fp),  POINTER    :: CSPEC_FULL(:,:,:,:)

      !=================================================================
      ! READ_CSPEC_FILE begins here!
      !=================================================================

      ! Assume success
      RC         = GIGC_SUCCESS

      ! Initialize some variables
      !TMP(:,:,:) = 0e0
      ! use 1e-30 as min (hotp 2/25/09)
      TMP(:,:,:) = 1e-30

      ! Initialize GEOS-Chem species array from Chemistry State object
      ! [# molecules/cm3] (mpayer, 12/6/12)
      CSPEC_FULL => State_Chm%Species

      !=================================================================
      ! Open checkpoint file and read top-of-file header
      !=================================================================

      ! Find a free file LUN
      IU_RST   = findFreeLUN()
     
      ! Copy the output restart file name into a local variable
      FILENAME = TRIM( Input_Opt%SPEC_RST_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to name
      !FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )

      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME )
      ENDIF
 100  FORMAT( '     - READ_CSPEC_FILE: Reading ', a )
 
      ! Check to see if cspec restart file exists
      IT_EXISTS = FILE_EXISTS( FILENAME )
      IF ( .not. IT_EXISTS ) THEN
         NULLIFY( CSPEC_FULL )
         RETURN
      ENDIF 

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'C S P E C   R E S T A R T   F I L E'
     &                    // '   I N P U T  ( BPCH )'
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( 'READ_CSPEC_FILE: Reading ', a )

      ! Echo more output
      WRITE( 6, 120 )
 120  FORMAT( /, 'Min and Max of each species in file',
     &           ' [molec/cm3/box]:' )

      ! Loop over all species in the chemical mechanism
      DO N = 1, NTSPEC(1)

         ! The IFAM array denotes the index # of species in the CSPEC
         ! array that are "fake" ND65 prod/loss families.  If we 
         ! encounter one of these, then we can skip it.
         IF ( ANY( IFAM == N ) ) THEN
            CYCLE
         ENDIF

         ! Read the values of CSPEC
         READ( IU_RST, IOSTAT=IOS )
     &       MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_cspec_file:13' )

         ! Read data block header
         READ( IU_RST, IOSTAT=IOS )
     &         CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &         NTL,      NN,       NL,   IFIRST, JFIRST, LFIRST,
     &         NSKIP

         ! Error check
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_cspec_file:14' )

         ! Read data block
         READ( IU_RST, IOSTAT=IOS )
     &       ( ( ( TMP(I,J,L), I= 1, NTL), J=1,NN ), L = 1, NL)

         ! Error check
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_cspec_file:16' )

         !==============================================================
         ! Assign data from the TMP array to CSPEC
         !==============================================================

         ! Only process checkpoint data 
         IF ( CATEGORY(1:8) == 'IJ-CHK-$' .and.
     &        NTL           == ILONG      .and. 
     &        NN            == ILAT       .and. 
     &        NL            == IPVERT            ) THEN

            ! Also make sure we have the proper species
            IF ( TRIM( NAMEGAS(N) ) /= TRIM( RESERVED ) ) THEN
               IF ( am_I_Root ) THEN
                  WRITE( 6, 200 ) TRIM( NAMEGAS(N) )
               ENDIF
 200           FORMAT( 'Species mismatch for ', a )
               CALL ERROR_STOP( 'Species mismatch!', 
     &                          'READ_CSPEC_FILE (restart_mod.f)' )
            ENDIF

            ! Read data from temporary array into CSPEC_FULL
            ! NOTE: Using N instead of NTRACER will skip over
            ! the "fake" ND65 prod/loss families (bmy, 7/12/11)
            CSPEC_FULL(:,:,:,N) = TMP(:,:,:)

            ! Print the min & max of each species
            WRITE( 6, 300 ) N, TRIM( NAMEGAS(N) ), 
     &                      MINVAL(TMP), MAXVAL(TMP)
 300        FORMAT( 'Species ', i3, ', ', a9, ': Min = ', es15.9, 
     &              '  Max = ',  es15.9 )

         ELSE
            CALL ERROR_STOP(' Restart data is not correct ', 
     &                   ' reading CSPEC, restart_mod')

         ENDIF
 
      ENDDO ! N 

 555  CONTINUE

      ! Close file
      CLOSE( IU_RST )      

      ! Free pointer
      NULLIFY( CSPEC_FULL )

      !### Debug
      IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN
         CALL DEBUG_MSG( '### READ_CSPEC_FILE: read file' )
      ENDIF

      END SUBROUTINE READ_CSPEC_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: make_cspec_file
!
! !DESCRIPTION: Subroutine MAKE\_CSPEC\_FILE writes GEOS-Chem chemical
!  species concentrations into a checkpoint file (binary punch file format).
!  The chemical species are saved into the CSPEC array, which is used by
!  the SMVGEAR and KPP solvers.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MAKE_CSPEC_FILE( am_I_Root, Input_Opt, YYYYMMDD, 
     &                            HHMMSS,    State_Chm, RC        )
!
! !USES:
!
      USE BPCH2_MOD
      USE COMODE_LOOP_MOD
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE COMODE_MOD,         ONLY : JLOP
      USE FILE_MOD,           ONLY : IOERROR
      USE GRID_MOD,           ONLY : GET_XOFFSET
      USE GRID_MOD,           ONLY : GET_YOFFSET
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TIME_MOD,           ONLY : GET_TAU
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)    :: YYYYMMDD    ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)    :: HHMMSS      ! hh:mm:ss   GMT time
      TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  27 Aug 2004 - D. Henze    - Initial version, based on MAKE_RESTART
!  11 Jul 2011 - R. Yantosca - Now skip over ND65 families when writing
!                              species to the restart.cspec.YYYYMMDDhh file!
!  11 Jul 2011 - R. Yantosca - Added ProTeX headers
!  12 Jul 2011 - R. Yantosca - Save species name to restart file using the
!                              RESERVED field of the bpch file
!  21 Jul 2011 - M. Long     - Now include F77_CMN_SIZE instead of CMN_SIZE
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: IU_RST
      INTEGER              :: I,    I0, IOS, J,  J0, L, N, JLOOP
      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH
      CHARACTER(LEN=255)   :: FILENAME

      ! Temporary storage arrays for checkpointed variables
      REAL*4               :: TMP(ILONG, ILAT, IPVERT)

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      ! make HALFPOLAR variable (hotp 2/25/09)
      !INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180 = 1

      INTEGER 		   :: MAX_nitr_max
      INTEGER 		   :: NSOFAR
      
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 1/24/12)
      REAL(fp), POINTER      :: CSPEC_FULL(:,:,:,:)

      !=================================================================
      ! MAKE_CSPEC_FILE begins here!
      !=================================================================

      ! Assume success
      RC       =  GIGC_SUCCESS

      ! Clear some arrays 
      ! use minimum value instead of zero hotp 2/25/09
      !TMP(:,:,:)   = 0e0
      TMP(:,:,:)   = 1e-30

      ! Initialize GEOS-Chem species array from Chemistry State object
      ! [# molecules/cm3] (mpayer, 12/6/12)
      CSPEC_FULL => State_Chm%Species

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Checkpoint File: ' // 
     &           'Instantaneous Species Concentrations (#/cm3)'
      CATEGORY = 'IJ-CHK-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      ! get value of HALFPOLAR hotp 2/25/09
      HALFPOLAR = GET_HALFPOLAR()

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      ! Find a free file LUN
      IU_RST   = findFreeLUN()

      ! Copy the output restart filename template into a local variable
      FILENAME = TRIM( Input_Opt%SPEC_RST_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to filename
      !FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'C S P E C   R E S T A R T   F I L E'
     &                    // '   O U T P U T  ( BPCH )'
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( 'MAKE_CSPEC_FILE: Writing ', a )
      
      ! Echo more output
      WRITE( 6, 120 )
 120  FORMAT( /, 'Min and Max of each species [molec/cm3/box]:' )

      !=================================================================
      ! Write each checkpointed quantity to the checkpoint file
      !=================================================================

      ! Checkpt additional values for full chem simulation
  
      ! Write the final species concetrations after full chemistry
      UNIT = 'molec/cm3/box'
          
      ! Loop over the total # of species.  This also includes the "fake"
      ! prod/loss family species for the ND65 diagnostic (bmy, 7/11/11)
      DO N = 1 , NTSPEC(1)
         
         ! The IFAM array denotes the index # of species in the CSPEC
         ! array that are "fake" ND65 prod/loss families.  If we 
         ! encounter one of these, then we can skip it (bmy, 7/11/11)
         IF ( ANY( IFAM == N ) ) THEN
            CYCLE
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, IPVERT
         DO J = 1, ILAT
         DO I = 1, ILONG
        
            ! fix recommended by pls, save CSPEC_FULL (hotp 2/18/09)
            ! now set conc no smaller than 1d-30
            IF ( CSPEC_FULL(I,J,L,N) .LE. 1e-30_fp ) THEN
                 TMP(I,J,L) = 1e-30
            ELSE
                 TMP(I,J,L) = CSPEC_FULL(I,J,L,N)
            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
 
         ! Save the species name to the RESERVED slot of the bpch file
         RESERVED = TRIM( NAMEGAS(N) )

         ! Write the data block to the CSPEC checkpoint file
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               ILONG,     ILAT,      IPVERT,    I0+1,
     &               J0+1,      1,         TMP(:,:,:) )
 
         ! Print the min & max of each tracer written to file
         WRITE( 6, 130 ) N, TRIM( NAMEGAS(N) ), 
     &                   MINVAL( TMP ), MAXVAL( TMP )
 130     FORMAT( 'Species ', i3, ', ', a9, ': Min = ', es15.9, 
     &           '  Max = ', 
     &           es15.9 )

      ENDDO

      ! Mark end of section in log
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

      ! Close file
      CLOSE( IU_RST )

      ! Free pointer
      NULLIFY( CSPEC_FULL )

      !### Debug
      IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN
         CALL DEBUG_MSG( '### MAKE_CSPEC_FILE: wrote file' )
      ENDIF

      END SUBROUTINE MAKE_CSPEC_FILE
!EOC
      END MODULE RESTART_MOD
