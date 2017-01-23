!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: restart_mod.F
!
! !DESCRIPTION: Module RESTART\_MOD contains variables and routines which 
!  are used to read and write restart files for GEOS-Chem species in units
!  of [v/v] dry mixing ratio, which represents mol species per mol dry air.
!\\
!\\
! !INTERFACE:
!
      MODULE RESTART_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE inquireMod,    ONLY : findFreeLUN
      USE PhysConstants
      USE Precision_Mod
      USE UnitConv_Mod

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS
!
      PUBLIC  :: READ_GC_RESTART
      PUBLIC  :: WRITE_GC_RESTART
      PUBLIC  :: INIT_GC_RESTART
      PUBLIC  :: COMPARE_TRC_SPC ! temporary for tracer removal debugging
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CHECK_RST_DIMENSIONS
!
! !DEFINED PARAMETERS:
!
      ! Prefix of GEOS-Chem output restart filename
      CHARACTER(LEN=255), PARAMETER :: Output_GC_Rst_Prefix = 
     &                                 'GEOSChem_restart'
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
!  11 Jul 2016 - E. Lundgren - Store only species in restart files
!  12 Jul 2016 - E. Lundgren - Remove bpch restart file format option and
!                              associated routines
!  09 Aug 2016 - E. Lundgren - Remove obsolete routines convert_tracer_to_vv
!                              and set_restart and associated module vars;
!                              replace all tracer language to species; rename 
!                              OUTPUT_RESTART_FILE to OUTPUT_GC_RST_PREFIX
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
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
! !IROUTINE: init_gc_restart
!
! !DESCRIPTION: Subroutine INIT\_GC\_RESTART initializes the GEOS-Chem
!  restart file collection and populates it with containers. This collection
!  includes restart concentrations for each species as well as additional
!  data if running the mercury simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GC_RESTART( am_I_Root, Input_Opt, 
     &                            State_Chm, RC         )
!
! !USES:
!
      USE GC_GRID_MOD,        ONLY : AREA_M2
      USE HCO_Error_Mod
      USE HCO_Diagn_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS: 
!
      LOGICAL,          INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput),   INTENT(INOUT) :: Input_Opt   ! Input Options object
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  09 Feb 2016 - E. Lundgren - Initial version
!  20 Apr 2016 - E. Lundgren - Implement ocean and snow Hg variables
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with SpcInfo%Name from species database
!  11 Jul 2016 - E. Lundgren - Remove tracer container; now only use species
!  12 Jul 2016 - E. Lundgren - Rename routine from init_gc_restart_nc
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local Variables  
      INTEGER                :: GC_RST_Collection, N, M, cID
      INTEGER                :: DeltaYMD, DeltaHMS 
      REAL(sp)               :: TS
      REAL(fp), POINTER      :: GridAreas_ptr(:,:)
      CHARACTER(LEN=60)      :: ContainerName, Prefix
      CHARACTER(LEN=255)     :: MSG
      CHARACTER(LEN=255)     :: LOC='INIT_GC_RESTART (restart_mod.F90)'
      TYPE(Species), POINTER :: SpcInfo

      ! For Hg simulation
      INTEGER                   :: Num_Hg_Categories
      CHARACTER(LEN=60)         :: HgSpc
      CHARACTER(LEN=4), POINTER :: Hg_Cat_Name(:)

      !=================================================================
      ! INIT_GC_RESTART begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

      ! Initialize pointers
      SpcInfo       => NULL()
      Hg_Cat_Name   => NULL()

      ! Define collection variables
      GridAreas_ptr => AREA_M2(:,:,1)
      TS            =  GET_TS_CHEM() * 60.0_sp

      ! Create restart collection for GEOS-Chem.  This will keep the
      ! GEOS-Chem restart separate from other output.
#if defined ( ESMF_ )
      deltaYMD = 0 
      deltaHMS = 1
#else
      deltaYMD = 99999999
      deltaHMS = 999999
#endif     
      ! Need to adjust this for restart file
      CALL DiagnCollection_Create( am_I_Root, HcoState%Diagn, 
     &                             NX           = IIPAR,           
     &                             NY           = JJPAR,           
     &                             NZ           = LLPAR,           
     &                             TS           = TS,              
     &                             AM2          = GridAreas_ptr,       
     &                             COL          = GC_RST_Collection,    
     &                             PREFIX       = Output_GC_RST_Prefix,     
     &                             deltaYMD     = deltaYMD,       
     &                             deltaHMS     = deltaHMS,       
     &                             OutTimeStamp = HcoDiagnEnd,     
     &                             RC           = RC                )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP( 'Error creating GEOS-Chem restart collection',
     &                    LOC ) 
      ENDIF
      
      ! Cleanup
      GridAreas_ptr => NULL()
      
      ! Save GEOS-Chem restart collection ID in Input_Opt
      Input_Opt%GC_RST_COLLECTION = GC_RST_Collection

      !=================================================================
      ! If running Hg simulation, set Hg-specific local variables
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Set the # of tagged Hg categories to State_Chm value
         Num_Hg_Categories = State_Chm%N_Hg_CATS

         ! Point to Hg category names in State_Chm
         Hg_Cat_Name => State_Chm%Hg_Cat_Name
         
      ENDIF

      !=================================================================
      ! Create species restart containers
      !=================================================================

      ! Initialize container id
      cID = 0

       ! Loop over the total # of species
      DO N = 1 , State_Chm%nSpecies
      
         ! Get info about this species from the species database
         SpcInfo => State_Chm%SpcData(N)%Info

         ! Define container name and ID
         ContainerName = 'SPC_' // TRIM( SpcInfo%Name )
         cID = cID + 1

         ! Create container
         CALL Diagn_Create( am_I_Root, HcoState, 
     &                      Col       = GC_RST_Collection,         
     &                      cId       = cID,                 
     &                      cName     = TRIM( ContainerName ), 
     &                      SpaceDim  =  3,                
     &                      OutUnit   = 'mol mol-1',           
     &                      OutOper   = 'Instantaneous', 
     &                      RC        = RC )

         ! Free pointer
         SpcInfo => NULL()
         
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot create restart species container: ' 
     &            // TRIM(ContainerName)
            CALL ERROR_STOP( MSG, LOC ) 
         ENDIF  
      ENDDO

      !=================================================================
      ! If Hg simulation, create Hg restart containers
      !=================================================================

      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         !--------------------------------------------------------------
         ! Total Hg(0), Hg(II), and Hg(P) in ocean
         !--------------------------------------------------------------
         DO M = 1, 3
         
            SELECT CASE( M )
               CASE ( 1 )
                  HgSpc = 'Hg0' 
               CASE ( 2 )
                  HgSpc = 'Hg2'
               CASE ( 3 )
                  HgSpc = 'HgP' 
            END SELECT
         
            ! Define container name and ID
            ContainerName = 'Ocean_' // TRIM( HgSpc )
            cID = cID + 1
            
            ! Create container
            CALL Diagn_Create( am_I_Root, HcoState, 
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
                  Prefix = 'Snow_Hg_Ocean_Hg0'  ! Reducible on ocean
               CASE ( 2 )
                  Prefix = 'Snow_HgN_Ocean_Hg0' ! Non-reducible on ocean
               CASE ( 3 )
                  Prefix = 'Snow_Hg_Land_Hg0'   ! Reducible on land
               CASE ( 4 )
                  Prefix = 'Snow_HgN_Land_Hg0'  ! Non-reducible on land
            END SELECT
         
            ! Loop over all categories to include tagged Hg
            DO N = 1, Num_Hg_Categories

               ! Define container name and id
               IF ( N == 1 ) THEN
                  ContainerName = TRIM( Prefix )
               ELSE
                  ! Append category name
                  ContainerName = TRIM( Prefix         ) // '_' //
     &                            TRIM( Hg_Cat_Name(N) )
               ENDIF
               cID = cId + 1
            
               ! Create container
               CALL Diagn_Create( am_I_Root, HcoState, 
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
                     HgSpc = 'Hg0' 
                  CASE ( 2 )
                     HgSpc = 'Hg2'
                  CASE ( 3 )
                     HgSpc = 'HgP' 
               END SELECT
            
               DO N = 2, Num_Hg_Categories
               
                  ! Define container name and ID
                  ContainerName = 'Ocean_' // TRIM( HgSpc ) //
     &                            '_'      // TRIM( Hg_Cat_Name(N) )
                  cID = cId + 1
               
                  ! Create container
                  CALL Diagn_Create( am_I_Root, HcoState, 
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

         ! Free Hg index pointers
         Hg_Cat_Name => NULL()

      ENDIF

      END SUBROUTINE INIT_GC_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gc_restart
!
! !DESCRIPTION: Subroutine READ\_GC\_RESTART reads species concentrations 
!  [mol/mol] from the GEOS-Chem restart file and uses them to initialize 
!  species concentrations in [kg/kg dry]. If species data are missing from
!  the restart file, pre-configured background values are used. If using the
!  mercury simulation, additional restart data are read from file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_GC_RESTART( am_I_Root,  YYYYMMDD,  HHMMSS,
     &                            Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!     
      USE m_netcdf_io_open                         ! netCDF file open
      USE m_netcdf_io_read                         ! netCDF read
      USE m_netcdf_io_readattr                     ! netCDF attribute reads
      USE m_netcdf_io_get_dimlen                   ! netCDF get dimensions
      USE m_netcdf_io_close                        ! netCDF file close
      USE PHYSCONSTANTS,      ONLY : BOLTZ
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : EXPAND_DATE

      ! For Hg simulation restart file
      USE OCEAN_MERCURY_MOD,  ONLY : Hg0aq, Hg2aq, HgPaq, Hgaq_tot
      USE OCEAN_MERCURY_MOD,  ONLY : CHECK_OCEAN_MERCURY
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_OC,  SNOW_HG_STORED_OC
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_LN,  SNOW_HG_STORED_LN
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
!  20 Apr 2016 - E. Lundgren - Implement ocean and snow Hg variables
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  31 May 2016 - E. Lundgren - Replace Input_Opt%TRACER_MW_G with species
!                              database field emMW_g (emitted species g/mol)
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with SpcInfo%Name from species database
!  22 Jun 2016 - R. Yantosca - Now refer to Hg0_Id_List, Hg2_Id_List, and
!                              HgP_Id_List fields of State_Chm
!  11 Jul 2016 - E. Lundgren - Remove tracers and read only species
!  12 Jul 2016 - E. Lundgren - Rename from read_gc_restart_nc
!  18 Jul 2016 - M. Sulprizio- Remove special handling of ISOPN, MMN, CFCX, and
!                              HCFCX. Family tracers have been eliminated.
!  25 Jul 2016 - E. Lundgren - Store whether species in rst file in species db
!                              rather than module-level variable
!  03 Aug 2016 - E. Lundgren - Remove tracers; now only use species
!  11 Aug 2016 - E. Lundgren - Move source of background values to spc database
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, L, M, N      ! lon, lat, lev, indexes
      LOGICAL              :: FOUND              ! ocn Hg var in restart file
      CHARACTER(LEN=60)    :: Prefix             ! utility string
      CHARACTER(LEN=255)   :: LOC                ! routine location
      CHARACTER(LEN=255)   :: MSG                ! message 
      CHARACTER(LEN=255)   :: MyName             ! local name 
      REAL(fp)             :: MW_g               ! species molecular weight
      REAL(fp)             :: SMALL_NUM          ! small number threshold

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
      INTEGER              :: st4d(4), ct4d(4)   ! start and count values for
      INTEGER              :: st3d(3), ct3d(3)   ! 3D and 4D variable arrays
                                                 ! (includes time as dim)
      ! Temporary arrays and pointers
      REAL*4,  TARGET      :: Temp2D(IIPAR,JJPAR) 
      REAL*4,  TARGET      :: Temp3D(IIPAR,JJPAR,LLPAR)
      REAL*4,  POINTER     :: Ptr2D(:,:  )
      REAL*4,  POINTER     :: Ptr3D(:,:,:)

      ! For Hg simulation
      INTEGER                   :: Num_Hg_Categories 
      INTEGER                   :: Total_Hg_Id
      CHARACTER(LEN=60)         :: HgSpc
      CHARACTER(LEN=4), POINTER :: Hg_Cat_Name(:)

      ! Passive tracers
      LOGICAL         :: IsPassive
      REAL(fp)        :: Background_VV

      ! Objects
      TYPE(Species),    POINTER :: SpcInfo

      !=================================================================
      ! READ_GC_RESTART begins here!
      !=================================================================

      ! Assume success
      RC        = GC_SUCCESS

      ! Initialize pointers
      Ptr2D       => NULL()
      Ptr3D       => NULL()
      SpcInfo     => NULL()
      Hg_Cat_Name => NULL()

      ! Name of this routine
      LOC = 'READ_GC_RESTART (restart_mod.F)'

      ! Set minimum value threshold for [mol/mol]
      SMALL_NUM = 1.0e-30_fp

      ! Copy input filename to a local variable
      nc_file = TRIM( Input_Opt%IN_RST_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

      !=================================================================
      ! If running Hg simulation, set Hg-specific local variables
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Set the # of tagHg categories from State_Chm
         Num_Hg_Categories   =  State_Chm%N_Hg_CATS

         ! Set variable storing names for each of the Hg categories
         Hg_Cat_Name => State_Chm%Hg_Cat_Name
    
         ! Set Hg species index corresponding to a given Hg category number;
         ! total is always the first category
         Total_Hg_Id   =  State_Chm%Hg0_Id_List(1)

      ENDIF

      !=================================================================
      ! Open GEOS-Chem restart file
      !=================================================================

      ! Open file for read
      CALL Ncop_Rd( fID, TRIM( nc_file ) )

      ! Read the dimensions from the netCDF file
      CALL NcGet_DimLen( fID, 'lon',   NX )
      CALL NcGet_DimLen( fID, 'lat',   NY )
      CALL NcGet_DimLen( fID, 'lev',   NZ )
      CALL NcGet_DimLen( fID, 'time',  NT )

      ! Make sure the dimensions of the file are valid
      CALL Check_Rst_Dimensions( lon=NX, lat=NY, lev=NZ, time=NT,  
     &                           time_expected=1, LOC=LOC )

      ! Set 3D start and count indices for dims: lon, lat, time
      st3d   = (/ 1, 1, 1 /)
      ct3d   = (/ IIPAR, JJPAR, 1 /)

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
      ! Read species concentrations from NetCDF or use default 
      ! background [mol/mol]; store in State_Chm%Species in [kg/kg dry]
      !=================================================================

      ! IMPORTANT NOTE: the unit conversion from mol/mol to kg/kg uses
      ! the molecular weight stored in the species database which is
      ! a meaningful value for advected species but is a bad value (-1)
      ! for all others. Non-advected species should NOT be used when 
      ! State_Chm%Species units are in mass mixing ratio. Current
      ! units can be determined at any point by looking at 
      ! State_Chm%Spc_Units. (ewl, 8/11/16)

      ! Print header for min/max concentration to log
      WRITE( 6, 110 )
 110  FORMAT( /,'Min and Max of each species in file [mol/mol]:' )

      ! Initialize species to all zeroes
      State_Chm%Species = 0.e+0_fp

      ! Loop over species
      DO N = 1, State_Chm%nSpecies

         ! Get info about this species from the species database
         SpcInfo => State_Chm%SpcData(N)%Info
         MW_g    = State_Chm%SpcData(N)%Info%emMW_g

         ! Define variable name
         v_name = 'SPC_' // TRIM( SpcInfo%Name )

         ! Check if species data is in file
         ierr = Nf_Inq_Varid (fId, v_name, v_id)
         IF ( ierr /= NF_NOERR ) THEN
            SpcInfo%Is_InRestart = .FALSE.
         ELSE
            SpcInfo%Is_InRestart = .TRUE.
         ENDIF

         ! If data is in file, read in as [mol/mol] and convert to 
         ! [kg/kg dry]. Otherwise, set to background value [mol/mol]
         ! either stored in species database (advected species all levels and
         ! non-advected species levels up to LLCHEM) or a small number
         ! (non-advected species levels above LLCHEM) converted to 
         ! [kg/kg dry]
         IF ( SpcInfo%Is_InRestart ) THEN

            ! Initialize temporary array for this species and point to it
            Temp3D = 0.0e-0_fp
            Ptr3D => Temp3D

            ! Read variable from file and store in local array
            CALL NcRd( Ptr3D, fId, TRIM(v_name), st4d, ct4d )

            ! Print the min & max of each species as it is read from 
            ! the restart file in mol/mol
            IF ( Input_Opt%LPRT ) THEN
               WRITE( 6, 120 ) N, TRIM( SpcInfo%Name ),
     &                         MINVAL( Ptr3D ), MAXVAL( Ptr3D )
 120           FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9,
     &                 '  Max = ',es15.9)
            ENDIF

            ! Convert file value [mol/mol] to [kg/kg dry] for storage
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               ! Apply minimum value threshold where input conc is very low
               IF ( Ptr3D(I,J,L) < SMALL_NUM ) THEN
                  Ptr3D(I,J,L) = SMALL_NUM
               ENDIF
               State_Chm%Species(I,J,L,N) = Ptr3D(I,J,L) * MW_g / AIRMW
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE

            ! Set species to the background value converted to [kg/kg dry] 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            ! Loop over all grid boxes
            DO L = 1, LLPAR 
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               
               ! Special handling for MOH
               IF ( TRIM( SpcInfo%Name ) == 'MOH' ) THEN

                  !----------------------------------------------------
                  ! For methanol (MOH), use different initial
                  ! background concentrations for different regions of
                  ! the atmosphere:
                  !
                  ! (a) 2.0 ppbv MOH -- continental boundary layer
                  ! (b) 0.9 ppbv MOH -- marine boundary layer
                  ! (c) 0.6 ppbv MOH -- free troposphere
                  !
                  ! The concentrations listed above are from Heikes et
                  ! al, "Atmospheric methanol budget and ocean
                  ! implication", _Global Biogeochem. Cycles_, 2002.
                  ! These represent the best estimates for the methanol
                  ! conc.'s in the troposphere based on various
                  ! measurements.
                  !
                  ! MOH is an inactive chemical species in GEOS-CHEM,
                  ! so these initial concentrations will never change.
                  ! However, MOH acts as a sink for OH, and therefore
                  ! will affect both the OH concentration and the
                  ! methylchloroform lifetime.
                  !
                  ! We specify the MOH concentration as ppbv, but then
                  ! we need to multiply by CONV_FACTOR in order to
                  ! convert to [molec/cm3].  (bdf, bmy, 2/22/02)
                  !----------------------------------------------------
                  
                  ! Test for altitude (L < 9 is always in the trop)
                  IF ( L <= 9 ) THEN
                     ! Test for ocean/land boxes
                     IF ( State_Met%FRCLND(I,J) >= 0.5 ) THEN
                        ! Continental boundary layer: 2 ppbv MOH
                        State_Chm%Species(I,J,L,N) = 2.000e-9_fp 
     &                                               * MW_g / AIRMW
                     ELSE
                        ! Marine boundary layer: 0.9 ppbv MOH
                        State_Chm%Species(I,J,L,N) = 0.900e-9_fp 
     &                                               * MW_g / AIRMW
                     ENDIF
                  ELSE
                     ! Test for troposphere
                     IF ( ITS_IN_THE_TROP(I,J,L,State_Met) ) THEN
                        ! Free troposphere: 0.6 ppbv MOH
                        State_Chm%Species(I,J,L,N) = 0.600e-9_fp 
     &                                               * MW_g / AIRMW
                     ELSE
                        ! Strat/mesosphere:
                        State_Chm%Species(I,J,L,N) = 
     &                                     SMALL_NUM * MW_g / AIRMW
                     ENDIF
                  ENDIF

                  ! Print to log if debugging is on
                  IF ( Input_Opt%LPRT .AND. I == 1 .AND. J == 1 
     &                 .AND. L == 1 ) THEN
                     WRITE( 6, 130 ) N, TRIM( SpcInfo%Name )
 130                 FORMAT('Species ', i3, ', ', a9, 
     &                   ': see READ_GC_RESTART for special MOH values')
                  ENDIF

               ! For non-advected species at levels above LLCHEM, use a 
               ! small number for background
               ELSEIF ( L > LLCHEM .AND. 
     &                  ( .NOT. SpcInfo%Is_Advected ) ) THEN

                  State_Chm%Species(I,J,L,N) = 
     &                                SMALL_NUM * MW_g / AIRMW 

               ! For all other cases except MOH, use the background value  
               ! stored in the species database
               ELSE

                  State_Chm%Species(I,J,L,N) = SpcInfo%BackgroundVV 
     &                                         * MW_g / AIRMW

                  ! Print to log if debugging is on
                  IF ( Input_Opt%LPRT .AND. I == 1 .AND. J == 1
     &                 .AND. L == 1 ) THEN
                     WRITE( 6, 140 ) N, TRIM( SpcInfo%Name ),
     &                               SpcInfo%BackgroundVV
 140                 FORMAT('Species ', i3, ', ', a9, 
     &                    ': Use background = ', es15.9)
                  ENDIF

               ENDIF

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ENDIF

         ! Free pointer
         SpcInfo => NULL()

      ENDDO

      ! Set species units
      State_Chm%Spc_Units = 'kg/kg dry'

      ! If in debug mode, print out species min and max in [molec/cm3]
      IF ( am_I_Root .and. Input_Opt%LPRT ) THEN
         PRINT *, " "
         PRINT *, "Species min and max in molec/cm3"
         CALL ConvertSpc_KgKgDry_to_MND( am_I_Root, State_Met,
     &                                   State_Chm, RC )
         DO N = 1, State_Chm%nSpecies
            SpcInfo => State_Chm%SpcData(N)%Info
            WRITE(6,150) N, TRIM( SpcInfo%Name ), 
     &                   MINVAL( State_Chm%Species(:,:,:,N) ), 
     &                   MAXVAL( State_Chm%Species(:,:,:,N) ) 
 150        FORMAT( 'Species ', i3, ', ', a9, 
     &              ': Min = ', es15.9, ', Max = ', es15.9 )
            SpcInfo => NULL()
         ENDDO
         CALL ConvertSpc_MND_to_KgKgDry( am_I_Root, State_Met,
     &                                   State_Chm, RC )
      ENDIF

      !=================================================================
      ! Read ocean mercury variables
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Print total mass to log
         WRITE( 6, 160 )
 160     FORMAT(/, 'Total mass of each ocean and snow Hg species:')

         !--------------------------------------------------------------
         ! Total Hg in ocean
         !--------------------------------------------------------------
         DO M = 1, 3

            ! Define variable name
            SELECT CASE( M )
               CASE ( 1 )
                  HgSpc    = 'Hg0' 
               CASE ( 2 )
                  HgSpc    = 'Hg2'
               CASE ( 3 )
                  HgSpc    = 'HgP' 
            END SELECT
            v_name = 'Ocean_' // TRIM( HgSpc )

            ! Zero temporary array and set pointer
            Temp2D = 0.0d4 
            Ptr2D => Temp2D

            ! Check if variable is in file. If not found, initialize
            ! to zero (currently the case for tagged Hg) (ewl, 4/19/16)
            ierr = Nf_Inq_Varid (fId, TRIM(v_name), v_id)
            IF ( ierr /= NF_NOERR ) THEN
               FOUND = .FALSE.
               WRITE( 6, 170 ) TRIM( v_name )
            ELSE
               FOUND = .TRUE.
               CALL NcRd( Ptr2D, fId, v_name, st3d, ct3d )
            ENDIF

            ! Check for negative concentrations (jaf, 7/6/11)
            DO I = 1, IIPAR
            DO J = 1, JJPAR
               IF ( Temp2D(I,J) < 0.0d4 ) THEN
                  Temp2D(I,J) = 0.0d4
               ENDIF
            ENDDO
            ENDDO

            ! Assign ocean mercury data and write total mass to log file
            SELECT CASE( M )
               CASE ( 1 )
                  Hg0aq(:,:,Total_Hg_Id) = Ptr2D
                  IF ( FOUND ) THEN
                     WRITE( 6, 180 ) v_name, 
     &                               SUM( Hg0aq(:,:,Total_Hg_Id) ), 'kg'
                  ENDIF
               CASE ( 2 )
                  Hg2aq(:,:,Total_Hg_Id) = Ptr2D
                  IF ( FOUND ) THEN
                     WRITE( 6, 180 ) v_name, 
     &                               SUM( Hg2aq(:,:,Total_Hg_Id) ), 'kg'
                  ENDIF
               CASE ( 3 )
                  HgPaq(:,:,Total_Hg_Id) = Ptr2D
                  IF ( FOUND ) THEN
                     WRITE( 6, 180 ) v_name, 
     &                               SUM( HgPaq(:,:,Total_Hg_Id) ), 'kg'
                  ENDIF
            END SELECT

            ! Nullify pointer
            Ptr2D => NULL()
         ENDDO

         !--------------------------------------------------------------
         ! Additional tagged ocean Hg species
         !--------------------------------------------------------------
         IF ( Input_Opt%LSPLIT ) THEN
            DO M = 1, 3
               DO N = 2, Num_Hg_Categories

                  ! Define variable name. Include appended region.
                  SELECT CASE( M )
                     CASE ( 1 )
                        HgSpc = 'Hg0' 
                     CASE ( 2 )
                        HgSpc = 'Hg2'
                     CASE ( 3 )
                        HgSpc = 'HgP' 
                  END SELECT
                  v_name = 'Ocean_' // TRIM( HgSpc ) //
     &                     '_'      // TRIM( Hg_Cat_Name(N) )

                  ! Zero temporary array and set pointer
                  Temp2D = 0.0d4
                  Ptr2D => Temp2D
                  
                  ! Check if variable is in file. If not found, initialize
                  ! to zero (currently the case for tagged Hg) (ewl, 4/19/16)
                  ierr = Nf_Inq_Varid (fId, v_name, v_id)
                  IF ( ierr /= NF_NOERR ) THEN
                     FOUND = .FALSE.
                     WRITE( 6, 170 ) TRIM( v_name )
                  ELSE
                     FOUND = .TRUE.
                     CALL NcRd( Ptr2D, fId, v_name, st3d, ct3d )
                  ENDIF

                  ! Assign ocean mercury data and write total mass to log
                  SELECT CASE( M )
                     CASE ( 1 )
                        Hg0aq(:,:,N) = Ptr2D
                        IF ( FOUND ) THEN
                           WRITE( 6, 180 ) v_name, 
     &                                     SUM( Hg0aq(:,:,N) ), 'kg'
                        ENDIF
                     CASE ( 2 )
                        Hg2aq(:,:,N) = Ptr2D
                        IF ( FOUND ) THEN
                           WRITE( 6, 180 ) v_name, 
     &                                     SUM( Hg2aq(:,:,N) ), 'kg'
                        ENDIF
                     CASE ( 3 )
                        HgPaq(:,:,N) = Ptr2D
                        IF ( FOUND ) THEN
                           WRITE( 6, 180 ) v_name, 
     &                                     SUM( HgPaq(:,:,N) ), 'kg'
                        ENDIF
                  END SELECT

                  ! Nullify pointer
                  Ptr2D => NULL()
               ENDDO
            ENDDO

            ! Make sure tagged & total species sum up
            IF ( Input_Opt%USE_CHECKS ) THEN
               CALL CHECK_OCEAN_MERCURY( 'end of READ_GC_RESTART' )
            ENDIF
         ENDIF

         !--------------------------------------------------------------
         ! Hg snowpack on land and ocean
         !--------------------------------------------------------------
         DO M = 1, 4               
            DO N = 1, Num_Hg_Categories               

               ! Define variable name prefix
               SELECT CASE( M )
                  CASE ( 1 )
                     Prefix = 'Snow_Hg_Ocean_Hg0'  ! Reducible on ocean
                  CASE ( 2 )
                     Prefix = 'Snow_HgN_Ocean_Hg0' ! Non-reducible on ocean
                  CASE ( 3 )
                     Prefix = 'Snow_Hg_Land_Hg0'   ! Reducible on land
                  CASE ( 4 )
                     Prefix = 'Snow_HgN_Land_Hg0'  ! Non-reducible on land
               END SELECT

               IF ( N == 1 ) THEN
                  v_name = TRIM( Prefix )
               ELSE
                  ! Append category name if tagged
                  v_name = TRIM( Prefix         ) // '_' //
     &                     TRIM( Hg_Cat_Name(N) ) 
               ENDIF

               ! Zero temporary array and set pointer
               Temp2D = 0.0d4 
               Ptr2D => Temp2D

               ! Check if variable is in file. If not found, initialize
               ! to zero (currently the case for tagged Hg) (ewl, 4/19/16)
               ierr = Nf_Inq_Varid (fId, v_name, v_id)
               IF ( ierr /= NF_NOERR ) THEN
                  FOUND = .FALSE.
                  WRITE( 6, 170 ) TRIM( v_name )
               ELSE
                  FOUND = .TRUE.
                  CALL NcRd( Ptr2D, fId, v_name, st3d, ct3d )
               ENDIF

               ! Assign ocean mercury data and write total mass to file
               SELECT CASE( M )
                  CASE ( 1 )
                     SNOW_HG_OC(:,:,N) = Ptr2D
                     IF ( FOUND ) THEN
                        WRITE( 6, 180 ) v_name, 
     &                        SUM( SNOW_HG_OC(:,:,N) ), 'kg'
                     ENDIF
                  CASE ( 2 )
                     SNOW_HG_STORED_OC(:,:,N) = Ptr2D
                     IF ( FOUND ) THEN
                        WRITE( 6, 180 ) v_name, 
     &                        SUM( SNOW_HG_STORED_OC(:,:,N) ), 'kg'
                     ENDIF
                  CASE ( 3 )
                     SNOW_HG_LN(:,:,N) = Ptr2D
                     IF ( FOUND ) THEN
                        WRITE( 6, 180 ) v_name, 
     &                     SUM( SNOW_HG_LN(:,:,N) ), 'kg'
                     ENDIF
                  CASE ( 4 )
                     SNOW_HG_STORED_LN(:,:,N) = Ptr2D
                     IF ( FOUND ) THEN
                        WRITE( 6, 180 ) v_name, 
     &                        SUM( SNOW_HG_STORED_LN(:,:,N) ), 'kg'
                     ENDIF
               END SELECT

               ! Nullify pointer
               Ptr2D => NULL()

            ENDDO
         ENDDO

         ! Format strings
 170     FORMAT('Variable not found in restart file: ', a24)
 180     FORMAT( a24, ':   ', es15.9, 1x, a4)

         ! Print note that variables are initialized to zero if not 
         ! found (currently only happens in tagged Hg simulation)
         IF ( Input_Opt%LSPLIT ) THEN
            WRITE( 6, 190 )
 190        FORMAT( /, 'NOTE: all variables not found are initialized ',
     &                 'to zero') 
         ENDIF

         ! Free pointers for Hg indexing
         Hg_Cat_Name => NULL()

      ENDIF

      !=================================================================
      ! Clean up
      !=================================================================
      
      ! Mark end of section in log
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Close netCDF file
      CALL NcCl( fId )
      
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
! !IROUTINE: write_gc_restart
!
! !DESCRIPTION: Subroutine WRITE\_GC\_RESTART writes species concentrations 
!  [mol/mol] to the GEOS-Chem restart file in NetCDF format. If using the
!  mercury simulation, additional data such as ocean mercury values 
!  are written to the restart file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_GC_RESTART( am_I_Root, YYYYMMDD,  HHMMSS,
     &                             Input_Opt, State_Chm, State_Met,
     &                             RC )
!
! !USES:
!
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : EXPAND_DATE

      ! HEMCO
      USE HCO_Error_Mod
      USE HCO_Diagn_Mod,      ONLY : Diagn_Update
      USE HCOIO_Diagn_Mod,    ONLY : HCOIO_Diagn_WriteOut

      ! For Hg simulation
      USE OCEAN_MERCURY_MOD,  ONLY : Hg0aq, Hg2aq, HgPaq, Hgaq_tot 
      USE OCEAN_MERCURY_MOD,  ONLY : CHECK_OCEAN_MERCURY     
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_OC,  SNOW_HG_STORED_OC
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_LN,  SNOW_HG_STORED_LN
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      INTEGER,        INTENT(IN)    :: YYYYMMDD    ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)    :: HHMMSS      ! hh:mm:ss GMT time
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  09 Feb 2016 - E. Lundgren - Initial version
!  20 Apr 2016 - E. Lundgren - Implement ocean and snow Hg variables
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with SpcInfo%Name from species database
!  22 Jun 2016 - R. Yantosca - Now refer to Hg0_Id_List, Hg2_Id_List, and
!                              HgP_Id_List fields of State_Chm
!  11 Jul 2016 - E. Lundgren - Remove tracers and write only species
!  12 Jul 2016 - E. Lundgren - Rename from write_gc_restart_nc
!  10 Aug 2016 - E. Lundgren - Input species concentrations are now kg/kg dry
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local Variables      
      INTEGER              :: M, N, HCRC
      INTEGER              :: GC_RST_Collection
      REAL(fp)             :: MW_g
      REAL(fp)             :: SMALL_NUM
      CHARACTER(LEN=60)    :: ContainerName, Prefix
      CHARACTER(LEN=255)   :: MSG, FILENAME, LOC = 
     &                       'WRITE_GC_RESTART (restart_mod.F90)'

      ! Temporary arrays and pointers
      REAL*4,            TARGET :: Temp3D( IIPAR, JJPAR, LLPAR )
      REAL*4,            TARGET :: Temp2D( IIPAR, JJPAR )
      REAL*4,           POINTER :: Ptr3D(:,:,:)
      REAL*4,           POINTER :: Ptr2D(:,:  )

      ! For Hg simulation
      INTEGER                   :: Num_Hg_Categories
      INTEGER                   :: Total_Hg_Id
      CHARACTER(LEN=60)         :: HgSpc
      CHARACTER(LEN=4), POINTER :: Hg_Cat_Name(:)

      ! Objects
      TYPE(Species),    POINTER :: SpcInfo

      !=================================================================
      ! WRITE_GC_RESTART begins here!
      !=================================================================

      ! Assume success
      RC          = GC_SUCCESS

      ! Initialize pointers
      Ptr3D       => NULL()
      Ptr2D       => NULL()
      Hg_Cat_Name => NULL()
      SpcInfo     => NULL()

      ! Set the GEOS-Chem restart collection ID
      GC_RST_Collection = Input_Opt%GC_RST_COLLECTION

      ! Set minimum value to write to file [mol/mol]
      SMALL_NUM = 1.0e-30_fp

      !=================================================================
      ! If running Hg simulation, set Hg-specific local variables
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN 
         
         ! Set the # of tagged Hg categories to the State_Chm value
         Num_Hg_Categories = State_Chm%N_Hg_CATS

         ! Point to names for each of the Hg categories
         Hg_Cat_Name => State_Chm%Hg_Cat_Name

         ! Store Hg species index corresponding to a given Hg category
         ! number; total is always the first category
         Total_Hg_Id = State_Chm%Hg0_Id_List(1)  

      ENDIF

      !=================================================================
      ! Construct output restart filename
      !=================================================================

      ! Define output restart filename with date tokens
      FILENAME = TRIM( Output_GC_RST_Prefix ) // '.YYYYMMDDhhmm'

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Write to log
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   O U T P U T ' 
     &                   // '( NetCDF )'
      WRITE( 6, 200 ) TRIM( FILENAME )
 200  FORMAT( 'WRITE_GC_RESTART: Writing ', a )
      
      !=================================================================
      ! Populate species concentration restart containers
      !=================================================================

      ! Print min/max species concentrations in molec/cm3 to log if debugging
      IF ( Input_Opt%LPRT ) THEN
         WRITE( 6, 210 )
 210     FORMAT( /, 'Min and Max of each species before unit ',
     &              'conversion [molec/cm3]:' )
         CALL ConvertSpc_KgKgDry_to_MND( am_I_Root, State_Met,
     &                                   State_Chm, RC )
         DO N = 1, State_Chm%nSpecies
            SpcInfo => State_Chm%SpcData(N)%Info
            WRITE(6,220) N, TRIM( SpcInfo%Name ), 
     &                   MINVAL( State_Chm%Species(:,:,1:LLCHEM,N) ), 
     &                   MAXVAL( State_Chm%Species(:,:,1:LLCHEM,N) ) 
 220        FORMAT( 'Species ', i3, ', ', a9, 
     &              ': Min = ', es15.9, ', Max = ', es15.9 )
            SpcInfo => NULL()
         ENDDO
         CALL ConvertSpc_MND_to_KgKgDry( am_I_Root, State_Met,
     &                                   State_Chm, RC )

         ! Write header for printing [mol/mol] concentrations later on
         IF ( Input_Opt%LPRT ) THEN
            WRITE( 6, 230 )
 230        FORMAT( /, 'Min and Max of each species written to file',
     &                 ' [mol/mol]:' )
         ENDIF

      ENDIF

      ! Initialize temporary array to all zeros
      Temp3D = 0.e+0_sp

      ! For each species, populate the temporary array with 
      ! concentrations in [mol/mol]
      DO N = 1, State_Chm%nSpecies

         ! Get info about this species from the species database
         SpcInfo => State_Chm%SpcData(N)%Info
         MW_g    = State_Chm%SpcData(N)%Info%emMW_g

         ! Convert kg/kg to mol/mol
         Temp3D = State_Chm%Species(:,:,:,N) *
     &                      AIRMW / MW_g
         
         ! Apply minimum value threshold where output concentration
         ! is very low
         WHERE (Temp3D < SMALL_NUM) Temp3D = SMALL_NUM
 
         ! Point to the array
         Ptr3D => Temp3D
         
         ! Define container name
         ContainerName = 'SPC_' // TRIM( SpcInfo%Name )
        
         ! Populate container
         CALL Diagn_Update( am_I_Root, HcoState,
     &                      cName      = TRIM( ContainerName ),
     &                      COL        = GC_RST_Collection,
     &                      Array3D    = Ptr3D,
     &                      RC         = HCRC )

         ! Free the pointer
         Ptr3D => NULL()

         ! Write min and max written to file if debugging
         IF ( Input_Opt%LPRT ) THEN
            WRITE( 6, 240 ) N, TRIM( SpcInfo%Name ),
     &                      MINVAL( Temp3D(:,:,:) ),
     &                      MAXVAL( Temp3D(:,:,:) )
 240        FORMAT( 'Species ', i3, ', ', a9, ': Min = ', es15.9,
     &              '  Max = ', es15.9 )
         ENDIF

         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot update restart species container: '
     &            // TRIM(ContainerName)
            CALL ERROR_STOP( MSG, LOC )
         ENDIF

         ! Free pointer
         SpcInfo => NULL()        
 
      ENDDO

      !=================================================================
      ! If Hg simulation, populate ocean mercury restart containers
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Print total mass to log
         WRITE( 6, 250 )
 250     FORMAT( /, 'Total mass of each ocean and snow Hg species: ' ) 

         !--------------------------------------------------------------
         ! Total Hg in ocean
         !--------------------------------------------------------------
         DO M = 1, 3

            ! Define container name, point to data, and write total 
            ! mass to log
            Temp2D = 0.0d4
            SELECT CASE( M )
               CASE ( 1 )
                  ContainerName = 'Ocean_Hg0' 
                  Temp2D        = Hg0aq(:,:,Total_Hg_Id)
                  WRITE( 6, 260 ) ContainerName, 
     &                            SUM( Hg0aq(:,:,Total_Hg_Id) ), 'kg'
               CASE ( 2 )
                  ContainerName = 'Ocean_Hg2'
                  Temp2D        = Hg2aq(:,:,Total_Hg_Id)
                  WRITE( 6, 260 ) ContainerName, 
     &                            SUM( Hg2aq(:,:,Total_Hg_Id) ), 'kg'
               CASE ( 3 )
                  ContainerName = 'Ocean_HgP' 
                  Temp2D        = HgPaq(:,:,Total_Hg_Id)
                  WRITE( 6, 260 ) ContainerName, 
     &                            SUM( HgPaq(:,:,Total_Hg_Id) ), 'kg'
            END SELECT
            Ptr2D => Temp2D

            ! Populate container with data
            CALL Diagn_Update( am_I_Root, HcoState, 
     &                         cName      = TRIM( ContainerName ),
     &                         COL        = GC_RST_Collection,
     &                         Array2D    = Ptr2D,
     &                         RC         = HCRC )

            ! Free the pointer
            Ptr2D => NULL()

            IF ( HCRC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot update restart species container: ' 
     &               // TRIM( ContainerName )
               CALL ERROR_STOP( MSG, LOC ) 
            ENDIF
         ENDDO

         !--------------------------------------------------------------
         ! Additional tagged Hg species in ocean
         !--------------------------------------------------------------
         IF ( Input_Opt%LSPLIT ) THEN
            DO M = 1, 3
               DO N = 2, Num_Hg_Categories

                  ! Define container name, point to data, and print 
                  ! total mass to log
                  Temp2D = 0.0d4
                  SELECT CASE( M )
                     CASE ( 1 )
                        HgSpc    = 'Hg0' 
                        Temp2D   = Hg0aq(:,:,N)
                     CASE ( 2 )
                        HgSpc    = 'Hg2'
                        Temp2D   = Hg2aq(:,:,N)
                     CASE ( 3 )
                        HgSpc    = 'HgP' 
                        Temp2D   = HgPaq(:,:,N)
                  END SELECT
                  Ptr2D => Temp2D
                  ContainerName = 'Ocean_' // Trim( HgSpc ) // 
     &                            '_'      // TRIM( Hg_Cat_Name(N) )

                  ! Populate container with data
                  CALL Diagn_Update( am_I_Root, HcoState, 
     &                               cName      = TRIM( ContainerName ),
     &                               COL        = GC_RST_Collection,
     &                               Array2D    = Ptr2D,
     &                               RC         = HCRC )
                  
                  ! Free the pointer
                  Ptr2D => NULL()
                  
                  ! Make sure tagged & total species sum up
                  IF ( Input_Opt%USE_CHECKS ) THEN
                     CALL CHECK_OCEAN_MERCURY(
     &                          'end of MAKE_OCEAN_Hg_RESTART' )
                  ENDIF
            
                  IF ( HCRC /= HCO_SUCCESS ) THEN
                     MSG = 'Cannot update restart species container: ' 
     &                     // TRIM( ContainerName )
                     CALL ERROR_STOP( MSG, LOC ) 
                  ENDIF

                  ! Print total mass to log
                  SELECT CASE( M )
                     CASE ( 1 )
                        WRITE( 6, 260 ) ContainerName, 
     &                            SUM( Hg0aq(:,:,N) ), 'kg'
                     CASE ( 2 )
                        WRITE( 6, 260 ) ContainerName, 
     &                            SUM( Hg2aq(:,:,N) ), 'kg'
                     CASE ( 3 )
                        WRITE( 6, 260 ) ContainerName, 
     &                            SUM( HgPaq(:,:,N) ), 'kg'
                  END SELECT

               ENDDO
            ENDDO
         ENDIF

         !--------------------------------------------------------------
         ! Hg snowpack on land and ocean
         !--------------------------------------------------------------
         DO M = 1, 4               
            DO N = 1, Num_Hg_Categories

               ! Define data, point to it, and print total mass to log
               Temp2D = 0.0d4
               SELECT CASE( M )
                  CASE ( 1 )
                     Prefix = 'Snow_Hg_Ocean_Hg0'  ! Reducible on ocean
                     Temp2D = SNOW_HG_OC(:,:,N)
                  CASE ( 2 )
                     Prefix = 'Snow_HgN_Ocean_Hg0' ! Non-reducible on ocean
                     Temp2D = SNOW_HG_STORED_OC(:,:,N)
                  CASE ( 3 )
                     Prefix = 'Snow_Hg_Land_Hg0'  ! Reducible on land
                     Temp2D = SNOW_HG_LN(:,:,N)
                  CASE ( 4 )
                     Prefix = 'Snow_HgN_Land_Hg0' ! Non-reducible on land
                     Temp2D = SNOW_HG_STORED_LN(:,:,N)
               END SELECT
               Ptr2D => Temp2D

               ! Define container name
               IF ( N == 1 ) THEN
                  ContainerName = TRIM( Prefix )
               ELSE
                  ! Append category name if tagged
                  ContainerName = TRIM( Prefix         ) // '_' //
     &                            TRIM( Hg_Cat_Name(N) )
               ENDIF
               
               ! Populate container with data
               CALL Diagn_Update( am_I_Root, HcoState,
     &                            cName      = TRIM( ContainerName ),
     &                            COL        = GC_RST_Collection,
     &                            Array2D    = Ptr2D,
     &                            RC         = HCRC )
               
               ! Free the pointer
               Ptr2D => NULL()
               
               IF ( HCRC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot update restart species container: ' 
     &                  // TRIM( ContainerName )
                  CALL ERROR_STOP( MSG, LOC ) 
               ENDIF

               ! Print total mass to log
               SELECT CASE( M )
                  CASE ( 1 )
                     WRITE( 6, 260 ) ContainerName, 
     &                         SUM( SNOW_HG_OC(:,:,N) ), 'kg'
                  CASE ( 2 )
                     WRITE( 6, 260 ) ContainerName, 
     &                         SUM( SNOW_HG_STORED_OC(:,:,N) ), 'kg'
                  CASE ( 3 )
                     WRITE( 6, 260 ) ContainerName, 
     &                         SUM( SNOW_HG_LN(:,:,N) ), 'kg'
                  CASE ( 4 )
                     WRITE( 6, 260 ) ContainerName, 
     &                         SUM( SNOW_HG_STORED_LN(:,:,N) ), 'kg'
               END SELECT

            ENDDO
         ENDDO

         ! Format string for printing Hg masses to log
 260     FORMAT( a24, ':   ', es15.9, 1x, a4)

         ! Free pointers for Hg indexing
         Hg_Cat_Name => NULL()

      ENDIF

      !=================================================================
      ! Write output NetCDF restart file using HEMCO 
      !=================================================================
      CALL HCOIO_Diagn_WriteOut( am_I_Root, HcoState, 
     &                           ForceWrite  = .TRUE.,
     &                           UsePrevTime = .FALSE.,
     &                           COL         = GC_RST_Collection,  
     &                           RC          = RC                )

      ! End restart write section of log file
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

      !### Debug
      IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN 
         CALL DEBUG_MSG('### WRITE_GC_RESTART_NC: wrote file')
      ENDIF

      END SUBROUTINE WRITE_GC_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_rst_dimensions
!
! !DESCRIPTION: Subroutine CHECK\_RST\_DIMENSIONS checks to see if dimensions 
!  read from the netCDF file match the defined GEOS-Chem dimensions.  If not, 
!  then it will stop the GEOS-Chem simulation with an error message.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Check_Rst_Dimensions( lon, lat, lev, time, 
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
!     CALL CHECK_RST_DIMENSION( lon=X,  lat=Y,           lev=Z,         &
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
!  01 Apr 2016 - E. Lundgren - Rename this subroutine from check_dimensions to
!                              check_rst_dimensions for global clarity
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Error message string
      CHARACTER(LEN=255) :: errMsg                  
      
      !=================================================================
      ! CHECK_RST_DIMENSIONS begins here!
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

      END SUBROUTINE Check_Rst_Dimensions
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compare_trc_spc
!
! !DESCRIPTION: Subroutine COMPARE\_TRC\_SPC prints out
!  species values to log in kg/kg dry for advected species, and
!  prints out non-advected species molec/cm3. This can only be used
!  where tracers were kg/kg dry.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COMPARE_TRC_SPC( am_I_Root, Input_Opt,
     &                            State_Met, State_Chm, LOC, RC ) 
!
! !USES:
!     
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,            INTENT(IN) :: am_I_Root  ! Are we on the root CPU?
      TYPE(OptInput),     INTENT(IN) :: Input_Opt   ! Input Options object
      TYPE(MetState),     INTENT(IN) :: State_Met   ! Meteorology State object
      CHARACTER(LEN=*),   INTENT(IN) :: LOC
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
!  08 Aug 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, L, N      ! lon, lat, lev, species
      REAL(fp)             :: MW_g            ! species molecular weight
      LOGICAL              :: INITIAL_SPC_MND

      ! Objects
      TYPE(Species),    POINTER :: SpcInfo

      !=================================================================
      ! COMPARE_TRC_SPC begins here!
      !=================================================================

      ! Assume success
      RC        = GC_SUCCESS

      ! Only print to log if in debug mode
      IF ( Input_Opt%LPRT ) THEN
      
         ! Write read message to log
         WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
         WRITE( 6, '(a,/)' ) '*** Compare Tracers to Species *** '
     &                       // '( tracer removal debugging )'
               
         ! Print header for min/max tracer and species [kg/kg dry]
         WRITE( 6, 910 )
 910     FORMAT( /,'Min and Max of advected [kg/kg dry]:' )
         WRITE( 6, 960 ) ADJUSTL( LOC )
         
         ! If species is in molec/cm3, convert to [kg/kg dry]
         INITIAL_SPC_MND = .FALSE.
         IF ( State_Chm%Spc_Units == 'molec/cm3' ) THEN
            CALL ConvertSpc_MND_to_KgKgDry( am_I_Root, State_Met,
     &                                      State_Chm, RC )
            INITIAL_SPC_MND = .TRUE.
         ENDIF
         
         ! Loop over advected species
         DO N = 1, State_Chm%nAdvect
         
            ! Get info about this species from the species database
            SpcInfo => State_Chm%SpcData(N)%Info
         
            ! Print the min & max of each species in kg/kg dry
            WRITE( 6, 920 ) N, TRIM( SpcInfo%Name ),
     &                      MINVAL( State_Chm%Species(:,:,:,N) ),
     &                      MAXVAL( State_Chm%Species(:,:,:,N) )
 920        FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9,
     &              '  Max = ',es15.9)
         
         ENDDO
         
         ! Convert species from kg/kg dry to molec/cm3 for next prints
         CALL ConvertSpc_KgKgDry_to_MND( am_I_Root, State_Met,
     &                                   State_Chm, RC )
         
         ! Print header for min/max of all species in  [molec/cm3]
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )      
         WRITE( 6, 940 )
 940     FORMAT( /,'Min and Max of all species [molec/cm3]:' )
         WRITE( 6, 960 ) ADJUSTL( LOC )
         
         ! Loop over all species
         DO N = 1, State_Chm%nSpecies
         
            ! Get info about this species from the species database
            SpcInfo => State_Chm%SpcData(N)%Info
         
            ! Print the min & max of each species in molec/cm3
            WRITE( 6, 950 ) N, TRIM( SpcInfo%Name ),
     &                      MINVAL( State_Chm%Species(:,:,:,N) ),
     &                      MAXVAL( State_Chm%Species(:,:,:,N) )
 950        FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9,
     &              '  Max = ',es15.9)
                  
         ENDDO
         WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

         ! If species were not initially molec/cm3, convert to kg/kg dry
         IF ( .NOT. INITIAL_SPC_MND ) THEN
            CALL ConvertSpc_MND_to_KgKgDry( am_I_Root, State_Met,
     &                                      State_Chm, RC )
         ENDIF

 960     FORMAT( /, 'Location: ', a60 )

      ENDIF

      END SUBROUTINE COMPARE_TRC_SPC
!EOC
      END MODULE RESTART_MOD
