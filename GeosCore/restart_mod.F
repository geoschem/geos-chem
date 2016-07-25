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
      PUBLIC  :: INIT_GC_RESTART
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CHECK_RST_DIMENSIONS
      PRIVATE :: CONVERT_TRACER_TO_VV
!
! !DEFINED PARAMETERS:
!
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
!  11 Jul 2016 - E. Lundgren - Store only species in restart files
!  12 Jul 2016 - E. Lundgren - Remove bpch restart file format option and
!                              associated routines
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
!  includes restart data that traditionally would be found in three different 
!  binary-punch format files: tracers, species, and ocean mercury.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GC_RESTART( am_I_Root, Input_Opt, 
     &                            State_Chm, RC         )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GRID_MOD,           ONLY : AREA_M2
      USE Species_Mod,        ONLY : Species
      USE TIME_MOD,           ONLY : GET_TS_CHEM

      USE HCO_Error_Mod
      USE HCO_Diagn_Mod
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
!                              NAMEGAS with ThisSpc%Name from species database
!  11 Jul 2016 - E. Lundgren - Remove tracer container; now only use species
!  12 Jul 2016 - E. Lundgren - Rename routine from init_gc_restart_nc
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
      REAL(fp), POINTER  :: GridAreas_ptr(:,:)
      CHARACTER(LEN=60)  :: ContainerName, Prefix, Tracer
      CHARACTER(LEN=255) :: MSG, LOC = 
     &                      'INIT_GC_RESTART (restart_mod.F90)'

      ! For Hg indexing 
      INTEGER                   :: N_Hg_CATS
      CHARACTER(LEN=4), POINTER :: Hg_Cat_Name(:) => NULL()

      ! Objects
      TYPE(Species), POINTER :: ThisSpc

      !=================================================================
      ! INIT_GC_RESTART begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Initialize pointer
      ThisSpc => NULL()

      ! Define collection variables
      GridAreas_ptr => AREA_M2(:,:,1)
      TS            =  GET_TS_CHEM() * 60.0_sp

      ! Create restart collection for GEOS-Chem.  This will keep the
      ! GEOS-Chem restart separate from other output.
      !
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

      !=================================================================
      ! Now handle indexing of Hg tracers locally.  These variables
      ! were in tracerid_mod.F, which will be removed (bmy, 4/26/16)
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Now get the # of tagged Hg categories from State_Chm
         N_Hg_CATS = State_Chm%N_Hg_CATS

         ! Now get Hg category names from State_Chm
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
         ThisSpc => State_Chm%SpcData(N)%Info

         ! Define container name and ID
         ContainerName = 'SPC_' // TRIM( ThisSpc%Name )
         cID = cID + 1

         ! Create container
         CALL Diagn_Create( am_I_Root,                     
     &                      Col       = GC_RST_Collection,         
     &                      cId       = cID,                 
     &                      cName     = TRIM( ContainerName ), 
     &                      SpaceDim  =  3,                
     &                      OutUnit   = 'mol mol-1',           
     &                      OutOper   = 'Instantaneous', 
     &                      RC        = RC )

         ! Free pointer
         ThisSpc => NULL()
         
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
                  Tracer = 'Hg0' 
               CASE ( 2 )
                  Tracer = 'Hg2'
               CASE ( 3 )
                  Tracer = 'HgP' 
            END SELECT
         
            ! Define container name and ID
            ContainerName = 'Ocean_' // TRIM( Tracer )
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
                  Prefix = 'Snow_Hg_Ocean_Hg0'  ! Reducible on ocean
               CASE ( 2 )
                  Prefix = 'Snow_HgN_Ocean_Hg0' ! Non-reducible on ocean
               CASE ( 3 )
                  Prefix = 'Snow_Hg_Land_Hg0'   ! Reducible on land
               CASE ( 4 )
                  Prefix = 'Snow_HgN_Land_Hg0'  ! Non-reducible on land
            END SELECT
         
            ! Loop over all categories to include tagged Hg
            DO N = 1, N_Hg_CATS

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
                     Tracer = 'Hg0' 
                  CASE ( 2 )
                     Tracer = 'Hg2'
                  CASE ( 3 )
                     Tracer = 'HgP' 
               END SELECT
            
               DO N = 2, N_Hg_CATS
               
                  ! Define container name and ID
                  ContainerName = 'Ocean_' // TRIM( TRACER         ) //
     &                            '_'      // TRIM( Hg_Cat_Name(N) )
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
      CHARACTER(LEN=255) :: INFILE    ! A
dvected tracer input  restart file
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

      ! For GEOS-Chem netcdf restart file, use format that mirrors
      ! HEMCO restart filename. NOTE: first part of this string is
      ! also hard-coded when creating GEOS-Chem restart collection.
      ! If one is changed, then the other must be changed as well.
      OUTPUT_RESTART_FILE = "GEOSChem_restart.YYYYMMDDhhmm"

      END SUBROUTINE SET_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gc_restart
!
! !DESCRIPTION: Subroutine READ\_GC\_RESTART initializes GEOS-Chem 
!  concentrations of species from a restart file in NetCDF format. Default
!  initial values are used for species missing from the restart file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_GC_RESTART( am_I_Root,  YYYYMMDD,  HHMMSS,
     &                            Input_Opt, State_Met, State_Chm, 
     &                            RC ) 
!
! !USES:
!     
      USE m_netcdf_io_open                         ! netCDF file open
      USE m_netcdf_io_read                         ! netCDF read
      USE m_netcdf_io_readattr                     ! netCDF attribute reads
      USE m_netcdf_io_get_dimlen                   ! netCDF get dimensions
      USE m_netcdf_io_close                        ! netCDF file close
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PHYSCONSTANTS,      ONLY : BOLTZ
      USE Species_Mod,        ONLY : Species
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
!                              NAMEGAS with ThisSpc%Name from species database
!  22 Jun 2016 - R. Yantosca - Now refer to Hg0_Id_List, Hg2_Id_List, and
!                              HgP_Id_List fields of State_Chm
!  11 Jul 2016 - E. Lundgren - Remove tracers and read only species
!  12 Jul 2016 - E. Lundgren - Rename from read_gc_restart_nc
!  25 Jul 2016 - E. Lundgren - Store whether species in rst file in species db
!                              rather than module-level variable
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, L, M, N      ! lon, lat, lev, indexes
      INTEGER              :: Tracer_id          ! ocean Hg tracer ids
      CHARACTER(LEN=10)    :: UNIT
      CHARACTER(LEN=60)    :: Prefix, Tracer
      REAL(fp)             :: Background_VV      ! background species conc
      REAL(fp)             :: CONV_FACTOR        ! convert mol/mol -> molec/cm3
      LOGICAL              :: FOUND              ! ocn Hg var in restart file

      !*******************************************************************
      ! TEMPORARY PENDING FAMILIES REMOVAL
      LOGICAL              :: FAMILIES_SPC       ! families species
      REAL(fp)             :: TRC_TO_SPC_FACTOR  ! tracer conc -> spc conc
      !*******************************************************************
      
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

      ! Other
      CHARACTER(LEN=255)   :: LOC

      ! For Hg species indexing (bmy, 4/26/16)
      INTEGER                   :: N_Hg_CATS 
      INTEGER                   :: Id_Hg_tot
      INTEGER,          POINTER :: Hg0_Id_List(:)
      INTEGER,          POINTER :: Hg2_Id_List(:)
      INTEGER,          POINTER :: HgP_Id_List(:)
      CHARACTER(LEN=4), POINTER :: Hg_Cat_Name(:)

      ! Objects
      TYPE(Species),    POINTER :: ThisSpc

      !=================================================================
      ! READ_GC_RESTART begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Initialize pointers
      Ptr2D       => NULL()
      Ptr3D       => NULL()
      Hg0_Id_List => NULL()
      Hg2_Id_List => NULL()
      HgP_Id_List => NULL()
      Hg_Cat_Name => NULL()
      ThisSpc     => NULL()

      ! Name of this routine
      LOC       = 'READ_GC_RESTART (restart_mod.F)'

      ! Copy input filename to a local variable
      nc_file   = TRIM( INPUT_RESTART_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

      !=================================================================
      ! Now handle indexing for Hg species locally (bmy, 4/26/16)
      ! These were formerly in tracerid_mod.F, which is being removed
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Now get the # of tagHg categories from State_Chm
         N_Hg_CATS   =  State_Chm%N_Hg_CATS

         ! Names for each of the Hg categories
         Hg_Cat_Name => State_Chm%Hg_Cat_Name
    
         ! Hg species index corresponding to a given Hg category number
         Hg0_Id_List => State_Chm%Hg0_Id_List
         Hg2_Id_List => State_Chm%Hg2_Id_List
         HgP_Id_List => State_Chm%HgP_Id_List

         ! ID_Hg_Tot is always the first category
         Id_Hg_Tot   =  Hg0_Id_List(1)

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
     &                       time_expected=1, LOC=LOC )

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
      ! Read species concentrations if present, else set to default vals
      !=================================================================

      ! Print header for min/max concentration to log
      WRITE( 6, 110 )
 110  FORMAT( /,'Min and Max of each species in file:' )

      ! Initialize species to all zeroes
      State_Chm%Species = 0.e+0_fp

      ! Initially set species units with the assumption that all species
      ! are in the restart file
      State_Chm%Spc_Units = 'molec/cm3'

      ! Loop over species
      DO N = 1, State_Chm%nSpecies

         ! Get info about this species from the species database
         ThisSpc => State_Chm%SpcData(N)%Info

         ! Define variable name
         v_name = 'SPC_' // TRIM( ThisSpc%Name )

         ! Check if species data is in file
         ierr = Nf_Inq_Varid (fId, v_name, v_id)
         IF ( ierr /= NF_NOERR ) THEN
            ThisSpc%Is_InRestart = .FALSE.
         ELSE
            ThisSpc%Is_InRestart = .TRUE.
         ENDIF

         !************************************************************
         ! TEMPORARY CODE PENDING FAMILIES REMOVAL (ewl, 7/11/16)
         ! If species is not in the restart file, check if is is part of 
         ! a family.,If yes, then the legacy family tracer should be in
         ! the initial restart file. Set the variable name to the family
         ! tracer name and set a factor for converting tracer concentration
         ! to species concentration (e.g. ISOPN * 0.5 for ISOPND and ISOPNB).
         IF ( .NOT. ThisSpc%Is_InRestart ) THEN

            SELECT CASE ( TRIM( ThisSpc%Name ) )
            
               CASE ('ISOPND', 'ISOPNB')
                  v_name = 'SPC_ISOPN'
                  TRC_TO_SPC_FACTOR = 0.5e+0_fp
            
               CASE ('MVKN', 'MACRN')
                  v_name = 'SPC_ISOPN'
                  TRC_TO_SPC_FACTOR = 0.5e+0_fp
            
               CASE( 'HCFC123', 'HCFC141b', 'HCFC142b' )
                  v_name = 'SPC_HCFCX'
                  TRC_TO_SPC_FACTOR = 1.0e+0_fp / 3.0e+0_fp
            
               CASE( 'CFC113', 'CFC114', 'CFC115' )
                  v_name = 'SPC_CFCX'
                  TRC_TO_SPC_FACTOR = 1.0e+0_fp / 3.0e+0_fp
            
               CASE DEFAULT
                  TRC_TO_SPC_FACTOR = 1.0e+0_fp
            
            END SELECT
   
            IF ( TRC_TO_SPC_FACTOR > 1.0e+0_fp ) THEN
               FAMILIES_SPC = .TRUE.
            ENDIF

         ENDIF
         !************************************************************

         ! If data in file, read and convert to molec/cm2
         IF ( ThisSpc%Is_InRestart .OR. FAMILIES_SPC ) THEN

            ! Initialize temporary array for this tracer and point to it
            ! use 1e-30 as min (hotp 2/25/09)
            ! NOTE: the 1.0e-30 minimum is for [molec/cm3] in legacy
            Temp3D = 1.0e-30_fp
            Ptr3D => Temp3D

            ! Read variable from file and store in local array
            CALL NcRd( Ptr3D, fId, TRIM(v_name), st4d, ct4d )

            !************************************************************
            ! TEMPORARY CODE PENDING TRACER REMOVAL
            ! For the first run, all species in the restart file are
            ! advected tracers and so State_Chm%Tracers and units must be set
            
            ! NOTE: Do not use this temporary code with restart files
            !       generated by a previous run.

            ! Print the min & max of each tracer as it is read from the file
            IF ( Input_Opt%LPRT ) THEN
               WRITE( 6, 120 ) N, TRIM( Input_Opt%TRACER_NAME(N) ),
     &                         MINVAL( Ptr3D ), MAXVAL( Ptr3D )
 120           FORMAT( 'Tracer ', i3, ', ', a8, ': Min = ', es15.9,
     &                 '  Max = ',es15.9)
            ENDIF
            
            ! Convert tracer concentration from [v/v dry] to [kg/kg dry] and 
            ! store in State_Chm%Tracers
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               State_Chm%Tracers(I,J,L,N) = Ptr3D(I,J,L) / 
     &                                      Input_Opt%TCVV(N) 
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
      
            ! Set State_Chm units
            State_Chm%Trac_Units = 'kg/kg dry'
            !************************************************************

            ! Convert [mol/mol] to [molec/cm3] to store in State_Met%Species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, CONV_FACTOR )
            ! Prior to 7/11/16: Loop over all chemistry grid boxes
            ! Tracer removal update: Now loop over all boxes
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
            
               ! Set box-dependent unit conversion factor 
               CONV_FACTOR = 
     &                 State_Met%PMID_DRY(I,J,L) * 1000e+0_fp /
     &                 ( State_Met%T(I,J,L) * BOLTZ * 1e+7_fp )
 
               ! Convert units from [mol/mol] to [molec/cm3/box]
               State_Chm%Species(I,J,L,N) = 
     &                 Ptr3D(I,J,L) * CONV_FACTOR
                                    
               !*********************************************************
               ! TEMPORARY CODE PENDING FAMILIES REMOVAL
               ! Adjust conversion factor to set family species
               ! concentration equal to the tracer concentration divided 
               ! by the number of species in the family. 
               IF ( FAMILIES_SPC ) THEN
                  State_Chm%Species(I,J,L,N) =
     &                    State_Chm%Species(I,J,L,N) * TRC_TO_SPC_FACTOR
               ENDIF
               !*********************************************************
          
               ! Make a small number if concentration is very low
               State_Chm%Species(I,J,L,N) = 
     &              MAX( State_Chm%Species(I,J,L,N), 1.0e-30_fp )
            
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Print the min and max of this species in [molec/cm3]
            IF ( Input_Opt%LPRT ) THEN
               WRITE( 6, 180 ) N, TRIM( ThisSpc%Name ),
     &                   MINVAL( State_Chm%Species(:,:,1:LLCHEM,N) ),
     &                   MAXVAL( State_Chm%Species(:,:,1:LLCHEM,N) )
 180           FORMAT('Species (molec/cm3) ', i3, ', ', a9, ': MIN = ', 
     &                 es15.9, ' MAX = ', es15.9)
            ENDIF

            ! Nullify pointer
            Ptr3D => NULL()

         ELSE

            !**********************************************************
            ! TRACER REMOVAL NOTE: (ewl, 7/11/16)
            ! Prior to tracer removal work ALL species were initialized
            ! using initial values for the first run of a simulation.
            ! With tracers removed, advected species will be initialized
            ! from the legacy tracer restart file values rather than the
            ! default initial values. This will impact output values.
            !**********************************************************

            ! If species data is not in the restart file, then assign
            ! default background values in units of [mol/mol].
            ! These values were traditionally stored in globchem.dat
            ! and assigned along with a unit conversion to [molec/cm3]
            ! in gasconc called from chemdr. Now init_flexchem converts
            ! these initial [mol/mol] values to [molec/cm3]. 
            ! Eventually store these values in the species database.
            SELECT CASE ( TRIM( ThisSpc%Name ) )
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
            IF ( ThisSpc%Name(1:3) == 'DRY') THEN
               Background_VV = 0.0e+0_fp
            ENDIF

            ! Put this array into State_Chm%Species [mol/mol], but
            ! only up to LLCHEM levels
            State_Chm%Species(:,:,1:LLCHEM,N) = Background_VV
            
            ! Print the background value of this species if debugging
            IF ( Input_Opt%LPRT ) THEN
               WRITE( 6, 190 ) N, TRIM( ThisSpc%Name ),
     &                         Background_VV
 190           FORMAT('Species ', i3, ', ', a9, ': Background = ', 
     &                 es15.9)
            ENDIF

            ! NOTE: Keep State_Chm%Species background values as
            ! mol/mol here to preserve location of unit conversion
            ! in pre-FlexChem code. This is temporary to minimize
            ! differences with FlexChem.

            ! Set State_Chm species units to reflect that conversion
            ! to molec/cm3 is not yet complete. This will prevent
            ! the incorrect application of unitconv_mod species unit 
            ! conversions prior to the full conversions to molec/cm3
            ! in init_flexchem (ewl, 7/22/16)
            State_Chm%Spc_Units = 'some species in v/v'

            ! Free pointer
            ThisSpc => NULL()

         ENDIF
      ENDDO

      !=================================================================
      ! Read ocean mercury variables
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Print total mass to log
         WRITE( 6, 200 )
 200     FORMAT(/, 'Total mass of each ocean and snow Hg tracer:')

         !--------------------------------------------------------------
         ! Total Hg in ocean
         !--------------------------------------------------------------
         DO M = 1, 3

            ! Define variable name
            SELECT CASE( M )
               CASE ( 1 )
                  Tracer = 'Hg0' 
                  Tracer_id = Hg0_Id_List(1)
               CASE ( 2 )
                  Tracer = 'Hg2'
                  Tracer_id = Hg2_Id_List(1)
               CASE ( 3 )
                  Tracer = 'HgP' 
                  Tracer_id = HgP_Id_List(1)
            END SELECT
            v_name = 'Ocean_' // TRIM( Tracer )

            ! Zero temporary array and set pointer
            Temp2D = 0.0d4 
            Ptr2D => Temp2D

            ! Check if variable is in file. If not found, initialize
            ! to zero (currently the case for tagged Hg) (ewl, 4/19/16)
            ierr = Nf_Inq_Varid (fId, TRIM(v_name), v_id)
            IF ( ierr /= NF_NOERR ) THEN
               FOUND = .FALSE.
               WRITE( 6, 210 ) TRIM( v_name )
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
                  Hg0aq(:,:,ID_Hg_tot) = Ptr2D
                  IF ( FOUND ) THEN
                     WRITE( 6, 220 ) v_name, 
     &                               SUM( Hg0aq(:,:,ID_Hg_tot) ), 'kg'
                  ENDIF
               CASE ( 2 )
                  Hg2aq(:,:,ID_Hg_tot) = Ptr2D
                  IF ( FOUND ) THEN
                     WRITE( 6, 220 ) v_name, 
     &                               SUM( Hg2aq(:,:,ID_Hg_tot) ), 'kg'
                  ENDIF
               CASE ( 3 )
                  HgPaq(:,:,ID_Hg_tot) = Ptr2D
                  IF ( FOUND ) THEN
                     WRITE( 6, 220 ) v_name, 
     &                               SUM( HgPaq(:,:,ID_Hg_tot) ), 'kg'
                  ENDIF
            END SELECT

            ! Nullify pointer
            Ptr2D => NULL()
         ENDDO

         !--------------------------------------------------------------
         ! Additional tagged ocean Hg tracers
         !--------------------------------------------------------------
         IF ( Input_Opt%LSPLIT ) THEN
            DO M = 1, 3
               DO N = 2, N_Hg_CATS

                  ! Define variable name. Include appended region.
                  SELECT CASE( M )
                     CASE ( 1 )
                        Tracer = 'Hg0' 
                        Tracer_id = Hg0_Id_List(N)
                     CASE ( 2 )
                        Tracer = 'Hg2'
                        Tracer_id = Hg2_Id_List(N)
                     CASE ( 3 )
                        Tracer = 'HgP' 
                        Tracer_id = HgP_Id_List(N)
                  END SELECT
                  v_name = 'Ocean_' // TRIM( Tracer         ) //
     &                     '_'      // TRIM( Hg_Cat_Name(N) )

                  ! Zero temporary array and set pointer
                  Temp2D = 0.0d4
                  Ptr2D => Temp2D
                  
                  ! Check if variable is in file. If not found, initialize
                  ! to zero (currently the case for tagged Hg) (ewl, 4/19/16)
                  ierr = Nf_Inq_Varid (fId, v_name, v_id)
                  IF ( ierr /= NF_NOERR ) THEN
                     FOUND = .FALSE.
                     WRITE( 6, 210 ) TRIM( v_name )
                  ELSE
                     FOUND = .TRUE.
                     CALL NcRd( Ptr2D, fId, v_name, st3d, ct3d )
                  ENDIF

                  ! Assign ocean mercury data and write total mass to log
                  SELECT CASE( M )
                     CASE ( 1 )
                        Hg0aq(:,:,N) = Ptr2D
                        IF ( FOUND ) THEN
                           WRITE( 6, 220 ) v_name, 
     &                                     SUM( Hg0aq(:,:,N) ), 'kg'
                        ENDIF
                     CASE ( 2 )
                        Hg2aq(:,:,N) = Ptr2D
                        IF ( FOUND ) THEN
                           WRITE( 6, 220 ) v_name, 
     &                                     SUM( Hg2aq(:,:,N) ), 'kg'
                        ENDIF
                     CASE ( 3 )
                        HgPaq(:,:,N) = Ptr2D
                        IF ( FOUND ) THEN
                           WRITE( 6, 220 ) v_name, 
     &                                     SUM( HgPaq(:,:,N) ), 'kg'
                        ENDIF
                  END SELECT

                  ! Nullify pointer
                  Ptr2D => NULL()
               ENDDO
            ENDDO

            ! Make sure tagged & total tracers sum up
            IF ( Input_Opt%USE_CHECKS ) THEN
               CALL CHECK_OCEAN_MERCURY( 'end of READ_GC_RESTART' )
            ENDIF
         ENDIF

         !--------------------------------------------------------------
         ! Hg snowpack on land and ocean
         !--------------------------------------------------------------
         DO M = 1, 4               
            DO N = 1, N_Hg_CATS               

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
                  WRITE( 6, 210 ) TRIM( v_name )
               ELSE
                  FOUND = .TRUE.
                  CALL NcRd( Ptr2D, fId, v_name, st3d, ct3d )
               ENDIF

               ! Assign ocean mercury data and write total mass to file
               SELECT CASE( M )
                  CASE ( 1 )
                     SNOW_HG_OC(:,:,N) = Ptr2D
                     IF ( FOUND ) THEN
                        WRITE( 6, 220 ) v_name, 
     &                        SUM( SNOW_HG_OC(:,:,N) ), 'kg'
                     ENDIF
                  CASE ( 2 )
                     SNOW_HG_STORED_OC(:,:,N) = Ptr2D
                     IF ( FOUND ) THEN
                        WRITE( 6, 220 ) v_name, 
     &                        SUM( SNOW_HG_STORED_OC(:,:,N) ), 'kg'
                     ENDIF
                  CASE ( 3 )
                     SNOW_HG_LN(:,:,N) = Ptr2D
                     IF ( FOUND ) THEN
                        WRITE( 6, 220 ) v_name, 
     &                     SUM( SNOW_HG_LN(:,:,N) ), 'kg'
                     ENDIF
                  CASE ( 4 )
                     SNOW_HG_STORED_LN(:,:,N) = Ptr2D
                     IF ( FOUND ) THEN
                        WRITE( 6, 220 ) v_name, 
     &                        SUM( SNOW_HG_STORED_LN(:,:,N) ), 'kg'
                     ENDIF
               END SELECT

               ! Nullify pointer
               Ptr2D => NULL()

            ENDDO
         ENDDO

         ! Format strings
 210     FORMAT('Variable not found in restart file: ', a24)
 220     FORMAT( a24, ':   ', es15.9, 1x, a4)

         ! Print note that variables are initialized to zero if not 
         ! found (currently only happens in tagged Hg simulation)
         IF ( Input_Opt%LSPLIT ) THEN
            WRITE( 6, 230 )
 230        FORMAT( /, 'NOTE: all variables not found are initialized ',
     &                 'to zero') 
         ENDIF

         ! Free pointers for Hg indexing
         Hg0_Id_List => NULL()
         Hg2_Id_List => NULL()
         HgP_Id_List => NULL()
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
! !DESCRIPTION: Subroutine WRITE\_GC\_RESTART creates restart files for 
!  GEOS-Chem advected tracers [units: vol tracer /vol dry air] in NetCDF
!  format.
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
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PHYSCONSTANTS,      ONLY : BOLTZ
      USE Species_Mod,        ONLY : Species
      USE TIME_MOD,           ONLY : EXPAND_DATE

      ! HEMCO
      USE HCO_STATE_MOD,      ONLY : HCO_STATE
      USE HCO_Error_Mod
      USE HCO_Diagn_Mod,      ONLY : Diagn_Update
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState
      USE HCOIO_Diagn_Mod,    ONLY : HCOIO_Diagn_WriteOut

      ! For Hg simulation
      USE OCEAN_MERCURY_MOD,  ONLY : Hg0aq, Hg2aq, HgPaq, Hgaq_tot 
      USE OCEAN_MERCURY_MOD,  ONLY : CHECK_OCEAN_MERCURY     
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_OC,  SNOW_HG_STORED_OC
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_LN,  SNOW_HG_STORED_LN
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      INTEGER,        INTENT(IN)  :: YYYYMMDD    ! YYYY/MM/DD GMT date
      INTEGER,        INTENT(IN)  :: HHMMSS      ! hh:mm:ss GMT time
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
!  20 Apr 2016 - E. Lundgren - Implement ocean and snow Hg variables
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with ThisSpc%Name from species database
!  22 Jun 2016 - R. Yantosca - Now refer to Hg0_Id_List, Hg2_Id_List, and
!                              HgP_Id_List fields of State_Chm
!  11 Jul 2016 - E. Lundgren - Remove tracers and write only species
!  12 Jul 2016 - E. Lundgren - Rename from write_gc_restart_nc
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local Variables      
      INTEGER              :: I, J, L, M, N, HCRC
      INTEGER              :: GC_RST_Collection
      INTEGER              :: Tracer_id    ! id for ocean Hg tracers
      REAL(fp)             :: CONV_FACTOR  ! convert molec/cm3 -> mol/mol
      REAL*4               :: SMALL_NUM
      CHARACTER(LEN=60)    :: ContainerName, Tracer, Prefix
      CHARACTER(LEN=255)   :: MSG, FILENAME, LOC = 
     &                       'WRITE_GC_RESTART (restart_mod.F90)'

      ! Temporary arrays and pointers
      REAL*4,            TARGET :: Temp3D( IIPAR, JJPAR, LLPAR )
      REAL*4,            TARGET :: Temp2D( IIPAR, JJPAR )
      REAL*4,           POINTER :: Ptr3D(:,:,:)
      REAL*4,           POINTER :: Ptr2D(:,:  )

      ! For indexing of Hg species
      INTEGER                   :: N_Hg_CATS
      INTEGER                   :: Id_Hg_tot
      INTEGER,          POINTER :: Hg0_Id_List(:)
      INTEGER,          POINTER :: Hg2_Id_List(:)
      INTEGER,          POINTER :: HgP_Id_List(:)
      CHARACTER(LEN=4), POINTER :: Hg_Cat_Name(:)

      ! Pointer to HEMCO state (is this needed?)
      TYPE(HCO_STATE),  POINTER :: HcoState

      ! Objects
      TYPE(Species),    POINTER :: ThisSpc

      !=================================================================
      ! WRITE_GC_RESTART begins here!
      !=================================================================

      ! Assume success
      RC          = GIGC_SUCCESS

      ! Initialize pointers
      Ptr3D       => NULL()
      Ptr2D       => NULL()
      Hg0_Id_List => NULL()
      Hg2_Id_List => NULL()
      HgP_Id_List => NULL()
      Hg_Cat_Name => NULL()
      HcoState    => NULL()
      ThisSpc     => NULL()

      ! Set the GEOS-Chem restart collection ID
      GC_RST_Collection = Input_Opt%GC_RST_COLLECTION

      ! Set minimum value of species concentration [molec/cm3/box]
      SMALL_NUM = 1.e-30_sp

      !=================================================================
      ! Initialize local Hg indexing variables (bmy, 4/26/16)
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN 
         
         ! Now get the # of tagHg categories from State_Chm
         N_Hg_CATS   =  State_Chm%N_Hg_CATS

         ! Names for each of the Hg categories
         Hg_Cat_Name => State_Chm%Hg_Cat_Name

         ! Hg species index corresponding to a given Hg category number
         Hg0_Id_List => State_Chm%Hg0_Id_List
         Hg2_Id_List => State_Chm%Hg2_Id_List
         HgP_Id_List => State_Chm%HgP_Id_List

         ! Total Hg is always the first category
         ID_Hg_Tot   =  Hg0_Id_List(1)

      ENDIF

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
 100  FORMAT( 'WRITE_GC_RESTART: Writing ', a )
      
      !=================================================================
      ! Populate species concentration restart containers
      !=================================================================

      ! Print min/max species conc to log in [molec/cm3]
      IF ( Input_Opt%LPRT ) THEN
         WRITE( 6, 130 )
 130     FORMAT( /, 'Min and Max of each species before unit ',
     &              'conversion [molec/cm3]:' )
      ENDIF

      ! Initialize Temp3D to all zeros
      Temp3D = 0.e+0_sp

      ! Populate Temp3D for all levels up to LLCHEM
      DO N = 1, State_Chm%nSpecies

         ! Get info about this species from the species database
         ThisSpc => State_Chm%SpcData(N)%Info

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLCHEM
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

         ! If debugging, print out species values up to LLCHEM.
         ! Above this level, all values are zero.
         IF ( Input_Opt%LPRT ) THEN
            WRITE( 6, 140 ) N, TRIM( ThisSpc%Name ),
     &                      MINVAL( Temp3D(:,:,1:LLCHEM) ), 
     &                      MAXVAL( Temp3D(:,:,1:LLCHEM) )
 140        FORMAT( 'Species ', i3, ', ', a9, ': Min = ', es15.9, 
     &              '  Max = ', es15.9 )
         ENDIF

         ! Free pointer
         ThisSpc => NULL()

      ENDDO

      ! Write header for post-conversion write to log (if debugging)
      IF ( Input_Opt%LPRT ) THEN
         WRITE( 6, 150 )
 150     FORMAT( /, 'Min and Max of each species written to file',
     &              ' [mol/mol]:' )
      ENDIF

      ! Convert species values locally to [mol/mol] for output
      ! in NetCDF. Print min and max values to log.
      DO N = 1, State_Chm%nSpecies
      
         ! Remove check for ND65 prod/loss families pending work on ND65 
         ! diagnostics and incorporation of FlexChem (ewl, 2/4/16)
         !! If the species is included in the IFAM array, which denotes 
         !! the index # of species in the CSPEC array that are "fake" ND65 
         !! prod/loss families, then skip. 
         !IF ( ANY( IFAM == N ) ) THEN
         !   CYCLE
         !ENDIF
      
         ! Get info about this species from the species database
         ThisSpc => State_Chm%SpcData(N)%Info

         ! Clear the species concentration array to all zeroes before
         ! populating array up to LLCHEM
         Temp3D(:,:,:) = 0.e+0_sp
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, CONV_FACTOR )
         DO L = 1, LLCHEM
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! NOTE: species is output to the binary-punch "spec_rst"
            ! restart file as molec/cm3 but is output to the NetCDF
            ! restart file as mol/mol. This results in small differences
            ! in the initial molec/cm3 species values at the beginning
            ! of the next run that uses the restart file. The difference
            ! is due to the pressure-dependency of the unit conversion 
            ! between molec/cm3 and mol/mol. 

            ! Set box-dependent unit conversion factor to convert
            ! from [molec/cm3] to [mol/mol]
            CONV_FACTOR = State_Met%T(I,J,L) * BOLTZ * 1e+7_fp /
     &                    ( State_Met%PMID_DRY(I,J,L) * 1000e+0_fp )

            ! Store species conc in [mol/mol] in local array
            IF ( State_Chm%Species(I,J,L,N) > SMALL_NUM ) THEN
               Temp3D(I,J,L) = State_Chm%Species(I,J,L,N) *
     &                         CONV_FACTOR
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
         ContainerName = 'SPC_' // TRIM( ThisSpc%Name )
         
         ! Populate container
         CALL Diagn_Update( am_I_Root,                   
     &                      cName   = TRIM( ContainerName ),
     &                      COL     = GC_RST_Collection,
     &                      Array3D = Ptr3D,             
     &                      RC      = HCRC )
         
         ! Free the pointer
         Ptr3D => NULL()

         ! Print the min & max of each species written to file
         IF ( Input_Opt%LPRT ) THEN
            WRITE( 6, 160 ) N, TRIM( ThisSpc%Name ),
     &                      MINVAL( Temp3D(:,:,1:LLCHEM) ), 
     &                      MAXVAL( Temp3D(:,:,1:LLCHEM) )
 160        FORMAT( 'Species ', i3, ', ', a9, ': Min = ', es15.9, 
     &                 '  Max = ', es15.9 )
         ENDIF
         
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot update restart species container: ' 
     &            // TRIM(ContainerName)
            CALL ERROR_STOP( MSG, LOC ) 
         ENDIF  

         ! Free pointer
         ThisSpc => NULL()

      ENDDO

      !=================================================================
      ! If Hg simulation, populate ocean mercury restart containers
      !=================================================================
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Print total mass to log
         WRITE( 6, 170 )
 170     FORMAT( /, 'Total mass of each ocean and snow Hg tracer: ' ) 

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
                  Tracer_id = Hg0_Id_List(1)
                  Temp2D = Hg0aq(:,:,ID_Hg_tot)
                  WRITE( 6, 180 ) ContainerName, 
     &                            SUM( Hg0aq(:,:,ID_Hg_tot) ), 'kg'
               CASE ( 2 )
                  ContainerName = 'Ocean_Hg2'
                  Tracer_id = Hg2_Id_List(1)
                  Temp2D = Hg2aq(:,:,ID_Hg_tot)
                  WRITE( 6, 180 ) ContainerName, 
     &                            SUM( Hg2aq(:,:,ID_Hg_tot) ), 'kg'
               CASE ( 3 )
                  ContainerName = 'Ocean_HgP' 
                  Tracer_id = HgP_Id_List(1)
                  Temp2D = HgPaq(:,:,ID_Hg_tot)
                  WRITE( 6, 180 ) ContainerName, 
     &                            SUM( HgPaq(:,:,ID_Hg_tot) ), 'kg'
            END SELECT
            Ptr2D => Temp2D

            ! Populate container with data
            CALL Diagn_Update( am_I_Root,                   
     &                      cName   = TRIM( ContainerName ),
     &                      COL     = GC_RST_Collection,
     &                      Array2D = Ptr2D,             
     &                      RC      = HCRC )

            ! Free the pointer
            Ptr2D => NULL()

            IF ( HCRC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot update restart tracer container: ' 
     &               // TRIM( ContainerName )
               CALL ERROR_STOP( MSG, LOC ) 
            ENDIF
         ENDDO

         !--------------------------------------------------------------
         ! Additional tagged Hg tracers in ocean
         !--------------------------------------------------------------
         IF ( Input_Opt%LSPLIT ) THEN
            DO M = 1, 3
               DO N = 2, N_Hg_CATS

                  ! Define container name, point to data, and print 
                  ! total mass to log
                  Temp2D = 0.0d4
                  SELECT CASE( M )
                     CASE ( 1 )
                        Tracer = 'Hg0' 
                        Tracer_id = Hg0_Id_List(N)
                        Temp2D = Hg0aq(:,:,N)
                     CASE ( 2 )
                        Tracer = 'Hg2'
                        Tracer_id = Hg2_Id_List(N)
                        Temp2D = Hg2aq(:,:,N)
                     CASE ( 3 )
                        Tracer = 'HgP' 
                        Tracer_id = HgP_Id_List(N)
                        Temp2D = HgPaq(:,:,N)
                  END SELECT
                  Ptr2D => Temp2D
                  ContainerName = 'Ocean_' // Trim( Tracer         ) // 
     &                            '_'      // TRIM( Hg_Cat_Name(N) )

                  ! Populate container with data
                  CALL Diagn_Update( am_I_Root,                   
     &                               cName   = TRIM( ContainerName ),
     &                               COL     = GC_RST_Collection,
     &                               Array2D = Ptr2D,             
     &                               RC      = HCRC )
                  
                  ! Free the pointer
                  Ptr2D => NULL()
                  
                  ! Make sure tagged & total tracers sum up
                  IF ( Input_Opt%USE_CHECKS ) THEN
                     CALL CHECK_OCEAN_MERCURY(
     &                          'end of MAKE_OCEAN_Hg_RESTART' )
                  ENDIF
            
                  IF ( HCRC /= HCO_SUCCESS ) THEN
                     MSG = 'Cannot update restart tracer container: ' 
     &                     // TRIM( ContainerName )
                     CALL ERROR_STOP( MSG, LOC ) 
                  ENDIF

                  ! Print total mass to log
                  SELECT CASE( M )
                     CASE ( 1 )
                        WRITE( 6, 180 ) ContainerName, 
     &                            SUM( Hg0aq(:,:,N) ), 'kg'
                     CASE ( 2 )
                        WRITE( 6, 180 ) ContainerName, 
     &                            SUM( Hg2aq(:,:,N) ), 'kg'
                     CASE ( 3 )
                        WRITE( 6, 180 ) ContainerName, 
     &                            SUM( HgPaq(:,:,N) ), 'kg'
                  END SELECT

               ENDDO
            ENDDO
         ENDIF

         !--------------------------------------------------------------
         ! Hg snowpack on land and ocean
         !--------------------------------------------------------------
         DO M = 1, 4               
            DO N = 1, N_Hg_CATS

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
               CALL Diagn_Update( am_I_Root,                   
     &                            cName   = TRIM( ContainerName ),
     &                            COL     = GC_RST_Collection,
     &                            Array2D = Ptr2D,             
     &                            RC      = HCRC )
               
               ! Free the pointer
               Ptr2D => NULL()
               
               IF ( HCRC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot update restart tracer container: ' 
     &                  // TRIM( ContainerName )
                  CALL ERROR_STOP( MSG, LOC ) 
               ENDIF

               ! Print total mass to log
               SELECT CASE( M )
                  CASE ( 1 )
                     WRITE( 6, 180 ) ContainerName, 
     &                         SUM( SNOW_HG_OC(:,:,N) ), 'kg'
                  CASE ( 2 )
                     WRITE( 6, 180 ) ContainerName, 
     &                         SUM( SNOW_HG_STORED_OC(:,:,N) ), 'kg'
                  CASE ( 3 )
                     WRITE( 6, 180 ) ContainerName, 
     &                         SUM( SNOW_HG_LN(:,:,N) ), 'kg'
                  CASE ( 4 )
                     WRITE( 6, 180 ) ContainerName, 
     &                         SUM( SNOW_HG_STORED_LN(:,:,N) ), 'kg'
               END SELECT

            ENDDO
         ENDDO

         ! Format string for printing Hg masses to log
 180     FORMAT( a24, ':   ', es15.9, 1x, a4)

         ! Free pointers for Hg indexing
         Hg_Cat_Name => NULL()
         Hg0_Id_List => NULL()
         Hg2_Id_List => NULL() 
         HgP_Id_List => NULL() 

      ENDIF

      !=================================================================
      ! Write output NetCDF restart file using HEMCO 
      !=================================================================

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

      ! End restart write section of log file
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

      !### Debug
      IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN 
         CALL DEBUG_MSG('### WRITE_GC_RESTART: wrote file')
      ENDIF

      END SUBROUTINE WRITE_GC_RESTART
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
      END MODULE RESTART_MOD
