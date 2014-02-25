!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gamap_mod
!
! !DESCRIPTION: Module GAMAP\_MOD contains routines to create GAMAP 
!  "tracerinfo.dat" and "diaginfo.dat" files which are customized to each
!   particular GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE GAMAP_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD                    ! Dimensions of arrays
      USE CMN_DIAG_MOD                    ! Diagnostic parameters
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DO_GAMAP
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: CREATE_DINFO
      PRIVATE :: CREATE_TINFO
      PRIVATE :: WRITE_TINFO
      PRIVATE :: WRITE_SEPARATOR
      PRIVATE :: INIT_DIAGINFO
      PRIVATE :: INIT_TRACERINFO
      PRIVATE :: INIT_GAMAP
      PRIVATE :: CLEANUP_GAMAP
!
! !REMARKS:
!  For more information, please see the GAMAP Online Users' Manual:
!   http://acmg.seas.harvard.edu/gamap/doc/index.html
!
! !REVISION HISTORY: 
!  03 May 2005 - R. Yantosca - Initial version
!  (1 ) Minor bug fix for Rn/Pb/Be simulations (bmy, 5/11/05)
!  (2 ) Added ND09 diagnostic for HCN/CH3CN simulation. (bmy, 6/30/05)
!  (3 ) Added ND04 diagnostic for CO2 simulation (bmy, 7/25/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Add MBO to ND46 diagnostic (tmf, bmy, 10/20/05)
!  (6 ) Updated for tagged Hg simulation (cdh, bmy, 4/6/06)
!  (7 ) Updated for ND56 lightning flash diagnostics (ltm, bmy, 5/5/06)
!  (8 ) Updated for ND42 SOA concentration diagnostics (dkh, bmy, 5/22/06)
!  (9 ) Updated for ND36 CH3I simulation diagnostics (bmy, 7/25/06)
!  (10) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (11) Add routines INIT_DIAGINFO, INIT_TRACERINFO for clarity.  Added new
!        entries for biomass burning (ND28) and time in tropopshere (ND54)
!        in INIT_DIAGINFO and INIT_TRACERINFO. (phs, bmy, 10/17/06)
!  (12) Now write GPROD & APROD info to diaginfo.dat, tracerinfo.dat files,
!        for the SOA restart files (tmf, havala, bmy, 2/6/07)
!  (13) Added ND10 diagnostic for H2/HD simulation. (phs, 9/18/07)
!  (14) Change category name for ND31 diagnostic (bmy, 11/16/07)
!  (15) Add to tracerinfo.dat file for timeseries and Rn-Pb-Be (bmy, 2/22/08)
!  (16) Added ND52 diagnostic for gamma HO2 (jaegle 02/26/09)
!  (17) Add gamap info for dicarbonyl simulation (tmf, 3/10/09)
!  (18) Add C2H4 in ND46 (ccc, 3/10/09)
!  (19) Add EFLUX to ND67 (lin, ccc, 5/29/09)
!  (20) Minor bug fixes (bmy, phs, 10/9/09)
!  (20) Minor bug fixes (dkh, bmy, 11/19/09)
!  (21) Include second satellite overpass diagnostic.  Adjust AOD name to 550 
!        nm from 400 nm.  Add additional dust AOD bins.  Output values to 
!        hdf_mod. (amv, bmy, 12/1/09)
!  (22) Increase MAXTRACER from 120 to 325 (win, 6/25/09) 
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!  03 Aug 2010 - R. Yantosca - Now move the #include "CMN_SIZE" and
!                              #include "CMN_DIAG" to the top of module
!  13 Aug 2010 - R. Yantosca - Added modifications for MERRA
!  21 Sep 2010 - R. Yantosca - Removed duplicates in INIT_DIAGINFO
!  21 Oct 2010 - R. Yantosca - Bug fix in INIT_DIAGINFO
!  09 Dec 2010 - C. Carouge  - Modify MAXTRACER definition to account for 
!  07 Feb 2012 - E. Corbitt  - Added diagnostic info for tagged Hg simulation.
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA +
!                              semivolatile POA simulations (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! For "diaginfo.dat"
      INTEGER,           PARAMETER   :: MAXCAT    = 150
      INTEGER,           PARAMETER   :: SPACING   = 1000
      INTEGER                        :: NCATS
      INTEGER,           ALLOCATABLE :: OFFSET(:)
      CHARACTER(LEN=40), ALLOCATABLE :: CATEGORY(:)
      CHARACTER(LEN=40), ALLOCATABLE :: DESCRIPT(:)
      CHARACTER(LEN=255)             :: DFILE

      ! For "tracerinfo.dat"
      INTEGER,           PARAMETER   :: MAXDIAG   = MAX_DIAG
      INTEGER,           PARAMETER   :: MAXTRACER = 2 * MAX_TRACER
      INTEGER,           ALLOCATABLE :: NTRAC(:)
      INTEGER,           ALLOCATABLE :: INDEX(:,:)
      INTEGER,           ALLOCATABLE :: MOLC(:,:)
      REAL*4,            ALLOCATABLE :: MWT(:,:)
      REAL*4,            ALLOCATABLE :: SCALE(:,:)
      CHARACTER(LEN=40), ALLOCATABLE :: NAME(:,:)
      CHARACTER(LEN=40), ALLOCATABLE :: FNAME(:,:)
      CHARACTER(LEN=40), ALLOCATABLE :: UNIT(:,:)
      CHARACTER(LEN=255)             :: TFILE
      
      ! Other variables
      CHARACTER(LEN=16)              :: STAMP
      CHARACTER(LEN=40)              :: SIM_NAME

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_gamap
!
! !DESCRIPTION: Subroutine DO\_GAMAP is the driver program for creating 
!  the customized GAMAP files "diaginfo.dat" and "tracerinfo.dat". 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_GAMAP( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : SYSTEM_TIMESTAMP
      USE TRACER_MOD,         ONLY : GET_SIM_NAME 
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  03 May 2005 - R. Yantosca - Initial version
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! DO_GAMAP begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate and initialize variables
      CALL INIT_GAMAP( am_I_Root, Input_Opt, RC )

      ! Create a timestamp with the system date & time
      STAMP    = SYSTEM_TIMESTAMP()

      ! Get simulation name
      SIM_NAME = GET_SIM_NAME()

      ! Write "diaginfo.dat" file
      CALL CREATE_DINFO()

      ! Write "tracerinfo.dat" file
      CALL CREATE_TINFO( am_I_Root, Input_Opt, RC )

      ! Deallocate variables
      CALL CLEANUP_GAMAP

      END SUBROUTINE DO_GAMAP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: create_dinfo
!
! !DESCRIPTION: Subroutine CREATE\_DINFO writes information about diagnostic 
!  categories to a customized "diaginfo.dat" file. (bmy, 5/3/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CREATE_DINFO
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
! 
! !REVISION HISTORY: 
!  03 May 2005 - R. Yantosca - Initial version
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: IOS, N
      INTEGER :: IU_FILE

      !=================================================================
      ! CREATE_DINFO begins here!
      !=================================================================

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open "diaginfo.dat" file for output
      OPEN( IU_FILE, FILE=TRIM( DFILE ), STATUS='UNKNOWN', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'create_dinfo:1' )

      ! Write file header
      WRITE( IU_FILE, '(a)' ) '#' // REPEAT( '=', 78 )
      WRITE( IU_FILE,  100  ) STAMP
      WRITE( IU_FILE,  105  ) TRIM( SIM_NAME )
      WRITE( IU_FILE,  110  )
      WRITE( IU_FILE,  115  )
      WRITE( IU_FILE, '(a)' ) '# File Format:'
      WRITE( IU_FILE, '(a)' ) '# ' // REPEAT( '-', 77 )
      WRITE( IU_FILE,  120  )
      WRITE( IU_FILE,  125  )
      WRITE( IU_FILE,  130  )
      WRITE( IU_FILE,  135  )
      WRITE( IU_FILE,  125  )
      WRITE( IU_FILE, '(a)' ) '#'
      WRITE( IU_FILE,  140  ) SPACING
      WRITE( IU_FILE, '(a)' ) '#' // REPEAT( '=', 78 )

      ! Loop over categories
      DO N = 1, NCATS

         ! Write one line to "diaginfo.dat" file
         WRITE( IU_FILE, 145, IOSTAT=IOS ) 
     &        OFFSET(N), ADJUSTL( CATEGORY(N) ), ADJUSTL( DESCRIPT(N) )

         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'create_dinfo:1' )
      ENDDO

      ! FORMAT strings
 100  FORMAT( '# diaginfo.dat: Created by GEOS-CHEM at ', a,     /,'#' )
 105  FORMAT( '# ****** CUSTOMIZED FOR ', a, ' SIMULATION *****',/,'#' )
 110  FORMAT( '# This file contains category names and the offsets',
     &        ' which they are stored'                                 ) 
 115  FORMAT( '# in file "tracerinfo.dat".  This file is read into',
     &        ' GAMAP by routine', /, '# "ctm_diaginfo.pro".',   /,'#' )
 120  FORMAT( '# OFFSET    (I8 )  Constant to add to tracer numbers', 
     &        ' in order to distinguish', /, '#', 18x, 'for the given',
     &        ' diagnostic category, as stored in file', /, '#', 18x, 
     &        '"tracerinfo.dat".  OFFSET may be up to 8 digits long.'  )
 125  FORMAT( '#  --       (1X )  1-character spacer'                  )    
 130  FORMAT( '# CATEGORY  (A40)  Category name for CTM diagnostics.', 
     &        ' NOTE: The category name', /, '#', 18x, 
     &        'can be up to 40 chars long, but historically the', 
     &        ' GEOS-CHEM', /,'#', 18x, 'and GISS models have used an', 
     &        ' 8-character category name.'                            )
 135  FORMAT( '# COMMENT   (A  )  Descriptive comment string',   /,'#' )
 140  FORMAT('##### SPACING BETWEEN DIAGNOSTIC CATEGORY OFFSETS = ',i8 )
 145  FORMAT( i8, 1x, a40, 1x, a )

      ! Close file
      CLOSE( IU_FILE )

      END SUBROUTINE CREATE_DINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: create_tinfo
!
! !DESCRIPTION: Subroutine CREATE\_TINFO writes information about tracers to 
!  a customized tracerinfo.dat" file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CREATE_TINFO( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE FILE_MOD,    ONLY : IOERROR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  21 Apr 2005 - R. Yantosca - Initial version
!  (1 ) Now write out tracers in ug/m3 (dkh, bmy, 5/22/06)
!  (2 ) Now write out GPROD & APROD info (tmf, havala, bmy, 2/6/07)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  31 Oct 2012 - R. Yantosca - Now save out info about soil NOx restart file
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: D, IOS, N, T
      INTEGER           :: IU_FILE
      REAL*4            :: SCALE_NEW
      CHARACTER(LEN=2)  :: C
      CHARACTER(LEN=40) :: UNIT_NEW, NAME_NEW

      !=================================================================
      ! CREATE_TINFO begins here!
      !=================================================================

      ! Assume success
      RC      =  GIGC_SUCCESS

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open "tracerinfo.dat" file for output
      OPEN( IU_FILE, FILE=TRIM( TFILE ), STATUS='UNKNOWN', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'create_tinfo:1' )

      ! Write file header
      WRITE( IU_FILE, '(a)' ) '#' // REPEAT( '=', 78 )
      WRITE( IU_FILE,  100  ) STAMP
      WRITE( IU_FILE,  105  ) TRIM( SIM_NAME )
      WRITE( IU_FILE,  110  ) 
      WRITE( IU_FILE,  115  )
      WRITE( IU_FILE, '(a)' ) '# File Format:' 
      WRITE( IU_FILE, '(a)' ) '# ' // REPEAT( '-', 77 )    
      WRITE( IU_FILE,  120  )
      WRITE( IU_FILE,  125  )
      WRITE( IU_FILE,  130  )
      WRITE( IU_FILE,  135  ) 
      WRITE( IU_FILE,  140  )
      WRITE( IU_FILE,  145  )
      WRITE( IU_FILE,  150  )
      WRITE( IU_FILE,  125  )
      WRITE( IU_FILE,  155  )

      ! FORMAT strings
 100  FORMAT( '# tracerinfo.dat: Created by GEOS-CHEM at ', a,   /,'#' )
 105  FORMAT( '# ****** CUSTOMIZED FOR ', a, ' SIMULATION *****',/,'#' )
 110  FORMAT( '# This file contains name weight and index', 
     &        ' information about GEOS-CHEM'                           )
 115  FORMAT( '# tracers.  It is read by routine ',
     &        '"ctm_tracerinfo.pro" of the GAMAP package.',      /,'#' )
 120  FORMAT( '# NAME     (A8   )  Tracer name (up to 8 chars)'        )
 125  FORMAT( '#  --      (1X   )  1-character spacer'                 )
 130  FORMAT( '# FULLNAME (A30  )  Full tracer name (up to 30 chars)'  )
 135  FORMAT( '# MOLWT    (E10.0)  Molecular weight (kg/mole)'         )
 140  FORMAT( '# C        (I3   )  For HC''s: # moles C/moles tracer;',
     &        ' otherwise set=1'                                       )
 145  FORMAT( '# TRACER   (I9   )  Tracer number (up to 9 digits)'     )
 150  FORMAT( '# SCALE    (E10.3)  Standard scale factor to convert', 
     &        ' to unit given below'                                   )
 155  FORMAT( '# UNIT     (A40  )  Unit string',                 /,'#' )

      !-------------------------------------
      ! 0: Tracers [ppbv]
      !-------------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( IU_FILE, 0 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)

         ! GAMAP tracer number
         N = ( SPACING * 0 ) + T

         ! Write tracers [ppbv] to "tracerinfo.dat" file
         CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE(T,45), UNIT(T,45), N )
      ENDDO

      !-------------------------------------
      ! SPACING*1: Tracers [molec/cm2/s]
      !-------------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( IU_FILE, 100 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)
         
         ! GAMAP tracer number
         N = ( SPACING * 1 ) + T

         ! New scale
         SCALE_NEW = 1.0e0
         
         ! New unit
         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
            UNIT_NEW = 'atoms C/cm2/s'
         ELSE
            UNIT_NEW = 'molec/cm2/s'
         ENDIF

         ! Write tracers [molec/cm2/s] to "tracerinfo.dat"
         CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
      ENDDO

      !-------------------------------------
      ! SPACING*2: Tracers [molec/cm2]
      !-------------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( IU_FILE, 200 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)
         
         ! GAMAP tracer number
         N = ( SPACING * 2 ) + T

         ! New scale
         SCALE_NEW = 1.0e0
         
         ! New unit
         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
            UNIT_NEW = 'atoms C/cm2'
         ELSE
            UNIT_NEW = 'molec/cm2'
         ENDIF

         ! Write tracers [molec/cm2] to "tracerinfo.dat"
         CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
      ENDDO

      !-------------------------------------
      ! SPACING*3: Tracers [kg/s]
      !-------------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( IU_FILE, 300 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)

         ! GAMAP tracer number
         N = ( SPACING * 3 ) + T

         ! New scale
         SCALE_NEW = 1.0e0

         ! New unit
         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
            UNIT_NEW = 'kg C/s'
         ELSE
            UNIT_NEW = 'kg/s'
         ENDIF

         ! Write tracers [kg/s] to "tracerinfo.dat"
         CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
      ENDDO

      !-------------------------------------
      ! SPACING*4: Tracers [kg]
      !-------------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( IU_FILE, 400 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)

         ! GAMAP tracer number
         N = ( SPACING * 4 ) + T

         ! New scale
         SCALE_NEW = 1.0e0

         ! New unit
         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
            UNIT_NEW = 'kg C'
         ELSE
            UNIT_NEW = 'kg'
         ENDIF

         ! Write tracers [kg] to "tracerinfo.dat"
         CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
      ENDDO

! Add this later on (bmy, 5/22/06)
!      !-----------------------------------
!      ! SPACING*5: Tracers [ug/m3]
!      !-----------------------------------
!
!      ! Write separator line
!      CALL WRITE_SEPARATOR( IU_FILE, 500 )
!
!      ! Loop over tracers
!      DO T = 1, NTRAC(45)
!
!         ! GAMAP tracer number
!         N = ( SPACING * 5 ) + T
!
!         ! New scale
!         SCALE_NEW = 1.0e0
!
!         ! New unit
!         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
!            UNIT_NEW = 'ug C/m3'
!         ELSE
!            UNIT_NEW = 'ug/m3'
!         ENDIF
!
!         ! Write tracers [kg] to "tracerinfo.dat"
!         CALL WRITE_TINFO( IU_FILE, NAME(T,45), FNAME(T,45), MWT(T,45), 
!     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
!      ENDDO

!     !-------------------------------------
!     ! SPACING*6: GPROD & APROD [kg/kg]
!     !-------------------------------------
!     IF ( Input_Opt%LSOA ) THEN
!
!        ! Write separator line
!        CALL WRITE_SEPARATOR( IU_FILE, 600 )
!
!        ! Loop over tracers
!        ! New tracers for SOA simulation
!--- Previous (ccc, 3/17/10)
!         DO T = 1, 18
!        DO T = 1, 28
!           
!           ! GAMAP tracer number
!           N = ( SPACING * 6 ) + T
!           
!           ! Make a character 
!           WRITE( C, '(i2.2)' ) T
!
!           ! Tracer name
!           NAME_NEW = 'PROD' // C
!           
!           ! Write tracers [kg] to "tracerinfo.dat"
!           CALL WRITE_TINFO( IU_FILE, NAME_NEW, NAME_NEW,
!    &           REAL(1e0,4), 1, REAL(1e0,4), 'kg/kg', N )
!        ENDDO
!
!     ENDIF

      !-------------------------------------      
      ! SPACING*60: Soil NOx restart file
      !-------------------------------------   
      IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
         
         ! Write separator line
         CALL WRITE_SEPARATOR( IU_FILE, 999 )

         ! Dry period [h]
         NAME_NEW = 'P_DRY'
         UNIT_NEW = 'hours'
         CALL WRITE_TINFO(  IU_FILE,     NAME_NEW, 
     &                      NAME_NEW,    REAL(1e0,4), 1,            
     &                      REAL(1e0,4), UNIT_NEW,    60001  )

         ! Pulse factor [unitless]
         NAME_NEW = 'PFACTOR'
         UNIT_NEW = 'unitless'
         CALL WRITE_TINFO(  IU_FILE,     NAME_NEW, 
     &                      NAME_NEW,    REAL(1e0,4), 1,            
     &                      REAL(1e0,4), UNIT_NEW,    60002  )

         ! GWET from previous timestep [unitless]
         NAME_NEW = 'GWET_PRV'
         UNIT_NEW = 'unitless'
         CALL WRITE_TINFO(  IU_FILE,     NAME_NEW, 
     &                      NAME_NEW,    REAL(1e0,4), 1,            
     &                      REAL(1e0,4), UNIT_NEW,    60003  )

         ! Reservoir of deposited N [ng N/m2]
         NAME_NEW = 'N_RESVR'
         UNIT_NEW = 'ng N/m2'
         CALL WRITE_TINFO(  IU_FILE,     NAME_NEW, 
     &                      NAME_NEW,    REAL(1e0,4), 1,            
     &                      REAL(1e0,4), UNIT_NEW,    60004  )

      ENDIF

      !-------------------------------------
      ! All other diagnostics 
      !-------------------------------------
      DO D = 1, MAXDIAG

         ! If tracers are defined then...
         IF ( NTRAC(D) > 0 ) THEN

            ! Skip ND45, we already wrote tracers above
            IF ( D == 45 ) CYCLE

            ! Write separator
            CALL WRITE_SEPARATOR( IU_FILE, D )

            ! Write tracers to file
            DO T = 1, NTRAC(D)
            CALL WRITE_TINFO( IU_FILE, NAME(T,D), FNAME(T,D), MWT(T,D), 
     &                           MOLC(T,D), SCALE(T,D), UNIT(T,D),
     &                           INDEX(T,D) )
            ENDDO
         ENDIF
      ENDDO

      END SUBROUTINE CREATE_TINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_tinfo
!
! !DESCRIPTION: Subroutine WRITE\_TINFO writes one line to the customized 
!  "tracerinfo.dat" file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_TINFO( IU_FILE, NAME, FNAME, 
     &                        MWT, MOLC, SCALE, UNIT, N )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN) :: IU_FILE ! Logical unit number
      CHARACTER(LEN=*), INTENT(IN) :: NAME    ! GAMAP short tracer name
      CHARACTER(LEN=*), INTENT(IN) :: FNAME   ! GAMAP long tracer name
      REAL*4,           INTENT(IN) :: MWT     ! Molecular weight [kg/mole]
      INTEGER,          INTENT(IN) :: MOLC    ! Moles C/mole tracer (for HC's)
      INTEGER,          INTENT(IN) :: N       ! Tracer number
      REAL*4,           INTENT(IN) :: SCALE   ! GAMAP scale factor
      CHARACTER(LEN=*), INTENT(IN) :: UNIT    ! Unit string
! 
! !REVISION HISTORY: 
!  03 May 2005 - R. Yantosca - Initial version
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: IOS

      !=================================================================
      ! WRITE_TINFO begins here!
      !=================================================================

      ! Write one line to "tracerinfo.dat" file
      WRITE( IU_FILE, 100, IOSTAT=IOS ) 
     &     ADJUSTL( NAME ), ADJUSTL( FNAME ), MWT, 
     &     MOLC,            N,                SCALE, TRIM( UNIT )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_tinfo:1' )

      ! FORMAT string
 100  FORMAT( a8, 1x, a30, es10.3, i3, i9, es10.3, 1x, a )

      END SUBROUTINE WRITE_TINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_separator

!
! !DESCRIPTION: Subroutine WRITE\_SEPARATOR writes a separator block to the 
!  customized "tracerinfo.dat" file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_SEPARATOR( IU_FILE, DIAG )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: IU_FILE
      INTEGER, INTENT(IN) :: DIAG   ! GEOS-Chem diagnostic number
! 
! !REVISION HISTORY: 
!  03 May 2005 - R. Yantosca - Initial version
!  06 Feb 2007 - R. Yantosca - Added new header for GPROD & APROD info
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!  31 Oct 2012 - R. Yantosca - Add write statement for soil NOx restart file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: IOS
      CHARACTER(LEN=79) :: SEPARATOR

      !=================================================================
      ! WRITE_SEPARATOR begins here! 
      !=================================================================

      ! Create separator string
      SEPARATOR = '#' // REPEAT( '=', 78 )

      ! Write separator string
      WRITE( IU_FILE, '(a)', IOSTAT=IOS ) SEPARATOR
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_separator:1' )

      ! Write the appropriate message
      SELECT CASE( DIAG )
         CASE( 0 ) 
            WRITE( IU_FILE, 100, IOSTAT=IOS )
         CASE( 100 ) 
            WRITE( IU_FILE, 110, IOSTAT=IOS )
         CASE( 200 ) 
            WRITE( IU_FILE, 120, IOSTAT=IOS )
         CASE( 300 ) 
            WRITE( IU_FILE, 130, IOSTAT=IOS )
         CASE( 400 ) 
            WRITE( IU_FILE, 140, IOSTAT=IOS )
         CASE( 500 ) 
            WRITE( IU_FILE, 150, IOSTAT=IOS )
         CASE( 600 ) 
            WRITE( IU_FILE, 160, IOSTAT=IOS )
         CASE( 700 )          !(win, 7/9/09)
            WRITE( IU_FILE, 190, IOSTAT=IOS ) !sfarina - this can stay
         CASE( 999 ) 
            WRITE( IU_FILE, 170, IOSTAT=IOS )
         CASE DEFAULT
            WRITE( IU_FILE, 180, IOSTAT=IOS ) DIAG
      END SELECT

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_separator:2' )

      ! Write separator string
      WRITE( IU_FILE, '(a)', IOSTAT=IOS ) SEPARATOR
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_separator:3' )

      ! FORMAT strings
 100  FORMAT( '# GEOS-CHEM tracers [ppbv]'             )
 110  FORMAT( '# GEOS-CHEM tracers [molec/cm2/s]'      )
 120  FORMAT( '# GEOS-CHEM tracers [molec/cm2]'        )
 130  FORMAT( '# GEOS-CHEM tracers [kg/s]'             )
 140  FORMAT( '# GEOS-CHEM tracers [kg]'               )
 150  FORMAT( '# GEOS-CHEM tracers [ug/m3]'            )
 160  FORMAT( '# SOA GPROD & APROD [kg/kg]'            )
 170  FORMAT( '# Quantities for Soil NOx restart file' )
 180  FORMAT( '# ND', i2.2, ' diagnostic quantities'   )
 190  FORMAT( '# TOMAS aerosol rate [cm-3 s-1]'      )   !(win, 7/9/09)
      
      END SUBROUTINE WRITE_SEPARATOR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diaginfo
!
! !DESCRIPTION: Subroutine INIT\_DIAGINFO initializes the CATEGORY, DESCRIPT, 
!  and OFFSET variables, which are used to define the "diaginfo.dat" file 
!  for GAMAP.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAGINFO
! 
! !REVISION HISTORY: 
!  17 Oct 1996 - R. Yantosca - Initial version
!  (1 ) Split this code off from INIT_GAMAP, for clarity.  Now declare biomass
!        burning emissions w/ offset of 45000.  Now declare time in the 
!        troposphere diagnostic with offset of 46000. (phs, bmy, 10/17/06)
!  (2 ) Now add IJ-GPROD & IJ-APROD w/ offset of SPACING*6, for the SOA
!        GPROD & APROD restart file. (tmf, havala, bmy, 2/6/07)
!  (3 ) Now declare H2-HD sources w/ offset of 48000. Now declare H2-HD
!        production/loss w/ offset of 47000. (phs, 9/18/07)
!  (4 ) Change diagnostic category for ND31 diagnostic from "PS-PTOP" 
!        to "PEDGE-$" (bmy, 11/16/07)
!  (5 ) Add categories CH4-LOSS, CH4-EMISS and WET-FRAC (kjw, 8/18/09)
!  (6 ) Add potential temperature category. (fp, 2/26/10)
!  21 May 2010 - C. Carouge  - Add diagnostic for mercury simulation
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!  21 Sep 2010 - R. Yantosca - Remove duplicate definitions of CV-FLX-$,
!                              TURBMC-$, EW-FLX-$, NS-FLX-$, UP-FLX-$
!  21 Oct 2010 - R. Yantosca - Bug fix: MC-FRC-$ should have an offset of 
!                              SPACING*3 since it has units of kg/s.
!  31 Oct 2012 - R. Yantosca - Save soil NOx restart quantities as SPACING*60
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA +
!                              semivolatile POA simulations (H. Pye)
!  03 Sep 2013 - R. Yantosca - Restored PS-PTOP to avoid issues when using
!                              diaginfo.dat, tracerinfo.dat w/ REGRIDH_RESTART
!  30 Jan 2014 - R. Yantosca - Add ND60 TOMAS diagnostic catetories
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER :: N

      !=================================================================
      ! INIT_DIAGINFO begins here!
      !=================================================================

      N           = 1
      CATEGORY(N) = 'IJ-AVG-$'
      DESCRIPT(N) = 'Tracer concentration'
      OFFSET(N)   = SPACING * 0

      N           = N + 1
      CATEGORY(N) = 'IJ-24H-$'
      DESCRIPT(N) = '24-hr avg tracer conc.'
      OFFSET(N)   = SPACING * 0
         
      N           = N + 1
      CATEGORY(N) = 'INST-MAP'
      DESCRIPT(N) = 'Instantaneous tracer'
      OFFSET(N)   = SPACING * 0

      N           = N + 1
      CATEGORY(N) = 'ANTHSRCE'
      DESCRIPT(N) = 'Anthropogenic emissions'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'BIOFSRCE'
      DESCRIPT(N) = 'Biofuel emissions'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NO-AC-$'
      DESCRIPT(N) = 'Aircraft NO'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NO-AN-$'
      DESCRIPT(N) = 'Anthropogenic NO'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NO-BIOB'
      DESCRIPT(N) = 'Biomass NO'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NO-BIOF'
      DESCRIPT(N) = 'Biofuel NO'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NO-LI-$'
      DESCRIPT(N) = 'Lightning NO'      
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NO-SOIL'
      DESCRIPT(N) = 'Soil NO'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NO-FERT'
      DESCRIPT(N) = 'Fertilizer NO'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NO-STRT'
      DESCRIPT(N) = 'Stratopsheric NO'
      OFFSET(N)   = SPACING * 1 

      N           = N + 1
      CATEGORY(N) = 'INST_COL'
      DESCRIPT(N) = 'Instantaneous columns'
      OFFSET(N)   = SPACING * 2

      N           = N + 1
      CATEGORY(N) = 'CV-FLX-$'
      DESCRIPT(N) = 'Convective mass flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'TURBMC-$'
      DESCRIPT(N) = 'PBL mixing mass flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'EW-FLX-$'
      DESCRIPT(N) = 'E/W transport flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'NS-FLX-$'
      DESCRIPT(N) = 'N/S transport flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'UP-FLX-$'
      DESCRIPT(N) = 'Up/down transport flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'STRT-FLX'
      DESCRIPT(N) = 'Flux from stratosphere'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'RN--SRCE'
      DESCRIPT(N) = 'Rn-Pb-Be source'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'RN-DECAY'
      DESCRIPT(N) = 'Rn-Pb-Be loss'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'WETDCV-$'
      DESCRIPT(N) = 'Conv wet scavenging'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'WETDLS-$'
      DESCRIPT(N) = 'Wet deposition'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'MC-FRC-$'
      DESCRIPT(N) = 'Moist conv fraction'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'DMS-BIOG'
      DESCRIPT(N) = 'Biogenic DMS'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'DUSTSRCE'
      DESCRIPT(N) = 'Dust emission'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NVOCSRCE'
      DESCRIPT(N) = 'NVOC emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SALTSRCE'
      DESCRIPT(N) = 'Seasalt emission'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-AC-$'
      DESCRIPT(N) = 'Aircraft SO2 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-AN-$'
      DESCRIPT(N) = 'Anthro SO2 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-BIOB'
      DESCRIPT(N) = 'Biomass SO2 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-BIOF'
      DESCRIPT(N) = 'Biofuel SO2 emissions'
      OFFSET(N)   = SPACING * 4 

      N           = N + 1
      CATEGORY(N) = 'SO2-EV-$'
      DESCRIPT(N) = 'Erup. Volcano SO2'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-NV-$'
      DESCRIPT(N) = 'Non-Erup. Volcano SO2'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-SHIP'
      DESCRIPT(N) = 'SO2 from ship exhaust'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO4-AN-$'
      DESCRIPT(N) = 'Anthro SO4 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO4-BIOF'
      DESCRIPT(N) = 'Biofuel SO4 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NH3-ANTH'
      DESCRIPT(N) = 'Anthro NH3 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NH3-NATU'
      DESCRIPT(N) = 'Natural NH3 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NH3-BIOB'
      DESCRIPT(N) = 'Biomass NH3 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NH3-BIOF'
      DESCRIPT(N) = 'Biofuel NH3 emissions'
      OFFSET(N)   = SPACING * 4

#if   defined( TOMAS )
      N           = N + 1               !(win, 7/9/09)
      CATEGORY(N) = 'NK-EMISS'
      DESCRIPT(N) = 'Size-res aerosol number emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1               !(win, 7/9/09)
      CATEGORY(N) = 'SF-EMISS'
      DESCRIPT(N) = 'Size-res sulfate mass emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1               !(win, 7/9/09)
      CATEGORY(N) = 'SS-EMISS'
      DESCRIPT(N) = 'Size-res sea-salt mass emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1               !(win, 7/9/09)
      CATEGORY(N) = 'ECIL-SRC'
      DESCRIPT(N) = 'Size-res H-phillic EC mass emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1               !(win, 7/9/09)
      CATEGORY(N) = 'ECOB-SRC'
      DESCRIPT(N) = 'Size-res H-phobic EC mass emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1               !(win, 7/9/09)
      CATEGORY(N) = 'OCIL-SRC'
      DESCRIPT(N) = 'Size-res H-phillic OC mass emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1               !(win, 7/9/09)
      CATEGORY(N) = 'OCOB-SRC'
      DESCRIPT(N) = 'Size-res H-phobic OC mass emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1               !(win, 7/9/09)
      CATEGORY(N) = 'DUST-SRC'
      DESCRIPT(N) = 'Size-res dust mass emissions'
      OFFSET(N)   = SPACING * 4

#endif

      N           = N + 1
      CATEGORY(N) = 'TROPO-AV'
      DESCRIPT(N) = 'Trop avg''d tracer'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'TCMASS-$'
      DESCRIPT(N) = 'Tracer mass (kg)'
      OFFSET(N)   = SPACING * 4

#if   defined( TOMAS ) 

      N           = N + 1
      CATEGORY(N) = 'TMS-COND'
      DESCRIPT(N) = 'Condensation rate'
      OFFSET(N)   = SPACING * 6

      N           = N + 1
      CATEGORY(N) = 'TMS-COAG'
      DESCRIPT(N) = 'Condensation rate'
      OFFSET(N)   = SPACING * 6

      N           = N + 1
      CATEGORY(N) = 'TMS-NUCL'
      DESCRIPT(N) = 'Condensation rate'
      OFFSET(N)   = SPACING * 6

      N           = N + 1
      CATEGORY(N) = 'TMS-AQOX'
      DESCRIPT(N) = 'Condensation rate'
      OFFSET(N)   = SPACING * 6

      N           = N + 1
      CATEGORY(N) = 'AERO-FIX'
      DESCRIPT(N) = 'Condensation rate'
      OFFSET(N)   = SPACING * 6

      N           = N + 1
      CATEGORY(N) = 'TMS-SOA'
      DESCRIPT(N) = 'Condensation rate'
      OFFSET(N)   = SPACING * 6

      N           = N + 1           !(win, 7/9/09)
      CATEGORY(N) = 'TOMAS-3D'
      DESCRIPT(N) = '3-D aerosol rate of a selected microphysics'
      OFFSET(N)   = SPACING * 7
#endif

      N           = N + 1
      CATEGORY(N) = 'PEDGE-$'
      DESCRIPT(N) = 'Pressure at level edges'
      OFFSET(N)   = SPACING * 10 

      N           = N + 1
      CATEGORY(N) = 'PS-PTOP'
      DESCRIPT(N) = 'Pressure level edges'
      OFFSET(N)   = SPACING * 10 

      N           = N + 1
      CATEGORY(N) = 'DAO-FLDS'
      DESCRIPT(N) = 'GMAO 2-D met fields'
      OFFSET(N)   = SPACING * 11

      N           = N + 1
      CATEGORY(N) = 'DAO-3D-$'
      DESCRIPT(N) = 'GMAO 3-D met fields'
      OFFSET(N)   = SPACING * 12

      N           = N + 1
      CATEGORY(N) = 'JV-MAP-$'
      DESCRIPT(N) = 'Photolysis rates'
      OFFSET(N)   = SPACING * 13

      N           = N + 1
      CATEGORY(N) = 'OD-MAP-$'
      DESCRIPT(N) = 'Optical Depths'
      OFFSET(N)   = SPACING * 14

      N           = N + 1
      CATEGORY(N) = 'LANDMAP'
      DESCRIPT(N) = 'Land type map'
      OFFSET(N)   = SPACING * 15

      N           = N + 1
      CATEGORY(N) = 'CHEM-L=$'
      DESCRIPT(N) = 'Chemical Prod/Loss'
      OFFSET(N)   = SPACING * 16

      N           = N + 1
      CATEGORY(N) = 'PORL-L=$'
      DESCRIPT(N) = 'ND65 P/L family diagnostics'
      OFFSET(N)   = SPACING * 17

      N           = N + 1
      CATEGORY(N) = 'PL-SUL=$'
      DESCRIPT(N) = 'P/L of sulfur species'
      OFFSET(N)   = SPACING * 18

      N           = N + 1
      CATEGORY(N) = 'TIME-SER'
      DESCRIPT(N) = 'Timeseries quantities'
      OFFSET(N)   = SPACING * 19

      N           = N + 1
      CATEGORY(N) = 'CO--SRCE'
      DESCRIPT(N) = 'CO Source diagnostic'
      OFFSET(N)   = SPACING * 20

      N           = N + 1
      CATEGORY(N) = 'BIOGSRCE'
      DESCRIPT(N) = 'Biogenic emissions'
      OFFSET(N)   = SPACING * 21

      N           = N + 1
      CATEGORY(N) = 'ACETSRCE'
      DESCRIPT(N) = 'Acetone emissions'
      OFFSET(N)   = SPACING * 22

      N           = N + 1
      CATEGORY(N) = 'EMDIS-BL'
      DESCRIPT(N) = 'Emissions in PBL'
      OFFSET(N)   = SPACING * 23

      N           = N + 1
      CATEGORY(N) = 'BXHGHT-$'
      DESCRIPT(N) = 'Boxheight, airmass, etc'
      OFFSET(N)   = SPACING * 24

      N           = N + 1
      CATEGORY(N) = 'DXYP'
      DESCRIPT(N) = 'Surface area'
      OFFSET(N)   = SPACING * 25

      N           = N + 1
      CATEGORY(N) = 'TR-PAUSE'
      DESCRIPT(N) = 'Annual mean tropopause'
      OFFSET(N)   = SPACING * 26

      N           = N + 1
      CATEGORY(N) = 'PBLDEPTH'
      DESCRIPT(N) = 'Afternoon PBL height'
      OFFSET(N)   = SPACING * 27

      N           = N + 1
      CATEGORY(N) = 'WD-FRC-$'
      DESCRIPT(N) = 'Wet dep fraction'
      OFFSET(N)   = SPACING * 28

      N           = N + 1
      CATEGORY(N) = 'WD-LSR-$'
      DESCRIPT(N) = 'Large-scale rainout'
      OFFSET(N)   = SPACING * 29

      N           = N + 1
      CATEGORY(N) = 'WD-CVR-$'
      DESCRIPT(N) = 'Convective rainout'
      OFFSET(N)   = SPACING * 29

      N           = N + 1
      CATEGORY(N) = 'WD-LSW-$'
      DESCRIPT(N) = 'Large-scale washout'
      OFFSET(N)   = SPACING * 29

      N           = N + 1
      CATEGORY(N) = 'WD-CVW-$'
      DESCRIPT(N) = 'Convective washout'
      OFFSET(N)   = SPACING * 29

      N           = N + 1
      CATEGORY(N) = 'COBUDGET'
      DESCRIPT(N) = 'bnd CO-OH budget'
      OFFSET(N)   = SPACING * 31

      N           = N + 1
      CATEGORY(N) = 'CH4-LOSS'
      DESCRIPT(N) = 'CH4 Loss by OH'
      OFFSET(N)   = SPACING * 32

      N           = N + 1
      CATEGORY(N) = 'BC-ANTH'
      DESCRIPT(N) = 'Anthro BC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'BC-BIOB'
      DESCRIPT(N) = 'Biomass BC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'BC-BIOF'
      DESCRIPT(N) = 'Biofuel BC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-ANTH'
      DESCRIPT(N) = 'Anthro OC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-BIOB'
      DESCRIPT(N) = 'Biomass OC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-BIOF'
      DESCRIPT(N) = 'Biofuel OC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-BIOG'
      DESCRIPT(N) = 'Biogenic OC emission'
      OFFSET(N)   = SPACING * 33

      ! SOAupdate: update for new mtp (hotp 5/25/10)
      N           = N + 1
      CATEGORY(N) = 'OC-MTPA'
      DESCRIPT(N) = 'Biogenic MTPA emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-LIMO'
      DESCRIPT(N) = 'Biogenic ALPH emission'
      OFFSET(N)   = SPACING * 33

      ! SOAupdate: update for new mtp (hotp 5/25/10)
      N           = N + 1
      CATEGORY(N) = 'OC-MTPO'
      DESCRIPT(N) = 'Biogenic MTPO emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-SESQ'
      DESCRIPT(N) = 'Biogenic SESQ emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'PL-BC=$'
      DESCRIPT(N) = 'H-philic from H-phobic BC'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'PL-OC=$'
      ! SOAupdate (hotp 5/25/09)
      ! updated name for SOA (hotp 7/31/08)
      !DESCRIPT(N) = 'H-philic from H-phobic OC'
      DESCRIPT(N) = 'Philic OC prod/SOA prod'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'SALT-SR$'
      DESCRIPT(N) = 'Sea salt emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'BLKC-SR$'
      DESCRIPT(N) = 'Black carbon emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'ORGC-SR$'
      DESCRIPT(N) = 'Organic Carbon emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'SOAGM=$'
      DESCRIPT(N) = 'Dicarbonyl SOA production'
      OFFSET(N)   = SPACING * 33

!- eds 9/1/10 ----------------------------------------------
      N           = N + 1
      CATEGORY(N) = 'HG-SRCE'
      DESCRIPT(N) = 'Hg emissions'
      OFFSET(N)   = SPACING * 34

      N           = N + 1
      CATEGORY(N) = 'HG0-ANTH'
      DESCRIPT(N) = 'Hg(0) Anthro Emissions'
      OFFSET(N)   = SPACING * 61

      N           = N + 1
      CATEGORY(N) = 'HG0-AQUA'
      DESCRIPT(N) = 'Hg(0) Ocean Mass'
      OFFSET(N)   = SPACING * 62

      N           = N + 1
      CATEGORY(N) = 'HG0-RECY'
      DESCRIPT(N) = 'Hg(0) Land Prompt Recycling'
      OFFSET(N)   = SPACING * 63

      N           = N + 1
      CATEGORY(N) = 'HGNET-OC'
      DESCRIPT(N) = 'Hg(0) Net Ocean Emissions'
      OFFSET(N)   = SPACING * 64

      N           = N + 1
      CATEGORY(N) = 'HG0-GEOG'
      DESCRIPT(N) = 'Hg(0) Geogenic Emissions'
      OFFSET(N)   = SPACING * 65

      N           = N + 1
      CATEGORY(N) = 'HG2-ANTH'
      DESCRIPT(N) = 'Hg(II) Anthro Emissions'
      OFFSET(N)   = SPACING * 66

      N           = N + 1
      CATEGORY(N) = 'HG2-AQUA'
      DESCRIPT(N) = 'Hg(II) Ocean Mass'
      OFFSET(N)   = SPACING * 67

      N           = N + 1
      CATEGORY(N) = 'HG2-SINK'
      DESCRIPT(N) = 'Hg(II) Sinking'
      OFFSET(N)   = SPACING * 68

      N           = N + 1
      CATEGORY(N) = 'HGP-ANTH'
      DESCRIPT(N) = 'Hg(P) Anthro Emissions'
      OFFSET(N)   = SPACING * 69

      N           = N + 1
      CATEGORY(N) = 'HGT-AQUA'
      DESCRIPT(N) = 'Hg Total Ocean Mass'
      OFFSET(N)   = SPACING * 70

      N           = N + 1
      CATEGORY(N) = 'HGP-AQUA'
      DESCRIPT(N) = 'Hg(P) Ocean Mass'
      OFFSET(N)   = SPACING * 71

      N           = N + 1
      CATEGORY(N) = 'HGP-CONV'
      DESCRIPT(N) = 'Hg(II) Aqeous Conversion to Particle'
      OFFSET(N)   = SPACING * 72

      N           = N + 1
      CATEGORY(N) = 'HG0-BURN'
      DESCRIPT(N) = 'Hg(0) Biomass Burning Emissions'
      OFFSET(N)   = SPACING * 73

      N           = N + 1
      CATEGORY(N) = 'HG0-VEGT'
      DESCRIPT(N) = 'Hg(0) Veg Transp. Emissions'
      OFFSET(N)   = SPACING * 74

      N           = N + 1
      CATEGORY(N) = 'HG0-SOIL'
      DESCRIPT(N) = 'Hg(0) Soil Emissions'
      OFFSET(N)   = SPACING * 75

      N           = N + 1
      CATEGORY(N) = 'HG0-FXUP'
      DESCRIPT(N) = 'Hg(0) Gross Ocean Up Flux'
      OFFSET(N)   = SPACING * 76

      N           = N + 1
      CATEGORY(N) = 'HG0-FXDN'
      DESCRIPT(N) = 'Hg(0) Gross Ocean Down Flux'
      OFFSET(N)   = SPACING * 77

      N           = N + 1
      CATEGORY(N) = 'HG0-SNOW'
      DESCRIPT(N) = 'Hg(0) Snow Emissions'
      OFFSET(N)   = SPACING * 78

      N           = N + 1
      CATEGORY(N) = 'HG-NETOX'
      DESCRIPT(N) = 'Hg(0) to Hg(II) Net Oxidation'
      OFFSET(N)   = SPACING * 79

      N           = N + 1
      CATEGORY(N) = 'HG2-OXOH'
      DESCRIPT(N) = 'Hg(0) Oxidation to Hg(II) by OH'
      OFFSET(N)   = SPACING * 80

      N           = N + 1
      CATEGORY(N) = 'HG2-OXO3'
      DESCRIPT(N) = 'Hg(0) Oxidation to Hg(II) by O3'
      OFFSET(N)   = SPACING * 81

      N           = N + 1
      CATEGORY(N) = 'HG2-SALT'
      DESCRIPT(N) = 'Hg(II) Loss by Sea Salt'
      OFFSET(N)   = SPACING * 82

      N           = N + 1
      CATEGORY(N) = 'HG2-SSRX'
      DESCRIPT(N) = 'Hg(0) Sea Salt Reaction Rate'
      OFFSET(N)   = SPACING * 83

      N           = N + 1
      CATEGORY(N) = 'HG2-OXBR'
      DESCRIPT(N) = 'Hg(0) Oxidation to Hg(II) by Br'
      OFFSET(N)   = SPACING * 84

!---------------------------------------------------------

      N           = N + 1
      CATEGORY(N) = 'PL-HG2-$'
      DESCRIPT(N) = 'Prod / loss of Hg2'
      OFFSET(N)   = SPACING * 35

      N           = N + 1
      CATEGORY(N) = 'DRYD-FLX'
      DESCRIPT(N) = 'Drydep fluxes'
      OFFSET(N)   = SPACING * 36

      N           = N + 1
      CATEGORY(N) = 'DRYD-VEL'
      DESCRIPT(N) = 'Drydep velocities'
      OFFSET(N)   = SPACING * 37

      N           = N + 1
      CATEGORY(N) = 'HCN-PL-$'
      DESCRIPT(N) = 'HCN & CH3CN sinks'
      OFFSET(N)   = SPACING * 38

      N           = N + 1
      CATEGORY(N) = 'HCN-SRCE'
      DESCRIPT(N) = 'HCN & CH3CN sources'
      OFFSET(N)   = SPACING * 39

      N           = N + 1
      CATEGORY(N) = 'CO2-SRCE'
      DESCRIPT(N) = 'CO2 fluxes'
      OFFSET(N)   = SPACING * 40

      N           = N + 1
      CATEGORY(N) = 'OCEAN-HG'
      DESCRIPT(N) = 'Oceanic Hg emissions'
      OFFSET(N)   = SPACING * 41

      N           = N + 1
      CATEGORY(N) = 'LFLASH-$'
      DESCRIPT(N) = 'Lightning flash rates'
      OFFSET(N)   = SPACING * 42

      N           = N + 1
      CATEGORY(N) = 'IJ-SOA-$'
      DESCRIPT(N) = 'SOA concentrations'
      OFFSET(N)   = SPACING * 43

      N           = N + 1
      CATEGORY(N) = 'CH3ISRCE'
      DESCRIPT(N) = 'CH3I emissions'
      OFFSET(N)   = SPACING * 44

      N           = N + 1
      CATEGORY(N) = 'BIOBSRCE'
      DESCRIPT(N) = 'Biomass emissions'
      OFFSET(N)   = SPACING * 45

      N           = N + 1
      CATEGORY(N) = 'TIME-TPS'
      DESCRIPT(N) = 'Fraction of time in troposphere'
      OFFSET(N)   = SPACING * 46

      N           = N + 1
      CATEGORY(N) = 'PL-H2HD-'
      DESCRIPT(N) = 'Prod / loss of H2-HD'
      OFFSET(N)   = SPACING * 47

      N           = N + 1
      CATEGORY(N) = 'H2HD-SRC'
      DESCRIPT(N) = 'H2 HD emissions'
      OFFSET(N)   = SPACING * 48

      ! New ND52 diagnostic for gamma(HO2). (jaegle, 2/29/09)
      N           = N + 1
      CATEGORY(N) = 'GAMMA'
      DESCRIPT(N) = 'gamma HO2'
      OFFSET(N)   = SPACING * 49

      ! For mercury simulation only so we can use same spacing. (ccc, 5/21/10)
      N           = N + 1
      CATEGORY(N) = 'SNOW-HG'
      DESCRIPT(N) = 'Hg mass in snow and ice'
      OFFSET(N)   = SPACING * 49

!(FP 06/19/2009)
      N           = N + 1
      CATEGORY(N) = 'THETA-$'
      DESCRIPT(N) = 'Potential temperature'
      OFFSET(N)   = SPACING * 57

      N           = N + 1
      CATEGORY(N) = 'CH4-EMIS'
      DESCRIPT(N) = 'CH4 Emissions'
      OFFSET(N)   = SPACING * 50

      N           = N + 1
      CATEGORY(N) = 'WET-FRAC'
      DESCRIPT(N) = 'Wetland Fraction'
      OFFSET(N)   = SPACING * 51

      N           = N + 1
      CATEGORY(N) = 'SHIP-$$$'
      DESCRIPT(N) = 'Ship diagnostics'
      OFFSET(N)   = SPACING * 59

      N           = N + 1
      CATEGORY(N) = 'RST-SOIL'
      DESCRIPT(N) = 'Soil NOx restart file'
      OFFSET(N)   = SPACING * 60

      ! Added for POPs (clf, 3/2/2011)
      N           = N + 1
      CATEGORY(N) = 'PG-SRCE'
      DESCRIPT(N) = 'POPs Emissions'
      OFFSET(N)   = SPACING * 52

      N           = N + 1
      CATEGORY(N) = 'PG-PP-$'
      DESCRIPT(N) = 'POP partitioning/oxidation'
      OFFSET(N)   = SPACING * 53

! NEED TO DO THAT AT ONE POINT      
!      ! New for CSPEC
!      N           = N + 1
!      CATEGORY(N) = 'IJ-CHK-$'
!      DESCRIPT(N) = 'species concentrations after chemistry'
!      OFFSET(N)   = SPACING * 50
!      
      ! Number of categories
      NCATS = N
      
      END SUBROUTINE INIT_DIAGINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tracerinfo
!
! !DESCRIPTION: Subroutine INIT\_TRACERINFO initializes the NAME, FNAME, 
!  MWT, MOLC, INDEX, MOLC, UNIT arrays which are used to define the 
!  "tracerinfo.dat" file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TRACERINFO( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIAG03_MOD,         ONLY : ND03, PD03, PD03_PL
      USE DIAG04_MOD,         ONLY : ND04
      USE DIAG41_MOD,         ONLY : ND41
      USE DIAG42_MOD,         ONLY : ND42
      USE DIAG48_MOD,         ONLY : DO_SAVE_DIAG48
      USE DIAG49_MOD,         ONLY : DO_SAVE_DIAG49
      USE DIAG50_MOD,         ONLY : DO_SAVE_DIAG50
      USE DIAG51_MOD,         ONLY : DO_SAVE_DIAG51
      USE DIAG51b_MOD,        ONLY : DO_SAVE_DIAG51b
      USE DIAG53_MOD,         ONLY : ND53, PD53
      USE DIAG56_MOD,         ONLY : ND56
      USE DIAG63_MOD,         ONLY:  DO_SAVE_DIAG63
      USE DIAG_PL_MOD,        ONLY : DO_SAVE_PL,  GET_NFAM
      USE DIAG_PL_MOD,        ONLY : GET_FAM_MWT, GET_FAM_NAME
      USE DRYDEP_MOD,         ONLY : DEPNAME,     NUMDEP,    NTRAIND
      USE GET_POPSINFO_MOD,   ONLY : GET_POP_XMW
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TRACERID_MOD
      USE WETSCAV_MOD,        ONLY : GET_WETDEP_IDWETD
      USE WETSCAV_MOD,        ONLY : GET_WETDEP_NSOL
#if   defined( TOMAS )
      USE TOMAS_MOD,          ONLY : IBINS, ICOMP,   IDIAG  !(win, 7/14/09)
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  17 Oct 1996 - R. Yantosca & P. Le Sager - Initial version
!  (1 ) Split this code off from INIT_GAMAP, for clarity.  Also now declare
!        biomass burning emissions w/ offset of 45000.  Bug fix: write out
!        26 tracers for ND48, ND49, ND50, ND51 timeseries.  Also define
!        ND54 diagnostic with offset of 46000. (bmy, 10/17/06)
!  (2 ) Modifications for H2/HD in ND10, ND44 diagnostics (phs, 9/18/07)
!  (3 ) Now write out PBLDEPTH diagnostic information to "tracerinfo.dat" if 
!        any of ND41, ND48, ND49, ND50, ND51 are turned on.  Also set the
!        unit to "kg/s" for the Rn-Pb-Be ND44 drydep diag. (cdh, bmy, 2/22/08)
!  (4 ) Added C2H4 in ND46 (ccc, 2/2/09)
!  (5 ) Add EFLUX to ND67 (lin, ccc, 5/29/08)
!  (6 ) Bug fix in ND28: ALD2 should have 2 carbons, not 3.  Also bug fix
!        in ND66 to print out the name of ZMMU correctly. (dbm, bmy, 10/9/09)
!  (7 ) Previous bug fix was erroneous; now corrected (dkh, bmy, 11/19/09)
!  (8 ) Include second satellite overpass diagnostic.  Adjust AOD name to 550 
!        nm from 400 nm.  Add additional dust AOD bins (amv, bmy, 12/18/09)
!  (9 ) Add iniformation for ND61 (win, 7/9/09)
!  (10) Manually add info for ND44 when TOMAS aerosols dry deposition because
!        only numbers are real drydep species while the mass species just 
!        tag along but we need them to show up in diagnostic too.  Reference  
!        to TOMAS_MOD and use IDTNK1 from TRACERID_MOD (win, 7/14/09)
!  20 Jul 2010 - C. Carouge  - Modifications to ND03 for mercury.
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  02 Sep 2010 - R. Yantosca - In ND28: Omit SOA tracers if LSOA = .FALSE.
!  09-Dec-2010 - H. Amos     - Added RGM and PBM  tracers for the
!                              mercury simulation
!  09-Dec-2010 - H. Amos     - fix spacing and #s for PL-HG2-$ diagnostics
!  12 Nov 2010 - R. Yantosca - Need to save out surface pressure line to 
!                              tracerinfo.dat for the timeseries diagnostics
!  22 Mar 2011 - C. Friedman - Added POPs
!  24 Jan 2012 - M. Payer    - Change scale factors for Rn-Pb-Be simulation so
!                              units are in mBq/SCM
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  19 Mar 2012 - M. Payer    - Remove ACETdl and ACETgr from ND11 diagnostic.
!                              Acetone from dry leaf matter and grasses is now
!                              included in the direct emissions (ACETbg).
!                              (E. Fischer)
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox with NO and O3 as part of
!                              removal of NOx-Ox partitioning. Removed code for
!                              storing pure O3 as N_TRACERS+1.
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA +
!                              semivolatile POA simulations (H. Pye)
!  04 Sep 2013 - R. Yantosca - Fix ND44 tracer listing
!  04 Sep 2013 - R. Yantosca - Now also reference everything in tracerid_mod.F
!  04 Sep 2013 - R. Yantosca - Bug fix: # of ND46 tracers should be 20
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  31 Jan 2014 - R. Yantosca - Now save out tracerinfo.dat properly for TOMAS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N, NN, NYMDb, NHMSb, T
      LOGICAL :: DO_TIMESERIES

      REAL*8,     SAVE      :: POP_XMW
      REAL*8                :: DUM

      !=================================================================
      ! INIT_TRACERINFO begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      RC = GIGC_SUCCESS

      DUM = 1.0

      ! Set a flag if any timeseries diagnostics are turned on
      DO_TIMESERIES = ( DO_SAVE_DIAG48 .or. DO_SAVE_DIAG49 .or.
     &                  DO_SAVE_DIAG50 .or. DO_SAVE_DIAG51 .or.
     &                  DO_SAVE_DIAG51b )

      !----------------------------------
      ! General tracer information
      !----------------------------------
      DO T = 1, Input_Opt%N_TRACERS

         ! Store quantities for each tracer
         NAME (T,45) = Input_Opt%TRACER_NAME(T)
         FNAME(T,45) = TRIM( NAME(T,45) ) // ' tracer'
         MOLC (T,45) = INT( Input_Opt%TRACER_COEFF(T,1) )
         SCALE(T,45) = 1.0e+9
         INDEX(T,45) = T   !changed from N (phs, 3/19/03)
         MWT  (T,45) = Input_Opt%TRACER_MW_KG(T) 

         ! Special handling for hydrocarbons
         IF ( MOLC(T,45) > 1 ) THEN
            UNIT(T,45) = 'ppbC'
         ELSE
            UNIT(T,45) = 'ppbv'
         ENDIF

         ! Special handling for Rn-Pb-Be simulation (bmy, 5/11/05)
         IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN
            SELECT CASE( T ) 
               CASE( 1 )
                  UNIT (T,45) = 'mBq/SCM'
                  SCALE(T,45) = 5.6397E+22
               CASE( 2 ) 
                  UNIT (T,45) = 'mBq/SCM'
                  SCALE(T,45) = 2.6141E+19
               CASE( 3 )
                  UNIT (T,45) = 'mBq/SCM'
                  SCALE(T,45) = 4.0513E+21
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDIF

         ! For mercury, print as pptv (bmy, 1/24/06)
         IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
            UNIT(T,45)  = 'pptv'
            SCALE(T,45) = 1.0e+12
         ENDIF
      ENDDO

      ! Number of ND45 tracers 
      NTRAC(45) = Input_Opt%N_TRACERS

      !-------------------------------------
      ! Hg source, production & loss (ND03)
      !
      ! Updated for tagged Hg simulation
      ! (cdh, bmy, 1/9/06)
      ! Modified to avoid over-sized arrays.
      ! (ccc, 7/20/10)
      !-------------------------------------

      !=====================================
      ! Standard Hg simulation
      !=====================================
      IF ( (ND03 > 0) .AND. .NOT. (Input_Opt%LSPLIT) ) THEN !eds 9/7/10
         
         ! Number of tracers:
         ! PD03     = max # of tracers for HG-SRCE
         ! PD03_PL  = max # of tracers for PL-HG2-$

         ! Loop over tracers for HG-SRCE
!         DO T = 1, NTRAC(03)
         DO T = 1, PD03

            ! Define quantities
            UNIT (T,03) = 'kg'
            MOLC (T,03) = 1
            MWT  (T,03) = 201e-3
            SCALE(T,03) = 1e0

            ! Get name, long-name, index, and new units
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,03) = 'Hg0_an'
                  FNAME(T,03) = 'Anthro elemental Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 2  )
                  NAME (T,03) = 'Hg0_aq'
                  FNAME(T,03) = 'Ocean mass of elemental Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 3  )
                  NAME (T,03) = 'Hg0_oc'
                  FNAME(T,03) = 'Ocean-emitted elemental Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 4  )
                  NAME (T,03) = 'Hg0_ln'
                  FNAME(T,03) = 'Land re-emitted elemental Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 5  )
                  NAME (T,03) = 'Hg0_na'
                  FNAME(T,03) = 'Natural land source'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 6  )
                  NAME (T,03) = 'Hg2_an'
                  FNAME(T,03) = 'Anthro divalent Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 7  )
                  NAME (T,03) = 'Hg2_aq'
                  FNAME(T,03) = 'Ocean mass of divalent Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 8  )
                  NAME (T,03) = 'Hg2_sk'
                  FNAME(T,03) = 'Mass of Hg2 sunk in ocean'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 9  )
                  NAME (T,03) = 'HgP_an'
                  FNAME(T,03) = 'Anthro particulate Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 10 )
                  NAME (T,03) = 'Hgaq_tot'
                  FNAME(T,03) = 'Total oceanic mercury'
                  UNIT (T,03) = 'kg'
                  INDEX(T,03) = T + ( SPACING * 34 )
                  MWT  (T,03) = 0e0
               CASE( 11 )
                  NAME (T,03) = 'HgP_aq'
                  FNAME(T,03) = 'Ocean mass of particulate Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 12 )
                  NAME (T,03) = 'JorgC'
                  FNAME(T,03) = 'Mass of orgC sunk in ocean'
                  UNIT (T,03) = 'kg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 13 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg biomass burn emis'
                  UNIT (T,03) = 'kg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 14 )
                  NAME (T,03) = 'Hg_vg'
                  FNAME(T,03) = 'Hg vegetation emissions'
                  UNIT (T,03) = 'kg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 15 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg soil emissions'
                  UNIT (T,03) = 'kg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 16 )
                  NAME (T,03) = 'Hg_up'
                  FNAME(T,03) = 'Hg ocean up flux'
                  UNIT (T,03) = 'kg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 17 )
                  NAME (T,03) = 'Hg_down'
                  FNAME(T,03) = 'Hg ocean downflux'
                  UNIT (T,03) = 'kg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 18 )
                  NAME (T,03) = 'Hg0_snow'
                  FNAME(T,03) = 'Snow emission of Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO
         
         ! Loop over tracers for PL-HG2-$
         DO N = 1, PD03_PL
            T = N + PD03

            ! Define quantities
            UNIT (T,03) = 'kg'
            MOLC (T,03) = 1
            MWT  (T,03) = 201e-3
            SCALE(T,03) = 1e0


            ! Get name, long-name, index, and new units
            SELECT CASE( N )
               CASE( 1 )
                  NAME (T,03) = 'Hg2_Hg0'
                  FNAME(T,03) = 'Prod of Hg2 from Hg0'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 )
               CASE( 2 )
                  NAME (T,03) = 'Hg2_OH'
                  FNAME(T,03) = 'Prod of Hg2 from OH'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 ) 
               CASE( 3 )
                  NAME (T,03) = 'Hg2_O3'
                  FNAME(T,03) = 'Prod of Hg2 from O3'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 )
               CASE( 4 )
                  NAME (T,03) = 'Hg2_SS'
                  FNAME(T,03) = 'Loss of Hg2 from sea salt'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 )
               CASE( 5 )
                  NAME (T,03) = 'Hg2_SSR'
                  FNAME(T,03) = 'Loss rate Hg2 from sea salt'
                  UNIT (T,03) = '/s'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 )
               CASE( 6 )
                  NAME (T,03) = 'Hg2_Br'
                  FNAME(T,03) = 'Prod of Hg2 from Br'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 )
               CASE( 7 )
                  NAME (T,03) = 'Br'
                  FNAME(T,03) = 'Br concentration'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 )
               CASE( 8 )
                  NAME (T,03) = 'BrO'
                  FNAME(T,03) = 'BrO concentration'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 )
               CASE( 9 )
                  NAME (T,03) = 'PBM'
                  FNAME(T,03) = 'Particle bound mercury'
                  UNIT (T,03) = 'pptv'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 )      
               CASE( 10 )
                  NAME (T,03) = 'RGM'
                  FNAME(T,03) = 'Reactive gaseous mercury'
                  UNIT (T,03) = 'pptv'
                  INDEX(T,03) = ( N ) + ( SPACING * 35 ) 
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO

         ! Loop over tracers for OCEAN-HG
         DO N = 1, Input_Opt%N_TRACERS
         
            T = N + PD03 + PD03_PL

            ! Define quantities
            UNIT (T,03) = 'kg'
            MOLC (T,03) = 1
            MWT  (T,03) = 201e-3
            SCALE(T,03) = 1e0


            NAME (T,03) = Input_Opt%TRACER_NAME(N) 

!--- Keep HgP for tracer 3 in the ocean now on. (ccc, 7/20/10)
!                  ! Tracer 3 should be "HgC" instead of "HgP"
!                  IF ( TRIM( NAME(T,03) ) == 'HgP' ) THEN
!                     NAME(T,03) = 'HgC'
!                  ENDIF

            FNAME(T,03) = 'Oceanic ' // TRIM( NAME(T,03) )
            INDEX(T,03) = ( N ) + ( SPACING * 41 )
         ENDDO

         ! Loop over tracers for SNOW-HG
         DO N = 1, 1
         
            T = N + PD03 + PD03_PL + Input_Opt%N_TRACERS

            ! Define quantities
            UNIT (T,03) = 'kg'
            MOLC (T,03) = 1
            MWT  (T,03) = 201e-3
            SCALE(T,03) = 1e0


            NAME (T,03) = 'Snow_Hg'
            FNAME(T,03) = 'Tot. Hg in snowpack'
            INDEX(T,03) = ( N ) + ( SPACING * 49 )
         ENDDO
         
         ! Total max number of tracers for DIAG 03.
         NTRAC(03) = T
         
      ENDIF

      !=====================================
      ! Tagged Hg simulation
      !=====================================
      IF ( (ND03 > 0) .AND. (Input_Opt%LSPLIT) ) THEN !eds 9/7/10

         ! Define quantities
         UNIT (T,03) = 'kg'
         MOLC (T,03) = 1
         MWT  (T,03) = 201e-3
         SCALE(T,03) = 1e0

       	 !-------------! 
         ! HG0-ANTH    !
       	 !-------------! 

         T = 1
         SELECT CASE( T )
           CASE( 1 )
               NAME (T,03) = 'Hg0'
               FNAME(T,03) = 'Hg0'
               INDEX(T,03) = T + ( SPACING * 61 )  !eds 9/4/10
         END SELECT  


         ! Number of tracers
!         NTRAC(03) = N_Hg_CATS
         NTRAC(03) =  480

       	 !-------------! 
         ! HG0-AQUA    !
       	 !-------------! 

         ! Loop over tracers
!         DO N = 1, NTRAC(03)
         DO N = 1, N_Hg_CATS

         T = N + 1

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 2 )
                  NAME (T,03) = 'Hg0'
                  FNAME(T,03) = 'Hg0'
                  INDEX(T,03) = N + ( SPACING * 62 )  !eds 9/3/10 
              CASE( 3 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
       	       	  INDEX(T,03) =	N + ( SPACING *	62 )  
              CASE( 4 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 62 ) 
              CASE( 5 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 62 )  
              CASE( 6 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 7 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 8 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 9 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 10 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 11 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 12 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 13 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 14 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 15 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 16 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 17 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 18 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 19 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 20 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 21 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 22 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 23 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 24 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 25 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 26 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 27 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 28 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 29 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE( 30 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 62 )
              CASE DEFAULT
            END SELECT
         ENDDO

         !-------------!
         ! HG0-RECY    !
         !-------------!

         ! Loop over tracers
!         DO N = 1, NTRAC(03)
         DO N = 1, N_Hg_CATS

         T = N + 1 + N_Hg_CATS

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 31 )
                  NAME (T,03) = 'Hg0'
                  FNAME(T,03) = 'Hg0'
                  INDEX(T,03) = N + ( SPACING * 63 ) 
              CASE( 32 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 33 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 63 ) 
              CASE( 34 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 63 )  
              CASE( 35 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 36 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 37 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 38 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 39 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 40 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 41 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 42 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 43 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 44 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 45 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 46 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 47 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 48 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 49 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 50 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 51 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 52 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 53 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 54 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 55 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 56 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 57 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 58 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE( 59 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 63 )
              CASE DEFAULT
            END SELECT
         ENDDO

         !-------------!
         ! HGNET-OC    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 1 + (2d0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 60 )
                  NAME (T,03) = 'Hg0'
                  FNAME(T,03) = 'Hg0'
                  INDEX(T,03) = N + ( SPACING * 64 ) 
              CASE( 61 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 62 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 64 ) 
              CASE( 63 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 64 )  
              CASE( 64 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 65 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 66 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 67 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 68 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 69 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 70 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 71 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 72 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 73 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 74 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 75 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 76 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 77 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 78 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 79 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 80 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 81 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 82 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 83 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 84 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 85 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 86 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 87 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE( 88 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 64 )
              CASE DEFAULT
            END SELECT
         ENDDO

        

         !-------------!
         ! HG0-GEOG    !
         !-------------!

         T = 89
         SELECT CASE( T )
           CASE( 89 )
               NAME (T,03) = 'Hg0'
               FNAME(T,03) = 'Hg0'
               INDEX(T,03) = 1 + ( SPACING * 65 )          
         END SELECT


         !-------------!
         ! HG2-ANTH    !
         !-------------!

         T = 90
         N = 1
         SELECT CASE( T )
           CASE( 90 )
               NAME (T,03) = 'Hg2'
               FNAME(T,03) = 'Hg2'
               INDEX(T,03) = N + ( SPACING * 66 )    
         END SELECT


         !-------------!
         ! HG2-AQUA    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 3 + (3D0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 91 )
                  NAME (T,03) = 'Hg2'
                  FNAME(T,03) = 'Hg2'
                  INDEX(T,03) = N + ( SPACING * 67 ) 
              CASE( 92 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 93 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 67 ) 
              CASE( 94 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 67 )  
              CASE( 95 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 96 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 97 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 98 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 99 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 100 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 101 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 102 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 103 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 104 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 105 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 106 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 107 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 108 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 109 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 110 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 111 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 112 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 113 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 114 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 115 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 116 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 117 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 118 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE( 119 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 67 )
              CASE DEFAULT
            END SELECT
         ENDDO


         !-------------!
         ! HG2-SINK    !
         !-------------!

         T = 120
         N = 1
         SELECT CASE( T )
           CASE( 120 )
               NAME (T,03) = 'Hg2'
               FNAME(T,03) = 'Hg2'
               INDEX(T,03) = N + ( SPACING * 68 )          
         END SELECT


         !-------------!
         ! HGP-ANTH    !
         !-------------!

         T = 121
         N = 1
         SELECT CASE( T )
           CASE( 121 )
               NAME (T,03) = 'HgP'
               FNAME(T,03) = 'HgP'
               INDEX(T,03) = N + ( SPACING * 69 )          
         END SELECT


         !-------------!
         ! HGT-AQUA    !
         !-------------!

         T = 122
         N = 1
         SELECT CASE( T )
           CASE( 122 )
               NAME (T,03) = 'HgTOT'
               FNAME(T,03) = 'Hg Total'
               INDEX(T,03) = N + ( SPACING * 70 )          
         END SELECT



         !-------------!
         ! HGP-AQUA    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 6 + (4D0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 123 )
                  NAME (T,03) = 'HgP'
                  FNAME(T,03) = 'HgP'
                  INDEX(T,03) = N + ( SPACING * 71 ) 
              CASE( 124 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 125 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 71 ) 
              CASE( 126 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 71 )  
              CASE( 127 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 128 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 129 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 130 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 131 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 132 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 133 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 134 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 135 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 136 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 137 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 138 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 139 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 140 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 141 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 142 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 143 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 144 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 145 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 146 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 147 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 148 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 149 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 150 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE( 151 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 71 )
              CASE DEFAULT
            END SELECT
         ENDDO


         !-------------!
         ! HGP-CONV    !
         !-------------!

         T = 152
         N = 1
         SELECT CASE( T )
           CASE( 152 )
               NAME (T,03) = 'HgP'
               FNAME(T,03) = 'HgP'
               INDEX(T,03) = N + ( SPACING * 72 )          
         END SELECT


         !-------------!
         ! HG0-BURN    !
         !-------------!

         T = 153
         N = 1
         SELECT CASE( T )
           CASE( 153 )
               NAME (T,03) = 'Hg0'
               FNAME(T,03) = 'Hg0'
               INDEX(T,03) = N + ( SPACING * 73 )          
         END SELECT


         !-------------!
         ! HG0-VEGT    !
         !-------------!

         T = 154
         N = 1
         SELECT CASE( T )
           CASE( 154 )
               NAME (T,03) = 'Hg0'
               FNAME(T,03) = 'Hg0'
               INDEX(T,03) = N + ( SPACING * 74 )          
         END SELECT


         !-------------!
         ! HG0-SOIL    !
         !-------------!

         T = 155
         N = 1
         SELECT CASE( T )
           CASE( 155 )
               NAME (T,03) = 'Hg0'
               FNAME(T,03) = 'Hg0'
               INDEX(T,03) = N + ( SPACING * 75 )          
         END SELECT


         !-------------!
         ! HG0-FXUP    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 10 + (5D0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 156 )
                  NAME (T,03) = 'Hg0'
                  FNAME(T,03) = 'Hg0'
                  INDEX(T,03) = N + ( SPACING * 76 ) 
              CASE( 157 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 158 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 76 ) 
              CASE( 159 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 76 )  
              CASE( 160 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 161 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 162 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 163 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 164 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 165 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 166 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 167 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 168 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 169 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 170 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 171 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 172 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 173 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 174 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 175 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 176 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 177 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 178 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 179 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 180 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 181 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 182 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 183 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE( 184 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 76 )
              CASE DEFAULT
            END SELECT
         ENDDO


         !-------------!
         ! HG0-FXDN    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 10 + (6D0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 185 )
                  NAME (T,03) = 'Hg0'
                  FNAME(T,03) = 'Hg0'
                  INDEX(T,03) = N + ( SPACING * 77 ) 
              CASE( 186 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 187 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 77 ) 
              CASE( 188 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 77 )  
              CASE( 189 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 190 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 191 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 192 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 193 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 194 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 195 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 196 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 197 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 198 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 199 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 200 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 201 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 202 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 203 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 204 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 205 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 206 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 207 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 208 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 209 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 210 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 211 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 212 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE( 213 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 77 )
              CASE DEFAULT
            END SELECT
         ENDDO


         !-------------!
         ! HG0-SNOW    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 10 + (7D0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 214 )
                  NAME (T,03) = 'Hg0'
                  FNAME(T,03) = 'Hg0'
                  INDEX(T,03) = N + ( SPACING * 78 ) 
              CASE( 215 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 216 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 78 ) 
              CASE( 217 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 78 )  
              CASE( 218 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 219 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 220 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 221 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 222 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 223 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 224 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 225 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 226 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 227 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 228 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 229 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 230 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 231 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 232 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 233 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 234 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 235 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 236 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 237 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 238 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 239 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 240 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 241 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE( 242 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 78 )
              CASE DEFAULT
            END SELECT
         ENDDO


         !-------------!
         ! HG-NETOX    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 10 + (8D0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 243 )
                  NAME (T,03) = 'Hg0'
                  FNAME(T,03) = 'Hg0'
                  INDEX(T,03) = N + ( SPACING * 79 ) 
              CASE( 244 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 245 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 79 ) 
              CASE( 246 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 79 )  
              CASE( 247 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 248 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 249 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 250 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 251 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 252 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 253 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 254 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 255 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 256 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 257 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 258 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 259 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 260 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 261 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 262 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 263 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 264 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 265 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 266 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 267 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 268 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 269 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 270 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE( 271 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 79 )
              CASE DEFAULT
            END SELECT
         ENDDO


         !-------------!
         ! HG2-OXOH    !
         !-------------!

         T = 272
         N = 1
         SELECT CASE( T )
           CASE( 272 )
               NAME (T,03) = 'Hg2'
               FNAME(T,03) = 'Hg2'
               INDEX(T,03) = N + ( SPACING * 80 )          
         END SELECT



         !-------------!
         ! HG2-OXO3    !
         !-------------!

         T = 273
         N = 1
         SELECT CASE( T )
           CASE( 273 )
               NAME (T,03) = 'Hg2'
               FNAME(T,03) = 'Hg2'
               INDEX(T,03) = N + ( SPACING * 81 )          
         END SELECT



         !-------------!
         ! HG2-SALT    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 12 + (9D0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 274 )
                  NAME (T,03) = 'Hg2'
                  FNAME(T,03) = 'Hg2'
                  INDEX(T,03) = N + ( SPACING * 82 ) 
              CASE( 275 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 276 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 82 ) 
              CASE( 277 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 82 )  
              CASE( 278 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 279 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 280 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 281 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 282 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 283 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 284 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 285 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 286 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 287 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 288 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 289 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 290 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 291 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 292 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 293 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 294 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 295 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 296 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 297 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 298 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 299 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 300 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 301 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE( 302 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 82 )
              CASE DEFAULT
            END SELECT
         ENDDO


         !-------------!
         ! HG2-SSRX    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 12 + (10D0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 303 )
                  NAME (T,03) = 'Hg2'
                  FNAME(T,03) = 'Hg2'
                  INDEX(T,03) = N + ( SPACING * 83 ) 
              CASE( 304 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 305 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 83 ) 
              CASE( 306 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 83 )  
              CASE( 307 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 308 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 309 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 310 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 311 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 312 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 313 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 314 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 315 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 316 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 317 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 318 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 319 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 320 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 321 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 322 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 323 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 324 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 325 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 326 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 327 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 328 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 329 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 330 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE( 331 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 83 )
              CASE DEFAULT
            END SELECT
         ENDDO


         !-------------!
         ! HG2-OXBR    !
         !-------------!

         ! Loop over tracers
         DO N = 1, N_Hg_CATS

         T = N + 12 + (11D0 * N_Hg_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 332 )
                  NAME (T,03) = 'Hg2'
                  FNAME(T,03) = 'Hg2'
                  INDEX(T,03) = N + ( SPACING * 84 ) 
              CASE( 333 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 334 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 84 ) 
              CASE( 335 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 84 )  
              CASE( 336 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 337 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 338 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 339 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 340 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 341 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 342 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 343 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 344 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 345 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 346 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 347 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 348 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 349 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 350 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 351 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 352 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 353 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 354 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 355 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 356 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 357 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 358 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 359 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE( 360 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 84 )
              CASE DEFAULT
            END SELECT
         ENDDO


         !-------------!
         ! PL-HG2-$    !
         !-------------!

         T = 361
         N = 1
         SELECT CASE( T )
           CASE( 361 )
               NAME (T,03) = 'Br'
               FNAME(T,03) = 'Br Concentration'
               INDEX(T,03) = N + ( SPACING * 35 )
         END SELECT


         !-------------!
         ! PL-HG2-$    !
         !-------------!

         T = 362
         N = 2
         SELECT CASE( T )
           CASE( 362 )
               NAME (T,03) = 'BrO'
               FNAME(T,03) = 'BrO Concentration'
               INDEX(T,03) = N + ( SPACING * 35 )
         END SELECT


       	 !-------------! 
         ! OCEAN-HG    !
       	 !-------------! 

         ! Loop over tracers
         DO N = 1, 87

         T = N + 14 + (12D0 * N_HG_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 363 )
                  NAME (T,03) = 'Hg0'
                  FNAME(T,03) = 'Hg0'
                  INDEX(T,03) = N + ( SPACING * 41 )  !eds 9/9/10 
              CASE( 364 )
                  NAME (T,03) = 'Hg2'
                  FNAME(T,03) = 'Hg2'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 365 )
                  NAME (T,03) = 'HgP'
                  FNAME(T,03) = 'HgP'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 366 )
                  NAME (T,03) = 'Hg0_can'
                  FNAME(T,03) = 'Hg0 Canada'
       	       	  INDEX(T,03) =	N + ( SPACING *	41 )  
              CASE( 367 )
                  NAME (T,03) = 'Hg0_usa'
                  FNAME(T,03) = 'Hg0 USA'
                  INDEX(T,03) = N + ( SPACING * 41 ) 
              CASE( 368 )
                  NAME (T,03) = 'Hg0_cam'
                  FNAME(T,03) = 'Hg0 Central America'
                  INDEX(T,03) = N + ( SPACING * 41 )  
              CASE( 369 )
                  NAME (T,03) = 'Hg0_sam'
                  FNAME(T,03) = 'Hg0 South American'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 370 )
                  NAME (T,03) = 'Hg0_waf'
                  FNAME(T,03) = 'Hg0 West Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 371 )
                  NAME (T,03) = 'Hg0_eaf'
                  FNAME(T,03) = 'Hg0 East Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 372 )
                  NAME (T,03) = 'Hg0_saf'
                  FNAME(T,03) = 'Hg0 South Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 373 )
                  NAME (T,03) = 'Hg0_naf'
                  FNAME(T,03) = 'Hg0 North Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 374 )
                  NAME (T,03) = 'Hg0_eur'
                  FNAME(T,03) = 'Hg0 OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 375 )
                  NAME (T,03) = 'Hg0_eeu'
                  FNAME(T,03) = 'Hg0 Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 376 )
                  NAME (T,03) = 'Hg0_sov'
                  FNAME(T,03) = 'Hg0 Former USSR'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 377 )
                  NAME (T,03) = 'Hg0_mde'
                  FNAME(T,03) = 'Hg0 Middle East'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 378 )
                  NAME (T,03) = 'Hg0_sas'
                  FNAME(T,03) = 'Hg0 South Asia'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 379 )
                  NAME (T,03) = 'Hg0_eas'
                  FNAME(T,03) = 'Hg0 East Asia'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 380 )
                  NAME (T,03) = 'Hg0_sea'
                  FNAME(T,03) = 'Hg0 Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 381 )
                  NAME (T,03) = 'Hg0_jpn'
                  FNAME(T,03) = 'Hg0 Japan'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 382 )
                  NAME (T,03) = 'Hg0_ocn'
                  FNAME(T,03) = 'Hg0 Oceania'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 383 )
                  NAME (T,03) = 'Hg0_so'
                  FNAME(T,03) = 'Hg0 Soil'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 384 )
                  NAME (T,03) = 'Hg0_bb'
                  FNAME(T,03) = 'Hg0 Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 385 )
                  NAME (T,03) = 'Hg0_geo'
                  FNAME(T,03) = 'Hg0 Geogenic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 386 )
                  NAME (T,03) = 'Hg0_atl'
                  FNAME(T,03) = 'Hg0 Atlantic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 387 )
                  NAME (T,03) = 'Hg0_nat'
                  FNAME(T,03) = 'Hg0 North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 388 )
                  NAME (T,03) = 'Hg0_sat'
                  FNAME(T,03) = 'Hg0 South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 389 )
                  NAME (T,03) = 'Hg0_npa'
                  FNAME(T,03) = 'Hg0 North Pacific'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 390 )
                  NAME (T,03) = 'Hg0_arc'
                  FNAME(T,03) = 'Hg0 Arctic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 391 )
                  NAME (T,03) = 'Hg0_ant'
                  FNAME(T,03) = 'Hg0 Antarctica'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 392 )
                  NAME (T,03) = 'Hg0_ocn'
                  FNAME(T,03) = 'Hg0 Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 393 )
                  NAME (T,03) = 'Hg0_str'
                  FNAME(T,03) = 'Hg0 Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 394 )
                  NAME (T,03) = 'Hg2_can'
                  FNAME(T,03) = 'Hg2 Canada'
       	       	  INDEX(T,03) =	N + ( SPACING *	41 )  
              CASE( 395 )
                  NAME (T,03) = 'Hg2_usa'
                  FNAME(T,03) = 'Hg2 USA'
                  INDEX(T,03) = N + ( SPACING * 41 ) 
              CASE( 396 )
                  NAME (T,03) = 'Hg2_cam'
                  FNAME(T,03) = 'Hg2 Central America'
                  INDEX(T,03) = N + ( SPACING * 41 )  
              CASE( 397 )
                  NAME (T,03) = 'Hg2_sam'
                  FNAME(T,03) = 'Hg2 South American'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 398 )
                  NAME (T,03) = 'Hg2_waf'
                  FNAME(T,03) = 'Hg2 West Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 399 )
                  NAME (T,03) = 'Hg2_eaf'
                  FNAME(T,03) = 'Hg2 East Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 400 )
                  NAME (T,03) = 'Hg2_saf'
                  FNAME(T,03) = 'Hg2 South Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 401 )
                  NAME (T,03) = 'Hg2_naf'
                  FNAME(T,03) = 'Hg2 North Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 402 )
                  NAME (T,03) = 'Hg2_eur'
                  FNAME(T,03) = 'Hg2 OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 403 )
                  NAME (T,03) = 'Hg2_eeu'
                  FNAME(T,03) = 'Hg2 Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 404 )
                  NAME (T,03) = 'Hg2_sov'
                  FNAME(T,03) = 'Hg2 Former USSR'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 405 )
                  NAME (T,03) = 'Hg2_mde'
                  FNAME(T,03) = 'Hg2 Middle East'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 406 )
                  NAME (T,03) = 'Hg2_sas'
                  FNAME(T,03) = 'Hg2 South Asia'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 407 )
                  NAME (T,03) = 'Hg2_eas'
                  FNAME(T,03) = 'Hg2 East Asia'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 408 )
                  NAME (T,03) = 'Hg2_sea'
                  FNAME(T,03) = 'Hg2 Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 409 )
                  NAME (T,03) = 'Hg2_jpn'
                  FNAME(T,03) = 'Hg2 Japan'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 410 )
                  NAME (T,03) = 'Hg2_ocn'
                  FNAME(T,03) = 'Hg2 Oceania'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 411 )
                  NAME (T,03) = 'Hg2_so'
                  FNAME(T,03) = 'Hg2 Soil'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 412 )
                  NAME (T,03) = 'Hg2_bb'
                  FNAME(T,03) = 'Hg2 Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 413 )
                  NAME (T,03) = 'Hg2_geo'
                  FNAME(T,03) = 'Hg2 Geogenic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 414 )
                  NAME (T,03) = 'Hg2_atl'
                  FNAME(T,03) = 'Hg2 Atlantic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 415 )
                  NAME (T,03) = 'Hg2_nat'
                  FNAME(T,03) = 'Hg2 North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 416 )
                  NAME (T,03) = 'Hg2_sat'
                  FNAME(T,03) = 'Hg2 South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 417 )
                  NAME (T,03) = 'Hg2_npa'
                  FNAME(T,03) = 'Hg2 North Pacific'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 418 )
                  NAME (T,03) = 'Hg2_arc'
                  FNAME(T,03) = 'Hg2 Arctic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 419 )
                  NAME (T,03) = 'Hg2_ant'
                  FNAME(T,03) = 'Hg2 Antarctica'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 420 )
                  NAME (T,03) = 'Hg2_ocn'
                  FNAME(T,03) = 'Hg2 Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 421 )
                  NAME (T,03) = 'Hg2_str'
                  FNAME(T,03) = 'Hg2 Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 422 )
                  NAME (T,03) = 'HgP_can'
                  FNAME(T,03) = 'HgP Canada'
       	       	  INDEX(T,03) =	N + ( SPACING *	41 )  
              CASE( 423 )
                  NAME (T,03) = 'HgP_usa'
                  FNAME(T,03) = 'HgP USA'
                  INDEX(T,03) = N + ( SPACING * 41 ) 
              CASE( 424 )
                  NAME (T,03) = 'HgP_cam'
                  FNAME(T,03) = 'HgP Central America'
                  INDEX(T,03) = N + ( SPACING * 41 )  
              CASE( 425 )
                  NAME (T,03) = 'HgP_sam'
                  FNAME(T,03) = 'HgP South American'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 426 )
                  NAME (T,03) = 'HgP_waf'
                  FNAME(T,03) = 'HgP West Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 427 )
                  NAME (T,03) = 'HgP_eaf'
                  FNAME(T,03) = 'HgP East Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 428 )
                  NAME (T,03) = 'HgP_saf'
                  FNAME(T,03) = 'HgP South Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 429 )
                  NAME (T,03) = 'HgP_naf'
                  FNAME(T,03) = 'HgP North Africa'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 430 )
                  NAME (T,03) = 'HgP_eur'
                  FNAME(T,03) = 'HgP OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 431 )
                  NAME (T,03) = 'HgP_eeu'
                  FNAME(T,03) = 'HgP Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 432 )
                  NAME (T,03) = 'HgP_sov'
                  FNAME(T,03) = 'HgP Former USSR'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 433 )
                  NAME (T,03) = 'HgP_mde'
                  FNAME(T,03) = 'HgP Middle East'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 434 )
                  NAME (T,03) = 'HgP_sas'
                  FNAME(T,03) = 'HgP South Asia'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 435 )
                  NAME (T,03) = 'HgP_eas'
                  FNAME(T,03) = 'HgP East Asia'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 436 )
                  NAME (T,03) = 'HgP_sea'
                  FNAME(T,03) = 'HgP Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 437 )
                  NAME (T,03) = 'HgP_jpn'
                  FNAME(T,03) = 'HgP Japan'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 438 )
                  NAME (T,03) = 'HgP_ocn'
                  FNAME(T,03) = 'HgP Oceania'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 439 )
                  NAME (T,03) = 'HgP_so'
                  FNAME(T,03) = 'HgP Soil'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 440 )
                  NAME (T,03) = 'HgP_bb'
                  FNAME(T,03) = 'HgP Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 441 )
                  NAME (T,03) = 'HgP_geo'
                  FNAME(T,03) = 'HgP Geogenic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 442 )
                  NAME (T,03) = 'HgP_atl'
                  FNAME(T,03) = 'HgP Atlantic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 443 )
                  NAME (T,03) = 'HgP_nat'
                  FNAME(T,03) = 'HgP North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 444 )
                  NAME (T,03) = 'HgP_sat'
                  FNAME(T,03) = 'HgP South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 445 )
                  NAME (T,03) = 'HgP_npa'
                  FNAME(T,03) = 'HgP North Pacific'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 446 )
                  NAME (T,03) = 'HgP_arc'
                  FNAME(T,03) = 'HgP Arctic'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 447 )
                  NAME (T,03) = 'HgP_ant'
                  FNAME(T,03) = 'HgP Antarctica'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 448 )
                  NAME (T,03) = 'HgP_ocn'
                  FNAME(T,03) = 'HgP Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE( 449 )
                  NAME (T,03) = 'HgP_str'
                  FNAME(T,03) = 'HgP Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 41 )
              CASE DEFAULT
            END SELECT
         ENDDO


       	 !-------------! 
         ! SNOW-HG     !
       	 !-------------! 

         ! Loop over tracers
         DO N = 1, 31

         T = N + 14 + (15D0 * N_HG_CATS)

            ! Get name and long name for each tracer
            SELECT CASE( T )
              CASE( 450 )
                  NAME (T,03) = 'Hg0'
                  FNAME(T,03) = 'Hg0'
                  INDEX(T,03) = N + ( SPACING * 49 )  !eds 9/9/10 
              CASE( 451 )
                  NAME (T,03) = 'Hg2'
                  FNAME(T,03) = 'Hg2'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 452 )
                  NAME (T,03) = 'HgP'
                  FNAME(T,03) = 'HgP'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 453 )
                  NAME (T,03) = 'Hg_can'
                  FNAME(T,03) = 'Hg Canada'
       	       	  INDEX(T,03) =	N + ( SPACING *	49 )  
              CASE( 454 )
                  NAME (T,03) = 'Hg_usa'
                  FNAME(T,03) = 'Hg USA'
                  INDEX(T,03) = N + ( SPACING * 49 ) 
              CASE( 455 )
                  NAME (T,03) = 'Hg_cam'
                  FNAME(T,03) = 'Hg Central America'
                  INDEX(T,03) = N + ( SPACING * 49 )  
              CASE( 456 )
                  NAME (T,03) = 'Hg_sam'
                  FNAME(T,03) = 'Hg South American'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 457 )
                  NAME (T,03) = 'Hg_waf'
                  FNAME(T,03) = 'Hg West Africa'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 458 )
                  NAME (T,03) = 'Hg_eaf'
                  FNAME(T,03) = 'Hg East Africa'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 459 )
                  NAME (T,03) = 'Hg_saf'
                  FNAME(T,03) = 'Hg South Africa'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 460 )
                  NAME (T,03) = 'Hg_naf'
                  FNAME(T,03) = 'Hg North Africa'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 461 )
                  NAME (T,03) = 'Hg_eur'
                  FNAME(T,03) = 'Hg OECD Europe'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 462 )
                  NAME (T,03) = 'Hg_eeu'
                  FNAME(T,03) = 'Hg Eastern Europe'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 463 )
                  NAME (T,03) = 'Hg_sov'
                  FNAME(T,03) = 'Hg Former USSR'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 464 )
                  NAME (T,03) = 'Hg_mde'
                  FNAME(T,03) = 'Hg Middle East'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 465 )
                  NAME (T,03) = 'Hg_sas'
                  FNAME(T,03) = 'Hg South Asia'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 466 )
                  NAME (T,03) = 'Hg_eas'
                  FNAME(T,03) = 'Hg East Asia'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 467 )
                  NAME (T,03) = 'Hg_sea'
                  FNAME(T,03) = 'Hg Southeast Asia'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 468 )
                  NAME (T,03) = 'Hg_jpn'
                  FNAME(T,03) = 'Hg Japan'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 469 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Oceania'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 470 )
                  NAME (T,03) = 'Hg_so'
                  FNAME(T,03) = 'Hg Soil'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 471 )
                  NAME (T,03) = 'Hg_bb'
                  FNAME(T,03) = 'Hg Biomass Burning'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 472 )
                  NAME (T,03) = 'Hg_geo'
                  FNAME(T,03) = 'Hg Geogenic'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 473 )
                  NAME (T,03) = 'Hg_atl'
                  FNAME(T,03) = 'Hg Atlantic'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 474 )
                  NAME (T,03) = 'Hg_nat'
                  FNAME(T,03) = 'Hg North Atlantic'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 475 )
                  NAME (T,03) = 'Hg_sat'
                  FNAME(T,03) = 'Hg South Atlantic'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 476 )
                  NAME (T,03) = 'Hg_npa'
                  FNAME(T,03) = 'Hg North Pacific'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 477 )
                  NAME (T,03) = 'Hg_arc'
                  FNAME(T,03) = 'Hg Arctic'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 478 )
                  NAME (T,03) = 'Hg_ant'
                  FNAME(T,03) = 'Hg Antarctica'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 479 )
                  NAME (T,03) = 'Hg_ocn'
                  FNAME(T,03) = 'Hg Other Ocean'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE( 480 )
                  NAME (T,03) = 'Hg_str'
                  FNAME(T,03) = 'Hg Residual Stratosphere'
                  INDEX(T,03) = N + ( SPACING * 49 )
              CASE DEFAULT
            END SELECT
         ENDDO

      ENDIF ! LSPLIT


      !-------------------------------------
      ! CO2 sources & fluxes
      !-------------------------------------
      IF ( ND04 > 0 ) THEN
         
         ! Number of tracers
         ! Increased from 6 to 10. (rnassar, 30/08/10)
!         NTRAC(04) = 6
         NTRAC(04) = 10

         ! Loop over tracers: CO2-SRCE
         DO T = 1, NTRAC(04)

            ! Define quantities
            UNIT (T,04) = 'atoms C/cm2/s'
            INDEX(T,04) = T + ( SPACING * 40 )
            MOLC (T,04) = 1
            MWT  (T,04) = 12e-3
            SCALE(T,04) = 1e0

            ! Get name, long-name
            ! Names updated (rnassar, 30/08/10)
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,04) = 'CO2ff'
                  FNAME(T,04) = 'CO2 fossil fuel emiss'
               CASE( 2  )
!                  NAME (T,04) = 'CO2ocn'
!                  FNAME(T,04) = 'CO2 ocean exchange'
                  NAME (T,04) = 'CO2oc'
                  FNAME(T,04) = 'CO2 ocean emissions'
               CASE( 3  )
                  NAME (T,04) = 'CO2bal'
                  FNAME(T,04) = 'CO2 balanced biosphere'
               CASE( 4  )
                  NAME (T,04) = 'CO2bb'
                  FNAME(T,04) = 'CO2 biomass burning emiss'
               CASE( 5  )
                  NAME (T,04) = 'CO2bf'
                  FNAME(T,04) = 'CO2 biofuel emission'
               CASE( 6  )
!                  NAME (T,04) = 'CO2net'
                  NAME (T,04) = 'CO2nte'
                  FNAME(T,04) = 'CO2 net terr exchange'
               CASE( 7  )
                  NAME (T,04) = 'CO2shp'
                  FNAME(T,04) = 'CO2 ship emissions'
               CASE( 8  )
                  NAME (T,04) = 'CO2pln'
                  FNAME(T,04) = 'CO2 aircraft emissions'
               CASE( 9  )
                  NAME (T,04) = 'CO2che'
                  FNAME(T,04) = 'CO2 chemical oxidation'
               CASE( 10  )
                  NAME (T,04) = 'CO2sur'
                  FNAME(T,04) = 'CO2 chem surf correction'
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO
      ENDIF

      !-------------------------------------
      ! Sulfur production & loss (ND05)
      !-------------------------------------
      IF ( ND05 > 0 ) THEN

         ! Number of tracers
         NTRAC(05) = 10

         ! Loop over tracers
         DO T = 1, NTRAC(05)

            ! Get name, long-name, unit, molecular weight
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,05) = 'SO2dms'
                  FNAME(T,05) = 'P(SO2) from DMS+OH'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 2  )
                  NAME (T,05) = 'SO2no3'
                  FNAME(T,05) = 'P(SO2) from DMS+NO3'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 3  )
                  NAME (T,05) = 'SO2tot'
                  FNAME(T,05) = 'Total P(SO2)'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 4  )
                  NAME (T,05) = 'MSAdms'
                  FNAME(T,05) = 'P(MSA) from DMS'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 5  )
                  NAME (T,05) = 'SO4gas'
                  FNAME(T,05) = 'P(SO4) by gas phase'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 6  )
                  NAME (T,05) = 'SO4h2o2'
                  FNAME(T,05) = 'P(SO4) by cloud H2O2'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 7  )
                  NAME (T,05) = 'SO4o3'
                  FNAME(T,05) = 'P(SO4) by cloud O3'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 8  )
                  NAME (T,05) = 'SO4ss'
                  FNAME(T,05) = 'P(SO4) by seasalt O3'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 9  )
                  NAME (T,05) = 'LOH'
                  FNAME(T,05) = 'L(OH) by DMS'
                  UNIT (T,05) = 'kg OH'
                  MWT  (T,05) = 17e-3
               CASE( 10 )
                  NAME (T,05) = 'LNO3'
                  FNAME(T,05) = 'L(NO3) by DMS'
                  UNIT (T,05) = 'kg NO3'
                  MWT  (T,05) = 62e-3
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            INDEX(T,05) = T + ( SPACING * 18 )
            MOLC (T,05) = 1
            SCALE(T,05) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------
      ! Carbon aerosol emissions (ND07)
      !-------------------------------------
      IF ( ND07 > 0 ) THEN


         ! Number of tracers
         IF ( Input_Opt%LSOA ) THEN
            ! hotp updated for SOA (7/31/08)
            NTRAC(07) = 15
         ELSE
            NTRAC(07) = 2
         ENDIF

         ! Loop over tracers
         DO T = 1, NTRAC(07)

            ! Define the default quantities
            UNIT (T,07) = 'kg'
            MOLC (T,07) = 1
            MWT  (T,07) = 12e-3
            SCALE(T,07) = 1e0

            ! Get name, long-name, tracern umber
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,07) = 'BLKC'
                  FNAME(T,07) = 'Black (Elemental) Carbon'
                  INDEX(T,07) = IDTBCPI + ( SPACING * 33 )
                  UNIT (T,07) = 'kgC' ! SOAupdate (hotp 7/29/10)
               CASE( 2 )
                  NAME (T,07) = 'ORGC'
                  FNAME(T,07) = 'Organic Carbon'
                  UNIT (T,07) = 'kgC' ! SOAupdate (hotp 7/29/10)
                  ! SOAupdate: semivolpoa: replace OCPI with POA (hotp 2/17/09)
                  IF ( IDTOCPI > 0 )
     &               INDEX(T,07) = IDTOCPI + ( SPACING * 33 )
                  IF ( IDTPOA1  > 0 )
     &               INDEX(T,07) = IDTPOA1 + ( SPACING * 33 )
               CASE( 3 ) ! SOAupdate (hotp 7/29/10)
                  NAME (T,07) = 'MTPA'
                  FNAME(T,07) = 'lumped ALPH'
                  INDEX(T,07) = IDTMTPA + ( SPACING * 33 )
                  MWT  (T,07) = 136e-3
               CASE( 4 )
                  NAME (T,07) = 'LIMO'
                  FNAME(T,07) = 'Limonene'
                  MWT  (T,07) = 136e-3
                  INDEX(T,07) = IDTLIMO + ( SPACING * 33 )
               CASE( 5 ) ! SOAupdate (hotp 7/29/10)
                  NAME (T,07) = 'MTPO'
                  FNAME(T,07) = 'lumped other mtp'
                  INDEX(T,07) = IDTMTPO + ( SPACING * 33 )
                  MWT  (T,07) = 136e-3
               CASE( 6 ) ! SOAupdate hotp
                  NAME (T,07) = 'SESQ'
                  FNAME(T,07) = 'Sesqterpene'
                  MWT  (T,07) = 204e-3
                  INDEX(T,07) = IDTLIMO + 3 + ( SPACING * 33 )
               CASE( 7 ) ! SOAupdate hotp
                  NAME (T,07) = 'TSOA'
                  FNAME(T,07) = 'mono+sesq aerosol'
                  INDEX(T,07) = IDTTSOA1 + ( SPACING * 33 )
               CASE( 8 ) ! SOAupdate hotp
                  NAME (T,07) = 'ISOA'
                  FNAME(T,07) = 'isoprene aerosol'
                  INDEX(T,07) = IDTISOA1 + ( SPACING * 33 )
               CASE( 9 ) ! SOAupdate hotp
                  NAME (T,07) = 'ASOA'
                  FNAME(T,07) = 'BTX + IVOC aerosol'
                  INDEX(T,07) = IDTASOA1 + ( SPACING * 33 )
               CASE( 10 ) ! SOAupdate update 5/27/10 hotp
                  NAME (T,07) = 'OPOA'
                  FNAME(T,07) = 'OPOA production'
                  IF ( IDTOPOA1 > 0 ) THEN ! hotp 6/8/10
                     INDEX(T,07) = IDTOPOA1 + ( SPACING * 33 )
                  ELSE
                     INDEX(T,07) = 80 + (SPACING * 33 )
                  ENDIF
                  UNIT (T,07) = 'kgC'
               CASE( 11 ) ! SOAupdate update 5/27/10 hotp
                  NAME (T,07) = 'OPOG'
                  FNAME(T,07) = 'POG + OH reaction'
                  UNIT (T,07) = 'kgC'
                  IF ( IDTOPOA1 > 0 ) THEN ! hotp 6/8/10
                     INDEX(T,07) = IDTOPOG1 + ( SPACING * 33 )
                  ELSE
                     INDEX(T,07) = 81 + (SPACING * 33 )
                  ENDIF
               CASE( 12 )
                  NAME (T,07) = 'SOAG'
                  FNAME(T,07) = 'SOAG production in aq. aerosol'
                  MWT  (T,07) = 58e-3
                  INDEX(T,07) = IDTSOAG + ( SPACING * 33 )
               CASE( 13 )
                  NAME (T,07) = 'SOAM'
                  FNAME(T,07) = 'SOAM production in aq. aerosol'
                  MWT  (T,07) = 72e-3
                  INDEX(T,07) = IDTSOAM + ( SPACING * 33 )
               CASE( 14 )
                  NAME (T,07) = 'SOAG'
                  FNAME(T,07) = 'SOAG production in clouds'
                  MWT  (T,07) = 58e-3
                  INDEX(T,07) = 91 + ( SPACING * 33 )
               CASE( 15 )
                  NAME (T,07) = 'SOAM'
                  FNAME(T,07) = 'SOAM production in clouds'
                  MWT  (T,07) = 72e-3
                  INDEX(T,07) = 92 + ( SPACING * 33 )

               CASE DEFAULT
                  ! Nothing
            END SELECT

         ENDDO
      ENDIF

      !-------------------------------------
      ! HCN & CH3CN sources & sinks (ND09)
      !-------------------------------------
      IF ( ND09 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(09) = Input_Opt%N_TRACERS + 6

         ! Loop over tracers: HCN-PL-$
         DO T = 1, NTRAC(09)-6
            NAME (T,09) = Input_Opt%TRACER_NAME(T)
            FNAME(T,09) = TRIM( Input_Opt%TRACER_NAME(T) ) 
     &                    // ' lost via OH rxn'
            UNIT (T,09) = 'kg'
            INDEX(T,09) = T + ( SPACING * 38 )
            MOLC (T,09) = 1
            MWT  (T,09) = Input_Opt%TRACER_MW_KG(T)
            SCALE(T,09) = 1e0
         ENDDO

         ! Loop over tracers: HCN-SRCE
         DO T = NTRAC(09)-6+1, NTRAC(09)

            ! Get tracer number, C number, scale factor, unit
            INDEX(T,09) = ( T - Input_Opt%N_TRACERS ) + ( SPACING * 39 )
            MOLC (T,09) = 1
            SCALE(T,09) = 1e0
            UNIT (T,09) = 'molec/cm2/s'
              
            ! Get name, longname, unit, mwt
            IF ( T == Input_Opt%N_TRACERS+1 ) THEN
               NAME (T,09) = 'HCNbb'
               FNAME(T,09) = 'HCN biomass emissions'
               MWT  (T,09) = 27d-3
            ELSE IF ( T == Input_Opt%N_TRACERS+2 ) THEN
               NAME (T,09) = 'CH3CNbb'
               FNAME(T,09) = 'CH3CN biomass emissions'
               MWT  (T,09) = 41d-3
            ELSE IF ( T == Input_Opt%N_TRACERS+3 ) THEN
               NAME (T,09) = 'HCNdf'
               FNAME(T,09) = 'HCN domestic FF emiss'
               MWT  (T,09) = 27d-3
            ELSE IF ( T == Input_Opt%N_TRACERS+4 ) THEN
               NAME (T,09) = 'CH3CNdf'
               FNAME(T,09) = 'CH3CN domestic FF emiss'
               MWT  (T,09) = 41d-3
            ELSE IF ( T == Input_Opt%N_TRACERS+5 ) THEN
               NAME (T,09) = 'HCNoc'
               FNAME(T,09) = 'HCN ocean uptake'
               MWT  (T,09) = 27d-3
            ELSE IF ( T == Input_Opt%N_TRACERS+6 ) THEN
               NAME (T,09) = 'CH3CNoc'
               FNAME(T,09) = 'CH3CN ocean uptake'
               MWT  (T,09) = 41d-3
            ENDIF
         ENDDO
      ENDIF
    
      !-------------------------------------      
      ! H2-HD source, production & loss (ND10)
      !-------------------------------------      
      IF ( ND10 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(10) = 20

         ! Loop over tracers 
         DO T = 1, NTRAC(10)

            ! Define quantities
            UNIT (T,10) = 'molec/cm3/s' ! overwrite below when needed
            MOLC (T,10) = 1   
            SCALE(T,10) = 1e0

            ! Get name, long-name, index, and new units
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,10) = 'H2oh'
                  FNAME(T,10) = 'Loss of H2 by OH'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 2  )
                  NAME (T,10) = 'H2iso'
                  FNAME(T,10) = 'Prod of H2 from Isoprene'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 3  )
                  NAME (T,10) = 'H2ch4'
                  FNAME(T,10) = 'Prod of H2 from CH4'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 4  )
                  NAME (T,10) = 'H2ch3oh'
                  FNAME(T,10) = 'Prod of H2 from CH3OH'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 5  )
                  NAME (T,10) = 'H2mono'
                  FNAME(T,10) = 'Prod of H2 from MONO'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 6  )
                  NAME (T,10) = 'H2acet'
                  FNAME(T,10) = 'Prod of H2 from ACET'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 7  )
                  NAME (T,10) = 'H2o1d'
                  FNAME(T,10) = 'Loss of H2 by strat O1D'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 8  )
                  NAME (T,10) = 'HDoh'
                  FNAME(T,10) = 'Loss of HD by OH'
                  UNIT (T,10) = 's-1'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 9  )
                  NAME (T,10) = 'HDiso'
                  FNAME(T,10) = 'Prod of HD from Isoprene'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 10 )
                  NAME (T,10) = 'HDch4'
                  FNAME(T,10) = 'Prod of HD from CH4'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 11 )
                  NAME (T,10) = 'HDch3oh'
                  FNAME(T,10) = 'Prod of HD from CH3OH'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 12 )
                  NAME (T,10) = 'HDmono'
                  FNAME(T,10) = 'Prod of HD from MONO'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 13 )
                  NAME (T,10) = 'HDacet'
                  FNAME(T,10) = 'Prod of HD from ACET'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 14 )
                  NAME (T,10) = 'HDo1d'
                  FNAME(T,10) = 'Loss of HD by strat O1D'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 15 )
                  NAME (T,10) = 'Alpha'
                  FNAME(T,10) = 'Ratio of OH k rates kHD/kH2'
                  UNIT (T,10) = 'unitless'
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 16 )
                  NAME (T,10) = 'H2anth'
                  FNAME(T,10) = 'Fossil Fuel H2 (anthropogenic)'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )
               CASE( 17 )
                  NAME (T,10) = 'H2bb'
                  FNAME(T,10) = 'Biomass Burning of H2'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )

               CASE( 18 )
                  NAME (T,10) = 'H2bf'
                  FNAME(T,10) = 'Biofuel Burning of H2'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )

               CASE( 19 )
                  NAME (T,10) = 'H2ocean'
                  FNAME(T,10) = 'Ocean emissions of H2'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )
               CASE( 20 )
                  NAME (T,10) = 'HDocean'
                  FNAME(T,10) = 'Ocean emissions of HD'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )

               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO
      ENDIF

      !-------------------------------------
      ! Acetone sources & sinks (ND11)
      !-------------------------------------
      IF ( ND11 > 0 ) THEN

         ! Number of tracers
         NTRAC(11) = 5

         ! Loop over tracers
         DO T = 1, NTRAC(11)

            ! Get name and long-name
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,11) = 'ACETmo'
                  FNAME(T,11) = 'ACET from monoterpenes'
               CASE( 2 )
                  NAME (T,11) = 'ACETmb'
                  FNAME(T,11) = 'ACET from methyl butenol'
               CASE( 3 )
                  NAME (T,11) = 'ACETbg'
                  FNAME(T,11) = 'ACET from direct emissions'
               CASE( 4 )
                  NAME (T,11) = 'ACETop'
                  FNAME(T,11) = 'ACET from ocean source'
               CASE( 5 )
                  NAME (T,11) = 'ACETol'
                  FNAME(T,11) = 'ACET from ocean sink'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,11) = 'atoms C/cm2/s'
            INDEX(T,11) = T + ( SPACING * 22 )
            MOLC (T,11) = 3
            MWT  (T,11) = 12e-3
            SCALE(T,11) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Distrib of emissions in PBL (ND12)
      !-------------------------------------      
      IF ( ND12 > 0 ) THEN
         T           = 1
         NTRAC(12)   = T
         NAME (T,12) = 'BLfrac'
         FNAME(T,12) = 'Boundary layer fraction'
         UNIT (T,12) = 'unitless'
         INDEX(T,12) = T + ( SPACING * 23 )
         MOLC (T,12) = 1
         MWT  (T,12) = 0e0
         SCALE(T,12) = 1e0
      ENDIF

      !-------------------------------------      
      ! Fraction of grid box experiencing
      ! LS or conv precipitation (ND16)
      !-------------------------------------   
      IF ( ND16 > 0 ) THEN

         ! Number of tracers
         NTRAC(16) = 2

         ! Loop over tracers
         DO T = 1, NTRAC(16)

            ! Get name, longname
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,16) = 'LS-f'
                  FNAME(T,16) = 'Large scale fraction'
               CASE( 2 )
                  NAME (T,16) = 'CV-f'
                  FNAME(T,16) = 'Convective fraction' 
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,16) = 'unitless'
            INDEX(T,16) = T + ( SPACING * 28 )
            MOLC (T,16) = 1
            MWT  (T,16) = 0e0
            SCALE(T,16) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! F/f ratios for LS & convective
      ! rainout (ND17) & washout (ND18)
      !-------------------------------------
      IF ( ND17 + ND18 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(17) = GET_WETDEP_NSOL()

         ! Loop over tracers
         DO N = 1, NTRAC(17)

            ! Tracer number for Nth wetdep species
            T = GET_WETDEP_IDWETD(N)
          
            ! Define quantities
            NAME (N,17) = 'Ff_' // TRIM( Input_Opt%TRACER_NAME(T) )
            FNAME(N,17) = 'F/f ratio for ' // 
     &                     TRIM( Input_Opt%TRACER_NAME(T) )
            UNIT (N,17) = 'unitless'
            INDEX(N,17) = T + ( SPACING * 29 )
            MWT  (N,17) = Input_Opt%TRACER_MW_KG(T)
            MOLC (N,17) = INT( Input_Opt%TRACER_COEFF(T,1) )
            SCALE(N,17) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! CH4 Loss due to reaction w/ OH
      !-------------------------------------
      IF ( ND19 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(19) = 1
         T = 1
          
         ! Define quantities
         NAME (T,19) = 'CH4Loss'
         FNAME(T,19) = 'CH4 Loss from Reaction with OH'
         UNIT (T,19) = 'kg'
         INDEX(T,19) = T + ( SPACING * 32 )
         MWT  (T,19) = 1.6e-3
         MOLC (T,19) = 1
         SCALE(T,19) = 1e0

      ENDIF

      !-------------------------------------      
      ! Optical depths (ND21)
      !-------------------------------------
      IF ( ND21 > 0 .or. DO_TIMESERIES ) THEN

         ! Number of tracers
         NTRAC(21) = 27

         ! Loop over tracers
         DO T = 1, NTRAC(21)

            ! Define quantities
            UNIT (T,21) = 'unitless'
            INDEX(T,21) = T + ( SPACING * 14 )
            MOLC (T,21) = 1
            MWT  (T,21) = 0e0
            SCALE(T,21) = 1e0

            ! Get name long-name (and sometimes, unit)
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,21) = 'OPTD'
                  FNAME(T,21) = 'Cloud optical depth'
               CASE( 2  )
                  ! GEOS-3, GEOS-4: CLDTOT
                  NAME (T,21) = 'CLDTOT'
                  FNAME(T,21) = '3-D cloud frc'
               CASE( 4  )
                  NAME (T,21) = 'OPD'
                  FNAME(T,21) = 'Mineral dust opt depth'
               CASE( 5  )
                  NAME (T,21) = 'SD'
                  FNAME(T,21) = 'Mineral dust surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 6  )
                  NAME (T,21) = 'OPSO4'
                  FNAME(T,21) = 'Sulfate opt depth (550 nm)'
               CASE( 7  )
                  NAME (T,21) = 'HGSO4'
                  FNAME(T,21) = 'Hygr growth of SO4'
               CASE( 8  )
                  NAME (T,21) = 'SSO4'
                  FNAME(T,21) = 'Sulfate surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 9  )
                  NAME (T,21) = 'OPBC'
                  FNAME(T,21) = 'Black carbon opt depth (550 nm)'
               CASE( 10 )
                  NAME (T,21) = 'HGBC'
                  FNAME(T,21) = 'Hygr growth of BC'
               CASE( 11 )
                  NAME (T,21) = 'SBC'
                  FNAME(T,21) = 'BC surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 12 )
                  NAME (T,21) = 'OPOC'
                  FNAME(T,21) = 'Org carbon opt depth (550 nm)'
               CASE( 13 )
                  NAME (T,21) = 'HGOC'
                  FNAME(T,21) = 'Hygr growth of OC'
               CASE( 14 )
                  NAME (T,21) = 'SOC'
                  FNAME(T,21) = 'OC surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 15 )
                  NAME (T,21) = 'OPSSa'
                  FNAME(T,21) = 'Seasalt (accum) opt depth (550 nm)'
               CASE( 16 ) 
                  NAME (T,21) = 'HGSSa'
                  FNAME(T,21) = 'Hygr growth of seasalt (accum)'
               CASE( 17 )
                  NAME (T,21) = 'SSSa'
                  FNAME(T,21) = 'Seasalt (accum) surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 18 )
                  NAME (T,21) = 'OPSSc'
                  FNAME(T,21) = 'Seasalt (coarse) opt depth (550 nm)'
               CASE( 19 ) 
                  NAME (T,21) = 'HGSSc'
                  FNAME(T,21) = 'Hygr growth of seasalt (coarse)'
               CASE( 20 )
                  NAME (T,21) = 'SSSc'
                  FNAME(T,21) = 'Seasalt (coarse) surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 21 )
                  NAME (T,21) = 'OPD1'
                  FNAME(T,21) = 'Dust bin 1 AOD (550 nm)'
               CASE( 22 )
                  NAME (T,21) = 'OPD2'
                  FNAME(T,21) = 'Dust bin 2 AOD (550 nm)'
               CASE( 23 )
                  NAME (T,21) = 'OPD3'
                  FNAME(T,21) = 'Dust bin 3 AOD (550 nm)'
               CASE( 24 )
                  NAME (T,21) = 'OPD4'
                  FNAME(T,21) = 'Dust bin 4 AOD (550 nm)'
               CASE( 25 )
                  NAME (T,21) = 'OPD5'
                  FNAME(T,21) = 'Dust bin 5 AOD (550 nm)'
               CASE( 26 )
                  NAME (T,21) = 'OPD6'
                  FNAME(T,21) = 'Dust bin 6 AOD (550 nm)'
               CASE( 27 )
                  NAME (T,21) = 'OPD7'
                  FNAME(T,21) = 'Dust bin 7 AOD (550 nm)'
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Optical depths (ND22)
      !-------------------------------------
      IF ( ND22 > 0 ) THEN

         IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN 

            !------------------
            ! Full-chem run
            !------------------

            ! Number of tracers
            NTRAC(22) = 14 ! Increase for bromine (jpp, mpayer, 12/28/11)

            ! Loop over tracers
            DO T = 1, NTRAC(22)

               ! Get name
               SELECT CASE( T )
                  CASE( 1  )
                     NAME(T,22) = 'JNO2'
                  CASE( 2  )
                     NAME(T,22) = 'JHNO3'
                  CASE( 3  )
                     NAME(T,22) = 'JH2O2'
                  CASE( 4  )
                     NAME(T,22) = 'JCH2O'
                  CASE( 5  )
                     NAME(T,22) = 'JO3'
                  CASE( 6  )
                     NAME(T,22) = 'POH'
                  CASE( 7  )
                     NAME(T,22) = 'JGLYX'
                  CASE( 8  )
                     NAME(T,22) = 'JMGLY'
                  CASE( 9  )
                     NAME(T,22) = 'JBrO'
                  CASE( 10 )
                     NAME(T,22) = 'JHOBr'
                  CASE( 11 )
                     NAME(T,22) = 'JBrNO2'
                  CASE( 12 )
                     NAME(T,22) = 'JBrNO3'
                  CASE( 13 )
                     NAME(T,22) = 'JCHBr3'
                  CASE( 14 )
                     NAME(T,22) = 'JBr2'
                 CASE DEFAULT
                  ! Nothing
               END SELECT

               ! Define quantities
               FNAME(T,22) = TRIM( NAME(T,22) ) // ' photolysis rate'
               UNIT (T,22) = 's-1'
               INDEX(T,22) = T + ( SPACING * 13 )
               MOLC (T,22) = 1
               MWT  (T,22) = 0e0
               SCALE(T,22) = 1e0
            ENDDO
         
         ELSE IF ( Input_Opt%ITS_A_CH3I_SIM ) THEN

            !------------------
            ! CH3I run
            !------------------
            NTRAC(22)   = 1 
            NAME (T,22) = 'JCH3I'
            FNAME(T,22) = TRIM( NAME(T,22) ) // ' photolysis rate'
            UNIT (T,22) = 's-1'
            INDEX(T,22) = T + ( SPACING * 13 )
            MOLC (T,22) = 1
            MWT  (T,22) = 0e0
            SCALE(T,22) = 1e0            

         ELSE IF ( Input_Opt%ITS_A_HCN_SIM ) THEN

            !------------------
            ! HCN run
            !------------------
            NTRAC(22)   = 1 
            NAME (T,22) = 'JHCN'
            FNAME(T,22) = TRIM( NAME(T,22) ) // ' photolysis rate'
            UNIT (T,22) = 's-1'
            INDEX(T,22) = T + ( SPACING * 13 )
            MOLC (T,22) = 1
            MWT  (T,22) = 0e0
            SCALE(T,22) = 1e0    

         ENDIF
      ENDIF

      !-------------------------------------      
      ! Biomass emissions (ND28)
      !-------------------------------------
      IF ( ND28 > 0 ) THEN

         IF ( Input_Opt%ITS_A_CO2_SIM ) THEN

            !---------------------------
            ! CO2 simulation only
            !---------------------------

            ! Number of tracers
            NTRAC(28) = 1
            
            ! Define quantities
            NAME (T,28) = 'CO2'
            FNAME(T,28) = TRIM( NAME(T,28) ) // ' biomass'
            INDEX(T,28) = 1 + ( SPACING * 45 )
            MWT  (T,28) = 44e-3
            MOLC (T,28) = 1
            UNIT (T,28) = 'molec/cm2/s'
            SCALE(T,28) = 1e0

         ELSE

            !---------------------------
            ! Full-chemistry simulation
            !---------------------------

            ! Number of tracers
            ! Add GLYX, MGLY, GLYC, HAC, C2H2,BENZ, TOLU, XYLE, C2H4
            !-----------------------------------------------------------------
            ! Bug fix: omit SOA tracers if not needed. (bmy, 9/2/10)
            !NTRAC(28) = 23
            !-----------------------------------------------------------------
            IF ( Input_Opt%LSOA ) THEN
               ! SOAupdate: plus one for NAP (hotp 7/25/10)
               NTRAC(28) = 24   ! Traditional SOA simulation
            ELSE
               NTRAC(28) = 14   ! Standard simulation
            ENDIF

            ! Loop over tracers
            DO T = 1, NTRAC(28)

               ! removed hardwires (fp, 6/2009)
               ! Case statement
               ! Replaced NOx emissions with NO emissions (mpayer, 3/14/13)
               IF ( T == IDBNO ) THEN
                  NAME (T,28) = 'NO'
                  INDEX(T,28) = IDTNO + ( SPACING * 45 )
                  MWT  (T,28) = 14e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'molec/cm2/s'
                  
               ELSE IF ( T == IDBCO ) THEN
                  NAME (T,28) = 'CO'
                  INDEX(T,28) = IDTCO + ( SPACING * 45 )
                  MWT  (T,28) = 28e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'molec/cm2/s'

               ELSE IF ( T == IDBALK4 ) THEN
                  NAME (T,28) = 'ALK4'
                  INDEX(T,28) = IDTALK4 + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 4
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBACET ) THEN
                  NAME (T,28) = 'ACET'
                  INDEX(T,28) = IDTACET + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 3
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBMEK ) THEN
                  NAME (T,28) = 'MEK'
                  INDEX(T,28) = IDTMEK + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 4
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBALD2 ) THEN
                  NAME (T,28) = 'ALD2'
                  INDEX(T,28) = IDTALD2 + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 3
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBPRPE ) THEN
                  NAME (T,28) = 'PRPE'
                  INDEX(T,28) = IDTPRPE + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 3
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBC3H8 ) THEN
                  NAME (T,28) = 'C3H8'
                  INDEX(T,28) = IDTC3H8 + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 3
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBCH2O ) THEN
                  NAME (T,28) = 'CH2O'
                  INDEX(T,28) = IDTCH2O + ( SPACING * 45 )
                  MWT  (T,28) = 30e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'molec/cm2/s'

               ELSE IF ( T == IDBC2H6 ) THEN
                  NAME (T,28) = 'C2H6'
                  INDEX(T,28) = IDTC2H6 + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 2
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBSO2 ) THEN
                  NAME (T,28) = 'SO2'
                  INDEX(T,28) = IDTSO2 + ( SPACING * 45 )
                  MWT  (T,28) = 32e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'atoms S/cm2/s'

               ELSE IF ( T == IDBNH3 ) THEN
                  NAME (T,28) = 'NH3'
                  INDEX(T,28) = IDTNH3 + ( SPACING * 45 )
                  MWT  (T,28) = 17e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'molec/cm2/s'

               ELSE IF ( T == IDBBC ) THEN
                  NAME (T,28) = 'BC'
                  INDEX(T,28) = IDTBCPI + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBOC ) THEN
                  NAME (T,28) = 'OC'
                  ! SOAupdate: allow non-vol or semi-vol POA (hotp 8/24/09)
                  IF ( IDTOCPI > 0 ) 
     &               INDEX(T,28) = IDTOCPI + ( SPACING * 45 )
                  IF ( IDTPOA1  > 0 )
     &               INDEX(T,28) = IDTPOA1 + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBGLYX ) THEN
                  NAME (T,28) = 'GLYX'
                  INDEX(T,28) = IDTGLYX + ( SPACING * 45 )
                  MWT  (T,28) = 58e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'molec/cm2/s'

               ELSE IF ( T == IDBMGLY ) THEN
                  NAME (T,28) = 'MGLY'
                  INDEX(T,28) = IDTMGLY + ( SPACING * 45 )
                  MWT  (T,28) = 72e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'molec/cm2/s'

               ELSE IF ( T == IDBBENZ ) THEN
                  NAME (T,28) = 'BENZ'
                  INDEX(T,28) = IDTBENZ + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 6
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBTOLU ) THEN
                  NAME (T,28) = 'TOLU'
                  INDEX(T,28) = IDTTOLU + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 7
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBXYLE ) THEN
                  NAME (T,28) = 'XYLE'
                  INDEX(T,28) = IDTXYLE + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 8
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBC2H4 ) THEN
                  NAME (T,28) = 'C2H4'
                  INDEX(T,28) = IDTC2H4 + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 2
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBC2H2 ) THEN
                  NAME (T,28) = 'C2H2'
                  INDEX(T,28) = IDTC2H2 + ( SPACING * 45 )
                  MWT  (T,28) = 12e-3
                  MOLC (T,28) = 2
                  UNIT (T,28) = 'atoms C/cm2/s'

               ELSE IF ( T == IDBGLYC ) THEN
                  NAME (T,28) = 'GLYC'
                  INDEX(T,28) = IDTGLYC + ( SPACING * 45 )
                  MWT  (T,28) = 60e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'molec/cm2/s'

               ELSE IF ( T == IDBHAC ) THEN
                  NAME (T,28) = 'HAC'
                  INDEX(T,28) = IDTHAC + ( SPACING * 45 )
                  MWT  (T,28) = 74e-3
                  MOLC (T,28) = 1
                  UNIT (T,28) = 'molec/cm2/s'

               ! SOAupdate: add NAPSOA (hotp 8/23/09)
               ELSEIF( T == IDBNAP ) THEN
                     NAME (T,28) = 'NAP'
                     INDEX(T,28) = IDTNAP + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'atoms C/cm2/s' 

               ENDIF

               ! Define other quantities
               FNAME(T,28) = TRIM( NAME(T,28) ) // ' biomass'
               SCALE(T,28) = 1e0
            ENDDO
         ENDIF
      ENDIF

      !-------------------------------------      
      ! CO emissions (ND29)
      !-------------------------------------
      IF ( ND29 > 0 ) THEN

         ! Number of tracers
         NTRAC(29) = 5

         ! Loop over tracers
         DO T = 1, NTRAC(29)

            ! Get name, longname
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,29) = 'COanth'
                  FNAME(T,29) = 'Anthropogenic CO'
               CASE( 2 )
                  NAME (T,29) = 'CObb'
                  FNAME(T,29) = 'Biomass CO'
               CASE( 3 )
                  NAME (T,29) = 'CObf'
                  FNAME(T,29) = 'Biofuel CO'
               CASE( 4 )
                  NAME (T,29) = 'COmeth'
                  FNAME(T,29) = 'CO from methanol'
               CASE( 5 )
                  NAME (T,29) = 'COmono'
                  FNAME(T,29) = 'CO from monoterpenes'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,29) = 'molec/cm2/s'
            INDEX(T,29) = T + ( SPACING * 20 )
            MOLC (T,29) = 1
            MWT  (T,29) = 28e-3
            SCALE(T,29) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Land map (ND30)
      !-------------------------------------      
      IF ( ND30 > 0 ) THEN
         T           = 1
         NTRAC(30)   = T
         NAME (T,30) = 'LWI'
         FNAME(T,30) = 'GMAO Land-water indices'
         UNIT (T,30) = 'unitless'
         INDEX(T,30) = T + ( SPACING * 15 )
         MOLC (T,30) = 1
         MWT  (T,30) = 0e0
         SCALE(T,30) = 1e0
      ENDIF

      !-------------------------------------      
      ! Surface pressure (ND31)
      !-------------------------------------   
      IF ( ND31 > 0       .or. DO_SAVE_DIAG48 .or. 
     &     DO_SAVE_DIAG49 .or. DO_SAVE_DIAG50 .or. 
     &     DO_SAVE_DIAG51 .or. DO_SAVE_DIAG51b     ) THEN
         T           = 1
         NTRAC(31)   = T
         NAME (T,31) = 'PSURF'
         FNAME(T,31) = 'Surface pressure'
         UNIT (T,31) = 'hPa'
         INDEX(T,31) = T + ( SPACING * 10 )
         MOLC (T,31) = 1
         MWT  (T,31) = 0e0
         SCALE(T,31) = 1e0
      ENDIF

      !-------------------------------------      
      ! CH3I emissions sources (ND36)
      !-------------------------------------      
      IF ( ND36 > 0 .and. Input_Opt%ITS_A_CH3I_SIM ) THEN

         ! Number of tracers
         NTRAC(36) = NEMANTHRO

         ! The first few are tracers
         DO T = 1, NTRAC(36)
            IF ( T <= Input_Opt%N_TRACERS ) THEN
               NAME (T,36) = Input_Opt%TRACER_NAME(T)
               FNAME(T,36) = TRIM( NAME(T,36) ) // ' tracer'
            ELSE IF ( T == 6 ) THEN
               NAME (T,36) = 'CH3Iof'
               FNAME(T,36) = 'CH3I ocean flux??'
            ELSE IF ( T == 7 ) THEN
               NAME (T,36) = 'CH3Igf'
               FNAME(T,36) = 'CH3I gas flux??'
            ELSE IF ( T == 8 ) THEN
               NAME (T,36) = 'CH3Itf'
               FNAME(T,36) = 'CH3I total flux??'
            ENDIF

            ! Define other quantities
            MOLC (T,36) = 1
            SCALE(T,36) = 1e0
            INDEX(T,36) = T + ( SPACING * 44 )
            MWT  (T,36) = 141.9e-3
            UNIT (T,36) = 'ng/m2/s'
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Afternoon-average boundary 
      ! layer heights (ND41) + timeseries
      !-------------------------------------      
      IF ( ND41 > 0 .or. DO_TIMESERIES ) THEN

         ! Number of tracers
         NTRAC(41) = 2

         ! Loop over tracers
         DO T = 1, NTRAC(41)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME(T,41) = 'PBL-M'
                  UNIT(T,41) = 'm'
               CASE( 2 )
                  NAME(T,41) = 'PBL-L'
                  UNIT(T,41) = 'levels'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            FNAME(T,41) = 'PBL depth'
            INDEX(T,41) = T + ( SPACING * 27 )
            MOLC (T,41) = 1
            MWT  (T,41) = 0e0
            SCALE(T,41) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! SOA concentrations (ND42)
      !-------------------------------------    
      IF ( ND42 > 0 ) THEN

         ! Number of tracers
         ! Add new diagnostics with SOA5. (ccc, 3/10/10)
         !NTRAC(42) = 8
         NTRAC(42) = 24

         ! Loop over tracers
         DO T = 1, NTRAC(42)

            ! Get name, long name, unit, and mol wt for each field
            SELECT CASE( T )
               ! SOAupdate: updated (hotp 7/28/10)
               CASE( 1 )
                  NAME (T,42) = 'TSOA'
                  FNAME(T,42) = 'Aer prods of terpene ox'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 150e-3
               CASE( 2 )
                  NAME (T,42) = 'ISOA'
                  FNAME(T,42) = 'Aer prods of isoprene ox'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 150e-3
               CASE( 3 )
                  NAME (T,42) = 'ASOA'
                  FNAME(T,42) = 'Aer prods of arom+IVOC ox'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 150e-3
               CASE( 4 )
                  NAME (T,42) = 'POA'
                  FNAME(T,42) = 'Aer from SVOCs'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 12e-3
               CASE( 5 )
                  NAME (T,42) = 'OPOA'
                  FNAME(T,42) = 'Aer prods of pog ox'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 12e-3
               CASE( 6 )
                  NAME (T,42) = 'sumOA'
                  FNAME(T,42) = 'all organic aerosol'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 150e-3
               CASE( 7 )
                  NAME (T,42) = 'sumOC'
                  FNAME(T,42) = 'all organic carbon'
                  UNIT (T,42) = 'ugC/m3'
                  MWT  (T,42) = 12e-3
               CASE( 8 )
                  NAME (T,42) = 'biogOA'
                  FNAME(T,42) = 'all biogenic OA'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 1e-6
               CASE( 9 )
                  NAME (T,42) = 'BNO'
                  FNAME(T,42) = 'NO branching ratio: betaNO'
                  UNIT (T,42) = 'dimless'
                  MWT  (T,42) = 1e-6
               CASE( 10 )
                  NAME (T,42) = 'POA'
                  FNAME(T,42) = 'Aer from SVOCs'
                  UNIT (T,42) = 'ugC/m3'
                  MWT  (T,42) = 12e-3
               CASE( 11 )
                  NAME (T,42) = 'OPOA'
                  FNAME(T,42) = 'Aer prods of pog ox'
                  UNIT (T,42) = 'ugC/m3'
                  MWT  (T,42) = 12e-3
               ! end hotp
               CASE( 12 )
                  NAME (T,42) = 'OC'
                  FNAME(T,42) = 'OCPI + OCPO'
                  UNIT (T,42) = 'ug C/m3'
                  MWT  (T,42) = 12e-3
               CASE( 13 )
                  NAME (T,42) = 'BC'
                  FNAME(T,42) = 'BCPI + BCPO'
                  UNIT (T,42) = 'ug C/m3'
                  MWT  (T,42) = 12e-3
               CASE( 14 )
                  NAME (T,42) = 'SO4'
                  FNAME(T,42) = 'sulfate'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 96e-3
               CASE( 15 )
                  NAME (T,42) = 'NH4'
                  FNAME(T,42) = 'ammonium'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 18e-3
               CASE( 16 )
                  NAME (T,42) = 'NIT'
                  FNAME(T,42) = 'nitrate'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 62e-3
               CASE( 17 )
                  NAME (T,42) = 'SAL'
                  FNAME(T,42) = 'SALA + SALC'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 36e-3
               CASE( 19 )
                  NAME (T,42) = 'SOAG'
                  FNAME(T,42) = 'SOAG'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 1e-6
               CASE( 20 )
                  NAME (T,42) = 'SOAM'
                  FNAME(T,42) = 'SOAM'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 1e-6
                  
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            INDEX(T,42) = T + ( SPACING * 43 )
            MOLC (T,42) = 1
            SCALE(T,42) = 1e0
         ENDDO
      ENDIF
      
      !-------------------------------------      
      ! Chemically-produced OH, etc (ND43)
      !-------------------------------------      
      IF ( ND43 > 0 ) THEN 

         ! Number of tracers
         ! update for arom (dkh, 06/22/07)  
         !NTRAC(43) = 5
         NTRAC(43) = 11

         ! Loop over tracers
         DO T = 1, NTRAC(43)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME(T,43) = 'OH'
                  UNIT(T,43) = 'molec/cm3'
                  MWT (T,43) = 17e-3
               CASE( 3 )
                  NAME(T,43) = 'HO2'
                  UNIT(T,43) = 'v/v'
                  MWT (T,43) = 33e-3
               ! update for arom (dkh, 06/22/07)  
               CASE( 6 )
                  NAME(T,43) = 'LBRO2H'
                  UNIT(T,43) = 'molec/cm3'
                  MWT (T,43) = 127e-3
               CASE( 7 )
                  NAME(T,43) = 'LBRO2N'
                  UNIT(T,43) = 'molec/cm3'
                  MWT (T,43) = 127e-3
               CASE( 8 )
                  NAME(T,43) = 'LTRO2H'
                  UNIT(T,43) = 'molec/cm3'
                  MWT (T,43) = 141e-3
               CASE( 9 )
                  NAME(T,43) = 'LTRO2N'
                  UNIT(T,43) = 'molec/cm3'
                  MWT (T,43) = 141e-3
               CASE( 10)
                  NAME(T,43) = 'LXRO2H'
                  UNIT(T,43) = 'molec/cm3'
                  MWT (T,43) = 155e-3
               CASE( 11)
                  NAME(T,43) = 'LXRO2N'
                  UNIT(T,43) = 'molec/cm3'
                  MWT (T,43) = 155e-3
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            FNAME(T,43) = 'Chemically produced ' // TRIM( NAME(T,43) )
            INDEX(T,43) = T + ( SPACING * 16 )
            MOLC (T,43) = 1
            SCALE(T,43) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Dry deposition fluxes     (ND44)
      ! Dry deposition velocities (ND44)
      !-------------------------------------      
      IF ( ND44 > 0 ) THEN

         ! Special handling for tagged simulations
         IF ( Input_Opt%ITS_A_TAGOX_SIM   .or. 
     &        Input_Opt%ITS_A_MERCURY_SIM .or. 
     &        Input_Opt%ITS_A_H2HD_SIM  ) THEN

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%             FOR TAGGED SIMULATIONS ONLY             %%%
            !%%%                                                     %%%
            !%%%          Save drydep flux for all tracers           %%%
            !%%%     but drydep velocity for only NUMDEP species     %%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! Drydep flux (all tracers)
            DO T = 1, Input_Opt%N_TRACERS
               NAME (T,44)  = TRIM( Input_Opt%TRACER_NAME(T) ) // 'df'        
               FNAME(T,44)  = TRIM( Input_Opt%TRACER_NAME(T) ) // 
     &                        ' drydep flux'
               UNIT (T,44)  = 'molec/cm2/s'
               MWT  (T,44)  = Input_Opt%TRACER_MW_KG(T) 
               MOLC (T,44)  = INT( Input_Opt%TRACER_COEFF(T,1) )
               SCALE(T,44)  = 1.0e0
               INDEX(T,44)  = T + ( SPACING * 36 )
               NTRAC(44)    = NTRAC(44) + 1
            ENDDO

            ! Drydep velocity (only deposited species)
            DO N = 1, NUMDEP
               T            = NTRAIND(N)
               NN           = N + Input_Opt%N_TRACERS
               NAME (NN,44) = TRIM( DEPNAME(N) ) // 'dv'        
               FNAME(NN,44) = TRIM( DEPNAME(N) ) // 
     &                        ' drydep velocity'
               UNIT (NN,44) = 'cm/s'
               MWT  (NN,44) = Input_Opt%TRACER_MW_KG(T) 
               MOLC (NN,44) = INT( Input_Opt%TRACER_COEFF(T,1) )
               SCALE(NN,44) = 1.0e0
               INDEX(NN,44) = T + ( SPACING * 37 )
               NTRAC(44)    = NTRAC(44) + 1
            ENDDO

         ELSE

#if   defined( TOMAS )
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%        FOR TOMAS MICROPHYSICS SIMULATIONS           %%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! Loop over # of dry depositing species
            DO N = 1, NUMDEP + (ICOMP-IDIAG)* IBINS
               
               IF ( N <= NUMDEP ) THEN 

                  ! Handle drydep flux first
                  SELECT CASE( TRIM( Input_Opt%DEPNAME(N) ) )

                     CASE( 'ISOPNB' )

                        !-----------------------------------------------
                        ! Drydep flux: ISOPNB.  Use the same info as
                        ! for ISOPND but save as w/ index N_TRACERS+1.
                        !-----------------------------------------------
                        T             = IDTISOPN
                        NAME (N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                                  'df'
                        FNAME(N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                                  ' drydep flux'
                        MWT  (N,44)   = Input_Opt%TRACER_MW_KG(T) 
                        MOLC (N,44)   = INT(Input_Opt%TRACER_COEFF(T,1))
                        SCALE(N,44)   = 1.0e0
                        UNIT (N,44)   = 'molec/cm2/s'
                        INDEX(N,44)   = ( Input_Opt%N_TRACERS + 1  )
     &                                + ( SPACING             * 36 )
                        NTRAC(44)     = NTRAC(44) + 1

                        !-----------------------------------------------
                        ! Drydep velocity: ISOPNB.  Use the same info as
                        ! for ISOPND but save as w/ index N_TRACERS+1.
                        !-----------------------------------------------
                        T             = IDTISOPN
                        NN            = N + Input_Opt%NUMDEP
                        NAME (NN,44)  = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                                 'dv'
                        FNAME(NN,44)  = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                                  ' drydep velocity'
                        MWT  (NN,44)  = Input_Opt%TRACER_MW_KG(T) 
                        MOLC (NN,44)  = INT(Input_Opt%TRACER_COEFF(T,1))
                        SCALE(NN,44)  = 1.0e0
                        UNIT (NN,44)  = 'cm/s'
                        INDEX(NN,44)  = ( Input_Opt%N_TRACERS + 1  )
     &                                + ( SPACING             * 37 )
                        NTRAC(44)     = NTRAC(44) + 1

                     CASE( 'MVKN' )

                        !-----------------------------------------------
                        ! Drydep flux: MVKN.  Use the same info as
                        ! for MMN but save as w/ index N_TRACERS+2.
                        !-----------------------------------------------
                        T             = IDTMMN
                        NAME (N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                       'df'
                        FNAME(N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                       ' drydep flux'
                        MWT  (N,44)   = Input_Opt%TRACER_MW_KG(T) 
                        MOLC (N,44)   = INT(Input_Opt%TRACER_COEFF(T,1))
                        SCALE(N,44)   = 1.0e0
                        UNIT (N,44)   = 'molec/cm2/s'
                        INDEX(N,44)   = ( Input_Opt%N_TRACERS + 2  )
     &                                + ( SPACING             * 36 )
                        NTRAC(44)     = NTRAC(44) + 1
                        
                        !--------------------------------------------------
                        ! Drydep velocity: MVKN.  Use the same info as
                        ! for MMN but save as w/ index N_TRACERS+2.
                        !--------------------------------------------------
                        T             = IDTMMN
                        NN            = N + Input_Opt%NUMDEP
                        NAME (NN,44)  = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                                  'dv'
                        FNAME(NN,44)  = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                                  ' drydep flux'
                        MWT  (NN,44)  = Input_Opt%TRACER_MW_KG(T) 
                        MOLC (NN,44)  = INT(Input_Opt%TRACER_COEFF(T,1))
                        SCALE(NN,44)  = 1.0e0
                        UNIT (NN,44)  = 'cm/s'
                        INDEX(NN,44)  = ( Input_Opt%N_TRACERS + 2  )
     &                                + ( SPACING             * 37 )
                        NTRAC(44)     = NTRAC(44) + 1

                     CASE DEFAULT

                        !-----------------------------------------------
                        ! Drydep flux: all other G-C tracers
                        !-----------------------------------------------
                        T             = NTRAIND(N)         
                        NAME (N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                                  'df'
                        FNAME(N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                                 ' drydep flux'
                        MWT  (N,44)   = Input_Opt%TRACER_MW_KG(T) 
                        MOLC (N,44)   = INT(Input_Opt%TRACER_COEFF(T,1))
                        SCALE(N,44)   = 1.0e0
                        INDEX(N,44)   = T + ( SPACING * 36 )
                        NTRAC(44)     = NTRAC(44) + 1

                        ! For the Rn simulation, unit is kg/s (bmy, 2/22/08)
                        IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN
                           UNIT(N,44) = 'kg/s'
                        ELSE
                           UNIT(N,44) = 'molec/cm2/s'
                        ENDIF

                        !-----------------------------------------------
                        ! Drydep velocity all other G-C tracers
                        !-----------------------------------------------
                        NN            = N + NUMDEP  +(ICOMP-IDIAG)*IBINS
                        NAME (NN,44)  = TRIM( DEPNAME(N) ) // 'dv'        
                        FNAME(NN,44)  = TRIM( DEPNAME(N) ) // 
     &                                  ' drydep velocity'
                        MWT  (NN,44)  = Input_Opt%TRACER_MW_KG(T) 
                        MOLC (NN,44)  = INT(Input_Opt%TRACER_COEFF(T,1))
                        UNIT (NN,44)  = 'cm/s'                     
                        SCALE(NN,44)  = 1.0e0
                        INDEX(NN,44)  = T + ( SPACING * 37 )
                        NTRAC(44)     = NTRAC(44) + 1

                  END SELECT

               ELSE

                  !-----------------------------------------------
                  ! Drydep flux: TOMAS aerosol tracers
                  !-----------------------------------------------
                  T = IDTNK1 + IBINS - 1 + ( N - NUMDEP )

                  ! Drydep flux (extended deposited species)
                  NAME (N,44)  = TRIM( Input_Opt%TRACER_NAME(T) ) // 
     &                          'df'
                  FNAME(N,44)  = TRIM( Input_Opt%TRACER_NAME(T) ) // 
     &                          ' drydep flux'
                  UNIT (N,44)  = 'molec/cm2/s'
                  MWT  (N,44)  = Input_Opt%TRACER_MW_KG(T) 
                  MOLC (N,44)  = INT( Input_Opt%TRACER_COEFF(T,1) )
                  SCALE(N,44)  = 1.0e0
                  INDEX(N,44)  = T + ( SPACING * 36 )
                  NTRAC(44)    = NTRAC(44) + 1
                  
                  !-----------------------------------------------
                  ! Drydep velocity all other G-C tracers
                  !-----------------------------------------------
                  NN            = N + NUMDEP + (ICOMP-IDIAG)* IBINS
                  NAME (NN,44)  = TRIM( Input_Opt%TRACER_NAME(T) ) // 
     &                            'dv'        
                  FNAME(NN,44)  = TRIM( Input_Opt%TRACER_NAME(T) ) //
     &                           ' drydep velocity'
                  MWT  (NN,44)  = Input_Opt%TRACER_MW_KG(T) 
                  MOLC (NN,44)  = INT(Input_Opt%TRACER_COEFF(T,1))
                  UNIT (NN,44)  = 'cm/s'                     
                  SCALE(NN,44)  = 1.0e0
                  INDEX(NN,44)  = T + ( SPACING * 37 )
                  NTRAC(44)     = NTRAC(44) + 1

               ENDIF
            ENDDO
#else
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%       FOR REGULAR FULL-CHEMISTRY SIMULATIONS        %%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! Loop over # of depositing species
            DO N = 1, Input_Opt%NUMDEP

               ! Handle drydep flux first
               SELECT CASE( TRIM( Input_Opt%DEPNAME(N) ) )

                  CASE( 'ISOPNB' )

                     !--------------------------------------------------
                     ! Drydep flux: ISOPNB.  Use the same info as
                     ! for ISOPND but save as w/ index N_TRACERS+1.
                     !--------------------------------------------------
                     T             = IDTISOPN
                     NAME (N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                               'df'
                     FNAME(N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                               ' drydep flux'
                     MWT  (N,44)   = Input_Opt%TRACER_MW_KG(T) 
                     MOLC (N,44)   = INT( Input_Opt%TRACER_COEFF(T,1) )
                     SCALE(N,44)   = 1.0e0
                     UNIT (N,44)   = 'molec/cm2/s'
                     INDEX(N,44)   = ( Input_Opt%N_TRACERS + 1  )
     &                             + ( SPACING             * 36 )
                     NTRAC(44)     = NTRAC(44) + 1

                     !--------------------------------------------------
                     ! Drydep velocity: ISOPNB.  Use the same info as
                     ! for ISOPND but save as w/ index N_TRACERS+1.
                     !--------------------------------------------------
                     T             = IDTISOPN
                     NN            = N + Input_Opt%NUMDEP
                     NAME (NN,44)  = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                              'dv'
                     FNAME(NN,44)  = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                             ' drydep velocity'
                     MWT  (NN,44)  = Input_Opt%TRACER_MW_KG(T) 
                     MOLC (NN,44)  = INT( Input_Opt%TRACER_COEFF(T,1) )
                     SCALE(NN,44)  = 1.0e0
                     UNIT (NN,44)  = 'cm/s'
                     INDEX(NN,44)  = ( Input_Opt%N_TRACERS + 1  )
     &                             + ( SPACING             * 37 )
                     NTRAC(44)     = NTRAC(44) + 1

                  CASE( 'MVKN' )

                     !--------------------------------------------------
                     ! Drydep flux: MVKN.  Use the same info as
                     ! for MMN but save as w/ index N_TRACERS+2.
                     !--------------------------------------------------
                     T             = IDTMMN
                     NAME (N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                               'df'
                     FNAME(N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                              ' drydep flux'
                     MWT  (N,44)   = Input_Opt%TRACER_MW_KG(T) 
                     MOLC (N,44)   = INT( Input_Opt%TRACER_COEFF(T,1) )
                     SCALE(N,44)   = 1.0e0
                     UNIT (N,44)   = 'molec/cm2/s'
                     INDEX(N,44)   = ( Input_Opt%N_TRACERS + 2  )
     &                             + ( SPACING             * 36 )
                     NTRAC(44)     = NTRAC(44) + 1
                    
                     !--------------------------------------------------
                     ! Drydep velocity: MVKN.  Use the same info as
                     ! for MMN but save as w/ index N_TRACERS+2.
                     !--------------------------------------------------
                     T             = IDTMMN
                     NN            = N + Input_Opt%NUMDEP
                     NAME (NN,44)  = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                              'dv'
                     FNAME(NN,44)  = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                             ' drydep velocity'
                     MWT  (NN,44)  = Input_Opt%TRACER_MW_KG(T) 
                     MOLC (NN,44)  = INT( Input_Opt%TRACER_COEFF(T,1) )
                     SCALE(NN,44)  = 1.0e0
                     UNIT (NN,44)  = 'cm/s'
                     INDEX(NN,44)  = ( Input_Opt%N_TRACERS + 2  )
     &                             + ( SPACING             * 37 )
                     NTRAC(44)     = NTRAC(44) + 1

                  CASE DEFAULT

                     !--------------------------------------------------
                     ! Drydep flux: all other G-C tracers
                     !--------------------------------------------------
                     T             = NTRAIND(N)         
                     NAME (N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                               'df'
                     FNAME(N,44)   = TRIM( Input_Opt%DEPNAME(N) ) // 
     &                              ' drydep flux'
                     MWT  (N,44)   = Input_Opt%TRACER_MW_KG(T) 
                     MOLC (N,44)   = INT( Input_Opt%TRACER_COEFF(T,1) )
                     SCALE(N,44)   = 1.0e0
                     INDEX(N,44)   = T + ( SPACING * 36 )
                     NTRAC(44)     = NTRAC(44) + 1

                     ! For the Rn simulation, unit is kg/s (bmy, 2/22/08)
                     IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN
                        UNIT(N,44) = 'kg/s'
                     ELSE
                        UNIT(N,44) = 'molec/cm2/s'
                     ENDIF

                     !--------------------------------------------------
                     ! Drydep velocity all other G-C tracers
                     !--------------------------------------------------
                     NN            = N + Input_Opt%NUMDEP
                     NAME (NN,44)  = TRIM( DEPNAME(N) ) // 'dv'        
                     FNAME(NN,44)  = TRIM( DEPNAME(N) ) // 
     &                               ' drydep velocity'
                     MWT  (NN,44)  = Input_Opt%TRACER_MW_KG(T) 
                     MOLC (NN,44)  = INT( Input_Opt%TRACER_COEFF(T,1) )
                     UNIT (NN,44)  = 'cm/s'                     
                     SCALE(NN,44)  = 1.0e0
                     INDEX(NN,44)  = T + ( SPACING * 37 )
                     NTRAC(44)     = NTRAC(44) + 1

               END SELECT

            ENDDO
#endif
         ENDIF
      ENDIF

      !-------------------------------------      
      ! Biogenic emissions (ND46)
      !-------------------------------------      
      IF ( ND46 > 0 ) THEN 

         ! Number of tracers
         NTRAC(46) = 20

         ! Loop over tracers
         DO T = 1, NTRAC(46)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME(T,46) = 'ISOP'
                  MOLC(T,46) = 5
               CASE( 2 )
                  NAME(T,46) = 'ACET'
                  MOLC(T,46) = 3
               CASE( 3 )
                  NAME(T,46) = 'PRPE'
                  MOLC(T,46) = 3
               CASE( 4 )
                  NAME(T,46) = 'MONOT'
                  MOLC(T,46) = 10
               CASE( 5 )
                  NAME(T,46) = 'MBO'
                  MOLC(T,46) = 5
               CASE( 6 )
                  NAME(T,46) = 'C2H4'
                  MOLC(T,46) = 2
               CASE( 7 )
                  NAME(T,46) = 'APINE'
                  MOLC(T,46) = 10
               CASE( 8 )
                  NAME(T,46) = 'BPINE'
                  MOLC(T,46) = 10
               CASE( 9 )
                  NAME(T,46) = 'LIMON'
                  MOLC(T,46) = 10
               CASE( 10)
                  NAME(T,46) = 'SABIN'
                  MOLC(T,46) = 10
               CASE( 11)
                  NAME(T,46) = 'MYRCN'
                  MOLC(T,46) = 10
               CASE( 12)
                  NAME(T,46) = 'CAREN'
                  MOLC(T,46) = 10
               CASE( 13)
                  NAME(T,46) = 'OCIMN'
                  MOLC(T,46) = 10
               ! SOAupdate: Biogenic emissions from sesq hotp 3/1/10
               ! MOLC is #C/molecule (10 for mono, 15 for sesq) hotp 5/10/10
               CASE( 14)
                  NAME(T,46) = 'FARN'
                  MOLC(T,46) = 15 ! [C/molecule]
               CASE( 15)              
                  NAME(T,46) = 'BCAR'
                  MOLC(T,46) = 15 ! [C/molecule]
               CASE( 16)              
                  NAME(T,46) = 'OSQT'
                  MOLC(T,46) = 15 ! [C/molecule]
               CASE( 17)              
                  NAME(T,46) = 'OMTP'
                  MOLC(T,46) = 10 ! [C/molecule]
               !jpp, 9/6/07
               CASE( 18 )
                  INDEX(T,46)  = T + ( SPACING * 21 )
                  NAME  (T,46) = 'CHBr3'
                  FNAME (T,46) = 'Ocean bromoform emissions'
                  UNIT  (T,46) = 'kg/m2/s'
                  MWT   (T,46) = 2.53e-1   !in [kg/mol]
                  SCALE(T,46) = 1e0
               CASE( 19 )
                  INDEX (T,46)  = T + ( SPACING * 21 )
                  NAME  (T,46) = 'CH2Br2'
                  FNAME (T,46) = 'Ocean dibromomethane emissions'
                  UNIT  (T,46) = 'kg/m2/s'
                  MWT   (T,46) = 1.74e-1   !in [kg/mol]
                  SCALE (T,46) = 1e0
               CASE( 20 )
                  INDEX (T,46)  = T + ( SPACING * 21 )
                  NAME  (T,46) = 'SSBr2'
                  FNAME (T,46) = 'sea-salt Br2 emissions'
                  UNIT  (T,46) = 'kg/m2/s'
                  MWT   (T,46) = 1.60e-1   !in [kg/mol]
                  SCALE (T,46) = 1e0
               CASE DEFAULT
                  ! Nothing
            END SELECT

            !jpp, 9/6/07: preventing these default quantities
            !             from overwriting CHBr3's info.
            IF ( T < 18 ) THEN
               ! Define the rest of the quantities
               INDEX(T,46) = T + ( SPACING * 21 )
               FNAME(T,46) = TRIM( NAME(T,46) ) // ' emissions'
               UNIT (T,46) = 'atoms C/cm2/s'
               MWT  (T,46) = 12e-3
               SCALE(T,46) = 1e0
            ENDIF
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Time series diags (ND{48,49,50,51})
      !-------------------------------------      
      IF ( DO_TIMESERIES ) THEN 

         ! Number of tracers
         ! Increased from 26 to 32 (mpb,2009)
         NTRAC(48) = 32

         ! Loop over tracers
         DO T = 1, NTRAC(48)

            ! Define quantities
            INDEX(T,48) = T + ( SPACING * 19 )
            MOLC (T,48) = 1
            MWT  (T,48) = 0e0
            SCALE(T,48) = 1e0

            ! Get name, long-name (and others where necessary)
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,48) = 'O3'
                  FNAME(T,48) = 'O3 time series'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 48e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 2  )
                  NAME (T,48) = 'OH'
                  FNAME(T,48) = 'OH time series'
                  UNIT (T,48) = 'molec/cm3'
               CASE( 3  )
                  NAME (T,48) = 'NOy'
                  FNAME(T,48) = 'NOy time series'
                  UNIT (T,48) = 'ppbv'
                  SCALE(T,48) = 1.0e+9
               CASE( 4  )
                  NAME (T,48) = 'NO2dv'
                  FNAME(T,48) = 'NO2 drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 5  )
                  NAME (T,48) = 'O3dv'
                  FNAME(T,48) = 'O3 drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 6  )
                  NAME (T,48) = 'PANdv'
                  FNAME(T,48) = 'PAN drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 7  )
                  NAME (T,48) = 'HNO3dv'
                  FNAME(T,48) = 'HNO3 drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 8  )
                  NAME (T,48) = 'H2O2dv'
                  FNAME(T,48) = 'H2O2 drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 9  )
                  NAME (T,48) = 'NO'
                  FNAME(T,48) = 'NO time series'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 30e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 10  )
                  NAME (T,48) = 'PBLTOP'
                  FNAME(T,48) = 'PBL top pressure'
                  UNIT (T,48) = 'hPa'
               CASE( 11  )
                  NAME (T,48) = 'LT'
                  FNAME(T,48) = 'Station local time'
                  UNIT (T,48) = 'h'
               CASE( 12 )
                  NAME (T,48) = 'ISOPf'
                  FNAME(T,48) = 'ISOP emission flux'
                  UNIT (T,48) = 'molec/cm2/s'
                  MWT  (T,48) = 12e-3
               CASE( 13 )
                  NAME (T,48) = 'MOLEN'
                  FNAME(T,48) = 'Monin-Obhukov length'
                  UNIT (T,48) = 'm'
               CASE( 14 )
                  NAME (T,48) = 'NO3'
                  FNAME(T,48) = 'NO3 time series'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 62e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 15 )
                  NAME (T,48) = 'PS'
                  FNAME(T,48) = 'Surface pressure'
                  UNIT (T,48) = 'hPa'
               CASE( 16 )
                  NAME (T,48) = 'AVGW'
                  FNAME(T,48) = 'Water vapor mixing ratio'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 33e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 17 )
                  NAME (T,48) = 'RH'
                  FNAME(T,48) = 'Relative humidity'
                  UNIT (T,48) = '%'
               CASE( 18 )
                  NAME (T,48) = 'HCNcol'
                  FNAME(T,48) = 'HCN column density'
                  UNIT (T,48) = 'molec/cm2'
               CASE( 19 )
                  NAME (T,48) = 'CF'
                  FNAME(T,48) = 'Cloud fraction'
                  UNIT (T,48) = 'unitless'
               CASE( 20 )
                  NAME (T,48) = 'ODCOL'
                  FNAME(T,48) = 'Cloud optical depth'
                  UNIT (T,48) = 'unitless'
               CASE( 21 )
                  NAME (T,48) = 'CThgt'
                  FNAME(T,48) = 'Cloud top height'
                  UNIT (T,48) = 'hPa'
               CASE( 22 )
                  NAME (T,48) = 'AIRDEN'
                  FNAME(T,48) = 'Air density'
                  UNIT (T,48) = 'molec/cm3'
               CASE( 23 )
                  NAME (T,48) = 'TOTDUST'
                  FNAME(T,48) = 'Total dust opt depth'
                  UNIT (T,48) = 'unitless'
               CASE( 24 )
                  NAME (T,48) = 'TOTSALT'
                  FNAME(T,48) = 'Total seasalt tracer'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 36e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 25 )
                  NAME (T,48) = 'NO2'
                  FNAME(T,48) = 'NO2 concentration'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 46e-3
                  SCALE(T,48) = 1.0e+9                  
               CASE( 26 )
                  NAME (T,48) = 'TOTAOD'
                  FNAME(T,48) = 'Total aerosol opt depth'
                  UNIT (T,48) = 'unitless'
                  ! +++++++++++++++++++++++++++++++++++++
                  ! Add in gamma factors (mpb,2008)
                  ! +++++++++++++++++++++++++++++++++++++
               CASE( 27 )
                  NAME (T,48) = 'G_T'
                  FNAME(T,48) = 'Temperature activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 28 )
                  NAME (T,48) = 'G_PAR'
                  FNAME(T,48) = 'Light activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 29 )              
                  NAME (T,48) = 'G_AGE'
                  FNAME(T,48) = 'Leaf-age activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 30 )
                  NAME (T,48) = 'G_LAI'
                  FNAME(T,48) = 'LAI activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 31 )
                  NAME (T,48) = 'G_SM'
                  FNAME(T,48) = 'Soil-moisture activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 32 )
                  NAME (T,48) = 'D_LAI'
                  FNAME(T,48) = 'Daily LAI'
                  UNIT (T,48) = 'm2/m2'
               CASE DEFAULT
                  ! Nothing
            END SELECT

         ENDDO
      ENDIF


      !-----------------------------------
      ! gamma HO2 (ND52) jaegle (02/26/09)
      !-----------------------------------
      IF ( ND52 > 0 ) THEN
         T           = 1
         NTRAC(52)   = T
         NAME (T,52) = 'GAMMAHO2'
         FNAME(T,52) = 'Gamma HO2'
         UNIT (T,52) = 'unitless'
         INDEX(T,52) = T + ( SPACING * 49 )
         MOLC (T,52) = 1
         MWT  (T,52) = 0e0
         SCALE(T,52) = 1e0
      ENDIF


      !-----------------------------------
      ! POPS (ND53) clf (03/22/11)
      !-----------------------------------
      IF ( ND53 > 0 ) THEN
         NTRAC(53)  = PD53

         ! Loop over tracers
         DO T = 1, PD53

            ! Define quantities
            UNIT (T,53) = 'kg'
            MOLC (T,53) = 1
            MWT  (T,53) = GET_POP_XMW(DUM)
            SCALE(T,53) = 1e0

            ! Get name, long-name, index, and new units
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,53) = 'POPT'
                  FNAME(T,53) = 'POP total emissions'
                  INDEX(T,53) = T + ( SPACING * 52 )
               CASE( 2  )
                  NAME (T,53) = 'POPPOC'
                  FNAME(T,53) = 'POPPOC emissions'
                  INDEX(T,53) = T + ( SPACING * 52 )
               CASE( 3  )
                  NAME (T,53) = 'POPPBC'
                  FNAME(T,53) = 'POPPBC emissions'
                  INDEX(T,53) = T + ( SPACING * 52 )
               CASE( 4  )
                  NAME (T,53) = 'POPG'
                  FNAME(T,53) = 'POPG emissions'
                  INDEX(T,53) = T + ( SPACING * 52 )
               CASE( 5  )
                  NAME (T,53) = 'POPPOC'
                  FNAME(T,53) = 'Gross POPPOC lost to gas'
                  INDEX(T,53) = T + ( SPACING * 53 )
               CASE( 6  )
                  NAME (T,53) = 'POPPOC'
                  FNAME(T,53) = 'Gross POPPOC formed from gas'
                  INDEX(T,53) = T + ( SPACING * 53 )
               CASE( 7  )
                  NAME (T,53) = 'POPPBC'
                  FNAME(T,53) = 'Gross POPPBC lost to gas'
                  INDEX(T,53) = T + ( SPACING * 53 )
               CASE( 8  )
                  NAME (T,53) = 'POPPBC'
                  FNAME(T,53) = 'Gross POPPBC formed from gas'
                  INDEX(T,53) = T + ( SPACING * 53 )
               CASE( 9  )
                  NAME (T,53) = 'POPG'
                  FNAME(T,53) = 'POPG oxidized by OH'
                  INDEX(T,53) = T + ( SPACING * 53 )
               CASE( 10  )
                  NAME (T,53) = 'POPPOC'
                  FNAME(T,53) = 'POPPOC oxidized by O3'
                  INDEX(T,53) = T + ( SPACING * 53 )
               CASE( 11  )
                  NAME (T,53) = 'POPPBC'
                  FNAME(T,53) = 'POPPBC oxidized by O3'
                  INDEX(T,53) = T + ( SPACING * 53 )
            END SELECT
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Time in the troposphere (ND54)
      !-------------------------------------      
      IF ( ND54 > 0 ) THEN
         T           = 1
         NTRAC(54)   = T
         NAME (T,54) = 'TIMETROP'
         FNAME(T,54) = 'Time in the troposphere'
         UNIT (T,54) = 'unitless'
         INDEX(T,54) = T + ( SPACING * 46 )
         MOLC (T,54) = 1
         MWT  (T,54) = 0e0
         SCALE(T,54) = 1e0
      ENDIF

      !-------------------------------------      
      ! Tropopause quantities (ND55)
      !-------------------------------------      
      IF ( ND55 > 0 ) THEN 

         ! Number of tracers
         NTRAC(55) = 3

         ! Loop over tracers
         DO T = 1, NTRAC(55)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,55) = 'TP-LEVEL'
                  FNAME(T,55) = 'Tropopause level'
                  UNIT (T,55) = 'level'
               CASE( 2 )
                  NAME (T,55) = 'TP-HGHT'
                  FNAME(T,55) = 'Tropopause height'
                  UNIT (T,55) = 'km'
               CASE( 3 )
                  NAME (T,55) = 'TP-PRESS'
                  FNAME(T,55) = 'Tropopause pressure'
                  UNIT (T,55) = 'hPa'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            INDEX(T,55) = T + ( SPACING * 26 )
            MOLC (T,55) = 1
            MWT  (T,55) = 0e0
            SCALE(T,55) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Lightning flash rates (ND56)
      !-------------------------------------      
      IF ( ND56 > 0 ) THEN 

         ! Number of tracers
         NTRAC(56) = 3

         ! Loop over tracers
         DO T = 1, NTRAC(56)

            ! Get name and long name for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,56) = 'L-RATE'
                  FNAME(T,56) = 'Total lightning flash rate'
               CASE( 2 )
                  NAME (T,56) = 'IC-RATE'
                  FNAME(T,56) = 'Intra-cloud flash rate'
               CASE( 3 )
                  NAME (T,56) = 'CG-RATE'
                  FNAME(T,56) = 'Cloud-ground flash rate'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,56) = 'flashes/min/km2'
            INDEX(T,56) = T + ( SPACING * 42 )
            MOLC (T,56) = 1
            MWT  (T,56) = 0e0
            SCALE(T,56) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! THETA (ND57) (FP 6/2009)
      !-------------------------------------      
      IF ( ND57 > 0 ) THEN
         T           = 1
         NTRAC(57)   = T
         NAME (T,57) = 'THETA'
         FNAME(T,57) = 'Potential Temperature'
         UNIT (T,57) = 'K'
         INDEX(T,57) = T + ( SPACING * 57 )
         MOLC (T,57) = 1
         MWT  (T,57) = 0e0
         SCALE(T,57) = 1e0
      ENDIF

      !-------------------------------------
      ! CH4 Emissions
      !-------------------------------------
      IF ( ND58 > 0 ) THEN

         ! Number of tracers
         NTRAC(58) = 12

         ! Loop over tracers
         DO T = 1, NTRAC(58)

            ! Get name and long name for each tracer
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,58) = 'CH4-TOT'
                  FNAME(T,58) = 'CH4 Emissions total (w/o sab)'
               CASE( 2 )
                  NAME (T,58) = 'CH4-GAO'
                  FNAME(T,58) = 'CH4 Emissions gas & oil'
               CASE( 3 )
                  NAME (T,58) = 'CH4-COL'
                  FNAME(T,58) = 'CH4 Emissions coal'
               CASE( 4 )
                  NAME (T,58) = 'CH4-LIV'
                  FNAME(T,58) = 'CH4 Emissions livestock'
               CASE( 5 )
                  NAME (T,58) = 'CH4-WST'
                  FNAME(T,58) = 'CH4 Emissions waste'
               CASE( 6 )
                  NAME (T,58) = 'CH4-BFL'
                  FNAME(T,58) = 'CH4 Emissions biofuel'
               CASE( 7 )
                  NAME (T,58) = 'CH4-RIC'
                  FNAME(T,58) = 'CH4 Emissions rice'
               CASE( 8 )
                  NAME (T,58) = 'CH4-OTA'
                  FNAME(T,58) = 'CH4 Emissions other anthro'
               CASE( 9 )
                  NAME (T,58) = 'CH4-BBN'
                  FNAME(T,58) = 'CH4 Emissions bioburn'
               CASE( 10 )
                  NAME (T,58) = 'CH4-WTL'
                  FNAME(T,58) = 'CH4 Emissions wetlands'
               CASE( 11 )
                  NAME (T,58) = 'CH4-SAB'
                  FNAME(T,58) = 'CH4 Emissions soil abs'
               CASE( 12 )
                  NAME (T,58) = 'CH4-OTN'
                  FNAME(T,58) = 'CH4 Emissions other natural'
               CASE( 13 )
                  NAME (T,58) = 'CH4-Loss'
                  FNAME(T,58) = 'CH4 Loss by Reaction w/ OH'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,58) = 'kg'
            INDEX(T,58) = T + ( SPACING * 50 )
            MOLC (T,58) = 1
            MWT  (T,58) = 1.6e-3
            SCALE(T,58) = 1e0
         ENDDO
      ENDIF

#if   defined( TOMAS )

      !-------------------------------------
      ! TOMAS 3-D rate (ND60)
      !-------------------------------------
      IF ( ND60 > 0 ) THEN
         
         ! # of tracers 
         NTRAC(60) = TMAX(60) 
         
         ! Loop over tracers
         DO T = 1, NTRAC(60)

            IF ( ( T .ge. IDTNK1 )       .and. 
     &           ( T .lt. IDTNK1+IBINS )      ) THEN
               UNIT(T,60) = 'no.'
            ELSE
               UNIT(T,60) = 'kg'
            ENDIF

            ! Define the rest of the quantities
            NAME (T,60) = Input_Opt%TRACER_NAME(T)
            FNAME(T,60) = 'Rates from ND60 TOMAS diagnostic'
            INDEX(T,60) = T + ( SPACING * 6 )
            MOLC (T,60) = 1
            MWT  (T,60) = 1e0
            SCALE(T,60) = 1e0
         ENDDO

      ENDIF

      !-------------------------------------
      ! TOMAS 3-D rate (ND61)                !(win, 7/9/09)
      !-------------------------------------
      IF ( ND61 > 0 ) THEN
         NTRAC(61) = 2

         ! Loop over tracers
         DO T = 1, NTRAC(61)
            
            ! Get name and long name for each tracer
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,61) = 'NUC10NM'
                  FNAME(T,61) = 'Particle Dp>10nm formation rate'
               CASE( 2 )
                  NAME (T,61) = 'NUC-1NM'
                  FNAME(T,61) = 'Nucleation rate (cluster size)'
               CASE DEFAULT
                  !Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,61) = 'cm-3s-1'
            INDEX(T,61) = T + ( SPACING * 7 )
            MOLC (T,61) = 1
            MWT  (T,61) = 1e0
            SCALE(T,61) = 1e0
         ENDDO
      ENDIF

#else

      !-------------------------------------      
      ! Wetland Fraction (ND60)
      !-------------------------------------      
      IF ( ND60 > 0 ) THEN

         ! Number of tracers
         NTRAC(60) = 1
         T         = 1

         NAME (T,60) = 'WETFRAC'
         FNAME(T,60) = 'Wetland Fraction'
         INDEX(T,60) = T + ( SPACING * 51 )
         MOLC (T,60) = 0 
         MWT  (T,60) = 0e0
         SCALE(T,60) = 1e0
         UNIT(T,60)  = 'unitless'
      ENDIF

#endif

      !-------------------------------------      
      ! SHIP DIAGNOSTICS (ND63)
      !-------------------------------------      
      IF ( DO_SAVE_DIAG63 ) THEN 

         ! Number of tracers
         NTRAC(63) = 5

         ! Loop over tracers
         DO T = 1, NTRAC(63)

            ! Define quantities
            INDEX(T,63) = T + ( SPACING * 59 )
            MOLC (T,63) = 0
            MWT  (T,63) = 0e0
            SCALE(T,63) = 1e0
            UNIT (T,63) = 'unitless'

            ! Get name, long-name (and others where necessary)
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,63) = 'FracNOx'
                  FNAME(T,63) = 'Fraction of NOx remaining'
               CASE( 2  )
                  NAME (T,63) = 'IntOPE'
                  FNAME(T,63) = 'Integrated OPE'
               CASE( 3  )
                  NAME (T,63) = 'PlumeNOx'
                  FNAME(T,63) = 'Ship plume NOx remaining'
                  UNIT (T,63) = 'molec/cm2/s'
               CASE( 4  )
                  NAME (T,63) = 'PlumeO3'
                  FNAME(T,63) = 'Ship plume Ox production'
                  UNIT (T,63) = 'molec/cm2/s'
               CASE( 5  )
                  NAME (T,63) = 'ShipNOx'
                  FNAME(T,63) = 'Ship NOx emissions'
                  UNIT (T,63) = 'molec/cm2/s'
               CASE DEFAULT
                  ! Nothing
            END SELECT

         ENDDO
      ENDIF

      !-------------------------------------      
      ! Family prod & loss (ND65)
      !-------------------------------------      
      IF ( DO_SAVE_PL ) THEN

         ! Number of P/L families
         NTRAC(65) = GET_NFAM()

         ! Loop over each P/L family
         DO T = 1, NTRAC(65)
            NAME (T,65) = GET_FAM_NAME( T ) 
            FNAME(T,65) = TRIM( NAME(T,65) ) // ' P/L family'
            INDEX(T,65) = T + ( SPACING * 17 )
            MOLC (T,65) = 1 
            MWT  (T,65) = GET_FAM_MWT( T )
            SCALE(T,65) = 1e0
            
            ! Unit for Tag Ox is kg/s, otherwise molec/cm3/s 
            IF ( Input_Opt%ITS_A_TAGOX_SIM ) THEN
               UNIT(T,65) = 'kg/s'
            ELSE
               UNIT(T,65) = 'molec/cm3/s'
            ENDIF


#if   defined( TOMAS )
            ! For sulfur tracking species, use 'kg S' unit and
            ! give a descriptive FNAME. (win, 8/9/09)
            IF ( NAME(T,65) == 'PSO21' ) THEN
               FNAME(T,65) = 'P(SO2) from DMS+OH (Online)'
               UNIT(T,65)  = 'kg S'
               MWT(T,65)   = 64e-3
            ELSE IF ( NAME(T,65) == 'PSO22' ) THEN
               FNAME(T,65) = 'P(SO2) from DMS+NO3 (Online)'
               UNIT(T,65)  = 'kg S'
               MWT(T,65)   = 64e-3
            ELSE IF ( NAME(T,65) == 'PSO4' ) THEN
               FNAME(T,65) = 'P(H2SO4) from SO2+OH (Online)'
               UNIT(T,65)  = 'kg S'
            ELSE IF ( NAME(T,65) == 'PMSA' ) THEN
               FNAME(T,65) = 'P(MSA) from DMS (Online)'
               UNIT(T,65) = 'kg S'
            ENDIF
#endif

         ENDDO
      ENDIF

      !-------------------------------------      
      ! 3-D GMAO met fields (ND66)
      ! also for timeseries diagnostics
      !-------------------------------------      
      IF ( ND66 > 0 .or. DO_TIMESERIES ) THEN

         ! Number of tracers
         NTRAC(66) = 6

         ! Loop over tracers
         DO T = 1, NTRAC(66)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1  )
                  NAME(T,66) = 'UWND'
                  UNIT(T,66) = 'm/s'
               CASE( 2  )
                  NAME(T,66) = 'VWND'
                  UNIT(T,66) = 'm/s'
               CASE( 3  )
                  NAME(T,66) = 'TMPU'
                  UNIT(T,66) = 'K'
               CASE( 4 )
                  NAME(T,66) = 'SPHU'
                  UNIT(T,66) = 'g/kg'
               CASE( 5 )
#if   defined( GEOS_4 )
                  NAME(T,66) = 'ZMMU'
                  UNIT(T,66) = 'Pa/s'
#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
                  NAME(T,66) = 'CMFMC'
                  UNIT(T,66) = 'kg/m2/s'
#else
                  NAME(T,66) = 'CLDMAS'
                  UNIT(T,66) = 'kg/m2/s'
#endif
               CASE( 6 )
                  NAME(T,66) = 'DTRAIN'
                  UNIT(T,66) = 'kg/m2/s'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            FNAME(T,66) = 'GMAO ' // TRIM( NAME(T,66) ) // ' field'
            INDEX(T,66) = T + ( SPACING * 12 )
            MOLC (T,66) = 1
            MWT  (T,66) = 0e0
            SCALE(T,66) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! 2-D GMAO met fields (ND67)
      !-------------------------------------      
      IF ( ND67 > 0 .or. DO_TIMESERIES ) THEN 

         ! Number of tracers
         NTRAC(67) = 23 ! (Lin, 03/31/09)

         ! Loop over tracers
         DO T = 1, NTRAC(67)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1  )
                  NAME(T,67) = 'HFLUX'
                  UNIT(T,67) = 'W/m2'
               CASE( 2  )
                  NAME(T,67) = 'RADSWG'
                  UNIT(T,67) = 'W/m2'
               CASE( 3  )
                  NAME(T,67) = 'PREACC'
                  UNIT(T,67) = 'W/m2'
               CASE( 4  )
                  NAME(T,67) = 'PRECON'
                  UNIT(T,67) = 'W/m2'
               CASE( 5  )
                  NAME(T,67) = 'TS'
                  UNIT(T,67) = 'K'
               CASE( 6  )
                  NAME(T,67) = 'RADSWT'
                  UNIT(T,67) = 'W/m2'
               CASE( 7  )
                  NAME(T,67) = 'USTAR'
                  UNIT(T,67) = 'm/s'
               CASE( 8  )
                  NAME(T,67) = 'Z0'
                  UNIT(T,67) = 'm'
               CASE( 9  )
                  NAME(T,67) = 'PBL'
#if   defined( GEOS_4 )
                  UNIT(T,67) = 'm'
#else
                  UNIT(T,67) = 'hPa'
#endif
               CASE( 10 )
                  NAME(T,67) = 'CLDFRC'
                  UNIT(T,67) = 'unitless'
               CASE( 11 )
                  NAME(T,67) = 'U10M'
                  UNIT(T,67) = 'm/s'
               CASE( 12 )
                  NAME(T,67) = 'V10M'
                  UNIT(T,67) = 'm/s'
               CASE( 13 )
                  NAME(T,67) = 'PS-PBL'
                  UNIT(T,67) = 'hPa'
               CASE( 14 )
                  NAME(T,67) = 'ALBD'
                  UNIT(T,67) = 'unitless'
               CASE( 15 )
                  NAME(T,67) = 'PHIS'
                  UNIT(T,67) = 'm'
               CASE( 16 )
                  NAME(T,67) = 'CLDTOP'
                  UNIT(T,67) = 'level'
               CASE( 17 )
                  NAME(T,67) = 'TROPP'
                  UNIT(T,67) = 'hPa'
               CASE( 18 )
                  NAME(T,67) = 'SLP'
                  UNIT(T,67) = 'hPa'
               CASE( 19 )
                  NAME(T,67) = 'TSKIN'
                  UNIT(T,67) = 'K'
               CASE( 20 )
                  NAME(T,67) = 'PARDF'
                  UNIT(T,67) = 'W/m2'
               CASE( 21 )
                  NAME(T,67) = 'PARDR'
                  UNIT(T,67) = 'W/m2'
               CASE( 22 )
                  NAME(T,67) = 'GWET'
                  UNIT(T,67) = 'unitless'
               CASE( 23 )
                  ! Add EFLUX (Lin, 05/16/08)
                  NAME(T,67) = 'EFLUX'
                  UNIT(T,67) = 'W/m2'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            FNAME(T,67) = 'GMAO ' // TRIM( NAME(T,67) ) // ' field'
            INDEX(T,67) = T + ( SPACING * 11 )
            MOLC (T,67) = 1
            MWT  (T,67) = 0e0
            SCALE(T,67) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Grid box heights and related 
      ! quantities (ND68) + timeseries
      !-------------------------------------      
      IF ( ND68 > 0 .or. DO_TIMESERIES ) THEN

         ! Number of tracers
         NTRAC(68) = 4

         ! Loop over tracers
         DO T = 1, NTRAC(68)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,68) = 'BXHEIGHT'
                  FNAME(T,68) = 'Grid box height'
                  UNIT (T,68) = 'm'
               CASE( 2 )
                  NAME (T,68) = 'AD'
                  FNAME(T,68) = 'Air mass in grid box'
                  UNIT (T,68) = 'kg'
               CASE( 3 )
                  NAME (T,68) = 'AVGW'
                  FNAME(T,68) = 'Mixing ratio of H2O vapor'
                  UNIT (T,68) = 'v/v'
               CASE( 4 )
                  NAME (T,68) = 'N(AIR)'
                  FNAME(T,68) = 'Number density of air'
                  UNIT (T,68) = 'molec/m3'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            INDEX(T,68) = T + ( SPACING * 24 )
            MOLC (T,68) = 1
            MWT  (T,68) = 0e0
            SCALE(T,68) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Grid box surface area (ND69)
      !-------------------------------------      
      IF ( ND69 > 0 ) THEN 
         T           = 1
         NTRAC(69)   = T
         NAME (T,69) = 'DXYP'
         FNAME(T,69) = 'Grid box surface area'
         UNIT (T,69) = 'm2'
         INDEX(T,69) = T + ( SPACING * 25 )
         MOLC (T,69) = 1
         MWT  (T,69) = 0e0
         SCALE(T,69) = 1e0
      ENDIF

      END SUBROUTINE INIT_TRACERINFO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gamap
!
! !DESCRIPTION: Subroutine INIT\_GAMAP allocates and initializes all 
!  module variables.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GAMAP( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE TIME_MOD,    ONLY : EXPAND_DATE, GET_NHMSb, GET_NYMDb
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE HDF_MOD,     ONLY : INIT_HDF
      USE HDF_MOD,     ONLY : HDFCATEGORY
      USE HDF_MOD,     ONLY : HDFDESCRIPT
      USE HDF_MOD,     ONLY : HDFNAME
      USE HDF_MOD,     ONLY : HDFFNAME
      USE HDF_MOD,     ONLY : HDFUNIT
      USE HDF_MOD,     ONLY : HDFMOLC
      USE HDF_MOD,     ONLY : HDFMWT
      USE HDF_MOD,     ONLY : HDFSCALE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  22 Apr 2005 - R. Yantosca - Initial version 
!  (1 ) Now add proper UNIT & SCALE for Rn/Pb/Be simulations (bmy, 5/11/05)
!  (2 ) Added HCN & CH3CN source & sink info for ND09 (bmy, 6/27/05)
!  (3 ) Bug fix: removed duplicate category names.  Updated for CO2-SRCE 
!        diagnostic.  Now references ND04 from "diag04_mod.f. 
!        (pns, bmy, 7/25/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now save MBO as tracer #5 for ND46 (tmf, bmy, 10/20/05)
!  (6 ) Now add categories CV-FLX-$, TURBMC-$, EW-FLX-$, NS-FLX-$, UP-FLX-$
!        which had been inadvertently omitted.  Also add OCEAN-HG category.
!        Rewrote do loop and case statement to add new diagnostics to ND03. 
!        Now make units of Hg tracers "pptv", not "ppbv".  Now remove 
!        restriction on printing out cloud mass flux in GEOS-4 for the ND66 
!        diagnostic.  Added new sea salt category. (cdh, eck, bmy, 4/6/06)
!  (7 ) Now references ND56 from "diag56_mod.f" (ltm, bmy, 5/5/06) 
!  (8 ) Now references ND42 from "diag42_mod.f".  Also updated for extra SOA
!        tracers in ND07 diagnostic. (dkh, bmy, 5/22/06)
!  (9 ) Updated ND36 for CH3I simulation (bmy, 7/25/06)
!  (10) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (11) Split into INIT_DIAGINFO, INIT_TRACERINFO for clarity (bmy, 9/28/06)
!  (12) Save output to HDF_MOD (amv, bmy, 12/18/09)
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS, NYMDb, NHMSb

      !=================================================================
      ! INIT_GAMAP begins here!
      !=================================================================

      ! Save from arguments to module variables
      DFILE = Input_Opt%GAMAP_DIAGINFO
      TFILE = Input_Opt%GAMAP_TRACERINFO

      ! Get starting date & time
      NYMDb = GET_NYMDb()
      NHMSb = GET_NHMSb()
 
      ! Replace any date/time tokens in the file names
      CALL EXPAND_DATE( DFILE, NYMDb, NHMSb )
      CALL EXPAND_DATE( TFILE, NYMDb, NHMSb )

      !=================================================================
      ! Allocate module arrays
      !=================================================================

      ALLOCATE( OFFSET( MAXCAT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OFFSET' )
      OFFSET = 0

      ALLOCATE( CATEGORY( MAXCAT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CATEGORY' )
      CATEGORY = ''

      ALLOCATE( DESCRIPT( MAXCAT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DESCRIPT' )
      DESCRIPT = ''

      ALLOCATE( NTRAC( MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NTRAC' )
      NTRAC = 0

      ALLOCATE( INDEX( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'INDEX' )      
      INDEX = 0

      ALLOCATE( MOLC( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOLC' )      
      MOLC = 0

      ALLOCATE( MWT( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MWT' )      
      MWT = 0.0

      ALLOCATE( SCALE( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SCALE' )      
      SCALE = 0.0

      ALLOCATE( NAME( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NAME' )      
      NAME = ''

      ALLOCATE( FNAME( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FNAME' )      
      FNAME = ''

      ALLOCATE( UNIT( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UNIT' )      
      UNIT = ''

      !=================================================================
      ! Initialize arrays for "diaginfo.dat" & "tracerinfo.dat" files
      !=================================================================

      ! Initialize arrays for "diaginfo.dat"
      CALL INIT_DIAGINFO

      ! Initialize arrays for "tracerinfo.dat"
      CALL INIT_TRACERINFO( am_I_Root, Input_Opt, RC )

      ! Store values for hdf output
      IF ( Input_Opt%LND50_HDF    .or. 
     &     Input_Opt%LND51_HDF    .or. 
     &     Input_Opt%LND51b_HDF ) THEN

         CALL INIT_HDF(MAXCAT, MAXTRACER, MAXDIAG)

         HDFCATEGORY(:) = CATEGORY(:)
         HDFDESCRIPT(:) = DESCRIPT(:)
         HDFNAME(:,:)   = NAME(:,:)
         HDFFNAME(:,:)  = FNAME(:,:)
         HDFUNIT(:,:)   = UNIT(:,:)
         HDFMOLC(:,:)   = MOLC(:,:)
         HDFMWT(:,:)    = MWT(:,:)
         HDFSCALE(:,:)  = SCALE(:,:)
         
      ENDIF

      END SUBROUTINE INIT_GAMAP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_gamap
!
! !DESCRIPTION: Subroutine CLEANUP\_GAMAP deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GAMAP
! 
! !REVISION HISTORY: 
!  25 Apr 2005 - R. Yantosca - Initial version
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !===================================================================
      ! CLEANUP_GAMAP begins here!
      !===================================================================
      IF ( ALLOCATED( CATEGORY ) ) DEALLOCATE( CATEGORY )
      IF ( ALLOCATED( DESCRIPT ) ) DEALLOCATE( DESCRIPT )
      IF ( ALLOCATED( FNAME    ) ) DEALLOCATE( FNAME    )
      IF ( ALLOCATED( INDEX    ) ) DEALLOCATE( INDEX    )
      IF ( ALLOCATED( MOLC     ) ) DEALLOCATE( MOLC     )
      IF ( ALLOCATED( MWT      ) ) DEALLOCATE( MWT      )
      IF ( ALLOCATED( NAME     ) ) DEALLOCATE( NAME     )
      IF ( ALLOCATED( NTRAC    ) ) DEALLOCATE( NTRAC    )
      IF ( ALLOCATED( OFFSET   ) ) DEALLOCATE( OFFSET   )
      IF ( ALLOCATED( SCALE    ) ) DEALLOCATE( SCALE    )
      IF ( ALLOCATED( UNIT     ) ) DEALLOCATE( UNIT     )

      END SUBROUTINE CLEANUP_GAMAP
!EOC
      END MODULE GAMAP_MOD
