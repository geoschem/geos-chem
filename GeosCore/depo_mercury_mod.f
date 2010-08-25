!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: DEPO_MERCURY_MOD
!
! !DESCRIPTION: Module DEPO\_MERCURY\_MOD contains routines to handle
!  deposition fluxes for mercury. 
!
! !INTERFACE: 
!
      MODULE DEPO_MERCURY_MOD
!
! !USES:
! 
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: ADD_Hg2_DD
      PUBLIC :: ADD_Hg2_WD
      PUBLIC :: ADD_HgP_DD
      PUBLIC :: ADD_HgP_WD
      PUBLIC :: ADD_HG2_SNOWPACK
      PUBLIC :: RESET_HG_DEP_ARRAYS
      PUBLIC :: CHECK_DIMENSIONS
      PUBLIC :: READ_GTMM_RESTART
      PUBLIC :: MAKE_GTMM_RESTART
      PUBLIC :: UPDATE_DEP
      PUBLIC :: INIT_DEPO_MERCURY
      PUBLIC :: CLEANUP_DEPO_MERCURY
!
! !PUBLIC DATA MEMBERS:
!  
      PUBLIC :: DD_HG2, DD_HGP, WD_HG2, WD_HGP
      PUBLIC :: HG2mth_wd, HG0mth_dd, HG2mth_dd
      PUBLIC :: SNOW_HG
      PUBLIC :: LHGSNOW
      REAL*8,  ALLOCATABLE :: DD_Hg2(:,:,:)
      REAL*8,  ALLOCATABLE :: DD_HgP(:,:,:)
      REAL*8,  ALLOCATABLE :: WD_Hg2(:,:,:)
      REAL*8,  ALLOCATABLE :: WD_HgP(:,:,:)
      REAL*8,  ALLOCATABLE :: HG0mth_dd(:,:)
      REAL*8,  ALLOCATABLE :: HG2mth_dd(:,:)
      REAL*8,  ALLOCATABLE :: HG2mth_wd(:,:)
      REAL*8,  ALLOCATABLE :: SNOW_HG(:,:,:) !CDH Hg stored in snow+ice
      REAL*8,  ALLOCATABLE :: Hg0dryGEOS(:,:), HgIIdryGEOS(:,:), 
     &                        HgIIwetGEOS(:,:)
!
! !PRIVATE DATA MEMBERS:
!
      CHARACTER(LEN=255)   :: GTMM_RST_FILE

      LOGICAL :: LHGSNOW
!
! !REVISION HISTORY:
! 23 Apr 2010 - C. Carouge  - Initial version
!
!EOP
!------------------------------------------------------------------------------

      CONTAINS


!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_Hg2_DD
!
! !DESCRIPTION: Subroutine ADD_Hg2_DD computes the amount of Hg(II) dry deposited 
!  out of the atmosphere into the column array DD_Hg2. 
!  (sas, cdh, bmy, 1/19/05, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_Hg2_DD( I, J, N, DRY_Hg2)
!
! !USES
!
      USE TRACERID_MOD, ONLY : GET_Hg2_CAT
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I, J, N   ! GEOS-Chem long, lat and tracer index
      REAL*8,  INTENT(IN)   :: DRY_Hg2   ! Hg(II) dry deposited out of the 
                                         ! atmosphere [kg]
!
! !REVISION HISTORY:
!  (1 ) DD_Hg2 is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_Hg2_CAT to return the Hg category #. (cdh, bmy, 3/28/06)
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: NN
      
      !=================================================================
      ! ADD_Hg2_DD begins here!
      !=================================================================

      ! Get the index for DD_Hg2 based on the tracer number
      NN = GET_Hg2_CAT( N )

      ! Store dry deposited Hg(II) into DD_Hg2 array
      IF ( NN > 0 ) THEN
         DD_Hg2(I,J,NN) = DD_Hg2(I,J,NN) + DRY_Hg2
        
      ENDIF
      
     
      ! Return to calling program
      END SUBROUTINE ADD_Hg2_DD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_Hg2_WD
!
! !DESCRIPTION: Subroutine ADD_Hg2_WD computes the amount of Hg(II) wet scavenged 
!  out of the atmosphere into the column array WD_Hg2. 
!  (sas, cdh, bmy, 1/19/05, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_Hg2_WD( I, J, N, WET_Hg2 )
!
! !USES
!
      USE TRACERID_MOD, ONLY : GET_Hg2_CAT
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I, J, N   ! GEOS-Chem long, lat and tracer index
      REAL*8,  INTENT(IN)   :: WET_Hg2   ! Hg(II) scavenged out of the 
                                         ! atmosphere [kg]
!
! !REVISION HISTORY:
!  (1 ) WD_Hg2 is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_Hg2_CAT to return the Hg category #. (cdh, bmy, 3/28/06)
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: NN

      !=================================================================
      ! ADD_Hg2_WD begins here!
      !=================================================================

      ! Get Hg2 category number
      NN = GET_Hg2_CAT( N ) 
     
      ! Store wet deposited Hg(II) into WD_Hg2 array
      IF ( NN > 0 ) THEN
         WD_Hg2(I,J,NN) = WD_Hg2(I,J,NN) + WET_Hg2
         
      ENDIF

      ! Return to calling program
      END SUBROUTINE ADD_Hg2_WD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_HgP_DD
!
! !DESCRIPTION: Subroutine ADD_HgP_DD computes the amount of HgP dry deposited 
!  out of the atmosphere into the column array DD_HgP. 
!  (sas, cdh, bmy, 1/19/05, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_HgP_DD( I, J, N, DRY_HgP )
!
! !USES
!
      USE TRACERID_MOD, ONLY : GET_HgP_CAT
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I, J, N   ! GEOS-Chem long, lat and tracer index
      REAL*8,  INTENT(IN)   :: DRY_HgP   ! HgP dry deposited out of the 
                                         ! atmosphere [kg]
!
! !REVISION HISTORY:
!  (1 ) DD_HgP is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_HgP_CAT to return the Hg category #. (cdh, bmy, 3/28/06)
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      INTEGER               :: NN

      !=================================================================
      ! ADD_HgP_DD begins here!
      !=================================================================
      
      ! Get the index for DD_Hg2 based on the tracer number
      NN = GET_HgP_CAT( N )

      ! Store dry deposited Hg(II) into DD_Hg2 array
      IF ( NN > 0 ) THEN
         DD_HgP(I,J,NN) = DD_HgP(I,J,NN) + DRY_HgP
        
      ENDIF

      ! Return to calling program
      END SUBROUTINE ADD_HgP_DD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_HgP_WD
!
! !DESCRIPTION: Subroutine ADD_HgP_WD computes the amount of HgP wet scavenged
!  out of the atmosphere into the column array WD_HgP. 
!  (sas, cdh, bmy, 1/19/05, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_HgP_WD( I, J, N, WET_HgP )
!
! !USES
!
      USE TRACERID_MOD, ONLY : GET_HgP_CAT
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I, J, N   ! GEOS-Chem long, lat and tracer index
      REAL*8,  INTENT(IN)   :: WET_HgP   ! HgP scavenged out of the 
                                         ! atmosphere [kg]
!
! !REVISION HISTORY:
!  (1 ) WD_HgP is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_HgP_CAT to return the Hg category #. (cdh, bmy, 3/28/06)
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: NN

      !=================================================================
      ! ADD_Hg2_WD begins here!
      !=================================================================
      
      ! Get Hg2 category number
      NN = GET_HgP_CAT( N ) 

       ! Store wet deposited HgP into WD_HgP array
      IF ( NN > 0 ) THEN
         WD_HgP(I,J,NN) = WD_HgP(I,J,NN) + WET_HgP
        
      ENDIF
      
      ! Return to calling program
      END SUBROUTINE ADD_HgP_WD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_HG2_SNOWPACK
!
! !DESCRIPTION: Subroutine RESET_HG_DEP_ARRAYS adds Hg2 deposition to snowpack.
!  (cdh, 9/2/08, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_HG2_SNOWPACK( I, J, N, DEP_Hg2 )
!
! !USES:
!
      USE DAO_MOD,           ONLY : SNOW, SNOMAS 
      USE DAO_MOD,           ONLY : IS_ICE
      USE TRACERID_MOD,      ONLY : GET_Hg2_CAT, GET_HgP_CAT
      USE TRACERID_MOD,      ONLY : IS_Hg2, IS_HgP

#     include 'define.h'
!
! !INPUT PARAMETERS:
!
      ! Arguments as input
      INTEGER, INTENT(IN)   :: I, J, N
      REAL*8,  INTENT(IN)   :: Dep_Hg2
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge  - Moved from mercury_mod.f to depo_mercury_mod.f
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      REAL*8                :: SNOW_HT
      INTEGER               :: NN

      !=================================================================
      ! ADD_HG2_SNOWPACK begins here!
      !=================================================================
      
      ! Return if snowpack model is disabled
      IF (.NOT. LHGSNOW) RETURN

      IF ( IS_Hg2( N ) ) THEN
         ! Get Hg2 category number
         NN = GET_Hg2_CAT( N ) 
      ELSE IF ( IS_HgP( N ) ) THEN
         ! Get HgP category number
         NN = GET_HgP_CAT( N ) 
      ENDIF

#if   defined( GEOS_5 ) || defined( MERRA )
      ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
      SNOW_HT = SNOMAS(I,J)
#else
      ! GEOS1-4 snow heigt (water equivalent) in mm
      SNOW_HT = SNOW(I,J)
#endif 

      ! Check if there is snow on the ground, or if this is sea ice
      IF ( (SNOW_HT > 1d0) .OR. (IS_ICE(I,J)) ) THEN
    
         IF (DEP_HG2<0d0) THEN
            WRITE(6,'(3I6,2G12.4)') I,J,NN,DEP_HG2,SNOW_HG(I,J,NN)
         ENDIF

         SNOW_HG(I,J,NN) = SNOW_HG(I,J,NN) + MAX( DEP_HG2, 0D0 )

      ENDIF

      ! Return to calling program
      END SUBROUTINE ADD_HG2_SNOWPACK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RESET_HG_DEP_ARRAYS
!
! !DESCRIPTION: Subroutine RESET_HG_DEP_ARRAYS resets the wet and dry 
!  deposition arrays for Hg(II) and Hg(p) to zero. This allows us to call 
!  OCEAN_MERCURY_FLUX and LAND_MERCURY_FLUX in any order in MERCURY_MOD. 
!  (cdh, 9/2/08, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE RESET_HG_DEP_ARRAYS
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! Reset deposition arrays.
      DD_Hg2 = 0d0
      WD_Hg2 = 0d0
      DD_HgP = 0d0
      WD_HgP = 0d0

      END SUBROUTINE RESET_HG_DEP_ARRAYS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: MAKE\_GTMM\_RESTART
!
! !DESCRIPTION: MAKE\_GTMM\_RESTART writes a GTMM restart file with deposition
!  fluxes and store deposition fluxes for continuous runs. (ccc, 9/15/09) 
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE MAKE_GTMM_RESTART( NYMD, NHMS, TAU )
! 
! !USES:
!
      USE BPCH2_MOD
      USE DIAG_MOD,      ONLY : AD39, AD44, AD38
      USE DIRECTORY_MOD, ONLY : RUN_DIR
      USE FILE_MOD,      ONLY : IU_FILE
      USE GRID_MOD,      ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE TRACERID_MOD,  ONLY : ID_Hg0, ID_Hg2, ID_Hg_tot
      USE TIME_MOD,      ONLY : GET_CT_DYN, GET_CT_CHEM
    
#     include "CMN_SIZE"          ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: NYMD    ! Year-Month-Date
      INTEGER, INTENT(IN)   :: NHMS    ! and Hour-Min-Sec for which to create 
                                       ! a restart file 
      REAL*8,  INTENT(IN)   :: TAU     ! GEOS-CHEM TAU value corresponding to 
                                       ! NYMD, NHMS
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: HALFPOLAR, CENTER180
      INTEGER               :: IFIRST,    JFIRST,   LFIRST
      INTEGER               :: N,         NN
      REAL*8                :: TS_DYN,    TS_CHEM
      REAL*4                :: LONRES,    LATRES,   ARRAY(IGLOB,JGLOB,1)
      CHARACTER(LEN=20)     :: MODELNAME
      CHARACTER(LEN=40)     :: CATEGORY,  UNIT,     RESERVED
      CHARACTER(LEN=255)    :: FILENAME

      !=================================================================
      ! MAKE_GTMM_RESTART begins here!
      !=================================================================
    
      ! Initialize values
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LFIRST    = 1
      HALFPOLAR = GET_HALFPOLAR()
      CENTER180 = 1
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      CATEGORY  = 'DRYD-FLX'
      RESERVED  = ''
      UNIT      = 'molec/cm2/s'

      ! Expand date in filename
      FILENAME  = TRIM( RUN_DIR ) // TRIM( GTMM_RST_FILE )
      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_RESTART_FILE: Writing ', a )

      ! Open BPCH file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME )

      !---------------------------
      ! Total Hg(0) dry deposition
      !---------------------------
      ARRAY(:,:,1) = HG0mth_dd
    
      CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,         
     &     HALFPOLAR, CENTER180, CATEGORY, N,                     
     &     UNIT,      TAU,       TAU,      RESERVED,              
     &     IIPAR,     JJPAR,     1,        IFIRST,                
     &     JFIRST,    LFIRST,    ARRAY(:,:,1) )

      !---------------------------
      ! Hg(II) dry deposition
      !---------------------------
      ARRAY(:,:,1) = HG2mth_dd

      CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,         
     &     HALFPOLAR, CENTER180, CATEGORY, N,                     
     &     UNIT,      TAU,       TAU,      RESERVED,              
     &     IIPAR,     JJPAR,     1,        IFIRST,                
     &    JFIRST,    LFIRST,    ARRAY(:,:,1) )

      !---------------------------
      ! Hg(II) wet deposition
      !---------------------------
      CATEGORY  = 'WETDLS-$'
      ARRAY(:,:,1) = HG2mth_wd

      CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,        
     &     HALFPOLAR, CENTER180, CATEGORY, N,                    
     &     UNIT,      TAU,       TAU,      RESERVED,             
     &     IIPAR,     JJPAR,     1,        IFIRST,               
     &     JFIRST,    LFIRST,    ARRAY(:,:,1) )
    
      ! Close file
      CLOSE( IU_FILE )
    
      ! Return to calling program
      END SUBROUTINE MAKE_GTMM_RESTART
  
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: READ\_GTMM\_RESTART
!
! !DESCRIPTION: Subroutine READ\_GTMM\_RESTART reads dry and wet deposition 
!  for mercury from GTMM restart. (ccc, 9/15/09)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE READ_GTMM_RESTART( YYYYMMDD, HHMMSS, 
     &                            Hg0dryGEOS, HgIIdryGEOS, HgIIwetGEOS )
! 
! !USES:
!
      USE BPCH2_MOD,     ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_MOD, ONLY : RUN_DIR
      USE ERROR_MOD,     ONLY : DEBUG_MSG
      USE FILE_MOD,      ONLY : IU_FILE,     IOERROR
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE TRACER_MOD,    ONLY : STT,         TRACER_NAME, TRACER_MW_G
      USE TRACERID_MOD,  ONLY : GET_Hg0_CAT, GET_Hg2_CAT, N_Hg_CATS
      USE TRACERID_MOD,  ONLY : ID_Hg0,      ID_Hg2

#     include "CMN_SIZE"
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)   :: YYYYMMDD, HHMMSS
!
! !OUTPUT PARAMETERS:
!
      REAL*8, DIMENSION(IIPAR, JJPAR)   :: Hg0dryGEOS
      REAL*8, DIMENSION(IIPAR, JJPAR)   :: HgIIdryGEOS
      REAL*8, DIMENSION(IIPAR, JJPAR)   :: HgIIwetGEOS     
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: IOS, I, J, L, NN, N_gtmm
      INTEGER               :: NCOUNT(NNPAR) 
      REAL*4                :: FLUX(IIPAR,JJPAR)
      CHARACTER(LEN=255)    :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER               :: NI,        NJ,      NL
      INTEGER               :: IFIRST,    JFIRST,  LFIRST
      INTEGER               :: NTRACER,   NSKIP
      INTEGER               :: HALFPOLAR, CENTER180
      REAL*4                :: LONRES,    LATRES
      REAL*8                :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)     :: MODELNAME
      CHARACTER(LEN=40)     :: CATEGORY
      CHARACTER(LEN=40)     :: UNIT     
      CHARACTER(LEN=40)     :: RESERVED
  
      !=================================================================
      ! READ_GTMM_RESTART begins here!
      !=================================================================

      ! Copy input file name to a local variable
      FILENAME = TRIM( RUN_DIR ) // TRIM( GTMM_RST_FILE )
    
      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
    
      ! Echo some input to the screen
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100   ) 
      WRITE( 6, 110   ) TRIM( FILENAME )
 100  FORMAT( 'G T M M  H g   R E S T A R T   F I L E   I N P U T' )
 110  FORMAT( /, 'READ_GTMM_RESTART: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )
    
      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      DO 
         READ( IU_FILE, IOSTAT=IOS )                              
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
       
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT
       
         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_rst:1' )
       
         READ( IU_FILE, IOSTAT=IOS )                               
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,  
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,    
     &        NSKIP
       
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_rst:2' )
       
         READ( IU_FILE, IOSTAT=IOS )                               
     &        ( ( FLUX(I,J), I=1,NI ), J=1,NJ )
       
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_rst:3' )
       
         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================
       
         ! Process dry deposition data 
         IF ( CATEGORY(1:8) == 'DRYD-FLX' ) THEN 
          
            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            CALL CHECK_DIMENSIONS( NI, NJ, NL )
          
            ! Save into arrays
            IF ( ANY( ID_Hg0 == NTRACER ) ) THEN
             
               !----------
               ! Hg(0)
               !----------
             
               ! Get the Hg category #
               NN              = GET_Hg0_CAT( NTRACER )
             
               ! Store ocean Hg(0) in Hg0aq array
               Hg0dryGEOS(:,:)   = FLUX(:,:)
             
               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1
             
            ELSE IF ( ANY( ID_Hg2 == NTRACER ) ) THEN
             
               !----------
               ! Hg(II)
               !----------
             
               ! Get the Hg category #
               NN              = GET_Hg2_CAT( NTRACER )
             
               ! Store ocean Hg(II) in Hg2_aq array
               HgIIdryGEOS(:,:)   = FLUX(:,:)
             
               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1
             
            ENDIF
         ENDIF
       
         ! Process wet deposition data
         IF ( CATEGORY(1:8) == 'WETDLS-$' ) THEN 

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! %%%% KLUDGE: CHECK_DIMENSIONS only works for NL=1 !!!!
            ! And we are only interested by the surface flux...
            NL = 1
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            CALL CHECK_DIMENSIONS( NI, NJ, NL )
          
            IF ( ANY( ID_Hg2 == NTRACER ) ) THEN
             
               !----------
               ! Hg(II)
               !----------
             
               ! Get the Hg category #
               NN              = GET_Hg2_CAT( NTRACER )
             
               ! Store ocean Hg(II) in Hg2_aq array
               HgIIwetGEOS(:,:)   = FLUX(:,:)
             
               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1
             
            ENDIF
         ENDIF
      ENDDO
    
      ! Close file
      CLOSE( IU_FILE )      
    
      END SUBROUTINE READ_GTMM_RESTART
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: UPDATE_DEP
!
! !DESCRIPTION: Subroutine UPDATE\_DEP update the monthly average for wet and 
!  dry deposition of Hg0 and Hg2 for mercury from GTMM restart. (ccc, 6/4/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE UPDATE_DEP( NN )
! 
! !USES:
!
      USE DIAG_MOD,     ONLY : AD38,   AD39,   AD44
      USE LOGICAL_MOD,  ONLY : LGTMM
      USE TIME_MOD,     ONLY : GET_CT_DYN,  GET_CT_CHEM 
      USE TRACERID_MOD, ONLY : IDTHg0, IDTHg2
!
! !INPUT PARAMETERS:
!
      INTEGER :: NN    ! Hg2 ID for wet deposition
!
! !REVISION HISTORY:
! 
!  4 June 2010  - C. Carouge  - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N
      REAL*8  :: SCALEDYN, SCALECHEM

      !=================================================================
      ! UPDATE_DEP begins here!
      !=================================================================

      ! counter variables 
      SCALEDYN   = DBLE( GET_CT_DYN()  ) + 1d-32
      SCALECHEM  = DBLE( GET_CT_CHEM() ) + 1d-32

      ! Hg2 total wet deposition at the surface
      HG2mth_wd = HG2mth_wd + ( SUM(AD38(:,:,:,NN), DIM=3) + 
     &                          SUM(AD39(:,:,:,NN), DIM=3) )
     &                      / SCALEDYN

      ! Hg0 total dry deposition at the surface
      N = IDTHg0
      HG0mth_dd = HG0mth_dd + AD44(:,:,N,1) / SCALECHEM

      ! Hg2 total dry deposition at the surface
      N = IDTHg2
      HG2mth_dd = HG2mth_dd + AD44(:,:,N,1) / SCALECHEM

      ! Return to calling program
      END SUBROUTINE UPDATE_DEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: CHECK\_DIMENSIONS
!
! !DESCRIPTION: Subroutine CHECK\_DIMENSIONS makes sure that the dimensions of 
!  the Hg restart file extend to cover the entire grid.
!  (sas, cdh, bmy, 3/28/06)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE CHECK_DIMENSIONS( NI, NJ, NL ) 
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
#     include "CMN_SIZE"
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: NI, NJ, NL
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! CHECK_DIMENSIONS begins here!
      !=================================================================

      ! Error check longitude dimension: NI must equal IIPAR
      IF ( NI /= IIPAR ) THEN
         WRITE( 6, 100 ) 
 100     FORMAT( 'ERROR reading in Hg restart file', /
     &           'Wrong number of longitudes encountered', /
     &           'STOP in CHECK_DIMENSIONS ("ocean_mercury_mod.f")' )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Error check latitude dimension: NJ must equal JJPAR
      IF ( NJ /= JJPAR ) THEN
         WRITE( 6, 110 ) 
 110     FORMAT( 'ERROR reading in Hg restart file', /
     &           'Wrong number of longitudes encountered', /
     &           'STOP in CHECK_DIMENSIONS ("ocean_mercury_mod.f")' )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF
      
      ! Error check vertical dimension: NL must equal LLPAR
      IF ( NL /= 1 ) THEN
         WRITE( 6, 120 ) 
 120     FORMAT( 'ERROR reading in Hg restart file', /
     &           'Wrong number of longitudes encountered', /
     &           'STOP in CHECK_DIMENSIONS ("ocean_mercury_mod.f")' )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK_DIMENSIONS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_DEPO_MERCURY
!
! !DESCRIPTION: Subroutine INIT\_DEPO\_MERCURY initialize deposition arrays for
!  mercury (ccc, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE INIT_DEPO_MERCURY( THIS_Hg_RST_FILE )
!
! !USES
!
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE LOGICAL_MOD,  ONLY : LGTMM
      USE TRACERID_MOD, ONLY : N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters
!
! !INPUT PARAMETERS:
!
      ! Name of the GTMM restart file
      CHARACTER(LEN=*), INTENT(IN) :: THIS_Hg_RST_FILE
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge   - Moved arrays allocation from ocean_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                      :: AS

      ! GTMM restart file name
      GTMM_RST_FILE = THIS_Hg_RST_FILE

      ! Allocate arrays
      ALLOCATE( DD_Hg2( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DD_Hg2' )
      DD_Hg2 = 0d0

      ALLOCATE( DD_HgP( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DD_HgP' )
      DD_HgP = 0d0

      ALLOCATE( WD_Hg2( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WD_Hg2' )
      WD_Hg2 = 0d0

      ALLOCATE( WD_HgP( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WD_HgP' )
      WD_HgP = 0d0

      ! CDH for snowpack
      ALLOCATE( SNOW_HG( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SNOW_HG' )
      SNOW_HG = 0d0

      IF ( LGTMM ) THEN
         ALLOCATE( HG0mth_dd( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'HG0mth_dd' )
         HG0mth_dd = 0d0
         
         ALLOCATE( HG2mth_dd( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'HG2mth_dd' )
         HG2mth_dd = 0d0
         
         ALLOCATE( HG2mth_wd( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'HG2mth_wd' )
         HG2mth_wd = 0d0
      ENDIF

      END SUBROUTINE INIT_DEPO_MERCURY
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_DEPO_MERCURY
!
! !DESCRIPTION: Subroutine CLEANUP\_DEPO\_MERCURY deallocate all arrays
!  (ccc, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE CLEANUP_DEPO_MERCURY
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge   - Moved from ocean_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IF ( ALLOCATED( DD_Hg2      ) ) DEALLOCATE( DD_Hg2      )
      IF ( ALLOCATED( DD_HgP      ) ) DEALLOCATE( DD_HgP      )
      IF ( ALLOCATED( WD_Hg2      ) ) DEALLOCATE( WD_Hg2      )
      IF ( ALLOCATED( WD_HgP      ) ) DEALLOCATE( WD_HgP      )
      IF ( ALLOCATED( SNOW_HG     ) ) DEALLOCATE( SNOW_HG     )!CDH for snowpack
      IF ( ALLOCATED( HG0mth_dd   ) ) DEALLOCATE( HG0mth_dd   )
      IF ( ALLOCATED( HG2mth_dd   ) ) DEALLOCATE( HG2mth_dd   )
      IF ( ALLOCATED( HG2mth_wd   ) ) DEALLOCATE( HG2mth_wd   )
      
      END SUBROUTINE CLEANUP_DEPO_MERCURY

      END MODULE DEPO_MERCURY_MOD
!EOC
