!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_Br_mod
!
! !DESCRIPTION: Module GLOBAL\_Br\_MOD contains variables and routines for 
!  reading the global monthly mean Br concentration from disk. 
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_Br_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to store global monthly mean BR field
      REAL(fp), PUBLIC, ALLOCATABLE :: BR_TROP(:,:,:)
      REAL(fp), PUBLIC, ALLOCATABLE :: BR_STRAT(:,:,:)
      REAL(fp), PUBLIC, ALLOCATABLE :: BR_MERGE(:,:,:)

      ! Array to store global monthly mean BrO field
      REAL(fp), PUBLIC, ALLOCATABLE :: BRO_TROP(:,:,:)
      REAL(fp), PUBLIC, ALLOCATABLE :: BRO_STRAT(:,:,:)
      REAL(fp), PUBLIC, ALLOCATABLE :: BRO_MERGE(:,:,:)

      ! Array to store global monthly J-BrO field
      REAL(fp), PUBLIC, ALLOCATABLE :: J_BRO(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      ! Remove obsolete routine
      !PUBLIC :: GET_GLOBAL_Br_NEW   
      PUBLIC :: GET_GLOBAL_Br   
      PUBLIC :: INIT_GLOBAL_Br
      PUBLIC :: CLEANUP_GLOBAL_Br

! !REFERENCES
!  (1 ) Holmes, C. D., et al. (2006), Global lifetime of elemental mercury
!       against oxidation by atomic bromine in the free troposphere, Geophys.
!       Res. Lett., 33(20).
!  (2 ) Holmes, C.D., et al. (2010) Global atmospheric model for mercury 
!       including oxidation by bromine atoms, AC&P, 10, 12,037-12,057. 
!  (3 ) Parrella, J. et al. (2012), Tropospheric bromine chemistry:     
!       implications for present and pre-industrial ozone and mercury, ACP.
!
! !REVISION HISTORY:
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  19 Apr 2012 - E.S. Corbitt - Added LGCBROMINE to use GEOS-Chem bromine.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  14 Nov 2014 - M. Yannetti - Added PRECISION_MOD
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
! !IROUTINE: get_global_Br
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_Br reads global Br from 
!  binary punch files stored in the /data/ctm/GEOS\_MEAN directory.  This Br 
!  data is needed as oxidant for mercury chemistry.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_Br( am_I_Root, Input_Opt, 
     &                          State_Met, THISMONTH, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT
      USE BPCH2_MOD,          ONLY : GET_RES_EXT
      USE BPCH2_MOD,          ONLY : GET_TAU0
      USE BPCH2_MOD,          ONLY : READ_BPCH2
      USE CMN_SIZE_MOD
      USE CHEMGRID_MOD,       ONLY : GET_TPAUSE_LEVEL
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D
      USE OCEAN_MERCURY_MOD,  ONLY : LGCBROMINE     !eds 4/19/12
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)  :: THISMONTH   ! Current month
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  THIS IS A NEW VERSION OF THIS SUBROUTINE WHICH COMBINES Br CONCENTRATIONS
!  FROM MULTIPLE DATA SOURCES
!                                                                             .
!  NOTE: The bpch file input will eventually be replaced by netCDF input.
! 
! !REVISION HISTORY: 
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  (1 ) GET_GLOBAL_BR assumes that we are reading global BR data that occupies
!        all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) for IDL
!        regridding code which will produce the appropriate BR files.
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  06 Nov 2014 - R. Yantosca - Replace TRANSFER_3D_TROP with direct casts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                      :: I, J, L
      REAL*4                       :: ARRAY(IIPAR,JJPAR,LGLOB)
      REAL*4                       :: ARRAY2(IIPAR,JJPAR,LLCHEM)
      REAL(fp)                       :: XTAU
!      REAL*4                       :: TPAUSE(IIPAR,JJPAR)
      CHARACTER(LEN=255)           :: FILENAME
      INTEGER                      :: TPL
 
      ! Location of archived Br
      CHARACTER(LEN=255)           :: BR_DIR

      ! First time flag
      LOGICAL, SAVE                :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_BR begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate BR array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_BR
         FIRST = .FALSE.
      ENDIF

      ! EDS NOTE (4/19/12): need to move Br to directory, note the years used,
      ! and regrid for other resolutions
      BR_DIR = TRIM( Input_Opt%DATA_DIR ) // 'mercury_201007/'

      IF ( LGCBROMINE ) THEN

         !-----------------------------------------------------------------
         ! Read Br from GEOS-Chem
         !-----------------------------------------------------------------
         
         ! Filename
#if defined( GRID05x0666 ) && defined( NESTED_NA )
         ! Need to read emission file at nested-resolution
         FILENAME = TRIM( Input_Opt%DATA_DIR )   // 
     &              'mercury_201203/BrOx.GC.'    //
     &              GET_NAME_EXT()               // '.' // GET_RES_EXT()
#else
         FILENAME = TRIM( BR_DIR ) // 'BrOx.GC.' //
     &              GET_NAME_EXT() // '.'        // GET_RES_EXT()
#endif
   
         ! Get the TAU0 value for the start of the given month
         XTAU = GET_TAU0( THISMONTH, 1, 2007 ) 
   
         ! Echo some information to the standard output
         IF ( am_I_Root ) THEN
            WRITE( 6, 120 ) TRIM( FILENAME )
 120        FORMAT( '     - GET_GLOBAL_BR: Reading BR from: ', a )
         ENDIF
            
         CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 45,
     &                    XTAU,	    IIPAR,      JJPAR,
     &                    LGLOB,    ARRAY,      QUIET=.FALSE. )
   
         ! Assign data from ARRAY to the module variable BR
         CALL TRANSFER_3D( ARRAY, BR_TROP )

         ! Convert	ppbv ->	pptv
         BR_TROP  = BR_TROP * 1e+3_fp                         
   
         !-----------------------------------------------------------------
         ! Read BrO from GEOS-Chem
         !-----------------------------------------------------------------
   
         CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 46,
     &                    XTAU,	    IIPAR,      JJPAR,
     &                    LGLOB,    ARRAY,      QUIET=.TRUE. )
   
         ! Assign data from ARRAY to the module variable BRO
         CALL TRANSFER_3D( ARRAY, BRO_TROP )
       
         ! Convert ppbv -> pptv
         BRO_TROP  = BRO_TROP * 1e+3_fp 

      ELSE

         !-----------------------------------------------------------------
         ! Read Br from pTOMCAT biogenic bromocarbons
         !-----------------------------------------------------------------

         ! Filename
#if defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( Input_Opt%DATA_DIR )             //
     &              'mercury_201203/BrOx.TOMCAT_org.'      //
     &              GET_NAME_EXT() // '.' // GET_RES_EXT()
#else
         FILENAME = TRIM( BR_DIR ) // 'BrOx.TOMCAT_org.'   //
     &              GET_NAME_EXT() // '.' // GET_RES_EXT() 
#endif
         
         ! Get the TAU0 value for the start of the given month
         XTAU = GET_TAU0( THISMONTH, 1, 1985 ) 
         
         ! Echo some information to the standard output
         IF ( am_I_Root ) THEN
            WRITE( 6, 130 ) TRIM( FILENAME )
 130        FORMAT( '     - GET_GLOBAL_BR: Reading BR from: ', a )
         ENDIF

         ! Read BR data from the binary punch file
         CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 6,     
     &                    XTAU,      IIPAR,     JJPAR,      
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )
         
         ! Assign data from ARRAY to the module variable BR
         CALL TRANSFER_3D( ARRAY, BR_TROP )
         
         !-----------------------------------------------------------------
         ! Read BrO from pTOMCAT biogenic bromocarbons
         !-----------------------------------------------------------------
         
         ! Read BRO data from the binary punch file
         CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 7,     
     &                    XTAU,      IIPAR,     JJPAR,      
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )
         
         ! Assign data from ARRAY to the module variable BRO
         CALL TRANSFER_3D( ARRAY, BRO_TROP )

      ENDIF

      !-----------------------------------------------------------------
      ! Read Br from GMI for stratosphere
      !-----------------------------------------------------------------

      ! Filename
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = TRIM( Input_Opt%DATA_DIR )    //
     &           'mercury_201203/BrOx.GMI.'    //
     &           GET_NAME_EXT() // '.'         // GET_RES_EXT()
#else
      FILENAME = TRIM( BR_DIR ) // 'BrOx.GMI.' //
     &           GET_NAME_EXT() // '.'         // GET_RES_EXT()
#endif

      ! Get the TAU0 value for the start of the given month
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Echo some information to the standard output
      IF ( am_I_Root ) THEN
         WRITE( 6, 120 ) TRIM( FILENAME )
      ENDIF

      ! Read BR data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 6,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable BR
      CALL TRANSFER_3D( ARRAY, BR_STRAT )

      !-----------------------------------------------------------------
      ! Read BrO from GMI for stratosphere
      !-----------------------------------------------------------------

      ! Read BR data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 7,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable BRO
      CALL TRANSFER_3D( ARRAY, BRO_STRAT )

      !-----------------------------------------------------------------
      ! Use pTOMCAT or GEOS-Chem exclusively in the troposphere.
      ! In the stratosphere, use the greater value from either COMBO or
      ! the tropospheric model. COMBO source gases include CH3Br and 
      ! halons, while pTOMCAT and GEOS-Chem includes CH3Br and 
      ! shorter-lived gases.
      !-----------------------------------------------------------------

      BR_MERGE  = BR_TROP
      BRO_MERGE = BRO_TROP

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, TPL )
      DO I=1, IIPAR
      DO J=1, JJPAR      
         
         ! First layer in the stratosphere
         TPL = GET_TPAUSE_LEVEL( I, J, State_Met )

         BR_MERGE(I,J,TPL:LLPAR) = MERGE(
     &        BR_STRAT(I,J,TPL:LLPAR), 
     &        BR_TROP(I,J,TPL:LLPAR), 
     &        MASK=BR_STRAT(I,J,TPL:LLPAR)>BR_TROP(I,J,TPL:LLPAR) )

         BRO_MERGE(I,J,TPL:LLPAR) = MERGE(
     &        BRO_STRAT(I,J,TPL:LLPAR), 
     &        BRO_TROP(I,J,TPL:LLPAR), 
     &        MASK=BR_STRAT(I,J,TPL:LLPAR)>BR_TROP(I,J,TPL:LLPAR) )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !----------------------------------------------------------------
      ! OLD VERSION (cdh, 3/4/09)
      ! Previously switched at the tropopause, but pTOMCAT includes
      ! short-lived source gases that continue to dominate Br-x in 
      ! the lower stratosphere. So we want to use pTOMCAT in the lower
      ! strat
      !----------------------------------------------------------------
      !!-----------------------------------------------------------------
      !! Switch between COMBO and pTOMCAT at either the instantaneous
      !! tropopause or at the monthly-mean tropopause
      !!-----------------------------------------------------------------
      !
      !IF ( .NOT. LVARTROP ) THEN
      !   
      !   ! Read the monthly mean tropopause level
      !   ! Filename
      !   FILENAME = TRIM( BR_DIR ) // 'TPause.' // GET_NAME_EXT() // 
      ! &        '.'           // GET_RES_EXT()
      !
      !   ! Echo some information to the standard output
      !   WRITE( 6, 120 ) TRIM( FILENAME )
      !
      !
      !   ! Read BR data from the binary punch file
      !   CALL READ_BPCH2( FILENAME, 'TR-PAUSE', 1,     
      ! &        XTAU,      IIPAR,     JJPAR,      
      ! &        1,         TPAUSE,    QUIET=.TRUE. )
      !
      !
      !   ! Merge GMI above the tropopause with TOMCAT below it
      !   BR_MERGE = BR_TROP
      !   
      !   DO I=1, IIPAR
      !   DO J=1, JJPAR
      !
      !      ! First layer in the stratosphere
      !      TPL = INT( CEILING( TPAUSE(I,J) ) )
      !
      !      BR_MERGE( I, J, TPL:LLPAR ) = BR_STRAT( I, J, TPL:LLPAR )
      !
      !   ENDDO
      !   ENDDO
      !   
      !ENDIF
      !
      !----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Read J-BrO 
      ! Originally archived from GEOS-Chem Br-y simulation (prelim version)
      ! without Br-y wet deposition, using GEOS4 2002 meteorology.
      ! J-values are saved only for the troposphere. Values are daytime 
      ! averages.
      !-----------------------------------------------------------------

      ! Get the TAU0 value for the start of the given month
      XTAU = GET_TAU0( THISMONTH, 1, 2002 )

      ! Filename
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = TRIM( Input_Opt%DATA_DIR )             //
     &           'mercury_201203/jBrO.daytime.'         //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()
#else
      FILENAME = TRIM( BR_DIR ) // 'jBrO.daytime.'      //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT() 
#endif

      ! Echo some information to the standard output
      IF ( am_I_Root ) THEN
         WRITE( 6, 120 ) TRIM( FILENAME )
      ENDIF

      ! Read BR data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'LIFE-T=$', 19,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LLCHEM,    ARRAY2,    QUIET=.TRUE. )

      ! Assign data from ARRAY2 to the module variable BR
      J_BrO(:,:,1:LLCHEM_FIX) = ARRAY2(:,:,1:LLCHEM_FIX)

      END SUBROUTINE GET_GLOBAL_Br
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_Br
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_Br allocates and zeroes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_Br
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD 
! 
! !REVISION HISTORY: 
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_GLOBAL_BR begins here!
      !=================================================================

      !-------------------------------------
      ! Br Arrays
      !-------------------------------------

      ! Allocate BR_TROP array
      ALLOCATE( BR_TROP( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BR_TROP' )
      BR_TROP = 0e+0_fp

      ! Allocate BR_STRAT array
      ALLOCATE( BR_STRAT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BR_STRAT' )
      BR_STRAT = 0e+0_fp

      ! Allocate BR_MERGE array
      ALLOCATE( BR_MERGE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BR_MERGE' )
      BR_MERGE = 0e+0_fp


      !-------------------------------------
      ! BrO Arrays
      !-------------------------------------

      ! Allocate J_BrO array
      ALLOCATE( J_BrO( IIPAR, JJPAR, LLCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'J_BrO' )
      J_BrO = 0e+0_fp

      ! Allocate BrO_TROP array
      ALLOCATE( BrO_TROP( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BrO_TROP' )
      BrO_TROP = 0e+0_fp

      ! Allocate BrO_STRAT array
      ALLOCATE( BrO_STRAT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BrO_STRAT' )
      BrO_STRAT = 0e+0_fp

      ! Allocate BrO_MERGE array
      ALLOCATE( BrO_MERGE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BrO_MERGE' )
      BrO_MERGE = 0e+0_fp

      END SUBROUTINE INIT_GLOBAL_BR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_Br
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_Br deallocates module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_Br
! 
! !REVISION HISTORY: 
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_GLOBAL_BR begins here!
      !=================================================================
      IF ( ALLOCATED( BR_TROP   ) ) DEALLOCATE( BR_TROP   ) 
      IF ( ALLOCATED( BR_STRAT  ) ) DEALLOCATE( BR_STRAT  )  
      IF ( ALLOCATED( BR_MERGE  ) ) DEALLOCATE( BR_MERGE  ) 
      IF ( ALLOCATED( J_BrO     ) ) DEALLOCATE( J_BrO     ) 
      IF ( ALLOCATED( BrO_TROP  ) ) DEALLOCATE( BrO_TROP  ) 
      IF ( ALLOCATED( BrO_STRAT ) ) DEALLOCATE( BrO_STRAT )      
      IF ( ALLOCATED( BrO_MERGE ) ) DEALLOCATE( BrO_MERGE ) 

      END SUBROUTINE CLEANUP_GLOBAL_Br
!EOC
      END MODULE GLOBAL_Br_MOD
