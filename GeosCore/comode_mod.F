!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: comode_mod
!
! !DESCRIPTION: Module COMODE\_MOD contains allocatable arrays for SMVGEAR 
!  that were previously contained in common blocks in header file "comode.h".
!\\
!\\
! !INTERFACE: 
!
      MODULE COMODE_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      !======================================================================
      ! ABSHUM    : array for absolute humidity [H2O molec/cm3]
      ! AIRDENS   : array for air density [molec/cm3]
      ! CSPEC     : array of chemical species concentration [molec/cm3]
      ! CSPEC_FULL: array of chemical species for full potential troposphere
      ! CSUMA     : array for time of sunrise/sunset, from midnight [s]
      ! CSUMC     : array for temporary storage 
      ! ERADIUS   : array for aerosol or dust radii [cm]
      ! ERRMX2    : array for storing stiffness values 
      ! IXSAVE    : array of grid box longitude indices
      ! IYSAVE    : array of grid box latitude indices
      ! IZSAVE    : array of grid box altitude indices
      ! JLOP      : array of 1-D grid box indices
      ! PRESS3    : array for grid box pressure [mb]
      ! REMIS     : array for emissions from GEOS-CHEM [molec/cm3] 
      ! T3        : array for grid box temperature [K]
      ! TAREA     : array for surface area of aerosol or dust [cm2/cm3]
      ! VOLUME    : array for grid box volume [cm3]
      !======================================================================
      REAL*8,  ALLOCATABLE, PUBLIC :: ABSHUM(:) 
      REAL*8,  ALLOCATABLE, PUBLIC :: AIRDENS(:) 
      REAL*8,  ALLOCATABLE, PUBLIC :: CSPEC(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: CSUMA(:) 
      REAL*8,  ALLOCATABLE, PUBLIC :: CSUMC(:) 
      REAL*8,  ALLOCATABLE, PUBLIC :: ERADIUS(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: ERRMX2(:) 
      INTEGER, ALLOCATABLE, PUBLIC :: IXSAVE(:)
      INTEGER, ALLOCATABLE, PUBLIC :: IYSAVE(:)
      INTEGER, ALLOCATABLE, PUBLIC :: IZSAVE(:)
      INTEGER, ALLOCATABLE, PUBLIC :: JLOP(:,:,:)
      INTEGER, ALLOCATABLE, PUBLIC :: JLOP_PREVIOUS(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PRESS3(:)      
      REAL*8,  ALLOCATABLE, PUBLIC :: REMIS(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: T3(:)      
      REAL*8,  ALLOCATABLE, PUBLIC :: TAREA(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: VOLUME(:)      
      REAL*8,  ALLOCATABLE, PUBLIC :: WTAREA(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: WERADIUS(:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC :: CLEANUP_COMODE
      PUBLIC :: INIT_COMODE
!
! !REMARKS:
!  In case you were wondering, "comode" stands for:
!     "COMmon blocks: Ordinary Differential Equations"
!
! !REVISION HISTORY:
!  31 Aug 2000 - R. Yantosca - Initial version
!  (1 ) Now zero CSPEC after allocating memory (bmy, 9/8/00)
!  (2 ) Now declare more SMVGEAR arrays allocatable (bmy, 10/19/00)
!  (3 ) Updated comments (bmy, 9/4/01)
!  (4 ) Now make ERADIUS, TAREA 2-D arrays, for het chem (bmy, 11/15/01)
!  (5 ) DARSFCA is now obsolete, remove it.  Now allocate ERADIUS and
!        TAREA arrays to be of size (ITLOOP,NDUST+NAER).  (rvm, bmy, 2/27/02)
!  (5 ) Removed obsolete code from 2/02 (bmy, 4/15/02)
!  (6 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (7 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (8 ) Now add CSUMA, CSUMC, ERRMX2 arrays for SMVGEAR II (bmy, 7/18/03)
!  (9 ) Now also references "tracer_mod.f" (bmy, 9/28/04)
!  (10) Add WTAREA and WERADIUS variables. 
!       For SOA production from reactive uptake of dicarbonyls, 
!       archived WTAREA and WERADIUS should include dusts, 
!       but excludes BCPO and OCPO (tmf, ccc, 1/7/09)
!  (11) Added 3 *_KPP arrays (phs,ks,dhk, 09/15/09)
!  (12) Removed 3 *_KPP arrays (phs, 09/16/09)
!  21 Dec 2010 - R. Yantosca - Added ProTeX headers
!  04 Apr 2013 - R. Yantosca - Removed CSPEC_FULL (now in State_Chm)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete 
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
! !IROUTINE: init_comode
!
! !DESCRIPTION: Subroutine INIT\_COMODE allocates memory for allocatable 
!  arrays that were previously contained in common blocks in "comode.h".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_COMODE( am_I_Root, Input_Opt, RC )
!
! !USES:
!    
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD 
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  31 Aug 2000 - R. Yantosca - Initial version
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Cosmetic chagnes (bmy, 2/27/03)
!  (3 ) Now allocate CSUMA, CSUMC, ERRMX2; cosmetic changes (bmy, 7/18/03)
!  (4 ) Now allocate certain arrays for offline aerosol sim (bmy, 9/28/04)
!  21 Dec 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  13 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC arguments for GIGC
!  04 Apr 2013 - R. Yantosca - Removed CSPEC_FULL (now in State_Chm)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS
      LOGICAL :: IT_IS_AN_AEROSOL_SIM
      LOGICAL :: IT_IS_A_FULLCHEM_SIM

      !=================================================================
      ! INIT_COMODE begins here!
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 )
      ENDIF
 100  FORMAT( '     - INIT_COMODE: Allocating arrays for SMVGEAR...' )

      ! Assume success
      RC                   = GIGC_SUCCESS

      ! Copy logical fields from INPUT_OPT to local variables for use below
      IT_IS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM
      IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM

      !----------------------------------
      ! FULL CHEMISTRY SIMULATION
      !----------------------------------
      IF ( IT_IS_A_FULLCHEM_SIM ) THEN
      
         ALLOCATE( ABSHUM( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ABSHUM' )
         ABSHUM = 0d0
      
         ALLOCATE( AIRDENS( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIRDENS' )
         AIRDENS = 0d0      

         ALLOCATE( CSPEC( ITLOOP, IGAS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC' )
         CSPEC = 0d0

         ALLOCATE( CSUMA( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSUMA' )
         CSUMA = 0d0
      
         ALLOCATE( CSUMC( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSUMC' )
         CSUMC = 0d0

         ALLOCATE( ERADIUS( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ERADIUS' )
         ERADIUS = 0d0      

         ALLOCATE( ERRMX2( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ERRMX2' )
         ERRMX2 = 0d0
           
         ALLOCATE( IXSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IXSAVE' )
         IXSAVE = 0
      
         ALLOCATE( IYSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IYSAVE' )
         IYSAVE = 0
      
         ALLOCATE( IZSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IZSAVE' )
         IZSAVE = 0
      
         ALLOCATE( JLOP( ILONG, ILAT, IPVERT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'JLOP' )
         JLOP = 0
      
         !vartrop fix (dkh, 05/08/11)
         ALLOCATE( JLOP_PREVIOUS( ILONG, ILAT, IPVERT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'JLOP_PREVIOUS' )
         JLOP_PREVIOUS = 0

         ALLOCATE( PRESS3( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRESS3' )
         PRESS3 = 0d0
      
         ALLOCATE( REMIS( ITLOOP, MAXGL3 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'REMIS' )
         REMIS = 0d0
      
         ALLOCATE( T3( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'T3' )
         T3 = 0d0

         ALLOCATE( TAREA( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAREA' )
         TAREA = 0d0      
      
         ALLOCATE( VOLUME( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'VOLUME' )
         VOLUME = 0d0

         ALLOCATE( WTAREA( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'WTAREA' )
         WTAREA = 0d0      

         ALLOCATE( WERADIUS( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'WERADIUS' )
         WERADIUS = 0d0      

      ENDIF

      !----------------------------------
      ! OFFLINE AEROSOL SIMULATION
      !----------------------------------
      IF ( IT_IS_AN_AEROSOL_SIM ) THEN

         ALLOCATE( ABSHUM( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ABSHUM' )
         ABSHUM = 0d0
      
         ALLOCATE( AIRDENS( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIRDENS' )
         AIRDENS = 0d0      

         ALLOCATE( ERADIUS( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ERADIUS' )
         ERADIUS = 0d0      

         ALLOCATE( IXSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IXSAVE' )
         IXSAVE = 0
      
         ALLOCATE( IYSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IYSAVE' )
         IYSAVE = 0
      
         ALLOCATE( IZSAVE( ITLOOP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'IZSAVE' )
         IZSAVE = 0
      
         ALLOCATE( JLOP( ILONG, ILAT, IPVERT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'JLOP' )
         JLOP = 0

         !vartrop fix (dkh, 05/08/11)
         ALLOCATE( JLOP_PREVIOUS(ILONG, ILAT, IPVERT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'JLOP_PREVIOUS' )
         JLOP_PREVIOUS = 0

         ALLOCATE( TAREA( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAREA' )
         TAREA = 0d0      

         ALLOCATE( WTAREA( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'WTAREA' )
         WTAREA = 0d0

         ALLOCATE( WERADIUS( ITLOOP, NDUST+NAER ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'WERADIUS' )
         WERADIUS = 0d0      
         
      ENDIF

      END SUBROUTINE INIT_COMODE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_comode
!
! !DESCRIPTION: Subroutine CLEANUP\_COMODE deallocates memory from allocatable 
!  arrays that were previously contained in common blocks in "comode.h" 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_COMODE()
! 
! !REVISION HISTORY: 
!  31 Aug 2000 - R. Yantosca - Initial version
!  (1 ) Now deallocate CSPEC, CSUMA, ERRMX2; cosmetic changes (bmy, 7/18/03)
!  21 Dec 2010 - R. Yantosca - Added ProTeX headers
!  13 Nov 2012 - R. Yantosca - Cosmetic changes
!   4 Apr 2013 - R. Yantosca - Removed CSPEC_FULL (now in State_Chm)
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_COMODE begins here!
      !=================================================================
      IF ( ALLOCATED( ABSHUM        ) ) DEALLOCATE( ABSHUM        )
      IF ( ALLOCATED( AIRDENS       ) ) DEALLOCATE( AIRDENS       )
      IF ( ALLOCATED( CSPEC         ) ) DEALLOCATE( CSPEC         )
      IF ( ALLOCATED( CSUMA         ) ) DEALLOCATE( CSUMA         )
      IF ( ALLOCATED( CSUMC         ) ) DEALLOCATE( CSUMC         )
      IF ( ALLOCATED( ERADIUS       ) ) DEALLOCATE( ERADIUS       )
      IF ( ALLOCATED( ERRMX2        ) ) DEALLOCATE( ERRMX2        )
      IF ( ALLOCATED( IXSAVE        ) ) DEALLOCATE( IXSAVE        )
      IF ( ALLOCATED( IYSAVE        ) ) DEALLOCATE( IYSAVE        )
      IF ( ALLOCATED( IZSAVE        ) ) DEALLOCATE( IZSAVE        )
      IF ( ALLOCATED( JLOP          ) ) DEALLOCATE( JLOP          )
      IF ( ALLOCATED( JLOP_PREVIOUS ) ) DEALLOCATE( JLOP_PREVIOUS )
      IF ( ALLOCATED( PRESS3        ) ) DEALLOCATE( PRESS3        )     
      IF ( ALLOCATED( REMIS         ) ) DEALLOCATE( REMIS         )
      IF ( ALLOCATED( T3            ) ) DEALLOCATE( T3            )     
      IF ( ALLOCATED( TAREA         ) ) DEALLOCATE( TAREA         )
      IF ( ALLOCATED( VOLUME        ) ) DEALLOCATE( VOLUME        )  
      IF ( ALLOCATED( WTAREA        ) ) DEALLOCATE( WTAREA        )
      IF ( ALLOCATED( WERADIUS      ) ) DEALLOCATE( WERADIUS      )

      END SUBROUTINE CLEANUP_COMODE
!EOC
      END MODULE COMODE_MOD

