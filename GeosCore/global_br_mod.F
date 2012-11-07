!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to store global monthly mean BR field
      REAL*8, PUBLIC, ALLOCATABLE :: BR_TROP(:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: BR_STRAT(:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: BR_MERGE(:,:,:)

      ! Array to store global monthly mean BrO field
      REAL*8, PUBLIC, ALLOCATABLE :: BRO_TROP(:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: BRO_STRAT(:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: BRO_MERGE(:,:,:)

      ! Array to store global monthly J-BrO field
      REAL*8, PUBLIC, ALLOCATABLE :: J_BRO(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      ! Remove obsolete routine
      !PUBLIC :: GET_GLOBAL_Br_NEW   
      PUBLIC :: GET_GLOBAL_Br   
      PUBLIC :: INIT_GLOBAL_Br
      PUBLIC :: CLEANUP_GLOBAL_Br
!
! !REVISION HISTORY:
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
#if defined( DEVEL )
      SUBROUTINE GET_GLOBAL_Br( THISMONTH, State_Met )
#else
      SUBROUTINE GET_GLOBAL_Br( THISMONTH )
#endif
!
! !USES:
!
      !USE LOGICAL_MOD,    ONLY : LVARTROP     ! Comment this out for now
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT
      USE BPCH2_MOD,      ONLY : GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0
      USE BPCH2_MOD,      ONLY : READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR ! cdh
      USE TRANSFER_MOD,   ONLY : TRANSFER_3D
      USE TRANSFER_MOD,   ONLY : TRANSFER_3D_TROP
      USE TROPOPAUSE_MOD, ONLY : GET_TPAUSE_LEVEL
#if defined( DEVEL )
      USE GIGC_State_Met_Mod, ONLY : MetState
#endif

      USE CMN_SIZE_MOD                 ! Size parameters

!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: THISMONTH   ! Current month
#if defined( DEVEL )
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
#endif
!
! !REMARKS:
!  THIS IS A NEW VERSION OF THIS SUBROUTINE WHICH COMBINES Br CONCENTRATIONS
!  FROM MULTIPLE DATA SOURCES
! 
! !REVISION HISTORY: 
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  (1 ) GET_GLOBAL_BR assumes that we are reading global BR data that occupies
!        all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) for IDL
!        regridding code which will produce the appropriate BR files.
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                      :: I, J, L
      REAL*4                       :: ARRAY(IIPAR,JJPAR,LGLOB)
      REAL*4                       :: ARRAY2(IIPAR,JJPAR,LLTROP)
      REAL*8                       :: XTAU
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

      ! Allocate BR array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_BR
         FIRST = .FALSE.
      ENDIF

      ! Get the TAU0 value for the start of the given month
      XTAU = GET_TAU0( THISMONTH, 1, 1985 ) !eds;  cdh-is this OK?

      BR_DIR = TRIM( DATA_DIR ) // 'mercury_201007/'

      !-----------------------------------------------------------------
      ! Read Br from pTOMCAT biogenic bromocarbons
      !-----------------------------------------------------------------

      ! Filename
      FILENAME = TRIM( BR_DIR ) // 'BrOx.TOMCAT_org.' //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT() 

      ! Echo some information to the standard output
      WRITE( 6, 120 ) TRIM( FILENAME )
 120  FORMAT( '     - GET_GLOBAL_BR: Reading BR from: ', a )

      ! Read BR data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 6,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable BR
      CALL TRANSFER_3D( ARRAY, BR_TROP )

      !-----------------------------------------------------------------
      ! Read BrO from pTOMCAT biogenic bromocarbons
      !-----------------------------------------------------------------

      ! Read BRO data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 7,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable BRO
      CALL TRANSFER_3D( ARRAY, BRO_TROP )

      !-----------------------------------------------------------------
      ! Read Br from GMI for stratosphere
      !-----------------------------------------------------------------

      ! Filename
      FILENAME = TRIM( BR_DIR ) // 'BrOx.GMI.' //
     &           GET_NAME_EXT() // '.'         // GET_RES_EXT()

      ! Echo some information to the standard output
      WRITE( 6, 120 ) TRIM( FILENAME )

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
      ! Use pTOMCAT exclusively in the troposphere.
      ! In the stratosphere, use the greater value from either COMBO or
      ! pTOMCAT. COMBO source gases include CH3Br and halons, while pTOMCAT 
      ! includes CH3Br and shorter-lived gases.
      !-----------------------------------------------------------------

      BR_MERGE  = BR_TROP
      BRO_MERGE = BRO_TROP

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, TPL )
      DO I=1, IIPAR
      DO J=1, JJPAR      
         
         ! First layer in the stratosphere
#if defined( DEVEL )
         TPL = GET_TPAUSE_LEVEL( I, J, State_Met )
#else
         TPL = GET_TPAUSE_LEVEL(I,J)
#endif

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
      FILENAME = TRIM( BR_DIR ) // 'jBrO.daytime.' //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT() 

      ! Echo some information to the standard output
      WRITE( 6, 120 ) TRIM( FILENAME )

      ! Read BR data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'LIFE-T=$', 19,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LLTROP,    ARRAY2,    QUIET=.TRUE. )

      ! Assign data from ARRAY2 to the module variable BR
      CALL TRANSFER_3D_TROP( ARRAY2, J_BrO )

      END SUBROUTINE GET_GLOBAL_Br
!EOC
! This routine is now obsolete
!!!!------------------------------------------------------------------------------
!!!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!!!------------------------------------------------------------------------------
!!!!BOP
!!!!
!!!! !IROUTINE: get_global_Br
!!!!
!!!! !DESCRIPTION: Subroutine GET\_GLOBAL\_Br reads global Br from binary 
!!!!  punch files stored in the /data/ctm/GEOS\_MEAN directory.  This Br data 
!!!!  is needed as oxidant for mercury chemistry.
!!!!\\
!!!!\\
!!!! !INTERFACE:
!!!!
!!!      SUBROUTINE GET_GLOBAL_Br( THISMONTH )
!!!!
!!!! !USES:
!!!!
!!!      USE BPCH2_MOD,    ONLY : GET_NAME_EXT
!!!      USE BPCH2_MOD,    ONLY : GET_RES_EXT
!!!      USE BPCH2_MOD,    ONLY : GET_TAU0
!!!      USE BPCH2_MOD,    ONLY : READ_BPCH2
!!!      USE TRANSFER_MOD, ONLY : TRANSFER_3D
!!!
!!!#     include "CMN_SIZE"                 ! Size parameters
!!!
!!!!
!!!! !INPUT PARAMETERS: 
!!!!
!!!      INTEGER, INTENT(IN) :: THISMONTH   ! Current month
!!!! 
!!!! !REVISION HISTORY: 
!!!!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!!!!  (1 ) GET_GLOBAL_BR assumes that we are reading global Br data that occupies
!!!!        all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) for IDL
!!!!        regridding code which will produce the appropriate Br files.
!!!!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!!!!EOP
!!!!------------------------------------------------------------------------------
!!!!BOC
!!!!
!!!! !LOCAL VARIABLES:
!!!!
!!!      INTEGER                      :: I, J, L
!!!      REAL*4                       :: ARRAY(IGLOB,JGLOB,LGLOB)
!!!      REAL*8                       :: XTAU
!!!      CHARACTER(LEN=255)           :: FILENAME
!!!
!!!      ! Location of archived Br
!!!      CHARACTER(LEN=30), PARAMETER :: BR_DIR = 
!!!     &                                '/home/cdh/GC/Archived-Br/'
!!!
!!!      ! First time flag
!!!      LOGICAL,           SAVE      :: FIRST = .TRUE. 
!!!
!!!      !=================================================================
!!!      ! GET_GLOBAL_BR begins here!
!!!      !=================================================================
!!!
!!!      ! Allocate BR array, if this is the first call
!!!      IF ( FIRST ) THEN
!!!         CALL INIT_GLOBAL_BR
!!!         FIRST = .FALSE.
!!!      ENDIF
!!!
!!!      ! Filename
!!!      FILENAME = TRIM( BR_DIR ) // 'Br_3Dglobal.' 
!!!     &           // GET_NAME_EXT() // '.' // GET_RES_EXT()
!!!
!!!      ! Echo some information to the standard output
!!!      WRITE( 6, 110 ) TRIM( FILENAME )
!!! 110  FORMAT( '     - GET_GLOBAL_BR: Reading Br, BrO from: ', a )
!!!
!!!      ! Time stamp for data
!!!      XTAU = GET_TAU0( THISMONTH, 1, 1985 ) 
!!!
!!!      ! Read BR data from the binary punch file
!!!      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 6,     
!!!     &                 XTAU,      IGLOB,     JGLOB,      
!!!     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )
!!!
!!!      ! Assign data from ARRAY to the module variable BR
!!!      CALL TRANSFER_3D( ARRAY, BR_MERGE )
!!!
!!!      END SUBROUTINE GET_GLOBAL_Br
!!!!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      BR_TROP = 0d0

      ! Allocate BR_STRAT array
      ALLOCATE( BR_STRAT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BR_STRAT' )
      BR_STRAT = 0d0

      ! Allocate BR_MERGE array
      ALLOCATE( BR_MERGE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BR_MERGE' )
      BR_MERGE = 0d0


      !-------------------------------------
      ! BrO Arrays
      !-------------------------------------

      ! Allocate J_BrO array
      ALLOCATE( J_BrO( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'J_BrO' )
      J_BrO = 0d0

      ! Allocate BrO_TROP array
      ALLOCATE( BrO_TROP( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BrO_TROP' )
      BrO_TROP = 0d0

      ! Allocate BrO_STRAT array
      ALLOCATE( BrO_STRAT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BrO_STRAT' )
      BrO_STRAT = 0d0

      ! Allocate BrO_MERGE array
      ALLOCATE( BrO_MERGE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BrO_MERGE' )
      BrO_MERGE = 0d0

      END SUBROUTINE INIT_GLOBAL_BR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
