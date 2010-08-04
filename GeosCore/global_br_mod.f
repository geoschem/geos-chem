
      MODULE GLOBAL_BR_MOD
!
!******************************************************************************
!  Module GLOBAL_BR_MOD contains variables and routines for reading the
!  global monthly mean Br concentration from disk. 
!  (bmy, cdh, 7/28/00, 10/3/05, 7/6/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) BR (REAL*8)       : stores global monthly mean BR field
!  
!  Module Routines:
!  ============================================================================
!  (1 ) GET_BR            : Wrapper for GET_GLOBAL_BR
!  (2 ) GET_GLOBAL_BR     : Reads global monthly mean BR from disk
!  (3 ) INIT_GLOBAL_BR    : Allocates & initializes the BR array
!  (4 ) CLEANUP_GLOBAL_BR : Deallocates the BR array
!
!  GEOS-CHEM modules referenced by global_nox_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f  : Module containing routines for binary punch file I/O
!  (2 ) error_mod.f  : Module containing NaN and other error-check routines
!
!  NOTES:
!  (1 ) Copied from global_oh_mod.f (cdh, 7/5/06)
!******************************************************************************
!     
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Array to store global monthly mean BR field
      REAL*8, ALLOCATABLE :: BR_TROP(:,:,:)
      REAL*8, ALLOCATABLE :: BR_STRAT(:,:,:)
      REAL*8, ALLOCATABLE :: BR_MERGE(:,:,:)

      ! Array to store global monthly mean BrO field
      REAL*8, ALLOCATABLE :: BRO_TROP(:,:,:)
      REAL*8, ALLOCATABLE :: BRO_STRAT(:,:,:)
      REAL*8, ALLOCATABLE :: BRO_MERGE(:,:,:)

      ! Array to store global monthly J-BrO field
      REAL*8, ALLOCATABLE :: J_BRO(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_GLOBAL_BR_NEW( THISMONTH )
!
! THIS IS A NEW VERSION OF THIS SUBROUTINE WHICH COMBINES BR CONCENTRATIONS
! FROM MULTIPLE DATA SOURCES
!
!******************************************************************************
!  Subroutine GET_GLOBAL_BR reads global BR from binary punch files stored
!  in the /data/ctm/GEOS_MEAN directory.  This BR data is needed as oxidant
!  for mercury chemistry  
!  (bmy, cdh 7/28/00, 10/3/05, 7/5/06)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) GET_GLOBAL_BR assumes that we are reading global BR data that occupies
!        all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) for IDL
!        regridding code which will produce the appropriate BR files.
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D,  TRANSFER_3D_TROP
      USE TROPOPAUSE_MOD,ONLY : GET_TPAUSE_LEVEL
!      USE LOGICAL_MOD,   ONLY : LVARTROP

      IMPLICIT NONE

#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: THISMONTH

      ! Local variables
      INTEGER              :: I, J, L
      REAL*4               :: ARRAY(IGLOB,JGLOB,LGLOB)
      REAL*4               :: ARRAY2(IGLOB,JGLOB,LLTROP)
      REAL*8               :: XTAU
!      REAL*4               :: TPAUSE(IGLOB,JGLOB)
      CHARACTER(LEN=255)   :: FILENAME
      INTEGER              :: TPL

      ! Location of archived Br
      CHARACTER(LEN=30), PARAMETER :: BR_DIR = 
     &                                '/home/cdh/GC/Archived-Br/'

      ! First time flag
      LOGICAL, SAVE        :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_BR_NEW begins here!
      !=================================================================

      ! Allocate BR array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_BR
         FIRST = .FALSE.
      ENDIF

      ! Get the TAU0 value for the start of the given month
      XTAU = GET_TAU0( THISMONTH, 1, 1985 ) !eds;  cdh-is this OK?

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
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable BR
      CALL TRANSFER_3D( ARRAY, BR_TROP )

      !-----------------------------------------------------------------
      ! Read BrO from pTOMCAT biogenic bromocarbons
      !-----------------------------------------------------------------

      ! Read BRO data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 7,     
     &                 XTAU,      IGLOB,     JGLOB,      
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
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable BR
      CALL TRANSFER_3D( ARRAY, BR_STRAT )

      !-----------------------------------------------------------------
      ! Read BrO from GMI for stratosphere
      !-----------------------------------------------------------------

      ! Read BR data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 7,     
     &                 XTAU,      IGLOB,     JGLOB,      
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
         TPL = GET_TPAUSE_LEVEL(I,J)

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
      ! &        XTAU,      IGLOB,     JGLOB,      
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
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LLTROP,    ARRAY2,    QUIET=.TRUE. )

      ! Assign data from ARRAY2 to the module variable BR
      CALL TRANSFER_3D_TROP( ARRAY2, J_BrO )



      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_BR_NEW

!------------------------------------------------------------------------------

      SUBROUTINE GET_GLOBAL_BR( THISMONTH )
!
!******************************************************************************
!  Subroutine GET_GLOBAL_BR reads global BR from binary punch files stored
!  in the /data/ctm/GEOS_MEAN directory.  This BR data is needed as oxidant
!  for mercury chemistry  
!  (bmy, cdh 7/28/00, 10/3/05, 7/5/06)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) GET_GLOBAL_BR assumes that we are reading global BR data that occupies
!        all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) for IDL
!        regridding code which will produce the appropriate BR files.
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      IMPLICIT NONE

#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: THISMONTH

      ! Local variables
      INTEGER              :: I, J, L
      REAL*4               :: ARRAY(IGLOB,JGLOB,LGLOB)
      REAL*8               :: XTAU
      CHARACTER(LEN=255)   :: FILENAME

      ! Location of archived Br
      CHARACTER(LEN=30), PARAMETER :: BR_DIR = 
     &                                '/home/cdh/GC/Archived-Br/'

      ! First time flag
      LOGICAL, SAVE        :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_BR begins here!
      !=================================================================

      ! Allocate BR array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_BR
         FIRST = .FALSE.
      ENDIF

      ! Filename
      FILENAME = TRIM( BR_DIR ) // 'Br_3Dglobal.' 
     &           // GET_NAME_EXT() // '.' // GET_RES_EXT()

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_BR: Reading Br, BrO from: ', a )

      ! Time stamp for data
      XTAU = GET_TAU0( THISMONTH, 1, 1985 ) 

      ! Read BR data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 6,     
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable BR
      CALL TRANSFER_3D( ARRAY, BR_MERGE )


      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_BR

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GLOBAL_BR
!
!******************************************************************************
!  Subroutine INIT_GLOBAL_BR allocates and zeroes the BR array, which holds 
!  global monthly mean BR concentrations. (bmy, cdh, 7/28/00, 5/4/04, 7/6/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE" 

      ! Local variables
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

      ! Return to calling program
      END SUBROUTINE INIT_GLOBAL_BR
      
!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GLOBAL_BR
!
!******************************************************************************
!  Subroutine CLEANUP_GLOBAL_BR deallocates the BR array. 
!  (bmy, cdh, 7/28/00, 7/6/06)
!
!  NOTES:
!******************************************************************************
!        
      !=================================================================
      ! CLEANUP_GLOBAL_BR begins here!
      !=================================================================
      IF ( ALLOCATED( BR_TROP  ) ) DEALLOCATE( BR_TROP  ) 
      IF ( ALLOCATED( BR_STRAT ) ) DEALLOCATE( BR_STRAT ) 
      IF ( ALLOCATED( BR_MERGE ) ) DEALLOCATE( BR_MERGE ) 

      IF ( ALLOCATED( J_BrO    ) ) DEALLOCATE( J_BrO  ) 

      IF ( ALLOCATED( BrO_TROP  ) ) DEALLOCATE( BrO_TROP  ) 
      IF ( ALLOCATED( BrO_STRAT ) ) DEALLOCATE( BrO_STRAT )      
      IF ( ALLOCATED( BrO_MERGE ) ) DEALLOCATE( BrO_MERGE ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_GLOBAL_BR

!------------------------------------------------------------------------------

      END MODULE GLOBAL_BR_MOD
