!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_bc_mod.f
!
! !DESCRIPTION:  Module GLOBAL_BC_MOD contains variables and routines for reading the
! global monthly mean OC concentration from disk. Based on module GLOBAL_OH_MOD.
! (clf, 1/19/2011).
!\\
!\\
! !INTERFACE:
!
      MODULE GLOBAL_BC_MOD
! 
! !USES:
!
      IMPLICIT NONE

!
! !PUBLIC TYPES:
!

! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: GET_GLOBAL_BC
      PUBLIC :: CLEANUP_GLOBAL_BC
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: BC
!
! !REVISION HISTORY:
!  19 January 2011 - C.L. Friedman - Initial Version
!
! !REMARKS:
! Under construction
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header:
!
!  Module GLOBAL_BC_MOD contains variables and routines for reading the global monthly
!  mean BC concentration from disk. (clf, 1/19/2011)
!
!  Module Variables:
!  ===========================================================================
!  (1 ) BC    (REAL*8)    : Array to store global monthly mean BC field
!
!  Module Routines:
!  ===========================================================================
!  (1 ) GET_BC            : Wrapper for GET_GLOBAL_BC
!  (2 ) GET_GLOBAL_BC     : Reads global monthly mean BC from disk
!  (3 ) INIT_GLOBAL_BC    : Allocates and initializes the BC array
!  (4 ) CLEANUP_GLOBAL_BC : Deallocates the BC array
!
!  GEOS-CHEM modules referenced by global_bc_mod.f
!  ===========================================================================
!  (1 ) bpch2_mod.f    : Module containing routines for binary punch file I/O
!  (2 ) error_mod.f    : Module containing NaN and other error-check routines
!
!  Notes:
!  ============================================================================
!  (1) 19 January 2011 C.L. Friedman - Initial version
!
!******************************************************************************
!
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Array to store global monthly mean BC field
      REAL*8, ALLOCATABLE :: BC(:,:,:)
!      REAL*8, ALLOCATABLE :: BCPHIL(:,:,:)
!      REAL*8, ALLOCATABLE :: BCPHOB(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GET_GLOBAL_BC
!
! !DESCRIPTION: GET_GLOBAL_BC reads global BC from binary punch files stored on disk.
!  BC data is needed for partitioning of gas phase organics onto BC particles (e.g., 
!  POPs). (clf, 1/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_BC( THISMONTH, THISYEAR )
!
! !INPUT PARAMETERS: 
! 
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
!      References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : BC_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D
      USE LOGICAL_MOD,   ONLY : LFUTURE

      USE CMN_SIZE_MOD        ! Size parameters ! H. Amos, 7 Mar 2013

      IMPLICIT NONE
!H. Amos, 7 Mar 2013
!#     include "CMN_SIZE"    ! Size parameters

!
! !REVISION HISTORY: 
!  19 January 2011 - C.L. Friedman - Initial Version
!
! !REMARKS:
! Under construction
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  GET_GLOBAL_BC reads global BC from binary punch files stored on disk.
!  BC data is needed for partitioning of gas phase organics onto BC particles (e.g., 
!  POPs). (clf, 1/19/2011)
!
!  Arguments as Input:
!  ============================================================================
!  (1 )  THISMONTH (INTEGER) : Current month number (1-12)
!
!
!  Local variables:
!  ============================================================================
!  (1 ) I, J, L                      INTEGER
!  (2 ) ARRAY(IGLOB, JGLOB, LGLOB)   REAL*4
!  (3 ) XTAU                         REAL*8
!  (4 ) FILENAME                     CHARACTER


!  NOTES:
!  (1 ) Based on GET_GLOBAL_OH subroutine in GLOBAL_OH_MOD
!******************************************************************************

!BOC

      ! Arguments
      INTEGER, INTENT(IN)  :: THISMONTH, THISYEAR

      ! Local variables
      INTEGER              :: I, J, L
      REAL*4               :: ARRAY(IGLOB,JGLOB,LGLOB)
      REAL*8               :: XTAU
      CHARACTER(LEN=255)   :: FILENAME

      ! First time flag
      LOGICAL, SAVE        :: FIRST = .TRUE. 

 
      !=================================================================
      ! GET_GLOBAL_BC begins here!
      !=================================================================

      ! Allocate BC array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_BC
         FIRST = .FALSE.
      ENDIF

! Default: GEOS-5 @ 4x5 reduced vertical layers (H. Amos, 15Mar2013)
      FILENAME = '/home/hamos/GC/input_files/POPs/PAHs_2004_GEOS' //
     &           '/4x5/CLF_fullchem/BCPO_FC_4x5.bpch'

      ! If GCAP:
#if defined( GCAP)

!      IF ( LFUTURE ) THEN

!      WRITE(6,*) 'LFUTURE=', LFUTURE

         ! Friedman GCAP OC/BC aerosol simulations w/future anthropogenic emissions:
         FILENAME = '/net/fs03/d0/geosdata/data/GCAP_4x5/PAHs_2004/' //
     &              'GCAP_OCBC/PCPE_BCPO.bpch'

!      ELSE

         ! Friedman full-chem OC/BC aerosol simulations:
!         FILENAME = '/net/fs03/d0/geosdata/data/GCAP_4x5/PAHs_2004/' //
!     &              'GCAP_OCBC/presem_BCPO_GCAP.bpch'

!      ENDIF


#elif defined ( GRIDREDUCED ) && defined( GEOS_5 )
      ! My full-chem OC/BC aerosol simulations:
!      LFUTURE = .FALSE.

      FILENAME = '/net/fs03/d0/geosdata/data/GEOS_4x5/PAHs_2004/' //
     &           '4x5/CLF_fullchem/BCPO_FC_4x5.bpch'

!      If we want to read in from the input file:
!      FILENAME = TRIM( BC_DIR ) // 'BC_3Dglobal.' // GET_NAME_EXT() // 
!     &                              '.'           // GET_RES_EXT()

#endif

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_BC: Reading BC from: ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      XTAU = GET_TAU0( THISMONTH, 1, THISYEAR )

      ! From Qiaoqiao's or CLF's run:
      CALL READ_BPCH2( FILENAME, 'IJ-24H-$', 13,     
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LGLOB,     ARRAY,     QUIET=.FALSE. )

      ! Assign data from ARRAY to the module variable BCPHOB
      CALL TRANSFER_3D( ARRAY, BC )

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_BC
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  INIT_GLOBAL_BC
!
! !DESCRIPTION: Subroutine INIT_GLOBAL_BC allocates and zeroes the BC array, 
! which holds global monthly mean BC concentrations. (clf, 1/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_BC
!
! !INPUT PARAMETERS: 
!
! 
! !INPUT/OUTPUT PARAMETERS:
!
!
! !OUTPUT PARAMETERS:

      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE CMN_SIZE_MOD        ! Size parameters ! H. Amos, 7 Mar 2013

!#     include "CMN_SIZE" 
!
!
! !REVISION HISTORY: 
!  19 January 2011 - C.L. Friedman - Initial Version
!
! !REMARKS:
! Under construction
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine INIT_GLOBAL_BC deallocates and zeroes the BC array, which holdes
!  global monthly mean BC concentrations. (clf, 1/19/2011)
!
!  NOTES:
!  (1 ) Based on INIT_GLOBAL_OH subroutine in GLOBAL_OH_MOD
!******************************************************************************

!BOC

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_GLOBAL_BC begins here!
      !=================================================================

      ! Allocate BC array
      ALLOCATE( BC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )

!      ALLOCATE( BCPHIL( IIPAR, JJPAR, LLPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPHIL' )

!      ALLOCATE( BCHPOB( IIPAR, JJPAR, LLPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPHOB' )

      ! Zero BC array
      BC = 0d0
!      BCPHIL = 0d0
!      BCPHOB = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_GLOBAL_BC    

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CLEANUP_GLOBAL_BC
!
! !DESCRIPTION: Subroutine CLEANUP_GLOBAL_BC deallocates the BC array. (clf, 1/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_BC
!
! !INPUT PARAMETERS: 
!
! 
! !INPUT/OUTPUT PARAMETERS:
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  19 January 2011 - C.L. Friedman - Initial Version
!
! !REMARKS:
! Under construction
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine CLEANUP_GLOBAL_BC deallocates the BC array. (clf, 1/19/2011)
!
!  NOTES:
!  (1 ) Based on CLEANUP_GLOBAL_OH subroutine in GLOBAL_OH_MOD
!******************************************************************************

!BOC

      !=================================================================
      ! CLEANUP_GLOBAL_BC begins here!
      !=================================================================
      IF ( ALLOCATED( BC ) ) DEALLOCATE( BC ) 
!      IF ( ALLOCATED( BCPHIL ) ) DEALLOCATE( BCPHIL ) 
!      IF ( ALLOCATED( BCPHOB ) ) DEALLOCATE( BCPHOB )
     
      ! Return to calling program
      END SUBROUTINE CLEANUP_GLOBAL_BC

!EOC
!------------------------------------------------------------------------------
      
      END MODULE GLOBAL_BC_MOD
