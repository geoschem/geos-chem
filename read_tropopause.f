! $Id: read_tropopause.f,v 1.3 2005/06/30 18:55:30 bmy Exp $
      SUBROUTINE READ_TROPOPAUSE
!
!******************************************************************************
!  Subroutine READ_TROPOPAUSE reads in the annual mean tropopause. 
!  (qli, bmy, 12/13/99, 7/20/04)
!
!  NOTES:
!  (1 ) Call READ_BPCH2 to read in the annual mean tropopause data
!        which is stored in binary punch file format. (bmy, 12/13/99)
!  (2 ) Now also read integer flags for ND27 diagnostic -- these determine
!        how to sum fluxes from boxes adjacent to the annual mean tropoause.
!        (qli, bmy, 1/7/00)
!  (3 ) Cosmetic changes (bmy, 3/17/00)
!  (4 ) Reference F90 module "bpch2_mod" which contains routine "read_bpch2"
!        for reading data from binary punch files (bmy, 6/28/00)
!  (5 ) Call TRANSFER_2D from "transfer_mod.f" to cast data from REAL*4 to
!        INTEGER and also to resize to (IIPAR,JJPAR).  ARRAY needs to be of 
!        size (IGLOB,JGLOB).  Also updated comments and made cosmetic changes. 
!        Removed obsolete variables.(bmy, 9/26/01)
!  (6 ) Removed obsolete code from 9/01 (bmy, 10/26/01)
!  (7 ) Now read annual mean tropopause files from the ann_mean_trop_200202/
!        subdirectory of DATA_DIR (bmy, 1/24/02)
!  (8 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (9 ) Now write file name to stdout (bmy, 4/3/02)
!  (10) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (11) Now call READ_BPCH2 with QUIET=.TRUE. to suppress printing of extra
!        info to stdout.  Also updated FORMAT strings. (bmy, 3/14/03)
!  (12) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPAUSE, IFLX
!---------------------------------------------
! Prior to 7/20/04:
!#     include "CMN_SETUP" ! DATA_DIR
!---------------------------------------------

      ! Local Variables
      INTEGER             :: I, J, LMAX, LMIN, COUNT
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      CHARACTER (LEN=255) :: FILENAME

      !=================================================================
      ! READ_TROPOPAUSE begins here!
      !
      ! Read the annual mean tropopause from disk (binary punch file 
      ! format).  Transfer data into an array of size (IIPAR,JJPAR).
      !=================================================================
      
      ! Create filename
      FILENAME = TRIM( DATA_DIR )                      // 
     &           'ann_mean_trop_200202/ann_mean_trop.' //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()

      ! Write file name to stdout
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - READ_TROPOPAUSE: Reading ', a )

      ! Annual mean tropopause is tracer #1  
      CALL READ_BPCH2( FILENAME, 'TR-PAUSE', 1, 0d0, 
     &                 IGLOB,     JGLOB,     1, ARRAY, QUIET=.TRUE. )

      ! Copy from REAL*4 to INTEGER and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( ARRAY(:,:,1), LPAUSE )

      ! Integer flags for ND27 diagnostic is tracer #4
      CALL READ_BPCH2( FILENAME, 'TR-PAUSE', 4, 0d0, 
     &                 IGLOB,     JGLOB,     1, ARRAY, QUIET=.TRUE. )

      ! Copy from REAL*4 to INTEGER and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( ARRAY(:,:,1), IFLX )

      !=================================================================
      ! L <  LPAUSE(I,J) are tropospheric boxes  --> call SMVGEAR 
      ! L >= LPAUSE(I,J) are stratospheric boxes --> call SCHEM   
      !
      ! LMIN = level where the minimum tropospheric extent occurs
      ! LMAX = level where the maximum tropospheric extent occurs.
      !
      ! Write LMAX and LMIN to the standard output.
      !
      ! Also make sure that LMAX does not exceed LLTROP, since LLTROP 
      ! is used to dimension the chemistry arrays in "comode.h". 
      !=================================================================
      LMIN = MINVAL( LPAUSE ) - 1
      LMAX = MAXVAL( LPAUSE ) - 1
      
      WRITE( 6, 120 ) LMIN 
 120  FORMAT( '     - READ_TROPOPAUSE: Minimum tropospheric extent,',
     &        ' L=1 to L=', i3 )

      WRITE( 6, 130 ) LMAX 
 130  FORMAT( '     - READ_TROPOPAUSE: Maximum tropospheric extent,',
     &        ' L=1 to L=', i3 )
    
      IF ( LMAX > LLTROP ) THEN
         WRITE( 6, '(a)' ) 'READ_TROPOPAUSE: LLTROP is set too low!' 
         WRITE( 6, 131   ) LMAX, LLTROP
 131     FORMAT( 'LMAX = ', i3, '  LLTROP = ', i3 )
         WRITE( 6, '(a)' ) 'STOP in READ_TROPOPAUSE.F!!!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Also write the number of tropopsheric and stratospheric boxes
      ! to the standard output.
      !
      ! Recall that tropospheric boxes extend up to LPAUSE - 1, so 
      ! we have to subtract 1 for each box from the sum of LPAUSE.
      !=================================================================
      COUNT = SUM( LPAUSE ) - ( IIPAR * JJPAR )

      WRITE( 6, 140 ) COUNT
 140  FORMAT( '     - READ_TROPOPAUSE: # of tropopsheric boxes:  ', i8 )
      
      WRITE( 6, 150 ) ( IIPAR * JJPAR * LLPAR ) - COUNT
 150  FORMAT( '     - READ_TROPOPAUSE: # of stratospheric boxes: ', i8 )

      ! Return to calling program
      END SUBROUTINE READ_TROPOPAUSE
