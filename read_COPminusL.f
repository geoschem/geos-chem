C $Id: read_COPminusL.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      SUBROUTINE READ_COPminusL( MCOSTRAT ) 
!
!******************************************************************************
!  Subroutine READ_COPminusL reads production and destruction rates 
!  for CO in the stratosphere. (bnd, bmy, 1999, 9/26/01, 4/2/02)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) MCOSTRAT (INTEGER) : Current month (1-12)
!
!  NOTES:
!  (1 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (2 ) Use TRANSFER_ZONAL from "transfer_mod.f" to copy an array of size 
!        (1,JGLOB,LGLOB) to (1,JJPAR,LLPAR).  Eliminate obsolete variables.
!        Now use 3 arguments in call to GET_TAU0.  Also updated comments
!        and made some cosmetic changes. (bmy, 9/26/01)
!  (3 ) Deleted obsolete, commented-out code from 9/01 (bmy, 11/26/01)
!  (4 ) Now read COprod and COloss files directly from the
!        DATA_DIR/pco_lco_200203/ subdirectory (bmy, 4/2/02)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TRANSFER_MOD, ONLY : TRANSFER_ZONAL

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_SETUP"  ! DATA_DIR
#     include "CMN_CO"     ! BAIRDENS, STT2GCO, CO arrays

      ! Arguments
      INTEGER, INTENT(IN) :: MCOSTRAT

      ! Local variables
      CHARACTER(LEN=255)  :: FILENAME
      REAL*4              :: ARRAY(1,JJPAR,LLPAR) 
      REAL*8              :: XTAU

      print*,'Reading in CO_prod.bpch and CO_loss.bpch in SR chemco.'

      !=================================================================
      ! READ_COPminusL begins here!
      !
      ! Initialize some variables
      !=================================================================
      ARRAY(:,:,:) = 0e0
      CO_prod(:,:) = 0d0
      CO_loss(:,:) = 0d0

      ! TAU value at the beginning of this month
      ! Use "generic" year 1985
      XTAU = GET_TAU0( MCOSTRAT, 1, 1985 )

      !=================================================================
      ! Read in CO production rates
      !=================================================================
!-----------------------------------------------------------------------
! Prior to 4/2/02:
!      FILENAME = TRIM( DATA_DIR ) // 'COprod.'  //
!     &           GET_NAME_EXT()   // '.'         // GET_RES_EXT()
!-----------------------------------------------------------------------
      FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/COprod.' //
     &           GET_NAME_EXT()   // '.'                      // 
     &           GET_RES_EXT()

      ! Read data
      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 9,     XTAU, 
     &                 1,         JGLOB,     LGLOB, ARRAY )

      ! Copy and resize
      CALL TRANSFER_ZONAL( ARRAY(1,:,:), CO_prod )
         
      !=================================================================
      ! Read in CO loss rates
      !=================================================================
!-----------------------------------------------------------------------
! Prior to 4/2/02:
!      FILENAME = TRIM( DATA_DIR ) // 'COloss.'  //
!     &           GET_NAME_EXT()   // '.'        // GET_RES_EXT()
!-----------------------------------------------------------------------
      FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/COloss.'  //
     &           GET_NAME_EXT()   // '.'                       // 
     &           GET_RES_EXT()

      ! Read data
      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 10,    XTAU,    
     &                 1,         JGLOB,     LGLOB, ARRAY )

      ! Copy and resize
      CALL TRANSFER_ZONAL( ARRAY(1,:,:), CO_loss )

      ! Close restart file and return to calling program
      CLOSE( 1 )

      ! Return to calling program
      END SUBROUTINE READ_COPminusL




