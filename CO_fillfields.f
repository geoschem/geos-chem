C $Id: CO_fillfields.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
       SUBROUTINE CO_FILLFIELDS( LMN, NCLIMATOLOGY )
!
!*****************************************************************************
! Created by bnd/bey (12/98). -- updated (bmy, 7/20/00)
!  SR CO_FILLFIELDS unzips mean punch files the directory
!  /data/ctm/GEOS_MEAN, and then calls subroutine READFIELDS to
!  read the NOx, HC, and OH fields stored in these punch files.
!
! NOTES:
! (1) Now use function GET_TAU0 (from "bpch2_mod.f") to return the TAU0 
!     value used to index the binary punch file. (bmy, 7/20/00)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : GET_TAU0

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT
#     include "CMN_SETUP" ! TEMP_DIR, DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: LMN, NCLIMATOLOGY 

      ! Local variables
      REAL*8              :: TAU_MONTH

      CHARACTER(LEN=255)  :: TEMPO,TEMPO2
      CHARACTER(LEN=255)  :: CHARCP,CHARRM
      CHARACTER(LEN=255)  :: FIELD_DIR

      ! Month names for 1994
!      CHARACTER(LEN=5)    :: BMONTH(12) =                        (/
!     &                          'jan94', 'feb94', 'mar94', 'apr94',
!     &                          'may94', 'jun94', 'jul94', 'aug94',
!     &                          'sep94', 'oct94', 'nov94', 'dec94' /)
!
!*****************************************************************************
!  CO_FILLFIELDS begins here!
!
!  Unzip and copy to the punch file to the temporary directory.
!*****************************************************************************
!
      ! Read the OH file from the data directory (bmy, 5/17/00)
cbnd      FIELD_DIR = TRIM( DATA_DIR )
cbnd       FIELD_DIR = '/users/ctm/bnd/rundao/'
       FIELD_DIR = '/r/amalthea/He/data/ctm/GEOS_MEAN/1994/'
cbnd
      ! temporary file name
      TEMPO     = 'tempo'
      TEMPO2     = 'rm -f '//TRIM( TEMP_DIR )//'tempo'
      CALL SYSTEM( TRIM ( TEMPO2 ) )

      ! String argument for the SYSTEM call
cbnd      CHARCP = TRIM( UNZIP_CMD  ) // SPACE(1:1)           //
cbnd     &         TRIM( FIELD_DIR  ) // TRIM( BMONTH( LMN ) ) //
cbnd     &         TRIM( ZIP_SUFFIX ) // TRIM( REDIRECT     ) //

      IF ( LMN.NE.12 ) THEN

      CHARCP = TRIM( 'cat'  ) // SPACE(1:1)           //
     &      TRIM( FIELD_DIR ) // TRIM( 'ctm.bpch.1994' ) //
     &      TRIM( REDIRECT  ) //
     &      SPACE(1:1)        // TRIM( TEMP_DIR     ) //
     &      TRIM( TEMPO     )

      ELSE

      CHARCP = TRIM( 'cat'  ) // SPACE(1:1)           //
     &      TRIM( FIELD_DIR ) // TRIM( 'ctm.bpch.1994dec' ) //
     &      TRIM( REDIRECT  ) //
     &      SPACE(1:1)        // TRIM( TEMP_DIR     ) //
     &      TRIM( TEMPO     )

      ENDIF

      WRITE( 6, '(''=========================================='')' )
      WRITE( 6, '(''Inside FILLFIELDS''                         )' )
      WRITE( 6, '(a)') TRIM( CHARCP )

      CALL SYSTEM( TRIM( CHARCP ) ) 
!
!*****************************************************************************
!  Call CO_READFIELDS to open the file and read the NOx, HC, and OH fields
!  stored in the punch file that we just copied.
!
!  OH values will be returned in array BIJ.
!*****************************************************************************
!
      ! Select the proper month (bmy, 4/14/00)
      TAU_MONTH = GET_TAU0( LMN, 1994 )
       print*,'TAU_MONTH=',TAU_MONTH

      CALL CO_READFIELDS( TEMPO,NCLIMATOLOGY,TAU_MONTH )
!
!*****************************************************************************
!  Remove the temporary file, after we have read it
!*****************************************************************************
!
      CHARRM = TRIM( REMOVE_CMD ) // SPACE(1:1)   //
     &         TRIM( TEMP_DIR   ) // TRIM( TEMPO )

      CALL SYSTEM( TRIM( CHARRM ) )
      CALL SYSTEM( TRIM ( TEMPO2 ) )
!
!*****************************************************************************
!  Return to calling program
!*****************************************************************************
!
      END SUBROUTINE CO_FILLFIELDS
