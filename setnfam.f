! $Id: setnfam.f,v 1.2 2003/07/08 15:31:06 bmy Exp $
      SUBROUTINE SETNFAM 
!
!*****************************************************************************
!  Subroutine SETNFAM reads in number of families for production and loss 
!  info.  This number is needed when setting up diagnostic arrays. 
!  (bey, bmy, 3/16/00, 4/21/03)
!
!  SETNFAM is called from NDXX_SETUP, when ND65 is turned on. 
!
!  NOTES:
!  (1 ) Original code is from Loretta Mickley and Isabelle Bey
!  (2 ) Trap I/O errors with the ERR statement 
!        Also made cosmetic changes, updated comments. (bmy, 3/16/00)
!  (3 ) Now reference ERROR_STOP from "error_mod.f".  Updated comments and
!        made some cosmetic changes. (bmy, 10/15/02)
!  (4 ) Replace 65 with IU_FILE from "file_mod.f" (bmy, 4/21/03)
!*****************************************************************************
!
      ! Reference sto F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IU_FILE

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN"
      
      ! Local variables
      CHARACTER(LEN=7) :: JUNK
      INTEGER          :: ICOUNT

      !=================================================================
      ! SETNFAM begins here! 
      !
      ! Read the "prodloss.dat" file and keep track of the number of 
      ! families for which to archive production/loss diagnostics.
      !=================================================================

      ! Initialize counter
      ICOUNT = 0

      ! Open the "prodloss.dat" file 
      OPEN( IU_FILE, FILE='prodloss.dat', STATUS='OLD',
     &               FORM='FORMATTED',    ERR=700 )

      ! Read the first 7 characters from prodloss.dat
 99   READ( IU_FILE, '(a7)', END=999, ERR=800 ) JUNK

      ! If the string is "*family", then increment ICOUNT
      IF ( JUNK == '*family' ) THEN
         ICOUNT = ICOUNT + 1
      ENDIF

      ! Read next line
      GO TO 99

 999  CONTINUE

      ! Save the number of families in NFAM and echo to std output
      NFAM = ICOUNT
      
      WRITE( 6, 90 ) NFAM
 90   FORMAT( '     - SETNFAM: # of families for ND65 diag: ', i3 )

      ! No families found!  Print error message.
      IF ( NFAM == 0 ) THEN
         CALL ERROR_STOP( 'No families found!', 'setnfam.f' )
      ENDIF

      ! Close the file
      CLOSE( IU_FILE )

      ! Exit
      RETURN

      ! Trap file open errors
 700  CONTINUE
      CALL ERROR_STOP( 'Open error in "prodloss.dat"!', 'setnfam.f' )

      ! Trap file read errors
 800  CONTINUE
      CALL ERROR_STOP( 'Read error in "prodloss.dat"!', 'setnfam.f' )

      ! Return to calling program
      END SUBROUTINE SETNFAM
