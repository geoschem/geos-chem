C $Id: debug_msg.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE DEBUG_MSG( MESSAGE )
!
!******************************************************************************
!  Subroutine DEBUG_MSG prints a message and flushes to stdout.  This is 
!  useful for determining the exact location where errors occur (bmy, 1/7/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MESSAGE (CHAR*(*)) : Descriptive message string
!
!  NOTES: 
!  (1 ) Now just write the message and flush the buffer (bmy, 7/5/01)
!  (2 ) Renamed from "paftop.f" to "debug_msg.f" (bmy, 1/7/02)
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: MESSAGE

      !=================================================================
      ! DEBUG_MSG begins here!
      !
      ! Write label and flush to stdout buffer
      !================================================================= 
      WRITE( 6, '(a)' ) TRIM( MESSAGE )
      CALL FLUSH( 6 )

      ! Return to calling program
      END SUBROUTINE DEBUG_MSG
