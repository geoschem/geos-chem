! $Id: setjfam.f,v 1.1 2003/06/30 20:26:02 bmy Exp $
      SUBROUTINE SETJFAM( NACTIVE, NINAC )
!
!******************************************************************************
!  Subroutine SETJFAM stores reads the "prodloss.dat" file for the ND65
!  chemical family production/loss diagnostics (ljm, bmy, 1999, 4/21/03)
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) NACTIVE (INTEGER) : Number of active SMVGEAR species
!  (2 ) NINAC   (INTEGER) : Number of inactive SMVGEAR species
!
!  NOTES:
!  (1 ) Replace NAMESPEC with NAMEGAS for SMVGEAR II.  Added comment header
!        and updated comments.  Now references IU_FILE and IOERROR from
!        F90 module "file_mod.f".  Now trap I/O errors using routine IOERROR.
!        Make DEFMR a parameter for safety's sake.   Need to increment NACTIVE
!        for SMVGEAR II or else the last species will be overwritten w/ the 
!        first ND65 family.  Set NCS = NCSURBAN, since we have defined our 
!        GEOS-CHEM mechanism in the urban slot of SMVGEAR II.(bmy, 4/21/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IOERROR, IU_FILE

      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "comode.h"     ! SMVGEAR II arrays

      ! Arguments
      INTEGER, INTENT(INOUT) :: NACTIVE, NINAC
      
      ! Local variables
      INTEGER                :: I, J, JGAS0, JGAS, IOS
      REAL*8, PARAMETER      :: DEFMR = 0d0
      CHARACTER(LEN=7)       :: JUNK

      !=================================================================
      ! SETJFAM begins here!
      !=================================================================
      
      ! Need increment NACTIVE for SMVGEAR II or else the last species
      ! will be overwritten w/ the first ND65 family (bmy, 4/18/03)
      NACTIVE = NACTIVE + 1
      JGAS0   = NACTIVE 

      ! Set NCS = NCSURBAN, since we have defined our GEOS-CHEM 
      ! mechanism in the urban slot of SMVGEAR II. (bmy, 4/21/03)
      NCS     = NCSURBAN

      !=================================================================
      ! Read in family names for prod and loss. Assume these
      ! families are active.  Assume initial mixing ratio defmr. 
      ! Note that when setjfam is called, nactive = active species +1.
      !=================================================================

      ! Open file
      OPEN( IU_FILE,      FILE='prodloss.dat', 
     &      STATUS='OLD', FORM='FORMATTED', IOSTAT=IOS )

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'setjfam:1' )

      DO 

         ! Read each line
         READ( IU_FILE, '(a7)', IOSTAT=IOS ) JUNK

         ! IOS < 0 is EOF; otherwise it's an error
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'setjfam:2' )

         ! If the line says *family, then there is data to read...
         IF ( JUNK == '*family' ) THEN

            ! Read the tracer name
            READ( IU_FILE, '(22x,6a8)', IOSTAT=IOS ) XINP(1)
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'setjfam:3' )

            ! Update variables
            JGAS              = NACTIVE
            NTSPEC(NCS)       = NACTIVE + IGAS - NINAC
            NAMEGAS(JGAS)     = XINP(1)
            QBKCHEM(JGAS,NCS) = DEFMR
            NACTIVE           = NACTIVE + 1
         ENDIF

      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      !=================================================================
      ! Write out family names to "smv2.log" file
      !=================================================================
      WRITE( IO93, '(/,a)'      ) REPEAT( '=', 79 )
      WRITE( IO93, '(a)'        ) 'Families for prod or loss output:'
      WRITE( IO93, '(a,/)'      ) REPEAT( '=', 79 )
      WRITE( IO93, '(10(a7,1x))' ) ( TRIM( NAMEGAS(J) ), J=JGAS0,JGAS )

      ! Return to calling program
      END SUBROUTINE SETJFAM
