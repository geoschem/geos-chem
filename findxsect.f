! $Id: findxsect.f,v 1.1 2003/06/30 20:26:06 bmy Exp $
      SUBROUTINE FINDXSECT
!
!******************************************************************************
!  Subroutine FINDXSECT finds the position of each photolysis species in
!  the SLOW-J cross-section array XSECT (lwh, jyl, gmg, djj, 1994; bmy, 4/1/03)
!
!  NOTES:
!  (1 ) Now references routines ERROR_STOP and GEOS_CHEM_STOP from 
!        "error_mod.f".  NAMESPEC is now NAMEGAS for SMVGEAR II.  Updated 
!         comments. (bdf, bmy, 4/1/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP, GEOS_CHEM_STOP

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! SMVGEAR II arrays
#     include "comsol.h"  ! SLOW-J arrays

      ! Local variables
      INTEGER :: I, NK, K, IFNC, J

      !=================================================================
      ! FINDXSECT begins here!
      !=================================================================

      ! Stop the run safely if NCS /= 1 
      IF ( NCS /= 1 ) THEN
         CALL ERROR_STOP( 'NCS /= 1', 'findxsect.f' )
      ENDIF

      DO 260 I = 1, NPHOT
         NK = NRATES(NCS) + I
         J  = IRM(1,NK,NCS)

         ! Write species name to "smv2.log" file
         WRITE( IO93, * ) NAMEGAS(J), 'PHO# ',I

         INAME(I) = 0

         ! Find photolysis index for this species name
         DO 250 K= 1, MXSPE
            IF ( CNAME(K) == NAMEGAS(J) ) INAME(I) = K
 250     CONTINUE

         ! Test if this species has a cross-section or not
         IF ( INAME(I) == 0 ) THEN
            WRITE(IO93,*) NAMEGAS(J),' has no C.S.'

            ! check IFNC to see if it's OK to have no C.S. for this species
            IFNC = DEFPRAT(NK,NCS)
            IF (IFNC < 24 .OR. IFNC > 26 ) THEN
               WRITE(*,*) 'Error: ',NAMEGAS(J),' has no C.S.'
               WRITE(*,*) '... see 8col.dat'
               STOP 40
            ENDIF

            !for some species, use the following default cross-sections:
            IF ( IFNC == 24 ) WRITE(IO93,*) 'using MP C.S.'
            IF ( IFNC == 25 ) WRITE(IO93,*) 'using MNO3 C.S.'
            IF ( IFNC == 26 ) WRITE(IO93,*) 'using ACET C.S.'
         ENDIF
 260  CONTINUE
      
      ! Return to calling program
      RETURN
      END SUBROUTINE FINDXSECT
