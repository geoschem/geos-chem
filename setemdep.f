! $Id: setemdep.f,v 1.1 2003/06/30 20:26:06 bmy Exp $
      SUBROUTINE SETEMDEP( NTRACER )
!
!******************************************************************************
!  Subroutine SETEMDEP stores SMVGEAR reaction numbers (listed in "chem.dat")
!  corresponding to GEOS-CHEM tracers which emit and dry deposit into the
!  NTEMIS and NTDEP index arrays.  (lwh, jyl, gmg, djj, 1994; bmy, 4/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTRACER (INTEGER) : Number of GEOS-CHEM tracers
!
!  NOTES:
!  (1 ) Now references "drydep_mod.f" and "tracerid_mod.f".  Updated comments
!        and made cosmetic changes. (bmy, 12/5/02)
!  (2 ) Cosmetic changes (bmy, 3/14/03)
!  (3 ) Updated for SMVGEAR II (gcc, bdf, bmy, 4/21/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPNAME, NUMDEP
      USE TRACERID_MOD, ONLY : IDEMIS,  IDTRMB, NEMANTHRO, NEMBIOG

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! SMVGEAR II arrays

      ! Arguments
      INTEGER, INTENT(IN) :: NTRACER

      ! Local variables
      INTEGER             :: I, N, NK, NCS_TEMP
      CHARACTER(LEN=14)   :: NAME1

      !=================================================================
      ! SETEMDEP begins here!
      !=================================================================      

      ! Write header to "smv2.log"
      WRITE( IO93, '(/,a)' ) REPEAT( '=', 79 )
      WRITE( IO93, '(a)'   ) 'SETEMDEP: Emission & deposition species'
      WRITE( IO93, '(a,/)' ) REPEAT( '=', 79 )

      !=================================================================
      ! Flag EMISSION REACTIONS in "globchem.dat" for GEOS-CHEM tracers
      !=================================================================

      ! Loop over different kinds of chemistry
      DO NCS = 1, NCSGAS

         ! Loop over GEOS-CHEM tracers
         DO I = 1, NTRACER 

            ! Rxn # for Ith GEOS-CHEM tracer
            NTEMIS(I,NCS) = 0

            ! Loop over emission species from "globchem.dat"
            DO N = 1, NEMIS(NCS)

               ! Rxn # for Nth emission species in "globchem.dat"
               NK = NKEMIS(N,NCS)

               ! Match "chem.dat" rxn number w/ GEOS-CHEM tracer number
               ! IRM is the species # for the first product of the NKth rxn
               ! IDTRMB is the species # of the GEOS-CHEM tracer which emits
               IF ( IDEMIS(I) /= 0 ) THEN
                  IF ( IRM(NPRODLO,NK,NCS) == IDTRMB(I,IDEMIS(I)) ) THEN
                     NTEMIS(I,NCS) = NK
                  ENDIF
               ENDIF
            ENDDO

            ! Write warning to "smv2.log" if no emission rxn was found
            IF ( NTEMIS(I,NCS) == 0 ) THEN
               WRITE( IO93, 100 ) I
 100           FORMAT( 'WARNING: no emission reaction for tracer ', i3 )
            ENDIF
         ENDDO

         ! The total # of emission species will be NEMANTHRO [anthro] + 
         ! NEMBIOG [bio], so reset NEMIS accordingly
         NEMIS(NCS) = NEMANTHRO + NEMBIOG

         ! Echo output to stdout
         WRITE( 6, 110 ) NEMIS(NCS)
 110     FORMAT( '     - SETEMDEP: Number of emitted '
     &           'species in "globchem.dat":', i3 )
      ENDDO                     

      !=================================================================
      ! Flag DRYDEP REACTIONS from "chem.dat" for each GEOS-CHEM tracer
      !=================================================================

      ! There is only drydep in the surface layer, which
      ! is accounted for in the "URBAN" chemistry slot
      NCS = NCSURBAN

      ! Loop over GEOS_CHEM drydep tracers
      DO I = 1, NUMDEP

         ! Rxn # of the Ith GEOS-CHEM drydep tracer
         NTDEP(I) = 0

         ! Loop over drydep species from "globchem.dat"
         DO N = 1, NDRYDEP(NCS)

            ! Rxn number and name of Nth drydep species in "globchem.dat" 
            NK    = NKDRY(N,NCS)
            NAME1 = NAMEGAS(IRM(1,NK,NCS))

            ! If we can match NAME1 against the GEOS-CHEM drydep tracer
            ! names in DEPNAME, then store the rxn number in NTDEP
            IF ( DEPNAME(I) == NAME1 ) THEN
               NTDEP(I) = NK
               EXIT
            ENDIF
         ENDDO

         ! Write warning to "smv2.log" if no drydep rxn was found
         IF ( NTDEP(I) == 0 ) THEN
            WRITE( IO93, 120 ) DEPNAME(I)
 120        FORMAT( 'WARNING: no deposition reaction for ', a )
         ENDIF
      ENDDO

      ! Echo output to stdout
      WRITE( 6, 130 ) NDRYDEP(1)
 130  FORMAT( '     - SETEMDEP: Number of drydep  species '
     &        'in "globchem.dat":', i3 )

      WRITE( 6, 140 ) NUMDEP
 140  FORMAT( '     - SETEMDEP: Number of all GEOS-CHEM '
     &        'drydep species     :', i3 )

      ! Reset NCS = NCSURBAN, since we have defined our GEOS-CHEM
      ! mechanism in the urban slot of SMVGEAR II (bmy, 4/21/03)
      NCS = NCSURBAN

      ! Return to calling program
      END SUBROUTINE SETEMDEP
