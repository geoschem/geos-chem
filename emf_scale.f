! $Id: emf_scale.f,v 1.5 2009/05/06 14:14:46 ccarouge Exp $
      SUBROUTINE EMF_SCALE( I,    J,    N,     NN, 
     &                      IREF, JREF, JSCEN, XEMISR, XEMISRN )       
!
!******************************************************************************
!  Subroutine EMF_SCALE (bmy, 4/2/98, 10/3/05) does the following:
! 
!  (1) Saves original values of EMISR, EMISRN, EMISPN
!      so that they can be restored later (after scaling)
!
!  (2) Scales emissions to weekend or weekday usage (using scale factors
!      stored in the SCNR89 array)
!
!  NOTES:
!  (1 ) Use F90 syntax for declarations, etc. (bmy, 4/14/99)
!  (2 ) Now test with N instead of NN.  N is the emission species, and can 
!        be equal to zero, which denotes that the species is not emitted.
!        This is necessary now, since IDEOX always = 0, but IDTOX is always
!        nonzero. (bmy, 4/19/99)
!  (3 ) Commented out special cases via ICASE.  Also made a few cosmetic
!        changes and updated comments. (bmy, 1/2/01)
!  (4 ) Remove old obsolete commented-out code (bmy, 4/20/01)
!  (5 ) Now references "tracerid_mod.f" (bmy, 11/6/02)
!  (6 ) Now references LFFNOX from "logical_mod.f" (bmy, 7/20/04)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Modified to add weekday/weekend scaling to aromatics,
!        C2H4, C2H2 (tmf, 1/7/09)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : IDTALK4, IDTC3H8, IDTISOP, IDTCO
      USE TRACERID_MOD, ONLY : IDTNOX,  IDTOX,   IDTPRPE
      USE TRACERID_MOD, ONLY : IDTMEK,  IDTC2H2, IDTC2H4, IDTACET
      USE TRACERID_MOD, ONLY : IDTBENZ, IDTTOLU, IDTXYLE, IDTC2H6

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN_O3"
#     include "comode.h"
      
      ! Arguments
      INTEGER, INTENT(IN)    :: I, J, N, NN, IREF, JREF, JSCEN
      REAL*8,  INTENT(INOUT) :: XEMISR, XEMISRN(NOXLEVELS)

      ! Local Variables
      INTEGER                :: LL
      REAL*8                 :: SFAC89, PCRE
!
!*****************************************************************************
!  EMF_SCALE begins here!
!
!  Define PCRE, PCUE, PCPE scale factors
!*****************************************************************************
!
      PCRE = .64d0
!
!*****************************************************************************
!  Save original values in temp variables so that they can be restored later
!  Use the appropriate multi-level arrays for NOx emissions
!*****************************************************************************
!
      IF ( NN == IDTNOX ) THEN
         XEMISRN(1:NOXLEVELS) = EMISRN(IREF,JREF,1:NOXLEVELS)
      ELSE
         XEMISR = EMISR(IREF,JREF,N)
      ENDIF
!
!*****************************************************************************
!  Scale emissions by weekend/weekday:
!  Saturday: JSCEN = 1; Sunday JSCEN=2; Weekday: JSCEN = 3
!*****************************************************************************
!
      ! NOx weekday/weekend emissions: Use SCNR89(1,JSCEN) as scale factor 
      IF ( NN == IDTNOX ) THEN   
         SFAC89 = SCNR89(1,JSCEN) 

         EMISRN(IREF,JREF,1:NOXLEVELS) = 
     &        EMISRN(IREF,JREF,1:NOXLEVELS) * SFAC89

      ! Ox weekday/weekend emissions: Use SCNR89(1,JSCEN) as scale factor
      ! CO weekday/weekend emissions: Use SCNR89(2,JSCEN) as scale factor
      ! HC weekday/weekend emissions: Use SCNR89(3,JSCEN) as scale factor
      ! Otherwise:                    Use 1d0             as scale factor
      ELSE 
         IF ( NN == IDTOX ) THEN 
            SFAC89 = SCNR89(1,JSCEN)

         ELSE IF ( NN == IDTCO ) THEN
            SFAC89 = SCNR89(2,JSCEN)

         ELSE IF ( NN == IDTALK4 .or. NN == IDTC2H2 .or.
     &             NN == IDTPRPE .or. NN == IDTC2H4 .or.
     &             NN == IDTC3H8 .or. NN == IDTTOLU .or. 
     &             NN == IDTXYLE ) THEN
            SFAC89 = SCNR89(3,JSCEN) 

         ELSE
            SFAC89 = 1d0

         ENDIF

         EMISR(IREF,JREF,N) = EMISR(IREF,JREF,N) * SFAC89
      ENDIF

      ! Return to calling program
      END SUBROUTINE EMF_SCALE


      



