! $Id: emf_scale.f,v 1.1 2003/06/30 20:26:09 bmy Exp $
      SUBROUTINE EMF_SCALE( I,    J,    N,     NN, 
     &                      IREF, JREF, JSCEN, XEMISR, XEMISRN )       
!
!*****************************************************************************
!  Subroutine EMF_SCALE (bmy, 4/2/98, 11/6/02) does the following:
! 
!  (1) Saves original values of EMISR, EMISRN, EMISPN
!      so that they can be restored later (after scaling)
!
!  (2) If LFFNOX=F then set anthro emissions of NOx and Ox = 0
!
!  (3) Scales emissions to weekend or weekday usage (using scale factors
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
!*****************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD
      
      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_O3"
!-------------------------------
!#     include "comtrid.h"
!-------------------------------
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
!  If LFFNOX=F then set fossil fuel NOx, Ox emissions = 0
!*****************************************************************************
!
      IF ( .not. LFFNOX ) THEN
         IF ( N == IDENOX ) THEN
            EMISRN(IREF,JREF,1:NOXLEVELS) = 0d0
         ELSE IF ( N == IDTOX ) THEN
            EMISR(IREF,JREF,N) = 0d0
         ENDIF
      ENDIF
!
!*****************************************************************************
!  Scale emissions by weekend/weekday:
!  Saturday: JSCEN = 1; Sunday JSCEN=2; Weekday: JSCEN = 3
!*****************************************************************************
!
      ! NOx weekday/weekend emissions: Use SCNR89(1,JSCEN) as scale factor 
      IF ( N == IDENOX ) THEN   
         SFAC89 = SCNR89(1,JSCEN) 

         EMISRN(IREF,JREF,1:NOXLEVELS) = 
     &        EMISRN(IREF,JREF,1:NOXLEVELS) * SFAC89

      ! Ox weekday/weekend emissions: Use SCNR89(1,JSCEN) as scale factor
      ! CO weekday/weekend emissions: Use SCNR89(2,JSCEN) as scale factor
      ! HC weekday/weekend emissions: Use SCNR89(3,JSCEN) as scale factor
      ! Otherwise:                    Use 1d0             as scale factor
      ELSE 
         IF ( N == IDEOX ) THEN 
            SFAC89 = SCNR89(1,JSCEN)

         ELSE IF ( N == IDECO ) THEN
            SFAC89 = SCNR89(2,JSCEN)

         ELSE IF ( N == IDEALK4 .or. N == IDEISOP .or. 
     &             N == IDPRPE  .or. N == IDEC3H8 ) THEN
            SFAC89 = SCNR89(3,JSCEN) 

         ELSE
            SFAC89 = 1d0

         ENDIF

            EMISR(IREF,JREF,N) = EMISR(IREF,JREF,N) * SFAC89
         ENDIF

      ! Return to calling program
      END SUBROUTINE EMF_SCALE


      



