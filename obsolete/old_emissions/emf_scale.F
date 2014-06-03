!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: emf_scale
!
! !DESCRIPTION: Subroutine EMF\_SCALE does the following:
!
!  \begin{itemize}
!  \item Saves original values of EMISR, EMISRN, EMISPN so that they can be 
!        restored later (after scaling)
!  \item Scales emissions to weekend or weekday usage (using scale factors
!        stored in the SCNR89 array)
!  \end{itemize}
!
! !INTERFACE:
!
      SUBROUTINE EMF_SCALE( I,    J,    N,     NN, 
     &                      IREF, JREF, JSCEN, XEMISR, XEMISRN )       
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTALK4, IDTC3H8, IDTISOP, IDTCO
      USE TRACERID_MOD, ONLY : IDTNO,   IDTO3,   IDTPRPE
      USE TRACERID_MOD, ONLY : IDTMEK,  IDTC2H2, IDTC2H4, IDTACET
      USE TRACERID_MOD, ONLY : IDTBENZ, IDTTOLU, IDTXYLE, IDTC2H6

      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD
      USE CMN_O3_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: I       ! GEOS-Chem longitude index
      INTEGER, INTENT(IN)    :: J       ! GEOS-Chem latitude index
      INTEGER, INTENT(IN)    :: N       ! GEOS-Chem emission species index
      INTEGER, INTENT(IN)    :: NN      ! GEOS-Chem advected tracer index
      INTEGER, INTENT(IN)    :: IREF    ! Offset index I+I0
      INTEGER, INTENT(IN)    :: JREF    ! Offset index J+J0
      INTEGER, INTENT(IN)    :: JSCEN   ! Day index (Sat=1, Sun=2, Weekday=3)

!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,  INTENT(INOUT) :: XEMISR               ! HC emissions, scaled
      REAL*8,  INTENT(INOUT) :: XEMISRN(NOXLEVELS)   ! NOx emissions, scaled
!
! !REMARKS:
!  This is historical baggage...we need to clean this up one of these days. 
! 
! !REVISION HISTORY: 
!  02 Apr 1998 - R. Yantosca - Initial version
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
!  19 Nov 2010 - R. Yantosca - Added ProTeX headers
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox emissions with NO and O3
!                              emissions as part of removal of NOx-Ox
!                              partitioning
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: LL
      REAL*8  :: SFAC89, PCRE

      !=================================================================
      ! EMF_SCALE begins here!
      !
      ! Define PCRE, PCUE, PCPE scale factors
      !=================================================================
      PCRE = .64d0

      !=================================================================
      ! Save original values in temp variables so that they can be 
      ! restored later.   Use the appropriate multi-level arrays for 
      ! NOx emissions
      !=================================================================
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( NN == IDTNO ) THEN
         XEMISRN(1:NOXLEVELS) = EMISRN(IREF,JREF,1:NOXLEVELS)
      ELSE
         XEMISR = EMISR(IREF,JREF,N)
      ENDIF

      !=================================================================
      ! Scale emissions by weekend/weekday:
      !  Saturday: JSCEN = 1; Sunday JSCEN=2; Weekday: JSCEN = 3
      !=================================================================

      ! NOx weekday/weekend emissions: Use SCNR89(1,JSCEN) as scale factor 
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( NN == IDTNO ) THEN   
         SFAC89 = SCNR89(1,JSCEN) 

         EMISRN(IREF,JREF,1:NOXLEVELS) = 
     &        EMISRN(IREF,JREF,1:NOXLEVELS) * SFAC89

      ! Ox weekday/weekend emissions: Use SCNR89(1,JSCEN) as scale factor
      ! CO weekday/weekend emissions: Use SCNR89(2,JSCEN) as scale factor
      ! HC weekday/weekend emissions: Use SCNR89(3,JSCEN) as scale factor
      ! Otherwise:                    Use 1d0             as scale factor
      ELSE 
         ! Replaced Ox with O3 (mpayer, 3/14/13)
         IF ( NN == IDTO3 ) THEN 
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

      END SUBROUTINE EMF_SCALE
!EOC

      



