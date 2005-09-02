! $Id: chem.f,v 1.5 2005/09/02 15:16:59 bmy Exp $
      SUBROUTINE CHEM( FIRSTCHEM, NPTS,     SUNCOS,  SUNCOSB, 
     &                 CLOUDS,    ALT,      SURFALT, TOTO3, 
     &                 IDXAIR,    IDXO3,    OPTD,    UVALBEDO  )
!
!******************************************************************************
!  Subroutine CHEM is the driver for the photolysis code (either the ancient
!  SLOW-J code or the newer FAST-J code) and the SMVGEAR chemistry solver. 
!  (lwh, jyl, gmg, djj, 1990's; bmy, 4/1/03, 6/23/05)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) FIRSTCHEM (LOGICAL) : FIRSTCHEM=T if we are on the first chem timestep
!  (2 ) NPTS      (INTEGER) : # of grid boxes passed to SMVGEAR solver
!  (3 ) SUNCOS    (REAL*8 ) : Array for COS( Solar Zenith Angle ) [unitless]
!  (4 ) SUNCOSB   (REAL*8 ) : Array for SUNCOS at time 1 hr from now [unitless]
!  (5 ) CLOUDS    (REAL*8 ) : Albedos at 2-km intervals from 0 to 20-km
!  (6 ) ALT       (REAL*8 ) : 3-D Grid box midpoint altitudes [cm]
!  (7 ) SURFALT   (REAL*8 ) : Surface grid box altitudes [m]
!  (8 ) TOTO3     (REAL*8 ) : Total ozone column [molec/cm3]
!  (9 ) IDXAIR    (INTEGER) : Index for standard temperature profile
!  (10) IDXO3     (INTEGER) : Index for standard ozone profile
!  (11) OPTD      (REAL*8 ) : Grid box optical depths [unitless]
!  (12) UVALBEDO  (REAL*8 ) : UV albedoes at the surface [unitless]
!
!  ALT, SURFALT, TOTO3, IDXAIR, IDXO3 are needed for SLOW-J photolysis
!  OPTD, UVALBEDO                     are needed for FAST-J photolysis
! 
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!  (2 ) Remove LSAMERAD argument, it's obsolete.  Add SUNCOSB to the arg list.
!        Now remove test for GETIFSUN, since FAST-J has an internal test for
!        daytime/nighttime.  Pass SUNCOSB to "physproc.f". (gcc, bmy, 7/30/03)
!  (3 ) SLOW-J is now obsolete; remove LSLOWJ #ifdef blocks (bmy, 6/23/05)
!******************************************************************************
!
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "comode.h"     ! SMVGEAR II arrays

      ! Arguments
      LOGICAL, INTENT(INOUT) :: FIRSTCHEM
      INTEGER, INTENT(IN)    :: NPTS
      INTEGER, INTENT(IN)    :: IDXAIR(NLAT)
      INTEGER, INTENT(IN)    :: IDXO3(NLAT)
      REAL*8,  INTENT(IN)    :: SUNCOS(MAXIJ)
      REAL*8,  INTENT(IN)    :: SUNCOSB(MAXIJ)
      REAL*8,  INTENT(IN)    :: CLOUDS(MAXIJ,11)
      REAL*8,  INTENT(IN)    :: ALT(MAXIJ,NPVERT)
      REAL*8,  INTENT(IN)    :: SURFALT(MAXIJ)
      REAL*8,  INTENT(IN)    :: TOTO3(NLAT)
      REAL*8,  INTENT(IN)    :: UVALBEDO(IIPAR,JJPAR)
      REAL*8,  INTENT(IN)    :: OPTD(LLPAR,IIPAR,JJPAR)                   

      !=================================================================
      ! CHEM begins here!
      !=================================================================
      
      ! At present, we are only doing tropospheric chemistry, which 
      ! for the moment we are storing in SMVGEAR II's "urban" slot
      NCS = NCSURBAN

      !=================================================================
      ! Call photolysis routine to compute J-Values
      !=================================================================

      ! FAST-J only: call FAST_J to compute J-values
      CALL FAST_J( SUNCOS, OPTD, UVALBEDO )              

      !================================================================
      ! Call chemistry routines
      !================================================================

      ! PHYSPROC calls both CALCRATE, which computes rxn rates 
      ! and SMVGEAR, which is the chemistry solver
      CALL PHYSPROC( SUNCOS, SUNCOSB )
      
      ! SCHEM applies a simplified strat chemistry in order
      ! to prevent stuff from building up in the stratosphere
      CALL SCHEM

      ! We have gone thru one chemistry timestep already
      FIRSTCHEM = .FALSE.

      ! Return to calling program
      END SUBROUTINE CHEM
