! $Id: chem.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      SUBROUTINE CHEM( FIRSTCHEM, LSAMERAD, NPTS,    SUNCOS, 
     &                 CLOUDS,    ALT,      SURFALT, TOTO3, 
     &                 IDXAIR,    IDXO3,    OPTD,    UVALBEDO ) 
!
!******************************************************************************
!  Subroutine CHEM is the driver for the photolysis code (either the ancient
!  SLOW-J code or the newer FAST-J code) and the SMVGEAR chemistry solver. 
!  (lwh, jyl, gmg, djj, 1990's; bmy, 4/1/03)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) FIRSTCHEM (LOGICAL) : FIRSTCHEM=T if we are on the first chem timestep
!  (2 ) LSAMERAD  (LOGICAL) : *** Now Obsolete *** 
!  (3 ) NPTS      (INTEGER) : # of grid boxes passed to SMVGEAR solver
!  (4 ) SUNCOS    (REAL*8 ) : Array for COS( Solar Zenith Angle ) [unitless]
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
!******************************************************************************
!
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "comode.h"     ! SMVGEAR II arrays
#if   defined( LSLOWJ )
#     include "comsol.h"     ! SLOW-J arrays
#endif

      ! Arguments
      LOGICAL, INTENT(INOUT) :: FIRSTCHEM
      LOGICAL, INTENT(IN)    :: LSAMERAD
      INTEGER, INTENT(IN)    :: NPTS
      INTEGER, INTENT(IN)    :: IDXAIR(NLAT)
      INTEGER, INTENT(IN)    :: IDXO3(NLAT)
      REAL*8,  INTENT(IN)    :: SUNCOS(MAXIJ)
      REAL*8,  INTENT(IN)    :: CLOUDS(MAXIJ,11)
      REAL*8,  INTENT(IN)    :: ALT(MAXIJ,NPVERT)
      REAL*8,  INTENT(IN)    :: SURFALT(MAXIJ)
      REAL*8,  INTENT(IN)    :: TOTO3(NLAT)
      REAL*8,  INTENT(IN)    :: UVALBEDO(IIPAR,JJPAR)
      REAL*8,  INTENT(IN)    :: OPTD(LLPAR,IIPAR,JJPAR)                   

      ! External functions
      INTEGER, EXTERNAL      :: GETIFSUN

      !=================================================================
      ! CHEM begins here!
      !=================================================================

      ! Do gas phase calculation only at this stage, NCS=1
      NCS   = 1

      ! Test if we have sunlight boxes
      IFSUN = GETIFSUN( SUNCOS )

      !=================================================================
      ! Call photolysis routine to compute J-Values
      !=================================================================
      IF ( ( IFSUN == 1 .AND. .NOT. LSAMERAD ) .OR. FIRSTCHEM ) THEN

#if   defined( LFASTJ )

         ! FAST-J only: call FAST_J to compute J-values
         CALL FAST_J( SUNCOS, OPTD, UVALBEDO )              

#elif defined( LSLOWJ )                                            

         ! SLOW-J only: call SOL to compute column densities and J-values
         CALL SOL( NPTS,  SUNCOS,  ALT,   SURFALT,
     &             TOTO3, CLOUDS, IDXAIR, IDXO3 )

#endif

      ENDIF

#if   defined( LSLOWJ )
      
      ! SLOW-J only: Call FINDXSECT which gets index INAME, which gives the
      ! position of each photolysis species in the cross-section array XSECT  
      IF ( FIRSTCHEM ) CALL FINDXSECT

#endif

      !================================================================
      ! Call chemistry routines
      !================================================================

      ! PHYSPROC calls both CALCRATE, which computes rxn rates 
      ! and SMVGEAR, which is the chemistry solver
      CALL PHYSPROC( SUNCOS )

      ! SCHEM applies a simplified strat chemistry in order
      ! to prevent stuff from building up in the stratosphere
      CALL SCHEM

      ! We have gone thru one chemistry timestep already
      FIRSTCHEM = .FALSE.

      ! Return to calling program
      END SUBROUTINE CHEM
