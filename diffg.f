C $Id: diffg.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      REAL*8 FUNCTION DIFFG(TK,PRESS,XM)

      IMPLICIT NONE
C===========================================================================
C  This function calculates the molecular diffusivity (m2 s-1) in air for a
C  gas X of molecular weight XM (kg) at temperature TK (K) and 
C  pressure PRESS (Pa).
C
C  We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
C  radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
C  radius of air is given in a Table on p. 479 of Levine [1988].  The Table
C  also gives radii for some other molecules.  Rather than requesting the user
C  to supply a molecular radius we specify here a generic value of 2.E-10 m for
C  all molecules, which is good enough in terms of calculating the diffusivity
C  as long as molecule is not too big.
C
C============================================================================
      REAL*8 TK,PRESS,XM,XMAIR,RADAIR,PI,RADX,RGAS,AVOGAD,AIRDEN,Z,DIAM
      REAL*8 FRPATH,SPEED
      DATA XMAIR/28.8E-3/, RADAIR/1.2E-10/, PI/3.1415926535897932/ 
      DATA RADX/1.5E-10/, RGAS/8.32/, AVOGAD/6.023E23/
C
C* Calculate air density AIRDEN
      AIRDEN = PRESS*AVOGAD/(RGAS*TK)
C
C* Calculate the mean free path for gas X in air: eq. 8.5 of Seinfeld [1986];
C   DIAM is the collision diameter for gas X with air.
      Z = XM/XMAIR
      DIAM = RADX+RADAIR
      FRPATH = 1./(PI*SQRT(1.+Z)*AIRDEN*(DIAM**2.))
C* Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      SPEED = SQRT(8*RGAS*TK/(PI*XM))
C* Calculate diffusion coefficient of gas X in air; eq. 8.9 of Seinfeld
C* [1986]
      DIFFG = (3.*PI/32.)*(1.+Z)*FRPATH*SPEED

      RETURN
      END
