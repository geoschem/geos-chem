C $Id: funbas.f,v 1.1 2003/06/30 20:26:02 bmy Exp $
      REAL*8 FUNCTION FUNBAS(I,ISPEC,TT)

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "comsol.h"

C**********************************************************************
C                                                                     *
C  HARVARD TROPOSPHERIC CHEMISTRY MODULE FOR 3-D APPLICATIONS         *
C  by Larry Horowitz, Jinyou Liang, Gerry Gardner, Prof. Daniel Jacob *
C                           (Release V2.0)                            *
C**********************************************************************
C
C***************************************************************************C
C   Calculates temperature dependence of OZONE cross-sections from
C   2675-3425 based on RJS's analysis of Bass data. Formulation is:
C       Q(L) = ABASS(L) + BBASS(L)*T(DEG.C) + CBASS(L)*T(DEG.C)**2
C     WAVELENGTHS/2675,2725,2775,2825,2875,2925,2975,3O25,
C                 3075,3125,3175,3225,3275,3325,3375,3425/
C***************************************************************************C
C
      INTEGER I,ISPEC,IB
      REAL*8 TT,TB
      REAL*8 ABASS(16),BBASS(16),CBASS(16)
C
      DATA ABASS/ 8.6348E-18, 7.9559E-18, 4.7925E-18, 3.0765E-18,
     2            1.8787E-18, 1.0338E-18, 5.3570E-19, 2.7430E-19,
     3            1.3900E-19, 7.0600E-20, 3.5700E-20, 1.7423E-20,
     4            8.3145E-21, 4.0000E-21, 1.9576E-21, 8.8154E-22/
      DATA BBASS/ 2.9036E-22, 1.4341E-21, 1.1244E-21, 1.1463E-21,
     2            1.1953E-21, 9.3861E-22, 7.5860E-22, 4.8831E-22,
     3            3.2141E-22, 1.9707E-22, 1.1714E-22, 7.1821E-23,
     4            4.1105E-23, 2.4335E-23, 1.5765E-23, 1.0052E-23/
      DATA CBASS/-3.7953E-24, 1.1119E-23,-2.5398E-24,-3.9890E-24,
     2            2.7978E-24, 2.2221E-24, 3.6727E-24, 1.9728E-24,
     3            1.4252E-24, 1.0259E-24, 6.4958E-25, 3.7252E-25,
     4            2.0089E-25, 1.3979E-25, 1.0669E-25, 6.1010E-26/
C
C
      IF (I.LT.23 .OR. I.GT.38) THEN
         FUNBAS = XSECT(I,ISPEC,1,1)
      ELSE
         TB = TT - 273.15
         IB = I-22
         FUNBAS = ABASS(IB) + BBASS(IB)*TB + CBASS(IB)*TB*TB
      ENDIF
      RETURN
      END
