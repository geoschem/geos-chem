C $Id: MIESCT.f,v 1.1 2009/09/16 14:06:47 bmy Exp $
      SUBROUTINE MIESCT
C-----------------------------------------------------------------------
C   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
C     Prather, 1974, Astrophys. J. 192, 787-792.
C         Sol'n of inhomogeneous Rayleigh scattering atmosphere. 
C         (original Rayleigh w/ polarization)
C     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
C         Raman scattering in the atmospheres of the major planets.
C         (first use of anisotropic code)
C     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
C         Chemistry of a polluted cloudy boundary layer,
C         (documentation of extension to anisotropic scattering)
C
C    takes atmospheric structure and source terms from std J-code
C    ALSO limited to 4 Gauss points, only calculates mean field!
C
C   mean rad. field ONLY (M=1)
C   initialize variables FIXED/UNUSED in this special version:
C   FTOP = 1.0 = astrophysical flux (unit of pi) at SZA, -ZU0, use for scaling
C   FBOT = 0.0 = external isotropic flux on lower boundary 
C   SISOTP = 0.0 = Specific Intensity of isotropic radiation incident from top
C
C   SUBROUTINES:  MIESCT              needs 'jv_mie.cmn'
C                 BLKSLV              needs 'jv_mie.cmn'
C                 GEN (ID)            needs 'jv_mie.cmn'
C                 LEGND0 (X,PL,N)
C                 MATIN4 (A)
C                 GAUSSP (N,XPT,XWT)
C-----------------------------------------------------------------------
      USE JV_MIE_MOD

      IMPLICIT NONE
#     include "define.h"
#     include "define.h"


      integer i, id, im
      real*8  cmeq1
C-----------------------------------------------------------------------
C---fix scattering to 4 Gauss pts = 8-stream
      CALL GAUSSP (N,EMU,WT)
C---solve eqn of R.T. only for first-order M=1
C      ZFLUX = (ZU0*FZ(ND)*ZREFL+FBOT)/(1.0d0+ZREFL)
      ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      M=1
      DO I=1,N
        CALL LEGND0 (EMU(I),PM0,MFIT)
        DO IM=M,MFIT
          PM(I,IM) = PM0(IM)
        ENDDO
      ENDDO
C
      CMEQ1 = 0.25D0
      CALL LEGND0 (-ZU0,PM0,MFIT)
      DO IM=M,MFIT
        PM0(IM) = CMEQ1*PM0(IM)
      ENDDO
C
      CALL BLKSLV
C
      DO ID=1,ND,2
        FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
      ENDDO

      RETURN
      END
