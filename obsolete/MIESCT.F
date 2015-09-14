! $Id: MIESCT.f,v 1.1 2009/09/16 14:06:47 bmy Exp $
      SUBROUTINE MIESCT
!-----------------------------------------------------------------------
!   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Sol'n of inhomogeneous Rayleigh scattering atmosphere. 
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field!
!
!   mean rad. field ONLY (M=1)
!   initialize variables FIXED/UNUSED in this special version:
!   FTOP = 1.0 = astrophysical flux (unit of pi) at SZA, -ZU0, use for scaling
!   FBOT = 0.0 = external isotropic flux on lower boundary 
!   SISOTP = 0.0 = Specific Intensity of isotropic radiation incident from top
!
!   SUBROUTINES:  MIESCT              needs 'jv_mie.cmn'
!                 BLKSLV              needs 'jv_mie.cmn'
!                 GEN (ID)            needs 'jv_mie.cmn'
!                 LEGND0 (X,PL,N)
!                 MATIN4 (A)
!                 GAUSSP (N,XPT,XWT)
!-----------------------------------------------------------------------
      USE JV_MIE_MOD

      IMPLICIT NONE


      integer i, id, im
      real*8  cmeq1
!-----------------------------------------------------------------------
!---fix scattering to 4 Gauss pts = 8-stream
      CALL GAUSSP (N,EMU,WT)
!---solve eqn of R.T. only for first-order M=1
!      ZFLUX = (ZU0*FZ(ND)*ZREFL+FBOT)/(1.0d0+ZREFL)
      ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      M=1
      DO I=1,N
        CALL LEGND0 (EMU(I),PM0,MFIT)
        DO IM=M,MFIT
          PM(I,IM) = PM0(IM)
        ENDDO
      ENDDO
!
      CMEQ1 = 0.25D0
      CALL LEGND0 (-ZU0,PM0,MFIT)
      DO IM=M,MFIT
        PM0(IM) = CMEQ1*PM0(IM)
      ENDDO
!
      CALL BLKSLV
!
      DO ID=1,ND,2
        FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
      ENDDO

      RETURN
      END
