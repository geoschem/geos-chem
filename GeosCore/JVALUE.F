! $Id: JVALUE.f,v 1.1 2009/09/16 14:06:47 bmy Exp $
      SUBROUTINE JVALUE( SA )
!-----------------------------------------------------------------------
c  Calculate the actinic flux at each level for the current SZA value.
!        quit when SZA > 98.0 deg ==> tangent height = 63 km
!             or         99.                           80 km
!-----------------------------------------------------------------------
!  Add the following input variables for CTM interface (bmy, 9/13/99)
!
!  Variable  Type    Dimensn Units   Description
!  --------  ----    ------- -----   -----------
!  SA        dble    -       -       Surface Albedo
!-----------------------------------------------------------------------
c
c     AVGF   Attenuation of beam at each level for each wavelength
c     FFF    Actinic flux at each desired level
c     WAVE   Effective wavelength of each wavelength bin
c     XQO2   Absorption cross-section of O2
c     XQO3   Absorption cross-section of O3
c
!-----------------------------------------------------------------------
#     include "define.h"

      USE CMN_FJ_MOD, ONLY : IPAR, JPAR, LPAR, JPMAX, JPPJ, JPNL
      USE JV_CMN_MOD

      IMPLICIT NONE

!=============== INPUT PARAMETERS ======================================
      REAL*8, INTENT(IN) :: SA

!=============== LOCAL VARIABLES =======================================
      integer j, k
      real*8  wave, xseco3, xseco2
      real*8  AVGF(lpar),XQO3(NB),XQO2(NB)
!
      do J=1,jpnl
        do K=NW1,NW2
          FFF(K,J) = 0.d0
        enddo
      enddo
c
c---SZA check
c      write(6,1000) SZA, RFLECT, (OD(nslon,nslat,j),j=1,lpar)
      if(SZA.gt.szamax) GOTO 99
c
!---Calculate spherical weighting functions
      CALL SPHERE
c
!---Loop over all wavelength bins
      do K=NW1,NW2
        WAVE = WL(K)
        do J=1,NB
          XQO3(J) = XSECO3(K,dble(TJ(J)))
        enddo
        do J=1,NB
          XQO2(J) = XSECO2(K,dble(TJ(J)))
        enddo
!-----------------------------------------
        CALL OPMIE(K,WAVE,XQO2,XQO3,AVGF)
!-----------------------------------------
        do J=1,jpnl
          FFF(K,J) = FFF(K,J) + FL(K)*AVGF(J)
        enddo
      enddo
c
   99 continue
 1000 format('  SZA=',f6.1,' Reflectvty=',f6.3,' OD=',10(1pe10.3))
      return
      end
