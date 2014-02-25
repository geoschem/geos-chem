!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: jvalue 
!
! !DESCRIPTION: subroutine JVALUE calculates the actinic flux at each level 
!  for the current solar zenith angle.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE JVALUE( SA, am_I_Root )
!
! !USES:
!
      USE CMN_FJ_MOD,   ONLY : JPMAX, JPPJ, JPNL
      USE CMN_SIZE_MOD, ONLY : LLPAR
      USE JV_CMN_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: SA          ! Surface albedo [unitless]
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REMARKS:
!-----------------------------------------------------------------------
!        quit when SZA > 98.0 deg ==> tangent height = 63 km
!             or         99.                           80 km
!-----------------------------------------------------------------------
!  AVGF   Attenuation of beam at each level for each wavelength
!  FFF    Actinic flux at each desired level
!  WAVE   Effective wavelength of each wavelength bin
!  XQO2   Absorption cross-section of O2
!  XQO3   Absorption cross-section of O3
! 
! !REVISION HISTORY: 
!         1997 - O. Wild     - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  31 Jul 2012 - R. Yantosca - Added ProTeX headers
!  10 Aug 2012 - R. Yantosca - Replace LPAR with LLPAR
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer j, k
      real*8  wave, xseco3, xseco2
      real*8  AVGF(LLPAR),XQO3(NB),XQO2(NB)
!
      do J=1,jpnl
        do K=NW1,NW2
          FFF(K,J) = 0.d0
        enddo
      enddo
c
c---SZA check
!      IF(  am_I_Root ) THEN
!         write(6,1000) SZA, RFLECT, (OD(nslon,nslat,j),j=1,LLPAR)
!      ENDIF
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
      END SUBROUTINE JVALUE
!EOC
