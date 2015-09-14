!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sphere
!
! !DESCRIPTION: Subroutine SPHERE calculates spherical geometry; derives 
!  tangent heights, slant path lengths and air mass factor for each layer. 
!  Not called when SZA > 98 degrees.  Beyond 90 degrees, include treatment 
!  of emergent beam (where tangent height is below altitude J-value desired at).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SPHERE
!
! !USES:
!
      USE CMN_FJ_MOD,   ONLY : JPMAX, JPPJ
      USE CMN_SIZE_MOD, ONLY : LLPAR
      USE JV_CMN_MOD

      IMPLICIT NONE
!
! !REMARKS:
!  GMU     MU, cos(solar zenith angle)
!  RZ      Distance from centre of Earth to each point (cm)
!  RQ      Square of radius ratios
!  TANHT   Tangent height for the current SZA
!  XL      Slant path between points
!  AMF     Air mass factor for slab between level and level above
! 
! 
! !REVISION HISTORY: 
!         1997 - O. Wild     - Initial version
!  31 Jul 2012 - R. Yantosca - Added ProTeX headers
!  10 Aug 2012 - R. Yantosca - Replace LPAR with LLPAR
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer i, j, k, ii
      real*8 airmas, gmu, xmu1, xmu2, xl, diff
      REAL*8 Ux,H,RZ(NB),RQ(NB),ZBYR
c
c  Inlined air mass factor function for top of atmosphere
      AIRMAS(Ux,H) = (1.0d0+H)/SQRT(Ux*Ux+2.0d0*H*(1.0d0-
     $         0.6817d0*EXP(-57.3d0*ABS(Ux)/SQRT(1.0d0+5500.d0*H))/
     $                                             (1.0d0+0.625d0*H)))
c
      GMU = U0
      RZ(1)=RAD+Z(1)
      ZBYR = ZZHT/RAD
      DO 2 II=2,NB
        RZ(II) = RAD + Z(II)
        RQ(II-1) = (RZ(II-1)/RZ(II))**2
    2 CONTINUE
      IF (GMU.LT.0.0D0) THEN
        TANHT = RZ(nlbatm)/DSQRT(1.0D0-GMU**2)
      ELSE
        TANHT = RZ(nlbatm)
      ENDIF
c
c  Go up from the surface calculating the slant paths between each level
c  and the level above, and deriving the appropriate Air Mass Factor
      DO 16 J=1,NB
        DO K=1,NB
          AMF(K,J)=0.D0
        ENDDO
c
c  Air Mass Factors all zero if below the tangent height
        IF (RZ(J).LT.TANHT) GOTO 16
c  Ascend from layer J calculating AMFs
        XMU1=ABS(GMU)
        DO 12 I=J,LLPAR
          XMU2=DSQRT(1.0D0-RQ(I)*(1.0D0-XMU1**2))
          XL=RZ(I+1)*XMU2-RZ(I)*XMU1
          AMF(I,J)=XL/(RZ(I+1)-RZ(I))
          XMU1=XMU2
   12   CONTINUE
c  Use function and scale height to provide AMF above top of model
        AMF(NB,J)=AIRMAS(XMU1,ZBYR)
c
c  Twilight case - Emergent Beam
        IF (GMU.GE.0.0D0) GOTO 16
        XMU1=ABS(GMU)
c  Descend from layer J
        DO 14 II=J-1,1,-1
          DIFF=RZ(II+1)*DSQRT(1.0D0-XMU1**2)-RZ(II)
          if(II.eq.1) DIFF=max(DIFF,0.d0)   ! filter
c  Tangent height below current level - beam passes through twice
          IF (DIFF.LT.0.0D0) THEN
            XMU2=DSQRT(1.0D0-(1.0D0-XMU1**2)/RQ(II))
            XL=ABS(RZ(II+1)*XMU1-RZ(II)*XMU2)
            AMF(II,J)=2.d0*XL/(RZ(II+1)-RZ(II))
            XMU1=XMU2
c  Lowest level intersected by emergent beam
          ELSE
            XL=RZ(II+1)*XMU1*2.0D0
c            WTING=DIFF/(RZ(II+1)-RZ(II))
c            AMF(II,J)=(1.0D0-WTING)*2.D0**XL/(RZ(II+1)-RZ(II))
            AMF(II,J)=XL/(RZ(II+1)-RZ(II))
            GOTO 16
          ENDIF
   14   CONTINUE
c
   16 CONTINUE
      RETURN
      END SUBROUTINE SPHERE
!EOC
