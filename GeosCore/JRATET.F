C $Id: JRATET.f,v 1.1 2009/09/16 14:06:48 bmy Exp $
      SUBROUTINE JRATET( T, IDAY )
C-----------------------------------------------------------------------
c  Calculate and print J-values. Note that the loop in this routine 
c  only covers the jpnl levels actually needed by the CTM.
C-----------------------------------------------------------------------
C  Add the following input variables for CTM interface (bmy, 9/7/99)
C
C  Variable  Type    Dimensn Units   Description
C  --------  ----    ------- -----   -----------
C  T         dble     [LPAR]  [K]    Vertical temperature profile
C  IDAY      int        -      -     Day of Year (0-365 or 0-366)
C-----------------------------------------------------------------------
c
c     FFF    Actinic flux at each level for each wavelength bin
c     QQQ    Cross sections for species (read in in RD_TJPL)
c     SOLF   Solar distance factor, for scaling; normally given by:
c                      1.0-(0.034*cos(real(iday-172)*2.0*pi/365.))
c     TQQ    Temperatures at which QQQ cross sections supplied
c
c  NOTES
c  (1 ) Added a pressure-dependancy function selector 'pdepf' 
c        in 'jv_spec.dat'. (tmf, 1/7/09)
c  (2 ) Added pressure dependency for MGLY. (tmf, 1/7/09)
c  (3 ) Updated pressure dependency algorithm for ACET. (tmf, 1/7/09)
c  (4 ) Added pressure dependancy for MeCOVi, EtCOMe, MeCOCHO. Rewritten 
c        pressure dependancy for Acetone according to FAST-JX v6.4.
c        See more detailed documentation for Acetone in fjx_acet_mod.f.
c        (ccc, 4/20/09)
C-----------------------------------------------------------------------

      USE FJX_ACET_MOD

#     include "define.h"
      
      USE CMN_FJ_MOD, ONLY : IPAR, JPAR, LPAR, JPMAX, JPPJ, JPNL
      USE JV_CMN_MOD

      IMPLICIT NONE
#     include "define.h"

C=============== INPUT PARAMETERS ======================================
      REAL*8,  INTENT(IN) :: T(LPAR)
      INTEGER, INTENT(IN) :: IDAY

C=============== LOCAL VARIABLES =======================================

C     Add Pressure dependancy function selector PF. (tmf, 1/7/09) 
      integer i, j, k, l, PF
      real*8 qptemp

C     For new pressure-dependency algorithm: (tmf, 1/7/09) 
      real*8 xp, xa, xb, xc

C     For new pressure dependency algo. for acetone
C     All variables "*_F" are results from external functions from
C     fjx_acet_mod.f (ccc, 4/20/09)
      real*8 TFACA
      real*8 TFAC0
      real*8 TFAC1, TFAC2
      real*8 QQQA , QQ1A , QQ1B
      real*8 QQ2

      real*8 qo2tot, qo3tot, qo31d, qo33p, qqqt
      real*8 xseco2, xseco3, xsec1d, solf, tfact

C     Parameters for Solar distance compensation
      real*8  PI, TWOPI
      PARAMETER (PI=3.14159265358979324D0,TWOPI=2.*PI)

C     Physical constants
      REAL*8  Na, R
      PARAMETER (Na=6.02217d23, R=8.3143d0)

C     Scale actinic flux (FFF) by Solar distance factor (SOLF)
      solf=1.d0-(0.034d0*cos(dble(iday-172)*2.d0*pi/365.d0))
C----------------------------------------------------------------------
C If you want to set SOLF = 1.0 for testing, uncomment the next line
C      SOLF = 1d0
C----------------------------------------------------------------------
C
      do I=1,jpnl
       VALJ(1) = 0.d0
       VALJ(2) = 0.d0
       VALJ(3) = 0.d0
       do K=NW1,NW2                       ! Using model 'T's here
         QO2TOT= XSECO2(K,dble(T(I))) 
         VALJ(1) = VALJ(1) + QO2TOT*FFF(K,I)
         QO3TOT= XSECO3(K,dble(T(I)))
         QO31D = XSEC1D(K,dble(T(I)))*QO3TOT
         QO33P = QO3TOT - QO31D
         VALJ(2) = VALJ(2) + QO33P*FFF(K,I)
         VALJ(3) = VALJ(3) + QO31D*FFF(K,I)
       enddo
C------Calculate remaining J-values with T-dep X-sections 
       do J=4,NJVAL
         VALJ(J) = 0.d0
         TFACT = 0.d0
         L = jpdep(J)

C        To choose different forms of pres. dependancy. (ccc, 4/20/09)
         if ( L.ne.0 ) PF = pdepf(L)

         if(TQQ(2,J).gt.TQQ(1,J)) TFACT = max(0.d0,min(1.d0,
     $        (T(I)-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)) ))

C        FAST_JX introduces a new pres. dependancy for acetone (ccc, 4/20/09)
C        Special calculations for the temperature interpolation factors
         if ( PF.eq.2 ) then
            TFACA=TFACA_F(dble(T(I)), J      )
            TFAC0=TFAC0_F(dble(T(I)), J+1    )
            TFAC1=TFAC_F (dble(T(I)), NJVAL+1)
            TFAC2=TFAC_F (dble(T(I)), NJVAL+2)
         else if ( PF.eq.3 ) then
            TFACA=TFACA_F(dble(T(I)), J-1    )
            TFAC0=TFAC0_F(dble(T(I)), J      )
         endif

         do K=NW1,NW2
           QQQT = QQQ(K,1,J-3) + (QQQ(K,2,J-3) - QQQ(K,1,J-3))*TFACT
           if(L.eq.0) then
             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)
           else

              ! Select pressure dependancy function (tmf, 1/31/06)
              if (PF .eq. 1) then
C----------------------------------------------------------------------
C Prior to 9/17/99
C Original form for acetaldehyde P-dep -- believed to be incorrect (pjc)
C             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)*
C     $                   (1.d0+zpdep(K,L)*(pj(i)+pj(i+1))*0.5d0)
C----------------------------------------------------------------------
C Essentially the change is the replacement of the factor
C
C   (1 + a P)     with               1
C                           ---------------------
C                             (1 + b density)
C
C where a and b are constants, P is pressure, and density is the 
C density of air in molec-cm(-3)   (pjc, 9/17/99)
C----------------------------------------------------------------------
              VALJ(J)=VALJ(J)+QQQT*FFF(K,I)/(1 + 
     $                 (zpdep(K,L)*Na*1d-6 /(R*T(I))) * 
     $                 (pj(i)+pj(i+1))*0.5d0*1d2)

             else if ( PF .eq. 4 ) then
C-----------------------------------------------------------------------
C For MGLY
C       y = a + ( b * exp(-p/c) )
C    where y is the ratio between Omega(p) / Omega(p=0);
C          x is the atmospheric pressure [Pa]
C          a,b,c are MGLYPDEP(:,1), MGLYPDEP(:,2), MGLYPDEP(:,3)
C-----------------------------------------------------------------------
                 xp = (pj(i)+pj(i+1))*0.5d0*1.d2   ! pressure [Pa]
                 xa = mglypdep( K, 1 )
                 xb = mglypdep( K, 2 )
                 xc = mglypdep( K, 3 )
                 qptemp = 1.0d0

                 if ( abs( xc ) .ge. 1.d-10 ) then
                    qptemp = xa + ( xb * exp(-xp/xc) )
                 endif

                 VALJ(J) = VALJ(J) + QQQT*FFF(K,I)*qptemp

              else if ( PF.eq.2 ) then
C             Acetone pressure dependency from FAST-JX (ccc, 4/20/09)
C             J1(acetone-a) ==> CH3CO + CH3
C             Special values for Xsect
                 QQQA = QQ1_F (TFACA, J      , K            )
                 QQ2  = QQ2_F (TFAC0, J+1    , K, dble(T(I)))
                 QQ1A = QQ1_F (TFAC1, NJVAL+1, K            )
                 QQ1B = QQ1_F (TFAC2, NJVAL+2, K            ) * 4.d-20

                 VALJ(J) = VALJ(J) + FFF(K,L)*QQQA *
     &            (1.d0-QQ2)/(QQ1A + (QQ1B*Na*1d-6 /(R*T(I))) * 
     $            (pj(i)+pj(i+1))*0.5d0*1d2)
              else if ( PF.eq.3 ) then
C             Second acetone pressure dependency from FAST-JX (ccc, 4/20/09)
C             J2(acetone-b) ==> CH3 + CO + CH3
C             Special values for Xsect
                 QQQA = QQ1_F (TFACA, J-1    , K            )
                 QQ2  = QQ2_F (TFAC0, J      , K, dble(T(I)))

                 VALJ(J) = VALJ(J) + FFF(K,L)*QQQA*QQ2
              endif
           endif
         enddo
       enddo
       do j=1,jppj
         zj(i,j)=VALJ(jind(j))*jfacta(j)*SOLF
       enddo
cc       write(6,'(I5,1P,7E10.3/(5X,7E10.3))') I, (VALJ(J), J=1,NJVAL)
      enddo
      return
      end
