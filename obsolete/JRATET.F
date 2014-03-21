!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: jratet
!
! !DESCRIPTION: Subroutine JRATET calculates and prints J-values. Note that 
!  the loop in this routine only covers the jpnl levels actually needed by 
!  the CTM.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE JRATET( T, IDAY )
!
! !USES:
!
      USE FJX_ACET_MOD

      USE CMN_FJ_MOD,   ONLY : JPMAX, JPPJ, JPNL
      USE CMN_SIZE_MOD, ONLY : LLPAR
      USE JV_CMN_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: T(LLPAR)   ! Temperature [K]
      INTEGER, INTENT(IN) :: IDAY       ! Day of year (0-365 or 0-366)
!
! !REMARKS:
!  FFF    Actinic flux at each level for each wavelength bin
!  QQQ    Cross sections for species (read in in RD_TJPL)
!  SOLF   Solar distance factor, for scaling; normally given by:
!                   1.0-(0.034*cos(real(iday-172)*2.0*pi/365.))
!  TQQ    Temperatures at which QQQ cross sections supplied
!
! !REVISION HISTORY: 
!         1997 - O. Wild     - Initial version
!  (1 ) Added a pressure-dependancy function selector 'pdepf' 
!        in 'jv_spec.dat'. (tmf, 1/7/09)
!  (2 ) Added pressure dependency for MGLY. (tmf, 1/7/09)
!  (3 ) Updated pressure dependency algorithm for ACET. (tmf, 1/7/09)
!  (4 ) Added pressure dependancy for MeCOVi, EtCOMe, MeCOCHO. Rewritten 
!        pressure dependancy for Acetone according to FAST-JX v6.4.
!        See more detailed documentation for Acetone in fjx_acet_mod.f.
!        (ccc, 4/20/09)
!  25 Aug 2011 - R. Yantosca - Rewrite IF statement to prevent PF from
!                              never being initialized.
!  31 Jul 2012 - R. Yantosca - Added ProTeX headers
!  10 Aug 2012 - R. Yantosca - Replace LPAR with LLPAR
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!     Add Pressure dependancy function selector PF. (tmf, 1/7/09) 
      integer i, j, k, l, PF
      real*8 qptemp

!     For new pressure-dependency algorithm: (tmf, 1/7/09) 
      real*8 xp, xa, xb, xc

!     For new pressure dependency algo. for acetone
!     All variables "*_F" are results from external functions from
!     fjx_acet_mod.f (ccc, 4/20/09)
      real*8 TFACA
      real*8 TFAC0
      real*8 TFAC1, TFAC2
      real*8 QQQA , QQ1A , QQ1B
      real*8 QQ2

      real*8 qo2tot, qo3tot, qo31d, qo33p, qqqt
      real*8 xseco2, xseco3, xsec1d, solf, tfact

!     Parameters for Solar distance compensation
      real*8  PI, TWOPI
      PARAMETER (PI=3.14159265358979324D0,TWOPI=2.*PI)

!     Physical constants
      REAL*8  Na, R
      PARAMETER (Na=6.02217d23, R=8.3143d0)

!     Scale actinic flux (FFF) by Solar distance factor (SOLF)
      solf=1.d0-(0.034d0*cos(dble(iday-172)*2.d0*pi/365.d0))
!----------------------------------------------------------------------
! If you want to set SOLF = 1.0 for testing, uncomment the next line
!      SOLF = 1d0
!----------------------------------------------------------------------
!
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
!------Calculate remaining J-values with T-dep X-sections 
       do J=4,NJVAL
         VALJ(J) = 0.d0
         TFACT = 0.d0
         L = jpdep(J)

!        To choose different forms of pres. dependancy. (ccc, 4/20/09)
         if ( L.ne.0 ) then
            PF = pdepf(L)
         else
            PF = -1
         endif

         if(TQQ(2,J).gt.TQQ(1,J)) TFACT = max(0.d0,min(1.d0,
     $        (T(I)-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)) ))

!        FAST_JX introduces a new pres. dependancy for acetone (ccc, 4/20/09)
!        Special calculations for the temperature interpolation factors
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
!----------------------------------------------------------------------
! Prior to 9/17/99
! Original form for acetaldehyde P-dep -- believed to be incorrect (pjc)
!             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)*
!     $                   (1.d0+zpdep(K,L)*(pj(i)+pj(i+1))*0.5d0)
!----------------------------------------------------------------------
! Essentially the change is the replacement of the factor
!
!   (1 + a P)     with               1
!                           ---------------------
!                             (1 + b density)
!
! where a and b are constants, P is pressure, and density is the 
! density of air in molec-cm(-3)   (pjc, 9/17/99)
!----------------------------------------------------------------------
              VALJ(J)=VALJ(J)+QQQT*FFF(K,I)/(1 + 
     $                 (zpdep(K,L)*Na*1d-6 /(R*T(I))) * 
     $                 (pj(i)+pj(i+1))*0.5d0*1d2)

             else if ( PF .eq. 4 ) then
!-----------------------------------------------------------------------
! For MGLY
!       y = a + ( b * exp(-p/c) )
!    where y is the ratio between Omega(p) / Omega(p=0);
!          x is the atmospheric pressure [Pa]
!          a,b,c are MGLYPDEP(:,1), MGLYPDEP(:,2), MGLYPDEP(:,3)
!-----------------------------------------------------------------------
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
!             Acetone pressure dependency from FAST-JX (ccc, 4/20/09)
!             J1(acetone-a) ==> CH3CO + CH3
!             Special values for Xsect
                 QQQA = QQ1_F (TFACA, J      , K            )
                 QQ2  = QQ2_F (TFAC0, J+1    , K, dble(T(I)))
                 QQ1A = QQ1_F (TFAC1, NJVAL+1, K            )
                 QQ1B = QQ1_F (TFAC2, NJVAL+2, K            ) * 4.d-20

                 VALJ(J) = VALJ(J) + FFF(K,L)*QQQA *
     &            (1.d0-QQ2)/(QQ1A + (QQ1B*Na*1d-6 /(R*T(I))) * 
     $            (pj(i)+pj(i+1))*0.5d0*1d2)
              else if ( PF.eq.3 ) then
!             Second acetone pressure dependency from FAST-JX (ccc, 4/20/09)
!             J2(acetone-b) ==> CH3 + CO + CH3
!             Special values for Xsect
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
      END SUBROUTINE JRATET
!EOC
