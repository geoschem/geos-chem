      subroutine CO_set_prof( Pcol, MONTH, YLAT, EDGES, PTOPCO )
C-----------------------------------------------------------------------
c  Routine to read and integrate atmospheric profiles required by Fast-J.
c  O3 and T are taken from the supplied climatology and integrated to 
c  the CTM levels.  (Modified from SR set_prof.f, bnd)
C-----------------------------------------------------------------------
C  Add the following input variables for CTM interface (bmy, 9/13/99)
C
C  Variable  Type    Dimensn Units   Description
C  --------  ----    ------- -----   -----------
C  P         dble       -     [mb]   Surface Pressure - Model top pressure
C  T         dble     [LPAR]  [K]    Vertical temperature profile
C
C-----------------------------------------------------------------------
c
c     PJ       Pressure at boundaries of model levels (hPa)
c     masfac   Conversion factor for pressure to column density
c
c     TJ       Temperature profile on model grid
c     DM       Air density for each model level (molecules.cm-2)
c     DO3      Ozone profile for each model level (molecules.cm-2)
c     PSTD     Approximate pressures of levels for supplied climatology
c
C-----------------------------------------------------------------------
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"

C=============== INPUT PARAMETERS ======================================
      INTEGER, INTENT(IN)    :: MONTH
      REAL*8,  INTENT(IN)    :: Pcol, YLAT, EDGES(NB), PTOPCO

C=============== LOCAL VARIABLES =======================================
      integer i, k, l, m
      real*8  dlogp,F0,T0,PB,PC,XC,masfac
      real*8  pstd(52),oref2(51),tref2(51)
c
      print*,'Calculate above O3 column using fastj columns'
      print*,'in SR CO_set_prof.'
c
c  The PJ's are the pressures at each sigma edge
c  Replace Oliver's ETAA-ETAB system with a pure sigma system (bmy, 9/13/99)
      DO I = 1, NB
         PJ(I) = PTOPCO + ( EDGES(I) * Pcol )
      ENDDO
      PJ(NB+1) = 0.d0
c
c  Set up pressure levels for O3/T climatology - assume that value
c  given for each 2 km z* level applies from 1 km below to 1 km above,
c  so select pressures at these boundaries. Surface level values at
c  1000 mb are assumed to extend down to the actual P(nslon,nslat).
c
      pstd(1) = max(PJ(1),1000.d0)
      pstd(2) = 1000.d0*10.d0**(-1.d0/16.d0)
      dlogp = 10.d0**(-2.d0/16.d0)
      do i=3,51
        pstd(i) = pstd(i-1)*dlogp
      enddo
      pstd(52) = 0.d0
c
c  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
      masfac=100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)
c
c  Select appropriate MONTHly and latitudinal profiles
c  Now use YLAT instead of Oliver's YDGRD(NSLAT) (bmy, 9/13/99) 
      m = max(1,min(12,MONTH))
      l = max(1,min(18,(int(ylat)+99)/10))
c
c  Temporary arrays for climatology data
      do i=1,51
        oref2(i)=oref(i,l,m)
        tref2(i)=tref(i,l,m)
      enddo
c
c  Apportion O3 and T on supplied climatology z* levels onto CTM levels 
c  with mass (pressure) weighting, assuming constant mixing ratio and
c  temperature half a layer on either side of the point supplied.
c
      do i = 1,NB
        F0 = 0.d0
        T0 = 0.d0
        do k = 1,51
          PC = min(PJ(i),pstd(k))
          PB = max(PJ(i+1),pstd(k+1))
          if(PC.gt.PB) then
            XC = (PC-PB)/(PJ(i)-PJ(i+1))
            F0 = F0 + oref2(k)*XC
            T0 = T0 + tref2(k)*XC
          endif
        enddo
        TJ(i) = T0
        DO3(i)= F0*1.d-6
      enddo
c
c  Calculate column quantities for Fast-J
c
      do i=1,NB
        DM(i)  = (PJ(i)-PJ(i+1))*masfac
        DO3(i) = DO3(i)*DM(i)
      enddo
c
      return
      end
