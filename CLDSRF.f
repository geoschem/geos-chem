C $Id: CLDSRF.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE CLDSRF( ODCOL, SA )
C-----------------------------------------------------------------------
c  Routine to set cloud and surface properties
C-----------------------------------------------------------------------
C  Add the following input variables for CTM interface (bmy, 9/13/99)
C
C  Variable  Type    Dimensn Units   Description
C  --------  ----    ------- -----   -----------
C  ODCOL     dble     [LPAR]   -     Vertical optical depth profile
C  SA        dble       -      -     Surface Albedo
C-----------------------------------------------------------------------
c     rflect   Surface albedo (Lambertian)
c     odmax    Maximum allowed optical depth, above which they are scaled
c     odcol    Optical depth at each model level
c     odsum    Column optical depth
c     nlbatm   Level of lower photolysis boundary - usually surface
C-----------------------------------------------------------------------
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"

C=============== INPUT PARAMETERS ======================================
      REAL*8, INTENT(INOUT) :: ODCOL(LPAR)
      REAL*8, INTENT(IN)    :: SA

C=============== LOCAL VARIABLES =======================================
      integer i, j, k
      real*8  odsum, odmax, odtot
c
c Default lower photolysis boundary as bottom of level 1
      nlbatm = 1
c
c Set surface albedo
      RFLECT = dble(SA)
      RFLECT = max(0.d0,min(1.d0,RFLECT))
c
c Zero aerosol column
      do k=1,MX
        do i=1,NB
          AER(k,i) = 0.d0
        enddo
      enddo
c
c Scale optical depths as appropriate - limit column to 'odmax'
      odmax = 200.d0
      odsum =   0.d0
      do i=1,lpar
        odcol(i) = dble(odcol(i))
        odsum    = odsum + odcol(i)
      enddo
      if(odsum.gt.odmax) then
        odsum = odmax/odsum
        do i=1,lpar
          odcol(i) = odcol(i)*odsum
        enddo
        odsum = odmax
      endif
c
c  Use clear-sky conditions
c      do i=1,jpnl
c        odcol(i)=0.d0
c      enddo
c
c Set sub-division switch if appropriate
      odtot=0.d0
      jadsub(nb)=0
      jadsub(nb-1)=0
      do i=nb-1,1,-1
        k=2*i
        jadsub(k)=0
        jadsub(k-1)=0
        odtot=odtot+odcol(i)
        if(odtot.gt.0.d0.and.odcol(i).ne.0.d0.and.
     $                                     dtausub.gt.0.d0) then
          if(odtot.le.dtausub) then
            jadsub(k)=1
            jadsub(k-1)=1
          elseif(odtot.gt.dtausub) then
            jadsub(k)=1
            jadsub(k-1)=0
            do j=1,2*(i-1)
              jadsub(j)=0
            enddo
            go to 20
          endif
        endif
      enddo
 20   continue
c
      return
      end
