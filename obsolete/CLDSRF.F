! $Id: CLDSRF.f,v 1.1 2009/09/16 14:06:49 bmy Exp $
      SUBROUTINE CLDSRF( ODCOL, SA )
!-----------------------------------------------------------------------
c  Routine to set cloud and surface properties
!-----------------------------------------------------------------------
!  Add the following input variables for CTM interface (bmy, 9/13/99)
!
!  Variable  Type    Dimensn Units   Description
!  --------  ----    ------- -----   -----------
!  ODCOL     dble     [LPAR]   -     Vertical optical depth profile
!  SA        dble       -      -     Surface Albedo
!-----------------------------------------------------------------------
c     rflect   Surface albedo (Lambertian)
c     odmax    Maximum allowed optical depth, above which they are scaled
c     odcol    Optical depth at each model level
c     odsum    Column optical depth
c     nlbatm   Level of lower photolysis boundary - usually surface
!  10 Aug 2012 - R. Yantosca - Replace LPAR with LLPAR
!  20 Aug 2013 - R. Yantosca - Remove "define.h", it's now obsolete
!-----------------------------------------------------------------------

      USE CMN_FJ_MOD, ONLY : LLPAR, JPMAX, JPPJ
      USE JV_CMN_MOD

      IMPLICIT NONE

!=============== INPUT PARAMETERS ======================================
      REAL*8, INTENT(INOUT) :: ODCOL(LLPAR)
      REAL*8, INTENT(IN)    :: SA

!=============== LOCAL VARIABLES =======================================
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
      do i=1,LLPAR
        odcol(i) = dble(odcol(i))
        odsum    = odsum + odcol(i)
      enddo
      if(odsum.gt.odmax) then
        odsum = odmax/odsum
        do i=1,LLPAR
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
