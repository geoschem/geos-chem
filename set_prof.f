! $Id: set_prof.f,v 1.1 2003/06/30 20:26:07 bmy Exp $
      subroutine set_prof( P,     T,    ODCOL,   SA, 
     &                     MONTH, YLAT, OPTDUST, OPTAER )
C-----------------------------------------------------------------------
c  Routine to set up atmospheric profiles required by Fast-J using a
c  doubled version of the level scheme used in the CTM. First pressure
c  and z* altitude are defined, then O3 and T are taken from the supplied
c  climatology and integrated to the CTM levels (may be overwritten with
c  values directly from the CTM, if desired) and then black carbon and
c  aerosol profiles are constructed.
c                                                     Oliver (04/07/99)
c  Mineral dust profiles are also constructed (rvm, 06/04/00)
c  and other aerosol profiles are also constructed (rvm, bmy, 2/27/02)
C-----------------------------------------------------------------------
C  Add the following input variables for CTM interface (bmy, rvm, 9/30/00)
C
C  Variable  Type    Dimension   Units   Description
C  --------  ----    ---------   -----   -----------
C  P         dble       -         [mb]   Surface Pressure - Model top pressure
C  T         dble     [LPAR]      [K]    Vertical temperature profile
C  ODCOL     dble     [LPAR]       -     Vertical optical depth profile
C  SA        dble       -          -     Surface Albedo
C  OPTDUST   dble     [LPAR,NDUST] -     Mineral dust optical depths
C  OPTAER    dble  [LPAR,NAER*NRH] -     Other aerosols
C
C  Also note: since we parallelize over columns, make T and ODCOL
C  1-D vectors.  In the original code from O. Wild, these were 3-D
C  arrays (bmy, 9/13/99)
C
C  OPTDUST(LGLOB,NDUST) is now OPTDUST(LPAR,NDUST), where LPAR = LLPAR.
C  (bmy, 9/26/01)
C-----------------------------------------------------------------------
c
c     pj       Pressure at boundaries of model levels (hPa)
c     z        Altitude of boundaries of model levels (cm)
c     odcol    Optical depth at each model level
c     masfac   Conversion factor for pressure to column density
c
c     TJ       Temperature profile on model grid
c     DM       Air column for each model level (molecules.cm-2)
c     DO3      Ozone column for each model level (molecules.cm-2)
c     DBC      Mass of Black Carbon at each model level (g.cm-3)  !  .....!
c     PSTD     Approximate pressures of levels for supplied climatology
c
C-----------------------------------------------------------------------
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"  ! NDUST, NAER

C=============== INPUT PARAMETERS ======================================
      INTEGER, INTENT(IN)    :: MONTH
      REAL*8,  INTENT(IN)    :: P, T(LPAR), SA, YLAT
      REAL*8,  INTENT(INOUT) :: ODCOL(LPAR)
      REAL*8,  INTENT(IN)    :: OPTDUST(LPAR,NDUST)  
      REAL*8,  INTENT(IN)    :: OPTAER(LPAR,NAER*NRH)  

C=============== LOCAL VARIABLES =======================================
      integer i, k, l, m, n
      real*8  dlogp,f0,t0,b0,pb,pc,xc,masfac,scaleh
      real*8  pstd(52),oref2(51),tref2(51),bref2(51)
c
      ! Now use the hybrid pressure formulation (bmy, 8/22/02)
      DO I = 1, NB
         PJ(I) = ETAA(I) + ( ETAB(I) * P )
      ENDDO
      pj(NB+1) = 0.d0
c
c  Set up cloud and surface properties
      call CLDSRF( ODCOL, SA )
c
c  Set up pressure levels for O3/T climatology - assume that value
c  given for each 2 km z* level applies from 1 km below to 1 km above,
c  so select pressures at these boundaries. Surface level values at
c  1000 mb are assumed to extend down to the actual P(nslon,nslat).
c
      pstd(1) = max(pj(1),1000.d0)
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
c  Select appropriate monthly and latitudinal profiles
c  Now use YLAT instead of Oliver's YDGRD(NSLAT) (bmy, 9/13/99) 
      m = max(1,min(12,month))
      l = max(1,min(18,(int(ylat)+99)/10))
c
c  Temporary arrays for climatology data
      do i=1,51
        oref2(i)=oref(i,l,m)
        tref2(i)=tref(i,l,m)
        bref2(i)=bref(i)
      enddo
c
c  Apportion O3 and T on supplied climatology z* levels onto CTM levels 
c  with mass (pressure) weighting, assuming constant mixing ratio and
c  temperature half a layer on either side of the point supplied.
c
      do i = 1,NB
        F0 = 0.d0
        T0 = 0.d0
        B0 = 0.d0
        do k = 1,51
          PC = min(pj(i),pstd(k))
          PB = max(pj(i+1),pstd(k+1))
          if(PC.gt.PB) then
            XC = (PC-PB)/(pj(i)-pj(i+1))
            F0 = F0 + oref2(k)*XC
            T0 = T0 + tref2(k)*XC
            B0 = B0 + bref2(k)*XC
          endif
        enddo
        TJ(i) = T0
        DO3(i)= F0*1.d-6
        DBC(i)= B0
      enddo
c
c  Insert model values here to replace or supplement climatology.
c  Note that CTM temperature is always used in x-section calculations
c  (see JRATET); TJ is used in actinic flux calculation only.
c
c      do i=1,lpar
c        DO3(i) = my_ozone(i)        ! Volume Mixing Ratio
c        TJ(i)  = T(I)                ! Kelvin
c      enddo
c      DO3(lpar+1) = my_ozone*exp()  ! Above top of model (or use climatology)
c      TJ(lpar+1)  = my_temp(lpar)   ! Above top of model (or use climatology)
c
c
c  Calculate effective altitudes using scale height at each level
      z(1) = 0.d0
      do i=1,lpar
        scaleh=1.3806d-19*masfac*TJ(i)
        z(i+1) = z(i)-(log(pj(i+1)/pj(i))*scaleh)
      enddo
c
c  Add Aerosol Column - include aerosol types here. Currently use soot
c  water and ice; assume black carbon x-section of 10 m2/g, independent
c  of wavelength; assume limiting temperature for ice of -40 deg C.
c
      do i=1,lpar
        AER(1,i) = DBC(i)*10.d0*(z(i+1)-z(i))

        ! Turn off uniform black carbon profile (rvm, bmy, 2/27/02)
        AER(1,i) = 0d0

        if(T(I).gt.233.d0) then
          AER(2,i) = odcol(i)
          AER(3,i) = 0.d0
        else
          AER(2,i) = 0.d0
          AER(3,i) = odcol(i)
        endif   

        ! Also add in aerosol optical depth columns (rvm, bmy, 9/30/00)
        do n=1,ndust
          AER(3+n,i) = optdust(i,n)	
        enddo
        
        ! Also add in other aerosol optical depth columns (rvm, bmy, 2/27/02)
        do n = 1, NAER*NRH
           AER(3+n+NDUST,i) = optaer(i,n)
        enddo


      enddo
      do k=1,MX
        AER(k,lpar+1) = 0.d0
      enddo
c
c  Calculate column quantities for Fast-J
      do i=1,NB
        DM(i)  = (PJ(i)-PJ(i+1))*masfac
        DO3(i) = DO3(i)*DM(i)
      enddo
c
      return
      end
