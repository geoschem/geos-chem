#ifdef APM
!************************************************************************
! This computer software was prepared by Battelle Memorial Institute,
! hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with
! the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
! CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
! LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! MOSAIC module: see chem/module_mosaic_driver.F for references and terms
! of use
!************************************************************************
!
!      Calculate maxsat at given updraft velocity and particle size dist
!      modified from relevant WRF_Chem code for APM
!      by F. Yu, UAlbany, 2012/11/19

!      Reference:
!      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
!      3. Sectional representation. J. Geophys. Res., 107, 2002.

MODULE apm_mixactivate
PRIVATE
PUBLIC activate
CONTAINS

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! 06-nov-2005 rce - grid_id & ktau added to arg list
      subroutine activate(maxsatout, wbar, tair, pres,&
                      maxd_atype, ntype_aer, maxd_asize, nsize_aer,    &
                      na, am, hygro)

      implicit none

      integer,intent(in) :: maxd_atype      ! dimension of types
      integer,intent(in) :: maxd_asize      ! dimension of sizes
      integer,intent(in) :: ntype_aer       ! number of types
      integer,intent(in) :: nsize_aer(maxd_atype) ! number of sizes for type
      real*8,intent(in) :: wbar          ! grid cell mean vertical velocity (m/s)
      real*8,intent(in) :: tair          ! air temperature (K)
      real*8,intent(in) :: pres         ! pressure (Pa)
      real*8,intent(in) :: na(maxd_asize,maxd_atype)     ! aerosol number concentration (/m3)
      real*8,intent(in) :: am(maxd_asize,maxd_atype)   ! aerosol section dry radius (m)
      real*8,intent(in) :: hygro(maxd_asize,maxd_atype)  ! bulk hygroscopicity of aerosol mode

      real*8,intent(inout) :: maxsatout


      real*8,parameter :: g=9.81,rhowater=1000.,XLV=2.5E6,r_d=287.,r_v=461.6,mwdry=28.966
      real*8,parameter :: cp=7.*r_d/2.,rvovrd=r_v/r_d,SVP1=0.6112,SVP2=17.67,SVP3=29.65
      real*8,parameter :: EP_2=R_d/R_v

      real*8, parameter :: surften = 0.076 ! surface tension of water w/respect to air (N/m)
      real*8, parameter :: p0 = 1013.25e2  ! reference pressure (Pa)
      real*8, parameter :: t0 = 273.15     ! reference temperature (K)

      real*8 :: rhoair        ! air density (kg/m3)
      real*8 diff0 ! diffusivity (m2/s)
      real*8 conduct0 ! thermal conductivity (Joule/m/sec/deg)
      real*8 es ! saturation vapor pressure
      real*8 qs ! water vapor saturation mixing ratio
      real*8 dqsdt ! change in qs with temperature
      real*8 gg ! thermodynamic function (m2/s)
      real*8 sqrtg ! sqrt(gg)
      real*8 sm(maxd_asize,maxd_atype) ! critical supersaturation for number mode radius
      real*8 zeta, eta
      real*8 alpha
      real*8 gamma
      real*8 beta
      real*8 totn ! total aerosol number concentration
      real*8 aten ! surface tension parameter
      real*8 gmsm ! critical supersaturation at radius gmrad
      real*8 sumns

      real*8 alw,sqrtalw
      real*8 smax
      integer m,n

!      mathematical constants
      real*8 third, twothird, sixth, zero, one, two, three

      real*8, parameter :: sq2  = 1.4142135624
      real*8, parameter :: sqpi = 1.7724538509
      real*8, parameter :: pi   = 3.1415926536
      real*8, parameter :: nsmall = 1.0e-20    ! aer number conc in #/m3
      real*8, parameter :: amsmall = 5.0e-9    ! aer dry radius in m

      zero = 0.0
      one = 1.0
      two = 2.0
      three = 3.0
      third = 1.0/3.0
      twothird = 2.0/3.0 !wig, 1-Mar-2009: Corrected value from 2/6
      sixth = 1.0/6.0

      rhoair = pres/(r_d*tair)

      diff0=0.211e-4*(p0/pres)*(tair/t0)**1.94
      conduct0=(5.69+0.017*(tair-t0))*4.186e2*1.e-5 ! convert to J/m/s/deg
      es=1000.*svp1*exp( svp2*(tair-t0)/(tair-svp3) )
      qs=ep_2*es/(pres-es)
      dqsdt=xlv/(r_v*tair*tair)*qs
      alpha=g*(xlv/(cp*r_v*tair*tair)-1./(r_d*tair))
      gamma=(1+xlv/cp*dqsdt)/(rhoair*qs)
      gg=1./(rhowater/(diff0*rhoair*qs)+xlv*rhowater/(conduct0*tair)*(xlv/(r_v*tair)-1.))
      sqrtg=sqrt(gg)
      beta=4.*pi*rhowater*gg*gamma
      aten=2.*surften/(r_v*tair*rhowater)

      totn=1.d-30
      sumns=1.d-30
      do n=1,ntype_aer
      do m=1,nsize_aer(n)
!      internal mixture of aerosols

         if (am(m,n).gt.amsmall .and. na(m,n).gt.nsmall) then
!            sectional model.
!            need to use bulk properties because parameterization doesn't
!            work well for narrow bins.
         totn=totn+na(m,n)

         if(hygro(m,n).gt.1.d-10)then
            sm(m,n)=2.d0*aten/(3.*am(m,n))*sqrt(aten/(3.*hygro(m,n)*am(m,n)))
         else
            sm(m,n)=100.d0
         endif
!        write(6,*)'sm,hygro,am=',sm(m,n),hygro(m,n),am(m,n)
         else
            sm(m,n)=1.d0
         endif
         sumns=sumns+na(m,n)/sm(m,n)**twothird
      end do ! size
      end do ! type

      gmsm=totn/sumns
      gmsm=gmsm*gmsm*gmsm
      gmsm=sqrt(gmsm)

!         write(6,*)'uniform updraft =',wbar

      alw=alpha*wbar
      sqrtalw=sqrt(alw)
      zeta=2.*sqrtalw*aten/(3.*sqrtg)

!     sectional model.
!     use bulk properties
      if(totn.gt.1.d-10)then
        eta=2.d0*alw*sqrtalw/(totn*beta*sqrtg)
      else
        eta=1.d10
      endif

      call maxsat(zeta,eta,gmsm,smax)

      maxsatout=smax

      return
      end subroutine activate

!----------------------------------------------------------------------
      subroutine maxsat(zeta,eta,gmsm,smax)

!      Calculates maximum supersaturation for multiple competing aerosol
!      sections/types.

      implicit none

      real*8, intent(in)  :: gmsm ! critical supersaturation for number mode radius
      real*8, intent(in)  :: zeta, eta
      real*8, intent(out) :: smax ! maximum supersaturation

      real*8 :: g1, g2
      real*8 thesum
      integer n ! type index

      if(zeta.gt.1.d5*eta .or. &
        gmsm*gmsm.gt.1.d5*eta)then
!       weak forcing. essentially none activated
        smax=1.d-20
      else
!       significant activation of this mode. calc activation all modes.
        go to 1
      endif

      return

  1   continue

      thesum=0.d0
      if(eta.gt.1.d-20)then
        g1=sqrt(zeta/eta)
        g1=g1*g1*g1
        g2=gmsm/sqrt(eta+3.d0*zeta)
        g2=sqrt(g2)
        g2=g2*g2*g2
        thesum=thesum + (0.5d0*g1 + g2)/(gmsm*gmsm)
      else
        thesum=1.d20
      endif

      smax=1./sqrt(thesum)

      return

      end subroutine maxsat
!----------------------------------------------------------------------

END MODULE apm_mixactivate
#endif
