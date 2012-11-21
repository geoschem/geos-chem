!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: vdiff_mod
!
! !DESCRIPTION: Module VDIFF\_MOD includes all routines for the non-local PBL
!  mixing scheme.
!\\
!\\
! !INTERFACE: 
!
MODULE VDIFF_MOD
! 
! !USES:
!
  USE TRACER_MOD,    ONLY : pcnst => N_TRACERS
  USE LOGICAL_MOD,   ONLY : LPRT
  USE ERROR_MOD,     ONLY : DEBUG_MSG
  USE VDIFF_PRE_MOD, ONLY : plev  => LLPAR
  USE CMN_SIZE_MOD,  ONLY : IIPAR, JJPAR, LLPAR

  IMPLICIT NONE
# include "define.h"
  
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: DO_PBL_MIX_2 
!
! !PRIVATE DATA MEMBERS:
!  
  save
  

  integer :: plevp
  
  real*8, parameter ::          &
       rearth = 6.37122d6,      & ! radius earth (m)
       cpwv   = 1.81d3,         &
       cpair  = 1004.64d0,      &
       rair   = 287.04d0,       &
       rh2o   = 461.d0,         &
       zvir   = rh2o/rair - 1., &
       gravit = 9.80616d0,      &
       ra     = 1.d0/rearth,    &
       epsilo = 0.622d0,        &
       latvap = 2.5104d06,      &
       latice = 3.336d5,        &
       cappa  = rair/cpair,     &
       rhoh2o = 1.d3,           &
       r_g    = rair / gravit,  &
       tfh2o  = 273.16d0

!-----------------------------------------------------------------------
! 	... pbl constants
!-----------------------------------------------------------------------

  ! These are constants, so use PARAMETER tag
  real*8, parameter ::   &
       betam  = 15.d0,   & ! constant in wind gradient expression
       betas  =  5.d0,   & ! constant in surface layer gradient expression
       betah  = 15.d0,   & ! constant in temperature gradient expression 
       fak    =  8.5d0,  & ! constant in surface temperature excess         
       fakn   =  7.2d0,  & ! constant in turbulent prandtl number
       ricr   =   .3d0,  & ! critical richardson number
       sffrac =   .1d0,  & ! surface layer fraction of boundary layer
       vk     =   .4d0     ! von karmans constant

  ! These are assigned later, so we can't use the PARAMETER tag
  real*8 ::              & 
       g,                & ! gravitational acceleration
       onet,             & ! 1/3 power in wind gradient expression
       ccon,             & ! fak * sffrac * vk
       binm,             & ! betam * sffrac
       binh                ! betah * sffrac

!-----------------------------------------------------------------------
! 	... constants used in vertical diffusion and pbl
!-----------------------------------------------------------------------
  real*8 :: &
       zkmin            ! minimum kneutral*f(ri)
!#if defined( DEVEL )
  real*8, allocatable :: ml2(:)   ! mixing lengths squaredB
!#else
!  real*8 :: ml2(plevp)   ! mixing lengths squared
!#endif
  real*8, allocatable :: qmincg(:)   ! min. constituent concentration 
                                     !  counter-gradient term
  
  integer :: &
       ntopfl, &        ! top level to which vertical diffusion is applied.
       npbl             ! maximum number of levels in pbl from surface

  logical, parameter :: divdiff = .true. , arvdiff = .false.
  
  logical, parameter :: pblh_ar = .true.
  
  logical, parameter :: pbl_mean_drydep = .false.  ! use mean concentration 
                                                   !  within the  PBL for 
                                                   ! calculating drydep fluxes
  logical, parameter :: drydep_back_cons = .false. ! backward consistency 
                                                   !  with previous GEOS-Chem 
                                                   !  drydep budgets 
                                                   !-- useless when 
                                                   !   pbl_mean_drydep=.false.
!
! !REMARKS:
!  The non-local PBL mixing routine VDIFF modifies the specific humidity,
!  (State_Met%SPHU) field.  Therefore, we must pass State_Met as an argument
!  to DO_PBL_MIX_2 and VDIFFDR with INTENT(INOUT).
!
! !REVISION HISTORY:
!  (1 ) This code is modified from mo_vdiff.F90 in MOZART-2.4. (lin, 5/14/09)
!  07 Oct 2009 - R. Yantosca - Added CVS Id Tag
!  24 Sep 2010 - J. Lin      - Modified ND15 to account for all mixing
!                              processes but not dry deposition and emissions.
!  17 Dec 2010 - R. Yantosca - Declare constants w/ the PARAMETER attribute
!  20 Dec 2010 - R. Yantosca - Bug fixes for the parallelization
!  02 Mar 2011 - R. Yantosca - Bug fixes for PGI compiler: these mostly
!                              involve explicitly using "D" exponents
!  25 Mar 2011 - R. Yantosca - Corrected bug fixes noted by Jintai Lin
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.2 met
!  22 Jun 2012 - R. Yantosca - Now use pointers to flip arrays in vertical
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: pbinti
!
! !DESCRIPTION: Subroutine PBINTI initializes time independent variables 
!  of pbl package
!\\
!\\
! !INTERFACE:
!
  subroutine pbinti( gravx )

!
! !USES:
!
    implicit none
!
! !INPUT PARAMETERS: 
!
    real*8, intent(in) :: gravx     !  acceleration of gravity
!
! !REVISION HISTORY: 
!  02 Mar 2011 - R. Yantosca - Bug fixes for PGI compiler: these mostly
!                              involve explicitly using "D" exponents
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! pbinti begins here!
    !=================================================================
!-----------------------------------------------------------------------
! 	... basic constants
!-----------------------------------------------------------------------
    plevp = plev+1
    g     = gravx
    onet  = 1d0/3.d0
    
!-----------------------------------------------------------------------
! 	... derived constants
!-----------------------------------------------------------------------
    ccon = fak*sffrac*vk
    binm = betam*sffrac
    binh = betah*sffrac
    
  end subroutine pbinti
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!
! !IROUTINE: vdiff
!
! !DESCRIPTION:
!  Subroutine vdiff is the driver routine to compute vertical diffusion of 
!  momentum, moisture, trace constituents and potential temperature. 
!\\
!\\
! !INTERFACE:
!
  subroutine vdiff( lat,        ip,        uwnd,       vwnd,        &
                    tadv,       pmid,      pint,       rpdel_arg,   &
                    rpdeli_arg, ztodt,     zm_arg,     shflx_arg,   &
                    sflx,       thp_arg,   as2,        pblh_arg,    &
                    kvh_arg,    kvm_arg,   tpert_arg,  qpert_arg,   &
                    cgs_arg,    shp,       wvflx_arg,  plonl,       &
                    State_Met,  taux_arg,   tauy_arg,  ustar_arg )
!
! !USES:
!
    USE DIAG_MOD,           ONLY : TURBFLUP
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE TRACER_MOD,         ONLY : TCVV
    USE VDIFF_PRE_MOD,      ONLY : ND15

    implicit none
!
! !INPUT PARAMETERS: 
!
    integer, intent(in) :: lat, ip ! latitude index, long tile index
    integer, intent(in) :: plonl   ! number of local longitudes
    real*8,  intent(in) ::   &
         ztodt                     ! 2 delta-t
    real*8,  intent(in) ::   &
         uwnd(:,:,:),        &     ! u wind input
         vwnd(:,:,:),        &     ! v wind input
         tadv(:,:,:),        &     ! temperature input
         pmid(:,:,:),        &     ! midpoint pressures
         pint(:,:,:),        &     ! interface pressures
         rpdel_arg(:,:,:),   &     ! 1./pdel  (thickness bet interfaces)
         rpdeli_arg(:,:,:),  &     ! 1./pdeli (thickness bet midpoints)
         zm_arg(:,:,:),      &     ! midpoint geoptl height above sfc
         shflx_arg(:,:),     &     ! surface sensible heat flux (w/m2)
         sflx(:,:,:),        &     ! surface constituent flux (kg/m2/s)
         wvflx_arg(:,:)            ! water vapor flux (kg/m2/s)
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    real*8, intent(inout) :: &
         as2(:,:,:,:),       &     ! moist, tracers after vert. diff
         shp(:,:,:),         &     ! specific humidity (kg/kg)
         thp_arg(:,:,:)            ! pot temp after vert. diffusion
!
! !OUTPUT PARAMETERS: 
!
    real*8, intent(out) ::   &
         kvh_arg(:,:,:),     &     ! coefficient for heat and tracers
         kvm_arg(:,:,:),     &     ! coefficient for momentum
         tpert_arg(:,:),     &     ! convective temperature excess
         qpert_arg(:,:),     &     ! convective humidity excess
         cgs_arg(:,:,:)            ! counter-grad star (cg/flux)

    real*8, optional, intent(inout) :: &
         taux_arg(:,:),      &     ! x surface stress (n)
         tauy_arg(:,:),      &     ! y surface stress (n)
         ustar_arg(:,:)            ! surface friction velocity

    real*8, intent(inout) :: pblh_arg(:,:) ! boundary-layer height [m]

!
! !REMARKS:
!  Free atmosphere diffusivities are computed first; then modified by the 
!  boundary layer scheme; then passed to individual parameterizations mvdiff, 
!  qvdiff.
!
!  The free atmosphere diffusivities are based on standard mixing length forms 
!  for the neutral diffusivity multiplied by functions of Richardson number. 
!  k = l^2 * |dv/dz| * f(ri). The same functions are used for momentum, 
!  potential temperature, and constitutents.
!
!  The stable Richardson num function (ri>0) is taken from Holtslag and 
!  Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*ri*(1 + 8*ri)).
!  The unstable richardson number function (ri<0) is taken from ccm1.
!  f = sqrt(1 - 18*ri)
!
! !REVISION HISTORY:
! (1 ) All arguments are full arrays now. Latitude slices are copied in local
!      variables. (ccc, 11/19/09)
!  24 Sep 2010 - J. Lin      - Moved call to ND15 at the end of vdiff.
!                              Modified to account for all mixing processes.
!  02 Mar 2011 - R. Yantosca - Bug fixes for PGI compiler: these mostly
!                              involve explicitly using "D" exponents
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: &
         i, &                   ! longitude index
         k, l, &                ! vertical index
         m                      ! constituent index
    integer :: &
         indx(plonl), &         ! array of indices of potential q<0
         nval, &                ! num of values which meet criteria
         ii                     ! longitude index of found points
    real*8 :: &
         dvdz2 , &                 ! (du/dz)**2 + (dv/dz)**2
         dz , &                    ! delta z between midpoints
         fstab , &                 ! stable f(ri)
         funst , &                 ! unstable f(ri)
         rinub, &                  ! richardson no=(g/theta)(dtheta/dz)/
                                   !               (du/dz**2+dv/dz**2)
         sstab, &                  ! static stability = g/th  * dth/dz
         kvn, &                    ! neutral kv
         tmp2, &                   ! temporary storage
         rcpair, &                 ! 1./cpair
         ztodtgor, &               ! ztodt*gravit/rair
         gorsq                     ! (gravit/rair)**2
    real*8 :: &
         cah(plonl,plev), &        ! -upper diag for heat and constituts
         cam(plonl,plev), &        ! -upper diagonal for momentum
         cch(plonl,plev), &        ! -lower diag for heat and constits
         ccm(plonl,plev), &        ! -lower diagonal for momentum
         cgh(plonl,plevp), &       ! countergradient term for heat
         cgq(plonl,plevp,pcnst), & ! countergrad term for constituent
         cgsh(plonl,plevp), &      ! countergrad term for sh
         kvf(plonl,plevp)          ! free atmosphere kv at interfaces
    real*8 :: &
         potbar(plonl,plevp), &    ! pintm1(k)/(.5*(tm1(k)+tm1(k-1))
         tmp1(plonl), &            ! temporary storage
         dubot(plonl), &           ! lowest layer u change from stress
         dvbot(plonl), &           ! lowest layer v change from stress
         dtbot(plonl), &           ! lowest layer t change from heat flx
         dqbot(plonl,pcnst), &     ! lowest layer q change from const flx
         dshbot(plonl), &          ! lowest layer sh change from wvflx
         thx(plonl,plev), &        ! temperature input + counter gradient
         thv(plonl,plev), &        ! virtual potential temperature
         qmx(plonl,plev,pcnst), &  ! constituents input + counter grad
         shmx(plonl,plev), &       ! sh input + counter grad
         zeh(plonl,plev), &        ! term in tri-diag. matrix system (t & q)
         zem(plonl,plev), &        ! term in tri-diag. matrix system (momentum)
         termh(plonl,plev), &      ! 1./(1.+cah(k) + cch(k) - cch(k)*zeh(k-1))
         termm(plonl,plev)         ! 1./(1.+cam(k) + ccm(k) - ccm(k)*zem(k-1))
    logical :: adjust(plonl)

    real*8 :: &
         um1(plonl,plev), &        ! u wind input
         vm1(plonl,plev), &        ! v wind input
         tm1(plonl,plev), &        ! temperature input
         pmidm1(plonl,plev), &     ! midpoint pressures
         pintm1(plonl,plevp), &    ! interface pressures
         rpdel(plonl,plev), &      ! 1./pdel  (thickness bet interfaces)
         rpdeli(plonl,plev), &     ! 1./pdeli (thickness bet midpoints)
         zm(plonl,plev), &         ! midpoint geoptl height above sfc
         shflx(plonl), &           ! surface sensible heat flux (w/m2)
         cflx(plonl,pcnst), &      ! surface constituent flux (kg/m2/s)
         wvflx(plonl)              ! water vapor flux (kg/m2/s)
    real*8 :: &
         qp1(plonl,plev,pcnst), &  ! moist, tracers after vert. diff
         shp1(plonl,plev), &       ! specific humidity (kg/kg)
         thp(plonl,plev)           ! pot temp after vert. diffusion
    real*8 :: &
         kvh(plonl,plevp), &       ! coefficient for heat and tracers
         kvm(plonl,plevp), &       ! coefficient for momentum
         tpert(plonl), &           ! convective temperature excess
         qpert(plonl), &           ! convective humidity excess
         cgs(plonl,plevp)          ! counter-grad star (cg/flux)

    real*8 :: &
        taux(plonl), &             ! x surface stress (n)
        tauy(plonl), &             ! y surface stress (n)
        ustar(plonl)               ! surface friction velocity

    real*8 :: pblh(plonl)          ! boundary-layer height [m]

    real*8 :: qp0(plonl,plev,pcnst) ! To store initial concentration values
                                    ! (as2)

    !=================================================================
    ! vdiff begins here!
    !=================================================================

    !### Debug
    IF ( LPRT .and. ip < 5 .and. lat < 5 ) &
         CALL DEBUG_MSG( '### VDIFF: vdiff begins' )

    !Populate local variables with values from arguments.(ccc, 11/17/09)
    um1    = uwnd(:,lat,:)
    vm1    = vwnd(:,lat,:)
    tm1    = tadv(:,lat,:)    
    pmidm1 = pmid(:,lat,:)  
    pintm1 = pint(:,lat,:)      
    rpdel  = rpdel_arg(:,lat,:)  
    rpdeli = rpdeli_arg(:,lat,:)
    zm     = zm_arg(:,lat,:)
    shflx  = shflx_arg(:,lat)
    cflx   = sflx(:,lat,:)
    wvflx  = wvflx_arg(:,lat)
    qp1    = as2(:,lat,:,:)
    qp0    = as2(:,lat,:,:)
    shp1   = shp(:,lat,:)
    thp    = thp_arg(:,lat,:)
    kvh    = kvh_arg(:,lat,:)
    kvm    = kvm_arg(:,lat,:)
    tpert  = tpert_arg(:,lat)
    qpert  = qpert_arg(:,lat)
    cgs    = cgs_arg(:,lat,:)
    pblh   = pblh_arg(:,lat)
        
    IF (PRESENT(taux_arg )) taux  = taux_arg(:,lat)
    IF (PRESENT(tauy_arg )) tauy  = tauy_arg(:,lat)
    IF (PRESENT(ustar_arg)) ustar = ustar_arg(:,lat)

!-----------------------------------------------------------------------
! 	... convert the surface fluxes to lowest level tendencies
!-----------------------------------------------------------------------
    rcpair = 1.d0/cpair
    do i = 1,plonl
       tmp1(i)      = ztodt*gravit*rpdel(i,plev)
       ! simplified treatment -- dubot and dvbot are not used under current PBL scheme, anyway
!ccc       if (present(taux) .and. present(tauy)) then 
       if (present(taux_arg) .and. present(tauy_arg)) then 
          dubot(i)     = taux(i)*tmp1(i)
          dvbot(i)     = tauy(i)*tmp1(i)
       endif
       dshbot(i)    = wvflx(i)*tmp1(i)
       dtbot(i)     = shflx(i)*tmp1(i)*rcpair
       kvf(i,plevp) = 0.d0
    end do
    do m = 1,pcnst
       dqbot(:plonl,m) = cflx(:plonl,m)*tmp1(:plonl)
    end do
    
!      !### Debug
    IF ( LPRT .and. ip < 5 .and. lat < 5 ) &
         CALL DEBUG_MSG( '### VDIFF: diffusion begins' )

!-----------------------------------------------------------------------
! 	... set the vertical diffusion coefficient above the top diffusion level
!-----------------------------------------------------------------------
    do k = 1,ntopfl
       kvf(:plonl,k) = 0.d0
    end do

!-----------------------------------------------------------------------
! 	... compute virtual potential temperature for use in static stability 
!           calculation.  0.61 is 1. - r(water vapor)/r(dry air).  use 0.61 instead
!           of a computed variable in order to obtain an identical simulation to
!           case 414.
!-----------------------------------------------------------------------
!      call virtem( thp, shp1, thv, plonl )
    do k = 1,plev
       thv(:,k) = thp(:,k)*(1. + zvir*shp1(:,k))
    end do
    
!      !### Debug
    IF ( LPRT .and. ip < 5 .and. lat < 5 ) &
         CALL DEBUG_MSG( '### VDIFF: compute free atmos. diffusion' )

!-----------------------------------------------------------------------
! 	... compute the free atmosphere vertical diffusion coefficients
!           kvh = kvq = kvm. 
!-----------------------------------------------------------------------
    do k = ntopfl,plev-1
       do i = 1,plonl
!-----------------------------------------------------------------------
! 	... vertical shear squared, min value of (delta v)**2 prevents zero shear.
!-----------------------------------------------------------------------
          dvdz2 = (um1(i,k) - um1(i,k+1))**2 + (vm1(i,k) - vm1(i,k+1))**2
          dvdz2 = max( dvdz2,1.d-36 )
          dz    = zm(i,k) - zm(i,k+1)
          dvdz2 = dvdz2/(dz**2)
!-----------------------------------------------------------------------
! 	... static stability (use virtual potential temperature)
!-----------------------------------------------------------------------
          sstab = gravit*2.d0*(thv(i,k) - thv(i,k+1))/((thv(i,k) + thv(i,k+1))*dz)
!-----------------------------------------------------------------------
! 	... richardson number, stable and unstable modifying functions
!-----------------------------------------------------------------------
          rinub = sstab/dvdz2
          fstab = 1.0d0/(1.0d0 + 10.0d0*rinub*(1.0d0 + 8.0d0*rinub))
          funst = max( 1.d0 - 18.d0*rinub,0.d0 )
!-----------------------------------------------------------------------
! 	... select the appropriate function of the richardson number
!-----------------------------------------------------------------------
          if( rinub < 0.d0 ) then
             fstab = sqrt( funst )
          end if
!-----------------------------------------------------------------------
! 	... neutral diffusion coefficient
!           compute mixing length (z), where z is the interface height estimated
!           with an 8 km scale height.
!-----------------------------------------------------------------------
          kvn = ml2(k)*sqrt( dvdz2 )
!-----------------------------------------------------------------------
! 	... full diffusion coefficient (modified by f(ri)),
!-----------------------------------------------------------------------
          kvf(i,k+1) = max( zkmin,kvn*fstab )
       end do
    end do
    
    !### Debug
    IF ( LPRT .and. ip < 5 .and. lat < 5 ) &
         CALL DEBUG_MSG( '### VDIFF: pbldif begins' )

!-----------------------------------------------------------------------
! 	... determine the boundary layer kvh (=kvq), kvm, 
!           counter gradient terms (cgh, cgq, cgs)
!           boundary layer height (pblh) and 
!           the perturbation temperature and moisture (tpert and qpert)
!           the free atmosphere kv is returned above the boundary layer top.
!-----------------------------------------------------------------------
    ! ustar must always be inputted 
!ccc    if (present(taux) .and. present(tauy)) then
    if (present(taux_arg) .and. present(tauy_arg)) then
       call pbldif( thp, shp1, zm, um1, vm1, &
                    tm1, pmidm1, kvf, cflx, shflx, &
                    kvm, kvh, &
                    cgh, cgq, cgs, pblh, tpert, qpert, &
                    wvflx, cgsh, plonl, & 
                    taux=taux, tauy=tauy, ustar=ustar)
    else
       call pbldif( thp, shp1, zm, um1, vm1, &
                    tm1, pmidm1, kvf, cflx, shflx, &
                    kvm, kvh, &
                    cgh, cgq, cgs, pblh, tpert, qpert, &
                    wvflx, cgsh, plonl, ustar=ustar )
    endif
    
    !### Debug
    IF ( LPRT .and. ip < 5 .and. lat < 5 ) &
         CALL DEBUG_MSG( '### VDIFF: after pbldif' )

!-----------------------------------------------------------------------
! 	... add the counter grad terms to potential temp, specific humidity
!           and other constituents in the bdry layer. note, npbl gives the max
!           num of levels which are permitted to be within the boundary layer.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 	... first set values above boundary layer
!-----------------------------------------------------------------------
    do k = 1,plev-npbl
       do i = 1,plonl
          thx(i,k)  = thp(i,k)
          shmx(i,k) = shp1(i,k)
       end do
       do m = 1,pcnst
          do i = 1,plonl
             qmx(i,k,m) = qp1(i,k,m)
          end do
       end do
    end do
    do k = 2,plev
       do i = 1,plonl
          potbar(i,k) = pintm1(i,k)/(.5d0*(tm1(i,k) + tm1(i,k-1)))
       end do
    end do
    do i = 1,plonl
       potbar(i,plevp) = pintm1(i,plevp)/tm1(i,plev)
    end do

!-----------------------------------------------------------------------
! 	... now focus on the boundary layer
!-----------------------------------------------------------------------
    ztodtgor = ztodt*gravit/rair
    do k = plev-npbl+1,plev
       do i = 1,plonl
          tmp1(i) = ztodtgor*rpdel(i,k)
          thx(i,k) = thp(i,k) + tmp1(i) &
                     *(potbar(i,k+1)*kvh(i,k+1)*cgh(i,k+1) &
                     - potbar(i,k)*kvh(i,k)*cgh(i,k))
          shmx(i,k) = shp1(i,k) + tmp1(i) &
                      *(potbar(i,k+1)*kvh(i,k+1)*cgsh(i,k+1) &
                      - potbar(i,k)*kvh(i,k)*cgsh(i,k))
       end do
       do m = 1,pcnst
          do i = 1,plonl
             qmx(i,k,m) = qp1(i,k,m) + tmp1(i) &
                          *(potbar(i,k+1)*kvh(i,k+1)*cgq(i,k+1,m) &
                          - potbar(i,k)*kvh(i,k)*cgq(i,k,m))
          end do
       end do
    end do

!-----------------------------------------------------------------------
! 	... check for neg qs in each constituent and put the original vertical
!           profile back if a neg value is found. a neg value implies that the
!           quasi-equilibrium conditions assumed for the countergradient term are
!           strongly violated.
!-----------------------------------------------------------------------
    do m = 1,pcnst
       adjust(:plonl) = .false.
       do k = plev-npbl+1,plev
          do i = 1,plonl
             if( qmx(i,k,m) < qmincg(m) ) then
                adjust(i) = .true.
             end if
          end do
       end do
!-----------------------------------------------------------------------
! 	... find long indices of those columns for which negatives were found
!-----------------------------------------------------------------------
       nval = count( adjust(:plonl) )
!-----------------------------------------------------------------------
! 	... replace those columns with original values
!-----------------------------------------------------------------------
       if( nval > 0 ) then
          do k = plev-npbl+1,plev
             where( adjust(:plonl) )
                qmx(:plonl,k,m) = qp1(:plonl,k,m)
             endwhere
          end do
       end if
    end do

!-----------------------------------------------------------------------
! 	... repeat above for sh
!-----------------------------------------------------------------------
    adjust(:plonl) = .false.
    do k = plev-npbl+1,plev
       do i = 1,plonl
!-----------------------------------------------------------------------
! 	... 1.e-12 is the value of qmin (=qmincg) used in ccm2.
!-----------------------------------------------------------------------
          if( shmx(i,k) < 1.d-12 ) then
             adjust(i) = .true.
          end if
       end do
    end do
!-----------------------------------------------------------------------
! 	... find long indices of those columns for which negatives were found
!-----------------------------------------------------------------------
    nval = count( adjust(:plonl) )
!-----------------------------------------------------------------------
! 	... replace those columns with original values
!-----------------------------------------------------------------------
    if( nval > 0 ) then
       do k = plev-npbl+1,plev
          where( adjust(:plonl) )
             shmx(:plonl,k) = shp1(:plonl,k)
          endwhere
       end do
    end if

!-----------------------------------------------------------------------
! 	... determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs
!           of the tridiagonal diffusion matrix. the diagonal elements are a
!           combination of ca and cc; they are not explicitly provided to the 
!           solver
!-----------------------------------------------------------------------
    gorsq = (gravit/rair)**2
    do k = ntopfl,plev-1
       do i = 1,plonl
          tmp2 = ztodt*gorsq*rpdeli(i,k)*(potbar(i,k+1)**2)
          cah(i,k  ) = kvh(i,k+1)*tmp2*rpdel(i,k  )
          cam(i,k  ) = kvm(i,k+1)*tmp2*rpdel(i,k  )
          cch(i,k+1) = kvh(i,k+1)*tmp2*rpdel(i,k+1)
          ccm(i,k+1) = kvm(i,k+1)*tmp2*rpdel(i,k+1)
       end do
    end do
!-----------------------------------------------------------------------
! 	... the last element of the upper diagonal is zero.
!-----------------------------------------------------------------------
    do i = 1,plonl
       cah(i,plev) = 0.d0
       cam(i,plev) = 0.d0
    end do
!-----------------------------------------------------------------------
! 	... calculate e(k) for heat & momentum vertical diffusion.  this term is 
!           required in solution of tridiagonal matrix defined by implicit diffusion eqn.
!-----------------------------------------------------------------------
    do i = 1,plonl
       termh(i,ntopfl) = 1.d0/(1.d0 + cah(i,ntopfl))
       termm(i,ntopfl) = 1.d0/(1.d0 + cam(i,ntopfl))
       zeh(i,ntopfl)   = cah(i,ntopfl)*termh(i,ntopfl)
       zem(i,ntopfl)   = cam(i,ntopfl)*termm(i,ntopfl)
    end do
    do k = ntopfl+1,plev-1
       do i = 1,plonl
          termh(i,k) = 1.d0/(1.d0 + cah(i,k) + cch(i,k) - cch(i,k)*zeh(i,k-1))
          termm(i,k) = 1.d0/(1.d0 + cam(i,k) + ccm(i,k) - ccm(i,k)*zem(i,k-1))
          zeh(i,k)   = cah(i,k)*termh(i,k)
          zem(i,k)   = cam(i,k)*termm(i,k)
       end do
    end do
    
    !### Debug
    IF ( LPRT .and. ip < 5 .and. lat < 5 ) &
         CALL DEBUG_MSG( '### VDIFF: starting diffusion' )

!-----------------------------------------------------------------------
! 	... diffuse constituents
!-----------------------------------------------------------------------
    call qvdiff( pcnst, qmx, dqbot, cch, zeh, &
	         termh, qp1, plonl )

!-----------------------------------------------------------------------
! 	... identify and correct constituents exceeding user defined bounds
!-----------------------------------------------------------------------
!      call qneg3( 'vdiff   ', lat, qp1, plonl )
! just use a simplified treatment
    where (qp1 < 0.d0)
       qp1 = 0.d0
    endwhere

!-----------------------------------------------------------------------
! 	... diffuse sh
!-----------------------------------------------------------------------
    call qvdiff( 1, shmx, dshbot, cch, zeh, &
	         termh, shp1, plonl )

!-----------------------------------------------------------------------
! 	... correct sh
!-----------------------------------------------------------------------
!      call shneg( 'vdiff:sh', lat, shp1, plonl )
! just use a simplified treatment
    where (shp1 < 1.d-12)
       shp1 = 0.d0
    endwhere

!-----------------------------------------------------------------------
! 	... diffuse potential temperature
!-----------------------------------------------------------------------
    call qvdiff( 1, thx, dtbot, cch, zeh, &
	         termh, thp, plonl )

    !Output values from local variables to arguments.(ccc, 11/17/09)
    as2(:,lat,:,:)   = qp1 
    shp(:,lat,:)     = shp1
    thp_arg(:,lat,:) = thp
    kvh_arg(:,lat,:) = kvh    
    kvm_arg(:,lat,:) = kvm    
    tpert_arg(:,lat) = tpert  
    qpert_arg(:,lat) = qpert  
    cgs_arg(:,lat,:)   = cgs    
    pblh_arg(:,lat)  = pblh   
        
    IF (PRESENT(taux_arg )) taux_arg(:,lat)  = taux   
    IF (PRESENT(tauy_arg )) tauy_arg(:,lat)  = tauy   
    IF (PRESENT(ustar_arg)) ustar_arg(:,lat) = ustar  

    !=======================================================
    ! ND15 Diagnostic: 
    ! mass change due to mixing in the boundary layer
    ! ND15 diagnostic moved here to not count the emissions
    ! and dry deposition in the Turbulent Flux.
    ! Needs to call qvdiff with emis+dep = 0 (dqbot) to 
    ! account for all mixing. (ccc, 9/24/10)
    !=======================================================
    IF ( ND15 > 0 ) THEN

       dqbot = 0d0
       call qvdiff( pcnst, qmx, dqbot, cch, zeh, &
            termh, qp1, plonl )

       DO M = 1, pcnst
       DO L = 1, plev 
       do I = 1, plonl
          ! Arrays in vdiff are upside-down
          K = plev - L + 1
          ! qp1 and qp0 are volume mixing ratio
          TURBFLUP(I,lat,k,M) = TURBFLUP(I,lat,k,M) &
                              + (qp1(I,L,M) - qp0(I,L,M)) &
                              * State_Met%AD(I,lat,k) &
                              / ( TCVV(M) * ztodt )
       enddo
       enddo
       ENDDO

    ENDIF

  end subroutine vdiff
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: pbldif
!
! !DESCRIPTION: Subroutine PBLDIF computes the atmospheric boundary layer.
!  The nonlocal scheme determines eddy diffusivities based on a diagnosed 
!  boundary layer height and a turbulent velocity scale. Also, countergradient 
!  effects for heat and moisture, and constituents are included, along with 
!  temperature and humidity perturbations which measure the strength of 
!  convective thermals in the lower part of the atmospheric boundary layer.
!\\
!\\
! References:
!
!  \begin{enumerate}
!  \item Holtslag, A. A. M., and B. A. Boville, 1993: \emph{Local versus 
!         nonlocal boundary-layer diffusion in a global climate model}, 
!         \underline{J. Clim.}, \textbf{6}, 1825-1842.
!  \end{enumerate}
!
! !INTERFACE:
!
  subroutine pbldif( th      ,q       ,z       ,u       ,v, &
                     t       ,pmid    ,kvf     ,cflx    ,shflx, &
                     kvm     ,kvh, &
                     cgh     ,cgq     ,cgs     ,pblh    ,tpert, &
                     qpert   ,wvflx   ,cgsh    ,plonl, &
                     taux    ,tauy    ,ustar)
!
! !USES:
! 
    implicit none
!
! !INPUT PARAMETERS: 
!
    integer, intent(in) :: &
	 plonl
    real*8, intent(in) :: &
         th(plonl,plev), &          ! potential temperature [k]
         q(plonl,plev), &           ! specific humidity [kg/kg]
         z(plonl,plev), &           ! height above surface [m]
         u(plonl,plev), &           ! windspeed x-direction [m/s]
         v(plonl,plev), &           ! windspeed y-direction [m/s]
         t(plonl,plev), &           ! temperature (used for density)
         pmid(plonl,plev), &        ! midpoint pressures
         kvf(plonl,plevp), &        ! free atmospheric eddy diffsvty [m2/s]
         cflx(plonl,pcnst), &       ! surface constituent flux (kg/m2/s)
         wvflx(plonl), &            ! water vapor flux (kg/m2/s)
         shflx(plonl)               ! surface heat flux (w/m2)
!
! !INPUT/OUTPUT PARAMETERS: 
!
    real*8, optional, intent(inout) :: &
         taux(plonl), &            ! x surface stress (n)
         tauy(plonl), &            ! y surface stress (n)
         ustar(plonl)              ! surface friction velocity

    real*8, intent(inout) :: pblh(plonl)       ! boundary-layer height [m]
!
! !OUTPUT PARAMETERS:
!
    real*8, intent(out) :: &
         kvm(plonl,plevp), &        ! eddy diffusivity for momentum [m2/s]
         kvh(plonl,plevp), &        ! eddy diffusivity for heat [m2/s]
         cgh(plonl,plevp), &        ! counter-gradient term for heat [k/m]
         cgq(plonl,plevp,pcnst), &  ! counter-gradient term for constituents
         cgsh(plonl,plevp), &       ! counter-gradient term for sh
         cgs(plonl,plevp), &        ! counter-gradient star (cg/flux)
         tpert(plonl), &            ! convective temperature excess
         qpert(plonl)               ! convective humidity excess
!
! !REVISION HISTORY: 
!  02 Mar 2011 - R. Yantosca - Bug fixes for PGI compiler: these mostly
!                              involve explicitly using "D" exponents
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real*8, parameter :: tiny = 1.d-36      ! lower bound for wind magnitude
    
    integer :: &
         i, &                 ! longitude index
         k, &                 ! level index
         m                    ! constituent index
    logical :: &
         unstbl(plonl), &     ! pts w/unstbl pbl (positive virtual ht flx)
         stblev(plonl), &     ! stable pbl with levels within pbl
         unslev(plonl), &     ! unstbl pbl with levels within pbl
         unssrf(plonl), &     ! unstb pbl w/lvls within srf pbl lyr
         unsout(plonl), &     ! unstb pbl w/lvls in outer pbl lyr
         check(plonl)         ! true=>chk if richardson no.>critcal
    real*8 :: &
         heatv(plonl), &         ! surface virtual heat flux
         thvsrf(plonl), &        ! sfc (bottom) level virtual temperature
         thvref(plonl), &        ! reference level virtual temperature
         tkv, &                  ! model level potential temperature
         therm(plonl), &         ! thermal virtual temperature excess
         phiminv(plonl), &       ! inverse phi function for momentum
         phihinv(plonl), &       ! inverse phi function for heat 
         wm(plonl), &            ! turbulent velocity scale for momentum
         vvk, &                  ! velocity magnitude squared
         zm(plonl), &            ! current level height
         zp(plonl), &            ! current level height + one level up
         khfs(plonl), &          ! surface kinematic heat flux [mk/s]
         kqfs(plonl,pcnst), &    ! sfc kinematic constituent flux [m/s]
         kshfs(plonl), &         ! sfc kinematic moisture flux [m/s]
         zmzp                    ! level height halfway between zm and zp
    real*8 :: &
         rino(plonl,plev), &     ! bulk richardson no. from level to ref lev
         tlv(plonl), &           ! ref. level pot tmp + tmp excess
         fak1(plonl), &          ! k*ustar*pblh
         fak2(plonl), &          ! k*wm*pblh
         fak3(plonl), &          ! fakn*wstr/wm 
         pblk(plonl), &          ! level eddy diffusivity for momentum
         pr(plonl), &            ! prandtl number for eddy diffusivities
         zl(plonl), &            ! zmzp / obukhov length
         zh(plonl), &            ! zmzp / pblh
         zzh(plonl), &           ! (1-(zmzp/pblh))**2
         wstr(plonl), &          ! w*, convective velocity scale
         rrho(plonl), &          ! 1./bottom level density (temporary)
         obklen(plonl), &        ! obukhov length
         ustr                    ! unbounded ustar
    real*8 :: &
         term, &                 ! intermediate calculation
         fac, &                  ! interpolation factor
         pblmin                  ! min pbl height due to mechanical mixing
    
    !=================================================================
    ! pbldif begins here!
    !=================================================================
    
!------------------------------------------------------------------------
! 	... compute kinematic surface fluxes
!------------------------------------------------------------------------
    do i = 1,plonl
       rrho(i)  = rair*t(i,plev)/pmid(i,plev)
       if (present(taux) .and. present(tauy)) then
          ustr     = sqrt( sqrt( taux(i)**2 + tauy(i)**2 )*rrho(i) )
          ustar(i) = max( ustr,.01d0 )
       endif
       khfs(i)  = shflx(i)*rrho(i)/cpair
       kshfs(i) = wvflx(i)*rrho(i)
    end do
    do m = 1,pcnst
       kqfs(:plonl,m) = cflx(:plonl,m)*rrho(:plonl)
    end do

!------------------------------------------------------------------------
! 	... initialize output arrays with free atmosphere values
!------------------------------------------------------------------------
    do k = 1,plevp
       kvm(:,k)  = kvf(:,k)
       kvh(:,k)  = kvf(:,k)
       cgh(:,k)  = 0.d0
       cgsh(:,k) = 0.d0
       cgs(:,k)  = 0.d0
    end do
    do m = 1,pcnst
       do k = 1,plevp
          cgq(:,k,m) = 0.d0
       end do
    end do

!------------------------------------------------------------------------
! 	... compute various arrays for use later:
!------------------------------------------------------------------------
    do i = 1,plonl
       thvsrf(i) = th(i,plev)*(1.0d0 + 0.61d0*q(i,plev))
       heatv(i)  = khfs(i) + 0.61d0*th(i,plev)*kshfs(i)
       wm(i)     = 0.d0
       therm(i)  = 0.d0
       qpert(i)  = 0.d0
       tpert(i)  = 0.d0
       fak3(i)   = 0.d0  
       zh(i)     = 0.d0  
       obklen(i) = -thvsrf(i)*ustar(i)**3 &
                   /(g*vk*(heatv(i) + sign( 1.d-10,heatv(i) )))
    end do
    
    if (pblh_ar) then  ! use archived PBLH
       
       do i = 1,plonl
          if( heatv(i) > 0.d0 ) then
             unstbl(i) = .true.
          else
             unstbl(i) = .false.
          end if
          thvref(i) = th(i,plev)*(1.0d0 + 0.61d0*q(i,plev))
       end do

    else ! use derived PBLH

!------------------------------------------------------------------------
! 	... define first a new factor fac=100 for use in richarson number
!           calculate virtual potential temperature first level
!           and initialize pbl height to z1
!------------------------------------------------------------------------
       fac = 100.
       do i = 1,plonl
          thvref(i) = th(i,plev)*(1.0d0 + 0.61d0*q(i,plev))
          pblh(i)   = z(i,plev)
          check(i)  = .true.
!------------------------------------------------------------------------
! 	... initialization of lowest level ri number 
!           (neglected in initial holtslag implementation)
!------------------------------------------------------------------------
          rino(i,plev) = 0.d0
       end do

!------------------------------------------------------------------------
! 	... pbl height calculation:
!           search for level of pbl. scan upward until the richardson number between
!           the first level and the current level exceeds the "critical" value.
!------------------------------------------------------------------------
       do k = plev-1,plev-npbl+1,-1
          do i = 1,plonl
             if( check(i) ) then
                vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2 + &
                      fac*ustar(i)**2
                vvk = max( vvk,tiny )
                tkv = th(i,k)*(1. + .61d0*q(i,k))
                rino(i,k) = g*(tkv - thvref(i))*(z(i,k)-z(i,plev))/ &
                            (thvref(i)*vvk)
                if( rino(i,k) >= ricr ) then
                   pblh(i) = z(i,k+1) &
                             + (ricr - rino(i,k+1)) &
                             /(rino(i,k) - rino(i,k+1))*(z(i,k) - z(i,k+1))
                   check(i) = .false.
                end if
             end if
          end do
       end do
       
!------------------------------------------------------------------------
! 	... set pbl height to maximum value where computation exceeds number of
!           layers allowed
!------------------------------------------------------------------------
       do i = 1,plonl
          if( check(i) ) then
             pblh(i) = z(i,plevp-npbl)
          end if
       end do
       
!------------------------------------------------------------------------
! 	... improve estimate of pbl height for the unstable points.
!           find unstable points (virtual heat flux is positive):
!------------------------------------------------------------------------
       do i = 1,plonl
          if( heatv(i) > 0.d0 ) then
             unstbl(i) = .true.
             check(i)  = .true.
          else
             unstbl(i) = .false.
             check(i)  = .false.
          end if
       end do

!------------------------------------------------------------------------
! 	... for the unstable case, compute velocity scale and the
!           convective temperature excess:
!------------------------------------------------------------------------
       do i = 1,plonl
          if( check(i) ) then
             phiminv(i)   = (1.d0 - binm*pblh(i)/obklen(i))**onet
             wm(i)        = ustar(i)*phiminv(i)
             therm(i)     = heatv(i)*fak/wm(i)       
             rino(i,plev) = 0.d0
             tlv(i)       = thvref(i) + therm(i)
          end if
       end do

!------------------------------------------------------------------------
! 	... improve pblh estimate for unstable conditions using the
!           convective temperature excess:
!------------------------------------------------------------------------
       do k = plev-1,plev-npbl+1,-1
          do i = 1,plonl
             if( check(i) ) then
                vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2 &
                      + fac*ustar(i)**2
                vvk = max( vvk,tiny )
                tkv = th(i,k)*(1. + 0.61d0*q(i,k))
                rino(i,k) = g*(tkv - tlv(i))*(z(i,k)-z(i,plev)) &
                            /(thvref(i)*vvk)
                if( rino(i,k) >= ricr ) then
                   pblh(i) = z(i,k+1) + (ricr - rino(i,k+1)) &
                             /(rino(i,k) - rino(i,k+1))*(z(i,k) - z(i,k+1))
                   check(i) = .false.
                end if
             end if
          end do
       end do

!------------------------------------------------------------------------
! 	... points for which pblh exceeds number of pbl layers allowed;
!           set to maximum
!------------------------------------------------------------------------
       do i = 1,plonl
          if( check(i) ) then
             pblh(i) = z(i,plevp-npbl)
          end if
       end do

!------------------------------------------------------------------------
! pbl height must be greater than some minimum mechanical mixing depth
! several investigators have proposed minimum mechanical mixing depth
! relationships as a function of the local friction velocity, u*.  we 
! make use of a linear relationship of the form h = c u* where c=700.
! the scaling arguments that give rise to this relationship most often 
! represent the coefficient c as some constant over the local coriolis
! parameter.  here we make use of the experimental results of koracin 
! and berkowicz (1988) [blm, vol 43] for wich they recommend 0.07/f
! where f was evaluated at 39.5 n and 52 n.  thus we use a typical mid
! latitude value for f so that c = 0.07/f = 700.
!------------------------------------------------------------------------
       do i = 1,plonl
          pblmin  = 700.d0*ustar(i)
          pblh(i) = max( pblh(i),pblmin )
       end do
       
    endif ! if pblh_ar


!------------------------------------------------------------------------
! 	... pblh is now available; do preparation for diffusivity calculation:
!------------------------------------------------------------------------
    do i = 1,plonl
       pblk(i) = 0.d0
       fak1(i) = ustar(i)*pblh(i)*vk
!------------------------------------------------------------------------
! 	... do additional preparation for unstable cases only, set temperature
!           and moisture perturbations depending on stability.
!------------------------------------------------------------------------
       if( unstbl(i) ) then
          phiminv(i) = (1.d0 - binm*pblh(i)/obklen(i))**onet
          phihinv(i) = sqrt(1.d0 - binh*pblh(i)/obklen(i))
          wm(i)      = ustar(i)*phiminv(i)
          fak2(i)    = wm(i)*pblh(i)*vk
          wstr(i)    = (heatv(i)*g*pblh(i)/thvref(i))**onet 
          fak3(i)    = fakn*wstr(i)/wm(i)
          tpert(i)   = max( khfs(i)*fak/wm(i),0.d0 )   
          qpert(i)   = max( kshfs(i)*fak/wm(i),0.d0 )    
       else
          tpert(i)   = max( khfs(i)*fak/ustar(i),0.d0 ) 
          qpert(i)   = max( kshfs(i)*fak/ustar(i),0.d0 ) 
       end if
    end do

!------------------------------------------------------------------------
! 	... main level loop to compute the diffusivities and counter-gradient terms
!------------------------------------------------------------------------
    do k = plev,plev-npbl+2,-1
!------------------------------------------------------------------------
! 	... find levels within boundary layer
!------------------------------------------------------------------------
       do i = 1,plonl
          unslev(i) = .false.
          stblev(i) = .false.
          zm(i) = z(i,k)
          zp(i) = z(i,k-1)
          if( zkmin == 0.d0 .and. zp(i) > pblh(i) ) then
             zp(i) = pblh(i)
          end if
          if( zm(i) < pblh(i) ) then
             zmzp = 0.5d0*(zm(i) + zp(i))
             zh(i) = zmzp/pblh(i)
             zl(i) = zmzp/obklen(i)
             if( zh(i) <= 1.d0 ) then
                zzh(i) = (1.d0 - zh(i))**2
             else
                zzh(i) = 0.d0
             end if
!------------------------------------------------------------------------
! 	... stblev for points zm < plbh and stable and neutral
!           unslev for points zm < plbh and unstable
!------------------------------------------------------------------------
             if( unstbl(i) ) then
                unslev(i) = .true.
             else
                stblev(i) = .true.
             end if
          end if
       end do
!------------------------------------------------------------------------
! 	... stable and neutral points; set diffusivities; counter-gradient
!           terms zero for stable case:
!------------------------------------------------------------------------
       do i = 1,plonl
          if( stblev(i) ) then
             if( zl(i) <= 1.d0 ) then
                pblk(i) = fak1(i)*zh(i)*zzh(i)/(1. + betas*zl(i))
             else
                pblk(i) = fak1(i)*zh(i)*zzh(i)/(betas + zl(i))
             end if
             kvm(i,k) = max( pblk(i),kvf(i,k) )
             kvh(i,k) = kvm(i,k)
          end if
       end do
!------------------------------------------------------------------------
! 	... unssrf, unstable within surface layer of pbl
!           unsout, unstable within outer   layer of pbl
!------------------------------------------------------------------------
       do i = 1,plonl
          unssrf(i) = .false.
          unsout(i) = .false.
          if( unslev(i) ) then
             if( zh(i) < sffrac ) then
                unssrf(i) = .true.
             else
                unsout(i) = .true.
             end if
          end if
       end do
!------------------------------------------------------------------------
! 	... unstable for surface layer; counter-gradient terms zero
!------------------------------------------------------------------------
       do i = 1,plonl
          if( unssrf(i) ) then
             term    = (1.d0 - betam*zl(i))**onet
             pblk(i) = fak1(i)*zh(i)*zzh(i)*term
             pr(i)   = term/sqrt(1.d0 - betah*zl(i))
          end if
       end do
!------------------------------------------------------------------------
! 	... unstable for outer layer; counter-gradient terms non-zero:
!------------------------------------------------------------------------
       do i = 1,plonl
          if( unsout(i) ) then
             pblk(i)   = fak2(i)*zh(i)*zzh(i)
             cgs(i,k)  = fak3(i)/(pblh(i)*wm(i))
             cgh(i,k)  = khfs(i)*cgs(i,k)
             pr(i)     = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
             cgsh(i,k) = kshfs(i)*cgs(i,k)
          end if
       end do
       do m = 1,pcnst
          do i = 1,plonl
             if( unsout(i) ) then
                cgq(i,k,m) = kqfs(i,m)*cgs(i,k)
             end if
          end do
       end do
!------------------------------------------------------------------------
! 	... for all unstable layers, set diffusivities
!------------------------------------------------------------------------
       do i = 1,plonl
          if( unslev(i) ) then
             kvm(i,k) = max( pblk(i),kvf(i,k) )
             kvh(i,k) = max( pblk(i)/pr(i),kvf(i,k) )
          end if
       end do

    end do
    
  end subroutine pbldif
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: qvdiff
!
! !DESCRIPTION: Subroutine QVDIFF solve vertical diffusion eqtn for constituent
!  with explicit srfc flux.
!\\
!\\ 
! !INTERFACE:
!
  subroutine qvdiff( ncnst, qm1, qflx, cc, ze, &
	             term, qp1, plonl )
! 
! !USES:
!
    implicit none
!
! !INPUT PARAMETERS: 
!
    integer, intent(in) :: &
         plonl
    integer, intent(in) :: &
         ncnst                    ! num of constituents being diffused
    real*8, intent(in) :: &
         qm1(plonl,plev,ncnst), & ! initial constituent
         qflx(plonl,ncnst), &     ! sfc q flux into lowest model level
         cc(plonl,plev), &        ! -lower diag coeff.of tri-diag matrix
         term(plonl,plev)         ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
!
! !INPUT/OUTPUT PARAMETERS: 
!
    real*8, intent(inout) :: &
         ze(plonl,plev)           ! term in tri-diag. matrix system
!
! !OUTPUT PARAMETERS: 
!
    real*8, intent(out) :: &
         qp1(plonl,plev,ncnst)    ! final constituent
!
! !REMARKS:
!  Procedure for solution of the implicit equation follows :
!  Richtmyer and Morton (1967,pp 198-199)
!
! !REVISION HISTORY: 
!  02 Mar 2011 - R. Yantosca - Bug fixes for PGI compiler: these mostly
!                              involve explicitly using "D" exponents
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real*8 :: &
         zfq(plonl,plev,pcnst), & ! terms appear in soln of tri-diag sys
         tmp1d(plonl)             ! temporary workspace (1d array)
    integer :: &
         i, k, &               ! longitude,vertical indices
         m                     ! constituent index

    !=================================================================
    ! qvdiff begins here!
    !=================================================================
!-----------------------------------------------------------------------
! 	... calculate fq(k).  terms fq(k) and e(k) are required in solution of 
!           tridiagonal matrix defined by implicit diffusion eqn.
!           note that only levels ntopfl through plev need be solved for.
!           no vertical diffusion is applied above this level
!-----------------------------------------------------------------------
    do m = 1,ncnst
       do i = 1,plonl
          zfq(i,ntopfl,m) = qm1(i,ntopfl,m)*term(i,ntopfl)
       end do
       do k = ntopfl+1,plev-1
          do i = 1,plonl
             zfq(i,k,m) = (qm1(i,k,m) + cc(i,k)*zfq(i,k-1,m))*term(i,k)
          end do
       end do
    end do
!-----------------------------------------------------------------------
! 	... bottom level: (includes  surface fluxes)
!-----------------------------------------------------------------------
    do i = 1,plonl
       tmp1d(i) = 1.d0/(1.d0 + cc(i,plev) - cc(i,plev)*ze(i,plev-1))
       ze(i,plev) = 0.
    end do
    do m = 1,ncnst
       do i = 1,plonl
          zfq(i,plev,m) = (qm1(i,plev,m) + qflx(i,m) &
                          + cc(i,plev)*zfq(i,plev-1,m))*tmp1d(i)
       end do
    end do
!-----------------------------------------------------------------------
! 	... perform back substitution
!-----------------------------------------------------------------------
    do m = 1,ncnst
       do i = 1,plonl
          qp1(i,plev,m) = zfq(i,plev,m)
       end do
       do k = plev-1,ntopfl,-1
          do i = 1,plonl
             qp1(i,k,m) = zfq(i,k,m) + ze(i,k)*qp1(i,k+1,m)
          end do
       end do
    end do

  end subroutine qvdiff
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vdiffar
!
! !DESCRIPTION: Subroutine VDIFFAR is the driver routine to compute vertical 
!  diffusion of trace constituents using archived coefficients for cgs and kvh.
!  This is a gutted version of vdiff.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE VDIFFAR( lat   ,tadv , &
                      pmid  ,pint ,rpdel_arg ,rpdeli_arg  ,ztodt, &
                      sflx  ,as2  ,kvh_arg   ,cgs_arg     ,plonl )
!
! !USES:
! 
    implicit none
!
! !INPUT PARAMETERS: 
!
    integer, intent(in) :: lat     ! latitude index
    integer, intent(in) :: plonl   ! lon tile dim
    real*8, intent(in) :: &
         ztodt , &                 ! 2 delta-t
         tadv(:,:,:), &        ! temperature input
         pmid(:,:,:), &     ! midpoint pressures
         pint(:,:,:), &    ! interface pressures
         rpdel_arg(:,:,:), &      ! 1./pdel  (thickness bet interfaces)
         rpdeli_arg(:,:,:), &     ! 1./pdeli (thickness bet midpoints)
         sflx(:,:,:), &      ! surface constituent flux (kg/m2/s)
         kvh_arg(:,:,:), &       ! coefficient for heat and tracers
         cgs_arg(:,:,:)          ! counter-grad star (cg/flux)
!
! !INPUT/OUTPUT PARAMETERS: 
!
    real*8, intent(inout) :: &
         as2(:,:,:,:)     ! moist, tracers after vert. diff
!
! !REVISION HISTORY: 
!  02 Mar 2011 - R. Yantosca - Bug fixes for PGI compiler: these mostly
!                              involve explicitly using "D" exponents
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real*8 :: &
         cah(plonl,plev), &       ! -upper diag for heat and constituts
         cch(plonl,plev), &       ! -lower diag for heat and constits
         cgq(plonl,plevp,pcnst), &! countergrad term for constituent
         potbar(plonl,plevp), &   ! pintm1(k)/(.5*(tm1(k)+tm1(k-1))
         tmp1(plonl), &           ! temporary storage
         tmp2, &                  ! temporary storage
         ztodtgor, &              ! ztodt*gravit/rair
         gorsq, &                 ! (gravit/rair)**2
         dqbot(plonl,pcnst), &    ! lowest layer q change from const flx
         qmx(plonl,plev,pcnst), & ! constituents input + counter grad
         zeh(plonl,plev), &       ! term in tri-diag. matrix system (t & q)
         termh(plonl,plev)        ! 1./(1. + cah(k) + cch(k) - cch(k)*zeh(k-1))
    integer :: &
         indx(plonl), &        ! array of indices of potential q<0
         ilogic(plonl), &      ! 1 => adjust vertical profile
         nval, &               ! num of values which meet criteria
         ii                    ! longitude index of found points
    integer :: &
         i, &                  ! longitude index
         k, &                  ! vertical index
         m                     ! constituent index

    real*8 :: &
         tm1(plonl,plev), &        ! temperature input
         pmidm1(plonl,plev), &     ! midpoint pressures
         pintm1(plonl,plevp), &    ! interface pressures
         rpdel(plonl,plev), &      ! 1./pdel  (thickness bet interfaces)
         rpdeli(plonl,plev), &     ! 1./pdeli (thickness bet midpoints)
         cflx(plonl,pcnst), &      ! surface constituent flux (kg/m2/s)
         kvh(plonl,plevp), &       ! coefficient for heat and tracers
         cgs(plonl,plevp)          ! counter-grad star (cg/flux)
    real*8 :: &
         qp1(plonl,plev,pcnst)     ! moist, tracers after vert. diff
    !=================================================================
    ! vdiffar begins here!
    !=================================================================
    
    !Populate local variables with values from arguments (ccc, 11/17/09)
    tm1    = tadv(:,lat,:)
    pmidm1 = pmid(:,lat,:)
    pintm1 = pint(:,lat,:)
    rpdel  = rpdel_arg(:,lat,:)
    rpdeli = rpdeli_arg(:,lat,:)
    cflx   = sflx(:,lat,:)
    kvh    = kvh_arg(:,lat,:)
    cgs    = cgs_arg(:,lat,:)
    qp1    = as2(:,lat,:,:)
    

!-----------------------------------------------------------------------
! 	... convert the surface fluxes to lowest level tendencies
!-----------------------------------------------------------------------
    do i = 1,plonl
       tmp1(i) = ztodt*gravit*rpdel(i,plev)
    end do
    do m = 1,pcnst
       do i = 1,plonl
          dqbot(i,m) = cflx(i,m)*tmp1(i)
       end do
    end do

!-----------------------------------------------------------------------
! 	... counter gradient terms:
!-----------------------------------------------------------------------
    call pbldifar( tm1, pmidm1, cflx, cgs, cgq, plonl )

!-----------------------------------------------------------------------
! 	... add the counter grad terms to potential temp, specific humidity
!           and other constituents in the bdry layer. note, npbl gives the max
!           num of levels which are permitted to be within the boundary layer.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 	... first set values above boundary layer
!-----------------------------------------------------------------------
    do k = 1,plev-npbl
       do m = 1,pcnst
          qmx(:,k,m) = qp1(:,k,m)
       end do
    end do
    do k = 2,plev
       potbar(:,k) = pintm1(:,k)/(.5d0*(tm1(:,k) + tm1(:,k-1)))
    end do
    potbar(:,plevp) = pintm1(:,plevp)/tm1(:,plev)

!-----------------------------------------------------------------------
! 	... now focus on the boundary layer
!-----------------------------------------------------------------------
    ztodtgor = ztodt*gravit/rair
    do k = plev-npbl+1,plev
       do i = 1,plonl
          tmp1(i) = ztodtgor*rpdel(i,k)
       end do
       do m = 1,pcnst
          do i = 1,plonl
             qmx(i,k,m) = qp1(i,k,m) + tmp1(i)*(potbar(i,k+1)*kvh(i,k+1)* &
                          cgq(i,k+1,m) - potbar(i,k)*kvh(i,k)*cgq(i,k,m))
          end do
       end do
    end do

!-----------------------------------------------------------------------
! 	... check for neg qs in each constituent and put the original vertical
!           profile back if a neg value is found. a neg value implies that the
!           quasi-equilibrium conditions assumed for the countergradient term are
!           strongly violated.
!           original code rewritten by rosinski 7/8/91 to vectorize in longitude.
!-----------------------------------------------------------------------
    do m = 1,pcnst
       ilogic(:plonl) = 0
       do k = plev-npbl+1,plev
          do i = 1,plonl
             if( qmx(i,k,m) < qmincg(m) ) then
                ilogic(i) = 1
             end if
          end do
       end do
!-----------------------------------------------------------------------
! 	... find long indices of those columns for which negatives were found
!-----------------------------------------------------------------------
       nval = count( ilogic(:plonl) == 1 )

!-----------------------------------------------------------------------
! 	... replace those columns with original values
!-----------------------------------------------------------------------
       if( nval > 0 ) then
          do k = plev-npbl+1,plev
             where( ilogic(:plonl) == 1 )
                qmx(:plonl,k,m) = qp1(:plonl,k,m)
             endwhere
          end do
       end if
    end do

!-----------------------------------------------------------------------
! 	... determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs
!           of the tridiagonal diffusion matrix. the diagonal elements are a
!           combination of ca and cc; they are not explicitly provided to the 
!           solver
!-----------------------------------------------------------------------
    gorsq = (gravit/rair)**2
    do k = ntopfl,plev-1
       do i = 1,plonl
          tmp2 = ztodt*gorsq*rpdeli(i,k)*(potbar(i,k+1)**2)
          cah(i,k  ) = kvh(i,k+1)*tmp2*rpdel(i,k  )
          cch(i,k+1) = kvh(i,k+1)*tmp2*rpdel(i,k+1)
       end do
    end do
!-----------------------------------------------------------------------
! 	... the last element of the upper diagonal is zero.
!-----------------------------------------------------------------------
    do i = 1,plonl
       cah(i,plev) = 0.d0
    end do
!-----------------------------------------------------------------------
! 	... calculate e(k) for heat vertical diffusion.  this term is 
!           required in solution of tridiagonal matrix defined by implicit 
!           diffusion eqn.
!-----------------------------------------------------------------------
    do i = 1,plonl
       termh(i,ntopfl) = 1.d0/(1.d0 + cah(i,ntopfl))
       zeh(i,ntopfl) = cah(i,ntopfl)*termh(i,ntopfl)
    end do
    do k = ntopfl+1,plev-1
       do i = 1,plonl
          termh(i,k) = 1.d0/(1.d0 + cah(i,k) + cch(i,k) - cch(i,k)*zeh(i,k-1))
          zeh(i,k) = cah(i,k)*termh(i,k)
       end do
    end do
!-----------------------------------------------------------------------
! 	... diffuse constituents
!-----------------------------------------------------------------------
    call qvdiff( pcnst, qmx, dqbot, cch, zeh, &
	         termh, qp1, plonl )
!-----------------------------------------------------------------------
! 	... identify and correct constituents exceeding user defined bounds
!-----------------------------------------------------------------------
!      call qneg3( 'vdiff   ', lat, qp1(1,1,1), plonl )
!     simplified treatment
    where (qp1 < 0.d0)
       qp1 = 0.d0
    endwhere
    
    !Output values from local variables to arguments.(ccc, 11/17/09)
    as2(:,lat,:,:) = qp1
    
  END SUBROUTINE VDIFFAR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: pbldifar
!
! !DESCRIPTION: Subroutine PBLDIFAR is a modified version of pbldif which only
!  calculates cgq given cgs.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PBLDIFAR( t, pmid, cflx, cgs, cgq, plonl )
! 
! !USES:
!
    implicit none
!
! !INPUT PARAMETERS: 
!
    integer, intent(in) :: &
         plonl
    real*8, intent(in) :: &
         t(plonl,plev), &        ! temperature (used for density)
         pmid(plonl,plev), &     ! midpoint pressures
         cflx(plonl,pcnst), &     ! surface constituent flux (kg/m2/s)
         cgs(plonl,plevp)        ! counter-gradient star (cg/flux)
!
! !OUTPUT PARAMETERS: 
!
    real*8, intent(out) :: &
         cgq(plonl,plevp,pcnst)  ! counter-gradient term for constituents
!
! !REVISION HISTORY: 
!  02 Mar 2011 - R. Yantosca - Bug fixes for PGI compiler: these mostly
!                              involve explicitly using "D" exponents
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: &
         i, &                 ! longitude index
         k, &                 ! level index
         m                    ! constituent index
    real*8 :: &
         rrho(plonl), &          ! 1./bottom level density
         kqfs(plonl,pcnst)       ! sfc kinematic constituent flux [m/s]

    !=================================================================
    ! pbldifar begins here!
    !=================================================================

!------------------------------------------------------------------------
! 	... compute kinematic surface fluxes
!------------------------------------------------------------------------
    rrho(:) = rair*t(:,plev)/pmid(:,plev)
    do m = 1,pcnst
       kqfs(:,m) = cflx(:,m)*rrho(:)
    end do
!------------------------------------------------------------------------
! 	... initialize output arrays with free atmosphere values
!------------------------------------------------------------------------
    do m = 1,pcnst
       do k = 1,plevp
          cgq(:,k,m) = 0.d0
       end do
    end do
!------------------------------------------------------------------------
! 	... compute the counter-gradient terms:
!------------------------------------------------------------------------
    do k = plev,plev-npbl+2,-1
       do m = 1,pcnst
          cgq(:,k,m) = kqfs(:,m)*cgs(:,k)
       end do
    end do

  END SUBROUTINE PBLDIFAR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vdinti
!
! !DESCRIPTION: Subroutine VDINTI initializes time independent fields for 
!  vertical diffusion. Calls initialization routine for boundary layer scheme.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE VDINTI
!
! !USES:
! 
    USE PRESSURE_MOD, ONLY : GET_AP, GET_BP
    USE ERROR_MOD,    ONLY : ALLOC_ERR
    
    implicit none
!
! !REVISION HISTORY: 
!  02 Mar 2011 - R. Yantosca - Bug fixes for PGI compiler: these mostly
!                              involve explicitly using "D" exponents
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real*8, parameter :: pbl_press = 400.e2         ! pressure cap for pbl (pa)
    integer :: k, &                               ! vertical loop index
               m
    
    integer :: AS
    
!#if defined( DEVEL )
    real*8, allocatable :: ref_pmid(:)
!#else
!    real*8 :: ref_pmid(LLPAR)
!#endif

    !=================================================================
    ! vdinti begins here!
    !=================================================================

!#if defined( DEVEL )
    ALLOCATE( ref_pmid(LLPAR), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'ref_pmid' )
    ref_pmid = 0.d0
    plevp = plev+1
!#endif
!-----------------------------------------------------------------------
! 	... hard-wired numbers.
!           zkmin = minimum k = kneutral*f(ri)
!-----------------------------------------------------------------------
    zkmin = .01d0

!-----------------------------------------------------------------------
! 	... set physical constants for vertical diffusion and pbl
!-----------------------------------------------------------------------

    ! REF_PMID is indexed with K=1 being the atm. top and K=PLEV being 
    ! the surface.  Eliminate call to UPSIDEDOWN (bmy, 12/21/10)
    do k = 1, plev
       ref_pmid(plev-k+1) = 0.5d0*( GET_AP(k  )*100.d0 + GET_BP(k  )*1.d5 + &
                                    GET_AP(k+1)*100.d0 + GET_BP(k+1)*1.d5 )
    enddo

!-----------------------------------------------------------------------
! 	... derived constants
!           ntopfl = top level to which v-diff is applied
!           npbl   = max number of levels (from bottom) in pbl
!-----------------------------------------------------------------------
    do k = plev,1,-1
       if( ref_pmid(k) < pbl_press ) then
          exit
       end if
    end do
    npbl = max( 1,plev - k )
    write(*,*) 'vdinti: pbl height will be limited to bottom ',npbl, &
               ' model levels. top is ',1.d-2*ref_pmid(plevp-npbl),' hpa'
    if( plev == 1 ) then
       ntopfl = 0
    else
       ntopfl = 1
    end if

!-----------------------------------------------------------------------
! 	... set the square of the mixing lengths
!-----------------------------------------------------------------------
!#if defined( DEVEL )
    ALLOCATE( ml2(plevp), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'ml2' )
!#endif

    ml2(1) = 0.d0
    do k = 2,plev
       ml2(k) = (30.d0)**2
    end do
    ml2(plevp) = 0.d0
!-----------------------------------------------------------------------
! 	... set the minimum mixing ratio for the counter-gradient term.
!           normally this should be the same as qmin.
!-----------------------------------------------------------------------
    
    ALLOCATE( qmincg(pcnst), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'DETRAINE' )
    
    do m = 1,pcnst
       qmincg(m) = 0.d0
    end do

!-----------------------------------------------------------------------
! 	... initialize pbl variables
!-----------------------------------------------------------------------
    call pbinti( gravit)

  END SUBROUTINE VDINTI
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vdiffdr
!
! !DESCRIPTION: Subroutine VDIFFDR calculates the vertical diffusion on a 
!  latitude slice of data.
!
!  \begin{enumerate}
!  \item The dummy argument as2 is in v\/v. (lin, 06/04/08)
!  \item TCVV and TRACER\_MW\_KG assume 12 g/mol for all HCs. Thus, when using 
!         them to convert units of HCs to be the inputs for vdiffdr, the
!         converted units are NOT kg/kg for concentrations and kg/m2/s for 
!         surface flux. However, since the units for both inputs are 
!         consistent, there should not be any problem. (lin, 06/04/08)
!  \end{enumerate}
!
! !INTERFACE:
!
  SUBROUTINE VDIFFDR( as2, State_Met )
!
! !USES:
! 
    USE COMODE_MOD,         ONLY : JLOP,      REMIS,   VOLUME
    USE DAO_MOD,            ONLY : USTAR
    USE DAO_MOD,            ONLY : IS_ICE, IS_LAND
    USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_DD, ADD_HgP_DD
    USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_SNOWPACK
    USE DIAG_MOD,           ONLY : AD44
    USE DRYDEP_MOD,         ONLY : DEPNAME, NUMDEP, NTRAIND, DEPSAV
    USE DRYDEP_MOD,         ONLY : DRYHg0, DRYHg2, DRYHgP !cdh
    USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GRID_MOD,           ONLY : GET_AREA_M2
    USE LOGICAL_MOD,        ONLY : LDYNOCEAN, LGTMM !cdh
    USE LOGICAL_MOD,        ONLY : LSOILNOX
    USE OCEAN_MERCURY_MOD,  ONLY : Fp, Fg !hma
    USE OCEAN_MERCURY_MOD,  ONLY : LHg2HalfAerosol !cdh
    USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_m, COMPUTE_PBL_HEIGHT, &
                                   GET_PBL_MAX_L, GET_FRAC_UNDER_PBLTOP
    USE PRESSURE_MOD,       ONLY : GET_PEDGE, GET_PCENTER
    USE TIME_MOD,           ONLY : GET_TS_CONV, GET_TS_EMIS
    USE TRACER_MOD,         ONLY : N_TRACERS,  TRACER_MW_KG, TCVV, &
                                   ID_EMITTED, TRACER_COEFF, TRACER_COEFF, &
                                   TRACER_NAME
    USE TRACER_MOD,         ONLY : ITS_A_TAGOX_SIM, ITS_A_TAGCO_SIM
    USE TRACER_MOD,         ONLY : ITS_A_CH4_SIM
    USE TRACER_MOD,         ONLY : ITS_A_MERCURY_SIM ! (cdh 8/28/09)
    USE TRACER_MOD,         ONLY : ITS_A_FULLCHEM_SIM  !bmy
    USE TRACERID_MOD,       ONLY : IS_Hg0, IS_Hg2, IS_HgP
    USE VDIFF_PRE_MOD,      ONLY : IIPAR, JJPAR, IDEMS, NEMIS, NCS, ND44, &
                                   NDRYDEP, emis_save

#   include "define.h"

    implicit none
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Meteorology State object
    TYPE(MetState), INTENT(INOUT)         :: State_Met   

    ! Advected species
    REAL*8,         intent(inout), TARGET :: as2(IIPAR,JJPAR,LLPAR,N_TRACERS) 
!
! !REMARKS:
!  Need to declare the Meteorology State object (State_MET) with
!  INTENT(INOUT).  This is because VDIFF will modify the specific
!  humidity field. (bmy, 11/21/12)
!
! !REVISION HISTORY:
! (1 ) Calls to vdiff and vdiffar are now done with full arrays as arguments.
!       (ccc, 11/19/09)
!  04 Jun 2010 - C. Carouge  - Updates for mercury simulations with GTMM 
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  24 Sep 2010 - J. Lin      - Move ND15 to vdiff.  
!  21 Dec 2010 - R. Yantosca - Add logical flags for different sim types
!  21 Dec 2010 - R. Yantosca - Now call ITS_A_FULLCHEM_SIM instead of
!                              relying on NCS == 0
!  22 Dec 2010 - C. Carouge  - Combine array flipping w/ unit conversion 
!                              to save on operations
!  02 Mar 2011 - R. Yantosca - Bug fixes for PGI compiler: these mostly
!                              involve explicitly using "D" exponents
!  26 Apr 2011 - J. Fisher   - Use MERRA land fraction information
!  25 Oct 2011 - H. Amos     - bring Hg2 gas-particle partitioning code into
!                              v9-01-02
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.2 in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  22 Jun 2012 - R. Yantosca - Now use pointers to flip arrays in vertical
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: I,J,L,JLOOP,N,NN

!    REAL*8                :: SNOW_HT !cdh - obsolete
    REAL*8                :: FRAC_NO_HG0_DEP !jaf 
    LOGICAL               :: ZERO_HG0_DEP !jaf 

    real*8, TARGET, dimension(IIPAR,JJPAR,LLPAR) :: pmid, rpdel, rpdeli, zm
    real*8, TARGET, dimension(IIPAR,JJPAR,LLPAR+1) :: pint
    real*8, TARGET, dimension(IIPAR,JJPAR,N_TRACERS) :: sflx
    real*8, TARGET, dimension(IIPAR,JJPAR,N_TRACERS) :: eflx, dflx ! surface flux
    real*8, TARGET, dimension(IIPAR,JJPAR,LLPAR+1) :: cgs, kvh, kvm
    real*8, TARGET, dimension(IIPAR,JJPAR) :: pblh, tpert, qpert
    real*8, TARGET, dimension(IIPAR,JJPAR,LLPAR) :: thp         ! potential temperature
    real*8, TARGET, dimension(IIPAR,JJPAR) :: shflx    ! water vapor flux
    real*8, TARGET, dimension(IIPAR,JJPAR,LLPAR) :: t1
    real*8, TARGET, dimension(IIPAR,JJPAR,LLPAR,N_TRACERS) :: as ! save tracer MR 
                                                         ! before vdiffdr
    real*8 :: vtemp
    real*8 :: p0 = 1.d5
    real*8 :: dtime
    real*8 :: wk1, wk2
    real*8 :: soilflux
    integer :: pbl_top
      
    REAL*8  :: DEP_KG !(cdh, 8/28/09)

    ! Array to store a single level of the AS2 array,
    ! so as not to blow up the parallelization (ccc, 12/22.10)
    REAL*8, dimension(IIPAR, JJPAR, N_TRACERS)  :: as2_scal

    ! Add flags
    LOGICAL :: IS_CH4, IS_FULLCHEM, IS_Hg, IS_TAGOx, IS_TAGCO

    ! Pointers 
    REAL*8,  POINTER :: p_um1   (:,:,:  )
    REAL*8,  POINTER :: p_vm1   (:,:,:  )
    REAL*8,  POINTER :: p_tadv  (:,:,:  )
    REAL*8,  POINTER :: p_hflux (:,:    )
    REAL*8,  POINTER :: p_ustar (:,:    )
    REAL*8,  POINTER :: p_pmid  (:,:,:  )
    REAL*8,  POINTER :: p_pint  (:,:,:  )
    REAL*8,  POINTER :: p_rpdel (:,:,:  ) 
    REAL*8,  POINTER :: p_rpdeli(:,:,:  )
    REAL*8,  POINTER :: p_zm    (:,:,:  )
    REAL*8,  POINTER :: p_thp   (:,:,:  )
    REAL*8,  POINTER :: p_kvh   (:,:,:  )
    REAL*8,  POINTER :: p_kvm   (:,:,:  )
    REAL*8,  POINTER :: p_cgs   (:,:,:  )
    REAL*8,  POINTER :: p_shp   (:,:,:  )
    REAL*8,  POINTER :: p_t1    (:,:,:  )
    REAL*8,  POINTER :: p_as2   (:,:,:,:)

    !=================================================================
    ! vdiffdr begins here!
    !=================================================================

    !### Debug
    IF ( LPRT ) CALL DEBUG_MSG( '### VDIFFDR: VDIFFDR begins' )
    
    ! Initialize local arrays. (ccc, 12/21/10)
    pmid    = 0d0
    rpdel   = 0d0
    rpdeli  = 0d0
    zm      = 0d0
    pint    = 0d0
    sflx    = 0d0
    eflx    = 0d0
    dflx    = 0d0
    soilflux = 0d0
    cgs     = 0d0
    kvh     = 0d0
    kvm     = 0d0
    pblh    = 0d0
    tpert   = 0d0
    qpert   = 0d0
    thp     = 0d0
    shflx   = 0d0
    t1      = 0d0
    as2_scal= 0d0

    ! Test for different types of simulations and save in local variables/
    ! These are used in the parallel DO loops below (bmy, 12/21/10)
    IS_CH4      = ITS_A_CH4_SIM()
    IS_FULLCHEM = ITS_A_FULLCHEM_SIM()
    IS_Hg       = ITS_A_MERCURY_SIM()
    IS_TAGCO    = ITS_A_TAGCO_SIM()
    IS_TAGOX    = ITS_A_TAGOX_SIM()

    dtime = GET_TS_CONV()*60d0 ! min -> second
    
    shflx = State_Met%EFLUX / latvap ! latent heat -> water vapor flux

! (Turn off parallelization for now, skim 6/20/12)
    
!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L )
    do J = 1, JJPAR
    do I = 1, IIPAR

    ! calculate variables related to pressure
    do L = 1, LLPAR
       pmid(I,J,L) = GET_PCENTER(I,J,L)*100.d0 ! hPa -> Pa
       pint(I,J,L) = GET_PEDGE(I,J,L)*100.d0   ! hPa -> Pa
       ! calculate potential temperature
       thp(I,J,L) = State_Met%T(I,J,L)*(p0/pmid(I,J,L))**cappa
    enddo
    pint(I,J,LLPAR+1) = GET_PEDGE(I,J,LLPAR+1)
    
    enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L )
    do J = 1, JJPAR
    do I = 1, IIPAR
    do L = 1, LLPAR
    ! Corrected calculation of zm.
    ! Use temperature instead of virtual temperature to be consistent with 
    ! the calculation of BXHEIGHT. (lin, 06/02/08)
    !zm(I,J,L) = sum(BXHEIGHT(I,J,1:L))
       zm(I,J,L) = sum( State_Met%BXHEIGHT(I,J,1:L)) &
                 - log( pmid(I,J,L)/pint(I,J,L+1) )  &
                 * r_g * State_Met%T(I,J,L)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    ! Have to separate the calculations of pmid/pint and rpdel/rpdeli.
    ! (Lin, 06/02/08)
!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L )
    do J = 1, JJPAR
    do I = 1, IIPAR
       do L = 1, LLPAR
          rpdel(I,J,L) = 1.d0 / (pint(I,J,L) - pint(I,J,L+1))
       enddo

       !rpdeli(I,J,1) = 1.d0 / (PS(I,J) - pmid(I,J,1))
       rpdeli(I,J,1) = 0.d0 ! follow mozart setup (shown in mo_physlic.F90) 

       do L = 2, LLPAR
          rpdeli(I,J,L) = 1.d0 / (pmid(I,J,L-1) - pmid(I,J,L))
       enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !!! calculate surface flux = emissions - dry deposition !!!

    ! Define slice of AS2, so as not to blow up the parallelization
    ! (ccc, bmy, 12/20/10)
    as2_scal = as2(:,:,1,:)

!$OMP PARALLEL DO       &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( I, J, L, N, NN, JLOOP, wk1, wk2, pbl_top, DEP_KG ) &
!!$OMP PRIVATE( SNOW_HT ) &
!$OMP PRIVATE( FRAC_NO_HG0_DEP, ZERO_HG0_DEP )
    do J = 1, JJPAR
    do I = 1, IIPAR

       !----------------------------------------------------------------
       ! Add emissions for full-chemistry simulation
       !----------------------------------------------------------------
       IF ( IS_FULLCHEM .and. NCS > 0 ) THEN

          do N = 1, NEMIS(NCS)
             NN = IDEMS(N)

             ! for emissions in the lowest model layer only
             IF ( NN > 0 ) THEN
                JLOOP = JLOP(I,J,1)
                eflx(I,J,NN) = REMIS(JLOOP,N) * TRACER_MW_KG(NN)
             ENDIF

          enddo

          ! additional step to convert from molec spec/cm3/s to kg/m2/s
          eflx(I,J,:) = eflx(I,J,:) * State_Met%BXHEIGHT(I,J,1) / &
                        6.022d23    * 1.d6

          !add the tracer coef. (i.e., one ISOP molecule has five carbon atoms)
          ! (lin, 06/07/08) 
          do N = 1, N_TRACERS
             if ( ID_EMITTED(N) .le. 0 ) cycle
             eflx(I,J,N) = eflx(I,J,N) * TRACER_COEFF(N,ID_EMITTED(N))
          enddo

          ! add surface emis of aerosols 
          ! (after converting kg/box/timestep to kg/m2/s)
          ! Should NOT use ID_EMITTED here, since it is only for gases 
          ! for SMVGEAR. (Lin, 06/10/08)
          do N = 1, N_TRACERS
             eflx(I,J,N) = eflx(I,J,N) + emis_save(I,J,N)       &
                                       / GET_AREA_M2( I, J, 1 ) &
                                       / GET_TS_EMIS() / 60.d0
          enddo

       ENDIF

       !----------------------------------------------------------------
       ! Zero emissions for tagged CO simulation
       !
       ! CO emis are considered in tagged_co_mod.f.  This over-
       ! simplified treatment may be inconsistent with the full 
       ! chemistry simulation. Hopefully this simplification wouldn't 
       ! cause too much problem, since the std. tagged_co simulation 
       ! is also approximate, anyway. (Lin, 06/20/09) 
       !----------------------------------------------------------------
       IF ( IS_TAGCO ) THEN
          eflx(I,J,:) = 0d0 
       ENDIF

       !----------------------------------------------------------------
       ! Add emissions for offline CH4 simulation
       !----------------------------------------------------------------
       IF ( IS_CH4 ) THEN
          ! add surface emis
          ! (after converting kg/box/timestep to kg/m2/s)
          ! Should NOT use ID_EMITTED here, since it is only for gases 
          ! for SMVGEAR. (Lin, 06/10/08)
          do N = 1, N_TRACERS
             eflx(I,J,N) = eflx(I,J,N) + emis_save(I,J,N) &
                         / GET_AREA_M2( I, J, 1 )         &
                         / GET_TS_EMIS() / 60.d0
          enddo
       endif

       !----------------------------------------------------------------
       ! Add emissions for offline mercury simulation
       !----------------------------------------------------------------
       IF ( IS_Hg ) THEN
          do N = 1, N_TRACERS
             eflx(I,J,N) = eflx(I,J,N) + emis_save(I,J,N) & 
                         / GET_AREA_M2( I, J, 1 )         &
                         / GET_TS_EMIS() / 60.d0
          enddo
       ENDIF

       !----------------------------------------------------------------
       ! Apply dry deposition frequencies
       !----------------------------------------------------------------
       do N = 1, NUMDEP ! NUMDEP includes all gases/aerosols
          ! Now include sea salt dry deposition (jaegle 5/11/11)
          IF (TRIM( DEPNAME(N) ) == 'DST1'.OR. &
              TRIM( DEPNAME(N) ) == 'DST2'.OR. &
              TRIM( DEPNAME(N) ) == 'DST3'.OR. &
              TRIM( DEPNAME(N) ) == 'DST4') CYCLE
              !TRIM( DEPNAME(N) ) == 'SALA'.OR. &
              !TRIM( DEPNAME(N) ) == 'SALC') CYCLE

          ! gases + aerosols for full chemistry 
          NN   = NTRAIND(N)
          if (NN == 0) CYCLE

          ! adding the backward consistency with previous GEOS-Chem drydep 
          ! calculation. (Lin, 06/04/2008) 
          ! given that as2 is in v/v
          !dflx(I,J,NN) = DEPSAV(I,J,N) * as2(I,J,1,NN) / TCVV(NN) 
          
          ! use mean concentration within the PBL for calculating drydep 
          ! fluxes
          if (pbl_mean_drydep) then 
             wk1 = 0.d0
             wk2 = 0.d0
             pbl_top = GET_PBL_MAX_L() ! the highest layer the PBL reaches, 
             ! globally
             do L = 1, pbl_top
                wk1 = wk1 + as2(I,J,L,NN) * State_Met%AD(I,J,L)* &
                      GET_FRAC_UNDER_PBLTOP(I,J,L)
                wk2 = wk2 + State_Met%AD(I,J,L) * &
                      GET_FRAC_UNDER_PBLTOP(I,J,L)
             enddo
             ! since we only use the ratio of wk1 / wk2, there should not be
             ! a problem even if the PBL top is lower than the top of the 
             ! first (lowest) model layer
             ! given that as2 is in v/v
             dflx(I,J,NN) = DEPSAV(I,J,N) * (wk1/(wk2+1.d-30)) / TCVV(NN)
             
             ! consistency with the standard GEOS-Chem setup (Lin, 07/14/08)
             if (drydep_back_cons) then 
                dflx(I,J,NN) = dflx(I,J,NN) * (wk2+1.d-30) / &
                               State_Met%AD(I,J,1)         * &
                               State_Met%BXHEIGHT(I,J,1)   / &
                               GET_PBL_TOP_m(I,J)
             endif
          else

             ! only use the lowest model layer for calculating drydep fluxes
             ! given that as2 is in v/v
             ! NOTE: Now use as2_scal(I,J,NN), instead of as2(I,J,1,NN) to 
             ! avoid seg faults in parallelization (ccarouge, bmy, 12/20/10)
             dflx(I,J,NN) = DEPSAV(I,J,N) * as2_scal(I,J,NN) / TCVV(NN)

             !------------------------------------------------------------------
             !Prior to 25 Oct 2011, H Amos
             !! If flag is set to treat Hg2 as half aerosol, half gas, then
             !! use average deposition velocity (cdh, 9/01/09)
             !IF ( LHG2HALFAEROSOL .AND. IS_HG2(NN) ) THEN
             !
             !   ! NOTE: Now use as2_scal(I,J,NN), instead of as2(I,J,1,NN) to 
             !   ! avoid seg faults in parallelization (ccarouge, bmy, 12/20/10)
             !   dflx(I,J,NN) =  &
             !        ( DEPSAV(I,J,DRYHg2) +  DEPSAV(I,J,DRYHgP) ) / 2D0 * &
             !        as2_scal(I,J,NN) / TCVV(NN) 
             !   
             !ENDIF
             !
             IF ( IS_HG2(NN) ) THEN

                IF ( LHG2HALFAEROSOL ) THEN
                   ! NOTE: Now use as2_scal(I,J,NN), instead of as2(I,J,1,NN) to 
                   ! avoid seg faults in parallelization (ccarouge, bmy, 12/20/10)

                   ! partition Hg2 50/50 gas/particle
                   dflx(I,J,NN) =  &
                        ( DEPSAV(I,J,DRYHg2) +  DEPSAV(I,J,DRYHgP) ) / 2D0 * &
                        as2_scal(I,J,NN) / TCVV(NN) 
                ELSE
                   
                   ! temperature-dependent Hg2 partitioning
                   dflx(I,J,NN) = ( DEPSAV(I,J,DRYHg2)*Fg(I,J,1) + &
                                    DEPSAV(I,J,DRYHgP)*Fp(I,J,1) ) * &
                                    as2_scal(I,J,NN) / TCVV(NN) 
                ENDIF
                   
             ENDIF
             !------------------------------------------------------------------
          endif
          
!--jaf.start          
          ! Restructure to use fractional land cover info available in
          ! MERRA (jaf,4/26/11)
!          ! Hg(0) exchange with the ocean is handled by ocean_mercury_mod
!          ! so disable deposition over water here.
!          ! LWI is exactly zero for water (includes some fresh water) 
!          ! in GEOS-5. Note LWI is not exactly zero for earlier GEOS versions,
!          ! but non-local PBL mixing is only defined for GEOS-5.
!          ! ocean_mercury_mod defines ocean based on fraction 
!          ! land, albedo and mixed layer depth. The difference with LWI is
!          ! small. (cdh, 8/28/09) 
!          IF ( IS_Hg .AND. IS_HG0(NN) .AND. LWI(I,J) == 0 ) THEN
!             DFLX(I,J,NN) = 0D0
!          ENDIF
!
!          ! CDH (9/11/09)
!          ! Turn off Hg(0) deposition to snow and ice because we haven't yet
!          ! included emission from these surfaces and most field studies
!          ! suggest Hg(0) emissions exceed deposition during sunlit hours.
!#if   defined( GEOS_5 ) || defined( MERRA )
!          ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
!          SNOW_HT = SNOMAS(I,J)
!#else
!          ! GEOS1-4 snow heigt (water equivalent) in mm
!          SNOW_HT = SNOW(I,J)
!#endif 
!
!          IF ( IS_Hg .AND. IS_HG0(NN) .AND. &
!               ( IS_ICE(I,J) .OR. (IS_LAND(I,J) .AND. SNOW_HT > 10d0) ) ) THEN
!             DFLX(I,J,NN) = 0D0
!          ENDIF

          ! Hg(0) exchange with the ocean is handled by ocean_mercury_mod
          ! so disable deposition over water here.
          ! Turn off Hg(0) deposition to snow and ice because we haven't yet
          ! included emission from these surfaces and most field studies
          ! suggest Hg(0) emissions exceed deposition during sunlit hours.

          ! Except in MERRA, we assume entire grid box is water or ice
          ! if conditions are met (jaf, 4/26/11)
          FRAC_NO_HG0_DEP = 1d0

#if   defined( MERRA ) || defined( GEOS_57 )
          FRAC_NO_HG0_DEP = MIN( State_Met%FROCEAN(I,J) + &
                                 State_Met%FRSNO(I,J)   + &
                                 State_Met%FRLANDIC(I,J), 1d0)
          ZERO_HG0_DEP    = ( FRAC_NO_HG0_DEP > 0d0 )

#elif defined( GEOS_5 )
          ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
          ZERO_HG0_DEP = (( State_Met%LWI(I,J) == 0      )  .OR.  &
                          ( IS_ICE ( I, J, State_Met     )) .OR.  &
                          ( IS_LAND( I, J, State_Met     )  .AND. &
                            State_Met%SNOMAS(I,J) > 10d0 ))

#else
          ! GEOS1-4 snow heigt (water equivalent) in mm
          ZERO_HG0_DEP = (( State_Met%LWI(I,J) == 0      )  .OR.  &
                          ( IS_ICE ( I, J, State_Met     )) .OR.  &
                          ( IS_LAND( I, J, State_Met     )  .AND. &
                            State_Met%SNOW(I,J)   > 10d0 ))
#endif
          
          IF ( IS_Hg .AND. IS_HG0(NN) ) THEN
             IF ( ZERO_HG0_DEP ) THEN
                DFLX(I,J,NN) = DFLX(I,J,NN) * &
                               MAX(1d0 - FRAC_NO_HG0_DEP,0d0)
             ENDIF
          ENDIF

!--jaf.end

       enddo

       !----------------------------------------------------------------
       ! Apply dry deposition frequencies for Tagged Ox simulation
       ! (Jintai Lin, 06/21/08)
       !----------------------------------------------------------------
       IF ( IS_TAGOX ) THEN
          do N = 2, N_TRACERS ! the first species, Ox, has been done above
             if (pbl_mean_drydep) then
                wk1 = 0.d0
                wk2 = 0.d0
                pbl_top = GET_PBL_MAX_L() ! the highest layer the PBL reaches,
                                          ! globally
                do L = 1, pbl_top
                   wk1 = wk1 + as2(I,J,L,N) * State_Met%AD(I,J,L) * &
                               GET_FRAC_UNDER_PBLTOP(I,J,L)
                   wk2 = wk2 + State_Met%AD(I,J,L) * &
                               GET_FRAC_UNDER_PBLTOP(I,J,L)
                enddo
                ! since we only use the ratio of wk1 / wk2, there should not be
                ! a problem even if the PBL top is lower than the top of the 
                ! first (lowest) model layer
                ! given that as2 is in v/v
                dflx(I,J,N) = DEPSAV(I,J,1) * (wk1/(wk2+1.d-30)) / TCVV(1)

                ! Consistent with the standard GEOS-Chem setup.(Lin, 07/14/08) 
                if (drydep_back_cons) then 
                   dflx(I,J,N) = dflx(I,J,N) * (wk2+1.d-30) / &
                                 State_Met%AD(I,J,1)        * &
                                 State_Met%BXHEIGHT(I,J,1)  / &
                                 GET_PBL_TOP_m(I,J)
                endif
             else 
                ! only use the lowest model layer for calculating drydep fluxes
                ! given that as2 is in v/v
                ! NOTE: Now use as2_scal(I,J,NN), instead of as2(I,J,1,NN) to 
                ! avoid seg faults in parallelization (ccarouge, bmy, 12/20/10)
                dflx(I,J,N) = DEPSAV(I,J,1) * as2_scal(I,J,N) / TCVV(1) 
             endif
          enddo
       endif

       ! virtual temperature in the lowest model layer
       ! vtemp = tadv(I,J,1)*(1. + zvir*shp(I,J,1))
       ! for deposition: additional step to convert from s-1 to kg/m2/s
       ! dflx(I,J,:) = dflx(I,J,:) * pmid(I,J,1) / rair / vtemp * BXHEIGHT(I,J,1)
       ! alternate method to convert from s-1 to kg/m2/s
       dflx(I,J,:) = dflx(I,J,:) * State_Met%AD(I,J,1) / &
                     GET_AREA_M2( I, J, 1 ) 

       ! surface flux = emissions - dry deposition
       sflx(I,J,:) = eflx(I,J,:) - dflx(I,J,:) ! kg/m2/s


       !----------------------------------------------------------------
       ! Archive Hg deposition for surface reservoirs (cdh, 08/28/09)
       !----------------------------------------------------------------
       IF ( IS_Hg ) THEN
          
          ! Loop over # of drydep species
          DO N = 1, NUMDEP
             
             ! GEOS_Chem tracer number
             NN = NTRAIND(N)
             
             ! Deposition mass, kg
             DEP_KG = dflx( I, J, NN ) * GET_AREA_M2( I, J, 1 ) &
                    * GET_TS_CONV()    * 60d0

             IF ( IS_Hg2(NN) ) THEN 
                
                CALL ADD_HG2_DD( I, J, NN, DEP_KG )
                CALL ADD_Hg2_SNOWPACK( I, J, NN, DEP_KG, State_Met )

             ELSE IF ( IS_HgP( NN ) ) THEN
                
                CALL ADD_HGP_DD( I, J, NN, DEP_KG )
                CALL ADD_Hg2_SNOWPACK( I, J, NN, DEP_KG, State_Met )

             ENDIF

          ENDDO
       ENDIF

    enddo
    enddo
!$OMP END PARALLEL DO

    ! drydep fluxes diag. for SMVGEAR mechanism 
    ! for gases -- moved from DRYFLX in drydep_mod.f to here
    ! for aerosols -- 
    if (ND44 > 0 .or. LGTMM .or. LSOILNOX) then

       do N = 1, NUMDEP
          SELECT CASE ( DEPNAME(N) )
             ! non gases + aerosols for fully chemistry 
             !CASE ( 'DST1', 'DST2', 'DST3', 'DST4', 'SALA', &
             !       'SALC' )
	     ! now include sea salt dry deposition (jaegle 5/11/11)
             CASE ( 'DST1', 'DST2', 'DST3', 'DST4')
                CYCLE
             CASE DEFAULT
                ! Locate position of each tracer in DEPSAV
                NN   = NTRAIND(N)
                if (NN == 0) CYCLE
                ! only for the lowest model layer
                ! Convert : kg/m2/s -> molec/cm2/s
                ! consider timestep difference between convection and emissions
		IF(ND44 > 0 .or. LGTMM) THEN                
		AD44(:,:,N,1) = AD44(:,:,N,1) + dflx(:,:,NN) &
                                / TRACER_MW_KG(NN) * 6.022d23 * 1.d-4 &
                                * GET_TS_CONV() / GET_TS_EMIS() 
		ENDIF

		IF(LSOILNOX) THEN
			soilflux = 0d0
			DO J = 1, JJPAR
		 	DO I = 1, IIPAR			
				soilflux = dflx(I,J,NN) &
		                        / TRACER_MW_KG(NN) * 6.022d23 * 1.d-4 &
		                        * GET_TS_CONV() / GET_TS_EMIS() 

		      	        CALL SOIL_DRYDEP ( I, J, 1, NN, soilflux)
			ENDDO
			ENDDO

		ENDIF

          END SELECT
       enddo

       ! Add ITS_A_TAGOX_SIM (Lin, 06/21/08)
       IF ( IS_TAGOX ) THEN
          ! The first species, Ox, has been done above
          do N = 2, N_TRACERS 
             ! Convert : kg/m2/s -> molec/cm2/s
             ! Consider timestep difference between convection and emissions
             AD44(:,:,N,1) = AD44(:,:,N,1) + dflx(:,:,N) &
                             / TRACER_MW_KG(1) * 6.022d23 * 1.d-4 &
                             * GET_TS_CONV() / GET_TS_EMIS()
             AD44(:,:,N,2) = AD44(:,:,1,2) ! drydep velocity
          enddo
       endif


    endif

	!Maasa, Add SoilNOx deposition to allow SN code to work with NLPBL on.

    !### Debug
    IF ( LPRT ) CALL DEBUG_MSG( '### VDIFFDR: after emis. and depdrp' )

    if( divdiff ) then
      
       if ( pblh_ar ) then
       do J = 1, JJPAR
       do I = 1, IIPAR
         pblh(I,J) = GET_PBL_TOP_m(I,J) ! obtain archived PBLH
       enddo
       enddo
       endif

       !-------------------------------------------------------------------
       ! Now use pointers to flip arrays in the vertical (bmy, 6/22/15)
       !-------------------------------------------------------------------

       ! 3-D fields on level centers
       p_um1              => State_Met%U    ( :, :, LLPAR  :1:-1    )   
       p_vm1              => State_Met%V    ( :, :, LLPAR  :1:-1    )
       p_tadv             => State_Met%T    ( :, :, LLPAR  :1:-1    )
       p_hflux            => State_Met%HFLUX
       p_ustar            => State_Met%USTAR
       p_pmid             => pmid           ( :, :, LLPAR  :1:-1    )
       p_rpdel            => rpdel          ( :, :, LLPAR  :1:-1    )
       p_rpdeli           => rpdeli         ( :, :, LLPAR  :1:-1    )
       p_zm               => zm             ( :, :, LLPAR  :1:-1    )
       p_thp              => thp            ( :, :, LLPAR  :1:-1    )
!#if defined( DEVEL )
       p_shp              => State_Met%SPHU ( :, :, LLPAR  :1:-1    )
!!#else
!       p_shp              => shp           ( :, :, LLPAR  :1:-1    )
!#endif
       ! 3-D fields on level edges
       p_pint             => pint           ( :, :, LLPAR+1:1:-1    )
       p_kvh              => kvh            ( :, :, LLPAR+1:1:-1    )
       p_kvm              => kvm            ( :, :, LLPAR+1:1:-1    )
       p_cgs              => cgs            ( :, :, LLPAR+1:1:-1    )

       ! Tracer concentration fields
       p_as2              => as2   ( :, :, LLPAR  :1:-1, : )

       ! Convert v/v -> m/m (i.e., kg/kg)
       DO N = 1, N_TRACERS
          p_as2(:,:,:,N)  =  p_as2(:,:,:,N) / TCVV(N) 
       ENDDO

       ! Convert g/kg -> kg/kg
       p_shp              =  p_shp * 1.d-3 

       !### Debug
       IF ( LPRT ) CALL DEBUG_MSG( '### VDIFFDR: before vdiff' )

!$OMP PARALLEL DO DEFAULT( SHARED )      &
!$OMP PRIVATE( J )     
       do J = 1, JJPAR
          call vdiff( J,        1,      p_um1,  p_vm1,           &
                      p_tadv,   p_pmid, p_pint, p_rpdel,         &
                      p_rpdeli, dtime,  p_zm,   p_hflux,         &
                      sflx,     p_thp,  p_as2,  pblh,            &
                      p_kvh,    p_kvm,  tpert,  qpert,           &
                      p_cgs,    p_shp,  shflx,  IIPAR,           &
                      State_Met,        ustar_arg=ustar )
       enddo
!$OMP END PARALLEL DO

       !### Debug
       IF ( LPRT ) CALL DEBUG_MSG( '### VDIFFDR: after vdiff' )

       ! Convert kg/kg -> v/v
       DO N = 1, N_TRACERS
          p_as2(:,:,:,N) = p_as2(:,:,:,N) * TCVV(N)
       ENDDO

       ! Convert kg/kg -> g/kg
       p_shp    = p_shp * 1.d3

       ! Free pointers
       NULLIFY( p_um1,   p_vm1,    p_tadv, p_pmid, p_pint )
       NULLIFY( p_rpdel, p_rpdeli, p_zm,   p_thp,  p_cgs  )
       NULLIFY( p_kvh,   p_kvm,    p_shp,  p_as2,  p_hflux)

    else if( arvdiff ) then
!-----------------------------------------------------------------------
!  	... vertical diffusion using archived values of cgs and kvh.
!
!       %%% NOTE: THIS SECTION IS NORMALLY NOT EXECUTED %%%
!       %%% BECAUSE ARVDIFF IS SET TO .FALSE. ABOVE     %%% 
!-----------------------------------------------------------------------

       !-------------------------------------------------------------------
       ! Now use pointers to flip arrays in the vertical (bmy, 6/22/15)
       !-------------------------------------------------------------------

       ! INPUTS: 3-D fields on level centers
       p_tadv   => State_Met%T( :, :, LLPAR  :1:-1   )
       p_pmid   => pmid       ( :, :, LLPAR  :1:-1   )
       p_rpdel  => rpdel      ( :, :, LLPAR  :1:-1   )
       p_rpdeli => rpdeli     ( :, :, LLPAR  :1:-1   )

       ! INPUTS: 3-D fields on level edges
       p_pint   => pint       ( :, :, LLPAR+1:1:-1   )
       p_kvh    => kvh        ( :, :, LLPAR+1:1:-1   )
       p_cgs    => cgs        ( :, :, LLPAR+1:1:-1   )

       ! INPUTS: Tracer concentration fields
       p_as2    => as2        ( :, :, LLPAR:1:-1,  : )

       ! Convert from v/v -> m/m (i.e., kg/kg)
       do N = 1, N_TRACERS
          p_as2(:,:,:,N) = p_as2(:,:,:,N) / TCVV(N) 
       enddo

       !### Debug
       IF ( LPRT ) CALL DEBUG_MSG( '### VDIFFDR: before vdiffar' )

!!$OMP PARALLEL DO DEFAULT( SHARED )   &
!!$OMP PRIVATE( J )
       do J = 1, JJPAR
          call vdiffar( J,      p_tadv, p_pmid, p_pint, p_rpdel, p_rpdeli,  &
                        dtime,  sflx,   p_as2,  p_kvh,  p_cgs,   IIPAR     )
      enddo
!!$OMP END PARALLEL DO

       !### Debug
       IF ( LPRT ) CALL DEBUG_MSG( '### VDIFFDR: after vdiffar' )

       ! Convert from m/m (i.e. kg/kg) -> v/v
       do N = 1, N_TRACERS
          p_as2(:,:,:,N) = p_as2(:,:,:,N) * TCVV(N) 
       enddo

       ! Free pointers
       NULLIFY( p_tadv, p_pmid, p_rpdel, p_rpdeli )
       NULLIFY( p_pint, p_kvh,  p_cgs,   p_as2    )

    end if

    !-------------------------------------------------------------------
    ! re-compute PBL variables wrt derived pblh (in m)
    !-------------------------------------------------------------------
    if (.not. pblh_ar) then

       ! PBL is in m 
       State_Met%PBLH = pblh 

       CALL COMPUTE_PBL_HEIGHT( State_Met )
    endif

!      !### Debug
    IF ( LPRT ) CALL DEBUG_MSG( '### VDIFFDR: VDIFFDR finished' )

  END SUBROUTINE VDIFFDR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_pbl_mix_2
!
! !DESCRIPTION: Subroutine DO\_PBL\_MIX\_2 is the driver routine for planetary 
!  boundary layer mixing. The PBL layer height and related quantities are 
!  always computed.   Mixing of tracers underneath the PBL top is toggled 
!  by the DO\_TURBDAY switch. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_PBL_MIX_2( DO_TURBDAY, State_Met )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE LOGICAL_MOD,        ONLY : LTURB, LPRT
    USE PBL_MIX_MOD,        ONLY : INIT_PBL_MIX, COMPUTE_PBL_HEIGHT
    USE TIME_MOD,           ONLY : ITS_TIME_FOR_EMIS
    USE TRACER_MOD,         ONLY : N_TRACERS, STT, TCVV, ITS_A_FULLCHEM_SIM
    USE VDIFF_PRE_MOD,      ONLY : EMISRR, EMISRRN

    IMPLICIT NONE
#   include "define.h"
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN) :: DO_TURBDAY  ! Switch which turns on PBL
                                              !  mixing of tracers
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY: 
!  11 Feb 2005 - R. Yantosca - Initial version
!  21 Dec 2010 - R. Yantosca - Now only call SETEMIS for fullchem simulations
!  22 Dec 2010 - R. Yantosca - Bug fix: print debug output only if LPRT=T
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE :: FIRST = .TRUE.

    !=================================================================
    ! DO_PBL_MIX2 begins here!
    !=================================================================
    !call flush(6)
    
    ! First-time initialization
    IF ( FIRST ) THEN
       CALL INIT_PBL_MIX
       call vdinti
       FIRST = .FALSE.
    ENDIF

    ! Compute PBL height and related quantities
    CALL COMPUTE_PBL_HEIGHT( State_Met )

    !=================================================================
    ! For full-chemistry simulations, call routine SETEMIS
    ! which sets up the emission rates array REMIS
    !=================================================================
    IF ( ITS_A_FULLCHEM_SIM() ) THEN

       ! If it's time to do emissions, call SETEMIS
       IF ( ITS_TIME_FOR_EMIS() ) THEN 
          CALL SETEMIS( EMISRR, EMISRRN, .TRUE., State_Met )
          IF ( LPRT ) CALL DEBUG_MSG( '### DO_PBL_MIX_2: aft SETEMIS' )
       ENDIF

    ENDIF

    ! Do mixing of tracers in the PBL (if necessary)
    IF ( DO_TURBDAY ) THEN 
       CALL VDIFFDR( STT, State_Met )
       IF( LPRT ) CALL DEBUG_MSG( '### DO_PBL_MIX_2: after VDIFFDR' )
    ENDIF


  END SUBROUTINE DO_PBL_MIX_2
!EOC  
END MODULE vdiff_mod

