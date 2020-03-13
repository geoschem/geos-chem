!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: vdiff_mod.F90
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
  USE ERROR_MOD,     ONLY : DEBUG_MSG              ! Routine for debug output
  USE PhysConstants                                ! Physical constants
  USE PRECISION_MOD                                ! For GEOS-Chem Precision(fp)
  USE VDIFF_PRE_MOD, ONLY : PCNST                  ! N_TRACERS
  USE VDIFF_PRE_MOD, ONLY : prtDebug               ! Debug print?
  USE VDIFF_PRE_MOD, ONLY : LTURB                  ! Do PBL mixing?

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: DO_PBL_MIX_2
!
! !PRIVATE DATA MEMBERS:
!
  save

  integer :: plev
  integer :: plevp

  real(fp), parameter ::            &
       cpair  = 1004.64e+0_fp,      &
       latvap = 2.5104e+06_fp,      &
       rhoh2o = 1.e+3_fp,           &
       tfh2o  = 273.16e+0_fp,       &
       rair   = Rd,                 &
       rh2o   = Rv,                 &
       gravit = g0,                 &
       zvir   = rh2o/rair - 1.,     &
       cappa  = Rd/cpair,           &
       r_g    = Rd / g0

  ! Obsolete variables and variables that are now defined with global params
  ! (ewl, 1/7/16)
!  real(fp), parameter ::          &
!       rearth = 6.37122e+6_fp,      & ! not used
!       cpwv   = 1.81e+3_fp,         & ! not used
!       ra     = 1.e+0_fp/rearth,    & ! not used
!       epsilo = 0.622e+0_fp,        & ! not used
!       latice = 3.336e+5_fp,        & ! not used
!       rair   = 287.04e+0_fp,       & ! now use global Rd
!       rh2o   = 461.e+0_fp,         & ! now use global Rv
!       gravit = 9.80616e+0_fp,      & ! now use global g0
!       cappa  = rair/cpair,     &     ! now use global Rd
!       r_g    = rair / gravit,  &     ! now use global Rd and g0

!-----------------------------------------------------------------------
! 	... pbl constants
!-----------------------------------------------------------------------

  ! These are constants, so use PARAMETER tag
  real(fp), parameter ::   &
       betam  = 15.e+0_fp,   & ! constant in wind gradient expression
       betas  =  5.e+0_fp,   & ! constant in surface layer gradient expression
       betah  = 15.e+0_fp,   & ! constant in temperature gradient expression
       fak    =  8.5e+0_fp,  & ! constant in surface temperature excess
       fakn   =  7.2e+0_fp,  & ! constant in turbulent prandtl number
       ricr   =   .3e+0_fp,  & ! critical richardson number
       sffrac =   .1e+0_fp,  & ! surface layer fraction of boundary layer
       vk     =   .4e+0_fp     ! von karmans constant

  ! These are assigned later, so we can't use the PARAMETER tag
  real(fp) ::              &
       g,                & ! gravitational acceleration
       onet,             & ! 1/3 power in wind gradient expression
       ccon,             & ! fak * sffrac * vk
       binm,             & ! betam * sffrac
       binh                ! betah * sffrac

!-----------------------------------------------------------------------
! 	... constants used in vertical diffusion and pbl
!-----------------------------------------------------------------------
  real(fp) :: &
       zkmin            ! minimum kneutral*f(ri)
  real(fp), allocatable :: ml2(:)   ! mixing lengths squaredB
  real(fp), allocatable :: qmincg(:)   ! min. constituent concentration
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

!----------------------------------------------------------------------
! Diagnostic quantities
!-----------------------------------------------------------------------

  LOGICAL :: Do_ND44           = .FALSE. ! Do ND44 bpch diagnostic
!
! !REMARKS:
!  The non-local PBL mixing routine VDIFF modifies the specific humidity,
!  (State_Met%SPHU) field.  Therefore, we must pass State_Met as an argument
!  to DO_PBL_MIX_2 and VDIFFDR with INTENT(INOUT).
!                                                                             .
!  Because logical_mod.F and tracer_mod.F have been superseded by Input_Opt,
!  we now use VDIFF_PRE_MOD to supply values
!
! !REVISION HISTORY:
!  (1 ) This code is modified from mo_vdiff.F90 in MOZART-2.4. (lin, 5/14/09)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
contains
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    real(fp), intent(in) :: gravx     !  acceleration of gravity
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
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
    onet  = 1e+0_fp/3.e+0_fp

!-----------------------------------------------------------------------
! 	... derived constants
!-----------------------------------------------------------------------
    ccon = fak*sffrac*vk
    binm = betam*sffrac
    binh = betah*sffrac

  end subroutine pbinti
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
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
  subroutine vdiff( lat,        ip,        uwnd,       vwnd,                 &
                    tadv,       pmid,      pint,       rpdel_arg,            &
                    rpdeli_arg, ztodt,     zm_arg,     shflx_arg,            &
                    sflx,       thp_arg,   as2,        pblh_arg,             &
                    kvh_arg,    kvm_arg,   tpert_arg,  qpert_arg,            &
                    cgs_arg,    shp,       wvflx_arg,  plonl,                &
                    Input_Opt,  State_Met, State_Chm,  State_Diag,           &
                    taux_arg,   tauy_arg,  ustar_arg,  RC                   )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : IS_SAFE_DIV
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Met_Mod,      ONLY : MetState

    implicit none
!
! !INPUT PARAMETERS:
!
    integer, intent(in) :: lat, ip ! latitude index, long tile index
    integer, intent(in) :: plonl   ! number of local longitudes
    real(fp),  intent(in) ::   &
         ztodt                     ! 2 delta-t
    real(fp),  intent(in) ::   &
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
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    real(fp), intent(inout) :: &
         as2(:,:,:,:),       &     ! moist, tracers after vert. diff
         shp(:,:,:),         &     ! specific humidity (kg/kg)
         thp_arg(:,:,:)            ! pot temp after vert. diffusion
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object

    real(fp), optional, intent(inout) :: &
         taux_arg(:,:),      &     ! x surface stress (n)
         tauy_arg(:,:),      &     ! y surface stress (n)
         ustar_arg(:,:)            ! surface friction velocity

    real(fp), intent(inout) :: pblh_arg(:,:) ! boundary-layer height [m]

!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: RC
    real(fp), intent(out) :: &
         kvh_arg(:,:,:),     &     ! coefficient for heat and tracers
         kvm_arg(:,:,:),     &     ! coefficient for momentum
         tpert_arg(:,:),     &     ! convective temperature excess
         qpert_arg(:,:),     &     ! convective humidity excess
         cgs_arg(:,:,:)            ! counter-grad star (cg/flux)
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    real(fp) :: &
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
    real(fp) :: &
         cah(plonl,plev), &        ! -upper diag for heat and constituts
         cam(plonl,plev), &        ! -upper diagonal for momentum
         cch(plonl,plev), &        ! -lower diag for heat and constits
         ccm(plonl,plev), &        ! -lower diagonal for momentum
         cgh(plonl,plevp), &       ! countergradient term for heat
         cgq(plonl,plevp,pcnst), & ! countergrad term for constituent
         cgsh(plonl,plevp), &      ! countergrad term for sh
         kvf(plonl,plevp)          ! free atmosphere kv at interfaces
    real(fp) :: &
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

    real(fp) :: &
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
    real(fp) :: &
         qp1(plonl,plev,pcnst), &  ! moist, tracers after vert. diff
         shp1(plonl,plev), &       ! specific humidity (kg/kg)
         thp(plonl,plev)           ! pot temp after vert. diffusion
    real(fp) :: &
         kvh(plonl,plevp), &       ! coefficient for heat and tracers
         kvm(plonl,plevp), &       ! coefficient for momentum
         tpert(plonl), &           ! convective temperature excess
         qpert(plonl), &           ! convective humidity excess
         cgs(plonl,plevp)          ! counter-grad star (cg/flux)

    real(fp) :: &
        taux(plonl), &             ! x surface stress (n)
        tauy(plonl), &             ! y surface stress (n)
        ustar(plonl)               ! surface friction velocity

    real(fp) :: pblh(plonl)          ! boundary-layer height [m]

    real(fp) :: qp0(plonl,plev,pcnst) ! To store initial concentration values
                                    ! (as2)

    real(fp) :: sum_qp0, sum_qp1    ! Jintai Lin 20180809

    !=================================================================
    ! vdiff begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Initialize
    sum_qp0  = 0.0_fp
    sum_qp1  = 0.0_fp

    !### Debug
    IF ( prtDebug .and. ip < 5 .and. lat < 5 ) &
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
    rcpair = 1.e+0_fp/cpair
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
       kvf(i,plevp) = 0.e+0_fp
    end do
    do m = 1,pcnst
       dqbot(:plonl,m) = cflx(:plonl,m)*tmp1(:plonl)
    end do

!      !### Debug
    IF ( prtDebug .and. ip < 5 .and. lat < 5 ) &
         CALL DEBUG_MSG( '### VDIFF: diffusion begins' )

!-----------------------------------------------------------------------
! 	... set the vertical diffusion coefficient above the top diffusion level
!-----------------------------------------------------------------------
    do k = 1,ntopfl
       kvf(:plonl,k) = 0.e+0_fp
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
    IF ( prtDebug .and. ip < 5 .and. lat < 5 ) &
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
          dvdz2 = (um1(i,k) - um1(i,k+1))**2 + &
                 (vm1(i,k) - vm1(i,k+1))**2
          dvdz2 = max( dvdz2,1.e-36_fp )
          dz    = zm(i,k) - zm(i,k+1)
          dvdz2 = dvdz2/(dz**2)
!-----------------------------------------------------------------------
! 	... static stability (use virtual potential temperature)
!-----------------------------------------------------------------------
          sstab = gravit*2.e+0_fp*(thv(i,k) - thv(i,k+1))/((thv(i,k) &
                 + thv(i,k+1))*dz)
!-----------------------------------------------------------------------
! 	... richardson number, stable and unstable modifying functions
!-----------------------------------------------------------------------
          rinub = sstab/dvdz2
          fstab = 1.0e+0_fp/(1.0e+0_fp + 10.0e+0_fp*rinub*(1.0e+0_fp &
                 + 8.0e+0_fp*rinub))
          funst = max( 1.e+0_fp - 18.e+0_fp*rinub,0.e+0_fp )
!-----------------------------------------------------------------------
! 	... select the appropriate function of the richardson number
!-----------------------------------------------------------------------
          if( rinub < 0.e+0_fp ) then
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
    IF ( prtDebug .and. ip < 5 .and. lat < 5 ) &
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
                    taux=taux, tauy=tauy, ustar=ustar )
    else
       call pbldif( thp, shp1, zm, um1, vm1, &
                    tm1, pmidm1, kvf, cflx, shflx, &
                    kvm, kvh, &
                    cgh, cgq, cgs, pblh, tpert, qpert, &
                    wvflx, cgsh, plonl, ustar=ustar )
    endif

    !### Debug
    IF ( prtDebug .and. ip < 5 .and. lat < 5 ) &
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
          potbar(i,k) = pintm1(i,k)/(.5e+0_fp*(tm1(i,k) + tm1(i,k-1)))
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
       cah(i,plev) = 0.e+0_fp
       cam(i,plev) = 0.e+0_fp
    end do

!-----------------------------------------------------------------------
! 	... calculate e(k) for heat & momentum vertical diffusion.  this term is
!           required in solution of tridiagonal matrix defined by implicit diffusion eqn.
!-----------------------------------------------------------------------
    do i = 1,plonl
       termh(i,ntopfl) = 1.e+0_fp/(1.e+0_fp + cah(i,ntopfl))
       termm(i,ntopfl) = 1.e+0_fp/(1.e+0_fp + cam(i,ntopfl))
       zeh(i,ntopfl)   = cah(i,ntopfl)*termh(i,ntopfl)
       zem(i,ntopfl)   = cam(i,ntopfl)*termm(i,ntopfl)
    end do
    do k = ntopfl+1,plev-1
       do i = 1,plonl
          termh(i,k) = 1.e+0_fp/(1.e+0_fp + cah(i,k) + cch(i,k) &
                      - cch(i,k)*zeh(i,k-1))
          termm(i,k) = 1.e+0_fp/(1.e+0_fp + cam(i,k) + ccm(i,k) &
                      - ccm(i,k)*zem(i,k-1))
          zeh(i,k)   = cah(i,k)*termh(i,k)
          zem(i,k)   = cam(i,k)*termm(i,k)
       end do
    end do

    !### Debug
    IF ( prtDebug .and. ip < 5 .and. lat < 5 ) &
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
    where (qp1 < 0.e+0_fp)
       qp1 = 0.e+0_fp
    endwhere

!-----------------------------------------------------------------------
! Simple bug fix to ensure mass conservation - Jintai Lin 20180809
!   Without this fix, mass is almost but not completed conserved,
!   which is OK for full chemistry simulations but a big problem
!   for long lived species such as CH4 and CO2
!-----------------------------------------------------------------------
    DO M = 1, pcnst
    do I = 1, plonl

       ! total mass in the PBL (ignoring the v/v -> m/m conversion)
       !   including pre-mixing mass and surface flux (emis+drydep)
       sum_qp0 = sum(qp0(I,ntopfl:plev,M) * &
                 State_Met%AD(I,lat,plev-ntopfl+1:1:-1)) &
               + dqbot(I,M) * State_Met%AD(I,lat,1)

       ! total mass in the PBL (ignoring the v/v -> m/m conversion)
       sum_qp1 = sum(qp1(I,ntopfl:plev,M) * &
                 State_Met%AD(I,lat,plev-ntopfl+1:1:-1))

       IF ( IS_SAFE_DIV( sum_qp0, sum_qp1 ) ) THEN
          qp1(I,ntopfl:plev,M) = qp1(I,ntopfl:plev,M) * &
                                 sum_qp0 / sum_qp1
       ENDIF

    enddo
    ENDDO

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
       shp1 = 0.e+0_fp
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

  end subroutine vdiff
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
                     taux    ,tauy    ,ustar )
!
! !USES:
!
    implicit none
!
! !INPUT PARAMETERS:
!
    integer, intent(in) :: &
	 plonl
    real(fp), intent(in) :: &
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
    real(fp), optional, intent(inout) :: &
         taux(plonl), &            ! x surface stress (n)
         tauy(plonl), &            ! y surface stress (n)
         ustar(plonl)              ! surface friction velocity

    real(fp), intent(inout) :: pblh(plonl)       ! boundary-layer height [m]
!
! !OUTPUT PARAMETERS:
!
    real(fp), intent(out) :: &
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp), parameter :: tiny = 1.e-36_fp      ! lower bound for wind magnitude

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
    real(fp) :: &
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
    real(fp) :: &
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
    real(fp) :: &
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
          ustar(i) = max( ustr,.01e+0_fp )
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
       cgh(:,k)  = 0.e+0_fp
       cgsh(:,k) = 0.e+0_fp
       cgs(:,k)  = 0.e+0_fp
    end do
    do m = 1,pcnst
       do k = 1,plevp
          cgq(:,k,m) = 0.e+0_fp
       end do
    end do

!------------------------------------------------------------------------
! 	... compute various arrays for use later:
!------------------------------------------------------------------------
    do i = 1,plonl
       thvsrf(i) = th(i,plev)*(1.0e+0_fp + 0.61e+0_fp*q(i,plev))
       heatv(i)  = khfs(i) + 0.61e+0_fp*th(i,plev)*kshfs(i)
       wm(i)     = 0.e+0_fp
       therm(i)  = 0.e+0_fp
       qpert(i)  = 0.e+0_fp
       tpert(i)  = 0.e+0_fp
       fak3(i)   = 0.e+0_fp
       zh(i)     = 0.e+0_fp
       obklen(i) = -thvsrf(i)*ustar(i)**3 &
                   /(g*vk*(heatv(i) + sign( 1.d-10,heatv(i) )))
    end do

    if (pblh_ar) then  ! use archived PBLH

       do i = 1,plonl
          if( heatv(i) > 0.e+0_fp ) then
             unstbl(i) = .true.
          else
             unstbl(i) = .false.
          end if
          thvref(i) = th(i,plev)*(1.0e+0_fp + 0.61e+0_fp*q(i,plev))
       end do

    else ! use derived PBLH

!------------------------------------------------------------------------
! 	... define first a new factor fac=100 for use in richarson number
!           calculate virtual potential temperature first level
!           and initialize pbl height to z1
!------------------------------------------------------------------------
       fac = 100.
       do i = 1,plonl
          thvref(i) = th(i,plev)*(1.0e+0_fp + 0.61e+0_fp*q(i,plev))
          pblh(i)   = z(i,plev)
          check(i)  = .true.
!------------------------------------------------------------------------
! 	... initialization of lowest level ri number
!           (neglected in initial holtslag implementation)
!------------------------------------------------------------------------
          rino(i,plev) = 0.e+0_fp
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
                tkv = th(i,k)*(1. + .61e+0_fp*q(i,k))
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
          if( heatv(i) > 0.e+0_fp ) then
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
             phiminv(i)   = (1.e+0_fp - binm*pblh(i)/obklen(i))**onet
             wm(i)        = ustar(i)*phiminv(i)
             therm(i)     = heatv(i)*fak/wm(i)
             rino(i,plev) = 0.e+0_fp
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
                tkv = th(i,k)*(1. + 0.61e+0_fp*q(i,k))
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
          pblmin  = 700.e+0_fp*ustar(i)
          pblh(i) = max( pblh(i),pblmin )
       end do

    endif ! if pblh_ar

!------------------------------------------------------------------------
! 	... pblh is now available; do preparation for diffusivity calculation:
!------------------------------------------------------------------------
    do i = 1,plonl
       pblk(i) = 0.e+0_fp
       fak1(i) = ustar(i)*pblh(i)*vk
!------------------------------------------------------------------------
! 	... do additional preparation for unstable cases only, set temperature
!           and moisture perturbations depending on stability.
!------------------------------------------------------------------------
       if( unstbl(i) ) then
          phiminv(i) = (1.e+0_fp - binm*pblh(i)/obklen(i))**onet
          phihinv(i) = sqrt(1.e+0_fp - binh*pblh(i)/obklen(i))
          wm(i)      = ustar(i)*phiminv(i)
          fak2(i)    = wm(i)*pblh(i)*vk
          wstr(i)    = (heatv(i)*g*pblh(i)/thvref(i))**onet
          fak3(i)    = fakn*wstr(i)/wm(i)
          tpert(i)   = max( khfs(i)*fak/wm(i),0.e+0_fp )
          qpert(i)   = max( kshfs(i)*fak/wm(i),0.e+0_fp )
       else
          tpert(i)   = max( khfs(i)*fak/ustar(i),0.e+0_fp )
          qpert(i)   = max( kshfs(i)*fak/ustar(i),0.e+0_fp )
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
          if( zkmin == 0.e+0_fp .and. zp(i) > pblh(i) ) then
             zp(i) = pblh(i)
          end if
          if( zm(i) < pblh(i) ) then
             zmzp = 0.5e+0_fp*(zm(i) + zp(i))
             zh(i) = zmzp/pblh(i)
             zl(i) = zmzp/obklen(i)
             if( zh(i) <= 1.e+0_fp ) then
                zzh(i) = (1.e+0_fp - zh(i))**2
             else
                zzh(i) = 0.e+0_fp
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
             if( zl(i) <= 1.e+0_fp ) then
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
             term    = (1.e+0_fp - betam*zl(i))**onet
             pblk(i) = fak1(i)*zh(i)*zzh(i)*term
             pr(i)   = term/sqrt(1.e+0_fp - betah*zl(i))
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    real(fp), intent(in) :: &
         qm1(plonl,plev,ncnst), & ! initial constituent
         qflx(plonl,ncnst), &     ! sfc q flux into lowest model level
         cc(plonl,plev), &        ! -lower diag coeff.of tri-diag matrix
         term(plonl,plev)         ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
!
! !INPUT/OUTPUT PARAMETERS:
!
    real(fp), intent(inout) :: &
         ze(plonl,plev)           ! term in tri-diag. matrix system
!
! !OUTPUT PARAMETERS:
!
    real(fp), intent(out) :: &
         qp1(plonl,plev,ncnst)    ! final constituent
!
! !REMARKS:
!  Procedure for solution of the implicit equation follows :
!  Richtmyer and Morton (1967,pp 198-199)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp) :: &
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
       tmp1d(i) = 1.e+0_fp/(1.e+0_fp + cc(i,plev) - cc(i,plev)*ze(i,plev-1))
       ze(i,plev) = 0.e+0_fp
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    real(fp), intent(in) :: &
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
    real(fp), intent(inout) :: &
         as2(:,:,:,:)     ! moist, tracers after vert. diff
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp) :: &
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

    real(fp) :: &
         tm1(plonl,plev), &        ! temperature input
         pmidm1(plonl,plev), &     ! midpoint pressures
         pintm1(plonl,plevp), &    ! interface pressures
         rpdel(plonl,plev), &      ! 1./pdel  (thickness bet interfaces)
         rpdeli(plonl,plev), &     ! 1./pdeli (thickness bet midpoints)
         cflx(plonl,pcnst), &      ! surface constituent flux (kg/m2/s)
         kvh(plonl,plevp), &       ! coefficient for heat and tracers
         cgs(plonl,plevp)          ! counter-grad star (cg/flux)
    real(fp) :: &
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
       potbar(:,k) = pintm1(:,k)/(.5e+0_fp*(tm1(:,k) + tm1(:,k-1)))
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
       cah(i,plev) = 0.e+0_fp
    end do
!-----------------------------------------------------------------------
! 	... calculate e(k) for heat vertical diffusion.  this term is
!           required in solution of tridiagonal matrix defined by implicit
!           diffusion eqn.
!-----------------------------------------------------------------------
    do i = 1,plonl
       termh(i,ntopfl) = 1.e+0_fp/(1.e+0_fp + cah(i,ntopfl))
       zeh(i,ntopfl) = cah(i,ntopfl)*termh(i,ntopfl)
    end do
    do k = ntopfl+1,plev-1
       do i = 1,plonl
          termh(i,k) = 1.e+0_fp/(1.e+0_fp + cah(i,k) + cch(i,k) &
                      - cch(i,k)*zeh(i,k-1))
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
    where (qp1 < 0.e+0_fp)
       qp1 = 0.e+0_fp
    endwhere

    !Output values from local variables to arguments.(ccc, 11/17/09)
    as2(:,lat,:,:) = qp1

  END SUBROUTINE VDIFFAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    real(fp), intent(in) :: &
         t(plonl,plev), &        ! temperature (used for density)
         pmid(plonl,plev), &     ! midpoint pressures
         cflx(plonl,pcnst), &     ! surface constituent flux (kg/m2/s)
         cgs(plonl,plevp)        ! counter-gradient star (cg/flux)
!
! !OUTPUT PARAMETERS:
!
    real(fp), intent(out) :: &
         cgq(plonl,plevp,pcnst)  ! counter-gradient term for constituents
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
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
    real(fp) :: &
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
          cgq(:,k,m) = 0.e+0_fp
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE VDINTI( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PRESSURE_MOD,   ONLY : GET_AP, GET_BP
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(GrdState), INTENT(IN ) :: State_Grid ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp), parameter :: pbl_press = 400.e2     ! pressure cap for pbl (pa)
    integer :: k, &                               ! vertical loop index
               m

    real(fp)  :: ref_pmid(State_Grid%NZ)

    !=================================================================
    ! vdinti begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ref_pmid = 0.e+0_fp
    plev  = State_Grid%NZ
    plevp = plev+1
!-----------------------------------------------------------------------
! 	... hard-wired numbers.
!           zkmin = minimum k = kneutral*f(ri)
!-----------------------------------------------------------------------
    zkmin = .01e+0_fp

!-----------------------------------------------------------------------
! 	... set physical constants for vertical diffusion and pbl
!-----------------------------------------------------------------------

    ! REF_PMID is indexed with K=1 being the atm. top and K=PLEV being
    ! the surface.  Eliminate call to UPSIDEDOWN (bmy, 12/21/10)
    do k = 1, plev
       ref_pmid(plev-k+1) = 0.5e+0_fp*( GET_AP(k  )*100.e+0_fp &
           + GET_BP(k  )*1.e+5_fp + &
           GET_AP(k+1)*100.e+0_fp + GET_BP(k+1)*1.e+5_fp )
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
               ' model levels. top is ',1.e-2_fp*ref_pmid(plevp-npbl),' hpa'
    if( plev == 1 ) then
       ntopfl = 0
    else
       ntopfl = 1
    end if

!-----------------------------------------------------------------------
! 	... set the square of the mixing lengths
!-----------------------------------------------------------------------
    ALLOCATE( ml2(plevp), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod:ML2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ml2(1) = 0.e+0_fp
    do k = 2,plev
       ml2(k) = (30.e+0_fp)**2
    end do
    ml2(plevp) = 0.e+0_fp
!-----------------------------------------------------------------------
! 	... set the minimum mixing ratio for the counter-gradient term.
!           normally this should be the same as qmin.
!-----------------------------------------------------------------------

    ALLOCATE( qmincg(pcnst), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod:QMINCG', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    do m = 1,pcnst
       qmincg(m) = 0.e+0_fp
    end do

!-----------------------------------------------------------------------
! 	... initialize pbl variables
!-----------------------------------------------------------------------
    call pbinti( gravit)

  END SUBROUTINE VDINTI
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE VDIFFDR( Input_Opt,  State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
    USE GLOBAL_CH4_MOD,     ONLY : CH4_EMIS
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PBL_MIX_MOD,        ONLY : COMPUTE_PBL_HEIGHT
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CONV, GET_TS_EMIS, GET_TS_CHEM
    USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_DD, ADD_HgP_DD
    USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_SNOWPACK
    USE MERCURY_MOD,        ONLY : HG_EMIS
    USE OCEAN_MERCURY_MOD,  ONLY : Fg !hma
    USE OCEAN_MERCURY_MOD,  ONLY : OMMFP => Fp
    USE OCEAN_MERCURY_MOD,  ONLY : LHg2HalfAerosol !cdh

    ! HEMCO update
    USE HCO_INTERFACE_MOD,  ONLY : GetHcoID, GetHcoVal, GetHcoDiagn

    implicit none
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met    ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REMARKS:
!  (1) Need to declare the Meteorology State object (State_MET) with
!      INTENT(INOUT).  This is because VDIFF will modify the specific
!      humidity field. (bmy, 11/21/12)
!                                                                            .
!  (2) VDIFF also archives drydep fluxes to the soil NOx emissions module
!      (by calling routine SOIL_DRYDEP) and to the ND44 diagnostic.
!                                                                            .
!  (3) As of July 2016, we assume that all of the advected species are listed
!      first in the species database.  This is the easiest way to pass a slab
!      to the TPCORE routine.  This may change in the future. (bmy, 7/13/16)

! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: I, J, L, N, NN, D, NA, nAdvect, ND, nDryDep

    REAL(fp)                :: FRAC_NO_HG0_DEP !jaf
    LOGICAL                 :: ZERO_HG0_DEP !jaf

    real(fp), TARGET, dimension(State_Grid%NX,State_Grid%NY,State_Grid%NZ)   :: pmid, rpdel, rpdeli, zm
    real(fp), TARGET, dimension(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1) :: pint
    real(fp), TARGET, dimension(State_Grid%NX,State_Grid%NY,State_Chm%nAdvect) :: sflx
    real(fp), TARGET, dimension(State_Grid%NX,State_Grid%NY,State_Chm%nAdvect) :: eflx, dflx
                                                              ! surface flux
    real(fp), TARGET, dimension(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1) :: cgs, kvh, kvm
    real(fp), TARGET, dimension(State_Grid%NX,State_Grid%NY)         :: pblh, tpert, qpert
    real(fp), TARGET, dimension(State_Grid%NX,State_Grid%NY,State_Grid%NZ)   :: thp   ! potential temp
    real(fp), TARGET, dimension(State_Grid%NX,State_Grid%NY)         :: shflx ! water vapor flux
    real(fp), TARGET, dimension(State_Grid%NX,State_Grid%NY,State_Grid%NZ)   :: t1
                                                         ! save tracer MR
                                                         ! before vdiffdr
    real(fp) :: vtemp
    real(fp) :: p0 = 1.e+5_fp
    real(fp) :: dtime
    real(fp) :: wk1, wk2
    real(fp) :: soilflux
    integer  :: pbl_top

    REAL(fp) :: DEP_KG !(cdh, 8/28/09)

    ! Array to store a single level of the AS2 array,
    ! so as not to blow up the parallelization (ccc, 12/22.10)
    REAL(fp), dimension(State_Grid%NX, State_Grid%NY, State_Chm%nAdvect)  :: as2_scal

    ! Pointers
    REAL(fp),  POINTER :: p_um1   (:,:,:  )
    REAL(fp),  POINTER :: p_vm1   (:,:,:  )
    REAL(fp),  POINTER :: p_tadv  (:,:,:  )
    REAL(fp),  POINTER :: p_hflux (:,:    )
    REAL(fp),  POINTER :: p_ustar (:,:    )
    REAL(fp),  POINTER :: p_pmid  (:,:,:  )
    REAL(fp),  POINTER :: p_pint  (:,:,:  )
    REAL(fp),  POINTER :: p_rpdel (:,:,:  )
    REAL(fp),  POINTER :: p_rpdeli(:,:,:  )
    REAL(fp),  POINTER :: p_zm    (:,:,:  )
    REAL(fp),  POINTER :: p_thp   (:,:,:  )
    REAL(fp),  POINTER :: p_kvh   (:,:,:  )
    REAL(fp),  POINTER :: p_kvm   (:,:,:  )
    REAL(fp),  POINTER :: p_cgs   (:,:,:  )
    REAL(fp),  POINTER :: p_shp   (:,:,:  )
    REAL(fp),  POINTER :: p_t1    (:,:,:  )
    REAL(fp),  POINTER :: p_as2   (:,:,:,:)

    REAL(fp),  POINTER :: DEPSAV  (:,:,:  )  ! IM, JM, nDryDep

    ! For values from Input_Opt
    LOGICAL            :: IS_CH4,    IS_FULLCHEM, IS_Hg
    LOGICAL            :: IS_TAGCO,  IS_AEROSOL,  IS_RnPbBe, LDYNOCEAN
    LOGICAL            :: LGTMM,     LSOILNOX

    ! HEMCO update
    LOGICAL            :: FND
    REAL(fp)           :: TMPFLX, EMIS, DEP
    INTEGER            :: HCRC,   TOPMIX

    ! PARANOX loss fluxes (kg/m2/s), imported from
    ! HEMCO PARANOX extension module (ckeller, 4/15/2015)
    REAL(f4), POINTER, SAVE :: PNOXLOSS_O3  (:,:) => NULL()
    REAL(f4), POINTER, SAVE :: PNOXLOSS_HNO3(:,:) => NULL()

    ! SAVEd scalars
    LOGICAL,           SAVE :: FIRST = .TRUE.
    INTEGER,           SAVE :: id_O3
    INTEGER,           SAVE :: id_HNO3

    ! For pointing to the species database
    TYPE(Species),  POINTER :: SpcInfo
    INTEGER                 :: Hg_Cat

    ! For diagnostics
    REAL(fp)                :: EmMW_kg

    ! For error trapping
    CHARACTER(LEN=255)      :: ErrMsg, ThisLoc

    !=================================================================
    ! vdiffdr begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS               ! Assume success
    nAdvect = State_Chm%nAdvect        ! # of advected species
    ErrMsg  = ''
    ThisLoc = ' -> at VDIFF (in module GeosCore/vdiff_mod.F90)'

    !### Debug
    IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: VDIFFDR begins' )

    ! NOTE: The prior behavior of the code assumed that NUMDEP=0 when
    ! drydep was shut off.  NUMDEP=0 prevented the main drydep loops
    ! from occurring, thus avoiding seg faults.  We now replicate the
    ! same behavior here. (bmy, 7/27/16)
    IF ( Input_Opt%LDRYD ) THEN
       nDryDep = State_Chm%nDryDep
    ELSE
       nDryDep = 0
    ENDIF

    ! Initialize pointers
    SpcInfo => NULL()
    DEPSAV  => State_Chm%DryDepSav

    ! Initialize local arrays. (ccc, 12/21/10)
    pmid    = 0e+0_fp
    rpdel   = 0e+0_fp
    rpdeli  = 0e+0_fp
    zm      = 0e+0_fp
    pint    = 0e+0_fp
    sflx    = 0e+0_fp
    eflx    = 0e+0_fp
    dflx    = 0e+0_fp
    soilflux = 0e+0_fp
    cgs     = 0e+0_fp
    kvh     = 0e+0_fp
    kvm     = 0e+0_fp
    pblh    = 0e+0_fp
    tpert   = 0e+0_fp
    qpert   = 0e+0_fp
    thp     = 0e+0_fp
    shflx   = 0e+0_fp
    t1      = 0e+0_fp
    as2_scal= 0e+0_fp

    ! Copy values from Input_Opt (bmy, 8/1/13)
    IS_CH4       = Input_Opt%ITS_A_CH4_SIM
    IS_Hg        = Input_Opt%ITS_A_MERCURY_SIM
    IS_AEROSOL   = Input_Opt%ITS_AN_AEROSOL_SIM
    LDYNOCEAN    = Input_Opt%LDYNOCEAN
    LGTMM        = Input_Opt%LGTMM
    LSOILNOX     = Input_Opt%LSOILNOX

    dtime = GET_TS_CONV() ! second

    shflx = State_Met%EFLUX / latvap ! latent heat -> water vapor flux

    ! First-time setup
    IF ( FIRST ) THEN

       ! Get species IDs
       id_O3   = Ind_('O3'  )
       id_HNO3 = Ind_('HNO3')

       ! On first call, get pointers to the PARANOX loss fluxes. These are
       ! stored in diagnostics 'PARANOX_O3_DEPOSITION_FLUX' and
       ! 'PARANOX_HNO3_DEPOSITION_FLUX'. The call below links pointers
       ! PNOXLOSS_O3 and PNOXLOSS_HNO3 to the data values stored in the
       ! respective diagnostics. The pointers will remain unassociated if
       ! the diagnostics do not exist (ckeller, 4/10/2015).
       CALL GetHcoDiagn( 'PARANOX_O3_DEPOSITION_FLUX'  , &
                         .FALSE.,   HCRC, Ptr2D = PNOXLOSS_O3 )
       CALL GetHcoDiagn( 'PARANOX_HNO3_DEPOSITION_FLUX', &
                         .FALSE.,   HCRC, Ptr2D = PNOXLOSS_HNO3 )

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

! (Turn off parallelization for now, skim 6/20/12)

!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L )
    do J = 1, State_Grid%NY
    do I = 1, State_Grid%NX

    ! calculate variables related to pressure
    do L = 1, State_Grid%NZ
       pmid(I,J,L) = State_Met%PMID(I,J,L)*100.e+0_fp  ! hPa -> Pa
       pint(I,J,L) = State_Met%PEDGE(I,J,L)*100.e+0_fp ! hPa -> Pa
       ! calculate potential temperature
       thp(I,J,L) = State_Met%T(I,J,L)*(p0/pmid(I,J,L))**cappa
    enddo
    pint(I,J,State_Grid%NZ+1) = State_Met%PEDGE(I,J,State_Grid%NZ+1)

    enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L )
    do J = 1, State_Grid%NY
    do I = 1, State_Grid%NX
    do L = 1, State_Grid%NZ
       ! Corrected calculation of zm.
       ! Box height calculation now uses virtual temperature.
       ! Therefore, use virtual temperature in hypsometric equation.
       ! (ewl, 3/3/15)
       zm(I,J,L) = sum( State_Met%BXHEIGHT(I,J,1:L)) &
                 - log( pmid(I,J,L)/pint(I,J,L+1) )  &
                 * r_g * State_Met%TV(I,J,L)

    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    ! Have to separate the calculations of pmid/pint and rpdel/rpdeli.
    ! (Lin, 06/02/08)
!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L )
    do J = 1, State_Grid%NY
    do I = 1, State_Grid%NX
       do L = 1, State_Grid%NZ
          rpdel(I,J,L) = 1.e+0_fp / (pint(I,J,L) - pint(I,J,L+1))
       enddo

       !rpdeli(I,J,1) = 1.e+0_fp / (PS(I,J) - pmid(I,J,1))
       rpdeli(I,J,1) = 0.e+0_fp ! follow mozart setup (shown in mo_physlic.F90)

       do L = 2, State_Grid%NZ
          rpdeli(I,J,L) = 1.e+0_fp / (pmid(I,J,L-1) - pmid(I,J,L))
       enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    !!! calculate surface flux = emissions - dry deposition !!!

    ! Define slice of AS2, so as not to blow up the parallelization
    ! (ccc, bmy, 12/20/10)
    !
    ! NOTE: For now, so as to avoid having to rewrite the internals
    ! of the TPCORE routines, just point to 1:nAdvect entries of
    ! State_Chm%Species.  This is OK for now, as of July 2016, all of
    ! the advected species are listed first.  This may change in the
    ! future, but we'll worry about that later. (bmy, 7/13/16)
    as2_scal = State_Chm%Species(:,:,1,1:nAdvect)

!$OMP PARALLEL DO                                                         &
!$OMP DEFAULT( SHARED )                                                   &
!$OMP PRIVATE( I,      J,               L,            N                 ) &
!$OMP PRIVATE( WK1,    WK2,             PBL_TOP,      DEP_KG,  TOPMIX   ) &
!$OMP PRIVATE( fnd,    emis,            dep,          NA,      ND       ) &
!$OMP PRIVATE( TMPFLX, FRAC_NO_HG0_DEP, ZERO_HG0_DEP, SpcInfo, Hg_Cat   )
    do J = 1, State_Grid%NY
    do I = 1, State_Grid%NX

       ! PBL top level [integral model levels]
       topmix      = MAX( 1, FLOOR( State_Met%PBL_TOP_L(I,J) ) )

       !--------------------------------------------------------------------
       ! Add emissions & deposition values calculated in HEMCO.
       ! Here we only consider emissions below the PBL top.
       !
       ! For the full-chemistry simulations, emissions above the PBL
       ! top will be applied in routine SETEMIS, which occurs just
       ! before the SMVGEAR/KPP solvers are invoked.
       !
       ! For the specialty simulations, emissions above the PBL top
       ! will be applied in the chemistry routines for each
       ! specialty simulation.
       !
       ! For more information, please see this wiki page:
       ! http://wiki.geos-chem.org/Distributing_emissions_in_the_PBL
       !--------------------------------------------------------------------
       DO NA = 1, nAdvect

          ! Add total emissions in the PBL to the EFLX array
          ! which tracks emission fluxes.  Units are [kg/m2/s].
          tmpflx = 0.0e+0_fp
          DO L = 1, TOPMIX
             CALL GetHcoVal ( NA, I, J, L, fnd, emis=emis )
             IF ( .NOT. fnd ) EXIT
             tmpflx = tmpflx + emis
          ENDDO
          eflx(I,J,NA) = eflx(I,J,NA) + tmpflx

          ! Also add drydep frequencies calculated by HEMCO to the DFLX
          ! array. These values are stored in 1/s. They are added in the
          ! same manner as the DEPSAV values from drydep_mod.F90.
          ! DFLX will be converted to kg/m2/s lateron. (ckeller, 04/01/2014)
          CALL GetHcoVal ( NA, I, J, 1, fnd, dep=dep )
          IF ( fnd ) THEN
             dflx(I,J,NA) = dflx(I,J,NA) + ( dep * as2_scal(I,J,NA)  &
                            / ( AIRMW                                &
                                / State_Chm%SpcData(NA)%Info%emMW_g ) )
          ENDIF
       ENDDO

       !--------------------------------------------------------------------
       ! Overwrite emissions for offline CH4 simulation.
       ! CH4 emissions become stored in CH4_EMIS in global_ch4_mod.F90.
       ! We use CH4_EMIS here instead of the HEMCO internal emissions
       ! only to make sure that total CH4 emissions are properly defined
       ! in a multi-tracer CH4 simulation. For a single-tracer simulation
       ! and/or all other source types, we could use the HEMCO internal
       ! values set above and would not need the code below.
       ! Units are already in kg/m2/s. (ckeller, 10/21/2014)
       !--------------------------------------------------------------------
       IF ( IS_CH4 ) THEN
          do NA = 1, nAdvect
             eflx(I,J,NA) = CH4_EMIS(I,J,NA)
          enddo
       ENDIF

       !--------------------------------------------------------------------
       ! Overwrite emissions for offline mercury simulation
       ! HG emissions become stored in HG_EMIS in mercury_mod.F90.
       ! This is a workaround to ensure backwards compatibility.
       ! Units are already in kg/m2/s. (ckeller, 10/21/2014)
       !--------------------------------------------------------------------
       IF ( IS_Hg ) THEN
          do NA = 1, nAdvect
             eflx(I,J,NA) = HG_EMIS(I,J,NA)
          enddo
       ENDIF

       !--------------------------------------------------------------------
       ! Apply dry deposition frequencies
       ! These are the frequencies calculated in drydep_mod.F90
       ! The HEMCO drydep frequencies (from air-sea exchange and
       ! PARANOX) were already added above.
       !
       ! NOTES:
       ! (1) Loops over only the drydep species
       ! (2) If drydep is turned off, nDryDep=0 and the loop won't execute
       ! (3) Tagged species are included in this loop. via species database
       !--------------------------------------------------------------------
       DO ND = 1, nDryDep

          ! Get the species ID from the drydep ID
          N = State_Chm%Map_DryDep(ND)

          IF ( N <= 0 ) CYCLE

          ! Point to the corresponding Species Database entry
          SpcInfo => State_Chm%SpcData(N)%Info

          ! use mean concentration within the PBL for calculating drydep
          ! fluxes
          if (pbl_mean_drydep) then
             wk1 = 0.e+0_fp
             wk2 = 0.e+0_fp
             pbl_top = State_Met%PBL_MAX_L ! the highest layer the PBL reaches,
                                           ! globally
             do L = 1, pbl_top
                wk1 = wk1 + State_Chm%Species    (I,J,L,N) * &
                            State_Met%AD         (I,J,L  ) * &
                            State_Met%F_UNDER_PBLTOP(I,J,L)

                wk2 = wk2 + State_Met%AD         (I,J,L)   * &
                            State_Met%F_UNDER_PBLTOP(I,J,L)
             enddo
             ! since we only use the ratio of wk1 / wk2, there should not be
             ! a problem even if the PBL top is lower than the top of the
             ! first (lowest) model layer
             ! given that as2 is in v/v
             ! Now add to existing dflx (ckeller, 10/16/2014).
             dflx(I,J,NN) = dflx(I,J,NN) &
                          + DEPSAV(I,J,N) * (wk1/(wk2+1.e-30_fp))  &
                          / ( AIRMW / SpcInfo%emMW_g )


             ! consistency with the standard GEOS-Chem setup (Lin, 07/14/08)
             if (drydep_back_cons) then
                dflx(I,J,N) = dflx(I,J,N) * (wk2+1.e-30_fp) / &
                               State_Met%AD(I,J,1)         * &
                               State_Met%BXHEIGHT(I,J,1)   / &
                               State_Met%PBL_TOP_m(I,J)
             endif
          else

             ! only use the lowest model layer for calculating drydep fluxes
             ! given that as2 is in v/v
             ! NOTE: Now use as2_scal(I,J,NN), instead of as2(I,J,1,NN) to
             ! avoid seg faults in parallelization (ccarouge, bmy, 12/20/10)
             ! Now add to existing dflx (ckeller, 10/16/2014).
             dflx(I,J,N) = dflx(I,J,N) &
                         + DEPSAV(I,J,ND) * as2_scal(I,J,N) /   &
                         ( AIRMW / SpcInfo%emMW_g )
          endif

          ! Hg(0) exchange with the ocean is handled by ocean_mercury_mod
          ! so disable deposition over water here.
          ! Turn off Hg(0) deposition to snow and ice because we haven't yet
          ! included emission from these surfaces and most field studies
          ! suggest Hg(0) emissions exceed deposition during sunlit hours.
          FRAC_NO_HG0_DEP = MIN( State_Met%FROCEAN(I,J) + &
                                 State_Met%FRSNO(I,J)   + &
                                 State_Met%FRLANDIC(I,J), 1e+0_fp)
          ZERO_HG0_DEP    = ( FRAC_NO_HG0_DEP > 0e+0_fp )

          IF ( IS_Hg .AND. SpcInfo%Is_Hg0 ) THEN
             IF ( ZERO_HG0_DEP ) THEN
                DFLX(I,J,N) = DFLX(I,J,N) * &
                               MAX(1e+0_fp - FRAC_NO_HG0_DEP, 0e+0_fp)
             ENDIF
          ENDIF

          ! Free species database pointer
          SpcInfo => NULL()
       enddo

       ! virtual temperature in the lowest model layer
       ! vtemp = tadv(I,J,1)*(1. + zvir*shp(I,J,1))
       ! for deposition: additional step to convert from s-1 to kg/m2/s
       ! dflx(I,J,:) = dflx(I,J,:) * pmid(I,J,1) / rair / vtemp * BXHEIGHT(I,J,1)
       ! alternate method to convert from s-1 to kg/m2/s
       dflx(I,J,:) = dflx(I,J,:) * State_Met%AD(I,J,1) / &
                     State_Grid%Area_M2(I,J)

       ! Now that dflx is in kg/m2/s, add PARANOX loss to this term. The PARANOX
       ! loss term is already in kg/m2/s. PARANOX loss (deposition) is calculated
       ! for O3 and HNO3 by the PARANOX module, and data is exchanged via the
       ! HEMCO diagnostics. The data pointers PNOXLOSS_O3 and PNOXLOSS_HNO3 have
       ! been linked to these diagnostics at the beginning of this routine
       ! (ckeller, 4/10/15).
       IF ( ASSOCIATED( PNOXLOSS_O3 ) .AND. id_O3 > 0 ) THEN
          dflx(I,J,id_O3) = dflx(I,J,id_O3) + PNOXLOSS_O3(I,J)
       ENDIF
       IF ( ASSOCIATED( PNOXLOSS_HNO3 ) .AND. id_HNO3 > 0 ) THEN
          dflx(I,J,id_HNO3) = dflx(I,J,id_HNO3) + PNOXLOSS_HNO3(I,J)
       ENDIF

       ! surface flux = emissions - dry deposition
       sflx(I,J,:) = eflx(I,J,:) - dflx(I,J,:) ! kg/m2/s

       !--------------------------------------------------------------------
       ! Archive Hg deposition for surface reservoirs (cdh, 08/28/09)
       !--------------------------------------------------------------------
       IF ( IS_Hg ) THEN

          ! Loop over only the drydep species
          ! If drydep is turned off, nDryDep=0 and the loop won't execute
          DO ND = 1, nDryDep

             ! Get the species ID from the drydep ID
             N = State_Chm%Map_DryDep(ND)

             ! Point to the Species Database entry for tracer N
             SpcInfo => State_Chm%SpcData(N)%Info

             ! Deposition mass, kg
             DEP_KG = dflx( I, J, N ) * State_Grid%Area_M2(I,J) &
                    * GET_TS_CONV()

             IF ( SpcInfo%Is_Hg2 ) THEN

                ! Get the category number for this Hg2 tracer
                Hg_Cat = SpcInfo%Hg_Cat

                ! Archive dry-deposited Hg2
                CALL ADD_Hg2_DD      ( I, J, Hg_Cat,    DEP_KG              )
                CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat,    DEP_KG,              &
                                       State_Met, State_Chm, State_Diag     )

             ELSE IF ( SpcInfo%Is_HgP ) THEN

                ! Get the category number for this HgP tracer
                Hg_Cat = SpcInfo%Hg_Cat

                ! Archive dry-deposited HgP
                CALL ADD_HgP_DD      ( I, J, Hg_Cat,    DEP_KG              )
                CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat,    DEP_KG,              &
                                       State_Met, State_Chm, State_Diag     )

             ENDIF

             ! Free pointer
             SpcInfo => NULL()
          ENDDO
       ENDIF

    enddo
    enddo
!$OMP END PARALLEL DO

    !IF ( prtDebug ) THEN
    !   ! Print debug output for the advected species
    !   write(*,*) 'eflx and dflx values HEMCO [kg/m2/s]'
    !   do NA = 1, nAdvect
    !      write(*,*) 'eflx TRACER ', NA, ': ', SUM(eflx(:,:,NA))
    !      write(*,*) 'dflx TRACER ', NA, ': ', SUM(dflx(:,:,NA))
    !   enddo
    !ENDIF

    !=======================================================================
    ! DIAGNOSTICS: Compute drydep flux loss due to mixing [molec/cm2/s]
    !
    ! NOTE: Dry deposition of "tagged" species (e.g. in tagO3, tagCO, tagHg
    ! specialty simulations) are accounted for in species 1..nDrydep,
    ! so we don't need to do any further special handling.
    !=======================================================================
    if ( Do_ND44 .or. LGTMM .or. LSOILNOX .or.  &
         State_Diag%Archive_DryDepMix     .or.  &
         State_Diag%Archive_DryDep       ) then

       ! Loop over only the drydep species
       ! If drydep is turned off, nDryDep=0 and the loop won't execute
       DO ND = 1, nDryDep

          ! Get the species ID from the drydep ID
          N = State_Chm%Map_DryDep(ND)

          ! Skip if not a valid species
          IF ( N <= 0 ) CYCLE

          ! Point to the Species Database entry for this tracer
          ! NOTE: Assumes a 1:1 tracer index to species index mapping
          SpcInfo => State_Chm%SpcData(N)%Info

          ! Get the (emitted) molecular weight of the species in kg
          EmMW_kg = SpcInfo%emMW_g * 1.e-3_fp

          !-----------------------------------------------------------------
          ! HISTORY: Update dry deposition flux loss [molec/cm2/s]
          !
          ! DFLX is in kg/m2/s.  We convert to molec/cm2/s by:
          !
          ! (1) multiplying by 1e-4 cm2/m2        => kg/cm2/s
          ! (2) multiplying by ( AVO / EmMW_KG )  => molec/cm2/s
          !
          ! The term AVO/EmMW_kg = (molec/mol) / (kg/mol) = molec/kg
          !
          ! NOTE: we don't need to multiply by the ratio of TS_CONV /
          ! TS_CHEM, as the updating frequency for HISTORY is determined
          ! by the "frequency" setting in the "HISTORY.rc"input file.
          ! The old bpch diagnostics archived the drydep due to chemistry
          ! every chemistry timestep = 2X the dynamic timestep.  So in
          ! order to avoid double-counting the drydep flux from mixing,
          ! you had to multiply by TS_CONV / TS_CHEM.
          !
          ! ALSO NOTE: When comparing History output to bpch output,
          ! you must use an updating frequency equal to the dynamic
          ! timestep so that the drydep fluxes due to mixing will
          ! be equivalent w/ the bpch output.  It is also recommended to
          ! turn off chemistry so as to be able to compare the drydep
          ! fluxes due to mixing in bpch vs. History as an "apples-to-
          ! apples" comparison.
          !
          !    -- Bob Yantosca (yantosca@seas.harvard.edu)
          !-----------------------------------------------------------------
          IF ( State_Diag%Archive_DryDepMix .OR. &
               State_Diag%Archive_DryDep ) THEN
             State_Diag%DryDepMix(:,:,ND) = Dflx(:,:,N)               &
                                            * 1.0e-4_fp               &
                                            * ( AVO / EmMW_kg  )
          ENDIF

          !-----------------------------------------------------------------
          ! If Soil NOx is turned on, then call SOIL_DRYDEP to
          ! archive dry deposition fluxes for nitrogen species
          ! (SOIL_DRYDEP will exit if it can't find a match.
          !
          ! Locate position of each tracer in DEPSAV
          ! (get tracer id)
          ! NOTE: trc_id was previously NN and drydep id D
          ! was previously N. This is changed for convention consistency
          ! within subroutine (ewl, 1/25/16)
          !-----------------------------------------------------------------
	  IF ( LSOILNOX ) THEN
             soilflux = 0e+0_fp
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX
                soilflux = dflx(I,J,N) &
	          / ( SpcInfo%emMW_g * 1.e-3_fp ) &
                  * AVO * 1.e-4_fp &
                  * GET_TS_CONV() / GET_TS_EMIS()

                CALL SOIL_DRYDEP ( I, J, 1, N, soilflux, State_Chm )
             ENDDO
             ENDDO
	  ENDIF

          ! Free species database pointer
          SpcInfo => NULL()

       enddo ! D

    endif

	!Maasa, Add SoilNOx deposition to allow SN code to work with NLPBL on.

    !### Debug
    IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: after emis. and depdrp' )

    if( divdiff ) then

       if ( pblh_ar ) then
       do J = 1, State_Grid%NY
       do I = 1, State_Grid%NX
         pblh(I,J) = State_Met%PBL_TOP_m(I,J) ! obtain archived PBLH
       enddo
       enddo
       endif

       !--------------------------------------------------------------------
       ! Now use pointers to flip arrays in the vertical (bmy, 6/22/15)
       !--------------------------------------------------------------------

       ! 3-D fields on level centers
       p_um1              => State_Met%U    ( :, :, State_Grid%NZ  :1:-1    )
       p_vm1              => State_Met%V    ( :, :, State_Grid%NZ  :1:-1    )
       p_tadv             => State_Met%T    ( :, :, State_Grid%NZ  :1:-1    )
       p_hflux            => State_Met%HFLUX
       p_ustar            => State_Met%USTAR
       p_pmid             => pmid           ( :, :, State_Grid%NZ  :1:-1    )
       p_rpdel            => rpdel          ( :, :, State_Grid%NZ  :1:-1    )
       p_rpdeli           => rpdeli         ( :, :, State_Grid%NZ  :1:-1    )
       p_zm               => zm             ( :, :, State_Grid%NZ  :1:-1    )
       p_thp              => thp            ( :, :, State_Grid%NZ  :1:-1    )
       p_shp              => State_Met%SPHU ( :, :, State_Grid%NZ  :1:-1    )

       ! 3-D fields on level edges
       p_pint             => pint           ( :, :, State_Grid%NZ+1:1:-1    )
       p_kvh              => kvh            ( :, :, State_Grid%NZ+1:1:-1    )
       p_kvm              => kvm            ( :, :, State_Grid%NZ+1:1:-1    )
       p_cgs              => cgs            ( :, :, State_Grid%NZ+1:1:-1    )

       ! Species concentration fields
       !
       ! NOTE: For now, so as to avoid having to rewrite the internals
       ! of the VDIFF routines, just point to 1:nAdvect entries of
       ! State_Chm%Species.  This is OK for now, as of July 2016, all of
       ! the advected species are listed first.  This may change in the
       ! future, but we'll worry about that later. (bmy, 7/13/16)
       p_as2              => State_Chm%Species( :, :, State_Grid%NZ:1:-1, 1:nAdvect )

       ! Convert v/v -> m/m (i.e., kg/kg)
       DO NA = 1, nAdvect
          p_as2(:,:,:,NA) =  p_as2(:,:,:,NA) / ( AIRMW       &
                             / State_Chm%SpcData(NA)%Info%emMW_g )
       ENDDO

       ! Convert g/kg -> kg/kg
       p_shp              =  p_shp * 1.e-3_fp

       !### Debug
       IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: before vdiff' )

!$OMP PARALLEL DO DEFAULT( SHARED )      &
!$OMP PRIVATE( J )
       do J = 1, State_Grid%NY
          call vdiff( J,         1,         p_um1,      p_vm1,               &
                      p_tadv,    p_pmid,    p_pint,     p_rpdel,             &
                      p_rpdeli,  dtime,     p_zm,       p_hflux,             &
                      sflx,      p_thp,     p_as2,      pblh,                &
                      p_kvh,     p_kvm,     tpert,      qpert,               &
                      p_cgs,     p_shp,     shflx,      State_Grid%NX,       &
                      Input_Opt, State_Met, State_Chm,  State_Diag,          &
                      ustar_arg=p_ustar,    RC=RC                           )
       enddo
!$OMP END PARALLEL DO

       !### Debug
       IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: after vdiff' )

       ! Convert kg/kg -> v/v
       DO NA = 1, nAdvect
          p_as2(:,:,:,NA) = p_as2(:,:,:,NA) * ( AIRMW          &
                            / State_Chm%SpcData(NA)%Info%emMW_g )
       ENDDO

       ! Convert kg/kg -> g/kg
       p_shp    = p_shp * 1.e+3_fp

       ! Free pointers
       NULLIFY( p_um1,   p_vm1,    p_tadv, p_pmid, p_pint )
       NULLIFY( p_rpdel, p_rpdeli, p_zm,   p_thp,  p_cgs  )
       NULLIFY( p_kvh,   p_kvm,    p_shp,  p_as2,  p_hflux)

    else if( arvdiff ) then
!---------------------------------------------------------------------------
!  	... vertical diffusion using archived values of cgs and kvh.
!
!       %%% NOTE: THIS SECTION IS NORMALLY NOT EXECUTED %%%
!       %%% BECAUSE ARVDIFF IS SET TO .FALSE. ABOVE     %%%
!--------------------------------------------------------------------------

       !-------------------------------------------------------------------
       ! Now use pointers to flip arrays in the vertical (bmy, 6/22/15)
       !-------------------------------------------------------------------

       ! INPUTS: 3-D fields on level centers
       p_tadv   => State_Met%T( :, :, State_Grid%NZ  :1:-1   )
       p_pmid   => pmid       ( :, :, State_Grid%NZ  :1:-1   )
       p_rpdel  => rpdel      ( :, :, State_Grid%NZ  :1:-1   )
       p_rpdeli => rpdeli     ( :, :, State_Grid%NZ  :1:-1   )

       ! INPUTS: 3-D fields on level edges
       p_pint   => pint       ( :, :, State_Grid%NZ+1:1:-1   )
       p_kvh    => kvh        ( :, :, State_Grid%NZ+1:1:-1   )
       p_cgs    => cgs        ( :, :, State_Grid%NZ+1:1:-1   )

       ! INPUTS: Tracer concentration fields
       ! NOTE: For now, so as to avoid having to rewrite the internals
       ! of the VDIFF routines, just point to 1:nAdvect entries of
       ! State_Chm%Species.  This is OK for now, as of July 2016, all of
       ! the advected species are listed first.  This may change in the
       ! future, but we'll worry about that later. (bmy, 7/13/16)
       p_as2    => State_Chm%Species( :, :, State_Grid%NZ:1:-1, 1:nAdvect )

       ! Convert from v/v -> m/m (i.e., kg/kg)
       do NA = 1, nAdvect
          p_as2(:,:,:,NA) = p_as2(:,:,:,NA) / ( AIRMW           &
                            / State_Chm%SpcData(NA)%Info%emMW_g )
       enddo

       !### Debug
       IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: before vdiffar' )

!!$OMP PARALLEL DO DEFAULT( SHARED )   &
!!$OMP PRIVATE( J )
       do J = 1, State_Grid%NY
          call vdiffar( J,      p_tadv, p_pmid, p_pint, p_rpdel, p_rpdeli,     &
                        dtime,  sflx,   p_as2,  p_kvh,  p_cgs,   State_Grid%NX )
      enddo
!!$OMP END PARALLEL DO

       !### Debug
       IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: after vdiffar' )

       ! Convert from m/m (i.e. kg/kg) -> v/v
       do NA = 1, nAdvect
          p_as2(:,:,:,NA) = p_as2(:,:,:,NA) * ( AIRMW /       &
                            State_Chm%SpcData(NA)%Info%emMW_g )
       enddo

       ! Free pointers
       NULLIFY( p_tadv, p_pmid, p_rpdel, p_rpdeli )
       NULLIFY( p_pint, p_kvh,  p_cgs,   p_as2    )

    end if

    !-----------------------------------------------------------------------
    ! re-compute PBL variables wrt derived pblh (in m)
    !-----------------------------------------------------------------------
    if (.not. pblh_ar) then

       ! PBL is in m
       State_Met%PBLH = pblh

       ! Compute PBL quantities
       CALL COMPUTE_PBL_HEIGHT( State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "COMPUTE_PBL_HEIGHT"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    endif

    !### Debug
    IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: VDIFFDR finished' )

    ! Nullify pointers
    NULLIFY( DEPSAV )

  END SUBROUTINE VDIFFDR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE DO_PBL_MIX_2( DO_VDIFF,   Input_Opt,  State_Chm, &
                           State_Diag, State_Grid, State_Met, RC )
!
! !USES:
!
    USE Calc_Met_Mod,       ONLY : AIRQNT
    USE Diagnostics_Mod,    ONLY : Compute_Column_Mass
    USE Diagnostics_Mod,    ONLY : Compute_Budget_Diagnostics
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PBL_MIX_MOD,        ONLY : INIT_PBL_MIX
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : ITS_TIME_FOR_EMIS
    USE Time_Mod,           ONLY : Get_Ts_Dyn
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
#ifdef TOMAS
#ifdef BPCH_DIAG
    !=======================================================================
    ! These are only needed if GEOS-Chem is compiled for TOMAS
    !=======================================================================
    USE CMN_DIAG_MOD,       ONLY : ND44
#endif
#endif

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: DO_VDIFF     ! Switch which turns on PBL
                                                  !  mixing of tracers
    TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met    ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  11 Feb 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    LOGICAL            :: prtDebug

    ! Strings
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    REAL(fp)           :: DT_Dyn

    !=======================================================================
    ! DO_PBL_MIX_2 begins here!
    !=======================================================================

    ! Initialize
    RC       =  GC_SUCCESS                                ! Assume success
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) ! Print debug output?
    ErrMsg   = ''
    ThisLoc  = ' -> at DO_PBL_MIX_2 (in module GeosCore/vdiff_mod.F90)'

    !-----------------------------------------------------------------------
    ! First-time initialization
    ! NOTE: Should really move this into the init stage
    !-----------------------------------------------------------------------
    IF ( FIRST ) THEN

       ! Make sure the various PBL arrays are allocated
       CALL INIT_PBL_MIX( Input_Opt, State_Grid, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "INIT_PBL_MIX"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Perform more internal initializations
       CALL Vdinti( Input_Opt, State_Grid, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "VDINTI"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

#ifdef TOMAS
#ifdef BPCH_DIAG
       ! Set a flag to denote we should archive ND44 bpch diagnostic
       ! NOTE: this will only be valid if BPCH_DIAG=y
       Do_ND44 = ( ND44 > 0 )
#endif
#endif

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !-----------------------------------------------------------------------
    ! Do mixing of species in the PBL (if necessary)
    !-----------------------------------------------------------------------
    IF ( DO_VDIFF ) THEN

       !------------------------------------------------------
       ! Non-local PBL mixing budget diagnostics - Part 1 of 2
       !
       ! WARNING: The mixing budget diagnostic includes the application
       ! of species tendencies (emissions fluxes and dry deposition
       ! rates) below the PBL when using non-local PBL mixing. This is
       ! done for all species with defined emissions / dry deposition
       ! rates, including dust. These tendencies below the PBL are
       ! therefore not included in the emissions/dry deposition budget
       ! diagnostics when using non-local PBL mixing. (ewl, 9/26/18)
       !------------------------------------------------------
       IF ( State_Diag%Archive_BudgetMixing ) THEN
          ! Get initial column masses
          CALL Compute_Column_Mass( Input_Opt,                           &
                                    State_Chm, State_Grid, State_Met,    &
                                    State_Chm%Map_Advect,                &
                                    State_Diag%Archive_BudgetMixingFull, &
                                    State_Diag%Archive_BudgetMixingTrop, &
                                    State_Diag%Archive_BudgetMixingPBL,  &
                                    State_Diag%BudgetMass1,              &
                                    RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Mixing budget diagnostics error 1 (non-local mixing)'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       !----------------------------------------
       ! Unit conversion #1
       !----------------------------------------

       ! Convert species concentration to v/v dry
       CALL Convert_Spc_Units( Input_Opt, State_Chm, State_grid, State_Met, &
                               'v/v dry', RC, OrigUnit=OrigUnit )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "Convert_Spc_Units" (to v/v dry)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !----------------------------------------
       ! PBL mixing
       !----------------------------------------

       ! Do non-local PBL mixing
       CALL VDIFFDR( Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "VDIFFDR"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Debug print
       IF( prtDebug ) THEN
          CALL DEBUG_MSG( '### DO_PBL_MIX_2: after VDIFFDR' )
       ENDIF

       !----------------------------------------
       ! Update airmass etc. and mixing ratios
       !----------------------------------------

       ! Update air quantities and species concentrations with updated
       ! specific humidity (ewl, 10/28/15)
       !
       ! NOTE: Prior to October 2015, air quantities were not updated
       ! with specific humidity modified in VDIFFDR at this point in
       ! the model
       CALL AIRQNT( Input_Opt, State_Chm, State_Grid, State_Met, &
                    RC, Update_Mixing_Ratio=.TRUE. )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "AIRQNT"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Debug print
       IF( prtDebug ) THEN
          CALL DEBUG_MSG( '### DO_PBL_MIX_2: after AIRQNT' )
       ENDIF

       !----------------------------------------
       ! Unit conversion #2
       !----------------------------------------

       ! Convert species back to the original units
       CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                               OrigUnit,  RC )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "Convert_Spc_Units" (from v/v dry)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !------------------------------------------------------
       ! Non-local PBL mixing budget diagnostics - Part 2 of 2
       !------------------------------------------------------
       IF ( State_Diag%Archive_BudgetMixing ) THEN

          ! Get dynamics timestep [s]
          DT_Dyn = Get_Ts_Dyn()

          ! Get final column masses and compute diagnostics
          CALL Compute_Column_Mass( Input_Opt,                              &
                                    State_Chm, State_Grid, State_Met,       &
                                    State_Chm%Map_Advect,                   &
                                    State_Diag%Archive_BudgetMixingFull,    &
                                    State_Diag%Archive_BudgetMixingTrop,    &
                                    State_Diag%Archive_BudgetMixingPBL,     &
                                    State_Diag%BudgetMass2,                 &
                                    RC )
          CALL Compute_Budget_Diagnostics( State_Grid,                      &
                                       State_Chm%Map_Advect,                &
                                       DT_Dyn,                              &
                                       State_Diag%Archive_BudgetMixingFull, &
                                       State_Diag%Archive_BudgetMixingTrop, &
                                       State_Diag%Archive_BudgetMixingPBL,  &
                                       State_Diag%BudgetMixingFull,         &
                                       State_Diag%BudgetMixingTrop,         &
                                       State_Diag%BudgetMixingPBL,          &
                                       State_Diag%BudgetMass1,              &
                                       State_Diag%BudgetMass2,              &
                                       RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Mixing budget diagnostics error 2 (non-local mixing)'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

    ENDIF

  END SUBROUTINE DO_PBL_MIX_2
!EOC
END MODULE vdiff_mod
