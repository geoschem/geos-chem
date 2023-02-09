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
MODULE Vdiff_Mod
!
! !USES:
!
  USE Error_Mod,     ONLY : Debug_Msg
  USE PhysConstants, ONLY : AIRMW, AVO, g0, Rd, Rv, Rdg0, VON_KARMAN
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Cleanup_Vdiff
  PUBLIC :: Do_Vdiff
  PUBLIC :: Init_Vdiff
  PUBLIC :: Max_PblHt_For_Vdiff
!
! !REMARKS:
!  The non-local PBL mixing routine VDIFF modifies the specific humidity,
!  (State_Met%SPHU) field.  Therefore, we must pass State_Met as an argument
!  to DO_PBL_MIX_2 and VDIFFDR with INTENT(INOUT).
!                                                                             .
! !REVISION HISTORY:
!  (1 ) This code is modified from mo_vdiff.F90 in MOZART-2.4. (lin, 5/14/09)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Physical constants
  REAL(fp), PARAMETER   :: cpair  = 1004.64_fp
  REAL(fp), PARAMETER   :: latvap = 2.5104e+06_fp
  REAL(fp), PARAMETER   :: rhoh2o = 1.0e+3_fp
  REAL(fp), PARAMETER   :: tfh2o  = 273.16_fp
  REAL(fp), PARAMETER   :: rair   = Rd
  REAL(fp), PARAMETER   :: rh2o   = Rv
  REAL(fp), PARAMETER   :: g      = G0
  REAL(fp), PARAMETER   :: gravit = g0
  REAL(fp), PARAMETER   :: zvir   = rh2o/rair - 1.0_fp
  REAL(fp), PARAMETER   :: cappa  = Rd/cpair
  REAL(fp), PARAMETER   :: r_g    = Rdg0

  ! PBL constants
  REAL(fp), PARAMETER   :: betam  = 15.0_fp  ! For wind gradient expression
  REAL(fp), PARAMETER   :: betas  = 5.0_fp   ! For surface layer gradient
  REAL(fp), PARAMETER   :: betah  = 15.0_fp  ! For temperature gradient
  REAL(fp), PARAMETER   :: fak    = 8.50_fp  ! For surface temperature excess
  REAL(fp), PARAMETER   :: fakn   = 7.20_fp  ! For turbulent prandtl number
  REAL(fp), PARAMETER   :: ricr   = 0.3_fp   ! For critical richardson number
  REAL(fp), PARAMETER   :: sffrac = 0.1_fp   ! For surface layer fraction of BL
  REAL(fp), PARAMETER   :: vk     = VON_KARMAN

  ! Derived constants
  REAL(fp), PARAMETER   :: ccon   = fak    * sffrac * vk
  REAL(fp), PARAMETER   :: binm   = betam  * sffrac
  REAL(fp), PARAMETER   :: binh   = betah  * sffrac
  REAL(fp), PARAMETER   :: onet   = 1.0_fp / 3.0_fp

  ! Options
  LOGICAL,  PARAMETER   :: divdiff = .TRUE.
  LOGICAL,  PARAMETER   :: arvdiff = .FALSE.
  LOGICAL,  PARAMETER   :: pblh_ar = .TRUE.
!
! !PRIVATE DATA MEMBERS:
!
  LOGICAL               :: prtDebug          ! Should we print debug info?
  INTEGER               :: nspcmix           ! # of species for mixing
  INTEGER               :: plev              ! # of levels
  INTEGER               :: plevp             ! # of level edges
  INTEGER               :: ntopfl            ! top level to which vertical
                                             !  diffusion is applied.
  INTEGER               :: npbl              ! max # of levels in pbl
  REAL(fp)              :: zkmin             ! minimum kneutral*f(ri)
  REAL(fp), ALLOCATABLE :: ml2(:)            ! mixing lengths squared
  REAL(fp), ALLOCATABLE :: qmincg(:)         ! min. constituent concentration

CONTAINS
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
                    sflx,       thp_arg,   pblh_arg,                         &
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
         cgq(plonl,plevp,nspcmix),&! countergrad term for constituent
         cgsh(plonl,plevp), &      ! countergrad term for sh
         kvf(plonl,plevp)          ! free atmosphere kv at interfaces
    real(fp) :: &
         potbar(plonl,plevp), &    ! pintm1(k)/(.5*(tm1(k)+tm1(k-1))
         tmp1(plonl), &            ! temporary storage
         dubot(plonl), &           ! lowest layer u change from stress
         dvbot(plonl), &           ! lowest layer v change from stress
         dtbot(plonl), &           ! lowest layer t change from heat flx
         dqbot(plonl,nspcmix), &   ! lowest layer q change from const flx
         dshbot(plonl), &          ! lowest layer sh change from wvflx
         thx(plonl,plev), &        ! temperature input + counter gradient
         thv(plonl,plev), &        ! virtual potential temperature
         qmx(plonl,plev,nspcmix), &! constituents input + counter grad
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
         cflx(plonl,nspcmix), &    ! surface constituent flux (kg/m2/s)
         wvflx(plonl)              ! water vapor flux (kg/m2/s)
    real(fp) :: &
         qp1(plonl,plev,nspcmix), &! moist, tracers after vert. diff
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

    real(fp) :: pblh(plonl)             ! boundary-layer height [m]

    real(fp) :: qp0(plonl,plev,nspcmix) ! To store initial concentration values

    real(fp) :: sum_qp0, sum_qp1        ! Jintai Lin 20180809

    !=================================================================
    ! vdiff begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Zero/initialize local variables for safety's sake
    indx     = 0
    adjust   = .FALSE.
    sum_qp0  = 0.0_fp
    sum_qp1  = 0.0_fp
    cah      = 0.0_fp
    cam      = 0.0_fp
    cch      = 0.0_fp
    ccm      = 0.0_fp
    cgh      = 0.0_fp
    cgq      = 0.0_fp
    cgs      = 0.0_fp
    cgsh     = 0.0_fp
    kvf      = 0.0_fp
    potbar   = 0.0_fp
    tmp1     = 0.0_fp
    dubot    = 0.0_fp
    dvbot    = 0.0_fp
    dtbot    = 0.0_fp
    dqbot    = 0.0_fp
    dshbot   = 0.0_fp
    thx      = 0.0_fp
    thv      = 0.0_fp
    qmx      = 0.0_fp
    shmx     = 0.0_fp
    zeh      = 0.0_fp
    zem      = 0.0_fp
    termh    = 0.0_fp
    termm    = 0.0_fp
    taux     = 0.0_fp
    tauy     = 0.0_fp
    ustar    = 0.0_fp
    qp0      = 0.0_fp
    um1      = uwnd(:,lat,:)
    vm1      = vwnd(:,lat,:)
    tm1      = tadv(:,lat,:)
    pmidm1   = pmid(:,lat,:)
    pintm1   = pint(:,lat,:)
    rpdel    = rpdel_arg(:,lat,:)
    rpdeli   = rpdeli_arg(:,lat,:)
    zm       = zm_arg(:,lat,:)
    shflx    = shflx_arg(:,lat)
    cflx     = sflx(:,lat,:)
    wvflx    = wvflx_arg(:,lat)
!ewl    qp1      = as2(:,lat,:,:)
!ewl    qp0      = as2(:,lat,:,:)
    shp1     = shp(:,lat,:)
    thp      = thp_arg(:,lat,:)
    kvh      = kvh_arg(:,lat,:)
    kvm      = kvm_arg(:,lat,:)
    tpert    = tpert_arg(:,lat)
    qpert    = qpert_arg(:,lat)
    cgs      = cgs_arg(:,lat,:)
    pblh     = pblh_arg(:,lat)

    !### Debug
    IF ( prtDebug .and. ip < 5 .and. lat < 5 ) &
         CALL DEBUG_MSG( '### VDIFF: vdiff begins' )

    IF (PRESENT(taux_arg )) taux  = taux_arg(:,lat)
    IF (PRESENT(tauy_arg )) tauy  = tauy_arg(:,lat)
    IF (PRESENT(ustar_arg)) ustar = ustar_arg(:,lat)

    ! Set initial species concentrations
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( M, I, L ) 
    DO M = 1, nspcmix
    DO I = 1, plonl
    DO L = 1, plev
       qp1(I,L,M) = State_Chm%Species(M)%Conc(I,lat,plev-L+1)
       qp0(I,L,M) = State_Chm%Species(M)%Conc(I,lat,plev-L+1)
    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

! resume...

!-----------------------------------------------------------------------
! 	... convert the surface fluxes to lowest level tendencies
!-----------------------------------------------------------------------
    rcpair = 1.0_fp/cpair
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
       kvf(i,plevp) = 0.0_fp
    end do
    do m = 1,nspcmix
       dqbot(:plonl,m) = cflx(:plonl,m)*tmp1(:plonl)
    end do

!      !### Debug
    IF ( prtDebug .and. ip < 5 .and. lat < 5 ) &
         CALL DEBUG_MSG( '### VDIFF: diffusion begins' )

!-----------------------------------------------------------------------
! 	... set the vertical diffusion coefficient above the top diffusion level
!-----------------------------------------------------------------------
    do k = 1,ntopfl
       kvf(:plonl,k) = 0.0_fp
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
          sstab = gravit*2.0_fp*(thv(i,k) - thv(i,k+1))/((thv(i,k) &
                 + thv(i,k+1))*dz)
!-----------------------------------------------------------------------
! 	... richardson number, stable and unstable modifying functions
!-----------------------------------------------------------------------
          rinub = sstab/dvdz2
          fstab = 1.0_fp/(1.0_fp + 10.0_fp*rinub*(1.0_fp &
                                 + 8.0_fp*rinub))
          funst = max( 1.0_fp - 18.0_fp*rinub,0.0_fp )
!-----------------------------------------------------------------------
! 	... select the appropriate function of the richardson number
!-----------------------------------------------------------------------
          if( rinub < 0.0_fp ) then
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
       do m = 1,nspcmix
          do i = 1,plonl
             qmx(i,k,m) = qp1(i,k,m)
          end do
       end do
    end do
    do k = 2,plev
       do i = 1,plonl
          potbar(i,k) = pintm1(i,k)/(0.5_fp*(tm1(i,k) + tm1(i,k-1)))
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
       do m = 1,nspcmix
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
    do m = 1,nspcmix
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
          if( shmx(i,k) < 1.0e-12_fp ) then
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
       cah(i,plev) = 0.0_fp
       cam(i,plev) = 0.0_fp
    end do

!-----------------------------------------------------------------------
! 	... calculate e(k) for heat & momentum vertical diffusion.  this term is
!           required in solution of tridiagonal matrix defined by implicit diffusion eqn.
!-----------------------------------------------------------------------
    do i = 1,plonl
       termh(i,ntopfl) = 1.0_fp/(1.0_fp + cah(i,ntopfl))
       termm(i,ntopfl) = 1.0_fp/(1.0_fp + cam(i,ntopfl))
       zeh(i,ntopfl)   = cah(i,ntopfl)*termh(i,ntopfl)
       zem(i,ntopfl)   = cam(i,ntopfl)*termm(i,ntopfl)
    end do
    do k = ntopfl+1,plev-1
       do i = 1,plonl
! Suspect that these lines lead to numerical instability
!         termh(i,k) = 1.e+0_fp/(1.e+0_fp + cah(i,k) + cch(i,k) &
!                     - cch(i,k)*zeh(i,k-1))
!         termm(i,k) = 1.e+0_fp/(1.e+0_fp + cam(i,k) + ccm(i,k) &
!                     - ccm(i,k)*zem(i,k-1))
          termh(i,k) =                                                       &
            1.0_fp / ( 1.0_fp + cah(i,k) + cch(i,k)*( 1.0_fp - zeh(i,k-1) ) )
          termm(i,k) =                                                       &
            1.0_fp / ( 1.0_fp + cam(i,k) + ccm(i,k)*( 1.0_fp - zem(i,k-1) ) )
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

    call qvdiff( nspcmix, qmx, dqbot, cch, zeh, &
	         termh, qp1, plonl )

!-----------------------------------------------------------------------
! 	... identify and correct constituents exceeding user defined bounds
!-----------------------------------------------------------------------
!      call qneg3( 'vdiff   ', lat, qp1, plonl )
! just use a simplified treatment
    where (qp1 < 0.0_fp)
       qp1 = 0.0_fp
    endwhere

!-----------------------------------------------------------------------
! Simple bug fix to ensure mass conservation - Jintai Lin 20180809
!   Without this fix, mass is almost but not completed conserved,
!   which is OK for full chemistry simulations but a big problem
!   for long lived species such as CH4 and CO2
!-----------------------------------------------------------------------
    DO M = 1, nspcmix
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
       shp1 = 0.0_fp
    endwhere

!-----------------------------------------------------------------------
! 	... diffuse potential temperature
!-----------------------------------------------------------------------
    call qvdiff( 1, thx, dtbot, cch, zeh, &
	         termh, thp, plonl )

    !Output values from local variables to arguments.(ccc, 11/17/09)
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

    ! Set species concentrations
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( M, I, L ) 
    DO M = 1, nspcmix
    DO I = 1, plonl
    DO L = 1, plev
       State_Chm%Species(M)%Conc(I,lat,L) = qp1(I,plev-L+1,M)
    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

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
         cflx(plonl,nspcmix), &     ! surface constituent flux (kg/m2/s)
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
         cgq(plonl,plevp,nspcmix), &! counter-gradient term for constituents
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
    real(fp), parameter :: tiny = 1.0e-20_fp  !e-36    ! lower bound for wind magnitude

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
         kqfs(plonl,nspcmix), &  ! sfc kinematic constituent flux [m/s]
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

    ! Zero/initialize local variables if they are not initialized below
    unstbl  = .FALSE.
    stblev  = .FALSE.
    unslev  = .FALSE.
    unssrf  = .FALSE.
    unsout  = .FALSE.
    check   = .FALSE.
    thvref  = 0.0_fp
    tkv     = 0.0_fp
    phiminv = 0.0_fp
    phihinv = 0.0_fp
    vvk     = 0.0_fp
    zm      = 0.0_fp
    zp      = 0.0_fp
    khfs    = 0.0_fp
    kqfs    = 0.0_fp
    kshfs   = 0.0_fp
    zmzp    = 0.0_fp
    rino    = 0.0_fp
    tlv     = 0.0_fp
    fak1    = 0.0_fp
    fak2    = 0.0_fp
    pblk    = 0.0_fp
    pr      = 0.0_fp
    zl      = 0.0_fp
    zzh     = 0.0_fp
    wstr    = 0.0_fp
    rrho    = 0.0_fp
    ustr    = 0.0_fp
    term    = 0.0_fp
    fac     = 0.0_fp
    pblmin  = 0.0_fp

!------------------------------------------------------------------------
! 	... compute kinematic surface fluxes
!------------------------------------------------------------------------
    do i = 1,plonl
       rrho(i)  = rair*t(i,plev)/pmid(i,plev)
       if (present(taux) .and. present(tauy)) then
          ustr     = sqrt( sqrt( taux(i)**2 + tauy(i)**2 )*rrho(i) )
          ustar(i) = max( ustr, 0.01_fp )
       endif
       khfs(i)  = shflx(i)*rrho(i)/cpair
       kshfs(i) = wvflx(i)*rrho(i)
    end do
    do m = 1,nspcmix
       kqfs(:plonl,m) = cflx(:plonl,m)*rrho(:plonl)
    end do

!------------------------------------------------------------------------
! 	... initialize output arrays with free atmosphere values
!------------------------------------------------------------------------
    do k = 1,plevp
       kvm(:,k)  = kvf(:,k)
       kvh(:,k)  = kvf(:,k)
       cgh(:,k)  = 0.0_fp
       cgsh(:,k) = 0.0_fp
       cgs(:,k)  = 0.0_fp
    end do
    do m = 1,nspcmix
       do k = 1,plevp
          cgq(:,k,m) = 0.0_fp
       end do
    end do

!------------------------------------------------------------------------
! 	... compute various arrays for use later:
!------------------------------------------------------------------------
    do i = 1,plonl
       thvsrf(i) = th(i,plev)*(1.0_fp + 0.61_fp*q(i,plev))
       heatv(i)  = khfs(i) + 0.61_fp*th(i,plev)*kshfs(i)
       wm(i)     = 0.0_fp
       therm(i)  = 0.0_fp
       qpert(i)  = 0.0_fp
       tpert(i)  = 0.0_fp
       fak3(i)   = 0.0_fp
       zh(i)     = 0.0_fp
       obklen(i) = -thvsrf(i)*ustar(i)**3 &
                   /(g*vk*(heatv(i) + sign( 1.0e-10_fp, heatv(i) )))
    end do

    if (pblh_ar) then  ! use archived PBLH

       do i = 1,plonl
          if( heatv(i) > 0.0_fp ) then
             unstbl(i) = .true.
          else
             unstbl(i) = .false.
          end if
          thvref(i) = th(i,plev)*(1.0_fp + 0.61_fp*q(i,plev))
       end do

    else ! use derived PBLH

!------------------------------------------------------------------------
! 	... define first a new factor fac=100 for use in richarson number
!           calculate virtual potential temperature first level
!           and initialize pbl height to z1
!------------------------------------------------------------------------
       fac = 100.0_fp
       do i = 1,plonl
          thvref(i) = th(i,plev)*(1.0_fp + 0.61_fp*q(i,plev))
          pblh(i)   = z(i,plev)
          check(i)  = .true.
!------------------------------------------------------------------------
! 	... initialization of lowest level ri number
!           (neglected in initial holtslag implementation)
!------------------------------------------------------------------------
          rino(i,plev) = 0.0_fp
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
                tkv = th(i,k)*(1.0_fp + 0.61_fp*q(i,k))
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
          if( heatv(i) > 0.0_fp ) then
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
             phiminv(i)   = (1.0_fp - binm*pblh(i)/obklen(i))**onet
             wm(i)        = ustar(i)*phiminv(i)
             therm(i)     = heatv(i)*fak/wm(i)
             rino(i,plev) = 0.0_fp
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
                tkv = th(i,k)*(1. + 0.61_fp*q(i,k))
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
          pblmin  = 700.0_fp*ustar(i)
          pblh(i) = max( pblh(i),pblmin )
       end do

    endif ! if pblh_ar

!------------------------------------------------------------------------
! 	... pblh is now available; do preparation for diffusivity calculation:
!------------------------------------------------------------------------
    do i = 1,plonl
       pblk(i) = 0.0_fp
       fak1(i) = ustar(i)*pblh(i)*vk
!------------------------------------------------------------------------
! 	... do additional preparation for unstable cases only, set temperature
!           and moisture perturbations depending on stability.
!------------------------------------------------------------------------
       if( unstbl(i) ) then
          phiminv(i) = (1.0_fp - binm*pblh(i)/obklen(i))**onet
          phihinv(i) = sqrt(1.0_fp - binh*pblh(i)/obklen(i))
          wm(i)      = ustar(i)*phiminv(i)
          fak2(i)    = wm(i)*pblh(i)*vk
          wstr(i)    = (heatv(i)*g*pblh(i)/thvref(i))**onet
          fak3(i)    = fakn*wstr(i)/wm(i)
          tpert(i)   = max( khfs(i)*fak/wm(i),0.0_fp )
          qpert(i)   = max( kshfs(i)*fak/wm(i),0.0_fp )
       else
          tpert(i)   = max( khfs(i)*fak/ustar(i),0.0_fp )
          qpert(i)   = max( kshfs(i)*fak/ustar(i),0.0_fp )
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
! NOTE: Do not test for floating-point equality, which due to roundoff
! may never occur. 
!          if( zkmin == 0.0_fp .and. zp(i) > pblh(i) ) then
          if ( ( .not. ABS( zkmin ) > 0.0_fp ) .and. ( zp(i) > pblh(i) ) ) then
             zp(i) = pblh(i)
          end if
          if( zm(i) < pblh(i) ) then
             zmzp = 0.5_fp*(zm(i) + zp(i))
             zh(i) = zmzp/pblh(i)
             zl(i) = zmzp/obklen(i)
             if( zh(i) <= 1.0_fp ) then
                zzh(i) = (1.0_fp - zh(i))**2
             else
                zzh(i) = 0.0_fp
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
             if( zl(i) <= 1.0_fp ) then
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
             term    = (1.0_fp - betam*zl(i))**onet
             pblk(i) = fak1(i)*zh(i)*zzh(i)*term
             pr(i)   = term/sqrt(1.0_fp - betah*zl(i))
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
       do m = 1,nspcmix
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
	             term, qp1, plonl, lat )
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
    integer, optional :: lat
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
         zfq(plonl,plev,nspcmix), & ! terms appear in soln of tri-diag sys
         tmp1d(plonl)               ! temporary workspace (1d array)
    integer :: &
         i, k, &               ! longitude,vertical indices
         m                     ! constituent index
    logical :: doPrt

    !=================================================================
    ! qvdiff begins here!
    !=================================================================
    doPrt = .FALSE.
    if ( present( lat ) ) doPrt = ( lat == 6 )

    ! Zero output arguments for safety's sake
    qp1   = 0.0_fp

    ! Zero local variables for safety's sake
    zfq   = 0.0_fp
    tmp1d = 0.0_fp

!-----------------------------------------------------------------------
! 	... calculate fq(k).  terms fq(k) and e(k) are required in solution of
!           tridiagonal matrix defined by implicit diffusion eqn.
!           note that only levels ntopfl through plev need be solved for.
!           no vertical diffusion is applied above this level
!-----------------------------------------------------------------------
    do m = 1,ncnst
       do i = 1,plonl
          zfq(i,ntopfl,m) = qm1(i,ntopfl,m) * term(i,ntopfl)
       end do
       do k = ntopfl+1,plev-1
          do i = 1,plonl
             zfq(i,k,m) = (qm1(i,k,m) + ( cc(i,k) * zfq(i,k-1,m) ) )         &
                        * term(i,k)
          end do
       end do
    end do
!-----------------------------------------------------------------------
! 	... bottom level: (includes  surface fluxes)
!-----------------------------------------------------------------------
    do i = 1,plonl
! Suspect that this line leads to numerical instability
!       tmp1d(i) = 1.0_fp /                                                   &
!                 (1.0_fp + cc(i,plev) - ( cc(i,plev) * ze(i,plev-1) ) )
       tmp1d(i) = 1.0_fp / ( 1.0_fp + cc(i,plev) * ( 1.0_fp - ze(i,plev-1) ) )
       ze(i,plev) = 0.0_fp
    end do
    do m = 1,ncnst
       do i = 1,plonl
          zfq(i,plev,m) =                                                    &
               (qm1(i,plev,m) + qflx(i,m) + cc(i,plev)*zfq(i,plev-1,m))      &
               * tmp1d(i)
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
             qp1(i,k,m) = zfq(i,k,m) + ( ze(i,k) * qp1(i,k+1,m) )
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
         cah(plonl,plev), &        ! -upper diag for heat and constituts
         cch(plonl,plev), &        ! -lower diag for heat and constits
         cgq(plonl,plevp,nspcmix),&! countergrad term for constituent
         potbar(plonl,plevp), &    ! pintm1(k)/(.5*(tm1(k)+tm1(k-1))
         tmp1(plonl), &            ! temporary storage
         tmp2, &                   ! temporary storage
         ztodtgor, &               ! ztodt*gravit/rair
         gorsq, &                  ! (gravit/rair)**2
         dqbot(plonl,nspcmix), &   ! lowest layer q change from const flx
         qmx(plonl,plev,nspcmix),& ! constituents input + counter grad
         zeh(plonl,plev), &        ! term in tri-diag. matrix system (t & q)
         termh(plonl,plev)         ! 1./(1. + cah(k) + cch(k) - cch(k)*zeh(k-1))
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
         cflx(plonl,nspcmix), &    ! surface constituent flux (kg/m2/s)
         kvh(plonl,plevp), &       ! coefficient for heat and tracers
         cgs(plonl,plevp)          ! counter-grad star (cg/flux)
    real(fp) :: &
         qp1(plonl,plev,nspcmix)   ! moist, tracers after vert. diff
    !=================================================================
    ! vdiffar begins here!
    !=================================================================

    ! Zero/initialize local variables for safety's sake
    indx   = 0
    ilogic = 0
    cah    = 0.0_fp
    cch    = 0.0_fp
    cgq    = 0.0_fp
    potbar = 0.0_fp
    tmp1   = 0.0_fp
    dqbot  = 0.0_fp
    qmx    = 0.0_fp
    zeh    = 0.0_fp
    termh  = 0.0_fp
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
    do m = 1,nspcmix
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
       do m = 1,nspcmix
          qmx(:,k,m) = qp1(:,k,m)
       end do
    end do
    do k = 2,plev
       potbar(:,k) = pintm1(:,k)/(0.5_fp*(tm1(:,k) + tm1(:,k-1)))
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
       do m = 1,nspcmix
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
    do m = 1,nspcmix
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
       cah(i,plev) = 0.0_fp
    end do
!-----------------------------------------------------------------------
! 	... calculate e(k) for heat vertical diffusion.  this term is
!           required in solution of tridiagonal matrix defined by implicit
!           diffusion eqn.
!-----------------------------------------------------------------------
    do i = 1,plonl
       termh(i,ntopfl) = 1.0_fp/(1.0_fp + cah(i,ntopfl))
       zeh(i,ntopfl) = cah(i,ntopfl)*termh(i,ntopfl)
    end do
    do k = ntopfl+1,plev-1
       do i = 1,plonl
          termh(i,k) = 1.0_fp/(1.0_fp + cah(i,k) + cch(i,k) &
                      - cch(i,k)*zeh(i,k-1))
          zeh(i,k) = cah(i,k)*termh(i,k)
       end do
    end do
!-----------------------------------------------------------------------
! 	... diffuse constituents
!-----------------------------------------------------------------------
    call qvdiff( nspcmix, qmx, dqbot, cch, zeh, &
	         termh, qp1, plonl )
!-----------------------------------------------------------------------
! 	... identify and correct constituents exceeding user defined bounds
!-----------------------------------------------------------------------
!      call qneg3( 'vdiff   ', lat, qp1(1,1,1), plonl )
!     simplified treatment
    where (qp1 < 0.0_fp)
       qp1 = 0.0_fp
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
         cflx(plonl,nspcmix), &  ! surface constituent flux (kg/m2/s)
         cgs(plonl,plevp)        ! counter-gradient star (cg/flux)
!
! !OUTPUT PARAMETERS:
!
    real(fp), intent(out) :: &
         cgq(plonl,plevp,nspcmix) ! counter-gradient term for constituents
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
         rrho(plonl), &       ! 1./bottom level density
         kqfs(plonl,nspcmix)  ! sfc kinematic constituent flux [m/s]

    !=================================================================
    ! pbldifar begins here!
    !=================================================================

    ! Zero output arguments for safety's sake
    cgq  = 0.0_fp

    ! Zero local variables for safety's sake
    rrho = 0.0_fp
    kqfs = 0.0_fp

!------------------------------------------------------------------------
! 	... compute kinematic surface fluxes
!------------------------------------------------------------------------
    rrho(:) = rair*t(:,plev)/pmid(:,plev)
    do m = 1,nspcmix
       kqfs(:,m) = cflx(:,m)*rrho(:)
    end do
!------------------------------------------------------------------------
! 	... initialize output arrays with free atmosphere values
!------------------------------------------------------------------------
    do m = 1,nspcmix
       do k = 1,plevp
          cgq(:,k,m) = 0.0_fp
       end do
    end do
!------------------------------------------------------------------------
! 	... compute the counter-gradient terms:
!------------------------------------------------------------------------
    do k = plev,plev-npbl+2,-1
       do m = 1,nspcmix
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
! !IROUTINE: Init_Vdiff
!
! !DESCRIPTION: Initializes fields used by the non-local BL mixing scheme.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Vdiff( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chem State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Init_Vdiff begins here!
    !=================================================================

    ! Assume success
    RC       = GC_SUCCESS

    ! Exit if this is a dry-run simulation
    IF ( Input_Opt%DryRun ) RETURN

    ! Initialize
    plev     = State_Grid%NZ           ! # of levels
    plevp    = State_Grid%NZ + 1       ! # of level edges
    nspcmix  = State_Chm%nAdvect       ! # of species
    zkmin    = 0.01_fp                 ! = minimum k = kneutral * f(ri)

    ! Set the square of the mixing lengths
    ALLOCATE( ml2(plevp), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod:ML2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ml2(1     ) = 0.0_fp
    ml2(2:plev) = 900.0_fp ! = (30.0_fp)**2
    ml2(plevp ) = 0.0_fp

    ! Set the minimum mixing ratio for the counter-gradient term.
    ! normally this should be the same as qmin.
    ALLOCATE( qmincg(nspcmix), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod:QMINCG', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    qmincg = 0.0_fp

  END SUBROUTINE Init_Vdiff
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Max_PblHt_For_Vdiff
!
! !DESCRIPTION: Computes the maximum boundary layer height variables
!  for the non-local mixing.  This was rewritten to avoid assuming the
!  a specific grid configuration.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Max_PblHt_For_Vdiff( Input_Opt, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt
    TYPE(GrdState), INTENT(IN)  :: State_Grid
    TYPE(MetState), INTENT(IN)  :: State_Met
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC
!
! !REMARKS:
!  This routine contains code that was originally in Init_Vdiff.  But it
!  has to be separated so that it can be called after the initial met fields
!  are read from disk.  This allows us to use the surface pressure field
!  instead of referencing the Ap and Bp parameters directly.
!
! !REVISION HISTORY:
!  18 May 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: pbl_press = 400.0_fp   ! pressure cap for pbl (Pa)
!
! !LOCAL VARIABLES
!
    ! Scalars
    INTEGER             :: K
    REAL(fp)            :: cells_per_layer

    ! Arrays
    REAL(fp)            :: ref_pmid(State_Grid%NZ)

    ! Strings
    CHARACTER(LEN=255)  :: errMsg
    CHARACTER(LEN=255)  :: thisLoc

    !===================================================================
    ! Max_PblHt_for_Vdiff begins here!
    !===================================================================

    ! Assume success
    RC              = GC_SUCCESS

    ! Exit if the non-local PBL mixing is not being used
    IF ( .not. ( Input_Opt%LTURB .and. Input_Opt%LNLPBL ) ) RETURN

    ! Initialize
    ref_pmid        = 0.0_fp
    cells_per_layer = DBLE( State_Grid%NX * State_Grid%NY )
    errMsg          = ''
    thisLoc         = &
     ' -> at Max_PblHt_for_Vdiff (in module GeosCore/vdiff_mod.F90)'

    ! Now use the average initial surface pressure per layer instead of
    ! having to reference the Ap and Bp.  This should make it easier
    ! to interface to external models such as CESM.
    DO k = 1, plev
       ref_pmid(plev-k+1) = SUM( State_Met%PMid(:,:,K) ) / cells_per_layer
    ENDDO

    ! Derived constants
    ! ntopfl = top level to which v-diff is applied
    ! npbl   = max number of levels (from bottom) in pbl
    DO k = plev,1,-1
       IF( ref_pmid(k) < pbl_press ) then
          EXIT
       ENDIF
    ENDDO

    ! Compute the number of PBL levels
    npbl = MAX( 1, plev - k )
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE(6,*) 'Init_Vdiff: pbl height will be limited to bottom ',npbl,  &
            ' model levels.'
       WRITE(6,*) 'Top is ',ref_pmid(plevp-npbl),' hpa'
    ENDIF

    ! Set the ntopfl
    IF ( plev == 1 ) THEN
       ntopfl = 0
    ELSE
       ntopfl = 1
    ENDIF

  END SUBROUTINE Max_PblHt_For_Vdiff
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
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PBL_MIX_MOD,        ONLY : COMPUTE_PBL_HEIGHT
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CONV, GET_TS_EMIS, GET_TS_CHEM
    USE OCEAN_MERCURY_MOD,  ONLY : Fg !hma
    USE OCEAN_MERCURY_MOD,  ONLY : OMMFP => Fp
    USE OCEAN_MERCURY_MOD,  ONLY : LHg2HalfAerosol !cdh

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
!  (2) As of July 2016, we assume that all of the advected species are listed
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
    ! Scalars
    INTEGER             :: I, J, L, N, NA, nAdvect, EC
    REAL(fp)            :: dtime

    ! Arrays
    REAL(fp), TARGET    :: tpert (State_Grid%NX,State_Grid%NY)
    REAL(fp), TARGET    :: qpert (State_Grid%NX,State_Grid%NY)
    REAL(fp), TARGET    :: shflx (State_Grid%NX,State_Grid%NY)
    REAL(fp), TARGET    :: pmid  (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp), TARGET    :: rpdel (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp), TARGET    :: rpdeli(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp), TARGET    :: thp   (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp), TARGET    :: t1    (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp), TARGET    :: zm    (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp), TARGET    :: cgs   (State_Grid%NX, State_Grid%NY,State_Grid%NZ+1)
    REAL(fp), TARGET    :: pint  (State_Grid%NX, State_Grid%NY,State_Grid%NZ+1)
    REAL(fp), TARGET    :: kvh   (State_Grid%NX, State_Grid%NY,State_Grid%NZ+1)
    REAL(fp), TARGET    :: kvm   (State_Grid%NX, State_Grid%NY,State_Grid%NZ+1)

    ! Pointers
    REAL(fp), POINTER   :: p_pblh  (:,:    )
    REAL(fp), POINTER   :: p_sflx  (:,:,:  )
    REAL(fp), POINTER   :: p_um1   (:,:,:  )
    REAL(fp), POINTER   :: p_vm1   (:,:,:  )
    REAL(fp), POINTER   :: p_tadv  (:,:,:  )
    REAL(fp), POINTER   :: p_hflux (:,:    )
    REAL(fp), POINTER   :: p_ustar (:,:    )
    REAL(fp), POINTER   :: p_pmid  (:,:,:  )
    REAL(fp), POINTER   :: p_pint  (:,:,:  )
    REAL(fp), POINTER   :: p_rpdel (:,:,:  )
    REAL(fp), POINTER   :: p_rpdeli(:,:,:  )
    REAL(fp), POINTER   :: p_zm    (:,:,:  )
    REAL(fp), POINTER   :: p_thp   (:,:,:  )
    REAL(fp), POINTER   :: p_kvh   (:,:,:  )
    REAL(fp), POINTER   :: p_kvm   (:,:,:  )
    REAL(fp), POINTER   :: p_cgs   (:,:,:  )
    REAL(fp), POINTER   :: p_shp   (:,:,:  )
    REAL(fp), POINTER   :: p_t1    (:,:,:  )

    ! For error trapping
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: p0 = 1.e+5_fp

    !=================================================================
    ! Vdiffdr begins here!
    !=================================================================

    !### Debug info
    IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: VDIFFDR begins' )

    ! Initialize
    RC       =  GC_SUCCESS
    nAdvect  =  State_Chm%nAdvect
    ErrMsg   = ''
    ThisLoc  = ' -> at VDIFF (in module GeosCore/vdiff_mod.F90)'
    pmid     =  0.0_fp
    rpdel    =  0.0_fp
    rpdeli   =  0.0_fp
    zm       =  0.0_fp
    pint     =  0.0_fp
    cgs      =  0.0_fp
    kvh      =  0.0_fp
    kvm      =  0.0_fp
    tpert    =  0.0_fp
    qpert    =  0.0_fp
    thp      =  0.0_fp
    shflx    =  0.0_fp
    t1       =  0.0_fp
    dtime    =  GET_TS_CONV()            ! second
    shflx    =  State_Met%EFLUX / latvap ! latent heat -> water vapor flux

!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( I, J, L )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       DO L = 1, State_Grid%NZ

          ! Convert PMID and PEDGE from hPa to Pa
          pmid(I,J,L) = State_Met%PMID(I,J,L)  * 100.0_fp
          pint(I,J,L) = State_Met%PEDGE(I,J,L) * 100.0_fp

          ! Potential temperature [K]
          thp(I,J,L) = State_Met%T(I,J,L)*(p0/pmid(I,J,L))**cappa
       ENDDO

       ! PEDGE at the top of the atmosphere
       pint(I,J,State_Grid%NZ+1) = State_Met%PEDGE(I,J,State_Grid%NZ+1)

       ! Corrected calculation of zm.
       ! Box height calculation now uses virtual temperature.
       ! Therefore, use virtual temperature in hypsometric equation.
       ! (ewl, 3/3/15)
       do L = 1, State_Grid%NZ
          zm(I,J,L) = SUM( State_Met%BXHEIGHT(I,J,1:L))                     &
                    - log( pmid(I,J,L)/pint(I,J,L+1) )                      &
                    * r_g * State_Met%TV(I,J,L)

          rpdel(I,J,L) = 1.0_fp / (pint(I,J,L) - pint(I,J,L+1))
       enddo

       !rpdeli(I,J,1) = 1.e+0_fp / (PS(I,J) - pmid(I,J,1))
       rpdeli(I,J,1) = 0.0_fp ! follow mozart setup (shown in mo_physlic.F90)
       do L = 2, State_Grid%NZ
          rpdeli(I,J,L) = 1.0_fp / (pmid(I,J,L-1) - pmid(I,J,L))
       enddo

    enddo
    enddo
!$OMP END PARALLEL DO

    !### Debug
    IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: after emis. and depdrp' )

    !--------------------------------------------------------------------
    ! Now use pointers to flip arrays in the vertical (bmy, 6/22/15)
    !--------------------------------------------------------------------
    p_sflx   => State_Chm%SurfaceFlux
    p_hflux  => State_Met%HFLUX
    p_ustar  => State_Met%USTAR
    p_pblh   => State_Met%PBL_TOP_m
    p_um1    => State_Met%U        (:, :, State_Grid%NZ  :1:-1           )
    p_vm1    => State_Met%V        (:, :, State_Grid%NZ  :1:-1           )
    p_tadv   => State_Met%T        (:, :, State_Grid%NZ  :1:-1           )
    p_pmid   => pmid               (:, :, State_Grid%NZ  :1:-1           )
    p_rpdel  => rpdel              (:, :, State_Grid%NZ  :1:-1           )
    p_rpdeli => rpdeli             (:, :, State_Grid%NZ  :1:-1           )
    p_zm     => zm                 (:, :, State_Grid%NZ  :1:-1           )
    p_thp    => thp                (:, :, State_Grid%NZ  :1:-1           )
    p_shp    => State_Met%SPHU     (:, :, State_Grid%NZ  :1:-1           )
    p_pint   => pint               (:, :, State_Grid%NZ+1:1:-1           )
    p_kvh    => kvh                (:, :, State_Grid%NZ+1:1:-1           )
    p_kvm    => kvm                (:, :, State_Grid%NZ+1:1:-1           )
    p_cgs    => cgs                (:, :, State_Grid%NZ+1:1:-1           )

    ! Convert v/v -> m/m (i.e., kg/kg)
    DO NA = 1, nAdvect
       State_Chm%Species(NA)%Conc(:,:,:) = State_Chm%Species(NA)%Conc(:,:,:)                                     &
                       / ( AIRMW / State_Chm%SpcData(NA)%Info%MW_g )
    ENDDO

    ! Convert g/kg -> kg/kg
    p_shp              =  p_shp * 1.e-3_fp

    !### Debug
    IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: before vdiff' )

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( J, EC  )
    DO J = 1, State_Grid%NY
       CALL Vdiff( J,                 1,         p_um1,     p_vm1,           &
                   p_tadv,            p_pmid,    p_pint,    p_rpdel,         &
                   p_rpdeli,          dtime,     p_zm,      p_hflux,         &
                   p_sflx,            p_thp,     p_pblh,                     &
                   p_kvh,             p_kvm,     tpert,     qpert,           &
                   p_cgs,             p_shp,     shflx,     State_Grid%NX,   &
                   Input_Opt,         State_Met, State_Chm, State_Diag,      &
                   ustar_arg=p_ustar, RC=EC                                 )
    ENDDO
    !$OMP END PARALLEL DO

    !### Debug
    IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: after vdiff' )

    ! Convert kg/kg -> v/v
    DO NA = 1, nAdvect
       State_Chm%Species(NA)%Conc(:,:,:) = State_Chm%Species(NA)%Conc(:,:,:)                                     &
                       * ( AIRMW / State_Chm%SpcData(NA)%Info%MW_g )
    ENDDO

    ! Convert kg/kg -> g/kg
    p_shp    = p_shp * 1.0e+3_fp

    ! Nullify pointers
    p_sflx   => NULL()
    p_um1    => NULL()
    p_vm1    => NULL()
    p_tadv   => NULL()
    p_hflux  => NULL()
    p_ustar  => NULL()
    p_pmid   => NULL()
    p_rpdel  => NULL()
    p_rpdeli => NULL()
    p_zm     => NULL()
    p_thp    => NULL()
    p_shp    => NULL()
    p_pint   => NULL()
    p_kvh    => NULL()
    p_kvm    => NULL()
    p_cgs    => NULL()

    !### Debug
    IF ( prtDebug ) CALL DEBUG_MSG( '### VDIFFDR: VDIFFDR finished' )

  END SUBROUTINE VDIFFDR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Vdiff
!
! !DESCRIPTION: Subroutine DO\_PBL\_MIX\_2 is the driver routine for planetary
!  boundary layer mixing. The PBL layer height and related quantities are
!  always computed.   Mixing of tracers underneath the PBL top is toggled
!  by the DO\_TURBDAY switch.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Vdiff( Input_Opt,  State_Chm, State_Diag,                    &
                       State_Grid, State_Met, RC                            )
!
! !USES:
!
    USE Calc_Met_Mod,       ONLY : AirQnt
    USE Diagnostics_Mod,    ONLY : Compute_Budget_Diagnostics
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : ITS_TIME_FOR_EMIS
    USE Time_Mod,           ONLY : Get_Ts_Dyn
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units

    IMPLICIT NONE
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
! !REVISION HISTORY:
!  11 Feb 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: prtDebug
    INTEGER            :: TS_Dyn
    REAL(fp)           :: DT_Dyn

    ! Strings
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: errMsg
    CHARACTER(LEN=255) :: thisLoc

    !=======================================================================
    ! DO_VDIFF begins here!
    !=======================================================================

    ! Initialize
    RC       =  GC_SUCCESS                                ! Assume success
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) ! Print debug output?
    ErrMsg   = ''
    ThisLoc  = ' -> at DO_PBL_MIX_2 (in module GeosCore/vdiff_mod.F90)'

    !=======================================================================
    ! Non-local PBL mixing budget diagnostics - Part 1 of 2
    !
    ! WARNING: The mixing budget diagnostic includes the application
    ! of species tendencies (emissions fluxes and dry deposition
    ! rates) below the PBL when using non-local PBL mixing. This is
    ! done for all species with defined emissions / dry deposition
    ! rates, including dust. These tendencies below the PBL are
    ! therefore not included in the emissions/dry deposition budget
    ! diagnostics when using non-local PBL mixing. (ewl, 9/26/18)
    !=======================================================================
    IF ( State_Diag%Archive_BudgetMixing ) THEN

       ! Get initial column masses (full, trop, PBL)
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = State_Diag%Archive_BudgetMixingFull,               &
            diagFull    = NULL(),                                            &
            mapDataFull = State_Diag%Map_BudgetMixingFull,                   &
            isTrop      = State_Diag%Archive_BudgetMixingTrop,               &
            diagTrop    = NULL(),                                            &
            mapDataTrop = State_Diag%Map_BudgetMixingTrop,                   &
            isPBL       = State_Diag%Archive_BudgetMixingPBL,                &
            diagPBL     = NULL(),                                            &
            mapDataPBL  = State_Diag%Map_BudgetMixingPBL,                    &
            colMass     = State_Diag%BudgetColumnMass,                       &
            before_op   = .TRUE.,                                            &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Non-local mixing budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Unit conversion #1
    !=======================================================================

    ! Convert species concentration to v/v dry
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid,                &
                            State_Met, 'v/v dry', RC,                        &
                            OrigUnit=OrigUnit                               )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "Convert_Spc_Units" (to v/v dry)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! PBL mixing
    !=======================================================================

    ! Do non-local PBL mixing
    CALL VDIFFDR( Input_Opt,  State_Chm, State_Diag,                         &
                  State_Grid, State_Met, RC                                 )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "VDIFFDR"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Debug print
    IF( prtDebug ) THEN
       CALL DEBUG_MSG( '### DO_PBL_MIX_2: after VDIFFDR' )
    ENDIF

    !=======================================================================
    ! Update airmass etc. and mixing ratios
    !=======================================================================

    ! Update air quantities and species concentrations with updated
    ! specific humidity (ewl, 10/28/15)
    !
    ! NOTE: Prior to October 2015, air quantities were not updated
    ! with specific humidity modified in VDIFFDR at this point in
    ! the model
    CALL AirQnt( Input_Opt, State_Chm, State_Grid,                           &
                 State_Met, RC,        Update_Mixing_Ratio=.TRUE.           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "AIRQNT"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Debug print
    IF( prtDebug ) THEN
       CALL DEBUG_MSG( '### DO_PBL_MIX_2: after AIRQNT' )
    ENDIF

    !=======================================================================
    ! Unit conversion #2
    !=======================================================================

    ! Convert species back to the original units
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid,                &
                            State_Met, OrigUnit,  RC                        )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "Convert_Spc_Units"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Non-local PBL mixing budget diagnostics - Part 2 of 2
    !=======================================================================
    IF ( State_Diag%Archive_BudgetMixing ) THEN

       ! Get dynamics timestep [s]
       TS_Dyn = Get_Ts_Dyn()
       DT_Dyn = DBLE( Ts_Dyn )

       ! Compute change in column masses (after mixing - before mixing)
       ! and store in diagnostic arrays.  Units are [kg/s].
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = State_Diag%Archive_BudgetMixingFull,               &
            diagFull    = State_Diag%BudgetMixingFull,                       &
            mapDataFull = State_Diag%Map_BudgetMixingFull,                   &
            isTrop      = State_Diag%Archive_BudgetMixingTrop,               &
            diagTrop    = State_Diag%BudgetMixingTrop,                       &
            mapDataTrop = State_Diag%Map_BudgetMixingTrop,                   &
            isPBL       = State_Diag%Archive_BudgetMixingPBL,                &
            diagPBL     = State_Diag%BudgetMixingPBL,                        &
            mapDataPBL  = State_Diag%Map_BudgetMixingPBL,                    &
            colMass     = State_Diag%BudgetColumnMass,                       &
            timeStep    = DT_Dyn,                                            &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Non-local mixing budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Do_Vdiff
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Vdiff
!
! !DESCRIPTION: Deallocates module arrays
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_VDiff( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialize
    RC = GC_SUCCESS

    ! Deallocate arrays
    IF ( ALLOCATED( ml2 ) ) THEN
       DEALLOCATE( ml2, STAT=RC )
       CALL GC_CheckVar( 'vdiff_mod.F90:ML2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( qmincg ) ) THEN
       DEALLOCATE( qmincg, STAT=RC )
       CALL GC_CheckVar( 'vdiff_mod.F90:QMINCG', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Cleanup_Vdiff
!EOC
END MODULE Vdiff_Mod
