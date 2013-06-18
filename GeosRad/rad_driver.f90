!------------------------------------------------------------------------------
!             Atmospheric and Environmental Research (AER Corp.)
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: 
!
!\\
!\\
! !INTERFACE:
!
SUBROUTINE rad_driver( IIPAR, JJPAR,   NLAY, NCOL, ICLD_GC,ISEED,T,TSKIN,  &
                       AVGW,  SUNCOS,  DOY,   PEDGE, PCENTER,  &
                       O3VMR, CH4VMR,N2OVMR,&
                       cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
                       TAUCLD,CLDFR,CLIQWP,CICEWP, &
                       RELIQ,REICE,&
                       ASDIR,ASDIF,ALDIR,ALDIF, &
                       EMIS,LLWRAD,LSWRAD,    &
                       TAUAER_LW,TAUAER_SW,SSAAER,ASMAER,        &
                       lw_uflux,lw_dflux,sw_uflux,sw_dflux, &
                       lw_ufluxc,lw_dfluxc,sw_ufluxc,sw_dfluxc) 
!
! !USES:

        use parkind, only : im=>kind_im, rb=>kind_rb
        use rrlw_con, only: gascon, avogad
        use parrrtm, only : nbndlw, ngptlw
        use parrrsw, only : nbndsw, ngptsw,naerec
        use rrtmg_lw_rad, only : rrtmg_lw
        use rrtmg_sw_rad, only : rrtmg_sw
        use mcica_subcol_gen_lw, only : mcica_subcol_lw
        use mcica_subcol_gen_sw, only : mcica_subcol_sw
!        USE LOGICAL_MOD,  ONLY : LLWRAD,LSWRAD

        implicit none

!

!
! !INPUT PARAMETERS: 
!
  ! Scalars
  REAL*8,  INTENT(IN) :: DOY
  INTEGER, INTENT(IN) :: IIPAR
  INTEGER, INTENT(IN) :: JJPAR
  INTEGER, INTENT(IN) :: NCOL      !! NCOL = IIPAR*JJPAR
  INTEGER, INTENT(IN) :: NLAY      !! NLAY = LLPAR
  INTEGER, INTENT(IN) :: icld_gc
  LOGICAL, INTENT(IN) :: LLWRAD    !from LOGICAL_MOD so we know
  LOGICAL, INTENT(IN) :: LSWRAD    !whether to do LW and SW
  INTEGER, INTENT(IN) :: iseed

  ! 2-D Arrays
  REAL*8,  INTENT(IN) :: TSKIN  (IIPAR,JJPAR      )
  REAL*8,  INTENT(IN) :: SUNCOS (IIPAR*JJPAR      )
  real (kind=rb), intent(in)   :: o3vmr(ncol,nlay)
  real (kind=rb), intent(in)   :: ch4vmr(ncol,nlay)
  real (kind=rb), intent(inout):: n2ovmr(ncol,nlay)
  real (kind=rb), intent(in)  :: cfc11vmr(ncol,nlay)
  real (kind=rb), intent(in) :: cfc12vmr(ncol,nlay)
  real (kind=rb), intent(in)  :: cfc22vmr(ncol,nlay)
  real (kind=rb), intent(in)  :: ccl4vmr(ncol,nlay)
  real (kind=rb), intent(in)   :: taucld(ncol,nlay)
  real (kind=rb), intent(in)   :: cldfr(ncol,nlay)
 !setting as inout so can fix for testing DAR
  real (kind=rb), intent(inout)   :: asdir(ncol)
  real (kind=rb), intent(inout)   :: asdif(ncol)
  real (kind=rb), intent(inout)   :: aldir(ncol)
  real (kind=rb), intent(inout)   :: aldif(ncol)
  real (kind=rb), intent(inout)   :: emis(ncol,nbndlw)
  

! Common cloud variables
  real (kind=rb), intent(in)  :: cicewp(ncol,nlay)
  real (kind=rb), intent(in) :: cliqwp(ncol,nlay)
  real (kind=rb), intent(in) :: reice(ncol,nlay)
  real (kind=rb), intent(in) :: reliq(ncol,nlay)

! Aerosol variables
!  LW
  real (kind=rb),intent(inout)  :: tauaer_lw(ncol,nlay,nbndlw)

! SW
   real(kind=rb),intent(inout)   :: tauaer_sw(ncol,nlay,nbndsw)  ! Aerosol optical depth (iaer=10 only)
                                                      !    Dimensions: (ncol,nlay,nbndsw)
                                                      ! (non-delta scaled)      
   real(kind=rb),intent(inout)   :: ssaaer(ncol,nlay,nbndsw)     ! Aerosol single scattering albedo (iaer=10 only)
                                                      !    Dimensions: (ncol,nlay,nbndsw)
                                                      ! (non-delta scaled)      

   real(kind=rb),intent(inout) :: asmaer(ncol,nlay,nbndsw)       ! Aerosol asymmetry parameter (iaer=10 only)
                                                      !    Dimensions: !    (ncol,nlay,nbndsw)
  ! 3-D arrays
  REAL*8,  INTENT(IN) :: AVGW   (IIPAR,JJPAR,NLAY)
  REAL*8,  INTENT(IN) :: PEDGE  (IIPAR,JJPAR,NLAY)
  REAL*8,  INTENT(IN) :: PCENTER(IIPAR,JJPAR,NLAY)
  REAL*8,  INTENT(IN) :: T      (IIPAR,JJPAR,NLAY)
!

! !RETURN VALUE:
! !OUTPUT PARAMETERS:
  REAL*8, INTENT(INOUT)    :: LW_UFLUX(ncol,nlay+1)
  REAL*8, INTENT(INOUT)    :: LW_DFLUX(ncol,nlay+1)
  REAL*8, INTENT(INOUT)    :: SW_UFLUX(ncol,nlay+1)
  REAL*8, INTENT(INOUT)    :: SW_DFLUX(ncol,nlay+1)
  REAL*8, INTENT(INOUT)    :: LW_UFLUXC(ncol,nlay+1)
  REAL*8, INTENT(INOUT)    :: LW_DFLUXC(ncol,nlay+1)
  REAL*8, INTENT(INOUT)    :: SW_UFLUXC(ncol,nlay+1)
  REAL*8, INTENT(INOUT)    :: SW_DFLUXC(ncol,nlay+1)

!

!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  01 Oct 1995 - R. Yantosca - Initial version
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!


  ! Flags and dimensions
  integer (kind=im) :: icld, idrv
  integer (kind=im) :: inflglw, iceflglw,liqflglw
  integer (kind=im) :: inflgsw, iceflgsw,liqflgsw


! Profile variables
  real (kind=rb)  :: play(ncol,nlay)
  real (kind=rb)  :: plev(ncol,nlay+1)
  real (kind=rb)  :: tlay(ncol,nlay)
  real (kind=rb)  :: tlev(ncol,nlay+1)
  real (kind=rb)  :: h2ovmr(ncol,nlay)
  real (kind=rb)  :: co2vmr(ncol,nlay)
  real (kind=rb)  :: o2vmr(ncol,nlay)
  
  
! LW Surface variables
  real (kind=rb)  :: tsfc(ncol)

!  LW Cloud variables
  real (kind=rb)  :: taucld_lw(nbndlw,ncol,nlay)

! SW Surface variables - NOW PASSED FROM GEOSCHEM
!  real(kind=rb)  :: asdir(ncol)          ! UV/vis surface albedo direct rad
                                         !    Dimensions: (ncol)
!  real(kind=rb) :: aldir(ncol)           ! Near-IR surface albedo direct rad
                                         !    Dimensions: (ncol)             
!  real(kind=rb) :: asdif(ncol)           ! UV/vis surface albedo: diffuse rad
                                         !    Dimensions: (ncol)
!  real(kind=rb) :: aldif(ncol)           ! Near-IR surface albedo: diffuse rad
                                         !    Dimensions: (ncol)
! SW solar variables
  integer(kind=im) :: dyofyr               ! Day of the year (used to get Earth/Sun
                                               !  distance if adjflx not provided)
  real(kind=rb)  :: adjes=1.0              ! Flux adjustment for Earth/Sun distance
  real(kind=rb)  :: coszen(ncol)           ! Cosine of solar zenith angle
                                               !    Dimensions: (ncol)
  real(kind=rb)  :: scon=1368.22           ! Solar constant (W/m2)

! SW cloud variables

   real(kind=rb)  :: taucld_sw(nbndsw,ncol,nlay)   ! In-cloud optical depth
                                                       !    Dimensions: (nbndsw,ncol,nlay)
   real(kind=rb)  :: ssacld(nbndsw,ncol,nlay)      ! In-cloud single scattering albedo
                                                      !    Dimensions: (nbndsw,ncol,nlay)
   real(kind=rb)  :: asmcld(nbndsw,ncol,nlay)      ! In-cloud asymmetry parameter
                                                      !    Dimensions: (nbndsw,ncol,nlay)
   real(kind=rb)  :: fsfcld(nbndsw,ncol,nlay)      ! In-cloud forward scattering fraction
                                                      !    Dimensions: (nbndsw,ncol,nlay)

   real(kind=rb) :: ecaer(ncol,nlay,naerec)        ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                                      !    Dimensions: !    (ncol,nlay,naerec)


! Longwave Flux variables
  real(kind=rb)  :: uflx(ncol,nlay+1)         ! Total sky longwave upward flux (W/m2)
                                                      !    Dimensions:
                                                      !    (ncol,nlay+1)

  real(kind=rb) :: dflx(ncol,nlay+1) ! Total sky longwave downward flux (W/m2)
                                                      !    Dimensions:
                                                      !    (ncol,nlay+1)
  real(kind=rb) :: hr(ncol,nlay)   ! Total sky longwave radiative heating rate (K/d)
                                                      !    Dimensions:
                                                      !    (ncol,nlay)
  real(kind=rb) :: uflxc(ncol,nlay+1)! Clear sky longwave upward flux (W/m2)
                                                      !    Dimensions:
                                                      !    (ncol,nlay+1)
  real(kind=rb) :: dflxc(ncol,nlay+1)! Clear sky longwave downward flux (W/m2)
                                                      !    Dimensions:
                                                      !    (ncol,nlay+1)
  real(kind=rb) :: hrc(ncol,nlay) ! Clear sky longwave radiative heating rate (K/d)
                                                      !    Dimensions:
                                                      !    (ncol,nlay)

!- Optional Output
  real(kind=rb) :: duflx_dt(ncol,nlay)
                                                      ! change in upward  longwave flux  (w/m2/k)
                                                      ! with respect to  surface  temperature
                                                      !    Dimensions:     (ncol,nlay)
  real(kind=rb) :: duflxc_dt(ncol,nlay)
                                                      ! change in clear  sky upward  longwave flux  (w/m2/k)
                                                      ! with respect to  surface  temperature
                                                      !    Dimensions:
                                                      !    (ncol,nlay)

! Shortwave flux variables
! ----- Output -----

  real(kind=rb)  :: swuflx(ncol,nlay+1)       ! Total sky shortwave upward flux (W/m2)
                                                      !    Dimensions:
                                                      !    (ncol,nlay+1)
  real(kind=rb) :: swdflx(ncol,nlay+1)       ! Total sky shortwave downward flux (W/m2)
                                                      !    Dimensions:
                                                      !    (ncol,nlay+1)
  real(kind=rb) :: swhr(ncol,nlay)         ! Total sky shortwave radiative heating rate (K/d)
                                                      !    Dimensions:
                                                      !    (ncol,nlay)
  real(kind=rb) :: swuflxc(ncol,nlay+1)      ! Clear sky shortwave upward flux (W/m2)
                                                      !    Dimensions:
                                                      !    (ncol,nlay+1)
  real(kind=rb) :: swdflxc(ncol,nlay+1)      ! Clear sky shortwave downward flux (W/m2)
                                                      !    Dimensions:
                                                      !    (ncol,nlay+1)
  real(kind=rb) :: swhrc(ncol,nlay)        ! Clear sky shortwave radiative heating rate (K/d)
                                                      !    Dimensions:
                                                      !    (ncol,nlay)


! Local variables
  integer  :: i, j, l, ib, jloop, ijloop
  real*8 :: gcair 
  real*8 :: rhoa, rhob, rhosum
  real*8 :: hr_temp
  real (kind=rb)   :: cldfr_sw(ncol,nlay)

! Pressure related variables
  real*8   :: plev_temp(ncol,nlay)

!McICA variables
  integer(kind=im)      :: iplon
  integer(kind=im),save :: permuteseed
  integer(kind=im)      :: irng=1  ! Mersenne Twister random number generator

  real(kind=rb)         :: relqmcl(ncol,nlay)
  real(kind=rb)         :: reicmcl(ncol,nlay)

! McICA LW specific
  real(kind=rb)         :: cldfmcl(ngptlw,ncol,nlay)
  real(kind=rb)         :: ciwpmcl(ngptlw,ncol,nlay)
  real(kind=rb)         :: clwpmcl(ngptlw,ncol,nlay)
  real(kind=rb)         :: taucmcl_lw(ngptlw,ncol,nlay)

! McICA SW specific
  real(kind=rb)         :: cldfmcl_sw(ngptsw,ncol,nlay)
  real(kind=rb)         :: ciwpmcl_sw(ngptsw,ncol,nlay)
  real(kind=rb)         :: clwpmcl_sw(ngptsw,ncol,nlay)
  real(kind=rb)         :: taucmcl_sw(ngptsw,ncol,nlay)
  real(kind=rb)         :: ssacmcl(ngptsw,ncol,nlay)
  real(kind=rb)         :: asmcmcl(ngptsw,ncol,nlay)
  real(kind=rb)         :: fsfcmcl(ngptsw,ncol,nlay)

   
!Shape array for reshaping arrays from 3D to 2D (nlon,nlat,nlay) to (ncol,nlay)
  integer :: ishape1(1) 
  integer :: ishape2(2) 

  ishape1(1) = ncol
  ishape2(1) = ncol
  ishape2(2) = nlay

!  print *,'rad_driver'
! Gridding is somewhat difficult
! GC provides variables at the center of grid boxes (including T and P)
! I will assume these are layer values (in RRTMG language) and calculate
! the values at the edges and call them level values (sort of backwards
! from what we usually do)



! Reshape arrays from GC to reduce to 2D
       h2ovmr = reshape(avgw,ishape2)
       tlay = reshape(t,ishape2)
       play = reshape(pcenter,ishape2)
       plev_temp = reshape(pedge,ishape2)
       tsfc = reshape(tskin,ishape1)

! Get level values
       gcair = 1.0e-3*gascon/avogad
       do i=1,ncol
          plev(i,1) = plev_temp(i,1)     ! set lowest level to surface pressure
          tlev(i,1) = tlay(i,1)     ! set lowest level to layer temperature  (KLUDGE)
          plev(i,nlay+1) = play(i,nlay)
          tlev(i,nlay+1) = tlay(i,nlay)
          do l=2,nlay
             rhoa = play(i,l-1)/(gcair*tlay(i,l-1))
             rhob = play(i,l)/(gcair*tlay(i,l))
             rhosum = rhoa+rhob
             plev(i,l) = (rhoa*play(i,l-1)+rhob*play(i,l))/rhosum
             tlev(i,l) = (rhoa*tlay(i,l-1)+rhob*tlay(i,l))/rhosum
          end do
       end do

! Fill co2, n2o and o2 arrays with reasonable atmospheric values   

       co2vmr(:,:) = 3.90e-4
!       n2ovmr(:,:) = 3.20e-7
       o2vmr(:,:) =  0.209

       ssacld(:,:,:) = 0.0
       asmcld(:,:,:) = 0.0
       fsfcld(:,:,:) = 0.0


       select case (icld_gc)
!  Cloud setup for clear
       case (0)
          idrv = 0
          icld = 0
          inflglw = 0
          inflgsw = 0
          taucmcl_lw(:,:,:) = 0.0
          taucmcl_sw(:,:,:) = 0.0
          iceflglw = 0       
          liqflglw = 0       
          iceflgsw = 0       
          liqflgsw = 0       
          !print *,'clear'

!  Cloud setup for McICA cloud (only option now)
       case (1)
          idrv = 0
          icld = 2                  !maximum random overlap
          inflglw = 2
          inflgsw = 2
          taucld_lw(:,:,:) = 0.0    ! taucld not used 
          taucld_sw(:,:,:) = 0.0
          taucmcl_lw(:,:,:) = 0.0      ! used only as a check
          taucmcl_sw(:,:,:) = 0.0 
          iceflglw = 2       !Streamer
          liqflglw = 1       !Hu and Stamnes
          iceflgsw = 2       !Streamer
          liqflgsw = 1       !Hu and Stamnes
!          print *,'McICA'
      end select


! Set solar variables
       coszen = suncos
       dyofyr = doy



! Initialize fluxes to avoid nasty surprises
        uflx(:,:) = 0.0 
        dflx(:,:) = 0.0 
        hr(:,:) = 0.0 
        uflxc(:,:) = 0.0 
        dflxc(:,:) = 0.0 
        hrc(:,:) = 0.0 
        duflx_dt(:,:) = 0.0 
        duflxc_dt(:,:) = 0.0 

        swuflx(:,:) = 0.0 
        swdflx(:,:) = 0.0 
        swhr(:,:) = 0.0 
        swuflxc(:,:) = 0.0 
        swdflxc(:,:) = 0.0 
        swhrc(:,:) = 0.0 

        IF (LLWRAD) THEN
!        print *,'will call rrtmg_lw'
        permuteseed=iseed+ngptsw+1


        call mcica_subcol_lw ( iplon, ncol, nlay, icld, permuteseed, irng,play,&
                       cldfr, cicewp, cliqwp, reice, reliq, taucld_lw, cldfmcl,&
                       ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl_lw)
        call rrtmg_lw &
            (ncol    ,nlay    ,icld    ,idrv    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr, &
             cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
             inflglw ,iceflglw,liqflglw,cldfmcl   , &
             taucmcl_lw ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl   , &
             tauaer_lw  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, &
             duflx_dt,duflxc_dt )
    ENDIF
        
        IF (LSWRAD) THEN

!        print *,'will call rrtmg_sw'
        permuteseed=permuteseed+ngptlw+1
        call mcica_subcol_sw(iplon, ncol, nlay, icld, permuteseed, irng, play, &
                       cldfr, cicewp, cliqwp, reice, reliq, taucld_sw, ssacld, asmcld, fsfcld, &
                       cldfmcl_sw, ciwpmcl_sw, clwpmcl_sw, reicmcl, relqmcl, &
                       taucmcl_sw, ssacmcl, asmcmcl, fsfcmcl)
        
        call rrtmg_sw &
            (ncol    ,nlay    ,icld    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr, &
             asdir   ,asdif   ,aldir   ,aldif   , &  
             coszen  ,adjes   ,dyofyr  ,scon    , &  
             inflgsw ,iceflgsw,liqflgsw,cldfmcl_sw , &  
             taucmcl_sw  ,ssacmcl  ,asmcmcl  ,fsfcmcl  , &  
             ciwpmcl_sw  ,clwpmcl_sw  ,reicmcl   ,relqmcl   , &  
             tauaer_sw  ,ssaaer  ,asmaer  ,ecaer   , &  
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc)
        ENDIF
!

             lw_uflux = uflx
             lw_dflux = dflx
             sw_uflux = swuflx
             sw_dflux = swdflx
             lw_ufluxc = uflxc
             lw_dfluxc = dflxc
             sw_ufluxc = swuflxc
             sw_dfluxc = swdflxc

!         
!       

100          format (i10,3f10.2,8f15.3)


        return

        end subroutine rad_driver
