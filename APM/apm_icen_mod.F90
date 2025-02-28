#ifdef APM
!program test
!
!use APM_ICEN_MOD, only : nucleati
!
!!
!  parameter (naer_all = 3)
!
!  real*8   wbar                ! grid cell mean vertical velocity (m/s)
!  real*8 :: tair                ! temperature (K)
!  real*8 :: relhum              ! relative humidity with respective to liquid
!  real*8 :: cldn                ! new value of cloud fraction    (fraction)
!  real*8 :: qc                  ! liquid water mixing ratio (kg/kg)
!!  real*8 :: nfice               ! ice mass fraction
!  real*8 :: na(naer_all)        ! aerosol number concentration (/m3)
!!
!!
!  real*8 :: nuci               ! ice number nucleated (#/kg)
!  real*8 :: onihf              ! nucleated number from homogeneous freezing of so4
!  real*8 :: oniimm             ! nucleated number from immersion freezing
!  real*8 :: onidep             ! nucleated number from deposition nucleation
!  real*8 :: onimey             ! nucleated number from deposition nucleation  (meyers: mixed phase)
!
!
!  wbar = 0.25 
!  relhum = 0.8
!  cldn = 0.3
!  qc   = 1.d-4
!!  nfice = 0.1   ! not used
!
!! already in #/cm3
!  na(1) = 10.0 !SP
!  na(2) = 1.0d-2  !dust
!  na(3) = 10.0  !soot
!
!  DO I=1,50
!     tair = 220.0 + float(I) 
!     call nucleati(wbar, tair, relhum, cldn, qc,     &
!       na, nuci, onihf, oniimm, onidep, onimey)
!
!     write(6,100)tair,nuci, onihf, oniimm, onidep, onimey
!100  FORMAT(100(1PE10.3))
!  ENDDO
!
!end program
!



MODULE APM_ICEN_MOD

!---------------------------------------------------------------------------------
! Purpose:
!   CAM Interface for aerosol activation


 implicit none
 private

 public :: nucleati 

contains

subroutine nucleati(wbar, tair, relhum, cldn, qc,    &
       na, nuci, onihf, oniimm, onidep, onimey)
 
!---------------------------------------------------------------
! Purpose:
!  The parameterization of ice nucleation.
!
! Method: The current method is based on Liu & Penner (2005)
!  It related the ice nucleation with the aerosol number, temperature and the
!  updraft velocity. It includes homogeneous freezing of sulfate, immersion
!  freezing of soot, and Meyers et al. (1992) deposition nucleation
!
! Authors: Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
!----------------------------------------------------------------

!-----------------------------------------------------
! Input Arguments
!
  real*8, intent(in) :: wbar                ! grid cell mean vertical velocity (m/s)
  real*8, intent(in) :: tair                ! temperature (K)
  real*8, intent(in) :: relhum              ! relative humidity with respective to liquid
  real*8, intent(in) :: cldn                ! new value of cloud fraction    (fraction)
  real*8, intent(in) :: qc                  ! liquid water mixing ratio (kg/kg)
!  real*8, intent(in) :: nfice               ! ice mass fraction
!  real*8, intent(in) :: na(naer_all)        ! aerosol number concentration (/m3) Yu --> in #/cm3
  real*8, intent(in) :: na(3)        ! aerosol number concentration (/m3) Yu --> in #/cm3


!
! Output Arguments
!
  real*8, intent(out) :: nuci               ! ice number nucleated (#/kg) -->#/cm3 (Yu)
  real*8, intent(out) :: onihf              ! nucleated number from homogeneous freezing of so4
  real*8, intent(out) :: oniimm             ! nucleated number from immersion freezing
  real*8, intent(out) :: onidep             ! nucleated number from deposition nucleation
  real*8, intent(out) :: onimey             ! nucleated number from deposition nucleation  (meyers: mixed phase)

!
! Local workspace
!
  real*8  so4_num                                      ! so4 aerosol number (#/cm^3)
  real*8  soot_num                                     ! soot (hydrophilic) aerosol number (#/cm^3)
!  real*8  dst1_num,dst2_num,dst3_num,dst4_num          ! dust aerosol number (#/cm^3)
  real*8  dst_num                                      ! total dust aerosol number (#/cm^3)
  real*8  nihf                                         ! nucleated number from homogeneous freezing of so4
  real*8  niimm                                        ! nucleated number from immersion freezing
  real*8  nidep                                        ! nucleated number from deposition nucleation
  real*8  nimey                                        ! nucleated number from deposition nucleation (meyers)
  real*8  n1,ni                                        ! nucleated number
  real*8  tc,A,B,C,regm,RHw                            ! work variable
  real*8  esl,esi,deles                                ! work variable
  real*8  dst_scale
  real*8  subgrid
  real*8 dmc,ssmc         ! variables for modal scheme.

    so4_num=0.0
    soot_num=0.0
    dst_num=0.0
!    dst1_num = 0.0
!    dst2_num = 0.0
!    dst3_num = 0.0
!    dst4_num = 0.0     

!For modal aerosols, assume for the upper troposphere:
! soot = accumulation mode
! sulfate = aiken mode
! dust = coarse mode
! since modal has internal mixtures.

!    if(idxsul .gt. 0) then 
!       so4_num=na(idxsul)*1.0e-6 ! #/cm^3
!    end if
!
!continue above philosophy here....

!    if(idxbcphi .gt. 0) then 
!      soot_num=na(idxbcphi)*1.0e-6 !#/cm^3
!    end if
!
!    if(idxdst1 .gt. 0) then 
!       dst1_num=na(idxdst1)  *1.0e-6 !#/cm^3
!    end if
!
!    if(idxdst2 .gt. 0) then 
!       dst2_num=na(idxdst2)*1.0e-6 !#/cm^3
!    end if
!
!    if(idxdst3 .gt. 0) then 
!       dst3_num=na(idxdst3)*1.0e-6 !#/cm^3
!    end if
!
!    if(idxdst4 .gt. 0) then 
!       dst4_num=na(idxdst4)*1.0e-6 !#/cm^3
!    end if
!
!    dst_num =dst1_num+dst2_num+dst3_num+dst4_num

    so4_num = na(1)
    dst_num = na(2)
    soot_num = na(3)

! no soot nucleation 
!    soot_num=0.0


    ni=0.
    tc=tair-273.15

    ! initialize
    niimm=0.
    nidep=0.
    nihf=0.

    if(so4_num.ge.1.0e-10 .and. (soot_num+dst_num).ge.1.0e-10 .and. cldn.gt.0.) then

!-----------------------------
! RHw parameterization for heterogeneous immersion nucleation
    A = 0.0073
    B = 1.477
    C = 131.74
    RHw=(A*tc*tc+B*tc+C)*0.01   ! RHi ~ 120-130%

    subgrid = 1.2

    if((tc.le.-35.0) .and. ((relhum*polysvp(tair,0)/polysvp(tair,1)*subgrid).ge.1.2)) then ! use higher RHi threshold

       A = -1.4938 * log(soot_num+dst_num) + 12.884
       B = -10.41  * log(soot_num+dst_num) - 67.69
       regm = A * log(wbar) + B

!       WRITE(6,*)regm

       if(tc.gt.regm) then    ! heterogeneous nucleation only
         if(tc.lt.-40. .and. wbar.gt.1.) then ! exclude T<-40 & W>1m/s from hetero. nucleation
           call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
           niimm=0.
           nidep=0.
           n1=nihf
         else
           call hetero(tc,wbar,soot_num+dst_num,niimm,nidep)
           nihf=0.
           n1=niimm+nidep
         endif
       elseif (tc.lt.regm-5.) then ! homogeneous nucleation only
         call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
         niimm=0.
         nidep=0.
         n1=nihf
       else        ! transition between homogeneous and heterogeneous: interpolate in-between
         if(tc.lt.-40. .and. wbar.gt.1.) then ! exclude T<-40 & W>1m/s from hetero. nucleation
           call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
           niimm=0.
           nidep=0.
           n1=nihf
         else

           call hf(regm-5.,wbar,relhum,subgrid,so4_num,nihf)
           call hetero(regm,wbar,soot_num+dst_num,niimm,nidep)

           if(nihf.le.(niimm+nidep)) then
             n1=nihf
           else
             n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5.)
           endif
         endif
       endif

       ni=n1

    endif
    endif
1100  continue

! deposition/condensation nucleation in mixed clouds (-40<T<0C) (Meyers, 1992)
    if(tc.lt.0. .and. tc.gt.-37. .and. qc.gt.1.e-12) then
      esl = polysvp(tair,0)     ! over water in mixed clouds
      esi = polysvp(tair,1)     ! over ice
      deles = (esl - esi)
      nimey=1.e-3*exp(12.96*deles/esi - 0.639) 
    else
      nimey=0.
    endif

    nuci=ni+nimey
    if(nuci.gt.9999..or.nuci.lt.0.) then
       write(6, *) 'Warning: incorrect ice nucleation number (nuci reset =0)'
       write(6, *) ni, tair, relhum, wbar, nihf, niimm, nidep,deles,esi,dst_num
       nuci=0.
    endif

!    nuci=nuci*1.e+6/rhoair    ! change unit from #/cm3 to #/kg
!    onimey=nimey*1.e+6/rhoair
!    onidep=nidep*1.e+6/rhoair
!    oniimm=niimm*1.e+6/rhoair
!    onihf=nihf*1.e+6/rhoair

!Yu +

    nuci=nuci
    onihf=nihf
    onidep=nidep
    oniimm=niimm
    onimey=nimey

  return
  end subroutine nucleati

  subroutine hetero(T,ww,Ns,Nis,Nid)

    real*8, intent(in)  :: T, ww, Ns
    real*8, intent(out) :: Nis, Nid

    real*8 A11,A12,A21,A22,B11,B12,B21,B22
    real*8 A,B,C

!---------------------------------------------------------------------
! parameters

      A11 = 0.0263
      A12 = -0.0185
      A21 = 2.758
      A22 = 1.3221
      B11 = -0.008
      B12 = -0.0468
      B21 = -0.2667
      B22 = -1.4588

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0    ! don't include deposition nucleation for cirrus clouds when T<-37C

      return
  end subroutine hetero


  subroutine hf(T,ww,RH,subgrid,Na,Ni)

      real*8, intent(in)  :: T, ww, RH, subgrid, Na
      real*8, intent(out) :: Ni

      real*8    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
      real*8    A2_fast,B2_fast
      real*8    C1_fast,C2_fast,k1_fast,k2_fast
      real*8    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
      real*8    C1_slow,C2_slow,k1_slow,k2_slow
      real*8    regm
      real*8    A,B,C
      real*8    RHw

!---------------------------------------------------------------------
! parameters

      A1_fast  =0.0231
      A21_fast =-1.6387  !(T>-64 deg)
      A22_fast =-6.045   !(T<=-64 deg)
      B1_fast  =-0.008
      B21_fast =-0.042   !(T>-64 deg)
      B22_fast =-0.112   !(T<=-64 deg)
      C1_fast  =0.0739
      C2_fast  =1.2372

      A1_slow  =-0.3949
      A2_slow  =1.282
      B1_slow  =-0.0156
      B2_slow  =0.0111
      B3_slow  =0.0217
      C1_slow  =0.120
      C2_slow  =2.312

      Ni = 0.0

!----------------------------
!RHw parameters
      A = 6.0e-4*log(ww)+6.6e-3
      B = 6.0e-2*log(ww)+1.052
      C = 1.68  *log(ww)+129.35
      RHw=(A*T*T+B*T+C)*0.01

!      WRITE(6,*)"RHw=", RHw

      if((T.le.-37.0) .and. ((RH*subgrid).ge.RHw)) then

        regm = 6.07*log(ww)-55.0

        if(T.ge.regm) then    ! fast-growth regime

          if(T.gt.-64.0) then
            A2_fast=A21_fast
            B2_fast=B21_fast
          else
            A2_fast=A22_fast
            B2_fast=B22_fast
          endif

          k1_fast = exp(A2_fast + B2_fast*T + C2_fast*log(ww))
          k2_fast = A1_fast+B1_fast*T+C1_fast*log(ww)

          Ni = k1_fast*Na**(k2_fast)
          Ni = min(Ni,Na)

        else       ! slow-growth regime

          k1_slow = exp(A2_slow + (B2_slow+B3_slow*log(ww))*T + C2_slow*log(ww))
          k2_slow = A1_slow+B1_slow*T+C1_slow*log(ww)

          Ni = k1_slow*Na**(k2_slow)
          Ni = min(Ni,Na)

        endif

      end if

      return
  end subroutine hf

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! error function in single precision
!
!    Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
!    You may use, copy, modify this code for any purpose and 
!    without fee. You may distribute this ORIGINAL package.

      function derf(x)
      implicit real (a - h, o - z)
      real*8 a,b,x
      dimension a(0 : 64), b(0 : 64)
      integer i,k
      data (a(i), i = 0, 12) / & 
         0.00000000005958930743d0, -0.00000000113739022964d0, & 
         0.00000001466005199839d0, -0.00000016350354461960d0, &
         0.00000164610044809620d0, -0.00001492559551950604d0, &
         0.00012055331122299265d0, -0.00085483269811296660d0, &
         0.00522397762482322257d0, -0.02686617064507733420d0, &
         0.11283791670954881569d0, -0.37612638903183748117d0, &
         1.12837916709551257377d0 / 
      data (a(i), i = 13, 25) / &
         0.00000000002372510631d0, -0.00000000045493253732d0, &
         0.00000000590362766598d0, -0.00000006642090827576d0, &
         0.00000067595634268133d0, -0.00000621188515924000d0, &
         0.00005103883009709690d0, -0.00037015410692956173d0, &
         0.00233307631218880978d0, -0.01254988477182192210d0, &
         0.05657061146827041994d0, -0.21379664776456006580d0, &
         0.84270079294971486929d0 / 
      data (a(i), i = 26, 38) / &
         0.00000000000949905026d0, -0.00000000018310229805d0, &
         0.00000000239463074000d0, -0.00000002721444369609d0, &
         0.00000028045522331686d0, -0.00000261830022482897d0, &
         0.00002195455056768781d0, -0.00016358986921372656d0, &
         0.00107052153564110318d0, -0.00608284718113590151d0, &
         0.02986978465246258244d0, -0.13055593046562267625d0, &
         0.67493323603965504676d0 / 
      data (a(i), i = 39, 51) / &
         0.00000000000382722073d0, -0.00000000007421598602d0, &
         0.00000000097930574080d0, -0.00000001126008898854d0, &
         0.00000011775134830784d0, -0.00000111992758382650d0, &
         0.00000962023443095201d0, -0.00007404402135070773d0, &
         0.00050689993654144881d0, -0.00307553051439272889d0, &
         0.01668977892553165586d0, -0.08548534594781312114d0, &
         0.56909076642393639985d0 / 
      data (a(i), i = 52, 64) / &
         0.00000000000155296588d0, -0.00000000003032205868d0, &
         0.00000000040424830707d0, -0.00000000471135111493d0, &
         0.00000005011915876293d0, -0.00000048722516178974d0, &
         0.00000430683284629395d0, -0.00003445026145385764d0, &
         0.00024879276133931664d0, -0.00162940941748079288d0, &
         0.00988786373932350462d0, -0.05962426839442303805d0, &
         0.49766113250947636708d0 / 
      data (b(i), i = 0, 12) / &
         -0.00000000029734388465d0, 0.00000000269776334046d0, &
         -0.00000000640788827665d0, -0.00000001667820132100d0, &
         -0.00000021854388148686d0, 0.00000266246030457984d0, &
         0.00001612722157047886d0, -0.00025616361025506629d0, &
         0.00015380842432375365d0, 0.00815533022524927908d0, &
         -0.01402283663896319337d0, -0.19746892495383021487d0,& 
         0.71511720328842845913d0 / 
      data (b(i), i = 13, 25) / &
         -0.00000000001951073787d0, -0.00000000032302692214d0, &
         0.00000000522461866919d0, 0.00000000342940918551d0, &
         -0.00000035772874310272d0, 0.00000019999935792654d0, &
         0.00002687044575042908d0, -0.00011843240273775776d0, &
         -0.00080991728956032271d0, 0.00661062970502241174d0, &
         0.00909530922354827295d0, -0.20160072778491013140d0, &
         0.51169696718727644908d0 / 
      data (b(i), i = 26, 38) / &
         0.00000000003147682272d0, -0.00000000048465972408d0, &
         0.00000000063675740242d0, 0.00000003377623323271d0, &
         -0.00000015451139637086d0, -0.00000203340624738438d0,& 
         0.00001947204525295057d0, 0.00002854147231653228d0, &
         -0.00101565063152200272d0, 0.00271187003520095655d0, &
         0.02328095035422810727d0, -0.16725021123116877197d0, &
         0.32490054966649436974d0 / 
      data (b(i), i = 39, 51) / &
         0.00000000002319363370d0, -0.00000000006303206648d0, &
         -0.00000000264888267434d0, 0.00000002050708040581d0, &
         0.00000011371857327578d0, -0.00000211211337219663d0, &
         0.00000368797328322935d0, 0.00009823686253424796d0, &
         -0.00065860243990455368d0, -0.00075285814895230877d0,& 
         0.02585434424202960464d0, -0.11637092784486193258d0, &
         0.18267336775296612024d0 / 
      data (b(i), i = 52, 64) / &
         -0.00000000000367789363d0, 0.00000000020876046746d0, &
         -0.00000000193319027226d0, -0.00000000435953392472d0, &
         0.00000018006992266137d0, -0.00000078441223763969d0, &
         -0.00000675407647949153d0, 0.00008428418334440096d0, &
         -0.00017604388937031815d0, -0.00239729611435071610d0, &
         0.02064129023876022970d0, -0.06905562880005864105d0, &
         0.09084526782065478489d0 / 
      w = abs(x)
      if (w .lt. 2.2d0) then
          t = w * w
          k = int(t)
          t = t - k
          k = k * 13
          y = ((((((((((((a(k) * t + a(k + 1)) * t + &
             a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t + &
             a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t + &
             a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t + &
             a(k + 11)) * t + a(k + 12)) * w
      else if (w .lt. 6.9d0) then
          k = int(w)
          t = w - k
          k = 13 * (k - 2)
          y = (((((((((((b(k) * t + b(k + 1)) * t + &
             b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
             b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
             b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
             b(k + 11)) * t + b(k + 12)
          y = y * y
          y = y * y
          y = y * y
          y = 1 - y * y
      else
          y = 1
      end if
      if (x .lt. 0) y = -y
      derf = y
      end function derf
!



      function polysvp (T,type)
!  Compute saturation vapor pressure by using
! function from Goff and Gatch (1946)

!  Polysvp returned in units of pa.
!  T is input in units of K.
!  type refers to saturation with respect to liquid (0) or ice (1)

      real*8 dum

      real*8 T,polysvp

      integer type

! ice

      if (type.eq.1) then

! Goff Gatch equation (good down to -100 C)

         polysvp = 10.**(-9.09718*(273.16/t-1.)-3.56654* &
          log10(273.16/t)+0.876793*(1.-t/273.16)+ &
          log10(6.1071))*100.

      end if

! Goff Gatch equation, uncertain below -70 C

      if (type.eq.0) then
         polysvp = 10.**(-7.90298*(373.16/t-1.)+ &
             5.02808*log10(373.16/t)- &
             1.3816e-7*(10.**(11.344*(1.-t/373.16))-1.)+ &
             8.1328e-3*(10.**(-3.49149*(373.16/t-1.))-1.)+ &
             log10(1013.246))*100.
         end if


      end function polysvp
!--xl

end module APM_ICEN_MOD
#endif
