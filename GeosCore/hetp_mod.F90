!######### HETP: AEROSOL THERMODYNAMIC EQUILIBRIUM OF THE ###################
!#########       NH4-SO4-NO3-Na-Cl-Ca-Mg-K SYSTEM         ###################
!
!Copyright (C) 2023  Stefan Miller, Environment and Climate Change Canada
!                    Contact: Stefan.Miller (at) ec.gc.ca
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!############################################################################
! ## HETP Code
! ## Module file for HETP code
! ## Contains: (1) parameters related to ITP/activity coefficient calculations
! ##           (2) ZSR data arrays 
! ##           (3) other parameters 
!
! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
! ## Written by Stefan Miller
!############################################################################
module hetp_mod

    IMPLICIT NONE
    PRIVATE

    PUBLIC  :: mach_hetp_main_15cases
    PRIVATE :: mach_hetp_calca2
    PRIVATE :: mach_hetp_calcb4
    PRIVATE :: mach_hetp_calcc2
    PRIVATE :: mach_hetp_calcd3
    PRIVATE :: mach_hetp_calce4
    PRIVATE :: mach_hetp_calcf2
    PRIVATE :: mach_hetp_calcg5
    PRIVATE :: mach_hetp_calch6
    PRIVATE :: mach_hetp_calci6
    PRIVATE :: mach_hetp_calcj3
    PRIVATE :: mach_hetp_calco7
    PRIVATE :: mach_hetp_calcm8
    PRIVATE :: mach_hetp_calcp13
    PRIVATE :: mach_hetp_calcl9
    PRIVATE :: mach_hetp_calck4
    PRIVATE :: mach_hetp_calcact1
    PRIVATE :: mach_hetp_calcact1b
    PRIVATE :: mach_hetp_calcact2
    PRIVATE :: mach_hetp_calcact2b
    PRIVATE :: mach_hetp_calcact3
    PRIVATE :: mach_hetp_calcact3b
    PRIVATE :: mach_hetp_calcact4
    PRIVATE :: mach_hetp_calcact4b
    PRIVATE :: mach_hetp_poly3v 
    PRIVATE :: mach_hetp_poly
    PRIVATE :: mach_hetp_adjust

    CONTAINS

    !###################################################################################
    ! ## HETP Code
    ! ## Branch/subcase determination for metastable state 
    ! AEROSOL THERMODYNAMIC EQUILIBRIUM OF THE NH4-SO4-NO3-Na-Cl-Ca-Mg-K system
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !###################################################################################
    subroutine mach_hetp_main_15cases(TS, TA, TN, TNa, TCl, TCa, TK, TMg, temp, rh,    &
                                      so4, hso4, caso4, nh4, nh3, no3, hno3, cl, hcl,  &
                                      na, ca, k, mg, h, oh, lwc, frna, frca, frk, frmg,&
                                      frso4, case_number)
    !
       use mach_hetp_mod     
       implicit none
    !
    !  INPUT VARIABLES (total gas + aerosol, input as mol/m3 air)
       real(kind=8),    intent(in)  :: TS    ! Total sulfate
       real(kind=8),    intent(in)  :: TA    ! Total ammonium
       real(kind=8),    intent(in)  :: TN    ! Total nitrate
       real(kind=8),    intent(in)  :: TNa   ! Total sodium
       real(kind=8),    intent(in)  :: TCl   ! Total chloride
       real(kind=8),    intent(in)  :: TCa   ! Total calcium
       real(kind=8),    intent(in)  :: TK    ! Total potassium
       real(kind=8),    intent(in)  :: TMg   ! Total magnesium
       real(kind=8),    intent(in)  :: temp  ! Air temperature (K)
       real(kind=8),    intent(in)  :: rh    ! Relative humidity (0-1 scale; zero allowed)
    !
    !  OUTPUT VARIABLES (output as mol/m3 air)
       real(kind=8),    intent(out) :: so4     ! SO4-- (aq)
       real(kind=8),    intent(out) :: hso4    ! HSO4- (aq)
       real(kind=8),    intent(out) :: caso4   ! CaSO4 (s)
       real(kind=8),    intent(out) :: nh4     ! NH4+  (aq)
       real(kind=8),    intent(out) :: nh3     ! NH3   (g)
       real(kind=8),    intent(out) :: no3     ! NO3-  (aq)
       real(kind=8),    intent(out) :: hno3    ! HNO3  (g)
       real(kind=8),    intent(out) :: cl      ! Cl-   (aq)
       real(kind=8),    intent(out) :: hcl     ! HCl   (g)
       real(kind=8),    intent(out) :: na      ! Na+   (aq)
       real(kind=8),    intent(out) :: ca      ! Ca2+  (aq)  
       real(kind=8),    intent(out) :: k       ! K+    (aq)
       real(kind=8),    intent(out) :: mg      ! Mg2+  (aq)
       real(kind=8),    intent(out) :: h       ! H+    (aq)
       real(kind=8),    intent(out) :: oh      ! OH-   (aq)
       real(kind=8),    intent(out) :: lwc     ! Aerosol liquid water content
       real(kind=8),    intent(out) :: frna    ! Free Na  
       real(kind=8),    intent(out) :: frca    ! Free Ca   
       real(kind=8),    intent(out) :: frk     ! Free K  
       real(kind=8),    intent(out) :: frmg    ! Free Mg
       real(kind=8),    intent(out) :: frso4   ! Free SO4
       real(kind=4),    intent(out) :: case_number 
    !!if_off
    !
    !  ## Local variables ##
       integer(kind=4), parameter  :: nr   = 15       !number of eq reactions
    !
       real(kind=8), dimension(nr) :: k0 = (/3569119.63798333d0, 0.983992783234808d0, 57.639d0,               &
                                             1.805d-5, 1.015d-2, 141.166162470542d0,                          &
                                             1.71700241471167d0, 30.9720858820060d0, 5.746d-17,               &
                                             1.010d-14, 1.971d6, 37.6654800230020d0,                          &
                                             12.2216391328847d0, 0.474106585740804d0, 2.511d6/)
    !
       real(kind=8), dimension(nr) :: p1 = (/29.4737118977339d0, 1.613671606774050d-2, 13.79d0,               &
                                             -1.50d0, 8.85d0, -2.88443799710923d0,                            &
                                             -2.63431889805915d0, -5.19198839479649d0,-74.38d0,               &
                                             -22.52d0, 30.20d0, -1.56284095116099d0,                          &
                                             -8.20148594143076d0, 0.992408038166253d0, 29.17d0/)
    !
       real(kind=8), dimension(nr) :: p2 = (/16.8318498917489d0, 0.0d0, -5.393d0, 26.920d0,                   &
                                             25.140d0, 15.8287226365167d0, 38.5722877074814d0,                &
                                             54.4022131344720d0, 6.120d0, 26.920d0,                           &
                                             19.910d0, 16.8992061582872d0, 16.0067356266538d0,                &
                                             39.4996391628578d0, 16.830d0/)
    !
    !
       real(kind=8) :: sulrat, sodrat, so4rat, crnarat, crrat, rest
       real(kind=8) :: ccaso4i, cafri, ccano32i, cnacli, cmgno32i, ccacl2i, rest1
       real(kind=8) :: cna2so4i, frso4i, nafri, cnano3i, no3fr, rest2, cmgso4i
       real(kind=8) :: frmgi, no3fri, cmgcl2i, clfri, rest3
       real(kind=8) :: frnh4, frno3, frcl, frca_init, frna_init, frk_init, frmg_init
    !
    !
    !  On input, set output aqueous phase species to the total gas + aerosol concentration; 
    !  the output speciation will be partitioned and these variables reused 
       so4 = TS
       nh4 = TA
       no3 = TN
       na  = TNa
       cl  = TCl
       ca  = TCa
       k   = TK
       mg  = TMg
    !
    !  Set other output species to 0.0d0
       hso4  = 0.0d0
       caso4 = 0.0d0
       nh3   = 0.0d0
       hno3  = 0.0d0
       hcl   = 0.0d0
       h     = 0.0d0
       oh    = 0.0d0
       lwc   = 0.0d0
       frna  = 0.0d0
       frca  = 0.0d0
       frk   = 0.0d0
       frmg  = 0.0d0
    !
       frso4     = 0.0d0
       frnh4     = 0.0d0
       frno3     = 0.0d0
       frcl      = 0.0d0
       frca_init = 0.0d0
       frna_init = 0.0d0
       frk_init  = 0.0d0
       frmg_init = 0.0d0
    !
    !  Sulfate ratios; initialize to 0.0d0
       sulrat  = 0.0d0
       sodrat  = 0.0d0
       so4rat  = 0.0d0
       crnarat = 0.0d0
       crrat   = 0.0d0
    !
    !  #######################
    !  ### SORTING SECTION ###
    !  #######################
    !
    !  ### Branch 0 ###
    !  ### No cases, all species < tiny: do nothing 
       if (TS + TA + TN + TNa + TCl + TCa + TK + TMg <= tiny) then
          case_number = 0.0d0
    !
    !  ### Branch 1 ###
    !  ISRP1F: Only sulfate and ammonium 
       else if (TN + TNa + TCl + TCa + TK + TMg <= tiny) then
          so4 = max(so4, tiny)
          nh4 = max(nh4, tiny)
          no3 = 0.0d0
          na  = 0.0d0
          cl  = 0.0d0
          ca  = 0.0d0
          k   = 0.0d0
          mg  = 0.0d0
    !
    !  ## Calculate sulfate ratio
          sulrat = nh4 / so4
    !
    !  ## 1. Sulfate poor (case: calca2)
          if (2.0 <= sulrat) then
             case_number = 1.0d0
             call mach_hetp_calca2(so4, nh4, nh3, hso4, h, oh,          &
                                   lwc, rh, temp, k0, p1, p2, nr)
    !
    !  ## 2. Sulfate rich, no acid (case: calcb4)
          else if (1.0 <= sulrat .and. sulrat < 2.0) then
             case_number = 2.0d0
             call mach_hetp_calcb4(so4, nh4, nh3, hso4, h, lwc,         &
                                   rh, temp, k0, p1, p2, nr)
    !
    !  ## 3. Sulfate rich, free acid (case: calcc2)
          else if (sulrat < 1.0) then
             case_number = 3.0d0
             call mach_hetp_calcc2(so4, nh4, nh3, hso4, h, lwc,         &
                                   rh, temp, k0, p1, p2, nr)           
          end if
    !
    !  ### Branch 2 ###
    !  ISRP2F: Only sulfate, ammonium and nitrate 
       else if (TNa + TCl + TCa + TK + TMg <= tiny) then
          so4 = max(so4, tiny)
          nh4 = max(nh4, tiny)
          no3 = max(no3, tiny)
          na  = 0.0d0
          cl  = 0.0d0
          ca  = 0.0d0
          k   = 0.0d0
          mg  = 0.0d0
    !
    !     ## Calculate sulfate ratio
          sulrat  = nh4 / so4
    !
    !  ## 4. Sulfate poor (case: calcd3)
          if (2.0 <= sulrat) then
             case_number = 4.0d0
             call mach_hetp_calcd3(so4, nh4, hno3, nh3, hso4, h, no3, &
                                   lwc, rh, temp, k0, p1, p2, nr)
    !
    !  ## 5. Sulfate rich, no acid (case: calce4)
          else if (1.0 <= sulrat .and. sulrat < 2.0) then
             case_number = 5.0d0
             call mach_hetp_calce4(so4, nh4, hno3, hso4, h, no3,       &
                                   lwc, rh, temp, k0, p1, p2, nr)
    !
    !  ## 6. Sulfate rich, free acid (case: calcf2)
          elseif (sulrat < 1.0) then
             case_number = 6.0d0
             call mach_hetp_calcf2(so4, nh4, hno3, hso4, h, no3,       &
                                   lwc, rh, temp, k0, p1, p2, nr)
          end if
    !
    !  ### Branch 3 ###
    !  ISRP3F: Only sulfate, ammonium, nitrate, chloride and calcium 
       else if (TCa + TK + TMg <= tiny) then
    !  ## Adjust for too little ammonium and chloride            
          nh4 = max(nh4, tiny)     ! In ISORROPIA II nh4 and cl are set to 1.0d-10 (not tiny)
          cl  = max(cl,  tiny)     ! creating 'excess mass', HETP uses tiny for mass conservation
    !
    !  ## Adjust for too little sodium, sulfate and nitrate combined
          if (na + so4 + no3 <= tiny) then
             na  = 1.0d-20  ! In ISORROPIA II na and so4 are set to 1.0d-10 (not tiny)
             so4 = 1.0d-20   ! creating 'excess mass', HETP uses tiny for mass conservation
          end if
    !
          na  = max(na,  tiny)
          so4 = max(so4, tiny)
          no3 = max(no3, tiny)
          ca  = 0.0d0
          k   = 0.0d0
          mg  = 0.0d0
    !
    !  ## Check if too much sodium, if too much sodium adjust
    !  Keep track of free amount that results after adjustment
          if (na > 2.0d0*so4 + no3 + cl) then
             frna_init = na - (1.0 - 1.0d-6)*(2.0d0*so4 + no3 + cl)
             na        = (1.0 - 1.0d-6)*(2.0d0*so4 + no3 + cl)
    !        write(*,*), 'Warning error: Na adjusted'
          end if
    !
    !  ## Calculate sulfate ratios
          sulrat = (na + nh4) / so4
          sodrat = na / so4
    !
    !  ## 7. Sulfate poor and sodium poor (case: calcg5)
          if (2.0 <= sulrat .and. sodrat < 2.0) then
             case_number = 7.0d0
             call mach_hetp_calcg5(so4, nh4, nh3, hno3, hcl, hso4,      &
                                   na, cl, no3, h, lwc, rh, temp, k0,   &
                                   p1, p2, nr)
    !
    !  ## 8. Sulfate poor and sodium rich (case: calch6)
          else if (sulrat >= 2.0 .and. sodrat >= 2.0) then
             case_number = 8.0d0
             call mach_hetp_calch6(so4, nh4, nh3, hno3, hcl, hso4,      &
                                   na, cl, no3, h, lwc, frna, rh,       &
                                   temp, k0, p1, p2, nr)
    !
    !  ## 9. Sulfate rich, no acid (case: calci6)
          else if (1.0 <= sulrat .and. sulrat < 2.0) then
             case_number = 9.0d0
             call mach_hetp_calci6(so4, nh4, nh3, hno3, hcl, hso4,      &
                                   na, cl, no3, h, lwc, frna, rh,       &
                                   temp, k0, p1, p2, nr)
    !
    !  ## 10. Sulfate rich, free acid (case: calcj3)
          else if (sulrat < 1.0) then
             case_number = 10.0d0
             call mach_hetp_calcj3(so4, nh4, nh3, hno3, hcl, hso4,      &
                                   na, cl, no3, h, lwc, rh, temp, k0,   &
                                   p1, p2, nr)
          end if
    !
    !  ### Branch 4 ###
    !  ISRP4F: All species are possibly present, at least one crustal species 
    !  (Ca, K, Mg) is present; Na and Cl are not necessarily present 
       else
    !  ######In HETP; do not perform the below mass adjustments############
    !  ####################################################################
    !   ## Adjust for too little ammonium and chloride
    !      nh4 = max(nh4, 1.0d-10)
    !      cl  = max(cl,  1.0d-10)
    !
    !   ## Adjust for too little sodium, sulfate and nitrate combined 
    !      if (na + so4 + no3 <= 1.0d-10) then
    !         na  = 1.0d-10
    !         so4 = 1.0d-10
    !      end if   
    !  ####################################################################
    !
          na  = max(na,  tiny)
          so4 = max(so4, tiny)
          nh4 = max(nh4, tiny)
          no3 = max(no3, tiny)
          cl  = max(cl,  tiny)
          ca  = max(ca,  tiny)
          k   = max(k,   tiny)
          mg  = max(mg,  tiny)
    !
    !  ## Check if too much sodium + crustals; if too much adjust
          rest = 2.0*so4 + no3 + cl
          if (na + ca + k + mg > rest) then
             ccaso4i  = min(so4, ca)
             frso4i   = max(so4 - ccaso4i, 0.0d0)
             cafri    = max(ca  - ccaso4i, 0.0d0)
             ccano32i = min(cafri, 0.5d0*no3)
             cafri    = max(cafri - ccano32i, 0.0d0)
             no3fri   = max(no3 - 2.0d0*ccano32i, 0.0d0)
             ccacl2i  = min(cafri, 0.5d0*cl)
             clfri    = max(cl - 2.0d0*ccacl2i, 0.0d0)
             rest1    = 2.0d0*frso4i + no3fri + clfri
    !
             cna2so4i = min(frso4i, 0.5d0*na)
             frso4i   = max(frso4i - cna2so4i, 0.0d0)
             nafri    = max(na - 2.0d0*cna2so4i, 0.0d0)
             cnacli   = min(nafri, clfri)
             nafri    = max(nafri - cnacli, 0.0d0)
             clfri    = max(clfri - cnacli, 0.0d0)
             cnano3i  = min(nafri, no3fri)
             no3fr    = max(no3fri - cnano3i, 0.0d0)
             rest2    = 2.0d0*frso4i + no3fri + clfri
    !
             cmgso4i  = min(frso4i, mg)
             frmgi    = max(mg - cmgso4i, 0.0d0)
             frso4i   = max(frso4i - cmgso4i, 0.0d0)
             cmgno32i = min(frmgi, 0.5d0*no3fri)
             frmgi    = max(frmgi - cmgno32i, 0.0d0)
             no3fri   = max(no3fri - 2.0d0*cmgno32i, 0.0d0)
             cmgcl2i  = min(frmgi, 0.5d0*clfri)
             clfri    = max(clfri - 2.0d0*cmgcl2i, 0.0d0)
             rest3    = 2.0d0*frso4i + no3fri + clfri
    !
             if (ca > rest) then         ! Ca > 2*SO4 + Cl + NO3?
                frca_init = max(ca - (1.0 - 1.0d-6)*rest, 0.0d0)
                frna_init = na
                frk_init  = k
                frmg_init = mg
                ca        = (1.0 - 1.0d-6)*rest
                na        = 0.0d0
                k         = 0.0d0
                mg        = 0.0d0
    !            write(*,*), 'Warning error: Ca, Na, K, Mg in excess'
    !
             else if (na > rest1) then  ! Na > 2*frso4 + frcl + frno3?
                frna_init = max(na - (1.0 - 1.0d-6)*rest1, 0.0d0)
                frk_init  = k
                frmg_init = mg
                na        = (1.0 - 1.0d-6)*rest1
                k         = 0.0d0
                mg        = 0.0d0
    !            write(*,*), 'Warning error: Na, K, Mg in excess'
    !
             else if (mg > rest2) then  ! Mg > 2*frso4 + frcl + frno3?
                frk_init  = k
                frmg_init = max(mg - (1.0 - 1.0d-6)*rest2, 0.0d0)
                mg        = (1.0 - 1.0d-6)*rest2
                k         = 0.0d0
    !            write(*,*), 'Warning error: K, Mg in excess'
    !
             else if (k > rest3) then   ! K > 2*frso4 + frcl + frno3?
                frk_init  = max(k - (1.0 - 1.0d-6)*rest3, 0.0d0)
                k         = (1.0 - 1.0d-6)*rest3
    !            write(*,*), 'Warning error: K in excess'
             end if
          end if
    !
    !  ## Calculate sulfate ratios
          so4rat  = (na + nh4 + ca + k + mg) / so4
          crnarat = (na + ca + k + mg) / so4
          crrat   = (ca + k + mg) / so4
    !
    !  ## 11. Sulfate poor and dust + sodium poor
          if (2.0 <= so4rat .and. crnarat < 2.0) then
             case_number = 11.0d0
             call mach_hetp_calco7(so4, nh4, nh3, hno3, hcl, hso4,      &
                                   na, cl, no3, h, lwc, ca, k, mg,      &
                                   caso4, frmg, frna, frca, frk, rh,    &          
                                   temp, k0, p1, p2, nr)
    !
    !  ## 12-13. Sulfate poor and dust + sodium rich
          else if (so4rat >= 2.0 .and. crnarat >= 2.0) then
             if (crrat <= 2.0) then
                case_number = 12.0d0
                call mach_hetp_calcm8(so4, nh4, nh3, hno3, hcl, hso4,  &
                                      na, cl, no3, h, lwc, ca, k, mg,  &
                                      caso4, frca, frmg, frk, frna,    &
                                      rh, temp, k0, p1, p2, nr)
    !
             else if (crrat > 2.0) then
                case_number = 13.0d0
                call mach_hetp_calcp13(so4, nh4, nh3, hno3, hcl, hso4, &
                                       na, cl, no3, h, lwc, ca, k, mg, &
                                       caso4, frso4, frmg, frk, frca,  &
                                       frna, rh, temp, k0, p1, p2, nr)
             end if
    !
    !  ## 14. Sulfate rich (no acid)
          else if (1.0 <= so4rat .and. so4rat < 2.0) then
             case_number = 14.0d0
             call mach_hetp_calcl9(so4, nh4, nh3, hno3, hcl, hso4,     &
                                   na, cl, no3, h, lwc, ca, k, mg,     &
                                   caso4, frmg, frk, frca, frna,       &
                                   frso4, rh, temp, k0, p1, p2, nr)
    !
    !  ## 15. Sulfate super rich (free acid)
          else if (so4rat < 1.0) then
             case_number = 15.0d0
             call mach_hetp_calck4(so4, nh4, nh3, hno3, hcl, hso4,     &
                                   na, cl, no3, h, lwc, ca, k, mg,     &
                                   caso4, rh, temp, k0, p1, p2, nr)
          end if
       end if
    !
    !  ### Return output ###
    !  Add up free amounts after chemical partitioning at thermodynamic equilibrium
       frca  = frca + frca_init
       frna  = frna + frna_init
       frk   = frk  + frk_init
       frmg  = frmg + frmg_init
       lwc   = lwc * 1.0d3 / 18.0d0   ! kg/m3 air -> mol/m3 air
    !
       return
    end subroutine mach_hetp_main_15cases
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: A2; Sulfate poor
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calca2(so4_i, nh4_i, nh3g_i, hso4_i, h_i, oh_i,                &
                                lwn_i, rh, temp, k0, p1, p2, nr) 
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr   
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: oh_i    
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables
       real(kind=8)     :: so4, nh4, hso4, gnh3, h, oh, lwn, t, aw, so4_t, nh4_t
       real(kind=8)     :: khso4, knh3, kh2o, errin, errouloc, errinlocb
       real(kind=8)     :: omebe, omehi, y1, y2, y3, x3, c1, k1, k2, k3, k4, ya, yb, xa, xb
       real(kind=8)     :: dx, tt1, tt2, lwnsq, gmax, loccondition
       real(kind=8)     :: nh, sigma, xt, xf, xh, delta, rr, gx, gx2, u1
       integer(kind=4)  :: j, k, rooteval, irh, nmax
       logical(kind=4)  :: condition, noroot, earlyexit, soln, frst, calain, calou
       real(kind=8), dimension(13) :: gama, gamin, gamou
    !
       so4   = so4_i
       nh4   = nh4_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       gnh3  = 0.0d0
       h     = 0.0d0
       oh    = 0.0d0
       lwn   = 0.0d0
       so4_t = 0.0d0
       nh4_t = nh4
       noroot    =.false.
       earlyexit = .false.
       soln      = .false.
       frst      = .true.
       calain    = .true.
       calou     = .true.
       loccondition = 0.0d0
       gama  = 0.1d0
       gamou = 0.1d0
       gamin = 1.0d10
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code 
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt1 = tstd / t - 1.0d0
       tt2 = 1.0d0 + log(tstd / t) - tstd / t
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                           &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    ! 
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !  
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Aerosol liquid water content 
       so4_t = so4
       lwn   = max(so4_t/awas(irh), tiny)
    !
    !  ## Constant values
       lwnsq = aw*lwn*lwn 
       c1    = r*t
       k1    = kh2o*lwnsq
       k2    = knh3*c1
       k3    = khso4*lwn
       k4    = k2/kh2o
    !
    !
    !  ### STAGE 1: Root tracking ###
    !  ## Find a subinterval [xa,xb] on the larger interval [I1,I2] where a sign change occurs
       rooteval = 0
       condition = .true.
       do while (rooteval < 2 .or. (condition .and. rooteval < ndiv+1))
          rooteval = rooteval + 1
    !
    !  ## Set low limit and high limit
          if (rooteval == 1) then
             omehi = 2.0d0*so4   ! high limit: from NH4+ -> NH3(g) + H+(aq)
             y1 = 1.0d0
          end if 
    !      
    !  ## Begin search on subinterval 
          if (rooteval == 2) then
             soln = .false.
             if (abs(y2) <= eps) then
                earlyexit = .true.
             end if 
             y1 = y2
    !
             if (earlyexit) then
                dx = 0.0d0
             else
                dx = (omehi - tiny)/float(ndiv) 
             end if 
    !
             omehi = max(omehi - dx, tiny) 
          end if 
    !
    !  ## Continue search
          if (rooteval > 2) then
             if (loccondition < 0.0d0) then
    !  ## 1. Root has been found on the subinterval; save x values for ITP search
                y1    = y1
                omebe = omebe
                omehi = omehi
             else
    !  ## 2. No root has been found; continue searching in the next subinterval
                y1    = y2
                omebe = omehi
                omehi = max(omehi - dx, tiny)
             end if 
          end if 
    !      
    !  ## Solve the system of equations
          frst   = .true.
          calain = .true.
    !
          if (( .not. soln)) then
             h     = omehi
             so4_t = max(so4 - (so4/((k3/gama(7)*(gama(8)/gama(7))**2.0)/omehi + 1.0d0)), tiny) 
             nh4_t = max(nh4/(1.0d0/(k4*(gama(8)/gama(9))**2.0)/omehi + 1.0d0), 2.0d0*so4_t)
             hso4  = so4/((k3/gama(7)*(gama(8)/gama(7))**2.0)/omehi + 1.0d0)
             gnh3  = max(nh4 - nh4_t, tiny)
             oh    = k1/omehi
          end if 
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)         
             k = k + 1
    
    !  ## Reset gamou
             if (( .not. soln) .and. frst) then
                gamou(4)  = gama(4)
                gamou(7)  = gama(7)
                gamou(8)  = gama(8)
                gamou(9)  = gama(9)
                gamou(13) = gama(13)
             end if 
    
    !  ## Reset gamin
             if (( .not. soln)) then
                gamin(4)  = gama(4)
                gamin(7)  = gama(7)
                gamin(8)  = gama(8)
                gamin(9)  = gama(9)
                gamin(13) = gama(13)
             end if 
    !
             call mach_hetp_calcact1b(h, nh4_t, so4_t, hso4, lwn, gama, t, soln,   &
                                      frst, calain, calou)
    !
             if (frst) then
                errouloc = 0
                errouloc = max(errouloc, abs(gamou(7)  - gama(7))  / gamou(7))
                errouloc = max(errouloc, abs(gamou(8)  - gama(8))  / gamou(8))
                errouloc = max(errouloc, abs(gamou(9)  - gama(9))  / gamou(9))
                errouloc = max(errouloc, abs(gamou(4)  - gama(4))  / gamou(4))
                errouloc = max(errouloc, abs(gamou(13) - gama(13)) / gamou(13))
                calou = errouloc .ge. epsact
                frst  = .false.
             end if 
    !
             errinlocb = 0
             errinlocb = max(errinlocb, abs(gamin(7)  - gama(7))  / gamin(7))
             errinlocb = max(errinlocb, abs(gamin(8)  - gama(8))  / gamin(8))
             errinlocb = max(errinlocb, abs(gamin(9)  - gama(9))  / gamin(9))
             errinlocb = max(errinlocb, abs(gamin(4)  - gama(4))  / gamin(4))
             errinlocb = max(errinlocb, abs(gamin(13) - gama(13)) / gamin(13))
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients 
             errin = max(errin, 0.0d0)
             errin = max(errin, abs((gamin(7) - gama(7))  / gamin(7)))
             errin = max(errin, abs((gamin(8) - gama(8))  / gamin(8)))
             errin = max(errin, abs((gamin(9) - gama(9))  / gamin(9)))
             errin = max(errin, abs((gamin(4) - gama(4))  / gamin(4)))
             errin = max(errin, abs((gamin(13)- gama(13)) / gamin(13)))
    !
    !  ## Solve system of equations, using new activity coefficients 
             if (( .not. soln)) then
                h     = omehi
                so4_t = max(so4 - (so4/((k3/gama(7)*(gama(8)/gama(7))**2.0)/omehi + 1.0d0)), tiny)  
                nh4_t = max(nh4/(1.0d0/(k4*(gama(8)/gama(9))**2.0)/omehi + 1.0d0), 2.0d0*so4_t)
                hso4  = so4/((k3/gama(7)*(gama(8)/gama(7))**2.0)/omehi + 1.0d0)
                gnh3  = max(nh4 - nh4_t, tiny)
                oh    = k1/omehi
             end if
          end do
    !
          y2 = (nh4_t/(2.0d0*so4_t + hso4) - 1.0d0) + h/(2.0d0*so4_t + hso4)
    !
    !  ## Check for criteria to exit root tracking 
          loccondition = sign(1.0d0,y1)*sign(1.0d0,y2)
    !
    !  ## AFTER iterating through ALL ndiv subdivided intervals
          if (rooteval == ndiv + 1) then
    !  ## (1) No solution; reset x-value to 'tiny' and use to solve system
             if (loccondition > 0.0d0 .and. abs(y2) > eps .and. ( .not. earlyexit)) then
                noroot = .true.
    !           write(*,*), 'Warning in CALCA2: no solution found; no interval with sign change'
                omehi = tiny   ! Reset to tiny 
                omebe = omehi
             else if (loccondition > 0.0d0 .and. abs(y2) <= eps .and. ( .not. earlyexit)) then
    !  ## (2) Solution is assumed and ITP is not required
                noroot = .true.
                soln   = .true.
             else if (earlyexit) then
    !  ## (3) Solution is assumed and ITP is not required
                noroot = .true.
                omehi = tiny    ! Reset to tiny 
                omebe = omehi
             end if
          end if 
    !
    !  ## Test for criteria to exit root tracking 
          condition = .false.
          if (loccondition > 0.0d0 .and. ( .not. noroot)) then
             condition = .true.
          else
             soln = .true.
          end if 
       end do
    !
    !  ## The root tracking did not iterate until ndiv
    !  ## Check if an x-value of tiny was a solution, if so reset x-value to tiny
       if (rooteval < ndiv + 1 .and. earlyexit) then
    !  Solution has been assumed 
          noroot = .true.
          omehi  = tiny    ! Reset to tiny 
          omebe  = omehi
          soln   = .true.
       end if 
    !
       if (( .not. earlyexit) .and. ( .not. noroot)) then
          soln = .false.
       end if
    !
    !
    !  ### STAGE 2: modified bisection search (using ITP algorithm) ###
    !  ## Initialize static ITP variables 
       if (( .not. noroot)) then
          yb = y1
          ya = y2
    !
          xa = omehi
          xb = omebe
          x3 = omebe
    !
          if (xa == xb) then
             noroot = .true.
             gx = tiny        
          else
             gx = xb - xa
          end if 
    !
          gx2  = (xa+xb)*0.5d0
          u1   = 0.2d0/gx
          nh   = log10(abs(gx/(2.0d0*eps*gx2))) / log10of2
          nmax = int(nh) + 2 
       else
          x3   = omehi
          soln = .true.
       end if 
    !
    !  ## Start search
       j          = 0
       condition = .true.
       do while (j < maxit .and. condition)
    !  ## Set dynamic ITP variables 
          if (( .not. noroot) .and. ( .not. soln)) then
             gx   = xb - xa
             xh   = 0.5d0*(xa + xb)
             rr   = max(gx2*eps*2.d0**(real(nmax - j)) - 0.5d0*gx, 0.0d0)
             delta= u1*(max(gx, 0.0d0))**2.0d0   
             xf   = max((yb*xa - ya*xb) / (yb - ya), 0.0d0)
    !
             sigma = sign(1.0d0, xh - xf)
             if (delta <= abs(xh - xf)) then
                xt = xf + sigma*delta               
             else
                xt = xh
             end if
    !
             if (abs(xt - xh) <= rr) then
                x3 = xt
             else
                x3 = xh - sigma*rr
             end if 
          end if       
    !
          j = j + 1
    !      
          if (( .not. soln)) then
             gmax = 0.1d0
             gmax = max(gmax, gama(7))
             gmax = max(gmax, gama(8))
             gmax = max(gmax, gama(9))
             gmax = max(gmax, gama(4))
             gmax = max(gmax, gama(13))
          end if 
    !
    !  ## Reinitialize activity coefficients if gmax > 100.0d0
          if (gmax > 100.0d0 .and. ( .not. soln)) then
             gama(7)  = 0.1d0
             gama(8)  = 0.1d0
             gama(9)  = 0.1d0
             gama(4)  = 0.1d0
             gama(13) = 0.1d0
             gamin(7) = 1.0d10
             gamin(8) = 1.0d10
             gamin(9) = 1.0d10
             gamin(4) = 1.0d10
             gamin(13)= 1.0d10
             calou    = .true.
             frst     = .true.
          end if 
    !    
    !  ## Solve system of equations
          frst    = .true.
          calain  = .true.
          if (( .not. soln)) then
             h     = x3
             so4_t = max(so4 - (so4/((k3/gama(7)*(gama(8)/gama(7))**2.0)/x3 + 1.0d0)), tiny) 
             nh4_t = max(nh4/(1.0d0/(k4*(gama(8)/gama(9))**2.0)/x3 + 1.0d0), 2.0d0*so4_t)
             hso4  = so4/((k3/gama(7)*(gama(8)/gama(7))**2.0)/x3 + 1.0d0)
             gnh3  = max(nh4 - nh4_t, tiny)
             oh    = k1/x3
          end if
    !
    !  ## Iterate until convergence of activity coefficients 
          k = 0
          errin = 1.0d0
          do while ( k < nsweep-1 .and.  errin >= epsact)
             k = k + 1
    !
    !  ## Reset gamou
             if (( .not. soln) .and. frst) then
                gamou(4)  = gama(4)
                gamou(7)  = gama(7)
                gamou(8)  = gama(8)
                gamou(9)  = gama(9)
                gamou(13) = gama(13)
             end if 
    !
    !  ## Reset gamin
             if (( .not. soln)) then
                gamin(4)  = gama(4)
                gamin(7)  = gama(7)
                gamin(8)  = gama(8)
                gamin(9)  = gama(9)
                gamin(13) = gama(13)
             end if 
    !
             call mach_hetp_calcact1b(h, nh4_t, so4_t, hso4, lwn, gama, t, soln,   &
                                      frst, calain, calou)
    !
             if (frst) then
                errouloc = 0
                errouloc = max(errouloc, abs(gamou(7)  - gama(7))  / gamou(7))
                errouloc = max(errouloc, abs(gamou(8)  - gama(8))  / gamou(8))
                errouloc = max(errouloc, abs(gamou(9)  - gama(9))  / gamou(9))
                errouloc = max(errouloc, abs(gamou(4)  - gama(4))  / gamou(4))
                errouloc = max(errouloc, abs(gamou(13) - gama(13)) / gamou(13))
                calou    = errouloc .ge. epsact
                frst     = .false.
             end if 
    !
             errinlocb = 0
             errinlocb = max(errinlocb, abs(gamin(7)  - gama(7))  / gamin(7))
             errinlocb = max(errinlocb, abs(gamin(8)  - gama(8))  / gamin(8))
             errinlocb = max(errinlocb, abs(gamin(9)  - gama(9))  / gamin(9))
             errinlocb = max(errinlocb, abs(gamin(4)  - gama(4))  / gamin(4))
             errinlocb = max(errinlocb, abs(gamin(13) - gama(13)) / gamin(13))
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients 
             errin = max(errin, 0.0d0)
             errin = max(errin, abs((gamin(7) - gama(7))  / gamin(7)))
             errin = max(errin, abs((gamin(8) - gama(8))  / gamin(8)))
             errin = max(errin, abs((gamin(9) - gama(9))  / gamin(9)))
             errin = max(errin, abs((gamin(4) - gama(4))  / gamin(4)))
             errin = max(errin, abs((gamin(13)- gama(13)) / gamin(13)))
    !
    !  ## Solve system of equations, with new activity coefficients 
             if (( .not. soln)) then
                h     = x3
                so4_t = max(so4 - (so4/((k3/gama(7)*(gama(8)/gama(7))**2.0)/x3 + 1.0d0)), tiny) 
                nh4_t = max(nh4/(1.0d0/(k4*(gama(8)/gama(9))**2.0)/x3 + 1.0d0), 2.0d0*so4_t)
                hso4  = so4/((k3/gama(7)*(gama(8)/gama(7))**2.0)/x3 + 1.0d0)
                gnh3  = max(nh4 - nh4_t, tiny)
                oh    = k1/x3
             end if 
          end do
    !     
          y3 = (nh4_t/(2.0d0*so4_t + hso4) - 1.0d0) + h/(2.0d0*so4_t + hso4) 
    !     
          condition = .false.
          if (noroot) then
    !  ## If no root on interval then do not perform ITP
             xa = x3
             xb = x3
          else if (y3 > 0.0d0 .and. ( .not. soln)) then
             xb = x3
             yb = y3
          else if (y3 < 0.0d0 .and. ( .not. soln)) then
             xa = x3
             ya = y3
          else if (( .not. soln)) then
             xa = x3
             xb = x3
          end if
    !
    !  ## Check for convergence criteria to exit ITP:
          if (abs(xb - xa) > abs(xa*eps) .and. ( .not. noroot)) then
             condition = .true.
             soln   = .false.
          else
             soln   = .true.
          end if
       end do
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       hso4_i  = hso4
       h_i     = h
       nh3g_i  = gnh3
       lwn_i   = lwn
       oh_i    = oh
    !
       return
    end subroutine mach_hetp_calca2
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: B4; Sulfate rich, no free acid
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcb4(so4_i, nh4_i, nh3g_i, hso4_i, h_i, lwn_i,                &
                                rh, temp, k0, p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables:
       real(kind=8)     :: so4, nh4, hso4, gnh3, h, lwn, t, aw, tt1, tt2, c1, so4_t, nh4_t
       real(kind=8)     :: khso4, knh3, kh2o, bb, cc, dd, hh, ff, v, ak1, errin, tt0
       real(kind=8)     :: m4, m9, m13, a3, c2, c3, c4
       integer(kind=4)  :: j, irh
       logical(kind=4)  :: gg
       real(kind=8), dimension(13) :: gama, gamin 
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       gnh3  = 0.0d0
       h     = 0.0d0
       lwn   = 0.0d0
       so4_t = 0.0d0
       nh4_t = nh4
       gama  = 0.1d0
       gamin = 1.0d10
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*          &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    ! 
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Calculate dry composition; used to calculate initial aerosol liquid water
       gg = 2.0*so4 - nh4 <= nh4 - so4 
    !
       if (gg) then   
    !  NH4HSO4 >= (NH4)2SO4
          m4  = 2.0d0*nh4 - 3.0d0*so4   ! clc
          m9  = 0.0d0                   ! nh4hs4
          m13 = 2.0d0*so4 - nh4         ! nh42s4
       else 
    !  NH4HSO4 < (NH4)2SO4
          m4  = 0.0d0                   ! clc
          m9  = 3.0d0*so4 - 2.0d0*nh4   ! nh4hs4
          m13 = nh4 - so4               ! nh42s4
       end if
    !
    !  ## Save ZSR parameters as variables to limit indirect addressing in loops
       c2   = awlc(irh)
       c3   = awab(irh)
       c4   = awas(irh)
    !
    !  ## Initial aerosol liquid water content
       lwn  = max(m13/c2 + m9/c3 + m4/c4, tiny)
       c1   = khso4*(lwn/0.1d0)
    ! 
    !
    !  ### MAJOR SYSTEM H+/HSO4-/SO42- ###
    !  ## Setup initial conditions
    !  ## 1. Calculate initial concentrations of H+/HSO4-/SO42- 
       bb  = so4 + c1 - nh4  
       cc  = -c1*so4
    !
       if (bb /= 0.d0) then
          dd = cc / (bb*bb)
          v  = 4.d0*dd
       else
          v  = 1.d+03
       end if       
    !
       if (abs(v) <= smrt .and. bb /= 0.d0) then
          hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
       else
          hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))                            !positive root
       end if 
    !
    !  ## 2. Speciation  
       so4_t = max(tiny, min(hh, so4))                  
       hso4  = max(tiny, min(so4 - so4_t, so4)) 
       h     = max(tiny, min(c1*hso4/so4_t, so4))
    !
    !  ## 3. Aerosol liquid water content  
    !  ## Correct for HSO4 dissociation      
       gg = so4_t - h < hso4 + h
       if (gg) then
          m4   = 0.0d0
          m13  = so4_t - h
          m9   = max((hso4 + h) - (so4_t - h), 0.0d0)
       else
          m13  = hso4 + h
          m9   = 0.0d0
          m4   = max((so4_t - h) - (hso4 + h), 0.0d0)
       end if 
       lwn  = max(m13/c2 + m9/c3 + m4/c4, tiny)
    !
    !  ## Iterative search for solution with convergence of activity coefficients
       j     = 0
       errin = 1.0d0
       do while (j < nsweep-1 .and. errin >= epsact)
          j = j + 1
    !
    !  ## Reset gamin
          gamin(4)  = gama(4)
          gamin(7)  = gama(7)
          gamin(8)  = gama(8)
          gamin(9)  = gama(9)
          gamin(13) = gama(13)
    !
          call mach_hetp_calcact1(h, nh4_t, so4_t, hso4, lwn, gama, t)
    !
    !  ## Test convergence criterion: max change in the activity coefficients
          errin = 0.0d0
          errin = max(errin, abs((gamin(7) - gama(7))  / gamin(7)))
          errin = max(errin, abs((gamin(8) - gama(8))  / gamin(8)))
          errin = max(errin, abs((gamin(9) - gama(9))  / gamin(9)))
          errin = max(errin, abs((gamin(4) - gama(4))  / gamin(4)))
          errin = max(errin, abs((gamin(13)- gama(13)) / gamin(13)))
    
    !  ## 1. Solve system of equations
          ak1 = khso4*((gama(8)/gama(7))*(gama(8)/gama(7)))*(lwn/gama(7))
          bb  = so4 + ak1 - nh4_t
          cc  = -ak1*so4 
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v  = 4.d0*dd
          else
             v = 1.d+03
          end if       
    !
          if (abs(v) <= smrt .and. bb /= 0.d0) then
             hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
          else
             hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))   !positive root
          end if 
    !
    !  ## 2. Speciation
          so4_t = max(tiny, min(hh, so4))
          hso4  = max(tiny, min(so4 - so4_t, so4))
          h     = max(tiny, min(ak1*hso4/so4_t, so4))
    !
    !  ## 3. Aerosol liquid water content
    !  ## Correct for HSO4 dissociation 
          gg = so4_t - h < hso4 + h
          if (gg) then   
             m4  = 0.0d0
             m13 = so4_t - h
             m9  = max((hso4 + h) - (so4_t - h), 0.0d0)  ! NH4HSO4
          else 
             m13 = hso4 + h
             m9  = 0.0d0
             m4  = max((so4_t - h) - (hso4 + h), 0.0d0)  ! (NH4)2SO4
          end if 
    !
          lwn = max(m13/c2 + m9/c3 + m4/c4, tiny)
       end do
    !
    !
    !  ### MINOR SYSTEM: NH4+/NH3/H+ ###
    !  Ammonia in the gas phase is assumed a minor species that does not significantly perturb 
    !  the aerosol equilibrium: NH3(g) + H+(aq) <==> NH4+(aq); H+ determined from the major
    !  system (solved above) is used initially.  
       if (lwn > tiny) then
    !  ## 1. Calculate NH3 sublimation
          ff = 0.0d0
          a3 = (knh3/kh2o)*r*t*(gama(10)/gama(5))**2.0
          bb = h + 1.0d0/a3    ! bb always > 0
          cc = -nh4_t/a3
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d+03
          end if  
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             ff = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             ff = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
          end if 
    !
    !  ## 2. Speciation 
    !  ## Due to round-off, ff may be .gt. nh4_t giving negative nh4_t
    !  ## If this condition is true, then set ff = nh4_t
          ff    = max(tiny, min(ff, nh4_t))
          gnh3  = ff
          nh4_t = max(nh4_t - ff, 0.0d0)
          h     = h + ff
       end if
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       hso4_i  = hso4
       h_i     = h
       nh3g_i  = gnh3
       lwn_i   = lwn
    !
       return
    end subroutine mach_hetp_calcb4
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: C2; Sulfate rich, free acid
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcc2(so4_i, nh4_i, nh3g_i, hso4_i, h_i, lwn_i,               &
                                rh, temp, k0, p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables
       real(kind=8)      :: so4, nh4, hso4, gnh3, h, lwn, frh2so4, t, aw
       real(kind=8)      :: khso4, knh3, kh2o, bb, cc, dd, hh, v, errin, tt0
       real(kind=8)      :: so4_t, nh4_t, tt1, tt2, c1, c2, c3, ff, a3
       integer(kind=4)   :: j, irh
       real(kind=8), dimension(13) :: gama, gamin
    !
    !  ### Initialize variables ### 
       so4   = so4_i     
       nh4   = nh4_i   
       aw    = rh         
       t     = temp      
       hso4  = 0.0d0      
       gnh3  = 0.0d0      
       h     = 0.0d0      
       lwn   = 1.0d-20   
       so4_t = 0.0d0      
       nh4_t = 0.0d0      
       gama  = 0.1d0
       gamin = 1.0d10
       frh2so4 = 0.0d0
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    ! 
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Aerosol liquid water content 
       c3  = so4 - nh4
       lwn  = max(max(c3, 0.0d0)/awsa(irh) + nh4/awab(irh), tiny)
    !
    !  ## Constant values
       c1  = khso4*1.0d-19    
       c2  = khso4*lwn
    !
    !
    !  ### MAJOR SYSTEM H+/HSO4-/SO42- ###
    !  ## Setup initial conditions
    !  ## 1. Calculate initial concentrations of H+/HSO4-/SO42-
       bb = c3 + c1
       cc = -c1*so4
    !
       if (bb /= 0.d0) then
          dd = cc / (bb*bb)
          v  = 4.d0*dd
       else
          v  = 1.d+03
       end if       
    !
       if (abs(v) <= smrt .and. bb /= 0.d0) then
          hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
       else
          hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))                            !positive root
       end if     
    !
    !  ## 2. Speciation  
       nh4_t   = nh4
       so4_t   = hh                
       hso4    = max(so4 - hh, tiny) 
       h       = c3 + hh
       frh2so4 = max(so4_t + hso4 - nh4_t, 0.0d0)   ! Free H2SO4
    !  
    !  ## Iterative search for solution with convergence of activity coefficients
       j     = 0
       errin = 1.0d0
       do while (j < nsweep-1 .and. errin >= epsact)
          j = j + 1
    !
    !  ## Reset gamin
          gamin(4)  = gama(4)
          gamin(7)  = gama(7)
          gamin(8)  = gama(8)
          gamin(9)  = gama(9)
          gamin(13) = gama(13)
    !
          call mach_hetp_calcact1(h, nh4_t, so4_t, hso4, lwn, gama, t)
    !
    !  ## Test convergence criterion: max change in the activity coefficients
          errin = 0.0d0
          errin = max(errin, abs((gamin(7) - gama(7))  / gamin(7)))
          errin = max(errin, abs((gamin(8) - gama(8))  / gamin(8)))
          errin = max(errin, abs((gamin(9) - gama(9))  / gamin(9)))
          errin = max(errin, abs((gamin(4) - gama(4))  / gamin(4)))
          errin = max(errin, abs((gamin(13)- gama(13)) / gamin(13)))
    !
    !  ## 1. Solve system of equations
          bb = c3 + c2/gama(7)*(gama(8)/gama(7))**2.0
          cc = -(c2/gama(7)*(gama(8)/gama(7))**2.0)*so4
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v  = 4.d0*dd
          else
             v = 1.d+03
          end if       
    !
          if (abs(v) <= smrt .and. bb /= 0.d0) then
             hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
          else
             hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))   !positive root
          end if 
    !
    !  ## 2. Speciation
          so4_t   = hh                
          hso4    = max(so4 - hh, tiny) 
          h       = c3 + hh
          frh2so4 = max(so4_t + hso4 - nh4_t, 0.0d0)    
       end do
    !
    !
    !  ### MINOR SYSTEM: NH4+/NH3/H+ ###
    !  Ammonia in the gas phase is assumed a minor species that does not significantly perturb 
    !  the aerosol equilibrium: NH3(g) + H+(aq) <==> NH4+(aq); H+ determined from the major
    !  system (solved above) is used initially.
       if (lwn > tiny) then
    !  Calculate NH3 sublimation
          ff = 0.0d0
          a3 = (knh3/kh2o)*r*t*(gama(10)/gama(5))**2.0
          bb = h + 1.0d0/a3    ! bb always > 0
          cc = -nh4_t/a3
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d+03
          end if  
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             ff = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             ff = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
          end if 
    !
    !  Speciation 
    !  Due to round-off, ff may be .gt. nh4_t giving negative nh4_t
    !  If this condition is true, then set ff = nh4_t
          ff    = max(tiny, min(ff, nh4_t))
          gnh3  = ff
          nh4_t = max(nh4_t - ff, 0.0d0)
          h     = h + ff
       end if
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       hso4_i  = hso4
       h_i     = h
       nh3g_i  = gnh3
       lwn_i   = lwn
    !
       return
    end subroutine mach_hetp_calcc2
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: D3; Sulfate poor
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcd3(so4_i, nh4_i, hno3g_i, nh3g_i, hso4_i, h_i, no3_i,              &
                                lwn_i, rh, temp, k0, p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables
       real(kind=8)     :: so4, nh4, hso4, ghno3, gnh3, gnh3_i, ghno3_i, no3, h, lwn, t, aw, nh4no3, nh42s4
       real(kind=8)     :: knh3, khso4, kh2o, khno3, knh4no3, k1, k2, k3, c1, c1a, c1b, c2, c3, c4, c10, c11
       real(kind=8)     :: bb, cc, dd, hh, v, ak1, errin, tt0, a3, a7, psi3, denm, m1, m2
       real(kind=8)     :: so4_t, nh4_t, no3_t, tt1, tt2, x3, gmax, y3, psilo, psihi, k4, k5, k6
       real(kind=8)     :: omehi, omebe, y1, y2, dx, a4, ylo, yhi, pshi, ya, yb, xa, xb
       real(kind=8)     :: nh, sigma, xt, xf, xh, delta, rr, gx, gx2, u1, errouloc, errinlocb
       real(kind=8)     :: cl, cl_t, ghcl, caso4, gama10sq, lwnsq
       integer(kind=4)  :: ii, j, k, rooteval, count_rt, irh, nmax
       logical(kind=4)  :: condition, soln, earlyexit, frst, calain, calou, noroot, redo, earlye
       real(kind=8),  dimension(13) :: gama, gamin, gamou, gama_min
       real(kind=8)     :: y3_min, x3_min, y3_lastiter, x3_lastiter
       real(kind=8)     :: no3_min, nh4_min, h_min, lwn_min, gnh3_min, ghno3_min
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       ghno3 = 0.0d0
       gnh3  = 0.0d0
       gnh3_i= 0.0d0
       ghno3_i=0.0d0
       h     = 0.0d0
       lwn   = 0.0d0
       nh4no3= 0.0d0
       nh42s4= 0.0d0
       so4_t = 0.0d0
       nh4_t = nh4
       no3_t = 0.0d0
       cl    = 0.0d0
       cl_t  = 0.0d0
       ghcl  = 0.0d0
       caso4 = 0.0d0
       soln  = .false.
       frst  = .true.
       calain= .true.
       calou = .true.
       redo  = .false.
       noroot= .false.
       count_rt = 1
       earlyexit = .false.
       gama  = 0.1d0
       gamou = 0.1d0
       gamin = 1.0d10
       earlye = .false.
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    !  
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !   
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3= k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. NH4NO3(s) <==> NH3(g) + HNO3(g)                            (xk10)
       knh4no3 = k0(9) * exp(p1(9)*tt1 + p2(9)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Calculate NH4NO3 that volatizes (ISORROPIA II: subroutine 'CALCD1A')
       nh42s4 = so4
       c1     = max(0.0d0, min(nh4 - 2.0d0*nh42s4, no3))
       c1a    = max(no3 - c1, 0.0d0)
       c1b    = max(nh4 - c1 - 2.0d0*nh42s4, 0.0d0)
       c2     = c1a + c1b
       c3     = sqrt(c2*c2 + 4.0d0*(knh4no3/(r*t)/(r*t)))
       c4     = min(c1, 0.5d0*(-c2 + c3))
    !
    !  ## Initial speciation
       nh4no3 = max(c1 - c4, 0.0d0)   ! NH4NO3 (s)
       gnh3   = c1b + c4              ! NH3 (g)
       ghno3  = c1a + c4              ! HNO3 (g)
       gnh3_i = gnh3
       ghno3_i= ghno3
       so4_t  = nh42s4
       hso4   = 0.0d0
       nh4_t  = nh4no3
       no3_t  = nh4no3
    ! 
    !  ## Initial aerosol liquid water content 
       m1  = (so4_t + hso4)       ! (NH4)2SO4
       m2  = (nh4_t - 2.0d0*m1)   ! free NH4 
    ! 
    !  ## Save ZSR allocations to limit indirect addressing  
       c10 = awas(irh) 
       c11 = awan(irh)
       lwn = max(m1/c10 + max(min(m2, no3_t), 0.0d0)/c11, tiny)
    !
    !  ## Constant values
       k1  = khno3*r*t
       k2  = (knh3/kh2o)*r*t
       k3  = kh2o*aw 
       pshi= gnh3
       k4  = nh4no3*2.0d0*nh42s4 + nh4no3*nh4no3
       k5  = 2.0d0*nh42s4 + nh4no3        
    !
    !
    !  ### STAGE 1: Root tracking ###
    !  ## Find a subinterval [xa,xb] on the larger interval [I1,I2] where a sign change occurs
       rooteval = 0
       condition = .true.
       do while (rooteval < 2 .or. (condition .and. rooteval < ndiv + 1))   ! Begin outer loop for root tracking
          rooteval = rooteval + 1
    !
    !  ## Set high limit for root tracking (i.e. lower bound)
          if (rooteval == 1 .and. count_rt == 1) then
            omehi = tiny   
            y1    = 1.0d0
          end if
    !
    !  ## Begin search on subinterval 
          if (rooteval == 2) then
             soln = .false.
             if (count_rt == 1 .or. redo) then
                ylo = y2
             end if 
    !
             gmax = 0.1d0
             gmax = max(gmax, gama( 4))
             gmax = max(gmax, gama( 5))
             gmax = max(gmax, gama( 7))
             gmax = max(gmax, gama( 8))
             gmax = max(gmax, gama( 9))
             gmax = max(gmax, gama(10))
             gmax = max(gmax, gama(13))
    !
    !  ## Reinitialize activity coefficients if gmax > 100.0d0
             if (gmax > 100.0d0 .and. ( .not. soln)) then
                gama  = 0.1d0
                gamin = 1.0d10
                gamou = 1.0d10
                calou = .true.
                frst  = .true.
             end if
    !
             if (abs(y2) <= eps .and. (count_rt == 1 .or. redo)) then
                earlyexit = .true.
             end if 
    !
             y1 = y2
    !
             if (earlyexit .or. (count_rt > 1 .and. ( .not. redo))) then
                dx = 0.0d0
             else if (count_rt > 1 .and. redo) then
                dx = (max(pshi-omehi, tiny))/float(ndiv)
                omebe = omehi               ! Lower bound of subinterval (xa)
             elseif (count_rt == 1) then
                dx = (pshi-tiny)/float(ndiv)
                omebe = omehi               ! Lower bound of subinterval (xa)
             else
                dx = 0.0d0
             end if
    !
             omehi = omehi + dx             ! Upper bound of subinterval (xb)
          end if    
    !
    !  ## Continue search 
          if (rooteval > 2) then
             if (( .not. soln)) then
                gmax = 0.1d0
                gmax = max(gmax, gama( 4))
                gmax = max(gmax, gama( 5))
                gmax = max(gmax, gama( 7))
                gmax = max(gmax, gama( 8))
                gmax = max(gmax, gama( 9))
                gmax = max(gmax, gama(10))
                gmax = max(gmax, gama(13))
             end if
    !
    !  ## Reinitialize activity coefficients if gmax > 100.0d0
             if (gmax > 100.0d0 .and. ( .not. soln)) then
                gama  = 0.1d0
                gamin = 1.0d10
                gamou = 1.0d10
                calou = .true.
                frst  = .true.
             end if
    
             if (sign(1.0d0,y1)*sign(1.0d0,y2) < 0.0d0 .and. (count_rt == 1 .or. redo)) then
    !  ## 1. Root has been found on the subinterval; save x values for ITP search
                y1    = y1
                omebe = omebe
                omehi = omehi
             elseif (count_rt == 1 .or. redo) then
    !  ## 2. No root has been found, continue searching in the next subinterval
                y1    = y2
                omebe = omehi
                omehi = omehi + dx
             end if 
          end if 
    !
    !  ## Solve the system of equations 
          frst   = .true.
          calain = .true.
          if (( .not. soln) .and. (count_rt == 1 .or. redo)) then
             lwnsq    = lwn*lwn
             gama10sq = gama(10)*gama(10)
             a3       = k1*lwnsq/gama10sq 
             a4       = k2*gama10sq/(gama(5)*gama(5))  
             a7       = k3*lwnsq                      
    !
    !  ## 1. Calculate dissociation quantities 
             k6   = a3*a4*(gnh3_i - omehi)
             psi3 = k6*ghno3_i - k4 - nh4no3*omehi
             psi3 = psi3 / (k6 + k5 + omehi)        
             psi3 = min(max(psi3, 0.0d0), ghno3_i)
    !
    !  ## 2. Calculate H+
             bb = omehi - psi3
    !  Use ISORROPIA II current scheme 
             denm = bb + sqrt(bb*bb + 4.0d0*a7)
             if (denm <= tiny) then
                denm = (bb + abs(bb)) + 2.0d0*a7/abs(bb)
             end if
             h = 2.0d0*a7/denm
    !
    !  ## 3. Speciation 
             nh4_t = omehi + k5    
             so4_t = nh42s4
             no3_t = psi3 + nh4no3
             ghno3 = max(ghno3_i - psi3, tiny2)
             gnh3  = max(gnh3_i - omehi, tiny2)
    !       
    !  ## 4. Aerosol liquid Water content
             m1  = so4_t + hso4       ! (NH4)2SO4
             m2  = nh4_t - 2.0d0*m1   ! free NH4   
             lwn = max(m1/c10 + max(min(m2, no3_t), 0.0d0)/c11, tiny)
          end if
    !
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)         
             k = k + 1
    !
    !  ## Reset gamou
             if (( .not. soln) .and. frst .and. (count_rt == 1 .or. redo)) then
                gamou( 4)  = gama( 4)
                gamou( 5)  = gama( 5)
                gamou( 7)  = gama( 7)
                gamou( 8)  = gama( 8)
                gamou( 9)  = gama( 9)
                gamou(10)  = gama(10)
                gamou(13)  = gama(13)
             end if 
    !
    !  ## Reset gamin
             if (( .not. soln) .and. (count_rt == 1 .or. redo)) then
                gamin( 4)  = gama( 4)
                gamin( 5)  = gama( 5)
                gamin( 7)  = gama( 7)
                gamin( 8)  = gama( 8)
                gamin( 9)  = gama( 9)
                gamin(10)  = gama(10)
                gamin(13)  = gama(13)
             end if 
    !
             call mach_hetp_calcact2b(h, nh4_t, so4_t, hso4, no3_t, lwn, gama, t, soln,   &
                                      frst, calain, calou)
    !
             if (frst) then
                errouloc = 0
                errouloc = max(errouloc, abs(gamou(4 ) - gama(4 )) / gamou(4 ))
                errouloc = max(errouloc, abs(gamou(5 ) - gama(5 )) / gamou(5 ))
                errouloc = max(errouloc, abs(gamou(7 ) - gama(7 )) / gamou(7 ))
                errouloc = max(errouloc, abs(gamou(8 ) - gama(8 )) / gamou(8 ))
                errouloc = max(errouloc, abs(gamou(9 ) - gama(9 )) / gamou(9 ))
                errouloc = max(errouloc, abs(gamou(10) - gama(10)) / gamou(10))
                errouloc = max(errouloc, abs(gamou(13) - gama(13)) / gamou(13))
                calou    = errouloc .ge. epsact
                frst     = .false.
             end if 
    !
             errinlocb = 0
             errinlocb = max(errinlocb, abs(gamin(4 ) - gama(4 )) / gamin(4 ))
             errinlocb = max(errinlocb, abs(gamin(5 ) - gama(5 )) / gamin(5 ))
             errinlocb = max(errinlocb, abs(gamin(7 ) - gama(7 )) / gamin(7 ))
             errinlocb = max(errinlocb, abs(gamin(8 ) - gama(8 )) / gamin(8 ))
             errinlocb = max(errinlocb, abs(gamin(9 ) - gama(9 )) / gamin(9 ))
             errinlocb = max(errinlocb, abs(gamin(10) - gama(10)) / gamin(10))
             errinlocb = max(errinlocb, abs(gamin(13) - gama(13)) / gamin(13))
             calain    = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients 
             do ii = 1, 13
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do            
    !
    !  ## Solve system of equations, with new activity coefficients
             if (( .not. soln) .and. (count_rt == 1 .or. redo)) then
                lwnsq    = lwn*lwn
                gama10sq = gama(10)*gama(10)
                a3       = k1*lwnsq/gama10sq
                a4       = k2*gama10sq/(gama(5)*gama(5)) 
                a7       = k3*lwnsq                      
    !
    !  ## 1. Calculate dissociation quantities 
                k6   = a3*a4*(gnh3_i - omehi)
                psi3 = k6*ghno3_i - k4 - nh4no3*omehi
                psi3 = psi3 / (k6 + k5 + omehi)         
                psi3 = min(max(psi3, 0.0d0), ghno3_i)
    !
    !  ## 2. Calculate H+
                bb = (omehi - psi3)
    !  Use ISORROPIA II current scheme 
                denm = bb + sqrt(bb*bb + 4.0d0*a7)
                if (denm <= tiny) then
                   denm = (bb + abs(bb)) + 2.0d0*a7/abs(bb)
                end if
                h = 2.0d0*a7/denm
    !
    !  ## 3. Speciation 
                nh4_t = omehi + k5     
                so4_t = nh42s4
                no3_t = psi3 + nh4no3
                ghno3 = max(ghno3_i - psi3, tiny2)
                gnh3  = max(gnh3_i - omehi, tiny2)
    !       
    !  ## 4. Aerosol liquid water content
                m1  = (so4_t + hso4)       ! (NH4)2SO4
                m2  = (nh4_t - 2.0d0*m1)   ! free NH4   
                lwn = max(m1/c10 + max(min(m2, no3_t), 0.0d0)/c11, tiny)
             end if
          end do
    !
          y2 = nh4_t/h/max(gnh3, tiny)/a4 - 1.0d0  !Function value
    !
    !  ## Check for criteria to exit root tracking 
          condition = .false.
          if (sign(1.0d0,y1)*sign(1.0d0,y2) > 0.0d0 .and. ( .not. noroot) .and. abs(y2) > eps) then
             condition = .true.         
          elseif (sign(1.0d0,y1)*sign(1.0d0,y2) < 0.0d0 .and. abs(y2) > eps) then
    !  ## Interval had been found where sign change occurs; exit root tracking and proceed to ITP
             soln = .true. 
          else
    !  ## abs(y2) <= eps; solution is assumed; exit root tracking and proceed to minor system (no ITP)
             soln   = .true.
             noroot = .true.
             earlye = .true.
          end if 
    !
          if (rooteval == ndiv+1) then
             count_rt = count_rt + 1
    !
             if (redo .or. count_rt == 2) then
                yhi  = y2
                redo = .false.
             end if
    !
             if (abs(y2) < eps) then
                redo = .false.
                noroot = .true.
             else if (ylo < 0.0d0 .and. yhi < 0.0d0) then
                redo = .false.
                omehi = tiny
                noroot = .true.
             else if (ylo > 0.0d0 .and. yhi > 0.0d0) then
                if (count_rt == 2) then
                   pshi  = tiny
                   omebe = tiny
                   omehi = tiny - 0.1d0*(nh4no3 + nh42s4)   !psi4lo
                   psihi = omebe
                   psilo = omehi
                elseif (count_rt > 2) then
                   psihi = psilo
                   psilo = psihi - 0.1d0*(nh4no3 + nh42s4)
                   pshi  = psihi
                   omebe = psihi
                   omehi = psilo
                end if 
    ! 
                if (omehi < -1.0d0*(nh4no3 + nh42s4)) then
    !   ## No solution
                   omehi = tiny 
                   omebe = tiny
                   redo = .false.
                   noroot = .true.
                else 
    !   ## Include sulfate in initial calculation of liquid water content 
                   so4_t = nh42s4
                   hso4  = 0.0d0
                   nh4_t = nh4no3
                   no3_t = nh4no3
    !
                   m1  = (so4_t + hso4)      ! (NH4)2SO4
                   m2  = (nh4_t - 2.0d0*m1)  ! free NH4   
                   lwn = max(m1/c10 + max(min(m2, no3_t), 0.0d0)/c11, tiny) ! Water content 
    !
                   redo     = .true.
                   rooteval = 0 
                end if 
             else 
                redo = .false.
             end if 
          end if
       end do
    !
    !
    !  ### STAGE 2: modified bisection search (using ITP algorithm) ###
    !  ## Initialize static ITP variables 
       if (( .not. noroot)) then
          ya = y1
          yb = y2
          xa = omebe
          xb = omehi
          x3 = omehi
    !
          if (xa == xb) then
             noroot = .true.
             gx = tiny
          else
             gx = xb - xa
          end if 
    !
          gx2  = (xa+xb)*0.5d0
          u1   = 0.2d0/gx
          nh   = log10(abs(gx/(2.0d0*eps*gx2))) / log10of2    
          nmax = int(nh) + 2                              
       else
          x3 = omehi
       end if 
    !
    !  ## Start search
       y3_lastiter = 0.0d0
       x3_lastiter = 0.0d0
       y3_min = 1.0d20
       x3_min = 0.0d0
    !
       j = 0
       condition = .true.
    !
       if (earlye) then
          soln = .true.
       else
          soln = .false.
       end if
    !
       do while (j < maxit .and. condition)   ! Begin outer loop for ITP search
    !  ## Set dynamic ITP variables 
    !  1. Track x3 and y3 of the previous iteration
          if (j > 0) then
             y3_lastiter = y3
             x3_lastiter = x3
          end if 
    !    
    
    !  2. Track the minimum y3 that is found before ending
    !     This is for a post-convergence correction, in case of an oscillatory solution
          if (abs(0.0d0 - y3_min) > abs(0.0d0 - y3_lastiter) .and. j > 0) then
             y3_min    = y3_lastiter
             x3_min    = x3_lastiter
             ghno3_min = ghno3
             no3_min   = no3_t
             h_min     = h
             gnh3_min  = gnh3
             nh4_min   = nh4_t
             lwn_min   = lwn
             gama_min  = gama
          end if
    !
          if (( .not. noroot) .and. ( .not. soln)) then
             if (yb - ya == 0.0d0) then
                write(*,*), '######        ABORT       ######'
                write(*,*), 'Zero divide in ITP reset: CALCD3'
                write(*,*), 'SO4 in = ', so4
                write(*,*), 'NH4 in = ', nh4
                write(*,*), 'NO3 in = ', no3
                write(*,*), 'Temp in= ', t
                write(*,*), 'RH in  = ', aw
                return
             end if  
    !
             gx   = xb - xa
             xh   = 0.5d0*(xa + xb)
             rr   = max(gx2*eps*2.d0**(real(nmax - j)) - 0.5d0*gx, 0.0d0)
             delta= u1*(max(gx, 0.0d0))**2.0d0   
             xf   = max((yb*xa - ya*xb) / (yb - ya), 0.0d0)
    !
             sigma = sign(1.0d0, xh - xf)
             if (delta <= abs(xh - xf)) then   
                xt = xf + sigma*delta           
             else
                xt = xh
             end if
    !
             if (abs(xt - xh) <= rr) then
                x3 = xt
             else
                x3 = xh - sigma*rr
             end if 
          end if 
    !
          j = j + 1
    !
          if (( .not. soln)) then
             gmax = 0.1d0
             gmax = max(gmax, gama( 4))
             gmax = max(gmax, gama( 5))
             gmax = max(gmax, gama( 7))
             gmax = max(gmax, gama( 8))
             gmax = max(gmax, gama( 9))
             gmax = max(gmax, gama(10))
             gmax = max(gmax, gama(13))
          end if
    !
    !  ## Reinitialize activity coefficients if gmax > 100.0d0
          if (gmax > 100.0d0 .and. ( .not. soln)) then
             gama  = 0.1d0
             gamin = 1.0d10
             gamou = 1.0d10
             calou = .true.
             frst  = .true.
          end if
    !
    !  ## Solve system of equations
          frst    = .true.
          calain  = .true.
          if (( .not. soln)) then
             lwnsq    = lwn*lwn
             gama10sq = gama(10)*gama(10)
             a3       = k1*lwnsq/gama10sq 
             a4       = k2*gama10sq/(gama(5)*gama(5))  
             a7       = k3*lwnsq                      
    !
    !  ## 1. Calculate dissociation quantities
             k6   = a3*a4*(gnh3_i - x3) 
             psi3 = k6*ghno3_i  - k4 - nh4no3*x3
             psi3 = psi3 / (k6 + k5 + x3)         
             psi3 = min(max(psi3, 0.0d0), ghno3_i)
    !
    !  ## 2. Calculate H+
             bb = x3 - psi3
    !  Use ISORROPIA II current scheme 
             denm = bb + sqrt(bb*bb + 4.0d0*a7)
             if (denm <= tiny) then
                denm = (bb + abs(bb)) + 2.0d0*a7/abs(bb)
             end if
             h = 2.0d0*a7/denm
    !
    !  ## 3. Speciation 
             nh4_t = x3 + k5     
             so4_t = nh42s4
             no3_t = psi3 + nh4no3
             ghno3 = max(ghno3_i - psi3, tiny2)
             gnh3  = max(gnh3_i - x3, tiny2)
    !     
    !  ## 4. Aerosol liquid water content
             m1  = (so4_t + hso4)       ! (NH4)2SO4
             m2  = (nh4_t - 2.0d0*m1)   ! free NH4   
             lwn = max(m1/c10 + max(min(m2, no3_t), 0.0d0)/c11, tiny)
          end if
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)         
             k = k + 1
    !
    !  ## Reset gamin and gamou
             if (( .not. soln) .and. frst) then
                gamou( 4)  = gama( 4)
                gamou( 5)  = gama( 5)
                gamou( 7)  = gama( 7)
                gamou( 8)  = gama( 8)
                gamou( 9)  = gama( 9)
                gamou(10)  = gama(10)
                gamou(13)  = gama(13)
             end if 
    !
             if (( .not. soln)) then
                gamin( 4)  = gama( 4)
                gamin( 5)  = gama( 5)
                gamin( 7)  = gama( 7)
                gamin( 8)  = gama( 8)
                gamin( 9)  = gama( 9)
                gamin(10)  = gama(10)
                gamin(13)  = gama(13)
             end if 
    !
             call mach_hetp_calcact2b(h, nh4_t, so4_t, hso4, no3_t, lwn, gama, t, soln,    &
                                      frst, calain, calou)
    !
             if (frst) then
                errouloc = 0.0d0
                errouloc = max(errouloc, abs(gamou(4 ) - gama(4 )) / gamou(4 ))
                errouloc = max(errouloc, abs(gamou(5 ) - gama(5 )) / gamou(5 ))
                errouloc = max(errouloc, abs(gamou(7 ) - gama(7 )) / gamou(7 ))
                errouloc = max(errouloc, abs(gamou(8 ) - gama(8 )) / gamou(8 ))
                errouloc = max(errouloc, abs(gamou(9 ) - gama(9 )) / gamou(9 ))
                errouloc = max(errouloc, abs(gamou(10) - gama(10)) / gamou(10))
                errouloc = max(errouloc, abs(gamou(13) - gama(13)) / gamou(13))
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    !
             errinlocb = 0.0d0
             errinlocb = max(errinlocb, abs(gamin(4 ) - gama(4 )) / gamin(4 ))
             errinlocb = max(errinlocb, abs(gamin(5 ) - gama(5 )) / gamin(5 ))
             errinlocb = max(errinlocb, abs(gamin(7 ) - gama(7 )) / gamin(7 ))
             errinlocb = max(errinlocb, abs(gamin(8 ) - gama(8 )) / gamin(8 ))
             errinlocb = max(errinlocb, abs(gamin(9 ) - gama(9 )) / gamin(9 ))
             errinlocb = max(errinlocb, abs(gamin(10) - gama(10)) / gamin(10))
             errinlocb = max(errinlocb, abs(gamin(13) - gama(13)) / gamin(13))
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients 
             do ii = 1, 13
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do        
    !
    !  ## Solve system of equations, with new activity coefficients
             if (( .not. soln)) then
                lwnsq    = lwn*lwn
                gama10sq = gama(10)*gama(10)
                a3       = k1*lwnsq/gama10sq 
                a4       = k2*gama10sq/(gama(5)*gama(5))  
                a7       = k3*lwnsq                      
    !
    !  ## 1. Calculate dissociation quantities 
                k6   = a3*a4*(gnh3_i - x3)
                psi3 = k6*ghno3_i  - k4 - nh4no3*x3
                psi3 = psi3 / (k6 + k5 + x3)         
                psi3 = min(max(psi3, 0.0d0), ghno3_i)
    !
    !  ## 2. Calculate H+
                bb = x3 - psi3
    !  Use ISORROPIA II current scheme 
                denm = bb + sqrt(bb*bb + 4.0d0*a7)
                if (denm <= tiny) then
                   denm = (bb + abs(bb)) + 2.0d0*a7/abs(bb)
                end if
                h = 2.0d0*a7/denm
    !
    !  ## 3. Speciation 
                nh4_t = x3 + k5     
                so4_t = nh42s4
                no3_t = psi3 + nh4no3
                ghno3 = max(ghno3_i - psi3, tiny2)
                gnh3  = max(gnh3_i - x3, tiny2)
    !       
    !  ## 4. Aerosol liquid Water content
                m1  = (so4_t + hso4)       ! (NH4)2SO4
                m2  = (nh4_t - 2.0d0*m1)   ! free NH4   
                lwn = max(m1/c10 + max(min(m2, no3_t), 0.0d0)/c11, tiny)
             end if
          end do
    !     
          y3 = nh4_t/h/max(gnh3, tiny)/a4 - 1.0d0
    !
          condition = .false.
          if (noroot) then
    !  ## If no root on interval, then do not perform ITP
             xa = x3
             xb = x3
          else if (y3 > 0.0d0 .and. ( .not. soln)) then
             if (ya < 0.0d0) then
                yb = y3
                xb = x3
             elseif (ya > 0.0d0) then
                ya = y3
                xa = x3
             end if 
          else if (y3 < 0.0d0 .and. ( .not. soln)) then
             if (ya < 0.0d0) then
                ya = y3
                xa = x3
             elseif (ya > 0.0d0) then
                yb = y3
                xb = x3
             end if 
          else if (( .not. soln)) then
             xa = x3
             xb = x3
          end if
    !
    !  ## Check for convergence criteria to exit ITP:
          if (xb - xa > abs(xa*eps) .and. ( .not. noroot)) then
             condition = .true.
             soln   = .false.
          else
             soln   = .true.
          end if
    !
    ! ## Exit ITP if the function being bisected evaluates to a value <= eps; solution is assumed
          if (abs(y3) <= eps .and. ( .not. noroot)) then
             soln = .true.
             condition = .false.
          end if
    !
    ! ## Post-convergence correction: 
    ! ## If the calculated y-value (y3) after ITP is further from zero than an earlier iteration 
    ! ## then reset to the x-value/concentrations/activity coefficients that were found to minimize
    ! ## the objective function (i.e., y3); in this case, this is chosen as the solution
          if (( .not. condition) .and. ( .not. noroot) .and. abs(y3) > 0.1d0) then     
             if (abs(0.0d0 - y3_min) < abs(0.0d0 - y3) .and. abs(y3_min - y3) > 1.0d-1) then
    !            write(*,*), 'Warning: oscillatory behavior; possibility of no valid solution!'
                x3 = x3_min
                nh4_t = nh4_min
                no3_t = no3_min
                h     = h_min
                lwn   = lwn_min
                ghno3 = ghno3_min
                gnh3  = gnh3_min
    ! 
    !  ## Reset activity coefficients 
                gama = gama_min
                a4   = k2*(gama(10)/gama(5))*(gama(10)/gama(5))
                y3   = nh4_t/h/max(gnh3, tiny)/a4 - 1.0d0
             end if
          end if 
       end do
    ! 
    !
    !  ### MINOR SYSTEM: HSO4-/SO42-/H+ ###
      hso4 = 0.0d0
      if (h > tiny) then
         ak1 = khso4*lwn/gama(7)*(gama(8)/gama(7))*(gama(8)/gama(7))
         bb  =-(h + so4_t + ak1)
         cc  = h*so4_t 
         dd  = bb*bb - 4.d0*cc
    ! 
         if (dd >= 0.0d0) then
            if (bb /= 0.d0) then
               dd = cc/(bb*bb)
               v  = 4.d0*dd
            else
               v = 1.d3
            end if
    !
            if (abs(v) <= smrt .and. bb /= 0.d0) then
               hh = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
            else
               hh = 0.5d0*(-bb - sqrt(bb*bb - 4.d0*cc))
            end if 
    !
            hh    = max(tiny, min(hh, min(h, so4_t)))     ! To avoid negative H+   (i.e., if hh > h or hh > so4_t)
            h     = max(h - hh, 0.0d0)
            so4_t = max(so4_t - hh, 0.0d0)
            hso4  = hh
         end if 
      end if  
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       h_i     = h
       hno3g_i = ghno3
       nh3g_i  = gnh3
       lwn_i   = lwn
    !
       return
    end subroutine mach_hetp_calcd3
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: E4; Sulfate rich, no acid
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calce4(so4_i, nh4_i, hno3g_i, hso4_i, h_i, no3_i,              &
                                lwn_i, rh, temp, k0, p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: hno3g_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables
       real(kind=8)      :: so4, nh4, hso4, ghno3, no3, h, lwn, t, aw
       real(kind=8)      :: khso4, kh2o, khno3, bb, cc, dd, hh, v, ak1, errin, tt0
       real(kind=8)      :: so4_t, nh4_t, no3_t, m4, m9, m13, tt1, tt2, c1, c2
       real(kind=8)      :: gnh3, cl, cl_t, ghcl, caso4, c3, c4, c5
       integer(kind=4)   :: j, irh
       logical(kind=4)   :: gg
       real(kind=8), dimension(13) :: gama, gamin
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       ghno3 = 0.0d0
       h     = 0.0d0
       lwn   = 0.0d0
       so4_t = 0.0d0
       nh4_t = nh4
       no3_t = 0.0d0
       gnh3  = 0.0d0
       cl    = 0.0d0
       cl_t  = 0.0d0
       ghcl  = 0.0d0
       caso4 = 0.0d0
       gama  = 0.1d0
       gamin = 1.0d10
    !
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !  
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !   
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3 = k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Calculate dry composition 
       gg = 2.0*so4 - nh4 <= nh4 - so4
    !
       if (gg) then   
    !  NH4HSO4 >= (NH4)2SO4
          m4  = 2.0d0*nh4 - 3.0d0*so4
          m9  = 0.0d0
          m13 = 2.0d0*so4 - nh4
       else 
    !  NH4HSO4 < (NH4)2SO4
          m4  = 0.0d0
          m9  = 3.0d0*so4 - 2.0d0*nh4
          m13 = nh4 - so4
       end if
    !
    !  ## Aerosol liquid water content (initial)
    !  Save ZSR arrays as variables to limit indirect addressing
       c3 = awlc(irh)
       c4 = awab(irh)
       c5 = awas(irh)
       lwn  = m13/c3 + m9/c4 + m4/c5
    !
    !  Constant variables 
       c1   = khno3*r*t
       c2   = khso4*(lwn/0.1d0)
    !
    !
    !  ### MAJOR SYSTEM H+/HSO4-/SO42- ###
    !  ## Setup initial conditions
       bb = so4 + c2 - nh4  
       cc = -c2*so4
    !
       if (bb /= 0.d0) then
          dd = cc / (bb*bb)
          v  = 4.d0*dd
       else
          v  = 1.d+03
       end if       
    !
       if (abs(v) <= smrt .and. bb /= 0.d0) then
          hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
       else
          hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))                            !positive root
       end if 
    !
    !  ## 2. Speciation  
       so4_t = max(tiny, min(hh, so4))                  
       hso4  = max(tiny, min(so4 - so4_t, so4)) 
       h     = max(tiny, min(c2*hso4/so4_t, so4))
    !  
    !  ## 3. Aerosol liquid water content       
       gg = so4_t - h < hso4 + h
       if (gg) then
          m4   = 0.0d0
          m13  = so4_t - h
          m9   = max((hso4 + h) - (so4_t - h), 0.0d0)
       else
          m13  = hso4 + h
          m9   = 0.0d0
          m4   = max((so4_t - h) - (hso4 + h), 0.0d0)
       end if 
       lwn  = max(m13/c3 + m9/c4 + m4/c5, tiny)
    !  
    !
    !  ## Iterative search for solution with convergence of activity coefficients
       j     = 0
       errin = 1.0d0
       do while (j < nsweep-1 .and. errin >= epsact)
          j = j + 1
    
    !  ## Reset gamin
          gamin(4)  = gama(4)
          gamin(5)  = gama(5)
          gamin(7)  = gama(7)
          gamin(8)  = gama(8)
          gamin(9)  = gama(9)
          gamin(10) = gama(10)
          gamin(13) = gama(13)
    !     
          call mach_hetp_calcact2(h, nh4_t, so4_t, hso4, no3_t, lwn, gama, t)
    !
    !  ## Test convergence criterion: max change in the activity coefficients
          errin = 0.0d0
          errin = max(errin, abs((gamin(7) - gama(7))  / gamin(7)))
          errin = max(errin, abs((gamin(8) - gama(8))  / gamin(8)))
          errin = max(errin, abs((gamin(9) - gama(9))  / gamin(9)))
          errin = max(errin, abs((gamin(10)- gama(10)) / gamin(10)))
          errin = max(errin, abs((gamin(4) - gama(4))  / gamin(4)))
          errin = max(errin, abs((gamin(5) - gama(5))  / gamin(5)))
          errin = max(errin, abs((gamin(13)- gama(13)) / gamin(13)))
    
    !  ## 1. Solve system of equations 
          ak1 = khso4*((gama(8)/gama(7))*(gama(8)/gama(7)))*(lwn/gama(7))
          bb  = so4 + ak1 - nh4_t
          cc  = -ak1*so4 
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v = 4.d0*dd
          else
             v = 1.d+03
          end if       
    !
          if (abs(v) <= smrt .and. bb /= 0.d0) then
             hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
          else
             hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))   !positive root
          end if 
    !
    !  ## 2. Speciation
          so4_t = max(tiny, min(hh, so4))
          hso4  = max(tiny, min(so4 - so4_t, so4))
          h     = max(tiny, min(ak1*hso4/so4_t, so4))
    !
    !  ## 3. Aerosol liquid water content
          gg = so4_t - h < hso4 + h
          if (gg) then  
             m4  = 0.0d0 
             m13 = so4_t - h
             m9  = max((hso4 + h) - (so4_t - h), 0.0d0)  ! NH4HSO4
          else 
             m9  = 0.0d0
             m13 = hso4 + h
             m4  = max((so4_t - h) - (hso4 + h), 0.0d0)  ! (NH4)2SO4
          end if 
          lwn = max(m13/c3 + m9/c4 + m4/c5, tiny)
       end do
    !
    !
    ! ### MINOR SYSTEM: NO3-/HNO3/H+ ###
    ! Nitric acid in the liquid phase assumed a minor species; nitric acid is dissolved
    ! from the [HNO3(g)] -> [H+] + [NO3-] equilibrium using the [H+] from the major system
      hh = 0.0d0
      if (lwn > tiny) then
          bb = c1*(lwn/gama(10))*(lwn/gama(10)) + h    !!bb always > 0
          cc = c1*(lwn/gama(10))*(lwn/gama(10))*no3
    
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v = 4.d0*dd
          else
             v = 1.d+03
          end if       
    !
          if (abs(v) <= smrt .and. bb /= 0.d0) then
             hh = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
          else
             hh = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
          end if 
       end if 
    !
       ghno3 = max(no3 - hh, 0.0d0)  
       no3_t = hh
       h     = h + hh  
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       hso4_i  = hso4
       no3_i   = no3_t
       h_i     = h
       hno3g_i = ghno3
       lwn_i   = lwn
    !
       return
    end subroutine mach_hetp_calce4
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: F2; Sulfate rich, free acid
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcf2(so4_i, nh4_i, hno3g_i, hso4_i, h_i, no3_i,              &
                                lwn_i, rh, temp, k0, p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: hno3g_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables
       real(kind=8)     :: so4, nh4, no3, hso4, ghno3, h, lwn, frh2so4, t, aw
       real(kind=8)     :: khso4, kh2o, khno3, bb, cc, dd, hh, v, errin, tt0
       real(kind=8)     :: so4_t, nh4_t, no3_t, tt1, tt2, c1, c2, c3, c4
       real(kind=8)     :: gnh3, cl, cl_t, ghcl, caso4
       integer(kind=4)  :: j, irh
       real(kind=8), dimension(13) :: gama, gamin
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       ghno3 = 0.0d0
       h     = 0.0d0
       lwn   = 1.0d-20
       so4_t = 0.0d0
       nh4_t = 0.0d0
       no3_t = 0.0d0
       frh2so4 = 0.0d0
       gnh3  = 0.0d0
       cl    = 0.0d0
       cl_t  = 0.0d0
       ghcl  = 0.0d0
       caso4 = 0.0d0
       gama  = 0.1d0
       gamin = 1.0d10
    !
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !   
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3 = k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh  = max(min(int(aw*100+0.5), 100), 1)
    
    !  ## Aerosol liquid water content 
       c3   = so4 - nh4
       lwn  = max(max(c3, 0.0d0)/awsa(irh) + nh4/awab(irh), tiny)
    !
    !  ## Constant values
       c1   = khno3*r*t
       c2   = khso4*1.0d-19
       c4   = khso4*lwn
    !
    !
    !  ### MAJOR SYSTEM H+/HSO4-/SO42- ###
    !  ## Setup initial conditions
       bb = c3 + c2 
       cc = -c2*so4
    !
       if (bb /= 0.d0) then
          dd = cc / (bb*bb)
          v  = 4.d0*dd
       else
          v  = 1.d+03
       end if       
    !
       if (abs(v) <= smrt .and. bb /= 0.d0) then
          hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
       else
          hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))                            !positive root
       end if      
    !
    !  ## 2. Speciation  
       nh4_t   = nh4
       so4_t   = hh                
       hso4    = max(so4-hh, tiny) 
       h       = c3 + hh
       frh2so4 = max(so4_t + hso4 - nh4_t, 0.0d0)
    !  
    !
    !  ## Iterative search for solution with convergence of activity coefficients
       j     = 0
       errin = 1.0d0
       do while (j < nsweep-1 .and. errin >= epsact)
          j = j + 1
    !
    !  ## Reset gamin
          gamin(4)  = gama(4)
          gamin(5)  = gama(5)
          gamin(7)  = gama(7)
          gamin(8)  = gama(8)
          gamin(9)  = gama(9)
          gamin(10) = gama(10)
          gamin(13) = gama(13)
    !
          call mach_hetp_calcact2(h, nh4_t, so4_t, hso4, no3_t, lwn, gama, t)
    !
    !  ## Test convergence criterion: max change in the activity coefficients
          errin = 0.0d0
          errin = max(errin, abs((gamin(7) - gama(7))  / gamin(7)))
          errin = max(errin, abs((gamin(8) - gama(8))  / gamin(8)))
          errin = max(errin, abs((gamin(9) - gama(9))  / gamin(9)))
          errin = max(errin, abs((gamin(10)- gama(10)) / gamin(10)))
          errin = max(errin, abs((gamin(4) - gama(4))  / gamin(4)))
          errin = max(errin, abs((gamin(5) - gama(5))  / gamin(5)))
          errin = max(errin, abs((gamin(13)- gama(13)) / gamin(13)))
    !
    !  ## 1. Solve system of equations 
          bb = c3 + c4/gama(7)*(gama(8)/gama(7))**2.0 
          cc = -(c4/gama(7)*(gama(8)/gama(7))**2.0)*so4
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d+03
          end if       
    !
          if (abs(v) <= smrt .and. bb /= 0.d0) then
             hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
          else
             hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))   !positive root
          end if 
    !
    !  ## 2. Speciation
          so4_t   = hh                
          hso4    = max(so4 - hh, tiny) 
          h       = c3 + hh
          frh2so4 = max(so4_t + hso4 - nh4_t, 0.0d0)
       end do
    !
    !
    ! ### MINOR SYSTEM: NO3-/HNO3/H+ ###
    ! Nitric acid in the liquid phase assumed a minor species; nitric acid is dissolved
    ! from the [HNO3(g)] -> [H+] + [NO3-] equilibrium using the [H+] from the major system
      hh = 0.0d0
      if (lwn > tiny) then
          bb = c1*(lwn/gama(10))*(lwn/gama(10)) + h    !!bb always > 0
          cc = c1*(lwn/gama(10))*(lwn/gama(10))*no3
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v = 4.d0*dd
          else
             v = 1.d+03
          end if       
    !
          if (abs(v) <= smrt .and. bb /= 0.d0) then
             hh = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
          else
             hh = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
          end if 
       end if 
    !
       ghno3 = max(no3 - hh, 0.0d0)  
       no3_t = hh
       h     = h + hh  
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       h_i     = h
       hno3g_i = ghno3
       lwn_i   = lwn
    !
       return
    end subroutine mach_hetp_calcf2
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: G5; Sulfate poor; sodium poor
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcg5(so4_i, nh4_i, nh3g_i, hno3g_i, hclg_i, hso4_i, na_i,    &
                                cl_i, no3_i, h_i, lwn_i, rh, temp, k0, p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: na_i    
       real(kind=8),    intent(inout) :: cl_i    
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: hclg_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables
       real(kind=8)      :: so4, nh4, hso4, gnh3, h, lwn, no3, cl, na, ghno3, ghcl, caso4, t, aw
       real(kind=8)      :: khso4, knh3, kh2o, khno3, khcl, bb, cc, dd, hh, v, ak1, frnh4, tt0
       real(kind=8)      :: errouloc, errinlocb, errin, ohi, smin, scon, m4, m5, m6
       real(kind=8)      :: a4, a5, psi4, psi5, nh, sigma, xt, xf, xh, delta, rr, gx, gx2
       real(kind=8)      :: u1, ya, yb, xa, xb, c1, c2, c3, c3a, c4, c5, k1, k2, k3, k4, loccon
       real(kind=8)      :: omehi, omebe, y1, y2, y3, x3, dx, a6, gmax, zsr2
       real(kind=8)      :: so4_t, nh4_t, no3_t, na_t, cl_t, tt1, tt2, zsr3, zsr4
       real(kind=8)      :: y3_min, x3_min, y3_lastiter, x3_lastiter, lwnsq, gama10sq
       real(kind=8)      :: no3_min, nh4_min, cl_min, h_min, lwn_min, ghcl_min, gnh3_min, ghno3_min
       integer(kind=4)   :: j, k, ii, rooteval, irh, nmax
       logical(kind=4)   :: condition, noroot, soln, frst, calain, calou, earlyexit
       real(kind=8), dimension(13) :: gama, gamin, gamou, gama_min
    ! 
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       na    = na_i
       cl    = cl_i
       aw    = rh
       t     = temp
       caso4 = 0.0d0
       hso4  = 0.0d0
       gnh3  = 0.0d0
       ghno3 = 0.0d0
       ghcl  = 0.0d0
       h     = 0.0d0
       lwn   = 0.0d0
       so4_t = 0.0d0
       nh4_t = 0.0d0
       no3_t = 0.0d0
       na_t  = 0.0d0
       cl_t  = 0.0d0
       noroot=.false.
       soln  =.false.
       gmax  = 0.0d0
       calou =.true.
       gama  = 0.1d0
       gamin = 1.0d10
       gamou = 0.1d0
       earlyexit = .false.
    !
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    ! 
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !   
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3 = k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. HCl(g) <==> H+(aq) + Cl-(aq)                               (xk3)
       khcl = k0(11) * exp(p1(11)*tt1 + p2(11)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Create variables to hold ZSR variables (to limit indirect addressing)
       zsr2 = awas(irh)
       zsr3 = awan(irh)
       zsr4 = awac(irh)
    !
    !  ## Constant values
       c1 = r*t
       c2 = 0.5d0*na                !na2so4
       c5 = c2/awss(irh)
       c3 = max(so4 - c2, 0.0d0)    !frso4
       c3a= 2.0d0*c3                
       c4 = max(nh4 - c3a, 0.0d0)   !frnh4
       k1 = knh3/kh2o*c1
       k2 = khno3*c1
       k3 = khcl*c1
       k4 = kh2o*aw
    !
    !  ## Initial speciation
       na_t  = na
       so4_t = c3 + c2
    !
    !  ## Initial aerosol liquid water content 
       lwn = c3/zsr2 + c5
    !
    !
    !  ### STAGE 1: Root tracking ###
    !  ## Find a subinterval [xa,xb] on the larger interval [I1,I2] where a sign change occurs
       rooteval  = 0
       condition = .true.
       do while (rooteval < 2 .or. (condition .and. rooteval < ndiv + 1))   ! Begin outer loop for root tracking
          rooteval = rooteval + 1
    !
    !  ## Set high limit for root tracking (i.e. lower bound)
          if (rooteval == 1) then
            omehi = tiny   
            y1    = 0.0d0
          end if
    !      
    !  ## Begin search on subinterval 
          if (rooteval == 2) then
             soln = .false.
             y1   = y2
    !
             dx = (cl-tiny-tiny)/float(ndiv) 
    !
             omebe = omehi            ! Lower bound of subinterval (xa)
             omehi = omehi + dx       ! Upper bound of subinterval (xb)
          end if 
    !  
    !  ## Continue search 
          if (rooteval > 2) then
             if (loccon < 0.0d0) then
    !  ## 1. Root has been found on the subinterval; save x values for ITP search
                y1    = y1
                omebe = omebe
                omehi = omehi
             else
    !  ## 2. No root has been found, continue searching in the next subinterval
                y1    = y2
                omebe = omehi
                omehi = omehi + dx
             end if 
          end if 
    !      
    !  ## Solve the system of equations 
          frst = .true.
          calain = .true.
          if (( .not. soln)) then
             lwnsq = lwn*lwn
             gama10sq = gama(10)*gama(10)
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))   
             a5 = k2*lwnsq/gama10sq               
             a6 = k3*lwnsq/(gama(11)*gama(11))    
    ! 
    !  ## 1. Calculate dissociation quantities 
             if (no3 >= tiny) then
                psi5 = omehi*no3/(a6/a5*(cl-omehi) + omehi)
             else 
                psi5 = tiny
             end if 
    !
    !  ## 2. Account for NH3 evaporation
             if (so4 > tiny) then
                bb = -(c4 + omehi + psi5 + 1.0d0/a4)
                cc = c4*(psi5 + omehi) - c3a/a4
    !               
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                else
                   psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if 
             else
                psi4 = tiny
             end if 
    !
    !  ## 3. Speciation
             nh4_t = c3a + psi4 
             cl_t  = omehi        
             no3_t = psi5   
             gnh3  = max(c4  - psi4,  tiny2)
             ghno3 = max(no3 - psi5,  tiny2)
             ghcl  = max(cl  - omehi, tiny2)
    !        
    !  ## 4. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t
             scon = k4*lwnsq  
    !
             if (smin > tiny) then                        ! H+ in excess
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                end if 
    !  
                h    = max(hh,sqrt(scon))
    !           ohi  = scon/h
             else                                        ! OH- in excess
                bb = smin
                cc =-scon
    ! 
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                end if 
    !
                ohi  = max(hh,sqrt(scon))
                h = scon/ohi 
             end if
    !
    !  ## 5. Aerosol liquid water content
             m4     = max(so4_t - c2, 0.0d0)		        ! (NH4)2SO4
             frnh4  = max(nh4_t - 2.0d0*m4, 0.0d0)
             m5     = min(no3_t, frnh4)				! NH4NO3
             frnh4  = max(frnh4 - m5, 0.0d0)
             m6     = min(cl_t, frnh4)				! NH4Cl
             lwn    = max(c5 + m4/zsr2 + m5/zsr3 + m6/zsr4, tiny)
          end if 
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)        
             k = k + 1
    !        
    !  ## Reset gamou
             if (( .not. soln) .and. frst) then
                gamou(1)  = gama(1)
                gamou(2)  = gama(2)
                gamou(3)  = gama(3)
                gamou(4)  = gama(4)
                gamou(5)  = gama(5)
                gamou(6)  = gama(6)
                gamou(7)  = gama(7)
                gamou(8)  = gama(8)
                gamou(9)  = gama(9)
                gamou(10) = gama(10)
                gamou(11) = gama(11)
                gamou(12) = gama(12)
                gamou(13) = gama(13)
             end if 
    !
    !  ## Reset gamin
             if (( .not. soln)) then
                gamin(1)  = gama(1)
                gamin(2)  = gama(2)
                gamin(3)  = gama(3)
                gamin(4)  = gama(4)
                gamin(5)  = gama(5)
                gamin(6)  = gama(6)
                gamin(7)  = gama(7)
                gamin(8)  = gama(8)
                gamin(9)  = gama(9)
                gamin(10) = gama(10)
                gamin(11) = gama(11)
                gamin(12) = gama(12)
                gamin(13) = gama(13)
             end if 
    !
             call mach_hetp_calcact3b(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                     lwn, gama, t, soln, frst, calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 13
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    !
             errinlocb = 0
             do ii = 1, 13
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact  
    !  
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients
             do ii = 1, 13
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do
    !
    !  ## Solve system of equations, using new activity coefficients 
            if (( .not. soln)) then
                lwnsq = lwn*lwn
                gama10sq = gama(10)*gama(10)
    !
                a4 = k1*gama10sq/(gama(5)*gama(5))   
                a5 = k2*lwnsq/gama10sq               
                a6 = k3*lwnsq/(gama(11)*gama(11))     
    !
    !  ## 1. Calculate dissociation quantities 
                if (no3 >= tiny) then
                   psi5 = omehi*no3/(a6/a5*(cl - omehi) + omehi)
                else 
                   psi5 = tiny
                end if 
    !
    !  ## 2. Account for NH3 evaporation
                if (so4 > tiny) then
                   bb = -(c4 + omehi + psi5 + 1.0d0/a4)
                   cc = c4*(psi5 + omehi) - c3a/a4
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                   else
                      psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if 
                else
                   psi4 = tiny
                end if 
    !
    !  ## 3. Speciation
                nh4_t = c3a + psi4 
                cl_t  = omehi        
                no3_t = psi5
                gnh3  = max(c4  - psi4,  tiny2)
                ghno3 = max(no3 - psi5,  tiny2)
                ghcl  = max(cl  - omehi, tiny2)
    !
    !  ## 4. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t
                scon = k4*lwnsq   
    !
                if (smin > tiny) then                        ! H+ in excess        
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                   end if 
    !
                   h = max(hh,sqrt(scon))
    !               ohi  = scon/h
                else                                        ! OH- in excess
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                   end if 
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if
    !       
    !  ## 5. Aerosol liquid water content 
                m4     = max(so4_t - c2, 0.0d0)		                 ! (NH4)2SO4
                frnh4  = max(nh4_t - 2.0d0*m4, 0.0d0)
                m5     = min(no3_t, frnh4)				         ! NH4NO3
                frnh4  = max(frnh4 - m5, 0.0d0)
                m6     = min(cl_t, frnh4)				         ! NH4Cl
                lwn    = max(c5 + m4/zsr2 + m5/zsr3 + m6/zsr4, tiny)
             end if 
          end do
    !     
          y2 = h*cl_t/ghcl/a6 - 1.0d0  ! Function value
    !
    !  ## Check for criteria to exit root tracking before ndiv iterations of root tracking
          condition = .false.
          loccon = sign(1.0d0,y1)*sign(1.0d0,y2)
          if (loccon > 0.0d0 .and. ( .not. noroot) .and. abs(y2) > eps) then
             condition = .true.         
          elseif (loccon < 0.0d0 .and. abs(y2) > eps) then
    !  ## Interval had been found where sign change occurs; exit root tracking and proceed to ITP
             soln = .true. 
          else
    !  ## abs(y2) <= eps; solution is assumed; exit root tracking and proceed to minor system (no ITP)
             if (rooteval >= 2) then
                ! ## Root has been found on the end of a subinterval
                soln   = .true.
                noroot = .true.
                earlyexit = .true.
             else
                ! ## Root is found at the start of larger interval, but still proceed
                condition = .true.
             end if
          end if 
    !
    !  ## Too little Cl; reset x-value to tiny and exit root tracking 
          if (cl <= tiny) then
             condition = .false.
             rooteval  = 2       ! Increment rooteval by 1 to force exit from root tracking
             noroot    = .true.
             omehi     = tiny    ! Reset x-value to tiny, and solve system with this value
             omebe     = omehi
          end if 
    !
    !  ### AFTER iterating through ALL ndiv subdivided intervals
          if (rooteval == ndiv + 1) then
             if (loccon > 0.0d0 .and. abs(y2) > eps) then
    !  ## (1) No solution
                noroot = .true.
                omehi  = tiny    ! Reset to tiny 
                omebe  = omehi
    !            write(*,*), 'Warning in CALCG5: no solution found'
             else if (loccon > 0.0d0 .and. abs(y2) <= eps) then
    !  ## (2) Solution is assumed and ITP is not required
                noroot = .true.
                soln   = .true.
             end if 
          end if  
       end do  !End outer loop for root tracking
    !
    !
    !
    !  ### STAGE 2: modified bisection search (using ITP algorithm) ###
    !  ## Initialize static ITP variables 
       if (( .not. noroot)) then
          ya = y1
          yb = y2
          xa = omebe
          xb = omehi
          x3 = omehi
    !
          if (xa == xb) then
             noroot = .true.
             gx = tiny
          else
             gx = xb - xa
          end if 
    !
          gx2  = (xa+xb)*0.5d0
          u1   = 0.2d0/gx
          nh   = log10(abs(gx/(2.0d0*eps*gx2))) / log10of2   
          nmax = int(nh) + 2                               
       else
          x3 = omehi
       end if 
    ! 
    !  ## Start search
    !  1. Store value of previous y3 for later use
       y3_lastiter = 0.0d0
       x3_lastiter = 0.0d0
       y3_min = 1.0d20
       x3_min = 0.0d0
    !
       j = 0
       condition = .true.
    !
       if (( .not. earlyexit)) then
          soln = .false.
       end if 
    !
       do while (j < maxit .and. condition)   
    !  ## Set dynamic ITP variables 
    !  1. Track x3 and y3 of the previous iteration
          if (j > 0) then
             y3_lastiter = y3
             x3_lastiter = x3
          end if 
    !
    !  2. Track the minimum y3 that is found before ending
          if (abs(0.0d0 - y3_min) > abs(0.0d0 - y3_lastiter) .and. j > 0) then
             y3_min    = y3_lastiter
             x3_min    = x3_lastiter
             ghno3_min = ghno3
             no3_min   = no3_t
             h_min     = h
             ghcl_min  = ghcl
             cl_min    = cl_t
             gnh3_min  = gnh3
             nh4_min   = nh4_t
             lwn_min   = lwn
             gama_min  = gama
          end if
    !
          if (( .not. noroot) .and. ( .not. soln)) then
             if (yb - ya == 0.0d0) then
                write(*,*), '######        ABORT       ######'
                write(*,*), 'Zero divide in ITP reset: CALCG5'
                write(*,*), 'SO4 in = ', so4
                write(*,*), 'NH4 in = ', nh4
                write(*,*), 'NO3 in = ', no3
                write(*,*), 'Na in  = ', na
                write(*,*), 'Cl in  = ', cl
                write(*,*), 'Temp in= ', t
                write(*,*), 'RH in  = ', aw
                return
             end if  
    !
             gx   = xb - xa
             xh   = 0.5d0*(xa + xb)
             rr   = max(gx2*eps*2.d0**(real(nmax - j)) - 0.5d0*gx, 0.0d0)
             delta= u1*(max(gx, 0.0d0))**2.0d0   
             xf   = max((yb*xa - ya*xb) / (yb - ya), 0.0d0)
    !
             sigma = sign(1.0d0, xh - xf)
             if (delta <= abs(xh - xf)) then
                xt = xf + sigma*delta           
             else
                xt = xh
             end if
    !
             if (abs(xt - xh) <= rr) then
                x3 = xt
             else
                x3 = xh - sigma*rr
             end if
          end if 
    !
          j = j + 1
    !
          if (( .not. soln)) then
             gmax = 0.1d0
             gmax = max(gmax, gama(1))
             gmax = max(gmax, gama(2))
             gmax = max(gmax, gama(3))
             gmax = max(gmax, gama(4))
             gmax = max(gmax, gama(5))
             gmax = max(gmax, gama(6))
             gmax = max(gmax, gama(7))
             gmax = max(gmax, gama(8))
             gmax = max(gmax, gama(9))
             gmax = max(gmax, gama(10))
             gmax = max(gmax, gama(11))
             gmax = max(gmax, gama(12))
             gmax = max(gmax, gama(13))
          end if
    !
    !  ## Reinitialize activity coefficients if gmax > 100.0d0
          if (gmax > 100.0d0 .and. ( .not. soln)) then
             gama  = 0.1d0
             gamin = 1.0d10
             gamou = 1.0d10
             calou = .true.
             frst  = .true.
          end if
    !
    !  ## Solve system of equations 
          frst    = .true.
          calain  = .true.
          if (( .not. soln)) then
             lwnsq    = lwn*lwn
             gama10sq = gama(10)*gama(10)
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))  
             a5 = k2*lwnsq/gama10sq              
             a6 = k3*lwnsq/(gama(11)*gama(11))   
    ! 
    !  ## 1. Calculate dissociation quantities 
             if (no3 >= tiny) then
                psi5 = x3*no3/(a6/a5*(cl-x3) + x3)
             else 
                psi5 = tiny
             end if 
    !
    !  ## 2. Account for NH3 evaporation
             if (so4 > tiny) then
                bb = -(c4 + x3 + psi5 + 1.0d0/a4)
                cc = c4*(psi5 + x3) - c3a/a4
    !               
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                else
                   psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if 
             else
                psi4 = tiny
             end if 
    !
    !  ## 3. Speciation
             nh4_t = c3a + psi4 
             cl_t  = x3        
             no3_t = psi5   
             gnh3  = max(c4  - psi4, tiny2)
             ghno3 = max(no3 - psi5, tiny2)
             ghcl  = max(cl  - x3,   tiny2)
    !
    !  ## 4. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t
             scon = k4*lwnsq  
    !
             if (smin > tiny) then                        ! H+ in excess
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                end if 
    !
                h   = max(hh,sqrt(scon))
    !            ohi = scon/h
             else                                        ! OH- in excess    
                bb = smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                end if 
    !
                ohi = max(hh,sqrt(scon))
                h   = scon/ohi
             end if 
    !        
    !  ## 5. Aerosol liquid water content
             m4     = max(so4_t - c2, 0.0d0)		           ! (NH4)2SO4
             frnh4  = max(nh4_t - 2.0d0*m4, 0.0d0)   
             m5     = min(no3_t, frnh4)				   ! NH4NO3
             frnh4  = max(frnh4 - m5, 0.0d0)
             m6     = min(cl_t, frnh4)				   ! NH4Cl
             lwn    = max(c5 + m4/zsr2 + m5/zsr3 + m6/zsr4, tiny)
          end if 
    !
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)         
             k = k + 1
    !  ## Reset gamin and gamou
             if (( .not. soln) .and. frst) then
                gamou(1)  = gama(1)
                gamou(2)  = gama(2)
                gamou(3)  = gama(3)
                gamou(4)  = gama(4)
                gamou(5)  = gama(5)
                gamou(6)  = gama(6)
                gamou(7)  = gama(7)
                gamou(8)  = gama(8)
                gamou(9)  = gama(9)
                gamou(10) = gama(10)
                gamou(11) = gama(11)
                gamou(12) = gama(12)
                gamou(13) = gama(13)
             end if 
    !
             if (( .not. soln)) then
                gamin(1)  = gama(1)
                gamin(2)  = gama(2)
                gamin(3)  = gama(3)
                gamin(4)  = gama(4)
                gamin(5)  = gama(5)
                gamin(6)  = gama(6)
                gamin(7)  = gama(7)
                gamin(8)  = gama(8)
                gamin(9)  = gama(9)
                gamin(10) = gama(10)
                gamin(11) = gama(11)
                gamin(12) = gama(12)
                gamin(13) = gama(13)
             end if 
    !
             call mach_hetp_calcact3b(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                      lwn, gama, t, soln, frst, calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 13
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    !
             errinlocb = 0
             do ii = 1, 13
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients 
             do ii = 1, 13
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do        
    !
    !  ## Solve system of equations, with new activity coefficients 
             if (( .not. soln)) then
                lwnsq    = lwn*lwn
                gama10sq = gama(10)*gama(10)
    !
                a4 = k1*gama10sq/(gama(5)*gama(5))   
                a5 = k2*lwnsq/gama10sq               
                a6 = k3*lwnsq/(gama(11)*gama(11))    
    !
    !  ## 1. Calculate dissociation quantities 
                if (no3 >= tiny) then
                   psi5 = x3*no3/(a6/a5*(cl-x3) + x3)
                else 
                   psi5 = tiny
                end if 
    !
    !  ## 2. Account for NH3 evaporation
                if (so4 > tiny) then
                   bb = -(c4 + x3 + psi5 + 1.0d0/a4)
                   cc = c4*(psi5 + x3) - c3a/a4
    !               
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                   else
                      psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if 
                else
                   psi4 = tiny
                end if 
    !
    !  ## 3. Speciation
                nh4_t = c3a + psi4 
                cl_t  = x3        
                no3_t = psi5
                gnh3  = max(c4  - psi4, tiny2)
                ghno3 = max(no3 - psi5, tiny2)
                ghcl  = max(cl  - x3,   tiny2)
    !
    !  ## 4. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t
                scon = k4*lwnsq  
    !
                if (smin > tiny) then                        ! H+ in excess  
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                   end if 
    !
                   h   = max(hh,sqrt(scon))
    !              ohi = scon/h 
                else                                        ! OH- in excess   
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                   end if 
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if
    !       
    !  ## 5. Aerosol liquid water content 
                m4     = max(so4_t - c2, 0.0d0)		        ! (NH4)2SO4
                frnh4  = max(nh4_t - 2.0d0*m4, 0.0d0)
                m5     = min(no3_t, frnh4)				! NH4NO3
                frnh4  = max(frnh4 - m5, 0.0d0)
                m6     = min(cl_t, frnh4)				! NH4Cl
                lwn    = max(c5 + m4/zsr2 + m5/zsr3 + m6/zsr4, tiny)
             end if
          end do
    !     
          y3 = h*cl_t/ghcl/a6 - 1.0d0
    !
          condition = .false.
          if (noroot) then
    !  ## If no root on interval then do not perform ITP
             xa = x3
             xb = x3
          else if (y3 > 0.0d0 .and. ( .not. soln)) then
             xb = x3
             yb = y3
          else if (y3 < 0.0d0 .and. ( .not. soln)) then
             xa = x3
             ya = y3
          else if (( .not. soln)) then
             xa = x3
             xb = x3
          end if
    !
    !  ## Check for convergence criteria to exit ITP:
          if (xb - xa > abs(xa*eps) .and. ( .not. noroot)) then
             condition = .true.
             soln      = .false.
          else
             soln = .true.
          end if
    !
    ! ## Exit ITP if the function being bisected evaluates to a value <= eps; solution is assumed
          if (abs(y3) <= eps .and. ( .not. noroot)) then
             soln = .true.
             condition = .false.
          end if
    !
    ! ## Post-convergence correction: 
    ! ## If the calculated y-value (y3) after ITP is further from zero than an earlier iteration 
    ! ## then reset to the x-value/concentrations/activity coefficients that were found to minimize
    ! ## the objective function (i.e., y3); in this case, this is chosen as the solution
          if (( .not. condition) .and. ( .not. noroot) .and. abs(y3) > 0.1d0) then     
             if (abs(0.0d0 - y3_min) < abs(0.0d0 - y3) .and. abs(y3_min - y3) > 1.0d-1) then 
                x3    = x3_min
                cl_t  = cl_min
                nh4_t = nh4_min
                no3_t = no3_min
                h     = h_min
                lwn   = lwn_min
                ghcl  = ghcl_min
                ghno3 = ghno3_min
                gnh3  = gnh3_min
    ! 
    !  ## Reset activity coefficients 
                gama  = gama_min
                a6    = k3*(lwn/gama(11))*(lwn/gama(11))
                y3    = h*cl_t/ghcl/a6 - 1.0d0
             end if
          end if
       end do ! End outer loop of ITP search
    !
    !
    !  ### MINOR SYSTEM: HSO4-/SO42-/H+ ###
       hso4 = 0.0d0
       if (h > tiny .and. so4_t > tiny .and. lwn >= 1.0d-19) then
          ak1 = khso4*lwn/gama(7)*(gama(8)/gama(7))*(gama(8)/gama(7))
          bb  =-(h + so4_t + ak1)
          cc  = h*so4_t 
          dd  = bb*bb - 4.d0*cc
    ! 
          if (dd >= 0.0d0) then
             if (bb /= 0.d0) then
                dd = cc/(bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d3
             end if
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                hh = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
             else
                hh = 0.5d0*(-bb - sqrt(bb*bb - 4.d0*cc))
             end if 
    !
             hh    = max(tiny, min(hh, min(h, so4_t)))     ! To avoid negative H+   (i.e., if hh > h or hh > so4_t)
             h     = max(h - hh, 0.0d0)
             so4_t = max(so4_t - hh, 0.0d0)
             hso4  = hh
          end if 
       end if
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       na_i    = na_t
       cl_i    = cl_t
       nh3g_i  = gnh3
       hno3g_i = ghno3
       hclg_i  = ghcl
       h_i     = h
       lwn_i   = lwn
    ! 
       return
    end subroutine mach_hetp_calcg5
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: H6; Sulfate poor; sodium rich
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calch6(so4_i, nh4_i, nh3g_i, hno3g_i, hclg_i, hso4_i,          &
                                na_i, cl_i, no3_i, h_i, lwn_i, frna_d, rh, temp, k0,    &
                                p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: na_i    
       real(kind=8),    intent(inout) :: cl_i    
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: hclg_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i  
       real(kind=8),    intent  (out) :: frna_d   
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    !
    ! ## Local variables
       real(kind=8)     :: so4, nh4, hso4, gnh3, h, lwn, no3, cl, na, ghno3, ghcl, caso4
       real(kind=8)     :: t, aw, khso4, knh3, kh2o, khno3, khcl, frnh4, frno3, frcl, tt0
       real(kind=8)     :: bb, cc, dd, hh, v, ak1, errin, ohi, smin, scon, a5, psi4, psi5
       real(kind=8)     :: m5, m6, omehi, omebe, y1, y2, y3, x3, dx, a6, a4, u1, ya, yb
       real(kind=8)     :: xa, xb, so4_t, nh4_t, no3_t, na_t, cl_t, c5, c6, c7, sumzsr
       real(kind=8)     :: zsr1, zsr2, tt1, tt2, gmax, c1, c2, c3, k1, k2, k3, k4, loccon
       real(kind=8)     :: nh, sigma, xt, xf, xh, delta, rr, gx, gx2, errouloc, errinlocb
       real(kind=8)     :: frno3_d, frcl_d, lwnsq, gama10sq
       integer(kind=4)  :: irh, nmax, j, k, ii, rooteval
       logical(kind=4)  :: condition, noroot, earlyexit, soln, frst, calain, calou
       real(kind=8), dimension(13) :: gama, gamin, gamou, gama_min
       real(kind=8)      :: y3_min, x3_min, y3_lastiter, x3_lastiter
       real(kind=8)      :: no3_min, nh4_min, cl_min, h_min, lwn_min, ghcl_min, gnh3_min, ghno3_min
    !
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       na    = na_i
       cl    = cl_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       caso4 = 0.0d0
       gnh3  = 0.0d0
       ghno3 = 0.0d0
       ghcl  = 0.0d0
       h     = 0.0d0
       lwn   = tiny
       so4_t = so4
       nh4_t = 0.0d0
       no3_t = 0.0d0
       na_t  = 0.0d0
       cl_t  = 0.0d0
       frna_d= 0.0d0
       frcl_d= 0.0d0
       frno3_d= 0.0d0
       noroot=.false.
       soln  =.false.
       gama  = 0.1d0
       gamin = 1.0d10
       gamou = 1.0d10
       gmax  = 0.0d0
       earlyexit = .false.
       calou   = .true.
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*(tt2))
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*(tt2)))*                  &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    !
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*(tt2))
    !
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3 = k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. HCl(g) <==> H+(aq) + Cl-(aq)                               (xk3)
       khcl = k0(11) * exp(p1(11)*tt1 + p2(11)*(tt2))
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    
    !  ## Calculate constant parameters
       c1     = r*t
       c3     = 2.0d0*so4               !na2so4
       frna_d = max(na - c3, 0.0d0)                      
       c2     = min(frna_d, no3)        !nano3
       c5     = max(no3-c2, 0.0d0)      !frno3
       frno3_d= c5
       frna_d = max(frna_d - c2, 0.0d0)
       c7     = min(frna_d, cl)         !nacl
       c6     = max(cl - c7, 0.0d0)     !frcl
       frcl_d = c6
       frna_d = max(frna_d - c7, 0.0d0)   
       k1     = knh3/kh2o*c1
       k2     = khno3*c1
       k3     = khcl*c1
       k4     = kh2o*aw
    !
    !  ## Constant ZSR parameter
       sumzsr = c7/awsc(irh) + so4/awss(irh) + c2/awsn(irh) 
       zsr1   = awan(irh)
       zsr2   = awac(irh)
    !
    !  ## Initial speciation
       na_t  = c2 + c7 + c3
    !
    !
    !  ### STAGE 1: Root tracking ###
    !  ## Find a subinterval [xa,xb] on the larger interval [I1,I2] where a sign change occurs
       rooteval = 0
       condition = .true.
       do while (rooteval < 2 .or. (condition .and. rooteval < ndiv + 1))   ! Begin outer loop for root tracking
          rooteval = rooteval + 1
    !
    !  ## Set high limit for root tracking (i.e. lower bound)
          if (rooteval == 1) then
            omehi = tiny
            y1    = 1.0d0
          end if
    !
    !  ## Begin search on subinterval 
          if (rooteval == 2) then
             soln = .false.
             if (abs(y2) <= eps) then
                earlyexit = .true.
             end if
             y1 = y2
    !
             if (earlyexit) then
                dx = 0.0d0
             else
                dx = (c6-tiny-tiny)/float(ndiv)
             end if
             omebe = omehi            ! Lower bound of subinterval (xa)
             omehi = omehi + dx       ! Upper bound of subinterval (xb)
          end if
    !
    !  ## Continue search 
          if (rooteval > 2) then
             if (loccon < 0.0d0) then
    !  ## 1. Root has been found on the subinterval; save x values for ITP search
                y1    = y1
                omebe = omebe
                omehi = omehi
             else
    !  ## 2. No root has been found, continue searching in the next subinterval
                y1    = y2
                omebe = omehi
                omehi = omehi + dx
             end if
          end if
    !
    !  ## Solve the system of equations 
          frst   = .true.
          calain = .true.
          if (( .not. soln)) then
             lwnsq = lwn*lwn
             gama10sq = gama(10)*gama(10)
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))   
             a5 = k2*lwnsq/gama10sq                
             a6 = k3*lwnsq/(gama(11)*gama(11))     
    !
    !  ## 1. Calculate dissociation quantities
             psi5 = c5*(omehi + c7) - a6/a5*c2*(c6-omehi)
             psi5 = max(psi5 / (a6/a5*(c6 - omehi) + omehi + c7), tiny)
    !
             if (nh4 > tiny .and. lwn > tiny) then
                bb = -(nh4 + omehi + psi5 + 1.0d0/a4)
                cc = nh4*(psi5 + omehi)
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                else
                   psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                psi4 = min(psi4, nh4)
             else
                psi4 = tiny
             end if
    !
    !  ## 2. Speciation
             nh4_t = psi4
             cl_t  = omehi + c7
             no3_t = psi5  + c2
             gnh3  = max(nh4 - psi4,  tiny2)
             ghno3 = max(c5  - psi5,  tiny2)
             ghcl  = max(c6  - omehi, tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t
             scon = k4*lwnsq
    !
             if (smin > tiny) then                        ! H+ in excess
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                h   = max(hh,sqrt(scon))
    !           ohi = scon/h
             else                                        ! OH- in excess
                bb = smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                ohi = max(hh,sqrt(scon))
                h   = scon/ohi
             end if
    !
    !  ## 4. Aerosol liquid water content
             frno3  = max(no3_t - c2, 0.0d0)
             frcl   = max(cl_t  - c7, 0.0d0)
             m5     = min(nh4_t, frno3)
             frnh4  = max(nh4_t - m5 , 0.0d0)
             m6     = min(frcl, frnh4)
             lwn    = max(sumzsr + m5/zsr1 + m6/zsr2, tiny)
          end if 
    !
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)
             k = k + 1
    !
    !  ## Reset gamou
             if (( .not. soln) .and. frst) then
                gamou(1)  = gama(1)
                gamou(2)  = gama(2)
                gamou(3)  = gama(3)
                gamou(4)  = gama(4)
                gamou(5)  = gama(5)
                gamou(6)  = gama(6)
                gamou(7)  = gama(7)
                gamou(8)  = gama(8)
                gamou(9)  = gama(9)
                gamou(10) = gama(10)
                gamou(11) = gama(11)
                gamou(12) = gama(12)
                gamou(13) = gama(13)
             end if 
    !
    !  ## Reset gamin
             if (( .not. soln)) then
                gamin(1)  = gama(1)
                gamin(2)  = gama(2)
                gamin(3)  = gama(3)
                gamin(4)  = gama(4)
                gamin(5)  = gama(5)
                gamin(6)  = gama(6)
                gamin(7)  = gama(7)
                gamin(8)  = gama(8)
                gamin(9)  = gama(9)
                gamin(10) = gama(10)
                gamin(11) = gama(11)
                gamin(12) = gama(12)
                gamin(13) = gama(13)
             end if 
    !
             call mach_hetp_calcact3b(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                     lwn, gama, t, soln, frst, calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 13
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    !
             errinlocb = 0
             do ii = 1, 13
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients
             do ii = 1, 13
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do
    !
    !  ## Solve system of equations, using new activity coefficients 
             if (( .not. soln)) then
                lwnsq    = lwn*lwn 
                gama10sq = gama(10)*gama(10)
    !
                a4 = k1*gama10sq/(gama(5)*gama(5))    
                a5 = k2*lwnsq/gama10sq                
                a6 = k3*lwnsq/(gama(11)*gama(11))              
    !
    !  ## 1. Calculate dissociation quantities
                psi5 = c5*(omehi + c7) - a6/a5*c2*(c6-omehi)
                psi5 = max(psi5 / (a6/a5*(c6 - omehi) + omehi + c7), tiny)
    !           
                bb = -(nh4 + omehi + psi5 + 1.0d0/a4)
                cc = nh4*(psi5 + omehi)
                if (nh4 > tiny .and. lwn > tiny) then
                   bb = -(nh4 + omehi + psi5 + 1.0d0/a4)
                   cc = nh4*(psi5+ omehi)
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                   else
                      psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   psi4 = min(psi4, nh4)
                else
                   psi4 = tiny
                end if
    !
    !  ## 2. Speciation
                nh4_t = psi4
                cl_t  = omehi + c7
                no3_t = psi5 + c2
                gnh3  = max(nh4 - psi4,  tiny2)
                ghno3 = max(c5  - psi5,  tiny2)
                ghcl  = max(c6  - omehi, tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t
                scon = k4*lwnsq  
    !
                if (smin > tiny) then                        ! H+ in excess
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   h   = max(hh,sqrt(scon))
    !              ohi = scon/h
                else                                        ! OH- in excess
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if
    !
    !  ## 4. Aerosol liquid water content
                frno3  = max(no3_t - c2, 0.0d0)
                frcl   = max(cl_t  - c7, 0.0d0)
                m5     = min(nh4_t, frno3)
                frnh4  = max(nh4_t - m5 , 0.0d0)
                m6     = min(frcl, frnh4)
                lwn    = max(sumzsr + m5/zsr1 + m6/zsr2, tiny)
             end if
          end do
    !
          y2 = nh4_t*cl_t/ghcl/gnh3/a6/a4 - 1.0d0    !Function value
    !
    !  ## Check for criteria to exit root tracking 
          condition = .false.
          loccon = sign(1.0d0,y1)*sign(1.0d0,y2)
          if (loccon > 0.0d0 .and. ( .not. noroot) .and. abs(y2) > eps .and. c6 > tiny) then
             condition = .true.
          elseif (loccon < 0.0d0 .and. abs(y2) > eps) then
    !  ## Interval had been found where sign change occurs; exit root tracking and proceed to ITP
             soln = .true.
          else
    !  ## abs(y2) <= eps; solution is assumed; exit root tracking and proceed to minor system (no ITP)
             soln   = .true.
             noroot = .true.
          end if 
    !
    !  ## Too little frcl; reset x-value to tiny and exit root tracking 
          if (c6 <= tiny) then
             condition = .false.
             rooteval  = 2       ! Increment rooteval by 1 to force exit from root tracking
             noroot    = .true.
             omehi     = tiny    ! Reset x-value to tiny, and solve system with this value
             omebe     = omehi
          end if 
    !
    !  ### AFTER iterating through ALL ndiv subdivided intervals
          if (rooteval == ndiv + 1) then
             if (loccon > 0.0d0 .and. abs(y2) > eps .and. c6 > tiny .and. ( .not. earlyexit)) then
    !  ## (1) No solution
                noroot = .true.
                omehi = tiny    ! Reset to tiny
                omebe = omehi
    !            write(*,*), 'Warning in CALCH6: no solution found'
             else if (loccon > 0.0d0 .and. abs(y2) <= eps .and. c6 > tiny .and. ( .not. earlyexit)) then
    !  ## (2) Solution is assumed and ITP is not required
                noroot = .true.
                !soln = .true.
             end if
          end if
       end do  !End outer loop for root tracking
    !
    !
    !
    !  ### STAGE 2: modified bisection search (using ITP algorithm) ###
    !  ## Initialize static ITP variables 
       if (( .not. noroot)) then
          ya = y1
          yb = y2
          xa = omebe
          xb = omehi
          x3 = omehi
    !
          if (xa == xb) then
             noroot = .true.
             gx = tiny
          else
             gx = xb - xa
          end if 
    !
          gx2  = (xa+xb)*0.5d0
          u1   = 0.2d0/gx
          nh   = log10(abs(gx/(2.0d0*eps*gx2))) / log10of2   
          nmax = int(nh) + 2                               
       else
          x3 = omehi
       end if
    !
    !  ## Start search
    !  1. Store value of previous y3 for later use
       y3_lastiter = 0.0d0
       x3_lastiter = 0.0d0
       y3_min = 1.0d20
       x3_min = 0.0d0
    !
       j = 0
       condition = .true.
       soln      = .false.
       do while (j < maxit .and. condition)   ! Begin outer loop for ITP search
    !  ## Set dynamic ITP variables 
    !  1. Track x3 and y3 of the previous iteration
          if (j > 0) then
             y3_lastiter = y3
             x3_lastiter = x3
          end if
    !
    !  2. Track the minimum y3 that is found before ending
          if (abs(0.0d0 - y3_min) > abs(0.0d0 - y3_lastiter) .and. j > 0) then
             y3_min    = y3_lastiter
             x3_min    = x3_lastiter
             ghno3_min = ghno3
             no3_min   = no3_t
             h_min     = h
             ghcl_min  = ghcl
             cl_min    = cl_t
             gnh3_min  = gnh3
             nh4_min   = nh4_t
             lwn_min   = lwn
             gama_min  = gama
          end if
    !
          if ((.not. noroot) .and. (.not. soln)) then
             if (yb - ya == 0.0d0) then
                write(*,*), '######        ABORT       ######'
                write(*,*), 'Zero divide in ITP reset: CALCH6'
                write(*,*), 'SO4 in = ', so4
                write(*,*), 'NH4 in = ', nh4
                write(*,*), 'NO3 in = ', no3
                write(*,*), 'Na in  = ', na
                write(*,*), 'Cl in  = ', cl
                write(*,*), 'Temp in= ', t
                write(*,*), 'RH in  = ', aw
                return
             end if  
    !
             gx   = xb - xa
             xh   = 0.5d0*(xa + xb)
             rr   = max(gx2*eps*2.d0**(real(nmax - j)) - 0.5d0*gx, 0.0d0)
             delta= u1*(max(gx, 0.0d0))**2.0d0   
             xf   = max((yb*xa - ya*xb) / (yb - ya), 0.0d0)
    !
             sigma = sign(1.0d0, xh - xf)
             if (delta <= abs(xh - xf)) then
                xt = xf + sigma*delta           
             else
                xt = xh
             end if
    !
             if (abs(xt - xh) <= rr) then
                x3 = xt
             else
                x3 = xh - sigma*rr
             end if 
          end if 
    !
          j = j + 1
    !
          if (( .not. soln)) then
             gmax = 0.1d0
             gmax = max(gmax, gama(1))
             gmax = max(gmax, gama(2))
             gmax = max(gmax, gama(3))
             gmax = max(gmax, gama(4))
             gmax = max(gmax, gama(5))
             gmax = max(gmax, gama(6))
             gmax = max(gmax, gama(7))
             gmax = max(gmax, gama(8))
             gmax = max(gmax, gama(9))
             gmax = max(gmax, gama(10))
             gmax = max(gmax, gama(11))
             gmax = max(gmax, gama(12))
             gmax = max(gmax, gama(13))
          end if
    !
    !  ## Reinitialize activity coefficients if gmax > 100.0d0
          if ((gmax > 100.0d0) .and. (( .not. soln))) then
             gama  = 0.1d0
             gamin = 1.0d10
             gamou = 1.0d10
             calou = .true.
             frst  = .true.
          end if
    !
    !  ## Solve system of equations
          frst    = .true.
          calain  = .true.
          if (( .not. soln)) then
             lwnsq    = lwn*lwn
             gama10sq = gama(10)*gama(10)
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))   
             a5 = k2*lwnsq/gama10sq                
             a6 = k3*lwnsq/(gama(11)*gama(11))     
    !
    !  ## 1. Calculate dissociation quantities
             psi5 = c5*(x3 + c7) - a6/a5*c2*(c6-x3)
             psi5 = max(psi5 / (a6/a5*(c6 - x3) + x3 + c7), tiny)
    !
             if (nh4 > tiny .and. lwn > tiny) then
                bb = -(nh4 + x3 + psi5 + 1.0d0/a4)
                cc = nh4*(psi5+ x3)
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                else
                   psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                psi4 = min(psi4, nh4)
             else
                psi4 = tiny
             end if
    !
    !  ## 2. Speciation
             nh4_t = psi4
             cl_t  = x3 + c7
             no3_t = psi5 + c2
             gnh3  = max(nh4 - psi4,  tiny2)
             ghno3 = max(c5  - psi5,  tiny2)
             ghcl  = max(c6  - x3, tiny2)
    
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t
             scon = k4*lwnsq  
    !
             if (smin > tiny) then                        ! H+ in excess
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                h   = max(hh,sqrt(scon))
    !           ohi = scon/h
             else                                        ! OH- in excess
                bb = smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                ohi = max(hh,sqrt(scon))
                h   = scon/ohi
             end if
    !
    !  ## 4. Aerosol liquid water content
             frno3  = max(no3_t - c2, 0.0d0)
             frcl   = max(cl_t  - c7, 0.0d0)
             m5     = min(nh4_t, frno3)
             frnh4  = max(nh4_t - m5 , 0.0d0)
             m6     = min(frcl, frnh4)
             lwn    = max(sumzsr + m5/zsr1 + m6/zsr2, tiny)
          end if 
    !
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)
             k = k + 1
    !
    !  ## Reset gamin and gamou
             if ((( .not. soln)) .and. (frst)) then
                gamou(1)  = gama(1)
                gamou(2)  = gama(2)
                gamou(3)  = gama(3)
                gamou(4)  = gama(4)
                gamou(5)  = gama(5)
                gamou(6)  = gama(6)
                gamou(7)  = gama(7)
                gamou(8)  = gama(8)
                gamou(9)  = gama(9)
                gamou(10) = gama(10)
                gamou(11) = gama(11)
                gamou(12) = gama(12)
                gamou(13) = gama(13)
             end if 
    !
             if (( .not. soln)) then
                gamin(1)  = gama(1)
                gamin(2)  = gama(2)
                gamin(3)  = gama(3)
                gamin(4)  = gama(4)
                gamin(5)  = gama(5)
                gamin(6)  = gama(6)
                gamin(7)  = gama(7)
                gamin(8)  = gama(8)
                gamin(9)  = gama(9)
                gamin(10) = gama(10)
                gamin(11) = gama(11)
                gamin(12) = gama(12)
                gamin(13) = gama(13)
             end if 
    !
             call mach_hetp_calcact3b(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                     lwn, gama, t, soln, frst, calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 13
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    !
             errinlocb = 0
             do ii = 1, 13
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients
             do ii = 1, 13
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do
    !
    !  ## Solve system of equations, with new activity coefficients 
             if (( .not. soln)) then
                lwnsq    = lwn*lwn
                gama10sq = gama(10)*gama(10)
    ! 
                a4 = k1*gama10sq/(gama(5)*gama(5))    
                a5 = k2*lwnsq/gama10sq                
                a6 = k3*lwnsq/(gama(11)*gama(11))     
    !
    !  ## 1. Calculate dissociation quantities
                psi5 = c5*(x3 + c7) - a6/a5*c2*(c6-x3)
                psi5 = max(psi5 / (a6/a5*(c6 - x3) + x3 + c7), tiny)
    !
                if (nh4 > tiny .and. lwn > tiny) then
                   bb = -(nh4 + x3 + psi5 + 1.0d0/a4)
                   cc = nh4*(psi5 + x3)
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                   else
                      psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   psi4 = min(psi4, nh4)
                else
                   psi4 = tiny
                end if
    !
    !  ## 2. Speciation
                nh4_t = psi4
                cl_t  = x3 + c7
                no3_t = psi5 + c2
                gnh3  = max(nh4 - psi4, tiny2)
                ghno3 = max(c5  - psi5, tiny2)
                ghcl  = max(c6  - x3,   tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t
                scon = k4*lwnsq  
    !
                if (smin > tiny) then                        ! H+ in excess
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   h   = max(hh,sqrt(scon))
    !              ohi = scon/h
                else                                        ! OH- in excess
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if
    !
    !  ## 4. Aerosol liquid water content
                frno3  = max(no3_t - c2, 0.0d0)
                frcl   = max(cl_t  - c7, 0.0d0)
                m5     = min(nh4_t, frno3)
                frnh4  = max(nh4_t - m5 , 0.0d0)
                m6     = min(frcl, frnh4)
                lwn    = max(sumzsr + m5/zsr1 + m6/zsr2, tiny)
             end if   
          end do
    !
          y3 = nh4_t*cl_t/ghcl/gnh3/a6/a4 - 1.0d0    !Function value
    !
          condition = .false.
          if (noroot) then
    !  ## If no root on interval then do not perform ITP
             xa = x3
             xb = x3
          else if ((y3 > 0.0d0) .and. (( .not. soln))) then
             xb = x3
             yb = y3
          else if ((y3 < 0.0d0) .and. (( .not. soln))) then
             xa = x3
             ya = y3
          else if (( .not. soln)) then
             xa = x3
             xb = x3
          end if
    !
    !  ## Check for convergence criteria to exit ITP:
          if ((xb - xa > abs(xa*eps)) .and. (( .not. noroot))) then
             condition = .true.
             soln      = .false.
          else
             soln      = .true.
          end if
    !
    ! ## Exit ITP if the function being bisected evaluates to a value <= eps; solution is assumed
          if ((abs(y3) <= eps) .and. (( .not. noroot))) then
             soln = .true.
             condition = .false.
          end if
    !
    ! ## Post-convergence correction: 
    ! ## If the calculated y-value (y3) after ITP is further from zero than an earlier iteration 
    ! ## then reset to the x-value/concentrations/activity coefficients that were found to minimize
    ! ## the objective function (i.e., y3); in this case, this is chosen as the solution
          if ((( .not. condition)) .and. (( .not. noroot)) .and. (abs(y3) > 0.1d0)) then     
             if (abs(0.0d0 - y3_min) < abs(0.0d0 - y3) .and. abs(y3_min - y3) > 1.0d-1) then
                x3    = x3_min
                cl_t  = cl_min
                nh4_t = nh4_min
                no3_t = no3_min
                h     = h_min
                lwn   = lwn_min
                ghcl  = ghcl_min
                ghno3 = ghno3_min
                gnh3  = gnh3_min
    ! 
    !  ## Reset activity coefficients 
                gama = gama_min
                a4   = k1*(gama(10)/gama(5))*(gama(10)/gama(5))
                a6   = k3*(lwn/gama(11))*(lwn/gama(11))
                y3   = nh4_t*cl_t/ghcl/gnh3/a6/a4 - 1.0d0
             end if
          end if
       end do ! End outer loop of ITP search
    !
    !
    !  ### MINOR SYSTEM: HSO4-/SO42-/H+ ###
       hso4 = 0.0d0
       if (h > tiny .and. so4_t > tiny .and. lwn >= 1.0d-19) then
          ak1 = khso4*lwn/gama(7)*(gama(8)/gama(7))*(gama(8)/gama(7))
          bb  =-(h + so4_t + ak1)
          cc  = h*so4_t 
          dd  = bb*bb - 4.d0*cc
    ! 
          if (dd >= 0.0d0) then
             if (bb /= 0.d0) then
                dd = cc/(bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d3
             end if
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                hh = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
             else
                hh = 0.5d0*(-bb - sqrt(bb*bb - 4.d0*cc))
             end if 
    !
             hh    = max(tiny, min(hh, min(h, so4_t)))     ! To avoid negative H+   (i.e., if hh > h or hh > so4_t)
             h     = max(h - hh, 0.0d0)
             so4_t = max(so4_t - hh, 0.0d0)
             hso4  = hh
          end if 
       end if
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       na_i    = na_t
       cl_i    = cl_t
       nh3g_i  = gnh3
       hno3g_i = ghno3
       hclg_i  = ghcl
       h_i     = h
       lwn_i   = lwn
    !
       return
    end subroutine mach_hetp_calch6
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: I6; Sulfate rich, no free acid
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calci6(so4_i, nh4_i, nh3g_i, hno3g_i, hclg_i, hso4_i,          &
                                na_i, cl_i, no3_i, h_i, lwn_i, frna, rh, temp, k0,      &
                                p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: na_i    
       real(kind=8),    intent(inout) :: cl_i    
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: hclg_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i  
       real(kind=8),    intent  (out) :: frna 
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables
       real(kind=8)     :: so4, nh4, hso4, gnh3, h, lwn
       real(kind=8)     :: no3, cl, na, ghno3, ghcl, caso4, t, aw
       real(kind=8)     :: so4_t, nh4_t, no3_t, na_t, cl_t, m1, m2, m3, delcl
       real(kind=8)     :: na2so4, nh42s4, nh4hs4, nahso4, slc, tt1, tt2, c1, c2, c3
       real(kind=8)     :: knh3, kh2o, khno3, khcl, khso4, bb, cc, dd, ff, v, errin
       real(kind=8)     :: a3, a4, a6, psi6, frnh4, frso4, delno, tt0
       integer(kind=4)  :: j, ii, islv, irh
       logical(kind=4)  :: ispoly
       real(kind=8), dimension(13) :: gama, gamin
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       na    = na_i
       cl    = cl_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       gnh3  = 0.0d0
       ghno3 = 0.0d0
       ghcl  = 0.0d0
       h     = 0.0d0
       caso4 = 0.0d0
       lwn   = tiny
       so4_t = 0.0d0
       nh4_t = 0.0d0
       no3_t = 0.0d0
       na_t  = 0.0d0
       cl_t  = 0.0d0
       delcl = 0.0d0
       na2so4= 0.0d0 
       nh42s4= 0.0d0
       nh4hs4= 0.0d0
       nahso4= 0.0d0
       slc   = 0.0d0
       gama  = 0.1d0
       gamin = 1.0d10
    !
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
                 (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    ! 
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !   
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3 = k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. HCl(g) <==> H+(aq) + Cl-(aq)                               (xk3)
       khcl = k0(11) * exp(p1(11)*tt1 + p2(11)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1) 
    !
    !  ## Calculate non-volatile solids (subroutine 'CALCI1A' in ISORROPIA II)
       na2so4  = 0.5d0*na                                  !na2so4
       frna    = na - 2.0d0*na2so4
       frso4   = max(so4 - na2so4, 0.0d0)
       slc     = min(nh4/3.0d0, frso4*0.5d0)               !(nh4)3h(so4)2
       frso4   = max(frso4 - 2.0d0*slc, 0.0d0)
       frnh4   = max(nh4 - 3.0d0*slc, 0.0d0)
    !      
       if (frso4 <= tiny) then
          slc     = max(slc - frnh4, 0.0d0)        
          nh42s4  = 2.0d0*frnh4                            !(nh4)2so4
       else if (frnh4 <= tiny .and. na2so4 <= tiny) then
          nh4hs4  = 3.0d0*min(frso4, slc)                  !nh4hso4
          slc     = max(slc - frso4, 0.0d0)
       else if (frnh4 <= tiny .and. na2so4 > tiny) then
          nh4hs4  = 3.0d0*min(frso4, slc)
          slc     = max(slc - frso4, 0.0d0)
          frso4   = max(frso4-nh4hs4/3.0d0, 0.0d0)
          nahso4  = 2.0d0*frso4                            !nahso4
          na2so4  = max(na2so4 - frso4, 0.0d0)
          frna    = max(na - nahso4 - 2.0d0*na2so4, 0.0d0)
       end if 
    !
    !  ## Calculate gaseous species 
       ghno3 = no3
       ghcl  = cl 
    
    !  ## Calculate constants
       c1 = slc + nahso4 + nh4hs4
       c2 = slc + na2so4 + nh42s4
    !
    !  ### MAJOR SYSTEM H+/HSO4-/SO42- ###
    !  ## Setup initial conditions
    !  ## 1. Calculate dissociation quantities
       a6 = khso4*1.0d-19
       bb = c2 + a6  
       cc = -a6*c1
    !
       if (bb /= 0.d0) then
          dd = cc/(bb*bb)
          v  = 4.d0*dd
       else
          v  = 1.d3
       end if
    !
       if (abs(v) <= smrt .and. bb > 0.d0) then
          psi6 = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
       else
          psi6 = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
       end if 
    !
    !  ## 2. Speciation 
       h     = psi6
       na_t  = 2.0d0*na2so4 + nahso4
       nh4_t = 3.0d0*slc + 2.0d0*nh42s4 + nh4hs4
    !
    !  Due to round-off, psi6 may be .gt. lc + nahso4 + nh4hs4 giving negative hso4
    !  If this condition is true, reset psi6 to [c1] and proceed
       psi6  = min(psi6, c1)
       so4_t = c2 + psi6
       hso4  = max(c1 - psi6, 0.0d0)
    !
    !  ## 3. Aerosol liquid water content 
       lwn = max(nh42s4/awas(irh) + na2so4/awss(irh) + nh4hs4/awab(irh) + nahso4/awsb(irh) +  &
                 slc/awlc(irh), tiny)
       c3  = khso4*lwn
    !
    !
    !  ## Iterative search for solution with convergence of activity coefficients
       j     = 0
       errin = 1.0d0
       do while (j < nsweep-1 .and. errin >= epsact)
          j = j + 1
    !
    !  ## Reset gamin
          gamin(1)  = gama(1)
          gamin(2)  = gama(2)
          gamin(3)  = gama(3)
          gamin(4)  = gama(4)
          gamin(5)  = gama(5)
          gamin(6)  = gama(6)
          gamin(7)  = gama(7)
          gamin(8)  = gama(8)
          gamin(9)  = gama(9)
          gamin(10) = gama(10)
          gamin(11) = gama(11)
          gamin(12) = gama(12)
          gamin(13) = gama(13)
    !
          call mach_hetp_calcact3(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                  lwn, gama, t)
    !
          errin = 0.0d0
    !  ## Test for convergence of activity coefficients 
          do ii = 1, 13
             errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
          end do   
    !
    !  ## 1. Calculate dissociation quantities
          a6 = c3/gama(7)*(gama(8)/gama(7))*(gama(8)/gama(7))
          bb = c2 + a6   
          cc = -a6*c1
    !
          if (bb /= 0.d0) then
             dd = cc/(bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d3
          end if
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             psi6 = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             psi6 = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
          end if 
    !
    !  ## 2. Speciation 
          h     = psi6
          psi6  = min(psi6, c1)
          so4_t = c2 + psi6
          hso4  = max(c1 - psi6, 0.0d0)
       end do
    !
    !
    !  ### MINOR SYSTEM #1: Cl-/HCl/NO3-/HNO3/H+ ###
    !  ## Calculate dissolution of HCl, HNO3 in the presence of (H,SO4)
    !  ## HCl, HNO3 are considered minor species 
       ispoly = .false.
    !
    !  ## 1. Special case: aerosol liquid water content (lwn) = 0.0d0
       if (lwn <= tiny) then
    !  Gaseous species
          ghcl  = max(cl  - cl_t,  0.0d0)
          ghno3 = max(no3 - no3_t, 0.0d0) 
          m1    = 0.0d0
          m2    = 0.0d0
          m3    = 0.0d0
    !
    !  ## 2. Special case: HCl = HNO3 = 0.0d0
       else if (cl <= tiny .and. no3 <= tiny) then
          m1  = 0.0d0
          m2  = 0.0d0
          m3  = 0.0d0
    !
    !  ## 3. Special case: HCl = 0.0d0
       else if (cl <= tiny) then
    !  Nitric acid in the liquid phase is assumed a minor species
    !  HNO3(g) <--> (H+) + (NO3-), using (H+) from the sulfates
          ff = 0.0d0
          m1 = 0.0d0
          m2 = 0.0d0
          m3 = 0.0d0
          if (lwn > tiny) then 
             bb = h + khno3*r*t*(lwn/gama(10))**2.0
             cc = (khno3*r*t*(lwn/gama(10))**2.0)*no3
    !
             if (bb /= 0.d0) then
                dd = cc / (bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d+03
             end if       
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                ff = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
             else
                ff = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
             end if 
          end if 
    !
    !  Speciation  
          ghno3 = max(no3 - ff, 0.0d0)
          no3_t = ff
          h     = h + ff
    !
    !  ## 4. Special case: HNO3 = 0.0d0
       else if (no3 <= tiny) then
    !  Hydrochloric acid in the liquid phase is assumed a minor species
    !  HCl (g) <--> (H+) + (Cl-) using (H+) from the sulfates
          m1 = 0.0d0
          m2 = 0.0d0
          m3 = 0.0d0
          ff = 0.0d0
          if (lwn > tiny) then
             bb = h + khcl*r*t*(lwn/gama(11))**2.0
             cc = (khcl*r*t*(lwn/gama(11))**2.0)*cl
    !
             if (bb /= 0.d0) then
                dd = cc / (bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d+03
             end if       
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                ff = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
             else
                ff = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
             end if 
          end if 
    !
    !  Speciation        
          ghcl = max(cl - ff, 0.0d0) 
          cl_t = ff 
          h    = h + ff
    !
    !  ## 5. All else: not a special case
       else
          ispoly = .true.
          a3 = khno3*r*t*(lwn/gama(10))**2.0  !HNO3
          a4 = khcl*r*t*(lwn/gama(11))**2.0   !HCl
    !
    !  Calculate cubic equation coefficients
          m1 = (a3*no3 + a4*cl + (h + a4)*(a3 - a4))/(a3 - a4)
          m2 = ((h + a4)*a4*cl - a4*(a3 - a4)*cl)/(a3 - a4)
          m3 = -a4*a4*cl*cl/(a3 - a4)
       end if 
    !
    !  ## Calculate roots 
    !    call mach_hetp_poly3v(m1, m2, m3, delcl, islv)     !# Original ISORROPIA subroutine 
        call mach_hetp_poly(m1, m2, m3, cl, delcl, islv)  
    !   
       delno = 0.0d0
       if (ispoly) then
          if (islv /= 0) then
    !  Tiny amounts of HCl are assumed when there is no root
             delcl = tiny                                        !Change in Cl-
          end if
    !
          a3    = khno3*r*t*(lwn/gama(10))**2.0                  !HNO3
          a4    = khcl*r*t*(lwn/gama(11))**2.0                   !HCl
          delno = min(a3*no3*delcl/(a4*cl + (a3 - a4)*delcl), no3) !Change in NO3-
    !
          if (delcl < 0.0d0 .or. delno < 0.0d0 .or. delcl > cl .or. delno > no3) then
             delcl = tiny                                        !Change in Cl-
             delno = tiny                                        !Change in NO3-
          end if 
    !
    !  ## Effect on liquid phase
          h     = h + delno + delcl       ! H+   change
          cl_t  = cl_t  + delcl           ! Cl-  change
          no3_t = no3_t + delno           ! NO3- change
    !
    !  ## Gaseous species
          ghcl  = max(cl  - cl_t,  0.0d0)
          ghno3 = max(no3 - no3_t, 0.0d0) 
       end if 
    !
    !
    !  ### MINOR SYSTEM #2: NH4+/NH3/H+ ###
    !  Ammonia in the gas phase is assumed a minor species that does not significantly perturb 
    !  the aerosol equilibrium: NH3(g) + H+(aq) <==> NH4+(aq); H+ determined from the major
    !  system (solved above) is used initially.  
       if (lwn > tiny) then
    !  ## Calculate NH3 sublimation
          ff = 0.0d0
          a3 = (knh3/kh2o)*r*t*(gama(10)/gama(5))**2.0
          bb = h + 1.0d0/a3    ! bb always > 0
          cc = -nh4_t/a3
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d+03
          end if  
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             ff = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             ff = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
          end if 
    !
    !  ## Speciation 
    !  Due to round-off, ff may be .gt. nh4_t giving negative nh4_t
    !  If this condition is true, then set ff = nh4_t
          ff    = max(tiny, min(ff, nh4_t))
          gnh3  = ff
          nh4_t = max(nh4_t - ff, 0.0d0)
          h     = h + ff
       end if
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       na_i    = na_t
       cl_i    = cl_t
       nh3g_i  = gnh3
       hno3g_i = ghno3
       hclg_i  = ghcl
       h_i     = h
       lwn_i   = lwn
    !
       return
    end subroutine mach_hetp_calci6
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: J3; Sulfate rich, free acid
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcj3(so4_i, nh4_i, nh3g_i, hno3g_i, hclg_i, hso4_i,          &
                                na_i, cl_i, no3_i, h_i, lwn_i, rh, temp, k0,            &
                                p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: na_i    
       real(kind=8),    intent(inout) :: cl_i    
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: hclg_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables
       real(kind=8)     :: so4, nh4, hso4, gnh3, h, lwn
       real(kind=8)     :: no3, cl, na, ghno3, ghcl, caso4, t, aw
       real(kind=8)     :: knh3, khso4, kh2o, khno3, khcl
       real(kind=8)     :: bb, cc, dd, ff, v, errin, a3, a4, psi6, delno, tt0, c3, c5, c6
       real(kind=8)     :: delcl, m1, m2, m3, tt1, tt2, c1, c2, so4_t, nh4_t, no3_t, na_t, cl_t
       integer(kind=4)  :: j, ii, islv, irh
       logical(kind=4)  :: ispoly
       real(kind=8), dimension(13) :: gama, gamin
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       na    = na_i
       cl    = cl_i
       aw    = rh
       t     = temp
    !
       m1    = 0.0d0
       m2    = 0.0d0
       m3    = 0.0d0
       hso4  = 0.0d0
       gnh3  = 0.0d0
       ghno3 = 0.0d0
       ghcl  = 0.0d0
       h     = 0.0d0
       caso4 = 0.0d0
       lwn   = tiny
       so4_t = 0.0d0
       nh4_t = 0.0d0
       no3_t = 0.0d0
       na_t  = 0.0d0
       cl_t  = 0.0d0
       delcl = 0.0d0
       gama  = 0.1d0
       gamin = 1.0d10
    !
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
                 (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    ! 
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !   
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3 = k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. HCl(g) <==> H+(aq) + Cl-(aq)                               (xk3)
       khcl = k0(11) * exp(p1(11)*tt1 + p2(11)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Constant values
       c1 = max(so4 - nh4 - na, tiny)
       c2 = c1 + na + nh4 
       c3 = nh4/awab(irh) + na/awsb(irh)
       c5 = nh4 + na
       c6 = awsa(irh)
    !
    !
    !  ### MAJOR SYSTEM H+/HSO4-/SO42- ###
    !  ## Setup initial conditions
    !  ## 1. Calculate dissociation quantities
       a3 = khso4*1.0d-19
       bb = a3 + c1             ! bb always > 0
       cc = -a3*c2
    !
       if (bb /= 0.d0) then
          dd = cc/(bb*bb)
          v  = 4.d0*dd
       else
          v  = 1.d3
       end if
    !
       if (abs(v) <= smrt .and. bb > 0.d0) then
          psi6 = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
       else
          psi6 = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
       end if 
    !
    !  ## 2. Speciation
    !  Due to round-off, subtraction of psi6 may give negative hso4; if this condition is true reset psi6
       psi6  = min(psi6, c2)
       h     = c1 + psi6 
       na_t  = na
       nh4_t = nh4
       so4_t = psi6
       hso4  = max(c2 - psi6, 0.0d0)
    !
    !  ## 3. Aerosol liquid water content
       lwn = max(c3 + max((so4_t + hso4 - c5), 0.0d0)/c6, tiny)
    !
    !
    !  ## Iterative search for solution with convergence of activity coefficients
       j     = 0
       errin = 1.0d0
       do while (j < nsweep-1 .and. errin >= epsact)
          j = j + 1
    !
    !  ## Reset gamin
          gamin(1)  = gama(1)
          gamin(2)  = gama(2)
          gamin(3)  = gama(3)
          gamin(4)  = gama(4)
          gamin(5)  = gama(5)
          gamin(6)  = gama(6)
          gamin(7)  = gama(7)
          gamin(8)  = gama(8)
          gamin(9)  = gama(9)
          gamin(10)  = gama(10)
          gamin(11)  = gama(11)
          gamin(12)  = gama(12)
          gamin(13)  = gama(13)
    !
          call mach_hetp_calcact3(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t, lwn, gama, t)
    !
          errin = 0.0d0
    !  ## Test for convergence of activity coefficients 
          do ii = 1, 13
             errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
          end do   
    !
    !  ## 1. Calculate dissociation quantities
          a3   = khso4*lwn/gama(7)*(gama(8)/gama(7))*(gama(8)/gama(7))
          bb   = a3 + c1                    !! bb always > 0
          cc   = -a3*c2
    !
          if (bb /= 0.d0) then
             dd = cc/(bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d3
          end if
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             psi6 = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             psi6 = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
          end if 
    !
    !  ## 2. Speciation 
          psi6  = min(psi6, c2)
          h     = c1 + psi6
          so4_t = psi6
          hso4  = max(c2 - psi6, 0.0d0)
    
    !  ## 3. Aerosol liquid water content
          lwn = max(c3 + max((so4_t + hso4 - c5), 0.0d0)/c6, tiny) 
       end do
    !
    !
    !  ### MINOR SYSTEM #1: Cl-/HCl/NO3-/HNO3/H+ ###
    !  ## Calculate dissolution of HCl, HNO3 in the presence of (H,SO4)
    !  ## HCl, HNO3 are considered minor species 
       ispoly = .false.
    !
    !  ## 1. Special case: aerosol liquid water content (lwn) = 0.0d0
       if (lwn <= tiny) then
    !  Gaseous species
          ghcl  = max(cl  - cl_t,  0.0d0)
          ghno3 = max(no3 - no3_t, 0.0d0) 
          m1    = 0.0d0
          m2    = 0.0d0
          m3    = 0.0d0
    !
    !  ## 2. Special case: HCl = HNO3 = 0.0d0
       else if (cl <= tiny .and. no3 <= tiny) then
          m1  = 0.0d0
          m2  = 0.0d0
          m3  = 0.0d0
    !
    !  ## 3. Special case: HCl = 0.0d0
       else if (cl <= tiny) then
    !  Nitric acid in the liquid phase is assumed a minor species
    !  HNO3(g) <--> (H+) + (NO3-), using H+ from the sulfates
          ff = 0.0d0
          m1 = 0.0d0
          m2 = 0.0d0
          m3 = 0.0d0
          if (lwn > tiny) then 
             bb = h + khno3*r*t*(lwn/gama(10))**2.0
             cc = (khno3*r*t*(lwn/gama(10))**2.0)*no3
    !
             if (bb /= 0.d0) then
                dd = cc / (bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d+03
             end if       
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                ff = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
             else
                ff = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
             end if 
          end if 
    !
    !  Speciation  
          ghno3 = max(no3 - ff, 0.0d0)
          no3_t = ff
          h     = h + ff
    !
    !  ## 4. Special case: HNO3 = 0.0d0
       else if (no3 <= tiny) then
    !  Hydrochloric acid in the liquid phase is assumed a minor species
    !  HCl (g) <--> (H+) + (Cl-) using (H+) from the sulfates
          m1 = 0.0d0
          m2 = 0.0d0
          m3 = 0.0d0
          ff = 0.0d0
          if (lwn > tiny) then
             bb = h + khcl*r*t*(lwn/gama(11))**2.0
             cc = (khcl*r*t*(lwn/gama(11))**2.0)*cl
    !
             if (bb /= 0.d0) then
                dd = cc / (bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d+03
             end if       
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                ff = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
             else
                ff = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
             end if 
          end if 
    !
    !  Speciation        
          ghcl = max(cl - ff, 0.0d0) 
          cl_t = ff 
          h    = h + ff
    !
    !  ## 5. All else: not a special case
       else
          ispoly = .true.
          a3 = khno3*r*t*(lwn/gama(10))**2.0  !HNO3
          a4 = khcl*r*t*(lwn/gama(11))**2.0   !HCl
    !
    !  Calculate cubic equation coefficients
          m1 = (a3*no3 + a4*cl + (h + a4)*(a3 - a4))/(a3 - a4)
          m2 = ((h + a4)*a4*cl - a4*(a3 - a4)*cl)/(a3 - a4)
          m3 = -a4*a4*cl*cl/(a3 - a4)
       end if 
    !
    !  ## Calculate roots 
    !   call mach_hetp_poly3v(m1, m2, m3, delcl, islv)     # Original ISORROPIA subroutine  
       call mach_hetp_poly(m1, m2, m3, cl, delcl, islv)
    !   
       delno = 0.0d0
       if (ispoly) then
          if (islv /= 0) then
    !  Tiny amounts of HCl are assumed when there is no root
             delcl = tiny                                          !Change in Cl-
          end if
    !
          a3    = khno3*r*t*(lwn/gama(10))**2.0                    !HNO3
          a4    = khcl*r*t*(lwn/gama(11))**2.0                     !HCl
          delno = min(a3*no3*delcl/(a4*cl + (a3 - a4)*delcl), no3) !Change in NO3-
    !
          if (delcl < 0.0d0 .or. delno < 0.0d0 .or. delcl > cl .or. delno > no3) then
             delcl = tiny                                          !Change in Cl-
             delno = tiny                                          !Change in NO3-
          end if 
    !
    !  ## Effect on liquid phase
          h     = h + delno + delcl       ! H+   change
          cl_t  = cl_t  + delcl           ! Cl-  change
          no3_t = no3_t + delno           ! NO3- change
    !
    !  ## Gaseous species
          ghcl  = max(cl  - cl_t,  0.0d0)
          ghno3 = max(no3 - no3_t, 0.0d0) 
       end if 
    !
    !
    !  ### MINOR SYSTEM #2: NH4+/NH3/H+ ###
    !  Ammonia in the gas phase is assumed a minor species that does not significantly perturb 
    !  the aerosol equilibrium: NH3(g) + H+(aq) <==> NH4+(aq); H+ determined from the major
    !  system (solved above) is used initially.  
       if (lwn > tiny) then
    !  ## Calculate NH3 sublimation
          ff = 0.0d0
          a3 = (knh3/kh2o)*r*t*(gama(10)/gama(5))**2.0
          bb = h + 1.0d0/a3    ! bb always > 0
          cc = -nh4_t/a3
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d+03
          end if  
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             ff = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             ff = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
          end if 
    !
    !  ## Speciation 
    !  ## Due to round-off, ff may be .gt. nh4_t giving negative nh4_t
    !  ## If this condition is true, then set ff = nh4_t
          ff    = max(tiny, min(ff, nh4_t))
          gnh3  = ff
          nh4_t = max(nh4_t - ff, 0.0d0)
          h     = h + ff
       end if
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t 
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       na_i    = na_t
       cl_i    = cl_t
       nh3g_i  = gnh3
       hno3g_i = ghno3
       hclg_i  = ghcl
       h_i     = h
       lwn_i   = lwn
    !
       return
    end subroutine mach_hetp_calcj3
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: O7; Sulfate poor; dust and sodium poor
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calco7(so4_i, nh4_i, nh3g_i, hno3g_i, hclg_i, hso4_i,          &
                                na_i, cl_i, no3_i, h_i, lwn_i, ca_i, k_i, mg_i,         &
                                caso4_i, frmg, frna, frca, frk, rh, temp, k0,           &
                                p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: na_i    
       real(kind=8),    intent(inout) :: cl_i    
       real(kind=8),    intent(inout) :: ca_i    
       real(kind=8),    intent(inout) :: k_i     
       real(kind=8),    intent(inout) :: mg_i    
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: hclg_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent(inout) :: caso4_i 
       real(kind=8),    intent(inout) :: frna
       real(kind=8),    intent(inout) :: frmg    
       real(kind=8),    intent(inout) :: frk   
       real(kind=8),    intent(inout) :: frca 
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    !
    ! ## Local variables
       real(kind=8)     :: so4, nh4, hso4, gnh3, h, lwn, caso4, no3, cl, na, ca, pk, mg, ghno3, ghcl
       real(kind=8)     :: t, aw, knh3, khso4, kh2o, khno3, khcl, bb, cc, dd, hh, v, errin
       real(kind=8)     :: smin, ohi, frnh4, a4, a5, psi4, psi5, m4, m5, m6, ak1, scon
       real(kind=8)     :: tt0, errinlocb, errouloc, loccon
       real(kind=8)     :: omehi, omebe, y1, y2, y3, x3, dx, c1, c2, c2a, c3, gmax, ya, yb, xa, xb
       real(kind=8)     :: so4_t, nh4_t, no3_t, na_t, cl_t, ca_t, pk_t, mg_t, k1, k2, k3, k4
       real(kind=8)     :: so4fr, a6, na2so4, k2so4, mgso4, tt1, tt2, c6a, c7, c8, c9, c10
       real(kind=8)     :: nh, sigma, xt, xf, xh, delta, rr, gx, gx2, u1, lwnsq, gama10sq
       integer(kind=4)  :: j, k, ii, rooteval, irh, nmax
       logical(kind=4)  :: condition, noroot, soln, frst, calain, calou, earlye
       real(kind=8), dimension(23) :: gama, gamin, gamou, gama_min
       real(kind=8)      :: y3_min, x3_min, y3_lastiter, x3_lastiter
       real(kind=8)      :: no3_min, nh4_min, cl_min, h_min, lwn_min, ghcl_min, gnh3_min, ghno3_min
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       na    = na_i
       cl    = cl_i
       ca    = ca_i
       pk    = k_i
       mg    = mg_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       gnh3  = 0.0d0
       ghno3 = 0.0d0
       ghcl  = 0.0d0
       h     = 0.0d0
       lwn   = tiny
       so4_t = 0.0d0
       nh4_t = 0.0d0
       no3_t = 0.0d0
       na_t  = 0.0d0
       cl_t  = 0.0d0
       ca_t  = 0.0d0
       pk_t  = 0.0d0
       mg_t  = 0.0d0
       caso4 = 0.0d0
       so4fr = 0.0d0
       na2so4= 0.0d0
       k2so4 = 0.0d0
       mgso4 = 0.0d0
       noroot=.false.
       frk   = 0.0d0
       frmg  = 0.0d0
       frca  = 0.0d0
       frna  = 0.0d0
       soln  = .false.
       calou = .true.
       gama  = 0.1d0
       gamin = 1.0d10
       gamou = 0.1d0
       earlye = .false.
    ! 
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0    = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  Net NH3: k2*k3                                                   (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    !
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3= k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. HCl(g) <==> H+(aq) + Cl-(aq)                               (xk3)
       khcl = k0(11) * exp(p1(11)*tt1 + p2(11)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Constants
       c1 = r*t
       k1 = (knh3/kh2o)*c1
       k2 = khno3*c1
       k3 = khcl*c1
       k4 = kh2o*aw
    !
    !  ## Calculate dry salt composition 
    !  ## Salts assumed to have completely dissolved (MgSO4, Na2SO4, K2SO4)
       caso4 = min(ca, so4)                           
       so4fr = max(so4 - caso4, 0.0d0)                
       frca  = max(ca  - caso4, 0.0d0)                
       k2so4 = min(0.5d0*pk, so4fr)                   
       frk   = max(pk  - 2.0d0*k2so4, 0.0d0)         
       so4fr = max(so4fr  - k2so4, 0.0d0)             
       na2so4= min(0.5d0*na, so4fr)                   
       frna  = max(na  - 2.0d0*na2so4, 0.0d0)       
       so4fr = max(so4fr  - na2so4, 0.0d0)         
       mgso4 = min(mg, so4fr)                         
       frmg  = max(mg  - mgso4, 0.0d0)               
       so4fr = max(so4fr  - mgso4, 0.0d0)          
    !
    !  ## Initial speciation 
       na_t  = 2.0d0*na2so4
       so4_t = na2so4 +  max(so4fr, 0.0d0) + k2so4 + mgso4        
       pk_t  = 2.0d0*k2so4
       mg_t  = mgso4
    !
    !  ## Constants
       c2  = 0.5d0*na_t - 0.5d0*pk_t - mg_t
       c2a = 2.0d0*max(so4fr, 0.0d0)
       c3  = max(nh4 - c2a,0.0d0)             ! free nh4 == nh3 dry
    !
    !  ## Save ZSR arrays as variables to limit indirect addressing 
       c6a = na2so4/awss(irh) + k2so4/awps(irh) + mgso4/awms(irh)
       c7  = awas(irh)
       c8  = awan(irh)
       c9  = awac(irh)
       c10 = na2so4 + k2so4 + mgso4
    !
    !  ## Initial aerosol liquid water content from dry salt composition 
       lwn = max(max(so4fr, 0.0d0)/c7 + c6a, tiny)
    !
    !
    !  ### STAGE 1: Root tracking ###
    !  ## Find a subinterval [xa,xb] on the larger interval [I1,I2] where a sign change occurs
       rooteval = 0
       condition = .true.
       do while (rooteval < 2 .or. (condition .and. rooteval < ndiv + 1))   ! Begin outer loop for root tracking
          rooteval = rooteval + 1
    !
    ! ## Set high limit for root tracking (i.e. lower bound)
          if (rooteval == 1) then
            omehi = tiny
            y1    = 1.0d0
          end if
    !
    !  ## Begin search on subinterval 
          if (rooteval == 2) then
             soln = .false.
             y1    = y2
    !
             dx = (cl-tiny-tiny)/float(ndiv)
    !
             omebe = omehi            ! Lower bound of subinterval (xa)
             omehi = omehi + dx       ! Upper bound of subinterval (xb)
          end if
    !
    !  ## Continue search 
          if (rooteval > 2) then
             if (sign(1.0d0,y1)*sign(1.0d0,y2) < 0.0d0) then
    !  ## 1. Root has been found on the subinterval; save x values for ITP search
                y1    = y1
                omebe = omebe
                omehi = omehi
             else
    !  ## 2. No root has been found, continue searching in the next subinterval
                y1    = y2
                omebe = omehi
                omehi = omehi + dx
             end if
          end if
    !
    !  ## Solve the system of equations 
          frst   = .true.
          calain = .true.
    !
          if (( .not. soln)) then
             lwnsq    = lwn*lwn
             gama10sq = gama(10)*gama(10)
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))      
             a5 = k2*lwnsq/gama10sq                  
             a6 = k3*lwnsq/(gama(11)*gama(11))       
    !
             if (no3 >= tiny) then
                psi5 = min(omehi*no3/(a6/a5*(cl - omehi) + omehi), no3)
             else
                psi5 = tiny
             end if
    !
    !  ## 1. Account for NH3 evaporation
             if (so4 > tiny) then
                bb = -(c3 + omehi + psi5 + 1.0d0/a4)
                cc = c3*(psi5 + omehi) - c2a/a4
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
    
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                else
                   psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                psi4 = max(min(psi4, c3) , 0.0d0)
             else
                psi4 = tiny
             end if
    !
    !  ## 2. Speciation
             nh4_t = c2a + psi4
             cl_t  = omehi
             no3_t = psi5
             gnh3  = max(c3  - psi4,  tiny2)
             ghno3 = max(no3 - psi5,  tiny2)
             ghcl  = max(cl  - omehi, tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t - pk_t - 2.0d0*mg_t
             scon = k4*lwnsq     
    !
             if (smin > tiny) then                        ! H+ in excess
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                h   = max(hh,sqrt(scon))
    !           ohi = scon/h
             else                                        ! OH- in excess
                bb = smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                ohi = max(hh,sqrt(scon))
                h   = scon/ohi
             end if 
    !
    !  ## 4. Aerosol liquid water content
             m4     = max(so4_t - c10, 0.0d0)
             frnh4  = max(nh4_t - 2.0d0*m4, 0.0d0)
             m5     = min(no3_t, frnh4)
             frnh4  = max(frnh4 - m5, 0.0d0)
             m6     = min(cl_t, frnh4)
             lwn    = max(c6a + m4/c7 + m5/c8 + m6/c9, tiny)
          end if 
    !
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)
             k = k + 1
    !
    !  ## Reset gamou
             if ((.not. soln) .and. (frst)) then
                gamou(1)   = gama(1)
                gamou(2)   = gama(2)
                gamou(3)   = gama(3)
                gamou(4)   = gama(4)
                gamou(5)   = gama(5)
                gamou(6)   = gama(6)
                gamou(7)   = gama(7)
                gamou(8)   = gama(8)
                gamou(9)   = gama(9)
                gamou(10)  = gama(10)
                gamou(11)  = gama(11)
                gamou(12)  = gama(12)
                gamou(13)  = gama(13)
                gamou(14)  = gama(14)
                gamou(15)  = gama(15)
                gamou(16)  = gama(16)
                gamou(17)  = gama(17)
                gamou(18)  = gama(18)
                gamou(19)  = gama(19)
                gamou(20)  = gama(20)
                gamou(21)  = gama(21)
                gamou(22)  = gama(22)
                gamou(23)  = gama(23)
             end if 
    !
    !  ## Reset gamin
             if (( .not. soln)) then
                gamin(1)   = gama(1)
                gamin(2)   = gama(2)
                gamin(3)   = gama(3)
                gamin(4)   = gama(4)
                gamin(5)   = gama(5)
                gamin(6)   = gama(6)
                gamin(7)   = gama(7)
                gamin(8)   = gama(8)
                gamin(9)   = gama(9)
                gamin(10)  = gama(10)
                gamin(11)  = gama(11)
                gamin(12)  = gama(12)
                gamin(13)  = gama(13)
                gamin(14)  = gama(14)
                gamin(15)  = gama(15)
                gamin(16)  = gama(16)
                gamin(17)  = gama(17)
                gamin(18)  = gama(18)
                gamin(19)  = gama(19)
                gamin(20)  = gama(20)
                gamin(21)  = gama(21)
                gamin(22)  = gama(22)
                gamin(23)  = gama(23)
             end if 
    !
             call mach_hetp_calcact4b(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                      ca_t, pk_t, mg_t, lwn, gama, t, soln, frst,   &
                                      calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 23
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    
             errinlocb = 0
             do ii = 1, 23
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients
             do ii = 1, 23
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do
    !
    !  ## Solve system of equations, using new activity coefficients 
             if (( .not. soln)) then
                lwnsq    = lwn*lwn
                gama10sq = gama(10)*gama(10)
    !
                a4 = k1*gama10sq/(gama(5)*gama(5))     
                a5 = k2*lwnsq/gama10sq                 
                a6 = k3*lwnsq/(gama(11)*gama(11))      
    !
                if (no3 >= tiny) then
                   psi5 = min(omehi*no3/(a6/a5*(cl - omehi) + omehi), no3)
                else
                   psi5 = tiny
                end if
    !
    !  ## 1. Account for NH3 evaporation
                if (so4 > tiny) then
                   bb = -(c3 + omehi + psi5 + 1.0d0/a4)
                   cc = c3*(psi5 + omehi) - c2a/a4
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                   else
                      psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   psi4 = max(min(psi4, c3) , 0.0d0)
                else
                   psi4 = tiny
                end if
    !
    !  ## 2. Speciation
                nh4_t = c2a + psi4
                cl_t  = omehi
                no3_t = psi5
                gnh3  = max(c3  - psi4,  tiny2)
                ghno3 = max(no3 - psi5,  tiny2)
                ghcl  = max(cl  - omehi, tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t - pk_t - 2.0d0*mg_t
                scon = k4*lwnsq   
    !
                if (smin > tiny) then                        ! H+ in excess
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   h   = max(hh,sqrt(scon))
    !              ohi = scon/h
                else                                        ! OH- in excess
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if 
    !
    !  ## 4. Aerosol liquid water content
                m4     = max(so4_t - c10, 0.0d0)
                frnh4  = max(nh4_t - 2.0d0*m4, 0.0d0)
                m5     = min(no3_t, frnh4)
                frnh4  = max(frnh4 - m5, 0.0d0)
                m6     = min(cl_t, frnh4)
                lwn    = max(c6a + m4/c7 + m5/c8 + m6/c9, tiny)
             end if     
          end do
    !
          y2 = h*cl_t/ghcl/a6 - 1.0d0  !Function value
    !
    !
    !  ## Check for criteria to exit root tracking 
          condition = .false.
          loccon = sign(1.0d0,y1)*sign(1.0d0,y2)
          if ((loccon > 0.0d0) .and. (( .not. noroot)) .and. (abs(y2) > eps) .and. (cl > tiny)) then
             condition = .true.         
          elseif (loccon < 0.0d0 .and. abs(y2) > eps) then
    !  ## Interval had been found where sign change occurs; exit root tracking and proceed to ITP
             soln = .true. 
          else
    !  ## abs(y2) <= eps; solution is assumed; exit root tracking and proceed to minor system (no ITP)
             if (rooteval >= 2) then
                soln   = .true.
                noroot = .true.
                earlye = .true.
             else
                condition = .true.
             end if
          end if 
    !
    !  ## Too little frcl; reset x-value to tiny and exit root tracking 
          if (cl <= tiny) then
             condition = .false.
             rooteval  = 2       ! Increment rooteval by 1 to force exit from root tracking
             noroot    = .true.
             omehi     = tiny    ! Reset x-value to tiny, and solve system with this value
             omebe     = omehi
          end if 
    !
    !  ### AFTER iterating through ALL ndiv subdivided intervals
          if (rooteval == ndiv + 1) then
             loccon = sign(1.0d0,y1)*sign(1.0d0,y2)
             if (loccon > 0.0d0 .and. abs(y2) > eps) then
    !  ## (1) No solution
                noroot = .true.
                omehi  = tiny    ! Reset to tiny
                omebe  = omehi
    !           write(*,*), 'Warning in CALCO7: no solution found'
             else if (loccon > 0.0d0 .and. abs(y2) <= eps) then
    !  ## (2) Solution is assumed and ITP is not required
                noroot = .true.
             end if
          end if
       end do  !End outer loop for root tracking
    !
    !
    !
    !  ### STAGE 2: modified bisection search (using ITP algorithm) ###
    !  ## Initialize static ITP variables 
       if (( .not. noroot)) then
          ya = y1
          yb = y2
          xa = omebe
          xb = omehi
          x3 = omehi
    !
          if (xa == xb) then
             noroot = .true.
             gx = tiny
          else
             gx = xb - xa
          end if 
    !
          gx2  = (xa+xb)*0.5d0
          u1   = 0.2d0/gx
          nh   = log10(abs(gx/(2.0d0*eps*gx2))) / log10of2
          nmax = int(nh) + 2
       else
          x3 = omehi
       end if 
    !
    !  ## Start search
       y3_lastiter = 0.0d0
       x3_lastiter = 0.0d0
       y3_min = 1.0d20
       x3_min = 0.0d0
    !
       j = 0
       condition = .true.
    !
       if (earlye) then
          soln = .true.
       else
          soln = .false.
       end if
    !
       do while (j < maxit .and. condition)   ! Begin outer loop for ITP search
    !  1. Track x3 and y3 of the previous iteration
          if (j > 0) then
             y3_lastiter = y3
             x3_lastiter = x3
          end if
    !
    !  2. Track the minimum y3 that is found before ending
          if (abs(0.0d0 - y3_min) > abs(0.0d0 - y3_lastiter) .and. j > 0) then
             y3_min    = y3_lastiter
             x3_min    = x3_lastiter
             ghno3_min = ghno3
             no3_min   = no3_t
             h_min     = h
             ghcl_min  = ghcl
             cl_min    = cl_t
             gnh3_min  = gnh3
             nh4_min   = nh4_t
             lwn_min   = lwn
             gama_min  = gama
          end if
    !
          if ((.not. noroot) .and. (.not. soln)) then
    !  ## Set dynamic ITP variables 
             if (yb - ya == 0.0d0) then
                write(*,*), '######        ABORT       ######'
                write(*,*), 'Zero divide in ITP reset: CALCO7'
                write(*,*), 'SO4 in = ', so4
                write(*,*), 'NH4 in = ', nh4
                write(*,*), 'NO3 in = ', no3
                write(*,*), 'Na in  = ', na
                write(*,*), 'Cl in  = ', cl
                write(*,*), 'Ca in  = ', ca
                write(*,*), 'Mg in  = ', mg
                write(*,*), 'K in   = ', pk
                write(*,*), 'Temp in= ', t
                write(*,*), 'RH in  = ', aw
                return
             end if  
    !
             gx   = xb - xa
             xh   = 0.5d0*(xa + xb)
             rr   = max(gx2*eps*2.d0**(real(nmax - j)) - 0.5d0*gx, 0.0d0)
             delta= u1*(max(gx, 0.0d0))**2.0d0   
             xf   = max((yb*xa - ya*xb) / (yb - ya), 0.0d0)
    !
             sigma = sign(1.0d0, xh - xf)
             if (delta <= abs(xh - xf)) then 
                xt = xf + sigma*delta             
             else
                xt = xh
             end if
    !
             if (abs(xt - xh) <= rr) then
                x3 = xt
             else
                x3 = xh - sigma*rr
             end if 
          end if 
    !
          j = j + 1
    !
          if (( .not. soln)) then
             gmax = 0.1d0
             gmax = max(gmax, gama(1))
             gmax = max(gmax, gama(2))
             gmax = max(gmax, gama(3))
             gmax = max(gmax, gama(4))
             gmax = max(gmax, gama(5))
             gmax = max(gmax, gama(6))
             gmax = max(gmax, gama(7))
             gmax = max(gmax, gama(8))
             gmax = max(gmax, gama(9))
             gmax = max(gmax, gama(10))
             gmax = max(gmax, gama(11))
             gmax = max(gmax, gama(12))
             gmax = max(gmax, gama(13))
             gmax = max(gmax, gama(14))
             gmax = max(gmax, gama(15))
             gmax = max(gmax, gama(16))
             gmax = max(gmax, gama(17))
             gmax = max(gmax, gama(18))
             gmax = max(gmax, gama(19))
             gmax = max(gmax, gama(20))
             gmax = max(gmax, gama(21))
             gmax = max(gmax, gama(22))
             gmax = max(gmax, gama(23))
          end if
    !
    !  ## Reinitialize activity coefficients if gmax > 100.0d0
          if (gmax > 100.0d0 .and. ( .not. soln)) then
             gama  = 0.1d0
             gamin = 1.0d10
             gamou = 1.0d10
             calou = .true.
             frst  = .true.
          end if
    !
          frst   = .true.
          calain = .true.
          if (( .not. soln)) then
             lwnsq    = lwn*lwn
             gama10sq = gama(10)*gama(10)
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))   
             a5 = k2*lwnsq/gama10sq               
             a6 = k3*lwnsq/(gama(11)*gama(11))    
    !
             if (no3 >= tiny) then
                psi5 = min(x3*no3/(a6/a5*(cl - x3) + x3), no3)
             else
                psi5 = tiny
             end if
    !
    !  ## 1. Account for NH3 evaporation
             if (so4 > tiny) then
                bb = -(c3 + x3 + psi5 + 1.0d0/a4)
                cc = c3*(psi5 + x3) - c2a/a4
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                else
                   psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                psi4 = max(min(psi4, c3) , 0.0d0)
             else
                psi4 = tiny
             end if
    !
    !  ## 2. Speciation
             nh4_t = c2a + psi4
             cl_t  = x3
             no3_t = psi5
             gnh3  = max(c3  - psi4, tiny2)
             ghno3 = max(no3 - psi5, tiny2)
             ghcl  = max(cl  - x3,   tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t - pk_t - 2.0d0*mg_t
             scon = k4*lwnsq     
    !
             if (smin > tiny) then                        ! H+ in excess
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                h   = max(hh,sqrt(scon))
    !           ohi = scon/h
             else                                        ! OH- in excess
                bb = smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                ohi = max(hh,sqrt(scon))
                h   = scon/ohi
             end if
    !
    !  ## 4. Aerosol liquid water content
             m4     = max(so4_t - c10, 0.0d0)
             frnh4  = max(nh4_t - 2.0d0*m4, 0.0d0)
             m5     = min(no3_t, frnh4)
             frnh4  = max(frnh4 - m5, 0.0d0)
             m6     = min(cl_t, frnh4)
             lwn    = max(c6a + m4/c7 + m5/c8 + m6/c9, tiny)
          end if
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)
             k = k + 1
    !
    !  ## Reset gamin and gamou
             if (( .not. soln) .and. frst) then
                gamou(1)   = gama(1)
                gamou(2)   = gama(2)
                gamou(3)   = gama(3)
                gamou(4)   = gama(4)
                gamou(5)   = gama(5)
                gamou(6)   = gama(6)
                gamou(7)   = gama(7)
                gamou(8)   = gama(8)
                gamou(9)   = gama(9)
                gamou(10)  = gama(10)
                gamou(11)  = gama(11)
                gamou(12)  = gama(12)
                gamou(13)  = gama(13)
                gamou(14)  = gama(14)
                gamou(15)  = gama(15)
                gamou(16)  = gama(16)
                gamou(17)  = gama(17)
                gamou(18)  = gama(18)
                gamou(19)  = gama(19)
                gamou(20)  = gama(20)
                gamou(21)  = gama(21)
                gamou(22)  = gama(22)
                gamou(23)  = gama(23)
             end if 
    !
             if (( .not. soln)) then
                gamin(1)   = gama(1)
                gamin(2)   = gama(2)
                gamin(3)   = gama(3)
                gamin(4)   = gama(4)
                gamin(5)   = gama(5)
                gamin(6)   = gama(6)
                gamin(7)   = gama(7)
                gamin(8)   = gama(8)
                gamin(9)   = gama(9)
                gamin(10)  = gama(10)
                gamin(11)  = gama(11)
                gamin(12)  = gama(12)
                gamin(13)  = gama(13)
                gamin(14)  = gama(14)
                gamin(15)  = gama(15)
                gamin(16)  = gama(16)
                gamin(17)  = gama(17)
                gamin(18)  = gama(18)
                gamin(19)  = gama(19)
                gamin(20)  = gama(20)
                gamin(21)  = gama(21)
                gamin(22)  = gama(22)
                gamin(23)  = gama(23)
             end if 
    !
             call mach_hetp_calcact4b(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                      ca_t, pk_t, mg_t, lwn, gama, t, soln, frst,   &
                                      calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 23
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    
             errinlocb = 0
             do ii = 1, 23
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients
             do ii = 1, 23
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do
    !
    !  ## Solve system of equations, with new activity coefficients
             if (( .not. soln)) then
                lwnsq    = lwn*lwn
                gama10sq = gama(10)*gama(10)
    !
                a4 = k1*gama10sq/(gama(5)*gama(5))     
                a5 = k2*lwnsq/gama10sq                 
                a6 = k3*lwnsq/(gama(11)*gama(11))    
    !
                if (no3 >= tiny) then
                   psi5 = min(x3*no3/(a6/a5*(cl - x3) + x3), no3)
                else
                   psi5 = tiny
                end if
    !
    !  ## 1. Account for NH3 evaporation
                if (so4 > tiny) then
                   bb = -(c3 + x3 + psi5 + 1.0d0/a4)
                   cc = c3*(psi5 + x3) - c2a/a4
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      psi4 = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                   else
                      psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   psi4 = max(min(psi4, c3) , 0.0d0)
                else
                   psi4 = tiny
                end if
    !
    !  ## 2. Speciation
                nh4_t = c2a + psi4
                cl_t  = x3
                no3_t = psi5
                gnh3  = max(c3  - psi4, tiny2)
                ghno3 = max(no3 - psi5, tiny2)
                ghcl  = max(cl  - x3,   tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t - pk_t - 2.0d0*mg_t
                scon = k4*lwnsq       
    !
                if (smin > tiny) then                        ! H+ in excess
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   h = max(hh,sqrt(scon))
    !               ohi  = scon/h
                else                                        ! OH- in excess
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if
    !
    !  ## 4. Aerosol liquid water content
                m4     = max(so4_t - c10, 0.0d0)
                frnh4  = max(nh4_t - 2.0d0*m4, 0.0d0)
                m5     = min(no3_t, frnh4)
                frnh4  = max(frnh4 - m5, 0.0d0)
                m6     = min(cl_t, frnh4)
                lwn    = max(c6a + m4/c7 + m5/c8 + m6/c9, tiny)
             end if      
          end do
    !
          y3 = h*cl_t/ghcl/a6 - 1.0d0    !Function value
    !
          condition = .false.
          if (noroot) then
    !  ## If no root on interval then do not perform ITP
             xa = x3
             xb = x3
          else if (y3 > 0.0d0 .and. ( .not. soln)) then
             xb = x3
             yb = y3
          else if (y3 < 0.0d0 .and. ( .not. soln)) then
             xa = x3
             ya = y3
          else if (( .not. soln)) then
             xa = x3
             xb = x3
          end if
    !
    !  ## Test for convergence criteria to exit ITP:
          if (xb - xa > abs(xa*eps) .and. cl > tiny .and. ( .not. noroot)) then
             condition = .true.
             soln   = .false.
          else
             soln   = .true.
          end if
    !
    ! ## Exit ITP if the function being bisected evaluates to a value <= eps; solution is assumed
          if (abs(y3) <= eps .and. ( .not. noroot)) then
             soln = .true.
             condition = .false.
          end if
    !
    ! ## Post-convergence correction: 
    ! ## If the calculated y-value (y3) after ITP is further from zero than an earlier iteration 
    ! ## then reset to the x-value/concentrations/activity coefficients that were found to minimize
    ! ## the objective function (i.e., y3); in this case, this is chosen as the solution
          if (( .not. condition) .and. ( .not. noroot) .and. abs(y3) > 0.1d0) then     
             if (abs(0.0d0 - y3_min) < abs(0.0d0 - y3) .and. abs(y3_min - y3) > 1.0d-1) then
                x3    = x3_min
                cl_t  = cl_min
                nh4_t = nh4_min
                no3_t = no3_min
                h     = h_min
                lwn   = lwn_min
                ghcl  = ghcl_min
                ghno3 = ghno3_min
                gnh3  = gnh3_min
    ! 
    !  ## Reset activity coefficients 
                gama  = gama_min
                a6    = k3*(lwn/gama(11))*(lwn/gama(11))
                y3    = h*cl_t/ghcl/a6 - 1.0d0
             end if
          end if
       end do ! End outer loop of ITP search
    !
    !
    !  ### MINOR SYSTEM: HSO4-/SO42-/H+ ###
       hso4 = 0.0d0
       if (h > tiny .and. so4_t > tiny .and. lwn >= 1.0d-19) then
          ak1 = khso4*lwn/gama(7)*(gama(8)/gama(7))*(gama(8)/gama(7))
          bb =-(h + so4_t + ak1)
          cc = h*so4_t
          dd = bb*bb - 4.d0*cc
    !
          if (dd >= 0.0d0) then
             if (bb /= 0.d0) then
                dd = cc/(bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d3
             end if
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                hh = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
             else
                hh = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
             end if
    !
             hh    = max(tiny, min(hh, min(h, so4_t)))     ! To avoid negative H+   (i.e., if hh > h or hh > so4_t)
             h     = max(h     - hh, 0.0d0)
             so4_t = max(so4_t - hh, 0.0d0)
             hso4  = hh
          end if
       end if
    !
    !
    !  ### Perform mass adjustment to machine precision, if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
          so4_i   = so4_t
          nh4_i   = nh4_t
          no3_i   = no3_t
          hso4_i  = hso4
          na_i    = na_t
          cl_i    = cl_t
          ca_i    = ca_t
          k_i     = pk_t
          mg_i    = mg_t
          nh3g_i  = gnh3
          hno3g_i = ghno3
          hclg_i  = ghcl
          h_i     = h
          lwn_i   = lwn
          caso4_i = caso4
    !
       return
    end subroutine mach_hetp_calco7
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: M8; Sulfate poor; dust and sodium rich
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcm8(so4_i, nh4_i, nh3g_i, hno3g_i, hclg_i, hso4_i,          &
                                na_i, cl_i, no3_i, h_i, lwn_i, ca_i, k_i, mg_i,         &
                                caso4_i, frca, frmg, frk, frna, rh, temp,        &
                                k0, p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: na_i    
       real(kind=8),    intent(inout) :: cl_i    
       real(kind=8),    intent(inout) :: ca_i    
       real(kind=8),    intent(inout) :: k_i     
       real(kind=8),    intent(inout) :: mg_i    
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: hclg_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent(inout) :: caso4_i 
       real(kind=8),    intent(inout) :: frna 
       real(kind=8),    intent(inout) :: frmg
       real(kind=8),    intent(inout) :: frk 
       real(kind=8),    intent(inout) :: frca
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    ! ## Local variables
       real(kind=8)     :: so4, nh4, hso4, gnh3, h, lwn, caso4, no3, cl, na, ca, pk, mg, ghno3, ghcl
       real(kind=8)     :: t, aw, knh3, khso4, kh2o, khno3, khcl, bb, cc, dd, hh, v, errin, smin, ohi
       real(kind=8)     :: a4, a5, psi4, psi5, m5, m6, ak1, scon, tt0
       real(kind=8)     :: errouloc, errinlocb, c1, gmax, ztot, loccon
       real(kind=8)     :: omehi, omebe, y1, y2, y3, x3, dx, a6, tt1, tt2, k1, k2, ya, yb, xa, xb
       real(kind=8)     :: so4_t, nh4_t, no3_t, na_t, cl_t, ca_t, pk_t, mg_t, k3, k4
       real(kind=8)     :: na2so4, chi4, chi5, chi6, nacl, nano3, k2so4, mgso4, c5, c6
       real(kind=8)     :: frcl_d, frno3_d, frnh4_d, frno3, frnh4, frcl, frso4, lwnsq, gama10sq
       real(kind=8)     :: nh, sigma, xt, xf, xh, delta, rr, gx, gx2, u1
       integer(kind=4)  :: j, k, ii, rooteval, irh, nmax
       logical(kind=4)  :: condition, noroot, earlyexit, soln, frst, calain, calou, earlye
       real(kind=8), dimension(23) :: gama, gamin, gamou, gama_min
       real(kind=8)      :: y3_min, x3_min, y3_lastiter, x3_lastiter
       real(kind=8)      :: no3_min, nh4_min, cl_min, h_min, lwn_min, ghcl_min, gnh3_min, ghno3_min
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       na    = na_i
       cl    = cl_i
       ca    = ca_i
       pk    = k_i
       mg    = mg_i
       aw    = rh
       t     = temp
       na2so4= 0.0d0
       chi4  = 0.0d0
       chi5  = 0.0d0
       chi6  = 0.0d0
       nacl  = 0.0d0
       nano3 = 0.0d0
       k2so4 = 0.0d0
       mgso4 = 0.0d0
       hso4  = 0.0d0
       gnh3  = 0.0d0
       ghno3 = 0.0d0
       ghcl  = 0.0d0
       h     = 0.0d0
       lwn   = tiny
       so4_t = 0.0d0
       nh4_t = 0.0d0
       no3_t = 0.0d0
       na_t  = 0.0d0
       cl_t  = 0.0d0
       ca_t  = 0.0d0
       pk_t  = 0.0d0
       mg_t  = 0.0d0
       caso4 = 0.0d0
       frmg  = 0.0d0
       frk   = 0.0d0
       frca  = 0.0d0
       frso4 = 0.0d0
       frna  = 0.0d0
       frnh4 = 0.0d0
       frcl  = 0.0d0
       frno3 = 0.0d0
       frno3_d=0.0d0
       frnh4_d=0.0d0
       frcl_d =0.0d0
       earlyexit = .false.
       noroot =.false.
       soln   = .false.
       calou  = .true.
       gmax   = 0.0d0
       gama   = 0.1d0
       gamou  = 0.1d0
       gamin  = 1.0d10
       earlye = .false.
    ! 
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    ! 
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !   
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3= k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. HCl(g) <==> H+(aq) + Cl-(aq)                               (xk3)
       khcl = k0(11) * exp(p1(11)*tt1 + p2(11)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Constants
       c1 = r*t
       k1 = (knh3/kh2o)*c1
       k2 = khno3*c1
       k3 = khcl*c1
       k4 = kh2o*aw                
    !
    !  ## Calculate dry salt composition; keep track of 'free' amounts
    !  ## Salts assumed to have completely dissolved (NaNO3, NaCl, MgSO4, Na2SO4, K2SO4)
       caso4   = min(ca, so4)                               
       frso4   = max(so4 - caso4, 0.0d0)                          
       frca    = max(ca - caso4, 0.0d0)             
       k2so4   = min(0.5d0*pk, frso4)                         
       frk     = max(pk - 2.0d0*k2so4, 0.0d0)   
       frso4   = max(frso4 - k2so4, 0.0d0) 
       mgso4   = min(mg, frso4)                              
       frmg    = max(mg - mgso4, 0.0d0)           
       frso4   = max(frso4 - mgso4, 0.0d0)        
       na2so4  = max(frso4, 0.0d0)                               
       frna    = max(na - 2.0d0*na2so4, 0.0d0) 
       frso4   = max(frso4 - na2so4, 0.0d0)   
       nano3   = min(frna, no3)
       frna    = max(frna - nano3, 0.0d0)     
       frno3_d = max(no3 - nano3, 0.0d0)                          
       chi4    = nh4                          
    !   frnh4_d = max(nh4 - chi4, 0.0d0)                                      
       chi5    = max(no3 - nano3, 0.0d0)     
    !   frno3_d = max(frno3_d - chi5, 0.0d0)   
       nacl    = min(frna, cl)  
    !   frcl_d  = max(cl - nacl, 0.0d0)       
       frna    = max(frna - nacl, 0.0d0)      
       chi6    = max(cl - nacl, 0.0d0)        
    !   frcl_d  = max(frcl_d - chi6, 0.0d0)              
    !
    !  ## Initial aqueous speciation 
       na_t  = nano3 + nacl + 2.0d0*na2so4
       so4_t = na2so4 + k2so4 + mgso4
       pk_t  = 2.0d0*k2so4
       mg_t  = mgso4
    !
    !  ## Constant ZSR 
       ztot = nacl/awsc(irh) + na2so4/awss(irh) + nano3/awsn(irh) + &
              k2so4/awps(irh) + mgso4/awms(irh)
       c5   = awan(irh)
       c6   = awac(irh)
    !
    !
    !  ### STAGE 1: Root tracking ###
    !  ## Find a subinterval [xa,xb] on the larger interval [I1,I2] where a sign change occurs
       rooteval = 0
       condition = .true.
       do while (rooteval < 2 .or. (condition .and. rooteval < ndiv + 1))   ! Begin outer loop for root tracking
          rooteval = rooteval + 1
    !
    !  ## Set high limit for root tracking (i.e. lower bound)
          if (rooteval == 1) then
             omehi = tiny   
             y1    = 1.0d0
          end if
    !      
    !  ## Begin search on subinterval 
          if (rooteval == 2) then
             soln = .false.
             if (abs(y2) <= eps) then
                earlyexit = .true.
                dx = 0.0d0
             else
                dx = (chi6-tiny-tiny)/float(ndiv)
             end if 
    !
             y1 = y2
    !
             omebe = omehi            ! Lower bound of subinterval (xa)
             omehi = omehi + dx       ! Upper bound of subinterval (xb)
          end if 
    !
    !  ## Continue search 
          if (rooteval > 2) then
             if (sign(1.0d0,y1)*sign(1.0d0,y2) < 0.0d0) then
    !  ## 1. Root has been found on the subinterval; save x values for ITP search
                y1    = y1
                omebe = omebe
                omehi = omehi
             else
    !  ## 2. No root has been found, continue searching in the next subinterval
                y1    = y2
                omebe = omehi
                omehi = omehi + dx
             end if 
          end if 
    !      
    !  ## Solve the system of equations 
          frst   = .true.
          calain = .true.
          if (( .not. soln)) then
             lwnsq = lwn*lwn
             gama10sq = gama(10)*gama(10)
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))    
             a5 = k2*lwnsq/gama10sq                
             a6 = k3*lwnsq/(gama(11)*gama(11))     
    !
    !  ## 1. Calculate dissociation quantities 
             psi5 = chi5*(omehi + nacl) - a6/a5*nano3*(chi6 - omehi)
             psi5 = psi5/(a6/a5*(chi6 - omehi) + omehi + nacl)
             psi5 = min(max(psi5, tiny), chi5)
    !
             if (nh4 > tiny .and. lwn > tiny) then
                bb   = -(nh4 + omehi + psi5 + 1.0d0/a4)
                cc   = nh4*(psi5 + omehi)
                psi4 = min(max(0.5d0*(-bb - sqrt(max(bb*bb - 4.0d0*cc, 0.0d0))), 0.0d0), nh4)
             else
                psi4 = tiny
             end if
    !
    !  ## 2. Speciation 
             nh4_t = psi4
             cl_t  = omehi + nacl
             no3_t = psi5 + nano3
             gnh3  = max(nh4  - psi4,  tiny2) 
             ghno3 = max(chi5 - psi5,  tiny2)
             ghcl  = max(chi6 - omehi, tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t - pk_t - 2.0d0*mg_t
             scon = k4*lwnsq   
    !
             if (smin > tiny) then                        ! H+ in excess
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                end if 
    !
                h   = max(hh,sqrt(scon))
    !           ohi = scon/h
             else                                        ! OH- in excess
                bb = smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                end if 
    !
                ohi = max(hh,sqrt(scon))
                h   = scon/ohi
             end if
    !       
    !  ## 4. Aerosol liquid water content
             frno3  = max(no3_t - nano3, 0.0d0)
             frcl   = max(cl_t  - nacl, 0.0d0)
             m5     = min(nh4_t, frno3)
             frnh4  = max(nh4_t - m5, 0.0d0)
             m6     = min(frcl, frnh4)   
             lwn    = max(ztot + m5/c5 + m6/c6, tiny)
          end if
    !
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep-1 .and.  errin >= epsact)         
             k = k + 1
    !
    !  ## Reset gamou
             if (( .not. soln) .and. frst) then
                gamou(1)   = gama(1)
                gamou(2)   = gama(2)
                gamou(3)   = gama(3)
                gamou(4)   = gama(4)
                gamou(5)   = gama(5)
                gamou(6)   = gama(6)
                gamou(7)   = gama(7)
                gamou(8)   = gama(8)
                gamou(9)   = gama(9)
                gamou(10)  = gama(10)
                gamou(11)  = gama(11)
                gamou(12)  = gama(12)
                gamou(13)  = gama(13)
                gamou(14)  = gama(14)
                gamou(15)  = gama(15)
                gamou(16)  = gama(16)
                gamou(17)  = gama(17)
                gamou(18)  = gama(18)
                gamou(19)  = gama(19)
                gamou(20)  = gama(20)
                gamou(21)  = gama(21)
                gamou(22)  = gama(22)
                gamou(23)  = gama(23)
             end if
    !
    !  ## Reset gamin
             if (( .not. soln)) then
                gamin(1)   = gama(1)
                gamin(2)   = gama(2)
                gamin(3)   = gama(3)
                gamin(4)   = gama(4)
                gamin(5)   = gama(5)
                gamin(6)   = gama(6)
                gamin(7)   = gama(7)
                gamin(8)   = gama(8)
                gamin(9)   = gama(9)
                gamin(10)  = gama(10)
                gamin(11)  = gama(11)
                gamin(12)  = gama(12)
                gamin(13)  = gama(13)
                gamin(14)  = gama(14)
                gamin(15)  = gama(15)
                gamin(16)  = gama(16)
                gamin(17)  = gama(17)
                gamin(18)  = gama(18)
                gamin(19)  = gama(19)
                gamin(20)  = gama(20)
                gamin(21)  = gama(21)
                gamin(22)  = gama(22)
                gamin(23)  = gama(23)
             end if 
    !
             call mach_hetp_calcact4b(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                      ca_t, pk_t, mg_t, lwn, gama, t, soln, frst,   &
                                      calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 23
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    !
             errinlocb = 0
             do ii = 1, 23
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients 
             do ii = 1, 23
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do       
    !
    !  ## Solve system of equations, with new activity coefficients
             if (( .not. soln)) then
                lwnsq = lwn*lwn
                gama10sq = gama(10)*gama(10)
    !
                a4 = k1*gama10sq/(gama(5)*gama(5))    
                a5 = k2*lwnsq/gama10sq                
                a6 = k3*lwnsq/(gama(11)*gama(11))     
    !
    !  ## 1. Calculate dissociation quantities 
                psi5 = chi5*(omehi + nacl) - a6/a5*nano3*(chi6 - omehi)
                psi5 = psi5/(a6/a5*(chi6 - omehi) + omehi + nacl)
                psi5 = min(max(psi5, tiny), chi5)
    !
                if (nh4 > tiny .and. lwn > tiny) then
                   bb   = -(nh4 + omehi + psi5 + 1.0d0/a4)
                   cc   = nh4*(psi5 + omehi)
                   psi4 = min(max(0.5d0*(-bb - sqrt(max(bb*bb - 4.0d0*cc, 0.0d0))), 0.0d0), nh4)
                else
                   psi4 = tiny
                end if
    !
    !  ## 2. Speciation 
                nh4_t = psi4
                cl_t  = omehi + nacl
                no3_t = psi5 + nano3
                gnh3  = max(nh4  - psi4,  tiny2) 
                ghno3 = max(chi5 - psi5,  tiny2)
                ghcl  = max(chi6 - omehi, tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t - pk_t - 2.0d0*mg_t
                scon = k4*lwnsq   
    !
                if (smin > tiny) then                        ! H+ in excess        
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else 
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                   end if 
    !
                   h   = max(hh,sqrt(scon))
    !              ohi = scon/h
                else                                        ! OH- in excess
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                   end if 
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if
    !          
    !  ## 4. Aerosol liquid water content
                frno3  = max(no3_t - nano3, 0.0d0)
                frcl   = max(cl_t  - nacl, 0.0d0)
                m5     = min(nh4_t, frno3)
                frnh4  = max(nh4_t - m5, 0.0d0)
                m6     = min(frcl, frnh4)   
                lwn    = max(ztot + m5/c5 + m6/c6, tiny)
             end if 
          end do
    !     
          y2 = h*cl_t/ghcl/a6 - 1.0d0  !Function value 
    !
    !  ## Check for criteria to exit root tracking 
          condition = .false.
    !
         loccon = sign(1.0d0,y1)*sign(1.0d0,y2)
         if (loccon > 0.0d0 .and. ( .not. noroot) .and. abs(y2) > eps .and. chi6 > tiny) then
             condition = .true.         
          elseif (loccon < 0.0d0 .and. abs(y2) > eps) then
    !  ## Interval had been found where sign change occurs; exit root tracking and proceed to ITP
             soln = .true. 
          else
    !  ## abs(y2) <= eps; solution is assumed; exit root tracking and proceed to minor system (no ITP)
             soln   = .true.
             noroot = .true.
             earlye = .true.
          end if 
    !
    !  ## Too little frcl, or 'tiny' is a root; reset x-value to tiny and exit root tracking 
          if (chi6 <= tiny .or. earlyexit) then
             condition = .false.
             rooteval  = 2       ! Increment rooteval by 1 to force exit from root tracking
             noroot    = .true.
             omehi     = tiny    ! Reset x-value to tiny, and solve system with this value
             omebe     = omehi
          end if 
    !
    !  ### AFTER iterating through ALL ndiv subdivided intervals
          if (rooteval == ndiv + 1) then
             loccon = sign(1.0d0,y1)*sign(1.0d0,y2)
             if (loccon > 0.0d0 .and. abs(y2) > eps) then
    !  ## (1) No solution
                noroot = .true.
                omehi  = tiny
                omebe  = omehi
    !           write(*,*), 'Warning in CALCM8: no solution found'
             else if (loccon > 0.0d0 .and. abs(y2) <= eps) then
    !  ## (2) Solution is assumed and ITP is not required
                noroot = .true.
             end if 
          end if 
       end do  !End outer loop for root tracking
    !
    !
    !
    !  ### STAGE 2: modified bisection search (using ITP algorithm) ###
    !  ## Initialize static ITP variables 
       if (( .not. noroot)) then
          ya = y1
          yb = y2
          xa = omebe
          xb = omehi
          x3 = omehi
    !
          if (xa == xb) then
             noroot = .true.
             gx     = tiny
          else
             gx     = xb - xa
          end if 
    !
          gx2  = (xa+xb)*0.5d0
          u1   = 0.2d0/gx
          nh   = log10(abs(gx/(2.0d0*eps*gx2))) / log10of2  
          nmax = int(nh) + 2                               
       else
          x3 = omehi
       end if 
    !
    !  ## Start search
    !  1. Store value of previous y3 for later use
       y3_lastiter = 0.0d0
       x3_lastiter = 0.0d0
       y3_min = 1.0d20
       x3_min = 0.0d0
    !
       j = 0
       condition = .true.
    !
       if (earlye) then
          soln = .true.
       else
          soln = .false.
       end if
    !
       do while (j < maxit .and. condition)   ! Begin outer loop for ITP search
    !  1. Track x3 and y3 of the previous iteration
          if (j > 0) then
             y3_lastiter = y3
             x3_lastiter = x3
          end if
    !
    !  2. Track the minimum y3 that is found before ending
          if (abs(0.0d0 - y3_min) > abs(0.0d0 - y3_lastiter) .and. j > 0) then
             y3_min    = y3_lastiter
             x3_min    = x3_lastiter
             ghno3_min = ghno3
             no3_min   = no3_t
             h_min     = h
             ghcl_min  = ghcl
             cl_min    = cl_t
             gnh3_min  = gnh3
             nh4_min   = nh4_t
             lwn_min   = lwn
             gama_min  = gama
          end if
    !
          if (( .not. noroot) .and. ( .not. soln)) then
    !  ## Set dynamic ITP variables 
             if (yb - ya == 0.0d0) then
                write(*,*), '######        ABORT       ######'
                write(*,*), 'Zero divide in ITP reset: CALCM8'
                write(*,*), 'SO4 in = ', so4
                write(*,*), 'NH4 in = ', nh4
                write(*,*), 'NO3 in = ', no3
                write(*,*), 'Na in  = ', na
                write(*,*), 'Cl in  = ', cl
                write(*,*), 'Ca in  = ', ca
                write(*,*), 'Mg in  = ', mg
                write(*,*), 'K in   = ', pk
                write(*,*), 'Temp in= ', t
                write(*,*), 'RH in  = ', aw
                return
             end if  
    !
             gx   = xb - xa
             xh   = 0.5d0*(xa + xb)
             rr   = max(gx2*eps*2.d0**(real(nmax - j)) - 0.5d0*gx, 0.0d0)
             delta= u1*(max(gx, 0.0d0))**2.0d0   
             xf   = max((yb*xa - ya*xb) / (yb - ya), 0.0d0)
    !
             sigma = sign(1.0d0, xh - xf)
             if (delta <= abs(xh - xf)) then        
                xt = xf + sigma*delta          
             else
                xt = xh
             end if
    !
             if (abs(xt - xh) <= rr) then
                x3 = xt
             else
                x3 = xh - sigma*rr
             end if 
          end if 
    !
          j = j + 1
    !
          if (( .not. soln)) then
             gmax = 0.1d0
             gmax = max(gmax, gama(1))
             gmax = max(gmax, gama(2))
             gmax = max(gmax, gama(3))
             gmax = max(gmax, gama(4))
             gmax = max(gmax, gama(5))
             gmax = max(gmax, gama(6))
             gmax = max(gmax, gama(7))
             gmax = max(gmax, gama(8))
             gmax = max(gmax, gama(9))
             gmax = max(gmax, gama(10))
             gmax = max(gmax, gama(11))
             gmax = max(gmax, gama(12))
             gmax = max(gmax, gama(13))
             gmax = max(gmax, gama(14))
             gmax = max(gmax, gama(15))
             gmax = max(gmax, gama(16))
             gmax = max(gmax, gama(17))
             gmax = max(gmax, gama(18))
             gmax = max(gmax, gama(19))
             gmax = max(gmax, gama(20))
             gmax = max(gmax, gama(21))
             gmax = max(gmax, gama(22))
             gmax = max(gmax, gama(23))
          end if
    !
    !  ## Reinitialize activity coefficients if gmax > 100.0d0
          if (gmax > 100.0d0 .and. ( .not. soln)) then
             gama = 0.1d0
             gamin = 1.0d10
             gamou = 1.0d10
             calou = .true.
             frst  = .true.
          end if
    !
    !  ## Solve system of equations
          frst    = .true.
          calain  = .true.
          if (( .not. soln)) then
             lwnsq = lwn*lwn
             gama10sq = gama(10)*gama(10)
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))    
             a5 = k2*lwnsq/gama10sq                
             a6 = k3*lwnsq/(gama(11)*gama(11))     
    !
    !  ## 1. Calculate dissociation quantities 
             psi5 = chi5*(x3 + nacl) - a6/a5*nano3*(chi6 - x3)
             psi5 = psi5/(a6/a5*(chi6 - x3) + x3 + nacl)
             psi5 = min(max(psi5, tiny), chi5)
    !
             if (nh4 > tiny .and. lwn > tiny) then
                bb   = -(nh4 + x3 + psi5 + 1.0d0/a4)
                cc   = nh4*(psi5 + x3)
                psi4 = min(max(0.5d0*(-bb - sqrt(max(bb*bb - 4.0d0*cc, 0.0d0))), 0.0d0), nh4)
             else
                psi4 = tiny
             end if
    !
    !  ## 2. Speciation 
             nh4_t = psi4
             cl_t  = x3 + nacl
             no3_t = psi5 + nano3   
             gnh3  = max(nh4  - psi4, tiny2) 
             ghno3 = max(chi5 - psi5, tiny2)
             ghcl  = max(chi6 - x3,   tiny2)           
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t - pk_t - 2.0d0*mg_t
             scon = k4*lwnsq 
    !
             if (smin > tiny) then                        ! H+ in excess    
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                end if 
    !
                h   = max(hh,sqrt(scon))
    !           ohi = scon/h 
             else                                        ! OH- in excess   
                bb = smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                end if 
    !
                ohi = max(hh,sqrt(scon))
                h   = scon/ohi
             end if 
    !
    !  ## 4. Aerosol liquid water content
             frno3  = max(no3_t - nano3, 0.0d0)
             frcl   = max(cl_t  - nacl, 0.0d0)
             m5     = min(nh4_t, frno3)
             frnh4  = max(nh4_t - m5, 0.0d0)
             m6     = min(frcl, frnh4)   
             lwn    = max(ztot + m5/c5 + m6/c6, tiny)
          end if 
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 0
          do while ( k < nsweep .and.  errin >= epsact)         
             k = k + 1
    !
    !  ## Reset gamin and gamou
             if (( .not. soln) .and. frst) then
                gamou(1)   = gama(1)
                gamou(2)   = gama(2)
                gamou(3)   = gama(3)
                gamou(4)   = gama(4)
                gamou(5)   = gama(5)
                gamou(6)   = gama(6)
                gamou(7)   = gama(7)
                gamou(8)   = gama(8)
                gamou(9)   = gama(9)
                gamou(10)  = gama(10)
                gamou(11)  = gama(11)
                gamou(12)  = gama(12)
                gamou(13)  = gama(13)
                gamou(14)  = gama(14)
                gamou(15)  = gama(15)
                gamou(16)  = gama(16)
                gamou(17)  = gama(17)
                gamou(18)  = gama(18)
                gamou(19)  = gama(19)
                gamou(20)  = gama(20)
                gamou(21)  = gama(21)
                gamou(22)  = gama(22)
                gamou(23)  = gama(23)
             end if 
    !
             if (( .not. soln)) then
                gamin(1)   = gama(1)
                gamin(2)   = gama(2)
                gamin(3)   = gama(3)
                gamin(4)   = gama(4)
                gamin(5)   = gama(5)
                gamin(6)   = gama(6)
                gamin(7)   = gama(7)
                gamin(8)   = gama(8)
                gamin(9)   = gama(9)
                gamin(10)  = gama(10)
                gamin(11)  = gama(11)
                gamin(12)  = gama(12)
                gamin(13)  = gama(13)
                gamin(14)  = gama(14)
                gamin(15)  = gama(15)
                gamin(16)  = gama(16)
                gamin(17)  = gama(17)
                gamin(18)  = gama(18)
                gamin(19)  = gama(19)
                gamin(20)  = gama(20)
                gamin(21)  = gama(21)
                gamin(22)  = gama(22)
                gamin(23)  = gama(23)
             end if 
    !
             call mach_hetp_calcact4b(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                      ca_t, pk_t, mg_t, lwn, gama, t, soln,  frst,  &
                                      calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 23
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst  = .false.
             end if 
    !
             errinlocb = 0
             do ii = 1, 23
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients 
             do ii = 1, 23
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do    
    !
    !  ## Solve system of equations, with new activity coefficients
             if (( .not. soln)) then
                lwnsq    = lwn*lwn
                gama10sq = gama(10)*gama(10)
    !
                a4 = k1*gama10sq/(gama(5)*gama(5))    
                a5 = k2*lwnsq/gama10sq                
                a6 = k3*lwnsq/(gama(11)*gama(11))     
    !
    !  ## 1. Calculate dissociation quantities 
                psi5 = chi5*(x3 + nacl) - a6/a5*nano3*(chi6 - x3)
                psi5 = psi5/(a6/a5*(chi6 - x3) + x3 + nacl)
                psi5 = min(max(psi5, tiny), chi5)
    !
                if (nh4 > tiny .and. lwn > tiny) then
                   bb   = -(nh4 + x3 + psi5 + 1.0d0/a4)
                   cc   = nh4*(psi5 + x3)
                   psi4 = min(max(0.5d0*(-bb - sqrt(max(bb*bb - 4.0d0*cc, 0.0d0))), 0.0d0), nh4)
                else
                   psi4 = tiny
                end if
    !
    !  ## 2. Speciation 
                nh4_t = psi4
                cl_t  = x3 + nacl
                no3_t = psi5 + nano3 
                gnh3  = max(nh4  - psi4, tiny2) 
                ghno3 = max(chi5 - psi5, tiny2)
                ghcl  = max(chi6 - x3,   tiny2)                                          
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t  &
                       - pk_t - 2.0d0*mg_t
                scon = k4*lwnsq    
    !
                if (smin > tiny) then                        ! H+ in excess    
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                   end if 
    !
                   h   = max(hh,sqrt(scon))
    !              ohi = scon/h
                else                                        ! OH- in excess 
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
                   end if 
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if
    !
    !  ## 4. Aerosol liquid water content
                frno3  = max(no3_t - nano3, 0.0d0)
                frcl   = max(cl_t  - nacl, 0.0d0)
                m5     = min(nh4_t, frno3)
                frnh4  = max(nh4_t - m5, 0.0d0)
                m6     = min(frcl, frnh4)   
                lwn    = max(ztot + m5/c5 + m6/c6, tiny)
             end if 
          end do
    !     
          y3 = h*cl_t/ghcl/a6 - 1.0d0  !Function value    
    !
          condition = .false.
          if (noroot) then
    !  ## If no root on interval then do not perform ITP
             xa = x3
             xb = x3
          else if (y3 > 0.0d0 .and. ( .not. soln)) then
             xb = x3
             yb = y3
          else if (y3 < 0.0d0 .and. ( .not. soln)) then
             xa = x3
             ya = y3
          else if (( .not. soln)) then
             xa = x3
             xb = x3
          end if
    !
    !  ## Check for convergence criteria to exit ITP
          if (xb - xa > abs(xa*eps) .and. chi6 > tiny .and. ( .not. noroot)) then
             condition = .true.
             soln   = .false.
          else
             soln   = .true.
          end if
    !
    ! ## Exit ITP if the function being bisected evaluates to a value <= eps; solution is assumed
          if (abs(y3) <= eps .and. ( .not. noroot)) then
             soln = .true.
             condition = .false.
          end if
    !
    ! ## Post-convergence correction: 
    ! ## If the calculated y-value (y3) after ITP is further from zero than an earlier iteration 
    ! ## then reset to the x-value/concentrations/activity coefficients that were found to minimize
    ! ## the objective function (i.e., y3); in this case, this is chosen as the solution
          if (( .not. condition) .and. ( .not. noroot) .and. abs(y3) > 0.1d0) then     
             if (abs(0.0d0 - y3_min) < abs(0.0d0 - y3) .and. abs(y3_min - y3) > 1.0d-1) then
                x3    = x3_min
                cl_t  = cl_min
                nh4_t = nh4_min
                no3_t = no3_min
                h     = h_min
                lwn   = lwn_min
                ghcl  = ghcl_min
                ghno3 = ghno3_min
                gnh3  = gnh3_min
    ! 
    !  ## Reset activity coefficients 
                gama  = gama_min
                a6    = k3*(lwn/gama(11))*(lwn/gama(11))
                y3    = h*cl_t/ghcl/a6 - 1.0d0
             end if
          end if
       end do ! End outer loop of ITP search
    !
    !
    !  ### MINOR SYSTEM: HSO4-/SO42-/H+ ###
       hso4 = 0.0d0
       if (h > tiny .and. so4_t > tiny .and. lwn >= 1.0d-19) then
          ak1 = khso4*lwn/gama(7)*(gama(8)/gama(7))*(gama(8)/gama(7))
          bb  =-(h + so4_t + ak1)
          cc  = h*so4_t 
          dd  = bb*bb - 4.d0*cc
    ! 
          if (dd >= 0.0d0) then
             if (bb /= 0.d0) then
                dd = cc/(bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d3
             end if
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                hh = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
             else
                hh = 0.5d0*(-bb - sqrt(bb*bb - 4.d0*cc))
             end if 
    !
             hh    = max(tiny, min(hh, min(h, so4_t)))     ! To avoid negative H+   (i.e., if hh > h or hh > so4_t)
             h     = max(h - hh, 0.0d0)
             so4_t = max(so4_t - hh, 0.0d0)
             hso4  = hh
          end if 
       end if
    !
    !
    !  ### Perform mass adjustment if excess exists, to balance mass to machine precision ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       na_i    = na_t
       cl_i    = cl_t
       ca_i    = ca_t
       k_i     = pk_t
       mg_i    = mg_t
       nh3g_i  = gnh3
       hno3g_i = ghno3
       hclg_i  = ghcl
       h_i     = h
       lwn_i   = lwn
       caso4_i = caso4
    !
       return
    end subroutine mach_hetp_calcm8
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: P13; Sulfate poor; dust and sodium rich
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcp13(so4_i, nh4_i, nh3g_i, hno3g_i, hclg_i, hso4_i,         &
                                 na_i, cl_i, no3_i, h_i, lwn_i, ca_i, k_i, mg_i,        &
                                 caso4_i, frso4, frmg, frk, frca, frna, rh, temp,       &
                                 k0, p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: na_i    
       real(kind=8),    intent(inout) :: cl_i    
       real(kind=8),    intent(inout) :: ca_i    
       real(kind=8),    intent(inout) :: k_i     
       real(kind=8),    intent(inout) :: mg_i    
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: hclg_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent(inout) :: caso4_i 
       real(kind=8),    intent(inout) :: frso4
       real(kind=8),    intent(inout) :: frmg  
       real(kind=8),    intent(inout) :: frk    
       real(kind=8),    intent(inout) :: frca   
       real(kind=8),    intent(inout) :: frna 
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    !
    ! ## Local variables
       real(kind=8)     :: so4, nh4, hso4, gnh3, h, lwn, caso4, no3, cl, na, ca, pk, mg, ghno3, ghcl
       real(kind=8)     :: t, aw, knh3, khso4, kh2o, khno3, khcl, a4, a5, psi4, psi5, ak1, m5, m6, tt0
       real(kind=8)     :: bb, cc, dd, hh, v, errin, smin, ohi, scon, errouloc, errinlocb
       real(kind=8)     :: frcl, frno3, frnh4, cano32, kno3, kcl, loccon
       real(kind=8)     :: omehi, omebe, y1, y2, y3, x3, dx, a6, tt1, tt2, k1, k2, ya, yb, xa, xb
       real(kind=8)     :: so4_t, nh4_t, no3_t, na_t, cl_t, ca_t, pk_t, mg_t, k3, k4
       real(kind=8)     :: chi5, chi6, nacl, nano3, k2so4, mgso4, c1, c2, c3, c4, c5
       real(kind=8)     :: mgno32, mgcl2, cacl2, gmax, ztot, lwnsq, gama10sq
       real(kind=8)     :: nh, sigma, xt, xf, xh, delta, rr, gx, gx2, u1
       integer(kind=4)  :: j, k, ii, rooteval, irh, nmax
       logical(kind=4)  :: condition, noroot, earlyexit, soln, frst, calain, calou, earlye
       real(kind=8), dimension(23) :: gama, gamin, gamou, gama_min
       real(kind=8)      :: y3_min, x3_min, y3_lastiter, x3_lastiter
       real(kind=8)      :: no3_min, nh4_min, cl_min, h_min, lwn_min, ghcl_min, gnh3_min, ghno3_min
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       na    = na_i
       cl    = cl_i
       ca    = ca_i
       pk    = k_i
       mg    = mg_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       gnh3  = 0.0d0
       ghno3 = 0.0d0
       ghcl  = 0.0d0
       h     = 0.0d0
       lwn   = tiny
       so4_t = 0.0d0
       nh4_t = 0.0d0
       no3_t = 0.0d0
       na_t  = 0.0d0
       cl_t  = 0.0d0
       ca_t  = 0.0d0
       pk_t  = 0.0d0
       mg_t  = 0.0d0
       caso4 = 0.0d0
       noroot=.false.
       earlyexit = .false.
       soln  = .false.
       calou = .true.
       gmax  = 0.0d0
       gama  = 0.1d0
       gamin = 1.0d10
       gamou = 0.1d0
       earlye = .false.
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    !
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3= k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. HCl(g) <==> H+(aq) + Cl-(aq)                               (xk3)
       khcl = k0(11) * exp(p1(11)*tt1 + p2(11)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Calculate dry salt composition 
    !  ## Salts assumed to have completely dissolved
    !  ## Dissolved salts include: Ca(NO3)2, CaCl2, K2SO4, KNO3, KCl, MgSO4, Mg(NO3)2, MgCl2, NaNO3, NaCl
       caso4 = min(so4, ca)                        
       frca  = max(ca - caso4, 0.0d0)              
       frso4 = max(so4 - caso4, 0.0d0)             
       k2so4 = min(frso4, 0.5d0*pk)                   
       frk   = max(pk - 2.0d0*k2so4, 0.0d0)        
       frso4 = max(frso4 - k2so4, 0.0d0)             
       mgso4 = min(frso4, mg)                                     
       frmg  = max(mg - mgso4, 0.0d0)      
       frso4 = max(frso4 - mgso4, 0.0d0)     
       nacl  = min(na, cl)                         
       frna  = max(na - nacl, 0.0d0)               
       frcl  = max(cl - nacl, 0.0d0)               
       cano32= min(frca, 0.5d0*no3)                   
       frca  = max(frca - cano32, 0.0d0)              
       frno3 = max(no3 - 2.0d0*cano32, 0.0d0)      
       cacl2 = min(frca, 0.5d0*frcl)                     
       frca  = max(frca - cacl2, 0.0d0)               
       frcl  = max(frcl - 2.0d0*cacl2, 0.0d0)         
       mgno32= min(frmg, 0.5d0*frno3)                    
       frmg  = max(frmg - mgno32, 0.0d0)             
       frno3 = max(frno3 - 2.0d0*mgno32, 0.0d0)       
       mgcl2 = min(frmg, 0.5d0*frcl)                     
       frmg  = max(frmg - mgcl2, 0.0d0)               
       frcl  = max(frcl - 2.0d0*mgcl2, 0.0d0)        
       nano3 = min(frna, frno3)                          
       frna  = max(frna - nano3, 0.0d0)               
       frno3 = max(frno3 - nano3, 0.0d0)             
       kcl   = min(frk, frcl)                            
       frk   = max(frk - kcl, 0.0d0)               
       frcl  = max(frcl - kcl, 0.0d0)                 
       kno3  = min(frk, frno3)                          
       frk   = max(frk - kno3, 0.0d0)               
       frno3 = max(frno3 - kno3, 0.0d0)             
       chi5  = frno3                                     ! HNO3(g)
       chi6  = frcl                                      ! HCl(g)
       frno3 = 0.0d0
       frcl  = 0.0d0
    !
    !  ## Initial speciation
       na_t  = nacl + nano3
       so4_t = k2so4 + mgso4
       ca_t  = cano32 + cacl2
       pk_t  = 2.0d0*k2so4 + kno3 + kcl
       mg_t  = mgso4 + mgno32 + mgcl2
    !
    !  ## Calculate constant parameters
       c1 = r*t
       c2 = nacl + kcl + 2.0d0*mgcl2 + 2.0d0*cacl2
       c3 = nano3 + 2.0d0*cano32 + kno3 + 2.0d0*mgno32
       c4 = -1.0d0*nano3 - 2.0d0*cano32 - kno3 - 2.0d0*mgno32
       c5 = -1.0d0*nacl - 2.0d0*cacl2 - kcl - 2.0d0*mgcl2
       k1 = (knh3/kh2o)*c1
       k2 = khno3*c1
       k3 = khcl*c1
       k4 = kh2o*aw
    !
    !  ## Constant ZSR parameters
       ztot = nacl/awsc(irh)   + nano3/awsn(irh) + cano32/awcn(irh) + &
              cacl2/awcc(irh)  + kno3/awpn(irh)  + kcl/awpc(irh)    + &
              mgno32/awmn(irh) + mgcl2/awmc(irh) + k2so4/awps(irh)  + &
              mgso4/awms(irh)
    !
    !
    !  ### STAGE 1: Root tracking ###
    !  ## Find a subinterval [xa,xb] on the larger interval [I1,I2] where a sign change occurs
       rooteval = 0
       condition = .true.
       do while (rooteval < 2 .or. (condition .and. rooteval < ndiv + 1))   ! Begin outer loop for root tracking
          rooteval = rooteval + 1
    !
    !  ## Set high limit for root tracking (i.e. lower bound)
          if (rooteval == 1) then
            omehi = tiny
            y1    = 1.0d0
          end if
    !
    !  ## Begin search on subinterval 
          if (rooteval == 2) then
             soln = .false.
             if (abs(y2) <= eps) then
                earlyexit = .true.
             end if
    !
             y1    = y2
    !
             if (earlyexit) then
                dx = 0.0d0
             else
                dx = (chi6-tiny-tiny)/float(ndiv)
             end if
    !
             omebe = omehi            ! Lower bound of subinterval (xa)
             omehi = omehi + dx       ! Upper bound of subinterval (xb)
          end if
    !
    !  ## Continue search 
          if (rooteval > 2) then
             if (sign(1.0d0,y1)*sign(1.0d0,y2) < 0.0d0) then
    !  ## 1. Root has been found on the subinterval; save x values for ITP search
                y1    = y1
                omebe = omebe
                omehi = omehi
             else
    !  ## 2. No root has been found, continue searching in the next subinterval
                y1    = y2
                omebe = omehi
                omehi = omehi + dx
             end if
          end if
    !
    !  ## Solve the system of equations 
          frst   = .true.
          calain = .true.
          if (( .not. soln)) then
             lwnsq    = lwn*lwn
             gama10sq = gama(10)*gama(10) 
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))    
             a5 = k2*lwnsq/gama10sq                
             a6 = k3*lwnsq/(gama(11)*gama(11))     
    !
    !  ## 1. Calculate dissociation quantities
             psi5 = chi5*(omehi + c2) - a6/a5*(c3)*(chi6 - omehi)
             psi5 = psi5/(a6/a5*(chi6 - omehi) + omehi + c2)
             psi5 = min(max(psi5, tiny), chi5)
    
             if (nh4 > tiny .and. lwn > tiny) then
                bb = -(nh4 + omehi + psi5 + 1.0d0/a4)
                cc = nh4*(psi5 + omehi)
    
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   psi4 = -0.5*bb - 0.5*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                else
                   psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    
                psi4 = min(max(psi4, 0.0d0), nh4)
             else
                psi4 = tiny
             end if
    !
    !  ## 2. Speciation
             nh4_t = psi4
             cl_t  = omehi + c2
             no3_t = psi5 + c3
             gnh3  = max(nh4  - psi4,  tiny2)
             ghno3 = max(chi5 - psi5,  tiny2)
             ghcl  = max(chi6 - omehi, tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t  &
                    - pk_t - 2.0d0*mg_t - 2.0d0*ca_t
             scon = k4*lwnsq   
    !
             if (smin > tiny) then                        ! H+ in excess
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                h   = max(hh,sqrt(scon))
    !           ohi = scon/h
             else                                        ! OH- in excess
                bb = smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                ohi = max(hh,sqrt(scon))
                h   = scon/ohi
             end if
    !
    !  ## 4. Aerosol liquid water content
             frno3  = max(no3_t + c4, 0.0d0)
             frcl   = max(cl_t  + c5, 0.0d0)
             m5     = min(nh4_t, frno3)
             frnh4  = max(nh4_t - m5, 0.0d0)
             m6     = min(frcl, frnh4)
             lwn    = max(ztot + m5/awan(irh) + m6/awac(irh), tiny)
          end if 
    !
    !  ## Iterate until convergence of activity coefficients
          errin = 1.0d0
          k     = 1
          do while ( k < nsweep .and.  errin >= epsact)
             k = k + 1
    !
    !  ## Reset gamou
             if (( .not. soln) .and. frst) then
                gamou(1)   = gama(1)
                gamou(2)   = gama(2)
                gamou(3)   = gama(3)
                gamou(4)   = gama(4)
                gamou(5)   = gama(5)
                gamou(6)   = gama(6)
                gamou(7)   = gama(7)
                gamou(8)   = gama(8)
                gamou(9)   = gama(9)
                gamou(10)  = gama(10)
                gamou(11)  = gama(11)
                gamou(12)  = gama(12)
                gamou(13)  = gama(13)
                gamou(14)  = gama(14)
                gamou(15)  = gama(15)
                gamou(16)  = gama(16)
                gamou(17)  = gama(17)
                gamou(18)  = gama(18)
                gamou(19)  = gama(19)
                gamou(20)  = gama(20)
                gamou(21)  = gama(21)
                gamou(22)  = gama(22)
                gamou(23)  = gama(23)
             end if 
    !
    !  ## Reset gamin
             if (( .not. soln)) then
                gamin(1)   = gama(1)
                gamin(2)   = gama(2)
                gamin(3)   = gama(3)
                gamin(4)   = gama(4)
                gamin(5)   = gama(5)
                gamin(6)   = gama(6)
                gamin(7)   = gama(7)
                gamin(8)   = gama(8)
                gamin(9)   = gama(9)
                gamin(10)  = gama(10)
                gamin(11)  = gama(11)
                gamin(12)  = gama(12)
                gamin(13)  = gama(13)
                gamin(14)  = gama(14)
                gamin(15)  = gama(15)
                gamin(16)  = gama(16)
                gamin(17)  = gama(17)
                gamin(18)  = gama(18)
                gamin(19)  = gama(19)
                gamin(20)  = gama(20)
                gamin(21)  = gama(21)
                gamin(22)  = gama(22)
                gamin(23)  = gama(23)
             end if 
    !
             call mach_hetp_calcact4b( h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,    &
                                      ca_t, pk_t, mg_t, lwn, gama, t, soln, frst,   &
                                      calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 23
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    
             errinlocb = 0
             do ii = 1, 23
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients
             do ii = 1, 23
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do
    !
    !  ## Solve system of equations, using new activity coefficients 
             if (( .not. soln)) then
                lwnsq    = lwn*lwn 
                gama10sq = gama(10)*gama(10)
    !
                a4 = k1*gama10sq/(gama(5)*gama(5))   
                a5 = k2*lwnsq/gama10sq               
                a6 = k3*lwnsq/(gama(11)*gama(11))    
    !
    !  ## 1. Calculate dissociation quantities
                psi5 = chi5*(omehi + c2) - a6/a5*(c3)*(chi6 - omehi)
                psi5 = psi5/(a6/a5*(chi6 - omehi) + omehi + c2)
                psi5 = min(max(psi5, tiny), chi5)
    
                if (nh4 > tiny .and. lwn > tiny) then
                   bb = -(nh4 + omehi + psi5 + 1.0d0/a4)
                   cc = nh4*(psi5 + omehi)
    
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      psi4 = -0.5*bb - 0.5*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                   else
                      psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    
                   psi4 = min(max(psi4, 0.0d0), nh4)
                else
                   psi4 = tiny
                end if
    !
    !  ## 2. Speciation
                nh4_t = psi4
                cl_t  = omehi + c2
                no3_t = psi5 + c3
                gnh3  = max(nh4  - psi4,  tiny2)
                ghno3 = max(chi5 - psi5,  tiny2)
                ghcl  = max(chi6 - omehi, tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t  &
                       - pk_t - 2.0d0*mg_t - 2.0d0*ca_t
                scon = k4*lwnsq   
    !
                if (smin > tiny) then                        ! H+ in excess
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   h   = max(hh,sqrt(scon))
    !              ohi = scon/h
                else                                        ! OH- in excess
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if
    !
    !  ## 4. Aerosol liquid water content
                frno3  = max(no3_t + c4, 0.0d0)
                frcl   = max(cl_t  + c5, 0.0d0)
                m5     = min(nh4_t, frno3)
                frnh4  = max(nh4_t - m5, 0.0d0)
                m6     = min(frcl, frnh4)
                lwn    = max(ztot + m5/awan(irh) + m6/awac(irh), tiny)
             end if 
          end do
    !
          y2 = h*cl_t/ghcl/a6 - 1.0d0  !Function value
    !
    !  ## Check for criteria to exit root tracking 
          condition = .false.
          loccon = sign(1.0d0,y1)*sign(1.0d0,y2)
          if (loccon > 0.0d0 .and. ( .not. noroot) .and. abs(y2) > eps .and. chi6 > tiny) then
             condition = .true.         
          elseif (loccon < 0.0d0 .and. abs(y2) > eps) then
    !  ## Interval had been found where sign change occurs; exit root tracking and proceed to ITP
             soln = .true. 
          else
    !  ## abs(y2) <= eps; solution is assumed; exit root tracking and proceed to minor system (no ITP)
             soln   = .true.
             noroot = .true.
             earlye = .true.
          end if 
    !
    !  ## Too little frcl, or 'tiny' is a root; reset x-value to tiny and exit root tracking 
          if (chi6 <= tiny .or. earlyexit) then
             condition = .false.
             rooteval  = 2       ! Increment rooteval by 1 to force exit from root tracking
             noroot    = .true.
             omehi     = tiny    ! Reset x-value to tiny, and solve system with this value
             omebe     = omehi
          end if
    !
    !  ### AFTER iterating through ALL ndiv subdivided intervals
          if (rooteval == ndiv + 1) then
             loccon = sign(1.0d0,y1)*sign(1.0d0,y2)
             if (loccon > 0.0d0 .and. abs(y2) > eps) then
    !  ## (1) No solution
                noroot = .true.
                omehi  = tiny
                omebe  = omehi
    !           write(*,*), 'Warning in CALCP13: no solution found'
             else if (loccon > 0.0d0 .and. abs(y2) <= eps) then
    !  ## (2) Solution is assumed and ITP is not required
                noroot = .true.
             end if
          end if
       end do  !End outer loop for root tracking
    !
    !
    !
    !  ### STAGE 2: modified bisection search (using ITP algorithm) ###
    !  ## Initialize static ITP variables   
       if (( .not. noroot)) then
          ya = y1
          yb = y2
          xa = omebe
          xb = omehi
          x3 = omehi
    !
          if (xa == xb) then
             noroot = .true.
             gx = tiny
          else
             gx = xb - xa
          end if 
    !
          gx2  = (xa+xb)*0.5d0
          u1   = 0.2d0/gx
          nh   = log10(abs(gx/(2.0d0*eps*gx2))) / log10of2 
          nmax = int(nh) + 2                              
       else 
          x3 = omehi
       end if                                      
    !
    !  ## Start search
       y3_lastiter = 0.0d0
       x3_lastiter = 0.0d0
       y3_min = 1.0d20
       x3_min = 0.0d0
    !
       j = 0
       condition = .true.
    !
       if (earlye) then
          soln = .true.
       else
          soln = .false.
       end if
    !
       do while (j < maxit .and. condition)   ! Begin outer loop for ITP search
    !  ## Set dynamic ITP variables 
    !  1. Track x3 and y3 of the previous iteration
          if (j > 0) then
             y3_lastiter = y3
             x3_lastiter = x3
          end if 
    !
    !  2. Track the minimum y3 that is found before ending
          if (abs(0.0d0 - y3_min) > abs(0.0d0 - y3_lastiter) .and. j > 0) then
             y3_min    = y3_lastiter
             x3_min    = x3_lastiter
             ghno3_min = ghno3
             no3_min   = no3_t
             h_min     = h
             ghcl_min  = ghcl
             cl_min    = cl_t
             gnh3_min  = gnh3
             nh4_min   = nh4_t
             lwn_min   = lwn
             gama_min  = gama
          end if
    !
          if (( .not. noroot) .and. ( .not. soln)) then
             if (yb - ya == 0.0d0) then
                write(*,*), '######        ABORT       ######'
                write(*,*), 'Zero divide in ITP reset: CALCP13'
                write(*,*), 'SO4 in = ', so4
                write(*,*), 'NH4 in = ', nh4
                write(*,*), 'NO3 in = ', no3
                write(*,*), 'Na in  = ', na
                write(*,*), 'Cl in  = ', cl
                write(*,*), 'Temp in= ', t
                write(*,*), 'RH in  = ', aw
                return
             end if  
    !
             gx   = xb - xa
             xh   = 0.5d0*(xa + xb)
             rr   = max(gx2*eps*2.d0**(real(nmax - j)) - 0.5d0*gx, 0.0d0)
             delta= u1*(max(gx, 0.0d0))**2.0d0   
             xf   = max((yb*xa - ya*xb) / (yb - ya), 0.0d0)
    !
             sigma = sign(1.0d0, xh - xf)
             if (delta <= abs(xh - xf)) then
                xt = xf + sigma*delta              
             else
                xt = xh
             end if
    !
             if (abs(xt - xh) <= rr) then
                x3 = xt
             else
                x3 = xh - sigma*rr
             end if 
          end if 
    !
          j = j + 1
    !
          if (( .not. soln)) then
             gmax = 0.1d0
             gmax = max(gmax, gama(1))
             gmax = max(gmax, gama(2))
             gmax = max(gmax, gama(3))
             gmax = max(gmax, gama(4))
             gmax = max(gmax, gama(5))
             gmax = max(gmax, gama(6))
             gmax = max(gmax, gama(7))
             gmax = max(gmax, gama(8))
             gmax = max(gmax, gama(9))
             gmax = max(gmax, gama(10))
             gmax = max(gmax, gama(11))
             gmax = max(gmax, gama(12))
             gmax = max(gmax, gama(13))
             gmax = max(gmax, gama(14))
             gmax = max(gmax, gama(15))
             gmax = max(gmax, gama(16))
             gmax = max(gmax, gama(17))
             gmax = max(gmax, gama(18))
             gmax = max(gmax, gama(19))
             gmax = max(gmax, gama(20))
             gmax = max(gmax, gama(21))
             gmax = max(gmax, gama(22))
             gmax = max(gmax, gama(23))
          end if
    !
    !  ## Reinitialize activity coefficients if gmax > 100.0d0
          if (gmax > 100.0d0 .and. ( .not. soln)) then
             gama  = 0.1d0
             gamin = 1.0d10
             gamou = 1.0d10
             calou = .true.
             frst  = .true.
          end if
    !
    !  ## Solve system of equations
          frst    = .true.
          calain  = .true.
            if (( .not. soln)) then
             lwnsq    = lwn*lwn 
             gama10sq = gama(10)*gama(10)
    !
             a4 = k1*gama10sq/(gama(5)*gama(5))   
             a5 = k2*lwnsq/gama10sq               
             a6 = k3*lwnsq/(gama(11)*gama(11))         
    !
    !  ## 1. Calculate dissociation quantities
             psi5 = chi5*(x3 + c2) - a6/a5*(c3)*(chi6 - x3)
             psi5 = psi5/(a6/a5*(chi6 - x3) + x3 + c2)
             psi5 = min(max(psi5, tiny), chi5)
    
             if (nh4 > tiny .and. lwn > tiny) then
                bb = -(nh4 + x3 + psi5 + 1.0d0/a4)
                cc = nh4*(psi5 + x3)
    
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   psi4 = -0.5*bb - 0.5*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                else
                   psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    
                psi4 = min(max(psi4, 0.0d0), nh4)
             else
                psi4 = tiny
             end if
    !
    !  ## 2. Speciation
             nh4_t = psi4
             cl_t  = x3 + c2
             no3_t = psi5  + c3
             gnh3  = max(nh4  - psi4, tiny2)
             ghno3 = max(chi5 - psi5, tiny2)
             ghcl  = max(chi6 - x3,   tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
             smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t  &
                    - pk_t - 2.0d0*mg_t - 2.0d0*ca_t
             scon = k4*lwnsq  
    !
             if (smin > tiny) then                        ! H+ in excess
                bb =-smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                h   = max(hh,sqrt(scon))
    !           ohi = scon/h
             else                                        ! OH- in excess
                bb = smin
                cc =-scon
    !
                if (bb /= 0.d0) then
                   dd = cc/(bb*bb)
                   v  = 4.d0*dd
                else
                   v  = 1.d3
                end if
    !
                if (abs(v) <= smrt .and. bb /= 0.d0) then
                   hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                else
                   hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                end if
    !
                ohi = max(hh,sqrt(scon))
                h   = scon/ohi
             end if
    !
    !  ## 4. Aerosol liquid water content
             frno3  = max(no3_t + c4, 0.0d0)
             frcl   = max(cl_t  + c5, 0.0d0)
             m5     = min(nh4_t, frno3)
             frnh4  = max(nh4_t - m5, 0.0d0)
             m6     = min(frcl, frnh4)
             lwn    = max(ztot + m5/awan(irh) + m6/awac(irh), tiny)
          end if 
    !
    !  ## Iterate until convergence of activity coefficients 
          errin = 1.0d0
          k     = 1
          do while ( k < nsweep .and.  errin >= epsact)
             k = k + 1
    !
    !  ## Reset gamin and gamou
             if (( .not. soln) .and. frst) then
                gamou(1)   = gama(1)
                gamou(2)   = gama(2)
                gamou(3)   = gama(3)
                gamou(4)   = gama(4)
                gamou(5)   = gama(5)
                gamou(6)   = gama(6)
                gamou(7)   = gama(7)
                gamou(8)   = gama(8)
                gamou(9)   = gama(9)
                gamou(10)  = gama(10)
                gamou(11)  = gama(11)
                gamou(12)  = gama(12)
                gamou(13)  = gama(13)
                gamou(14)  = gama(14)
                gamou(15)  = gama(15)
                gamou(16)  = gama(16)
                gamou(17)  = gama(17)
                gamou(18)  = gama(18)
                gamou(19)  = gama(19)
                gamou(20)  = gama(20)
                gamou(21)  = gama(21)
                gamou(22)  = gama(22)
                gamou(23)  = gama(23)
             end if 
    !
             if (( .not. soln)) then
                gamin(1)   = gama(1)
                gamin(2)   = gama(2)
                gamin(3)   = gama(3)
                gamin(4)   = gama(4)
                gamin(5)   = gama(5)
                gamin(6)   = gama(6)
                gamin(7)   = gama(7)
                gamin(8)   = gama(8)
                gamin(9)   = gama(9)
                gamin(10)  = gama(10)
                gamin(11)  = gama(11)
                gamin(12)  = gama(12)
                gamin(13)  = gama(13)
                gamin(14)  = gama(14)
                gamin(15)  = gama(15)
                gamin(16)  = gama(16)
                gamin(17)  = gama(17)
                gamin(18)  = gama(18)
                gamin(19)  = gama(19)
                gamin(20)  = gama(20)
                gamin(21)  = gama(21)
                gamin(22)  = gama(22)
                gamin(23)  = gama(23)
             end if 
    !
             call mach_hetp_calcact4b(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,    &
                                      ca_t, pk_t, mg_t, lwn, gama, t, soln, frst,  &
                                      calain, calou)
    !
             if (frst) then
                errouloc = 0
                do ii = 1, 23
                   errouloc = max(errouloc, abs(gamou(ii) - gama(ii)) / gamou(ii))
                end do
                calou = errouloc .ge. epsact
                frst   = .false.
             end if 
    !
             errinlocb = 0
             do ii = 1, 23
                errinlocb = max(errinlocb, abs(gamin(ii) - gama(ii)) / gamin(ii))
             end do
             calain = errinlocb .ge. epsact
    !
             errin = 0.0d0
    !  ## Test for convergence of activity coefficients
             do ii = 1, 23
                errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
             end do
    !
    !  ## Solve system of equations, with new activity coefficients
             if (( .not. soln)) then
                lwnsq    = lwn*lwn
                gama10sq = gama(10)*gama(10) 
    !
                a4 = k1*gama10sq/(gama(5)*gama(5))     
                a5 = k2*lwnsq/gama10sq                 
                a6 = k3*lwnsq/(gama(11)*gama(11))      
    !
    !  ## 1. Calculate dissociation quantities
                psi5 = chi5*(x3 + c2) - a6/a5*(c3)*(chi6 - x3)
                psi5 = psi5/(a6/a5*(chi6 - x3) + x3 + c2)
                psi5 = min(max(psi5, tiny), chi5)
    !
                if (nh4 > tiny .and. lwn > tiny) then
                   bb   = -(nh4 + x3 + psi5 + 1.0d0/a4)
                   cc   = nh4*(psi5 + x3)
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      psi4 = -0.5*bb - 0.5*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
                   else
                      psi4 = 0.5d0*(-bb - sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   psi4 = min(max(psi4, 0.0d0), nh4)
                else
                   psi4 = tiny
                end if
    !
    !  ## 2. Speciation
                nh4_t = psi4
                cl_t  = x3 + c2
                no3_t = psi5  + c3
                gnh3  = max(nh4  - psi4, tiny2)
                ghno3 = max(chi5 - psi5, tiny2)
                ghcl  = max(chi6 - x3,   tiny2)
    !
    !  ## 3. Calculate H+; smin = (negative charge) - (positive charge)
                smin = 2.0d0*so4_t + no3_t + cl_t - na_t - nh4_t  &
                       - pk_t - 2.0d0*mg_t - 2.0d0*ca_t
                scon = k4*lwn*lwn
    !
                if (smin > tiny) then                        ! H+ in excess
                   bb =-smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   h   = max(hh,sqrt(scon))
    !              ohi = scon/h
                else                                        ! OH- in excess
                   bb = smin
                   cc =-scon
    !
                   if (bb /= 0.d0) then
                      dd = cc/(bb*bb)
                      v  = 4.d0*dd
                   else
                      v  = 1.d3
                   end if
    !
                   if (abs(v) <= smrt .and. bb /= 0.d0) then
                      hh = -0.5d0*bb + 0.5d0*abs(bb) - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)
                   else
                      hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
                   end if
    !
                   ohi = max(hh,sqrt(scon))
                   h   = scon/ohi
                end if
    !
    !  ## 4. Aerosol liquid water content
                frno3  = max(no3_t + c4, 0.0d0)
                frcl   = max(cl_t  + c5, 0.0d0)
                m5     = min(nh4_t, frno3)
                frnh4  = max(nh4_t - m5, 0.0d0)
                m6     = min(frcl, frnh4)
                lwn    = max(ztot + m5/awan(irh) + m6/awac(irh), tiny)
             end if 
          end do
    !
          y3 = h*cl_t/ghcl/a6 - 1.0d0  !Function value
    !
          condition = .false.
          if (noroot) then
    !  ## If no root on interval then do not perform ITP
             xa = x3
             xb = x3
          else if (y3 > 0.0d0 .and. ( .not. soln)) then
             xb = x3
             yb = y3
          else if (y3 < 0.0d0 .and. ( .not. soln)) then
             xa = x3
             ya = y3
          else if (( .not. soln)) then
             xa = x3
             xb = x3
          end if
    !
    !  ## Check for convergence criteria to exit ITP:
          if (xb - xa > abs(xa*eps) .and. chi6 > tiny .and. ( .not. noroot)) then
             condition = .true.
             soln   = .false.
          else
             soln   = .true.
          end if
    !
    ! ## Exit ITP if the function being bisected evaluates to a value <= eps; solution is assumed
          if (abs(y3) <= eps .and. ( .not. noroot)) then
             soln = .true.
             condition = .false.
          end if
    !
    ! ## Post-convergence correction: 
    ! ## If the calculated y-value (y3) after ITP is further from zero than an earlier iteration 
    ! ## then reset to the x-value/concentrations/activity coefficients that were found to minimize
    ! ## the objective function (i.e., y3); in this case, this is chosen as the solution
          if (( .not. condition) .and. ( .not. noroot) .and. abs(y3) > 0.1d0) then     
             if (abs(0.0d0 - y3_min) < abs(0.0d0 - y3) .and. abs(y3_min - y3) > 1.0d-1) then
                x3    = x3_min
                cl_t  = cl_min
                nh4_t = nh4_min
                no3_t = no3_min
                h     = h_min
                lwn   = lwn_min
                ghcl  = ghcl_min
                ghno3 = ghno3_min
                gnh3  = gnh3_min
    ! 
    !  ## Reset activity coefficients 
                gama  = gama_min
                a6    = k3*(lwn/gama(11))*(lwn/gama(11))
                y3    = h*cl_t/ghcl/a6 - 1.0d0
             end if
          end if
       end do ! End outer loop of ITP search
    !
    !
    !  ### MINOR SYSTEM: HSO4-/SO42-/H+ ###
       hso4 = 0.0d0
       if (h > tiny .and. so4_t > tiny .and. lwn >= 1.0d-19) then
          ak1 = khso4*lwn/gama(7)*(gama(8)/gama(7))*(gama(8)/gama(7))
          bb  =-(h + so4_t + ak1)
          cc  = h*so4_t 
          dd  = bb*bb - 4.d0*cc
    ! 
          if (dd >= 0.0d0) then
             if (bb /= 0.d0) then
                dd = cc/(bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d3
             end if
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                hh = -0.5d0*bb - 0.5d0*abs(bb) + ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/abs(bb)   ! Negative root
             else
                hh = 0.5d0*(-bb - sqrt(bb*bb - 4.d0*cc))
             end if 
    !
             hh    = max(tiny, min(hh, min(h, so4_t)))     ! To avoid negative H+   (i.e., if hh > h or hh > so4_t)
             h     = max(h - hh, 0.0d0)
             so4_t = max(so4_t - hh, 0.0d0)
             hso4  = hh
          end if 
       end if
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)   
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       na_i    = na_t
       cl_i    = cl_t
       ca_i    = ca_t
       k_i     = pk_t
       mg_i    = mg_t
       nh3g_i  = gnh3
       hno3g_i = ghno3
       hclg_i  = ghcl
       h_i     = h
       lwn_i   = lwn
       caso4_i = caso4
    !
       return
    end subroutine mach_hetp_calcp13
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: L9; sulfate rich; no free acid 
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcl9(so4_i, nh4_i, nh3g_i, hno3g_i, hclg_i, hso4_i,          &
                                na_i, cl_i, no3_i, h_i, lwn_i, ca_i, k_i, mg_i,         &
                                caso4_i, frmg, frk, frca, frna, frso4, rh, temp, k0,    &
                                p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: hso4_i  
       real(kind=8),    intent(inout) :: na_i    
       real(kind=8),    intent(inout) :: cl_i    
       real(kind=8),    intent(inout) :: ca_i    
       real(kind=8),    intent(inout) :: k_i     
       real(kind=8),    intent(inout) :: mg_i    
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: hclg_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent(inout) :: caso4_i 
       real(kind=8),    intent(out)   :: frmg
       real(kind=8),    intent(out)   :: frk
       real(kind=8),    intent(out)   :: frca
       real(kind=8),    intent(out)   :: frna
       real(kind=8),    intent(out)   :: frso4
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    !
    !  ## Local variables
       real(kind=8)     :: so4, nh4, hso4, gnh3, h, lwn, caso4
       real(kind=8)     :: no3, cl, na, ca, pk, mg, ghno3, ghcl, t, aw
       real(kind=8)     :: knh3, khso4, kh2o, khno3, khcl, frso4_2, frnh4, frnh4_2, tt0
       real(kind=8)     :: bb, cc, dd, ff, hh, v, errin, a3, a4, a9, delno
       real(kind=8)     :: delcl, m1, m2, m3, c1, c2, tt1, tt2
       real(kind=8)     :: so4_t, nh4_t, no3_t, na_t, cl_t, ca_t, pk_t, mg_t
       real(kind=8)     :: clc, nh4hs4, nahso4, na2so4, nh42s4, k2so4, mgso4, kkhso4, exnh4
       real(kind=8)     :: frcl, frno3
       integer(kind=4)  :: j, ii, islv, irh
       logical(kind=4)  :: ispoly
       real(kind=8), dimension(23) :: gama, gamin
    !
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       na    = na_i
       cl    = cl_i
       ca    = ca_i
       pk    = k_i
       mg    = mg_i
       aw    = rh
       t     = temp
       caso4 = 0.0d0
       hso4  = 0.0d0
       gnh3  = 0.0d0
       ghno3 = 0.0d0
       ghcl  = 0.0d0
       h     = 0.0d0
       lwn   = tiny
       so4_t = 0.0d0
       nh4_t = 0.0d0
       no3_t = 0.0d0
       na_t  = 0.0d0
       cl_t  = 0.0d0
       ca_t  = 0.0d0
       pk_t  = 0.0d0
       mg_t  = 0.0d0
       delcl = 0.0d0
       na2so4= 0.0d0
       nahso4= 0.0d0
       clc   = 0.0d0
       nh42s4= 0.0d0
       nh4hs4= 0.0d0
       kkhso4= 0.0d0
       k2so4 = 0.0d0
       mgso4 = 0.0d0
       frmg  = 0.0d0
       frna  = 0.0d0
       frk   = 0.0d0
       frca  = 0.0d0
       frcl  = 0.0d0
       frno3 = 0.0d0
       frso4 = 0.0d0
       frnh4 = 0.0d0
       exnh4 = 0.0d0
       gama  = 0.1d0
       gamin = 1.0d10
    !
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
                 (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    !
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3 = k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. HCl(g) <==> H+(aq) + Cl-(aq)                               (xk3)
       khcl = k0(11) * exp(p1(11)*tt1 + p2(11)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Find dry composition: all solids except CaSO4 assumed to be completely dissolved
    !  ## Subroutine 'CALCL1A' in ISORROPIA II, but modified here for mass balance
       caso4  = min(ca, so4)                     
       frso4  = max(so4 - caso4, 0.0d0)           
       frca   = max(ca  - caso4, 0.0d0)
    !
       k2so4  = min(0.5d0*pk, frso4)    
       frk    = max(pk - 2.0d0*k2so4, 0.0d0)      
       frso4  = max(frso4 - k2so4, 0.0d0)
    !             
       na2so4 = min(0.5d0*na, frso4)                  
       frna   = max(na - 2.0d0*na2so4, 0.0d0)      
       frso4  = max(frso4 - na2so4, 0.0d0)      
    !
       mgso4  = min(mg, frso4)                        
       frmg   = max(mg - mgso4, 0.0d0)            
       frso4  = max(frso4 - mgso4, 0.0d0)     
    !
       clc    = min(nh4/3.0d0, frso4*0.5d0)             ! (NH4)3H(SO4)2(s)
       frso4  = max(frso4  - 2.0d0*clc, 0.0d0)        
       frnh4  = max(nh4 - 3.0d0*clc, 0.0d0)        
    !
       if (frso4 <= tiny) then
          frso4   = frso4 + 2.0d0*clc - 2.0d0*max(clc-frnh4, 0.0d0)
          frnh4_2 = frnh4 + 3.0d0*clc - 3.0d0*max(clc-frnh4, 0.0d0)
          clc     = max(clc - frnh4, 0.0d0)          
          nh42s4  = min(0.5d0*frnh4_2, frso4)                         
          frnh4   = max(frnh4_2 - 2.0d0*nh42s4, 0.0d0)                
          frso4   = max(frso4 - nh42s4, 0.0d0)                        
       else if (frnh4 <= tiny) then
          nh4hs4  = 3.0d0*min(frso4, clc)             
          clc     = max(clc - frso4, 0.0d0)  
          frso4   = max(frso4 - nh4hs4/3.0d0, 0.0d0)  
    !
          if (na2so4 > tiny) then
             frso4_2 = frso4
             frso4   = frso4 + (na2so4 - max(na2so4 - frso4, 0.0d0))
             na2so4  = max(na2so4 - frso4_2, 0.0d0)    
             frna    = max(na - 2.0d0*na2so4, 0.0d0)
             nahso4  = min(frna, frso4)
             frna    = max(frna - nahso4, 0.0d0)
             frso4   = max(frso4 - nahso4, 0.0d0)
          end if
    !
          if (k2so4 > tiny) then
             frso4_2 = frso4
             frso4   = frso4 + k2so4 - max(k2so4 - frso4, 0.0d0)
             k2so4   = max(k2so4 - frso4_2, 0.0d0)
             frk     = max(pk - 2.0d0*k2so4, 0.0d0)
             kkhso4  = min(frk, frso4)
             frso4   = max(frso4 - kkhso4, 0.0d0)
             frk     = max(frk - kkhso4, 0.0d0)
          end if
       end if
    !
    !  ## Speciation
       na_t  = 2.0d0*na2so4 + nahso4
       nh4_t = 3.0d0*clc + 2.0d0*nh42s4 + nh4hs4
       pk_t  = kkhso4 + 2.0d0*k2so4
       mg_t  = mgso4
    !
    !  ## Gaseous species
       ghno3 = no3
       ghcl  = cl
    !
    !  ## Constant values
       c1 = clc + na2so4 + nh42s4 + k2so4 + mgso4
       c2 = kkhso4 + nh4hs4 + clc + nahso4
    !
    !
    !  ### MAJOR SYSTEM H+/HSO4-/SO42- ###
    !  ## Setup initial conditions
    !  ## 1. Calculate dissociation quantities
       a9 = khso4*1.0d-19
       bb = c1 + a9               ! bb always > 0
       cc = -a9*c2
    !
       if (bb /= 0.d0) then
          dd = cc/(bb*bb)
          v  = 4.d0*dd
       else
          v  = 1.d3
       end if
    !
       if (abs(v) <= smrt .and. bb > 0.d0) then
          hh = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
       else
          hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
       end if
    !
    !  ## 2. Speciation
       hh    = max(hh, tiny)   ! Avoid negative hh
       hh    = min(hh, c2)     ! Avoid negative HSO4
       h     = hh
       so4_t = c1 + hh
       hso4  = max(c2 - hh, 0.0d0)
    !
    !  ## 3. Aerosol liquid water content
       lwn = max(nh42s4/awas(irh) + na2so4/awss(irh) + nh4hs4/awab(irh) +    &
                 nahso4/awsb(irh) + clc/awlc(irh)    + k2so4/awps(irh)  +    &
                 kkhso4/awpb(irh) + mgso4/awms(irh), tiny)
    !
    !  ## Iterative search for solution with convergence of activity coefficients
       errin = 1.0d0
       j     = 0
       do while (j < nsweep-1 .and. errin >= epsact)
          j = j + 1
    !
    !  ## Reset gamin
          gamin(1)  = gama(1)
          gamin(2)  = gama(2)
          gamin(3)  = gama(3)
          gamin(4)  = gama(4)
          gamin(5)  = gama(5)
          gamin(6)  = gama(6)
          gamin(7)  = gama(7)
          gamin(8)  = gama(8)
          gamin(9)  = gama(9)
          gamin(10) = gama(10)
          gamin(11) = gama(11)
          gamin(12) = gama(12)
          gamin(13) = gama(13)
          gamin(14) = gama(14)
          gamin(15) = gama(15)
          gamin(16) = gama(16)
          gamin(17) = gama(17)
          gamin(18) = gama(18)
          gamin(19) = gama(19)
          gamin(20) = gama(20)
          gamin(21) = gama(21)
          gamin(22) = gama(22)
          gamin(23) = gama(23)
    !
          call mach_hetp_calcact4(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                  ca_t, pk_t, mg_t, lwn, gama, t)
    !
    !  ## Test for convergence of activity coefficients 
          errin = 0.0d0
          do ii = 1, 23
             errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
          end do
    !
    !  ## Solve system with new set of activity coefficients
    !  ## 1. Calculate dissociation quantities
          a9 = khso4*lwn/gama(7)*(gama(8)/gama(7))**2.0
          bb = c1 + a9        ! bb always > 0
          cc = -a9*c2
    !
          if (bb /= 0.d0) then
             dd = cc/(bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d3
          end if
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             hh = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc, 0.0d0)))
          end if
    !
    !  ## 2. Speciation
          hh    = max(hh, tiny)   ! Avoid negative hh
          hh    = min(hh, c2)     ! Avoid negative HSO4
          h     = hh
          so4_t = c1 + hh
          hso4  = max(c2 - hh, 0.0d0)
       end do
    !
    !
    !  ### MINOR SYSTEM #1: Cl-/HCl/NO3-/HNO3/H+ ###
    !  ## Calculate dissolution of HCl, HNO3 in the presence of (H,SO4)
    !  ## HCl, HNO3 are considered minor species 
       ispoly = .false.
    !
    !  ## 1. Special case: aerosol liquid water content (lwn) = 0.0d0
       if (lwn <= tiny) then
    !  Gaseous species
          ghcl  = max(cl  - cl_t,  0.0d0)
          ghno3 = max(no3 - no3_t, 0.0d0) 
          m1    = 0.0d0
          m2    = 0.0d0
          m3    = 0.0d0
    !
    !  ## 2. Special case: HCl = HNO3 = 0.0d0
       else if (cl <= tiny .and. no3 <= tiny) then
          m1  = 0.0d0
          m2  = 0.0d0
          m3  = 0.0d0
    !
    !  ## 3. Special case: HCl = 0.0d0
       else if (cl <= tiny) then
    !  Nitric acid in the liquid phase is assumed a minor species
    !  HNO3(g) <--> (H+) + (NO3-), using (H+) from the sulfates
          ff = 0.0d0
          m1 = 0.0d0
          m2 = 0.0d0
          m3 = 0.0d0
          if (lwn > tiny) then 
             bb = h + khno3*r*t*(lwn/gama(10))**2.0
             cc = (khno3*r*t*(lwn/gama(10))**2.0)*no3
    !
             if (bb /= 0.d0) then
                dd = cc / (bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d+03
             end if       
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                ff = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
             else
                ff = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
             end if 
          end if 
    !
    !  Speciation  
          ghno3 = max(no3 - ff, 0.0d0)
          no3_t = ff
          h     = h + ff
    !
    !  ## 4. Special case: HNO3 = 0.0d0
       else if (no3 <= tiny) then
    !  Hydrochloric acid in the liquid phase is assumed a minor species
    !  HCl (g) <--> (H+) + (Cl-) using (H+) from the sulfates
          m1 = 0.0d0
          m2 = 0.0d0
          m3 = 0.0d0
          ff = 0.0d0
          if (lwn > tiny) then
             bb = h + khcl*r*t*(lwn/gama(11))**2.0
             cc = (khcl*r*t*(lwn/gama(11))**2.0)*cl
    !
             if (bb /= 0.d0) then
                dd = cc / (bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d+03
             end if       
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                ff = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
             else
                ff = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
             end if 
          end if 
    !
    !  Speciation        
          ghcl = max(cl - ff, 0.0d0) 
          cl_t = ff 
          h    = h + ff
    !
    !  ## 5. All else: not a special case
       else
          ispoly = .true.
          a3 = khno3*r*t*(lwn/gama(10))**2.0  !HNO3
          a4 = khcl*r*t*(lwn/gama(11))**2.0   !HCl
    !
    !  Calculate cubic equation coefficients
          m1 = (a3*no3 + a4*cl + (h + a4)*(a3 - a4))/(a3 - a4)
          m2 = ((h + a4)*a4*cl - a4*(a3 - a4)*cl)/(a3 - a4)
          m3 = -a4*a4*cl*cl/(a3 - a4)
       end if 
    !
    !  ## Calculate roots 
    !   call mach_hetp_poly3v(m1, m2, m3, delcl, islv)       # Original ISORROPIA subroutine
       call mach_hetp_poly(m1, m2, m3, cl, delcl, islv)
    !   
       delno = 0.0d0
       if (ispoly) then
          if (islv /= 0) then
    !  Tiny amounts of HCl are assumed when there is no root
             delcl = tiny                                          !Change in Cl-
          end if
    !
          a3    = khno3*r*t*(lwn/gama(10))**2.0                    !HNO3
          a4    = khcl*r*t*(lwn/gama(11))**2.0                     !HCl
          delno = min(a3*no3*delcl/(a4*cl + (a3 - a4)*delcl), no3) !Change in NO3-
    !
          if (delcl < 0.0d0 .or. delno < 0.0d0 .or. delcl > cl .or. delno > no3) then
             delcl = tiny                                          !Change in Cl-
             delno = tiny                                          !Change in NO3-
          end if 
    !
    !  ## Effect on liquid phase
          h     = h + delno + delcl       ! H+   change
          cl_t  = cl_t  + delcl           ! Cl-  change
          no3_t = no3_t + delno           ! NO3- change
    !
    !  ## Gaseous species
          ghcl  = max(cl  - cl_t,  0.0d0)
          ghno3 = max(no3 - no3_t, 0.0d0) 
       end if 
    !
    !
    !  ### MINOR SYSTEM #2: NH4+/NH3/H+ ###
    !  Ammonia in the gas phase is assumed a minor species that does not significantly perturb 
    !  the aerosol equilibrium: NH3(g) + H+(aq) <==> NH4+(aq); H+ determined from the major
    !  system (solved above) is used initially. 
       if (lwn > tiny) then
    !  Calculate NH3 sublimation
          ff = 0.0d0
          a3 = (knh3/kh2o)*r*t*(gama(10)/gama(5))**2.0
          bb = h + 1.0d0/a3    ! bb always > 0
          cc = -nh4_t/a3
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d+03
          end if  
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             ff = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             ff = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
          end if 
    !
    !  Speciation 
    !  Due to round-off, ff may be .gt. nh4_t giving negative nh4_t
    !  If this condition is true, then set ff = nh4_t
          ff    = max(tiny2, min(ff, nh4_t))
          gnh3  = ff
          nh4_t = max(nh4_t - ff, 0.0d0)
          h     = h + ff
       end if
    !
    !  Add any free nh4 back to the gas phase
       gnh3 = gnh3 + frnh4
       frnh4 = 0.0d0
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       na_i    = na_t
       cl_i    = cl_t
       ca_i    = ca_t
       k_i     = pk_t
       mg_i    = mg_t
       nh3g_i  = gnh3
       hno3g_i = ghno3
       hclg_i  = ghcl
       h_i     = h
       lwn_i   = lwn
       caso4_i = caso4
    !
       return
    end subroutine mach_hetp_calcl9
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Subcase: K4; sulfate rich; free acid
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calck4(so4_i, nh4_i, nh3g_i, hno3g_i, hclg_i, hso4_i,          &
                                na_i, cl_i, no3_i, h_i, lwn_i, ca_i, k_i, mg_i,         &
                                caso4_i, rh, temp, k0, p1, p2, nr)
    !
       use mach_hetp_mod
       implicit none
    !
       integer(kind=4), intent   (in) :: nr
       real(kind=8),    intent   (in) :: k0      (nr)
       real(kind=8),    intent   (in) :: p1      (nr)
       real(kind=8),    intent   (in) :: p2      (nr)
       real(kind=8),    intent(inout) :: so4_i   
       real(kind=8),    intent(inout) :: nh4_i   
       real(kind=8),    intent(inout) :: no3_i   
       real(kind=8),    intent(inout) :: hso4_i 
       real(kind=8),    intent(inout) :: na_i    
       real(kind=8),    intent(inout) :: cl_i    
       real(kind=8),    intent(inout) :: ca_i    
       real(kind=8),    intent(inout) :: k_i     
       real(kind=8),    intent(inout) :: mg_i    
       real(kind=8),    intent(inout) :: nh3g_i  
       real(kind=8),    intent(inout) :: hno3g_i 
       real(kind=8),    intent(inout) :: hclg_i  
       real(kind=8),    intent(inout) :: h_i     
       real(kind=8),    intent(inout) :: lwn_i   
       real(kind=8),    intent(inout) :: caso4_i 
       real(kind=8),    intent   (in) :: rh      
       real(kind=8),    intent   (in) :: temp    
    ! 
    !  ## Local variables
       real(kind=8)      :: so4, nh4, hso4, gnh3, h, lwn, caso4
       real(kind=8)      :: no3, cl, na, ca, pk, mg, ghno3, ghcl, t, aw
       real(kind=8)      :: knh3, khso4, kh2o, khno3, khcl
       real(kind=8)      :: bb, cc, dd, ff, v, hh, errin, a3, a4, delno, tt0
       real(kind=8)      :: delcl, m1, m2, m3, c1, c2, c3, c4, c5, c6, c7, tt1, tt2
       real(kind=8)      :: so4_t, nh4_t, no3_t, na_t, cl_t, ca_t, pk_t, mg_t
       integer(kind=4)   :: j, ii, islv, irh
       logical(kind=4)   :: ispoly
       real(kind=8), dimension(23) :: gama, gamin
    !
    !
    !  ### Initialize variables ### 
       so4   = so4_i
       nh4   = nh4_i
       no3   = no3_i
       na    = na_i
       cl    = cl_i
       ca    = ca_i
       pk    = k_i
       mg    = mg_i
       aw    = rh
       t     = temp
       hso4  = 0.0d0
       gnh3  = 0.0d0
       ghno3 = 0.0d0
       ghcl  = 0.0d0
       h     = 0.0d0
       lwn   = tiny
       so4_t = 0.0d0
       nh4_t = 0.0d0
       no3_t = 0.0d0
       na_t  = 0.0d0
       cl_t  = 0.0d0
       ca_t  = 0.0d0
       pk_t  = 0.0d0
       mg_t  = 0.0d0
       caso4 = 0.0d0
       delcl = 0.0d0
       gama  = 0.1d0
       gamin = 1.0d10
    !
    !
    !  ### Calculate equilibrium constants and other static variables ###
    !  ## Set RH to a range between 0.5% and 99.5%: RH = 0.00d0 will 
    !  ## cause division by zero, aborting the code
       aw = max(aw, 0.005d0)
       aw = min(aw, 0.995d0)
    !
       tt0 = tstd / t
       tt1 = tt0 - 1.0d0
       tt2 = 1.0d0 + log(tt0) - tt0
    !
    !  ## 1. HSO4(aq) <==> H+(aq) + SO4=(aq)                            (xk1)
       khso4 = k0(5) * exp(p1(5)*tt1 + p2(5)*tt2)
    !
    !  ## 2. k2 = NH3(g) <==> NH3(aq)                                   (xk21)
    !  ## 3. k3 = NH3(aq) + H2O(aq) <==> NH4+(aq) + OH-(aq)             (xk22)
    !  ## Net NH3: k2*k3                                                (xk2)
       knh3 = (k0(3) * exp(p1(3)*tt1 + p2(3)*tt2))*                  &
              (k0(4) * exp(p1(4)*tt1 + p2(4)*tt2))
    ! 
    !  ## 4. H2O(aq) <==> H+(aq) + OH-(aq)                              (xkw)
       kh2o = k0(10) * exp(p1(10)*tt1 + p2(10)*tt2)
    !   
    !  ## 5. HNO3(g) <==> H+(aq) + NO3-(aq)                             (xk4)
       khno3 = k0(15) * exp(p1(15)*tt1 + p2(15)*tt2)
    !
    !  ## 6. HCl(g) <==> H+(aq) + Cl-(aq)                               (xk3)
       khcl = k0(11) * exp(p1(11)*tt1 + p2(11)*tt2)
    !
    !  ## Calculate ZSR position parameter
       irh = max(min(int(aw*100+0.5), 100), 1)
    !
    !  ## Constants
       c1 = max(so4 - nh4 - na - ca - pk - mg, tiny)  ! Free H2SO4
       c2 = c1 + pk + na + nh4
       c3 = c1 + mg
       c4 = c1*mg
    !
    !  ## ZSR constants
       c5 = nh4/awab(irh) + na/awsb(irh) + pk/awpb(irh) + mg/awms(irh)
       c6 = nh4 + na + pk + mg
       c7 = awsa(irh)
    !
    !
    !  ### MAJOR SYSTEM H+/HSO4-/SO42- ###
    !  ## Setup initial conditions
    !  ## 1. Calculate dissociation quantities 
       a4 = khso4*1.0d-19
       bb = a4 + c3             ! bb always > 0
       cc = -a4*c2 + c4
    !
       if (bb /= 0.d0) then
          dd = cc/(bb*bb)
          v  = 4.d0*dd
       else
          v  = 1.d3
       end if
    !
       if (abs(v) <= smrt .and. bb > 0.d0) then
          hh = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
       else
          hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc,0.0d0)))
       end if 
       hh = min(c2, hh)
    !
    !  ## 2. Speciation 
       h     = c1 + hh
       na_t  = na
       nh4_t = nh4
       so4_t = hh + mg
       hso4  = max(c2 - hh, 0.0d0)
       pk_t  = pk
       mg_t  = mg
       caso4 = ca
    !
    !  ## 3. Aerosol liquid water content
       lwn = max(c5 + max(so4_t + hso4 - c6, 0.0d0)/c7, tiny)
    !
    !  ## Iterative search for solution with convergence of activity coefficients
       errin = 1.0d0
       j     = 0
       do while (j < nsweep-1 .and. errin >= epsact)
          j = j + 1
    !
    !  ## Reset gamin
          gamin(1)  = gama(1)
          gamin(2)  = gama(2)
          gamin(3)  = gama(3)
          gamin(4)  = gama(4)
          gamin(5)  = gama(5)
          gamin(6)  = gama(6)
          gamin(7)  = gama(7)
          gamin(8)  = gama(8)
          gamin(9)  = gama(9)
          gamin(10) = gama(10)
          gamin(11) = gama(11)
          gamin(12) = gama(12)
          gamin(13) = gama(13)
          gamin(14) = gama(14)
          gamin(15) = gama(15)
          gamin(16) = gama(16)
          gamin(17) = gama(17)
          gamin(18) = gama(18)
          gamin(19) = gama(19)
          gamin(20) = gama(20)
          gamin(21) = gama(21)
          gamin(22) = gama(22)
          gamin(23) = gama(23)
    !
          call mach_hetp_calcact4(h, nh4_t, so4_t, hso4, no3_t, cl_t, na_t,     &
                                  ca_t, pk_t, mg_t, lwn, gama, t)
    !
    !  ## Test for convergence of activity coefficients 
          errin = 0.0d0
          do ii = 1, 23
             errin = max(max(errin, abs((gamin(ii) - gama(ii)) / gamin(ii))), 0.0d0)
          end do   
    !
    !  ## Solve system with new set of activity coefficients
    !  ## 1. Calculate dissociation quantities 
          a4 = khso4*lwn/gama(7)*(gama(8)/gama(7))**2.0
          bb = a4 + c3         
          cc = -a4*c2 + c4
    !
          if (bb /= 0.d0) then
             dd = cc/(bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d3
          end if
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             hh = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             hh = 0.5d0*(-bb + sqrt(max(bb*bb - 4.d0*cc,0.0d0)))
          end if 
          hh = min(c2, hh)
    !
    !  ## 2. Speciation 
          h     = c1 + hh
          so4_t = hh + mg
          hso4  = max(c2 - hh, 0.0d0)
    !
    !  ## 3. Aerosol liquid water content
          lwn = max(c5 + max(so4_t + hso4 - c6, 0.0d0)/c7, tiny)
       end do
    !
    !
    !  ### MINOR SYSTEM #1: Cl-/HCl/NO3-/HNO3/H+ ###
    !  ## Calculate dissolution of HCl, HNO3 in the presence of (H,SO4)
    !  ## HCl, HNO3 are considered minor species 
       ispoly = .false.
    !
    !  ## 1. Special case: aerosol liquid water content (lwn) = 0.0d0
       if (lwn <= tiny) then
    !  Gaseous species
          ghcl  = max(cl  - cl_t,  0.0d0)
          ghno3 = max(no3 - no3_t, 0.0d0) 
          m1    = 0.0d0
          m2    = 0.0d0
          m3    = 0.0d0
    !
    !  ## 2. Special case: HCl = HNO3 = 0.0d0
       else if (cl <= tiny .and. no3 <= tiny) then
          m1  = 0.0d0
          m2  = 0.0d0
          m3  = 0.0d0
    !
    !  ## 3. Special case: HCl = 0.0d0
       else if (cl <= tiny) then
    !  Nitric acid in the liquid phase is assumed a minor species
    !  HNO3(g) <--> (H+) + (NO3-), using (H+) from the sulfates
          ff = 0.0d0
          m1 = 0.0d0
          m2 = 0.0d0
          m3 = 0.0d0
          if (lwn > tiny) then 
             bb = h + khno3*r*t*(lwn/gama(10))**2.0
             cc = (khno3*r*t*(lwn/gama(10))**2.0)*no3
    !
             if (bb /= 0.d0) then
                dd = cc / (bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d+03
             end if       
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                ff = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
             else
                ff = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
             end if 
          end if 
    !
    !  Speciation  
          ghno3 = max(no3 - ff, 0.0d0)
          no3_t = ff
          h     = h + ff
    !
    !  ## 4. Special case: HNO3 = 0.0d0
       else if (no3 <= tiny) then
    !  Hydrochloric acid in the liquid phase is assumed a minor species
    !  HCl (g) <--> (H+) + (Cl-) using (H+) from the sulfates
          m1 = 0.0d0
          m2 = 0.0d0
          m3 = 0.0d0
          ff = 0.0d0
          if (lwn > tiny) then
             bb = h + khcl*r*t*(lwn/gama(11))**2.0
             cc = (khcl*r*t*(lwn/gama(11))**2.0)*cl
    !
             if (bb /= 0.d0) then
                dd = cc / (bb*bb)
                v  = 4.d0*dd
             else
                v  = 1.d+03
             end if       
    !
             if (abs(v) <= smrt .and. bb /= 0.d0) then
                ff = -0.5d0*bb + 0.5d0*abs(bb) + (1.d0-dd*(1.d0-dd*(2.d0-dd*(5.d0-14.d0*dd))))*cc/abs(bb)
             else
                ff = 0.5d0*(-bb + sqrt(bb*bb + 4.d0*cc))   
             end if 
          end if 
    !
    !  Speciation        
          ghcl = max(cl - ff, 0.0d0) 
          cl_t = ff 
          h    = h + ff
    !
    !  ## 5. All else: not a special case
       else
          ispoly = .true.
          a3 = khno3*r*t*(lwn/gama(10))**2.0  !HNO3
          a4 = khcl*r*t*(lwn/gama(11))**2.0   !HCl
    !
    !  Calculate cubic equation coefficients
          m1 = (a3*no3 + a4*cl + (h + a4)*(a3 - a4))/(a3 - a4)
          m2 = ((h + a4)*a4*cl - a4*(a3 - a4)*cl)/(a3 - a4)
          m3 = -a4*a4*cl*cl/(a3 - a4)
       end if 
    !
    !  ## Calculate roots 
       call mach_hetp_poly(m1, m2, m3, cl, delcl, islv)
    !   call mach_hetp_poly3v(m1, m2, m3, delcl, islv)      # Original ISORROPIA subroutine 
    !   
       delno = 0.0d0
       if (ispoly) then
          if (islv /= 0) then
    !  Tiny amounts of HCl are assumed when there is no root
             delcl = tiny                                          !Change in Cl-
          end if
    !
          a3    = khno3*r*t*(lwn/gama(10))**2.0                    !HNO3
          a4    = khcl*r*t*(lwn/gama(11))**2.0                     !HCl
          delno = min(a3*no3*delcl/(a4*cl + (a3 - a4)*delcl), no3) !Change in NO3-
    !
          if (delcl < 0.0d0 .or. delno < 0.0d0 .or. delcl > cl .or. delno > no3) then
             delcl = tiny                                         !Change in Cl-
             delno = tiny                                         !Change in NO3-
          end if 
    !
    !  ## Effect on liquid phase
          h     = h + delno + delcl       ! H+   change
          cl_t  = cl_t  + delcl           ! Cl-  change
          no3_t = no3_t + delno           ! NO3- change
    !
    !  ## Gaseous species
          ghcl  = max(cl  - cl_t,  0.0d0)
          ghno3 = max(no3 - no3_t, 0.0d0) 
       end if 
    !
    !
    !  ### MINOR SYSTEM #2: NH4+/NH3/H+ ###
    !  Ammonia in the gas phase is assumed a minor species that does not significantly perturb 
    !  the aerosol equilibrium: NH3(g) + H+(aq) <==> NH4+(aq); H+ determined from the major
    !  system (solved above) is used initially.  
       if (lwn > tiny) then
    !  Calculate NH3 sublimation
          ff = 0.0d0
          a3 = (knh3/kh2o)*r*t*(gama(10)/gama(5))**2.0
          bb = h + 1.0d0/a3    ! bb always > 0
          cc = -nh4_t/a3
    !
          if (bb /= 0.d0) then
             dd = cc / (bb*bb)
             v  = 4.d0*dd
          else
             v  = 1.d+03
          end if  
    !
          if (abs(v) <= smrt .and. bb > 0.d0) then
             ff = - ((((14.d0*dd + 5.d0)*dd + 2.d0)*dd + 1.d0)*dd + 1.d0)*cc/bb
          else
             ff = 0.5d0*(-bb + sqrt(bb*bb - 4.d0*cc))
          end if 
    !
    !  Speciation 
    !  Due to round-off, ff may be .gt. nh4_t giving negative nh4_t
    !  If this condition is true, then set ff = nh4_t
          ff    = max(tiny, min(ff, nh4_t))
          gnh3  = ff
          nh4_t = max(nh4_t - ff, 0.0d0)
          h     = h + ff
       end if
    !
    !
    !  ### Perform mass adjustment if excess exists ###
       call mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                             nh4_t, gnh3, cl_t, ghcl, caso4)
    !
    !
    !  ### Save result and return ###
       so4_i   = so4_t 
       nh4_i   = nh4_t
       no3_i   = no3_t
       hso4_i  = hso4
       na_i    = na_t
       cl_i    = cl_t
       ca_i    = ca_t
       k_i     = pk_t
       mg_i    = mg_t
       nh3g_i  = gnh3
       hno3g_i = ghno3
       hclg_i  = ghcl
       h_i     = h
       lwn_i   = lwn
       caso4_i = caso4
    !
       return
    end subroutine mach_hetp_calck4
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Calculates multi-component activity coefficients from Bromley's method
    ! ## of ammonium-sulfate aerosol system for case B4, C2
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcact1(h, nh4_t, so4_t, hso4, lwn, gama, t)
    !
       use mach_hetp_mod,         only: tiny
       implicit none
    !
       real(kind=8),       intent(in) :: nh4_t   
       real(kind=8),       intent(in) :: so4_t  
       real(kind=8),       intent(in) :: hso4    
       real(kind=8),       intent(in) :: lwn    
       real(kind=8),       intent(in) :: h       
       real(kind=8),       intent(in) :: t       
       real(kind=8),    intent(inout) :: gama    (13)
    !
    !  ## Local variables
       real(kind=8)    :: ionic, sion, tc, c, xx, xx2, c1, c2, hh
       real(kind=8)    :: c3, c4, c5, c5a
       real(kind=8)    :: ff11, ff13, f2a2, f2a3, ch1, ch2
       real(kind=8)    :: h2, k1, k2, f11, f12
       real(kind=8)    :: g04, g06, g07, g08, g09, g11
    !
    !  ## Calculate ionic strength of solution
       ionic = h + nh4_t + 4.0d0*so4_t + hso4
       ionic = max(min(0.5d0*ionic/lwn, 100.d0), tiny)
    !
    !  ## Calculate the binary activity coefficients
       sion = sqrt(ionic)
       c3   = exp(-0.023d0*ionic*ionic*ionic)
       c4   = -0.5107d0*sion
       c5   = 1.0d0 + 0.1d0*ionic
       c5a  = 1.0d0 + sion
    !
    !  ## Coefficients at 25C
       if (ionic < 6.0d0) then
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + 2.0d0*xx
    
          c   = 1.0d0 + -0.0055d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + 2.0d0*xx
    
          c   = 1.0d0 + 0.44d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
    
          c   = 1.0d0 + 0.0451d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
    
          c   = 1.0d0 + 0.33d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       else
          xx  = c4/(1.0d0 + sion)
          xx2 = 2.0d0*xx
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + xx2
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + xx2
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       end if 
    !
    !  ## Correct for temperature other than 298 K
       tc = abs(t - 298.0d0)
       if (abs(tc) > 1.0) then
          c1  = 1.125d0 - 0.005d0*(t - 273.0d0)
          c2  = (0.125d0 - 0.005d0*(t - 273.0d0))*(0.039d0*ionic**0.92 -   &
                 0.41d0*sion/c5a)
          c3  = 2.0d0*c2
          g04 = c1*g04 - c3                 
          g06 = c1*g06 - c2         
          g07 = c1*g07 - c3         
          g08 = c1*g08 - c2                 
          g11 = c1*g11 - c2    
       end if
    !
    !  ## Correction: g9 is g09 which is not calculated in calcact1
    !                 use g09 from calcact3 to represent g09 (slc.2.2012)
       g09 = g06 + g08 - g11 
    !
    !  ## Calculate multicomponent activity coefficients
       hh   = (0.511d0*(298.0d0/t)**1.5)*sion/c5a
       f11  = so4_t/lwn
       f12  = hso4/lwn
       ch1  = 2.25d0/ionic
       ch2  = 1.0d0/ionic
       h2   = 2.0d0*hh
       k1   = ch1*f11
       k2   = ch2*f12
       ff11 = k1*(g07 + h2) + k2*(g08 + hh) 
       ff13 = k1*(g04 + h2) + k2*(g09 + hh)  
    
       f11  = h/lwn
       f12  = nh4_t/lwn
       f2a2 = ((ch1*f11)*(g07 + h2) + (ch1*f12)*(g04 + h2))*0.5d0 
       f2a3 = ch2*f11*(g08 + hh) + ch2*f12*(g09 + hh)
    
    !  ## log10 of activity coefficients
       gama(4)  = ((ff13 + f2a2) / 3.0d0 - hh)*2.0d0  ! (NH4)2SO4
       gama(7)  = ((ff11 + f2a2) / 3.0d0 - hh)*2.0d0  ! 2H-SO4
       gama(8)  = ((ff11 + f2a3) * 0.5d0 - hh)        ! H-HSO4
       gama(9)  = ((ff13 + f2a3) * 0.5d0 - hh)        ! NH4HSO4
       gama(13) = 0.6d0*gama(4) + 0.4d0*gama(9)      ! lc; scape
    !
    !  ## Convert log(gama) coefficients to gama
       gama(4)  = max(-5.0d0, min(gama(4), 5.0d0))
       gama(4)  = 10.0d0**gama(4)
       gama(13) = max(-5.0d0, min(gama(13), 5.0d0))
       gama(13) = 10.0d0**gama(13)
       gama(7)  = max(-5.0d0, min(gama(7), 5.0d0))
       gama(7)  = 10.0d0**gama(7)
       gama(8)  = max(-5.0d0, min(gama(8), 5.0d0))
       gama(8)  = 10.0d0**gama(8)
       gama(9)  = max(-5.0d0, min(gama(9), 5.0d0))
       gama(9)  = 10.0d0**gama(9)
    !
       return
    end subroutine mach_hetp_calcact1
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Calculates multi-component activity coefficients from Bromley's method
    ! ## of ammonium-sulfate aerosol system for case A2
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcact1b(h, nh4_t, so4_t, hso4, lwn, gama, t, soln,  &
                                   frst, calain, calou)
    !
       use mach_hetp_mod,         only: tiny
       implicit none
    !
       real(kind=8),       intent(in) :: nh4_t   
       real(kind=8),       intent(in) :: so4_t  
       real(kind=8),       intent(in) :: hso4    
       real(kind=8),       intent(in) :: lwn    
       real(kind=8),       intent(in) :: h       
       real(kind=8),       intent(in) :: t       
       real(kind=8),    intent(inout) :: gama    (13)
       logical(kind=4),    intent(in) :: soln   
       logical(kind=4),    intent(in) :: frst    
       logical(kind=4),    intent(in) :: calain  
       logical(kind=4),    intent(in) :: calou   
    !
    !  ## Local variables
       real(kind=8)    :: ionic, sion, tc, c, xx, xx2, c1, c2, hh
       real(kind=8)    :: c3, c4, c5, c5a
       real(kind=8)    :: ff11, ff13, f2a2, f2a3, ch1, ch2
       real(kind=8)    :: h2, k1, k2, f11, f12
       real(kind=8)    :: g04, g06, g07, g08, g09, g11
    !
      if ((( .not. soln)).and.  &
         ((frst  .and. calou)  .or.   &
          (( .not. frst) .and. calain))) then
    !  ## Calculate ionic strength of solution
       ionic = h + nh4_t + 4.0d0*so4_t + hso4
       ionic = max(min(0.5d0*ionic/lwn, 100.d0), tiny)
    !
    !  ## Calculate the binary activity coefficients
       sion = sqrt(ionic)
       c3   = exp(-0.023d0*ionic*ionic*ionic)
       c4   = -0.5107d0*sion
       c5   = 1.0d0 + 0.1d0*ionic
       c5a  = 1.0d0 + sion
    !
    !  ## Coefficients at 25C
       if (ionic < 6.0d0) then
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + 2.0d0*xx
    
          c   = 1.0d0 + -0.0055d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + 2.0d0*xx
    
          c   = 1.0d0 + 0.44d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
    
          c   = 1.0d0 + 0.0451d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
    
          c   = 1.0d0 + 0.33d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       else
          xx  = c4/(1.0d0 + sion)
          xx2 = 2.0d0*xx
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + xx2
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + xx2
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       end if 
    !
    !  ## Correct for temperature other than 298 K
       tc = abs(t - 298.0d0)
       if (abs(tc) > 1.0) then
          c1  = 1.125d0 - 0.005d0*(t - 273.0d0)
          c2  = (0.125d0 - 0.005d0*(t - 273.0d0))*(0.039d0*ionic**0.92 -   &
                 0.41d0*sion/c5a)
          c3  = 2.0d0*c2
          g04 = c1*g04 - c3                 
          g06 = c1*g06 - c2         
          g07 = c1*g07 - c3         
          g08 = c1*g08 - c2                 
          g11 = c1*g11 - c2    
       end if
    !
    !  ## Correction: g9 is g09 which is not calculated in calcact1
    !                 use g09 from calcact3 to represent g09 (slc.2.2012)
       g09 = g06 + g08 - g11 
    !
    !  ## Calculate multicomponent activity coefficients
       hh   = (0.511d0*(298.0d0/t)**1.5)*sion/c5a
       f11  = so4_t/lwn
       f12  = hso4/lwn
       ch1  = 2.25d0/ionic
       ch2  = 1.0d0/ionic
       h2   = 2.0d0*hh
       k1   = ch1*f11
       k2   = ch2*f12
       ff11 = k1*(g07 + h2) + k2*(g08 + hh) 
       ff13 = k1*(g04 + h2) + k2*(g09 + hh)  
    
       f11  = h/lwn
       f12  = nh4_t/lwn
       f2a2 = ((ch1*f11)*(g07 + h2) + (ch1*f12)*(g04 + h2))*0.5d0 
       f2a3 = ch2*f11*(g08 + hh) + ch2*f12*(g09 + hh)
    
    !  ## log10 of activity coefficients
       gama(4) = ((ff13 + f2a2) / 3.0d0 - hh)*2.0d0  ! (NH4)2SO4
       gama(7) = ((ff11 + f2a2) / 3.0d0 - hh)*2.0d0  ! 2H-SO4
       gama(8) = ((ff11 + f2a3) * 0.5d0 - hh)        ! H-HSO4
       gama(9) = ((ff13 + f2a3) * 0.5d0 - hh)        ! NH4HSO4
       gama(13)= 0.6d0*gama(4) + 0.4d0*gama(9)       ! lc; scape
    !
    !  ## Convert log(gama) coefficients to gama
       gama(4)  = max(-5.0d0, min(gama(4), 5.0d0))
       gama(4)  = 10.0d0**gama(4)
       gama(13) = max(-5.0d0, min(gama(13), 5.0d0))
       gama(13) = 10.0d0**gama(13)
       gama(7)  = max(-5.0d0, min(gama(7), 5.0d0))
       gama(7)  = 10.0d0**gama(7)
       gama(8)  = max(-5.0d0, min(gama(8), 5.0d0))
       gama(8)  = 10.0d0**gama(8)
       gama(9)  = max(-5.0d0, min(gama(9), 5.0d0))
       gama(9)  = 10.0d0**gama(9)
       end if 
    !
       return
    end subroutine mach_hetp_calcact1b
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Calculates multi-component activity coefficients from Bromley's method
    ! ## of ammonium-sulfate-nitrate aerosol system for case E4 and F2
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcact2(h, nh4_t, so4_t, hso4, no3, lwn, gama, t)
    !
       use mach_hetp_mod,         only: tiny
       implicit none
    !
       real(kind=8),       intent(in) :: nh4_t   
       real(kind=8),       intent(in) :: so4_t  
       real(kind=8),       intent(in) :: hso4    
       real(kind=8),       intent(in) :: no3     
       real(kind=8),       intent(in) :: lwn   
       real(kind=8),       intent(in) :: h       
       real(kind=8),       intent(in) :: t       
       real(kind=8),    intent(inout) :: gama(13)
    !
    !  ## Local variables:
       real(kind=8)    :: ionic, sion, tc, c, xx, xx2, c1, c2, hh
       real(kind=8)    :: c3, c4, c5, c5a
       real(kind=8)    :: ff11, ff13, f2a2, f2a3, f2a4, ch1, ch2
       real(kind=8)    :: h2, k1, k2, k3, f11, f12, f13
       real(kind=8)    :: g04, g05, g06, g07, g08, g09, g10, g11
    !
    !  ## Calculate ionic strength of solution
       ionic = h + nh4_t + 4.0d0*so4_t + hso4 + no3
       ionic = max(min(0.5d0*ionic/lwn, 100.d0), tiny)
    !
    !  ## Calculate the binary activity coefficients
       sion = sqrt(ionic)
       c3   = exp(-0.023d0*ionic*ionic*ionic)
       c4   = -0.5107d0*sion
       c5   = 1.0d0 + 0.1d0*ionic
       c5a  = 1.0d0 + sion
    !
    !  ## Coefficients at 25C
       if (ionic < 6.0d0) then
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + 2.0d0*xx
    
          c   = 1.0d0 + -0.06325d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
    
          c   = 1.0d0 + -0.0055d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + 2.0d0*xx
    
          c   = 1.0d0 + 0.44d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
    
          c   = 1.0d0 + 0.0451d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
    
          c   = 1.0d0 + 0.143d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
    
          c   = 1.0d0 + 0.33d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       else
          xx  = c4/(1.0d0 + sion)
          xx2 = 2.0d0*xx
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + xx2
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + xx2
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       end if 
    !
    !  ## Correct for temperature other than 298 K
       tc = abs(t - 298.0d0)
       if (tc > 1.0) then
          c1  = 1.125d0  - 0.005d0*(t - 273.0d0)
          c2  = (0.125d0 - 0.005d0*(t - 273.0d0))*(0.039d0*ionic**0.92 - 0.41d0*sion/c5a)
          c3  = 2.0d0*c2
          g04 = c1*g04 - c3         
          g05 = c1*g05 - c2         
          g06 = c1*g06 - c2         
          g07 = c1*g07 - c3         
          g08 = c1*g08 - c2         
          g10 = c1*g10 - c2         
          g11 = c1*g11 - c2         
       end if
    !
    !  ## Correction: g9 is g09 which is not calculated in calcact1
    !                 use g09 from calcact3 to represent g09 (slc.2.2012)
       g09 = g06 + g08 - g11 
    !
    !  ## Calculate multicomponent activity coefficients
       hh = (0.511*(298.0/t)**1.5)*sion/c5a
    !
       f11   = so4_t/lwn
       f12   = hso4/lwn
       f13   = no3/lwn
       ch1   = 2.25d0/ionic
       ch2   = 1.0d0/ionic
       h2    = 2.0d0*hh
       k1    = ch1*f11
       k2    = ch2*f12
       k3    = ch2*f13
       ff11  = k1*(g07 + h2) + k2*(g08 + hh) + k3*(g10 + hh)
       ff13  = k1*(g04 + h2) + k2*(g09 + hh) + k3*(g05 + hh)
    !
       f11   = h/lwn
       f12   = nh4_t/lwn
       k1    = ch2*f11
       k2    = ch2*f12
       f2a2  = ((ch1*f11)*(g07 + h2) + (ch1*f12)*(g04 + h2))*0.5d0 
       f2a3  = k1*(g08 + hh) + k2*(g09 + hh)
       f2a4  = k1*(g10 + hh) + k2*(g05 + hh)
    !
    !  ## log10 of activity coefficients
       gama(4)  = ((ff13 + f2a2) / 3.0d0 - hh)*2.0d0  ! (NH4)2SO4
       gama(5)  = ((ff13 + f2a4) * 0.5d0 - hh)        ! NH4NO3
       gama(7)  = ((ff11 + f2a2) / 3.0d0 - hh)*2.0d0  ! 2H-SO4
       gama(8)  = ((ff11 + f2a3) * 0.5d0 - hh)        ! H-HSO4
       gama(9)  = ((ff13 + f2a3) * 0.5d0 - hh)        ! NH4HSO4
       gama(10) = ((ff11 + f2a4) * 0.5d0 - hh)        ! HNO3
       gama(13) = 0.6d0*gama(4) + 0.4d0*gama(9)       ! lc; scape
    !
    !  ## Convert log(gama) coefficients to gama
       gama(4)  = max(-5.0d0, min(gama(4), 5.0d0))
       gama(4)  = 10.0d0**gama(4)
       gama(5)  = max(-5.0d0, min(gama(5), 5.0d0))
       gama(5)  = 10.0d0**gama(5)
       gama(13) = max(-5.0d0, min(gama(13), 5.0d0))
       gama(13) = 10.0d0**gama(13)
       gama(7)  = max(-5.0d0, min(gama(7), 5.0d0))
       gama(7)  = 10.0d0**gama(7)
       gama(8)  = max(-5.0d0, min(gama(8), 5.0d0))
       gama(8)  = 10.0d0**gama(8)
       gama(9)  = max(-5.0d0, min(gama(9), 5.0d0))
       gama(9)  = 10.0d0**gama(9)
       gama(10) = max(-5.0d0, min(gama(10), 5.0d0))
       gama(10) = 10.0d0**gama(10)
       return
    end subroutine mach_hetp_calcact2
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Calculates multi-component activity coefficients from Bromley's method
    ! ## of ammonium-sulfate-nitrate aerosol system for case D3
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcact2b(h, nh4_t, so4_t, hso4, no3, lwn, gama, t, soln,   &
                                   frst, calain, calou)
    !
       use mach_hetp_mod,         only: tiny
       implicit none
    !
       real(kind=8),       intent(in) :: nh4_t   
       real(kind=8),       intent(in) :: so4_t  
       real(kind=8),       intent(in) :: hso4    
       real(kind=8),       intent(in) :: no3     
       real(kind=8),       intent(in) :: lwn   
       real(kind=8),       intent(in) :: h       
       real(kind=8),       intent(in) :: t       
       real(kind=8),    intent(inout) :: gama(13)
       logical(kind=4),    intent(in) :: soln    
       logical(kind=4),    intent(in) :: frst    
       logical(kind=4),    intent(in) :: calain  
       logical(kind=4),    intent(in) :: calou   
    !
    !  ## Local variables:
       real(kind=8)    :: ionic, sion, tc, c, xx, xx2, c1, c2, hh
       real(kind=8)    :: c3, c4, c5, c5a
       real(kind=8)    :: ff11, ff13, f2a2, f2a3, f2a4, ch1, ch2
       real(kind=8)    :: h2, k1, k2, k3, f11, f12, f13
       real(kind=8)    :: g04, g05, g06, g07, g08, g09, g10, g11
    !
       if ((( .not. soln)) .and.  &
          ((frst  .and. calou) .or.    &
           (( .not. frst) .and. calain)) ) then
    !  ## Calculate ionic strength of solution
       ionic = h + nh4_t + 4.0d0*so4_t + hso4 + no3
       ionic = max(min(0.5d0*ionic/lwn, 100.d0), tiny)
    !
    !  ## Calculate the binary activity coefficients
       sion = sqrt(ionic)
       c3   = exp(-0.023d0*ionic*ionic*ionic)
       c4   = -0.5107d0*sion
       c5   = 1.0d0 + 0.1d0*ionic
       c5a  = 1.0d0 + sion
    !
    !  ## Coefficients at 25C
       if (ionic < 6.0d0) then
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + 2.0d0*xx
    
          c   = 1.0d0 + -0.06325d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
    
          c   = 1.0d0 + -0.0055d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + 2.0d0*xx
    
          c   = 1.0d0 + 0.44d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
    
          c   = 1.0d0 + 0.0451d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
    
          c   = 1.0d0 + 0.143d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
    
          c   = 1.0d0 + 0.33d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       else
          xx  = c4/(1.0d0 + sion)
          xx2 = 2.0d0*xx
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + xx2
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + xx2
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       end if 
    !
    !  ## Correct for temperature other than 298 K
       tc = abs(t - 298.0d0)
       if (tc > 1.0) then
          c1  = 1.125d0  - 0.005d0*(t - 273.0d0)
          c2  = (0.125d0 - 0.005d0*(t - 273.0d0))*(0.039d0*ionic**0.92 - 0.41d0*sion/c5a)
          c3  = 2.0d0*c2
          g04 = c1*g04 - c3         
          g05 = c1*g05 - c2         
          g06 = c1*g06 - c2         
          g07 = c1*g07 - c3         
          g08 = c1*g08 - c2         
          g10 = c1*g10 - c2         
          g11 = c1*g11 - c2         
       end if
    !
    !  ## Correction: g9 is g09 which is not calculated in calcact1
    !                 use g09 from calcact3 to represent g09 (slc.2.2012)
       g09 = g06 + g08 - g11 
    !
    !  ## Calculate multicomponent activity coefficients
       hh = (0.511*(298.0/t)**1.5)*sion/c5a
    !
       f11   = so4_t/lwn
       f12   = hso4/lwn
       f13   = no3/lwn
       ch1   = 2.25d0/ionic
       ch2   = 1.0d0/ionic
       h2    = 2.0d0*hh
       k1    = ch1*f11
       k2    = ch2*f12
       k3    = ch2*f13
       ff11  = k1*(g07 + h2) + k2*(g08 + hh) + k3*(g10 + hh)
       ff13  = k1*(g04 + h2) + k2*(g09 + hh) + k3*(g05 + hh)
    !
       f11   = h/lwn
       f12   = nh4_t/lwn
       k1    = ch2*f11
       k2    = ch2*f12
       f2a2  = ((ch1*f11)*(g07 + h2) + (ch1*f12)*(g04 + h2))*0.5d0 
       f2a3  = k1*(g08 + hh) + k2*(g09 + hh)
       f2a4  = k1*(g10 + hh) + k2*(g05 + hh)
    !
    !  ## log10 of activity coefficients
       gama(4)  = ((ff13 + f2a2) / 3.0d0 - hh)*2.0d0  ! (NH4)2SO4
       gama(5)  = ((ff13 + f2a4) * 0.5d0 - hh)        ! NH4NO3
       gama(7)  = ((ff11 + f2a2) / 3.0d0 - hh)*2.0d0  ! 2H-SO4
       gama(8)  = ((ff11 + f2a3) * 0.5d0 - hh)        ! H-HSO4
       gama(9)  = ((ff13 + f2a3) * 0.5d0 - hh)        ! NH4HSO4
       gama(10) = ((ff11 + f2a4) * 0.5d0 - hh)        ! HNO3
       gama(13) = 0.6d0*gama(4) + 0.4d0*gama(9)       ! lc; scape
    !
    !  ## Convert log(gama) coefficients to gama
       gama(4)  = max(-5.0d0, min(gama(4), 5.0d0))
       gama(4)  = 10.0d0**gama(4)
       gama(5)  = max(-5.0d0, min(gama(5), 5.0d0))
       gama(5)  = 10.0d0**gama(5)
       gama(13) = max(-5.0d0, min(gama(13), 5.0d0))
       gama(13) = 10.0d0**gama(13)
       gama(7)  = max(-5.0d0, min(gama(7), 5.0d0))
       gama(7)  = 10.0d0**gama(7)
       gama(8)  = max(-5.0d0, min(gama(8), 5.0d0))
       gama(8)  = 10.0d0**gama(8)
       gama(9)  = max(-5.0d0, min(gama(9), 5.0d0))
       gama(9)  = 10.0d0**gama(9)
       gama(10) = max(-5.0d0, min(gama(10), 5.0d0))
       gama(10) = 10.0d0**gama(10)
       end if 
    !
       return
    end subroutine mach_hetp_calcact2b
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Calculates multi-component activity coefficients from Bromley's method
    ! ## of ammonium-sulfate-nitrate-sodium-chloride aerosol system for case 
    ! ## I6 and J3
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcact3(h, nh4_t, so4_t, hso4, no3, cl, na, lwn, gama, t)
    !
       use mach_hetp_mod,         only: tiny
       implicit none
    !
       real(kind=8),       intent(in) :: nh4_t
       real(kind=8),       intent(in) :: so4_t
       real(kind=8),       intent(in) :: hso4 
       real(kind=8),       intent(in) :: no3  
       real(kind=8),       intent(in) :: cl  
       real(kind=8),       intent(in) :: na   
       real(kind=8),       intent(in) :: lwn   
       real(kind=8),       intent(in) :: h      
       real(kind=8),       intent(in) :: t       
       real(kind=8),    intent(inout) :: gama(13)
    !
    !  ## Local variables:
       real(kind=8)    :: ionic, sion, c, xx, xx2, c1, c2, hh
       real(kind=8)    :: c3, c4, c5, c5a
       real(kind=8)    :: f11, f12, f13, f14, k1, k2, k3, k4, h2
       real(kind=8)    :: f21, f22, f23, ch1, ch2
       real(kind=8)    :: f2a1, f2a2, f2a3, f2a4, ff11, ff12, ff13
       real(kind=8)    :: g01, g02, g03, g04, g05, g06, g07, g08, g09
       real(kind=8)    :: g10, g11, g12
       logical(kind=4) :: tc
    !
    !  ## Calculate ionic strength of solution
       ionic = h + na + nh4_t + cl + 4.0d0*so4_t + hso4 + no3
       ionic = max(min(0.5d0*ionic/lwn, 100.d0), tiny)
    !
    !  ## Calculate the binary activity coefficients
       sion = sqrt(ionic)
       c3   = exp(-0.023d0*ionic*ionic*ionic)
       c4   = -0.5107d0*sion
       c5   = 1.0d0 + 0.1d0*ionic
       c5a  = 1.0d0 + sion
    !
    !  ## Coefficients at 25C
       if (ionic < 6.0d0) then
          c   = 1.0d0 + 0.12265d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g01 = 0.39495d0 + 0.60505d0*c5**2.23          
          g01 = log10(g01) + xx
    
          c   = 1.0d0 + -0.01045d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g02 = 0.23765d0 + 0.76235d0*c5**(-0.19)       
          g02 = 2.0d0*log10(g02) + 2.0d0*xx
    
          c   = 1.0d0 + -0.02145d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g03 = 0.22465d0 + 0.77535d0*c5**(-0.39)       
          g03 = log10(g03) + xx
    
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + 2.0d0*xx
    
          c   = 1.0d0 + -0.06325d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
    
          c   = 1.0d0 + -0.0055d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + 2.0d0*xx
    
          c   = 1.0d0 + 0.44d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
     
          c   = 1.0d0 + 0.0451d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
    
          c   = 1.0d0 + 0.143d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
    
          c   = 1.0d0 + 0.33d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       else
          xx  = c4/(1.0d0 + sion)
          xx2 = 2.0d0*xx
          g01 = 0.39495d0 + 0.60505d0*c5**2.23          
          g01 = log10(g01) + xx
          g02 = 0.23765d0 + 0.76235d0*c5**(-0.19)       
          g02 = 2.0d0*log10(g02) + xx2
          g03 = 0.22465d0 + 0.77535d0*c5**(-0.39)       
          g03 = log10(g03) + xx
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + xx2
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + xx2
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       end if 
    !
    !  ## Correct coefficients for temperature other than 298 K
       tc = abs(t - 298.0d0) > 1.0d0
       if (tc) then
          c1  = 1.125d0  - 0.005d0*(t - 273.0d0)
          c2  = (0.125d0 - 0.005d0*(t - 273.0d0))*(0.039d0*ionic**0.92d0 - 0.41d0*sion/c5a)
          c3  = 2.0d0*c2
          g01 = c1*g01 - c2         !g01
          g02 = c1*g02 - c3         !g02
          g03 = c1*g03 - c2         !g03
          g04 = c1*g04 - c3         !g04
          g05 = c1*g05 - c2         !g05
          g06 = c1*g06 - c2         !g06
          g07 = c1*g07 - c3         !g07
          g08 = c1*g08 - c2         !g08
          g10 = c1*g10 - c2         !g10
          g11 = c1*g11 - c2         !g11
       end if
    !
       g09 = g06 + g08 - g11  
       g12 = g01 + g08 - g11
    !
    !  ## Calculate multicomponent activity coefficients
       hh = (0.511d0*(298.0d0/t)**1.5)*sion/c5a
    !
       f11  = cl/lwn
       f12  = so4_t/lwn
       f13  = hso4/lwn
       f14  = no3/lwn
       ch1  = 1.0d0/ionic
       ch2  = 2.25d0/ionic
       k1   = ch1*f11
       k2   = ch2*f12
       k3   = ch1*f13
       k4   = ch1*f14
       h2   = 2.0d0*hh
    !
       ff11 = k1*(g11 + hh) + k2*(g07 + h2) + k3*(g08 + hh) + k4*(g10 + hh)
       ff12 = k1*(g01 + hh) + k2*(g02 + h2) + k3*(g12 + hh) + k4*(g03 + hh)
       ff13 = k1*(g06 + hh) + k2*(g04 + h2) + k3*(g09 + hh) + k4*(g05 + hh)
    !
       f21  = h/lwn
       f22  = na/lwn
       f23  = nh4_t/lwn
       k1   = ch1*f21
       k2   = ch1*f22
       k3   = ch1*f23
       f2a1 = k1*(g11 + hh) + k2*(g01 + hh) + k3*(g06 + hh)
       f2a2 = ((ch2*f21)*(g07 + h2) + (ch2*f22)*(g02 + h2) + (ch2*f23)*(g04 + h2))*0.5d0
       f2a3 = k1*(g08 + hh) + k2*(g12 + hh) + k3*(g09 + hh)
       f2a4 = k1*(g10 + hh) + k2*(g03 + hh) + k3*(g05 + hh)
    !
    !  ## log10 of activity coefficients
       gama(1)  = ((ff12 + f2a1) * 0.5d0 - hh)        ! NaCl
       gama(2)  = ((ff12 + f2a2) / 3.0d0 - hh)*2.0d0  ! Na2SO4
       gama(3)  = ((ff12 + f2a4) * 0.5d0 - hh)        ! NaNO3
       gama(4)  = ((ff13 + f2a2) / 3.0d0 - hh)*2.0d0  ! (NH4)2SO4
       gama(5)  = ((ff13 + f2a4) * 0.5d0 - hh)        ! NH4NO3
       gama(6)  = ((ff13 + f2a1) * 0.5d0 - hh)        ! NH4Cl
       gama(7)  = ((ff11 + f2a2) / 3.0d0 - hh)*2.0d0  ! 2H-SO4
       gama(8)  = ((ff11 + f2a3) * 0.5d0 - hh)        ! H-HSO4
       gama(9)  = ((ff13 + f2a3) * 0.5d0 - hh)        ! NH4HSO4
       gama(10) = ((ff11 + f2a4) * 0.5d0 - hh)        ! HNO3
       gama(11) = ((ff11 + f2a1) * 0.5d0 - hh)        ! HCl
       gama(12) = ((ff12 + f2a3) * 0.5d0 - hh)        ! NaHSO4
       gama(13) = 0.6d0*gama(4) + 0.4d0*gama(9)       ! lc; scape
    
    !  ## Convert log(gama) coefficients to gama
       gama(1)  = max(-5.0d0, min(gama(1), 5.0d0))
       gama(1)  = 10.0**gama(1)
       gama(2)  = max(-5.0d0, min(gama(2), 5.0d0))
       gama(2)  = 10.0**gama(2)
       gama(3)  = max(-5.0d0, min(gama(3), 5.0d0))
       gama(3)  = 10.0**gama(3)
       gama(4)  = max(-5.0d0, min(gama(4), 5.0d0))
       gama(4)  = 10.0**gama(4)
       gama(5)  = max(-5.0d0, min(gama(5), 5.0d0))
       gama(5)  = 10.0**gama(5)
       gama(6)  = max(-5.0d0, min(gama(6), 5.0d0))
       gama(6)  = 10.0**gama(6)
       gama(7)  = max(-5.0d0, min(gama(7), 5.0d0))
       gama(7)  = 10.0**gama(7)
       gama(8)  = max(-5.0d0, min(gama(8), 5.0d0))
       gama(8)  = 10.0**gama(8)
       gama(9)  = max(-5.0d0, min(gama(9), 5.0d0))
       gama(9)  = 10.0**gama(9)
       gama(10) = max(-5.0d0, min(gama(10), 5.0d0))
       gama(10) = 10.0**gama(10)
       gama(11) = max(-5.0d0, min(gama(11), 5.0d0))
       gama(11) = 10.0**gama(11)
       gama(12) = max(-5.0d0, min(gama(12), 5.0d0))
       gama(12) = 10.0**gama(12)
       gama(13) = max(-5.0d0, min(gama(13), 5.0d0))
       gama(13) = 10.0**gama(13)
    !
       return
    end subroutine mach_hetp_calcact3
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Calculates multi-component activity coefficients from Bromley's method
    ! ## of ammonium-sulfate-nitrate-sodium-chloride aerosol system for case 
    ! ## G5 and H6
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcact3b(h, nh4_t, so4_t, hso4, no3, cl, na,     &
                                  lwn, gama, t, soln, frst, calain, calou)
    !
       use mach_hetp_mod,         only: tiny
       implicit none
    !
       real(kind=8),       intent(in) :: nh4_t
       real(kind=8),       intent(in) :: so4_t
       real(kind=8),       intent(in) :: hso4 
       real(kind=8),       intent(in) :: no3  
       real(kind=8),       intent(in) :: cl  
       real(kind=8),       intent(in) :: na   
       real(kind=8),       intent(in) :: lwn   
       real(kind=8),       intent(in) :: h      
       real(kind=8),       intent(in) :: t       
       logical(kind=4),    intent(in) :: soln 
       real(kind=8),    intent(inout) :: gama    (13)
       logical(kind=4),    intent(in) :: frst    
       logical(kind=4),    intent(in) :: calain  
       logical(kind=4),    intent(in) :: calou   
    !
    !  ## Local variables:
       real(kind=8)    :: ionic, sion, c, xx, xx2, c1, c2, hh
       real(kind=8)    :: c3, c4, c5, c5a
       real(kind=8)    :: f11, f12, f13, f14, k1, k2, k3, k4, h2
       real(kind=8)    :: f21, f22, f23, ch1, ch2
       real(kind=8)    :: f2a1, f2a2, f2a3, f2a4, ff11, ff12, ff13
       real(kind=8)    :: g01, g02, g03, g04, g05, g06, g07, g08, g09
       real(kind=8)    :: g10, g11, g12
       logical(kind=4) :: tc
    !
    !
       if ((( .not. soln)) .and.                          &
          ((frst  .and. calou) .or.    &
           (( .not. frst) .and. calain)) ) then
    !  ## Calculate ionic strength of solution
       ionic = h + na + nh4_t + cl + 4.0d0*so4_t + hso4 + no3
       ionic = max(min(0.5d0*ionic/lwn, 100.d0), tiny)
    !
    !  ## Calculate the binary activity coefficients
       sion = sqrt(ionic)
       c3   = exp(-0.023d0*ionic*ionic*ionic)
       c4   = -0.5107d0*sion
       c5   = 1.0d0 + 0.1d0*ionic
       c5a  = 1.0d0 + sion
    !
    !  ## Coefficients at 25C
       if (ionic < 6.0d0) then
          c   = 1.0d0 + 0.12265d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g01 = 0.39495d0 + 0.60505d0*c5**2.23          
          g01 = log10(g01) + xx
    
          c   = 1.0d0 + -0.01045d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g02 = 0.23765d0 + 0.76235d0*c5**(-0.19)       
          g02 = 2.0d0*log10(g02) + 2.0d0*xx
    
          c   = 1.0d0 + -0.02145d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g03 = 0.22465d0 + 0.77535d0*c5**(-0.39)       
          g03 = log10(g03) + xx
    
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + 2.0d0*xx
    
          c   = 1.0d0 + -0.06325d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
    
          c   = 1.0d0 + -0.0055d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + 2.0d0*xx
    
          c   = 1.0d0 + 0.44d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
     
          c   = 1.0d0 + 0.0451d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
    
          c   = 1.0d0 + 0.143d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
    
          c   = 1.0d0 + 0.33d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       else
          xx  = c4/(1.0d0 + sion)
          xx2 = 2.0d0*xx
          g01 = 0.39495d0 + 0.60505d0*c5**2.23          
          g01 = log10(g01) + xx
          g02 = 0.23765d0 + 0.76235d0*c5**(-0.19)       
          g02 = 2.0d0*log10(g02) + xx2
          g03 = 0.22465d0 + 0.77535d0*c5**(-0.39)       
          g03 = log10(g03) + xx
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + xx2
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + xx2
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
       end if 
    !
    !  ## Correct coefficients for temperature other than 298 K
       tc = abs(t - 298.0d0) > 1.0d0
       if (tc) then
          c1  = 1.125d0  - 0.005d0*(t - 273.0d0)
          c2  = (0.125d0 - 0.005d0*(t - 273.0d0))*(0.039d0*ionic**0.92d0 - 0.41d0*sion/c5a)
          c3  = 2.0d0*c2
          g01 = c1*g01 - c2         !g01
          g02 = c1*g02 - c3         !g02
          g03 = c1*g03 - c2         !g03
          g04 = c1*g04 - c3         !g04
          g05 = c1*g05 - c2         !g05
          g06 = c1*g06 - c2         !g06
          g07 = c1*g07 - c3         !g07
          g08 = c1*g08 - c2         !g08
          g10 = c1*g10 - c2         !g10
          g11 = c1*g11 - c2         !g11
       end if
    !
       g09 = g06 + g08 - g11  
       g12 = g01 + g08 - g11
    !
    !  ## Calculate multicomponent activity coefficients
       hh = (0.511d0*(298.0d0/t)**1.5)*sion/c5a
    !
       f11  = cl/lwn
       f12  = so4_t/lwn
       f13  = hso4/lwn
       f14  = no3/lwn
       ch1  = 1.0d0/ionic
       ch2  = 2.25d0/ionic
       k1   = ch1*f11
       k2   = ch2*f12
       k3   = ch1*f13
       k4   = ch1*f14
       h2   = 2.0d0*hh
    !
       ff11 = k1*(g11 + hh) + k2*(g07 + h2) + k3*(g08 + hh) + k4*(g10 + hh)
       ff12 = k1*(g01 + hh) + k2*(g02 + h2) + k3*(g12 + hh) + k4*(g03 + hh)
       ff13 = k1*(g06 + hh) + k2*(g04 + h2) + k3*(g09 + hh) + k4*(g05 + hh)
    !
       f21  = h/lwn
       f22  = na/lwn
       f23  = nh4_t/lwn
       k1   = ch1*f21
       k2   = ch1*f22
       k3   = ch1*f23
       f2a1 = k1*(g11 + hh) + k2*(g01 + hh) + k3*(g06 + hh)
       f2a2 = ((ch2*f21)*(g07 + h2) + (ch2*f22)*(g02 + h2) + (ch2*f23)*(g04 + h2))*0.5d0
       f2a3 = k1*(g08 + hh) + k2*(g12 + hh) + k3*(g09 + hh)
       f2a4 = k1*(g10 + hh) + k2*(g03 + hh) + k3*(g05 + hh)
    !
    !  ## log10 of activity coefficients
       gama(1)  = ((ff12 + f2a1) * 0.5d0 - hh)        ! NaCl
       gama(2)  = ((ff12 + f2a2) / 3.0d0 - hh)*2.0d0  ! Na2SO4
       gama(3)  = ((ff12 + f2a4) * 0.5d0 - hh)        ! NaNO3
       gama(4)  = ((ff13 + f2a2) / 3.0d0 - hh)*2.0d0  ! (NH4)2SO4
       gama(5)  = ((ff13 + f2a4) * 0.5d0 - hh)        ! NH4NO3
       gama(6)  = ((ff13 + f2a1) * 0.5d0 - hh)        ! NH4Cl
       gama(7)  = ((ff11 + f2a2) / 3.0d0 - hh)*2.0d0  ! 2H-SO4
       gama(8)  = ((ff11 + f2a3) * 0.5d0 - hh)        ! H-HSO4
       gama(9)  = ((ff13 + f2a3) * 0.5d0 - hh)        ! NH4HSO4
       gama(10) = ((ff11 + f2a4) * 0.5d0 - hh)        ! HNO3
       gama(11) = ((ff11 + f2a1) * 0.5d0 - hh)        ! HCl
       gama(12) = ((ff12 + f2a3) * 0.5d0 - hh)        ! NaHSO4
       gama(13) = 0.6d0*gama(4) + 0.4d0*gama(9)       ! lc; scape
    
    !  ## Convert log(gama) coefficients to gama
       gama(1)  = max(-5.0d0, min(gama(1), 5.0d0))
       gama(1)  = 10.0**gama(1)
       gama(2)  = max(-5.0d0, min(gama(2), 5.0d0))
       gama(2)  = 10.0**gama(2)
       gama(3)  = max(-5.0d0, min(gama(3), 5.0d0))
       gama(3)  = 10.0**gama(3)
       gama(4)  = max(-5.0d0, min(gama(4), 5.0d0))
       gama(4)  = 10.0**gama(4)
       gama(5)  = max(-5.0d0, min(gama(5), 5.0d0))
       gama(5)  = 10.0**gama(5)
       gama(6)  = max(-5.0d0, min(gama(6), 5.0d0))
       gama(6)  = 10.0**gama(6)
       gama(7)  = max(-5.0d0, min(gama(7), 5.0d0))
       gama(7)  = 10.0**gama(7)
       gama(8)  = max(-5.0d0, min(gama(8), 5.0d0))
       gama(8)  = 10.0**gama(8)
       gama(9)  = max(-5.0d0, min(gama(9), 5.0d0))
       gama(9)  = 10.0**gama(9)
       gama(10) = max(-5.0d0, min(gama(10), 5.0d0))
       gama(10) = 10.0**gama(10)
       gama(11) = max(-5.0d0, min(gama(11), 5.0d0))
       gama(11) = 10.0**gama(11)
       gama(12) = max(-5.0d0, min(gama(12), 5.0d0))
       gama(12) = 10.0**gama(12)
       gama(13) = max(-5.0d0, min(gama(13), 5.0d0))
       gama(13) = 10.0**gama(13)
       end if 
    !
       return
    end subroutine mach_hetp_calcact3b
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Calculates multi-component activity coefficients from Bromley's method
    ! ## of ammonium-sulfate-nitrate-sodium-chloride-calcium-potassium-magnesium
    ! ## aerosol system for case K4 and L9
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcact4(h, nh4_t, so4_t, hso4, no3, cl, na,     &
                                  ca, pk, mg, lwn, gama, t)
    !
       use mach_hetp_mod,         only: tiny
       implicit none
    !
       real(kind=8),    intent   (in) :: nh4_t 
       real(kind=8),    intent   (in) :: so4_t 
       real(kind=8),    intent   (in) :: hso4   
       real(kind=8),    intent   (in) :: no3    
       real(kind=8),    intent   (in) :: cl     
       real(kind=8),    intent   (in) :: na    
       real(kind=8),    intent   (in) :: ca    
       real(kind=8),    intent   (in) :: pk      
       real(kind=8),    intent   (in) :: mg      
       real(kind=8),    intent   (in) :: lwn     
       real(kind=8),    intent   (in) :: h      
       real(kind=8),    intent   (in) :: t       
       real(kind=8),    intent(inout) :: gama    (23)
    !
    !  ## Local variables:
    !
       real(kind=8)    :: ionic, sion, tc, tx, c, xx, xx2, c1, c2, hh
       real(kind=8)    :: c3, c4, c5, c5a
       real(kind=8)    :: f11, f12, f13, f14, k1, k2, k3, k4, h2, h3
       real(kind=8)    :: f21, f22, f23, ch1, ch2, ch3
       real(kind=8)    :: f2a1, f2a2, f2a3, f2a4, f2b1, f2b2, f2b3, f2b4
       real(kind=8)    :: ff11, ff12, ff13, ff14, ff15, ff16
       real(kind=8)    :: g01, g02, g03, g04, g05, g06, g07, g08, g09
       real(kind=8)    :: g10, g11, g12, g15, g16, g17, g18, g19, g20
       real(kind=8)    :: g21, g22, g23
    !
    !  ## Calculate ionic strength of solution
       ionic = h + na + nh4_t + cl + 4.0d0*so4_t + hso4  &
               + no3 + 4.0d0*ca + pk + 4.0d0*mg
       ionic = max(min(0.5d0*ionic/lwn, 100.d0), tiny)
    
    !  ## Calculate the binary activity coefficients
       sion = sqrt(ionic)
       c3   = exp(-0.023d0*ionic*ionic*ionic)
       c4   = -0.5107d0*sion
       c5   = 1.0d0 + 0.1d0*ionic
       c5a  = 1.0d0 + sion
    !
    !  ## Coefficients at 25C
       if (ionic < 6.0d0) then
          c   = 1.0d0 + 0.12265d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g01 = 0.39495d0 + 0.60505d0*c5**2.23          
          g01 = log10(g01) + xx
    
          c   = 1.0d0 + -0.01045d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g02 = 0.23765d0 + 0.76235d0*c5**(-0.19)       
          g02 = 2.0d0*log10(g02) + 2.0d0*xx
    
          c   = 1.0d0 + -0.02145d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g03 = 0.22465d0 + 0.77535d0*c5**(-0.39)       
          g03 = log10(g03) + xx
    
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + 2.0d0*xx
    
          c   = 1.0d0 + -0.06325d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
    
          c   = 1.0d0 + -0.0055d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + 2.0d0*xx
    
          c   = 1.0d0 + 0.44d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
     
          c   = 1.0d0 + 0.0451d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
    
          c   = 1.0d0 + 0.143d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
    
          c   = 1.0d0 + 0.33d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
    
          c   = 1.0d0 + 0.05115d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g15 = 0.31045d0 + 0.68955d0*c5**0.93  
          g15 = 2.0d0*log10(g15) + 2.0d0*xx
    
          c   = 1.0d0 + 0.132d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g16 = 0.406d0 + 0.594d0*c5**2.40               
          g16 = 2.0d0*log10(g16) + 2.0d0*xx
    
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g17 = 0.23375d0 + 0.76625d0*c5**(-0.25) 
          g17 = 2.0d0*log10(g17) + 2.0d0*xx
    
          c   = 1.0d0 + -0.12815d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g19 = 0.09855d0 + 0.90145d0*c5**(-2.33)        
          g19 = log10(g19) + xx
    
          c   = 1.0d0 + 0.0506d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g20 = 0.3098d0 + 0.6902d0*c5**0.92             
          g20 = log10(g20) + xx
    
          c   = 1.0d0 + 0.00825d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g21 = 0.25975d0 + 0.74025d0*c5**0.15            
          g21 = 4.0d0*log10(g21) + 4.0d0*xx
    
          c   = 1.0d0 + 0.1276d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g22 = 0.4008d0 + 0.5992d0*c5**2.32           
          g22 = 2.0d0*log10(g22) + 2.0d0*xx
    
          c   = 1.0d0 + 0.1595d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g23 = 0.4385d0 + 0.5615d0*c5**2.90           
          g23 = 2.0d0*log10(g23) + 2.0d0*xx
       else
          xx  = c4/(1.0d0 + sion)
          xx2 = 2.0d0*xx
          g01 = 0.39495d0 + 0.60505d0*c5**2.23          
          g01 = log10(g01) + xx
          g02 = 0.23765d0 + 0.76235d0*c5**(-0.19)       
          g02 = 2.0d0*log10(g02) + xx2
          g03 = 0.22465d0 + 0.77535d0*c5**(-0.39)       
          g03 = log10(g03) + xx
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + xx2
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + xx2
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
          g15 = 0.31045d0 + 0.68955d0*c5**0.93  
          g15 = 2.0d0*log10(g15) + xx2
          g16 = 0.406d0 + 0.594d0*c5**2.40               
          g16 = 2.0d0*log10(g16) + xx2
          g17 = 0.23375d0 + 0.76625d0*c5**(-0.25) 
          g17 = 2.0d0*log10(g17) + xx2
          g19 = 0.09855d0 + 0.90145d0*c5**(-2.33)        
          g19 = log10(g19) + xx
          g20 = 0.3098d0 + 0.6902d0*c5**0.92             
          g20 = log10(g20) + xx
          g21 = 0.25975d0 + 0.74025d0*c5**0.15            
          g21 = 4.0d0*log10(g21) + 4.0d0*xx
          g22 = 0.4008d0 + 0.5992d0*c5**2.32            
          g22 = 2.0d0*log10(g22) + xx2
          g23 = 0.4385d0 + 0.5615d0*c5**2.90              
          g23 = 2.0d0*log10(g23) + xx2
       end if 
    !
    !  ## Correct for temperature other than 298 K
       tc = abs(t - 298.0d0)
       tx = 0.005d0*(t - 273.0d0)
       if (tc > 1.0d0) then
          c1  = 1.125d0  - tx
          c2  = (0.125d0 - tx)*(0.039d0*ionic**0.92 - 0.41d0*sion/c5a)
          c3  = c2*2.0d0
          g01 = c1*g01 - c2         !g01
          g02 = c1*g02 - c3         !g02
          g03 = c1*g03 - c2         !g03
          g04 = c1*g04 - c3         !g04
          g05 = c1*g05 - c2         !g05
          g06 = c1*g06 - c2         !g06
          g07 = c1*g07 - c3         !g07
          g08 = c1*g08 - c2         !g08
          g10 = c1*g10 - c2         !g10
          g11 = c1*g11 - c2         !g11
          g15 = c1*g15 - c3         !g15
          g16 = c1*g16 - c3         !g16
          g17 = c1*g17 - c3         !g17
          g19 = c1*g19 - c2         !g19
          g20 = c1*g20 - c2         !g20
          g21 = c1*g21 - c2*4.0d0   !g21
          g22 = c1*g22 - c3         !g22
          g23 = c1*g23 - c3         !g23
       end if
    !
       g09 = g06 + g08 - g11  
       g12 = g01 + g08 - g11   
       g18 = g08 + g20 - g11   
    !
    !  ## Calculate multicomponent activity coefficients
       hh   = (0.511d0*(298.0d0/t)**1.5)*sion/c5a
       f11  = cl/lwn
       f12  = so4_t/lwn
       f13  = hso4/lwn
       f14  = no3/lwn
       ch1  = 1.0d0/ionic
       ch2  = 2.25d0/ionic
       k1   = ch1*f11
       k2   = ch2*f12
       k3   = ch1*f13
       k4   = ch1*f14
       h2   = 2.0d0*hh
    !
       ff11 = k1*(g11 + hh) + k2*(g07 + h2) + k3*(g08 + hh) + k4*(g10 + hh)
       ff12 = k1*(g01 + hh) + k2*(g02 + h2) + k3*(g12 + hh) + k4*(g03 + hh)
       ff13 = k1*(g06 + hh) + k2*(g04 + h2) + k3*(g09 + hh) + k4*(g05 + hh)
    !
       f21   = h/lwn
       f22   = na/lwn
       f23   = nh4_t/lwn
       k1    = ch1*f21
       k2    = ch1*f22
       k3    = ch1*f23
       f2a1  = k1*(g11 + hh) + k2*(g01 + hh) + k3*(g06 + hh)
       f2a2  = ((ch2*f21)*(g07 + h2) + (ch2*f22)*(g02 + h2) + (ch2*f23)*(g04 + h2))*0.5d0
       f2a3  = k1*(g08 + hh) + k2*(g12 + hh) + k3*(g09 + hh)
       f2a4  = k1*(g10 + hh) + k2*(g03 + hh) + k3*(g05 + hh)
    !
       ch1 = 2.25d0/ionic
       ch2 = 4.0d0/ionic
       ch3 = 1.0d0/ionic
       h3  = 4.0d0*hh
       f11 = cl/lwn
       f12 = so4_t/lwn
       f13 = no3/lwn
       f14 = hso4/lwn
       k1  = ch1*f11
       k2  = ch2*f12
       k3  = ch1*f13
    !
       ff14 = (k1*(g16 + h2) + k2*h3 + k3*(g15 + h2))*0.5d0
       ff15 = (ch3*f11)*(g20 + hh) + (ch1*f12)*(g17 + h2) + (ch3*f14)*(g18 + hh) +  &
                 (ch3*f13)*(g19 + hh)
       ff16 = (k1*(g23 + h2) + k2*(g21 + h3) + k3*(g22 + h2))*0.5d0
    !
       f21  = ca/lwn
       f22  = pk/lwn
       f23  = mg/lwn
       k1   = ch1*f21
       k2   = ch1*f23
       k3   = ch3*f22
       f2b1 = k1*(g16 + h2) + k3*(g20 + hh) + k2*(g23 + h2)
       f2b2 = ((ch2*f21)*h3 + (ch1*f22)*(g17 + h2) + (ch2*f23)*(g21 + h3))*0.5d0
       f2b3 = (ch3*f22)*(g18 + hh)
       f2b4 = k1*(g15 + h2) + k3*(g19 + hh) + k2*(g22 + h2)
    !
    !  ## log10 of activity coefficients
       gama(1)   = ((ff12 + f2a1) * 0.5d0 - hh)        ! NaCl
       gama(2)   = ((ff12 + f2a2) / 3.0d0 - hh)*2.0d0  ! Na2SO4
       gama(3)   = ((ff12 + f2a4) * 0.5d0 - hh)        ! NaNO3
       gama(4)   = ((ff13 + f2a2) / 3.0d0 - hh)*2.0d0  ! (NH4)2SO4
       gama(5)   = ((ff13 + f2a4) * 0.5d0 - hh)        ! NH4NO3
       gama(6)   = ((ff13 + f2a1) * 0.5d0 - hh)        ! NH4Cl
       gama(7)   = ((ff11 + f2a2) / 3.0d0 - hh)*2.0d0  ! 2H-SO4
       gama(8)   = ((ff11 + f2a3) * 0.5d0 - hh)        ! H-HSO4
       gama(9)   = ((ff13 + f2a3) * 0.5d0 - hh)        ! NH4HSO4
       gama(10)  = ((ff11 + f2a4) * 0.5d0 - hh)        ! HNO3
       gama(11)  = ((ff11 + f2a1) * 0.5d0 - hh)        ! HCl
       gama(12)  = ((ff12 + f2a3) * 0.5d0 - hh)        ! NaHSO4
       gama(13)  = 0.6d0*gama(4) + 0.4d0*gama(9)       ! lc; scape
       gama(14)  = 0.0d0                               ! CaSO4
       gama(15)  = ((ff14 + f2b4) / 3.0d0 - hh)*2.0d0  ! Ca(NO3)2
       gama(16)  = ((ff14 + f2b1) / 3.0d0 - hh)*2.0d0  ! CaCl2
       gama(17)  = ((ff15 + f2b2) / 3.0d0 - hh)*2.0d0  ! K2SO4
       gama(18)  = ((ff15 + f2b3) * 0.5d0 - hh)        ! KHSO4
       gama(19)  = ((ff15 + f2b4) * 0.5d0 - hh)        ! KNO3
       gama(20)  = ((ff15 + f2b1) * 0.5d0 - hh)        ! KCl
       gama(21)  = ((ff16 + f2b2) * 0.25d0- hh)*4.0d0  ! MgSO4
       gama(22)  = ((ff16 + f2b4) / 3.0d0 - hh)*2.0d0  ! Mg(NO3)2
       gama(23)  = ((ff16 + f2b1) / 3.0d0 - hh)*2.0d0  ! MgCl2
    
    !  ## Convert log(gama) coefficients to gama
       gama(1)  = max(-5.0d0, min(gama(1), 5.0d0))
       gama(1)  = 10.0**gama(1)
       gama(2)  = max(-5.0d0, min(gama(2), 5.0d0))
       gama(2)  = 10.0**gama(2)
       gama(3)  = max(-5.0d0, min(gama(3), 5.0d0))
       gama(3)  = 10.0**gama(3)
       gama(4)  = max(-5.0d0, min(gama(4), 5.0d0))
       gama(4)  = 10.0**gama(4)
       gama(5)  = max(-5.0d0, min(gama(5), 5.0d0))
       gama(5)  = 10.0**gama(5)
       gama(6)  = max(-5.0d0, min(gama(6), 5.0d0))
       gama(6)  = 10.0**gama(6)
       gama(7)  = max(-5.0d0, min(gama(7), 5.0d0))
       gama(7)  = 10.0**gama(7)
       gama(8)  = max(-5.0d0, min(gama(8), 5.0d0))
       gama(8)  = 10.0**gama(8)
       gama(9)  = max(-5.0d0, min(gama(9), 5.0d0))
       gama(9)  = 10.0**gama(9)
       gama(10) = max(-5.0d0, min(gama(10), 5.0d0))
       gama(10) = 10.0**gama(10)
       gama(11) = max(-5.0d0, min(gama(11), 5.0d0))
       gama(11) = 10.0**gama(11)
       gama(12) = max(-5.0d0, min(gama(12), 5.0d0))
       gama(12) = 10.0**gama(12)
       gama(13) = max(-5.0d0, min(gama(13), 5.0d0))
       gama(13) = 10.0**gama(13)
       gama(14) = 1.0d0
       gama(15) = max(-5.0d0, min(gama(15), 5.0d0))
       gama(15) = 10.0**gama(15)
       gama(16) = max(-5.0d0, min(gama(16), 5.0d0))
       gama(16) = 10.0**gama(16)
       gama(17) = max(-5.0d0, min(gama(17), 5.0d0))
       gama(17) = 10.0**gama(17)
       gama(18) = max(-5.0d0, min(gama(18), 5.0d0))
       gama(18) = 10.0**gama(18)
       gama(19) = max(-5.0d0, min(gama(19), 5.0d0))
       gama(19) = 10.0**gama(19)
       gama(20) = max(-5.0d0, min(gama(20), 5.0d0))
       gama(20) = 10.0**gama(20)
       gama(21) = max(-5.0d0, min(gama(21), 5.0d0))
       gama(21) = 10.0**gama(21)
       gama(22) = max(-5.0d0, min(gama(22), 5.0d0))
       gama(22) = 10.0**gama(22)
       gama(23) = max(-5.0d0, min(gama(23), 5.0d0))
       gama(23) = 10.0**gama(23)
    !
       return
    end subroutine mach_hetp_calcact4
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Calculates multi-component activity coefficients from Bromley's method
    ! ## of ammonium-sulfate-nitrate-sodium-chloride-calcium-potassium-magnesium
    ! ## aerosol system for case O7, M8 and P13
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_calcact4b(h, nh4_t, so4_t, hso4, no3, cl, na,     &
                                   ca, pk, mg, lwn, gama, t, soln, frst,   &
                                   calain, calou)
    !
       use mach_hetp_mod,         only: tiny
       implicit none
    !
       real(kind=8),    intent   (in) :: nh4_t 
       real(kind=8),    intent   (in) :: so4_t 
       real(kind=8),    intent   (in) :: hso4   
       real(kind=8),    intent   (in) :: no3    
       real(kind=8),    intent   (in) :: cl     
       real(kind=8),    intent   (in) :: na    
       real(kind=8),    intent   (in) :: ca    
       real(kind=8),    intent   (in) :: pk      
       real(kind=8),    intent   (in) :: mg      
       real(kind=8),    intent   (in) :: lwn     
       real(kind=8),    intent   (in) :: h      
       real(kind=8),    intent   (in) :: t       
       real(kind=8),    intent(inout) :: gama    (23)
       logical(kind=4), intent   (in) :: soln 
       logical(kind=4), intent   (in) :: frst    
       logical(kind=4), intent   (in) :: calain  
       logical(kind=4), intent   (in) :: calou   
    !
    !  ## Local variables:
    !
       real(kind=8)    :: ionic, sion, tc, tx, c, xx, xx2, c1, c2, hh
       real(kind=8)    :: c3, c4, c5, c5a
       real(kind=8)    :: f11, f12, f13, f14, k1, k2, k3, k4, h2, h3
       real(kind=8)    :: f21, f22, f23, ch1, ch2, ch3
       real(kind=8)    :: f2a1, f2a2, f2a3, f2a4, f2b1, f2b2, f2b3, f2b4
       real(kind=8)    :: ff11, ff12, ff13, ff14, ff15, ff16
       real(kind=8)    :: g01, g02, g03, g04, g05, g06, g07, g08, g09
       real(kind=8)    :: g10, g11, g12, g15, g16, g17, g18, g19, g20
       real(kind=8)    :: g21, g22, g23
    !
    !
       if ((( .not. soln)) .and.                          &
          ((frst  .and. calou ) .or.    &
           (( .not. frst) .and. calain)) ) then
    !
    !  ## Calculate ionic strength of solution
       ionic = h + na + nh4_t + cl + 4.0d0*so4_t + hso4  &
               + no3 + 4.0d0*ca + pk + 4.0d0*mg
       ionic = max(min(0.5d0*ionic/lwn, 100.d0), tiny)
    
    !  ## Calculate the binary activity coefficients
       sion = sqrt(ionic)
       c3   = exp(-0.023d0*ionic*ionic*ionic)
       c4   = -0.5107d0*sion
       c5   = 1.0d0 + 0.1d0*ionic
       c5a  = 1.0d0 + sion
    !
    !  ## Coefficients at 25C
       if (ionic < 6.0d0) then
          c   = 1.0d0 + 0.12265d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g01 = 0.39495d0 + 0.60505d0*c5**2.23          
          g01 = log10(g01) + xx
    
          c   = 1.0d0 + -0.01045d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g02 = 0.23765d0 + 0.76235d0*c5**(-0.19)       
          g02 = 2.0d0*log10(g02) + 2.0d0*xx
    
          c   = 1.0d0 + -0.02145d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g03 = 0.22465d0 + 0.77535d0*c5**(-0.39)       
          g03 = log10(g03) + xx
    
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + 2.0d0*xx
    
          c   = 1.0d0 + -0.06325d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
    
          c   = 1.0d0 + -0.0055d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + 2.0d0*xx
    
          c   = 1.0d0 + 0.44d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
     
          c   = 1.0d0 + 0.0451d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
    
          c   = 1.0d0 + 0.143d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
    
          c   = 1.0d0 + 0.33d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
    
          c   = 1.0d0 + 0.05115d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g15 = 0.31045d0 + 0.68955d0*c5**0.93  
          g15 = 2.0d0*log10(g15) + 2.0d0*xx
    
          c   = 1.0d0 + 0.132d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g16 = 0.406d0 + 0.594d0*c5**2.40               
          g16 = 2.0d0*log10(g16) + 2.0d0*xx
    
          c   = 1.0d0 + -0.01375d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g17 = 0.23375d0 + 0.76625d0*c5**(-0.25) 
          g17 = 2.0d0*log10(g17) + 2.0d0*xx
    
          c   = 1.0d0 + -0.12815d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g19 = 0.09855d0 + 0.90145d0*c5**(-2.33)        
          g19 = log10(g19) + xx
    
          c   = 1.0d0 + 0.0506d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g20 = 0.3098d0 + 0.6902d0*c5**0.92             
          g20 = log10(g20) + xx
    
          c   = 1.0d0 + 0.00825d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g21 = 0.25975d0 + 0.74025d0*c5**0.15            
          g21 = 4.0d0*log10(g21) + 4.0d0*xx
    
          c   = 1.0d0 + 0.1276d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g22 = 0.4008d0 + 0.5992d0*c5**2.32           
          g22 = 2.0d0*log10(g22) + 2.0d0*xx
    
          c   = 1.0d0 + 0.1595d0*c3
          xx  = c4/(1.0d0 + c*sion)
          g23 = 0.4385d0 + 0.5615d0*c5**2.90           
          g23 = 2.0d0*log10(g23) + 2.0d0*xx
       else
          xx  = c4/(1.0d0 + sion)
          xx2 = 2.0d0*xx
          g01 = 0.39495d0 + 0.60505d0*c5**2.23          
          g01 = log10(g01) + xx
          g02 = 0.23765d0 + 0.76235d0*c5**(-0.19)       
          g02 = 2.0d0*log10(g02) + xx2
          g03 = 0.22465d0 + 0.77535d0*c5**(-0.39)       
          g03 = log10(g03) + xx
          g04 = 0.23375d0 + 0.76625d0*c5**(-0.25)  
          g04 = 2.0d0*log10(g04) + xx2
          g05 = 0.17525d0 + 0.82475d0*c5**(-1.15)         
          g05 = log10(g05) + xx
          g07 = 0.2435d0 + 0.7565d0*c5**(-0.1)           
          g07 = 2.0d0*log10(g07) + xx2
          g08 = 0.77d0 + 0.23d0*c5*c5*c5*c5*c5*c5*c5*c5
          g08 = log10(g08) + xx
          g06 = 0.3033d0 + 0.6967d0*c5**0.82           
          g06 = log10(g06) + xx
          g10 = 0.419d0 + 0.581d0*c5**2.60           
          g10 = log10(g10) + xx
          g11 = 0.64d0 + 0.36d0*c5*c5*c5*c5*c5*c5
          g11 = log10(g11) + xx
          g15 = 0.31045d0 + 0.68955d0*c5**0.93  
          g15 = 2.0d0*log10(g15) + xx2
          g16 = 0.406d0 + 0.594d0*c5**2.40               
          g16 = 2.0d0*log10(g16) + xx2
          g17 = 0.23375d0 + 0.76625d0*c5**(-0.25) 
          g17 = 2.0d0*log10(g17) + xx2
          g19 = 0.09855d0 + 0.90145d0*c5**(-2.33)        
          g19 = log10(g19) + xx
          g20 = 0.3098d0 + 0.6902d0*c5**0.92             
          g20 = log10(g20) + xx
          g21 = 0.25975d0 + 0.74025d0*c5**0.15            
          g21 = 4.0d0*log10(g21) + 4.0d0*xx
          g22 = 0.4008d0 + 0.5992d0*c5**2.32            
          g22 = 2.0d0*log10(g22) + xx2
          g23 = 0.4385d0 + 0.5615d0*c5**2.90              
          g23 = 2.0d0*log10(g23) + xx2
       end if 
    !
    !  ## Correct for temperature other than 298 K
       tc = abs(t - 298.0d0)
       tx = 0.005d0*(t - 273.0d0)
       if (tc > 1.0d0) then
          c1  = 1.125d0  - tx
          c2  = (0.125d0 - tx)*(0.039d0*ionic**0.92 - 0.41d0*sion/c5a)
          c3  = c2*2.0d0
          g01 = c1*g01 - c2         !g01
          g02 = c1*g02 - c3         !g02
          g03 = c1*g03 - c2         !g03
          g04 = c1*g04 - c3         !g04
          g05 = c1*g05 - c2         !g05
          g06 = c1*g06 - c2         !g06
          g07 = c1*g07 - c3         !g07
          g08 = c1*g08 - c2         !g08
          g10 = c1*g10 - c2         !g10
          g11 = c1*g11 - c2         !g11
          g15 = c1*g15 - c3         !g15
          g16 = c1*g16 - c3         !g16
          g17 = c1*g17 - c3         !g17
          g19 = c1*g19 - c2         !g19
          g20 = c1*g20 - c2         !g20
          g21 = c1*g21 - c2*4.0d0   !g21
          g22 = c1*g22 - c3         !g22
          g23 = c1*g23 - c3         !g23
       end if
    !
       g09 = g06 + g08 - g11  
       g12 = g01 + g08 - g11   
       g18 = g08 + g20 - g11   
    !
    !  ## Calculate multicomponent activity coefficients
       hh   = (0.511d0*(298.0d0/t)**1.5)*sion/c5a
       f11  = cl/lwn
       f12  = so4_t/lwn
       f13  = hso4/lwn
       f14  = no3/lwn
       ch1  = 1.0d0/ionic
       ch2  = 2.25d0/ionic
       k1   = ch1*f11
       k2   = ch2*f12
       k3   = ch1*f13
       k4   = ch1*f14
       h2   = 2.0d0*hh
    !
       ff11 = k1*(g11 + hh) + k2*(g07 + h2) + k3*(g08 + hh) + k4*(g10 + hh)
       ff12 = k1*(g01 + hh) + k2*(g02 + h2) + k3*(g12 + hh) + k4*(g03 + hh)
       ff13 = k1*(g06 + hh) + k2*(g04 + h2) + k3*(g09 + hh) + k4*(g05 + hh)
    !
       f21   = h/lwn
       f22   = na/lwn
       f23   = nh4_t/lwn
       k1    = ch1*f21
       k2    = ch1*f22
       k3    = ch1*f23
       f2a1  = k1*(g11 + hh) + k2*(g01 + hh) + k3*(g06 + hh)
       f2a2  = ((ch2*f21)*(g07 + h2) + (ch2*f22)*(g02 + h2) + (ch2*f23)*(g04 + h2))*0.5d0
       f2a3  = k1*(g08 + hh) + k2*(g12 + hh) + k3*(g09 + hh)
       f2a4  = k1*(g10 + hh) + k2*(g03 + hh) + k3*(g05 + hh)
    !
       ch1 = 2.25d0/ionic
       ch2 = 4.0d0/ionic
       ch3 = 1.0d0/ionic
       h3  = 4.0d0*hh
       f11 = cl/lwn
       f12 = so4_t/lwn
       f13 = no3/lwn
       f14 = hso4/lwn
       k1  = ch1*f11
       k2  = ch2*f12
       k3  = ch1*f13
    !
       ff14 = (k1*(g16 + h2) + k2*h3 + k3*(g15 + h2))*0.5d0
       ff15 = (ch3*f11)*(g20 + hh) + (ch1*f12)*(g17 + h2) + (ch3*f14)*(g18 + hh) +  &
                 (ch3*f13)*(g19 + hh)
       ff16 = (k1*(g23 + h2) + k2*(g21 + h3) + k3*(g22 + h2))*0.5d0
    !
       f21  = ca/lwn
       f22  = pk/lwn
       f23  = mg/lwn
       k1   = ch1*f21
       k2   = ch1*f23
       k3   = ch3*f22
       f2b1 = k1*(g16 + h2) + k3*(g20 + hh) + k2*(g23 + h2)
       f2b2 = ((ch2*f21)*h3 + (ch1*f22)*(g17 + h2) + (ch2*f23)*(g21 + h3))*0.5d0
       f2b3 = (ch3*f22)*(g18 + hh)
       f2b4 = k1*(g15 + h2) + k3*(g19 + hh) + k2*(g22 + h2)
    !
    !  ## log10 of activity coefficients
       gama(1)   = ((ff12 + f2a1) * 0.5d0 - hh)        ! NaCl
       gama(2)   = ((ff12 + f2a2) / 3.0d0 - hh)*2.0d0  ! Na2SO4
       gama(3)   = ((ff12 + f2a4) * 0.5d0 - hh)        ! NaNO3
       gama(4)   = ((ff13 + f2a2) / 3.0d0 - hh)*2.0d0  ! (NH4)2SO4
       gama(5)   = ((ff13 + f2a4) * 0.5d0 - hh)        ! NH4NO3
       gama(6)   = ((ff13 + f2a1) * 0.5d0 - hh)        ! NH4Cl
       gama(7)   = ((ff11 + f2a2) / 3.0d0 - hh)*2.0d0  ! 2H-SO4
       gama(8)   = ((ff11 + f2a3) * 0.5d0 - hh)        ! H-HSO4
       gama(9)   = ((ff13 + f2a3) * 0.5d0 - hh)        ! NH4HSO4
       gama(10)  = ((ff11 + f2a4) * 0.5d0 - hh)        ! HNO3
       gama(11)  = ((ff11 + f2a1) * 0.5d0 - hh)        ! HCl
       gama(12)  = ((ff12 + f2a3) * 0.5d0 - hh)        ! NaHSO4
       gama(13)  = 0.6d0*gama(4) + 0.4d0*gama(9)       ! lc; scape
       gama(14)  = 0.0d0                               ! CaSO4
       gama(15)  = ((ff14 + f2b4) / 3.0d0 - hh)*2.0d0  ! Ca(NO3)2
       gama(16)  = ((ff14 + f2b1) / 3.0d0 - hh)*2.0d0  ! CaCl2
       gama(17)  = ((ff15 + f2b2) / 3.0d0 - hh)*2.0d0  ! K2SO4
       gama(18)  = ((ff15 + f2b3) * 0.5d0 - hh)        ! KHSO4
       gama(19)  = ((ff15 + f2b4) * 0.5d0 - hh)        ! KNO3
       gama(20)  = ((ff15 + f2b1) * 0.5d0 - hh)        ! KCl
       gama(21)  = ((ff16 + f2b2) * 0.25d0- hh)*4.0d0  ! MgSO4
       gama(22)  = ((ff16 + f2b4) / 3.0d0 - hh)*2.0d0  ! Mg(NO3)2
       gama(23)  = ((ff16 + f2b1) / 3.0d0 - hh)*2.0d0  ! MgCl2
    
    !  ## Convert log(gama) coefficients to gama
       gama(1)  = max(-5.0d0, min(gama(1), 5.0d0))
       gama(1)  = 10.0**gama(1)
       gama(2)  = max(-5.0d0, min(gama(2), 5.0d0))
       gama(2)  = 10.0**gama(2)
       gama(3)  = max(-5.0d0, min(gama(3), 5.0d0))
       gama(3)  = 10.0**gama(3)
       gama(4)  = max(-5.0d0, min(gama(4), 5.0d0))
       gama(4)  = 10.0**gama(4)
       gama(5)  = max(-5.0d0, min(gama(5), 5.0d0))
       gama(5)  = 10.0**gama(5)
       gama(6)  = max(-5.0d0, min(gama(6), 5.0d0))
       gama(6)  = 10.0**gama(6)
       gama(7)  = max(-5.0d0, min(gama(7), 5.0d0))
       gama(7)  = 10.0**gama(7)
       gama(8)  = max(-5.0d0, min(gama(8), 5.0d0))
       gama(8)  = 10.0**gama(8)
       gama(9)  = max(-5.0d0, min(gama(9), 5.0d0))
       gama(9)  = 10.0**gama(9)
       gama(10) = max(-5.0d0, min(gama(10), 5.0d0))
       gama(10) = 10.0**gama(10)
       gama(11) = max(-5.0d0, min(gama(11), 5.0d0))
       gama(11) = 10.0**gama(11)
       gama(12) = max(-5.0d0, min(gama(12), 5.0d0))
       gama(12) = 10.0**gama(12)
       gama(13) = max(-5.0d0, min(gama(13), 5.0d0))
       gama(13) = 10.0**gama(13)
       gama(14) = 1.0d0
       gama(15) = max(-5.0d0, min(gama(15), 5.0d0))
       gama(15) = 10.0**gama(15)
       gama(16) = max(-5.0d0, min(gama(16), 5.0d0))
       gama(16) = 10.0**gama(16)
       gama(17) = max(-5.0d0, min(gama(17), 5.0d0))
       gama(17) = 10.0**gama(17)
       gama(18) = max(-5.0d0, min(gama(18), 5.0d0))
       gama(18) = 10.0**gama(18)
       gama(19) = max(-5.0d0, min(gama(19), 5.0d0))
       gama(19) = 10.0**gama(19)
       gama(20) = max(-5.0d0, min(gama(20), 5.0d0))
       gama(20) = 10.0**gama(20)
       gama(21) = max(-5.0d0, min(gama(21), 5.0d0))
       gama(21) = 10.0**gama(21)
       gama(22) = max(-5.0d0, min(gama(22), 5.0d0))
       gama(22) = 10.0**gama(22)
       gama(23) = max(-5.0d0, min(gama(23), 5.0d0))
       gama(23) = 10.0**gama(23)
       end if
    !
       return
    end subroutine mach_hetp_calcact4b
    
    
    
    !############################################################################
    ! ## HETP Code  
    ! ## Finds the smallest positive real root of a cubic equation analytically
    !
    ! ## Equation of the form: 
    ! ## x*x*x + a1*x*x + a2*x + a3 = 0.0d0
    ! ## Input:  a1, a2, a3
    ! ## Output: root, islv  ! Minimum positive real root 
    !                        ! root set to 1.0d30 if no root is found
    !                        ! islv > 0 if no root is found
    !
    ! ## Special case: quadratic equation solved separately
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_poly3v (a1, a2, a3, root, islv)
    !
       implicit none
    !
       real(kind=8),    intent (in) :: a1
       real(kind=8),    intent (in) :: a2
       real(kind=8),    intent (in) :: a3 
       real(kind=8),    intent(out) :: root
       integer(kind=4), intent(out) :: islv
    !
    !  ## Local variables
    !
       integer(kind=4)            :: j, ix
       real(kind=8)               :: thet, coef, s, sqd, ssig, tsig, t, d, q, r, pi
       real(kind=8)               :: x(3)
       real(kind=8), parameter    :: expon = 1.0d0/3.0d0, zero = 0.0d0,                &
                                     thet1 = 120.0d0/180.0d0, thet2 = 240.0d0/180.0d0, &
                                     eps = 1.d-50
    !
       pi = acos(-1.0d0)
       islv = 1
       root = 0.0d0
       x    = 0.0d0
    !
    !  #### 1. Quadratic equation
       d = a1*a1 - 4.0d0*a2
    
       if (abs(a3) <= eps) then
          ix   = 1
          x(1) = 0.0d0
       end if
    !
       if (abs(a3) <= eps .and. d >= zero) then
          ix   = 3
          x(2) = 0.5d0*(-a1 + sqrt(d))
          x(3) = 0.5d0*(-a1 - sqrt(d))
       end if
    !
    !  #### 2. Cubic equation
       q = (3.0d0*a2 - a1*a1)/9.0d0
       r = (9.0d0*a1*a2 - 27.0d0*a3 - 2.0d0*a1*a1*a1)/54.0d0
       d = q*q*q + r*r
    !
    !  ## Calculate roots *
    !  ## i. d < 0, three real roots
       if (abs(a3) > eps .and. d < -eps) then
          ix   = 3
          thet = expon*acos(r / sqrt(-q*q*q))
          coef = 2.0d0*sqrt(-q)
          x(1) = coef*cos(thet) - expon*a1
          x(2) = coef*cos(thet + thet1*pi) - expon*a1
          x(3) = coef*cos(thet + thet2*pi) - expon*a1
       end if
    !
    !  ## ii. d = 0, three real (one double) roots
       if (abs(a3) > eps .and. d >= -eps .and. d <= eps) then
          ix   = 2
          s    = (abs(r)**expon)*sign(1.0e0, real(r))
          x(1) = 2.0d0*s - expon*a1
          x(2) = -s - expon*a1
       end if
    !
    !  ## iii. d > 0, one real root
       if (abs(a3) > eps .and. d > eps) then
          ix   = 1
          sqd  = sqrt(d)
          ssig = sign(1.0e0, real(r + sqd))
          tsig = sign(1.0e0, real(r - sqd))
          s    = ssig*(abs(r + sqd))**expon
          t    = tsig*(abs(r - sqd))**expon
          x(1) = s + t - expon*a1
       end if
    
    !  #### Select the appropriate root
    !  ## Note that islv == 1 if there are no positive roots with a magnitude > eps 
       islv = 1
       do j = 1, 3
          islv = min(islv, int(0.5d0*(1.0d0 + sign(1.0d0, (eps - x(j))))))
       end do
    ! 
    !  ## Set roots with magnitude <= eps to 1.0d+30; then determine the smallest root 
       root = 1.0d30
       do j = 1, 3
          x(j) = max(x(j), 1.0d30*(0.5d0*(1.0d0 - sign(1.0d0, x(j) - eps))))
       end do
    !
       do j = 1, 3
          root = min(root, x(j))
       end do
    return
    end
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Finds the smallest positive real root of a cubic equation analytically;
    ! ## if no root is found analytically, an ITP search is then performed over
    ! ## the range (tiny2, cl)
    !
    ! ## Equation of the form:
    ! ## x*x*x + a1*x*x + a2*x + a3 = 0.0d0
    ! ## Input:  a1, a2, a3, cl
    ! ## Output: root, islv  ! Minimum positive real root 
    !                        ! root set to 1.0d30 if no root is found
    !                        ! islv > 0 if no root is found
    !
    ! ## Special case: quadratic equation solved separately
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_poly(a1, a2, a3, cl, root, islv)
    !
       implicit none
    !
       real(kind=8),    intent (in) :: a1
       real(kind=8),    intent (in) :: a2
       real(kind=8),    intent (in) :: a3 
       real(kind=8),    intent (in) :: cl
       real(kind=8),    intent(out) :: root
       integer(kind=4), intent(out) :: islv
    !
    ! Local variables
    !
       integer(kind=4), parameter :: jmax = 60
       integer(kind=4)            :: j, ix
       real(kind=8)               :: thet, coef, s, sqd, ssig, tsig, t, d, q, r, pi
       real(kind=8)               :: x(3)
       real(kind=8), parameter    :: expon = 1.0d0/3.0d0, zero = 0.0d0,                &
                                     thet1 = 120.0d0/180.0d0, thet2 = 240.0d0/180.0d0, &
                                     eps = 1.d-50
       real(kind=8)               :: y1, x1, y2, x2, nh, nmax, eps2!, k1, k2, 
       real(kind=8)               :: xf, xh, xt, sigma, delta, rr, x3, y3, gx2, gx, u1
       integer(kind=4)            :: k
       real(kind=8)    :: rooteval, ndiv, dx
       logical(kind=4) :: condition, ITPsearch, bisect
    !
       bisect = .true.
    !
       pi   = acos(-1.0d0)
       islv = 1
       root = 0.0d0
       x    = 0.0d0
    !
    !  #### 1. Quadratic equation
       d = a1*a1 - 4.0d0*a2
    !
       if (abs(a3) <= eps) then
          ix   = 1
          x(1) = 0.0d0
       end if
    !
       if (abs(a3) <= eps .and. d >= zero) then
          ix   = 3
          x(2) = 0.5d0*(-a1 + sqrt(d))
          x(3) = 0.5d0*(-a1 - sqrt(d))
       end if
    !
    !  #### 2. Cubic equation
       q = (3.0d0*a2 - a1*a1)/9.0d0
       r = (9.0d0*a1*a2 - 27.0d0*a3 - 2.0d0*a1*a1*a1)/54.0d0
       d = q*q*q + r*r
    !
    !  ## Calculate roots 
    !  ## i. d < 0, three real roots
       if (abs(a3) > eps .and. d < -eps) then
          ix   = 3
          thet = expon*acos(r / sqrt(-q*q*q))
          coef = 2.0d0*sqrt(-q)
          x(1) = coef*cos(thet) - expon*a1
          x(2) = coef*cos(thet + thet1*pi) - expon*a1
          x(3) = coef*cos(thet + thet2*pi) - expon*a1
       end if
    !
    !  ## ii. d = 0, three real (one double) roots
       if (abs(a3) > eps .and. d >= -eps .and. d <= eps) then
          ix   = 2
          s    = (abs(r)**expon)*sign(1.0e0, real(r))
          x(1) = 2.0d0*s - expon*a1
          x(2) = -s - expon*a1
       end if
    !
    !  ## iii. d > 0, one real root
       if (abs(a3) > eps .and. d > eps) then
          ix   = 1
          sqd  = sqrt(d)
          ssig = sign(1.0e0, real(r + sqd))
          tsig = sign(1.0e0, real(r - sqd))
          s    = ssig*(abs(r + sqd))**expon
          t    = tsig*(abs(r - sqd))**expon
          x(1) = s + t - expon*a1
       end if
    
    !  #### Select the appropriate root
    !  ## Note that islv == 1 if there are no positive roots with a magnitude > eps 
       islv = 1
       do j = 1, 3
          islv = min(islv, int(0.5d0*(1.0d0 + sign(1.0d0, (eps - x(j))))))
       end do
    ! 
    !  ## Set roots with magnitude <= eps to 1.0d+30; then determine the smallest root 
       root = 1.0d30
       do j = 1, 3
          x(j) = max(x(j), 1.0d30*(0.5d0*(1.0d0 - sign(1.0d0, x(j) - eps))))
       end do
    !
       do j = 1, 3
          root = min(root, x(j))
       end do
    !
       ITPsearch = .false.
    !  ## No root less than 1.0d29 was found; perform ITP search instead
       if (root > 1.0d29) then  
          rooteval = 0
          ndiv     = 5
          x1       = 1.0d-28           !Lower bound is tiny2
          y1       = x1*x1*x1 + a1*x1*x1 + a2*x1 + a3
          dx       = cl/real(ndiv-1)   !Upper bound is total chlorine 
          condition= .false.
    !
    ! ## Subdivision search for a root
          do while (rooteval < 5 .and. ( .not. condition))
             x2 = x1 + dx 
             y2 = x2*x2*x2 + a1*x2*x2 + a2*x2 + a3
    !
    ! ## Test for a sign change
             if (y1 == 0.0d0) then
    ! ## 1. x1 is a root
                condition = .true.
                ITPsearch = .false.
                root      = x1
                islv      = 0
             elseif (y2 == 0.0d0) then 
    ! ## 2. x2 is a root
                condition = .true.
                ITPsearch = .false.
                root      = x2
                islv      = 0
             elseif (y1*y2 < 0.0d0) then
    ! ## 3. Interval has a sign change
                condition = .true.
                ITPsearch = .true.
             elseif (y1*y2 > 0.0d0) then
    ! ## 4. No sign change; continue onto next interval
                x1 = x2
                y1 = y2
                ITPsearch = .false.
             end if
    !
             rooteval = rooteval + 1
          end do
       end if 
    !
    !   k1 = 0.1d0
    !   k2 = 2.0d0
       eps2 = 1.0d-9
    !
    !  #### Perform an ITP search to refine the root, if a root exists
       if (bisect .and. ITPsearch) then
          gx   = x2 - x1
          gx2  = (x1+x2)*0.5d0
          u1   = 0.2d0/gx
          nh   = log10(abs(gx/(2.0d0*gx2*eps2))) / 0.301029995663981195d0
          nmax = nint(nh) + 2
          k    = 0 
    !
          do while (abs(x2-x1) > x2*eps2 .and. k < jmax)
             gx   = x2 - x1
             xh   = 0.5d0*(x1 + x2)
             rr   = max(gx2*eps*2.d0**(real(nmax - k)) - 0.5d0*gx, 0.0d0)
             delta= u1*(max(gx, 0.0d0))**2.0d0   
             xf   = max((y2*x1 - y1*x2) / (y2 - y1), 0.0d0)
             sigma = sign(1.0d0, xh - xf)
    !
             if (delta <= abs(xh-xf)) then
                xt = xf + sigma*delta
             else
                xt = xh
             end if 
    !
             if (abs(xt - xh) <= rr) then
                x3 = xt
             else
                x3 = xh - sigma*rr
             end if
    !
             y3 = x3*x3*x3 + a1*x3*x3 + a2*x3 + a3
    !
             if (y3 > 0.0d0) then
                x2 = x3
                y2 = y3
             elseif (y3 < 0.0d0) then
                x1 = x3
                y1 = y3
             else
                x1 = x3
                x2 = x3
             end if
    !    
             k = k + 1
          end do
    !
          root = x3
          islv = 0
       end if
       return
    end subroutine mach_hetp_poly
    
    
    
    !############################################################################
    ! ## HETP Code 
    ! ## Adjusts to force mass balance for volatile species and sulfate
    ! ## Calculate EXCESS mass only, remove excess, first from aerosol (liquid) 
    ! ## phase, second from the gas phase, and third from the solid phase 
    !
    ! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
    ! ## Written by Stefan Miller
    !
    ! ## Code is based on ISORROPIA II, obtained from the CMAQ air-quality
    ! ## model (https://github.com/USEPA/CMAQ/tree/main/CCTM/src/aero/aero6)
    !############################################################################
    subroutine mach_hetp_adjust(so4, no3, nh4, cl, so4_t, hso4, no3_t, ghno3,   &
                                nh4_t, gnh3, cl_t, ghcl, caso4)
    !
       use mach_hetp_mod,         only: tiny2
       implicit none
    !
       real(kind=8),    intent   (in) :: so4   
       real(kind=8),    intent   (in) :: nh4   
       real(kind=8),    intent   (in) :: no3   
       real(kind=8),    intent   (in) :: cl    
       real(kind=8),    intent(inout) :: so4_t 
       real(kind=8),    intent(inout) :: hso4  
       real(kind=8),    intent(inout) :: caso4 
       real(kind=8),    intent(inout) :: no3_t 
       real(kind=8),    intent(inout) :: ghno3 
       real(kind=8),    intent(inout) :: nh4_t 
       real(kind=8),    intent(inout) :: gnh3  
       real(kind=8),    intent(inout) :: cl_t  
       real(kind=8),    intent(inout) :: ghcl  
    !
    !  ## Local variables:
       real(kind=8) :: exnh4, exno3, exs4, excl
    !
       exnh4 = 0.0d0
       exno3 = 0.0d0
       exs4  = 0.0d0
       excl  = 0.0d0
    !
    !  ## Calculate excess as: solution - input (units mol/m3 for all species)
    !  1. Expected: NH4+(aq) + NH3(g)  - TA = 0.0d0
          exnh4  = max(nh4_t + gnh3  - nh4, 0.0d0)
    !
    !  2. Expected: NO3-(aq) + HNO3(g) - TN = 0.0d0
          exno3  = max(no3_t + ghno3 - no3, 0.0d0)
    !
    !  3. Expected: Cl-(aq)  + HCl(g)  - TCl  = 0.0d0
          excl   = max(cl_t  + ghcl  - cl , 0.0d0)
    !
    !  4. Expected: SO4--(aq) + HSO4-(aq) + CaSO4(s) - TS = 0.0d0
          exs4   = max(so4_t + hso4  + caso4 - so4, 0.0d0)
    !
    !
    !  ## 1. Adjust ammonium
       if (exnh4 >= tiny2) then
          if (nh4_t > exnh4) then
             nh4_t = max(nh4_t - exnh4, 0.0d0)
          else
             exnh4 = max(exnh4 - nh4_t, 0.0d0)
             nh4_t = 0.0d0
             gnh3  = max(gnh3  - exnh4, 0.0d0)
          end if
       end if
    !
    !  ## 2. Adjust nitrate
       if (exno3 >= tiny2) then
          if (no3_t > exno3) then
             no3_t = max(no3_t - exno3, 0.0d0)
          else
             exno3 = max(exno3 - no3_t, 0.0d0)
             no3_t = 0.0d0
             ghno3 = max(ghno3 - exno3, 0.0d0)
          end if
       end if
    !
    !  ## 3. Adjust chloride
       if (excl >= tiny2) then
          if (cl_t > excl) then
             cl_t = max(cl_t - excl, 0.0d0)
          else
             excl = max(excl - cl_t, 0.0d0)
             cl_t = 0.0d0
             ghcl = max(ghcl - excl, 0.0d0)
          end if
       end if
    !
    !  ## 4. Adjust sulfate
       if (exs4 >= tiny2) then
          if (hso4 > exs4) then
             hso4 = max(hso4 - exs4, 0.0d0)
          else
             exs4 = max(exs4 - hso4, 0.0d0)
             hso4 = 0.0d0
    !
             if (so4_t > exs4) then
                so4_t = max(so4_t - exs4, 0.0d0)
             else
                exs4  = max(exs4 - so4_t, 0.0d0)
                so4_t = 0.0d0
    !
                if (caso4 > exs4) then  
                   caso4 = max(caso4 - exs4, 0.0d0)
                else
                   exs4  = max(exs4 - caso4, 0.0d0)
                   caso4 = 0.0d0
                end if
             end if
          end if
       end if
    !
       return
    end subroutine mach_hetp_adjust
end module hetp_mod
