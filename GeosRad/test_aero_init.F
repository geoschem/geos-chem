      program test_AERO_INIT
    
        parameter (nwv=61,nr=7,nspec=6)
        parameter (nbndlw=16,nbndsw=14)
        real*8    :: lambda(nwv)
        real*8   :: rrtmg_lmb(nbndlw+nbndsw)

! RRTMG locations for aerosol OD 
        REAL*8  :: aero_band_sel(nbndlw+nbndsw) 
        INTEGER  :: id_aer_lmb0 (nbndlw+nbndsw)
        INTEGER  :: id_aer_lmb1 (nbndlw+nbndsw)

       ! Local variables
       integer  :: ib, iflag

      data aero_band_sel /180.,425.,565.,665.,760.,900.,1030.,1130.,
     &                   1285.,1435.,1640.,1940.,2165.,2315.,2490.,
     &                   2925., 2903.,3601.,4310.,4892.,5623.,6872.,
     &                   7872., 10590.,14420.,18970.,25015.,30390.,
     &                   43507.,1412./

      data lambda /250, 300, 350, 400, 450, 500, 550, 600, 650, 700,
     &  750, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3200, 
     & 3390, 3500, 3750, 4000, 4500, 5000, 5500, 6000, 6200, 6500, 7200, 
     & 7900, 8200, 8500, 8700, 9000, 9200, 9500, 9800, 10000, 10600, 
     &11000, 11500, 12500, 13000, 14000, 14800, 15000, 16400, 17200, 
     & 18000, 18500, 20000, 21300, 22500, 25000, 27900, 
     & 30000, -30536, -25536/
    
       rrtmg_lmb = 1.0e7/aero_band_sel   ! convert from cm-1 to nm

       ! Find bracketing GC aerosol indices for each RRTMG band
       do ib=nbndlw+nbndsw,1,-1
          iflag = 0 
          if(rrtmg_lmb(ib).lt.lambda(1)) iflag = -1
          if(rrtmg_lmb(ib).gt.lambda(59)) iflag = 1 
          select case (iflag)
       !Handle short wavelength points
             case (-1)
                id_aer_lmb0(ib) = 1 
                id_aer_lmb1(ib) = 2 

       !Handle long wavelength points
             case (1) 
                id_aer_lmb0(ib) = 58
                id_aer_lmb1(ib) = 59
    
             case default
                ip = 1 
                do while (rrtmg_lmb(ib).gt. lambda(ip)) 
                   ip=ip+1
                end do
                id_aer_lmb0(ib) =  ip-1
                id_aer_lmb1(ib) =  ip  
             end select
       end do

       do ib=1,nbndlw+nbndsw
          print *, rrtmg_lmb(ib),lambda(id_aer_lmb0(ib)),
     &             lambda(id_aer_lmb1(ib))
       end do
   
       end
