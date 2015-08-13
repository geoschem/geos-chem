!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: processData
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE processData
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  
  IMPLICIT NONE
!
! !REVISION HISTORY:
!  09 Jul 2010 - C. Caroauge - Initial version
!  01 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  integer   :: i
  integer, dimension(26)    :: st
  character(len=f_len_output+14) :: filename1

  st(:)=0.0e+0_fp
!!!!!monthly sum of total fluxes for all veg. pixels
!!!!! UNITS = g
  NPPmonthly(yr,mo)    = SUM(gridAreab(:,1)*(frac_tree(:,1)+frac_herb(:,1))*NPP(:,mo))
  respmonthly(yr,mo)   = SUM(gridAreab(:,1)*((frac_tree(:,1)*wresp(:,1))+(frac_herb(:,1)*hresp(:,1))))
  combmonthly(yr,mo)   = SUM(gridAreab(:,1)*((frac_tree(:,1)*wcomb(:,1))+(frac_herb(:,1)*hresp(:,1))))
  herbmonthly(yr,mo)   = SUM(gridAreab(:,1)*((frac_tree(:,1)*wherb(:,1))+(frac_herb(:,1)*hherb(:,1))))
  biofmonthly(yr,mo)   = SUM(gridAreab(:,1)*wbiof(:,1)*frac_tree(:,1))
  
  
  respmonthly_hg(yr,mo)= SUM(gridAreab(:,1)*((frac_tree(:,1)*wresp_hg(:,1))+(frac_herb(:,1)*hresp_hg(:,1))))
  combmonthly_hg(yr,mo)=SUM(gridAreab(:,1)*((frac_tree(:,1)*wcomb_hg(:,1))+(frac_herb(:,1)*hcomb_hg(:,1))))
  reemmonthly_hg(yr,mo)= SUM(gridAreab(:,1)*reemitted(:,1))
  photmonthly_hg(yr,mo)= SUM(gridAreab(:,1)*photoreduced(:,1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!monthly fluxes of respiration by each soil class
      ! these have already been weighted by f_tree and f_herb (see doTreeCarbonHg and doHerbCarbonHg)
  resp_surfstr(:,mo)=resp_surfstr(:,mo)*gridAreab(:,1) !units = g
  resp_surfmet(:,mo)=resp_surfmet(:,mo)*gridAreab(:,1)
  resp_surfmic(:,mo)=resp_surfmic(:,mo)*gridAreab(:,1)
  resp_soilstr(:,mo)=resp_soilstr(:,mo)*gridAreab(:,1)
  resp_soilmet(:,mo)=resp_soilmet(:,mo)*gridAreab(:,1)
  resp_soilmic(:,mo)=resp_soilmic(:,mo)*gridAreab(:,1)
  resp_slow(:,mo)=resp_slow(:,mo)*gridAreab(:,1)
  resp_armored(:,mo)=resp_armored(:,mo)*gridAreab(:,1)
  
  resp_surfstr_hg(:,mo)=resp_surfstr_hg(:,mo)*gridAreab(:,1) !units = g
  resp_surfmet_hg(:,mo)=resp_surfmet_hg(:,mo)*gridAreab(:,1)
  resp_surfmic_hg(:,mo)=resp_surfmic_hg(:,mo)*gridAreab(:,1)
  resp_soilstr_hg(:,mo)=resp_soilstr_hg(:,mo)*gridAreab(:,1)
  resp_soilmet_hg(:,mo)=resp_soilmet_hg(:,mo)*gridAreab(:,1)
  resp_soilmic_hg(:,mo)=resp_soilmic_hg(:,mo)*gridAreab(:,1)
  resp_slow_hg(:,mo)=resp_slow_hg(:,mo)*gridAreab(:,1)
  resp_armored_hg(:,mo)=resp_armored_hg(:,mo)*gridAreab(:,1)
  
  Hg_pool_fluxes1(:,mo)=Hg_pool_fluxes1(:,mo)*gridAreab(:,1) !units = g
  Hg_pool_fluxes2(:,mo)=Hg_pool_fluxes2(:,mo)*gridAreab(:,1) !units = g
  Hg_pool_fluxes3(:,mo)=Hg_pool_fluxes3(:,mo)*gridAreab(:,1) !units = g
  Hg_pool_fluxes4(:,mo)=Hg_pool_fluxes4(:,mo)*gridAreab(:,1) !units = g
  Hg_pool_fluxes5(:,mo)=Hg_pool_fluxes5(:,mo)*gridAreab(:,1) !units = g
  Hg_pool_fluxes6(:,mo)=Hg_pool_fluxes6(:,mo)*gridAreab(:,1) !units = g
  
  
!!!!!!!!!!Monthly size of pools, units=g
  slowmonthly(yr,mo)   = SUM(gridAreab(:,1)*(frac_tree(:,1)*slowpool(:,1)+frac_herb(:,1)*hslowpool(:,1)))
  armoredmonthly(yr,mo)= SUM(gridAreab(:,1)*(frac_tree(:,1)*armoredpool(:,1)+frac_herb(:,1)*harmoredpool(:,1)))
  surfstrmonthly(yr,mo)= SUM(gridAreab(:,1)*(frac_tree(:,1)*surfstrpool(:,1)+frac_herb(:,1)*hsurfstrpool(:,1)))
  surfmetmonthly(yr,mo)= SUM(gridAreab(:,1)*(frac_tree(:,1)*surfmetpool(:,1)+frac_herb(:,1)*hsurfmetpool(:,1)))
  surfmicmonthly(yr,mo)= SUM(gridAreab(:,1)*(frac_tree(:,1)*surfmicpool(:,1)+frac_herb(:,1)*hsurfmicpool(:,1)))
  soilstrmonthly(yr,mo)= SUM(gridAreab(:,1)*(frac_tree(:,1)*soilstrpool(:,1)+frac_herb(:,1)*hsoilstrpool(:,1)))
  soilmetmonthly(yr,mo)= SUM(gridAreab(:,1)*(frac_tree(:,1)*soilmetpool(:,1)+frac_herb(:,1)*hsoilmetpool(:,1)))
  soilmicmonthly(yr,mo)= SUM(gridAreab(:,1)*(frac_tree(:,1)*soilmicpool(:,1)+frac_herb(:,1)*hsoilmicpool(:,1)))
  leafmonthly(yr,mo)   = SUM(gridAreab(:,1)*(frac_tree(:,1)*leafpool(:,1)+frac_herb(:,1)*hleafpool(:,1)))   
!!!!!!!Monthly size of Hg pools, units=g
  slowmonthly_hg(yr,mo)   =SUM(gridAreab(:,1)*(frac_tree(:,1)*slowpool_hg(:,1)&
       & +frac_herb(:,1)*hslowpool_hg(:,1)))
  armoredmonthly_hg(yr,mo)=SUM(gridAreab(:,1)*(frac_tree(:,1)*armoredpool_hg(:,1)&
       & +frac_herb(:,1)*harmoredpool_hg(:,1)))
  surfstrmonthly_hg(yr,mo)=SUM(gridAreab(:,1)*(frac_tree(:,1)*surfstrpool_hg(:,1)&
       & +frac_herb(:,1)*hsurfstrpool_hg(:,1)))
  surfmetmonthly_hg(yr,mo)=SUM(gridAreab(:,1)*(frac_tree(:,1)*surfmetpool_hg(:,1)&
       & +frac_herb(:,1)*hsurfmetpool_hg(:,1)))
  surfmicmonthly_hg(yr,mo)=SUM(gridAreab(:,1)*(frac_tree(:,1)*surfmicpool_hg(:,1)&
       & +frac_herb(:,1)*hsurfmicpool_hg(:,1)))
  soilstrmonthly_hg(yr,mo)=SUM(gridAreab(:,1)*(frac_tree(:,1)*soilstrpool_hg(:,1)&
       & +frac_herb(:,1)*hsoilstrpool_hg(:,1)))
  soilmetmonthly_hg(yr,mo)=SUM(gridAreab(:,1)*(frac_tree(:,1)*soilmetpool_hg(:,1)&
       & +frac_herb(:,1)*hsoilmetpool_hg(:,1)))
  soilmicmonthly_hg(yr,mo)= SUM(gridAreab(:,1)*(frac_tree(:,1)*soilmicpool_hg(:,1)&
       & +frac_herb(:,1)*hsoilmicpool_hg(:,1)))
  leafmonthly_hg(yr,mo)   = SUM(gridAreab(:,1)*(frac_tree(:,1)*leafpool_hg(:,1)&
       & +frac_herb(:,1)*hleafpool_hg(:,1)))   
  HgAqmonthly(yr,mo)      = SUM(gridAreab(:,1)*(frac_tree(:,1)*HgAq(:,1)&
       & +frac_herb(:,1)*hHgAq(:,1)))
  total_tree_hg(:,1)=slowpool_hg(:,1)+armoredpool_hg(:,1)+surfstrpool_hg(:,1)&
       &+surfmetpool_hg(:,1)+surfmicpool_hg(:,1)+soilstrpool_hg(:,1)+soilmetpool_hg(:,1)&
       &+soilmicpool_hg(:,1)
  total_herb_hg(:,1)=hslowpool_hg(:,1)+harmoredpool_hg(:,1)+hsurfstrpool_hg(:,1)&
       &+hsurfmetpool_hg(:,1)+hsurfmicpool_hg(:,1)+hsoilstrpool_hg(:,1)+hsoilmetpool_hg(:,1)&
       &+hsoilmicpool_hg(:,1)
  
  IF (mo .eq. 12) THEN
     DO i=1, n_veg
        IF (veg1(i,1) .eq. 1) THEN
           biomeAnnual_Hg(yr,1)=biomeAnnual_Hg(yr,1)&
                &+(((frac_tree(i,1)*total_tree_hg(i,1))+(frac_herb(i,1)*total_herb_hg(i,1)))*gridAreab(i,1))
        ENDIF
        IF (veg1(i,1) .eq. 3) THEN
           biomeAnnual_Hg(yr,2)=biomeAnnual_Hg(yr,2)&
                &+(((frac_tree(i,1)*total_tree_hg(i,1))+(frac_herb(i,1)*total_herb_hg(i,1)))*gridAreab(i,1))
        ENDIF
        IF (veg1(i,1) .eq. 4) THEN
           biomeAnnual_Hg(yr,3)=biomeAnnual_Hg(yr,3)&
                &+(((frac_tree(i,1)*total_tree_hg(i,1))+(frac_herb(i,1)*total_herb_hg(i,1)))*gridAreab(i,1))
        ENDIF
        IF (veg1(i,1) .eq. 6) THEN
           biomeAnnual_Hg(yr,4)=biomeAnnual_Hg(yr,4)&
                &+(((frac_tree(i,1)*total_tree_hg(i,1))+(frac_herb(i,1)*total_herb_hg(i,1)))*gridAreab(i,1))
        ENDIF
        IF (veg1(i,1) .eq. 10) THEN
           biomeAnnual_Hg(yr,5)=biomeAnnual_Hg(yr,5)&
                &+(((frac_tree(i,1)*total_tree_hg(i,1))+(frac_herb(i,1)*total_herb_hg(i,1)))*gridAreab(i,1))
        ENDIF
     END DO
  ENDIF
  
  IF ((yr .eq. HgPOOLSequilibriumYear) .or. (yr .eq. preindYear) .or. (yr .eq. indYear))THEN
!$OMP PARALLEL         &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
     respEQ(:,mo)=(frac_tree(:,1)*wresp(:,1)+frac_herb(:,1)*hresp(:,1))
     combEQ(:,mo)=(frac_tree(:,1)*wcomb(:,1)+frac_herb(:,1)*hcomb(:,1))
     herbEQ(:,mo)=(frac_tree(:,1)*wherb(:,1)+frac_herb(:,1)*hherb(:,1))
     biofEQ(:,mo)=wbiof(:,1)!*frac_tree(:,1)
     
     respEQ_hg(:,mo)=(frac_tree(:,1)*wresp_hg(:,1))+(frac_herb(:,1)*hresp_hg(:,1))
     combEQ_hg(:,mo)=(frac_tree(:,1)*wcomb_hg(:,1))+(frac_herb(:,1)*hcomb_hg(:,1))
     reemitEQ_hg(:,mo)=reemitted(:,1)
     photoEQ_hg(:,mo)=photoreduced(:,1)
     leafpoolEQ(:,mo)   =(frac_tree(:,1)*leafpool(:,1))    + (frac_herb(:,1)*hleafpool(:,1))
     slowpoolEQ(:,mo)   =(frac_tree(:,1)*slowpool(:,1))    + (frac_herb(:,1)*hslowpool(:,1))
     armoredpoolEQ(:,mo)=(frac_tree(:,1)*armoredpool(:,1)) + (frac_herb(:,1)*harmoredpool(:,1))
     surfstrpoolEQ(:,mo)=(frac_tree(:,1)*surfstrpool(:,1)) + (frac_herb(:,1)*hsurfstrpool(:,1))
     soilstrpoolEQ(:,mo)=(frac_tree(:,1)*soilstrpool(:,1)) + (frac_herb(:,1)*hsoilstrpool(:,1))
     surfmetpoolEQ(:,mo)=(frac_tree(:,1)*surfmetpool(:,1)) + (frac_herb(:,1)*hsurfmetpool(:,1))
     soilmetpoolEQ(:,mo)=(frac_tree(:,1)*soilmetpool(:,1)) + (frac_herb(:,1)*hsoilmetpool(:,1))
     surfmicpoolEQ(:,mo)=(frac_tree(:,1)*surfmicpool(:,1)) + (frac_herb(:,1)*hsurfmicpool(:,1))
     soilmicpoolEQ(:,mo)=(frac_tree(:,1)*soilmicpool(:,1)) + (frac_herb(:,1)*hsoilmicpool(:,1))
     
     leafpoolEQ_hg(:,mo)   =(frac_tree(:,1)*leafpool_Hg(:,1))    + (frac_herb(:,1)*hleafpool_Hg(:,1))
     slowpoolEQ_hg(:,mo)   =(frac_tree(:,1)*slowpool_Hg(:,1))    + (frac_herb(:,1)*hslowpool_Hg(:,1))
     armoredpoolEQ_hg(:,mo)=(frac_tree(:,1)*armoredpool_Hg(:,1)) + (frac_herb(:,1)*harmoredpool_Hg(:,1))
     surfstrpoolEQ_hg(:,mo)=(frac_tree(:,1)*surfstrpool_Hg(:,1)) + (frac_herb(:,1)*hsurfstrpool_Hg(:,1))
     soilstrpoolEQ_hg(:,mo)=(frac_tree(:,1)*soilstrpool_Hg(:,1)) + (frac_herb(:,1)*hsoilstrpool_Hg(:,1))
     surfmetpoolEQ_hg(:,mo)=(frac_tree(:,1)*surfmetpool_Hg(:,1)) + (frac_herb(:,1)*hsurfmetpool_Hg(:,1))
     soilmetpoolEQ_hg(:,mo)=(frac_tree(:,1)*soilmetpool_Hg(:,1)) + (frac_herb(:,1)*hsoilmetpool_Hg(:,1))
     surfmicpoolEQ_hg(:,mo)=(frac_tree(:,1)*surfmicpool_Hg(:,1)) + (frac_herb(:,1)*hsurfmicpool_Hg(:,1))
     soilmicpoolEQ_hg(:,mo)=(frac_tree(:,1)*soilmicpool_Hg(:,1)) + (frac_herb(:,1)*hsoilmicpool_Hg(:,1))
     HgAqEQ_hg(:,mo)          =(frac_tree(:,1)*HgAq(:,1))           + (frac_herb(:,1)*hHgAq(:,1))
!$OMP END WORKSHARE
!$OMP END PARALLEL

     IF (mo .eq. 12) THEN
        ! write data to files
        filename1(1:f_len_output)=outputpath
        filename1(f_len_output+1:f_len_output+10)='gridareaea'
        IF (yr .eq. preindYear) THEN
           filename1(f_len_output+11:f_len_output+14)='_pre'
        ENDIF
        IF (yr .eq. HgPoolsequilibriumYear) THEN
           filename1(f_len_output+11:f_len_output+14)='_pst'
        ENDIF
        IF (yr .eq. indYear) THEN
           filename1(f_len_output+11:f_len_output+14)='_ind'
        ENDIF
        OPEN(UNIT=5, file=filename1)
        WRITE(5,*) gridareab
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='leafpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) leafpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='heafpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) hleafpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='slowpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) slowpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='hlowpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) hslowpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='armdpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) armoredpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='hrmdpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) harmoredpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='tUSTpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfstrpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='hUSTpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) hsurfstrpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='tOSTpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilstrpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='hOSTpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) hsoilstrpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='tUMTpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmetpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='hUMTpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) hsurfmetpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='tOMTpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmetpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='hOMTpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) hsoilmetpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='tUMIpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmicpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='hUMIpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) hsurfmicpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='tOMIpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmicpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='hOMIpoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) hsoilmicpool_Hg
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='tHgApoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) HgAq
        CLOSE(5) 
        filename1(f_len_output+1:f_len_output+10)='hHgApoolHg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) hHgAq
        CLOSE(5) 
        
        filename1(f_len_output+1:f_len_output+10)='NPPmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) NPPmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='RESmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) respmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='COMmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) combmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='HERmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) herbmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='BIFmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) biofmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='RES_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) respmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='COM_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) combmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='REM_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) reemmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='PHT_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) photmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='SLWmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) slowmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='ARMmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) armoredmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='USTmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfstrmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UMEmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmetmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UMImonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmicmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OSTmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilstrmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OMEmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmetmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OMImonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmicmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='LEAmonthly'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) leafmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='SLW_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) slowmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='ARM_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) armoredmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UST_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfstrmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UME_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmetmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UMI_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmicmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OST_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilstrmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OME_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmetmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OMI_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmicmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='LEA_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) leafmonthly_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='HgA_Hg_mon'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) HgAqmonthly
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UST_Hg_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_surfstr_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UME_Hg_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_surfmet_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UMI_Hg_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_surfmic_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OST_Hg_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_soilstr_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OME_Hg_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_soilmet_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OMI_Hg_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_soilmic_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='SLW_Hg_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_slow_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='ARM_Hg_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_armored_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UST_CC_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_surfstr
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UME_CC_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_surfmet
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UMI_CC_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_surfmic
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OST_CC_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_soilstr
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OME_CC_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_soilmet
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OMI_CC_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_soilmic
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='SLW_CC_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_slow
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='ARM_CC_res'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) resp_armored
        CLOSE(5)
        
        filename1(f_len_output+1:f_len_output+10)='RESP_EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) respEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='COMB_EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) combEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='HERB_EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) herbEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='BIOF_EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) biofEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='RESP_EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) respEQ_Hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='COMB_EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) combEQ_Hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='REEM_EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) reemitEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='PHOT_EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) photoEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='LEAF_EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) leafpoolEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='SLOW_EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) slowpoolEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='ARMD_EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) armoredpoolEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UST__EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfstrpoolEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OST__EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilstrpoolEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UME__EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmetpoolEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OME__EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmetpoolEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UMI__EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmicpoolEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OMI__EQ_CC'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmicpoolEQ
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='LEAF_EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) leafpoolEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='SLOW_EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) slowpoolEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='ARMD_EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) armoredpoolEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UST__EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfstrpoolEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OST__EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilstrpoolEQ_Hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UME__EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmetpoolEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OME__EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmetpoolEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='UMI__EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) surfmicpoolEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='OMI__EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) soilmicpoolEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='HgAq_EQ_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) HgAqEQ_hg
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='Pool1_flux'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) Hg_pool_fluxes1
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='Pool2_flux'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) Hg_pool_fluxes2
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='Pool3_flux'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) Hg_pool_fluxes3
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='Pool4_flux'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) Hg_pool_fluxes4
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='Pool5_flux'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) Hg_pool_fluxes5
        CLOSE(5)
        filename1(f_len_output+1:f_len_output+10)='Pool6_flux'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) Hg_pool_fluxes6
        CLOSE(5)
        
        filename1(f_len_output+1:f_len_output+10)='biome_an_Hg'
        OPEN(UNIT=5, file=filename1, STATUS="NEW")
        WRITE(5,*) biomeAnnual_Hg
        CLOSE(5)
        
     ENDIF
  ENDIF
END SUBROUTINE processData
!EOC
