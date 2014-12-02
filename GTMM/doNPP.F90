!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doNPP
!
! !DESCRIPTION: Subroutin doNPP calculate net primary production (NPP)
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doNPP
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  
  implicit none
!
! !REVISION HISTORY:
!  09 Jul 2010 - C. Carouge  - Parallelization
!  01 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER :: i
  character(len=f_len_output+4) :: filename3
  
  filename3(1:f_len_output)=outputpath
  
  !T1 is the temperature effect relating for each site the 
  !optimal temperature for highest NDVI (TOPT) to the temp
  !considered to be generally the most optimal for 
  !photosynthesis, 20C.  T2 is the effect of variable air
  !temperature throughout the year around TOPT

!$OMP PARALLEL         &
!$OMP DEFAULT(SHARED)  
!$OMP DO PRIVATE(i)
  DO i=1, n_veg
     T1(i,1)=0.800e+0_fp+(0.0200e+0_fp*topt(i,1))-(0.0005e+0_fp &
        *(topt(i,1)**2.000e+0_fp))
     T2low(i,1)=1.00e+0_fp/(1.00e+0_fp+exp(0.200e+0_fp*(topt(i,1) &
        -10.000e+0_fp-airt1(i,mo))))
     T2high(i,1)=1.00e+0_fp/(1.00e+0_fp+exp(0.300e+0_fp*(airt1(i,mo) &
        -10.000e+0_fp-topt(i,1))))
     NPPtemp(i,1)=T1(i,1)*(1.1919e+0_fp*T2low(i,1)*T2high(i,1))
     IF (T1(i,1)     .lt. 0e+0_fp .or.                                          &
         T2low(i,1)  .lt. 0e+0_fp .or.                                          &
         T2high(i,1) .lt. 0e+0_fp) THEN 
        NPPtemp(i,1)=0.00e+0_fp
     ELSE IF (airt1(i,mo) .lt. -10e+0_fp) THEN
        NPPtemp(i,1)=0.00e+0_fp
     END IF

     epsilona(i,1)=EMAX*NPPtemp(i,1)*NPPmoist(i,mo)
     IPAR(i,1)=FPAR(i,mo)*solrad1(i,mo)*1.3148333e+0_fp !solarconversion
     NPP(i,mo)=epsilona(i,1)*IPAR(i,1)
  END DO
!$OMP END DO

!$OMP WORKSHARE
  !derive abiotic effect for each timestep in the cycle
  bgtemp(:,mo)=Q10**((airt1(:,mo)-30.000e+0_fp)/10.000e+0_fp)
  abiotic(:,mo)=bgmoist(:,mo)*bgtemp(:,mo)
  
  WHERE (abiotic(:,mo) > 1e+0_fp)
     abiotic(:,mo)=1.0e+0_fp
  END WHERE

!$OMP END WORKSHARE
!$OMP END PARALLEL


  IF ( yr == NPPequilibriumYear .AND. mo == 12 ) THEN
     filename3 = outputpath // 'NPP'
     OPEN(UNIT=4, FILE=filename3, FORM="FORMATTED", STATUS="NEW")
     WRITE(4,*) NPP
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4) = 'ABIO'
     OPEN(UNIT=4, FILE=filename3, FORM="FORMATTED", STATUS="NEW")
     WRITE(4,*) abiotic
     CLOSE(4)
     
  ENDIF
  
  IF (yr .eq. 1 .and. mo .eq. 12) THEN
     filename3(f_len_output+1:f_len_output+4)='bgtp'
     OPEN(UNIT=4, FILE=filename3, FORM="FORMATTED",&
          STATUS="NEW")
     WRITE(4,*) bgtemp
     CLOSE(4)
  ENDIF
END SUBROUTINE doNPP
!EOC
