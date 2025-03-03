!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getSoilMoistParams
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE getSoilMoistParams
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
!  09 Jul 2010 - C. Carouge  - Parallelization
!  01 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP  
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
  REAL(fp) :: wp(n_veg, 1)
  REAL(fp) :: fc(n_veg, 1)
  REAL(fp) :: ssc(n_veg, 1)
  REAL(fp) :: alpha(n_veg, 1)
  REAL(fp) :: beta(n_veg, 1)
  REAL(fp) :: text(n_veg, 1)
  INTEGER :: i, j
  
  text=soiltext1
  
!$OMP PARALLEL DO     &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i)
  DO i=1,n_veg
     IF (soiltext1(i,1) .eq. 7) THEN
        text(i,1)=6
     ENDIF
     !for grassland biomes (7, 10, 12)
     IF (veg1(i,1) == 7 .or. veg1(i,1) == 10 .or. veg1(i,1) == 12) THEN
        IF (text(i,1) == 1) THEN
           wp(i,1)=129.000e+0_fp
           fc(i,1)=232.000e+0_fp
           ssc(i,1)=455.000e+0_fp
           alpha(i,1)=0.0020e+0_fp
           beta(i,1)=-6.540e+0_fp
        ELSE IF (text(i,1) == 2) THEN
           wp(i,1)=81.000e+0_fp
           fc(i,1)=202.000e+0_fp
           ssc(i,1)=394.000e+0_fp
           alpha(i,1)=0.002e+0_fp
           beta(i,1)=-5.480e+0_fp
        ELSE IF (text(i,1) == 3) THEN
           wp(i,1)=129.000e+0_fp
           fc(i,1)=232.000e+0_fp
           ssc(i,1)=455.000e+0_fp
           alpha(i,1)=0.002e+0_fp
           beta(i,1)=-6.540e+0_fp
        ELSE IF (text(i,1) == 4) THEN
           wp(i,1)=220.000e+0_fp
           fc(i,1)=393.000e+0_fp
           ssc(i,1)=642.000e+0_fp
           alpha(i,1)=0.013e+0_fp
           beta(i,1)=-6.570e+0_fp
        ELSE IF (text(i,1) == 5) THEN
           wp(i,1)=270.000e+0_fp
           fc(i,1)=405.000e+0_fp
           ssc(i,1)=527.000e+0_fp
           alpha(i,1)=0.006e+0_fp
           beta(i,1)=-9.470e+0_fp
        ELSE IF (text(i,1) == 6) THEN
           wp(i,1)=275.000e+0_fp
           fc(i,1)=363.000e+0_fp
           ssc(i,1)=387.000e+0_fp
           alpha(i,1)=0.004e+0_fp
           beta(i,1)=-13.800e+0_fp
        END IF
     ELSE !tree biomes
        IF (text(i,1) == 1) THEN
           wp(i,1)=258.000e+0_fp
           fc(i,1)=463.000e+0_fp
           ssc(i,1)=909.000e+0_fp
           alpha(i,1)=0.002e+0_fp
           beta(i,1)=-6.540e+0_fp
        ELSE IF (text(i,1) == 2) THEN
           wp(i,1)=203.000e+0_fp
           fc(i,1)=506.000e+0_fp
           ssc(i,1)=984.000e+0_fp
           alpha(i,1)=0.002e+0_fp
           beta(i,1)=-5.480e+0_fp
        ELSE IF (text(i,1) == 3) THEN
           wp(i,1)=258.000e+0_fp
           fc(i,1)=463.000e+0_fp
           ssc(i,1)=909.000e+0_fp
           alpha(i,1)=0.002e+0_fp
           beta(i,1)=-6.540e+0_fp
        ELSE IF (text(i,1) == 4) THEN
           wp(i,1)=338.000e+0_fp
           fc(i,1)=604.000e+0_fp
           ssc(i,1)=987.000e+0_fp
           alpha(i,1)=0.013e+0_fp
           beta(i,1)=-6.570e+0_fp
        ELSE IF (text(i,1) == 5) THEN
           wp(i,1)=433.000e+0_fp
           fc(i,1)=647.000e+0_fp
           ssc(i,1)=843.000e+0_fp
           alpha(i,1)=0.006e+0_fp
           beta(i,1)=-9.470e+0_fp
        ELSE IF (text(i,1) == 6) THEN
           wp(i,1)=472.000e+0_fp
           fc(i,1)=622.000e+0_fp
           ssc(i,1)=663.000e+0_fp
           alpha(i,1)=0.004e+0_fp
           beta(i,1)=-13.800e+0_fp
        END IF
     END IF
     SMparams(i,1)=wp(i,1)
     SMparams(i,2)=fc(i,1)
     SMparams(i,3)=ssc(i,1)
     SMparams(i,4)=alpha(i,1)
     SMparams(i,5)=beta(i,1)
     last_soilm(i,1)=wp(i,1) !begin soilmoisture set at wilting point
      
  END DO
!$OMP END PARALLEL DO

END SUBROUTINE getSoilMoistParams
!EOC
