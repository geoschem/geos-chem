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
  
  IMPLICIT NONE
!
! !REVISION HISTORY:
!  09 July 2010 - C. Carouge  - Parallelization
!EOP  
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
  REAL*8 :: wp(n_veg, 1)
  REAL*8 :: fc(n_veg, 1)
  REAL*8 :: ssc(n_veg, 1)
  REAL*8 :: alpha(n_veg, 1)
  REAL*8 :: beta(n_veg, 1)
  REAL*8 :: text(n_veg, 1)
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
           wp(i,1)=129.000d0
           fc(i,1)=232.000d0
           ssc(i,1)=455.000d0
           alpha(i,1)=0.0020d0
           beta(i,1)=-6.540d0
        ELSE IF (text(i,1) == 2) THEN
           wp(i,1)=81.000d0
           fc(i,1)=202.000d0
           ssc(i,1)=394.000d0
           alpha(i,1)=0.002d0
           beta(i,1)=-5.480d0
        ELSE IF (text(i,1) == 3) THEN
           wp(i,1)=129.000d0
           fc(i,1)=232.000d0
           ssc(i,1)=455.000d0
           alpha(i,1)=0.002d0
           beta(i,1)=-6.540d0
        ELSE IF (text(i,1) == 4) THEN
           wp(i,1)=220.000d0
           fc(i,1)=393.000d0
           ssc(i,1)=642.000d0
           alpha(i,1)=0.013d0
           beta(i,1)=-6.570d0
        ELSE IF (text(i,1) == 5) THEN
           wp(i,1)=270.000d0
           fc(i,1)=405.000d0
           ssc(i,1)=527.000d0
           alpha(i,1)=0.006d0
           beta(i,1)=-9.470d0
        ELSE IF (text(i,1) == 6) THEN
           wp(i,1)=275.000d0
           fc(i,1)=363.000d0
           ssc(i,1)=387.000d0
           alpha(i,1)=0.004d0
           beta(i,1)=-13.800d0
        END IF
     ELSE !tree biomes
        IF (text(i,1) == 1) THEN
           wp(i,1)=258.000d0
           fc(i,1)=463.000d0
           ssc(i,1)=909.000d0
           alpha(i,1)=0.002d0
           beta(i,1)=-6.540d0
        ELSE IF (text(i,1) == 2) THEN
           wp(i,1)=203.000d0
           fc(i,1)=506.000d0
           ssc(i,1)=984.000d0
           alpha(i,1)=0.002d0
           beta(i,1)=-5.480d0
        ELSE IF (text(i,1) == 3) THEN
           wp(i,1)=258.000d0
           fc(i,1)=463.000d0
           ssc(i,1)=909.000d0
           alpha(i,1)=0.002d0
           beta(i,1)=-6.540d0
        ELSE IF (text(i,1) == 4) THEN
           wp(i,1)=338.000d0
           fc(i,1)=604.000d0
           ssc(i,1)=987.000d0
           alpha(i,1)=0.013d0
           beta(i,1)=-6.570d0
        ELSE IF (text(i,1) == 5) THEN
           wp(i,1)=433.000d0
           fc(i,1)=647.000d0
           ssc(i,1)=843.000d0
           alpha(i,1)=0.006d0
           beta(i,1)=-9.470d0
        ELSE IF (text(i,1) == 6) THEN
           wp(i,1)=472.000d0
           fc(i,1)=622.000d0
           ssc(i,1)=663.000d0
           alpha(i,1)=0.004d0
           beta(i,1)=-13.800d0
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
