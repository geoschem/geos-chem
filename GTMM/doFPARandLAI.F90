!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doFPARandLAI
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doFPARandLAI(FIRST)
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays

  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  implicit none
!
! !INPUT PARAMETERS
!
  LOGICAL, INTENT(IN) :: FIRST
!
! !REVISION HISTORY:
!  09 Jul 2010 - C. Carouge  - Adapted to restart simulations. Parallelization
!  25 Nov 2014 - M. Yannetti - Added PRECISION_MOD 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!  
  INTEGER :: i 
  character(len=f_len_output+4) :: filename3
  
  filename3(1:f_len_output)=outputpath

!--- Previous to (ccc, 11/4/09)      
!      IF (yr .eq. 1 .and. mo .eq. 1) THEN
  IF ( FIRST ) THEN
     maxallLAI=0.0e+0_fp
!$OMP PARALLEL DO      &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i)
     DO i=1,n_veg
        IF (veg1(i,1) .eq. 1) THEN
           srmax(i,1)=SRMAX1
           LAImax(i,1)=LAIMAX1
        ENDIF
        IF (veg1(i,1) .eq. 2) THEN
           srmax(i,1)=SRMAX2
           LAImax(i,1)=LAIMAX2
        ENDIF
        IF (veg1(i,1) .eq. 3) THEN
           srmax(i,1)=SRMAX3
           LAImax(i,1)=LAIMAX3
        ENDIF
        IF (veg1(i,1) .eq. 4) THEN
           srmax(i,1)=SRMAX4
           LAImax(i,1)=LAIMAX4
        ENDIF
        IF (veg1(i,1) .eq. 5) THEN
           srmax(i,1)=SRMAX5
           LAImax(i,1)=LAIMAX5
        ENDIF
        IF (veg1(i,1) .eq. 6) THEN
           srmax(i,1)=SRMAX6
           LAImax(i,1)=LAIMAX6
        ENDIF
        IF (veg1(i,1) .eq. 7) THEN
           srmax(i,1)=SRMAX7
           LAImax(i,1)=LAIMAX7
        ENDIF
        IF (veg1(i,1) .eq. 8) THEN
           srmax(i,1)=SRMAX8
           LAImax(i,1)=LAIMAX8
        ENDIF
        IF (veg1(i,1) .eq. 9) THEN
           srmax(i,1)=SRMAX9
           LAImax(i,1)=LAIMAX9
        ENDIF
        IF (veg1(i,1) .eq. 10) THEN
           srmax(i,1)=SRMAX10
           LAImax(i,1)=LAIMAX10
        ENDIF
        IF (veg1(i,1) .eq. 11) THEN
           srmax(i,1)=SRMAX11
           LAImax(i,1)=LAIMAX11
        ENDIF
        IF (veg1(i,1) .eq. 12) THEN
           srmax(i,1)=SRMAX12
           LAImax(i,1)=LAIMAX12
        ENDIF
     END DO
!$OMP END PARALLEL DO
  END IF
      
  !begin to calculate FPAR
      
!$OMP PARALLEL         &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
  sr(:,1)=(1.000e+0_fp+NDVI1(:,mo))/(1.000e+0_fp-NDVI1(:,mo))


  FPAR(:,mo)=(((sr(:,1)-SRMIN)/(srmax(:,1)-SRMIN))*(FPARMAX-FPARMIN))+  &
             FPARMIN
!$OMP END WORKSHARE

!!$OMP PARALLEL         &
!!$OMP DEFAULT(SHARED)
!$OMP DO PRIVATE(i)
  DO i=1,n_veg
     IF (FPAR(i,mo) .gt. FPARMAX) THEN
        FPAR(i,mo)=FPARMAX
     ENDIF
     IF (FPAR(i,mo) .lt. FPARMIN) THEN
        FPAR(i,mo)=FPARMIN
     ENDIF
  END DO
!$OMP END DO

  ! Begin to calculate LAI
  !clustered canopies (veg = 4, 5, 9)

!$OMP DO PRIVATE(i)
  DO i=1, n_veg
     IF (veg1(i,1) .eq. 4 .or. veg1(i,1) .eq. 5 .or. veg1(i,1) .eq.9) THEN
        LAI_temp(i,1)=LAImax(i,1)*(FPAR(i,mo)/FPARMAX)
     END IF
         
     !mix of clustered and homogeneous canopy
     
     IF (veg1(i,1) .eq. 3) THEN
        LAI_temp(i,1)=(LAImax(i,1)*(log(1e+0_fp-FPAR(i,mo))/ &
                      log(1e+0_fp-FPARMAX)+ &
                      FPAR(i,mo)/FPARMAX))/2e+0_fp
     ENDIF

     !homogeneous canopy
     
     IF (veg1(i,1) .ne. 3 .and.      &
         veg1(i,1) .ne. 4 .and.      &
         veg1(i,1) .ne. 5 .and.      &
         veg1(i,1) .ne. 9) THEN
        LAI_temp(i,1)=LAImax(i,1)*(log(1e+0_fp-FPAR(i,mo))/ &
                      log(1e+0_fp-FPARMAX))
     ENDIF
         
  END DO
!$OMP END DO

!$OMP DO PRIVATE(i)      
  DO i=1, n_veg
     IF (LAI_temp(i,1) .lt. 0) THEN
        LAI_temp(i,1)=0.0e+0_fp
     ENDIF
  END DO
!$OMP END DO
!$OMP END PARALLEL

  maxallLAI=8.0e+0_fp
  LAI(:,mo)=LAI_temp(:,1)
  filename3(f_len_output+1:f_len_output+4)='flai'
  IF (mo .eq. 12 .and. yr .eq. NPPequilibriumYear ) THEN
     OPEN(UNIT=4, file=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4,FMT="(12F9.2)") LAI
     CLOSE(4)
     
  ENDIF
  
END SUBROUTINE doFPARandLAI
!EOC
                 
              
