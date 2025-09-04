!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doHerbivory
!
! !DESCRIPTION: Subroutine doHerbivory calculate herbivory analog to 
!  McNaughton (Science, 1989) as fraction of foliage NPP.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doHerbivory
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays

  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
  
  IMPLICIT NONE
!
! !REMARKS:
!  Herbivory analog to McNaughton is computed as:
!                                                                             .
!     log C = 2.04*(log NFP)-4.8   -->  C = NFP^2.04*10^(-4.8)
!                                                                             .
!  where C= consumption, NFP = Net foliage production (NPP
!  delivered to leaves)  units kJ/m2/yr
!
! !REVISION HISTORY:
!  09 Jul 2010 - C. Carouge  - Parallelization.
!  25 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
  INTEGER :: i
  REAL(fp)  :: herb(n_veg, 1)
  character(len=f_len_output+4) :: filename3
  
  filename3(1:f_len_output)=outputpath
      
  !converting kJ/m2/yr to gC/m2/yr
  !
  !NPP(j/m2/yr)=NPP(gC/m2/yr)*energy content/carbon content
  !where energy content = 1.6*10^4
  !and carbon content = 0.45
  !
  !so 1 gC/m2/yr = 35.5 kJ/m2/yr
  !! 1/35.5 = 0.028
  !

  herb(:,1)=0.0e+0_fp
!$OMP PARALLEL DO     &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i)
  DO i=1, n_veg
     herb(i,1)=0.028e+0_fp*(35.5e+0_fp*sum(NPP(i,1:12))/ &
         2e+0_fp)**2.04e+0_fp
     herb(i,1)=herb(i,1)*(10e+0_fp**(-4.8e+0_fp))
     grass_herbivory(i,1)=herb(i,1)
     herb(i,1)=0.0e+0_fp
     herb(i,1)=0.028e+0_fp*(35.5e+0_fp*sum(NPP(i,1:12))/ &
         3e+0_fp)**2.04e+0_fp
     herb(i,1)=herb(i,1)*(10.00e+0_fp**(-4.8e+0_fp))
     trees_herbivory(i,1)=herb(i,1)
  END DO
!$OMP END PARALLEL DO

  !Seasonality in herbivory is calculated as in Randerson et al
  !(GBC 1996) scaling linearly with npp (66%) with a non zero
  !intercept (33%) representing a minimum consumption limit
  !outside the growing season - scalar is equal for C3 and C4
  !NPP
  herb(1,:)=0.0e+0_fp
!$OMP PARALLEL DO     &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i)
  DO i=1, n_veg
     IF (sum(NPP(i,1:12)) .eq. 0e+0_fp) THEN
        herb(i,1)=(0.08333333333e+0_fp)
     ELSE
        herb(i,1)=0.666667e+0_fp*(NPP(i,mo)/sum(NPP(i,1:12)))+ &
             0.33333e+0_fp*0.08333e+0_fp
     END IF
  END DO
!$OMP END PARALLEL DO
  herb_seasonality(:,mo)=herb(:,1)
  
  
  IF (yr .eq. NPPequilibriumYear .and. mo .eq. 12) THEN
     filename3(f_len_output+1:f_len_output+4)='fhsn'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, *) herb_seasonality
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4)='fgrh'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, *) grass_herbivory
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4)='ftrh'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, *) trees_herbivory
     CLOSE(4)
  ENDIF
END SUBROUTINE doHerbivory
!EOC
