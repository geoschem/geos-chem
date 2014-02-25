!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doLeafRootShedding
!
! !DESCRIPTION: \subsection*{Overview}
!  Subroutine doLeafRootShedding define the scalars that predict 
!  the seasonality of leaf shedding and root decay based on changes in LAI
!  This needs improvement.
!
!\subsection*{References}
! Randerson, Thompson, Malmstrom, Field and Fung 1996, GBC 10(4) p585
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doLeafRootShedding
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
  
  implicit none
!
! !REVISION HISTORY:
!  09 July 2010 - C. Carouge  - Parallelization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER :: i,j
  character(len=f_len_output+4) :: filename3

  filename3(1:f_len_output)=outputpath
  
  
  MINLAI(:,1)=12.000d0

  SUMLAI(:,1)=sum(lais(:,2:13),2)
  
  SUMLAInew(:,mo)=SUMLAI(:,1)
  LTVARSUM(:,1)=0.0d0

  DO i=13,2,-1
!$OMP PARALLEL DO        &
!$OMP DEFAULT(SHARED)    &
!$OMP PRIVATE(j)
     DO j=1, n_veg
        IF (lais(j,i) .gt. lais(j,i-1)) THEN
           LTVARSUM(j,1)=LTVARSUM(j,1)+(lais(j,i)-lais(j,i-1))
        ENDIF
        IF (lais(j,i) .gt. MAXLAI(j,1)) THEN
           MAXLAI(j,1)=lais(j,i)
        ENDIF
        IF (lais(j,i) .lt. MINLAI(j,1)) THEN
           MINLAI(j,1)=lais(j,i)
        ENDIF
     END DO
!$OMP END PARALLEL DO
  END DO
      
  AVELAI(:,1)=SUMLAI(:,1)/12.00000d0

!$OMP PARALLEL DO        &
!$OMP DEFAULT(SHARED)    &
!$OMP PRIVATE(i)
  DO i=1, n_veg
     IF (AVELAI(i,1) .gt. 0d0) THEN
        LTCON(i,1)=MINLAI(i,1)/AVELAI(i,1)
     ELSE
        LTCON(i,1)=0.0d0
     ENDIF

     IF (lais(i,2)-lais(i,1) .gt. 0d0) THEN
        LTVAR(i,1)=(lais(i,2)-lais(i,1))
     ELSE
        LTVAR(i,1)=0.000d0
     ENDIF

     IF (LTVARSUM(i,1) .gt. 0d0) THEN
        litterscalar(i,mo)=(LTCON(i,1)/12.000d0)+(1.000d0-LTCON(i,1))*(LTVAR(i,1)/LTVARSUM(i,1))
     ELSE
        litterscalar(i,mo)=0.000d0
     ENDIF

     IF(SUMLAI(i,1) .gt. 0d0) THEN
        rootlitscalar(i,mo)=0.7d0*((litterscalar(i,mo)+LAI(i,mo)/SUMLAI(i,1))/2.000d0)+(0.3d0/12.000d0)
     ELSE
        rootlitscalar(i,mo)=0.000d0
     ENDIF

     hlitterscalar(i,mo)=0.500d0+(0.500d0-abiotic(i,mo)/2.00d0)
  END DO
!$OMP END PARALLEL DO
  
  IF (yr .eq. NPPequilibriumYear .and. mo .eq. 12) THEN
     filename3(f_len_output+1:f_len_output+4)='fltc'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, *) LTCON
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4)='flvr'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, *) LTVARSUM
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4)='fl_i'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, FMT="(12F9.2)") LAI
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4)='frls'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, FMT="(12F9.2)") rootlitscalar
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4)='f_ls'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, FMT="(12F9.2)") litterscalar
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4)='fhls'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, FMT="(12F9.2)") hlitterscalar
     CLOSE(4)
     
     filename3(f_len_output+1:f_len_output+4)='fali'
     OPEN(UNIT=4, FILE=filename3, STATUS="NEW", &
          FORM="FORMATTED")
     WRITE(4, *) AVELAI
     CLOSE(4)
  ENDIF
  
END SUBROUTINE doLeafRootShedding
!EOC
