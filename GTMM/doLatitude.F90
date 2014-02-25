!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doLatitude
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doLatitude
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER :: i, j
  REAL*8  :: startLatitude 
  
  startLatitude=90.000d0-(180.000d0/rows)/2.000d0
  
  DO i=1,rows
     latitude(:,i)=startLatitude
     startLatitude=startLatitude-180.000d0/rows
  END DO

  !all latitudes north of 50N become 50, all south of 50S
  !become -50
  
  DO i=1,rows
     IF (latitude(1,i) .gt. 50d0) THEN 
        latitude(1:columns, i)=50d0
     ELSE IF (latitude(i,1) .lt. -50d0) THEN
        latitude(1:columns, i)=-50d0
     ENDIF
  END DO

  latitude1=maskfile(latitude, mask2)
  
END SUBROUTINE doLatitude
!EOC
