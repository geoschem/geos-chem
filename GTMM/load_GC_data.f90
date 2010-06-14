SUBROUTINE load_GC_data(month, TS, PREACC, RADSWG)

  USE defineConstants
  USE loadCASAinput
  USE CasaRegridModule
  
  implicit none

  ! ARGUMENTS
  INTEGER, INTENT(IN) :: month
  REAL*8, INTENT(IN), DIMENSION(I4x5, J4x5)  :: TS, &
       PREACC, RADSWG
  real*8, dimension(columns, rows)      ::   geos1x1
  

  CALL regrid4x5to1x1(1, TS, geos1x1)
  airt(:,:,month)=geos1x1
  airt(:,:,month)=airt(:,:,month)-273.15d0
  

  CALL regrid4x5to1x1(1, PREACC, geos1x1)
  ppt(:,:,month)=geos1x1
  ppt_mo(1,month)=sum(ppt(:,:,month))
  

  CALL regrid4x5to1x1(1, RADSWG, geos1x1)
  solrad(:,:,month)=geos1x1
  
END SUBROUTINE load_GC_data
