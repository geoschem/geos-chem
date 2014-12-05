!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: load_GC_data
!
! !DESCRIPTION: Subroutine load\_GC\_data is only used with GTMM coupled to
!  GEOS-Chem. The subroutine regrid the temperature, precipation and
!  radiation fields to 1x1. The met fields are read in GTMM\_DR in GEOS\_Chem.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE load_GC_data(month, TS, PREACC, RADSWG)
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE CasaRegridModule
  
  implicit none
!
! !INPUT PARAMETERS:
!
  INTEGER, INTENT(IN) :: month
  REAL*8, INTENT(IN), DIMENSION(I4x5, J4x5)  :: TS, &
       PREACC, RADSWG
!
! !REVISION HISTORY:
!  09 July 2010 - C. Carouge  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
!EOC
