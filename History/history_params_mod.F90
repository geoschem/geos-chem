!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: History_Params_Mod 
!
! !DESCRIPTION: Contains defined parameters for the GEOS-Chem
!  History Component.
!\\
!\\
! !INTERFACE:
!
MODULE History_Params_Mod
!
! !USES:
!
  USE Precision_Mod
  PRIVATE
!
! !DEFINED PARAMETERS:
!
  ! Missing data value for netCDF
  REAL(f4), PARAMETER, PUBLIC :: UNDEFINED = -9.999e-30_f4
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
END MODULE History_Params_Mod
!EOC
