!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: history_params_mod.F90
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

  ! Operation codes:
  ! 0 = Copy       data from source pointer to the HISTORY ITEM data array
  ! 1 = Accumulate data from source pointer to the HISTORY ITEM data array
  INTEGER,  PARAMETER, PUBLIC :: COPY_FROM_SOURCE  = 0
  INTEGER,  PARAMETER, PUBLIC :: ACCUM_FROM_SOURCE = 1
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!  03 Aug 2017 - R. Yantosca - Add operation code parameters
!  04 Aug 2017 - R. Yantosca - Rename operation code parameters
!EOP
!------------------------------------------------------------------------------
!BOC
END MODULE History_Params_Mod
!EOC
