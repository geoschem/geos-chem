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
  ! Missing data values
  INTEGER,          PARAMETER, PUBLIC :: UNDEFINED_INT = -999
  REAL(f4),         PARAMETER, PUBLIC :: UNDEFINED     = -1.0e+31_f4
  CHARACTER(LEN=9), PARAMETER, PUBLIC :: UNDEFINED_STR = 'not found'

  ! Operation codes:
  ! 0 = Copy       data from source pointer to the HISTORY ITEM data array
  ! 1 = Accumulate data from source pointer to the HISTORY ITEM data array
  INTEGER,          PARAMETER, PUBLIC :: COPY_FROM_SOURCE  = 0
  INTEGER,          PARAMETER, PUBLIC :: ACCUM_FROM_SOURCE = 1

  ! Action codes (used for testing if it is time to do something)
  ! 0 = Update (aka archive) data from the source into the HISTORY ITEM
  ! 1 = Write data to the netCDF file
  ! 2 = Close the current netCDF file and open the one for the next interval
  INTEGER,          PARAMETER, PUBLIC :: ACTION_UPDATE     = 0
  INTEGER,          PARAMETER, PUBLIC :: ACTION_FILE_WRITE = 1
  INTEGER,          PARAMETER, PUBLIC :: ACTION_FILE_CLOSE = 2
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!  03 Aug 2017 - R. Yantosca - Add operation code parameters
!  04 Aug 2017 - R. Yantosca - Rename operation code parameters
!  09 Aug 2017 - R. Yantosca - Add UNDEFINED_INT
!  15 Aug 2017 - R. Yantosca - Add UNDEFINED_STR
!  16 Aug 2017 - R. Yantosca - Added ACTION_* parameters
!EOP
!------------------------------------------------------------------------------
!BOC
END MODULE History_Params_Mod
!EOC
