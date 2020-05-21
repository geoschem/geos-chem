!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_vertgrid_mod.F90
!
! !DESCRIPTION: Module HCO\_VERTGRID\_Mod contains routines and
! variables for the definition of the vertical grid and related
! calculations.
! \\
! !INTERFACE:
!
MODULE HCO_VertGrid_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Arr_Mod
  USE HCO_Types_Mod, ONLY : VertGrid

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCO_VertGrid_Init
  PUBLIC :: HCO_VertGrid_Define
  PUBLIC :: HCO_VertGrid_Cleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !PARAMETERS:
!
  INTEGER,  PARAMETER, PUBLIC  :: HCO_ZTYPE_HYBSIG = 1

  ! Ap [Pa] for 47 levels (48 edges)
  REAL(hp), PARAMETER          :: Ap47(48) = (/                       &
              0.000000e+00_hp, 4.804826e-00_hp, 6.593752e+02_hp, 1.313480e+03_hp, &
              1.961311e+03_hp, 2.609201e+03_hp, 3.257081e+03_hp, 3.898201e+03_hp, &
              4.533901e+03_hp, 5.169611e+03_hp, 5.805321e+03_hp, 6.436264e+03_hp, &
              7.062198e+03_hp, 7.883422e+03_hp, 8.909992e+03_hp, 9.936521e+03_hp, &
              1.091817e+04_hp, 1.189586e+04_hp, 1.286959e+04_hp, 1.429100e+04_hp, &
              1.562600e+04_hp, 1.696090e+04_hp, 1.816190e+04_hp, 1.930970e+04_hp, &
              2.032590e+04_hp, 2.121500e+04_hp, 2.187760e+04_hp, 2.238980e+04_hp, &
              2.243630e+04_hp, 2.168650e+04_hp, 2.011920e+04_hp, 1.769300e+04_hp, &
              1.503930e+04_hp, 1.278370e+04_hp, 1.086630e+04_hp, 9.236572e+03_hp, &
              7.851231e+03_hp, 5.638791e+03_hp, 4.017541e+03_hp, 2.836781e+03_hp, &
              1.979160e+03_hp, 9.292942e+02_hp, 4.076571e+02_hp, 1.650790e+02_hp, &
              6.167791e+01_hp, 2.113490e+01_hp, 6.600001e+00_hp, 1.000000e+00_hp  &
                                                                      /)

  ! Bp [unitless] for 47 levels (48 edges)
  REAL(hp), PARAMETER          :: Bp47(48) = (/                       &
              1.000000e+00_hp, 9.849520e-01_hp, 9.634060e-01_hp, 9.418650e-01_hp, &
              9.203870e-01_hp, 8.989080e-01_hp, 8.774290e-01_hp, 8.560180e-01_hp, &
              8.346609e-01_hp, 8.133039e-01_hp, 7.919469e-01_hp, 7.706375e-01_hp, &
              7.493782e-01_hp, 7.211660e-01_hp, 6.858999e-01_hp, 6.506349e-01_hp, &
              6.158184e-01_hp, 5.810415e-01_hp, 5.463042e-01_hp, 4.945902e-01_hp, &
              4.437402e-01_hp, 3.928911e-01_hp, 3.433811e-01_hp, 2.944031e-01_hp, &
              2.467411e-01_hp, 2.003501e-01_hp, 1.562241e-01_hp, 1.136021e-01_hp, &
              6.372006e-02_hp, 2.801004e-02_hp, 6.960025e-03_hp, 8.175413e-09_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp  &
                                                                      /)

  ! Ap [Pa] for 72 levels (73 edges)
  REAL(hp), PARAMETER          :: Ap72(73) = (/                       &
              0.000000e+00_hp, 4.804826e+00_hp, 6.593752e+02_hp, 1.313480e+03_hp, &
              1.961311e+03_hp, 2.609201e+03_hp, 3.257081e+03_hp, 3.898201e+03_hp, &
              4.533901e+03_hp, 5.169611e+03_hp, 5.805321e+03_hp, 6.436264e+03_hp, &
              7.062198e+03_hp, 7.883422e+03_hp, 8.909992e+03_hp, 9.936521e+03_hp, &
              1.091817e+04_hp, 1.189586e+04_hp, 1.286959e+04_hp, 1.429100e+04_hp, &
              1.562600e+04_hp, 1.696090e+04_hp, 1.816190e+04_hp, 1.930970e+04_hp, &
              2.032590e+04_hp, 2.121500e+04_hp, 2.187760e+04_hp, 2.238980e+04_hp, &
              2.243630e+04_hp, 2.168650e+04_hp, 2.011920e+04_hp, 1.769300e+04_hp, &
              1.503930e+04_hp, 1.278370e+04_hp, 1.086630e+04_hp, 9.236572e+03_hp, &
              7.851231e+03_hp, 6.660341e+03_hp, 5.638791e+03_hp, 4.764391e+03_hp, &
              4.017541e+03_hp, 3.381001e+03_hp, 2.836781e+03_hp, 2.373041e+03_hp, &
              1.979160e+03_hp, 1.645710e+03_hp, 1.364340e+03_hp, 1.127690e+03_hp, &
              9.292942e+02_hp, 7.619842e+02_hp, 6.216801e+02_hp, 5.046801e+02_hp, &
              4.076571e+02_hp, 3.276431e+02_hp, 2.620211e+02_hp, 2.084970e+02_hp, &
              1.650790e+02_hp, 1.300510e+02_hp, 1.019440e+02_hp, 7.951341e+01_hp, &
              6.167791e+01_hp, 4.758061e+01_hp, 3.650411e+01_hp, 2.785261e+01_hp, &
              2.113490e+01_hp, 1.594950e+01_hp, 1.197030e+01_hp, 8.934502e+00_hp, &
              6.600001e+00_hp, 4.758501e+00_hp, 3.270000e+00_hp, 2.000000e+00_hp, &
              1.000000e+00_hp /)

  ! Bp [unitless] for 72 levels (73 edges)
  REAL(hp), PARAMETER          :: Bp72(73) = (/                       &
              1.000000e+00_hp, 9.849520e-01_hp, 9.634060e-01_hp, 9.418650e-01_hp, &
              9.203870e-01_hp, 8.989080e-01_hp, 8.774290e-01_hp, 8.560180e-01_hp, &
              8.346609e-01_hp, 8.133039e-01_hp, 7.919469e-01_hp, 7.706375e-01_hp, &
              7.493782e-01_hp, 7.211660e-01_hp, 6.858999e-01_hp, 6.506349e-01_hp, &
              6.158184e-01_hp, 5.810415e-01_hp, 5.463042e-01_hp, 4.945902e-01_hp, &
              4.437402e-01_hp, 3.928911e-01_hp, 3.433811e-01_hp, 2.944031e-01_hp, &
              2.467411e-01_hp, 2.003501e-01_hp, 1.562241e-01_hp, 1.136021e-01_hp, &
              6.372006e-02_hp, 2.801004e-02_hp, 6.960025e-03_hp, 8.175413e-09_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, 0.000000e+00_hp, &
              0.000000e+00_hp /)
!
! PUBLIC TYPES:
!
!
! !REVISION HISTORY:
!  28 Sep 2015 - C. Keller   - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_VertGrid_Init
!
! !DESCRIPTION: Function HCO\_VertGrid\_Init initializes the vertical
!  grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_VertGrid_Init( zGrid, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(VertGrid), POINTER        :: zGrid     ! vertical grid
    INTEGER,        INTENT(INOUT)  :: RC        ! Return code
!
! !REVISION HISTORY:
!  28 Sep 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    !=====================================================================
    ! HCO_VertGrid_Init begins here!
    !=====================================================================

    ! Initialize
    IF ( .NOT. ASSOCIATED(zGrid) ) ALLOCATE( zGrid )

    ! Initialize vertical grid type. For now, always assume to be hybrid
    ! sigma coordinates
    zGrid%ZTYPE =  HCO_ZTYPE_HYBSIG
    zGrid%Ap    => NULL()
    zGrid%Bp    => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_VertGrid_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_VertGrid_Define
!
! !DESCRIPTION: Function HCO\_VertGrid\_Define initializes the vertical
!  grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_VertGrid_Define( HcoConfig, zGrid, nz, Ap, Bp, RC )
!
! !USES:
!
    USE HCO_TYPES_MOD,    ONLY : ConfigObj
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN   )            :: nz        ! # of vertical levels
    REAL(hp),       INTENT(IN   ), OPTIONAL  :: Ap(nz+1)  ! Ap values
    REAL(hp),       INTENT(IN   ), OPTIONAL  :: Bp(nz+1)  ! Bp values
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ConfigObj),POINTER                  :: HcoConfig ! HEMCO config obj
    TYPE(VertGrid), POINTER                  :: zGrid     ! vertical grid
    INTEGER,        INTENT(INOUT)            :: RC        ! Return code
!
! !REVISION HISTORY:
!  28 Sep 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, AS
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = 'HCO_VertGrid_Define (hco_vertgrid_mod.F90)'

    !=====================================================================
    ! HCO_VertGrid_Define begins here!
    !=====================================================================

    ! Allocate AP and BP
    ALLOCATE(zGrid%Ap(nz+1), zGrid%Bp(nz+1), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoConfig%Err, 'Cannot allocate Ap / Bp', RC, THISLOC=LOC )
       RETURN
    ENDIF
    zGrid%Ap = 0.0_hp
    zGrid%Bp = 0.0_hp

    ! Set Ap
    IF ( PRESENT(Ap) ) THEN
       zGrid%Ap(:) = Ap(:)
    ELSE
       IF ( nz > 72 ) THEN
          WRITE(MSG,*) 'Vertical grid has more than 72 vertical levels', &
                       '- please provide Ap values in configuration file.'
          CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ELSEIF ( nz > 47 ) THEN
          zGrid%Ap(:) = Ap72(1:(nz+1))
       ELSE
          zGrid%Ap(:) = Ap47(1:(nz+1))
       ENDIF
    ENDIF

    ! Set Bp
    IF ( PRESENT(Bp) ) THEN
       zGrid%Bp(:) = Bp(:)
    ELSE
       IF ( nz > 72 ) THEN
          WRITE(MSG,*) 'Vertical grid has more than 72 vertical levels', &
                       '- please provide Bp values in configuration file.'
          CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ELSEIF ( nz > 47 ) THEN
          zGrid%Bp(:) = Bp72(1:(nz+1))
       ELSE
          zGrid%Bp(:) = Bp47(1:(nz+1))
       ENDIF
    ENDIF

    ! Verbose
    IF ( HcoConfig%amIRoot .AND. HCO_IsVerb(HcoConfig%Err,1) ) THEN
       WRITE(MSG,*) ' HEMCO vertical sigma-hybrid coordinates: '
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,*) 'Ap [Pa]       (first and last): ', zGrid%Ap(1), zGrid%Ap(nz+1)
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,*) 'Bp [unitless] (first and last): ', zGrid%Bp(1), zGrid%Bp(nz+1)
       CALL HCO_MSG(HcoConfig%Err,MSG,SEP2='-')
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_VertGrid_Define
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_VertGrid_Cleanup
!
! !DESCRIPTION: Function HCO\_VertGrid\_Cleanup cleans up the vertical
!  grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_VertGrid_Cleanup( zGrid )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(VertGrid), POINTER        :: zGrid     ! vertical grid
!
! !REVISION HISTORY:
!  28 Sep 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Eventually deallocate Ap/Bp vectors
    IF( ASSOCIATED(zGrid%Ap) ) THEN
       DEALLOCATE(zGrid%Ap)
       zGrid%Ap => NULL()
    ENDIF
    IF( ASSOCIATED(zGrid%Bp) ) THEN
       DEALLOCATE(zGrid%Bp)
       zGrid%Bp => NULL()
    ENDIF

  END SUBROUTINE HCO_VertGrid_Cleanup
!EOC
END MODULE HCO_VertGrid_Mod
!EOM
