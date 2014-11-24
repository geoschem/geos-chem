!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_grid_mod.F90
!
! !DESCRIPTION: Module HCO\_GRID\_MOD contains routines to specify 
!  parameters of the HEMCO emissions grid.
!\\  
!\\
! !INTERFACE: 
!
MODULE HCO_Grid_Mod
! 
! !USES:
!
  USE HCO_ERROR_MOD

  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCO_SET_GRID
!
! !REVISION HISTORY:
!  23 Aug 2013 - C. Keller   - Initial version, Based on grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GRID_SET
!
! !DESCRIPTION: Subroutine HCO\_SET\_GRID defines the HEMCO grid. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_SET_GRID( HcoState, I,     J,     L,         &
                             XMID,      YMID,  XEDGE, YSIN,      &
                             AREA_M2,   BXHEIGHT_M,   OnSimGrid, &
                             RC ) 
!
! !USES:
!
    USE HCO_TYPE_MOD,       ONLY : HCO_State
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState         ! HEMCO State object
!
! !INPUT PARAMETERS:
!  
    INTEGER, INTENT(IN)            :: I                 ! # of lons on this CPU 
    INTEGER, INTENT(IN)            :: J                 ! # of lats on this CPU 
    INTEGER, INTENT(IN)            :: L                 ! # of levs on this CPU 
    REAL*8,  INTENT(IN), TARGET    :: XMID(I,J,L)       ! lon midpoints
    REAL*8,  INTENT(IN), TARGET    :: YMID(I,J,L)       ! lat midpoints
    REAL*8,  INTENT(IN), TARGET    :: XEDGE(I+1,J,L)    ! lon edges
    REAL*8,  INTENT(IN), TARGET    :: YSIN(I,J+1,L)     ! lat edges (sin)
    REAL*8,  INTENT(IN), TARGET    :: AREA_M2(I,J,L)    ! grid box areas (m2)
    REAL*8,  INTENT(IN), TARGET    :: BXHEIGHT_M(I,J,L) ! grid box heights (m)
    LOGICAL, INTENT(IN)            :: OnSimGrid         ! Emission grid = Simulation grid?
!
! !OUTPUT PARAMETERS:
!  
    INTEGER, INTENT(INOUT)         :: RC                ! Success or failure?
!
! !REMARKS:
!  23 Aug 2013 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 

    !======================================================================
    ! HCO_Set_Grid begins here 
    !======================================================================

    ! Copy grid dimensions
    HcoState%ISIZE = I
    HcoState%JSIZE = J
    HcoState%LSIZE = L

    ! Assign pointers to grid arrays
    HcoState%XMID       => XMID 
    HcoState%YMID       => YMID 
    HcoState%XEDGE      => XEDGE 
    HcoState%YSIN       => YSIN
    HcoState%AREA_M2    => AREA_M2 
    HcoState%BXHEIGHT_M => BXHEIGHT_M 

    ! Is emission grid equal to simulation grid?
    HcoState%OnSimGrid = OnSimGrid

    ! Return w/ success
    RC = HCO_SUCCESS

END SUBROUTINE HCO_SET_GRID
!EOC
END MODULE HCO_Grid_Mod
