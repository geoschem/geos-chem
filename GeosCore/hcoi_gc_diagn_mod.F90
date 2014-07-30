!BOM
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_gc_diagn_mod.F90
!
! !DESCRIPTION: Module HCOi\_GC\_Diagn\_Mod.F90 is the GEOS-Chem interface 
! module for the HEMCO diagnostics.
! \\
! !INTERFACE:
!
MODULE HCOI_GC_DIAGN_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_DIAGN_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOI_GC_DIAGN_INIT
!
! !REMARKS:
!  HEMCO diagnostics are still in testing mode. We will fully activate them
!  at a later time.  They will be turned on when debugging & unit testing.
!
! !REVISION HISTORY:
!  04 May 2014 - C. Keller   - Initial version. 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  28 Jul 2014 - C. Keller   - Split off from hcoio_diagn_mod.F90 and moved
!                              from HEMCO/Interface to GeosCore.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_DiagN_Init
!
! !DESCRIPTION: Subroutine HCOI\_GC\_Diagn\_Init initializes the HEMCO 
! diagnostics in GEOS-Chem. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Diagn_Init( am_I_Root, HcoState, RC ) 
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_GetHcoID, HCO_State
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                         :: I, N, AS
    CHARACTER(LEN=255)              :: MSG, LOC

    !=================================================================
    ! HCOI_GC_DIAGN_INIT begins here!
    !=================================================================

    ! Init 
    LOC = 'HCOI_GC_DIAGN_INIT (hcoi_gc_diagn_mod.F90)'
    RC  = HCO_SUCCESS

       !=================================================================
       ! Define manual diagnostics (no AutoFill)
       !=================================================================

       !=================================================================
       ! Define automatic diagnostics (AutoFill)
       !=================================================================

       ! Total NO
       I = HCO_GetHcoID( 'NO', HcoState )
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'NOtotal', &
                           ExtNr    = -1, &
                           Cat      = -1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           LevIDx   = -1, &
                           OutUnit  = 'kg/m2/s', &
                           WriteFreq = 'Hourly',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from anthropogenic emissions 
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'NO_anthr', &
                           ExtNr    = 0,  &
                           Cat      = 1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           OutUnit  = 'kg/m2/s', &
                           WriteFreq = 'Hourly',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from aircrafts (AEIC) 
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName     = 'NO_aircraft', &
                           ExtNr     = 0,  &
                           Cat       = 20, &
                           Hier      = -1, &
                           HcoID     = I, &
                           SpaceDim  = 2, &
                           LevIdx    = -1, &
                           OutUnit   = 'kg/m2/s', &
                           WriteFreq = 'Hourly',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from biomass burning 
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'NO_biomass', &
                           ExtNr    = 111, &
                           Cat      = -1,  &
                           Hier     = -1,  &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           OutUnit  = 'kg/m2/s', &
                           WriteFreq = 'Hourly',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from ships (ParaNox) 
       I = HCO_GetHcoID( 'NO', HcoState )
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'ParaNOx', &
                           ExtNr    = 102,  &
                           Cat      = -1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           LevIDx   = -1, &
                           OutUnit  = 'kg/m2/s', &
                           WriteFreq = 'Hourly',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from soils 
       I = HCO_GetHcoID( 'NO', HcoState )
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'SoilNOx', &
                           ExtNr    = 104,  &
                           Cat      = -1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           OutUnit  = 'kg/m2/s', &
                           WriteFreq = 'Hourly',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from lightning 
       I = HCO_GetHcoID( 'NO', HcoState )
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'LightNOx', &
                           ExtNr    = 103,  &
                           Cat      = -1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           LevIDx   = -1, &
                           OutUnit  = 'kg/m2/s', &
                           WriteFreq = 'Hourly',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

!       ! Total SO2
!       I = HCO_GetHcoID( 'SO2', HcoState )
!       CALL Diagn_Create ( am_I_Root, &
!                           HcoState,  &
!                           cName    = 'SO2total', &
!                           ExtNr    = -1, &
!                           Cat      = -1, &
!                           Hier     = -1, &
!                           HcoID    = I, &
!                           SpaceDim = 2, &
!                           LevIDx   = -1, &
!                           OutUnit  = 'kg/m2/s', &
!                           WriteFreq = 'Hourly',  &
!                           AutoFill  = 1, &
!                           cID       = N, & 
!                           RC        = RC ) 
!       IF ( RC /= HCO_SUCCESS ) RETURN 

    ! Leave w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCOI_GC_Diagn_Init
!EOC
END MODULE HCOI_GC_Diagn_Mod
