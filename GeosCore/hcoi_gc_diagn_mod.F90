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
 
  ! GEOS-Chem diagnostic switches and arrays
  USE CMN_SIZE_MOD
  USE CMN_DIAG_MOD
  USE DIAG_MOD

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
  SUBROUTINE HCOI_GC_Diagn_Init( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,   ONLY : HCO_GetHcoID, HCO_State
    USE HCOX_State_Mod,  ONLY : Ext_State
    USE HCO_ExtList_Mod, ONLY : GetExtNr
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
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
    INTEGER                         :: I, ID1, ID2, N, AS
    INTEGER                         :: ExtNr, Cat, Hier
    CHARACTER(LEN=255)              :: MSG, LOC
    CHARACTER(LEN=1)                :: ISTR
    CHARACTER(LEN=4)                :: SpcName
    CHARACTER(LEN=31)               :: DiagnName 

    !=================================================================
    ! HCOI_GC_DIAGN_INIT begins here!
    !=================================================================

    ! Init 
    LOC = 'HCOI_GC_DIAGN_INIT (hcoi_gc_diagn_mod.F90)'
    RC  = HCO_SUCCESS

       !=================================================================
       ! Define GEOS-Chem diagnostics (ADXX)
       !=================================================================

       !-----------------------------------------------------------------
       ! DUST (AD06)
       !-----------------------------------------------------------------
       IF ( ( ExtState%DustDead .OR. ExtState%DustGinoux ) .AND. &
            ND06 > 0 ) THEN

          ! Get Ext. Nr of used extension
          ExtNr = GetExtNr( 'DustDead' )
          IF ( ExtNr <= 0 ) ExtNr = GetExtNr( 'DustGinoux' )

          ! Do for each dust bin
          DO I = 1, NDSTBIN

             ! Get species name (e.g. DST1) and define diagnostics name
             ! (e.g. AD06_DST1).
             WRITE(ISTR,'(i1.1)') I
             SpcName   = 'DST' // TRIM(ISTR)
             DiagnName = 'AD06_' // TRIM(SpcName)         

             ! Get HEMCO species ID 
             ID1 = HCO_GetHcoID( TRIM(SpcName), HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN      
             ENDIF
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = -1,              &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg',            &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 
          ENDDO !Loop over dust species
       ENDIF

       !-----------------------------------------------------------------
       ! CARBON AEROSOLS (AD07)
       !-----------------------------------------------------------------
       IF ( ND07 > 0 .AND. Input_Opt%LCARB ) THEN

          ! Get HEMCO species names
          SpcName = 'BC'
          ID1 = HCO_GetHcoID( TRIM(SpcName), HcoState )
          IF ( ID1 <= 0 ) THEN
             MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN      
          ENDIF

          SpcName = 'OC'
          ID2 = HCO_GetHcoID( TRIM(SpcName), HcoState )
          IF ( ID2 <= 0 ) THEN
             MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN      
          ENDIF

          ! Anthropogenic carbon
          DiagnName = 'AD07_BC_ANTHRO'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = 0,               &
                              Cat       = 1,               &
                              Hier      = -1,              &
                              HcoID     = ID1,             &
                              SpaceDim  = 2,               &
                              LevIDx    = -1,              &
                              OutUnit   = 'kg',            &
                              WriteFreq = 'Manual',        &
                              AutoFill  = 1,               &
                              cID       = N,               & 
                              RC        = RC                ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

          DiagnName = 'AD07_OC_ANTHRO'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = 0,               &
                              Cat       = 1,               &
                              Hier      = -1,              &
                              HcoID     = ID2,             &
                              SpaceDim  = 2,               &
                              LevIDx    = -1,              &
                              OutUnit   = 'kg',            &
                              WriteFreq = 'Manual',        &
                              AutoFill  = 1,               &
                              cID       = N,               & 
                              RC        = RC                ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

          ! Biofuel
          DiagnName = 'AD07_BC_BIOFUEL'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = 0,               &
                              Cat       = 2,               &
                              Hier      = -1,              &
                              HcoID     = ID1,             &
                              SpaceDim  = 2,               &
                              LevIDx    = -1,              &
                              OutUnit   = 'kg',            &
                              WriteFreq = 'Manual',        &
                              AutoFill  = 1,               &
                              cID       = N,               & 
                              RC        = RC                ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

          DiagnName = 'AD07_OC_BIOFUEL'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = 0,               &
                              Cat       = 2,               &
                              Hier      = -1,              &
                              HcoID     = ID2,             &
                              SpaceDim  = 2,               &
                              LevIDx    = -1,              &
                              OutUnit   = 'kg',            &
                              WriteFreq = 'Manual',        &
                              AutoFill  = 1,               &
                              cID       = N,               & 
                              RC        = RC                ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

          ! Biomass burning
          IF ( ExtState%GFED3 .OR. ExtState%FINN ) THEN
         
             ! Get Ext. Nr of used extension
             ExtNr = GetExtNr( 'GFED3' )
             IF ( ExtNr <= 0 ) ExtNr = GetExtNr( 'FINN' )
             IF ( ExtNr <= 0 ) THEN
                MSG = 'Cannot find GFED3 or FINN extension!'
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN      
             ENDIF
             Cat   = -1
          ELSE
             ExtNr = 0
             Cat   = 3
          ENDIF 

          DiagnName = 'AD07_BC_BIOMASS' 
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
                              Cat       = Cat,             &
                              Hier      = -1,              &
                              HcoID     = ID1,             &
                              SpaceDim  = 2,               &
                              LevIDx    = -1,              &
                              OutUnit   = 'kg',            &
                              WriteFreq = 'Manual',        &
                              AutoFill  = 1,               &
                              cID       = N,               & 
                              RC        = RC                ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

          DiagnName = 'AD07_OC_BIOMASS'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
                              Cat       = Cat,             &
                              Hier      = -1,              &
                              HcoID     = ID2,             &
                              SpaceDim  = 2,               &
                              LevIDx    = -1,              &
                              OutUnit   = 'kg',            &
                              WriteFreq = 'Manual',        &
                              AutoFill  = 1,               &
                              cID       = N,               & 
                              RC        = RC                ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

          ! BIOGENIC
          ExtNr = GetExtNr('MEGAN_Mono')
          IF ( ExtNr > 0 ) THEN
             DiagnName = 'AD07_OC_BIOGENIC'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = -1,              &
                                 Hier      = -1,              &
                                 HcoID     = ID2,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg',            &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 
          ENDIF

          ! SOA: NVOC BIOGENIC
          IF ( Input_Opt%LSOA ) THEN
             ExtNr = GetExtNr('MEGAN_SOA')
             IF ( ExtNr < 0 ) THEN
                MSG = 'MEGAN SOA emissions are not turned on!'
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN      
             ENDIF
            
             ! MTPA 
             ID1 = HCO_GetHcoID( 'MTPA', HcoState )
             IF ( ID1 <= 0 ) THEN
                CALL HCO_ERROR ( 'MTPA is not a species', RC, THISLOC=LOC )
                RETURN      
             ENDIF
             CALL Diagn_Create ( am_I_Root,                        & 
                                 HcoState,                         &
                                 cName     = 'AD07_MTPA_BIOGENIC', &
                                 ExtNr     = ExtNr,                &
                                 Cat       = -1,                   &
                                 Hier      = -1,                   &
                                 HcoID     = ID1,                  &
                                 SpaceDim  = 2,                    &
                                 LevIDx    = -1,                   &
                                 OutUnit   = 'kg',                 &
                                 WriteFreq = 'Manual',             &
                                 AutoFill  = 1,                    &
                                 cID       = N,                    & 
                                 RC        = RC                     ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 

             ! MTPO 
             ID1 = HCO_GetHcoID( 'MTPO', HcoState )
             IF ( ID1 <= 0 ) THEN
                CALL HCO_ERROR ( 'MTPO is not a species', RC, THISLOC=LOC )
                RETURN      
             ENDIF
             CALL Diagn_Create ( am_I_Root,                        & 
                                 HcoState,                         &
                                 cName     = 'AD07_MTPO_BIOGENIC', &
                                 ExtNr     = ExtNr,                &
                                 Cat       = -1,                   &
                                 Hier      = -1,                   &
                                 HcoID     = ID1,                  &
                                 SpaceDim  = 2,                    &
                                 LevIDx    = -1,                   &
                                 OutUnit   = 'kg',                 &
                                 WriteFreq = 'Manual',             &
                                 AutoFill  = 1,                    &
                                 cID       = N,                    & 
                                 RC        = RC                     ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 

             ! LIMO 
             ID1 = HCO_GetHcoID( 'LIMO', HcoState )
             IF ( ID1 <= 0 ) THEN
                CALL HCO_ERROR ( 'LIMO is not a species', RC, THISLOC=LOC )
                RETURN      
             ENDIF
             CALL Diagn_Create ( am_I_Root,                        & 
                                 HcoState,                         &
                                 cName     = 'AD07_LIMO_BIOGENIC', &
                                 ExtNr     = ExtNr,                &
                                 Cat       = -1,                   &
                                 Hier      = -1,                   &
                                 HcoID     = ID1,                  &
                                 SpaceDim  = 2,                    &
                                 LevIDx    = -1,                   &
                                 OutUnit   = 'kg',                 &
                                 WriteFreq = 'Manual',             &
                                 AutoFill  = 1,                    &
                                 cID       = N,                    & 
                                 RC        = RC                     ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 

             ! SESQ 
             ID1 = HCO_GetHcoID( 'SESQ', HcoState )
             IF ( ID1 <= 0 ) THEN
                CALL HCO_ERROR ( 'SESQ is not a species', RC, THISLOC=LOC )
                RETURN      
             ENDIF
             CALL Diagn_Create ( am_I_Root,                        & 
                                 HcoState,                         &
                                 cName     = 'AD07_SESQ_BIOGENIC', &
                                 ExtNr     = ExtNr,                &
                                 Cat       = -1,                   &
                                 Hier      = -1,                   &
                                 HcoID     = ID1,                  &
                                 SpaceDim  = 2,                    &
                                 LevIDx    = -1,                   &
                                 OutUnit   = 'kg',                 &
                                 WriteFreq = 'Manual',             &
                                 AutoFill  = 1,                    &
                                 cID       = N,                    & 
                                 RC        = RC                     ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 
          ENDIF
       ENDIF ! CARBON


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
