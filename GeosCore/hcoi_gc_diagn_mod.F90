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
! module for the HEMCO diagnostics. For every GEOS-Chem emissions diagnostics,
! a corresponding HEMCO diagnostics is created. The HEMCO diagnostics become
! (automatically) filled and updated when calling HEMCO. They are passed
! back to GEOS-Chem when writing the diagnostics (e.g. in diag3.F).
!
! Notes:
! \begin{itemize}
! \item The category specific diagnostics (anthropogenic, aircraft, etc.)
!  explicitly assume certain category numbers in the HEMCO configuration 
!  file (e.g. Cat=1 for anthropogenic, Cat=20 for aircraft, etc.).
!  Diagnostics will not represent what they should if these category numbers
!  get changed!
! \item Most biofuel emissions are included in the anthropogenic inventories
!  and hence not distinguishable from those. For most compounds, no biofuel 
!  diagnostics are written.
! \item In HEMCO, ocean sinks are treated as drydep and the calculated 
!  deposition velocities are passed to drydep\_mod.F. Hence, no Acetone
!  ocean sink is calculated by HEMCO and the DMS diagnostics only includes
!  the ocean flux (this is NOT the net flux!!). 
!  If needed, we can build a simple wrapper in hcoi\_gc\_main\_mod.F90 that
!  explicitly calculates oceanic fluxes.
! \end{itemize}
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
          IF ( ExtNr <= 0 ) THEN
             CALL HCO_ERROR ( 'Cannot find dust extension', RC, THISLOC=LOC )
             RETURN      
          ENDIF

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
       ! ==> BIOMASS diagnostics are defined below (ND28)
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

       !-----------------------------------------------------------------
       ! Sea salt aerosols (AD08)
       !-----------------------------------------------------------------
       IF ( ND08 > 0 .AND. Input_Opt%LSSALT .AND. ExtState%SeaSalt ) THEN

          ! Get Ext. Nr of used extension
          ExtNr = GetExtNr( 'SeaSalt' )
          IF ( ExtNr <= 0 ) THEN
             CALL HCO_ERROR ( 'Cannot find extension SeaSalt', RC, THISLOC=LOC )
             RETURN      
          ENDIF

          ! Do for both sea salt aerosol modes
          DO I = 1, 2

             ! Get species name and define diagnostics name
             IF ( I == 1 ) THEN
                SpcName   = 'SALA'
                DiagnName = 'AD08_SALA'
             ELSEIF ( I == 2 ) THEN
                SpcName   = 'SALC'
                DiagnName = 'AD08_SALC'
             ENDIF

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
          ENDDO !Loop over sea salt species
       ENDIF

       !-----------------------------------------------------------------
       ! ACETONE (AD11)
       !-----------------------------------------------------------------
       ! --> 3 manually defined diagnostics in MEGAN
       ! --> 1 automatically filled diagnostics in Seaflux
       ! --> Ocean sink is passed to drydep and not explicitly written out!!
       IF ( ND11 > 0 ) THEN 

          ! Get HEMCO species ID 
          SpcName   = 'ACET'
          ID1 = HCO_GetHcoID( TRIM(SpcName), HcoState )
          IF ( ID1 <= 0 ) THEN
             MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN      
          ENDIF

          ! Define diagnostics in MEGAN 
          ExtNr = GetExtNr( 'MEGAN' )
          IF ( ExtNr <= 0 ) THEN
             CALL HCO_WARNING ( 'Cannot find extension MEGAN', RC, THISLOC=LOC )
          ELSE
             ! There are three manual diagnostics in MEGAN
             DO I = 1,3

                ! Define diagnostics names. These names have to match the
                ! names called in hcox_megan_mod.F90.
                IF ( I == 1 ) THEN
                   DiagnName = 'MEGAN_ACET_MONO'
                ELSEIF ( I == 2 ) THEN
                   DiagnName = 'MEGAN_ACET_MBO'
                ELSEIF ( I == 3 ) THEN
                   DiagnName = 'MEGAN_ACET_DIRECT'
                ENDIF
   
                ! Create diagnostics. Don't use AutoFill here since the 
                ! diagnostics update calls are explicitly called in 
                ! hcox_megan_mod.F90.
                CALL Diagn_Create ( am_I_Root,                   & 
                                    HcoState,                    &
                                    cName     = TRIM(DiagnName), &
                                    ExtNr     = ExtNr,           &
                                    Cat       = -1,              &
                                    Hier      = -1,              &
                                    HcoID     = ID1,             &
                                    SpaceDim  = 2,               &
                                    LevIDx    = -1,              &
                                    OutUnit   = 'kg/m2/s',       &
                                    WriteFreq = 'Manual',        &
                                    AutoFill  = 0,               &
                                    cID       = N,               & 
                                    RC        = RC                ) 
                IF ( RC /= HCO_SUCCESS ) RETURN 
             ENDDO
          ENDIF

          ! Define diagnostics in SeaFlux 
          ExtNr = GetExtNr( 'SeaFlux' )
          IF ( ExtNr <= 0 ) THEN
             CALL HCO_WARNING ( 'Cannot find extension SeaFlux', RC, THISLOC=LOC )
          ELSE
             ! Create diagnostics. Set AutoFill to on for this diagnostics.
             DiagnName = 'AD11_OCEAN_SOURCE'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = -1,              &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 
          ENDIF
       ENDIF

       !-----------------------------------------------------------------
       ! Sulfur emissions (ND13) 
       !-----------------------------------------------------------------
       ! --> For DMS, only positive flux is diagnosed
       ! --> Don't diagnose biofuel as most inventory include it w/ anthro
       ! --> Volcano emissions are lumped (eruptive + noneruptive)
       ! ==> BIOMASS diagnostics are defined below (ND28)
       !-----------------------------------------------------------------

       IF (   ND13 > 0 .AND. &
            ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. &
              Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN

          ! DMS emissions
          ! As for acetone, only keep track of flux from ocean. Deposition
          ! from atmosphere is handled by drydep.
          ExtNr = GetExtNr( 'SeaFlux' )
          IF ( ExtNr <= 0 ) THEN
             MSG = 'Cannot find SeaFlux extension - ' // &
                   'no DMS diagnostics will be written!'
             CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
          ELSE
             ID1 = HCO_GetHcoID( 'DMS', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN      
             ENDIF
             DiagnName = 'AD13_DMS_OCEAN_SOURCE'
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
          ENDIF 

          ! SO2 emissions ...
          ExtNr = 0
          ID1   = HCO_GetHcoID( 'SO2', HcoState )
          IF ( ID1 <= 0 ) THEN
             MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN      
          ENDIF

          ! ... from aircrafts ...
          DiagnName = 'AD13_SO2_AIRCRAFT'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
                              Cat       = 20,              &
                              Hier      = -1,              &
                              HcoID     = ID1,             &
                              SpaceDim  = 3,               &
                              LevIDx    = -1,              &
                              OutUnit   = 'kg',            &
                              WriteFreq = 'Manual',        &
                              AutoFill  = 1,               &
                              cID       = N,               & 
                              RC        = RC                )
          IF ( RC /= HCO_SUCCESS ) RETURN 

          ! ... anthropogenic ...
          DiagnName = 'AD13_SO2_ANTHROPOGENIC'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
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

          ! ... biofuel ...
          DiagnName = 'AD13_SO2_BIOFUEL'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
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

          ! ... from volcanoes ...
          DiagnName = 'AD13_SO2_VOLCANO'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
                              Cat       = 50,              &
                              Hier      = -1,              &
                              HcoID     = ID1,             &
                              SpaceDim  = 3,               &
                              LevIDx    = -1,              &
                              OutUnit   = 'kg',            &
                              WriteFreq = 'Manual',        &
                              AutoFill  = 1,               &
                              cID       = N,               & 
                              RC        = RC                )
          IF ( RC /= HCO_SUCCESS ) RETURN 

          ! ... from ships ...
          DiagnName = 'AD13_SO2_SHIP'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
                              Cat       = 10,              &
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

          ! NH3 emissions ...
          ID1 = HCO_GetHcoID( 'NH3', HcoState )
          IF ( ID1 <= 0 ) THEN
             MSG = 'This is not a HEMCO species: ' // TRIM(SpcName) // &
                   ' No diagnostics written for this compound!'
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          ! ... anthropogenic ...
          ExtNr     = 0
          DiagnName = 'AD13_NH3_ANTHROPOGENIC'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
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

          ! ... biofuel ...
          ExtNr     = 0
          DiagnName = 'AD13_NH3_BIOFUEL'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
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

          ! ... natural
          ExtNr     = 0
          DiagnName = 'AD13_NH3_NATURAL'
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
                              Cat       = 3,               &
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

          ! SO4 emissions ...
          ID1 = HCO_GetHcoID( 'SO4', HcoState )
          IF ( ID1 <= 0 ) THEN
             MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          ! ... anthropogenic ...
          ExtNr     = 0
          DiagnName = 'AD13_SO4_ANTHROPOGENIC' 
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
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

          ! ... biofuel ...
          ExtNr     = 0
          DiagnName = 'AD13_SO4_BIOFUEL' 
          CALL Diagn_Create ( am_I_Root,                   & 
                              HcoState,                    &
                              cName     = TRIM(DiagnName), &
                              ExtNr     = ExtNr,           &
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
       ENDIF

       !-----------------------------------------------------------------
       ! Biomass burning emissions (ND07, ND13, ND28, ND29, ND32) 
       !-----------------------------------------------------------------
       ! ==> write one single biomass burning diagnostics per species.
       ! ==> Biomass buring comes from GFED or FINN inventory. If none
       !     of those inventories is used, it is assumed that biomass
       !     burning has category 3.
       ! ==> Diagnostics are returned in kg/m2/s.
       !-----------------------------------------------------------------
       ! Get ext. Nr of used extension
       Cat   = -1
       ExtNr = GetExtNr( 'GFED3' )
       IF ( ExtNr <= 0 ) ExtNr = GetExtNr( 'FINN' )
       IF ( ExtNr <= 0 ) THEN
          ExtNr = 0
          Cat   = 3
       ENDIF
 
       ! For ND28 diagnostics, biomass burning emissions must be
       ! turned on! 
       IF ( ExtNr <= 0 .AND. ND28 > 0 ) THEN
          MSG = 'ND28 error: No biomass burning emissions defined'
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       IF ( ExtNr > 0 ) THEN
          
          ! ND28 only:
          IF ( ND28 > 0 ) THEN
             ID1 = HCO_GetHcoID( 'ALK4', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_ALK4'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

             ID1 = HCO_GetHcoID( 'ACET', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_ACET'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

             ID1 = HCO_GetHcoID( 'MEK', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_MEK'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

             ID1 = HCO_GetHcoID( 'ALD2', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_ALD2'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

             ID1 = HCO_GetHcoID( 'PRPE', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_PRPE'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

             ID1 = HCO_GetHcoID( 'C3H8', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_C3H8'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

             ID1 = HCO_GetHcoID( 'CH2O', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_CH2O'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

             ID1 = HCO_GetHcoID( 'C2H6', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_C2H6'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF ! ND28 only

          IF ( ND28 > 0 .OR.                                         &
               ( ND13 > 0 .AND. (Input_Opt%ITS_A_FULLCHEM_SIM .OR.   &
                                 Input_Opt%ITS_AN_AEROSOL_SIM     )) &
              ) THEN
             ID1 = HCO_GetHcoID( 'SO2', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_SO2'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

             ID1 = HCO_GetHcoID( 'NH3', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_NH3'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          IF ( ND28 > 0 .OR. ( ND07 > 0 .AND. Input_Opt%LCARB ) ) THEN

             ID1 = HCO_GetHcoID( 'OC', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_OC'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

             ID1 = HCO_GetHcoID( 'BC', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_BC'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          ! CO (ND28 and/or ND29)
          IF ( ND28 > 0 .OR. ND29 > 0 ) THEN
             ID1 = HCO_GetHcoID( 'CO', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_CO'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          ! NO (ND28 and/or ND32)
          IF ( ND28 > 0 .OR. ND32 > 0 ) THEN
             ID1 = HCO_GetHcoID( 'NO', HcoState )
             IF ( ID1 <= 0 ) THEN
                MSG = 'This is not a HEMCO species: ' // TRIM(SpcName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             DiagnName = 'BIOMASS_NO'
             CALL Diagn_Create ( am_I_Root,                   & 
                                 HcoState,                    &
                                 cName     = TRIM(DiagnName), &
                                 ExtNr     = ExtNr,           &
                                 Cat       = Cat,             &
                                 Hier      = -1,              &
                                 HcoID     = ID1,             &
                                 SpaceDim  = 2,               &
                                 LevIDx    = -1,              &
                                 OutUnit   = 'kg/m2/s',       &
                                 WriteFreq = 'Manual',        &
                                 AutoFill  = 1,               &
                                 cID       = N,               & 
                                 RC        = RC                ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF
       ENDIF ! Biomass burning


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

    ! Leave w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCOI_GC_Diagn_Init
!EOC
END MODULE HCOI_GC_Diagn_Mod
