#if defined( MODEL_GEOS ) || defined( MODEL_GCHP )
#include "MAPL_Generic.h"
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tracer_mod.F90
!
! !DESCRIPTION: Module TRACER\_MOD is used to implement passive tracers, 
! typically used in the TransportTracers simulation.
!\\
!\\
! !INTERFACE:
!
MODULE Tracer_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Tracer_Source_Phase
  PUBLIC  :: Tracer_Sink_Phase
!
! !PRIVATE MEMBER FUNCTIONS:
!

!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tracer_source_phase
!
! !DESCRIPTION: Subroutine TRACER\_SOURCE\_PHASE
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tracer_Source_Phase( Input_Opt, State_Chm, State_Grid, &
                                  State_Met, RC                     )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_State_GC_Mod, ONLY : HcoState
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE Input_Opt_Mod,    ONLY : OptInput
    USE PhysConstants,    ONLY : AVO
    USE Species_Mod,      ONLY : Species
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
    USE Time_Mod,         ONLY : Get_Ts_Chem
    USE UnitConv_Mod,     ONLY : Convert_Spc_Units

#if defined( MODEL_GEOS ) || defined( MODEL_GCHP )
    USE ESMF
    USE MAPL
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN   ) :: Input_Opt   ! Input options object
    TYPE(GrdState), INTENT(IN   ) :: State_Grid  ! Grid state object
    TYPE(MetState), INTENT(IN   ) :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY:
!  02 Mar 2023 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: I, J,  L
    INTEGER                :: N, DT
    REAL(fp)               :: Local_Tally
    REAL(fp)               :: Total_Area
    REAL(fp)               :: Total_Spc

    ! SAVEd scalars
    LOGICAL,  SAVE         :: First = .TRUE.

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg,  ThisLoc
    CHARACTER(LEN=63)      :: OrigUnit

    ! Arrays
    REAL(fp)               :: Flux(State_Grid%NX,State_Grid%NY)
    REAL(fp)               :: Mask(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Pointers to fields in the HEMCO data structure
    REAL(f4),      POINTER :: O3(:,:,:)   => NULL()

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

#if defined( MODEL_GEOS ) || defined( MODEL_GCHP )
    INTEGER       :: status
    TYPE(ESMF_VM) :: vm
#endif

    !========================================================================
    ! Tracer_Source_Phase begins here!
    !========================================================================

    ! Initialize
    RC          = GC_SUCCESS
    ErrMsg      = ''
    ThisLoc     = &
       ' -> at Tracer_Source_Phase (in module GeosCore/tracer_mod.F90)'
    DT          = GET_TS_CHEM()
    Local_Tally = 0.0_fp
    Total_Area  = 0.0_fp
    Total_Spc   = 0.0_fp
    Mask        = 1.0_fp
    Flux        = 0.0_fp

#if defined( MODEL_GEOS ) || defined( MODEL_GCHP )
    call ESMF_VmGetCurrent(vm, rc=status)
    _VERIFY(status)
#endif

    !=======================================================================
    ! Convert species units to v/v dry
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'v/v dry', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (kg/kg dry -> v/v dry)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Apply tracer source
    !========================================================================

    ! Loop over species
    DO N = 1, State_Chm%nAdvect

       ! Point to the Species Database entry for species N
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Skip Rn-Pb-Be tracers for now
       IF ( SpcInfo%Is_RadioNuclide ) CYCLE

       ! Skip if this species does not have a source or if the 
       ! source is handled by HEMCO
       IF ( TRIM(SpcInfo%Src_Mode) == 'none'   .or. &
            TRIM(SpcInfo%Src_Mode) == 'HEMCO' ) CYCLE

       IF ( First ) THEN
 
          !------------------------------------------------------------------
          ! Convert Src_Value to v/v or to seconds for aoa tracers
          !------------------------------------------------------------------
          IF(TRIM(SpcInfo%Src_Units) ==   'ppmv'    ) THEN

             SpcInfo%Src_Value = SpcInfo%Src_Value * 1.0E-06

          ELSE IF(TRIM(SpcInfo%Src_Units) ==   'ppbv'    ) THEN

             SpcInfo%Src_Value = SpcInfo%Src_Value * 1.0E-09

          ELSE IF(TRIM(SpcInfo%Src_Units) ==   'pptv'    ) THEN

             SpcInfo%Src_Value = SpcInfo%Src_Value * 1.0E-12

          ELSE IF ( TRIM(SpcInfo%Src_Units) == 'timestep' ) THEN

             IF ( TRIM(SpcInfo%Units) == 'hours' ) THEN

                SpcInfo%Src_Value = SpcInfo%Src_Value * DT / ( 3600_fp )

             ELSE IF(TRIM(SpcInfo%Units) == 'days'  ) THEN

                SpcInfo%Src_Value = SpcInfo%Src_Value * DT / &
                                    ( 3600_fp * 24.0_fp )

             ELSE IF(TRIM(SpcInfo%Units) == 'years' ) THEN

                SpcInfo%Src_Value = SpcInfo%Src_Value * DT / &
                                    ( 3600_fp * 24.0_fp * 365.25_fp )

             ELSE
                ErrMsg = TRIM( SpcInfo%Name ) // ': Src_Units: timestep' // &
                         ' requires species Units hours, days or years.'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          ENDIF

       ENDIF

       !---------------------------------------------------------------------
       ! Get model fields (if needed)
       !---------------------------------------------------------------------
       IF ( TRIM(SpcInfo%Src_Mode) == 'model_field' ) THEN

          ! For stOX tracer get O3 [v/v] from HEMCO
          IF ( TRIM(SpcInfo%Name ) == 'stOX' ) THEN
             CALL HCO_GetPtr( HcoState, 'GLOBAL_O3', O3, RC )
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Cannot get pointer to GLOBAL_O3!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          ELSE

             ErrMsg = 'Src_Mode: model_field not currently supported for' // &
                   TRIM( SpcInfo%Name ) // &
                   'Please modify species_database.yml or add this '      // &
                   'capability in tracer_mod.F90.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN

          ENDIF

       ENDIF

       !---------------------------------------------------------------------
       ! Determine mask
       !---------------------------------------------------------------------
       IF ( TRIM(SpcInfo%Src_Horiz) == 'lat_zone' ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J                                               )&
          !$OMP COLLAPSE( 2                                                 )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Set mask to zero outside of latitude zone
             IF ( State_Grid%YMid(I,J) < SpcInfo%Src_LatMin .and. &
                  State_Grid%YMid(I,J) > SpcInfo%Src_LatMax ) THEN
                Mask(I,J,:) = 0.0_fp
             ENDIF

          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ENDIF

       IF ( TRIM(SpcInfo%Src_Vert) == 'surface' ) THEN

          ! Set mask to zero above the surface
          MASK(:,:,2:State_Grid%NZ) = 0.0_fp

       ELSE IF ( TRIM(SpcInfo%Src_Vert) == 'pressures' ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, L                                            )&
          !$OMP COLLAPSE( 3                                                 )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Set mask to zero outside of pressure levels
             IF ( State_Met%PEDGE(I,J,L+1) < SpcInfo%Src_PresMin   .and. &
                  State_Met%PEDGE(I,J,L)   > SpcInfo%Src_PresMax ) THEN
                Mask(I,J,L) = 0.0_fp
             ENDIF

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ELSE IF ( TRIM(SpcInfo%Src_Vert) == 'troposphere' ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, L                                            )&
          !$OMP COLLAPSE( 3                                                 )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Set mask to zero outside of troposphere
             IF ( .not. State_Met%InTroposphere(I,J,L) ) &
                Mask(I,J,L) = 0.0_fp

          ENDDO
          ENDDO
          ENDDO

       ELSE IF ( TRIM(SpcInfo%Src_Vert) == 'stratosphere' ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, L                                            )&
          !$OMP COLLAPSE( 3                                                 )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Set mask to zero outside of stratosphere
             IF ( .not. State_Met%InStratosphere(I,J,L) ) &
                Mask(I,J,L) = 0.0_fp

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ENDIF

       !--------------------------------------------------------------------
       ! Apply source
       !--------------------------------------------------------------------
       IF ( TRIM(SpcInfo%Src_Mode) == 'constant' ) THEN
       
          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, L                                            )&
          !$OMP COLLAPSE( 3                                                 )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             IF ( SpcInfo%Src_Add ) THEN

                ! Add source
                State_Chm%Species(N)%Conc(I,J,L) =      &
                   State_Chm%Species(N)%Conc(I,J,L) +   &
                   ( SpcInfo%Src_Value * Mask(I,J,L) )

             ELSE

                ! Replace value
                State_Chm%Species(N)%Conc(I,J,L) =      &
                   ( SpcInfo%Src_Value * Mask(I,J,L) )

             ENDIF

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ELSE IF ( TRIM(SpcInfo%Src_Mode) == 'model_field' ) THEN

          ! Special handling for stOX tracer
          IF ( TRIM(SpcInfo%Name) == 'stOX' ) THEN

             !$OMP PARALLEL DO                                               &
             !$OMP DEFAULT( SHARED                                          )&
             !$OMP PRIVATE( I, J, L                                         )&
             !$OMP COLLAPSE( 3                                              )
             DO L = 1, State_Grid%NZ
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX
                ! Replace value
                State_Chm%Species(N)%Conc(I,J,L) = O3(I,J,L) * Mask(I,J,L)
             ENDDO
             ENDDO
             ENDDO
             !$OMP END PARALLEL DO
          ENDIF

       ELSE IF ( TRIM(SpcInfo%Src_Mode) == 'maintain_mixing_ratio' ) THEN

          ! To distrubute tracer uniformly on the surface, compute the
          ! total area
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Accumulate the total area [m2]
             Local_Tally = Local_Tally + &
                           ( State_Grid%Area_M2(I,J) * Mask(I,J,1) )

          ENDDO
          ENDDO

#if defined( MODEL_GCHP ) || defined( MODEL_GEOS )
          ! Sum across all nodes
          call MAPL_CommsAllReduceSum(vm, sendbuf=Local_Tally, recvbuf=Total_Area, cnt=1, RC=status)
#else
          Total_Area = Local_Tally
#endif

          ! Reinitialize
          Local_Tally = 0.0_fp
          
          ! Compute mol of tracer needed to achieve the desired value
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             Local_Tally = Local_Tally &
                   + ( SpcInfo%Src_Value - State_Chm%Species(N)%Conc(I,J,L) ) &
                   * ( State_Met%AIRNUMDEN(I,J,L) / AVO )                     &
                   *  State_Met%AIRVOL(I,J,L)
          ENDDO
          ENDDO
          ENDDO

#if defined( MODEL_GCHP ) || defined( MODEL_GEOS )
          ! Sum across all nodes
          call MAPL_CommsAllReduceSum(vm, sendbuf=Local_Tally, recvbuf=Total_Spc, cnt=1, __RC__)
#else
          Total_Spc = Local_Tally
#endif

          ! Compute flux [mol/m2]
          Flux(:,:) =  ( Total_Spc / Total_Area ) * Mask(:,:,1)

          ! Update species concentrations at surface [mol/mol]
          State_Chm%Species(N)%Conc(:,:,1) = State_Chm%Species(N)%Conc(:,:,1) &
               + ( Flux(:,:) * AVO            )                               &
               / ( State_Met%BXHEIGHT(:,:,1)                                  &
               *   State_Met%AIRNUMDEN(:,:,1) )

       ELSE
          ErrMsg = TRIM( SpcInfo%Name ) // ': Src_Mode '                    // &
                   TRIM( SpcInfo%Src_Mode ) // ' not currently supported. ' // &
                   'Please modify species_database.yml or add this '        // &
                   'capability in tracer_mod.F90.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Free pointers
       SpcInfo => NULL()

    ENDDO

    !=======================================================================
    ! Convert species units back to original unit
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Reset after the first time
    IF ( First ) First = .FALSE.

  END SUBROUTINE Tracer_Source_Phase
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tracer_sink_phase
!
! !DESCRIPTION: Subroutine TRACER\_SINK\_PHASE performs loss chemistry
!  on passive species with finite atmospheric lifetimes.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tracer_Sink_Phase( Input_Opt, State_Chm, State_Grid, &
                                State_Met, RC                     )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_State_GC_Mod, ONLY : HcoState
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE Input_Opt_Mod,    ONLY : OptInput
    USE PhysConstants,    ONLY : AVO
    USE Species_Mod,      ONLY : Species
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
    USE Time_Mod,         ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN   ) :: Input_Opt   ! Input options object
    TYPE(GrdState), INTENT(IN   ) :: State_Grid  ! Grid state object
    TYPE(MetState), INTENT(IN   ) :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY:
!  04 Sep 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: I, J, L, N
    REAL(fp)               :: DT
    REAL(fp)               :: DecayConstant
    REAL(fp)               :: DecayRate
    REAL(fp)               :: LO3_kg

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg,  ThisLoc

    ! Arrays
    REAL(fp)               :: Mask(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Pointers to fields in the HEMCO data structure
    REAL(f4),    POINTER   :: LOSS_O3(:,:,:)   => NULL()

    ! Parameters
    REAL(fp),    PARAMETER :: ln2 = 0.693147181E+00_fp
    REAL(fp),    PARAMETER :: Day2Sec = 86400_fp

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !========================================================================
    ! Tracer_Sink_Phase begins here!
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at Tracer_Sink_Phase (in module GeosCore/tracer_mod.F90)'
    DT       = GET_TS_CHEM()
    Mask     = 1.0_fp

    !========================================================================
    ! Apply tracer sink
    !========================================================================

    ! Loop over species
    DO N = 1, State_Chm%nAdvect

       ! Point to the Species Database entry for species N
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Skip Rn-Pb-Be tracers for now
       IF ( SpcInfo%Is_RadioNuclide ) CYCLE

       ! Skip if no sink
       IF ( TRIM(SpcInfo%Snk_Mode) == 'none' ) CYCLE

       !---------------------------------------------------------------------
       ! Compute the decay rate (if needed)
       !---------------------------------------------------------------------
       IF ( TRIM(SpcInfo%Snk_Mode) == 'efolding' .or. &
            TRIM(SpcInfo%Snk_Mode) == 'halflife' ) THEN

          IF ( TRIM(SpcInfo%Snk_Mode) == 'efolding' ) THEN

             ! Compute the decay constant (s-1) from the e-folding time
             !  ln(N/No) = ln(1/e) = (-1) =  -decayConstant * e-folding time
             DecayConstant = 1.0 / ( SpcInfo%Snk_Period * Day2Sec )

          ELSE IF ( TRIM(SpcInfo%Snk_Mode) == 'halflife' ) THEN

             ! Compute the decay constant (s-1) from the half-life:
             !  ln(N/No) = ln(1/2) = -decayConstant * halfLife
             DecayConstant = ln2 / ( SpcInfo%Snk_Period * Day2Sec )

          ENDIF

          DecayRate  = EXP( - DecayConstant * DT )

          !### Debug output
          IF ( Input_Opt%Verbose ) THEN
             WRITE( 6,100 ) ADJUSTL( SpcInfo%Name ), DecayRate
 100         FORMAT( '     -  Species name, decay rate: ', a15, es13.6 )
          ENDIF

       ENDIF

       !--------------------------------------------------------------------
       ! Get chemical loss (if needed)
       !--------------------------------------------------------------------
       IF ( TRIM(SpcInfo%Snk_Mode) == 'chemical_loss' ) THEN

          ! For stOX tracer get O3 loss [molec/cm3/s] from HEMCO
          IF ( TRIM(SpcInfo%Name ) == 'stOX' ) THEN

             CALL HCO_GetPtr( HcoState, 'O3_LOSS', LOSS_O3,   RC )
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Cannot get pointer to O3_LOSS!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          ELSE

             ErrMsg = 'Snk_Mode: chemical_loss not currently supported '  // &
                   'for ' // TRIM( SpcInfo%Name )  // &
                   'Please modify species_database.yml or add this '      // &
                   'capability in tracer_mod.F90.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN

          ENDIF

       ENDIF

       !---------------------------------------------------------------------
       ! Determine mask
       !---------------------------------------------------------------------
       IF ( TRIM(SpcInfo%Snk_Horiz) == 'lat_zone' ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J                                               )&
          !$OMP COLLAPSE( 2                                                 )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Set mask to zero outside of latitude zone
             IF ( State_Grid%YMid(I,J) < SpcInfo%Snk_LatMin .and. &
                  State_Grid%YMid(I,J) > SpcInfo%Snk_LatMax ) THEN
                Mask(I,J,:) = 0.0_fp
             ENDIF

          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ENDIF

       IF ( TRIM(SpcInfo%Snk_Vert) == 'surface' ) THEN

          ! Set mask to zero above the surface
          MASK(:,:,2:State_Grid%NZ) = 0.0_fp

       ELSE IF ( TRIM(SpcInfo%Snk_Vert) == 'boundary_layer' ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, L                                            )&
          !$OMP COLLAPSE( 3                                                 )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Set mask to zero outside the boundary layer
             IF ( L > State_Met%PBL_TOP_L(I,J) ) THEN
                Mask(I,J,L) = 0.0_fp
             ENDIF

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ELSE IF ( TRIM(SpcInfo%Snk_Vert) == 'troposphere' ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, L                                            )&
          !$OMP COLLAPSE( 3                                                 )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Set mask to zero outside of troposphere
             IF ( .not. State_Met%InTroposphere(I,J,L) ) &
                Mask(I,J,L) = 0.0_fp

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ELSE IF ( TRIM(SpcInfo%Snk_Vert) == 'stratosphere' ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, L                                            )&
          !$OMP COLLAPSE( 3                                                 )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Set mask to zero outside of stratosphere
             IF ( .not. State_Met%InStratosphere(I,J,L) ) &
                Mask(I,J,L) = 0.0_fp

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ENDIF

       !---------------------------------------------------------------------
       ! Apply loss
       !---------------------------------------------------------------------

       IF ( TRIM(SpcInfo%Snk_Mode) == 'efolding' .or. &
            TRIM(SpcInfo%Snk_Mode) == 'halflife' ) THEN

          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, L                                            )&
          !$OMP COLLAPSE( 3                                                 )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             State_Chm%Species(N)%Conc(I,J,L) =                &
               State_Chm%Species(N)%Conc(I,J,L) * ( DecayRate * Mask(I,J,L) )
          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ELSE IF ( TRIM(SpcInfo%Snk_Mode) == 'constant' ) THEN

          WHERE( Mask > 0.0_fp )
             State_Chm%Species(N)%Conc =  SpcInfo%Snk_Value
          ENDWHERE

       ELSE IF ( TRIM(SpcInfo%Snk_Mode) == 'chemical_loss' ) THEN

          ! Special handling for stOX tracer
          IF ( TRIM(SpcInfo%Name ) == 'stOX' ) THEN

             !$OMP PARALLEL DO                                               &
             !$OMP DEFAULT( SHARED                                          )&
             !$OMP PRIVATE( I, J, L, LO3_kg                                 )&
             !$OMP COLLAPSE( 3                                              )
             DO L = 1, State_Grid%NZ
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX

                ! Convert L(O3) from [molec/cm3/s] to [kg]
                LO3_kg = LOSS_O3(I,J,L)               & ! in molec/cm3/s
                 *  (SpcInfo%MW_g * 1000.0_fp) / AVO  & ! molec/cm3/s -> kg/m3/s
                 * State_Met%AIRVOL(I,J,L)            & ! kg/m3/s     -> kg/s
                 * DT                                   ! kg/s        -> kg

                State_Chm%Species(N)%Conc(I,J,L) = &
                   State_Chm%Species(N)%Conc(I,J,L) - ( LO3_kg * MASK(I,J,L) )

             ENDDO
             ENDDO
             ENDDO
             !$OMP END PARALLEL DO

          ENDIF

       ENDIF

    ENDDO

  END SUBROUTINE Tracer_Sink_Phase
!EOC
END MODULE Tracer_Mod