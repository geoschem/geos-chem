#ifdef MODEL_CLASSIC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_boundary_condition_mod.F90
!
! !DESCRIPTION: Module SET\_BOUNDARY\_CONDITION\_MOD sets boundary conditions
!  for the GEOS-Chem "Classic" nested-grid model.
!\\
!\\
! !INTERFACE:

MODULE Set_Boundary_Conditions_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Set_Boundary_Conditions
!
! !REMARKS:
!  This module was split for two purposes:
!  (1) Avoid subroutine creep in HCO_Utilities_GC_Mod as this is
!      purely GC code.
!  (2) Allow for future extension if handling of boundary conditions
!      will change (for example introducing rate-of-change)
!
! !REVISION HISTORY:
!  28 Jul 2023 - H.P. Lin   - Initial version
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
! !IROUTINE: set_boundary_conditions
!
! !DESCRIPTION: Subroutine SET\_BOUNDARY\_CONDITIONS sets the boundary
!  conditions using the boundary conditions read from HEMCO for nested-grid
!  simulations.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Set_Boundary_Conditions( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
   USE ErrCode_Mod,      ONLY : GC_SUCCESS, GC_Error
   USE Input_Opt_Mod,    ONLY : OptInput
   USE Species_Mod,      ONLY : Species, SpcConc
   USE State_Chm_Mod,    ONLY : ChmState
   USE State_Grid_Mod,   ONLY : GrdState
   USE Time_Mod,         ONLY : TIMESTAMP_STRING
   USE PhysConstants,    ONLY : AIRMW
   USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput),   INTENT(IN)    :: Input_Opt  ! Input options
   TYPE(GrdState),   INTENT(IN)    :: State_Grid ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState),   INTENT(INOUT) :: State_Chm  ! Chemistry State
!
! !OUTPUT PARAMETERS:
!
   INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REMARKS:
!  Split off from HEMCO code (Get\_Boundary\_Conditions) in order to be called
!  more frequently throughout timesteps.
!
! !REVISION HISTORY:
!  28 Jul 2023 - H.P. Lin    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER              :: I, J, L, N, NA     ! lon, lat, lev, spc indexes
   LOGICAL              :: Perturb_CH4_BC
   REAL(fp)             :: MW_g_CH4           ! CH4 molecular weight

   ! Strings
   CHARACTER(LEN=16)    :: STAMP
   CHARACTER(LEN=255)   :: errMsg, thisLoc

   !=================================================================
   ! SET_BOUNDARY_CONDITIONS begins here!
   !=================================================================

   ! Initialize
   RC      = GC_SUCCESS
   errMsg  = ''
   thisLoc = &
 ' -> at Set_Boundary_Conditions (in GeosCore/set_boundary_conditions_mod.F90)'


   ! We only need to get boundary conditions if this is a nested-grid
   ! simulation.  Otherwise the BoundaryCond field won't be allocated.
   IF ( .not. State_Grid%NestedGrid ) RETURN

   ! Ensure species array is in kg/kg as State_Chm%BoundaryCond is in kg/kg dry
   IF ( State_Chm%Spc_Units /= KG_SPECIES_PER_KG_DRY_AIR ) THEN
      IF ( Input_Opt%amIRoot ) THEN
          WRITE(6, '(a)') 'Unit check failure: Current units are '        // &
               UNIT_STR(State_Chm%Spc_Units)                              // &
               ', expected kg/kg dry'
      ENDIF
      errMsg = 'Unit check failure: Cannot apply nested-grid boundary '   // &
               'conditions if units are not kg/kg dry. Your run may '     // &
               ' have failed previous to this error.'
      CALL GC_Error( errMsg, RC, thisLoc )
      RETURN
   ENDIF

   !=========================================================================
   ! Loop over grid boxes and apply BCs to the specified buffer zone
   !=========================================================================
   !$OMP PARALLEL DO                                                         &
   !$OMP DEFAULT( SHARED                                                    )&
   !$OMP PRIVATE( I, J, L, N                                                )&
   !$OMP COLLAPSE( 2                                                        )
   DO NA = 1, State_Chm%nAdvect
   DO L  = 1, State_Grid%NZ

      ! Get the species ID from the advected species ID
      N = State_Chm%Map_Advect(NA)

      ! Optionally perturb the CH4 boundary conditions
      ! Use ppb values specified in geoschem_config.yml
      ! Convert to [kg/kg dry] (nbalasus, 8/31/2023)
      Perturb_CH4_BC = ( State_Chm%SpcData(N)%Info%Name == "CH4" .AND. &
                         ( Input_Opt%ITS_A_CH4_SIM .OR. Input_Opt%ITS_A_CARBON_SIM ) .AND. &
                         Input_Opt%DoPerturbCH4BoundaryConditions .AND. &
                         ( .NOT. State_Chm%IsCH4BCPerturbed ) )
      MW_g_CH4       =   State_Chm%SpcData(N)%Info%MW_g

      ! First loop over all latitudes of the nested domain
      DO J = 1, State_Grid%NY

         ! West BC
         DO I = 1, State_Grid%WestBuffer
            IF ( Perturb_CH4_BC ) THEN
               State_Chm%BoundaryCond(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N) + &
                                                 Input_Opt%CH4BoundaryConditionIncreaseWest * 1.0e-9_fp * MW_g_CH4 / AIRMW
            ENDIF
            State_Chm%Species(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
         ENDDO

         ! East BC
         DO I = (State_Grid%NX-State_Grid%EastBuffer)+1, State_Grid%NX
            IF ( Perturb_CH4_BC ) THEN
               State_Chm%BoundaryCond(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N) + &
                                                 Input_Opt%CH4BoundaryConditionIncreaseEast * 1.0e-9_fp * MW_g_CH4 / AIRMW
            ENDIF
            State_Chm%Species(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
         ENDDO

      ENDDO

      ! Then loop over the longitudes of the nested domain
      DO I = 1+State_Grid%WestBuffer,(State_Grid%NX-State_Grid%EastBuffer)

         ! South BC
         DO J = 1, State_Grid%SouthBuffer
            IF ( Perturb_CH4_BC ) THEN
               State_Chm%BoundaryCond(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N) + &
                                                 Input_Opt%CH4BoundaryConditionIncreaseSouth * 1.0e-9_fp * MW_g_CH4 / AIRMW
            ENDIF
            State_Chm%Species(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
         ENDDO

         ! North BC
         DO J = (State_Grid%NY-State_Grid%NorthBuffer)+1, State_Grid%NY
            IF ( Perturb_CH4_BC ) THEN
               State_Chm%BoundaryCond(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N) + &
                                                 Input_Opt%CH4BoundaryConditionIncreaseNorth * 1.0e-9_fp * MW_g_CH4 / AIRMW
            ENDIF
            State_Chm%Species(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
         ENDDO
      ENDDO

   ENDDO
   ENDDO
   !OMP END PARALLEL DO

   ! If the boundary conditions have already been perturbed, don't do it again
   IF ( Perturb_CH4_BC ) THEN
      State_Chm%IsCH4BCPerturbed = .TRUE.
   ENDIF

   ! Echo output. This will be at every time step,
   ! so comment this out when unnecessary.
   IF ( Input_Opt%amIRoot .and. Input_Opt%Verbose ) THEN
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, * ) 'SET_BOUNDARY_CONDITIONS: Done applying BCs at ', STAMP
   ENDIF

 END SUBROUTINE Set_Boundary_Conditions
!EOC
END MODULE Set_Boundary_Conditions_Mod
#endif
