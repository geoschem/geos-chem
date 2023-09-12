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
   INTEGER              :: I, J, L, N, NA     ! lon, lat, lev, spc indexes
   CHARACTER(LEN=255)   :: LOC                ! routine location
   CHARACTER(LEN=16)    :: STAMP

   !=================================================================
   ! SET_BOUNDARY_CONDITIONS begins here!
   !=================================================================

   ! Name of this routine
   LOC = &
 ' -> at Set_Boundary_Conditions (in GeosCore/set_boundary_conditions_mod.F90)'

   ! Assume success
   RC        = GC_SUCCESS

   ! We only need to get boundary conditions if this is a nested-grid
   ! simulation.  Otherwise the BoundaryCond field won't be allocated.
   IF ( .not. State_Grid%NestedGrid ) RETURN

   ! Ensure species array is in kg/kg as State_Chm%BoundaryCond is in kg/kg dry
   IF ( TRIM(State_Chm%Spc_Units) /= 'kg/kg dry' ) THEN
      IF ( Input_Opt%amIRoot ) THEN
          WRITE(6,*) 'Unit check failure: Current units are ', &
               State_Chm%Spc_Units, ', expected kg/kg dry'
      ENDIF
      CALL GC_Error( 'Unit check failure: Cannot apply nested-grid boundary conditions if units are not kg/kg dry. Your run may have failed previous to this error.', RC, LOC )
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

      ! First loop over all latitudes of the nested domain
      DO J = 1, State_Grid%NY

         ! West BC
         DO I = 1, State_Grid%WestBuffer
            State_Chm%Species(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
         ENDDO

         ! East BC
         DO I = (State_Grid%NX-State_Grid%EastBuffer)+1, State_Grid%NX
            State_Chm%Species(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
         ENDDO

      ENDDO

      ! Then loop over the longitudes of the nested domain
      DO I = 1+State_Grid%WestBuffer,(State_Grid%NX-State_Grid%EastBuffer)

         ! South BC
         DO J = 1, State_Grid%SouthBuffer
            State_Chm%Species(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
         ENDDO

         ! North BC
         DO J = (State_Grid%NY-State_Grid%NorthBuffer)+1, State_Grid%NY
            State_Chm%Species(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
         ENDDO
      ENDDO

   ENDDO
   ENDDO
   !OMP END PARALLEL DO

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
