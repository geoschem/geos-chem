!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_utilities_gc_mod.F90
!
! !DESCRIPTION: Module hco\_utilities\_gc\_mod.F90 is the high-level component
!  used by GEOS-Chem routines to receive emissions and generally data from
!  HEMCO. It interacts with HEMCO routines directly but also GEOS-Chem data
!  structures.
!\\
!\\
!  Please note the following guidelines for inclusion of routines in this module
!  to prevent feature creep like the former HCOI\_GC\_Main\_Mod:
!  - The routines must use both HEMCO and GEOS-Chem data structures. If not,
!    they may be better suited for EMISSIONS\_MOD or FLEXGRID\_READ\_MOD,
!    or implemented as common tools in HEMCO/Interfaces/hco\_interface\_common.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_Utilities_GC_Mod
!
! !USES:
!
  USE Precision_Mod
  USE HCO_Error_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GetHcoValEmis
  PUBLIC :: GetHcoValDep
  PUBLIC :: Compute_Sflx_For_Vdiff

#if defined( MODEL_CLASSIC )
  !=========================================================================
  ! These are only needed for GEOS-Chem "Classic"
  !=========================================================================
  PUBLIC :: Get_GC_Restart
  PUBLIC :: Get_Met_Fields
  PUBLIC :: Get_Boundary_Conditions
#endif

!
! !REMARKS:
!  Mostly wrapper functions migrated from the former HCO_Interface_Mod
!
! !REVISION HISTORY:
!  12 Mar 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoValEmis
!
! !DESCRIPTION: Subroutine GetHcoVal is a wrapper routine to return an
! emission (kg/m2/s) or deposition (1/s) value from the HEMCO state object
! for a given GEOS-Chem tracer at position I, J, L.
! A value of zero is returned if no HEMCO species is defined for the given
! tracer, and the output parameter Found is set to false.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetHcoValEmis ( TrcID, I, J, L, Found, Emis )
!
! !USES:
!
    USE HCO_Interface_Common, ONLY : GetHcoVal
    USE HCO_State_GC_Mod,     ONLY : ExtState
    USE HCO_State_GC_Mod,     ONLY : HcoState
!
! !INPUT ARGUMENTS:
!
    INTEGER,            INTENT(IN   )  :: TrcID   ! GEOS-Chem tracer ID
    INTEGER,            INTENT(IN   )  :: I, J, L ! Position
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,            INTENT(  OUT)  :: Found   ! Was this tracer ID found?
    REAL(hp),           INTENT(  OUT)  :: Emis    ! Emissions  [kg/m2/s]
!
! !REMARKS:
!  This subroutine is currently just a stub to call the equivalent in HEMCO
!  utilities. This assumes that HEMCO and GEOS-Chem grids match.
!  When the intermediate grid is implemented, a regridding routine will live
!  here and regrid data on-demand.
!
!  This also assumes that TrcID matches HEMCO tracer ID. If not, a mapping
!  needs to be performed here.
!
! !REVISION HISTORY:
!  20 Oct 2014 - C. Keller - Initial Version
!  12 Mar 2020 - H.P. Lin  - Now wrapper around common utilities
!EOP
!------------------------------------------------------------------------------
!BOC
    CALL GetHcoVal( HcoState, ExtState, TrcID, I, J, L, Found, Emis=Emis )

  END SUBROUTINE GetHcoValEmis
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoValEmis
!
! !DESCRIPTION: Subroutine GetHcoVal is a wrapper routine to return an
! emission (kg/m2/s) or deposition (1/s) value from the HEMCO state object
! for a given GEOS-Chem tracer at position I, J, L.
! A value of zero is returned if no HEMCO species is defined for the given
! tracer, and the output parameter Found is set to false.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetHcoValDep ( TrcID, I, J, L, Found, Dep )
!
! !USES:
!
    USE HCO_Interface_Common, ONLY : GetHcoVal
    USE HCO_State_GC_Mod,     ONLY : ExtState
    USE HCO_State_GC_Mod,     ONLY : HcoState
!
! !INPUT ARGUMENTS:
!
    INTEGER,            INTENT(IN   )  :: TrcID   ! GEOS-Chem tracer ID
    INTEGER,            INTENT(IN   )  :: I, J, L ! Position
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,            INTENT(  OUT)  :: Found   ! Was this tracer ID found?
    REAL(hp),           INTENT(  OUT)  :: Dep     ! Deposition [1/s]
!
! !REMARKS:
!  This subroutine is currently just a stub to call the equivalent in HEMCO
!  utilities. This assumes that HEMCO and GEOS-Chem grids match.
!  When the intermediate grid is implemented, a regridding routine will live
!  here and regrid data on-demand.
!
!  This also assumes that TrcID matches HEMCO tracer ID. If not, a mapping
!  needs to be performed here.
!
! !REVISION HISTORY:
!  20 Oct 2014 - C. Keller - Initial Version
!  12 Mar 2020 - H.P. Lin  - Now wrapper around common utilities
!EOP
!------------------------------------------------------------------------------
!BOC
    CALL GetHcoVal ( HcoState, ExtState, TrcID, I, J, L, Found, Dep=Dep )

  END SUBROUTINE GetHcoValDep
!EOC
#if defined ( MODEL_CLASSIC )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_met_fields
!
! !DESCRIPTION: Subroutine GET\_MET\_FIELDS calls the various routines to get
! met fields from HEMCO.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Get_Met_Fields( Input_Opt, State_Chm, State_Grid, &
                            State_Met, Phase, RC )
!
! ! USES:
!
   USE Calc_Met_Mod
   USE ErrCode_Mod
   USE FlexGrid_Read_Mod 
   USE HCO_State_GC_Mod,  ONLY : ExtState
   USE HCO_State_GC_Mod,  ONLY : HcoState
   USE HCO_EmisList_Mod,  ONLY : HCO_GetPtr
   USE Input_Opt_Mod,     ONLY : OptInput
   USE Pressure_Mod,      ONLY : Set_Floating_Pressures
   USE State_Chm_Mod,     ONLY : ChmState
   USE State_Grid_Mod,    ONLY : GrdState
   USE State_Met_Mod,     ONLY : MetState
   USE Time_Mod
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput),   INTENT(IN   )          :: Input_Opt  ! Input options
   TYPE(GrdState),   INTENT(IN   )          :: State_Grid ! Grid State
   INTEGER,          INTENT(IN   )          :: Phase      ! Run phase
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(MetState),   INTENT(INOUT)          :: State_Met  ! Meteorology State
   TYPE(ChmState),   INTENT(INOUT)          :: State_Chm  ! Chemistry State
   INTEGER,          INTENT(INOUT)          :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY:
!  07 Feb 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER              :: N_DYN              ! Dynamic timestep in seconds
   INTEGER              :: D(2)               ! Variable for date and time
   LOGICAL              :: FOUND              ! Found in restart file?
   LOGICAL              :: Update_MR          ! Update species mixing ratio?
   CHARACTER(LEN=255)   :: v_name             ! Variable name

   ! Pointers
   REAL*4,  POINTER     :: Ptr2D(:,:)
   REAL*4,  POINTER     :: Ptr3D(:,:,:)

   !=================================================================
   !    *****  R E A D   M E T   F I E L D S    *****
   !    *****  At the start of the GEOS-Chem simulation  *****
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize pointers
   Ptr2D       => NULL()
   Ptr3D       => NULL()

   !----------------------------------
   ! Read time-invariant data (Phase 0 only)
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      CALL FlexGrid_Read_CN( Input_Opt, State_Grid, State_Met )
   ENDIF

   !----------------------------------
   ! Read 1-hr time-averaged data
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      D = GET_FIRST_A1_TIME()
   ELSE
      D = GET_A1_TIME()
   ENDIF
   IF ( PHASE == 0 .or. ITS_TIME_FOR_A1() .and. &
        .not. ITS_TIME_FOR_EXIT() ) THEN
      CALL FlexGrid_Read_A1( D(1), D(2), Input_Opt, State_Grid, State_Met )
   ENDIF

   !----------------------------------
   ! Read 3-hr time averaged data
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      D = GET_FIRST_A3_TIME()
   ELSE
      D = GET_A3_TIME()
   ENDIF
   IF ( PHASE == 0 .or. ITS_TIME_FOR_A3() .and. &
        .not. ITS_TIME_FOR_EXIT() ) THEN
      CALL FlexGrid_Read_A3( D(1), D(2), Input_Opt, State_Grid, State_Met )
   ENDIF

   !----------------------------------
   ! Read 3-hr instantanous data
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      D = GET_FIRST_I3_TIME()
      CALL FlexGrid_Read_I3_1( D(1), D(2), Input_Opt, State_Grid, State_Met )

      ! On first call, attempt to get instantaneous met fields for prior
      ! timestep from the GEOS-Chem restart file. Otherwise, initialize
      ! to met fields for this timestep.

      !-------------
      ! TMPU
      !-------------

      ! Define variable name
      v_name = 'TMPU1'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%TMPU1 = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize TMPU1    from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'TMPU1    not found in restart, keep as value at t=0'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      !-------------
      ! SPHU
      !-------------

      ! Define variable name
      v_name = 'SPHU1'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%SPHU1 = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize SPHU1    from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'SPHU1    not found in restart, keep as value at t=0'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      !-------------
      ! PS1_WET
      !-------------

      ! Define variable name
      v_name = 'PS1WET'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%PS1_WET = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize PS1_WET  from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'PS1_WET  not found in restart, keep as value at t=0'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

      !-------------
      ! PS1_DRY
      !-------------

      ! Define variable name
      v_name = 'PS1DRY'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%PS1_DRY = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize PS1_DRY  from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'PS1_DRY  not found in restart, keep as value at t=0'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

      !-------------
      ! DELP_DRY
      !-------------

      ! Define variable name
      v_name = 'DELPDRY'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%DELP_DRY = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize DELP_DRY from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'DELP_DRY not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
      ! and compute avg dry pressure near polar caps
      CALL Set_Dry_Surface_Pressure( State_Grid, State_Met, 1 )
      CALL AvgPole( State_Grid, State_Met%PS1_DRY )

      ! Compute avg moist pressure near polar caps
      CALL AvgPole( State_Grid, State_Met%PS1_WET )

      ! Initialize surface pressures prior to interpolation
      ! to allow initialization of floating pressures
      State_Met%PSC2_WET = State_Met%PS1_WET
      State_Met%PSC2_DRY = State_Met%PS1_DRY
      CALL Set_Floating_Pressures( State_Grid, State_Met, RC )

      ! Call AIRQNT to compute initial air mass quantities
      ! Do not update initial tracer concentrations since not read
      ! from restart file yet (ewl, 10/28/15)
      CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, &
                   RC, update_mixing_ratio=.FALSE. )

   ENDIF

   ! Read in I3 fields at t+3hours for this timestep
   IF ( ITS_TIME_FOR_I3() .and. .not. ITS_TIME_FOR_EXIT() ) THEN

      D = GET_I3_TIME()
      CALL FlexGrid_Read_I3_2( D(1), D(2), Input_Opt, State_Grid, State_Met )
      
      ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
      ! and compute avg dry pressure near polar caps
      CALL Set_Dry_Surface_Pressure( State_Grid, State_Met, 2 )
      CALL AvgPole( State_Grid, State_Met%PS2_DRY )

      ! Compute avg moist pressure near polar caps
      CALL AvgPole( State_Grid, State_Met%PS2_WET )

   ENDIF

 END SUBROUTINE Get_Met_Fields
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gc_restart
!
! !DESCRIPTION: Subroutine GET\_GC\_RESTART reads species concentrations
!  [mol/mol] from the GEOS-Chem restart file and uses them to initialize
!  species concentrations in [kg/kg dry]. If species data are missing from
!  the restart file, pre-configured background values are used. If using the
!  mercury simulation, additional restart data are read from file.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Get_GC_Restart( Input_Opt, State_Chm, State_Grid, &
                            State_Met, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Error_Mod
   USE HCO_State_GC_Mod,  ONLY : HcoState
   USE HCO_EmisList_Mod,  ONLY : HCO_GetPtr
   USE PhysConstants,     ONLY : AIRMW
   USE Input_Opt_Mod,     ONLY : OptInput
   USE Species_Mod,       ONLY : Species
   USE State_Chm_Mod,     ONLY : ChmState
   USE State_Grid_Mod,    ONLY : GrdState
   USE State_Met_Mod,     ONLY : MetState
   USE Time_Mod,          ONLY : Expand_Date
   USE UnitConv_Mod,      ONLY : Convert_Spc_Units
#ifdef APM
   USE APM_Init_Mod,      ONLY : APMIDS
#endif
   USE Ocean_Mercury_Mod, ONLY : Check_Ocean_Mercury
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
   TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!
!  09 Feb 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER              :: I, J, L, M, N      ! lon, lat, lev, indexes
   LOGICAL              :: FOUND              ! Found in restart file?
   CHARACTER(LEN=60)    :: Prefix             ! utility string
   CHARACTER(LEN=255)   :: LOC                ! routine location
   CHARACTER(LEN=255)   :: MSG                ! message
   CHARACTER(LEN=255)   :: v_name             ! variable name
   REAL(fp)             :: MW_g               ! species molecular weight
   REAL(fp)             :: SMALL_NUM          ! small number threshold
   CHARACTER(LEN=63)    :: OrigUnit

   ! Temporary arrays and pointers
   REAL*4,  TARGET           :: Temp2D(State_Grid%NX,State_Grid%NY)
   REAL*4,  TARGET           :: Temp3D(State_Grid%NX,State_Grid%NY, &
                                       State_Grid%NZ)
   REAL*4,  POINTER          :: Ptr2D(:,:  )
   REAL*4,  POINTER          :: Ptr3D(:,:,:)

   ! For Hg simulation
   INTEGER                   :: Num_Hg_Categories
   INTEGER                   :: Total_Hg_Id
   CHARACTER(LEN=60)         :: HgSpc
   CHARACTER(LEN=4), POINTER :: Hg_Cat_Name(:)

   ! Default background concentration
   REAL(fp)                  :: Background_VV

   ! Objects
   TYPE(Species),    POINTER :: SpcInfo

   !=================================================================
   ! READ_GC_RESTART begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize pointers
   Ptr2D       => NULL()
   Ptr3D       => NULL()
   SpcInfo     => NULL()
   Hg_Cat_Name => NULL()

   ! Name of this routine
   LOC = ' -> at Get_GC_Restart (in GeosCore/hco_utilities_gc_mod.F90)'

   ! Set minimum value threshold for [mol/mol]
   SMALL_NUM = 1.0e-30_fp

   !=================================================================
   ! If running Hg simulation, set Hg-specific local variables
   !=================================================================
   IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

      ! Set the # of tagHg categories from State_Chm
      Num_Hg_Categories   =  State_Chm%N_Hg_CATS

      ! Set variable storing names for each of the Hg categories
      Hg_Cat_Name => State_Chm%Hg_Cat_Name

      ! Set Hg species index corresponding to a given Hg category number;
      ! total is always the first category
      Total_Hg_Id   =  State_Chm%Hg0_Id_List(1)

   ENDIF

   !=================================================================
   ! Open GEOS-Chem restart file
   !=================================================================

   ! Write read message to log
   WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
   WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'

   !=================================================================
   ! Read species concentrations from NetCDF or use default
   ! background [mol/mol]; store in State_Chm%Species in [kg/kg dry]
   !=================================================================

   ! Print header for min/max concentration to log
   WRITE( 6, 110 )
110 FORMAT( 'Min and Max of each species in restart file [mol/mol]:' )

   ! Initialize species to all zeroes
   State_Chm%Species = 0.e+0_fp

   ! Loop over species
   DO N = 1, State_Chm%nSpecies

      ! Get info about this species from the species database
      SpcInfo => State_Chm%SpcData(N)%Info
      MW_g    =  SpcInfo%MW_g

      ! Define variable name
      v_name = 'SPC_' // TRIM( SpcInfo%Name )

      ! Initialize temporary array for this species and point to it
      Temp3D = 0.0_fp
      Ptr3D => Temp3D

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), &
                       Ptr3D,    RC,       FOUND=FOUND )

      ! Check if species data is in file
      IF ( FOUND ) THEN
         SpcInfo%Is_InRestart = .TRUE.
      ELSE
         SpcInfo%Is_InRestart = .FALSE.
      ENDIF

      ! If data is in file, read in as [mol/mol] and convert to
      ! [kg/kg dry]. Otherwise, set to background value [mol/mol]
      ! either stored in species database (advected species all levels and
      ! non-advected species levels in the chemistry grid) or a small number
      ! (non-advected species levels above the chemistry grid) converted to
      ! [kg/kg dry]
      IF ( SpcInfo%Is_InRestart ) THEN

         ! Print the min & max of each species as it is read from
         ! the restart file in mol/mol
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 120 ) N, TRIM( SpcInfo%Name ), &
                            MINVAL( Ptr3D ), MAXVAL( Ptr3D )
120         FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9, &
                    '  Max = ',es15.9)
         ENDIF

         ! Convert file value [mol/mol] to [kg/kg dry] for storage
         !$OMP PARALLEL DO       &
         !$OMP DEFAULT( SHARED ) &
         !$OMP PRIVATE( I, J, L )
         DO L = 1, State_Grid%NZ
         DO J = 1, State_Grid%NY
         DO I = 1, State_Grid%NX
            State_Chm%Species(I,J,L,N) = Ptr3D(I,J,L) * MW_g / AIRMW
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO

      ELSE

         ! Set species to the background value converted to [kg/kg dry]
         !$OMP PARALLEL DO       &
         !$OMP DEFAULT( SHARED ) &
         !$OMP PRIVATE( I, J, L )
         ! Loop over all grid boxes
         DO L = 1, State_Grid%NZ
         DO J = 1, State_Grid%NY
         DO I = 1, State_Grid%NX

            ! For non-advected species at levels above chemistry grid,
            ! use a small number for background
            IF ( L > State_Grid%MaxChemLev .and. &
                     .NOT. SpcInfo%Is_Advected ) THEN

               State_Chm%Species(I,J,L,N) = SMALL_NUM * MW_g / AIRMW

            ! For all other cases, use the background value
            ! stored in the species database
            ELSE

               State_Chm%Species(I,J,L,N) = SpcInfo%BackgroundVV &
                                            * MW_g / AIRMW

               ! Print to log if debugging is on
               IF ( Input_Opt%amIRoot .AND. &
                    I == 1 .AND. J == 1 .AND. L == 1 ) THEN
                  WRITE( 6, 140 ) N, TRIM( SpcInfo%Name ), SpcInfo%BackgroundVV
140               FORMAT('Species ', i3, ', ', a9, &
                         ': Use background = ', es15.9)
               ENDIF


            ENDIF

         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO

#ifdef APM
         !================================================================
         ! APM MICROPHYSICS
         !================================================================
         WRITE(*,*)'APM run does not find '// TRIM( SpcInfo%Name ),N
         IF(SpcInfo%Name(1:9)=='APMSPBIN2')THEN
            !$OMP PARALLEL DO        &
            !$OMP DEFAULT( SHARED  ) &
            !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               State_Chm%Species(I,J,L,N) = &
               State_Chm%Species(I,J,L,APMIDS%id_SO4)/20.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:9)=='APMSPBIN3')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               State_Chm%Species(I,J,L,N) = &
               State_Chm%Species(I,J,L,APMIDS%id_SO4)/20.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           !GanLuotest
           IF(SpcInfo%Name(1:10)=='APMSEABIN0')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               State_Chm%Species(I,J,L,N) = &
               State_Chm%Species(I,J,L,APMIDS%id_SALA)/9.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:10)=='APMSEABIN1')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               State_Chm%Species(I,J,L,N) = &
               State_Chm%Species(I,J,L,APMIDS%id_SALC)/10.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:10)=='APMDSTBIN1')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very low
               State_Chm%Species(I,J,L,N) = &
               SUM(State_Chm%Species(I,J,L,APMIDS%id_DST1:APMIDS%id_DST4))/6.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:8)=='APMBCBIN')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very low
               State_Chm%Species(I,J,L,N) = 1.D-30
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:8)=='APMOCBIN')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very low
               State_Chm%Species(I,J,L,N) = 1.D-30
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
#endif
      ENDIF

      ! Free pointer
      SpcInfo => NULL()

   ENDDO

   ! Set species units
   State_Chm%Spc_Units = 'kg/kg dry'

   ! If in debug mode, print out species min and max in [molec/cm3]
   IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) THEN

      ! Convert units
      PRINT *, " "
      PRINT *, "Species min and max in molec/cm3"
      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                              'molec/cm3', RC, OrigUnit=OrigUnit )

      ! Trap error
      IF ( RC /= GC_SUCCESS ) THEN
         Msg = 'Error returned from Convert_Spc_Units, call #1!'
         CALL GC_Error( Msg, RC, Loc )
         RETURN
      ENDIF

      ! Print values
      DO N = 1, State_Chm%nSpecies
         SpcInfo => State_Chm%SpcData(N)%Info
         WRITE(6,150) N, TRIM( SpcInfo%Name ),                 &
                         MINVAL( State_Chm%Species(:,:,:,N) ), &
                         MAXVAL( State_Chm%Species(:,:,:,N) )
150      FORMAT( 'Species ', i3, ', ', a9,                     &
                 ': Min = ', es15.9, ', Max = ', es15.9 )
         SpcInfo => NULL()
      ENDDO

      ! Convert units back
      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                              OrigUnit,  RC )

      ! Trap error
      IF ( RC /= GC_SUCCESS ) THEN
         Msg = 'Error returned from Convert_Spc_Units, call #2!'
         CALL GC_Error( Msg, RC, Loc )
         RETURN
      ENDIF

   ENDIF

   !=================================================================
   ! Get variables for FlexChem
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM .and. Input_Opt%LCHEM ) THEN

      ! Define variable name
      v_name = 'KPP_HVALUE'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%KPPHvalue = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize KPP H-value from restart file'
            WRITE(6,160) MINVAL( State_Chm%KPPHvalue(:,:,:) ), &
                         MAXVAL( State_Chm%KPPHvalue(:,:,:) )
160         FORMAT( 'KPP_HVALUE: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%KPPHvalue = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'KPP_HVALUE     not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

   ENDIF

   !=================================================================
   ! Get variables for Soil NOx emissions
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      ! Define variable name
      v_name = 'WETDEP_N'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%WetDepNitrogen = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize wet deposited nitrogen from restart file'
            WRITE(6,170) MINVAL( State_Chm%WetDepNitrogen(:,:) ), &
                         MAXVAL( State_Chm%WetDepNitrogen(:,:) )
170         FORMAT( 12x, '  WETDEP_N: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%WetDepNitrogen = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'WETDEP_N       not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

      ! Define variable name
      v_name = 'DRYDEP_N'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%DryDepNitrogen = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize dry deposited nitrogen from restart file'
            WRITE(6,180) MINVAL( State_Chm%DryDepNitrogen(:,:) ), &
                         MAXVAL( State_Chm%DryDepNitrogen(:,:) )
180         FORMAT( 12x, '  DRYDEP_N: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%DryDepNitrogen = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'DRYDEP_N       not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

   ENDIF

   !=================================================================
   ! Read variables for sulfate chemistry
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. &
        Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      ! Define variable name
      v_name = 'H2O2_AFTERCHEM'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%H2O2AfterChem = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize H2O2 from restart file'
            WRITE(6,190) MINVAL( State_Chm%H2O2AfterChem(:,:,:) ), &
                         MAXVAL( State_Chm%H2O2AfterChem(:,:,:) )
190         FORMAT( 12x, 'H2O2_AChem: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%H2O2AfterChem = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'H2O2_AFTERCHEM not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      ! Define variable name
      v_name = 'SO2_AFTERCHEM'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%SO2AfterChem = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize dry deposited nitrogen from restart file'
            WRITE(6,200) MINVAL( State_Chm%SO2AfterChem(:,:,:) ), &
                         MAXVAL( State_Chm%SO2AfterChem(:,:,:) )
200         FORMAT( 12x, ' SO2_AChem: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%SO2AfterChem = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'SO2_AFTERCHEM  not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

   ENDIF

   !=================================================================
   ! Read variables for UCX
   !=================================================================
   IF ( Input_Opt%LUCX ) THEN

      ! Define variable name
      v_name = 'STATE_PSC'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%STATE_PSC = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize PSC from restart for UCX'
            WRITE(6,210) MINVAL( State_Chm%STATE_PSC(:,:,:) ), &
                         MAXVAL( State_Chm%STATE_PSC(:,:,:) )
210         FORMAT( 12x, ' STATE_PSC: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         IF ( Input_OPt%amIRoot ) THEN
#ifdef ESMF_
            ! ExtData and HEMCO behave ambiguously - if the file was found
            ! but was full of zeros throughout the domain of interest, it
            ! will result in the same output from ExtData as if the field
            ! was missing from the file. As such, HEMCO cannot distinguish
            ! between a missing file and a field of zeros
            WRITE(6,*) 'PSC restart either all zeros in the '
            WRITE(6,*) 'root domain, or the restart file did '
            WRITE(6,*) 'not contain STATE_PSC. Root domain '
            WRITE(6,*) 'will be initialized PSC-free'
         ENDIF
#else
            WRITE(6,*) 'STATE_PSC      not found in restart, initialize PSC-free'
         ENDIF
#endif
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

   ENDIF

   !=================================================================
   ! Read ocean mercury variables
   !=================================================================
   IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

      ! Print total mass to log
      WRITE( 6, 220 )
220   FORMAT(/, 'Total mass of each ocean and snow Hg species:')

      !--------------------------------------------------------------
      ! Total Hg in ocean
      !--------------------------------------------------------------
      DO M = 1, 3

         ! Define variable name
         SELECT CASE( M )
         CASE ( 1 )
            HgSpc    = 'Hg0'
         CASE ( 2 )
            HgSpc    = 'Hg2'
         CASE ( 3 )
            HgSpc    = 'HgP'
         END SELECT
         v_name = 'OCEAN_' // TRIM( HgSpc )

         ! Get variable from HEMCO and store in local array
         CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

         ! Check if variable is in file
         IF ( FOUND ) THEN

            ! Check for negative concentrations (jaf, 7/6/11)
            DO I = 1, State_Grid%NX
            DO J = 1, State_Grid%NY
               IF ( Ptr2D(I,J) < 0.0d4 ) THEN
                  Ptr2D(I,J) = 0.0d4
               ENDIF
            ENDDO
            ENDDO

            ! Assign ocean mercury data and write total mass to log file
            SELECT CASE( M )
            CASE ( 1 )
               State_Chm%OceanHg0(:,:,Total_Hg_Id) = Ptr2D
               WRITE( 6, 240 ) TRIM( v_name ), &
                            SUM( State_Chm%OceanHg0(:,:,Total_Hg_Id) ), 'kg'
            CASE ( 2 )
               State_Chm%OceanHg2(:,:,Total_Hg_Id) = Ptr2D
               WRITE( 6, 240 ) TRIM( v_name ),  &
                            SUM( State_Chm%OceanHg2(:,:,Total_Hg_Id) ), 'kg'
            CASE ( 3 )
               State_Chm%OceanHgP(:,:,Total_Hg_Id) = Ptr2D
               WRITE( 6, 240 ) TRIM( v_name ),  &
                            SUM( State_Chm%OceanHgP(:,:,Total_Hg_Id) ), 'kg'
            END SELECT

         ELSE
            WRITE( 6, 230 ) TRIM( v_name )
         ENDIF

         ! Nullify pointer
         Ptr2D => NULL()

      ENDDO

      !--------------------------------------------------------------
      ! Additional tagged ocean Hg species
      !--------------------------------------------------------------
      IF ( Input_Opt%LSPLIT ) THEN
         DO M = 1, 3
            DO N = 2, Num_Hg_Categories

               ! Define variable name. Include appended region.
               SELECT CASE( M )
               CASE ( 1 )
                  HgSpc = 'Hg0'
               CASE ( 2 )
                  HgSpc = 'Hg2'
               CASE ( 3 )
                  HgSpc = 'HgP'
               END SELECT
               v_name = 'OCEAN_' // TRIM( HgSpc ) //  &
                        '_'      // TRIM( Hg_Cat_Name(N) )

               ! Get variable from HEMCO and store in local array
               CALL HCO_GetPtr( HcoState, TRIM(v_name), &
                                Ptr2D, RC, FOUND=FOUND )

               ! Check if variable is in file
               IF ( FOUND ) THEN

                  ! Assign ocean mercury data and write total mass to log
                  SELECT CASE( M )
                  CASE ( 1 )
                     State_Chm%OceanHg0(:,:,N) = Ptr2D
                     WRITE( 6, 240 ) TRIM( v_name ),  &
                                     SUM( State_Chm%OceanHg0(:,:,N) ), 'kg'
                  CASE ( 2 )
                     State_Chm%OceanHg2(:,:,N) = Ptr2D
                     WRITE( 6, 240 ) TRIM( v_name ),  &
                                     SUM( State_Chm%OceanHg2(:,:,N) ), 'kg'
                  CASE ( 3 )
                     State_Chm%OceanHgP(:,:,N) = Ptr2D
                     WRITE( 6, 240 ) TRIM( v_name ),  &
                                     SUM( State_Chm%OceanHgP(:,:,N) ), 'kg'
                  END SELECT

               ELSE
                  WRITE( 6, 230 ) TRIM( v_name )
               ENDIF

               ! Nullify pointer
               Ptr2D => NULL()

            ENDDO
         ENDDO

         ! Make sure tagged & total species sum up
         IF ( Input_Opt%USE_CHECKS ) THEN
            CALL CHECK_OCEAN_MERCURY( State_Chm, State_Grid, &
                                      'end of READ_GC_RESTART' )
         ENDIF
      ENDIF

      !--------------------------------------------------------------
      ! Hg snowpack on land and ocean
      !--------------------------------------------------------------
      DO M = 1, 4
         DO N = 1, Num_Hg_Categories

            ! Define variable name prefix
            SELECT CASE( M )
            CASE ( 1 )
               Prefix = 'SNOW_HG_OCEAN'        ! Reducible on ocean
            CASE ( 2 )
               Prefix = 'SNOW_HG_OCEAN_STORED' ! Non-reducible on ocean
            CASE ( 3 )
               Prefix = 'SNOW_HG_LAND'         ! Reducible on land
            CASE ( 4 )
               Prefix = 'SNOW_HG_LAND_STORED'  ! Non-reducible on land
            END SELECT

            IF ( N == 1 ) THEN
               v_name = TRIM( Prefix )
            ELSE
               ! Append category name if tagged
               v_name = TRIM( Prefix         ) // '_' // &
                        TRIM( Hg_Cat_Name(N) )
            ENDIF

            ! Get variable from HEMCO and store in local array
            CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

            ! Check if variable is in file
            IF ( FOUND ) THEN

               ! Assign ocean mercury data and write total mass to file
               SELECT CASE( M )
               CASE ( 1 )
                  State_Chm%SnowHgOcean(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgOcean(:,:,N) ), 'kg'
               CASE ( 2 )
                  State_Chm%SnowHgOceanStored(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgOceanStored(:,:,N) ),'kg'
               CASE ( 3 )
                  State_Chm%SnowHgLand(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgLand(:,:,N) ), 'kg'
               CASE ( 4 )
                  State_Chm%SnowHgLandStored(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgLandStored(:,:,N) ), 'kg'
               END SELECT

            ELSE
               WRITE( 6, 230 ) TRIM( v_name )
            ENDIF

            ! Nullify pointer
            Ptr2D => NULL()

         ENDDO
      ENDDO

      ! Format strings
230   FORMAT( a24, ' not found in restart file, set to zero')
240   FORMAT( a24, ':   ', es15.9, 1x, a4)

      ! Print note that variables are initialized to zero if not
      ! found (currently only happens in tagged Hg simulation)
      IF ( Input_Opt%LSPLIT ) THEN
         WRITE( 6, 250 )
250      FORMAT( /, 'NOTE: all variables not found in restart ', &
                    'are initialized to zero')
      ENDIF

      ! Free pointers for Hg indexing
      Hg_Cat_Name => NULL()

   ENDIF

   !=================================================================
   ! Clean up
   !=================================================================

   ! Mark end of section in log
   IF ( Input_Opt%LPRT .AND. Input_Opt%amIRoot ) THEN
      CALL DEBUG_MSG('### DONE GET_GC_RESTART')
   ENDIF
   WRITE( 6, '(a)' ) REPEAT( '=', 79 )

 END SUBROUTINE Get_GC_Restart
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_boundary_conditions
!
! !DESCRIPTION: Subroutine GET\_BOUNDARY\_CONDITIONS calls the various routines
! to get boundary conditions from HEMCO for nested grid simulations.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Get_Boundary_Conditions( Input_Opt, State_Chm, State_Grid, &
                                     State_Met, YYYYMMDD,  HHMMSS, RC )
!
! ! USES:
!
   USE ErrCode_Mod
   USE HCO_State_GC_Mod, ONLY : HcoState
   USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
   USE Input_Opt_Mod,    ONLY : OptInput
   USE PhysConstants,    ONLY : AIRMW
   USE Species_Mod,      ONLY : Species
   USE State_Chm_Mod,    ONLY : ChmState
   USE State_Grid_Mod,   ONLY : GrdState
   USE State_Met_Mod,    ONLY : MetState
   USE Time_Mod
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput),   INTENT(IN   )          :: Input_Opt  ! Input options
   TYPE(GrdState),   INTENT(IN   )          :: State_Grid ! Grid State
   INTEGER,          INTENT(IN   )          :: YYYYMMDD   ! GMT date
   INTEGER,          INTENT(IN   )          :: HHMMSS     ! GMT time
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(MetState),   INTENT(INOUT)          :: State_Met  ! Meteorology State
   TYPE(ChmState),   INTENT(INOUT)          :: State_Chm  ! Chemistry State
   INTEGER,          INTENT(INOUT)          :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Apr 2019 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER              :: I, J, L, N, NA     ! lon, lat, lev, spc indexes
   INTEGER              :: t_index            ! Time index
   LOGICAL              :: FOUND              ! Found in restart file?
   LOGICAL, SAVE        :: FIRST = .TRUE.     ! Is this the first routine call?
   CHARACTER(LEN=60)    :: Prefix             ! utility string
   CHARACTER(LEN=255)   :: LOC                ! routine location
   CHARACTER(LEN=255)   :: MSG                ! message
   CHARACTER(LEN=255)   :: v_name             ! variable name
   REAL(fp)             :: MW_g               ! species molecular weight
   CHARACTER(LEN=16)    :: STAMP

   ! Temporary arrays and pointers
   REAL*4,  TARGET      :: Temp3D(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
   REAL*4,  POINTER     :: Ptr3D(:,:,:)
   REAL(fp), POINTER    :: Spc(:,:,:,:)

   ! Objects
   TYPE(Species), POINTER :: SpcInfo

   !=================================================================
   ! GET_BOUNDARY_CONDITIONS begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! We only need to get boundary conditions if this is a nested-grid
   ! simulation.  Otherwise the BoundaryCond field won't be allocated.
   IF ( .not. State_Grid%NestedGrid ) RETURN

   ! Initialize pointers
   Ptr3D     => NULL()
   SpcInfo   => NULL()

   ! Point to species array [kg/kg]
   Spc       => State_Chm%Species

   ! Name of this routine
   LOC = ' -> at Get_Boundary_Conditions (in GeosCore/hco_utilities_gc_mod.F90)'

   ! Find the proper time-slice to read from disk
   t_index = ( HHMMSS / 030000 ) + 1

   ! Stop w/ error if the time index is invalid
   IF ( t_index < 1 .or. t_index > 8 ) THEN
      WRITE( MSG, 100 ) t_index
100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
      CALL GC_Error( MSG, RC, LOC)
      RETURN
   ENDIF

   !=================================================================
   ! Read species concentrations from NetCDF [mol/mol] and
   ! store in State_Chm%BoundaryCond in [kg/kg dry]
   !=================================================================

   ! Print header for min/max concentration to log
   IF ( Input_Opt%amIRoot ) THEN
      WRITE( 6, 110 )
110   FORMAT( 'Min and Max of each species in BC file [mol/mol]:' )
   ENDIF

   ! Initialize BCs to all zeroes
   State_Chm%BoundaryCond = 0.e+0_fp

   ! Loop over advected species
   DO NA = 1, State_Chm%nAdvect

      ! Get the species ID from the advected species ID
      N = State_Chm%Map_Advect(NA)

      ! Get info about this species from the species database
      SpcInfo => State_Chm%SpcData(N)%Info
      MW_g    =  SpcInfo%MW_g

      ! Define variable name
      v_name = 'BC_' // TRIM( SpcInfo%Name )

      ! Initialize temporary array for this species and point to it
      Temp3D = 0.0_fp
      Ptr3D => Temp3D

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, &
                       TIDX=t_index, FOUND=FOUND )

      ! Check if BCs are found
      IF ( FOUND ) THEN

         ! Print the min & max of each species as it is read from
         ! the BC file in mol/mol if debug is turned on in input.geos
         IF ( Input_Opt%amIRoot ) THEN
            IF ( FIRST .or. Input_Opt%LPRT ) THEN
               WRITE( 6, 120 ) N, TRIM( SpcInfo%Name ), &
                               MINVAL( Ptr3D ), MAXVAL( Ptr3D )
120            FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9, &
                       '  Max = ',es15.9)
            ENDIF
         ENDIF

         ! Copy data from file to State_Chm%BoundaryCond
         ! and convert from [mol/mol] to [kg/kg dry]
         State_Chm%BoundaryCond(:,:,:,N) = Ptr3D(:,:,:) * MW_g / AIRMW

      ELSE

         ! Print to log if debug is turned on in input.geos
         IF ( Input_Opt%amIRoot ) THEN
            IF ( FIRST .or. Input_Opt%LPRT ) THEN
               WRITE( 6, 130 ) N, TRIM( SpcInfo%Name ), SpcInfo%BackgroundVV
130            FORMAT('Species ', i3, ', ', a9, ': Use background = ', es15.9)
            ENDIF
         ENDIF

         ! Use the background value stored in the species database
         State_Chm%BoundaryCond(:,:,:,N) = SpcInfo%BackgroundVV &
                                            * MW_g / AIRMW

      ENDIF


      ! Loop over grid boxes and apply BCs to the specified buffer zone
      !$OMP PARALLEL DO       &
      !$OMP DEFAULT( SHARED ) &
      !$OMP PRIVATE( I, J, L )
      DO L = 1, State_Grid%NZ

         ! First loop over all latitudes of the nested domain
         DO J = 1, State_Grid%NY

            ! West BC
            DO I = 1, State_Grid%WestBuffer
               State_Chm%Species(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N)
            ENDDO

            ! East BC
            DO I = (State_Grid%NX-State_Grid%EastBuffer)+1, State_Grid%NX
               State_Chm%Species(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N)
            ENDDO

         ENDDO

         ! Then loop over the longitudes of the nested domain
         DO I = 1+State_Grid%WestBuffer,(State_Grid%NX-State_Grid%EastBuffer)

            ! South BC
            DO J = 1, State_Grid%SouthBuffer
               Spc(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N)
            ENDDO

            ! North BC
            DO J = (State_Grid%NY-State_Grid%NorthBuffer)+1, State_Grid%NY
               Spc(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N)
            ENDDO
         ENDDO

      ENDDO
      !OMP END PARALLEL DO

      ! Free pointer
      SpcInfo => NULL()

   ENDDO

   ! Reset FIRST flag
   FIRST = .FALSE.

   ! Echo output
   IF ( Input_Opt%amIRoot ) THEN
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 140 ) STAMP
140   FORMAT( 'GET_BOUNDARY_CONDITIONS: Done applying BCs at ', a )
   ENDIF

 END SUBROUTINE Get_Boundary_Conditions
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Compute_Sflx_for_Vdiff
!
! !DESCRIPTION: Computes the surface flux (\= emissions - drydep) for the
!  non-local PBL mixing.  This code was removed from within the non-local
!  PBL mixing driver routine VDIFFDR.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Sflx_for_Vdiff( Input_Opt,  State_Chm, State_Diag,      &
                                     State_Grid, State_Met, RC              )
!
! ! USES:
!
    USE Depo_Mercury_Mod,     ONLY : Add_Hg2_DD
    USE Depo_Mercury_Mod,     ONLY : Add_HgP_DD
    USE Depo_Mercury_Mod,     ONLY : Add_Hg2_SnowPack
    USE ErrCode_Mod
    USE Get_Ndep_Mod,         ONLY : Soil_Drydep
    USE Global_CH4_Mod,       ONLY : CH4_Emis
    USE HCO_Interface_Common, ONLY : GetHcoDiagn
    USE HCO_EmisList_Mod,     ONLY : HCO_GetPtr
    USE HCO_State_GC_Mod,     ONLY : ExtState
    USE HCO_State_GC_Mod,     ONLY : HcoState   
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Mercury_Mod,          ONLY : HG_Emis
    USE PhysConstants
    USE Species_Mod,          ONLY : Species
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Chm_Mod,        ONLY : Ind_
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE Time_Mod,             ONLY : Get_Ts_Conv
    USE Time_Mod,             ONLY : Get_Ts_Emis
    USE UnitConv_Mod,         ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(GrdState),   INTENT(IN)    :: State_Grid  ! Grid State
    TYPE(MetState),   INTENT(IN)    :: State_Met   ! Meteorology State
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag  ! Diagnostics State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  18 May 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE           :: first    = .TRUE.
    INTEGER, SAVE           :: id_O3    = -1
    INTEGER, SAVE           :: id_HNO3  = -1

    ! Scalars
    LOGICAL                 :: found,   zeroHg0Dep
    INTEGER                 :: I,       J
    INTEGER                 :: L,       NA
    INTEGER                 :: ND,      N
    INTEGER                 :: Hg_Cat,  topMix
    INTEGER                 :: S
    REAL(fp)                :: dep,     emis
    REAL(fp)                :: MW_kg,   fracNoHg0Dep
    REAL(fp)                :: tmpFlx

    ! Strings
    CHARACTER(LEN=63)       :: origUnit
    CHARACTER(LEN=255)      :: errMsg,  thisLoc

    ! Arrays
    REAL(fp), TARGET        :: eflx(State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Chm%nAdvect                       )
    REAL(fp), TARGET        :: dflx(State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Chm%nAdvect                       )

    ! Pointers and Objects
    REAL(fp),       POINTER :: spc(:,:,:)
    REAL(f4), SAVE, POINTER :: PNOxLoss_O3  (:,:) => NULL()
    REAL(f4), SAVE, POINTER :: PNoxLoss_HNO3(:,:) => NULL()
    TYPE(Species),  POINTER :: ThisSpc

    !=======================================================================
    ! Compute_Sflx_For_Vdiff begins here!
    !=======================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    dflx    =  0.0_fp
    eflx    =  0.0_fp
    spc     => State_Chm%Species(:,:,1,1:State_Chm%nAdvect)
    ThisSpc => NULL()
    errMsg  = ''
    thisLoc = &
    ' -> at Compute_Sflx_for_Vdiff (in module GeosCore/hco_utilities_gc_mod.F90)'

    !=======================================================================
    ! Convert units to v/v dry
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met,     &
                            'v/v dry', RC,        OrigUnit=OrigUnit         )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "Convert_Spc_Units" (to v/v dry)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! First-time setup: Get pointers to the PARANOX loss fluxes.
    ! These are stored in diagnostics 'PARANOX_O3_DEPOSITION_FLUX' and
    ! 'PARANOX_HNO3_DEPOSITION_FLUX'. The call below links pointers
    ! PNOXLOSS_O3 and PNOXLOSS_HNO3 to the data values stored in the
    ! respective diagnostics. The pointers will remain unassociated if
    ! the diagnostics do not exist.
    !=======================================================================
    IF ( FIRST ) THEN

       ! Get species IDs
       id_O3   = Ind_('O3'  )
       id_HNO3 = Ind_('HNO3')

#if !defined( MODEL_CESM )
       IF ( id_O3 > 0 ) THEN
          CALL GetHcoDiagn(                                                  &
               HcoState       = HcoState,                                    &
               ExtState       = ExtState,                                    &
               DiagnName      = 'PARANOX_O3_DEPOSITION_FLUX',                &
               StopIfNotFound = .FALSE.,                                     &
               Ptr2D          = PNOxLoss_O3,                                 &
               RC             = RC                                          )
       ENDIF

       IF ( id_HNO3 > 0 ) THEN
          CALL GetHcoDiagn(                                                  &
               HcoState       = HcoState,                                    &
               ExtState       = ExtState,                                    &
               DiagnName      = 'PARANOX_HNO3_DEPOSITION_FLUX',              &
               StopIfNotFound = .FALSE.,                                     &
               Ptr2D          = PNOxLoss_HNO3,                               &
               RC             = RC                                          )
       ENDIF
#endif

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !=======================================================================
    ! Add emissions & deposition values calculated in HEMCO.
    ! Here we only consider emissions below the PBL top.
    !
    ! For the full-chemistry simulations, emissions above the PBL
    ! top will be applied in routine SETEMIS, which occurs just
    ! before the SMVGEAR/KPP solvers are invoked.
    !
    ! For the specialty simulations, emissions above the PBL top
    ! will be applied in the chemistry routines for each
    ! specialty simulation.
    !
    ! For more information, please see this wiki page:
    ! http://wiki.geos-chem.org/Distributing_emissions_in_the_PBL
    !========================================================================
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED )                                                  &
    !$OMP PRIVATE( I,       J,            topMix,     NA,     N             )&
    !$OMP PRIVATE( thisSpc, tmpflx,       found,      emis,   dep           )&
    !$OMP PRIVATE( ND,      fracNoHg0Dep, zeroHg0Dep, Hg_cat                )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

#if !defined( MODEL_CESM )
       ! PBL top level [integral model levels]
       topMix = MAX( 1, FLOOR( State_Met%PBL_TOP_L(I,J) ) )

       ! Loop over advected species
       DO NA = 1, State_Chm%nAdvect

          ! Get the modelId
          N = State_Chm%Map_Advect(NA)

          ! Point to the corresponding entry in the species database
          ThisSpc => State_Chm%SpcData(N)%Info

          !------------------------------------------------------------------
          ! Add total emissions in the PBL to the EFLX array
          ! which tracks emission fluxes.  Units are [kg/m2/s].
          !------------------------------------------------------------------
          IF ( Input_Opt%ITS_A_CH4_SIM ) THEN

             ! CH4 emissions become stored in CH4_EMIS in global_ch4_mod.F90.
             ! We use CH4_EMIS here instead of the HEMCO internal emissions
             ! only to make sure that total CH4 emissions are properly defined
             ! in a multi-tracer CH4 simulation. For a single-tracer simulation
             ! and/or all other source types, we could use the HEMCO internal
             ! values set above and would not need the code below.
             ! Units are already in kg/m2/s. (ckeller, 10/21/2014)
             !
             !%%% NOTE: MAYBE THIS CAN BE REMOVED SOON (bmy, 5/18/19)%%%
             eflx(I,J,NA) = CH4_EMIS(I,J,NA)

          ELSE IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

             ! HG emissions become stored in HG_EMIS in mercury_mod.F90.
             ! This is a workaround to ensure backwards compatibility.
             ! Units are already in kg/m2/s. (ckeller, 10/21/2014)
             !
             !%%% NOTE: MAYBE THIS CAN BE REMOVED SOON (bmy, 5/18/19)%%%
             eflx(I,J,NA) = HG_EMIS(I,J,NA)

          ELSE

             ! Compute emissions for all other simulation
             tmpFlx = 0.0_fp
             DO L = 1, topMix
                CALL GetHcoValEmis( NA, I, J, L, found, emis )
                IF ( .NOT. found ) EXIT
                tmpFlx = tmpFlx + emis
             ENDDO
             eflx(I,J,NA) = eflx(I,J,NA) + tmpFlx

          ENDIF

          !------------------------------------------------------------------
          ! Also add drydep frequencies calculated by HEMCO (e.g. from the
          ! air-sea exchange module) to DFLX.  These values are stored
          ! in 1/s.  They are added in the same manner as the drydep freq values
          ! from drydep_mod.F90.  DFLX will be converted to kg/m2/s later.
          ! (ckeller, 04/01/2014)
          !------------------------------------------------------------------
          CALL GetHcoValDep( NA, I, J, L, found, dep )
          IF ( found ) THEN
             dflx(I,J,NA) = dflx(I,J,NA)                                     &
                          + ( dep * spc(I,J,NA) / (AIRMW / ThisSpc%MW_g)  )
          ENDIF

          ! Free pointers
          ThisSpc => NULL()
       ENDDO
#endif

       !=====================================================================
       ! Apply dry deposition frequencies
       ! These are the frequencies calculated in drydep_mod.F90
       ! The HEMCO drydep frequencies (from air-sea exchange and
       ! PARANOX) were already added above.
       !
       ! NOTES:
       ! (1) Loops over only the drydep species
       ! (2) If drydep is turned off, nDryDep=0 and the loop won't execute
       ! (3) Tagged species are included in this loop. via species database
       !=====================================================================
       DO ND = 1, State_Chm%nDryDep

          ! Get the species ID from the drydep ID
          N = State_Chm%Map_DryDep(ND)

          IF ( N <= 0 ) CYCLE

          ! Point to the corresponding Species Database entry
          ThisSpc => State_Chm%SpcData(N)%Info

          ! only use the lowest model layer for calculating drydep fluxes
          ! given that spc is in v/v
          dflx(I,J,N) = dflx(I,J,N) &
                      + State_Chm%DryDepFreq(I,J,ND) * spc(I,J,N)             &
                      /  ( AIRMW                    / ThisSpc%MW_g        )

          IF ( Input_Opt%ITS_A_MERCURY_SIM .and. ThisSpc%Is_Hg0 ) THEN

             ! Hg(0) exchange with the ocean is handled by ocean_mercury_mod
             ! so disable deposition over water here.
             ! Turn off Hg(0) deposition to snow and ice because we haven't yet
             ! included emission from these surfaces and most field studies
             ! suggest Hg(0) emissions exceed deposition during sunlit hours.
             fracNoHg0Dep = MIN( State_Met%FROCEAN(I,J) + &
                                 State_Met%FRSNO(I,J)   + &
                                 State_Met%FRLANDIC(I,J), 1e+0_fp)
             zeroHg0Dep   = ( fracNoHg0Dep > 0e+0_fp )

             IF ( zeroHg0Dep ) THEN
                dflx(I,J,N) = dflx(I,J,N) * MAX( 1.0_fp-fracNoHg0Dep, 0.0_fp )
             ENDIF
          ENDIF

          ! Free species database pointer
          ThisSpc => NULL()
       ENDDO

       !=====================================================================
       ! Convert DFLX from 1/s to kg/m2/s
       !
       ! If applicable, add PARANOX loss to this term. The PARANOX
       ! loss term is already in kg/m2/s. PARANOX loss (deposition) is
       ! calculated for O3 and HNO3 by the PARANOX module, and data is
       ! exchanged via the HEMCO diagnostics.  The data pointers PNOXLOSS_O3
       ! and PNOXLOSS_HNO3 have been linked to these diagnostics at the
       ! beginning of this routine (ckeller, 4/10/15).
       !=====================================================================
       dflx(I,J,:) = dflx(I,J,:) * State_Met%AD(I,J,1)                        &
                                 / State_Grid%Area_M2(I,J)

       IF ( ASSOCIATED( PNOxLoss_O3 ) .AND. id_O3 > 0 ) THEN
          dflx(I,J,id_O3) = dflx(I,J,id_O3) + PNOxLoss_O3(I,J)
       ENDIF

       IF ( ASSOCIATED( PNOXLOSS_HNO3 ) .AND. id_HNO3 > 0 ) THEN
          dflx(I,J,id_HNO3) = dflx(I,J,id_HNO3) + PNOxLOss_HNO3(I,J)
       ENDIF

       !=====================================================================
       ! Surface flux (SFLX) = emissions (EFLX) - dry deposition (DFLX)
       !
       ! SFLX is what we need to pass into routine VDIFF
       !=====================================================================
       State_Chm%SurfaceFlux(I,J,:) = eflx(I,J,:) - dflx(I,J,:) ! kg/m2/s

       !=====================================================================
       ! Archive Hg deposition for surface reservoirs (cdh, 08/28/09)
       !=====================================================================
       IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

          ! Loop over only the drydep species
          ! If drydep is turned off, nDryDep=0 and the loop won't execute
          DO ND = 1, State_Chm%nDryDep

             ! Get the species ID from the drydep ID
             N = State_Chm%Map_DryDep(ND)

             ! Point to the Species Database entry for tracer N
             ThisSpc => State_Chm%SpcData(N)%Info

             ! Deposition mass, kg
             dep = dflx(I,J,N) * State_Grid%Area_M2(I,J) * GET_TS_CONV()

             IF ( ThisSpc%Is_Hg2 ) THEN

                ! Get the category number for this Hg2 tracer
                Hg_Cat = ThisSpc%Hg_Cat

                ! Archive dry-deposited Hg2
                CALL ADD_Hg2_DD      ( I, J, Hg_Cat, dep                    )
                CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat, dep,                    &
                                       State_Met, State_Chm, State_Diag     )

             ELSE IF ( ThisSpc%Is_HgP ) THEN

                ! Get the category number for this HgP tracer
                Hg_Cat = ThisSpc%Hg_Cat

                ! Archive dry-deposited HgP
                CALL ADD_HgP_DD      ( I, J, Hg_Cat, dep                    )
                CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat, dep,                    &
                                       State_Met, State_Chm, State_Diag     )

             ENDIF

             ! Free pointer
             ThisSpc => NULL()
          ENDDO
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !### Uncomment for debug output
    !WRITE( 6, '(a)' ) 'eflx and dflx values HEMCO [kg/m2/s]'
    !DO NA = 1, State_Chm%nAdvect
    !   WRITE(6,*) 'eflx TRACER ', NA, ': ', SUM(eflx(:,:,NA))
    !   WRITE(6,*) 'dflx TRACER ', NA, ': ', SUM(dflx(:,:,NA))
    !   WRITE(6,*) 'sflx TRACER ', NA, ': ', SUM(State_Chm%SurfaceFlux(:,:,NA))
    !ENDDO

    !=======================================================================
    ! DIAGNOSTICS: Compute drydep flux loss due to mixing [molec/cm2/s]
    !
    ! NOTE: Dry deposition of "tagged" species (e.g. in tagO3, tagCO, tagHg
    ! specialty simulations) are accounted for in species 1..nDrydep,
    ! so we don't need to do any further special handling.
    !=======================================================================
    IF ( Input_Opt%LGTMM              .or. Input_Opt%LSOILNOX          .or.  &
         State_Diag%Archive_DryDepMix .or. State_Diag%Archive_DryDep ) THEN

       ! Loop over only the drydep species
       ! If drydep is turned off, nDryDep=0 and the loop won't execute
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( ND, N, ThisSpc, MW_kg, tmpFlx, S                      )
       DO ND = 1, State_Chm%nDryDep

          ! Get the species ID from the drydep ID
          N = State_Chm%Map_DryDep(ND)

          ! Skip if not a valid species
          IF ( N <= 0 ) CYCLE

          ! Point to the Species Database entry for this tracer
          ! NOTE: Assumes a 1:1 tracer index to species index mapping
          ThisSpc => State_Chm%SpcData(N)%Info

          ! Get the molecular weight of the species in kg
          MW_kg = ThisSpc%MW_g * 1.e-3_fp

          !-----------------------------------------------------------------
          ! HISTORY: Update dry deposition flux loss [molec/cm2/s]
          !
          ! DFLX is in kg/m2/s.  We convert to molec/cm2/s by:
          !
          ! (1) multiplying by 1e-4 cm2/m2        => kg/cm2/s
          ! (2) multiplying by ( AVO / MW_KG )    => molec/cm2/s
          !
          ! The term AVO/MW_kg = (molec/mol) / (kg/mol) = molec/kg
          !
          ! NOTE: we don't need to multiply by the ratio of TS_CONV /
          ! TS_CHEM, as the updating frequency for HISTORY is determined
          ! by the "frequency" setting in the "HISTORY.rc"input file.
          ! The old bpch diagnostics archived the drydep due to chemistry
          ! every chemistry timestep = 2X the dynamic timestep.  So in
          ! order to avoid double-counting the drydep flux from mixing,
          ! you had to multiply by TS_CONV / TS_CHEM.
          !
          ! ALSO NOTE: When comparing History output to bpch output,
          ! you must use an updating frequency equal to the dynamic
          ! timestep so that the drydep fluxes due to mixing will
          ! be equivalent w/ the bpch output.  It is also recommended to
          ! turn off chemistry so as to be able to compare the drydep
          ! fluxes due to mixing in bpch vs. History as an "apples-to-
          ! apples" comparison.
          !
          !    -- Bob Yantosca (yantosca@seas.harvard.edu)
          !-----------------------------------------------------------------
          IF ( State_Diag%Archive_DryDepMix   .or.                           &
               State_Diag%Archive_DryDep    ) THEN
             S = State_Diag%Map_DryDepMix%id2slot(ND)
             IF ( S > 0 ) THEN
                State_Diag%DryDepMix(:,:,S) = Dflx(:,:,N)                    &
                                            * 1.0e-4_fp                      &
                                            * ( AVO / MW_kg  )
             ENDIF
          ENDIF

          !-----------------------------------------------------------------
          ! If Soil NOx is turned on, then call SOIL_DRYDEP to
          ! archive dry deposition fluxes for nitrogen species
          ! (SOIL_DRYDEP will exit if it can't find a match.
          !-----------------------------------------------------------------
	  IF ( Input_Opt%LSOILNOX ) THEN
             tmpFlx = 0.0_fp
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX
                tmpFlx = dflx(I,J,N)                                         &
	               / MW_kg                                               &
                       * AVO           * 1.e-4_fp                            &
                       * GET_TS_CONV() / GET_TS_EMIS()

                CALL Soil_DryDep( I, J, 1, N, tmpFlx, State_Chm )
             ENDDO
             ENDDO
	  ENDIF

          ! Free species database pointer
          ThisSpc => NULL()
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !=======================================================================
    ! Unit conversion #2: Convert back to the original units
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid,                &
                            State_Met, OrigUnit,  RC                        )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "Convert_Spc_Units" (from v/v dry)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Compute_Sflx_For_Vdiff
!EOC
END MODULE HCO_Utilities_GC_Mod
