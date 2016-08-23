!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: flexchem_setup_mod
!
! !DESCRIPTION: Module FLEXCHEM\_SETUP\_MOD 
!\\
!\\
! !INTERFACE: 
!
MODULE FLEXCHEM_SETUP_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
  USE GIGC_Input_Opt_Mod, ONLY : OptInput
  USE GIGC_State_Chm_Mod, ONLY : ChmState
  USE GIGC_State_Met_Mod, ONLY : MetState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: CSPECTOKPP
  PUBLIC :: INIT_FLEXCHEM
  PUBLIC :: HSAVE_KPP
!    
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  15 Jun 2016 - M. Sulprizio- Remove STTTOCSPEC mapping array. Species and
!                              tracers have a 1:1 mapping currently so mapping
!                              is not required
!  18 Jul 2016 - M. Sulprizio- Remove FAMILIES_KLUDGE routine. Family tracers
!                              have been eliminated.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Index arrays
  INTEGER,   ALLOCATABLE :: CSPECTOKPP(:)

  ! Save the H value for the Rosenbrock solver
  REAL(fp),  ALLOCATABLE :: HSAVE_KPP(:,:,:) 

CONTAINS
!EOC  
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_FLEXCHEM
!
! !DESCRIPTION: Subroutine INIT\_FLEXCHEM is used to allocate
!               arrays for the KPP solver. 
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE INIT_FLEXCHEM( Input_Opt, State_Met, State_Chm, am_I_Root, RC)  
!
! !USES:
!
    USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
    USE CMN_SIZE_MOD
    USE PHYSCONSTANTS,      ONLY : BOLTZ
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_DIAGN_MOD,      ONLY : Diagn_Create
    USE HCO_ERROR_MOD
    USE GIGC_ErrCode_Mod
    USE GIGC_State_Chm_Mod, ONLY : IND_
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE gckpp_Global,       ONLY : NSPEC, NREACT
    USE gckpp_Monitor,      ONLY : SPC_NAMES, EQN_NAMES
    USE PRECISION_MOD
    USE Species_Mod,        ONLY : Species
    USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
!
! !INPUT PARAMETERS:
!    
    TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
    TYPE(MetState), INTENT(IN)    :: State_Met ! Met State object
    LOGICAL,        INTENT(IN)    :: am_I_Root ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!    
    INTEGER,        INTENT(OUT)   :: RC        ! Success or failure?
!    
! !REVISION HISTORY:
!  16 Sep 2009 - P. Le Sager - init
!  09 Dec 2009 - C. Carouge  - R_KPP and CSPEC_FOR_KPP are not used anymore
!  02 Aug 2012 - R. Yantosca - Now use am_I_Root to print on root CPU
!  22 Dec 2015 - M. Sulprizio- Use State_Met%AIRNUMDEN to convert initial
!                              species concentrations from v/v to molec/cm3
!  29 Jan 2016 - M. Sulprizio- Add calls to Register_Tracer and Register_Species
!                              to populate Tracer_Name, Tracer_Id, Species_Name,
!                              and Species_ID fields in State_Chm
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with SpcInfo%Name from species database
!  06 Jun 2016 - M. Sulprizio- Replace Get_Indx with Spc_GetIndx to use the
!                              fast-species lookup from the species database
!  06 Jun 2016 - M. Sulprizio- Remove calls to Register_Tracer and
!                              Register_Species; these routines were made
!                              obsolete by the species database
!  14 Jun 2016 - M. Sulprizio- Replace Spc_GetIndx with Ind_ (M. Long)
!  25 Jul 2016 - E. Lundgren - Add check that species was not in restart file
!                              prior to v/v -> molec/cm3 conversion
!  02 Aug 2016 - E. Lundgren - Move unit conversion of species background
!                              values to restart_mod
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: AS, D, N, N1, FOUND
    INTEGER            :: cId, Collection
    CHARACTER(LEN=15)  :: OutOper,  WriteFreq
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'INIT_FLEXCHEM (flexchem_setup_mod.F90)' 

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! INIT_FLEXCHEM begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Initialize pointer
    SpcInfo => NULL()

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = 'Instantaneous'
    WriteFreq  = 'Always'

    IF ( am_I_Root ) WRITE( 6, 100 )
100 FORMAT( '     - INIT_FLEXCHEM: Allocating arrays for FLEX_CHEMISTRY' )

    ALLOCATE( CSPECTOKPP(NSPEC), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CSPECTOKPP = 0

    ALLOCATE( HSAVE_KPP( IIPAR, JJPAR, LLCHEM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    HSAVE_KPP = 0.d0

    ! Loop over GEOS-Chem species
    DO N = 1, State_Chm%nSpecies

       ! Initialize
       FOUND=0

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Loop over KPP species
       DO N1 = 1, NSPEC

          ! Create vector to map
          ! State_Chm%Species order defined in species_database_mod.F90
          ! to KPP order defined in gckpp_Parameters.F90
          IF ( ADJUSTL( TRIM( SPC_NAMES(N1) ) ) == &
               ADJUSTL( TRIM( SpcInfo%Name  ) ) ) THEN
             CSPECTOKPP(N1)=N
             FOUND=1
             EXIT
          ENDIF

       ENDDO

       ! Print info to log
       IF (FOUND .NE. 1) THEN
          WRITE (6,'(a8,a17)') TRIM( SpcInfo%Name ),   &
             ' NOT found in KPP'
       ENDIF
       IF (FOUND .EQ. 1) THEN
          WRITE (6,'(a8,a29,I4)') TRIM( SpcInfo%Name ), &
             ' was found in KPP with index ', N1
       ENDIF

       ! Free pointer
       SpcInfo => NULL()

    ENDDO

    IF (am_I_Root) THEN
       write(*,'(a)') ' KPP Reaction Reference '
       DO D = 1,NREACT
          write(*,'(i,a3,a85)') D,' | ',EQN_NAMES(D)
       END DO
    ENDIF

    ! Return to calling program
    RETURN

  END SUBROUTINE INIT_FLEXCHEM
!EOC
END MODULE FLEXCHEM_SETUP_MOD

