!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: passive_species_mod.F90
!
! !DESCRIPTION: Module passive\_species\_mod.F90 contains variables and routines
!  for using passive species in GEOS-Chem. Passive species are species that are 
!  passively transported by GEOS-Chem, with a simple first order loss rate 
!  being applied to each species. The number of passive species, corresponding 
!  species properties as well as loss rates and default initial concentrations 
!  can be specified by the user via the GEOS-Chem input file (input.geos).
! 
!  The passive species module is designed to work in combination with any
!  existing GEOS-Chem simulation type, even though it has only been tested
!  with the Radon simulation at this point.
!
! !REMARKS:  
!  Passive species are defined in input.geos in the PASSIVE SPECIES MENU. For
!  instance, to use the Rn-Pb-Be simulation with two passive species (PASV1
!  and PASV2) with atmospheric lifetimes of 1 hour and 3.8 days, respectively, 
!  add the following entries to input.geos:
!
!  %%% PASSIVE SPECIES MENU %%%:
!  Number of passive spec. : 2
!  Passive species #1      : PASV1 50.0 3600.0 1.0e-12
!  Passive species #2      : PASV2 250.0 328320.0 1.0e-09
!
!  For each passive species you must define the molecular weight (g/mol),
!  atmospheric lifetime (s), and default initial atmospheric concentration
!  (v/v). The default initial concentration is only used if the passive species
!  are not included in the restart file.
!
!  There must be a matching entry in the advected species menu for every
!  passive species defined in the passive species menu:
!
!  %%% ADVECTED SPECIES MENU %%%:
!  Type of simulation      : 1
!  Number of Advected Spec.: 5
!  Species Entries ------->: Name
!  Species #1              : Rn
!  Species #2              : Pb
!  Species #3              : Be7
!  Species #4              : PASV1
!  Species #5              : PASV2
! 
!  In this example, species 1-3 are the default species for this simulation type
!  while species 4-5 are the user-specific passive species.
!
!  As for regular GEOS-Chem species, emissions can be assigned to passive 
!  species via the HEMCO configuration file. For example, to assign a uniform 
!  flux of 1e-3 and 1e-9 kg/m2/s to PASV1 and PASV2, respectively, add the
!  following line to section base emissions of your HEMCO_Config.rc:
!
!  0 PASV1_Flux 1.0e-3  - - - xy kg/m2/s PASV1 - 1 1
!  0 PASV2_Flux 1.0e-9  - - - xy kg/m2/s PASV2 - 1 1
!
!  There is no default diagnostics for these emissions but you can easily
!  define a HEMCO diagnostics for each passive species by creating a
!  corresponding entry in the HEMCO diagnostics file (HEMCO_diagn.rc):
!
!  # Name       Spec  ExtNr Cat Hier Dim Unit
!  PASV1_TOTAL  PASV1 -1    -1  -1   2   kg/m2/s
!  PASV2_TOTAL  PASV2 -1    -1  -1   2   kg/m2/s
!
!  To activate these diagnostics you need to link to the HEMCO diagnostics file
!  in HEMCO_Config.rc, as well as define the diagnostics output frequency (in
!  the settings section of the HEMCO configuration file). For example to write
!  the diagnostics at 30-minute intervals:
!
!  DiagnFile:                 HEMCO_Diagn.rc
!  DiagnFreq:                 00000000 003000
!\\
!\\
! !INTERFACE:
!
MODULE PASSIVE_SPECIES_MOD 
!
! !USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: INIT_PASSIVE_SPECIES
  PUBLIC :: ADD_PASSIVE_SPECIES
  PUBLIC :: PASSIVE_SPECIES_GETRATE
  PUBLIC :: PASSIVE_SPECIES_INQUIRE
  PUBLIC :: CLEANUP_PASSIVE_SPECIES
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !REVISION HISTORY:
!  04 Sep 2015 - C. Keller   - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  INTEGER,                        PUBLIC :: NPASSIVE = 0
  CHARACTER(LEN=63), ALLOCATABLE, PUBLIC :: PASSIVE_NAME    (:)
  INTEGER,  ALLOCATABLE                  :: PASSIVE_ID      (:)
  REAL(fp), ALLOCATABLE                  :: PASSIVE_MW      (:)
  REAL(fp), ALLOCATABLE                  :: PASSIVE_TAU     (:)
  REAL(fp), ALLOCATABLE                  :: PASSIVE_INITCONC(:)

  !============================================================================
  ! MODULE ROUTINES -- follow below the "CONTAINS" statement
  !============================================================================
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Passive_Species
!
! !DESCRIPTION: Subroutine INIT\_PASSIVE\_SPECIES initializes the passive
! species arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_PASSIVE_SPECIES( am_I_Root, nSpc, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    INTEGER,          INTENT(IN   )  :: nSpc       ! # of passive species 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER         :: AS

    !=================================================================
    ! INIT_PASSIVE_SPECIES begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS
  
    ! Set module variable 
    NPASSIVE = nSpc
    
    ! Allocate arrays
    IF ( NPASSIVE > 0 ) THEN
       ALLOCATE( PASSIVE_ID(NPASSIVE),       PASSIVE_TAU(NPASSIVE), &
                 PASSIVE_INITCONC(NPASSIVE), PASSIVE_MW(NPASSIVE),  &
                 PASSIVE_NAME(NPASSIVE),     STAT=AS )
       IF ( AS /= 0 ) THEN
          WRITE(*,*) 'Cannot allocate passive species arrays'
          RC = GC_FAILURE
          RETURN
       ENDIF

       PASSIVE_ID(:)       = -999
       PASSIVE_NAME(:)     = ''
       PASSIVE_MW(:)       = 0.0_fp
       PASSIVE_TAU(:)      = 0.0_fp
       PASSIVE_INITCONC(:) = 0.0_fp
    ENDIF

  END SUBROUTINE INIT_PASSIVE_SPECIES
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Add_Passive_Species
!
! !DESCRIPTION: Subroutine ADD\_PASSIVE\_SPECIES registers a passive species 
!  based on the passed input arguments.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ADD_PASSIVE_SPECIES( am_I_Root,   SpcName, SpcTau, &
                                  SpcInitConc, SpcMW,   RC ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root   ! root CPU?
    CHARACTER(LEN=*), INTENT(IN   )  :: SpcName     ! Species name
    REAL(fp),         INTENT(IN   )  :: SpcTau      ! Species lifetime (s) 
    REAL(fp),         INTENT(IN   )  :: SpcInitConc ! Species default init conc (v/v) 
    REAL(fp),         INTENT(IN   )  :: SpcMW       ! Species molec. weight (g/mol) 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, IDX, SPCID
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = 'ADD_PASSIVE_SPECIES (passive_species_mod.F90)'

    !=================================================================
    ! ADD_PASSIVE_SPECIES begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Error check: cannot define passive species if # of passive species is 0.
    IF ( NPASSIVE <= 0 ) THEN
       WRITE(*,*) 'Cannot add passive species ', TRIM(SpcName), &
             ': # of passive species is smaller than 1!'
       RC = GC_FAILURE
       RETURN 
    ENDIF 

    ! Find empty slot
    IDX = -1
    DO I = 1, NPASSIVE
       IF ( PASSIVE_ID(I) < 0 ) THEN
          IDX = I
          EXIT
       ENDIF
    ENDDO 
    
    IF ( IDX <= 0 ) THEN
       WRITE(*,*) 'Cannot add passive species ', TRIM(SpcName), &
                    ': all ', NPASSIVE, ' species slots are already being used.'
       RC = GC_FAILURE
       RETURN 
    ENDIF 

    ! Register species
    PASSIVE_ID(IDX)       = IDX
    PASSIVE_NAME(IDX)     = TRIM(SpcName)
    PASSIVE_TAU(IDX)      = SpcTau 
    PASSIVE_INITCONC(IDX) = SpcInitConc 
    PASSIVE_MW(IDX)       = SpcMW

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE( 6, '(a)' ) 'Added passive species: '
       WRITE( 6, 110   ) ' - Species name                 : ', PASSIVE_NAME(IDX) 
       WRITE( 6, 120   ) ' - Molec. weight [g/mol]       : ', PASSIVE_MW(IDX)
       WRITE( 6, 120   ) ' - Lifetime [s]                : ', PASSIVE_TAU(IDX)
       WRITE( 6, 130   ) ' - Default concentration [v/v] : ', PASSIVE_INITCONC(IDX)
    ENDIF

110 FORMAT( A, A )
120 FORMAT( A, F10.2  )
130 FORMAT( A, ES10.2 )
 
    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE ADD_PASSIVE_SPECIES
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Passive_Species_GetRate
!
! !DESCRIPTION: Subroutine PASSIVE\_SPECIES\_GETRATE returns the unitless decay
!  rate for the given species and chemistry time step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PASSIVE_SPECIES_GETRATE( am_I_Root, SpcName, DT, Rate, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    CHARACTER(LEN=*), INTENT(IN   )  :: SpcName    ! Passive species name 
    REAL(fp),         INTENT(IN   )  :: DT         ! Time step in s
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),         INTENT(  OUT)  :: Rate       ! Decay rate (unitless 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: N, ID
    REAL(fp)            :: Decay
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = "PASSIVE_SPECIES_GETRATE (PASSIVE_SPECIES_MOD.F90)"

    REAL(fp), PARAMETER :: ln2 = 0.693147181E+00_fp

    !=================================================================
    ! PASSIVE_SPECIES_GETRATE begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Find the index
    ID = -1
    DO N = 1, NPASSIVE

       IF ( TRIM(SpcName) == TRIM(PASSIVE_NAME(N)) ) THEN
          ID = N
          EXIT
       ENDIF
    ENDDO

    ! Error check
    IF ( ID <= 0 ) THEN
       WRITE(*,*) 'This is not a passive species ', TRIM(SpcName)
       RC = GC_FAILURE
       RETURN 
    ENDIF 

    ! No loss needed if tau is zero or negative
    IF ( PASSIVE_TAU(N) <= 0.0_fp ) THEN
       Rate = 1.0

    ! Calculate decay rate (unitless)
    ELSE 
       Decay    = ln2 / PASSIVE_TAU(N)
       Rate     = EXP( - DT * Decay )

    ENDIF

    ! Return w/ success 
    RC = GC_SUCCESS
 
  END SUBROUTINE PASSIVE_SPECIES_GETRATE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Passive_Species_Inquire
!
! !DESCRIPTION: Function PASSIVE\_SPECIES\_INQUIRE is a wrapper routine to 
!  inquire information about a passive species. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PASSIVE_SPECIES_INQUIRE( SpcName, IsPassive, MW, InitConc ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),   INTENT(IN   )  :: SpcName    ! GC species name 
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,  OPTIONAL, INTENT(  OUT)  :: IsPassive  ! Is SpcID a passive spec.?
    REAL(fp), OPTIONAL, INTENT(  OUT)  :: MW         ! Molecular weight (g/mol) 
    REAL(fp), OPTIONAL, INTENT(  OUT)  :: InitConc   ! Initial conc. (v/v)
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!  09 Mar 2017 - C. Keller    - Bug fix: initialize molw to default value of 0.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER    :: N, ID
    LOGICAL    :: IsPass
    REAL(fp)   :: InConc, molw

    !=================================================================
    ! PASSIVE_SPECIES_INQUIRE begins here!
    !=================================================================

    ! Init
    IsPass = .FALSE.
    InConc = 0.0_fp 
    molw   = 0.0_fp 

    ! Nothing to do if no passive species defined
    IF ( NPASSIVE > 0 ) THEN 

       ! Loop over all passive species
       DO N = 1, NPASSIVE
          IF ( TRIM(SpcName) == TRIM(PASSIVE_NAME(N)) ) THEN
             IsPass = .TRUE.
             InConc = PASSIVE_INITCONC(N)
             molw   = PASSIVE_MW(N)
             EXIT
          ENDIF
       ENDDO
    ENDIF

    ! Pass to output
    IF ( PRESENT(IsPassive) ) IsPassive = IsPass
    IF ( PRESENT(InitConc ) ) InitConc  = InConc
    IF ( PRESENT(MW       ) ) MW        = molw    

  END SUBROUTINE PASSIVE_SPECIES_INQUIRE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Passive_Species
!
! !DESCRIPTION: Subroutine CLEANUP\_PASSIVE\_SPECIES finalizes the passive
! species arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_PASSIVE_SPECIES( am_I_Root, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! CLEANUP_PASSIVE_SPECIES begins here!
    !=================================================================

    ! Deallocate arrays
    IF ( ALLOCATED( PASSIVE_ID       ) ) DEALLOCATE( PASSIVE_ID       )
    IF ( ALLOCATED( PASSIVE_TAU      ) ) DEALLOCATE( PASSIVE_TAU      )
    IF ( ALLOCATED( PASSIVE_INITCONC ) ) DEALLOCATE( PASSIVE_INITCONC )
    IF ( ALLOCATED( PASSIVE_MW       ) ) DEALLOCATE( PASSIVE_MW       )
    IF ( ALLOCATED( PASSIVE_NAME     ) ) DEALLOCATE( PASSIVE_NAME     )

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE CLEANUP_PASSIVE_SPECIES
!EOC
END MODULE PASSIVE_SPECIES_MOD 
