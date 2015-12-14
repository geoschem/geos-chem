!$ID$
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE:
!
! !DESCRIPTION:
!              
!\\
!\\
! !INTERFACE: 
!
MODULE FLEXCHEM_SETUP_MOD
! 
! !USES:
!  
  
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: STTTOCSPEC, CSPECTOKPP
  PUBLIC :: INIT_FLEXCHEM
  PUBLIC :: HSAVE_KPP
  
  INTEGER,   ALLOCATABLE :: STTTOCSPEC(:), CSPECTOKPP(:)
  REAL(fp),  ALLOCATABLE :: HSAVE_KPP(:,:,:) 
!    
! !REVISION HISTORY:
!   
!EOP
!------------------------------------------------------------------------------
    
CONTAINS
  
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_GCKPP_COMODE
!
! !DESCRIPTION: Subroutine INIT\_GCKPP\_COMODE is used to allocate
!               arrays for the KPP solver. 
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE INIT_FLEXCHEM( Input_Opt, State_Chm, am_I_Root, RC)  

    USE CMN_SIZE_MOD
    USE COMODE_LOOP_MOD,    ONLY : QBKGAS, NGAS, NAMEGAS
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_DIAGN_MOD,      ONLY : Diagn_Create
    USE HCO_ERROR_MOD
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE gckpp_Global,       ONLY : NSPEC, NREACT
    USE gckpp_Monitor,      ONLY : SPC_NAMES, EQN_NAMES
    USE RESTART_MOD,        ONLY : READ_CSPEC_FILE 
    USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
!
! !INPUT PARAMETERS:
!    
    TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
    LOGICAL,        INTENT(IN)    :: am_I_Root ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!    
    INTEGER,        INTENT(OUT):: RC
!    
! !REVISION HISTORY:
!  16 Sep 2009 - P. Le Sager - init
!  09 Dec 2009 - C. Carouge  - R_KPP and CSPEC_FOR_KPP are not used anymore
!  02 Aug 2012 - R. Yantosca - Now use am_I_Root to print on root CPU
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS, D, I, N, N1, FOUND
    INTEGER            :: cId, Collection
    CHARACTER(LEN=15)  :: OutOper,  WriteFreq
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'INIT_FLEXCHEM (flexchem_setup_mod.F90)' 
    LOGICAL            :: IT_EXISTS

    RC=0

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND45_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND45_OUTPUT_FREQ

    IF ( am_I_Root ) WRITE( 6, 100 )
    ALLOCATE( STTTOCSPEC(Input_Opt%N_TRACERS), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    IF ( am_I_Root ) WRITE( 6, 100 )
    ALLOCATE( CSPECTOKPP(NSPEC), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    IF ( am_I_Root ) WRITE( 6, 100 )
    ALLOCATE( HSAVE_KPP( IIPAR, JJPAR, LLCHEM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    HSAVE_KPP = 0.d0

!     mje Find link between tracers and species 
!     mje This is easier as we don't advect families any more
    DO N=1,Input_Opt%N_TRACERS
       FOUND=0
       STTTOCSPEC(N)=0
       DO N1=1,NSPEC
          IF (ADJUSTL(TRIM(SPC_NAMES(N1))) == &
               ADJUSTL(TRIM(Input_Opt%TRACER_NAME(N)))) THEN
             STTTOCSPEC(N)=N1
             FOUND=1
          ENDIF
       ENDDO
       IF (FOUND .NE. 1) THEN
          WRITE (6,*) TRIM(Input_Opt%TRACER_NAME(N)),' Not a species'
       ENDIF
    ENDDO

!   MSL - Create vector to map species in CSPEC order, defined
!         in READCHEM, to KPP's order, as seen in gckpp_Parameters.F90
    CSPECTOKPP=0
    DO N=1,Input_Opt%N_SPECIES
       FOUND=0
       DO N1=1,NSPEC
          IF (ADJUSTL(TRIM(SPC_NAMES(N1))) == &
               ADJUSTL(TRIM(NAMEGAS(N)))) THEN
             CSPECTOKPP(N1)=N
             FOUND=1
             EXIT
          ENDIF
       ENDDO
       IF (FOUND .NE. 1) THEN
          WRITE (6,*) TRIM(NAMEGAS(N)),' NOT FOUND IN KPP'
       ENDIF
       IF (FOUND .EQ. 1) THEN
          WRITE (6,*) TRIM(NAMEGAS(N)),' FOUND IN KPP WITH INDEX ', NSPEC
       ENDIF
    ENDDO

!------------------------------------------------------------------------------
! Initialize Initial Species Concentrations
! -- If read from a cspec restart file, otherwise
! -- Set species to their assigned background concentrations
!
! -- For the time being, we still read in globchem.dat thru READCHEM,
!    but this should be squashed for simplicity and stream-lining.
!    MSL - 12/2/2015
!------------------------------------------------------------------------------

    IF ( Input_Opt%LSVCSPEC ) THEN
       ! Read the CSPEC restart file.  If not found, then 
       ! return IT_EXISTS = .FALSE.
       CALL READ_CSPEC_FILE( am_I_Root,  Input_Opt, GET_NYMD(),  &
                             GET_NHMS(), IT_EXISTS, State_Chm,  RC  )
       
       IF ( .not. IT_EXISTS ) THEN 
          
          !-----------------------------------------------------------
          ! (2) Initialize CSPEC with default values from globchem.dat.
          !-----------------------------------------------------------
          IF ( am_I_Root ) THEN
             WRITE(6,*) '    - CHEMDR: CSPEC restart not found, use background values'
          ENDIF
          
          DO N=1,Input_Opt%N_SPECIES
             FOUND=0
             DO N1=1,NSPEC
                IF (ADJUSTL(TRIM(SPC_NAMES(N1))) == &
                     ADJUSTL(TRIM(NAMEGAS(N)))) THEN
                   FOUND=1
                   State_Chm%Species(:,:,:,N) = QBKGAS(N) ! Is this mol/mol?
                ENDIF
                IF (FOUND .EQ. 1) THEN
                   WRITE (6,*) TRIM(SPC_NAMES(N1)),' initialized to', QBKGAS(N), &
                        TRIM(NAMEGAS(N))
                   EXIT
                ENDIF
             ENDDO
             IF (FOUND .NE. 1) &
                WRITE (6,*) TRIM(SPC_NAMES(N1)),' NOT initialized'
          ENDDO
       ENDIF
    ELSE ! No CSPEC file specified in input.geod. Set background values
       DO N=1,Input_Opt%N_SPECIES
          FOUND=0
          IF ( (ADJUSTL(TRIM(NAMEGAS(N)))) .eq. "" ) CYCLE
          DO N1=1,NSPEC
             IF (ADJUSTL(TRIM(SPC_NAMES(N1))) == &
                  ADJUSTL(TRIM(NAMEGAS(N)))) THEN
                FOUND=1
                State_Chm%Species(:,:,:,N) = QBKGAS(N) ! Is this mol/mol?
             ENDIF
             IF (FOUND .EQ. 1) THEN
                WRITE (6,'(2a,e8.1)') TRIM(SPC_NAMES(N1)),' initialized to ', QBKGAS(N)
                EXIT
             ENDIF
          ENDDO
          IF (FOUND .NE. 1) &
               WRITE (6,*) TRIM(NAMEGAS(N)),' NOT initialized'
       ENDDO
    ENDIF
    
!------------------------------------------------------------------------------
! Initialize Diagnostics Per the new diagnostic package | M. Long 1-21-15
!------------------------------------------------------------------------------
    
      DO D = 1, Input_Opt%N_TRACERS
         
            !-----------------------------------------------------------
            ! Create containers for tracer concentration
            !-----------------------------------------------------------

            ! Diagnostic name
            DiagnName = TRIM( Input_Opt%TRACER_NAME(D) )

            ! Create container
            CALL Diagn_Create( am_I_Root,                     &
                 Col       = Collection,        & 
                 cID       = cID,               &
                 cName     = TRIM( DiagnName ), &
                 AutoFill  = 0,                 &
                 ExtNr     = -1,                &
                 Cat       = -1,                &
                 Hier      = -1,                &
                 HcoID     = -1,                &
                 SpaceDim  =  3,                &
                 LevIDx    = -1,                &
                 OutUnit   = 'cm-3',            &
                 OutOper   = TRIM( OutOper   ), &
                 OkIfExist = .TRUE.,            &
                 RC        = RC )
            
            IF ( RC /= HCO_SUCCESS ) THEN
               write(*,*) 'Collection: ', collection, RC, trim(writefreq),trim(outoper)
               MSG = ': Cannot create diagnostics: ' // TRIM(DiagnName)
               CALL ERROR_STOP( MSG, LOC ) 
            ENDIF
      ENDDO

      IF (am_I_Root) THEN
         write(*,'(a)') ' KPP Reaction Reference '
         DO D = 1,NREACT
            write(*,'(i,a3,a85)') D,' | ',EQN_NAMES(D)
         END DO
      ENDIF

100 FORMAT( '     - INIT_FLEXCHEM: Allocating arrays for FLEX_CHEMISTRY' )
    ! Return to calling program
    RETURN
END SUBROUTINE INIT_FLEXCHEM
END MODULE FLEXCHEM_SETUP_MOD

