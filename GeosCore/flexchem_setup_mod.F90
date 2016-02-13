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

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: STTTOCSPEC, CSPECTOKPP
  PUBLIC :: INIT_FLEXCHEM
  PUBLIC :: HSAVE_KPP
!    
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  INTEGER,   ALLOCATABLE :: STTTOCSPEC(:)
  INTEGER,   ALLOCATABLE :: CSPECTOKPP(:)
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

    USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
    USE CMN_SIZE_MOD
    USE COMODE_LOOP_MOD,    ONLY : QBKGAS, NGAS, NAMEGAS, BK, SMAL2
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_DIAGN_MOD,      ONLY : Diagn_Create
    USE HCO_ERROR_MOD
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState, Register_Tracer, Register_Species
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE gckpp_Global,       ONLY : NSPEC, NREACT
    USE gckpp_Monitor,      ONLY : SPC_NAMES, EQN_NAMES
    USE PRECISION_MOD
    USE RESTART_MOD,        ONLY : READ_CSPEC_FILE 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: AS, D, I, J, L, N, N1, FOUND
    INTEGER            :: cId, Collection
    CHARACTER(LEN=15)  :: OutOper,  WriteFreq
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'INIT_FLEXCHEM (flexchem_setup_mod.F90)' 
    LOGICAL            :: IT_EXISTS
    REAL(fp)           :: AIRNUMDEN

    ! Assume success
    RC = 0

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND45_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND45_OUTPUT_FREQ

    IF ( am_I_Root ) WRITE( 6, 100 )
    ALLOCATE( STTTOCSPEC(Input_Opt%N_TRACERS), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    STTTOCSPEC = 0

    ALLOCATE( CSPECTOKPP(NSPEC), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    CSPECTOKPP = 0

    ALLOCATE( HSAVE_KPP( IIPAR, JJPAR, LLCHEM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    HSAVE_KPP = 0.d0

    ! Loop over GEOS-Chem tracers
    DO N=1,Input_Opt%N_TRACERS

       ! Register tracers in the State_Chm object
       ! NOTE: This is here to populate the Trac_Id and Trac_Name fields in
       ! State_Chm. This can be removed when we fully utilize the species
       ! database (mps, 1/19/16)
       CALL Register_Tracer( Name      = Input_Opt%TRACER_NAME(N),  &
                             ID        = Input_Opt%ID_TRACER(N),    &
                             State_Chm = State_Chm,                 &
                             Status    = RC                    )
       IF ( am_I_Root .and. RC > 0 ) THEN
          WRITE( 6, 200 ) Input_Opt%TRACER_NAME(N), RC
 200      FORMAT( 'Registered Tracer : ', a14, i5 )
       ENDIF

       ! mje Find link between tracers and species 
       ! mje This is easier as we don't advect families any more
       FOUND=0
       DO N1=1,Input_Opt%N_SPECIES
          IF ( ADJUSTL(TRIM(NAMEGAS(N1))) == &
               ADJUSTL(TRIM(Input_Opt%TRACER_NAME(N)))) THEN
             STTTOCSPEC(N)=N1
             FOUND=1
          ENDIF
       ENDDO
       IF (FOUND .NE. 1) THEN
          WRITE (6,'(a8,a17)') TRIM(Input_Opt%TRACER_NAME(N)), &
             ' is not a species'
       ENDIF

    ENDDO

    ! Loop over GEOS-Chem species
    DO N=1,Input_Opt%N_SPECIES

       ! Register species (active + inactive) in the State_Chm object
       ! NOTE: This is here to populate the Spec_Id and Spec_Name fields in
       ! State_Chm. This can be removed when we fully utilize the species
       ! database (mps, 1/19/16)
       CALL Register_Species( NAME      = NAMEGAS(N),  &
                              ID        = N,           &
                              State_Chm = State_Chm,   &
                              Status    = RC        )
       IF ( am_I_Root .and. RC > 0 ) THEN
          WRITE( 6, 205 ) NAMEGAS(N), RC
205       FORMAT( 'Registered Species : ', a14, i5 )
       ENDIF

       ! MSL - Create vector to map species in State_Chm%Species order, defined
       !       in READCHEM, to KPP's order, as seen in gckpp_Parameters.F90
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
          WRITE (6,'(a8,a17)') TRIM(NAMEGAS(N)),   &
             ' NOT found in KPP'
       ENDIF
       IF (FOUND .EQ. 1) THEN
          WRITE (6,'(a8,a29,I4)') TRIM(NAMEGAS(N)), &
             ' was found in KPP with index ', N1
       ENDIF
    ENDDO

!------------------------------------------------------------------------------
! Initialize Initial Species Concentrations
! -- If read from a species restart file, otherwise
! -- Set species to their assigned background concentrations
!
! -- For the time being, we still read in globchem.dat thru READCHEM,
!    but this should be squashed for simplicity and stream-lining.
!    MSL - 12/2/2015
!------------------------------------------------------------------------------

    ! Read CSPEC restart file if option is turned on in input.geos
    IF ( Input_Opt%LSVCSPEC ) THEN

       ! Read the species restart file.  If not found, then 
       ! return IT_EXISTS = .FALSE.
       CALL READ_CSPEC_FILE( am_I_Root,  Input_Opt, GET_NYMD(),  &
                             GET_NHMS(), IT_EXISTS, State_Chm,  RC  )

    ENDIF

    !--------------------------------------------------------------------------
    ! If no species restart file is found or if not using species restart file
    ! then initialize species with default values from globchem.dat
    !--------------------------------------------------------------------------
    IF ( ( Input_Opt%LSVCSPEC .and. .not. IT_EXISTS ) .or. &
         ( .not. Input_Opt%LSVCSPEC ) ) THEN 

       IF ( am_I_Root ) THEN
          WRITE(6,*) '    - FLEXCHEM_SETUP: Not using species restart. ', &
             'Initialize species to background values from globchem.dat.'
       ENDIF

       ! Loop over species
       DO N=1,Input_Opt%N_SPECIES
          FOUND=0
          IF ( (ADJUSTL(TRIM(NAMEGAS(N)))) .eq. "" ) CYCLE
          DO N1=1,NSPEC

             IF ( ADJUSTL(TRIM(SPC_NAMES(N1))) == &
                  ADJUSTL(TRIM(NAMEGAS(N)   )) ) THEN

                FOUND=1

                ! Loop over grid boxes
                DO L = 1, LLCHEM 
                DO J = 1, JJPAR
                DO I = 1, IIPAR

                   ! Set grid box dry air density [molec/cm3]
                   AIRNUMDEN = State_Met%AIRNUMDEN(I,J,L)

                   !========================================================
                   ! For methanol (MOH), now use different initial background
                   ! concentrations for different regions of the atmosphere:
                   !
                   ! (a) 2.0 ppbv MOH -- continental boundary layer
                   ! (b) 0.9 ppbv MOH -- marine boundary layer
                   ! (c) 0.6 ppbv MOH -- free troposphere
                   !
                   ! The concentrations listed above are from Heikes et al,
                   ! "Atmospheric methanol budget and ocean implication",
                   ! _Global Biogeochem. Cycles_, submitted, 2002.  These
                   ! represent the best estimates for the methanol conc.'s
                   ! in the troposphere based on various measurements.
                   !
                   ! MOH is an inactive chemical species in GEOS-CHEM, so
                   ! these initial concentrations will never change. However,
                   ! MOH acts as a sink for OH, and therefore will affect
                   ! both the OH concentration and the methylchloroform
                   ! lifetime.
                   !
                   ! We specify the MOH concentration as ppbv, but then we
                   ! need to multiply by PRESS3(JLOOP) / ( T3(JLOOP) * BK )
                   ! in order to convert to [molec/cm3].  (bdf, bmy, 2/22/02)
                   !========================================================
                   IF ( NAMEGAS(N) == 'MOH' ) THEN

                      !------------------------------
                      ! Test for altitude
                      ! L < 9 is always in the trop.
                      !------------------------------
                      IF ( L <= 9 ) THEN

                         !---------------------------
                         ! Test for ocean/land boxes
                         !---------------------------
                         IF ( State_Met%FRCLND(I,J) >= 0.5 ) THEN

                            ! Continental boundary layer: 2 ppbv MOH
                            State_Chm%Species(I,J,L,N) = 2.000e-9_fp * &
                                                         AIRNUMDEN

                            ! Make sure MOH conc. is not negative
                            ! (SMAL2 = 1d-99)
                            State_Chm%Species(I,J,L,N) = &
                                    MAX(State_Chm%Species(I,J,L,N),SMAL2)

                         ELSE

                            ! Marine boundary layer: 0.9 ppbv MOH
                            State_Chm%Species(I,J,L,N) = 0.900e-9_fp * &
                                                         AIRNUMDEN

                            ! Make sure MOH conc. is not negative
                            ! (SMAL2 = 1d-99)
                            State_Chm%Species(I,J,L,N) = &
                                    MAX(State_Chm%Species(I,J,L,N),SMAL2)
                         ENDIF

                      ELSE

                         !---------------------------
                         ! Test for troposphere
                         !---------------------------
                         IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) THEN

                            ! Free troposphere: 0.6 ppbv MOH
                            State_Chm%Species(I,J,L,N) = 0.600e-9_fp * &
                                                         AIRNUMDEN

                            ! Make sure MOH conc. is not negative
                            ! (SMAL2 = 1d-99)
                            State_Chm%Species(I,J,L,N) = &
                                    MAX(State_Chm%Species(I,J,L,N),SMAL2)

                         ELSE

                            ! Strat/mesosphere: set MOH conc. to
                            ! SMAL2 = 1d-99
                            State_Chm%Species(I,J,L,N) = SMAL2

                         ENDIF

                      ENDIF

                      ! Debug
                      IF (I .eq. 10 .and. J .eq. 10 .and. L .eq. 10) THEN
                         WRITE(6,*) 'Initialized MOH by region'
                      ENDIF

                   !========================================================
                   ! For species other than MOH
                   ! Copy default background conc. from globchem.dat
                   !========================================================
                   ELSE

                      ! Convert from v/v to molec/cm3
                      State_Chm%Species(I,J,L,N) = QBKGAS(N) * AIRNUMDEN

                      ! Make sure concentration is not negative
                      ! (SMAL2 = 1d-99)
                      State_Chm%Species(I,J,L,N) = &
                              MAX(State_Chm%Species(I,J,L,N),SMAL2)

                   ENDIF

                ENDDO
                ENDDO
                ENDDO

             ENDIF
             IF (FOUND .EQ. 1) THEN
                WRITE (6,'(a8,a16,e10.3,a4)') TRIM(SPC_NAMES(N1)), &
                               ' initialized to ', QBKGAS(N), ' v/v'
                EXIT
             ENDIF

          ENDDO

          IF (FOUND .NE. 1) &
               WRITE (6,'(a8,a16)') TRIM(NAMEGAS(N)),' NOT initialized'

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
!EOC
END MODULE FLEXCHEM_SETUP_MOD

