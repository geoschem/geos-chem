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
  PUBLIC :: STTTOCSPEC, CSPECTOKPP
  PUBLIC :: INIT_FLEXCHEM
  PUBLIC :: HSAVE_KPP
  PUBLIC :: FAMILIES_KLUDGE
!    
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Index arrays
  INTEGER,   ALLOCATABLE :: STTTOCSPEC(:)
  INTEGER,   ALLOCATABLE :: CSPECTOKPP(:)

  ! Save the H value for the Rosenbrock solver
  REAL(fp),  ALLOCATABLE :: HSAVE_KPP(:,:,:) 

  ! Tracer and species indices for use in FAMILIES_KLUDGE
  INTEGER                :: T_MMN,   S_MVKN,    S_MACRN
  INTEGER                :: T_ISOPN, S_ISOPND,  S_ISOPNB
  INTEGER                :: T_CFCX,  S_CFC113,  S_CFC114,   S_CFC115
  INTEGER                :: T_HCFCX, S_HCFC123, S_HCFC141b, S_HCFC142b

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
    USE GIGC_State_Chm_Mod, ONLY : Register_Tracer
    USE GIGC_State_Chm_Mod, ONLY : Register_Species
    USE GIGC_State_Chm_Mod, ONLY : Get_Indx
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

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !% IMPORTANT:                                                  %
       !% N_SPECIES = NTSPEC(NCS) = active + inactive species         %
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                   AIRNUMDEN = State_Met%PMID_DRY(I,J,L) * 1000e+0_fp / &
                       ( State_Met%T(I,J,L) * BK ) !State_Met%AIRNUMDEN(I,J,L)

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
                WRITE (6,'(i4,a10,a16,e10.3,a4)') N, TRIM(SPC_NAMES(N1)), &
                               ' initialized to ', QBKGAS(N), ' v/v'
                EXIT
             ENDIF

          ENDDO

          IF (FOUND .NE. 1) &
               WRITE (6,'(i4,a10,a16)') N, TRIM(NAMEGAS(N)),' NOT initialized'

       ENDDO
    ENDIF

#if defined( DEVEL )
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
#endif

    IF (am_I_Root) THEN
       write(*,'(a)') ' KPP Reaction Reference '
       DO D = 1,NREACT
          write(*,'(i,a3,a85)') D,' | ',EQN_NAMES(D)
       END DO
    ENDIF

100 FORMAT( '     - INIT_FLEXCHEM: Allocating arrays for FLEX_CHEMISTRY' )

    !=====================================================================
    ! Save species indices during the Init phase to avoid
    ! repeated costly string matching operations
    !=====================================================================

    ! MMN family
    T_MMN      = Get_Indx( 'MMN',      State_Chm%Trac_ID, State_Chm%Trac_Name )
    S_MVKN     = Get_Indx( 'MVKN',     State_Chm%Spec_ID, State_Chm%Spec_Name )
    S_MACRN    = Get_Indx( 'MACRN',    State_Chm%Spec_ID, State_Chm%Spec_Name )

    ! ISOPN family
    T_ISOPN    = Get_Indx( 'ISOPN',    State_Chm%Trac_ID, State_Chm%Trac_Name )
    S_ISOPND   = Get_Indx( 'ISOPND',   State_Chm%Spec_ID, State_Chm%Spec_Name )
    S_ISOPNB   = Get_Indx( 'ISOPNB',   State_Chm%Spec_ID, State_Chm%Spec_Name )

    ! CFCX family
    T_CFCX     = Get_Indx( 'CFCX',     State_Chm%Trac_ID, State_Chm%Trac_Name )
    S_CFC113   = Get_Indx( 'CFC113',   State_Chm%Spec_ID, State_Chm%Spec_Name )
    S_CFC114   = Get_Indx( 'CFC114',   State_Chm%Spec_ID, State_Chm%Spec_Name )
    S_CFC115   = Get_Indx( 'CFC115',   State_Chm%Spec_ID, State_Chm%Spec_Name )

    ! HCFCX family
    T_HCFCX    = Get_Indx( 'HCFCX',    State_Chm%Trac_ID, State_Chm%Trac_Name )
    S_HCFC123  = Get_Indx( 'HCFC123',  State_Chm%Spec_ID, State_Chm%Spec_Name )
    S_HCFC141b = Get_Indx( 'HCFC141b', State_Chm%Spec_ID, State_Chm%Spec_Name )
    S_HCFC142b = Get_Indx( 'HCFC142b', State_Chm%Spec_ID, State_Chm%Spec_Name )

    ! Return to calling program
    RETURN

  END SUBROUTINE INIT_FLEXCHEM

  SUBROUTINE FAMILIES_KLUDGE(am_I_Root,STT,IO,SC,OPT,RC)
    USE CMN_SIZE_MOD,         ONLY : IIPAR, JJPAR, LLPAR
    REAL(kind=fp)  :: STT(:,:,:,:)
    TYPE(OptInput) :: IO ! Short-hand for Input_Opt
    TYPE(ChmState) :: SC ! Short-hand for State_Chem
    INTEGER        :: RC,OPT
    LOGICAL        :: am_I_Root

    INTEGER N,S1,S2,S3
    REAL(kind=fp)  ::  QSUM(IIPAR,JJPAR,LLPAR)
    REAL(kind=fp)  :: QTEMP(IIPAR,JJPAR,LLPAR)
    
    
      ! K -- L -- U -- D -- G -- E -- C -- O -- D -- E -- ! -- ! -- ! -- !
      ! THIS IS A TEMPORARY FIX AND NEEDS TO BE RESOLVED THROUGHOUT
      ! THE MODEL AS SOON AS IT APPEARS POSSIBLE TO DO SO!
      ! -- The following code ensures that the two remaining tracer
      !    families are dealt with appropriately.
      !    The Tracer MMN is a sum of the species MVKN and MACRN
      !    The Tracer ISOPN is a sum of ISOPND and ISOPNB
      !       The tracer restart file for the pre-flex GEOS-Chem, includes
      !    these families, but the removal of the routines
      !    "lump" and "partition" killed the code resposible for them.
      !       Here, we want to install a hard-code fix with the complete
      !    expectation that we will no longer use species families. Thus,
      !    when possible, the two remaining families need to be removed, 
      !    and this kludge disabled.
      !   M.S.L. - Jan., 5 2016
      !
      ! Prevent div-by-zero statements.  Renamed SUM to QSUM to avoid
      ! conflicts with the Fortran intrinsic function SUM.
      !  (bmy, 3/28/16)
      !

    IF (OPT .eq. 1) THEN
      ! Part 1: From Tracers to Species

      ! -- Process MMN family
      N  = T_MMN     ! MMN   species index
      S1 = S_MVKN    ! MVKN  species index
      S2 = S_MACRN   ! MACRN species index

      STT(:,:,:,N) = STT(:,:,:,N)*IO%TRACER_COEFF(N,1) ! A correction...
      ! ... to account for the division by TRACER_COEFF above.
      QSUM   = (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) &
             + (SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2))

      ! -- -- First, do MVKN
      WHERE( QSUM > 0.0_fp ) 
         QTEMP = (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) / QSUM
      ELSEWHERE
         QTEMP = 0.0_fp
      ENDWHERE
      SC%Species(:,:,:,S1) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )

      ! -- -- Then, do MACRN
      WHERE( QSUM > 0.0_fp ) 
         QTEMP = (SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2)) / QSUM
      ELSEWHERE
         QTEMP = 0.0_fp
      ENDWHERE
      SC%Species(:,:,:,S2) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )

      ! -- Process ISOPN family
      ! -- Get the tracer and species indices
      N  = T_ISOPN    ! ISOPN  tracer  index
      S1 = S_ISOPND   ! ISOPND species index
      S2 = S_ISOPNB   ! ISOPNB species index

      STT(:,:,:,N) = STT(:,:,:,N)*IO%TRACER_COEFF(N,1) ! A correction...
      ! ... to account for the division by TRACER_COEFF above.
      QSUM   = (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) &
             + (SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2))

      ! -- -- First, do ISOPND
      WHERE( QSUM > 0.0_fp ) 
         QTEMP = (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) / QSUM
      ELSEWHERE
         QTEMP = 0.0_fp
      ENDWHERE
      SC%Species(:,:,:,S1) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )

      ! -- -- Then, do ISOPNB
      WHERE( QSUM > 0.0_fp ) 
         QTEMP = (SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2)) / QSUM
      ELSEWHERE
         QTEMP = 0.0_fp
      ENDWHERE
      SC%Species(:,:,:,S2) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )

#if defined( UCX )
      IF ( IO%LUCX ) THEN
         ! -- Process CFCX family
         N  = T_CFCX     ! CFCX   tracer  index
         S1 = S_CFC113   ! CFC113 species index
         S2 = S_CFC114   ! CFC114 species index
         S3 = S_CFC115   ! CFC115 species index

         STT(:,:,:,N) = STT(:,:,:,N)*IO%TRACER_COEFF(N,1) ! A correction...
         ! ... to account for the division by TRACER_COEFF above.
         QSUM  = (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) &
               + (SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2)) &
               + (SC%Species(:,:,:,S3)*IO%TRACER_COEFF(N,3))

         ! -- -- First, do CFC113
         WHERE( QSUM > 0.0_fp )
            QTEMP = (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) / QSUM
         ELSEWHERE
            QTEMP = 0.0_fp
         ENDWHERE
         SC%Species(:,:,:,S1) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )

         ! -- -- Then, do CFC114
         WHERE( QSUM > 0.0_fp )
            QTEMP = (SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2)) / QSUM
         ELSEWHERE
            QTEMP = 0.0_fp
         ENDWHERE
         SC%Species(:,:,:,S2) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )

         ! -- -- Then, do CFC115
         WHERE( QSUM > 0.0_fp )
            QTEMP = (SC%Species(:,:,:,S3)*IO%TRACER_COEFF(N,3)) / QSUM
         ELSEWHERE
            QTEMP = 0.0_fp
         ENDWHERE
         SC%Species(:,:,:,S3) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )

         ! -- Process HCFCX family
         N  = T_HCFCX      ! HCFCX    tracer  index
         S1 = S_HCFC123    ! HCFC123  species index
         S2 = S_HCFC141b   ! HCFC141b species index
         S3 = S_HCFC142b   ! HCFC142b species index

         STT(:,:,:,N) = STT(:,:,:,N)*IO%TRACER_COEFF(N,1) ! A correction...
         ! ... to account for the division by TRACER_COEFF above.
         QSUM  = (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) &
               + (SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2)) &
               + (SC%Species(:,:,:,S3)*IO%TRACER_COEFF(N,3))

         ! -- -- First, do HCFC123
         WHERE( QSUM > 0.0_fp )
            QTEMP = (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) / QSUM
         ELSEWHERE
            QTEMP = 0.0_fp
         ENDWHERE
         SC%Species(:,:,:,S1) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )

         ! -- -- Then, do HCFC141b
         WHERE( QSUM > 0.0_fp )
            QTEMP = (SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2)) / QSUM
         ELSEWHERE
            QTEMP = 0.0_fp
         ENDWHERE
         SC%Species(:,:,:,S2) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )

         ! -- -- Then, do HCFC142b
         WHERE( QSUM > 0.0_fp )
            QTEMP = (SC%Species(:,:,:,S3)*IO%TRACER_COEFF(N,3)) / QSUM
         ELSEWHERE
            QTEMP = 0.0_fp
         ENDWHERE
         SC%Species(:,:,:,S3) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )
      ENDIF !IF UCX
#endif

      ! E -- N -- D -- O -- F -- K -- L -- U -- D -- G -- E -- -- P -- T -- 1
      ELSEIF (OPT .eq. 2) THEN
         ! K -- L -- U -- D -- G -- E -- -- P -- T -- 2
         ! Part 2: From Species to Tracers

         ! -- Process MMN family
         N  = T_MMN     ! MMN   species index
         S1 = S_MVKN    ! MVKN  species index
         S2 = S_MACRN   ! MACRN species index
         STT(:,:,:,N) = &
               (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) &
              +(SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2))

         ! -- Process ISOPN family
         ! -- Get the tracer and species indices
         N  = T_ISOPN    ! ISOPN  tracer  index
         S1 = S_ISOPND   ! ISOPND species index
         S2 = S_ISOPNB   ! ISOPNB species index
         STT(:,:,:,N) = &
               (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) &
              +(SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2))
         
#if defined( UCX )
         ! -- Process CFCX family
         N  = T_CFCX     ! CFCX   tracer  index
         S1 = S_CFC113   ! CFC113 species index
         S2 = S_CFC114   ! CFC114 species index
         S3 = S_CFC115   ! CFC115 species index

         STT(:,:,:,N) = &
               (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) &
              +(SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2)) &
              +(SC%Species(:,:,:,S3)*IO%TRACER_COEFF(N,3))

         ! -- Process HCFCX family
         N  = T_HCFCX      ! HCFCX    tracer  index
         S1 = S_HCFC123    ! HCFC123  species index
         S2 = S_HCFC141b   ! HCFC141b species index
         S3 = S_HCFC142b   ! HCFC142b species index

         STT(:,:,:,N) = &
               (SC%Species(:,:,:,S1)*IO%TRACER_COEFF(N,1)) &
              +(SC%Species(:,:,:,S2)*IO%TRACER_COEFF(N,2)) &
              +(SC%Species(:,:,:,S3)*IO%TRACER_COEFF(N,3))
#endif
      ! E -- N -- D -- O -- F -- K -- L -- U -- D -- G -- E -- -- P -- T -- 2

      ELSE
         write(*,*) 'OPT NOT SET TO 1 OR 2 IN FlexChem_Setup_Mod::FAMILIES_KLUDGE'
      ENDIF
    
    RETURN

  END SUBROUTINE FAMILIES_KLUDGE
!EOC
END MODULE FLEXCHEM_SETUP_MOD

