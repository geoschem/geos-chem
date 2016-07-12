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
  PUBLIC :: FAMILIES_KLUDGE
!    
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  15 Jun 2016 - M. Sulprizio- Remove STTTOCSPEC mapping array. Species and
!                              tracers have a 1:1 mapping currently so mapping
!                              is not required
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
    USE RESTART_MOD,        ONLY : NONADV_SPC_IN_RESTART
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
!                              NAMEGAS with ThisSpc%Name from species database
!  06 Jun 2016 - M. Sulprizio- Replace Get_Indx with Spc_GetIndx to use the
!                              fast-species lookup from the species database
!  06 Jun 2016 - M. Sulprizio- Remove calls to Register_Tracer and
!                              Register_Species; these routines were made
!                              obsolete by the species database
!  14 Jun 2016 - M. Sulprizio- Replace Spc_GetIndx with Ind_ (M. Long)
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

    ! For species unit conversion code
    REAL(fp)           :: CONV_FACTOR        ! [mol/mol] -> [molec/cm3]

    ! Objects
    TYPE(Species), POINTER :: ThisSpc

    !=================================================================
    ! INIT_FLEXCHEM begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Initialize pointer
    ThisSpc => NULL()

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
       ThisSpc => State_Chm%SpcData(N)%Info

       ! Loop over KPP species
       DO N1 = 1, NSPEC

          ! Create vector to map
          ! State_Chm%Species order defined in species_database_mod.F90
          ! to KPP order defined in gckpp_Parameters.F90
          IF ( ADJUSTL( TRIM( SPC_NAMES(N1) ) ) == &
               ADJUSTL( TRIM( ThisSpc%Name  ) ) ) THEN
             CSPECTOKPP(N1)=N
             FOUND=1
             EXIT
          ENDIF

       ENDDO

       ! Print info to log
       IF (FOUND .NE. 1) THEN
          WRITE (6,'(a8,a17)') TRIM( ThisSpc%Name ),   &
             ' NOT found in KPP'
       ENDIF
       IF (FOUND .EQ. 1) THEN
          WRITE (6,'(a8,a29,I4)') TRIM( ThisSpc%Name ), &
             ' was found in KPP with index ', N1
       ENDIF

       ! Free pointer
       ThisSpc => NULL()

    ENDDO

    ! If non-advected species data was not in the NetCDF restart file then 
    ! convert the default background values set in restart_mod.F from 
    ! [mol/mol] to [molec/cm3]. Otherwise, do nothing since concentrations 
    ! for these species were converted to [molec/cm3] after being read in.
    ! Note that advected species concentrations are now always read in 
    ! from the restart file but non-advected species concentrations are
    ! only read in if using a restart file output from a previous run.
    ! (ewl, 7/12/16)
    IF ( .NOT. NONADV_SPC_IN_RESTART ) THEN

       ! Write message to log
       IF ( am_I_Root .and. Input_Opt%LPRT ) THEN
          WRITE(6,110) 
110       FORMAT('     - INIT_FLEXCHEM: converting species background', &
                 ' from [mol mol-1] to [molec/cm3/box]', //, &
                 'Initial value of each species after unit conversion ', &
                 '[molec/cm3]:')
       ENDIF

       ! Loop over GEOS-Chem species
       DO N = 1, State_Chm%nSpecies

          ! Get info about this species from the species database
          ThisSpc => State_Chm%SpcData(N)%Info

          ! Skip non-KPP species
          IF ( .not. ThisSpc%Is_Kpp ) CYCLE

          !$OMP PARALLEL DO &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, L, CONV_FACTOR )
          ! Loop over all chemistry grid boxes
          DO L = 1, LLCHEM 
          DO J = 1, JJPAR
          DO I = 1, IIPAR

             ! Set box-dependent unit conversion factor 
             !CONV_FACTOR = State_Met%AIRNUMDEN(I,J,L)
             CONV_FACTOR = State_Met%PMID_DRY(I,J,L) * 1000e+0_fp / &
                         ( State_Met%T(I,J,L) * ( BOLTZ * 1e+7_fp ) )

             IF ( TRIM( ThisSpc%Name ) /= 'MOH' ) THEN

                ! Convert units from [mol/mol] to [molec/cm3/box]
                State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N) * &
                                             CONV_FACTOR

             ELSE

                !----------------------------------------------------
                ! For methanol (MOH), now use different initial
                ! background concentrations for different regions of
                ! the atmosphere:
                !
                ! (a) 2.0 ppbv MOH -- continental boundary layer
                ! (b) 0.9 ppbv MOH -- marine boundary layer
                ! (c) 0.6 ppbv MOH -- free troposphere
                !
                ! The concentrations listed above are from Heikes et
                ! al, "Atmospheric methanol budget and ocean
                ! implication", _Global Biogeochem. Cycles_, 2002.
                ! These represent the best estimates for the methanol
                ! conc.'s in the troposphere based on various
                ! measurements.
                !
                ! MOH is an inactive chemical species in GEOS-CHEM,
                ! so these initial concentrations will never change.
                ! However, MOH acts as a sink for OH, and therefore
                ! will affect both the OH concentration and the
                ! methylchloroform lifetime.
                !
                ! We specify the MOH concentration as ppbv, but then
                ! we need to multiply by CONV_FACTOR in order to
                ! convert to [molec/cm3].  (bdf, bmy, 2/22/02)
                !----------------------------------------------------

                ! Test for altitude (L < 9 is always in the trop)
                IF ( L <= 9 ) THEN
                   ! Test for ocean/land boxes
                   IF ( State_Met%FRCLND(I,J) >= 0.5 ) THEN
                      ! Continental boundary layer: 2 ppbv MOH
                      State_Chm%Species(I,J,L,N) = 2.000e-9_fp * CONV_FACTOR
                   ELSE
                      ! Marine boundary layer: 0.9 ppbv MOH
                      State_Chm%Species(I,J,L,N) = 0.900e-9_fp * CONV_FACTOR
                   ENDIF
                ELSE
                   ! Test for troposphere
                   IF ( ITS_IN_THE_TROP(I,J,L,State_Met) ) THEN
                      ! Free troposphere: 0.6 ppbv MOH
                      State_Chm%Species(I,J,L,N) = 0.600e-9_fp * CONV_FACTOR
                   ELSE
                      ! Strat/mesosphere:
                      State_Chm%Species(I,J,L,N) = 0.0e+0_fp
                   ENDIF
                ENDIF
             ENDIF

             ! Make a small number if concentration is neg or zero
             State_Chm%SPECIES(I,J,L,N) = &
                        MAX( State_Chm%Species(I,J,L,N), 1d-99 )

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

          ! Print the min and max of each species in [molec/cm3/box]
          ! after conversion if in debug mode (ewl, 3/1/16)
          IF ( Input_Opt%LPRT ) THEN
             WRITE(6,120) N, &
                          TRIM( ThisSpc%Name ), &
                          MINVAL( State_Chm%Species(:,:,1:LLCHEM,N) ), &
                          MAXVAL( State_Chm%Species(:,:,1:LLCHEM,N) ) 
120          FORMAT( 'Species ', i3, ', ', a9, ': Min = ', &
                      es15.9, ', Max = ', es15.9 )
          ENDIF

          ! Free pointer
          ThisSpc => NULL()

       ENDDO ! State_Chm%nSpecies

       ! Mark end of unit conversion in log
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )

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

    !=====================================================================
    ! Save species indices during the Init phase to avoid
    ! repeated costly string matching operations
    !=====================================================================

    ! MMN family
    T_MMN      = IND_( 'MMN',  'T')
    S_MVKN     = IND_( 'MVKN'     )
    S_MACRN    = IND_( 'MACRN'    )

    ! ISOPN family
    T_ISOPN    = IND_( 'ISOPN','T')
    S_ISOPND   = IND_( 'ISOPND'   )
    S_ISOPNB   = IND_( 'ISOPNB'   )

    ! CFCX family
    T_CFCX     = IND_( 'CFCX', 'T')
    S_CFC113   = IND_( 'CFC113'   )
    S_CFC114   = IND_( 'CFC114'   )
    S_CFC115   = IND_( 'CFC115'   )

    ! HCFCX family
    T_HCFCX    = IND_( 'HCFCX','T')
    S_HCFC123  = IND_( 'HCFC123'  )
    S_HCFC141b = IND_( 'HCFC141b' )
    S_HCFC142b = IND_( 'HCFC142b' )

    ! Return to calling program
    RETURN

  END SUBROUTINE INIT_FLEXCHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FAMILIES_KLUDGE
!
! !DESCRIPTION: Subroutine FAMILIES\_KLUDGE is a temporary fix to account for
!  family tracers in FlexChem
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE FAMILIES_KLUDGE(am_I_Root,STT,IO,SC,SM,OPT,RC)
!
! !USES:
!
    USE CMN_SIZE_MOD,         ONLY : IIPAR, JJPAR, LLPAR
    USE GIGC_ErrCode_Mod
    USE PHYSCONSTANTS,        ONLY : AVO
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: IO           ! Short-hand for Input_Opt
    TYPE(MetState), INTENT(IN)    :: SM           ! Short-hand for State_Met
    INTEGER,        INTENT(IN)    :: OPT          ! 1=Trc->Spc, 2=Spc->trc
    LOGICAL,        INTENT(IN)    :: am_I_Root    ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: SC           ! Short-hand for State_Chem
    REAL(kind=fp),  INTENT(INOUT) :: STT(:,:,:,:) ! GC tracer array
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REMARKS:
!  K -- L -- U -- D -- G -- E -- C -- O -- D -- E -- ! -- ! -- ! -- !
!  THIS IS A TEMPORARY FIX AND NEEDS TO BE RESOLVED THROUGHOUT
!  THE MODEL AS SOON AS IT APPEARS POSSIBLE TO DO SO!
!  -- The following code ensures that the two remaining tracer
!     families are dealt with appropriately.
!     The Tracer MMN is a sum of the species MVKN and MACRN
!     The Tracer ISOPN is a sum of ISOPND and ISOPNB
!        The tracer restart file for the pre-flex GEOS-Chem, includes
!     these families, but the removal of the routines
!     "lump" and "partition" killed the code resposible for them.
!        Here, we want to install a hard-code fix with the complete
!     expectation that we will no longer use species families. Thus,
!     when possible, the two remaining families need to be removed, 
!     and this kludge disabled.
!    M.S.L. - Jan., 5 2016
!   
! !REVISION HISTORY:
!  05 Jan 2016 - M. Long     - Initial version
!  28 Mar 2016 - R. Yantosca - Prevent div-by-zero statements. Renamed SUM to
!                              QSUM to avoid conflicts with the Fortran
!                              intrinsic function SUM.
!  14 Jun 2016 - M. Sulprizio- Add ProTeX headers
!  14 Jun 2016 - M. Sulprizio- Handle unit conversions between kg and molec/cm3
!                              for family tracers in this routine instead of in
!                              flex_chemdr.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER        :: N,S1,S2,S3,I,J,L
    REAL(kind=fp)  :: QSUM(IIPAR,JJPAR,LLPAR)
    REAL(kind=fp)  :: QTEMP(IIPAR,JJPAR,LLPAR)
    REAL(kind=fp)  :: STTTEMP(IIPAR,JJPAR,LLPAR)
    REAL(fp)       :: MW_kg


    ! Assume success
    RC = GIGC_SUCCESS

    !-------------------------------------------------------------------------
    ! Part 1: From Tracers to Species
    !-------------------------------------------------------------------------
    IF (OPT .eq. 1) THEN

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%% Process MMN family
       !%%%%
       !%%%% NOTE: MMN = 1.0 MVKN + 1.0 MACRN
       !%%%% so TRACER_COEFF is always 1 in this case.
       !%%%% This will let us get rid of IO%TRACER_COEFF in the code.
       !%%%% (bmy, 5/17/16)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       N  = T_MMN     ! MMN   tracer  index
       S1 = S_MVKN    ! MVKN  species index
       S2 = S_MACRN   ! MACRN species index

       ! Moles C per moles species
       MW_kg      = SC%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Convert concentrations from kg to molec/cm3 for SC%Species
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          STTTEMP(I,J,L) = STT(I,J,L,N) * &
                           ( AVO / MW_KG ) / ( SM%AIRVOL(I,J,L) * 1e+6_fp )
       ENDDO
       ENDDO
       ENDDO

       ! Sum of constituents
       QSUM   = SC%Species(:,:,:,S1) + SC%Species(:,:,:,S2)

       ! -- -- First, do MVKN
       WHERE( QSUM > 0.0_fp ) 
          QTEMP = SC%Species(:,:,:,S1) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S1) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )

       ! -- -- Then, do MACRN
       WHERE( QSUM > 0.0_fp ) 
          QTEMP = SC%Species(:,:,:,S2) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S2) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%% Process ISOPN family
       !%%%%
       !%%%% NOTE: ISOPN = 1.0 ISOPND + 1.0 ISOPNB
       !%%%% so TRACER_COEFF is always 1 in this case.
       !%%%% This will let us get rid of IO%TRACER_COEFF in the code.
       !%%%% (bmy, 5/17/16)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       N  = T_ISOPN    ! ISOPN  tracer  index
       S1 = S_ISOPND   ! ISOPND species index
       S2 = S_ISOPNB   ! ISOPNB species index

       ! Moles C per moles species
       MW_kg      = SC%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Convert concentrations from kg to molec/cm3 for SC%Species
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          STTTEMP(I,J,L) = STT(I,J,L,N) * &
                           ( AVO / MW_KG ) / ( SM%AIRVOL(I,J,L) * 1e+6_fp )
       ENDDO
       ENDDO
       ENDDO

       ! Sum of constitutents
       QSUM = SC%Species(:,:,:,S1) + SC%Species(:,:,:,S2)

       ! -- -- First, do ISOPND
       WHERE( QSUM > 0.0_fp ) 
          QTEMP = ( SC%Species(:,:,:,S1) ) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S1) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )

       ! -- -- Then, do ISOPNB
       WHERE( QSUM > 0.0_fp ) 
          QTEMP = ( SC%Species(:,:,:,S2) ) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S2) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )

#if defined( UCX )
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%% Process CFCX family
       !%%%%
       !%%%% NOTE: CFCX = 1.0 CFC113 + 1.0 CFC114 + 1.0 CFC115
       !%%%% so TRACER_COEFF is always 1 in this case.
       !%%%% This will let us get rid of IO%TRACER_COEFF in the code.
       !%%%% (bmy, 5/17/16)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       N  = T_CFCX     ! CFCX   tracer  index
       S1 = S_CFC113   ! CFC113 species index
       S2 = S_CFC114   ! CFC114 species index
       S3 = S_CFC115   ! CFC115 species index

       ! Moles C per moles species
       MW_kg      = SC%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Convert concentrations from kg to molec/cm3 for SC%Species
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          STTTEMP(I,J,L) = STT(I,J,L,N) * &
                           ( AVO / MW_KG ) / ( SM%AIRVOL(I,J,L) * 1e+6_fp )
       ENDDO
       ENDDO
       ENDDO

       ! Sum of constituents
       QSUM  = SC%Species(:,:,:,S1) &
             + SC%Species(:,:,:,S2) &
             + SC%Species(:,:,:,S3) 

       ! -- -- First, do CFC113
       WHERE( QSUM > 0.0_fp )
          QTEMP = SC%Species(:,:,:,S1) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S1) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )

       ! -- -- Then, do CFC114
       WHERE( QSUM > 0.0_fp )
          QTEMP = SC%Species(:,:,:,S2) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S2) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )

       ! -- -- Then, do CFC115
       WHERE( QSUM > 0.0_fp )
          QTEMP = SC%Species(:,:,:,S3) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S3) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%% Process HCFCX family
       !%%%%
       !%%%% NOTE: HCFCX = 1.0 HCFC123 + 1.0 HCFC141b + 1.0 HFCC142b
       !%%%% so TRACER_COEFF is always 1 in this case.
       !%%%% This will let us get rid of IO%TRACER_COEFF in the code.
       !%%%% (bmy, 5/17/16)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       N  = T_HCFCX      ! HCFCX    tracer  index
       S1 = S_HCFC123    ! HCFC123  species index
       S2 = S_HCFC141b   ! HCFC141b species index
       S3 = S_HCFC142b   ! HCFC142b species index

       ! Moles C per moles species
       MW_kg      = SC%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Convert concentrations from kg to molec/cm3 for SC%Species
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          STTTEMP(I,J,L) = STT(I,J,L,N) * &
                           ( AVO / MW_KG ) / ( SM%AIRVOL(I,J,L) * 1e+6_fp )
       ENDDO
       ENDDO
       ENDDO
      
       ! Sum of constituents
       QSUM  = SC%Species(:,:,:,S1) &
             + SC%Species(:,:,:,S2) &
             + SC%Species(:,:,:,S3)

       ! -- -- First, do HCFC123
       WHERE( QSUM > 0.0_fp )
          QTEMP = SC%Species(:,:,:,S1) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S1) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )

       ! -- -- Then, do HCFC141b
       WHERE( QSUM > 0.0_fp )
          QTEMP = SC%Species(:,:,:,S2) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S2) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )

       ! -- -- Then, do HCFC142b
       WHERE( QSUM > 0.0_fp )
          QTEMP = SC%Species(:,:,:,S3) / QSUM
       ELSEWHERE
          QTEMP = 0.0_fp
       ENDWHERE
       SC%Species(:,:,:,S3) = MAX( QTEMP*STTTEMP(:,:,:), 1.E-99_fp )
#endif

    !------------------------------------------------------------------------
    ! Part 2: From Species to Tracers
    !------------------------------------------------------------------------
    ELSEIF (OPT .eq. 2) THEN

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%% Process MMN family
       !%%%%
       !%%%% NOTE: MMN = 1.0 MVKN + 1.0 MACRN
       !%%%% so TRACER_COEFF is always 1 in this case.
       !%%%% This will let us get rid of IO%TRACER_COEFF in the code.
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       N  = T_MMN     ! MMN   species index
       S1 = S_MVKN    ! MVKN  species index
       S2 = S_MACRN   ! MACRN species index

       STT(:,:,:,N) = SC%Species(:,:,:,S1) + SC%Species(:,:,:,S2)

       ! Moles C per moles species
       MW_kg      = SC%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Convert concentrations from molec/cm3 to kg for STT
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          STT(I,J,L,N) = STT(I,J,L,N) / &
                         ( AVO / MW_KG ) * ( SM%AIRVOL(I,J,L) * 1e+6_fp )
       ENDDO
       ENDDO
       ENDDO

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%% Process ISOPN family
       !%%%%
       !%%%% NOTE: ISOPN = 1.0 ISOPND + 1.0 ISOPNB
       !%%%% so TRACER_COEFF is always 1 in this case.
       !%%%% This will let us get rid of IO%TRACER_COEFF in the code.
       !%%%% (bmy, 5/17/16)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       N  = T_ISOPN    ! ISOPN  tracer  index
       S1 = S_ISOPND   ! ISOPND species index
       S2 = S_ISOPNB   ! ISOPNB species index

       STT(:,:,:,N) = SC%Species(:,:,:,S1) + SC%Species(:,:,:,S2)
         
       ! Moles C per moles species
       MW_kg      = SC%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Convert concentrations from molec/cm3 to kg for STT
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          STT(I,J,L,N) = STT(I,J,L,N) / &
                         ( AVO / MW_KG ) * ( SM%AIRVOL(I,J,L) * 1e+6_fp )
       ENDDO
       ENDDO
       ENDDO

#if defined( UCX )
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%% Process CFCX family
       !%%%%
       !%%%% NOTE: CFCX = 1.0 CFC113 + 1.0 CFC114 + 1.0 CFC115
       !%%%% so TRACER_COEFF is always 1 in this case.
       !%%%% This will let us get rid of IO%TRACER_COEFF in the code.
       !%%%% (bmy, 5/17/16)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       N  = T_CFCX     ! CFCX   tracer  index
       S1 = S_CFC113   ! CFC113 species index
       S2 = S_CFC114   ! CFC114 species index
       S3 = S_CFC115   ! CFC115 species index

       STT(:,:,:,N) = SC%Species(:,:,:,S1) &
                    + SC%Species(:,:,:,S2) &
                    + SC%Species(:,:,:,S3)

       ! Moles C per moles species
       MW_kg      = SC%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Convert concentrations from molec/cm3 to kg for STT
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          STT(I,J,L,N) = STT(I,J,L,N) / &
                         ( AVO / MW_KG ) * ( SM%AIRVOL(I,J,L) * 1e+6_fp )
       ENDDO
       ENDDO
       ENDDO

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%% Process HCFCX family
       !%%%%
       !%%%% NOTE: HCFCX = 1.0 HCFC123 + 1.0 HCFC141b + 1.0 HFCC142b
       !%%%% so TRACER_COEFF is always 1 in this case.
       !%%%% This will let us get rid of IO%TRACER_COEFF in the code.
       !%%%% (bmy, 5/17/16)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       N  = T_HCFCX      ! HCFCX    tracer  index
       S1 = S_HCFC123    ! HCFC123  species index
       S2 = S_HCFC141b   ! HCFC141b species index
       S3 = S_HCFC142b   ! HCFC142b species index

       STT(:,:,:,N) = SC%Species(:,:,:,S1) &
                    + SC%Species(:,:,:,S2) &
                    + SC%Species(:,:,:,S3)

       ! Moles C per moles species
       MW_kg      = SC%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Convert concentrations from molec/cm3 to kg for STT
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          STT(I,J,L,N) = STT(I,J,L,N) / &
                         ( AVO / MW_KG ) * ( SM%AIRVOL(I,J,L) * 1e+6_fp )
       ENDDO
       ENDDO
       ENDDO
#endif

    ELSE
       write(*,*) 'OPT NOT SET TO 1 OR 2 IN FlexChem_Setup_Mod::FAMILIES_KLUDGE'
    ENDIF
    
    RETURN

  END SUBROUTINE FAMILIES_KLUDGE
!EOC
END MODULE FLEXCHEM_SETUP_MOD

