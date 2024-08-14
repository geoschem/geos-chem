!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: photolysis_mod.F90
!
! !DESCRIPTION: Module PHOTOLYSIS\_MOD contains routines and variables
!  for GEOS-Chem photolysis.
!\\
!\\
! !INTERFACE:
!
MODULE PHOTOLYSIS_MOD
!
! !USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_Photolysis
  PUBLIC  :: Do_Photolysis
  PUBLIC  :: PHOTRATE_ADJ
!
! !PRIVATE MEMBER FUNCTIONS:
!
!ewl  PRIVATE :: RD_AOD
!ewl  PRIVATE :: CALC_AOD
  PRIVATE :: SET_AER
  PRIVATE :: RD_PROF_NC
!
! !REVISION HISTORY:
!  20 Mar 2023 - E. Lundgren - initial version, adapted from fast_jx_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
   ! Species ID flags
   INTEGER :: id_NIT, id_NITs, id_SALA, id_SALC

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_photolysis
!
! !DESCRIPTION: Subroutine INIT\_PHOTOLYSIS initializes GEOS-Chem photolysis.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_PHOTOLYSIS( Input_Opt, State_Grid, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE Charpak_Mod,    ONLY : CSTRIP
#ifdef FASTJX
    USE CMN_FJX_Mod,    ONLY : JVN_, W_, JLABEL, RNAMES, WL, JFACTA
#else
    USE Cldj_Cmn_Mod,   ONLY : JVN_, W_, JLABEL, RNAMES, WL, JFACTA
#endif
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : Planck, CConst
    USE State_Chm_Mod,  ONLY : ChmState, Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
#ifdef FASTJX
    USE Fjx_Interface_Mod,  ONLY : Init_FastJX
#else
    USE Cldj_Interface_Mod, ONLY : Init_CloudJ
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(DgnState), INTENT(IN)    :: State_Diag  ! Diagnostics State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 Mar 2023 - E. Lundgren - initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=50 ) :: TEXT
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Scalars
    LOGICAL            :: notDryRun
    INTEGER            :: J, K
    REAL(fp)           :: ND64MULT

    ! Pointers
    INTEGER,   POINTER :: GC_Photo_ID(:)

    !=======================================================================
    ! INIT_PHOTOLYSIS begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
      ' -> at Init_Photolysis  (in module GeosCore/photolysis_mod.F90)'

    ! Set pointers
    GC_Photo_ID => State_Chm%Phot%GC_Photo_ID

    !--------------------------------------------------------------------
    ! Initialize Fast-JX if photolysis enabled
    !
    ! NOTE: we need to call this for a dry-run so that we can get
    ! a list of all of the lookup tables etc read by Fast-JX
    !--------------------------------------------------------------------
    IF ( Input_Opt%Do_Photolysis ) THEN
#ifdef FASTJX
       IF ( TRIM(Input_Opt%Fast_JX_Dir) == MISSING_STR ) THEN
          ErrMsg = 'Init_Photolysis: Fast-JX directory missing in ' // &
                   'in geoschem_config.yml!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CALL Init_FastJX( Input_Opt, State_Diag, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_FastJX"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
#else
       IF ( TRIM(Input_Opt%CloudJ_Dir) == MISSING_STR ) THEN
          ErrMsg = 'Init_Photolysis: Cloud-J directory missing in ' // &
                   'geoschem_config.yml!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CALL Init_CloudJ( Input_Opt, State_Grid, State_Diag, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_CloudJ"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
#endif
    ENDIF
    
!ewl    !------------------------------------------------------------------------
!ewl    ! Read in AOD data even if photolysis disabled
!ewl    ! (or just print file name if in dry-run mode)
!ewl    !------------------------------------------------------------------------
!ewl    CALL RD_AOD( Input_Opt, State_Chm, RC )
!ewl    IF ( RC /= GC_SUCCESS ) THEN
!ewl       ErrMsg = 'Error encountered in routine "RD_AOD"!'
!ewl       CALL GC_Error( ErrMsg, RC, ThisLoc )
!ewl       RETURN
!ewl    ENDIF

    !--------------------------------------------------------------------
    ! Read in T & O3 climatology to fill e.g. upper layers or if O3 not calc.
    !--------------------------------------------------------------------
    ! NOTE: Cloud-J reads in an ascii file with this data during initialization
    ! and uses it prior to calling Cloud_JX within the Cloud-J standalone. In
    ! GEOS-Chem we read a netcdf file instead and use the data within
    ! subroutine Set_Prof_Fjx if using Fast-JX and Set_Prof_CloudJ if using
    ! Cloud-J. The data is stored in State_Chm%Phot%TREF/%OREF. Cloud-J
    ! globals variables TREF and OREF are only used for Cloud-J standalone.
    CALL RD_PROF_NC( Input_Opt, State_Grid, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Rd_Prof_Nc"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Exit without doing any computations if we are doing a dry-run
    !------------------------------------------------------------------------
    IF ( Input_Opt%DryRun ) RETURN

!ewl    !------------------------------------------------------------------------
!ewl    ! Compute the required wavelengths in the LUT to calculate requested AOD
!ewl    !------------------------------------------------------------------------
!ewl    IF (Input_Opt%amIRoot) WRITE(6,*) 'Wavelength optics read successfully'
!ewl    CALL CALC_AOD( Input_Opt, State_Chm, RC )

    !------------------------------------------------------------------------
    ! Exit if photolysis disabled (zero J-values)
    !------------------------------------------------------------------------
    IF ( .NOT. Input_Opt%Do_Photolysis ) RETURN

    !--------------------------------------------------------------------
    ! Set up MIEDX array to interpret between GC and FJX aerosol indexing
    !------------------------------------------------------------------------
    CALL SET_AER( Input_Opt, State_Chm, RC )

    !========================================================================
    ! Flag special reactions that will be later adjusted by
    ! routine PHOTRATE_ADJ (called from FlexChem)
    !========================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       ! Loop over all photolysis reactions
       DO K = 1, State_Chm%Phot%nPhotRxns

          ! Strip all blanks from the reactants and products list
          TEXT = JLABEL(K)
          CALL CSTRIP( TEXT )

          !IF ( Input_Opt%amIRoot ) THEN
          !   WRITE(*,*) K, TRIM( TEXT )
          !ENDIF

          ! Look for certain reactions
          SELECT CASE( TRIM( TEXT ) )
             CASE( 'O2PHOTONOO' )
                State_Chm%Phot%RXN_O2 = K     ! O2 + hv -> O + O
             CASE( 'O3PHOTONO2O' )
                State_Chm%Phot%RXN_O3_1 = K   ! O3 + hv -> O2 + O
             CASE( 'O3PHOTONO2O(1D)' )
                State_Chm%Phot%RXN_O3_2 = K   ! O3 + hv -> O2 + O(1D)
             CASE( 'SO4PHOTONSO2OHOH' )
                State_Chm%Phot%RXN_H2SO4 = K  ! SO4 + hv -> SO2 + OH + OH
             CASE( 'NO2PHOTONNOO' )
                State_Chm%Phot%RXN_NO2 = K    ! NO2 + hv -> NO + O
             CASE( 'NOPHOTONNO' )
                State_Chm%Phot%RXN_NO = K     ! NO + hv -> N + O
             CASE( 'NO3PHOTONNO2O' )
                State_Chm%Phot%RXN_NO3 = K    ! NO3 + hv -> NO2 + O
             CASE( 'N2OPHOTONN2O' )
                State_Chm%Phot%RXN_N2O = K    ! N2O + hv -> N2 + O
             CASE( 'NITsPHOTONHNO2' )
                State_Chm%Phot%RXN_JNITSa = K ! NITs + hv -> HNO2
             CASE( 'NITsPHOTONNO2' )
                State_Chm%Phot%RXN_JNITSb = K ! NITs + hv -> NO2
             CASE( 'NITPHOTONHNO2' )
                State_Chm%Phot%RXN_JNITa = K  ! NIT + hv -> HNO2
             CASE( 'NITPHOTONNO2' )
                State_Chm%Phot%RXN_JNITb = K  ! NIT + hv -> NO2
             CASE( 'HNO3PHOTONNO2OH' )
                State_Chm%Phot%RXN_JHNO3 = K  ! HNO3 + hv = OH + NO2
             CASE DEFAULT
                ! Nothing
          END SELECT
       ENDDO

       !---------------------------------------------------------------------
       ! Error check the various rxn flags
       !---------------------------------------------------------------------
       IF ( State_Chm%Phot%RXN_O2 < 0 ) THEN
          ErrMsg = 'Could not find rxn O2 + hv -> O + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_O3_1 < 0 ) THEN
          ErrMsg = 'Could not find rxn O3 + hv -> O2 + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_O3_2 < 0 ) THEN
          ErrMsg = 'Could not find rxn O3 + hv -> O2 + O(1D)'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
       ENDIF

       IF ( State_Chm%Phot%RXN_NO2 < 0 ) THEN
          ErrMsg = 'Could not find rxn NO2 + hv -> NO + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_NO2 < 0 ) THEN
          ErrMsg = 'Could not find rxn NO2 + hv -> NO + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_JNITSa < 0 ) THEN
          ErrMsg = 'Could not find rxn NITS + hv -> HNO2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_JNITSb < 0 ) THEN
          ErrMsg = 'Could not find rxn NITS + hv -> NO2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_JNITa < 0 ) THEN
          ErrMsg = 'Could not find rxn NIT + hv -> HNO2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_JNITb < 0 ) THEN
          ErrMsg = 'Could not find rxn NIT + hv -> NO2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_H2SO4  < 0 ) THEN
          ErrMsg = 'Could not find rxn SO4 + hv -> SO2 + OH + OH!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_NO3 < 0 ) THEN
          ErrMsg = 'Could not find rxn NO3 + hv -> NO2 + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_NO < 0 ) THEN
          ErrMsg = 'Could not find rxn NO + hv -> O + N'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Phot%RXN_N2O < 0 ) THEN
          ErrMsg = 'Could not find rxn N2O + hv -> N2 + O(1D)'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Print out saved rxn flags for fullchem simulations
       !---------------------------------------------------------------------
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 100 ) REPEAT( '=', 79 )
          WRITE( 6, 110 )
          WRITE( 6, 120 ) State_Chm%Phot%RXN_O2
          WRITE( 6, 130 ) State_Chm%Phot%RXN_O3_1
          WRITE( 6, 140 ) State_Chm%Phot%RXN_O3_2
          WRITE( 6, 180 ) State_Chm%Phot%RXN_JNITSa
          WRITE( 6, 190 ) State_Chm%Phot%RXN_JNITSb
          WRITE( 6, 200 ) State_Chm%Phot%RXN_JNITa
          WRITE( 6, 210 ) State_Chm%Phot%RXN_JNITb
          WRITE( 6, 160 ) State_Chm%Phot%RXN_H2SO4
          WRITE( 6, 170 ) State_Chm%Phot%RXN_NO2
          WRITE( 6, 100 ) REPEAT( '=', 79 )
       ENDIF
    ENDIF

    !========================================================================
    ! Flag reactions for diagnostics (only in Hg chem)
    !========================================================================
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
        ! Loop over all photolysis reactions
        DO K = 1, State_Chm%Phot%nPhotRxns

           ! Strip all blanks from the reactants and products list
           TEXT = JLABEL(K)
           CALL CSTRIP( TEXT )

           ! Look for certain reactions
           SELECT CASE( TRIM( TEXT ) )
              CASE( 'O3PHOTONO2O' )
                 State_Chm%Phot%RXN_O3_1 = K ! O3 + hv -> O2 + O
              CASE( 'O3PHOTONO2O(1D)' )
                 State_Chm%Phot%RXN_O3_2 = K ! O3 + hv -> O2 + O(1D)
              CASE( 'NO2PHOTONNOO' )
                 State_Chm%Phot%RXN_NO2 = K  ! NO2 + hv -> NO + O
              CASE( 'BrOPHOTONBrO' )
                 State_Chm%Phot%RXN_BrO = K  ! BrO + hv -> Br + O
              CASE( 'ClOPHOTONClO' )
                 State_Chm%Phot%RXN_ClO = K  ! ClO + hv -> Cl + O
              CASE DEFAULT
                 ! Nothing
           END SELECT
        ENDDO

        !--------------------------------------------------------------------
        ! Error check the various rxn flags
        !--------------------------------------------------------------------
        IF ( State_Chm%Phot%RXN_O3_1 < 0 ) THEN
           ErrMsg = 'Could not find rxn O3 + hv -> O2 + O'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

        IF ( State_Chm%Phot%RXN_O3_2 < 0 ) THEN
           ErrMsg = 'Could not find rxn O3 + hv -> O2 + O(1D) #1'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

        IF ( State_Chm%Phot%RXN_NO2 < 0 ) THEN
           ErrMsg = 'Could not find rxn NO2 + hv -> NO + O'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

        IF ( State_Chm%Phot%RXN_BrO < 0 ) THEN
           ErrMsg = 'Could not find rxn BrO + hv -> Br + O'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

        IF ( State_Chm%Phot%RXN_ClO < 0 ) THEN
           ErrMsg = 'Could not find rxn ClO + hv -> Cl + O'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

       !---------------------------------------------------------------------
       ! Print out saved rxn flags for Hg simulation
       !---------------------------------------------------------------------
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 100 ) REPEAT( '=', 79 )
          WRITE( 6, 110 )
          WRITE( 6, 130 ) State_Chm%Phot%RXN_O3_1
          WRITE( 6, 140 ) State_Chm%Phot%RXN_O3_2
          WRITE( 6, 170 ) State_Chm%Phot%RXN_NO2
          WRITE( 6, 220 ) State_Chm%Phot%RXN_BrO
          WRITE( 6, 230 ) State_Chm%Phot%RXN_ClO
          WRITE( 6, 100 ) REPEAT( '=', 79 )
       ENDIF
    ENDIF

    ! Skip further processing if we are in dry-run mode
    IF (  .not. Input_Opt%DryRun ) THEN

       ! Define species IDs
       id_NIT  = IND_('NIT')
       id_NITs = IND_('NITs')
       id_SALA = IND_('SALA')
       id_SALC = IND_('SALC')

       ! Get the GEOS-Chem photolysis index for each of the 1...JVN_ entries
       ! in the FJX_j2j.dat file.  We'll use this for the diagnostics.
       DO J = 1, State_Chm%Phot%nMaxPhotRxns

          IF ( J == State_Chm%Phot%Rxn_O3_2 ) THEN

             !------------------------------------------------------------
             ! O3 + hv = O + O(1D)
             !
             ! Save this as JO3_O1D in the nPhotol+1 slot
             !------------------------------------------------------------
             GC_Photo_Id(J) = State_Chm%nPhotol + 1

          ELSE IF ( J == State_Chm%Phot%Rxn_O3_1 ) THEN

             !------------------------------------------------------------
             ! O3 + hv -> O + O
             !
             ! Save this as JO3_O3P in the nPhotol+2 slot
             !-------------------------------------------------------------
             GC_Photo_Id(J) = State_Chm%nPhotol + 2

          ELSE

             !------------------------------------------------------------
             ! Everything else
             !
             ! Find the matching GEOS-Chem photolysis species number
             !------------------------------------------------------------
             GC_Photo_Id(J) = Ind_( RNAMES(J), 'P' )

          ENDIF

          ! Print the mapping
          IF ( Input_Opt%amIRoot ) THEN
             IF ( GC_Photo_Id(J) > 0 ) THEN
                WRITE(6, 240) RNAMES(J), J, GC_Photo_Id(J), JFACTA(J)
240             FORMAT( a10, ':', i7, 2x, i7, 2x, f7.4 )
             ENDIF
          ENDIF
       ENDDO

       !=====================================================================
       ! Compute factors for UV flux diagnostics if turned on
       !=====================================================================
       IF ( State_Diag%Archive_UVFluxNet      .or. &
            State_Diag%Archive_UVFluxDirect   .or. &
            State_Diag%Archive_UVFluxDiffuse ) THEN
          ND64MULT  = PLANCK * CCONST * 1.0e+13_fp
          State_Chm%Phot%UVXFACTOR = 0e+0_fp
          DO J = 1, W_
             State_Chm%Phot%UVXFACTOR(J) = ND64MULT/WL(J)
          ENDDO
       ENDIF
    ENDIF

    ! Free pointers
    GC_Photo_ID => NULL()
    
    ! FORMAT statements
100 FORMAT( a                                                 )
110 FORMAT( 'Photo rxn flags saved for use in PHOTRATE_ADJ:', / )
120 FORMAT( 'RXN_O2     [ O2   + hv -> O + O         ]  =  ', i5 )
130 FORMAT( 'RXN_O3_1   [ O3   + hv -> O2 + O        ]  =  ', i5 )
140 FORMAT( 'RXN_O3_2a  [ O3   + hv -> O2 + O(1D) #1 ]  =  ', i5 )
150 FORMAT( 'RXN_O3_2b  [ O3   + hv -> O2 + O(1D) #2 ]  =  ', i5 )
160 FORMAT( 'RXN_H2SO4  [ SO4  + hv -> SO2 + OH + OH ]  =  ', i5 )
170 FORMAT( 'RXN_NO2    [ NO2  + hv -> NO + O        ]  =  ', i5 )
180 FORMAT( 'RXN_JNITSa [ NITS + hv -> HNO2          ]  =  ', i5 )
190 FORMAT( 'RXN_JNITSb [ NITS + hv -> NO2           ]  =  ', i5 )
200 FORMAT( 'RXN_JNITa  [ NIT  + hv -> HNO2          ]  =  ', i5 )
210 FORMAT( 'RXN_JNITb  [ NIT  + hv -> NO2           ]  =  ', i5 )
220 FORMAT( 'RXN_BrO    [ BrO  + hv -> Br + O        ]  =  ', i5 )
230 FORMAT( 'RXN_ClO    [ ClO  + hv -> Cl + O        ]  =  ', i5 )

  END SUBROUTINE Init_Photolysis
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_photolysis
!
! !DESCRIPTION: Subroutine DO\_PHOTOLYSIS loops over longitude and latitude,
!  and computes J-Values for each column at every chemistry time-step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Photolysis( Input_Opt, State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
#ifdef FASTJX
    USE Fjx_Interface_Mod,  ONLY : Run_FastJX
#else
    USE Cldj_Interface_Mod, ONLY : Run_CloudJ
#endif

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
!
! !REVISION HISTORY:
!  20 Mar 2023 - E. Lundgren - initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: Wavelength
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! DO_PHOTOLYSIS begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at DO_PHOTOLYSIS (in module GeosCore/photolysis_mod.F90)'
    WAVELENGTH = 0

#ifdef FASTJX
    CALL Run_FastJX( Wavelength, Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC )
#else
    CALL Run_CloudJ( Input_Opt, State_Chm, State_Diag, State_Grid, State_Met, RC )
#endif

    IF ( RC /= GC_SUCCESS ) THEN
#ifdef FASTJX
       ErrMsg = 'Error encountered in subroutine Run_FastJX!'
#else
       ErrMsg = 'Error encountered in subroutine Run_CloudJ!'
#endif
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Do_Photolysis
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photrate_adj
!
! !DESCRIPTION: Subroutine PHOTRATE\_ADJ adjusts certain photolysis rates
!  for chemistry.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PHOTRATE_ADJ( Input_Opt, State_Chm,  State_Diag, State_Met,     &
                           I, J, L, FRAC, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
    USE Error_Mod,      ONLY : SAFE_DIV
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input_Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm  ! Chemistry State object
    TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
    INTEGER,        INTENT(IN)    :: I, J, L    ! Lon, lat, lev indices
    REAL(fp),       INTENT(IN)    :: FRAC       ! Result of SO4_PHOTFRAC,
                                                !  called from DO_FLEXCHEM
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
!
! !REMARKS:
!  NOTE: The netCDF diagnostics are attached in DO_FLEXCHEM so that we have
!  access to the adjusted rates.
!    -- Bob Yantosca, 19 Dec 2017
!
!  %%%% NOTE: WE SHOULD UPDATE THE COMMENTS TO MAKE SURE THAT WE DO      %%%%
!  %%%% NOT KEEP ANY CONFLICTING OR INCORRECT INFORMATION (bmy, 3/28/16) %%%%
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: RXN_JNITSa, RXN_JNITSb, RXN_JNITa, RXN_JNITb
    INTEGER  :: RXN_JHNO3,  RXN_H2SO4,  RXN_O3_1,  RXN_O3_2
    REAL(fp) :: JscaleNITs, JscaleNIT,  JNITChanA, JNITChanB
    REAL(fp) :: C_O2,     C_N2, C_H2,   ITEMPK, RO1DplH2O
    REAL(fp) :: RO1DplH2, RO1D, NUMDEN, TEMP,   C_H2O
    REAL(fp) :: C_NIT, C_NITs, C_SALA, C_SALC, FAC

    ! Pointers
    REAL*8, POINTER :: ZPJ(:,:,:,:)

    !=================================================================
    ! PHOTRATE_ADJ begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    TEMP    = State_Met%T(I,J,L)                                 ! K
    NUMDEN  = State_Met%AIRNUMDEN(I,J,L)                         ! molec/cm3
    C_H2O   = State_Met%AVGW(I,J,L) * State_Met%AIRNUMDEN(I,J,L) ! molec/cm3
 
    ! Reaction flags
    RXN_JNITSa = State_Chm%Phot%RXN_JNITSa
    RXN_JNITSb = State_Chm%Phot%RXN_JNITSb
    RXN_JNITa  = State_Chm%Phot%RXN_JNITa
    RXN_JNITb  = State_Chm%Phot%RXN_JNITb
    RXN_JHNO3  = State_Chm%Phot%RXN_JHNO3
    RXN_H2SO4  = State_Chm%Phot%RXN_H2SO4
    RXN_O3_1   = State_Chm%Phot%RXN_O3_1
    RXN_O3_2   = State_Chm%Phot%RXN_O3_2

    ! Pointers
    ZPJ => State_Chm%Phot%ZPJ

    ! For all mechanisms. Set the photolysis rate of NITs and NIT to a
    ! scaled value of JHNO3. NOTE: this is set in geoschem_config.yml
    ! Allow particulate nitrate photolysis in the troposphere only
    IF ( Input_Opt%hvAerNIT .and.               &
         State_Met%InTroposphere(I,J,L) ) THEN

       ! Get NIT and NITs concentrations [molec cm-3]
       C_NIT      = State_Chm%Species(id_NIT)%Conc(I,J,L)
       C_NITs     = State_Chm%Species(id_NITs)%Conc(I,J,L)


       ! Get sea-salt concentrations [molec cm-3]
       C_SALA     = State_Chm%Species(id_SALA)%Conc(I,J,L)
       C_SALC     = State_Chm%Species(id_SALC)%Conc(I,J,L)

       ! Scaling factor for J(NIT)
       FAC        = SAFE_DIV( C_SALA, C_SALA + C_NIT, 1e+0_fp )

       ! Set FRAC_NIT to a minimum of 0.1
       FAC        = MAX( 0.1e+0_fp, FAC )

       JscaleNITs = Input_Opt%hvAerNIT_JNITs
       JscaleNIT  = Input_Opt%hvAerNIT_JNIT

       ! convert reaction channel % to a fraction
       JNITChanA  = Input_Opt%JNITChanA
       JNITChanB  = Input_Opt%JNITChanB
       JNITChanA  = JNITChanA / 100.0_fp
       JNITChanB  = JNITChanB / 100.0_fp

       ! Set the photolysis rate of NITs
       ZPJ(L,RXN_JNITSa,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNITs
       ZPJ(L,RXN_JNITSb,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNITs

       ! Set the photolysis rate of NIT
       ZPJ(L,RXN_JNITa,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNIT * FAC
       ZPJ(L,RXN_JNITb,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNIT * FAC

       ! NOTE: channel scaling is 1 in FJX_j2j.dat, then updated here
       ZPJ(L,RXN_JNITSa,I,J) = ZPJ(L,RXN_JNITSa,I,J) * JNITChanA
       ZPJ(L,RXN_JNITa,I,J) = ZPJ(L,RXN_JNITa,I,J) * JNITChanA
       ZPJ(L,RXN_JNITSb,I,J) = ZPJ(L,RXN_JNITSb,I,J) * JNITChanB
       ZPJ(L,RXN_JNITb,I,J) = ZPJ(L,RXN_JNITb,I,J) * JNITChanB

    ! Gotcha to set JNIT and JNITs to zero if hvAerNIT switch is off
    ELSE

       ! Set the photolysis rate of NITs to zero
       ZPJ(L,RXN_JNITSa,I,J) = 0.0_fp
       ZPJ(L,RXN_JNITSb,I,J) = 0.0_fp
       ! Set the photolysis rate of NIT to zero
       ZPJ(L,RXN_JNITa,I,J) = 0.0_fp
       ZPJ(L,RXN_JNITb,I,J) = 0.0_fp

    ENDIF

    !==============================================================
    ! SPECIAL TREATMENT FOR H2SO4+hv -> SO2 + 2OH
    !
    ! Only allow photolysis of H2SO4 when gaseous (SDE 04/11/13)
    !==============================================================

    ! Calculate if H2SO4 expected to be gaseous or aqueous
    ! Only allow photolysis above 6 hPa
    ! RXN_H2SO4 specifies SO4 + hv -> SO2 + OH + OH
    ZPJ(L,RXN_H2SO4,I,J) = ZPJ(L,RXN_H2SO4,I,J) * FRAC

    !==============================================================
    ! SPECIAL TREATMENT FOR O3+hv -> O+O2
    !
    ! [O1D]ss=J[O3]/(k[H2O]+k[N2]+k[O2])
    ! SO, THE EFFECTIVE J-VALUE IS J*k[H2O]/(k[H2O]+k[N2]+k[O2])
    !
    ! We don't want to do this if strat-chem is in use, as all
    ! the intermediate reactions are included - this would be
    ! double-counting (SDE 04/01/13)
    !==============================================================

    ! Need to subtract O3->O1D from rate
    ! RXN_O3_1 specifies: O3 + hv -> O2 + O
    ! RXN_O3_2 specifies: O3 + hv -> O2 + O(1D)
    ZPJ(L,RXN_O3_1,I,J) = ZPJ(L,RXN_O3_1,I,J) &
                        - ZPJ(L,RXN_O3_2,I,J)

    ! Free pointers
    ZPJ => NULL()

  END SUBROUTINE PHOTRATE_ADJ
!EOC
!32l!------------------------------------------------------------------------------
!32l!                  GEOS-Chem Global Chemical Transport Model                  !
!32l!------------------------------------------------------------------------------
!32l!BOP
!32l!
!32l! !IROUTINE: rd_aod
!32l!
!32l! !DESCRIPTION: Subroutine RD\_AOD reads aerosol phase functions that are
!32l!  used to scale diagnostic output to an arbitrary wavelengh.  This
!32l!  facilitates comparing with satellite observations.
!32l!\\
!32l!\\
!32l! !INTERFACE:
!32l!
!32l  SUBROUTINE RD_AOD( Input_Opt, State_Chm, RC )
!32l!
!32l! !USES:
!32l!
!32l    USE ErrCode_Mod
!32l    USE Input_Opt_Mod, ONLY : OptInput
!32l    USE InquireMod,    ONLY : FindFreeLUN
!32l    USE State_Chm_Mod, ONLY : ChmState
!32l#if defined( MODEL_CESM )
!32l    USE UNITS,         ONLY : freeUnit
!32l#endif
!32l!
!32l! !INPUT PARAMETERS:
!32l!
!32l    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
!32l!
!32l! !INPUT/OUTPUT PARAMETERS:
!32l!
!32l    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!32l!
!32l! !OUTPUT PARAMETERS:
!32l!
!32l    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
!32l!
!32l! !REMARKS:
!32l!  The .dat files for each species contain the optical properties
!32l!  at multiple wavelengths to be used in the online calculation of the aerosol
!32l!  optical depth diagnostics.
!32l!  These properties have been calculated using the same size and optical
!32l!  properties as the FJX_spec.dat file used for the FAST-J photolysis
!32l!  calculations (which is now redundant for aerosols, the values in the .dat
!32l!  files here are now used). The file currently contains 11 wavelengths
!32l!  for Fast-J and other commonly used wavelengths for satellite and
!32l!  AERONET retrievals. 30 wavelengths follow that map onto RRTMG
!32l!  wavebands for radiaitive flux calculations (not used if RRTMG is off).
!32l!  A complete set of optical properties from 250-2000 nm for aerosols is
!32l!  available at:
!32l!  ftp://ftp.as.harvard.edu/geos-chem/data/aerosol_optics/hi_spectral_res
!32l!                                                                             .
!32l!     -- Colette L. Heald, 05/10/10)
!32l!     -- David A. Ridley, 05/10/13 (update for new optics files)
!32l!
!32l! !REVISION HISTORY:
!32l!  10 May 2010 - C. Heald      - Initial version
!32l!  See https://github.com/geoschem/geos-chem for complete history
!32l!EOP
!32l!------------------------------------------------------------------------------
!32l!BOC
!32l!
!32l! !LOCAL VARIABLES
!32l!
!32l    ! Scalars
!32l    INTEGER            :: I, J, K, N, g
!32l    INTEGER            :: IOS, NJ1
!32l    LOGICAL            :: LBRC, FileExists
!32l
!32l    ! Strings
!32l    CHARACTER(LEN=78 ) :: TITLE0
!32l    CHARACTER(LEN=255) :: DATA_DIR
!32l    CHARACTER(LEN=255) :: THISFILE
!32l    CHARACTER(LEN=255) :: FileMsg
!32l    CHARACTER(LEN=255) :: ErrMsg
!32l    CHARACTER(LEN=255) :: ThisLoc
!32l
!32l    ! String arrays
!32l    CHARACTER(LEN=30)  :: SPECFIL(8)
!32l
!32l    ! Pointers
!32l    REAL*8, POINTER :: WVAA  (:,:)    
!32l    REAL*8, POINTER :: RHAA  (:,:)    
!32l    REAL*8, POINTER :: RDAA  (:,:,:)    
!32l    REAL*8, POINTER :: RWAA  (:,:,:)    
!32l    REAL*8, POINTER :: SGAA  (:,:)    
!32l    REAL*8, POINTER :: REAA  (:,:,:)    
!32l    REAL*8, POINTER :: NCMAA (:,:,:)  
!32l    REAL*8, POINTER :: NRLAA (:,:,:)  
!32l    REAL*8, POINTER :: QQAA  (:,:,:,:)  
!32l    REAL*8, POINTER :: ALPHAA(:,:,:,:)  
!32l    REAL*8, POINTER :: SSAA  (:,:,:,:)  
!32l    REAL*8, POINTER :: ASYMAA(:,:,:,:)  
!32l    REAL*8, POINTER :: PHAA  (:,:,:,:,:)
!32l
!32l    !================================================================
!32l    ! RD_AOD begins here!
!32l    !================================================================
!32l
!32l    ! Initialize
!32l    RC       = GC_SUCCESS
!32l    ErrMsg   = ''
!32l    ThisLoc  = ' -> at RD_AOD (in module GeosCore/photolysis_mod.F90)'
!32l    LBRC     = Input_Opt%LBRC
!32l    DATA_DIR = TRIM( Input_Opt%AER_OPTICS_DIR )
!32l
!32l    ! Set Pointers
!32l    WVAA   => State_Chm%Phot%WVAA
!32l    RHAA   => State_Chm%Phot%RHAA
!32l    RDAA   => State_Chm%Phot%RDAA
!32l    RWAA   => State_Chm%Phot%RWAA
!32l    SGAA   => State_Chm%Phot%SGAA
!32l    REAA   => State_Chm%Phot%REAA
!32l    NRLAA  => State_Chm%Phot%NRLAA
!32l    NCMAA  => State_Chm%Phot%NCMAA
!32l    QQAA   => State_Chm%Phot%QQAA
!32l    ALPHAA => State_Chm%Phot%ALPHAA
!32l    SSAA   => State_Chm%Phot%SSAA
!32l    ASYMAA => State_Chm%Phot%ASYMAA
!32l    PHAA   => State_Chm%Phot%PHAA
!32l
!32l    ! Get a free LUN
!32l    NJ1 = findFreeLUN()
!32l
!32l    ! IMPORTANT: aerosol_mod.F and dust_mod.F expect aerosols in this order
!32l    !
!32l    ! Treating strat sulfate with GADS data but modified to match
!32l    ! the old Fast-J values size (r=0.09um, sg=0.6) - I think there's
!32l    ! evidence that this is too smale and narrow e.g. Deshler et al. 2003
!32l    ! NAT should really be associated with something like cirrus cloud
!32l    ! but for now we are just treating the NAT like the sulfate... limited
!32l    ! info but ref index is similar e.g. Scarchilli et al. (2005)
!32l    !(DAR 05/2015)
!32l    DATA SPECFIL /"so4.dat","soot.dat","org.dat", &
!32l                  "ssa.dat","ssc.dat",            &
!32l                  "h2so4.dat","h2so4.dat",        &
!32l                  "dust.dat"/
!32l
!32l    ! Loop over the array of filenames
!32l    DO k = 1, State_Chm%Phot%NSPAA
!32l
!32l       ! Choose different set of input files for standard (trop+strat chenm)
!32l       ! and tropchem (trop-only chem) simulations
!32l       THISFILE = TRIM( DATA_DIR ) // '/' // TRIM( SPECFIL(k) )
!32l
!32l       !--------------------------------------------------------------
!32l       ! In dry-run mode, print file path to dryrun log and cycle.
!32l       ! Otherwise, print file path to stdout and continue.
!32l       !--------------------------------------------------------------
!32l
!32l       ! Test if the file exists
!32l       INQUIRE( FILE=TRIM( ThisFile ), EXIST=FileExists )
!32l
!32l       ! Test if the file exists and define an output string
!32l       IF ( FileExists ) THEN
!32l          FileMsg = 'PHOTOLYSIS (RD_AOD): Opening'
!32l       ELSE
!32l          FileMsg = 'PHOTOLYSIS (RD_AOD): REQUIRED FILE NOT FOUND'
!32l       ENDIF
!32l
!32l       ! Write to stdout for both regular and dry-run simulations
!32l       IF ( Input_Opt%amIRoot ) THEN
!32l          WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
!32l300       FORMAT( a, ' ', a )
!32l       ENDIF
!32l
!32l       ! For dry-run simulations, cycle to next file.
!32l       ! For regular simulations, throw an error if we can't find the file.
!32l       IF ( Input_Opt%DryRun ) THEN
!32l          CYCLE
!32l       ELSE
!32l          IF ( .not. FileExists ) THEN
!32l             WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
!32l             CALL GC_Error( ErrMsg, RC, ThisLoc )
!32l             RETURN
!32l          ENDIF
!32l       ENDIF
!32l
!32l       !--------------------------------------------------------------
!32l       ! If not a dry-run, read data from each species file
!32l       !--------------------------------------------------------------
!32l
!32l       ! Open file
!32l       OPEN( NJ1, FILE=TRIM( THISFILE ), STATUS='OLD', IOSTAT=RC )
!32l
!32l       ! Error check
!32l       IF ( RC /= 0 ) THEN
!32l          ErrMsg = 'Error opening file: ' // TRIM( ThisFile )
!32l          CALL GC_Error( ErrMsg, RC, ThisLoc )
!32l          RETURN
!32l       ENDIF
!32l
!32l       ! Read header lines
!32l       READ(  NJ1, '(A)' ) TITLE0
!32l       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0
!32l
!32l       ! Second header line added for more info
!32l       READ(  NJ1, '(A)' ) TITLE0
!32l       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0
!32l
!32l       READ(  NJ1, '(A)' ) TITLE0
!32l110    FORMAT( 3x, a20 )
!32l
!32l       IF (k == 1 .OR. k == 3) THEN
!32l       ! for SO4 and ORGANICS, dry aerosol size varies, therefore all 
!32l       ! opt properties vary. 
!32l       DO g = 1, State_Chm%Phot%NDRg
!32l       DO i = 1, State_Chm%Phot%NRAA
!32l       DO j = 1, State_Chm%Phot%NWVAA
!32l
!32l          READ(NJ1,*) WVAA(j,k),RHAA(i,k),NRLAA(j,i,k),NCMAA(j,i,k), &
!32l                      RDAA(i,k,g),RWAA(i,k,g),SGAA(i,k),QQAA(j,i,k,g),   &
!32l                      ALPHAA(j,i,k,g),REAA(i,k,g),SSAA(j,i,k,g),         &
!32l                      ASYMAA(j,i,k,g),(PHAA(j,i,k,n,g),n=1,8)
!32l
!32l          ! make note of where 1000nm is for FAST-J calcs
!32l          IF (WVAA(j,k).EQ.1000.0) State_Chm%Phot%IWV1000=J
!32l
!32l       ENDDO
!32l       ENDDO
!32l       ENDDO
!32l
!32l       ELSE
!32l       ! For other species, keep g = default Rg (DRg) 
!32l       g = State_Chm%Phot%DRg
!32l       DO i = 1, State_Chm%Phot%NRAA
!32l       DO j = 1, State_Chm%Phot%NWVAA
!32l
!32l          READ(NJ1,*) WVAA(j,k),RHAA(i,k),NRLAA(j,i,k),NCMAA(j,i,k), &
!32l                      RDAA(i,k,g),RWAA(i,k,g),SGAA(i,k),QQAA(j,i,k,g),   &
!32l                      ALPHAA(j,i,k,g),REAA(i,k,g),SSAA(j,i,k,g),         &
!32l                      ASYMAA(j,i,k,g),(PHAA(j,i,k,n,g),n=1,8)
!32l
!32l          ! make note of where 1000nm is for FAST-J calcs
!32l          IF (WVAA(j,k).EQ.1000.0) State_Chm%Phot%IWV1000=J
!32l
!32l       ENDDO
!32l       ENDDO
!32l
!32l       ENDIF
!32l
!32l       ! Close file
!32l       CLOSE( NJ1 )
!32l
!32l    ENDDO
!32l    
!32l#if defined( MODEL_CESM )
!32l   CALL freeUnit(NJ1)
!32l#endif
!32l
!32l  ! Free pointers
!32l    WVAA   => NULL()
!32l    RHAA   => NULL()
!32l    RDAA   => NULL()
!32l    RWAA   => NULL()
!32l    SGAA   => NULL()
!32l    REAA   => NULL()
!32l    NCMAA  => NULL()
!32l    NRLAA  => NULL()
!32l    QQAA   => NULL()
!32l    ALPHAA => NULL()
!32l    SSAA   => NULL()
!32l    ASYMAA => NULL()
!32l    PHAA   => NULL()
!32l
!32l  END SUBROUTINE RD_AOD
!32l!EOC
!32l!------------------------------------------------------------------------------
!32l!                  GEOS-Chem Global Chemical Transport Model                  !
!32l!------------------------------------------------------------------------------
!32l!BOP
!32l!
!32l! !IROUTINE: calc_aod
!32l!
!32l! !DESCRIPTION: Subroutine CALC\_AOD works out the closest tie points
!32l! in the optics LUT wavelengths and the coefficients required to
!32l! calculate the angstrom exponent for interpolating optics to the requested
!32l! wavelength. If the wavelength requested matches a standard wavelength
!32l! in the LUT then we skip the interpolation (DAR 09/2013)
!32l!\\
!32l!\\
!32l! !INTERFACE:
!32l!
!32l  SUBROUTINE CALC_AOD( Input_Opt, State_Chm, RC )
!32l!
!32l! !USES:
!32l!
!32l    USE Input_Opt_Mod, ONLY : OptInput
!32l#ifdef RRTMG
!32l    USE PARRRTM,       ONLY : NBNDLW
!32l#endif
!32l    USE State_Chm_Mod, ONLY : ChmState
!32l!
!32l! !INPUT PARAMETERS:
!32l!
!32l    TYPE(OptInput), INTENT(IN)    :: Input_Opt
!32l!
!32l! !INPUT/OUTPUT PARAMETERS:
!32l!
!32l    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!32l!
!32l! !OUTPUT PARAMETERS:
!32l!
!32l    INTEGER,        INTENT(IN)    :: RC
!32l!
!32l! !REMARKS:
!32l!  Now the user is able to select any 3 wavelengths for optics
!32l!  output in the geoschem_config.yml file we need to be able to interpolate
!32l!  to those wavelengths based on what is available in the optics
!32l!  look-up table.
!32l!                                                                             .
!32l!  The standard lookup table currently has values for
!32l!  11 common wavelengths followed by 30 that are required by RRTMG.
!32l!  Only those required to interpolate to user requested
!32l!  wavelengths are selected from the standard wavelengths. RRTMG
!32l!  wavelengths are not used in the interpolation for AOD output
!32l!  (DAR 10/2013)
!32l!                                                                             .
!32l!   UPDATE: because the RT optics output doesnt have access to the
!32l!   standard wavelengths we now calculate two sets of values: one
!32l!   for the ND21 and diag3 outputs that use the standard wavelengths
!32l!   and one for RRTMG diagnostics that interpolate the optics from RRTMG
!32l!   wavelengths. Perhaps a switch needs adding to switch off the RT
!32l!   optics output (and interpolation) if this ends up costing too
!32l!   much and is not used, but it is ideal to have an optics output
!32l!   that matches exactly what RRTMG uses to calculate the fluxes
!32l!
!32l! !REVISION HISTORY:
!32l!  18 Jun 2013 - D. Ridley   - Initial version
!32l!  See https://github.com/geoschem/geos-chem for complete history
!32l!EOP
!32l!------------------------------------------------------------------------------
!32l!BOC
!32l!
!32l! !LOCAL VARIABLES
!32l!
!32l    INTEGER             :: MINWV, MAXWV, N, N0, N1, W, NSTEP
!32l    INTEGER             :: NWVAA, NWVAA0, NWVREQUIRED, NRTWVREQUIRED
!32l    REAL(fp)            :: WVDIF
!32l
!32l    ! Pointers
!32l    INTEGER, POINTER :: IWVREQUIRED  (:)    
!32l    INTEGER, POINTER :: IRTWVREQUIRED(:)
!32l    INTEGER, POINTER :: IWVSELECT    (:,:)
!32l    INTEGER, POINTER :: IRTWVSELECT  (:,:)
!32l    REAL*8,  POINTER :: ACOEF_WV     (:)
!32l    REAL*8,  POINTER :: BCOEF_WV     (:)
!32l    REAL*8,  POINTER :: CCOEF_WV     (:)
!32l    REAL*8,  POINTER :: ACOEF_RTWV   (:)
!32l    REAL*8,  POINTER :: BCOEF_RTWV   (:)
!32l    REAL*8,  POINTER :: CCOEF_RTWV   (:)
!32l    REAL*8,  POINTER :: WVAA         (:,:)
!32l
!32l    !================================================================
!32l    ! CALC_AOD begins here!
!32l    !================================================================
!32l
!32l    ! Constants State_Chm%Phot
!32l    NWVAA         = State_Chm%Phot%NWVAA
!32l    NWVAA0        = State_Chm%Phot%NWVAA0
!32l
!32l    ! Scalars in State_Chm%Phot that will be set in this subroutine
!32l    NWVREQUIRED   = State_Chm%Phot%NWVREQUIRED
!32l    NRTWVREQUIRED = State_Chm%Phot%NRTWVREQUIRED
!32l
!32l    ! Set pointers
!32l    IWVREQUIRED   => State_Chm%Phot%IWVREQUIRED
!32l    IRTWVREQUIRED => State_Chm%Phot%IRTWVREQUIRED
!32l    IWVSELECT     => State_Chm%Phot%IWVSELECT
!32l    IRTWVSELECT   => State_Chm%Phot%IRTWVSELECT
!32l    ACOEF_WV      => State_Chm%Phot%ACOEF_WV
!32l    BCOEF_WV      => State_Chm%Phot%BCOEF_WV
!32l    CCOEF_WV      => State_Chm%Phot%CCOEF_WV
!32l    ACOEF_RTWV    => State_Chm%Phot%ACOEF_RTWV
!32l    BCOEF_RTWV    => State_Chm%Phot%BCOEF_RTWV
!32l    CCOEF_RTWV    => State_Chm%Phot%CCOEF_RTWV
!32l    WVAA          => State_Chm%Phot%WVAA
!32l
!32l    !cycle over standard wavelengths
!32l    N0=1
!32l    N1=NWVAA0
!32l    NSTEP=1
!32l    NWVREQUIRED=0
!32l    DO W=1,Input_Opt%NWVSELECT
!32l       MINWV     = -999
!32l       MAXWV     =  999
!32l       DO N=N0,N1,NSTEP ! 1 to 11
!32l          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
!32l          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
!32l             MINWV = WVDIF
!32l             IWVSELECT(1,W)=N             
!32l          ENDIF
!32l          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
!32l             MAXWV = WVDIF
!32l             IWVSELECT(2,W)=N
!32l          ENDIF
!32l       ENDDO
!32l       IF (IWVSELECT(2,W).EQ.IWVSELECT(1,W)) THEN
!32l          !we have a match!
!32l          MINWV=0
!32l          MAXWV=0
!32l          !add this wavelength to those for output
!32l          NWVREQUIRED=NWVREQUIRED+1
!32l          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
!32l       ELSE
!32l          !we are going to have to interpolate to the requested wavelength
!32l          NWVREQUIRED=NWVREQUIRED+1
!32l          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
!32l          NWVREQUIRED=NWVREQUIRED+1
!32l          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(2,W)
!32l       ENDIF
!32l
!32l       !Error check - ensure we have a match or requested wavelength
!32l       !falls within two LUT tie points
!32l       IF (MINWV.EQ.-999) THEN
!32l          ! requested wavelength is shorter than min wv in LUT
!32l          ! set to min
!32l          write(6,*) 'ERROR requested wavelength is too short!!'
!32l          write(6,*) 'Defaulting to LUT min: ',WVAA(1,1)
!32l          IWVSELECT(1,W)=1
!32l          IWVSELECT(2,W)=1 !300nm
!32l          NWVREQUIRED=NWVREQUIRED-1
!32l          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
!32l       ENDIF
!32l       IF (MAXWV.EQ.999) THEN
!32l          ! requested wavelength is longer than min wv in LUT
!32l          ! set to max
!32l          write(6,*) 'ERROR requested wavelength is too long!!'
!32l          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0,1)
!32l          IWVSELECT(1,W)=NWVAA0
!32l          IWVSELECT(2,W)=NWVAA0 !1020nm
!32l          NWVREQUIRED=NWVREQUIRED-1
!32l          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
!32l       ENDIF
!32l
!32l       !now calcualte the angstrom exponent coefs for interpolation -
!32l       !this is done here to save time and repetition in aerosol_mod.F
!32l       IF (IWVSELECT(1,W).NE.IWVSELECT(2,W)) THEN
!32l          ACOEF_WV(W) = WVAA(IWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
!32l          BCOEF_WV(W) =1.0d0/(LOG(WVAA(IWVSELECT(2,W),1)/ &
!32l                                  WVAA(IWVSELECT(1,W),1)))
!32l          !relative location of selected wavelength between tie points
!32l          !for interpolating SSA and ASYM for output in aerosol_mod.F and
!32l          !dust_mod.F
!32l          CCOEF_WV(W) =(Input_Opt%WVSELECT(W)-WVAA(IWVSELECT(1,W),1))/ &
!32l                      (WVAA(IWVSELECT(2,W),1)-WVAA(IWVSELECT(1,W),1))
!32l       ENDIF
!32l       IF ( Input_Opt%amIRoot ) THEN
!32l          write(6,*) 'N WAVELENGTHS: ',Input_Opt%NWVSELECT
!32l          write(6,*) 'WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
!32l          write(6,*) 'WAVELENGTH REQUIRED:', NWVREQUIRED
!32l          !write(6,*) IWVSELECT(1,W),WVAA(IWVSELECT(1,W),1)
!32l          !write(6,*) IWVSELECT(2,W),WVAA(IWVSELECT(2,W),1)
!32l          !write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
!32l          write(6,*) '*********************************'
!32l       ENDIF
!32l    ENDDO !Input_Opt%NWVSELECT
!32l#ifdef RRTMG
!32l    !repeat for RRTMG wavelengths to get the closest wavelength
!32l    !indices and the interpolation coefficients
!32l    !Indices are relative to all wavelengths in the LUT i.e. the RRTMG
!32l    !wavelengths start at NWVAA0+1
!32l    N0=NWVAA0+1
!32l    N1=NWVAA
!32l    NSTEP=1
!32l    NRTWVREQUIRED=0
!32l    DO W=1,Input_Opt%NWVSELECT
!32l       MINWV     = -999
!32l       MAXWV     =  999
!32l       DO N=N0,N1,NSTEP
!32l          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
!32l          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
!32l             MINWV = WVDIF
!32l             IRTWVSELECT(1,W)=N
!32l          ENDIF
!32l          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
!32l             MAXWV = WVDIF
!32l             IRTWVSELECT(2,W)=N
!32l          ENDIF
!32l       ENDDO
!32l       IF (IRTWVSELECT(2,W).EQ.IRTWVSELECT(1,W)) THEN
!32l          !we have a match!
!32l          MINWV=0
!32l          MAXWV=0
!32l          !add this wavelength to those for output
!32l          NRTWVREQUIRED=NRTWVREQUIRED+1
!32l          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
!32l       ELSE
!32l          !we are going to have to interpolate to the requested
!32l          !wavelength
!32l          NRTWVREQUIRED=NRTWVREQUIRED+1
!32l          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
!32l          NRTWVREQUIRED=NRTWVREQUIRED+1
!32l          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(2,W)
!32l       ENDIF
!32l
!32l       !Error check - ensure we have a match or requested wavelength
!32l       !falls within two LUT tie points
!32l       IF (MINWV.EQ.-999) THEN
!32l          ! requested wavelength is shorter than min wv in LUT
!32l          ! set to min
!32l          write(6,*) 'ERROR requested wavelength is too short!!'
!32l          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA-1,1)
!32l          IRTWVSELECT(1,W)=NWVAA-1
!32l          IRTWVSELECT(2,W)=NWVAA-1
!32l          NRTWVREQUIRED=NRTWVREQUIRED-1
!32l          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
!32l       ENDIF
!32l       IF (MAXWV.EQ.999) THEN
!32l          ! requested wavelength is longer than min wv in LUT
!32l          ! set to max
!32l          write(6,*) 'ERROR requested wavelength is too long!!'
!32l          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0+1,1)
!32l          IRTWVSELECT(1,W)=NWVAA0+1
!32l          IRTWVSELECT(2,W)=NWVAA0+1
!32l          NRTWVREQUIRED=NRTWVREQUIRED-1
!32l          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
!32l       ENDIF
!32l
!32l       !now calcualte the angstrom exponent coefs for interpolation -
!32l       !this is done here to save time and repetition in aerosol_mod.F
!32l       IF (IRTWVSELECT(1,W).NE.IRTWVSELECT(2,W)) THEN
!32l          ACOEF_RTWV(W) = WVAA(IRTWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
!32l          BCOEF_RTWV(W) =1.0d0/(LOG(WVAA(IRTWVSELECT(2,W),1)/ &
!32l                                    WVAA(IRTWVSELECT(1,W),1)))
!32l          !relative location of selected wavelength between tie points
!32l          !for interpolating SSA and ASYM for output in aerosol_mod.F and
!32l          !dust_mod.F
!32l          CCOEF_RTWV(W) =(Input_Opt%WVSELECT(W)-WVAA(IRTWVSELECT(1,W),1))/ &
!32l                      (WVAA(IRTWVSELECT(2,W),1)-WVAA(IRTWVSELECT(1,W),1))
!32l       ENDIF
!32l       !convert wavelength index to that required by rrtmg_rad_transfer
!32l       !i.e. without the standard and LW wavelengths
!32l       IRTWVSELECT(1,W) = IRTWVSELECT(1,W) - NWVAA0 - NBNDLW
!32l       IRTWVSELECT(2,W) = IRTWVSELECT(2,W) - NWVAA0 - NBNDLW
!32l       IF ( Input_Opt%amIRoot ) THEN
!32l          write(6,*) 'N RT WAVELENGTHS: ',Input_Opt%NWVSELECT
!32l          write(6,*) 'RT WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
!32l          write(6,*) 'RT WAVELENGTH REQUIRED:', NRTWVREQUIRED
!32l          write(6,*) IRTWVSELECT(1,W),WVAA(IRTWVSELECT(1,W)+NWVAA0+NBNDLW,1)
!32l          write(6,*) IRTWVSELECT(2,W),WVAA(IRTWVSELECT(2,W)+NWVAA0+NBNDLW,1)
!32l          write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
!32l          write(6,*) '*********************************'
!32l       ENDIF
!32l    ENDDO !Input_Opt%NWVSELECT
!32l#endif
!32l
!32l    ! Copy values back into State_Chm
!32l    State_Chm%Phot%NWVREQUIRED   = NWVREQUIRED
!32l    State_Chm%Phot%NRTWVREQUIRED = NRTWVREQUIRED
!32l
!32l    ! Free pointers
!32l    IWVREQUIRED   => NULL() 
!32l    IRTWVREQUIRED => NULL()
!32l    IWVSELECT     => NULL()
!32l    IRTWVSELECT   => NULL()
!32l    ACOEF_WV      => NULL()
!32l    BCOEF_WV      => NULL()
!32l    CCOEF_WV      => NULL()
!32l    ACOEF_RTWV    => NULL()
!32l    BCOEF_RTWV    => NULL()
!32l    CCOEF_RTWV    => NULL()
!32l    WVAA          => NULL()
!32l
!32l  END SUBROUTINE CALC_AOD
!32l!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_aer
!
! !DESCRIPTION: Subroutine SET\_AER fills out the array MIEDX.
!  Each entry connects a GEOS-Chem aerosol to its Fast-JX counterpart:
!  MIEDX(Fast-JX index) = (GC index)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_AER( Input_Opt, State_Chm, RC )
!
! !USES:
!
#ifdef FASTJX
    USE CMN_FJX_Mod,    ONLY : AN_, NAA, TITLAA
#else
    USE Cldj_Cmn_Mod,   ONLY : AN_, NAA, TITLAA
#endif
    USE CMN_SIZE_Mod,   ONLY : NRHAER, NRH
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input options
!
! !INPUT/OUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  31 Mar 2013 - S. D. Eastham - Adapted from J. Mao FJX v6.2 implementation
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    INTEGER            :: I, J, K
    INTEGER            :: IND(NRHAER)
    INTEGER,   POINTER :: MIEDX(:)

    !=================================================================
    ! SER_AER begins here!
    !=================================================================

    ! Initialize
    RC = GC_SUCCESS
    ErrMsg = ''
    ThisLoc = ' -> at Set_Aer (in module GeosCore/photolysis_mod.F90)'


    ! Set pointer
    MIEDX => State_Chm%Phot%MIEDX

    ! Taken from aerosol_mod.F
    IND = (/22,29,36,43,50/)

    DO I=1,AN_
       MIEDX(I) = 0
    ENDDO

    ! Select Aerosol/Cloud types to be used - define types here
    ! Each of these types must be listed in the order used by OPMIE.F

    ! Clouds
    MIEDX(1)  =  3   !  Black carbon absorber
    MIEDX(2)  = 10   !  Water Cloud (Deirmenjian 8 micron)
    MIEDX(3)  = 14   !  Irregular Ice Cloud (Mishchenko)

    ! Dust
    MIEDX(4)  = 15   !  Mineral Dust  .15 micron    (rvm, 9/30/00)
    MIEDX(5)  = 16   !  Mineral Dust  .25 micron    (rvm, 9/30/00)
    MIEDX(6)  = 17   !  Mineral Dust  .4  micron    (rvm, 9/30/00)
    MIEDX(7)  = 18   !  Mineral Dust  .8  micron    (rvm, 9/30/00)
    MIEDX(8)  = 19   !  Mineral Dust 1.5  micron    (rvm, 9/30/00)
    MIEDX(9)  = 20   !  Mineral Dust 2.5  micron    (rvm, 9/30/00)
    MIEDX(10) = 21   !  Mineral Dust 4.0  micron    (rvm, 9/30/00)

    ! Aerosols
    DO I=1,NRHAER
       DO J=1,NRH
          MIEDX(10+((I-1)*NRH)+J)=IND(I)+J-1
       ENDDO
    ENDDO

    ! Stratospheric aerosols - SSA/STS and solid PSCs
#ifdef FASTJX
    MIEDX(10+(NRHAER*NRH)+1) = 4  ! SSA/LBS/STS
#else
    MIEDX(10+(NRHAER*NRH)+1) = 1  ! SSA/LBS/STS
#endif
    MIEDX(10+(NRHAER*NRH)+2) = 14 ! NAT/ice PSCs

    ! Ensure all 'AN_' types are valid selections
    do i=1,AN_
       IF (Input_Opt%amIRoot) write(6,1000) MIEDX(i),TITLAA(MIEDX(i))
       if (MIEDX(i).gt.NAA.or.MIEDX(i).le.0) then
          if (Input_Opt%amIRoot) then
             write(6,1200) MIEDX(i),NAA
          endif
          ErrMsg = 'Bad MIEDX value in "Set_AER"!'
          call GC_Error( ErrMsg, RC, ThisLoc )
          return
       endif
    enddo

    ! Free pointer
    MIEDX => NULL()

1000 format('Using Aerosol type: ',i3,1x,a)
1200 format('Aerosol type ',i3,' unsuitable; supplied values must be ', &
            'between 1 and ',i3)

  END SUBROUTINE SET_AER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_prof_nc
!
! !DESCRIPTION: Subroutine RD\_PROF\_NC reads in the reference climatology
!  from a NetCDF file rather than an ASCII .dat.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_PROF_NC( Input_Opt, State_Grid, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Chm_Mod,  ONLY : ChmState

#if defined( MODEL_CESM )
    USE CAM_PIO_UTILS,     ONLY : CAM_PIO_OPENFILE
    USE IOFILEMOD,         ONLY : GETFIL
    USE PIO,               ONLY : PIO_CLOSEFILE
    USE PIO,               ONLY : PIO_INQ_DIMID
    USE PIO,               ONLY : PIO_INQ_DIMLEN
    USE PIO,               ONLY : PIO_INQ_VARID
    USE PIO,               ONLY : PIO_GET_VAR
    USE PIO,               ONLY : PIO_NOERR
    USE PIO,               ONLY : PIO_NOWRITE
    USE PIO,               ONLY : FILE_DESC_T
#else
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This file was automatically generated by the Perl scripts in the
!  NcdfUtilities package (which ships w/ GEOS-Chem) and was subsequently
!  hand-edited.
!
! !REVISION HISTORY:
!  19 Apr 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists          ! Does input file exist?
    INTEGER            :: fId                 ! netCDF file ID

    ! Strings
    CHARACTER(LEN=255) :: nc_dir              ! netCDF directory name
    CHARACTER(LEN=255) :: nc_file             ! netCDF file name
    CHARACTER(LEN=255) :: nc_path             ! netCDF path name
    CHARACTER(LEN=255) :: v_name              ! netCDF variable name
    CHARACTER(LEN=255) :: a_name              ! netCDF attribute name
    CHARACTER(LEN=255) :: a_val               ! netCDF attribute value
    CHARACTER(LEN=255) :: FileMsg
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! Arrays
    INTEGER            :: st3d(3), ct3d(3)    ! For 3D arrays

#if defined( MODEL_CESM )
    type(FILE_DESC_T)  :: ncid
    INTEGER            :: vId, iret
#endif

    ! Pointers
    REAL(fp), POINTER :: OREF(:,:,:)
    REAL(fp), POINTER :: TREF(:,:,:)

    !=================================================================
    ! RD_PROF_NC begins here!
    !=================================================================

    ! Initialize
    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at RD_PROF_NC (in module GeosCore/photolysis_mod.F90)'

#if defined( MODEL_CESM )
    ! In the CESM model, only read on the root chunk, but on all CPUs (hplin, 7/3/24)
    IF ( State_Grid%CPU_Subdomain_ID .ne. State_Grid%CPU_Subdomain_FirstID ) RETURN
#endif

    ! Set pointers
    OREF => State_Chm%Phot%OREF
    TREF => State_Chm%Phot%TREF

    ! Directory and file names
    nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // '/' // 'FastJ_201204' // '/'
    nc_file = 'fastj.jv_atms_dat.nc'
    nc_path = TRIM( nc_dir ) // TRIM( nc_file )

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Test if the file exists
    INQUIRE( FILE=TRIM( nc_path ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'PHOTOLYSIS (RD_PROF_NC): Opening'
    ELSE
       FileMsg = 'PHOTOLYSIS (RD_PROF_NC): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( nc_path )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( Input_Opt%DryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( nc_path )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=========================================================================
    ! Open and read data from the netCDF file
    !=========================================================================

    ! Open netCDF file
#if defined( MODEL_CESM )
    ! Note: In CESM environment, PIO_OPENFILE is a collective operation and must
    ! be called by all CPUs. (hplin, 7/3/24)
    CALL CAM_PIO_OPENFILE( ncid, TRIM(nc_path), PIO_NOWRITE )
#else
    CALL Ncop_Rd( fId, TRIM(nc_path) )
#endif

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) REPEAT( '%', 79 )
       WRITE( 6, 110 ) TRIM(nc_file)
       WRITE( 6, 120 ) TRIM(nc_dir)
    ENDIF

    !----------------------------------------
    ! VARIABLE: T
    !----------------------------------------

    ! Variable name
    v_name = "T"

    ! Read T from file
    st3d   = (/  1,  1,  1 /)
    ct3d   = (/ 51, 18, 12 /)
#if defined( MODEL_CESM )
    iret = PIO_INQ_VARID( ncid, trim(v_name), vid  )
    iret = PIO_GET_VAR( ncid, vid, st3d, ct3d, TREF )
#else
    CALL NcRd( TREF, fId, TRIM(v_name), st3d, ct3d )

    ! Read the T:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF
#endif

    !----------------------------------------
    ! VARIABLE: O3
    !----------------------------------------

    ! Variable name
    v_name = "O3"

    ! Read O3 from file
    st3d   = (/  1,  1,  1 /)
    ct3d   = (/ 51, 18, 12 /)
#if defined( MODEL_CESM )
    iret = PIO_INQ_VARID( ncid, trim(v_name), vid  )
    iret = PIO_GET_VAR( ncid, vid, st3d, ct3d, OREF )
#else
    CALL NcRd( OREF, fId, TRIM(v_name), st3d, ct3d )

    ! Read the O3:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF
#endif

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Close netCDF file
#if defined( MODEL_CESM )
    CALL PIO_CLOSEFILE( ncid )
#else
    CALL NcCl( fId )
#endif

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 140 )
       WRITE( 6, 100 ) REPEAT( '%', 79 )
    ENDIF

    ! Free pointers
    OREF => NULL()
    TREF => NULL()

    ! FORMAT statements
100 FORMAT( a                                              )
110 FORMAT( '%% Opening file  : ',         a               )
120 FORMAT( '%%  in directory : ',         a, / , '%%'     )
130 FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
140 FORMAT( '%% Successfully closed file!'                 )

  END SUBROUTINE RD_PROF_NC
!EOC

END MODULE PHOTOLYSIS_MOD
