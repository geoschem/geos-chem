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
    ! Exit if photolysis disabled (zero J-values)
    !------------------------------------------------------------------------
    IF ( Input_Opt%DryRun .OR. .NOT. Input_Opt%Do_Photolysis ) RETURN

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
