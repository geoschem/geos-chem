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
  PRIVATE :: RD_AOD
  PRIVATE :: CALC_AOD
  PRIVATE :: SET_AER
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
  SUBROUTINE INIT_PHOTOLYSIS( Input_Opt, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE Charpak_Mod,    ONLY : CSTRIP
    USE CMN_FJX_Mod,    ONLY : JVN_, W_, JLABEL, RNAMES, WL, JFACTA
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : Planck, CConst
    USE State_Chm_Mod,  ONLY : ChmState, Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE Fjx_Interface_Mod,  ONLY : Init_FastJX
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
    ENDIF

    !------------------------------------------------------------------------
    ! Read in AOD data even if photolysis disabled
    ! (or just print file name if in dry-run mode)
    !------------------------------------------------------------------------
    CALL RD_AOD( Input_Opt, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_AOD"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Exit without doing any computations if we are doing a dry-run
    !------------------------------------------------------------------------
    IF ( Input_Opt%DryRun ) RETURN

    !------------------------------------------------------------------------
    ! Compute the required wavelengths in the LUT to calculate requested AOD
    !------------------------------------------------------------------------
    IF (Input_Opt%amIRoot) WRITE(6,*) 'Wavelength optics read successfully'
    CALL CALC_AOD( Input_Opt, State_Chm, RC )

    !------------------------------------------------------------------------
    ! Exit if photolysis disabled (zero J-values)
    !------------------------------------------------------------------------
    IF ( .NOT. Input_Opt%Do_Photolysis ) RETURN

    !------------------------------------------------------------------------
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
    USE Fjx_Interface_Mod,  ONLY : Run_FastJX

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

    CALL Run_FastJX( Wavelength, Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) THEN
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
!  access to the adjusted rates.  Only the bpch diagnostics are updated
!  here.
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
! !IROUTINE: rd_aod
!
! !DESCRIPTION: Subroutine RD\_AOD reads aerosol phase functions that are
!  used to scale diagnostic output to an arbitrary wavelengh.  This
!  facilitates comparing with satellite observations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_AOD( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE InquireMod,    ONLY : FindFreeLUN
    USE State_Chm_Mod, ONLY : ChmState
#if defined( MODEL_CESM )
    USE UNITS,         ONLY : freeUnit
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  The .dat files for each species contain the optical properties
!  at multiple wavelengths to be used in the online calculation of the aerosol
!  optical depth diagnostics.
!  These properties have been calculated using the same size and optical
!  properties as the FJX_spec.dat file used for the FAST-J photolysis
!  calculations (which is now redundant for aerosols, the values in the .dat
!  files here are now used). The file currently contains 11 wavelengths
!  for Fast-J and other commonly used wavelengths for satellite and
!  AERONET retrievals. 30 wavelengths follow that map onto RRTMG
!  wavebands for radiaitive flux calculations (not used if RRTMG is off).
!  A complete set of optical properties from 250-2000 nm for aerosols is
!  available at:
!  ftp://ftp.as.harvard.edu/geos-chem/data/aerosol_optics/hi_spectral_res
!                                                                             .
!     -- Colette L. Heald, 05/10/10)
!     -- David A. Ridley, 05/10/13 (update for new optics files)
!
! !REVISION HISTORY:
!  10 May 2010 - C. Heald      - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Scalars
    INTEGER            :: I, J, K, N
    INTEGER            :: IOS, NJ1
    LOGICAL            :: LBRC, FileExists

    ! Strings
    CHARACTER(LEN=78 ) :: TITLE0
    CHARACTER(LEN=255) :: DATA_DIR
    CHARACTER(LEN=255) :: THISFILE
    CHARACTER(LEN=255) :: FileMsg
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! String arrays
    CHARACTER(LEN=30)  :: SPECFIL(8)

    ! Pointers
    REAL*8, POINTER :: WVAA  (:,:)    
    REAL*8, POINTER :: RHAA  (:,:)    
    REAL*8, POINTER :: RDAA  (:,:)    
    REAL*8, POINTER :: RWAA  (:,:)    
    REAL*8, POINTER :: SGAA  (:,:)    
    REAL*8, POINTER :: REAA  (:,:)    
    REAL*8, POINTER :: NCMAA (:,:,:)  
    REAL*8, POINTER :: NRLAA (:,:,:)  
    REAL*8, POINTER :: QQAA  (:,:,:)  
    REAL*8, POINTER :: ALPHAA(:,:,:)  
    REAL*8, POINTER :: SSAA  (:,:,:)  
    REAL*8, POINTER :: ASYMAA(:,:,:)  
    REAL*8, POINTER :: PHAA  (:,:,:,:)

    !================================================================
    ! RD_AOD begins here!
    !================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at RD_AOD (in module GeosCore/fast_jx_mod.F90)'
    LBRC     = Input_Opt%LBRC
    DATA_DIR = TRIM( Input_Opt%FAST_JX_DIR )

    ! Set Pointers
    WVAA   => State_Chm%Phot%WVAA
    RHAA   => State_Chm%Phot%RHAA
    RDAA   => State_Chm%Phot%RDAA
    RWAA   => State_Chm%Phot%RWAA
    SGAA   => State_Chm%Phot%SGAA
    REAA   => State_Chm%Phot%REAA
    NRLAA  => State_Chm%Phot%NRLAA
    NCMAA  => State_Chm%Phot%NCMAA
    QQAA   => State_Chm%Phot%QQAA
    ALPHAA => State_Chm%Phot%ALPHAA
    SSAA   => State_Chm%Phot%SSAA
    ASYMAA => State_Chm%Phot%ASYMAA
    PHAA   => State_Chm%Phot%PHAA

#if defined( MODEL_CESM )
    IF ( Input_Opt%amIRoot ) THEN
       NJ1 = findFreeLUN()
    ELSE
       NJ1 = 0
    ENDIF
#else
    ! Get a free LUN
    NJ1 = findFreeLUN()
#endif

    ! IMPORTANT: aerosol_mod.F and dust_mod.F expect aerosols in this order
    !
    ! Treating strat sulfate with GADS data but modified to match
    ! the old Fast-J values size (r=0.09um, sg=0.6) - I think there's
    ! evidence that this is too smale and narrow e.g. Deshler et al. 2003
    ! NAT should really be associated with something like cirrus cloud
    ! but for now we are just treating the NAT like the sulfate... limited
    ! info but ref index is similar e.g. Scarchilli et al. (2005)
    !(DAR 05/2015)
    DATA SPECFIL /"so4.dat","soot.dat","org.dat", &
                  "ssa.dat","ssc.dat",            &
                  "h2so4.dat","h2so4.dat",        &
                  "dust.dat"/

    ! Loop over the array of filenames
    DO k = 1, State_Chm%Phot%NSPAA

       ! Choose different set of input files for standard (trop+strat chenm)
       ! and tropchem (trop-only chem) simulations
       THISFILE = TRIM( DATA_DIR ) // TRIM( SPECFIL(k) )

       !--------------------------------------------------------------
       ! In dry-run mode, print file path to dryrun log and cycle.
       ! Otherwise, print file path to stdout and continue.
       !--------------------------------------------------------------

       ! Test if the file exists
       INQUIRE( FILE=TRIM( ThisFile ), EXIST=FileExists )

       ! Test if the file exists and define an output string
       IF ( FileExists ) THEN
          FileMsg = 'FAST-JX (RD_AOD): Opening'
       ELSE
          FileMsg = 'FAST-JX (RD_AOD): REQUIRED FILE NOT FOUND'
       ENDIF

       ! Write to stdout for both regular and dry-run simulations
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
300       FORMAT( a, ' ', a )
       ENDIF

       ! For dry-run simulations, cycle to next file.
       ! For regular simulations, throw an error if we can't find the file.
       IF ( Input_Opt%DryRun ) THEN
          CYCLE
       ELSE
          IF ( .not. FileExists ) THEN
             WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       !--------------------------------------------------------------
       ! If not a dry-run, read data from each species file
       !--------------------------------------------------------------

#if defined( MODEL_CESM )
       ! Only read file on root thread if using CESM
       IF ( Input_Opt%amIRoot ) THEN
#endif

       ! Open file
       OPEN( NJ1, FILE=TRIM( THISFILE ), STATUS='OLD', IOSTAT=RC )

       ! Error check
       IF ( RC /= 0 ) THEN
          ErrMsg = 'Error opening file: ' // TRIM( ThisFile )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Read header lines
       READ(  NJ1, '(A)' ) TITLE0
       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0

       ! Second header line added for more info
       READ(  NJ1, '(A)' ) TITLE0
       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0

       READ(  NJ1, '(A)' ) TITLE0
110    FORMAT( 3x, a20 )

       DO i = 1, State_Chm%Phot%NRAA
       DO j = 1, State_Chm%Phot%NWVAA

          READ(NJ1,*) WVAA(j,k),RHAA(i,k),NRLAA(j,i,k),NCMAA(j,i,k), &
                      RDAA(i,k),RWAA(i,k),SGAA(i,k),QQAA(j,i,k),   &
                      ALPHAA(j,i,k),REAA(i,k),SSAA(j,i,k),         &
                      ASYMAA(j,i,k),(PHAA(j,i,k,n),n=1,8)

          ! make note of where 1000nm is for FAST-J calcs
          IF (WVAA(j,k).EQ.1000.0) State_Chm%Phot%IWV1000=J

       ENDDO
       ENDDO

       ! Close file
       CLOSE( NJ1 )

#if defined( MODEL_CESM )
       ENDIF
#endif

    ENDDO
    
#if defined( MODEL_CESM )
    IF ( Input_Opt%amIRoot ) CALL freeUnit(NJ1)
#endif

  ! Free pointers
    WVAA   => NULL()
    RHAA   => NULL()
    RDAA   => NULL()
    RWAA   => NULL()
    SGAA   => NULL()
    REAA   => NULL()
    NCMAA  => NULL()
    NRLAA  => NULL()
    QQAA   => NULL()
    ALPHAA => NULL()
    SSAA   => NULL()
    ASYMAA => NULL()
    PHAA   => NULL()

  END SUBROUTINE RD_AOD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_aod
!
! !DESCRIPTION: Subroutine CALC\_AOD works out the closest tie points
! in the optics LUT wavelengths and the coefficients required to
! calculate the angstrom exponent for interpolating optics to the requested
! wavelength. If the wavelength requested matches a standard wavelength
! in the LUT then we skip the interpolation (DAR 09/2013)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_AOD( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
#ifdef RRTMG
    USE PARRRTM,       ONLY : NBNDLW
#endif
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: RC
!
! !REMARKS:
!  Now the user is able to select any 3 wavelengths for optics
!  output in the geoschem_config.yml file we need to be able to interpolate
!  to those wavelengths based on what is available in the optics
!  look-up table.
!                                                                             .
!  The standard lookup table currently has values for
!  11 common wavelengths followed by 30 that are required by RRTMG.
!  Only those required to interpolate to user requested
!  wavelengths are selected from the standard wavelengths. RRTMG
!  wavelengths are not used in the interpolation for AOD output
!  (DAR 10/2013)
!                                                                             .
!   UPDATE: because the RT optics output doesnt have access to the
!   standard wavelengths we now calculate two sets of values: one
!   for the ND21 and diag3 outputs that use the standard wavelengths
!   and one for RRTMG diagnostics that interpolate the optics from RRTMG
!   wavelengths. Perhaps a switch needs adding to switch off the RT
!   optics output (and interpolation) if this ends up costing too
!   much and is not used, but it is ideal to have an optics output
!   that matches exactly what RRTMG uses to calculate the fluxes
!
! !REVISION HISTORY:
!  18 Jun 2013 - D. Ridley   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER             :: MINWV, MAXWV, N, N0, N1, W, NSTEP
    INTEGER             :: NWVAA, NWVAA0, NWVREQUIRED, NRTWVREQUIRED
    REAL(fp)            :: WVDIF

    ! Pointers
    INTEGER, POINTER :: IWVREQUIRED  (:)    
    INTEGER, POINTER :: IRTWVREQUIRED(:)
    INTEGER, POINTER :: IWVSELECT    (:,:)
    INTEGER, POINTER :: IRTWVSELECT  (:,:)
    REAL*8,  POINTER :: ACOEF_WV     (:)
    REAL*8,  POINTER :: BCOEF_WV     (:)
    REAL*8,  POINTER :: CCOEF_WV     (:)
    REAL*8,  POINTER :: ACOEF_RTWV   (:)
    REAL*8,  POINTER :: BCOEF_RTWV   (:)
    REAL*8,  POINTER :: CCOEF_RTWV   (:)
    REAL*8,  POINTER :: WVAA         (:,:)

    !================================================================
    ! CALC_AOD begins here!
    !================================================================

    ! Constants State_Chm%Phot
    NWVAA         = State_Chm%Phot%NWVAA
    NWVAA0        = State_Chm%Phot%NWVAA0

    ! Scalars in State_Chm%Phot that will be set in this subroutine
    NWVREQUIRED   = State_Chm%Phot%NWVREQUIRED
    NRTWVREQUIRED = State_Chm%Phot%NRTWVREQUIRED

    ! Set pointers
    IWVREQUIRED   => State_Chm%Phot%IWVREQUIRED
    IRTWVREQUIRED => State_Chm%Phot%IRTWVREQUIRED
    IWVSELECT     => State_Chm%Phot%IWVSELECT
    IRTWVSELECT   => State_Chm%Phot%IRTWVSELECT
    ACOEF_WV      => State_Chm%Phot%ACOEF_WV
    BCOEF_WV      => State_Chm%Phot%BCOEF_WV
    CCOEF_WV      => State_Chm%Phot%CCOEF_WV
    ACOEF_RTWV    => State_Chm%Phot%ACOEF_RTWV
    BCOEF_RTWV    => State_Chm%Phot%BCOEF_RTWV
    CCOEF_RTWV    => State_Chm%Phot%CCOEF_RTWV
    WVAA          => State_Chm%Phot%WVAA

    !cycle over standard wavelengths
    N0=1
    N1=NWVAA0
    NSTEP=1
    NWVREQUIRED=0
    DO W=1,Input_Opt%NWVSELECT
       MINWV     = -999
       MAXWV     =  999
       DO N=N0,N1,NSTEP ! 1 to 11
          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
             MINWV = WVDIF
             IWVSELECT(1,W)=N             
          ENDIF
          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
             MAXWV = WVDIF
             IWVSELECT(2,W)=N
          ENDIF
       ENDDO
       IF (IWVSELECT(2,W).EQ.IWVSELECT(1,W)) THEN
          !we have a match!
          MINWV=0
          MAXWV=0
          !add this wavelength to those for output
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ELSE
          !we are going to have to interpolate to the requested wavelength
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(2,W)
       ENDIF

       !Error check - ensure we have a match or requested wavelength
       !falls within two LUT tie points
       IF (MINWV.EQ.-999) THEN
          ! requested wavelength is shorter than min wv in LUT
          ! set to min
          write(6,*) 'ERROR requested wavelength is too short!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(1,1)
          IWVSELECT(1,W)=1
          IWVSELECT(2,W)=1 !300nm
          NWVREQUIRED=NWVREQUIRED-1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ENDIF
       IF (MAXWV.EQ.999) THEN
          ! requested wavelength is longer than min wv in LUT
          ! set to max
          write(6,*) 'ERROR requested wavelength is too long!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0,1)
          IWVSELECT(1,W)=NWVAA0
          IWVSELECT(2,W)=NWVAA0 !1020nm
          NWVREQUIRED=NWVREQUIRED-1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ENDIF

       !now calcualte the angstrom exponent coefs for interpolation -
       !this is done here to save time and repetition in aerosol_mod.F
       IF (IWVSELECT(1,W).NE.IWVSELECT(2,W)) THEN
          ACOEF_WV(W) = WVAA(IWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
          BCOEF_WV(W) =1.0d0/(LOG(WVAA(IWVSELECT(2,W),1)/ &
                                  WVAA(IWVSELECT(1,W),1)))
          !relative location of selected wavelength between tie points
          !for interpolating SSA and ASYM for output in aerosol_mod.F and
          !dust_mod.F
          CCOEF_WV(W) =(Input_Opt%WVSELECT(W)-WVAA(IWVSELECT(1,W),1))/ &
                      (WVAA(IWVSELECT(2,W),1)-WVAA(IWVSELECT(1,W),1))
       ENDIF
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) 'N WAVELENGTHS: ',Input_Opt%NWVSELECT
          write(6,*) 'WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
          write(6,*) 'WAVELENGTH REQUIRED:', NWVREQUIRED
          !write(6,*) IWVSELECT(1,W),WVAA(IWVSELECT(1,W),1)
          !write(6,*) IWVSELECT(2,W),WVAA(IWVSELECT(2,W),1)
          !write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
          write(6,*) '*********************************'
       ENDIF
    ENDDO !Input_Opt%NWVSELECT
#ifdef RRTMG
    !repeat for RRTMG wavelengths to get the closest wavelength
    !indices and the interpolation coefficients
    !Indices are relative to all wavelengths in the LUT i.e. the RRTMG
    !wavelengths start at NWVAA0+1
    N0=NWVAA0+1
    N1=NWVAA
    NSTEP=1
    NRTWVREQUIRED=0
    DO W=1,Input_Opt%NWVSELECT
       MINWV     = -999
       MAXWV     =  999
       DO N=N0,N1,NSTEP
          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
             MINWV = WVDIF
             IRTWVSELECT(1,W)=N
          ENDIF
          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
             MAXWV = WVDIF
             IRTWVSELECT(2,W)=N
          ENDIF
       ENDDO
       IF (IRTWVSELECT(2,W).EQ.IRTWVSELECT(1,W)) THEN
          !we have a match!
          MINWV=0
          MAXWV=0
          !add this wavelength to those for output
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ELSE
          !we are going to have to interpolate to the requested
          !wavelength
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(2,W)
       ENDIF

       !Error check - ensure we have a match or requested wavelength
       !falls within two LUT tie points
       IF (MINWV.EQ.-999) THEN
          ! requested wavelength is shorter than min wv in LUT
          ! set to min
          write(6,*) 'ERROR requested wavelength is too short!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA-1,1)
          IRTWVSELECT(1,W)=NWVAA-1
          IRTWVSELECT(2,W)=NWVAA-1
          NRTWVREQUIRED=NRTWVREQUIRED-1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ENDIF
       IF (MAXWV.EQ.999) THEN
          ! requested wavelength is longer than min wv in LUT
          ! set to max
          write(6,*) 'ERROR requested wavelength is too long!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0+1,1)
          IRTWVSELECT(1,W)=NWVAA0+1
          IRTWVSELECT(2,W)=NWVAA0+1
          NRTWVREQUIRED=NRTWVREQUIRED-1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ENDIF

       !now calcualte the angstrom exponent coefs for interpolation -
       !this is done here to save time and repetition in aerosol_mod.F
       IF (IRTWVSELECT(1,W).NE.IRTWVSELECT(2,W)) THEN
          ACOEF_RTWV(W) = WVAA(IRTWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
          BCOEF_RTWV(W) =1.0d0/(LOG(WVAA(IRTWVSELECT(2,W),1)/ &
                                    WVAA(IRTWVSELECT(1,W),1)))
          !relative location of selected wavelength between tie points
          !for interpolating SSA and ASYM for output in aerosol_mod.F and
          !dust_mod.F
          CCOEF_RTWV(W) =(Input_Opt%WVSELECT(W)-WVAA(IRTWVSELECT(1,W),1))/ &
                      (WVAA(IRTWVSELECT(2,W),1)-WVAA(IRTWVSELECT(1,W),1))
       ENDIF
       !convert wavelength index to that required by rrtmg_rad_transfer
       !i.e. without the standard and LW wavelengths
       IRTWVSELECT(1,W) = IRTWVSELECT(1,W) - NWVAA0 - NBNDLW
       IRTWVSELECT(2,W) = IRTWVSELECT(2,W) - NWVAA0 - NBNDLW
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) 'N RT WAVELENGTHS: ',Input_Opt%NWVSELECT
          write(6,*) 'RT WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
          write(6,*) 'RT WAVELENGTH REQUIRED:', NRTWVREQUIRED
          write(6,*) IRTWVSELECT(1,W),WVAA(IRTWVSELECT(1,W)+NWVAA0+NBNDLW,1)
          write(6,*) IRTWVSELECT(2,W),WVAA(IRTWVSELECT(2,W)+NWVAA0+NBNDLW,1)
          write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
          write(6,*) '*********************************'
       ENDIF
    ENDDO !Input_Opt%NWVSELECT
#endif

    ! Copy values back into State_Chm
    State_Chm%Phot%NWVREQUIRED   = NWVREQUIRED
    State_Chm%Phot%NRTWVREQUIRED = NRTWVREQUIRED

    ! Free pointers
    IWVREQUIRED   => NULL() 
    IRTWVREQUIRED => NULL()
    IWVSELECT     => NULL()
    IRTWVSELECT   => NULL()
    ACOEF_WV      => NULL()
    BCOEF_WV      => NULL()
    CCOEF_WV      => NULL()
    ACOEF_RTWV    => NULL()
    BCOEF_RTWV    => NULL()
    CCOEF_RTWV    => NULL()
    WVAA          => NULL()

  END SUBROUTINE CALC_AOD
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
    USE CMN_FJX_Mod, ONLY : AN_, NAA, TITLAA
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
    MIEDX(10+(NRHAER*NRH)+1) = 4  ! SSA/LBS/STS
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

END MODULE PHOTOLYSIS_MOD
