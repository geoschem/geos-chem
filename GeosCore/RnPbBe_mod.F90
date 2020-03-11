!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: RnPbBe_mod.F90
!
! !DESCRIPTION: Module RnPbBe\_MOD contains variables and routines used
!  for the 222Rn-210Pb-7Be simulation. (hyl, swu, bmy, 6/14/01, 8/4/06)
!\\
!\\
! !INTERFACE:
!
MODULE RnPbBe_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CHEMRnPbBe
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Liu,H., D.Jacob, I.Bey, and R.M.Yantosca, Constraints from 210Pb
!        and 7Be on wet deposition and transport in a global three-dimensional
!        chemical tracer model driven by assimilated meteorological fields,
!        JGR, 106, D11, 12,109-12,128, 2001.
!  (2 ) Jacob et al.,Evaluation and intercomparison of global atmospheric
!        transport models using Rn-222 and other short-lived tracers,
!        JGR, 1997 (102):5953-5970
!  (3 ) Dorothy Koch, JGR 101, D13, 18651, 1996.
!  (4 ) Lal, D., and B. Peters, Cosmic ray produced radioactivity on the
!        Earth. Handbuch der Physik, 46/2, 551-612, edited by K. Sitte,
!        Springer-Verlag, New York, 1967.
!  (5)  GMI tracer suite, https://gmi.gsfc.nasa.gov/uploads/files/gmi_tracersuite.pdf
!  (6 ) Koch and Rind, Beryllium 10/beryllium 7 as a tracer of stratospheric
!        transport, JGR, 103, D4, 3907-3917, 1998.
!
! !REVISION HISTORY:
!  14 Jun 2001 - H. Liu - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
  ! Species ID flags
  INTEGER  :: id_Rn222
  INTEGER  :: id_Pb210, id_Pb210Strat
  INTEGER  :: id_Be7,   id_Be7Strat
  INTEGER  :: id_Be10,  id_Be10Strat

  ! Exponential terms
  REAL(fp) :: EXP_Rn, EXP_Pb, EXP_Be7, EXP_Be10

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemRnPbBe
!
! !DESCRIPTION: Subroutine CHEMRnPbBe performs loss chemistry on 222Rn,
!  210Pb, 7Be, and 10Be.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEMRnPbBe( Input_Opt,  State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  31 Oct 1999 - H. Liu - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: I, J, L
    REAL(fp)            :: ADD_Pb
    REAL(fp)            :: Decay
    REAL(fp)            :: DTCHEM
    REAL(fp)            :: Pb_LOST,   PbStrat_LOST
    REAL(fp)            :: Be7_LOST,  Be7Strat_LOST
    REAL(fp)            :: Be10_LOST, Be10Strat_LOST

    ! SAVEd scalars
    LOGICAL, SAVE       :: FIRSTCHEM = .TRUE.

    ! Strings
    CHARACTER(LEN=255)  :: ErrMsg
    CHARACTER(LEN=255)  :: ThisLoc

    ! Arrays
    REAL(fp)            :: Rn_LOST(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Pointers
    REAL(fp), POINTER   :: Spc(:,:,:,:)
!
! !DEFINED PARAMETERS
!
    ! Ratio of molecular weights of 210Pb/222Rn
    REAL(fp), PARAMETER :: Pb_Rn_RATIO = 210e+0_fp / 222e+0_fp

    ! Ln 2
    REAL(fp), PARAMETER :: ln2         = 0.693147181E+00_fp

    ! Half-life in days
    REAL(fp), PARAMETER :: RnTau       = 3.83E+00_fp ! Liu et al. (2001)
    REAL(fp), PARAMETER :: Pb210Tau    = 8.25E+03_fp ! Liu et al. (2001)
    REAL(fp), PARAMETER :: Be7Tau      = 53.3E+00_fp ! Liu et al. (2001)
    REAL(fp), PARAMETER :: Be10Tau     = 5.84E+08_fp ! GMI "tracer" mechanism

    !=================================================================
    ! CHEMRnPbBe begins here!
    !=================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    ErrMsg  =  ''
    ThisLoc =  ' -> at ChemRnPbBe (in module GeosCore/RnPbBe_mod.F90)'

    ! Chemistry timestep [s]
    DTCHEM  =  GET_TS_CHEM()

    ! Point to the species array
    Spc     => State_Chm%Species

    !-----------------------------------------------------------------
    ! Pre-compute exponential terms and do other first-time setup
    !-----------------------------------------------------------------
    IF ( FIRSTCHEM ) THEN

       ! Fraction of Rn222 left after radioactive decay
       Decay     = ln2 / ( RnTau * 24.E+00_fp * 3600.E+00_fp )
       EXP_Rn    = EXP( -DTCHEM * Decay        )

       ! Fraction of Pb210 left after radioactive decay
       !Decay    = ln2 / ( PbTau * 24.E+00_fp * 3600.E+00_fp )
       EXP_Pb    = EXP( -DTCHEM * 9.725E-10_fp )

       ! Fraction of Be7 left after radioactive decay
       !Decay    = ln2 / ( Be7Tau * 24.E+00_fp * 3600.E+00_fp )
       EXP_Be7   = EXP( -DTCHEM * 1.506E-7_fp  )

       ! Fraction of Be10 left after radioactive decay
       Decay     = ln2 / ( Be10Tau * 24.E+00_fp * 3600.E+00_fp )
       EXP_Be10  = EXP( -DTCHEM * Decay        )

       ! Species ID flags
       id_Rn222      = Ind_('Rn222'      )
       id_Pb210      = Ind_('Pb210'      )
       id_Pb210Strat = Ind_('Pb210Strat' )
       id_Be7        = Ind_('Be7'        )
       id_Be7Strat   = Ind_('Be7Strat'   )
       id_Be10       = Ind_('Be10'       )
       id_Be10Strat  = Ind_('Be10Strat'  )

       ! Reset FIRSTCHEM flag
       FIRSTCHEM = .FALSE.

       ! testing only
       IF ( Input_Opt%amIRoot ) THEN
          write(*,*) ''
          write(*,*) '### GEOS-Chem Radon simulation ###'
          write(*,*) '    Timestep (secs)   : ', DTCHEM
          write(*,*) '    Rn lifetime (days): ', RnTau
          write(*,*) '    Rn decadence      : ', EXP_Rn
          write(*,*) ''
       ENDIF
    ENDIF

    !=================================================================
    ! Radioactive decay of Rn222
    !=================================================================

    ! Make sure Rn222 is a defined species
    IF ( id_Rn222 > 0 ) THEN
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Rn_LOST = amount of 222Rn lost to decay [kg]
          Rn_LOST(I,J,L) = Spc(I,J,L,id_Rn222) * ( 1.0_fp - EXP_Rn )

          !-----------------------------------------------------------
          ! HISTORY (aka netCDF diagnostics)
          !
          ! Rn222 lost to radioactive decay
          !-----------------------------------------------------------

          ! Units: [kg/s], but consider eventually changing to [kg/m2/s]
          IF ( State_Diag%Archive_RadDecay ) THEN
             State_Diag%RadDecay(I,J,L,1) = Rn_LOST(I,J,L) / DTCHEM
          ENDIF

          ! Subtract Rn_LOST from Spc [kg]
          Spc(I,J,L,id_Rn222) = Spc(I,J,L,id_Rn222) - Rn_LOST(I,J,L)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !=================================================================
    ! Radioactive decay of Pb210
    !=================================================================

    ! Make sure Pb210 is a defined species
    IF ( id_Pb210 > 0 .or. id_Pb210Strat > 0 ) THEN
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, ADD_Pb, Pb_LOST, PbStrat_LOST )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! ADD_Pb = Amount of 210Pb gained by decay from 222Rn [kg]
          ADD_Pb = Rn_LOST(I,J,L) * Pb_Rn_RATIO

          !-----------------------------------------------------------
          ! HISTORY (aka netCDF diagnostics)
          !
          ! Pb210 emission from Rn222 decay
          !-----------------------------------------------------------

          ! Units: [kg/s], but consider eventually changing to [kg/m2/s]
          IF ( State_Diag%Archive_PbFromRnDecay ) THEN
             State_Diag%PbFromRnDecay(I,J,L) = ( ADD_Pb / DTCHEM )
          ENDIF

          ! Add 210Pb gained by decay from 222Rn into Spc [kg]
          Spc(I,J,L,id_Pb210) = Spc(I,J,L,id_Pb210) + ADD_Pb

          ! Update stratospheric 210Pb [kg]
          IF ( State_Met%InStratosphere(I,J,L) .and. id_Pb210Strat > 0 ) THEN
             Spc(I,J,L,id_Pb210Strat) = Spc(I,J,L,id_Pb210Strat) + ADD_Pb
          ENDIF

          ! Amount of 210Pb lost to radioactive decay [kg]
          ! NOTE: we've already added in the 210Pb gained from 222Rn
          Pb_LOST = Spc(I,J,L,id_Pb210) * ( 1.0_fp - EXP_Pb )
          IF ( id_Pb210Strat > 0 ) THEN
             PbStrat_LOST = Spc(I,J,L,id_Pb210Strat) * ( 1.0_fp - EXP_Pb)
          ENDIF

          !-----------------------------------------------------------
          ! HISTORY (aka netCDF diagnostics)
          !
          ! Pb210 lost to radioactive decay
          !-----------------------------------------------------------

          ! Units: [kg/s], but consider eventually changing to [kg/m2/s]
          IF ( State_Diag%Archive_RadDecay ) THEN
             State_Diag%RadDecay(I,J,L,2) = ( Pb_LOST / DTCHEM )
             State_Diag%RadDecay(I,J,L,3) = ( PbStrat_LOST /DTCHEM )
          ENDIF

          ! Subtract 210Pb lost to decay from Spc [kg]
          Spc(I,J,L,id_Pb210) = Spc(I,J,L,id_Pb210) - Pb_LOST

          ! Update stratospheric 210Pb [kg]
          IF ( id_Pb210Strat > 0 ) THEN
             Spc(I,J,L,id_Pb210Strat) = Spc(I,J,L,id_Pb210Strat) - PbStrat_LOST
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !=================================================================
    ! Radioactive decay of Be7
    !=================================================================

    ! Make sure Be7 is a defined species
    IF ( id_Be7 > 0 .or. id_Be7Strat > 0 ) THEN
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, Be7_LOST, Be7Strat_LOST )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Amount of 7Be lost to decay [kg]
          Be7_LOST = Spc(I,J,L,id_Be7) * ( 1d0 - EXP_Be7 )
          IF ( id_Be7Strat > 0 ) THEN
             Be7Strat_LOST = Spc(I,J,L,id_Be7Strat) * ( 1d0 - EXP_Be7 )
          ENDIF

          !-----------------------------------------------------------
          ! HISTORY (aka netCDF diagnostics)
          !
          ! Be7 lost to radioactive decay
          !-----------------------------------------------------------

          ! Units: [kg/s], but consider eventually changing to [kg/m2/s]
          IF ( State_Diag%Archive_RadDecay ) THEN
             State_Diag%RadDecay(I,J,L,4) = ( Be7_LOST      / DTCHEM )
             State_Diag%RadDecay(I,J,L,5) = ( Be7Strat_LOST / DTCHEM )
          ENDIF

          ! Subtract amount of 7Be lost to decay from Spc [kg]
          Spc(I,J,L,id_Be7) = Spc(I,J,L,id_Be7) - Be7_LOST

          ! Update stratospheric 7Be [kg]
          IF ( id_Be7Strat > 0 ) THEN
             Spc(I,J,L,id_Be7Strat) = Spc(I,J,L,id_Be7Strat) - Be7Strat_LOST
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !=================================================================
    ! Radioactive decay of Be10
    !=================================================================

    ! Make sure Be10 is a defined species
    IF ( id_Be10 > 0 .or. id_Be10Strat > 0 ) THEN
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, Be10_LOST, Be10Strat_LOST )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Amount of 10Be lost to decay [kg]
          Be10_LOST = Spc(I,J,L,id_Be10) * ( 1d0 - EXP_Be10 )
          IF ( id_Be10Strat > 0 ) THEN
             Be10Strat_LOST = Spc(I,J,L,id_Be10Strat) * ( 1d0 - EXP_Be10 )
          ENDIF

          !-----------------------------------------------------------
          ! HISTORY (aka netCDF diagnostics)
          !
          ! Be10 lost to radioactive decay
          !-----------------------------------------------------------

          ! Units: [kg/s], but consider eventually changing to [kg/m2/s]
          IF ( State_Diag%Archive_RadDecay ) THEN
             State_Diag%RadDecay(I,J,L,6) = ( Be10_LOST      / DTCHEM)
             State_Diag%RadDecay(I,J,L,7) = ( Be10Strat_LOST / DTCHEM)
          ENDIF

          ! Subtract amount of 10Be lost to decay from Spc [kg]
          Spc(I,J,L,id_Be10) = Spc(I,J,L,id_Be10) - Be10_LOST

          ! Update stratospheric 10Be [kg]
          IF ( id_Be10Strat > 0 ) THEN
             Spc(I,J,L,id_Be10Strat) = Spc(I,J,L,id_Be10Strat) - Be10Strat_LOST
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CHEMRnPbBe
!EOC
END MODULE RnPbBe_MOD
