!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: RnPbBe_mod.F
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
!  14 Jun 2001 - H. Liu      - Initial version  
!  (1 ) Added existing routines to this module (bmy, 6/14/01)
!  (2 ) Updated comments (bmy, 9/4/01)
!  (3 ) Eliminate AVGF; redimensioned XTRA2 (bmy, 9/25/01)
!  (4 ) Replace references to PW(I,J) with P(I,J) (bmy, 10/3/01)
!  (5 ) Remove obsolete code from 9/01 and 10/01 (bmy, 10/23/01)
!  (6 ) Removed duplicate variable declarations (bmy, 11/15/01)
!  (7 ) Now read files from DATA_DIR/RnPbBe_200203/ directory.  
!        Also updated comments. (bmy, 3/29/02)
!  (8 ) Incorporated latest changes from Hongyu Liu.  Also split off the
!        code to read in the 7Be emissions into a separate routine. 
!        Add parallel DO-loops in several places.  Cleaned up DRYFLXRnPbBe,
!        and now make sure ND44 accurately represents the drydep fluxes
!        of 210Pb and 7Be. (hyl, bmy, 8/7/02)
!  (9 ) Now reference AD from "dao_mod.f".  Now references "error_mod.f".
!        Moved routine DRYFLXRnPbBe into "drydep_mod.f".  (bmy, 1/27/03)
!  (10) Now references the new "time_mod.f" (bmy, 2/11/03)
!  (11) Bug fix in EMISSRnPbBe -- take abs( lat) for 7Be emiss. (bmy, 6/10/03)
!  (12) Bug fix in EMISSRnPbBe -- shut off 222Rn emissions in polar regions
!        (swu, bmy, 10/28/03)
!  (13) Now references "directory_mod.f", "logical_mod.f", and "tracer_mod.f"
!        (bmy, 7/20/04)
!  (14) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 5/24/05)
!  (15) Now references "tropopause_mod.f"
!  (16) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  19 Nov 2010 - R. Yantosca - Added ProTeX headers
!  08 Nov 2011 - R. Yantosca - Prevent out-of-bounds errors in diagnostics
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use routines from the new grid_mod.F90
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  07 Jul 2014 - R. Yantosca - Removed routines now orphaned by HEMCO
!  22 Aug 2014 - R. Yantosca - Removed LATSOU, PRESOU, BESOU arrays, these
!                              are now defined in the HEMCO code.
!  22 Aug 2014 - R. Yantosca - Remove XNUMOL_* parameters; these are obsolete
!  04 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  30 Jan 2015 - E. Lundgren - Add new diagnostics stuctures for netCDF output.
!  02 May 2016 - R. Yantosca - Now define IDTRn, IDTPb, IDTBe locally
!  06 Dec 2018 - M. Sulprizio- Add Be10 to mechanism following GMI tracer suite
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
      SUBROUTINE CHEMRnPbBe( am_I_Root, Input_Opt,  State_Met,
     &                       State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,       ONLY : AD01
      USE DIAG_MOD,       ONLY : AD02
#endif
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Chm_Mod,  ONLY : Ind_
      USE State_Diag_Mod, ONLY : DgnState
      USE State_Met_Mod,  ONLY : MetState
      USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  (1 ) Now use F90 syntax (bmy, hyl, 3/22/99)
!  (2 ) Add FIRSTCHEM as an argument.  Only compute the exponential terms
!        when FIRSTCHEM = .TRUE., and save the values for later use
!        (bmy, 3/24/99)
!  (3 ) Cosmetic changes (bmy, 10/13/99)
!  (4 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (5 ) Cosmetic changes (bmy, 7/12/00)
!  (6 ) Added to module "RnPbBe_mod.f".  Also updated comments 
!        and made cosmetic changes. (bmy, 6/14/01)
!  (7 ) Add diagnostics for Rn/Be emissions.  Also cleaned up some old code
!        and added parallel DO-loops.  Updated comments. (hyl, 8/6/02)
!  (8 ) Now make FIRSTCHEM a local SAVEd variable.  (bmy, 1/27/03)
!  (9 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 2/11/03)
!  (10) Now references STT and N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (11) Remove reference to CMN; it's obsolete.  Now use inquiry functions 
!        from "tropopause_mod.f" to diagnose strat boxes. (bmy, 8/15/05)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  08 Nov 2011 - R. Yantosca - Prevent out-of-bounds errors in diagnostics
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  22 Aug 2014 - R. Yantosca - Copy emissions above the PBL to tracer array
!  22 Aug 2014 - R. Yantosca - Cosmetic changes, for clarity
!  04 Sep 2014 - R. Yantosca - Add minor changes for efficiency
!  04 Nov 2014 - M. Yannetti - Changed REAL*8 to REAL(fp) 
!  12 Dec 2014 - M. Yannetti - Changed Hemco REAL*8 to REAL(hp)
!   2 May 2016 - R. Yantosca - Now define IDTRn, IDTPb, IDTBe locally
!  22 Jun 2016 - R. Yantosca - Now use Ind_() to define species ID's
!  22 Jun 2016 - R. Yantosca - Rename species ID's to id_Rn, id_Pb, id_Be7
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  10 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
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
      REAL(fp)            :: Rn_LOST(IIPAR,JJPAR,LLPAR)

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
      ThisLoc =  ' -> at ChemRnPbBe (in module GeosCore/RnPbBe_mod.F)'

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
         IF ( am_I_Root ) THEN
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
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Rn_LOST = amount of 222Rn lost to decay [kg]
            Rn_LOST(I,J,L) = Spc(I,J,L,id_Rn222) * ( 1.0_fp - EXP_Rn )

#if defined( BPCH_DIAG )
            !-----------------------------------------------------------
            ! ND02 (bpch) diagnostic:
            !
            ! Rn222 lost to radioactive decay
            !-----------------------------------------------------------

            ! Units: [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
               AD02(I,J,L,1)   = AD02(I,J,L,1) 
     &                         + ( Rn_LOST(I,J,L) / DTCHEM )
            ENDIF
#endif

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
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, ADD_Pb, Pb_LOST, PbStrat_LOST )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
           
            ! ADD_Pb = Amount of 210Pb gained by decay from 222Rn [kg]
            ADD_Pb = Rn_LOST(I,J,L) * Pb_Rn_RATIO 

#if defined( BPCH_DIAG )
            !-----------------------------------------------------------
            ! ND01 (bpch) diagnostic
            !
            ! Pb210 emission from Rn222 decay
            !-----------------------------------------------------------

            ! Units: [kg/s]
            IF ( ND01 > 0 .and. L <= LD01 ) THEN
               AD01(I,J,L,1) = AD01(I,J,L,1) + ( ADD_Pb / DTCHEM )

               IF ( State_Met%InStratosphere(I,J,L) ) THEN
                  AD01(I,J,L,2) = AD01(I,J,L,2) + ( ADD_Pb / DTCHEM )
               ENDIF
            ENDIF
#endif

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
            IF ( State_Met%InStratosphere(I,J,L) .and.
     &           id_Pb210Strat > 0 ) THEN
               Spc(I,J,L,id_Pb210Strat) = Spc(I,J,L,id_Pb210Strat) +
     &                                    ADD_Pb
            ENDIF

            ! Amount of 210Pb lost to radioactive decay [kg]
            ! NOTE: we've already added in the 210Pb gained from 222Rn
            Pb_LOST = Spc(I,J,L,id_Pb210) * ( 1.0_fp - EXP_Pb )
            IF ( id_Pb210Strat > 0 ) THEN
               PbStrat_LOST = Spc(I,J,L,id_Pb210Strat) *
     &                        ( 1.0_fp - EXP_Pb)
            ENDIF

#if defined( BPCH_DIAG )
            !-----------------------------------------------------------
            ! ND02 (bpch) diagnostic:
            !
            ! Pb210 lost to radioactive decay
            !-----------------------------------------------------------

            ! Units: [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
               AD02(I,J,L,2) = AD02(I,J,L,2) + ( Pb_LOST      / DTCHEM )
               AD02(I,J,L,3) = AD02(I,J,L,3) + ( PbStrat_LOST / DTCHEM )
            ENDIF
#endif

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
               Spc(I,J,L,id_Pb210Strat) = Spc(I,J,L,id_Pb210Strat) -
     &                                    PbStrat_LOST
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
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, Be7_LOST, Be7Strat_LOST )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Amount of 7Be lost to decay [kg]
            Be7_LOST = Spc(I,J,L,id_Be7) * ( 1d0 - EXP_Be7 )
            IF ( id_Be7Strat > 0 ) THEN
               Be7Strat_LOST = Spc(I,J,L,id_Be7Strat) *
     &                         ( 1d0 - EXP_Be7 )
            ENDIF

#if defined( BPCH_DIAG )
            !-----------------------------------------------------------
            ! ND02 (bpch) diagnostic:
            !
            ! Be7 lost to radioactive decay
            !-----------------------------------------------------------

            ! Units: [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
               AD02(I,J,L,4) = AD02(I,J,L,4) + ( Be7_LOST      / DTCHEM)
               AD02(I,J,L,5) = AD02(I,J,L,5) + ( Be7Strat_LOST / DTCHEM)
            ENDIF
#endif

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
               Spc(I,J,L,id_Be7Strat) = Spc(I,J,L,id_Be7Strat) -
     &                                  Be7Strat_LOST
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
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, Be10_LOST, Be10Strat_LOST )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Amount of 10Be lost to decay [kg]
            Be10_LOST = Spc(I,J,L,id_Be10) * ( 1d0 - EXP_Be10 )
            IF ( id_Be10Strat > 0 ) THEN
               Be10Strat_LOST = Spc(I,J,L,id_Be10Strat) *
     &                          ( 1d0 - EXP_Be10 )
            ENDIF

#if defined( BPCH_DIAG )
            !-----------------------------------------------------------
            ! ND02 (bpch) diagnostic:
            !
            ! Be10 lost to radioactive decay
            !-----------------------------------------------------------

            ! Units: [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
               AD02(I,J,L,6) = AD02(I,J,L,6) + ( Be10_LOST      /DTCHEM)
               AD02(I,J,L,7) = AD02(I,J,L,7) + ( Be10Strat_LOST /DTCHEM)
            ENDIF
#endif

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
               Spc(I,J,L,id_Be10Strat) = Spc(I,J,L,id_Be10Strat) -
     &                                   Be10Strat_LOST
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


