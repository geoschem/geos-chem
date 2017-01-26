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
! !PRIVATE MEMBER FUNCTIONS:
!
#if defined( NC_DIAG )
      PRIVATE :: Diagn_RnPbBe_Update
#endif

      ! NOTE: WE NEED TO KEEP THIS HERE FOR NOW SINCE IT IS USED
      ! BY CODE IN mercury_mod.F.  REMOVE IT LATER. (bmy, 7/7/14)
      ! Removed. (mdy, 12/4/14)
!      PUBLIC  :: SLQ
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Species ID flags
      INTEGER :: id_Rn, id_Pb, id_Be7

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
!  210Pb, and 7Be.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMRnPbBe( am_I_Root, Input_Opt, 
     &                       State_Met, State_Chm, RC )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRATMESO
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE DIAG_MOD,           ONLY : AD01
      USE DIAG_MOD,           ONLY : AD02
#endif
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                  :: I,        J,         L
      INTEGER                  :: N !,        NA,        nAdvect
      REAL(fp)                 :: ADD_Pb,   Be_LOST,   DTCHEM   
      REAL(fp)                 :: Pb_LOST,  Decay

      ! SAVEd scalars
      LOGICAL,       SAVE      :: FIRSTCHEM = .TRUE.
      REAL(fp),      SAVE      :: EXP_Rn,   EXP_Pb,    EXP_Be

      ! Arrays
      REAL(fp)                 :: Rn_LOST(IIPAR,JJPAR,LLPAR)

      ! To calculate decay rate

      ! For netCDF diagnostics via HEMCO
      CHARACTER(LEN=255)       :: MSG
      CHARACTER(LEN=255)       :: LOC = 'CHEMRNPBBE (RnPbBe_mod.F)'
      CHARACTER(LEN=31)        :: Spc_Name
      REAL(fp),      TARGET    :: PB_EMISS(IIPAR, JJPAR, LLPAR)
      REAL(fp),      TARGET    :: RN_DECAY(IIPAR, JJPAR, LLPAR)
      REAL(fp),      TARGET    :: PB_DECAY(IIPAR, JJPAR, LLPAR) 
      REAL(fp),      TARGET    :: BE_DECAY(IIPAR, JJPAR, LLPAR)  

      ! Pointers
      REAL(fp),      POINTER   :: Spc(:,:,:,:)
!
! !DEFINED PARAMETERS
!
      ! Ratio of molecular weights of 210Pb/222Rn
      REAL(fp),      PARAMETER :: Pb_Rn_RATIO = 210e+0_fp / 222e+0_fp

      ! Ln 2
      REAL(fp),      PARAMETER :: ln2   = 0.693147181E+00_fp

      ! Lifetime of Radon-222 in days
      REAL(fp),      PARAMETER :: RnTau = 3.83E+00_fp

      !=================================================================
      ! CHEMRnPbBe begins here!
      !=================================================================

      ! Assume success
      RC     =  GC_SUCCESS

      ! Chemistry timestep [s]
      DTCHEM =  GET_TS_CHEM() * 60e+0_fp

      ! Point to the species array
      Spc    => State_Chm%Species

      ! Pre-compute exponential terms only on first timestep
      IF ( FIRSTCHEM ) THEN 

         ! Fraction of (222Rn, 210Pb, 7Be) left after radioactive decay
         Decay     = ln2 / ( RnTau * 24.E+00_fp * 3600.E+00_fp ) 
         EXP_Rn    = EXP( -DTCHEM * Decay        )
         EXP_Pb    = EXP( -DTCHEM * 9.725E-10_fp ) 
         EXP_Be    = EXP( -DTCHEM * 1.506E-7_fp  )

         ! Species ID flags
         id_Rn     = Ind_('Rn'  )
         id_Pb     = Ind_('Pb'  )
         id_Be7    = Ind_('Be7' )

         ! Reset FIRSTCHEM flag
         FIRSTCHEM = .FALSE.

         ! testing only 
         IF(am_I_Root) THEN
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
      
      ! Species index for Rn222
      N = id_Rn

      ! Make sure Rn222 is a defined species
      IF ( N > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Rn_LOST = amount of 222Rn lost to decay [kg]
            Rn_LOST(I,J,L)     = Spc(I,J,L,N) * ( 1.0_fp - EXP_Rn )

            ! ND02 diag: 222Rn lost to decay [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
#if defined( BPCH_DIAG )
               AD02(I,J,L,1)   = AD02(I,J,L,1) 
     &                         + ( Rn_LOST(I,J,L) / DTCHEM )
#endif
#if defined( NC_DIAG )
               RN_DECAY(I,J,L) = RN_LOST(I,J,L) / DTCHEM
#endif
            ENDIF

            ! Subtract Rn_LOST from STT [kg]
            Spc(I,J,L,N)       = Spc(I,J,L,N) - Rn_LOST(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

#if defined( NC_DIAG )
         ! Get species name
         Spc_Name = State_Chm%SpcData(N)%Info%Name

         ! Make sure ND02 is turned on
         IF ( Input_Opt%ND02 > 0 ) THEN

            ! Update ND02 diagnostic: Rn lost to decay [kg/s]
            CALL Diagn_RnPbBe_Update( am_I_Root,
     &                                'ND02', 
     &                                Spc_Name, 
     &                                RN_DECAY, 
     &                                INput_Opt%DIAG_COLLECTION )
         ENDIF
#endif
         
      ENDIF

      !=================================================================
      ! Radioactive decay of Pb210
      !=================================================================
      
      ! Species index for Pb210
      N = id_Pb

      ! Make sure Pb210 is a defined species
      IF ( N > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, ADD_Pb, Pb_LOST )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
           
            ! ADD_Pb = Amount of 210Pb gained by decay from 222Rn [kg]
            ADD_Pb                     = Rn_LOST(I,J,L) * Pb_Rn_RATIO 

#if defined( GCAP )
            !%%%%% FOR GCAP SIMULATION: Divide the 210Pb source by 3.5
            !%%%%% to correct for the strat-trop exchange.  A similar
            !%%%%% line of code is done in the HEMCO extension where
            !%%%%% the Be7 source is computed.  This replicates the
            !%%%%% pre-HEMCO code. (bmy, 9/4/14)
            IF ( ITS_IN_THE_STRATMESO( I, J, L, State_Met ) ) THEN
               ADD_Pb          = ADD_Pb / 3.5e+0_fp
            ENDIF
#endif

            ! ND01 diag: 210Pb emission from 222Rn decay [kg/s]
            IF ( ND01 > 0 .and. L <= LD01 ) THEN
#if defined( BPCH_DIAG )
               AD01(I,J,L)     = AD01(I,J,L) + ( ADD_Pb / DTCHEM )
#endif
#if defined ( NC_DIAG )
               PB_EMISS(I,J,L) = ADD_Pb / DTCHEM
#endif
            ENDIF

            ! Add 210Pb gained by decay from 222Rn into STT [kg]
            Spc(I,J,L,N)       = Spc(I,J,L,N) + ADD_Pb          

            ! Amount of 210Pb lost to radioactive decay [kg]
            ! NOTE: we've already added in the 210Pb gained from 222Rn
            Pb_LOST            = Spc(I,J,L,N) * ( 1.0_fp - EXP_Pb )

            ! ND02 diag: 210Pb lost to decay [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
#if defined( BPCH_DIAG )
               AD02(I,J,L,2)   = AD02(I,J,L,2) + ( Pb_LOST / DTCHEM )
#endif
#if defined( NC_DIAG )
               PB_DECAY(I,J,L) = Pb_LOST / DTCHEM
#endif
            ENDIF

            ! Subtract 210Pb lost to decay from STT [kg]
            Spc(I,J,L,N)       = Spc(I,J,L,N) - Pb_LOST
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

#if defined( NC_DIAG )
         ! Get species name
         Spc_Name = State_Chm%SpcData(N)%Info%Name

         ! Make sure ND02 is turned on
         IF ( Input_Opt%ND02 > 0 ) THEN

            ! Update ND02 diagnostic: Pb lost to decay [kg/s]
            CALL Diagn_RnPbBe_Update( am_I_Root,
     &                                'ND02', 
     &                                Spc_Name, 
     &                                PB_DECAY, 
     &                                Input_Opt%DIAG_COLLECTION )
         ENDIF

         ! Make sure ND01 is turned on
         IF ( Input_Opt%ND01 > 0 ) THEN

            ! Update ND01 diagnostic: Pb emission [kg/s]
            CALL Diagn_RnPbBe_Update( am_I_Root,
     &                                'ND01', 
     &                                Spc_Name, 
     &                                PB_EMISS, 
     &                                Input_Opt%DIAG_COLLECTION )
         ENDIF

#endif

      ENDIF

      !=================================================================
      ! Radioactive decay of Be7
      !=================================================================

      ! Species index for Be7
      N = id_Be7

      ! Make sure Be7 is a defined species
      IF ( N > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, Be_LOST )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Amount of 7Be lost to decay [kg]
            Be_LOST            = Spc(I,J,L,N) * ( 1d0 - EXP_Be )

            ! ND02 diag: 7Be lost to decay [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
#if defined( BPCH_DIAG )
               AD02(I,J,L,3)   = AD02(I,J,L,3) + ( Be_LOST / DTCHEM )
#endif
#if defined( NC_DIAG )
               BE_DECAY(I,J,L) = Be_LOST / DTCHEM
#endif
            ENDIF

            ! Subtract amount of 7Be lost to decay from STT [kg]
            Spc(I,J,L,N)       = Spc(I,J,L,N) - Be_LOST
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

#if defined( NC_DIAG )
         ! Get species name
         Spc_Name = State_Chm%SpcData(N)%Info%Name

         ! Make sure ND02 is turned on
         IF ( Input_Opt%ND02 > 0 ) THEN

            ! Update ND02 diagnostic: Be lost to decay [kg/s]
            CALL Diagn_RnPbBe_Update( am_I_Root,
     &                                'ND02', 
     &                                Spc_Name, 
     &                                BE_DECAY, 
     &                                Input_Opt%DIAG_COLLECTION )
         ENDIF
#endif

      ENDIF

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE CHEMRnPbBe
!EOC

#if defined( NC_DIAG )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: diagn_rnpbbe_update
!
! !DESCRIPTION: Subroutine DIAGN\_RNPBBE\_UPDATE updates the Rn/Pb/Be7
!  diagnostics ND01 (emissions) and ND02 (loss due to decay) that are 
!  written to netCDF. ND01 Rn and Be7 emissions diagnostics are included
!  in HEMCO emissions and are therefore not updated in this subroutine.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Diagn_RnPbBe_Update( am_I_Root, NDxx, SpeciesName, 
     &                                DiagArray, ColNo )
!
! !USES:
!
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
      USE HCO_Error_Mod
      USE Error_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL,            INTENT(IN)  :: am_I_Root    ! Are we on the root CPU?
      CHARACTER(LEN=4),   INTENT(IN)  :: NDxx         ! Diagnostic id
      CHARACTER(LEN=*),   INTENT(IN)  :: SpeciesName  ! Species name
      REAL(fp), TARGET,   INTENT(IN)  :: DiagArray(:,:,:)
      INTEGER,            INTENT(IN)  :: ColNo        ! Diag collection #
!
! !REVISION HISTORY:
!  30 Jan 2015 - E. Lundgren   - Initial version
! 
! !REMARKS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)  :: DiagnName
      CHARACTER(LEN=255)  :: MSG
      CHARACTER(LEN=255)  :: LOC = 'DIAGN_RNPBBE_UPDATE (RnPbBe_mod.F)'
      INTEGER             :: HCRC
      REAL(fp), POINTER   :: Ptr3D(:,:,:) => NULL()

      !=================================================================
      ! DIAGN_RNPBBE_UPDATE begins here!
      !=================================================================
      
      ! Define the diagnostic name
      SELECT CASE ( NDxx )
         CASE ( 'ND01' )
            DiagnName = 'EMISS_' // TRIM( SpeciesName )
         CASE ( 'ND02' )
            DiagnName = 'DECAY_' // TRIM( SpeciesName )
         CASE DEFAULT
            MSG = 'Undefined Rn/Pb/Be7 diagnostic: ' // NDxx
            CALL ERROR_STOP( TRIM( MSG ), LOC )
      END SELECT

      ! Point to the array
      Ptr3D => DiagArray

      ! Update container for diagnostic
      CALL Diagn_Update( am_I_Root, HcoState,
     &                   cName   = TRIM( DiagnName ),
     &                   Array3D = Ptr3D,
     &                   COL     = ColNo,
     &                   RC      = HCRC )
 
      ! Free the pointer before error handling
      Ptr3D => NULL()
 
      ! Stop with error if the diagnostic update was unsuccessful
      IF ( HCRC /= HCO_SUCCESS ) THEN
         MSG = 'Cannot update ' // NDxx // ' diagnostic: ' 
     &         // TRIM( DiagnName )
         CALL ERROR_STOP( MSG, LOC)
      ENDIF

      END SUBROUTINE Diagn_RnPbBe_Update 
!EOC
#endif
      END MODULE RnPbBe_MOD


