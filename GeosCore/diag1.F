!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag1.F
!
! !DESCRIPTION: Subroutine DIAG1 accumulates diagnostic quantities on every 
!  dynamic timestep.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG1( am_I_Root,  Input_Opt, State_Chm,
     &                  State_Grid, State_Met, RC )
!
! !USES:
!
      ! References to F90 modules
      USE CMN_SIZE_MOD
      USE DAO_MOD,            ONLY : IS_ICE, IS_WATER, IS_LAND
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD33, AD35, AD54 
      USE DIAG_MOD,           ONLY : AD47
      USE DIAG_MOD,           ONLY : AD30, AD31, AD45, AD57 
      USE DIAG_MOD,           ONLY : AD67, AD68, AD69, LTOTH
      USE DIAG_MOD,           ONLY : AD71
      USE DIAG_MOD,           ONLY : AD71_DAY,   AD71_HR
      USE DIAG_MOD,           ONLY : AD71_LDAY,  AD71_LHR
      USE DIAG_MOD,           ONLY : AD71_COUNT, AD71_HRCT
      USE DIAG03_MOD,         ONLY : AD03_RGM, AD03_PBM, ND03   
#endif
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
      USE HCO_ERROR_MOD
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants
      USE PRECISION_MOD
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Grid_Mod,     ONLY : GrdState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : ITS_TIME_FOR_CHEM
      USE TIME_MOD,           ONLY : GET_DAY
      USE TIME_MOD,           ONLY : GET_HOUR
      USE HCO_INTERFACE_MOD,  ONLY : HcoState

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
! !REMARKS:
!  For a complete list of GEOS-Chem diagnostics, please see this web page:
!     http://acmg.seas.harvard.edu/geos/doc/man/appendix_5.html
!
! !REVISION HISTORY:
!  (1 ) This subroutine was reconstructed from gmg's version of (10/10/97)
!  (2 ) GISS-specific code has been eliminated (bmy, 3/15/99)
!  (3 ) UWND, VWND, WW no longer needs to be passed (bmy, 4/7/99)
!  (4 ) Use F90 syntax for declarations, etc (bmy, 4/7/99)
!  (5 ) Remove counter KWACC...this is now redundant (bmy, 11/5/99)
!  (6 ) ND31, ND33, ND35, ND67, and ND69 now use dynamically 
!        allocatable arrays declared in "diag_mod.f". (bmy, 3/9/00)
!  (7 ) LTOTH is now an allocatable array in "diag_mod.f". (bmy, 3/17/00)
!  (8 ) Add parallel loops over tracer where expedient (bmy, 5/4/00)
!  (9 ) Updated comments and diagnostics list.  Also add more parallel
!        loops for ND31 and ND68.  (bmy, 6/21/00)
!  (10) Use NTRACE to dimension STT_VV instead of NNPAR (bmy, 10/17/00)
!  (11) Removed obsolete code from 10/17/00 (bmy, 12/21/00)
!  (12) Updated diagnostic list & comments, cosmetic changes (bmy, 6/19/01)
!  (13) Updated diagnostic list & comments (bmy, 9/4/01)
!  (14) Now reference AVGW from "dao_mod.f", and make sure it is allocated
!        before we reference it in the ND68 diagnostic.  Also reference PBL, 
!        PS, AIRDEN from "dao_mod.f". (bmy, 9/25/01)
!  (15) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (16) Renamed ND33 to "ATMOSPHERIC COLUMN SUM OF TRACER", since this is
!        a sum over all levels and not just in the troposphere.  Also
!        removed more obsolete code from 9/01.  Now use P(I,J)+PTOP instead
!        of PS, since that is the way to ensure that we use will be used
!        consistently.  Remove reference to PS from "dao_mod.f"(bmy, 4/11/02)
!  (17) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE.  Also removed
!        obsolete, commented-out code.  Also now replaced reference to
!        P(IREF,JREF) with P(I,J). (bmy, 6/25/02)
!  (18) Replaced references to P(I,J) with call to GET_PEDGE(I,J,1) from
!        "pressure_mod.f"  Eliminated obsolete commented-out code from
!        6/02. (dsa, bdf, bmy, 8/20/02)
!  (19) Now reference AD, and BXHEIGHT from "dao_mod.f".  Removed obsolete 
!        code.  Now refEerence IDTOX from "tracerid_mod.f". (bmy, 11/6/02)
!  (20) Now replace DXYP(J) with routine GET_AREA_M2 from "grid_mod.f"
!        (bmy, 2/4/03)
!  (21) Now compute PBL top for ND67 for GEOS-4/fvDAS.  Also now include
!        SCALE_HEIGHT from header file "CMN_GCTM". (bmy, 6/23/03)
!  (22) Now references N_TRACERS, STT, and ITS_A_FULLCHEM_SIM from
!        "tracer_mod.f" (bmy, 7/20/04)
!  (23) Fixed ND67 PS-PBL for GCAP and GEOS-5 met fields (swu, bmy, 6/9/05)
!  (24) Now archive ND30 diagnostic for land/water/ice flags (bmy, 8/18/05)
!  (25) Now reference XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (26) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (27) Added count for time in the troposphere - array AD54 (phs, 9/22/06)
!  (28) Now only archive O3 in ND45 and ND47 at chem timsteps (phs, 1/24/07)
!  (29) Bug fix: Update ND30 for both GEOS-3 and otherwise.  Also now save
!        3-D pressure edges in ND31 instead of PS-PTOP.  Revert to the !
!        pre-near-land ND30 diagnostic algorithm. (bmy, 1/28/04)
!  (30) Use LTO3 for O3 in ND45. (ccc, 7/20/09)
!  (31) Add potential temperature diagnostic in ND57 (fp, 2/3/10)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  15 Feb 2011 - R. Yantosca - Added modifications for APM from G. Luo
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 as part of removal of
!                              NOx-Ox partitioning
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  29 Aug 2013 - R. Yantosca - In ND57, we need to make ND57 !$OMP PRIVATE
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  07 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  20 Jan 2015 - R. Yantosca - Added new netCDF diagnostics
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE and
!                              GET_PCENTER with State_Met%PMID.
!  24 Mar 2015 - E. Lundgren - Remove dependency on tracer_mod
!  25 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg
!  16 Apr 2015 - E. Lundgren - Add new State_Met variables to ND68
!  19 Oct 2015 - C. Keller   - Now use Input_Opt%ND68 instead of ND68;
!                              Rename AIRDEN diagnostics to AIRDENSITY to
!                              avoid name conflict with GEOS-5 model.
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  02 May 2016 - R. Yantosca - Now declare IDTPASV locally
!  31 May 2016 - E. Lundgren - Replace input_opt%TRACER_MW_G with species
!                              database field emMW_g (emitted species g/mol)
!  31 May 2016 - E. Lundgren - Remove usage of TCVV; replace with AIRMW/emMW_g
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with ThisSpc%Name from species database
!  22 Jun 2016 - R. Yantosca - Now use Ind_() to define id_PASV species ID
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  16 Sep 2016 - E. Lundgren - Remove passive species that are per total air
!                              since moisture fix corrects v/v dry
!  20 Sep 2016 - E. Lundgren - Simplify met fields included in ND68 (8 total)
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER       :: I,          J,          K
      INTEGER       :: L,          N,          NA
      REAL(fp)      :: EmMW_g,     P0,         Spc_VV

      ! SAVEd scalars
      LOGICAL, SAVE :: FIRST   = .TRUE.
      LOGICAL, SAVE :: Do_ND03 = .FALSE.
      INTEGER, SAVE :: id_Hg2  = -1
      INTEGER, SAVE :: id_HgP  = -1
      REAL(fp),SAVE :: EmMW_g_Hg2 = -1.0_fp
      REAL(fp),SAVE :: EmMW_g_HgP = -1.0_fp

      !=================================================================
      ! DIAG1 begins here!
      !=================================================================

      ! Initialize
      RC      =  GC_SUCCESS

#if defined( BPCH_DIAG )
      ! Are certain diagnostics turned on?
      IF ( .not. Input_Opt%ND03 > 0 ) RETURN

      !=================================================================
      ! First-time setup
      !=================================================================
      IF ( FIRST ) THEN 

         ! Define species ID's  on the first call
         id_Hg2  = Ind_( 'Hg2' )
         id_HgP  = Ind_( 'HgP' )

         ! Pre-save the molecular weight of Hg2 [g]
         IF ( id_Hg2 > 0 ) THEN
            EmMW_g_Hg2 = State_Chm%SpcData(id_Hg2)%Info%EmMW_g
         ENDIF

         ! Pre-save the molecuar weight of HgP [s]
         IF ( id_HgP > 0 ) THEN
            EmMW_g_HgP = State_Chm%SpcData(id_HgP)%Info%EmMW_g
         ENDIF

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Archive diagnostics.  For better efficiency, place everything
      ! within a single OpenMP parallel loop, so that we can enclose
      ! more work within the parallel region than the prior code.
      !=================================================================
      IF ( Input_Opt%ND03 > 0 ) THEN

      ! Loop over advected species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED                             )
!$OMP+PRIVATE( I, J, L, N, NA, EmMw_g, P0, Spc_VV )
!$OMP+SCHEDULE( DYNAMIC, 4                        )
      DO NA = 1, State_Chm%nAdvect

         ! Initialize
         N      = State_Chm%Map_Advect(NA)          ! Species ID #
         EmMW_g = State_Chm%SpcData(N)%Info%EmMW_g  ! Emitted MW (g)

         ! Loop over grid boxes
         DO L = 1, State_Grid%NZ
         DO J = 1, State_Grid%NY
         DO I = 1, State_Grid%NX

            ! Initialize
            P0      = 0.0_fp
            Spc_VV  = 0.0_fp

            !===========================================================
            ! ND03: Hg speciality simulation diagnostics
            !===========================================================
            IF ( NA == 1 ) THEN

               !--------------------------------------------------------
               ! Reactive gaseous mercury [pptv]
               !--------------------------------------------------------
               IF ( id_Hg2 > 0 ) THEN
                  AD03_RGM(I,J,L) = AD03_RGM(I,J,L) 
     &                            + State_Chm%Species(I,J,L,id_Hg2)
     &                            * ( AIRMW / EmMw_g_Hg2 * 1e+12_fp )
               ENDIF

               !--------------------------------------------------------
               ! Particulate bound mercury [pptv]
               !--------------------------------------------------------
               IF ( id_HgP > 0 ) THEN
                  AD03_PBM(I,J,L) = AD03_PBM(I,J,L) 
     &                            + State_Chm%Species(I,J,L,id_HgP)
     &                            * ( AIRMW / EmMW_g_HgP * 1e+12_fp )
               ENDIF
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      ENDIF
#endif
      END SUBROUTINE DIAG1
!EOC
