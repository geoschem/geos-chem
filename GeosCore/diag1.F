!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: diag1
!
! !DESCRIPTION: Subroutine DIAG1 accumulates diagnostic quantities on every 
!  dynamic timestep.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG1( am_I_Root, Input_Opt,
     &                  State_Met, State_Chm, RC )
!
! !USES:
!
      ! References to F90 modules
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
      USE DAO_MOD,            ONLY : IS_ICE, IS_WATER, IS_LAND
      USE DIAG_MOD,           ONLY : AD33, AD35, AD54 
      USE DIAG_MOD,           ONLY : AD47
#if defined( BPCH_DIAG )
      USE DIAG_MOD,           ONLY : AD30, AD31, AD45, AD57 
      USE DIAG_MOD,           ONLY : AD67, AD68, AD69, LTOTH
#endif
      USE DIAG_MOD,           ONLY : AD71
      USE DIAG_MOD,           ONLY : AD71_DAY,   AD71_HR
      USE DIAG_MOD,           ONLY : AD71_LDAY,  AD71_LHR
      USE DIAG_MOD,           ONLY : AD71_COUNT, AD71_HRCT
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Chm_Mod, ONLY : Ind_
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE Species_Mod,        ONLY : Species
      USE TIME_MOD,           ONLY : ITS_TIME_FOR_CHEM
      USE DIAG03_MOD,         ONLY : AD03_RGM, AD03_PBM, ND03   

      ! For ND71
      USE TIME_MOD,           ONLY : GET_DAY
      USE TIME_MOD,           ONLY : GET_HOUR
      USE ERROR_MOD,          ONLY : SAFE_DIV

#if defined( APM )
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif

      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE PHYSCONSTANTS
      USE PRECISION_MOD
      USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
      USE HCO_ERROR_MOD
      USE CHEMGRID_MOD,       ONLY : GET_TPAUSE_LEVEL
#if defined( NC_DIAG ) 
      USE DIAGNOSTICS_MOD,    ONLY : DiagnSpec, DiagnUpdate_Met
#endif

      IMPLICIT NONE
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
!  22 Jun 2016 - R. Yantosca - Now use IND_() to define id_PASV species ID
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: AVGW_ALLOCATED, IS_FULLCHEM, IS_CHEM
      INTEGER            :: I, J, K, L, N, NN, IREF, JREF, LN45, N_TRC
      REAL(fp)           :: FDTT, XLOCTM, AREA_M2
#if   defined( APM ) 
      REAL(fp), TARGET   :: STT_VV(IIPAR,JJPAR,LLPAR,
     &                             Input_Opt%N_TRACERS+N_APMTRA)
#else 
      REAL(fp), TARGET   :: STT_VV(IIPAR,JJPAR,LLPAR,
     &                             Input_Opt%N_TRACERS+2)
#endif 

      ! For FP potential temperature (hotp 7/31/09)
      REAL(fp)          :: P0

      ! For fields from Input_Opt
      LOGICAL            :: ITS_A_FULLCHEM_SIM
      INTEGER            :: N_TRACERS

      ! For netCDF diagnostics via HEMCO
      REAL(fp), TARGET   :: Landmap(IIPAR, JJPAR)
      INTEGER            :: TPL
      REAL(fp)           :: PRESS
      REAL(fp), POINTER  :: Ptr3D(:,:,:)
      REAL(fp), POINTER  :: Ptr2D(:,:)
      REAL(fp)           :: Arr3D(IIPAR,JJPAR,LLPAR)
      REAL(fp)           :: Arr2D(IIPAR,JJPAR)
      CHARACTER(LEN=255) :: DiagnName

      ! For species ID's
      TYPE(Species), POINTER :: ThisSpc
      LOGICAL,       SAVE    :: FIRST = .TRUE.
      INTEGER,       SAVE    :: id_PASV
      INTEGER,       SAVE    :: id_Hg2
      INTEGER,       SAVE    :: id_HgP

      !=================================================================
      ! DIAG1 begins here!
      !=================================================================
      
      ! Assume success
      RC                 =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      ITS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
      N_TRACERS          = Input_Opt%N_TRACERS

      ! Initialize pointers
      Ptr3D              => NULL()
      Ptr2D              => NULL()
      ThisSpc            => NULL()

      ! Pick the proper # of tracers for diagnostics
#if   defined( APM )
      N_TRC = N_TRACERS + N_APMTRA    ! G-C advected tracers + APM tracers
#else 
      N_TRC = N_TRACERS               ! G-C advected tracers
#endif

      ! Is it a fullchem run?
      IS_FULLCHEM = ITS_A_FULLCHEM_SIM
      IS_CHEM     = ITS_TIME_FOR_CHEM()
      
      !---------------------------------------------
      ! First-time setup
      !---------------------------------------------
      IF ( FIRST ) THEN 

         ! Define species ID's  on the first call
         id_Hg2  = Ind_('Hg2' )
         id_HgP  = Ind_('HgP' )
         id_PASV = Ind_('PASV')

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Compute conc. in mixing ratio for ND35, ND45, ND47 diagnostics 
      IF ( ND35 > 0 .or. ND45 > 0 .or. ND47 > 0 .or. ND71 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, ThisSpc )
         DO N = 1, N_TRC
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Get info about this tracer from the database
            ! NOTE: Assumes 1:1 mapping of tracer and species indexes
            ThisSpc => State_Chm%SpcData(N)%Info

            ! Tracer conc. in v/v dry air
            STT_VV(I,J,L,N) = 
     &         MAX( State_Chm%Tracers(I,J,L,N) * AIRMW / ThisSpc%emMW_g,
     &              0e+0_fp )

            IF ( N == id_PASV ) THEN

               ! Also express tracer conc. in v/v total air (mps, 11/18/15)
               STT_VV(I,J,L,N+1) =
     &            MAX( State_Chm%Tracers(I,J,L,N) *
     &                 ( 1.0e+0_fp - State_Met%SPHU(I,J,L) * 1.e-3_fp ) 
     &                 * State_Met%MOISTMW(I,J,L) 
     &                 / ThisSpc%emMW_g,
     &                 0e+0_fp )


               ! Also express tracer conc. in kg/kg total air (mps, 11/18/15)
               STT_VV(I,J,L,N+2) =
     &            MAX( State_Chm%Tracers(I,J,L,N) *
     &                 ( 1.0e+0_fp - State_Met%SPHU(I,J,L) * 1.e-3_fp ),
     &                 0e+0_fp )

            ENDIF

            ! Free pointer
            ThisSpc => NULL()
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

! H Amos, 25 Oct 2011--------------------------------------------------

#if defined( BPCH_DIAG )
      !-----------------------------------------------!
      ! Write ND03 diagnostics for RGM and PBM (pptv) !
      !-----------------------------------------------!
      IF(ND03 > 0) THEN

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L)
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Reactive gaseous mercury, RGM
            IF ( id_Hg2 > 0 ) THEN
               AD03_RGM(I,J,L) = AD03_RGM(I,J,L) +
     &              State_Chm%Tracers(I,J,L,id_Hg2) * AIRMW /
     &              State_Chm%SpcData(id_Hg2)%Info%emMW_g * 
     &              1e+12_fp
            ENDIF

            ! Reactive particulate mercury, PBM
            IF ( id_HgP > 0 ) THEN
               AD03_PBM(I,J,L) = AD03_PBM(I,J,L) + 
     &              State_Chm%Tracers(I,J,L,id_HgP) * AIRMW /
     &              State_Chm%SpcData(id_HgP)%Info%emMW_g *
     &              1e+12_fp
            ENDIF
                                                       
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO 
      ENDIF

      !================================================================= 
      ! ND30: Land/water/ice flags
      !=================================================================
      IF ( ND30 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( IS_WATER( I, J, State_Met ) )
     &           AD30(I,J) = AD30(I,J) + 0e0
            IF ( IS_LAND ( I, J, State_Met ) )
     &           AD30(I,J) = AD30(I,J) + 1e0
            IF ( IS_ICE  ( I, J, State_Met ) )
     &           AD30(I,J) = AD30(I,J) + 2e0
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !================================================================= 
      ! ND31: Surface pressure diagnostic (PS - PTOP) in hPa
      !=================================================================
      IF ( ND31 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LD31
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD31(I,J,L) = AD31(I,J,L) + State_Met%PEDGE( I, J, L ) 
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
#endif
      !================================================================= 
      ! ND33: Atmospheric column sum of tracer [kg]
      !=================================================================  
      IF ( ND33 > 0 ) THEN   

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRC
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD33(I,J,N) = AD33(I,J,N) + State_Chm%Tracers(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
         ENDDO  
!$OMP END PARALLEL DO

      ENDIF

      !================================================================= 
      ! ND35: 500 HPa fields.  
      !
      ! NOTES: 
      ! (1 ) Use level 9 for both GEOS-1 and GEOS-STRAT.
      !       They are both close to 500 hPa (bmy, 4/7/99)  
      !================================================================= 
      IF ( ND35 > 0 ) THEN
         L = 9
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, N )
         DO N = 1, N_TRACERS
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD35(I,J,N) = AD35(I,J,N) + STT_VV(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

#if defined( BPCH_DIAG )
      !================================================================= 
      ! ND45: Tracer (V/V) at level 1 to level LD45, averaged over 
      !       the time period OTH_HR1 to OTH_HR2.  
      !
      !  NOTES: 
      !  (1 ) AD45 array replaces the AIJ array for this diagnostic 
      !        (bmy, 3/22/99)
      !  (2 ) Add parallel loop over tracers (bmy, 5/4/00)
      !  (3 ) Use LTO3 and not LTOTH for O3. (ccc, 7/20/09)
      !  (4 ) Removed code for storing pure O3 as N_TRACERS+1 because O3
      !       is now a tracer (mpayer, 3/14/13)
      !================================================================= 
      IF ( ND45 > 0 ) THEN
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRC
            DO L = 1, LD45  
            DO J = 1, JJPAR 
            DO I = 1, IIPAR
               AD45(I,J,L,N) = AD45(I,J,L,N) + 
     &                         STT_VV(I,J,L,N) * LTOTH(I,J)

               IF ( N == id_PASV ) THEN

                  ! Also save passive tracer in v/v total air (mps, 11/18/15)
                  AD45(I,J,L,N+1) = AD45(I,J,L,N+1) + 
     &                              STT_VV(I,J,L,N+1) * LTOTH(I,J)

                  ! Also save passive tracer in kg/kg total air (mps, 11/18/15)
                  AD45(I,J,L,N+2) = AD45(I,J,L,N+2) + 
     &                              STT_VV(I,J,L,N+2) * LTOTH(I,J)

               ENDIF

            ENDDO
            ENDDO
            ENDDO
         ENDDO   
!$OMP END PARALLEL DO
      ENDIF
#endif

      !================================================================= 
      ! ND47: Tracer (V/V) at level 1 to level LD45, 
      ! averaged from 0-24 hours
      !
      !  (1 ) Added parallel loop over tracers (bmy, 5/4/00)
      !  (2 ) Removed code for storing pure O3 as N_TRACERS+1 because O3
      !       is now a tracer (mpayer, 3/14/13)
      !================================================================= 
      IF ( ND47 > 0 ) THEN
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRC
            DO L = 1, LD47
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               AD47(I,J,L,N) = AD47(I,J,L,N) + STT_VV(I,J,L,N) 

               IF ( N == id_PASV ) THEN

                  ! Also save passive tracer in v/v total air (mps, 11/18/15)
                  AD47(I,J,L,N+1) = AD47(I,J,L,N+1) + STT_VV(I,J,L,N+1)

                  ! Also save passive tracer in kg/kg total air (mps, 11/18/15)
                  AD47(I,J,L,N+2) = AD47(I,J,L,N+2) + STT_VV(I,J,L,N+2)

               ENDIF

            ENDDO   
            ENDDO
            ENDDO
         ENDDO 
!$OMP END PARALLEL DO   
      ENDIF

      !================================================================= 
      ! ND54: Count time the box was tropospheric
      !================================================================= 
      IF ( ND54 > 0 ) THEN

            DO L = 1, LD54
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) )
     &              AD54(I,J,L) = AD54(I,J,L) + 1.
            ENDDO   
            ENDDO
            ENDDO
            
      ENDIF
#if defined( BPCH_DIAG )
      !================================================================= 
      ! ND57: Potential temperature
      !================================================================= 
      IF ( ND57 > 0 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, P0 )
         DO L = 1, LD57
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            P0          = State_Met%PEDGE(I,J,1)

            AD57(I,J,L) = AD57(I,J,L) + State_Met%T(I,J,L) *
     &                    ( P0 / State_Met%PMID( I, J, L ) )**0.286

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
#endif
#if defined( BPCH_DIAG )
      !=================================================================  
      ! ND67: Store PBL top pressure [hPa]
      !=================================================================  
      IF ( ND67 > 0 ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! PBL is in [m], use hydrostatic law to get [hPa]
            AD67(I,J,13) = AD67(I,J,13) + 
     &                     ( State_Met%PEDGE(I,J,1) 
     &                     * EXP( -State_Met%PBLH(I,J)
     &                     / SCALE_HEIGHT ) )
         ENDDO
         ENDDO
      ENDIF
#endif
#if defined( BPCH_DIAG )
      !================================================================= 
      ! ND68: 1:  BXHEIGHT     - box height [m] 
      !       2:  AD           - dry air mass [kg]
      !       3:  AVGW         - mol H2O / mol dry air [v/v]
      !       4:  N_DRYAIR     - dry air number density [molecules/m3]
      !       5:  MAIRDEN      - moist air mass density [kg/m3]
      !       6:  T            - temperature [K]
      !       7:  TV           - virtual temperature [K]
      !       8:  PMID         - arithmetic mean pressure [hPa]
      !       9:  PMEAN        - altitude-weighted mean pressure [hPa]
      !       10: PEDGE        - level pressure (bottom edge) [hPa]
      !       11: RH           - relative humidity [%]
      !       12: PMID_DRY     - dry air partial pressure at PMID [hPa]
      !       13: PMEAN_DRY    - dry air partial pressure at PMEAN [hPa]
      !================================================================= 
      IF ( Input_Opt%ND68 > 0 ) THEN

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LD68
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD68(I,J,L,1)  = AD68(I,J,L,1)  + State_Met%BXHEIGHT(I,J,L)
            AD68(I,J,L,2)  = AD68(I,J,L,2)  + State_Met%AD(I,J,L)
            AD68(I,J,L,3)  = AD68(I,J,L,3)  + State_Met%AVGW(I,J,L)
            AD68(I,J,L,4)  = AD68(I,J,L,4)  + State_Met%AIRDEN(I,J,L) 
     &                                        * XNUMOLAIR
            AD68(I,J,L,5)  = AD68(I,J,L,5)  + State_Met%MAIRDEN(I,J,L) 
            AD68(I,J,L,6)  = AD68(I,J,L,6)  + State_Met%T(I,J,L)
            AD68(I,J,L,7)  = AD68(I,J,L,7)  + State_Met%TV(I,J,L)
            AD68(I,J,L,8)  = AD68(I,J,L,8)  + State_Met%PMID(I,J,L)
            AD68(I,J,L,9)  = AD68(I,J,L,9)  + State_Met%PMEAN(I,J,L)
            AD68(I,J,L,10) = AD68(I,J,L,10) + State_Met%PEDGE(I,J,L)
            AD68(I,J,L,11) = AD68(I,J,L,11) + State_Met%RH(I,J,L)
            AD68(I,J,L,12) = AD68(I,J,L,12) + State_Met%PMID_DRY(I,J,L)
            AD68(I,J,L,13) = AD68(I,J,L,13) + State_Met%PMEAN_DRY(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
#endif
#if defined( BPCH_DIAG )
      !================================================================= 
      ! ND69: Grid box surface areas [m2]
      !
      ! NOTE: Only save areas on the first timestep since the
      !       grid box surface areas are a time-invariant field.
      !================================================================= 
      IF ( ND69 > 0 ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD69(I,J,1) = GET_AREA_M2( I, J, 1 )
         ENDDO
         ENDDO
      ENDIF
#endif
      
      !================================================================= 
      ! ND71: Tracer (V/V) at level 1, daily maximum hourly avg.
      !================================================================= 
      IF ( ND71 > 0 ) THEN

         IF ((GET_DAY().ne.AD71_LDAY).and.(AD71_LDAY.ge.0)) THEN
            ! It's a new day, with no diagnostics written
            AD71(:,:,:) = AD71(:,:,:) + AD71_DAY(:,:,:)
            AD71_DAY(:,:,:) = 0e0
            AD71_COUNT = AD71_COUNT + 1
         ENDIF

         IF ((GET_HOUR().ne.AD71_LHR).and.(AD71_LHR.ge.0)) THEN
            ! New hour - get average of last data
!            AD71_HR = AD71_HR / AD71_HRCT
            ! Set daily max
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
            DO N = 1, N_TRC
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  AD71_DAY(I,J,N) = MAX(AD71_DAY(I,J,N),
     &                     (AD71_HR(I,J,N)/(AD71_HRCT)))
               ENDDO   
               ENDDO
            ENDDO 
!$OMP END PARALLEL DO
            AD71_HR = 0e0
            AD71_HRCT = 0
         ENDIF

         ! Set the tracking variables
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRC
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               AD71_HR(I,J,N) = AD71_HR(I,J,N) + STT_VV(I,J,1,N)
            ENDDO   
            ENDDO
         ENDDO
!$OMP END PARALLEL DO   

         AD71_LDAY = GET_DAY()
         AD71_LHR = GET_HOUR()
         AD71_HRCT = AD71_HRCT + 1
      ENDIF

!##############################################################################
!###                                                                        ###
!###   NEW GEOS-Chem DIAGNOSTICS VIA HEMCO LINKED LIST --> netCDF output    ###
!###                                                                        ###
!##############################################################################

#if defined( NC_DIAG )

      !=================================================================
      ! ND30: Landmap diagnostics (netCDF diagnostic structure)
      !=================================================================
      IF ( Input_Opt%ND30 > 0 ) THEN

         ! Construct diagnostic name
         ! (better to have function to get it since already defined)
         DiagnName = 'LANDMAP'

         ! Compute landmap
         IF ( ND30 > 0 ) THEN  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               IF ( IS_WATER( I, J, State_Met ) )
     &           Landmap(I,J) = 0e0
               IF ( IS_LAND ( I, J, State_Met ) )
     &           Landmap(I,J) = 1e0
               IF ( IS_ICE  ( I, J, State_Met ) )
     &           Landmap(I,J) = 2e0
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF

         ! Point to the landmap array
         Ptr2D => Landmap
    
         ! Update the landmap diagnostic
         CALL Diagn_Update( am_I_Root,
     &                      cName   = TRIM( DiagnName ),
     &                      Array2D = Ptr2D,
     &                      RC      = RC )
    
         ! Free the pointer before error handling
         Ptr2D => NULL()
    
         ! Stop with error if the diagnostics were unsuccessful
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL ERROR_STOP(
     &           'Cannot update landmap diagnostic (ND30)',
     &           'DIAG1 (diag1.F)')
         ENDIF
 
      ENDIF

      !=================================================================
      ! ND45: Tracer concentrations (netCDF diagnostic structure)
      !=================================================================
      IF ( Input_Opt%ND45 > 0 ) THEN

         ! Loop over all tracers
         DO N = 1, Input_Opt%N_TRACERS

            ! If this tracer number NN is scheduled for output in input.geos,
            ! then archive the latest concentration into the diagnostic structure
            IF ( ANY( Input_Opt%TINDEX(45,:) == N ) ) THEN

               IF ( N == id_PASV ) THEN

                  !-------------------------------------
                  ! Passive tracer in v/v dry air
                  !-------------------------------------

                  ! Construct diagnostic name
                  DiagnName = 'TRACER_CONC_' // 
     &                         TRIM( Input_Opt%TRACER_NAME(N) ) //
     &                        '_VVDry'

                  ! Point to the proper slice of the DEPSAV array
                  Ptr3D => STT_VV(:,:,:,N)

                  ! Update the deposition velocity diagostics
                  CALL Diagn_Update( am_I_Root,
     &                               cName   = TRIM( DiagnName ),
     &                               Array3D = Ptr3D,
     &                               RC      = RC )

                  !-------------------------------------
                  ! Passive tracer in v/v total air
                  !-------------------------------------

                  ! Construct diagnostic name
                  DiagnName = 'TRACER_CONC_' // 
     &                         TRIM( Input_Opt%TRACER_NAME(N) ) //
     &                        '_VVTotal'

                  ! Point to the proper slice of the DEPSAV array
                  Ptr3D => STT_VV(:,:,:,N+1)

                  ! Update the deposition velocity diagostics
                  CALL Diagn_Update( am_I_Root,
     &                               cName   = TRIM( DiagnName ),
     &                               Array3D = Ptr3D,
     &                               RC      = RC )

                  !-------------------------------------
                  ! Passive tracer in kg/kg total air
                  !-------------------------------------

                  ! Construct diagnostic name
                  DiagnName = 'TRACER_CONC_' // 
     &                         TRIM( Input_Opt%TRACER_NAME(N) ) //
     &                        '_KgKgTotal'

                  ! Point to the proper slice of the DEPSAV array
                  Ptr3D => STT_VV(:,:,:,N+2)

                  ! Update the deposition velocity diagostics
                  CALL Diagn_Update( am_I_Root,
     &                               cName   = TRIM( DiagnName ),
     &                               Array3D = Ptr3D,
     &                               RC      = RC )

               ELSE

                  ! Construct diagnostic name
                  DiagnName = 'TRACER_CONC_' // 
     &                         TRIM( Input_Opt%TRACER_NAME(N) )

                  ! Point to the proper slice of the DEPSAV array
                  Ptr3D => STT_VV(:,:,:,N)

                  ! Update the deposition velocity diagostics
                  CALL Diagn_Update( am_I_Root,
     &                               cName   = TRIM( DiagnName ),
     &                               Array3D = Ptr3D,
     &                               RC      = RC )

               ENDIF

               ! Free the pointer
               Ptr3D => NULL()

               ! STOP wieth error if the diagnostics were unsuccessful
               IF ( RC /= HCO_SUCCESS ) THEN
                  CALL ERROR_STOP(
     &                 'Cannot update TRACER_CONC diag diagnostics',
     &                 'DIAG1 (diag1.F)')
               ENDIF
            ENDIF
         ENDDO

         IF ( DiagnSpec ) THEN

            ! Loop over all species
            DO N = 1, State_Chm%nSpecies

               ! Get info about this species from the species database
               ThisSpc => State_Chm%SpcData(N)%Info
   
               ! Construct diagnostic name
               IF ( TRIM( ThisSpc%Name ) == '' ) CYCLE
               DiagnName = 'SPECIES_CONC_' // TRIM( ThisSpc%Name ) 
         
               ! Get value in #/cm3
               Ptr3D => State_Chm%Species(:,:,:,N)
   
               ! Update the deposition velocity diagostics
               CALL Diagn_Update( am_I_Root,
     &                            cName   = TRIM( DiagnName ),
     &                            Array3D = Ptr3D,
     &                            RC      = RC )
   
               ! Free the pointers
               Ptr3D   => NULL()
               ThisSpc => NULL()

               ! STOP wieth error if the diagnostics were unsuccessful
               IF ( RC /= HCO_SUCCESS ) THEN
                  CALL ERROR_STOP(
     &                    'Cannot update SPECIES_CONC diag diagnostics',
     &                    'DIAG1 (diag1.F)')
               ENDIF
            ENDDO
   
         ENDIF ! DiagnSpec
      ENDIF !ND45
   

      !=================================================================
      ! Met diagnostics: ND31, ND55, ND56, ND57, ND66, ND67, ND68, ND69
      !=================================================================

      ! For now, use ND68 as indicator that Met diagnostic group is on
      IF ( Input_Opt%ND68 > 0 ) THEN
         CALL DiagnUpdate_Met( am_I_Root, Input_Opt, State_Met, RC ) 
      ENDIF

      ! Formerly, netcdf ND68 diagnostic included:
!      !================================================================= 
!      ! ND68: 1:  BXHEIGHT     - box height [m] 
!      !       2:  AD           - dry air mass [kg]
!      !       3:  AVGW         - mol H2O / mol dry air [v/v]
!      !       4:  AIRDENSITY   - dry air number mass density [kg/m3]
!      !                          (Note this is different units than bpch)
!      !       5:  MAIRDENSITY  - moist air mass density [kg/m3]
!      !       6:  T            - temperature [K]
!      !       7:  TV           - virtual temperature [K]
!      !       8:  PMID         - arithmetic mean pressure [hPa]
!      !       9:  PMEAN        - altitude-weighted mean pressure [hPa]
!      !       10: PEDGE        - level pressure (bottom edge) [hPa]
!      !       11: RH           - relative humidity [%]
!      !       12: PMID_DRY     - dry air partial pressure at PMID [hPa]
!      !       13: PMEAN_DRY    - dry air partial pressure at PMEAN [hPa]
!      !=================================================================   
!         ! There are 13 diagnostics
!         DO I = 1,13
!
!            ! Construct diagnostic name and get proper array
!            SELECT CASE ( I ) 
!               CASE ( 1 ) 
!                  DiagnName = 'BOXHEIGHT'
!                  Ptr3D => State_Met%BXHEIGHT
!               CASE ( 2 ) 
!                  DiagnName = 'AIRMASS'
!                  Ptr3D => State_Met%AD
!               CASE ( 3 ) 
!                  DiagnName = 'AVGW'
!                  Ptr3D => State_Met%AVGW
!               CASE ( 4 ) 
!                  DiagnName = 'AIRDENSITY'
!                  Ptr3D => State_Met%AIRDEN
!               CASE ( 5 ) 
!                  DiagnName = 'MAIRDENSITY'
!                  Ptr3D => State_Met%MAIRDEN
!               CASE ( 6 ) 
!                  DiagnName = 'T'
!                  Ptr3D => State_Met%T
!               CASE ( 7 ) 
!                  DiagnName = 'TV'
!                  Ptr3D => State_Met%TV
!               CASE ( 8 ) 
!                  DiagnName = 'PMID'
!                  Ptr3D => State_Met%PMID
!               CASE ( 9 ) 
!                  DiagnName = 'PMEAN'
!                  Ptr3D => State_Met%PMEAN
!               CASE ( 10 ) 
!                  DiagnName = 'PEDGE'
!                  Ptr3D => State_Met%PEDGE
!               CASE ( 11 ) 
!                  DiagnName = 'RH'
!                  Ptr3D => State_Met%RH
!                CASE ( 12 ) 
!                  DiagnName = 'PMID_DRY'
!                  Ptr3D => State_Met%PMID_DRY
!               CASE ( 13 ) 
!                  DiagnName = 'PMEAN_DRY'
!                  Ptr3D => State_Met%PMEAN_DRY
!            END SELECT
! 
!            ! Update diagnostics
!            IF ( ASSOCIATED(Ptr3D) ) THEN 
!               CALL Diagn_Update( am_I_Root,
!     &                            cName   = TRIM( DiagnName ),
!     &                            Array3D = Ptr3D,
!     &                            RC      = RC )
!
!               ! Free the pointer
!               Ptr3D => NULL()
!   
!               ! STOP with error if the diagnostics were unsuccessful
!               IF ( RC /= HCO_SUCCESS ) THEN
!                  CALL ERROR_STOP(
!     &                 'Cannot update diagnostics '//TRIM(DiagnName),  
!     &                 'DIAG1 (diag1.F)')
!               ENDIF
!            ENDIF   
!         ENDDO
!
!         ! Now also add TROPP_PRESSURE and TROPP_LEVEL
!         DiagnName = 'TROPP_PRESSURE'
!         Ptr2D => State_Met%TROPP
!         CALL Diagn_Update( am_I_Root, cName=DiagnName, 
!     &      Array2D=Ptr2D, RC=RC )
!         Ptr2D => NULL()
!         
!         DiagnName = 'TROPP_LEVEL'
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, TPL )
!         DO J=1,JJPAR
!         DO I=1,IIPAR
!            TPL        = GET_TPAUSE_LEVEL( I, J, State_Met )
!            Arr2D(I,J) = TPL
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO   
!         CALL Diagn_Update( am_I_Root, cName=DiagnName, 
!     &      Array2D=Arr2D, RC=RC )
! 
!      ENDIF ! ND68
#endif

      END SUBROUTINE DIAG1
!EOC
