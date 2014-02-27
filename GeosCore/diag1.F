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
      USE DAO_MOD,            ONLY : IS_ICE, IS_WATER, IS_LAND
      USE DIAG_MOD,           ONLY : AD30, AD31, AD33, AD35, AD45, AD54 
      USE DIAG_MOD,           ONLY : AD47, AD67, AD68, AD69, LTOTH
      USE DIAG_MOD,           ONLY : AD57
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE TIME_MOD,           ONLY : ITS_TIME_FOR_CHEM
      USE TRACER_MOD,         ONLY : XNUMOLAIR
      USE TRACERID_MOD,       ONLY : IDTO3
      USE TRACERID_MOD,       ONLY : ID_HG2, ID_HGP, ID_Hg_TOT
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
      USE DIAG03_MOD,         ONLY : AD03_RGM, AD03_PBM, ND03   

#if   defined( APM )
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! Diagnostic arrays & parameters
      USE CMN_GCTM_MOD             ! Physical constants

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
!
! !REMARKS:
!  For a complete list of GEOS-Chem diagnostics, please see this web page:
!     http://acmg.seas.harvard.edu/geos/doc/man/appendix_5.html
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: AVGW_ALLOCATED, IS_FULLCHEM, IS_CHEM
      INTEGER :: I, J, K, L, N, NN, IREF, JREF, LN45, N_TRC
      REAL*8  :: FDTT, XLOCTM, AREA_M2
#if   defined( APM ) 
      REAL*8  :: STT_VV(IIPAR,JJPAR,LLPAR,
     &                  Input_Opt%N_TRACERS+N_APMTRA)
#else
      REAL*8  :: STT_VV(IIPAR,JJPAR,LLPAR,Input_Opt%N_TRACERS)
#endif 

      !FOR ND57
      ! For FP potential temperature (hotp 7/31/09)
      REAL*8  :: P0

      ! For fields from Input_Opt
      LOGICAL :: ITS_A_FULLCHEM_SIM
      INTEGER :: N_TRACERS
      REAL*8  :: TCVV(Input_Opt%N_TRACERS)

      !=================================================================
      ! DIAG1 begins here!
      !=================================================================
      
      ! Assume success
      RC                 =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      ITS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
      N_TRACERS          = Input_Opt%N_TRACERS
      TCVV               = Input_Opt%TCVV(1:N_TRACERS)

      ! Pick the proper # of tracers for diagnostics
#if   defined( APM )
      N_TRC = N_TRACERS + N_APMTRA    ! G-C advected tracers + APM tracers
#else 
      N_TRC = N_TRACERS               ! G-C advected tracers
#endif

      ! Is it a fullchem run?
      IS_FULLCHEM = ITS_A_FULLCHEM_SIM
      IS_CHEM     = ITS_TIME_FOR_CHEM()
      
      ! Compute conc. in mixing ratio for ND35, ND45, ND47 diagnostics 
      IF ( ND35 > 0 .or. ND45 > 0 .or. ND47 > 0 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRC
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            STT_VV(I,J,L,N) = 
     &         MAX( State_Chm%Tracers(I,J,L,N) * TCVV(N) /
     &              State_Met%AD(I,J,L), 0d0 )
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

! H Amos, 25 Oct 2011--------------------------------------------------

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
          AD03_RGM(I,J,L) = AD03_RGM(I,J,L)+
     &                      State_Chm%Tracers(I,J,L,ID_Hg2(ID_Hg_tot)) *
     &                      TCVV(ID_Hg2(ID_Hg_tot)) /
     &                      State_Met%AD(I,J,L) *1D12 

          ! Reactive particulate mercury, PBM
          AD03_PBM(I,J,L) = AD03_PBM(I,J,L) + 
     &                      State_Chm%Tracers(I,J,L,ID_HgP(ID_Hg_tot)) *
     &                      TCVV(ID_HgP(ID_Hg_tot)) / 
     &                      State_Met%AD(I,J,L) *1D12
                                                       
        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO 

      ENDIF
!-------------------------------------------------------------------------

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
            AD31(I,J,L) = AD31(I,J,L) + GET_PEDGE( I, J, L ) 
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

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
     &              STT_VV(I,J,L,N) * LTOTH(I,J)
            ENDDO
            ENDDO
            ENDDO
         ENDDO   
!$OMP END PARALLEL DO 
      ENDIF

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

            P0          = GET_PEDGE(I,J,1)

            AD57(I,J,L) = AD57(I,J,L) + State_Met%T(I,J,L) *
     &                  ( P0 / GET_PCENTER( I, J, L ) )**0.286

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================  
      ! ND67: Store PBL top pressure [hPa]
      !=================================================================  
      IF ( ND67 > 0 ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! PBL is in [m], use hydrostatic law to get [hPa]
            AD67(I,J,13) = AD67(I,J,13) + 
     &                   ( GET_PEDGE(I,J,1) * EXP( -State_Met%PBLH(I,J)
     &                   / SCALE_HEIGHT ) )

         ENDDO
         ENDDO
      ENDIF

      !================================================================= 
      ! ND68: Quantity 1: BXHEIGHT(I,J,L) in meters 
      !       Quantity 2: AD(I,J,L)       in kg
      !       Quantity 3: AVGW(I,J,L)     in v/v
      !       Quantity 4: N_AIR(I,J,L)    in molecules air / m3
      !
      ! NOTE: AVGW is now an allocatable array from "dao_mod.f"
      !================================================================= 
      IF ( ND68 > 0 ) THEN

         ! Set a flag for whether AVGW is allocated or not
         AVGW_ALLOCATED = ASSOCIATED( State_Met%AVGW )

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LD68
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD68(I,J,L,1) = AD68(I,J,L,1) + State_Met%BXHEIGHT(I,J,L)
            AD68(I,J,L,2) = AD68(I,J,L,2) + State_Met%AD(I,J,L)
            AD68(I,J,L,4) = AD68(I,J,L,4) + State_Met%AIRDEN(L,I,J) *
     &                      XNUMOLAIR

            ! Make sure AVGW is now allocated (bmy, 9/25/01)
            IF ( AVGW_ALLOCATED ) THEN
               AD68(I,J,L,3) = AD68(I,J,L,3) + State_Met%AVGW(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

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
      
      END SUBROUTINE DIAG1
!EOC
