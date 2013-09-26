! $Id: emissdr.f,v 1.6 2010/03/15 19:33:24 ccarouge Exp $
      SUBROUTINE EMISSDR( am_I_Root, Input_Opt, 
     &                    State_Met, State_Chm, STATUS )
!
!******************************************************************************
!  Subroutine EMISSDR computes emissions for the full chemistry simulation
!  Emissions are stored in various arrays, which are then passed to the
!  SMVGEAR solver via routine "setemis.f". (bmy, 10/8/98, 2/25/10)
!
!  NOTES:
!  (1 ) Now accounts for seasonal NOx emissions, and multi-level NOx 
!        emissions. (bmy, 10/8/98)
!  (2 ) Surface NOx and 100m NOx are now placed into the correct sigma level.
!        (bmy, 10/8/98)
!  (3 ) Eliminate GISS-Specific code (bmy, 3/15/99)
!  (4 ) Now includes monoterpenes for ACETONE emissions (bdf, 4/8/99)
!  (5 ) Now uses allocatable arrays for ND29, ND36, and ND46 diagnostics.
!        Also made some cosmetic changes, and updated comments (bmy, 3/16/00)
!  (6 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (7 ) Add reference to BIOBURN in "biomass_mod.f" and to BIOFUEL_BURN
!        in "biofuel_mod.f".  Also remove references to BURNEMIS and TWOODIJ.
!        (bmy, 9/12/00)
!  (8 ) Remove reference to "biomass.h" -- that is replaced by F90 module
!        "biomass_mod.f". (bmy, 9/25/00)
!  (9 ) Remove obsolete code from 9/12/00 and 9/25/00 (bmy, 12/21/00)
!  (10) Add CO source from monoterpenes (bnd, bmy, 12/21/00)
!  (11) Add monoterpene source to the ND46 diagnostic.  Renamed EMXX to
!        EMIS (for Isoprene), and commented out Larry Horowitz's "special
!        cases".  Also made some cosmetic changes. (bmy, 1/2/01)
!  (12) Added CO source from CH3OH oxidation (bmy, 1/3/01)
!  (13) Removed obsolete code from 1/2/01 (bmy, 3/15/01)
!  (14) Now initialize GEMISNOX2.  Also updated comments. (bdf, bmy, 6/15/01)
!  (15) Now references routines from "acetone_mod.f" for the biogenic
!        emission of acetone into the SMVGEAR arrays.  Now use 
!        EMISRR(I,J,IDEACET) to archive ND46, since the biogenic
!        acetone emissions are now computed in this array.  Also define
!        XNUMOL_C so as not to rely on IDTISOP being defined.  Also add
!        LASTMONTH variable to flag when we change month. (bmy, 9/4/01)
!  (16) Now reference AIREMISS from "aircraft_nox_mod.f" (bmy, 2/14/02)
!  (17) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE.  Also removed
!        obsolete, commented-out code. (bmy, 6/25/02)
!  (18) Now references IDTNOX, etc. from "tracerid_mod.f".  Now references
!        SUNCOS from "dao_mod.f".  Now make FIRSTEMISS a local SAVEd variable
!        instead of an argument. (bmy, 11/15/02)
!  (19) Now replaced DXYP(JREF)*1d4 with GET_AREA_CM2 from "grid_mod.f".
!        Now remove MONTH from call to BIOBURN.  Now use functions GET_MONTH,
!        GET_LOCALTIME, GET_ELAPSED_MIN, GET_TS_EMIS, GET_LOCALTIME from 
!        "time_mod.f".  Now use functions GET_XOFFSET and GET_YOFFSET from 
!        "grid_mod.f". (bmy, 2/11/03)
!  (20) Now pass I, J to EMISOP, EMISOP_GRASS, EMISOP_MB (bmy, 12/9/03)
!  (21) Now references EMLIGHTNING from "lightning_nox_mod.f" (bmy, 4/14/04)
!  (22) Now references "logical_mod.f".  Now replaced LFOSSIL with LANTHRO.
!        (bmy, 7/20/04)
!  (23) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (24) Now can use MEGAN inventory for biogenic VOCs.  Now references
!        "megan_mod.f" (bmy, tmf, 10/25/05)
!  (25) Now call EMLIGHTNING_NL from "lightning_nox_nl_mod.f" for GEOS-4 so
!        that we can use the new near-land lightning formulation. 
!        (ltm, bmy, 5/11/06)
!  (26) Added switch for BIOGENIC emissions.  Now revert to single 
!        lightning_nox_mod.f.  Remove reference to old GEMISNOX array; this
!        has been replaced by module arrays EMIS_LI_NOx and EMIS_AC_NOx, in 
!        order to avoid common block errors. (ltm, bmy, phs, 10/3/07)
!  (27) Add biogenic emission of MONX and C2H4. (tmf, 1/20/09)
!  (28) Add a switch, LMEGANMONO, for monoterpene and MBO emission to choose
!       between MEGAN and GEIA inventory indenpendantly from Isoprene.
!      This is because although we did, at one time, got monoterpene
!      emission factors from Alex Guenther, he never published it.
!      So there is no reference for it. Use of those emission factors have
!      caused much confusion among users, so we should probably not use it 
!      for the time being. (tmf, 1/20/09)
!  (29) Update to MEGAN v2.1 (mpb, 11/20/09)
!  (29) Move XLTMMP to module MEGANUT_MOD (ccc, 11/20/09)
!  (30) Change arguments order in some MEGAN functions for coherence
!       (ccc, 11/30/09)
!  (31) Remove reference to obsolete embedded chemistry stuff in "CMN" 
!       (bmy, 2/25/10)
!  (32) Add a NOx fertilizer switch and a scaling factor for ISOP emissions
!       (fp, 6/09)
!  05 Oct 2011 - R. Yantosca - Now use SUNCOS30 array, which has the cos(SZA)
!                              defined at 30 mins after the GMT hour, in order
!                              to be consistent w/ the chemistry.
!  07 Oct 2011 - R. Yantosca - Rename SUNCOS30 to SUNCOS_MID, which is the
!                              cos(SZA) at the midpt of the chemistry timestep
!  06 Dec 2011 - E. Fischer  - Change biogenic acetone emissions to MEGAN
!  08 Dec 2011 - M. Payer    - Remove GEIA biogenic emissions option
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  22 Oct 2012 - R. Yantosca - Prevent out-of-bounds errors in ND46 if bromine
!                              emissions are turned off (e.g. for tagged CO)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  27 Nov 2012 - R. Yantosca - Replace SUNCOS_MID with State_Met%SUNCOSmid
!  13 Dec 2012 - R. Yantosca - Now remove reference to obsolete CMN_DEP_mod.F
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  22 Jul 2013 - M. Sulprizio- Now copy LRCPAIR from Input_Opt
!  31 Jul 2013 - M. Sulprizio- Now use aircraft emissions from the AEIC
!                              inventory (S. Eastham)
!  13 Aug 2013 - M. Sulprizio- Update MEGAN biogenics to include sesquiterpenes
!                              for updated SOA simulation (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      ! References to F90 modules
      USE ACETONE_MOD,        ONLY : EMISS_BIOACET, OCEAN_SOURCE_ACET
      USE AEIC_MOD,           ONLY : EMIT_AEIC
      USE BIOFUEL_MOD,        ONLY : BIOFUEL_BURN
      USE BROMOCARB_MOD,      ONLY : EMIS_CHBr3, EMIS_CH2Br2
      USE BROMOCARB_MOD,      ONLY : Br_SCALING
      USE CANOPY_NOX_MOD,     ONLY : GET_CANOPY_NOX
      USE DIAG_MOD,           ONLY : AD29,          AD46
      USE DIAG_MOD,           ONLY : AD32_fe,       AD32_so
      USE EMISSIONS_MOD,      ONLY : ISOP_SCALING
      USE GET_NDEP_MOD,       ONLY : GET_DEP_N 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE GRID_MOD,           ONLY : GET_YMID
      USE GRID_MOD,           ONLY : GET_XMID
      USE GRID_MOD,           ONLY : GET_XOFFSET,   GET_YOFFSET
      USE LIGHTNING_NOX_MOD,  ONLY : EMLIGHTNING
      USE MEGAN_MOD,          ONLY : GET_EMISOP_MEGAN
      USE MEGAN_MOD,          ONLY : GET_EMMBO_MEGAN
      USE MEGAN_MOD,          ONLY : GET_EMMONOT_MEGAN
      USE MEGAN_MOD,          ONLY : GET_EMTERP_MEGAN
      USE MEGANUT_MOD,        ONLY : XLTMMP
      USE MODIS_LAI_MOD
      USE RCP_MOD,            ONLY : RCP_AIREMISS! (cdh, 10/14/11)
      USE SOIL_NOX_MOD,       ONLY : SOIL_NOX_EMISSION   
      USE TIME_MOD,           ONLY : GET_MONTH,     GET_TAU
      USE TIME_MOD,           ONLY : GET_TS_EMIS,   GET_LOCALTIME
      USE TIME_MOD,           ONLY : GET_DAY_OF_YEAR
      USE TRACERID_MOD,       ONLY : IDEACET,       IDTISOP,   IDEISOP   
      USE TRACERID_MOD,       ONLY : IDECO,         IDEPRPE,   NEMANTHRO 
      USE TRACERID_MOD,       ONLY : IDEMONX, IDEC2H4
      USE TRACERID_MOD,       ONLY : IDTMONX, IDTC2H4
      USE TRACERID_MOD,       ONLY : IDECHBr3, IDECH2Br2
 
      USE CMN_SIZE_MOD             ! Size parameters
      USE COMODE_LOOP_MOD          ! IVERT?
      USE CMN_O3_MOD               ! Emissions arrays
      USE CMN_DIAG_MOD             ! Diagnostic arrays and switches
      USE CMN_NOX_MOD              ! GEMISNOX2
      USE COMMSOIL_MOD             ! Soil wetness variables

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
      INTEGER,        INTENT(OUT)   :: STATUS      ! Success or failure?
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE. 
      LOGICAL                :: NO_TAGCO
      INTEGER                :: I, J, L, N, IJLOOP 
      INTEGER                :: I0, J0, IOFF, JOFF, IREF, JREF
      INTEGER                :: NDAY,JSCEN,NN
      INTEGER, SAVE          :: LASTMONTH = -99
      REAL*8                 :: DTSRCE,  XLOCTM,   EMIS,  TMMP
      REAL*8                 :: BIOSCAL, AREA_CM2, EMMO,  ACETSCAL
      REAL*8                 :: TMPVAL,  EMMB,     GRASS, BIO_ACET
      REAL*8                 :: CONVERT(NVEGTYPE), GMONOT(NVEGTYPE)
      REAL*8                 :: SC, PDF, PDR
      INTEGER                :: DOY, MONTH, RC, M
      REAL*8                 :: FERTDIAG
      REAL*8                 :: TS_EMIS
      REAL*8                 :: TOTALPULSE
      REAL*8                 :: SOILFRT
      REAL*8                 :: DEP_FERT, TEMP_FERT

      ! Add biogenic emission scale factor for ethene (tmf, 1/13/06)
      REAL*8                 :: BIOSCALEC2H4

      ! Molecules C / kg C
      REAL*8,  PARAMETER     :: XNUMOL_C = 6.022d+23 / 12d-3 

      REAL*8                 :: AREA_M2 !jpp 9/10/07

      ! Molecular weight of CHBr3 [kg/mol], (jpp, 9/10/07)
      REAL*8,  PARAMETER     :: MWT_CHBr3  = 2.53d-1

      ! Molecular weight of CH2Br2 [kg/mol], (jpp, 9/16/07)
      REAL*8,  PARAMETER     :: MWT_CH2Br2 = 1.74d-1

      ! Add in new monoterpene species (mpb,2009)
      REAL*8                 :: APINE  , BPINE , LIMON , SABIN 
      REAL*8                 :: MYRCN  , CAREN , OCIMN

      ! SOAupdate: Sesquiterpenes (hotp 3/1/10)
      REAL*8                 :: FARN, BCAR, OSQT
      REAL*8                 :: OMTP ! other monoterpenes 3/5/10

      ! For fields from Input_Opt
      LOGICAL                :: LANTHRO
      LOGICAL                :: LLIGHTNOX
      LOGICAL                :: LSOILNOX
      LOGICAL                :: LFERTILIZERNOX
      LOGICAL                :: LBIOMASS
      LOGICAL                :: LBIONOX
      LOGICAL                :: LBIOFUEL
      LOGICAL                :: LWOODCO   
      LOGICAL                :: LMEGAN
      LOGICAL                :: LMEGANMONO
      LOGICAL                :: LBIOGENIC
      LOGICAL                :: LNLPBL
      LOGICAL                :: LRCPAIR
      LOGICAL                :: LAEIC
      LOGICAL                :: ITS_A_TAGCO_SIM
      INTEGER                :: N_TRACERS
!
!******************************************************************************
!  EMISSDR begins here!
!
!  DTSRCE = emission timestep in seconds
!
!  Call subroutines to set up ISOP and monoterpene emission (first time only!)
!******************************************************************************
!
      ! Assume success
      STATUS           = GIGC_SUCCESS

      ! Copy values from Input_Opt
      LNLPBL           = Input_Opt%LNLPBL
      LANTHRO          = Input_Opt%LANTHRO
      LLIGHTNOX        = Input_Opt%LLIGHTNOX
      LSOILNOX         = Input_Opt%LSOILNOX
      LFERTILIZERNOX   = Input_Opt%LFERTILIZERNOX
      LBIOMASS         = Input_Opt%LBIOMASS
      LBIONOX          = Input_Opt%LBIOMASS
      LBIOFUEL         = Input_Opt%LBIOFUEL
      LWOODCO          = Input_Opt%LBIOFUEL
      LMEGAN           = Input_Opt%LMEGAN
      LMEGANMONO       = Input_Opt%LMEGANMONO
      LBIOGENIC        = Input_Opt%LBIOGENIC
      LNLPBL           = Input_Opt%LNLPBL
      LRCPAIR          = Input_Opt%LRCPAIR
      LAEIC            = Input_Opt%LAEIC
      ITS_A_TAGCO_SIM  = Input_Opt%ITS_A_TAGCO_SIM
      N_TRACERS        = Input_Opt%N_TRACERS

      ! This is not a tagged CO simulation
      NO_TAGCO = ( .not. ITS_A_TAGCO_SIM )

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SOIL EMISSIONS NOX [molecules/cm3/s]

      IF ( LSOILNOX ) THEN           

	 !Zero the Diag49 arrays
         IF ( FIRSTEMISS ) THEN
            INST_SOIL(:,:) = 0.0
            INST_FERT(:,:) = 0.0

            FIRSTEMISS = .FALSE.

         ENDIF

         !***************************************************************
         !First zero the arrays in which emissions will be stored
         !GEMISNOX2(I,J)   = Array which stores NOx emissions from soils.
         !                   Units are [molec NOx/cm3/s].
         GEMISNOX2   = 0d0
         !***************************************************************

         ! Now need to call GET_CANOPY_NOX to break ugly dependency between
         ! drydep and soil NOx emissions. (bmy, 6/22/09) 
         ! Now a function of the new MODIS/Koppen biome map (J.D. Maasakkers)
         CALL GET_CANOPY_NOX( State_Met )
               
         ! Get day of year/month and timestep information
         DOY     = GET_DAY_OF_YEAR() 
         TS_EMIS = GET_TS_EMIS()           
         MONTH   = GET_MONTH()

         ! Loop over each land grid-box, removed loop over landpoints
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, DEP_FERT, SOILFRT, FERTDIAG, RC, IJLOOP)
 	 DO J = 1, JJPAR 
	 DO I = 1, IIPAR

 	    ! Do not calculate soil NOx emissions if there is no soil in
            ! the gridbox
 	    IF (LAND2(I,J,1) .EQ. 1) CYCLE

            IJLOOP = ( (J-1) * IIPAR ) + I 

            ! Get Deposited Fertilizer
            CALL GET_DEP_N( I, J, TS_EMIS, DEP_FERT )
            CALL FLUSH(6)

            ! Get N fertilizer reservoir associated with chemical and
            ! manure fertilizer
            SOILFRT =  SOILFERT(I,J,GET_DAY_OF_YEAR())  

            ! Put in constraint if dry period gt 1 yr, keep at 1yr to
            ! avoid unrealistic pulse
            IF ( DRYPERIOD(I,J) .GT. 8760.) DRYPERIOD(I,J) = 8760.
 
            ! Return NOx emission from soils [molec/cm2/s]
            CALL SOIL_NOX_EMISSION( TS_EMIS,        
     &                              State_Met%TS(I,J),
     &                              State_Met%GWETTOP(I,J),
     &                              SOILFRT,  
     &                              GWET_PREV(I,J), DRYPERIOD(I,J),
     &                              PFACTOR(I,J),   GEMISNOX2(I,J),
     &                              DEP_FERT,       FERTDIAG, 
     &                              CLIM(I,J,:),    LAND2(I,J,:),
     &                              RC,             
     &                              State_Met%SUNCOSmid(I,J),
     &				    State_Met%U10M(I,J),
     &                              State_Met%V10M(I,J), 
     &                              CANOPYNOX(IJLOOP,:), 
     &                              GC_LAI(I,J) )

            ! Archive emissions
            IF ( ND32 > 0 ) THEN
               AD32_so(I,J) = AD32_so(I,J) + GEMISNOX2(I,J)
               AD32_fe(I,J) = AD32_fe(I,J) + FERTDIAG
            ENDIF
            
            INST_SOIL(I,J) = GEMISNOX2(I,J)
            INST_FERT(I,J) = FERTDIAG

            GEMISNOX2(I,J) = GEMISNOX2(I,J)  / 
     &                     ( State_Met%BXHEIGHT(I,J,1) * 100d0 )

         ENDDO            
         ENDDO
!$OMP END PARALLEL DO

         !Zero arrays. This code assumes that ts_chem = ts_emiss
         DRY_HNO3 = 0.d0
         DRY_NH4  = 0.d0
         DRY_NH3  = 0.d0
         DRY_NIT  = 0.d0
         DRY_NO2  = 0.d0
         DRY_PAN  = 0.d0

         WET_HNO3 = 0.d0
         WET_NH4  = 0.d0
         WET_NH3  = 0.d0
         WET_NIT  = 0.d0

      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!******************************************************************************
!  BE CAREFUL WITH USING A WINDOW RELATIVE TO THE EMISSIONS WINDOW
!  NEED SPECIFY THE OFFSET OF THE SUB-WINDOW
!
!  first zero the arrays in which emissions will be stored
!
!  EMISRRN(I,J,L)   = Emission rate of NOx (N = IDENOX ) into box (I,J,L).
!                     Units are [molec NOx/box/s].
!
!  EMISRR(I,J,N)    = Emission rate of tracer N into surface box (I,J,1). 
!                     Units are [molec tracer/box/s].
!
!  GEMISNOX(I,J,L)  = Array which stores NOx emissions from aircraft,
!                     and lightning.  Units are [molec NOx/cm3/s].
!
!  GEMISNOX2(I,J)   = Array which stores NOx emissions from soils.
!                     Units are [molec NOx/cm3/s].
!
!  NOTE: Now use F90 array initialization syntax (bmy, 3/15/99)
!******************************************************************************
!
      ! These need to be initialized on every call
      EMISRRN     = 0d0
      EMISRR      = 0d0
      !GEMISNOX2   = 0d0
      EMISS_BVOC  = 0d0  ! Add BVOC emissions (mpb,2009)
      
      ! Loop over latitudes
      IJLOOP = 0
      DO J = 1, JJPAR
         JREF = J + J0

         ! Loop over longitues
         DO I = 1, IIPAR
            IREF   = I + I0
            IJLOOP = IJLOOP + 1

            ! Compute surface area of grid boxes in cm^2 
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! jpp, 9/10/07: wanted surface area in m^2
            AREA_M2  = AREA_CM2 / (100.d0 ** 2)

            ! Zero biogenic acetone (bmy, 9/14/01)
            BIO_ACET  = 0d0

            ! Use function GET_LOCALTIME to get the local time at lon I
            ! Middle of time step is between 10pm-2am when IHOUR = 1
            IHOUR = NINT( ( GET_LOCALTIME( I, J, 1 ) ) / 4 ) + 1
            IF ( IHOUR .EQ. 7 ) IHOUR = 1

            !=================================================================
            ! attenuate emissions on the weekend ---
            ! scale factors for Saturday/Sunday/Weekday must average out to 1!
            !    JSCEN = 1 Saturday 
            !    JSCEN = 2 Sunday
            !    JSCEN = 3 Weekday 
            !
            ! 1 Jan 1980 and 1 Jan 1985 were both Tuesdays, so NDAY mod 7 = 4 
            ! is a Saturday and NDAY mod 7 = 5 is a Sunday (bmy, 3/23/98)
            !=================================================================
            NDAY = ( GET_TAU() / 24d0 ) 
            IF ( MOD( NDAY, 7 ) .eq. 4 ) THEN
               JSCEN = 1
            ELSE IF ( MOD( NDAY, 7 ) .eq. 5 ) THEN
               JSCEN = 2
            ELSE
               JSCEN = 3
            ENDIF

            ! Fossil Fuel emissions (kg / Grid-Box / Time-Step)
            ! NN = tracer number corresponding to emission species N
            IF ( LANTHRO ) THEN 
               DO N = 1, NEMANTHRO
                  NN = IDEMS(N)
                  IF ( NN /= 0 ) THEN
                     CALL EMFOSSIL( I, J, N, NN, IREF, JREF, JSCEN,
     &                              Input_Opt )
                  ENDIF
               ENDDO
            ENDIF

!-----------------------------------------------------------------------------
! LIGHTNING EMISSIONS NOX [molecules/cm3/s]
!
            IF ( LLIGHTNOX .AND. I == 1 .AND. J == 1 ) THEN
               CALL EMLIGHTNING( State_Met, State_Chm )
            ENDIF
!-----------------------------------------------------------------------------
! AIRCRAFT emissions [kg/box]
!
            IF ( LAEIC .AND. I == 1 .AND. J == 1 ) THEN
               ! Choose between RCP inventory and default
               ! (cdh, 6/14/12)
               IF (LRCPAIR) THEN
                  CALL RCP_AIREMISS( State_Met )
               ELSE
                  CALL EMIT_AEIC( Input_Opt, State_Met, State_Chm )
               ENDIF
            ENDIF
!-----------------------------------------------------------------------------
! NOx AND CO from biofuel combustion [kg/box]  
!
            IF ( LWOODCO .AND. I == 1 .AND. J == 1 ) THEN
               CALL BIOFUEL_BURN
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF
!----------------------------------------------------------------------------
! BIOGENIC EMISSIONS OF VARIOUS QUANTITIES [Atoms C/box/time step]
!
            IF ( LBIOGENIC ) THEN

               ! Temperature
               TMMP  = XLTMMP( I, J, State_Met%TS, IJLOOP )
            
               ! Modified to choose MEGAN inventory indenpendantly 
               ! for ISOP and MONX/MBO. (ccc, 1/20/09)
               IF ( LMEGAN ) THEN

                  !------------------
                  ! MEGAN Isoprene
                  !------------------

                  ! Cosine of solar zenith angle
                  SC   = State_Met%SUNCOSmid(I,J)

                  ! Diffuse and direct PAR
                  PDR  = State_Met%PARDR(I,J)
                  PDF  = State_Met%PARDF(I,J)

                  ! Isoprene         
                  EMIS = GET_EMISOP_MEGAN(  I, J,     SC, TMMP, 
     &                                      PDR, PDF, XNUMOL_C )

               ENDIF
 
               !FP_ISOP (6/2009)
               EMIS = ISOP_SCALING * EMIS

               IF ( LMEGANMONO ) THEN

                  !------------------
                  ! MEGAN biogenics
                  !------------------

                  ! Methyl butenol
                  EMMB = GET_EMMBO_MEGAN(   I, J,     SC, TMMP,
     &                                      PDR, PDF, XNUMOL_C )

                  ! ------------------------------------------
                  ! Alpha Pinene emissions
                  APINE = GET_EMTERP_MEGAN( I , J , SC, TMMP ,  
     &                       PDR , PDF , XNUMOL_C , 'APINE' )
                  ! ------------------------------------------
                  ! Beta Pinene emissions
                  BPINE = GET_EMTERP_MEGAN( I , J , SC, TMMP ,  
     &                       PDR , PDF , XNUMOL_C , 'BPINE' )
                  ! ------------------------------------------
                  ! Limonene emissions
                  LIMON = GET_EMTERP_MEGAN( I , J , SC, TMMP , 
     &                       PDR , PDF , XNUMOL_C , 'LIMON' )
                  ! ------------------------------------------
                  ! Sabinene emissions
                  SABIN = GET_EMTERP_MEGAN( I , J , SC, TMMP ,  
     &                       PDR , PDF , XNUMOL_C , 'SABIN' )               
                  ! ------------------------------------------
                  ! Mycrene emissions
                  MYRCN = GET_EMTERP_MEGAN( I , J , SC, TMMP ,  
     &                       PDR , PDF , XNUMOL_C , 'MYRCN' )
                  ! ------------------------------------------
                  ! 3-Carene emissions
                  CAREN = GET_EMTERP_MEGAN( I , J , SC, TMMP ,  
     &                       PDR , PDF , XNUMOL_C , 'CAREN' )
                  ! ------------------------------------------
                  ! Ocimene emissions,
                  OCIMN = GET_EMTERP_MEGAN( I , J , SC, TMMP ,   
     &                       PDR , PDF , XNUMOL_C , 'OCIMN' )
                  ! ------------------------------------------

                  ! Total monoterpenes = sum of individual
                  EMMO = APINE + BPINE + LIMON + SABIN + 
     &                   MYRCN + CAREN + OCIMN

                  ! SOAupdate: Sesquiterpenes, OMTP (hotp 3/1/10) 
                  ! (not actually used here)
                  ! ------------------------------------------
                  ! a-Farnesene emissions
                  FARN = GET_EMTERP_MEGAN( I , J , SC, TMMP ,  
     &                      PDR , PDF , XNUMOL_C , 'FARNE' )
                  ! ------------------------------------------
                  ! b-Caryophyllene emissions
                  BCAR = GET_EMTERP_MEGAN( I , J , SC, TMMP ,  
     &                      PDR , PDF , XNUMOL_C , 'BCARE' )
                  ! ------------------------------------------
                  ! Other sesquiterpene emissions
                  OSQT = GET_EMTERP_MEGAN( I , J , SC, TMMP ,   
     &                      PDR , PDF , XNUMOL_C , 'OSQTE' )
                  ! ------------------------------------------
                  ! Other monoterpene emissions
                  OMTP = GET_EMTERP_MEGAN( I , J , SC, TMMP ,   
     &                      PDR , PDF , XNUMOL_C , 'OMTPE' )
                  ! ------------------------------------------

               ENDIF

               ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               ! Store BVOC emissions for diagnositics 48,49,51 (mpb,2009)
               !
               ! Species            | Order 
               ! ----------------------------
               ! Isoprene           = 1
               ! Total Monoterpenes = 2 
               ! MBO                = 3 
               ! Alpha-Pinene       = 4
               ! Beta-Pinene        = 5
               ! Limonene           = 6
               ! Sabinene           = 7
               ! Mycrene            = 8
               ! 3-Carene           = 9
               ! Ocimene            = 10
               ! a-Farnesene        = 11
               ! b-Caryophyllene    = 12
               ! Other sesquiterp.  = 13
               ! Other monoterpenes = 14

               EMISS_BVOC( I , J , 1 ) =  EMIS   / AREA_CM2 / DTSRCE     
               EMISS_BVOC( I , J , 2 ) =  EMMO   / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 3 ) =  EMMB   / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 4 ) =  APINE  / AREA_CM2 / DTSRCE     
               EMISS_BVOC( I , J , 5 ) =  BPINE  / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 6 ) =  LIMON  / AREA_CM2 / DTSRCE  
               EMISS_BVOC( I , J , 7 ) =  SABIN  / AREA_CM2 / DTSRCE     
               EMISS_BVOC( I , J , 8 ) =  MYRCN  / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 9 ) =  CAREN  / AREA_CM2 / DTSRCE  
               EMISS_BVOC( I , J , 10) =  OCIMN  / AREA_CM2 / DTSRCE
               ! SOAupdate: sesquiterpenes (hotp 3/1/10), other monoterps (3/5)
               EMISS_BVOC( I , J , 11) =  FARN   / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 12) =  BCAR   / AREA_CM2 / DTSRCE  
               EMISS_BVOC( I , J , 13) =  OSQT   / AREA_CM2 / DTSRCE 
               EMISS_BVOC( I , J , 14) =  OMTP   / AREA_CM2 / DTSRCE
               

!-----------------------------------------------------------------------------
! BIOGENIC ACETONE EMISSIONS
!
               IF ( IDEACET /= 0 ) THEN

                  ! evf, edits to use MEGAN biogenic acetone emiss (5/25/2011)
                  ! Compute biogenic acetone emissions [atoms C/box/s]
                  CALL EMISS_BIOACET( I,    J,    TMMP,  EMMO, 
     &                                SC, PDR, PDF, XNUMOL_C,
     &                                EMIS, EMMB, GRASS, BIO_ACET )

               
                  ! Also add ocean source of acetone [atoms C/box/s]
                  CALL OCEAN_SOURCE_ACET( I, J, BIO_ACET, State_Met )

                  ! Add biogenic acetone to anthro source [atoms C/box/s]
                  ! NOTE: Don't save into EMISRR for the tagged CO 
                  ! simulation (jaf, mak, bmy, 2/14/08)
                  IF ( NO_TAGCO ) THEN
                    EMISRR(I,J,IDEACET) = EMISRR(I,J,IDEACET) + BIO_ACET
                  ENDIF
               ENDIF
!-----------------------------------------------------------------------------

               !==============================================================
               ! save biogenic isoprene emission for later use
               ! EMISRR has units [atoms C/box/s] 
               !==============================================================
               IF ( IDTISOP /= 0 ) THEN

                  ! NOTE: Don't save into EMISRR for the tagged CO 
                  ! simulation (jaf, mak, bmy, 2/14/08)
                  IF ( NO_TAGCO ) THEN
                     EMISRR(I,J,IDEISOP) = EMISRR(I,J,IDEISOP) + 
     &                                     ( EMIS / DTSRCE )
                  ENDIF
               ENDIF

            !=================================================================
            ! save biogenic monoterpene emission for later use
            ! EMISRR has units [atoms C/box/s] (tmf, 4/10/06)
            !=================================================================
            IF ( IDTMONX /= 0 ) THEN
               EMISRR(I,J,IDEMONX) = EMISRR(I,J,IDEMONX) + 
     &                               ( EMMO / DTSRCE )
            ENDIF

!------------------------------------------------------------------------------
!
!******************************************************************************
!  Biogenic source of CO -- from oxidation of METHANOL and MONOTERPENES
!
!  CO from METHANOL oxidation -- scaled from ISOPRENE (bnd, 1/2/01)
!
!    We need to scale the Isoprene flux to get the CH3OH (methanol) flux.
!    Currently, the annual isoprene flux in GEOS-CHEM is ~ 397 Tg C.
!
!    Daniel Jacob recommends a flux of 100 Tg/yr CO from CH3OH oxidation  
!    based on Singh et al. 2000 [JGR 105, 3795-3805] who estimate a global 
!    methanol source of 122 Tg yr-1, of which most (75 Tg yr-1) is 
!    "primary biogenic".  He also recommends for now that the CO flux 
!    from CH3OH oxidation be scaled to monthly mean isoprene flux.
!
!    To get CO from METHANOL oxidation, we must therefore multiply
!    the ISOPRENE flux by the following scale factor:
!      ( 100 Tg CO from CH3OH Oxidation  / 397 Tg C from Isoprene Flux ) *
!      (  12 g C/mole                    / 28 g CO/mole                )
!
!  CO from MONOTERPENE oxidation (bnd, bmy, 1/2/01)
!
!    Assume the production of CO from monoterpenes is instantaneous even 
!    though the lifetime of intermediate species may be on the order of hours 
!    or days.  This assumption will likely cause CO from monoterpene oxidation
!    to be too high in the box in which the monoterpene is emitted.
!
!    The CO yield here is taken from:
!
!    Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
!      "The ultimate yield of CO from the tropospheric oxidation of terpenes 
!       (including both O3 and OH reactions) was estimated to be 20% on the 
!       carbon number basis."  They studied ALPHA- & BETA-pinene.
!
!    Vinckier et al. Fresenius Env. Bull., Vol. 7, p.361-368 (1998)
!      "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
!******************************************************************************
!
               !=====================================================
               ! CO from MONOTERPENE oxidation [molec CO/box/s] 
               !====================================================
 
               ! NOTE: Don't save into EMISRR for the tagged CO simulation.  
               ! Also for tagged CO we don't use monoterpenes  ??????
               ! (jaf, mak, bmy, 2/14/08)
               IF ( NO_TAGCO ) THEN
                  TMPVAL            = ( EMMO / DTSRCE ) * 0.2d0
                  EMISRR(I,J,IDECO) = EMISRR(I,J,IDECO) + TMPVAL

                  ! ND29: CO-source from monoterpenes [molec/cm2/s]
                  IF ( ND29 > 0 ) THEN
                     AD29(I,J,5) = AD29(I,J,5) + ( TMPVAL / AREA_CM2 )
                  ENDIF
               ENDIF
!
!******************************************************************************
!  Biogenic source of PRPE -- scaled to ISOPRENE
!
!  Also, add biogenic emissions of alkenes. We do this by scaling to
!  isoprene emissions (probably OK for summertime conditions). The
!  scaling factor is based on work by Allen Goldstein. His values
!  indicate emission ratios of ethene:propene:butene=4:2:1 (on a
!  per molecule basis), with total emissions approx. equal to
!  10% of isoprene emissions (again, on molecule basis).
!  BIOSCAL is in units of atoms C (alkenes) / atoms C (isoprene)
!******************************************************************************
! Change this factor to exclude ethene (bey, ljm)
!    (10 molec alkenes / 100 molec isop) * (1 molec isop / 5 atoms C isop)
!    *(3 molec butene + propene / 7 molec total alkenes)
!    *(3.3333 atoms C but+prop mix/ 1 molec but+prop mix)
!    = 0.0286 atoms C butene+propene / atom C isop
! Note that 3.3333 atoms C/molecule is the weighted average for this mix.
!******************************************************************************
!
               BIOSCAL = 0.0286d0 ! new factor, (ljm, bey, 9/28/98)

               IF ( IDEPRPE /= 0 ) THEN

                  ! NOTE: Don't save into EMISRR for the tagged 
                  ! CO simulation. (jaf, mak, bmy, 2/14/08)
                  IF ( NO_TAGCO ) THEN 
                     EMISRR(I,J,IDEPRPE) = EMISRR(I,J,IDEPRPE) +
     &                                     ( EMIS / DTSRCE ) * BIOSCAL
                  ENDIF
               ENDIF

!=======================================================================
! Add biogenic emission of ethene (C2H4) --> scaled to isoprene
!
! Scale factor BIOSCALEC2H4 =
!   ( 10 molec alkenes / 100 molec isop ) * ( 1 molec isop / 5 atoms C )
!   * ( 4 molec ethene / 7 molec alkenes )
!   * ( 2 atoms C / 1 molec ethene )
!   = 0.022857d0  [atoms C / atoms C isop]
! (tmf, 1/13/06)

            BIOSCALEC2H4 = 0.022857d0

            IF ( IDEC2H4 /= 0 ) THEN
               EMISRR(I,J,IDEC2H4) = EMISRR(I,J,IDEC2H4) +
     &            ( EMIS / DTSRCE ) * BIOSCALEC2H4
            ENDIF
!=======================================================================


!
!-----------------------------------------------------------------------------
! BIOGENIC CHBr3 and CH2Br2 EMISSIONS, (jpp, 8/24/07)
!  Both in [molecules/box/s] as expected by setemiss.f
!
!

            IF ( IDECHBr3 /= 0 ) THEN
               EMISRR(I,J,IDECHBr3) = EMIS_CHBr3( I, J, State_Met )

               ! Scale emissions
               EMISRR(I,J,IDECHBr3) = EMISRR(I,J,IDECHBr3) * Br_SCALING

            ENDIF

            IF ( IDECH2Br2 /= 0 ) THEN
               EMISRR(I,J,IDECH2Br2) = EMIS_CH2Br2(I,J)

               ! Scale emissions
               EMISRR(I,J,IDECH2Br2) = EMISRR(I,J,IDECH2Br2) * 
     &                                 Br_SCALING

            ENDIF

!
!******************************************************************************
!  ND46 diagnostic: Biogenic emissions 
!
!     AD46(:,:,1)  = Total biogenic ISOP     emissions [atoms C/cm2/s]
!     AD46(:,:,2)  = Total biogenic ACET     emissions [atoms C/cm2/s]
!     AD46(:,:,3)  = Total biogenic PRPE     emissions [atoms C/cm2/s]
!     AD46(:,:,4)  = Total biogenic MONOT    emissions [atoms C/cm2/s]
!     AD46(:,:,5)  = Total biogenic MBO      emissions [atoms C/cm2/s]
!     AD46(:,:,6)  = Total biogenic C2H4     emissions [atoms C/cm2/s]
!     AD46(:,:,7)  = Total biogenic a-pinene emissions [atoms C/cm2/s]
!     AD46(:,:,8)  = Total biogenic b-pinene emissions [atoms C/cm2/s]
!     AD46(:,:,9)  = Total biogenic limonene emissions [atoms C/cm2/s]
!     AD46(:,:,10) = Total biogenic sabinene emissions [atoms C/cm2/s]
!     AD46(:,:,11) = Total biogenic mycrene  emissions [atoms C/cm2/s]
!     AD46(:,:,12) = Total biogenic 3-carene emissions [atoms C/cm2/s]
!     AD46(:,:,13) = Total biogenic ocimene  emissions [atoms C/cm2/s]
!     AD46(:,:,14) = Total biogenic FARN     emissions [atoms C/cm2/s]
!     AD46(:,:,15) = Total biogenic BCAR     emissions [atoms C/cm2/s]
!     AD46(:,:,16) = Total biogenic OSQT     emissions [atoms C/cm2/s]
!     AD46(:,:,17) = Total biogenic OMTP     emissions [atoms C/cm2/s]
!     AD46(:,:,18) = Total biogenic CHBr3    emissions [kg/m2/s]
!     AD46(:,:,19) = Total biogenic CH2Br2   emissions [kg/m2/s]
!     AD46(:,:,20) = Total biogenic SSBr2    emissions [kg/m2/s]
!
!  NOTES: 
!  (1 ) Now make ACET tracer #2 and PRPE tracer #3 (bmy, 9/13/01)
!  (2 ) Now archive ND46 as [atoms C/cm2/s] here (bmy, 9/13/01)
!  (3 ) Added MBO emission diagnostics [atoms C/cm2/s] (bmy, tmf, 10/20/05)
!  22 Oct 2012 - R. Yantosca - Prevent out-of-bounds errors in ND46 if bromine
!                              emissions are turned off (e.g. for tagged CO)
!  13 Aug 2013 - M. Sulprizio- Add sesquiterpenes for updated SOA (H.Pye)
!******************************************************************************
!
               IF ( ND46 > 0 ) THEN

                  ! ISOP emissions [atoms C/cm2/s] -- tracer #1
                  AD46(I,J,1) = AD46(I,J,1) + ( EMIS / AREA_CM2 /DTSRCE)

                  ! ACET emissions [atoms C/cm2/s] -- tracer #2
                  AD46(I,J,2) = AD46(I,J,2) + ( BIO_ACET / AREA_CM2 )

                  ! PRPE emissions [atoms C/cm2/s] -- tracer #3
                  AD46(I,J,3) = AD46(I,J,3) + 
     &                          ( EMIS * BIOSCAL / AREA_CM2 / DTSRCE )

                  ! Monoterpene emissions [atoms C/cm2/s] -- tracer #4
                  AD46(I,J,4) = AD46(I,J,4) + ( EMMO / AREA_CM2 /DTSRCE) 

                  ! MBO emissions [atoms C/cm2/s] -- tracer #5
                  AD46(I,J,5) = AD46(I,J,5) + ( EMMB / AREA_CM2 /DTSRCE) 

                  ! C2H4 emissions [atoms C/cm2/s] -- tracer #6 (tmf, 1/13/06)
                  AD46(I,J,6) = AD46(I,J,6) + 
     &            ( EMIS * BIOSCALEC2H4 / AREA_CM2 / DTSRCE )

                  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++
                  !              MEGAN v2.1 - (mpb,2009)               !
                  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++

                  ! Aplha Pinene emissions [atoms C/cm2/s] -- tracer #7
                  AD46(I,J,7)  = AD46(I,J,7)  +
     &                         ( APINE        / AREA_CM2 / DTSRCE ) 
 
                  ! Beta Pinene emissions [atoms C/cm2/s] -- tracer #8
                  AD46(I,J,8)  = AD46(I,J,8)  +
     &                         ( BPINE        / AREA_CM2 / DTSRCE ) 

                  ! Limonene emissions [atoms C/cm2/s] -- tracer #9
                  AD46(I,J,9)  = AD46(I,J,9)  +
     &                         ( LIMON        / AREA_CM2 / DTSRCE ) 
 
                  ! Sabinene emissions [atoms C/cm2/s] -- tracer #10
                  AD46(I,J,10) = AD46(I,J,10) + 
     &                         ( SABIN        / AREA_CM2 / DTSRCE) 

                  ! Mycrene emissions [atoms C/cm2/s] -- tracer #11
                  AD46(I,J,11) = AD46(I,J,11) +
     &                         ( MYRCN        / AREA_CM2 / DTSRCE) 

                  ! 3-Carene emissions [atoms C/cm2/s] -- tracer #12
                  AD46(I,J,12) = AD46(I,J,12) +
     &                         ( CAREN        / AREA_CM2 / DTSRCE) 

                  ! Ocimene emissions [atoms C/cm2/s] -- tracer #13
                  AD46(I,J,13) = AD46(I,J,13) +
     &                         ( OCIMN        / AREA_CM2 / DTSRCE) 

                  ! SOAupdate (hotp 7/29/10)
                  ! Sesquiterpenes: 14-16 (hotp 3/1/10)
                  ! Farnesene emissions [atoms C/cm2/s] -- tracer #14
                  AD46(I,J,14) = AD46(I,J,14) + 
     &                         ( FARN         / AREA_CM2 / DTSRCE)

                  ! b-caryophyllene emissions [atoms C/cm2/s] -- tracer #15
                  AD46(I,J,15) = AD46(I,J,15) + 
     &                         ( BCAR         / AREA_CM2 / DTSRCE)

                  ! Other SQT emissions [atoms C/cm2/s] -- tracer #16
                  AD46(I,J,16) = AD46(I,J,16) + 
     &                         ( OSQT         / AREA_CM2 / DTSRCE)
 
                  ! Other MTP emissions [atoms C/cm2/s] -- tracer #17
                  AD46(I,J,17) = AD46(I,J,17) + 
     &                         ( OMTP         / AREA_CM2 / DTSRCE)
                  ! end hotp

                  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++

                  ! CHBr3 emissions [kg/m2/s] -- tracer #18
                  IF ( IDECHBr3 > 0 ) THEN
                     AD46(I,J,18) = AD46(I,J,18) 
     &                            + ( EMISRR(I,J,IDECHBr3)  / AREA_M2 ) 
     &                            * ( MWT_CHBr3             / AVG     )
                  ENDIF

                  ! CH2Br2 emissions [kg/m2/s] -- tracer #19
                  IF ( IDECH2Br2 > 0 ) THEN
                     AD46(I,J,19) = AD46(I,J,19) 
     &                            + ( EMISRR(I,J,IDECH2Br2) / AREA_M2 ) 
     &                            * ( MWT_CH2Br2            / AVG     )
                  ENDIF

                  ! Note: IDEBr2 is added to AD46(I,J,20) in ssa_bromine_mod.F
               
               ENDIF ! ND46
            ENDIF    ! LBIOGENIC
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE EMISSDR                                                

