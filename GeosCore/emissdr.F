! $Id: emissdr.f,v 1.6 2010/03/15 19:33:24 ccarouge Exp $
#if defined( DEVEL )
      SUBROUTINE EMISSDR( State_Met, EMISS_TEND, BIOEMISS_TEND )
#else
      SUBROUTINE EMISSDR( State_Met )
#endif
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
!******************************************************************************
!
      ! References to F90 modules
      USE ACETONE_MOD,        ONLY : EMISS_BIOACET, OCEAN_SOURCE_ACET
      USE AIRCRAFT_NOX_MOD,   ONLY : AIREMISS
      USE BIOFUEL_MOD,        ONLY : BIOFUEL_BURN
      USE BROMOCARB_MOD,      ONLY : EMIS_CHBr3, EMIS_CH2Br2
      USE BROMOCARB_MOD,      ONLY : Br_SCALING
      USE CANOPY_NOX_MOD,     ONLY : GET_CANOPY_NOX
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD
      USE CMN_O3_MOD
      USE CMN_DIAG_MOD
      USE CMN_NOX_MOD
      USE DIAG_MOD,           ONLY : AD29,          AD46
      USE DIAG_MOD,           ONLY : AD32_fe,       AD32_so
      USE EMISSIONS_MOD,      ONLY : ISOP_SCALING
      USE GET_NDEP_MOD,       ONLY : GET_DEP_N 
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE GRID_MOD,           ONLY : GET_YMID
      USE GRID_MOD,           ONLY : GET_XMID
      USE GRID_MOD,           ONLY : GET_XOFFSET,   GET_YOFFSET
      USE LIGHTNING_NOX_MOD,  ONLY : EMLIGHTNING
      USE LOGICAL_MOD,        ONLY : LANTHRO,       LLIGHTNOX, LSOILNOX
      USE LOGICAL_MOD,        ONLY : LFERTILIZERNOX
      USE LOGICAL_MOD,        ONLY : LAIRNOX,       LBIONOX,   LWOODCO   
      USE LOGICAL_MOD,        ONLY : LMEGAN, LMEGANMONO, LBIOGENIC
      USE MEGAN_MOD,          ONLY : GET_MEGAN_EMISSIONS 
      USE MEGANUT_MOD,        ONLY : XLTMMP
      USE MODIS_LAI_MOD
      USE SOIL_NOX_MOD,       ONLY : SOIL_NOX_EMISSION   
      USE TIME_MOD,           ONLY : GET_MONTH,     GET_TAU
      USE TIME_MOD,           ONLY : GET_TS_EMIS,   GET_LOCALTIME
      USE TIME_MOD,           ONLY : GET_DAY_OF_YEAR
      USE TRACER_MOD,         ONLY : ITS_A_TAGCO_SIM
      USE TRACERID_MOD,       ONLY : IDEACET,       IDTISOP,   IDEISOP   
      USE TRACERID_MOD,       ONLY : IDECO,         IDEPRPE,   NEMANTHRO 
      USE TRACERID_MOD,       ONLY : IDEMONX, IDEC2H4
      USE TRACERID_MOD,       ONLY : IDTMONX, IDTC2H4
      USE TRACERID_MOD,       ONLY : IDECHBr3, IDECH2Br2

#if defined( DEVEL )
      USE TRACERID_MOD,       ONLY : IDTNOX
      USE TRACER_MOD,         ONLY : N_TRACERS
#endif
 
      USE CMN_SIZE_MOD             ! Size parameters
      USE COMODE_LOOP_MOD          ! IVERT?
      USE CMN_O3_MOD               ! Emissions arrays
      USE CMN_DIAG_MOD             ! Diagnostic arrays and switches
      USE CMN_NOX_MOD              ! GEMISNOX2
      USE COMMSOIL_MOD             ! Soil wetness variables

      IMPLICIT NONE
#     include "define.h"

      !Arguments
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object

#if defined( DEVEL )
      REAL*8,         INTENT(INOUT) :: 
     &     EMISS_TEND(IIPAR,JJPAR,LLPAR,N_TRACERS+1),
     &     BIOEMISS_TEND(IIPAR,JJPAR,LLPAR,N_TRACERS+1)
#endif

      ! Local variables
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE. 
      LOGICAL                :: NO_TAGCO
      INTEGER                :: I, J, L, N, IJLOOP 
      INTEGER                :: I0, J0, IOFF, JOFF, IREF, JREF
      INTEGER                :: NDAY,JSCEN,NN
      INTEGER, SAVE          :: LASTMONTH = -99
      REAL*8                 :: DTSRCE,  XLOCTM,   TMMP
      REAL*8                 :: EMIS_ISOP, EMIS_PRPE
      REAL*8                 :: EMIS_MONO, EMIS_MBOX, EMIS_C2H4
      REAL*8                 :: EMIS_APIN, EMIS_BPIN, EMIS_LIMO
      REAL*8                 :: EMIS_SABI, EMIS_MYRC, EMIS_CARE
      REAL*8                 :: EMIS_OCIM, EMIS_FAXX, EMIS_AAXX
      REAL*8                 :: EMIS_ALD2, EMIS_OMON, EMIS_MOHX
      REAL*8                 :: EMIS_ETOH
      REAL*8                 :: AREA_CM2
      REAL*8                 :: TMPVAL,    GRASS, BIO_ACET
      REAL*8                 :: CONVERT(NVEGTYPE), GMONOT(NVEGTYPE)
      REAL*8                 :: SC, PDF, PDR
      INTEGER                :: DOY, MONTH, RC, M
      REAL*8                 :: FERTDIAG
      REAL*8                 :: TS_EMIS
      REAL*8                 :: TOTALPULSE
      REAL*8                 :: SOILFRT
      REAL*8                 :: DEP_FERT, TEMP_FERT


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
      ! Add in gamma activity factors for isoprene ONLY (mpb,2009)
      REAL*8                 :: GAMMA_LAI
      REAL*8                 :: GAMMA_LEAF_AGE
      REAL*8                 :: GAMMA_P
      REAL*8                 :: GAMMA_T
      REAL*8                 :: GAMMA_SM

      ! External functions
      REAL*8, EXTERNAL       :: BOXVL
!
!******************************************************************************
!  EMISSDR begins here!
!
!  DTSRCE = emission timestep in seconds
!
!  Call subroutines to set up ISOP and monoterpene emission (first time only!)
!******************************************************************************
!
      ! This is not a tagged CO simulation
      NO_TAGCO = ( .not. ITS_A_TAGCO_SIM() )

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
            ! Zero other biogenic species as well
            EMIS_ISOP = 0d0
            EMIS_PRPE = 0d0
            EMIS_MONO = 0d0
            EMIS_MBOX = 0d0
            EMIS_C2H4 = 0d0
            EMIS_APIN = 0d0
            EMIS_BPIN = 0d0
            EMIS_LIMO = 0d0
            EMIS_SABI = 0d0
            EMIS_MYRC = 0d0
            EMIS_CARE = 0d0
            EMIS_OCIM = 0d0
            EMIS_FAXX = 0d0
            EMIS_AAXX = 0d0
            EMIS_ALD2 = 0d0
            EMIS_OMON = 0d0
            EMIS_MOHX = 0d0
            EMIS_ETOH = 0d0

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
                     CALL EMFOSSIL( I, J, N, NN, IREF, JREF, JSCEN )
                  ENDIF
               ENDDO
            ENDIF

!-----------------------------------------------------------------------------
! LIGHTNING EMISSIONS NOX [molecules/cm3/s]
!
#if defined( DEVEL )
            IF ( LLIGHTNOX .AND. I == 1 .AND. J == 1 ) 
     &           CALL EMLIGHTNING( State_Met, EMISS_TEND(:,:,:,IDTNOX) )
#else
            IF ( LLIGHTNOX ) CALL EMLIGHTNING( I, J, State_Met )
#endif
!-----------------------------------------------------------------------------
! AIRCRAFT emissions NOx [molecules/cm3/s]
!
#if defined( DEVEL )
            IF ( LAIRNOX .AND. I == 1 .AND. J == 1 ) 
     &           CALL AIREMISS( State_Met, EMISS_TEND(:,:,:,IDTNOX) )
#else
            IF ( LAIRNOX .AND. I == 1 .AND. J == 1 ) 
     &           CALL AIREMISS( State_Met )
#endif
!-----------------------------------------------------------------------------
! NOx AND CO from biofuel combustion [kg/box]  
!
#if defined( DEVEL )
            IF ( LWOODCO .AND. I == 1 .AND. J == 1 ) 
     &           CALL BIOFUEL_BURN( State_Met, BIOEMISS_TEND(:,:,1,:) )
#else
            IF ( LWOODCO .AND. I == 1 .AND. J == 1 ) 
     &           CALL BIOFUEL_BURN( State_Met )
#endif
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

                  ! Isoprene         
                  ! Updated to new MEGAN routine (dbm, 12/2012)
                  EMIS_ISOP = GET_MEGAN_EMISSIONS( I,         J, 
     &                                             State_Met, 'ISOP' )

               ENDIF
 
               !FP_ISOP (6/2009)
               EMIS_ISOP = ISOP_SCALING * EMIS_ISOP

               IF ( LMEGANMONO ) THEN

                  !------------------
                  ! MEGAN biogenics
                  !------------------
                  ! Updated the calls below for new MEGAN routine (dbm, 12/2012)

                  ! Methyl butenol
                  EMIS_MBOX = GET_MEGAN_EMISSIONS( I,         J, 
     &                                             State_Met, 'MBOX' )

                  ! ------------------------------------------
                  ! Aplha Pinene emissions
                  EMIS_APIN = GET_MEGAN_EMISSIONS( I,         J, 
     &                                             State_Met, 'APIN' )
                  ! ------------------------------------------
                  ! Beta Pinene emissions
                  EMIS_BPIN = GET_MEGAN_EMISSIONS( I,         J, 
     &                                             State_Met, 'BPIN' )
                  ! ------------------------------------------
                  ! Limonene emissions 
                  EMIS_LIMO = GET_MEGAN_EMISSIONS( I,         J,
     &                                             State_Met, 'LIMO' )
                  ! ------------------------------------------
                  ! Sabinene emissions
                  EMIS_SABI = GET_MEGAN_EMISSIONS( I,         J, 
     &                                             State_Met, 'SABI' )         
                  ! ------------------------------------------
                  ! Mycrene emissions
                  EMIS_MYRC = GET_MEGAN_EMISSIONS( I,         J, 
     &                                             State_Met, 'MYRC' )
                  ! ------------------------------------------
                  ! 3-Carene emissions 
                  EMIS_CARE = GET_MEGAN_EMISSIONS( I,         J 
     &                                             State_Met, 'CARE' )
                  ! ------------------------------------------
                  ! Ocimene emissions
                  EMIS_OCIM = GET_MEGAN_EMISSIONS( I,         J, 
     &                                             State_Met, 'OCIM' )
                  ! ------------------------------------------
                  ! Other monoterpenes
                  ! (added 12/2012; dbm)
                  EMIS_OMON = GET_MEGAN_EMISSIONS( I,         J, 
     &                                             State_Met, 'OMON' )
                  ! ------------------------------------------

                  ! Total monoterpenes = sum of individual
                  ! dbm, now add other lumped monoterpenes (11/2012)
                  EMIS_MONO = EMIS_APIN + EMIS_BPIN + EMIS_LIMO + 
     &                        EMIS_SABI + EMIS_MYRC + EMIS_CARE + 
     &                        EMIS_OCIM + EMIS_OMON

               ENDIF

               !------------------
               ! MEGAN ALD2
               !------------------
               IF ( IDEALD2 /= 0 ) THEN

                  EMIS_ALD2 = GET_MEGAN_EMISSIONS( I,          J,  
     &                                             State_Met, 'ALD2' )

                  ! NOTE: Don't save into EMISRR for the tagged CO 
                  ! simulation (jaf, mak, bmy, 2/14/08)
                  IF ( NO_TAGCO ) THEN
                     EMISRR(I,J,IDEALD2) = EMISRR(I,J,IDEALD2) + 
     &                                     ( EMIS_ALD2 / DTSRCE )
                  ENDIF

               ENDIF

               !------------------
               ! MEGAN Acetone
               !------------------
               ! Computed within acetone_mod.f

               ! ----------------------------
               ! Other biogenics. 
               ! Calls included here for future incorporation or
               ! specialized simulations. (dbm, 12/2012)
               ! ----------------------------

               !------------------
               ! MEGAN Methanol
               !------------------
               ! Uncomment this to compute biogenic methanol emissions
!              EMIS_MOHX = GET_MEGAN_EMISSIONS( I,         J,  
!     &                                         State_Met, 'MOHX' )

               !------------------
               ! MEGAN Ethanol
               !------------------
               ! Uncomment this to compute biogenic ethanol emissions
!               EMIS_ETOH = GET_MEGAN_EMISSIONS( I,         J,
!     &                                          State_Met, 'ETOH' )

               !------------------
               ! MEGAN Formic acid
               !------------------
               ! Uncomment this to compute biogenic formic acid emissions
!               EMIS_FAXX = GET_MEGAN_EMISSIONS( I,         J,
!     &                                          State_Met, 'FAXX' ) 

               !------------------
               ! MEGAN Acetic acid
               !------------------
               ! Uncomment this to compute biogenic acetic acid emissions
!               EMIS_AAXX = GET_MEGAN_EMISSIONS( I,         J,
!     &                                          State_Met, 'AAXX' )

!-----------------------------------------------------------------------------
! BIOGENIC ACETONE EMISSIONS
!
               IF ( IDEACET /= 0 ) THEN

!                  ! evf, edits to use MEGAN biogenic acetone emiss (5/25/2011)
!                  ! Compute biogenic acetone emissions [atoms C/box/s]
!                  CALL EMISS_BIOACET( I,    J,    TMMP,  EMMO, 
!     &                                SC, PDR, PDF, XNUMOL_C,
!     &                                EMIS, EMMB, GRASS, BIO_ACET )

                  ! evf, edits to use MEGAN biogenic acetone emiss (5/25/2011)
                  ! Compute biogenic acetone emissions [atoms C/box/s]
                  CALL EMISS_BIOACET( I,    J,    TMMP,  EMIS_MONO, 
     &                                SC, PDR, PDF, XNUMOL_C,
     &                                EMIS_ISOP, EMIS_MBOX, 
     &                                GRASS, BIO_ACET )
               
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
     &                                     ( EMIS_ISOP / DTSRCE )
                  ENDIF
               ENDIF

               !=================================================================
               ! save biogenic monoterpene emission for later use
               ! EMISRR has units [atoms C/box/s] (tmf, 4/10/06)
               !=================================================================
               IF ( IDTMONX /= 0 ) THEN
                  EMISRR(I,J,IDEMONX) = EMISRR(I,J,IDEMONX) + 
     &                 ( EMIS_MONO / DTSRCE )
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
                  TMPVAL            = ( EMIS_MONO / DTSRCE ) * 0.2d0
                  EMISRR(I,J,IDECO) = EMISRR(I,J,IDECO) + TMPVAL

                  ! ND29: CO-source from monoterpenes [molec/cm2/s]
                  IF ( ND29 > 0 ) THEN
                     AD29(I,J,5) = AD29(I,J,5) + ( TMPVAL / AREA_CM2 )
                  ENDIF
               ENDIF
!
!******************************************************************************
!  Biogenic source of PRPE

!  Now uses MEGAN2.1 (dbm, 12/2012)
!******************************************************************************
!
               IF ( IDEPRPE /= 0 ) THEN
                  
                  EMIS_PRPE = GET_MEGAN_EMISSIONS(   I, J,     SC, TMMP,
     &                 PDR, PDF, 'PRPE' )

                  ! NOTE: Don't save into EMISRR for the tagged 
                  ! CO simulation. (jaf, mak, bmy, 2/14/08)
                  IF ( NO_TAGCO ) THEN 
                     EMISRR(I,J,IDEPRPE) = EMISRR(I,J,IDEPRPE) +
     &                    ( EMIS_PRPE / DTSRCE )
                  ENDIF
               ENDIF

!=======================================================================
! Add biogenic emission of ethene (C2H4)
!
! Now uses MEGAN2.1 (dbm, 12/2012

               IF ( IDEC2H4 /= 0 ) THEN

                  EMIS_C2H4 = GET_MEGAN_EMISSIONS(   I, J,     SC, TMMP,
     &                 PDR, PDF, 'C2H4' )

                  EMISRR(I,J,IDEC2H4) = EMISRR(I,J,IDEC2H4) +
     &                 ( EMIS_C2H4 / DTSRCE ) 
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
                  EMISRR(I,J,IDECHBr3) = 
     &                 EMISRR(I,J,IDECHBr3) * Br_SCALING

               ENDIF

               IF ( IDECH2Br2 /= 0 ) THEN
                  EMISRR(I,J,IDECH2Br2) = EMIS_CH2Br2(I,J)

                  ! Scale emissions
                  EMISRR(I,J,IDECH2Br2) = 
     &                 EMISRR(I,J,IDECH2Br2) * Br_SCALING

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
!     AD46(:,:,14) = Total biogenic HCOOH    emissions [molec/cm2/s]
!     AD46(:,:,15) = Total biogenic ACTA     emissions [molec/cm2/s]
!     AD46(:,:,16) = Total biogenic ALD2     emissions [atoms C/cm2/s]
!     AD46(:,:,17) = Total other monoterpene emissions [atoms C/cm2/s]
!     AD46(:,:,18) = Total biogenic methanol emissions [molec/cm2/s]
!     AD46(:,:,19) = Total biogenic ethanol  emissions [atoms C/cm2/s]
!     AD46(:,:,20) = Total biogenic CHBr3    emissions [kg/m2/s]
!     AD46(:,:,21) = Total biogenic CH2Br2   emissions [kg/m2/s]
!     AD46(:,:,22) = Total biogenic SSBr2    emissions [kg/m2/s]
!
!  NOTES: 
!  (1 ) Now make ACET tracer #2 and PRPE tracer #3 (bmy, 9/13/01)
!  (2 ) Now archive ND46 as [atoms C/cm2/s] here (bmy, 9/13/01)
!  (3 ) Added MBO emission diagnostics [atoms C/cm2/s] (bmy, tmf, 10/20/05)
!  22 Oct 2012 - R. Yantosca - Prevent out-of-bounds errors in ND46 if bromine
!                              emissions are turned off (e.g. for tagged CO)
!  01 Dec 2012 - D. Millet   - Updated tracers (dbm, 12/2012)
!******************************************************************************
!
               IF ( ND46 > 0 ) THEN

                  ! ISOP emissions [atoms C/cm2/s] -- tracer #1
                  AD46(I,J,1) = AD46(I,J,1) + 
     &                 ( EMIS_ISOP / AREA_CM2 / DTSRCE )

                  ! ACET emissions [atoms C/cm2/s] -- tracer #2
                  AD46(I,J,2) = AD46(I,J,2) + ( BIO_ACET / AREA_CM2 )

                  ! PRPE emissions [atoms C/cm2/s] -- tracer #3
                  AD46(I,J,3) = AD46(I,J,3) + 
     &                 ( EMIS_PRPE / AREA_CM2 / DTSRCE )

                  ! Total monoterpene emissions [atoms C/cm2/s] -- tracer #4
                  AD46(I,J,4) = AD46(I,J,4) + 
     &                 ( EMIS_MONO / AREA_CM2 / DTSRCE ) 

                  ! MBO emissions [atoms C/cm2/s] -- tracer #5
                  AD46(I,J,5) = AD46(I,J,5) + 
     &                 ( EMIS_MBOX / AREA_CM2 / DTSRCE ) 

                  ! C2H4 emissions [atoms C/cm2/s] -- tracer #6 (tmf, 1/13/06)
                  AD46(I,J,6) = AD46(I,J,6) + 
     &                 ( EMIS_C2H4 / AREA_CM2 / DTSRCE )

                  ! Aplha Pinene emissions [atoms C/cm2/s] -- tracer #7
                  AD46(I,J,7) = AD46(I,J,7) + 
     &                 ( EMIS_APIN / AREA_CM2 / DTSRCE ) 

                  ! Beta Pinene emissions [atoms C/cm2/s] -- tracer #8
                  AD46(I,J,8) = AD46(I,J,8) + 
     &                 ( EMIS_BPIN / AREA_CM2 / DTSRCE ) 

                  ! Limonene emissions [atoms C/cm2/s] -- tracer #9
                  AD46(I,J,9) = AD46(I,J,9) + 
     &                 ( EMIS_LIMO / AREA_CM2 / DTSRCE ) 
 
                  ! Sabinene emissions [atoms C/cm2/s] -- tracer #10
                  AD46(I,J,10) = AD46(I,J,10) + 
     &                 ( EMIS_SABI / AREA_CM2 / DTSRCE ) 

                  ! Mycrene emissions [atoms C/cm2/s] -- tracer #11
                  AD46(I,J,11) = AD46(I,J,11) + 
     &                 ( EMIS_MYRC / AREA_CM2 / DTSRCE ) 

                 ! 3-Carene emissions [atoms C/cm2/s] -- tracer #12
                  AD46(I,J,12) = AD46(I,J,12) + 
     &                 ( EMIS_CARE / AREA_CM2 / DTSRCE ) 

                  ! Ocimene emissions [atoms C/cm2/s] -- tracer #13
                  AD46(I,J,13) = AD46(I,J,13) + 
     &                 ( EMIS_OCIM / AREA_CM2 / DTSRCE ) 

                  ! HCOOH [atoms /cm2/s] -- tracer #14
                  AD46(I,J,14) = AD46(I,J,14) + 
     &                 ( EMIS_FAXX / AREA_CM2 / DTSRCE ) 

                  ! ACTA [atoms /cm2/s] -- tracer #15
                  AD46(I,J,15) = AD46(I,J,15) + 
     &                 ( EMIS_AAXX / AREA_CM2 / DTSRCE ) 

                  ! ALD2 [atoms C/cm2/s] -- tracer #16
                  AD46(I,J,16) = AD46(I,J,16) + 
     &                 ( EMIS_ALD2 / AREA_CM2 / DTSRCE ) 

                  ! Other monoterpenes [atoms C/cm2/s] -- tracer #17
                  AD46(I,J,17) = AD46(I,J,17) + 
     &                 ( EMIS_OMON / AREA_CM2 / DTSRCE ) 

                  ! Methanol [molec/cm2/s] -- tracer #18
                  AD46(I,J,18) = AD46(I,J,18) + 
     &                 ( EMIS_MOHX / AREA_CM2 / DTSRCE ) 

                  ! Ethanol [atomsC/cm2/s] -- tracer #19
                  AD46(I,J,19) = AD46(I,J,19) + 
     &                 ( EMIS_ETOH / AREA_CM2 / DTSRCE ) 

                  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++
                  
                  ! CHBr3 emissions [kg/m2/s] -- tracer #20
                  AD46(I,J,20) = AD46(I,J,20) + 
     &                 ( EMISRR(I,J,IDECHBr3) / AREA_M2) * 
     &                 ( MWT_CHBr3 / AVG )

                  ! CH2Br2 emissions [kg/m2/s] -- tracer #21
                  AD46(I,J,21) = AD46(I,J,21) + 
     &                 ( EMISRR(I,J,IDECH2Br2) / AREA_M2) * 
     &                 ( MWT_CH2Br2 / AVG )

                  ! Note: IDEBr2 is added to AD46(I,J,16) in ssa_bromine_mod.F
               
               ENDIF  

               ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               ! Store instantaneous BVOC emissions for diagnostics 48,49,51 (mpb,2009)
               ! Added compounds and changed order to be consistent with
               ! AD46 (dbm, 12/2012)
               !
               ! Species            | Order 
               ! ----------------------------
               ! Isoprene           = 1
               ! Acetone            = 2
               ! Propene            = 3
               ! Total Monoterpenes = 4 
               ! MBO                = 5 
               ! Ethene             = 6
               ! Alpha-Pinene       = 7
               ! Beta-Pinene        = 8
               ! Limonene           = 9
               ! Sabinene           = 10
               ! Mycrene            = 11
               ! 3-Carene           = 12
               ! Ocimene            = 13
               ! Formic acid        = 14
               ! Acetic acid        = 15
               ! Acetaldehyde       = 16           
               ! Other monoterpenes = 17
               ! Methanol           = 18
               ! Ethanol            = 19
               ! UNITS: atomsC/cm2/s or molec/cm2/s
               ! BIO_ACET in atomsC/box/s. All other EMIS_**** variables
               ! are in atomsC/box/step or molec/box/step
               ! modify to have same order as AD46 (dbm, 11/2012)
               EMISS_BVOC( I , J , 1 ) =  EMIS_ISOP  / AREA_CM2 / DTSRCE   
               ! Don't need to divide BIO_ACET by DTSRCE
               EMISS_BVOC( I , J , 2 ) =  BIO_ACET / AREA_CM2
               EMISS_BVOC( I , J , 3 ) =  EMIS_PRPE  / AREA_CM2 / DTSRCE          
               EMISS_BVOC( I , J , 4 ) =  EMIS_MONO  / AREA_CM2 / DTSRCE                 
               EMISS_BVOC( I , J , 5 ) =  EMIS_MBOX  / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 6 ) =  EMIS_C2H4  / AREA_CM2 / DTSRCE     
               EMISS_BVOC( I , J , 7 ) =  EMIS_APIN  / AREA_CM2 / DTSRCE     
               EMISS_BVOC( I , J , 8 ) =  EMIS_BPIN  / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 9 ) =  EMIS_LIMO  / AREA_CM2 / DTSRCE  
               EMISS_BVOC( I , J , 10) =  EMIS_SABI  / AREA_CM2 / DTSRCE     
               EMISS_BVOC( I , J , 11) =  EMIS_MYRC  / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 12) =  EMIS_CARE  / AREA_CM2 / DTSRCE  
               EMISS_BVOC( I , J , 13) =  EMIS_OCIM  / AREA_CM2 / DTSRCE 
               EMISS_BVOC( I , J , 14) =  EMIS_FAXX  / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 15) =  EMIS_AAXX  / AREA_CM2 / DTSRCE     
               EMISS_BVOC( I , J , 16) =  EMIS_ALD2  / AREA_CM2 / DTSRCE     
               EMISS_BVOC( I , J , 17) =  EMIS_OMON  / AREA_CM2 / DTSRCE
               EMISS_BVOC( I , J , 18) =  EMIS_MOHX  / AREA_CM2 / DTSRCE     
               EMISS_BVOC( I , J , 19) =  EMIS_ETOH  / AREA_CM2 / DTSRCE     

               ! CHBr3 emissions [kg/m2/s] -- tracer #14
               IF ( IDECHBr3 > 0 ) THEN
                  AD46(I,J,14) = AD46(I,J,14) 
     &                         + ( EMISRR(I,J,IDECHBr3)  / AREA_M2 ) 
     &                         * ( MWT_CHBr3             / AVG     )
               ENDIF

               ! CH2Br2 emissions [kg/m2/s] -- tracer #15
               IF ( IDECHBr3 > 0 ) THEN
                  AD46(I,J,15) = AD46(I,J,15) 
     &                         + ( EMISRR(I,J,IDECH2Br2) / AREA_M2 ) 
     &                         * ( MWT_CH2Br2            / AVG     )
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE EMISSDR                                                

