!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: anthroems
!
! !DESCRIPTION: Subroutine ANTHROEMS reads anthropogenic tracers for each 
!  season.  NOx emissions at levels other than the surface are now accounted 
!  for.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ANTHROEMS( NSEASON, Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod,  ONLY : OptInput
      USE FUTURE_EMISSIONS_MOD,ONLY : GET_FUTURE_SCALE_ALK4ff
      USE FUTURE_EMISSIONS_MOD,ONLY : GET_FUTURE_SCALE_C2H6ff
      USE FUTURE_EMISSIONS_MOD,ONLY : GET_FUTURE_SCALE_C3H8ff
      USE FUTURE_EMISSIONS_MOD,ONLY : GET_FUTURE_SCALE_COff  
      USE FUTURE_EMISSIONS_MOD,ONLY : GET_FUTURE_SCALE_NOxff 
      USE FUTURE_EMISSIONS_MOD,ONLY : GET_FUTURE_SCALE_PRPEff
      USE FUTURE_EMISSIONS_MOD,ONLY : GET_FUTURE_SCALE_TONEff
      USE GEIA_MOD,            ONLY : READ_GEIA,    READ_C3H8_C2H6_NGAS
      USE GEIA_MOD,            ONLY : READ_LIQCO2,  READ_TODX
      USE GEIA_MOD,            ONLY : READ_TOTCO2,  TOTAL_FOSSIL_TG
      USE GRID_MOD,            ONLY : GET_AREA_CM2, GET_XOFFSET
      USE GRID_MOD,            ONLY : GET_YOFFSET
      USE TIME_MOD,            ONLY : GET_TS_EMIS,  GET_YEAR
      USE TIME_MOD,            ONLY : GET_SEASON
      USE TRACERID_MOD,        ONLY : IDEACET,      IDEALK4
      USE TRACERID_MOD,        ONLY : IDEC2H6,      IDEC3H8
      USE TRACERID_MOD,        ONLY : IDECO,        IDEMEK
      USE TRACERID_MOD,        ONLY : IDENO,       IDEPRPE
      USE TRACERID_MOD,        ONLY : NEMANTHRO
      USE TRACERID_MOD,        ONLY : IDEBENZ,     IDETOLU,    IDEXYLE
      USE TRACERID_MOD,        ONLY : IDEC2H4,      IDEC2H2
      USE TRACERID_MOD,        ONLY : IDTBENZ,     IDTTOLU,    IDTXYLE
      USE TRACERID_MOD,        ONLY : IDTC2H4,      IDTC2H2,   IDTNO2
      USE SCALE_ANTHRO_MOD,    ONLY : GET_ANNUAL_SCALAR
      USE SCALE_ANTHRO_MOD,    ONLY : GET_ANNUAL_SCALAR_05x0666_NESTED
      USE EDGAR_MOD,           ONLY : READ_AROMATICS, READ_C2H4
      USE EDGAR_MOD,           ONLY : READ_C2H2
      USE EDGAR_MOD,           ONLY : READ_AROMATICS_05x0666
      USE EDGAR_MOD,           ONLY : READ_C2H4_05x0666
      USE EDGAR_MOD,           ONLY : READ_C2H2_05x0666
      ! SOAupdate: for gas phase NAP chemistry (hotp 6/24/09)
      ! get location of NAP emissions in array
      USE TRACERID_MOD,        ONLY : IDENAP

      USE CMN_SIZE_MOD              ! Size parameters
      USE COMODE_LOOP_MOD           ! IDEMS
      USE CMN_O3_MOD                ! EMIST, EMISR, EMISRR, etc.

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: NSEASON     ! Current season (1-4)
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
!
! !REMARKS:
!  NSEASON:  is the seasonal index for NOx emissions:
!    NSEASON=1 --> winter (Dec, Jan, Feb)
!    NSEASON=2 --> spring (Mar, Apr, May)
!    NSEASON=3 --> summer (Jun, Jul, Aug)
!    NSEASON=4 --> autumn (Sep, Oct, Nov)
!                                                                             .
!  Passed Via CMN_O3_mod:
!  ===========================================================================
!  Fossil Fuel arrays: EMISTNOX,  EMISTCO,   EMISTETHE, EMISTPRPE,
!                      EMISTC2H6, EMISTC3H8, EMISTALK4, EMISTACET,
!                      EMISTMEK,  EMISTSOX
!                                                                             .
!  Emissions arrays:   EMIST, EMISTN, EMISR, EMISRN, EMISRR, EMISRRN
!                                                                             .
!  REFERENCES:
!  ===========================================================================
!  (1 ) Zhang, Y.X., and S. Tao, "Global atmospheric emission inventory of
!        polycyclic aromatic hydrocarbons (PAHs) for 2004", Atmos. Environ., 
!        Vol 43, 812-819, doi:10.1016/J.ATMOSENV.2008.10.050, 2009.
! 
! !REVISION HISTORY: 
!  04 Jun 1998 - R. Yantosca - Initial version
!  (1 ) We now read the new merge file, created for SASS. (bey, 2/99)
!  (2 ) ANTHROEMS should be called each time the season changes, since
!        the GEIA NOx emissions are seasonal.
!  (3 ) NOx emissions are stored separately in EMISTN, EMISRN, EMISRRN.  
!        This is because the NOx emissions can be located across several 
!        sigma levels, whereas the other tracers are only emitted into 
!        the surface level.
!  (4 ) NO2 is no longer emitted as the emission species for Ox.
!        (bey, bmy, 4/14/99)
!  (5 ) There are 3 different types of scale factors for anthro emissions:
!        (a) Yearly since 1985: done in anthroems.f
!        (b) Weekday/weekend:   done in emf_scale.f
!        (c) Time of day:       done in emfossil.f 
!  (6 ) At present NEMANTHRO = Total number of emitted tracers 
!        (set in tracerid.f).  We no longer use moments in emissions.
!        ORDER = NOx, CO, PRPE, C3H8, ALK4, C2H6, ALD2.
!  (7 ) NOx is assumed to be the first tracer (N=1).  The first usable 
!        row for tracers other than NOx in EMIST(I,J,N), etc. is N=2. 
!  (8 ) Need to offset EMISR, which has global dimensions.  
!        EMIST has window dimensions.
!  (9 ) Now trap I/O errors and stop gracefully if file open or read
!        errors are encountered.  Print an error message to alert user
!        which file triggered the I/O error. (bmy, 4/14/99)
!  (10) Eliminate GISS-specific code and PLUMES code (bmy, 4/14/99)
!  (11) Now use F90 syntax where expedient. (bmy, 4/14/99)
!  (12) Cosmetic changes, added comments (bmy, 3/17/00)
!  (13) Do not let SCALYEAR go higher than 1996, since right now we don't
!        have FF scaling data beyond 1996.  Also cosmetic changes and
!        updated comments. (bmy, 4/6/01)
!  (14) Now reference routines from GEIA_MOD for reading scale factor and
!        other emissions data from disk. (bmy, 4/23/01)
!  (15) Now read fossil-fuel emissions from a binary punch file (bmy, 4/23/01)
!  (16) CO and hydrocarbons are read from disk once per year.  Fossil fuel
!        scale factors are also applied once per 
!  (17) Now comment out LNAPAPNOX.  Also total fossil fuel emissions
!        and echo to std output. (bmy, 4/27/01)
!  (18) Bug fix: Now convert units for CO, Hydrocarbon tracers only once
!        per year.  Convert units for NOx once per season. (bmy, 6/7/01)
!  (19) Bug fix: Now index CH26 correctly when totaling it (bmy, 8/30/01)
!  (20) Now take C3H8 and C2H6 emissions as scaled from natural gas.  Read
!        these in subroutine READ_C3H8_C2H6_NGAS.  Also scale anthropogenic
!        ACET by 0.82 in order to match the acetone paper (bdf, bmy, 9/10/01)
!  (21) Removed obsolete, commented-out code from 6/01 (bmy, 11/26/01)
!  (22) Eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (23) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (24) Now reference IDTNOX, IDENOX, etc. from "tracerid_mod.f".  Also
!        do not let SCALEYEAR exceed 1998. (bmy, 1/13/03)
!  (25) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 from "grid_mod.f"
!        Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now I0 and J0 are local variables.  Now use functions GET_TS_EMIS,
!        GET_YEAR, GET_SEASON from "time_mod.f". (bmy, 2/11/03)
!  (26) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (27) Now replace FMOL with TRACER_MW_KG (bmy, 10/25/05)
!  (28) Modified for IPCC future emissions scale factors (swu, bmy, 5/30/06)
!  (29) Extend max value for FSCALYR to 2002 (bmy, 7/18/06)
!  (30) Use updated int'annual scale factors for 1985-2003 (amv, 08/24/07)
!  (31) As default, use EDGARv2.0 emission (fossil fuel + industry) 
!        for year 1985, scale to target year with CO2 from liquid fuel, 
!        for aromatics, C2H4, and C2H2. (tmf, 6/13/07)
!  (32) GET_ANNUAL_SCALAR_05x0666_NESTED_CH renamed to 
!        GET_ANNUAL_SCALAR_05x0666_NESTED (amv, bmy, 12/18/09)
!  19 Nov 2010 - R. Yantosca - Added ProTeX headers
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  25 Mar 2013 - R. Yantosca - Now use logical fields from Input_Opt
!  13 Aug 2013 - M. Sulprizio- Add NAP emissions for SOA + semivolatile POA
!                              simulation (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Remove SEAC4RS C-preprocessor switch
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE                 :: FIRST = .TRUE.
      INTEGER                       :: SCALEYEAR
      INTEGER, SAVE                 :: LASTYEAR
      INTEGER                       :: I, I0,  IREF, J, J0, JREF
      INTEGER                       :: K, L,   LL,   N, NN 
      REAL*8                        :: DTSRCE, AREA_CM2

      REAL*4 :: E_BENZ(IIPAR,JJPAR), E_TOLU(IIPAR,JJPAR),
     &          E_XYLE(IIPAR,JJPAR)
      REAL*4 :: E_C2H4(IIPAR,JJPAR), E_C2H2(IIPAR,JJPAR)
      REAL*8 :: GEOS1x1(I1x1,J1x1,1)
      REAL*8 :: TEMP(IIPAR,JJPAR)

      ! SOAupdate hotp
      ! for gas phase NAP chemistry (hotp 6/24/09)
      ! Get anthropogenic (FF) NAP emissions by scaling BENZ emissions with the
      ! following factor. Factor is ratio of TgC NAP to TgC BENZ emissions
      ! or equivalently, molec C NAP to molec C BENZ. Scaling should produce
      ! about 0.09 TgC NAP/year, consistent with non-BB,BF emissions predicted
      ! by Zhang and Tao 2009 Atm Env
      ! Based on year 2000 1x1 inv  (hotp 11/14/09)
      REAL*8, PARAMETER    :: NAPTOBENZSCALE = 0.06861d0

      ! NAPSOA: factor to scale total NAP emissions to POA (hotp 7/24/09)
      REAL*8, PARAMETER    :: NAPTOTALSCALE = 66.09027d0

      ! For values from Input_Opt
      LOGICAL :: LFUTURE
      INTEGER :: N_TRACERS
      REAL*8  :: TRACER_MW_KG(Input_Opt%N_TRACERS)
      REAL*8  :: NAPEMISS

      !=================================================================
      ! ANTHROEMS begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LFUTURE      = Input_Opt%LFUTURE
      N_TRACERS    = Input_Opt%N_TRACERS
      TRACER_MW_KG = Input_Opt%TRACER_MW_KG(1:N_TRACERS)
      NAPEMISS     = Input_Opt%NAPEMISS

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
      WRITE( 6, '(a)' ) 'A N T H R O P O G E N I C   E M I S S I O N S'
      WRITE( 6, '(a)' )
      WRITE( 6, 110   ) GET_YEAR(), GET_SEASON()
 110  FORMAT( 'ANTHROEMS: NYEAR, NSEASON = ', i4, 1x, i2 )

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Get nested-grid offsets
      I0     = GET_XOFFSET()
      J0     = GET_YOFFSET()

      ! As of March 2009, the GEIA input files for GEOS_5, 0.5X0.666, 
      ! NESTED_CHINA are already cropped to the nested grid domain.  
      ! So there is no need for offsetting the emission data. 
      ! Reset I0 = 0, J0 = 0. (tmf, 3/5/09)
#if   defined( GRID05x0666 ) && defined( NESTED_CH )

      I0 = 0
      J0 = 0

#elif defined( GRID025x03125 ) && defined( NESTED_NA )

      I0 = 0
      J0 = 0

#elif defined( GRID025x03125 ) && defined( NESTED_EU )

      I0 = 0
      J0 = 0

#endif

      !=================================================================
      ! If FSCALYR < 0 then use this year (JYEAR) for the scaling 
      ! factors.  Otherwise, use the value of FSCALYR as specified in 
      ! 'input.ctm'.
      !
      ! Do not let SCALEYEAR exceed 1998 for now, since this is the 
      ! latest year for which we have data from CDIAC. (bmy, 1/13/03)
      !
      ! Do not limit default SCALEYEAR - this is done in 
      !  GET_ANNUAL_SCALAR. Allow users to force the scaling year
      !  as before with a value GT 0 in input.geos (phs, 3/11/08)
      !=================================================================
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF

      !=================================================================      
      ! Do the following only on the very first call to ANTHROEMS...
      !=================================================================      
      IF ( FIRST ) THEN 
         
         ! Zero emission arrays
         EMISTNOX  = 0e0
         EMISTCO   = 0e0
         EMISTALK4 = 0e0
         EMISTACET = 0e0
         EMISTMEK  = 0e0
         EMISTPRPE = 0e0
         EMISTC3H8 = 0e0
         EMISTC2H6 = 0e0
         EMISTETHE = 0e0
         EMISTSOX  = 0e0
         EMISTBENZ = 0e0
         EMISTTOLU = 0e0
         EMISTXYLE = 0e0
         EMISTC2H4 = 0e0
         EMISTC2H2 = 0e0
         ! SOAupdate: for gas phase NAP chemistry (hotp 6/24/09)
         ! zero the NAP emission array
         EMISTNAP  = 0e0

         ! Zero arrays for holding CO & Hydrocarbons
         EMIST = 0d0
         EMISR = 0d0

         ! Read time-of-day scale factors (TODN, TODH, TODB)
         ! and weekday-weekend scale factors (SCNR89)
         CALL READ_TODX( TODN, TODH, TODB, SCNR89 )
      
         ! Read emissions from binary punch file format for entire year:
         ! NOx [molec NOx/cm2/s], CO [molec CO/cm2/s], HC's [atoms C/cm2/s]
         ! NOTE: We don't read in ETHE or SOx for our chemistry mechanism.
         CALL READ_GEIA( E_NOX  = EMISTNOX,  E_CO   = EMISTCO,   
     &                   E_ALK4 = EMISTALK4, E_ACET = EMISTACET,
     &                   E_MEK  = EMISTMEK,  E_PRPE = EMISTPRPE )

         ! Read C3H8 and C2H6 emissions, scaled from Natural Gas emissions
         ! as computed by Yaping Xiao (xyp@io.harvard.edu)
         CALL READ_C3H8_C2H6_NGAS( E_C3H8=EMISTC3H8, E_C2H6=EMISTC2H6 )

         !================================================================
         ! Read EDGARv2 aromatics emission for 1985  (tmf, 7/30/08)
         !================================================================
         IF ( IDTBENZ /= 0 .AND. IDTTOLU /= 0 .AND. IDTXYLE /= 0 ) THEN
#if defined(GRID05x0666) && defined( NESTED_CH )
         CALL READ_AROMATICS_05x0666( E_BENZ, E_TOLU, E_XYLE )
#else
         CALL READ_AROMATICS( E_BENZ, E_TOLU, E_XYLE )
#endif

            EMISTBENZ = E_BENZ
            EMISTTOLU = E_TOLU
            EMISTXYLE = E_XYLE
         ENDIF

         ! SOAupdate: for gas phase NAP chemistry (hotp 6/24/09)
         ! Get NAP FF emissions by scaling BENZ
         EMISTNAP  = E_BENZ * NAPTOBENZSCALE
         ! NAPSOA: scale up global total (hotp 7/24/09)
         EMISTNAP = EMISTNAP * NAPTOTALSCALE
         ! Now set according to input.geos (hopt 8/9/09)
         EMISTNAP = EMISTNAP * NAPEMISS

         !================================================================
         ! Read EDGARv2 C2H4 emission for 1985   (tmf, 7/30/08)
         !================================================================
         IF ( IDTC2H4 /= 0 ) THEN
#if defined(GRID05x0666) && defined( NESTED_CH )
         CALL READ_C2H4_05x0666( E_C2H4 )
#else
         CALL READ_C2H4( E_C2H4 )
#endif

            EMISTC2H4 = E_C2H4
         ENDIF
         !================================================================
         ! Read EDGARv2 C2H2 emission for 1985  (tmf, 7/30/08)
         !================================================================
         IF ( IDTC2H2 /= 0 ) THEN
#if defined(GRID05x0666) && defined( NESTED_CH )
         CALL READ_C2H2_05x0666( E_C2H2 )
#else
         CALL READ_C2H2( E_C2H2 )
#endif

            EMISTC2H2 = E_C2H2
         ENDIF
         !============================================================== 
         ! Apply IPCC future scale factors to emissions (if necessary)
         !==============================================================
         IF ( LFUTURE ) THEN

            ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Future CO [molec/cm2/s]
               EMISTCO(I,J)      = EMISTCO(I,J)                    *
     &                             GET_FUTURE_SCALE_COff( I, J )

               ! Future C2H6 [atoms C/cm2/s]
               EMISTC2H6(I,J)    = EMISTC2H6(I,J)                  *
     &                             GET_FUTURE_SCALE_C2H6ff( I, J )

               ! Future C3H8 emissions [atoms C/cm2/s]
               EMISTC3H8(I,J)    = EMISTC3H8(I,J)                  * 
     &                             GET_FUTURE_SCALE_C3H8ff( I, J )

               ! Future ALK4 [atoms C/cm2/s]
               EMISTALK4(I,J)    = EMISTALK4(I,J)                  * 
     &                             GET_FUTURE_SCALE_ALK4ff( I, J )

               ! Future PRPE [atoms C/cm2/s]
               EMISTPRPE(I,J)    = EMISTPRPE(I,J)                  *
     &                             GET_FUTURE_SCALE_PRPEff( I, J )

               ! Future ACET [atoms C/cm2/s]
               EMISTACET(I,J)    = EMISTACET(I,J)                  * 
     &                             GET_FUTURE_SCALE_TONEff( I, J )

               ! Future MEK [atoms C/cm2/s]
               EMISTMEK(I,J)     = EMISTMEK(I,J)                   * 
     &                             GET_FUTURE_SCALE_TONEff( I, J )

               ! Future NOx [molec/cm2/s]
               EMISTNOX(I,J,:,:) = EMISTNOX(I,J,:,:)               *
     &                             GET_FUTURE_SCALE_NOxff( I, J )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
      ENDIF

      !================================================================= 
      ! Do the following on the first call to ANTHROEMS,
      ! or whenever we enter into a new year:
      !================================================================= 
      IF ( FIRST .or. SCALEYEAR /= LASTYEAR ) THEN

         WRITE( 6, * )

         ! now use updated scalars (amv, phs, 3/11/08)
         CALL GET_ANNUAL_SCALAR( 71, 1985, SCALEYEAR, FTOTCO2 )
         CALL GET_ANNUAL_SCALAR( 72, 1985, SCALEYEAR, FLIQCO2 )

         ! Set SCALEYEAR to this YEAR
         LASTYEAR = SCALEYEAR

         !==============================================================
         ! Apply scale factors to CO and Hydrocarbon emission species
         ! These are aseasonal, so we only have to do this the first
         ! time that anthroems.f is called.  
         !
         ! EMIST(I,J,N) contains CO and hydrocarbon emission species
         ! in units of [molec (C)/cm2/s] 
         !
         ! NOTE: We always assume NOx is the first tracer (N=1), so the 
         ! first valid entry in EMIST(I,J,N) will be the N=2 row.  
         !==============================================================

         ! CO: Scale by liquid CO2 scale factors 
         IF ( IDECO /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDECO) = EMISTCO(IREF,JREF) * 
     &                               FLIQCO2(IREF,JREF)
               ENDDO
            ENDDO

            ! Print total in Tg CO
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDECO), IIPAR, JJPAR, 
     &                            1,                28d-3, 'CO' )
         ENDIF

         ! ALK4: scale by liquid fuel CO scale factors
         IF ( IDEALK4 /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDEALK4) = EMISTALK4(IREF,JREF) * 
     &                                 FLIQCO2(IREF,JREF)
               ENDDO
            ENDDO

            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEALK4), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'ALK4' )
         ENDIF

         ! ACET: scale by liquid fuel CO scale factors
         IF ( IDEACET /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  
                  ! Also multiply by 0.82 in order to match the
                  ! a posteriori acetone source (bdf, bmy, 9/5/01)
                  EMIST(I,J,IDEACET) = EMISTACET(IREF,JREF) * 
     &                                 FLIQCO2(IREF,JREF)   * 0.82d0
               ENDDO
            ENDDO

            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEACET), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'ACET' )
         ENDIF

         ! MEK: scale by liquid fuel CO scale factors 
         IF ( IDEMEK /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDEMEK) = EMISTMEK(IREF,JREF) * 
     &                                FLIQCO2(IREF,JREF)
               ENDDO
            ENDDO
         
            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEMEK), IIPAR, JJPAR, 
     &                            1,                 12d-3, 'MEK' )
         ENDIF

         ! PRPE: Scale by liquid CO2 scale factors 
         IF ( IDEPRPE /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDEPRPE) = EMISTPRPE(IREF,JREF) * 
     &                                 FLIQCO2(IREF,JREF)
               ENDDO
            ENDDO
         
            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEPRPE), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'PRPE' )
         ENDIF

         ! C3H8: scale by liquid fuel CO scale factors
         IF ( IDEC3H8 /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDEC3H8) = EMISTC3H8(IREF,JREF) * 
     &                                 FLIQCO2(IREF,JREF)
               ENDDO
            ENDDO
         
            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEC3H8), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'C3H8' )
         ENDIF

         ! C2H6: scale by liquid fuel CO scale factors
         IF ( IDEC2H6 /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDEC2H6) = EMISTC2H6(IREF,JREF) * 
     &                                 FLIQCO2(IREF,JREF)
               ENDDO 
            ENDDO

            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEC2H6), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'C2H6' )
         ENDIF         

         !=============================================================
         ! Default emissions for BENZ, TOLU, XYLE, C2H2, C2H4 
         !  are for year 1985 only.  Scale to target year
         !=============================================================
         ! BENZ: for year 1985
         IF ( IDEBENZ /= 0 .AND. IDTBENZ /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDEBENZ) = EMISTBENZ(IREF,JREF) * 
     &               FLIQCO2(IREF, JREF) 
               ENDDO
            ENDDO
         
            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEBENZ), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'BENZ' )
         ENDIF

         ! TOLU: for year 1985
         IF ( IDETOLU /= 0 .AND. IDTTOLU /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDETOLU) = EMISTTOLU(IREF,JREF) *  
     &               FLIQCO2(IREF, JREF) 
               ENDDO
            ENDDO
         
            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDETOLU), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'TOLU' )
         ENDIF

         ! XYLE: for year 1985
         IF ( IDEXYLE /= 0 .AND. IDTXYLE /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDEXYLE) = EMISTXYLE(IREF,JREF) * 
     &               FLIQCO2(IREF, JREF) 
               ENDDO
            ENDDO
         
            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEXYLE), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'XYLE' )
         ENDIF

         ! C2H4: for year 1985
         IF ( IDEC2H4 /= 0 .AND. IDTC2H4 /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDEC2H4) = EMISTC2H4(IREF,JREF) *
     &               FLIQCO2(IREF, JREF) 

               ENDDO
            ENDDO
         
            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEC2H4), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'C2H4' )
         ENDIF

         ! C2H2: for year 1985
         IF ( IDEC2H2 /= 0 .AND. IDTC2H2 /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMIST(I,J,IDEC2H2) = EMISTC2H2(IREF,JREF) *
     &               FLIQCO2(IREF, JREF) 

               ENDDO
            ENDDO
         
            ! Print total in Tg C
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDEC2H2), IIPAR, JJPAR, 
     &                            1,                  12d-3, 'C2H2' )
         ENDIF

         ! SOAupdate: for gas phase NAP chemsitry (hotp 6/24/09)
         ! NAP 
         IF ( IDENAP /= 0 ) THEN
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  ! scale this now since BENZ is also 1985 (hotp 7/28/10)
                  EMIST(I,J,IDENAP) = EMISTNAP(IREF,JREF) *
     &               FLIQCO2(IREF, JREF)
               ENDDO
            ENDDO
            
            ! Print total in Tg C 
            CALL TOTAL_FOSSIL_TG( EMIST(:,:,IDENAP), IIPAR, JJPAR, 
     &                            1,                 12d-3, 'NAP' )
         ENDIF

         !==============================================================
         ! Convert CO and hydrocarbon emissions from [molec (C)/cm2/s] 
         ! to [kg (C)/box/emission timestep].  Store in array EMISR.
         !==============================================================

         ! Loop over the anthropogenic tracers 
         DO N = 1, NEMANTHRO

            ! NN is the actual CTM tracer # 
            ! corresponding to emissions species N
            NN = IDEMS(N)

            ! Skip NOx
            ! Replaced NOx with NO (mpayer, 3/14/13)
            IF ( N == IDENO ) CYCLE

            ! Skip if some tracer is not present because there are more 
            ! anthro. tracers for dicarbonyl chemistry. (ccc, 4/16/09)
            IF ( NN == 0 ) CYCLE

            ! Convert units
            DO J = 1, JJPAR
               JREF = J + J0
               
               DO I = 1, IIPAR
                  IREF = I + I0
                  
                  ! Grid box surface area [cm2]
                  AREA_CM2 = GET_AREA_CM2( I, J, 1 )

                  EMIST(I,J,N) = EMIST(I,J,N) * TRACER_MW_KG(NN) * 
     &                           DTSRCE * AREA_CM2 / 6.023d23

                  EMISR(IREF,JREF,N) = EMIST(I,J,N)
               ENDDO    
            ENDDO       
         ENDDO  

      ENDIF   ! FIRST or SCALEYEAR /= LASTYEAR

      !==============================================================
      ! Apply total fuel CO2 scale factors to NOx emissions
      ! This has to be done once per season (4x/year); 
      ! that is, every time that ANTHROEMS is called. 
      !==============================================================

      ! Zero NOx emission arrays
      EMISTN = 0d0
      EMISRN = 0d0

      ! NOX: scale by total CO2 scale factors
      ! Now loop over 2 levels only (phs, 10/19/09)
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( IDENO > 0 ) THEN
         DO LL = 1, 2 !NOXLEVELS
            DO J = 1, JJPAR
               JREF = J + J0
               DO I = 1, IIPAR
                  IREF = I + I0
                  EMISTN(I,J,LL) = EMISTNOX(IREF,JREF,NSEASON,LL) * 
     &                             FTOTCO2(IREF,JREF)
               ENDDO
            ENDDO
         ENDDO
         
         ! Print total in Tg N
         ! Now loop over 2 levels only (phs, 10/19/09)
         CALL TOTAL_FOSSIL_TG( EMISTN,    IIPAR, JJPAR,  
     &                         2,         14d-3, 'NO', NSEASON )
      ENDIF   

      !=================================================================
      ! Convert all emission species from [molec (C)/cm2/s] to 
      ! [kg/box/emission timestep] and store in EMISRN, EMISR arrays.
      !=================================================================

      ! Loop over the anthropogenic tracers
      DO N = 1, NEMANTHRO

         ! NN is the actual CTM tracer # 
         ! corresponding to emissions species N
         NN = IDEMS(N)

         ! Do unit conversion for NOx separately, since it is multi-level
         ! Now loop over 2 levels only (phs, 10/19/09)
         ! Replaced NOx with NO (mpayer, 3/14/13)
         IF ( N == IDENO ) THEN
            DO LL = 1, 2 !NOXLEVELS
               DO J = 1, JJPAR
                  JREF = J + J0
               
                  DO I = 1, IIPAR
                     IREF = I + I0

                     ! Grid box surface area [cm2]
                     AREA_CM2 = GET_AREA_CM2( I, J, 1 )
               
                     ! Convert NOx using molecular weight of NO2 since original
                     ! NOx emissions are in mass units of NO2 (mpayer, 4/18/13)
                     EMISTN(I,J,LL) = EMISTN(I,J,LL) *
     &                                TRACER_MW_KG(IDTNO2) * 
     &                                DTSRCE * AREA_CM2  / 6.023d23
               
                     EMISRN(IREF,JREF,LL) = EMISTN(I,J,LL)
                  ENDDO
               ENDDO
            ENDDO

            ! Exit from the loop over anthropogenic tracers
            EXIT
         ENDIF
      ENDDO  

      !================================================================
      ! Cleanup and quit
      !================================================================

      ! Set first-time-flag FALSE for next iteration
      FIRST = .FALSE.

      ! Pretty output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      END SUBROUTINE ANTHROEMS
!EOC

