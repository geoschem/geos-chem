!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: emfossil
!
! !DESCRIPTION: Subroutine EMFOSSIL emits fossil fuels into the EMISRR and 
!  EMISRRN arrays, which are then passed to SMVGEAR. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMFOSSIL( I, J, N, NN, IREF, JREF, JSCEN, Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE BRAVO_MOD,          ONLY : GET_BRAVO_ANTHRO, GET_BRAVO_MASK
      USE CAC_ANTHRO_MOD,     ONLY : GET_CANADA_MASK,  GET_CAC_ANTHRO
      !USE DAO_MOD,            ONLY : IS_WATER
      USE DIAG_MOD,           ONLY : AD29,   AD32_an,  AD36
      USE EDGAR_MOD,          ONLY : GET_EDGAR_CO,     GET_EDGAR_NOx
      USE EDGAR_MOD,          ONLY : GET_EDGAR_TODN
      USE EMEP_MOD,           ONLY : GET_EMEP_ANTHRO,  GET_EUROPE_MASK
      USE EPA_NEI_MOD,        ONLY : GET_EPA_ANTHRO,   GET_USA_MASK
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE RETRO_MOD,          ONLY : GET_RETRO_ANTHRO
      USE C2H6_MOD,           ONLY : GET_C2H6_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : GET_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : NEI05_MASK => USA_MASK
      USE LOGICAL_MOD,        ONLY : LICOADSSHIP !(cklee, 6/30/09)
      USE RCP_MOD,            ONLY : GET_RCP_EMISSION
                              
      USE STREETS_ANTHRO_MOD, ONLY : GET_SE_ASIA_MASK
      USE STREETS_ANTHRO_MOD, ONLY : GET_STREETS_ANTHRO
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TIME_MOD,           ONLY : GET_DAY_OF_WEEK_LT
      USE TIME_MOD,           ONLY : GET_HOUR
      USE TRACERID_MOD,       ONLY : IDENO,    IDEO3,    IDEHNO3
      USE TRACERID_MOD,       ONLY : IDTCO,    IDTHNO3
      USE TRACERID_MOD,       ONLY : IDTC2H6,  IDTNO2
      USE VISTAS_ANTHRO_MOD,  ONLY : GET_VISTAS_ANTHRO
      USE ICOADS_SHIP_MOD,    ONLY : GET_ICOADS_SHIP !(cklee, 7/09/09)

      USE CMN_SIZE_MOD             ! Size parameters
      USE COMODE_LOOP_MOD          ! IHOUR
      USE CMN_O3_MOD               ! EMISR, EMISRR, etc...
      USE CMN_DIAG_MOD             ! Diagnostic switches & arrays

      ! 10/24/12, ckeller: NOX diurnal scale factors fix:
      USE TIME_MOD,           ONLY : GET_LOCALTIME
 
      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I       ! GEOS-Chem longitude index
      INTEGER, INTENT(IN)   :: J       ! GEOS-Chem latitude index
      INTEGER, INTENT(IN)   :: N       ! GEOS-Chem emission species index
      INTEGER, INTENT(IN)   :: NN      ! GEOS-Chem advected tracer index
      INTEGER, INTENT(IN)   :: IREF    ! Offset index I+I0
      INTEGER, INTENT(IN)   :: JREF    ! Offset index J+J0
      INTEGER, INTENT(IN)   :: JSCEN   ! Day index (Sat=1, Sun=2, Weekday=3)
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
!
! !REMARKS:
!  In most cases, I0=J0=0, so IREF=I and JREF=J.  The offsets I0 and J0 are
!  mostly historical baggage.
!                                                                             .
!  NOTE: The source code for ship emissions has been commented out, but left
!  in place.  Although PARANOX computes ship emissions in CALCRATE, we may
!  have to disable this for the Grid-Independent model (because we have to
!  have a clean separation between emissions & chemistry).
! 
! !REVISION HISTORY: 
!  19 Apr 1999 - R. Yantosca - Initial version
!  (1 ) Uses the correct seasonal NOx and multi-level NOx (anthroems.f)
!  (2 ) Uses anthro scale factors for years since 1985 (from anthroems.f)
!  (3 ) Scales emissions based on weekday/weekend (emf_scale.f)
!  (4 ) Preserves old sensitivity study cases (emf_scale.f, emissdr.f)
!  (5 ) Scales emissions based on time of day (emfossil.f)
!  (6 ) Get rid of all GISS and PLUMES code (bmy, 4/19/99)
!  (7 ) Now use F90 syntax for declarations, etc. (bmy, 4/19/99)
!  (8 ) Now use allocatable arrays for ND29 and ND36 diagnostics.
!        Also made minor cosmetic changes & updated comments. (bmy, 3/16/00)
!  (9 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (10) Enhance anthropogenic CO emission by 8%, to account for CO production
!        from oxidation of anthropogenic VOC's (bnd, bmy, 1/2/01)
!  (11) Comment out scaling by 1.08 for anthro CO (bmy, 2/12/01)
!  (12) Eliminate obsolete commented-out code (bmy, 4/20/01)
!  (13) Now use 2% as the enhancment factor for CO instead of 1.08,
!        according to new jal numbers (bmy, 4/26/01)
!  (14) Now references "tracerid_mod.f" (bmy, 11/6/02)
!  (15) Now replaced DXYP(JREF)*1d4 with GET_AREA_CM2(J).  Now use function
!        GET_TS_EMIS() from "time_mod.f" (bmy, 2/11/03)
!  (16) Now can overwrite existing emissions with EPA/NEI data over the 
!        continental USA if LNEI99=T.  Now reference LNEI99 from F90
!        module "logical_mod.f".  Now reference GET_EPA_ANTHRO and
!        GET_USA_MASK from "epa_nei_mod.f". (rch, rjp, bmy, 11/5/04)
!  (17) Now references GET_DAY_OF_WEEK from "time_mod.f" to correctly figure
!        out if this is a weekday or weekend. (bmy, 7/6/05)
!  (18) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (19) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (20) Now apply EMEP European emissions if necessary.  Remove reference
!        to CMN, it's now obsolete. (bdf, bmy, 11/1/05)
!  (21) Rewrite IF statements to avoid seg fault errors when LEMEP and LNEI99 
!        are turned off. (bmy, 2/1/06)
!  (22) Now apply BRAVO Mexican emissions if necessary (rjp, kfb, bmy, 6/26/06)
!  (23) Now apply EDGAR emissions if necessary.  Also now only do the the 
!        EDGAR, EPA, EMEP, and BRAVO function calls in the LL=1 block.
!        (avd, bmy, 7/10/06)
!  (24) Now do BRAVO emissions before EPA/NEI99 emissions in order to avoid 
!        zero emissions in some boxes.  Now add David Streets emissions for 
!        NOx over SE Asia and CO over just China (yxw, bmy, 8/17/06)
!  (25) Bug fix: Now only execute EDGAR CO block if the tracer is CO.
!        Also, David Streets' CO is now applied over SE ASIA. (bmy, 9/8/06)
!  (26) Now references ITS_A_TAGCO_SIM from "tracer_mod.f".  Enhance CO prod
!        by 18.5% for tagged CO sim here instead of in "tagged_co_mod.f".
!        (bmy, 2/14/08)
!  (27) Use more robust test to only screen out "missing" values in EMEP,
!        BRAVO, and David Streets emissions. (avd, phs, bmy, 11/19/08)   
!  (28) Ship NOx is emitted as HNO3+10*O3  (phs, 3/4/08)
!  (29) Apply spatially-varying diurnal scalars for NOx (amv, 08/24/07)
!  (30) Now apply CAC Canadian emissions if necessary (amv, 01/09/08)
!  (31) Moved down BRAVO parts and add BRAVO and EPA emissions where they 
!        overlap (phs, 5/7/08)
!  (32) Now overwrite USA NOx with VISTAS if necessary (amv, 12/02/08)
!  (33) Modified CO scaling (jaf, 2/25/09)
!  (34) Add a test on existing emissions for EPA/NEI. (hotp, ccc, 5/29/09)
!  (35) Updated ship treatment (phs, 7/0/09)
!  (36) Add NEI2005 (amv, phs, 10/20/09)
!  (37) Bug fix for tagged CO and 0.5 x 0.666 Nested Grid (yxw, bmy, 11/23/09)
!  (38) Bug fix for array EMISRR, if emissions are already present in this 
!        array (e.g. ship O3 or HNO3) they no longer get overwritten.
!        (gvinken, 11/16/10)
!  19 Nov 2010 - R. Yantosca - Added ProTeX headers
!  24 Feb 2012 - M. Payer    - Commented out ship emissions, which has been
!                              moved to calcrate.F for PARANOX.  Left original
!                              code in place, but commented out.
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  22 Mar 2012 - M. Payer    - Update for C2H6. RETRO ethane emissions are too
!                              low so we will use Yaping Xiao's offline emiss.
!  24 Oct 2012 - C. Keller   - Use localtime instead of UTC to derive NOx 
!                              diurnal scale factors.
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox emissions with NO and O3
!                              emissions as part of removal of NOx-Ox
!                              partitioning
!  18 Apr 2013 - M. Payer    - Convert NOx using XNUMOL for NO2 since original
!                              NOx emissions are in mass units of NO2
!  25 Mar 2013 - R. Yantosca - Now use logical fields from Input_Opt
!  14 Jun 2013 - R. Yantosca - Now determine weekday/weekend with respect to
!                              the local time at each grid box.  (Formerly,
!                              this had been done w/r/t the GMT time.
!  22 Jul 2013 - M. Sulprizio- Now copy LRCP and LRCPSHIP from Input_Opt
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: WEEKDAY
      INTEGER :: L, LL, K, DOW_LT, HOUR
      REAL*8  :: TODX, DTSRCE, AREA_CM2           
      REAL*8  :: EMX(NOXLEVELS)  
      REAL*8  :: XEMISR
      REAL*8  :: XEMISRN(NOXLEVELS)
      REAL*8  :: BRAVO, EPA_NEI, EMEP,   EDGAR, STREETS
      REAL*8  :: CAC,   SHIP,    VISTAS, NEI05
      REAL*8  :: RETRO
      REAL*8  :: C2H6_ANTHRO
      ! 10/24/12, ckeller: NOX diurnal scale factors fix:
      INTEGER :: NOXHOUR
      REAL*8  :: RCP

      ! For fields from Input_Opt
      LOGICAL :: LBRAVO
      LOGICAL :: LEMEP
      LOGICAL :: LNEI99
      LOGICAL :: LEDGARNOx
      LOGICAL :: LEDGARCO
      LOGICAL :: LSTREETS
      LOGICAL :: LCAC
      LOGICAL :: LEDGARSHIP
      LOGICAL :: LARCSHIP
      LOGICAL :: LEMEPSHIP
      LOGICAL :: LVISTAS
      LOGICAL :: LICARTT
      LOGICAL :: LNEI05
      LOGICAL :: LRETRO
      LOGICAL :: LHIST
      LOGICAL :: LRCP
      LOGICAL :: LRCPSHIP
      LOGICAL :: ITS_A_TAGCO_SIM
      INTEGER :: N_TRACERS
      REAL*8  :: XNUMOL(Input_Opt%N_TRACERS)

      !=================================================================
      ! EMFOSSIL begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LBRAVO          = Input_Opt%LBRAVO
      LEMEP           = Input_Opt%LEMEP 
      LNEI99          = Input_Opt%LNEI99
      LEDGARNOx       = Input_Opt%LEDGARNOx
      LEDGARCO        = Input_Opt%LEDGARCO
      LSTREETS        = Input_Opt%LSTREETS  
      LCAC            = Input_Opt%LCAC
      LEDGARSHIP      = Input_Opt%LEDGARSHIP
      LARCSHIP        = Input_Opt%LARCSHIP
      LEMEPSHIP       = Input_Opt%LEMEPSHIP
      LVISTAS         = Input_Opt%LVISTAS
      LICARTT         = Input_Opt%LICARTT
      LNEI05          = Input_Opt%LNEI05
      LRETRO          = Input_Opt%LRETRO
      LHIST           = Input_Opt%LHIST
      LRCP            = Input_Opt%LRCP
      LRCPSHIP        = Input_Opt%LRCPSHIP
      ITS_A_TAGCO_SIM = Input_Opt%ITS_A_TAGCO_SIM
      N_TRACERS       = Input_Opt%N_TRACERS
      XNUMOL          = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Emission timestep [s]
      DTSRCE          = GET_TS_EMIS() * 60d0

      ! Surface area of grid box
      AREA_CM2        = GET_AREA_CM2( I, J, 1 )

      ! GMT hour of day
      HOUR            = GET_HOUR()

      ! Determine if we should use weekday or weekend NEI
      ! emissions at grid box (I,J,L).  Since NEI is over
      ! the US, then weekend is Sat/Sun.
      DOW_LT          = GET_DAY_OF_WEEK_LT( I, J, 1 )
      WEEKDAY         = ( DOW_LT > 0 .and. DOW_LT < 6 )

      !=================================================================
      ! Call EMF_SCALE to do the following:
      ! (1) Save original values of EMISR, EMISRN
      ! (2) If LFFNOX=F, turn off NOx, Ox emissions
      ! (3) Scale emissions to weekend/weekday usage
      !=================================================================
      CALL EMF_SCALE( I, J, N, NN, IREF, JREF, JSCEN, XEMISR, XEMISRN )       

      !=================================================================    
      ! ADD ANTHROPOGENIC EMISSIONS TO TRACER TOTALS
      ! NOTE APPROPRIATE TIME-OF-DAY FACTOR (TOD) MUST BE
      ! ESTABLISHED FOR EACH TRACER; 
      ! WITH IHOUR = 1-6 (1 = 10pm-2am)
      ! and tracer index distinguishing NOx-HC- BIO
      !
      ! NOx only: account for all NOx levels (LL=1,NOXLEVELS)
      !=================================================================    
      ! Replaced NOx emissions with NO emissions (mpayer, 3/14/13)
      IF ( N == IDENO ) THEN

         ! Initialize work variables
         EMX(:)  = 0d0

         ! 10/24/12, ckeller: fix for EDGAR diurnal scale factors:
         NOXHOUR = MIN( 23, NINT( GET_LOCALTIME ( I, J, 1 ) ) )
         TODX = GET_EDGAR_TODN(I,J,NOXHOUR)

         ! Use spatially varying diurnal scale factors 
         ! from EDGAR (amv, phs, 3/10/08)
!         TODX = GET_EDGAR_TODN(I,J,HOUR)

         ! Loop over all of the emission levels for NOx (e.g. surface, 100m)
         DO LL = 1, NOXLEVELS
            EMX(LL)  = TODX * EMISRN(IREF,JREF,LL)

            !-----------------------------------------------------------
            ! Get NOx from the EDGAR or RCP inventory (global)
            !-----------------------------------------------------------

            ! If we are using EDGAR emissions
            IF ( LEDGARNOx ) THEN

               ! Put all emissions into 1st level
               IF ( LL == 1 ) THEN

                  ! Get EDGAR emissions for NOx [molec/cm2/s]
                  EDGAR   = GET_EDGAR_NOx( I, J, MOLEC_CM2_S=.TRUE. )

                  ! Apply EDGAR time-of-day factor
                  EDGAR   = EDGAR * TODX

                  ! Replace GEIA with EPA/NEI emissions at surface
                  EMX(LL) = EDGAR * ( DTSRCE * AREA_CM2 ) /
     &                      XNUMOL(IDTNO2) 
            
               ELSE 
            
                  ! Zero EDGAR emissions in the 2nd level 
                  EMX(LL) = 0d0                   
                  
               ENDIF

            ! If we are using RCP emissions (cdh, 10/14/11)
            ELSEIF ( LRCP ) THEN

               ! Put all emissions into 1st level
               IF ( LL == 1 ) THEN
                  
                  ! Get RCP emissions for NOx [molec/cm2/s]
                  RCP   = GET_RCP_EMISSION( I, J, NN, 
     &                                      LAND=.TRUE., SHIP=.FALSE. )

                  ! Apply EDGAR time-of-day factor
                  RCP   = RCP * TODX

                  ! Replace GEIA with RCP emissions at surface
                  EMX(LL) = RCP * ( DTSRCE * AREA_CM2 ) / XNUMOL(IDTNO2) 

               ELSE 
            
                  ! Zero RCP emissions in the 2nd level 
                  EMX(LL) = 0d0                   
                  
               ENDIF                  
            ENDIF

            !-----------------------------------------------------------
            ! Get NOx from EMEP inventory over Europe 
            !-----------------------------------------------------------

            ! If we are using EMEP ...
![eml Don't apply regional overwrites of NOx over global historical emissions
            IF ( LEMEP .and. .not. LHIST ) THEN
!eml]
               ! If we are over the European region ...
               IF ( GET_EUROPE_MASK( I, J ) > 0d0 ) THEN 

                  IF ( LL == 1 ) THEN
            
                     ! Get EMEP emissions for NOx 
                     EMEP    = GET_EMEP_ANTHRO( I, J, NN, KG_S=.FALSE. )

                     ! Apply time-of-day factor
                     EMEP    = EMEP * TODX

                     ! Replace GEIA with EMEP emissions at surface
                     EMX(LL) = EMEP * ( DTSRCE * AREA_CM2 ) /
     &                         XNUMOL(IDTNO2) 
                  
                  ELSE 
            
                     ! Zero GEIA emissions in the 2nd level 
                     ! where the EMEP emissions are nonzero
                     EMX(LL) = 0d0                   
                  
                  ENDIF

               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! Get NOx from EPA/NEI or VISTAS inventory over the USA 
            !-----------------------------------------------------------

            ! If we are using EPA/NEI emissions
![eml 
            IF ( LNEI99 .and. .not. LHIST ) THEN
!eml]

               ! If we are over the USA ...               
               IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN 
            
                  IF ( LL == 1 ) THEN
                  
                     ! Get EPA emissions for NOx 
                     EPA_NEI = GET_EPA_ANTHRO( I, J, NN, WEEKDAY )

                     ! Apply time-of-day factor
                     EPA_NEI = EPA_NEI * TODX

                     ! Replace GEIA with EPA/NEI emissions at surface
                     EMX(LL) = EPA_NEI * 
     &                         ( DTSRCE * AREA_CM2 ) / XNUMOL(IDTNO2) 
            
                  ELSE 
            
                     ! Zero GEIA emissions in the 2nd level 
                     ! where the EPA/NEI emissions are nonzero
                     EMX(LL) = 0d0                   
                  
                  ENDIF
               ENDIF
            ENDIF

![eml
            IF ( LVISTAS .and. .not. LHIST ) THEN
!eml]
               ! If we are over the USA ...
               IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN

                  IF ( LL == 1 ) THEN

                     ! Get VISTAS emissions for NOx 
                     VISTAS = GET_VISTAS_ANTHRO( I, J, NN, WEEKDAY )

                     ! Apply time-of-day factor
                     VISTAS = VISTAS * TODX

                     ! Replace with VISTAS emissions at surface
                     EMX(LL) = VISTAS * 
     &                         ( DTSRCE * AREA_CM2 ) / XNUMOL(IDTNO2)

                  ELSE 

                     EMX(LL) = 0d0

                  ENDIF
               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! Get NOx from the David Streets' inventory (SE Asia)
            !-----------------------------------------------------------

            ! If we are using David Streets' emissions
![eml
            IF ( LSTREETS .and. .not. LHIST ) THEN
!eml]
               ! If we are over the SE Asia region
               IF ( GET_SE_ASIA_MASK( I, J ) > 0d0 ) THEN

                  ! Put all emissions into 1st level
                  IF ( LL == 1 ) THEN

                     ! Get David Streets' emissions for NOx [molec/cm2/s]
                     STREETS = GET_STREETS_ANTHRO( I, J, NN, 
     &                                             MOLEC_CM2_S=.TRUE. )

                     ! Apply time-of-day factor
                     STREETS = STREETS * TODX

                     ! Replace base emissions with STREETS
                     EMX(LL) = STREETS * 
     &                         ( DTSRCE * AREA_CM2 ) / XNUMOL(IDTNO2) 
            
                  ELSE 
            
                     ! Zero EDGAR emissions in the 2nd level 
                     EMX(LL) = 0d0                   
                  
                  ENDIF
               ENDIF
            ENDIF


            !-----------------------------------------------------------
            ! Get NOx from BRAVO inventory over MEXICO
            !-----------------------------------------------------------

            ! If we are using BRAVO ...
![eml
            IF ( LBRAVO .and. .not. LHIST ) THEN
!eml]
               ! If we are over the Mexican region ...
               IF ( GET_BRAVO_MASK( I, J ) > 0d0 ) THEN 

                  IF ( LL == 1 ) THEN

                     ! Get BRAVO emissions for NOx 
                     ! (and apply time-of-day factor)
                     BRAVO   = GET_BRAVO_ANTHRO( I, J, NN ) * TODX

                     ! Replace GEIA with BRAVO emissions at surface
                     ! Now, if on border, add to NEI99 emissions (phs, 5/7/08)
                     IF ( LNEI99 ) THEN 
                        IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN

                           EMX(LL) = EMX(LL) + BRAVO * ( DTSRCE *
     &                               AREA_CM2 ) / XNUMOL(IDTNO2) 
                        ENDIF                           
                     ELSE
                        EMX(LL) = BRAVO * ( DTSRCE*AREA_CM2 ) /
     &                            XNUMOL(IDTNO2) 

                     ENDIF

                  ELSE 
            
                     ! Zero GEIA emissions in the 2nd level 
                     ! where the BRAVO emissions are nonzero
                     EMX(LL) = 0d0                   
                  
                  ENDIF
               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! Get NOx from the CAC inventory (Canada)
            !-----------------------------------------------------------

            ! If we are using CAC emissions
![eml
            IF ( LCAC .and. .not. LHIST ) THEN
!eml]
               ! If we are over the SE Asia region
               IF ( GET_CANADA_MASK( I, J ) > 0d0 ) THEN

                  ! Put all emissions into 1st level
                  IF ( LL == 1 ) THEN

                     ! Get CAC emissions for NOx [molec/cm2/s]
                     CAC = GET_CAC_ANTHRO( I, J, NN,
     &                                             MOLEC_CM2_S=.TRUE. )

                     ! Apply time-of-day factor
                     CAC = CAC * TODX

                     IF ( LNEI99 ) THEN

                        ! If on border, add to NEI99 emissions (which has
                        ! no Canadian component)
                        IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
                           EMX(LL) = EMX(LL) + CAC * ( DTSRCE * 
     &                               AREA_CM2 ) / XNUMOL(IDTNO2)
                        ELSE
                           EMX(LL) = CAC * ( DTSRCE * AREA_CM2 ) /
     &                               XNUMOL(IDTNO2)
                        ENDIF
                     ELSE

                        ! Replace base emissions with CAC
                        EMX(LL) = CAC * ( DTSRCE * AREA_CM2 ) /
     &                            XNUMOL(IDTNO2)
                     ENDIF

                  ELSE

                     ! Zero CAC emissions in the 2nd level
                     EMX(LL) = 0d0

                  ENDIF
               ENDIF
            ENDIF

![eml
            IF ( LNEI05 .and. .not. LHIST ) THEN
!eml]

               ! If we are over the USA and CAN/MEX
               IF ( NEI05_MASK( I, J ) > 0d0 ) THEN

                  ! Get EPA emissions for NOx
                  NEI05 = GET_NEI2005_ANTHRO( I, J, LL, NN, WEEKDAY,
     &                                        MOLEC_CM2_s = .TRUE.)

                  ! Apply time-of-day factor
                  NEI05 = NEI05 * TODX

                  ! Replace GEIA with EPA/NEI emissions at surface
                  EMX(LL) = NEI05 *
     &                      ( DTSRCE * AREA_CM2 ) / XNUMOL(IDTNO2)

               ENDIF

            ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 2/24/12:
! Comment this out and leave this here, in case we have to re-install it
! in the context of the Grid-Independent model (bmy, 2/24/12) 
!
!            !SHIP EMISSIONS NOW DONE IN CALCRATE.F (gvinken, 07/27/10)
!
!            !-----------------------------------------------------------
!            ! Add ship emissions emitted as HNO3 and 10*O3,
!            ! i.e., ozone production efficiency (OPE)=10.
!            ! See : Chen, G., et al. (2005), An investigation of the 
!            ! chemistry of ship emission plumes during ITCT 2002, 
!            ! J. Geophys. Res., 110, D10S90, doi:10.1029/2004JD005236.
!            ! (djj, phs, 3/4/08)
!            ! Now also process EMEP NOx ship emissions, available
!            ! from 1990 with EMEP 2005 (phs, 6/08)
!            ! Correctly handle LEMEPSHIP=.TRUE. (phs 7/9/09) 
!            !-----------------------------------------------------------
!
!
!            ! DO it only once (1st level)
!            IF ( LL == 1 ) THEN
!
!               ! Reset
!               SHIP = 0D0
!
!               ! handle global inventory first
!               IF ( LEDGARSHIP ) THEN 
!
!                  ! Get SHIP EDGAR emissions for NOx [molec/cm2/s]
!                  SHIP = GET_EDGAR_NOx( I, J, 
!     &                                  MOLEC_CM2_S=.TRUE., SHIP=.TRUE.)
!                  
!               ! ICOADS ship emissions (cklee,7/09/09)
!               ELSE IF ( LICOADSSHIP ) THEN
!
!                  ! Get ICOADS  emissions for NOx [molec/cm2/s]
!                  SHIP = GET_ICOADS_SHIP( I, J, NN, MOLEC_CM2_S=.TRUE. )
!                  
!               ENDIF
!
!               ! Overwrite Europe
!               IF ( LEMEPSHIP ) THEN
!                 
!                  IF ( GET_EUROPE_MASK( I, J ) > 0d0 )
!
!                  ! Get SHIP EMEP emissions for NOx [molec/cm2/s]
!     &            SHIP = GET_EMEP_ANTHRO( I, J, NN, SHIP=.TRUE.)
!
!               ENDIF
!                  
!               ! Store as HNO3 and O3
!               ! Convert molec/cm2/s to molec/box/s (cdh, 10/20/2011)
!               EMISRR(I,J,IDEHNO3) = SHIP * AREA_CM2
!               EMISRR(I,J,IDEOX)   = 10D0 * SHIP * AREA_CM2 
!               
!               ! ND36 = Anthro source diagnostic...store as [molec/cm2]
!               ! and convert to [molec/cm2/s] in DIAG3.F
!               IF ( ND36 > 0 ) THEN
!
!                  AD36(I,J,IDEHNO3) = AD36(I,J,IDEHNO3) + SHIP * DTSRCE
!
!                  AD36(I,J,IDEOX) = AD36(I,J,IDEOX) + 10D0 * SHIP * 
!     &                              DTSRCE
!                  
!               ENDIF
!                  
!            ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !-----------------------------------------------------------
            ! Store in EMISRRN array and archive diagnostics
            !-----------------------------------------------------------

            ! EMISRRN [molec/box/s] is referenced by LL
            EMISRRN(I,J,LL) = EMISRRN(I,J,LL) +
     &                       ( EMX(LL) * XNUMOL(IDTNO2) / DTSRCE )

            ! ND32 = save anthro NOx for levels L=1,NOXEXTENT [molec/cm2/s]
            IF ( ND32 > 0 ) THEN
               AD32_an(I,J,LL) = AD32_an(I,J,LL) + 
     &              ( EMX(LL) * XNUMOL(IDTNO2) / ( DTSRCE * AREA_CM2 ) )   
            ENDIF

            ! ND36 = save anthro emissions in [molec/cm2]
            ! and then convert to [molec/cm2/s] in DIAG3.F
            IF ( ND36 > 0 ) THEN
               AD36(I,J,N) = AD36(I,J,N) +
     &            ( EMX(LL) * XNUMOL(IDTNO2) / AREA_CM2 )
            ENDIF    
         ENDDO     

      !=================================================================    
      ! All other emitted tracers except NOx! 
      !=================================================================    
      ELSE

         ! Initialize work variables
         EMX(:)  = 0d0

         ! Use appropriate scale factor for time of day
         ! Replaced Ox emissions with O3 emissions (mpayer, 3/14/13)
         IF ( N == IDEO3 ) THEN
            TODX = TODN(IHOUR)
         ELSE
            TODX = TODH(IHOUR)
         ENDIF

         EMX(1) = TODX * EMISR(IREF,JREF,N) 
 
         !--------------------------------------------------------------
         ! Get CO emissions from the EDGAR inventory (global)
         !--------------------------------------------------------------

         ! If we are using EDGAR CO ...
         IF ( NN == IDTCO .and. LEDGARCO ) THEN
         
            ! Get EDGAR CO
            EDGAR  = GET_EDGAR_CO( I, J, MOLEC_CM2_S=.TRUE. )

            ! Apply time of day factor
            EDGAR  = EDGAR * TODX
         
            ! Convert from molec/cm2/s to kg/box/timestep in order
            ! to be in the proper units for EMISRR array
            EMX(1) = EDGAR * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 

         ENDIF

         !--------------------------------------------------------------
         ! Get hydrocarbon emissions from RETRO inventory
         !--------------------------------------------------------------
        
         ! If we are using RETRO emissions ...
         IF ( LRETRO ) THEN

	    ! Get RETRO emissions
            RETRO = GET_RETRO_ANTHRO( I, J, NN )

	    ! -1 indicates tracer NN does not have RETRO emissions
            IF ( .not. ( RETRO < 0d0 ) ) THEN

	       ! Apply time-of-day factor
	       RETRO = RETRO * TODX

               ! Convert from molec/cm2/s to kg/box/timestep in order
               ! to be in the proper units for EMISRR array
	       EMX(1) = RETRO * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN)	    

	    ENDIF  
         ENDIF

         !--------------------------------------------------------------
         ! Get C2H6 emissions from Yaping Xiao's inventory (mpayer, 3/22/12)
         !--------------------------------------------------------------

         ! If C2H6 is a defined tracer ...
         IF ( NN == IDTC2H6 ) THEN

               C2H6_ANTHRO = GET_C2H6_ANTHRO( I, J, NN )

               ! Apply time-of-day factor
	       C2H6_ANTHRO = C2H6_ANTHRO * TODX

               ! Convert from molC/cm2/s to kg/box/timestep in order
               ! to be in the proper units for EMISRR array
               EMX(1) = C2H6_ANTHRO * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN)

         ENDIF
  
         !--------------------------------------------------------------
         ! Get CO & hydrocarbons emissions from RCP inventory
         ! (cdh, 10/14/11)
         !--------------------------------------------------------------

         IF ( LRCP ) THEN

	    ! Get RCP emissions
            IF (NN==IDTCO) THEN
               ! Ship CO emissions are handled separately below
               RCP = GET_RCP_EMISSION( I, J, NN, 
     &              LAND=.TRUE., SHIP=.FALSE. )
            ELSE
               ! Land and ship emissions for all hydrocarbons
               RCP = GET_RCP_EMISSION( I, J, NN, 
     &              LAND=.TRUE., SHIP=.TRUE. )
            ENDIF

	    ! -1 means tracer NN does not have RCP emissions
            IF ( RCP >= 0d0 ) THEN

	       ! Apply time-of-day factor
	       RCP = RCP * TODX

               ! Convert from molec/cm2/s to kg/box/timestep in order
               ! to be in the proper units for EMISRR array
               EMX(1) = RCP * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN)	    

	    ENDIF  

         ENDIF

         !--------------------------------------------------------------
         ! Get CO & Hydrocarbons from EMEP inventory over Europe
         !--------------------------------------------------------------

         ! If we are using EMEP emissions ...
         IF ( LEMEP ) THEN

            ! If we are over the European region ...
            IF ( GET_EUROPE_MASK( I, J ) > 0d0 ) THEN
         
               ! Get EMEP emissions 
               EMEP = GET_EMEP_ANTHRO( I, J, NN )

               ! -1 indicates tracer NN does not have EMEP emissions
               IF ( .not. ( EMEP < 0d0 ) ) THEN

                  ! Apply time-of-day factor
                  EMEP   = EMEP * TODX

                  ! Convert from molec/cm2/s to kg/box/timestep in order
                  ! to be in the proper units for EMISRR array
                  EMX(1) = EMEP * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 

               ENDIF
            ENDIF
         ENDIF


         !--------------------------------------------------------------
         ! Get CO & Hydrocarbons from EPA/NEI inventory over the USA 
         !--------------------------------------------------------------

         ! If we are using EPA/NEI99 emissions ...
         IF ( LNEI99 ) THEN

            ! If we are over the USA ...
            IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
         
               ! Get EPA/NEI emissions (and apply time-of-day factor)
               EPA_NEI = GET_EPA_ANTHRO( I, J, NN, WEEKDAY )

               ! -1 indicates tracer NN does not have EPA/NEI emissions
               IF ( .not. ( EPA_NEI < 0d0 ) ) THEN

                  ! Apply time-of-day factor
                  EPA_NEI = EPA_NEI * TODX
         
                  ! Convert from molec/cm2/s to kg/box/timestep in order
                  ! to be in the proper units for EMISRR array
                  EMX(1)  = EPA_NEI * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 

               ENDIF
            ENDIF
         ENDIF

         !--------------------------------------------------------------
         ! Get CO from David Streets' inventory over Europe
         !--------------------------------------------------------------

         ! If we are using David Streets' emissions ...
         IF ( LSTREETS ) THEN

            ! If we are over the China region ...
            IF ( GET_SE_ASIA_MASK( I, J ) > 0d0 ) THEN
         
               ! Get STREETS emissions 
               STREETS = GET_STREETS_ANTHRO( I, J, NN, 
     &                                       MOLEC_CM2_S=.TRUE. )
         
               ! -1 indicates tracer NN does not have BRAVO emissions
               IF ( .not. ( STREETS < 0d0 ) ) THEN

                  ! Apply time-of-day factor
                  STREETS = STREETS * TODX

                  ! Convert from molec/cm2/s to kg/box/timestep in order
                  ! to be in the proper units for EMISRR array
                  EMX(1) = STREETS * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 

               ENDIF
            ENDIF
         ENDIF


         !--------------------------------------------------------------
         ! Get CO from BRAVO inventory over MEXICO
         !--------------------------------------------------------------

         ! If we are using BRAVO emissions ...
         IF ( LBRAVO ) THEN

            ! If we are over the Mexican region ...
            IF ( GET_BRAVO_MASK( I, J ) > 0d0 ) THEN
         
               ! Get BRAVO emissions 
               BRAVO = GET_BRAVO_ANTHRO( I, J, NN )
         
               ! -1 indicates tracer NN does not have BRAVO emissions
               IF ( .not. ( BRAVO < 0d0 ) ) THEN

                  ! Apply time-of-day factor
                  BRAVO  = BRAVO * TODX

                  ! Convert from molec/cm2/s to kg/box/timestep in order
                  ! to be in the proper units for EMISRR array.
                  ! Now, if on border, add to NEI99 emissions (phs, 5/7/08)
                  IF ( LNEI99 ) THEN 
                     IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN

                        EMX(1) = EMX(1) + 
     &                       BRAVO * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 
                     ENDIF
                  ELSE

                     EMX(1) = BRAVO * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN)

                  ENDIF

               ENDIF
            ENDIF
         ENDIF


         !--------------------------------------------------------------
         ! Get CAC other emissions over Canada
         !--------------------------------------------------------------

         ! If we are using CAC emissions ...
         IF ( LCAC ) THEN

            ! If we are over the China region ...
            IF ( GET_CANADA_MASK( I, J ) > 0d0 ) THEN

               ! Get CAC emissions 
               CAC = GET_CAC_ANTHRO( I, J, NN, MOLEC_CM2_S=.TRUE. )
         
               ! -1 indicates tracer NN does not have CAC emissions
               IF ( .not. ( CAC < 0d0 ) ) THEN

                  ! Apply time-of-day factor
                  CAC = CAC * TODX

                  IF ( LNEI99 ) THEN
                     ! If on border, add to NEI99 emissions (which contain
                     ! no Canadian component)
                     IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
                        EMX(1) = EMX(1) + CAC *
     &                           ( DTSRCE * AREA_CM2 ) / XNUMOL(NN)
                     ELSE
                        EMX(1) = CAC *
     &                           ( DTSRCE * AREA_CM2 ) / XNUMOL(NN)
                     ENDIF
                  ELSE

                    ! Else replace base emissions with CAC
                     EMX(1) = CAC *
     &                        ( DTSRCE * AREA_CM2 ) / XNUMOL(NN)
                  ENDIF

               ENDIF
            ENDIF
         ENDIF

         ! If we are using EPA/NEI2005 emissions ...
         IF ( LNEI05 ) THEN

            ! If we are over the USA ...
            IF ( NEI05_MASK( I, J ) > 0d0 ) THEN

               NEI05 = 0D0
               
               ! Loop over all of the emission levels
               ! For now lump levels together (phs, 10/20/09)
               DO LL = 1, NOXLEVELS 
               
                  ! Get EPA/NEI emissions               
                  EPA_NEI = GET_NEI2005_ANTHRO( I, J, LL, NN,
     &                 WEEKDAY, MOLEC_CM2_S=.TRUE. )

                  ! -1 indicates tracer NN does not have EPA/NEI emissions
                  IF ( EPA_NEI < 0d0 )  EXIT

                  NEI05 = NEI05 + EPA_NEI
                  
               ENDDO

               
               IF ( EPA_NEI > -1d0 ) THEN

                  ! Apply time-of-day factor
                  NEI05 = NEI05 * TODX
  
                  ! Convert from molec/cm2/s to kg/box/timestep in order
                  ! to be in the proper units for EMISRR array
                  EMX(1)  = NEI05 * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN)

               ENDIF
            ENDIF              
         ENDIF

         ! Account for CO production from anthropogenic VOC's
         ! -> For Tagged CO, enhance CO production by 18.5%
         ! -> For full-chem, enhance CO production by 2%
         ! (bnd, bmy, 4/26/01; jaf, mak, bmy, 2/14/08)
         ! Scaling factor is now correctly applied after 
         ! calculating emissions. (jaf, ccc, 2/25/09)
         ! Modifications of the scaling using Rynda GRL 2008.
         ! (jaf, ccc, 2/25/09)
         ! Added a nested if (phs, 7/9/09)
         IF ( ITS_A_TAGCO_SIM ) THEN           
            IF ( LICARTT ) THEN
               IF ( GET_USA_MASK(I,J) > 0.d0 ) THEN
                  IF ( NN == IDTCO ) EMX(1) = EMX(1) * 1.39d0
               ELSE
                  IF ( NN == IDTCO ) EMX(1) = EMX(1) * 1.19d0
               ENDIF
            ELSE
               IF ( NN == IDTCO ) EMX(1) = EMX(1) * 1.19d0
            ENDIF
         ELSE
            IF ( NN == IDTCO ) EMX(1) = EMX(1) * 1.02d0
         ENDIF


         !--------------------------------------------------------------
         ! Add ship emissions for CO (phs, 7/9/09)
         !--------------------------------------------------------------
         SHIP = 0D0
         
         IF ( NN == IDTCO ) THEN

            ! get global inventory first
            IF ( LEDGARSHIP ) THEN
               
               SHIP = GET_EDGAR_CO( I, J, MOLEC_CM2_S=.TRUE.,
     $                              SHIP=.TRUE.)
            
            ELSE IF ( LICOADSSHIP ) THEN

               SHIP = GET_ICOADS_SHIP( I, J, NN, MOLEC_CM2_S=.TRUE. )

            ELSE IF ( LRCPSHIP ) THEN

               SHIP = GET_RCP_EMISSION( I, J, NN, 
     &              LAND=.FALSE., SHIP=.TRUE. )

            ENDIF

            ! overwrite Europe
            IF ( LEMEPSHIP ) THEN
            
               ! If we are over the European region ...
               ! (bug fix, now overwrite only EMEP region, cdh 1/2/2013)
               IF ( GET_EUROPE_MASK( I, J ) > 0d0 ) THEN
                  SHIP =  GET_EMEP_ANTHRO( I, J, NN, SHIP=.TRUE.)
               ENDIF

            ENDIF

            ! Convert to same units as EMX(1), and add
            SHIP = SHIP * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 

            EMX(1) = EMX(1) + SHIP
            
         ENDIF
     
         !--------------------------------------------------------------
         ! Store in EMISRR array and archive diagnostics
         !--------------------------------------------------------------
         EMISRR(I,J,N) = EMISRR(I,J,N) + EMX(1) * XNUMOL(NN) / DTSRCE

         ! ND29 = CO source diagnostic... 
         ! store as [molec/cm2/s] in AD29(:,:,1)
         IF ( ND29 > 0 .and. NN == IDTCO ) THEN
            AD29(I,J,1) = AD29(I,J,1) +
     &           ( EMX(1) * XNUMOL(NN) / ( DTSRCE * AREA_CM2 ) ) 
         ENDIF

         ! ND36 = Anthro source diagnostic...store as [molec/cm2] 
         ! and convert to [molec/cm2/s] in DIAG3.F
         IF ( ND36 > 0 ) THEN
            AD36(I,J,N) = AD36(I,J,N) + 
     &           ( EMX(1) * XNUMOL(NN) / AREA_CM2 ) 
         ENDIF    

      ENDIF

      !=================================================================
      ! Restore EMISR, EMISRN to original values
      !=================================================================
      IF ( N == IDENO ) THEN
         EMISRN(IREF,JREF,1:NOXLEVELS) = XEMISRN(1:NOXLEVELS)
      ELSE
         EMISR(IREF,JREF,N) = XEMISR
      ENDIF  

      END SUBROUTINE EMFOSSIL
!EOC
