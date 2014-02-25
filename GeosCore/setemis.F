!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: setemis
!
! !DESCRIPTION: Subroutine SETEMIS places emissions computed from GEOS-Chem
!  subroutines into arrays for SMVGEAR II chemistry. 
!\\
!\\
!  SETEMIS converts from units of [molec tracer/box/s] to units of
!  [molec chemical species/cm3/s], and stores in the REMIS array.  For
!  hydrocarbons that are carried through the GEOS-CHEM model as [molec C], 
!  these are converted back to [molec hydrocarbon], and then stored in REMIS.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETEMIS( EMISRR,    EMISRRN,   am_I_Root,
     &                    Input_Opt, State_Met, RC         )
!
! !USES:
!
      USE BIOFUEL_MOD,        ONLY : BIOFUEL,   BFTRACE, NBFTRACE
      USE BIOMASS_MOD,        ONLY : BIOMASS,   BIOTRCE
      USE BIOMASS_MOD,        ONLY : BIOBGAS
      USE COMODE_MOD,         ONLY : JLOP,      REMIS,   VOLUME
      USE COMODE_MOD,         ONLY : IXSAVE,    IYSAVE,  IZSAVE
      USE DIAG_MOD,           ONLY : AD12
      USE EMISSIONS_MOD,      ONLY : NOx_SCALING
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE LIGHTNING_NOX_MOD,  ONLY : EMIS_LI_NOx
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_L
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE TRACERID_MOD,       ONLY : CTRMB,     IDEMIS,  IDENO
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT

      ! Ship plume emissions moved from calcrate.F (cdh, 4/24/2013)
      USE DRYDEP_MOD,         ONLY : SHIPO3DEP
      USE ICOADS_SHIP_MOD,    ONLY : GET_ICOADS_SHIP 
      USE EDGAR_MOD,          ONLY : GET_EDGAR_NOx
      USE EMEP_MOD,           ONLY : GET_EMEP_ANTHRO, GET_EUROPE_MASK
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TRACERID_MOD,       ONLY : IDTNO,  IDTO3,  IDTHNO3
      USE TRACERID_MOD,       ONLY : IDNO,   IDO3,   IDHNO3, IDNO2
      USE TRACERID_MOD,       ONLY : IDENO,  IDEO3,  IDEHNO3
      USE DIAG_MOD,           ONLY : AD32_SHIP, AD32_SHIP_COUNT 
      USE DIAG_MOD,           ONLY : AD36_SHIP, AD36_SHIP_COUNT
      USE DIAG_MOD,           only : AD63,      AD63_COUNT
      USE DIAG63_MOD,         ONLY : DO_SAVE_DIAG63
      USE PARANOX_MOD,        ONLY : INTERPOLATE_LUT2
      USE COMODE_MOD,         ONLY : CSPEC, AIRDENS
      USE RCP_MOD,            ONLY : GET_RCP_EMISSION, RCP_AC_NOx
      USE CMN_SIZE_MOD             ! Size parameters
      USE COMODE_LOOP_MOD          ! IDEMS, NEMIS
      USE CMN_DIAG_MOD             ! Diagnostic flags
      USE CMN_NOX_MOD              ! GEMISNOX2

      IMPLICIT NONE

!
! !INPUT PARAMETERS:
!
      ! CO, hydrocarbon emission   [molec tracer/box/s]
      REAL*8,         INTENT(IN)  :: EMISRR(IIPAR,JJPAR,NEMPARA+NEMPARB)

      ! Multi-level NOx emissions  [molec NOx/box/s]
      REAL*8,         INTENT(IN)  :: EMISRRN(IIPAR,JJPAR,NOXEXTENT)  

      ! Is this the root CPU?
      LOGICAL,        INTENT(IN)  :: am_I_Root

      ! Input Options object
      TYPE(OptInput), INTENT(IN)  :: Input_Opt

      ! Meteorology State object
      TYPE(MetState), INTENT(IN)  :: State_Met
!
! !OUTPUT PARAMETERS:
!
      ! Success or failure?
      INTEGER,        INTENT(OUT) :: RC 
!
! !REMARKS:
!  Developers: lwh, jyl, gmg, djj, bdf, bmy, 6/8/98, 6/11/08
!  (lwh, jyl, gmg, djj, bdf, bmy, 6/8/98, 6/11/08)
!                                                                             .
!  REMIS(JLOOP,N) = emis. rate of species corr. to tracer N in box JLOOP
!                   (reaction number NTEMIS(N))
!
! !REVISION HISTORY: 
!  (1 ) Original code from Harvard Tropospheric Chemistry Module for 3-D 
!        applications by Larry Horowitz, Jinyou Liang, Gerry Gardner, 
!        Prof. Daniel Jacob of Harvard University (Release V2.0)  
!  (2 ) New version 3.0 by Bob Yantosca to place NOx emissions into boxes  
!        above the surface. (bmy, 6/8/98)     
!  (3 ) Also now do chemistry up to the location of the annual mean         
!         tropopause (bmy, 12/9/99)                                         
!  (4 ) BURNEMIS is now dynamically allocatable and is contained in F90       
!        module "biomass_mod.f".  BIOTRCE and NBIOTRCE are also contained
!        in "biomass_mod.f".  (bmy, 9/12/00)                     
!  (5 ) BIOFUEL is now dynamically allocatable and is contained in F90
!        module "biofuel_mod.f".  BFTRACE and NBFTRACE are also contained
!        in "biofuel_mod.f" (bmy, 9/12/00, 4/17/01)
!  (6 ) BURNEMIS and BIOFUEL are now treated as true global arrays,  
!        and need to be referenced by the global offset variables          
!        IREF = I + I0 and JREF = J + J0 (bmy, 9/12/00)                    
!  (7 ) Now reference JLOP, REMIS, VOLUME from F90 module "comode_mod.f",  
!        in order to save memory (bmy, 10/19/00)                           
!  (8 ) Now add in up to NBFTRACE biofuel species (bmy, 4/17/01) 
!  (9 ) Add new subroutine header, updated comments, cosmetic changes.
!        (bmy, 4/17/01)  
!  (10) Updated comments -- GEMISNOX is [molec/cm3/s]. (bdf, bmy, 6/7/01)     
!  (11) For GEOS-3, we now distributes surface emissions throughout the 
!        boundary layer.  This is necessary since the first couple of GEOS-3 
!        surface layers are very thin.  Piling up of emissions into a small 
!        layer will cause SMVGEAR to choke.  (bdf, bmy, 6/15/01)
!  (12) Also now reference BFTRACE and NBFTRACE from "biofuel_mod.f", 
!        and reference AD12 from "diag_mod.f". (bdf, bmy, 6/15/01)
!  (13) For GEOS-1, GEOS-STRAT, emit into the surface layer, as we did
!        in prior versions. (bmy, 6/26/01)
!  (14) Bug fix: corrected a typo for the biofuel emissions (bmy, 7/10/01)
!  (15) Bug fix: make sure BIOMASS and BIOFUEL, and SOIL NOx emissions have 
!        units of [molec/box/s] before distributing thru the boundary layer.  
!        This involves multiplication by VOLUME(JLOOP1) and division by
!        VOLUME(JLOOP). (bmy, 7/16/01)
!  (16) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  BIOFUEL(:,IREF,JREF) is now
!        BIOFUEL(:,I,J). BURNEMIS(:,IREF,JREF) is now BURNEMIS(:,I,J).
!        Replace PW(I,J) with P(I,J). (bmy, 9/28/01)
!  (17) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (18) Now references GET_PEDGE from "pressure_mod.f", to compute P at 
!        the bottom edge of grid box (I,J,L).  (dsa, bdf, bmy, 8/21/02)
!  (19) Now reference IDTNOX, IDENOX, etc from "tracerid_mod.f" (bmy, 11/6/02)
!  (20) Remove references to IREF, JREF (bmy, 2/11/03)
!  (21) NEMIS is now NEMIS(NCS) for SMVGEAR II (gcc, bdf, bmy, 4/1/03)
!  (22) Added parallel loop over N.  Also directly substituted JLOP(I,J,1) 
!        for all instances of JLOOP1.  Updated comments. (hamid, bmy, 3/19/04)
!  (23) Bug fix for COMPAQ compiler...do not use EXIT from w/in parallel loop.
!        (auvray, bmy, 11/29/04)
!  (24) Now replace XTRA2 with GET_PBL_TOP_L in "pbl_mix_mod.f".  Now remove
!        reference to CMN, it's obsolete.  Now references GET_TPAUSE_LEVEL
!        from "tropopause_mod.f" (bmy, 8/22/05)
!  (25) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (26) Now updated for new "biomass_mod.f" (bmy, 4/5/06)
!  (27) Now account for the different definition of tropopause in case 
!        of variable tropopause.   The BIOMASS array from "biomass_mod.f" is 
!        now in units of [molec CO/cm2/s].  Adjust unit conversion accordingly.
!        Also replace NBIOMAX with NBIOMAX_GAS, since aerosol biomass is
!        handled elsewhere.  (bdf, phs, bmy, 9/27/06)
!  (28) Now replace GEMISNOX array (from CMN_NOX) with module arrays
!        EMIS_LI_NOx and EMIS_AC_NOx (ltm, bmy, 10/3/07)
!  (29) Bug fix: resize EMISRR to be consistent w/ CMN_O3 (bmy, jaf, 6/11/08)
!  (30) Limit emissions into the surface level only (lin, 5/29/09)
!  (31) Bug fix: cycle if IDEMIS(NN) <= 0 to avoid array-out-of-bounds
!        errors (bmy, 8/6/09)
!  (32) Check for emissions above PBL -anthro NOx only for now- (phs, 10/27/09)
!  (33) Modify selection of biomass burning emissions (hotp, 8/3/09)
!  (34) Moved NOx scaling to improve parallelization. (ccc, 11/10/10)
!  16 Dec 2010 - R. Yantosca - Removed obsolete, commented-out code
!  16 Dec 2010 - R. Yantosca - Added ProTeX headers
!  21 Dec 2010 - R. Yantosca - Now set REMIS=0d0.  Also updated comments.
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now reference IXSAVE, IZSAVE from comode_mod.F
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  02 Aug 2012 - R. Yantosca - Add error trap for JLOOP=0 under DEVEL tag
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LPRT and Input_Opt%LNLPBL
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  17 Jun 2013 - R. Yantosca - Bug fix: declare FRACTION_NOX, INT_OPE with
!                              REAL*4 explicitly to avoid numerical errors
!                              when OMP=yes.
!  18 Jun 2013 - R. Yantosca - Bug fix: always add into REMIS to better ensure
!                              that identical results will be obtained when
!                              compiling with OMP=no vs. OMP=yes.
!  18 Jun 2013 - R. Yantosca - Remove COEF11 variable, that is always =1 for
!                              O3, HNO3, and NO.
!  22 Jul 2013 - M. Sulprizio- Now copy LRCPSHIP from Input_Opt
!  31 Jul 2013 - M. Sulprizio- Now only add aircraft emissions to REMIS if using
!                              RCP aircraft emissions. Emissions are added
!                              directly to STT in aeic_mod.F when using AEIC
!                              aircraft emissions.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!      
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL          :: IS_LI_NOx, IS_AC_NOx
      INTEGER          :: I, J,  JLOOP, JLOOP1, LTROP
      INTEGER          :: L, LL, N, NN,  NBB, NBF, TOP, TOPMIX
      REAL*8           :: COEF1,   TOTPRES, DELTPRES
      REAL*8           :: EMIS_BL, NOXTOT,  TOTAL, A_CM2

      ! Added by G. Vinken, C. Holmes
      REAL*8           :: SHIP,   TOTO3, JNO2, JO1D
      REAL*8           :: DTSRCE, AREA_CM2
!      REAL*4           :: F_NOX,        I_OPE 
      REAL*4           :: FRACTION_NOX, INT_OPE 
      INTEGER          :: NK
      CHARACTER*8      :: SPECNAME
      REAL*8, EXTERNAL :: FJFUNC
      
      ! For fields from Input_Opt
      LOGICAL          :: LDRYD,       LNLPBL,     LPRT
      LOGICAL          :: LICOADSSHIP, LEDGARSHIP, LEMEPSHIP     
      LOGICAL          :: LRCPSHIP,    LRCPAIR

      !%%% BMY KLUDGE to only print out warning message about 
      !%%% non-emitted tracers on the first timestep (bmy, 5/28/13)
      LOGICAL, SAVE    :: FIRST = .TRUE.

      !=================================================================
      ! SETEMIS begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS
      
      ! some ajdustments for non-local PBL (Lin, 03/31/09)
      ! NOTE: If NCS=0 and we are NOT doing a full-chemistry run, then
      ! setting REMIS=0 can lead to a seg-fault error.  Therefore, you should
      ! add an error check in the calling routine to prevent SETEMIS from
      ! being called if we are not doing full-chemistry. (bmy, 12/21/10)
      IF ( NCS == 0 ) THEN
         REMIS = 0d0
         RETURN
      ENDIF

      ! Test if the EMIS_LI_NOx array is allocated
      IS_LI_NOx   = ALLOCATED( EMIS_LI_NOx )
      
      ! Now use logical fields from Input_Opt
      LDRYD       = Input_Opt%LDRYD
      LNLPBL      = Input_Opt%LNLPBL
      LPRT        = Input_Opt%LPRT
      LICOADSSHIP = Input_Opt%LICOADSSHIP
      LEDGARSHIP  = Input_Opt%LEDGARSHIP 
      LEMEPSHIP   = Input_Opt%LEMEPSHIP
      LRCPSHIP    = Input_Opt%LRCPSHIP
      LRCPAIR     = Input_Opt%LRCPAIR

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      IF ( LPRT .and. am_I_Root ) THEN
         WRITE(6,*) 'SCALE NOx emissions by', NOx_SCALING
      ENDIF

      ! Initialize REMIS(*,N) -- the emission rate array
      REMIS = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( N,       NN,     NBB,      NBF,          I       )
!$OMP+PRIVATE( J,       L,      JLOOP,    COEF1,        TOP     )
!$OMP+PRIVATE( TOTPRES, NOXTOT, DELTPRES, EMIS_BL,      A_CM2   )
!$OMP+PRIVATE( TOTO3,   SHIP,   AREA_CM2, FRACTION_NOX, INT_OPE )
!$OMP+PRIVATE( NK,      JNO2,   JO1D,     SPECNAME,     TOPMIX  )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over emission species
      DO N = 1, NEMIS(NCS)

         ! Get CTM tracer number NN corresponding to emission species N
         NN = IDEMS(N)
         IF ( NN == 0 ) CYCLE

         ! We have to search for the biomass burning species in 
         ! BIOTRCE with the same CTM tracer number NN as in IDEMS
         NBB = 0
         IF ( ALLOCATED( BIOMASS ) ) THEN 
            ! Now loop over NBIOMAX and select gas phase species (hotp 8/3/09)
            !DO I = 1, NBIOMAX_GAS
            DO I = 1, NBIOMAX
               ! Only select gas phase chemistry species (ignore
               ! aerosols) (hotp 8/3/09)
               IF ( BIOBGAS(I) ) THEN

                  IF ( BIOTRCE(I) == NN ) THEN 
                     NBB = I
                     EXIT
                  ENDIF
               ENDIF ! end BIOBGAS (hotp 8/3/09)
            ENDDO
         ENDIF

         ! We have to search for the biofuel burning species in 
         ! BFTRACE with the same CTM tracer number NN as in IDEMS
         NBF = 0
         IF ( ALLOCATED( BIOFUEL ) ) THEN
            DO I = 1, NBFTRACE
               IF ( BFTRACE(I) == NN ) THEN
                  NBF = I
                  EXIT
              ENDIF
            ENDDO
         ENDIF            

         ! Skip to next tracer if IDEMIS(NN) is not defined in 
         ! order to avoid array-out-of-bounds errors (bmy, 8/6/04)
         IF ( IDEMIS(NN) <= 0 ) THEN
            IF ( FIRST .and. am_I_Root ) THEN
               PRINT*, 'Tracer ', NN, ' is not an emitted species!'
            ENDIF
            CYCLE
         ENDIF

         ! COEF1 = molecules of emission species / molecules of tracer
         COEF1 = 1.0 + CTRMB(NN, IDEMIS(NN))         

         ! Loop over Lat and Long boxes
         DO J = 1, NLAT
         DO I = 1, NLONG

            !===========================================================
            ! For GEOS-3: distribute surface emissions throughout the
            ! entire boundary layer.  Define some variables here.
            ! (bdf, 6/15/01)
            !===========================================================

            ! Top level of the boundary layer
            ! guard for b.l. being in first level.
            TOP = FLOOR( GET_PBL_TOP_L( I, J ) )
            IF ( TOP == 0 ) TOP = 1

            ! Check for possible emissions above PBL (phs, 27/10/09)
            TOPMIX = MIN( TOP, NOXEXTENT )
            
            ! Add option for non-local PBL (Lin, 03/31/09)
            IF (LNLPBL) TOP = 1

            ! Pressure thickness of entire boundary layer [hPa]
            TOTPRES = GET_PEDGE(I,J,1) - GET_PEDGE(I,J,TOP+1)

            ! For NOx only....
            ! Replaced NOx with NO (mpayer, 3/14/13)
            IF ( N == IDENO ) THEN

               !========================================================
               ! Anthropogenic NOx emissions [molec/box/s]
               ! Distribute emissions thru the entire boundary layer
               !========================================================

               ! Sum anthro NOx emissions over all levels [molec NOx/box/s]
               NOXTOT = 0d0
               DO L = 1, TOPMIX
                  NOXTOT = NOXTOT + EMISRRN(I,J,L)
               ENDDO
               
               ! Loop over the boundary layer
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0
                  
                  IF ( JLOOP /= 0 ) THEN
                     
                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
                     
                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES
                     
                     ! Of the total anthro NOx, the fraction DELTPRES/TOTPRES
                     ! goes into level L, since that is the fraction of the
                     ! boundary layer occupied by level L.  Also divide NOx 
                     ! by COEF1 to convert from [molec NOx/box/s] to 
                     ! [molec NO/box/s], which is really what gets emitted.
                     EMIS_BL        = ( NOXTOT   / COEF1   ) *
     &                    ( DELTPRES / TOTPRES )
                     
                     ! Convert anthro NOx emissions from [molec NO/box/s]
                     ! to [molec NO/cm3/s] and store in the REMIS array
                     REMIS(JLOOP,N) = EMIS_BL / VOLUME(JLOOP)

                  ENDIF
               ENDDO

               ! put any emissions above PBL into corresponding
               ! REMIS (phs, 27/10/09)
               IF ( NOXEXTENT > TOPMIX ) THEN
                  DO L = TOPMIX+1, NOXEXTENT
                     JLOOP          = JLOP(I,J,L)
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
                     ! Add this error trap to prevent out of bounds error
                     ! but we should benchmark first before adding to
                     ! the std G-C code (bmy, 8/2/12)
                     IF ( JLOOP == 0 ) CYCLE
#endif
                     REMIS(JLOOP,N) = EMISRRN(I,J,L) /
     $                    ( VOLUME(JLOOP) * COEF1 )
                  ENDDO
               ENDIF   
                     
               
               !========================================================
               ! Soil Nox emissions [molec/cm3/s] 
               ! Distribute emissions thru the entire boundary layer
               !========================================================
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
                     
                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES

                     ! Soil NOx is in [molec/cm3/s], so we need to multiply
                     ! by VOLUME(JLOP(I,J,1)) to convert it to [molec/box/s],
                     ! VOLUME(JLOP(I,J,1)) is the volume in cm3 of the surface
                     ! grid box (I,J,1).  Then we need to divide that by 
                     ! COEF1 to convert from [molec NOx/box/s] to 
                     ! [molec NO/box/s], which is really what gets emitted.  
                     ! Of the total soil NOx, the fraction DELTPRES/TOTPRES 
                     ! goes into level L, since that is the fraction of the 
                     ! boundary layer occupied by level L.  Store in EMIS_BL.
                     EMIS_BL    = ( GEMISNOX2(I,J)                    *
     &                              VOLUME( JLOP(I,J,1) ) / COEF1   ) * 
     &                              ( DELTPRES            / TOTPRES )

                     ! Since EMIS_BL is in [molec NO/box/s], we have to
                     ! divide by VOLUME(JLOOP), which is the volume of the
                     ! grid box (I,J,L) to convert back to [molec/cm3/s].
                     ! Store in the REMIS array for SMVGEAR.
                     REMIS(JLOOP,N) = REMIS(JLOOP,N) + 
     &                                ( EMIS_BL / VOLUME(JLOOP) ) 

                  ENDIF
               ENDDO

               !========================================================
               ! Aircraft and Lightning NOx [molec/cm3/s]
               ! Distribute emissions in the troposphere
               !========================================================
               DO L = 1, LLTROP

                  ! Now use ITS_IN_THE_STRAT to test if we have passed from
                  ! the trop into the strat.  We only want to partition
                  ! emissions in the troposphere. (bmy, 9/10/10)
                  IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) EXIT
                  
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     !-----------------
                     ! Lightning NOx
                     !-----------------
                     IF ( IS_LI_NOx ) THEN
 
                        ! Divide lightning NOx by COEF1 to convert
                        ! from [molec NOx/cm3/s] to [molec NO/cm3/s], since
                        ! NO is the actual emission species for NOx.
                        EMIS_BL        = EMIS_LI_NOx(I,J,L) / COEF1 

                        ! Save lightning NOx [molec NO/cm3/s] in REMIS
                        REMIS(JLOOP,N) = REMIS(JLOOP,N) + EMIS_BL

                     ENDIF

                     !-----------------
                     ! Aircraft NOx
                     !-----------------

                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     !% NOTE: Aircraft emissions from the AEIC     %
                     !% inventory are added directly into the STT  %
                     !% array in aeic_mod.F (mpayer, 7/31/13)      %
                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                     ! If using aircraft emissions from RCP scenarios
                     IF ( LRCPAIR ) THEN

                        ! Divide aircraft NOx by COEF1 to convert
                        ! from [molec NOx/cm3/s] to [molec NO/cm3/s], since
                        ! NO is the actual emission species for NOx.
                        EMIS_BL        = RCP_AC_NOx(I,J,L) / COEF1 

                        ! Save aircraft NOx [molec NO/cm3/s] in REMIS
                        REMIS(JLOOP,N) = REMIS(JLOOP,N) + EMIS_BL

                     ENDIF

                     ! NOx scaling moved here for optimisation. (ccc, 11/10/10)
                     ! Add possible scaling of NOx emissions
                     REMIS(JLOOP,N) = NOx_SCALING * REMIS(JLOOP,N)             

                  ENDIF
               ENDDO

               !========================================================
               ! Ship NOx (emitted as NO, O3, HNO3 after plume evolution)
               ! [molec/cm3/s]
               ! Distribute emissions thru the entire boundary layer
               !========================================================
                
               ! Update only if ship emissions are turned on
               IF ( LICOADSSHIP .or. LEDGARSHIP .or. LEMEPSHIP .or. 
     &              LRCPSHIP ) THEN

                  ! Surface area of grid box
                  AREA_CM2 = GET_AREA_CM2( I, J, 1 )
   
                  ! Reset
                  SHIP = 0D0
                  
                  ! handle global inventory first
                  IF ( LEDGARSHIP ) THEN 

                      ! Get SHIP EDGAR emissions for NOx [molec/cm2/s]
                     SHIP = GET_EDGAR_NOx( I, J, 
     &                      MOLEC_CM2_S=.TRUE., SHIP=.TRUE.)
                         
                  
                  ! ICOADS ship emissions (cklee,7/09/09)
                  ELSE IF ( LICOADSSHIP ) THEN

                     ! Get ICOADS  emissions for NOx [molec/cm2/s]
                     SHIP = GET_ICOADS_SHIP( I, J, IDTNO, 
     &                      MOLEC_CM2_S=.TRUE. )
                  
                  ! RCP ship emissions (cdh, 10/14/11)
                  ELSE IF ( LRCPSHIP ) THEN

                     ! Get SHIP RCP emissions for NOx [molec/cm2/s]
                     SHIP = GET_RCP_EMISSION( I, J, IDTNO,  
     &                      LAND=.FALSE., SHIP=.TRUE.)

                  ENDIF

                  ! Overwrite Europe
                  IF ( LEMEPSHIP ) THEN
                 
                     ! Prevent overwriting ICOADS data with a 0 value from EMEP 
                     IF (GET_EUROPE_MASK( I, J ) > 0d0) THEN
           
                        ! Get SHIP EMEP emissions for NOx [molec/cm2/s]
                        SHIP = GET_EMEP_ANTHRO( I, J, IDTNO, 
     &                         SHIP=.TRUE. )

                     ENDIF

                  ENDIF
                  
                  ! Add possible scaling of NOx emissions
                  SHIP = SHIP * NOx_SCALING
                  
                  !----------------------------------
                  ! Get J-Values for J(NO2) and J(O3)
                  !----------------------------------

                  ! Check if sun is up
                  IF ( State_Met%SUNCOSmid(I,J) > 0d0 ) THEN

                     ! Loop over photolysis reactions to find NO2, O3
                     DO L = 1, JPHOTRAT(NCS)
                        
                        ! Reaction number
                        NK  = NRATES(NCS) + L
                        
                        ! Name of species being photolyzed
                        SPECNAME = NAMEGAS(IRM(1,NK,NCS))

                        ! Check if this is NO2 or O3, store values, 1/s
                        SELECT CASE ( TRIM( SPECNAME ) )
                        CASE ( 'NO2' )
                           JNO2 = FJFUNC(I,J,1,L,1,SPECNAME)
                        CASE ( 'O3' )
                           JO1D = FJFUNC(I,J,1,L,1,SPECNAME)
                        CASE DEFAULT
                        END SELECT
                        
                     ENDDO

                  ELSE

                     ! J-values are zero when sun is down
                     JNO2 = 0d0
                     JO1D = 0d0

                  ENDIF

                  ! Determine fraction of NOx remaining and integrated Ozone
                  ! Production Efficiency for ship emiss (gvinken,mpayer,2/7/12)
                  ! Uses surface-layer concentrations of O3, NO, NO2 (molec/cm3)
                  ! from CSPEC and air density (molec/cm3)
                  JLOOP   = JLOP(I,J,1)
                  CALL INTERPOLATE_LUT2( I, J,  
     &                 CSPEC(JLOOP,IDO3),  CSPEC(JLOOP,IDNO),
     &                 CSPEC(JLOOP,IDNO2), AIRDENS(JLOOP),
     &                 JO1D,               JNO2, 
     &                 FRACTION_NOx,       INT_OPE,
     &                 State_Met           )

                  !-----------------
                  ! Ship NO
                  !-----------------
 
                  ! Loop over the boundary layer
                  DO L = 1, TOP
                     JLOOP   = JLOP(I,J,L)
                     EMIS_BL = 0d0

                     IF ( JLOOP == 0 ) CYCLE
                     
                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
                     
                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES

                     ! Of the total ship NOx, the fraction FRACTION_NOX
                     ! survives after plume dilution and chemistry. 
                     ! The fraction DELTPRES/TOTPRES goes into level L,
                     ! since that is the fraction of the boundary layer
                     ! mass occupied by level L. 
                     ! Final units: molec(NO)/cm2/s
                     EMIS_BL = ( SHIP     * FRACTION_NOX )
     &                       * ( DELTPRES / TOTPRES      )

                     ! Convert anthro NOx emissions from [molec NO/cm2/s]
                     ! to [molec NO/cm3/s] and store in the REMIS array
                     REMIS(JLOOP,IDENO) = REMIS(JLOOP,IDENO) 
     &                                  + ( EMIS_BL * AREA_CM2 ) 
     &                                  / VOLUME(JLOOP)

                  ENDDO

                  !-----------------
                  ! Ship HNO3
                  !-----------------

                  ! Loop over the boundary layer
                  DO L = 1, TOP
                     JLOOP   = JLOP(I,J,L)
                     EMIS_BL = 0d0
                  
                     IF ( JLOOP == 0 ) CYCLE
                     
                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
                     
                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES
                     
                     ! Of the total ship NOx, the fraction (1d0-FRACTION_NOX)
                     ! is converted to HNO3 during plume dilution and chemistry. 
                     ! The fraction DELTPRES/TOTPRES goes into level L,
                     ! since that is the fraction of the boundary layer
                     ! mass occupied by level L. Also, divide by COEF11 
                     ! to convert from [molec tracer/cm2/s] to 
                     ! [molec chem species/cm2/s].
                     ! Final units: molec(HNO3)/cm2/s
                     EMIS_BL = SHIP * ( 1D0      - FRACTION_NOX ) * 
     &                                ( DELTPRES / TOTPRES      )

                     ! Convert ship HNO3 emissions from [molec HNO3/cm2/s]
                     ! to [molec HNO3/cm3/s] and store in the REMIS array
                     REMIS(JLOOP,IDEHNO3) = REMIS(JLOOP,IDEHNO3) 
     &                                    + ( EMIS_BL * AREA_CM2 ) 
     &                                    / VOLUME(JLOOP)
                        
                  ENDDO

                  !-----------------
                  ! Ship O3
                  !-----------------

                  ! Total emissions over the boundary layer
                  ! Of the total ship NOx, the fraction 
                  ! (1d0-FRACTION_NOX)*INT_OPE is converted to O3 during 
                  ! plume dilution and chemistry. 
                  ! Final units: molec(O3)/cm2/s
                  EMIS_BL = SHIP * (1D0 - FRACTION_NOX) * INT_OPE

                  ! Ship plume chemistry can create or destroy net O3
                  ! Treat positive O3 production as emissions 
                  ! Treat O3 *destruction* in plume as dry deposition
                  ! (cdh, 3/21/2013)
                  IF (EMIS_BL >= 0d0) THEN

                     ! No need for deposition when net O3 production from 
                     ! ships is positive
!$OMP CRITICAL
                     IF (LDRYD) SHIPO3DEP(I,J) = 0d0
!$OMP END CRITICAL

                     ! Loop over the boundary layer
                     DO L = 1, TOP
                        JLOOP   = JLOP(I,J,L)
                  
                        IF ( JLOOP == 0 ) CYCLE
                        
                        ! Thickness of level L [mb]
                        DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
                     
                        ! Add option for non-local PBL (Lin, 03/31/09)
                        IF (LNLPBL) DELTPRES = TOTPRES

                        ! The fraction DELTPRES/TOTPRES goes into level L,
                        ! since that is the fraction of the boundary layer
                        ! mass occupied by level L. 
                        ! The ratio ( AREA_CM2 / VOLUME(JLOOP) ) converts
                        ! molec/cm2 to molec/cm3
                        ! Final units: molec(O3)/cm3/s
                        REMIS(JLOOP,IDEO3) = REMIS(JLOOP,IDEO3) 
     &                                     + EMIS_BL 
     &                                     * ( DELTPRES / TOTPRES ) 
     &                                     * AREA_CM2 / VOLUME(JLOOP) 
         
                     ENDDO
                     
                  ELSE

                     ! No change in REMIS

                     ! Initialize
                     TOTO3 = 0d0

                     ! Loop over the boundary layer
                     ! If LNLPBL, then TOP=1
                     DO L = 1, TOP
                        JLOOP   = JLOP(I,J,L)

                        IF ( JLOOP == 0 ) CYCLE

                        ! Cumulative O3 in the PBL, molec(O3)/cm2
                        TOTO3 = TOTO3 + ( State_Met%BXHEIGHT(I,J,L) *
     &                          1d2 ) * CSPEC(JLOOP,IDO3)
                        
                     ENDDO

                     ! Effective O3 deposition frequency
                     ! EMIS_BL (negative) is total O3 destruction in the 
                     ! column, molec(O3)/cm2/s
                     ! TOTO3 is the total O3 content of the PBL, molec(O3)/cm2
                     ! Final units: 1/s
!$OMP CRITICAL
                     IF (LDRYD) SHIPO3DEP(I,J) = ABS( EMIS_BL / TOTO3 )
!$OMP END CRITICAL

                  ENDIF

                  !-----------------
                  ! Ship diagnostics
                  !-----------------

                  !ND36 = Anthro source diagnostic [molec/cm2/s]
                  ! We need separate arrays for ship emissions because 
                  ! other emissions are calculated at different model times
                  IF ( ND36 > 0 ) THEN

                     AD36_SHIP(I,J,IDEHNO3) = AD36_SHIP(I,J,IDEHNO3) + 
     &                    SHIP * (1D0 - FRACTION_NOx)
                     
                     AD36_SHIP(I,J,IDEO3)  = AD36_SHIP(I,J,IDEO3) + 
     &                    SHIP * (1D0 - FRACTION_NOx) * INT_OPE
    
                     AD36_SHIP(I,J,IDENO) = AD36_SHIP(I,J,IDENO) + 
     &                    SHIP * FRACTION_NOx

                     ! Increment counter, only once per time step
                     IF ( (I == 1) .and. (J==1) ) THEN
                        AD36_SHIP_COUNT = AD36_SHIP_COUNT + 1
                     ENDIF

                  ENDIF

                  ! ND32 = save anthro NOx for levels L=1,NOXEXTENT
                  ! [molec/cm2/s]
                  IF ( ND32 > 0 ) THEN

                     AD32_ship(I,J) = AD32_ship(I,J) + 
     &                    ( SHIP * FRACTION_NOx )

                     ! Increment counter, only once per time step
                     IF ( (I == 1) .and. (J==1) ) THEN
                        AD32_SHIP_COUNT = AD32_SHIP_COUNT + 1
                     ENDIF

                  ENDIF

                  ! ND63 = save the fraction of NOx remaining, integrated OPE
                  ! and total emissions of NOx and O3
                  IF ( DO_SAVE_DIAG63 ) THEN

                     ! Surviving fraction of emitted NOx, mol/mol
                     AD63(I,J,1) = AD63(I,J,1) + FRACTION_NOx

                     ! Integraed OPE, mol/mol
                     AD63(I,J,2) = AD63(I,J,2) + INT_OPE

                     ! Surviving NOx, molec/cm2 
                     ! (converted to molec/cm2/s in diag63_mod)
                     AD63(I,J,3) = AD63(I,J,3) + 
     &                    SHIP * FRACTION_NOx 

                     ! Ozone produced, molec/cm2
                     ! (converted to molec/cm2/s in diag63_mod)
                     AD63(I,J,4) = AD63(I,J,4) + 
     &                    SHIP * (1D0 - FRACTION_NOx) * INT_OPE

                     ! NOx emitted at ship stack, molec/cm2
                     ! (converted to molec/cm2/s in diag63_mod)
                     AD63(I,J,5) = AD63(I,J,5) + SHIP 

                     ! Increment counter, only once per time step
                     IF ( (I == 1) .and. (J==1) ) THEN
                        AD63_COUNT = AD63_COUNT + 1d0
                     ENDIF

                  ENDIF

               ENDIF
              
            ELSE

               !========================================================
               ! Anthropogenic tracers other than NOx [molec/box/s]
               ! Distribute emissions thru the entire boundary layer
               !========================================================
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN 

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES

                     ! Of the total tracer, the fraction DELTPRES/TOTPRES
                     ! goes into level L, since that is the fraction of the
                     ! boundary layer occupied by level L.  Also divide the
                     ! tracer by COEF1 to convert from [molec tracer/box/s] 
                     ! to [molec species/box/s].  For example, ISOPRENE is
                     ! carried by GEOS-CHEM as 5 carbons, so you would divide
                     ! by 5 to get "effective molecules" of ISOPRENE.
                     EMIS_BL        = ( EMISRR(I,J,N) / COEF1   ) *
     &                                ( DELTPRES      / TOTPRES )

                     ! Convert emissions from [molec species/box/s] to 
                     ! [molec species/cm3/s] and store in the REMIS array
                     REMIS(JLOOP,N) = REMIS(JLOOP,N) 
     &                              + ( EMIS_BL / VOLUME(JLOOP) )

                  ENDIF
               ENDDO
            ENDIF

            !===========================================================
            ! Add biomass burning source [molec/cm3/s]
            ! Distribute emissions thru the entire boundary layer
            !===========================================================
            IF ( NBB /= 0 ) THEN
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES
 
                     ! Grid box area [cm2]
                     A_CM2    = GET_AREA_CM2( IXSAVE(JLOOP),
     &                                        IYSAVE(JLOOP),
     &                                        IZSAVE(JLOOP)  )
               
                     ! Biomass burning is in [molec/cm2/s], so we need to 
                     ! multiply by A_CM2 to convert it to [molec/box/s].
                     ! Then we need to divide that by COEF1 to convert from 
                     ! [molec tracer/box/s] to [molec species/box/s].
                     ! Of the total biomass burning emissions, the fraction 
                     ! DELTPRES/TOTPRES goes into level L, since that is the 
                     ! fraction of the boundary layer occupied by level L.  
                     ! Store in EMIS_BL.
                     EMIS_BL  = ( BIOMASS(I,J,NBB) * A_CM2 / COEF1   ) *
     &                          ( DELTPRES                 / TOTPRES )

                     ! Since EMIS_BL is in [molec species/box/s], we 
                     ! have to divide by VOLUME(JLOOP), which is the 
                     ! volume of the grid box (I,J,L) to convert back to 
                     ! [molec species/cm3/s].  Store in the REMIS array.
                     REMIS(JLOOP,N) = REMIS(JLOOP,N) + 
     &                                ( EMIS_BL / VOLUME(JLOOP) )
                  ENDIF
               ENDDO
            ENDIF

            !===========================================================
            ! Add biofuel burning source [molec/cm3/s]
            ! Distribute emissions thru the entire boundary layer 
            !===========================================================
            IF ( NBF /= 0 ) THEN
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES

                     ! Biofuel burning is in [molec/cm3/s], so we need to 
                     ! multiply by VOLUME(JLOP(I,J,1)) to convert it to 
                     ! [molec/box/s], VOLUME(JLOP(I,J,1)) is the volume in cm3 
                     ! of the surface grid box (I,J,1).  Then we need to 
                     ! divide that by COEF1 to convert from 
                     ! [molec tracer/box/s] to [molec species/box/s].
                     ! Of the total biofuel burning emissions, the fraction 
                     ! DELTPRES/TOTPRES goes into level L, since that is the 
                     ! fraction of the boundary layer occupied by level L.  
                     ! Store in EMIS_BL.
                     EMIS_BL  = ( BIOFUEL(NBF,I,J) *
     &                            VOLUME( JLOP(I,J,1) ) / COEF1 ) *
     &                            ( DELTPRES / TOTPRES )

                     ! Since EMIS_BL is in [molec species/box/s], we 
                     ! have to divide by VOLUME(JLOOP), which is the 
                     ! volume of the grid box (I,J,L) to convert back to 
                     ! [molec species/cm3/s].  Store in the REMIS array.
                     REMIS(JLOOP,N) = REMIS(JLOOP,N) + 
     &                                ( EMIS_BL / VOLUME(JLOOP) )
                  ENDIF
               ENDDO
            ENDIF

            !===========================================================
            ! ND12 Diagnostic: Save the fraction of the boundary layer
            ! occupied by level L into the AD12 diagnostic array.
            !===========================================================
            IF ( N == 1 .and. ND12 > 0 ) THEN
               DO L = 1, MIN( FLOOR( GET_PBL_TOP_L( I, J ) ), LD12 )

                  ! Thickness of layer L [mb]
                  DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
 
                 ! Save boundary layer fraction into AD12
                  AD12(I,J,L) = AD12(I,J,L) + ( DELTPRES / TOTPRES )
               ENDDO
            ENDIF

         ENDDO  ! I
         ENDDO  ! J
      ENDDO     ! N
!$OMP END PARALLEL DO

      !%%% BMY KLUDGE to only print out warning message about 
      !%%% non-emitted tracers on the first timestep (bmy, 5/28/13)
      FIRST = .FALSE.

      END SUBROUTINE SETEMIS
!EOC
