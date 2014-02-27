!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: aerosol_mod
!
! !DESCRIPTION: Module AEROSOL\_MOD contains variables and routines for
!  computing optical properties for aerosols which are needed for both the
!  FAST-J photolysis and ND21 optical depth diagnostics. (bmy, 7/20/04,
!  2/10/09)
!\\
!\\
! !INTERFACE: 
!
      MODULE AEROSOL_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: AEROSOL_CONC
      PUBLIC :: AEROSOL_RURALBOX
      PUBLIC :: CLEANUP_AEROSOL
      PUBLIC :: INIT_AEROSOL
      PUBLIC :: RDAER
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: SOILDUST
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Pye, H.O.T., and J.H. Seinfeld, "A global perspective on aerosol from
!        low-volatility organic compounds", Atmos. Chem. & Phys., Vol 10, pp
!        4377-4401, 2010.
!
! !REVISION HISTORY:
!  (1 ) Added AEROSOL_RURALBOX routine (bmy, 9/28/04)
!  (2 ) Now convert ABSHUM from absolute humidity to relative humidity in
!         AEROSOL_RURALBOX, using the same algorithm as in "gasconc.f".
!         (bmy, 1/27/05)
!  (3 ) Now references "tropopause_mod.f" (bmy, 8/22/05)
!  (4 ) Now add contribution of SOA4 into Hydrophilic OC (dkh, bmy, 5/18/06)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Add support for variable tropopause (bdf, phs, 9/14/06)
!  (7 ) Now set OCF=2.1 in AEROSOL_CONC for consistency w/ carbon_mod.f
!       (tmf, 2/10/09)
!  (8 ) Add WTAREA and WERADIUS for dicarbonyl SOA production.  
!       WTAREA is the same as TAREA, but excludes dry dust, BCPO and OCPO; 
!       use same units as TAREA.
!       WERADIUS is same as ERADIUS, but excludes dry dust, BCPO and OCPO;
!       use same units as ERADIUS. (tmf, 3/2/09)
!  (9 ) Add SOAG and SOAM species. (tmf, ccc, 3/2/09)
!  (10) Modify AOD output to wavelength specified in jv_spec_aod.dat 
!       (clh, 05/07/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  05 Mar 2013 - R. Yantosca - Now make INIT_AEROSOL a public routine
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !========================================================================
      ! MODULE VARIABLES:
      !
      ! BCPI        : Hydrophilic black carbon aerosol   [kg/m3]
      ! BCPO        : Hydrophobic black carbon aerosol   [kg/m3]
      ! OCPI        : Hydrophilic organic carbon aerosol [kg/m3]
      ! OCPO        : Hydrophilic organic carbon aerosol [kg/m3]
      ! SALA        : Accumulation mode seasalt aerosol  [kg/m3] 
      ! SALC        : Coarse mode seasalt aerosol        [kg/m3]
      ! SO4_NH4_NIT : Lumped SO4-NH4-NIT aerosol         [kg/m3]
      ! SOILDUST    : Mineral dust aerosol from soils    [kg/m3]
      !========================================================================
      REAL*8, ALLOCATABLE :: BCPI(:,:,:)
      REAL*8, ALLOCATABLE :: BCPO(:,:,:)
      REAL*8, ALLOCATABLE :: OCPI(:,:,:)
      REAL*8, ALLOCATABLE :: OCPO(:,:,:)
      REAL*8, ALLOCATABLE :: SALA(:,:,:)
      REAL*8, ALLOCATABLE :: SALC(:,:,:)
      REAL*8, ALLOCATABLE :: SO4_NH4_NIT(:,:,:)
      REAL*8, ALLOCATABLE :: SOILDUST(:,:,:,:)

      ! Mass of hydrophobic aerosol from Mian Chin
      REAL*8, ALLOCATABLE, SAVE   :: DAERSL(:,:,:,:)

      ! Mass of hydrophilic aerosol from Mian Chin
      REAL*8, ALLOCATABLE, SAVE   :: WAERSL(:,:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerosol_ruralbox
!
! ! DESCRIPTION: Subroutine AEROSOL\_RURALBOX computes quantities that are
!   needed by RDAER. This mimics the call to RURALBOX, which is only done for
!   fullchem runs. (bmy, 9/28/04, 9/14/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AEROSOL_RURALBOX( N_TROP,    am_I_Root, 
     &                             Input_Opt, State_Met )
!
! !USES:
!
      USE COMODE_MOD,         ONLY : ABSHUM, AIRDENS, IXSAVE   
      USE COMODE_MOD,         ONLY : IYSAVE, IZSAVE,  JLOP       
      USE DAO_MOD,            ONLY : MAKE_AVGW
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP

      USE CMN_SIZE_MOD             ! Size parameters
      USE COMODE_LOOP_MOD          ! AD, AVG, WTAIR, etc.
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root ! Is this the root CPU?
!
! !INPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input Options object
      TYPE(MetState), INTENT(INOUT) :: State_Met ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: N_TROP    ! # tropospheric boxes
!
! !REVISION HISTORY: 
!  (1 ) Now convert ABSHUM from absolute humidity to relative humidity in
!        AEROSOL_RURALBOX, using the same algorithm as in "gasconc.f".
!        (bmy, 1/27/05)
!  (2 ) Now references ITS_IN_THE_TROP from "tropopause_mod.f" to diagnose
!        boxes w/in the troposphere. (bmy, 8/22/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Modified for variable tropopause (phs, bdf, 9/14/06)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  08 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER, SAVE        :: N_TROP_BOXES
      INTEGER              :: I, J, L, JLOOP
      REAL*8               :: CONSEXP, TK, VPRESH2O
      REAL*8,  EXTERNAL    :: BOXVL
      
      !=================================================================
      ! AEROSOL_RURALBOX begins here!
      !=================================================================

      ! Initialize
      NLONG  = IIPAR
      NLAT   = JJPAR
      NVERT  = IVERT
      NPVERT = NVERT

      ! Create AVGW field -- mixing ratio of water [v/v]
      CALL MAKE_AVGW( State_Met )

      !=================================================================
      ! Pre-save SMVGEAR loop indices on the first call
      !=================================================================

      ! bdf-phs: must do it everytime with a variable tropopause
      IF ( FIRST .or. Input_Opt%LVARTROP ) THEN

         ! Initialize 1-D index
         JLOOP  = 0

         ! Loop over grid boxes
         DO L = 1, NVERT
         DO J = 1, NLAT
         DO I = 1, NLONG

            ! JLOP is the 1-D grid box loop index
            JLOP(I,J,L) = 0

            !----------------------------------
            ! Boxes w/in ann mean tropopause
            !----------------------------------
            IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) THEN

               ! Increment JLOOP for trop boxes
               JLOOP          = JLOOP + 1

               ! Save JLOOP in SMVGEAR array JLOP
               JLOP(I,J,L)    = JLOOP

               ! These translate JLOOP back to an (I,J,L) triplet
               IXSAVE(JLOOP)  = I
               IYSAVE(JLOOP)  = J
               IZSAVE(JLOOP)  = L                                 

            ENDIF
         ENDDO
         ENDDO
         ENDDO

         ! JLOOP is now the number of boxes w/in GEOS-CHEM's annual mean 
         ! tropopause.  Copy to SAVEd variable N_TROP_BOXES.
         IF ( am_I_Root ) THEN
            WRITE( 6, 100 ) jloop
 100        FORMAT( '  in aerosol ruralbox, val of trop boxes: ', I10)
         ENDIF

         N_TROP_BOXES = JLOOP

         ! Set NTLOOP, NTTLOOP here.  Howeve, we will have to reset these
         ! after the call to READER, since READER redefines these. 
         NTLOOP       = JLOOP
         NTTLOOP      = JLOOP

         ! Reset first-time flag
         FIRST        = .FALSE.
      ENDIF

      ! Copy N_TROP_BOXES to NTROP for passing back to calling program
      N_TROP = N_TROP_BOXES
      
      !=================================================================
      ! Compute AIRDENS and ABSHUM at every timestep
      !
      ! NOTE: In the full-chemistry simulation, SMVGEAR uses the ABSHUM 
      ! array for both absolute humidity [molec H2O/cm3] and relative 
      ! humidity [fraction].  This conversion is done within subroutine 
      ! "gasconc.f", which is called from "chemdr.f".
      !
      ! The computation of aerosol optical depths is done in routine
      ! RDAER of "aerosol_mod.f".  In the full-chemistry simulation,
      ! RDAER is called after "gasconc.f".  At the time when routine
      ! RDAER is called, ABSHUM has already been converted to relative 
      ! humidity.
      !
      ! For the offline aerosol simulation, we must also convert ABSHUM
      ! from absolute humidity to relative humidity using the same
      ! algorithm from "gasconc.f" (see code below).  This will ensure 
      ! that aerosol optical depths in the offline aerosol simulation 
      ! will be  computed in the same way as in the full chemistry 
      ! simulation. (bmy, 1/27/05)
      !=================================================================

      ! Initialize 1-D index
      JLOOP = 0

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP, TK, CONSEXP, VPRESH2O )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, NVERT
      DO J = 1, NLAT
      DO I = 1, NLONG

         ! Get 1-D loop index
         JLOOP = JLOP(I,J,L)

         !----------------------------------
         ! Only process tropospheric boxes
         !----------------------------------
         IF ( JLOOP > 0 ) THEN

            ! Air density in [molec/cm3]
            AIRDENS(JLOOP) = State_Met%AD(I,J,L)    * 1000.d0 / 
     &                       BOXVL(I,J,L,State_Met) * AVG  / WTAIR

            ! ABSHUM = absolute humidity [molec H2O/cm3 air]
            ABSHUM(JLOOP)  = State_Met%AVGW(I,J,L)  * AIRDENS(JLOOP)

            ! Convert ABSHUM to relative humidity [fraction]
            ! using the same algorithm as in "gasconc.f"
            TK             = State_Met%T(I,J,L)

            CONSEXP        = 17.2693882d0      * 
     &                       ( TK - 273.16d0 ) / ( TK - 35.86d0 )
            VPRESH2O       = CONSVAP * EXP( CONSEXP ) / TK
            ABSHUM(JLOOP) = ABSHUM(JLOOP) / VPRESH2O

         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE AEROSOL_RURALBOX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerosol_conc
!
! !DESCRIPTION: Subroutine AEROSOL\_CONC computes aerosol concentrations in
!  kg/m3 from the tracer mass in kg in the STT array.  These are needed to
!  compute optical properties for photolysis and for the optical depth
!  diagnostics. (bmy, 7/20/04, 2/10/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AEROSOL_CONC( am_I_Root, Input_Opt, 
     &                         State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTBCPI, IDTBCPO, IDTDST1, IDTDST2
      USE TRACERID_MOD,       ONLY : IDTDST3, IDTDST4, IDTNH4,  IDTNIT  
      USE TRACERID_MOD,       ONLY : IDTOCPO, IDTOCPI, IDTSALA, IDTSALC 
      USE TRACERID_MOD,       ONLY : IDTSO4,   IDTSOAG, IDTSOAM
      ! SOAupdate: update SOA tracers (hotp 7/25/10)
      USE TRACERID_MOD,       ONLY : IDTPOA1,  IDTPOA2
      USE TRACERID_MOD,       ONLY : IDTOPOA1, IDTOPOA2
      ! SOAupdate: new mtp replaces SOA1-4 (hotp 5/20/10
      USE TRACERID_MOD,       ONLY : IDTTSOA1, IDTTSOA2 
      USE TRACERID_MOD,       ONLY : IDTTSOA3, IDTTSOA0
      USE TRACERID_MOD,       ONLY : IDTISOA1, IDTISOA2, IDTISOA3
      ! SOAupdate: replace SOA5 with lumped aromatic (hotp 5/12/10)
      USE TRACERID_MOD,       ONLY : IDTASOAN, IDTASOA1
      USE TRACERID_MOD,       ONLY : IDTASOA2, IDTASOA3
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  This code was originally included in "chemdr.f", but the same computation 
!  also needs to be done for offline aerosol simulations.  Therefore, we have 
!  split this code off into a separate subroutine which can be called by both 
!  fullchem and offline aerosol simulations.
!
! !REVISION HISTORY: 
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now add contribution from SOA4 into Hydrophilic OC (dkh, bmy, 5/18/06)
!  (3 ) Now set OCF=2.1 to be consistent w/ "carbon_mod.f" (tmf, 2/10/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  13 Nov 2012 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  05 Mar 2013 - R. Yantosca - Remove call to INIT_AEROSOL, this is now done
!                              in the initialization stage
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! We carry carbon mass only in OC and here need to multiply by
      ! 1.4 to account for the mass of the other chemical components
      ! (rjp, bmy, 7/15/04)
      !REAL*8,  PARAMETER :: OCF     = 2.1d0
      ! SOAupdate: semivolpoa: reduce OCF to 1.4 for POA (hotp 2/17/09)
      REAL*8,  PARAMETER :: OCF     = 1.4d0

      ! SOAupdate: semivolpoa4opoa: add OPOA mass (hotp 3/17/09)
      REAL*8,  PARAMETER :: OCFOPOA = 2.1d0  ! = 1.4*1.5

      ! For SOAG, assume the total aerosol mass/glyoxal mass = 1.d0 
      ! for now (tmf, 1/7/09)
      REAL*8,  PARAMETER :: OCFG = 1.d0

      ! For SOAM, assume the total aerosol mass/methylglyoxal mass = 1.d0 
      ! for now (tmf, 1/7/09)
      REAL*8,  PARAMETER :: OCFM = 1.d0
!
! !LOCAL VARIABLES:
!
      ! SAVEd variables
      LOGICAL, SAVE      :: FIRST = .TRUE. 

      ! Non-SAVEd variables
      INTEGER            :: I, J, L, N
      LOGICAL            :: LCARB
      LOGICAL            :: LDUST
      LOGICAL            :: LSOA
      LOGICAL            :: LSSALT
      LOGICAL            :: LSULF

      ! Pointers
      REAL*8, POINTER    :: STT(:,:,:,:)
      REAL*8, POINTER    :: AIRVOL(:,:,:)

      !=================================================================
      ! AEROSOL_CONC begins here!
      !=================================================================

      ! Assume success
      RC     = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LCARB  = Input_Opt%LCARB
      LDUST  = Input_Opt%LDUST
      LSOA   = Input_Opt%LSOA
      LSSALT = Input_Opt%LSSALT
      LSULF  = Input_Opt%LSULF

      ! Initialize pointers
      STT    => State_Chm%Tracers
      AIRVOL => State_Met%AIRVOL

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      
         !==============================================================
         ! S U L F A T E   A E R O S O L S
         !
         ! Dump hydrophilic aerosols into one array that will be passed 
         ! to RDAER and then used for heterogeneous chemistry as well 
         ! as photolysis rate calculations interatively. 
         !
         ! For the full-chemistry run, If LSULF=F, then we read these 
         ! aerosol data from Mian's simulation.  If LSULF=T then we use 
         ! the online tracers.
         !
         ! Now assume that all sulfate, ammonium, and nitrate are 
         ! hydrophilic but sooner or later we can pass only hydrophilic 
         ! aerosols from the thermodynamic calculations for this 
         ! purpose.  This dumping should be done before calling INITGAS, 
         ! which converts the unit of STT from kg/box to molec/cm3.
         !
         ! Units of SO4_NH4_NIT are [kg/m3].  (rjp, bmy, 3/23/03)
         !==============================================================
         IF ( LSULF ) THEN

            ! Compute SO4 aerosol concentration [kg/m3]
            SO4_NH4_NIT(I,J,L) = ( STT(I,J,L,IDTSO4)     + 
     &                             STT(I,J,L,IDTNH4)     +
     &                             STT(I,J,L,IDTNIT) )   /
     &                             AIRVOL(I,J,L)

         ENDIF

         !==============================================================
         ! C A R B O N  &  2 n d A R Y   O R G A N I C   A E R O S O L S
         !
         ! Compute hydrophilic and hydrophobic BC and OC in [kg/m3]
         ! Also add online 2ndary organics if necessary
         !==============================================================
         IF ( LCARB ) THEN

            ! Hydrophilic BC [kg/m3]
            BCPI(I,J,L) = STT(I,J,L,IDTBCPI) / AIRVOL(I,J,L)

            ! Hydrophobic BC [kg/m3]
            BCPO(I,J,L) = STT(I,J,L,IDTBCPO) / AIRVOL(I,J,L)

            ! Hydrophobic OC [kg/m3]
            ! SOAupdate: Treat either OCPO (x2.1) or POA (x1.4) 
            ! POA now 1 and 2 (hotp 10/11/09)
            IF ( IDTPOA1 > 0 .and. IDTPOA2 > 0 ) THEN
               OCPO(I,J,L) = STT(I,J,L,IDTPOA1) * OCF / AIRVOL(I,J,L) +
     &                       STT(I,J,L,IDTPOA2) * OCF / AIRVOL(I,J,L)
            ELSE
               OCPO(I,J,L) = STT(I,J,L,IDTOCPO) * OCFOPOA / 
     &                       AIRVOL(I,J,L)
            ENDIF

            IF ( LSOA ) THEN

               ! Hydrophilic primary OC plus SOA [kg/m3].  We need
               ! to multiply by OCF to account for the mass of other 
               ! components which are attached to the OC aerosol.
               ! (rjp, bmy, 7/15/04)
               ! SOAupdate: Update traditional SOA (hotp 7/21/10)
               ! for new mtp and isop (hotp 5/20/10)
               OCPI(I,J,L) = ( STT(I,J,L,IDTTSOA1)
     &                      +  STT(I,J,L,IDTTSOA2)
     &                      +  STT(I,J,L,IDTTSOA3)
     &                      +  STT(I,J,L,IDTTSOA0)
     &                      +  STT(I,J,L,IDTISOA1)
     &                      +  STT(I,J,L,IDTISOA2)
     &                      +  STT(I,J,L,IDTISOA3) )
     &                      /  AIRVOL(I,J,L) 

               IF ( IDTOPOA1 > 0 .and. IDTOPOA2 > 0 ) THEN  ! hotp 7/28/10
                  OCPI(I,J,L) = OCPI(I,J,L) 
     &                        + ( STT(I,J,L,IDTOPOA1) * OCFOPOA
     &                        / AIRVOL(I,J,L) )
     &                        + ( STT(I,J,L,IDTOPOA2) * OCFOPOA
     &                        / AIRVOL(I,J,L) )
               ENDIF

               IF ( IDTOCPI > 0 ) THEN  ! hotp 7/28/10
                  OCPI(I,J,L) = OCPI(I,J,L)                                  
     &                        + ( STT(I,J,L,IDTOCPI) * OCFOPOA
     &                        / AIRVOL(I,J,L) )
               ENDIF

               ! LUMPAROMIVOC (always present) hotp 7/28/10
               OCPI(I,J,L) = OCPI(I,J,L) +
     &                  ( (  STT(I,J,L,IDTASOAN) 
     &                     + STT(I,J,L,IDTASOA1)
     &                     + STT(I,J,L,IDTASOA2)
     &                     + STT(I,J,L,IDTASOA3) )
     &                     / AIRVOL(I,J,L) )

               ! Check to see if we are simulating SOAG and SOAM (tmf, 1/7/09)
               IF ( IDTSOAG > 0 ) THEN
                 OCPI(I,J,L) = OCPI(I,J,L) + 
     &                         STT(I,J,L,IDTSOAG) * OCFG /
     &                         AIRVOL(I,J,L)
               ENDIF

               IF ( IDTSOAM > 0 ) THEN
                 OCPI(I,J,L) = OCPI(I,J,L) + 
     &                         STT(I,J,L,IDTSOAM) * OCFM /
     &                         AIRVOL(I,J,L)
               ENDIF

            ELSE

               ! Hydrophilic primary and SOA OC [kg/m3].   We need
               ! to multiply by OCF to account for the mass of other 
               ! components which are attached to the OC aerosol.
               ! (rjp, bmy, 7/15/04)
               ! SOAupdate: use 2.1 (OCFOPOA) (hotp 7/21/10)
               OCPI(I,J,L) = STT(I,J,L,IDTOCPI) * OCFOPOA /
     &                       AIRVOL(I,J,L)
                  
            ENDIF ! LSOA

            ! Now avoid division by zero (bmy, 4/20/04)
            BCPI(I,J,L) = MAX( BCPI(I,J,L), 1d-35 )
            OCPI(I,J,L) = MAX( OCPI(I,J,L), 1d-35 )
            BCPO(I,J,L) = MAX( BCPO(I,J,L), 1d-35 )
            OCPO(I,J,L) = MAX( OCPO(I,J,L), 1d-35 )
         
         ENDIF ! LCARB
            
         !===========================================================
         ! M I N E R A L   D U S T   A E R O S O L S
         !
         ! NOTE: We can do better than this! Currently we carry 4 
         ! dust tracers...but het. chem and fast-j use 7 dust bins 
         ! hardwired from Ginoux.
         !
         ! Now, I apportion the first dust tracer into four smallest 
         ! dust bins equally in mass for het. chem and fast-j. 
         ! 
         ! Maybe we need to think about chaning our fast-j and het. 
         ! chem to use just four dust bins or more flexible 
         ! calculations depending on the number of dust bins. 
         ! (rjp, 03/27/04)
         !
         ! Now splitting mass into bins in fractions derived from
         ! Highwood et al. (2003).  Data is from log-normal fit of
         ! PCASP measurements of Saharan dust (Solid line in Fig.4b)
         ! (dar, 04/25/10)
         !===========================================================
         IF ( LDUST ) THEN

            ! Lump 1st dust tracer for het chem
            SOILDUST(I,J,L,1) = 
     &            0.06d0 * STT(I,J,L,IDTDST1) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,2) =
     &            0.12d0 * STT(I,J,L,IDTDST1) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,3) =
     &            0.24d0 * STT(I,J,L,IDTDST1) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,4) =
     &            0.58d0 * STT(I,J,L,IDTDST1) / AIRVOL(I,J,L)

            ! Other hetchem bins
            SOILDUST(I,J,L,5) = STT(I,J,L,IDTDST2) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,6) = STT(I,J,L,IDTDST3) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,7) = STT(I,J,L,IDTDST4) / AIRVOL(I,J,L)
            
         ENDIF

         !===========================================================
         ! S E A S A L T   A E R O S O L S
         !
         ! Compute accumulation & coarse mode concentration [kg/m3]
         !===========================================================
         IF ( LSSALT ) THEN

            ! Accumulation mode seasalt aerosol [kg/m3]
            SALA(I,J,L) = STT(I,J,L,IDTSALA) / AIRVOL(I,J,L)
            
            ! Coarse mode seasalt aerosol [kg/m3]
            SALC(I,J,L) = STT(I,J,L,IDTSALC) / AIRVOL(I,J,L)

            ! Avoid division by zero
            SALA(I,J,L) = MAX( SALA(I,J,L), 1d-35 )
            SALC(I,J,L) = MAX( SALC(I,J,L), 1d-35 )

         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      NULLIFY( STT    )
      NULLIFY( AIRVOL )

      END SUBROUTINE AEROSOL_CONC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdaer
!
! !DESCRIPTION: Subroutine RDAER reads global aerosol concentrations as
!  determined by Mian Chin.  Calculates optical depth at each level for
!  "set\_prof.f". Also calculates surface area for heterogeneous chemistry. It
!  uses aerosol parameters in FAST-J input file "jv\_spec.dat" for these
!  calculations. (rvm, rjp, tdf, bmy, 11/04/01, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RDAER( am_I_Root, Input_Opt, State_Met,  RC,
     &                  MONTH,     YEAR,      WAVELENGTH      )
!
! !USES:
!
      USE CMN_FJ_MOD,         ONLY : JPMAX
      USE CMN_FJ_MOD,         ONLY : JPPJ   
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE COMODE_LOOP_MOD
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT
      USE BPCH2_MOD,          ONLY : GET_RES_EXT
      USE BPCH2_MOD,          ONLY : GET_TAU0
      USE BPCH2_MOD,          ONLY : READ_BPCH2
      USE COMODE_MOD,         ONLY : ABSHUM
      USE COMODE_MOD,         ONLY : ERADIUS
      USE COMODE_MOD,         ONLY : IXSAVE
      USE COMODE_MOD,         ONLY : IYSAVE
      USE COMODE_MOD,         ONLY : IZSAVE
      USE COMODE_MOD,         ONLY : TAREA 
      USE COMODE_MOD,         ONLY : WTAREA
      USE COMODE_MOD,         ONLY : WERADIUS
      USE DIAG_MOD,           ONLY : AD21
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE JV_CMN_MOD
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root     ! Are we on root CPU?
      TYPE(MetState), INTENT(IN)  :: State_Met     ! Meteorology State
      INTEGER,        OPTIONAL    :: MONTH         ! # of current month
      INTEGER,        OPTIONAL    :: YEAR          ! 4-digit year
      INTEGER,        OPTIONAL    :: WAVELENGTH    ! Logical indicator
                                                   !  = 0: AOD computed 
                                                   !       at 999 nm
                                                   !  = 1: AOD computed 
                                                   !       at wavelength 
                                                   !       in jv_spec_aod.dat
      TYPE(OptInput), INTENT(IN)  :: Input_Opt     ! Input Options object
!
! !OUTPUT PARAMETERS:                             
!
      INTEGER,        INTENT(OUT) :: RC            ! Success or failure?
!
! !REVISION HISTORY: 
!  (1 ) At the point in which "rdaer.f" is called, ABSHUM is actually
!        absolute humidity and not relative humidity (rvm, bmy, 2/28/02)
!  (2 ) Now force double-precision arithmetic by using the "D" exponent.
!        (bmy, 2/28/02)
!  (3 ) At present aerosol growth is capped at 90% RH.  The data
!        in jv_spec.dat could be used to allow a particle to grow to
!        99% RH if desired. (rvm, 3/15/02)
!  (4 ) Bug fix: TEMP2 needs to be sized (IIPAR,JJPAR,LLPAR) (bmy, 5/30/02)
!  (5 ) Now reference BXHEIGHT from "dao_mod.f".  Also references ERROR_STOP
!        from "error_mod.f".  Delete local declaration of TIME, since that
!        is also declared w/in comode.h -- this causes compile-time errors
!        on the ALPHA platform. (gcc, bmy, 11/6/02)
!  (6 ) Now use the online SO4, NH4, NIT aerosol, taken from the STT array, 
!        and passed via SO4_NH4_NIT argument if sulfate chemistry is turned on.
!        Otherwise, read monthly mean sulfate from disk.  (rjp, bmy, 3/23/03)
!  (7 ) Now call READ_BPCH2 with QUIET=.TRUE., which prevents info from being
!        printed to stdout.  Also made cosmetic changes. (bmy, 3/27/03)
!  (8 ) Add BCPI, BCPO, OCPI, OCPO to the arg list.  Bug fix: for online
!        sulfate & carbon aerosol tracers, now make sure these get updated
!        every timestep.  Now references "time_mod.f".  Now echo info about
!        which online/offline aerosols we are using.  Updated comments.
!        (bmy, 4/9/04)
!  (9 ) Add SALA, SALC to the arg list (rjp, bec, bmy, 4/20/04)
!  (10) Now references DATA_DIR from "directory_mod.f".  Now references LSULF,
!        LCARB, LSSALT from "logical_mod.f".  Added minor bug fix for 
!        conducting the appropriate scaling for optical depth for ND21
!        diagnostic.  Now make MONTH and YEAR optional arguments.  Now bundled
!        into "aerosol_mod.f".  (rvm, aad, clh, bmy, 7/20/04)
!  (11) Now remove FWET from extinction efficiency computation (avd, 8/3/10)
!  (12) Include third input argument to determine the wavelength at which
!        the AOD should be computed. This will set the optical properties
!        that are used for the calculation of the AOD. The ND21 diagnostic
!        should only be updated when WAVELENGTH = 1. (skim, 02/03/11)
!  09 Mar 2011 - R. Yantosca - Set MSDENS(2) = 1800 for APM (G. Luo)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  13 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC arguments for GIGC
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL             :: FIRST = .TRUE.
      LOGICAL             :: IT_IS_A_FULLCHEM_SIM
      LOGICAL             :: DO_READ_DATA
      CHARACTER(LEN=255)  :: FILENAME
      INTEGER             :: THISMONTH, THISYEAR
      INTEGER             :: I, J, L, N, R, JLOOP, IRH, IRHN
      INTEGER, SAVE       :: MONTH_LAST = -999
      REAL*4              :: TEMP(IIPAR,JJPAR,LGLOB)
      REAL*8              :: TEMP2(IIPAR,JJPAR,LLPAR)
      REAL*8              :: MSDENS(NAER), XTAU, DRYAREA

      ! Effective radius at RH bins read in from "jv_spec.dat"
      REAL*8		  :: RW(NRH)	

      ! Effective radius at RH after interpolation
      REAL*8		  :: REFF       

      ! Q at different RH bins read in from "jv_spec.dat"
      REAL*8		  :: QW(NRH)	
 
      ! Used to interpolate between sizes
      REAL*8		  :: FRAC       
 
      ! Change in Q (extinction efficiency)
      REAL*8		  :: SCALEQ     

      ! Change in Radius with RH
      REAL*8		  :: SCALER     

      ! Chnge in Optical Depth vs RH
      REAL*8		  :: SCALEOD(IIPAR,JJPAR,LLPAR,NRH) 

      ! Change in Vol vs RH 
      REAL*8		  :: SCALEVOL(IIPAR,JJPAR,LLPAR)  

      ! Relative Humidities
      REAL*8,  SAVE       :: RH(NRH)   = (/0d0,0.5d0,0.7d0,0.8d0,0.9d0/)

      ! Index to aerosol types in jv_spec.dat
      ! The following are ordered according to the mass densities below
      INTEGER, SAVE	  :: IND(NAER) = (/22, 29, 36, 43, 50/)

      ! Local variables for quantities from Input_Opt
      CHARACTER(LEN=255)  :: DATA_DIR
      LOGICAL             :: LCARB
      LOGICAL             :: LSSALT
      LOGICAL             :: LSULF

      ! Pointers
      REAL*8, POINTER     :: BXHEIGHT(:,:,:)

      !=================================================================
      ! RDAER begins here!
      !=================================================================

      ! Assume success
      RC                   = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      DATA_DIR             = Input_Opt%DATA_DIR
      IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
      LCARB                = Input_Opt%LCARB
      LSSALT               = Input_Opt%LSSALT
      LSULF                = Input_Opt%LSULF

      ! Initialize pointers
      BXHEIGHT            => State_Met%BXHEIGHT

      ! Copy MONTH argument to local variable THISMONTH
      IF ( PRESENT( MONTH ) ) THEN
         THISMONTH = MONTH
      ELSE
         THISMONTH = 0
      ENDIF

      ! Copy YEAR argument to local variable THISYEAR
      IF ( PRESENT( YEAR ) ) THEN
         THISYEAR = YEAR
      ELSE
         THISYEAR = 0
      ENDIF

      ! Set a logical flag if we have to read data from disk
      ! (once per month, for full-chemistry simulations)
      DO_READ_DATA = ( IT_IS_A_FULLCHEM_SIM .and. ITS_A_NEW_MONTH() )

      !=================================================================
      ! For full-chemistry runs w/ offline fields, define filename
      !=================================================================
      IF ( DO_READ_DATA ) THEN
   
         ! Filename
         FILENAME = TRIM( DATA_DIR ) // 'aerosol_200106/aerosol.' // 
     &              GET_NAME_EXT()   // '.' // GET_RES_EXT()

         ! Use the "generic" year 1996
         XTAU = GET_TAU0( THISMONTH, 1, 1996 )

      ENDIF

      !=================================================================
      ! S U L F A T E   A E R O S O L S
      !
      ! If LSULF = TRUE, then take the lumped SO4, NH4, NIT 
      ! concentrations [kg/m3] computed by AEROSOL_CONC, and save 
      ! into WAERSL(:,:,:,1) for use w/ FAST-J and hetchem.  This is 
      ! updated every timestep.  (For fullchem and offline runs)
      !
      ! If LSULF = FALSE, then read monthly mean offline sulfate aerosol   
      ! concentrations [kg/m3] from disk at the start of each month.
      ! (For fullchem simulations only)
      !=================================================================
      IF ( LSULF ) THEN 

         !-----------------------------------
         ! Use online aerosol concentrations
         !-----------------------------------
         IF ( FIRST ) THEN
            IF ( am_I_Root ) WRITE( 6, 100 ) 
 100        FORMAT( '     - RDAER: Using online SO4 NH4 NIT!' )
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            WAERSL(I,J,L,1) = SO4_NH4_NIT(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE 

         !-----------------------------------
         ! Read from disk -- fullchem only
         !-----------------------------------
         IF ( DO_READ_DATA ) THEN
            
            ! Print filename
            IF ( am_I_Root ) WRITE( 6, 105 ) TRIM( FILENAME )
 105        FORMAT( '     - RDAER: Reading SULFATE from ', a )

            ! Read data
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 1,     
     &                       XTAU,      IIPAR,     JJPAR,     
     &                       LGLOB,     TEMP,      QUIET=.TRUE. )

            ! Cast to REAL*8 and resize
            CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,1) )
         ENDIF
      ENDIF

      !=================================================================
      ! C A R B O N  &  2 n d A R Y   O R G A N I C   A E R O S O L S
      !
      ! If LCARB = TRUE, then take Hydrophilic OC, Hydrophobic OC,
      ! Hydropilic BC, and Hydrophobic BC, and 2ndary organic aerosol
      ! concentrations [kg/m3] that have been computed by AEROSOL_CONC.   
      ! Save these into DAERSL and WAERSL for use w/ FAST-J and hetchem.  
      ! These fields are updated every chemistry timestep.
      ! (For both fullchem and offline simulations)
      !
      ! If LCARB = FALSE, then read monthly mean carbon aerosol
      ! concentrations [kg/m3] from disk at the start of each month.
      ! (For full chemistry simulations only)
      !=================================================================
      IF ( LCARB ) THEN

         !-----------------------------------
         ! Use online aerosol concentrations
         !-----------------------------------
         IF ( FIRST ) THEN
            IF ( am_I_Root ) WRITE( 6, 110 ) 
 110        FORMAT( '     - RDAER: Using online BCPI OCPI BCPO OCPO!' )
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Hydrophilic BC (a.k.a EC) [kg/m3]
            WAERSL(I,J,L,2) = BCPI(I,J,L)

            ! Hydrophilic OC [kg/m3]
            WAERSL(I,J,L,3) = OCPI(I,J,L)

            ! Hydrophobic BC (a.k.a EC) [kg/m3]
            DAERSL(I,J,L,1) = BCPO(I,J,L)

            ! Hydrophobic OC [kg/m3]
            DAERSL(I,J,L,2) = OCPO(I,J,L)

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE 

         !-----------------------------------
         ! Read from disk -- fullchem only
         !-----------------------------------
         IF ( DO_READ_DATA ) THEN

            ! Print filename
            IF ( am_I_Root ) WRITE( 6, 115 ) TRIM( FILENAME )
 115        FORMAT( '     - RDAER: Reading BC and OC from ', a )

            !--------------------------------
            ! Read Hydrophobic BC
            !--------------------------------
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 2,     
     &                       XTAU,      IIPAR,     JJPAR,     
     &                       LGLOB,     TEMP,      QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, DAERSL(:,:,:,1) )

            !---------------------------------
            ! Read Hydrophilic BC
            !---------------------------------
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 3,     
     &                       XTAU,      IIPAR,     JJPAR,     
     &                       LGLOB,     TEMP,      QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,2) )

            !---------------------------------
            ! Read Hydrophobic OC
            !---------------------------------
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 4,     
     &                       XTAU,      IIPAR,     JJPAR,     
     &                       LGLOB,     TEMP,      QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, DAERSL(:,:,:,2) )

            !---------------------------------
            ! Read Hydrophilic OC
            !---------------------------------
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 5,     
     &                       XTAU,      IIPAR,     JJPAR,     
     &                       LGLOB,     TEMP,      QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,3) )
         ENDIF
      ENDIF

      !=================================================================
      ! S E A S A L T   A E R O S O L S
      !
      ! If LSSALT = TRUE, then take accumulation and coarse mode
      ! seasalt aerosol concentrations [kg/m3] that are passed from
      ! "chemdr.f".  Save these into WAERSL for use w/ FAST-J and
      ! hetchem.  These fields are updated every chemistry timestep.
      ! (For both fullchem and offline simulations)
      !
      ! If LSSALT = FALSE, then read monthly-mean coarse sea-salt 
      ! aerosol concentrations [kg/m3] from the binary punch file.  
      ! Also merge the coarse sea salt aerosols into a combined bin 
      ! rather than carrying them separately.
      ! (For fullchem simulations only)
      !=================================================================
      IF ( LSSALT ) THEN

         !-----------------------------------
         ! Use online aerosol concentrations
         !-----------------------------------
         IF ( FIRST ) THEN
            IF ( am_I_Root ) WRITE( 6, 120 ) 
 120        FORMAT( '     - RDAER: Using online SALA SALC' )
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Accumulation mode seasalt aerosol [kg/m3]
            WAERSL(I,J,L,4) = SALA(I,J,L)

            ! Coarse mode seasalt aerosol [kg/m3]
            WAERSL(I,J,L,5) = SALC(I,J,L)

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE
            
         !-----------------------------------
         ! Read from disk -- fullchem only
         !-----------------------------------
         IF ( DO_READ_DATA ) THEN

            ! Print filename
            IF ( am_I_Root ) WRITE( 6, 125 ) TRIM( FILENAME )
 125        FORMAT( '     - RDAER: Reading SEASALT from ', a )

            !----------------------------------
            ! Offline -- read Sea Salt (accum)
            !----------------------------------
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 6,     
     &                       XTAU,      IIPAR,     JJPAR,     
     &                       LGLOB,     TEMP,      QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,4) )

            !----------------------------------
            ! Offline -- read Sea Salt (coarse)
            !----------------------------------
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 7,     
     &                       XTAU,      IIPAR,     JJPAR,     
     &                       LGLOB,     TEMP,      QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,5) )

            !----------------------------------
            ! Offline -- read Sea Salt (coarse)
            !----------------------------------
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 8,     
     &                       XTAU,      IIPAR,     JJPAR,     
     &                       LGLOB,     TEMP,      QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, TEMP2 )

            ! Accumulate into one size bin
            WAERSL(:,:,:,5) = WAERSL(:,:,:,5) + TEMP2 

            !----------------------------------
            ! Offline -- read Sea Salt (coarse)
            !----------------------------------
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 9,     
     &                       XTAU,      IIPAR,     JJPAR,     
     &                       LGLOB,     TEMP,      QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, TEMP2 )

            ! Accumulate into one size bin
            WAERSL(:,:,:,5) = WAERSL(:,:,:,5) + TEMP2 

         ENDIF   
      ENDIF 

      !=================================================================
      ! Calculate optical depth and surface area at each timestep
      ! to account for the change in relative humidity
      !
      ! For the optical depth calculation, this involves carrying the 
      ! optical depth at each RH as separate aerosols since OPMIE.f 
      ! treats the phase functions and single scattering albedos 
      ! separately. (An alternative would be to rewrite OPMIE.f)
      !
      ! Scaling is sufficient for the surface area calculation
      !=================================================================
      MSDENS(1) = 1700.0d0    !SO4
#if   defined( APM )
      MSDENS(2) = 1800.0d0    !BC 
#else
      MSDENS(2) = 1000.0d0    !BC 
#endif
      MSDENS(3) = 1800.0d0    !OC 
      MSDENS(4) = 2200.0d0    !SS (accum)
      MSDENS(5) = 2200.0d0    !SS (coarse)

      ! Loop over types of aerosol
      DO N = 1, NAER

         ! Zero array
         SCALEOD(:,:,:,:) = 0d0
         
         !==============================================================
         ! Determine aerosol growth rates from the relative 
         ! humidity in each box
         !
         ! The optical depth scales with the radius and Q alone
         ! since SCALEDENS cancels as follows
         ! 
         !    SCALER 	= RW / RDRY
         !    SCALEDENS = DENSWET / DENSDRY
         !    SCALEM 	= SCALEDENS * SCALER**3
         !    SCALEOD 	= (SCALEQ * SCALEM) / (SCALEDENS * SCALER)
         !          	= SCALEQ * SCALER**2
         !
         ! Cap aerosol values at 90% relative humidity since
         ! aerosol growth at that point becomes highly nonlinear and 
         ! relative humidities above this value essentially mean
         ! there is a cloud in that grid box
         !
         ! Q is the extinction efficiency
         !
         ! Each grid box (I,J,L) will fall into one of the RH bins, 
         ! since each grid box will have a different RH value.  So,
         ! for SCALEOD(I,J,L,:), only one of the IRH bins will contain
         ! nonzero data, while the other IRH bins will all be zero.
         !==============================================================
 
         ! Loop over relative humidity bins
         DO R = 1, NRH

            IF ( WAVELENGTH > 0 ) THEN

               ! Now use the optical properties at the wavelength specified
               !  in jv_spec_aod.dat (skim, 02/03/11)

               ! Wet radius in "jv_spec_aod.dat"
               RW(R) = RAA_AOD(IND(N)+R-1)
               
               ! Extinction efficiency for Q for each RH bin
               QW(R) = QAA_AOD(IND(N)+R-1)

            ELSE 

               ! Wet radius in "jv_spec.dat"
               RW(R) = RAA(4,IND(N)+R-1)	

               ! Extinction efficiency Q for each RH bin
               QW(R) = QAA(4,IND(N)+R-1)

            ENDIF

         ENDDO

         ! Loop over SMVGEAR grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, IRH, JLOOP, SCALEQ, SCALER, REFF, FRAC )
!$OMP+SCHEDULE( DYNAMIC )
         DO JLOOP = 1, NTLOOP

            ! Get 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)

            ! Sort into relative humidity bins
            IF (      ABSHUM(JLOOP) <= RH(2) ) THEN   
               IRH = 1
            ELSE IF ( ABSHUM(JLOOP) <= RH(3) ) THEN
               IRH = 2
            ELSE IF ( ABSHUM(JLOOP) <= RH(4) ) THEN
               IRH = 3
            ELSE IF ( ABSHUM(JLOOP) <= RH(5) ) THEN
               IRH = 4
            ELSE 
               IRH = 5
            ENDIF

            ! For the NRHth bin, we don't have to interpolate
            ! For the other bins, we have to interpolate 
            IF ( IRH == NRH ) THEN
               SCALEQ = QW(NRH) / QW(1)  !QW(1) is dry extinction eff.
               REFF   = RW(NRH) 

            ELSE                

               ! Interpolate between different RH
               FRAC = (ABSHUM(JLOOP)-RH(IRH)) / (RH(IRH+1)-RH(IRH))
               IF ( FRAC > 1.0d0 ) FRAC = 1.0d0
               
               SCALEQ = (FRAC*QW(IRH+1) + (1.d0-FRAC)*QW(IRH)) / QW(1)
               REFF   = FRAC*RW(IRH+1)  + (1.d0-FRAC)*RW(IRH)

            ENDIF

            SCALER                 = REFF / RW(1)
            SCALEOD(I,J,L,IRH)     = SCALEQ * SCALER * SCALER
            SCALEVOL(I,J,L)        = SCALER**3
            ERADIUS(JLOOP,NDUST+N) = 1.0D-4 * REFF

            !==============================================================
            ! ND21 Diagnostic: 
            !
            ! Computed here:
            ! --------------
            ! #7  Hygroscopic growth of SO4                [unitless]
            ! #10 Hygroscopic growth of Black Carbon       [unitless]
            ! #13 Hygroscopic growth of Organic Carbon     [unitless]
            ! #16 Hygroscopic growth of Sea Salt (accum)   [unitless]
            ! #19 Hygroscopic growth of Sea Salt (coarse)  [unitless]
            !==============================================================
            IF ( ND21 > 0 .and. L <= LD21 .and. WAVELENGTH > 0 ) THEN
               AD21(I,J,L,4+3*N) = AD21(I,J,L,4+3*N) +SCALEOD(I,J,L,IRH)
            ENDIF

         ENDDO
!$OMP END PARALLEL DO
      
         !==============================================================
         ! Convert concentration [kg/m3] to optical depth [unitless].
         !
         ! ODAER = ( 0.75 * BXHEIGHT * AERSL * QAA ) / 
         !         ( MSDENS * RAA * 1e-6 )
         ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
         !
         ! Units ==> AERSL    [ kg/m3    ]
         !           MSDENS   [ kg/m3    ]
         !           RAA      [ um       ]  
         !           BXHEIGHT [ m        ]
         !           QAA      [ unitless ]
         !           ODAER    [ unitless ]
         !
         ! NOTES: 
         ! (1 ) Do the calculation at QAA(4,:) (i.e. 999 nm).          
         ! (2 ) RAA is the 'effective radius', Hansen and Travis, 1974
         ! (3 ) Report at the more relevant QAA(2,:) (i.e. 400 nm)   
         !       Although SCALEOD would be slightly different at 400nm 
         !       than at 1000nm as done here, FAST-J does currently 
         !       allow one to provide different input optical depths at 
         !       different wavelengths.  Therefore the reported value at
         !       determined with QAA(2,:) is as used in FAST-J. 
         ! (4 ) Now use explicit indices in parallel DO-loops, since
         !       some compilers may not like array masks in parallel
         !       regions (bmy, 2/28/02)
         !==============================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, R, IRHN ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO R = 1, NRH
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Bin for aerosol type and relative humidity
            IRHN = ( (N-1) * NRH ) + R

            IF ( WAVELENGTH > 0 ) THEN

            ! Save aerosol optical depth for each combination 
            ! of aerosol type and relative humidity into ODAER, 
            ! which will get passed to the FAST-J routines
            ODAER(I,J,L,IRHN) = SCALEOD(I,J,L,R) * 0.75d0
     &                        * BXHEIGHT(I,J,L) 
     &                        * WAERSL(I,J,L,N)  * QAA_AOD(IND(N)) / 
     &                        ( MSDENS(N) * RAA_AOD(IND(N)) * 1.0D-6 )

            ELSE

            ! Save aerosol optical depth for each combination 
            ! of aerosol type and relative humidity into ODAER, 
            ! which will get passed to the FAST-J routines
            ODAER(I,J,L,IRHN) = SCALEOD(I,J,L,R) * 0.75d0 
     &                        * BXHEIGHT(I,J,L) 
     &                        * WAERSL(I,J,L,N)  * QAA(4,IND(N)) / 
     &                        ( MSDENS(N) * RAA(4,IND(N)) * 1.0D-6 )

            ENDIF

         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         !  Calculate Aerosol Surface Area
         !
         !  Units ==> AERSL    [ kg aerosol m^-3 air ]
         !            MSDENS   [ kg aerosol m^-3 aerosol ]
         !            ERADIUS  [ cm      ]
         !            TAREA    [ cm^2 dry aerosol/cm^3 air ]
         !
         !  Note: first find volume of aerosol (cm^3 arsl/cm^3 air), then
         !        multiply by 3/radius to convert to surface area in cm^2
         !
         !  Wet Volume = AERSL * SCALER**3 / MSDENS
         !  Wet Surface Area = 3 * (Wet Volume) / ERADIUS 
         !
         !  Effective radius for surface area and optical depths 
         !  are identical.
         !==============================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, JLOOP ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO JLOOP = 1, NTLOOP

            ! Get 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)

            !========================================================
            ! NOTES:          
            !    WAERSL   [ kg dry mass of wet aerosol m^-3 air ]
            !    ERADIUS  [ cm wet aerosol radius ]
            !    MSDENS   [ kg dry mass of aerosol m^-3 dry volume of aerosol ]
            !    TAREA    [ cm^2 wet sfc area of aerosol cm^-3 air ]   
            !    WTAREA   : same as TAREA, but excludes dry dust, BCPO and OCPO 
            !               use same units as TAREA    (tmf, 4/18/07) 
            !    WERADIUS : same as ERADIUS, but excludes dry dust, BCPO and OCPO
            !               use same units as ERADIUS  (tmf, 4/18/07)
            ! Wet dust WTAREA and WERADIUS are archived in dust_mod.f.
            !========================================================
          
            ! Store aerosol surface areas in TAREA, and be sure
            ! to list them following the dust surface areas
            TAREA(JLOOP,N+NDUST) = 3.D0                     * 
     &                             WAERSL(I,J,L,N)          *  
     &                             SCALEVOL(I,J,L)          / 
     &                             ( ERADIUS(JLOOP,NDUST+N) * 
     &                               MSDENS(N) )  

            WTAREA(JLOOP, N)   = TAREA(JLOOP, N+NDUST)
            WERADIUS(JLOOP, N) = ERADIUS(JLOOP, N+NDUST)


         ENDDO
!$OMP END PARALLEL DO

      ENDDO  !Loop over NAER

      !==============================================================
      ! Account for hydrophobic aerosols (BC and OC), N=2 and N=3
      !==============================================================
      DO N = 2, 3

         ! Index for combination of aerosol type and RH
         IRHN = ( (N-1) * NRH ) + 1

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF ( WAVELENGTH > 0 ) THEN

            ! Aerosol optical depth
            ODAER(I,J,L,IRHN) = ODAER(I,J,L,IRHN) + 0.75d0 *
     &                          BXHEIGHT(I,J,L)  * 
     &                          DAERSL(I,J,L,N-1) * QAA_AOD(IND(N))   / 
     &                          ( MSDENS(N) * RAA_AOD(IND(N)) * 1.0D-6 )

            ELSE

            ! Aerosol optical depth
            ODAER(I,J,L,IRHN) = ODAER(I,J,L,IRHN) + 0.75d0 *
     &                          BXHEIGHT(I,J,L)  * 
     &                          DAERSL(I,J,L,N-1) * QAA(4,IND(N))   / 
     &                          ( MSDENS(N) * RAA(4,IND(N)) * 1.0D-6 )

            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Effective radius
         IF ( WAVELENGTH > 0 ) THEN
            REFF = 1.0D-4 * RAA_AOD(IND(N))
         ELSE
            REFF = 1.0D-4 * RAA(4,IND(N))
         ENDIF

         ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, JLOOP, DRYAREA ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO JLOOP = 1, NTLOOP

            ! Get 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)

            ! Dry surface area
            DRYAREA = 3.D0 * DAERSL(I,J,L,N-1) / ( REFF * MSDENS(N) )  

            ! Add surface area to TAREA array
            TAREA(JLOOP,N+NDUST) = TAREA(JLOOP,N+NDUST) + DRYAREA

            ! Define a new effective radius that accounts 
            ! for the hydrophobic aerosol 
            ERADIUS(JLOOP,NDUST+N) = ( ERADIUS(JLOOP,NDUST+N) * 
     &                                  TAREA(JLOOP,N+NDUST)  +
     &                                  REFF * DRYAREA)       / 
     &                               ( TAREA(JLOOP,N+NDUST) + DRYAREA )

         ENDDO
!$OMP END PARALLEL DO

      ENDDO
       
      !==============================================================
      ! ND21 Diagnostic: Aerosol OD's, Growth Rates, Surface Areas
      !
      ! Computed in other routines:
      ! ---------------------------------
      ! #1: Cloud optical depths (1000 nm) --> from "optdepth_mod.f"       
      ! #2: Max Overlap Cld Frac           --> from "optdepth_mod.f" 
      ! #3: Random Overlap Cld Frac        --> from "optdepth_mod.f" 
      ! #4: Dust optical depths            --> from "rdust.f"
      ! #5: Dust surface areas             --> from "rdust.f"
      !
      ! Computed previously in "rdaer.f":
      ! ---------------------------------
      ! #7  Hygroscopic growth of SO4                [unitless]
      ! #10 Hygroscopic growth of Black Carbon       [unitless]
      ! #13 Hygroscopic growth of Organic Carbon     [unitless]
      ! #16 Hygroscopic growth of Sea Salt (accum)   [unitless]
      ! #19 Hygroscopic growth of Sea Salt (coarse)  [unitless]
      !
      ! Computed here:
      ! ---------------------------------
      ! #6  Sulfate Optical Depth (400 nm)           [unitless]
      ! #8  Sulfate Surface Area                     [cm2/cm3 ]
      ! #9  Black Carbon Optical Depth (400 nm)      [unitless]
      ! #11 Black Carbon Surface Area                [cm2/cm3 ]
      ! #12 Organic Carbon Optical Depth (400 nm)    [unitless]
      ! #14 Organic Carbon Surface Area              [cm2/cm3 ]
      ! #15 Sea Salt (accum) Opt Depth (400 nm)      [unitless]
      ! #17 Sea Salt (accum) Surface Area            [cm2/cm3 ]
      ! #18 Sea Salt (coarse) Opt Depth(400 nm)      [unitless]
      ! #20 Sea Salt (coarse) Surface Area           [cm2/cm3 ]
      !
      ! #21: Dust optical depths (0.15 um)     --> from "rdust.f"
      ! #22: Dust optical depths (0.25 um)     --> from "rdust.f"
      ! #23: Dust optical depths (0.4 um)     --> from "rdust.f"
      ! #24: Dust optical depths (0.8 um)     --> from "rdust.f"
      ! #25: Dust optical depths (1.5 um)     --> from "rdust.f"
      ! #26: Dust optical depths (2.5 um)     --> from "rdust.f"
      ! #27: Dust optical depths (4.0 um)     --> from "rdust.f"

      ! NOTE: The cloud optical depths are actually recorded at
      !       1000 nm, but vary little with wavelength.
      !==============================================================
      IF ( ND21 > 0 .and. WAVELENGTH > 0 ) THEN

         ! Loop over aerosol types (dust handled in dust_mod.f)
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, IRHN, J, JLOOP, L, N, R ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO N = 1, NAER

            !------------------------------------
            ! Aerosol Optical Depths [uhitless]
            ! Scale of optical depths w/ RH 
            !------------------------------------
            DO R = 1, NRH
            DO L = 1, LD21
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Index for type of aerosol and RH value
               IRHN = ( (N-1) * NRH ) + R

               ! Optical Depths (scaled to jv_spec_aod.dat wavelength, clh)
               ! Scaling based on wavelength no long necessary (skim)
               AD21(I,J,L,3+3*N) = AD21(I,J,L,3+3*N) + ODAER(I,J,L,IRHN)

            ENDDO
            ENDDO
            ENDDO
            ENDDO

            !------------------------------------
            ! Aerosol Surface Areas [cm2/cm3]
            !------------------------------------
            DO JLOOP = 1, NTLOOP

               ! Get 3-D grid box indices
               I = IXSAVE(JLOOP)
               J = IYSAVE(JLOOP)
               L = IZSAVE(JLOOP)

               ! Add aerosol surface areas 
               IF ( L <= LD21 ) THEN
                  AD21(I,J,L,5+3*N) = AD21(I,J,L,5+3*N) + 
     &                                TAREA(JLOOP,N+NDUST)
               ENDIF 
            ENDDO

         ENDDO
!$OMP END PARALLEL DO

      ENDIF 

      !=================================================================
      ! To turn off the radiative effects of different aerososl
      ! uncomment the following lines
      !=================================================================
      !DO R = 1,NRH
      !  ODAER(:,:,:,R)       = 0.d0  !sulfate
      !  ODAER(:,:,:,R+NRH)   = 0.d0  !BC
      !  ODAER(:,:,:,R+2*NRH) = 0.d0  !OC
      !  ODAER(:,:,:,R+3*NRH) = 0.d0  !SS(accum)
      !  ODAER(:,:,:,R+4*NRH) = 0.d0  !SS(coarse)
      !ENDDO

      !=================================================================
      ! To turn off heterogeneous chemistry on different aerosols
      ! uncomment the following lines
      !=================================================================
      !TAREA(:,NDUST+1) = 0.d0	!Sulfate
      !TAREA(:,NDUST+2) = 0.d0	!BC 
      !TAREA(:,NDUST+3) = 0.d0	!OC 
      !TAREA(:,NDUST+4) = 0.d0	!SS (accum)
      !TAREA(:,NDUST+5) = 0.d0	!SS (coarse)

      ! Free pointers
      NULLIFY( BXHEIGHT )

      ! Reset first-time flag
      FIRST = .FALSE.

      END SUBROUTINE RDAER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_aerosol
!
! !DESCRIPTION: Subroutine INIT\_AEROSOL allocates and zeroes module arrays
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_AEROSOL( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  05 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_AEROSOL begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate arrays
      ALLOCATE( BCPI( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPI' )
      BCPI = 0d0

      ALLOCATE( BCPO( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPO' )
      BCPO = 0d0

      ALLOCATE( OCPI( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPI' )
      OCPI = 0d0

      ALLOCATE( OCPO( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPO' )
      OCPO = 0d0

      ALLOCATE( SALA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SALA' )
      SALA = 0d0

      ALLOCATE( SALC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SALC' )
      SALC = 0d0

      ALLOCATE( SO4_NH4_NIT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_NH4_NIT' )
      SO4_NH4_NIT = 0d0

      ALLOCATE( SOILDUST( IIPAR, JJPAR, LLPAR, NDUST ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOILDUST' )
      SOILDUST = 0d0

      ! Mass of hydrophobic aerosol from Mian Chin
      ALLOCATE( DAERSL(IIPAR,JJPAR,LLPAR,2), STAT=AS)       
      IF ( AS /= 0) CALL ALLOC_ERR( 'DAERSL' )
      DAERSL = 0d0

      ! Mass of hydrophilic aerosol from Mian Chin
      ALLOCATE( WAERSL(IIPAR,JJPAR,LLPAR,NAER), STAT=AS)       
      IF ( AS /= 0) CALL ALLOC_ERR( 'WAERSL' )
      DAERSL = 0d0

      ! Return to calling program

      END SUBROUTINE INIT_AEROSOL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_aerosol
!
! !DESCRIPTION: Subroutine CLEANUP\_AEROSOL deallocates all module arrays
!  (bmy, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_AEROSOL
! 
! !REVISION HISTORY:
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!

      !=================================================================
      ! CLEANUP_AEROSOL begins here!
      !=================================================================
      IF ( ALLOCATED( BCPI        ) ) DEALLOCATE( BCPI        )
      IF ( ALLOCATED( BCPO        ) ) DEALLOCATE( BCPO        )
      IF ( ALLOCATED( OCPI        ) ) DEALLOCATE( OCPI        )
      IF ( ALLOCATED( OCPO        ) ) DEALLOCATE( OCPO        )
      IF ( ALLOCATED( SALA        ) ) DEALLOCATE( SALA        )
      IF ( ALLOCATED( SALC        ) ) DEALLOCATE( SALC        )
      IF ( ALLOCATED( SO4_NH4_NIT ) ) DEALLOCATE( SO4_NH4_NIT )
      IF ( ALLOCATED( SOILDUST    ) ) DEALLOCATE( SOILDUST    )

      END SUBROUTINE CLEANUP_AEROSOL
!EOC
      END MODULE AEROSOL_MOD
