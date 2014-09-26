!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: dust_mod
!
! !DESCRIPTION: Module DUST\_MOD contains routines for computing dust aerosol 
!  emissions, chemistry, and optical depths.
!\\
!\\
! !INTERFACE: 
!
      MODULE DUST_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CHEMDUST      
#if   defined( TOMAS )
      PUBLIC  :: SETTLEDUST
#endif
      PUBLIC  :: RDUST_ONLINE
      PUBLIC  :: RDUST_OFFLINE
      PUBLIC  :: INIT_DUST
      PUBLIC  :: CLEANUP_DUST
!     
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: DRY_SETTLING
      PRIVATE :: DRY_DEPOSITION
!
!  !REVISION HISTORY:
!  30 Mar 2004 - T. D. Fairlie - Initial version
!  (1 ) Bug fix in SRC_DUST_DEAD (bmy, 4/14/04)
!  (2 ) Now references "logical_mod.f", "directory_mod.f", and "tracer_mod.f"
!        Added comments. (bmy, 7/2/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) Bug fix in snow height computation (bmy, 11/18/05)
!  (6 ) Now only do drydep if LDRYD=T (bmy, 5/23/06)
!  (7 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (8 ) Updated output print statement in SRC_DUST_DEAD (bmy, 1/23/07)
!  (9 ) Modifications for GEOS-5 (bmy, 1/24/07)
!  (10) Modified to archive only hydrophilic aerosol/aqueous dust surface area
!        (excluding BCPO and OCPO) for aqueous chemistry calculations 
!        Dust surfaces are considered aqueous only when RH > 35% (tmf, 3/6/09)
!  (11) Add AOD output for all dust size bins (clh, 5/7/10)
!  (12) Modify AOD output to wavelength specified in jv_spec_aod.dat 
!       (clh, 05/07/10)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  03 Sep 2010 - R. Yantosca - Bug fix in SRC_DUST_DEAD
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x
!  01 Mar 2012 - R. Yantosca - Now reference the new grid_mod.F90
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  14 Nov 2012 - R. Yantosca - Add modifications for GIGC
!  04 Mar 2013 - R. Yantosca - Now call INIT_DUST from the init stage
!                              which facilitates connection to GEOS-5 GCM
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  20 Jun 2014 - R. Yantosca - Remove obsolete emissions code; we now use HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      INTEGER              :: DRYDST1, DRYDST2, DRYDST3, DRYDST4
      REAL*8,  ALLOCATABLE :: FRAC_S(:)
      REAL*8,  ALLOCATABLE :: SRCE_FUNC(:,:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemdust
!
! !DESCRIPTION: Subroutine CHEMDUST is the interface between the GEOS-Chem 
!  main program and the dust chemistry routines that mostly calculates dust
!  dry deposition.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMDUST( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTDST1
      USE TRACERID_MOD,       ONLY : IDTDST2
      USE TRACERID_MOD,       ONLY : IDTDST3
      USE TRACERID_MOD,       ONLY : IDTDST4
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
!  30 Mar 2004 - T. D. Fairlie - Initial version
!  (1 ) Now references STT from "tracer_mod.f" and LDUST from "logical_mod.f"
!        (bmy, 7/20/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now only do dry deposition if LDRYD = T (bmy, 5/23/06)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Now pass State_Met as an argument
!  05 Mar 2013 - R. Yantosca - Add ND70 debug print output
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL, SAVE :: FIRST = .TRUE.

      ! Non-SAVEd scalars
      LOGICAL       :: LDRYD
      LOGICAL       :: LDUST
      LOGICAL       :: prtDebug
      INTEGER       :: N

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEMDUST begins here!
      !=================================================================

      ! Assume success
      RC    = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LDRYD = Input_Opt%LDRYD
      LDUST = Input_Opt%LDUST

      ! Set a flag for debug output
      prtDebug = ( Input_Opt%LPrt .and. am_I_Root )

      ! Execute on first call only
      IF ( FIRST ) THEN
 
         ! Stop w/ error if dust tracer flags are undefined
         IF ( IDTDST1 + IDTDST2 + IDTDST3 + IDTDST4 == 0 ) THEN
            IF ( LDUST ) THEN 
               CALL ERROR_STOP( 
     &              'LDUST=T but dust tracers are undefined!',
     &              'EMISSDUST ("dust_mod.f")' )
            ENDIF
         ENDIF

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Do dust settling & deposition
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      !----------------------
      ! Dust settling
      !----------------------
      CALL DRY_SETTLING(  am_I_Root, Input_Opt,  
     &                    State_Met, STT(:,:,:,IDTDST1:IDTDST4), RC )

      ! Debug print
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDUST: a DRY_SETTLING' )
      ENDIF

      !----------------------
      ! Dust deposition
      !----------------------
      IF ( LDRYD ) THEN
         CALL DRY_DEPOSITION( am_I_Root,                  Input_Opt, 
     &                        STT(:,:,:,IDTDST1:IDTDST4), RC         )

         ! Debug print
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMDUST: a DRY_DEPOSITION' )
         ENDIF
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEMDUST
!EOC


#if   defined( TOMAS )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: settledust 
!
! !DESCRIPTION: Subroutine SETTLEDUST is the interface between the 
!  size-resolved dry deposition subroutine AERO\_DRYDEP and the dust module. 
!  This is to call only gravitational settling and deals with removal of 
!  aerosol number with the dust mass.  (win, 7/17/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETTLEDUST( am_I_Root, Input_Opt, 
     &                       State_Met, State_Chm, RC )
!
! !USES:
! 
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : NUMDEP, NTRAIND
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TOMAS_MOD,          ONLY : IBINS, Xk, SRTDUST
      USE TRACERID_MOD,       ONLY : IDTDUST1, IDTNK1

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND44
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
!  17 Jul 2009 - W. Trivitayanurak - Initial version
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  13 Dec 2012 - M. Payer    - Add am_I_Root, Input_Opt, RC as arguments
!  31 May 2013 - R. Yantosca - Now pass State_Chm via the arg list
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL, SAVE :: FIRST = .TRUE.

      ! Non-SAVEd scalars
      INTEGER       :: N, BIN, I, J, L, NN
      REAL*8        :: DU0(IIPAR,JJPAR,LLPAR,IBINS)
      REAL*8        :: DIF, FLUXN, FLUXD
      REAL*8        :: DT_SETTL, AREA_CM2

      !debug
      integer       :: ii, jj , ix, jx, bb
      data ii,jj, ix, jx, bb /37, 24, 58, 34, 30 /

      ! For values from Input_Opt
      INTEGER         :: N_TRACERS
      REAL*8          :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER   :: STT(:,:,:,:)

      !=================================================================
      ! SETTLEDUST begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy values from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT       => State_Chm%Tracers

      !=================================================================
      ! Do dust settling & deposition
      !=================================================================

      ! Dust settling timestep [s]
      DT_SETTL = GET_TS_CHEM() * 60d0

      ! Save initial dust mass
      DO BIN = 1, IBINS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         DU0(I,J,L,BIN) = STT(I,J,L,IDTDUST1-1+BIN) 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ! Dust settling
      CALL DRY_SETTLING( am_I_Root,
     &                   Input_Opt,  
     &                   State_Met,
     &                   STT(:,:,:,IDTDUST1:IDTDUST1-1+IBINS),
     &                   RC )

      ! Calculate change in number to correspond with dust redistribution
      ! by gravitational settling
      DO BIN = 1, IBINS
      NN = NUMDEP + (SRTDUST-1)*IBINS + BIN

      DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )
         
            FLUXD = 0d0
            FLUXN = 0d0
!debug            if(i==ii .and. j==jj .and. bin==bb) 
!     &          print *,'L    DU0(',I,J,L,BIN,')   DIF    ',
!     &           'FLUXD  AD44' 
!            if(i==ix .and. j==jx .and. bin==bb) 
!     &          print *,'L    DU0(',I,J,L,BIN,')   DIF    ',
!     &           'FLUXD  AD44' 
!debug-----
            DO L = 1, LLPAR
               DIF = DU0(I,J,L,BIN) - STT(I,J,L,IDTDUST1-1+BIN)

               STT(I,J,L,IDTNK1-1+BIN) = STT(I,J,L,IDTNK1-1+BIN) - 
     &                                DIF/(SQRT( Xk(BIN)*Xk(BIN+1)))
            
               ! Convert flux from [kg/s] to [molec/cm2/s]
               FLUXD = FLUXD + 
     &              DIF / DT_SETTL * XNUMOL(IDTDUST1-1+BIN) / AREA_CM2
               FLUXN = FLUXN + DIF/(SQRT( Xk(BIN)*Xk(BIN+1))) /
     &                      DT_SETTL * XNUMOL(IDTNK1-1+BIN) / AREA_CM2

!debug               if(i==ii .and. j==jj .and. bin==bb) then
!                  print *,L, DU0(I,J,L,BIN), DIF , FLUXD, AD44(I,J,NN,1)
!               endif
!               if(i==ix .and. j==jx .and. bin==bb) then
!                  print *,L, DU0(I,J,L,BIN), DIF , FLUXD, AD44(I,J,NN,1)
!               endif
!debug-----
            ENDDO
         !========================================================      
         ! ND44: Dry deposition diagnostic [#/cm2/s]
         !========================================================
            IF ( ND44 > 0 ) THEN
               AD44(I,J,Input_Opt%IDDEP(BIN),1) = 
     &         AD44(I,J,Input_Opt%IDDEP(BIN),1) + FLUXN
               AD44(I,J,NN,1) = AD44(I,J,NN,1) + FLUXD
            ENDIF

         ENDDO
      ENDDO
      ENDDO

      ! Free pointer memory
      NULLIFY( STT )

      END SUBROUTINE SETTLEDUST
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissdust
!
! !DESCRIPTION: Subroutine EMISSDUST is the driver routine for the dust 
!  emission module.  You may call either the GINOUX or the DEAD dust source 
!  function. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISSDUST( am_I_Root, Input_Opt, 
     &                      State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD             ! ND59
      USE CMN_SIZE_MOD             ! Size parameters
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTDST1
      USE TRACERID_MOD,       ONLY : IDTDST2
      USE TRACERID_MOD,       ONLY : IDTDST3
      USE TRACERID_MOD,       ONLY : IDTDST4
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TRACERID_MOD,       ONLY : IDTDUST1, IDTNK1      !(win, 7/17/09)
      USE TOMAS_MOD,          ONLY : IBINS, XK             !(win, 7/17/09)
      USE DIAG_MOD,           ONLY : AD59_DUST, AD59_NUMB  !(win, 7/17/09)
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
!  30 Mar 2004 - T. D. Fairlie - Initial version
!  (1 ) Now reference LDEAD, LDUST, LPRT from "logical_mod.f".  Now reference!
!        STT from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Add if condition for selecting between emitting 4-bin or 30-bin 
!       dust.  Add emission diagnostic calculation for 30bin dust(win, 7/17/09)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  26 Nov 2012 - R. Yantosca - Now pass am_I_Root, Input_Opt, State_Met as args
!  26 Feb 2013 - R. Yantosca - Now pass Input_Opt to SRC_DUST_GINOUX
!  26 Feb 2013 - R. Yantosca - Changed INPUT_OPT to INTENT(IN), since we are
!                              now no longer calling INIT_DUST from here,
!                              it is now called in the init stage
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL, SAVE :: FIRST = .TRUE.
      
      ! Non-SAVEd scalars
      INTEGER       :: N_TRACERS
      LOGICAL       :: LDEAD
      LOGICAL       :: LDUST
      LOGICAL       :: LPRT
      INTEGER       :: I, J, K                            !(win, 7/17/09)
      REAL*8        :: MEMIS, MINIT(IIPAR,JJPAR,1,IBINS)  !(win, 7/17/09)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:) => NULL()

      !=================================================================
      ! EMISSDUST begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      N_TRACERS = Input_Opt%N_TRACERS
      LDEAD     = Input_Opt%LDEAD
      LDUST     = Input_Opt%LDUST
      LPRT      = Input_Opt%LPRT

      ! Execute on first-call only
      IF ( FIRST ) THEN

         ! Return if dust ID flags are not defined
         IF ( IDTDUST1 == 0 ) THEN
            IF ( LDUST ) THEN 
               CALL ERROR_STOP( 
     &              'LDUST=T but dust tracers are undefined!',
     &              'EMISSDUST ("dust_mod.f")' )
            ENDIF
         ENDIF
 
         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Call appropriate emissions routine
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      !=================================================================
      ! For TOMAS microphysics
      !=================================================================
      IF ( IDTNK1 > 0 .and. IDTDUST1 > 0 ) THEN

         MINIT(:,:,1,1:IBINS) = STT(:,:,1,IDTDUST1:IDTDUST1-1+IBINS)

         IF ( LDEAD ) THEN
            ! still didn't figure out why run would crash w/ this option (win, 7/17/09)
            print *,'Currently the TOMAS code does not work with ',
     &           'dust DEAD emission yet!  Switch to GINOUX for now'
            stop

!            ! Use Zender's DEAD dust source function
!            CALL SRC_DUST_DEAD( STT(:,:,:,IDTDUST1:IDTDUST1-1+IBINS),
!     &                          State_Met )
!            !### Debug
!            IF ( LPRT ) 
!     &           CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_DEAD' )
         ELSE 
!            ! Use Paul Ginoux's dust source function
!            CALL SRC_DUST_GINOUX( STT(:,:,:,IDTDUST1:IDTDUST1-1+IBINS),
     &                            State_Met, Input_Opt )
            !### Debug
            IF ( LPRT ) 
     &           CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_GINOUX')
         ENDIF

         IF ( ND59 > 0 ) THEN 
            DO K = 1, IBINS
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               MEMIS = STT(I,J,1,IDTDUST1-1+K) - MINIT(I,J,1,K)
               IF ( MEMIS == 0.d0 ) GOTO 10

               AD59_DUST(I,J,1,K) = AD59_DUST(I,J,1,K) + MEMIS ! kg ????
               STT(I,J,1,IDTNK1-1+K) = STT(I,J,1,IDTNK1-1+K) +
     &                                 MEMIS / (sqrt(Xk(K)*Xk(K+1)))

               AD59_NUMB(I,J,1,K) = AD59_NUMB(I,J,1,K) + 
     &                              MEMIS / (sqrt(Xk(K)*Xk(K+1)))
 10            CONTINUE
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ENDIF
!#else
!      !=================================================================
!      ! For standard GEOS-Chem simulations w/o TOMAS:
!      !=================================================================
!      IF ( LDEAD ) THEN
!
!         ! Use Zender's DEAD dust source function
!         CALL SRC_DUST_DEAD( STT(:,:,:,IDTDST1:IDTDST4), State_Met )
!
!         !### Debug
!         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_DEAD' )
!         
!      ELSE
!
!         ! Use Paul Ginoux's dust source function
!         CALL SRC_DUST_GINOUX( STT(:,:,:,IDTDST1:IDTDST4 ), 
!     &                         State_Met,
!     &                         Input_Opt )
!
!         !### Debug
!         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_GINOUX')
!
!      ENDIF    
!#endif

      ! Free pointers
      STT => NULL()

      END SUBROUTINE EMISSDUST
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dry_settling 
!
! !DESCRIPTION: Subroutine DRY\_SETTLING computes the dry settling of 
!  dust tracers.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DRY_SETTLING( am_I_Root, Input_Opt, State_Met, TC, RC )
!
! !USES:
!
      USE CMN_GCTM_MOD
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTDST1
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)   :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)   :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)   :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      REAL*8,        INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN) ! Dust [kg]
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)  :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  30 Mar 2004 - T. D. Fairlie - Initial version
!  (1 ) Updated comments, cosmetic changes (bmy, 3/30/04)
!  (2 ) Remove reference to CMN, it's not needed (bmy, 7/20/04)
!  (3 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%XNUMOL(1:N_TRACERS)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: I, J, L, N, N_TRACERS
      REAL*8                :: DT_SETTL, DELZ,  DELZ1
      REAL*8                :: REFF,     DEN,   CONST   
      REAL*8                :: NUM,      LAMDA, FLUX
      REAL*8                :: AREA_CM2, TC0(LLPAR)
      REAL*8                :: TOT1,     TOT2

      ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa      
      REAL*8                :: P 

      ! Diameter of aerosol [um]
      REAL*8                :: Dp

      ! Pressure * DP
      REAL*8                :: PDp 

      ! Temperature (K)    
      REAL*8                :: TEMP        

      ! Slip correction factor
      REAL*8                :: Slip        

      ! Viscosity of air (Pa s)
      REAL*8                :: Visc   

      ! Settling velocity of particle (m/s)
      REAL*8                :: VTS(LLPAR)  

      ! Molecules tracer / kg tracer, from Input Options
      REAL*8                :: XNUMOL(Input_Opt%N_TRACERS)
!
! !DEFINED PARAMETERS:
!      
      REAL*8,  PARAMETER    :: C1 =  0.7674D0
      REAL*8,  PARAMETER    :: C2 =  3.079d0 
      REAL*8,  PARAMETER    :: C3 =  2.573D-11
      REAL*8,  PARAMETER    :: C4 = -1.424d0

      !=================================================================
      ! DRY_SETTLING begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Dust settling timestep [s]
      DT_SETTL  = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,        L,    N,     DEN,  REFF, DP    )
!$OMP+PRIVATE( CONST, AREA_CM2, VTS,  TEMP,  P,    PDP,  SLIP  )
!$OMP+PRIVATE( VISC,  TC0,      DELZ, DELZ1, TOT1, TOT2, FLUX  )

      ! Loop over dust bins
      DO N = 1, NDSTBIN

         ! Initialize
         DEN   = Input_Opt%DUSTDEN(N)
         REFF  = Input_Opt%DUSTREFF(N)
         DP    = 2D0 * REFF * 1.D6              ! Dp [um] = particle diameter
         CONST = 2D0 * DEN * REFF**2 * G0 / 9D0

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, 1 )

               ! Initialize settling velocity
               DO L = 1, LLPAR
                  VTS(L) = 0d0
               ENDDO

               ! Loop over levels
               DO L = 1, LLPAR

                  ! Get P [kPa], T [K], and P*DP
                  P    = GET_PCENTER(I,J,L) * 0.1d0
                  TEMP = State_Met%T(I,J,L)
                  PDP  = P * DP

                  !=====================================================
                  ! # air molecule number density
                  ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
                  !
                  ! # gas mean free path
                  ! lamda = 1.d6 / 
                  !     &   ( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
                  !
                  ! # Slip correction
                  ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * 
                  !      &  exp( -1.1 * Dp / (2. * lamda))) / Dp
                  !=====================================================
                  ! NOTE, Slip correction factor calculations following 
                  !       Seinfeld, pp464 which is thought to be more 
                  !       accurate but more computation required.
                  !=====================================================

                  ! Slip correction factor as function of (P*dp)
                  SLIP = 1d0 + 
     &                   ( 15.60d0 + 7.0d0 * EXP(-0.059d0*PDP) ) / PDP
            
                  !=====================================================
                  ! NOTE, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
                  ! which produce slip correction factor with small 
                  ! error compared to the above with less computation.
                  !=====================================================

                  ! Viscosity [Pa s] of air as a function of temp (K)
                  VISC = 1.458d-6 * (TEMP)**(1.5d0) / ( TEMP + 110.4d0 )

                  ! Settling velocity [m/s]
                  VTS(L) = CONST * SLIP / VISC

               ENDDO

               ! Method is to solve bidiagonal matrix 
               ! which is implicit and first order accurate in Z
               DO L = 1, LLPAR
                  TC0(L) = TC(I,J,L,N)
               ENDDO

               ! We know the boundary condition at the model top
               L           = LLCHEM
               DELZ        = State_Met%BXHEIGHT(I,J,L)
               TC(I,J,L,N) = TC(I,J,L,N) / 
     &                       ( 1.d0 + DT_SETTL * VTS(L) / DELZ )

               DO L = LLCHEM-1, 1, -1
                  DELZ        = State_Met%BXHEIGHT(I,J,L)
                  DELZ1       = State_Met%BXHEIGHT(I,J,L+1)
                  TC(I,J,L,N) = 1.d0 / 
     &                          ( 1.d0 + DT_SETTL * VTS(L)   / DELZ )
     &                 * (TC(I,J,L,N)  + DT_SETTL * VTS(L+1) / DELZ1
     &                 *  TC(I,J,L+1,N) )
               ENDDO

               !========================================================      
               ! ND44: Dry deposition diagnostic [#/cm2/s]
               !========================================================
               IF ( ND44 > 0 ) THEN

                  ! Initialize
                  TOT1 = 0d0
                  TOT2 = 0d0
            
                  ! Compute column totals of TCO(:) and TC(I,J,:,N)
                  DO L = 1, LLPAR
                     TOT1 = TOT1 + TC0(L)
                     TOT2 = TOT2 + TC(I,J,L,N)
                  ENDDO

                  ! Convert dust flux from [kg/s] to [#/cm2/s]
                  FLUX = ( TOT1 - TOT2 ) / DT_SETTL  
                  FLUX = FLUX * XNUMOL(IDTDST1) / AREA_CM2 

                  ! Save in AD44
                  AD44(I,J,Input_Opt%IDDEP(N),1) = 
     &            AD44(I,J,Input_OPt%IDDEP(N),1) + FLUX
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE DRY_SETTLING
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dry_deposition
!
! !DESCRIPTION: Subroutine DRY\_DEPOSITION computes the loss of dust due to 
!  dry deposition at the surface using an implicit method. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DRY_DEPOSITION( am_I_Root, Input_Opt, TC, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTDST1
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      REAL*8,       INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NDSTBIN)  ! Dust [kg]
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  30 Mar 2004 - T. D. Fairlie - Initial version
!  (1 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%XNUMOL(1:N_TRACERS)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER :: I,   J,   L,      N,    N_TRACERS
      REAL*8  :: OLD, NEW, DTCHEM, FLUX, AREA_CM2

      ! Arrays
      REAL*8  :: XNUMOL(Input_Opt%N_TRACERS)

      !=================================================================
      ! DRY_DEPOSITION begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM    = Input_Opt%TS_CHEM * 60d0

      ! Loop over dust bins
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, OLD, NEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over dust bins
      DO N = 1, NDSTBIN    

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, 1 )

               ! Original dust concentration at surface
               OLD = TC(I,J,1,N)

               ! Dust left after dry deposition
               NEW = OLD * 
     &               EXP( -DEPSAV(I,J,Input_Opt%IDDEP(N)) * DTCHEM  )

               !========================================================
               ! ND44 diagnostic: dust drydep loss [#/cm2/s]
               !========================================================
               IF ( ND44 > 0 ) THEN

                  ! Convert drydep flux from [kg/s] to [#/cm2/s]
                  FLUX     = ( OLD - NEW ) / DTCHEM 
                  FLUX     = FLUX * XNUMOL(IDTDST1) / AREA_CM2 
               
                  ! Store in AD44
                  AD44(I,J,Input_Opt%IDDEP(N),1) = 
     &            AD44(I,J,Input_Opt%IDDEP(N),1) + FLUX
               ENDIF

               ! Save back into STT
               TC(I,J,1,N) = NEW 
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE DRY_DEPOSITION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdust_online
!
! !DESCRIPTION: Subroutine RDUST\_ONLINE reads global mineral dust 
!  concentrations as determined by P. Ginoux.  Calculates dust optical 
!  depth at each level for the FAST-J routine "set\_prof.f". 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RDUST_ONLINE( am_I_Root, Input_Opt,  State_Met, 
     &                         DUST,      WAVELENGTH, RC         )
!
! !USES:
!
      USE CMN_FJX_MOD,        ONLY : ODMDUST, IND999, NDUST
      USE CMN_FJX_MOD,        ONLY : QAA, QAA_AOD, RAA, RAA_AOD
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE COMODE_MOD,         ONLY : ERADIUS, IXSAVE, IYSAVE 
      USE COMODE_MOD,         ONLY : IZSAVE,  JLOP,   TAREA
      USE COMODE_MOD,         ONLY : WTAREA, WERADIUS
      USE COMODE_LOOP_MOD
      USE DIAG_MOD,           ONLY : AD21
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D
!
! !INPUT PARAMETERS: 
! 
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
      REAL*8,         INTENT(IN)  :: DUST(IIPAR,JJPAR,LLPAR,NDUST) !Dust [kg/m3]
      INTEGER,        INTENT(IN)  :: WAVELENGTH
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC       ! Success or failure?      
! 
! !REVISION HISTORY: 
!  01 Apr 2004 - R. Martin, R. Park - Initial version
!  (1 ) Bundled into "dust_mod.f" (bmy, 4/1/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f".  Now parallelize over
!        the L-dimension for ND21 diagnostics. (bmy, 7/20/04)
!  (3 ) Archive only hydrophilic aerosol/aqueous dust surface area 
!       (excluding BCPO and OCPO), WTAREA and WERADIUS. (tmf, 3/6/09)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  03 Feb 2011 - S. Kim.     - Include wavelength argument to determine the 
!                              wavelength at which the AOD should be computed. 
!                              This will set the optical properties that are
!                              used for the calculation of the AOD. The ND21 
!                              diagnostic should only be updated when 
!                              WAVELENGTH = 1. (skim, 02/03/11)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, JLOOP, L, N
      REAL*8             :: MSDENS(NDUST), XTAU
      ! Added to calculate aqueous dust surface area (WTAREA, WERADIUS)
      ! (tmf, 3/6/09)
      REAL*8             :: XRH
      REAL*8             :: CRITRH      ! Critical RH [%], above which 
                                        !  heteorogeneous chem takes place

      !=================================================================
      ! RDUST_ONLINE begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Dust density 
      MSDENS(1) = 2500.0d0
      MSDENS(2) = 2500.0d0
      MSDENS(3) = 2500.0d0
      MSDENS(4) = 2500.0d0
      MSDENS(5) = 2650.0d0
      MSDENS(6) = 2650.0d0
      MSDENS(7) = 2650.0d0

      ! Critical RH, above which heteorogeneous chem takes place (tmf, 6/14/07)
      CRITRH = 35.0d0   ! [%]

      !=================================================================     
      ! Convert concentration [kg/m3] to optical depth [unitless].
      !
      ! ODMDUST = ( 0.75 * BXHEIGHT * CONC * QAA ) / 
      !           ( MSDENS * RAA * 1e-6 )
      ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
      !
      !  Units ==> DUST     [ kg/m3    ]
      !            MSDENS   [ kg/m3    ]
      !            RAA      [ um       ]
      !            BXHEIGHT [ m        ]
      !            QAA      [ unitless ]
      !            ODMDUST  [ unitless ]
      !
      ! NOTES: 
      ! (1) Do the calculation at QAA(IND999,:) (i.e. 999 nm).          
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, NDUST
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( WAVELENGTH > 0 ) THEN

            ODMDUST(I,J,L,N) = 0.75d0 * State_Met%BXHEIGHT(I,J,L) * 
     &                         DUST(I,J,L,N) * QAA_AOD(14+N)     / 
     &                        ( MSDENS(N) * RAA_AOD(14+N) * 1.0D-6 )

         ELSE

            ODMDUST(I,J,L,N) = 0.75d0 * State_Met%BXHEIGHT(I,J,L) * 
     &                         DUST(I,J,L,N) * QAA(IND999,14+N)     / 
     &                        ( MSDENS(N) * RAA(IND999,14+N) * 1.0D-6 )

         ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !==============================================================
      ! Calculate Dust Surface Area
      !
      ! Units ==> DUST     [ kg dust/m^3 air    ]
      !           MSDENS   [ kg dust/m^3 dust   ]
      !           RAA      [ um                 ]
      !           TAREA    [ cm^2 dust/cm^3 air ]
      !           ERADIUS  [ cm                 ]
      !
      ! NOTE: first find volume of dust (cm3 dust/cm3 air), then 
      !       multiply by 3/radius to convert to surface area in cm2
      !  
      ! TAREA(:,1:NDUST) and ERADIUS(:,1:NDUST) are for 
      ! the NDUST FAST-J dust wavelength bins (read into DUST)
      !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L, N, XRH )
      DO N     = 1, NDUST
      DO JLOOP = 1, NTTLOOP

         ! Compute 3-D grid box indices
         I = IXSAVE(JLOOP)
         J = IYSAVE(JLOOP)
         L = IZSAVE(JLOOP)

         IF ( WAVELENGTH > 0 ) THEN

            ERADIUS(JLOOP,N) = RAA_AOD(14+N) * 1.0D-4

         ELSE

            ERADIUS(JLOOP,N) = RAA(IND999,14+N) * 1.0D-4

         ENDIF

         TAREA(JLOOP,N)   = 3.D0 / ERADIUS(JLOOP,N) *
     &                      DUST(I,J,L,N) / MSDENS(N)  

         ! Archive WTAREA and WERADIUS when RH > 35%  (tmf, 6/13/07) 
         ! Get RH  
         XRH                = State_Met%RH( I, J, L )  ! [%]
         WTAREA(JLOOP, N)   = 0.d0
         WERADIUS(JLOOP, N) = 0.d0

         IF ( XRH >= CRITRH ) THEN
            WTAREA(JLOOP, N)   = TAREA(JLOOP, N)
            WERADIUS(JLOOP, N) = ERADIUS(JLOOP, N)
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! ND21 Diagnostic: 
      !
      ! Tracer #1: Cloud optical depths    (from "optdepth_mod.f")
      ! Tracer #2: Max Overlap Cld Frac    (from "optdepth_mod.f")
      ! Tracer #3: Random Overlap Cld Frac (from "optdepth_mod.f")
      ! Tracer #4: Dust optical depths (total all size bins)
      ! Tracer #5: Dust surface areas (from all size bins)
      ! Tracers #21-27: Dust AOD for each size bin
      !==============================================================
      IF ( ND21 > 0 .and. WAVELENGTH > 0 ) THEN

         DO N = 1, NDUST

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L )
            DO L = 1, LD21
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               !--------------------------------------
               ! ND21 tracer #4: Dust optical depths (clh)
               !--------------------------------------
               AD21(I,J,L,4) = AD21(I,J,L,4) + ODMDUST(I,J,L,N)

               !--------------------------------------
               ! ND21 tracer #21-27: size-resolved dust optical depths
               !--------------------------------------
               AD21(I,J,L,21+(N-1)) = AD21(I,J,L,21+(N-1)) + 
     &                                ODMDUST(I,J,L,N)
  

               !--------------------------------------
               ! ND21 tracer #5: Dust surface areas
               !--------------------------------------
               IF ( L <= LLCHEM ) THEN
               
                  ! Convert 3-D indices to 1-D index
                  ! JLOP is only defined in the tropopause
                  JLOOP = JLOP(I,J,L)
             
                  ! Add to AD21
                  IF ( JLOOP > 0 ) THEN
                     AD21(I,J,L,5) = AD21(I,J,L,5) + TAREA(JLOOP,N)
                  ENDIF
               ENDIF
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDDO
      ENDIF 

      END SUBROUTINE RDUST_ONLINE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdust_offline
!
! !DESCRIPTION: Subroutine RDUST\_OFFLINE reads global mineral dust 
!  concentrations as determined by P. Ginoux.  Calculates dust optical 
!  depth at each level for the FAST-J routine "set\_prof.f". 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RDUST_OFFLINE( am_I_Root, Input_Opt, State_Met,
     &                          THISMONTH, THISYEAR, WAVELENGTH, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,          ONLY : GET_TAU0,     READ_BPCH2
      USE CMN_FJX_MOD,        ONLY : IND999
      USE CMN_FJX_MOD,        ONLY : ODMDUST, QAA, RAA, RAA_AOD, QAA_AOD
      USE COMODE_MOD,         ONLY : ERADIUS, IXSAVE, IYSAVE
      USE COMODE_MOD,         ONLY : IZSAVE,  JLOP,   TAREA
      USE COMODE_MOD,         ONLY : WTAREA, WERADIUS
      USE COMODE_LOOP_MOD
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD21
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root  ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met  ! Meteorology State object
      INTEGER,        INTENT(IN)  :: THISMONTH  ! Current month (1-12)
      INTEGER,        INTENT(IN)  :: THISYEAR   ! Current year  (YYYY format)
      INTEGER,        INTENT(IN)  :: WAVELENGTH ! Determine which wavelength to
                                                !  use for optical properties
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY: 
!  (1 ) RDUST was patterned after rdaerosol.f (rvm, 9/30/00)
!  (2 ) Don't worry about rewinding the binary file...reading from
!        binary files is pretty fast.  And it's only done once a month.
!  (3 ) Now references punch file utility routines from F90 module
!        "bpch2_mod.f".  Also reference variable DATA_DIR from the
!         header file "CMN_SETUP". (bmy, 9/30/00) 
!  (4 ) Now selects proper GEOS-STRAT dust field for 1996 or 1997.
!        Also need to pass THISYEAR thru the arg list. (rvm, bmy, 11/21/00)
!  (5 ) CONC is now declared as REAL*8 (rvm, bmy, 12/15/00)
!  (6 ) Removed obsolete code from 12/15/00 (bmy, 12/21/00)
!  (7 ) CONC(IIPAR,JJPAR,LGLOB,NDUST) is now CONC(IIPAR,JJPAR,LLPAR,NDUST).
!        Now use routine TRANSFER_3D from "transfer_mod.f" to cast from REAL*4
!        to REAL*8 and also to convert from {IJL}GLOB to IIPAR,JJPAR,LLPAR 
!        space.  Use 3 arguments in call to GET_TAU0.  Updated comments.
!        (bmy, 9/26/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (9 ) Now reference ERADIUS, IXSAVE, IYSAVE, IZSAVE, TAREA from 
!        "comode_mod.f".  Compute ERADIUS and TAREA for the NDUST dust
!        size bins from FAST-J.  Renamed CONC to DUST to avoid conflicts.
!        Also reference NTTLOOP from "comode.h".  Also added parallel
!        DO-loops.  Also renamed MONTH and YEAR to THISMONTH and THISYEAR
!        to avoid conflicts w/ other variables. (bmy, 11/15/01)
!  (10) Bug fix: Make sure to use 1996 dust data for Dec 1995 for the
!        GEOS-STRAT met field dataset.  Set off CASE statement with an
!        #if defined( GEOS_STRAT ) block. (rvm, bmy, 1/2/02)
!  (11) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (12) Now report dust optical depths in ND21 diagnostic at 400 nm.  Now
!       report dust optical depths as one combined diagnostic field instead 
!        of 7 separate fields.  Now reference JLOP from "comode_mod.f".  
!        Now save aerosol surface areas as tracer #5 of the ND21 diagnostic.  
!        (rvm, bmy, 2/28/02)
!  (13) Remove declaration for TIME, since that is also defined in the
!        header file "comode.h" (bmy, 3/20/02)
!  (14) Now read mineral dust files directly from the DATA_DIR/dust_200203/
!        subdirectory (bmy, 4/2/02)
!  (15) Now reference BXHEIGHT from "dao_mod.f".  Also reference ERROR_STOP
!        from "error_mod.f". (bmy, 10/15/02)
!  (16) Now call READ_BPCH2 with QUIET=TRUE to suppress extra informational
!        output from being printed.  Added cosmetic changes. (bmy, 3/14/03)
!  (17) Since December 1997 dust data does not exist, use November 1997 dust
!        data as a proxy. (bnd, bmy, 6/30/03)
!  (18) Bundled into "dust_mod.f" and renamed to RDUST_OFFLINE. (bmy, 4/1/04)
!  (19) Now references DATA_DIR from "directory_mod.f".  Now parallelize over 
!        the L-dimension for ND21 diagnostic. (bmy, 7/20/04)
!  (20) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (21) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (22) Archive only hydrophilic aerosol/aqueous dust surface area 
!       (excluding BCPO and OCPO), WTAREA and WERADIUS. (tmf, 3/6/09)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  03 Feb 2011 - S. Kim      - Include third input argument to determine the 
!                              wavelength at which the AOD should be computed. 
!                              This will set the optical properties that are 
!                              used for the calculation of the AOD.  The ND21 
!                              diagnostic should only be updated when 
!                              WAVELENGTH = 1.
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, JLOOP, L, N
      INTEGER, SAVE       :: MONTH_LAST = -999
      REAL*4              :: TEMP(IIPAR,JJPAR,LGLOB)
      REAL*8              :: DUST(IIPAR,JJPAR,LLPAR,NDUST)
      REAL*8              :: MSDENS(NDUST), XTAU
      CHARACTER (LEN=255) :: FILENAME

      ! Added to calculate aqueous dust surface area (WTAREA, WERADIUS)
      ! (tmf, 3/6/09)
      REAL*8             :: XRH
      REAL*8             :: CRITRH      ! Critical RH [%], above which 
                                        !  heteorogeneous chem takes place

      !=================================================================
      ! RDUST begins here!
      !
      ! Read aerosol data from the binary punch file during the first 
      ! chemistry timestep and, after that, at the start of each month.
      !=================================================================
      IF ( THISMONTH /= MONTH_LAST ) THEN   
         
         ! Save the current month
         MONTH_LAST = THISMONTH

         ! Get TAU0 value used to index the punch file
         ! Use the "generic" year 1985
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
         
         ! Select proper dust file name for GEOS-1, GEOS-3, or GEOS-4
         FILENAME = TRIM( Input_Opt%DATA_DIR ) // 'dust_200203/dust.' //
     &              GET_NAME_EXT()             // '.'                 // 
     &              GET_RES_EXT()

         ! Echo filename
         IF ( am_I_Root ) THEN
            WRITE( 6, 100 ) TRIM( FILENAME )
         ENDIF
 100     FORMAT( '     - RDUST: Reading ', a )

         ! Read aerosol concentrations [kg/m3] for each 
         ! dust type from the binary punch file
         DO N = 1, NDUST 
            CALL READ_BPCH2( FILENAME, 'MDUST-$', N,     XTAU,
     &                       IIPAR,     JJPAR,    LGLOB, TEMP, 
     &                       QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, DUST(:,:,:,N) )
         ENDDO

         !==============================================================
         ! Convert concentration [kg/m3] to optical depth [unitless].
         !
         ! ODMDUST = ( 0.75 * BXHEIGHT * CONC * QAA ) / 
         !           ( MSDENS * RAA * 1e-6 )
         ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
         !
         !  Units ==> DUST     [ kg/m3    ]
         !            MSDENS   [ kg/m3    ]
         !            RAA      [ um       ]
         !            BXHEIGHT [ m        ]
         !            QAA      [ unitless ]
         !            ODMDUST  [ unitless ]
         !
         ! NOTES: 
         ! (1) Do the calculation at QAA(IND999,:) (i.e. 999 nm).          
         !==============================================================
         MSDENS(1) = 2500.0d0
         MSDENS(2) = 2500.0d0
         MSDENS(3) = 2500.0d0
         MSDENS(4) = 2500.0d0
         MSDENS(5) = 2650.0d0
         MSDENS(6) = 2650.0d0
         MSDENS(7) = 2650.0d0

      ! Critical RH, above which heteorogeneous chem takes place (tmf, 6/14/07)
      CRITRH = 35.0d0   ! [%]

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, NDUST
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF ( WAVELENGTH > 0 ) THEN

               ODMDUST(I,J,L,N) = 0.75d0        *
     &                            State_Met%BXHEIGHT(I,J,L) * 
     &                            DUST(I,J,L,N) * QAA_AOD(14+N)     / 
     &                           ( MSDENS(N) * RAA_AOD(14+N) * 1.0D-6 )

            ELSE

               ODMDUST(I,J,L,N) = 0.75d0        *
     &                            State_Met%BXHEIGHT(I,J,L) * 
     &                            DUST(I,J,L,N) * QAA(IND999,14+N)  / 
     &                            ( MSDENS(N) * RAA(IND999,14+N) *
     &                            1.0D-6 )

            ENDIF

         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Echo information
         IF ( am_I_Root ) THEN
            WRITE( 6, 110 ) 
         ENDIF
 110     FORMAT( '     - RDUST: Finished computing optical depths' )

         !==============================================================
         ! Calculate Dust Surface Area
         !
         ! Units ==> DUST     [ kg dust/m^3 air    ]
         !           MSDENS   [ kg dust/m^3 dust   ]
         !           RAA      [ um                 ]
         !           TAREA    [ cm^2 dust/cm^3 air ]
         !           ERADIUS  [ cm                 ]
         !
         ! NOTE: first find volume of dust (cm3 dust/cm3 air), then 
         !       multiply by 3/radius to convert to surface area in cm2
         !  
         ! TAREA(:,1:NDUST) and ERADIUS(:,1:NDUST) are for 
         ! the NDUST FAST-J dust wavelength bins (read into DUST)
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L, N, XRH )
         DO N     = 1, NDUST
         DO JLOOP = 1, NTTLOOP

            ! Compute 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)

            IF ( WAVELENGTH > 0 ) THEN

               ERADIUS(JLOOP,N) = RAA_AOD(14+N) * 1.0D-4

            ELSE

               ERADIUS(JLOOP,N) = RAA(IND999,14+N) * 1.0D-4

            ENDIF

            TAREA(JLOOP,N)   = 3.D0 / ERADIUS(JLOOP,N) *
     &                         DUST(I,J,L,N) / MSDENS(N)  

         ! Archive WTAREA and WERADIUS when RH > 35% (tmf, 6/13/07) 
         ! Get RH  
         XRH                = State_Met%RH( I, J, L )  ! [%]
         WTAREA(JLOOP, N)   = 0.d0
         WERADIUS(JLOOP, N) = 0.d0

         IF ( XRH >= CRITRH ) THEN
            WTAREA(JLOOP, N)   = TAREA(JLOOP, N)
            WERADIUS(JLOOP, N) = ERADIUS(JLOOP, N)
         ENDIF

         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! ND21 Diagnostic: 
         !
         ! Tracer #1: Cloud optical depths    (from "optdepth_mod.f")
         ! Tracer #2: Max Overlap Cld Frac    (from "optdepth_mod.f")
         ! Tracer #3: Random Overlap Cld Frac (from "optdepth_mod.f")
         ! Tracer #4: Dust optical depths at  (from all size bins)
         ! Tracer #5: Dust surface areas (from all size bins)
         !==============================================================
         IF ( ND21 > 0 .and. WAVELENGTH > 0 ) THEN
            DO N = 1, NDUST

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L ) 
               DO L = 1, LD21
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  !--------------------------------------
                  ! ND21 tracer #4: Dust optical depths
                  !--------------------------------------
                  AD21(I,J,L,4) = AD21(I,J,L,4) + ODMDUST(I,J,L,N)


                  !--------------------------------------
                  ! ND21 tracer #21-27: size-resolved dust optical depths
                  !--------------------------------------

                  ! Added the optical depths for each size bin
                  ! (skim, 02/03/11)

                  AD21(I,J,L,21+(N-1)) = AD21(I,J,L,21+(N-1)) + 
     &                                   ODMDUST(I,J,L,N)
  
                  !--------------------------------------
                  ! ND21 tracer #5: Dust surface areas
                  !--------------------------------------
                  IF ( L <= LLCHEM ) THEN

                     ! Convert 3-D indices to 1-D index
                     ! JLOP is only defined in the tropopause
                     JLOOP = JLOP(I,J,L)
             
                     ! Add to AD21
                     IF ( JLOOP > 0 ) THEN
                        AD21(I,J,L,5) = AD21(I,J,L,5) + TAREA(JLOOP,N)
                     ENDIF
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ENDDO
         ENDIF 
      ENDIF

      END SUBROUTINE RDUST_OFFLINE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_dust
!
! !DESCRIPTION: Subroutine INIT\_DUST allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DUST( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DRYDEP_MOD,         ONLY : NTRAIND
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
#if   defined( TOMAS )
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE TOMAS_MOD,          ONLY : IBINS, Xk
      USE TRACERID_MOD,       ONLY : IDTNK1
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  30 Mar 2004 - R. Yantosca - Initial version
!  (1 ) Now references LDEAD from "logical_mod.f" (bmy, 7/20/04)
!  (2 ) Modify to work with 30bin dust. Reference to IBINS from tomas_mod
!       for number of total bi}n = 30 bins.  (win, 7/17/09)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  26 Feb 2013 - M. Long     - Now use fields from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Non-SAVEd scalars
      LOGICAL            :: LDEAD
      INTEGER            :: AS, N
#if   defined( TOMAS )
      INTEGER            :: I
      INTEGER            :: BIN
      REAL*8,  PARAMETER :: PI = 3.14159265d0  !(win, 7/17/09)
      INTEGER            :: ACTMODEBINS        !(scf, 2/27/13)

# if   defined( TOMAS40 )
      ACTMODEBINS = 10  !10 activation mode size bins
# elif defined( TOMAS15 )
      ACTMODEBINS = 3   !3  activation mode size bins
# else
      ACTMODEBINS = 0   !no bins < 10nm
# endif

#endif

      !=================================================================
      ! INIT_DUST begins here!
      !=================================================================

      ! Assume success
      RC    = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LDEAD = Input_Opt%LDEAD

#if defined( TOMAS )
      !=================================================================
      ! INITIALIZATION SECTION FOR TOMAS MICROPHYSICS
      !=================================================================

      !-------------------------------
      ! Set up dust density
      !-------------------------------
      IF ( NDSTBIN == 4 ) THEN
         Input_Opt%DUSTDEN(1:4) = (/ 2500.d0, 2650.d0, 
     &                               2650.d0, 2650.d0 /)

      ELSE IF ( NDSTBIN == IBINS ) THEN

         DO I = 1, IBINS
            IF (Xk(I) < 4.0D-15 ) THEN
               Input_Opt%DUSTDEN(I)  = 2500.d0
            ELSE
               Input_Opt%DUSTDEN(I)  = 2650.d0
            ENDIF
         ENDDO

      ENDIF

      !--------------------------------
      ! Set up dust effective radius
      !--------------------------------
      IF ( NDSTBIN == 4 ) THEN
         Input_Opt%DUSTREFF(1) = 0.73d-6
         Input_Opt%DUSTREFF(2) =  1.4d-6
         Input_Opt%DUSTREFF(3) =  2.4d-6
         Input_Opt%DUSTREFF(4) =  4.5d-6

      ELSE IF ( NDSTBIN == IBINS ) THEN

         ! TOMAS dust Reff (win, 7/17/09)
         DO I = 1, IBINS
            Input_Opt%DUSTREFF(I) = 
     &           0.5d0 * ( SQRT(Xk(I) * Xk(I+1)) /
     &           Input_Opt%DUSTDEN(I) * 6.d0/PI  )**( 0.333d0 )
         ENDDO 
        
      ENDIF

      !----------------------------------
      ! Set up FRAC_S (only for Ginoux)
      !----------------------------------
      IF ( .not. LDEAD ) THEN

         ! Allocate the FRAC_S array
         ALLOCATE( FRAC_S( NDSTBIN ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRAC_S' )
         
         IF ( NDSTBIN == 4 ) THEN

            ! 4 dust bins
            FRAC_S(1) = 0.095d0
            FRAC_S(2) =   0.3d0
            FRAC_S(3) =   0.3d0
            FRAC_S(4) =   0.3d0

         ELSE IF ( NDSTBIN == IBINS ) THEN

            DO I = 1, IBINS
               FRAC_S( I )  = 0.00d-00 !scf initialize
                                       !first few bins in TOMAS15/40 frac_s = 0.0d0
            ENDDO

# if  defined( TOMAS12 ) || defined( TOMAS15 )

            FRAC_S(ACTMODEBINS + 1  )  = 7.33E-10
            FRAC_S(ACTMODEBINS + 2  )  = 2.032E-08
            FRAC_S(ACTMODEBINS + 3  )  = 3.849E-07
            FRAC_S(ACTMODEBINS + 4  )  = 5.01E-06
            FRAC_S(ACTMODEBINS + 5  )  = 4.45E-05
            FRAC_S(ACTMODEBINS + 6  )  = 2.714E-04
            FRAC_S(ACTMODEBINS + 7  )  = 1.133E-03
            FRAC_S(ACTMODEBINS + 8  )  = 3.27E-03
            FRAC_S(ACTMODEBINS + 9  )  = 6.81E-03
            FRAC_S(ACTMODEBINS + 10 )  = 1.276E-02
            FRAC_S(ACTMODEBINS + 11 )  = 2.155E-01
            FRAC_S(ACTMODEBINS + 12 )  = 6.085E-01

# else
! else TOMAS30 or 40

            FRAC_S( ACTMODEBINS +  1 )  = 1.05d-10
            FRAC_S( ACTMODEBINS +  2 )  = 6.28d-10
            FRAC_S( ACTMODEBINS +  3 )  = 3.42d-09
            FRAC_S( ACTMODEBINS +  4 )  = 1.69d-08
            FRAC_S( ACTMODEBINS +  5 )  = 7.59d-08
            FRAC_S( ACTMODEBINS +  6 )  = 3.09d-07
            FRAC_S( ACTMODEBINS +  7 )  = 1.15d-06
            FRAC_S( ACTMODEBINS +  8 )  = 3.86d-06
            FRAC_S( ACTMODEBINS +  9 )  = 1.18d-05
            FRAC_S( ACTMODEBINS + 10 )  = 3.27d-05
            FRAC_S( ACTMODEBINS + 11 )  = 8.24d-05
            FRAC_S( ACTMODEBINS + 12 )  = 1.89d-04
            FRAC_S( ACTMODEBINS + 13 )  = 3.92d-04
            FRAC_S( ACTMODEBINS + 14 )  = 7.41d-04
            FRAC_S( ACTMODEBINS + 15 )  = 1.27d-03
            FRAC_S( ACTMODEBINS + 16 )  = 2.00d-03
            FRAC_S( ACTMODEBINS + 17 )  = 2.89d-03
            FRAC_S( ACTMODEBINS + 18 )  = 3.92d-03
            FRAC_S( ACTMODEBINS + 19 )  = 5.26d-03
            FRAC_S( ACTMODEBINS + 20 )  = 7.50d-03
            FRAC_S( ACTMODEBINS + 21 )  = 1.20d-02
            FRAC_S( ACTMODEBINS + 22 )  = 2.08d-02
            FRAC_S( ACTMODEBINS + 23 )  = 3.62d-02
            FRAC_S( ACTMODEBINS + 24 )  = 5.91d-02
            FRAC_S( ACTMODEBINS + 25 )  = 8.74d-02
            FRAC_S( ACTMODEBINS + 26 )  = 1.15d-01
            FRAC_S( ACTMODEBINS + 27 )  = 1.34d-01
            FRAC_S( ACTMODEBINS + 28 )  = 1.37d-01
            FRAC_S( ACTMODEBINS + 29 )  = 1.24d-01
            FRAC_S( ACTMODEBINS + 30 )  = 9.85d-02

# endif

         ENDIF

         ! Match dry depositing tracer in drydep list with NK1-NK30
         DO BIN = 1, IBINS
         DO N   = 1, Input_Opt%NUMDEP
            IF ( NTRAIND(N) == ( IDTNK1-1+BIN ) )THEN
               Input_Opt%IDDEP(BIN) = N
               GOTO 100
            ENDIF            
         ENDDO
 100     CONTINUE
         ENDDO

      ENDIF

#else

      !=================================================================
      ! INITIALIZATION SECTION FOR STANDARD GEOS-CHEM SIMULATION
      !=================================================================

      ! Dust radii
      Input_Opt%DUSTREFF(1:NDSTBIN) = (/ 0.73d-6, 1.4d-6,  
     &                                   2.4d-6,  4.5d-6  /)

      ! Dust density
      Input_Opt%DUSTDEN(1:NDSTBIN) = (/ 2500.d0, 2650.d0, 
     &                                  2650.d0, 2650.d0 /)

      ! These only have to be allocated for the Ginoux source function
      IF ( .not. LDEAD ) THEN
         ALLOCATE( FRAC_S( NDSTBIN ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRAC_S' )
         FRAC_S(1:NDSTBIN)   = (/ 0.095d0, 0.3d0, 0.3d0, 0.3d0 /)
      ENDIF

      ! Find drydep species in DEPSAV
      IF( Input_Opt%LDRYD ) THEN
         DO N = 1, Input_Opt%NUMDEP
            SELECT CASE ( TRIM( Input_Opt%DEPNAME(N) ) )
               CASE ( 'DST1' )
                  DRYDST1 = N
               CASE ( 'DST2' )
                  DRYDST2 = N
               CASE ( 'DST3' )
                  DRYDST3 = N
               CASE ( 'DST4' )
                  DRYDST4 = N
               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO

         ! This may lead to out of bounds errors
         Input_Opt%IDDEP(1) = DRYDST1
         Input_Opt%IDDEP(2) = DRYDST2
         Input_Opt%IDDEP(3) = DRYDST3
         Input_Opt%IDDEP(4) = DRYDST4
      ENDIF
#endif 

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ) || defined( TOMAS )
      ! EXPERIMENTAL: For archiving the dust source for GIGC
      ALLOCATE( SRCE_FUNC(IIPAR,JJPAR,3), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SRCE_FUNC' )
      SRCE_FUNC = 0.d0
#endif

      END SUBROUTINE INIT_DUST
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_dust
!
! !DESCRIPTION: Subroutine CLEANUP\_DUST deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DUST
! 
! !REVISION HISTORY: 
!  30 Mar 2004 - R. Yantosca - Initial version
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  26 Feb 2013 - R. Yantosca - Now use Input_Opt instead of local arrays
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DUST begins here!
      !=================================================================
      IF ( ALLOCATED( FRAC_S   ) ) DEALLOCATE( FRAC_S   )

      END SUBROUTINE CLEANUP_DUST
!EOC
      END MODULE DUST_MOD
