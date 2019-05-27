!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: dust_mod.F
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
      USE inquireMod,    ONLY : findFreeLUN
      USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

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
      PUBLIC  :: GET_DUST_ALK
      PUBLIC  :: INIT_DUST
      PUBLIC  :: CLEANUP_DUST
!     
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: DRY_SETTLING
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
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  20 Jun 2014 - R. Yantosca - Remove obsolete emissions code; we now use HEMCO
!  13 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  01 Apr 2015 - R. Yantosca - Remove obsolete DUSTMIX, DRY_DEPOSITION routines
!  16 Jun 2016 - J. Kaiser   - Move tracerIDs into variable names
!  20 Jun 2016 - R. Yantosca - Rename IDT* species ID's to id_*
!  20 Jun 2016 - R. Yantosca - Now only define species ID's in the INIT phase
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  17 Mar 2017 - R. Yantosca - Remove reference to Input_Opt fields
!  02 Nov 2017 - R. Yantosca - Add modifications for netCDF diagnostics
!  15 Nov 2017 - E. Lundgren - Implement dust AOD and surf area nc diagnostics
!  03 Jan 2018 - M. Sulprizio- Replace UCX CPP switch with Input_Opt%LUCX
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Species ID flags
      INTEGER               :: id_DST1,  id_DST2,  id_DST3,   id_DST4  
      INTEGER               :: id_DAL1,  id_DAL2,  id_DAL3,   id_DAL4  
      INTEGER               :: id_DUST1, id_NK1

      ! Arrays
      REAL(fp), ALLOCATABLE :: FRAC_S(:)
      REAL(fp), ALLOCATABLE :: SRCE_FUNC(:,:,:)

#if defined( TOMAS )
      ! To replicate the obsolete Input_Opt%IDDEP field
      ! Set to a large placeholder value
      INTEGER               :: TOMAS_IDDEP(1000)
#endif
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
      SUBROUTINE CHEMDUST( am_I_Root, Input_Opt,  State_Met,
     &                     State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Met_Mod,      ONLY : MetState
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
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
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
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
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  01 Apr 2015 - R. Yantosca - Remove call to DRY_DEPOSITION, this is now
!                              done in mixing_mod.F90.
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  10 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  16 Mar 2017 - R. Yantosca - Clean up calls to DRY_SETTLING
!  02 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL, SAVE      :: FIRST = .TRUE.

      ! Non-SAVEd scalars
      LOGICAL            :: LDRYD
      LOGICAL            :: LDUST
      LOGICAL            :: LDSTUP
      LOGICAL            :: prtDebug

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      !=================================================================
      ! CHEMDUST begins here!
      !=================================================================

      ! Assume success
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at CHEMDUST (in module GeosCore/dust_mod.F)'

      ! Copy fields from INPUT_OPT to local variables for use below
      LDRYD  = Input_Opt%LDRYD
      LDUST  = Input_Opt%LDUST
      LDSTUP = Input_Opt%LDSTUP

      ! Set a flag for debug output
      prtDebug = ( Input_Opt%LPrt .and. am_I_Root )

      ! Execute on first call only
      IF ( FIRST ) THEN
 
         ! Stop w/ error if dust species flags are undefined
         ! NOTE: Ind_() returns -1 for species not found, so for the
         ! algorithm below to work, we need to make sure all id's are
         ! not less than zero (bmy, 7/5/16)
         IF ( MAX( id_DST1, 0 )  + 
     &        MAX( id_DST2, 0 )  + 
     &        MAX( id_DST3, 0 )  + 
     &        MAX( id_DST4, 0 ) == 0 ) THEN 
            IF ( LDUST ) THEN 
               ErrMsg = 'LDUST=T but dust species are undefined!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
         ENDIF

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Do dust settling & deposition
      !=================================================================

      !-------------------------------------------------------------
      ! Dust settling for regular dust species
      !
      ! NOTE: Assumes that id_DST1, id_DST2, id_DST3, id_DST4 are 
      ! contiguous in the species list.  This is usually correct,
      ! but we may have to update this later. (bmy, 7/5/16)
      !-------------------------------------------------------------
      CALL DRY_SETTLING( am_I_Root, Input_Opt, State_Met, State_Chm,
     &                   State_Diag, id_DST1,  id_DST4,   RC         )
   
      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN 
         ErrMsg = 'Error encountered in call to "DRY_SETTLING"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !-------------------------------------------------------------
      ! Dust settling for acid uptake on dust species
      !
      ! NOTE: Assumes that id_DAL1, id_DAL2, id_DAL3, id_DAL4 are 
      ! contiguous in the species list.  This is usually correct,
      ! but we may have to update this later. (bmy, 7/5/16)
      !
      ! ALSO NOTE: Dry settling of DSTNIT, DSTSO4 occurs in 
      ! GRAV_SETTLING, called from CHEMSULFATE in SULFATE_MOD.
      !-------------------------------------------------------------
      IF ( LDSTUP ) THEN
         CALL DRY_SETTLING( am_I_Root,  Input_Opt, State_Met, State_Chm,
     &                      State_Diag, id_DAL1,   id_DAL4,   RC       )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN 
            ErrMsg = 'Error encountered in call to "DRY_SETTLING"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

      ! Debug print
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDUST: a DRY_SETTLING' )
      ENDIF

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
      SUBROUTINE SETTLEDUST( am_I_Root, Input_Opt,  State_Met,
     &                       State_Chm, State_Diag, RC )
!
! !USES:
! 
      USE CMN_SIZE_MOD             ! Size parameters
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD             ! ND44
      USE DIAG_MOD,           ONLY : AD44
#endif
      USE ErrCode_Mod
      USE ERROR_MOD
      USE GC_GRID_MOD,        ONLY : GET_AREA_CM2
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : AVO
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TOMAS_MOD,          ONLY : IBINS, Xk, SRTDUST
      USE Species_Mod,        ONLY : Species
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
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
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
!  16 Jul 2015 - R. Yantosca - Now pass INDEX and IDISP to DRY_SETTLING
!  31 May 2016 - E. Lundgren - Replace Input_Opt%XNUMOL with AVO/(emMW_g*1e-3)
!                              where emMW_g is from species database
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  10 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  17 Mar 2017 - R. Yantosca - Update subroutine call to DRY_SETTLING
!  02 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL, SAVE          :: FIRST = .TRUE.

      ! Non-SAVEd scalars
      INTEGER                :: N, BIN, I, J, L, NN, IDISP
      REAL(fp)               :: DU0(IIPAR,JJPAR,LLPAR,IBINS)
      REAL(fp)               :: DIF, FLUXN, FLUXD
      REAL(fp)               :: DT_SETTL, AREA_CM2

      !debug
      integer                :: ii, jj , ix, jx, bb
      data                      ii,jj, ix, jx, bb /37, 24, 58, 34, 30 /
      CHARACTER(LEN=255)     :: MSG, LOC

      ! Pointers
      REAL(fp),      POINTER :: Spc(:,:,:,:)
      TYPE(Species), POINTER :: ThisSpc

      !=================================================================
      ! SETTLEDUST begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS

      ! Check that species units are in [kg] (ewl, 8/13/15)
      IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
         MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
         LOC = 'Routine SETTLEDUST in dust_mod.F'
         CALL GC_Error( MSG, RC, LOC )
      ENDIF

      ! Set pointers
      Spc     => State_Chm%Species  ! [kg], for now
      ThisSpc => NULL()

      !=================================================================
      ! Do dust settling & deposition
      !=================================================================

      ! Dust settling timestep [s]
      DT_SETTL = GET_TS_CHEM()

      ! Save initial dust mass
      DO BIN = 1, IBINS

         ! Species ID
         N = id_DUST1 - 1 + BIN

         DO L   = 1, LLPAR
         DO J   = 1, JJPAR
         DO I   = 1, IIPAR
            DU0(I,J,L,BIN) = Spc(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
      ENDDO

      ! Dust settling
      CALL DRY_SETTLING( am_I_Root,        Input_Opt,
     &                   State_Met,        State_Chm,
     &                   State_Diag,       id_DUST1,
     &                   id_DUST1-1+IBINS, RC         )
   

      ! Calculate change in number to correspond with dust redistribution
      ! by gravitational settling
      DO BIN = 1, IBINS

         ! Species ID
         N       =  id_DUST1 - 1 + BIN

         ! Look up this species in the species database
         ThisSpc => State_Chm%SpcData(N)%Info

         ! Drydep index (??? -- replace w/ species DB info?)
         NN      = State_Chm%nDryDep + (SRTDUST-1)*IBINS + BIN

         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )
         
            FLUXD = 0e+0_fp
            FLUXN = 0e+0_fp
!debug            if(i==ii .and. j==jj .and. bin==bb) 
!     &          print *,'L    DU0(',I,J,L,BIN,')   DIF    ',
!     &           'FLUXD  AD44' 
!            if(i==ix .and. j==jx .and. bin==bb) 
!     &          print *,'L    DU0(',I,J,L,BIN,')   DIF    ',
!     &           'FLUXD  AD44' 
!debug-----
            DO L = 1, LLPAR
               DIF = DU0(I,J,L,BIN) - Spc(I,J,L,id_DUST1-1+BIN)

               Spc(I,J,L,id_NK1-1+BIN) = Spc(I,J,L,id_NK1-1+BIN) - 
     &                                DIF/(SQRT( Xk(BIN)*Xk(BIN+1)))
            
               ! Convert flux from [kg/s] to [molec/cm2/s]
               FLUXD = FLUXD +                                                 
     &                 DIF / DT_SETTL * AVO / 
     &                 ( 1.e-3_fp * ThisSpc%emMW_g ) /        
     &                 AREA_CM2                                                
           
               
               FLUXN = FLUXN + DIF/(SQRT( Xk(BIN)*Xk(BIN+1))) /
     &                 DT_SETTL * AVO / 
     &                 ( 1.e-3_fp * ThisSpc%emMW_g ) / 
     &                 AREA_CM2

!debug               if(i==ii .and. j==jj .and. bin==bb) then
!                  print *,L, DU0(I,J,L,BIN), DIF , FLUXD, AD44(I,J,NN,1)
!               endif
!               if(i==ix .and. j==jx .and. bin==bb) then
!                  print *,L, DU0(I,J,L,BIN), DIF , FLUXD, AD44(I,J,NN,1)
!               endif
!debug-----
            ENDDO

#if defined( BPCH_DIAG )
            !=========================================================== 
            !  Dry deposition diagnostic [#/cm2/s] (bpch)
            !===========================================================
            IF ( ND44 > 0 ) THEN

               !--------------------------------------------------------
               ! ND44 DIAGNOSTIC (bpch) 
               ! Dry deposition flux loss [#/cm2/s]
               !
               ! NOTE: Bpch diagnostics are being phased out.
               !--------------------------------------------------------
               !%%% NOTE: Now use TOMAS_IDDEP, which replicates the
               !%%% since-removed Input_Opt%IDDEP field (bmy, 3/17/17)
               AD44(I,J,TOMAS_IDDEP(BIN),1) = 
     &         AD44(I,J,TOMAS_IDDEP(BIN),1) + FLUXN
               AD44(I,J,NN,1) = AD44(I,J,NN,1) + FLUXD
            ENDIF
#endif

            ! Add to State_Diag later when TOMAS diagnostics are implemented

         ENDDO
         ENDDO

         ! Free pointer
         ThisSpc => NULL()

      ENDDO

      ! Free pointer memory
      Spc => NULL()

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
      SUBROUTINE EMISSDUST( am_I_Root, Input_Opt,  State_Met,
     &                      State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD             ! Size parameters
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD             ! ND59
      USE DIAG_MOD,           ONLY : AD59_DUST, AD59_NUMB  !(win, 7/17/09)
#endif
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
      USE State_Met_Mod,      ONLY : MetState
      USE TOMAS_MOD,          ONLY : IBINS, XK             !(win, 7/17/09)
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
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
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
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  02 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL,  SAVE    :: FIRST = .TRUE.
      
      ! Non-SAVEd scalars
      LOGICAL           :: LDEAD
      LOGICAL           :: LDUST
      LOGICAL           :: LPRT
      LOGICAL           :: LINTERP 
      INTEGER           :: I, J, K
      REAL(fp)          :: MEMIS, MINIT(IIPAR,JJPAR,1,IBINS)

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      !=================================================================
      ! EMISSDUST begins here!
      !=================================================================
     
      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at EMISSDUST (in module GeosCore/dust_mod.F)'

      ! Initialize
      LDEAD     =  Input_Opt%LDEAD
      LDUST     =  Input_Opt%LDUST
      LPRT      =  Input_Opt%LPRT

      ! Execute on first-call only
      IF ( FIRST ) THEN

         ! Return if dust ID flags are not defined
         IF ( id_DUST1 == 0 ) THEN
            IF ( LDUST ) THEN 
               ErrMsg = 'LDUST=T but dust species are undefined!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
         ENDIF
 
         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Call appropriate emissions routine
      !=================================================================

      ! Point to chemical species array [kg]
      Spc => State_Chm%Species

      !=================================================================
      ! For TOMAS microphysics
      !=================================================================
      IF ( id_NK1 > 0 .and. id_DUST1 > 0 ) THEN

         MINIT(:,:,1,1:IBINS) = Spc(:,:,1,id_DUST1:id_DUST1-1+IBINS)

         IF ( LDEAD ) THEN
            ! still didn't figure out why run would crash w/ this option (win, 7/17/09)
            print *,'Currently the TOMAS code does not work with ',
     &           'dust DEAD emission yet!  Switch to GINOUX for now'
            stop

         ELSE 

            !### Debug
            IF ( LPRT ) 
     &           CALL DEBUG_MSG( '### EMISSDUST: a SRC_DUST_GINOUX')
         ENDIF

#if defined( BPCH_DIAG )
         IF ( ND59 > 0 ) THEN 
            DO K = 1, IBINS
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               MEMIS = Spc(I,J,1,id_DUST1-1+K) - MINIT(I,J,1,K)
               IF ( MEMIS == 0.e+0_fp ) GOTO 10

               AD59_DUST(I,J,1,K) = AD59_DUST(I,J,1,K) + MEMIS ! kg ????
               Spc(I,J,1,id_NK1-1+K) = Spc(I,J,1,id_NK1-1+K) +
     &                                 MEMIS / (sqrt(Xk(K)*Xk(K+1)))

               AD59_NUMB(I,J,1,K) = AD59_NUMB(I,J,1,K) + 
     &                              MEMIS / (sqrt(Xk(K)*Xk(K+1)))
 10            CONTINUE
            ENDDO
            ENDDO
            ENDDO
         ENDIF
#endif
      ENDIF

      ! Free pointers
      Spc => NULL()

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
!  dust species.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DRY_SETTLING( am_I_Root, Input_Opt,  State_Met,
     &                         State_Chm, State_Diag, Ind0,
     &                         Ind1,      RC                     )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD44
#endif
      USE ErrCode_Mod
      USE GC_GRID_MOD,        ONLY : GET_AREA_CM2
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      INTEGER,        INTENT(IN)    :: Ind0        ! Starting dust species #
      INTEGER,        INTENT(IN)    :: Ind1        ! Ending   dust species #
!
! !INPUT/OUTPUT PARAMETERS
! 
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
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
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  26 Feb 2015 - E. Lundgren - Replace GET_PCENTER with State_Met%PMID and
!                              remove dependency on pressure_mod
!  25 Jan 2016 - E. Lundgren - Update netcdf drydep flux diagnostic
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  31 May 2016 - E. Lundgren - Use TCVV instead of XNUMOL for molecular weights
!  22 Jun 2016 - M. Yannetti - Replace TCVV with species db MW and phys constant
!  16 Mar 2017 - R. Yantosca - Remove TC, IDISP; Loop over slots Ind0:Ind1
!                              of the State_Chm%Species array.
!  02 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I,        J,     L
      INTEGER                :: N,        NA,    ND
      REAL(fp)               :: DT_SETTL, DELZ,  DELZ1
      REAL(fp)               :: REFF,     DEN,   CONST   
      REAL(fp)               :: NUM,      LAMDA, FLUX
      REAL(fp)               :: AREA_CM2, TC0(LLPAR)
      REAL(fp)               :: TOT1,     TOT2

      ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa      
      REAL(fp)               :: P 

      ! Diameter of aerosol [um]
      REAL(fp)               :: Dp

      ! Pressure * DP
      REAL(fp)               :: PDp 

      ! Temperature (K)    
      REAL(fp)               :: TEMP        

      ! Slip correction factor
      REAL(fp)               :: Slip        

      ! Viscosity of air (Pa s)
      REAL(fp)               :: Visc   

      ! Settling velocity of particle (m/s)
      REAL(fp)               :: VTS(LLPAR)  

      ! Pointers
      TYPE(Species), POINTER :: ThisSpc
      REAL(fp),      POINTER :: TC(:,:,:,:)

!
! !DEFINED PARAMETERS:
!      
      REAL(fp),  PARAMETER   :: C1 =  0.7674e+0_fp
      REAL(fp),  PARAMETER   :: C2 =  3.079e+0_fp 
      REAL(fp),  PARAMETER   :: C3 =  2.573e-11_fp
      REAL(fp),  PARAMETER   :: C4 = -1.424e+0_fp 

      !=================================================================
      ! DRY_SETTLING begins here!
      !=================================================================

      ! Initialize
      RC        =  GC_SUCCESS
      ThisSpc   => NULL()
      TC        => State_Chm%Species

      ! Dust settling timestep [s]
      DT_SETTL  = GET_TS_CHEM()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,        L,    N,     DEN,  REFF, DP    )
!$OMP+PRIVATE( CONST, AREA_CM2, VTS,  TEMP,  P,    PDP,  SLIP  )
!$OMP+PRIVATE( VISC,  TC0,      DELZ, DELZ1, TOT1, TOT2, FLUX  )
!$OMP+PRIVATE( NA,    ThisSpc,  ND                             )

      ! Loop over only the advected dust species
      DO NA = Ind0, Ind1

         ! Look up this species in the species database
         ! NOTE: The advected species are listed first in the master
         ! species list, so it's OK to use the advected species ID
         ! to query the Species Database (bmy, 3/16/17)
         ThisSpc => State_Chm%SpcData(NA)%Info
         
         ! Get the drydep ID corresponding to this species
         ND      =  ThisSpc%DryDepId

         ! Density [kg/m3] and radius [m]
         DEN     =  ThisSpc%Density
         REFF    =  ThisSpc%Radius

         ! Dp [um] = particle diameter
         DP      =  2e+0_fp * REFF * 1.e+6_fp    
         CONST   =  2e+0_fp * DEN * REFF**2 * g0 / 9e+0_fp

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, 1 )

               ! Initialize settling velocity
               DO L = 1, LLPAR
                  VTS(L) = 0e+0_fp
               ENDDO

               ! Loop over levels
               DO L = 1, LLPAR

                  ! Get P [kPa], T [K], and P*DP
                  ! Use moist air pressure for mean free path (ewl, 3/2/15)
                  P    = State_Met%PMID(I,J,L) * 0.1e+0_fp
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
                  SLIP = 1e+0_fp + ( 15.60e+0_fp + 7.0e+0_fp * 
     &                   EXP(-0.059e+0_fp*PDP) ) / PDP
            
                  !=====================================================
                  ! NOTE, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
                  ! which produce slip correction factor with small 
                  ! error compared to the above with less computation.
                  !=====================================================

                  ! Viscosity [Pa s] of air as a function of temp (K)
                  VISC = 1.458e-6_fp * (TEMP)**(1.5e+0_fp) / 
     &                 ( TEMP + 110.4e+0_fp )

                  ! Settling velocity [m/s]
                  VTS(L) = CONST * SLIP / VISC

               ENDDO

               ! Method is to solve bidiagonal matrix 
               ! which is implicit and first order accurate in Z
               DO L = 1, LLPAR
                  TC0(L) = TC(I,J,L,NA)
               ENDDO

               ! We know the boundary condition at the model top
               L           = LLCHEM
               DELZ        = State_Met%BXHEIGHT(I,J,L)
               TC(I,J,L,NA) = TC(I,J,L,NA) / 
     &                        ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )

               DO L = LLCHEM-1, 1, -1
                  DELZ         = State_Met%BXHEIGHT(I,J,L)
                  DELZ1        = State_Met%BXHEIGHT(I,J,L+1)
                  TC(I,J,L,NA) = 1.e+0_fp / 
     &                 ( 1.e+0_fp + DT_SETTL * VTS(L)   / DELZ )
     &                 * (TC(I,J,L,NA)  + DT_SETTL * VTS(L+1) / DELZ1
     &                 *  TC(I,J,L+1,NA) )
               ENDDO

               !========================================================      
               ! Dry deposition diagnostic [molec/cm2/s]
               !========================================================
#if defined( BPCH_DIAG )
               !-----------------------------------------------------
               ! ND44 DIAGNOSTIC (bpch)
               ! Dry deposition flux loss [molec/cm2/s]
               !-----------------------------------------------------
               IF ( ND44 > 0 ) THEN

                  ! Initialize
                  TOT1 = 0e+0_fp
                  TOT2 = 0e+0_fp
                  
                  ! Compute column totals of TCO(:) and TC(I,J,:,N)
                  DO L = 1, LLPAR
                     TOT1 = TOT1 + TC0(L)
                     TOT2 = TOT2 + TC(I,J,L,NA)
                  ENDDO
                  
                  ! Convert dust flux from [kg/s] to [molec/cm2/s]
                  FLUX = ( TOT1 - TOT2 ) / DT_SETTL  
                  FLUX = FLUX * AVO * ( AIRMW / 
     &                    ThisSpc%emMw_g    ) /
     &                   ( AIRMW * 1.e-3_fp ) / AREA_CM2 
                  
                  ! For bpch diagnostic, save in global array AD44
                  AD44(I,J,ND,1) = AD44(I,J,ND,1) + FLUX
               ENDIF
#endif

               IF ( State_Diag%Archive_DryDepChm .OR.
     &              State_Diag%Archive_DryDep        ) THEN
                  !-----------------------------------------------------
                  ! HISTORY (aka netCDF diagnostics)
                  ! Dry deposition flux loss [molec/cm2/s]
                  !
                  ! NOTE: Eventually think about converting this
                  ! diagnostic to more standard units [kg/m2/s]
                  !-----------------------------------------------------
                  ! Initialize
                  TOT1 = 0e+0_fp
                  TOT2 = 0e+0_fp
            
                  ! Compute column totals of TCO(:) and TC(I,J,:,N)
                  DO L = 1, LLPAR
                     TOT1 = TOT1 + TC0(L)
                     TOT2 = TOT2 + TC(I,J,L,NA)
                  ENDDO

                  ! Convert dust flux from [kg/s] to [molec/cm2/s]
                  FLUX = ( TOT1 - TOT2 ) / DT_SETTL  
                  FLUX = FLUX * AVO * ( AIRMW / 
     &                    ThisSpc%emMw_g    ) /
     &                   ( AIRMW * 1.e-3_fp ) / AREA_CM2 

                  ! Drydep flux in chemistry only
                  State_Diag%DryDepChm(I,J,ND) = FLUX
               ENDIF
            ENDDO 
         ENDDO    

         ! Nullify pointer
         ThisSpc => NULL()
      ENDDO       
!$OMP END PARALLEL DO

      ! Free pointers
      ThisSpc => NULL()
      TC      => NULL()

      END SUBROUTINE DRY_SETTLING
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
     &                         State_Chm, State_Diag, DUST,
     &                         ODSWITCH,  RC                     )
!
! !USES:
!
      USE CMN_FJX_MOD
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG ) || defined( MODEL_GEOS )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD21
#endif
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D
!
! !INPUT PARAMETERS: 
! 
      LOGICAL,        INTENT(IN)  :: am_I_Root     ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt     ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met     ! Meteorology State object
      REAL(fp),       INTENT(IN)  :: DUST(IIPAR,JJPAR,LLPAR,NDUST) !Dust [kg/m3]
      INTEGER,        INTENT(IN)  :: ODSWITCH
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC            ! Success or failure?
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
!  04 Sep 2012 - D. Ridley   - WAVELENGTH now ODSWITCH for clarity now
!                              that multiple wavelengths can be
!                              calculated at once.
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  12 May 2016 - M. Sulprizio- Remove 1D arrays that depend on JLOOP. ERADIUS,
!                              TAREA, WERADIUS, WTAREA are now pointers that
!                              point to 3D fields in State_Chm.
!  16 May 2016 - M. Sulprizio- Remove JLOOP entirely and loop over LLPAR, JJPAR,
!                              IIPAR instead.
!  27 Mar 2017 - R. Yantosca - Fixed bug in formula for NOUT diagnostic index
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  15 Nov 2017 - E. Lundgren - Implement dust AOD and surf area nc diagnostics
!  23 Mar 2018 - E. Lundgren - Get number of wavelengths from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: LINTERP
      INTEGER            :: I, J, L, N, NOUT, W
      INTEGER            :: IWV, IIWV, NWVS, IDST
      REAL(fp)           :: MSDENS(NDUST), XTAU

      ! Added to calculate aqueous dust surface area (WTAREA, WERADIUS)
      ! (tmf, 3/6/09)
      REAL(fp)           :: XRH
      REAL(fp)           :: CRITRH      ! Critical RH [%], above which 
                                        !  heteorogeneous chem takes place

      INTEGER            :: IWV1,  IWV2
      REAL(fp)           :: ACOEF, BCOEF, LOGTERM

      ! Pointers
      REAL(fp), POINTER   :: ERADIUS(:,:,:,:)
      REAL(fp), POINTER   :: TAREA(:,:,:,:)
      REAL(fp), POINTER   :: WERADIUS(:,:,:,:)
      REAL(fp), POINTER   :: WTAREA(:,:,:,:)

      ! For diagnostics
      REAL(fp)            :: tempOD(IIPAR,JJPAR,LLPAR,NDUST,3)
      LOGICAL             :: LINTERPARR(Input_Opt%NWVSELECT)

      !=================================================================
      ! RDUST_ONLINE begins here!
      !=================================================================

      ! Assume success
      RC        =  GC_SUCCESS

      ! Initialize pointers
      ERADIUS   => State_Chm%AeroRadi     ! Aerosol Radius     [cm]
      TAREA     => State_Chm%AeroArea     ! Aerosol Area       [cm2/cm3]
      WERADIUS  => State_Chm%WetAeroRadi  ! Wet Aerosol Radius [cm]
      WTAREA    => State_Chm%WetAeroArea  ! Wet Aerosol Area   [cm2/cm3]

      ! Index for dust in ODAER and LUT arrays
      ! Dust properties are saved to different indices in RD_AOD for
      ! UCX vs tropchem simulations
      IF ( Input_Opt%LUCX ) THEN
         IDST   = 8
      ELSE
         IDST   = 6
      ENDIF

      ! Dust density 
      MSDENS(1) = 2500.0e+0_fp
      MSDENS(2) = 2500.0e+0_fp
      MSDENS(3) = 2500.0e+0_fp
      MSDENS(4) = 2500.0e+0_fp
      MSDENS(5) = 2650.0e+0_fp
      MSDENS(6) = 2650.0e+0_fp
      MSDENS(7) = 2650.0e+0_fp

      ! Critical RH, above which heteorogeneous chem takes place (tmf, 6/14/07)
      CRITRH = 35.0e+0_fp   ! [%]

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
      ! DAR Oct 2012
      ! if the call is from chemistry (ODSWITCH=0) only need one wavelength
      ! if radiation on (LRAD=TRUE) then cycle over all LUT wavelengths
      ! if radiation off but call is for AOD then cycle over the
      ! wavelength required or the wavelengths between which to
      ! interpolate

      IF (ODSWITCH .EQ. 0) THEN
         NWVS   = 1
      ELSE
         IF ( Input_Opt%LRAD ) THEN !Loop over all RT wavelengths (30)
            ! plus any required for calculating the AOD
            NWVS   = NWVAART+NWVREQUIRED
         ELSE                       !Loop over wavelengths needed (from RD_AOD)
            NWVS   = NWVREQUIRED
         ENDIF
      ENDIF
      
      DO IIWV = 1, NWVS
         ! get current wavelength index
         ! (specified by user or all wavelengths for RRTMG)
         ! (1000nm always used for FAST-J and stored in first index)
         IF (ODSWITCH .EQ. 0) THEN
            ! only doing for 1000nm i.e. IWV=10 in LUT
            ! N.B. NWVS is fixed to 1 above - only one wavelength
            IWV=IWV1000
         ELSE
            IF ( Input_Opt%LRAD ) THEN
            ! RRTMG wavelengths begin after NWVAA0 standard wavelengths
            ! but add on any others required
               IF (IIWV.LE.NWVAART) THEN
                  !index of RRTMG wavelengths starts after the standard NWVAA0
                  !(currently NWVAA0=11, set in CMN_FJX_mod based on the .dat
                  !LUT)
                  IWV = IIWV+NWVAA0
               ELSE
                  !now we calculate at wvs for the requested AOD
                  !offset index by NWVAART i.e. start from 1 
                  IWV = IWVREQUIRED(IIWV-NWVAART)
               ENDIF
            ELSE
               ! IWVREQUIRED lists the index of requires standard wavelengths
               IWV = IWVREQUIRED(IIWV)
            ENDIF
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, NDUST
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! dust stored in the IDST species bin of LUT variables
            ODMDUST(I,J,L,IWV,N) = 0.75e+0_fp *
     &                             State_Met%BXHEIGHT(I,J,L) * 
     &                             DUST(I,J,L,N) * QQAA(IWV,N,IDST)  / 
     &                           ( MSDENS(N) * RDAA(N,IDST) * 1.0e-6_fp)

#if defined ( RRTMG )
            !add dust optics to the RT code arrays
            !SSA and ASYM copying seems a little redundant...
            !will keep this way for uniformity for now but
            !possibly could deal with SSA and ASYM in RT module
            RTODAER(I,J,L,IWV,NAER+2+N) = ODMDUST(I,J,L,IWV,N)
            RTSSAER(I,J,L,IWV,NAER+2+N) = SSAA(IWV,N,IDST)
            RTASYMAER(I,J,L,IWV,NAER+2+N) = ASYMAA(IWV,N,IDST)
#endif

         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDDO !wavelength loop

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
!$OMP+PRIVATE( I, J, L, N, XRH )
      DO N = 1, NDUST
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip non-chemistry boxes
         IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

         ERADIUS(I,J,L,N) = RDAA(N,IDST) * 1.0e-4_fp
         
         TAREA(I,J,L,N)   = 3.e+0_fp / ERADIUS(I,J,L,N) *
     &                      DUST(I,J,L,N) / MSDENS(N)  

         ! Archive WTAREA and WERADIUS when RH > 35%  (tmf, 6/13/07) 
         ! Get RH  
         XRH                = State_Met%RH( I, J, L )  ! [%]
         WTAREA(I,J,L, N)   = 0.e+0_fp
         WERADIUS(I,J,L, N) = 0.e+0_fp

         IF ( XRH >= CRITRH ) THEN
            WTAREA(I,J,L, N)   = TAREA(I,J,L, N)
            WERADIUS(I,J,L, N) = ERADIUS(I,J,L, N)
         ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#if defined( BPCH_DIAG ) || defined( MODEL_GEOS )
      !=================================================================
      ! ND21 Diagnostic: 
      !
      ! Tracer #1: Cloud optical depths    (from "optdepth_mod.f")
      ! Tracer #2: Max Overlap Cld Frac    (from "optdepth_mod.f")
      ! Tracer #3: Random Overlap Cld Frac (from "optdepth_mod.f")
      ! Tracer #4: Dust optical depths (total all size bins)
      ! Tracer #5: Dust surface areas (from all size bins)
      ! Tracers #31-37: Dust AOD for each size bin
      ! Tracers #38-44: Dust SSA for each size bin
      ! Tracers #45-51: Dust ASYM for each size bin
      !==============================================================
      IF ( Input_Opt%ND21 > 0 .and. ODSWITCH .EQ. 1 ) THEN

!#############################################################################
!### Prior to 4/6/15
!### NOTE: Activating this parallel loop causes differences in the ND21
!### diagnostic outputs for some reason.  We have not been able to diagnose
!### this.  The cheap solution is to disable the parallelization, and then
!### we get identical results between single and multiprocessor runs.
!### (bmy, 4/6/15)
!###
!###!$OMP PARALLEL DO
!###!$OMP+DEFAULT( SHARED )
!###!$OMP+PRIVATE( I, J, L, N, W, NOUT, LINTERP )
!#############################################################################
         DO N=1, NDUST
         DO W = 1, Input_Opt%NWVSELECT !loop over number of wavelengths
            IF (IWVSELECT(1,W).EQ.IWVSELECT(2,W)) THEN
               LINTERP=.FALSE.
            ELSE
               LINTERP=.TRUE.
            ENDIF
            NOUT = 5 + (NRHAER*3) + ((NRHAER+NDUST)*(W-1)) + N
            ! equivalent to getting 21-27 in the loop for W = 1
            ! etc for W = 2 and W=3 (fills OD dust bins 1-7 for each lambda)

            DO L = 1, Input_Opt%LD21
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               IF ( .not. LINTERP ) THEN

                  !--------------------------------------
                  ! ND21 tracer #4: Dust optical depths (clh)
                  !--------------------------------------
                  AD21(I,J,L,4) = AD21(I,J,L,4) + 
     &                            ODMDUST(I,J,L,IWVSELECT(1,W),N)

                  !--------------------------------------
                  ! ND21 tracer #21-27: size-resolved dust optical depths
                  !--------------------------------------
                  AD21(I,J,L,NOUT) = AD21(I,J,L,NOUT) + 
     &                               ODMDUST(I,J,L,IWVSELECT(1,W),N)

               ELSE

                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  ! AOD sometimes zero (if Q zero), must catch this
                  IF ((ODMDUST(I,J,L,IWVSELECT(1,W),N).GT.0).AND.
     &                (ODMDUST(I,J,L,IWVSELECT(2,W),N).GT.0)) THEN

                     !--------------------------------------
                     ! ND21 tracer #4: Dust optical depths (clh)
                     !--------------------------------------
                     AD21(I,J,L,4) = AD21(I,J,L,4) +
     &                             ODMDUST(I,J,L,IWVSELECT(2,W),N)*
     &                             ACOEF_WV(W)**(BCOEF_WV(W)*LOG(
     &                             ODMDUST(I,J,L,IWVSELECT(1,W),N)/
     &                             ODMDUST(I,J,L,IWVSELECT(2,W),N)))

                     AD21(I,J,L,NOUT) = AD21(I,J,L,NOUT) +
     &                                ODMDUST(I,J,L,IWVSELECT(2,W),N)*
     &                                ACOEF_WV(W)**(BCOEF_WV(W)*LOG(
     &                                ODMDUST(I,J,L,IWVSELECT(1,W),N)/
     &                                ODMDUST(I,J,L,IWVSELECT(2,W),N)))


                  ENDIF
               ENDIF

               !--------------------------------------
               ! ND21 tracer #5: Dust surface areas
               !--------------------------------------
               IF ( L <= LLCHEM ) THEN
                           
                  ! Add to AD21
                  AD21(I,J,L,5) = AD21(I,J,L,5) + TAREA(I,J,L,N)

               ENDIF
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
!#############################################################################
!### Prior to 4/6/15
!### NOTE: Activating this parallel loop causes differences in the ND21
!### diagnostic outputs for some reason.  We have not been able to diagnose
!### this.  The cheap solution is to disable the parallelization, and then
!### we get identical results between single and multiprocessor runs.
!### (bmy, 4/6/15)
!###
!###!$OMP END PARALLEL DO
!#############################################################################
      ENDIF
#endif

      IF ( State_Diag%Archive_AOD .AND. ODSWITCH .EQ. 1 ) THEN
         
         ! Initialize local arrays
         tempOD = 0.0_fp
         LINTERPARR(:) = .FALSE.

         ! Interpolate?
         DO W = 1, Input_Opt%NWVSELECT
            IF (IWVSELECT(1,W).EQ.IWVSELECT(2,W)) THEN
               LINTERPARR(W) =.FALSE.
            ELSE
               LINTERPARR(W) =.TRUE.
            ENDIF
         ENDDO

         ! Loop over dust bins, # of wavelengths, and all grid cells
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, W, NOUT, LINTERP )
         DO W = 1, Input_Opt%NWVSELECT 
         DO N = 1, NDUST
         DO L = 1, LLCHEM
         DO J = 1, JJPAR
         DO I = 1, IIPAR               

            IF ( .not. LINTERPARR(W) ) THEN

               tempOD(I,J,L,N,W) = ODMDUST(I,J,L,IWVSELECT(1,W),N)
            
            ELSEIF ((ODMDUST(I,J,L,IWVSELECT(1,W),N).GT.0).AND.
     &             (ODMDUST(I,J,L,IWVSELECT(2,W),N).GT.0)) THEN
            
               ! Interpolate using angstrom exponent between
               ! closest available wavelengths (coefs pre-calculated 
               ! in CALC_AOD)
               ! AOD sometimes zero (if Q zero), must catch this!
               tempOD(I,J,L,N,W) = 
     &                  ODMDUST(I,J,L,IWVSELECT(2,W),N)*
     &                  ACOEF_WV(W)**(BCOEF_WV(W)*LOG(
     &                  ODMDUST(I,J,L,IWVSELECT(1,W),N)/
     &                  ODMDUST(I,J,L,IWVSELECT(2,W),N)))

            ENDIF

            !---------------------------------------------------
            ! Set size-resolved dust optical depth diagnostic
            !---------------------------------------------------
            IF ( State_Diag%Archive_AODDustWL1 .AND. W == 1 ) THEN
               State_Diag%AODDustWL1(I,J,L,N) = tempOD(I,J,L,N,1)
            ENDIF
            IF ( State_Diag%Archive_AODDustWL2 .AND. W == 2 ) THEN
               State_Diag%AODDustWL2(I,J,L,N) = tempOD(I,J,L,N,2)
            ENDIF
            IF ( State_Diag%Archive_AODDustWL3 .AND. W == 3 ) THEN
               State_Diag%AODDustWL3(I,J,L,N) = tempOD(I,J,L,N,3)
            ENDIF


         ENDDO ! longitude
         ENDDO ! latitude
         ENDDO ! level
         ENDDO ! bin
         ENDDO ! wavelength
!$OMP END PARALLEL DO

         !---------------------------------------------------
         ! Set dust optical depth diagnostic
         !---------------------------------------------------
         IF ( State_Diag%Archive_AODDust ) THEN
!$OMP PARALLEL DO

!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLCHEM
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               State_Diag%AODDust(I,J,L) = SUM( tempOD(I,J,L,:,:) )
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF

      ENDIF

      ! Archive total dust surface area (sum across all bins)
      IF ( State_Diag%Archive_AerSurfAreaDust ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLCHEM
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            State_Diag%AerSurfAreaDust(I,J,L) =  
     &                                  SUM( TAREA(I,J,L,1:NDUST) )
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Free pointers
      ERADIUS  => NULL()
      TAREA    => NULL()
      WERADIUS => NULL()
      WTAREA   => NULL()

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
      SUBROUTINE RDUST_OFFLINE( am_I_Root, Input_Opt,  State_Met,
     &                          State_Chm, State_Diag, THISMONTH, 
     &                          THISYEAR,  ODSWITCH,   RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_FJX_MOD
#if defined( BPCH_DIAG ) || defined( MODEL_GEOS )
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,     READ_BPCH2
      USE CMN_DIAG_MOD
      USE DIAG_MOD,       ONLY : AD21
#endif
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Diag_Mod, ONLY : DgnState
      USE State_Met_Mod,  ONLY : MetState
      USE TRANSFER_MOD,   ONLY : TRANSFER_3D

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root  ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
      INTEGER,        INTENT(IN)    :: THISMONTH  ! Current month (1-12)
      INTEGER,        INTENT(IN)    :: THISYEAR   ! Current year  (YYYY format)
      INTEGER,        INTENT(IN)    :: ODSWITCH   ! Determine which wavelength 
                                                  !  to use for optical 
                                                  !  properties
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
!
! !REMARKS:
!  ##########################################################################
!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!  #####    WILL EVENTUALLY BE READ IN FROM netCDF FILES VIA HEMCO!     #####
!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!  ##########################################################################
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
!  (5 ) CONC is now declared as REAL(fp) (rvm, bmy, 12/15/00)
!  (6 ) Removed obsolete code from 12/15/00 (bmy, 12/21/00)
!  (7 ) CONC(IIPAR,JJPAR,LGLOB,NDUST) is now CONC(IIPAR,JJPAR,LLPAR,NDUST).
!        Now use routine TRANSFER_3D from "transfer_mod.f" to cast from REAL*4
!        to REAL(fp) and also to convert from {IJL}GLOB to IIPAR,JJPAR,LLPAR 
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
!  04 Sep 2012 - D. Ridley   - WAVELENGTH now ODSWITCH for clarity now
!                              that multiple wavelengths can be
!                              calculated at once.
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path
!  23 Oct 2015 - E. Lundgren - Fix index for dust in ODAER/LUT arrays if UCX on
!  12 May 2016 - M. Sulprizio- Remove 1D arrays that depend on JLOOP. ERADIUS,
!                              TAREA, WERADIUS, WTAREA are now pointers that
!                              point to 3D fields in State_Chm.
!  16 May 2016 - M. Sulprizio- Remove JLOOP entirely and loop over LLPAR, JJPAR,
!                              IIPAR instead.
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  25 Jan 2018 - E. Lundgren - Implement dust AOD and surf area nc diagnostics
!  23 Mar 2018 - E. Lundgren - Get number of wavelengths from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL             :: LINTERP
      INTEGER             :: I, J, L, N, NOUT, W
      INTEGER             :: IWV, IIWV, NWVS, IDST
      INTEGER, SAVE       :: MONTH_LAST = -999
      LOGICAL, SAVE       :: DO_DIAGN = .FALSE.
      REAL(f4)            :: TEMP(IIPAR,JJPAR,LGLOB)
      REAL(fp)            :: DUST(IIPAR,JJPAR,LLPAR,NDUST)
      REAL(fp)            :: MSDENS(NDUST)
      REAL(f8)            :: XTAU
      CHARACTER (LEN=255) :: FILENAME

      ! Added to calculate aqueous dust surface area (WTAREA, WERADIUS)
      ! (tmf, 3/6/09)
      REAL(fp)            :: XRH
      REAL(fp)            :: CRITRH      ! Critical RH [%], above which 
                                         !  heteorogeneous chem takes place

      ! Pointers
      REAL(fp), POINTER   :: ERADIUS(:,:,:,:)
      REAL(fp), POINTER   :: TAREA(:,:,:,:)
      REAL(fp), POINTER   :: WERADIUS(:,:,:,:)
      REAL(fp), POINTER   :: WTAREA(:,:,:,:)

      ! For diagnostics
      REAL(fp)            :: tempOD(IIPAR,JJPAR,LLPAR,NDUST,3)
      LOGICAL             :: LINTERPARR(Input_Opt%NWVSELECT)

      !=================================================================
      ! RDUST begins here!
      !=================================================================

      ! Assume success
      RC        =  GC_SUCCESS

#if defined( BPCH_DIAG )
      !#################################################################
      !####  WE HAVE BLOCKED OUT THIS CODE SO AS TO AVOID COMPILING
      !####  BPCH SPECIFIC CODE IN GCHP.  THIS MAY BREAK SOME THINGS
      !####  GOING FORWARD.  WORRY ABOUT IT LATER.  BPCH IS BEING
      !####  TOTALLY PHASED OUT (bmy, 1/12/18)
      !#################################################################
 

      ! Initialize pointers
      ERADIUS   => State_Chm%AeroRadi     ! Aerosol Radius     [cm]
      TAREA     => State_Chm%AeroArea     ! Aerosol Area       [cm2/cm3]
      WERADIUS  => State_Chm%WetAeroRadi  ! Wet Aerosol Radius [cm]
      WTAREA    => State_Chm%WetAeroArea  ! Wet Aerosol Area   [cm2/cm3]

      ! Index for dust in ODAER and LUT arrays
      ! Dust properties are saved to different indices in RD_AOD for
      ! UCX vs tropchem simulations
      IF ( Input_Opt%LUCX ) THEN
         IDST   = 8
      ELSE
         IDST   = 6
      ENDIF

      ! Read aerosol data from the binary punch file during the first 
      ! chemistry timestep and, after that, at the start of each month.
      IF ( THISMONTH /= MONTH_LAST .or. DO_DIAGN ) THEN
         
         ! Get TAU0 value used to index the punch file
         ! Use the "generic" year 1985
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
         
         ! Select proper dust file name for GEOS-1, GEOS-3, or GEOS-4
         FILENAME = TRIM( Input_Opt%DATA_DIR ) // 
     &              TRIM( Input_Opt%RES_DIR  ) // 'dust_200203/dust.' //
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
      MSDENS(1) = 2500.0e+0_fp
      MSDENS(2) = 2500.0e+0_fp
      MSDENS(3) = 2500.0e+0_fp
      MSDENS(4) = 2500.0e+0_fp
      MSDENS(5) = 2650.0e+0_fp
      MSDENS(6) = 2650.0e+0_fp
      MSDENS(7) = 2650.0e+0_fp

      ! Critical RH, above which heteorogeneous chem takes place (tmf, 6/14/07)
      CRITRH    = 35.0e+0_fp  ! [%]

      ! DAR 09/2013
      ! There are two ways RDUST_ can be called:
      ! (1) When Fast-J requires aerosol optics at 1000nm (ODSWITCH=0)
      ! (2) Before diags are accumulated, from RECOMPUTE_AOD (ODSWITCH=1)
      ! for (1) we just need the optics stored at a single wavelength
      !     not for all user specified wavelengths, hence NWVS=1, IWV=IWV1000
      ! for (2) we need to determine if RRTMG is switched on
      !     if LRAD=true, calculation is over total minus standard wavelengths
      !     in optics dat files (NWVAA-NWVAA0)
      !     if LRAD=false, calculation is for the wavelengths required for
      !     user-requested wavelength output. These are determined in CALC_AOD
      !     within RD_AOD.F and stored in IWVREQUIRED. The coefficients to
      !     interpolate from the LUT wavelengths to the user-requested
      !     waveelenths (in CALC_AOD) are used here.

      ! Select number of wavelengths required to loop over
      IF (ODSWITCH .EQ. 0) THEN
         NWVS   = 1
      ELSE
         IF ( Input_Opt%LRAD ) THEN !Loop over all RT wavelengths (30)
            ! plus any required for calculating the AOD
            NWVS   = NWVAA-NWVAA0+NWVREQUIRED
         ELSE                       !Loop over wavelengths needed (from RD_AOD)
            NWVS   = NWVREQUIRED
         ENDIF
      ENDIF

      DO IIWV = 1, NWVS
         !now select the correct LUT wavelength
         IF (ODSWITCH .EQ. 0) THEN
            ! only doing for 1000nm (IWV1000 is set in RD_AOD)
            ! N.B. NWVS is fixed to 1 above - only one wavelength
            IWV=IWV1000
         ELSE
            IF ( Input_Opt%LRAD ) THEN
               ! RRTMG wavelengths begin after NWVAA0 standard wavelengths
               ! but add on any others required
               IF (IIWV.LE.30) THEN
                  !index of RRTMG wavelengths starts after the standard NWVAA0
                  !(currently NWVAA0=11, set in CMN_FJX_mod based on the .dat
                  !LUT)
                  IWV = IIWV+NWVAA0
               ELSE
                  !now we calculate at wvs for the requested AOD
                  IWV = IWVREQUIRED(IIWV-30)
               ENDIF
            ELSE
               ! IWVREQUIRED lists the index of requires standard wavelengths
               IWV = IWVREQUIRED(IIWV)
            ENDIF
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, NDUST
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! dust stored in the IDST species bin of LUT variables
            ODMDUST(I,J,L,IWV,N) = 0.75e+0_fp *
     &                             State_Met%BXHEIGHT(I,J,L) *
     &                             DUST(I,J,L,N) * QQAA(IWV,N,IDST)  /
     &                           ( MSDENS(N) * RDAA(N,IDST) * 
     &                             1.0e-6_fp )

#if defined ( RRTMG )
            !add dust optics to the RT code arrays
            !SSA and ASYM copying seems a little redundant...
            !will keep this way for uniformity for now but
            !possibly could deal with SSA and ASYM in RT module
            RTODAER(I,J,L,IWV,NAER+4+N) = ODMDUST(I,J,L,IWV,N)
            RTSSAER(I,J,L,IWV,NAER+4+N) = SSAA(IWV,N,IDST)
            RTASYMAER(I,J,L,IWV,NAER+4+N) = ASYMAA(IWV,N,IDST)
#endif

         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDDO !wavelength loop

      ! Echo information
      IF ( am_I_Root ) THEN
         WRITE( 6, 110 )
      ENDIF
 110  FORMAT( '     - RDUST: Finished computing optical depths' )

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
!$OMP+PRIVATE( I, J, L, N, XRH )
      DO N = 1, NDUST
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip non-chemistry boxes
         IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

         ERADIUS(I,J,L,N) = RDAA(N,6) * 1.0e-4_fp

         TAREA(I,J,L,N)   = 3.e+0_fp / ERADIUS(I,J,L,N) *
     &                      DUST(I,J,L,N) / MSDENS(N)  

         ! Archive WTAREA and WERADIUS when RH > 35% (tmf, 6/13/07) 
         ! Get RH  
         XRH               = State_Met%RH(I,J,L)  ! [%]
         WTAREA(I,J,L,N)   = 0.e+0_fp
         WERADIUS(I,J,L,N) = 0.e+0_fp

         IF ( XRH >= CRITRH ) THEN
            WTAREA(I,J,L,N)   = TAREA(I,J,L,N)
            WERADIUS(I,J,L,N) = ERADIUS(I,J,L,N)
         ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#if defined( BPCH_DIAG ) || defined( MODEL_GEOS )
      !==============================================================
      ! ND21 Diagnostic: 
      !
      ! Tracer #1: Cloud optical depths    (from "optdepth_mod.f")
      ! Tracer #2: Max Overlap Cld Frac    (from "optdepth_mod.f")
      ! Tracer #3: Random Overlap Cld Frac (from "optdepth_mod.f")
      ! Tracer #4: Dust optical depths at  (from all size bins)
      ! Tracer #5: Dust surface areas (from all size bins)
      !==============================================================
      IF ( ND21 > 0 .and. ODSWITCH .EQ. 1 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, W, NOUT, LINTERP )
         DO N=1, NDUST
         DO W = 1, Input_Opt%NWVSELECT !different to NWVSELECT if interp
            IF (IWVSELECT(1,W).EQ.IWVSELECT(2,W)) THEN
               LINTERP=.FALSE.
            ELSE
               LINTERP=.TRUE.
            ENDIF
            SELECT CASE(W)
               CASE(1)
                  NOUT=6+3*NAER+(N-1)
               CASE(2)
                  NOUT=6+3*NAER+NAER+NDUST+(N-1)
               CASE(3)
                  NOUT=6+3*NAER+2*(NAER+NDUST)+(N-1)
               CASE DEFAULT
                  NOUT=6+3*NAER+(N-1)
            END SELECT

            DO L = 1, LD21
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               IF ( .not. LINTERP ) THEN

                  !--------------------------------------
                  ! ND21 tracer #4: Dust optical depths (clh, dar)
                  !--------------------------------------
                  AD21(I,J,L,4) = AD21(I,J,L,4) +
     &                            ODMDUST(I,J,L,IWVSELECT(1,W),N)

                  !--------------------------------------
                  ! ND21 tracer #21-27: size-resolved dust optical depths
                  !--------------------------------------
                  AD21(I,J,L,NOUT) = AD21(I,J,L,NOUT) +
     &                               ODMDUST(I,J,L,IWVSELECT(1,W),N)

               ELSE

                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  ! AOD sometimes zero (if Q zero), must catch this
                  IF ((ODMDUST(I,J,L,IWVSELECT(1,W),N).EQ.0).OR.
     &                (ODMDUST(I,J,L,IWVSELECT(2,W),N).EQ.0)) THEN
                     !Dont add anything to diagnostic
                  ELSE

                     !--------------------------------------
                     ! ND21 tracer #4: Dust optical depths (clh, dar)
                     !--------------------------------------
                     AD21(I,J,L,4) = AD21(I,J,L,4) +
     &                             ODMDUST(I,J,L,IWVSELECT(2,W),N)*
     &                             ACOEF_WV(W)**(BCOEF_WV(W)*LOG(
     &                             ODMDUST(I,J,L,IWVSELECT(1,W),N)/
     &                             ODMDUST(I,J,L,IWVSELECT(2,W),N)))

                     AD21(I,J,L,NOUT) = AD21(I,J,L,NOUT) +
     &                                ODMDUST(I,J,L,IWVSELECT(2,W),N)*
     &                                ACOEF_WV(W)**(BCOEF_WV(W)*LOG(
     &                                ODMDUST(I,J,L,IWVSELECT(1,W),N)/
     &                                ODMDUST(I,J,L,IWVSELECT(2,W),N)))

                  ENDIF
               ENDIF

               !--------------------------------------
               ! ND21 tracer #5: Dust surface areas
               !--------------------------------------
               IF ( L <= LLTROP ) THEN

                  ! Add to AD21
                  AD21(I,J,L,5) = AD21(I,J,L,5) + TAREA(I,J,L,N)

               ENDIF
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
#endif

      IF ( State_Diag%Archive_AOD .AND. ODSWITCH .EQ. 1 ) THEN
         
         ! Initialize local arrays
         tempOD = 0.0_fp
         LINTERPARR(:) = .FALSE.

         ! Interpolate?
         DO W = 1, Input_Opt%NWVSELECT
            IF (IWVSELECT(1,W).EQ.IWVSELECT(2,W)) THEN
               LINTERPARR(W) =.FALSE.
            ELSE
               LINTERPARR(W) =.TRUE.
            ENDIF
         ENDDO

         ! Loop over dust bins, # of wavelengths, and all grid cells
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, W, NOUT, LINTERP )
         DO W = 1, Input_Opt%NWVSELECT
         DO N = 1, NDUST
         DO L = 1, LLCHEM
         DO J = 1, JJPAR
         DO I = 1, IIPAR               
            IF ( .not. LINTERPARR(W) ) THEN
               tempOD(I,J,L,W,N) = ODMDUST(I,J,L,IWVSELECT(1,W),N)
            ELSEIF ((ODMDUST(I,J,L,IWVSELECT(1,W),N).GT.0).AND.
     &             (ODMDUST(I,J,L,IWVSELECT(2,W),N).GT.0)) THEN
                !Dont add anything to diagnostic
            ELSE
               ! Interpolate using angstrom exponent between
               ! closest available wavelengths (coefs pre-calculated 
               ! in CALC_AOD)
               ! AOD sometimes zero (if Q zero), must catch this!
               tempOD(I,J,L,W,N) = 
     &                  ODMDUST(I,J,L,IWVSELECT(2,W),N)*
     &                  ACOEF_WV(W)**(BCOEF_WV(W)*LOG(
     &                  ODMDUST(I,J,L,IWVSELECT(1,W),N)/
     &                  ODMDUST(I,J,L,IWVSELECT(2,W),N)))
            ENDIF

            !---------------------------------------------------
            ! Set size-resolved dust optical depth diagnostic
            !---------------------------------------------------
            IF ( State_Diag%Archive_AODDustWL1 .AND. W == 1 ) THEN
               State_Diag%AODDustWL1(I,J,L,N) = tempOD(I,J,L,N,1)
            ENDIF
            IF ( State_Diag%Archive_AODDustWL2 .AND. W == 2 ) THEN
               State_Diag%AODDustWL2(I,J,L,N) = tempOD(I,J,L,N,2)
            ENDIF
            IF ( State_Diag%Archive_AODDustWL3 .AND. W == 3 ) THEN
               State_Diag%AODDustWL3(I,J,L,N) = tempOD(I,J,L,N,3)
            ENDIF

         ENDDO ! longitude
         ENDDO ! latitude
         ENDDO ! level
         ENDDO ! bin
         ENDDO ! wavelength
!$OMP END PARALLEL DO

         !---------------------------------------------------
         ! Set dust optical depth diagnostic
         !---------------------------------------------------
         IF ( State_Diag%Archive_AODDust ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLCHEM
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               State_Diag%AODDust(I,J,L) = SUM( tempOD(I,J,L,:,:) )
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF

      ENDIF

      !------------------------------------------------------
      ! Archive total dust surface area (sum across all bins)
      !------------------------------------------------------
      IF ( State_Diag%Archive_AerSurfAreaDust ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLTROP
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            State_Diag%AerSurfAreaDust(I,J,L) =  
     &                                  SUM( TAREA(I,J,L,1:NDUST) )
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

         ! Save the current month
         ! Set diagnostic to aggregate next call time. 
         IF ( THISMONTH /= MONTH_LAST ) THEN
            MONTH_LAST = THISMONTH
            DO_DIAGN = .TRUE.
         ELSE
            DO_DIAGN = .FALSE.
         ENDIF
         
      ENDIF
      
      ! Free pointers
      ERADIUS  => NULL()
      TAREA    => NULL()
      WERADIUS => NULL()
      WTAREA   => NULL()

#endif

      END SUBROUTINE RDUST_OFFLINE
!EOC
!-----------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_dust_alk
!
! !DESCRIPTION: Subroutine GET\_DUST\_ALK returns: (1) dust alkalinity,
!  ALK\_d(NDSTBIN) [v/v], (2) rate coefficients, KTS(NDSTBIN), KTN(NDSTBIN),
!  for uptake of SO2 and HNO3 on dust for use in sulfate\_mod.f for chemistry
!  on dust aerosols, (3) fraction, KTH(NDSTBIN), of the size-weighted total
!  area of aerosols in the grid box. GET\_DUST\_ALK is analogous to GET\_ALK
!  for seasalt (bec, 12/7/04; tdf 04/08/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_DUST_ALK( I, J, L, ALK_d, KTS, KTN, KTH,
     &                         Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE CMN_SIZE_MOD             ! Size parameters
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : PI, AIRMW
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)    :: I, J, L        ! Grid box indices
      TYPE(OptInput), INTENT(IN)    :: Input_Opt      ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met      ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm      ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(OUT)   :: ALK_d(NDSTBIN) ! Dust alkalinity [v/v]
      REAL(fp),       INTENT(OUT)   :: KTS  (NDSTBIN) ! Rate coef for uptake of
                                                      ! SO2 on dust [s-1]
      REAL(fp),       INTENT(OUT)   :: KTN  (NDSTBIN) ! Rate coef for uptake of
                                                      ! HNO3 on dust [s-1]
      REAL(fp),       INTENT(OUT)   :: KTH  (NDSTBIN) ! Fraction of the size-
                                                      ! weighted total area
                                                      ! of aerosols in grid box
!
! !REVISION HISTORY:
!  08 Apr 2008 - T.D. Fairlie- Initial version
!  16 Sep 2013 - M. Sulprizio- Added ProTeX headers
!  17 Sep 2013 - M. Sulprizio- Now pass Input_Opt, State_Met, and State_Chm
!                              as arguments
!  14 Nov 2013 - M. Sulprizio- Bug fix: Avoid div-by-zero in calculation of
!                              gas-to-particle rate constants for SO2 and HNO3
!  05 Jan 2016 - E. Lundgren - Use global PI
!  12 May 2016 - M. Sulprizio- Remove 1D arrays that depend on JLOOP. ERADIUS,
!                              TAREA, are now pointers that point to 3D fields
!                              fields in State_Chm.
!  22 Jun 2016 - M. Yannetti - Replace TCVV with species db MW and phys constant
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER :: MINDAT = 1.d-20

      ! Need this for dust
      !REAL(fp), PARAMETER :: GAMMA_SO2 = 0.05d0  !(Song & Carmichael, 2001)
      ! 200 times smaller 8/28/2K9
      REAL(fp), PARAMETER :: GAMMA_SO2 = 2.5d-4

      !tdf V9 4/1/2K9 Applying Song et al.(2007) reduced value 
      REAL(fp), PARAMETER :: GAMMA_H2SO4 = 1.d0  

      !  Need this for dust
      !REAL(fp), PARAMETER :: GAMMA_HNO3 = 0.1d0 ! (Song & Carmichael, 2001)
      ! 200 times smaller 8/28/2K9
      REAL(fp), PARAMETER :: GAMMA_HNO3 = 5.0d-4 

      REAL(fp), PARAMETER :: DG = 0.2d0 ! gas phase diffusion coeff. [cm2/s]
      REAL(fp), PARAMETER :: v = 3.0d4  ! cm/s
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER             :: IRH
      INTEGER             :: IBIN, ISBIN
      REAL(fp)            :: N1, KT1, KT1N
      REAL(fp)            :: AREA, HGF
      REAL(fp)            :: CONST1, CONST2, CONST
      REAL(fp)            :: AIR_DENS
      REAL(fp)            :: A1, B1, A1N, B1N, T1, R1
      REAL(fp)            :: DD, RD, DM
      REAL(fp)            :: TOTAL_AREA, DF_TOTAL_AREA
      REAL(fp)            :: DST_d (NDSTBIN), ALK, GAMS, GAMN
      REAL(fp)            :: SULF_AREA, BC_AREA, OC_AREA
      REAL(fp)            :: SULF_RAD, BC_RAD, OC_RAD
      REAL(fp)            :: SULF_FAC, BC_FAC, OC_FAC
      REAL(fp)            :: SSA_AREA, SSC_AREA
      REAL(fp)            :: SSA_RAD, SSC_RAD
      REAL(fp)            :: SSA_FAC, SSC_FAC
      LOGICAL, SAVE       :: FIRST = .TRUE.

      ! Arrays
      ! Dust Surface Areas                                ! tdf 08/20/09
      REAL(fp)            :: AREA_d(NDSTBIN)              ! [cm^2/cm^3]

      ! Dust Surface Areas within sub-bins 1-4 of BIN 1   ! tdf 08/20/09
      REAL(fp)            :: AREA_sd1(4)                  ! [cm^2/cm^3]

      ! Dust Effective Radius                             ! tdf 08/20/09
      REAL(fp)            :: RD_d(NDSTBIN)                ! [cm]

      ! Dust Effective Radii for sub-bins 1-4 of BIN 1    ! tdf 08/20/09
      REAL(fp)            :: RD_sd1(4)                    ! [cm]

      ! Dust size-weighted Surface Areas                  ! tdf 08/20/09
      REAL(fp)            :: DF_AREA_d(NDSTBIN)           ! [1/s]

      ! Dust size-weighted Surface Areas for sub-bins 1-4 ! tdf 08/20/09
      REAL(fp)            :: DF_AREA_sd1(4)               ! [1/s]

      ! Molecular weights
      REAL(fp)            :: MW_DST1, MW_DST2, MW_DST3, MW_DST4

      ! Pointers
      REAL(fp), POINTER   :: Spc(:,:,:,:)
      REAL(fp), POINTER   :: ERADIUS(:,:,:,:)
      REAL(fp), POINTER   :: TAREA(:,:,:,:)

      !=================================================================
      ! GET_DUST_ALK begins here!
      !=================================================================

      ! Initialize pointers
      Spc      => State_Chm%Species  ! GEOS-Chem species array [v/v dry]
      ERADIUS  => State_Chm%AeroRadi ! Aerosol Radius [cm]
      TAREA    => State_Chm%AeroArea ! Aerosol Area [cm2/cm3]

      ! Get MWs from species database 
      MW_DST1 = State_Chm%SpcData(id_DST1)%Info%emMW_g
      MW_DST2 = State_Chm%SpcData(id_DST2)%Info%emMW_g
      MW_DST3 = State_Chm%SpcData(id_DST3)%Info%emMW_g
      MW_DST4 = State_Chm%SpcData(id_DST4)%Info%emMW_g
     
      ! Zero variables
      DO IBIN = 1, NDSTBIN
         ALK_d (IBIN)  = 0.0e+0_fp
         KTS   (IBIN)  = 0.0e+0_fp
         KTN   (IBIN)  = 0.0e+0_fp
         KTH   (IBIN)  = 0.0e+0_fp
         AREA_d (IBIN) = 0.0e+0_fp
         RD_d (IBIN)   = 0.0e+0_fp
      END DO

      ! Air density [kg/m3]
      AIR_DENS = State_Met%AD(I,J,L) / State_Met%AIRVOL(I,J,L) 

      ! Retrieve Dust Alkalinity [v/v dry from Spc array
      ALK_d(1) = Spc(I,J,L,id_DAL1)
      ALK_d(2) = Spc(I,J,L,id_DAL2)
      ALK_d(3) = Spc(I,J,L,id_DAL3)
      ALK_d(4) = Spc(I,J,L,id_DAL4)

      ! Dust [kg/m3] from Spc, used to compute dust surface area
      ! Units: (moles/mole).(kg(air)/m3).(kg(dust)/mole)/(kg(air)/mole)
      DST_d(1) = Spc(I,J,L,id_DST1) * AIR_DENS / ( AIRMW / MW_DST1 )
      DST_d(2) = Spc(I,J,L,id_DST2) * AIR_DENS / ( AIRMW / MW_DST2 )
      DST_d(3) = Spc(I,J,L,id_DST3) * AIR_DENS / ( AIRMW / MW_DST3 )
      DST_d(4) = Spc(I,J,L,id_DST4) * AIR_DENS / ( AIRMW / MW_DST4 )

      ! tdf Now get aerosol surface area from TAREA (cm2/cm3)
      SULF_AREA = TAREA(I,J,L,NDUST+1)
      BC_AREA   = TAREA(I,J,L,NDUST+2)
      OC_AREA   = TAREA(I,J,L,NDUST+3)
      SSA_AREA  = TAREA(I,J,L,NDUST+4)
      SSC_AREA  = TAREA(I,J,L,NDUST+5)

      ! tdf Now get aerosol effective radius from ERADIUS (cm)
      SULF_RAD  = ERADIUS(I,J,L,NDUST+1)
      BC_RAD    = ERADIUS(I,J,L,NDUST+2)
      OC_RAD    = ERADIUS(I,J,L,NDUST+3)
      SSA_RAD   = ERADIUS(I,J,L,NDUST+4)
      SSC_RAD   = ERADIUS(I,J,L,NDUST+5)

      ! tdf Quotients [s/cm] used to weight surface area for H2SO4 uptake
      SULF_FAC = (SULF_RAD / DG + 4.e+0_fp/(V*GAMMA_H2SO4) )
      BC_FAC   = (  BC_RAD / DG + 4.e+0_fp/(V*GAMMA_H2SO4) )
      OC_FAC   = (  OC_RAD / DG + 4.e+0_fp/(V*GAMMA_H2SO4) )
      SSA_FAC  = ( SSA_RAD / DG + 4.e+0_fp/(V*GAMMA_H2SO4) )
      SSC_FAC  = ( SSC_RAD / DG + 4.e+0_fp/(V*GAMMA_H2SO4) )

      !tdf Surface areas and effective radii for sub-bins 1-4 of dust bin 1
      DO ISBIN = 1, 4
         T1 = TAREA  (I,J,L,ISBIN)
         R1 = ERADIUS(I,J,L,ISBIN)
         AREA_sd1    (ISBIN) = T1
         RD_sd1      (ISBIN) = R1
         !tdf surface area for sub bins 1-4 in bin 1, weighted by gas-phase
         !tdf diffusion and collision limitations
         !tdf used to compute proportionate uptake of H2SO4 only  [1/s]
         DF_AREA_sd1 (ISBIN) = T1
     &                       / (R1/DG + 4.0e+0_fp/(V*GAMMA_H2SO4))
      END DO

      !-----------------------------------------------------------------------
      ! Very Simple Formulation: For each size bin (i)   ! tdf 8/20/09
      ! Dust Area density = 3 * Dust Mass density  / (REFF(i) * DUSTDEN)
      ! TAREA computed   in RDUST_ONLINE - Units: cm^2(dust) / cm^3(air) 
      ! ERADIUS computed in RDUST_ONLINE - Units: cm
      ! NB: I am now subdividing the submicron dust size bin 
      !     using TAREA (I,J,L,1->4), and ERADIUS (I,J,L,1->4).
      !-----------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Find Dust surface area density in grid-box, AREA_d [cm^2/cm^3].
      !  Also find the size-weighted surface area density, DF_AREA_d [1/s].
      !  The latter represents the gas-phase diffusion and surface 
      !  limited weighting and is used to determine the fraction of H2SO4 
      !  taken up on dust, versus taken up on other aerosols.
      !                                                 tdf 08/21/09
      !-------------------------------------------------------------------------

      ! tdf Loop over size bins  (NDSTBIN = 4)
      DO IBIN = 1, NDSTBIN

         ! Dust Area density in grid box,      AREA_d [cm^2/cm^3]    tdf 8/21/09
         ! Dust weighted surface area density, DF_AREA_d [1/s]       tdf 8/21/09

         IF (IBIN .EQ. 1) THEN
            ! For Dust size bin 1, sum over the 4 size sub bins  tdf 8/21/09
!tdf 11/25/2013            AREA_d   (IBIN) = AREA_sd1(1) + AREA_sd1(1)        ! [cm^2/cm^3]
            AREA_d   (IBIN) = AREA_sd1(1) + AREA_sd1(2)        ![cm^2/cm^3]
     &                      + AREA_sd1(3) + AREA_sd1(4) 
            DF_AREA_d(IBIN) = DF_AREA_sd1(1) + DF_AREA_sd1(3)  ! [1/s]
     &                      + DF_AREA_sd1(2) + DF_AREA_sd1(4)
         ELSE
            T1 = TAREA(I,J,L,3+IBIN)      ! [cm^2/cm^3]
            R1 = ERADIUS(I,J,L,3+IBIN)    ! [cm]
!tdf 11/25/2013
            RD_d     (IBIN) = R1
!tdf 11/25/2013
            AREA_d   (IBIN) = T1          ! [cm^2/cm^3]
            DF_AREA_d(IBIN) = T1          ! [1/s]
     &                       / (R1/DG + 4.0D0/(V*GAMMA_H2SO4))
         ENDIF

      END DO

      ! tdf
      !if (I .eq. 1 .and. J .eq. 63 .and. L .eq. 6) then
      !   print *,' GET_DUST_ALK: AREA'
      !   print *,' SULF_AREA ',SULF_AREA,'cm2/cm3'
      !   print *,'   BC_AREA ',BC_AREA,'cm2/cm3'
      !   print *,'   OC_AREA ',OC_AREA,'cm2/cm3'
      !   print *,'  SSA_AREA ',SSA_AREA,'cm2/cm3'
      !   print *,'  SSC_AREA ',SSC_AREA,'cm2/cm3'
      !   print *,'  AREA_d(1) ',AREA_d(1),'cm2/cm3'
      !   print *,'  AREA_d(2) ',AREA_d(2),'cm2/cm3'
      !   print *,'  AREA_d(3) ',AREA_d(3),'cm2/cm3'
      !   print *,'  AREA_d(4) ',AREA_d(4),'cm2/cm3'
      !endif

      ! tdf total aerosol surface area  [cm^2/cm^3]
      TOTAL_AREA = SULF_AREA + BC_AREA + OC_AREA
     &           + SSA_AREA  + SSC_AREA
     &           + AREA_d(1) + AREA_d(2) + AREA_d(3) + AREA_d(4)

      ! tdf total surface area weighted by gas-phase diffusion limitation [1/s]
      DF_TOTAL_AREA = SULF_AREA / SULF_FAC + BC_AREA  / BC_FAC
     &              + OC_AREA   / OC_FAC   + SSA_AREA / SSA_FAC
     &              + SSC_AREA  / SSC_FAC
     &              + DF_AREA_d(1) 
     &              + DF_AREA_d(2) 
     &              + DF_AREA_d(3) 
     &              + DF_AREA_d(4) 

      ! tdf
      !if (I .eq. 1 .and. J .eq. 63 .and. L .eq. 6) then
      !   print *,' GET_DUST_ALK: DF_AREA'
      !   print *,' DF_SULF_AREA  ',SULF_AREA/SULF_FAC,'1/s'
      !   print *,'   DF_BC_AREA  ',BC_AREA/BC_FAC,'1/s'
      !   print *,'   DF_OC_AREA  ',OC_AREA/OC_FAC,'1/s'
      !   print *,'  DF_SSA_AREA  ',SSA_AREA/SSA_FAC,'1/s'
      !   print *,'  DF_SSC_AREA  ',SSC_AREA/SSC_FAC,'1/s'
      !   print *,'  DF_AREA_d(1) ',DF_AREA_d(1),'1/s'
      !   print *,'  DF_AREA_d(2) ',DF_AREA_d(2),'1/s'
      !   print *,'  DF_AREA_d(3) ',DF_AREA_d(3),'1/s'
      !   print *,'  DF_AREA_d(4) ',DF_AREA_d(4),'1/s'
      !endif

      ! tdf Total Dust Alkalinity
      ALK = ALK_d(1) + ALK_d(2) + ALK_d(3) + ALK_d(4)  ! [v/v]

      ! set humidity index IRH as a percent
      IRH = State_Met%RH(I,J,L)
      IRH = MAX(  1, IRH )
      IRH = MIN( 99, IRH )

      ! hygroscopic growth factor for dust: Set to NO GROWTH for now
      IF ( IRH < 100 ) HGF = 1.0e+0_fp

      ! tdf Loop over size bins (NDSTBIN = 4)
      DO IBIN = 1, NDSTBIN

         !----------------------------------
         ! SO2 uptake onto particles
         !----------------------------------

         !tdf 2/11/2K9
         !tdf Following relative uptake rates of Preszler-Prince et al.(2007)
         IF (IRH >= 90 ) THEN     
            GAMS = GAMMA_SO2 * 2.e+0_fp
         ELSE IF (IRH >= 84 ) THEN  
            GAMS = GAMMA_SO2 * (0.5e+0_fp  + 1.5e+0_fp*(IRH-84.e+0_fp)/
     &             (90.e+0_fp-84.e+0_fp))
         ELSE IF (IRH >= 76 ) THEN 
            GAMS = GAMMA_SO2 * (0.16e+0_fp + 0.34e+0_fp*(IRH-76.e+0_fp)/
     &             (84.e+0_fp-76.e+0_fp))
         ELSE IF (IRH >= 33 ) THEN  
            GAMS = GAMMA_SO2 * (0.03e+0_fp + 0.13e+0_fp*(IRH-33.e+0_fp)/
     &             (76.e+0_fp-33.e+0_fp))
         ELSE IF (IRH >= 20 ) THEN 
            GAMS = GAMMA_SO2*0.03e+0_fp
         ELSE                        ! 0.0 below 20%
            GAMS = GAMMA_SO2*0.0e+0_fp
         ENDIF

         ! Check for sufficient alkalinity          tdf 3/28/2K8
         IF ( ALK > MINDAT ) THEN

            ! calculate gas-to-particle rate constant for uptake of
            ! SO2 onto dust aerosols [Jacob, 2000] analytical solution
            ! Corrected based on discussions with Becky     tdf 07/14/08
            KT1    = 0.0e+0_fp

            IF (IBIN .EQ. 1) THEN
               ! tdf Sum over the 1-4 sub-bins for bin 1      ! tdf 08/21/2K9
               DO ISBIN = 1, 4
                  RD     = RD_sd1 (ISBIN)        ! effective radius [cm]
                  AREA   = AREA_sd1 (ISBIN)      ! Dust Surface Area [cm^2/cm^3]

                  ! Prevent divide by zero if GAMS = 0 (tdf, mps, 11/14/13)
                  IF ( GAMS > 0.e+0_fp ) THEN
                     CONST1 = 4.e+0_fp/(V*GAMS)  ! Collision [s/cm]
                     CONST2 = RD/DG              ! Diffusion [s/cm]
                     CONST  = CONST1 + CONST2
                     KT1    = KT1 + AREA / CONST ! [cm^2/cm^3] * [cm/s] = [1/s]
                  ELSE
                     KT1    = KT1                ! [cm^2/cm^3] * [cm/s] = [1/s]
                  ENDIF

               END DO
            ELSE
               RD     = RD_d (IBIN)              ! effective radius [cm]
               AREA   = AREA_d (IBIN)            ! Dust Surface Area [cm^2/cm^3]

               ! Prevent divide by zero if GAMS = 0 (tdf, mps, 11/14/13)
               IF ( GAMS > 0.e+0_fp ) THEN
                  CONST1 = 4.e+0_fp/(V*GAMS)     ! Collision [s/cm]
                  CONST2 = RD/DG                 ! Diffusion [s/cm]
                  CONST  = CONST1 + CONST2
                  KT1    = AREA / CONST          ! [cm^2/cm^3] * [cm/s] = [1/s]
               ELSE
                  KT1    = 0.0e+0_fp             ! [cm^2/cm^3] * [cm/s] = [1/s]
               ENDIF

            ENDIF

            KTS(IBIN) = KT1

         ELSE

            ! If no alkalinity, set rate coefficients to zero
            !tdf
            KTS(IBIN)  = 0.0e+0_fp

         ENDIF

         !----------------------------------
         ! HNO3 uptake onto particles
         !----------------------------------

         !tdf 2/11/2K9
         !tdf Following uptake coefficients of Liu et al.(2007)
         IF (IRH >= 80 ) THEN         
            GAMN = GAMMA_HNO3 * 2.1e+0_fp
         ELSE IF (IRH >= 70 ) THEN   
            GAMN = GAMMA_HNO3 * (1.3e+0_fp  + 0.7e+0_fp   * 
     &             (IRH - 70.e+0_fp)/ 10.e+0_fp)
         ELSE IF (IRH >= 60 ) THEN 
            GAMN = GAMMA_HNO3 * (1.0e+0_fp  + 0.3e+0_fp   * 
     &             (IRH - 60.e+0_fp)/  10.e+0_fp)
         ELSE IF (IRH >= 50 ) THEN
            GAMN = GAMMA_HNO3 * (0.7e+0_fp  + 0.3e+0_fp   * 
     &             (IRH - 50.e+0_fp)/ 10.e+0_fp)
         ELSE IF (IRH >= 30 ) THEN  
            GAMN = GAMMA_HNO3 * (0.19e+0_fp + 0.255e+0_fp * 
     &             (IRH - 30.e+0_fp)/ 10.e+0_fp)
         ELSE IF (IRH >= 10 ) THEN  
            GAMN = GAMMA_HNO3 * (0.03e+0_fp + 0.08e+0_fp  * 
     &             (IRH - 10.e+0_fp)/ 10.e+0_fp)
         ELSE                        ! 0.0 below 10%
            GAMN = GAMMA_HNO3*0.0e+0_fp
         ENDIF

         ! Check for sufficient alkalinity      tdf 3/28/2K8
         IF ( ALK > MINDAT ) THEN

            ! calculate gas-to-particle rate constant for uptake of
            ! HNO3 onto dust aerosols [Jacob, 2000] analytical solution
            ! Corrected based on discussions with Becky     tdf 07/14/08
            KT1    = 0.0e+0_fp

            IF (IBIN .EQ. 1) THEN
               ! tdf Sum over the 1-4 sub-bins for bin 1      ! tdf 08/21/2K9
               DO ISBIN = 1, 4
                  RD = RD_sd1 (ISBIN)            ! effective radius [cm]
                  AREA = AREA_sd1 (ISBIN)        ! Dust Surface Area [cm^2/cm^3]

                  ! Prevent divide by zero if GAMN = 0 (tdf, mps, 11/14/13)
                  IF ( GAMN > 0.e+0_fp ) THEN
                     CONST1 = 4.e+0_fp/(V*GAMN)  ! Collision [s/cm]
                     CONST2 = RD/DG              ! Diffusion [s/cm]
                     CONST  = CONST1 + CONST2
                     KT1    = KT1 + AREA / CONST ! [cm^2/cm^3] * [cm/s] = [1/s]
                  ELSE
                     KT1    = KT1                ! [cm^2/cm^3] * [cm/s] = [1/s]
                  ENDIF

               END DO
            ELSE
               RD     = RD_d (IBIN)              ! effective radius [cm]
               AREA   = AREA_d (IBIN)            ! Dust Surface Area [cm^2/cm^3]

               ! Prevent divide by zero if GAMN = 0 (tdf, mps, 11/14/13)
               IF ( GAMN > 0.e+0_fp ) THEN
                  CONST1 = 4.e+0_fp/(V*GAMN)     ! Collision [s/cm]
                  CONST2 = RD/DG                 ! Diffusion [s/cm]
                  CONST  = CONST1 + CONST2
                  KT1    = AREA / CONST          ! [cm^2/cm^3] * [cm/s] = [1/s]
               ELSE
                  KT1    = 0.0e+0_fp             ! [cm^2/cm^3] * [cm/s] = [1/s]
               ENDIF

            ENDIF

            KTN(IBIN) = KT1

         ELSE

            ! If no alkalinity, set rate coefficients to zero
            !tdf
            KTN(IBIN)  = 0.0e+0_fp

         ENDIF

         !----------------------------------
         ! H2SO4 uptake onto particles
         !----------------------------------

         ! Uptake not limited by dust alkalinity      tdf 3/02/2K9

         !tdf As of 08/20/09, we use AREA and size weighted uptake
         !tdf where now KTH is a fractional uptake for each size bin
         !tdf with respect to total aerosol surface area.

         KT1    = DF_AREA_d(IBIN) / DF_TOTAL_AREA   ! Fraction

         KTH(IBIN) = KT1

      ! tdf End Loop over size bins
      END DO

      ! Free pointers
      Spc     => NULL()
      ERADIUS => NULL()
      TAREA   => NULL()

      END SUBROUTINE GET_DUST_ALK
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
      SUBROUTINE INIT_DUST( am_I_Root, Input_Opt,
     &                      State_Chm, State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : PI, AIRMW
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Diag_Mod,     ONLY : DgnState
#if   defined( TOMAS )
      USE TOMAS_MOD,          ONLY : IBINS, Xk
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  30 Mar 2004 - R. Yantosca - Initial version
!  (1 ) Now references LDEAD from "logical_mod.f" (bmy, 7/20/04)
!  (2 ) Modify to work with 30bin dust. Reference to IBINS from tomas_mod
!       for number of total bi}n = 30 bins.  (win, 7/17/09)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  26 Feb 2013 - M. Long     - Now use fields from Input_Opt
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  23 Sep 2015 - R. Yantosca - Now accept State_Chm as an argument so that
!                              we can use the species database
!  30 Sep 2015 - R. Yantosca - DD_A_Density is now renamed to Density
!  30 Sep 2015 - R. Yantosca - DD_A_Radius is now renamed to Radius
!  05 Jan 2016 - E. Lundgren - Use global parameter PI
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  16 Jun 2016 - J. Kaiser   - Now define species ID's with Ind_() function
!  20 Jun 2016 - R. Yantosca - Only define species ID's in the INIT phase
!  24 Jun 2016 - R. Yantosca - Add error checks for dust uptake species
!  05 Oct 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Non-SAVEd scalars
      INTEGER                :: AS, N
#if   defined( TOMAS )
      INTEGER                :: I
      INTEGER                :: BIN
#endif

      ! Strings
      CHARACTER(LEN=255)     :: ErrMsg
      CHARACTER(LEN=255)     :: ThisLoc

      ! Objects
      TYPE(Species), POINTER :: SpcInfo

      !=================================================================
      ! INIT_DUST begins here!
      !=================================================================

      ! Initialize
      RC      =  GC_SUCCESS
      ErrMsg  =  ''
      ThisLoc =  '-> at INIT_DUST (in module GeosCore/dust_mod.F)'
      SpcInfo => NULL()

      !=================================================================
      ! INITIALIZATION SECTION FOR STANDARD GEOS-CHEM SIMULATION
      !=================================================================
     
      !-----------------------------------------------------------------
      ! DST1 - DST4 species (most fullchem + aerosol simulations)
      !
      ! NOTE: We can consider removing DUSTDEN and DUSTREFF from 
      ! the Input_Opt object at a later point.  These can be directly
      ! obtained from the species database object now. (bmy, 6/20/16)
      !-----------------------------------------------------------------

      id_DST1                  =  Ind_('DST1'    )
      id_DST2                  =  Ind_('DST2'    )
      id_DST3                  =  Ind_('DST3'    )
      id_DST4                  =  Ind_('DST4'    )

      !-----------------------------------------------------------------
      ! DAL1 - DAL4 and SO4d1 - SO4d4 species (acid uptake sims only)
      !-----------------------------------------------------------------

      id_DAL1                  =  Ind_('DSTAL1'  )
      id_DAL2                  =  Ind_('DSTAL2'  )
      id_DAL3                  =  Ind_('DSTAL3'  )
      id_DAL4                  =  Ind_('DSTAL4'  )

      ! Error check the dust uptake species
      IF ( Input_Opt%LDSTUP ) THEN
         IF ( id_DAL1 < 0 .or. id_DAL2 < 0   .or.  
     &        id_DAL3 < 0 .or. id_DAL4 < 0 ) THEN
            ErrMsg = 'Dust uptake species are undefined!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

#if defined( TOMAS )
      !=================================================================
      ! INITIALIZATION SECTION FOR TOMAS MICROPHYSICS
      !=================================================================

      ! TOMAS species ID flags
      id_DUST1 = Ind_('DUST1')
      id_NK1   = Ind_('NK1'  ) 

      !----------------------------------
      ! Set up FRAC_S (only for Ginoux)
      !----------------------------------
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !% NOTE: LDEAD has not been set by HCOI_GC_INIT yet. This code needs %
      !% to be moved or modified accordingly (mps, 4/9/15)                 %    
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !% NOTE: This code replicates the obsolete Input_Opt%IDDEP field,    %
      !% which has since been removed.  TOMAS TEAM please take note.       %
      !% (bmy, 3/17/17)                                                    %
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF ( .not. Input_Opt%LDEAD ) THEN

         ! Initialize
         TOMAS_IDDEP = 0
         
         DO BIN = 1, IBINS
         DO N   = 1, State_Chm%nDryDep
            IF ( State_Chm%Map_DryDep(N) == ( id_NK1-1+BIN ) )THEN
               TOMAS_IDDEP(BIN) = N
               GOTO 100
            ENDIF            
         ENDDO
 100     CONTINUE
         ENDDO

      ENDIF
#endif 

#if defined( ESMF_ ) || defined( TOMAS )
      ! EXPERIMENTAL: For archiving the dust source for GCHP
      !
      ! Changed to use the ESMF_ flag and not EXTERNAL_GRID/EXTERNAL_FORCING, as
      ! WRF-GC which uses these flags does not require SRCE_FUNC (hplin, 1/22/19)
      ALLOCATE( SRCE_FUNC(IIPAR,JJPAR,3), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SRCE_FUNC' )
      SRCE_FUNC = 0.e+0_fp
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
      IF ( ALLOCATED( FRAC_S    ) ) DEALLOCATE( FRAC_S    )
      IF ( ALLOCATED( SRCE_FUNC ) ) DEALLOCATE( SRCE_FUNC )

      END SUBROUTINE CLEANUP_DUST
!EOC
      END MODULE DUST_MOD

