!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: convection_mod
!
! !DESCRIPTION: Module CONVECTION\_MOD contains routines which select the 
!  proper convection code for GEOS-3, GEOS-4, GEOS-5, MERRA, or GCAP met 
!  field data sets. 
!\\
!\\
! !INTERFACE: 
!
      MODULE CONVECTION_MOD
! 
! !USES:
!
      USE GC_TYPE_MOD
      USE GIGC_ERRCODE_MOD
      USE SMV_PHYSCONST_MOD

      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DO_CONVECTION
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: DO_GEOS4_CONVECT
      PRIVATE :: DO_GCAP_CONVECT
      PRIVATE :: NFCLDMX
      PRIVATE :: DO_MERRA_CONVECTION
!
! !REVISION HISTORY:
!  27 Jan 2004 - R. Yantosca - Initial version
!  (1 ) Contains new updates for GEOS-4/fvDAS convection.  Also now references
!        "error_mod.f".  Now make F in routine NFCLDMX a 4-D array to avoid
!        memory problems on the Altix. (bmy, 1/27/04)
!  (2 ) Bug fix: Now pass NTRACE elements of TCVV to FVDAS_CONVECT in routine 
!        DO_CONVECTION (bmy, 2/23/04)  
!  (3 ) Now references "logical_mod.f" and "tracer_mod.f" (bmy, 7/20/04)
!  (4 ) Now also references "ocean_mercury_mod.f" and "tracerid_mod.f" 
!        (sas, bmy, 1/19/05)
!  (5 ) Now added routines DO_GEOS4_CONVECT and DO_GCAP_CONVECT by breaking 
!        off code from DO_CONVECTION, in order to implement GCAP convection
!        in a much cleaner way. (swu, bmy, 5/25/05)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) Shut off scavenging in shallow convection for GCAP (swu, bmy, 11/1/05)
!  (8 ) Modified for tagged Hg simulation (cdh, bmy, 1/6/06)
!  (9 ) Bug fix: now only call ADD_Hg2_WD if LDYNOCEAN=T (phs, 2/8/07)
!  (10) Fix for GEOS-5 met fields in routine NFCLDMX (swu, 8/15/07)
!  (11) Resize DTCSUM array in NFCLDMX to save memory (bmy, 1/31/08)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as for GEOS-5
!  29 Sep 2010 - R. Yantosca - Added modifications for MERRA
!  05 Oct 2010 - R. Yantosca - Added ND14 and ND38 diagnostics to 
!                              DO_MERRA_CONVECTION routine
!  16 Aug 2011 - J. Fisher   - Minor bug fixes in DO_MERRA_CONVECTION
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  02 Mar 2012 - R. Yantosca - Now reference the new grid_mod.F90
!  22 Oct 2012 - R. Yantosca - Now reference Headers/gigc_errcode_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_convection
!
! !DESCRIPTION: Subroutine DO\_CONVECTION calls the appropriate convection 
!  driver program for different met field data sets.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_CONVECTION( am_I_Root, Input_Opt,
     &                          State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : CONVFLUP
      USE DIAG_MOD,           ONLY : AD38
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE TRACERID_MOD,       ONLY : IDTHg2
      USE TRACERID_MOD,       ONLY : IDTHgP
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE WETSCAV_MOD,        ONLY : COMPUTE_F
      USE WETSCAV_MOD,        ONLY : H2O2s
      USE WETSCAV_MOD,        ONLY : SO2s
#if   defined( APM )
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif
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
!  25 May 2005 - S. Wu       - Initial version
!  08 Feb 2007 - R. Yantosca - Now reference "CMN_SIZE".  Now references 
!                              CLDMAS, CMFMC, DTRAIN from "dao_mod.f" so that 
!                              we can pass either GEOS-5 or GEOS-3 meteorology
!                              to NFCLDMX. 
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as for GEOS-5
!  29 Sep 2010 - R. Yantosca - Now call DO_MERRA_CONVECTION for MERRA met
!  05 Oct 2010 - R. Yantosca - Now attach diagnostics to MERRA conv routine
!  06 Oct 2010 - R. Yantosca - Parallelized call to DO_MERRA_CONVECTION
!  15 Oct 2010 - H. Amos     - Now get BXHEIGHT, T from dao_mod.f
!  15 Oct 2010 - R. Yantosca - Now get LDYNOCEAN from logical_mod.f
!  15 Oct 2010 - R. Yantosca - Now get ITS_A_MERCURY_SIM from tracer_mod.f
!  15 Oct 2010 - R. Yantosca - Now get IDTHg2, IDTHgP from tracerid_mod.f
!  15 Oct 2010 - R. Yantosca - Now get H2O2s, SO2s from wetscav_mod.f
!  15 Oct 2010 - H. Amos     - Now pass BXHEIGHT, T, to DO_MERRA_CONVECTION
!  15 Oct 2010 - R. Yantosca - Now pass H2O2s, SO2s to DO_MERRA_CONVECTION
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  29 Aug 2011 - R. Yantosca - Bug fix: reposition #if defined(APM) statement
!  09 Feb 2012 - R. Yantosca - For GEOS-5.7, PFICU and PFLCU fields are
!                              defined on level edges.  Pass the top edge
!                              of each level to DO_MERRA_CONVECTION
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  21 Jun 2012 - R. Yantosca - Now use poiners to pass array slices to routines
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  04 Feb 2013 - S. Kim      - Bug fix: H2O2s, SO2s, STT are not in State_Met
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  31 May 2013 - R. Yantosca - Now pass objects to NFCLDMX
!  03 Jun 2013 - R. Yantosca - Bug fix: pass State_Chm to DO_MERRA_CONVECTION
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N_TOT_TRC               ! # of tracers to convect
 
      ! For fields from Input_Opt
      LOGICAL         :: LDYNOCEAN
      LOGICAL         :: ITS_A_MERCURY_SIM
      INTEGER         :: N_TRACERS
      REAL*8          :: TCVV(Input_Opt%N_TRACERS)
      REAL*8          :: TRACER_MW_KG(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! DO_CONVECT begins here!
      !=================================================================

#if   defined( GCAP ) 

      !=================================================================
      ! GCAP met fields
      !=================================================================

      ! Assume success
      RC                =  GIGC_SUCCESS

      ! Call GCAP driver routine
      CALL DO_GCAP_CONVECT
     &   ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#elif defined( GEOS_4 )

      !=================================================================
      ! GEOS-4 met fields
      !=================================================================

      ! Assume success
      RC                =  GIGC_SUCCESS

      ! Call GEOS-4 driver routine
      CALL DO_GEOS4_CONVECT
     &   ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#elif defined( GEOS_5 )

      !=================================================================
      ! GEOS-5 met fields
      !=================================================================

      ! Assume success
      RC                =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      ITS_A_MERCURY_SIM = Input_Opt%ITS_A_MERCURY_SIM
      N_TRACERS         = Input_Opt%N_TRACERS

      ! Pick the proper # of tracers
#if   defined( APM )  
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! If APM is used, then convect
                                         !  the G-C advected tracers +
#else                                    
      N_TOT_TRC = N_TRACERS              !  the APM aerosol tracers
#endif

      ! Call the S-J Lin convection routine for GEOS-5
      CALL NFCLDMX( am_I_Root, 
     &              Input_Opt, 
     &              State_Met%CMFMC(:,:,2:LLPAR+1),
     &              State_Met%DTRAIN,
     &              State_Met, 
     &              State_Chm, 
     &              RC                               )

#elif defined( MERRA ) || defined( GEOS_FP )

      !=================================================================
      ! MERRA met fields
      !=================================================================
      
      ! Objects
      TYPE(GC_IDENT)     :: IDENT
      TYPE(GC_DIMS )     :: DIMINFO
      TYPE(SPEC_2_TRAC)  :: COEF
      TYPE(GC_OPTIONS)   :: OPTIONS
      TYPE(ID_TRAC)      :: IDT
                   
      ! Scalars          
      INTEGER            :: I, J, L, N, NN, TS_DYN
      REAL*8             :: AREA_M2, DT
                         
      ! Arrays           
      REAL*8             :: DIAG14 ( LLPAR,   Input_Opt%N_TRACERS )
      REAL*8             :: DIAG38 ( LLPAR,   Input_Opt%N_TRACERS )
      REAL*8             :: F      ( LLPAR,   Input_Opt%N_TRACERS )
      REAL*8,  TARGET    :: FSOL   ( IIPAR, JJPAR, LLPAR,   
     &                                        Input_Opt%N_TRACERS )
      INTEGER            :: ISOL   (          Input_Opt%N_TRACERS )
      REAL*8             :: PEDGE  ( LLPAR+1                      )

      ! Pointers
      REAL*8,  POINTER   :: p_AD      (:    )   
      REAL*8,  POINTER   :: p_BXHEIGHT(:    )
      REAL*8,  POINTER   :: p_CMFMC   (:    )
      REAL*8,  POINTER   :: p_DQRCU   (:    )
      REAL*8,  POINTER   :: p_DTRAIN  (:    )
      REAL*8,  POINTER   :: p_H2O2s   (:    )
      REAL*8,  POINTER   :: p_FSOL    (:,:,:)
      REAL*8,  POINTER   :: p_PFICU   (:    )   
      REAL*8,  POINTER   :: p_PFLCU   (:    )
      REAL*8,  POINTER   :: p_REEVAPCN(:    )
      REAL*8,  POINTER   :: p_SO2s    (:    )
      REAL*8,  POINTER   :: p_STT     (:,:  )
      REAL*8,  POINTER   :: p_T       (:    )
      
      ! Parameters
      LOGICAL, PARAMETER :: FIRST = .TRUE.

      !-----------------------------------------------------------------
      ! Initialization
      !-----------------------------------------------------------------

      ! Assume success
      RC                =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LDYNOCEAN         = Input_Opt%LDYNOCEAN
      ITS_A_MERCURY_SIM = Input_Opt%ITS_A_MERCURY_SIM
      N_TRACERS         = Input_Opt%N_TRACERS
      TCVV              = Input_Opt%TCVV(1:N_TRACERS)
      TRACER_MW_KG      = Input_Opt%TRACER_MW_KG(1:N_TRACERS)

      ! Pick the proper # of tracers
#if   defined( APM )  
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! If APM is used, then convect
#else                                    !  the G-C advected tracers +
      N_TOT_TRC = N_TRACERS              !  the APM aerosol tracers
#endif

      ! Define variables and object fields
      TS_DYN                  = GET_TS_DYN()         ! Dynamic timestep [min]
      DT                      = DBLE( TS_DYN )       ! Dynamic timestep [min]
      DIMINFO%L_COLUMN        = LLPAR                ! # of vertical levels
      DIMINFO%N_TRACERS       = N_TRACERS            ! # of tracers
      IDT%Hg2                 = IDTHg2               ! Tracer ID flag for Hg2
      IDT%HgP                 = IDTHgP               ! Tracer ID flag for HgP
      OPTIONS%USE_Hg          = ITS_A_MERCURY_SIM    ! Mercury simulation?
      OPTIONS%USE_Hg_DYNOCEAN = LDYNOCEAN            !  with dynamic ocean?
      OPTIONS%USE_DIAG14      = ( ND14 > 0 )         ! Use ND14 diagnostic?
      OPTIONS%USE_DIAG38      = ( ND38 > 0 )         ! Use ND38 diagnostic?

      ! Allocate the pointer array
      ALLOCATE( COEF%MOLWT_KG( N_TRACERS ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      COEF%MOLWT_KG = 0d0

      ! Loop over advected tracers
      DO N = 1, N_TRACERS

         ! Molecular weight [kg]
         COEF%MOLWT_KG(N) = TRACER_MW_KG(N)

         ! Set a pointer to each 3D slice of FSOL
         p_FSOL => FSOL(:,:,:,N)

         ! Fraction of soluble tracer
         CALL COMPUTE_F( N, p_FSOL, ISOL(N), State_Met, State_Chm ) 

         ! Nullify the pointer
         NULLIFY( p_FSOL )

      ENDDO
      
      !-----------------------------------------------------------------
      ! Do convection column by column
      !-----------------------------------------------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J,          AREA_M2,  I,       IDENT,      L       )
!$OMP+PRIVATE( PEDGE,      N,        F,       DIAG14,     DIAG38  )
!$OMP+PRIVATE( RC,         NN,       p_AD,    p_BXHEIGHT, p_CMFMC )
!$OMP+PRIVATE( p_DQRCU,    p_DTRAIN, p_H2O2s, p_PFICU,    p_PFLCU )
!$OMP+PRIVATE( P_REEVAPCN, p_SO2s,   p_STT,   p_T                 )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]
         AREA_M2           =  GET_AREA_M2( I, J, 1 ) 

         ! Initialize the IDENT object for error I/O
         IDENT%STDOUT_FILE =  ''                     ! Write to
         IDENT%STDOUT_LUN  =  6                      !  stdout
         IDENT%PET         =  0                      ! CPU #
         IDENT%I_AM(1)     =  'MAIN'                 ! Main routine
         IDENT%I_AM(2)     =  'DO_CONVECTION'        ! This routine
         IDENT%LEV         =  2                      ! This level
         IDENT%ERRMSG      =  ''                     ! Error msg
        !IDENT%VERBOSE     =  ( I==23 .and. J==34 )  ! Debug box
         IDENT%VERBOSE     =  ( I==67 .and. J==31 )  ! Debug box

         ! Pressure edges
         DO L = 1, LLPAR+1
            PEDGE(L)       =  GET_PEDGE( I, J, L )
         ENDDO

         ! For some reason, the parallelization chokes if we pass an array 
         ! slice of FSOL to DO_MERRA_CONVECTION.  It works if we manually 
         ! copy the data into a local array and then pass that to the routine.
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
            F(L,N)         =  FSOL(I,J,L,N)
         ENDDO
         ENDDO

#if   defined( GEOS_FP )

         ! For GEOS-5.7, PFICU and PFLCU are defined on level edges
         p_PFICU           => State_Met%PFICU   ( I, J, 2:LLPAR+1   )
         p_PFLCU           => State_Met%PFLCU   ( I, J, 2:LLPAR+1   )

#elif defined( MERRA   )

         ! For MERRA, PFICU and PFLCU are defined on level centers
         p_PFICU           => State_Met%PFICU   ( I, J, 1:LLPAR     )
         p_PFLCU           => State_Met%PFLCU   ( I, J, 1:LLPAR     )

#endif

         ! Pointers to slices of larger arrays
         p_AD              => State_Met%AD      ( I, J,  :          )
         p_BXHEIGHT        => State_Met%BXHEIGHT( I, J,  :          )
         p_CMFMC           => State_Met%CMFMC   ( I, J, 2:LLPAR+1   )
         p_DQRCU           => State_Met%DQRCU   ( I, J,  :          )
         p_DTRAIN          => State_Met%DTRAIN  ( I, J,  :          )
         p_H2O2s           =>           H2O2s   ( I, J,  :          )
         p_REEVAPCN        => State_Met%REEVAPCN( I, J,  :          )
         p_SO2s            =>           SO2s    ( I, J,  :          )
         p_STT             => State_Chm%Tracers ( I, J,  :,       : )
         p_T               => State_Met%T       ( I, J,  :          )
         
         !--------------------------
         ! Do the cloud convection
         !--------------------------
         CALL DO_MERRA_CONVECTION( IDENT     = IDENT,
     &                             DIMINFO   = DIMINFO,
     &                             COEF      = COEF,
     &                             IDT       = IDT,
     &                             OPTIONS   = OPTIONS,
     &                             AREA_M2   = AREA_M2,
     &                             PEDGE     = PEDGE,
     &                             TS_DYN    = DT, 
     &                             F         = F,
     &                             AD        = p_AD,
     &                             BXHEIGHT  = p_BXHEIGHT,
     &                             CMFMC     = p_CMFMC,
     &                             DQRCU     = p_DQRCU,
     &                             DTRAIN    = p_DTRAIN,
     &                             H2O2s     = p_H2O2s,
     &                             PFICU     = p_PFICU,
     &                             PFLCU     = p_PFLCU,
     &                             REEVAPCN  = p_REEVAPCN,            
     &                             SO2s      = p_SO2s,
     &                             T         = p_T,
     &                             Q         = p_STT,
     &                             DIAG14    = DIAG14,
     &                             DIAG38    = DIAG38,
     &                             I         = I,
     &                             J         = J,
     &                             State_Met = State_Met,
     &                             State_Chm = State_Chm,
     &                             RC        = RC         )

         ! Stop if error is encountered
         IF ( RC /= GIGC_SUCCESS ) THEN
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Free pointer memory
         NULLIFY( p_AD,       p_BXHEIGHT, p_CMFMC, p_DQRCU )
         NULLIFY( p_DTRAIN,   p_H2O2s,    p_PFICU, p_PFLCU )
         NULLIFY( p_REEVAPCN, p_SO2s,     p_STT,   p_T     )

         !--------------------------
         ! ND14 diagnostic [kg/s]
         !--------------------------
         IF ( ND14 > 0 ) THEN
            DO N = 1, N_TRACERS
            DO L = 1, LD14
               CONVFLUP(I,J,L,N) = CONVFLUP(I,J,L,N) + DIAG14(L,N)
            ENDDO
            ENDDO
         ENDIF

         !--------------------------
         ! ND38 diagnostic [kg/s]
         !--------------------------
         IF ( ND38 > 0 ) THEN
            DO N = 1, N_TRACERS

               ! Wetdep species for each tracer N
               NN = ISOL(N)

               ! Only archive wet-deposited species
               IF ( NN > 0 ) THEN
                  DO L = 1, LD38
                     AD38(I,J,L,NN) = AD38(I,J,L,NN) + DIAG38(L,N)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free the memory
      DEALLOCATE( COEF%MOLWT_KG )

#endif

      END SUBROUTINE DO_CONVECTION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_geos4_convect
!
! !DESCRIPTION: Subroutine DO\_GEOS4\_CONVECT is a wrapper for the 
!  GEOS-4/fvDAS convection code.  This was broken off from the old 
!  DO\_CONVECTION routine above.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_GEOS4_CONVECT( am_I_Root, Input_Opt,
     &                             State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD37
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE FVDAS_CONVECT_MOD,  ONLY : INIT_FVDAS_CONVECT, FVDAS_CONVECT
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE TIME_MOD,           ONLY : GET_TS_CONV
      USE WETSCAV_MOD,        ONLY : COMPUTE_F
#if   defined( APM ) 
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif

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
!  25 May 2005 - S. Wu       - Initial version
!  (1 ) Now use array masks to flip arrays vertically in call to FVDAS_CONVECT
!        (bmy, 5/25/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Add a check to set negative values in STT to TINY (ccc, 4/15/09)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE   :: FIRST = .TRUE.
      INTEGER         :: I, ISOL, J, L, L2, N, NSTEP
      !INTEGER         :: INDEXSOL(N_TRACERS)
      INTEGER         :: CONVDT
      !REAL*8          :: F(IIPAR,JJPAR,LLPAR,N_TRACERS)
      REAL*8          :: RPDEL(IIPAR,JJPAR,LLPAR)
      REAL*8          :: DP(IIPAR,JJPAR,LLPAR)
      REAL*8          :: P1, P2, TDT
      INTEGER         :: N_TOT_TRC
#if   defined( APM )  
      INTEGER         :: INDEXSOL(Input_Opt%N_TRACERS+N_APMTRA)
      REAL*8          :: F       (IIPAR,JJPAR,LLPAR,
     &                            Input_Opt%N_TRACERS+N_APMTRA)
#else                 
      INTEGER         :: INDEXSOL(Input_Opt%N_TRACERS)
      REAL*8          :: F       (IIPAR,JJPAR,LLPAR,
     &                            Input_Opt%N_TRACERS)
#endif

      ! For fields from Input_Opt
      LOGICAL         :: LPRT
      INTEGER         :: N_TRACERS
      REAL*8          :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL*8, POINTER :: STT     (:,:,:,:)   ! Ptr to Tracer concentrations
      REAL*8, POINTER :: p_HKETA (:,:,:  )   ! Ptr to HKETA  met field
      REAL*8, POINTER :: p_HKBETA(:,:,:  )   ! Ptr to HKBETA met field
      REAL*8, POINTER :: p_ZMMU  (:,:,:  )   ! Ptr to ZMMU   met field
      REAL*8, POINTER :: p_ZMMD  (:,:,:  )   ! Ptr to ZMMD   met field
      REAL*8, POINTER :: p_ZMEU  (:,:,:  )   ! Ptr to ZMEU   met field

      !=================================================================
      ! DO_GEOS4_CONVECT begins here!
      !=================================================================
      
      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LPRT      = Input_Opt%LPRT
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Pick the proper # of tracers for convection
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA  
#else
      N_TOT_TRC = N_TRACERS
#endif

      ! Convection timestep [s]
      CONVDT = GET_TS_CONV() * 60d0 
       
      ! NSTEP is the # of internal convection timesteps.  According to
      ! notes in the old convection code, 300s works well. (swu, 12/12/03)
      NSTEP  = CONVDT / 300    
      NSTEP  = MAX( NSTEP, 1 ) 

      ! TIMESTEP*2; will be divided by 2 before passing to CONVTRAN 
      TDT    = DBLE( CONVDT ) * 2.0D0 / DBLE( NSTEP )

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_FVDAS_CONVECT
         FIRST = .FALSE.
      ENDIF
         
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV: a INIT_FV' )

      !=================================================================
      ! Before calling convection, compute the fraction of insoluble
      ! tracer (Finsoluble) lost in updrafts.  Finsoluble = 1-Fsoluble.
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, ISOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, N_TOT_TRC

         ! Get fraction of tracer scavenged and the soluble tracer 
         ! index (ISOL). For non-soluble tracers, F=0 and ISOL=0.
         CALL COMPUTE_F( N, F(:,:,:,N), ISOL, State_Met, State_Chm ) 
         
         ! Store ISOL in an array for later use
         INDEXSOL(N) = ISOL

         ! Loop over grid boxes
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! ND37 diagnostic: store fraction of tracer 
            ! lost in moist convective updrafts ("MC-FRC-$")
            IF ( ND37 > 0 .and. ISOL > 0 .and. L <= LD37 ) THEN
               AD37(I,J,L,ISOL) = AD37(I,J,L,ISOL) + F(I,J,L,N) 
            ENDIF

            ! GEOS-4 convection routines need the insoluble fraction
            F(I,J,L,N) = 1d0 - F(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
       
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV: a COMPUTE_F' )

      !=================================================================
      ! Compute pressure thickness arrays DP and RPDEL
      ! These arrays are indexed from atm top --> surface
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, L2, P1, P2 )
      DO L = 1, LLPAR

         ! L2 runs from the atm top down to the surface
         L2 = LLPAR - L + 1

         ! Loop over surface grid boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR
               
            ! Pressure at bottom and top edges of grid box [hPa]
            P1 = GET_PEDGE(I,J,L)
            P2 = GET_PEDGE(I,J,L+1)

            ! DP = Pressure difference between top & bottom edges [Pa]
            DP(I,J,L2) = ( P1 - P2 ) * 100.0d0

            ! RPDEL = reciprocal of DP [1/hPa]
            RPDEL(I,J,L2) = 100.0d0 / DP(I,J,L2) 
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV: a DP, RPDEL' )
   
      !=================================================================
      ! Flip arrays in the vertical and call FVDAS_CONVECT
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT      => State_Chm%Tracers

      ! Point to met fields and flip in vertical
      p_HKETA  => State_Met%HKETA (:,:,LLPAR:1:-1  )
      p_HKBETA => State_Met%HKBETA(:,:,LLPAR:1:-1  )
      p_ZMMU   => State_Met%ZMMU  (:,:,LLPAR:1:-1  )
      p_ZMMD   => State_Met%ZMMD  (:,:,LLPAR:1:-1  )
      p_ZMEU   => State_Met%ZMEU  (:,:,LLPAR:1:-1  )

      ! Call the fvDAS convection routines (originally from NCAR!)
      CALL FVDAS_CONVECT( TDT,   
     &                    N_TOT_TRC, 
     &                    STT   (:,:,LLPAR:1:-1,:),    
     &                    RPDEL,         
     &                    p_HKETA,
     &                    p_HKBETA,
     &                    p_ZMMU,
     &                    p_ZMMD,
     &                    p_ZMEU,
     &                    DP,     
     &                    NSTEP,    
     &                    F     (:,:,LLPAR:1:-1,:),         
     &                    TCVV,   
     &                    INDEXSOL,
     &                    State_Met )

      

      ! Add a check to set negative values in STT to TINY
      ! (ccc, 4/15/09)
!$OMP PARALLEL DO        
!$OMP+DEFAULT( SHARED   )
!$OMP+PRIVATE( N, L, J, I )
      DO N = 1,N_TOT_TRC
      DO L = 1,LLPAR
      DO J = 1,JJPAR
      DO I = 1,IIPAR
         STT(I,J,L,N) = MAX(STT(I,J,L,N),TINY(1d0))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      NULLIFY( STT                                       )
      NULLIFY( p_HKETA, p_HKBETA, p_ZMMU, p_ZMMD, p_ZMEU )

      !### Debug! 
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV: a FVDAS_CONVECT' )

      END SUBROUTINE DO_GEOS4_CONVECT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_gcap_convect
!
! !DESCRIPTION: Subroutine DO\_GCAP\_CONVECT is a wrapper for the GCAP 
!  convection code.  This was broken off from the old DO\_CONVECTION routine 
!  above.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_GCAP_CONVECT( am_I_Root, Input_Opt,
     &                            State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD37
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GCAP_CONVECT_MOD,   ONLY : GCAP_CONVECT
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PEDGE, GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_CONV
      USE WETSCAV_MOD,        ONLY : COMPUTE_F
#if   defined( APM )
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif
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
!  25 May 2005 - S. Wu       - Initial version
!  (1 ) Now use array masks to flip arrays vertically in call to GCAP_CONVECT
!        (bmy, 5/25/05)
!  (2 ) Shut off scavenging in shallow convection for GCAP below 700 hPa
!        (swu, bmy, 11/1/05)
!  (3 ) Add a check to set negative values in STT to TINY (ccc, 4/15/09)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  04 Nov 2013 - M. Sulprizio- Now use pointer variables to flip met fields in
!                              the vertical
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I, ISOL, J, L, L2, N, NSTEP  
      !INTEGER       :: INDEXSOL(N_TRACERS) 
      INTEGER       :: CONVDT    
      !REAL*8        :: F(IIPAR,JJPAR,LLPAR,N_TRACERS)
      REAL*8        :: DP(IIPAR,JJPAR,LLPAR)
      REAL*8        :: P1, P2, TDT   
      REAL*8        :: GAINMASS
      INTEGER       :: N_TOT_TRC      
#if   defined( APM )
      INTEGER       :: INDEXSOL(Input_Opt%N_TRACERS+N_APMTRA) 
      REAL*8        :: F       (IIPAR,JJPAR,LLPAR,
     &                          Input_Opt%N_TRACERS+N_APMTRA)
#else
      INTEGER       :: INDEXSOL(Input_Opt%N_TRACERS) 
      REAL*8        :: F       (IIPAR,JJPAR,LLPAR,
     &                          Input_Opt%N_TRACERS)
#endif

      ! For fields from Input_Opt
      LOGICAL       :: LPRT
      INTEGER       :: N_TRACERS
      REAL*8        :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL*8, POINTER :: STT        (:,:,:,:)   ! Ptr to Tracer concentrations
      REAL*8, POINTER :: p_UPDE     (:,:,:  )   ! Ptr to UPDE     met field
      REAL*8, POINTER :: p_DNDE     (:,:,:  )   ! Ptr to DNDE     met field
      REAL*8, POINTER :: p_ENTRAIN  (:,:,:  )   ! Ptr to ENTRAIN  met field
      REAL*8, POINTER :: p_DETRAINE (:,:,:  )   ! Ptr to DETRAINE met field
      REAL*8, POINTER :: p_UPDN     (:,:,:  )   ! Ptr to UPDN     met field
      REAL*8, POINTER :: p_DNDN     (:,:,:  )   ! Ptr to DNDN     met field
      REAL*8, POINTER :: p_DETRAINN (:,:,:  )   ! Ptr to DETRAINN met field

      !=================================================================
      ! DO_GCAP_CONVECT begins here!
      !=================================================================
      
      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LPRT      = Input_Opt%LPRT
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Pick the proper # of tracers for convection
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! G-C tracers + APM tracers
#else
      N_TOT_TRC = N_TRACERS              ! G-C tracers only
#endif

      ! Test??
      GAINMASS = 0d0

      ! Convection timestep [s]
      CONVDT = GET_TS_CONV() * 60d0 
       
      ! NSTEP is the # of internal convection timesteps.  According to
      ! notes in the old convection code, 300s works well. (swu, 12/12/03)
      NSTEP  = CONVDT / 300    
      NSTEP  = MAX( NSTEP, 1 ) 

      ! TIMESTEP*2; will be divided by 2 before passing to CONVTRAN 
      TDT    = DBLE( CONVDT ) * 2.0D0 / DBLE( NSTEP )
         
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GCAP_CONV: a INIT_FV' )

      !=================================================================
      ! Before calling convection, compute the fraction of insoluble
      ! tracer (Finsoluble) lost in updrafts.  Finsoluble = 1-Fsoluble.
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, ISOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, N_TOT_TRC

         ! Get fraction of tracer scavenged and the soluble tracer 
         ! index (ISOL). For non-soluble tracers, F=0 and ISOL=0.
         CALL COMPUTE_F( N, F(:,:,:,N), ISOL, State_Met, State_Chm ) 
         
         ! Store ISOL in an array for later use
         INDEXSOL(N) = ISOL

         ! Loop over grid boxes
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Shut off scavenging in shallow convection for GCAP
            ! (swu, bmy, 11/1/05)
            IF ( GET_PCENTER( I, J, L ) > 700d0 ) F(I,J,L,N) = 0d0

            ! ND37 diagnostic: store fraction of tracer 
            ! lost in moist convective updrafts ("MC-FRC-$")
            IF ( ND37 > 0 .and. ISOL > 0 .and. L <= LD37 ) THEN
               AD37(I,J,L,ISOL) = AD37(I,J,L,ISOL) + F(I,J,L,N) 
            ENDIF

            ! GEOS-4 convection routines need the insoluble fraction
            F(I,J,L,N) = 1d0 - F(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
       
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GCAP_CONV: a COMPUTE_F' )

      !=================================================================
      ! Compute pressure thickness arrays DP and RPDEL
      ! These arrays are indexed from atm top --> surface
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, L2, P1, P2 )
      DO L = 1, LLPAR

         ! L2 runs from the atm top down to the surface
         L2 = LLPAR - L + 1

         ! Loop over surface grid boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR
               
            ! Pressure at bottom and top edges of grid box [hPa]
            P1 = GET_PEDGE(I,J,L)
            P2 = GET_PEDGE(I,J,L+1)

            ! DP = Pressure difference between top & bottom edges [Pa]
            DP(I,J,L2) = ( P1 - P2 ) * 100.0d0
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GCAP_CONV: a DP, RPDEL' )
   
      !=================================================================
      ! Flip arrays in the vertical and call FVDAS_CONVECT
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Point to met fields and flip in vertical
      p_UPDE     => State_Met%UPDE     (:,:,LLPAR:1:-1  )
      p_DNDE     => State_Met%DNDE     (:,:,LLPAR:1:-1  )
      p_ENTRAIN  => State_Met%ENTRAIN  (:,:,LLPAR:1:-1  )
      p_DETRAINE => State_Met%DETRAINE (:,:,LLPAR:1:-1  )
      p_UPDN     => State_Met%UPDN     (:,:,LLPAR:1:-1  )
      p_DNDN     => State_Met%DNDN     (:,:,LLPAR:1:-1  )
      p_DETRAINN => State_Met%DETRAINN (:,:,LLPAR:1:-1  )

      ! Call the GCAP convection routines 
      CALL GCAP_CONVECT( TDT, 
     &                   STT      (:,:,LLPAR:1:-1,:), 
     &                   N_TOT_TRC,
     &                   DP, ! I think this is the correct way (bmy, 5/25/05)
     &                   NSTEP, 
     &                   F        (:,:,LLPAR:1:-1,:),     
     &                   TCVV,   
     &                   INDEXSOL, 
     &                   p_UPDE,
     &                   p_DNDE,
     &                   p_ENTRAIN,
     &                   p_DETRAINE,
     &                   p_UPDN,
     &                   p_DNDN,
     &                   p_DETRAINN )

      ! Add a check to set negative values in STT to TINY
      ! (ccc, 4/15/09)
!$OMP PARALLEL DO        
!$OMP+DEFAULT( SHARED   )
!$OMP+PRIVATE( N, L, J, I )
      DO N = 1,N_TOT_TRC
      DO L = 1,LLPAR
      DO J = 1,JJPAR
      DO I = 1,IIPAR
         STT(I,J,L,N) = MAX(STT(I,J,L,N),TINY(1d0))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT                                   )
      NULLIFY( p_UPDE, p_DNDE, p_ENTRAIN, p_DETRAINE )
      NULLIFY( p_UPDN, p_DNDN, p_DETRAINN            )

      !### Debug! 
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GCAP_CONV: a GCAP_CONVECT' )

      END SUBROUTINE DO_GCAP_CONVECT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nfcldmx
!
! !DESCRIPTION: Subroutine NFCLDMX is S-J Lin's cumulus transport module for 
!  3D GSFC-CTM, modified for the GEOS-Chem model.  The "NF" stands for "no 
!  flipping", and denotes that you don't have to flip the tracer array Q in 
!  the main program before passing it to NFCLDMX.
!\\
!\\
!  NOTE: NFCLDMX can be used with GEOS-1, GEOS-STRAT, and GEOS-3 met fields.
!  For GEOS-4/fVDAS, you must use the routines in "fvdas\_convect\_mod.f"
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NFCLDMX( am_I_Root, Input_Opt, CLDMAS, 
     &                    DTRN,      State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_WD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_HgP_WD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_SNOWPACK
      USE DIAG_MOD,           ONLY : AD37
      USE DIAG_MOD,           ONLY : AD38
      USE DIAG_MOD,           ONLY : CONVFLUP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE PRESSURE_MOD,       ONLY : GET_BP
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE TIME_MOD,           ONLY : GET_TS_CONV
      USE TRACER_MOD,         ONLY : ITS_A_MERCURY_SIM
      USE TRACERID_MOD,       ONLY : IS_Hg2
      USE TRACERID_MOD,       ONLY : IS_HgP 
      USE WETSCAV_MOD,        ONLY : COMPUTE_F
      USE ERROR_MOD,          ONLY : ALLOC_ERR

      IMPLICIT NONE

!
! !INPUT PARAMETERS: 
!
      ! Are we on the root CPU
      LOGICAL,        INTENT(IN)    :: am_I_Root

      ! CLDMAS : Cloud mass flux (at upper edges of each level) [kg/m2/s]
      REAL*8,         INTENT(IN)    :: CLDMAS(IIPAR,JJPAR,LLPAR)

      ! Detrainment mass flux [kg/m2/s]
      REAL*8,         INTENT(IN)    :: DTRN(IIPAR,JJPAR,LLPAR)

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

! !REMARKS:
!  (1) The "NF" stands for "no flipping", and denotes that you don't have to 
!  flip the tracer array Q in the main program before passing it to NFCLDMX.
!  (bmy, 2/12/97, 1/31/08)
!
!  (2) This version has been customized to work with GEOS-5 met fields.
!
!  Reference:
!  ============================================================================
!  Lin, SJ.  "Description of the parameterization of cumulus transport
!     in the 3D Goddard Chemistry Transport Model, NASA/GSFC, 1996.
!
!  Vertical indexing:
!  ============================================================================
!  The indexing of the vertical sigma levels has been changed from 
!  SJ-Lin's original code:
!
!                 Old Method          New Method
!                  (SJ Lin)    
!
!               ------------------------------------- Top of Atm.
!                  k = 1               k = NLAY
!               ===================================== Max Extent 
!                  k = 2               k = NLAY-1      of Clouds
!               -------------------------------------
!              
!                   ...                 ...             
!
!               -------------------------------------
!                  k = NLAY-3          k = 4
!               -------------------------------------
!                  k = NLAY-2          k = 3
!               ------------------------------------- Cloud base
!                  k = NLAY-1          k = 2
!               -    -    -    -    -    -    -    - 
!                  k = NLAY            k = 1
!               ===================================== Ground
!
!      which means that:
!
!                 Old Method                      New Method
!                  (SJ Lin)
!
!            k-1      ^                      k+1      ^
!            ---------|---------             ---------|---------
!                     |                               |
!                  CMFMC(k)                       CMFMC(k)
!             
!                                 becomes
!            k     DTRAIN(k),                k     DTRAIN(k),       
!                 QC(k), Q(k)                     QC(k), Q(k)   
!          
!                     ^                               ^
!            ---------|---------             ---------|---------
!                     |                               |   
!            k+1   CMFMC(k+1)                k-1   CMFMC(k-1)
!
!
!      i.e., the lowest level    used to be  NLAY  but is now  1
!            the level below k   used to be  k+1   but is now  k-1.
!            the level above k   used to be  k-1   but is now  k+1
!            the top of the atm. used to be  1     but is now  NLAY.
!
!  The old method required that the vertical dimensions of the CMFMC, DTRAIN, 
!  and Q arrays had to be flipped before and after calling CLDMX.  Also, 
!  diagnostic arrays generated within CLDMX also had to be flipped.  The new 
!  indexing eliminates this requirement (and also saves on array operations).  
!
!  Major Modifications:
!  ============================================================================
!  Original Author:   Shian-Jiann Lin, Code 910.3, NASA/GSFC
!  Original Release:  12 February 1997
!                     Version 3, Detrainment and Entrainment are considered.
!                     The algorithm reduces to that of version 2 if Dtrn = 0.
!                                                                             .
!  Modified By:       Bob Yantosca, for Harvard Atmospheric Sciences
!  Modified Release:  27 January 1998
!                     Version 3.11, contains features of V.3 but also 
!                     scavenges soluble tracer in wet convective updrafts.
!                                                                             .
!                     28 April 1998
!                     Version 3.12, now includes mass flux diagnostic
!                                                                             .
!                     11 November 1999
!                     Added mass-flux diagnostics
!                                                                             .
!                     04 January 2000
!                     Updated scavenging constant AS2
!                                                                             .
!                     14 March 2000
!                     Added new wet scavenging code and diagnostics
!                     based on the GMI algorithm
!                                                                             .
!                     02 May 2000
!                     Added parallel loop over tracers! 
! 
! !REVISION HISTORY: 
!  12 Feb 1997 - M. Prather  - Initial version
!  (1 ) NFCLDMX is written in Fixed-Form Fortran 90.
!  (2 ) Added TCVV to the argument list.  Also cleaned up argument
!        and local variable declarations. (bey, bmy, 11/10/99)
!  (3 ) AD38 and CONVFLUP are now declared allocatable in "diag_mod.f". 
!        (bmy, 11/29/99)
!  (4 ) Bug fix for tagged CO tracer run (bey, bmy, 1/4/00)
!  (5 ) Add new routines for computing scavenging coefficients,
!        as well as adding the AD37 diagnostic array. (bmy, 3/14/00)
!  (6 ) Updated comments (bmy, 10/2/01)
!  (7 ) Now print a header to stdout on the first call, to confirm that 
!        NFCLDMX has been called (bmy, 4/15/02)
!  (8 ) Remove PZ from the arg list -- it isn't used! (bmy, 8/22/02)
!  (9 ) Fixed ND38 diagnostic so that it now reports correctly (must divide
!        by DNS).  Updatec comments, cosmetic changes. (bmy, 1/27/03)
!  (10) Bug fix: remove duplicate K from PRIVATE declaration (bmy, 3/23/03)
!  (11) Now removed all arguments except NC, TCVV, Q from the arg list -- the
!        other arguments can be supplied via F90 modules.  Now references
!        "dao_mod.f", "grid_mod.f", "pressure_mod.f", and "time_mod.f".
!        (bmy, 3/27/03)
!  (12) Bundled into "convection_mod.f" (bmy, 6/26/03)
!  (13) Make sure K does not go out of bounds in ND38 diagnostic.  Now make 
!        F a 4-D array in order to avoid memory problems on the Altix.  
!        (bmy, 1/27/04)
!  (14) Now references both "ocean_mercury_mod.f" and "tracerid_mod.f".
!        Now call ADD_Hg2_WD from "ocean_mercury_mod.f" to pass the amt of Hg2
!        lost by wet scavenging (sas, bmy, 1/19/05)
!  (15) Now references IS_Hg2 from "tracerid_mod.f".  Now pass tracer # IC
!        to ADD_Hg2_WD. (cdh, bmy, 1/6/06)
!  (16) Bug fix: now only call ADD_Hg2_WD if LDYNOCEAN=T (phs, 2/8/07)
!  (17) Now make CLDMAS, DTRN as arguments, so that we can pass either
!        GEOS-3 or GEOS-3 met data.  Redimension DTCSUM with NC instead of 
!        NNPAR.  In many cases, NC is less than NNPAR and this will help to 
!        save memory especially when running at 2x25 or greater resolution 
!        (bmy, 1/31/08)
!  (18) Add a check to set negative values in Q to TINY (ccc, 4/15/09)
!  (19) Updates for mercury simulation (ccc, 5/17/10)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  29 May 2013 - R. Yantosca - Now set TINY = 1d-60 only for TOMAS code
!  31 May 2013 - R. Yantosca - Now pass State_Chm and then have Q point to 
!                              State_Chm%Tracers.  This is for TOMAS.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  04 Feb 2014 - R. Yantosca - Bug fix for TOMAS: call COMPUTE_F In its own
!                              separate parallel loop.  Also save the values
!                              of ISOL in the ISOL_SAVE array so that we can
!                              pass them to the parallel tracer loop.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE          :: FIRST = .TRUE.
      LOGICAL, SAVE          :: IS_Hg = .TRUE.
      INTEGER                :: I, J, K, KTOP, L, N, NDT
      INTEGER                :: IC, ISTEP, JUMP, JS, JN, NS
      INTEGER                :: IMR, JNP, NLAY, AS
      REAL*8                 :: SDT, CMOUT, ENTRN, DQ, AREA_M2
      REAL*8                 :: T0, T1, T2, T3, T4, TSUM, DELQ
      REAL*8                 :: DTCSUM(IIPAR,JJPAR,LLPAR,
     &                          Input_Opt%N_TRACERS)

      REAL*8, ALLOCATABLE, SAVE :: DSIG(:)

      ! F is the fraction of tracer lost to wet scavenging in updrafts
      REAL*8                 :: F(IIPAR,JJPAR,LLPAR,
     &                            Input_Opt%N_TRACERS)

      ! Local Work arrays
      REAL*8                 :: BMASS(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: QB(IIPAR,JJPAR)
      REAL*8                 :: MB(IIPAR,JJPAR)
      REAL*8                 :: QC(IIPAR,JJPAR) 

      ! TINY = a very small number
#if defined( TOMAS )
      REAL*8, PARAMETER      :: TINY  = 1d-60   ! Set to 1d-60 for TOMAS code
#else
      REAL*8, PARAMETER      :: TINY  = 1d-14   ! Leave as 1d-14 for std GEOS-Chem
#endif
      REAL*8, PARAMETER      :: TINY2 = 1d-30 

      ! ISOL is an index for the diagnostic arrays
      INTEGER                :: ISOL
      INTEGER                :: ISOL_SAVE(Input_Opt%N_TRACERS)

      ! QC_PRES and QC_SCAV are the amounts of tracer 
      ! preserved against and lost to wet scavenging
      REAL*8                 :: QC_PRES, QC_SCAV 

      ! DNS is the double precision value for NS
      REAL*8                 :: DNS
     
      ! Amt of Hg2, HgP scavenged out of the column (sas, bmy, 1/19/05)
      REAL*8                 :: WET_Hg2
      REAL*8                 :: WET_HgP 

      ! For fields from Input_Opt
      LOGICAL                :: LDYNOCEAN
      LOGICAL                :: LGTMM
      INTEGER                :: NC
      REAL*8                 :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER        :: Q(:,:,:,:)

      !=================================================================
      ! NFCLDMX begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy values from Input_Opt
      LDYNOCEAN = Input_Opt%LDYNOCEAN
      LGTMM     = Input_Opt%LGTMM
      NC        = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:NC)

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      Q         => State_Chm%Tracers

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'N F C L D M X  -- by S-J Lin'
         WRITE( 6, '(a)' ) 'Modified for GEOS-CHEM by Bob Yantosca'
         WRITE( 6, '(a)' ) 'Last Modification Date: 1/27/04'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

#if defined( GEOS_5 ) || defined( GEOS_FP ) || defined( MERRA )
         ! Do nothing, these are all hybrid grids
#else
         ! NOTE: We don't need to do this for GEOS-5 (bmy, 6/27/07)
         ! DSIG is the sigma-level thickness (NOTE: this assumes that
         ! we are using a pure-sigma grid.  Use new routine for fvDAS.)
         ALLOCATE( DSIG(LLPAR), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'DSIG_CONVECTION' )

         DO L = 1, LLPAR        
            DSIG(L) = GET_BP(L) - GET_BP(L+1)
         ENDDO
#endif

         ! Flag to denote if this is a mercury simulation (sas, bmy, 1/19/05)
         IS_Hg = ( Input_Opt%ITS_A_MERCURY_SIM .and. LDYNOCEAN )

         ! Reset first time flag
         FIRST = .FALSE.
      ENDIF
      
      ! Define dimensions
      IMR  = IIPAR
      JNP  = JJPAR
      NLAY = LLPAR

      ! Convection timestep [s]
      NDT  = GET_TS_CONV() * 60d0

      !=================================================================
      ! Define active convective region, from J = JS(outh) to 
      ! J = JN(orth), and to level K = KTOP. 
      !
      ! Polar regions are too cold to have moist convection.
      ! (Dry convection should be done elsewhere.)
      !
      ! We initialize the ND14 diagnostic each time we start a new 
      ! time step loop.  Only initialize DTCSUM array if the ND14 
      ! diagnostic is turned on.  This saves a quite a bit of time. 
      ! (bmy, 12/15/99)       
      !=================================================================
      IF ( ND14 > 0 ) DTCSUM = 0d0

      KTOP = NLAY - 1
      JUMP = (JNP-1) / 20
      JS   = 1 + JUMP
      JN   = JNP - JS + 1

      !=================================================================
      ! Internal time step for convective mixing is 300 sec.
      ! Doug Rotman (LLNL) says that 450 sec works just as well.       
      !=================================================================
      NS  = NDT / 300
      NS  = MAX(NS,1)
      SDT = FLOAT(NDT) / FLOAT(NS)
      DNS = DBLE( NS )

!=============================================================================
!  BMASS has units of kg/m^2 and is equivalent to AD(I,J,L) / AREA_M2
!
!   Ps - Pt (mb)| P2 - P1 | 100 Pa |  s^2  | 1  |  1 kg        kg
!  -------------+---------+--------+-------+----+--------  =  -----
!               | Ps - Pt |   mb   | 9.8 m | Pa | m^2 s^2      m^2
!
!  This is done to keep BMASS in the same units as CLDMAS * SDT
!
!  We can parallelize over levels here.  The only quantities that need to 
!  be held local are the loop counters (I, IC, J, JREF, K). (bmy, 5/2/00)
!
!  Now use routine GET_AREA_M2 from "grid_mod.f" to get surface area of
!  grid boxes in m2. (bmy, 2/4/03)
!=============================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_M2, K )
!$OMP+SCHEDULE( DYNAMIC )
      DO K = 1, NLAY
      DO J = 1, JJPAR
      DO I = 1, IMR
         AREA_M2      = GET_AREA_M2( I, J, K )
         BMASS(I,J,K) = State_Met%AD(I,J,K) / AREA_M2
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! (1)  S c a v e n g i n g   i n   C l o u d   U p d r a f t s
      !
      ! Call COMPUTE_F to compute the fraction of tracer scavenged
      ! in convective cloud updrafts.  COMPUTE_F works for both full
      ! chemistry (NSRCX == 3) and Rn-Pb-Be chemistry (NSRCX == 1)
      ! simulations.  It is best to compute the fraction of tracer
      ! scavenged at this point, outside the internal time step loop.
      ! This will avoid having to repeat the entire calculation of F
      ! for NS times.
      !
      ! ISOL, which is returned from COMPUTE_F, is the tracer index
      ! used for diagnostic arrays AD37 and AD38.  ISOL = 0 for all
      ! non-soluble tracers.
      !
      ! NOTE: We now move this into a separate parallel loop so as to
      ! ensure that TOMAS simulations result in identical results
      ! when parallelization is turned on compared to when it is
      ! turned off. (bmy, sfarina, 2/4/14)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED   )
!$OMP+PRIVATE( IC, ISOL, I, J, K ) 
!$OMP+SCHEDULE( DYNAMIC )
      DO IC = 1, NC

         ! Compute the fraction of tracer IC scavenged in updrafts
         CALL COMPUTE_F( IC, F(:,:,:,IC), ISOL, State_Met, State_Chm ) 

         ! Save ISOL values in an array for use below
         ISOL_SAVE(IC) = ISOL

         ! ND37 diagnostic -- store F only for soluble tracers
         IF ( ND37 > 0 .and. ISOL > 0 ) THEN
            DO K = 1, LD37
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               AD37(I,J,K,ISOL) = AD37(I,J,K,ISOL) + F(I,J,K,IC) 
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! (2)  T r a c e r   L o o p 
      !
      ! We now parallelize over tracers, since tracers are independent 
      ! of each other.  The parallel loop only takes effect if you 
      ! compile with the f90 "-mp" switch.  Otherwise the compiler will 
      ! interpret the parallel-processing directives as comments, and 
      ! the loop will execute on a single thread.
      !
      ! The following types of quantities must be held local for 
      ! parallelization:
      ! (1) Loop counters ( I, IC, ISTEP, J, K )
      ! (2) Scalars that are assigned values inside the tracer loop: 
      !     ( CMOUT, DELQ, ENTRN, ISOL, QC_PRES, etc. )
      ! (3) Arrays independent of tracer ( F, MB, QB, QC )
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( CMOUT, DELQ, ENTRN, I, IC, ISTEP, ISOL, J, AREA_M2, K  ) 
!$OMP+PRIVATE( MB, QB, QC, QC_PRES, QC_SCAV, T0, T1, T2, T3, T4, TSUM )
!$OMP+PRIVATE( WET_Hg2, WET_HgP                                       )
!$OMP+SCHEDULE( DYNAMIC )
      DO IC = 1, NC

         ! Get previously archived value of ISOL (bmy, 2/4/14)
         ISOL = ISOL_SAVE(IC)

         !==============================================================
         ! (3)  I n t e r n a l   T i m e   S t e p   L o o p
         !
         ! The internal time step is currently 300 seconds.  
         !==============================================================
         DO ISTEP = 1, NS

            !===========================================================
            ! (4)  B e l o w   C l o u d   B a s e  (K < 3)
            !
            ! Loop over longitude and latitude (I,J), and consider what
            ! is below the cloud base (i.e. below the 2nd sigma level).  
            !===========================================================
            DO J = JS, JN
            DO I = 1, IMR

               !========================================================
               ! (4.1) If Cloud Mass Flux exists at (I,J,2), 
               !       then compute QB.
               !
               ! QB is "weighted average" mixing ratio below the cloud 
               ! base.  QB is used to compute QC, which is the mixing 
               ! ratio of the air that moved in cumulus transport up 
               ! to the next level.  MB is the total mass of air below 
               ! the cloud base.
               !========================================================

               IF ( CLDMAS(I,J,2) .gt. TINY ) THEN 
#if   defined( GEOS_5 )

                  ! Need to replace DSIG w/ the difference 
                  ! of pressure edges for GEOS-5 (bmy, 6/27/07)
                  QB(I,J) = 
     &       ( Q(I,J,1,IC) * ( GET_PEDGE(I,J,1) - GET_PEDGE(I,J,2) )   + 
     &         Q(I,J,2,IC) * ( GET_PEDGE(I,J,2) - GET_PEDGE(I,J,3) ) ) / 
     &                       ( GET_PEDGE(I,J,1) - GET_PEDGE(I,J,3) )

#else
                  ! for GEOS-3
                  QB(I,J) = 
     &               ( Q(I,J,1,IC) * DSIG(1)   + 
     &                 Q(I,J,2,IC) * DSIG(2) ) / ( DSIG(1) + DSIG(2) )

#endif

                  ! total mass of air below the cloud base
                  MB(I,J) = BMASS(I,J,1) + BMASS(I,J,2)

!=============================================================================
!         Total mass of tracer below cloud base (i.e. MB  * QB ) +   
!         Subsidence into cloud base from above (i.e. SDT * C(2) * Q(3) ) 
!  QC =  -----------------------------------------------------------------
!             Total air mass below cloud base (i.e. MB + C(2)*Q(3) )
!=============================================================================
                  QC(I,J) = 
     &               ( MB(I,J)       * QB(I,J)       + 
     &                 CLDMAS(I,J,2) * Q(I,J,3,IC)   * SDT ) /
     &               ( MB(I,J)       + CLDMAS(I,J,2) * SDT )
     
                  !=====================================================
                  ! DQ = QB - QC is the total mass to be transported 
                  ! out of the cloud base.  Changes below cloud base 
                  ! are proportional to the background mass.
                  ! 
                  ! Subtract DQ from Q(*,*,K=1,*) and from Q(*,*,K=2,*), 
                  ! but do not make Q(*,*,K=1,*) or Q(*,*,K=2,*) < 0.
                  !=====================================================
!-----------------------------------------------------------------------------
!  Prior to 1/4/00:
!  This ensures additivity when using tagged CO tracers (bey, bmy, 1/4/00)
!                  DQ = QB(I,J) - QC(I,J)
!                  
!                  IF ( DQ .GT. Q(I,J,1,IC)   .OR. 
!     &                 DQ .GT. Q(I,J,2,IC) ) THEN
!                     Q(I,J,2,IC) = QC(I,J)
!                     Q(I,J,1,IC) = QC(I,J)
!                  ELSE
!                     Q(I,J,2,IC) = Q(I,J,2,IC) - DQ
!                     Q(I,J,1,IC) = Q(I,J,1,IC) - DQ
!                  ENDIF
!-----------------------------------------------------------------------------
                  Q(I,J,2,IC) = QC(I,J)
                  Q(I,J,1,IC) = QC(I,J)

               !========================================================
               ! If there is no Cloud mass flux, set QC = Q(K=3) 
               ! at this I,J location
               !========================================================
               ELSE
                  QC(I,J) = Q(I,J,3,IC)
               ENDIF
            ENDDO !I
            ENDDO !J

            !===========================================================
            ! (5)  A b o v e   C l o u d   B a s e
            !
            ! Loop over sigma levels K and longitude and latitude (I,J) 
            !
            ! Recall that K+1 is the level BELOW level K, and that K-1 
            ! is the level ABOVE level K.  This is the way that SJ Lin 
            ! indexes the sigma levels.
            !===========================================================
            DO K = 3,  KTOP
            DO J = JS, JN
            DO I = 1,  IMR

               ! Grid box surface area [m2]
               AREA_M2 = GET_AREA_M2( I, J, K )

!==============================================================================
!  (5.1)  M a s s   B a l a n c e   i n   C l o u d  ===>  QC(I,J)
!
!  QC_PRES = QC(I,J,K-1) * AP = amt of Qc preserved against wet scavenging
!  QC_SCAV = QC(I,J,K-1) * AS = amt of Qc lost to wet scavenging
!  CMOUT   = air mass flowing out of cloud at level K 
!  ENTRN   = Entrainment: air mass flowing into cloud at level K   
!
!  If ENTRN > 0 then compute the new value of QC(I,J):
!
!                CLDMAS(K-1) * QC_PRES  +  ENTRN(K) * Q(K)
!    QC(I,J) =  -------------------------------------------
!                      CLDMAS(I,J,K) + DTRN(I,J,K)
!
!            =   tracer mass coming in from below      (i.e. level K-1) + 
!                tracer mass coming in from this level (i.e. level K)
!               -----------------------------------------------------------
!                             total mass coming into cloud
!
!  Otherwise, preserve the previous value of QC(I,J).  This will ensure
!  that TERM1 - TERM2 is not a negative quantity (see below).
!  
!  Entrainment must be >= 0 (since we cannot have a negative flux of air
!  into the cloud).  This condition is strong enough to ensure that
!  CMOUT > 0 and will prevent floating-point exception.
!==============================================================================

                  IF ( CLDMAS(I,J,K-1) .gt. TINY ) THEN
                     CMOUT    = CLDMAS(I,J,K) + DTRN(I,J,K)
                     ENTRN    = CMOUT         - CLDMAS(I,J,K-1)
                     QC_PRES  = QC(I,J) * ( 1d0 - F(I,J,K,IC) )
                     QC_SCAV  = QC(I,J) * F(I,J,K,IC)
!#if defined ( DEBUG )
!       IF ( IC == 129 ) THEN
!          print *, I, J, K, QC_PRES, QC_SCAV, QC(I,J), F(I,J,K,IC)
!       ENDIF
!#endif

                     IF ( ENTRN .ge. 0 ) THEN
                        QC(I,J) = ( CLDMAS(I,J,K-1) * QC_PRES       + 
     &                              ENTRN           * Q(I,J,K,IC) ) / 
     &                              CMOUT
                     ENDIF   

!==============================================================================
!  (5.2)  M a s s   B a l a n c e   i n   L e v e l  ===>  Q(I,J,K,IC)
!
!  The cumulus transport above the cloud base is done as follows:
!     C_k-1  = cloud air mass flux from level k-1 to level k
!     C_k    = cloud air mass flux from level k   to level k+1
!     QC_k-1 = mixing ratio of tracer INSIDE CLOUD at level k-1
!     QC_k   = mixing ratio of tracer INSIDE CLOUD at level k
!     Q_k    = mixing ratio of tracer in level k
!     Q_k+1  = mixing ratio of tracer in level k+1
!                                
!                       |                    |
!      k+1     ^        |       Cloud        |3)      C_k * Q_k+1
!              |        |         ^          |            |
!      --------|--------+---------|----------+------------|--------
!              |        |         |          |            V
!      k     C_k        |2)   C_k * QC_k     |
!                       |                    |
!                       |                    |
!                       |         ^          |4)    C_k-1 * Q_k           
!              ^        |         |          |            |
!      --------|--------+---------|----------+------------|----------
!              |        |         |          |            |
!      k-1   C_k-1      |1) C_k-1 * QC_k-1   |            V
!                       |         * AP       |
!
!  There are 4 terms that contribute to mass flow in and out of level k:
!
!  1) C_k-1 * QC_PRES = tracer convected from level k-1 to level k 
!  2) C_k   * QC_k    = tracer convected from level k   to level k+1 
!  3) C_k   * Q_k+1   = tracer subsiding from level k+1 to level k 
!  4) C_k-1 * Q_k     = tracer subsiding from level k   to level k-1 
!
!  Therefore the change in tracer concentration is given by
!     DELQ = (Term 1) - (Term 2) + (Term 3) - (Term 4)
!
!  and Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ.  
!
!  The term T0 is the amount of tracer that is scavenged out of the box.
!  Compute that term here for the ND38 diagnostic below. (bmy, 1/27/03)
!==============================================================================
                     T0   =  CLDMAS(I,J,K-1) * QC_SCAV   
                     T1   =  CLDMAS(I,J,K-1) * QC_PRES
                     T2   = -CLDMAS(I,J,K  ) * QC(I,J       )
                     T3   =  CLDMAS(I,J,K  ) * Q (I,J,K+1,IC)
                     T4   = -CLDMAS(I,J,K-1) * Q (I,J,K,  IC)

                     TSUM = T1 + T2 + T3 + T4 

                     DELQ = ( SDT / BMASS(I,J,K) ) * TSUM

                     ! If DELQ > Q then do not make Q negative!!!
                     IF ( Q(I,J,K,IC) + DELQ < 0 ) THEN
                        DELQ = -Q(I,J,K,IC)
                     ENDIF

                     Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ
                                
                     !==================================================
                     ! ND14 Diagnostic: Upward mass flux due to wet 
                     ! convection.  The diagnostic ND14 works only for 
                     ! levels >= than 3.  This is ok for now since I 
                     ! want to have a closed budget with a box starting 
                     ! from the ground.
                     !
                     ! DTCSUM(I,J,K,IC) is the flux (kg/box/sec) in 
                     ! the box (I,J), for the tracer IC going out of 
                     ! the top of the layer K to the layer above (K+1) 
                     ! (bey, 11/10/99).                     
                     !==================================================
                     IF ( ND14 > 0 ) THEN
                        DTCSUM(I,J,K,IC) = DTCSUM(I,J,K,IC) +  
     &                       (-T2-T3) * AREA_M2 / TCVV(IC)  
                     ENDIF

                     !==================================================
                     ! ND38 Diagnostic: loss of soluble tracer to wet
                     ! scavenging in cloud updrafts [kg/s].  We must 
                     ! divide by DNS, the # of internal timesteps.
                     !==================================================
                     IF ( ( ND38 > 0 .or. LGTMM )
     &                    .and. ISOL > 0 .and. K <= LD38 ) THEN
!#if defined ( DEBUG )
!                        if(TCVV(IC) <= 0.d0) then
!                           print *, I, J, K, IC, ISOL, TCVV(IC)
!                        endif
!#endif
                        AD38(I,J,K,ISOL) = AD38(I,J,K,ISOL) +
     &                       ( T0 * AREA_M2 / ( TCVV(IC) * DNS ) )
                     ENDIF
  
                     !=================================================
                     ! Pass the amount of Hg2 and HgP lost in wet 
                     ! scavenging [kg] to "ocean_mercury_mod.f" via 
                     ! ADD_Hg2_WET and ADD_HgP_WET. We must also divide 
                     ! by DNS, the # of internal timesteps. 
                     ! (sas, bmy, eck, eds, 1/19/05, 1/6/06, 7/30/08)
                     !=================================================
                     IF ( IS_Hg .and. IS_Hg2( IC ) ) THEN

                        ! Wet scavenged Hg(II) in [kg/s]
                        WET_Hg2 = ( T0 * AREA_M2 ) / ( TCVV(IC) * DNS )

                        ! Convert [kg/s] to [kg]
                        WET_Hg2 = WET_Hg2 * NDT 

                        ! Pass to "ocean_mercury_mod.f"
                        CALL ADD_Hg2_WD( I, J, IC, WET_Hg2 )
                        CALL ADD_Hg2_SNOWPACK( I, J, IC, WET_Hg2,
     &                                         State_Met )
                     ENDIF

                     IF ( IS_Hg .and. IS_HgP( IC ) ) THEN

                        ! Wet scavenged Hg(P) in [kg/s]
                        WET_HgP = ( T0 * AREA_M2 ) / ( TCVV(IC) * DNS )

                        ! Convert [kg/s] to [kg]
                        WET_HgP = WET_HgP * NDT 

                        ! Pass to "ocean_mercury_mod.f"
                        CALL ADD_HgP_WD( I, J, IC, WET_HgP )
                        CALL ADD_Hg2_SNOWPACK( I, J, IC, WET_HgP,
     &                                         State_Met )
                     ENDIF

                  !=====================================================
                  ! No cloud transport if cloud mass flux < TINY; 
                  ! Change Qc to q
                  !=====================================================
                  ELSE                     
                     QC(I,J) = Q(I,J,K,IC)

#if   defined( GEOS_5 ) 
                     !--------------------------------------------------
                     ! FIX FOR GEOS-5 MET FIELDS!
                     !
                     ! Bug fix for the cloud base layer, which is not 
                     ! necessarily in the boundary layer, and for 
                     ! GEOS-5, there could be "secondary convection 
                     ! plumes - one in the PBL and another one not.
                     !
                     ! NOTE: T2 and T3 are the same terms as described
                     ! in the above section.
                     !
                     ! (swu, 08/13/2007)
                     !--------------------------------------------------
                     IF ( CLDMAS(I,J,K) > TINY ) THEN 

                        ! Tracer convected from K -> K+1 
                        T2   = -CLDMAS(I,J,K  ) * QC(I,J)

                        ! Tracer subsiding from K+1 -> K 
                        T3   =  CLDMAS(I,J,K  ) * Q (I,J,K+1,IC)

                        ! Change in tracer concentration
                        DELQ = ( SDT / BMASS(I,J,K) ) * (T2 + T3)

                        ! If DELQ > Q then do not make Q negative!!!
                        IF ( Q(I,J,K,IC) + DELQ < 0.0d0 ) THEN 
                            DELQ = -Q(I,J,K,IC)
                        ENDIF
  
                        ! Add change in tracer to Q array
                        Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ

                     ENDIF
#endif
                  ENDIF
               ENDDO  !I
            ENDDO     !J
            ENDDO     !K
         ENDDO        !NSTEP
      ENDDO           !IC
!$OMP END PARALLEL DO

      !=================================================================
      ! ND14 Diagnostic: Store into the CONVFLUP array.  
      ! Also divide by the number of internal timesteps (DNS).
      !================================================================= 
      IF ( ND14 > 0 ) THEN
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, IC, J, K ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO IC = 1,  NC
         DO K  = 3,  KTOP
         DO J  = JS, JN
         DO I  = 1,  IMR
            CONVFLUP(I,J,K,IC) = CONVFLUP(I,J,K,IC) +
     &                           ( DTCSUM(I,J,K,IC) / DNS )
         ENDDO
         ENDDO
         ENDDO
         ENDDO 
!$OMP END PARALLEL DO
      ENDIF

      ! Add a check to set negative values in Q to TINY2
      ! (ccc, 4/15/09)
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED   )
!$OMP+PRIVATE( N, L, J, I )
      DO N = 1,NC
      DO L = 1,LLPAR
      DO J = 1,JJPAR
      DO I = 1,IIPAR
         Q(I,J,L,N) = MAX(Q(I,J,L,N),TINY2)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer memory
      NULLIFY( Q )

      END SUBROUTINE NFCLDMX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_merra_convection
!
! !DESCRIPTION: Subroutine DO\_MERRA\_CONVECTION (formerly called NFCLDMX)
!  is S-J Lin's cumulus transport module for 3D GSFC-CTM, modified for the 
!  GEOS-Chem model.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_MERRA_CONVECTION( IDENT,    DIMINFO,  COEF, 
     &                                IDT,      OPTIONS,  AD,      
     &                                AREA_M2,  BXHEIGHT, CMFMC,   
     &                                DQRCU,    DTRAIN,   F,       
     &                                PEDGE,    PFICU,    PFLCU,   
     &                                REEVAPCN, T,        TS_DYN,   
     &                                Q,        DIAG14,   DIAG38,   
     &                                H2O2s,    SO2s,     I,        
     &                                J,        RC,       State_Met,
     &                                State_Chm                      )
!
! !USES:
!
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_SNOWPACK
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_WD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_HgP_WD
      USE ERROR_MOD,          ONLY : IT_IS_NAN, IT_IS_FINITE
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP  ! hma Nov 3, debug
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IS_Hg2
      USE TRACERID_MOD,       ONLY : IS_HgP
      USE WETSCAV_MOD,        ONLY : WASHOUT
      USE WETSCAV_MOD,        ONLY : LS_K_RAIN
      USE WETSCAV_MOD,        ONLY : LS_F_PRIME
!
! !INPUT PARAMETERS:
!  
      TYPE(SPEC_2_TRAC), INTENT(IN) :: COEF        ! Obj w/ spec <-> trac map
      TYPE(GC_DIMS),     INTENT(IN) :: DIMINFO     ! Obj w/ array dimensions
      TYPE(ID_TRAC),     INTENT(IN) :: IDT         ! Obj w/ tracer ID flags
      TYPE(GC_OPTIONS),  INTENT(IN) :: OPTIONS     ! Obj w/ logical switches
      TYPE(MetState),    INTENT(IN) :: State_Met   ! Meteorology State object
      REAL*8,            INTENT(IN) :: AD(:)       ! Air mass [kg]    
      REAL*8,            INTENT(IN) :: AREA_M2     ! Surface area [m2]
      REAL*8,            INTENT(IN) :: BXHEIGHT(:) ! Box height [m]
      REAL*8,            INTENT(IN) :: CMFMC(:)    ! Cloud mass flux [kg/m2/s] 
      REAL*8,            INTENT(IN) :: DQRCU(:)    ! Precip production rate:
                                                   !  convective [kg/kg/s]
      REAL*8,            INTENT(IN) :: DTRAIN(:)   ! Detrainment flux [kg/m2/s]
      REAL*8,            INTENT(IN) :: F(:,:)      ! Fraction of soluble tracer
                                                   !  for updraft scavenging 
                                                   !  [unitless].  ! This is 
                                                   !  computed by routine 
                                                   !  COMPUTE_UPDRAFT_FSOL
      REAL*8,            INTENT(IN) :: PEDGE(:)    ! P @ level box edges [hPa]
      REAL*8,            INTENT(IN) :: PFICU(:)    ! Dwnwd flux of convective
                                                   !  ice precip [kg/m2/s]
      REAL*8,            INTENT(IN) :: PFLCU(:)    ! Dwnwd flux of convective
                                                   !  liquid precip [kg/m2/s]
      REAL*8,            INTENT(IN) :: REEVAPCN(:) ! Evap of precip'ing conv.
                                                   !  condensate [kg/kg/s]
      REAL*8,            INTENT(IN) :: T(:)        ! air temperature [K]
      REAL*8,            INTENT(IN) :: TS_DYN      ! Dynamic timestep [min]
      INTEGER,           INTENT(IN) :: I, J        ! Lon & lat indices
!                                                  
! !INPUT/OUTPUT PARAMETERS:                        
!        
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(GC_IDENT), INTENT(INOUT) :: IDENT       ! Obj w/ info from ESMF etc.
      REAL*8,         INTENT(INOUT) :: H2O2s(:)
      REAL*8,         INTENT(INOUT) :: SO2s(:)  
      REAL*8,         INTENT(INOUT) :: Q(:,:)      ! Tracer conc. [mol/mol]  
!                                                  
! !OUTPUT PARAMETERS:                              
!                     
      REAL*8,           INTENT(OUT) :: DIAG14(:,:) ! Array for ND14 diagnostic
      REAL*8,           INTENT(OUT) :: DIAG38(:,:) ! Array for ND38 diagnostic
      INTEGER,          INTENT(OUT) :: RC          ! Return code
!
! !REMARKS:
!  Reference:
!  ============================================================================
!  Lin, SJ.  "Description of the parameterization of cumulus transport
!     in the 3D Goddard Chemistry Transport Model, NASA/GSFC, 1996.
!                                                                             .
!  Unit conversion for BMASS:
!
!      Ps - Pt (mb)| P2 - P1 | 100 Pa |  s^2  | 1  |  1 kg        kg
!     -------------+---------+--------+-------+----+--------  =  -----
!                  | Ps - Pt |   mb   | 9.8 m | Pa | m^2 s^2      m^2
!
!                                                                             .
!  NOTE: We are passing I & J down to this routine so that it can call the
!  proper code from "mercury_mod.f".  Normally, we wouldn't pass I & J as
!  arguments to columnized code.  This prevents rewriting the mercury_mod.f
!  routines ADD_Hg2_
!
! !REVISION HISTORY:
!  15 Jul 2009 - R. Yantosca - Columnized and cleaned up.
!                            - CLDMAS renamed to CMFMC and DTRN renamed
!                              to DTRAIN for consistency w/ GEOS-5.
!  17 Jul 2009 - R. Yantosca - Now do unit conversion of Q array from
!                              [kg] --> [v/v] and vice versa internally
!  14 Dec 2009 - R. Yantosca - Now remove internal unit conversion, since
!                              Q now comes in as [mol/mol] (=[v/v]) from the
!                              calling routine.
!  14 Dec 2009 - R. Yantosca - Remove COEF from the argument list
!  06 May 2010 - R. Yantosca - Now add IDENT via the argument list
!  29 Sep 2010 - R. Yantosca - Modified for MERRA met fields
!  05 Oct 2010 - R. Yantosca - Now pass COEF via the argument list
!  05 Oct 2010 - R. Yantosca - Attach ND14 and ND38 diagnostics
!  15 Oct 2010 - H. Amos     - Added BXHEIGHT and T as arguments
!  15 Oct 2010 - R. Yantosca - Added I, J, H2O2s and SO2s as arguments
!  15 Oct 2010 - H. Amos     - Added scavenging below cloud base
!  06 Apr 2011 - M.Fu, H.Amos- Bug fix: make sure washout adheres to the same
!                              algorithm as in the wet deposition code.
!  27 Jul 2011 - R. Yantosca - Declare CLDBASE as INTEGER to avoid PGI errors
!  16 Aug 2011 - J. Fisher   - Bug fix: use IS_Hg2() and IS_HgP to test if 
!                              a tracer is Hg2 or HgP (for tagged species)
!  16 Aug 2011 - J. Fisher   - Now use WETLOSS instead of T0_SUM in the ND38 
!                              diagnostic below the cloud.  Using T0_SUM leads 
!                              us to over-count the tracer scavenged out of 
!                              the column.
!  22 Oct 2012 - R. Yantosca - Now reference Headers/gigc_errcode_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  31 May 2013 - R. Yantosca - Now pass State_Chm to WASHOUT
!  05 Sep 2013 - R. Yantosca - Bug fix: DT is apparently undefined, but still
!                              passed to WASHOUT.  Use SDT instead.  This
!                              avoids a floating-point error.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER :: TINYNUM = 1d-14 
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL           :: AER,         IS_Hg 
      INTEGER           :: IC,          ISTEP,    K,       KTOP
      INTEGER           :: NC,          NDT,      NLAY,    NS    
      INTEGER           :: CLDBASE
      REAL*8            :: CMFMC_BELOW, ALPHA,    ALPHA2
      REAL*8            :: CMOUT,       DELQ,     DQ,      DNS
      REAL*8            :: ENTRN,       QC,       QC_PRES, QC_SCAV 
      REAL*8            :: SDT,         T0,       T0_SUM,  T1    
      REAL*8            :: T2,          T3,       T4,      TCVV   
      REAL*8            :: TCVV_DNS,    TSUM,     DT_DNS
      REAL*8            :: LOST,        GAINED,   WETLOSS, MASS_WASH
      REAL*8            :: MASS_NOWASH, QDOWN,    DT,      F_WASHOUT
      REAL*8            :: K_RAIN,      WASHFRAC, WET_Hg2, WET_HgP 
      REAL*8            :: MB,          QB

      ! Arrays
      REAL*8            :: BMASS ( DIMINFO%L_COLUMN  )   
      REAL*8            :: PDOWN ( DIMINFO%L_COLUMN  )
      REAL*8            :: QB_NUM( DIMINFO%N_TRACERS )

      !========================================================================
      ! (0)  I n i t i a l i z a t i o n
      !========================================================================
      
      ! Put this routine name on error trace stack
      IDENT%LEV               = IDENT%LEV + 1
      IDENT%I_AM( IDENT%LEV ) = 'DO_MERRA_CONVECTION'

      ! # of levels and # of tracers
      NLAY = DIMINFO%L_COLUMN
      NC   = DIMINFO%N_TRACERS

      ! Safety check #1: Invalid NLAY w/r/t Q
      IF ( NLAY < 1 .or. NLAY > SIZE( Q, 1 ) ) THEN
         WRITE( IDENT%ERRMSG, 200 ) NLAY, SIZE( Q, 1 )
         RC = GIGC_FAILURE
         RETURN
      ENDIF

      ! Safety check #2: Invalid NC w/r/t Q
      IF ( NC < 1 .or. NC > SIZE( Q, 2 ) ) THEN
         WRITE( IDENT%ERRMSG, 200 ) NC, SIZE( Q, 2 )
         RC = GIGC_FAILURE
         RETURN
      ENDIF

      ! Safety check #3: Invalid NLAY & NC w/r/t F
      IF ( NLAY > SIZE( F, 1 ) .or. NC > SIZE( F, 2 ) ) THEN
         WRITE( IDENT%ERRMSG, 220 ) NLAY, SIZE( F, 1 ),
     &                              NC,   SIZE( F, 2 )
         RC = GIGC_FAILURE
         RETURN
      ENDIF

      ! Formats
 200  FORMAT( 'Error: NLAY = ',  i5, ' but SIZE( Q,1 ) = ', i5 )
 210  FORMAT( 'Error: NC = ',    i5, ' but SIZE( Q,2 ) = ', i5 )
 220  FORMAT( 'Error: NLAY = ',  i5, ' but SIZE( F,1 ) = ', i5,
     &        '  or NC = ',      i5, ' but SIZE( F,2 ) = ', i5 )

      ! Top level for convection
      KTOP   = NLAY - 1

      ! Convection timestep [s]
      NDT    = TS_DYN * 60d0

      ! Internal time step for convective mixing is 300 sec.
      ! Doug Rotman (LLNL) says that 450 sec works just as well.       
      NS     = NDT / 300
      NS     = MAX( NS, 1 )
      SDT    = DBLE( NDT ) / DBLE( NS )
      DNS    = DBLE( NS )
      DT_DNS = NDT * DNS
      
      !-----------------------------------------------------------------
      ! Determine location of the cloud base, which is the level where
      ! we start to have non-zero convective precipitation formation
      !-----------------------------------------------------------------

      ! Minimum value of cloud base is the surface level
      CLDBASE = 1

      ! Find the cloud base
      DO K = 1, NLAY
         IF ( DQRCU(K) > 0d0 ) THEN
            CLDBASE = K
            EXIT
         ENDIF
      ENDDO

      !-----------------------------------------------------------------
      ! Compute PDOWN and BMASS
      !-----------------------------------------------------------------
      DO K = 1, NLAY

         ! PDOWN is the convective precipitation leaving each
         ! box [cm3 H2O/cm2 air/s].This accounts for the 
         ! contribution from both liquid & ice precip
         PDOWN(K) = ( ( PFLCU(K) / 1000d0 ) 
     &            +   ( PFICU(K) /  917d0 ) ) * 100d0

         ! BMASS has units of kg/m^2 and is equivalent to AD(K) / AREA_M2
         ! This is done to keep BMASS in the same units as CMFMC * SDT
         BMASS(K) = AD(K)    / AREA_M2

      ENDDO

      !-----------------------------------------------------------------
      ! Compute MB (the total mass of air below the cloud base)
      ! and QB_NUM (the numerator of QB, the weighted avg mixing ratio
      ! below the cloud base)   
      !-----------------------------------------------------------------

      ! Zero variables
      MB      = 0d0
      QB_NUM  = 0d0

      ! Loop over up to just below the cloud base
      DO K  = 1, CLDBASE-1
         
         ! Total mass of air below the cloud base
         MB        = MB + BMASS(K)

         ! numerator of QB
         QB_NUM(:) = QB_NUM(:) + Q(K,:) * ( PEDGE(K) - PEDGE(K+1) )

      ENDDO

      ! Is this a Hg simulation?
      IS_Hg = ( OPTIONS%USE_Hg .and. OPTIONS%USE_Hg_DYNOCEAN )

      !========================================================================
      ! (1)  T r a c e r   L o o p 
      !========================================================================
      DO IC = 1, NC

         ! Initialize
         DIAG14(:,IC) = 0d0                            ! ND14 diag array      
         DIAG38(:,IC) = 0d0                            ! ND38 diag array
         TCVV         = MW_AIR / COEF%MOLWT_KG(IC)     ! Air MW/tracer MW
         TCVV_DNS     = TCVV * DNS                     ! TCVV * DNS

         !=====================================================================
         ! (2)  I n t e r n a l   T i m e   S t e p   L o o p
         !=====================================================================
         DO ISTEP = 1, NS

            ! Initialize
            QC     = 0d0                      
            T0_SUM = 0d0  

            !----------------------------------------------------------
            ! B e l o w   C l o u d   B a s e   (K < CLDBASE) 
            !
            ! QB is the "weighted avg" mixing ratio below the cloud
            ! base. QB is used to compute QC, which is the mixing 
            ! ratio of the air that moved in cumulus transport up to
            ! the next level. MB is the total mass of air below the
            ! the cloud base
            !-----------------------------------------------------------

            ! We need to make this a nested IF statement so that we don't
            ! get an out-of-bounds error when CLDBASE=1 (bmy, 11/18/10)
            IF ( CLDBASE > 1 ) THEN
               
               IF ( CMFMC(CLDBASE-1) > TINYNUM ) THEN
 
                  !-----------------------------------------------------
                  ! %%% Non-negligible Cloud mass flux %%% 
                  !-----------------------------------------------------

                  ! Compute the weighted avg mixing ratio below
                  ! the cloud base
                  QB = QB_NUM(IC) / ( PEDGE(1) - PEDGE(CLDBASE) )

                  !        Total mass of tracer below cloud base  +   
                  !        Subsidence into cloud base from above  
                  ! QC =  --------------------------------------------
                  !            Total air mass below cloud base 
                  !
                  QC = ( MB*QB + CMFMC(CLDBASE-1) * 
     &                           Q(CLDBASE,IC)    * SDT  ) /
     &                 ( MB    + CMFMC(CLDBASE-1) * SDT  ) 


                  ! Copy QC to all levels of the tracer array Q
                  ! that are below the cloud base level
                  Q(1:CLDBASE-1,IC) = QC

               ELSE                  

                  !-----------------------------------------------------
                  ! %%% Negligible cloud mass flux %%% 
                  !-----------------------------------------------------

                  ! When CMFMC is negligible, then set to the tracer
                  ! concentration at the cloud base level
                  QC = Q(CLDBASE,IC)
   
               ENDIF

            ELSE
               
               !-----------------------------------------------------
               ! If the cloud base happens at level 1, then just
               ! set QC to the tracer at the surface level
               !-----------------------------------------------------
               QC = Q(CLDBASE,IC) 

            ENDIF

            !==================================================================
            ! (3)  A b o v e   C l o u d   B a s e
            !==================================================================
            DO K = CLDBASE, KTOP
               
               ! Initialize
               ALPHA   = 0d0
               ALPHA2  = 0d0
               CMOUT   = 0d0
               ENTRN   = 0d0
               QC_PRES = 0d0
               QC_SCAV = 0d0
               
               ! CMFMC_BELOW is the air mass [kg/m2/s] coming into the
               ! grid box (K) from the box immediately below (K-1).
               IF ( K == 1 ) THEN
                  CMFMC_BELOW = 0d0
               ELSE
                  CMFMC_BELOW = CMFMC(K-1)
               ENDIF

               ! If we have a nonzero air mass flux coming from 
               ! grid box (K-1) into (K) ...
               IF ( CMFMC_BELOW > TINYNUM ) THEN

                  !------------------------------------------------------------
                  ! (3.1)  M a s s   B a l a n c e   i n   C l o u d
                  !
                  ! F(K,IC) = fraction of tracer IC in level K that is 
                  !           available for wet-scavenging by cloud updrafts.  
                  !
                  ! If ENTRN > 0 then compute the new value of QC:
                  !
                  !      tracer mass from below      (i.e. level K-1) + 
                  !      tracer mass from this level (i.e. level K)
                  !  = -----------------------------------------------------
                  !             total mass coming into cloud
                  !
                  ! Otherwise, preserve the previous value of QC.  This will 
                  ! ensure that TERM1 - TERM2 is not a negative quantity (see 
                  ! below).
                  !  
                  ! Entrainment must be >= 0 (since we cannot have a negative 
                  ! flux of air into the cloud).  This condition is strong 
                  ! enough to ensure that CMOUT > 0 and will prevent floating-
                  ! point exception.
                  !------------------------------------------------------------

                  ! Air mass flowing out of cloud at grid box (K)
                  CMOUT   = CMFMC(K) + DTRAIN(K)

                  ! Air mass flowing into cloud at grid box (K)
                  ENTRN   = CMOUT - CMFMC_BELOW

                  ! Amount of QC preserved against scavenging
                  QC_PRES = QC * ( 1d0 - F(K,IC) )

                  ! Amount of QC lost to scavenging
                  ! QC_SCAV = 0 for non-soluble tracer
                  QC_SCAV = QC * F(K,IC)

                  ! - - - - - - - - FOR SOLUBLE TRACERS ONLY - - - - - - - - - 
                  IF ( QC_SCAV > 0d0 ) THEN 

                     ! The fraction ALPHA is the fraction of raindrops that 
                     ! will re-evaporate soluble tracer while falling from 
                     ! grid box K+1 down to grid box K.  Avoid div-by-zero.

                     ! Initialize 
                     ALPHA = 0d0
                     
                     IF ( PDOWN(K+1)  > TINYNUM ) THEN 

                        ! %%%% CASE 1 %%%%
                        ! Partial re-evaporation. Less precip is leaving 
                        ! the grid box then entered from above.
                        IF ( PDOWN(K+1) > PDOWN(K) .AND.
     &                       PDOWN(K)   > TINYNUM        ) THEN

                           ! Define ALPHA, the fraction of raindrops that 
                           ! re-evaporate when falling from grid box
                           ! (I,J,L+1) to (I,J,L)
                           ALPHA = ( REEVAPCN(K) * AD(K)           )
     &                           / ( PDOWN(K+1)  * AREA_M2 * 10d0  )

                           ! Restrict ALPHA to be less than 1 
                           ! (>1 is unphysical)  (hma, 24-Dec-2010)  
                           IF ( ALPHA > 1d0 ) THEN 
                              ALPHA = 1d0 
                           ENDIF

                           ! We assume that 1/2 of the soluble tracer w/in 
                           ! the raindrops actually gets resuspended into 
                           ! the atmosphere
                           ALPHA2   = ALPHA * 0.5d0

                        ENDIF

                        ! %%%% CASE 2 %%%%
                        ! Total re-evaporation. Precip entered from above, 
                        ! but no precip is leaving grid box (ALPHA = 2 so 
                        ! that  ALPHA2 = 1)
                        IF ( PDOWN(K) < TINYNUM ) THEN
                           ALPHA2 = 1d0
                        ENDIF
                        
                     ENDIF
                    
                     ! The resuspension takes 1/2 the amount of the scavenged
                     ! aerosol (QC_SCAV) and adds that back to QC_PRES ...
                     QC_PRES  = QC_PRES + ( ALPHA2 * QC_SCAV )
                    
                     ! ... then we decrement QC_SCAV accordingly
                     QC_SCAV  = QC_SCAV * ( 1d0    - ALPHA2     )

                  ENDIF
                  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

                  ! Update QC taking entrainment into account
                  ! Prevent div by zero condition
                  IF ( ENTRN >= 0d0 .and. CMOUT > 0d0 ) THEN
                     QC   = ( CMFMC_BELOW * QC_PRES   + 
     &                       ENTRN       * Q(K,IC) ) / CMOUT
                  ENDIF

                  !------------------------------------------------------------
                  ! (3.2)  M a s s   B a l a n c e   i n   L e v e l  ==> Q
                  !
                  ! Terminology:
                  !
                  !  C_k-1   = cloud air mass flux from level k-1 to level k
                  !  C_k     = cloud air mass flux from level k   to level k+1
                  !  QC_k-1  = mixing ratio of tracer INSIDE CLOUD at level k-1
                  !  QC_k    = mixing ratio of tracer INSIDE CLOUD at level k
                  !  Q_k     = mixing ratio of tracer in level k
                  !  Q_k+1   = mixing ratio of tracer in level k+1
                  ! 
                  ! For convenience we denote:
                  !
                  !  QC_SCAV = Amount of tracer wet-scavenged in updrafts
                  !          = QC_k-1 * F(k,IC)
                  !
                  !  QC_PRES = Amount of tracer preserved against
                  !            wet-scavenging in updrafts
                  !          = QC_k-1 * ( 1 - F(k,IC) )
                  !
                  ! Where F(k,IC) is the fraction of tracer IC in level k
                  ! that is available for wet-scavenging by cloud updrafts.
                  ! F(k,IC) is computed by routine COMPUTE_UPDRAFT_FSOL
                  ! and passed to this routine as an argument.
                  !
                  ! The cumulus transport above the cloud base is done as 
                  ! follows:
                  !  
                  !                 ||///////////////////||
                  !                 ||//// C L O U D ////||
                  !                 ||                   ||
                  !   k+1     ^     ||         ^         ||3)   C_k * Q_k+1
                  !           |     ||         |         ||         |
                  !   --------|-----++---------|---------++---------|--------
                  !           |     ||         |         ||         |
                  !   k      C_k    ||2)   C_k * QC_k    ||         V
                  !                 ||                   ||
                  !                 ||                   ||
                  !           ^     ||         ^         ||4)   C_k-1 * Q_k 
                  !           |     ||         |         ||         |
                  !   --------|-----++---------|---------++---------|--------
                  !           |     ||         |         ||         |
                  !   k-1   C_k-1   ||1) C_k-1 * QC_k-1  ||         V
                  !                 ||         * (1 - F) ||
                  !                 ||                   ||
                  !                 ||//// C L O U D ////||
                  !                 ||///////////////////||
                  !
                  ! There are 4 terms that contribute to mass flow in 
                  ! and out of level k:
                  !
                  ! 1) C_k-1 * QC_PRES = tracer convected from k-1 to k 
                  ! 2) C_k   * QC_k    = tracer convected from k   to k+1 
                  ! 3) C_k   * Q_k+1   = tracer subsiding from k+1 to k 
                  ! 4) C_k-1 * Q_k     = tracer subsiding from k   to k-1 
                  !
                  ! Therefore the change in tracer concentration is given by
                  !
                  !    DELQ = (Term 1) - (Term 2) + (Term 3) - (Term 4)
                  !
                  ! and Q(K,IC) = Q(K,IC) + DELQ.  
                  !
                  ! The term T0 is the amount of tracer that is scavenged 
                  ! out of the box.
                  !------------------------------------------------------------
                  T0      =  CMFMC_BELOW * QC_SCAV   
                  T1      =  CMFMC_BELOW * QC_PRES
                  T2      = -CMFMC(K  )  * QC
                  T3      =  CMFMC(K  )  * Q(K+1,IC)
                  T4      = -CMFMC_BELOW * Q(K,  IC)
                 
                  TSUM    = T1 + T2 + T3 + T4 
                  
                  DELQ    = ( SDT / BMASS(K) ) * TSUM

                  ! If DELQ > Q then do not make Q negative!!!
                  IF ( Q(K,IC) + DELQ < 0 ) THEN
                     DELQ = -Q(K,IC)
                  ENDIF

                  ! Increment the tracer array [v/v]
                  Q(K,IC) = Q(K,IC) + DELQ

                  ! check for infinity
                  IF ( .not. IT_IS_FINITE( Q(K,IC) ) ) THEN
                     CALL GEOS_CHEM_STOP
                  ENDIF

                  ! Return if we encounter NaN
                  IF ( IT_IS_NAN( Q(K,IC) ) ) THEN 
                     WRITE( 6, 250 )
                     WRITE( 6, 255 ) K, IC, Q(K,IC) 
 250                 FORMAT( 'NaN encountered in DO_MERRA_CONVECTION!' )
 255                 FORMAT( 'K, IC, Q(K,IC): ', 2i4, 1x, es13.6 )
                     RC = GIGC_FAILURE
                     RETURN
                  ENDIF

                  ! Archive T0 for use in the next section
                  T0_SUM  = T0_SUM + T0
 
                  !------------------------------------------------------------
                  ! (3.3)  N D 1 4   D i a g n o s t i c
                  !
                  ! Archive upward mass flux due to wet convection.  
                  ! DTCSUM(K,IC) is the flux [kg/sec] in the box (I,J), 
                  ! for the tracer IC going out of the top of the layer K 
                  ! to the layer above (K+1)  (bey, 11/10/99). 
                  !------------------------------------------------------------
                  IF ( OPTIONS%USE_DIAG14 ) THEN
                     DIAG14(K,IC) = DIAG14(K,IC) 
     &                            + ( ( -T2-T3 ) * AREA_M2 / TCVV_DNS )
                  ENDIF

                  !------------------------------------------------------------
                  ! (3.4)  N D 3 8   D i a g n o s t i c
                  !
                  ! Archive the loss of soluble tracer to wet scavenging in 
                  ! cloud updrafts [kg/s].  We must divide by DNS, the # of 
                  ! internal timesteps.
                  !------------------------------------------------------------
                  IF ( OPTIONS%USE_DIAG38 .and. F(K,IC) > 0d0 ) THEN
                     DIAG38(K,IC) = DIAG38(K,IC)
     &                            + ( T0 * AREA_M2 / TCVV_DNS )

                     ! check for infinity (added by hma, 20101117)
                     IF ( .not. IT_IS_FINITE( DIAG38(K,IC) ) ) THEN
                        PRINT*, 'DIAG38 IS INFINITY at K = ', K
                        CALL GEOS_CHEM_STOP
                     ENDIF
                  ENDIF
                  
               ELSE

                  !------------------------------------------------------------
                  ! (3.5)  N o   C l o u d   M a s s   F l u x   B e l o w 
                  !------------------------------------------------------------

                  ! If there is no cloud mass flux coming from below
                  ! just set QC to the tracer concentration at this level
                  QC = Q(K,IC)
                  
                  ! Bug fix for the cloud base layer, which is not necessarily
                  ! in the boundary layer, and for MERRA, there could be 
                  ! "secondary convection" plumes - one in the PBL and another 
                  ! one not.  NOTE: T2 and T3 are the same terms as described 
                  ! in the above section.  (swu, 08/13/2007)
                  IF ( CMFMC(K) > TINYNUM ) THEN 

                     ! Tracer convected from K -> K+1 
                     T2   = -CMFMC(K) * QC

                     ! Tracer subsiding from K+1 -> K 
                     T3   =  CMFMC(K) * Q(K+1,IC)

                     ! Change in tracer concentration
                     DELQ = ( SDT / BMASS(K) ) * (T2 + T3)

                     ! If DELQ > Q then do not make Q negative!!!
                     IF ( Q(K,IC) + DELQ < 0.0d0 ) THEN 
                        DELQ = -Q(K,IC)
                     ENDIF
  
                     ! Add change in tracer to Q array
                     Q(K,IC) = Q(K,IC) + DELQ

                  ENDIF
               ENDIF
            ENDDO 
            
            !==================================================================
            ! (4)  B e l o w   C l o u d   B a s e
            !==================================================================
            DO K = CLDBASE-1, 1, -1

               ! Initialize
               QDOWN       = 0d0
               F_WASHOUT   = 0d0
               WASHFRAC    = 0d0
               ALPHA       = 0d0
               ALPHA2      = 0d0
               GAINED      = 0d0
               WETLOSS     = 0d0
               LOST        = 0d0
               MASS_WASH   = 0d0
               MASS_NOWASH = 0d0
               K_RAIN      = 0d0
               
               ! Check if...
               ! (1) there is precip coming into box (I,J,K) from (I,J,K+1)
               ! (2) there is re-evaporation happening in grid box (I,J,K)
               ! (3) there is tracer to re-evaporate
               IF ( PDOWN(K+1)  > 0 .and.   
     &              REEVAPCN(K) > 0 .and. 
     &              T0_SUM      > 0        ) THEN

                  ! Compute F_WASHOUT, the fraction of grid box (I,J,L)
                  ! experiencing washout 
               
                  ! Convert PDOWN the downward flux of precip leaving grid
                  ! box (K+1) from [cm3 H20/cm2 area/s] to [cm3 H20/cm3 air/s]
                  QDOWN = PDOWN(K+1) / ( BXHEIGHT(K+1) * 100d0  )

                  ! Compute K_RAIN and F_WASHOUT based on the flux of precip 
                  ! leaving grid box (K+1). 
                  K_RAIN   = LS_K_RAIN( QDOWN )
                  F_WASHOUT= LS_F_PRIME( QDOWN, K_RAIN ) 
                  
                  ! Call WASHOUT to compute the fraction of tracer lost
                  ! to washout in grid box (I,J,K)
                  CALL WASHOUT( I, J, K,   IC,       BXHEIGHT(K), 
     &                          T(K),      QDOWN,    SDT,   
     &                          F_WASHOUT, H2O2s(K), SO2s(K), 
     &                          WASHFRAC,  AER,      State_Met,
     &                          State_Chm                         )

                  ! Check if the tracer is an aerosol or not
                  IF ( AER == .TRUE. ) THEN
 
                     !---------------------------------------------------------
                     ! Washout of aerosol tracers
                     ! This is modeled as a kinetic process
                     !---------------------------------------------------------

                     ! Define ALPHA, the fraction of raindrops that 
                     ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
                     ALPHA   = ( REEVAPCN(K) * AD(K)          )
     &                       / ( PDOWN(K+1)  * AREA_M2 * 10d0 )

                     ! ALPHA2 is the fraction of the rained-out aerosols
                     ! that gets resuspended in grid box (I,J,L)
                     ALPHA2  = 0.5d0 * ALPHA

                     ! GAINED is the rained out aerosol coming down from 
                     ! grid box (I,J,L+1) that will evaporate and re-enter 
                     ! the atmosphere in the gas phase in grid box (I,J,L).
                     GAINED  = T0_SUM * ALPHA2

                     ! Amount of aerosol lost to washout in grid box
                     WETLOSS = Q(K,IC) * WASHFRAC - GAINED

                     ! LOST is the rained out aerosol coming down from
                     ! grid box (I,J,L+1) that will remain in the liquid
                     ! phase in grid box (I,J,L) and will NOT re-evaporate.
                     LOST    = T0_SUM - GAINED

                     ! Update T0_SUM, the total amount of scavenged
                     ! tracer that will be passed to the grid box below
                     T0_SUM = T0_SUM + WETLOSS  


                  ELSE

                     !---------------------------------------------------------
                     ! Washout of non-aerosol tracers
                     ! This is modeled as an equilibrium process
                     !---------------------------------------------------------

                     ! MASS_NOWASH is the amount of non-aerosol tracer in 
                     ! grid box (I,J,L) that is NOT available for washout.
                     MASS_NOWASH = ( 1d0 - F_WASHOUT ) * Q(K,IC)
               
                     ! MASS_WASH is the total amount of non-aerosol tracer
                     ! that is available for washout in grid box (I,J,L).
                     ! It consists of the mass in the precipitating
                     ! part of box (I,J,L), plus the previously rained-out
                     ! tracer coming down from grid box (I,J,L+1).
                     ! (Eq. 15, Jacob et al, 2000).
                     MASS_WASH   = ( F_WASHOUT * Q(K,IC) ) + T0_SUM

                     ! WETLOSS is the amount of tracer mass in 
                     ! grid box (I,J,L) that is lost to washout.
                     ! (Eq. 16, Jacob et al, 2000)
                     WETLOSS     = MASS_WASH * WASHFRAC - T0_SUM 
                     
                     ! The tracer left in grid box (I,J,L) is what was
                     ! in originally in the non-precipitating fraction 
                     ! of the box, plus MASS_WASH, less WETLOSS. 
                     Q(K,IC)     = Q(K,IC) - WETLOSS
                     
                     ! Updated T0_SUM, the total scavenged tracer
                     ! that will be passed to the grid box below
                     T0_SUM      = T0_SUM + WETLOSS 

                  ENDIF

                  !------------------------------------------------------------
                  ! N D 1 4   D i a g n o s t i c
                  !
                  ! Archive upward mass flux due to wet convection.  
                  ! DTCSUM(K,IC) is the flux [kg/sec] in the box (I,J), 
                  ! for the tracer IC going out of the top of the layer K 
                  ! to the layer above (K+1)  (bey, 11/10/99). 
                  !------------------------------------------------------------
                  IF ( OPTIONS%USE_DIAG14 ) THEN
                     DIAG14(K,IC) = DIAG14(K,IC) 
     &                            + ( ( -T2-T3 ) * AREA_M2 / TCVV_DNS )
                  ENDIF

                  !------------------------------------------------------------
                  !  N D 3 8   D i a g n o s t i c
                  !
                  ! Archive the loss of soluble tracer to wet scavenging in 
                  ! cloud updrafts [kg/s].  We must divide by DNS, the # of 
                  ! internal timesteps.
                  !------------------------------------------------------------
                  IF ( OPTIONS%USE_DIAG38 .and. F(K,IC) > 0d0 ) THEN
                     DIAG38(K,IC) = DIAG38(K,IC)
     &                            + ( WETLOSS * AREA_M2 / TCVV_DNS )
                  ENDIF
                
                  ! CHECK for infinity (added by hma, 20101117)
                  IF ( .not. IT_IS_FINITE( DIAG38(K,IC) ) ) THEN
                     PRINT*, 'DIAG38 IS INFINITY at K = ', K
                     CALL GEOS_CHEM_STOP
                  ENDIF
               ENDIF
            ENDDO  

            !==================================================================
            ! (5)  M e r c u r y   O c e a n   M o d e l   A r c h i v a l
            !
            ! Pass the amount of Hg2 and HgP lost in wet  scavenging [kg] to 
            ! "ocean_mercury_mod.f" via ADD_Hg2_WET and ADD_HgP_WET.   We must 
            ! also divide  by DNS, the # of internal timesteps. 
            ! (sas, bmy, eck, eds, 1/19/05, 1/6/06, 7/30/08)
            !==================================================================
 
            !--------------------------------------
            ! Hg2 
            !--------------------------------------
            IF ( IS_Hg .and. IS_Hg2( IC ) ) THEN

               ! Wet scavenged Hg(II) in [kg/s]
               WET_Hg2 = ( T0_SUM * AREA_M2 / TCVV_DNS )
               
               ! Convert [kg/s] to [kg]
               WET_Hg2 = WET_Hg2 * NDT 

               ! Pass to "ocean_mercury_mod.f"
               CALL ADD_Hg2_WD      ( I, J, IC, WET_Hg2 )
               CALL ADD_Hg2_SNOWPACK( I, J, IC, WET_Hg2, State_Met )
            ENDIF

            !--------------------------------------
            ! HgP
            !--------------------------------------
            IF ( IS_Hg .and. IS_HgP( IC ) ) THEN

               ! Wet scavenged Hg(P) in [kg/s]
               WET_HgP = ( T0_SUM * AREA_M2 / TCVV_DNS )

               ! Convert [kg/s] to [kg]
               WET_HgP = WET_HgP * NDT 

               ! Pass to "ocean_mercury_mod.f"
               CALL ADD_HgP_WD      ( I, J, IC, WET_HgP )
               CALL ADD_Hg2_SNOWPACK( I, J, IC, WET_HgP, State_Met )
            ENDIF

         ENDDO 
      ENDDO      

      !================================================================
      ! Succesful return!
      !================================================================

      ! Set error code to success
      RC                      = GIGC_SUCCESS

      ! Remove this routine name from error trace stack
      IDENT%I_AM( IDENT%LEV ) = ''
      IDENT%LEV               = IDENT%LEV - 1

      END SUBROUTINE DO_MERRA_CONVECTION
!EOC
      END MODULE CONVECTION_MOD
