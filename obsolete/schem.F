!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: schem
!
! !DESCRIPTION: Subroutine SCHEM performs simplified stratospheric chemistry, 
!  which means only reactions with OH and photolysis are considered.  The 
!  production and loss of CO and NOy in the stratosphere are taken from Dylan 
!  Jones' 2-D model. 
!\\
!\\
! !INTERFACE:
!
#if defined( DEVEL )
      SUBROUTINE SCHEM( PEDGE )
#else
      SUBROUTINE SCHEM
#endif
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT,     GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,         READ_BPCH2
      USE DAO_MOD,        ONLY : AD, T
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE ERROR_MOD,      ONLY : ALLOC_ERR
      USE TIME_MOD,       ONLY : GET_MONTH,        GET_TAU
      USE TIME_MOD,       ONLY : GET_TS_CHEM,      TIMESTAMP_STRING
      USE TRACER_MOD,     ONLY : N_TRACERS,        STT
      USE TRACER_MOD,     ONLY : TRACER_MW_KG,     XNUMOLAIR
      USE TRACERID_MOD,   ONLY : IDTACET, IDTALD2, IDTALK4, IDTC2H6
      USE TRACERID_MOD,   ONLY : IDTC3H8, IDTCH2O, IDTH2O2, IDTHNO4
      USE TRACERID_MOD,   ONLY : IDTISOP, IDTMACR, IDTMEK,  IDTMP  
      USE TRACERID_MOD,   ONLY : IDTMVK,  IDTPMN,  IDTPRPE, IDTR4N2
      USE TRACERID_MOD,   ONLY : IDTRCHO
      USE TRANSFER_MOD,   ONLY : TRANSFER_ZONAL
      USE TROPOPAUSE_MOD, ONLY : GET_MIN_TPAUSE_LEVEL, ITS_IN_THE_STRAT
      ! jpp, 7/14/09: putting in the bromocarbon tracer id's so
      !               a simple decay term can be applied to avoid buildup.
      USE TRACERID_MOD,   ONLY : IDTCH3BR, IDTCH2Br2, IDTCHBr3
      USE DAO_MOD,        ONLY : SUNCOS
      USE TRACER_MOD,     ONLY : TCVV
      USE TRANSFER_MOD,   ONLY : TRANSFER_3D

      USE CMN_SIZE_MOD         ! Size parameters

#if defined( DEVEL )
      USE GC_TYPE_MOD,    ONLY : GC_MET_LOCAL
      USE GC_TYPE2_MOD,   ONLY : EXT_STRATOH, EXT_SJVALUE
      USE GC_TYPE2_MOD,   ONLY : EXT_COLOSS,  EXT_COPROD
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% KLUDGE FOR NESTED GRIDS (bmy, 7/3/12)
!%%% Use the BR_SCALING variable to determine if we have to shut off
!%%% the bromine species. 
!%%%
      USE BROMOCARB_MOD,  ONLY : BR_SCALING
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IMPLICIT NONE
#     include "define.h"
! 
! !REVISION HISTORY: 
!  01 Oct 1995 - M. Prather  - Initial version
!  20 Nov 1999 - Q. Li - Initial version

!  (1 ) Now read all inputs (stratospheric OH, monthly mean J-values,  
!        P(CO) rates, and L(CO) rates) from binary punch file format. 
!        (bmy, 12/10/99) 
!  (2 ) Uses READ_BPCH2 to read from binary file format (bmy, 12/10/99)
!  (3 ) Make sure the DO-loops go in the order N-L-J-I to avoid disk
!        swapping problems (bmy, 12/10/99)
!  (4 ) Remove reactions for HNO3 photolysis and HNO3 + OH.  The HNO3
!        concentrations that we read in from disk are from Dylan's 2-D
!        model, where chemistry is already taken into account. 
!        (qli, bmy, 12/23/99)
!  (5 ) Remove obsolete code from 12/23/99. (bmy, 4/18/00)
!  (6 ) Bug fixes: Cap RDLOSS so that it does not exceed 1.0.
!        Now declare RDLOSS, T1L, RC, K0, K1, K2, K3, M as REAL*8 
!        Cosmetic changes & update comments (bmy, 5/4/00)
!  (7 ) Reference F90 module "bpch2_mod" which contains routine "read_bpch2"
!        for reading data from binary punch files (bmy, 6/28/00)
!  (8 ) Now all monthly mean J-values are in the same file (bmy, 6/30/00)
!  (9 ) Now use function GET_TAU0 (from "bpch2_mod.f") to return the TAU0 
!        value used to index the binary punch file. (bmy, 7/20/00)
!  (10) Declared arrays for reading data from disk to be both ALLOCATABLE
!        and SAVE.  Also cosmetic changes & some cleanup. (bmy, 9/8/00) 
!  (11) Activated parallel DO-loops (bmy, 12/12/00)
!  (12) Now use 3 arguments (M/D/Y) in call to GET_TAU0.  ARRAY needs to be 
!        of size (IIPAR,JJPAR).  Use JJPAR,LGLOB in calls to READ_BPCH2.
!        Use TRANSFER_ZONAL (from "transfer_mod.f") to cast from REAL*4 to 
!        REAL*8 and resize arrays to (JJPAR,LLPAR).  Updated comments, 
!        made cosmetic changes. (bmy, 9/27/01)
!  (13) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (14) Now read COprod and COloss files directly from the
!        DATA_DIR/pco_lco_200203/ subdirectory.  Also read stratOH files
!        directly from the DATA_DIR/stratOH_200203/ subdirectory.  Also 
!        read stratjv files directly from the DATA_DIR/stratjv_200203/ 
!        subdirectory. (bmy, 4/2/02)
!  (15) Now reference AD and T from "dao_mod.f".  Also reference routine
!        ALLOC_ERR from "error_mod.f".  Now reference IDTOX, IDTNOX, etc.
!        from "tracerid_mod.f". (bmy, 11/6/02)
!  (16) Now use functions GET_TS_CHEM, GET_MONTH and GET_TAU, and 
!        TIMESTAMP_STRING from the new "time_mod.f".   Also call READ_BPCH2 
!        with QUIET=.TRUE., which prevents info from being printed to the 
!        log file. (bmy, 3/14/03)
!  (17) LINUX has a problem putting a function call w/in a WRITE statement.  
!        Now save output from TIMESTAMP_STRING to STAMP and print that.
!        (bmy, 9/29/03)
!  (18) Now reference STT and TRACER_MW_KG from "tracer_mod.f".  Now reference
!        DATA_DIR from "directory_mod.f".  Bug fix: now loop over N_TRACERS
!        and not NNPAR.  NNPAR is the max # of tracers but may not be the
!        actual number of tracers. (bmy, 7/20/04)
!  (19) Now references GET_MIN_TPAUSE_LEVEL and ITS_IN_THE_STRAT from
!        "tropopause_mod.f".  Now remove reference to CMN, it's obsolete.
!        (bmy, 8/22/05)
!  (20) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (21) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  07 Mar 2012 - R. Yantosca - Start DO loops from L=1 instead of L=LMIN
!EOP
!------------------------------------------------------------------------------
!BOC
      REAL*8, INTENT(IN) :: PEDGE(IIPAR, JJPAR, LLPAR)

!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE             :: FIRST = .TRUE.

      INTEGER                   :: I, IOS, J, L, N, NN, LMIN
      INTEGER, SAVE             :: MONTHSAVE = 0 
      
      ! Number of photolysis species (currently is 13)
      INTEGER, PARAMETER        :: NSPHOTO = 13  

      ! Tracers that undergo photolysis loss in the stratosphere
      INTEGER                   :: SPHOTOID(NSPHOTO) = (/ 
     &                               3,  8,  9, 10, 11, 12, 13, 
     &                              14, 17, 20, 22, 23, 24/)
!jp_to_continue     &                              51, 52 /) ! jpp, added CHBr3 & CH2Br2

      ! Character variables
      CHARACTER(LEN=16 )        :: STAMP
      CHARACTER(LEN=255)        :: FILENAME, dayfile, nightfile

      ! REAL*4 arrays -- for reading from binary data files
      REAL*4                    :: ARRAY(1,JJPAR,LGLOB) 
      REAL*4, ALLOCATABLE, SAVE :: STRATOH(:,:)
      REAL*4, ALLOCATABLE, SAVE :: SJVALUE(:,:,:) 
      REAL*4, ALLOCATABLE, SAVE :: COPROD(:,:)
      REAL*4, ALLOCATABLE, SAVE :: COLOSS(:,:)

      ! REAL*8 variables
      REAL*8                    :: k0,     k1,     k2,  k3, XTAU
      REAL*8                    :: DTCHEM, RDLOSS, T1L, M,  TK, RC 

      ! jpp, 6/27/2011: adding an array to store the
      ! inorganic bromine species for stratosphere.
      REAL*4  :: Bry_temp(IGLOB,JGLOB,LGLOB)
      REAL*8  :: Bry_day(IIPAR, JJPAR, LLPAR, 6)
      REAL*8  :: Bry_night(IIPAR, JJPAR, LLPAR, 6)
      INTEGER :: br_nos(6)=(/44, 45, 46, 47, 48, 50/) ! Br2, Br, BrO, HOBr, HBr, BrNO3
      INTEGER :: IJWINDOW, IB
      REAL*8  :: GMU

      ! External functions
      REAL*8, EXTERNAL          :: BOXVL

      !=================================================================
      ! SCHEM begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - SCHEM: Strat chemistry at ', a )

      !=================================================================
      ! If it is the first call to SCHEM, allocate arrays for reading 
      ! data. These arrays are declared SAVE so they will be preserved 
      ! between calls. 
      !=================================================================
      IF ( FIRST ) THEN 
         ALLOCATE( STRATOH( JJPAR, LLPAR ), STAT=IOS )
         IF ( IOS /= 0 ) CALL ALLOC_ERR( 'STRATOH' )
         STRATOH = 0e0

         ALLOCATE( SJVALUE( JJPAR, LLPAR, NSPHOTO ), STAT=IOS )
         IF ( IOS /= 0 ) CALL ALLOC_ERR( 'SJVALUE' )
         SJVALUE = 0e0

         ALLOCATE( COPROD( JJPAR, LLPAR ), STAT=IOS )
         IF ( IOS /= 0 ) CALL ALLOC_ERR( 'COPROD' )
         COPROD = 0e0

         ALLOCATE( COLOSS( JJPAR, LLPAR ), STAT=IOS )
         IF ( IOS /= 0 ) CALL ALLOC_ERR( 'COLOSS' )
         COLOSS = 0e0
      ENDIF

#if defined( EXTERNAL_FORCING )

         STRATOH = EXT_STRATOH
         SJVALUE = EXT_SJVALUE
         COPROD  = EXT_COPROD
         COLOSS  = EXT_COLOSS

#else
      !=================================================================
      ! If it is a new month (or the first call to SCHEM), 
      ! do the following:
      !
      ! (1) Read archived J-values and store in SJVALUE
      ! (2) Read archived CO production rates and store in COPROD
      ! (3) Read archived CO loss rates and store in COLOSS
      !
      ! NOTES
      ! (a) All of the above-mentioned data are stored in binary punch 
      !     files, for ease of use.  
      !
      ! (b) STRATOH, SJVALUE, CO_PROD, and CO_LOSS are now declared 
      !     as both ALLOCATABLE and SAVE.  If SCHEM is called, then 
      !     data will be declared for these arrays, and the values in 
      !     these arrays will be preserved between calls.  
      !
      ! (c) If SCHEM is never called (i.e. if you are running another 
      !     type of chemistry simulation), then memory never gets 
      !     allocated to STRATOH, SJVALUE, CO_PROD, and CO_LOSS.  
      !     This saves on computational resources.       
      !=================================================================
      IF ( GET_MONTH() /= MONTHSAVE .or. FIRST ) THEN
         MONTHSAVE = GET_MONTH()
      
         ! TAU value at the beginning of this month
         XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )

         !==============================================================
         ! Read this month's OH 
         !==============================================================
         FILENAME = TRIM( DATA_DIR ) // 'stratOH_200203/stratOH.' // 
     &              GET_NAME_EXT()   // '.'                       // 
     &              GET_RES_EXT()

         ! Read data
         CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 1,     
     &                    XTAU,      1,         JJPAR,     
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% SCHEM PATCH: Do not call TRANSFER_ZONAL when using the GEOS-5 grid
!%%% since the SCHEM data are saved at the reduced resolution.  This will
!%%% be fixed in a subsequent version. (ltm, bmy, 6/2/10)
!%%%---------------------------------------------------------------------
!%%%         ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
!%%%         CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATOH )
!%%%---------------------------------------------------------------------
#if   defined( GEOS_5 )  || defined( MERRA ) || defined( GEOS_57 )
         ! Cast from REAL*4 to REAL*8
         STRATOH = ARRAY(1,:,:)
#else
         ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
         CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATOH )
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
         !==============================================================
         ! Read in monthly mean archived J-values
         !==============================================================
         FILENAME = TRIM( DATA_DIR ) // 'stratjv_200203/stratjv.' //
     &              GET_NAME_EXT()   // '.'                       // 
     &              GET_RES_EXT()

         DO NN = 1, NSPHOTO
            N = SPHOTOID(NN)

            ! Read data
            CALL READ_BPCH2( FILENAME, 'JV-MAP-$', N,     
     &                       XTAU,      1,         JJPAR,     
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% SCHEM PATCH: Do not call TRANSFER_ZONAL when using the GEOS-5 grid
!%%% since the SCHEM data are saved at the reduced resolution.  This will
!%%% be fixed in a subsequent version. (ltm, bmy, 6/2/10)
!%%%---------------------------------------------------------------------
!%%% Original code here:
!%%%            ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
!%%%            CALL TRANSFER_ZONAL( ARRAY(1,:,:), SJVALUE(:,:,NN) )
!%%%---------------------------------------------------------------------
#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )
            ! Cast from REAL*4 to REAL*8
            SJVALUE(:,:,NN) = ARRAY(1,:,:)
#else
            ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), SJVALUE(:,:,NN) )
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ENDDO

         !==============================================================
         ! Read in CO production rates
         !==============================================================
         FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/COprod.' //
     &              GET_NAME_EXT()   // '.'                      // 
     &              GET_RES_EXT()

         ! Read data
         CALL READ_BPCH2( FILENAME, 'PORL-L=$', 9,     
     &                    XTAU,      1,         JJPAR,     
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% SCHEM PATCH: Do not call TRANSFER_ZONAL when using the GEOS-5 grid
!%%% since the SCHEM data are saved at the reduced resolution.  This will
!%%% be fixed in a subsequent version. (ltm, bmy, 6/2/10)
!%%%---------------------------------------------------------------------
!%%% Original code here:
!%%%         ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
!%%%         CALL TRANSFER_ZONAL( ARRAY(1,:,:), COPROD )
!%%%---------------------------------------------------------------------
#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )
         ! Cast from REAL*4 to REAL*8
         COPROD = ARRAY(1,:,:)
#else
         ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
         CALL TRANSFER_ZONAL( ARRAY(1,:,:), COPROD )
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
         !==============================================================
         ! Read in CO loss rates
         !==============================================================
         FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/COloss.' //
     &              GET_NAME_EXT()   // '.'                      // 
     &              GET_RES_EXT()

         ! Read data
         CALL READ_BPCH2( FILENAME, 'PORL-L=$', 10,    
     &                    XTAU,      1,         JJPAR,     
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% SCHEM PATCH: Do not call TRANSFER_ZONAL when using the GEOS-5 grid
!%%% since the SCHEM data are saved at the reduced resolution.  This will
!%%% be fixed in a subsequent version. (ltm, bmy, 6/2/10)
!%%%---------------------------------------------------------------------
!%%% Original code:
!%%%         ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
!%%%         CALL TRANSFER_ZONAL( ARRAY(1,:,:), COLOSS )
!%%%---------------------------------------------------------------------
#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )
         ! Cast from REAL*4 to REAL*8
         COLOSS = ARRAY(1,:,:)
#else
         ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
         CALL TRANSFER_ZONAL( ARRAY(1,:,:), COLOSS )
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !==============================================================
         ! Read in stored Bry species concentrations for stratosphere.
         ! Stored by Q. Liang using the GEOS CCM. (jpp, 6/27/2011)
         !==============================================================

         ! the daytime concentrations
         dayfile = TRIM( DATA_DIR ) // 'bromine_201205/' //
     &             'CCM_stratosphere_Bry/Bry_Stratosphere_day.bpch.'//
     &             GET_NAME_EXT()   // '.' // GET_RES_EXT()

         ! the nighttime concentrations
         nightfile = TRIM(DATA_DIR) // 'bromine_201205/' //
     &             'CCM_stratosphere_Bry/Bry_Stratosphere_night.bpch.'//
     &             GET_NAME_EXT()   // '.' // GET_RES_EXT()

         DO NN = 1, 6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% KLUDGE FOR NESTED GRIDS (bmy, 7/3/12)
!%%% Right now there is no good way to shut off bromine.  We don't have
!%%% nested grid data for BRY_TEMP and BRY_DAY yet.  Therefore, use the
!%%% kludge that if BR_SCALING = 0 then don't read any of these files.
!%%% Figure out a more consistent way to implement this later.
!%%% (bmy, 7/3/12)
!%%%
            IF ( ABS( Br_SCALING ) > 0d0 ) THEN

               ! 1. Read daytime data
               CALL READ_BPCH2( DAYFILE, 'IJ-AVG-$', br_nos(NN),    
     &                          XTAU,    IGLOB,      JGLOB,     
     &                          LGLOB,   Bry_temp,   QUIET=.TRUE. )

               ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
               CALL TRANSFER_3D( Bry_temp(:,:,:), Bry_day(:,:,:,NN) )

               !==============================================================
               ! now convert units from [ppt] to [kg]
               !  STT(I,J,L,N) [v/v] * AD(I,J,L) / TCVV(N) = STT(I,J,L,N) [kg]
               !==============================================================

               ! 2. Read nighttime data
               CALL READ_BPCH2( NIGHTFILE, 'IJ-AVG-$', br_nos(NN),    
     &                          XTAU,      IGLOB,      JGLOB,     
     &                          LGLOB,     Bry_temp,   QUIET=.TRUE. )

               ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
               CALL TRANSFER_3D( Bry_temp(:,:,:), Bry_night(:,:,:,NN) )

            ELSE
               ! Echo msg
               WRITE( 6, 777 )
               WRITE( 6, 778 )
         
               ! Set BRY_DAY and BRY_NIGHT to low, positive values
               ! so that SMVGEAR won't crash
               Bry_day   = 1e-30
               Bry_night = 1e-30
            ENDIF

 777      FORMAT( '===> KLUDGE: Turn off Bromine in SCHEM       <====' )
 778      FORMAT( '===> Set BRY_DAY and BRY_NIGHT to low values <====' )
!%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !==============================================================
            ! now convert units from [ppt] to [kg]
            !  STT(I,J,L,N) [v/v] * AD(I,J,L) / TCVV(N) = STT(I,J,L,N) [kg]
            !==============================================================
!jpt !$OMP PARALLEL DO
!jpt !$OMP+DEFAULT( SHARED )
!jpt !$OMP+PRIVATE( IJWINDOW, I, J, L )
!jpt !$OMP+SCHEDULE( DYNAMIC )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               ! 1. the daytime values
               bry_day(I,J,L,NN) = bry_day(I,J,L,NN)
     &                           * 1.d-12    ! convert from [ppt]
     &                           * AD(I,J,L) / TCVV(br_nos(NN))

               ! 2. the nighttime values
               bry_night(I,J,L,NN) = bry_night(I,J,L,NN)
     &                             * 1.d-12    ! convert from [ppt]
     &                             * AD(I,J,L) / TCVV(br_nos(NN))

               ! Only proceed for stratospheric boxes
               IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

                  ! Set the Bry boundary conditions. Simulated
                  ! output from the GEOS5 CCM stratosphere.
                  ! (jpp, 6/27/2011)
                  IJWINDOW   = (J-1)*IIPAR + I
!                  GMU        = SUNCOS(IJWINDOW)

                  ! daytime
!                  if ( GMU > 0.d0 ) then
                  IF (SUNCOS(IJWINDOW) > 0.d0) THEN
                     STT(I,J,L, br_nos(NN) ) = Bry_day(I,J,L,NN)
                  ELSE
                     ! nighttime
                     STT(I,J,L, br_nos(NN) ) = Bry_night(I,J,L,NN)
                  ENDIF

               ENDIF

            ENDDO     
            ENDDO
            ENDDO
!jpt !$OMP END PARALLEL DO

         ENDDO
      ENDIF

#endif

      !=================================================================
      ! Do photolysis for selected tracers with this 
      ! month's archived J-values
      !=================================================================

      ! Get the minimum level extent of the ann mean tropopause
      LMIN = GET_MIN_TPAUSE_LEVEL()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, NN )
!$OMP+SCHEDULE( DYNAMIC )
      DO NN = 1, NSPHOTO
         N = SPHOTOID(NN)

         DO L = 1,    LLPAR
         DO J = 1,    JJPAR
         DO I = 1,    IIPAR

            ! Only proceed for stratospheric boxes
            IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

               ! Compute photolysis loss 
               STT(I,J,L,N) = STT(I,J,L,N) * 
     &                        EXP( -SJVALUE(J,L,NN) * DTCHEM )
            ENDIF

         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !print*, 'In schem, done with photolysis'

      !=================================================================
      ! CO is special -- 
      ! use archived P, L rates for CO chemistry in stratosphere
      !=================================================================
      CALL CO_STRAT_PL( COPROD, COLOSS, PEDGE )

      !=================================================================
      ! Reaction with OH -- compute rate constants for each tracer
      !=================================================================
      !print*, 'In schem, before reaction with OH'

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, M, TK, RC, k0, k1, RDLOSS, T1L )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1,    N_TRACERS
      DO L = 1,    LLPAR
      DO J = 1,    JJPAR
      DO I = 1,    IIPAR

         ! Only proceed for stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

            ! Density of air at grid box (I,J,L) in molec/cm3
            M = AD(I,J,L) / BOXVL(I,J,L) * XNUMOLAIR

            ! Temperature at grid box (I,J,L) in K
            TK = T(I,J,L)

            ! Select proper reaction rate w/ OH for the given tracer
            ! Some rates are temperature or density dependent
            IF ( N == IDTALK4 ) THEN
               RC = 8.20D-12 * EXP(  -300.D0 / TK )
         
            ELSE IF ( N == IDTISOP ) THEN
               RC = 2.55D-11 * EXP(   410.D0 / TK )

            ELSE IF ( N == IDTH2O2 ) THEN 
               RC = 2.90D-12 * EXP(  -160.D0 / TK )
               
            ELSE IF ( N == IDTACET ) THEN
               RC = 1.70D-12 * EXP(  -600.D0 / TK )
               
            ELSE IF ( N == IDTMEK  ) THEN 
               RC = 2.92D-13 * EXP(   414.D0 / TK )
            
            ELSE IF ( N == IDTALD2 ) THEN 
               RC = 1.40D-12 * EXP( -1860.D0 / TK )
               
            ELSE IF ( N == IDTRCHO ) THEN 
               RC = 2.00D-11
               
            ELSE IF ( N == IDTMVK  ) THEN 
               RC = 4.13D-12 * EXP(   452.D0 / TK )
                  
            ELSE IF ( N == IDTMACR ) THEN 
               RC = 1.86D-11 * EXP(  -175.D0 / TK )
            
            ELSE IF ( N == IDTPMN  ) THEN 
               RC = 3.60D-12

            ELSE IF ( N == IDTR4N2 ) THEN
               RC = 1.30D-12
               
            ELSE IF ( N == IDTPRPE ) THEN 
               k0 = 8.0D-27 * ( 300.D0 / TK )**3.5
               k1 = 3.0D-11

               RC = k1 * k0 * M / ( k1 + k0*M )
               RC = RC * 0.5 ** (1 / ( 1 + LOG10( k0*M/k1 )**2 ) )

            ELSE IF ( N == IDTC3H8 ) THEN
               RC = 8.00D-12 * EXP(  -590.D0 / TK )
               
            ELSE IF ( N == IDTCH2O ) THEN
               RC = 1.00D-12
               
            ELSE IF ( N == IDTC2H6 ) THEN
               RC =  7.9D-12 * EXP( -1030.D0 / TK )
               
            ELSE IF ( N == IDTHNO4 ) THEN
               RC = 1.30D-12 * EXP(   380.D0 / TK )
               
            ELSE IF ( N == IDTMP ) THEN
               RC = 1.14D-12 * EXP(   200.D0 / TK )

            ! +++++++++++++++++++++++++++++++++++++++++++
            ! jpp, 7/14/09: adding decay by OH only
            !   for stratospheric loss of bromocarbons.
            !   FLAG: CHBr3 should have a simple decay
            !       in addition to OH to approximate hv.
            ! +++++++++++++++++++++++++++++++++++++++++++
            ELSE IF ( N == IDTCH2Br2 ) THEN
               
               RC = 2.0d-12 * EXP ( - 840.d0 / TK )

            ELSE IF ( N == IDTCH3Br ) THEN

               RC = 2.35d-12 * EXP ( - 1300.d0 / TK )

            ELSE IF ( N == IDTCHBr3 ) THEN

               RC = 1.35d-12 * EXP ( - 600.d0 / TK )
            ! +++++++++++++++++++++++++++++++++++++++++++

            ELSE
               RC = 0d0
               
            ENDIF

            ! Compute loss with OH based on the rate constants from above
            ! Cap RDLOSS so that it does not exceed 1.0 (bmy, 5/4/00)
            RDLOSS       = RC * STRATOH(J,L) * DTCHEM
            RDLOSS       = MIN( RDLOSS, 1d0 )

            ! T1L is the absolute amount of STT lost to rxn with OH
            ! Subtract T1L from STT 
            T1L          = STT(I,J,L,N) * RDLOSS
            STT(I,J,L,N) = STT(I,J,L,N) - T1L
         
            ! Oxidation of PRPE as source of ACET with 80% yield
            IF ( N == IDTPRPE ) THEN
               STT(I,J,L,IDTACET) = STT(I,J,L,IDTACET) +
     &              0.8d0 * T1L * 
     &              TRACER_MW_KG(IDTACET) / TRACER_MW_KG(IDTPRPE)
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Set FIRST = .FALSE. -- we have been thru SCHEM at least once now
      FIRST = .FALSE.

      END SUBROUTINE SCHEM
!EOC
