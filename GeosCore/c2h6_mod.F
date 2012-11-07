!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: c2h6_mod
!
! !DESCRIPTION: Module C2H6\_MOD contains variables and routines used for the
!  tagged C2H6 (ethane) simulation. (xyp, qli, bmy, 7/28/01, 4/5/06)
!\\
!\\
! !INTERFACE:
!
      MODULE C2H6_MOD
!
! !USES:
!
      IMPLICIT NONE
#     include "define.h" 
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: GET_C2H6_ANTHRO
      PUBLIC  :: EMISSC2H6
      PUBLIC  :: CHEMC2H6
      PUBLIC  :: CLEANUP_C2H6
!
! PRIVATE DATA MEMBERS:
!
      PRIVATE :: NGASC2H6
      PRIVATE :: FMOL_C2H6
      PRIVATE :: XNUMOL_C2H6
!
! !REMARKS:
!  Setting LSPLIT = T in "input.geos" will run with the following tracers:
!     (1) Total C2H6
!     (2) C2H6 from biomass burning
!     (3) C2H6 from biofuel burning
!     (4) C2H6 from natural gas leaking/venting (e.g. "anthro" C2H6)
!                                                                             .
!  Setting LSPLIT = F in "input.geos" will run w/ the following tracers:
!     (1) Total C2H6
!
! !REVISION HISTORY:
!  (1 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (2 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (3 ) Now reference BXHEIGHT and T from "dao_mod.f".  Also references
!        "error_mod.f".  Removed obsolete code.  Now references F90 module
!         tracerid_mod.f". (bmy, 11/15/02)
!  (4 ) Now references "grid_mod.f" and the new "time_mod.f" (bmy, 2/11/03)
!  (5 ) Now references "directory_mod.f", "logical_mod.f", and "tracer_mod.f".
!        (bmy, 7/20/04)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) Now modified
!  22 Mar 2012 - M. Payer    - Add function GET_C2H6_ANTHRO
!  22 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Molecular weight of C2H6 [kg/mole]
      REAL*8, PARAMETER   :: FMOL_C2H6 = 30d-3

      ! Ratio of [molec C2H6/kg C2H6]
      REAL*8, PARAMETER   :: XNUMOL_C2H6 = 6.022d+23/FMOL_C2H6
!
! PRIVATE TYPES:
!
      ! Global monthly mean natural gas C2H6 emissions
      REAL*8, ALLOCATABLE :: NGASC2H6(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissc2h6
!
! !DESCRIPTION: Subroutine EMISSC2H6 reads in C2H6 emissions for the Tagged
!  C2H6 run. (xyp, qli, bmy, 7/21/00, 4/5/06)
!\\
!\\
! !INTERFACE:
!
#if defined( DEVEL )      
      SUBROUTINE EMISSC2H6( State_Met, EMISS, BIO_EMISS)
#else
      SUBROUTINE EMISSC2H6
#endif
!
! !USES:
!
      USE BIOMASS_MOD,   ONLY : BIOMASS
      USE TRACERID_MOD,  ONLY : IDBC2H6      
      USE BIOFUEL_MOD,   ONLY : BIOFUEL, BIOFUEL_BURN
      USE DIAG_MOD,      ONLY : AD36
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE GEIA_MOD,      ONLY : READ_C3H8_C2H6_NGAS, TOTAL_FOSSIL_TG
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,   ONLY : LSPLIT, LBIOMASS, LBIOFUEL, LANTHRO
      USE TIME_MOD,      ONLY : GET_MONTH, GET_TS_EMIS
      USE TRACER_MOD,    ONLY : STT
      USE TRACER_MOD,    ONLY : ITS_A_C2H6_SIM
      USE TRACERID_MOD,  ONLY : IDBFC2H6, IDEC2H6, IDTC2H6
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
#if defined( DEVEL )
      USE GIGC_State_Met_Mod, ONLY : MetState
#endif

      USE CMN_SIZE_MOD        ! Size parameters
      USE CMN_O3_MOD          ! EMISTC2H6
      USE CMN_DIAG_MOD        ! Diagnostic arrays & switches
      USE CMN_MOD             ! STT, etc.
!
! !INPUT PARAMETERS:
!
#if defined( DEVEL )
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
#if defined( DEVEL )
      REAL*8, INTENT(INOUT)  :: EMISS(IIPAR,JJPAR)
      REAL*8, INTENT(INOUT)  :: BIO_EMISS(:,:,:)
#endif
!
! !REVISION HISTORY:
!  (1 ) BURNEMIS and BIOFUEL are now dimensioned with IIPAR,JJPAR instead of
!        IIPAR,JJPAR.  Remove BXHEIGHT from the arg list, since ND28 and ND36
!        diags are archived in BIOBURN and BIOFUEL_BURN.  Now use routine
!        TRANSFER_2D from "transfer_mod.f" to cast from REAL*4 to REAL*8.
!        Now print emission totals for C2H6 emissions to stdout. (bmy, 1/25/02)
!  (2 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (3 ) Now references IDBC2H6 etc from "tracerid_mod.f".  Now make FIRSTEMISS
!        a local SAVEd variable instead of an argument. (bmy, 11/15/02)
!  (4 ) Now use GET_AREA_CM2 from "grid_mod.f" to get grid box surface
!        area in cm2.  Remove references to DXYP.  Use routines GET_MONTH
!        and GET_TS_EMIS from "time_mod.f".  Remove MONTH from call to
!        BIOBURN. (bmy, 2/11/03)
!  (5 ) Now replace CMN_SETUP w/ references from "logical_mod.f" and
!        "directory_mod.f".  Now references STT from "tracer_mod.f".
!        Replace LFOSSIL with LANTHRO (bmy, 7/20/04)
!  (6 ) Now make sure all USE statements are USE, ONLY.  Also eliminate 
!        reference to BPCH2_MOD, it's obsolete. (bmy, 10/3/05)
!  (7 ) Now modified for new "biomass_mod.f" (bmy, 4/5/06)
!  (8 ) BIOMASS(:,:,IDBCO) from "biomass_mod.f" is now in units of 
!        [atoms C/cm2/s].  Adjust unit conversion accordingly. (bmy, 9/27/06)
!  (9 ) Now IDBC2H6 is defined in TRACERID_MOD. (fp, hotp , 7/31/09)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  22 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE.
      INTEGER, SAVE          :: LASTMONTH  = -99
      INTEGER                :: I, J, L, AS
      REAL*4                 :: ARRAY(IIPAR,JJPAR)
      REAL*8                 :: AREA_CM2,    XTAU  
      REAL*8                 :: E_C2H6_BB,   E_C2H6_BF
      REAL*8                 :: E_C2H6_NGAS, DTSRCE
      CHARACTER(LEN=255)     :: FILENAME
!
! !EXTERNAL FUNCTIONS:
!
      REAL*8, EXTERNAL       :: BOXVL

      !=================================================================
      ! EMISS_C2H6 begins here!
      !=================================================================
      IF ( FIRSTEMISS ) THEN 

         ! Allocate NGASC2H6 array, if this is the first emission 
         CALL INIT_C2H6

         ! Set first-time flag to false
         FIRSTEMISS = .FALSE.
      ENDIF

      ! DTSRCE is the number of seconds per emission timestep
      DTSRCE = GET_TS_EMIS() * 60d0

      !=================================================================
      ! Process biomass C2H6 emissions ored in BURNEMIS(IDBC2H6,:,:) 
      ! in [molec C/cm3/s].  Convert to [kg C2H6] and store in STT.
      !=================================================================
      IF ( LBIOMASS ) THEN

         ! Only process biomass C2H6 emissions if offline C2H6 sim.
         ! For fullchem sim, C2H6 biomass emissions are read from GFED2/3. 
         ! (mpayer, 3/22/12)
         IF ( ITS_A_C2H6_SIM() ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, E_C2H6_BB )
            ! Loop over lats & lons
            DO J = 1, JJPAR
            DO I = 1, IIPAR

                  ! Grid box area [cm2]
                  AREA_CM2 = GET_AREA_CM2( I, J, 1 )
            
                  ! Convert [atoms C/cm2/s] to [kg C2H6] and store in E_C2H6
                  E_C2H6_BB = BIOMASS(I,J,IDBC2H6)   / 2.0d0  / 
     &                        XNUMOL_C2H6 * AREA_CM2 * DTSRCE  

                  ! Add BB C2H6 to tracer #1 -- total C2H6 [kg C2H6]
                  STT(I,J,1,1) = STT(I,J,1,1) + E_C2H6_BB  

                  ! Add BB C2H6 to tracer #2 -- BB C2H6
                  IF ( LSPLIT ) THEN
                     STT(I,J,1,2) = STT(I,J,1,2) + E_C2H6_BB 
                  ENDIF

            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
      ENDIF
      
      !=================================================================
      ! Process biofuel C2H6 emissions stored in BIOFUEL(IDBFC2H6,:,:) 
      ! in [molec C/cm3/s.  Convert to [kg C2H6] and store in STT. 
      !=================================================================
      IF ( LBIOFUEL ) THEN

         ! Only read biofuel burning emissions if offline C2H6 sim in order
         ! to avoid double counting of biofuel emissions in fullchem sim
         ! (mpayer, 3/22/12)
         IF ( ITS_A_C2H6_SIM() ) THEN

            ! Read biofuel burning emissions (and update ND34 diagnostic)
#if defined( DEVEL )
            CALL BIOFUEL_BURN( State_Met, BIO_EMISS )
#else
            CALL BIOFUEL_BURN
#endif

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, E_C2H6_BF )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
 
               ! Convert [molec C/cm3/s] to [kg C2H6] and store in E_C2H6
#if defined( DEVEL )
               E_C2H6_BF = BIOFUEL(IDBFC2H6,I,J) / 2.0d0 / 
     &                     XNUMOL_C2H6 * BOXVL(I,J,1,State_Met) * DTSRCE
#else
               E_C2H6_BF = BIOFUEL(IDBFC2H6,I,J) / 2.0d0 / 
     &                     XNUMOL_C2H6 * BOXVL(I,J,1) * DTSRCE  
#endif

               ! Add BF C2H6 to tracer #1 -- total C2H6 [kg C2H6]
               STT(I,J,1,1) = STT(I,J,1,1) + E_C2H6_BF

               ! Add BF C2H6 to tracer #3 -- BF C2H6
               IF ( LSPLIT ) THEN
                  STT(I,J,1,3) = STT(I,J,1,3) + E_C2H6_BF 
               ENDIF

            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
      ENDIF

      !=================================================================
      ! Process anthro (natural gas venting/leakage) C2H6 emissions
      ! This source is 6.3 Tg C/yr, following Wang et al. [1998].
      ! The distribution follows natural gas venting/leakage of CH4.
      ! Contact: Yaping Xiao (xyp@io.harvard.edu)
      !=================================================================
      IF ( LANTHRO ) THEN 

         ! Read C2H6 emissions only if it's a new month
         IF ( GET_MONTH() /= LASTMONTH ) THEN

            ! Fancy output...
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            WRITE( 6, '(a)' ) 'EMISSC2H6: Reading anthro C2H6!'

            ! Read C2H6 emissions [atoms C/cm2/s]
            CALL READ_C3H8_C2H6_NGAS( E_C2H6=ARRAY )
            
            ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY, NGASC2H6 )

            ! Print emission totals in Tg C
            CALL TOTAL_FOSSIL_TG( NGASC2H6, IIPAR, JJPAR, 
     &                            1,        12d-3, 'C2H6' )

            ! Fancy output...
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )

            ! Save current month in LASTMONTH
            LASTMONTH = GET_MONTH()
         ENDIF

         ! Only add anthro emissions to STT if offline C2H6 sim. For
         ! fullchem sim, this is done in emfossil (mpayer, 3/22/12)
         IF ( ITS_A_C2H6_SIM() ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, E_C2H6_NGAS  )
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid box surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, 1 )

               ! Convert NGAS C2H6 from [atoms C/cm2/s] to [kg C2H6]
               E_C2H6_NGAS = NGASC2H6(I,J) / 2.0d0  /
     &                       XNUMOL_C2H6 * AREA_CM2 * DTSRCE 

               ! Add NGAS C2H6 to tracer #1 -- total C2H6 [kg C2H6]
               STT(I,J,1,1) = STT(I,J,1,1) + E_C2H6_NGAS

               ! Add NGAS C2H6 to tracer #4 -- NGAS C2H6
               IF ( LSPLIT ) THEN
                  STT(I,J,1,4) = STT(I,J,1,4) + E_C2H6_NGAS 
               ENDIF

               ! ND36 = Anthro source diagnostic...store as [moleC/cm2]
               ! and convert to [moleC/cm2/s] in DIAG3.F
               IF ( ND36 > 0 ) THEN
                  AD36(I,J,IDEC2H6) = AD36(I,J,IDEC2H6) + 
     &                              ( NGASC2H6(I,J) * DTSRCE )
               ENDIF
  
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
      ENDIF
         
#if defined( DEVEL )
      EMISS(:,:) = EMISS(:,:) + NGASC2H6(:,:)
#endif

      END SUBROUTINE EMISSC2H6
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemc2h6
!
! !DESCRIPTION: Subroutine CHEMC2H6 performs C2H6 chemistry. Loss of C2H6 is
!  via reaction with OH. (xyp, qli, bmy, 10/19/99, 7/20/04)
!\\
!\\
! !INTERFACE:
!
#if defined( DEVEL )
      SUBROUTINE CHEMC2H6( State_Met )
#else
      SUBROUTINE CHEMC2H6
#endif
!
! !USES:
!
      USE DAO_MOD,       ONLY : T
      USE GLOBAL_OH_MOD, ONLY : OH,        GET_GLOBAL_OH
      USE LOGICAL_MOD,   ONLY : LSPLIT 
      USE TIME_MOD,      ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACER_MOD,    ONLY : N_TRACERS, STT
#if defined( DEVEL )
      USE GIGC_State_Met_Mod, ONLY : MetState
#endif

      USE CMN_SIZE_MOD        ! Size parameters
!
! !INPUT PARAMETERS:
!
#if defined( DEVEL )
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
#endif
! 
! !REVISION HISTORY:
!  (1 ) Now do chemistry all the way to the model top. 
!  (2 ) Use monthly mean OH fields for oxidation -- reference the monthly 
!        mean OH array and the routine which reads it from disk in 
!       "global_oh_mod.f" (bmy, 1/25/02)
!  (3 ) Now reference T from "dao_mod.f".  Also make FIRSTCHEM a local SAVEd
!        variable. (bmy, 11/15/02)
!  (4 ) Now use functions GET_MONTH and GET_TS_CHEM from "time_mod.f".
!  (5 ) Now reference STT & N_TRACERS from "tracer_mod.f".  Now reference 
!        LSPLIT from "logical_mod.f" (bmy, 7/20/04)
!  22 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE          :: LASTMONTH = -99
      INTEGER                :: I, J, L, N
      REAL*8                 :: DTCHEM, KRATE

      !=================================================================
      ! CHEMC2H6 begins here! 
      !=================================================================
      IF ( FIRSTCHEM ) THEN
         FIRSTCHEM = .FALSE.  ! save for future use?
      ENDIF

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      !=================================================================
      ! Read in the tropospheric OH fields of the (LMN)th month
      ! OH data will be saved into the OH array of "global_oh_mod.f" 
      !=================================================================
      IF ( GET_MONTH() /= LASTMONTH ) THEN
         CALL GET_GLOBAL_OH( GET_MONTH() )
         LASTMONTH = GET_MONTH()
      ENDIF

      !=================================================================
      ! Do C2H6 chemistry -- C2H6 Loss due to chemical reaction with OH
      !
      ! DECAY RATE: The decay rate (KRATE) is calculated by:
      !
      !    OH + C2H6 -> H2O + C2H5 (JPL '97)
      !    k = 8.7D-12 * exp(-1070/T)
      !
      ! KRATE has units of [ molec^2 C2H6 / cm6 / s ]^-1.
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, N, KRATE )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Loss rate [molec2 C2H6/cm6/s]^-1
#if defined( DEVEL )
         KRATE = 8.7d-12 * EXP( -1070.0d0 / State_Met%T(I,J,L) )
#else
         KRATE = 8.7d-12 * EXP( -1070.0d0 / T(I,J,L) )
#endif

         ! Apply loss to total C2H6 (tracer #1)
         STT(I,J,L,1) = STT(I,J,L,1) *
     &                  ( 1d0 - KRATE * OH(I,J,L) * DTCHEM )

         ! If we are running w/ tagged tracers,
         ! then also apply the loss to each of these
         IF ( LSPLIT ) THEN 
            DO N = 2, N_TRACERS
            
               ! Subtract loss of C2H6 by OH and store in STT [kg C2H6]
               ! Loss = k * [C2H6] * [OH] * dt
               STT(I,J,L,N) = STT(I,J,L,N) *
     &                        ( 1d0 - KRATE * OH(I,J,L) * DTCHEM )
            ENDDO
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE CHEMC2H6
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_c2h6_anthro
!
! !DESCRIPTION: Function GET\_C2H6\_ANTHRO returns the monthly average
!  anthropogenic C2H6 emissions at GEOS-Chem grid box (I,J).  Data will
!  be returned in units of [atoms C/cm2/s].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_C2H6_ANTHRO( I, J, N ) RESULT( C2H6_ANTHRO )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTC2H6
      USE CMN_SIZE_MOD       ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I   ! GEOS-Chem longitude index
      INTEGER, INTENT(IN)   :: J   ! GEOS-Chem latitude index
      INTEGER, INTENT(IN)   :: N   ! GEOS-Chem tracer index
!
! !RETURN VALUE:
!
      REAL*8                :: C2H6_ANTHRO
! 
! !REVISION HISTORY: 
!  22 Mar 2012 - M. Payer    - Initial version adapted from GET_RETRO_ANTHRO
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      !  GET_C2H6_ANTHRO begins here
      !=================================================================

      IF ( N == IDTC2H6 ) THEN
         C2H6_ANTHRO = NGASC2H6(I,J)
      ENDIF

      END FUNCTION GET_C2H6_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_c2h6
!
! !DESCRIPTION: Subroutine INIT\_C2H6 allocates and zeroes the NGASC2H6 array,
!  which holds global monthly mean natural gas C2H6 emissions. (qli, bmy,
!  1/1/01, 10/15/02)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_C2H6
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD
! 
! !REVISION HISTORY:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  22 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      !  INIT_C2H6 begins here!
      !=================================================================

      ! Allocate NGASC2H6 array
      ALLOCATE( NGASC2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NGASC2H6' )

      ! Zero NGASC2H6 array
      NGASC2H6 = 0d0

      END SUBROUTINE INIT_C2H6
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_c2h6
!
! !DESCRIPTION: Subroutine CLEANUP\_C2H6 deallocates the natural gas C2H6
!  emission array.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_C2H6
! 
! !REVISION HISTORY:
!  22 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_C2H6 begins here!
      !=================================================================

      IF ( ALLOCATED( NGASC2H6 ) ) DEALLOCATE( NGASC2H6 )

      END SUBROUTINE CLEANUP_C2H6
!EOC
      END MODULE C2H6_MOD
