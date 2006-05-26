! $Id: c2h6_mod.f,v 1.6 2006/05/26 17:45:14 bmy Exp $
      MODULE C2H6_MOD
!
!******************************************************************************
!  Module C2H6_MOD contains variables and routines used for the tagged 
!  C2H6 (ethane) simulation. (xyp, qli, bmy, 7/28/01, 4/5/06)
!
!  Setting LSPLIT = T in "input.geos" will run with the following tracers:
!     (1) Total C2H6
!     (2) C2H6 from biomass burning
!     (3) C2H6 from biofuel burning
!     (4) C2H6 from natural gas leaking/venting (e.g. "anthro" C2H6)
!
!  Setting LSPLIT = F in "input.geos" will run w/ the following tracers:
!     (1) Total C2H6
!
!  Module Variables:
!  ============================================================================
!  (1 ) NGASC2H6        : Array to store C2H6 emissions from natural gas
!  (2 ) FMOL_C2H6       : Molecular weight of C2H6 [kg/mole]
!  (3 ) XNUMOL_C2H6     : Ratio of molecules C2H6 per kg C2H6
!
!  Module Procedures:
!  ============================================================================
!  (1 ) EMISSC2H6       : Routine that performs emission of C2H6 tracers 
!  (2 ) CHEMC2H6        : Routine that performs chemistry for C2H6 tracers
!  (3 ) INIT_NGAS
! 
!  GEOS-Chem modules referenced by "c2h6_mod.f"
!  ============================================================================
!  (1 ) biofuel_mod.f   : Module containing routines to read biofuel emissions
!  (2 ) biomass_mod.f   : Module containing routines to read biomass emissions
!  (3 ) dao_mod.f       : Module containing arrays for DAO met fields!
!  (4 ) diag_mod.f      : Module containing GEOS-CHEM diagnostic arrays
!  (5 ) error_mod.f     : Module containing NaN and other error check routines
!  (6 ) geia_mod.f      : Module containing routines to read anthro emissions
!  (7 ) grid_mod.f      : Module containing horizontal grid information
!  (8 ) global_oh_mod.f : Module containing routines to read 3-D OH field
!  (9 ) time_mod.f      : Module containing routines to compute date & time
!  (10) tracerid_mod.f  : Module containing pointers to tracers and emissions
!  (11) transfer_mod.f  : Module containing routines to cast and resize arrays
!
!  NOTES:
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
!******************************************************************************
!
      IMPLICIT NONE 

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "c2h6_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE NGASC2H6
      PRIVATE FMOL_C2H6
      PRIVATE XNUMOL_C2H6

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Array to store global monthly mean natural gas C2H6 emissions
      REAL*8, ALLOCATABLE :: NGASC2H6(:,:)

      ! FMOL_C2H6: molecular weight of C2H6 [kg/mole]
      REAL*8, PARAMETER   :: FMOL_C2H6 = 30d-3

      ! XNUMOL_C2H6  ratio of [molec C2H6/kg C2H6]
      REAL*8, PARAMETER   :: XNUMOL_C2H6 = 6.022d+23/FMOL_C2H6

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------
      
      SUBROUTINE EMISSC2H6
!
!******************************************************************************
!  Subroutine EMISSC2H6 reads in C2H6 emissions for the Tagged C2H6 run.
!  (xyp, qli, bmy, 7/21/00, 4/5/06)
!
!  NOTES:
!  (1 ) BURNEMIS and BIOFUEL are now dimensioned with IIPAR,JJPAR instead of
!        IGLOB,JGLOB.  Remove BXHEIGHT from the arg list, since ND28 and ND36
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
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,   ONLY : BIOMASS, IDBC2H6
      USE BIOFUEL_MOD,   ONLY : BIOFUEL, BIOFUEL_BURN
      USE DIAG_MOD,      ONLY : AD36
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE GEIA_MOD,      ONLY : READ_C3H8_C2H6_NGAS, TOTAL_FOSSIL_TG
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,   ONLY : LSPLIT, LBIOMASS, LBIOFUEL, LANTHRO
      USE TIME_MOD,      ONLY : GET_MONTH, GET_TS_EMIS
      USE TRACER_MOD,    ONLY : STT
      USE TRACERID_MOD,  ONLY : IDBC2H6, IDBFC2H6, IDEC2H6
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! STT, etc.
#     include "CMN_O3"       ! EMISTC2H6
#     include "CMN_DIAG"     ! Diagnostic arrays & switches

      ! Local variables
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE.
      INTEGER, SAVE          :: LASTMONTH  = -99
      INTEGER                :: I, J, L, AS
      REAL*4                 :: ARRAY(IGLOB,JGLOB)
      REAL*8                 :: AREA_CM2,    XTAU  
      REAL*8                 :: E_C2H6_BB,   E_C2H6_BF
      REAL*8                 :: E_C2H6_NGAS, DTSRCE
      CHARACTER(LEN=255)     :: FILENAME
      
      ! External functions
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

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, E_C2H6_BB )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
 
            ! Convert [molec C/cm3/s] to [kg C2H6] and store in E_C2H6
            E_C2H6_BB = BIOMASS(I,J,IDBC2H6)       / 2.0d0  / 
     &                  XNUMOL_C2H6 * BOXVL(I,J,1) * DTSRCE  

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
      
      !=================================================================
      ! Process biofuel C2H6 emissions stored in BIOFUEL(IDBFC2H6,:,:) 
      ! in [molec C/cm3/s.  Convert to [kg C2H6] and store in STT. 
      !=================================================================
      IF ( LBIOFUEL ) THEN

         ! Read biofuel burning emissions (and update ND34 diagnostic)
         CALL BIOFUEL_BURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, E_C2H6_BF )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
 
            ! Convert [molec C/cm3/s] to [kg C2H6] and store in E_C2H6
            E_C2H6_BF = BIOFUEL(IDBFC2H6,I,J) / 2.0d0 / 
     &                  XNUMOL_C2H6 * BOXVL(I,J,1) * DTSRCE  

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
            CALL TOTAL_FOSSIL_TG( NGASC2H6, IGLOB, JGLOB, 
     &                            1,        12d-3, 'C2H6' )

            ! Fancy output...
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )

            ! Save current month in LASTMONTH
            LASTMONTH = GET_MONTH()
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, E_C2H6_NGAS  )
         DO J = 1, JJPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            DO I = 1, IIPAR

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
     &                                ( NGASC2H6(I,J) * DTSRCE )
               ENDIF  
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
         
      ! Return to calling program
      END SUBROUTINE EMISSC2H6

!------------------------------------------------------------------------------

      SUBROUTINE CHEMC2H6
!
!******************************************************************************
!  Subroutine CHEM_C2H6 performs C2H6 chemistry. Loss of C2H6 is via reaction 
!  with OH. (xyp, qli, bmy, 10/19/99, 7/20/04)
!
!  Arguments as input:
!  ==========================================================================
!  (1 ) FIRSTCHEM (LOGICAL) : First time flag for chemistry 
!
!  NOTES:
!  (1 ) Now do chemistry all the way to the model top. 
!  (2 ) Use monthly mean OH fields for oxidation -- reference the monthly 
!        mean OH array and the routine which reads it from disk in 
!       "global_oh_mod.f" (bmy, 1/25/02)
!  (3 ) Now reference T from "dao_mod.f".  Also make FIRSTCHEM a local SAVEd
!        variable. (bmy, 11/15/02)
!  (4 ) Now use functions GET_MONTH and GET_TS_CHEM from "time_mod.f".
!  (5 ) Now reference STT & N_TRACERS from "tracer_mod.f".  Now reference 
!        LSPLIT from "logical_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : AIRVOL,    T
      USE GLOBAL_OH_MOD, ONLY : OH,        GET_GLOBAL_OH
      USE LOGICAL_MOD,   ONLY : LSPLIT 
      USE TIME_MOD,      ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACER_MOD,    ONLY : N_TRACERS, STT

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE          :: LASTMONTH = -99
      INTEGER                :: I, J, L, N
      REAL*8                 :: DTCHEM, KRATE

      ! External functions
      REAL*8, EXTERNAL       :: BOXVL

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
         KRATE = 8.7d-12 * EXP( -1070.0d0 / T(I,J,L) )

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

      ! Return to calling program
      END SUBROUTINE CHEMC2H6

!------------------------------------------------------------------------------

      SUBROUTINE INIT_C2H6
!
!******************************************************************************
!  Subroutine INIT_C2H6 allocates and zeroes the NGASC2H6 array, which holds 
!  global monthly mean natural gas C2H6 emissions. (qli, bmy, 1/1/01, 10/15/02)
!
!  NOTES:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      ! Local variables
      INTEGER :: AS

      ! Allocate NGASC2H6 array
      ALLOCATE( NGASC2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NGASC2H6' )

      ! Zero NGASC2H6 array
      NGASC2H6 = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_C2H6

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_C2H6
!
!******************************************************************************
!  Subroutine CLEANUP_C2H6 deallocates the natural gas C2H6 emission array.
!
!  NOTES:
!******************************************************************************
!
      IF ( ALLOCATED( NGASC2H6 ) ) DEALLOCATE( NGASC2H6 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_C2H6

!------------------------------------------------------------------------------

      END MODULE C2H6_MOD
