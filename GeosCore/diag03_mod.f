! $Id: diag03_mod.f,v 1.1 2009/09/16 14:06:36 bmy Exp $
      MODULE DIAG03_MOD
!
!******************************************************************************
!  Module DIAG03_MOD contains arrays and routines for archiving the ND03
!  diagnostic -- Hg emissions, mass, and production. (bmy, 1/21/05, 9/5/06) 
!
!  Module Variables:
!  ============================================================================
!  (1 ) AD03         (REAL*4) : Array for Hg emissions & ocean masses 
!  (2 ) AD03_Hg2_Hg0 (REAL*4) : Array for Hg(II) produced from Hg(0)
!  (3 ) AD03_Hg2_OH  (REAL*4) : Array for Hg(II) produced from OH
!  (4 ) AD03_Hg2_O3  (REAL*4) : Array for Hg(II) produced from O3
!
!  Module Routines:
!  ============================================================================
!  (1 ) ZERO_DIAG03           : Sets all module arrays to zero
!  (2 ) WRITE_DIAG03          : Writes data in module arrays to bpch file
!  (3 ) INIT_DIAG03           : Allocates all module arrays
!  (4 ) CLEANUP_DIAG03        : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag03_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f           : Module w/ routines for binary pch file I/O
!  (2 ) error_mod.f           : Module w/ NaN and other error check routines
!  (3 ) file_mod.f            : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f            : Module w/ horizontal grid information
!  (5 ) time_mod.f            : Module w/ routines to compute date & time
!
!  Nomenclature: 
!  ============================================================================
!  (1 ) Hg(0)  a.k.a. Hg0     : Elemental   mercury
!  (2 ) Hg(II) a.k.a. Hg2     : Divalent    mercury
!  (3 ) HgP                   : Particulate mercury
!
!  NOTES:
!  (1 ) Updated for GCAP grid (bmy, 6/28/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Add 2 extra diagnostics to ND03. Set PD03=15.  (cdh, bmy, 12/15/05)
!  (4 ) Add loss of Hg2 by sea salt (eck, bmy, 4/6/06)
!  (5 ) Replace TINY(1d0) w/ 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (6 ) Updates to mercury simulation (ccc, 5/17/10)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag03_mod.f"
      !=================================================================

      ! Make everything PUBLIC
      PUBLIC 

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      INTEGER              :: ND03, LD03
      INTEGER, PARAMETER   :: PD03 = 25

      ! Arrays
      REAL*4,  ALLOCATABLE :: AD03(:,:,:)
      REAL*4,  ALLOCATABLE :: AD03_Hg2_Hg0(:,:,:)
      REAL*4,  ALLOCATABLE :: AD03_Hg2_Br(:,:,:) !cdh added diagnostic
      REAL*4,  ALLOCATABLE :: AD03_Hg2_OH(:,:,:)
      REAL*4,  ALLOCATABLE :: AD03_Hg2_O3(:,:,:)
      REAL*4,  ALLOCATABLE :: AD03_Hg2_SS(:,:,:)
      REAL*4,  ALLOCATABLE :: AD03_nat(:,:,:)
      REAL*4,  ALLOCATABLE :: AD03_Hg2_SSR(:,:) !CDH for sea salt loss rate
      REAL*4,  ALLOCATABLE :: AD03_Br(:,:,:,:) !CDH for bromine

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE ZERO_DIAG03
!
!******************************************************************************
!  Subroutine ZERO_DIAG03 zeroes the ND03 diagnostic arrays. 
!  (bmy, 1/21/05, 4/6/06)
!
!  NOTES:
!  (1 ) Now references N_Hg_CATS from "tracerid_mod.f".  Now zero AD03_Hg2_SS
!        array. (bmy, 4/6/06)
!  (2 ) Now use broadcast assignment and double precision 0D0 to zero arrays,
!        rather than nested DO loops and single precision 0E0. (cdh, 8/14/08)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : N_Hg_CATS

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER :: I, J, L, N

      !=================================================================
      ! ZERO_DIAG03 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND03 == 0 ) RETURN

      ! Zero arrays
      AD03         = 0D0
      AD03_Hg2_Hg0 = 0D0
      AD03_Hg2_Br  = 0D0 !cdh added diagnostic
      AD03_Hg2_OH  = 0D0
      AD03_Hg2_O3  = 0D0
      AD03_Hg2_SS  = 0D0
      AD03_Hg2_SSR = 0D0 !cdh added diagnostic
      AD03_nat     = 0D0 !cdh moved here from mercury_mod
      AD03_Br      = 0D0 !cdh for bromine

      ! Return to calling program
      END SUBROUTINE ZERO_DIAG03

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG03
!
!******************************************************************************
!  Subroutine WRITE_DIAG03 writes the ND03 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 1/21/05, 2/24/06)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) HG-SRCE  : Anthropogenic HG0 emission      : kg       : 1
!  (2 ) HG-SRCE  : Total mass of oceanic Hg0       : kg       : 1
!  (3 ) HG-SRCE  : Oceanic HgO emission            : kg       : 1
!  (4 ) HG-SRCE  : Land reemission                 : kg       : 1
!  (5 ) HG-SRCE  : Land natural emission           : kg       : 1
!  (6 ) HG-SRCE  : Anthropogenic Hg2 emission      : kg       : 1
!  (7 ) HG-SRCE  : Total mass of oceanic Hg2       : kg       : 1
!  (8 ) HG-SRCE  : Mass of Hg2 sunk in the ocean   : kg       : 1
!  (9 ) HG-SRCE  : Anthropogenic HgP emission      : kg       : 1
!  (10) HG-SRCE  : Henry's law piston velocity Kw  : cm/h     : em timesteps
!  (11) HG-SRCE  : Mass of Hg(C)                   : kg       : 1
!  (12) HG-SRCE  : Converted to Colloidal          : kg       : 1
!  (13) HG-SRCE  : Biomass burning emissions       : kg       : 1
!  (14) HG-SRCE  : Emissions from vegetation       : kg       : 1
!  (15) HG-SRCE  : Emissions from soils            : kg       : 1
!  (16) HG-SRCE  : Flux-up Hg0 volat from ocean    : kg       : 1
!  (17) HG-SRCE  : Flux-down Hg0 dry dep to ocean  : kg       : 1
!  (18) PL-HG2-$ : Production of Hg2 from Hg0      : kg       : 1
!  (19) PL-HG2-$ : Production of Hg2 from rxn w/OH : kg       : 1
!  (20) PL-HG2-$ : Production of Hg2 from rxn w/O3 : kg       : 1
!  (21) PL-HG2-$ : Loss of Hg2 from rxn w/ seasalt : kg       : 1 
!
!  NOTES:
!  (1 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Add HgC ocean mass and converted to colloidal to ND03 diagnostic.
!        The units of the Kw and conversion terms in ND03 should be kg
!        and not divided by the scale factor. (cdh, sas, bmy, 2/26/02)
!  (4 ) Replace TINY(1d0) w/ 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (5 ) Fixed tracer numbers (NN) for 'PL-HG2-$' diagnostic quantities.
!        (cdh, 8/13/08)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,     ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe
      USE TIME_MOD,     ONLY : GET_CT_CHEM ! CDH for sea salt loss rate
      USE TRACERID_MOD, ONLY : N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! TINDEX

      ! Local variables
      INTEGER               :: CENTER180, HALFPOLAR,   IFIRST
      INTEGER               :: JFIRST,    LFIRST,      LMAX
      INTEGER               :: M,         N,           NN
      REAL*4                :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*4                :: LONRES,    LATRES
      REAL*8                :: DIAGb,     DIAGe,       SCALE
      CHARACTER(LEN=20)     :: MODELNAME 
      CHARACTER(LEN=40)     :: CATEGORY,  RESERVED,    UNIT
      REAL*8                :: NCHEMSTEP !CDH for sea salt loss rate

      !=================================================================
      ! WRITE_DIAG03 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND03 == 0 ) RETURN

      ! Initialize
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = GET_HALFPOLAR()
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LATRES    = DJSIZE
      LFIRST    = 1
      LONRES    = DISIZE
      MODELNAME = GET_MODELNAME()
      RESERVED  = ''
      SCALE     = DBLE( GET_CT_EMIS() ) + 1d-32
      NCHEMSTEP = DBLE( GET_CT_CHEM() ) + TINY( 1d0 ) !CDH for sea salt loss rat         
      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! Loop over ND03 diagnostic tracers
      DO M = 1, TMAX(3)

         ! Get ND03 tracer #
         N = TINDEX(3,M)

         ! Pick the proper array & dimensions
         IF ( N == 1 .or. N == 3 .or. N == 4  .or.
     &        N == 5 .or. N == 6 .or. N == 9 .or. N==13 .or. N==14
     &        .or. N==15 ) THEN
               
            !--------------------------------
            ! #1,3,4,5,6,9,13,14,15: Hg emissions
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N)

         ELSE IF ( N == 2 .or. N == 7 ) THEN

            !--------------------------------
            ! #2,7: Hg0, Hg2 ocean masses
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N) / SCALE
               
         ELSE IF ( N == 8 ) THEN
            
            !--------------------------------
            ! #8: Hg2 sinking loss rate
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N)

         ELSE IF ( N == 10 ) THEN
               
            !--------------------------------
            ! #10: Kw (piston velocity)
            ! Divide by # of emiss timesteps
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'cm/h'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N) / SCALE

         ELSE IF ( N == 11 ) THEN

            !--------------------------------
            ! #11: Hg(C) ocean mass
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N) / SCALE

         ELSE IF ( N == 12 ) THEN

            !--------------------------------
            ! #12: Converted to colloidal
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N)

         ELSE IF (N == 16)  THEN
            
            !--------------------------------
            ! #16: Flux-up (Hg0 volat from ocean)
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N)
            
         ELSE IF (N == 17)  THEN

            !--------------------------------
            ! #17: Flux-down (Hg0 dry dep to ocean)
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N)

         ELSE IF ( N == 18 ) THEN

            !--------------------------------
            ! #18: Production of Hg2 from Hg0
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = 1
            ARRAY(:,:,1:LMAX) = AD03_Hg2_Hg0(:,:,1:LMAX)

         ELSE IF ( N == 19 ) THEN

            !--------------------------------
            ! #20: Prod of Hg(II) from rxn w/OH
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = 2
            ARRAY(:,:,1:LMAX) = AD03_Hg2_OH(:,:,1:LMAX)

         ELSE IF ( N == 20 ) THEN

            !--------------------------------
            ! #21: Prod of Hg(II) from rxn w/O3
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = 3
            ARRAY(:,:,1:LMAX) = AD03_Hg2_O3(:,:,1:LMAX)
      
         ELSE IF ( N == 21 ) THEN
            
            !--------------------------------
            ! #21: Loss of Hg2 from rxn w/sea salt
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = N_Hg_CATS
            NN                = 4
            ARRAY(:,:,1:LMAX) = AD03_Hg2_SS(:,:,1:LMAX)
      
         ELSE IF ( N == 22 ) THEN
            
            !--------------------------------
            ! #22: Loss of Hg2 from rxn w/sea salt
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = '/s'
            LMAX              = 1
            NN                = 5
            ARRAY(:,:,1)      = AD03_Hg2_SSR(:,:) / NCHEMSTEP

         ELSE IF ( N == 23 ) THEN

            !--------------------------------
            ! #23: Prod of Hg(II) from rxn w/Br 
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = 6
            ARRAY(:,:,1:LMAX) = AD03_Hg2_Br(:,:,1:LMAX)
      
         ELSE IF ( N == 24 ) THEN

            !--------------------------------
            ! #24: Br concentration
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'molec/cm3'
            LMAX              = LD03
            NN                = 7
            ARRAY(:,:,1:LMAX) = AD03_Br(:,:,1:LMAX,1) / NCHEMSTEP
      
         ELSE IF ( N == 25 ) THEN

            !--------------------------------
            ! #25: Br concentration
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'molec/cm3'
            LMAX              = LD03
            NN                = 8
            ARRAY(:,:,1:LMAX) = AD03_Br(:,:,1:LMAX,2) / NCHEMSTEP

         ELSE

            !--------------------------------
            ! Otherwise skip to next N
            !--------------------------------
            CYCLE

         ENDIF

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LMAX,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LMAX) )
      ENDDO

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG03

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG03
!
!******************************************************************************
!  Subroutine INIT_DIAG03 allocates all module arrays (bmy, 1/21/05, 4/6/06)
!
!  NOTES:
!  (1 ) Now allocates AD03_Hg2_SS (eck, bmy, 4/6/06)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE TRACERID_MOD, ONLY : N_Hg_CATS
   
#     include "CMN_SIZE" 

      ! Local variables
      INTEGER :: AS
      
      !=================================================================
      ! INIT_DIAG03 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND03 == 0 ) THEN
         LD03 = 0
         RETURN
      ENDIF

      ! Get number of levels for 3-D arrays
      LD03 = MIN( ND03, LLPAR )

      ! 2-D array ("HG-SRCE")
      ALLOCATE( AD03( IIPAR, JJPAR, PD03-3 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03' )

      ! 3-D arrays ("PL-HG2-$")
      ALLOCATE( AD03_Hg2_Hg0( IIPAR, JJPAR, LD03 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_Hg0' )

      ALLOCATE( AD03_Hg2_OH( IIPAR, JJPAR, LD03 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_OH' )

      !cdh added diagnostic
      ALLOCATE( AD03_Hg2_Br( IIPAR, JJPAR, LD03 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_Br' )

      ALLOCATE( AD03_Hg2_O3( IIPAR, JJPAR, LD03 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_O3' )

      ALLOCATE( AD03_Hg2_SS( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_SS' )

      ALLOCATE( AD03_nat( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_nat' )

      !begin CDH for sea salt loss rate
      ALLOCATE( AD03_Hg2_SSR( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_SSR' )
      !End CDH

      !begin CDH for bromine
      ALLOCATE( AD03_Br( IIPAR, JJPAR, LLPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_BR' )
      !End CDH

      ! Zero arrays
      CALL ZERO_DIAG03

      ! Return to calling program
      END SUBROUTINE INIT_DIAG03

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG03
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG03 deallocates all module arrays 
!  (bmy, 1/21/05, 4/6/06)
!
!  NOTES:
!  (1 ) Now deallocates AD03_Hg2_SS (eck, bmy, 4/6/06)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG03 begins here!
      !=================================================================
      IF ( ALLOCATED( AD03         ) ) DEALLOCATE( AD03         ) 
      IF ( ALLOCATED( AD03_Hg2_Hg0 ) ) DEALLOCATE( AD03_Hg2_Hg0 )
      IF ( ALLOCATED( AD03_Hg2_OH  ) ) DEALLOCATE( AD03_Hg2_OH  )
      IF ( ALLOCATED( AD03_Hg2_Br  ) ) DEALLOCATE( AD03_Hg2_Br  )!cdh new diag
      IF ( ALLOCATED( AD03_Hg2_O3  ) ) DEALLOCATE( AD03_Hg2_O3  )
      IF ( ALLOCATED( AD03_Hg2_SS  ) ) DEALLOCATE( AD03_Hg2_SS  )
      IF ( ALLOCATED( AD03_nat     ) ) DEALLOCATE( AD03_nat     )
      IF ( ALLOCATED( AD03_Hg2_SSR ) ) DEALLOCATE( AD03_Hg2_SSR ) !CDH sea salt rate
      IF ( ALLOCATED( AD03_Br      ) ) DEALLOCATE( AD03_Br ) !CDH Bromine

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG03

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG03_MOD
