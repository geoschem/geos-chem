!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag03_mod
!
! !DESCRIPTION:  Module DIAG03\_MOD contains arrays and routines for archiving 
!  the ND03 diagnostic -- Hg emissions, mass, and production.
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG03_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !DEFINED PARAMETERS:
!
      INTEGER, PUBLIC, PARAMETER   :: PD03    = 18         ! Dim of AD03 array
      INTEGER, PUBLIC, PARAMETER   :: PD03_PL = 10         ! # of PL-HG2 diags

!
! !PUBLIC DATA MEMBERS:
!
      ! Scalars
      INTEGER, PUBLIC              :: ND03                 ! NDO3 on/off flag
      INTEGER, PUBLIC              :: LD03                 ! # of levels

      ! Arrays
!      REAL*4,  PUBLIC, ALLOCATABLE :: AD03(:,:,:)          ! Diagnostic arrays
!      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_Hg0(:,:,:)  !  for the prod/loss
!      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_Br(:,:,:)   !  and mass of 
!      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_OH(:,:,:)   !  various Hg
!      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_O3(:,:,:)   !  species
!      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_SS(:,:,:)   !
!      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_nat(:,:,:)      !
!      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_SSR(:,:)    !
!- eds 8/31/10 ----------------------------------------------------------
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03(:,:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_Hg0(:,:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_Br(:,:,:,:) 
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_OH(:,:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_O3(:,:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_SS(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_nat(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Hg2_SSR(:,:,:)
!-------------------------------------------------------------------------------
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_Br(:,:,:,:)     !
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_RGM(:,:,:)      !
      REAL*4,  PUBLIC, ALLOCATABLE :: AD03_PBM(:,:,:)      !

!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC :: ZERO_DIAG03
      PUBLIC :: WRITE_DIAG03
      PUBLIC :: INIT_DIAG03
      PUBLIC :: CLEANUP_DIAG03
!
! !REMARKS:
!  Nomenclature: 
!  ============================================================================
!  (1 ) Hg(0)  a.k.a. Hg0     : Elemental   mercury
!  (2 ) Hg(II) a.k.a. Hg2     : Divalent    mercury
!  (3 ) RGM a.k.a. Hg(II)gas  : Reactive (oxidized) gaseous mercury
!  (4 ) PBM a.k.a. Hg(II)P    : Reactive (oxidized) particulate mercury
!
! !REVISION HISTORY:
!  21 Jan 2005 - R. Yantosca - Initial version
!  (1 ) Updated for GCAP grid (bmy, 6/28/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Add 2 extra diagnostics to ND03. Set PD03=15.  (cdh, bmy, 12/15/05)
!  (4 ) Add loss of Hg2 by sea salt (eck, bmy, 4/6/06)
!  (5 ) Replace TINY(1d0) w/ 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (6 ) Updates to mercury simulation (ccc, 5/17/10)
!  (7 ) Added mercury tagged tracers (eds 9/2/10)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  07 Feb 2012 - E. Corbitt  - Added new diagnostics for tagged simulation.
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
! !IROUTINE: zero_diag03
!
! !DESCRIPTION: Subroutine ZERO\_DIAG03 zeroes all module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ZERO_DIAG03
!
! !USES:
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : N_Hg_CATS

      USE CMN_SIZE_MOD  ! Size parameters
! 
! !REVISION HISTORY: 
!  21 Jan 2005 - R. Yantosca - Initial version
!  (1 ) Now references N_Hg_CATS from "tracerid_mod.f".  Now zero AD03_Hg2_SS
!        array. (bmy, 4/6/06)
!  (2 ) Now use broadcast assignment and double precision 0D0 to zero arrays,
!        rather than nested DO loops and single precision 0E0. (cdh, 8/14/08)
!  (3 ) Now zeros RGM and PBM diagnostics. (hma 20100219)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L, N

      !=================================================================
      ! ZERO_DIAG03 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND03 == 0 ) RETURN

      ! Zero arrays
      AD03         = 0D0
      AD03_Hg2_Hg0 = 0D0
      AD03_Hg2_Br  = 0D0 
      AD03_Hg2_OH  = 0D0
      AD03_Hg2_O3  = 0D0
      AD03_Hg2_SS  = 0D0
      AD03_Hg2_SSR = 0D0 
      AD03_nat     = 0D0
      AD03_Br      = 0D0  
      AD03_RGM     = 0D0     
      AD03_PBM     = 0D0

      END SUBROUTINE ZERO_DIAG03
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag03
!
! !DESCRIPTION: Subroutine WRITE\_DIAG03 writes the ND03 diagnostic arrays to 
!  the binary punch file at the proper time.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG03
!
! !USES:
!
      USE BPCH2_MOD,    ONLY : BPCH2
      USE BPCH2_MOD,    ONLY : GET_MODELNAME
      USE BPCH2_MOD,    ONLY : GET_HALFPOLAR
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GRID_MOD,     ONLY : GET_XOFFSET
      USE GRID_MOD,     ONLY : GET_YOFFSET
      USE TIME_MOD,     ONLY : GET_CT_EMIS
      USE TIME_MOD,     ONLY : GET_DIAGb
      USE TIME_MOD,     ONLY : GET_DIAGe
      USE TIME_MOD,     ONLY : GET_CT_CHEM  
      USE TRACERID_MOD, ONLY : N_Hg_CATS
      USE TIME_MOD,     ONLY : GET_CT_DIAG, GET_Hg2_DIAG !H Amos, 20100218
      USE LOGICAL_MOD,  ONLY : LSPLIT !eds 9/7/10

      USE CMN_SIZE_MOD     ! Size parameters
      USE CMN_DIAG_MOD     ! TINDEX
!
! !REMARKS:
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  The following list is outdated and not reliable (cdh, 7/5/11)
!  (1 ) HG-SRCE  : Anthropogenic HG0 emission      : kg       : 1
!  (2 ) HG-SRCE  : Total mass of oceanic Hg0       : kg       : 1
!  (3 ) HG-SRCE  : Oceanic HgO emission            : kg       : 1
!  (4 ) HG-SRCE  : Land reemission                 : kg       : 1
!  (5 ) HG-SRCE  : Land natural emission           : kg       : 1
!  (6 ) HG-SRCE  : Anthropogenic Hg2 emission      : kg       : 1
!  (7 ) HG-SRCE  : Total mass of oceanic Hg2       : kg       : 1
!  (8 ) HG-SRCE  : Mass of Hg2 sunk in the ocean   : kg       : 1
!  (9 ) HG-SRCE  : Anthropogenic HgP emission      : kg       : 1
!  (10) HG-SRCE  : Henry's law piston velocity Kw  : cm/h     : em timesteps  (anls, redo)
!  (11) HG-SRCE  : Mass of Hg(P)                   : kg       : 1
!  (12) HG-SRCE  : Converted to Particulate        : kg       : 1
!  (13) HG-SRCE  : Biomass burning emissions       : kg       : 1
!  (14) HG-SRCE  : Emissions from vegetation       : kg       : 1
!  (15) HG-SRCE  : Emissions from soils            : kg       : 1
!  (16) HG-SRCE  : Flux-up Hg0 volat from ocean    : kg       : 1
!  (17) HG-SRCE  : Flux-down Hg0 dry dep to ocean  : kg       : 1
!  (18) PL-HG2-$ : Production of Hg2 from Hg0      : kg       : 1
!  (19) PL-HG2-$ : Production of Hg2 from rxn w/OH : kg       : 1
!  (20) PL-HG2-$ : Production of Hg2 from rxn w/O3 : kg       : 1
!  (21) PL-HG2-$ : Loss of Hg2 from rxn w/ seasalt : kg       : 1 
!  (22) PL-HG2-$ : Prod of Hg2 form rxn w/ Br      : kg       : 1
!  (23) PL-HG2-$ : Br concentration                : molec/cm3: 1
!  (24) PL-HG2-$ : Br concentration                : molec/cm3: 1
!  (27) PL-HG2-$ : Reactive gaseous mercury        : pptv     : 1
!  (28) PL-HG2-$ : Reactice particule mercury      : pptv     : 1
!
!  NOTES:
! 
! !REVISION HISTORY: 
!  21 Jan 2005 - R. Yantosca - Initial version
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
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: CENTER180, HALFPOLAR,   IFIRST
      INTEGER           :: JFIRST,    LFIRST,      LMAX
      INTEGER           :: M,         N,           NN
!- eds 9/1/10 ---------------------------------------------------------------
      INTEGER           :: NT
      REAL*4            :: ARRAY(IIPAR,JJPAR,LLPAR)
!      REAL*4            :: ARRAY(IIPAR,JJPAR,LLPAR,N_Hg_CATS) 
!----------------------------------------------------------------------------
      REAL*4            :: LONRES,    LATRES
      REAL*8            :: DIAGb,     DIAGe,       SCALE
      CHARACTER(LEN=20) :: MODELNAME 
      CHARACTER(LEN=40) :: CATEGORY,  RESERVED,    UNIT
      REAL*8            :: NCHEMSTEP 
      REAL*8            :: NDIAGSTEP       ! hma, for RGM and PBM
      REAL*8            :: NDIAGSTEP_Hg2   ! hma, for Fg and Fp
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
      NDIAGSTEP = DBLE( GET_CT_DIAG() ) + TINY( 1d0 ) ! for RGM and PBM
      NDIAGSTEP_Hg2 = DBLE ( GET_Hg2_DIAG() )+ TINY( 1d0 ) ! for Fg and Fp


      !=================================================================
      ! Different diagnostics for total vs. tagged simulations (eds)
      !=================================================================
      IF (LSPLIT) THEN

         !=================================================================
         ! Write data to the bpch file
         !=================================================================

         ! Loop over ND03 diagnostic tracers
!         DO M = 1, TMAX(3) !eds 9/9/10 -> should equal PD03_PL as of tag3

         ! Get ND03 tracer #
!         N = TINDEX(3,M)
         DO N = 1, 26 !eds 9/9/10 tag6

         ! Loop over tag IDs
         DO NT = 1, N_HG_CATS

         ! Pick the proper array & dimensions


         IF ( N == 1 )  THEN
           
            IF ( NT == 1 ) THEN

            !--------------------------------
            ! #1 Hg(0) Anthropogenic
            !--------------------------------
            CATEGORY          = 'HG0-ANTH'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)

            ELSE
              CYCLE
            ENDIF


         ELSE IF ( N == 2 )  THEN

            !--------------------------------
            ! #1 Hg(0) Ocean Mass    
            !--------------------------------
       	    CATEGORY   	      =	'HG0-AQUA'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT) / SCALE


         ELSE IF ( N == 3 ) THEN

            !--------------------------------
            ! #3 Net Ocean
            !--------------------------------
            CATEGORY          = 'HGNET-OC'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N        
!            NT                = ID_Hg0(NTG)
            ARRAY(:,:,1)      = AD03(:,:,N,NT)


         ELSE IF ( N == 4 ) THEN

            !--------------------------------
            ! #4  Land Prompt Recycling
            !--------------------------------
            CATEGORY          = 'HG0-RECY'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)


         ELSE IF ( N == 5 ) THEN

       	    IF ( NT == 1 ) THEN

            !--------------------------------
            ! #5  Geogenic
            !--------------------------------
            CATEGORY          = 'HG0-GEOG'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)

       	    ELSE
       	      CYCLE
       	    ENDIF


         ELSE IF ( N == 6 ) THEN

       	    IF ( NT == 1 ) THEN

            !--------------------------------
            ! #6 Hg(II) Anthropogenic
            !--------------------------------
            CATEGORY          = 'HG2-ANTH'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)

       	    ELSE
       	      CYCLE
       	    ENDIF


         ELSE IF ( N == 7 ) THEN

            !--------------------------------
            ! #7 Hg(II) Ocean Mass
            !--------------------------------
            CATEGORY          = 'HG2-AQUA'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT) /SCALE


         ELSE IF ( N == 8 ) THEN

       	    IF ( NT == 1 ) THEN
 
            !--------------------------------
            ! #8 Hg(II) Sinking
            !--------------------------------
            CATEGORY          = 'HG2-SINK'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)

            ELSE
              CYCLE
            ENDIF


         ELSE IF ( N == 9 ) THEN

       	    IF ( NT == 1 ) THEN

            !--------------------------------
            ! #9 Hg(P) Anthropogenic
            !--------------------------------
            CATEGORY          = 'HGP-ANTH'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
!            NT                = ID_HgP(NTG)
            ARRAY(:,:,1)      = AD03(:,:,N,NT)

       	    ELSE
       	      CYCLE
       	    ENDIF


         ELSE IF ( N == 10 ) THEN

            IF ( NT == 1 ) THEN

            !--------------------------------
            ! #10 Total Aqueous
            !--------------------------------
            CATEGORY          = 'HGT-AQUA'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)

            ELSE
              CYCLE
            ENDIF


         ELSE IF ( N == 11 ) THEN

            !--------------------------------
            ! #11 Hg(P) Aqueous Mass
            !--------------------------------
            CATEGORY          = 'HGP-AQUA'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)


         ELSE IF ( N == 12 ) THEN

            !--------------------------------
            ! #12 Conversion to Particle
            !--------------------------------
            IF ( NT == 1 ) THEN

            CATEGORY          = 'HGP-CONV'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)
 
            ELSE
              CYCLE
            ENDIF


         ELSE IF ( N == 13 ) THEN

       	    IF ( NT == 1 ) THEN

            !--------------------------------
            ! #13 Biomass Burning
            !--------------------------------
            CATEGORY          = 'HG0-BURN'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)

       	    ELSE
       	      CYCLE
       	    ENDIF


         ELSE IF ( N == 14 ) THEN

       	    IF ( NT == 1 ) THEN

            !--------------------------------
            ! #14 Vegetation Transpiration
            !--------------------------------
            CATEGORY          = 'HG0-VEGT'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)

       	    ELSE
       	      CYCLE
       	    ENDIF


         ELSE IF ( N == 15 ) THEN

       	    IF ( NT == 1 ) THEN

            !--------------------------------
            ! #15 Soils
            !--------------------------------
            CATEGORY          = 'HG0-SOIL'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)

       	    ELSE
       	      CYCLE
       	    ENDIF


         ELSE IF ( N == 16 ) THEN

            !--------------------------------
            ! #16 Gross Ocean Up Flux
            !--------------------------------
            CATEGORY          = 'HG0-FXUP'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)


         ELSE IF ( N == 17 ) THEN

            !--------------------------------
            ! #17 Gross Ocean Down Flux
            !--------------------------------
            CATEGORY          = 'HG0-FXDN'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)


         ELSE IF ( N == 18 ) THEN

            !--------------------------------
            ! #18 Snow
            !--------------------------------
            CATEGORY          = 'HG0-SNOW'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,NT)


         ELSE IF ( N == 19 ) THEN

            !--------------------------------
            ! #19: Net Ox of Hg(0) to Hg(II)
            !--------------------------------
            CATEGORY          = 'HG-NETOX'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Hg2_Hg0(:,:,1:LMAX,NT)



         ELSE IF ( N == 20 ) THEN

            !--------------------------------
            ! #20: Prod of Hg(II) from rxn w/OH
            !--------------------------------

            IF ( NT == 1 ) THEN

            CATEGORY          = 'HG2-OXOH'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Hg2_OH(:,:,1:LMAX,NT)

            ELSE
              CYCLE
            ENDIF


         ELSE IF ( N == 21 ) THEN

            !--------------------------------
            ! #21: Prod of Hg(II) from rxn w/O3
            !--------------------------------

            IF ( NT == 1 ) THEN

            CATEGORY          = 'HG2-OXO3'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Hg2_O3(:,:,1:LMAX,NT)

            ELSE
              CYCLE
            ENDIF

      
         ELSE IF ( N == 22 ) THEN

            !--------------------------------
            ! #22: Loss of Hg2 from rxn w/sea salt
            !--------------------------------
            CATEGORY          = 'HG2-SALT'
            UNIT              = 'kg'
            NN                = N
            LMAX              = 1
            ARRAY(:,:,1)      = AD03_Hg2_SS(:,:,NT)

      
         ELSE IF ( N == 23 ) THEN

            !--------------------------------
            ! #23: Loss of Hg2 from rxn w/sea salt
            !--------------------------------
            CATEGORY          = 'HG2-SSRX'
            UNIT              = '/s'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1) = AD03_Hg2_SSR(:,:,NT) / NCHEMSTEP


         ELSE IF ( N == 24 ) THEN

            !--------------------------------
            ! #24: Prod of Hg(II) from rxn w/Br 
            !--------------------------------
            CATEGORY          = 'HG2-OXBR'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Hg2_Br(:,:,1:LMAX,NT)

      
         ELSE IF ( N == 25 ) THEN

            IF ( NT == 1 ) THEN

            !--------------------------------
            ! #25: Br concentration
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'molec/cm3'
            LMAX              = LD03
            NN                = 1
            ARRAY(:,:,1:LMAX) = AD03_Br(:,:,1:LMAX,1) / NCHEMSTEP

            ELSE
               CYCLE
            ENDIF

      
         ELSE IF ( N == 26 ) THEN

            IF ( NT == 2 ) THEN

            !--------------------------------
            ! #26: BrO concentration
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'molec/cm3'
            LMAX              = LD03
            NN                = 2
            ARRAY(:,:,1:LMAX) = AD03_Br(:,:,1:LMAX,2) / NCHEMSTEP

            ELSE
       	       CYCLE
            ENDIF


         ELSE

            !--------------------------------
            ! Otherwise skip to next N
            !--------------------------------
            CYCLE

         ENDIF

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NT,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LMAX,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LMAX) )
      ENDDO

      ENDDO 

      ELSE

      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! Loop over ND03 HG-SRCE diagnostic tracers
      DO M = 1, TMAX(3)

         ! Get ND03 tracer #
         N = TINDEX(3,M)

         ! Pick the proper array & dimensions
         IF ( N == 2 .or. N == 7 .or. N == 10 .or. N == 11 ) THEN

            !--------------------------------
            ! #2,7,10,11: Hg0, Hg2, Hg(P), Hg_tot ocean masses
            ! Divide by # of emiss timesteps
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,1) / SCALE
               
         ELSE IF ( N <= 18 ) THEN
            
            !--------------------------------
            ! #1,3,4,5,6,9,13,14,15,18: Hg emissions
            ! #8: Hg2_tot sinking
            ! #12: Carbon sinking               !anls
            ! #16: Flux-up (Hg0 volat from ocean)
            ! #17: Flux-down (Hg0 dry dep to ocean)
            !--------------------------------
            CATEGORY          = 'HG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03(:,:,N,1)

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


      ! Loop over ND03 PL-HG2-$ diagnostics
      DO N=1, PD03_PL


         ! Pick array and units
         IF ( N == 1 ) THEN

            !--------------------------------
            ! #1: Production of Hg2 from Hg0
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Hg2_Hg0(:,:,1:LMAX,1)

         ELSE IF ( N == 2 ) THEN

            !--------------------------------
            ! #2: Prod of Hg(II) from rxn w/OH
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Hg2_OH(:,:,1:LMAX,1)

         ELSE IF ( N == 3 ) THEN

            !--------------------------------
            ! #3: Prod of Hg(II) from rxn w/O3
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Hg2_O3(:,:,1:LMAX,1)
      
         ELSE IF ( N == 4 ) THEN
            
            !--------------------------------
            ! #4: Loss of Hg2 from rxn w/sea salt
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = N_Hg_CATS
            NN                = N
            ARRAY(:,:,1)      = AD03_Hg2_SS(:,:,1)
      
         ELSE IF ( N == 5 ) THEN
            
            !--------------------------------
            ! #5: Loss of Hg2 from rxn w/sea salt
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = '/s'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD03_Hg2_SSR(:,:,1) / NCHEMSTEP

         ELSE IF ( N == 6 ) THEN

            !--------------------------------
            ! #6: Prod of Hg(II) from rxn w/Br 
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'kg'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Hg2_Br(:,:,1:LMAX,1)
      
         ELSE IF ( N == 7 ) THEN

            !--------------------------------
            ! #7: Br concentration
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'molec/cm3'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Br(:,:,1:LMAX,1) / NCHEMSTEP
      
         ELSE IF ( N == 8 ) THEN

            !--------------------------------
            ! #8: BrO concentration
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'molec/cm3'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_Br(:,:,1:LMAX,2) / NCHEMSTEP

!This is a duplicate of N==7, must be a bug (H Amos, 26 Oct 2011)
!         ELSE IF ( N == 9 ) THEN
!
!            !--------------------------------
!            ! #8: Br concentration
!            !--------------------------------
!            CATEGORY          = 'PL-HG2-$'
!            UNIT              = 'molec/cm3'
!            LMAX              = LD03
!            NN                = N
!            ARRAY(:,:,1:LMAX) = AD03_Br(:,:,1:LMAX,2) / NCHEMSTEP

         ELSE IF ( N == 9 ) THEN
            !--------------------------------
            ! #9: Reactive particulate mercury
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'pptv'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_PBM(:,:,1:LMAX) / NDIAGSTEP

         ELSE IF ( N == 10 ) THEN
            !--------------------------------
            ! #10: Reactive gaseous mercury
            !--------------------------------
            CATEGORY          = 'PL-HG2-$'
            UNIT              = 'pptv'
            LMAX              = LD03
            NN                = N
            ARRAY(:,:,1:LMAX) = AD03_RGM(:,:,1:LMAX) / NDIAGSTEP 

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

      END IF

      END SUBROUTINE WRITE_DIAG03
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag03
!
! !DESCRIPTION: Subroutine INIT\_DIAG03 allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG03
!
! !USES:
!
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE TRACERID_MOD, ONLY : N_Hg_CATS
   
      USE CMN_SIZE_MOD 
! 
! !REVISION HISTORY: 
!  21 Jan 2005 - R. Yantosca - Initial version
!  (1 ) Now allocates AD03_Hg2_SS (eck, bmy, 4/6/06)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  (3 ) Adapted for tagged tracers (eds 2/7/12)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

!- eds 8/31/10 -----------------------------------------------------
!      ! 2-D array ("HG-SRCE")
!      ALLOCATE( AD03( IIPAR, JJPAR, PD03 ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03' )

      ! 3-D array ("HG-SRCE")
      ALLOCATE( AD03( IIPAR, JJPAR, PD03, N_HG_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03' )
!------------------------------------------------------------------

!- eds 8/31/10 ----------------------------------------------------
!      ! 3-D arrays ("PL-HG2-$")
!      ALLOCATE( AD03_Hg2_Hg0( IIPAR, JJPAR, LD03 ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_Hg0' )

      ! 4-D arrays ("PL-HG2-$")
      ALLOCATE( AD03_Hg2_Hg0( IIPAR, JJPAR, LD03, N_HG_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_Hg0' )
!-----------------------------------------------------------------

! eds 8/31/10 ----------------------------------------------------
!      ALLOCATE( AD03_Hg2_OH( IIPAR, JJPAR, LD03 ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_OH' )

      ALLOCATE( AD03_Hg2_OH( IIPAR, JJPAR, LD03, N_HG_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_OH' )
!-----------------------------------------------------------------

! eds 8/31/10 ----------------------------------------------------
!      ALLOCATE( AD03_Hg2_Br( IIPAR, JJPAR, LD03 ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_Br' )

      ALLOCATE( AD03_Hg2_Br( IIPAR, JJPAR, LD03, N_HG_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_Br' )
!-----------------------------------------------------------------

! eds 8/31/10 ----------------------------------------------------
!      ALLOCATE( AD03_Hg2_O3( IIPAR, JJPAR, LD03 ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_O3' )

      ALLOCATE( AD03_Hg2_O3( IIPAR, JJPAR, LD03, N_HG_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_O3' )
!-----------------------------------------------------------------

      ALLOCATE( AD03_Hg2_SS( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_SS' )

      ALLOCATE( AD03_nat( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_nat' )

! eds 8/31/10 ----------------------------------------------------
!      ALLOCATE( AD03_Hg2_SSR( IIPAR, JJPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_SSR' )

      ALLOCATE( AD03_Hg2_SSR( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_Hg2_SSR' )
!-----------------------------------------------------------------

      ALLOCATE( AD03_Br( IIPAR, JJPAR, LLPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_BR' )
      !End CDH

! add hma 20100216-------------------------------------------------------
      ! Reactive particulate mercury
      ALLOCATE( AD03_PBM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_PBM' )

      ! Reactive gaseous mercury
      ALLOCATE( AD03_RGM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03_RGM' )
!------------------------------------------------------------------------

      ! Zero arrays
      CALL ZERO_DIAG03

      END SUBROUTINE INIT_DIAG03
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag03
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG03 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG03
! 
! !REVISION HISTORY: 
!  21 Jan 2005 - R. Yantosca - Initial version
!  (1 ) Now deallocates AD03_Hg2_SS (eck, bmy, 4/6/06)
!  (2 ) Now deallocates AD03_PBM, AD03_RGM (hma 20100216)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DIAG03 begins here!
      !=================================================================
      IF ( ALLOCATED( AD03         ) ) DEALLOCATE( AD03         ) 
      IF ( ALLOCATED( AD03_Hg2_Hg0 ) ) DEALLOCATE( AD03_Hg2_Hg0 )
      IF ( ALLOCATED( AD03_Hg2_OH  ) ) DEALLOCATE( AD03_Hg2_OH  )
      IF ( ALLOCATED( AD03_Hg2_Br  ) ) DEALLOCATE( AD03_Hg2_Br  )
      IF ( ALLOCATED( AD03_Hg2_O3  ) ) DEALLOCATE( AD03_Hg2_O3  )
      IF ( ALLOCATED( AD03_Hg2_SS  ) ) DEALLOCATE( AD03_Hg2_SS  )
      IF ( ALLOCATED( AD03_nat     ) ) DEALLOCATE( AD03_nat     )
      IF ( ALLOCATED( AD03_Hg2_SSR ) ) DEALLOCATE( AD03_Hg2_SSR ) 
      IF ( ALLOCATED( AD03_Br      ) ) DEALLOCATE( AD03_Br      ) 
! added hma 20100216-----------------------------------------------------------
      IF ( ALLOCATED( AD03_PBM     ) ) DEALLOCATE( AD03_PBM     )
      IF ( ALLOCATED( AD03_RGM     ) ) DEALLOCATE( AD03_RGM     )
!------------------------------------------------------------------------------

      END SUBROUTINE CLEANUP_DIAG03
!EOC
      END MODULE DIAG03_MOD
