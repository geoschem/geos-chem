!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag04_mod.F
!
! !DESCRIPTION: Module DIAG04\_MOD contains arrays and routines for archiving
!  the ND04 diagnostic -- CO2 emissions and fluxes.
!\\
!\\
! !INTERFACE:
!
      MODULE DIAG04_MOD
!
! !USES:
!
      USE PRECISION_MOD     ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      ! Scalars
      INTEGER              :: ND04, LD04
      INTEGER, PARAMETER   :: PD04 = 1

      ! Arrays
      REAL*4,  ALLOCATABLE :: AD04_chem(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_DIAG04
      PUBLIC :: INIT_DIAG04
      PUBLIC :: WRITE_DIAG04
      PUBLIC :: ZERO_DIAG04
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%   BUYER BEWARE! Tagged CO2 tracers only work for 2 x 2.5 grid!   %%%
!  %%%   Someone will have to make this more general later on...        %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                             .
! !REVISION HISTORY:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (3 ) Modified for ship emissions (2-D), aircraft emissions (3-D) and
!       chemical source for CO2 (3-D) (RayNassar, 2009-12-23)
!  20 May 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  07 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
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
! !IROUTINE: zero_diag04
!
! !DESCRIPTION: Subroutine ZERO\_DIAG04 zeroes the ND04 diagnostic array.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ZERO_DIAG04
!
! !USES:
!
      USE CMN_SIZE_MOD  ! Size parameters
!
! !REVISION HISTORY:
!  26 Jul 2005 - R. Yantosca - Initial version
!  18 May 2010 - R. Nassar   - Also zero AD04_PLANE, AD04_CHEM arrays
!  18 May 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! ZERO_DIAG04 begins here!
      !=================================================================
#if defined( BPCH_DIAG )

      ! Exit if ND04 is turned off
      IF ( ND04 == 0 ) RETURN

      ! Zero 2-D array (for N=7 tracers) and 3-D plane and chem arrays
      AD04_chem(:,:,:)  = 0e0
#endif

      END SUBROUTINE ZERO_DIAG04
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag04
!
! !DESCRIPTION: Subroutine WRITE\_DIAG04 writes the ND04 diagnostic arrays
!  to the binary punch file at the proper time.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG04( Input_Opt, State_Grid )
!
! !USES:
!
#if defined( BPCH_DIAG )
      USE BPCH2_MOD,      ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
      USE FILE_MOD,       ONLY : IU_BPCH
#endif
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Grid_Mod, ONLY : GrdState
      USE TIME_MOD,       ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_DIAG_MOD  ! TINDEX
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input options
      TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !REMARKS:
!   # : Field     : Description                  : Units       : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) CO2-SRCE  : CO2 fossil fuel emissions    : molec/cm2/s : SCALE
!  (2 ) CO2-SRCE  : CO2 ocean emissions          : molec/cm2/s : SCALE
!  (3 ) CO2-SRCE  : CO2 balanced biosphere       : molec/cm2/s : SCALE
!  (4 ) CO2-SRCE  : CO2 biomass emissions        : molec/cm2/s : SCALE
!  (5 ) CO2-SRCE  : CO2 biofuel emissions        : molec/cm2/s : SCALE
!  (6 ) CO2-SRCE  : CO2 net terrestrial exchange : molec/cm2/s : SCALE
!  (7 ) CO2-SRCE  : CO2 ship emissions           : molec/cm2/s : SCALE
!  (8 ) CO2-SRCE  : CO2 aircraft emissions (3-D) : molec/cm2/s : SCALE
!  (9 ) CO2-SRCE  : CO2 chemical source (3-D)    : molec/cm2/s : SCALE
!  (10) CO2-SRCE  : CO2 chem source surf correct : molec/cm2/s : SCALE!
!
! !REVISION HISTORY:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!  18 May 2010 - R. Nassar   - Now write out AD04_PLANE, ADO4_CHEM
!  18 May 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG )
      INTEGER            :: CENTER180, HALFPOLAR, IFIRST, JFIRST
      INTEGER            :: LFIRST,    LMAX,      M,      N
      REAL(f4)           :: ARRAY(State_Grid%NX,State_Grid%NY,
     &                            State_Grid%NZ)
      REAL(f4)           :: LONRES,    LATRES
      REAL(f8)           :: DIAGb,     DIAGe
      REAL(fp)           :: SCALE
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY,  RESERVED,    UNIT

      !=================================================================
      ! WRITE_DIAG04 begins here!
      !=================================================================

      ! Exit if ND04 is turned off
      IF ( ND04 == 0 ) RETURN

      ! Initialize
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = GET_HALFPOLAR()
      IFIRST    = State_Grid%XMinOffset + 1
      JFIRST    = State_Grid%YMinOffset + 1
      LATRES    = State_Grid%DY
      LFIRST    = 1
      LONRES    = State_Grid%DX
      MODELNAME = GET_MODELNAME( Input_Opt, State_Grid )
      RESERVED  = ''
      SCALE     = DBLE( GET_CT_EMIS() ) + 1e-32_fp

      !=================================================================
      ! Write data to the bpch file
      ! Note: if any of the ARRAY or AD04* dimensions are wrong, the
      ! run will crash with "ERROR RUNNING GEOS-CHEM" at the end.
      !=================================================================

      ! Loop over ND04 diagnostic tracers
      DO M = 1, TMAX(4)

         ! Get quantities
         N            = TINDEX(4,M)

         IF ( N == 1 ) THEN

            CATEGORY     = 'CO2-SRCE'
            UNIT         = 'molec/cm2/s'
            LMAX = LD04
            ARRAY(:,:,1:LMAX) = AD04_chem(:,:,1:LMAX) / SCALE

         ELSE

            CYCLE

         ENDIF

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               State_Grid%NX, State_Grid%NY, LMAX, IFIRST,
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LMAX) )

      ENDDO
#endif

      END SUBROUTINE WRITE_DIAG04
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag04
!
! !DESCRIPTION: Subroutine INIT\_DIAG04 allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG04( State_Grid )
!
! !USES:
!
      USE ERROR_MOD,      ONLY : ALLOC_ERR
      USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
      TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  26 Jul 2005 - R. Yantosca - Initial version
!  18 May 2010 - R. Nassar   - Now initialize AD04_PLANE, AD04_CHEM
!  18 May 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_DIAG04 begins here!
      !=================================================================
#if defined( BPCH_DIAG )

      ! Exit if ND04 is turned off
      IF ( ND04 == 0 ) RETURN

      ! Get number of levels for 3-D arrays
      LD04 = MIN( ND04, State_Grid%NZ )

      ! 3-D arrays ("CO2-SRCE")
      ALLOCATE( AD04_chem(State_Grid%NX,State_Grid%NY,State_Grid%NZ),
     &          STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD04_chem' )

      ! Zero arrays
      CALL ZERO_DIAG04
#endif

      END SUBROUTINE INIT_DIAG04
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag04
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG04 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG04
!
! !REVISION HISTORY:
!  26 Jul 2005 - R. Yantosca - Initial version
!  18 May 2010 - R. Nassar   - Now ce
!  18 May 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DIAG04 begins here!
      !=================================================================
#if defined( BPCH_DIAG )
      IF ( ALLOCATED( AD04_chem  ) ) DEALLOCATE( AD04_chem  )
#endif

      END SUBROUTINE CLEANUP_DIAG04
!EOC
      END MODULE DIAG04_MOD
