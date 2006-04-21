! $Id: pressure_mod.f,v 1.10 2006/04/21 15:40:05 bmy Exp $
      MODULE PRESSURE_MOD
!
!******************************************************************************
!  Module PRESSURE_MOD contains variables and routines which specify the grid 
!  box pressures for both hybrid or pure-sigma models.  This is necessary
!  for running GEOS-CHEM with the new GEOS-4/fvDAS meteorological fields.
!  (dsa, bmy, 8/27/02, 5/24/05)
!
!  The Hybrid ETA-coordinate (dsa, 8/27/02, 4/14/04)
!  ============================================================================
!  Pressure at layer edges are defined as follows:
!  
!     P(I,J,L) = AP(L) + ( BP(L) * PS(i,j) )
!  
!  where
!     PS = Psfc - PTOP (GEOS-1, GEOS-STRAT, GEOS-3) or
!     PS = Psfc        (GEOS-4) 
!        here Psfc is the true surface pressure [hPa]
!
!     AP has the same units as PS [hPa]
!     BP is a unitless constant given at layer edges.
!     In all cases  BP(LLPAR+1) = 0., BP(1) = 1.
!     The pressure at the model top is PTOP = AP(LLPAR+1)
!  
!  For a pure sigma system (GEOS-1, GEOS-STRAT, GEOS-3), this reduces to:
!     AP(L) = PTOP for all L
!     BP(L) = SIGE(L) (sigma at edges)
!
!  Module Variables:
!  ============================================================================
!  (1 ) AP   (REAL*8)         : "A" term for hybrid ETA coordinate
!  (2 ) BP   (REAL*8)         : "B" term for hybrid ETA coordinate
!  (3 ) PFLT (REAL*8)         : "Floating" surface pressure field
!
!  Module Routines:
!  ============================================================================
!  (1 ) GET_AP                : Returns "A" term for hybrid ETA coordinate
!  (2 ) GET_BP                : Returns "B" term for hybrid ETA coordinate
!  (3 ) SET_FLOATING_PRESSURE : Initializes PFLT w/ Psurface from "main.f"
!  (3 ) GET_PEDGE             : Returns pressure at bottom edge of box (I,J L)
!  (4 ) GET_PCENTER           : Returns pressure at center of box (I,J,L)
!  (5 ) INIT_PRESSURE         : Allocates and zeroes all module arrays
!  (6 ) CLEANUP_PRESSURE      : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by pressure_mod.f
!  ============================================================================
!  (1 ) error_mod.f    : Module containing I/O error and NaN check routines
!
!  NOTES:
!  (1 ) Be sure to check PFLT for NaN or Infinities (bmy, 8/27/02)
!  (2 ) Updated comments (bmy, 5/8/03)
!  (3 ) Updated format string for fvDAS (bmy, 6/19/03)
!  (4 ) Bug fix: use PFLT instead of PFLT-PTOP for GEOS-4 (bmy, 10/24/03)
!  (5 ) Modifications for 30L and 55L GEOS-4 grids (bmy, 11/3/03)
!  (6 ) Added parallel DO-loop in SET_FLOATING_PRESSURE (bmy, 4/14/04)
!  (7 ) Modified for GCAP and GEOS-5 grids (swu, bmy, 5/24/05)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "pressure_mod.f"
      !=================================================================

      ! Make everything PUBLIC ...
      PUBLIC
      
      ! ... except these variables
      PRIVATE :: AP, BP, PFLT

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      REAL*8, ALLOCATABLE :: AP(:)
      REAL*8, ALLOCATABLE :: BP(:)
      REAL*8, ALLOCATABLE :: PFLT(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_AP( L ) RESULT( AP_TEMP )
!
!******************************************************************************
!  Function GET_AP returns the "A" term [hPa] for the hybrid ETA coordinate.
!  (dsa, bmy, 8/20/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) L (INTEGER) : AP will be returned at the bottom edge of level L
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"  ! Size parameters
      
      ! Arguments
      INTEGER, INTENT(IN) :: L 

      ! Local variables
      REAL*8 :: AP_TEMP

      !=================================================================
      ! GET_AP begins here!
      !=================================================================      
      AP_TEMP = AP(L)

      ! Return to calling program
      END FUNCTION GET_AP

!------------------------------------------------------------------------------

      FUNCTION GET_BP( L ) RESULT( BP_TEMP )
!
!******************************************************************************
!  Function GET_BP returns the "B" term [unitless] for the hybrid ETA 
!  coordinate (dsa, bmy, 8/20/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) L (INTEGER) : BP will be returned at the bottom edge of level L
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"
      
      INTEGER, INTENT(IN)  :: L !edge level

      REAL*8 :: BP_TEMP

      !=================================================================
      ! GET_BP begins here!
      !=================================================================
      BP_TEMP = BP(L)

      ! Return to calling program
      END FUNCTION GET_BP

!------------------------------------------------------------------------------

      SUBROUTINE SET_FLOATING_PRESSURE( PS )
!
!******************************************************************************
!  Subroutine SET_FLOATING_PRESSURE initializes the floating pressure field
!  PFLT with a pressure from the main program.  This is needed to initialize 
!  and reset PFLT after transport. (dsa, bdf, bmy, 8/27/02, 4/14/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) PS (REAL*8) :: Array containing pressure with which to initialize PFLT
!
!  NOTES:
!  (1 ) Now check PFLT for NaN or Infinities (bmy, 8/27/02)
!  (2 ) Added parallel DO-loop (bmy, 4/14/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : CHECK_VALUE

#     include "CMN_SIZE"
   
      ! Arguments
      REAL*8, INTENT(IN) :: PS(IIPAR,JJPAR)

      ! Local variables
      INTEGER            :: I, J

      !=================================================================
      ! SET_FLOATING_PRESSURE begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Assign into PFLT array
         PFLT(I,J) = PS(I,J) 

         ! Check for NaN or Infinities
         CALL CHECK_VALUE( PFLT(I,J), (/I,J,0,0/), 
     &                    'PFLT',   'set_floating_pressure:1' )
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SET_FLOATING_PRESSURE

!------------------------------------------------------------------------------

      FUNCTION GET_PEDGE( I, J, L ) RESULT( PEDGE )
!
!******************************************************************************
!  Function GET_PEDGE returns the pressure at the bottom edge of level L.
!  (dsa, bmy, 8/20/02, 10/24/03)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) P (REAL*8 ) : P_surface - P_top (PS-PTOP)
!  (2 ) L (INTEGER) : Pressure will be returned at the bottom edge of level L
!
!  NOTES:
!  (1 ) Bug fix: use PFLT instead of PFLT-PTOP for GEOS-4 (bmy, 10/24/03)
!******************************************************************************
!
#     include "CMN_SIZE"  ! PTOP

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L 
      
      ! Return value
      REAL*8              :: PEDGE 

      !=================================================================
      ! GET_PEDGE begins here!
      !=================================================================
#if   defined( GEOS_4 )

      ! Here Ap is in [hPa] and Bp is [unitless].  
      ! For GEOS-4, we need to have PFLT as true surface pressure, 
      ! since Ap(1)=0 and Bp(1)=1.0.  This ensures that the true
      ! surface pressure will be returned for L=1. (bmy, 10/24/03)
      PEDGE = AP(L) + ( BP(L) * PFLT(I,J) )

#else 

      ! Here Ap is just PTOP in [hPa] and Bp are the sigma edges [unitless]  
      ! PFLT is the true surface pressure, so subtract PTOP from it
      PEDGE = AP(L) + ( BP(L) * ( PFLT(I,J) - PTOP ) )

#endif     
 
      ! Return to calling program
      END FUNCTION GET_PEDGE 

!------------------------------------------------------------------------------

      FUNCTION GET_PCENTER( I, J, L ) RESULT( PCENTER )
!
!******************************************************************************
!  Function GET_PEDGE returns the pressure at the bottom edge of level L.
!  (dsa, bmy, 8/20/02, 6/19/03)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) P_BOT (REAL*8 ) : P_surface - P_top (PS-PTOP)
!  (2 ) L     (INTEGER) : Pressure will be returned at the center of level L
!
!  NOTES:
!  (1 ) Updated format string for fvDAS (bmy, 6/19/03)
!******************************************************************************
!
#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! SIG

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L     
      
      ! Return value
      REAL*8              :: PCENTER 

      !=================================================================
      ! GET_PCENTER begins here!
      !=================================================================

      ! The pressure at the center of a grid-box is found
      ! by averaging the pressures at the box's two edges
      PCENTER = 0.5d0 * ( GET_PEDGE(I,J,L) + GET_PEDGE(I,J,L+1) )

      ! Return to calling program
      END FUNCTION GET_PCENTER

!------------------------------------------------------------------------------

      SUBROUTINE INIT_PRESSURE
!
!******************************************************************************
!  Subroutine INIT_PRESSURE allocates and initializes the AP and BP arrays.
!  It must be called in "main.f", after SIGE is defined.  GEOS-4 uses fvDAS, 
!  which requires the hybrid pressure system specified by the listed values 
!  of AP and BP, while earlier versions of GEOS use a pure sigma pressure
!  system.  GCAP met fields (based on GISS) also use a hybrid system. 
!  (dsa, swu, bmy, 8/20/02, 5/24/05)
!
!  NOTES:
!  (1 ) Now reference ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Now echo Ap, Bp to std output (bmy, 3/14/03)
!  (3 ) Now print LLPAR+1 levels for Ap, Bp.  Remove reference to SIGE, it's
!        obsolete.  Also now use C-preprocessor switch GRID30LEV instead of
!        IF statements to define vertical coordinates. (bmy, 11/3/03)
!  (4 ) Now modified for both GCAP & GEOS-5 vertical grids (swu, bmy, 5/24/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! LLPAR, PTOP

      ! Local Variables
      INTEGER :: AS
      INTEGER :: L

      ! Temporary array for GCAP met fields (swu, bmy, 5/24/05)
      REAL*8  :: SIGE_GCAP(LLPAR+1)  

      !=================================================================
      ! INIT_PRESSURE begins here!
      !=================================================================
      ALLOCATE( PFLT( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PFLT' )
      PFLT = 0d0

      ALLOCATE( AP( LLPAR + 1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AP' )
      AP = 1d0

      ALLOCATE( BP( LLPAR + 1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BP' )
      BP = 0d0

#if   defined( GEOS_4 ) || defined( GEOS_5 )
      
      !=================================================================
      ! GEOS-4 or GEOS-5 vertical coordinates (30 or 55 levels)
      !=================================================================

#if   defined( GRID30LEV )

      !--------------------------------
      ! GEOS-4 or GEOS-5 30 level grid
      !--------------------------------

      ! Ap [hPa] for 30 levels (31 edges)
      AP = (/  0.000000d0,   0.000000d0,  12.704939d0,  35.465965d0, 
     &        66.098427d0, 101.671654d0, 138.744400d0, 173.403183d0, 
     &       198.737839d0, 215.417526d0, 223.884689d0, 224.362869d0, 
     &       216.864929d0, 201.192093d0, 176.929993d0, 150.393005d0, 
     &       127.837006d0, 108.663429d0,  92.365662d0,  78.512299d0, 
     &        56.387939d0,  40.175419d0,  28.367815d0,  19.791553d0, 
     &         9.292943d0,   4.076567d0,   1.650792d0,   0.616779d0, 
     &         0.211349d0,   0.066000d0,   0.010000d0 /)

      ! Bp [unitless] for 30 levels (31 edges)
      BP = (/  1.000000d0,   0.985110d0,   0.943290d0,   0.867830d0, 
     &         0.764920d0,   0.642710d0,   0.510460d0,   0.378440d0, 
     &         0.270330d0,   0.183300d0,   0.115030d0,   0.063720d0, 
     &         0.028010d0,   0.006960d0,   0.000000d0,   0.000000d0, 
     &         0.000000d0,   0.000000d0,   0.000000d0,   0.000000d0, 
     &         0.000000d0,   0.000000d0,   0.000000d0,   0.000000d0, 
     &         0.000000d0,   0.000000d0,   0.000000d0,   0.000000d0, 
     &         0.000000d0,   0.000000d0,   0.000000d0 /)

#else

      !--------------------------------
      ! GEOS-4 or GEOS-5 55 level grid
      !--------------------------------

      ! AP [hPa] for 55 levels (56 edges)
      AP = (/ 0.000000d0,   0.000000d0,  12.704939d0,  35.465965d0, 
     &       66.098427d0, 101.671654d0, 138.744400d0, 173.403183d0,
     &      198.737839d0, 215.417526d0, 223.884689d0, 224.362869d0,
     &      216.864929d0, 201.192093d0, 176.929993d0, 150.393005d0,
     &      127.837006d0, 108.663429d0,  92.365662d0,  78.512299d0, 
     &       66.603378d0,  56.387939d0,  47.643932d0,  40.175419d0, 
     &       33.809956d0,  28.367815d0,  23.730362d0,  19.791553d0, 
     &       16.457071d0,  13.643393d0,  11.276889d0,   9.292943d0,
     &        7.619839d0,   6.216800d0,   5.046805d0,   4.076567d0, 
     &        3.276433d0,   2.620212d0,   2.084972d0,   1.650792d0,
     &        1.300508d0,   1.019442d0,   0.795134d0,   0.616779d0, 
     &        0.475806d0,   0.365041d0,   0.278526d0,   0.211349d0, 
     &        0.159495d0,   0.119703d0,   0.089345d0,   0.066000d0, 
     &        0.047585d0,   0.032700d0,   0.020000d0,   0.010000d0 /)

      ! BP [unitless] for 55 levels (56 edges)
      BP = (/  1.000000d0,  0.985110d0,   0.943290d0,   0.867830d0,
     &         0.764920d0,  0.642710d0,   0.510460d0,   0.378440d0,
     &         0.270330d0,  0.183300d0,   0.115030d0,   0.063720d0,
     &         0.028010d0,  0.006960d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0 /)


#endif

#elif defined( GEOS_3 ) 

      !=================================================================
      ! GEOS-3 vertical coordinates (30 or 48 levels)
      !=================================================================

#if   defined( GRID30LEV )

      !--------------------------------
      ! GEOS-3 30 level grid
      !--------------------------------

      ! AP [hPa] is just PTOP for a pure-sigma grid
      AP = PTOP

      ! BP [unitless] is just SIGE for a pure-sigma grid
      BP = (/ 1.000000d0, 0.997095d0, 0.991200d0, 0.981500d0, 
     &        0.967100d0, 0.946800d0, 0.919500d0, 0.884000d0, 
     &        0.839000d0, 0.783000d0, 0.718200d0, 0.647600d0, 
     &        0.574100d0, 0.500000d0, 0.427800d0, 0.359500d0, 
     &        0.297050d0, 0.241950d0, 0.194640d0, 0.155000d0, 
     &        0.122680d0, 0.096900d0, 0.076480d0, 0.047610d0, 
     &        0.029600d0, 0.018380d0, 0.007040d0, 0.002530d0, 
     &        0.000765d0, 0.000155d0, 0.000000d0 /)

#else

      !--------------------------------
      ! GEOS-3 48 level grid
      !--------------------------------

      ! AP [hPa] is just PTOP for a pure-sigma grid
      AP = PTOP
         
      ! BP [unitless] is just SIGE for a pure-sigma grid
      BP = (/ 1.000000d0, 0.997095d0, 0.991200d0, 0.981500d0,    
     &        0.967100d0, 0.946800d0, 0.919500d0, 0.884000d0,    
     &        0.839000d0, 0.783000d0, 0.718200d0, 0.647600d0,    
     &        0.574100d0, 0.500000d0, 0.427800d0, 0.359500d0,    
     &        0.297050d0, 0.241950d0, 0.194640d0, 0.155000d0,    
     &        0.122680d0, 0.096900d0, 0.076480d0, 0.060350d0,   
     &        0.047610d0, 0.037540d0, 0.029600d0, 0.023330d0,   
     &        0.018380d0, 0.014480d0, 0.011405d0, 0.008975d0,  
     &        0.007040d0, 0.005500d0, 0.004280d0, 0.003300d0,  
     &        0.002530d0, 0.001900d0, 0.001440d0, 0.001060d0,  
     &        0.000765d0, 0.000540d0, 0.000370d0, 0.000245d0, 
     &        0.000155d0, 9.20000d-5, 4.75000d-5, 1.76800d-5, 
     &        0.000000d0 /)
       
#endif  

#elif defined( GEOS_STRAT )

      !=================================================================
      ! GEOS-STRAT vertical coordinates (26 levels)
      !=================================================================

      ! AP [hPa] is just PTOP for a pure-sigma grid
      AP = PTOP

      ! BP [unitless] is just SIGE for a pure-sigma grid
      BP = (/ 1.d0, 0.987871d0, 0.954730d0, 0.905120d0, 0.845000d0, 
     &              0.780000d0, 0.710000d0, 0.639000d0, 0.570000d0, 
     &              0.503000d0, 0.440000d0, 0.380000d0, 0.325000d0, 
     &              0.278000d0, 0.237954d0, 0.202593d0, 0.171495d0, 
     &              0.144267d0, 0.121347d0, 0.102098d0, 0.085972d0, 
     &              0.072493d0, 0.061252d0, 0.051896d0, 0.037692d0, 
     &              0.019958d0, 0.000000d0 /)

#elif defined( GEOS_1 )

      !=================================================================
      ! GEOS-1 vertical coordinates (20 levels)
      !=================================================================

      ! AP [hPa] is just PTOP for a pure-sigma model
      AP = PTOP

      ! BP [unitless] is just SIGE for a pure-sigma grid
      BP = (/ 1.d0, 0.987871d0, 0.954730d0, 0.905120d0, 0.843153d0, 
     &              0.772512d0, 0.696448d0, 0.617779d0, 0.539000d0, 
     &              0.462000d0, 0.387500d0, 0.316500d0, 0.251000d0, 
     &              0.194500d0, 0.149800d0, 0.114600d0, 0.085500d0, 
     &              0.060500d0, 0.039000d0, 0.019000d0, 0.000000d0 /)

#elif defined( GCAP )

      !=================================================================
      ! GCAP vertical coordinates (23 levels)
      !=================================================================

      ! Define SIGMA edges from GISS 23L model
      SIGE_GCAP = (/ 1.0d0,           0.9712230d0,     0.9340528d0,     
     &               0.8800959d0,     0.8021583d0,     0.6714628d0, 
     &               0.5035971403d0,  0.3297362030d0,  0.1966426820d0,  
     &               0.1139088720d0,  0.0503597111d0,  0.0000000000d0, 
     &              -0.0395683460d0, -0.0764988065d0, -0.1124700233d0, 
     &              -0.1419664323d0, -0.1585131884d0, -0.1678657085d0, 
     &              -0.1743045598d0, -0.1781055182d0, -0.1793033630d0, 
     &              -0.1796822548d0, -0.1798187047d0, -0.1798536479d0 /)

      ! Convert SIGMA to AP and BP coordinates
      DO  L = 1, LLTROP
         AP(L) = 150d0 + ( SIGE_GCAP(L) * ( PTOP - 150d0 ) )
         BP(L) = SIGE_GCAP(L)
      ENDDO

      DO L = LLTROP+1, LLPAR+1  
         AP(L) = 150d0 + ( SIGE_GCAP(L) * ( 984d0 - 150d0 ) )
         BP(L) = 0d0
      ENDDO

#endif
      
      ! Echo info to std output
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'V E R T I C A L   G R I D   S E T U P'
      WRITE( 6, '(a,/)' ) 'INIT_PRESSURE: Vertical coordinates!'
      WRITE( 6, '( ''Ap '', /, 6(f11.6,1x) )' ) AP(1:LLPAR+1)
      WRITE( 6, '(a)'   )
      WRITE( 6, '( ''Bp '', /, 6(f11.6,1x) )' ) BP(1:LLPAR+1)
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE INIT_PRESSURE

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_PRESSURE
!
!******************************************************************************
!  Subroutine CLEANUP_PRESSURE deallocates all allocated arrays at the
!  end of a GEOS-CHEM model run. (dsa, bmy, 8/20/02)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_PRESSURE begins here!
      !=================================================================
      IF ( ALLOCATED( AP   ) ) DEALLOCATE( AP   )
      IF ( ALLOCATED( BP   ) ) DEALLOCATE( BP   )
      IF ( ALLOCATED( PFLT ) ) DEALLOCATE( PFLT )

      ! Return to calling program
      END SUBROUTINE CLEANUP_PRESSURE

!------------------------------------------------------------------------------

      END MODULE PRESSURE_MOD
