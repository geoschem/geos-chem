! $Id: global_ch4_mod.f,v 1.6 2005/09/02 15:17:12 bmy Exp $
      MODULE GLOBAL_CH4_MOD
!
!******************************************************************************
!  Module GLOBAL_CH4_MOD contains variables and routines for simulating
!  CH4 chemistry in the troposphere (jsw, bnd, bmy, 1/17/01, 8/16/05)
!
!  Module Variables:
!  =========================================================================== 
!  (1 ) N_CH4      (INTEGER) : Number of budget items in TCH4
!  (2 ) BAIRDENS   (REAL*8 ) : Array for air density [molec/cm3]
!  (3 ) BOH        (REAL*8 ) : Array for OH values [molec/cm3]
!  (4 ) COPROD     (REAL*8 ) : Array for zonal mean P(CO) [v/v/s]
!  (5 ) PAVG       (REAL*8 ) : Array for 24-h avg surface pressure [mb]
!  (6 ) TAVG       (REAL*8 ) : Array for 24-h avg temperature [K]
!  (7 ) TCH4       (REAL*8 ) : Array for CH4 budget (N_CH4 items)
!  (8 ) NCMSALTS   (INTEGER) : # of altitudes for CMS climatological OH
!  (9 ) NCMSLATS   (INTEGER) : # of latitudes for CMS climatological OH
!  (10) CMSALTS    (REAL*8 ) : Altitude values for CMS climatological OH
!  (11) CMSLATS    (REAL*8 ) : Latitude values for CMS climatological OH
!  (12) AVGOH      (REAL*8 ) : Array for CMS climatological OH [molec/cm3]
!  (13) FMOL_CH4   (REAL*8 ) : Molecular weight of CH4 [kg/mole]
!  (14) XNUMOL_CH4 (REAL*8 ) : Molecules CH4 / kg CH4
!
!  Module Routines: 
!  =========================================================================== 
!  (1 ) GET_GLOBAL_CH4 : Computes latitudinal, yearly CH4 gradient
!  (2 ) CH4_AVGTP      : Computes 24-h average pressure & temperature
!  (3 ) EMISSCH4       : Handles CH4 emissions
!  (4 ) CHEMCH4        : Handles CH4 chemistry (various sinks)
!  (5 ) CLIMATOL_OH    : Reads in CMS climatological OH field
!  (6 ) INTERPOH       : Interpolates CMS clim. OH to GEOS-CHEM resolution
!  (7 ) CH4_DECAY      : Computes decay of CH4 w/ OH in the troposphere
!  (8 ) CO_OHSAVE      : Computes CH3CCl3 lifetime based on OH concentration
!  (9 ) CH4_STRAT      : Computes loss of CH4 in the stratosphere
!  (10) CH4_BUDGET     : Computes global CH4 budgets, sources & sinks
!  (11) SUM_CH4        : Sums a sub-region of the TCH4 budget array
!  (12) INIT_CH4       : Allocates and zeroes module arrays
!  (13) CLEANUP_CH4    : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by global_ch4_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module containing routines for binary punch file I/O
!  (2 ) diag_mod.f     : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) dao_mod.f      : Module containing arrays for DAO met fields
!  (4 ) diag_oh_mod.f  : Module containing arrays for mean OH & CH3CCl3 life
!  (4 ) error_mod.f    : Module containing NaN and other error check routines
!  (5 ) grid_mod.f     : Module containing horizontal grid information
!  (6 ) pressure_mod.f : Module containing routines to compute P(I,J,L) 
!  (7 ) time_mod.f     : Module containing routines to compute date & time
!
!  NOTES:
!  (1 ) Merged routines from jsw's CH4 code  into "global_ch4_mod.f" 
!        (bmy, 1/16/01)
!  (2 ) XNUMOL_CH4 and TCH4 have to be public - all other variables can
!        be made private, so as not to conflict with other common-block
!        definitions (bmy, 1/17/01)
!  (3 ) Minor fixes from jsw added (jsw, bmy, 2/17/01)
!  (4 ) Removed some F90 module references from EMISSCH4 (bmy, 3/20/01)
!  (5 ) Eliminate obsolete commented-out code (bmy, 4/20/01)
!  (6 ) Updated comments (bmy, 9/4/01)
!  (7 ) Fixes for binary punch file in READ_COPROD (bmy, 9/26/01)
!  (8 ) Removed obsolete code from READ_COPROD (bmy, 10/24/01)
!  (9 ) Minor bug fixes for compilation on ALPHA (bmy, 11/15/01)
!  (10) Eliminate obsolete code from 11/01 (bmy, 2/27/02)
!  (11) Now eliminate PS from the arg list to CH4_AVGTP (4/11/02)
!  (12) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (13) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (14) Now reference "file_mod.f".  Also removed obsolete code. (bmy, 6/27/02)
!  (15) Now references "pressure_mod.f" (bmy, 8/21/02)
!  (16) Now reference AD and T from "dao_mod.f".  Now reference "error_mod.f".
!        Remove obsolete code from various routines.  Remove reference to
!        header file "comtrid.h" -- it's not used. (bmy, 11/6/02)
!  (17) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (18) Now references "grid_mod.f" and "time_mod.f" (bmy, 3/27/03)
!  (19) Updates to GET_GLOBAL_CH4 (bmy, 7/1/03)
!  (20) Now references "directory_mod.f", "tracer_mod.f", and "diag_oh_mod.f"
!        (bmy, 7/20/04)
!  (21) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!******************************************************************************
!     
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "global_ch4_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE :: N_CH4,    BAIRDENS,  BOH    
      PRIVATE :: COPROD,   PAVG,      TAVG      
      PRIVATE :: NSEAS,    NCMSALTS,  NCMSLATS 
      PRIVATE :: CMSALTS,  CMSLATS,   AVGOH
      PRIVATE :: FMOL_CH4 
      
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Number of CH4 budget types 
      INTEGER, PARAMETER   :: N_CH4 = 12
                            
      ! Various arrays      
      REAL*8,  ALLOCATABLE :: BAIRDENS(:,:,:)
      REAL*8,  ALLOCATABLE :: BOH(:,:,:)
      REAL*8,  ALLOCATABLE :: COPROD(:,:,:)
      REAL*8,  ALLOCATABLE :: PAVG(:,:,:)
      REAL*8,  ALLOCATABLE :: TAVG(:,:,:)
      REAL*8,  ALLOCATABLE :: TCH4(:,:,:,:)

      ! For Clarisa's Climatological OH
      INTEGER, PARAMETER   :: NSEAS    = 4
      INTEGER, PARAMETER   :: NCMSALTS = 7
      INTEGER, PARAMETER   :: NCMSLATS = 24
      
      REAL*8               :: CMSALTS(NCMSALTS) =
     &    (/ 1000d0, 900d0, 800d0, 700d0, 500d0, 300d0, 200d0 /)

      REAL*8               :: CMSLATS(NCMSLATS) =
     &    (/ 90d0,  84d0,  76d0,  68d0,  60d0,  52d0,  44d0,  36d0, 
     &       28d0,  20d0,  12d0,   4d0,  -4d0, -12d0, -20d0, -28d0,
     &      -36d0, -44d0, -52d0, -60d0, -68d0, -76d0, -84d0, -90d0 /)

      REAL*8,  ALLOCATABLE  :: AVGOH(:,:,:)

      ! FMOL_CH4   =  kg CH4    / mole CH4
      ! XNUMOL_CH4 =  molec CH4 / kg CH4
      REAL*8, PARAMETER    :: FMOL_CH4   = 16d-3
      REAL*8, PARAMETER    :: XNUMOL_CH4 = 6.0221d23 / 16d-3

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_GLOBAL_CH4( THISYEAR, VARIABLE_CH4, 
     &                           A3090S, A0030S, A0030N, A3090N )
!
!******************************************************************************
!  Subroutine GET_GLOBAL_CH4 computes the latitudinal gradient in CH4
!  corresponding to year (jsw, bnd, bmy, 1/3/01, 7/1/03)
!
!  Arguments as Input:
!  ===========================================================================
!  (1) THISYEAR     (INTEGER) : Current month number (1-12)
!  (2) VARIABLE_CH4 (LOGICAL) : Flag for selecting variable or constant CH4
!
!  Arguments as Output:
!  ===========================================================================
!  (3 ) A3090S      (REAL*8 ) : CH4 concentration [ppbv], 90S - 30S lat
!  (4 ) A0030S      (REAL*8 ) : CH4 concentration [ppbv], 30S - 00  lat
!  (5 ) A0030N      (REAL*8 ) : CH4 concentration [ppbv], 00  - 30N lat
!  (6 ) A3090N      (REAL*8 ) : CH4 concentration [ppbv], 30N - 90N lat
!
!  NOTES:
!  (1 ) GET_GLOBAL_CH4 only has to be called at the start of the new year,
!        as long as A3090S, A0030S, A0030N, A3090N are saved in the
!        calling program (bmy, 1/3/01)
!  (2 ) Also need to compute yearly gradients for CH4 beyond 1997 --
!        will do this later (bmy, 1/3/01)
!  (3 ) Bug fix: add missing comma to FORMAT statement (bmy, 3/23/03)
!  (4 ) Place WRITE statments w/in an !$OMP CRITICAL block, so as to make
!        sure that only one processor at a time writes them.  Also now use
!        F90 REPEAT intrinsic function.  Also replaced old CH4 gradient values
!        with updated values for 1983-2001.  Use data for 2001 as a proxy for
!        years past 2001, since data for those years has not been reported
!        yet. (mje, bmy, 7/7/03)
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)  :: THISYEAR
      LOGICAL, INTENT(IN)  :: VARIABLE_CH4
      REAL*8,  INTENT(OUT) :: A3090S, A0030S, A0030N, A3090N

      !=================================================================
      ! New methane data from 1983-2001 (mje, bmy, 7/7/03)
      !
      ! Methane measurements are from CMDL website:
      ! ftp://140.172.192.211/ccg/ch4/flask/month
      ! 
      ! Measurements includes all sites other than:
      ! BAL BSC HUN MHD OXK TAP SEY IZO KUM MID ASK
      !
      ! Sites are separated into 4 latitude bands:
      !    (1) 90S - 30S;  (2) 30S - 00S;  
      !    (3) 00N - 30N;  (4) 30N - 90N
      ! 
      ! Bob Yantosca (bmy@io.harvard.edu) maintains the archive 
      ! of the IDL code needed to process the methane data.
      !=================================================================
      IF ( VARIABLE_CH4 ) THEN

         ! Select latitudinal CH4 gradient by year...
         SELECT CASE ( THISYEAR )

            CASE( 1983 )
               A3090S = 1559.89d0
               A0030S = 1575.68d0
               A0030N = 1627.04d0
               A3090N = 1682.40d0

            CASE( 1984 )
               A3090S = 1578.59d0
               A0030S = 1587.03d0
               A0030N = 1635.20d0
               A3090N = 1702.69d0

            CASE( 1985 )
               A3090S = 1588.78d0
               A0030S = 1600.98d0
               A0030N = 1648.02d0
               A3090N = 1716.23d0
               
            CASE( 1986 )
               A3090S = 1598.28d0  
               A0030S = 1612.76d0  
               A0030N = 1664.98d0  
               A3090N = 1731.23d0

            CASE( 1987 )
               A3090S = 1611.65d0  
               A0030S = 1622.34d0  
               A0030N = 1681.88d0  
               A3090N = 1741.44d0

            CASE( 1988 )
               A3090S = 1620.31d0  
               A0030S = 1634.43d0  
               A0030N = 1691.88d0  
               A3090N = 1753.92d0

            CASE( 1989 )
               A3090S = 1634.89d0  
               A0030S = 1647.15d0  
               A0030N = 1699.20d0  
               A3090N = 1759.64d0

            CASE( 1990 )
               A3090S = 1643.58d0  
               A0030S = 1653.97d0  
               A0030N = 1712.33d0  
               A3090N = 1769.97d0

            CASE( 1991 )
               A3090S = 1654.38d0  
               A0030S = 1665.13d0  
               A0030N = 1722.64d0  
               A3090N = 1779.76d0

            CASE( 1992 )
               A3090S = 1668.22d0  
               A0030S = 1673.40d0  
               A0030N = 1732.30d0  
               A3090N = 1786.76d0

            CASE( 1993 )
               A3090S = 1667.04d0  
               A0030S = 1677.26d0  
               A0030N = 1733.96d0  
               A3090N = 1790.82d0

            CASE( 1994 )
               A3090S = 1670.85d0  
               A0030S = 1681.07d0  
               A0030N = 1740.88d0  
               A3090N = 1797.05d0

            CASE( 1995 )
               A3090S = 1681.00d0  
               A0030S = 1689.19d0  
               A0030N = 1751.25d0  
               A3090N = 1802.51d0

            CASE( 1996 )
               A3090S = 1682.23d0  
               A0030S = 1690.72d0  
               A0030N = 1751.64d0  
               A3090N = 1805.18d0
            
            CASE( 1997 )
               A3090S = 1687.94d0  
               A0030S = 1693.35d0  
               A0030N = 1755.41d0  
               A3090N = 1805.92d0

            CASE( 1998 )
               A3090S = 1696.98d0  
               A0030S = 1703.54d0  
               A0030N = 1764.94d0  
               A3090N = 1820.58d0

            CASE( 1999 )
               A3090S = 1705.64d0  
               A0030S = 1714.18d0  
               A0030N = 1769.83d0  
               A3090N = 1823.48d0

            CASE( 2000 )
               A3090S = 1707.14d0  
               A0030S = 1715.63d0  
               A0030N = 1769.11d0  
               A3090N = 1822.85d0

            CASE( 2001 )
               A3090S = 1705.68d0  
               A0030S = 1709.52d0  
               A0030N = 1767.51d0  
               A3090N = 1822.53d0

            ! Use 2001 data as the default for years past 2001, until we 
            ! can get actual data for these years (bmy, 7/3/03)
            CASE DEFAULT
               A3090S = 1705.68d0  
               A0030S = 1709.52d0  
               A0030N = 1767.51d0  
               A3090N = 1822.53d0

         END SELECT

      ELSE
         
         ! ...otherwise assume constant global CH4
         A3090S = 1700.0d0
         A0030S = 1700.0d0
         A0030N = 1700.0d0
         A3090N = 1700.0d0
         
      ENDIF

      !=================================================================
      ! Print the latitudinal CH4 gradient for this year to stdout
      !=================================================================
!$OMP CRITICAL
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 105   ) THISYEAR
 105  FORMAT( 'GET_GLOBAL_CH4: YEAR = ', i4 )

      WRITE( 6, 110 ) A3090N, A0030N, A0030S, A3090S 
 110  FORMAT( 'CH4 (90N - 30N) : ', f7.1, ' [ppbv]', /,
     &        'CH4 (30N - 00 ) : ', f7.1, ' [ppbv]', /,
     &        'CH4 (00  - 30S) : ', f7.1, ' [ppbv]', /,
     &        'CH4 (30S - 90S) : ', f7.1, ' [ppbv]' )

      ! Indicate to the log file if we are using CH4 gradient data
      ! from 2001 as a proxy for years past 2001 (mje, bmy, 7/7/03)
      IF ( THISYEAR > 2001 ) THEN
         WRITE( 6, 115 ) 
 115     FORMAT( /, 'Using CH4 gradient data from 2001 as a proxy',
     &           /, 'since 2001 is the last year with reported data!' )
      ENDIF

      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!$OMP END CRITICAL

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_CH4

!------------------------------------------------------------------------------

      SUBROUTINE CH4_AVGTP
!
!******************************************************************************
!  Subroutine CH4_AVGTP gets the 24-h average surface pressure and temperature
!  needed for the CH4 simulation. (jsw, bnd, bmy, 1/16/01, 7/20/04)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry and
!        placed into module "global_ch4_mod.f" by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_AVGTP is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Removed duplicate definition for NTDT, NMIN (bmy, 11/15/01)
!  (4 ) Removed PS from argument list.  Now use P(I,J)+PTOP instead of
!        PS, this ensures that we have consistency between P and AD.
!        (bmy, 4/11/02)
!  (5 ) Removed obsolete code (bmy, 6/27/02)
!  (6 ) Now uses GET_PCENTER from "pressure_mod.f" to return the pressure
!        at the midpoint of the box (I,J,L).  Also added parallel DO-loops.
!        Updated comments. (dsa, bdf, bmy, 8/21/02)
!  (7 ) Now reference T from "dao_mod.f".  Now reference GEOS_CHEM_STOP from
!        "error_mod.f" (bmy, 10/15/02)
!  (8 ) Removed NTDT, NMIN from the arg list.  Now uses functions GET_TS_DYN,
!        GET_TS_CHEM, and GET_ELAPSED_MIN from "time_mod.f" (bmy, 3/27/03)
!  (9 ) Remove reference to CMN, it's not needed (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : T
      USE ERROR_MOD,    ONLY : GEOS_CHEM_STOP
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TIME_MOD,     ONLY : GET_TS_DYN, GET_TS_CHEM, GET_ELAPSED_MIN

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER             :: NTDT, NMIN
      INTEGER             :: I, J, L, NTIMES, MNDT, K, M, N
      INTEGER             :: NNEW, NNCOUNT 
      REAL*8              :: Ptemp(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CH4_AVGTP begins here!
      !=================================================================

      ! Get quantities from "time_mod.f"
      NTDT   = GET_TS_DYN() * 60
      NMIN   = GET_ELAPSED_MIN()
      MNDT   = NTDT  / 60 
      NTIMES = GET_TS_CHEM() / MNDT

      ! NTIMES is the number of dynamic timesteps in a chem timestep
      IF ( NMIN <= GET_TS_CHEM() ) NTIMES = NTIMES + 1

      ! At the start of the run...
      IF ( NMIN == 0 ) THEN

         ! Initialize NNEW
	 NNEW = 0

         ! Error check --  NCHEM has to be 1440 min
         IF ( GET_TS_CHEM() /= 1440 ) THEN   
            WRITE(*,*) ' '
            WRITE(*,*) 'CO-OH parameterization option (i.e., NSRCX=5)!' 
            WRITE(*,*) 'Use a chemistry time step = 24 hours'
            WRITE(*,*) '(i.e., NCHEM=1440 min.)'
            WRITE(*,*) ' '
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Error check -- need timestep to be divisible by 1440
         IF ( mod( GET_TS_CHEM(), MNDT ) /= 0 ) THEN   
            WRITE(*,*) ' '
            WRITE(*,*) 'CO-OH parameterization option (i.e., NSRCX=5)!'
            WRITE(*,*) 'The chemistry time step (i.e., 24 hours) is'
            WRITE(*,*) 'not evenly divisible by the meteorological'
            WRITE(*,*) 'data read-in time step (i.e., 6 hours).  This'
            WRITE(*,*) 'will mess up SR avgtp which calculates a 24-'
            WRITE(*,*) 'hour average temperature and pressure to be'
            WRITE(*,*) 'used by SR getinfo.'
            WRITE(*,*) ' '
            CALL GEOS_CHEM_STOP
         ENDIF

         ! If NCHEM < NTDT then stop program.
         IF ( GET_TS_CHEM() < MNDT ) THEN   
            WRITE(*,*) ' '
            WRITE(*,*) 'When using the CO-OH parameterization'
            WRITE(*,*) 'option (i.e., NSRCX=5), take a 24-hour'
            WRITE(*,*) 'time step (i.e., NCHEM=1440 min.) because'
            WRITE(*,*) 'the OH parameterization produces a 24-hour'
            WRITE(*,*) 'average [OH]'
            WRITE(*,*) ' '
            CALL GEOS_CHEM_STOP
         ENDIF
      ENDIF

      !=================================================================
      ! If a new 24-hr period, set Pavg = 0, and reset NNEW, NCOUNT
      !=================================================================
      IF ( NNEW == 0 ) THEN 
         Pavg(:,:,:) = 0d0
         Tavg(:,:,:) = 0d0
	 NNEW        = 1
	 NNCOUNT     = 0
      ENDIF

      !=================================================================
      ! Archive quantities
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PTEMP )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
                  
         ! Archive pressure
         Pavg(I,J,L) = Pavg(I,J,L) + GET_PCENTER(I,J,L)

         ! Archive temperature
         Tavg(I,J,L) = Tavg(I,J,L) + T(I,J,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !================================================================
      ! Keep track to see if at end of NCHEM time step.
      ! If so, divide PAVG & TAVG by the number of times archived.
      !=================================================================
      NNCOUNT = NNCOUNT + 1

      IF ( NNCOUNT == NTIMES ) THEN
         Pavg(:,:,1:LLPAR) = Pavg(:,:,1:LLPAR) / DBLE( NTIMES )
         Tavg(:,:,1:LLPAR) = Tavg(:,:,1:LLPAR) / DBLE( NTIMES )
         NNEW              = 0
      ENDIF

      ! Return to calling program
      END SUBROUTINE CH4_AVGTP

!------------------------------------------------------------------------------

      SUBROUTINE EMISSCH4
!
!******************************************************************************
!  Subroutine EMISSCH4 places emissions of CH4 [kg] into the STT array.
!  (jsw, bnd, bey, bmy, 1/16/01, 7/20/04)
!
!  I might want to eventually adjust my swamps (WS) source upward
!  because it is currently 39.1 Tg/yr whereas in Fung et al. [1991] it is
!  80 Tg/yr in their preferred scenario.  The source strength for many of
!  the other sources is also different in the emissions files from their
!  preferred scenario. (jsw, 1/16/01)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) EMISSCH4 is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) GLOBASEAEMIS, GLOBSEAEMIS are diagnostics by jsw.
!  (4 ) Do not multiply CO emissions by 1.28 anymore (jsw, bmy, 2/12/01)
!  (5 ) Renamed input files to CH4_monthly.geos.{RES} and 
!        CH4_aseasonal.geos.{RES}. (bmy, 2/12/01)
!  (6 ) Add reference to "CMN_SETUP" for the DATA_DIR variable (bmy, 2/13/01)
!  (7 ) Removed references to "biofuel_mod.f" and "biomass_mod.f"; these
!        weren't necessary (bmy, 3/20/01)
!  (8 ) Now reference IU_FILE and IOERROR from "file_mod.f".  Now use IU_FILE
!        instead of IUNIT as the file unit #. (bmy, 6/27/02)
!  (9 ) Now reference BXHEIGHT and SUNCOS from "dao_mod.f".  Remove reference 
!        to header file "comtrid.h" -- it's not used.  Make FIRSTEMISS a local
!        SAVEd variable.  Also use MONTH from "CMN" instead of the variable
!        LMN. (bmy, 11/15/02)
!  (10) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f".
!        Now use function GET_MONTH and GET_TS_EMIS from "time_mod.f". 
!        Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        I0 and J0 are now local variables. (bmy, 3/27/03)
!  (11) Now reference STT from "tracer_mod.f".  Now reference DATA_DIR from
!        "directory_mod.f". (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DAO_MOD,       ONLY : BXHEIGHT,     SUNCOS
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE,      IOERROR
      USE GRID_MOD,      ONLY : GET_AREA_CM2, GET_XOFFSET, GET_YOFFSET 
      USE TIME_MOD,      ONLY : GET_MONTH,    GET_TS_EMIS
      USE TRACER_MOD,    ONLY : STT

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches

      ! Local Variables
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE. 
      INTEGER                :: I, IREF, IJLOOP, J, JREF, L, N, K, LMN
      INTEGER                :: IOS, SCALEYEAR, I0, J0
      INTEGER, SAVE          :: LASTYEAR
      REAL*8                 :: DTSRCE,       AREA_CM2
      REAL*8                 :: GLOBASEAEMIS, GLOBSEAEMIS  
      REAL*8                 :: EMIS1(IGLOB,JGLOB)  
      REAL*8                 :: EMIS2(IGLOB,JGLOB)
      REAL*8                 :: EMIS3(IGLOB,JGLOB)    
      REAL*8                 :: EMIS4(IGLOB,JGLOB)
      REAL*8                 :: EMIS5(IGLOB,JGLOB)    
      REAL*8                 :: EMIS6(IGLOB,JGLOB)
      REAL*8                 :: EMIS7(IGLOB,JGLOB)
      REAL*8                 :: EMIS8(IGLOB,JGLOB)    
      REAL*8                 :: EMIS9(IGLOB,JGLOB)    
      REAL*8                 :: EMIS10(IGLOB,JGLOB,12)
      REAL*8                 :: EMIS11(IGLOB,JGLOB,12) 
      REAL*8                 :: EMIS12(IGLOB,JGLOB,12)
      REAL*8                 :: EMIS13(IGLOB,JGLOB,12) 
      REAL*8                 :: EMIS14(IGLOB,JGLOB,12)
      REAL*8                 :: EMIS15(IGLOB,JGLOB,12)
      CHARACTER(LEN=255)     :: FILENAME

      ! External functions
      REAL*8, EXTERNAL       :: BOXVL

      ! LMN is the current month
      LMN = GET_MONTH()

      ! Get nested-grid offsets
      I0  = GET_XOFFSET()
      J0  = GET_YOFFSET()

      !=================================================================
      ! jsw: I might want to use the following scaling code as a 
      ! template for scaling my yearly CH4 emissions for different 
      ! sources.  Question: Do these pre-existing scale factors contain 
      ! different values for different geographic regions?  If so, I 
      ! might want to eventually use the same scale factors for my 
      ! fossil fuel CH4 emissions, although I'd need to modify to get 
      ! solid fuel scale factors (coal) and gaseous fuel (natural gas).
      !
      ! Read aseasonal CH4 emissions 
      !=================================================================
!---------------------------------------------------------------------------
! Prior to 8/16/05:
!      FILENAME = TRIM( DATA_DIR )      // 
!     &           'CH4_aseasonal.geos.' // GET_RES_EXT()
!---------------------------------------------------------------------------
      FILENAME = TRIM( DATA_DIR )            // 
     &           'CH4_200202/CH4_aseasonal.' // GET_NAME_EXT_2D() // 
     &           '.'                         // GET_RES_EXT()

      OPEN( IU_FILE, FILE=TRIM( FILENAME ), IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissch4:1' )

      READ( IU_FILE, '(6(E14.6))', IOSTAT=IOS ) 
     &   EMIS1, EMIS2, EMIS3, EMIS4, EMIS5, EMIS6, EMIS7, EMIS8, EMIS9
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissch4:2' )

      CLOSE( IU_FILE )

      !=================================================================
      ! Read monthly-varying CH4 emissions
      !=================================================================
!---------------------------------------------------------------------------
! Prior to 8/16/05:
!      FILENAME = TRIM( DATA_DIR )    // 
!     &           'CH4_monthly.geos.' // GET_RES_EXT()
!---------------------------------------------------------------------------
      FILENAME = TRIM( DATA_DIR )          // 
     &           'CH4_200202/CH4_monthly.' // GET_NAME_EXT_2D() // 
     &           '.'                       // GET_RES_EXT()

      OPEN( IU_FILE, FILE=TRIM( FILENAME ), IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissch4:3' )      

      READ( IU_FILE, '(6(E14.6))', IOSTAT=IOS )
     &     EMIS10,EMIS11,EMIS12,EMIS13,EMIS14,EMIS15
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissch4:4' )      

      CLOSE( IU_FILE )

      ! Chemistry timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      !=================================================================  
      ! In the following, we modify the STT with emissions rates.
      !
      ! Aseasonal CH4 emissions. 
      ! EMIS1 through EMIS9 (I,J) in molec CH4/cm^2/s
      !
      ! NOTES: 
      ! (1) Test to see if low emissions are causing low [CH4].
      !      Multiply emissions by a factor of 1.28, excluding the
      !      soil sink (jsw, 11/18/00) 
      !=================================================================  

      !### Debug
      !print*,'Inside emissco, before aseasonal emissions'
      GLOBASEAEMIS = 0d0  
      GLOBSEAEMIS  = 0d0

      !jsw  J0 and I0 are global variables, both set = 0.
      DO J = 1, JJPAR
         JREF = J + J0     

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )
  
      DO I = 1, IIPAR
         IREF = I + I0
         
         !### Debug
         IF ( JREF == 300 ) THEN
            WRITE(6,*) 'Inside emissco   STT(IREF,30,1,1)='
            WRITE(6,*) STT(IREF,JREF,1,1) 
            
            WRITE(6,*) '(EMIS1+) XNUMOL_CH4*DTSRCE*DXYP(JREF)*10000'
            WRITE(6,*) ( EMIS1(IREF,JREF) + EMIS2(IREF,JREF) +
     &                   EMIS3(IREF,JREF) + EMIS4(IREF,JREF) +
     &                   EMIS5(IREF,JREF) + EMIS6(IREF,JREF) + 
     &                   EMIS7(IREF,JREF) + EMIS8(IREF,JREF) + 
     &                   EMIS9(IREF,JREF) ) / 
     &                 XNUMOL_CH4 * DTSRCE * AREA_CM2


            WRITE(6,*) '(EMIS1+):'
            WRITE(6,*) ( EMIS1(IREF,JREF) + EMIS2(IREF,JREF) + 
     &                   EMIS3(IREF,JREF) + EMIS4(IREF,JREF) +
     &                   EMIS5(IREF,JREF) + EMIS6(IREF,JREF) + 
     &                   EMIS7(IREF,JREF) + EMIS8(IREF,JREF) + 
     &                   EMIS9(IREF,JREF) )
         ENDIF

         ! Modify emisison rates
         !jsw I think 10000 (or 1d4) converts m^2 to cm^2
         STT(IREF,JREF,1,1) = STT(IREF,JREF,1,1) + 
     &        ( ( EMIS1(IREF,JREF) + EMIS2(IREF,JREF) + 
     &            EMIS3(IREF,JREF) + EMIS4(IREF,JREF) +
     &            EMIS5(IREF,JREF) + EMIS6(IREF,JREF) +
     &            EMIS7(IREF,JREF) + EMIS8(IREF,JREF) ) +
     &          EMIS9(IREF,JREF) ) / 
     &        XNUMOL_CH4 * DTSRCE * AREA_CM2

         !### Debug
         GLOBASEAEMIS = GLOBASEAEMIS + 
     &                  ( EMIS1(IREF,JREF) + EMIS2(IREF,JREF) + 
     &                    EMIS3(IREF,JREF) + EMIS4(IREF,JREF) +
     &                    EMIS5(IREF,JREF) + EMIS6(IREF,JREF) +
     &                    EMIS7(IREF,JREF) + EMIS8(IREF,JREF) + 
     &                    EMIS9(IREF,JREF) ) / 
     &                  XNUMOL_CH4 * DTSRCE * AREA_CM2
      ENDDO
      ENDDO

      !-----------------------------------------------------------------
      !### Debug
      !WRITE( 6, * ) 'Inside emissco  Global asea CH4 emis for this'
      !WRITE( 6, * ) 'time step in Tg ='
      !WRITE( 6, * ) globaseaemis/1.e9
      !WRITE( 6, * ) after aseasonal CH4 emissions'
      !-----------------------------------------------------------------

      !=================================================================  
      ! In the following, we modify the STT with emissions rates.
      !
      ! Monthly CH4 emissions.
      ! EMIS10 through EMIS15 in molec CH4/cm^2/s
      !=================================================================  

      !### Debug
      !WRITE (66666,*) 'LMN (month?) = ',LMN

      DO J = 1, JJPAR
         JREF = J + J0 

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

      DO I = 1, IIPAR
         IREF = I + I0

         !### Debug
         IF ( JREF == 300 ) THEN
            WRITE(6,*) 'Inside emissco SEAS   STT(IREF,30,1,1)='
            WRITE(6,*) STT(IREF,JREF,1,1) 

            !jsw I think 10000 (or 1d4) converts m^2 to cm^2
            WRITE(6,*) '(EMIS10+) XNUMOL_CH4*DTSRCE*DXYP(JREF)*10000'
            WRITE(6,*) ( EMIS10(IREF,JREF,LMN) + EMIS11(IREF,JREF,LMN) + 
     &                   EMIS12(IREF,JREF,LMN) + EMIS13(IREF,JREF,LMN) +
     &                   EMIS14(IREF,JREF,LMN) + EMIS15(IREF,JREF,LMN) )
     &                 / XNUMOL_CH4 * DTSRCE * AREA_CM2

            WRITE(*,*) '(EMIS10+):'
            WRITE(*,*) ( EMIS10(IREF,JREF,LMN) + EMIS11(IREF,JREF,LMN) + 
     &                   EMIS12(IREF,JREF,LMN) + EMIS13(IREF,JREF,LMN) +
     &                   EMIS14(IREF,JREF,LMN) + EMIS15(IREF,JREF,LMN) )
         ENDIF

         ! Modify STT with emission rates of CH4
         ! Multiply by 1.28 to test if low emissions are causing low [CH4]
         ! (jsw, 11/18/00)
         STT(IREF,JREF,1,1) = STT(IREF,JREF,1,1) + 
     &        ( EMIS10(IREF,JREF,LMN) + EMIS12(IREF,JREF,LMN) + 
     &          EMIS13(IREF,JREF,LMN) + EMIS14(IREF,JREF,LMN) + 
     &          EMIS15(IREF,JREF,LMN) ) 
     &        / XNUMOL_CH4 * DTSRCE * AREA_CM2

         !### Debug
         GLOBSEAEMIS = GLOBSEAEMIS + 
     &        ( EMIS10(IREF,JREF,LMN) + EMIS12(IREF,JREF,LMN) + 
     &          EMIS13(IREF,JREF,LMN) + EMIS14(IREF,JREF,LMN) + 
     &          EMIS15(IREF,JREF,LMN) )
     *        / XNUMOL_CH4 * DTSRCE * AREA_CM2

      ENDDO
      ENDDO

      !-----------------------------------------------------------------
      !### Debug (jsw)
      !IF (LMN.EQ.10) THEN
      !   write(*,*) 'Inside emissco  Global sea CH4 emis for this'
      !   write(*,*) 'time step in Tg ='
      !   write(*,*) globseaemis/1.e9
      !ENDIF
      !WRITE(6,*) 'Inside emissco, after seasonal CH4 emissions'
      !-----------------------------------------------------------------

      !=================================================================  
      ! Sum up CH4 budgets
      !
      ! The 1.28 in the following lines is for increasing the total 
      ! source strength from 453 Tg/yr to 581 Tg/yr. (jsw, 11/18/00)
      !
      ! DXYP(JREF) * 1d4 is grid box surface area in cm2
      !
      ! NOTE: Don't multiply by 1.28 anymore (jsw, bmy, 2/12/01)
      !=================================================================  
      DO J=1,JJPAR
         JREF = J + J0

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

      DO I = 1,IIPAR

         TCH4(I,J,1,5) = TCH4(I,J,1,5) +
     &        ( EMIS1(I,J) * AREA_CM2 * DTSRCE ) !* 1.28d0 )
            
         TCH4(I,J,1,6) = TCH4(I,J,1,6)+
     &        ( EMIS4(I,J) * AREA_CM2 * DTSRCE ) !* 1.28d0 )
            
         TCH4(I,J,1,7) = TCH4(I,J,1,7)+
     &        ( EMIS5(I,J) * AREA_CM2 * DTSRCE ) !* 1.28d0 )

         TCH4(I,J,1,8) = TCH4(I,J,1,8)+
     &        ( EMIS12(I,J,LMN) * AREA_CM2 * DTSRCE ) !* 1.28d0 )
         
         TCH4(I,J,1,9) = TCH4(I,J,1,9)+
     *        ( EMIS10(I,J,LMN) * AREA_CM2 * DTSRCE ) !* 1.28d0 )

         TCH4(I,J,1,10) = TCH4(I,J,1,10)+
     *        ( EMIS9(I,J) * AREA_CM2 * DTSRCE )

         TCH4(I,J,1,4) = TCH4(I,J,1,4) + 
     &        ( EMIS1(I,J)      + EMIS2(I,J)      + 
     &          EMIS3(I,J)      + EMIS4(I,J)      + 
     &          EMIS5(I,J)      + EMIS6(I,J)      + 
     &          EMIS7(I,J)      + EMIS8(I,J)      + 
     &          EMIS10(I,J,LMN) + EMIS12(I,J,LMN) + 
     &          EMIS13(I,J,LMN) + EMIS14(I,J,LMN) + 
     &          EMIS15(I,J,LMN) ) * 
     &        AREA_CM2 * DTSRCE !* 1.28d0 

      ENDDO
      ENDDO

      ! Set FIRSTEMISS to FALSE
      FIRSTEMISS = .FALSE.

      ! Return to calling program
      END SUBROUTINE EMISSCH4

!------------------------------------------------------------------------------

      SUBROUTINE CHEMCH4
!
!******************************************************************************
!  Subroutine CHEMCH4 computes the chemical loss of CH4 (sources - sinks).
!  (jsw, bnd, bmy, 6/8/00, 7/20/04)
!
!  CH4 SOURCES
!  ============================================================================
!  (1 ) Oxidation of methane, isoprene and monoterpenes (SRCO_fromHCs).
!  (2 ) Direct emissions of CO from fossil fuel combustion, biomass 
!        burning and wood (for fuel) burning (SR SETEMIS).
!  (3 ) Emissions.
!
!  CH4 SINKS:
!  ============================================================================
!  (1 ) Removal of CO by OH (SR OHparam & CO_decay).
!  (2 ) CO uptake by soils (neglected).
!  (3 ) Transport of CO to stratosphere from troposphere 
!        (in dynamical subroutines).
!  (4 ) Removal by OH (Clarissa's OH--climatol_OH.f and CO_decay.f)
!  (5 ) Transport of CH4 between troposphere and stratosphere, and 
!        destruction in strat (CH4_strat.f).
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CHEMCH4 is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/12/01)
!  (4 ) LD43 is already declared in CMN_DIAG; don't redefine it (bmy, 11/15/01)
!  (5 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (6 ) Now reference AD from "dao_mod.f".  Now reference GEOS_CHEM_STOP from
!        "error_mod.f"  Now make FIRSTCHEM a local SAVEd variable.  Now 
!        reference ALBD from "dao_mod.f".  Now use MONTH and JDATE from "CMN"
!        instead of LMN and LDY. (bmy, 11/15/02)
!  (7 ) Remove NYMDb, NYMDe from the arg list.  Now use functions GET_MONTH,
!        GET_NYMDb, GET_NYMDe, GET_MONTH, GET_DAY from the new "time_mod.f"
!        (bmy, 3/27/03) 
!  (8 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DAO_MOD,       ONLY : AD, ALBD
      USE DIAG_MOD,      ONLY : AD43
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,      ONLY : GET_DAY, GET_MONTH, GET_NYMDb, GET_NYMDe

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! LPAUSE
#     include "CMN_DIAG"     ! ND43, AD43

      ! Local variables
      LOGICAL                :: FIRSTCHEM = .TRUE.
      INTEGER                :: I, J, L, K, M, N
      INTEGER                :: IJ, JJ, NPART, III, JJJ
      INTEGER                :: NOHDO
      INTEGER                :: NCLIMATOLOGY, NCLIMATOLOGY2
      INTEGER, SAVE          :: NBUDGET, NTALDT

      CHARACTER(LEN=255)     :: FILENAME

      ! Number of days per month
      INTEGER                :: NODAYS(12) = (/ 31, 28, 31, 30, 
     &                                          31, 30, 31, 31, 
     &                                          30, 31, 30, 31 /)

      ! External functions 
      REAL*8 , EXTERNAL      :: BOXVL
  
      ! Weight of air (taken from "comode.h") 
      REAL*8, PARAMETER      :: WTAIR = 28.966d0

      !=================================================================
      ! CHEMCH4 begins here!
      !
      ! Call subroutine READER which opens and reads from data files   
      ! m.dat AND tracer.dat                                           
      !=================================================================
      WRITE( 6, '(a)' ) '--- ENTERING CHEMCH4! ---'

      CALL READER( FIRSTCHEM )

      !=================================================================
      ! (0) Calculate each box's air density [molec/cm3]
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         BAIRDENS(I,J,L) = AD(I,J,L) * 1000d0   / BOXVL(I,J,L) * 
     &                                 6.023D23 / WTAIR
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! (1) If the first time step ...
      !=================================================================
      IF ( FIRSTCHEM ) THEN

         ! Counter for total number of timesteps per month for CO budget.
         NTALDT = 1

         ! Flag for SR CO_budget
         NBUDGET = 1 

         ! Zero CO Production array
         COPROD(:,:,:) = 0d0

         ! Read zonally-averaged CO production [v/v/s]
         CALL READ_COPROD

         ! Added following line to increase strat. sink strength
         ! Hmm, the values I printed out above for COprod are very small, 
         ! all less than 1e-15. (jsw)
         COprod = COprod * 3d0
      ENDIF                   
      
      ! Increment counter of timesteps 
      NTALDT = NTALDT + 1

      !=================================================================
      ! (2) Calculate the production and destruction of CO from 
      !     gas-phase chemistry only.       
      !
      ! Concerning O3, there are 3 options:
      !   A) The OH parameterization is calculated using GEOS monthly 
      !      means (NCLIMATOLOGY=0) for the independent variable O3.  
      !      The O3 column above independent variable is determined 
      !      using jal's O3  climatologies for both the tropospheric 
      !      and stratospheric portions of the O3 column 
      !      (NCLIMATOLOGY2=1).  
      !
      !   B) The O3 variable is determined from jal's O3 climatolgies 
      !      (tropospheric portion) and the o3 column above variable 
      !      is determined from jal's O3 climatolgies (NCLIMATOLOGY=1 & 
      !      NCLIMATOLOGY2=1).
      !
      !   C) The O3 variable is determined by GEOS monthly means 
      !      (NCLIMATOLOGY=0) and the o3 column above variable is 
      !      determined from fastj climatologies for both the 
      !      tropospheric and stratospheric portion of the O3 column 
      !      (NCLIMATOLOGY=0 & NCLIMATOLOGY2=0).
      !=================================================================
      NCLIMATOLOGY  = 0 !jsw changed this value from 0
      NCLIMATOLOGY2 = 1

      ! Error check
      IF( NCLIMATOLOGY == 1 .AND. NCLIMATOLOGY2 == 0 ) THEN
         PRINT*,'Stopped in SR CHEMCO.'
         PRINT*,'This combination of options is not allowed!'
         PRINT*,'Reset NCLIMATOLOGY and/or NCLIMATOLOGY2.'
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! (3) get parameterized OH fields or monthly mean fields.
      !
      ! Variables of note:
      ! ----------------------------------------------------------------
      ! (1) BOH = storage array for OH fields.
      !
      ! (2) NOHDO = switch
      !       = 0 : Use OH field from full chemistry monthly avg field.
      !       = 1 : Get parameterized OH field.
      !       = 2 : Get Clarissa's climatological OH (jsw)
      !
      ! (3) LPAUSE =  the vertical level of the tropopause.  Above this
      !     level, no [OH] is calculated.  The user can feed this
      !     SR a high value for LPAUSE which effectively turns this 
      !     option off (i.e., LPAUSE > MVRTBX). If the [OH] = -999 
      !     then the [OH] was not calculated.
      !=================================================================
      BOH(:,:,:) = 0d0

      ! Change value of NOHDO as listed above
      NOHDO = 2

      SELECT CASE ( NOHDO )

         ! NOHDO = 0: Get full chem monthly avg OH field
         CASE ( 0 )
            ! Comment out for now (bmy, 1/16/01)
            !BOH(:,:,:) = BBIJ(:,:,:,NFIELDS2)
            
         ! NOHDO = 1: Get parameterized OH field 
         CASE ( 1 )
            ! Comment out for now (bmy, 1/16/01)
            !CALL GETINFO( LMN, ALBD   )
            !CALL OHPARAM( BOH, LPAUSE )

            !### Debug
            !write(*,*) 'Inside chemco  BOH(1,J,1)=',(BOH(1,J,1),J=1,46)
            !write(*,*) 'BOH(1,23,L)=',(BOH(1,23,L),L=1,20)

         ! NOHDO = 2: Get Clarisa's climatological OH
         CASE ( 2 )
            CALL CLIMATOL_OH( GET_MONTH(), GET_DAY() )

         ! Error
         CASE DEFAULT
            WRITE( 6, '(a)' ) 'Invalid selection for NOHDO!'
            WRITE( 6, '(a)' ) 'Halting execution in CHEMCH4!'
            CALL GEOS_CHEM_STOP
            
      END SELECT

      !=================================================================
      ! (3.1) ND43 diagnostics...save [OH] in molecules/cm3
      !=================================================================
      IF ( ND43 > 0 ) THEN
         DO L = 1, LD43
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( L < LPAUSE(I,J) ) THEN
               AD43(I,J,L,1) = AD43(I,J,L,1) + BOH(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! (4) Diagnostic for Methyl Chloroform (ND23)
      !=================================================================
      CALL CH4_OHSAVE

      !=================================================================
      ! (5) calculate rate of decay of CH4 by OH oxidation.
      !=================================================================
      CALL CH4_DECAY

      !=================================================================
      ! (6) do CH4 chemistry in layers above tropopause.
      !=================================================================
      CALL CH4_STRAT

      !=================================================================
      ! (7) write budget (i.e., monthly average fields).
      !
      ! Check to make sure the start and end times are on the
      ! first of a month.  If not the SR CO_budget will not
      ! work properly!
      !=================================================================
      NPART = GET_NYMDb() / 100 

      IF ( ( GET_NYMDb() - NPART*100 ) /= 1 ) THEN
         print*,'Start date not equal to 1st of month!!!'
         print*,'  Therefore, SR CO_budget will not work!!!'
         CALL GEOS_CHEM_STOP
      ENDIF

      NPART = GET_NYMDe() /100 

      IF ( ( GET_NYMDe() - NPART*100 ) /= 1 ) THEN      
         print*,'End date not equal to 1st of month!!!'
         print*,'  Therefore, SR CO_budget will not work!!!'
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Call CH4_budget
      IF ( NTALDT == NODAYS( GET_MONTH() ) ) THEN
         CALL CH4_BUDGET( FIRSTCHEM, NBUDGET )

         NBUDGET = 0
	 NTALDT  = 0
      ENDIF
      print*, 'after ch4_budget'
      call flush(6)

      ! Set FIRSTCHEM to FALSE
      FIRSTCHEM = .FALSE.

      ! Return to calling program
      END SUBROUTINE CHEMCH4

!------------------------------------------------------------------------------

      SUBROUTINE READ_COPROD
!
!*****************************************************************************
!  Subroutine READ_COPROD reads production and destruction rates for CO in 
!  the stratosphere. (bnd, bmy, 1/17/01, 8/16/05)
!
!  Module Variables:
!  ===========================================================================
!  (1) COPROD (REAL*8) : Array containing P(CO) for all 12 months [v/v/s]
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) READ_COPROD is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) ARRAY needs to be dimensioned (1,JGLOB,LGLOB) (bmy, 9/26/01)
!  (4 ) Remove obsolete code from 9/01 (bmy, 10/24/01)
!  (5 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Now reads data for both GEOS and GCAP grids (bmy, 8/16/05)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_ZONAL

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER            :: J, L, M
      REAL*4             :: ARRAY(1,JGLOB,LGLOB) 
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME
 
      !=================================================================
      ! READ_COPROD begins here!
      ! 
      ! Read P(CO) for all 12 months
      !=================================================================
      DO M = 1, 12 

         ! TAU value at the start of month M -- Use "generic" year 1985
         XTAU = GET_TAU0( M, 1, 1985 )

         ! Construct filename
!---------------------------------------------------------------------------
! Prior to 8/16/05:
!         FILENAME = TRIM( DATA_DIR ) // 'COprod.'  //
!     &              GET_NAME_EXT()   // '.'        // GET_RES_EXT()
!---------------------------------------------------------------------------
         FILENAME = TRIM( DATA_DIR )         // 
     &              'pco_lco_200203/COprod.' // GET_NAME_EXT_2D() //
     &              '.'                      // GET_RES_EXT()

         ! Read P(CO) in units of [v/v/s]
         CALL READ_BPCH2( FILENAME, 'PORL-L=$', 9,     
     &                    XTAU,      1,         JGLOB,     
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )
         
         ! Copy REAL*4 to REAL*8 data, and resize from (JGLOB,LGLOB) 
         ! to (JJPAR,LLPAR) -- vertically regrid if necessary
         CALL TRANSFER_ZONAL( ARRAY(1,:,:), COPROD(:,:,M) )
      ENDDO

      ! Return to calling program
      END SUBROUTINE READ_COPROD

!------------------------------------------------------------------------------

      SUBROUTINE CLIMATOL_OH( LMN, LDY )
!
!******************************************************************************
!  Subroutine CLIMATOL_OH gets climatological OH fields for use in 
!  calculating CH4 loss. (jsw, bnd, bmy, 1/16/00, 7/20/04)
!
!  OH from Clarissa Spivakovsky et al's paper accepted in 1999.
!  Written by James Wang, 7/24/00, based on OHparam.f written by Bryan Duncan.
!
!  Arguments as Input:
!  ===========================================================================
!  (1) LTPAUSE (INTEGER) : Array of tropopause heights [levels]
!  (2) LMN     (INTEGER) : Current month number (1-12)
!  (3) LDY     (INTEGER) : Current day number (1-31)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CLIMATOL_OH is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now reference IU_FILE and IOERROR from "file_mod.f".  Now use IU_FILE
!        instead of IUNIT as the file unit #. (bmy, 6/27/02) 
!  (4 ) Now use function GET_YMID of "grid_mod.f" to compute grid box
!        latitude (bmy, 2/3/03)
!  (5 ) Now references DATA_DIR from "directory_mod.f".  Also call routine
!        GET_SEASON from "time_mod.f". (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE GRID_MOD,      ONLY : GET_YMID
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE TIME_MOD,      ONLY : GET_SEASON

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPAUSE

      ! Arguments
      INTEGER, INTENT(IN) :: LMN, LDY

      ! Local variables
      INTEGER             :: I, J, L, IOS, SEASON
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! Read in Clarisa's climatological OH data from "avgOH.cms"
      !=================================================================
      !-----------------------------------------------------------------
      ! Prior to 8/16/05:
      !FILENAME = TRIM( DATA_DIR ) // 'avgOH.cms' 
      !-----------------------------------------------------------------
      FILENAME = TRIM( DATA_DIR ) // 'CH4_200202/avgOH.cms' 

      OPEN( IU_FILE, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'interpoh:1' )

      READ( IU_FILE, '(//)', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'interpoh:2' )      

      DO I = 1, NSEAS
         READ( IU_FILE, '(////)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'interpoh:3' )

         DO J = 1, NCMSLATS
            READ( IU_FILE, '(5x,7(F6.2))', IOSTAT=IOS ) 
     &          ( AVGOH(I,J,L), L=1,NCMSALTS )
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'interpoh:4' )
         ENDDO
      ENDDO

      CLOSE( IU_FILE )

      !=================================================================
      ! Interpolate Clarisa's OH to GEOS-CHEM grid
      !=================================================================

      ! Compute the proper season for indexing the OH file
      SEASON = GET_SEASON()

      ! Initialize BOH array
      BOH(:,:,:) = 0d0

      DO L = 1, MAXVAL( LPAUSE )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only consider tropospheric boxes
         IF ( L < LPAUSE(I,J) ) THEN 
            CALL INTERPOH( SEASON, GET_YMID(J), PAVG(I,J,L), BOH(I,J,L))     

            ! Error check 
            !IF ( BOH(I,J,L) < 1d-20 ) THEN
            !   PRINT*, 'Zero OH at : ', I, J, L
            !   !PRINT*, 'PAVG       : ', PAVG(I,J,L)
            !   !PRINT*, 'OHSEASON   : ', OH_SEASON 
            !   !PRINT*, 'OHLAT      : ', OH_LAT_ALL(J)
            !   !PRINT*,''
            !   !PRINT*,'OH not calculated for this box!'
            !   !PRINT*,'  BOH(',I,',',J,',',L,')'
            !   !PRINT*,'Stopped in SR OHPARAM'
            !   !PRINT*,''
            !   !STOP
            !ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      !### Debug (jsw)
      !IF (LDY.EQ.2) THEN 
      !   PRINT*, 'James Wang OH climatology subroutine climatol_OH.f'
      !   PRINT*, 'LDY = ',LDY,'LMN = ',LMN, ' BOH(1,J,2) = '
      !   write(*,*) (BOH(1,J,2), J=1,46)
      !   write(*,*) '(BOH(1,J,15):',(BOH(1,J,15), J=1,46)
      !   write(*,*) 'BOH(1,1,L):',(BOH(1,1,L), L=1,20)
      !   write(*,*) 'BOH(1,2,L):',(BOH(1,2,L), L=1,20)
      !ENDIF
  
      ! Return to calling program
      END SUBROUTINE CLIMATOL_OH

!------------------------------------------------------------------------------
 
      SUBROUTINE INTERPOH( OH_SEASON, OH_LAT, OH_PRESS, OH_VALUE )
!
!******************************************************************************
!  Subroutine INTERPOH interpolates Clarisa Spivakovsky's climatological OH
!  values to GEOS-CHEM grid resolution (jsw, bnd, bmy, 1/16/01, 10/15/02).
!
!  In subdomains where OH is very low or negligible due to low
!  sunlight (e.g., high latitudes in winter), concentrations of OH are set 
!  to climatological mean values as a function of latitude, altitude and 
!  season. This SR picks the appropriate average OH field described in 
!  Spivakovsky et al., "Three-dimensional climatological distribution of 
!  tropospheric OH: update and evaluation", accepted to JGR, 1999.  
!  The fields are stored in array avgOH and read in SR CLIMATOL_OH.
!
!  Arguments as Input:
!  ============================================================================
!  (1) OH_SEASON (INTEGER) : Number of the season for CMS climatology
!  (2) OH_LAT    (REAL*8 ) : Latitude value for CMS climatology
!  (3) OH_PRESS  (REAL*8 ) : Pressure value for CMS climatology
!  
!  Arguments as Output:
!  ============================================================================
!  (4) OH_VALUE  (REAL*8 ) : OH value at nearest GEOS-CHEM grid box
!
!  Module Variables:
!  ============================================================================
!  AVGOH         (REAL*8 ) : Array containing the climatological OH values.
!  NCMSALTS      (INTEGER) : Number of altitude levels of climatology
!  NCMSLATS      (INTEGER) : Number of latitude bands of climatology
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) INTERPOH is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Remove references to CMN_SIZE, CMN (bmy, 11/15/01)
!  (4 ) Eliminate obsolete code from 11/01 (bmy, 2/27/02)
!  (5 ) Now reference GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN)  :: OH_SEASON
      REAL*8,  INTENT(IN)  :: OH_LAT, OH_PRESS
      REAL*8,  INTENT(OUT) :: OH_VALUE

      ! Local variables
      INTEGER              :: I, J, K

      !=================================================================
      ! INTERPOH begins here!
      !
      ! The following do loops take the climatological [OH] from the 
      ! box closest to the GEOS-CHEM box of interest.  CMSLATS and 
      ! CMSALTS are the latitudes and altitudes of Clarissa's OH 
      ! climatology, which are defined as module variables. 
      ! (jsw, bmy, 1/16/01)  
      !=================================================================

      !### Debug (jsw)
      !write(*,*) 'OH_LAT:',OH_LAT,'OH_PRESS:',OH_PRESS

      I = OH_SEASON

      DO J = 1, NCMSLATS
         DO K = NCMSALTS, 1, -1

            IF( OH_LAT >= CMSLATS(J) ) THEN         
               IF( CMSALTS(K) >= OH_PRESS ) THEN
                  OH_VALUE = AVGOH(I,J,K) * 1d5

                  !### Debug (jsw)
                  !write(*,*) 'CMSLATS(J):',CMSLATS(J)
                  !write(*,*) 'CMSALTS(K):',CMSALTS(K)
                  !write(*,*) 'OH_VALUE=',OH_VALUE
                  GOTO 2
               ENDIF
           
               IF ( K == 1 ) THEN
                  OH_VALUE = AVGOH(I,J,K) * 1d5

                  !### Debug(jsw)
                  !write(*,*) 'CMSLATS(J):',CMSLATS(J)
                  !write(*,*) 'CMSALTS(K):',CMSALTS(K)
                  !write(*,*) 'OH_VALUE=',OH_VALUE
                  GOTO 2
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      !=================================================================
      ! Error Check -- stop if no GEOS-CHEM box is found
      !=================================================================
      WRITE( 6, '(a)'    ) REPEAT( '=', 79 )
      WRITE( 6, '(a)'     ) 'STOPPED IN SR interpOH!'
      WRITE( 6, '(a)'     ) 'Point lies nowhere!'
      WRITE( 6, '(a,3i5)' ) 'CMS Box (I,J,K): ', I, J, K
      WRITE( 6, '(a)'    ) REPEAT( '=', 79 )
      CALL GEOS_CHEM_STOP

 2    CONTINUE

      ! Return to calling program
      END SUBROUTINE INTERPOH

!------------------------------------------------------------------------------

      SUBROUTINE CH4_DECAY
!
!******************************************************************************
!  Subroutine CH4_DECAY calculates the decay rate of CH4 by OH.  OH is the 
!  only sink for CH4 considered here. (jsw, bnd, bmy, 1/16/01, 7/20/04)
!
!  The annual mean tropopause is stored in the LPAUSE array 
!  (from header file "CMN").  LPAUSE is defined such that: 
!
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!
!  We now use LPAUSE instead of NSKIPL to denote the strat/trop boundary. 
!  (bmy, 4/18/00)  
!
!  Monthly loss of CH4 is summed in TCO(3)
!     TCH4(3)  = CH4 sink by OH
!
!  Module Variables:
!  ============================================================================
!  (1) BOH        (REAL*8) : Array holding global OH concentrations
!  (2) XNUMOL_CH4 (REAL*8) : Molec CH4 / kg CH4
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_DECAY is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,    ONLY : AIRVOL
      USE TIME_MOD,   ONLY : GET_TS_CHEM
      USE TRACER_MOD, ONLY : STT

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE

      ! Local variables
      INTEGER          :: I, J, L, M, N
      REAL*8           :: DT, GCH4, STT2GCH4, KRATE

      ! External variables
      REAL*8, EXTERNAL :: BOXVL

      !=================================================================
      ! CH4_DECAY begins here!
      !=================================================================

      ! Chemistry timestep in seconds
      DT = GET_TS_CHEM() * 60d0

      !=================================================================
      ! Compute decay of CO by OH in the troposphere
      !
      ! The decay for CH4 is calculated by:
      !    OH + CH4 -> CH3 + H2O 
      !    k = 2.45E-12 exp(-1775/T)
      !
      ! This is from JPL '97. JPL '00 does not revise '97 value. (jsw)
      !=================================================================
      DO L = 1, MAXVAL( LPAUSE )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only consider tropospheric boxes
         IF ( L < LPAUSE(I,J) ) THEN 

            !jsw  Is it all right that I'm using 
            !     24-hr avg temperature to calc. rate coeff.?
            KRATE = 2.45d-12 * EXP( -1775d0 / Tavg(I,J,L) )  

            ! Conversion from [kg/box] --> [molec/cm3]
            ! [kg CH4/box] * [box/cm3] * XNUMOL_CH4 [molec CH4/kg CH4]
            STT2GCH4 = 1d0 / AIRVOL(I,J,L) / 1d6 * XNUMOL_CH4 

            ! CH4 in [molec/cm3]
            GCH4 = STT(I,J,L,1) * STT2GCH4

            ! Sum loss in TCH4(3) (molecules/box)
            TCH4(I,J,L,3) = TCH4(I,J,L,3)+ 
     &           ( GCH4 * BOXVL(I,J,L) * KRATE * BOH(I,J,L) * DT )

            ! Calculate new CH4 value: [CH4]=[CH4](1-k[OH]*delt) 
            GCH4 = GCH4 * ( 1d0 - KRATE * BOH(I,J,L) * DT )

            ! Convert back from [molec/cm3] --> [kg/box]
            STT(I,J,L,1) = GCH4 / STT2GCH4

            !### Debug(jsw)
            !IF ( ( J == 23 ) .AND. ( I == 100 ) ) THEN
            !   WRITE(*,*) 'James within CO_decay within LPAUSE loop.'
            !   WRITE(*,*) 'I=1, J=23, L=', L
            !   WRITE(*,*) 'krate(I,J,L)*BOH(I,J,L) =', krate*BOH(I,J,L)
            !   WRITE(*,*) 'krate(I,J,L) =', krate
            !ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CH4_DECAY

!------------------------------------------------------------------------------

      SUBROUTINE CH4_OHSAVE
! 
!*****************************************************************************
!  Subroutine CH4_OHSAVE archives the CH3CCl3 lifetime from the OH
!  used in the CH4 simulation. (bnd, jsw, bmy, 1/16/01, 7/20/04)
!
!  The annual mean tropopause is stored in the LPAUSE array 
!  (from header file "CMN").  LPAUSE is defined such that: 
!
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!
!  Module Variables
!  ===========================================================================
!  (1) BOH (REAL*8) : Array containing global OH field 
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_OHSAVE is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now call DO_DIAG_OH_CH4 to pass OH diagnostic info to the
!        "diag_oh_mod.f" (bmy, 7/20/04)
!*****************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,    ONLY : DIAGCHLORO
      USE DIAG_OH_MOD, ONLY : DO_DIAG_OH_CH4

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPAUSE
#     include "CMN_DIAG"  ! ND23 switch

      ! Local variables
      INTEGER          :: I, J, L
      REAL*8           :: KCLO, LOSS, OHMASS, MASST

      ! External functions
      REAL*8, EXTERNAL :: BOXVL

      !=================================================================
      ! CH4_OHSAVE begins here!
      !
      ! Archive quantities for CH3CCl3 lifetime diagnostic
      !=================================================================
      IF ( ND23 > 0 ) THEN
         DO L = 1, MAXVAL( LPAUSE )
         DO J = 1, JJPAR 
         DO I = 1, IIPAR 
            
            ! Only process tropospheric boxes (bmy, 4/17/00)
            IF ( L < LPAUSE(I,J) ) THEN
               OHMASS = BOH(I,J,L) * BAIRDENS(I,J,L) * BOXVL(I,J,L)
               MASST  = BAIRDENS(I,J,L) * BOXVL(I,J,L)

               ! Pass OH mass & total mass to "diag_oh_mod.f"
               CALL DO_DIAG_OH_CH4( I, J, L, OHMASS, MASST )
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CH4_OHSAVE

!------------------------------------------------------------------------------

      SUBROUTINE CH4_STRAT
!
!*****************************************************************************
!  Subroutine CH4_STRAT calculates uses production rates for CH4 to 
!  calculate loss of CH4 in above the tropopause. 
!  (jsw, bnd, bmy, 1/16/01, 7/20/04)
!
!  Production (mixing ratio/sec) rate provided by Dylan Jones.  
!  Only production by CH4 + OH is considered.
!  
!  The annual mean tropopause is stored in the LPAUSE array 
!  (from header file "CMN").  LPAUSE is defined such that: 
! 
!  Levels           1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!         LPAUSE(I,J) <= L <= LLPAR           are stratospheric (bmy, 4/18/00)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_STRAT is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Removed LMN from the arg list and made it a local variable.  Now use 
!        functions GET_MONTH and GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!*****************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,    ONLY : AIRVOL
      USE TIME_MOD,   ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACER_MOD, ONLY : STT

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE

      ! Local variables
      INTEGER             :: I, J, L, LMN
      REAL*8              :: DT, GCH4, STT2GCH4

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL

      !=================================================================
      ! CH4_STRAT begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DT  = GET_TS_CHEM() * 60d0

      ! Current month
      LMN = GET_MONTH()

      !=================================================================
      ! Loop over stratospheric boxes only
      !=================================================================
      DO L = MINVAL( LPAUSE ), LLPAR 
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( L >= LPAUSE(I,J) ) THEN

            ! Conversion factor [kg/box] --> [molec/cm3]
            ! [kg/box] / [AIRVOL * 1e6 cm3] * [XNUMOL_CH4 molec/mole]
            STT2GCH4 = 1d0 / AIRVOL(I,J,L) / 1d6 * XNUMOL_CH4
  
            ! CH4 in [molec/cm3]
            GCH4 = STT(I,J,L,1) * STT2GCH4

            ! Sum loss in TCH4(3) [molec CH4/box] in the stratosphere
            ! [molec/cm3] * [v/v/s] * [s] * [cm3/box] = [molec CH4/box]
            TCH4(I,J,L,3) = TCH4(I,J,L,3) + 
     &                      ( BAIRDENS(I,J,L) * COPROD(J,L,LMN) *
     &                        DT              * BOXVL(I,J,L)    )

            ! Calculate new CH4 value [molec CH4/cm3] in the stratosphere
            ! [v/v/s] * [s] * [molec/cm3] = [molec CH4/cm3]            
            GCH4 = GCH4 - ( COPROD(J,L,LMN) * DT * BAIRDENS(I,J,L) )

            ! Convert back from [molec CH4/cm3] --> [kg/box] 
            STT(I,J,L,1) = GCH4 / STT2GCH4
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CH4_STRAT

!------------------------------------------------------------------------------

      SUBROUTINE CH4_BUDGET( FIRSTCHEM, NBUDGET )
!
!******************************************************************************
!  Subroutine CH4_BUDGET calculates the budget of CH4.  This SR only works 
!  for monthly averages, so be sure to start on the first of the month 
!  and run to another first of the month!!!  (jsw, bnd, bmy, 1/16/01, 7/20/04)
!
!  Store the sources/sinks of CH4 in TCH4 in total molecules
!           ( 1) = Initial burden
!           ( 2) = Final burden
!  SINKS
!           ( 3) = CH4 sink by OH
!  SOURCES
!           ( 4) = Total Source
!           ( 5) = Animals
!           ( 6) = Gas Leakage
!           ( 7) = Coal
!           ( 8) = Bogs 
!           ( 9) = Rice
!           (10) = Soil absorption
!           (11) = Interhemispheric Exchange (+ = northward)
!           (12) = Bogs
!
!  Arguments as Input:
!  ============================================================================
!  (1) FIRSTCHEM (LOGICAL) : Flag, =T for very first chemistry timestep
!
!  Arguments as Output:
!  ============================================================================
!  (3) NBUDGET   (INTEGER) : Switch for computing budget (1=on)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_BUDGET is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/13/01)
!  (4 ) Renamed XLABEL to LABEL so as not to conflict w/ "CMN"
!  (5 ) Now use functions GET_MONTH, GET_YEAR, GET_DIAGb, and GET_CT_DYN from 
!        "time_mod.f".  Removed LMN from the arg list and made it a local 
!        variable.  Use functions GET_XOFFSET and GET_YOFFSET from 
!        "grid_mod.f".  (bmy, 3/27/03)
!  (6 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,   ONLY : GET_MONTH, GET_YEAR, GET_DIAGb, GET_CT_DYN
      USE TRACER_MOD, ONLY : STT

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE
 
      ! Arguments
      LOGICAL, INTENT(IN)    :: FIRSTCHEM
      INTEGER, INTENT(INOUT) :: NBUDGET

      ! Local variables
      INTEGER                :: I, J, K, L, M, NERROR, UD, LMN

      REAL*8                 :: STTCONV, TGS, SCALEDYN
      REAL*8                 :: NTP, NTQ, NTP2, NTQ2 
      REAL*8                 :: SOURCES, SINKS

      CHARACTER(LEN=16)      :: MERGE
      CHARACTER(LEN=13)      :: MERGE2

      ! For binary punch file, v. 2.0
      REAL*4                 :: ARRAY(IIPAR, JJPAR, LLPAR)
      REAL*4                 :: LONRES, LATRES

      INTEGER                :: IFIRST, JFIRST, LFIRST
      INTEGER, PARAMETER     :: HALFPOLAR = 1
      INTEGER, PARAMETER     :: CENTER180 = 1

      CHARACTER (LEN=20)     :: MODELNAME 
      CHARACTER (LEN=40)     :: UNIT
      CHARACTER (LEN=40)     :: RESERVED = ''
      CHARACTER (LEN=40)     :: CATEGORY 
      CHARACTER (LEN=80)     :: LABEL

      ! External functions 
      REAL*8, EXTERNAL       :: BOXVL 

      !=================================================================
      ! CH4_BUDGET begins here!
      !
      ! Initialize quantities 
      !=================================================================
      IFIRST    = GET_XOFFSET() + 1
      JFIRST    = GET_YOFFSET() + 1
      LFIRST    = 1
      LONRES    = DISIZE
      LATRES    = DJSIZE

      ! Current month
      LMN       = GET_MONTH()

      ! Make up a category name for GAMAP (use 8 characters)
      CATEGORY  = 'CH4BUDGT'

      ! Get the proper model name for the binary punch file
      MODELNAME = GET_MODELNAME()

      ! Descriptor string
      LABEL    = 'GEOS-CHEM -- CH4 Budget output (jsw, bmy, 1/16/01)'

      ! Unit of quantity being saved
      UNIT      = 'Tg'  !(NOTE: check w/ bnd to get the right units!!!)

      ! Scale factor for dynamic time steps
      SCALEDYN  = FLOAT( GET_CT_DYN() ) + 1D-20

      !=================================================================
      ! Convert initial burden of CO in TCO(1) from volume mixing
      ! ratio to total molecules/box for first month of simulation.
      !=================================================================
      IF( NBUDGET == 1 ) THEN
         DO L=1,LLPAR
         DO J=1,JJPAR
         DO I=1,IIPAR
            TCH4(I,J,L,1) = TCH4(I,J,L,1) * 
     &                      ( BAIRDENS(I,J,L) * BOXVL(I,J,L) )
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! Store the final burden of CH4 in TCH4(2) 
      ! Convert kg CH4/box to molecules/box.
      !=================================================================
      TCH4(:,:,:,2) = 0d0
      TCH4(:,:,:,2) = STT(:,:,:,1) * XNUMOL_CH4

      !=================================================================
      ! Write GLOBAL AVERAGES for all layers to ASCII file
      !=================================================================
      WRITE( MERGE, 2 ) GET_MONTH(), GET_YEAR()
 2    FORMAT( 'CObudget.', I2.2, '.',I4 )

      OPEN( 189, FILE=MERGE, STATUS='UNKNOWN' )
      REWIND( 189 )
      
      TGS     = 1.D-9
      STTCONV = XNUMOL_CH4/TGS
      SOURCES = 0.D0
      SINKS   = 0.D0
      NERROR  = 0
      
      WRITE(189,18)
      WRITE(189,1801)
 1801 FORMAT('*************************')
      WRITE(189,1800)
 1800 FORMAT('LAYERS 1 - 20')
      WRITE(189,1801)
      WRITE(189,18)

      WRITE(189,18)
      WRITE(189,38)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)
 1990 FORMAT('Tropospheric Burden')
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1)
      WRITE(189,20)NTP,NTP/STTCONV
      WRITE(189,21)NTP2,NTP2/STTCONV

      WRITE(189,18)
      WRITE(189,1991)
 1991 FORMAT('Stratospheric Burden')

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,21) NTP,NTP/STTCONV
      
      WRITE(189,18)
      WRITE(189,31)

c Sinks   jsw has checked correctness of code for sinks.
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,3,3,1)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,3,3,0)
      SINKS=NTP+NTQ

      WRITE(189,22) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      WRITE(189,220) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(189,29) 
      WRITE(189,34) SINKS,SINKS/STTCONV  !Just OH sink 
      WRITE(189,18)
      WRITE(189,30)

C Sources
!jsw      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,5,9,1)
!jsw      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,4,4,0)
!jsw      SOURCES=NTQ+NTP
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,4,4,1)
!jsw      WRITE(189,23) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      SOURCES=NTP
!jsw      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,4,4,0)
!jsw      WRITE(189,230) NTQ,NTQ/SOURCES*100.D0,NTQ/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,5,5,1)
      WRITE(189,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,9,9,1)
      WRITE(189,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,6,6,1)
      WRITE(189,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,7,7,1)
      WRITE(189,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,8,8,1)
      WRITE(189,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

Cjsw Following lines added by jsw.
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,10,10,1)
      WRITE(189,35) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      SINKS=SINKS-NTP  !Minus sign because soil absorption is negative.

      WRITE(189,29) 
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,18)
      WRITE(189,288) NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS,
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
 288  FORMAT('Initial-Final+Sources-Sinks=',E10.3,2x,F10.3)
      WRITE(189,18)
      WRITE(189,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV
 289  FORMAT('Net Gain          : ',E10.3,10x,F10.3)
      
 18   FORMAT()
 19   FORMAT('                    #Molecules               TG')
 20   FORMAT('  Start of Month  :',E10.3,10x,F10.3)
 21   FORMAT('  End of Month    :',E10.3,10x,F10.3)
 31   FORMAT('SINKS                            %Sink')
 22   FORMAT('  CH4 decay-trop   :',E10.3,2x,F6.1,2x,F10.3)
 220  FORMAT('  CH4 decay-strat  :',E10.3,2x,F6.1,2x,F10.3)
 34   FORMAT('Total Sinks       :',E10.3,10x,F10.3)
 30   FORMAT('SOURCES                          %Source')
 23   FORMAT('  Bogs            :',E10.3,2x,F6.1,2x,F10.3)
 !230  FORMAT('  CH4 Ox.-strat   :',E10.3,2x,F6.1,2x,F10.3)
 24   FORMAT('  Animals         :',E10.3,2x,F6.1,2x,F10.3)
 39   FORMAT('  Rice            :',E10.3,2x,F6.1,2x,F10.3)
 25   FORMAT('  Gas leakage     :',E10.3,2x,F6.1,2x,F10.3)
 26   FORMAT('  Coal            :',E10.3,2x,F6.1,2x,F10.3)
 27   FORMAT('  Bogs            :',E10.3,2x,F6.1,2x,F10.3)
 35   FORMAT('  Soil absorption :',E10.3,2x,F6.1,2x,F10.3)
 270  FORMAT('  N-S Ex.-trop    :',E10.3,2x,F6.1,2x,F10.3)
 2700 FORMAT('  N-S Ex.-strat   :',E10.3,2x,F6.1,2x,F10.3)
 29   FORMAT('                     ---------')
 28   FORMAT('Total Sources     :',E10.3,10x,F10.3)
      
      !=================================================================
      ! Write SOUTHERN HEMISPHERE averages to ASCII file
      ! jsw:  I have not modified the remaining code for CH4.
      !================================================================= 

      SOURCES = 0.D0
      SINKS   = 0.D0

      WRITE(189,18)
      WRITE(189,18)
      WRITE(189,36)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1)
      WRITE(189,21) NTP,NTP/STTCONV

      WRITE(189,18)
      WRITE(189,1991)
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0)
      WRITE(189,21) NTP,NTP/STTCONV
      WRITE(189,18)
      WRITE(189,31)

      ! Sinks
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF( NTP > 0d0) SINKS = SINKS + NTP

      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF( NTP2 > 0d0 ) SINKS = SINKS + NTP2      

      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,1)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,0)
      SINKS=SINKS+NTQ+NTQ2
      WRITE(189,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(189,220) NTQ2,NTQ2/SINKS*100.D0,NTQ2/STTCONV
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.GE.0.D0) THEN
         WRITE(189,270) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.GE.0.D0) THEN
         WRITE(189,2700) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      ENDIF
      
      WRITE(189,29)
      WRITE(189,34) SINKS,SINKS/STTCONV
      WRITE(189,18)
      WRITE(189,30)
      
      ! Sources
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,9,1)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0)
      SOURCES=NTQ+NTQ2
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         SOURCES=SOURCES-NTP
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         SOURCES=SOURCES-NTP
      ENDIF
      
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0)
      WRITE(189,23) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
!jsw      WRITE(189,230) NTP2,NTP2/SOURCES*100.D0,NTP2/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,5,5,1)
      WRITE(189,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,9,9,1)
      WRITE(189,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,6,6,1)
      WRITE(189,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,7,7,1)
      WRITE(189,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,8,8,1)
      WRITE(189,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,270) -NTP,-NTP/SOURCES*100.D0,-NTP/STTCONV
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,2700) -NTP,-NTP/SOURCES*100.D0,-NTP/STTCONV
      ENDIF
      
      WRITE(189,29)
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)
      
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0)
      WRITE(189,18)
      WRITE(189,288) (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS),
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(189,18)
      WRITE(189,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV

      !=================================================================
      ! Write NORTHERN HEMISPHERE averages to ASCII file 
      ! jsw:  I have not modified the remaining code for CH4.
      !================================================================= 

      SOURCES = 0.D0
      SINKS   = 0.D0

      WRITE(189,18)
      WRITE(189,18)
      WRITE(189,37)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1)
      WRITE(189,20) NTP,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1)
      WRITE(189,21) NTP,NTP/STTCONV
      
      WRITE(189,18)
      WRITE(189,1991)
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0)
      WRITE(189,20) NTP,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,21) NTP,NTP/STTCONV
      
      WRITE(189,18)
      WRITE(189,31)
c Sinks
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,1)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,0)
      SINKS=NTQ+NTQ2

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         SINKS=SINKS-NTP
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         SINKS=SINKS-NTP
      ENDIF

      WRITE(189,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(189,220) NTQ2,NTQ2/SINKS*100.D0,NTQ2/STTCONV
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,270) -NTP,-NTP/SINKS*100.D0,-NTP/STTCONV
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,2700) -NTP,-NTP/SINKS*100.D0,-NTP/STTCONV
      ENDIF
      WRITE(189,29)
      WRITE(189,34)SINKS,SINKS/STTCONV
      WRITE(189,18)
      WRITE(189,30)
C Sources
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,9,1)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0)
      SOURCES=NTQ+NTQ2

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.GE.0.D0) THEN
         SOURCES=SOURCES+NTP
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.GE.0.D0) THEN
         SOURCES=SOURCES+NTP
      ENDIF
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,1)
      WRITE(189,23) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0)
!jsw      WRITE(189,230) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,5,5,1)
      WRITE(189,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,9,9,1)
      WRITE(189,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,6,6,1)
      WRITE(189,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,7,7,1)
      WRITE(189,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,8,8,1)
      WRITE(189,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.GE.0.D0) THEN
         WRITE(189,270) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.GE.0.D0) THEN
         WRITE(189,2700) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      ENDIF
      WRITE(189,29)
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1)
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,18)
      WRITE(189,288) (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS),
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(189,18)
      WRITE(189,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV
 36   FORMAT('*****  Southern Hemisphere  *****')
 37   FORMAT('*****  Northern Hemisphere  *****')
 38   FORMAT('*****  Global  *****')
      
      CLOSE(189)

      !=================================================================
      ! Also save to binary punch file
      !=================================================================
      CALL BPCH2_HDR( 190, LABEL )

      DO K = 1, N_CH4
        
         ! Cast REAL*8 into REAL*4, convert from molec to Tg
         ARRAY(:,:,:) = TCH4(:,:,:,K) / STTCONV
        
         ! Save the data block 
         CALL BPCH2( 190,       MODELNAME,   LONRES,      LATRES,
     &               HALFPOLAR, CENTER180,   CATEGORY,    K,    
     &               UNIT,      GET_DIAGB(), GET_DIAGb(), RESERVED, 
     &               IIPAR,     JJPAR,       LLPAR,       IFIRST,  
     &               JFIRST,    LFIRST,      ARRAY )
      ENDDO

      CLOSE(190)

      !=================================================================
      ! Final burden at last of month equals initial burden
      ! of next month.  Also set TCH4 = 0 for next month.
      !=================================================================
      TCH4(:,:,:,1      ) = TCH4(:,:,:,2)
      TCH4(:,:,:,2:N_CH4) = 0d0

      ! Return to calling program
      END SUBROUTINE CH4_BUDGET

!------------------------------------------------------------------------------

      REAL*8 FUNCTION SUM_CH4( I1, I2, J1, J2, L1, L2, K1, K2, UPDOWN )
!
!******************************************************************************
!  Function SUM_CH4 sums a section of the TCH4 array bounded by the input
!  variables I1, I2, J1, J2, L1, L2, K1, K2.  SUM_CH4 is called by
!  module subroutine CH4_BUDGET. (jsw, bnd, bmy, 1/16/01)
!
!  Store the sources/sinks of CH4 in TCH4 in total molecules
!           ( 1) = Initial burden
!           ( 2) = Final burden
!  SINKS
!           ( 3) = CH4 sink by OH
!  SOURCES
!           ( 4) = Total Source
!           ( 5) = Animals
!           ( 6) = Gas Leakage
!           ( 7) = Coal
!           ( 8) = Bogs 
!           ( 9) = Rice
!           (10) = Soil absorption
!           (11) = Interhemispheric Exchange (+ = northward)
!           (12) = Bogs
!
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric (bmy, 4/17/00)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I1, I2 (INTEGER) : Min and max longitude indices of TCH4 to sum
!  (3-4) J1, J2 (INTEGER) : Min and max latitude  indices of TCH4 to sum 
!  (5-6) L1, L2 (INTEGER) : Min and max altitude  indices of TCH4 to sum 
!  (7-8) K1, K2 (INTEGER) : Min and max tracer    indices of TCH4 to sum
!  (9  ) UPDOWN (INTEGER) : Sum in troposphere (=1) or in stratosphere (=0)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_BUDGET is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/12/01)
!******************************************************************************
!    
#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! LPAUSE

      ! Arguments
      INTEGER, INTENT(IN) :: I1, I2, J1, J2, L1, L2
      INTEGER, INTENT(IN) :: K1, K2, UPDOWN
      
      ! Local variables
      INTEGER             :: I, J, K, L, LPAUSE_MIN, LPAUSE_MAX

      !=================================================================
      ! SUM_CH4 begins here!
      !=================================================================

      ! Compute the minimum value of LPAUSE once for use in
      ! the DO-loops below (bmy, 4/18/00)
      LPAUSE_MIN = MINVAL( LPAUSE )
      LPAUSE_MAX = MAXVAL( LPAUSE )

      !### Debug
      !print*,'LPAUSE MIN/MAX=',LPAUSE_MIN,LPAUSE_MAX  
      !print*,'L1,L2=',L1,L2
      
      ! Initialize SUM_CH4
      SUM_CH4 = 0d0

      ! Test on UPDOWN
      IF ( UPDOWN == 1 ) THEN

         !=============================================================
         ! UPDOWN = 1: Sum up from the surface to the tropopause
         !=============================================================
         DO K = K1, K2
         DO L = L1, LPAUSE_MAX
         DO J = J1, J2
         DO I = I1, I2
            IF ( L < LPAUSE(I,J) ) THEN 
               SUM_CH4 = SUM_CH4 + TCH4(I,J,L,K)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      ELSE

         !=============================================================
         ! UPDOWN = 0: Sum up from the tropopause to the atm top
         !=============================================================
         DO K = K1,         K2
         DO L = LPAUSE_MIN, L2
         DO J = J1,         J2
         DO I = I1,         I2
            IF ( L >= LPAUSE(I,J) ) THEN 
               SUM_CH4 = SUM_CH4 + TCH4(I,J,L,K)
            ENDIF            
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      
      ! Return to calling program
      END FUNCTION SUM_CH4

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GLOBAL_CH4
!
!******************************************************************************
!  Subroutine INIT_GLOBAL_CH4 allocates and zeroes module arrays. 
!  (bmy, 1/16/01, 10/15/02)
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

      ALLOCATE( AVGOH( NSEAS, NCMSLATS, NCMSALTS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVGOH' )
      AVGOH = 0d0

      ALLOCATE( BAIRDENS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BAIRDENS' )
      BAIRDENS = 0d0

      ALLOCATE( BOH( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BOH' )
      BOH = 0d0

      ALLOCATE( COPROD( JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COPROD' )
      COPROD = 0d0

      ALLOCATE( PAVG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PAVG' )
      PAVG = 0d0

      ALLOCATE( TAVG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAVG' )
      TAVG = 0d0

      ALLOCATE( TCH4( IIPAR, JJPAR, LLPAR, N_CH4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCH4' )
      TCH4 = 0d0      

      ! Return to calling program
      END SUBROUTINE INIT_GLOBAL_CH4

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GLOBAL_CH4
!
!******************************************************************************
!  Subroutine CLEANUP_GLOBAL_CH4 deallocates module arrays. (bmy, 1/16/01)
!******************************************************************************
! 
      IF ( ALLOCATED( BAIRDENS ) ) DEALLOCATE( BAIRDENS )
      IF ( ALLOCATED( BOH      ) ) DEALLOCATE( BOH      )
      IF ( ALLOCATED( COPROD   ) ) DEALLOCATE( COPROD   )
      IF ( ALLOCATED( TCH4     ) ) DEALLOCATE( TCH4     )

      END SUBROUTINE CLEANUP_GLOBAL_CH4

!------------------------------------------------------------------------------

      END MODULE GLOBAL_CH4_MOD
