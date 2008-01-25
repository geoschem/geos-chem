! $Id: get_global_ch4.f,v 1.3 2008/01/25 19:48:16 bmy Exp $
      SUBROUTINE GET_GLOBAL_CH4( THISYEAR, VARIABLE_CH4, 
     &                           A3090S, A0030S, A0030N, A3090N )
!
!******************************************************************************
!  Subroutine GET_GLOBAL_CH4 computes the latitudinal gradient in CH4
!  corresponding to year (jsw, bnd, bmy, 1/3/01, 1/25/08)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) THISYEAR     (INTEGER) : Current month number (1-12)
!  (2 ) VARIABLE_CH4 (LOGICAL) : Flag for selecting variable or constant CH4
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
!  (5 ) Split off from module "global_ch4_mod.f".  Updated for IPCC future
!        emissions scenarios. (swu, bmy, 5/30/06)     
!  (6 ) Add the preindustrial CH4 scenarios.  Also set 2001 as the default
!        in case we are running 2030 or 2050 met but present-day emissions.
!        (swu, havala, bmy, 1/25/08)
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCENARIO
      USE LOGICAL_MOD,          ONLY : LFUTURE

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)           :: THISYEAR
      LOGICAL, INTENT(IN)           :: VARIABLE_CH4
      REAL*8,  INTENT(OUT)          :: A3090S, A0030S, A0030N, A3090N

      ! Local variables
      CHARACTER(LEN=2)              :: FUTURE_SCENARIO

      !=================================================================
      ! GET_GLOBAL_CH4 begins here!
      !
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
      !
      ! Also add future emission scenarios for GCAP, as well as
      ! the preindustrial CH4 levels (swu, havala, bmy, 1/25/08)
      !=================================================================
      IF ( VARIABLE_CH4 ) THEN

         ! Get IPCC future scenario (e.g. A1, A2, B1, B2)
         IF ( LFUTURE ) THEN
            FUTURE_SCENARIO = GET_FUTURE_SCENARIO()
         ENDIF

         ! Select latitudinal CH4 gradient by year...
         SELECT CASE ( THISYEAR )

            ! Preindustrial years
            CASE ( :1750 )
               A3090S = 700.0d0
               A0030S = 700.0d0
               A0030N = 700.0d0
               A3090N = 700.0d0
               
            ! Modern-day years ...
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

            ! Future year 2030
            CASE( 2025:2035 )
            
               ! Pick the IPCC scenario.  If LFUTURE=F and FUTURE_SCENARIO
               ! are undefined, then we are running 2030 meteorology with 
               ! present-day emissions.  In this case, default to 2001 CH4 
               ! concentrations. (havala, 1/25/08)
               SELECT CASE( FUTURE_SCENARIO )
                  CASE( 'A1' )
                     A3090S = 2202.0d0 
                     A0030S = 2202.0d0
                     A0030N = 2202.0d0
                     A3090N = 2202.0d0 
                  CASE( 'B1' )
                     A3090S = 1927.0d0 
                     A0030S = 1927.0d0
                     A0030N = 1927.0d0
                     A3090N = 1927.0d0 
                  CASE( 'A2' )
                     ! Not defined yet
                  CASE( 'B2' )
                     ! Not defined yet
                  CASE DEFAULT
                     ! 2001 is the default
                     A3090S = 1705.68d0  
                     A0030S = 1709.52d0  
                     A0030N = 1767.51d0  
                     A3090N = 1822.53d0
               END SELECT

            ! Future year 2050
            CASE( 2045:2055 )

               ! Pick the IPCC scenario.  If LFUTURE=F and FUTURE_SCENARIO
               ! is undefined, then we are running 2050 meteorology with 
               ! present-day emissions.  In this case, default to 2001 CH4 
               ! concentrations. (havala, 1/25/08)
               SELECT CASE ( FUTURE_SCENARIO )
                  CASE ( 'A1' )
                     A3090S = 2400.0d0 
                     A0030S = 2400.0d0
                     A0030N = 2400.0d0
                     A3090N = 2400.0d0 
                  CASE ( 'B1' )
                     A3090S = 1881.0d0 
                     A0030S = 1881.0d0
                     A0030N = 1881.0d0
                     A3090N = 1881.0d0 
                  CASE ( 'A2' )
                     A3090S = 2562.0d0 
                     A0030S = 2562.0d0
                     A0030N = 2562.0d0
                     A3090N = 2562.0d0
                  CASE ( 'B2' )
                     A3090S = 2363.0d0 
                     A0030S = 2363.0d0
                     A0030N = 2363.0d0
                     A3090N = 2363.0d0
                  CASE DEFAULT
                     ! 2001 is the default
                     A3090S = 1705.68d0  
                     A0030S = 1709.52d0  
                     A0030N = 1767.51d0  
                     A3090N = 1822.53d0
               END SELECT

            ! Default is to use 2001 data for other years
            ! for which we do not yet have data (bmy, 5/30/06)
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
