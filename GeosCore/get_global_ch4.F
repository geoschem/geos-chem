! $Id: get_global_ch4.f,v 1.1 2009/09/16 14:06:26 bmy Exp $
      SUBROUTINE GET_GLOBAL_CH4( THISYEAR, VARIABLE_CH4, 
     &                           A3090S, A0030S, A0030N, A3090N,
     &                           am_I_Root, Input_Opt )
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
!  (7 ) Updated CH4 vales with version 2008-07-02 for 1983-2007. Also use
!        2007 for years past 2007 (jaf, 4/15/09) 
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LFUTURE instead of LFUTURE
!                              from logical_mod.F
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
! !USES:
!
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCENARIO
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: THISYEAR
      LOGICAL,        INTENT(IN)  :: VARIABLE_CH4
      LOGICAL,        INTENT(IN)  :: am_I_Root
      TYPE(OptInput), INTENT(IN)  :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(OUT) :: A3090S, A0030S, A0030N, A3090N
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=2)            :: FUTURE_SCENARIO

      !=================================================================
      ! GET_GLOBAL_CH4 begins here!
      !
      ! New methane data from 1983-2007 (mje, bmy, 7/7/03, jaf 4/15/09)
      ! version is 2008-07-02
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
         IF ( Input_Opt%LFUTURE ) THEN
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
               A3090S = 1583.48d0
               A0030S = 1598.24d0
               A0030N = 1644.37d0
               A3090N = 1706.48d0

            CASE( 1984 )
               A3090S = 1597.77d0
               A0030S = 1606.66d0
               A0030N = 1655.62d0
               A3090N = 1723.63d0

            CASE( 1985 )
               A3090S = 1608.08d0
               A0030S = 1620.43d0
               A0030N = 1668.11d0
               A3090N = 1736.78d0
               
            CASE( 1986 )
               A3090S = 1619.91d0  
               A0030S = 1632.24d0  
               A0030N = 1682.88d0  
               A3090N = 1752.71d0

            CASE( 1987 )
               A3090S = 1630.54d0  
               A0030S = 1640.54d0  
               A0030N = 1702.05d0  
               A3090N = 1763.03d0

            CASE( 1988 )
               A3090S = 1642.08d0  
               A0030S = 1651.60d0  
               A0030N = 1713.07d0  
               A3090N = 1775.66d0

            CASE( 1989 )
               A3090S = 1654.03d0  
               A0030S = 1666.12d0  
               A0030N = 1720.53d0  
               A3090N = 1781.83d0

            CASE( 1990 )
               A3090S = 1663.21d0  
               A0030S = 1672.45d0  
               A0030N = 1733.84d0  
               A3090N = 1791.92d0

            CASE( 1991 )
               A3090S = 1673.52d0  
               A0030S = 1683.87d0  
               A0030N = 1750.68d0  
               A3090N = 1800.90d0

            CASE( 1992 )
               A3090S = 1687.97d0  
               A0030S = 1692.97d0  
               A0030N = 1755.94d0  
               A3090N = 1807.16d0

            CASE( 1993 )
               A3090S = 1687.83d0  
               A0030S = 1696.48d0  
               A0030N = 1758.86d0  
               A3090N = 1810.99d0

            CASE( 1994 )
               A3090S = 1692.00d0  
               A0030S = 1701.41d0  
               A0030N = 1766.98d0  
               A3090N = 1817.12d0

            CASE( 1995 )
               A3090S = 1701.04d0  
               A0030S = 1709.07d0  
               A0030N = 1778.25d0  
               A3090N = 1822.04d0

            CASE( 1996 )
               A3090S = 1701.87d0  
               A0030S = 1711.01d0  
               A0030N = 1778.08d0  
               A3090N = 1825.23d0
            
            CASE( 1997 )
               A3090S = 1708.01d0  
               A0030S = 1713.91d0  
               A0030N = 1781.43d0  
               A3090N = 1825.15d0

            CASE( 1998 )
               A3090S = 1716.55d0  
               A0030S = 1724.57d0  
               A0030N = 1783.86d0  
               A3090N = 1839.72d0

            CASE( 1999 )
               A3090S = 1725.70d0  
               A0030S = 1734.06d0  
               A0030N = 1791.50d0  
               A3090N = 1842.59d0

            CASE( 2000 )
               A3090S = 1728.13d0  
               A0030S = 1737.70d0  
               A0030N = 1792.42d0  
               A3090N = 1840.83d0

            CASE( 2001 )
               A3090S = 1726.92d0  
               A0030S = 1730.72d0  
               A0030N = 1789.11d0  
               A3090N = 1841.85d0

            CASE( 2002 )
               A3090S = 1729.75d0
               A0030S = 1735.28d0
               A0030N = 1790.08d0
               A3090N = 1842.36d0

            CASE( 2003 )
               A3090S = 1729.64d0
               A0030S = 1735.49d0
               A0030N = 1795.89d0
               A3090N = 1853.97d0

            CASE( 2004 )
               A3090S = 1728.72d0
               A0030S = 1738.54d0
               A0030N = 1797.30d0
               A3090N = 1849.58d0

            CASE( 2005 )
               A3090S = 1727.10d0
               A0030S = 1734.65d0
               A0030N = 1795.73d0
               A3090N = 1849.79d0

            CASE( 2006 )
               A3090S = 1726.53d0
               A0030S = 1735.17d0
               A0030N = 1796.30d0
               A3090N = 1848.20d0

            CASE( 2007 )
               A3090S = 1732.52d0
               A0030S = 1741.68d0
               A0030N = 1801.38d0
               A3090N = 1855.55d0

            ! Future year 2030
            CASE( 2025:2035 )
            
               ! Pick the IPCC scenario.  If LFUTURE=F and FUTURE_SCENARIO
               ! are undefined, then we are running 2030 meteorology with 
               ! present-day emissions.  In this case, default to 2007 CH4 
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
                     ! 2007 is the default
                     A3090S = 1732.52d0
                     A0030S = 1741.68d0
                     A0030N = 1801.38d0
                     A3090N = 1855.55d0
               END SELECT

            ! Future year 2050
            CASE( 2045:2055 )

               ! Pick the IPCC scenario.  If LFUTURE=F and FUTURE_SCENARIO
               ! is undefined, then we are running 2050 meteorology with 
               ! present-day emissions.  In this case, default to 2007 CH4 
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
                     ! 2007 is the default
                     A3090S = 1732.52d0
                     A0030S = 1741.68d0
                     A0030N = 1801.38d0
                     A3090N = 1855.55d0
               END SELECT

            ! Default is to use 2007 data for other years
            ! for which we do not yet have data (jaf, 4/15/09)
            CASE DEFAULT
               A3090S = 1732.52d0
               A0030S = 1741.68d0
               A0030N = 1801.38d0
               A3090N = 1855.55d0

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
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 105   ) THISYEAR
 105     FORMAT( 'GET_GLOBAL_CH4: YEAR = ', i4 )

         WRITE( 6, 110 ) A3090N, A0030N, A0030S, A3090S 
 110     FORMAT( 'CH4 (90N - 30N) : ', f7.1, ' [ppbv]', /,
     &           'CH4 (30N - 00 ) : ', f7.1, ' [ppbv]', /,
     &           'CH4 (00  - 30S) : ', f7.1, ' [ppbv]', /,
     &           'CH4 (30S - 90S) : ', f7.1, ' [ppbv]' )

         ! Indicate to the log file if we are using CH4 gradient data
         ! from 2007 as a proxy for years past 2007 (mje, bmy, 7/7/03,
         ! jaf, 4/15/09)
         IF ( THISYEAR > 2007 ) THEN
            WRITE( 6, 115 ) 
 115        FORMAT( /, 'Using CH4 gradient data from 2007 as a proxy',
     &           /, 'since 2007 is the last year with reported data!' )
         ENDIF

         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      ENDIF
!$OMP END CRITICAL

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_CH4
