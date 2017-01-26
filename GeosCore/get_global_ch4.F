!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: get_global_ch4.F
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_CH4 computes the latitudinal gradient 
!  in CH4 corresponding to year.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_CH4( THISYEAR,  VARIABLE_CH4, 
     &                           A3090S,    A0030S, 
     &                           A0030N,    A3090N,
     &                           am_I_Root, Input_Opt )
!
! !USES:
!
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCENARIO
      USE Input_Opt_Mod,        ONLY : OptInput
      USE PRECISION_MOD 

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: THISYEAR      ! Current year
      LOGICAL,        INTENT(IN)  :: VARIABLE_CH4  ! =T: Use time-varying CH4
                                                   ! =F: Use constant CH4
      LOGICAL,        INTENT(IN)  :: am_I_Root     ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt     ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(OUT) :: A3090S        ! CH4 [ppbv], 90S - 30S lat
      REAL(fp),       INTENT(OUT) :: A0030S        ! CH4 [ppbv], 30S - 00  lat
      REAL(fp),       INTENT(OUT) :: A0030N        ! CH4 [ppbv], 00  - 30N lat
      REAL(fp),       INTENT(OUT) :: A3090N        ! CH4 [ppbv], 30N - 90N lat
!
! !REMARKS:
!  (1) 2007 was the prior default year for CH4 emissions.  This has now been 
!       changed to 2013.
!  (2) The FUTURE_SCENARIO option is probably now obsolete, as most scale
!       factors are now handled by HEMCO.  But we shall leave this intact 
!       until further notice.
!
! !REVISION HISTORY:
!  03 Jan 2001 - J. Wang     - Initial version
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
!  13 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  06 Jan 2015 - R. Yantosca - Initial version
!  12 May 2015 - K. Travis   - Updated CH4 values with Version: 2014-06-24 
!                              for 2008-2013.
!  12 May 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=2) :: FUTURE_SCENARIO

      !=================================================================
      ! GET_GLOBAL_CH4 begins here!
      !
      ! New methane data from 1983-2007 (mje, bmy, 7/7/03, jaf 4/15/09)
      ! version is 2008-07-02
      ! Added data from 2008-2013  (krt, 7/7/14) Version: 2014-06-24 
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
               A3090S = 700.0e+0_fp
               A0030S = 700.0e+0_fp
               A0030N = 700.0e+0_fp
               A3090N = 700.0e+0_fp
               
            ! Modern-day years ...
            CASE( 1983 )
               A3090S = 1583.48e+0_fp
               A0030S = 1598.24e+0_fp
               A0030N = 1644.37e+0_fp
               A3090N = 1706.48e+0_fp

            CASE( 1984 )
               A3090S = 1597.77e+0_fp
               A0030S = 1606.66e+0_fp
               A0030N = 1655.62e+0_fp
               A3090N = 1723.63e+0_fp

            CASE( 1985 )
               A3090S = 1608.08e+0_fp
               A0030S = 1620.43e+0_fp
               A0030N = 1668.11e+0_fp
               A3090N = 1736.78e+0_fp
               
            CASE( 1986 )
               A3090S = 1619.91e+0_fp  
               A0030S = 1632.24e+0_fp  
               A0030N = 1682.88e+0_fp  
               A3090N = 1752.71e+0_fp

            CASE( 1987 )
               A3090S = 1630.54e+0_fp  
               A0030S = 1640.54e+0_fp  
               A0030N = 1702.05e+0_fp  
               A3090N = 1763.03e+0_fp

            CASE( 1988 )
               A3090S = 1642.08e+0_fp  
               A0030S = 1651.60e+0_fp  
               A0030N = 1713.07e+0_fp  
               A3090N = 1775.66e+0_fp

            CASE( 1989 )
               A3090S = 1654.03e+0_fp  
               A0030S = 1666.12e+0_fp  
               A0030N = 1720.53e+0_fp  
               A3090N = 1781.83e+0_fp

            CASE( 1990 )
               A3090S = 1663.21e+0_fp  
               A0030S = 1672.45e+0_fp  
               A0030N = 1733.84e+0_fp  
               A3090N = 1791.92e+0_fp

            CASE( 1991 )
               A3090S = 1673.52e+0_fp  
               A0030S = 1683.87e+0_fp  
               A0030N = 1750.68e+0_fp  
               A3090N = 1800.90e+0_fp

            CASE( 1992 )
               A3090S = 1687.97e+0_fp  
               A0030S = 1692.97e+0_fp  
               A0030N = 1755.94e+0_fp  
               A3090N = 1807.16e+0_fp

            CASE( 1993 )
               A3090S = 1687.83e+0_fp  
               A0030S = 1696.48e+0_fp  
               A0030N = 1758.86e+0_fp  
               A3090N = 1810.99e+0_fp

            CASE( 1994 )
               A3090S = 1692.00e+0_fp  
               A0030S = 1701.41e+0_fp  
               A0030N = 1766.98e+0_fp  
               A3090N = 1817.12e+0_fp

            CASE( 1995 )
               A3090S = 1701.04e+0_fp  
               A0030S = 1709.07e+0_fp  
               A0030N = 1778.25e+0_fp  
               A3090N = 1822.04e+0_fp

            CASE( 1996 )
               A3090S = 1701.87e+0_fp  
               A0030S = 1711.01e+0_fp  
               A0030N = 1778.08e+0_fp  
               A3090N = 1825.23e+0_fp
            
            CASE( 1997 )
               A3090S = 1708.01e+0_fp  
               A0030S = 1713.91e+0_fp  
               A0030N = 1781.43e+0_fp  
               A3090N = 1825.15e+0_fp

            CASE( 1998 )
               A3090S = 1716.55e+0_fp  
               A0030S = 1724.57e+0_fp  
               A0030N = 1783.86e+0_fp  
               A3090N = 1839.72e+0_fp

            CASE( 1999 )
               A3090S = 1725.70e+0_fp  
               A0030S = 1734.06e+0_fp  
               A0030N = 1791.50e+0_fp  
               A3090N = 1842.59e+0_fp

            CASE( 2000 )
               A3090S = 1728.13e+0_fp  
               A0030S = 1737.70e+0_fp  
               A0030N = 1792.42e+0_fp  
               A3090N = 1840.83e+0_fp

            CASE( 2001 )
               A3090S = 1726.92e+0_fp  
               A0030S = 1730.72e+0_fp  
               A0030N = 1789.11e+0_fp  
               A3090N = 1841.85e+0_fp

            CASE( 2002 )
               A3090S = 1729.75e+0_fp
               A0030S = 1735.28e+0_fp
               A0030N = 1790.08e+0_fp
               A3090N = 1842.36e+0_fp

            CASE( 2003 )
               A3090S = 1729.64e+0_fp
               A0030S = 1735.49e+0_fp
               A0030N = 1795.89e+0_fp
               A3090N = 1853.97e+0_fp

            CASE( 2004 )
               A3090S = 1728.72e+0_fp
               A0030S = 1738.54e+0_fp
               A0030N = 1797.30e+0_fp
               A3090N = 1849.58e+0_fp

            CASE( 2005 )
               A3090S = 1727.10e+0_fp
               A0030S = 1734.65e+0_fp
               A0030N = 1795.73e+0_fp
               A3090N = 1849.79e+0_fp

            CASE( 2006 )
               A3090S = 1726.53e+0_fp
               A0030S = 1735.17e+0_fp
               A0030N = 1796.30e+0_fp
               A3090N = 1848.20e+0_fp

            CASE( 2007 )
               A3090S = 1732.52e+0_fp
               A0030S = 1741.68e+0_fp
               A0030N = 1801.38e+0_fp
               A3090N = 1855.55e+0_fp

            CASE( 2008 )
               A3090S = 1740.54e+0_fp
               A0030S = 1748.16e+0_fp
               A0030N = 1807.91e+0_fp
               A3090N = 1868.26e+0_fp 

            CASE( 2009 )
               A3090S = 1743.97e+0_fp
               A0030S = 1753.06e+0_fp 
               A0030N = 1816.44e+0_fp
               A3090N = 1872.97e+0_fp

            CASE( 2010 )
               A3090S = 1750.15e+0_fp
               A0030S = 1759.71e+0_fp 
               A0030N = 1820.60e+0_fp
               A3090N = 1878.29e+0_fp

            CASE( 2011 )
               A3090S = 1757.30e+0_fp
               A0030S = 1763.98e+0_fp
               A0030N = 1822.17e+0_fp
               A3090N = 1880.56e+0_fp

            CASE( 2012 )
               A3090S = 1761.27e+0_fp
               A0030S = 1767.83e+0_fp 
               A0030N = 1824.79e+0_fp
               A3090N = 1885.32e+0_fp

            CASE( 2013 )
               A3090S = 1765.66e+0_fp
               A0030S = 1773.37e+0_fp 
               A0030N = 1830.84e+0_fp
               A3090N = 1893.50e+0_fp 
               
            ! Future year 2030
            CASE( 2025:2035 )
            
               ! Pick the IPCC scenario.  If LFUTURE=F and FUTURE_SCENARIO
               ! are undefined, then we are running 2030 meteorology with 
               ! present-day emissions.  In this case, default to 2013 CH4 
               ! concentrations. (havala, 1/25/08)
               SELECT CASE( FUTURE_SCENARIO )
                  CASE( 'A1' )
                     A3090S = 2202.0e+0_fp 
                     A0030S = 2202.0e+0_fp
                     A0030N = 2202.0e+0_fp
                     A3090N = 2202.0e+0_fp 
                  CASE( 'B1' )
                     A3090S = 1927.0e+0_fp 
                     A0030S = 1927.0e+0_fp
                     A0030N = 1927.0e+0_fp
                     A3090N = 1927.0e+0_fp 
                  CASE( 'A2' )
                     ! Not defined yet
                  CASE( 'B2' )
                     ! Not defined yet
                  CASE DEFAULT
                     ! 2013 is the default
                     A3090S = 1765.66e+0_fp
                     A0030S = 1773.37e+0_fp 
                     A0030N = 1830.84e+0_fp
                     A3090N = 1893.50e+0_fp 
               END SELECT

            ! Future year 2050
            CASE( 2045:2055 )

               ! Pick the IPCC scenario.  If LFUTURE=F and FUTURE_SCENARIO
               ! is undefined, then we are running 2050 meteorology with 
               ! present-day emissions.  In this case, default to 2007 CH4 
               ! concentrations. (havala, 1/25/08)
               SELECT CASE ( FUTURE_SCENARIO )
                  CASE ( 'A1' )
                     A3090S = 2400.0e+0_fp 
                     A0030S = 2400.0e+0_fp
                     A0030N = 2400.0e+0_fp
                     A3090N = 2400.0e+0_fp 
                  CASE ( 'B1' )
                     A3090S = 1881.0e+0_fp 
                     A0030S = 1881.0e+0_fp
                     A0030N = 1881.0e+0_fp
                     A3090N = 1881.0e+0_fp 
                  CASE ( 'A2' )
                     A3090S = 2562.0e+0_fp 
                     A0030S = 2562.0e+0_fp
                     A0030N = 2562.0e+0_fp
                     A3090N = 2562.0e+0_fp
                  CASE ( 'B2' )
                     A3090S = 2363.0e+0_fp 
                     A0030S = 2363.0e+0_fp
                     A0030N = 2363.0e+0_fp
                     A3090N = 2363.0e+0_fp
                  CASE DEFAULT
                     ! 2013 is the default
                     A3090S = 1765.66e+0_fp
                     A0030S = 1773.37e+0_fp 
                     A0030N = 1830.84e+0_fp
                     A3090N = 1893.50e+0_fp 
               END SELECT

            ! Default is to use 2013 data for other years
            ! for which we do not yet have data (krt, 7/7/14)
            CASE DEFAULT
               A3090S = 1765.66e+0_fp
               A0030S = 1773.37e+0_fp 
               A0030N = 1830.84e+0_fp
               A3090N = 1893.50e+0_fp 
         END SELECT

      ELSE
         
         ! ...otherwise assume constant global CH4
         A3090S = 1700.0e+0_fp
         A0030S = 1700.0e+0_fp
         A0030N = 1700.0e+0_fp
         A3090N = 1700.0e+0_fp
         
      ENDIF

      !=================================================================
      ! Print the latitudinal CH4 gradient for this year to stdout
      !=================================================================
#if !defined ( ESMF_ ) 
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
         ! from 2013 as a proxy for years past 2013 (mje, bmy, 7/7/03,
         ! jaf, 4/15/09, krt, 7/7/14)
         IF ( THISYEAR > 2013 ) THEN
            WRITE( 6, 115 ) 
 115        FORMAT( /, 'Using CH4 gradient data from 2013 as a proxy',
     &           /, 'since 2013 is the last year with reported data!' )
         ENDIF

         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      ENDIF
!$OMP END CRITICAL
#endif

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_CH4
!EOC
