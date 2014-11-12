!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcox_gc_RnPbBe_include.H
!
! !DESCRIPTION: Include file containing 7Be emissions from Lal \& Peters 
!  (1967) on 33 pressure levels.  This data used to be read from a file,
!  but we have now hardwired it to facilitate I/O in the ESMF environment.
!\\
!\\
! !REMARKS:
!  (1) Reference: Lal, D., and B. Peters, Cosmic ray produced radioactivity 
!       on the Earth. Handbuch der Physik, 46/2, 551-612, edited by K. Sitte, 
!        Springer-Verlag, New York, 1967.
!
!  (2) In prior versions of GEOS-Chem, this routine was named READ_7BE, and
!      it read the ASCII file "7Be.Lal".   Because this data set is not placed
!      on a lat/lon grid, ESMF cannot regrid it.  To work around this, we now
!      hardwire this data in module arrays rather than read it from disk.
!                                                                             .
!  (3) Units of 7Be emissions are [stars/g air/s].  
!      Here, "stars" = # of nuclear disintegrations of cosmic rays
!                                                                             .
!  (4) Original data from Lal & Peters (1967), w/ these modifications:
!      (a) Replace data at (0hPa, 70S) following Koch 1996:
!          (i ) old value = 3000 
!          (ii) new value = 1900
!      (b) Copy data from 70S to 80S and 90S at all levels
!                                                                             .
! !REVISION HISTORY: 
!  07 Aug 2002 - H. Liu - Initial version
!  (1 ) This code was split off from routine EMISSRnPbBe below. (bmy, 8/7/02)
!  (2 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/19/04)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  02 Jul 2014 - R. Yantosca - Now hardwire the data instead of reading it
!                              from an ASCII file.  This facilitates ESMF I/O.
!  07 Jul 2014 - R. Yantosca - Now renamed to INIT_7Be_Emissions and added
!                              as a HEMCO extension
!  07 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!   8 Aug 2014 - R. Yantosca - Now split off into hcox_gc_RnPbBe_include.H
!  05 Nov 2014 - C. Keller   - Converted from double-precision to flexible
!                              (HEMCO) precision hp.
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Define latitudes [degrees North]
    LATSOU      = (/     0.0_hp,     10.0_hp,     20.0_hp,    30.0_hp,    40.0_hp,   &
                        50.0_hp,     60.0_hp,     70.0_hp,    80.0_hp,    90.0_hp  /)
 
    ! Define pressures [hPa]
    PRESOU      = (/     0.0_hp,     50.0_hp,     70.0_hp,    90.0_hp,   110.0_hp,   &
                       130.0_hp,    150.0_hp,    170.0_hp,   190.0_hp,   210.0_hp,   &
                       230.0_hp,    250.0_hp,    270.0_hp,   290.0_hp,   313.0_hp,   &
                       338.0_hp,    364.0_hp,    392.0_hp,   420.0_hp,   451.0_hp,   &
                       485.0_hp,    518.0_hp,    555.0_hp,   592.0_hp,   633.0_hp,   &
                       680.0_hp,    725.0_hp,    772.0_hp,   822.0_hp,   875.0_hp,   &
                       930.0_hp,    985.0_hp,   1030.0_hp                    /)

    ! Define 7Be emissions [stars/g air/s]
    ! 1 "star" = 1 nuclear disintegration via cosmic rays
    !
    ! NOTE: These statements were defined from printout of the file
    ! and need to be multiplied by 1d-5 below.
    BESOU(:,1)  = (/   150.0_hp,    156.0_hp,    188.0_hp,   285.0_hp,   500.0_hp,   &
                       910.0_hp,   1700.0_hp,   1900.0_hp,  1900.0_hp,  1900.0_hp  /)
                       
    BESOU(:,2)  = (/   280.0_hp,    310.0_hp,    390.0_hp,   590.0_hp,   880.0_hp,   &
                      1390.0_hp,   1800.0_hp,   1800.0_hp,  1800.0_hp,  1800.0_hp  /)
                       
    BESOU(:,3)  = (/   310.0_hp,    330.0_hp,    400.0_hp,   620.0_hp,   880.0_hp,   &
                      1280.0_hp,   1450.0_hp,   1450.0_hp,  1450.0_hp,  1450.0_hp  /)
                       
    BESOU(:,4)  = (/   285.0_hp,    310.0_hp,    375.0_hp,   570.0_hp,   780.0_hp,   &
                      1100.0_hp,   1180.0_hp,   1180.0_hp,  1180.0_hp,  1180.0_hp  /)
                       
    BESOU(:,5)  = (/   255.0_hp,    275.0_hp,    330.0_hp,   510.0_hp,   680.0_hp,   &
                       950.0_hp,   1000.0_hp,   1000.0_hp,  1000.0_hp,  1000.0_hp  /)
                       
    BESOU(:,6)  = (/   230.0_hp,    245.0_hp,    292.0_hp,   450.0_hp,   600.0_hp,   &
                       820.0_hp,    875.0_hp,    875.0_hp,   875.0_hp,   875.0_hp  /)
                       
    BESOU(:,7)  = (/   205.0_hp,    215.0_hp,    260.0_hp,   400.0_hp,   530.0_hp,   &
                       730.0_hp,    750.0_hp,    750.0_hp,   750.0_hp,   750.0_hp  /)
                       
    BESOU(:,8)  = (/   182.0_hp,    195.0_hp,    235.0_hp,   355.0_hp,   480.0_hp,   &
                       630.0_hp,    650.0_hp,    650.0_hp,   650.0_hp,   650.0_hp  /)
                       
    BESOU(:,9)  = (/   160.0_hp,    173.0_hp,    208.0_hp,   315.0_hp,   410.0_hp,   &
                       543.0_hp,    550.0_hp,    550.0_hp,   550.0_hp,   550.0_hp  /)
                       
    BESOU(:,10) = (/   148.0_hp,    152.0_hp,    185.0_hp,   280.0_hp,   370.0_hp,   &
                       480.0_hp,    500.0_hp,    500.0_hp,   500.0_hp,   500.0_hp  /)
                       
    BESOU(:,11) = (/   130.0_hp,    139.0_hp,    167.0_hp,   250.0_hp,   320.0_hp,   &
                       425.0_hp,    430.0_hp,    430.0_hp,   430.0_hp,   430.0_hp  /)
                       
    BESOU(:,12) = (/   116.0_hp,    123.0_hp,    148.0_hp,   215.0_hp,   285.0_hp,   &
                       365.0_hp,    375.0_hp,    375.0_hp,   375.0_hp,   375.0_hp  /)
                       
    BESOU(:,13) = (/   104.0_hp,    110.0_hp,    130.0_hp,   198.0_hp,   250.0_hp,   &
                       320.0_hp,    330.0_hp,    330.0_hp,   330.0_hp,   330.0_hp  /)
                       
    BESOU(:,14) = (/    93.0_hp,     99.0_hp,    118.0_hp,   170.0_hp,   222.0_hp,   &
                       280.0_hp,    288.0_hp,    288.0_hp,   288.0_hp,   288.0_hp  /)
                       
    BESOU(:,15) = (/    80.0_hp,     84.0_hp,    100.0_hp,   145.0_hp,   190.0_hp,   &
                       235.0_hp,    250.0_hp,    250.0_hp,   250.0_hp,   250.0_hp  /)
                       
    BESOU(:,16) = (/    72.0_hp,     74.0_hp,     88.0_hp,   129.0_hp,   168.0_hp,   &
                       210.0_hp,    218.0_hp,    218.0_hp,   218.0_hp,   218.0_hp  /)
                       
    BESOU(:,17) = (/    59.5_hp,     62.5_hp,     73.5_hp,   108.0_hp,   138.0_hp,   &
                       171.0_hp,    178.0_hp,    178.0_hp,   178.0_hp,   178.0_hp  /)
                       
    BESOU(:,18) = (/    50.0_hp,     53.0_hp,     64.0_hp,    90.0_hp,   115.0_hp,   &
                       148.0_hp,    150.0_hp,    150.0_hp,   150.0_hp,   150.0_hp  /)
                       
    BESOU(:,19) = (/    45.0_hp,   46.5_hp,   52.5_hp,    76.0_hp,    98.0_hp,   &
                       122.0_hp,    128.0_hp,    128.0_hp,   128.0_hp,   128.0_hp  /)
                       
    BESOU(:,20) = (/  36.5_hp,   37.5_hp,     45.0_hp,    61.0_hp,    77.0_hp,   &
                        98.0_hp,    102.0_hp,    102.0_hp,   102.0_hp,   102.0_hp  /)
                       
    BESOU(:,21) = (/  30.8_hp,     32.0_hp,   37.5_hp,  51.5_hp,    65.0_hp,   &
                        81.0_hp,     85.0_hp,     85.0_hp,    85.0_hp,    85.0_hp  /)

    BESOU(:,22) = (/  25.5_hp,   26.5_hp,     32.0_hp,  40.5_hp,    54_hp,   &
                      67.5_hp,   69.5_hp,   69.5_hp,  69.5_hp,  69.5_hp  /)

    BESOU(:,23) = (/  20.5_hp,   21.6_hp,   25.5_hp,    33.0_hp,    42.0_hp,   &
                      53.5_hp,     55.0_hp,     55.0_hp,    55.0_hp,    55.0_hp  /)

    BESOU(:,24) = (/  16.8_hp,   17.3_hp,     20.0_hp,    26.0_hp,  33.5_hp,   &
                        41.0_hp,     43.0_hp,     43.0_hp,    43.0_hp,    43.0_hp  /)

    BESOU(:,25) = (/    13.0_hp,   13.8_hp,   15.3_hp,  20.5_hp,  26.8_hp,   &
                      32.5_hp,   33.5_hp,   33.5_hp,  33.5_hp,  33.5_hp  /)

    BESOU(:,26) = (/  10.1_hp,   10.6_hp,   12.6_hp,  15.8_hp,    20.0_hp,   &
                      24.5_hp,   25.8_hp,   25.8_hp,  25.8_hp,  25.8_hp  /)

    BESOU(:,27) = (/   7.7_hp,   8.15_hp,    9.4_hp,  11.6_hp,  14.8_hp,   &
                      17.8_hp,   18.5_hp,   18.5_hp,  18.5_hp,  18.5_hp  /)
 
    BESOU(:,28) = (/   5.7_hp,   5.85_hp,   6.85_hp,  8.22_hp,    11.0_hp,   &
                      13.1_hp,   13.2_hp,   13.2_hp,  13.2_hp,  13.2_hp  /)

    BESOU(:,29) = (/   3.9_hp,    4.2_hp,   4.85_hp,     6.0_hp,   7.6_hp,   &
                         9.0_hp,    9.2_hp,    9.2_hp,   9.2_hp,   9.2_hp  /)

    BESOU(:,30) = (/     3.0_hp,   3.05_hp,   3.35_hp,   4.2_hp,   5.3_hp,   &
                       5.9_hp,   6.25_hp,   6.25_hp,  6.25_hp,  6.25_hp  /)

    BESOU(:,31) = (/  2.05_hp,    2.1_hp,   2.32_hp,   2.9_hp,   3.4_hp,   &
                       3.9_hp,    4.1_hp,    4.1_hp,   4.1_hp,   4.1_hp  /)

    BESOU(:,32) = (/  1.45_hp,   1.43_hp,   1.65_hp,  2.03_hp,   2.4_hp,   &
                      2.75_hp,   2.65_hp,   2.65_hp,  2.65_hp,  2.65_hp  /)

    BESOU(:,33) = (/  1.04_hp,   1.08_hp,   1.21_hp,   1.5_hp,  1.68_hp,   &
                       1.8_hp,    1.8_hp,    1.8_hp,   1.8_hp,   1.8_hp  /)

    ! All the numbers of BESOU need to be multiplied by 1e-5 in order to put 
    ! them into the correct data range.  NOTE: This multiplication statement
    ! needs to be preserved here in order to  ensure identical output to the
    ! prior code! (bmy, 7/7/14)
    BESOU = BESOU * 1.e-5_hp
!EOC
