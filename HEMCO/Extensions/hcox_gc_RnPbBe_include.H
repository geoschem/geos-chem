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

!EOP
!------------------------------------------------------------------------------
!BOC

    ! Define latitudes [degrees North]
    LATSOU      = (/     0d0,     10d0,     20d0,    30d0,    40d0,   &
                        50d0,     60d0,     70d0,    80d0,    90d0  /)
 
    ! Define pressures [hPa]
    PRESOU      = (/     0d0,     50d0,     70d0,    90d0,   110d0,   &
                       130d0,    150d0,    170d0,   190d0,   210d0,   &
                       230d0,    250d0,    270d0,   290d0,   313d0,   &
                       338d0,    364d0,    392d0,   420d0,   451d0,   &
                       485d0,    518d0,    555d0,   592d0,   633d0,   &
                       680d0,    725d0,    772d0,   822d0,   875d0,   &
                       930d0,    985d0,   1030d0                    /)

    ! Define 7Be emissions [stars/g air/s]
    ! 1 "star" = 1 nuclear disintegration via cosmic rays
    !
    ! NOTE: These statements were defined from printout of the file
    ! and need to be multiplied by 1d-5 below.
    BESOU(:,1)  = (/   150d0,    156d0,    188d0,   285d0,   500d0,   &
                       910d0,   1700d0,   1900d0,  1900d0,  1900d0  /)
                       
    BESOU(:,2)  = (/   280d0,    310d0,    390d0,   590d0,   880d0,   &
                      1390d0,   1800d0,   1800d0,  1800d0,  1800d0  /)
                       
    BESOU(:,3)  = (/   310d0,    330d0,    400d0,   620d0,   880d0,   &
                      1280d0,   1450d0,   1450d0,  1450d0,  1450d0  /)
                       
    BESOU(:,4)  = (/   285d0,    310d0,    375d0,   570d0,   780d0,   &
                      1100d0,   1180d0,   1180d0,  1180d0,  1180d0  /)
                       
    BESOU(:,5)  = (/   255d0,    275d0,    330d0,   510d0,   680d0,   &
                       950d0,   1000d0,   1000d0,  1000d0,  1000d0  /)
                       
    BESOU(:,6)  = (/   230d0,    245d0,    292d0,   450d0,   600d0,   &
                       820d0,    875d0,    875d0,   875d0,   875d0  /)
                       
    BESOU(:,7)  = (/   205d0,    215d0,    260d0,   400d0,   530d0,   &
                       730d0,    750d0,    750d0,   750d0,   750d0  /)
                       
    BESOU(:,8)  = (/   182d0,    195d0,    235d0,   355d0,   480d0,   &
                       630d0,    650d0,    650d0,   650d0,   650d0  /)
                       
    BESOU(:,9)  = (/   160d0,    173d0,    208d0,   315d0,   410d0,   &
                       543d0,    550d0,    550d0,   550d0,   550d0  /)
                       
    BESOU(:,10) = (/   148d0,    152d0,    185d0,   280d0,   370d0,   &
                       480d0,    500d0,    500d0,   500d0,   500d0  /)
                       
    BESOU(:,11) = (/   130d0,    139d0,    167d0,   250d0,   320d0,   &
                       425d0,    430d0,    430d0,   430d0,   430d0  /)
                       
    BESOU(:,12) = (/   116d0,    123d0,    148d0,   215d0,   285d0,   &
                       365d0,    375d0,    375d0,   375d0,   375d0  /)
                       
    BESOU(:,13) = (/   104d0,    110d0,    130d0,   198d0,   250d0,   &
                       320d0,    330d0,    330d0,   330d0,   330d0  /)
                       
    BESOU(:,14) = (/    93d0,     99d0,    118d0,   170d0,   222d0,   &
                       280d0,    288d0,    288d0,   288d0,   288d0  /)
                       
    BESOU(:,15) = (/    80d0,     84d0,    100d0,   145d0,   190d0,   &
                       235d0,    250d0,    250d0,   250d0,   250d0  /)
                       
    BESOU(:,16) = (/    72d0,   74.5d0,     88d0,   129d0,   168d0,   &
                       210d0,    218d0,    218d0,   218d0,   218d0  /)
                       
    BESOU(:,17) = (/  59.5d0,   62.5d0,   73.5d0,   108d0,   138d0,   &
                       171d0,    178d0,    178d0,   178d0,   178d0  /)
                       
    BESOU(:,18) = (/    50d0,     53d0,     64d0,    90d0,   115d0,   &
                       148d0,    150d0,    150d0,   150d0,   150d0  /)
                       
    BESOU(:,19) = (/    45d0,   46.5d0,   52.5d0,    76d0,    98d0,   &
                       122d0,    128d0,    128d0,   128d0,   128d0  /)
                       
    BESOU(:,20) = (/  36.5d0,   37.5d0,     45d0,    61d0,    77d0,   &
                        98d0,    102d0,    102d0,   102d0,   102d0  /)
                       
    BESOU(:,21) = (/  30.8d0,     32d0,   37.5d0,  51.5d0,    65d0,   &
                        81d0,     85d0,     85d0,    85d0,    85d0  /)

    BESOU(:,22) = (/  25.5d0,   26.5d0,     32d0,  40.5d0,    54d0,   &
                      67.5d0,   69.5d0,   69.5d0,  69.5d0,  69.5d0  /)

    BESOU(:,23) = (/  20.5d0,   21.6d0,   25.5d0,    33d0,    42d0,   &
                      53.5d0,     55d0,     55d0,    55d0,    55d0  /)

    BESOU(:,24) = (/  16.8d0,   17.3d0,     20d0,    26d0,  33.5d0,   &
                        41d0,     43d0,     43d0,    43d0,    43d0  /)

    BESOU(:,25) = (/    13d0,   13.8d0,   15.3d0,  20.5d0,  26.8d0,   &
                      32.5d0,   33.5d0,   33.5d0,  33.5d0,  33.5d0  /)

    BESOU(:,26) = (/  10.1d0,   10.6d0,   12.6d0,  15.8d0,    20d0,   &
                      24.5d0,   25.8d0,   25.8d0,  25.8d0,  25.8d0  /)

    BESOU(:,27) = (/   7.7d0,   8.15d0,    9.4d0,  11.6d0,  14.8d0,   &
                      17.8d0,   18.5d0,   18.5d0,  18.5d0,  18.5d0  /)
 
    BESOU(:,28) = (/   5.7d0,   5.85d0,   6.85d0,  8.22d0,    11d0,   &
                      13.1d0,   13.2d0,   13.2d0,  13.2d0,  13.2d0  /)

    BESOU(:,29) = (/   3.9d0,    4.2d0,   4.85d0,     6d0,   7.6d0,   &
                         9d0,    9.2d0,    9.2d0,   9.2d0,   9.2d0  /)

    BESOU(:,30) = (/     3d0,   3.05d0,   3.35d0,   4.2d0,   5.3d0,   &
                       5.9d0,   6.25d0,   6.25d0,  6.25d0,  6.25d0  /)

    BESOU(:,31) = (/  2.05d0,    2.1d0,   2.32d0,   2.9d0,   3.4d0,   &
                       3.9d0,    4.1d0,    4.1d0,   4.1d0,   4.1d0  /)

    BESOU(:,32) = (/  1.45d0,   1.43d0,   1.65d0,  2.03d0,   2.4d0,   &
                      2.75d0,   2.65d0,   2.65d0,  2.65d0,  2.65d0  /)

    BESOU(:,33) = (/  1.04d0,   1.08d0,   1.21d0,   1.5d0,  1.68d0,   &
                       1.8d0,    1.8d0,    1.8d0,   1.8d0,   1.8d0  /)

    ! All the numbers of BESOU need to be multiplied by 1e-5 in order to put 
    ! them into the correct data range.  NOTE: This multiplication statement
    ! needs to be preserved here in order to  ensure identical output to the
    ! prior code! (bmy, 7/7/14)
    BESOU = BESOU * 1d-5
!EOC
