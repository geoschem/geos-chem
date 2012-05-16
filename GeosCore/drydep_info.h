! $Id: drydep_info.h,v 1.1.1.1 2009/07/24 20:34:55 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: drydep_info.h
!
! !DESCRIPTION: Header file containing parameters and arrays for the
!  GEOS-Chem dry deposition routines.  These were taken from the input
!  files \texttt{drydep.coef} and \texttt{drydep.table}.
!\\
!\\
! !REMARKS:
!  NOTE: You should make sure that these values match up with those in 
!  "drydep.coef" and "drydep.table".  However, these values are unlikely 
!  to change anytime soon.
!
!  Dry deposition land types:
!  -------------------------------------
!  1  Snow/Ice (Wesely) - listed first.
!  2  Deciduous forest (Wesely)
!  3  Coniferous forest (Wesely)
!  4  Agricultural land (Wesely)
!  5  Shrub/grassland (Wesely)
!  6  Amazon forest (Jacob & Wofsy, JGR 1990)
!  7  Tundra (Jacob et al., JGR 1992)
!  8  Desert (Wesely)
!  9  Wetland (Wesely)
!  10 Urban (Wesely)
!  11 Water (Wesely)
! 
! !DEFINED PARAMETERS:
! 
      INTEGER, PARAMETER :: N_DLAND   = 11   ! # of drydep land types
      INTEGER, PARAMETER :: N_POLY    = 20   ! # of Baldocchi coefficients
      INTEGER, PARAMETER :: N_VEGTYPE = 74   ! # of Olson land types
      INTEGER, PARAMETER :: N_WATER   = 6    ! # of Olson types that are water
!
! !PUBLIC DATA MEMBERS:
! 
      !-----------------------------------------------------------------
      ! Baldocchi dry deposition coefficients (from "drydep.coef")
      !-----------------------------------------------------------------
      REAL*8, PARAMETER  :: DRYCOEFF(N_POLY) = 
     &     (/ -3.58d-01,  3.02d+00,  3.85d+00, -9.78d-02, -3.66d+00,  
     &         1.20d+01,  2.52d-01, -7.80d+00,  2.26d-01,  2.74d-01,  
     &         1.14d+00, -2.19d+00,  2.61d-01, -4.62d+00,  6.85d-01,
     &        -2.54d-01,  4.37d+00, -2.66d-01, -1.59d-01, -2.06d-01 /)

      !-----------------------------------------------------------------
      ! Resistances (R*) and max drydep velocities (Vsmax) for each of 
      ! the dry deposition land types (from "drydep.table")
      !-----------------------------------------------------------------

      ! Ri 
      INTEGER, PARAMETER :: IRI(N_DLAND)    = (/ 9999,  200,  400,  200,  
     &                                            200,  200,  200, 9999,  
     &                                            200, 9999, 9999      
     &                                        /) 
      ! Rlu    
      INTEGER, PARAMETER :: IRLU(N_DLAND)   = (/ 9999, 9000, 9000, 9000, 
     &                                           9000, 1000, 4000, 9999, 
     &                                           9000, 9999, 9999 
     &                                        /) 
              
      ! Rac   
      INTEGER, PARAMETER :: IRAC(N_DLAND)   = (/    0, 2000, 2000,  200,  
     &                                            100, 2000,    0,    0,  
     &                                            300,  100,    0 
     &                                        /) 
              
      ! Rgss  
      INTEGER, PARAMETER :: IRGSS(N_DLAND)  = (/  100,  500,  500,  150,  
     &                                            350,  200,  340, 1000,    
     &                                              0,  400,    0      
     &                                        /)
              
      ! Rgso  
      INTEGER, PARAMETER :: IRGSO(N_DLAND)  = (/ 3500,  200,  200,  150,  
     &                                            200,  200,  340,  400, 
     &                                           1000,  300, 2000      
     &                                        /)
              
      ! Rcls 
      INTEGER, PARAMETER :: IRCLS(N_DLAND)  = (/ 9999, 2000, 2000, 2000, 
     &                                           2000, 9999, 9999, 9999, 
     &                                           2500, 9999, 9999      
     &                                        /)
              
      ! Rclo  
      INTEGER, PARAMETER :: IRCLO(N_DLAND)  = (/ 1000, 1000, 1000, 1000, 
     &                                           1000, 9999, 9999, 9999, 
     &                                           1000, 9999, 9999      
     &                                        /) 
   
      ! Vsmax
      INTEGER, PARAMETER :: IVSMAX(N_DLAND) = (/  100,  100,  100,  100,  
     &                                            100,  100,  100,   10,  
     &                                            100,  100,   10 
     &                                        /)
    
      !-----------------------------------------------------------------
      ! Mapping between dry deposition land types <--> Olson land types
      ! (from "drydep.table")
      !-----------------------------------------------------------------

      ! Dry deposition landtypes corresponding to Olson landtypes
      INTEGER, PARAMETER :: IDEP(N_VEGTYPE) = (/ 11, 10,  5,  1,  1,  1,  
     &                                            2,  1,  8,  1,  1,  1,  
     &                                            1,  1,  1,  1,  5,  1,  
     &                                            1,  1,  3,  3,  3,  3,  
     &                                            2,  2,  2,  3,  2,  2,  
     &                                            4,  4,  2,  6,  1,  1,  
     &                                            9,  4,  4,  4,  5,  5,  
     &                                            5,  5,  5,  9,  5,  5,  
     &                                            5,  5,  8,  8,  5,  7,  
     &                                            6,  2,  2,  2,  2,  2,  
     &                                            3,  3,  3,  5,  5, 11, 
     &                                           11, 11, 11,  8,  1,  8,  
     &                                            9, 11 /)

      ! Indices of Olson Land Types that are water
      INTEGER, PARAMETER  :: IWATER(N_WATER) = (/  1, 66, 67, 
     &                                            68, 69, 74 /) 
!
! !REVISION HISTORY:
!  24 Jun 2009 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
