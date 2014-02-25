!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ef_mgn20_mod
!
! !DESCRIPTION: Module EF\_MGN20 contains emission factors for 20 MEGAN
!  species.  The values in this file have to be in the same order as in
!  SPC\_MGN.EXT
!\\
!\\
! !INTERFACE: 
!
      MODULE EF_MGN20_MOD
!
! !USES:
!
      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS: 
!
      ! Number of categories
      INTEGER, PARAMETER :: N_EF_SPC = 20
!
! !PUBLIC DATA MEMBERS:
!
      ! Speciation factor for broadleaf
      REAL*8 ::      EF_BT( N_EF_SPC )

      ! Speciation factor for needleleaf
      REAL*8 ::      EF_NT( N_EF_SPC )

      ! Speciation factor for shrub
      REAL*8 ::      EF_SB( N_EF_SPC )

      ! Speciation factor for herb/grass/crop
      REAL*8 ::      EF_HB( N_EF_SPC )

      DATA    EF_BT(  1)  , EF_NT(  1)  , EF_SB(  1)  , EF_HB(  1)
     &      / 13000.00    , 2000.00     , 11000.00    , 400.00       /
      DATA    EF_BT(  2)  , EF_NT(  2)  , EF_SB(  2)  , EF_HB(  2)
     &      / 20.00       , 75.00       , 22.00       , 0.30         /
      DATA    EF_BT(  3)  , EF_NT(  3)  , EF_SB(  3)  , EF_HB(  3)
     &      / 45.00       , 70.00       , 50.00       , 0.70         /
      DATA    EF_BT(  4)  , EF_NT(  4)  , EF_SB(  4)  , EF_HB(  4)
     &      / 45.00       , 100.00      , 52.00       , 0.70         /
      DATA    EF_BT(  5)  , EF_NT(  5)  , EF_SB(  5)  , EF_HB(  5)
     &      / 18.00       , 160.00      , 25.00       , 0.30         /
      DATA    EF_BT(  6)  , EF_NT(  6)  , EF_SB(  6)  , EF_HB(  6)
     &      / 90.00       , 60.00       , 85.00       , 1.00         /
      DATA    EF_BT(  7)  , EF_NT(  7)  , EF_SB(  7)  , EF_HB(  7)
     &      / 90.00       , 300.00      , 100.00      , 1.50         /
      DATA    EF_BT(  8)  , EF_NT(  8)  , EF_SB(  8)  , EF_HB(  8)
     &      / 180.00      , 450.00      , 200.00      , 2.00         /
      DATA    EF_BT(  9)  , EF_NT(  9)  , EF_SB(  9)  , EF_HB(  9)
     &      / 90.00       , 180.00      , 110.00      , 4.80         /
      DATA    EF_BT( 10)  , EF_NT( 10)  , EF_SB( 10)  , EF_HB( 10)
     &      / 35.00       , 30.00       , 30.00       , 0.50         /
      DATA    EF_BT( 11)  , EF_NT( 11)  , EF_SB( 11)  , EF_HB( 11)
     &      / 30.00       , 60.00       , 45.00       , 0.90         /
      DATA    EF_BT( 12)  , EF_NT( 12)  , EF_SB( 12)  , EF_HB( 12)
     &      / 75.00       , 110.00      , 85.00       , 1.4          /
      DATA    EF_BT( 13)  , EF_NT( 13)  , EF_SB( 13)  , EF_HB( 13)
     &      / 0.10        , 100.00      , 1.00        , 0.01         /
      DATA    EF_BT( 14)  , EF_NT( 14)  , EF_SB( 14)  , EF_HB( 14)
     &      / 800.00      , 800.00      , 800.00      , 800.00       /
      DATA    EF_BT( 15)  , EF_NT( 15)  , EF_SB( 15)  , EF_HB( 15)
     &      / 240.00      , 240.00      , 240.00      , 80.00        /
      DATA    EF_BT( 16)  , EF_NT( 16)  , EF_SB( 16)  , EF_HB( 16)
     &      / 30.00       , 30.00       , 30.00       , 30.00        /
      DATA    EF_BT( 17)  , EF_NT( 17)  , EF_SB( 17)  , EF_HB( 17)
     &      / 5.00        , 6.00        , 30.00       , 70.00        /
      DATA    EF_BT( 18)  , EF_NT( 18)  , EF_SB( 18)  , EF_HB( 18)
     &      / 240.00      , 240.00      , 240.00      , 80.00        /
      DATA    EF_BT( 19)  , EF_NT( 19)  , EF_SB( 19)  , EF_HB( 19)
     &      / 70.00       , 70.00       , 70.00       , 70.00        /
      DATA    EF_BT( 20)  , EF_NT( 20)  , EF_SB( 20)  , EF_HB( 20)
     &      / 1000.0      , 1000.0      , 1000.0      , 1000.0       /
!
! !REMARKS:
!  MEGAN v2.02
!  INPUT version 210
!
! !REVISION HISTORY:
!  Tan          12/02/06 - Creates this file
!  Guenther A.  08/11/07 - Creates this file again with updates and move
!                          from v2.0 to v2.02
!  Havala Pye   07/26/10 - EFs now double precision (replace REAL with
!                          REAL*8) (part of SOA + semivolatile POA)
!  15 Mar 2012 - M. Payer  - Converted to module from include file.
!                            Added ProTeX headers.
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE EF_MGN20_MOD
!EOC
