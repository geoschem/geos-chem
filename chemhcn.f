! $Id: chemhcn.f,v 1.3 2004/09/21 18:04:09 bmy Exp $
      SUBROUTINE CHEMHCN
!
!******************************************************************************
!  Subroutine CHEMHCN performs HCN chemistry. Loss is via reaction with OH, 
!  O(1D), photolysis and ocean uptake. (qli, bmy, 12/15/98, 7/20/04)
!
!  NOTES:
!  (1 ) Now use F90 syntax (bmy, 3/24/99)
!  (2 ) Now use double precision exponents, e.g. 1d0.  (bmy, 3/18/99)
!  (3 ) Add the stratospheric monthly mean O1D and OH fields from Hans
!        Schneider.  (qli, 3/31/99) 
!  (4 ) Add the simple photolysis loss mechenism of HCN.  (qli, 4/5/99)
!  (5 ) Now trap I/O errors (bmy, 4/12/99)
!  (6 ) Using SLOW-J in calculating the photolysis of HCN (bmy, qli, 4/16/99)
!  (7 ) Include the dependence of HCN air-sea transfer velocity on the Schmidt
!        number (qli, 4/16/99)
!  (8 ) Now use subroutine IOERROR to trap I/O errors (bmy, 5/27/99)
!  (9 ) OH is now passed back from "fillfields.f" via the argument list.
!       This saves a little memory in common blocks. (bmy, 10/20/99)
!  (10) Rearrange some comments (bmy, 10/20/99)
!  (11) AD22 is now declared allocatable in "diag_mod.f". (bmy, 11/29/99)
!  (12) ALBD is not used anymore; remove it from the argument list.
!        Also remove old code for trapping I/O errors.(bmy, 2/2/00)
!  (13) Add C-preprocessor switch LSLOWJ to bracket code for 
!        SLOW-J photolysis (bmy, 2/25/00)
!  (14) LTJV is now declared as allocatable in "diag_mod.f" (bmy, 4/10/00)
!  (15) Now reference AIRVOL, TS, U10M, and V10M from "dao_mod.f"
!        instead of from common block header files. (bmy, 6/23/00)
!  (16) Now reference the monthly mean OH array and the routine which reads
!        it from disk in "global_oh_mod.f" (bmy, 7/28/00)
!  (17) Removed obsolete code from 7/28/00 (bmy, 8/31/00)
!  (18) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (19) Now references AVGW from "dao_mod.f" (bmy, 9/24/01)
!  (20) Need to reference JLOP from "comode_mod.f" (bmy, 10/3/01)
!  (21) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (22) Now reference IU_FILE and IOERROR from "file_mod.f" (bmy, 6/26/02)
!  (23) Bug fixes for SLOW-J: reference ALBD from "dao_mod.f", and also
!        update argument list for the call to "ruralbox.f" (bmy, 7/31/02)
!  (24) Remove P and SIG from arg list to RURALBOX.  Deleted obsolete,
!        commented-out code from 6/02 and 7/02.  Updated comments, cosmetic
!        changes. (dsa, bdf, bmy, 9/23/02)
!  (25) Now reference AD, SUNCOS, T from "dao_mod.f".  Now make FIRSTCHEM a
!        local SAVEd variable.  Now set LMN = MONTH from "CMN". (bmy, 11/15/02)
!  (26) Now replace DXYP(J) with routine GET_AREA_M2 from "grid_mod.f"
!        Now use functions GET_MONTH, GET_TS_CHEM from "time_mod.f".
!        (bmy, 2/11/03)
!  (27) TWO_PI is now obsolete for SMVGEAR II (bdf, bmy, 4/1/03)
!  (28) Now references N_TRACERS and STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : AD,      AIRVOL, AVGW, ALBD, 
     &                          SUNCOS,  T,      TS,   U10M, V10M
      USE DIAG_MOD,      ONLY : AD22,    LTJV
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE COMODE_MOD,    ONLY : JLOP
      USE GLOBAL_OH_MOD, ONLY : OH,      GET_GLOBAL_OH
      USE GRID_MOD,      ONLY : GET_AREA_M2
      USE TIME_MOD,      ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACER_MOD,    ONLY : N_TRACERS, STT

      IMPLICIT NONE

#     include "CMN_SIZE"
!------------------------------
! Prior to 7/20/04:
!#     include "CMN"       
!------------------------------
#     include "CMN_DEP" 
#     include "CMN_HCN" 
#     include "CMN_DIAG"  ! For J-Value diagnostic
#     include "hcn.h"     ! Switches for HCN simulation

#if   defined( LSLOWJ )
#     include "comsol.h"
#     include "comode.h"
#endif
      
      ! Local variables
      LOGICAL, SAVE     :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE     :: LMN_LAST, IOS 
      INTEGER           :: I, J, K, L, N, J1, J2, LMN
      REAL*8            :: DTCHEM, RDLOSS, T1L  
      REAL*8            :: K0, K1, AIR_DENS 
      REAL*8            :: TEMP_11, TEMP_12, TEMP_1
      REAL*8            :: TEMP_21, TEMP_22, TEMP_2

      ! RC_HCN_OH : rate constant for the reaction of HCN + OH,
      ! (cm^3/molecule/s)
      REAL*8            :: RC_HCN_OH

      ! RC_HCN_O1D: rate constant for the reaction of HCN + O1D,
      ! (cm^3/molecule/s) taken from Cicerone(1983), the value is assumed
      REAL*8, PARAMETER :: RC_HCN_O1D = 1.0d-10 

      ! Some variables used in calculating HCN uptake by ocean, for their 
      ! meanings, see the explanation to the related code    
      REAL*8            :: H, H_, u, Tc, Sc, kl, kg, KH, KKl, Cg, F
      REAL*8, PARAMETER :: fu = 2.778d-6

      ! XNUMOL_AIR: (molecules air / kg air) 
      REAL*8, PARAMETER :: XNUMOL_AIR = 6.022d+23 / 0.02897d0 

      ! Rg: Universal gas constant, (m^3 atm / K / mol)
      REAL*8, PARAMETER :: Rg = 8.2d-5

      ! External functions
      REAL*8, EXTERNAL  :: BOXVL
                    
      ! Coefficients for fitting the Schmdit number 
      REAL*8, PARAMETER :: A0 =  1888.382d0
      REAL*8, PARAMETER :: A1 = -78.24071d0
      REAL*8, PARAMETER :: A2 =  1.267327d0
      REAL*8, PARAMETER :: A3 = -0.00810245d0

      ! Coefficients for fitting Henry's law constant of HCN, 
      ! taken from Edwards(1978)
      REAL*8, PARAMETER :: B1 = -49068.8d0
      REAL*8, PARAMETER :: B2 = -241.82d0 
      REAL*8, PARAMETER :: B3 =  0.315014d0  
      REAL*8, PARAMETER :: B4 =  1446.005d0

#if   defined( LSLOWJ )
      ! Necessary variables for slow-J photolysis (qli, bmy, 4/15/99)
      INTEGER      :: IJLOOP, IJWINDOW, NPTS
      INTEGER      :: IDXAIR(JJPAR),IDXO3(JJPAR)
     
      REAL*8       :: ALT(MAXIJ,LLPAR)
      REAL*8       :: SURFALT(MAXIJ)
      REAL*8       :: TOTO3(JJPAR)
      REAL*8       :: CLOUDS(MAXIJ,11)
      REAL*8       :: ROVMG, JVALUE 
      REAL*8       :: CLOUDREF(5)
      REAL*8, SAVE :: XSECT_HCN(MXWL)
#endif
     
      !=================================================================
      ! CHEMHCN begins here!
      !=================================================================

      ! Month
      LMN = GET_MONTH()

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      !=================================================================
      ! Do the following on the first chemistry timestep
      !=================================================================
      IF ( FIRSTCHEM ) THEN
         LMN_LAST  = 99

         ! Read the stratospheric OH field for levels 15 to 20, 
         ! totally 12 months, into STRAT_OH. NSKIPL is 14, 
         ! defined in CMN -- qli,3/31/99
         OPEN( IU_FILE, FILE = "OH.strat", IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'chemhcn:1' )
         
         DO I = 1, 12*JJPAR
            READ( IU_FILE, *, IOSTAT=IOS ) 
     &           ( STRAT_OH(I,J), J = 1,(LLPAR - LLTROP) )
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'chemhcn:2' )

         ENDDO 
  
         CLOSE( IU_FILE )

         ! Read the stratospheric O1D monthly mean field, 
         ! totally 12 months, into STRAT_O1D.  (qli,3/31/99)
         OPEN( IU_FILE, FILE = "O1D.strat", IOSTAT=IOS )
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'chemhcn:3' )

         DO I = 1, 12*JJPAR
            READ( IU_FILE, *, IOSTAT=IOS ) 
     &           ( STRAT_O1D(I,J), J = 1,(LLPAR - LLTROP + 2) )
            IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'chemhcn:4' )
         ENDDO 

         CLOSE( IU_FILE )

#if   defined( LSLOWJ )
         ! Read O3 & T profiles for SLOW-J photolysis (qli, bmy, 4/15/99)
         CALL JVALUEIN
         CALL READCSHCN( XSECT_HCN )
#endif

         ! Set FIRSTCHEM to false (bmy, 4/16/99)
         FIRSTCHEM = .FALSE.
      ENDIF

      !=================================================================
      ! Read in the tropospheric OH fields of the (LMN)th month into OH
      !=================================================================
      IF ( LMN /= LMN_LAST ) THEN 
         CALL GET_GLOBAL_OH( LMN )

         ! Copy from STRAT_OH  the stratospheric OH field into OH, 
         ! and copy from STRAT_O1D  the stratospheric O1D field into O1D. 
         ! ( qli, 3/31/99 ) 
         DO I = 1, IIPAR
            J1 = (LMN-1)*JJPAR + 1
            J2 = LMN*JJPAR
           
            OH(I,:,(LLTROP+1):LLPAR) = STRAT_OH(J1:J2,:)
            
            O1D(I,:,(LLTROP-1):LLPAR) = STRAT_O1D(J1:J2,:)  

         ENDDO
         
         LMN_LAST = LMN       
      ENDIF

      !=================================================================
      ! (1) HCN Loss due to chemical reaction with OH
      ! --------------------------------------------
      ! Calculate the rate constant for the reaction of HCN + OH, which 
      ! depends on height, i.e., varying total air number density, and 
      ! temperature in each grid box.
      !               K0 * [M]   {1 + [LOG(K0 * [M] / K1)] ^ 2 } ^ (-1) 
      !     K  =  ------------------- * 0.8
      !           1 + K0 * [M] / K1
      !
      !     K0 = 1.50E-31 * EXP( - 875 / T )
      !
      !     K1 = 1.20E-13 * EXP( - 400 / T )
      !
      ! where 
      !     K      : rate constant
      !     T      : temperature
      !     [M]    : total atmospheric density(cm^-3)
      !     K0, K1 : experimentally determined rate coefficients in the
      !               limit of low and high pressure, respectively.
      ! Reference  : Cicerone, R.J., and R. Zellner, The Atmospheric 
      !               Chemistry of Hydrogen Cyanide (HCN), J. Geophys. 
      !               Res., 88, 10,689-10,696, 1983       
      !
      ! -- qli, 1/24/1999
      !=================================================================
      N    = 1      ! Only one tracer, HCN
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         K0        = 1.50d-31  * EXP( - 875.d0 / T(I,J,L) )  

         K1        = 1.20d-13  * EXP( - 400.d0 / T(I,J,L) ) 

         ! AIR_DENS: density of air in grid box(I,J,L), in (molecule/cm^3); 
         AIR_DENS   = AD(I,J,L) / BOXVL(I,J,L) * XNUMOL_AIR

         TEMP_11    = K0  * AIR_DENS

         TEMP_12    = TEMP_11 / K1

         TEMP_1     = TEMP_11 / ( 1.d0 + TEMP_12 )

         TEMP_21    = LOG10( TEMP_12 )

         TEMP_22    = 1.d0 / ( 1.d0 + TEMP_21 * TEMP_21 )

         ! To save time: LOG(0.8) = -0.2231
         TEMP_2     = EXP( -0.2231d0 * TEMP_22 )

         RC_HCN_OH = TEMP_1 * TEMP_2
         
         RDLOSS    = RC_HCN_OH * OH(I,J,L) * DTCHEM
        
         T1L          = STT(I,J,L,N) * RDLOSS
         STT(I,J,L,N) = STT(I,J,L,N) - T1L
         HCN_LOSS(I,J,L,1) =   HCN_LOSS(I,J,L,1) + T1L
      
      ENDDO
      ENDDO
      ENDDO 
!
!*****************************************************************************
!(2) HCN Loss due to reaction with O(1D) 
!---------------------------------------
! The rate constant for HCN + O(1D)-> is assumed to be 1E-10 cm^3/molecule/s,
! independent of pressure and temperature and hence altitude.
!
! Reference  :  Cicerone, R.J., and R. Zellner, The Atmospheric Chemistry of
!               Hydrogen Cyanide (HCN), J. Geophys. Res., 88, 10,689-10,696,
!               1983       
!
! -- qli, 1/24/1999
!
! Note:  We take this loss into consideration only for levels 13 (=NSKIPL-2) 
! ====   to 20, because at levels lower than 13 this reaction is quite slow.
! 
! -- qli, 4/1/99
!*****************************************************************************
!   
      N    = 1
      DO L = LLTROP-1, LLPAR  
      DO J = 1, JJPAR
      DO I = 1, IIPAR
 
         RDLOSS       = RC_HCN_O1D   * O1D(I,J,L) * DTCHEM
        
         T1L          = STT(I,J,L,N) * RDLOSS
         STT(I,J,L,N) = STT(I,J,L,N) - T1L
         HCN_LOSS(I,J,L,2) =   HCN_LOSS(I,J,L,2) + T1L
     
      ENDDO
      ENDDO
      ENDDO
!
!*****************************************************************************
!(3) HCN Loss due to photolysis 
!-------------------------------
! The photolysis rate used here are those for HCl, because we do not have the 
! absorption cross sections for HCN at wavelength > 200 nm, which is the
! typical actinic flux wavelength in the lower stratopshere. We assumed that
! the photoabsorption and dissociation of HCN is similar to that of HCl.
!
!  The photo rate for each grid box is in (s^-1) so multiply this
!  by the number of seconds in the chemistry interval and use that
!  as the loss rate (i.e. the argument of the exponential).
!   
! -- bmy, qli, 4/15/99
!*****************************************************************************
!   
#ifdef HCNPHOTO

! Now we use C-preprocessor switch LSLOWJ (bmy, 2/25/00)
#if   defined( LSLOWJ )

      ! set up things from chemdr.f and ruralbox.f (qli, bmy, 4/15/99)
      NLONG   = IIPAR
      NLAT    = JJPAR
      NVERT   = IVERT
      NPVERT  = NVERT
      ! NIJLOOP = MAXIJ

      ! Get the altitude
      CALL GETALT  ( ALT, SURFALT, ROVMG  )

      ! Set up rural boxes
      NLOOP       = NLAT  * NLONG
      NTLOOP      = NLOOP * NVERT
 
      ! Call ruralbox to set up rural boxes
      CALL RURALBOX( AD,     T,      AVGW,   ALT,   ALBD,  
     &               SUNCOS, CLOUDS, LEMBED, IEBD1, IEBD2, 
     &               JEBD1,  JEBD2,  LPAUSE )

      ! Get the total ozone column
      CALL GETTOTO3( TOTO3, IDXAIR, IDXO3 )

      NTTLOOP     = NTLOOP
      NPTS        = NTTLOOP

      ! Call SOL to get the actinic flux (qli, bmy, 4/15/99)
      CALL SOL (NPTS,SUNCOS,ALT,SURFALT,TOTO3,CLOUDS,IDXAIR,IDXO3)

      ! Calculate the loss due to photolysis
      N = 1  ! Only one tracer, HCN 
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IJLOOP = JLOP(I,J,L)

         ! Only calculate J-values for daytime (bmy, 4/15/99)
         IF ( SUNCOS( (J-1)*IIPAR + I ) > 0 ) THEN
            JVALUE = 0d0
            
            DO K = 1, MXWL
               JVALUE = JVALUE + ( ACTFLX(K,IJLOOP) * XSECT_HCN(K) )
            ENDDO

            RDLOSS       = EXP( -JVALUE * DTCHEM )
            HCN_LOSS(I,J,L,3) = HCN_LOSS(I,J,L,3) + STT(I,J,L,N) *
     &                          (1d0 - RDLOSS)
            STT(I,J,L,N) = STT(I,J,L,N) * RDLOSS

            ! Add ND22 - J-Value diagnostic
            IF ( ND22 > 0 ) THEN
               IF ( LTJV(I,J) > 0 .and. L <= LD22 ) THEN
                  AD22(I,J,L,N) = AD22(I,J,L,N) + JVALUE
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
#endif
#endif

!
!***************************************************************************
! (4)HCN Loss due to ocean uptake 
! -------------------------------
! Which happens only in grid boxes in the first level, i.e., L = 1.
! The flux of gas (HCN) through the air-sea interface can be expressed as:
!
!  F = KKl * ( Cg / H - Cl)
!
! where
! (a)KKl    :  total transfer velocity expressed on a liquid phase (seawater) 
!               basis, 
!
!                 1      1        1
!               ----- = ---  + --------
!                KKl     kl     H * kg
!
!               where 
!                (a.1) kl : transfer velocity in liquid phase (seawater)
!
!                  Liss & Merlivat (1986):
!
!                  kl = fu * (0.17 * u) *(Sc/600)^(-2/3),       u < 3.6
!                       fu * (2.85 * u - 9.65)*(Sc/600)^(-1/2), 3.6 < u < 13
!                       fu * (5.90 * u - 49.3)*(Sc/600)^(-1/2), u > 13
!
!                  Wanninkhof (1992):
!                   
!                  kl = 0.31 *u*u *(Sc/660)^(-1/2)
!
!                  where u is the wind speed at the height of 10m. fu is the
!                  unit conversion factor(cm/hr to m/s). Sc is the Schmidt
!                  number for HCN:
!
!                    Sc = A1 * T^4 + A2 * T^3 + A3 * T^2 + A4 * T + A5
!
!                (a.2) kg : transfer velocity in gas phase
!                    
!                  kg = (5.2E-5 + 3.2E-3 * u) * SQRT(WTH2O / WTHCN)
!                  
!                  where WTH2O =18.0 and WTHCN =27.0 are the molecular weight 
!                  of H2O and HCN in (g/mole), respectively.
! 
! (b) Cg, Cl :  gas (HCN) concentration in gas and liquid phase, respectively 
! (c) H      :  dimensionless Henry's law constant
!
!                        1
!               H = ------------- 
!                    KH * Rg * T
!
!               where
!                (c.1) KH  :  Henry's law constant(moles / litr / atm) which 
!                             can be evaluated from the following equation:
!                      
!                      ln( H_ ) = B1 / T + B2 * ln(T) + B3 * T +B4 
!                      KH       = 1  / H_
!                       
!                (c.2) Rg  :  Universal gas constant
!                (c.3) T   :  Temperature
!
! (d)Reference:
!
!             Liss, P.S., and P.G. Slater, Flux of gases across the air-
!                sea interface, Nature, 247, 181-184, 1974
!             Liss, P.S., and L. Merlivat, Air-sea gas exchange rates: 
!                Introduction and synthesis, in The Role of Air-Sea Exchange 
!                in Geochemical Cycling, pp113-127, edited by P. Buat-Menard, 
!                D.Reidel, Norwell, MA, 1986
!             Edwards, T.J., G. Maurer, J. Newman, and J.M. Prausnitz, Vapor-
!                Liquid Eauilibria in Multicomponent Aqueous Solutions of
!                Volatile Weak Electrolytes, AIChE J., 24, 966-976, 1978    
!             Asher, W., The sea-surface microlayer and its effect on global
!                air-sea gas transfer, in The sea surface and global change,
!                edited by P.S. Liss and R.A. Duce, Cambridge Univ. Press, 
!                Cambridge, UK, 1997  
!             Wilke, C.R., and P. Chang, Correlation of diffusion coefficients
!                in dilute solutions, AIChE J., 1, 264-270, 1955
!             Jahne, B., G. Heinz, and W. Dietrich, Measurement of the 
!                diffusion coefficients of sparingly soluable gases in 
!                water with a modified Barrer method, J. Geophys. Res., 
!                92, 10,767-10,776, 1987                      
!             Wanninkhof,R., Relationship Between Wind Speed and Gas Exchnage 
!                Over the Ocean, J. Geophys. Res., 97, 7373-7382, 1992
!
! -- qli, 1/24/1999, 4/16/99
!*****************************************************************************
!
         N    = 1      ! Only one tracer, HCN
         L    = 1      ! Ocean uptake is effective in the 1st layer
         DO J = 1, JJPAR
         DO I = 1, IIPAR    

          ! TS(I,J) the surface temperature in grid box (I, J), in (K); 
          ! H_ is in ( kg atm / mole ); Convert KH from (mole/litr/atm) 
          ! to (mole/m^3/atm) by  the factor 1000; 
          H_  = B1 / TS(I,J) + B2 * LOG( TS(I,J) ) + B3 * TS(I,J) + B4
          
          H_  = EXP( H_ )
          
          KH  = 1000.d0 / H_ 
 
          H   = 1.d0 / ( KH * Rg * TS(I,J) )

           ! Wind speed at 10m can be evaluated from UV10M and V10M, in (m/s), 
           ! which are zonal and meridional wind at 10m, passed by CMN_UV10M; 
           u   = SQRT( U10M(I,J) * U10M(I,J) + V10M(I,J) * V10M(I,J) )

           ! Evaluate kl, kg, and KKl, in (m/s);
           ! Tc is the surface temperature in Celsius
           Tc = TS(I,J) - 273.15   

           Sc = A0 + A1* Tc + A2 * Tc**2 + A3 * Tc**3

!----------------------------------------------------------------------------
! Liss & Merlivat (1986) relationship 
!
#ifdef  LM86  
           IF ( u <= 3.6 ) THEN
             kl = fu * (0.17d0 * u         ) * (Sc/600.d0)**(-0.67d0)

           ELSEIF ( u > 13.0 ) THEN
             kl = fu * (5.90d0 * u - 49.3d0) * (Sc/600.d0)**(-0.50d0)

           ELSE
             kl = fu * (2.85d0 * u - 9.65d0) * (Sc/600.d0)**(-0.50d0)

           ENDIF

! Liss & Merlivat (1986) relationship enhanced by 1.6
#ifdef  EN_16 
           kl   = 1.6d0 * kl
#endif

#endif
!----------------------------------------------------------------------------
! WannIkohf (1992) relationship
!
#ifdef  W92   
           kl   = fu * 0.31 *u*u * (Sc/660.d0)**(-0.50d0)
#endif  
!----------------------------------------------------------------------------

           ! To save time: SQRT( 18. / 27. ) = 0.8165;
           kg    = (5.2d-5 + 3.2d-3 * u) * 0.8165d0 

           KKl   = kl * kg * H / ( kl + kg * H )

           ! Cg gas phase concentration of HCN, in (kg/m^3);
           Cg    = STT(I,J,L,N) / AIRVOL(I,J,L)  

           ! Assuming Cl, which is not available now, to be neglectible, 
           ! then F = KKl * Cg / H = KKl * Cg * KH * Rg * T, 
           ! invoking H = 1 / (KH * Rg * T);  F is in (kg/m^2/s);
           F     = KKl * KH * Rg * TS(I,J) * Cg 

           ! DXYP(J) is the Delta-X * Delta-Y surface area of grid 
           ! boxes (I,J,L=1),in (m^2); FRCLND(I,J) contains land fraction 
           ! of surface grid boxes (I,J,L=1)
           T1L   = F * GET_AREA_M2( J ) * (1.d0 - FRCLND(I,J)) * DTCHEM 

           STT(I,J,L,N) = STT(I,J,L,N) - T1L 
           HCN_LOSS(I,J,L,4) = HCN_LOSS(I,J,L,4) + T1L
 
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CHEMHCN



