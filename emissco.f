! $Id: emissco.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE EMISSCO( FIRSTEMISS, NSEASON, LMN, SUNCOS ) 
! 
!******************************************************************************
!  Subroutine EMISSCO reads CO emissions routine for the CO simulation
!  with parameterized OH.  EMISSCO is based on subroutine "emissdr.f".
!  (bey, 02/02/99, bnd & bmy, 9/13/00, 1/13/03)
!
!  NOTES:
!  (1 ) STT should be referenced (I,J,L,N) and not (IREF,JREF,L,N),
!        since it is of window dimension (IIPAR,JJPAR,LLPAR,NNPAR).
!  (2 ) Use the "D" exponent to explicitly force double precision.
!  (3 ) Now use subroutine IOERROR to trap I/O error conditions.
!  (4 ) TCVV_CO is now included in header file "CMN_CO_BUDGET" (bmy, 4/19/00)
!  (5 ) Added references to F90 module "biomass_mod.f" and "biofuel_mod.f".
!        Also no longer compute ND29 diagnostics here -- they are done in
!        BIOFUEL_BURN.  
!  (6 ) Now use IOS /= 0 to trap both errors and EOF. (bmy, 9/13/00)
!  (7 ) Removed obsolete code from 9/13/00 (bmy, 10/16/00)
!  (8 ) Enhance anthro, biomass, and biofuel CO emissions to account for
!        production of CO from VOC's (bnd, bmy, 1/3/01)
!  (9 ) Now make sure IDBCO is not zero -- to avoid subscript range
!        errors when indexing BURNEMIS.  Also make sure that IDBFCO is 
!        not zero, to avoid errors when indexing BIOFUEL. (bmy, 3/20/01)
!  (10) Now do not let SCALEYEAR go higher than 1996.  This is the last
!        year for which we have FF scale factor data. (bnd, bmy, 4/6/01)
!  (11) Eliminate obsolete commented-out code (bmy, 4/20/01)
!  (12) Now prompt user to check IDBCO and IDBFCO in "tracer.dat" if
!        these switches are turned off (bmy, 6/19/01)
!  (13) ND46 diagnostic is now archived as atoms C/cm2/s here, instead of
!        dividing by DTSRCE in "diag3.f" (bmy, 9/14/01)
!  (14) BIOFUEL (N,IREF,JREF) is now BIOFUEL(N,I,J).  BURNEMIS(N,IREF,JREF) is
!        now BURNEMIS(N,I,J).  (bmy, 9/28/01)
!  (15) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (16) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (17) Now reference IU_FILE and IOERROR from "file_mod.f".  Now use IU_FILE
!        as the file unit # instead of #65 and #66. (bmy, 6/26/02)
!  (18) Now reference BXHEIGHT from "dao_mod.f".  Also references IDBCO and
!        IDBFCO from "tracerid_mod.f".  Now do not let SCALEYEAR exceed 1998 
!        for the fossil fuel scaling. (bmy, 1/13/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BIOFUEL_MOD,  ONLY : BIOFUEL,  BIOFUEL_BURN
      USE BIOMASS_MOD,  ONLY : BURNEMIS, BIOBURN
      USE DIAG_MOD,     ONLY : AD29,     AD46
      USE FILE_MOD,     ONLY : IU_FILE,  IOERROR
      USE TRACERID_MOD, ONLY : IDBCO,    IDBFCO

      IMPLICIT NONE

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_O3"        ! EMIST, EMISRR, etc
#     include "CMN"           ! STT, many other arrays
#     include "CMN_DIAG"      ! Diagnostic switches
#     include "CMN_CO"        ! Arrays for CO parameterization
#     include "CMN_CO_BUDGET" ! TCO, XNUMOL_CO, TCVV_CO

      ! Arguments
      LOGICAL, INTENT(INOUT) :: FIRSTEMISS
      INTEGER, INTENT(IN)    :: LMN, NSEASON
      REAL*8,  INTENT(IN)    :: SUNCOS(MAXIJ)

      ! Local Variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      
      INTEGER, SAVE          :: LASTYEAR
      INTEGER                :: IREF, JREF, IJLOOP, SCALEYEAR
      INTEGER                :: I, J, L, N
      INTEGER                :: IOS
cbnd
      INTEGER                :: IPICK, CTM_lat, CTM_lon, NVLOOP
      REAL*8                 :: SCALEITBB, SCALEITFF,CONVERT_lon,REDISTR
      REAL*8                 :: SCALEITBF,COAfr
      REAL*8                 :: BOXPERCENT(10),TBOXPERCENT
cbnd
      REAL*8                 :: DTSRCE, BIO_CO, WOOD_CO
      REAL*8                 :: SCALE,  MOLRAT, BXHEIGHT_CM
      REAL*8                 :: TMMP,   EMXX,   EMMO
      REAL*8                 :: CONVERT(NVEGTYPE) 
      REAL*8                 :: GMONOT(NVEGTYPE)
      REAL*8                 :: EMX

      CHARACTER(LEN=4 )      :: CYEAR
      CHARACTER(LEN=25)      :: LIQCO2_FILE

      ! XNUMOL_C = atoms C / kg C
      REAL*8, PARAMETER      :: XNUMOL_C = 6.022d23 / 12d-3 
 
      ! External functions
      REAL*8, EXTERNAL       :: BOXVL
      REAL*8, EXTERNAL       :: XLTMMP
      REAL*8, EXTERNAL       :: EMISOP
      REAL*8, EXTERNAL       :: EMMONOT
!
cbnd
cbnd To increase/decrease biomass burning emissions and/or fossil fuel
cbnd emissions:
cbmy
cbmy The enhancement of biomass burning CO and biofuel burning CO is
cbmy now done in routines BIOBURN and BIOFUEL_BURN, respectively.
cbmy The enhancement needs to be done in these routines, in order
cbmy that the ND29 diagnostic archives the correct CO emissions.
cbmy
cbmy Therefore, in this routine, just increase CO by 20%.
cbmy Set SCALEITBB and SCALEITBF to 1d0, since the enhancement of 
cbmy BB CO and BF CO emissions has already been done. (bmy, 1/3/01)
cbmy

cbmy      SCALEITBB=1.1
cbmy      SCALEITFF=1.2
cbmy      SCALEITBF=1.1

      SCALEITBB = 1.00d0
      SCALEITFF = 1.20d0
      SCALEITBF = 1.00d0
cbnd
!
!******************************************************************************
!  Once per year, read in fossil fuel scale factors (these are ratios
!  of SCALEYEAR's emissions to 1985 emissions).
!  (1) FLIQCO2 = liquid fuel CO2 scale factors...used for CO, Hydrocarbons
!
!  If FSCALYR < 0 then use this year (JYEAR) for the scaling factors
!  Otherwise, use the value of FSCALYR as specified in 'input.ctm'.
!
!  If FSCALYR = 1985, then don't read the scale factors in, since they would
!  all be equal to 1 anyway.  This saves on extra multiplications, etc.
!******************************************************************************
!
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = JYEAR

         ! Put a cap on scaleyear for 1998 for now (bmy, 1/13/03)
         IF ( SCALEYEAR > 1998 ) SCALEYEAR = 1998
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF

      WRITE ( CYEAR, '(I4)' ) SCALEYEAR

      !-----------------------------------------------------------------------
      ! Debug output
      !print*,'ENTERING EMISSCO'
      !print*,'JYEAR=',JYEAR
      !print*,'SCALEYEAR=',SCALEYEAR
      !print*,'LASTYEAR=',LASTYEAR
      !print*,'FIRST=',FIRST
      !print*,'CYEAR=',CYEAR
      !print*,'FSCALYR=',FSCALYR
      !-----------------------------------------------------------------------

      IF ( FIRST .or. SCALEYEAR .ne. LASTYEAR ) THEN

         LIQCO2_FILE = 'scalefoss.liq.4x5.' // CYEAR

         IF ( SCALEYEAR .ne. 1985 ) THEN
            ! Read in data from the LIQUID fuel file scale factor 
            ! as REAL*4 data
            OPEN( IU_FILE, FILE = LIQCO2_FILE,  STATUS ='OLD',    
     &                     FORM ='UNFORMATTED', IOSTAT=IOS )
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissco:1' )

            ! Read in IREF, JREF dimensions
            READ( IU_FILE, IOSTAT=IOS ) IREF, JREF
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissco:2' )


            ! Read the data block
            READ( IU_FILE, IOSTAT=IOS ) ( ( FLIQCO2(I,J), I=1,IREF ), 
     &                                                    J=1,JREF )
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissco:3' )
            CLOSE( IU_FILE )
         ENDIF

         ! Set variables for next iteration
         FIRST    = .FALSE.
         LASTYEAR = SCALEYEAR
      ENDIF
!
!******************************************************************************
!  Echo JYEAR and NSEASON to the standard output
!******************************************************************************
!
!Comment out for now (bmy, 4/19/00)
!      WRITE ( 6, '(''--------------------------------------------'')' )
!      WRITE ( 6, '(''EMISSCO: NYEAR, NSEASON = '', i4,1x,i2)' )
!     &   JYEAR, NSEASON
!      WRITE (6, '(''Fossil Fuel Scale Year: '', a4) ') CYEAR
!
!******************************************************************************
!  Do the following on the very first time through EMISSCO
!******************************************************************************
!
      IF ( FIRSTEMISS ) THEN 

         ! Call some subroutines to set up ISOP emission (first time only!)
         CALL RDLIGHT
         CALL RDISOPT ( CONVERT )
         CALL RDMONOT ( GMONOT  )
         CALL SETBASE ( CONVERT, GMONOT )

         ! Initialize isoprene and montoterpene arrays
         SUMISOPCO(:,:) = 0.d0
         SUMMONOCO(:,:) = 0.d0

         ! Open the merge file
         WRITE( 6, '(a)' ) 'READING merge.4x5_CTM'
         
         ! NOTE: should upgrade to New fossil fuel file!
         OPEN( IU_FILE,  FILE='merge.4x5_CTM', STATUS='OLD', 
     &                   FORM='FORMATTED',     IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissco:4' ) 
         
         ! Read the fossil fuel emissions from the merge file
         READ( IU_FILE, '(7e10.3)', IOSTAT=IOS )
     &        EMISTNOX,  EMISTCO,   EMISTETHE, EMISTPRPE, EMISTC2H6,
     &        EMISTC3H8, EMISTALK4, EMISTACET, EMISTMEK,  EMISTSOX
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissco:5' )

         ! Close the merge file
         CLOSE ( IU_FILE )

         ! Set FIRSTEMISS = FALSE so we don't read this in again
         FIRSTEMISS = .FALSE.
      ENDIF
!
!******************************************************************************
!  NSRCE  =     the number of minutes between source events (GEOS-CTM)
!           OR  the number of hours   between source events (GISS-CTM)
!           AND is read in input.ctm (unit 5)
!  DTSRCE = NSRCE * 60 (GEOS-CTM) OR NSRCE * 3600 (GISS-CTM)
!******************************************************************************
!
      DTSRCE = 60.0d0 * NSRCE
!
!******************************************************************************
!  In the following, we modify the STT with emissions rates.
!
!  Anthropogenic CO emissions.
!  EMISTCO(IREF,JREF) in [molec CO/cm2/s]
!
!  DXYP(JREF)       is the surface area of the grid box in m^2.
!  DXYP(JREF) * 1e4 is the surface area of the grid box in cm^2.
!******************************************************************************
!
      ! Baseline year -- 1985
      IF ( SCALEYEAR == 1985 ) THEN
         DO J = 1, JJPAR
            JREF = J + J0
            DO I = 1, IIPAR
               IREF = I + I0
               
               ! Convert from [molec CO/cm2/s] to [kg CO] and store in STT
               STT(I,J,1,1) = STT(I,J,1,1) + 
     &                        ( EMISTCO(IREF,JREF) / XNUMOL_CO * 
cbnd     &                          DXYP(JREF) * 1d4 * DTSRCE )
     &                       DXYP(JREF) * 1d4 * DTSRCE )*SCALEITFF
cbnd

            ENDDO
         ENDDO

      ! Non-baseline year
      ELSE
         DO J = 1, JJPAR
            JREF = J + J0
            DO I = 1, IIPAR
               IREF = I + I0

               ! Convert from [molec CO/cm2/s] to [kg CO] and store in STT
               ! Also multiply by the FLIQCO2 scale factor array 
               STT(I,J,1,1) = STT(I,J,1,1) + 
     &                        ( EMISTCO(IREF,JREF) / XNUMOL_CO *
     &                          DXYP(JREF)         * 1d4       * 
cbnd     &                          FLIQCO2(IREF,JREF) * DTSRCE )
     &                       FLIQCO2(IREF,JREF) * DTSRCE )*SCALEITFF
cbnd

            ENDDO
         ENDDO
      ENDIF
!
!******************************************************************************
!  Biomass burning CO emissions 
!  stored in BURNEMIS(IDBCO,I,J) in [molec/cm3/s]
!******************************************************************************
!
      ! Make sure IDBCO > 0 
      IF ( IDBCO == 0 ) THEN
         WRITE( 6, '(a)' ) 'Biomass CO has been turned off!'
         WRITE( 6, '(a)' ) 'Check IDBCO in your "tracer.dat" file!'
         WRITE( 6, '(a)' ) 'STOP in emissco.f!'
         STOP
      ENDIF

      ! Read biomass burning emissions
      CALL BIOBURN( LMN )

      DO J = 1, JJPAR
         JREF = J + J0
         DO I = 1, IIPAR
            IREF = I + I0            
cbnd
            IPICK=0
cbnd
cbnd correct low emission factor for CO for N. Africa savannas.
            COAfr=1.
       IF(I.GE.32.AND.I.LE.50.AND.J.GE.24.AND.J.LE.29) COAfr=73./45.
cbnd

      !=================================================================
      ! Convert to degrees (1.25x1 degree grid (lon,lat)).
      !=================================================================
#if   defined( GRID4x5 )
      CONVERT_lon = ( DBLE(I) * 5.d0 ) * 1.d0 / 1.25d0
      CTM_lon     = INT( CONVERT_lon )
      CTM_lat     = ( J * 4 ) - 2

      IF(J.EQ.1    ) CTM_LAT = 2
      IF(J.EQ.JJPAR) CTM_LAT = 178

#elif defined( GRID2x25 )
      CONVERT_lon = ( DBLE(I) * 2.5d0 ) * 1.d0 / 1.25d0
      CTM_lon     = INT( CONVERT_lon )
      CTM_lat     = ( J * 2 ) - 1

      IF(J.EQ.1)     CTM_lat=1
      IF(J.EQ.JJPAR) CTM_lat=179

#elif defined( GRID1x1 )
      PRINT*, 'Need to compute CONVERT_LON for 1 x 1 grid!'
      PRINT*, 'STOP in TOMSAI (biomass_mod.f)'
      STOP

#endif

cbnd For the following regions, distribute biomass burning emissions into
cbnd the first 10 layers.  Most fires in these regions are intense creating
cbnd strong convection and lofting emissions into the free troposphere.

      ! Indonesia
      IF(CTM_lat.GE.83.and.CTM_lat.LE.99) THEN
         IF(CTM_lon.GE.221.and.CTM_lon.LE.269) THEN
            IPICK=1
         ENDIF
      ENDIF

      ! Canada and Alaska
      ! We have fire burn estimates for Canada, so we can use
      ! this data to fill in the TOMS data gap.
      IF(LMN.GE.5.AND.LMN.LE.9) THEN
         IF(CTM_lat.GE.141.and.CTM_lat.LE.161) THEN
            IF(CTM_lon.GE.16.and.CTM_lon.LE.96) THEN
               IPICK=1
            ENDIF
         ENDIF
      ENDIF

      ! Asiatic Russia
      IF(LMN.GE.5.AND.LMN.LE.9) THEN
         IF(CTM_lat.GE.136.and.CTM_lat.LE.161) THEN
            IF(CTM_lon.GE.211.and.CTM_lon.LE.291) THEN
               IPICK=1
            ENDIF
         ENDIF
      ENDIF
cbnd
cbnd Convert from [molec/cm3/s] to [kg] and store in STT.
cbnd Emit emissions only into layer 1.

      IF(IPICK.EQ.0) THEN
            STT(I,J,1,1) = STT(I,J,1,1) +
     &                     ( BURNEMIS(IDBCO,I,J) / XNUMOL_CO * 
cbnd     &                       BOXVL(I,J,1) * DTSRCE )
     &                       BOXVL(I,J,1) * DTSRCE )*SCALEITBB*COAfr

cbnd
      ELSE
cbnd Emit emissions into layers 1-10 assuming strong convection columns
cbnd  generated by forest fires.
cbnd  REDISTR ( molec / box 1 )
          REDISTR=BURNEMIS(IDBCO,I,J)*BOXVL(I,J,1)*DTSRCE
cbnd
cbnd Calculate the percent in each box of total volume in first 10 boxes. 
            TBOXPERCENT=0.      
          DO NVLOOP=1,10
            TBOXPERCENT=TBOXPERCENT+BOXVL(I,J,NVLOOP)      
          ENDDO
          DO NVLOOP=1,10
            BOXPERCENT(NVLOOP)=0.
            BOXPERCENT(NVLOOP)=BOXVL(I,J,NVLOOP)/TBOXPERCENT
          ENDDO
cbnd
cbnd Convert to kg.
          DO NVLOOP=1,10
            STT(I,J,NVLOOP,1) = STT(I,J,NVLOOP,1) + REDISTR *
     &        ( BOXPERCENT(NVLOOP) / XNUMOL_CO ) * SCALEITBB 
          ENDDO
cbnd
      ENDIF
            
         ENDDO
      ENDDO
!
!*****************************************************************************
!  CO Biofuel Emissions 
!  stored in BIOFUEL(N,IREF,JREF) in [molec/cm3/s]
!*****************************************************************************
!
      ! Make sure IDBFCO > 0 
      IF ( IDBFCO == 0 ) THEN
         WRITE( 6, '(a)' ) 'Biofuel CO has been turned off!'
         WRITE( 6, '(a)' ) 'Check IDBFCO in your "tracer.dat" file!'
         WRITE( 6, '(a)' ) 'STOP in emissco.f!'
         STOP
      ENDIF

      ! Read biofuel burning emissions
      CALL BIOFUEL_BURN

      DO J = 1, JJPAR
      DO I = 1, IIPAR
         STT(I,J,1,1) = STT(I,J,1,1) +
     &                   ( BIOFUEL(IDBFCO,I,J) / XNUMOL_CO * 
cbnd     &                       BOXVL(I,J,1)         * DTSRCE )
     &                       BOXVL(I,J,1)         * DTSRCE )*SCALEITBF
cbnd
      ENDDO
      ENDDO
!
!*****************************************************************************
!  ND29 diagnostics...save anthro CO, biomass burning CO, and wood 
!  burning CO in [molec CO/cm2/s].
!
!  NOTE: Now use AD29 array instead of AIJ (bmy, 4/16/00)
!*****************************************************************************
!
      IF ( ND29 > 0 ) THEN

         DO J = 1, JJPAR
            JREF = J + J0
            DO I = 1, IIPAR
               IREF = I + I0     

               ! BXHEIGHT is in meters, convert to cm
               BXHEIGHT_CM = BXHEIGHT(I,J,1) * 100d0

               ! Anthropogenic CO
               IF ( SCALEYEAR == 1985 ) THEN
cbnd                  AD29(I,J,1) = AD29(I,J,1) + EMISTCO(IREF,JREF)
               AD29(I,J,1) = AD29(I,J,1) + EMISTCO(IREF,JREF)*SCALEITFF
cbnd
               ELSE
                  AD29(I,J,1) = AD29(I,J,1) + ( EMISTCO(IREF,JREF) * 
cbnd     &                                          FLIQCO2(IREF,JREF) )
     &                         FLIQCO2(IREF,JREF) )*SCALEITFF
cbnd
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!*****************************************************************************
!  CO budget information
!*****************************************************************************
!
      DO J = 1, JJPAR
         JREF = J + J0
         DO I = 1, IIPAR
            IREF = I + I0
  
cbnd correct low emission factor for CO for N. Africa savannas.
            COAfr=1.
       IF(J.GE.32.AND.J.LE.50.AND.I.GE.24.AND.I.LE.29) COAfr=73./45.
cbnd
            ! Store biomass burning [molec CO] in TCO(I,J,1,6)
            TCO(I,J,1,6) = TCO(I,J,1,6) + 
     &                     ( BURNEMIS(IDBCO,IREF,JREF) * 
cbnd     &                       BOXVL(I,J,1) * DTSRCE )
     &                    BOXVL(I,J,1) * DTSRCE )*SCALEITBB*COAfr
cbnd

            ! Store biofuel burning [molec CO] in TCO(I,J,1,7)
            TCO(I,J,1,7) = TCO(I,J,1,7)+
!     &                     ( TWOODIJ(I,J) * BOXVL(I,J,1) * DTSRCE )
cbnd     &              ( BIOFUEL(2,IREF,JREF) * BOXVL(I,J,1) * DTSRCE )
     &    ( BIOFUEL(2,IREF,JREF) * BOXVL(I,J,1) * DTSRCE )*SCALEITBF
cbnd

        ! Store fossil fuels [molec CO]  in TCO(I,J,1,8)
        IF ( SCALEYEAR == 1985 ) THEN
            TCO(I,J,1,8) = TCO(I,J,1,8) +
     &                     ( EMISTCO(IREF,JREF) * DTSRCE * 
cbnd     &                DXYP(JREF)         * 1d4 )
     &                DXYP(JREF)         * 1d4 )*SCALEITFF
cbnd
        ! Non-baseline year
        ELSE

            TCO(I,J,1,8) = TCO(I,J,1,8) +
     &                     ( EMISTCO(IREF,JREF) * DTSRCE * 
cbnd     &                DXYP(JREF)  * 1d4 * FLIQCO2(IREF,JREF))
     &         DXYP(JREF)  * 1d4 * FLIQCO2(IREF,JREF))*SCALEITFF
cbnd
        ENDIF

        ENDDO
      ENDDO
!
!*****************************************************************************
!  Isoprene Emissions 
!  NOTE: Now use AD46 array instead of AIJ for diagnostics (bmy, 4/14/00)
!*****************************************************************************
!
      IJLOOP = 0

      DO J = 1, JJPAR
         JREF = J + J0
         DO I = 1, IIPAR
            
            ! 1-D loop index
            IJLOOP = IJLOOP + 1

            ! Surface temperature in K
            TMMP = XLTMMP(I,J,IJLOOP)

            ! EMXX = atoms C/box/time step from Isoprene
            EMXX = EMISOP( IJLOOP, SUNCOS, TMMP, XNUMOL_C ) 

            ! Sum isoprene emissions for use in SR CO_fromHCs.  
            ! Units are [atoms C/box/time step].
            SUMISOPCO(I,J) = SUMISOPCO(I,J) + EMXX

            ! Diagnostic in [atoms C/cm2/time step]
            ! Divide by DTSRCE in diag3.f to get [atoms C/cm2/s].
            IF ( ND46 > 0 ) THEN
               !EMX         = EMXX / ( DXYP(JREF) * 1D4 )
               AD46(I,J,1) = AD46(I,J,1) + EMX
            ENDIF
         ENDDO
      ENDDO
!
!*****************************************************************************
!  Monoterpene Emissions 
!  NOTE: Now use AD46 array instead of AIJ for diagnostics (bmy, 4/14/00)
!*****************************************************************************
!
      IJLOOP = 0
      DO J = 1, JJPAR
         JREF = J + J0
         DO I = 1, IIPAR

            ! 1-D loop index
            IJLOOP = IJLOOP + 1

            ! Surface temperature in K
            TMMP = XLTMMP(I,J,IJLOOP)

            ! EMMO = atoms C/box/time step from monoterpenes
            EMMO = EMMONOT( IJLOOP, TMMP, XNUMOL_C )

            ! Sum monoterpene emissions for use in SR CO_fromHCs.  
            ! Units are [atoms C/box/time step].
            SUMMONOCO(I,J) = SUMMONOCO(I,J) + EMMO

            ! ND46 Diagnostic in [atoms C/cm2/s] (bmy, 9/15/01)
            IF ( ND46 > 0 ) THEN
               EMX         = EMMO / ( DXYP(JREF) * 1D4 ) / DTSRCE
               AD46(I,J,2) = AD46(I,J,2) + EMX
            ENDIF
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE EMISSCO                                                 
 
