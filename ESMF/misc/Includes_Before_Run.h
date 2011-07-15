!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: Includes_Before_Run.h
!
! !DESCRIPTION: This include file contains the array assignments that need
!  to be made BEFORE the call to the Run method of the 
!  GEOSCHEMchem\_GridCompMod.F90 code.  These array assignments take data out 
!  of the import and internal states and saves them into locally-defined 
!  arrays (which pass the data down to the GEOS-Chem column chemistry code).
!\\
!\\
!  These assignments were placed into this separate include file to avoid 
!  bogging down the GEOSCHEMchem\_GridCompMod.F90 module.
!\\
!\\
! !REVISION HISTORY: 
!  16 Apr 2010 - R. Yantosca - Initial version
!  29 Apr 2010 - R. Yantosca - Add TO3 and FRCLAND to the MET object
!  06 May 2010 - R. Yantosca - Comment out dead species of CSPEC_1d
!  02 Jun 2010 - R. Yantosca - Now allocate Met%DQ*DTMST fields
!  02 Jun 2010 - R. Yantosca - Now allocate Met%TAUCL* fields
!  04 Jun 2010 - R. Yantosca - Bug fix: convert units of Met%AIRDENS to
!                              molec/m3 as is done in the GEOS-Chem
!  04 Jun 2010 - R. Yantosca - Now init Met%PBLH with the ZPBL state array
!  01 Jul 2010 - R. Yantosca - Add references to D_OH_MASS and D_AIR_MASS
!  02 Jul 2010 - R. Yantosca - Use OH_SCALE to avoid overflow for mean OH
!EOP
!------------------------------------------------------------------------------
!BOC
       !-------------------------------------------------------------------
       ! Get values from 2-D met field pointer arrays
       !-------------------------------------------------------------------
       Met%ALBD         =  ALBEDO    (I,J)               ! [m2]
       Met%CLDFRC       =  CLDFRC    (I,J)               ! [unitless]
       Met%FRCLND       =  FRLAND    (I,J)               ! [unitless]
       Met%HFLUX        =  SH        (I,J)               ! [W/m2]
       Met%LWI          =  LWI       (I,J)               ! [unitless]
       Met%PBLH         =  ZPBL      (I,J)               ! [m]
       Met%PRECCON      =  CN_PRCP   (I,J)               ! [kg/m2/s]
       Met%PRECTOT      =  AN_PRCP   (I,J) + &
                           CN_PRCP   (I,J) + &
                           LS_PRCP   (I,J)               ! [kg/m2/s]
       Met%RADSWG       =  RADSRF    (I,J)               ! [W/m2]
       Met%SST          =  TSEA      (I,J)               ! [K]
       Met%SUNCOS       =  COSZ      (I,J)               ! [unitless]  
       Met%TO3          =  TO3       (I,J)
       Met%TROPP        =  TROPP     (I,J)               ! [Pa]
       Met%TS           =  T2M       (I,J)               ! [K]
       Met%U10M         =  U10M      (I,J)               ! [m/s]
       Met%USTAR        =  USTAR     (I,J)               ! [m/s]
       Met%UVALBEDO     =  UVALBEDO  (I,J)               ! [unitless]
       Met%V10M         =  V10M      (I,J)               ! [m/s]
       Met%Z0           =  Z0        (I,J)               ! [m]       

       ! These are only needed for emissions, zero out
       Met%GWETTOP      = 0d0                            ! [unitless]
       Met%PARDR        = 0d0                            ! [W/m2]
       Met%PARDF        = 0d0                            ! [W/m2]

       ! Grid box surface areas [m2]
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%% NOTE: Assume a constant grid box size, therefore AREA_M2     %%%
       !%%% reduces to Re^2 * Delta-Lon * Delta-Lat (bmy, 4/12/10)       %%%
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Met%AREA_M2      = ( 2d0 * PI * Re**2 ) / DBLE( 72 )          &
                        * COS( latCtr(I,J) )                         &
                        * ( latEdg(I,11) - latEdg(I,10) )
                        !  NOTE the SIN term is constant, all boxes same size!

       ! Kludge, near poles we may get weird results due to the grid
       Met%AREA_M2      = ABS( Met%AREA_M2 )

       IF ( Ident%VERBOSE ) THEN
          WRITE( logLun, * ) '### AREA_M2 #1: ', Met%AREA_M2
       ENDIF

!       !### Hardwire AREA_M2 for testing the scientific part of the code
!       !### DEBUG only (bmy, 5/13/10)
!       Met%AREA_M2      = 1.839630e+11
!
!       IF ( Ident%VERBOSE ) THEN
!          WRITE( logLun, * ) '### AREA_M2 #2: ', Met%AREA_M2
!       ENDIF
       
       !-------------------------------------------------------------------
       ! Vector fields (copy instead of pointer assignment)  
       !
       ! NOTE: GEOS-Chem starts with level 1 = surface, which is the
       ! opposite of how GEOS-5 indexes fields.  Therefore we must flip
       ! the the vertical arrays.
       !
       ! ALSO NOTE: Met is dimensioned with REAL*8 pointer fields, where
       ! the pointers in the ESMF states are only REAL*4.  Therefore we
       ! will use a direct copy rather than a pointer reference.
       !-------------------------------------------------------------------
       
       ! Allocate the pointer array fields
       IF ( .not. ASSOCIATED ( Met%AD ) ) THEN 
          ALLOCATE( Met%AD( LM ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED ( Met%AIRDENS ) ) THEN
          ALLOCATE( Met%AIRDENS( LM ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%AIRVOL ) ) THEN
          ALLOCATE( Met%AIRVOL( LM ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%BXHEIGHT ) ) THEN
          ALLOCATE( Met%BXHEIGHT( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%CLDF ) ) THEN
          ALLOCATE( Met%CLDF    ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%CMFMC ) ) THEN
          ALLOCATE( Met%CMFMC   ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%DQIDTMST ) ) THEN
          ALLOCATE( Met%DQIDTMST( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%DQLDTMST ) ) THEN
          ALLOCATE( Met%DQLDTMST( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%DQVDTMST ) ) THEN
          ALLOCATE( Met%DQVDTMST( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%DTRAIN ) ) THEN
          ALLOCATE( Met%DTRAIN  ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%MOISTQ ) ) THEN
          ALLOCATE( Met%MOISTQ  ( LM   ), STAT=STATUS )  ! make from DQ*DTMST
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%OPTD ) ) THEN
          ALLOCATE( Met%OPTD    ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%PEDGE ) ) THEN
          ALLOCATE( Met%PEDGE   ( LM+1 ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF
       
       IF ( .not. ASSOCIATED( Met%PMID ) ) THEN
          ALLOCATE( Met%PMID    ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%RH ) ) THEN
          ALLOCATE( Met%RH      ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%SPHU ) ) THEN
          ALLOCATE( Met%SPHU    ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%T ) ) THEN
          ALLOCATE( Met%T       ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%TAUCLI ) ) THEN
          ALLOCATE( Met%TAUCLI  ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       IF ( .not. ASSOCIATED( Met%TAUCLW ) ) THEN
          ALLOCATE( Met%TAUCLW  ( LM   ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF


       ! Populate MET object
       Met%CLDF           = FCLDF  (I,J,LM:1:-1)     ! [unitless]
       Met%CMFMC          = CNV_FMC(I,J,LM:1:-1)     ! [kg/m2/s]
       Met%DTRAIN         = CNV_MFD(I,J,LM:1:-1)     ! [kg/m2/s]
       Met%MOISTQ         = MOISTQ (I,J,LM:1:-1)     ! [kg/kg/s]  !***
       Met%PEDGE(1:LM+1)  = PLE    (I,J,LM:0:-1)     ! [Pa]
       Met%PMID           = PL     (I,J,LM:1:-1)     ! [Pa]
       Met%T              = T      (I,J,LM:1:-1)     ! [K]
       Met%OPTD           = TAUCLI (I,J,LM:1:-1) + & ! [unitless]
                            TAUCLW (I,J,LM:1:-1) 
       Met%RH             = RH2    (I,J,LM:1:-1)     ! [unitless]       
       Met%SPHU           = Q      (I,J,LM:1:-1)     ! [kg/kg]
       Met%T              = T      (I,J,LM:1:-1)     ! [K]
       Met%TAUCLI         = TAUCLI (I,J,LM:1:-1)     ! [unitless]
       Met%TAUCLW         = TAUCLW (I,J,LM:1:-1)     ! [unitless]

       ! NOTE: These are not part of the import state yet
       ! But we can add these later (bmy, 6/2/10)
       Met%DQIDTMST       = 0d0                      ! [kg/kg/s]
       Met%DQLDTMST       = 0d0                      ! [kg/kg/s]
       Met%DQVDTMST       = 0d0                      ! [kg/kg/s]

       ! Make derived data fields
       DO L = 1, LM
          P1              = Met%PEDGE(L)             ! Bot edge prs [Pa]
          P2              = Met%PEDGE(L+1)           ! Top edge prs [Pa]

          Met%BXHEIGHT(L) = Rdg0                   & ! Grid box height [m]
                          * Met%T(L)               &           
                          * LOG( P1 / P2 )

          Met%AIRVOL  (L) = Met%BXHEIGHT(L)        & ! Grid box vol [m3]
                          * Met%AREA_M2
  
          Met%AD      (L) = ( P1 - P2 )            & ! Air mass [kg]
                          / g0                     &
                          * Met%AREA_M2

          Met%AIRDENS (L) = Met%AD(L)              & ! Air density [molec/m3]
                           * 1000.d0               & 
                           / ( Met%AIRVOL(L)*1d6 ) & 
                           * 6.02252d+23           &    
                           / 28.966d0              & 
                           * 1d6     
       ENDDO

       !-------------------------------------------------------------------
       ! Put strat Ox data for SCHEM from import state into local arrays
       ! %%% NOTE: This data is slated for replacement once we %%%
       ! %%%       implement the columnized LINOZ code!!!      %%%
       !-------------------------------------------------------------------

       ! Allocate the pointer array fields
       IF ( .not. ASSOCIATED( Schem%OH ) ) THEN
          ALLOCATE( Schem%OH( LM ), STAT=STATUS )
          VERIFY_(STATUS)       
       ENDIF

       IF ( .not. ASSOCIATED( Schem%PCO ) ) THEN
          ALLOCATE( Schem%PCO( LM ), STAT=STATUS )
          VERIFY_(STATUS)       
       ENDIF

       IF ( .not. ASSOCIATED( Schem%LCO ) ) THEN
          ALLOCATE( Schem%LCO( LM  ), STAT=STATUS )
          VERIFY_(STATUS)       
       ENDIF

       IF ( .not. ASSOCIATED( Schem%JVALUE ) ) THEN
          ALLOCATE( Schem%JVALUE( LM, 13 ), STAT=STATUS )
          VERIFY_(STATUS)
       ENDIF

       ! Strat OH, P(CO), L(CO)
       Schem%OH    ( 1:LM     )             = SOX_OH     ( I, J, LM:1:-1 )
       Schem%PCO   ( 1:LM     )             = SOX_PCO    ( I, J, LM:1:-1 )
       Schem%LCO   ( 1:LM     )             = SOX_LCO    ( I, J, LM:1:-1 )

       ! Strat J-values
       Schem%JVALUE( 1:LM, 1  )             = SOX_JV_NOx ( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 2  )             = SOX_JV_H2O2( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 3  )             = SOX_JV_ACET( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 4  )             = SOX_JV_MEK ( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 5  )             = SOX_JV_ALD2( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 6  )             = SOX_JV_RCHO( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 7  )             = SOX_JV_MVK ( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 8  )             = SOX_JV_MACR( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 9  )             = SOX_JV_R4N2( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 10 )             = SOX_JV_CH2O( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 11 )             = SOX_JV_N2O5( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 12 )             = SOX_JV_HNO4( I, J, LM:1:-1 )
       Schem%JVALUE( 1:LM, 13 )             = SOX_JV_MP  ( I, J, LM:1:-1 )

       !-------------------------------------------------------------------
       ! Put land type/LAI data from import state into local arrays
       !-------------------------------------------------------------------
       iReg_1d                              = IREG    (I,J)
       iLand_1d(1 )                         = ILAND_01(I,J)
       iLand_1d(2 )                         = ILAND_02(I,J)
       iLand_1d(3 )                         = ILAND_03(I,J)
       iLand_1d(4 )                         = ILAND_04(I,J)
       iLand_1d(5 )                         = ILAND_05(I,J)
       iLand_1d(6 )                         = ILAND_06(I,J)
       iLand_1d(7 )                         = ILAND_07(I,J)
       iLand_1d(8 )                         = ILAND_08(I,J)
       iLand_1d(9 )                         = ILAND_09(I,J)
       iLand_1d(10)                         = ILAND_10(I,J)
       iLand_1d(11)                         = ILAND_11(I,J)
       iLand_1d(12)                         = ILAND_12(I,J)
       iLand_1d(13)                         = ILAND_13(I,J)
       iLand_1d(14)                         = ILAND_14(I,J)
       iLand_1d(15)                         = ILAND_15(I,J)
       iUse_1d (1 )                         = IUSE_01 (I,J)
       iUse_1d (2 )                         = IUSE_02 (I,J)
       iUse_1d (3 )                         = IUSE_03 (I,J)
       iUse_1d (4 )                         = IUSE_04 (I,J)
       iUse_1d (5 )                         = IUSE_05 (I,J)
       iUse_1d (6 )                         = IUSE_06 (I,J)
       iUse_1d (7 )                         = IUSE_07 (I,J)
       iUse_1d (8 )                         = IUSE_08 (I,J)
       iUse_1d (9 )                         = IUSE_09 (I,J)
       iUse_1d (10)                         = IUSE_10 (I,J)
       iUse_1d (11)                         = IUSE_11 (I,J)
       iUse_1d (12)                         = IUSE_12 (I,J)
       iUse_1d (13)                         = IUSE_13 (I,J)
       iUse_1d (14)                         = IUSE_14 (I,J)
       iUse_1d (15)                         = IUSE_15 (I,J)
       lai_1d  (1 )                         = LAI_01  (I,J)
       lai_1d  (2 )                         = LAI_02  (I,J)
       lai_1d  (3 )                         = LAI_03  (I,J)
       lai_1d  (4 )                         = LAI_04  (I,J)
       lai_1d  (5 )                         = LAI_05  (I,J)
       lai_1d  (6 )                         = LAI_06  (I,J)
       lai_1d  (7 )                         = LAI_07  (I,J)
       lai_1d  (8 )                         = LAI_08  (I,J)
       lai_1d  (9 )                         = LAI_09  (I,J)
       lai_1d  (10)                         = LAI_10  (I,J)
       lai_1d  (11)                         = LAI_11  (I,J)
       lai_1d  (12)                         = LAI_12  (I,J)
       lai_1d  (13)                         = LAI_13  (I,J)
       lai_1d  (14)                         = LAI_14  (I,J)
       lai_1d  (15)                         = LAI_15  (I,J)  
       
       !-------------------------------------------------------------------
       ! Put emissions from import state into local arrays
       !-------------------------------------------------------------------
       emiss_1d( 1:LM, IdTracers%NOX  )     = EMISS_NOX ( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%Ox   )     = EMISS_O3  ( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%CO   )     = EMISS_CO  ( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%ALK4 )     = EMISS_ALK4( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%ISOP )     = EMISS_ISOP( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%HNO3 )     = EMISS_HNO3( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%ACET )     = EMISS_ACET( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%MEK  )     = EMISS_MEK ( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%ALD2 )     = EMISS_ALD2( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%PRPE )     = EMISS_PRPE( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%C3H8 )     = EMISS_C3H8( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%CH2O )     = EMISS_CH2O( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%C2H6 )     = EMISS_C2H6( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%DMS  )     = EMISS_DMS ( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%SO2  )     = EMISS_SO2 ( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%SO4  )     = EMISS_SO4 ( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%NH3  )     = EMISS_NH3 ( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%BCPI )     = EMISS_BCPI( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%OCPI )     = EMISS_OCPI( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%BCPO )     = EMISS_BCPO( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%OCPO )     = EMISS_OCPO( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%ALPH )     = EMISS_ALPH( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%LIMO )     = EMISS_LIMO( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%ALCO )     = EMISS_ALCO( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%DST1 )     = EMISS_DST1( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%DST2 )     = EMISS_DST2( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%DST3 )     = EMISS_DST3( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%DST4 )     = EMISS_DST4( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%SALA )     = EMISS_SALA( I, J, LM:1:-1 )
       emiss_1d( 1:LM, IdTracers%SALC )     = EMISS_SALC( I, J, LM:1:-1 )

       !-------------------------------------------------------------------
       ! Put SOA parameters from the internal state into local arrays
       !-------------------------------------------------------------------

       ! TERP and SESQ emissions
       terp_1d( 1:LM )                      = ORVC_TERP( I, J, LM:1:-1 ) 
       sesq_1d( 1:LM )                      = ORVC_SESQ( I, J, LM:1:-1 ) 
                                            
       ! GPROD parameters                   
       gProd_1d( 1:LM, 1, 1 )               = GSOAP11( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 2, 1 )               = GSOAP21( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 3, 1 )               = GSOAP31( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 1, 2 )               = GSOAP12( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 2, 2 )               = GSOAP22( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 3, 2 )               = GSOAP32( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 1, 3 )               = GSOAP13( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 2, 3 )               = GSOAP23( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 3, 3 )               = GSOAP33( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 1, 4 )               = GSOAP14( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 2, 4 )               = GSOAP24( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 3, 4 )               = GSOAP34( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 1, 5 )               = GSOAP15( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 2, 5 )               = GSOAP25( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 3, 5 )               = GSOAP35( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 1, 6 )               = GSOAP16( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 2, 6 )               = GSOAP26( I, J, LM:1:-1 )
       gProd_1d( 1:LM, 3, 6 )               = GSOAP36( I, J, LM:1:-1 )
                                            
       ! GPROD parameters                   
       aProd_1d( 1:LM, 1, 1 )               = ASOAP11( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 2, 1 )               = ASOAP21( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 3, 1 )               = ASOAP31( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 1, 2 )               = ASOAP12( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 2, 2 )               = ASOAP22( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 3, 2 )               = ASOAP32( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 1, 3 )               = ASOAP13( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 2, 3 )               = ASOAP23( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 3, 3 )               = ASOAP33( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 1, 4 )               = ASOAP14( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 2, 4 )               = ASOAP24( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 3, 4 )               = ASOAP34( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 1, 5 )               = ASOAP15( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 2, 5 )               = ASOAP25( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 3, 5 )               = ASOAP35( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 1, 6 )               = ASOAP16( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 2, 6 )               = ASOAP26( I, J, LM:1:-1 )
       aProd_1d( 1:LM, 3, 6 )               = ASOAP36( I, J, LM:1:-1 )

       !--------------------------------------------------------------------
       ! Put seasalt parameters from the internal state into local arrays
       !--------------------------------------------------------------------
       alkEmis_1d( 1:LM, 1 )                = ALK_EMIS_SALA( I, J, LM:1:-1 )
       alkEmis_1d( 1:LM, 2 )                = ALK_EMIS_SALC( I, J, LM:1:-1 )
       nDens_1d  ( 1:LM, 1 )                = N_DENS_SALA  ( I, J, LM:1:-1 )
       nDens_1d  ( 1:LM, 2 )                = N_DENS_SALC  ( I, J, LM:1:-1 )

       !--------------------------------------------------------------------
       ! Put wetdep parameters from the internal state into local arrays
       !--------------------------------------------------------------------
       h2o2s_1d( 1:LM )                     = H2O2s( I, J, LM:1:-1 )
       so2s_1d( 1:LM )                      = SO2s ( I, J, LM:1:-1 )

       !-------------------------------------------------------------------
       ! Put chemical species from internal state into CSPEC array
       !-------------------------------------------------------------------

       !%%% Active species in mechanism %%%
       cspec_1d( 1:LM, IdSpecies%A3O2     ) = A3O2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ACET     ) = ACET    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ALD2     ) = ALD2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ALK4     ) = ALK4    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ATO2     ) = ATO2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%B3O2     ) = B3O2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%C2H6     ) = C2H6    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%C3H8     ) = C3H8    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%CH2O     ) = CH2O    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%CO       ) = CO      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DMS      ) = DMS     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYCH2O  ) = DRYCH2O ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYH2O2  ) = DRYH2O2 ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYHNO3  ) = DRYHNO3 ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYN2O5  ) = DRYN2O5 ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYNO2   ) = DRYNO2  ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYO3    ) = DRYO3   ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYPAN   ) = DRYPAN  ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYPMN   ) = DRYPMN  ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYPPN   ) = DRYPPN  ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%DRYR4N2  ) = DRYR4N2 ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ETO2     ) = ETO2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ETP      ) = ETP     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%GCO3     ) = GCO3    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%GLYC     ) = GLYC    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%GP       ) = GP      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%GPAN     ) = GPAN    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%H2O2     ) = H2O2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%HAC      ) = HAC     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%HNO2     ) = HNO2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%HNO3     ) = HNO3    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%HNO4     ) = HNO4    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%HO2      ) = HO2     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%IALD     ) = IALD    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%IAO2     ) = IAO2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%IAP      ) = IAP     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%INO2     ) = INO2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%INPN     ) = INPN    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ISN1     ) = ISN1    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ISNP     ) = ISNP    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ISOP     ) = ISOP    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%KO2      ) = KO2     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%LISOPOH  ) = LISOPOH ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MACR     ) = MACR    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MAN2     ) = MAN2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MAO3     ) = MAO3    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MAOP     ) = MAOP    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MAP      ) = MAP     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MCO3     ) = MCO3    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MEK      ) = MEK     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MGLY     ) = MGLY    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MO2      ) = MO2     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MP       ) = MP      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MRO2     ) = MRO2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MRP      ) = MRP     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MSA      ) = SO4     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MVK      ) = MVK     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MVN2     ) = MVN2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%N2O5     ) = N2O5    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%NO       ) = NO      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%NO2      ) = NO2     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%NO3      ) = NO3     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%O3       ) = O3      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%OH       ) = OH      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%PAN      ) = PAN     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%PMN      ) = PMN     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%PO2      ) = PO2     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%PP       ) = PP      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%PRN1     ) = PRN1    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%PRPE     ) = PRPE    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%PRPN     ) = PRPN    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%R4N1     ) = R4N1    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%R4N2     ) = R4N2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%R4O2     ) = R4O2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%R4P      ) = R4P     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%RA3P     ) = RA3P    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%RB3P     ) = RB3P    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%RCHO     ) = RCHO    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%RCO3     ) = RCO3    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%RIO1     ) = RIO1    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%RIO2     ) = RIO2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%RIP      ) = RIP     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%RP       ) = RP      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%SO2      ) = SO2     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%SO4      ) = SO2     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%VRO2     ) = VRO2    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%VRP      ) = VRP     ( I, J, LM:1:-1 )

       !%%% Inactive species in mechanism %%%
       cspec_1d( 1:LM, IdSpecies%ACTA     ) = ACTA    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%CH4      ) = CH4     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%EMISSION ) = EMISSION( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%EOH      ) = EOH     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%GLCO3    ) = GLCO3   ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%GLP      ) = GLP     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%GLPAN    ) = GLPAN   ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%GLYX     ) = GLYX    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%H        ) = H       ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%H2       ) = H2      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%H2O      ) = H2O     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%HCOOH    ) = HCOOH   ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ISNO3    ) = ISNO3   ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%M        ) = M       ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MNO3     ) = MNO3    ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%MOH      ) = MOH     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%N2       ) = N2      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%NH2      ) = NH2     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%NH3      ) = NH3     ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%O        ) = O       ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%O2       ) = O2      ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%O2CH2OH  ) = O2CH2OH ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%RCOOH    ) = RCOOH   ( I, J, LM:1:-1 )
       cspec_1d( 1:LM, IdSpecies%ROH      ) = ROH     ( I, J, LM:1:-1 )

       !%%% Dead species -- comment out (bmy, 5/6/10) %%%
      !cspec_1d( 1:LM, IdSpecies%CO2      ) = CO2     ( I, J, LM:1:-1 )
      !cspec_1d( 1:LM, IdSpecies%DRYDEP   ) = DRYDEP  ( I, J, LM:1:-1 )
      !cspec_1d( 1:LM, IdSpecies%N2O      ) = N2O     ( I, J, LM:1:-1 )
      !cspec_1d( 1:LM, IdSpecies%O1D      ) = O1D     ( I, J, LM:1:-1 )
      !cspec_1d( 1:LM, IdSpecies%PPN      ) = PPN     ( I, J, LM:1:-1 )

       !-------------------------------------------------------------------
       ! Put advected tracers from internal state into TRACER array
       !-------------------------------------------------------------------
       tracer_1d( 1:LM, IdTracers%NOx  )    = TRC_NOx ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%Ox   )    = TRC_Ox  ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%PAN  )    = TRC_PAN ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%CO   )    = TRC_CO  ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%ALK4 )    = TRC_ALK4( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%ISOP )    = TRC_ISOP( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%HNO3 )    = TRC_HNO3( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%H2O2 )    = TRC_H2O2( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%ACET )    = TRC_ACET( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%MEK  )    = TRC_MEK ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%ALD2 )    = TRC_ALD2( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%RCHO )    = TRC_RCHO( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%MVK  )    = TRC_MVK ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%MACR )    = TRC_MACR( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%PMN  )    = TRC_PMN ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%PPN  )    = TRC_PPN ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%R4N2 )    = TRC_R4N2( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%PRPE )    = TRC_PRPE( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%C3H8 )    = TRC_C3H8( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%CH2O )    = TRC_CH2O( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%C2H6 )    = TRC_C2H6( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%N2O5 )    = TRC_N2O5( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%HNO4 )    = TRC_HNO4( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%MP   )    = TRC_MP  ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%DMS  )    = TRC_DMS ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SO2  )    = TRC_SO2 ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SO4  )    = TRC_SO4 ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SO4s )    = TRC_SO4s( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%MSA  )    = TRC_MSA ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%NH3  )    = TRC_NH3 ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%NH4  )    = TRC_NH4 ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%NIT  )    = TRC_NIT ( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%NITs )    = TRC_NITs( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%BCPI )    = TRC_BCPI( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%OCPI )    = TRC_OCPI( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%BCPO )    = TRC_BCPO( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%OCPO )    = TRC_OCPO( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%ALPH )    = TRC_ALPH( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%LIMO )    = TRC_LIMO( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%ALCO )    = TRC_ALCO( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SOG1 )    = TRC_SOG1( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SOG2 )    = TRC_SOG2( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SOG3 )    = TRC_SOG3( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SOG4 )    = TRC_SOG4( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SOA1 )    = TRC_SOA1( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SOA2 )    = TRC_SOA2( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SOA3 )    = TRC_SOA3( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SOA4 )    = TRC_SOA4( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%DST1 )    = TRC_DST1( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%DST2 )    = TRC_DST2( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%DST3 )    = TRC_DST3( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%DST4 )    = TRC_DST4( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SALA )    = TRC_SALA( I, J, LM:1:-1 )
       tracer_1d( 1:LM, IdTracers%SALC )    = TRC_SALC( I, J, LM:1:-1 )

       !-------------------------------------------------------------------
       ! Put mean OH parameters from the internal state into local arrays
       !
       ! NOTE: D_AIR_MASS and D_OH_MASS are REAL*4, so we need to make
       ! sure that we do not exceed 1e38, which is the max allowable value.
       ! Therefore, D_AIR_MASS and D_OH_MASS will store the values divided
       ! by OH_SCALE = 1d20 to prevent this overflow situation.
       !-------------------------------------------------------------------
       airMassDiag( 1:LM ) = DBLE( D_AIR_MASS( I, J, LM:1:-1 ) ) * OH_SCALE
       ohMassDiag ( 1:LM ) = DBLE( D_OH_MASS ( I, J, LM:1:-1 ) ) * OH_SCALE
