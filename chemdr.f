! $Id: chemdr.f,v 1.9 2004/07/15 18:17:44 bmy Exp $
      SUBROUTINE CHEMDR
!
!******************************************************************************
!  Subroutine CHEMDR is the driver subroutine for full chemistry w/ SMVGEAR.
!  Adapted from original code by lwh, jyl, gmg, djj. (bmy, 11/15/01, 7/15/04)
!
!  Important input variables from "dao_mod.f" and "uvalbedo_mod.f":
!  ============================================================================
!  ALBD        : DAO visible albedo                         [unitless]
!  AVGW        : Mixing ratio of water vapor                [v/v] 
!  BXHEIGHT    : Grid box heights                           [m]
!  OPTD        : DAO grid-box optical depths (for FAST-J)   [unitless]
!  SUNCOS      : Cosine of solar zenith angle               [unitless]
!  SUNCOSB     : Cosine of solar zenith angle 1 hr from now [unitless]
!  UVALBEDO    : TOMS UV albedo 340-380 nm (for FAST-J)     [unitless]
!
!  Important input variables from "comode.h" or "comode_mod.f":
!  ============================================================================
!  NPTS        : Number of points (grid-boxes) to calculate
!  REMIS       : Emission rates                             [molec/cm3/s-1]
!  RAERSOL     : Frequency of gas-aerosol collisions        [s-1]
!  PRESS       : Pressure                                   [Pa]
!  TMPK        : Temperature                                [K]
!  ABSHUM      : Absolute humidity                          [molec/cm3]
!  ALT         : Altitude                                   [cm]
!  SURFALT     : Surface altitude                           [m]
!  TOTO3       : Total ozone column                         [molec/cm3]
!  CLOUDS      : Albedos at 2-km intervals from 0 to 20-km  
!  IDXAIR      : Index for standard temperature profile
!  IDXO3       : Index for standard ozone profile
!  CSPEC       : Initial species concentrations             [molec/cm3]
!
!  Important output variables in "comode.h" etc.
!  ============================================================================
!  NAMESPEC    : Character array of species names
!  NNSPEC      : # of ACTIVE + INACTIVE (not DEAD) species
!  CSPEC       : Final species concentrations               [molec/cm3]
!
!  Other Important Variables
!  ============================================================================
!  MAXPTS      : Maximum number of points or grid-boxes (in "comsol.h")
!                (NPTS must be <= MAXPTS, for SLOW-J only)
!  MAXDEP      : Maximum number of deposition species (note # of
!                depositing species listed in tracer.dat must be <= MAXDEP)
!  IGAS        : Maximum number of gases, ACTIVE + INACTIVE
!  IO93        : I/O unit for output for "ctm.chem" file
!
!  Input files for SMVGEAR II:
!  ============================================================================
!   mglob.dat  : control switches                       (read in "reader.f")
!  tracer.dat  : list of tracers, emitting species      (read in "reader.f")
!                and depositing species
! globchem.dat : species list, reaction list,           (read in "chemset.f")
!                photolysis reaction list
!
!  Input files for FAST-J photolysis:
!  ============================================================================
!     ratj.d   : Lists photo species, branching ratios  (read in "rd_js.f")
! jv_atms.dat  : Climatology of T and O3                (read in "rd_prof.f")
! jv_spec.dat  : Cross-sections for each species        (read in "RD_TJPL.f")
!
!  Input files for SLOW-J photolysis:
!  ============================================================================
!  jvalue.dat  : Solar flux data, standard T and O3     (read in "jvaluein.f")
!                profiles, aerosol optical depths 
!    8col.dat  : SLOW-J cross-section data              (read in "jvaluein.f")
!  chemga.dat  : Aerosol data
!    o3du.dat  : O3 in Dobson units, cloud data         (read in "jvaluein.f")
!
!  NOTES:
!  (1 ) Cleaned up a lot of stuff.  SUNCOS, OPTD, ALBD, and AVGW are now 
!        referenced from dao_mod.f.  IREF and JREF are obsolete.  Also 
!        updated comments. (bmy, 9/27/01)
!  (2 ) Do not declare LPRT or set LPRT = .FALSE. in "chemdr.f".  LPRT is 
!        included via "CMN" and is defined in "main.f". (bmy, 10/9/01)
!  (3 ) Removed obsolete data from 9/01 (bmy, 10/23/01)
!  (4 ) ERADIUS(JLOOP) is now ERADIUS(JLOOP,1) and TAREA(JLOOP) is now
!        TAREA(JLOOP,1) for sulfate aerosol.  Updated comments. (bmy, 11/15/01)
!  (5 ) Renamed routine PAFTOP to DEBUG_MSG.  Also deleted obsolete code
!        from 11/01.  Enhanced debug output via DEBUG_MSG.  Also reference
!        the UVALBEDO array directly from "uvalbedo_mod.f".  Remove UVALBEDO
!        from the argument list.  Removed obsolete comments. (bmy, 1/15/02)
!  (6 ) Now pass LPAUSE to "initgas.f" via the arg list (bmy, 2/14/02)
!  (7 ) Now call "rdaer.f" instead of RDAEROSOL to read the aerosol and dust 
!        fields from disk.  Also, ignore hygroscopic growth for dust.  Now
!        pass SAVEHO2 and FRACNO2 arrays in the arg list to "ohsave.f"; these 
!        return HO2 conc.'s and NO2 fraction.  Delete NTRACE from call
!        to "ohsave.f", it's obsolete.  Delete reference to DARSFCA from
!        "comode_mod.f", it's obsolete. (rvm, bmy, 2/27/02)
!  (8 ) Removed obsolete code from 2/02. (bmy, 4/15/02)
!  (9 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (10) Now reference IU_CTMCHEM from "file_mod.f".  Assign the value of
!        IU_CTMCHEM (which =93) to IO93 for SMVGEAR routines.  Also open 
!        "ctm.chem" file on the first call as file unit # IO93.  Add
!        references to "CMN_DIAG" and "planeflight_mod.f".  Call routine
!        SETUP_PLANEFLIGHT to initialize the plane track diagnostic
!        after reading the "chem.dat" file.  (bmy, 7/2/02)
!  (11) Now reference AD, T and BXHEIGHT from "dao_mod.f".  Also removed 
!        obsolete commented out code in various sections below.  Now also
!        references "tracerid_mod.f".  Also remove reference to BIOTRCE, since
!        this is now obsolete.  Now make FIRSTCHEM a local SAVED variable
!        instead of an argument.  Now calls MAKE_AVGW, which was formerly
!        called in "main.f". (bmy, 11/15/02)
!  (12) Now reference "AIRVOL" from "dao_mod.f".  Now declare local array
!        SO4_NH4_NIT, which will contain lumped SO4, NH3, NIT aerosol.  Now
!        pass SO4_NH4_NIT to "rdaer.f" via the arg list if sulfate chemistry
!        is turned on.  Now also references CMN_SETUP. (rjp, bmy, 3/23/03)
!  (13) Removed ITAU from the arg list.  Removed reference to IHOUR.  Now use
!        functions GET_MONTH, GET_YEAR from "time_mod.f" (bmy, 3/27/03)
!  (14) Remove KYEAR and TWO_PI, these are now obsolete for SMVGEAR II.  Now 
!        open unit #93 and call READER in the same FIRSTCHEM if-block.  Now
!        Replace call to CHEMSET with call to READCHEM.  JPARSE is now called 
!        from w/in READCHEM.  Replace call to INITGAS w/ call to GASCONC.
!        Removed reference to "file_mod.f".  Remove call to SETPL, we now must
!        call this in "readchem.f" before the call to JSPARSE. 
!        (bdf, ljm, bmy, 5/8/03)
!  (15) Now reference routine GET_GLOBAL_CH4 from "global_ch4_mod.f".  Also
!        added CH4_YEAR as a SAVEd variable. (bnd, bmy, 7/1/03)
!  (16) Remove references to MONTHP, IMIN, ISEC; they are obsolete and not 
!        defined anywhere. (bmy, 7/16/03)
!  (17) Now reference SUNCOSB from "dao_mod.f".  Now pass SUNCOSB to "chem.f". 
!        Also remove LSAMERAD from call to CHEM, since it's obsolete. 
!        (gcc, bmy, 7/30/03)
!  (18) Added BCPO, BCPI, OCPO, OCPI, and SOILDUST arrays.  Now pass SOILDUST
!       to RDUST_ONLINE (in "dust_mod.f").  Now pass PIEC, POEC, PIOC, POOC to
!       "rdaer.f".  Now references "dust_mod.f". (rjp, tdf, bmy, 4/1/04)
!  (19) Added SALA and SALC arrays for passing seasalt to rdaer.f.  Now
!        rearranged the DO loop for computational efficiency. (bmy, 4/20/04)
!  (20) Added OCF parameter to account for the other chemical components that
!        are attached to OC.  Also now handle hydrophilic OC differently for
!        online & offline SOA. (rjp, bmy, 7/15/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,      ONLY : ABSHUM,   CSPEC,     ERADIUS, TAREA
      USE DAO_MOD,         ONLY : AD,       AIRVOL,    ALBD,    AVGW,    
     &                            BXHEIGHT, MAKE_AVGW, OPTD,    SUNCOS,  
     &                            SUNCOSB,  T
      USE DUST_MOD,        ONLY : RDUST_ONLINE, RDUST_OFFLINE
      USE ERROR_MOD,       ONLY : DEBUG_MSG
      USE GLOBAL_CH4_MOD,  ONLY : GET_GLOBAL_CH4
      USE PLANEFLIGHT_MOD, ONLY : SETUP_PLANEFLIGHT
      USE TIME_MOD,        ONLY : GET_MONTH, GET_YEAR
      USE TRACERID_MOD 
      USE UVALBEDO_MOD,    ONLY : UVALBEDO

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT, P, T, NSRCX, MONTH, JYEAR, etc.
#     include "CMN_O3"    ! EMISRRN, EMISRR
#     include "CMN_NOX"   ! SLBASE
#     include "comode.h"  ! SMVGEAR variables
#     include "CMN_DEP"   ! FRCLND
#     include "CMN_DIAG"  ! ND40
#     include "CMN_SETUP" ! LSULF, LCARB, LSOA, LDUST, LSSALT

      ! Local variables
      LOGICAL, SAVE       :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE       :: CH4_YEAR  = -1
      INTEGER             :: I, J, JLOOP, L, NPTS, N, MONTH, YEAR
      INTEGER             :: IDXAIR(JJPAR)
      INTEGER             :: IDXO3(JJPAR)
      REAL*8              :: XWETRAT, HUMEFF, ROVMG
      REAL*8              :: ALT(MAXIJ,IVERT) 
      REAL*8              :: SURFALT(MAXIJ)
      REAL*8              :: TOTO3(JJPAR)
      REAL*8              :: CLOUDS(MAXIJ,11)
      REAL*8              :: SO4_NH4_NIT(IIPAR,JJPAR,LLPAR)
      REAL*8              :: BCPI(IIPAR,JJPAR,LLPAR)
      REAL*8              :: BCPO(IIPAR,JJPAR,LLPAR)
      REAL*8              :: OCPI(IIPAR,JJPAR,LLPAR)
      REAL*8              :: OCPO(IIPAR,JJPAR,LLPAR)
      REAL*8              :: SALA(IIPAR,JJPAR,LLPAR)
      REAL*8              :: SALC(IIPAR,JJPAR,LLPAR)
      REAL*8              :: SOILDUST(IIPAR,JJPAR,LLPAR,NDUST)

      ! We carry carbon mass only in OC and here need to multiply by
      ! 1.4 to account for the mass of the other chemical components
      ! (rjp, bmy, 7/15/04)
      REAL*8, PARAMETER   :: OCF = 1.4d0

      !=================================================================
      ! CHEMDR begins here!
      !=================================================================

      ! Set some size variables
      NLAT   = JJPAR
      NLONG  = IIPAR
      NVERT  = IVERT 
      NPVERT = NVERT
      NPVERT = NVERT + IPLUME

      ! Get month and year
      MONTH = GET_MONTH()
      YEAR  = GET_YEAR()

      !=================================================================
      ! Compute AVGW, the mixing ratio of water vapor
      !=================================================================
      CALL MAKE_AVGW

      !=================================================================
      ! Open "smv2.log" output file and read chem mechanism switches
      !=================================================================
      IF ( FIRSTCHEM ) THEN
         
         ! Read from data files m.dat and tracer.dat
         CALL READER( FIRSTCHEM )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after READER' )

         ! Set NCS for urban chemistry only (since that is where we
         ! have defined the GEOS-CHEM mechanism) (bdf, bmy, 4/21/03)
         NCS = NCSURBAN
      ENDIF

      !=================================================================
      ! For SLOW-J photolysis only: 
      ! Call JVALUEIN which reads in parameters for J-value calculation.  
      ! JVALUEIN reads O3DU from the file "o3du.dat"
      !=================================================================
#if   defined( LSLOWJ )
      IF ( FIRSTCHEM ) CALL JVALUEIN            
      IF ( LPRT      ) CALL DEBUG_MSG( '### CHEMDR: after JVALUEIN' )
#endif

      !=================================================================      
      ! Call GETALT, which GETS THE ALTITUDE   
      !=================================================================
      CALL GETALT( ALT, SURFALT, ROVMG )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after GETALT' ) 

      !=================================================================      
      ! Call RURALBOX, which defines tropospheric boxes to be sent to
      ! the SMVGEAR solver, as well as setting up some SMVGEAR arrays.
      !=================================================================      
      NLOOP       = NLAT  * NLONG
      NTLOOP      = NLOOP * NVERT

      CALL RURALBOX( AD,     T,      AVGW,   ALT,   ALBD,  
     &               SUNCOS, CLOUDS, LEMBED, IEBD1, IEBD2, 
     &               JEBD1,  JEBD2,  LPAUSE )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after RURALBOX' ) 

      !=================================================================      
      ! For SLOW-J Photolysis only: 
      ! Call GETTOTO3, which gets the total ozone column 
      !=================================================================      
#if   defined( LSLOWJ )
      CALL GETTOTO3( TOTO3, IDXAIR, IDXO3 )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after GETTOTO3' )
#endif

      !=================================================================
      ! NTTLOOP is the total number of tropospheric grid boxes.
      ! SMVGEAR will do chemistry in all NTTLOOP grid boxes.
      !=================================================================
      NTTLOOP = NTLOOP

      !=================================================================
      ! Call SETMODEL which defines number of grid-blocks in calculation,
      ! and copies meteorological parameters into local variables 
      !=================================================================
      CALL SETMODEL

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after SETMODEL' )

      ! Do the following only on the first call ...
      IF ( FIRSTCHEM ) THEN

         ! Zero SO4_NH4_NIT on the first call only
         SO4_NH4_NIT = 0d0

         ! Zero EC/OC on the first call only
         BCPI        = 0d0
         BCPO        = 0d0
         OCPI        = 0d0
         OCPO        = 0d0
         SALA        = 0d0
         SALC        = 0d0

         ! Zero dust on the first call only
         SOILDUST    = 0d0

         !==============================================================
         ! Set up stuff for chemical prod/loss diagnostic (bey-02/99)
         !==============================================================
         IF ( NFAM > 0 ) THEN
            LFAMILY = .TRUE.
         ELSE
            LFAMILY = .FALSE.
         ENDIF

         NFAMILIES = NFAM
         write(6,*) 'NFAMILIES = NFAM = :', NFAMILIES, NFAM

         !==============================================================
         ! Call READCHEM -- the setup routine for gas-phase chemistry
         !==============================================================
         CALL READCHEM

         ! Set NCS=NCSURBAN here since we have defined our tropospheric
         ! chemistry mechanism in the urban slot of SMVGEAR II (bmy, 4/21/03)
         NCS = NCSURBAN

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after CHEMSET' )

         ! Set up plane flight diagnostic here -- test 
         IF ( ND40 > 0 ) CALL SETUP_PLANEFLIGHT

         !==============================================================
         ! If CH4 is a SMVGEAR II species, then call GET_GLOBAL_CH4
         ! to return the globally-varying CH4 conc. as a function of
         ! year and latitude bin.  (ICH4 is defined in READCHEM.)
         ! (bnd, bmy, 7/1/03)
         !==============================================================
         IF ( ICH4 > 0 .and. ( CH4_YEAR /= GET_YEAR() ) ) THEN

            ! Get CH4 [ppbv] in 4 latitude bins for each year
            CALL GET_GLOBAL_CH4( GET_YEAR(), .TRUE., C3090S, 
     &                           C0030S,     C0030N, C3090N )
            
            ! Save year for CH4 emissions
            CH4_YEAR = GET_YEAR()
         ENDIF

         !==============================================================
         ! If using Fast-J initialize it - must be after chemset 
         ! Pass ALL variables via argument list, except PTOP, which
         ! is declared as a parameter in "CMN_SIZE" (bmy, 2/10/00)
         !==============================================================
#if   defined( LFASTJ ) 
         CALL INPHOT( LLTROP, NPHOT ) 

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after INPHOT' )        
#endif

         !==============================================================
         ! Call SETTRACE which flags certain chemical species
         !==============================================================
         CALL SETTRACE( NTRACE )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after SETTRACE' )

         !==============================================================
         ! Call SETEMDEP which flags emission and dry deposition 
         ! reactions (read from chem.dat)
         !==============================================================
         CALL SETEMDEP( NTRACE )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after SETEMDEP' )

      ENDIF ! IF (FIRSTCHEM)

      ! Skip this section if all these are turned off
      IF ( LSULF .or. LCARB .or. LDUST .or. LSSALT ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !===========================================================
            ! Dump hydrophilic aerosols into one array that will be 
            ! passed to RDAER and then used for heterogeneous chemistry 
            ! as well as photolysis rate calculations interatively. 
            !
            ! If LSULF=F, then we read these aerosol data from Mian's 
            ! simulation.  If LSULF=T then we use the online tracers.
            !
            ! Now assume that all sulfate, ammonium, and nitrate are 
            ! hydrophilic but sooner or later we can pass only 
            ! hydrophilic aerosols from the thermodynamic calculations 
            ! for this purpose.  This dumping should be done before 
            ! calling INITGAS, which converts the unit of STT array 
            ! from kg/box to molec/cm3.
            !
            ! Units of SO4_NH4_NIT are [kg/m3].  (rjp, bmy, 3/23/03)
            !===========================================================
            IF ( LSULF ) THEN
               SO4_NH4_NIT(I,J,L) = ( STT(I,J,L,IDTSO4) + 
     &                                STT(I,J,L,IDTNH4) +
     &                                STT(I,J,L,IDTNIT) ) / 
     &                              AIRVOL(I,J,L)
            ENDIF

            !===========================================================
            ! Compute hydrophilic and hydrophobic BC and OC in [kg/m3]
            ! for passing to "rdaer.f" 
            !===========================================================
            !----------------------------------------------------------------
            ! Prior to 7/15/04:
            !IF ( LCARB ) THEN
            !
            !   ! Hydrophilic BC [kg/m3]
            !   BCPI(I,J,L) = STT(I,J,L,IDTBCPI) / AIRVOL(I,J,L)
            !
            !   ! Hydrophilic OC [kg/m3]
            !   OCPI(I,J,L) = STT(I,J,L,IDTOCPI) / AIRVOL(I,J,L)
            !
            !   ! Hydrophobic BC [kg/m3]
            !   BCPO(I,J,L) = STT(I,J,L,IDTBCPO) / AIRVOL(I,J,L)
            !
            !   ! Hydrophobic OC [kg/m3]                
            !   OCPO(I,J,L) = STT(I,J,L,IDTOCPO) / AIRVOL(I,J,L)
            !
            !   ! Now avoid division by zero (bmy, 4/20/04)
            !   BCPI(I,J,L) = MAX( BCPI(I,J,L), 1d-35 )
            !   OCPI(I,J,L) = MAX( OCPI(I,J,L), 1d-35 )
            !   BCPO(I,J,L) = MAX( BCPO(I,J,L), 1d-35 )
            !   OCPO(I,J,L) = MAX( OCPO(I,J,L), 1d-35 )
            !ENDIF
            !----------------------------------------------------------------
            IF ( LCARB ) THEN

               ! Hydrophilic BC [kg/m3]
               BCPI(I,J,L) = STT(I,J,L,IDTBCPI) / AIRVOL(I,J,L)

               ! Hydrophobic BC [kg/m3]
               BCPO(I,J,L) = STT(I,J,L,IDTBCPO) / AIRVOL(I,J,L)

               ! Hydrophobic OC [kg/m3] 
               OCPO(I,J,L) = STT(I,J,L,IDTOCPO) * OCF / AIRVOL(I,J,L)

               IF ( LSOA ) THEN

                  ! Hydrophilic primary OC plus SOA [kg/m3A.  We need
                  ! to multiply by OCF to account for the mass of other 
                  ! components which are attached to the OC aerosol.
                  ! (rjp, bmy, 7/15/04)
                  OCPI(I,J,L) = ( STT(I,J,L,IDTOCPI) * OCF 
     &                        +   STT(I,J,L,IDTSOA1)
     &                        +   STT(I,J,L,IDTSOA2)
     &                        +   STT(I,J,L,IDTSOA3)  ) / AIRVOL(I,J,L)
 
               ELSE

                  ! Hydrophilic primary and SOA OC [kg/m3].   We need
                  ! to multiply by OCF to account for the mass of other 
                  ! components which are attached to the OC aerosol.
                  ! (rjp, bmy, 7/15/04)
                  OCPI(I,J,L) = STT(I,J,L,IDTOCPI) * OCF / AIRVOL(I,J,L)
                  
               ENDIF

               ! Now avoid division by zero (bmy, 4/20/04)
               BCPI(I,J,L) = MAX( BCPI(I,J,L), 1d-35 )
               OCPI(I,J,L) = MAX( OCPI(I,J,L), 1d-35 )
               BCPO(I,J,L) = MAX( BCPO(I,J,L), 1d-35 )
               OCPO(I,J,L) = MAX( OCPO(I,J,L), 1d-35 )
            ENDIF

            !===========================================================
            ! Full chemistry with dust aerosol turned on.
            !
            ! Note, we can do better than this! Currently we carry 4 
            ! dust tracers...but het. chem and fast-j use 7 dust bins 
            ! hardwired from Ginoux.
            !
            ! Now, I apportion the first dust tracer into four smallest 
            ! dust bins equally in mass for het. chem and fast-j. 
            ! 
            ! Maybe we need to think about chaning our fast-j and het. 
            ! chem to use just four dust bins or more flexible 
            ! calculations depending on the number of dust bins. 
            ! (rjp, 03/27/04)
            !===========================================================
            IF ( LDUST ) THEN

               ! Lump 1st dust tracer for het chem
               DO N = 1, 4
                  SOILDUST(I,J,L,N) = 
     &                 0.25d0 * STT(I,J,L,IDTDST1) / AIRVOL(I,J,L)
               ENDDO

               ! Other hetchem bins
               SOILDUST(I,J,L,5) = STT(I,J,L,IDTDST2) / AIRVOL(I,J,L)
               SOILDUST(I,J,L,6) = STT(I,J,L,IDTDST3) / AIRVOL(I,J,L)
               SOILDUST(I,J,L,7) = STT(I,J,L,IDTDST4) / AIRVOL(I,J,L)
            ENDIF

            !===========================================================
            ! Sea salt tracers
            !===========================================================
            IF ( LSSALT ) THEN

               ! Accumulation mode seasalt aerosol [kg/m3]
               SALA(I,J,L) = STT(I,J,L,IDTSALA) / AIRVOL(I,J,L)
            
               ! Coarse mode seasalt aerosol [kg/m3]
               SALC(I,J,L) = STT(I,J,L,IDTSALC) / AIRVOL(I,J,L)

               ! Avoid division by zero
               SALA(I,J,L) = MAX( SALA(I,J,L), 1d-35 )
               SALC(I,J,L) = MAX( SALC(I,J,L), 1d-35 )
            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      !=================================================================
      ! Call GASCONC which initializes gas concentrations and sets 
      ! miscellaneous parameters.  GASCONC also calls PARTITION, which
      ! splits up family tracers like NOx and Ox into individual
      ! chemical species for SMVGEAR.
      !=================================================================
      CALL GASCONC( FIRSTCHEM, STT, NTRACE, XNUMOL, LPAUSE, FRCLND )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after GASCONC' )

      !=================================================================  
      ! Call SETEMIS which sets emission rates REMIS 
      ! SETEMIS now references BURNEMIS, BIOTRCE, NBIOTRCE,
      ! and BIOFUEL from F90 module "biomass_mod.f" (bmy, 9/12/00)
      !=================================================================
      CALL SETEMIS( EMISRR, EMISRRN )
      
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after SETEMIS' )

      !=================================================================
      ! Call RDAER -- computes aerosol optical depths
      !=================================================================
      !-----------------------------------------------------------------
      ! Prior to 4/20/04:
      !CALL RDAER( MONTH, YEAR, SO4_NH4_NIT, BCPI, BCPO, OCPI, OCPO )
      !-----------------------------------------------------------------
      CALL RDAER( MONTH, YEAR, SO4_NH4_NIT, 
     &            BCPI,  BCPO, OCPI, 
     &            OCPO,  SALA, SALC )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after RDAER' )

      !=================================================================
      ! If LDUST is turned on, then we have online dust aerosol in
      ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order to
      ! compute aerosol optical depth for FAST-J, etc.
      !
      ! If LDUST is turned off, then we do not have online dust aerosol
      ! in GEOS-CHEM...so read monthly-mean dust files from disk.
      ! (rjp, tdf, bmy, 4/1/04)
      !=================================================================
      IF ( LDUST ) THEN
         CALL RDUST_ONLINE( SOILDUST )
      ELSE
         CALL RDUST_OFFLINE( MONTH, YEAR )
      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after RDUST' )

      NPTS = NTTLOOP

      !=================================================================
      ! Call CHEM which PERFORMS GAS-PHASE CHEMISTRY. 
      !=================================================================
      CALL CHEM( FIRSTCHEM, NPTS,  SUNCOS, SUNCOSB, CLOUDS, ALT,
     &           SURFALT,   TOTO3, IDXAIR, IDXO3,   OPTD,   UVALBEDO )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after CHEM' )

      !=================================================================
      ! Call LUMP which LUMPS THE SPECIES TOGETHER.
      !=================================================================
      CALL LUMP( NTRACE, XNUMOL, STT )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after LUMP' )

      !=================================================================
      ! Call OHSAVE which saves info on OZONE, OH, AND NO concentrations 
      !=================================================================
      IF ( IDTNOX /= 0 .AND. IDTOX /= 0 ) THEN
         CALL OHSAVE( XNUMOL, STT,     FRACO3, FRACNO,  FRACNO2, 
     &                SAVEOH, SAVEHO2, SAVENO, SAVENO2, SAVENO3 )

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after OHSAVE' )
      ENDIF

      !=================================================================
      ! Call PLSAVE:  Chemical prod-loss diagnostic (bey-02/99)
      !=================================================================
      IF ( NFAM > 0 ) THEN
         CALL PLSAVE 

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMDR: after PLSAVE' )
      ENDIF
      
      !=================================================================
      ! Set FIRSTCHEM = .FALSE. -- we have gone thru one chem step
      !=================================================================
      FIRSTCHEM = .FALSE.

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### Now exiting CHEMDR!' )

      ! Return to calling program
      END SUBROUTINE CHEMDR




