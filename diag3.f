! $Id: diag3.f,v 1.5 2003/10/30 16:17:17 bmy Exp $
      SUBROUTINE DIAG3                                                      
! 
!******************************************************************************
!  Subroutine DIAG3 prints out I-J (Long-Lat) diagnostics to the BINARY
!  format punch file (bmy, bey, mgs, rvm, 5/27/99, 8/20/03)
!
!  The preferred file format is binary punch file format v. 2.0.  This
!  file format is very GAMAP-friendly.  GAMAP also supports the ASCII
!  style punch file, for backwards compatibility.
!
!  If LBPNCH = 1 call BPCH2 to print to the BINARY punch file 'ctm.bpch'
!
!  Also, in the "diag.dat" input file, one can specify which tracers will
!  be printed out for a given diagnostic.
!
!  At some point we will improve the diagnostic output, as much of this
!  is historical baggage. (bmy, 1/10/03)
!
!  NOTES: 
!  (1 ) Move the binary punch file code ahead of the ASCII punch file code
!        so that the multi-level diagnostics will not rewrite the title
!        strings for the binary punch file (bmy, 4/12/99)
!  (2 ) Rename "TITLE" to "CATEGORY" for consistency w/ binary punch
!        file format naming, and make it CHAR*40 (bmy, 5/27/99)
!  (3 ) Fixed bug for ND29.  Now test N against PD29. (bmy, 10/1/99)
!  (4 ) Cosmetic changes and some bug fixes (bmy, 10/15/99)
!  (5 ) Added ND14, ND15 diagnostics for mass flux, and updated
!        some comments.  Also freed up ND41 and ND42 diagnostics
!        since these were not being used.  (bey, bmy, 11/15/99)
!  (6 ) Fixed more bugs for ND32 binary punch file output (bmy, 11/18/99)
!  (7 ) Restrict adding TRCOFFSET to N for certain diagnostics
!        and/or simulations that don't require it.  Also updated 
!        some comments. (bmy, 11/23/99)
!  (8 ) Certain arrays for mass flux diagnostics and some NDxx diagnostics
!        are now declared allocatable in "diag_mod.f". (bmy, 11/30/99)
!  (9 ) Add ND55 diagnostics for tropopause (bmy, 11/30/99)
!  (10) Now make ND41 use allocatable array AD41 (bmy, 12/6/99)
!  (11) Bug fix for ND27 diagnostic (bmy, 2/7/00)
!  (12) ND31, ND33, ND35, ND67, and ND69 now use dynamically 
!        allocatable arrays declared in "diag_mod.f". (bmy, 2/17/00)
!  (13) Add air density to ND68 diagnostic (bmy, 3/10/00)
!  (14) As of 3/16/00, all diagnostics now use dynamically allocatable
!        arrays.  AIJ is no longer used (bmy, 3/16/00)
!  (15) Bug fix for ND47 diagnostic (bmy, 3/31/00)
!  (16) Fixes to ND43, ND46 diagnostics for CO run with parameterized OH 
!        (bmy, bnd, 4/18/00)
!  (17) For ND36 diagnostic for CH3I simulation, cycle when N > NEMANTHRO
!        instead of when N > NTRACE (mgs, hsu, bmy, 5/17/00)
!  (18) For DMS/SO2/SO4/MSA simulation, add TRCOFFSET to ND44 diagnostic
!        for dry deposition velocities & fluxes.  Also update the ND44 
!        error checking to avoid array-out-of-bounds errors. (bmy, 5/31/00)
!  (19) Now reference "bpch2_mod.f" which contains routines BPCH2 and 
!        GET_MODELNAME, which are used to write data to binary punch file 
!        format. (bmy, 6/22/00)
!  (20) Now reference LWI from "dao_mod.f" instead of from common block
!        header file "CMN_LWI". (bmy, 6/26/00)
!  (21) Added CLEVEL declaration for GEOS-3 data.  Also eliminated
!        obsolete code from 2/00 and 5/00. (bmy, 8/29/00)
!  (22) Added reference to BIOTRCE in "biomass_mod.f" and removed reference
!        to "bioburn.h".  Also added biofuel burning to ND32. (bmy, 9/12/00)
!  (23) Add TROPP and SLP to ND67 diagnostic.  Also changed the scale factor
!        for ND15 from SCALEDYN to SCALECONV (as it should be). 
!        (bmy, 10/18/00)
!  (24) Write out the proper number of levels to the punch file for ND14,
!        ND24, ND25, ND26, ND31, ND38, ND39. Also remove obsolete code from 
!        9/00. (bmy, 12/7/00)
!  (25) Don't scale the monthly mean AOD's by SCALECHEM (rvm, bmy, 12/8/00)
!  (26) Remove obsolete ASCII punch file format (bmy, 12/12/00)
!  (27) Added monoterpenes to ND46 diagnostic (bmy, 1/2/01)
!  (28) Added CO-sources from both isoprene and monoterpenes to ND29 
!        diagnostic (bmy, 1/2/01)
!  (29) Now archive biofuel emissions in ND34 diagnostic,  Also removed
!        obsolete code from 1/2/01. (bmy, 3/19/01)
!  (30) Now write biomass and biofuel tracers to the binary punch file
!        in the same order as they are listed in "diag.dat".  Also made
!        a few cosmetic changes. (bmy, 4/17/01)
!  (31) Use proper scale factor for ND21, depending on which met field
!        is being used.  Now use proper SCALEDYN for ND21 for CO simulation
!        with OH parameterization.  Also updated comments.(bmy, 4/23/01)
!  (32) Added code updates for multi-tracer Ox run (amf, bmy, 7/3/01)
!  (33) GEOS-2, GEOS-3 now use SCALECHEM for ND21 (bmy, 8/13/01)
!  (33) Now write acetone source/sink diagnostic to disk as ND11 (bmy, 9/4/01)
!  (34) Now add error check for N > NTRACE in ND27 diagnostic (bmy, 10/22/01)
!  (35) Deleted obsolete, commented-out code from 7/01 (bmy, 11/26/01)
!  (36) Replace LAIREMS w/ LLTROP since they are both equal (bmy, 2/14/02)
!  (37) Eliminated obsolete code from 1/02.  Monthly mean optical depths for 
!        ND21 diagnostic are now N=4, N=5; don't divide these by SCALEX.
!        Also added HO2, NO2 to ND43 chemical diagnostic.  Updated comments
!        for ND21 diagnostic w/ new aerosol tracers.  Now select correct
!        unit string for ND21 tracers. (rvm, bmy, 3/1/02)
!  (38) Now use category "COLUMN-T" instead of "TROPO-AV" for ND33, since
!        that is a column sum of tracer over the whole atmosphere, not just
!        the tropopshere. (bmy, 4/3/02)
!  (39) Now replace the file unit number "12" with a parameter IU_BPCH from
!        file mod (IU_BPCH is currently set to 12 for backwards compatibility)
!        Now merged "diag7.f" into "diag3.f" as the ND62 diagnostic.
!        Also added ND01 and ND02 diagnostics for Rn-Pb_Be simulation.
!        Also convert scale factors to REAL*8 precision.  Bug fix: set units
!        of ND44 drydep flux to kg/s for Rn-Pb-Be simulation. (bmy, 8/7/02)
!  (40) Bug fix: Save levels 1:LD13 for ND13 diagnostic for diagnostic
!        categories "SO2-AC-$" and "SO2-EV-$".  Now reference F90 module
!        "tracerid_mod.f".  Now reference NUMDEP from "drydep_mod.f".
!        Now save anthro, biofuel, biomass NH3 in ND13; also fixed ND13
!        tracer numbers.  For ND13, change scale factor from SCALESRCE to 1.
!        Now references "wetscav_mod.f".  Now also save true tracer numbers 
!        for ND38 and ND39 diagnostic.  Now also write out biomass SO2.
!        Now convert ND01, ND02, ND44 diagnostics for Rn/Pb/Be from kg to 
!        kg/s here. (bmy, 1/24/03)
!  (41) Now save out natural NH3 in ND13 as "NH3-NATU" (rjp, bmy, 3/23/03)
!  (42) Now replace DXYP(JREF) by routine GET_AREA_M2, GET_XOFFSET, and
!        GET_YOFFSET of "grid_mod.f".  Now references "time_mod.f".
!        DIAGb, DIAGe are now local variables.  Now remove obsolete statements
!        IF ( LBPNCH > 0 ).  Removed SCALE1, replaced with SCALEDYN. 
!        (bmy, 2/24/03)
!  (43) Added TSKIN, PARDF, PARDR, GWET to ND67 diagnostic.  For GEOS-4/fvDAS,
!        UWND, VWND, TMPU, SPHU are A-6 fields.  Adjust the ND66 scale factors 
!        accordingly.  Delete KZZ from ND66.  Updated comments. (bmy, 6/23/03)
!  (44) Bug fix: use LD68 instead of ND68 in DO-loop to avoid out-of-bounds 
!        error. (bec, bmy, 7/15/03)
!  (45) Now print out NTRACE drydep fluxes for tagged Ox.  Also tagged Ox 
!        now saves drydep in molec/cm2/s.  Now print out Kr85 prod/loss in 
!        ND03. (bmy, 8/20/03)
!******************************************************************************
! 
      ! References to F90 modules
      USE BPCH2_MOD
      USE BIOMASS_MOD, ONLY : BIOTRCE
      USE BIOFUEL_MOD, ONLY : NBFTRACE, BFTRACE
      USE DAO_MOD,     ONLY : LWI
      USE DIAG_MOD
      USE DRYDEP_MOD,  ONLY : NUMDEP,   NTRAIND
      USE FILE_MOD,    ONLY : IU_BPCH
      USE GRID_MOD,    ONLY : GET_AREA_M2, GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD
      USE TRACERID_MOD
      USE WETSCAV_MOD, ONLY : GET_WETDEP_NSOL, GET_WETDEP_IDWETD  

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! STT, P, T, NSRCX, etc...
#     include "CMN_DIAG"   ! Diagnostic switches & arrays
#     include "CMN_O3"     ! FMOL, XNUMOL
#     include "CMN_SETUP"  ! LSPLIT
#     include "comode.h"   ! IDEMS

      ! Local variables
      INTEGER            :: I, IREF, J, JREF, L, M, MM
      INTEGER            :: N, NN, NMAX, NTEST
      INTEGER            :: IE, IN, IS, IW, ITEMP(3)
      REAL*8             :: SCALE_TMP(IIPAR,JJPAR)
      REAL*8             :: SCALE_I6,  SCALE_A6,  SCALE_A3,  SCALED    
      REAL*8             :: SCALEDYN,  SCALECONV, SCALESRCE, SCALECHEM 
      REAL*8             :: SCALEX,    SECONDS,   PMASS,     PRESSX
      REAL*8             :: FDTT,      AREA_M2,   DIAGb,     DIAGe
      
      ! For binary punch file, version 2.0
      CHARACTER (LEN=40) :: CATEGORY 
      REAL*4             :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*4             :: LONRES, LATRES
      INTEGER            :: IFIRST, JFIRST, LFIRST
      INTEGER, PARAMETER :: HALFPOLAR = 1
      INTEGER, PARAMETER :: CENTER180 = 1
      CHARACTER (LEN=20) :: MODELNAME 
      CHARACTER (LEN=40) :: UNIT
      CHARACTER (LEN=40) :: RESERVED = ''
!
!******************************************************************************
!  DIAG3 begins here!
!
!  Define scale factors for division.  
!  Add a small number (e.g. 1d-20) to prevent division by zero errors.
!******************************************************************************
!
      ! Now use counter variables from "time_mod.f" (bmy, 3/27/03)
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      SECONDS   = ( DIAGe - DIAGb ) * 3600d0
      SCALED    = 1d0
      SCALEDYN  = DBLE( GET_CT_DYN()  ) + TINY( 1d0 )
      SCALECONV = DBLE( GET_CT_CONV() ) + TINY( 1d0 )
      SCALESRCE = DBLE( GET_CT_EMIS() ) + TINY( 1d0 )
      SCALECHEM = DBLE( GET_CT_CHEM() ) + TINY( 1d0 )
      SCALE_A3  = DBLE( GET_CT_A3()   ) + TINY( 1d0 )
      SCALE_A6  = DBLE( GET_CT_A6()   ) + TINY( 1d0 )
      SCALE_I6  = DBLE( GET_CT_I6()   ) + TINY( 1d0 )
!
!******************************************************************************
!  Setup for binary punch file:
!
!  IFIRST, JFIRST, LFIRST = I, J, L indices of the starting grid box 
!  LONRES                 = DISIZE, cast to REAL*4
!  LATRES                 = DJSIZE, cast to REAL*4
!******************************************************************************
!
      IFIRST = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LFIRST = 1
      LONRES = DISIZE
      LATRES = DJSIZE

      ! Get the proper model name for the binary punch file
      MODELNAME = GET_MODELNAME()
!
!******************************************************************************
!  ND01: Rn, Pb, Be emissions (Category: "RN--SRCE")
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  Rn222  : Emissions of 222Rn             : kg/s       : SCALESRCE
!  (2)  Pb210  : Emissions of 210Pb             : kg/s       : SCALECHEM
!  (3)  Be7    : Emissions of 7Be               : kg/s       : SCALESRCE
!******************************************************************************
!
      IF ( ND01 > 0 ) THEN
         CATEGORY = 'RN--SRCE'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(1)
            N  = TINDEX(1,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET
               
            ! Pb "emission" comes from chemical decay of Rn, which happens 
            ! in the chemistry routine, so use SCALECHEM (bmy, 1/27/03)
            IF ( N == IDTPB ) THEN
               SCALEX = SCALECHEM
            ELSE
               SCALEX = SCALESRCE
            ENDIF
                 
            ! Divide by # of emission timesteps
            DO L = 1, LD01
               ARRAY(:,:,L) = AD01(:,:,L,N) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD01,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD01) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND02: Rn, Pb, Be lost to radioactive decay (Category: "RN-DECAY")
!
!   # : Field  : Description                   : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  Rn222  : Loss of 222Rn                 : kg/s       : SCALECHEM
!  (2)  Pb210  : Loss of 210Pb                 : kg/s       : SCALECHEM
!  (3)  Be7    : Loss of 7Be                   : kg/s       : SCALECHEM
!******************************************************************************
!
      IF ( ND02 > 0 ) THEN
         CATEGORY = 'RN-DECAY'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(2)
            N  = TINDEX(2,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET

            ! Divide by # of chemistry timesteps
            DO L = 1, LD02
               ARRAY(:,:,L) = AD02(:,:,L,N) / SCALECHEM
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD02,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD02) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND03: Rn, Pb, Be emissions (Category: "RN--SRCE")
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  P(Kr)  : Emissions of Kr85              : kg         : SCALESRCE
!  (2)  L(Kr)  : Decay of Kr85                  : kg         : SCALECHEM
!******************************************************************************
!
      IF ( ND03 > 0 ) THEN
         CATEGORY = 'KRBUDGET'
         UNIT     = 'kg'

         DO M = 1, TMAX(3)
            N  = TINDEX(3,M)
            IF ( N > PD03 ) CYCLE
            NN = N
               
            ! Pick proper scale factor
            IF ( N == 1 ) THEN
               SCALEX = SCALESRCE
            ELSE
               SCALEX = SCALECHEM
            ENDIF

            ! Divide by # of emission timesteps
            DO L = 1, LD03
               ARRAY(:,:,L) = AD03(:,:,L,N) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD03,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD03) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND05: Production/Loss for coupled fullchem/aerosol runs (NSRCX==3) or
!        offline sulfate chemistry runs (NSRCX==10).      
!
!   # : Field  : Description                   : Units        : Scale factor
!  ----------------------------------------------------------------------------
!  (1 ) SO2dms : P(SO2) from DMS + OH          : kg S         : SCALEX
!  (2 ) SO2no3 : P(SO2) from DMS + NO3         : kg S         : SCALEX
!  (3 ) SO2    : Total P(SO2)                  : kg S         : SCALEX
!  (4 ) MSAdms : P(MSA) from DMS               : kg S         : SCALEX
!  (5 ) SO4gas : P(SO4) gas phase              : kg S         : SCALEX
!  (6 ) SO4aq  : P(SO4) aqueous phase          : kg S         : SCALEX
!  (7 ) SO4    : Total P(SO4)                  : kg S         : SCALEX
!  (8 ) LOH    : L(OH) by DMS                  : kg OH        : SCALEX
!  (9 ) LNO3   : L(NO3) by DMS                 : kg NO3       : SCALEX
!  (10) LH2O2  : L(H2O2)                       : kg H2O2      : SCALEX
!******************************************************************************
!
      IF ( ND05 > 0 ) THEN
         CATEGORY = 'PL-SUL=$'

         DO M = 1, TMAX(5)
            N = TINDEX(5,M)

            ! Tracers 8, 9, 10 are OH, NO3, H2O2,
            ! and are in [kg] instead of [kg S]
            IF ( N < 8 ) THEN 
               UNIT = 'kg S'
            ELSE
               UNIT = 'kg'
            ENDIF

            NN     = N + TRCOFFSET
            SCALEX = 1.d0

            DO L = 1, LD05
               ARRAY(:,:,L) = AD05(:,:,L,N) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD05,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD05) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND11: Acetone source & sink diagnostic (Category: "ACETSRCE")
!
!   # : Field  : Description                        : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  ACETmo : Acetone source from MONOTERPENES   : at C/cm2/s : SCALESRCE
!  (2)  ACETmb : Acetone source from METHYL BUTENOL : at C/cm2/s : SCALESRCE 
!  (3)  ACETbg : Acetone source from DIRECT EMISSION: at C/cm2/s : SCALESRCE 
!  (4)  ACETdl : Acetone source from DRY LEAF MATTER: at C/cm2/s : SCALESRCE 
!  (5)  ACETgr : Acetone source from GRASSLANDS     : at C/cm2/s : SCALESRCE 
!  (6)  ACETop : Acetone source from OCEANS         : at C/cm2/s : SCALESRCE 
!  (7)  ACETol : Acetone sink   from OCEANS         : at C/cm2/s : SCALECHEM
!******************************************************************************
!
      IF ( ND11 > 0 ) THEN
         CATEGORY = 'ACETSRCE'
         UNIT     = 'molC/cm2/s'

         DO M = 1, TMAX(11)
            N  = TINDEX(11,M)
            IF ( N > PD11 ) CYCLE
            NN = N 
               
            ! Acetone ocean sink is on the chemistry timestep
            ! but acetone sources are all on the emission timestep
            IF ( N == 7 ) THEN
               SCALEX = SCALECHEM
            ELSE
               SCALEX = SCALESRCE
            ENDIF

            ARRAY(:,:,1) = AD11(:,:,N) / SCALEX

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND12: distribution of suface emissions in the boundry layer: [fraction]
!
!   # : Field   : Description                         : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1) EMDIS-BL : Fraction of BL occupied by level L  : unitless : SCALECHEM
!******************************************************************************
!
      IF ( ND12 > 0 ) THEN
         UNIT     = 'unitless'
         CATEGORY = 'EMDIS-BL'

         DO L = 1, LD12
            ARRAY(:,:,L) = AD12(:,:,L) / SCALECHEM
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 1,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LLTROP,   IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD12) )
      ENDIF
!
!******************************************************************************
!  ND13: Sulfur emissions (for DMS/SO2/SO4/MSA/NH3/NH4/NIT chemistry)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) DMS-BIOG : Biogenic DMS emission           : kg S     : 1
!  (2 ) SO2-AC-$ : Aircraft SO2 emission           : kg S     : 1
!  (3 ) SO2-AN-$ : Anthropogenic SO2 emission      : kg S     : 1
!  (4 ) SO2-BIOB : Biomass SO2 emission            : kg S     : 1
!  (5 ) SO2-BIOF : Biofuel SO2 emission            : kg S     : 1
!  (6 ) SO2-NV-$ : Non-eruptive volcano SO2 em.    : kg S     : 1
!  (7 ) SO2-EV-$ : Eruptive volcano SO2 emissions  : kg S     : 1
!  (8 ) SO4-AN-$ : Anthropogenic SO4 emission      : kg S     : 1
!  (9 ) NH3-ANTH : Anthropogenic NH3 emission      : kg NH3   : 1
!  (10) NH3-NATU : Natural source NH3 emission     : kg NH3   : 1
!  (11) NH3-BIOB : Biomass burning NH3 emission    : kg NH3   : 1
!  (12) NH3-BIOF : Biofuel burning NH3 emission    : kg NH3   : 1
!******************************************************************************
!
      IF ( ND13 > 0 .and. ( NSRCX == 3 .or. NSRCX == 10 ) ) THEN
         UNIT = 'kg S'

         !==============================================================
         ! Biogenic DMS 
         !==============================================================
         CATEGORY     = 'DMS-BIOG'
         ARRAY(:,:,1) = AD13_DMS(:,:)
         N            = IDTDMS + TRCOFFSET
         
         CALL BPCH2( 12,        MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Aircraft SO2 
         !==============================================================
         CATEGORY = 'SO2-AC-$'
         N        = IDTSO2 + TRCOFFSET

         DO L = 1, LD13
            ARRAY(:,:,L) = AD13_SO2_ac(:,:,L)
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N, 
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD13,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD13) )

         !==============================================================
         ! Anthropogenic SO2 
         !==============================================================
         CATEGORY = 'SO2-AN-$'
         N        = IDTSO2 + TRCOFFSET
         
         DO L = 1, NOXEXTENT 
            ARRAY(:,:,L) = AD13_SO2_an(:,:,L)
         ENDDO
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &               IIPAR,     JJPAR,     NOXEXTENT, IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:NOXEXTENT) )

         !==============================================================
         ! Biomass SO2 
         !==============================================================
         CATEGORY     = 'SO2-BIOB'
         ARRAY(:,:,1) = AD13_SO2_bb(:,:)
         N            = IDTSO2 + TRCOFFSET

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )


         !==============================================================
         ! Biofuel SO2
         !==============================================================
         CATEGORY     = 'SO2-BIOF'
         ARRAY(:,:,1) = AD13_SO2_bf(:,:)
         N            = IDTSO2 + TRCOFFSET

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Eruptive volcano SO2 
         !==============================================================
         CATEGORY = 'SO2-EV-$'
         N        = IDTSO2 + TRCOFFSET

         DO L = 1, LD13
            ARRAY(:,:,L) = AD13_SO2_ev(:,:,L)
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD13,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD13) )

         !==============================================================
         ! Non-eruptive volcano SO2 
         !==============================================================
         CATEGORY = 'SO2-NV-$'
         N        = IDTSO2 + TRCOFFSET

         DO L = 1, LD13 
            ARRAY(:,:,L) = AD13_SO2_nv(:,:,L)
         ENDDO
               
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD13,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD13) )

         !==============================================================
         ! Anthropogenic SO4 
         !==============================================================
         CATEGORY = 'SO4-AN-$'
         N        = IDTSO4 + TRCOFFSET

         DO L = 1, NOXEXTENT 
            ARRAY(:,:,L) = AD13_SO4_an(:,:,L)
         ENDDO
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &               IIPAR,     JJPAR,     NOXEXTENT, IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:NOXEXTENT) )

         !==============================================================
         ! Anthropogenic NH3 
         !==============================================================
         UNIT         = 'kg'
         CATEGORY     = 'NH3-ANTH'
         ARRAY(:,:,1) = AD13_NH3_an(:,:) 
         N            = IDTNH3 + TRCOFFSET

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Natural source NH3
         !==============================================================
         CATEGORY     = 'NH3-NATU'
         ARRAY(:,:,1) = AD13_NH3_na(:,:) 
         N            = IDTNH3 + TRCOFFSET

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Biomass NH3
         !==============================================================
         CATEGORY     = 'NH3-BIOB'
         ARRAY(:,:,1) = AD13_NH3_bb(:,:)
         N            = IDTNH3 + TRCOFFSET

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Biofuel NH3 
         !==============================================================
         CATEGORY     = 'NH3-BIOF'
         ARRAY(:,:,1) = AD13_NH3_bf(:,:)
         N            = IDTNH3 + TRCOFFSET

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDIF    
!
!******************************************************************************
!  ND14: Upward mass flux from wet convection (NFCLDMX)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  CONVFLUP : Upward mass flux from wet conv  : kg/s     : SCALECONV
!
!  NOTES:
!  (1) Bug fix -- only write LD14 levels to the bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND14 > 0 ) THEN
         CATEGORY = 'CV-FLX-$'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(14)
            N  = TINDEX(14,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET
               
            ARRAY(:,:,1:LD14) = CONVFLUP(:,:,1:LD14,N) / SCALECONV

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD14,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD14) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND15: Upward mass flux from boundary layer mixing (TURBDAY)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  TURBFLUX : Upward mass flux from BL mixing : kg/s     : SCALECONV
!
!  NOTES:
!  (1) Bug fix -- only write LD15 levels to the bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND15 > 0 ) THEN
         CATEGORY = 'TURBMC-$'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(15)
            N  = TINDEX(15,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET

            ARRAY(:,:,1:LD15) = TURBFLUP(:,:,1:LD15,N) / SCALECONV
               
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD15,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD15) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND16: Fraction of grid box experiencing LS or convective precipitation
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WD-FLS-$ : LS precip fraction          : unitless  : CT16(:,:,:,1)
!  (2) WD-FCV-$ : Convective precip fraction  : unitless  : CT16(:,:,:,2)
!******************************************************************************
!
      IF ( ND16 > 0 ) THEN

         ! Large-scale area of precipitation
         CATEGORY = 'WD-FRC-$'
         UNIT     = 'unitless'

         DO M = 1, TMAX(16)
            N  = TINDEX(16,M)
            IF ( N > PD16 ) CYCLE
            NN = N 
              
            DO L = 1, LD16
               SCALE_TMP(:,:) = FLOAT( CT16(:,:,L,N) ) + 1d-20
               ARRAY(:,:,L)   = AD16(:,:,L,N) / SCALE_TMP(:,:)
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD16,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD16) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND17: Fraction of tracer lost rainout in LS and convective precip
!
!   #  Field    : Description                  : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WD-LSR-$ : Rainout fraction/LS Precip   : unitless  : CT17(:,:,:,1)
!  (2) WD-CVR-$ : Rainout fraction/conv precip : unitless  : CT17(:,:,:,2)
!******************************************************************************
!
      IF ( ND17 > 0 ) THEN
         UNIT = 'unitless'

         DO M = 1, TMAX(17)
            N = TINDEX(17,M)

            ! Rn-Pb-Be run has 2 soluble tracers
            ! Full chemistry run has 4 soluble tracers
            IF ( NSRCX == 1 .and. N > 2    ) CYCLE
            IF ( NSRCX == 3 .and. N > PD17 ) CYCLE

            NN = N + TRCOFFSET

            ! Large-scale rainout/washout fractions
            CATEGORY = 'WD-LSR-$'
               
            DO L = 1, LD17
               SCALE_TMP(:,:) = FLOAT( CT17(:,:,L,1) ) + 1d-20
               ARRAY(:,:,L)   = AD17(:,:,L,N,1) / SCALE_TMP(:,:) 
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD17,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD17) )


            ! Convective rainout/washout fractions
            CATEGORY = 'WD-CVR-$'

            DO L = 1, LD17
               SCALE_TMP(:,:) = FLOAT( CT17(:,:,L,2) ) + 1d-20
               ARRAY(:,:,L)   = AD17(:,:,L,N,2) / SCALE_TMP(:,:) 
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD17,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD17) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND18: Fraction of tracer lost to washout in LS or convective precip
!
!   #  Field    : Description                  : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WD-LSW-$ : Washout fraction/LS precip   : unitless  : CT18(:,:,:,1)
!  (2) WD-CVW-$ : Washout fraction/conv precip : unitless  : CT18(:,:,:,2)
!******************************************************************************
!
      IF ( ND18 > 0 ) THEN
         UNIT = 'unitless'

         DO M = 1, TMAX(18)
            N = TINDEX(18,M)

            ! Rn-Pb-Be run has 2 soluble tracers
            ! Full chemistry run has 4 soluble tracers
            IF ( NSRCX == 1 .and. N > 2    ) CYCLE
            IF ( NSRCX == 3 .and. N > PD18 ) CYCLE

            NN = N + TRCOFFSET

            ! Large-scale rainout/washout fractions
            CATEGORY = 'WD-LSW-$'

            DO L = 1, LD18
               SCALE_TMP(:,:) = FLOAT( CT18(:,:,L,1) ) + 1d-20
               ARRAY(:,:,L)   = AD18(:,:,L,N,1) / SCALE_TMP(:,:) 
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD18,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD18) )


            ! Convective washout fractions
            CATEGORY = 'WD-CVW-$'

            DO L = 1, LD18
               SCALE_TMP(:,:) = FLOAT( CT18(:,:,L,2) ) + 1d-20
               ARRAY(:,:,L)   = AD18(:,:,L,N,2) / SCALE_TMP(:,:) 
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD18,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD18) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND21: Optical depth diagnostics
!
!   # : Field : Description                        : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) OPTD   Cloud Optical Depth                 : unitless : SCALECHEM
!  (2 ) CLMO   Maximum Overlap Cloud Fraction      : unitless : SCALECHEM
!  (3 ) CLRO   Random  Overlap Cloud Fraction      : unitless : SCALECHEM
!  (4 ) OPD    Mineral Dust Optical Depth (400 nm) : unitless : none
!  (5 ) SD     Mineral Dust Surface Area           : cm2/cm3  : none
!  (6 ) OPSO4  Sulfate Optical Depth (400 nm)      : unitless : SCALECHEM
!  (7 ) HGSO4  Hygroscopic growth of SO4           : unitless : SCALECHEM
!  (8 ) SSO4   Sulfate Surface Area                : cm2/cm3  : SCALECHEM
!  (9 ) OPBC   Black Carbon Optical Depth (400 nm) : unitless : SCALECHEM
!  (10) HGBC   Hygroscopic growth of BC            : unitless : SCALECHEM
!  (11) SBC    Black Carbon Surface Area           : cm2/cm3  : SCALECHEM
!  (12) OPOC   Organic C Optical Depth (400 nm)    : unitless : SCALECHEM
!  (13) HGOC   Hygroscopic growth of OC            : unitless : SCALECHEM
!  (14) SOC    Organic Carbon Surface Area         : cm2/cm3  : SCALECHEM
!  (15) OPSSa  Sea Salt (accum) Opt Depth (400 nm) : unitless : SCALECHEM
!  (16) HGSSa  Hygroscopic growth of SSa           : unitless : SCALECHEM
!  (17) SSSa   Sea Salt (accum) Surface Area       : cm2/cm3  : SCALECHEM
!  (18) OPSSc  Sea Salt (coarse) Opt Depth(400 nm) : unitless : SCALECHEM
!  (19) HGSSc  Hygroscopic growth of SSc           : unitless : SCALECHEM
!  (20) SSSc   Sea Salt (coarse) Surface Area      : cm2/cm3  : SCALECHEM  
!
!  NOTES:
!  (1 ) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2 ) Don't divide monthly mean AOD by SCALECHEM (rvm, bmy, 12/8/00)
!  (3 ) Use SCALE_A6 for GEOS-2, GEOS-3 fields, since optical depths are read
!        in from disk every 6 hours.  Use SCALECHEM for GEOS-1, GEOS-STRAT
!        fields, since optical depths are computed every chemistry timestep.
!        Use SCALEDYN for CO-OH parameterization simulation. (bmy, 4/23/01)
!  (4 ) Now GEOS-2, GEOS-3 use SCALECHEM for ND21 (bmy, 8/13/01)
!  (5 ) Updated tracers for new aerosols from Mian Chin (rvm, bmy, 3/1/02)
!******************************************************************************
!
      IF ( ND21 > 0 ) THEN
         CATEGORY = 'OD-MAP-$'

         ! For CO-OH param run,   ND21 is updated every dynamic   timestep
         ! For other simulations, ND21 is updated every chemistry timestep
         IF ( NSRCX == 5 ) THEN
            SCALEX = SCALEDYN
         ELSE
            SCALEX = SCALECHEM
         ENDIF

         DO M = 1, TMAX(21)
            N  = TINDEX(21,M)
            IF ( N > PD21 ) CYCLE
            NN = N 
               
            ! Select proper unit string (cf list above)
            SELECT CASE( N ) 
               CASE ( 5, 8, 11, 14, 17, 20 )
                  UNIT = 'cm2/cm3'
               CASE DEFAULT
                  UNIT = 'unitless'
            END SELECT
               
            ! The monthly mean dust optical depths (N > 3)
            ! don't need to be scaled by SCALECHEM (rvm, bmy, 12/8/00)
            IF ( N > 3 .AND. N < 6 ) THEN
               ARRAY(:,:,1:LD21) = AD21(:,:,1:LD21,N)
            ELSE
               ARRAY(:,:,1:LD21) = AD21(:,:,1:LD21,N) / SCALEX
            ENDIF
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD21,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD21) )
         ENDDO    
      ENDIF
!
!******************************************************************************
!  ND22: J-value diagnostics
!
!   #  : Field : Description                   : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1  ) JNO2  : NO2   J-Value                 : s-1   : SCALE_JV
!  (2  ) JHNO3 : HNO3  J-Value                 : s-1   : SCALE_JV
!  (3  ) JH2O2 : H2O2  J-Value                 : s-1   : SCALE_JV
!  (4  ) JCH2O : CH2O  J-Value                 : s-1   : SCALE_JV
!  (5  ) JO3   : O3    J-Value                 : s-1   : SCALE_JV
!  (6  ) POH   : OH-source from O3 photolysis  : s-1   : SCALE_JV
!  (71 ) JCH3I : CH3I  J-value (s^-1)          : s-1   : SCALE_JV
!  (81 ) JHCN  : HCN   J-value (s^-1)          : s-1   : SCALE_JV
!
!  NOTES:
!  (1) We must add TRCOFFSET for CH3I and HCN runs, so that GAMAP can
!       recognize those photo rates as distinct from the NO2, HNO3,
!       H2O2, CH2O, O3, and POH photo rates.
!******************************************************************************
!
      IF ( ND22 > 0 ) THEN
         CATEGORY  = 'JV-MAP-$'
         SCALE_TMP = FLOAT( CTJV ) + 1d-20
         UNIT      = 's-1'

         DO M = 1, TMAX(22)
            N  = TINDEX(22,M)
            IF ( N > PD22 ) CYCLE
            
            ! Add TRCOFFSET to CH3I and HCN tracer numbers
            IF ( NSRCX == 2 .or. NSRCX == 4 ) THEN
               NN = N + TRCOFFSET
            ELSE
               NN = N
            ENDIF

            DO L = 1, LD22
               ARRAY(:,:,L) = AD22(:,:,L,N) / SCALE_TMP(:,:)
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD22,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD22) )
         ENDDO    
      ENDIF     
!
!******************************************************************************
!  ND24: Eastward mass flux from transport (TPCORE, XTP)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  MASSFLEW : Eastward mass flux - transport  : kg/s     : SCALEDYN
!
!  NOTES:
!  (1) MASSFLEW is REAL*8...store to ARRAY, which is REAL*4
!       before sending to BPCH or IJSCAL (bey, bmy, 4/23/99)
!  (2) Now only write LD24 levels out to the bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND24 > 0 ) THEN
         CATEGORY = 'EW-FLX-$'
         UNIT = 'kg/s'

         DO M = 1, TMAX(24)
            N  = TINDEX(24,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET
            
            ARRAY(:,:,1:LD24) = MASSFLEW(:,:,1:LD24,N) / SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD24,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD24) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND25: Northward mass flux from transport (TPCORE, YTP)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  MASSFLNS : Northward mass flux - transport : kg/s     : SCALEDYN
!
!  NOTES:
!  (1) MASSFLNS is REAL*8...store to ARRAY, which is REAL*4
!       before sending to BPCH or IJSCAL (bey, bmy, 4/23/99)
!  (2) Now only write LD25 levels out to the bpch file (bmy, 12/7/00)
!******************************************************************************
!  
      IF ( ND25 > 0 ) THEN
         CATEGORY = 'NS-FLX-$'
         UNIT = 'kg/s'

         DO M = 1, TMAX(25)
            N  = TINDEX(25,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET

            ARRAY(:,:,1:LD25) = MASSFLNS(:,:,1:LD25,N) / SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD25,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD25) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND26: Upward mass flux from transport (TPCORE, FZPPM)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  MASSFLUP : Upward mass flux - transport    : kg/s     : SCALEDYN
!
!  NOTES:
!  (1) MASSFLNS is REAL*8...store to ARRAY, which is REAL*4
!       before sending to BPCH or IJSCAL (bey, bmy, 4/23/99)
!  (2) Now only write LD26 levels to the bpch file (bmy, 12/7/00)
!******************************************************************************
!  
      IF ( ND26 > 0 ) THEN
         CATEGORY = 'UP-FLX-$'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(26)
            N  = TINDEX(26,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET
            
            ARRAY(:,:,1:LD26) = MASSFLUP(:,:,1:LD26,N) / SCALEDYN
               
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES, 
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD26,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD26) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND27: Cross-tropopause Stratospheric Influx of Ox 
!
!   #  : Field : Description                   : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1) : Ox    : Ox from the stratosphere      : kg/s  : SCALEDYN
!
!  NOTES:
!  (1) Only print out if we are doing a NOx-Ox-HC run (NSRCX == 3)
!       or a single tracer Ox run (NSRCX == 6). (bey, bmy, 11/10/99)
!  (2) Now consider the cross-tropopause stratospheric influx of ozone, 
!       which, in some grid boxes, includes horizontal influxes as well as 
!       up(down)ward flux. (qli, 1/5/2000) 
!  (3) Now error check for N > NTRACE (bmy, 10/23/01)
!******************************************************************************
!
      IF ( ND27 > 0 .and. IDTOX > 0 ) then
         IF ( NSRCX == 3 .or. NSRCX == 6 ) THEN
            CATEGORY = 'STRT-FLX'
            UNIT     = 'kg/s'

            ! Full chemistry   -- compute NOx, Ox, HNO3 fluxes
            ! Single tracer Ox -- compute Ox flux only, hardwire
            !                     to tracer = 1 (bmy, 2/7/00)
            IF ( NSRCX == 3 ) THEN
               ITEMP = (/ IDTNOX, IDTOX, IDTHNO3 /)
            ELSE
               ITEMP = (/ 1, 0, 0 /)
            ENDIF
            
            ! Loop over tracers
            DO M = 1, 3
               N = ITEMP(M)
               IF ( N == 0     ) CYCLE
               IF ( N > NTRACE ) CYCLE

               ! Loop over grid boxes
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Get the level of the tropopause
                  L = LPAUSE(I,J)

                  ! Initialize integer flags
                  IS = 0
                  IN = 0
                  IW = 0
                  IE = 0

                  ! Set integer flags based on the value of each bit of IFLX
                  IF ( BTEST( IFLX(I,J), 0 ) ) IS = 1
                  IF ( BTEST( IFLX(I,J), 1 ) ) IN = 1
                  IF ( BTEST( IFLX(I,J), 2 ) ) IW = 1
                  IF ( BTEST( IFLX(I,J), 3 ) ) IE = 1

                  ! Add fluxes from the top, south, and west
                  ARRAY(I,J,1) = MASSFLUP(I,J,L,N)          +
     &                           ( MASSFLNS(I,J,L,N) * IS ) +  
     &                           ( MASSFLEW(I,J,L,N) * IW ) 
               
                  ! Add fluxes from the north 
                  ! (take poles into account !)
                  IF ( J < JJPAR ) THEN
                     ARRAY(I,J,1) = ARRAY(I,J,1) -
     &                              ( MASSFLNS(I,J+1,L,N) * IN )
                  ELSE 
                     ARRAY(I,J,1) = ARRAY(I,J,1) -
     &                              ( MASSFLNS(I,  1,L,N) * IN )
                  ENDIF

                  ! Add fluxes from the east 
                  !(wrap around dateline if necessary)
                  IF ( I < IIPAR ) THEN
                     ARRAY(I,J,1) = ARRAY(I,J,1) -
     &                              ( MASSFLEW(I+1,J,L,N) * IE )
                  ELSE 
                     ARRAY(I,J,1) = ARRAY(I,J,1) -
     &                             ( MASSFLEW(  1,J,L,N) * IE )
                  ENDIF
               ENDDO
               ENDDO

               UNIT = 'kg/s'
               NN   = N + TRCOFFSET
                  
               ARRAY(:,:,1) = ARRAY(:,:,1) / SCALEDYN

               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     PD27,     IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )
            ENDDO
         ENDIF
      ENDIF
!
!******************************************************************************
!  ND28: Biomass burning diagnostic 
!
!   # : Field : Description   : Units            : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) NOx   : NOx           : molec NOx /cm2/s : SCALESRCE
!  (4 ) CO    : CO            : molec CO  /cm2/s : SCALESRCE
!  (9 ) ACET  : Acetone       : atoms C   /cm2/s : SCALESRCE
!  (10) MEK   : Ketones(>C3)  : atoms C   /cm2/s : SCALESRCE
!  (11) ALD2  : Acetaldehyde  : atoms C   /cm2/s : SCALESRCE
!  (18) PRPE  : Propene       : atoms C   /cm2/s : SCALESRCE
!  (19) C3H8  : Propane       : atoms C   /cm2/s : SCALESRCE
!  (20) C2HO  : Formaldehyde  : molec CH2O/cm2/s : SCALESRCE
!  (21) C2H6  : Ethane        : atoms C   /cm2/s : SCALESRCE
!
!  NOTES:
!  (1) Use the F90 intrinsic "ANY" function to make sure that N 
!       corresponds to actual biomass burning tracers (bmy, 4/8/99)
!  (2) ND28 now uses allocatable array AD28 instead of AIJ. (bmy, 3/16/00)
!  (3) Now write biofuel burning tracers to the punch file in the same order 
!       as they are listed in "diag.dat". (bmy, 4/17/01)
!******************************************************************************
!
      IF ( ND28 > 0 ) THEN
         CATEGORY = 'BIOBSRCE'
         UNIT     = ''

         DO M = 1, TMAX(28)
            N  = TINDEX(28,M)
            IF ( .not. ANY( BIOTRCE == N ) ) CYCLE
            NN = N + TRCOFFSET
            
            ARRAY(:,:,1) = AD28(:,:,M) / SCALESRCE
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND29: CO source diagnostics
!
!   #  Field  : Description             : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) COanth : CO from Anthro Sources  : mol/cm2/s : SCALESRCE
!  (2) CObb   : CO from Biomass Burning : mol/cm2/s : SCALESRCE
!  (3) CObf   : CO from Biofuel Burning : mol/cm2/s : SCALESRCE
!  (4) COmeth : CO from Methanol        : mol/cm2/s : SCALESRCE
!  (5) COmono : CO from Monoterpenes    : mol/cm2/s : SCALESRCE
!
!  NOTES:
!  (1) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2) ND29 now uses allocatable array AD29 instead of AIJ. (bmy, 3/16/00)
!  (3) Added CO-sources from isoprene and monoterpenes (bnd, bmy, 1/2/01)
!******************************************************************************
!
      IF ( ND29 > 0 ) THEN
         CATEGORY ='CO--SRCE'
         UNIT    = 'mol/cm2/s'

         DO M = 1, TMAX(29)
            N  = TINDEX(29,M)
            IF ( N > PD29 ) CYCLE
            NN = N 

            ARRAY(:,:,1) = AD29(:,:,N) / SCALESRCE

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND30: Land map diagnostic
!
!   #  Field : Description             : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) LWI   : GMAO Land-Water indices  : unitless  : SCALED 
!
!  NOTES: 
!  (1) Values are: 1=sea, 2=land, 3=land ice, 4=sea ice
!  (2) The antarctic sea ice pack will grow and retreat during the year. 
!  (3) Remove reference to the AAAIJ array (bmy, 12/12/00)
!******************************************************************************
!
      IF ( ND30 > 0 ) THEN
         CATEGORY = 'LANDMAP'
         UNIT     = 'unitless'
            
         ARRAY(:,:,1) = FLOAT( LWI(:,:) )
         NN           = 1 

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDIF
!
!******************************************************************************
!  ND31: Surface pressure diagnostic
!
!   #  Field : Description             : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) PS    : PS - PTOP               : mb        : SCALEDYN
!
!  NOTES: 
!  (1) The ASCII punch file was using SCALE2 instead of SCALE1.
!       This has now been fixed. (hyl, bmy, 12/21/99).
!  (2) Now use AD31 dynamically allocatable array (bmy, 2/17/00)
!  (3) Bug fix: write out 1 level to the bpch file (bmy, 12/7/00)
!  (4) Now remove SCALE1, replace with SCALEDYN (bmy, 2/24/03)
!******************************************************************************
!   
      IF ( ND31 > 0 ) THEN
         CATEGORY     = 'PS-PTOP'
         UNIT         = 'mb'
         ARRAY(:,:,1) = AD31(:,:,1) / SCALEDYN
         NN           = 1

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDIF
!
!******************************************************************************
!  ND32: NOx source diagnostic
!
!  Levels        : Field                  : Units       : Scale Factor
!  -------------------------------------------------------------------------
!  1 - LLTROP    : Aircraft NOx           : molec/cm2/s : SCALESRCE
!  1 - NOXEXTENT : Anthropogenic NOx      : molec/cm2/s : SCALESRCE
!  Surface       : Biomass Burning NOx    : molec/cm2/s : SCALESRCE
!  Surface       : Biofuel Burning NOx    : molec/cm2/s : SCALESRCE
!  Surface       : Fertilizer NOx         : molec/cm2/s : SCALESRCE
!  1 - LLCONVM   : Lightning NOx          : molec/cm2/s : SCALESRCE
!  Surface       : Soil NOx               : molec/cm2/s : SCALESRCE
!  Above TP      : NOx from upper boundary: molec/cm2/s : SCALEDYN
!
!  Print out all of the types of NOx, for all levels.
!
!  NOTES:
!  (1) Only print out ND32 if for an O3 chemistry run ( NSRCX == 3 ),
!       and if NOx is a defined tracer ( IDTNOX > 0 ). (bmy, 5/26/99)
!  (2) ND32 now uses allocatable arrays instead of AIJ. (bmy 3/16/00)
!  (3) Added biofuel burning to ND32 diagnostic (bmy, 9/12/00)
!******************************************************************************
!
      IF ( ND32 > 0 .and. IDTNOX > 0 .and. NSRCX == 3 ) THEN

         ! All categories of NOx are in molec/cm2/s
         UNIT = 'molec/cm2/s'

         !==============================================================
         ! Aircraft NOx
         !==============================================================
         CATEGORY = 'NOX-AC-$'
            
         DO L = 1, LLTROP
            ARRAY(:,:,L) = AD32_ac(:,:,L) / SCALESRCE               
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LLTROP,   IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LLTROP) )

         !==============================================================
         ! Anthropogenic NOx
         !==============================================================
         CATEGORY = 'NOX-AN-$'

         DO L = 1, NOXEXTENT 
            ARRAY(:,:,L) = AD32_an(:,:,L) / SCALESRCE
         ENDDO
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY,  IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &               IIPAR,     JJPAR,     NOXEXTENT, IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:NOXEXTENT) )

         !==============================================================
         ! Biomass Burning NOx
         !==============================================================
         CATEGORY     = 'NOX-BIOB'
         ARRAY(:,:,1) = AD32_bb(:,:) / SCALESRCE

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Binary punch file: NOx from Biofuel
         !==============================================================
         CATEGORY     = 'NOX-BIOF'
         ARRAY(:,:,1) = AD32_bf(:,:) / SCALESRCE

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &                  HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Fertilizer NOx
         !==============================================================
         CATEGORY     = 'NOX-FERT'
         ARRAY(:,:,1) = AD32_fe(:,:) / SCALESRCE
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Lightning NOx
         !==============================================================
         CATEGORY = 'NOX-LI-$'

         DO L = 1, LLCONVM 
            ARRAY(:,:,L) = AD32_li(:,:,L) / SCALESRCE
         ENDDO
               
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX, 
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LLCONVM,  IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LLCONVM) )

         !==============================================================
         ! Soil NOx
         !==============================================================
         CATEGORY     = 'NOX-SOIL'
         ARRAY(:,:,1) = AD32_so(:,:) / SCALESRCE
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Stratospheric NOx (boundary condition)
         !==============================================================
         CATEGORY     = 'NOX-STRT'
         ARRAY(:,:,1) = AD32_ub(:,:) / SCALEDYN
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDIF
!
!******************************************************************************
!  ND33: Atmospheric column sum of Tracer
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) COLUMN-T : Trop. Column Sum of Tracer  : kg        : SCALEDYN
!
!  NOTES: 
!  (1) Now use dynamically allocatable array AD33 (bmy, 2/17/00)
!  (2) Rename category to COLUMN-T, since this is a column sum of tracer over
!       the entire atmosphere, not just the troposphere. (bmy, 4/3/02)
!  (3) Now replace SCALE1 with SCALEDYN (bmy, 3/27/03)
!******************************************************************************
!
      IF ( ND33 > 0 ) THEN
         CATEGORY = 'COLUMN-T'
         UNIT     = 'kg'

         DO M = 1, TMAX(33)
            N  = TINDEX(33,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET 
            
            ARRAY(:,:,1) = AD33(:,:,N) / SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF  
!
!******************************************************************************
!  ND34: Biofuel burning diagnostic 
!
!   # : Field : Description         : Units            : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) NOx   : NOx                 : molec NOx /cm2/s : SCALESRCE
!  (4 ) CO    : CO                  : molec CO  /cm2/s : SCALESRCE
!  (5 ) ALK4  : Alkanes(>C4)        : atoms C   /cm2/s : SCALESRCE
!  (9 ) ACET  : Acetone             : atoms C   /cm2/s : SCALESRCE
!  (10) MEK   : Metyl Ethyl Ketone  : atoms C   /cm2/s : SCALESRCE
!  (11) ALD2  : Acetaldehyde        : atoms C   /cm2/s : SCALESRCE
!  (18) PRPE  : Alkenes(>=C3)       : atoms C   /cm2/s : SCALESRCE
!  (19) C3H8  : Propane             : atoms C   /cm2/s : SCALESRCE
!  (20) CH2O  : Formaldehyde        : molec CH2O/cm2/s : SCALESRCE
!  (21) C2H6  : Ethane              : atoms C   /cm2/s : SCALESRCE
!
!  NOTES:
!  (1) Use the F90 intrinsic "ANY" function to make sure that N 
!       corresponds to actual biofuel burning tracers (bmy, 3/15/01)
!  (3) Now write biofuel burning tracers to the punch file in the same order 
!       as they are listed in "diag.dat". (bmy, 4/17/01)
!******************************************************************************
!
      IF ( ND34 > 0 ) THEN
         CATEGORY = 'BIOFSRCE'
         UNIT     = ''
         
         DO M = 1, TMAX(34)
            N  = TINDEX(34,M)
            IF ( .not. ANY( BFTRACE == N ) ) CYCLE
            NN = N + TRCOFFSET

            ARRAY(:,:,1) = AD34(:,:,M) / SCALESRCE

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND35: Tracer concentration at 500 mb 
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) 500-AVRG : Tracer at 500 mb            : v/v       : SCALEDYN
!
!  NOTES:
!  (1) Now use dynamically allocatable array AD35 (bmy, 2/17/00)
!  (2) Now replace SCALE1 with SCALEDYN (bmy, 2/24/03)
!******************************************************************************
!
      IF ( ND35 > 0 ) THEN
         CATEGORY = '500-AVRG'        
         UNIT     = ''

         DO M = 1, TMAX(35)
            N  = TINDEX(35,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET
               
            ARRAY(:,:,1) = AD35(:,:,N) / SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF  
!
!******************************************************************************
!  ND36: Anthropogenic source diagnostic
!
!   #   Field  : Description                     : Units         : S. Factor
!  ---------------------------------------------------------------------------
!  (1 ) NOx    : NOx                             : mol/cm2/s     : SCALE3
!  (4 ) CO     : CO                              : mol/cm2/s     : SCALE3
!  (5 ) ALK4   : Alkanes(>C4)                    : atoms C/cm2/s : SCALE3
!  (9 ) ACET   : Acetone                         : atoms C/cm2/s : SCALE3
!  (10) MEK    : Ketones(>C3)                    : atoms C/cm2/s : SCALE3
!  (18) PRPE   : Propene                         : atoms C/cm2/s : SCALE3 
!  (19) C3H8   : Propane                         : atoms C/cm2/s : SCALE3 
!  (21) C2H6   : Ethane                          : atoms C/cm2/s : SCALE3
!  (71) CH3Ioc : Methyl Iodide (oceanic source)  : ng/m2/s       : SCALE3
!  (72) CH3Ibb : Methyl Iodide (biomass burning) : ng/m2/s       : SCALE3
!  (73) CH3Iwb : Methyl Iodide (wood burning)    : ng/m2/s       : SCALE3
!  (74) CH3Irc : Methyl Iodide (rice paddies)    : ng/m2/s       : SCALE3
!  (75) CH3Iwl : Methyl Iodide (wetlands)        : ng/m2/s       : SCALE3
!
!  NOTES:
!  (1) ND36 is also used for CH3I emissions diagnostics when NSRCX=2.
!  (2) For an O3 run (NSRCX = 3, the "default" run) make sure that the 
!       tracer number N matches an entry in the IDEMS emission index 
!       array (bmy, 4/9/99)  
!  (3) Write the tracers out to the punch file in the same order as
!       they are listed in the IDEMS array.  Thus, we have to re-assign
!       N = IDEMS(M) after we test to make sure it is a valid tracer
!       number (bmy, 4/16/99)
!  (4) For a CH3I run, make sure that the tracer number N is not larger
!       than NTRACE (bmy, 4/9/99) 
!  (5) ND36 now uses the AD36 array instead of AIJ. (bmy, 3/16/00)
!******************************************************************************
!                     
      IF ( ND36 > 0 ) THEN
         SELECT CASE ( NSRCX )
            CASE ( 2 ) 
               CATEGORY = 'CH3ISRCE'
               UNIT     = 'ng/m2/s'
            CASE DEFAULT
               CATEGORY = 'ANTHSRCE'
               UNIT     = ''
         END SELECT

         DO M = 1, TMAX(36)
            N  = TINDEX(36,M)
               
            SELECT CASE ( NSRCX )
               CASE ( 2 ) 
                  MM = N
                  IF ( N > NEMANTHRO ) CYCLE 
               CASE DEFAULT
                  MM = M
                  IF ( .not. ANY( IDEMS == N ) ) CYCLE               
                  N  = IDEMS(M)
            END SELECT

            NN = N + TRCOFFSET

            ARRAY(:,:,1) = AD36(:,:,MM) / SECONDS

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND37: Fraction of tracer scavenged in convective cloud updrafts
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WETCVF-$ : Scavenging fraction         : unitless  : SCALECONV
!******************************************************************************
!
      IF ( ND37 > 0 ) THEN
         CATEGORY = 'MC-FRC-$'
         UNIT     = 'unitless'

         DO M = 1, TMAX(37)
            N = TINDEX(37,M)
               
            ! Rn-Pb-Be run has 2 soluble tracers
            ! Full chemistry run has 4 soluble tracers
            IF ( NSRCX == 1 .and. N > 2    ) CYCLE
            IF ( NSRCX == 3 .and. N > PD37 ) CYCLE

            NN = N + TRCOFFSET

            DO L = 1, LD37
               ARRAY(:,:,L) = AD37(:,:,L,M) / SCALECONV
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD37,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD37) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND38: Rainout loss of tracer in convective updrafts
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WETDCV-$ : Rainout loss of tracer      : kg/s      : SCALECONV
!
!  NOTES:
!  (1) Now write only LD38 levels to bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND38 > 0 ) THEN
         CATEGORY = 'WETDCV-$'
         UNIT     = 'kg/s'

         ! Get actual # of soluble tracers
         M = GET_WETDEP_NSOL()

         ! Loop over soluble tracers
         DO N = 1, M

            ! Tracer number plus GAMAP offset
            NN = GET_WETDEP_IDWETD(N) + TRCOFFSET

            ! Divide by # of convective timesteps
            DO L = 1, LD38
               ARRAY(:,:,L) = AD38(:,:,L,N) / SCALECONV
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD38,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD38) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND39: Rainout loss of tracer in large scale rains 
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WETDLS-$ : Large-scale loss of tracer  : kg/s      : SCALEDYN
!******************************************************************************
!
      IF ( ND39 > 0 ) THEN
         CATEGORY = 'WETDLS-$'
         UNIT     = 'kg/s'

         ! Get actual # of soluble tracers
         M = GET_WETDEP_NSOL()
            
         ! Loop over soluble tracers
         DO N = 1, M
               
            ! Tracer number plus GAMAP offset
            NN = GET_WETDEP_IDWETD(N) + TRCOFFSET

            ! Divide by # of wetdep (= dynamic) timesteps
            DO L = 1, LD39
               ARRAY(:,:,L) = AD39(:,:,L,N) / SCALEDYN
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD39,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD39) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND40: Free diagnostic
!
!  ND41: Afternoon PBL depth (between 1200 and 1600 Local Time)
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) PBLDEPTH : Afternoon PBL heights       : m         : AFTTOT
!
!  NOTES:
!  (1) Bug fix: write one level to binary punch file (bmy, 12/12/00)
!******************************************************************************
!
      IF ( ND41 > 0 ) THEN 
         CATEGORY  = 'PBLDEPTH'
         SCALE_TMP = FLOAT( AFTTOT ) + 1d-20 

         DO M = 1, TMAX(41)
            N = TINDEX(41,M)
            IF ( N > PD41 ) CYCLE
            NN = N
               
            ! Select the proper unit string
            SELECT CASE ( N ) 
               CASE ( 1 )
                  UNIT = 'm'
               CASE ( 2 )
                  UNIT = 'level'
            END SELECT
                     
            ARRAY(:,:,1) = AD41(:,:,N) / SCALE_TMP

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF   
!
!******************************************************************************
!  ND42: Free diagnostic as of 11/24/99
!
!  ND43: Chemical production of OH and NO
!
!   # : Field : Description             : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1)  OH    : OH  Chemical Diagnostic : mol/cm3 : CTOH
!  (2)  NO    : NO  Chemical Diagnostic : v/v     : CTNO
!  (3)  HO2   : HO2 Chemical Diagnostic : v/v     : CTHO2
!  (4)  NO2   : NO2 Chemical Diagnostic : v/v     : CTNO2
!  (5)  NO3   : NO3 Chemical Diagnostic : v/v     : CTNO3
!
!  NOTES:
!  (1) Print output for either a NOx-Ox-HC run (NSRCX == 3), or a CO run
!       with parameterized OH (NSRCX == 5).  (bmy, 4/17/00)
!  (2) Add parentheses in IF test since .AND. has higher precedence
!       than .OR. (jsw, bmy, 12/5/00)
!  (3) Added HO2, NO2 to ND43 (rvm, bmy, 2/27/02)
!  (4) Added NO3 to ND43 (bmy, 1/16/03)
!******************************************************************************
!
      IF ( ND43 > 0 .and. ( NSRCX == 3 .or. NSRCX == 5 ) ) THEN 
         CATEGORY = 'CHEM-L=$' 

         DO M = 1, TMAX(43)
            N  = TINDEX(43,M)
            NN = N 

            SELECT CASE ( N )

               ! OH
               CASE ( 1 )
                  SCALE_TMP(:,:) = FLOAT( CTOH ) + 1d-20
                  UNIT           = 'mol/cm3'
                     
                  !================================================
                  ! For the CO-OH parameterization option, 24-hr 
                  ! chemistry time steps are taken.  As a result, 
                  ! the last day's OH is not calculated, so adjust 
                  ! SCALE_TMP accordingly. (bnd, bmy, 4/18/00)
                  !================================================
                  IF ( NSRCX == 5 ) THEN
                     SCALE_TMP(:,:) = SCALE_TMP(:,:) - 1
                  ENDIF

               ! NO
               CASE ( 2 ) 
                  SCALE_TMP(:,:) = FLOAT( CTNO ) + 1d-20
                  UNIT           = 'v/v'

               ! HO2 (rvm, bmy, 2/27/02)
               CASE ( 3 ) 
                  SCALE_TMP(:,:) = FLOAT( CTHO2 ) + 1d-20
                  UNIT           = 'v/v'

               ! NO2 (rvm, bmy, 2/27/02)
               CASE ( 4 ) 
                  SCALE_TMP(:,:) = FLOAT( CTNO2 ) + 1d-20
                  UNIT           = 'v/v'

               ! NO3 (rjp, bmy, 1/16/03)
               CASE ( 5 )
                  SCALE_TMP(:,:) = FLOAT( CTNO3 ) + 1d-20
                  UNIT           = 'v/v'

               CASE DEFAULT
                  CYCLE
            END SELECT

            DO L = 1, LD43 
               ARRAY(:,:,L) = AD43(:,:,L,N) / SCALE_TMP(:,:)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD43,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD43) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND44: Drydep flux (molec/cm2/s) and velocity (cm/s) diagnostics
!
!   #   : Field    : Quantity           : Units               : Scale factor
!  -------------------------------------------------------------------------
!  (1 ) : DRYD-FLX : drydep fluxes      : molec/cm2/s or kg/s : SCALECHEM
!  (2 ) : DRYD-VEL : drydep velocities  : cm/s                : SCALECHEM
!
!  NOTES: 
!  (1 ) Remove diagnostics for wetdep HNO3, H2O2 from ND44.
!  (2 ) For NSRCX == 1 (Rn-Pb-Be), save the actual tracer number 
!        instead of the dry deposition index.  Add TRCOFFSET to N.
!  (3 ) For NSRCX == 6 (single tracer Ox), drydep fluxes are in kg/s.
!  (4 ) ND44 now uses allocatable array AD44 instead of AIJ. (bmy, 3/16/00)
!  (5 ) Add code from amf for multi-tracer Ox (bmy, 7/3/01)
!  (6 ) Now divide by SCALECHEM since DRYFLX is only called after the
!        chemistry routines for all relevant simulations (bmy, 1/27/03)
!  (7 ) Now print out NTRACE drydep fluxes for tagged Ox.  Also tagged Ox 
!        now saves drydep in molec/cm2/s. (bmy, 8/19/03)
!******************************************************************************
!
      IF ( ND44 > 0 ) THEN
         
         !==============================================================
         ! Drydep fluxes
         !==============================================================
         IF ( NSRCX == 6 ) THEN
            M = NTRACE
         ELSE
            M = NUMDEP
         ENDIF

         DO N = 1, M

            ! Tracer number plus GAMAP offset
            IF ( NSRCX == 6 ) THEN
               NN = N + TRCOFFSET
            ELSE
               NN = NTRAIND(N) + TRCOFFSET
            ENDIF

            CATEGORY = 'DRYD-FLX'

            SELECT CASE ( NSRCX )

               ! Rn-Pb-Be
               CASE ( 1 )
                  UNIT         = 'kg/s'
                  ARRAY(:,:,1) = AD44(:,:,N,1) / SCALECHEM

               ! Default (fullchem)
               CASE DEFAULT
                  UNIT         = 'molec/cm2/s'
                  ARRAY(:,:,1) = AD44(:,:,N,1) / SCALECHEM 
                     
            END SELECT

            ! Write to file
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         ENDDO

         !==============================================================
         ! Drydep velocities
         !==============================================================
         DO N = 1, NUMDEP

            ! Tracer number plus GAMAP offset
            NN           = NTRAIND(N) + TRCOFFSET
            CATEGORY     = 'DRYD-VEL'
            UNIT         = 'cm/s'
            ARRAY(:,:,1) = AD44(:,:,N,2) / SCALESRCE

            ! Write to file
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         ENDDO
      ENDIF
!
!******************************************************************************
!  ND45: Tracer Mixing Ratio (v/v) for Levels L = 1, LD45
!        averaged between hours OTH_HR1 and OTH_HR2
!
!   # : Field   : Description            : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) IJ-AVG-$ : Tracer mixing ratio    : v/v     : CTOTH
!
!  NOTES:
!  (1) For NSRCX == 3 (NOx-Ox-HC run), store pure O3 with index NTRACE+1.
!  (2) Now store pure O3 as NNPAR+1 (now tracer #32). (bmy, 1/10/03)
!******************************************************************************
!
      IF ( ND45 > 0 ) THEN
         CATEGORY  = 'IJ-AVG-$'
         SCALE_TMP = FLOAT( CTOTH ) + 1d-20
         UNIT = ''   

         DO M = 1, TMAX(45)
            N  = TINDEX(45,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET

            DO L = 1, LD45
               ARRAY(:,:,L) = AD45(:,:,L,N) / SCALE_TMP(:,:)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD45,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD45) )

            ! Store pure O3 as NNPAR+1 (bmy, 1/10/03)
            IF ( NSRCX == 3 .and. NN == IDTOX ) THEN 
               DO L = 1, LD45
                  ARRAY(:,:,L) =
     &                 AD45(:,:,L,NTRACE+1) / SCALE_TMP(:,:)
               ENDDO
                  
               NN = NNPAR + 1 + TRCOFFSET   

               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     LD45,     IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD45) )
            ENDIF
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND46: Biogenic source diagnostic
!
!   # : Field : Description  : Units         : Scale Factor
!  ---------------------------------------------------------------------------
!  (1)  ISOP  : Isoprene     : atoms C/cm2/s : SCALE3
!  (2)  ACET  : Acetone      : atoms C/cm2/s : SCALE3
!  (3)  PRPE  : Propene      : atoms C/cm2/s : SCALE3
!  (4)  MONOT : Monoterpenes : atoms C/cm2/s : SCALE3
!  
!  NOTES:
!  (1) ND46 now uses allocatable array AD46 instead of AIJ (bmy, 3/16/00)
!  (2) Also write out PRPE for CO-OH run (NSRCX == 5), regardless of
!       the setting of IDTPRPE.  This is to print out monterpene 
!       diagnostics. (bnd, bmy, 4/18/00)
!  (3) Added monoterpenes as tracer #4.  This requires updated versions
!       of "tracerinfo.dat" and "diaginfo.dat" for GAMAP. (bmy, 1/2/01)
!******************************************************************************
!
      IF ( ND46 > 0 ) THEN
         CATEGORY = 'BIOGSRCE'
         UNIT     = ''

         DO M = 1, TMAX(46)
            N  = TINDEX(46,M)
            IF ( N > PD46 ) CYCLE
            NN = N

            ! Skip if ISOP, ACET, PRPE are not tracers
            IF ( N == 1 .and. IDTISOP == 0 ) CYCLE
            IF ( N == 2 .and. IDTACET == 0 ) CYCLE
            IF ( N == 3 .and. IDTPRPE == 0 ) CYCLE
            
            ARRAY(:,:,1) = AD46(:,:,N) / SCALESRCE
               
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND47: Tracer Mixing Ratio (v/v) for Levels L = 1, LD47
!        *always* averaged between 0000 and 2400 Local time.
!
!   # : Field   : Description                 : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) IJ-24H-$ : 24h avg Tracer mixing ratio : v/v     : SCALEDYN
!
!  NOTES:
!  (1) For NSRCX == 3 (NOx-Ox-HC run), store pure O3 with index NTRACE+1.
!  (2) Now store pure O3 as NNPAR+1 (now tracer #32). (bmy, 1/10/03) 
!  (3) Now replace SCALE1 with SCALEDYN
!******************************************************************************
!
      IF ( ND47 > 0 ) THEN
         CATEGORY = 'IJ-24H-$'
         UNIT     = ''   

         DO M = 1, TMAX(47)
            N  = TINDEX(47,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET
            
            DO L = 1, LD47
               ARRAY(:,:,L) = AD47(:,:,L,N) / SCALEDYN
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD47,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD47) )

            ! Store pure O3 as NNPAR+1 (bmy, 1/10/03)
            IF ( NSRCX == 3 .and. NN == IDTOX ) THEN 
               DO L = 1, LD47
                  ARRAY(:,:,L) = AD47(:,:,L,NTRACE+1) / SCALEDYN
               ENDDO
                     
               NN = NNPAR + 1 + TRCOFFSET   

               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     LD47,     IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD47) )
            ENDIF
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND55: Tropopause diagnostics
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) TP-LEVEL : Tropopause level            : unitless  : SCALEDYN
!  (2) TP-HGHT  : Tropopause height           : km        : SCALEDYN
!  (3) TP-PRESS : Tropopause pressure         : mb        : SCALEDYN
!******************************************************************************
!
      IF ( ND55 > 0 ) THEN
         CATEGORY = 'TR-PAUSE'

         DO M = 1, TMAX(55)
            N  = TINDEX(55,M)
            IF ( N > PD55 ) CYCLE
            NN = N

            ! Pick the appropriate unit string
            SELECT CASE ( N )
               CASE ( 1 )
                  UNIT = 'unitless'
               CASE ( 2 )
                  UNIT = 'km'
               CASE ( 3 )
                  UNIT = 'mb'
            END SELECT

            ARRAY(:,:,1) = AD55(:,:,N) / SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!*****************************************************************************
!  ND62: I-J Instantaneous Column Maps for Tracers (molec/cm^2)  
!
!  The unit conversion is as follows:
!
!    STT (kg) | 6.022e23 molec |   mole   | 1000 g |    1        |   m^2
!    ---------+----------------+----------+--------+-------------+----------
!             |     mole       | TCMASS g |  kg    | AREA_M2 m^2 | 10^4 cm^2
!
!
!  which is equivalent to
!
!   ( STT * 6.022e22 ) / ( TCMASS * AREA_M2 )
!*****************************************************************************
!
      IF ( ND62 > 0 ) THEN
         CATEGORY = 'INST-MAP'

         DO M = 1, TMAX(62)
            N  = TINDEX(62,M)
            IF ( N > NTRACE ) CYCLE
            NN = N + TRCOFFSET

            DO J = 1, JJPAR

               ! Grid box surface area [cm2]
               AREA_M2 = GET_AREA_M2( J )

               DO I = 1, IIPAR
                  ARRAY(I,J,1) = ( SUM( STT(I,J,:,N) ) * 6.022d22 )
     &                         / ( TCMASS(N)           * AREA_M2  )
               ENDDO
            ENDDO

            ! Write the proper unit string
            IF ( TCMASS(N) > 12d0 ) THEN
               UNIT = 'molec/cm2'
            ELSE
               UNIT = 'atoms C/cm2'
            ENDIF

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         ENDDO
      ENDIF
!
!******************************************************************************
!  ND65: Production/Loss of specified chemical families
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) PORL-L=$ : Chemical family P-L rates   : mol/cm3/s : SCALECHEM
!
!  NOTES:
!  (1 ) Make sure the units for NSRCX == 6 (single tracer O3) P-L 
!        coincide with those in "chemo3.f".  
!  (2 ) ND65 now uses allocatable array AD65 instead of AIJ. (bmy, 3/16/00)
!  (3 ) Add L(CH3I) to the ND65 diagnostic -- do not take the average 
!        but instead compute the total sum of L(CH3I) (nad, bmy, 3/20/01)
!  (4 ) Add updates for multi-tracer Ox run from amf (bmy, 7/3/01)
!******************************************************************************
!
      IF ( ND65 > 0 ) THEN     
         CATEGORY = 'PORL-L=$'

         DO M = 1, TMAX(65)
            N  = TINDEX(65,M)

            ! Cycle if we exceed the number of families
            IF ( N > NFAM ) CYCLE

            ! Don't add TRCOFFSET for single tracer Ox
            ! Also select proper unit string
            SELECT CASE ( NSRCX )
               CASE ( 2 )
                  NN     = N + TRCOFFSET
                  UNIT   = 'kg/s'
                  SCALEX = 1d0

               CASE ( 6 )
                  NN     = N + TRCOFFSET
                  UNIT   = 'kg/s'
                  SCALEX = SCALECHEM
                  
               CASE ( 10 )
                  NN     = N + TRCOFFSET
                  UNIT   = 'mol/cm3/s'
                  SCALEX = SCALECHEM
                  
               CASE DEFAULT 
                  NN     = N + TRCOFFSET
                  UNIT   = 'mol/cm3/s'
                  SCALEX = SCALECHEM
                  
            END SELECT

            DO L = 1, LD65
               ARRAY(:,:,L) = AD65(:,:,L,N) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD65,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD65) )
         ENDDO       
      ENDIF
!
!******************************************************************************
!  ND66: GMAO 3-D fields 
!
!   # : Field  : Description                       : Units   : Scale factor
!  --------------------------------------------------------------------------
!  (1)  UWND   : GMAO Zonal Winds                  : m/s     : SCALE_I6 or _A6
!  (2)  VWND   : GMAO Meridional Winds             : m/s     : SCALE_I6 or _A6
!  (3)  TMPU   : GMAO Temperatures                 : K       : SCALE_I6 or _A6
!  (4)  SPHU   : GMAO Specific Humidity            : g/kg    : SCALE_I6 or _A6
!  (5)  CLDMAS : GMAO Cloud Mass Flux              : kg/m2/s : SCALE_A6 or _A6
!
!  NOTES:
!  (1) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2) Add CLDMAS to ND66 diagnostic as field #6, but with tracer index
!       #7 (for compatibility with the existing GAMAP).  (rvm, bmy, 9/8/00)
!  (3) For GEOS-4/fvDAS, UWND, VWND, TMPU, SPHU are A-6 fields.  Adjust
!       the scale factors accordingly.  Also delete KZZ. (bmy, 6/23/03)
!******************************************************************************
!
      IF ( ND66 > 0 ) THEN
         CATEGORY = 'DAO-3D-$'

         DO M = 1, TMAX(66)
            N  = TINDEX(66,M)
            NN = N 
            
            SELECT CASE ( N )

               ! UWND, VWND
               CASE ( 1,2 )
#if   defined( GEOS_4 )
                  SCALEX = SCALE_A6
#else
                  SCALEX = SCALE_I6
#endif

                  UNIT   = 'm/s'

               ! TMPU
               CASE ( 3 )
#if   defined( GEOS_4 )
                  SCALEX = SCALE_A6
#else
                  SCALEX = SCALE_I6
#endif
                  UNIT   = 'K'

               ! SPHU
               CASE ( 4 )
#if   defined( GEOS_4 )
                  SCALEX = SCALE_A6
#else
                  SCALEX = SCALE_I6
#endif
                  UNIT   = 'g/kg'

               ! CLDMAS
               CASE( 5 )
                  NN     = 7
                  SCALEX = SCALE_A6
                  UNIT   = 'kg/m2/s'

               CASE DEFAULT
                  CYCLE
            END SELECT

            ARRAY(:,:,1:LD66) = AD66(:,:,1:LD66,N) / SCALEX

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD66,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD66) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND67: GMAO surface fields
!
!   # : Field  : Description                      : Units    : Scale factor
!  -----------------------------------------------------------------------
!  (1 ) HFLUX  : GMAO Sensible Heat Flux          : W/m2     : SCALE_A3
!  (2 ) RADSWG : GMAO Insolation @ Surface        : W/m2     : SCALE_A3
!  (3 ) PREACC : GMAO Accum Precip @ Surface      : mm/day   : SCALE_A3
!  (4 ) PRECON : GMAO Conv Precip @ Surface       : mm/day   : SCALE_A3
!  (5 ) TS     : GMAO Surface Air Temperature     : K        : SCALE_A3
!  (6 ) RADSWT : GMAO Insolation @ Top of Atm     : W/m2     : SCALE_A3
!  (7 ) USTAR  : GMAO Friction Velocity           : m/s      : SCALE_A3
!  (8 ) Z0     : GMAO Roughness Height            : m        : SCALE_A3
!  (9 ) PBL    : GMAO PBL depth                   : mb       : SCALE_A3
!  (10) CLDFRC : GMAO Cloud Fraction              : unitless : SCALE_A3
!  (11) U10M   : GMAO U-wind @ 10 m               : m/s      : SCALE_A3
!  (12) V10M   : GMAO V-wind @ 10 m               : m/s      : SCALE_A3
!  (13) PS-PBL : GMAO Boundary Layer Top Pressure : mb       : SCALEDYN
!  (14) ALBD   : GMAO Surface Albedo              : unitless : SCALE_I6 
!  (15) PHIS   : GMAO Geopotential Heights        : m        : SCALED 
!  (16) CLTOP  : GMAO Cloud Top Height            : levels   : SCALE_A6
!  (17) TROPP  : GMAO Tropopause pressure         : mb       : SCALE_I6
!  (18) SLP    : GMAO Sea Level pressure          : mb       : SCALE_I6
!  (19) TSKIN  : Ground/sea surface temp.         : hPa      : SCALE_A3
!  (20) PARDF  : Photosyn active diffuse rad.     : W/m2     : SCALE_A3
!  (21) PARDR  : Photosyn active direct  rad.     : W/m2     : SCALE_A3
!  (22) GWET   : Top soil wetness                 : unitless : SCALE_A3
!
!  NOTES:
!  (1 ) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2 ) Now use AD67 allocatable array (bmy, 2/17/00)
!  (3 ) Add TROPP as tracer #17 and SLP as tracer #18 (bmy, 10/11/00)
!  (4 ) Now replace SCALE1 with SCALEDYN (bmy, 3/27/03)
!  (5 ) Added TSKIN, PARDF, PARDR, GWET for GEOS-4 (bmy, 6/23/03)
!******************************************************************************
!
      IF ( ND67 > 0 ) THEN
         CATEGORY = 'DAO-FLDS'

         ! Binary punch file
         DO M = 1, TMAX(67)
            N  = TINDEX(67,M)
            NN = N 

            SELECT CASE ( N ) 
               CASE ( 1, 2, 6 )
                  SCALEX = SCALE_A3
                  UNIT   = 'W/m2'
               CASE ( 3, 4 )
                  SCALEX = SCALE_A3
                  UNIT   = 'mm/day'
               CASE ( 5 )
                  SCALEX = SCALE_A3
                  UNIT   = 'K'
               CASE ( 7, 11, 12 )
                  SCALEX = SCALE_A3
                  UNIT   = 'm/s'
               CASE ( 8 )
                  SCALEX = SCALE_A3
                  UNIT   = 'm'
               CASE ( 9 )
                  SCALEX = SCALE_A3
                  UNIT   = 'hPa'
               CASE ( 10 )
                  SCALEX = SCALE_A3
                  UNIT   = 'unitless'
               CASE ( 13 )
                  SCALEX = SCALEDYN
                  UNIT   = 'hPa'
               CASE ( 14 ) 
                  SCALEX = SCALE_I6 
                  UNIT   = 'unitless'
               CASE ( 15 )
                  SCALEX = SCALED
                  UNIT   = 'm'
               CASE ( 16 ) 
                  SCALEX = SCALE_A6
                  UNIT   = 'levels'
               CASE ( 17 ) 
                  SCALEX = SCALE_I6
                  UNIT   = 'hPa'
                  NN     = 19   ! for GAMAP
               CASE ( 18 ) 
                  SCALEX = SCALE_I6
                  UNIT   = 'hPa'
                  NN     = 21   ! for GAMAP
               CASE ( 19 ) 
                  SCALEX = SCALE_A3
                  UNIT   = 'K'
                  NN     = 29   ! for GAMAP
               CASE ( 20 ) 
                  SCALEX = SCALE_A3
                  UNIT   = 'W/m2'
                  NN     = 32   ! for GAMAP
               CASE ( 21 ) 
                  SCALEX = SCALE_A3
                  UNIT   = 'W/m2'
                  NN     = 33   ! for GAMAP
               CASE ( 22 ) 
                  SCALEX = SCALE_A3
                  UNIT   = 'unitless'
                  NN     = 27   ! for GAMAP
               CASE DEFAULT
                  CYCLE
            END SELECT
                  
            ARRAY(:,:,1) = AD67(:,:,N) / SCALEX

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND68: Grid box diagnostics
!
!   # : Field   : Description                       : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1) BXHEIGHT : Grid box height                   : m     : SCALEDYN
!  (2) AD       : Air mass in grid box              : kg    : SCALEDYN
!  (3) AVGW     : Mixing ratio of water vapor       : v/v   : SCALEDYN
!  (4) N(AIR)   : Number density of air             : m^-3  : SCALEDYN
!
!  NOTES:
!  (1) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2) Now replaced SCALE1 with SCALEDYN (bmy, 2/24/03)
!******************************************************************************
!
      IF ( ND68 > 0 ) THEN
         CATEGORY = 'BXHGHT-$'
         UNIT     = ''

         DO M = 1, TMAX(68)
            N  = TINDEX(68,M)
            IF ( N > PD68 ) CYCLE
            NN = N 

            DO L = 1, LD68
               ARRAY(:,:,L) = AD68(:,:,L,N) / SCALEDYN
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     ND68,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:ND68) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND69: Grid Box Surface Areas
!
!   # : Field : Description                       : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1) DXYP   : Surface area of grid box          : m^2   : SCALED = 1.0
!
!  NOTES:
!  (1) Only print DXYP for the first timestep, as it is an invariant field.
!  (2) We don't need to add TRCOFFSET to N.  This is not a CTM tracer.
!  (3) Now use the AD69 dynamically allocatable array (bmy, 2/17/00)
!******************************************************************************
!
      IF ( ND69 > 0 ) THEN 
         CATEGORY = 'DXYP'
         UNIT     = 'm2'

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 1,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    AD69(:,:,1) )

         ! Set ND69 = 0 so we won't print it out again
         ND69 = 0
      ENDIF

      ! Echo output
      WRITE( 6, '(a)' ) '     - DIAG3: Diagnostics written to bpch!'

      ! Return to calling program
      END SUBROUTINE DIAG3    
