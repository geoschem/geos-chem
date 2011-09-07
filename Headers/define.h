!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: define.h
!
! !DESCRIPTION: Include file "define.h" specifies C-preprocessor "switches" 
!  that are used to include or exclude certain sections of code.  
!\\
!\\
! !REMARKS:
!  List of "Switches"
!  ===========================================================================
!  (1 ) GCAP        : Enables code for GCAP   met fields & chemistry
!  (2 ) GEOS_3      : Enables code for GEOS-3 met fields & chemistry
!  (3 ) GEOS_4      : Enables code for GEOS-4 met fields & chemistry
!  (4 ) GEOS_5      : Enables code for GEOS-5 met fields & chemistry
!  (5 ) MERRA       : Enables code for MERRA  met fields & chemistry
!  (6 ) GRIDREDUCED : Enables code for reduced stratosphere grids
!  (7 ) GRID1x1     : Enables code for 1 x 1    GLOBAL        GRID
!  (8 ) NESTED_CH   : Enables code for CHINA  NESTED GRID
!  (9 ) NESTED_NA   : Enables code for N. AM. NESTED GRID
!  (10) NESTED_EUR  : Enables code for EUROPE NESTED GRID
!  (11) GRID1x125   : Enables code for 1 x 1.25 GLOBAL        GRID
!  (12) GRID2x25    : Enables code for 2 x 2.5  GLOBAL        GRID
!  (13) GRID4x5     : Enables code for 4 x 5    GLOBAL        GRID 
!  (14) IBM_AIX     : Enables code for IBM/AIX compiler
!  (15) IBM_XLF     : Enables code for IBM/XLF compiler
!  (16) LINUX_PGI   : Enables code for Linux w/ PGI compiler
!  (17) LINUX_IFORT : Enables code for Linux v8 or v9 "IFORT" compiler
!  (18) SPARC       : Enables code for Sun w/ SPARC or Sun Studio compiler
!  (19) GTMM_Hg     : Enables code for Hg simulation with GTMM
!                                                                            .
!  NOTES:
!  (1 ) "define.h" is #include'd at the top of CMN_SIZE.  All subroutines
!        that normally reference CMN_SIZE will also reference "define.h".
!  (2 ) Only define the "switches" that are *absolutely* needed for a
!        given implementation, as the criteria for code inclusion/exclusion 
!        is the #if defined() statement.  Undefined "switches" are "off".
!  (3 ) To turn off a switch, comment that line of code out.
!
! !REVISION HISTORY:
!  30 Nov 1999 - R. Yantosca - DO_MASSFLUX is obsolete, since the mass flux
!                              arrays are now declared allocatable in 
!                             "diag_mod.f".  
!  12 Apr 2000 - R. Yantosca - Eliminate DO_MASSB switch -- ND63 diagnostic 
!                              is now obsolete.
!  07 Jul 2000 - R. Yantosca - Add GEOS_3 and GRID1x1 switches for future use
!  03 Oct 2000 - R. Yantosca - Make sure that one of FULLCHEM, SMALLCHEM, or 
!                              LGEOSCO is turned on.  Also cosmetic changes. 
!  03 Sep 2001 - R. Yantosca - Added new switches "DEC_COMPAQ" and "SGI" 
!  16 Jul 2001 - R. Yantosca - Added new "LINUX" switch\
!  21 Nov 2001 - R. Yantosca - Added new "GEOS_4" switch for GEOS-4/fvDAS 
!                              met fields
!  20 Mar 2002 - R. Yantosca - Now enclose switch names in ' ', since the 
!                              PGI compiler chokes on barewords
!  25 Jun 2002 - R. Yantosca - Changed RCS ID tag comment character from "C" 
!                              to "!" to allow freeform compilation
!  23 Mar 2003 - R. Yantosca - Removed GEOS_2 switch; added GEOS_4 switch.  
!                              Also added SPARC switch to invoke Sun/Sparc 
!                              specific code.
!  27 Mar 2003 - R. Yantosca - Added IBM_AIX switch
!  21 Oct 2003 - R. Yantosca - Added INTEL_FC switch
!  31 Oct 2003 - R. Yantosca - GRID30LEV switch for 30L GEOS-3 or GEOS-4 grid
!  02 Dec 2003 - R. Yantosca - Renamed cpp switch "LINUX" to "LINUX_PGI".  
!                              Renamed cpp switch "INTEL_FC" to "LINUX_IFC".  
!                              Renamed cpp switch "SGI" to "SGI_MIPS".
!                              Added cpp switch "LINUX_EFC".  
!                              Removed cpp switch SMALLCHEM.
!  22 Mar 2004 - R. Yantosca - Added "A_LLK_03" switch to denote GEOS-4 
!                              "a_llk_03" met fields.  This will be temporary 
!                              since "a_llk_03" met fields will be replaced by
!                              a newer product.
!  01 Dec 2004 - R. Yantosca - Added NESTED_NA and NESTED_CH cpp switches.  
!                              Also add GRID1x125 cpp switch.
!  23 Jun 2005 - R. Yantosca - Removed obsolete A_LLK_03, LFASTJ, LSLOWJ, 
!                              FULLCHEM, LGEOSCO switches.  Also added extra 
!                              switches for GCAP and GEOS_5 met fields.  
!  18 Oct 2005 - R. Yantosca - Added LINUX_IFORT switch to delineate Intel 
!                              compilers v8 or v9 from v7.
!  04 Aug 2006 - R. Yantosca - Removed obsolete GEOS_1, GEOS_STRAT, LINUX_IFC, 
!                              LINUX_EFC switches.
!  07 Feb 2007 - R. Yantosca - Renamed GRID30LEV to GRIDREDUCED
!  06 Nov 2008 - R. Yantosca - Added IN_CLOUD_OD flag for reprocessed GEOS-5 
!                              met.  Added GRID05x0666 flag for GEOS-5 nested 
!                              grids (cf. yxw, dan, bmy, hyl)
!  08 Jul 2009 - R. Yantosca - Deleted support for old COMPAQ and SGI_MIPS 
!                              compilers.  Added switch for IBM XLF compiler. 
!  15 Oct 2009 - R. Yantosca - Remove IN_CLOUD_OD.  Added ProTex headers.
!  18 Dec 2009 - Aaron van D - Added NESTED_EU C-preprocessor switch
!  20 Jul 2010 - C. Carouge  - Added GTMM_Hg for mercury simulation.
!  12 Aug 2010 - R. Yantosca - Added MERRA switch for MERRA reanalysis met
!EOP
!------------------------------------------------------------------------------
!BOC
!
!==============================================================================
! Undefine all "switches" so that they cannot be accidentally reset  
!==============================================================================
#undef GCAP
#undef GEOS_3
#undef GEOS_4
#undef GEOS_5
#undef MERRA
#undef GRIDREDUCED
#undef GRID4x5
#undef GRID2x25  
#undef GRID1x125
#undef GRID1x1
#undef GRID05x0666
#undef NESTED_NA
#undef NESTED_CH
#undef NESTED_EU
#undef IBM_AIX
#undef IBM_XLF
#undef LINUX_PGI
#undef LINUX_IFORT
#undef SPARC
#undef GTMM_Hg

!==============================================================================
! Define the necessary "switches" for GEOS-CHEM. 
! Give each switch its own name as a value, since this will prevent
! the C-preprocessor overwriting the name everywhere in the code.
!==============================================================================

!----- Model types -----
!#define GCAP        'GCAP'
!#define GEOS_3      'GEOS_3'
!#define GEOS_4      'GEOS_4'
#define GEOS_5      'GEOS_5'
!#define MERRA       'MERRA'

!----- Grid sizes -----
!#define NESTED_CH   'NESTED_CH'
!#define NESTED_NA   'NESTED_NA'
!#define NESTED_EU   'NESTED_EU'
!#define GRID05x0666 'GRID05x0666'
!#define GRID1x1     'GRID1x1'
!#define GRID1x125   'GRID1x125'
!#define GRID2x25    'GRID2x25'
#define GRID4x5     'GRID4x5'
#define GRIDREDUCED 'GRIDREDUCED'

!----- Compilers -----
!#define IBM_AIX     'IBM_AIX'
!#define IBM_XLF     'IBM_XLF'
!#define LINUX_PGI   'LINUX_PGI'
#define LINUX_IFORT 'LINUX_IFORT'
!#define SPARC       'SPARC'

!----- Simulation type -----
!#define GTMM_Hg

!==============================================================================
! Force a compile error if GEOS_1, GEOS_STRAT, GEOS_3, GEOS_4 are undefined 
!==============================================================================
#if !defined(GEOS_3) && !defined(GEOS_4) && !defined(GEOS_5) && !defined(MERRA) && !defined(GCAP)
#error "ERROR: GEOS_STRAT, GEOS_3, GEOS_4, GEOS_5, and GCAP"
#error "are ALL und efined in header file define.h"
#endif

!==============================================================================
! Force a compile error if GRID1x1, GRID2x25, and GRID4x5 are all undefined 
!==============================================================================
#if !defined(GRID2x25) && !defined(GRID4x5) && !defined(GRID1x125) && !defined(GRID1x1) && !defined(GRID05x0666)
#error "ERROR: GRID4x5, GRID2x25, GRID1x125, GRID05x0666 and GRID1x1"
#error "are ALL undefined in header file define.h"
#endif

!==============================================================================
! Force a compile  error if all compiler switches are undefined
!==============================================================================
#if !defined(IBM_AIX) && !defined(IBM_XLF) && !defined(LINUX_PGI) && !defined(LINUX_IFORT) && !defined(SPARC)
#error "ERROR: One of IBM_AIX, IBL_XLF, LINUX_PGI, LINUX_IFORT,"
#error "SPARC must be defined in header file define.h"
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% SCHEM PATCH: Stop the run if we are running w/ GEOS-5 72-level grid
!%%% This prevents dimension mismatch when reading SCHEM data fields!
!%%% (ltm, bmy, 6/2/10)
#if defined( GEOS_5 ) && !defined( GRIDREDUCED )
#error "Cannot run GEOS-5 with the full vertical 72 level grid!"
#error "We are working on a patch to fix this soon!"
#endif 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% SCHEM PATCH: Stop the run if we are running w/ MERRA 72-level grid
!%%% This prevents dimension mismatch when reading SCHEM data fields!
!%%% (ltm, bmy, 6/2/10)
#if defined( MERRA ) && !defined( GRIDREDUCED )
#error "Cannot run MERRA with the full vertical 72 level grid!"
#error "We are working on a patch to fix this soon!"
#endif 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!EOC
