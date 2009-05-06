 ! $Id: define.h,v 1.68 2009/05/06 14:45:12 ccarouge Exp $
!
!******************************************************************************
!  Include file "define.h" specifies C-preprocessor "switches" that are 
!  used to include or exclude certain sections of code.  
!  (bmy, bdf, 1/30/98, 11/6/08)
!
!  List of "Switches"
!  ===========================================================================
!  (1 ) GCAP        : Enables code for GCAP   met fields & chemistry
!  (2 ) GEOS_3      : Enables code for GEOS-3 met fields & chemistry
!  (3 ) GEOS_4      : Enables code for GEOS-4 met fields & chemistry
!  (4 ) GEOS_5      : Enables code for GEOS-5 met fields & chemistry
!  (5 ) GRIDREDUCED : Enables code for reduced stratosphere grids
!  (6 ) GRID1x1     : Enables code for 1 x 1    GLOBAL        GRID
!  (7 ) NESTED_CH   : Enables code for 1 x 1    CHINA  NESTED GRID
!  (8 ) NESTED_NA   : Enables code for 1 x 1    N. AM. NESTED GRID
!  (9 ) GRID1x125   : Enables code for 1 x 1.25 GLOBAL        GRID
!  (10) GRID2x25    : Enables code for 2 x 2.5  GLOBAL        GRID
!  (11) GRID4x5     : Enables code for 4 x 5    GLOBAL        GRID 
!  (12) COMPAQ      : Enables code for Alpha w/ COMPAQ/HP Alpha compiler
!  (13) IBM_AIX     : Enables code for IBM/AIX compiler
!  (14) LINUX_PGI   : Enables code for Linux w/ PGI compiler
!  (15) LINUX_IFORT : Enables code for Linux v8 or v9 "IFORT" compiler
!  (16) SGI_MIPS    : Enables code for SGI Origin w/ MIPS compiler
!  (17) SPARC       : Enables code for Sun w/ SPARC or Sun Studio compiler
! 
!  NOTES:
!  (1 ) "define.h" is #include'd at the top of CMN_SIZE.  All subroutines
!        that normally reference CMN_SIZE will also reference "define.h".
!  (2 ) Only define the "switches" that are *absolutely* needed for a
!        given implementation, as the criteria for code inclusion/exclusion 
!        is the #if defined() statement.  Undefined "switches" are "off".
!  (3 ) To turn off a switch, comment that line of code out.
!  (4 ) As of 11/30/99, DO_MASSFLUX is obsolete, since the mass flux
!        arrays are now declared allocatable in "diag_mod.f".  
!  (5 ) Eliminate DO_MASSB switch -- ND63 diagnostic is now obsolete.
!        (bmy, 4/12/00)
!  (6 ) Add GEOS_3 and GRID1x1 switches for future use (bmy, 7/7/00)
!  (7 ) Make sure that one of FULLCHEM, SMALLCHEM, or LGEOSCO is turned on.
!        Also cosmetic changes. (bmy, 10/3/00)
!  (8 ) Added new switches "DEC_COMPAQ" and "SGI" (bmy, 3/9/01) 
!  (9 ) Added new "LINUX" switch (bmy, 7/16/01)
!  (10) Added new "GEOS_4" switch for GEOS-4/fvDAS met fields (bmy, 11/21/01)
!  (11) Now enclose switch names in ' ', since the PGI compiler chokes 
!        on barewords (bmy, 3/20/02)
!  (12) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!        compilation (bmy, 6/25/02)
!  (13) Removed GEOS_2 switch; added GEOS_4 switch.  Also added SPARC switch 
!        to invoke Sun/Sparc specific code. (bmy, 3/23/03)
!  (14) Added IBM_AIX switch (bmy, 6/27/03)
!  (15) Added INTEL_FC switch (bmy, 10/21/03)
!  (16) Added GRID30LEV switch for 30L GEOS-3 or GEOS-4 grid (bmy, 10/31/03)
!  (17) Renamed cpp switch "LINUX" to "LINUX_PGI".  Renamed cpp switch 
!        "INTEL_FC" to "LINUX_IFC".  Renamed cpp switch "SGI" to "SGI_MIPS".
!        Added cpp switch "LINUX_EFC".  Removed cpp switch SMALLCHEM.
!        (bmy, 12/2/03)
!  (18) Added "A_LLK_03" switch to denote GEOS-4 "a_llk_03" met fields.  This
!        will be temporary since "a_llk_03" met fields will be replaced by
!        a newer product.  (bmy, 3/22/04) 
!  (19) Added NESTED_NA and NESTED_CH cpp switches.  Also add GRID1x125
!        cpp switch. (bmy, 12/1/04)
!  (20) Removed obsolete A_LLK_03, LFASTJ, LSLOWJ, FULLCHEM, LGEOSCO switches.
!        Also added extra switches for GCAP and GEOS_5 met fields.  
!        (bmy, 6/23/05)
!  (21) Added LINUX_IFORT switch to delineate Intel compilers v8 or v9 
!        from v7. (bmy, 10/18/05)
!  (22) Removed GEOS_1, GEOS_STRAT, LINUX_IFC, LINUX_EFC (bmy, 8/4/06)
!  (23) Renamed GRID30LEV to GRIDREDUCED (bmy, 2/7/07)
!  (24) Added IN_CLOUD_OD flag for reprocessed GEOS-5 met.  Added GRID05x0666 
!        flag for GEOS-5 nested grids (yxw, dan, bmy, hyl, 11/6/08)
!******************************************************************************
!
!==============================================================================
! Undefine all "switches" so that they cannot be accidentally reset  
!==============================================================================
#undef GCAP
#undef GEOS_3
#undef GEOS_4
#undef GEOS_5
#undef GRIDREDUCED
#undef GRID4x5
#undef GRID2x25  
#undef GRID1x125
#undef GRID1x1
#undef GRID05x0666
#undef NESTED_NA
#undef NESTED_CH
#undef COMPAQ
#undef IBM_AIX
#undef LINUX_PGI
#undef LINUX_IFORT
#undef SGI_MIPS
#undef SPARC
#undef IN_CLOUD_OD

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

!----- Grid sizes -----
!#define NESTED_CH   'NESTED_CH'
!#define NESTED_NA   'NESTED_NA'
!#define GRID05x0666 'GRID05x0666'
!#define GRID1x1     'GRID1x1'
!#define GRID1x125   'GRID1x125'
!#define GRID2x25    'GRID2x25'
#define GRID4x5     'GRID4x5'
#define GRIDREDUCED 'GRIDREDUCED'

!----- Compilers -----
!#define COMPAQ      'COMPAQ'
!#define IBM_AIX     'IBM_AIX'
!#define LINUX_PGI   'LINUX_PGI'
#define LINUX_IFORT 'LINUX_IFORT'
!#define SGI_MIPS    'SGI_MIPS'
!#define SPARC       'SPARC'

!----- FOR GEOS-5 MET FIELDS ONLY -----
! NOTE: If you are using GEOS-5 met fields that were reprocessed to 
! correctly regrid the in-cloud optical depth and cloud fraction fields, 
! then be sure to uncomment the following line of code.  This will cause 
! FAST-J to interpret the optical depth correctly.  Leaving this option
! commented will cause a "quick fix" (i.e. multiplying the optical depth
! by the cloud fracton) to be applied, which should be a good enough fix
! in the meantime. (bmy, hyl, 10/24/08)
#define IN_CLOUD_OD  'IN_CLOUD_OD'

!==============================================================================
! Force a compile error if IN_CLOUD_OD is used with GEOS_3 or GEOS_4  
!==============================================================================
#if defined(GEOS_3) || defined(GEOS_4) || defined (GCAP)
#if defined(IN_CLOUD_OD)
#error "ERROR: IN_CLOUD_OD option set with GEOS_3, GEOS_4, or GCAP"
#endif
#endif

!==============================================================================
! Force a compile error if GEOS_1, GEOS_STRAT, GEOS_3, GEOS_4 are undefined 
!==============================================================================
#if !defined(GEOS_3) && !defined(GEOS_4) && !defined(GEOS_5) && !defined(GCAP)
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
#if !defined(COMPAQ) && !defined(IBM_AIX) && !defined(LINUX_PGI) && !defined(LINUX_IFORT) && !defined(SGI_MIPS) && !defined(SPARC)
#error "ERROR: One of COMPAQ, IBM_AIX, LINUX_PGI, LINUX_IFORT,"
#error "SGI_MIPS, SPARC must be defined in header file define.h"
#endif
