! $Id: define.h,v 1.24 2004/12/17 15:00:05 bmy Exp $
!
!******************************************************************************
!  Include file "define.h" specifies C-preprocessor "switches" that are 
!  used to include or exclude certain sections of code.  
!  (bmy, bdf, 1/30/98, 12/1/04)
!
!  List of "Switches"
!  ===========================================================================
!  (1 ) GEOS_1     : Enables code for GEOS-1 met fields & chemistry
!  (2 ) GEOS_STRAT : Enables code for GEOS-STRAT met fields & chemistry
!  (3 ) GEOS_3     : Enables code for GEOS-3 met fields & chemistry
!  (4 ) GEOS_4     : Enables code for GEOS-4 met fields & chemistry
!  (5 ) A_LLK_03   : Enables code for GEOS-4 Version 3 data ("a_llk_03")
!  (6 ) GRID30LEV  : Enables code for 30-level GEOS-3 or GEOS-4 grid
!  (7 ) GRID1x1    : Enables code for 1 x 1   GLOBAL        GRID
!  (8 ) NESTED_CH  : Enables code for 1 x 1   CHINA  NESTED GRID
!  (9 ) NESTED_NA  : Enables code for 1 x 1   N. AM. NESTED GRID
!  (10) GRID2x25   : Enables code for 2 x 2.5 GLOBAL        GRID
!  (11) GRID4x5    : Enables code for 4 x 5   GLOBAL        GRID 
!  (12) FULLCHEM   : Enables code for "Full" Chemistry (ISOP and NMHC)
!  (13) LGEOSCO    : Enables code for CO run w/ parameterized OH
!  (14) LFASTJ     : Enables code for FAST-J photolysis
!  (15) LSLOWJ     : Enables code for SLOW-J photolysis
!  (16) COMPAQ     : Enables code for Alpha w/ COMPAQ/HP Alpha compiler
!  (17) IBM_AIX    : Enables code for IBM/AIX compiler
!  (18) LINUX_PGI  : Enables code for Linux w/ PGI compiler
!  (19) LINUX_IFC  : Enables code for Linux w/ 32-bit Intel Fortran compiler
!  (20) LINUX_EFC  : Enables code for Linux w/ 64-bit Intel Fortran compiler
!  (21) SGI_MIPS   : Enables code for SGI Origin w/ MIPS compiler
!  (22) SPARC      : Enables code for Sun w/ SPARC compiler
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
!******************************************************************************
!
!==============================================================================
! Undefine all "switches" so that they cannot be accidentally reset  
!==============================================================================
#undef GEOS_1
#undef GEOS_STRAT
#undef GEOS_3
#undef GEOS_4
#undef A_LLK_03
#undef GRID30LEV
#undef GRID4x5
#undef GRID2x25  
#undef GRID1x125
#undef GRID1x1
#undef NESTED_NA
#undef NESTED_CH
#undef FULLCHEM  
#undef LGEOSCO
#undef LFASTJ
#undef LSLOWJ
#undef COMPAQ
#undef IBM_AIX
#undef LINUX_PGI
#undef LINUX_IFC
#undef LINUX_EFC
#undef SGI_MIPS
#undef SPARC

!==============================================================================
! Define the necessary "switches" for GEOS-CHEM. 
! Give each switch its own name as a value, since this will prevent
! the C-preprocessor overwriting the name everywhere in the code.
!==============================================================================

!----- Model types -----
!#define GEOS_1      'GEOS_1'       
!#define GEOS_STRAT  'GEOS_STRAT'
#define GEOS_3      'GEOS_3'
!#define GEOS_4      'GEOS_4'
!#define A_LLK_03    'A_LLK_03'

!----- Grid sizes -----
!#define GRID1x1     'GRID1x1'
!#define NESTED_CH   'NESTED_CH'
!#define NESTED_NA   'NESTED_NA'
!#define GRID1x125   'GRID1x125'
!#define GRID2x25    'GRID2x25'
#define GRID4x5     'GRID4x5'
#define GRID30LEV   'GRID30LEV'

!----- Chemistry -----
#define FULLCHEM    'FULLCHEM'
!#define LGEOSCO     'LGEOSCO'

!----- Photolysis -----
#define LFASTJ      'LFASTJ'
!#define LSLOWJ      'LSLOWJ'

!----- Compilers -----
!#define COMPAQ      'COMPAQ'
!#define IBM_AIX     'IBM_AIX'
!#define LINUX_PGI   'LINUX_PGI'
!#define LINUX_IFC   'LINUX_IFC'
!#define LINUX_EFC   'LINUX_EFC'
#define SGI_MIPS    'SGI_MIPS'
!#define SPARC       'SPARC'

!==============================================================================
! Force a compile error if GEOS_1, GEOS_STRAT, GEOS_3, GEOS_4 are undefined 
!==============================================================================
#if !defined( GEOS_1 ) && !defined( GEOS_STRAT ) && !defined( GEOS_3 ) && !defined( GEOS_4 )
#error "ERROR: GEOS_1, GEOS_STRAT, GEOS_3, and GEOS_4"
#error "are ALL undefined in header file define.h"
#endif

!==============================================================================
! Force a compile error if GRID1x1, GRID2x25, and GRID4x5 are all undefined 
!==============================================================================
#if !defined( GRID2x25 ) && !defined( GRID4x5 ) && !defined( GRID1x125 ) && !defined( GRID1x1 )
#error "ERROR: GRID4x5, GRID2x25, GRID1x125, and GRID1x1"
#error "are ALL undefined in header file define.h"
#endif

!==============================================================================
! Force a compile error if switches FULLCHEM, LGEOSCO are undefined
!==============================================================================
#if !defined( FULLCHEM ) && !defined( LGEOSCO )
#error "ERROR: One of FULLCHEM, LGEOSCO" 
#error "needs to be defined in header file define.h"
#endif

!==============================================================================
! Force a compile  error if all compiler switches are undefined
!==============================================================================
#if !defined(COMPAQ) && !defined(IBM_AIX) && !defined(LINUX_PGI) && !defined(LINUX_IFC) && !defined(LINUX_EFC) && !defined(SGI_MIPS) && ! defined(SPARC)
#error "ERROR: One of COMPAQ, IBM_AIX, LINUX_PGI,"
#error "LINUX_IFC, LINUX_EFC, SGI_MIPS, SPARC"
#error "needs to be defined in header file define.h"
#endif
