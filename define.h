! $Id: define.h,v 1.4 2003/11/06 21:07:18 bmy Exp $
!
!******************************************************************************
!  Include file "define.h" specifies C-preprocessor "switches" that are 
!  used to include or exclude certain sections of code.  
!  (bmy, bdf, 1/30/98, 10/31/03)
!
!  List of "Switches"
!  ===========================================================================
!  (1 ) GEOS_1     : Enables code for GEOS-1 met fields & chemistry
!  (2 ) GEOS_STRAT : Enables code for GEOS-STRAT met fields & chemistry
!  (3 ) GEOS_3     : Enables code for GEOS-3 met fields & chemistry
!  (4 ) GEOS_4     : Enables code for GEOS-4 met fields & chemistry
!  (5 ) GRID30LEV  : Enables code for 30-level GEOS-3 or GEOS-4 grid
!  (6 ) GRID1x1    : Enables code for 1 x 1   grid
!  (7 ) GRID2x25   : Enables code for 2 x 2.5 grid 
!  (8 ) GRID4x5    : Enables code for 4 x 5   grid 
!  (9 ) FULLCHEM   : Enables code for "Full" Chemistry (ISOP and NMHC)
!  (10) SMALLCHEM  : Enables code for "Small" chemistry (No ISOP, no NMHC)
!  (11) LGEOSCO    : Enables code for CO run w/ parameterized OH
!  (12) LFASTJ     : Enables code for FAST-J photolysis
!  (13) LSLOWJ     : Enables code for SLOW-J photolysis
!  (14) SGI        : Enables SGI specific code
!  (15) COMPAQ     : Enables COMPAQ Alpha specific code
!  (16) LINUX      : Enables Linux specific code
!  (17) SPARC      : Enables Sun/Sparc specific code
!  (18) IBM_AIX    : Enables IBM/AIX specific code
!  (19) INTEL_FC   : Enables INTEL FORTRAN COMPILER code
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
!******************************************************************************
!
!==============================================================================
! Undefine all "switches" so that they cannot be accidentally reset  
!==============================================================================
#undef GEOS_1
#undef GEOS_STRAT
#undef GEOS_3
#undef GEOS_4
#undef GRID30LEV
#undef GRID2x25  
#undef GRID4x5
#undef GRID1x1   
#undef FULLCHEM  
#undef SMALLCHEM 
#undef LGEOSCO
#undef LFASTJ
#undef LSLOWJ
#undef SGI
#undef COMPAQ
#undef LINUX
#undef SPARC
#undef IBM_AIX
#undef INTEL_FC

!==============================================================================
! Define the necessary "switches" for GEOS-CHEM. 
! Give each switch its own name as a value, since this will prevent
! the C-preprocessor overwriting the name everywhere in the code.
!==============================================================================
!#define GEOS_1      'GEOS_1'       
!#define GEOS_STRAT  'GEOS_STRAT'
#define GEOS_3      'GEOS_3'
!#define GEOS_4      'GEOS_4'

!#define GRID1x1     'GRID1x1'
!#define GRID2x25    'GRID2x25'
#define GRID4x5     'GRID4x5'
#define GRID30LEV   'GRID30LEV'

!#define SMALLCHEM   'SMALLCHEM'
#define FULLCHEM    'FULLCHEM'
!#define LGEOSCO     'LGEOSCO'

#define LFASTJ      'LFASTJ'
!#define LSLOWJ      'LSLOWJ'

!#define SGI         'SGI'
!#define COMPAQ      'COMPAQ'
!#define LINUX       'LINUX'
!#define SPARC       'SPARC'
!#define IBM_AIX     'IBM_AIX'
#define INTEL_FC    'INTEL_FC'

!==============================================================================
! Force a compile-time error if switches GEOS_1, GEOS_STRAT, 
! and GEOS_2, and GEOS_3 are all undefined. 
!==============================================================================
#if !defined( GEOS_1 ) && !defined( GEOS_STRAT ) && !defined( GEOS_3 ) && !defined( GEOS_4 )
#error "ERROR: GEOS_1, GEOS_STRAT, GEOS_3, and GEOS_4"
#error "are ALL undefined in header file define.h"
#endif

!==============================================================================
! Force a compile-time error if switches GRID1x1, GRID2x25,
! and GRID4x5 are all undefined. 
!==============================================================================
#if !defined( GRID2x25 ) && !defined( GRID4x5 ) && !defined( GRID1x1 )
#error "ERROR: GRID2x25, GRID4x5, and GRID1x1"
#error "are ALL undefined in header file define.h"
#endif

!==============================================================================
! Force a compile-time error if switches FULLCHEM, 
! SMALLCHEM, and LGEOSCO are all undefined
!==============================================================================
#if !defined( FULLCHEM ) && !defined( SMALLCHEM ) && !defined( LGEOSCO )
#error "ERROR: One of FULLCHEM, SMALLCHEM, LGEOSCO" 
#error "needs to be defined in header file define.h"
#endif

!==============================================================================
! Force a compile-time error if switches SGI, 
! DEC_COMPAQ, and LINUX are ALL undefined
!==============================================================================
#if !defined( SGI ) && !defined( COMPAQ ) && !defined( LINUX ) && !defined( SPARC ) && !defined( IBM_AIX ) && !defined( INTEL_FC )
#error "ERROR: One of SGI, DEC_COMPAQ, LINUX, SPARC, IBM_AIX, INTEL_FC"
#error "needs to be defined in header file define.h"
#endif
