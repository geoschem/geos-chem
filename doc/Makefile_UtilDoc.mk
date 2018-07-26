#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_UtilDoc.mk (in doc subdirectory)
#
# !DESCRIPTION: This Makefile fragment contains commands to build the 
#  documentation for the GEOS-Chem utility modules.  It is inlined into
#  the Makefile (in the doc subdirectory) by an "include" command.
#\\
#\\
# !REMARKS:
# To build the documentation, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# You must have the LaTeX utilities (latex, dvips, dvipdf) installed
# on your system in order to build the documentation.
#
# !REVISION HISTORY: 
#  14 Sep 2010 - R. Yantosca - Initial version, split off from Makefile
#  16 Dec 2010 - R. Yantosca - Renamed output files to "GC_Ref_Vol_2.*"
#  19 Jul 2011 - R. Yantosca - Changed *.f* to *.F* for ESMF compatibility
#  03 Apr 2012 - M. Payer    - Added *.F90 so that grid_mod.F90 and 
#                              regrid_a2a_mod.F90 are included
#  15 Jan 2014 - R. Yantosca - Now only create *.pdf output
#  10 Jul 2015 - R. Yantosca - Use ./protex to avoid problems on some systems
#  18 Nov 2016 - R. Yantosca - Now use protex -sfp to increment chapters
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files
SRC1 :=                     \
./util_intro.P              \
./util_grid.P               \
$(UTIL)/pressure_mod.F      \
$(UTIL)/gc_grid_mod.F90     \
$(UTIL)/regrid_a2a_mod.F90  \
$(UTIL)/bpch2_mod.F         \
$(UTIL)/transfer_mod.F      \
$(UTIL)/file_mod.F          \
./util_file.P               \
$(wildcard $(NC)/*.F*)      \
$(SRCNC)                    \
./util_time.P               \
$(UTIL)/time_mod.F          \
$(UTIL)/julday_mod.F        \
./util_units.P              \
$(UTIL)/unitconv_mod.F90    \
./util_error.P              \
$(UTIL)/error_mod.F         \
$(UTIL)/ifort_errmsg.F      \
./util_misc.P               \
$(UTIL)/geos_timers_mod.F90

# Output file names
TEX1 := GC_v11-02_Utility_Modules.tex
DVI1 := GC_v11-02_Utility_Modules.dvi
PDF1 := GC_v11-02_Utility_Modules.pdf

# Make commands
utildoc: 
	rm -f $(TEX1)
	./protex -sfp $(SRC1) > $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	dvipdf $(DVI1) $(PDF1)
	rm -f *.aux *.dvi *.log *.toc

#EOC
