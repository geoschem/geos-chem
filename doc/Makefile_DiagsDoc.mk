#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_DiagsDoc.mk (in doc subdirectory)
#
# !DESCRIPTION: This Makefile fragment contains commands to build the 
#  documentation for the GEOS-Chem diagnostics modules.  It is inlined into
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
#  18 Nov 2016 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files
SRC3 :=                               \
./gc_diags_intro.P                    \
./gc_diags.P                          \
$(CORE)/benchmark_mod.F               \
$(CORE)/diag1.F                       \
$(CORE)/diag_2pm.F                    \
$(CORE)/diag3.F                       \
$(CORE)/diag_mod.F                    \
$(CORE)/diag03_mod.F                  \
$(CORE)/diag04_mod.F                  \
$(CORE)/diag20_mod.F                  \
$(CORE)/diag41_mod.F                  \
$(CORE)/diag42_mod.F                  \
$(CORE)/diag49_mod.F                  \
$(CORE)/diag50_mod.F                  \
$(CORE)/diag51b_mod.F                 \
$(CORE)/diag53_mod.F                  \
$(CORE)/diag56_mod.F                  \
$(CORE)/diag63_mod.F                  \
$(CORE)/diagoh.F                      \
$(CORE)/diag_oh_mod.F                 \
$(CORE)/gamap_mod.F                   \
$(CORE)/initialize.F                  \
$(CORE)/ndxx_setup.F                  \
$(CORE)/ohsave.F                      \
$(CORE)/planeflight_mod.F             \

# Output file names
TEX3 := GC_v11-01_Diagnostic_Modules.tex
DVI3 := GC_v11-01_Diagnostic_Modules.dvi
PDF3 := GC_v11-01_Diagnostic_Modules.pdf

# Make commands
diagsdoc:
	rm -f $(TEX3)
	./protex -sfp $(SRC3) > $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	dvipdf $(DVI3) $(PDF3)
	rm -f *.aux *.dvi *.log *.toc

#EOC
