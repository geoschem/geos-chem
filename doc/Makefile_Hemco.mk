#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_Hemco.mk (in doc subdirectory)
#
# !DESCRIPTION: This Makefile fragment contains commands to build the
#  documentation for the HEMCO Source Code.  It is inlined into
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
#  08 Jul 2014 - R. Yantosca - Initial version
#  See https://github.com/geoschem/geos-chem for complete history
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files (order is important)
SRC4 :=                            \
./hemco_intro.P                    \
./hemco_core.P                     \
$(wildcard $(HCO)/Core/*.F*)       \
./hemco_extensions.P               \
$(wildcard $(HCO)/Extensions/*.F*) \
./hemco_interfaces.P               \
$(wildcard $(HCO)/Interfaces/*.F*)

# Output file names
TEX4 := GC_v11-02_HEMCO_Modules.tex
DVI4 := GC_v11-02_HEMCO_Modules.dvi
PDF4 := GC_v11-02_HEMCO_Modules.pdf

# Make commands
hemcodoc:
	rm -f $(TEX4)
	./protex -sfp $(SRC4) > $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	dvipdf $(DVI4) $(PDF4)
	rm -f *.aux *.dvi *.log *.toc

#EOC
