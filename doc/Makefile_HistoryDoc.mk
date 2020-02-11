#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_HistoryDoc.mk (in doc subdirectory)
#
# !DESCRIPTION: This Makefile fragment contains commands to build the
#  documentation for the GEOS-Chem History modules.  It is inlined into
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
#  12 Jul 2017 - R. Yantosca - Initial version
#  See https://github.com/geoschem/geos-chem for complete history
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files
SRC3 :=                           \
./gc_history_intro.P              \
./gc_history.P                    \
$(HIST)/history_params_mod.F90    \
$(HIST)/histitem_mod.F90          \
$(HIST)/metahistitem_mod.F90      \
$(HIST)/histcontainer_mod.F90     \
$(HIST)/metahistcontainer_mod.F90 \
$(HIST)/history_mod.F90

# Output file names
TEX3 := GC_v11-02_History_Modules.tex
DVI3 := GC_v11-02_History_Modules.dvi
PDF3 := GC_v11-02_History_Modules.pdf

# Make commands
historydoc:
	rm -f $(TEX3)
	./protex -sfp $(SRC3) > $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	dvipdf $(DVI3) $(PDF3)
	rm -f *.aux *.dvi *.log *.toc

#EOC
