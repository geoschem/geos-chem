#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_MakeDoc.mk (in doc subdirectory)
#
# !DESCRIPTION: This Makefile fragment contains commands to build the
#  documentation for the GEOS-Chem Makefiles  It is inlined into
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
#  See https://github.com/geoschem/geos-chem for complete history
#EOP
#------------------------------------------------------------------------------
#BOC


# List of source code files (order is important)
SRC5 :=                          \
make_intro.P                     \
$(ROOT)/Makefile                 \
$(ROOT)/Makefile_header.mk       \
$(HCO)/Makefile                  \
$(UTIL)/Makefile                 \
$(ISO)/Makefile                  \
$(CORE)/Makefile                 \
$(GTMM)/Makefile                 \
$(KPP)/Makefile                  \
$(KPP)/SOA_SVPOA/Makefile        \
$(KPP)/Standard/Makefile         \
$(KPP)/Tropchem/Makefile         \
$(DOC)/Makefile                  \
$(DOC)/Makefile_UtilDoc.mk       \
$(DOC)/Makefile_SrcDoc.mk        \
$(DOC)/Makefile_HistoryDoc.mk    \
$(DOC)/Makefile_MakeDoc.mk       \
$(DOC)/Makefile_Hemco.mk         \
$(HELP)/Makefile

# Output file names
TEX5 :=GC_v11-02_Makefiles.tex
DVI5 :=GC_v11-02_Makefiles.dvi
PDF5 :=GC_v11-02_Makefiles.pdf

# Make command
makedoc:
	rm -f $(TEX5)
	./protex -sfS $(SRC5) > $(TEX5)
	latex $(TEX5)
	latex $(TEX5)
	latex $(TEX5)
	dvipdf $(DVI5) $(PDF5)
	rm -f *.aux *.dvi *.log *.toc

#EOC
