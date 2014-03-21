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
#  16 Dec 2010 - R. Yantosca - Renamed output files to "GC_Ref_Vol_1.*"
#  15 Jan 2014 - R. Yantosca - Now only create *.pdf output
#  15 Jan 2014 - R. Yantosca - Now only prints prologues, avoids printing code
#EOP
#------------------------------------------------------------------------------
#BOC


# List of source code files (order is important)
SRC2 :=                        \
./intro.make                   \
$(ROOTDIR)/Makefile            \
$(ROOTDIR)/Makefile_header.mk  \
$(UTIL)/Makefile               \
$(ISO)/Makefile                \
$(CORE)/Makefile               \
$(KPP)/Makefile                \
$(KPP)/standard/Makefile       \
$(KPP)/SOA/Makefile            \
$(GTMM)/Makefile               \
$(DOC)/Makefile                \
$(DOC)/Makefile_SrcDoc.mk      \
$(DOC)/Makefile_UtilDoc.mk     \
$(DOC)/Makefile_GtmmDoc.mk     \
$(DOC)/Makefile_MakeDoc.mk     \
$(HELP)/Makefile


# Output file names
TEX2 := GC_Ref_Vol_1.tex
DVI2 := GC_Ref_Vol_1.dvi
PDF2 := GC_Ref_Vol_1.pdf


# Make command
makedoc: 
	rm -f $(TEX2)
	protex -sfS $(SRC2) > $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	dvipdf $(DVI2) $(PDF2)
	rm -f *.aux *.dvi *.log *.toc

#EOC
