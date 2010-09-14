#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
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
#
#   make TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# You must have the LaTeX utilities (latex, dvips, dvipdf) installed
# on your system in order to build the documentation.
#
# !REVISION HISTORY: 
#  14 Sep 2010 - R. Yantosca - Initial version, split off from Makefile
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files
SRC3 := ./intro.util $(wildcard $(UTIL)/*.f)


# Output file names
TEX3 := GEOS-Chem-Utility.tex
DVI3 := GEOS-Chem-Utility.dvi
PDF3 := GEOS-Chem-Utility.pdf
PS3  := GEOS-Chem-Utility.ps


# Make commands
utildoc: 
	rm -f $(TEX3)
	protex -sf $(SRC3) > $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	dvipdf $(DVI3) $(PDF3)
	dvips $(DVI3) -o $(PS3)
	rm -f *.aux *.dvi *.log *.toc

#EOC
