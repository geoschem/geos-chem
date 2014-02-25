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
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files
SRC3 :=                     \
./intro.util                \
$(wildcard $(UTIL)/*.F)     \
$(wildcard $(UTIL)/*.F90)


# Output file names
TEX3 := GC_Ref_Vol_2.tex
DVI3 := GC_Ref_Vol_2.dvi
PDF3 := GC_Ref_Vol_2.pdf


# Make commands
utildoc: 
	rm -f $(TEX3)
	protex -sf $(SRC3) > $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	dvipdf $(DVI3) $(PDF3)
	rm -f *.aux *.dvi *.log *.toc

#EOC
