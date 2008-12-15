# $Id: Makefile,v 1.1 2008/12/15 15:55:14 bmy Exp $
#==============================================================================
# Makefile: to create PDF and PS documentation from the ProTeX headers
#           contained in GEOS-Chem F90 module files. (bmy, 12/15/08)
#==============================================================================

# Make all documentation
all: tpcore

#---------------------------
# tpcore_fvdas_mod.f90
#---------------------------

# Documentation files
TPCORE_F90 = tpcore_fvdas.intro ../tpcore_fvdas_mod.f90 ../pjc_pfix_mod.f 
TPCORE_TEX = tpcore_fvdas.tex
TPCORE_DVI = tpcore_fvdas.dvi
TPCORE_PDF = tpcore_fvdas.pdf
TPCORE_PS  = tpcore_fvdas.ps

# Build LaTeX, *.pdf, and *.ps documentation
tpcore: 
	rm -f $(TPCORE_TEX)
	protex -sf $(TPCORE_F90) > $(TPCORE_TEX)
	latex $(TPCORE_TEX)
	latex $(TPCORE_TEX)
	latex $(TPCORE_TEX)
	dvipdf $(TPCORE_DVI) $(TPCORE_PDF)
	dvips $(TPCORE_DVI) -o $(TPCORE_PS)
	rm -f *.aux *.dvi *.log *.toc