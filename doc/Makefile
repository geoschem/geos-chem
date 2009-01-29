# $Id: Makefile,v 1.4 2009/01/29 15:35:49 bmy Exp $
#==============================================================================
# Makefile: to create PDF and PS documentation from the ProTeX headers
#           contained in GEOS-Chem F90 module files. (bmy, 12/17/08, 1/28/09)
#==============================================================================

# List of "phony" targets (i.e. targets not associated with files)
.PHONY: all docclean

# Make all documentation
all: tpcore emissions

# Remove all Latex, PS, PDF files
clean:
	rm -f *.tex *.ps *.pdf

#---------------------------
# tpcore_fvdas_mod.f90
#---------------------------

# Documentation files
TPCORE_F90 =             \
tpcore_fvdas.intro       \
../tpcore_fvdas_mod.f90  \
../pjc_pfix_mod.f 

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


#---------------------------
# Emissions document
#---------------------------

# Documentation files
EMISS_F90 =                  \
emissions.intro              \
../emissions_mod.f           \
../arctas_ship_emiss_mod.f   \
../cac_anthro_mod.f          \
../emep_mod.f                \
../scale_anthro_mod.f        \
../vistas_anthro_mod.f 

EMISS_TEX = emissions.tex
EMISS_DVI = emissions.dvi
EMISS_PDF = emissions.pdf
EMISS_PS  = emissions.ps

# Build LaTeX, *.pdf, and *.ps documentation
emissions: 
	rm -f $(EMISS_TEX)
	protex -sf $(EMISS_F90) > $(EMISS_TEX)
	latex $(EMISS_TEX)
	latex $(EMISS_TEX)
	latex $(EMISS_TEX)
	dvipdf $(EMISS_DVI) $(EMISS_PDF)
	dvips $(EMISS_DVI) -o $(EMISS_PS)
	rm -f *.aux *.dvi *.log *.toc
