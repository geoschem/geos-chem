# $Id: Makefile,v 1.8 2009/09/21 18:59:35 bmy Exp $
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in doc subdirectory)
#
# !DESCRIPTION: Makefile for building the documentation (in PDF and 
#  PostScript formats) for the GEOS-Chem model source code, makefiles,
#  and related scripts.
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
#  16 Sep 2009 - R. Yantosca - Initial version
#  21 Sep 2009 - R. Yantosca - Now call Makefile in help directory to
#                              display the help screen options
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
SHELL   = /bin/sh
ROOTDIR = ..
BINDIR  = $(ROOTDIR)/bin
CODEDIR = $(ROOTDIR)/GeosCore
HELPDIR = $(ROOTDIR)/help
KPPDIR  = $(ROOTDIR)/KPP

#==============================================================================
# Makefile targets
#==============================================================================

.PHONY: all docclean help

all: srcdoc makedoc

clean:
	rm -f *.tex *.ps *.pdf

help:
	@$(MAKE) -C $(HELPDIR)

#==============================================================================
# Build documentation for GEOS-Chem source code
#==============================================================================

SRC1 =                               \
./geos_chem.intro                    \
$(CODEDIR)/tpcore_fvdas_mod.f90      \
$(CODEDIR)/pjc_pfix_mod.f            \
$(CODEDIR)/emissions_mod.f           \
$(CODEDIR)/arctas_ship_emiss_mod.f   \
$(CODEDIR)/bravo_mod.f               \
$(CODEDIR)/cac_anthro_mod.f          \
$(CODEDIR)/emep_mod.f                \
$(CODEDIR)/scale_anthro_mod.f        \
$(CODEDIR)/vistas_anthro_mod.f 

TEX1 = GEOS-Chem.tex
DVI1 = GEOS-Chem.dvi
PDF1 = GEOS-Chem.pdf
PS1  = GEOS-chem.ps

srcdoc: 
	rm -f $(TEX1)
	protex -sf $(SRC1) > $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	dvipdf $(DVI1) $(PDF1)
	dvips $(DVI1) -o $(PS1)
	rm -f *.aux *.dvi *.log *.toc


#==============================================================================
# Build documentation for GEOS-Chem source code
#==============================================================================

SRC2 =                         \
./make.intro                   \
$(ROOTDIR)/Makefile            \
$(ROOTDIR)/Makefile_header.mk  \
$(CODEDIR)/Makefile            \
$(KPPDIR)/Makefile             \
$(KPPDIR)/43t/Makefile         \
$(KPPDIR)/54t/Makefile         \
$(KPPDIR)/int/Makefile         \
./Makefile                     \
$(HELPDIR)/Makefile

TEX2 = GEOS-Chem-Makefiles.tex
DVI2 = GEOS-Chem-Makefiles.dvi
PDF2 = GEOS-Chem-Makefiles.pdf
PS2  = GEOS-Chem-Makefiles.ps

makedoc: 
	rm -f $(TEX2)
	protex -fS $(SRC2) > $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	dvipdf $(DVI2) $(PDF2)
	dvips $(DVI2) -o $(PS2)
	rm -f *.aux *.dvi *.log *.toc

#EOC
