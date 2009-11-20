# $Id: Makefile,v 1.10 2009/11/20 21:43:00 bmy Exp $
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
#  19 Nov 2009 - R. Yantosca - Now build documentation for Makefiles in the
#                              GeosUtil and GeosCore subdirectories
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
BPCHDIR = $(ROOTDIR)/GeosBpch
CODEDIR = $(ROOTDIR)/GeosCore
HELPDIR = $(ROOTDIR)/help
KPPDIR  = $(ROOTDIR)/KPP
UTILDIR = $(ROOTDIR)/GeosUtil

#==============================================================================
# Makefile targets
#==============================================================================

.PHONY: all docclean help

all: srcdoc makedoc utildoc bpchdoc

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
# Build documentation for GEOS-Chem Makefiles
#==============================================================================

SRC2 =                         \
./make.intro                   \
$(ROOTDIR)/Makefile            \
$(ROOTDIR)/Makefile_header.mk  \
$(UTILDIR)/Makefile            \
$(BPCHDIR)/Makefile            \
$(CODEDIR)/Makefile            \
$(KPPDIR)/Makefile             \
$(KPPDIR)/43t/Makefile         \
$(KPPDIR)/54t/Makefile         \
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

#==============================================================================
# Build documentation for GEOS-Chem utility modules
#==============================================================================

SRC3 = ./util.intro         \
$(UTILDIR)/charpak_mod.f    \
$(UTILDIR)/directory_mod.f  \
$(UTILDIR)/error_mod.f      \
$(UTILDIR)/ifort_errmsg.f   \
$(UTILDIR)/file_mod.f       \
$(UTILDIR)/grid_mod.f       \
$(UTILDIR)/julday_mod.f     \
$(UTILDIR)/pressure_mod.f   \
$(UTILDIR)/regrid_1x1_mod.f \
$(UTILDIR)/time_mod.f       \
$(UTILDIR)/transfer_mod.f   \
$(UTILDIR)/unix_cmds_mod.f

TEX3 = GEOS-Chem-Utility.tex
DVI3 = GEOS-Chem-Utility.dvi
PDF3 = GEOS-Chem-Utility.pdf
PS3  = GEOS-Chem-Utility.ps

utildoc: 
	rm -f $(TEX3)
	protex -sf $(SRC3) > $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	dvipdf $(DVI3) $(PDF3)
	dvips $(DVI3) -o $(PS3)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for GEOS-Chem utility modules
#==============================================================================

SRC4 = ./bpch.intro    \
$(BPCHDIR)/bpch2_mod.f

TEX4 = GEOS-Chem-Bpch.tex
DVI4 = GEOS-Chem-Bpch.dvi
PDF4 = GEOS-Chem-Bpch.pdf
PS4  = GEOS-Chem-Bpch.ps

bpchdoc: 
	rm -f $(TEX4)
	protex -sf $(SRC4) > $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	dvipdf $(DVI4) $(PDF4)
	dvips $(DVI4) -o $(PS4)
	rm -f *.aux *.dvi *.log *.toc


#EOC
