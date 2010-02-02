# $Id: Makefile,v 1.12 2010/02/02 16:57:45 bmy Exp $
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
#  11 Dec 2009 - R. Yantosca - Now get SHELL from Makefile_header.mk
#  28 Jan 2010 - R. Yantosca - Now add the Makefiles in the GeosTomas and
#                              ISOROPIA subdirs to the documentation
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
ROOTDIR = ..
BIN     = $(ROOTDIR)/bin
BPCH    = $(ROOTDIR)/GeosBpch
CORE    = $(ROOTDIR)/GeosCore
HELP    = $(ROOTDIR)/help
ISO     = $(ROOTDIR)/ISOROPIA
KPP     = $(ROOTDIR)/KPP
TOM     = $(ROOTDIR)/GeosTomas
UTIL    = $(ROOTDIR)/GeosUtil

# Get the Unix shell in SHELL from the Makefile_header.mk
include $(ROOTDIR)/Makefile_header.mk

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
$(CORE)/tpcore_fvdas_mod.f90      \
$(CORE)/pjc_pfix_mod.f            \
$(CORE)/emissions_mod.f           \
$(CORE)/arctas_ship_emiss_mod.f   \
$(CORE)/bravo_mod.f               \
$(CORE)/cac_anthro_mod.f          \
$(CORE)/emep_mod.f                \
$(CORE)/scale_anthro_mod.f        \
$(CORE)/vistas_anthro_mod.f 

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
$(UTIL)/Makefile               \
$(ISO)/Makefile                \
$(CODE)/Makefile               \
$(KPP)/Makefile                \
$(KPP)/standard/Makefile       \
$(KPP)/SOA/Makefile            \
$(TOM)/Makefile                \
./Makefile                     \
$(HELP)/Makefile

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

SRC3 = ./util.intro      \
$(UTIL)/bpch2_mod.f      \
$(UTIL)/charpak_mod.f    \
$(UTIL)/directory_mod.f  \
$(UTIL)/error_mod.f      \
$(UTIL)/ifort_errmsg.f   \
$(UTIL)/file_mod.f       \
$(UTIL)/grid_mod.f       \
$(UTIL)/hdf_mod.f        \
$(UTIL)/julday_mod.f     \
$(UTIL)/pressure_mod.f   \
$(UTIL)/regrid_1x1_mod.f \
$(UTIL)/time_mod.f       \
$(UTIL)/transfer_mod.f   \
$(UTIL)/unix_cmds_mod.f

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

#EOC
