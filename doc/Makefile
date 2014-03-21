#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
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
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
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
#  27 Aug 2010 - R. Yantosca - Brought up to date for newly documented code
#  31 Aug 2010 - R. Yantosca - Add documentation for GTMM routines
#  14 Sep 2010 - R. Yantosca - Now split make commands into include files
#  20 Aug 2013 - R. Yantosca - Remove GeosTomas directory, that is obsolete
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
ROOTDIR := ..
BIN     := $(ROOTDIR)/bin
BPCH    := $(ROOTDIR)/GeosBpch
CORE    := $(ROOTDIR)/GeosCore
DOC     := $(ROOTDIR)/doc
ESMF    := $(ROOTDIR)/ESMF
GTMM    := $(ROOTDIR)/GTMM
HDR     := $(ROOTDIR)/Headers
HELP    := $(ROOTDIR)/help
ISO     := $(ROOTDIR)/ISOROPIA
KPP     := $(ROOTDIR)/KPP
UTIL    := $(ROOTDIR)/GeosUtil

# Get the Unix shell in SHELL from the Makefile_header.mk
include $(ROOTDIR)/Makefile_header.mk

#==============================================================================
# Makefile targets
#==============================================================================

.PHONY: doc docclean help

doc: srcdoc utildoc gtmmdoc makedoc

clean:
	rm -f *.tex *.ps *.pdf

help:
	@$(MAKE) -C $(HELPDIR)

#==============================================================================
# Build the GEOS-Chem documentation
#==============================================================================

# Commands to build the source code documentation
include ./Makefile_SrcDoc.mk

# Commands to build the utility module documentation
include ./Makefile_UtilDoc.mk

# Commands to build the makefile documentation
include ./Makefile_MakeDoc.mk

# Commands to build the GTMM documentation
include ./Makefile_GtmmDoc.mk


#EOC
