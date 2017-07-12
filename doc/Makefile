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
#  09 Jul 2014 - R. Yantosca - Now build HEMCO documentation
#  21 Jul 2014 - R. Yantosca - Now look for HEMCO modules in subdirectories
#  10 Jul 2015 - R. Yantosca - Added GCRT directory
#  18 Nov 2016 - R. Yantosca - Updated for v11-01
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Define variables
ROOT  :=..
BIN   :=$(ROOT)/bin
BPCH  :=$(ROOT)/GeosBpch
CORE  :=$(ROOT)/GeosCore
DOC   :=$(ROOT)/doc
ESMF  :=$(ROOT)/ESMF
GCRT  :=$(ROOT)/GCRT
GTMM  :=$(ROOT)/GTMM
HCO   :=$(ROOT)/HEMCO/
HDR   :=$(ROOT)/Headers
HELP  :=$(ROOT)/help
HIST  :=$(ROOT)/History
ISO   :=$(ROOT)/ISOROPIA
NC    :=$(ROOT)/NcdfUtil
KPP   :=$(ROOT)/KPP
UTIL  :=$(ROOT)/GeosUtil

# Get the Unix shell in SHELL from the Makefile_header.mk
include $(ROOT)/Makefile_header.mk

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: doc docclean help

doc: utildoc srcdoc diagsdoc hemcodoc makedoc

clean:
	@echo "===> Making clean in directory: doc <==="
	@rm -f *.tex *.ps *.pdf *.aux *.dvi *.log *.toc

help:
	@$(MAKE) -C $(HELP)

###############################################################################
###                                                                         ###
###  Build the GEOS-Chem documentation                                      ###
###                                                                         ###
###############################################################################

## Commands to build the utility module documentation
include ./Makefile_UtilDoc.mk

# Commands to build the "core" modules documentation
include ./Makefile_SrcDoc.mk

# Commands to build the diagnostic modules documentation
#include ./Makefile_DiagsDoc.mk
include ./Makefile_HistoryDoc.mk

# Commands to build the HEMCO documentation
include ./Makefile_Hemco.mk

# Commands to build the makefile documentation
include ./Makefile_MakeDoc.mk


#EOC
