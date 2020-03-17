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
#  See https://github.com/geoschem/geos-chem for complete history
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
CORE  :=$(ROOT)/GeosCore
DOC   :=$(ROOT)/doc
ESMF  :=$(ROOT)/ESMF
GTMM  :=$(ROOT)/GTMM
HCO   :=$(ROOT)/HEMCO/
HDR   :=$(ROOT)/Headers
HELP  :=$(ROOT)/help
HIST  :=$(ROOT)/History
ISO   :=$(ROOT)/ISORROPIA
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

doc: utildoc srcdoc historydoc hemcodoc makedoc

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

# Commands to build the History modules documentation
include ./Makefile_HistoryDoc.mk

# Commands to build the HEMCO documentation
include ./Makefile_Hemco.mk

# Commands to build the makefile documentation
include ./Makefile_MakeDoc.mk


#EOC
