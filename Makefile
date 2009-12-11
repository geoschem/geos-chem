# $Id: Makefile,v 1.8 2009/12/11 19:27:42 bmy Exp $
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (Main-level)
#
# !DESCRIPTION: This is a "router" makefile.  It calls the main GEOS-Chem 
# Makefile (in the GeosCore subdirectory) to direct the Unix "make" utility 
# how to build the GEOS-Chem source code.
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# GEOSDIR    Specifies the directory where GEOS-Chem "core" routines are found
#
# !REVISION HISTORY: 
#  16 Sep 2009 - R. Yantosca - Initial version
#  24 Nov 2009 - R. Yantosca - Now call libbpch and libcore targets in
#                              the Makefile in the GeosCore sub-directory
#  11 Dec 2009 - R. Yantosca - Now get SHELL from Makefile_header.mk
#EOP
#------------------------------------------------------------------------------
#BOC

# Get the Unix shell definition
include ./Makefile_header.mk

# Define variables
GEOSDIR = GeosCore

#=============================================================================
# Makefile targets: type "make help" for a complete list!
#=============================================================================

.PHONY: all lib libkpp libutil exe clean realclean doc docclean help

all:
	@$(MAKE) -C $(GEOSDIR) all

lib:
	@$(MAKE) -C $(GEOSDIR) lib

libcore:
	@$(MAKE) -C $(GEOSDIR) libcore

libkpp:
	@$(MAKE) -C $(GEOSDIR) libkpp

libutil:
	@$(MAKE) -C $(GEOSDIR) libutil

exe:
	@$(MAKE) -C $(GEOSDIR) exe

clean:
	@$(MAKE) -C $(GEOSDIR) clean

realclean:
	@$(MAKE) -C $(GEOSDIR) realclean

doc:
	@$(MAKE) -C $(GEOSDIR) doc

docclean: 
	@$(MAKE) -C $(GEOSDIR) docclean

help:
	@$(MAKE) -C $(GEOSDIR) help

#EOC



