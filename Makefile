# $Id: Makefile,v 1.1 2009/09/16 19:34:15 bmy Exp $
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
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# GEOSDIR    Specifies the directory where GEOS-Chem "core" routines are found
#
# !REVISION HISTORY: 
#  16 Sep 2009 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL   = /bin/sh
GEOSDIR = GeosCore

#=============================================================================
# Makefile targets: type "make help" for a complete list!
#=============================================================================

.PHONY: all lib exe clean realclean doc docclean help

all:
	@$(MAKE) -C $(GEOSDIR) all 

lib:
	@$(MAKE) -C $(GEOSDIR) lib

exe:
	@$(MAKE) -C $(GEOSDIR) all

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



