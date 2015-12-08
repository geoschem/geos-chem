#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the GeosRad/ subdirectory)
#
# !DESCRIPTION: This makefile compiles the RAD codes.  
#  the MOD directory. 
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in the directory   %%%
# %%% just above this one!                                                 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root-level directory of the GEOS-Chem code
# HDR        Specifies the directory where GEOS-Chem include files are found
# LIB        Specifies the directory where library files (*.a) are stored
# MOD        Specifies the directory where module files (*.mod) are stored
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
# R8         Interprets REAL as REAL*8
#
# !REVISION HISTORY: 
#  08 Aug 2012 - K. Cady-Pereira - Initial version
#  22 Aug 2012 - R. Yantosca     - Fixed typo in depend list for rrsw_kg24.o
#  22 Aug 2012 - R. Yantosca     - Added white space for readability
#  04 Jun 2015 - R. Yantosca     - Cosmetic changes
#  04 Jun 2015 - R. Yantosca     - Remove *.mod, *.a files with "make clean"
#  04 Jun 2015 - R. Yantosca     - Update debug target, remove help
#  04 Jun 2015 - R. Yantosca     - Manually remove mersennetwister.mod
#  07 Dec 2015 - R. Yantosca     - Restore fast "clean" command; add "slowclean"
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT    := ..
LIB     := $(ROOT)/lib
MOD     := $(ROOT)/mod

# Sometimes we may need to compile RRTMG with a lower optimization setting,
# so use the OPT variable if necessary
OPT     := -O1 -override-limits

# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source files
SOURCES := $(wildcard *.F90)

# List of object files (replace .F90 extensions with *.o)
OBJECTS := $(SOURCES:.F90=.o)

# List of module files.  Convert to lowercase, then prefix directory name.
MODULES :=$(OBJECTS:.o=.mod)
MODULES :=$(shell echo $(MODULES) | tr A-Z a-z)
MODULES :=$(foreach I,$(MODULES),$(MOD)/$(I))

# Bug fix for RRTMG: the file mcica_random_numbers.F90 contains a module 
# named MersenneTwister.  Make sure to also remove this module as well,
# by manually appending it to the modules list.
MODULES += $(MOD)/mersennetwister.mod

# Library file
LIBRARY :=librad.a

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug slowclean

lib: $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB)

clean:
	@echo "===> Making clean in directory: GeosRad <==="
	@rm -f *.o *.mod *.a *.x

slowclean:
	@echo "===> Making slowclean in directory: GeosRad <==="
	@rm -f $(OBJECTS) $(MODULES) $(LIBRARY) $(LIB)/$(LIBRARY)

debug:
	@echo "Targets : $(MAKECMDGOALS)"
	@echo "ROOT    : $(ROOT)"
	@echo "LIB     : $(LIB)"
	@echo "MOD     : $(MOD)"
	@echo "F90     : $(F90)"
	@echo "OBJECTS : $(OBJECTS)"
	@echo "MODULES : $(MODULES)"
	@echo "LIBRARY : $(LIBRARY)"

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###  (grep "USE " to get the list of module references!)                    ###
###                                                                         ###
###  From this list of dependencies, the "make" utility will figure out     ###
###  correct order of compilation (so we don't have to do that ourselves).  ###
###  This also allows us to compile on multiple processors with "make -j".  ###
###                                                                         ###
###  NOTES:                                                                 ###
###  (1) Only specify object-file dependencies that are within this         ###
###       directory.  Object files in other directories will be referenced  ### 
###       at link-time.                                                     ###
###  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
###       in this directory must have a listed dependency.                  ###
###                                                                         ###
###############################################################################

mcica_random_numbers.o   : mcica_random_numbers.F90                          \
                           parkind.o
	$(F90) $(R8) -c $(OPT) $<

mcica_subcol_gen_lw.o    : mcica_subcol_gen_lw.F90  mcica_random_numbers.o   \
                           parkind.o                parrrtm.o                \
                           rrlw_con.o               rrlw_wvn.o               \
                           rrlw_vsn.o
	$(F90) $(R8) -c $(OPT) $<


mcica_subcol_gen_sw.o    : mcica_subcol_gen_sw.F90  mcica_random_numbers.o   \
                           parkind.o                parrrsw.o                \
                           rrsw_con.o               rrsw_wvn.o               \
                           rrsw_vsn.o
	$(F90) $(R8) -c $(OPT) $<

parkind.o                : parkind.F90 
	$(F90) $(R8) -c $(OPT) $<


parrrtm.o                : parrrtm.F90             parkind.o 
	$(F90) $(R8) -c $(OPT) $<


#rad_driver.o            : rad_driver.F90                                    \
                           parkind.o               parrrtm.o                 \
                           parrrsw.o               rrlw_con.o                \
                           rrtmg_lw_rad.o          rrtmg_sw_rad.o
#	$(F90) $(R8) -c $(OPT) $<


rrlw_cld.o               : rrlw_cld.F90            parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_con.o               : rrlw_con.F90            parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg01.o              : rrlw_kg01.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg02.o              : rrlw_kg02.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg03.o              : rrlw_kg03.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg04.o              : rrlw_kg04.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg05.o               : rrlw_kg05.F90          parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg06.o               : rrlw_kg06.F90          parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg07.o               : rrlw_kg07.F90          parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg08.o              : rrlw_kg08.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg09.o              : rrlw_kg09.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg10.o              : rrlw_kg10.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg11.o              : rrlw_kg11.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg12.o              : rrlw_kg12.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg13.o              : rrlw_kg13.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg14.o              : rrlw_kg14.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg15.o              : rrlw_kg15.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_kg16.o              : rrlw_kg16.F90           parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_ncpar.o             : rrlw_ncpar.F90          parkind.o 
	$(F90) $(R8) -c $(OPT) $<

rrlw_ref.o               : rrlw_ref.F90            parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_tbl.o               : rrlw_tbl.F90            parkind.o 
	$(F90) $(R8) -c $(OPT) $<


rrlw_vsn.o               : rrlw_vsn.F90 
	$(F90) $(R8) -c $(OPT) $<


rrlw_wvn.o               : rrlw_wvn.F90                                      \
                           parkind.o               parrrtm.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_lw_cldprmc.o       : rrtmg_lw_cldprmc.F90                              \
                           parkind.o               parrrtm.o                 \
                           rrlw_cld.o              rrlw_wvn.o                \
                           rrlw_vsn.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_lw_init.o          : rrtmg_lw_init.F90                                 \
                           parkind.o               parrrtm.o                 \
                           rrlw_cld.o              rrlw_vsn.o                \
                           rrlw_wvn.o              rrlw_tbl.o                \
                           rrlw_con.o              rrtmg_lw_setcoef.o        \
                           rrlw_kg01.o             rrlw_kg02.o               \
                           rrlw_kg03.o             rrlw_kg04.o               \
                           rrlw_kg05.o             rrlw_kg06.o               \
                           rrlw_kg07.o             rrlw_kg08.o               \
                           rrlw_kg09.o             rrlw_kg10.o               \
                           rrlw_kg11.o             rrlw_kg12.o               \
                           rrlw_kg13.o             rrlw_kg14.o               \
                           rrlw_kg15.o             rrlw_kg16.o  
	$(F90) $(R8) -c $(OPT) $<


rrtmg_lw_k_g.o           : rrtmg_lw_k_g.F90                                  \
                           parkind.o               rrlw_vsn.o                \
                           rrlw_kg01.o             rrlw_kg02.o               \
                           rrlw_kg03.o             rrlw_kg04.o               \
                           rrlw_kg05.o             rrlw_kg06.o               \
                           rrlw_kg07.o             rrlw_kg08.o               \
                           rrlw_kg05.o             rrlw_kg06.o               \
                           rrlw_kg07.o             rrlw_kg08.o               \
                           rrlw_kg09.o             rrlw_kg10.o               \
                           rrlw_kg11.o             rrlw_kg12.o               \
                           rrlw_kg13.o             rrlw_kg14.o               \
                           rrlw_kg15.o             rrlw_kg16.o   
	$(F90) $(R8) -c $(OPT) $<


rrtmg_lw_rad.o           : rrtmg_lw_rad.F90                                  \
                           parkind.o               parrrtm.o                 \
                           rrlw_cld.o              rrlw_vsn.o                \
                           rrlw_wvn.o              rrlw_tbl.o                \
                           rrlw_con.o              rrtmg_lw_setcoef.o        \
                           rrtmg_lw_taumol.o       mcica_subcol_gen_lw.o     \
                           rrtmg_lw_rtrnmc.o       rrlw_wvn.o                \
                           rrlw_tbl.o              rrtmg_lw_cldprmc.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_lw_rtrnmc.o        : rrtmg_lw_rtrnmc.F90                               \
                           parkind.o               parrrtm.o                 \
                           rrlw_vsn.o              rrlw_wvn.o                \
                           rrlw_tbl.o              rrlw_con.o   
	$(F90) $(R8) -c $(OPT) $<


rrtmg_lw_setcoef.o       : rrtmg_lw_setcoef.F90                              \
                           parkind.o               parrrtm.o                 \
                           rrlw_vsn.o              rrlw_wvn.o                \
                           rrlw_ref.o     
	$(F90) $(R8) -c $(OPT) $<


rrtmg_lw_taumol.o        : rrtmg_lw_taumol.F90                               \
                           parkind.o               rrlw_vsn.o                \
                           rrlw_vsn.o              rrlw_wvn.o                \
                           rrlw_con.o              rrlw_ref.o                \
                           rrlw_kg01.o             rrlw_kg02.o               \
                           rrlw_kg03.o             rrlw_kg04.o               \
                           rrlw_kg05.o             rrlw_kg06.o               \
                           rrlw_kg07.o             rrlw_kg08.o               \
                           rrlw_kg09.o             rrlw_kg10.o               \
                           rrlw_kg11.o             rrlw_kg12.o               \
                           rrlw_kg13.o             rrlw_kg14.o               \
                           rrlw_kg15.o             rrlw_kg16.o  
	$(F90) $(R8) -c $(OPT) $<


parrrsw.o                : parrrsw.F90             parkind.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_aer.o               : rrsw_aer.F90                                      \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_cld.o               : rrsw_cld.F90            parkind.o
	$(F90) $(R8) -c $(OPT) $<



rrsw_con.o               : rrsw_con.F90            parkind.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg16.o              : rrsw_kg16.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg17.o              : rrsw_kg17.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg18.o              : rrsw_kg18.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg19.o              : rrsw_kg19.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg20.o              : rrsw_kg20.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg21.o              : rrsw_kg21.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg22.o              : rrsw_kg22.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg23.o              : rrsw_kg23.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg24.o              : rrsw_kg24.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg25.o              : rrsw_kg25.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg26.o              : rrsw_kg26.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg27.o              : rrsw_kg27.F90           parkind.o                 \
                           parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg28.o              : rrsw_kg28.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_kg29.o              : rrsw_kg29.F90                                     \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_ref.o               : rrsw_ref.F90            parkind.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_tbl.o               : rrsw_tbl.F90            parkind.o
	$(F90) $(R8) -c $(OPT) $<


rrsw_vsn.o               : rrsw_vsn.F90
	$(F90) $(R8) -c $(OPT) $<


rrsw_wvn.o               : rrsw_wvn.F90                                      \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_sw_cldprmc.o       : rrtmg_sw_cldprmc.F90                              \
                           parkind.o               parrrsw.o                 \
                           rrsw_cld.o              rrsw_wvn.o                \
                           rrsw_vsn.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_sw_init.o          : rrtmg_sw_init.F90                                 \
                           parkind.o               parrrsw.o                 \
                           rrsw_cld.o              rrsw_wvn.o                \
                           rrsw_vsn.o              rrsw_tbl.o                \
                           rrsw_con.o              rrsw_aer.o                \
                           rrsw_con.o              rrtmg_sw_setcoef.o        \
                           rrsw_kg16.o             rrsw_kg17.o               \
                           rrsw_kg18.o             rrsw_kg19.o               \
                           rrsw_kg20.o             rrsw_kg21.o               \
                           rrsw_kg22.o             rrsw_kg23.o               \
                           rrsw_kg24.o             rrsw_kg25.o               \
                           rrsw_kg26.o             rrsw_kg27.o               \
                           rrsw_kg28.o             rrsw_kg29.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_sw_k_g.o           : rrtmg_sw_k_g.F90                                  \
                           parkind.o               rrsw_kg16.o               \
                           rrsw_kg17.o             rrsw_kg18.o               \
                           rrsw_kg19.o             rrsw_kg20.o               \
                           rrsw_kg21.o             rrsw_kg22.o               \
                           rrsw_kg23.o             rrsw_kg24.o               \
                           rrsw_kg25.o             rrsw_kg26.o               \
                           rrsw_kg27.o             rrsw_kg28.o               \
                           rrsw_kg29.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_sw_rad.o           : rrtmg_sw_rad.F90                                  \
                           parkind.o               parrrsw.o                 \
                           rrsw_wvn.o              rrsw_vsn.o                \
                           rrsw_con.o              rrsw_aer.o                \
                           rrsw_con.o              rrtmg_sw_cldprmc.o        \
                           rrtmg_sw_init.o         rrtmg_sw_setcoef.o        \
                           rrtmg_sw_spcvmc.o       mcica_subcol_gen_sw.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_sw_reftra.o        : rrtmg_sw_reftra.F90                               \
                           parkind.o               rrsw_vsn.o                \
                           rrsw_tbl.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_sw_setcoef.o       : rrtmg_sw_setcoef.F90                              \
                           parkind.o               rrsw_vsn.o                \
                           rrsw_tbl.o              parrrsw.o                 \
                           rrsw_ref.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_sw_spcvmc.o        : rrtmg_sw_spcvmc.F90                               \
                           parkind.o               parrrsw.o                 \
                           rrsw_wvn.o              rrsw_vsn.o                \
                           rrsw_tbl.o              rrtmg_sw_reftra.o         \
                           rrtmg_sw_taumol.o       rrtmg_sw_vrtqdr.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_sw_taumol.o        : rrtmg_sw_taumol.F90                               \
                           parkind.o               parrrsw.o                 \
                           rrsw_con.o              rrsw_wvn.o                \
                           rrsw_vsn.o              rrsw_kg16.o               \
                           rrsw_kg17.o             rrsw_kg18.o  	     \
                           rrsw_kg19.o             rrsw_kg20.o   	     \
                           rrsw_kg21.o             rrsw_kg22.o  	     \
                           rrsw_kg23.o             rrsw_kg24.o   	     \
                           rrsw_kg25.o             rrsw_kg26.o  	     \
                           rrsw_kg27.o             rrsw_kg28.o               \
                           rrsw_kg29.o
	$(F90) $(R8) -c $(OPT) $<


rrtmg_sw_vrtqdr.o        : rrtmg_sw_vrtqdr.F90                               \
                           parkind.o               parrrsw.o
	$(F90) $(R8) -c $(OPT) $<


#EOC
