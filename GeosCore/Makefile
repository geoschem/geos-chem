# $Id: Makefile,v 1.1 2009/09/16 14:06:47 bmy Exp $
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the GeosCore subdirectory)
#
# !DESCRIPTION: This is the main GEOS-Chem makefile.  It compiles the 
# GEOS-Chem core source code files and bundles all of the object files (*.o) 
# into the libGeosCore.a library (located in the LIB directory).  Module files 
# (*.mod) are copied to the MOD directory.
#
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Most of the time this Makefile will be called automatically    %%%
# %%% from the router Makefile in the top-level directory.  However, if    %%%
# %%% you are in the ./GeosCore directory, then you can call this Makefile %%%
# %%% to build the GEOS-Chem source code, libraries, and executables.      %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# BIN        Specifies the directory where executable files are stored
# CC         Contains the default C compilation commands (for PGI only)
# DOC        Specifies the directory for generating documentation w/ ProTeX
# EXE        Specifies the name of the executable file
# F90        Contains the Fortran compilation commands
# FREEFORM   Contains the command to force F90 "free format" compilation
# HDR        Specifies the directory where include files are found
# LD         Contains the command to link to libraries & make executable
# LIB        Specifies the directory where library files (*.a) are stored
# KPP        Specifies the directory where th KPP solver files reside
# MOD        Specifies the directory where module files (*.mod) are stored
# OBJ        Specifies the list of object files (*.o) to be created.
# R8         Contains the command to force REAL -> REAL*8
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
#                                                                             .
# NOTE: CC, F90, FREEFORM, LD, R8 are included from "Makefile_header.mk".
#
# !REVISION HISTORY: 
#  27 Jul 2009 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL = /bin/sh
BIN   = ../bin
DOC   = ../doc
HDR   = ../Headers
EXE   = geos
LIB   = ../lib
LINK  = -lGeosCore -lKpp -lKppSolver
KPP   = ../KPP
MOD   = ../mod

# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include ../Makefile_header.mk

#=============================================================================
# List of files to compile (the order is important!).  We specify these as
# a list of object files (*.o).  For each object file, the "make" utility
# will find the corresponding source code file (*.f, *.f90) and compile it. 
#=============================================================================


# KPP files parameter (NTRAC)
NTRAC     = 43
KPPSOLVER = rosenbrock
#
.DEFAULT_GOAL := all
#-------------------------------------------------------------------

OBJ =                         \
ifort_errmsg.o                \
fjx_acet_mod.o                \
charpak_mod.o                 \
error_mod.o                   \
logical_mod.o                 \
directory_mod.o               \
unix_cmds_mod.o               \
tracer_mod.o                  \
julday_mod.o                  \
file_mod.o                    \
grid_mod.o                    \
time_mod.o                    \
bpch2_mod.o                   \
regrid_1x1_mod.o              \
pressure_mod.o                \
transfer_mod.o                \
future_emissions_mod.o        \
lai_mod.o                     \
tracerid_mod.o                \
benchmark_mod.o               \
comode_mod.o                  \
gckpp_Rates.o                 \
gckpp_Initialize.o            \
gckpp_Integrator.o            \
gckpp_Model.o	              \
diag_mod.o                    \
dao_mod.o                     \
tropopause_mod.o              \
pbl_mix_mod.o                 \
diag03_mod.o                  \
diag04_mod.o                  \
diag41_mod.o                  \
diag42_mod.o                  \
diag48_mod.o                  \
diag49_mod.o                  \
diag50_mod.o                  \
diag51_mod.o                  \
diag56_mod.o                  \
diag_oh_mod.o                 \
diag_pl_mod.o                 \
ocean_mercury_mod.o           \
drydep_mod.o                  \
vdiff_pre_mod.o               \
vdiff_mod.o                   \
scale_anthro_mod.o            \
edgar_mod.o                   \
bravo_mod.o                   \
emep_mod.o                    \
epa_nei_mod.o                 \
streets_anthro_mod.o          \
arctas_ship_emiss_mod.o       \
icoads_ship_mod.o             \
cac_anthro_mod.o              \
vistas_anthro_mod.o           \
geia_mod.o                    \
global_ch4_mod.o              \
global_hno3_mod.o             \
global_no3_mod.o              \
global_nox_mod.o              \
global_oh_mod.o               \
global_o1d_mod.o              \
global_o3_mod.o               \
uvalbedo_mod.o                \
RnPbBe_mod.o                  \
Kr85_mod.o                    \
acetone_mod.o                 \
aerosol_mod.o                 \
aircraft_nox_mod.o            \
biofuel_mod.o                 \
gc_biomass_mod.o              \
gfed2_biomass_mod.o           \
biomass_mod.o                 \
c2h6_mod.o                    \
ch3i_mod.o                    \
a3_read_mod.o                 \
a6_read_mod.o                 \
i6_read_mod.o                 \
gcap_read_mod.o               \
gwet_read_mod.o               \
xtra_read_mod.o               \
megan_mod.o                   \
carbon_mod.o                  \
lightning_nox_mod.o           \
optdepth_mod.o                \
planeflight_mod.o             \
restart_mod.o                 \
rpmares_mod.o                 \
isoropia_mod.o                \
wetscav_mod.o                 \
seasalt_mod.o                 \
sulfate_mod.o                 \
hcn_ch3cn_mod.o               \
tagged_co_mod.o               \
tagged_ox_mod.o               \
h2_hd_mod.o                   \
gcap_convect_mod.o            \
fvdas_convect_mod.o           \
convection_mod.o              \
pjc_pfix_mod.o                \
pjc_pfix_geos5_window_mod.o   \
dust_dead_mod.o               \
dust_mod.o                    \
co2_mod.o                     \
mercury_mod.o                 \
toms_mod.o                    \
tpcore_bc_mod.o               \
tpcore_fvdas_mod.o            \
tpcore_mod.o                  \
tpcore_window_mod.o           \
tpcore_geos5_window_mod.o     \
transport_mod.o               \
upbdflx_mod.o                 \
chemistry_mod.o               \
emissions_mod.o               \
gamap_mod.o                   \
input_mod.o                   \
CO_strat_pl.o                 \
airmas.o                      \
anthroems.o                   \
arsl1k.o                      \
backsub.o                     \
biofit.o                      \
boxvl.o                       \
calcrate.o                    \
chemdr.o                      \
cleanup.o                     \
decomp.o                      \
diag1.o                       \
diag3.o                       \
diag_2pm.o                    \
diagoh.o                      \
emf_scale.o                   \
emfossil.o                    \
emisop.o                      \
emisop_grass.o                \
emisop_mb.o                   \
emissdr.o                     \
emmonot.o                     \
fertadd.o                     \
findmon.o                     \
fyrno3.o                      \
fyhoro.o                      \
gasconc.o                     \
get_global_ch4.o              \
getifsun.o                    \
initialize.o                  \
jsparse.o                     \
ksparse.o                     \
lump.o                        \
main.o                        \
ndxx_setup.o                  \
ohsave.o                      \
partition.o                   \
pderiv.o                      \
physproc.o                    \
precipfrac.o                  \
pulsing.o                     \
rdisopt.o                     \
rdlai.o                       \
rdland.o                      \
rdlight.o                     \
rdmonot.o                     \
rdsoil.o                      \
readchem.o                    \
reader.o                      \
readlai.o                     \
ruralbox.o                    \
schem.o                       \
setbase.o                     \
setemdep.o                    \
setemis.o                     \
setmodel.o                    \
sfcwindsqr.o                  \
smvgear.o                     \
soilbase.o                    \
soilcrf.o                     \
soilnoxems.o                  \
soiltemp.o                    \
soiltype.o                    \
subfun.o                      \
sunparam.o                    \
tcorr.o                       \
tpcore_fvdas_mod.o            \
tropopause.o                  \
update.o                      \
xltmmp.o                      \
BLKSLV.o                      \
CLDSRF.o                      \
EFOLD.o                       \
FLINT.o                       \
GAUSSP.o                      \
GEN.o                         \
fjx_acet_mod.o                \
JRATET.o                      \
JVALUE.o                      \
LEGND0.o                      \
MATIN4.o                      \
MIESCT.o                      \
NOABS.o                       \
OPMIE.o                       \
RD_TJPL.o                     \
SPHERE.o                      \
XSEC1D.o                      \
XSECO2.o                      \
XSECO3.o                      \
fast_j.o                      \
fjfunc.o                      \
inphot.o                      \
jv_index.o                    \
mmran_16.o                    \
photoj.o                      \
rd_js.o                       \
rd_prof.o                     \
set_aer.o                     \
set_prof.o                    


#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: all clean realclean doc docclean help

all: lib libkpp exe

lib: $(OBJ)
	cp *.mod $(MOD)
	$(AR) crs libGeosCore.a $(OBJ)
	mv libGeosCore.a $(LIB)

libkpp: 
	@$(MAKE) -C $(KPP)

exe: main.o
	$(LD) main.o $(LINK) -o $(EXE)
	cp -f $(EXE) $(BIN)

clean:
	rm -f *.o *.mod geos

realclean:
	$(MAKE) clean
	$(MAKE) -C $(KPP) clean
	$(MAKE) docclean
	rm -f $(LIB)/*.a
	rm -f $(MOD)/*.mod

doc:
	@$(MAKE) -C $(DOC) all

docclean: 
	@$(MAKE) -C $(DOC) clean

help:
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'all       Build all GEOS-Chem libraries & executables'
	@echo 'lib       Only build GEOS-Chem objects & libraries in GeosCore subdir'
	@echo 'libkpp    Only build GEOS-Chem objects & libraries in KPP subdir'
	@echo 'exe       Build GEOS-Chem executable without KPP solver'
	@echo 'clean     Removes *.o, *.mod files in source code subdirs only'
	@echo 'realclean Removes *.o, *mod, *.lib *.a files source code, mod, lib, bin subdirs'
	@echo 'doc       Builds GEOS-Chem documentation (*.ps, *.pdf) with ProTeX in doc subdir'
	@echo 'docclean  Removes *.tex, *.pdf, *,ps from the doc subdir'
	@echo 'help      Displays this help screen'
	@echo ''
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'COMPILER=___   Options: ifort pgi sun xlf (default is ifort)'
	@echo 'DEBUG=yes      Compiles GEOS-Chem for use w/ a debugger'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debugging purposes)'
	@echo 'NTRAC=[43|54]  Specifies # of tracers for KPP chemical solver'
	@echo 'KPPSOLVER=___  Specifies the integrator used w/ KPP:'
	@echo '               Options: lsodes radau5 rosenbrock runge_kutta (default is lsodes)'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Original from Philippe 
#
##=============================================================================
##  Executable and Documentation
##=============================================================================
#.PHONY : all kppintegrator
#
#all: kppintegrator geos
#
#geos: $(MODS) $(OBJS) $(OBJSe) $(FJ)
#	$(F90)  $(MODS) $(OBJS)  $(OBJSe) $(FJ) -o geos
#
##=============================================================================
## Two rules for KPP : update solver only if needed, and set directory
##=============================================================================
#
## (1) update solver (S=Source, T=Target)
#SOLVER_SFILE=./KPP/int/gckpp_Integrator_$(KPPSOLVER).f90
#SOLVER_TFILE=./KPP/$(NTRAC)t/gckpp_Integrator.f90
#
#kppintegrator:
#	@diff -q $(SOLVER_SFILE) $(SOLVER_TFILE) ;\
#	if [ $$? == 1 ] ; then   \
#		echo " copy $(SOLVER_SFILE) --> $(SOLVER_TFILE)";\
#		cp $(SOLVER_SFILE) $(SOLVER_TFILE) ;    \
#	fi
#
## (2) specify files source location
#gckpp_%.o : ./KPP/$(NTRAC)t/gckpp_%.f90
#	$(F90) -c -free $< -o $@
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#==============================================================================
# Dependencies listing
#==============================================================================
tpcore_geos5_window_mod.o      : tpcore_geos5_window_mod.f90
	$(F90) -c $(R8) $*.f90
tpcore_mod.o                   : tpcore_mod.f
	$(F90) -c $(R8) $*.f
tpcore_window_mod.o            : tpcore_window_mod.f
	$(F90) -c $(R8) $*.f
