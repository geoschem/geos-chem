#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 17Oct2008  da Silva  Changed from distributed to centralized makefile.
#
#-------------------------------------------------------------------------

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # GMAO stuff
#  include ../Makefile_header.mk

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

THIS := $(shell basename `pwd`)
NAME := $(patsubst %_GridComp,%,$(THIS))
LIB_THIS := lib$(THIS).a

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRC_DIRS := .

INC_DIRS := . $(INC_MAPL_BASE)  $(INC_MFHDF3) ../Headers ../mod

MOD_DIRS = . $(INC_ESMF) $(INC_MPEU) $(INC_MPI) \
             $(INC_DIRS) $(INC_GEOS_SHARED) $(ESMAINC)/MAPL_Base \
             $(INC_CHEM_BASE) $(INC_CHEM_SHARED)

SRCS := $(foreach dir,$(SRC_DIRS), \
        $(wildcard $(dir)/*.[fFc]) $(wildcard $(dir)/*.[fF]90) )

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
ACGS := $(NAME)_ExportSpec___.h $(NAME)_GetPointer___.h $(NAME)_History___.rc
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

FREAL = $(FREAL4)

ROOTDIR = ..
HDR     = $(ROOTDIR)/Headers
MOD     = $(ROOTDIR)/mod
LIB     = $(ROOTDIR)/lib

ifeq ($(ESMF),yes)
USER_FFLAGS += -DESMF_
endif

ifeq ($(ESMF_TESTBED),yes)
USER_FFLAGS += -DESMF_TESTBED_
endif

###############################################################################
# Now add USER_FFLAGS setting from env variable (bmy, 4/20/10)
#USER_FFLAGS = $(BIG_ENDIAN) 
##############################################################################
USER_FFLAGS += $(BIG_ENDIAN) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

ACG_FLAGS += -F

vpath % /usr/include $(SRC_DIRS) $(INC_DIRS) $(MOD_DIRS)

$(LIB_THIS) lib : $(ACGS) $(DEPS) $(OBJS)

$(LIB_THIS) lib : $(	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS) 
$(LIB_THIS) esmf : $(ACGS) $(DEPS) $(OBJS)

$(LIB_THIS) esmf : $(	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS)
	mv $(LIB_THIS) $(LIB)

$(ACGS) : $(NAME)_Registry.rc $(ACG)
	@$(ACG) $(ACG_FLAGS) $(NAME)_Registry.rc

%.x : $(LIB) %.o
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB_THIS) \
	      $(LIB_CHEM_BASE) $(LIB_CHEM_SHARED) $(LIB_PILGRIM) \
              $(LIB_MAPL_BASE) $(LIB_CFIO) $(LIB_GFIO) $(LIB_MPEU) \
              $(LIB_ESMF) $(LIB_SDF) \
              $(LIB_SYS) $(LIB_MPI) $(ESMF_LDFLAGS)

test_rc: $(LIB)
	$(CP) $(ESMAETC)/GEOSCHEMchem_GridComp.rc     ./Tests

clean_rc:
	-$(RM) -f Tests/ExtData Tests/*_Registry.rc Tests/*_GridComp.rc 

clean realclean:
	-$(RM) -f *.[hoad]
	-$(RM) -f *.mod

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

# For parallel install
# --------------------
  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

# Dependencies
gc_esmf_drv.o   : gc_esmf_drv.F90 		gc_esmf_comp.o  \
		  gc_esmf_testbed_value_mod.o

gc_esmf_comp.o  : gc_esmf_comp.F90 		gc_esmf_chunkmod.o

gc_esmf_chunkmod.o : gc_esmf_chunkmod.F90

gc_esmf_testbed_value_mod.o : gc_esmf_testbed_value_mod.F90

gc_esmf_testbed.o : gc_esmf_testbed.F90     	gc_esmf_comp.o \
		    gc_esmf_testbed_value_mod.o	gc_esmf_drv.o

#.
