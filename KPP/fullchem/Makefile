#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User: Set here the F90 compiler and options
#       Pedefined compilers: INTEL, PGF, HPUX, LAHEY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#COMPILER = G95
#COMPILER = LAHEY
#COMPILER = INTEL
#COMPILER = PGF
#COMPILER = HPUX
COMPILER = GFORTRAN

FC_G95     = g95
FOPT_G95   = -cpp -O -pg -fbounds-check -fimplicit-none  -Wall -ftrace=full

FC_LAHEY   = lf95
# More aggressive for production runs:
#FOPT_LAHEY = -Cpp --pca -O
# More checking for debugging:
FOPT_LAHEY = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap --trace --chkglobal

FC_INTEL   = ifort 
# More aggressive for production runs:
#FOPT_INTEL = -cpp -O -fp-model precise -pc80 -prec_div
# More checking for debugging:
FOPT_INTEL = -cpp -O0 -fp-model strict -implicitnone -ftrapuv \
              -debug all -check all -warn all

FC_PGF     = pgf90
# More aggressive for production runs:
FOPT_PGF   = -Mpreprocess -O -fast -pc 80 -Kieee
# More checking for debugging:
#FOPT_PGF   = -Mpreprocess -O0 -Mbounds -Mchkfpstk -Mchkptr -Mchkstk \
#             -Ktrap=fp -pc 80 -Kieee

FC_HPUX    = f90
FOPT_HPUX  = -O -u +Oall +check=on

FC_GFORTRAN     = gfortran
#FOPT_GFORTRAN   = -cpp -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow #bounds
FOPT_GFORTRAN   = -cpp -O 

# define FULL_ALGEBRA for non-sparse integration
FC   = $(FC_$(COMPILER))
FOPT = $(FOPT_$(COMPILER)) # -DFULL_ALGEBRA

LIBS =
#LIBS = -llapack -lblas

# Command to create Matlab mex gateway routines 
# Note: use $(FC) as the mex Fortran compiler
MEX  = mex

GENSRC = gckpp_Precision.F90  \
	 gckpp_Parameters.F90     \
	 gckpp_Global.F90  

GENOBJ = gckpp_Precision.o    \
	 gckpp_Parameters.o       \
	 gckpp_Global.o     

FUNSRC = gckpp_Function.F90
FUNOBJ = gckpp_Function.o

JACSRC = gckpp_JacobianSP.F90  gckpp_Jacobian.F90
JACOBJ = gckpp_JacobianSP.o    gckpp_Jacobian.o

UTLSRC = gckpp_Rates.F90 gckpp_Util.F90 gckpp_Monitor.F90 fullchem_RateLawFuncs.F90 rateLawUtilFuncs.F90
UTLOBJ = gckpp_Rates.o   gckpp_Util.o   gckpp_Monitor.o fullchem_RateLawFuncs.o rateLawUtilFuncs.o

LASRC  = gckpp_LinearAlgebra.F90 
LAOBJ  = gckpp_LinearAlgebra.o   

STOCHSRC = gckpp_Stochastic.F90 
STOCHOBJ = gckpp_Stochastic.o 

MODSRC   = gckpp_Model.F90
MODOBJ   = gckpp_Model.o

INISRC   = gckpp_Initialize.F90
INIOBJ 	 = gckpp_Initialize.o

MAINSRC = kpp_standalone.F90   gckpp_Initialize.F90   gckpp_Integrator.F90 gckpp_Model.F90 
MAINOBJ = kpp_standalone.o     gckpp_Initialize.o     gckpp_Integrator.o


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User: modify the line below to include only the
#       objects needed by your application
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ALLOBJ = $(GENOBJ) $(JACOBJ) $(FUNOBJ)  $(HESOBJ) $(STMOBJ) \
	 $(UTLOBJ) $(LAOBJ) $(MODOBJ) $(INIOBJ) $(SFCOBJ) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User: modify the line below to include only the
#       executables needed by your application
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all:    exe

exe:	$(ALLOBJ) $(MAINOBJ) kpp_standalone_init.o
	$(FC) $(FOPT) kpp_standalone.F90 gckpp_Integrator.o kpp_standalone_init.o $(ALLOBJ) $(LIBS) -o kpp_standalone.exe


stochastic:$(ALLOBJ) $(STOCHOBJ) $(MAINOBJ)
	$(FC) $(FOPT) $(ALLOBJ) $(STOCHOBJ) $(MAINOBJ) $(LIBS) \
	-o gckpp_stochastic.exe

mex:    $(ALLOBJ)
	$(MEX) FC#$(FC) -fortran -O gckpp_mex_Fun.F90     $(ALLOBJ)
	$(MEX) FC#$(FC) -fortran -O gckpp_mex_Jac_SP.F90  $(ALLOBJ)
	$(MEX) FC#$(FC) -fortran -O gckpp_mex_Hessian.F90 $(ALLOBJ)

clean:
	rm -f *.o *.mod\
	gckpp*.dat kpp_standalone.exe gckpp*.mexglx \
	gckpp.map

distclean:
	rm -f *.o *.mod \
	gckpp*.dat kpp_standalone.exe gckpp.map \
	gckpp*.F90 gckpp_*.mexglx

gckpp_Precision.o: gckpp_Precision.F90 
	$(FC) $(FOPT) -c $<

gckpp_Parameters.o: gckpp_Parameters.F90 \
	            gckpp_Precision.o
	$(FC) $(FOPT) -c $<

gckpp_Monitor.o: gckpp_Monitor.F90 \
	             gckpp_Precision.o
	$(FC) $(FOPT) -c $<

gckpp_Global.o: gckpp_Global.F90 \
	            gckpp_Parameters.o gckpp_Precision.o
	$(FC) $(FOPT) -c $<

gckpp_Initialize.o: gckpp_Initialize.F90  $(GENOBJ) 
	$(FC) $(FOPT) -c $<

gckpp_Function.o: gckpp_Function.F90 $(GENOBJ) 
	$(FC) $(FOPT) -c $<

gckpp_Stochastic.o: gckpp_Stochastic.F90  $(GENOBJ) 
	$(FC) $(FOPT) -c $<

gckpp_JacobianSP.o: gckpp_JacobianSP.F90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

gckpp_Jacobian.o: gckpp_Jacobian.F90  $(GENOBJ) gckpp_JacobianSP.o
	$(FC) $(FOPT) -c $<

gckpp_LinearAlgebra.o: gckpp_LinearAlgebra.F90 $(GENOBJ) gckpp_JacobianSP.o
	$(FC) $(FOPT) -c $<

rateLawUtilFuncs.o: rateLawUtilFuncs.F90
	$(FC) $(FOPT) -c $<

fullchem_RateLawFuncs.o: fullchem_RateLawFuncs.F90 rateLawUtilFuncs.o
	$(FC) $(FOPT) -c $<

gckpp_Rates.o: gckpp_Rates.F90  $(GENOBJ) fullchem_RateLawFuncs.o
	$(FC) $(FOPT) -c $<

gckpp_HessianSP.o: gckpp_HessianSP.F90  $(GENOBJ)
	$(FC) $(FOPT) -c $<

gckpp_Hessian.o:  gckpp_Hessian.F90 $(GENOBJ) gckpp_HessianSP.o
	$(FC) $(FOPT) -c $<

gckpp_Util.o: gckpp_Util.F90  $(GENOBJ) gckpp_Monitor.o
	$(FC) $(FOPT) -c $<

gckpp_Main.o: gckpp_Main.F90  $(ALLOBJ) gckpp_Initialize.o gckpp_Model.o gckpp_Integrator.o
	$(FC) $(FOPT) -c $<

gckpp_Model.o: gckpp_Model.F90  $(ALLOBJ) gckpp_Integrator.o
	$(FC) $(FOPT) -c $<

gckpp_Integrator.o: gckpp_Integrator.F90  $(ALLOBJ)
	$(FC) $(FOPT) -c $<

kpp_standalone_init.o: kpp_standalone_init.F90 gckpp_Parameters.o
	$(FC) $(FOPT) -c $<

kpp_standalone.o: kpp_standalone.F90 kpp_standalone_init.o gckpp_Integrator.o $(ALLOBJ)
	$(FC) $(FOPT) -c $<




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
