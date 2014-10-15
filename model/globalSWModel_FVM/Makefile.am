SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../src/io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../src/io/libDQGModelIO.la \
	$(srcdir)/../../src/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

bin_PROGRAMS += globalSWM 
globalSWM_SOURCES =\
	SimParameters_mod.f90 \
	GridSet_mod.f90 \
	VariableSet_mod.f90 \
	GovernEquationSolver_mod.f90 \
	OutputData_mod.f90 \
	Exp_Williamson94_Case1.f90 \
	Exp_Williamson94_Case2.f90 \
	Exp_Williamson94_Case5.f90 \
	globalSWModel_main.f90

