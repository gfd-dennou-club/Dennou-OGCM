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


###############

bin_PROGRAMS = 

bin_PROGRAMS += globalSWM_DG
globalSWM_DG_SOURCES =\
	LagrangePolyn_mod.f90 \
	DGElement_mod.f90 \
	SimParameters_mod.f90 \
	GridSet_mod.f90 \
	DGHelper_mod.f90 \
	DGCalcusUtil_mod.f90 \
	VariableSet_mod.f90 \
	GovernEquationSolver_mod.f90 \
	OutputData_mod.f90 \
	experiment/W92_case1/Exp_Williamson92_Case1.f90 \
	experiment/W92_case2/Exp_Williamson92_Case2.f90 \
	experiment/W92_case5/Exp_Williamson92_Case5.f90 \
	experiment/WindDrivenCirc/Exp_WindDrivenCirc_mod.f90 \
	globalSWModel_main.f90

bin_PROGRAMS += lagrangeBasisTool 
lagrangeBasisTool_SOURCES =\
	LagrangeBasisTool.f90

#bin_PROGRAMS += rectDomainGridGen
#rectDomainGridGen_SOURCES =\
	RectDomainGridGen.f90

###########

TESTS = $(check_PROGRAMS)

check_PROGRAMS = 

check_PROGRAMS += test_DGBaseMod
test_DGBaseMod_SOURCES =\
	DGBaseMod_test.f90

test_DGBaseMod_LDADD =\
	LagrangePolyn_mod.o \
	DGElement_mod.o \
	$(LDADD)
