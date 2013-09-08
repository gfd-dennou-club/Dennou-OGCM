subdir =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../src/ \
        $(MODPATHFLAG)$(srcdir)/../src/math \
        $(MODPATHFLAG)$(srcdir)/../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../src/data \
        $(MODPATHFLAG)$(srcdir)/../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../src/io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../src/io/libDQGModelIO.la \
	$(srcdir)/../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

bin_PROGRAMS += scvgridGen
scvgridGen_SOURCES =\
	gridGen/gridGen.f90

bin_PROGRAMS += globalSWM 
globalSWM_SOURCES =\
	globalSWModel/SimParameters_mod.f90 \
	globalSWModel/GridSet_mod.f90 \
	globalSWModel/VariableSet_mod.f90 \
	globalSWModel/GovernEquationSolver_mod.f90 \
	globalSWModel/Exp_Williamson94_Case1.f90 \
	globalSWModel/Exp_Williamson94_Case2.f90 \
	globalSWModel/globalSWModel_main.f90

