SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../../src/io \
        $(MODPATHFLAG)$(srcdir)/../ \
        $(MODPATHFLAG)$(srcdir)/../experiment \
        $(MODPATHFLAG)$(srcdir)/../dynamics \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../util \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../experiment/libDSOGCMExp.la \
	$(srcdir)/../dynamics/libDSOGCMDynamics.la \
	$(srcdir)/../setup/libDSOGCMSetup.la \
	$(srcdir)/../util/libDSOGCMUtil.la \
	$(srcdir)/../util/libSpmlUtil.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS=diagVar
diagVar_SOURCES =\
	diagvar/DiagVarSet_mod.f90 \
	diagvar/DiagVarEval_mod.f90 \
	diagvar/DiagVarFileSet_mod.f90 \
	diagvar/BudgetAnalysis_mod.f90 \
	diagvar/SpectralAnalysis_mod.f90 \
	diagvar/diagVar_main.f90

