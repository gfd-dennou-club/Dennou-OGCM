SUBDIRS = util setup dynamics experiment unit_test

# util dynamics

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../src/io \
        $(MODPATHFLAG)$(srcdir)/../ogcm/util \
        $(MODPATHFLAG)$(srcdir)/setup \
        $(MODPATHFLAG)$(srcdir)/util \
        $(MODPATHFLAG)$(srcdir)/dynamics \
        $(MODPATHFLAG)$(srcdir)/experiment \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/experiment/libDSOGCMExp.la \
	$(srcdir)/dynamics/libDSOGCMDynamics.la \
	$(srcdir)/setup/libDSOGCMSetup.la \
	$(srcdir)/util/libDSOGCMUtil.la \
	$(srcdir)/../ogcm/util/libDOGCMUtil.la \
	$(srcdir)/../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

bin_PROGRAMS += dsogcm 
dsogcm_SOURCES =\
	GovernEqSolverDriver_mod.f90 \
	InitCond_mod.f90 \
	ogcm_main.f90

bin_PROGRAMS += diagVar
diagVar_SOURCES =\
	tool/DiagVarSet_mod.f90 \
	tool/DiagVarEval_mod.f90 \
	tool/DiagVarFileSet_mod.f90 \
	tool/diagVar_main.f90