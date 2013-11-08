SUBDIRS = util setup dynamics unit_test

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
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
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
	experiment/Exp_W94_Case2_mod.f90 \
	experiment/Exp_BarotRossbyWave_mod.f90 \
	experiment/Exp_InternalGravWave_mod.f90 \
	experiment/Exp_WindDrivenCirculation_mod.f90 \
	InitCond_mod.f90 \
	ogcm_main.f90

