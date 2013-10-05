SUBDIRS = setup util dynamics

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../src/io \
        $(MODPATHFLAG)$(srcdir)/setup \
        $(MODPATHFLAG)$(srcdir)/util \
        $(MODPATHFLAG)$(srcdir)/dynamics \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/dynamics/libDOGCMDynamics.la \
	$(srcdir)/util/libDOGCMUtil.la \
	$(srcdir)/setup/libDOGCMSetup.la \
	$(srcdir)/../../src/io/libDQGModelIO.la \
	$(srcdir)/../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

bin_PROGRAMS += dogcm 
dogcm_SOURCES =\
	GovernEqSolverDriver_mod.f90 \
	ogcm_main.f90

