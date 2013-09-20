SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../../src/io \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../util \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../../src/io/libDQGModelIO.la \
	$(srcdir)/../../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDOGCMDynamics.la 
nodist_libDOGCMDynamics_la_SOURCES= \
	HydroBouEqSolver_mod.f90


