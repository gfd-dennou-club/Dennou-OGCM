SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../include/ \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../util \
        $(MODPATHFLAG)$(srcdir)/../boundary_cond \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

noinst_LTLIBRARIES= libDSOGCMDynamics.la 
dist_libDSOGCMDynamics_la_SOURCES= \
	HydroBouEqSolverRHS_old_mod.f90 \
	HydroBouEqSolverRHS_v2_mod.f90 \
	HydroBouEqRHS_DynProc_mod.f90 \
	HydroBouEqSolverVImplProc_mod.f90


