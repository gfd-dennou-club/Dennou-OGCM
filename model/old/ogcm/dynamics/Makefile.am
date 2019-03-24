SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../include/ \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../util \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = TimeFilterParams

TimeFilterParams_SOURCES= \
	TimeFilterParams.f90

noinst_LTLIBRARIES= libDOGCMDynamics.la 
nodist_libDOGCMDynamics_la_SOURCES= \
	HydroBouEqSolverRHS_mod.f90 \
	BarotModeTimeFilter_mod.f90 \
	HydroBouEqSolver_mod.f90 \
	HydroBouEqSolverSelfStart_mod.f90


