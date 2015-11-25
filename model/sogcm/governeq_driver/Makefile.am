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
        $(MODPATHFLAG)$(srcdir)/../seaice \
        $(MODPATHFLAG)$(srcdir)/../boundary_cond \
        $(MODPATHFLAG)$(srcdir)/../io \
        $(MODPATHFLAG)$(srcdir)/../dynamics \
        $(MODPATHFLAG)$(srcdir)/../sgs_param \
	$(MODPATHFLAG)$(srcdir)/../physics_driver \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

noinst_LTLIBRARIES= libDSOGCMGovernEqDriver.la 
#nodist_libDSOGCMGovernEqDriver_la_SOURCES= \
	HydroBouEq_TimeInteg_old_mod.f90    \
	HydroBouEq_TimeInteg_v2_mod.f90 \
	SeaiceEq_TimeInteg_mod.f90      \
	GovernEqSolverDriver_mod.f90

nodist_libDSOGCMGovernEqDriver_la_SOURCES= \
	HydroBouEq_TimeInteg_v2_mod.f90 \
	HydroBouEqVorDivForm_TimeInteg_mod.f90 \
	SeaiceEq_TimeInteg_mod.f90      \
	GovernEqSolverDriver_mod.f90

