SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../../src/io \
        $(MODPATHFLAG)$(srcdir)/../seaice \
        $(MODPATHFLAG)$(srcdir)/../util \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDSOGCMSetup.la 
nodist_libDSOGCMSetup_la_SOURCES= \
	OptionParser_mod.f90 \
	Constants_mod.f90 \
	GridSet_mod.f90 \
	TemporalIntegSet_mod.f90 \
	BoundCondSet_mod.f90 \
	VariableSet_mod.f90 \
	GovernEqSet_mod.f90


