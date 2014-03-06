SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../../src/io \
        $(MODPATHFLAG)$(srcdir)/../util \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../../src/libDQGModelCore.la \
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
	DataFileSet_mod.f90 \
	RestartDataFileSet_mod.f90 \
	GovernEqSet_mod.f90


GridSet_mod.f90: \
	../util/SpmlUtil_mod.o

DataFileSet_mod.f90: \
	../util/EOSDriver_mod.o