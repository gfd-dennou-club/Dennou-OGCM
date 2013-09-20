SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../../src/io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../../src/io/libDQGModelIO.la \
	$(srcdir)/../../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDOGCMSetup.la 
nodist_libDOGCMSetup_la_SOURCES= \
	OptionParser_mod.f90 \
	Constants_mod.f90 \
	GridSet_mod.f90 \
	TemporalIntegSet_mod.f90 \
	VariableSet_mod.f90 \
	DataFileSet_mod.f90

