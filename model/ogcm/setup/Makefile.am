SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../include \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
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

