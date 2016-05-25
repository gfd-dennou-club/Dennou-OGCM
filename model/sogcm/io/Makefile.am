SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../include/ \
        $(MODPATHFLAG)$(srcdir)/../util \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../seaice \
	$(MODPATHFLAG)$(srcdir)/../boundary_cond \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDSOGCMIO.la 
dist_libDSOGCMIO_la_SOURCES= \
	DataFileSet_mod.f90 \
	RestartDataFileSet_mod.f90

