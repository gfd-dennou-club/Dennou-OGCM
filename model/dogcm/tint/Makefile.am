SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../libsrc/ \
        $(MODPATHFLAG)$(srcdir)/../util      \
        $(MODPATHFLAG)$(srcdir)/../admin     \
        $(MODPATHFLAG)$(srcdir)/../io        \
	$(MODPATHFLAG)$(srcdir)/../boundary  \
	$(MODPATHFLAG)$(srcdir)/../physics   \
        $(MODPATHFLAG)$(srcdir)/../dynamics  \
        $(MODPATHFLAG)$(srcdir)/../seaice    \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDOGCMTInt.la 
dist_libDOGCMTInt_la_SOURCES= \
	DOGCM_TInt_common_mod.f90 \
	DOGCM_TInt_driver_mod.f90

