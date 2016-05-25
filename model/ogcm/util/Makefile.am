SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../include \
        $(MODPATHFLAG)$(srcdir)/../setup \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDOGCMUtil.la 
nodist_libDOGCMUtil_la_SOURCES= \
	CGridFieldDataUtil_mod.f90 \
	VGridFieldDataUtil_mod.f90 \
	EqState_JM95_mod.f90 \
	EOSDriver_mod.f90


CGridFieldDataUtil_mod.f90:
	../setup/GridSet_mod.o

