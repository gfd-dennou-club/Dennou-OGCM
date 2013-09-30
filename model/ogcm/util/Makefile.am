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
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../../src/io/libDQGModelIO.la \
	$(srcdir)/../../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDOGCMUtil.la 
nodist_libDOGCMUtil_la_SOURCES= \
	CGridFieldDataUtil_mod.f90 \
	VGridFieldDataUtil_mod.f90 \
	EqState_JM95_mod.f90 \
	EOSDriver_mod.f90



