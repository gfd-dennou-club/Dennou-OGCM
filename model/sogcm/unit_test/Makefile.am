SUBDIRS =

test: check

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
        $(MODPATHFLAG)$(srcdir)/../dynamics \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../dynamics/libDSOGCMDynamics.la \
	$(srcdir)/../setup/libDSOGCMSetup.la \
	$(srcdir)/../util/libDSOGCMUtil.la \
	$(srcdir)/../../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

TESTS = $(check_PROGRAMS)

check_PROGRAMS = 

check_PROGRAMS += test_SpmlUtil_mod
test_SpmlUtil_mod_SOURCES =\
	SpmlUtil_mod_test.f90

check_PROGRAMS += test_HydroBousEqVDiffProc_mod
test_HydroBousEqVDiffProc_mod_SOURCES =\
	HydroBousEqVDiffProc_mod_test.f90

