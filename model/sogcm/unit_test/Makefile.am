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
        $(MODPATHFLAG)$(srcdir)/../physics_driver \
        $(MODPATHFLAG)$(srcdir)/../../ogcm/util \
        $(MODPATHFLAG)$(srcdir)/../experiment \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../InitCond_mod.o \
	$(srcdir)/../physics_driver/libDSOGCMPhysicsDriver.la \
	$(srcdir)/../experiment/libDSOGCMExp.la \
	$(srcdir)/../dynamics/libDSOGCMDynamics.la \
	$(srcdir)/../setup/libDSOGCMSetup.la \
	$(srcdir)/../util/libDSOGCMUtil.la \
	$(srcdir)/../util/libSpmlUtil.la \
	$(srcdir)/../../../src/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

TESTS = $(check_PROGRAMS)

check_PROGRAMS = 

check_PROGRAMS += test_SpmlUtil_mod
test_SpmlUtil_mod_SOURCES =\
	SpmlUtil_mod_test.f90

#check_PROGRAMS += test_HydroBousEqVDiffProc_mod
#test_HydroBousEqVDiffProc_mod_SOURCES =\
	HydroBousEqVDiffProc_mod_test.f90

check_PROGRAMS += test_SeaWaterProp_mod
test_SeaWaterProp_mod_SOURCES =\
	SeaWaterProp_mod_test.f90

check_PROGRAMS += test_EOS_mod
test_EOS_mod_SOURCES =\
	EOS_mod_test.f90

check_PROGRAMS += test_SGSConvAdjust_mod
test_SGSConvAdjust_mod_SOURCES = \
	SGSConvAdjust_test.f90

check_PROGRAMS += test_SGSSlowConvAdjust_mod
test_SGSSlowConvAdjust_mod_SOURCES = \
	SGSSlowConvAdjust_test.f90

check_PROGRAMS += test_GMScheme_mod
test_GMScheme_mod_SOURCES = \
	GMScheme_test.f90
