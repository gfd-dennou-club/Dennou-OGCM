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
	$(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../boudary_cond \
        $(MODPATHFLAG)$(srcdir)/../dynamics \
        $(MODPATHFLAG)$(srcdir)/../io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

noinst_LTLIBRARIES= \
	libDSOGCMSGSParamMesoEddy.la \
	libDSOGCMSGSParamConvection.la \
	libDSOGCMSGSParamVMixing.la

dist_libDSOGCMSGSParamMesoEddy_la_SOURCES= \
	mesoeddy/SGSEddyMixingHelper_mod.f90 \
	mesoeddy/SGSEddyMixing_mod.f90

dist_libDSOGCMSGSParamConvection_la_SOURCES= \
	convection/SGSConvAdjust_mod.f90 \
	convection/SGSSlowConvAdjust_mod.f90

dist_libDSOGCMSGSParamVMixing_la_SOURCES= \
	vmixing/SGSVMixingConst_mod.f90


