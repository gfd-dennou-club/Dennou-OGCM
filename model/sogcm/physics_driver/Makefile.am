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
        $(MODPATHFLAG)$(srcdir)/../sgs_param \
        $(MODPATHFLAG)$(srcdir)/../io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

noinst_LTLIBRARIES= \
	libDSOGCMPhysicsDriver.la

dist_libDSOGCMPhysicsDriver_la_SOURCES= \
	PhysicsDriver_mod.f90

libDSOGCMPhysicsDriver_la_LIBADD= \
	../sgs_param/libDSOGCMSGSParamConvection.la \
	../sgs_param/libDSOGCMSGSParamMesoEddy.la


