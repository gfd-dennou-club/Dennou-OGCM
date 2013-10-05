SUBDIRS = setup

# util dynamics

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../src/io \
        $(MODPATHFLAG)$(srcdir)/../ogcm/util \
        $(MODPATHFLAG)$(srcdir)/setup \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/setup/libDSOGCMSetup.la \
	$(srcdir)/../ogcm/util/libDOGCMUtil.la \
	$(srcdir)/../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

bin_PROGRAMS += dsogcm 
dsogcm_SOURCES =\
	ogcm_main.f90

