SUBDIRS = \
	util setup dynamics experiment governeq_driver \
	unit_test tool

# util dynamics

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../src/io \
        $(MODPATHFLAG)$(srcdir)/setup \
        $(MODPATHFLAG)$(srcdir)/util \
        $(MODPATHFLAG)$(srcdir)/dynamics \
        $(MODPATHFLAG)$(srcdir)/governeq_driver \
        $(MODPATHFLAG)$(srcdir)/experiment \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/governeq_driver/libDSOGCMGovernEqDriver.la \
	$(srcdir)/experiment/libDSOGCMExp.la \
	$(srcdir)/dynamics/libDSOGCMDynamics.la \
	$(srcdir)/setup/libDSOGCMSetup.la \
	$(srcdir)/util/libDSOGCMUtil.la \
	$(srcdir)/util/libSpmlUtil.la \
	$(srcdir)/../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

#####################
bin_PROGRAMS += dsogcm 
dsogcm_SOURCES =\
	InitCond_mod.f90 \
	ogcm_main.f90


#############################

