SUBDIRS = \
	util setup dynamics experiment physics_driver governeq_driver \
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
        $(MODPATHFLAG)$(srcdir)/physics_driver \
        $(MODPATHFLAG)$(srcdir)/experiment \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../src/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

#####################

if DSOGCM_3D

bin_PROGRAMS += dsogcm3d
dsogcm3d_SOURCES =\
	InitCond_mod.f90 \
	ogcm_main.f90
dsogcm3d_LDADD = \
	$(srcdir)/governeq_driver/libDSOGCMGovernEqDriver.la \
	$(srcdir)/physics_driver/libDSOGCMPhysicsDriver.la \
	$(srcdir)/experiment/libDSOGCMExp.la \
	$(srcdir)/dynamics/libDSOGCMDynamics.la \
	$(srcdir)/setup/libDSOGCMSetup.la \
	$(srcdir)/util/libDSOGCMUtil.la \
	$(srcdir)/util/libSpmlUtil.la \
	$(LDADD)

endif

#####################

if DSOGCM_AXISYM

bin_PROGRAMS += dsogcm_axisym
dsogcm_axisym_SOURCES =\
	InitCond_mod.f90 \
	ogcm_main.f90
dsogcm_axisym_LDADD = \
	$(srcdir)/governeq_driver/libDSOGCMGovernEqDriver.la \
	$(srcdir)/physics_driver/libDSOGCMPhysicsDriver.la \
	$(srcdir)/experiment/libDSOGCMExp.la \
	$(srcdir)/dynamics/libDSOGCMDynamics.la \
	$(srcdir)/setup/libDSOGCMSetup.la \
	$(srcdir)/util/libDSOGCMUtil.la \
	$(srcdir)/util/libSpmlUtil.la \
	$(LDADD)
endif

#############################

