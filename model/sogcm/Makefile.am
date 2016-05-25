SUBDIRS = \
	util setup seaice boundary_cond  io                      \
	dynamics sgs_param physics_driver governeq_driver        \
	experiment                                     \
	unit_test tool

################

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../include/ \
        $(MODPATHFLAG)$(srcdir)/setup \
        $(MODPATHFLAG)$(srcdir)/util \
        $(MODPATHFLAG)$(srcdir)/boundary_cond \
	$(MODPATHFLAG)$(srcdir)/dynamics \
        $(MODPATHFLAG)$(srcdir)/governeq_driver \
        $(MODPATHFLAG)$(srcdir)/physics_driver \
        $(MODPATHFLAG)$(srcdir)/seaice \
        $(MODPATHFLAG)$(srcdir)/io \
        $(MODPATHFLAG)$(srcdir)/experiment \
	$(FCFLAGS) $(SYSFFLAGS)


LDADD = \
	$(srcdir)/../../libsrc/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../../libsrc/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

######################

CLEANFILES = *$(MODEXT)
all: all-am copy_header

copy_header:
		@echo "Copy generated header files into $(srcdir)/include"
		@$(MKDIR_P) $(srcdir)/include
		@find . -name "*@MODEXT@" -exec $(install_sh_DATA) {} $(srcdir)/include \;

bin_PROGRAMS = 

#####################

if DSOGCM_3D

bin_PROGRAMS += dsogcm3d
dsogcm3d_SOURCES =\
	ogcm_main.f90

dsogcm3d_LDADD = \
	libDSOGCM.la \
	$(LDADD)

lib_LTLIBRARIES = libDSOGCM.la
libDSOGCM_la_SOURCES = DSOGCM_main_mod.f90
libDSOGCM_la_LIBADD = \
	$(srcdir)/experiment/libDSOGCMExp.la \
	$(srcdir)/governeq_driver/libDSOGCMGovernEqDriver.la \
	$(srcdir)/physics_driver/libDSOGCMPhysicsDriver.la \
	$(srcdir)/dynamics/libDSOGCMDynamics.la \
	$(srcdir)/io/libDSOGCMIO.la \
	$(srcdir)/boundary_cond/libDSOGCMBC.la \
	$(srcdir)/seaice/libDSOGCMSeaIce.la \
	$(srcdir)/setup/libDSOGCMSetup.la \
	$(srcdir)/util/libDSOGCMUtil.la \
	$(srcdir)/util/libSpmlUtil.la

endif

#####################

if DSOGCM_AXISYM

bin_PROGRAMS += dsogcm_axisym

dsogcm_axisym_SOURCES =\
	ogcm_main.f90
dsogcm_axisym_LDADD = \
	libDSOGCM.la \
	$(LDADD)

lib_LTLIBRARIES = libDSOGCM.la
libDSOGCM_la_SOURCES = DSOGCM_main_mod.f90
libDSOGCM_la_LIBADD = \
	$(srcdir)/experiment/libDSOGCMExp.la \
	$(srcdir)/governeq_driver/libDSOGCMGovernEqDriver.la \
	$(srcdir)/physics_driver/libDSOGCMPhysicsDriver.la \
	$(srcdir)/dynamics/libDSOGCMDynamics.la \
	$(srcdir)/io/libDSOGCMIO.la \
	$(srcdir)/boundary_cond/libDSOGCMBC.la \
	$(srcdir)/seaice/libDSOGCMSeaIce.la \
	$(srcdir)/setup/libDSOGCMSetup.la \
	$(srcdir)/util/libDSOGCMUtil.la \
	$(srcdir)/util/libSpmlUtil.la

endif

#############################

