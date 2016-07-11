SUBDIRS = \
	util io admin boundary seaice dynamics physics tint

#util setup seaice boundary_cond  io                      \
	dynamics sgs_param physics_driver governeq_driver        \
	experiment                                     \
	unit_test tool

################

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../include/ \
        $(MODPATHFLAG)$(srcdir)/admin          \
        $(MODPATHFLAG)$(srcdir)/util           \
        $(MODPATHFLAG)$(srcdir)/io             \
	$(MODPATHFLAG)$(srcdir)/boundary       \
	$(MODPATHFLAG)$(srcdir)/dynamics       \
        $(MODPATHFLAG)$(srcdir)/physics        \
        $(MODPATHFLAG)$(srcdir)/seaice         \
        $(MODPATHFLAG)$(srcdir)/tint           \
        $(MODPATHFLAG)$(srcdir)/exp            \
	$(FCFLAGS) $(SYSFFLAGS)


LDADD = \
	$(srcdir)/../../libsrc/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../../libsrc/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

#####################################################################

CLEANFILES = *$(MODEXT)
all: all-am copy_header

copy_header:
		@echo "Copy generated header files into $(srcdir)/include"
		@$(MKDIR_P) $(srcdir)/include
		@find . -name "*@MODEXT@" -exec $(install_sh_DATA) {} $(srcdir)/include \;

bin_PROGRAMS = 

######################################################################

#= Barotropic Rossby wave -
#exp_name = BarotRossby

#= Internal gravity wave --
#exp_name = IGW

#= WindDrivCirc -----------
#exp_name = WindDrivCirc

#= APEOGCirc --------------
#exp_name = APEOGCirc

#= APEOGCircSIce ----------
#exp_name = APEOGCircSIce

#= APECoupleClimate -------
exp_name = APECoupleClimate

#= EqJetAccel -------------
#exp_name = EqJetAccel

#######################################################################

if DSOGCM_3D

bin_PROGRAMS += dogcm3d
dogcm3d_SOURCES =\
	ogcm_main.f90
dogcm3d_LDADD = \
	libDOGCM.la \
	$(srcdir)/seaice/libDSIce.la              \
	$(LDADD)
endif

#######################################################################

if DSOGCM_AXISYM

bin_PROGRAMS += dogcm_axisym
dogcm_axisym_SOURCES =\
	ogcm_main.f90
dogcm_axisym_LDADD = \
	$(srcdir)/seaice/libDSIce.la              \
	libDOGCM.la \
	$(LDADD)
endif

#######################################################################

lib_LTLIBRARIES = libDOGCM.la
libDOGCM_la_SOURCES = \
	exp/$(exp_name)/DOGCM_Exp_$(exp_name)_mod.f90 \
	exp/DOGCM_Exp_driver_mod.f90                  \
	DOGCM_main_mod.f90

DOGCM_Exp_driver_mod.f90 : exp/$(exp_name)/DOGCM_Exp_$(exp_name)_mod.f90

libDOGCM_la_LIBADD = \
	$(srcdir)/tint/libDOGCMTInt.la            \
	$(srcdir)/physics/libDOGCMPhysics.la      \
	$(srcdir)/dynamics/libDOGCMDynamics.la    \
	$(srcdir)/io/libDOGCMIO.la                \
	$(srcdir)/boundary/libDOGCMBoundary.la    \
	$(srcdir)/admin/libDOGCMAdmin.la          \
	$(srcdir)/util/libDOGCMUtil.la            \
	$(srcdir)/util/libSpmlUtil.la


