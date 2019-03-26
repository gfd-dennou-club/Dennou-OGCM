SUBDIRS = src

################

bin_PROGRAMS = 
lib_LTLIBRARIES=
noinst_LTLIBRARIES=

AM_FCFLAGS = \
  $(MODPATHFLAG)$(srcdir)/../../include/ \
  $(MODPATHFLAG)$(srcdir)/include        \
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

######################################################################

#= Barotropic Rossby wave -
#exp_name = BarotRossby

#= Internal gravity wave --
#exp_name = IGW

#= WindDrivCirc -----------
#exp_name = WindDrivCirc

#= APEOGCirc --------------
exp_name = APEOGCirc

#= APEOGCircSIce ----------
#exp_name = APEOGCircSIce

#= APECoupleClimate -------
#exp_name = APECoupleClimate

#= EqJetAccel -------------
#exp_name = EqJetAccel

#= APEOGCircI98BC -------
#exp_name = APEOGCircI98BC


#######################################################################

if DOGCM_3D

bin_PROGRAMS += dogcm3d
dogcm3d_SOURCES =\
	ogcm_main.f90
dogcm3d_LDADD = \
	libDOGCM.la \
	$(srcdir)/seaice/libDSIce.la              \
	$(LDADD)
endif

if DOGCM_AXISYM

bin_PROGRAMS += dogcm_axisym
dogcm_axisym_SOURCES = src/ogcm_main.f90
dogcm_axisym_LDADD =\
	libDOGCM.la                                   \
	$(srcdir)/src/seaice/libDSIce.la              \
	$(LDADD)
endif

src/ogcm_main.f90: src/DOGCM_main_mod.o

#######################################################################

bin_PROGRAMS += ocndiag
ocndiag_SOURCES =\
	tool/ocndiag/OcnDiag_DiagUtil_mod.f90            \
	tool/ocndiag/OcnDiag_HeatTransport_mod.f90       \
	tool/ocndiag/ocndiag_main.f90
ocndiag_LDADD = \
	libDOGCM.la                              \
	$(srcdir)/src/seaice/libDSIce.la         \
	$(LDADD)

#######################################################################

lib_LTLIBRARIES += libDOGCM.la

libDOGCM_la_SOURCES = \
	src/DOGCM_main_mod.f90

libDOGCM_la_LIBADD = \
	src/tint/libDOGCMTInt.la            \
	src/physics/libDOGCMPhysics.la      \
	src/dynamics/libDOGCMDynamics.la    \
	src/io/libDOGCMIO.la                \
	src/boundary/libDOGCMBoundary.la    \
	src/admin/libDOGCMAdmin.la          \
	src/util/libDOGCMUtil.la            \
	src/util/libSpmlUtil.la             \
	libDOGCMExp.la


src/DOGCM_main_mod.f90 : \
	libDOGCMExp.la

######################################################################


noinst_LTLIBRARIES += libDOGCMExp.la

dist_libDOGCMExp_la_SOURCES = \
	exp/$(exp_name)/DOGCM_Exp_$(exp_name)_mod.f90 \
	src/DOGCM_Exp_driver_mod.f90

AM_FCFLAGS += "-DExp_$(exp_name)"

src/DOGCM_Exp_driver_mod.f90 : \
	exp/$(exp_name)/DOGCM_Exp_$(exp_name)_mod.o

