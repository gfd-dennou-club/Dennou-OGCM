SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../../src/io \
        $(MODPATHFLAG)$(srcdir)/../ \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../util \
        $(MODPATHFLAG)$(srcdir)/../boundary_cond \
        $(MODPATHFLAG)$(srcdir)/../seaice \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)


#################################################################

noinst_LTLIBRARIES= libDSOGCMExp.la 
dist_libDSOGCMExp_la_SOURCES= \
	InitCond_mod.f90 \
	exp_W94_Case2/Exp_W94_Case2_mod.f90 \
	exp_BarotRossbyWave/Exp_BarotRossbyWave_mod.f90 \
	exp_InternalGravWave/Exp_InternalGravWave_mod.f90 \
	exp_WindDrivenCirculation/Exp_WindDrivenCirculation_mod.f90 \
	exp_APEOGCirc/Exp_APEOGCirc_mod.f90 \
	exp_APEOGCircSeaice/Exp_APEOGCircSeaice_mod.f90 \
	exp_APECoupledAOGCMSeaice/Exp_APECoupledAOGCMSeaice_mod.f90


###########################################################

bin_PROGRAMS = \
	exp_WindDrivenCirculation/script/numsol_check_eq \
	exp_WindDrivenCirculation/check_numsol \
	exp_SeaIceThermTest/seaicetherm_S78Cases


################################################################

# Some programs to analyze the result of exp_WindDrivenCirculation

exp_WindDrivenCirculation_script_numsol_check_eq_SOURCES= \
	exp_WindDrivenCirculation/script/integrand_mod.f90   \
	exp_WindDrivenCirculation/script/numsol_check_eq.f90


exp_WindDrivenCirculation_check_numsol_SOURCES= \
	exp_WindDrivenCirculation/check_numsol.f90

exp_WindDrivenCirculation_check_numsol_LDADD= \
	$(srcdir)/libDSOGCMExp.la \
	$(srcdir)/../seaice/libDSOGCMSeaIce.la \
	$(srcdir)/../boundary_cond/libDSOGCMBC.la \
	$(srcdir)/../setup/libDSOGCMSetup.la \
	$(srcdir)/../util/libDSOGCMUtil.la \
	$(srcdir)/../util/libSpmlUtil.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

exp_WindDrivenCirculation_check_numsol: libDSOGCMExp.la

################################################################

# A program to validiate the implementin of sea ice thermdynamics model.

exp_SeaIceThermTest_seaicetherm_S78Cases_SOURCES= \
	exp_SeaIceThermTest/seaicetherm_S78Cases.f90

exp_SeaIceThermTest_seaicetherm_S78Cases_LDADD= \
	$(srcdir)/../seaice/libDSOGCMSeaIce.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)



