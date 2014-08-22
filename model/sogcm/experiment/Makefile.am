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
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)


bin_PROGRAMS = \
	exp_WindDrivenCirculation/script/numsol_check_eq \
	exp_WindDrivenCirculation/check_numsol


exp_WindDrivenCirculation_script_numsol_check_eq_SOURCES= \
	exp_WindDrivenCirculation/script/numsol_check_eq.f90


exp_WindDrivenCirculation_check_numsol_SOURCES= \
	exp_WindDrivenCirculation/check_numsol.f90

exp_WindDrivenCirculation/check_numsol.f90: ../InitCond_mod.o

exp_WindDrivenCirculation_check_numsol_LDADD= \
	$(srcdir)/../InitCond_mod.o \
	$(srcdir)/../experiment/libDSOGCMExp.la \
	$(srcdir)/../setup/libDSOGCMSetup.la \
	$(srcdir)/../util/libDSOGCMUtil.la \
	$(srcdir)/../util/libSpmlUtil.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)


noinst_LTLIBRARIES= libDSOGCMExp.la 
nodist_libDSOGCMExp_la_SOURCES= \
	exp_W94_Case2/Exp_W94_Case2_mod.f90 \
	exp_BarotRossbyWave/Exp_BarotRossbyWave_mod.f90 \
	exp_InternalGravWave/Exp_InternalGravWave_mod.f90 \
	exp_WindDrivenCirculation/Exp_WindDrivenCirculation_mod.f90 \
	exp_APEOGCirc/Exp_APEOGCirc_mod.f90

