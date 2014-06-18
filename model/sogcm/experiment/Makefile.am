SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../../src/io \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../util \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)


bin_PROGRAMS = \
	exp_WindDrivenCirculation/script/numsol_check_eq


exp_WindDrivenCirculation_script_numsol_check_eq_SOURCES= \
	exp_WindDrivenCirculation/script/numsol_check_eq.f90



noinst_LTLIBRARIES= libDSOGCMExp.la 
nodist_libDSOGCMExp_la_SOURCES= \
	exp_W94_Case2/Exp_W94_Case2_mod.f90 \
	exp_BarotRossbyWave/Exp_BarotRossbyWave_mod.f90 \
	exp_InternalGravWave/Exp_InternalGravWave_mod.f90 \
	exp_WindDrivenCirculation/Exp_WindDrivenCirculation_mod.f90
