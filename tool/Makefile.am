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
        $(MODPATHFLAG)$(srcdir)/../experiment \
        $(MODPATHFLAG)$(srcdir)/../governeq_driver \
        $(MODPATHFLAG)$(srcdir)/../seaice \
        $(MODPATHFLAG)$(srcdir)/../physics_driver \
        $(MODPATHFLAG)$(srcdir)/../dynamics \
        $(MODPATHFLAG)$(srcdir)/../sgs_param \
	$(MODPATHFLAG)$(srcdir)/../boundary_cond \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../util \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../experiment/libDSOGCMExp.la \
	$(srcdir)/../governeq_driver/libDSOGCMGovernEqDriver.la \
	$(srcdir)/../physics_driver/libDSOGCMPhysicsDriver.la \
	$(srcdir)/../boundary_cond/libDSOGCMBC.la \
	$(srcdir)/../seaice/libDSOGCMSeaIce.la \
	$(srcdir)/../dynamics/libDSOGCMDynamics.la \
	$(srcdir)/../setup/libDSOGCMSetup.la \
	$(srcdir)/../util/libDSOGCMUtil.la \
	$(srcdir)/../util/libSpmlUtil.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

###################################

if DSOGCM_3D

bin_PROGRAMS=diagVar3d
diagVar3d_SOURCES =\
	diagvar/DiagVarSet_mod.f90 \
	diagvar/DiagVarEval_mod.f90 \
	diagvar/DiagVarFileSet_mod.f90 \
	diagvar/EnergyBudgetAnalysis_mod.f90 \
	diagvar/EnergyBudgetAnalysisT2013_mod.f90 \
	diagvar/BudgetAnalysis_mod.f90    \
	diagvar/SpectralAnalysis_mod.f90  \
	diagvar/SGSEddyMixAnalysis_mod.f90 \
	diagvar/GMSchemeAnalysis_mod.f90 \
	diagvar/diagVar_main.f90

endif

if DSOGCM_AXISYM

bin_PROGRAMS=diagVar_axisym
diagVar_axisym_SOURCES =\
	diagvar/DiagVarSet_mod.f90 \
	diagvar/DiagVarEval_mod.f90 \
	diagvar/DiagVarFileSet_mod.f90 \
	diagvar/EnergyBudgetAnalysis_mod.f90 \
	diagvar/EnergyBudgetAnalysisT2013_mod.f90 \
	diagvar/BudgetAnalysis_mod.f90 \
	diagvar/SpectralAnalysis_mod.f90 \
	diagvar/GMSchemeAnalysis_mod.f90 \
	diagvar/SGSEddyMixAnalysis_mod.f90 \
	diagvar/diagVar_main.f90

endif
