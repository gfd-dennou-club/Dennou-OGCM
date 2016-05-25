SUBDIRS =

all: copy_header

copy_header:
		@echo "Copy generated header files into $(srcdir)/include"
		@$(MKDIR_P) $(srcdir)/include
		@find . -name "*@MODEXT@" -exec $(install_sh_DATA) {} $(srcdir)/include \;

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../include/ \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../libsrc/io/libDQGModelIO.la \
	$(srcdir)/../../libsrc/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../../libsrc/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

bin_PROGRAMS += globalSWM 
globalSWM_SOURCES =\
	SimParameters_mod.f90 \
	GridSet_mod.f90 \
	VariableSet_mod.f90 \
	GovernEquationSolver_mod.f90 \
	OutputData_mod.f90 \
	Exp_Williamson94_Case1.f90 \
	Exp_Williamson94_Case2.f90 \
	Exp_Williamson94_Case5.f90 \
	globalSWModel_main.f90

