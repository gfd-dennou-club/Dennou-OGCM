SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../include/        \
        $(MODPATHFLAG)$(srcdir)/../util            \
        $(MODPATHFLAG)$(srcdir)/../admin           \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

lib_LTLIBRARIES= libDSIce.la

dist_libDSIce_la_SOURCES = \
	admin/DSIce_Admin_Constants_mod.f90          \
	admin/DSIce_Admin_GridDef_mod.f90            \
	admin/DSIce_Admin_GaussSpmGrid_mod.f90       \
	admin/DSIce_Admin_GovernEq_mod.f90           \
	admin/DSIce_Admin_Grid_mod.f90               \
	admin/DSIce_Admin_TInteg_mod.f90             \
	admin/DSIce_Admin_Variable_mod.f90           \
	io/DSIce_IO_History_mod.f90                  \
	io/DSIce_IO_Restart_mod.f90                  \
	boundary/DSIce_Boundary_Vars_mod.f90         \
	boundary/DSIce_Boundary_common_mod.f90       \
	boundary/DSIce_Boundary_driver_mod.f90       \
	thermodyn/DSIce_ThermoDyn_Winton2000_mod.f90 \
	dynamics/DSIce_Dyn_fvm_mod.f90               \
	dynamics/DSIce_Dyn_driver_mod.f90            \
	tint/DSIce_TInt_common_mod.f90               \
	tint/DSIce_TInt_driver_mod.f90               \
	DSIce_main_mod.f90                     

io/DSIce_IO_History_mod.f90: \
	admin/DSIce_Admin_Grid_mod.o

admin/DSIce_Admin_Variable_mod.f90: \
	io/DSIce_IO_History_mod.o   \
	io/DSIce_IO_Restart_mod.o


#seaice_thermdyn_mom_SOURCES=\
	SeaIceThermDyn_mom_main.f90 \
	ice_thm.F90

#seaice_thermdyn_mom_LDADD=\
	$(LDADD) \
	libDSOGCMSEAICE.la
