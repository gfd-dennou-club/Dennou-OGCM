SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../../include/ \
	$(MODPATHFLAG)$(srcdir)/../setup \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= 

##########


noinst_LTLIBRARIES += libSpmlUtil.la

dist_libSpmlUtil_la_SOURCES= \
	math/at_module_omp.f90 \
	math/SpmlUtil_mod.f90

##########

noinst_LTLIBRARIES += libDOGCMUtil.la 
dist_libDOGCMUtil_la_SOURCES= \
	misc/OptionParser_mod.f90            \
	misc/UnitConversion_mod.f90          \
	misc/ProfUtil_mod.f90                \
	grid/GridIndex_mod.f90               \
	grid/DOGCM_GaussSpmGrid_mod.f90      \
	grid/DOGCM_GaussSpmVFvmGrid_mod.f90  \
	eos/EOS_Linear_mod.f90               \
	eos/EOS_SimpleNonLinear_mod.f90      \
	eos/SeaWaterProp_FM83_mod.f90        \
	eos/EOS_JM95_mod.f90                 \
	eos/EOSDriver_mod.f90                \
	math/TemporalIntegUtil_mod2.f90      \
	math/VSpmUtil_mod.f90                \
	math/VFvmUtil_mod.f90                \
	math/CalculusDriver_mod.f90            


# Dependency

math/CalculusDriver_mod.f90:      \
	math/VSpmUtil_mod.o       \
	math/VFvmUtil_mod.o     

grid/DOGCM_GaussSpmGrid_mod.f90: \
	../admin/DOGCM_Admin_Constants_mod.o

grid/DOGCM_GaussSpmVFvmGrid_mod.f90: \
	../admin/DOGCM_Admin_Constants_mod.o

eos/EOSDriver_mod.f90: \
	../admin/DOGCM_Admin_Constants_mod.o \
	./eos/EOS_Linear_mod.o \
	./eos/EOS_SimpleNonLinear_mod.o \
	./eos/EOS_JM95_mod.o



