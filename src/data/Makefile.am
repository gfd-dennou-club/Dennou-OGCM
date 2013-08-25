SUBDIRS = 

AM_FCFLAGS = $(SYSFFLAGS) \
	$(MODPATHFLAG)$(srcdir)/ \
	$(MODPATHFLAG)$(srcdir)/../math \
	$(MODPATHFLAG)$(srcdir)/../mesh

#AM_CXXFLAGS = $(SYSCXXFLAGS) \
#	-I $(srcdir)/data/ \
#	-I $(srcdir)/util/ \
#	-I $(srcdir)/io/ \
#	-I $(srcdir)/dynamics/ \
#	-I $(srcdir)/experiment/ \
#	-I $(srcdir)/

#LDADD = \
#	$(SYSLDFLAGS) $(SYSLDLIBS) \
#	$(srcdir)/libDIGMOBaseLib.la 

CLEANFILES = *$(MODEXT)

#bin_PROGRAMS = IGModel-NHRM

#IGModel_NHRM_SOURCES = \
	main.f90

noinst_LTLIBRARIES = libDQGModelData.la
nodist_libDQGModelData_la_SOURCES = \
	List_i_mod.f90 \
	List_d_mod.f90 \
	List_vec3d_mod.f90 \
	List_mod.f90 \
	volScalarField_mod.f90 \
	volVectorField_mod.f90 \
	surfaceScalarField_mod.f90 \
	surfaceVectorField_mod.f90 \
	pointScalarField_mod.f90 \
	pointVectorField_mod.f90 \
	GeometricField_mod.f90
