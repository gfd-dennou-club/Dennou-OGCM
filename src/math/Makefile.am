SUBDIRS = 

AM_FCFLAGS = $(SYSFFLAGS) \
	$(MODPATHFLAG)$(srcdir)/ 

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

noinst_LTLIBRARIES = libDQGModelMath.la
nodist_libDQGModelMath_la_SOURCES = \
	vector2d_mod.f90 \
	vector3d_mod.f90 \
	dotProd_vector2d_mod.f90 \
	dotProd_vector3d_mod.f90 \
	crossProd_mod.f90 \
	VectorSpace_mod.f90 \
	SphericalCoord_mod.f90 \
	geometry_mod.f90 

