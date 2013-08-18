SUBDIRS = 

AM_FCFLAGS = $(SYSFFLAGS) \
	$(MODPATHFLAG)$(srcdir)/ \
	$(MODPATHFLAG)$(srcdir)/../math 

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
	Field_i_mod.f90 \
	Field_d_mod.f90 \
	Field_vec3d_mod.f90



