SUBDIRS = 

AM_FCFLAGS = $(SYSFFLAGS) \
	$(MODPATHFLAG)$(srcdir)/ \
	$(MODPATHFLAG)$(srcdir)/../math \
	$(MODPATHFLAG)$(srcdir)/../mesh \
	$(MODPATHFLAG)$(srcdir)/../data \
	$(MODPATHFLAG)$(srcdir)/../fvm


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

lib_LTLIBRARIES = \
	libDQGModelIO.la

libDQGModelIO_la_SOURCES = \
	vtkDataWriter_mod.f90 \
	netcdfDataHelper_mod.f90 \
	netcdfDataWriter_mod.f90 \
	netcdfDataReader_mod.f90 



libDQGModelIO_la_LIBADD = 




