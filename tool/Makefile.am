SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../include/ \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../libsrc/io/libDQGModelIO.la \
	$(srcdir)/../libsrc/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../libsrc/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

bin_PROGRAMS += scvgridGen
#scvgridGen_SOURCES =\
	gridGen/ShoreLine_mod.f90 \
	gridGen/gridGen.f90
scvgridGen_SOURCES =\
	gridGen/ShoreLine_SphericalSegment_mod.f90 \
	gridGen/gridGen.f90

bin_PROGRAMS += nc2vtk
nc2vtk_SOURCES =\
	nc2vtkConvert/nc2vtk_main.f90

