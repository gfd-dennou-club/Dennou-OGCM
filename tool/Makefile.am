SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../src/ \
        $(MODPATHFLAG)$(srcdir)/../src/math \
        $(MODPATHFLAG)$(srcdir)/../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../src/data \
        $(MODPATHFLAG)$(srcdir)/../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../src/io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../src/io/libDQGModelIO.la \
	$(srcdir)/../src/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../src/libDQGModelCore.la \
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

