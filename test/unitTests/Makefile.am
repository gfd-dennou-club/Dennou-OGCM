SUBDIRS =

test: check

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../src/io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../src/io/libDQGModelIO.la \
	$(srcdir)/../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

TESTS = $(check_PROGRAMS)

check_PROGRAMS = 

check_PROGRAMS += test_geometricField_mod
test_geometricField_mod_SOURCES =\
	test_geometricField_mod.f90

check_PROGRAMS += test_io_mod
test_io_mod_SOURCES =\
	test_io_mod.f90
