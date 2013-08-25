subdir =

test: check

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../src/data \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS) \
	$(srcdir)/../../src/libDQGModelCore.la 

CLEANFILES = *$(MODEXT)

TESTS = $(check_PROGRAMS)

check_PROGRAMS = 

check_PROGRAMS += test_geometricField_mod
test_geometricField_mod_SOURCES =\
	temp.f90 \
	test_geometricField_mod.f90
