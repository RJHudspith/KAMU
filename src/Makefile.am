TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libKAMU.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/HEADERS/

CUBAFILES=./CUBA/hcubature.c

INTEGRANDFILES=./INTEGRANDS/integrands.c ./INTEGRANDS/shifted.c

IOFILES=./IO/cmdline.c

## all the source files apart from amu_integral*.c
libKAMU_a_SOURCES = ${CUBAFILES} ${INTEGRANDFILES} ${IOFILES}

## and their respective headers
include_HEADERS = \
	$(CUBAFILES:./CUBA/%.c=./HEADERS/%.h) \
	$(INTEGRANDFILES:./INTEGRANDS/%.c=./HEADERS/%.h) \
	$(IOFILES:./IO/%.c=./HEADERS/%.h) \
	../config.h

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = KAMU KINT

KAMU_SOURCES = amu_integral.c
KAMU_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/HEADERS/
KAMU_LDADD = libKAMU.a ${LDFLAGS}

KINT_SOURCES = amu_integrand.c
KINT_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/HEADERS/
KINT_LDADD = libKAMU.a ${LDFLAGS}

endif
