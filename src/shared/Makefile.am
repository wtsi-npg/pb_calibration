
noinst_LIBRARIES = libhelper_funcs.a

libhelper_funcs_a_SOURCES = aprintf.c die.c smalloc.c rts.c shared.c snp.c parse_bam.c cif.c \
                            aprintf.h die.h smalloc.h rts.h shared.h snp.h parse_bam.h cif.h

INCLUDES = @IO_LIB_CFLAGS@ @SAMTOOLS_CFLAGS@

test_rts_LDADD = -L../shared -lhelper_funcs
test_aprintf_LDADD = -L../shared -lhelper_funcs
test_shared_LDADD = -L../shared -lhelper_funcs
TESTS = test_rts test_aprintf test_shared
check_PROGRAMS = $(TESTS)

BUILT_SOURCES = version.h

version:

version.h: version
	echo '#ifndef PB_VERSION' > version.h
	echo '#define PB_VERSION' >> version.h
	echo -n const static char *version = \" >> version.h
	../../bin/gitver | tr -d "\n" >> version.h
	echo \"\; >> version.h
	echo '#endif' >> version.h
