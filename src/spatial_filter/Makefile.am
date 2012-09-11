bin_PROGRAMS = spatial_filter

spatial_filter_SOURCES = spatial_filter.c
BUILT_SOURCES = version.h

version:

version.h: version
	echo -n const char *version = \" > version.h
	../../bin/gitver | tr -d "\n" >> version.h
	echo \"\; >> version.h

spatial_filter_LDADD = @IO_LIB_LDFLAGS@ @SAMTOOLS_LDFLAGS@ -L../shared -lhelper_funcs

INCLUDES= @IO_LIB_CFLAGS@ @SAMTOOLS_CFLAGS@ -I@top_srcdir@/shared
