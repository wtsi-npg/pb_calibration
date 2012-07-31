bin_PROGRAMS = spatial_filter

spatial_filter_SOURCES = spatial_filter.c

spatial_filter_LDADD = @IO_LIB_LDFLAGS@ @SAMTOOLS_LDFLAGS@ -L../shared -lhelper_funcs

INCLUDES= @IO_LIB_CFLAGS@ @SAMTOOLS_CFLAGS@ -I@top_srcdir@/shared
