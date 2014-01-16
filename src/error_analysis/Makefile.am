bin_PROGRAMS = error_analysis

error_analysis_SOURCES = error_analysis.c

error_analysis_LDADD = @IO_LIB_LDFLAGS@ @SAMTOOLS_LDFLAGS@ -L../shared -lhelper_funcs

INCLUDES= @IO_LIB_CFLAGS@ @SAMTOOLS_CFLAGS@ -I@top_srcdir@/shared
