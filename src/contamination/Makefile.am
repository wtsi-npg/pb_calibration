bin_PROGRAMS = get_mixture get_thr_ploidy

get_mixture_SOURCES = get_mixture.c

get_mixture_LDADD = @IO_LIB_LDFLAGS@ @SAMTOOLS_LDFLAGS@ -L../shared -lhelper_funcs

get_thr_ploidy_SOURCES = get_thr_ploidy.c

get_thr_ploidy_LDADD = @IO_LIB_LDFLAGS@ @SAMTOOLS_LDFLAGS@ -L../shared -lhelper_funcs

AM_CPPFLAGS = @IO_LIB_CFLAGS@ @SAMTOOLS_CFLAGS@ -I@top_srcdir@/shared
