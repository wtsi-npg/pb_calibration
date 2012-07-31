bin_PROGRAMS = predictor_pu 

predictor_pu_SOURCES = predictor_pu.c

predictor_pu_LDADD = @IO_LIB_LDFLAGS@ @SAMTOOLS_LDFLAGS@ -L../shared -lhelper_funcs

INCLUDES= @IO_LIB_CFLAGS@ @SAMTOOLS_CFLAGS@ -I@top_srcdir@/shared
