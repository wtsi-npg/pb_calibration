bin_PROGRAMS = calibration_pu

calibration_pu_SOURCES = calibration_pu.c

calibration_pu_LDADD = @IO_LIB_LDFLAGS@ @SAMTOOLS_LDFLAGS@ -L../shared -lhelper_funcs

INCLUDES= @IO_LIB_CFLAGS@ @SAMTOOLS_CFLAGS@ -I@top_srcdir@/shared
