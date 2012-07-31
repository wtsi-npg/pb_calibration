check_PROGRAMS = fake_data check_srf

EXTRA_DIST = all_srf.pl srf_test.pl

INCLUDES= -I@top_srcdir@/shared -I@top_srcdir@/pb_calibration @IO_LIB_CFLAGS@

fake_data_SOURCES = fake_data.c ../pb_calibration/pb_calibration_util_zfp.c
fake_data_LDADD = -L../shared -lhelper_funcs

check_srf_SOURCES = check_srf.c
check_srf_LDADD = @IO_LIB_LDFLAGS@ -L../shared -lhelper_funcs

check-local:
	srcdir=${srcdir} ; export srcdir ; \
	perl ${srcdir}/all_srf.pl

clean-local:
	-rm -f test*.out
