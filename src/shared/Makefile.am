
noinst_LIBRARIES = libhelper_funcs.a

libhelper_funcs_a_SOURCES = aprintf.c die.c smalloc.c rts.c \
                            aprintf.h die.h smalloc.h rts.h

test_rts_LDADD = -L../shared -lhelper_funcs
test_aprintf_LDADD = -L../shared -lhelper_funcs
TESTS = test_rts test_aprintf
check_PROGRAMS = $(TESTS)
