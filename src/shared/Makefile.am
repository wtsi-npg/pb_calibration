noinst_LIBRARIES = libhelper_funcs.a

libhelper_funcs_a_SOURCES = aprintf.c die.c smalloc.c rts.c \
                            aprintf.h die.h smalloc.h rts.h

BUILT_SOURCES = version.h

version:

version.h: version
	echo -n const char *version = \" > version.h
	../../bin/gitver | tr -d "\n" >> version.h
	echo \"\; >> version.h
