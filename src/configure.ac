dnl Process this file with autoconf to produce a configure script.
AC_INIT(pb_cal, 8.0.0)
AM_CONFIG_HEADER(pb_config.h)
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

dnl Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB

dnl C language features
AX_CFLAGS_WARN_ALL
AC_C_INLINE

dnl Functions
AX_FUNC_ASPRINTF
AC_CHECK_FUNCS(popen)

dnl Large file support
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

dnl Look for Staden io_lib
AX_LIB_IO_LIB(1.12.1,
	      [have_io_lib=yes],
	      [AC_MSG_ERROR([Staden io_lib not found])])

AX_LIB_ZLIB(,[have_zlib=yes], [AC_MSG_ERROR([Abort: no zlib. Please rerun configure using the --with-zlib=DIR option.])])
AX_LIB_SAMTOOLS(0.1.3,[have_bam=yes],[AC_MSG_ERROR([Abort: no samtools. Please rerun configure using the ---with-samtools=DIR option.])])

AC_OUTPUT(Makefile shared/Makefile calibration_pu/Makefile predictor_pu/Makefile spatial_filter/Makefile)
