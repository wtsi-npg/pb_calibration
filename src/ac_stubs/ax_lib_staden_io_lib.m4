# SYNOPSIS
#
#   AX_LIB_IO_LIB([MINIMUM-VERSION], [ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of Staden io_lib library.
#   (http://staden.sourceforge.net/). It does this by running the
#   io_lib-config script, hence needing io_lib 1.10.0 or newer (where
#   that program first appeared).
#
#   The following output variables are set using AC_SUBST:
#
#     IO_LIB_CFLAGS
#     IO_LIB_LDFLAGS
#     IO_LIB_VERSION (if MINIMUM-VERSION is not "")
#
# LICENSE
#
#   Copyright (c) 2009 James Bonfield <jkb@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty
#   provided the copyright notice and this notice are preserved.
#   This file is offered as-is, without any warranty.

# AC_CHECK_IO_LIB([DEFAULT-ACTION], [MINIMUM-VERSION])
# Autoconf macro to find io_lib.
# If found it defines HAVE_IO_LIB and sets IO_LIB_CFLAGS and IO_LIB_LDFLAGS.
#
# DEFAULT-ACTION is the string "yes" or "no", defaulting to "yes" when not
# specified.
# Requires io_lib 1.10.0 or above (for the io_lib-config script).

AC_DEFUN([AX_LIB_IO_LIB],
[
  AC_ARG_WITH(io_lib, AC_HELP_STRING([--with-io_lib=DIR],
	                             [Look for io_lib root in DIR]),
	      [_io_lib_with=$withval],
	      [_io_lib_with=""])

  if test "x$_io_lib_with" = "x"
  then
    _iopath=$PATH
  else
    _iopath="$_io_lib_with/bin"
  fi

  AC_PATH_PROG([_io_lib_config], [staden-io_lib-config],,[$_iopath])
  if test "x$_io_lib_config" = "x"
  then
      AC_PATH_PROG([_io_lib_config], [io_lib-config],,[$_iopath])
  fi

  if test "x$_io_lib_config" != "x"
  then
    # Check version is sufficient; sneakily entirely in sh syntax
    if test "x$1" != "x"
    then
      AC_MSG_CHECKING([if io_lib version >= $1])

      _io_lib_version=`$_io_lib_config --version`
      SAVE_IFS=$IFS; IFS=.
      _val=0
      for v in $_io_lib_version; do _val=`expr $_val '*' 1000 + $v`; done
      _io_lib_version=$_val
      IFS=$SAVE_IFS

      _io_lib_wanted=`echo ifelse([$1],,[0],[$1])`
      SAVE_IFS=$IFS; IFS=.
      _val=0
      for v in $_io_lib_wanted; do _val=`expr $_val '*' 1000 + $v`; done
      _io_lib_wanted=$_val
      IFS=$SAVE_IFS

      if test $_io_lib_version -ge $_io_lib_wanted
      then
        AC_MSG_RESULT([yes])
        io_lib_version_ok=yes
      else
        AC_MSG_RESULT([no])
        io_lib_version_ok=no
      fi
      AC_SUBST([IO_LIB_VERSION])
    else
      io_lib_version_ok=yes;
    fi

    if test $io_lib_version_ok = "yes" 
    then    
      # Configure IO_LIB_CFLAGS and IO_LIB_LDFLAGS
      test x"$IO_LIB_CFLAGS" = "x" && IO_LIB_CFLAGS=`$_io_lib_config --cflags`
      test x"$IO_LIB_LDFLAGS"  = "x" && IO_LIB_LDFLAGS=`$_io_lib_config --libs`
      AC_DEFINE(HAVE_IO_LIB, 1, [Define to 1 if you have a working io_lib])
      AC_SUBST(IO_LIB_CFLAGS)
      AC_SUBST(IO_LIB_LDFLAGS)
    fi
  fi

  # Execute the conditional expressions
  if test "$io_lib_version_ok" = "yes"
  then
     # This is the IF-YES path
     ifelse([$2],,:,[$2])
  else
     # This is the IF-NO path
     ifelse([$3],,:,[$3])
  fi

])dnl
