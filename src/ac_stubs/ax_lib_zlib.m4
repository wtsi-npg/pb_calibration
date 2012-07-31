# SYNOPSIS
#
#   AX_LIB_ZLIB([MINIMUM-VERSION], [ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of zlib library.
#   It does this by checking for the header file zlib.h and the z library
#   object file. The location of these may be specified using the
#   --with-zlib=DIR command line option (eg --with-zlib=/usr/local),
#   using $DIR/include and $DIR/lib for the search path.
#
#   The following output variables are set using AC_SUBST:
#
#     ZLIB_CFLAGS
#     ZLIB_LDFLAGS
#     ZLIB_VERSION (if MINIMUM-VERSION is not "")
#
#   The C preprocessor symbol HAVE_ZLIB will be also defined with
#   AC_DEFINE if a functioning samtools is available.
#
# LICENSE
#
#   Copyright (c) 2009 James Bonfield <jkb@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty
#   provided the copyright notice and this notice are preserved.
#   This file is offered as-is, without any warranty.


AC_DEFUN([AX_LIB_ZLIB],
[
  AC_ARG_WITH(zlib,
	      AC_HELP_STRING([--with-zlib=DIR],[look for zlib in DIR]),
	      [_zlib_with=$withval],[_zlib_with="no"])

  ZLIB_ROOT=""
  if test "$_zlib_with" != "no"
  then
     if test -f "$_zlib_with/include/zlib.h"
     then
         ZLIB_ROOT=$_zlib_with
     fi
  fi

  # Check if it's a working library
  zlib_ok=no
  if test "x$ZLIB_ROOT" != "x"
  then
    _cppflags=$CPPFLAGS
    CPPFLAGS="$CPPFLAGS -I${ZLIB_ROOT}/include"
    _ldflags=$LDFLAGS
    LDFLAGS="$LFDLAGS -L${ZLIB_ROOT}/lib"
    AC_LANG_PUSH([C])
    AC_CHECK_LIB(z, inflateEnd,
	[AC_CHECK_HEADER(zlib.h, zlib_ok=yes, zlib_ok=no)])
    AC_LANG_POP([C])
    if test "$zlib_ok" != "yes"
    then
        # Backout and whinge
        CPPFLAGS=$_cppflags
        LDFLAGS=$_ldflags
        AC_MSG_WARN("--with-zlib specified, but non functioning")
    fi

  else
    # Maybe it works "out of the box"?
    AC_CHECK_LIB(z, inflateEnd,
	[AC_CHECK_HEADER(zlib.h, zlib_ok=yes, zlib_ok=no)])
  fi

  # Check version
  if test "x$1" != "x" && test "$zlib_ok" = "yes"
  then
      AC_MSG_CHECKING([if zlib version >= $1])

      for i in "$ZLIB_ROOT/include" "/usr/include" "/usr/share/include" "/usr/local/include"
      do
	  if test -f "$i/zlib.h"
	  then
	      ZLIB_VERSION=`sed -n 's/.*#define *ZLIB_VERSION *"\([^"]*\)"/\1/p' "$i/zlib.h"`
	      break
	  fi
      done

      v1=`expr "$1" : '\([[0-9]]*\)'`
      v2=`expr "$1" : '[[0-9]]*\.\([[0-9]]*\)'`
      v3=`expr "$1" : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
      want_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`

      v1=`expr "${ZLIB_VERSION:-}" : '\([[0-9]]*\)'`
      v2=`expr "${ZLIB_VERSION:-}" : '[[0-9]]*\.\([[0-9]]*\)'`
      v3=`expr "${ZLIB_VERSION:-}" : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
      have_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`
      if test `expr "$have_vers" ">=" "$want_vers"` = "1"
      then
          AC_MSG_RESULT([yes])
          AC_SUBST([ZLIB_VERSION])
      else
          AC_MSG_RESULT([no])
	  zlib_ok="no"
      fi
  fi

  # perform substitutions
  if test "$zlib_ok" = "yes"
  then
      AC_DEFINE(HAVE_ZLIB, 1,
         [Define to 1 if you have a functional libz.])
      if test "$ZLIB_ROOT" != ""
      then
          ZLIB_LDFLAGS="-L${ZLIB_ROOT}/lib -lz"
	  ZLIB_CFLAGS="-I${ZLIB_ROOT}/include"
      else
          ZLIB_LDFLAGS="-lz"
	  ZLIB_CFLAGS=
      fi
      AC_SUBST([ZLIB_LDFLAGS])
      AC_SUBST([ZLIB_CFLAGS])
  else
    AC_MSG_WARN("No functioning zlib found")
  fi

  # Not sure how many of these are needed, but it's belt-and-braces mode
  AH_TEMPLATE([HAVE_ZLIB], [Define if zlib is installed])
  AM_CONDITIONAL(HAVE_ZLIB, test "$zlib_ok" = "yes")


  # Execute the conditional expressions
  if test "$zlib_ok" = "yes"
  then
     # This is the IF-YES path
     ifelse([$2],,:,[$2])
  else
     # This is the IF-NO path
     ifelse([$3],,:,[$3])
  fi

  # Tidy up
  unset zlib_ok
  unset _cppflags
  unset _ldflags
])
