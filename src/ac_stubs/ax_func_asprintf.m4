AC_DEFUN([AX_FUNC_ASPRINTF],
[
AC_LANG_ASSERT([C])
AC_MSG_CHECKING([for asprintf])
AC_EGREP_CPP([asprintf],[#include <stdio.h>],
[_asprintf_needs_gnu_source=no],[
AC_EGREP_CPP([asprintf],[#define _GNU_SOURCE
#include <stdio.h>],
[_asprintf_needs_gnu_source=yes])
])
AS_IF([test "x$_asprintf_needs_gnu_source" = "xyes"],
[AC_DEFINE([ASPRINTF_NEEDS_GNU_SOURCE], [1], [asprintf needs _GNU_SOURCE])])
AC_LINK_IFELSE(AC_LANG_PROGRAM([[
#ifdef ASPRINTF_NEEDS_GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
]],
[[
char *str;
return asprintf(&str, "foo%d", 1);
]]),
[AC_MSG_RESULT([yes])
AC_DEFINE([HAVE_ASPRINTF], [1], [Have asprintf])
],
[AC_MSG_RESULT([no])
])
]
)
