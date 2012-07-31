#ifdef HAVE_CONFIG_H
#include "pb_config.h"
#endif

#ifdef ASPRINTF_NEEDS_GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdarg.h>

#ifdef HAVE_ASPRINTF
char * aprintf(const char *fmt, ...) {
  char *str = NULL;
  int res;
  va_list ap;

  va_start(ap, fmt);
  res = vasprintf(&str, fmt, ap);
  va_end(ap);

  if (res >= 0) return str;
  return NULL;
}
  
#else

/* Pretty much from the man page... */
char * aprintf(const char *fmt, ...) {
    int n, size = 100; /* Initial guess for size */
    char *p, *np;
    va_list ap;

    p = malloc(size);
    if (NULL == p) {
        perror("aprintf");
        exit(EXIT_FAILURE);
    }

    while (1) {
        /* Try to print in the allocated space. */
        va_start(ap, fmt);
        n = vsnprintf (p, size, fmt, ap);
        va_end(ap);
        /* If that worked, return the string. */
        if (n > -1 && n < size)
            return p;
        /* Else try again with more space. */
        if (n > -1)    /* glibc 2.1 */
            size = n+1; /* precisely what is needed */
        else           /* glibc 2.0 */
            size *= 2;  /* twice the old size */
        if ((np = realloc (p, size)) == NULL) {
            perror("aprintf");
            exit(EXIT_FAILURE);
        } else {
            p = np;
        }
    }
}

#endif /* HAVE_ASPRINTF */
