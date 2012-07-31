#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <smalloc.h>
#include <die.h>

void * _s_calloc(size_t nmemb, size_t size,
		 const char *file, unsigned int line, const char *func) {
  void *m = calloc(nmemb, size);

  if (NULL == m) {
    if (func) {
      die("Couldn't allocate %zd bytes in %s at %s line %u: %s\n",
	  nmemb * size, func, file, line, strerror(errno));
    } else {
      die("Couldn't allocate %zd bytes at %s line %u: %s\n",
	   nmemb * size, file, line, strerror(errno));
    }
  }
  return m;
}

void * _s_malloc(size_t size,
		 const char *file, unsigned int line, const char *func) {
  void *m = malloc(size);

  if (NULL == m) {
    if (func) {
      die("Couldn't allocate %zd bytes in %s at %s line %u: %s\n",
	  size, func, file, line, strerror(errno));
    } else {
      die("Couldn't allocate %zd bytes at %s line %u: %s\n",
	  size, file, line, strerror(errno));
    }
  }
  return m;
}

void * _s_realloc(void *ptr, size_t size,
		  const char *file, unsigned int line, const char *func) {
  void *m = realloc(ptr, size);

  if (NULL == m) {
    if (func) {
      die("Couldn't reallocate %zd bytes in %s at %s line %u: %s\n",
	  size, func, file, line, strerror(errno));
    } else {
      die("Couldn't reallocate %zd bytes at %s line %u: %s\n",
	  size, file, line, strerror(errno));
    }
  }
  return m;
}

char * _s_strdup(const char *s,
		 const char *file, unsigned int line, const char *func) {
  char *str = strdup(s);

  if (NULL == str) {
    if (func) {
      die("Couldn't duplicate string of %zd bytes in %s at %s line %u: %s\n",
	  strlen(str), func, file, line, strerror(errno));
    } else {
      die("Couldn't duplicate string of %zd bytes at %s line %u: %s\n",
	  strlen(str), file, line, strerror(errno));
    }
  }

  return str;
}
