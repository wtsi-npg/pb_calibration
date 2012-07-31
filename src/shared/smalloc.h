#ifndef SMALLOC_H_INCLUDED
#define SMALLOC_H_INCLUDED

void * _s_calloc(size_t nmemb, size_t size,
		 const char *file, unsigned int line, const char *func);

void * _s_malloc(size_t size,
		 const char *file, unsigned int line, const char *func);

void * _s_realloc(void *ptr, size_t size,
		  const char *file, unsigned int line, const char *func);

char * _s_strdup(const char *s,
		 const char *file, unsigned int line, const char *func);

/* If using C9X, we can have function names in out error messages */

#if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#  define __SMALLOC_FUNC __func__
#else
#  define __SMALLOC_FUNC ((const char *) 0)
#endif

#define scalloc(n, s)  _s_calloc((n), (s),  __FILE__, __LINE__, __SMALLOC_FUNC)
#define smalloc(s)     _s_malloc((s),       __FILE__, __LINE__, __SMALLOC_FUNC)
#define srealloc(p, s) _s_realloc((p), (s), __FILE__, __LINE__, __SMALLOC_FUNC)
#define sstrdup(s)     _s_strdup((s),       __FILE__, __LINE__, __SMALLOC_FUNC)


#endif /* SMALLOC_H_INCLUDED */
