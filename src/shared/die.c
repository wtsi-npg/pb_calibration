#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <die.h>

void display(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}

void die(const char *fmt, ...)
{
	va_list ap;
	va_start(ap,fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	exit(EXIT_FAILURE);
}
