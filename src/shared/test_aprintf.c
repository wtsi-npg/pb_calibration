#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "aprintf.h"

void checkString(char *desc, char *expected, char *got)
{
	if (strcmp(expected,got)) {
		fprintf(stderr,"%s\nExpected: %s\nGot     : %s\n", desc, expected, got);
		exit(1);
	}
//	fprintf(stderr,"%s\nExpected: %s\nGot     :%s\n", desc, expected, got);
}

int main(int argc, char *argv[])
{
	char *x = aprintf("Hello world");
	checkString("aprintf(str)", "Hello world", x);

	x = aprintf("pi = %4.2f", 3.14);
	checkString("aprintf(int)", "pi = 3.14", x);
	return 0;
}

