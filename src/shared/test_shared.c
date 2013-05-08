#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "shared.h"

void checkStr(char *desc, char *expected, char *got)
{
	if (strcmp(expected,got)) {
		fprintf(stderr,"%s\nExpected: %s\nGot     : %s\n", desc, expected, got);
		exit(1);
	}
//	fprintf(stderr,"%s\nExpected: %s\nGot     :%s\n", desc, expected, got);
}

void checkInt(char *desc, int expected, int got)
{
	if (expected != got) {
		fprintf(stderr,"%s\nExpected: %d\nGot     : %d\n", desc, expected, got);
		exit(1);
	}
//	fprintf(stderr,"%s\nExpected: %s\nGot     :%s\n", desc, expected, got);
}

int main(int argc, char *argv[])
{
	int ia[5] = {1,2,3,4,5};
	char hello[] = "Hello world";
	char seq[] = "ACGT";
        int word;
        char *str;

	checkInt("reverse_int(1)", 1, ia[0]);
	reverse_int(ia,5);
	checkInt("reverse_int(2)", 5, ia[0]);

	checkStr("reverse_seq(1)", "Hello world", hello);
	reverse_seq(hello);
	checkStr("reverse_seq(2)", "dlrow olleH", hello);

	word = str2word(seq, 4);
	checkInt("str2word(1)", 27, word);
	str = word2str(228, 4);
	checkStr("word2str(1)", "TGCA", str);

	return 0;
}

