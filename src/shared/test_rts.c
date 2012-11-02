#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "rts.h"

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
	Header hdr;
	Header h;
	int read;
	int tileArray[1];
	FILE *fp;

	tileArray[0] = 1;

	hdr.region_magic = REGION_MAGIC;
    hdr.coord_shift = COORD_SHIFT;
    hdr.coord_factor = COORD_FACTOR;
    hdr.ntiles = 1;
    hdr.tileArray = &tileArray;
    hdr.nreads = N_READS;
    hdr.region_size = 200;
    hdr.nregions = 100;
    hdr.nregions_x = 10;
    hdr.nregions_y = 11;
    for (read=0; read < hdr.nreads; read++)
        hdr.readLength[read] = 50;
    hdr.cmdLine = argv[0];
    hdr.ncomments = 0;
	addHeaderComment(&hdr, "Greetings universe");

	fp = fopen("x.filter", "w");
	writeHeader(fp,&hdr);
	fclose(fp);

	fp = fopen("x.filter", "r");
	readHeader(fp, &h);
	fclose(fp);

	checkStr("Region Magic", REGION_MAGIC, h.region_magic);
	checkInt("coord_shift", COORD_SHIFT, h.coord_shift);
	checkInt("coord_factor", COORD_FACTOR, h.coord_factor);
	checkInt("ntiles", 1, h.ntiles);
	checkInt("region_size", 200, h.region_size);
	checkInt("nregion", 100, h.nregions);
	checkInt("nregions_x", 10, h.nregions_x);
	checkInt("nregions_y", 11, h.nregions_y);
	checkStr("cmdLine", argv[0], h.cmdLine);
	return 0;
}

