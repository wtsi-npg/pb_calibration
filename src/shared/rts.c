#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <smalloc.h>
#include <die.h>
#include <rts.h>

static RegionTable *rts = NULL;
static int n_tiles = 0;
static int n_regions = 0;
static int n_reads = 0;
static int n_cycles = 0;

static inline int getIndex(int itile, int iregion, int iread, int icycle)
{
	return (itile * n_regions * n_reads * n_cycles) + (iregion * n_reads * n_cycles) + (iread * n_cycles) + icycle;
}

void initRTS(int tiles, int regions, int reads, int cycles) 
{
	n_tiles = tiles;
	n_regions = regions;
	n_reads = reads;
	n_cycles = cycles;
	rts = smalloc(n_tiles * n_regions * n_reads * n_cycles * sizeof(RegionTable));
}

void freeRTS(void)
{
	if (rts) free(rts);
	rts = NULL;
}

RegionTable *getRTS(int itile, int iregion, int iread, int icycle)
{
	if (!rts) die("RTS has not been initialised\n");

	if ((itile < 0) || (itile >= n_tiles)) die("Error in getRTS: itile is %d: must be in [0..%d]\n", itile, n_tiles-1);
	if ((iregion < 0) || (iregion >= n_regions)) die("Error in getRTS: iregion is %d: must be in [0..%d]\n", iregion, n_regions-1);
	if ((iread < 0) || (iread >= n_reads)) die("Error in getRTS: iread is %d: must be in [0..%d]\n", iread, n_reads-1);
	if ((icycle < 0) || (icycle >= n_cycles)) die("Error in getRTS: icycle is %d: must be in [0..%d]\n", icycle, n_cycles-1);

	int index = getIndex(itile,iregion,iread,icycle);
	return rts + index;
}

#ifdef TEST
int main(int argc, char *argv[]) 
{
	RegionTable *rt;
	int i;

	initRTS(2,3,4,5);

	rt = getRTS(0,0,0,0); i = getIndex(0,0,0,0);
	printf("getRTS(0,0,0,0) = %d (%p)\n", i, rt);
	rt = getRTS(0,0,0,1);
	printf("getRTS(0,0,0,1) = %d (%p)\n", i, rt);
	rt = getRTS(0,0,0,2);
	printf("getRTS(0,0,0,2) = %d (%p)\n", i, rt);
	rt = getRTS(0,0,0,3);
	printf("getRTS(0,0,0,3) = %d (%p)\n", i, rt);
	rt = getRTS(0,0,0,4);
	printf("getRTS(0,0,0,4) = %d (%p)\n", i, rt);
	rt = getRTS(0,0,1,0);
	printf("getRTS(0,0,1,0) = %d (%p)\n", i, rt);
	rt = getRTS(0,0,1,1);
	printf("getRTS(0,0,1,1) = %d (%p)\n", i, rt);
	rt = getRTS(1,2,3,4);
	printf("getRTS(1,2,3,4) = %d (%p)\n", i, rt);
}
#endif

