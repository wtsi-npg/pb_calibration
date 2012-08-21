#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <smalloc.h>
#include <die.h>
#include <rts.h>

#define N_READS 3

static RegionTable *rts[N_READS] = {NULL, NULL, NULL};
static int n_tiles = 0;
static int n_regions = 0;
static int n_cycles[N_READS] = {0, 0, 0};

  

static inline int getIndex(int itile, int iregion, int iread, int icycle)
{
	return (itile * n_regions * n_cycles[iread]) + (iregion * n_cycles[iread]) + icycle;
}

void initRTS(int tiles, int regions, int iread, int cycles) 
{
        if ((n_tiles > 0) && (tiles != n_tiles)) die("Error in initRTS: inconsistent ntiles: have %d expected %d\n", tiles, n_tiles);
        if ((n_regions > 0) && (regions != n_regions)) die("Error in initRTS: inconsistent nregions: have %d expected %d\n", regions, n_regions);
	if ((iread < 0) || (iread >= N_READS)) die("Error in initRTS: iread is %d: must be in [0..%d]\n", iread, N_READS-1);

        n_tiles = tiles;
	n_regions = regions;
	n_cycles[iread] = cycles;
	rts[iread] = smalloc(n_tiles * n_regions * n_cycles[iread] * sizeof(RegionTable));
}

void freeRTS(void)
{
        int iread;
        for (iread=0; iread<N_READS; iread++) {
                  if (rts[iread]) free(rts[iread]);
                  rts[iread] = NULL;
        }
}

RegionTable *getRTS(int itile, int iregion, int iread, int icycle)
{
	if (!rts) die("RTS has not been initialised\n");

	if ((itile < 0) || (itile >= n_tiles)) die("Error in getRTS: itile is %d: must be in [0..%d]\n", itile, n_tiles-1);
	if ((iregion < 0) || (iregion >= n_regions)) die("Error in getRTS: iregion is %d: must be in [0..%d]\n", iregion, n_regions-1);
	if ((iread < 0) || (iread >= N_READS)) die("Error in getRTS: iread is %d: must be in [0..%d]\n", iread, N_READS-1);
	if ((icycle < 0) || (icycle >= n_cycles[iread])) die("Error in getRTS: icycle is %d: must be in [0..%d]\n", icycle, n_cycles[iread]-1);

	int index = getIndex(itile,iregion,iread,icycle);
	return rts[iread] + index;
}

#ifdef TEST
int main(int argc, char *argv[]) 
{
	RegionTable *rt;
	int iread, i;

	for (iread=0; iread<N_READS; iread++)
                initRTS(2,3,iread,5);

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

