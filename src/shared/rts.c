#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <search.h>

#include <smalloc.h>
#include <die.h>
#include <rts.h>

static RegionTable *rts[N_READS] = {NULL, NULL, NULL};
static Header *_hdr;
static int n_tiles = 0;
static int n_regions = 0;
static int n_cycles[N_READS] = {0, 0, 0};
static char *filterData = NULL;

// The Perl 'chomp()' function is too useful not to recreate it here...
static void chomp(char *line)
{
    int n = strlen(line) - 1;
    if (line[n] == '\n') line[n] = 0;
    return;
}

//
// RTS Methods
//
static inline int getIndex(int itile, int iregion, int iread, int icycle)
{
	return (itile * n_regions * n_cycles[iread]) + (icycle * n_regions) + iregion;
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
	if ((itile < 0) || (itile >= n_tiles)) die("Error in getRTS: itile is %d: must be in [0..%d]\n", itile, n_tiles-1);
	if ((iregion < 0) || (iregion >= n_regions)) die("Error in getRTS: iregion is %d: must be in [0..%d]\n", iregion, n_regions-1);
	if ((iread < 0) || (iread >= N_READS)) die("Error in getRTS: iread is %d: must be in [0..%d]\n", iread, N_READS-1);
	if ((icycle < 0) || (icycle >= n_cycles[iread])) die("Error in getRTS: icycle is %d: must be in [0..%d]\n", icycle, n_cycles[iread]-1);

	if (!rts[iread]) die("RTS has not been initialised for iread %d\n", iread);

	int index = getIndex(itile,iregion,iread,icycle);
	return rts[iread] + index;
}

//
// Filter methods
//

void readHeader(FILE *fp, Header *hdr)
{
    int len = 1024;
    char line[1024];
    int i;

    fgets(line, len, fp); chomp(line); hdr->region_magic = strdup(line);
    fgets(line, len, fp); hdr->coord_shift = atoi(line);
    fgets(line, len, fp); hdr->coord_factor = atoi(line);
    fgets(line, len, fp); hdr->region_size = atoi(line);
    fgets(line, len, fp); hdr->ntiles = atoi(line);
	hdr->tileArray = smalloc(hdr->ntiles * sizeof(int));
	for (i=0; i < hdr->ntiles; i++) {
		fgets(line, len, fp); hdr->tileArray[i] = atoi(line);
	}
    fgets(line, len, fp); hdr->nregions = atoi(line);
    fgets(line, len, fp); hdr->nreads = atoi(line);
    for (i=0; i < hdr->nreads; i++) {
        fgets(line, len, fp); hdr->readLength[i] = atoi(line);
    }
    fgets(line, len, fp); chomp(line); hdr->cmdLine = strdup(line);
    fgets(line, len, fp); hdr->ncomments = atoi(line);
    for (i=0; i < hdr->ncomments; i++) {
        fgets(line, len, fp); chomp(line); hdr->comments[i] = strdup(line);
    }
}

void writeHeader(FILE *fp, Header *hdr)
{
    int i;
    fprintf(fp, "%s\n", hdr->region_magic);
    fprintf(fp, "%d\n", hdr->coord_shift);
    fprintf(fp, "%d\n", hdr->coord_factor);
    fprintf(fp, "%d\n", hdr->region_size);
    fprintf(fp, "%d\n", hdr->ntiles);
	for (i=0; i < hdr->ntiles; i++) {
		fprintf(fp, "%d\n", hdr->tileArray[i]);
	}
    fprintf(fp, "%d\n", hdr->nregions);
    fprintf(fp, "%d\n", hdr->nreads);
    for (i=0; i<hdr->nreads; i++)
        fprintf(fp, "%d\n", hdr->readLength[i]);
    fprintf(fp, "%s\n", hdr->cmdLine);
    fprintf(fp, "%d\n", hdr->ncomments);
    for (i=0; i<hdr->ncomments; i++)
        fprintf(fp, "%s\n", hdr->comments[i]);
}

void addHeaderComment(Header *hdr, char *comment)
{
    hdr->comments[hdr->ncomments++] = strdup(comment);
}

// Allocate memory for the filter file, then slurp the whole thing
void readFilterData(FILE *fp, Header *hdr)
{
	int i, n=0;
	for (i=0; i < hdr->nreads; i++) n += hdr->readLength[i];
	filterData = smalloc(hdr->ntiles * n * hdr->nregions);
	if (fread(filterData, hdr->ntiles * n * hdr->nregions, 1, fp) != 1) 
		die("Error reading filter file\n");

	_hdr = hdr;
}

static int keyComp(const void *k1, const void *k2)
{
	int key = *(int *)k1;
	int mem = *(int *)k2;
	return key == mem;
}

char getFilterData(int tile, int read, int cycle, int region)
{
	int itile, i, cycleLength=0;
	for (i=0; i < read; i++) cycleLength += _hdr->readLength[i];
	size_t nelem = _hdr->ntiles;
	itile = *(int *)lfind(&tile, _hdr->tileArray, &nelem, sizeof(int), &keyComp);
	return filterData[ (itile * cycleLength) + region];
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

