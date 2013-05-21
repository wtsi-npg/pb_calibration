#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <search.h>

#include <smalloc.h>
#include <die.h>
#include <rts.h>

static Header *_hdr = NULL;
static char *filterData = NULL;

// The Perl 'chomp()' function is too useful not to recreate it here...
static void chomp(char *line)
{
    int n = strlen(line) - 1;
    if (line[n] == '\n') line[n] = 0;
    return;
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
    fgets(line, len, fp); hdr->nregions_x = atoi(line);
    fgets(line, len, fp); hdr->nregions_y = atoi(line);
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
    fprintf(fp, "%d\n", hdr->nregions_x);
    fprintf(fp, "%d\n", hdr->nregions_y);
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
	if ((hdr->ntiles * n * hdr->nregions)) {
            filterData = smalloc(hdr->ntiles * n * hdr->nregions);
            if (fread(filterData, hdr->ntiles * n * hdr->nregions, 1, fp) != 1) 
                    die("Error reading filter file\n");
        }

	_hdr = hdr;
}

static int keyComp(const void *k1, const void *k2)
{
	int key = *(int *)k1;
	int mem = *(int *)k2;
	return key != mem;
}

char getFilterData(int tile, int read, int cycle, int region)
{
	int itile, i, cycleLength=0, totalCycleLength=0;
	void *pitile;
	for (i=0; i < read; i++) cycleLength += _hdr->readLength[i];
	for (i=0; i < _hdr->nreads; i++) totalCycleLength += _hdr->readLength[i];
	size_t nelem = _hdr->ntiles;
	pitile = lfind(&tile, _hdr->tileArray, &nelem, sizeof(int), &keyComp);
	if (!pitile) return 0;	// if tile not found in filter
	itile = ((int*)pitile - _hdr->tileArray);
	return filterData[ (itile * totalCycleLength * _hdr->nregions) + ((cycleLength + cycle) * _hdr->nregions) + region];
}

// which region is x in?
int x2region(int x, int region_size)
{
    int coord_shift = (_hdr ? _hdr->coord_shift : COORD_SHIFT);
    int coord_factor = (_hdr ? _hdr->coord_factor : COORD_FACTOR);
    float x_coord = (float)(x - coord_shift) / (float)coord_factor;
    return (int)(x_coord / region_size);
}

// which region is (x,y) in?
int xy2region(int x, int y, int region_size, int nregions_x, int nregions_y)
{
    return x2region(x, region_size) * nregions_y + x2region(y, region_size);
}

