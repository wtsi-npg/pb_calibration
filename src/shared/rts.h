#ifndef RTS_H_INCLUDED
#define RTS_H_INCLUDED

#include <stdio.h>

#define N_READS 3
#define N_COMMENTS 100

#define COORD_SHIFT   1000
#define COORD_FACTOR  10

#define REGION_MAGIC                "RGFL"
#define REGION_SIZE                 200

#define REGION_STATE_COVERAGE   (1<<1)
#define REGION_STATE_MISMATCH   (1<<2)
#define REGION_STATE_INSERTION  (1<<3)
#define REGION_STATE_DELETION   (1<<4)
#define REGION_STATE_SOFT_CLIP  (1<<5)
#define REGION_STATE_BAD        (1<<6)

// The header of the filter file
typedef struct {
    char *region_magic;
    int coord_shift;
    int coord_factor;
    int ntiles;
    int *tileArray;
    size_t *tileReadCountArray;
    int region_size;
    int nregions;
    int nregions_x;
    int nregions_y;
    int nreads;
    int readLength[N_READS];
    int totalReadLength;
    char *cmdLine;
    int ncomments;
    char *comments[N_COMMENTS];
} Header;

// An internal structure used to create the filter file
typedef struct {
        int align;
        int mismatch;
        int insertion;
        int deletion;
        int soft_clip;
        int known_snp;
        float quality;
	char state;
} RegionTable;

// Filter methods
void writeHeader(FILE *fp, Header *hdr);
void addHeaderComment(Header *hdr, char *comment);
void readHeader(FILE *fp, Header *hdr);
void readFilterData(FILE *fp, Header *hdr);
char* getFilterData(int tile, int read, int cycle, int region);

int x2region(int x, int region_size);
int xy2region(int x, int y, int region_size, int nregions_x, int nregions_y);

#endif /* RTS_H_INCLUDED */
