#ifndef RTS_H_INCLUDED
#define RTS_H_INCLUDED

typedef struct {
        int align;
        int mismatch;
        int insertion;
        int deletion;
        int soft_clip;
        int known_snp;
	char state;
} RegionTable;

void initRTS(int n_tiles, int n_regions, int iread, int n_cycles);
void freeRTS(void);
RegionTable *getRTS(int itile, int iregion, int iread, int icycle);

#endif /* DIE_H_INCLUDED */
