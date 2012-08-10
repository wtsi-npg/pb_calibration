#ifndef RTS_H_INCLUDED
#define RTS_H_INCLUDED

typedef struct {
    long align;		// FIXME do these really need to be 'long'? Can we make them 'int'?
    long mismatch;
    long insertion;
    long deletion;
    long soft_clip;
    int known_snp;	// FIXME can we get rid of this?
//    long no_call;
	char state;
} RegionTable;

void initRTS(int n_tiles, int n_regions, int n_reads, int n_cycles);
void freeRTS(void);
RegionTable *getRTS(int itile, int iregion, int iread, int icycle);

#endif /* DIE_H_INCLUDED */
