#ifndef CIF_H_INCLUDED
#define CIF_H_INCLUDED

#include <stdio.h>

typedef struct {
    long        lane;
    long        cycle;
    long        read;
    const char *dir;
} CifDir;

typedef union {
    int8_t  *i8;
    int16_t *i16;
    int32_t *i32;
} ChanData;

typedef struct {
    char *filename;
    int fd;
    int data_type;
    int num_channels;
    int cycle_number;
    off_t  cycle_start_pos;
    size_t num_entries;
    size_t chunk_num_entries;
    size_t chunk_start;
    ChanData *chan_data;
} CifCycleData;

typedef struct {
    CifCycleData *cycles;
    size_t        ncycles;
    size_t        size;
    size_t        num_spots;
} CifData;

CifDir *cif_dirs;
size_t n_cif_dirs;
size_t *cif_lane_index;
size_t n_cif_lanes;

void get_cif_dirs(char *intensity_dir);
void read_cif_chunk(CifCycleData *cycle, size_t spot_num);
int read_cif_file(char *name, int fd, CifData *cif_data);
CifData *load_cif_data(int lane, int tile, char *suffix);
void free_cif_data(CifData *cif_data);

#endif /* CIF_H_INCLUDED */
