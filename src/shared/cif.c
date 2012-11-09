/*
 * Copyright (c) 2006-2009, Genome Research Ltd (GRL).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *
 *     * Neither the name of the Genome Research Limited nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GRL ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL GRL BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifdef HAVE_PREAD
# define _XOPEN_SOURCE 500 // for pread
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <search.h>
#include <glob.h>
#include <regex.h>
#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <io_lib/misc.h>
#include <cif.h>
#include <smalloc.h>
#include <die.h>
#include <rts.h>

#define MAX_CIF_CHUNK_BYTES 4194304

#ifdef HAVE_PREAD
# warning "trying to use system pread"
# define pread_bytes pread
#else
# warning "using Rob's(?) pread"
/* Read count bytes at position offset in file fd.  It actually emulates
   pread(2) as the real thing doesn't seem to be any faster, and it
   may not be present. */
static ssize_t pread_bytes(int fd, void *buf, size_t count, off_t offset) 
{
	char *b = (char *) buf;
	ssize_t res;
	ssize_t total = 0;

	if (lseek(fd, offset, SEEK_SET) < 0) return -1;

	do {
		do {
			res = read(fd, b + total, count - total);
		} while (res < 0 && EINTR == errno);
		if (res > 0) total += res;
	} while (res > 0 && total < count);
	return res < 0 ? res : total;
}
#endif


static int cif_dir_compare(const void *va, const void *vb) 
{
	const CifDir *a = (const CifDir *) va;
	const CifDir *b = (const CifDir *) vb;

	if (a->lane  < b->lane)  return -1;
	if (a->lane  > b->lane)  return  1;
	if (a->read  < b->read)  return -1;
	if (a->read  > b->read)  return  1;
	if (a->cycle < b->cycle) return -1;
	if (a->cycle > b->cycle) return  1;
	return strcmp(a->dir, b->dir);
}

void get_cif_dirs(char *intensity_dir) 
{
	char msg[1024];
    char   *pattern = NULL;
    size_t  pattern_sz = 0;
    glob_t  glob_buf;
    regex_t    re;
    regmatch_t matches[4];
    size_t  i, l;
    int     res;

	n_cif_lanes = 0;
    pattern_sz = strlen(intensity_dir) + 30;
    pattern    = smalloc(pattern_sz);

    snprintf(pattern, pattern_sz, "%s/L*/C*", intensity_dir);

    memset(&glob_buf, 0, sizeof(glob_buf));
    memset(&re,       0, sizeof(re));

    res = glob(pattern, GLOB_ERR|GLOB_NOSORT|GLOB_NOESCAPE, NULL, &glob_buf);
    if (0 != res) {
        if (GLOB_NOMATCH != res) {
            perror("Looking for CIF directories");
        }
        globfree(&glob_buf);
        free(pattern);
        return;
    }

    cif_dirs = smalloc(glob_buf.gl_pathc * sizeof(CifDir));

    res = regcomp(&re, "/L([[:digit:]]+)/C([[:digit:]]+)\\.([[:digit:]]+)/*$", REG_EXTENDED);
    if (0 != res) {
        regerror(res, &re, msg, sizeof(msg));
        die("Regular expression compile failed in get_cif_dirs: %s\n", msg);
	}

    for (n_cif_dirs = 0, i = 0; i < glob_buf.gl_pathc; i++) {
        const char *p = glob_buf.gl_pathv[i];

        res = regexec(&re, p, sizeof(matches)/sizeof(matches[0]), matches, 0);
        if (REG_NOMATCH == res) continue;
        if (0 != res) {
			regerror(res, &re, msg, sizeof(msg));
			die("Regular expression match failed in get_cif_dirs: %s\n", msg);
		}
        cif_dirs[n_cif_dirs].lane  = strtol(p + matches[1].rm_so, NULL, 10);
        cif_dirs[n_cif_dirs].cycle = strtol(p + matches[2].rm_so, NULL, 10);
        cif_dirs[n_cif_dirs].read  = strtol(p + matches[3].rm_so, NULL, 10);
        cif_dirs[n_cif_dirs].dir   = strdup(p);
        if (NULL == cif_dirs[n_cif_dirs].dir) die("Out of memory in get_cif_dirs()");
        if (cif_dirs[n_cif_dirs].lane > n_cif_lanes) n_cif_lanes = cif_dirs[n_cif_dirs].lane;
        n_cif_dirs++;
    }
    regfree(&re);
    globfree(&glob_buf);

    if (0 == n_cif_dirs) {
        free(cif_dirs);
        return;
    }

    qsort(cif_dirs, n_cif_dirs, sizeof(*cif_dirs), cif_dir_compare);

    cif_lane_index = scalloc(n_cif_lanes + 2, sizeof(*cif_lane_index));

    for (i = 0, l = 0; i < n_cif_dirs; i++) {
        if (cif_dirs[i].lane > l) {
            long j;

            for (j = l + 1; j <= cif_dirs[i].lane; j++) {
                cif_lane_index[j] = i;
            }
            l = cif_dirs[i].lane;
        }
    }
    cif_lane_index[l + 1] = n_cif_dirs;

}


void read_cif_chunk(CifCycleData *cycle, size_t spot_num) 
{
    size_t num_entries;
    size_t num_bytes;
    size_t chan_size_bytes;
    ssize_t got;
    off_t  file_pos;
    int chan;

    assert(spot_num < cycle->num_entries);
    if (spot_num >= cycle->chunk_start
        && spot_num < cycle->chunk_start + cycle->chunk_num_entries) {
        return;
    }

    chan_size_bytes  = cycle->num_entries * cycle->data_type;
    num_entries  = (cycle->num_entries - spot_num < cycle->chunk_num_entries
                    ? cycle->num_entries - spot_num
                    : cycle->chunk_num_entries);
    num_bytes = num_entries * cycle->data_type;

    if (NULL == cycle->chan_data) {
        size_t chunk_size_bytes = cycle->chunk_num_entries * cycle->data_type;
        cycle->chan_data = smalloc(cycle->num_channels * sizeof(cycle->chan_data[0]));
        for (chan = 0; chan < cycle->num_channels; chan++) {
            cycle->chan_data[chan].i8 = smalloc(chunk_size_bytes);
        }
    }

    for (chan = 0; chan < cycle->num_channels; chan++) {
        file_pos = (cycle->cycle_start_pos
                    + chan_size_bytes * chan
                    + spot_num * cycle->data_type);
        got = pread_bytes(cycle->fd, cycle->chan_data[chan].i8, num_bytes, file_pos);
        if (got < 0) {
            die("Error reading %s: %s\n", cycle->filename, strerror(errno));
        }
        if (got < num_bytes) {
            die("Error: Did not get the expected amount of data from %s\n", cycle->filename);
        }
    }
    cycle->chunk_start = spot_num;

#ifdef WORDS_BIGENDIAN
    /* CIF files are little-endian, so we need to switch everything
       around on big-endian machines */
    switch (cycle->data_type) {
    case 1:
        break;
    case 2:
        for (chan = 0; chan < cycle->num_channels; chan++) {
            size_t i;
            int16_t *d = cycle->chan_data[chan].i16;
            for (i = 0; i < num_entries; i++) {
                d[i] = bswap_16(d[i]);
            }
        }
        break;
    case 4:
        for (chan = 0; chan < cycle->num_channels; chan++) {
            size_t i;
            int32_t *d = cycle->chan_data[chan].i32;
            for (i = 0; i < num_entries; i++) {
                d[i] = bswap_32(d[i]);
            }
        }
        break;
    default:
        die("Unexpected data type in CIF file %s\n", cycle->filename);
    }
#endif
}


int read_cif_file(char *name, int fd, CifData *cif_data) 
{
    uint8_t cif_header[13];
    int    num_cycles;
    int    first_cycle;
    int    num_channels = 4;
    size_t num_entries;
    size_t chunk_num_entries;
    size_t cycle_size_bytes;
    off_t  cycle_start_pos = sizeof(cif_header);
    int    data_type;
    CifCycleData *cycle;
    size_t i;

    if (sizeof(cif_header) != pread_bytes(fd, cif_header, sizeof(cif_header), 0)) return -1;
    if (0 != memcmp(cif_header, "CIF", 3))                 return -1;
    if (cif_header[3] != '\1')                             return -1;

    data_type   = cif_header[4];
    first_cycle = cif_header[5] | (cif_header[6] << 8);
    num_cycles  = cif_header[7] | (cif_header[8] << 8);
    num_entries = (cif_header[9]
                   | (cif_header[10] << 8)
                   | (cif_header[11] << 16)
                   | (cif_header[12] << 24));

    if (0 == cif_data->num_spots) {
        cif_data->num_spots = num_entries;
    } else if (cif_data->num_spots != num_entries) {
        fprintf(stderr, "Got unexpected number of entries in CIF file %s\n"
                "Expected: %zd; Got %zd\n",
                name, cif_data->num_spots, num_entries);
    }
    if (data_type < 1 || data_type > 4) {
        die("Unexpected data_type in CIF file %s\n", name);
    }

    chunk_num_entries = MAX_CIF_CHUNK_BYTES / (data_type * num_channels);
    if (chunk_num_entries > num_entries) chunk_num_entries = num_entries;
    cycle_size_bytes = num_entries * data_type * num_channels;

    for (i = 0; i < num_cycles; i++) {
        if (cif_data->ncycles == cif_data->size) {
            cif_data->size = cif_data->size > 0 ? cif_data->size * 2 : 128;
            cif_data->cycles = srealloc(cif_data->cycles,
                                        cif_data->size * sizeof(CifCycleData));
        }

        cycle = cif_data->cycles + cif_data->ncycles;

        cycle->filename     = sstrdup(name);
        cycle->fd           = fd;
        cycle->data_type    = data_type;
        cycle->num_channels = num_channels;
        cycle->cycle_number = i + first_cycle;
        cycle->cycle_start_pos = cycle_start_pos + i * cycle_size_bytes;
        cycle->num_entries  = num_entries;
        cycle->chunk_num_entries = chunk_num_entries;
        cycle->chunk_start  = num_entries; /* Will force first load */
        cycle->chan_data    = NULL; /* Allocate memory on first load */

        cif_data->ncycles++;
    }

    return 0;
}


CifData *load_cif_data(int lane, int tile, char *suffix) 
{
    size_t first_dir;
    size_t end_dir;
    size_t i;
    int cif = -1;
    CifData *cif_data = NULL;

	n_cif_lanes = 0;
    cif_data = scalloc(1, sizeof(CifData));

    first_dir = (lane <= n_cif_lanes
                 ? cif_lane_index[lane]
                 : cif_lane_index[n_cif_lanes + 1]);
    end_dir = (lane <= n_cif_lanes
               ? cif_lane_index[lane + 1]
               : cif_lane_index[n_cif_lanes + 1]);

    for (i = first_dir; i < end_dir; i++) {
		size_t cif_file_sz = cif_file_sz = strlen(cif_dirs[i].dir) + strlen(suffix) + 100;
		char *cif_file = smalloc(cif_file_sz);
        snprintf(cif_file, cif_file_sz, "%s/s_%d_%d.%s", cif_dirs[i].dir, lane, tile, suffix);
        cif = open(cif_file, O_RDONLY);
        if (cif < 0) {
            if (ENOENT == errno) continue;
            die("Couldn't open %s : %s\n", cif_file, strerror(errno));
        }

        if (0 != read_cif_file(cif_file, cif, cif_data)) {
            die("Error reading %s\n", cif_file);
        }
		free(cif_file);
		close(cif);
    }

    /* Check we have a full set of cycles */
    for (i = 0; i < cif_data->ncycles; i++) {
        if (cif_data->cycles[i].cycle_number > i + 1) {
            die("Error: Missing cycle %zd for lane %d tile %d from CIF files.\n", i+1, lane, tile);
        } else if (cif_data->cycles[i].cycle_number != i + 1) {
            die("Error: Unexpected cycle number %d for lane %d tile %d from CIF files.\n", cif_data->cycles[i].cycle_number, lane, tile);
        }
    }

    return cif_data;
}


void free_cif_data(CifData *cif_data) 
{
    size_t i;
    int c;

    if (NULL == cif_data->cycles)
        return;

    for (i = 0; i < cif_data->ncycles; i++) {
        if (NULL != cif_data->cycles[i].chan_data) {
            for (c = 0; c < cif_data->cycles[i].num_channels; c++) {
                if (NULL != cif_data->cycles[i].chan_data[c].i8)
                    free(cif_data->cycles[i].chan_data[c].i8);
            }
            free(cif_data->cycles[i].chan_data);
        }
        if (0 < cif_data->cycles[i].fd) {
            if (0 != close(cif_data->cycles[i].fd)) {
                if (NULL != cif_data->cycles[i].filename)
                    die("Error when closing %s: %s\n", cif_data->cycles[i].filename, strerror(errno));
                else
                    die("Error when closing cif file: %s\n", strerror(errno));
            }
        }
        if (NULL != cif_data->cycles[i].filename)
            free(cif_data->cycles[i].filename);
    }

    free(cif_data->cycles);

    free(cif_data);
}


