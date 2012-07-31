#include <pb_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include <die.h>
#include <smalloc.h>
#include <aprintf.h>

/* Hack to stop io_lib from trying to include its own config.h */
#ifdef HAVE_CONFIG_H
#undef HAVE_CONFIG_H
#endif
#include <io_lib/srf.h>
#include <io_lib/mFILE.h>
#include <io_lib/ztr.h>

typedef struct {
  char *bustard_dir;
  char *qseq_dir;
  char *qseq_suffix;
  int   filter_bad_reads;
  int   verbosity;
} Settings;

typedef struct {
  char    *regn_meta;
  size_t   regn_meta_sz;
  uint8_t *regn;
  size_t   regn_sz;
  char    *text;
  size_t   text_sz;
} Previous_data;

typedef struct {
  char *key;
  char *value;
} Key_value_pair;

typedef struct {
  char   *machine;
  size_t  machine_sz;
  int     run_id;
  int     lane;
  int     tile;
  int     x;
  int     y;
} Spot_info;

typedef struct {
  int lane;
  int tile;
  int nqseqs;
  char **qseq_names;
  FILE **qseqs;
  size_t seq_len;
  char *seq;
  char *qual;
} Trace_reader;

int key_compare(const void *va, const void *vb) {
  Key_value_pair *a = (Key_value_pair *) va;
  Key_value_pair *b = (Key_value_pair *) vb;

  return strcmp(a->key, b->key);
}

static inline char * find_text(Key_value_pair *kvp, size_t n_kvp, char *key) {
  Key_value_pair  to_find = { key, NULL };
  Key_value_pair *kv;

  kv = bsearch(&to_find, kvp, n_kvp, sizeof(*kvp), key_compare);

  return NULL != kv ? kv->value : NULL;
}

int check_file(Settings *opts, char *name, char *text) {
  char buffer[65536];
  char *path = NULL;
  FILE *f    = NULL;
  char *t = text;
  size_t len;
  size_t i;
  int res = 1;

  if (NULL == name) return 1;
  path = aprintf("%s/%s", opts->bustard_dir, name);
  f = fopen(path, "rb");
  if (NULL == f) {
    printf("Couldn't open %s: %s\n", path, strerror(errno));
    goto clean;
  }
  
  do {
    len = fread(buffer, 1, sizeof(buffer), f);
    for (i = 0; i < len; i++) {
      if (*t != buffer[i] || '\0' == *t) break;
      t++;
    }
    if (i < len) {
      printf("Difference found at char %zd reading file %s\n",
	      t - text, path);
      goto clean;
    }
  } while (len == sizeof(buffer));

  if (ferror(f)) {
    printf("Error reading %s: %s\n", path, strerror(errno));
    goto clean;
  }
  if (0 != fclose(f)) {
    printf("Error closing %s: %s\n", path, strerror(errno));
    f = NULL;
    goto clean;
  }
  f = NULL;

  if ('\0' != *t) {
    printf("EOF reading %s\n", path);
    goto clean;
  }

  res = 0;

 clean:
  if (NULL != f) fclose(f);
  if (NULL != path) free(path);

  return res;
}

int check_text(Settings *opts, Previous_data *seen) {
  static Key_value_pair *kvp = NULL;
  static size_t sz_kvp = 0;
  static int iolib_bug_warn = 0;

  size_t  n_kvp = 0;
  size_t  i;
  char   *start;
  char   *end = NULL;
  ssize_t remaining;
  int     res;
  char   *fn;

  if (seen->text_sz < 2) {
    printf("Zero-length TEXT chunk found\n");
    return -1;
  }
  if ('\0' != seen->text[0]) {
    printf("check_text given compressed TEXT data\n");
    return -1;
  }
  
  start     = seen->text + 1;
  remaining = seen->text_sz - 1;
  while (remaining > 0 && '\0' != *start) {
    if (n_kvp == sz_kvp) {
      sz_kvp = sz_kvp ? sz_kvp * 2 : 16;
      kvp = srealloc(kvp, sz_kvp * sizeof(*kvp));
    }
    kvp[n_kvp].key = start;
    end = memchr(start, 0, remaining);
    if (NULL == end) break;                    /* Unterminated key */
    remaining -= end - start + 1;
    if (remaining <= 0) { end = NULL; break; } /* Only key, no value */
    kvp[n_kvp].value = end + 1;
    end = memchr(kvp[n_kvp].value, 0, remaining);
    if (NULL == end) break;                    /* Unterminated value */
    remaining -= end - kvp[n_kvp].value + 1;
    n_kvp++;
    if (remaining <= 0) break;            /* Should be at least 1 byte left */
    start = end + 1;
  }
  if (remaining != 1 || '\0' != *start) {
    if (0 != remaining || NULL == end) {
      printf("Error decoding TEXT chunk\n");
      return -1;
    } else {
      /* Ignore bug in io_lib up to 1.12.1 */
      if (!iolib_bug_warn) {
	printf("\nWarning: Your version of io_lib doesn't include the final\n"
	       "NULL byte to indicate the end of a TEXT chunk.  This means\n"
	       "it strictly breaks the ZTR specification.  It'sbeen like\n"
	       "that for a long time though, and it isn't too serious so\n"
	       "we'll let it go for now.\n\n");
	iolib_bug_warn = 1;
      }
    }
  }

  qsort(kvp, n_kvp, sizeof(*kvp), key_compare);

  /* Look for repeat entries */
  for (i = 1; i < n_kvp; i++) {
    if (0 == strcmp(kvp[i - 1].key, kvp[i].key)) {
      printf("Duplicate identifier found in TEXT chunk\n");
      return -1;
    }
  }

  res  = check_file(opts, "config.xml",
		    find_text(kvp, n_kvp, "ILLUMINA_GA_BUSTARD_CONFIG"));
  res += check_file(opts, "../config.xml",
		    find_text(kvp, n_kvp, "ILLUMINA_GA_FIRECREST_CONFIG"));
  res += check_file(opts, "BustardSummary.xml",
		    find_text(kvp, n_kvp, "ILLUMINA_GA_BUSTARD_SUMMARY"));
  res += check_file(opts,
		    find_text(kvp, n_kvp, "ILLUMINA_GA_MATRIX_FWD_FILENAME"),
		    find_text(kvp, n_kvp, "ILLUMINA_GA_MATRIX_FWD"));
  res += check_file(opts,
		    find_text(kvp, n_kvp, "ILLUMINA_GA_PHASING_FWD_FILENAME"),
		    find_text(kvp, n_kvp, "ILLUMINA_GA_PHASING_FWD"));
  if (NULL != (fn = find_text(kvp, n_kvp, "ILLUMINA_GA_MATRIX_REV_FILENAME"))){
    res += check_file(opts, fn,
		      find_text(kvp, n_kvp, "ILLUMINA_GA_MATRIX_REV"));
  }
  if (NULL != (fn = find_text(kvp, n_kvp,"ILLUMINA_GA_PHASING_REV_FILENAME"))){
    res += check_file(opts, fn,
		      find_text(kvp, n_kvp, "ILLUMINA_GA_PHASING_REV"));
  }

  return res;
}

int check_trace_header(Settings *opts, srf_t *srf, Previous_data *seen) {
  int i;
  int j;
  int n;
  uint8_t  *data;
  uint32_t  val;
  ztr_chunk_t *chunk;

  if (NULL != srf->mf) {
    mfrecreate(srf->mf, NULL, 0);
  } else {
    srf->mf = mfcreate(NULL, 0);
  }
  if (NULL == srf->mf) die("%s\n", strerror(errno));

  if (srf->th.trace_hdr_size) {
    if (1 != mfwrite(srf->th.trace_hdr, srf->th.trace_hdr_size, 1, srf->mf)) {
      die("mfwrite failed: %s\n", strerror(errno));
    }
  }

  if (srf->ztr) delete_ztr(srf->ztr);

  mrewind(srf->mf);

  if (NULL != (srf->ztr = partial_decode_ztr(srf, srf->mf, NULL))) {
    srf->mf_pos = mftell(srf->mf);
  } else {
    /* We expect this to work for srfs made by illumina2srf */
    printf("partial_decode_ztr failed for trace header\n");
    srf->mf_pos = 0;
    return -1;
  }

  mfseek(srf->mf, 0, SEEK_END);
  srf->mf_end = mftell(srf->mf);

  /* Go through chunks */
  for (i = 0; i < srf->ztr->nchunks; i++) {
    chunk = &srf->ztr->chunk[i];
    if (opts->verbosity > 1) {
      printf(" Chunk %d type %.4s\n", i,
	      (char *) &chunk->type);
    }
    switch (chunk->type) {
    case ZTR_TYPE_HUFF:
      break;
    case ZTR_TYPE_BPOS:
      if (0 != uncompress_chunk(srf->ztr, chunk)) {
	printf("Couldn't uncompress BPOS chunk\n");
	return -1;
      }
      n = (chunk->dlength - 4) / 4;
      for (j = 0; j < n; j++) {
	data = (uint8_t *) &chunk->data[j * 4 + 4];
	val = (  ((uint32_t) data[0]) << 24
	       | ((uint32_t) data[1]) << 16
	       | ((uint32_t) data[2]) <<  8
	       | ((uint32_t) data[3]));
	if (val != j) {
	  printf("BPOS data misses cycles\n");
	  return -1;
	}
      }
      break;
    case ZTR_TYPE_REGN:
      if (0 != uncompress_chunk(srf->ztr, chunk)) {
	printf("Couldn't uncompress REGN chunk\n");
	return -1;
      }
      if (NULL == seen->regn) {
	/* Copy REGN chunk */
	seen->regn_meta_sz = chunk->mdlength;
	seen->regn_meta = smalloc(seen->regn_meta_sz);
	memcpy(seen->regn_meta, chunk->mdata, seen->regn_meta_sz);
	seen->regn_sz = chunk->dlength;
	seen->regn = smalloc(seen->regn_sz);
	memcpy(seen->regn, chunk->data, seen->regn_sz);
      } else {
	/* Compare with last copy */
	if (seen->regn_meta_sz != chunk->mdlength
	    || seen->regn_sz != chunk->dlength
	    || 0 != memcmp(seen->regn_meta, chunk->mdata, seen->regn_meta_sz)
	    || 0 != memcmp(seen->regn, chunk->data, seen->regn_sz)) {
	  printf("REGN chunk changed between header blocks\n");
	  return -1;
	}
      }
      break;
    case ZTR_TYPE_TEXT:
      if (0 != uncompress_chunk(srf->ztr, chunk)) {
	printf("Couldn't uncompress REGN chunk\n");
	return -1;
      }
      if (NULL == seen->text) {
	seen->text_sz = chunk->dlength;
	seen->text = smalloc(seen->text_sz);
	memcpy(seen->text, chunk->data, seen->text_sz);
	
	if (0 != check_text(opts, seen)) return -1;
      } else {
	if (seen->text_sz != chunk->dlength
	    || 0 != memcmp(seen->text, chunk->data, seen->text_sz)) {
	  printf("New TEXT chunk found\n");
	  seen->text_sz = chunk->dlength;
	  seen->text = srealloc(seen->text, seen->text_sz);
	  memcpy(seen->text, chunk->data, seen->text_sz);
	  
	  if (0 != check_text(opts, seen)) return -1;
	}
      }
      break;
    default:
      printf("Found unexpected chunk type in header block\n");
      return -1;
    }
  }

  return 0;
}

int decode_name(char *name, Spot_info *spot) {
  int u;

  u = strcspn(name, "_");
  if (name[u] != '_') return -1;

  if (5 != sscanf(name + u + 1, "%d:%d:%d:%d:%d",
		  &spot->run_id, &spot->lane,
		  &spot->tile, &spot->x, &spot->y)) {
    return -1;
  }

  if (spot->machine_sz < u + 1) {
    spot->machine = srealloc(spot->machine, u + 1);
    spot->machine_sz = u + 1;
  }

  memcpy(spot->machine, name, u);
  spot->machine[u] = '\0';

  return 0;
}

int close_qseq_files(Trace_reader *trc) {
  FILE *f;
  int   i;
  int   res = 0;

  for (i = 0; i < trc->nqseqs; i++) {
    if (NULL == (f = trc->qseqs[i])) continue;
    trc->qseqs[i] = NULL;
    if (0 != fclose(f)) {
      printf("Error closing %s: %s\n",
	      trc->qseq_names[i], strerror(errno));
      res = -1;
    }
    if (NULL != trc->qseq_names[i]) free(trc->qseq_names[i]);
    trc->qseq_names[i] = NULL;
  }

  return res;
}

int open_new_qseqs(Settings *opts, Trace_reader *trc, Spot_info *spot) {
  size_t  name_len;
  char   *name;
  int     i;

  if (NULL == trc->qseqs) {
    name_len = strlen(opts->qseq_dir) + strlen(opts->qseq_suffix) + 64;
    name = smalloc(name_len);

    for (i = 0; ; i++) {
      snprintf(name, name_len, "%s/s_%d_%d_%04d_%s",
	       opts->qseq_dir, spot->lane, i + 1,
	       spot->tile, opts->qseq_suffix);
      if (0 != access(name, F_OK)) {
	if (ENOENT == errno) break;
	printf("Error looking up %s: %s\n", name, strerror(errno));
	free(name);
	return -1;
      }
    }

    if (0 == i) {
      printf("Couldn't find qseq files\n");
      free(name);
      return -1;
    }

    trc->nqseqs     = i;
    trc->qseqs      = scalloc(i, sizeof(FILE *));
    trc->qseq_names = scalloc(i, sizeof(char *));
    free(name);

  } else {

    if (0 != close_qseq_files(trc)) return -1;

  }

  for (i = 0; i < trc->nqseqs; i++) {
    trc->qseq_names[i] = aprintf("%s/s_%d_%d_%04d_%s",
				 opts->qseq_dir, spot->lane, i + 1,
				 spot->tile, opts->qseq_suffix);
    trc->qseqs[i] = fopen(trc->qseq_names[i], "r");
    if (NULL == trc->qseqs[i]) {
      printf("Couldn't open %s: %s\n",
	      trc->qseq_names[i], strerror(errno));
      return -1;
    }
  }
  
  return 0;
}

static int read_qseq_line(FILE *f, char *fname, char *fields[11]) {
  static const char *spaces = " \t\n";
  static char       *buffer = NULL;
  static size_t      bufsz = 0;
  size_t             s;
  char              *b;
  char              *res;
  int                i;

  if (NULL == buffer) {
    bufsz = 1024;
    buffer = smalloc(bufsz);
  }

  b = buffer;
  s = bufsz;
  b[s - 1] = '*';
  while (NULL != (res = fgets(b, s, f))) {
    if (b[s - 1] == 0 && b[s - 2] != '\n') {
      bufsz *= 2;
      buffer = srealloc(buffer, bufsz);
      b = buffer + s - 1;
      s = bufsz - (b - buffer);
      b[s - 1] = '*';
    } else {
      break;
    }
  }

  if (NULL == res) {
    if (ferror(f)) {
      printf("Error reading %s: %s\n", fname, strerror(errno));
      return -1;
    } else { /* eof */
      return 1;
    }
  }

  /* (*line_no)++; */
  
  b = buffer;
  s = 0;
  for (i = 0; i < 11; i++) {
    fields[i] = buffer + s;
    s += strcspn(buffer + s, spaces);
    if (buffer[s] == '\0') {
      printf("Error reading qseq file: not enough fields found\n");
      return -1;
    }
    buffer[s++] = '\0';
    s += strspn(buffer + s, spaces);
  }

  return 0;
}

int read_qseqs(Settings *opts, Trace_reader *trc, Spot_info *spot,
	       size_t *skipped) {
  char  *fields[11];
  size_t count   = 0;
  size_t seq_len;
  size_t l;
  int    failed_read;
  int    i;
  
  *skipped = 0;
  do {
    failed_read = 1;
    count++;
    seq_len = 0;
    for (i = 0; i < trc->nqseqs; i++) {
      if (0 != read_qseq_line(trc->qseqs[i], trc->qseq_names[i], fields)) {
	return -1;
      }

      if (atoi(fields[7]) != i + 1) {
	printf("Error: Qseq file %s says it's read %s, expected %d\n",
		trc->qseq_names[i], fields[7], i + 1);
	return -1;
      }

      if (0 == i || trc->nqseqs - 1 == i) {
	failed_read = failed_read && *fields[10] != '1';
      }

      l = strlen(fields[8]);
      if (l != strlen(fields[9])) {
	fprintf(stderr,
		"Error reading qseq file %s: Seq and quality lengths differ\n",
		trc->qseq_names[i]);
	return -1;
      }
      if (trc->seq_len < seq_len + l) {
	trc->seq_len = seq_len + l;
	trc->seq  = srealloc(trc->seq,  trc->seq_len + 1);
	trc->qual = srealloc(trc->qual, trc->seq_len + 1);
      }
      strcpy(trc->seq + seq_len,  fields[8]);
      strcpy(trc->qual + seq_len, fields[9]);
      seq_len += l;
    }
  } while (opts->filter_bad_reads && failed_read);

  *skipped = count - 1;

  if (0 != strcmp(fields[0], spot->machine)
      || atoi(fields[1]) != spot->run_id
      || atoi(fields[2]) != spot->lane
      || atoi(fields[3]) != spot->tile
      || atoi(fields[4]) != spot->x
      || atoi(fields[5]) != spot->y) {
    printf("Out of synch. reading qseq files\n");
    return -1;
  }

  for (i = 0; i < trc->seq_len; i++) trc->qual[i] -= 64;
  return 0;
}

int read_trace_data(Settings *opts, Trace_reader *trc, Spot_info *spot) {
  size_t skipped;

  if (trc->lane != spot->lane || trc->tile != spot->tile) {
    if (0 != open_new_qseqs(opts, trc, spot)) return -1;
    trc->lane = spot->lane;
    trc->tile = spot->tile;
  }

  if (0 != read_qseqs(opts, trc, spot, &skipped)) return -1;
  return 0;
}

int check_slxi_trace(ztr_chunk_t *chunk, Trace_reader *trc) {
  size_t   i;
  uint8_t *data[4] = {
    (uint8_t *) chunk->data + 2,
    (uint8_t *) chunk->data + 2 + trc->seq_len * 2,
    (uint8_t *) chunk->data + 2 + trc->seq_len * 4,
    (uint8_t *) chunk->data + 2 + trc->seq_len * 6
  };
  uint16_t max = 0;
  uint16_t val = 0;
  int      max_base;
  int      base;

  for (i = 0; i < trc->seq_len; i++) {
    max_base = -1;
    max = 0;
    for (base = 0; base < 4; base++) {
      val = ((uint16_t) data[base][i * 2] << 8) | data[base][i * 2 + 1];
      if (max <= val) {
	max = val;
	max_base = base;
      }
    }
    if (trc->seq[i] != "CATG"[max_base]) {
      printf("SLXI trace does not match sequence\n");
      return -1;
    }
  }
  return 0;
}

int check_slxn_trace(ztr_t *ztr, ztr_chunk_t *chunk, Trace_reader *trc) {
  char    *offs;
  int      offset;
  size_t   i;
  int      base;
  uint16_t val;
  uint8_t *data;

  offs = ztr_lookup_mdata_value(ztr, chunk, "OFFS");
  if (NULL == offs) {
    printf("Couldn't get OFFS metadata from SLXN trace chunk\n");
    return -1;
  }

  offset = atoi(offs);

  for (base = 0; base < 4; base++) {
    data = (uint8_t *) chunk->data + base * trc->seq_len * 2 + 2;
    for (i = 0; i < trc->seq_len; i++) {
      val = (((uint16_t) data[i * 2] << 8) | data[i * 2 + 1]);
      if (val - offset != (((base & 1) == 0) ? 1 : -1)) {
	printf("Incorrect SLXN data\n");
	return -1;
      }
    }
  }
  return 0;
}

int check_proc_trace(ztr_chunk_t *chunk, Trace_reader *trc) {
  size_t   i;
  uint8_t *data[4] = {
    (uint8_t *) chunk->data + 2,
    (uint8_t *) chunk->data + 2 + trc->seq_len * 2,
    (uint8_t *) chunk->data + 2 + trc->seq_len * 4,
    (uint8_t *) chunk->data + 2 + trc->seq_len * 6
  };
  uint16_t max = 0;
  uint16_t val = 0;
  int      max_base;
  int      base;

  for (i = 0; i < trc->seq_len; i++) {
    max_base = -1;
    max = 0;
    for (base = 0; base < 4; base++) {
      val = ((uint16_t) data[base][i * 2] << 8) | data[base][i * 2 + 1];
      if (max <= val) {
	max = val;
	max_base = base;
      }
    }
    if (trc->seq[i] != "ACGT"[max_base]) {
      printf("PROC trace does not match sequence\n");
      return -1;
    }
  }
  return 0;
}

int check_trace(ztr_t *ztr, ztr_chunk_t *chunk, Trace_reader *trc) {
  char *type;
  
  if (chunk->dlength < trc->seq_len * 8 + 2) {
    fprintf(stderr,
	    "Trace is not long enough to account for all base calls\n");
    return -1;
  }

  type = ztr_lookup_mdata_value(ztr, chunk, "TYPE");
  if (NULL == type) {
    type = "PROC";
  }

  switch (type[0]) {
  case 'P':
    if (0 == strcmp(type + 1, "ROC")) {
      return check_proc_trace(chunk, trc);
    }
    break;
  case 'S':
    if ('L' == type[1] && 'X' == type[2]) {
      switch (type[3]) {
      case 'I':
	return check_slxi_trace(chunk, trc);
      case 'N':
	return check_slxn_trace(ztr, chunk, trc);
      default:
	break;
      }
    }
  default:
    break;
  }

  printf("Unexpected trace type '%s' found\n", type);
  return -1;
}

int get_second_calls(ztr_t *ztr, size_t nbases, int **indexes) {
  ztr_chunk_t *smp4_chunk = NULL;
  char        *type = NULL;
  size_t       i;
  int max_base;
  int second_base;
  int base;
  uint16_t     val;
  uint16_t     max;
  uint16_t     second;
  uint8_t     *data[4];

  for (i = 0; i < ztr->nchunks; i++) {
    if (ZTR_TYPE_SMP4 == ztr->chunk[i].type) {
      type = ztr_lookup_mdata_value(ztr, &ztr->chunk[i], "TYPE");
      if (NULL == type
	  || 0 == strcmp(type, "PROC")
	  || 0 == strcmp(type, "SLXI")) {
	smp4_chunk = &ztr->chunk[i];
	break;
      }
    }
  }

  if (NULL == smp4_chunk) return 1;

  if (0 != uncompress_chunk(ztr, smp4_chunk)) {
    printf("Couldn't uncompresss SMP4 chunk\n");
    return -1;
  }
  
  if (smp4_chunk->dlength != nbases * 8 + 2) {
    printf("Trace and basecalls have different number of samples\n");
    return -1;
  }

  *indexes = smalloc(nbases * sizeof(int));
  if (NULL == type || type[0] == 'P') {
    data[0] = (uint8_t *) smp4_chunk->data + 2;
    data[1] = (uint8_t *) smp4_chunk->data + 2 + nbases * 2;
    data[2] = (uint8_t *) smp4_chunk->data + 2 + nbases * 4;
    data[3] = (uint8_t *) smp4_chunk->data + 2 + nbases * 6;
  } else {
    data[1] = (uint8_t *) smp4_chunk->data + 2;
    data[0] = (uint8_t *) smp4_chunk->data + 2 + nbases * 2;
    data[3] = (uint8_t *) smp4_chunk->data + 2 + nbases * 4;
    data[2] = (uint8_t *) smp4_chunk->data + 2 + nbases * 6;
  }

  for (i = 0; i < nbases; i++) {
    max    = ((uint16_t) data[0][i * 2] << 8) | data[0][i * 2 + 1];
    second = ((uint16_t) data[1][i * 2] << 8) | data[1][i * 2 + 1];
    if (second > max) {
      val = max; max = second; second = val;
      max_base = 1;
      second_base = 0;
    } else {
      max_base = 0;
      second_base = 1;
    }
    for (base = 2; base < 4; base++) {
      val = ((uint16_t) data[base][i * 2] << 8) | data[base][i * 2 + 1];
      if (val > max) {
	second = max;
	second_base = max_base;
	max = val;
	max_base = base;
      } else if (val > second) {
	second = val;
	second_base = base;
      }
    }
    (*indexes)[i] = second_base < max_base ? second_base : second_base - 1;
  }

  return 0;
}

int check_cnf4(ztr_t *ztr, ztr_chunk_t *chunk, Trace_reader *trc) {
  size_t  i = 0;
  int    *second_calls = NULL;
  int     res;
  int     base;
  char    expected;
  char    val;

  if (trc->seq_len * 4 + 1 != chunk->dlength) return -1;
  if (0 != memcmp(trc->qual, chunk->data + 1, trc->seq_len)) return -1;

  res = get_second_calls(ztr, trc->seq_len, &second_calls);
  if (res < 0) {
    return -1;
  }

  if (0 == res) {
    for (i = 0; i < trc->seq_len; i++) {
      for (base = 0; base < 3; base++) {
	expected = base == second_calls[i] ? -chunk->data[i + 1] : -40;
	if (chunk->data[trc->seq_len + 1 + i * 3 + base] != expected)
	  goto quickexit;
      }
    }
  } else {
    for (i = 0; i < trc->seq_len; i++) {
      for (base = 0; base < 3; base++) {
	val = chunk->data[trc->seq_len + 1 + i * 3 + base];
	if (val != -40 && val != -chunk->data[i + 1]) goto quickexit;
      }
    }
  }
 quickexit:
  if (NULL != second_calls) free(second_calls);
  if (i != trc->seq_len) return -1;
  return 0;
}

int check_trace_body(Settings *opts, Trace_reader *trc, srf_t *srf) {
  char         name[512];
  Spot_info    spot = { NULL, 0, 0, 0, 0, 0, 0 };
  ztr_t       *ztr = NULL;
  ztr_chunk_t *chunk;
  int          start_chunk = 0;
  int          i;

  if (-1 == construct_trace_name(srf->th.id_prefix,
				 (unsigned char *)srf->tb.read_id,
				 srf->tb.read_id_length,
				 name, sizeof(name))) {
    printf("Couldn't construct read name\n");
    return -1;
  }

  if (0 != decode_name(name, &spot)) {
    printf("Couldn't decode read name.\n");
    return -1;
  }

  if (opts->verbosity > 0) {
    printf("%s\n", name);
  }

  if (0 != read_trace_data(opts, trc, &spot)) {
    printf("Couldn't get trace data for %s\n", name);
    return -1;
  }

  mfseek(srf->mf, srf->mf_end, SEEK_SET);
  if (srf->tb.trace_size) {
    mfwrite(srf->tb.trace, 1, srf->tb.trace_size, srf->mf);
    free(srf->tb.trace);
    srf->tb.trace = NULL;
  }
  mftruncate(srf->mf, mftell(srf->mf));
  mfseek(srf->mf, srf->mf_pos, SEEK_SET);
  
  if (srf->ztr) {
    start_chunk = srf->ztr->nchunks;
    ztr = ztr_dup(srf->ztr);
  }

  if (NULL == partial_decode_ztr(srf, srf->mf, ztr)) {
    printf("Couldn't decode ZTR data for %s\n", name);
    return -1;
  }

  for (i = start_chunk; i < ztr->nchunks; i++) {
    chunk = &ztr->chunk[i];
    if (opts->verbosity > 1) {
      printf(" Chunk %d type %.4s\n", i,
	      (char *) &chunk->type);
    }
    switch (chunk->type) {
    case ZTR_TYPE_BASE:
      if (0 != uncompress_chunk(ztr, chunk)) {
	printf("Couldn't uncompress BASE chunk\n");
	return -1;
      }
      if (trc->seq_len + 1 != chunk->dlength
	  || 0 != memcmp(trc->seq, chunk->data + 1, trc->seq_len)) {
	printf("Sequence differs for %s\n", name);
	return -1;
      }
      break;
    case ZTR_TYPE_CNF1:
      if (0 != uncompress_chunk(ztr, chunk)) {
	printf("Couldn't uncompress CNF1 chunk\n");
	return -1;
      }
      if (trc->seq_len + 1 != chunk->dlength
	  || 0 != memcmp(trc->qual, chunk->data + 1, trc->seq_len)) {
	printf("CNF1 confidence values differ for %s\n", name);
	return -1;
      }
      break;
    case ZTR_TYPE_CNF4:
      if (0 != uncompress_chunk(ztr, chunk)) {
	printf("Couldn't uncompress CNF4 chunk\n");
	return -1;
      }
      if (0 != check_cnf4(ztr, chunk, trc)) {
	printf("CNF4 confidence values differ for %s\n", name);
	return -1;
      }
      break;
    case ZTR_TYPE_SMP4:
      if (0 != uncompress_chunk(ztr, chunk)) {
	printf("Couldn't uncompress SMP4 chunk\n");
	return -1;
      }
      if (0 != check_trace(ztr, chunk, trc)) {
	printf("Trace doesn't match for %s\n", name);
	return -1;
      }
      break;
    default:
      printf("Found unexpected chunk type in trace body for %s\n",
	      name);
      return -1;
    }
  }

  delete_ztr(ztr);

  return 0;
}

int read_srf(srf_t *srf, char *srf_name, Settings *opts) {
  int type;
  off_t pos;
  Previous_data seen = { NULL, 0, NULL, 0, NULL, 0 };
  Trace_reader  trc = { 0, 0, 0, NULL, NULL, 0, NULL, NULL };

  do {
    if ((pos = ftello(srf->fp)) < 0) {
      printf("Couldn't ftello() on %s: %s\n",
	      srf_name, strerror(errno));
      return -1;
    }
    switch (type = srf_next_block_type(srf)) {
    case -1:
      /* EOF */
      return 0;
      
    case SRFB_CONTAINER:
      printf("Found container header at %zd\n", pos);
      if (0 != srf_read_cont_hdr(srf, &srf->ch)) {
	printf("srf_read_cont_hdr failed for %s", srf_name);
	return -1;
      }
      break;

    case SRFB_XML:
      printf("Found xml block at %zd\n", pos);
      if (0 != srf_read_xml(srf, &srf->xml)) {
	printf("srf_read_xml failed for %s\n", srf_name);
	return -1;
      }
      break;

    case SRFB_TRACE_HEADER:
      printf("Found trace header at %zd\n", pos);
      if (0 != srf_read_trace_hdr(srf, &srf->th)) {
	printf("srf_read_trace_header failed for %s\n", srf_name);
	return -1;
      }
      if (0 != check_trace_header(opts, srf, &seen)) {
	printf("Bad trace header found in %s\n", srf_name);
	return -1;
      }
      break;

    case SRFB_TRACE_BODY:
      if (opts->verbosity > 0) {
	printf("Found trace body at %zd\n", pos);
      }
      if (0 != srf_read_trace_body(srf, &srf->tb, 0)) {
	printf("srf_read_trace_body failed for %s\n", srf_name);
	return -1;
      }
      if (0 != check_trace_body(opts, &trc, srf)) {
	printf("Bad trace body found in %s\n", srf_name);
	return -1;
      }
      if (srf->tb.trace) {
	free(srf->tb.trace);
	srf->tb.trace = NULL;
      }
      break;

    case SRFB_INDEX:
      if (0 != srf_read_index_hdr(srf, &srf->hdr, 1)) {
	printf("srf_read_index_hdr failed for %s\n", srf_name);
	return -1;
      }
      /* Skip the index body */
      if (0 != fseeko(srf->fp, pos + srf->hdr.size, SEEK_SET)) {
	printf("Couldn't fseeko() on %s: %s\n",
		srf_name, strerror(errno));
	return -1;
      }
      break;

    case SRFB_NULL_INDEX: {
      uint64_t idx = 0;
      
      printf("Found null index at %zd\n", pos);
      if (1 != fread(&idx, 8, 1, srf->fp)) {
	printf("Couldn't read null index from %s: %s\n",
		srf_name, strerror(errno));
	return -1;
      }
      if (0 != idx) {
	printf("Null index does not contain zeros in %s\n",
		srf_name);
      }

      break;
    }
    default:
      printf("Block of unknown type '%d' found in %s\n\n",
	      type, srf_name);
      return -1;
    }
  } while (1);
}

int check_srf(char *srf_name, Settings *opts) {
  FILE *fp = fopen(srf_name, "rb");
  srf_t *srf;
  int res;
  
  if (NULL == fp) {
    printf("Couldn't open %s: %s\n", srf_name, strerror(errno));
    return -1;
  }

  srf = srf_create(fp);
  if (NULL == srf) {
    printf("srf_create failed for %s\n", srf_name);
    fclose(fp);
    return -1;
  }

  res = read_srf(srf, srf_name, opts);

  srf_destroy(srf, 0);
  if (0 != fclose(fp)) {
    printf("Error closing %s: %s\n", srf_name, strerror(errno));
    return -1;
  }

  return res;
}

int main(int argc, char **argv) {
  Settings opts = { NULL, NULL, NULL, 0, 0 };
  int res;
  char *usage = "Usage: %s -b <bustard_dir> [-f] <srf_file>\n";
  int exit_val = EXIT_SUCCESS;
  int i;
  
  setvbuf(stdout, (char *)NULL, _IOLBF, 0);

  while ((res = getopt(argc, argv, "b:q:s:fv")) >= 0) {
    switch (res) {
    case 'b':
      opts.bustard_dir = optarg;
      break;
    case 'q':
      opts.qseq_dir = optarg;
      break;
    case 's':
      opts.qseq_suffix = optarg;
      break;
    case 'f':
      opts.filter_bad_reads = 1;
      break;
    case 'v':
      opts.verbosity++;
      break;
    default:
      die(usage, argv[0]);
    }
  }

  if (NULL == opts.bustard_dir || optind == argc) die(usage, argv[0]);

  if (NULL == opts.qseq_dir) opts.qseq_dir = opts.bustard_dir;
  if (NULL == opts.qseq_suffix) {
    opts.qseq_suffix = (0 == strcmp(opts.qseq_dir, opts.bustard_dir)
			? "qseq.txt"
			: "custom_qseq.txt");
  }

  for (i = optind; i < argc; i++) {
    if (check_srf(argv[i], &opts)) exit_val = EXIT_FAILURE;
  }

  return exit_val;
}
