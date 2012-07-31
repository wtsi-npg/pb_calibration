#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <errno.h>
#include <dirent.h>
#include <stdint.h>
#include <assert.h>

#include <smalloc.h>
#include <aprintf.h>
#include <die.h>
#include <illumina_srf_util_zfp.h>

typedef enum {
  ot_1p3,
  ot_ipar,
  ot_rta,
  ot_ntypes,
} OutputType;

const static char * OutputTypeStrs[ot_ntypes] = {
  "1.3", "IPAR", "RTA"
};

typedef struct {
  OutputType  ot;
  int         ntiles;
  int         nlanes;
  int         nreads;
  int        *read_lengths;
  char       *out_dir;
  char       *machine_name;
  char       *run_id;
  int         tile_x;
  int         tile_y;
  int         density;
  int         prb;
} Settings;

typedef struct {
  struct tm   tm_now;
  char       *folder_dir;
  char       *config_dir;
  char       *data_dir;
  char       *signal_dir;    /* Firecrest */
  char       *basecall_dir;  /* Bustard   */
  char       *calibrate_dir; /* Gerald    */
  char       *matrix_dir;    /* Matrix files */
  char       *phasing_dir;   /* Phasing files */
  char      **lane_dirs;     /* RTA lane directories */
  char      **cycle_dirs;    /* RTA cycle directories */
  char       *signal_name;
  char       *basecall_name;
  char       *run_folder_date;
  char       *software_name;
  char       *bc_software;
  char       *software_version;
  int         cycles;
} State;

typedef struct {
  char **name;
  FILE **fp;
} Qseq_files;

typedef struct {
  char *int_name;
  zfp  *int_zfp;
  char *nse_name;
  zfp  *nse_zfp;
  char *sig_name;
  zfp  *sig_zfp;
  char *pos_name;
  FILE *pos_fp;
  size_t nspots;
  size_t spot;
  int16_t *tile_data;
  int16_t *proc_data;
} Trace_files;

void write_flowcell_id(Settings *opts, State *state) {
  char *name = aprintf("%s/FlowCellId.xml", state->config_dir);
  FILE *f = fopen(name, "w");

  if (NULL == f) {
    die("Couldn't open %s for writing : %s\n", name, strerror(errno));
  }
  fprintf(f,
	  "<?xml version=\"1.0\"?>\n"
	  "<FlowcellId />\n");
  if (0 != fclose(f)) die("Error writing to %s : %s\n", name, strerror(errno));
  free(name);
}

void write_firecrest_config(Settings *opts, State *state) {
  char *name = aprintf("%s/config.xml", state->signal_dir);
  FILE *f = fopen(name, "w");
  char start_time[128];
  int i;
  int cycle;
  
  if (NULL == f) {
    die("Couldn't open %s for writing : %s\n", name, strerror(errno));
  }

  strftime(start_time, sizeof(start_time), "%d-%m-%y %H:%M:%S %Z",
	   &state->tm_now);

  fprintf(f,
	  "<?xml version=\"1.0\"?>\n"
	  "<ImageAnalysis>\n");
  fprintf(f, "  <Run Name=\"%s\">\n", state->signal_name);
  fprintf(f, "    <Cycles First=\"1\" Last=\"%d\" Number=\"%d\" />\n",
	  state->cycles, state->cycles);
  fprintf(f,
	  "    <ImageParameters>\n"
	  "      <AutoOffsetFlag>0</AutoOffsetFlag>\n"
	  "      <AutoSizeFlag>0</AutoSizeFlag>\n"
	  "      <Fwhm>0</Fwhm>\n"
	  "      <RemappingDistance>0</RemappingDistance>\n"
	  "      <Threshold>0</Threshold>\n"
	  "    </ImageParameters>\n");
  fprintf(f,
	  "    <RunParameters>\n"
	  "      <AutoCycleFlag>0</AutoCycleFlag>\n"
	  "      <BasecallFlag>0</BasecallFlag>\n"
	  "      <Deblocked>0</Deblocked>\n"
	  "      <DebugFlag>0</DebugFlag>\n"
	  "      <FirstRunOnlyFlag>0</FirstRunOnlyFlag>\n");
  for (i = 0, cycle = 0; i < opts->nreads; cycle += opts->read_lengths[i++]) {
    fprintf(f,
	    "      <ImagingReads Index=\"%d\">\n"
	    "        <FirstCycle>%d</FirstCycle>\n"
	    "        <LastCycle>%d</LastCycle>\n"
	    "        <RunFolder>%s</RunFolder>\n"
	    "      </ImagingReads>\n",
	    i + 1, cycle + 1, cycle + opts->read_lengths[i],
	    state->folder_dir);
  }
  fprintf(f, "      <Instrument>%s</Instrument>\n", opts->machine_name);
  fprintf(f,
	  "      <IterativeMatrixFlag>0</IterativeMatrixFlag>\n"
	  "      <MakeFlag>1</MakeFlag>\n"
	  "      <MaxCycle>-1</MaxCycle>\n"
	  "      <MinCycle>-1</MinCycle>\n"
	  "      <Prb></Prb>\n");
  
  for (i = 0, cycle = 0; i < opts->nreads; cycle += opts->read_lengths[i++]) {
    fprintf(f,
	    "      <Reads Index=\"%d\">\n"
	    "        <FirstCycle>%d</FirstCycle>\n"
	    "        <LastCycle>%d</LastCycle>\n"
	    "        <RunFolder>%s</RunFolder>\n"
	    "      </Reads>\n",
	    i + 1, cycle + 1, cycle + opts->read_lengths[i],
	    state->folder_dir);
  }
  fprintf(f,
	  "      <RunFolder>%s</RunFolder>\n"
	  "      <RunFolderDate>%s</RunFolderDate>\n"
	  "      <RunFolderId>%s</RunFolderId>\n"
	  "      <SelectedFolders />\n"
	  "      <SelectedTiles />\n"
	  "      <Sig2></Sig2>\n",
	  state->folder_dir, state->run_folder_date, opts->run_id);
  fprintf(f, "    </RunParameters>\n");
  fprintf(f, "    <Software Name=\"%s\" Version=\"%s\" />\n",
	  state->software_name, state->software_version);
  fprintf(f, "    <TileSelection>\n");
  for (i = 0; i < opts->nlanes; i++) {
    fprintf(f, 
	    "      <Lane Index=\"%d\">\n"
	    "        <Sample>s</Sample>\n"
	    "        <TileRange Max=\"%d\" Min=\"1\" />\n"
	    "      </Lane>\n",
	    i + 1, opts->ntiles);
  }
  fprintf(f, "    </TileSelection>\n");
  fprintf(f,
	  "    <Time>\n"
	  "      <Start>%s</Start>\n"
	  "    </Time>\n",
	  start_time);
  fprintf(f, "    <UsedCycleFolders>C1.1");
  for (i = 1; i < state->cycles; i++) fprintf(f, ",C%d.1", i);
  fprintf(f, "</UsedCycleFolders>\n");
  fprintf(f, "    <User Name=\"%s\" />\n", "foo");
  fprintf(f, "  </Run>\n");
  fprintf(f, "</ImageAnalysis>\n");

  if (0 != fclose(f)) die("Error writing to %s : %s\n", name, strerror(errno));
  free(name);
}

void write_bustard_config(Settings *opts, State *state) {
  char *name = aprintf("%s/config.xml", state->basecall_dir);
  FILE *f = fopen(name, "w");
  char start_time[128];
  int i;
  int cycle;

  if (NULL == f) {
    die("Couldn't open %s for writing : %s\n", name, strerror(errno));
  }

  strftime(start_time, sizeof(start_time), "%d-%m-%y %H:%M:%S %Z",
	   &state->tm_now);

  fprintf(f,
	  "<?xml version=\"1.0\"?>\n"
	  "<BaseCallAnalysis>\n");
  fprintf(f, "  <Run Name=\"%s\">\n", state->basecall_name);
  fprintf(f,
	  "    <BaseCallParameters>\n"
	  "      <ChastityThreshold>0.600000</ChastityThreshold>\n");
  for (i = 0, cycle = 0; i < opts->nreads; cycle += opts->read_lengths[i++]) {
    fprintf(f,
	    "      <Matrix Path=\"\">\n"
	    "        <AutoFlag>0</AutoFlag>\n"
	    "        <AutoLane>0</AutoLane>\n"
	    "        <Cycle>2</Cycle>\n"
	    "        <CycleOffset>1</CycleOffset>\n"
	    "        <FirstCycle>%d</FirstCycle>\n"
	    "        <LastCycle>%d</LastCycle>\n"
	    "        <Read>%d</Read>\n"
	    "      </Matrix>\n",
	    cycle + 1, cycle + opts->read_lengths[i], i + 1);
  }
  fprintf(f, "      <MatrixElements />\n");
  for (i = 0, cycle = 0; i < opts->nreads; cycle += opts->read_lengths[i++]) {
    fprintf(f,
	    "      <Phasing Path=\"\">\n"
	    "        <AutoFlag>0</AutoFlag>\n"
	    "        <AutoLane>0</AutoLane>\n"
	    "        <Cycle>0</Cycle>\n"
	    "        <CycleOffset>0</CycleOffset>\n"
	    "        <FirstCycle>%d</FirstCycle>\n"
	    "        <LastCycle>%d</LastCycle>\n"
	    "        <Read>%d</Read>\n"
	    "      </Phasing>\n",
	    cycle + 1, cycle + opts->read_lengths[i], i + 1);
  }
  fprintf(f,
	  "      <PhasingRestarts />\n"
	  "      <PrbCompression>none</PrbCompression>\n"
	  "      <PureBases>25</PureBases>\n"
	  "      <QhgCompression>none</QhgCompression>\n"
	  "      <SeqCompression>none</SeqCompression>\n"
	  "      <Sig2Compression>gzip</Sig2Compression>\n"
	  "      <SmtFilter>failed-chastity</SmtFilter>\n"
	  "      <SmtRelation>le</SmtRelation>\n"
	  "      <SmtThreshold>1.000000</SmtThreshold>\n"
	  "    </BaseCallParameters>\n");
  fprintf(f,
	  "    <Cycles First=\"1\" Last=\"%d\" Number=\"%d\" />\n"
	  "    <Input Path=\"%s\" />\n",
	  state->cycles, state->cycles, state->signal_name);
  fprintf(f,
	  "    <RunParameters>\n"
	  "      <AutoCycleFlag>0</AutoCycleFlag>\n"
	  "      <BasecallFlag>0</BasecallFlag>\n"
	  "      <Compression>gzip</Compression>\n"
	  "      <CompressionSuffix>.gz</CompressionSuffix>\n"
	  "      <Deblocked>0</Deblocked>\n"
	  "      <DebugFlag>0</DebugFlag>\n"
	  "      <DirectionSuffix>.p</DirectionSuffix>\n"
	  "      <FirstRunOnlyFlag>0</FirstRunOnlyFlag>\n");
  for (i = 0, cycle = 0; i < opts->nreads; cycle += opts->read_lengths[i++]) {
    fprintf(f,
	    "      <ImagingReads Index=\"%d\">\n"
	    "        <FirstCycle>%d</FirstCycle>\n"
	    "        <LastCycle>%d</LastCycle>\n"
	    "        <RunFolder>%s</RunFolder>\n"
	    "      </ImagingReads>\n",
	    i + 1, cycle + 1, cycle + opts->read_lengths[i],
	    state->folder_dir);
  }
  fprintf(f,
	  "      <Instrument>%s</Instrument>\n"
	  "      <IterativeMatrixFlag>0</IterativeMatrixFlag>\n"
	  "      <MakeFlag>1</MakeFlag>\n"
	  "      <MaxCycle>-1</MaxCycle>\n"
	  "      <MinCycle>-1</MinCycle>\n"
	  "      <Prb></Prb>\n",
	  opts->machine_name);
  for (i = 0, cycle = 0; i < opts->nreads; cycle += opts->read_lengths[i++]) {
      fprintf(f,
	      "      <Reads Index=\"%d\">\n"
	      "        <FirstCycle>%d</FirstCycle>\n"
	      "        <LastCycle>%d</LastCycle>\n"
	      "        <RunFolder>%s</RunFolder>\n"
	      "      </Reads>\n",
	      i + 1, cycle + 1, cycle + opts->read_lengths[i],
	      state->folder_dir);
  }
  fprintf(f, 
	  "      <RunFolder>%s</RunFolder>\n"
	  "      <RunFolderDate>%s</RunFolderDate>\n"
	  "      <RunFolderId>%s</RunFolderId>\n"
	  "      <SelectedFolders />\n"
	  "      <SelectedTiles />\n"
	  "      <Sig2></Sig2>\n"
	  "    </RunParameters>\n",
	  state->folder_dir, state->run_folder_date, opts->run_id);
  fprintf(f, "    <Software Name=\"%s\" Version=\"%s\" />\n",
	  state->bc_software, state->software_version);
  fprintf(f, "    <TileSelection>\n");
  for (i = 0; i < opts->nlanes; i++) {
    fprintf(f,
	    "      <Lane Index=\"%d\">\n"
	    "        <Sample>s</Sample>\n"
	    "        <TileRange Max=\"%d\" Min=\"1\" />\n"
	    "      </Lane>\n",
	    i + 1, opts->ntiles);
  }
  fprintf(f,
	  "    </TileSelection>\n"
	  "    <Time>\n"
	  "      <Start>%s</Start>\n"
	  "    </Time>\n"
	  "    <User Name=\"foo\" />\n"
	  "  </Run>\n"
	  "</BaseCallAnalysis>\n",
	  start_time);
  
  if (0 != fclose(f)) die("Error writing to %s : %s\n", name, strerror(errno));
  free(name);	  
}

void write_bustard_summary(Settings *opts, State *state) {
  char *name = aprintf("%s/BustardSummary.xml", state->basecall_dir);
  FILE *f = fopen(name, "w");

  if (NULL == f) {
    die("Couldn't open %s for writing : %s\n", name, strerror(errno));
  }

  fprintf(f, "Testing testing 1...2...3...4...5   5...4...3...2...1\n");
  if (0 != fclose(f)) die("Error writing to %s : %s\n", name, strerror(errno));
  free(name);
}

void clean_dir(const char *path) {
  DIR           *dir;
  struct stat    st_buf;
  struct dirent *dirent;
  char          *buf = NULL;
  size_t         bufsz;
  size_t         pathsz;
  size_t         namesz;

  pathsz = strlen(path);
  bufsz = pathsz + 128;
  buf   = smalloc(bufsz);
  strncpy(buf, path, pathsz + 1);
  buf[pathsz++] = '/';

  dir = opendir(path);
  if (NULL == dir) die("Couldn't open directory %s : %s\n",
		       path, strerror(errno));
  while (NULL != (dirent = readdir(dir))) {
    namesz = strlen(dirent->d_name);
    if (bufsz < pathsz + namesz + 1) {
      bufsz = pathsz + namesz + 1;
      buf = srealloc(buf, bufsz);
    }
    strncpy(buf + pathsz, dirent->d_name, namesz + 1);
    
    if (0 != lstat(buf, &st_buf)) {
      die("Couldn't lookup %s : %s\n", buf, strerror(errno));
    }
    if (!S_ISDIR(st_buf.st_mode)) {
      if (0 != unlink(buf)) die("Couldn't remove %s : %s\n",
				buf, strerror(errno));
    }
  }

  if (0 != closedir(dir)) {
    die("Error on closing directory %s : %s\n", path, strerror(errno));
  }

  free(buf);
}

void make_dir(const char *path, const char *type, int clean) {
  struct stat st_buf;

  if (0 == stat(path, &st_buf)) {
    if (!S_ISDIR(st_buf.st_mode)) {
      die("Error: %s exists, but is not a directory.\n");
    }
    printf("Reused %s directory : %s\n", type, path);
    if (clean) clean_dir(path);
    return;
  }

  if (0 != mkdir(path, 0700)) {
    die("Couldn't mkdir %s : %s\n", path, strerror(errno));
  }

  printf("Made %s directory : %s\n", type, path);
}

void create_dirs(Settings *opts, State *state) {
  int i;

  make_dir(opts->out_dir, "root", 0);
  make_dir(state->folder_dir, "run folder", 1);
  make_dir(state->config_dir, "config", 1);
  make_dir(state->data_dir, "data", 1);
  make_dir(state->signal_dir, "signals", 1);
  make_dir(state->basecall_dir, "basecalls", 1);
  make_dir(state->matrix_dir, "matrix", 1);
  make_dir(state->phasing_dir, "phasing", 1);
  if (NULL != state->calibrate_dir) {
    make_dir(state->calibrate_dir, "recalibrated basecalls", 1);
  }
  if (NULL != state->lane_dirs) {
    for (i = 0; i < opts->nlanes; i++) {
      make_dir(state->lane_dirs[i], "RTA lane", 1);
    }
  }
  if (NULL != state->cycle_dirs) {
    for (i = 0; i < opts->nlanes * state->cycles; i++) {
      make_dir(state->cycle_dirs[i], "RTA cycle", 1);
    }
  }
}
 
void init_state(Settings *opts, State *state) {
  time_t     now;
  struct tm *tm;
  char       dd_mm_yyyy[12];
  int        i, j;

  state->cycles = 0;
  for (i = 0; i < opts->nreads; i++) {
    state->cycles += opts->read_lengths[i];
  }

  now = time(NULL);
  tm = localtime(&now);
  memcpy(&state->tm_now, tm, sizeof(struct tm));

  snprintf(dd_mm_yyyy, sizeof(dd_mm_yyyy), "%02d_%02d_%04d",
	   tm->tm_mday, tm->tm_mon + 1, tm->tm_year + 1900);

  state->run_folder_date = aprintf("%02d%02d%02d",
				   tm->tm_year % 100,
				   tm->tm_mon + 1,
				   tm->tm_mday);

  state->folder_dir = aprintf("%s/%s_%s_%s",
			      opts->out_dir,
			      state->run_folder_date,
			      opts->machine_name,
			      opts->run_id);
  state->config_dir = aprintf("%s/Config", state->folder_dir);
  state->data_dir = aprintf("%s/Data", state->folder_dir);

  switch (opts->ot) {
  case ot_1p3:
  case ot_ipar:
    state->software_name    = "Firecrest";
    state->bc_software      = "Bustard";
    state->software_version = ot_1p3 == opts->ot ? "1.3.4" : "1.4.0";
    state->signal_name = aprintf("C1-%d_Firecrest%s_%s_foo",
				 state->cycles,
				 state->software_version,
				 dd_mm_yyyy);
    state->basecall_name = aprintf("Bustard%s_%s_foo",
				   state->software_version,
				   dd_mm_yyyy);
    break;
  case ot_rta:
    state->software_name    = "RTA";
    state->bc_software      = "RTA";
    state->software_version = "1.4.15.0";
    state->signal_name         = aprintf("Intensities");
    state->basecall_name = aprintf("BaseCalls");
    break;
  default:
    die("This shouldn't happen\n");
    break;
  }

  state->signal_dir   = aprintf("%s/%s",
				state->data_dir, state->signal_name);
  state->basecall_dir = aprintf("%s/%s",
				state->signal_dir, state->basecall_name);
  state->matrix_dir = aprintf("%s/Matrix", state->basecall_dir);
  state->phasing_dir = aprintf("%s/Phasing", state->basecall_dir);
  if (opts->ot == ot_rta) {
    state->calibrate_dir = NULL;
  } else {
    state->calibrate_dir = aprintf("%s/GERALD_%s_foo",
				   state->basecall_dir, dd_mm_yyyy);
  }

  if (ot_rta == opts->ot) {
    state->lane_dirs = smalloc(opts->nlanes * sizeof(char *));
    state->cycle_dirs = smalloc(opts->nlanes * state->cycles * sizeof(char *));
    for (i = 0; i < opts->nlanes; i++) {
      state->lane_dirs[i] = aprintf("%s/L%03d", state->signal_dir, i + 1);
      for (j = 0; j < state->cycles; j++) {
	state->cycle_dirs[i * state->cycles + j]
	  = aprintf("%s/C%d.1", state->lane_dirs[i], j + 1);
      }
    }
  } else {
    state->lane_dirs  = NULL;
    state->cycle_dirs = NULL;
  }
}

void write_matrix_files(Settings *opts, State *state) {
  char *name;
  FILE *f;
  int l, r, c;

  static const char *fmt_rta =
    "%f\n%f\n%f\n%f\n"
    "%f\n%f\n%f\n%f\n"
    "%f\n%f\n%f\n%f\n"
    "%f\n%f\n%f\n%f\n";

  static const char *fmt_other =
    "%.2f %.2f %.2f %.2f\n"
    "%.2f %.2f %.2f %.2f\n"
    "%.2f %.2f %.2f %.2f\n"
    "%.2f %.2f %.2f %.2f\n";

  for (l = 0; l < opts->nlanes; l++) {
    for (r = 0, c = 0; r < opts->nreads; c += opts->read_lengths[r++]) {
      if (ot_rta == opts->ot) {
	name = aprintf("%s/s_%d_%d_matrix.txt",
		       state->matrix_dir, l + 1, r + 1);
      } else {
	name = aprintf("%s/s_%d_%02d_matrix.txt",
		       state->matrix_dir, l + 1, c + 2);
      }
      f = fopen(name, "w");
      if (NULL == f) {
	die("Couldn't open %s for writing: %s\n", name, strerror(errno));
      }
      printf("Writing %s\n", name);
      if (ot_rta != opts->ot) {
	fprintf(f,
		"# Auto-generated frequency response matrix\n"
		"> A\n"
		"> C\n"
		"> G\n"
		"> T\n");
      }
      fprintf(f, ot_rta == opts->ot ? fmt_rta : fmt_other,
	      0.0, 1.0, 0.0, 0.0,
	      1.0, 0.0, 0.0, 0.0,
	      0.0, 0.0, 0.0, 1.0,
	      0.0, 0.0, 1.0, 0.0);
      if (0 != fclose(f)) {
	die("Error closing %s: %s\n", name, strerror(errno));
      }
      free(name);
    }
  }
}

void write_phasing_files(Settings *opts, State *state) {
  char *name;
  FILE *f;
  int l, r, c;

  static const char *fmt_rta = "%f\n%f\n";
  static const char *fmt_other =
    "<Parameters>\n"
    "  <Phasing>%f</Phasing>\n"
    "  <Prephasing>%f</Prephasing>\n"
    "</Parameters>\n";

  for (l = 0; l < opts->nlanes; l++) {
    for (r = 0, c = 0; r < opts->nreads; c += opts->read_lengths[r++]) {
      if (ot_rta == opts->ot) {
	name = aprintf("%s/s_%d_%d_phasing.txt",
		       state->phasing_dir, l + 1, r + 1);
      } else {
	name = aprintf("%s/s_%d_%02d_phasing.xml",
		       state->phasing_dir, l + 1, c + 1);
      }
      f = fopen(name, "w");
      if (NULL == f) {
	die("Couldn't open %s for writing: %s\n", name, strerror(errno));
      }
      printf("Writing %s\n", name);
      
      fprintf(f, ot_rta == opts->ot ? fmt_rta : fmt_other, 0.0, 0.0);
      if (0 != fclose(f)) {
	die("Error closing %s: %s\n", name, strerror(errno));
      }
      free(name);
    }
  }
}

Qseq_files * open_qseq_files(char *path, char *suffix,
			     int lane, int tile, int nreads) {
  Qseq_files *qseqs = smalloc(sizeof(Qseq_files));
  int i;

  qseqs->name = smalloc(nreads * sizeof(char *));
  qseqs->fp   = smalloc(nreads * sizeof(FILE *));

  for (i = 0; i < nreads; i++) {
    qseqs->name[i] = aprintf("%s/s_%d_%d_%04d_%s",
			     path, lane + 1, i + 1, tile + 1, suffix);
    qseqs->fp[i] = fopen(qseqs->name[i], "w");
    if (NULL == qseqs->fp[i]) {
      die("Couldn't open %s for writing: %s\n",
	  qseqs->name[i], strerror(errno));
    }
  }

  printf("Writing %s/s_%d_*_%04d_%s\n", path, lane + 1, tile + 1, suffix);

  return qseqs;
}

void close_qseq_files(Qseq_files *qseqs, int nreads) {
  int i;

  for (i = 0; i < nreads; i++) {
    if (0 != fclose(qseqs->fp[i])) {
      die("Error closing %s: %s\n", qseqs->name[i], strerror(errno));
    }
    free(qseqs->name[i]);
  }
  free(qseqs->name);
  free(qseqs->fp);
  free(qseqs);
}

inline static int rnd2(void) {
  static long num = 0;
  static long bits = 0;
  int res;

  if (bits < 3) {
    num = random();
    bits = RAND_MAX;
  }

  res = num & 3;
  num  >>= 2;
  bits >>= 2;

  return res;
}

Trace_files * open_trace_files(int lane, int tile,
			       Settings *opts, State *state) {
  Trace_files *tf;

  tf = smalloc(sizeof(Trace_files));
  
  tf->int_name = tf->nse_name = tf->sig_name = tf->pos_name = NULL;
  tf->int_zfp = tf->nse_zfp = tf->sig_zfp = NULL;
  tf->pos_fp = NULL;
  tf->spot = tf->nspots = 0;
  tf->tile_data = NULL;
  tf->proc_data = NULL;
  
  switch (opts->ot) {
  case ot_1p3:
    tf->int_name = aprintf("%s/s_%d_%04d_int.txt",
			   state->signal_dir, lane + 1, tile + 1);
    tf->nse_name = aprintf("%s/s_%d_%04d_nse.txt",
			   state->signal_dir, lane + 1, tile + 1);
    tf->sig_name = aprintf("%s/s_%d_%04d_sig2.txt",
			   state->basecall_dir, lane + 1, tile + 1);
    break;
  case ot_ipar:
    tf->int_name = aprintf("%s/s_%d_%04d_int.txt.p",
			   state->signal_dir, lane + 1, tile + 1);
    tf->nse_name = aprintf("%s/s_%d_%04d_nse.txt.p",
			   state->signal_dir, lane + 1, tile + 1);
    tf->sig_name = aprintf("%s/s_%d_%04d_sig2.txt",
			   state->basecall_dir, lane + 1, tile + 1);
    break;
  default:
    break;
  }
  
  tf->pos_name = aprintf("%s/s_%d_%04d_pos.txt",
			 state->signal_dir, lane + 1, tile + 1);

  if (tf->int_name) {
    tf->int_zfp = zfopen_write(tf->int_name, COMP_GZIP);
    if (NULL == tf->int_zfp) {
      die("Couldn't open %s for writing: %s\n", tf->int_name, strerror(errno));
    }
  }
  if (tf->nse_name) {
    tf->nse_zfp = zfopen_write(tf->nse_name, COMP_GZIP);
    if (NULL == tf->nse_zfp) {
      die("Couldn't open %s for writing: %s\n", tf->nse_name, strerror(errno));
    }
  }
  if (tf->sig_name) {
    tf->sig_zfp = zfopen_write(tf->sig_name, COMP_GZIP);
    if (NULL == tf->sig_zfp) {
      die("Couldn't open %s for writing: %s\n", tf->sig_name, strerror(errno));
    }
  }
  if (tf->pos_name) {
    tf->pos_fp = fopen(tf->pos_name, "w");
    if (NULL == tf->pos_fp) {
      die("Couldn't open %s for writing: %s\n", tf->pos_name, strerror(errno));
    }
  }
  
  return tf;
}

void write_signal(Trace_files *tf, int lane, int tile, int x, int y,
		  Settings *opts, State *state, int16_t (*signal)[4]) {
  static char *line = NULL;
  static size_t len = 1024;
  int16_t *sigp;
  size_t p = 0, q = 0;
  int i;

  if (NULL == line) {
    line = smalloc(len);
  }

  switch (opts->ot) {
  case ot_1p3:
    p = q = snprintf(line, len, "%d\t%d\t%d\t%d",
		     lane + 1, tile + 1, x, y);
    for (i = 0; i < state->cycles; i++) {
      if (len - p < 256) {
	len *= 2;
	line = srealloc(line, len);
      }
      p += snprintf(line + p, len - p, "\t%d.0 %d.0 %d.0 %d.0",
		    signal[i][1], signal[i][0], signal[i][3], signal[i][2]);
    }
    line[p++] = '\n';
    line[p++] = '\0';
    zfputs(line, tf->int_zfp);

    p = q;
    for (i = 0; i < state->cycles; i++) {
      if (len - p < 256) {
	len *= 2;
	line = srealloc(line, len);
      }
      p += snprintf(line + p, len - p, "\t%d.0 %d.0 %d.0 %d.0",
		    1, -1, 1, -1);
    }
    line[p++] = '\n';
    line[p++] = '\0';
    zfputs(line, tf->nse_zfp);

    /* NB: falls through here */
  case ot_ipar:

    if (ot_ipar == opts->ot) {
      p = q = snprintf(line, len, "%d\t%d\t%d\t%d",
		     lane + 1, tile + 1, x, y);
    } else {
      p = q;
    }

    for (i = 0; i < state->cycles; i++) {
      if (len - p < 256) {
	len *= 2;
	line = srealloc(line, len);
      }
      p += snprintf(line + p, len - p, "\t%d.0 %d.0 %d.0 %d.0",
		    signal[i][0], signal[i][1], signal[i][2], signal[i][3]);
    }
    line[p++] = '\n';
    line[p++] = '\0';
    zfputs(line, tf->sig_zfp);

    break;

  default:
    break;
  }

  /* The IPAR int and nse files and all RTA files have formats that are
     difficult to write as we go.  In these cases we save the data
     and then write it all out when we get to the end of the tile */

  switch (opts->ot) {
  case ot_rta:
    assert(tf->spot < tf->nspots);
    sigp = tf->proc_data + tf->spot * state->cycles * 4;
    for (i = 0; i < state->cycles; i++) {
      *sigp++ = signal[i][0];
      *sigp++ = signal[i][1];
      *sigp++ = signal[i][2];
      *sigp++ = signal[i][3];
    }

    /* Note falls through */
  case ot_ipar:
    assert(tf->spot < tf->nspots);
    sigp = tf->tile_data + tf->spot * state->cycles * 4;
    for (i = 0; i < state->cycles; i++) {
      *sigp++ = signal[i][1];
      *sigp++ = signal[i][0];
      *sigp++ = signal[i][3];
      *sigp++ = signal[i][2];
    }
    
    tf->spot++;
    break;

  default:
    break;
  }
}

void close_trace_files(Trace_files *tf, Settings *opts, State *state) {
  char line[1024];
  size_t c, i;
  int16_t *sigp;

  if (ot_ipar == opts->ot) {
    snprintf(line, sizeof(line), "#CH4:OBJ%zd\n", tf->spot);
    zfputs(line, tf->int_zfp);
    zfputs(line, tf->nse_zfp);

    for (c = 0; c < state->cycles; c++) {
      for (i = 0; i < tf->spot; i++) {
	sigp = tf->tile_data + i * state->cycles * 4 + c * 4;
	
	snprintf(line, sizeof(line), "%d\t%d\t%d\t%d\n",
		 sigp[0], sigp[1], sigp[2], sigp[3]);
	zfputs(line, tf->int_zfp);
	zfputs("1\t-1\t1\t-1\n", tf->nse_zfp);
      }
      snprintf(line, sizeof(line), "#END CYCLE %zd\n", c);
      zfputs(line, tf->int_zfp);
      zfputs(line, tf->nse_zfp);
    }
  }

  if (NULL != tf->int_zfp) {
    if (0 != zfclose(tf->int_zfp)) {
      die("Error closing %s: %s\n", tf->int_name, strerror(errno));
    }
  }
  if (NULL != tf->int_name) free(tf->int_name);
  if (NULL != tf->nse_zfp) {
    if (0 != zfclose(tf->nse_zfp)) {
      die("Error closing %s: %s\n", tf->nse_name, strerror(errno));
    }
  }
  if (NULL != tf->nse_name) free(tf->nse_name);
  if (NULL != tf->sig_zfp) {
    if (0 != zfclose(tf->sig_zfp)) {
      die("Error closing %s: %s\n", tf->sig_name, strerror(errno));
    }
  }
  if (NULL != tf->sig_name) free(tf->sig_name);
  if (NULL != tf->pos_fp) {
    if (0 != fclose(tf->pos_fp)) {
      die("Error closing %s: %s\n", tf->pos_name, strerror(errno));
    }
  }
  if (NULL != tf->pos_name) free(tf->pos_name);
  
  if (NULL != tf->tile_data) free(tf->tile_data);
  if (NULL != tf->proc_data) free(tf->proc_data);
  free(tf);
}

void write_cif_files(int lane, int tile, Trace_files *tf,
		     Settings *opts, State *state) {
  uint8_t cif_header[13];
  uint8_t cnf_header[13];
  char *cif_name = NULL;
  char *cnf_name = NULL;
  char *dif_name = NULL;
  uint8_t *sig_data;
  FILE *cif;
  FILE *cnf;
  FILE *dif;
  int16_t *tdp;
  int   lane_idx = lane * state->cycles;
  int   cycle;
  size_t base;
  size_t s;

  memcpy(cif_header, "CIF\1\2\0\0\1\0", 9);
  cif_header[ 9] =  tf->spot        & 0xff;
  cif_header[10] = (tf->spot >>  8) & 0xff;
  cif_header[11] = (tf->spot >> 16) & 0xff;
  cif_header[12] = (tf->spot >> 24) & 0xff;
  memcpy(cnf_header, cif_header, sizeof(cnf_header));
  cnf_header[4] = 1;

  sig_data = smalloc(tf->spot * 2);

  for (cycle = 0; cycle < state->cycles; cycle++) {
    cif_name = aprintf("%s/s_%d_%d.cif",
		       state->cycle_dirs[lane_idx + cycle],
		       lane + 1, tile + 1);
    cnf_name = aprintf("%s/s_%d_%d.cnf",
		       state->cycle_dirs[lane_idx + cycle],
		       lane + 1, tile + 1);
    dif_name = aprintf("%s/s_%d_%d.dif",
		       state->cycle_dirs[lane_idx + cycle],
		       lane + 1, tile + 1);
    cif = fopen(cif_name, "wb");
    if (NULL == cif) {
      die("Couldn't open %s for writing: %s\n", cif_name, strerror(errno));
    }
    cnf = fopen(cnf_name, "wb");
    if (NULL == cnf) {
      die("Couldn't open %s for writing: %s\n", cnf_name, strerror(errno));
    }
    dif = fopen(dif_name, "wb");
    if (NULL == dif) {
      die("Couldn't open %s for writing: %s\n", dif_name, strerror(errno));
    }

    cif_header[5] = cnf_header[5] =  (cycle + 1)       & 0xff;
    cif_header[6] = cnf_header[6] = ((cycle + 1) >> 8) & 0xff;
    
    if (1 != fwrite(cif_header, sizeof(cif_header), 1, cif)) goto nowrite_cif;
    if (1 != fwrite(cnf_header, sizeof(cnf_header), 1, cnf)) goto nowrite_cnf;
    if (1 != fwrite(cif_header, sizeof(cif_header), 1, dif)) goto nowrite_dif;
    
    for (base = 0; base < 4; base++) {
      tdp = tf->tile_data + cycle * 4 + base;
      for (s = 0; s < tf->spot; s++, tdp += state->cycles * 4) {
	sig_data[s * 2]     =  (*tdp)       & 0xff;
	sig_data[s * 2 + 1] = ((*tdp) >> 8) & 0xff;
      }
      if (tf->spot != fwrite(sig_data, 2, tf->spot, cif)) goto nowrite_cif;

      tdp = tf->proc_data + cycle * 4 + base;
      for (s = 0; s < tf->spot; s++, tdp += state->cycles * 4) {
	sig_data[s * 2]     =  (*tdp)       & 0xff;
	sig_data[s * 2 + 1] = ((*tdp) >> 8) & 0xff;
      }
      if (tf->spot != fwrite(sig_data, 2, tf->spot, dif)) goto nowrite_dif;

      memset(sig_data, (base & 1) == 0 ? 1 : -1, tf->spot);
      if (tf->spot != fwrite(sig_data, 1, tf->spot, cnf)) goto nowrite_cnf;
    }
    if (0 != fclose(cif)) {
    nowrite_cif:
      die("Error writing to %s: %s\n", cif_name, strerror(errno));
    }
    if (0 != fclose(cnf)) {
    nowrite_cnf:
      die("Error writing to %s: %s\n", cnf_name, strerror(errno));
    }
    if (0 != fclose(dif)) {
    nowrite_dif:
      die("Error writing to %s: %s\n", dif_name, strerror(errno));
    }
  }
  free(sig_data);
  if (NULL != cif_name) free(cif_name);
  if (NULL != cnf_name) free(cnf_name);
  if (NULL != dif_name) free(dif_name);
}

void write_tile_data(int lane, int tile, Settings *opts, State *state) {
  Qseq_files *bc_qseqs  = NULL;
  Qseq_files *cal_qseqs = NULL;
  Trace_files *trace_files = NULL;
  unsigned char  bases[4] = { 0, 1, 2, 3 };
  int16_t (*signal)[4];
  int16_t (*prb)[4] = NULL;
  char *seq;
  char *qual;
  char *prb_line = NULL;
  char *prb_name = NULL;
  zfp  *prb_zfp  = NULL;
  size_t prb_len = state->cycles * 64;
  size_t prb_used;
  int x, y, r, c, b, j;
  int passed;
  char t;
  
  signal   = smalloc(state->cycles * sizeof(*signal));
  trace_files = open_trace_files(lane, tile, opts, state);

  for (r = 0, c = 0; r < opts->nreads; r++) {
    if (c < opts->read_lengths[r]) c = opts->read_lengths[r];
  }
  seq  = smalloc(c + 1);
  qual = smalloc(c + 1);

  bc_qseqs = open_qseq_files(state->basecall_dir, "qseq.txt",
			     lane, tile, opts->nreads);
  if (state->calibrate_dir) {
    cal_qseqs = open_qseq_files(state->calibrate_dir, "custom_qseq.txt",
				lane, tile, opts->nreads);
  }

  if (opts->prb) {
    prb      = smalloc(state->cycles * sizeof(*prb));
    prb_line = smalloc(prb_len);
    prb_name = aprintf("%s/s_%d_%04d_prb.txt",
		       state->basecall_dir, lane + 1, tile + 1);
    prb_zfp = zfopen_write(prb_name, COMP_GZIP);
    if (NULL == prb_zfp) {
      die("Couldn't open %s for writing: %s\n", prb_name, strerror(errno));
    }
  }

  if (opts->ot != ot_1p3) {
    trace_files->nspots = opts->tile_x * opts->tile_y / opts->density + 1;
    trace_files->tile_data = smalloc(trace_files->nspots * state->cycles * 4
				     * sizeof(int16_t));
    if (opts->ot == ot_rta) {
      trace_files->proc_data = smalloc(trace_files->nspots * state->cycles * 4
				       * sizeof(int16_t));
    }
  }

  for (y = 0, x = 0; y < opts->tile_y; y++) {
    for (; x < opts->tile_x; x += opts->density) {

      if (trace_files->pos_fp) {
	fprintf(trace_files->pos_fp, "%7.2f %7.2f\n", (double) x, (double) y);
      }

      for (r = 0, b = 0; r < opts->nreads; r++) {
	for (c = 0; c < opts->read_lengths[r]; c++, b++) {
	  /* Randomly shuffle bases */
	  j = rnd2(); t = bases[j]; bases[j] = bases[3]; bases[3] = t;
	  do { j = rnd2(); } while (j > 2);
	  t = bases[j]; bases[j] = bases[2]; bases[2] = t;
	  j = rnd2() & 1; t = bases[j]; bases[j] = bases[1]; bases[1] = t;

	  /* Make up some signal vals */
	  signal[b][bases[0]] = (random() & 0x7ff) + 0x100;
	  signal[b][bases[1]] = signal[b][bases[0]] / (rnd2() + 2);
	  signal[b][bases[2]] = signal[b][bases[1]] / (rnd2() + 2);
	  signal[b][bases[3]] = signal[b][bases[2]] / (rnd2() + 2);

	  if (opts->prb) {
	    prb[b][bases[0]] = (signal[b][bases[0]] - 0x100) >> 6;
	    prb[b][bases[1]] = -prb[b][bases[0]];
	    prb[b][bases[2]] = prb[b][bases[3]] = -40;
	  }

	  seq[c] = "ACGT"[bases[0]];
	  qual[c] = ((signal[b][bases[0]] - 0x100) >> 6) + 0x40;

	}
	seq[c] = qual[c] = '\0';

	passed = (random() & 0xff) > 8;
	fprintf(bc_qseqs->fp[r],
		"%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\n",
		opts->machine_name, opts->run_id, lane + 1, tile + 1,
		x, y, 0, r + 1, seq, qual, passed ? 1 : 0);
	if (cal_qseqs) {
	  for (c = 0; c < opts->read_lengths[r]; c++) {
	    qual[c]++;
	  }
	  fprintf(cal_qseqs->fp[r],
		  "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\n",
		  opts->machine_name, opts->run_id, lane + 1, tile + 1,
		  x, y, 0, r + 1, seq, qual, passed ? 1 : 0);
	}
      }

      if (opts->prb) {
	prb_used = 0;
	for (b = 0; b < state->cycles; b++) {
	  int n = snprintf(prb_line + prb_used, prb_len - prb_used,
			   "%s%d %d %d %d",
			   b > 0 ? "\t" : "",
			   prb[b][0], prb[b][1], prb[b][2], prb[b][3]);
	  assert(n > 0);
	  prb_used += n;
	}
	assert(prb_used < prb_len - 2);
	prb_line[prb_used++] = '\n';
	prb_line[prb_used++] = '\0';
	zfputs(prb_line, prb_zfp);
      }

      write_signal(trace_files, lane, tile, x, y, opts, state, signal);
    }
    x -= opts->tile_x;
  }

  if (ot_rta == opts->ot) {
    write_cif_files(lane, tile, trace_files, opts, state);
  }
  close_trace_files(trace_files, opts, state);
  close_qseq_files(bc_qseqs, opts->nreads);
  if (cal_qseqs) close_qseq_files(cal_qseqs, opts->nreads);
  if (NULL != prb_zfp && 0 != zfclose(prb_zfp)) {
    die("Error closing %s: %s\n", prb_name, strerror(errno));
  }
  if (NULL != prb_name) free(prb_name);
  
  free(signal);
  free(seq);
  free(qual);
  if (NULL != prb)      free(prb);
  if (NULL != prb_line) free(prb_line);
}

void write_lane_data(int lane, Settings *opts, State *state) {
  int tile;

  for (tile = 0; tile < opts->ntiles; tile++) {
    write_tile_data(lane, tile, opts, state);
  }
}

void write_data(Settings *opts, State *state) {
  int lane;

  for (lane = 0; lane < opts->nlanes; lane++) {
    write_lane_data(lane, opts, state);
  }
}

int main(int argc, char **argv) {
  Settings opts = {
    ot_1p3, 10, 8, 0, NULL, NULL, "XX1", "9876", 1000, 1000, 105, 0
  };
  State    state;
  int res;
  int i;
  char *c;
  char *usage =
    "Usage: %s -d <output_dir> -r <read_length> [-r <read_length> ...]\n";
  
  memset(&state, 0, sizeof(state));

  while ((res = getopt(argc, argv, "m:l:t:r:d:p")) >= 0) {
    switch (res) {
    case 'm':
      if (0 == strcmp(optarg, "1.3")) {
	opts.ot = ot_1p3;
      } else if (0 == strcmp(optarg, "ipar")) {
	opts.ot = ot_ipar;
      } else if (0 == strcmp(optarg, "rta")) {
	opts.ot = ot_rta;
      } else {
	fprintf(stderr, "Unknown mode: %s\n", optarg);
	return EXIT_FAILURE;
      }
      break;
    case 'l':
      opts.nlanes = atoi(optarg);
      break;
    case 't':
      opts.ntiles = atoi(optarg);
      break;
    case 'r':
      opts.nreads++;
      opts.read_lengths = srealloc(opts.read_lengths,
				   opts.nreads * sizeof(int));
      opts.read_lengths[opts.nreads - 1] = atoi(optarg);
      break;
    case 'd':
      opts.out_dir = optarg;
      break;
    case 'p':
      opts.prb = 1;
      break;
    default:
      die(usage, argv[0]);
    }
  }

  if (NULL == opts.out_dir || '\0' == *opts.out_dir) {
    fprintf(stderr, "No output directory specified!\n");
    die(usage, argv[0]);
  }
  
  c = opts.out_dir + strlen(opts.out_dir) - 1;
  while (c > opts.out_dir && *c == '/') {
    *c-- = '\0';
  }

  if (0 == strcmp(opts.out_dir, "/")) {
    die("Output directory '%s' specified.  Did you really want to do that?",
	opts.out_dir);
  }

  printf("Output to : %s\n", opts.out_dir);
  printf("Mode      : %s\n", OutputTypeStrs[opts.ot]);
  printf("No. lanes : %3d\n", opts.nlanes);
  printf("No. tiles : %3d\n", opts.ntiles);
  printf("No. reads : %3d\n", opts.nreads);
  for (i = 0; i < opts.nreads; i++) {
    printf("  Length of read %3d : %3d\n", i + 1, opts.read_lengths[i]);
  }
  printf("Machine name : %s\n", opts.machine_name);
  printf("Run Id       : %s\n", opts.run_id);
  printf("\n");

  init_state(&opts, &state);
  create_dirs(&opts, &state);
  write_flowcell_id(&opts, &state);
  write_firecrest_config(&opts, &state);
  write_bustard_config(&opts, &state);
  write_bustard_summary(&opts, &state);
  write_matrix_files(&opts, &state);
  write_phasing_files(&opts, &state);
  write_data(&opts, &state);

  return EXIT_SUCCESS;
}
