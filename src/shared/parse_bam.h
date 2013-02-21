#define BASE_ALIGN      (1<<0)
#define BASE_MISMATCH   (1<<1)
#define BASE_INSERTION  (1<<2)
#define BASE_DELETION   (1<<3)
#define BASE_SOFT_CLIP  (1<<4)
#define BASE_KNOWN_SNP  (1<<5)

int parse_bam_readinfo(samfile_t *fp, 
                         bam1_t *bam, 
                            int *bam_lane, 
                            int *bam_tile, 
                            int *bam_x, 
                            int *bam_y, 
                            int *bam_read, 
                         size_t *bam_offset);

int parse_bam_alignments(samfile_t *fp, 
                            bam1_t *bam, 
                              char *read_seq, 
                               int *read_qual, 
                              char *read_ref, 
                               int *read_mismatch, 
                         const int read_buff_size,
                         HashTable *snp_hash);

void bam_header_add_pg(       char *id,
                              char *pn,
                              char *ds,
                              char *cl,
                      bam_header_t *bam_header);


