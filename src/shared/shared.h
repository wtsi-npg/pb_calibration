#define MAXNH      4
#define LEN_SUBST  2
#define NUM_SUBST  1<<(2*LEN_SUBST)
#define LEN_CNTXT  3
#define NUM_CNTXT  1<<(2*LEN_CNTXT)

int *reverse_int(int *num, int n);
char *reverse_seq(char *seq);
char complement_base(char c);
char *complement_seq(char *seq);
const char *parse_next_int(const char *str, int *val, const char *sep);
char *append_int(char *cp, int i);
char *append_char(char *cp, char c);
void checked_chdir(const char *dir);
char *alloc_getcwd(void);
char *get_real_path_name(const char*); 
char *get_command_line(int argc, char **argv);
int str2word(char *seq, int NH);
char *word2str(int word, int NH);
int tile_sort(const void *t1, const void *t2);
