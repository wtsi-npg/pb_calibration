#define MAXNH       7

int *reverse_int(int *num, int n);
char *reverse_seq(char *seq);
char complement_base(char c);
char *complement_seq(char *seq);
const char *parse_next_int(const char *str, int *val, const char *sep);
char *append_char(char *cp, char c);
void checked_chdir(const char *dir);
char *alloc_getcwd(void);
char *get_real_path_name(const char*); 
char *get_command_line(int argc, char **argv);
int str2word(char *seq, int NH);
char *word2str(int word, int NH);
int int_cmp(const void *t1, const void *t2);
int int_sort(const void *t1, const void *t2);
int long_sort(const void *t1, const void *t2);
