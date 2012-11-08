#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <io_lib/hash_table.h>
#include "die.h"

HashTable *readSnpFile(char *snp_file)
{
    FILE *fp;
    HashTable *snp_hash;
    static const int line_size = 8192;
    char line[line_size];
    size_t count = 0;
    char last_chrom[100] = "";

    display("reading snp file %s\n", snp_file);

    fp = fopen(snp_file, "rb");
    if (NULL == fp) {
        die("ERROR: can't open known snp file %s: %s\n", snp_file, strerror(errno));
    }

    if (NULL == (snp_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS3))) {
        die("ERROR: creating snp hash table\n");
    }

    while (fgets(line, line_size, fp)) {
        char key[100];
        HashData hd;
        int bin, start, end;
        char chrom[100];

        if (4 != sscanf(line, "%d\t%s\t%d\t%d", &bin, chrom, &start, &end)) {
            die("ERROR: reading snp file\n%s\n", line);
        }

        /* N.B rod start is 0 based */
        snprintf(key, sizeof(key), "%s:%d", chrom, start);
        hd.i = 0;
        if (NULL == HashTableAdd(snp_hash, key, strlen(key), hd, NULL)) {
            die("ERROR: building snp hash table\n");
        }

        if (strcmp(chrom, last_chrom)) {
            strcpy(last_chrom, chrom);
            count = 0;
        }

        count++;
    }

    fclose(fp);

	return snp_hash;
}

