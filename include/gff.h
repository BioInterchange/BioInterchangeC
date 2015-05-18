/*
 * Copyright (c) 2015 CODAMONO, Ontario, Canada
 *
 * This Source Code Form is subject to the terms of the accompanying
 * LICENSE.txt file, or, available via the following URLs:
 *
 *   As text: http://www.codamono.com/license/biointerchange-l1.txt
 *
 *   As PDF:  http://www.codamono.com/license/biointerchange-l1.pdf
 *
 *   As DOCX: http://www.codamono.com/license/biointerchange-l1.docx
 */

#ifndef __biointerchange__gff__
#define __biointerchange__gff__

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/queue.h>

#include <document.h>

#include "fio.h"
#include "gen.h"

#ifdef __cplusplus
extern "C" {
#endif
    
#define BI_NFOUND -1

// Number of memory pages that should reside in memory when reading a genomics file:
#define BI_GEN_PG_MUL 4096

// Note: Needs to be at least 2, or otherwise, search operations will not overlap:
#define BI_FA_IDX_MUL 1024

// With the 80-character limit per line, this should be a safe bet about
// the maximum length of a sequence identifier (plus next line):
#define BI_FA_ID_MAXLEN 2048

#define GFF_FA_PFX_LEN 7

extern const char* GFF_FA;
extern const char* GFF_FA_PFX;
    
typedef enum
{
    BI_GEN_HEADER,
    BI_GEN_FEATURE,
    BI_GEN_FASTA
} bi_gen_state;

typedef enum
{
    BI_FASTA_REF
} bi_idx_ntry;
    
typedef enum
{
    BI_NKW = 0,
    /**
     * String value.
     */
    BI_VAL,
    /**
     * Comma separated values.
     */
    BI_CSEP,
    /**
     * CIGAR string in GFF3 format (needs reformatting to match actual CIGAR specification).
     */
    BI_XCIG,
    BI_STRUCT
} bi_attr;
    
void gen_init();
    
off_t gff_fnd_fa(int fd, gen_fstat* stat, off_t mx);

ldoc_trie_t* gff_idx_fa(int fd, gen_fstat* stat, off_t mx);

void gff_rd(int fd, off_t mx, ldoc_trie_t* idx);
    
char* gff_seq(int fd, off_t mx, ldoc_trie_t* idx, const char* id, off_t st, off_t en, bool rv);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* defined(__biointerchange__gff__) */
