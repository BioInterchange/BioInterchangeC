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

// Note: Needs to be at least 2, or otherwise, search operations will not overlap:
#define BI_FA_IDX_MUL 1024

// With the 80-character limit per line, this should be a safe bet about
// the maximum length of a sequence identifier (plus next line):
#define BI_FA_ID_MAXLEN 2048

#define GFF_FA_PFX_LEN 7

extern const char* GFF_C1;
extern const char* GFF_C2;
extern const char* GFF_C3;
extern const char* GFF_C4;
extern const char* GFF_C5;
extern const char* GFF_C6;
extern const char* GFF_C7;
extern const char* GFF_C8;
    
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
    
void gff_cbcks(gen_cbcks_t* cbcks);
    
off_t gff_fnd_fa(int fd, gen_fstat* stat, off_t mx);

ldoc_trie_t* gff_idx_fa(int fd, gen_fstat* stat, off_t mx);
    
char* gff_seq(int fd, off_t mx, ldoc_trie_t* idx, const char* id, off_t st, off_t en, bool rv);

char* gff_proc_cmt(char* ln, size_t lnlen);
    
void gff_proc_gbld(ldoc_nde_t* nde, char* val);
ldoc_nde_t* gff_proc_sregion(char* val, char** cid);
    
ldoc_struct_t gff_prgm_tpe(char* ky);
void gff_splt_attrs(ldoc_nde_t* ftr, ldoc_nde_t* usr, char* attrs);
    
ldoc_doc_t* gff_proc_ln(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt, gen_fstat* stat);
    
char* gff_proc_doc_ftr(ldoc_nde_t* ftr);
char* gff_proc_doc(ldoc_doc_t* doc);
    
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* defined(__biointerchange__gff__) */
