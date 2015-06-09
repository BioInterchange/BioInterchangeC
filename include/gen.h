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

#ifndef biointerchange_gen_h
#define biointerchange_gen_h

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <document.h>

#include "ext-python.h"
#include "fio.h"

#ifdef __cplusplus
extern "C" {
#endif
    
// Number of memory pages that should reside in memory when reading a genomics file:
#define BI_GEN_PG_MUL 4096

extern const char* JSONLD_CTX;
extern const char* JSONLD_GFF3;
extern const char* JSONLD_GTF;
extern const char* JSONLD_GVF;
extern const char* JSONLD_VCF;

extern const char* GEN_ATTRS;
extern const char* GEN_BUILD;
extern const char* GEN_COMMENT;
extern const char* GEN_END;
extern const char* GEN_LOCUS;
extern const char* GEN_REFERENCE;
extern const char* GEN_SEQUENCE;
extern const char* GEN_START;
extern const char* GEN_SOURCE;
extern const char* GEN_VARIANTS;

extern const char* GEN_EMPTY;
extern const char* GEN_NULL;
extern const char* GEN_TRUE;
extern const char* GEN_UNKNOWN;

extern ldoc_vis_nde_ord_t* json_vis_nde;;
extern ldoc_vis_ent_t* json_vis_ent;
    
// Letters in the alphabet; restriction for labelling alleles (A, B, C, etc.):
#define GEN_MAX_ALT       26
// Two characters needed to encode for the value '26':
#define GEN_MAX_ALT_CHARS  2
    
// Characters to step over when encoding for allele labels:
// Note -- assumes VCF_MAX_ALT <= 26, so that only one character is
//         required to encode for the allele.
#define GEN_STEP 3
    
extern char GEN_ALLELE[GEN_MAX_ALT * 2];
extern char GEN_ALLELES[GEN_STEP * (((GEN_MAX_ALT * (GEN_MAX_ALT + 1) / 2) + GEN_MAX_ALT + 1))];
    
typedef enum
{
    /**
     * Unknown keyword ("no keyword").
     */
    BI_NKW = 0,
    /**
     * Ignore. This key/value pair has been processed already.
     */
    BI_IGN,
    /**
     * String value.
     */
    BI_VAL,
    /**
     * Number value.
     */
    BI_NUM,
    /**
     * Comma separated values.
     */
    BI_CSEP,
    /**
     * Comma separated information about variants.
     */
    BI_CSEPVAR,
    /**
     * Information about the reference sequence.
     */
    BI_REFSEQ,
    /**
     * CIGAR string in GFF3 format (needs reformatting to match actual CIGAR specification).
     */
    BI_XCIG,
    BI_STRUCT
} bi_attr;
    
typedef enum
{
    GEN_FMT_GFF3,
    GEN_FMT_GTF,
    GEN_FMT_GVF,
    GEN_FMT_VCF
} gen_filetype_t;

typedef enum
{
    GEN_CTPE_SETUP,
    GEN_CTPE_CLEANUP,
    GEN_CTPE_PROCESS
} gen_ctpe_t;
    
typedef struct gen_fstat
{
    /**
     * Number of meta-information lines.
     */
    uint32_t meta;
    /**
     * Number of feature lines.
     */
    uint32_t ftrs;
    /**
     * Number of comment lines.
     */
    uint32_t comms;
    /**
     * Indicates whether the meta information was filtered.
     */
    bool cbk_meta_fltr;
    /**
     * Number of filtered feature lines.
     */
    uint32_t cbk_ftrs_fltr;
    /**
     * Offset to which the statistics are valid.
     */
    off_t off;
} gen_fstat;

typedef struct gen_fa_ntry_t
{
    off_t off;
    off_t seq;
    size_t llen;
    size_t lbrk;
} gen_fa_ntry_t;

/**
 * @brief Parser state.
 */
typedef struct gen_prsr_t
{
    bool fa_sct;
    bool vcf_ftr_sct;
    size_t vcf_col;
    char* vcf_hdr;
    char** vcf_hdr_off;
} gen_prsr_t;
    
typedef struct gen_cbcks_t
{
    ldoc_doc_t* (*proc_ln)(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt, gen_fstat* stat);
} gen_cbcks_t;

typedef struct gen_ctxt_t
{
    gen_filetype_t tpe;
    bool py;
    char* fgen;
    char* pycall;
    char* fname;
    FILE* fout;
    char* usr;
} gen_ctxt_t;
    
void gen_init();
    
void gen_nde_dsc_opt(ldoc_nde_t* nde, ldoc_nde_t* dsc, char* lbl);
    
char* gen_res_opt(ldoc_res_t* res);
char* gen_res_optx(ldoc_res_t* res);
char* gen_res_req(ldoc_res_t* res);
    
char* gen_term_crnl(char* s);
void gen_lwr(char* str);
bi_attr gen_kwd(char* str);
char gen_inv(char c);
    
void gen_xcig(char* str);

char* gen_escstr(char* str);
    
size_t gen_csplit(char* str, char c);

bool gen_join_attrs_ent(char* id, ldoc_ent_t* ent, char* attrs);
bool gen_join_attrs_nde(char* id, ldoc_nde_t* nde, char* attrs);
void gen_splt_attrs(ldoc_nde_t* ftr, ldoc_nde_t* usr, ldoc_nde_t* ref, ldoc_nde_t* vars, char* attrs);

ldoc_nde_t* gen_variants(char* seq, char sep, char** vseqs, size_t* vnum);

void gen_rd(int fd, off_t mx, ldoc_trie_t* idx, gen_cbcks_t* cbcks, gen_ctxt_t* ctxt);

void gen_ser(gen_ctxt_t* ctxt, gen_ctpe_t ctpe, ldoc_doc_t* doc, ldoc_doc_t* opt, gen_fstat* stat);
    
char* qk_alloc(size_t n);
void qk_free();
char* qk_heap_ptr();
char* qk_working_ptr();
bool qk_heap_empty();
void qk_purge();
char* qk_strdup(const char* s1);
char* qk_strndup(const char* s1, size_t n);
bool qk_strcat(const char* s1);
    
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
