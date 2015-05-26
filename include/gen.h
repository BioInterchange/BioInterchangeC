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
#include <unistd.h>

#include <document.h>

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
extern const char* GEN_SEQUENCE;
extern const char* GEN_START;
extern const char* GEN_SOURCE;

extern const char* GEN_NULL;
extern const char* GEN_TRUE;

extern ldoc_vis_nde_ord_t* json_vis_nde;;
extern ldoc_vis_ent_t* json_vis_ent;
    
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
     * Comma separated information about variants.
     */
    BI_CSEPVAR,
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
    void (*proc_ln)(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt);
} gen_cbcks_t;

void gen_init();
    
void gen_lwr(char* str);
bi_attr gen_kwd(char* str);
char gen_inv(char c);
    
void gen_xcig(char* str);

char* gen_escstr(char* str);
    
size_t gen_csplit(char* str, char c);

void gen_splt_attrs(ldoc_nde_t* ftr, ldoc_nde_t* usr, ldoc_nde_t* vars, char* attrs);

ldoc_nde_t* gen_variants(char* seq, char sep);

void gen_rd(int fd, off_t mx, ldoc_trie_t* idx, gen_cbcks_t* cbcks);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
