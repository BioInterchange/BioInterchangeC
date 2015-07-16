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
#include <json.h>

#include "ext-python.h"
#include "fio.h"

#ifdef __cplusplus
extern "C" {
#endif
    
// Error codes (exit codes):
#define MAIN_SUCCESS  0
    
// Never change these definitions or it will be very hard to interpret
// error codes across software versions!
#define MAIN_ERR_PGSZ 1
#define MAIN_ERR_HME1 2
#define MAIN_ERR_HME2 3
#define MAIN_ERR_PARA 4
// License file size wrong:
#define MAIN_ERR_LISZ 5
// License file format wrong:
#define MAIN_ERR_LICF 6
// License check network error:
#define MAIN_ERR_LICN 7
// License was rejected by the license server:
#define MAIN_ERR_LICV 8
#define MAIN_ERR_FNME 9
#define MAIN_ERR_FEXT 10
#define MAIN_ERR_FACC 11
// PEM_read_bio_X509 failed:
#define MAIN_ERR_SBIO 12
// X509_STORE_add_cert failed:
#define MAIN_ERR_SADD 13
    
// Ensure that no statistics are sent during license check:
#define GEN_STATS_PRIVATE 1
    
// Number of memory pages that should reside in memory when reading a genomics file:
#define BI_GEN_PG_MUL 4096

extern const char* JSONLD_CTX;

extern const char* JSONLD_GFF3_CTX1;
extern const char* JSONLD_GTF_CTX1;
extern const char* JSONLD_GVF_CTX1;
extern const char* JSONLD_VCF_CTX1;

extern const char* JSONLD_GFF3_X1;
extern const char* JSONLD_GTF_X1;
extern const char* JSONLD_GVF_X1;
extern const char* JSONLD_VCF_X1;

extern const char* JSONLD_GFF3_1;
extern const char* JSONLD_GTF_1;
extern const char* JSONLD_GVF_1;
extern const char* JSONLD_VCF_1;
    
extern const char* JSONLD_STAT_1;
    
extern const char* GEN_AFFECTED;
extern const char* GEN_AFFECTED_TPE;
extern const char* GEN_ALLELE_CNT;
extern const char* GEN_ALLELE_CNT_VCF;
extern const char* GEN_ALLELE_FRQ;
extern const char* GEN_ALLELE_FRQ_VCF;
extern const char* GEN_ALLELE_TTL;
extern const char* GEN_ALLELE_TTL_VCF;
extern const char* GEN_ALIGNMENT;
extern const char* GEN_ALIGNMENT_GFF3;
extern const char* GEN_ANNOTATIONS;
extern const char* GEN_ATTRIBUTE_MTHD;
extern const char* GEN_ATTRS;
extern const char* GEN_AVG_COVERAGE;
extern const char* GEN_AVG_COVERAGE_GVF;
extern const char* GEN_BUILD;
extern const char* GEN_BUILD_VAL;
extern const char* GEN_CIGAR;
extern const char* GEN_CIGAR_GFF3;
extern const char* GEN_CODON;
extern const char* GEN_COMMENT;
extern const char* GEN_DATA_SRC;
extern const char* GEN_DBXREF;
extern const char* GEN_DBXREF_GFF3;
extern const char* GEN_DEPTH;
extern const char* GEN_DEPTH_VCF;
extern const char* GEN_EFFECT;
extern const char* GEN_EFFECTS;
extern const char* GEN_END;
extern const char* GEN_GENOMIC_SRC;
extern const char* GEN_GLOBAL;
extern const char* GEN_ID;
extern const char* GEN_ID_GFF3;
extern const char* GEN_LANDMARK;
extern const char* GEN_LANDMARKS;
extern const char* GEN_LOCUS;
extern const char* GEN_ONT_ACCESSION;
extern const char* GEN_ONT_TERM;
extern const char* GEN_PHASED_GENO;
extern const char* GEN_PHENO_DESCR;
extern const char* GEN_POP;
extern const char* GEN_QUALITY_MAP;
extern const char* GEN_QUALITY_MAP0;
extern const char* GEN_QUALITY_RMS;
extern const char* GEN_READ_PAIR_SPAN;
extern const char* GEN_READ_PAIR_SPAN_GVF;
extern const char* GEN_REFERENCE;
extern const char* GEN_SAMPLES_DATA;
extern const char* GEN_SCORE_MTHD;
extern const char* GEN_SEQUENCE;
extern const char* GEN_SEQUENCES;
extern const char* GEN_START;
extern const char* GEN_STRAND;
extern const char* GEN_SOURCE;
extern const char* GEN_SOURCES;
extern const char* GEN_SOURCE_MTHD;
extern const char* GEN_TECHNOLOGY;
extern const char* GEN_TYPE;
extern const char* GEN_TYPE_GVF;
extern const char* GEN_TYPES;
extern const char* GEN_VARIANTS;
extern const char* GEN_VCFVERSION;
extern const char* GEN_VCFVERSION_VCF;

extern const char* GEN_EMPTY;
extern const char* GEN_NULL;
extern const char* GEN_TRUE;
extern const char* GEN_UNKNOWN;
extern const char* GEN_FORWARD;
    
// Statistics wrap-up:
extern const char* GEN_STATS;
extern const char* GEN_STAT_CMMS;
extern const char* GEN_STAT_META;
extern const char* GEN_STAT_META_FLTR;
extern const char* GEN_STAT_FTRS;
extern const char* GEN_STAT_FTRS_FLTR;
extern const char* GEN_STAT_RNTM;
extern const char* GEN_STAT_INVC;
extern const char* GEN_STAT_FNSH;
extern const char* GEN_STAT_LSEC;
    
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
     * Comma separated values of pairs (colon separated):
     */
    BI_CSEPCPAIR,
    /**
     * Comma separated information about variants.
     */
    BI_CSEPVAR,
    /**
     * Comma separated information about variants that is a number.
     */
    BI_CSEPVARN,
    /**
     * Comma separated information about variants & skip 8 characters of the key.
     */
    BI_CSEPVAR8,
    /**
     * GVF "Variant_effect" key.
     */
    BI_GVFEFFECT,
    /**
     * Information about the reference sequence.
     */
    BI_REFSEQ,
    /**
     * Information about the reference sequence that is a number.
     */
    BI_REFSEQN,
    /**
     * Information about the reference sequence & skip 10 characters of the key.
     */
    BI_REFSEQ10,
    /**
     * CIGAR string in GFF3 format (needs reformatting to match actual CIGAR specification).
     */
    BI_XCIG,
    /**
     * Target string in GFF3 format.
     */
    BI_TARGET,
    BI_STRUCT
} bi_attr;
    
typedef struct gen_attr_t
{
    bi_attr attr;
    const char* alt;
} gen_attr_t;

typedef enum
{
    GEN_FMT_GFF3,
    GEN_FMT_GTF,
    GEN_FMT_GVF,
    GEN_FMT_VCF,
    GEN_FMT_LDJ
} gen_filetype_t;

typedef enum
{
    GEN_FMT_CTX,
    GEN_FMT_INF,
    GEN_FMT_FTR,
    GEN_FMT_STT
} gen_doctype_t;
    
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
    char* ver;
} gen_ctxt_t;
    
void gen_init();
    
void gen_nde_dsc_opt(ldoc_nde_t* nde, ldoc_nde_t* dsc, char* lbl);
    
char* gen_res_opt(ldoc_res_t* res);
char* gen_res_optx(ldoc_res_t* res);
char* gen_res_req(ldoc_res_t* res);
    
char* gen_term_crnl(char* s);
void gen_ky(char* attr, char** val);
void gen_lwr(char* str);
void gen_lwrhyph(char* str);
char* gen_quoskp(char* str);
void gen_kwd(char* str, gen_attr_t* attr, bi_attr upfail);
char gen_inv(char c);
    
/**
 *  Returns the ordered-list node that was added to `dst`.
 */
ldoc_nde_t* gen_csep(ldoc_nde_t* dst, gen_attr_t kwd, char* ky, char* val);
ldoc_nde_t* gen_csep_dup(ldoc_nde_t* dst, gen_attr_t kwd, char* ky, char* val, bool dup);

void gen_xcig(char* str);
void gen_qk_revcig(char* str);
    
char* gen_escstr(char* str, gen_filetype_t tpe);
    
size_t gen_csplit(char* str, char c);

void gff_proc_tgt(ldoc_nde_t* nde, char* val);
    
ldoc_nde_t* gen_ctx(ldoc_nde_t* nde, bool* nw);
ldoc_nde_t* gen_find_nde(ldoc_nde_t* ctnr1, ldoc_nde_t* ctnr2, char* ky);
void gen_add_nw(ldoc_nde_t* nde, ldoc_nde_t* usr);

ldoc_content_t gen_smrt_tpe(char* val);
ldoc_content_t gen_smrt_flttpe(char* val);
    
bool gen_join_nde(ldoc_nde_t* nde);
bool gen_join_attrs_key(char* id, ldoc_nde_t* nde, ldoc_ent_t* ent, char* attrs);
bool gen_join_attrs_ent(char* id, ldoc_ent_t* ent, char* attrs);
bool gen_join_attrs_nde(char* id, ldoc_nde_t* nde, char* attrs);
void gen_splt_attrs(ldoc_nde_t* ftr, ldoc_nde_t* usr, ldoc_nde_t* ref, ldoc_nde_t* vars, char* attrs, bi_attr upfail);

ldoc_nde_t* gen_variants(char* seq, char sep, char** vseqs, size_t* vnum);

void gen_rd(int fd, off_t mx, ldoc_trie_t* idx, gen_cbcks_t* cbcks, gen_ctxt_t* ctxt);

void gen_rd_doc(int fd, off_t mx, gen_ctxt_t* ctxt);
    
void gen_ser(gen_ctxt_t* ctxt, gen_ctpe_t ctpe, ldoc_doc_t* doc, ldoc_doc_t* opt, gen_fstat* stat);
    
//
// Genomics file formatting
//

bool gen_proc_doc_usr(ldoc_nde_t* ftr);
bool gen_proc_nde(ldoc_nde_t* vars, char* attr, char* pre, char* astr, size_t vnum);
    
char* qk_alloc(size_t n);
void qk_free();
char* qk_heap_ptr();
char* qk_working_ptr();
bool qk_heap_empty();
void qk_purge();
char* qk_strdup(const char* s1);
char* qk_strndup(const char* s1, size_t n);
bool qk_strcat(const char* s1);
bool qk_strncat(const char* s1, size_t n);
    
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
