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

#ifndef biointerchange_gvf_h
#define biointerchange_gvf_h

#include "gen.h"
#include "gff.h"

#ifdef __cplusplus
extern "C" {
#endif
    
// Indices for structured pragma filters:
extern ldoc_trie_t* pgm_idx;

extern const char* GVF_C1;
extern const char* GVF_C2;
extern const char* GVF_C3;
extern const char* GVF_C4;
extern const char* GVF_C5;
extern const char* GVF_C6;
extern const char* GVF_C7;

typedef enum
{
    GVF_PRS_NONE = 0,
    GVF_PRS_REF,
    GVF_PRS_VAR
} gvf_prs_seq;
  
typedef enum
{
    GVF_PGM_TECHNOLOGY,
    GVF_PGM_DATA_SRC,
    GVF_PGM_SCR_MTD,
    GVF_PGM_SRC_MTD,
    GVF_PGM_ATTR_MTD,
    GVF_PGM_PHEN_DESC,
    GVF_PGM_PHSD_GT
} gvf_pgm_t;

void gvf_cbcks(gen_cbcks_t* cbcks);
    
//
// Internal use
//
    
void gvf_proc_effct(ldoc_nde_t* vars, char** val_cmp, bool* lend); // Called in gen.c
ldoc_doc_t* gvf_proc_ftr(int fd, off_t mx, ldoc_trie_t* idx, char* ln, size_t lnlen, char** cmt);
void gvf_tags(ldoc_nde_t* nde, char* strct, gvf_pgm_t tpe);
    
//
// Public API
//
    
void gvf_init();
    
ldoc_doc_t* gvf_proc_ln(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt, gen_fstat* stat);
    
char* gvf_proc_doc(ldoc_doc_t* doc, gen_doctype_t tpe);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
