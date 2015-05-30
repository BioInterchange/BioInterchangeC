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

typedef enum
{
    GVF_PRS_NONE = 0,
    GVF_PRS_REF,
    GVF_PRS_VAR
} gvf_prs_seq;
    
void gvf_proc_ln(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt);
    
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
