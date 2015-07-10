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

#ifndef __biointerchange__vcf__
#define __biointerchange__vcf__

#include <inttypes.h>
#include <stdio.h>

#include "document.h"

#include "gen.h"
#include "gff.h"

#ifdef __cplusplus
extern "C" {
#endif
    
typedef enum
{
    BI_VCF_OTHER,
    BI_VCF_GT
} bi_vcf_info;
    
void vcf_cbcks(gen_cbcks_t* cbcks);

char* vcf_proc_doc_ftr(ldoc_nde_t* ftr);
char* vcf_proc_doc_ftr_attrs(ldoc_nde_t* ftr);

ldoc_doc_t* vcf_proc_ln(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt, gen_fstat* stat);

char* vcf_proc_doc(ldoc_doc_t* doc, gen_doctype_t tpe);
    
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* defined(__biointerchange__vcf__) */
