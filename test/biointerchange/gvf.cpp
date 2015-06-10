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

#include <gtest/gtest.h>

#include <unistd.h>

#include "fio.h"
#include "gvf.h"

TEST(gvf, gvf_serialize_ftr)
{
    ldoc_doc_t* doc = ldoc_doc_new();
    
    // Line taken from: test-data/Saccharomyces_cerevisiae_incl_consequences.gvf
    char* ln = strdup("I\tSGRP\tSNV\t675\t675\t.\t+\t.\tID=76;variant_peptide=1 P YAL068W-A,0 P YAL068W-A;Variant_seq=G,T;reference_peptide=P;Variant_effect=upstream_gene_variant 1 transcript YAL067W-A,upstream_gene_variant 0 transcript YAL067W-A,downstream_gene_variant 1 transcript YAL069W,downstream_gene_variant 0 transcript YAL069W,synonymous_variant 1 mRNA YAL068W-A,synonymous_variant 0 mRNA YAL068W-A,downstream_gene_variant 1 transcript YAL068C,downstream_gene_variant 0 transcript YAL068C;Dbxref=SGRP:s01-675;Reference_seq=A");
    char* cmt = strdup("A \"comment\"!");
    doc = gvf_proc_ftr(0, 0, NULL, ln, strlen(ln), &cmt);
    
    char* qk = qk_alloc(10*1024*1024);
    
    gvf_proc_doc(doc);
    
    qk_free();
    
    free(cmt);
}

