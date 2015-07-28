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
#include "vcf.h"

// Lines taken from: test-data/mgp.v3.indels.rsIDdbSNPv137.vcf
// Command: head -n 69 test-data/mgp.v3.indels.rsIDdbSNPv137.vcf | tail -n 2 | tr $'\t' '^' | sed 's/\^/\\t/g'
static const char* vcf_ftr1_1 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t129P2\t129S1\t129S5\tAJ\tAKRJ\tBALBcJ\tC3HHeJ\tC57BL6NJ\tCASTEiJ\tCBAJ\tDBA2J\tFVBNJ\tLPJ\tNODShiLtJ\tNZOHlLtJ\tPWKPhJ\tSPRETEiJ\tWSBEiJ";
static const char* vcf_ftr1_2 = "1\t3000019\t.\tG\tGA\t40.49\tQual;MinAB;MinDP\tAC1=1;AC=36;AF1=1;AN=36;DP4=0,0,71,0;DP=75;INDEL;MQ=29;VDB=0.0006\tGT:GQ:DP:SP:PL:FI\t1/1:68:2:0:53,6,0:0\t1/1:99:6:0:103,18,0:1\t1/1:99:3:0:72,9,0:1\t1/1:99:3:0:69,9,0:0\t1/1:99:3:0:69,9,0:0\t1/1:78:3:0:58,9,0:0\t1/1:99:6:0:108,18,0:1\t1/1:99:5:0:98,15,0:1\t1/1:90:3:0:64,9,0:0\t1/1:99:6:0:104,18,0:1\t1/1:90:3:0:64,9,0:0\t1/1:34:2:0:36,6,0:0\t1/1:99:4:0:81,12,0:0\t1/1:99:4:0:75,12,0:0\t1/1:99:4:0:75,12,0:0\t1/1:99:4:0:85,12,0:0\t1/1:99:5:0:91,15,0:1\t1/1:99:5:0:95,15,0:1";

TEST(vcf, vcf_serialize_ftr)
{
    gen_init();
    char* qk = qk_alloc(10*1024*1024);
    
    ldoc_doc_t* doc;
    ldoc_doc_t* fdoc = ldoc_doc_new();
    
    gen_prsr_t st;
    gen_fstat stat;
    char* ln = strdup(vcf_ftr1_1);
    doc = vcf_proc_ln(0, strlen(vcf_ftr1_1), fdoc, NULL, ln, strlen(vcf_ftr1_1), &st, NULL, &stat);
    
    ln = strdup(vcf_ftr1_2);
    char* cmt = strdup("A \"comment\"!");
    doc = vcf_proc_ln(0, strlen(vcf_ftr1_2) + 1, fdoc, NULL, ln, strlen(vcf_ftr1_2) + 1, &st, &cmt, &stat);
    
    ldoc_ser_t* ser = ldoc_format(doc, json_vis_nde, json_vis_ent);
    printf("%s\n", ser->pld.str);
    
    qk_purge();
    
    // TODO: This fails because the document replies on the quickstack, which
    //       is purged here.
    // vcf_proc_doc(doc, GEN_FMT_FTR);
    
    qk_free();
    
    free(cmt);
}

