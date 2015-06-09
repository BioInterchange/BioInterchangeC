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
    
    ldoc_nde_t* ftr = doc->rt;
    
    ldoc_ent_t* ctx = ldoc_ent_new(LDOC_ENT_OR);
    ctx->pld.pair.anno.str = (char*)JSONLD_CTX;
    ctx->pld.pair.dtm.str = (char*)JSONLD_GFF3;
    
    ldoc_nde_t* attrs = ldoc_nde_new(LDOC_NDE_UA);
    attrs->mkup.anno.str = (char*)GEN_ATTRS;
    
    ldoc_ent_t* lm = ldoc_ent_new(LDOC_ENT_OR);
    lm->pld.pair.anno.str = (char*)GVF_C1;
    lm->pld.pair.dtm.str = (char*)"Chr1";
    
    ldoc_nde_t* lc = ldoc_nde_new(LDOC_NDE_UA);
    lc->mkup.anno.str = (char*)GEN_LOCUS;
    
    ldoc_ent_t* src = ldoc_ent_new(LDOC_ENT_OR);
    src->pld.pair.anno.str = (char*)GVF_C2;
    src->pld.pair.dtm.str = NULL;
    
    ldoc_ent_t* tpe = ldoc_ent_new(LDOC_ENT_OR);
    tpe->pld.pair.anno.str = (char*)GVF_C3;
    tpe->pld.pair.dtm.str = NULL;
    
    ldoc_ent_t* st = ldoc_ent_new(LDOC_ENT_NR);
    st->pld.pair.anno.str = (char*)GVF_C4;
    st->pld.pair.dtm.str = (char*)"1000";
    
    ldoc_ent_t* en = ldoc_ent_new(LDOC_ENT_NR);
    en->pld.pair.anno.str = (char*)GVF_C5;
    en->pld.pair.dtm.str = (char*)"1001";
    
    ldoc_ent_t* scr = ldoc_ent_new(LDOC_ENT_OR);
    scr->pld.pair.anno.str = (char*)GVF_C6;
    scr->pld.pair.dtm.str = (char*)"12.3";
    
    ldoc_ent_t* strnd = ldoc_ent_new(LDOC_ENT_OR);
    strnd->pld.pair.anno.str = (char*)GVF_C7;
    strnd->pld.pair.dtm.str = (char*)"+";
    
    char* cmt = strdup("A \"comment\"!");
    ldoc_ent_t* c = ldoc_ent_new(LDOC_ENT_OR);
    c->pld.pair.anno.str = (char*)GEN_COMMENT;
    c->pld.pair.dtm.str = gen_escstr(cmt);
    ldoc_nde_ent_push(ftr, c);
    
    ldoc_ent_t* sq = ldoc_ent_new(LDOC_ENT_OR);
    sq->pld.pair.anno.str = (char*)GEN_SEQUENCE;
    sq->pld.pair.dtm.str = (char*)"TATA";
    ldoc_nde_ent_push(ftr, sq);
    
    char* attrs_str = strdup("ID=1;Variant_seq=A;Variant_effect=upstream_gene_variant 0 transcript YAL067W-A,upstream_gene_variant 0 transcript YAL069W,upstream_gene_variant 0 transcript YAL068W-A,downstream_gene_variant 0 transcript YAL068C;Dbxref=SGRP:s01-84;Reference_seq=G");
    gen_splt_attrs(ftr, attrs, NULL, NULL, attrs_str);
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_nde_ent_push(ftr, ctx);
    
    // Source, type, score, strand:
    ldoc_nde_ent_push(ftr, src);
    ldoc_nde_ent_push(ftr, tpe);
    ldoc_nde_ent_push(ftr, scr);
    ldoc_nde_ent_push(ftr, strnd);
    
    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, lm);
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    
    ldoc_nde_dsc_push(ftr, lc);
    
    ldoc_nde_dsc_push(ftr, attrs);
    
    char* qk = qk_alloc(10*1024*1024);
    
    gvf_proc_doc(doc);
    
    qk_free();
    
    free(cmt);
    free(attrs_str);
}

