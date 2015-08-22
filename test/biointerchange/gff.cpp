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
#include "gff.h"

#define GFF_GFF_SMALL "../test-data/chromosome_BF.gff"
#define GFF_GFF_FA_OFF 72902
#define GFF_GFF_FA_REF1_ID 72903
#define GFF_GFF_FA_REF1_SEQ 72963
#define GFF_GFF_FA_REF1_LBRK 1
#define GFF_GFF_FA_REF1_LLEN 60
#define GFF_GFF_FA_NME1 "DDB0220052"

#define GFF_GFF_REGIONS "../test-data/regions.gff"
#define GFF_GFF_REGIONS_TMP "../test-data/regions.tmp"

/**
 * Pragma statements including the FASTA separator.
 */
#define GFF_GFF_FA_META 4

#define GFF_GFF_FA_FTRS 432

#define GFF_GFF_FA_COMMS 1

TEST(gff, fa_find)
{
    int fd = fio_opn(GFF_GFF_SMALL);
    size_t mx = fio_len(fd);
    
    gen_fstat stat = { 0, 0, 0, 0 };
    
    off_t off = gff_fnd_fa(fd, &stat, mx);
    
    // Was the right offset found?
    EXPECT_EQ(GFF_GFF_FA_OFF, off);
    
    // Are the number of seen meta data, comments and features right?
    EXPECT_EQ(GFF_GFF_FA_META, stat.meta);
    EXPECT_EQ(GFF_GFF_FA_FTRS, stat.ftrs);
    EXPECT_EQ(GFF_GFF_FA_COMMS, stat.comms);
    
    fio_cls(fd);
}

TEST(gff, fa_idx)
{
    int fd = fio_opn(GFF_GFF_SMALL);
    size_t mx = fio_len(fd);
    
    // Statistics have already been tested in `fa_find`.
    ldoc_trie_t* idx = gff_idx_fa(fd, NULL, mx);
    EXPECT_NE((ldoc_trie_t*)NULL, idx);

    fio_cls(fd);
}

TEST(gff, fa_idx_lookup)
{
    int fd = fio_opn(GFF_GFF_SMALL);
    size_t mx = fio_len(fd);
    
    // Statistics have already been tested in `fa_find`.
    ldoc_trie_t* idx = gff_idx_fa(fd, NULL, mx);
    EXPECT_NE((ldoc_trie_t*)NULL, idx);
    
    ldoc_trie_nde_t* res  = ldoc_trie_lookup(idx, GFF_GFF_FA_NME1, false);
    EXPECT_NE((ldoc_trie_nde_t*)NULL, res);
    EXPECT_EQ(res->anno.cat, BI_FASTA_REF);
    
    gen_fa_ntry_t* ntry = (gen_fa_ntry_t*)res->anno.pld;
    EXPECT_EQ(ntry->off, GFF_GFF_FA_REF1_ID);
    EXPECT_EQ(ntry->seq, GFF_GFF_FA_REF1_SEQ);
    EXPECT_EQ(ntry->lbrk, GFF_GFF_FA_REF1_LBRK);
    EXPECT_EQ(ntry->llen, GFF_GFF_FA_REF1_LLEN);
    
    fio_cls(fd);
}

TEST(gff, fa_idx_seq)
{
    int fd = fio_opn(GFF_GFF_SMALL);
    size_t mx = fio_len(fd);
    ldoc_trie_t* idx = gff_idx_fa(fd, NULL, mx);
    EXPECT_NE((ldoc_trie_t*)NULL, idx);
    
    char* seq = gff_seq(fd, mx, idx, GFF_GFF_FA_NME1, 140, 149, false);
    EXPECT_STREQ("TAGTCTACAA", seq);
    free(seq);

    seq = gff_seq(fd, mx, idx, GFF_GFF_FA_NME1, 140, 149, true);
    EXPECT_STREQ("TTGTAGACTA", seq);
    free(seq);
    
    fio_cls(fd);
}

TEST(gff, fa_idx_seq_lbrk)
{
    int fd = fio_opn(GFF_GFF_SMALL);
    size_t mx = fio_len(fd);
    ldoc_trie_t* idx = gff_idx_fa(fd, NULL, mx);
    EXPECT_NE((ldoc_trie_t*)NULL, idx);
    
    char* seq = gff_seq(fd, mx, idx, GFF_GFF_FA_NME1, 177, 182, false);
    EXPECT_STREQ("AGTATT", seq);
    free(seq);
    
    seq = gff_seq(fd, mx, idx, GFF_GFF_FA_NME1, 177, 182, true);
    EXPECT_STREQ("AATACT", seq);
    free(seq);
    
    fio_cls(fd);
}

TEST(gff, gff_rd)
{
    gen_init();
    
    int fd = fio_opn(GFF_GFF_SMALL);
    size_t mx = fio_len(fd);
    ldoc_trie_t* idx = gff_idx_fa(fd, NULL, mx);
    EXPECT_NE((ldoc_trie_t*)NULL, idx);
    
    gen_cbcks_t cbcks;
    gff_cbcks(&cbcks);
    
    gen_ctxt_t ctxt =
    {
        GEN_FMT_GFF3,
        false, // No Python API.
        false, // Not "quick".
        (char*)"unit-test.gff3", // Input filename; this file does not really exist.
        NULL,
        NULL,
        stdout,
        NULL, // User data.
        (char*)BIOINTERCHANGE_VERSION
    };
    
    gen_rd(fd, mx, idx, &cbcks, &ctxt);
    
    fio_cls(fd);
}

TEST(gff, gff_serialize_ftr)
{
    gen_init();
    
    ldoc_doc_t* doc = ldoc_doc_new();
    
    ldoc_nde_t* ftr = doc->rt;
    
    ldoc_ent_t* ctx = ldoc_ent_new(LDOC_ENT_OR);
    ctx->pld.pair.anno.str = (char*)JSONLD_CTX;
    ctx->pld.pair.dtm.str = (char*)JSONLD_GFF3_1;
    
    ldoc_nde_t* attrs = ldoc_nde_new(LDOC_NDE_UA);
    attrs->mkup.anno.str = (char*)GEN_ATTRS;
    
    ldoc_ent_t* lm = ldoc_ent_new(LDOC_ENT_OR);
    lm->pld.pair.anno.str = (char*)GFF_C1;
    lm->pld.pair.dtm.str = (char*)"Chr1";
    
    ldoc_nde_t* lc = ldoc_nde_new(LDOC_NDE_UA);
    lc->mkup.anno.str = (char*)GEN_LOCUS;
    
    ldoc_ent_t* src = ldoc_ent_new(LDOC_ENT_OR);
    src->pld.pair.anno.str = (char*)GFF_C2;
    src->pld.pair.dtm.str = NULL;
    
    ldoc_ent_t* tpe = ldoc_ent_new(LDOC_ENT_OR);
    tpe->pld.pair.anno.str = (char*)GFF_C3;
    tpe->pld.pair.dtm.str = NULL;
    
    ldoc_ent_t* st = ldoc_ent_new(LDOC_ENT_NR);
    st->pld.pair.anno.str = (char*)GFF_C4;
    st->pld.pair.dtm.str = (char*)"1000";
    
    ldoc_ent_t* en = ldoc_ent_new(LDOC_ENT_NR);
    en->pld.pair.anno.str = (char*)GFF_C5;
    en->pld.pair.dtm.str = (char*)"1001";
    
    ldoc_ent_t* scr = ldoc_ent_new(LDOC_ENT_OR);
    scr->pld.pair.anno.str = (char*)GFF_C6;
    scr->pld.pair.dtm.str = (char*)"12.3";
    
    ldoc_ent_t* strnd = ldoc_ent_new(LDOC_ENT_OR);
    strnd->pld.pair.anno.str = (char*)GFF_C7;
    strnd->pld.pair.dtm.str = (char*)"+";
    
    ldoc_ent_t* ph = ldoc_ent_new(LDOC_ENT_OR);
    ph->pld.pair.anno.str = (char*)GFF_C8;
    ph->pld.pair.dtm.str = (char*)"0";
    
    char* cmt = strdup("A \"comment\"!");
    ldoc_ent_t* c = ldoc_ent_new(LDOC_ENT_OR);
    c->pld.pair.anno.str = (char*)GEN_COMMENT;
    c->pld.pair.dtm.str = gen_escstr(cmt, GEN_FMT_GFF3);
    ldoc_nde_ent_push(ftr, c);
    
    ldoc_ent_t* sq = ldoc_ent_new(LDOC_ENT_OR);
    sq->pld.pair.anno.str = (char*)GEN_SEQUENCE;
    sq->pld.pair.dtm.str = (char*)"TATA";
    ldoc_nde_ent_push(ftr, sq);
    
    char* attrs_str = strdup("ID=DDB0220064;Parent=DDB_G0294346;Name=DDB0220064;description=BFNV1_0C0011_07644: Obtained from Geneid output run by Dictyostelium Genome Consortium at The Welcome Trust Sanger Institute;translation_start=1;Dbxref=Protein Accession Version:EAL60309.1,Inparanoid V. 5.1:DDB0220064,Protein Accession Number:EAL60309.1,Protein GI Number:60462053,UniProt:Q54AM3,Genome V. 2.0 ID:BFNV1_0C0011_07644;qualifier=Partial%2C 3%27 missing");
    gen_splt_attrs(ftr, attrs, NULL, NULL, attrs_str, BI_VAL);
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_nde_ent_push(ftr, ctx);
    
    // Source, type, score, strand, phase:
    ldoc_nde_ent_push(ftr, src);
    ldoc_nde_ent_push(ftr, tpe);
    ldoc_nde_ent_push(ftr, scr);
    ldoc_nde_ent_push(ftr, strnd);
    ldoc_nde_ent_push(ftr, ph);
    
    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, lm);
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    
    ldoc_nde_dsc_push(ftr, lc);
    
    ldoc_nde_dsc_push(ftr, attrs);
    
    char* qk = qk_alloc(10*1024*1024);
    
    ldoc_ser_t* ser = ldoc_format(doc, json_vis_nde, json_vis_ent);
    printf("%s\n", ser->pld.str);
    
    gff_proc_doc(doc, GEN_FMT_FTR);
    
    qk_free();
    
    free(cmt);
    free(attrs_str);
}

TEST(gff, gff_serialize_regions)
{
    gen_init();
    
    int fd = fio_opn(GFF_GFF_REGIONS);
    size_t mx = fio_len(fd);
    
    gen_cbcks_t cbcks;
    gff_cbcks(&cbcks);
    
    FILE* tmp = fopen(GFF_GFF_REGIONS_TMP, "w+");
    gen_ctxt_t ctxt =
    {
        GEN_FMT_GFF3,
        false, // No Python API.
        false, // Not "quick".
        (char*)"unit-test.gff", // Input filename; this file does not really exist.
        NULL,
        NULL,
        tmp,
        NULL, // User data.
        (char*)BIOINTERCHANGE_VERSION
    };
    
    gen_rd(fd, mx, NULL, &cbcks, &ctxt);
    
    fflush(tmp);
    fseek(tmp, 0, SEEK_SET);
    size_t lines = 0;
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line, &len, tmp)) != -1)
    {
        lines++;
        
        // Meta document:
        if (lines == 2)
        {
            // Make sure no "isolated" contig entries exist -- means that
            // the code did not pick up existing "contig" keys.
            EXPECT_EQ((char*)NULL, strstr(line, "\"contig\":{\"15\":{\"start\":1,\"end\":104043685}}"));
            
            // Ordering preserved by LibDocument:
            EXPECT_NE((char*)NULL, strstr(line, "\"contig\":{\"15\":{\"start\":1,\"end\":104043685},\"6\":{\"start\":1,\"end\":149736546},\"16\":{\"start\":1,\"end\":98207768},\"10\":{\"start\":1,\"end\":130694993},\"2\":{\"start\":1,\"end\":182113224},\"3\":{\"start\":1,\"end\":160039680},\"11\":{\"start\":1,\"end\":122082543},\"MT\":{\"start\":1,\"end\":16299},\"7\":{\"start\":1,\"end\":145441459},\"12\":{\"start\":1,\"end\":120129022},\"19\":{\"start\":1,\"end\":61431566},\"1\":{\"start\":1,\"end\":195471971},\"17\":{\"start\":1,\"end\":94987271},\"9\":{\"start\":1,\"end\":124595110},\"X\":{\"start\":1,\"end\":171031299},\"18\":{\"start\":1,\"end\":90702639},\"Y\":{\"start\":1,\"end\":91744698},\"4\":{\"start\":1,\"end\":156508116},\"14\":{\"start\":1,\"end\":124902244},\"13\":{\"start\":1,\"end\":120421639},\"5\":{\"start\":1,\"end\":151834684},\"8\":{\"start\":1,\"end\":129401213}}"));
        }
    }
    fclose(tmp);
    remove(GFF_GFF_REGIONS_TMP);
    
    // 1 context, 1 meta, 1 feature, 1 summary
    EXPECT_EQ(4, lines);
    
    fio_cls(fd);
}

