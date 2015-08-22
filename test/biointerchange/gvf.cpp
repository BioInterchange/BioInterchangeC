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

#define GVF_GVF_SMALL "../test-data/small.gvf"

#define GVF_GVF_REGIONS "../test-data/regions.gvf"
#define GVF_GVF_REGIONS_TMP "../test-data/regions.tmp"

// Line taken from: test-data/Saccharomyces_cerevisiae_incl_consequences.gvf
// Modified to include more information about variants.
static const char* gvf_ftr1 = "I\tSGRP\tSNV\t675\t675\t.\t+\t.\tID=76;variant_peptide=1 P YAL068W-A,0 P YAL068W-A;Variant_seq=G,T;reference_peptide=P;Variant_effect=upstream_gene_variant 1 transcript YAL067W-A,upstream_gene_variant 0 transcript YAL067W-A,downstream_gene_variant 1 transcript YAL069W,downstream_gene_variant 0 transcript YAL069W,synonymous_variant 1 mRNA YAL068W-A,synonymous_variant 0 mRNA YAL068W-A,downstream_gene_variant 1 transcript YAL068C,downstream_gene_variant 0 transcript YAL068C;Dbxref=SGRP:s01-675;Variant_codon=GAG,GGG;Reference_seq=A;Reference_codon=GAG";

static const char* gvf_ftr2 = "I\tSGRP\tSNV\t675\t675\t.\t+\t.\tVariant_seq=.,-;Reference_seq=-";

// Line taken from: http://www.sequenceontology.org/resources/gvf.html
static const char* gvf_techpfm1 = "Seqid=chr13;Source=mpileup;Type=SNV,nucleotide_insertion,nucleotide_deletion;";

// Made up example:
static const char* gvf_techpfm2 = "Type=SNV,nucleotide_insertion,nucleotide_deletion;Comment=A simple comment;example_key=Example Value";

TEST(gvf, gvf_rd)
{
    gen_init();
    gvf_init();
    
    int fd = fio_opn(GVF_GVF_SMALL);
    size_t mx = fio_len(fd);
    ldoc_trie_t* idx = gff_idx_fa(fd, NULL, mx);
    EXPECT_EQ((ldoc_trie_t*)NULL, idx);
    
    gen_cbcks_t cbcks;
    gvf_cbcks(&cbcks);
    
    gen_ctxt_t ctxt =
    {
        GEN_FMT_GVF,
        false, // No Python API.
        false, // Not "quick".
        (char*)"unit-test.gvf", // Input filename; this file does not really exist.
        NULL,
        NULL,
        stdout,
        NULL, // User data.
        (char*)BIOINTERCHANGE_VERSION
    };
    
    gen_rd(fd, mx, idx, &cbcks, &ctxt);
    
    fio_cls(fd);
}

TEST(gvf, gvf_serialize_ftr)
{
    ldoc_doc_t* doc = ldoc_doc_new();
    qk_alloc(10*1024*1024);
    
    // Line taken from: test-data/Saccharomyces_cerevisiae_incl_consequences.gvf
    // Modified to include more information about variants.
    char* ln = strdup(gvf_ftr1);
    char* cmt = strdup("A \"comment\"!");
    doc = gvf_proc_ftr(0, 0, NULL, ln, strlen(ln), &cmt);

    qk_purge();
    
    gvf_proc_doc(doc, GEN_FMT_FTR);
    
    qk_free();
    
    free(cmt);
}

TEST(gvf, gvf_serialize_ftr_empty_seqs)
{
    ldoc_doc_t* doc = ldoc_doc_new();
    qk_alloc(10*1024*1024);
    
    // Line taken from: test-data/Saccharomyces_cerevisiae_incl_consequences.gvf
    // Modified to include more information about variants.
    char* ln = strdup(gvf_ftr2);
    char* cmt = NULL;
    doc = gvf_proc_ftr(0, 0, NULL, ln, strlen(ln), &cmt);
    
    qk_purge();
    
    gvf_proc_doc(doc, GEN_FMT_FTR);
    
    qk_free();
}

TEST(gvf, gvf_serialize_regions)
{
    gen_init();
    gvf_init();
    
    int fd = fio_opn(GVF_GVF_REGIONS);
    size_t mx = fio_len(fd);
    
    gen_cbcks_t cbcks;
    gvf_cbcks(&cbcks);
    
    FILE* tmp = fopen(GVF_GVF_REGIONS_TMP, "w+");
    gen_ctxt_t ctxt =
    {
        GEN_FMT_GVF,
        false, // No Python API.
        false, // Not "quick".
        (char*)"unit-test.gvf", // Input filename; this file does not really exist.
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
    remove(GVF_GVF_REGIONS_TMP);
    
    // 1 context, 1 meta, 1 feature, 1 summary
    EXPECT_EQ(4, lines);
    
    fio_cls(fd);
}

TEST(gvf, gvf_tags_technology)
{
    gen_init();
    gvf_init();
    
    ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_UA);
    EXPECT_NE((ldoc_nde_t*)NULL, nde);

    char* pgm = strdup(gvf_techpfm1);
    EXPECT_NE((char*)NULL, pgm);
    
    gvf_tags(nde, pgm, GVF_PGM_TECHNOLOGY);

    free(pgm);
    
    pgm = strdup(gvf_techpfm2);
    EXPECT_NE((char*)NULL, pgm);
    
    gvf_tags(nde, pgm, GVF_PGM_TECHNOLOGY);

    // TODO Add more value checks!
    
    free(pgm);
}

TEST(gvf, gvf_tags_index_use)
{
    gen_init();
    gvf_init();
    
    // Index:
    
    ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_UA);
    EXPECT_NE((ldoc_nde_t*)NULL, nde);
    
    char* pgm = strdup(gvf_techpfm1);
    EXPECT_NE((char*)NULL, pgm);
    
    gvf_tags(nde, pgm, GVF_PGM_TECHNOLOGY);
    
    free(pgm);
    
    pgm = strdup(gvf_techpfm2);
    EXPECT_NE((char*)NULL, pgm);
    
    gvf_tags(nde, pgm, GVF_PGM_TECHNOLOGY);
    
    // Features:
    
    ldoc_doc_t* doc = ldoc_doc_new();
    qk_alloc(10*1024*1024);
    
    // Line taken from: test-data/Saccharomyces_cerevisiae_incl_consequences.gvf
    // Modified to include more information about variants.
    char* ln = strdup(gvf_ftr1);
    char* cmt = NULL;
    doc = gvf_proc_ftr(0, 0, NULL, ln, strlen(ln), &cmt);
    
    qk_purge();

    ldoc_ser_t* ser = ldoc_format(doc, json_vis_nde, json_vis_ent);
    printf("%s\n", ser->pld.str);

    qk_purge();
    
    gvf_proc_doc(doc, GEN_FMT_FTR);
    
    qk_free();

    free(pgm);
}

