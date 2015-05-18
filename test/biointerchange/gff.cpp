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
    
    gff_rd(fd, mx, idx);
    
    fio_cls(fd);
}
