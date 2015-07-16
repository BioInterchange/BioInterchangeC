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

#include "gen.h"

#define GFF3_CIGAR_1 "M8 D41 I1 M321"
#define GFF3_CIGAR_2 "M8D41I1M321"
#define GEN_CIGAR "8M41D1I321M"

TEST(gen, cigar_spaces)
{
    char* cig = strdup(GFF3_CIGAR_1);
    EXPECT_NE((char*)NULL, cig);
    
    gen_xcig(cig);
    EXPECT_STREQ(GEN_CIGAR, cig);
    
    free(cig);
}

TEST(gen, cigar_continuous)
{
    char* cig = strdup(GFF3_CIGAR_2);
    EXPECT_NE((char*)NULL, cig);
    
    gen_xcig(cig);
    EXPECT_STREQ(GEN_CIGAR, cig);
    
    free(cig);
}

TEST(gen, cigar_reverse)
{
    qk_alloc(1024);
    
    gen_qk_revcig((char*)GEN_CIGAR);
    EXPECT_STREQ(GFF3_CIGAR_1, qk_heap_ptr());
    
    qk_free();
}
