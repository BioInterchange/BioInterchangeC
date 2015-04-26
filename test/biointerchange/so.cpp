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

#include "so.h"

#define ONT_GAP_URI "http://purl.obolibrary.org/obo/SO_0000730"

TEST(so, trie_new)
{
    ldoc_trie_t* trie = ont_trie_new();
    
    EXPECT_NE((ldoc_trie_t*)NULL, trie);
    
    ldoc_trie_free(trie);
}

TEST(so, lookup)
{
    ldoc_trie_t* trie = ont_trie_new();
    
    EXPECT_NE((ldoc_trie_t*)NULL, trie);
    
    ldoc_trie_nde_t* nde = ldoc_trie_lookup(trie, "gap", false);

    EXPECT_NE((ldoc_trie_nde_t*)NULL, nde);
    EXPECT_STREQ(ONT_GAP_URI, (char*)nde->anno.pld);
    
    nde = ldoc_trie_lookup(trie, "gop123", false);
    
    EXPECT_EQ((ldoc_trie_nde_t*)NULL, nde);
    
    ldoc_trie_free(trie);
}
