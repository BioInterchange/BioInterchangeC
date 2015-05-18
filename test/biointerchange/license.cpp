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

#include "license.h"

#define LIC_EX_A "A8e3jba9r7C21"
#define LIC_EX_B "B8F3jba9r7C21"

TEST(lic, lic_local)
{
    lic_valid(LIC_EX_A, NULL);
}

TEST(lic, lic_http)
{
    lic_valid(LIC_EX_B, NULL);
}

TEST(lic, lic_genkey)
{
    unsigned char key1[LIC_KEYLEN];
    unsigned char key2[LIC_KEYLEN];
    
    int ret = RAND_bytes(key1, LIC_KEYLEN);
    EXPECT_EQ(1, ret);
    
    ret = RAND_bytes(key2, LIC_KEYLEN);
    EXPECT_EQ(1, ret);
    
    EXPECT_NE(0, memcmp(key1, key2, LIC_KEYLEN));
}

TEST(lic, lic_raw2escstr)
{
    unsigned char key[LIC_KEYLEN];
    
    int ret = RAND_bytes(key, LIC_KEYLEN);
    EXPECT_EQ(1, ret);
    
    char* str = lic_raw2escstr(key, LIC_KEYLEN);
    EXPECT_NE((char*)NULL, str);
    
    // printf("%s\n", str);
    
    free(str);
    
    str = lic_raw2escstr(key, LIC_KEYLEN / 2);
    EXPECT_NE((char*)NULL, str);
    
    // printf("%s\n", str);
    
    free(str);
}

TEST(lic, lic_genkey_enc)
{
    unsigned char key[LIC_KEYLEN];
    int ret = RAND_bytes(key, LIC_KEYLEN);
    EXPECT_EQ(1, ret);
    
    char* str;
    
    // For verifying key generation:
    str = lic_raw2escstr(key, LIC_KEYLEN);
    //printf("%s\n", str);
    free(str);
    
    uint8_t* cipher;
    size_t clen;
    lic_enc_t enc = lic_enc_sym((unsigned char*)EXE_MEMKEY, (unsigned char*)EXE_MEMIV, (char*)key, LIC_KEYLEN, &cipher, &clen);
    EXPECT_EQ(ENCODING_OK, enc);
    
    str = lic_raw2escstr(cipher, clen);
    EXPECT_NE((char*)NULL, str);

    // The actual encrypted key:
    //printf("%s\n", str);
    
    free(str);
}

TEST(lic, lic_genkey_dec)
{
    char* plainkey;
    size_t len;
    lic_dec_t dec = lic_dec_sym((unsigned char*)EXE_MEMKEY, (unsigned char*)EXE_MEMIV, EXE_SYMKEY1, EXE_SYMLEN1, &plainkey, &len);
    EXPECT_EQ(DECODING_OK, dec);
    
    char* str = lic_raw2escstr((unsigned char*)plainkey, LIC_KEYLEN);
    EXPECT_NE((char*)NULL, str);
    
    //printf("%s\n", str);
    
    free(str);
}

TEST(lic, lic_memkey_license)
{
    char* str = (char*)LIC_EX_B;
    
    uint8_t* cipher;
    size_t clen;
    lic_enc_t enc = lic_enc_sym((unsigned char*)EXE_MEMKEY, (unsigned char*)EXE_MEMIV, (char*)str, strlen(str), &cipher, &clen);
    EXPECT_EQ(ENCODING_OK, enc);
    
    str = lic_raw2escstr(cipher, clen);
    EXPECT_NE((char*)NULL, str);
    
    // The actual encrypted key:
    // printf("%s\n", str);
    
    free(str);
    free(cipher);

}