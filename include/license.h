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

#undef BIOINTERCHANGE_NOCRYPT

#ifndef biointerchange_license_h
#define biointerchange_license_h

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef BIOINTERCHANGE_CRYPT
#include <openssl/conf.h>
#include <openssl/err.h>
#include <openssl/evp.h>
#include <openssl/rand.h>
#endif // BIOINTERCHANGE_CRYPT

#define CURL_STATICLIB 1
#include <curl.h>
#include <openssl/ssl.h>

#include "gen.h"

#ifdef NDEBUG
#define LIC_URL     "https://www.codamono.com/license/"
#else
//#define LIC_URL     "http://localhost:8000/license/"
#define LIC_URL     "https://www.codamono.com/license/"
#endif // NDEBUG

#ifdef __APPLE__
#define EXE_SYMID   "mH2kuYvd0KYBpThg"
#else
#define EXE_SYMID   "9SDYi2ZXv3JzfWJn"
#endif

#ifdef BIOINTERCHANGE_CRYPT
#define LIC_KEYLEN  32
#define LIC_IVLEN   16

#define EXE_MEMKEY  "\x6e\x84\x8a\xd7\x0c\xcb\x58\x95\x5f\x67\x0d\x9b\x16\xf5\x56\x6e\x58\x1b\x37\x5e\xa2\x02\x0c\x98\x07\x5d\x05\x63\xff\x29\xe9\x4e"
#define EXE_MEMIV   "\x6e\x3f\x2d\xee\xee\x0c\xc2\x10\xbd\xa0\x7a\x53\xa2\x96\xe6\x7b"

// Both keys encrypted with mem-key and -iv:
#define EXE_SYMKEY1 "\x0c\x32\xda\x4a\x82\x04\xf3\x54\x8a\x73\x2c\x89\xa2\x85\xfc\x08\xea\x79\xa8\x26\x69\x95\xe0\x91\xcc\x5d\x65\x69\x06\x4e\x8f\x1c\xcd\xfd\x63\xa9\x18\xfb\x96\xa8\xa3\xe3\x28\xc3\x77\xa7\x75\xae"
#define EXE_SYMKEY2 "\xba\xb4\xbd\xe7\x55\x2a\x9c\x34\xde\x1b\x09\x70\xe0\xe6\x3c\x3d\x6c\x6a\x97\xcc\x8a\xcc\xc0\x73\x15\x09\x73\xcb\xae\xd7\x5a\x5e\x20\xf2\xbf\x5f\xff\x7f\x00\x00\xc3\x6c\x44\x62\xe6\x0b\x00\x37"

#define EXE_SYMLEN1 48
#define EXE_SYMLEN2 48
#endif // BIOINTERCHANGE_CRYPT

#ifdef __cplusplus
extern "C" {
#endif

/**
 *
 * Note: Do not change the order of these entries; only add to the end.
 *       Changing the order will break communication with the Python
 *       server implementation.
 */
typedef enum
{
    /**
     * License valid.
     */
    LICENSE_OK = 0,
    /**
     * Network/communiation error occurred.
     */
    LICENSE_NET,
    /**
     * Not an encoding token.
     */
    LICENSE_NENC,
    /**
     * License key not recorded.
     */
    LICENSE_NREC,
    /**
     * Invalid format.
     */
    LICENSE_INVFMT,
    /**
     * Internal data error.
     */
    LICENSE_INT,
    /**
     * License expired.
     */
    LICENSE_EXP,
    /**
     * License limit reached.
     */
    LICENSE_LMT,
    /**
     * Server response garbled.
     */
    LICENSE_SRV
} lic_status_t;

typedef enum
{
    ENCODING_OK,
    ENCODING_AES
} lic_enc_t;

typedef enum
{
    DECODING_OK,
    DECODING_AES
} lic_dec_t;

typedef struct lic_id
{
    lic_status_t status;
    char* lid;
} lic_id;
    
typedef struct lic_chksum_t
{
    uint8_t sum;
    uint16_t alt;
} lic_chksum_t;

lic_status_t lic_valid(char* lstr, gen_fstat* stat);
    
#ifdef BIOINTERCHANGE_CRYPT
lic_enc_t lic_enc_sym(unsigned char* key, unsigned char* iv, const char* plaintext, size_t len, uint8_t** ciphertext, size_t* cipherlen);
lic_dec_t lic_dec_sym(unsigned char* key, unsigned char* iv, const char* cipher, size_t len, char** plaintext, size_t* plainlen);
#endif // BIOINTERCHANGE_CRYPT

char* lic_raw2escstr(const unsigned char* raw, size_t len);

    
    
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
