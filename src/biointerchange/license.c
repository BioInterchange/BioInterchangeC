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

#include "license.h"

char* lic_raw2escstr(const unsigned char* raw, size_t len)
{
    char* str = (char*)malloc(len * 4 + 1);
    
    if (!str)
    {
        // TODO Error handling.
        
        return NULL;
    }

    char* ptr = str;
    while (len--)
    {
        snprintf(ptr, 5, "\\x%02x", *(raw++));
        
        ptr += 4;
    }
    *ptr = 0;
    
    return str;
}

static inline char* lic_cipher_raw2hex(uint8_t* cipher, size_t len)
{
    // TODO Does fail for len == 0.
    
    char* str = (char*)malloc(len * 2 + 1);
    
    if (!str)
    {
        // TODO Error handling.
    }
    
    char* dst = str;
    unsigned char* src = cipher;
    while (len--)
    {
        snprintf(dst, 3, "%02x", *src);
        
        src++;
        dst += 2;
    }
    
    return str;
}

static inline uint8_t* lic_cipher_hex2raw(char* str)
{
    // TODO Does fail for len == 0 and uneven strlen(str).
    
    size_t len = strlen(str) / 2;
    uint8_t* cipher = (uint8_t*)malloc(len);
    
    if (!cipher)
    {
        // TODO Error handling.
    }
    
    unsigned int val;
    char* src = str;
    uint8_t* dst = cipher;
    char buf[3] = { 0, 0, 0 };
    while (*src)
    {
        buf[0] = *(src++);
        buf[1] = *(src++);
        
        sscanf(buf, "%02x", &val);
        
        *(dst)++ = val;
    }
    
    return cipher;
}

static inline lic_status_t lic_chksum(char* str, lic_chksum_t* cs, char** off)
{
    cs->sum = 0;
    cs->alt = 0;
    char* s = str;
    while (*(s + 3))
    {
        if (!((*s >= '0' && *s <= '9') ||
              (*s >= 'a' && *s <= 'z')))
            return LICENSE_INVFMT;
        
        cs->sum += (uint8_t)*s;
        cs->alt ^= (uint8_t)*s;
        
        s++;
    }
    cs->sum = (cs->sum ^ 255) + 1;
    cs->alt = ((cs->alt & 0xF0) >> 4) ^ (cs->alt & 0x0F);
    
    if (off)
        *off = s;
    
    return LICENSE_OK;
}

lic_status_t lic_valid_fmt1(char* lstr, char** lcore)
{
    // A123456789ASS
    
    // Skip format:
    lstr++;
    
    size_t len = strlen(lstr);
    
    if (len != 12)
        return LICENSE_INVFMT;
    
    char* s;
    lic_chksum_t cs;
    lic_status_t chk = lic_chksum(lstr, &cs, &s);
    if (chk != LICENSE_OK)
        return chk;
    
    // Verify integrity:
    char ref[4];
    sprintf(ref, "%01X%02X", cs.alt, cs.sum);
    if (strcmp(ref, s))
        return LICENSE_INT;
    
    if (lcore)
    {
        *lcore = strndup(&lstr[0], 9);
        
        if (!*lcore)
        {
            // TODO Error handling.
        }
    }

    return LICENSE_OK;
}

size_t function(char* ptr, size_t size, size_t nmemb, void* userdata)
{
    size_t sz = size * nmemb;
    lic_status_t* status_ptr = (lic_status_t*)userdata;
    
    // Minimum length required for: {"valid":true}
    if (sz < 14)
    {
        *status_ptr = LICENSE_NET;
        
        return sz;
    }
    
    // Check whether the license is valid:
    if (!strncmp((char*)ptr, "{\"valid\":true}", 14) ||
        !strncmp((char*)ptr, "{\"valid\": true}", 15))
    {
        *status_ptr = LICENSE_OK;
        
        return sz;
    }
    
    
    return sz;
}

lic_status_t lic_valid_onln(char* lstr, gen_fstat* stat)
{
    CURL *curl = curl_easy_init();
    
    if (!curl)
    {
        // TODO Error handling.
        
    }
    
    unsigned char iv[LIC_IVLEN];
    if (RAND_bytes(iv, LIC_IVLEN) != 1)
    {
        // TODO Error handling.
    }
    
    char* symkey;
    size_t symlen;
    lic_dec_sym((unsigned char*)EXE_MEMKEY, (unsigned char*)EXE_MEMIV, (char*)EXE_SYMKEY1, EXE_SYMLEN1, &symkey, &symlen);
    // TODO Error handling.
    
    uint8_t* cipher;
    size_t clen;
    lic_enc_sym((unsigned char*)symkey, iv, (unsigned char*)lstr, strlen(lstr), &cipher, &clen);
    
    // TODO Error handling.
    
    char* chex = lic_cipher_raw2hex(cipher, clen);
    
    free(cipher);
    free(symkey);
    
    char* ivhex = lic_cipher_raw2hex(iv, LIC_IVLEN);
    
    if (!ivhex)
    {
        // TODO Error handling.
    }
    
    char* sstr;
    if (stat)
    {
        // 10 digits for 32 bit
        // 20 digits for 64 bit
        sstr = (char*)malloc(67 + 10 + 10 + 10 + 20 + 1);
        
        if (!sstr)
        {
            // TODO Error handling.
        }
        
        sprintf(sstr, "{\"stat-comments\":%u, \"stat-features\":%u, \"stat-meta\":%u, \"stat-offset\":%llu}", stat->comms, stat->ftrs, stat->meta, stat->off);

        lic_dec_sym((unsigned char*)EXE_MEMKEY, (unsigned char*)EXE_MEMIV, (char*)EXE_SYMKEY1, EXE_SYMLEN1, &symkey, &symlen);
        // TODO Error handling.
        
        uint8_t* scipher;
        size_t sclen;
        lic_enc_sym((unsigned char*)symkey, iv, sstr, strlen(sstr), &scipher, &sclen);
        
        // TODO Error handling.
        
        char* schex = lic_cipher_raw2hex(scipher, sclen);
        
        free(sstr);
        free(scipher);
        free(symkey);
        
        sstr = schex;
    }
    else
        sstr = "";
    
    // Assume stats are always included:
    char* post = (char*)malloc(strlen(EXE_SYMID) + strlen(ivhex) + strlen(chex) + strlen(sstr) + 59 + 1);
    
    if (!post)
    {
        // TODO Error handling.
    }
    
    if (stat)
    {
        sprintf(post, "{\"encoding\":\"%s\", \"iv\":\"%s\", \"license\":\"%s\", \"supplementary\":\"%s\"}", EXE_SYMID, ivhex, chex, sstr);
    }
    else
    {
        sprintf(post, "{\"encoding\":\"%s\", \"iv\":\"%s\", \"license\":\"%s\"}", EXE_SYMID, ivhex, chex);
    }
    
    lic_status_t status;
    CURLcode res;
    curl_easy_setopt(curl, CURLOPT_URL, "http://45.56.107.233/license/");
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, post);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, function);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &status);
    res = curl_easy_perform(curl);
    
    if (res != CURLE_OK)
    {
        // TODO Error handling.
        status = LICENSE_NET;
    }
    
    curl_easy_cleanup(curl);
    
    free(post);
    
    return status;
}

void lic_free(lic_id* lid)
{
    
}

lic_status_t lic_valid(char* lstr, gen_fstat* stat)
{
    if (!*lstr)
        return LICENSE_INVFMT;
    
    // Note: Do not bother checking online, if license is not
    //       properly formatted/encoded to start with.
    lic_status_t status;
    switch (*lstr)
    {
        case 'A':
            return lic_valid_fmt1(lstr, NULL);
        case 'B':
            status = lic_valid_fmt1(lstr, NULL);
            
            if (status != LICENSE_OK)
                return status;
            
            return lic_valid_onln(lstr, stat);
        default:
            return LICENSE_INVFMT;
    }
}

// Code below from: https://github.com/saju/misc/blob/master/misc/openssl_aes.c

/**
 AES encryption/decryption demo program using OpenSSL EVP apis
 gcc -Wall openssl_aes.c -lcrypto
 this is public domain code.
 Saju Pillai (saju.pillai@gmail.com)
 **/

/**
 * Create an 256 bit key and IV using the supplied key_data. salt can be added for taste.
 * Fills in the encryption and decryption ctx objects and returns 0 on success
 **/
int aes_init(unsigned char *key_data, int key_data_len, unsigned char *salt, EVP_CIPHER_CTX *e_ctx,
             EVP_CIPHER_CTX *d_ctx)
{
    int i, nrounds = 5;
    unsigned char key[32], iv[16];
    
    /*
     * Gen key & IV for AES 256 CBC mode. A SHA1 digest is used to hash the supplied key material.
     * nrounds is the number of times the we hash the material. More rounds are more secure but
     * slower.
     */
    i = EVP_BytesToKey(EVP_aes_256_cbc(), EVP_sha1(), salt, key_data, key_data_len, nrounds, key, iv);
    if (i != 32) {
        printf("Key size is %d bits - should be 256 bits\n", i);
        return -1;
    }
    
    EVP_CIPHER_CTX_init(e_ctx);
    EVP_EncryptInit_ex(e_ctx, EVP_aes_256_cbc(), NULL, key, iv);
    EVP_CIPHER_CTX_init(d_ctx);
    EVP_DecryptInit_ex(d_ctx, EVP_aes_256_cbc(), NULL, key, iv);
    
    printf("KEY: ");
    for (i = 0; i < 32; i++)
        printf("\\x%02x", key[i]);
    printf("\n");

    printf("IV : ");
    for (i = 0; i < 16; i++)
        printf("\\x%02x", iv[i]);
    printf("\n");
    
    return 0;
}

int aes_encrypt(unsigned char *plaintext, int plaintext_len, unsigned char *key, unsigned char *iv, unsigned char *ciphertext)
{
    EVP_CIPHER_CTX *ctx;
    
    int len;
    
    int ciphertext_len;
    
    /* Create and initialise the context */
    if(!(ctx = EVP_CIPHER_CTX_new()))
    {
        // TODO Error handling.
    }
    
    /* Initialise the encryption operation. IMPORTANT - ensure you use a key
     * and IV size appropriate for your cipher
     * In this example we are using 256 bit AES (i.e. a 256 bit key). The
     * IV size for *most* modes is the same as the block size. For AES this
     * is 128 bits */
    if(1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv))
    {
        // TODO Error handling.
    }
    
    /* Provide the message to be encrypted, and obtain the encrypted output.
     * EVP_EncryptUpdate can be called multiple times if necessary
     */
    if(1 != EVP_EncryptUpdate(ctx, ciphertext, &len, plaintext, plaintext_len))
    {
        // TODO Error handling.
    }
    
    ciphertext_len = len;
    
    /* Finalise the encryption. Further ciphertext bytes may be written at
     * this stage.
     */
    if(1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len))
    {
        // TODO Error handling.
    }
    
    ciphertext_len += len;
    
    /* Clean up */
    EVP_CIPHER_CTX_free(ctx);
    
    return ciphertext_len;
}

int aes_decrypt(unsigned char *ciphertext, int ciphertext_len, unsigned char *key, unsigned char *iv, unsigned char *plaintext)
{
    EVP_CIPHER_CTX *ctx;
    
    int len;
    
    int plaintext_len;
    
    /* Create and initialise the context */
    if(!(ctx = EVP_CIPHER_CTX_new()))
    {
        // TODO Error handling.
    }
    
    /* Initialise the decryption operation. IMPORTANT - ensure you use a key
     * and IV size appropriate for your cipher
     * In this example we are using 256 bit AES (i.e. a 256 bit key). The
     * IV size for *most* modes is the same as the block size. For AES this
     * is 128 bits */
    if(1 != EVP_DecryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv))
    {
        // TODO Error handling.
    }
    
    /* Provide the message to be decrypted, and obtain the plaintext output.
     * EVP_DecryptUpdate can be called multiple times if necessary
     */
    if(1 != EVP_DecryptUpdate(ctx, plaintext, &len, ciphertext, ciphertext_len))
    {
        // TODO Error handling.
    }
    
    plaintext_len = len;
    
    /* Finalise the decryption. Further plaintext bytes may be written at
     * this stage.
     */
    if(1 != EVP_DecryptFinal_ex(ctx, plaintext + len, &len))
    {
        // TODO Error handling.
    }
    
    plaintext_len += len;
    
    /* Clean up */
    EVP_CIPHER_CTX_free(ctx);
    
    return plaintext_len;
}

lic_enc_t lic_enc_sym(unsigned char* key, unsigned char* iv, const char* plaintext, size_t len, uint8_t** ciphertext, size_t* cipherlen)
{
    /* Buffer for ciphertext. Ensure the buffer is long enough for the
     * ciphertext which may be longer than the plaintext, dependant on the
     * algorithm and mode
     */
    unsigned char* cipher = (unsigned char*)malloc((len + 16) & ~15);
    
    int clen;
    
    // TODO I get a linker error when not commenting the following out:
    /* Initialise the library */
    //ERR_load_crypto_strings();
    //OpenSSL_add_all_algorithms();
    //OPENSSL_config(NULL);
    
    /* Encrypt the plaintext */
    clen = aes_encrypt(plaintext, len, key, iv, cipher);
    
    /* Do something useful with the ciphertext here */
    printf("Ciphertext is:\n");
    BIO_dump_fp(stdout, cipher, clen);
    
    /* Clean up */
    EVP_cleanup();
    ERR_free_strings();
    
    *ciphertext = cipher;
    *cipherlen = clen;
    
    return ENCODING_OK;
}

lic_dec_t lic_dec_sym(unsigned char* key, unsigned char* iv, const char* cipher, size_t len, char** plaintext, size_t* plainlen)
{
    // Slight overallocation:
    *plaintext = (char*)malloc(len + 1);
    
    /* Decrypt the ciphertext */
    *plainlen = aes_decrypt(cipher, len, key, iv, (unsigned char*)*plaintext);
    
    // In case this is a C string:
    (*plaintext)[*plainlen] = 0;
    
    /* Clean up */
    EVP_cleanup();
    ERR_free_strings();
    
    return DECODING_OK;
}
