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

#ifndef biointerchange_gen_h
#define biointerchange_gen_h

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef __cplusplus
extern "C" {
#endif

extern const char* JSONLD_CTX;
extern const char* JSONLD_GFF3;
extern const char* JSONLD_GTF;
extern const char* JSONLD_GVF;
extern const char* JSONLD_VCF;

typedef struct gen_fstat
{
    /**
     * Number of meta-information lines.
     */
    uint32_t meta;
    /**
     * Number of feature lines.
     */
    uint32_t ftrs;
    /**
     * Number of comment lines.
     */
    uint32_t comms;
    /**
     * Offset to which the statistics are valid.
     */
    off_t off;
} gen_fstat;

typedef struct gen_fa_ntry_t
{
    off_t off;
    off_t seq;
    size_t llen;
    size_t lbrk;
} gen_fa_ntry_t;

/**
 * @brief Parser state.
 */
typedef struct gen_prsr_t
{
    bool fa_sct;
} gen_prsr_t;

void gen_xcig(char* str);

char* gen_escstr(char* str);
    
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
