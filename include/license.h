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

#ifndef biointerchange_license_h
#define biointerchange_license_h

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <openssl/sha.h>

typedef enum
{
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
     * License valid.
     */
    LICENSE_OK
} lic_status_t;


typedef struct lic_id
{
    lic_status_t status;
    char* lid;
} lic_id;

#endif
