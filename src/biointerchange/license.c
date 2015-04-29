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

static lic_id LIC_INVFMT = { LICENSE_INVFMT, NULL };
static lic_id LIC_INT = { LICENSE_INT, NULL };

lic_id* lic_valid_fmt1(char* lstr)
{
    // Format: AXXXXXXXXASS
    
    // Skip format:
    lstr++;
    
    size_t len = strlen(lstr);
    
    if (len != 10)
        return &LIC_INVFMT;
    
    uint8_t sum = 0;
    uint16_t alt = 0;
    char* s = lstr;
    while (*(s + 3))
    {
        if (*s < '0' && *s > '9')
            return &LIC_INVFMT;
        
        sum += (uint8_t)*s;
        alt ^= (uint8_t)*s;
        
        s++;
    }
    sum = (sum ^ 255) + 1;
    alt = ((alt & 0xF0) >> 4) ^ (alt & 0x0F);
    
    // Verify integrity:
    char ref[4];
    sprintf(ref, "%01X%02X", alt, sum);
    if (!strcmp(ref, s))
        return &LIC_INT;
    
    lic_id* lid = (lic_id*)malloc(sizeof(lic_id));
    
    if (!lid)
    {
        // TODO Error handling.
    }
    
    lid->status = LICENSE_OK;
    lid->lid = strndup(&lstr[1], 8);
    
    return lid;
}

lic_id* lic_valid_onln(lic_id* lid)
{
    
    return lid;
}

void lic_free(lic_id* lid)
{
    
}

lic_id* lic_valid(char* lstr)
{
    if (!*lstr)
        return LICENSE_INVFMT;
    
    lic_id* lid;
    switch (*lstr)
    {
        case 'A':
            return lic_valid_fmt1(lstr);
        case 'B':
            lid = lic_valid_fmt1(lstr);
            
            if (lid->status != LICENSE_OK)
                return lid;
            
            return lic_valid_onln(lid);
        default:
            return &LIC_INVFMT;
    }
}
