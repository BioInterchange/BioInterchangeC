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

#include "gen.h"

const char* JSONLD_CTX = "@context";

/*
{
  "@context" :
  {
    "ID" :
    {
        "@ID" : "ID",
        "@type" : "http://www.biointerchange.org/gfvo#Identifier"
    }
  }
}
 */
const char* JSONLD_GFF3 = "http://www.biointerchange.org/jsonld/gff3.json";

const char* JSONLD_GTF = "http://www.biointerchange.org/jsonld/gtf.json";
const char* JSONLD_GVF = "http://www.biointerchange.org/jsonld/gvf.json";
const char* JSONLD_VCF = "http://www.biointerchange.org/jsonld/vcf.json";

void gen_xcig(char* str)
{
    char* wptr = str;
    
    char c = 0;
    while (*str)
    {
        if (*str >= '0' && *str <= '9')
        {
            *(wptr++) = *str;
        }
        else if (*str == ' ' && c)
        {
            *(wptr++) = c;
            c = 0;
        }
        else
        {
            if (c)
                *(wptr++) = c;
            
            c = *str;
        }
        
        str++;
    }
    
    if (c)
        *(wptr++) = c;
    
    *wptr = 0;
}

static inline bool gen_escchr(char* cptr)
{
    if (*cptr == '\n' ||
        *cptr == '"' ||
        *cptr == '\r' ||
        *cptr == '\t' ||
        *cptr == '\b' ||
        *cptr == '\f')
        return true;
    
    return false;
}

char* gen_escstr(char* str)
{
    size_t escchr = 0;
    
    // Count number of characters that need escaping:
    char* ptr = str;
    while (*ptr)
        if (gen_escchr(ptr++))
            escchr++;
    
    // Walk backwards and remove trailing newlines/carriage returns:
    while (ptr-- > str)
        if (*ptr == '\n' ||
            *ptr == '\r')
            *ptr = 0;
        else
            break;
    
    str = realloc(str, strlen(str) + escchr + 1);
    
    // TODO Error handling.
    
    // Escape characters:
    char c;
    ptr = str;
    while ((c = *ptr))
        if (gen_escchr(ptr))
        {
            *(ptr++) = '\\';
            
            switch (c)
            {
                case '"':
                    *(ptr++) = '"';
                    break;
                case '\n':
                    *(ptr++) = 'n';
                    break;
                case '\r':
                    *(ptr++) = 'r';
                    break;
                case '\t':
                    *(ptr++) = 't';
                    break;
                case '\b':
                    *(ptr++) = 'b';
                    break;
                case '\f':
                    *(ptr++) = 'f';
                    break;
                default:
                    // TODO Error handling. Internal error.
                    break;
            }
        }
        else
            *(ptr++) = c;
    *ptr = 0;
    
    return str;
}
