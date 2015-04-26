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

#include "vcf.h"

void vcf_splt_inf(char* inf)
{
    char* ntry = inf;
    while (*inf)
    {
        if (*inf == ':')
        {
            // Split to separate string:
            *inf = 0;
            
        }
        
        inf++;
    }
}
