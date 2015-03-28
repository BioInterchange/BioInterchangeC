//
//  vcf.c
//  biointerchange
//
//  Created by Joachim Baran on 2015-03-27.
//
//

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
