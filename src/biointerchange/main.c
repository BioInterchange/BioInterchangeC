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

#include <stdio.h>

#include "fio.h"
#include "gff.h"

int main()
{
    printf("Hello BioInterchange!\n");
    
    int fd = fio_opn("../BioInterchange/examples/chromosome_BF.gff");
    
    fio_cls(fd);
    
    printf("Bye!\n");
    
    return 0;
}
