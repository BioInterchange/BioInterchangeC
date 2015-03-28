//
//  gff.h
//  biointerchange
//
//  Created by Joachim Baran on 2015-03-27.
//
//

#ifndef __biointerchange__gff__
#define __biointerchange__gff__


#include <stdint.h>
#include <stdio.h>
#include <sys/mman.h>

#define BI_IDX_PG 8192

typedef enum
{
    BI_GFF_HEADER,
    BI_GFF_FEATURE,
    BI_GFF_FASTA
} bi_gff_state;

#endif /* defined(__biointerchange__gff__) */
