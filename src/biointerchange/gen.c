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