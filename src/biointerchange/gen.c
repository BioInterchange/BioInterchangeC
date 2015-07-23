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
#include "gvf.h"
#include "vcf.h"

const char* JSONLD_CTX = "@context";
const char* JSONLD_TPE = "@type";

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

const char* JSONLD_CLSS_CTX = "https://www.codamono.com/gfvo-squared#Context";
const char* JSONLD_CLSS_MTA = "https://www.codamono.com/gfvo-squared#Meta";
const char* JSONLD_CLSS_FTR = "https://www.codamono.com/gfvo-squared#Feature";
const char* JSONLD_CLSS_SUM = "https://www.codamono.com/gfvo-squared#Summary";

const char* JSONLD_GFF3_CTX1 = "https://www.codamono.com/jsonld/gff3-c1.json";
const char* JSONLD_GTF_CTX1 = "https://www.codamono.com/jsonld/gtf-c1.json";
const char* JSONLD_GVF_CTX1 = "https://www.codamono.com/jsonld/gvf-c1.json";
const char* JSONLD_VCF_CTX1 = "https://www.codamono.com/jsonld/vcf-c1.json";

const char* JSONLD_GFF3_X1 = "https://www.codamono.com/jsonld/gff3-x1.json";
const char* JSONLD_GTF_X1 = "https://www.codamono.com/jsonld/gtf-x1.json";
const char* JSONLD_GVF_X1 = "https://www.codamono.com/jsonld/gvf-x1.json";
const char* JSONLD_VCF_X1 = "https://www.codamono.com/jsonld/vcf-x1.json";

const char* JSONLD_GFF3_1 = "https://www.codamono.com/jsonld/gff3-f1.json";
const char* JSONLD_GTF_1 = "https://www.codamono.com/jsonld/gtf-f1.json";
const char* JSONLD_GVF_1 = "https://www.codamono.com/jsonld/gvf-f1.json";
const char* JSONLD_VCF_1 = "https://www.codamono.com/jsonld/vcf-f1.json";

const char* JSONLD_STAT_1 = "https://www.codamono.com/jsonld/biointerchange-s1.json";

const char* GEN_AFFECTED = "affected-features";         // (documented)
const char* GEN_AFFECTED_TPE = "affected-feature-type"; // (documented)
const char* GEN_ALIAS = "alias";                        // (documented)
const char* GEN_ALIAS_GFF3 = "Alias";                   // N/A
const char* GEN_ALIGNMENT = "alignment";                // (documented)
const char* GEN_ALIGNMENT_GFF3 = "Target";              // N/A
const char* GEN_ALLELE_CNT = "allele-count";            // (documented)
const char* GEN_ALLELE_CNT_VCF = "AC";                  // N/A
const char* GEN_ALLELE_CNTEXP = "allele-count-expected"; // (documented)
const char* GEN_ALLELE_CNTEXP_VCF = "EC";               // N/A
const char* GEN_ALLELE_FRQ = "allele-frequency";        // (documented)
const char* GEN_ALLELE_FRQ_VCF = "AF";                  // N/A
const char* GEN_ALLELE_TTL = "allele-total-number";     // (documented)
const char* GEN_ALLELE_TTL_VCF = "AN";                  // N/A
const char* GEN_ANCESTRAL_ALLELE = "ancestral-allele";  // (documented)
const char* GEN_ANCESTRAL_ALLELE_VCF = "AA";            // N/A
const char* GEN_ANNOTATIONS = "annotations";            // (documented)
const char* GEN_ANNOTATIONS_VCF = "FT";
const char* GEN_ATTRIBUTE_MTHD = "attribute-method";    // (documented)
const char* GEN_ATTRS = "user-defined";                 // (documented)
const char* GEN_AVG_COVERAGE = "average-coverage";      // (documented)
const char* GEN_AVG_COVERAGE_GVF = "Average_coverage";  // N/A
const char* GEN_BUILD = "genome-build";                 // (documented; takes value from line -- not this constant)
const char* GEN_BUILD_GFF3 = "genome-build";                 // N/A
const char* GEN_BUILD_VAL = "build";                    // (documented)
const char* GEN_CIGAR = "cigar-string";                 // (documented)
const char* GEN_CIGAR_GFF3 = "Gap";                     // N/A
const char* GEN_CIGAR_VCF = "CIGAR";                    // N/A
const char* GEN_CODON = "codon";                        // (documented)
const char* GEN_CODON_PHASE = "codon-phase";            // (documented)
const char* GEN_COMMENT = "comment";                    // (documented)
const char* GEN_COMMENT_GVF = "Comment";
const char* GEN_DATA_SRC = "data-source";               // (documented; takes value from line -- see `gvf_proc_prgm`)
const char* GEN_DATA_TPE = "data-type";                 // (documented; takes value from line -- see `gvf_tags`)
const char* GEN_DATA_TPE_GVF = "Data_type";
const char* GEN_DBXREF = "dbxref";                      // (documented)
const char* GEN_DBXREF_GFF3 = "Dbxref";                 // N/A
const char* GEN_DEPTH = "depth";                        // (documented)
const char* GEN_DEPTH_VCF = "DP";                       // N/A
const char* GEN_EFFECT = "effect";                      // (documented)
const char* GEN_EFFECTS = "effects";                    // (documented)
const char* GEN_END = "end";                            // (documented)
const char* GEN_END_VCF = "END";
const char* GEN_FILE_DATE = "file-date";                // (documented)
const char* GEN_FILE_DATE_VCF1 = "filedate";
const char* GEN_FILE_DATE_VCF2 = "fileDate";
const char* GEN_FMT_ALT = "allele-format";              // (documented)
const char* GEN_FMT_ALT_VCF = "ALT";
const char* GEN_FMT_INFO = "feature-format";            // (documented)
const char* GEN_FMT_INFO_VCF = "INFO";
const char* GEN_FMT_FLT = "annotation-format";          // (documented)
const char* GEN_FMT_FLT_VCF = "FILTER";
const char* GEN_FMT_GT = "genotype-format";             // (documented)
const char* GEN_FMT_GT_VCF = "FORMAT";
const char* GEN_FST_BKPNT = "breakpoint-fasta";         // (documented)
const char* GEN_FST_BKPNT_VCF = "assembly";
const char* GEN_FST_REF = "reference-fasta";            // (documented; code left untouched in GVF)
const char* GEN_FST_REF_GVF = "reference-fasta";
const char* GEN_FST_REF_VCF = "reference";
const char* GEN_GENOMIC_SRC = "genomic-source";         // (documented)
const char* GEN_GENOTYPE = "genotype";                  // (documented)
const char* GEN_GENOTYPE_VCF = "GT";
const char* GEN_GENOTYPE_LIKE = "genotype-likelihood";  // (documented)
const char* GEN_GENOTYPE_LIKE_VCF = "GL";
const char* GEN_GENOTYPE_LIKEP = "genotype-likelihood-phred-scaled"; // (documented)
const char* GEN_GENOTYPE_LIKEP_VCF = "PL";
const char* GEN_GENOTYPE_QUAL = "genotype-quality";     //
const char* GEN_GENOTYPE_QUAL_VCF = "GQ";               // N/A
const char* GEN_GENOTYPE_PROB = "genotype-probabilities-phred-scaled";
const char* GEN_GENOTYPE_PROB_VCF = "GP";
const char* GEN_GFFVERSION = "gff-version";             // (documented)
const char* GEN_GFFVERSION_GFF3 = "gff-version";
const char* GEN_GLOBAL = "global";                      // (documented)
const char* GEN_GVFVERSION = "gvf-version";             // (documented)
const char* GEN_GVFVERSION_GVF = "gvf-version";
const char* GEN_HAP_QUALITIES = "haplotype-qualities";
const char* GEN_HAP_QUALITIES_VCF = "HQ";
const char* GEN_ID = "id";                              // (documented)
const char* GEN_ID_GFF3 = "ID";                         // N/A
const char* GEN_INDIVIDUALS = "individuals";            // (documented)
const char* GEN_INDIVIDUALS_GVF1 = "individual-id";     // N/A
const char* GEN_INDIVIDUALS_GVF2 = "multi-individual";  // N/A
const char* GEN_LANDMARK = "landmark";                  // (documented)
const char* GEN_LANDMARKS = "landmarks";                // (documented)
const char* GEN_LANDMARKS_GVF = "Seqid";
const char* GEN_LOCUS = "locus";                        // (documented)
const char* GEN_MEMBER_1000G = "membership-1000G";      // (documented)
const char* GEN_MEMBER_1000G_VCF = "1000G";
const char* GEN_MEMBER_DBSNP = "membership-dbsnp";      // (documented)
const char* GEN_MEMBER_DBSNP_VCF = "DB";
const char* GEN_MEMBER_HM2 = "membership-hapmap-2";     // (documented)
const char* GEN_MEMBER_HM2_VCF = "H2";
const char* GEN_MEMBER_HM3 = "membership-hapmap-3";     // (documented)
const char* GEN_MEMBER_HM3_VCF = "H3";
const char* GEN_ONT_ACCESSION = "ontology-accession";   // not used?
const char* GEN_ONT_TERM = "ontology-term";             // (documented)
const char* GEN_PHASE_SET = "phase-set";                // (documented)
const char* GEN_PHASE_SET_VCF = "PS";
const char* GEN_PHASE_QUAL = "phasing-quality";         // (documented)
const char* GEN_PHASE_QUAL_VCF = "PQ";
const char* GEN_PHASED_GENO = "phased-genotypes";       // (documented)
const char* GEN_PHENO_DESCR = "phenotype-description";  // (documented)
const char* GEN_POP = "population";                     // (documented)
const char* GEN_QUALITY_MAP = "mapping-quality-rms";    // (documented)
const char* GEN_QUALITY_MAP_VCF = "MQ";
const char* GEN_QUALITY_MAP0 = "reads-with-zero-mapping-quality"; // (documented)
const char* GEN_QUALITY_MAP0_VCF = "MQ0";
const char* GEN_QUALITY_RMS = "base-quality-rms";       // (documented)
const char* GEN_QUALITY_RMS_VCF = "BQ";
const char* GEN_READ_PAIR_SPAN = "read-pair-span";      // (documented)
const char* GEN_READ_PAIR_SPAN_GVF = "Read_pair_span";
const char* GEN_REFERENCE = "reference";                // (documented)
const char* GEN_SAMPLES_DATA = "samples-with-data";     // (documented)
const char* GEN_SAMPLES_DATA_VCF = "NS";
const char* GEN_SCORE_MTHD = "score-method";            // (documented)
const char* GEN_SEQUENCE = "sequence";                  // (documented)
const char* GEN_SEQUENCE_REGION = "contig";             // (documented)
const char* GEN_SEQUENCE_REGION_GFF3 = "sequence-region";
const char* GEN_SEQUENCE_REGION_VCF = "contig";
const char* GEN_SEQUENCES = "sequences";                // (documented)
const char* GEN_SEX = "sex";                            // (documented)
const char* GEN_START = "start";                        // (documented)
const char* GEN_STRAND = "strand";                      // (documented)
const char* GEN_STRAND_BIAS = "strand-bias";            // (documented)
const char* GEN_STRAND_BIAS_VCF = "SB";
const char* GEN_SOMATIC = "somatic-mutation";           // (documented)
const char* GEN_SOMATIC_VCF = "SOMATIC";
const char* GEN_SOURCE = "source";                      // (documented)
const char* GEN_SOURCES = "sources";                    // (documented)
const char* GEN_SOURCES_GVF = "Source";
const char* GEN_SOURCE_MTHD = "source-method";          // (documented)
const char* GEN_TECHNOLOGY = "technology-platform";     // (documented; takes value from line -- see `gvf_proc_prgm`)
const char* GEN_TYPE = "type";                          // (documented)
const char* GEN_TYPE_GVF = "Type";
const char* GEN_TYPES = "types";                        // (documented)
const char* GEN_TYPES_GVF = "Type";
const char* GEN_VALIDATED = "experimentally-validated"; // (documented)
const char* GEN_VALIDATED_VCF = "VALIDATED";
const char* GEN_VARIANTS = "variants";                  // (documented)
const char* GEN_VCFVERSION = "vcf-version";             // (documented)
const char* GEN_VCFVERSION_VCF = "fileformat";          // N/A

const char* GEN_EMPTY = "";
const char* GEN_NULL = "null";
const char* GEN_TRUE = "true";
const char* GEN_UNKNOWN = ".";
const char* GEN_FORWARD = "+";

const char* GEN_COUNT = "count"; // TODO Check whether better replaced with GEN_ALLELE_CNT
const char* GEN_FREQUENCY = "frequency";

const char* GEN_SEQUENCE_GVF = "seq";

// Statistics wrap-up:
const char* GEN_STATS = "statistics";                   //
const char* GEN_STAT_CMMS = "comment-lines";            //
const char* GEN_STAT_META = "meta-lines";               //
const char* GEN_STAT_META_FLTR = "meta-lines-filtered"; //
const char* GEN_STAT_FTRS = "features";                 //
const char* GEN_STAT_FTRS_FLTR = "features-filtered";   //
const char* GEN_STAT_RNTM = "runtime";                  //
const char* GEN_STAT_INVC = "invocation";               //
const char* GEN_STAT_FNSH = "finish";                   //
const char* GEN_STAT_LSEC = "lapsed-seconds";           //

ldoc_vis_nde_ord_t* json_vis_nde;
ldoc_vis_ent_t* json_vis_ent;

// Assumes VCF_MAX_ALT <= 26; one character for an allele.
char GEN_ALLELE[GEN_MAX_ALT * 2];

// Assumes VCF_MAX_ALT <= 26; one character for an allele, so encoding of
// AA, AB, AC, is only VCF_STEP characters wide (allele 1, allel 2, null terminator).
char GEN_ALLELES[GEN_STEP * (((GEN_MAX_ALT * (GEN_MAX_ALT + 1) / 2) + GEN_MAX_ALT + 1))];

static inline void gen_allele_lbl(char* s, size_t j, size_t k)
{
    // Note: requires j <= k
    // Based on: https://samtools.github.io/hts-specs/VCFv4.2.pdf
    // Index formula: F(j/k) = (k*(k+1)/2)+j
    // Bi-allelic example: AA, AB, BB
    // Tri-allelic example: AA, AB, BB, AC, BC, CC
    off_t idx = GEN_STEP * (( k * ( k + 1 ) / 2 ) + j);
    GEN_ALLELES[idx] = 'A' + j;
    GEN_ALLELES[idx + 1] = 'A' + k;
    GEN_ALLELES[idx + 2] = 0;
}

void gen_err(int err, const char* s)
{
    switch (err)
    {
        case MAIN_ERR_JFMT_KY:
            fprintf(stderr, "Invalid JSON document structure.\n");
            if (s)
                fprintf(stderr, "Cannot find key: %s\n", s);
            break;
        case MAIN_ERR_JFMT_KYMAP:
            fprintf(stderr, "Invalid JSON document structure.\n");
            if (s)
                fprintf(stderr, "Cannot find key for genomic attribute: %s\n", s);
            break;
        case MAIN_ERR_JFMT_DEP:
            fprintf(stderr, "Invalid JSON document structure.\n");
            if (s)
                fprintf(stderr, "Cause: %s\n", s);
            break;
        default:
            break;
    }
    
    exit(err);
}

void gen_init()
{
    json_vis_nde = ldoc_vis_nde_ord_new();
    json_vis_nde->vis_setup = ldoc_vis_setup_json;
    json_vis_nde->vis_teardown = ldoc_vis_teardown_json;
    ldoc_vis_nde_uni(&(json_vis_nde->pre), ldoc_vis_nde_pre_json);
    ldoc_vis_nde_uni(&(json_vis_nde->infx), ldoc_vis_nde_infx_json);
    ldoc_vis_nde_uni(&(json_vis_nde->post), ldoc_vis_nde_post_json);
    
    json_vis_ent = ldoc_vis_ent_new();
    ldoc_vis_ent_uni(json_vis_ent, ldoc_vis_ent_json);
    
    char* s = GEN_ALLELES;
    
    size_t j = 0;
    for (; j < GEN_MAX_ALT; j++)
    {
        GEN_ALLELE[2 * j] = 'A' + j;
        GEN_ALLELE[2 * j + 1] = 0;
        
        size_t k = 0;
        for (; k < GEN_MAX_ALT; k++)
            if (j <= k)
                gen_allele_lbl(s, j, k);
    }
}

inline void gen_nde_dsc_opt(ldoc_nde_t* nde, ldoc_nde_t* dsc, char* lbl)
{
    if (dsc)
    {
        ldoc_nde_dsc_push(nde, dsc);
        
        return;
    }
    
    ldoc_ent_t* dsc_null = ldoc_ent_new(LDOC_ENT_OR);
    
    dsc_null->pld.pair.anno.str = lbl;
    dsc_null->pld.pair.dtm.str = NULL;
    
    ldoc_nde_ent_push(nde, dsc_null);
}

inline char* gen_term_crnl(char* s)
{
    char* b = s;
    
    while (*s)
    {
        if (*s == '\n' || *s == '\r')
        {
            *s = 0;
            
            return b;
        }
        
        s++;
    }
    
    return b;
}

inline void gen_ky(char* attr, char** val)
{
    while (*attr)
    {
        if (*attr == '=')
        {
            *attr = 0;
            *val = ++attr;
            
            return;
        }
        else
            attr++;
    }
    
    // No value assignment; just a key.
    *val = NULL;
}

inline char* gen_res_opt(ldoc_res_t* res)
{
    if (!res)
        return (char*)GEN_UNKNOWN;
    
    if (res->nde)
    {
        // TODO Internal error.
    }
    
    // TODO Should handle various entity types?
    if (res->info.ent->pld.pair.dtm.str)
        return res->info.ent->pld.pair.dtm.str;
    
    return (char*)GEN_UNKNOWN;
}

inline char* gen_res_optx(ldoc_res_t* res)
{
    if (!res)
        return (char*)GEN_EMPTY;
    
    return gen_res_req(res);
}

char* gen_res_req(ldoc_res_t* res)
{
    // TODO Does this cover all used entity types?
    if (!res || res->nde || res->info.ent->pld.pair.dtm.str)
    {
        // TODO Data error.
    }
    
    return res->info.ent->pld.pair.dtm.str;
}

inline void gen_lwr(char* str)
{
    while (*str)
    {
        if (*str >= 'A' && *str <= 'Z')
            *str = *str - 'A' + 'a';
        
        str++;
    }
}

inline void gen_lwrhyph(char* str)
{
    while (*str)
    {
        if (*str >= 'A' && *str <= 'Z')
            *str = *str - 'A' + 'a';
        else if (*str == '_')
            *str = '-';
        
        str++;
    }
}

inline char* gen_quoskp(char* str)
{
    // Does not start with a quote: bail out
    if (*str != '"')
        return str;

    // Seen a quote: proceed to end and remove last character,
    //               which is assumed to be a quote too!
    str++;
    char* c = str;
    while (*c)
        c++;
    if (c > str)
        *(c - 1) = 0;
    
    return str;
}

inline void gen_kwd(char* str, gen_attr_t* kwd, bi_attr upfail)
{
    if (*str >= 'A' && *str <= 'Z')
    {
        if (*str == 'A')
        {
            str++;
            if (*str == 'A')
            {
                str++;
                if (!*str)
                {
                    // VCF: AA, ancestral allele
                    kwd->attr = BI_VAL;
                    kwd->alt = GEN_ANCESTRAL_ALLELE;
                    return;
                }
            }
            if (*str == 'C')
            {
                str++;
                if (!*str)
                {
                    // VCF: AC, allele count
                    kwd->attr = BI_CSEPVARN;
                    kwd->alt = GEN_ALLELE_CNT;
                    return;
                }
            }
            else if (*str == 'F')
            {
                str++;
                if (!*str)
                {
                    // VCF: AF, allele frequency
                    kwd->attr = BI_CSEPVARN;
                    kwd->alt = GEN_ALLELE_FRQ;
                    return;
                }
            }
            else if (*str == 'N')
            {
                str++;
                if (!*str)
                {
                    // VCF: AN, total number of alleles in called
                    kwd->attr = BI_NUM;
                    kwd->alt = GEN_ALLELE_TTL;
                    return;
                }
            }
            else if (*str == 'l')
            {
                str++;
                if (*str == 'i')
                {
                    str++;
                    if (!strcmp(str, "as"))
                    {
                        // GFF3: Alias
                        kwd->attr = BI_CSEP;
                        kwd->alt = NULL;
                        return;
                    }
                }
            }
        }
        else if (*str == 'B')
        {
            str++;
            if (*str == 'Q')
            {
                str++;
                if (!*str)
                {
                    // VCF: BQ, RMS base quality
                    kwd->attr = BI_NUM;
                    kwd->alt = GEN_QUALITY_RMS;
                    return;
                }
            }
        }
        else if (*str == 'C')
        {
            str++;
            if (*str == 'I')
            {
                str++;
                if (*str == 'G')
                {
                    str++;
                    if (!strcmp(str, "AR"))
                    {
                        // VCF: CIGAR
                        kwd->attr = BI_XCIG;
                        kwd->alt = GEN_CIGAR;
                        return;
                    }
                }
            }
        }
        else if (*str == 'D')
        {
            str++;
            if (*str == 'B')
            {
                str++;
                if (!*str)
                {
                    // VCF: DB, dbSNP membership
                    kwd->attr = BI_BL;
                    kwd->alt = GEN_MEMBER_DBSNP;
                    return;
                }
            }
            if (*str == 'P')
            {
                str++;
                if (!*str)
                {
                    // VCF: DP, depth across samples
                    kwd->attr = BI_NUM;
                    kwd->alt = GEN_DEPTH;
                    return;
                }
            }
            else if (*str == 'b')
            {
                str++;
                if (*str == 'x')
                {
                    str++;
                    if (!strcmp(str, "ref"))
                    {
                        // GFF3/GVF: Dbxref
                        kwd->attr = BI_CSEPCPAIR;
                        kwd->alt = NULL;
                        return;
                    }
                }
            }
            else if (*str == 'e')
            {
                str++;
                if (*str == 'r')
                {
                    str++;
                    if (!strcmp(str, "ives_from"))
                    {
                        // GFF3:
                        kwd->attr = BI_VAL;
                        kwd->alt = NULL;
                        return;
                    }
                }
            }
        }
        else if (*str == 'E')
        {
            str++;
            if (*str == 'C')
            {
                str++;
                if (!*str)
                {
                    // VCF: EC, expected counts (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_ALLELE_CNTEXP;
                    return;
                }
            }
            else if (*str == 'N')
            {
                str++;
                if (*str == 'D')
                {
                    str++;
                    if (!*str)
                    {
                        // VCF: END, end position of variant description
                        kwd->attr = BI_NUM;
                        kwd->alt = NULL; // TODO Check whether this is correct.
                        return;
                    }
                }
            }
        }
        else if (*str == 'F')
        {
            str++;
            if (*str == 'T')
            {
                str++;
                if (!*str)
                {
                    // VCF: FT, filter (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_ANNOTATIONS;
                    return;
                }
            }
        }
        else if (*str == 'G')
        {
            str++;
            if (*str == 'L')
            {
                str++;
                if (!*str)
                {
                    // VCF: GL, genotype-likelihood (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_GENOTYPE_LIKE;
                    return;
                }
            }
            else if (*str == 'P')
            {
                str++;
                if (!*str)
                {
                    // VCF: GP, genotype-probabilities Phred-scaled (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_GENOTYPE_PROB;
                    return;
                }
            }
            else if (*str == 'Q')
            {
                str++;
                if (!*str)
                {
                    // VCF: GQ, genotype-quality (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_GENOTYPE_QUAL;
                    return;
                }
            }
            else if (*str == 'T')
            {
                str++;
                if (!*str)
                {
                    // VCF: GT, genotype (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_GENOTYPE;
                    return;
                }
            }
            else if (*str == 'a')
            {
                str++;
                if (*str == 'p')
                {
                    str++;
                    if (!*str)
                    {
                        // GFF3: Gap
                        kwd->attr = BI_XCIG;
                        kwd->alt = GEN_CIGAR;
                        return;
                    }
                }
            }
        }
        else if (*str == 'H')
        {
            str++;
            if (*str == '2')
            {
                str++;
                if (!*str)
                {
                    // VCF: H2, HapMap 2 membership
                    kwd->attr = BI_BL;
                    kwd->alt = GEN_MEMBER_HM2;
                    return;
                }
            }
            else if (*str == '3')
            {
                str++;
                if (!*str)
                {
                    // VCF: H3, HapMap 3 membership
                    kwd->attr = BI_BL;
                    kwd->alt = GEN_MEMBER_HM3;
                    return;
                }
            }
            else if (*str == 'Q')
            {
                str++;
                if (!*str)
                {
                    // VCF: HQ, haplotide quality (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_HAP_QUALITIES;
                    return;
                }
            }

        }
        else if (*str == 'I')
        {
            str++;
            if (*str == 'D')
            {
                str++;
                if (!*str)
                {
                    // GFF3/GVF: ID
                    kwd->attr = BI_VAL;
                    kwd->alt = GEN_ID;
                    return;
                }
            }
            else if(*str == 's')
            {
                str++;
                if (*str == '_')
                {
                    str++;
                    if (!strcmp(str, "circular"))
                    {
                        // GFF3: Is_circular
                        kwd->attr = BI_BL;
                        kwd->alt = NULL;
                        return;
                    }
                }
            }
        }
        else if (*str == 'M')
        {
            str++;
            if (*str == 'Q')
            {
                str++;
                if (!*str)
                {
                    // VCF: MQ, mapping quality
                    kwd->attr = BI_NUM;
                    kwd->alt = GEN_QUALITY_MAP;
                    return;
                }
                else if (*str == '0')
                {
                    str++;
                    if (!*str)
                    {
                        // VCF: MQ0, mapping quality == 0
                        kwd->attr = BI_NUM;
                        kwd->alt = GEN_QUALITY_MAP0;
                        return;
                    }
                }
            }
        }
        else if (*str == 'N')
        {
            str++;
            if (*str == 'S')
            {
                str++;
                if (!*str)
                {
                    // VCF: NS, number of samples with data
                    kwd->attr = BI_NUM;
                    kwd->alt = GEN_SAMPLES_DATA;
                    return;
                }
            }
            else if (*str == 'o')
            {
                str++;
                if (*str == 't')
                {
                    str++;
                    if (!strcmp(str, "e"))
                    {
                        // GFF3: Note
                        kwd->attr = BI_VAL;
                        kwd->alt = NULL;
                        return;
                    }
                }
            }
        }
        else if (*str == 'O')
        {
            str++;
            if (*str == 'n')
            {
                str++;
                if (*str == 't')
                {
                    str++;
                    if (!strcmp(str, "ology_term"))
                    {
                        // GFF3: Ontology_term
                        kwd->attr = BI_CSEP;
                        kwd->alt = GEN_ONT_TERM;
                        return;
                    }
                }
            }
        }
        else if (*str == 'P')
        {
            str++;
            if (*str == 'L')
            {
                str++;
                if (!*str)
                {
                    // VCF: PL, genotype-likelihood phred-scaled (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_GENOTYPE_LIKEP;
                    return;
                }
            }
            else if (*str == 'S')
            {
                str++;
                if (!*str)
                {
                    // VCF: PS, phase set (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_PHASE_SET;
                    return;
                }
            }
            else if (*str == 'Q')
            {
                str++;
                if (!*str)
                {
                    // VCF: PQ, phasing quality (genotype field)
                    kwd->attr = BI_NKW;
                    kwd->alt = (char*)GEN_PHASE_QUAL;
                    return;
                }
            }
            else if (*str == 'a')
            {
                str++;
                if (*str == 'r')
                {
                    str++;
                    if (!strcmp(str, "ent"))
                    {
                        // GFF3: Parent
                        kwd->attr = BI_CSEP;
                        kwd->alt = NULL;
                        return;
                    }
                }
            }
        }
        else if (*str == 'R')
        {
            str++;
            if (*str == 'e')
            {
                str++;
                if (*str == 'f')
                {
                    str++;
                    if (!strcmp(str, "erence_seq"))
                    {
                        // GVF: Reference_seq
                        kwd->attr = BI_IGN;
                        kwd->alt = NULL;
                        return;
                    }
                    else if (!strcmp(str, "erence_codon"))
                    {
                        // GVF: Reference_codon
                        kwd->attr = BI_REFSEQ10;
                        kwd->alt = NULL;
                        return;
                    }
                }
            }
        }
        else if (*str == 'S')
        {
            str++;
            if (*str == 'B')
            {
                str++;
                if (!*str)
                {
                    // VCF: SB, strand bias
                    kwd->attr = BI_NUM;
                    kwd->alt = GEN_STRAND_BIAS;
                    return;
                }
            }
            else if (*str == 'O')
            {
                str++;
                if (*str == 'M')
                {
                    str++;
                    if (!strcmp(str, "ATIC"))
                    {
                        // VCF: SOMATIC, somatic mutation
                        kwd->attr = BI_BL;
                        kwd->alt = GEN_SOMATIC;
                        return;
                    }
                }
            }
        }
        else if (*str == 'T')
        {
            str++;
            if (*str == 'a')
            {
                str++;
                if (*str == 'r')
                {
                    str++;
                    if (!strcmp(str, "get"))
                    {
                        // GFF3: Target
                        kwd->attr = BI_TARGET;
                        kwd->alt = NULL;
                        return;
                    }
                }
            }
        }
        else if (*str == 'V')
        {
            str++;
            if (*str == 'a')
            {
                str++;
                if (*str == 'r')
                {
                    str++;
                    if (!strcmp(str, "iant_seq"))
                    {
                        // GVF: Variant_seq
                        kwd->attr = BI_IGN;
                        kwd->alt = NULL;
                        return;
                    }
                    else if (!strcmp(str, "iant_effect"))
                    {
                        // GVF: Variant_effect
                        kwd->attr = BI_GVFEFFECT;
                        kwd->alt = NULL;
                        return;
                    }
                    else if (!strncmp(str, "iant_", 5))
                    {
                        // GVF: Variant_*
                        kwd->attr = BI_CSEPVAR8;
                        kwd->alt = NULL;
                        return;
                    }
                }
            }
            else if (*str == 'A')
            {
                str++;
                if (*str == 'L')
                {
                    str++;
                    if (!strcmp(str, "IDATED"))
                    {
                        // VCF: VALIDATED
                        kwd->attr = BI_BL;
                        kwd->alt = GEN_VALIDATED;
                        return;
                    }
                }
            }
        }
        
        kwd->attr = upfail;
        kwd->alt = NULL;
        return;
    }
    else if (*str == '1')
    {
        str++;
        if (*str == '0')
        {
            str++;
            if (*str == '0')
            {
                str++;
                if (!strcmp(str, "0G"))
                {
                    // VCF: 1000G, 1000-genomes membership
                    kwd->attr = BI_BL;
                    kwd->alt = GEN_MEMBER_1000G;
                    return;
                }
            }
        }
    }

    kwd->attr = BI_NKW;
    kwd->alt = NULL;
    return;
}

inline char gen_inv(char c)
{
    switch (c)
    {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'N':
            return 'N';
        case 'T':
            return 'A';
        default:
            // TODO Error.
            break;
    }
    
    return 0;
}

inline ldoc_content_t gen_smrt_tpe(char* val)
{
    if (!val)
        return LDOC_ENT_BR;
    
    if (*val == '.' && !*(val + 1))
        return LDOC_ENT_OR;
    
    bool dot = false;
    while (*val)
    {
        if (!(*val >= '0' && *val <= '9'))
        {
            if (*val == '.')
            {
                if (dot)
                    return LDOC_ENT_OR;
                else
                    dot = true;
            }
            else
                return LDOC_ENT_OR;
        }
        
        val++;
    }
    
    return LDOC_ENT_NR;
}

inline ldoc_content_t gen_smrt_flttpe(char* val)
{
    if (!val)
        return LDOC_ENT_BL;
    
    if (*val == '.' && !*(val + 1))
        return LDOC_ENT_TXT;
    
    bool dot = false;
    while (*val)
    {
        if (!(*val >= '0' && *val <= '9'))
        {
            if (*val == '.')
            {
                if (dot)
                    return LDOC_ENT_TXT;
                else
                    dot = true;
            }
            else
                return LDOC_ENT_TXT;
        }
        
        val++;
    }
    
    return LDOC_ENT_NUM;
}

void gen_xcig(char* str)
{
    // CIGAR format seems fine already (digits followed by operations):
    if (*str >= '0' && *str <= '9')
        return;
    
    char* wptr = str;
    
    // `c` stores the alignment operation ('M', 'D', etc.)
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

void gen_qk_revcig(char* str)
{
    // Pointer to where the count starts; used when an
    // alignment operation ('M', 'D', etc.) is found:
    char* n = str;
    uint8_t nlen = 0;
    while (*str)
    {
        if (*str >= '0' && *str <= '9')
            nlen++;
        else
        {
            // Cat: operation, then count, space as per GFF3 examples
            qk_strncat(str, 1);
            qk_strncat(n, nlen);
            if (*(str + 1))
                qk_strcat(" ");
            
            n = str + 1;
            nlen = 0;
        }
        
        str++;
    }
}

size_t gen_csplit(char* str, char c)
{
    size_t n = 0;
    
    while (*str)
    {
        if (*str == c)
        {
            *str = 0;
            n++;
        }
        
        str++;
    }
    
    return n;
}

inline ldoc_nde_t* gen_ctx(ldoc_nde_t* nde, bool* nw, const char* ctx_url)
{
    ldoc_nde_t* usr;
    
    if (!nde->ent_cnt)
    {
        ldoc_ent_t* ctx = ldoc_ent_new(LDOC_ENT_OR);
        
        // TODO Error handling.
        
        ctx->pld.pair.anno.str = (char*)JSONLD_CTX;
        ctx->pld.pair.dtm.str = (char*)ctx_url;
        ldoc_nde_ent_push(nde, ctx);

        ldoc_ent_t* tpe = ldoc_ent_new(LDOC_ENT_OR);
        
        // TODO Error handling.
        
        tpe->pld.pair.anno.str = (char*)JSONLD_TPE;
        tpe->pld.pair.dtm.str = (char*)JSONLD_CLSS_MTA;
        ldoc_nde_ent_push(nde, tpe);
        
        // User-defined pragmas:
        usr = ldoc_nde_new(LDOC_NDE_UA);
        usr->mkup.anno.str = strdup(GEN_ATTRS);
        
        if (nw)
            *nw = true;
    }
    else
    {
        const char* pth[1] = { GEN_ATTRS };
        ldoc_res_t* usr_ = ldoc_find_anno_nde(nde, (char**)pth, 1);
        
        if (!usr_)
        {
            // User-defined pragmas:
            usr = ldoc_nde_new(LDOC_NDE_UA);
            usr->mkup.anno.str = strdup(GEN_ATTRS);
          
            if (nw)
                *nw = true;
        }
        else
        {
            usr = usr_->info.nde;
            
            if (nw)
                *nw = false;
        }
    }
    
    return usr;
}

inline ldoc_nde_t* gen_find_nde(ldoc_nde_t* ctnr1, ldoc_nde_t* ctnr2, char* ky)
{
    ldoc_nde_t* nde = NULL;
    
    TAILQ_FOREACH(nde, &(ctnr1->dscs), ldoc_nde_entries)
    {
        if (!strcmp(nde->mkup.anno.str, ky))
            break;
    }
    
    if (!nde && ctnr2)
        TAILQ_FOREACH(nde, &(ctnr2->dscs), ldoc_nde_entries)
        {
            if (!strcmp(nde->mkup.anno.str, ky))
                break;
        }
    
    return nde;
}

inline void gen_add_nw(ldoc_nde_t* nde, ldoc_nde_t* usr)
{
    if (usr->dsc_cnt || usr->ent_cnt)
        ldoc_nde_dsc_push(nde, usr);
    else
    {
        free(usr->mkup.anno.str);
        ldoc_nde_free(usr);
    }
}

inline bool gen_join_attrs_key(char* id, ldoc_nde_t* nde, ldoc_ent_t* ent, char* attrs)
{
    char* attr_id;
    
    if (id)
        attr_id = id;
    else if (nde)
        attr_id = nde->mkup.anno.str;
    else
        attr_id = ent->pld.pair.anno.str;
    
    // Add attribute separator, if other attributes are present:
    if (qk_working_ptr() != attrs)
        qk_strcat(";");
    
    qk_strcat(attr_id);
    
    if (nde || (ent && ent->tpe != LDOC_ENT_BR))
        qk_strcat("=");
    
    return true;
}

inline bool gen_join_attrs_ent(char* id, ldoc_ent_t* ent, char* attrs)
{
    if (!ent)
        return true;
    
    gen_join_attrs_key(id, NULL, ent, attrs);
    
    char* s;
    bool sp = false;
    switch (ent->tpe)
    {
        case LDOC_ENT_OR:
            s = ent->pld.pair.dtm.str;
            while (*s)
            {
                if (*s == ' ' || *s == '\t')
                {
                    sp = true;
                    break;
                }
                
                s++;
            }
            
            if (sp)
                qk_strcat("\"");
            qk_strcat(ent->pld.pair.dtm.str);
            if (sp)
                qk_strcat("\"");
            break;
        case LDOC_ENT_NR:
            qk_strcat(ent->pld.pair.dtm.str);
            break;
        case LDOC_ENT_BR:
            // Nothing to do: handled by gen_join_attrs_key!
            break;
        default:
            qk_strcat(ent->pld.str);
            break;
    }
    
    return true;
}

inline bool gen_join_nde(ldoc_nde_t* nde)
{
    if (!nde)
        return true;
    
    bool fst = true;
    ldoc_ent_t* ent;
    TAILQ_FOREACH(ent, &(nde->ents), ldoc_ent_entries)
    {
        if (fst)
            fst = false;
        else
            qk_strcat(",");
        
        // TODO Account for more entity types?
        switch (ent->tpe)
        {
            case LDOC_ENT_TXT:
                qk_strcat(ent->pld.str);
                break;
            case LDOC_ENT_OR:
                qk_strcat(ent->pld.pair.anno.str);
                qk_strcat(":");
                qk_strcat(ent->pld.pair.dtm.str);
                break;
            default:
                qk_strcat(ent->pld.pair.dtm.str);
                break;
        }
    }
    
    return true;
}

inline bool gen_join_attrs_nde(char* id, ldoc_nde_t* nde, char* attrs)
{
    if (!nde)
        return true;
    
    gen_join_attrs_key(id, nde, NULL, attrs);
    
    return gen_join_nde(nde);
}

void gen_csepstr_dup(ldoc_nde_t* dst, char* val, bool dup)
{
    bool cnt = true;
    char* val_cmp = val;
    char* val_cmp_;
    ldoc_ent_t* kv_ent;
    do
    {
        // Check: this should always work, since strlen(val) > 0; but check!?
        val++;
        
        if (*val == ',' || !*val)
        {
            if (!*val)
                cnt = false;
            
            *val = 0;
            
            kv_ent = ldoc_ent_new(gen_smrt_flttpe(val_cmp));
            
            if (!kv_ent)
            {
                // TODO Error handling.
            }
            
            if (dup)
            {
                val_cmp_ = strdup(val_cmp);
                
                // TODO Error handling.
            }
            else
                val_cmp_ = val_cmp;
            
            // Remove quotes -- if those exist:
            val_cmp_ = gen_quoskp(val_cmp_);
            
            kv_ent->pld.str = val_cmp_;
            
            ldoc_nde_ent_push(dst, kv_ent);
            
            val_cmp = val + 1;
        }
    } while (cnt);
}

ldoc_nde_t* gen_csep(ldoc_nde_t* dst, gen_attr_t kwd, char* ky, char* val)
{
    return gen_csep_dup(dst, kwd, ky, val, false);
}

ldoc_nde_t* gen_csep_dup(ldoc_nde_t* dst, gen_attr_t kwd, char* ky, char* val, bool dup)
{
    ldoc_nde_t* kv_nde = ldoc_nde_new(LDOC_NDE_OL);
    
    if (!kv_nde)
    {
        // TODO Error handling.
    }
    
    if (kwd.alt)
        kv_nde->mkup.anno.str = (char*)kwd.alt;
    else
    {
        char* ky_;
        
        if (dup)
        {
            ky_ = strdup(ky);
            
            // TODO Error handling.
        }
        else
            ky_ = ky;
        
        kv_nde->mkup.anno.str = ky_;
    }
    
    // TODO Bug? `splt` is never set below!? -- Might be okay, since all tests are fine.
    bool cnt = true;
    char splt = 0;
    char* val_cmp = val;
    char* val_splt;
    char* val_cmp_;
    char* val_splt_;
    ldoc_ent_t* kv_ent;
    do
    {
        // Check: this should always work, since strlen(val) > 0; but check!?
        val++;
        
        if (splt && *val == splt)
        {
            // Key found in key/value pair (separated by splt):
            *(val++) = 0;
            val_splt = val;
        }
        else if (*val == ',' || !*val)
        {
            if (!*val)
                cnt = false;
            
            *val = 0;
            
            // Was: splt ? LDOC_ENT_OR : LDOC_ENT_TXT
            kv_ent = ldoc_ent_new(gen_smrt_flttpe(val_cmp));
            
            if (!kv_ent)
            {
                // TODO Error handling.
            }
            
            if (dup)
            {
                val_cmp_ = strdup(val_cmp);
                
                // TODO Error handling.
            }
            else
                val_cmp_ = val_cmp;
            
            if (splt)
            {
                if (dup)
                {
                    val_splt_ = strdup(val_splt);
                    
                    // TODO Error handling.
                }
                else
                    val_splt_ = val_splt;
                
                kv_ent->pld.pair.anno.str = val_cmp_;
                kv_ent->pld.pair.dtm.str = val_splt_;
            }
            else
            {
                // Remove quotes -- if those exist:
                val_cmp_ = gen_quoskp(val_cmp_);
                
                kv_ent->pld.str = val_cmp_;
            }
            
            ldoc_nde_ent_push(kv_nde, kv_ent);
            
            val_cmp = val + 1;
        }
    } while (cnt);
    
    ldoc_nde_dsc_push(dst, kv_nde);
    
    return kv_nde;
}

void gen_splt_attrs(ldoc_nde_t* ftr, ldoc_nde_t* usr, ldoc_nde_t* ref, ldoc_nde_t* vars, char* attrs, bi_attr upfail)
{
    /* Ruby prototype:
     attributes = {}
     hashes = attribute_string.split(';').map { |assignment|
     match = assignment.match(/([^=]+)=(.+)/) ;
     { match[1].strip => match[2].split(',').map { |value| URI.decode(value.strip) } }
     }
     hashes.map { |hash|
     hash.each_pair { |tag,list|
     attributes[tag] = list
     }
     }
     attributes
     */
    bool brk;
    char* val;
    char* attr = attrs;
    do
    {
        if (!*attrs || *attrs == ';')
        {
            brk = false;
            
            // Note: similar to gvf_tag!
            // If this is the end of the string, then make sure to bail out next:
            if (!*attrs)
                brk = true;
            
            // Isolate attribute assignment by making it an isolated string:
            *attrs = 0;
            
            // Isolate key/value:
            gen_ky(attr, &val);
            
            // Nothing to handle -- break condition 1:
            if (!*attr)
                return;
            
            // Key/value assignment, or, key-only handling.
            // Note: known keys become lower case, user-defined
            //       keys have their case preserved.
            ldoc_nde_t* dst;
            gen_attr_t kwd = { BI_NKW, (char*)NULL };
            gen_kwd(attr, &kwd, upfail);
            if (kwd.attr)
            {
                dst = ftr;
                
                gen_lwrhyph(attr);
                
                // Fall back to comma separated list handing, if
                // this is a comma separated variant list but no
                // variant node is given:
                if (kwd.attr == BI_CSEPVAR && !vars)
                    kwd.attr = BI_CSEP;
            }
            else
                dst = usr;
            
            bool lend = false;
            char* val_cmp;
            char splt = 0;
            // char* val_splt;
            ldoc_nde_t* kv_nde;
            ldoc_ent_t* kv_ent;
            size_t skp = 0;
            switch (kwd.attr)
            {
                case BI_CSEPCPAIR:
                    splt = ':';
                case BI_CSEP:
                    gen_csep(dst, kwd, attr, val);
                    
                    break;
                case BI_CSEPVAR8:
                    skp = 8;
                case BI_CSEPVAR:
                case BI_CSEPVARN:
                    val_cmp = val;
                    TAILQ_FOREACH(kv_nde, &(vars->dscs), ldoc_nde_entries)
                    {
                        if (lend)
                        {
                            // TODO Data format error.
                        }
                        
                        // Check: this should always work, since strlen(val) > 0; but check!?
                        while (*val && *val != ',')
                            val++;
                        
                        if (!*val)
                            lend = true;
                        else
                            *val = 0;
                        
                        kv_ent = ldoc_ent_new(kwd.attr == BI_CSEPVARN ? LDOC_ENT_NR : LDOC_ENT_OR);
                        
                        if (!kv_ent)
                        {
                            // TODO Error handling.
                        }
                        
                        if (kwd.alt)
                            kv_ent->pld.pair.anno.str = (char*)kwd.alt;
                        else
                            kv_ent->pld.pair.anno.str = attr + skp;
                        kv_ent->pld.pair.dtm.str = val_cmp;
                        ldoc_nde_ent_push(kv_nde, kv_ent);
                        
                        val_cmp = ++val;
                    }
                    
                    break;
                case BI_GVFEFFECT:
                    val_cmp = val;
                    // TODO Make this a function in gvf.{c,h}!
                    
                    gvf_proc_effct(vars, &val_cmp, &lend);
                    break;
                case BI_IGN:
                    // Ignore. Key/value pair has been processed already.
                    break;
                case BI_REFSEQ10:
                    skp = 10;
                case BI_REFSEQ:
                case BI_REFSEQN:
                    kv_ent = ldoc_ent_new(kwd.attr == BI_REFSEQN ? LDOC_ENT_NR : LDOC_ENT_OR);
                    
                    if (!kv_ent)
                    {
                        // TODO Error handling.
                    }
                    
                    if (kwd.alt)
                        kv_ent->pld.pair.anno.str = (char*)kwd.alt;
                    else
                        kv_ent->pld.pair.anno.str = &attr[skp];
                    kv_ent->pld.pair.dtm.str = val;
                    ldoc_nde_ent_push(ref, kv_ent);
                    
                    break;
                case BI_TARGET:
                case BI_XCIG:
                    {
                        // Check whether there is an alignment section already:
                        const char* almnt_pth[] = { GEN_ALIGNMENT };
                        ldoc_res_t* almnt = ldoc_find_anno_nde(ftr, (char**)almnt_pth, 1);
                        
                        ldoc_nde_t* almnt_nde;
                        if (almnt)
                        {
                            almnt_nde = almnt->info.nde;
                            ldoc_res_free(almnt);
                        }
                        else
                        {
                            almnt_nde = ldoc_nde_new(LDOC_NDE_UA);
                            
                            // TODO Error handling.
                            
                            almnt_nde->mkup.anno.str = (char*)GEN_ALIGNMENT;
                            
                            ldoc_nde_dsc_push(ftr, almnt_nde);
                        }
                        
                        if (kwd.attr == BI_XCIG)
                        {
                            gen_xcig(val);
                            
                            kv_ent = ldoc_ent_new(LDOC_ENT_OR);
                            kv_ent->pld.pair.anno.str = (char*)kwd.alt;
                            kv_ent->pld.pair.dtm.str = val;
                            
                            ldoc_nde_ent_push(almnt_nde, kv_ent);
                        }
                        else
                        {
                            gff_proc_tgt(almnt_nde, val);
                        }
                    }
                    break;
                default:
                    // Entity type:
                    //   LDOC_ENT_NR  -- if attribute is a number
                    //   LDOC_ENT_OR  -- if attribute is a string
                    //   gen_smrt_tpe -- otherwise
                    kv_ent = ldoc_ent_new(kwd.attr == BI_NUM ? LDOC_ENT_NR : (kwd.attr == BI_VAL ? LDOC_ENT_OR : (kwd.attr == BI_BL ? LDOC_ENT_BR : gen_smrt_tpe(val))));
                    
                    if (!kv_ent)
                    {
                        // TODO Error handling.
                    }
                    
                    // Assign alternative label, if given:
                    if (kwd.alt)
                        kv_ent->pld.pair.anno.str = (char*)kwd.alt;
                    else
                        kv_ent->pld.pair.anno.str = attr;
                    
                    // Assign value, or, truth value if only a keyword was seen:
                    if (kv_ent->tpe == LDOC_ENT_OR)
                    {
                        kv_ent->pld.pair.dtm.str = gen_quoskp(val);
                    }
                    else if (kv_ent->tpe != LDOC_ENT_BR)
                        kv_ent->pld.pair.dtm.str = val;
                    else if (!val)
                        kv_ent->pld.pair.dtm.bl = true;
                    else
                    {
                        gen_lwr(val);
                        
                        if (!strcmp(val, "true"))
                            kv_ent->pld.pair.dtm.bl = true;
                        else
                            kv_ent->pld.pair.dtm.bl = false;
                    }
                    
                    ldoc_nde_ent_push(dst, kv_ent);
                    
                    break;
            }
            
            // Attribute EOS reached earlier -- break condition 2:
            if (brk)
                return;
            
            // Next attribute will have to start here:
            attr = attrs + 1;
        }
        
        attrs++;
    } while (true); // See "break condition 1" and "break condition 2".
}

ldoc_nde_t* gen_variants(char* seq, char sep, char** vseqs, size_t* vnum)
{
    ldoc_nde_t* vars = ldoc_nde_new(LDOC_NDE_UA);
    
    // TODO Error handling.
    
    vars->mkup.anno.str = (char*)GEN_VARIANTS;
    
    if (!seq)
        return vars;
    
    char* v;
    bool cnt = *seq != 0;
    while (cnt)
    {
        v = seq;
        
        cnt = false;
        while (*seq)
        {
            if (*seq == sep)
            {
                *(seq++) = 0;
                cnt = true;
                break;
            }
            else if (*seq == ';')
            {
                // VCF attribute separator.
                *(seq++) = 0;
                break;
            }
            
            seq++;
        }
        
        ldoc_nde_t* var = ldoc_nde_new(LDOC_NDE_UA);
        
        // TODO Error handling.
        
        ldoc_ent_t* seq_i = ldoc_ent_new(LDOC_ENT_OR);

        // TODO Error handling.
        
        seq_i->pld.pair.anno.str = (char*)GEN_SEQUENCE;
        seq_i->pld.pair.dtm.str = v;
        
        vseqs[(*vnum)++] = v;
        
        var->mkup.anno.str = &GEN_ALLELE[2 * *vnum];
        
        ldoc_nde_ent_push(var, seq_i);
        
        ldoc_nde_dsc_push(vars, var);
    }
    
    return vars;
}

static inline uint8_t gen_escchr(char* cptr, gen_filetype_t tpe)
{
    if (*cptr == '\n' ||
        *cptr == '"' ||
        *cptr == '\r' ||
        *cptr == '\t' ||
        *cptr == '\b' ||
        *cptr == '\f')
        return 1;
    else if ((tpe == GEN_FMT_GFF3 || tpe == GEN_FMT_GVF) &&
             *cptr == '%')
    {
        cptr++;
        if ((*cptr >= '0' && *cptr <= '9') ||
            (*cptr >= 'a' && *cptr <= 'z') ||
            (*cptr >= 'A' && *cptr <= 'Z'))
        {
            cptr++;
            if ((*cptr >= '0' && *cptr <= '9') ||
                (*cptr >= 'a' && *cptr <= 'z') ||
                (*cptr >= 'A' && *cptr <= 'Z'))
                return 5;
        }
    }
    
    return 0;
}

char* gen_escstr(char* str, gen_filetype_t tpe)
{
    size_t escchr = 0;
    
    // Count number of characters that need escaping:
    char* ptr = str;
    while (*ptr)
        escchr += gen_escchr(ptr++, tpe);
    
    // Walk backwards and remove trailing newlines/carriage returns:
    while (ptr-- > str)
        if (*ptr == '\n' ||
            *ptr == '\r')
            *ptr = 0;
        else
            break;
    
    size_t nlen = strlen(str) + escchr + 1;
    str = realloc(str, nlen);
    
    // TODO Error handling.
    
    // Escape characters:
    char tmp[nlen];
    char c;
    ptr = str;
    char* tptr = tmp;
    while ((c = *ptr))
        if (gen_escchr(ptr++, tpe))
        {
            *(tptr++) = '\\';
            
            switch (c)
            {
                case '"':
                    *(tptr++) = '"';
                    break;
                case '\n':
                    *(tptr++) = 'n';
                    break;
                case '\r':
                    *(tptr++) = 'r';
                    break;
                case '\t':
                    *(tptr++) = 't';
                    break;
                case '\b':
                    *(tptr++) = 'b';
                    break;
                case '\f':
                    *(tptr++) = 'f';
                    break;
                case '%':
                    if (*ptr == '2')
                    {
                        if (*(ptr + 1) == '6')
                        {
                            *(tptr++) = '&';
                            ptr += 2;
                        }
                        else if (*(ptr + 1) == 'c' || *(ptr + 1) == 'C')
                        {
                            *(tptr++) = ',';
                            ptr += 2;
                        }
                    }
                    else if (*ptr == '3')
                    {
                        if (*(ptr + 1) == 'b' || *(ptr + 1) == 'B')
                        {
                            *(tptr++) = ';';
                            ptr += 2;
                        }
                        else if (*(ptr + 1) == 'd' || *(ptr + 1) == 'D')
                        {
                            *(tptr++) = '=';
                            ptr += 2;
                        }
                    }
                    else
                    {
                        *(tptr++) = 'u';
                        *(tptr++) = '0';
                        *(tptr++) = '0';
                        *(tptr++) = *(ptr++);
                        *(tptr++) = *(ptr++);
                    }
                    break;
                default:
                    // TODO Error handling. Internal error.
                    break;
            }
        }
        else
            *(tptr++) = c;
    *tptr = 0;
    
    memcpy(str, tmp, nlen);
    
    return str;
}

char* gen_rd_ln(fio_mem* mem, off_t mx, size_t llen, char* ln, size_t* ln_len, off_t off)
{
    if (!ln)
    {
        *ln_len = llen;
        ln = (char*)malloc(*ln_len + 1);
    }

    if (!ln)
    {
        // TODO Error handling. Initial alloc failed.
        exit(98);
    }
    
    if (llen > *ln_len)
    {
        ln = realloc(ln, llen + 1);
        *ln_len = llen;
    }
    
    if (!ln)
    {
        // TODO Error handling. Realloc failed.
        exit(99);
    }
    
    // Reached end-of-file:
    if (mem->mx < off + *ln_len)
        *ln_len = mem->mx - off;
    
    // Copy line (including newline characters):
    memcpy(ln, &((char*)mem->pg)[off - mem->off], llen);
    
    // Terminate string:
    ln[llen] = 0;
    
    return ln;
}

void gen_rd(int fd, off_t mx, ldoc_trie_t* idx, gen_cbcks_t* cbcks, gen_ctxt_t* ctxt)
{
    time_t tm_s = time(NULL);
    
    ldoc_doc_t* fdoc = ldoc_doc_new();
    ldoc_doc_t* ldoc;
    
    if (!fdoc)
    {
        // TODO Error handling.
    }
    
    off_t off = 0;
    off_t skp = 0;
    fio_mem* mem = NULL;
    int incr = 0;
    
    gen_prsr_t st;
    st.fa_sct = false;
    st.vcf_ftr_sct = false;
    st.vcf_col = 0;
    
    bool fdoc_purged = false;
    char* cmt = NULL;
    char* mem_cpy = NULL;
    size_t mem_len = 0;
    size_t lnlen;
    uint64_t ln_no = 0;
    gen_fstat stat = { 0, 0, 0, false, 0, 0 };
    while (!mem || mem->off + mem->ln < mx)
    {
        // Calculate offsets in case `off` is not landing on a page size:
        if (mem)
        {
            skp = off;
            
            off /= getpagesize();
            off *= getpagesize();
            
            skp -= off;
        }

        // Goto for increasing number of pages; note that `incr` is not set back, but kept on a high watermark:
    rd_incr_mem:
        mem = fio_mmap(NULL, fd, mx, getpagesize() * (BI_GEN_PG_MUL + incr), off);
        
        // TODO Error checking.
        
        // Figure out if (at least) one line can be read (based on current pointer):
        lnlen = fio_lnlen(mem, mem->off + skp);
        
        // Handle case where no line ending is visible:
        if (!lnlen && mem->off + mem->ln < mx)
        {
            incr++;
            goto rd_incr_mem;
        }
        
        // Still no line ending visible? Then read to the end of the buffer (this is implicit; check fio_mmap behavior):
        if (!lnlen)
            lnlen = mem->mx - mem->off;
        
        // Adjust offset in case of re-mapping on a non-page boundary:
        if (skp)
            off += skp;
        
        do
        {
            // Current line:
            mem_cpy = gen_rd_ln(mem, mx, lnlen, mem_cpy, &mem_len, off);

            ldoc = cbcks->proc_ln(fd, mx, fdoc, idx, mem_cpy, lnlen, &st, &cmt, &stat);
            
            if (ldoc)
            {
                if (!fdoc_purged)
                {
                    // Meta information about the file and its data contents:
                    
                    ldoc_doc_t* cdoc = ldoc_doc_new();

                    ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_OR);
                    ent->pld.pair.anno.str = strdup(JSONLD_CTX);
                    switch (ctxt->tpe)
                    {
                        case GEN_FMT_GFF3:
                            ent->pld.pair.dtm.str = strdup(JSONLD_GFF3_CTX1);
                            break;
                        case GEN_FMT_GTF:
                            ent->pld.pair.dtm.str = strdup(JSONLD_GTF_CTX1);
                            break;
                        case GEN_FMT_GVF:
                            ent->pld.pair.dtm.str = strdup(JSONLD_GVF_CTX1);
                            break;
                        case GEN_FMT_VCF:
                            ent->pld.pair.dtm.str = strdup(JSONLD_VCF_CTX1);
                            break;
                        default:
                            // TODO Internal error.
                            break;
                    }
                    ldoc_nde_ent_push(cdoc->rt, ent);

                    ldoc_ent_t* tpe = ldoc_ent_new(LDOC_ENT_OR);
                    tpe->pld.pair.anno.str = strdup(JSONLD_TPE);
                    tpe->pld.pair.dtm.str = strdup(JSONLD_CLSS_CTX);
                    ldoc_nde_ent_push(cdoc->rt, tpe);
                    
                    ent = ldoc_ent_new(LDOC_ENT_OR);
                    ent->pld.pair.anno.str = strdup("biointerchange-version");
                    ent->pld.pair.dtm.str = strdup(ctxt->ver);
                    ldoc_nde_ent_push(cdoc->rt, ent);
                    
                    ent = ldoc_ent_new(LDOC_ENT_OR);
                    ent->pld.pair.anno.str = strdup("input-file");
                    ent->pld.pair.dtm.str = strdup(ctxt->fgen);
                    ldoc_nde_ent_push(cdoc->rt, ent);
                    
                    ent = ldoc_ent_new(LDOC_ENT_OR);
                    ent->pld.pair.anno.str = strdup("output-file");
                    if (ctxt->fname)
                        ent->pld.pair.dtm.str = strdup(ctxt->fname);
                    else
                        ent->pld.pair.dtm.str = NULL;
                    ldoc_nde_ent_push(cdoc->rt, ent);

                    ent = ldoc_ent_new(LDOC_ENT_OR);
                    ent->pld.pair.anno.str = strdup("input-filetype");
                    switch (ctxt->tpe)
                    {
                        case GEN_FMT_GFF3:
                            ent->pld.pair.dtm.str = strdup("GFF3");
                            break;
                        case GEN_FMT_GTF:
                            ent->pld.pair.dtm.str = strdup("GTF");
                            break;
                        case GEN_FMT_GVF:
                            ent->pld.pair.dtm.str = strdup("GVF");
                            break;
                        case GEN_FMT_VCF:
                            ent->pld.pair.dtm.str = strdup("VCF");
                            break;
                        default:
                            // TODO Internal error.
                            break;
                    }
                    ldoc_nde_ent_push(cdoc->rt, ent);
                    
                    ent = ldoc_ent_new(LDOC_ENT_OR);
                    ent->pld.pair.anno.str = strdup("user-parameter");
                    if (ctxt->usr)
                        ent->pld.pair.dtm.str = strdup(ctxt->usr);
                    else
                        ent->pld.pair.dtm.str = NULL;
                    ldoc_nde_ent_push(cdoc->rt, ent);

                    ent = ldoc_ent_new(LDOC_ENT_OR);
                    ent->pld.pair.anno.str = strdup("python-callback");
                    if (ctxt->py)
                        ent->pld.pair.dtm.str = strdup(ctxt->pycall);
                    else
                        ent->pld.pair.dtm.str = NULL;
                    ldoc_nde_ent_push(cdoc->rt, ent);
                    
                    gen_ser(ctxt, GEN_CTPE_SETUP, cdoc, fdoc, &stat);
                    
                    ldoc_doc_free(cdoc);
                    
                    fdoc_purged = true;
                }
                
                gen_ser(ctxt, GEN_CTPE_PROCESS, ldoc, NULL, &stat);
                
                ldoc_doc_free(ldoc);
            }

            // Reset quick memory:
            qk_purge();
            
            // Next line:
            off += lnlen;
            ln_no++;
            lnlen = fio_lnlen(mem, off);
        } while (lnlen);
    }
    
    if (mem_cpy)
        free(mem_cpy);
    
    if (mem)
        fio_munmap(mem);
    
    // Wrap up:
    time_t tm_e = time(NULL);
    char* tm_sstr = strdup(ctime(&tm_s));
    char* tm_estr = strdup(ctime(&tm_e));
    
    ldoc_doc_t* sdoc = ldoc_doc_new();
    
    ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_OR);
    ent->pld.pair.anno.str = strdup(JSONLD_CTX);
    ent->pld.pair.dtm.str = strdup(JSONLD_STAT_1);
    ldoc_nde_ent_push(sdoc->rt, ent);
    
    ldoc_ent_t* stat_tpe = ldoc_ent_new(LDOC_ENT_OR);
    stat_tpe->pld.pair.anno.str = strdup(JSONLD_TPE);
    stat_tpe->pld.pair.dtm.str = strdup(JSONLD_CLSS_SUM);
    ldoc_nde_ent_push(sdoc->rt, stat_tpe);
    
    ldoc_nde_t* fstat = ldoc_nde_new(LDOC_NDE_UA);
    fstat->mkup.anno.str = strdup(GEN_STATS);
    ldoc_nde_dsc_push(sdoc->rt, fstat);

    char num[21]; // Space for a 64-bit number plus null-byte.
    
    sprintf(num, "%u", stat.comms);
    ldoc_ent_t* ent_stat = ldoc_ent_new(LDOC_ENT_NR);
    ent_stat->pld.pair.anno.str = strdup(GEN_STAT_CMMS);
    ent_stat->pld.pair.dtm.str = strdup(num);
    ldoc_nde_ent_push(fstat, ent_stat);

    sprintf(num, "%u", stat.meta);
    ent_stat = ldoc_ent_new(LDOC_ENT_NR);
    ent_stat->pld.pair.anno.str = strdup(GEN_STAT_META);
    ent_stat->pld.pair.dtm.str = strdup(num);
    ldoc_nde_ent_push(fstat, ent_stat);

    ent_stat = ldoc_ent_new(LDOC_ENT_BR);
    ent_stat->pld.pair.anno.str = strdup(GEN_STAT_META_FLTR);
    ent_stat->pld.pair.dtm.bl = stat.cbk_meta_fltr;
    ldoc_nde_ent_push(fstat, ent_stat);
    
    sprintf(num, "%u", stat.ftrs);
    ent_stat = ldoc_ent_new(LDOC_ENT_NR);
    ent_stat->pld.pair.anno.str = strdup(GEN_STAT_FTRS);
    ent_stat->pld.pair.dtm.str = strdup(num);
    ldoc_nde_ent_push(fstat, ent_stat);

    sprintf(num, "%u", stat.cbk_ftrs_fltr);
    ent_stat = ldoc_ent_new(LDOC_ENT_NR);
    ent_stat->pld.pair.anno.str = strdup(GEN_STAT_FTRS_FLTR);
    ent_stat->pld.pair.dtm.str = strdup(num);
    ldoc_nde_ent_push(fstat, ent_stat);
    
    ldoc_nde_t* rntm = ldoc_nde_new(LDOC_NDE_UA);
    rntm->mkup.anno.str = strdup(GEN_STAT_RNTM);
    ldoc_nde_dsc_push(sdoc->rt, rntm);
    
    ldoc_ent_t* ent_tm_s = ldoc_ent_new(LDOC_ENT_OR);
    ent_tm_s->pld.pair.anno.str = strdup(GEN_STAT_INVC);
    ent_tm_s->pld.pair.dtm.str = gen_term_crnl(tm_sstr);
    ldoc_nde_ent_push(rntm, ent_tm_s);

    ldoc_ent_t* ent_tm_e = ldoc_ent_new(LDOC_ENT_OR);
    ent_tm_e->pld.pair.anno.str = strdup(GEN_STAT_FNSH);
    ent_tm_e->pld.pair.dtm.str = gen_term_crnl(tm_estr);
    ldoc_nde_ent_push(rntm, ent_tm_e);

    sprintf(num, "%lu", tm_e - tm_s);
    ldoc_ent_t* ent_tm = ldoc_ent_new(LDOC_ENT_OR);
    ent_tm->pld.pair.anno.str = strdup(GEN_STAT_LSEC);
    ent_tm->pld.pair.dtm.str = strdup(num);
    ldoc_nde_ent_push(rntm, ent_tm);
    
    gen_ser(ctxt, GEN_CTPE_CLEANUP, sdoc, NULL, &stat);
    
    ldoc_doc_free(sdoc);
}

inline void gen_ser(gen_ctxt_t* ctxt, gen_ctpe_t ctpe, ldoc_doc_t* doc, ldoc_doc_t* opt, gen_fstat* stat)
{
    if (!doc)
        return;
    
    ldoc_ser_t* ser;

    if (ctxt->py)
    {
        ldoc_doc_t* ldoc_anno = NULL;
        
        switch (ctpe)
        {
            case GEN_CTPE_SETUP:
                // doc is the "context", which cannot be modifed:
                ser = ldoc_format(doc, json_vis_nde, json_vis_ent);
                fprintf(ctxt->fout, "%s\n", ser->pld.str);
                ldoc_ser_free(ser);
                
                // ldoc_anno will be a possibly modified version of opt:
                ldoc_anno = py_setup(doc, opt);
                ser = ldoc_format(doc, json_vis_nde, json_vis_ent);
                fprintf(ctxt->fout, "%s\n", ser->pld.str);
                ldoc_ser_free(ser);
                
                if (!ldoc_anno)
                    stat->cbk_meta_fltr = true;
                
                break;
            case GEN_CTPE_CLEANUP:
                ldoc_anno = py_cleanup(doc);
                break;
            case GEN_CTPE_PROCESS:
                ldoc_anno = py_process(doc);
                
                if (!ldoc_anno)
                    stat->cbk_ftrs_fltr++;
                break;
            default:
                // TODO Internal error.
                break;
        }
        
        if (ldoc_anno)
        {
            ser = ldoc_format(ldoc_anno, json_vis_nde, json_vis_ent);
            fprintf(ctxt->fout, "%s\n", ser->pld.str);
            ldoc_ser_free(ser);
        }
    }
    else
    {
        ser = ldoc_format(doc, json_vis_nde, json_vis_ent);
        fprintf(ctxt->fout, "%s\n", ser->pld.str);
        ldoc_ser_free(ser);
        
        if (opt)
        {
            ser = ldoc_format(opt, json_vis_nde, json_vis_ent);
            fprintf(ctxt->fout, "%s\n", ser->pld.str);
            ldoc_ser_free(ser);
        }
    }
}

/// Genomic file format serialization

char* gen_exp_ky(char* ky)
{
    char* str = ky;
    
    if (*str == 's')
    {
        str++;
        if (*str == 'e')
        {
            str++;
            if (*str == 'q')
            {
                str++;
                if (!strcmp(str, "uence"))
                {
                    ky = (char*)GEN_SEQUENCE_GVF;
                }
            }
        }
    }

    return ky;
}

inline void gen_proc_doc_prgm_kv(ldoc_nde_t* cntnr, char* ky, char* alt, char* sep)
{
    ldoc_res_t* ent = ldoc_find_anno_ent(cntnr, ky);
    
    if (ent && ent->info.ent->pld.pair.dtm.str)
    {
        if (!qk_heap_empty())
            qk_strcat("\n");
        
        qk_strcat("##");
        qk_strcat(alt);
        qk_strcat(sep);
        qk_strcat(ent->info.ent->pld.pair.dtm.str);
    }
}

bool gen_proc_doc_prgm(ldoc_nde_t* prgm, char* sep)
{
    const char* usr_id[] = { GEN_ATTRS };
    ldoc_res_t* usr = ldoc_find_anno_nde(prgm, (char**)usr_id, 1);
    
    // TODO Error handling.
    
    if (!usr)
        return true;
    
    ldoc_nde_t* dsc;
    TAILQ_FOREACH(dsc, &(usr->info.nde->dscs), ldoc_nde_entries)
    {
        ldoc_ent_t* ent;
        TAILQ_FOREACH(ent, &(dsc->ents), ldoc_ent_entries)
        {
            if (!qk_heap_empty())
                qk_strcat("\n");

            qk_strcat("##");
            qk_strcat(dsc->mkup.anno.str);
            qk_strcat(sep);
            qk_strcat(ent->pld.str);
        }
    }
    
    return true;
}

// TODO Obsolete?
bool gen_proc_doc_usr(ldoc_nde_t* ftr)
{
    const char* usr_id[] = { GEN_ATTRS };
    ldoc_res_t* usr = ldoc_find_anno_nde(ftr, (char**)usr_id, 1);
    
    // TODO Error handling.
    
    if (!usr)
        return true;
    
    // TODO Assumes that other attributes are present: so the
    //      semi-colon is always prepended. This might be a
    //      false assumption! Unlikely to happen -- but could happen!
    ldoc_ent_t* ent;
    TAILQ_FOREACH(ent, &(usr->info.nde->ents), ldoc_ent_entries)
    {
        qk_strcat(";");
        
        qk_strcat(ent->pld.pair.anno.str);
        qk_strcat("=");
        qk_strcat(ent->pld.pair.dtm.str);
    }
    
    return true;
}

// TODO Replace vstr with parametrizable quick-heap implementation.
bool gen_proc_nde(ldoc_nde_t* vars, char* attr, char* pre, char* astr, size_t vnum)
{
    // Attribute name might need adjusting:
    // (For example: "sequence" is shortened to "seq" in GVF
    char* attr_adj = gen_exp_ky(attr);
    
    // Attribute key:
    if (pre)
        strcat(astr, pre);
    strcat(astr, attr_adj);
    strcat(astr, "=");
    
    bool fst = true;
    ldoc_res_t* ent;
    ldoc_res_t* var;
    char* all_pth[1];
    size_t vrnt = 0;
    for (; vrnt < vnum; vrnt++)
    {
        if (fst)
            fst = false;
        else
            strcat(astr, ",");
        
        all_pth[0] = &GEN_ALLELE[(vrnt + 1) * 2];
        
        // TODO The order of the nodes will always be
        //      the same for following runs. This means
        //      that there is room for optimization.
        var = ldoc_find_anno_nde(vars, all_pth, 1);
        
        // TODO Error handling. Data error -- not supported.
        
        ent = ldoc_find_anno_ent(var->info.nde, attr);
        
        // TODO Error handling. Data error -- not supported.
        
        strcat(astr, ent->info.ent->pld.pair.dtm.str);
        
        // Entity has been processed, so remove it. Prevents that the
        // same information is serialized at a later stage again.
        ldoc_ent_rm(ent->info.ent);
        ldoc_ent_free(ent->info.ent);
    }

    return true;
}

/// JSON to genomics formats

void gen_rd_doc(int fd, off_t mx, gen_ctxt_t* ctxt)
{
    off_t err;
    off_t nxt;
    
    // NOTE: A lot of this code is derived from gen_rd. It might
    //       be possible to write a function for this.
    
    ldoc_doc_t* ldoc;
    
    off_t off = 0;
    off_t skp = 0;
    fio_mem* mem = NULL;
    int incr = 0;
    
    gen_prsr_t st;
    st.fa_sct = false;
    st.vcf_ftr_sct = false;
    st.vcf_col = 0;
    
    bool fdoc_purged = false;
    char* mem_cpy = NULL;
    size_t mem_len = 0;
    size_t lnlen;
    uint64_t ln_no = 0;
    gen_fstat stat = { 0, 0, 0, false, 0, 0 };
    // Temporary -- will be replaced with the actual output once '@context'
    // has been read:
    gen_filetype_t out = GEN_FMT_LDJ;
    gen_doctype_t tpe;
    while (!mem || mem->off + mem->ln < mx)
    {
        // Calculate offsets in case `off` is not landing on a page size:
        if (mem)
        {
            skp = off;
            
            off /= getpagesize();
            off *= getpagesize();
            
            skp -= off;
        }
        
        // Goto for increasing number of pages; note that `incr` is not set back, but kept on a high watermark:
    doc_incr_mem:
        mem = fio_mmap(NULL, fd, mx, getpagesize() * (BI_GEN_PG_MUL + incr), off);
        
        // TODO Error checking.
        
        // Figure out if (at least) one line can be read (based on current pointer):
        lnlen = fio_lnlen(mem, mem->off + skp);
        
        // Handle case where no line ending is visible:
        if (!lnlen && mem->off + mem->ln < mx)
        {
            incr++;
            goto doc_incr_mem;
        }
        
        // Still no line ending visible? Then read to the end of the buffer (this is implicit; check fio_mmap behavior):
        if (!lnlen)
            lnlen = mem->mx - mem->off;
        
        // Adjust offset in case of re-mapping on a non-page boundary:
        if (skp)
            off += skp;
        
        do
        {
            // Current line:
            mem_cpy = gen_rd_ln(mem, mx, lnlen, mem_cpy, &mem_len, off);
            
            // ldoc = ldoc_ldjson_read(mem_cpy, mem_len, &err, NULL);
            ldoc = ldoc_json_read(mem_cpy, mem_len, &err);

            // Establish document context:
            
            ldoc_res_t* json_ctx = ldoc_find_anno_ent(ldoc->rt, (char*)JSONLD_CTX);
            
            if (!json_ctx)
            {
                // TODO Error handling.
                exit(100);
            }
            
            if (json_ctx->nde)
            {
                // TODO Error handling.
                exit(100);
            }
            
            if (json_ctx->info.ent->tpe != LDOC_ENT_OR)
            {
                // TODO Error handling.
                exit(101);
            }
            
            if (out == GEN_FMT_LDJ)
            {
                if (!strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_GFF3_CTX1))
                    out = GEN_FMT_GFF3;
                else if (!strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_GTF_CTX1))
                    out = GEN_FMT_GTF;
                else if (!strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_GVF_CTX1))
                    out = GEN_FMT_GVF;
                else if (!strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_VCF_CTX1))
                    out = GEN_FMT_VCF;
                else
                {
                    // TODO Error handling.
                    exit(102);
                }
                
                tpe = GEN_FMT_CTX;
            }
            else if (!strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_GFF3_X1) ||
                     !strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_GTF_X1) ||
                     !strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_GVF_X1) ||
                     !strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_VCF_X1))
                tpe = GEN_FMT_INF;
            else if (!strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_GFF3_1) ||
                     !strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_GTF_1) ||
                     !strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_GVF_1) ||
                     !strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_VCF_1))
                tpe = GEN_FMT_FTR;
            else if (!strcmp(json_ctx->info.ent->pld.pair.dtm.str, JSONLD_STAT_1))
                tpe = GEN_FMT_STT;
    
            // Output data:
            
            char* lns = NULL;
            if (tpe == GEN_FMT_INF || tpe == GEN_FMT_FTR)
            {
                switch (out)
                {
                    case GEN_FMT_GFF3:
                        lns = gff_proc_doc(ldoc, tpe);
                        break;
                    case GEN_FMT_GTF:
                        // TODO
                        exit(1001);
                        break;
                    case GEN_FMT_GVF:
                        lns = gvf_proc_doc(ldoc, tpe);
                        break;
                    case GEN_FMT_VCF:
                        lns = vcf_proc_doc(ldoc, tpe);
                        break;
                    default:
                        // TODO Internal error.
                        exit(100);
                }
            }
            
            if (lns)
                fprintf(ctxt->fout, "%s\n", lns);
            
            // Wrap up:
            
            // ldoc = cbcks->proc_ln(fd, mx, fdoc, idx, mem_cpy, lnlen, &st, &cmt, &stat);
            
            //gen_ser(ctxt, GEN_CTPE_PROCESS, ldoc, NULL, &stat);
            ldoc_doc_free(ldoc);
            //}
            
            // Reset quick memory (resets lns):
            qk_purge();
            
            // Next line:
            off += lnlen;
            ln_no++;
            lnlen = fio_lnlen(mem, off);
        } while (lnlen);
    }
    
    if (mem_cpy)
        free(mem_cpy);
    
    if (mem)
        fio_munmap(mem);
}

/// Quick SINGLE THREADED string operations

static size_t qk_size;
static char* qk_heap;
static char* qk_ptr;

char* qk_alloc(size_t n)
{
    qk_heap = (char*)malloc(n);
    
    qk_ptr = qk_heap;
    qk_size = n;
    
    return qk_heap;
}

void qk_free()
{
    free(qk_heap);
}

inline char* qk_heap_ptr()
{
    return qk_heap;
}

inline char* qk_working_ptr()
{
    return qk_ptr;
}

inline bool qk_heap_empty()
{
    return qk_heap == qk_ptr;
}

inline void qk_purge()
{
    qk_ptr = qk_heap;
}

inline char* qk_strdup(const char* s1)
{
    size_t len = strlen(s1);
    
    // Out of quick memory:
    if (qk_ptr - qk_heap + len + 1 > qk_size)
        return NULL;
    
    off_t cpy = qk_ptr;
    memcpy(qk_ptr, s1, len + 1);
    
    qk_ptr += len + 1;
    
    return cpy;
}

inline char* qk_strndup(const char* s1, size_t n)
{
    // Out of quick memory:
    if (qk_ptr - qk_heap + n + 1 > qk_size)
        return NULL;
    
    off_t cpy = qk_ptr;
    memcpy(qk_ptr, s1, n);
    qk_ptr[n] = 0;
    
    qk_ptr += n + 1;
    
    return cpy;
}

inline bool qk_strcat(const char* s1)
{
    if (!s1)
        return true;
    
    size_t n = strlen(s1);
    
    // Out of quick memory (+1 required in case heap is empty):
    if (qk_ptr - qk_heap + n + 1 > qk_size)
        return false;
    
    // Reposition qk_ptr to null-byte of previous string:
    if (qk_ptr > qk_heap)
        qk_ptr--;

    memcpy(qk_ptr, s1, n + 1);
    
    qk_ptr += n + 1;
    
    return true;
}

inline bool qk_strncat(const char* s1, size_t n)
{
    if (!s1)
        return true;
    
    // Out of quick memory (+1 required in case heap is empty):
    if (qk_ptr - qk_heap + n + 1 > qk_size)
        return false;
    
    // Reposition qk_ptr to null-byte of previous string:
    if (qk_ptr > qk_heap)
        qk_ptr--;
    
    memcpy(qk_ptr, s1, n);
    
    qk_ptr += n + 1;
    
    *qk_ptr = 0;
    
    return true;
}

