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

const char* GEN_ATTRS = "user";
const char* GEN_BUILD = "build";
const char* GEN_COMMENT = "comment";
const char* GEN_END = "end";
const char* GEN_LOCUS = "locus";
const char* GEN_REFERENCE = "reference";
const char* GEN_SEQUENCE = "sequence";
const char* GEN_START = "start";
const char* GEN_SOURCE = "source";

const char* GEN_NULL = "null";
const char* GEN_TRUE = "true";

const char* GEN_COUNT = "count";
const char* GEN_DEPTH = "depth";
const char* GEN_FREQUENCY = "frequency";

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
    
    for (size_t j = 0; j < GEN_MAX_ALT; j++)
    {
        GEN_ALLELE[2 * j] = 'A' + j;
        GEN_ALLELE[2 * j + 1] = 0;
        
        for (size_t k = 0; k < GEN_MAX_ALT; k++)
            if (j <= k)
                gen_allele_lbl(s, j, k);
    }
}

static inline void gen_ky(char* attr, char** val)
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

inline void gen_lwr(char* str)
{
    while (*str)
    {
        if (*str >= 'A' && *str <= 'Z')
            *str = *str - 'A' + 'a';
        
        str++;
    }
}

inline bi_attr gen_kwd(char* str)
{
    if (*str >= 'A' && *str <= 'Z')
    {
        if (*str == 'A')
        {
            str++;
            if (*str == 'C')
            {
                str++;
                if (!*str)
                    return BI_CSEPVAR; // VCF: AC, allele count
            }
            else if (*str == 'F')
            {
                str++;
                if (!*str)
                    return BI_CSEPVAR; // VCF: AF, allele frequency
            }
            else if (*str == 'N')
            {
                str++;
                if (!*str)
                    return BI_NUM; // VCF: AN, total number of alleles in called genotypes
            }
            else if (*str == 'l')
            {
                str++;
                if (*str == 'i')
                {
                    str++;
                    if (!strcmp(str, "as")) // GFF3: Alias
                        return BI_CSEP;
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
                    return BI_NUM; // VCF: BQ, RMS base quality
            }
        }
        else if (*str == 'D')
        {
            str++;
            if (*str == 'P')
            {
                str++;
                if (!*str)
                    return BI_NUM; // VCF: DP, depth across samples
            }
            else if (*str == 'b')
            {
                str++;
                if (*str == 'x')
                {
                    str++;
                    if (!strcmp(str, "ref")) // GFF3: Dbxref
                        return BI_CSEP;
                }
            }
        }
        else if (*str == 'E')
        {
            str++;
            if (*str == 'N')
            {
                str++;
                if (*str == 'D')
                {
                    str++;
                    if (!*str)
                        return BI_NUM; // VCF: END, end position of variant description
                }
            }
        }
        else if (*str == 'G')
        {
            str++;
            if (*str == 'a')
            {
                str++;
                if (*str == 'p')
                {
                    str++;
                    if (!*str)
                        return BI_XCIG; // GFF3: Gap
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
                    return BI_NUM; // VCF: MQ, mapping quality
                else if (*str == '0')
                {
                    str++;
                    if (!*str)
                        return BI_NUM; // VCF: MQ0, mapping quality == 0
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
                    return BI_NUM; // VCF: NS, number of samples with data
            }
            else if (*str == 'o')
            {
                str++;
                if (*str == 't')
                {
                    str++;
                    if (!strcmp(str, "e"))
                        return BI_CSEP; // GFF3: Note
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
                        return BI_CSEP; // GFF3: Ontology_term
                }
            }
        }
        else if (*str == 'P')
        {
            str++;
            if (*str == 'a')
            {
                str++;
                if (*str == 'r')
                {
                    str++;
                    if (!strcmp(str, "ent"))
                        return BI_CSEP; // GFF3: Parent
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
                        return BI_IGN;
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
                        return BI_IGN; // GVF: Variant_seq
                }
            }
        }
        
        return BI_VAL;
    }
    
    return BI_NKW;
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

void gen_splt_attrs(ldoc_nde_t* ftr, ldoc_nde_t* usr, ldoc_nde_t* vars, char* attrs)
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
    char* attr = attrs;
    do
    {
        if (!*attrs || *attrs == ';')
        {
            bool brk = false;
            
            // If this is the end of the string, then make sure to bail out next:
            if (!*attrs)
                brk = true;
            
            // Isolate attribute assignment by making it an isolated string:
            *attrs = 0;
            
            // Isolate key/value:
            char* val;
            gen_ky(attr, &val);
            
            // Nothing to handle -- break condition 1:
            if (!*attr)
                return;
            
            // Key/value assignment, or, key-only handling:
            ldoc_nde_t* dst;
            bi_attr kind;
            // Yes, this is an assignment. Do not '=='!
            if ((kind = gen_kwd(attr)))
            {
                dst = ftr;
                
                gen_lwr(attr);
                
                // Fall back to comma separated list handing, if
                // this is a comma separated variant list but no
                // variant node is given:
                if (kind == BI_CSEPVAR && !vars)
                    kind = BI_CSEP;
            }
            else
            {
                dst = usr;
                
                if (!val)
                    val = (char*)GEN_TRUE;
            }
            
            bool lend = false;
            char* val_cmp;
            ldoc_nde_t* kv_nde;
            ldoc_ent_t* kv_ent;
            switch (kind)
            {
                case BI_CSEP:
                    kv_nde = ldoc_nde_new(LDOC_NDE_OL);
                    
                    if (!kv_nde)
                    {
                        // TODO Error handling.
                    }
                    
                    kv_nde->mkup.anno.str = attr;
                    
                    val_cmp = val;
                    do
                    {
                        // Check: this should always work, since strlen(val) > 0; but check!?
                        val++;
                        
                        if (*val == ',' || !*val)
                        {
                            *val = 0;
                            
                            // TODO Possible bug: LDOC_ENT_TXT but use of pld.pair!
                            
                            kv_ent = ldoc_ent_new(LDOC_ENT_TXT);
                            
                            if (!kv_ent)
                            {
                                // TODO Error handling.
                            }
                            
                            kv_ent->pld.pair.anno.str = attr;
                            kv_ent->pld.pair.dtm.str = val_cmp;
                            ldoc_nde_ent_push(kv_nde, kv_ent);
                            
                            val_cmp = val;
                        }
                    } while (*val);
                    
                    ldoc_nde_dsc_push(dst, kv_nde);
                    
                    break;
                case BI_CSEPVAR:
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
                        
                        kv_ent = ldoc_ent_new(LDOC_ENT_OR);
                        
                        if (!kv_ent)
                        {
                            // TODO Error handling.
                        }
                        
                        kv_ent->pld.pair.anno.str = attr;
                        kv_ent->pld.pair.dtm.str = val_cmp;
                        ldoc_nde_ent_push(kv_nde, kv_ent);
                        
                        val_cmp = ++val;
                    }
                    
                    break;
                case BI_IGN:
                    // Ignore. Key/value pair has been processed already.
                    break;
                default:
                    kv_ent = ldoc_ent_new(kind == BI_NUM ? LDOC_ENT_NR : LDOC_ENT_OR);
                    
                    if (!kv_ent)
                    {
                        // TODO Error handling.
                    }
                    
                    if (kind == BI_XCIG)
                        gen_xcig(val);
                    
                    kv_ent->pld.pair.anno.str = attr;
                    kv_ent->pld.pair.dtm.str = val;
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
    
    vars->mkup.anno.str = "variants";
    
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
    }
    
    if (llen > *ln_len)
    {
        ln = realloc(ln, llen + 1);
        *ln_len = llen;
    }
    
    if (!ln)
    {
        // TODO Error handling. Realloc failed.
    }
    
    // Copy line (including newline characters):
    memcpy(ln, mem->pg + (off - mem->off), *ln_len);
    
    // Terminate string:
    ln[*ln_len] = 0;
    
    return ln;
}

void gen_rd(int fd, off_t mx, ldoc_trie_t* idx, gen_cbcks_t* cbcks)
{
    ldoc_doc_t* fdoc = ldoc_doc_new();
    
    if (!fdoc)
    {
        // TODO Error handling.
    }
    
    off_t off = 0;
    fio_mem* mem = NULL;
    int incr = 0;
    
    gen_prsr_t st;
    st.fa_sct = false;
    st.vcf_ftr_sct = false;
    st.vcf_col = 0;
    
    char* cmt = NULL;
    char* mem_cpy = NULL;
    size_t mem_len = 0;
    size_t lnlen;
    uint64_t ln_no = 0;
    while (!mem || mem->off + mem->ln < mx)
    {
        // Goto for increasing number of pages; note that `incr` is not set back, but kept on a high watermark:
    gff_rd_incr_mem:
        mem = fio_mmap(NULL, fd, mx, getpagesize() * (BI_GEN_PG_MUL + incr), off);
        
        // TODO Error checking.
        
        // Figure out if (at least) one line can be read (based on current pointer):
        lnlen = fio_lnlen(mem, mem->off);
        
        // Handle case where no line ending is visible:
        if (!lnlen && mem->off + mem->ln < mx)
        {
            incr++;
            goto gff_rd_incr_mem;
        }
        
        // Still no line ending visible? Then read to the end of the buffer (this is implicit; check fio_mmap behavior):
        if (!lnlen)
            lnlen = mem->mx - mem->off;
        
        do {
            // Current line:
            mem_cpy = gen_rd_ln(mem, mx, lnlen, mem_cpy, &mem_len, off);
            
            cbcks->proc_ln(fd, mx, fdoc, idx, mem_cpy, lnlen, &st, &cmt);
            
            // Reset quick memory:
            qk_purge();
            
            // Next line:
            off += lnlen;
            ln_no++;
            lnlen = fio_lnlen(mem, off);
            
            //printf("%lu: %lu\n", ln_no, lnlen);
            //printf("%s", mem_cpy);
        } while (lnlen);
    }
    
    if (mem_cpy)
        free(mem_cpy);
    
    if (mem)
        fio_munmap(mem);
    
    // Meta information about the file and its data contents:
    
    ldoc_ser_t* ser = ldoc_format(fdoc, json_vis_nde, json_vis_ent);
    
    // Only output non-empty documents:
    if (fdoc->rt->dsc_cnt > 0 || fdoc->rt->ent_cnt > 0)
        printf("%s\n", ser->sclr.str);
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

void qk_purge()
{
    qk_ptr = qk_heap;
}

char* qk_strdup(const char* s1)
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

char* qk_strndup(const char* s1, size_t n)
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

