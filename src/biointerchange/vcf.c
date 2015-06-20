
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

static const char* VCF_C1   = "landmark";
static const char* VCF_C2_1 = "start";
static const char* VCF_C2_2 = "end";
static const char* VCF_C3   = "id";
static const char* VCF_C4   = "reference";
// static const char* VCF_C5 handled by generic implementation.
static const char* VCF_C6   = "score";
static const char* VCF_C7   = "filters";

static const char* VCF_ALLELES = "alleles";
static const char* VCF_GENOTYPE = "genotype";
static const char* VCF_GENOLHOOD = "genotype-likelihood";
static const char* VCF_GENOLHOODE = "gle"; // TODO Obsolete. (Discussed on VCF mailing list. Will affect VCF 4.3 and later.
static const char* VCF_PHREDLHOOD = "genotype-likelihood-phred-scaled";
static const char* VCF_PHASED = "phased";
static const char* VCF_SAMPLES = "samples";

void vcf_cbcks(gen_cbcks_t* cbcks)
{
    cbcks->proc_ln = &vcf_proc_ln;
}

static inline ldoc_struct_t vcf_prgm_tpe(char* ky)
{
    char* ky_ = ky;

    if (*ky == 'A')
    {
        ky++;
        if (*ky == 'L')
        {
            ky++;
            if (*ky == 'T')
            {
                ky++;
                if (!*ky) // ALT
                    return LDOC_NDE_OO;
            }
        }
    }
    else if (*ky == 'F')
    {
        ky++;
        if (*ky == 'I')
        {
            ky++;
            if (*ky == 'L')
            {
                ky++;
                if (!strcmp(ky, "TER")) // FILTER
                    return LDOC_NDE_OO;
            }
        }
        else if (*ky == 'O')
        {
            ky++;
            if (*ky == 'R')
            {
                ky++;
                if (!strcmp(ky, "MAT")) // FORMAT
                    return LDOC_NDE_OO;
            }
        }
    }
    else if (*ky == 'I')
    {
        ky++;
        if (*ky == 'N')
        {
            ky++;
            if (*ky == 'F')
            {
                ky++;
                if (!strcmp(ky, "O")) // INFO
                    return LDOC_NDE_OO;
            }
        }
    }
    else if (*ky == 'P')
    {
        ky++;
        if (*ky == 'E')
        {
            ky++;
            if (*ky == 'D')
            {
                ky++;
                if (!strcmp(ky, "IGREE")) // PEDIGREE
                    return LDOC_NDE_OO;
                else if (!strcmp(ky, "IGREEDB")) // PEDIGREEDB (see also pedigreeDB)
                    return LDOC_NDE_OO;
            }
        }
    }
    else if (*ky == 'S')
    {
        ky++;
        if (*ky == 'A')
        {
            ky++;
            if (*ky == 'M')
            {
                ky++;
                if (!strcmp(ky, "PLE")) // SAMPLE
                    return LDOC_NDE_OO;
            }
        }
    }
    else if (*ky == 'a')
    {
        ky++;
        if (*ky == 's')
        {
            ky++;
            if (*ky == 's')
            {
                ky++;
                if (!strcmp(ky, "embly")) // assembly
                    return LDOC_NDE_UA;
            }
        }
    }
    else if (*ky == 'c')
    {
        ky++;
        if (*ky == 'o')
        {
            ky++;
            if (*ky == 'n')
            {
                ky++;
                if (!strcmp(ky, "tig")) // contig
                    return LDOC_NDE_OO;
            }
        }
    }
    else if (*ky == 'f')
    {
        ky++;
        if (*ky == 'i')
        {
            ky++;
            if (*ky == 'l')
            {
                ky++;
                if (!strcmp(ky, "eDate") || !strcmp(ky, "edate")) // fileDate & filedate
                    return LDOC_NDE_UA;
                if (!strcmp(ky, "eformat")) // fileformat
                    return LDOC_NDE_UA;
            }
        } else if (*ky == 'f')
        {
            ky++;
            if (*ky == 'f')
            {
                ky++;
                if (!strcmp(ky, "-version"))
                    return LDOC_NDE_UA;
            }
        }
    }
    else if (*ky == 'p')
    {
        ky++;
        if (*ky == 'e')
        {
            ky++;
            if (*ky == 'd')
            {
                ky++;
                if (!strcmp(ky, "igreeDB")) // pedigreeDB (see also PEDIGREEDB)
                    return LDOC_NDE_OO;
            }
        }
        else if (*ky == 'h')
        {
            ky++;
            if (*ky == 'a')
            {
                ky++;
                if (!strcmp(ky, "sing")) // phasing
                    return LDOC_NDE_UA;
            }
        }
    }
    else if (*ky == 'r')
    {
        ky++;
        if (*ky == 'e')
        {
            ky++;
            if (*ky == 'f')
            {
                ky++;
                if (!strcmp(ky, "erence")) // reference
                    return LDOC_NDE_UA;
            }
        }
    }
    else if (*ky == 's')
    {
        ky++;
        if (*ky == 'o')
        {
            ky++;
            if (*ky == 'u')
            {
                ky++;
                if (!strcmp(ky, "rce")) // source
                    return LDOC_NDE_UA;
            }
        }
    }
    
    return LDOC_NDE_OL;
}

static inline void vcf_proc_brckt(ldoc_nde_t* cntnr, char* id, char* inf)
{
    if (*inf != '<')
    {
        // TODO Format error.
    }
    
    inf++;
    
    ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_UA);
    
    id = strdup(id);
    gen_lwr(id);
    nde->mkup.anno.str = id;
    
    // TODO Does not deal with escaped '"' within strings ('"' opened/closed)
    ldoc_ent_t* ent;
    bool cnt = true;
    char* ky = inf;
    char* val;
    bool strng = false;
    while (cnt)
    {
        if (!strng && *inf == '=')
        {
            *inf = 0;
            val = inf + 1;
            
            if (*val == '"')
                strng = true;
        }
        else if (strng && inf > val && *inf == '"')
        {
            strng = false;
        }
        else if (!strng && (*inf == ',' || *inf == '>'))
        {
            if (*inf == '>')
                cnt = false;
            
            *inf = 0;

            // Remove quotes for known keys:
            if (!strcmp(ky, "Description"))
            {
                // Make sure there are quotes -- and at least a character:
                if (*val == '"')
                {
                    char* tmp = val;
                    while (*(++tmp));
                    
                    *(--tmp) = 0;
                    
                    val++;
                }
            }
            
            // Determine entity type based on key:
            ldoc_content_t tpe;
            if (!strcmp(ky, "ID"))
                tpe = LDOC_ENT_OR;
            else
                tpe = gen_smrt_tpe(val);
            
            // Make key "JSON" like (lower case):
            ky = strdup(ky);
            gen_lwr(ky);
            
            ent = ldoc_ent_new(tpe);
            
            // TODO Error handling.

            ent->pld.pair.anno.str = ky;
            
            // Handle "unknown":
            if (*val == '.' && !*(val + 1))
                ent->pld.pair.dtm.str = NULL;
            else
                ent->pld.pair.dtm.str = strdup(val);
            
            ldoc_nde_ent_push(nde, ent);
            
            ky = inf + 1;
        }
        else if (!*inf)
            cnt = false;
        
        inf++;
    }
    
    ldoc_nde_dsc_push(cntnr, nde);
}

static inline ldoc_doc_t* vcf_proc_prgm(ldoc_doc_t* doc, char* ln, size_t lnlen, char** cmt)
{
    // Skip leading markup:
    ln += 2;
    
    // Remove trailing line breaks and white space (account for skip ahead above):
    off_t off = lnlen - 3;
    while ((ln[off] == '\n' || ln[off] == '\r' || ln[off] == ' ' || ln[off] == '\t') && off > 0)
        ln[off--] = 0;
    
    // Separate key and value:
    char* val = ln;
    while (*val &&
           *val != '=')
        val++;
    
    // Value empty? In any case: terminate key string:
    if (!*val)
    {
        val = NULL;
    }
    else
    {
        *val = 0;
        val++;
        
        // Skip trailing whitespace in values:
        while (*val == ' ' ||
               *val == '\t')
            val++;
    }
    
    ldoc_nde_t* stmt = NULL;
    TAILQ_FOREACH(stmt, &(doc->rt->dscs), ldoc_nde_entries)
    {
        if (!strcmp(stmt->mkup.anno.str, ln))
            break;
    }
    
    ldoc_struct_t tpe = vcf_prgm_tpe(ln);
    if (!stmt)
    {
        if (tpe == LDOC_NDE_UA)
        {
            ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_OR);
            ent->pld.pair.anno.str = strdup(ln);
            ent->pld.pair.dtm.str = strdup(val);
            
            ldoc_nde_ent_push(doc->rt, ent);
        }
        else
        {
            // TODO Use own types, so that this conversion is not necessary:
            stmt = ldoc_nde_new(tpe == LDOC_NDE_OO ? LDOC_NDE_OL : tpe);
            stmt->mkup.anno.str = strdup(ln);
            
            ldoc_nde_dsc_push(doc->rt, stmt);
        }
    }
    
    // Possible identifier for use in the comment section:
    char* id = NULL;
    
    if (tpe == LDOC_NDE_OO)
    {
        // Fine with string comparisons here, since this case
        // will (hopefully) not be true many times.
        if (!strcmp(ln, "FILTER"))
        {
            gen_lwr(ln);
            vcf_proc_brckt(doc->rt, ln, val);
        }
        else if (!strcmp(ln, "INFO"))
        {
            gen_lwr(ln);
            vcf_proc_brckt(doc->rt, ln, val);
        }
        else if (!strcmp(ln, "FORMAT"))
        {
            gen_lwr(ln);
            vcf_proc_brckt(doc->rt, ln, val);
        }
        else if (!strcmp(ln, "contig"))
        {
            vcf_proc_brckt(doc->rt, ln, val);
        }
        else if (!strcmp(ln, "SAMPLE"))
        {
            gen_lwr(ln);
            vcf_proc_brckt(doc->rt, ln, val);
        }
        else if (!strcmp(ln, "ALT"))
        {
            gen_lwr(ln);
            vcf_proc_brckt(doc->rt, ln, val);
        }
        else if (!strcmp(ln, "PEDIGREE"))
        {
            gen_lwr(ln);
            vcf_proc_brckt(doc->rt, ln, val);
        }
        else if (!strcmp(ln, "PEDIGREEDB") || !strcmp(ln, "pedigreeDB"))
        {
            
        }
        else
        {
            // TODO Internal error.
        }
    }
    else if (tpe == LDOC_NDE_OL)
    {
        if (!strcmp(ln, "genome-build"))
        {
            gff_proc_gbld(stmt, val);
        }
        else
        {
            // Unknown pragmas:
            ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_TXT);
            
            // TODO Error handling.
            
            ent->pld.str = strdup(val);
            
            // TODO Error handling.
            
            ldoc_nde_ent_push(stmt, ent);
        }
    }
    else if (tpe != LDOC_NDE_UA)
    {
        ldoc_ent_t* ntry = ldoc_ent_new(LDOC_ENT_OR);
        ntry->pld.str = strdup(val);
        
        ldoc_nde_ent_push(stmt, ntry);
    }
    
    // Add comments in a separate section:
    if (*cmt)
    {
        // Find previous comment section (if it exists):
        ldoc_nde_t* c = NULL;
        ldoc_nde_t* iter;
        TAILQ_FOREACH(iter, &(doc->rt->dscs), ldoc_nde_entries)
        {
            if (!strcmp(iter->mkup.anno.str, GEN_COMMENT))
            {
                c = iter;
                break;
            }
        }
        
        // If no comment section exists, create one now:
        if (!c)
        {
            c = ldoc_nde_new(LDOC_NDE_OO);
            c->mkup.anno.str = strdup(GEN_COMMENT);
            
            ldoc_nde_dsc_push(doc->rt, c);
        }
        
        // Find whether the meta information has been commented on already:
        ldoc_nde_t* cpgm = NULL;
        TAILQ_FOREACH(iter, &(c->dscs), ldoc_nde_entries)
        {
            if (!strcmp(iter->mkup.anno.str, ln))
            {
                cpgm = iter;
                break;
            }
        }
        
        // No comment for the meta information present yet? Create a new node:
        if (!cpgm)
        {
            cpgm = ldoc_nde_new(LDOC_NDE_OL);
            cpgm->mkup.anno.str = strdup(ln);
            
            ldoc_nde_dsc_push(c, cpgm);
            
        }
        
        // Now add the actual comment:
        ldoc_ent_t* cent = ldoc_ent_new(LDOC_ENT_OR);;
        if (id)
        {
            cent->pld.pair.anno.str = id;
            cent->pld.pair.dtm.str = gen_escstr(*cmt, GEN_FMT_VCF);
        }
        else
        {
            cent->pld.pair.anno.str = strdup(val);
            cent->pld.pair.dtm.str = gen_escstr(*cmt, GEN_FMT_VCF);
            
        }
        ldoc_nde_ent_push(cpgm, cent);
        
        // Erase comment, but it is not released (free'd) yet:
        *cmt = NULL;
    }
    
    return doc;
}

static inline ldoc_nde_t* vcf_proc_gle(char* val, size_t len, char* ref, char** vseqs, size_t vnum)
{
    /*
     On May 28, 2015 at 11:06:25 AM, Erik Garrison (erik.garrison@gmail.com) wrote:
     
     
     I contributed this and think it is broken (and unused). It is now mostly dealt with by the extension of the GL field to allow for arbitrary ploidy.
     */

    // Example:
    // GLE=0:-75.22,1:-223.42,0/0:-323.03,1/0:-99.29,1/1:-802.53
    
    ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_UA);
    
    nde->mkup.anno.str = (char*)VCF_GENOLHOODE;
    
    size_t n = 0;
    char* v = val;
    
    // Need to be large enough to account for the values of idx1 and idx2:
    uint16_t idx1_, idx2_;
    char idx1[GEN_MAX_ALT_CHARS + 1];
    char idx2[GEN_MAX_ALT_CHARS + 1];
    while (len--)
    {
        // TODO Check whether val - v is larger than VCF_MAX_ALT_CHARS - 1.
        
        if (!len || *val == ',')
        {
            // Adjust val pointer, if len reached:
            if (!len)
                val++;
            
            // Find colon and optional forward slash:
            char* cln = v;
            char* slsh = NULL;
            while (*cln != ':' && cln < val)
            {
                if (*cln == '/')
                {
                    // TODO Could check format by testing whether slsh is NULL.
                    slsh = cln;
                    
                    *slsh = 0;
                }
                
                cln++;
            }
            
            if (cln >= val)
            {
                // TODO Format error.
            }
            
            *cln = 0;
            
            // Convert indexes:
            size_t ilen = slsh ? slsh - v : cln - v;
            memcpy(idx1, v, ilen);
            idx1[ilen] = 0;
            idx1_ = 2 * strtol(idx1, NULL, 10);
            if (slsh)
            {
                ilen = cln - (slsh + 1);
                memcpy(idx2, slsh + 1, ilen);
                idx2[ilen] = 0;
                idx2_ = 2* strtol(idx2, NULL, 10);
            }
            else
                idx2_ = 1; // Null byte.
            
            char gt[3] = { &GEN_ALLELE[idx1_], &GEN_ALLELE[idx2_], 0 };
            
            ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_NR);
            
            // TODO Ignore this strdup for now. Since GLE most likely unused.
            ent->pld.pair.anno.str = strdup(gt);
            
            if (val - (cln + 1) == 1 && *(cln + 1) == '.')
                ent->pld.pair.dtm.str = NULL;
            else
                ent->pld.pair.dtm.str = strndup(v, val - (cln + 1)); // TODO Same, ignore for now. GLE unlikely to be used.
            
            ldoc_nde_ent_push(nde, ent);
            
            n++;
            v = val + 1;
        }
        
        val++;
    }
    
    return nde;
}

static inline ldoc_nde_t* vcf_proc_glpl(char* val, size_t len, bool pl)
{
    ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_UA);
    
    if (pl)
        nde->mkup.anno.str = (char*)VCF_PHREDLHOOD;
    else
        nde->mkup.anno.str = (char*)VCF_GENOLHOOD;
    
    size_t slen;
    size_t n = 0;
    char* v = val;
    while (len--)
    {
        // TODO Check whether val - v is larger than VCF_MAX_ALT_CHARS - 1.
        
        if (!len || *val == ',')
        {
            ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_NR);
            
            ent->pld.pair.anno.str = &GEN_ALLELES[n * GEN_STEP];
            
            slen = len ? val - v : val - v + 1;
            if (slen == 1 && *v == '.')
                ent->pld.pair.dtm.str = NULL;
            else
                ent->pld.pair.dtm.str = qk_strndup(v, slen);
            
            ldoc_nde_ent_push(nde, ent);
            
            n++;
            v = val + 1;
        }
        
        val++;
    }
    
    return nde;
}

static inline ldoc_nde_t* vcf_proc_gt(char* val, size_t len, char* ref, char** vseqs, size_t vnum)
{
    ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_UA);
    ldoc_nde_t* nde_seq = ldoc_nde_new(LDOC_NDE_OL);
    
    nde->mkup.anno.str = (char*)VCF_GENOTYPE;
    
    nde_seq->mkup.anno.str = (char*)GEN_SEQUENCE;
    
    char* v = val;
    uint16_t idx_; // Needs to be large enough to account for the value in idx.
    char idx[GEN_MAX_ALT_CHARS + 1];
    bool phased = false;
    char gt[len + 1]; // Genotype alleles. Size overapproximation.
    char* g = gt;
    while (len--)
    {
        // TODO Check whether val - v is larger than VCF_MAX_ALT_CHARS - 1.
        
        if (!len || *val == '/' || *val == '|')
        {
            ldoc_ent_t* ent_seq = ldoc_ent_new(LDOC_ENT_TXT);
            
            if (*val == '|')
                phased = true;
            
            if (!len)
                val++;
            
            // Isolate number:
            memcpy(idx, v, val - v);
            idx[val - v] = 0;
            
            // Convert number string to an actual number:
            idx_ = strtol(idx, NULL, 10);
            
            if (!idx_)
            {
                ent_seq->pld.str = ref;
                *(g++) = 'A';
            }
            else
            {
                ent_seq->pld.str = vseqs[idx_ - 1];
                *(g++) = 'A' + idx_;
            }
            
            ldoc_nde_ent_push(nde_seq, ent_seq);
            
            v = val + 1;
        }
        
        val++;
    }
    *g = 0;
    
    ldoc_nde_dsc_push(nde, nde_seq);

    ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_OR);
    
    ent->pld.pair.anno.str = (char*)VCF_ALLELES;
    ent->pld.pair.dtm.str = qk_strdup(gt);
    
    ldoc_nde_ent_push(nde, ent);
    
    ent = ldoc_ent_new(LDOC_ENT_OR);
    
    ent->pld.pair.anno.str = (char*)VCF_PHASED;
    ent->pld.pair.dtm.str = phased ? "true" : "false";
    
    ldoc_nde_ent_push(nde, ent);
    
    return nde;
}

static inline ldoc_ent_t* vcf_proc_num(char* val, size_t len, char* lbl, size_t lbllen)
{
    ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_NR);
    
    ent->pld.pair.anno.str = qk_strndup(lbl, lbllen);
    
    if (len == 1 && *val == '.')
        ent->pld.pair.dtm.str = NULL;
    else
        ent->pld.pair.dtm.str = qk_strndup(val, len);
    
    return ent;
}

static inline void vcf_proc_smpl(ldoc_nde_t* smpls, gen_prsr_t* stt, size_t i, char* fmt, char* smpl, char* ref, char** vseqs, size_t vnum)
{
    // Note: Going through fmt character-by-character, because the names
    //       are so short that it might be as efficient as looking up a
    //       fixed format that *needs* to be the same for the whole file.
    
    // TODO Assumes that fmt is not empty.
    
    ldoc_nde_t* s = ldoc_nde_new(LDOC_NDE_UA);
    
    s->mkup.anno.str = qk_strdup(stt->vcf_hdr_off[i]);
    
    // TODO Error handling.

    char* val;
    do
    {
        val = smpl;
        
        // Find label from format string:
        char* id = fmt;
        while (*fmt && *fmt != ':')
            fmt++;

        // Find value part:
        while (*smpl && *smpl != ':')
            smpl++;
        
        // TODO Optimization.
        //   1. could be optimized by not checking length multiple times.
        //   2. could be optimized by not checking prefixes multiple times.
        if (fmt - id == 2 && id[0] == 'G' && id[1] == 'T')
        {
            // GT: genotype
            ldoc_nde_t* gt = vcf_proc_gt(val, smpl - val, ref, vseqs, vnum);
            
            ldoc_nde_dsc_push(s, gt);
        }
        else if (fmt - id == 2 && id[0] == 'G' && id[1] == 'L')
        {
            // PL: phred-scaled genotype likelihoods
            ldoc_nde_t* gl = vcf_proc_glpl(val, smpl - val, false);
            
            ldoc_nde_dsc_push(s, gl);
        }
        else if (fmt - id == 3 && id[0] == 'G' && id[1] == 'L' && id[2] == 'E')
        {
            // GLE: genotype likelihoods
            ldoc_nde_t* gt = vcf_proc_gle(val, smpl - val, ref, vseqs, vnum);
            
            ldoc_nde_dsc_push(s, gt);
        }
        else if (fmt - id == 2 && id[0] == 'P' && id[1] == 'L')
        {
            // PL: phred-scaled genotype likelihoods
            ldoc_nde_t* pl = vcf_proc_glpl(val, smpl - val, true);
            
            ldoc_nde_dsc_push(s, pl);
        }
        else if ((fmt - id == 2 &&
                  ((id[0] == 'A' && id[1] == 'N') ||
                   (id[0] == 'B' && id[1] == 'Q') ||
                   (id[0] == 'D' && id[1] == 'P') ||
                   (id[0] == 'G' && id[1] == 'P') || // TODO Double check if this is not a comma separated list.
                   (id[0] == 'G' && id[1] == 'Q') ||
                   (id[0] == 'M' && id[1] == 'Q') ||
                   (id[0] == 'N' && id[1] == 'S') ||
                   (id[0] == 'P' && id[1] == 'Q'))) ||
                 (fmt - id == 3 &&
                  ((id[0] == 'E' && id[1] == 'N' && id[2] == 'D') ||
                   (id[0] == 'M' && id[1] == 'Q' && id[2] == '0'))))
        {
            // Numeric fields with a single value.
            ldoc_ent_t* num = vcf_proc_num(val, smpl - val, id, fmt - id);
            
            ldoc_nde_ent_push(s, num);
        }
        else
        {
            ldoc_ent_t* anno = ldoc_ent_new(LDOC_ENT_OR);
            
            // TODO Error handling.

            // Default: add simple key/value pair
            anno->pld.pair.anno.str = qk_strndup(id, fmt - id);
            anno->pld.pair.dtm.str = qk_strndup(val, smpl - val);
            
            ldoc_nde_ent_push(s, anno);
        }
        
        if (*smpl)
            smpl++;
    } while (*(fmt++));
    
    ldoc_nde_dsc_push(smpls, s);
}

static inline void vcf_proc_optlst(ldoc_nde_t* ftr, char* id, char* lst)
{
    ldoc_ent_t* ent;
    
    // Check whether there is no list:
    if (!lst || (lst[0] == '.' && !lst[1]))
    {
        ent = ldoc_ent_new(LDOC_ENT_OR);
        
        // TODO Error handling.
        
        ent->pld.pair.anno.str = (char*)id;
        ent->pld.pair.dtm.str = NULL;
        
        ldoc_nde_ent_push(ftr, ent);
        
        return;
    }
    
    ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_OL);
    
    // TODO Error handling.
    
    nde->mkup.anno.str = id;
    
    char* f;
    bool cnt = true;
    while (cnt)
    {
        f = lst;
        
        while (*lst && *lst != ';')
            lst++;

        if (!*lst)
            cnt = false;
        else
            *(lst++) = 0;
        
        ent = ldoc_ent_new(LDOC_ENT_TXT);
        
        // TODO Error handling.
        
        ent->pld.str = f;
        
        ldoc_nde_ent_push(nde, ent);
    }
    
    ldoc_nde_dsc_push(ftr, nde);
}

static inline ldoc_doc_t* vcf_proc_ftr(int fd, off_t mx, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* stt, char** cmt)
{
    ldoc_doc_t* doc = ldoc_doc_new();
    
    if (!doc)
    {
        // TODO Error handling.
    }
    
    // Remove trailing line breaks:
    off_t off = lnlen - 1;
    while ((ln[off] == '\n' || ln[off] == '\r') && off > 0)
        ln[off--] = 0;
    
    // Note: this can be optimized by creating the entities in the for-loop!
    // TODO Verify whether this creates stack overflow problems when there are
    //      lots of samples in a VCF file.
    off = 0;
    char* coff[stt->vcf_col];
    coff[0] = ln;
    for (uint8_t col = 1; col <= stt->vcf_col; col++)
    {
        while (ln[off] != '\t' && off < lnlen)
            off++;
        
        if (off >= lnlen)
        {
            // TODO Error handling.
        }
        
        ln[off] = 0;
        coff[col] = ln + off + 1;
        
        // Handling of "unknown" values ("." values):
        // Note: more efficient than using gff_proc_optvalx below.
        if (coff[col - 1][0] == '.' && (off + 1 >= lnlen || coff[col - 1][1] == 0))
            coff[col - 1] = NULL;
    }
    
    ldoc_nde_t* ftr = doc->rt;
    
    /* Reimplement later where only used "things" are defined:
     ldoc_nde_t* ctx = ldoc_nde_new(LDOC_NDE_UA);
     ctx->mkup.anno.str = (char*)JSONLD_CTX;
     */
    
    ldoc_ent_t* ctx = ldoc_ent_new(LDOC_ENT_OR);
    ctx->pld.pair.anno.str = (char*)JSONLD_CTX;
    ctx->pld.pair.dtm.str = (char*)JSONLD_VCF;
    
    ldoc_nde_t* lc = ldoc_nde_new(LDOC_NDE_UA);
    lc->mkup.anno.str = (char*)GEN_LOCUS;

    ldoc_ent_t* lm = ldoc_ent_new(LDOC_ENT_OR);
    lm->pld.pair.anno.str = (char*)VCF_C1;
    lm->pld.pair.dtm.str = coff[0];
    
    ldoc_ent_t* st = ldoc_ent_new(LDOC_ENT_NR);
    st->pld.pair.anno.str = (char*)VCF_C2_1;
    st->pld.pair.dtm.str = coff[1];
    
    ldoc_ent_t* en = ldoc_ent_new(LDOC_ENT_NR);
    en->pld.pair.anno.str = (char*)VCF_C2_2;
    en->pld.pair.dtm.str = coff[1];
    
    vcf_proc_optlst(ftr, (char*)VCF_C3, coff[2]);

    ldoc_nde_t* ref = ldoc_nde_new(LDOC_NDE_UA);
    ref->mkup.anno.str = (char*)VCF_C4;
    
    ldoc_ent_t* ref_seq = ldoc_ent_new(LDOC_ENT_OR);
    ref_seq->pld.pair.anno.str = (char*)GEN_SEQUENCE;
    ref_seq->pld.pair.dtm.str = coff[3];
    ldoc_nde_ent_push(ref, ref_seq);
    
    // Note: assume that no more than VCF_MAX_ALT variants are observed:
    size_t vnum = 0;
    char* vseqs[GEN_MAX_ALT];
    ldoc_nde_t* vars = gen_variants(coff[4], ',', vseqs, &vnum);
    
    ldoc_ent_t* scr = ldoc_ent_new(LDOC_ENT_NR);
    scr->pld.pair.anno.str = (char*)VCF_C6;
    scr->pld.pair.dtm.str = coff[5];
    
    vcf_proc_optlst(ftr, (char*)VCF_C7, coff[6]);
    
    // Add comment lines -- if available:
    if (*cmt)
    {
        ldoc_ent_t* c = ldoc_ent_new(LDOC_ENT_OR);
        c->pld.pair.anno.str = (char*)GEN_COMMENT;
        c->pld.pair.dtm.str = gen_escstr(*cmt, GEN_FMT_VCF);
        ldoc_nde_ent_push(ftr, c);
        
        // Erase comment, but it is not released (free'd) yet:
        *cmt = NULL;
    }
    
    ldoc_nde_t* attrs = ldoc_nde_new(LDOC_NDE_UA);
    attrs->mkup.anno.str = (char*)GEN_ATTRS;
    
    gen_splt_attrs(ftr, attrs, ref, vars, coff[7], BI_NKW);
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_nde_ent_push(ftr, ctx);
    
    // Score:
    ldoc_nde_ent_push(ftr, scr);
    
    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, lm);
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    
    ldoc_nde_dsc_push(ftr, lc);
    
    // Reference & variants:
    ldoc_nde_dsc_push(ftr, ref);
    ldoc_nde_dsc_push(ftr, vars);
    
    // Do not add user defined sub-tree if it is empty:
    if (attrs->dsc_cnt || attrs->ent_cnt)
        ldoc_nde_dsc_push(ftr, attrs);
    
    // Sample information -- if provided:
    if (stt->vcf_col > 8)
    {
        ldoc_nde_t* smpls = ldoc_nde_new(LDOC_NDE_UA);
        
        // TODO Error handling.
        
        smpls->mkup.anno.str = (char*)VCF_SAMPLES;
        
        for (size_t i = 9; i < stt->vcf_col; i++)
        {
            if (coff[i])
                vcf_proc_smpl(smpls, stt, i - 9, coff[8], coff[i], coff[3], vseqs, vnum);
            else
            {
                ldoc_ent_t* smpl_null = ldoc_ent_new(LDOC_ENT_OR);
                
                // TODO Error handling.
                
                smpl_null->pld.pair.anno.str = qk_strdup(stt->vcf_hdr_off[i - 9]);
                smpl_null->pld.pair.dtm.str = NULL;
                
                ldoc_nde_ent_push(smpls, smpl_null);
            }
        }
        
        ldoc_nde_dsc_push(ftr, smpls);
    }
    
    return doc;
}

ldoc_doc_t* vcf_proc_ln(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt, gen_fstat* stat)
{
    ldoc_doc_t* ldoc = NULL;
    
    if (lnlen > 3)
    {
        if (ln[0] == '#')
        {
            if (ln[1] == '#')
            {
                // VCF file have no FASTA section, so this is easier than
                // GFF3/GVF handling.
                
                // Meta line (pragma statement):
                vcf_proc_prgm(fdoc, ln, lnlen, cmt);
                stat->meta++;
            }
            else
            {
                if (lnlen > 7 && !strncmp(ln, "#CHROM\t", 7))
                {
                    // TODO If already 'true', raise error.
                    st->vcf_ftr_sct = true;
                    
                    // Save values for second parse that accumulates sample IDs:
                    char* ln_ = ln;
                    size_t lnlen_ = lnlen;
                    
                    // Note that this works, because we are still on the '##...' prefix:
                    st->vcf_col = 0;
                    while(*(ln++) && lnlen--)
                        if (*ln == '\t')
                            st->vcf_col++;
                    st->vcf_col++;
                    
                    if (st->vcf_col < 8)
                    {
                        // TODO Error handling. Too few columns for a VCF file.
                    }
                    
                    // Second parse: accumulate sample IDs:
                    
                    st->vcf_hdr = strdup(ln_);
                    
                    // TODO Error handling.
                    
                    ln_ = st->vcf_hdr;
                    
                    st->vcf_hdr_off = (char**)malloc(sizeof(char*) * st->vcf_col + 1 - 9);
                    
                    // TODO Error handling.
                    
                    // Note: Pre-increment possible, since we are still on '##...':
                    size_t col = 0;
                    while(*(++ln_) && --lnlen_)
                        if (*ln_ == '\t' || *ln_ == '\n' || !lnlen_)
                        {
                            *ln_ = 0;
                            
                            if (col >= 8)
                                st->vcf_hdr_off[col - 8] = ln_ + 1;
                            
                            col++;
                        }
                    *ln_ = 0;
                }
                else
                {
                    // Comment -- note that officially VCF files do not have comments, so re-use GFF3 implementation:
                    char* ccmt = gff_proc_cmt(ln, lnlen);
                    stat->comms++;
                    
                    // Append to existing comment, or, create a new one:
                    if (*cmt)
                    {
                        size_t ccmtlen = strlen(ccmt);
                        size_t cmtlen = strlen(*cmt) + ccmtlen;
                        *cmt = realloc(*cmt, cmtlen + 1);
                        
                        // TODO Error handling.
                        
                        strncat(*cmt, ccmt, ccmtlen);
                        
                        free(ccmt);
                    }
                    else
                        *cmt = ccmt;
                }
            }
        }
        else if (st->vcf_ftr_sct)
        {
            // Feature:
            ldoc = vcf_proc_ftr(fd, mx, idx, ln, lnlen, st, cmt);
            stat->ftrs++;
        }
        else
        {
            // TODO Error handling. Corrupt VCF file.
        }
    }
    
    return ldoc;
}

static inline void vcf_proc_doc_optlst(ldoc_nde_t* nde, char* id, char* empty)
{
    // If filters are an entity, then the assigned value should be NULL!
    ldoc_res_t* lst = ldoc_find_anno_ent(nde, id);
    if (lst)
    {
        if (lst->info.ent->pld.pair.dtm.str)
        {
            // TODO Data error.
        }
        
        qk_strcat(GEN_UNKNOWN);
    }
    else
    {
        const char* lst_id[] = { id };
        lst = ldoc_find_anno_nde(nde, (char**)lst_id, 1);
        
        if (lst)
        {
            if (!lst->info.nde->ent_cnt)
                qk_strcat(empty);
            else
            {
                bool fst = true;
                ldoc_ent_t* fltr;
                TAILQ_FOREACH(fltr, &(lst->info.nde->ents), ldoc_ent_entries)
                {
                    if (fst)
                        fst = false;
                    else
                        qk_strcat(";");
                    
                    // TODO Check entity type. Needs to be LDOC_ENT_TXT.
                    qk_strcat(fltr->pld.str);
                }
            }
        }
        else
            qk_strcat(GEN_UNKNOWN);
    }
}

static inline bool vcf_proc_doc_info(ldoc_nde_t* ftr, char* attrs)
{
    ldoc_res_t* res = ldoc_find_anno_ent(ftr, (char*)GEN_ALLELE_CNT);
    if (res)
        gen_join_attrs_ent((char*)GEN_ALLELE_CNT_VCF, res->info.ent, attrs);

    res = ldoc_find_anno_ent(ftr, (char*)GEN_ALLELE_FRQ);
    if (res)
        gen_join_attrs_ent((char*)GEN_ALLELE_FRQ_VCF, res->info.ent, attrs);

    res = ldoc_find_anno_ent(ftr, (char*)GEN_ALLELE_TTL);
    if (res)
        gen_join_attrs_ent((char*)GEN_ALLELE_TTL_VCF, res->info.ent, attrs);
    
    return true;
}

inline char* vcf_proc_doc_ftr(ldoc_nde_t* ftr)
{
    const char* lc_id[] = { GEN_LOCUS };
    ldoc_res_t* lc = ldoc_find_anno_nde(ftr, (char**)lc_id, 1);
    
    ldoc_res_t* lm = ldoc_find_anno_ent(lc->info.nde, (char*)VCF_C1);
    ldoc_res_t* st = ldoc_find_anno_ent(lc->info.nde, (char*)VCF_C2_1);
    ldoc_res_t* en = ldoc_find_anno_ent(lc->info.nde, (char*)VCF_C2_2);
    
    // TODO Check that st and en value are the same!
    
    const char* ref_id[] = { VCF_C4 };
    ldoc_res_t* ref = ldoc_find_anno_nde(ftr, (char**)ref_id, 1);
    ldoc_res_t* ref_seq = ldoc_find_anno_ent(ref->info.nde, (char*)GEN_SEQUENCE);
    
    qk_strcat(gen_res_req(lm));
    qk_strcat("\t");
    qk_strcat(gen_res_req(st));
    qk_strcat("\t");
    vcf_proc_doc_optlst(ftr, (char*)VCF_C3, (char*)GEN_UNKNOWN);
    qk_strcat("\t");
    qk_strcat(gen_res_opt(ref_seq));
    qk_strcat("\t");
    
    const char* vars_id[] = { GEN_VARIANTS };
    ldoc_res_t* vars = ldoc_find_anno_nde(ftr, (char**)vars_id, 1);
    off_t allele = 0;
    ldoc_res_t* var;
    ldoc_res_t* seq;
    char* str;
    while (true) // See "Break condition" below.
    {
        const char* var_id[] = { &GEN_ALLELE[(allele + 1) * 2] };
        var = ldoc_find_anno_nde(vars->info.nde, (char**)var_id, 1);

        if (!var) // Break condition.
            break;
        
        seq = ldoc_find_anno_ent(var->info.nde, (char*)GEN_SEQUENCE);
        
        if (!seq)
            str = ".";
        else
            str = seq->info.ent->pld.pair.dtm.str;
        
        if (allele)
            qk_strcat(",");
        
        qk_strcat(str);
        
        allele++;
    }
    qk_strcat("\t");
    
    // SCORE column:
    ldoc_res_t* scr = ldoc_find_anno_ent(ftr, (char*)VCF_C6);
    qk_strcat(gen_res_opt(scr));
    qk_strcat("\t");

    // FILTER column:
    vcf_proc_doc_optlst(ftr, (char*)VCF_C7, "PASS");
    qk_strcat("\t");

    // INFO -- standards part:
    char* attrs = qk_working_ptr();
    vcf_proc_doc_info(ftr, attrs);
    
    // INFO -- user-defined part:
    const char* usr_pth[] = { GEN_ATTRS };
    ldoc_res_t* usr = ldoc_find_anno_nde(ftr, (char**)usr_pth, 1);
    if (usr)
    {
        ldoc_ent_t* ent;
        TAILQ_FOREACH(ent, &(usr->info.nde->ents), ldoc_ent_entries)
        {
            gen_join_attrs_ent(NULL, ent, attrs);
        }
    }
    
    /*
    ldoc_res_t* nme = ldoc_find_anno_ent(ftr, "name");
    
    const char* dbxref_id[] = { "dbxref" };
    ldoc_res_t* dbxref = ldoc_find_anno_nde(ftr, (char**)dbxref_id, 1);
    
    const char* prnt_pth[] = { "parent" };
    ldoc_res_t* prnt = ldoc_find_anno_nde(ftr, (char**)prnt_pth, 1);
    
    char* attrs = qk_working_ptr();
    
    if (id && !id->nde)
        gen_join_attrs_ent("ID", id->info.ent, attrs);
    
    if (nme && !nme->nde)
        gen_join_attrs_ent("Name", nme->info.ent, attrs);
    
    if (prnt && prnt->nde)
        gen_join_attrs_nde("Parent", prnt->info.nde, attrs);
    
    if (dbxref && dbxref->nde)
        gen_join_attrs_nde("Dbxref", dbxref->info.nde, attrs);
    */
    
    return NULL;
}

inline char* vcf_proc_doc_ftr_attrs(ldoc_nde_t* ftr)
{
    return NULL;
}

char* vcf_proc_doc(ldoc_doc_t* doc)
{
    char* attr = vcf_proc_doc_ftr_attrs(doc->rt);
    
    vcf_proc_doc_ftr(doc->rt);
    
    if (attr && *attr)
        qk_strcat(";");
    qk_strcat(attr);
    
    free(attr);
}

