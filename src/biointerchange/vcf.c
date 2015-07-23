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

    if (!strcmp(id, GEN_FMT_ALT_VCF))
        id = strdup(GEN_FMT_ALT);
    else if (!strcmp(id, GEN_FMT_INFO_VCF))
        id = strdup(GEN_FMT_INFO);
    else if (!strcmp(id, GEN_FMT_FLT_VCF))
        id = strdup(GEN_FMT_FLT);
    else if (!strcmp(id, GEN_FMT_GT_VCF))
        id = strdup(GEN_FMT_GT);
    else
    {
        id = strdup(id);
        gen_lwr(id);
    }

    char* pth[] = { id };
    ldoc_res_t* res = ldoc_find_anno_nde(cntnr, pth, 1);
    
    ldoc_nde_t* nde;
    if (res)
        nde = res->info.nde;
    else
        nde = ldoc_nde_new(LDOC_NDE_UA);
    
    nde->mkup.anno.str = id;
    
    ldoc_nde_t* nde_brckt = ldoc_nde_new(LDOC_NDE_UA);
    
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
            {
                gen_attr_t kwd;
                gen_kwd(val, &kwd, BI_NKW);
                
                if (kwd.alt)
                    nde_brckt->mkup.anno.str = (char*)kwd.alt;
                else
                    nde_brckt->mkup.anno.str = strdup(val);
                
                ky = inf + 1;
            }
            else
            {
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
                
                ldoc_nde_ent_push(nde_brckt, ent);
                
                ky = inf + 1;
            }
        }
        else if (!*inf)
            cnt = false;
        
        inf++;
    }
    
    ldoc_nde_dsc_push(nde, nde_brckt);
    
    if (!res)
        ldoc_nde_dsc_push(cntnr, nde);
}

static inline void vcf_proc_optlst(ldoc_nde_t* ftr, char* id, char* lst)
{
    ldoc_ent_t* ent;
    
    // Check whether there is no list:
    if (!lst ||
        (lst[0] == '.' && !lst[1]) ||
        !strcmp(lst, "PASS"))
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

static inline void vcf_proc_idlst(ldoc_nde_t* ftr, char* lst)
{
    // Find first semi-colon (if exists):
    if (lst)
    {
        char* sep = lst;
        while (*sep && *sep != ';')
            sep++;
        
        // Split second, third, etc., IDs into aliases:
        if (*sep)
        {
            *sep = 0;
            
            vcf_proc_optlst(ftr, (char*)GEN_ALIAS, sep + 1);
        }
    }
    
    // The actual ID (first entry):
    ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_OR);
    ent->pld.pair.anno.str = (char*)GEN_ID;
    ent->pld.pair.dtm.str = lst;
    
    ldoc_nde_ent_push(ftr, ent);
}

static inline ldoc_doc_t* vcf_proc_prgm(ldoc_doc_t* doc, char* ln, size_t lnlen, char** cmt)
{
    bool usr_nw;
    ldoc_nde_t* usr = gen_ctx(doc->rt, &usr_nw, JSONLD_VCF_X1);
    
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
    
    ldoc_struct_t tpe = vcf_prgm_tpe(ln);

    // LDOC_NDE_OO entries are <> ("bracket") structures, which are
    // handled independently by `vcf_proc_brckt`:
    ldoc_nde_t* stmt;
    if (tpe != LDOC_NDE_OO)
    {
        stmt = gen_find_nde(doc->rt, usr, ln);
        
        if (!stmt)
        {
            if (tpe == LDOC_NDE_UA)
            {
                ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_OR);
                
                // TODO Error handling.
                
                if (!strcmp(ln, GEN_VCFVERSION_VCF))
                {
                    ln = (char*)GEN_VCFVERSION;
                    
                    // Skip "VCFv" part of VCF versions:
                    if (!strncmp(val, "VCFv", 4))
                        val += 4;
                }
                else if (!strcmp(ln, GEN_FILE_DATE_VCF1) ||
                         !strcmp(ln, GEN_FILE_DATE_VCF2))
                {
                    ln = (char*)GEN_FILE_DATE;
                }
                else if (!strcmp(ln, GEN_FST_REF_VCF))
                {
                    ln = (char*)GEN_FST_REF;
                }
                else if (!strcmp(ln, GEN_FST_BKPNT_VCF))
                {
                    ln = (char*)GEN_FST_BKPNT;
                }
                
                ent->pld.pair.anno.str = strdup(ln);
                ent->pld.pair.dtm.str = strdup(val);
                
                ldoc_nde_ent_push(doc->rt, ent);
            }
            else
            {
                ldoc_nde_t* dst;
                if (tpe != LDOC_NDE_OL ||
                    !strcmp(ln, "X example X -- unused in VCF"))
                    dst = doc->rt;
                else
                    dst = usr; // User defined pragmas.
                
                // TODO Use own types, so that this conversion is not necessary:
                stmt = ldoc_nde_new(tpe == LDOC_NDE_OO ? LDOC_NDE_OL : tpe);
                stmt->mkup.anno.str = strdup(ln);
                
                ldoc_nde_dsc_push(dst, stmt);
            }
        }
    }
    
    // Possible identifier for use in the comment section:
    char* id = NULL;
    
    if (tpe == LDOC_NDE_OO)
    {
        // Fine with string comparisons here, since this case
        // will (hopefully) not be true many times.
        if (!strcmp(ln, "FILTER"))
            vcf_proc_brckt(doc->rt, ln, val);
        else if (!strcmp(ln, "INFO"))
            vcf_proc_brckt(doc->rt, ln, val);
        else if (!strcmp(ln, "FORMAT"))
            vcf_proc_brckt(doc->rt, ln, val);
        else if (!strcmp(ln, "contig"))
            vcf_proc_brckt(doc->rt, ln, val);
        else if (!strcmp(ln, "SAMPLE"))
            vcf_proc_brckt(doc->rt, ln, val);
        else if (!strcmp(ln, "ALT"))
            vcf_proc_brckt(doc->rt, ln, val);
        else if (!strcmp(ln, "PEDIGREE"))
            vcf_proc_brckt(doc->rt, ln, val);
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
        if (!strcmp(ln, "X example X -- unused in VCF"))
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
    
    // Add user-defined pragmas -- if those exist:
    if (usr_nw)
        gen_add_nw(doc->rt, usr);
    
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
                idx2_ = 2 * strtol(idx2, NULL, 10);
            }
            else
                idx2_ = 1; // Null byte.
            
            char gt[3] = { GEN_ALLELE[idx1_], GEN_ALLELE[idx2_], 0 };
            
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

static inline void vcf_proc_ec(ldoc_nde_t* cntnr, char* val, size_t len)
{
    const char* pth[] = { GEN_VARIANTS };
    ldoc_res_t* res = ldoc_find_anno_nde(cntnr, (char**)pth, 1);
    
    ldoc_nde_t* vars;
    if (res)
        vars = res->info.nde;
    else
    {
        vars = ldoc_nde_new(LDOC_NDE_UA);
        
        // TODO Error handling.
        
        vars->mkup.anno.str = (char*)GEN_VARIANTS;
        
        ldoc_nde_dsc_push(cntnr, vars);
    }
    
    ldoc_nde_t* nde;
    size_t slen;
    size_t n = 1;
    char* v = val;
    while (len--)
    {
        // TODO Check whether val - v is larger than VCF_MAX_ALT_CHARS - 1.
        
        if (!len || *val == ',')
        {
            char* alt[] = { &GEN_ALLELE[n * 2] };
            res = ldoc_find_anno_nde(vars, alt, 1);
            
            if (res)
                nde = res->info.nde;
            else
            {
                nde = ldoc_nde_new(LDOC_NDE_UA);
                
                // TODO Error handling.
                
                nde->mkup.anno.str = &GEN_ALLELE[n * 2];
                
                ldoc_nde_dsc_push(vars, nde);
            }
            
            ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_NR);
            
            ent->pld.pair.anno.str = (char*)GEN_ALLELE_CNTEXP;
            
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
}

static inline void vcf_proc_ft(ldoc_nde_t* cntnr, char* val, size_t len)
{
    vcf_proc_optlst(cntnr, (char*)GEN_ANNOTATIONS, val);
}

static inline void vcf_proc_hq(ldoc_nde_t* cntnr, char* val, size_t len, bool pl)
{
    gen_attr_t kwd = { BI_NKW, (char*)NULL };
    gen_csep_dup(cntnr, kwd, (char*)GEN_HAP_QUALITIES, val, false);
}

static inline void vcf_proc_glgppl(ldoc_nde_t* cntnr, char* val, size_t len, gen_alt_t alt)
{
    ldoc_nde_t* nde;
    size_t slen;
    size_t n = 0;
    char* v = val;
    while (len--)
    {
        // TODO Check whether val - v is larger than VCF_MAX_ALT_CHARS - 1.
        
        if (!len || *val == ',')
        {
            char* pth[] = { &GEN_ALLELES[n * GEN_STEP] };
            ldoc_res_t* res = ldoc_find_anno_nde(cntnr, pth, 1);
            
            if (res)
                nde = res->info.nde;
            else
            {
                nde = ldoc_nde_new(LDOC_NDE_UA);
                
                // TODO Error handling.
                
                nde->mkup.anno.str = &GEN_ALLELES[n * GEN_STEP];
                
                ldoc_nde_dsc_push(cntnr, nde);
            }
            
            ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_NR);
            
            switch (alt)
            {
                case GEN_ALT_GL:
                    ent->pld.pair.anno.str = (char*)GEN_GENOTYPE_LIKE;
                    break;
                case GEN_ALT_GP:
                    ent->pld.pair.anno.str = (char*)GEN_GENOTYPE_PROB;
                    break;
                case GEN_ALT_PL:
                    ent->pld.pair.anno.str = (char*)GEN_GENOTYPE_LIKEP;
                    break;
                default:
                    // TODO Internal error.
                    exit(123);
                    break;
            }
            
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
}

static inline ldoc_nde_t* vcf_proc_gt(char* val, size_t len, char* ref, char** vseqs, size_t vnum)
{
    ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_UA);
    ldoc_nde_t* nde_seq = ldoc_nde_new(LDOC_NDE_OL);
    
    nde->mkup.anno.str = (char*)VCF_GENOTYPE;
    
    nde_seq->mkup.anno.str = (char*)GEN_SEQUENCES;
    
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
    
    ent = ldoc_ent_new(LDOC_ENT_BR);
    
    ent->pld.pair.anno.str = (char*)VCF_PHASED;
    ent->pld.pair.dtm.bl = phased;
    
    ldoc_nde_ent_push(nde, ent);
    
    return nde;
}

static inline ldoc_ent_t* vcf_proc_num(char* val, size_t len, char* lbl, size_t lbllen)
{
    ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_NR);
    
    if (lbllen == 2)
    {
        if (*lbl == 'A' && *(lbl + 1) == 'N') // AN: total number of alleles
            ent->pld.pair.anno.str = qk_strdup(GEN_ALLELE_TTL);
        else if (*lbl == 'B' && *(lbl + 1) == 'Q') // BQ: base quality
            ent->pld.pair.anno.str = qk_strdup(GEN_QUALITY_RMS);
        else if (*lbl == 'D' && *(lbl + 1) == 'P') // DP: depth
            ent->pld.pair.anno.str = qk_strdup(GEN_DEPTH);
        else if (*lbl == 'G' && *(lbl + 1) == 'Q') // GQ: genotype quality
            ent->pld.pair.anno.str = qk_strdup(GEN_GENOTYPE_QUAL);
        else if (*lbl == 'M' && *(lbl + 1) == 'Q') // MQ: mapping quality
            ent->pld.pair.anno.str = qk_strdup(GEN_QUALITY_MAP_VCF);
        else if (*lbl == 'N' && *(lbl + 1) == 'S') // NS: number of samples with data
            ent->pld.pair.anno.str = qk_strdup(GEN_SAMPLES_DATA);
        else if (*lbl == 'P' && *(lbl + 1) == 'S') // PS: phase set
            ent->pld.pair.anno.str = qk_strdup(GEN_PHASE_SET);
        else if (*lbl == 'P' && *(lbl + 1) == 'Q') // PQ: phasing quality
            ent->pld.pair.anno.str = qk_strdup(GEN_PHASE_QUAL);
        else
            ent->pld.pair.anno.str = qk_strndup(lbl, lbllen);
    }
    else if (lbllen == 3)
    {
        if (*lbl == 'E' && *(lbl + 1) == 'N' && *(lbl + 1) == 'D') // END
            ent->pld.pair.anno.str = qk_strdup(GEN_END);
        else if (*lbl == 'M' && *(lbl + 1) == 'Q' && *(lbl + 1) == '0') // MQ0
            ent->pld.pair.anno.str = qk_strdup(GEN_QUALITY_MAP0);
        else
            ent->pld.pair.anno.str = qk_strndup(lbl, lbllen);
    }
    else
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
        
        // GT comes first, because it is everywhere; then alphabetical order.
        // TODO Optimization.
        //   1. could be optimized by not checking length multiple times.
        //   2. could be optimized by not checking prefixes multiple times.
        if (fmt - id == 2 && id[0] == 'G' && id[1] == 'T')
        {
            // GT: genotype
            ldoc_nde_t* gt = vcf_proc_gt(val, smpl - val, ref, vseqs, vnum);
            
            ldoc_nde_dsc_push(s, gt);
        }
        else if (fmt - id == 2 && id[0] == 'E' && id[1] == 'C')
        {
            // PL: phred-scaled genotype likelihoods
            vcf_proc_ec(s, val, smpl - val);
        }
        else if (fmt - id == 2 && id[0] == 'F' && id[1] == 'T')
        {
            // FT: filter
            vcf_proc_ft(s, val, smpl - val);
        }
        else if (fmt - id == 2 && id[0] == 'G' && id[1] == 'L')
        {
            // GL: genotype likelihoods
            vcf_proc_glgppl(s, val, smpl - val, GEN_ALT_GL);
        }
        else if (fmt - id == 2 && id[0] == 'G' && id[1] == 'P')
        {
            // GP:  phred-scaled genotype posterior probabilities
            vcf_proc_glgppl(s, val, smpl - val, GEN_ALT_GP);
        }
        else if (fmt - id == 3 && id[0] == 'G' && id[1] == 'L' && id[2] == 'E')
        {
            // GLE: genotype likelihoods
            ldoc_nde_t* gt = vcf_proc_gle(val, smpl - val, ref, vseqs, vnum);
            
            ldoc_nde_dsc_push(s, gt);
        }
        else if (fmt - id == 2 && id[0] == 'H' && id[1] == 'Q')
        {
            // HQ: haplotype qualities
            vcf_proc_hq(s, val, smpl - val, false);
        }
        else if (fmt - id == 2 && id[0] == 'P' && id[1] == 'L')
        {
            // PL: phred-scaled genotype likelihoods
            vcf_proc_glgppl(s, val, smpl - val, GEN_ALT_PL);
        }
        else if ((fmt - id == 2 &&
                  ((id[0] == 'A' && id[1] == 'N') || // AN:
                   (id[0] == 'B' && id[1] == 'Q') || // BQ:
                   (id[0] == 'D' && id[1] == 'P') || // DP:
                   (id[0] == 'G' && id[1] == 'P') || // GP: ; TODO Double check if this is not a comma separated list.
                   (id[0] == 'G' && id[1] == 'Q') || // GQ:
                   (id[0] == 'M' && id[1] == 'Q') || // MQ:
                   (id[0] == 'N' && id[1] == 'S') || // NS:
                   (id[0] == 'P' && id[1] == 'S') || // PS:
                   (id[0] == 'P' && id[1] == 'Q'))) || // PQ:
                 (fmt - id == 3 &&
                  ((id[0] == 'E' && id[1] == 'N' && id[2] == 'D') || // END
                   (id[0] == 'M' && id[1] == 'Q' && id[2] == '0')))) // MQ0
        {
            // Numeric fields with a single value.
            ldoc_ent_t* num = vcf_proc_num(val, smpl - val, id, fmt - id);
            
            ldoc_nde_ent_push(s, num);
        }
        else
        {
            const char* pth[] = { GEN_ATTRS };
            ldoc_res_t* res = ldoc_find_anno_nde(s, (char**)pth, 1);
            
            ldoc_nde_t* usr;
            if (!res)
            {
                usr = ldoc_nde_new(LDOC_NDE_UA);
                
                // TODO Error handling.
                
                usr->mkup.anno.str = (char*)GEN_ATTRS;
                
                ldoc_nde_dsc_push(s, usr);
            }
            else
                usr = res->info.nde;
            
            char* qk_val = qk_strndup(val, smpl - val);
            ldoc_ent_t* anno = ldoc_ent_new(gen_smrt_tpe(qk_val));
            
            // TODO Error handling.

            // Default: add simple key/value pair
            anno->pld.pair.anno.str = qk_strndup(id, fmt - id);
            anno->pld.pair.dtm.str = qk_val;
            
            ldoc_nde_ent_push(usr, anno);
        }
        
        if (*smpl)
            smpl++;
    } while (*(fmt++));
    
    ldoc_nde_dsc_push(smpls, s);
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
    uint8_t col = 1;
    off = 0;
    char* coff[stt->vcf_col];
    coff[0] = ln;
    for (; col <= stt->vcf_col; col++)
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
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_ent_t* ctx = ldoc_ent_new(LDOC_ENT_OR);
    ctx->pld.pair.anno.str = (char*)JSONLD_CTX;
    ctx->pld.pair.dtm.str = (char*)JSONLD_VCF_1;
    ldoc_nde_ent_push(ftr, ctx);

    ldoc_ent_t* tpe = ldoc_ent_new(LDOC_ENT_OR);
    tpe->pld.pair.anno.str = (char*)JSONLD_TPE;
    tpe->pld.pair.dtm.str = (char*)JSONLD_CLSS_FTR;
    ldoc_nde_ent_push(ftr, tpe);
    
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
    
    // Assume that all features are reported on the forward strand:
    ldoc_ent_t* strnd = ldoc_ent_new(LDOC_ENT_OR);
    strnd->pld.pair.anno.str = (char*)GEN_STRAND;
    strnd->pld.pair.dtm.str = (char*)GEN_FORWARD;

    vcf_proc_idlst(ftr, coff[2]);

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
    
    vcf_proc_optlst(ftr, (char*)GEN_ANNOTATIONS, coff[6]);
    
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
    
    // Score:
    ldoc_nde_ent_push(ftr, scr);
    
    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, lm);
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    ldoc_nde_ent_push(lc, strnd);
    
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
        
        size_t i = 9;
        for (; i < stt->vcf_col; i++)
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

// JSON to VCF

static inline void vcf_proc_doc_strct(ldoc_nde_t* meta, char* ky, char* alt)
{
    char* pth[] = { ky };
    ldoc_res_t* res = ldoc_find_anno_nde(meta, pth, 1);
    
    if (!res || !res->info.nde->dsc_cnt)
        return;
    
    ldoc_nde_t* cntnr;
    TAILQ_FOREACH(cntnr, &(res->info.nde->dscs), ldoc_nde_entries)
    {
        if (cntnr->dsc_cnt)
        {
            // TODO Format error.
        }
        
        if (!qk_heap_empty())
            qk_strcat("\n");
        
        qk_strcat("##");
        qk_strcat(alt);
        qk_strcat("=<ID=");
        qk_strcat(cntnr->mkup.anno.str);
        
        ldoc_ent_t* ent;
        TAILQ_FOREACH(ent, &(cntnr->ents), ldoc_ent_entries)
        {
            qk_strcat(",");
        
            // Figure out whether `dtm.str` contains spaces, if so,
            // then it is necessary to add quotes:
            bool sp = false;
            char* s = ent->pld.pair.dtm.str;
            while (*s)
            {
                if (*s == ' ' || *s == '\t')
                {
                    sp = true;
                    break;
                }
                
                s++;
            }
            
            qk_strcat(ent->pld.pair.anno.str);
            qk_strcat("=");
            
            if (sp)
                qk_strcat("\"");
            qk_strcat(ent->pld.pair.dtm.str);
            if (sp)
                qk_strcat("\"");
        }
        
        qk_strcat(">");
    }
}

static inline void vcf_proc_doc_meta(ldoc_nde_t* meta)
{
    gen_proc_doc_prgm_kv(meta, (char*)GEN_VCFVERSION, (char*)GEN_VCFVERSION_VCF, "=VCFv");
    gen_proc_doc_prgm_kv(meta, (char*)GEN_FILE_DATE, (char*)GEN_FILE_DATE_VCF1, "=");
    gen_proc_doc_prgm_kv(meta, (char*)GEN_FST_REF, (char*)GEN_FST_REF_VCF, "=");
    gen_proc_doc_prgm_kv(meta, "phasing", "phasing", "=");
    gen_proc_doc_prgm_kv(meta, (char*)GEN_FST_BKPNT, (char*)GEN_FST_BKPNT_VCF, "=");
    
    vcf_proc_doc_strct(meta, (char*)GEN_SEQUENCE_REGION, (char*)GEN_SEQUENCE_REGION_VCF);
    vcf_proc_doc_strct(meta, (char*)GEN_FMT_INFO, (char*)GEN_FMT_INFO_VCF);
    vcf_proc_doc_strct(meta, (char*)GEN_FMT_FLT, (char*)GEN_FMT_FLT_VCF);
    vcf_proc_doc_strct(meta, (char*)GEN_FMT_GT, (char*)GEN_FMT_GT_VCF);
    vcf_proc_doc_strct(meta, (char*)GEN_FMT_ALT, (char*)GEN_FMT_ALT_VCF);
    
    gen_proc_doc_prgm(meta, "=");
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

static bool vcf_hdr = false;
static char vcf_smpls[1024 * 1024];
static char vcf_fmt[512];

static inline void vcf_proc_doc_smpl_fmt(ldoc_nde_t* smpl, const char* ky, const char* alt, bool nde)
{
    ldoc_res_t* fmt;
    
    if (nde)
    {
        if (!alt)
        {
            const char* fld_pth[] = { ky };
            fmt = ldoc_find_anno_nde(smpl, (char**)fld_pth, 1);
            
            if (fmt)
            {
                ldoc_ent_t* usr;
                TAILQ_FOREACH(usr, &(fmt->info.nde->ents), ldoc_ent_entries)
                {
                    if (*vcf_fmt)
                        strcat(vcf_fmt, ":");
                    
                    strcat(vcf_fmt, usr->pld.pair.anno.str);
                }
                return;
            }
        }
        if (!strcmp(ky, GEN_GENOTYPE_LIKE) ||
            !strcmp(ky, GEN_GENOTYPE_LIKEP))
        {
            const char* fld_pth[] = { "AA", ky };
            fmt = ldoc_find_anno_nde(smpl, (char**)fld_pth, 2);
        }
        else
        {
            const char* fld_pth[] = { ky };
            fmt = ldoc_find_anno_nde(smpl, (char**)fld_pth, 1);
        }
    }
    else
        fmt = ldoc_find_anno_ent(smpl, (char*)ky);
    
    if (fmt)
    {
        if (*vcf_fmt)
            strcat(vcf_fmt, ":");
        
        strcat(vcf_fmt, alt);
    }
}

static inline void vcf_proc_doc_smpl_hdrfmt(ldoc_nde_t* ftr)
{
    ldoc_nde_t* smpl;
    
    const char* smpl_pth[] = { "samples" };
    ldoc_res_t* res = ldoc_find_anno_nde(ftr, (char**)smpl_pth, 1);
    
    if (!vcf_hdr)
    {
        vcf_smpls[0] = 0;

        qk_strcat("##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        
        TAILQ_FOREACH(smpl, &(res->info.nde->dscs), ldoc_nde_entries)
        {
            char* sid = smpl->mkup.anno.str;
            
            if (!*vcf_smpls)
            {
                strcat(vcf_smpls, sid);
                vcf_smpls[strlen(vcf_smpls) + 1] = 0;
            }
            else
            {
                char* s = &vcf_smpls[strlen(vcf_smpls)];
                
                s++;
                strcat(s, sid);
                s[strlen(s) + 1] = 0;
            }
            
            qk_strcat("\t");
            qk_strcat(sid);
        }
        
        qk_strcat("\n");
    }
    
    vcf_fmt[0] = 0;
    
    TAILQ_FOREACH(smpl, &(res->info.nde->dscs), ldoc_nde_entries)
    {
        vcf_proc_doc_smpl_fmt(smpl, GEN_GENOTYPE, GEN_GENOTYPE_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_DEPTH, GEN_DEPTH_VCF, false);
        vcf_proc_doc_smpl_fmt(smpl, GEN_ANNOTATIONS, GEN_ANNOTATIONS_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_GENOTYPE_LIKE, GEN_GENOTYPE_LIKE_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_GENOTYPE_LIKEP, GEN_GENOTYPE_LIKEP_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_GENOTYPE_PROB, GEN_GENOTYPE_PROB_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_GENOTYPE_QUAL, GEN_GENOTYPE_QUAL_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_HAP_QUALITIES, GEN_HAP_QUALITIES_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_PHASE_SET, GEN_PHASE_SET_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_PHASE_QUAL, GEN_PHASE_QUAL_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_ALLELE_CNTEXP, GEN_ALLELE_CNTEXP_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_QUALITY_MAP, GEN_QUALITY_MAP_VCF, true);
        vcf_proc_doc_smpl_fmt(smpl, GEN_ATTRS, NULL, true);
    }
    
    vcf_hdr = true;
}

static inline void vcf_proc_doc_gt(ldoc_nde_t* smpl)
{
    const char* gt_pth[] = { GEN_GENOTYPE };
    ldoc_res_t* fld = ldoc_find_anno_nde(smpl, (char**)gt_pth, 1);
    
    if (!fld)
    {
        qk_strcat(".");
        return;
    }
    
    ldoc_res_t* ph = ldoc_find_anno_ent(fld->info.nde, "phased");
    bool phsd = ph ? (ph->info.ent->tpe == LDOC_ENT_BR ? ph->info.ent->pld.pair.dtm.bl : false) : false;
    
    ldoc_res_t* all = ldoc_find_anno_ent(fld->info.nde, "alleles");
    
    char* astr = all->info.ent->pld.pair.dtm.str;
    char gt[2] = { 0, 0 };
    while (*astr)
    {
        if (astr > all->info.ent->pld.pair.dtm.str)
            qk_strcat(phsd ? "|" : "/");
        
        gt[0] = *astr - 'A' + '0';
        qk_strcat(gt);
        
        astr++;
    }
}

static inline void vcf_proc_doc_ent(ldoc_nde_t* smpl, const char* ky)
{
    ldoc_res_t* fld = ldoc_find_anno_ent(smpl, (char*)ky);
    
    if (!fld)
    {
        qk_strcat(".");
        return;
    }

    qk_strcat(fld->info.ent->pld.pair.dtm.str);
}

static inline void vcf_proc_doc_lst(ldoc_nde_t* smpl, const char* ky, char* sep)
{
    const char* pth[] = { ky };
    ldoc_res_t* fld = ldoc_find_anno_nde(smpl, (char**)pth, 1);
    
    if (!fld)
    {
        qk_strcat(".");
        return;
    }
    
    ldoc_ent_t* ent;
    TAILQ_FOREACH(ent, &(fld->info.nde->ents), ldoc_ent_entries)
    {
        if (ent != TAILQ_FIRST(&(fld->info.nde->ents)))
            qk_strcat(sep);
        
        qk_strcat(ent->pld.str);
    }
}

static inline void vcf_proc_doc_glgppl(ldoc_nde_t* smpl, gen_alt_t tpe)
{
    ldoc_res_t* res;
    uint8_t n = 0;
    for (; n < GEN_MAX_ALT; n++)
    {
        char* pth[] = { &GEN_ALLELES[GEN_STEP * n] };
        res = ldoc_find_anno_nde(smpl, pth, 1);
        
        if (!res)
            return;
        
        ldoc_res_t* res_ent = NULL;
        switch (tpe)
        {
            case GEN_ALT_GL:
                res_ent = ldoc_find_anno_ent(res->info.nde, (char*)GEN_GENOTYPE_LIKE);
                break;
            case GEN_ALT_GP:
                res_ent = ldoc_find_anno_ent(res->info.nde, (char*)GEN_GENOTYPE_PROB);
                break;
            case GEN_ALT_PL:
                res_ent = ldoc_find_anno_ent(res->info.nde, (char*)GEN_GENOTYPE_LIKEP);
                break;
            default:
                // TODO: Internal error
                break;
        }
        
        if (!res_ent)
            return;
        
        if (n)
            qk_strcat(",");
        
        qk_strcat(res_ent->info.ent->pld.pair.dtm.str);

    }
}

static inline bool vcf_proc_doc_smpl(ldoc_nde_t* ftr, char* attrs)
{
    const char* smpl_pth[] = { "samples" };
    ldoc_res_t* res = ldoc_find_anno_nde(ftr, (char**)smpl_pth, 1);
    
    if (!res)
        return true;

    qk_strcat("\t");
    qk_strcat(vcf_fmt);
    
    char* soff;
    char* sid = vcf_smpls;
    while (*sid)
    {
        qk_strcat("\t");
        soff = qk_working_ptr();
        
        const char* pth[] = { sid };
        ldoc_res_t* smpl_res = ldoc_find_anno_nde(res->info.nde, (char**)pth, 1);

        if (!smpl_res)
            gen_err(MAIN_ERR_JFMT_KY, sid);
        
        char* fmt = strdup(vcf_fmt);
        char* fid = fmt;
        while (*fid)
        {
            // Separate fields:
            if (fid > fmt)
                qk_strcat(":");
            
            // Terminate this format ID:
            char* fid_ = fid;
            while (*fid_ && *fid_ != ':')
                fid_++;
            *fid_ = 0;
            
            if (!strcmp(fid, GEN_GENOTYPE_VCF))
                vcf_proc_doc_gt(smpl_res->info.nde); // GT: genotype
            else if (!strcmp(fid, GEN_DEPTH_VCF))
                vcf_proc_doc_ent(smpl_res->info.nde, GEN_DEPTH); // DP: depth
            else if (!strcmp(fid, GEN_ANNOTATIONS_VCF))
                vcf_proc_doc_lst(smpl_res->info.nde, GEN_ANNOTATIONS, ";"); // FT: filter
            else if (!strcmp(fid, GEN_GENOTYPE_LIKE_VCF))
                vcf_proc_doc_glgppl(smpl_res->info.nde, GEN_ALT_GL); // GL: genotype likelihood
            else if (!strcmp(fid, GEN_GENOTYPE_PROB_VCF))
                vcf_proc_doc_glgppl(smpl_res->info.nde, GEN_ALT_GP); // GP: genotype probabilities phred scaled
            else if (!strcmp(fid, GEN_GENOTYPE_LIKEP_VCF))
                vcf_proc_doc_glgppl(smpl_res->info.nde, GEN_ALT_PL); // PL: genotype likelihood phred scaled
            else if (!strcmp(fid, GEN_GENOTYPE_QUAL_VCF))
                vcf_proc_doc_ent(smpl_res->info.nde, GEN_GENOTYPE_QUAL); // GQ: genotype quality
            else if (!strcmp(fid, GEN_HAP_QUALITIES_VCF))
                vcf_proc_doc_lst(smpl_res->info.nde, GEN_HAP_QUALITIES, ","); // HQ: haplotype qualities
            else if (!strcmp(fid, GEN_PHASE_QUAL_VCF))
                vcf_proc_doc_ent(smpl_res->info.nde, GEN_PHASE_QUAL); // PQ: phasing quality
            else if (!strcmp(fid, GEN_QUALITY_MAP_VCF))
                vcf_proc_doc_ent(smpl_res->info.nde, GEN_QUALITY_MAP); // MQ: mapping quality
            else if (!strcmp(fid, GEN_PHASE_SET_VCF))
                vcf_proc_doc_ent(smpl_res->info.nde, GEN_PHASE_SET); // PS: phase set
            else
            {
                // User defined attributes:
                const char* usr_pth[] = { GEN_ATTRS };
                ldoc_res_t* res_usr = ldoc_find_anno_nde(smpl_res->info.nde, (char**)usr_pth, 1);
                
                if (!res_usr)
                    gen_err(MAIN_ERR_JFMT_KY, GEN_ATTRS);
                
                ldoc_res_t* res_fid = ldoc_find_anno_ent(res_usr->info.nde, fid);
                
                if (!res_fid)
                    gen_err(MAIN_ERR_JFMT_KY, fid);
                
                qk_strcat(res_fid->info.ent->pld.pair.dtm.str);
                
                ldoc_res_free(res_fid);
                ldoc_res_free(res_usr);
            }
            
            // Proceed to the next format ID:
            while (*fid)
                fid++;
            fid++;
        }
        
        free(fmt);
        ldoc_res_free(smpl_res);
        
        // Proceed to the next sample ID:
        while (*sid)
            sid++;
        sid++;
    }
    
    ldoc_res_free(res);
    
    return true;
}

static inline bool vcf_proc_doc_info(ldoc_nde_t* ftr, char* attrs)
{
    /*
    ldoc_res_t* res = ldoc_find_anno_ent(ftr, (char*)GEN_ALLELE_CNT);
    if (res)
        gen_join_attrs_ent((char*)GEN_ALLELE_CNT_VCF, res->info.ent, attrs);

    res = ldoc_find_anno_ent(ftr, (char*)GEN_ALLELE_FRQ);
    if (res)
        gen_join_attrs_ent((char*)GEN_ALLELE_FRQ_VCF, res->info.ent, attrs);
    */
    
    // Alleles:
    
    ldoc_res_t* res = ldoc_find_anno_ent(ftr, (char*)GEN_ALLELE_TTL);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_ALLELE_TTL_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_ANCESTRAL_ALLELE);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_ANCESTRAL_ALLELE_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    // Memberships/truth-tags:
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_MEMBER_1000G);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_MEMBER_1000G_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_MEMBER_DBSNP);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_MEMBER_DBSNP_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_MEMBER_HM2);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_MEMBER_HM2_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_MEMBER_HM3);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_MEMBER_HM3_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_SOMATIC);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_SOMATIC_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    // Qualities:
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_QUALITY_MAP);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_QUALITY_MAP_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_QUALITY_RMS);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_QUALITY_RMS_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_QUALITY_MAP0);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_QUALITY_MAP0_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    // Samples:
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_SAMPLES_DATA);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_SAMPLES_DATA_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    // Other:
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_DEPTH);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_DEPTH_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_STRAND_BIAS);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_STRAND_BIAS_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    res = ldoc_find_anno_ent(ftr, (char*)GEN_VALIDATED);
    if (res)
    {
        gen_join_attrs_ent((char*)GEN_VALIDATED_VCF, res->info.ent, attrs);
        ldoc_res_free(res);
    }
    
    const char* algn_pth[] = { GEN_ALIGNMENT };
    res = ldoc_find_anno_nde(ftr, (char**)algn_pth, 1);
    if (res && res->nde)
    {
        res = ldoc_find_anno_ent(res->info.nde, (char*)GEN_CIGAR);
        
        if (res && !res->nde)
        {
            // Add attribute separator, if other attributes are present:
            if (qk_working_ptr() != attrs)
                qk_strcat(";");

            qk_strcat(GEN_CIGAR_VCF);
            qk_strcat("=");
            qk_strcat(res->info.ent->pld.pair.dtm.str);
        }
    }
    
    if (res)
        ldoc_res_free(res);

    vcf_proc_doc_smpl(ftr, attrs);
    
    return true;
}

inline char* vcf_proc_doc_ftr(ldoc_nde_t* ftr)
{
    vcf_proc_doc_smpl_hdrfmt(ftr);
    
    const char* lc_id[] = { GEN_LOCUS };
    ldoc_res_t* lc = ldoc_find_anno_nde(ftr, (char**)lc_id, 1);
    
    ldoc_res_t* lm = ldoc_find_anno_ent(lc->info.nde, (char*)VCF_C1);
    ldoc_res_t* st = ldoc_find_anno_ent(lc->info.nde, (char*)VCF_C2_1);
    ldoc_res_t* en = ldoc_find_anno_ent(lc->info.nde, (char*)VCF_C2_2);
    
    ldoc_res_free(lc);
    
    // TODO Check that st and en value are the same!
    if (!st || !en)
    {
        // TODO Format error.
    }
    else
    {
        if (st->nde ||
            en->nde ||
            st->info.ent->tpe != LDOC_ENT_NR ||
            en->info.ent->tpe != LDOC_ENT_NR ||
            strcmp(st->info.ent->pld.str, en->info.ent->pld.str))
        {
            gen_err(MAIN_ERR_JFMT_DEP, "Start/end coordinate mismatched or malformatted");
        }
    }
    
    const char* ref_id[] = { VCF_C4 };
    ldoc_res_t* ref = ldoc_find_anno_nde(ftr, (char**)ref_id, 1);
    ldoc_res_t* ref_seq = ldoc_find_anno_ent(ref->info.nde, (char*)GEN_SEQUENCE);
    
    qk_strcat(gen_res_req(lm));
    qk_strcat("\t");
    qk_strcat(gen_res_req(st));
    qk_strcat("\t");
    
    ldoc_res_free(lm);
    ldoc_res_free(st);
    ldoc_res_free(en);
    
    // ID and aliases:
    ldoc_res_t* ftr_id = ldoc_find_anno_ent(ftr, (char*)VCF_C3);
    if (ftr_id && !ftr_id->nde)
    {
        if (ftr_id->info.ent->pld.str)
            qk_strcat(ftr_id->info.ent->pld.str);
        
        // TODO: If there is NO id, then there SHOULD NOT be any aliases (document format error).
        const char* alias_id[] = { GEN_ALIAS };
        ldoc_res_t* als = ldoc_find_anno_nde(ftr, (char**)alias_id, 1);
        if (als)
        {
            ldoc_ent_t* als_str;
            TAILQ_FOREACH(als_str, &(als->info.nde->ents), ldoc_ent_entries)
            {
                if (ftr_id->info.ent->pld.str)
                    qk_strcat(";");
                
                qk_strcat(als_str->pld.str);
            }
        }
    }
    
    if (ftr_id)
        ldoc_res_free(ftr_id);
    
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
        {
            str = seq->info.ent->pld.pair.dtm.str;
            ldoc_res_free(seq);
        }
        
        if (allele)
            qk_strcat(",");
        
        qk_strcat(str);
        
        ldoc_res_free(var);
        
        allele++;
    }
    qk_strcat("\t");
    
    ldoc_res_free(vars);
    
    // SCORE column:
    ldoc_res_t* scr = ldoc_find_anno_ent(ftr, (char*)VCF_C6);
    qk_strcat(gen_res_opt(scr));
    qk_strcat("\t");
    
    ldoc_res_free(scr);

    // FILTER column:
    vcf_proc_doc_optlst(ftr, (char*)GEN_ANNOTATIONS, "PASS");
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
        
        ldoc_res_free(usr);
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
    
    ldoc_res_free(ref);
    ldoc_res_free(ref_seq);
    
    return NULL;
}

inline char* vcf_proc_doc_ftr_attrs(ldoc_nde_t* ftr)
{
    return NULL;
}

char* vcf_proc_doc(ldoc_doc_t* doc, gen_doctype_t tpe)
{
    char* attr;
    
    switch (tpe)
    {
        case GEN_FMT_INF:
            vcf_proc_doc_meta(doc->rt);
            
            return qk_heap_ptr();
        case GEN_FMT_FTR:
            attr = vcf_proc_doc_ftr_attrs(doc->rt);
            
            vcf_proc_doc_ftr(doc->rt);
            
            if (attr && *attr)
                qk_strcat(";");
            qk_strcat(attr);
            
            free(attr);
            
            return qk_heap_ptr();
        default:
            // TODO Internal error.
            return NULL;
    }
}

