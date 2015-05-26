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
static const char* VCF_C7   = "filter";

static const char* VCF_SAMPLES = "samples";

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

static inline ldoc_struct_t vcf_prgm_tpe(char* ky)
{
    char* ky_ = ky;
    
    if (*ky == 'g')
    {
        ky++;
        if (*ky == 'e')
        {
            ky++;
            if (*ky == 'n')
            {
                ky++;
                if (!strcmp(ky, "omic-source"))
                    return LDOC_NDE_UA;
            }
        } else if (*ky == 'v')
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
    
    return LDOC_NDE_OL;
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
        if (!strcmp(ln, "sequence-region"))
        {
            ldoc_nde_t* nde = gff_proc_sregion(val, &id);
            
            ldoc_nde_dsc_push(stmt, nde);
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
            cent->pld.pair.dtm.str = gen_escstr(*cmt);
        }
        else
        {
            cent->pld.pair.anno.str = strdup(val);
            cent->pld.pair.dtm.str = gen_escstr(*cmt);
            
        }
        ldoc_nde_ent_push(cpgm, cent);
        
        // Erase comment, but it is not released (free'd) yet:
        *cmt = NULL;
    }
    
    return doc;
}

static inline void vcf_proc_smpl(ldoc_nde_t* smpls, gen_prsr_t* stt, size_t i, char* fmt, char* smpl)
{
    // Note: Going through fmt character-by-character, because the names
    //       are so short that it might be as efficient as looking up a
    //       fixed format that *needs* to be the same for the whole file.
    
    // TODO Assumes that fmt is not empty.
    
    ldoc_nde_t* s = ldoc_nde_new(LDOC_NDE_UA);
    
    s->mkup.anno.str = strdup(stt->vcf_hdr_off[i]);
    
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
        
        ldoc_ent_t* anno = ldoc_ent_new(LDOC_ENT_OR);
        
        // TODO Error handling.
        
        anno->pld.pair.anno.str = strndup(id, fmt - id);
        anno->pld.pair.dtm.str = strndup(val, smpl - val);
        
        ldoc_nde_ent_push(s, anno);
        
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
    
    ldoc_ent_t* lm = ldoc_ent_new(LDOC_ENT_OR);
    lm->pld.pair.anno.str = (char*)VCF_C1;
    lm->pld.pair.dtm.str = coff[0];
    
    ldoc_nde_t* lc = ldoc_nde_new(LDOC_NDE_UA);
    lc->mkup.anno.str = (char*)GEN_LOCUS;
    
    ldoc_ent_t* st = ldoc_ent_new(LDOC_ENT_OR);
    st->pld.pair.anno.str = (char*)VCF_C2_1;
    st->pld.pair.dtm.str = coff[1];
    
    ldoc_ent_t* en = ldoc_ent_new(LDOC_ENT_OR);
    en->pld.pair.anno.str = (char*)VCF_C2_2;
    en->pld.pair.dtm.str = coff[1];
    
    ldoc_ent_t* id = ldoc_ent_new(LDOC_ENT_NR);
    id->pld.pair.anno.str = (char*)VCF_C3;
    id->pld.pair.dtm.str = coff[2];
    
    ldoc_ent_t* ref = ldoc_ent_new(LDOC_ENT_NR);
    ref->pld.pair.anno.str = (char*)VCF_C4;
    ref->pld.pair.dtm.str = coff[3];
    
    ldoc_nde_t* vars = gen_variants(coff[4], ',');
    
    ldoc_ent_t* scr = ldoc_ent_new(LDOC_ENT_OR);
    scr->pld.pair.anno.str = (char*)VCF_C6;
    scr->pld.pair.dtm.str = coff[5];
    
    ldoc_ent_t* fltr = ldoc_ent_new(LDOC_ENT_OR);
    fltr->pld.pair.anno.str = (char*)VCF_C7;
    fltr->pld.pair.dtm.str = coff[6];
    
    // Add comment lines -- if available:
    if (*cmt)
    {
        ldoc_ent_t* c = ldoc_ent_new(LDOC_ENT_OR);
        c->pld.pair.anno.str = (char*)GEN_COMMENT;
        c->pld.pair.dtm.str = gen_escstr(*cmt);
        ldoc_nde_ent_push(ftr, c);
        
        // Erase comment, but it is not released (free'd) yet:
        *cmt = NULL;
    }
    
    ldoc_nde_t* attrs = ldoc_nde_new(LDOC_NDE_UA);
    attrs->mkup.anno.str = (char*)GEN_ATTRS;
    
    gen_splt_attrs(ftr, attrs, vars, coff[7]);
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_nde_ent_push(ftr, ctx);
    
    // Score:
    ldoc_nde_ent_push(ftr, scr);
    
    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    
    ldoc_nde_ent_push(ftr, lm);
    ldoc_nde_dsc_push(ftr, lc);
    
    // Reference & variants:
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
                vcf_proc_smpl(smpls, stt, i - 9, coff[8], coff[i]);
            else
            {
                ldoc_ent_t* smpl_null = ldoc_ent_new(LDOC_ENT_OR);
                
                // TODO Error handling.
                
                smpl_null->pld.pair.anno.str = strdup(stt->vcf_hdr_off[i - 9]);
                smpl_null->pld.pair.dtm.str = NULL;
                
                ldoc_nde_ent_push(smpls, smpl_null);
            }
        }
        
        ldoc_nde_dsc_push(ftr, smpls);
    }
    
    return doc;
}

void vcf_proc_ln(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt)
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
        }
        else
        {
            // TODO Error handling. Corrupt VCF file.
        }
    }
    
    if (ldoc)
    {
        ldoc_ser_t* ser = ldoc_format(ldoc, json_vis_nde, json_vis_ent);
        
        printf("%s\n", ser->sclr.str);
    }
}
