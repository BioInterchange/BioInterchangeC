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

#include "gvf.h"

static const char* GVF_C1  = "seqid";
static const char* GVF_C2  = "source";
static const char* GVF_C3  = "type";
static const char* GVF_C4  = "start";
static const char* GVF_C5  = "end";
static const char* GVF_C6  = "score";
static const char* GVF_C7  = "strand";

void gvf_cbcks(gen_cbcks_t* cbcks)
{
    cbcks->proc_ln = &gvf_proc_ln;
}

static inline ldoc_struct_t gvf_prgm_tpe(char* ky)
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
    
    // If not a GVF pragma, it might be a GFF3 pragma:
    return gff_prgm_tpe(ky_);
}

static inline ldoc_doc_t* gvf_proc_ftr(int fd, off_t mx, ldoc_trie_t* idx, char* ln, size_t lnlen, char** cmt)
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
    off = 0;
    char* coff[9];
    coff[0] = ln;
    for (uint8_t col = 1; col < 9; col++)
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
        if (coff[col - 1][0] == '.' && coff[col - 1][1] == 0)
            coff[col - 1] = NULL;
    }
    
    // See if landmark has a sequence:
    char* seq;
    if (idx && coff[3] && coff[4])
    {
        bool rv = coff[6][0] == '+' ? true : false;
        uint64_t fst = strtoull(coff[3], NULL, 10); // TODO Add second parameter to figure out errors.
        uint64_t fen = strtoull(coff[4], NULL, 10); // TODO Ditto.
        seq = gff_seq(fd, mx, idx, coff[0], fst, fen, rv);
    }
    else
    {
        seq = NULL;
    }
    
    ldoc_nde_t* ftr = doc->rt;
    
    ldoc_ent_t* ctx = ldoc_ent_new(LDOC_ENT_OR);
    ctx->pld.pair.anno.str = (char*)JSONLD_CTX;
    ctx->pld.pair.dtm.str = (char*)JSONLD_GVF;
    
    ldoc_nde_t* attrs = ldoc_nde_new(LDOC_NDE_UA);
    attrs->mkup.anno.str = (char*)GEN_ATTRS;
    
    ldoc_ent_t* lm = ldoc_ent_new(LDOC_ENT_OR);
    lm->pld.pair.anno.str = (char*)GVF_C1;
    lm->pld.pair.dtm.str = coff[0];
    
    ldoc_nde_t* lc = ldoc_nde_new(LDOC_NDE_UA);
    lc->mkup.anno.str = (char*)GEN_LOCUS;
    
    ldoc_ent_t* src = ldoc_ent_new(LDOC_ENT_OR);
    src->pld.pair.anno.str = (char*)GVF_C2;
    src->pld.pair.dtm.str = coff[1];
    
    ldoc_ent_t* tpe = ldoc_ent_new(LDOC_ENT_OR);
    tpe->pld.pair.anno.str = (char*)GVF_C3;
    tpe->pld.pair.dtm.str = coff[2];
    
    ldoc_ent_t* st = ldoc_ent_new(LDOC_ENT_NR);
    st->pld.pair.anno.str = (char*)GVF_C4;
    st->pld.pair.dtm.str = coff[3];
    
    ldoc_ent_t* en = ldoc_ent_new(LDOC_ENT_NR);
    en->pld.pair.anno.str = (char*)GVF_C5;
    en->pld.pair.dtm.str = coff[4];
    
    ldoc_ent_t* scr = ldoc_ent_new(LDOC_ENT_OR);
    scr->pld.pair.anno.str = (char*)GVF_C6;
    scr->pld.pair.dtm.str = coff[5];
    
    ldoc_ent_t* strnd = ldoc_ent_new(LDOC_ENT_OR);
    strnd->pld.pair.anno.str = (char*)GVF_C7;
    strnd->pld.pair.dtm.str = coff[6];
    
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
    
    // Add sequence information -- if available:
    if (seq)
    {
        ldoc_ent_t* sq = ldoc_ent_new(LDOC_ENT_OR);
        sq->pld.pair.anno.str = (char*)GEN_SEQUENCE;
        sq->pld.pair.dtm.str = seq;
        ldoc_nde_ent_push(ftr, sq);
    }
    
    // Find Reference_seq and Variant_seq and create respective nodes.
    // Note 1: Without X_seq, other X_* entries will be ignored.
    // Note 2: Assume that no more than VCF_MAX_ALT variants are observed.
    ldoc_nde_t* ref = NULL;
    size_t vnum = 0;
    char* vseqs[GEN_MAX_ALT];
    ldoc_nde_t* vars = NULL;
    char* xseq = coff[8];
    gvf_prs_seq prs = GVF_PRS_NONE;
    while (xseq - ln < lnlen)
    {
        if (*xseq == 'R' && *(xseq + 1) == 'e' && !strncmp(xseq + 2, "ference_seq=", 12))
        {
            xseq += 14;
            prs = GVF_PRS_REF;
        }
        else if (*xseq == 'V' && *(xseq + 1) == 'a' && !strncmp(xseq + 2, "riant_seq=", 10))
        {
            xseq += 12;
            prs = GVF_PRS_VAR;
        }
        
        if (prs)
        {
            char *xseq_ = xseq;
            while (xseq_ - ln < lnlen && *xseq_ != ';')
                xseq_++;
            
            char* val = qk_strndup(xseq, xseq_ - xseq);
            
            switch (prs)
            {
                case GVF_PRS_REF:
                    ref = ldoc_nde_new(LDOC_NDE_UA);
                    ref->mkup.anno.str = (char*)GEN_REFERENCE;
                    
                    ldoc_ent_t* ref_seq = ldoc_ent_new(LDOC_ENT_OR);
                    ref_seq->pld.pair.anno.str = (char*)GEN_SEQUENCE;
                    ref_seq->pld.pair.dtm.str = val;
                    ldoc_nde_ent_push(ref, ref_seq);
                    
                    break;
                case GVF_PRS_VAR:
                    vars = gen_variants(val, ',', vseqs, &vnum);
                    
                    break;
                default:
                    // TODO Internal error.
                    break;
            }
            
            prs = GVF_PRS_NONE;
        }
        
        
        if (ref && vars)
            break;
        
        xseq++;
    }
    
    // Generic implementation for parsing attributes:
    gen_splt_attrs(ftr, attrs, ref, vars, coff[8]);
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_nde_ent_push(ftr, ctx);
    
    // Source, type, score, strand:
    ldoc_nde_ent_push(ftr, src);
    ldoc_nde_ent_push(ftr, tpe);
    ldoc_nde_ent_push(ftr, scr);
    ldoc_nde_ent_push(ftr, strnd);
    
    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    
    ldoc_nde_ent_push(ftr, lm);
    ldoc_nde_dsc_push(ftr, lc);

    // Reference & variants:
    ldoc_nde_dsc_push(ftr, ref);
    ldoc_nde_dsc_push(ftr, vars);
    
    // Do not add user defined sub-tree if it is empty:
    if (attrs->dsc_cnt || attrs->ent_cnt)
        ldoc_nde_dsc_push(ftr, attrs);
    
    return doc;
}


static inline ldoc_doc_t* gvf_proc_prgm(ldoc_doc_t* doc, char* ln, size_t lnlen, char** cmt)
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
           *val != ' ' &&
           *val != '\t')
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
    
    ldoc_struct_t tpe = gvf_prgm_tpe(ln);
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

void gvf_proc_ln(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt)
{
    ldoc_doc_t* ldoc = NULL;
    
    if (lnlen > 3)
    {
        if (ln[0] == '#')
        {
            if (ln[1] == '#')
            {
                // Meta line (pragma statement):
                
                // Does not catch if FASTA section starts at the end of the file
                // without FASTA lines present -- but that case does not matter
                // anyway.
                if (lnlen > GFF_FA_PFX_LEN &&
                    !strncmp(ln, GFF_FA_PFX, GFF_FA_PFX_LEN) &&
                    (ln[GFF_FA_PFX_LEN] == '\n' ||
                     ln[GFF_FA_PFX_LEN] == '\r' ||
                     ln[GFF_FA_PFX_LEN] == '\t' ||
                     ln[GFF_FA_PFX_LEN] == ' '))
                {
                    st->fa_sct = true;
                }
                else
                {
                    // Nope, real meta line:
                    gvf_proc_prgm(fdoc, ln, lnlen, cmt);
                }
            }
            else
            {
                // Comment:
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
        else if (st->fa_sct)
        {
            
        }
        else
        {
            // Feature:
            ldoc = gvf_proc_ftr(fd, mx, idx, ln, lnlen, cmt);
        }
    }
    
    if (ldoc)
    {
        ldoc_ser_t* ser = ldoc_format(ldoc, json_vis_nde, json_vis_ent);
        
        printf("%s\n", ser->pld.str);
    }
}
