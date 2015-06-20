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

const char* GVF_C1  = "landmark";
const char* GVF_C2  = "source";
const char* GVF_C3  = "type";
const char* GVF_C4  = "start";
const char* GVF_C5  = "end";
const char* GVF_C6  = "score";
const char* GVF_C7  = "strand";

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

inline void gvf_proc_effct(ldoc_nde_t* vars, char** val_cmp, bool* lend)
{
    // TODO This function seems long and uses many variables. Can be shortened?
    
    // Remember where the processing started.
    char* val = *val_cmp;
    
    // Format: sequence_variant index feature_type feature_ID feature_ID
    char* val_spc;
    char* spc;
    off_t eff;
    ldoc_ent_t* kv_ent;
    ldoc_res_t* nde_res;
    ldoc_res_t* nde_eff;
    ldoc_nde_t* eff_lst;
    ldoc_nde_t* aff_lst;
    ldoc_nde_t* eff_ctnr;
    bool flst_done;
    while (!*lend)
    {
        // Check: this should always work, since strlen(val) > 0; but check!?
        while (*val && *val != ',')
            val++;
        
        if (!*val)
            *lend = true;
        else
            *val = 0;
        
        // Go through values that are separated by space:
        spc = *val_cmp;
        while (*spc && *spc != ' ')
            spc++;
        *spc = 0;
        
        // Now *val_cmp denotes the sequence variant SO type.
        // Find variant index:
        val_spc = ++spc;
        while (*spc && *spc != ' ')
            spc++;
        *spc = 0;
        eff = strtoll(val_spc, NULL, 10);
        // TODO Check that eff is a number and less than dsc_cnt!
        // Find the node (can the desclaration be moved up?):
        const char* vars_id[] = { &GEN_ALLELE[(eff + 1) * 2] };
        nde_res = ldoc_find_anno_nde(vars, (char**)vars_id, 1);
        
        if (!nde_res)
        {
            // TODO Error in data format.
        }
        
        const char* effs_id[] = { GEN_EFFECTS };
        nde_eff = ldoc_find_anno_nde(nde_res->info.nde, (char**)effs_id, 1);
        
        if (!nde_eff)
        {
            eff_lst = ldoc_nde_new(LDOC_NDE_OL);
            
            // TODO Error handling.
            
            eff_lst->mkup.anno.str = (char*)GEN_EFFECTS;
            ldoc_nde_dsc_push(nde_res->info.nde, eff_lst);
        }
        else
            eff_lst = nde_eff->info.nde;
        
        eff_ctnr = ldoc_nde_new(LDOC_NDE_UA);
        
        // TODO Error handling.
        
        ldoc_nde_dsc_push(eff_lst, eff_ctnr);
        
        // First space-separated value: effect
        kv_ent = ldoc_ent_new(LDOC_ENT_OR);
        
        if (!kv_ent)
        {
            // TODO Error handling.
        }
        
        kv_ent->pld.pair.anno.str = (char*)GEN_EFFECT;
        kv_ent->pld.pair.dtm.str = *val_cmp;
        ldoc_nde_ent_push(eff_ctnr, kv_ent);
        
        // Third space-separated value: affected feature type
        val_spc = ++spc;
        while (*spc && *spc != ' ')
            spc++;
        *spc = 0;
        
        kv_ent = ldoc_ent_new(LDOC_ENT_OR);
        
        if (!kv_ent)
        {
            // TODO Error handling.
        }
        
        kv_ent->pld.pair.anno.str = (char*)GEN_AFFECTED_TPE;
        kv_ent->pld.pair.dtm.str = val_spc;
        ldoc_nde_ent_push(eff_ctnr, kv_ent);
        
        // Remainder: feature IDs
        
        aff_lst = ldoc_nde_new(LDOC_NDE_OL);
        
        // TODO Error handling.
        
        aff_lst->mkup.anno.str = (char*)GEN_AFFECTED;
        ldoc_nde_dsc_push(eff_ctnr, aff_lst);
        
        flst_done = false;
        while (!flst_done)
        {
            val_spc = ++spc;
            while (*spc && *spc != ' ')
                spc++;
            if (!*spc)
                flst_done = true;
            else
                *spc = 0;
            
            kv_ent = ldoc_ent_new(LDOC_ENT_TXT);
            
            if (!kv_ent)
            {
                // TODO Error handling.
            }
            
            //kv_ent->pld.pair.anno.str = attr;
            //kv_ent->pld.pair.dtm.str = val_spc;
            kv_ent->pld.str = val_spc;
            
            ldoc_nde_ent_push(aff_lst, kv_ent);
        }
        
        *val_cmp = ++val;
    }
}

inline ldoc_doc_t* gvf_proc_ftr(int fd, off_t mx, ldoc_trie_t* idx, char* ln, size_t lnlen, char** cmt)
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
        c->pld.pair.dtm.str = gen_escstr(*cmt, GEN_FMT_GVF);
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
                    if (*val != '.')
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
    gen_splt_attrs(ftr, attrs, ref, vars, coff[8], BI_VAL);
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_nde_ent_push(ftr, ctx);
    
    // Source, type, score, strand:
    ldoc_nde_ent_push(ftr, src);
    ldoc_nde_ent_push(ftr, tpe);
    ldoc_nde_ent_push(ftr, scr);
    ldoc_nde_ent_push(ftr, strnd);
    
    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, lm);
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    
    ldoc_nde_dsc_push(ftr, lc);

    // Reference & variants:
    gen_nde_dsc_opt(ftr, ref, (char*)GEN_REFERENCE);
    gen_nde_dsc_opt(ftr, vars, (char*)GEN_VARIANTS);
    
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
            cent->pld.pair.dtm.str = gen_escstr(*cmt, GEN_FMT_GVF);
        }
        else
        {
            cent->pld.pair.anno.str = strdup(val);
            cent->pld.pair.dtm.str = gen_escstr(*cmt, GEN_FMT_GVF);
            
        }
        ldoc_nde_ent_push(cpgm, cent);
        
        // Erase comment, but it is not released (free'd) yet:
        *cmt = NULL;
    }
    
    return doc;
}

ldoc_doc_t* gvf_proc_ln(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt, gen_fstat* stat)
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
                    stat->meta++;
                }
            }
            else
            {
                // Comment:
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
        else if (st->fa_sct)
        {
            
        }
        else
        {
            // Feature:
            ldoc = gvf_proc_ftr(fd, mx, idx, ln, lnlen, cmt);
            stat->ftrs++;
        }
    }
    
    return ldoc;
}

void gvf_proc_attr_effct(ldoc_nde_t* effs, char* allele, char* astr)
{
    // Format:
    //   "effects":[
    //     {
    //       "effect":"downstream_gene_variant",
    //       "affected-feature-type":"transcript",
    //       "affected-features":[
    //                             "YAL056W"
    //                           ]
    //     },
    //   ...
    ldoc_nde_t* eff;
    ldoc_nde_t* ftrs;
    ldoc_ent_t* ftr;
    ldoc_res_t* ent;
    bool fst = true;
    while (effs->dsc_cnt)
    {
        if (fst)
            fst = false;
        else
            strcat(astr, ",");
        
        eff = effs->dscs.tqh_first;
        
        ent = ldoc_find_anno_ent(eff, "effect");
        // TODO Error handling. Data format error.
        
        strcat(astr, ent->info.ent->pld.pair.dtm.str);
        strcat(astr, " ");
        
        // Index.
        char* e = &astr[strlen(astr)];
        sprintf(e, "%u ", allele[0] - 'B');
        
        ent = ldoc_find_anno_ent(eff, "affected-feature-type");
        // TODO Error handling. Data format error.
        
        strcat(astr, ent->info.ent->pld.pair.dtm.str);
        
        // TODO There should be only one descendant.
        ftrs = eff->dscs.tqh_first;
        
        while (ftrs->ent_cnt)
        {
            strcat(astr, " ");

            ftr = ftrs->ents.tqh_first;
            
            strcat(astr, ftr->pld.str);
            
            ldoc_ent_rm(ftr);
            ldoc_ent_free(ftr);
        }
        
        ldoc_nde_rm(eff);
        ldoc_nde_free(eff);
    }
}

char* gvf_proc_doc_ftr_attrs(ldoc_nde_t* ftr)
{
    // Everything about the reference sequence:
    const char* ref_id[] = { GEN_REFERENCE };
    ldoc_res_t* ref = ldoc_find_anno_nde(ftr, (char**)ref_id, 1);
    
    // All variant alleles (B, C, D, ...):
    const char* vars_id[] = { GEN_VARIANTS };
    ldoc_res_t* vars = ldoc_find_anno_nde(ftr, (char**)vars_id, 1);
    
    // TODO Replace with some quick-heap implementation!
    char* astr = (char*)malloc(1*1024*1024);
    *astr = 0;

    bool fst = true;
    
    //
    // Reference:
    //
    
    ldoc_res_t* ref_seq = ldoc_find_anno_ent(ref->info.nde, (char*)GEN_SEQUENCE);
    strcat(astr, "Reference_sequence=");
    strcat(astr, ref_seq->info.ent->pld.pair.dtm.str);
    fst = false;
    
    const char* cdn_id[] = { GEN_CODON };
    ldoc_res_t* cdn = ldoc_find_anno_nde(ref->info.nde, (char**)cdn_id, 1);
    if (cdn)
    {
        if (!fst)
            strcat(astr, ";");
        strcat(astr, "Reference_codon=");
        strcat(astr, cdn->info.ent->pld.pair.dtm.str);
    }
    
    //
    // Variants:
    //
    
    size_t vnum = vars->info.nde->dsc_cnt; // Number of alleles ('B', 'C', etc.)
    
    if (vnum > 26)
    {
        // TODO Data error. Not supported.
    }
    
    // TODO Is this iteration right? Necessary to find 'B', 'C', etc. in order?
    // Note 1: this assumes that all allele nodes contain the
    //         the same named entities.
    // Note 2: alleles must be labeled 'B', 'C', etc.
    char* ent_nme;
    ldoc_ent_t* ent;
    while (vars->info.nde->dscs.tqh_first->ent_cnt)
    {
        if (fst)
            fst = false;
        else
            strcat(astr, ";");
        
        ent = vars->info.nde->dscs.tqh_first->ents.tqh_first;
        ent_nme = ent->pld.pair.anno.str;
        
        gen_proc_nde(vars->info.nde, ent_nme, "Variant_", astr, vnum);
    }
    
    if (fst)
        fst = false;
    else
        strcat(astr, ";");
    
    strcat(astr, "Variant_effect=");
    
    // Note: this assumes that the only descendant (if any) is "effects":
    char* nde_nme;
    ldoc_nde_t* allele;
    ldoc_nde_t* nde;
    TAILQ_FOREACH(allele, &(vars->info.nde->dscs), ldoc_nde_entries)
    {
        while (allele->dsc_cnt)
        {
            nde = allele->dscs.tqh_first;
            nde_nme = nde->mkup.anno.str;
            
            // TODO nde_nme will have to be "effects" here.
            
            gvf_proc_attr_effct(nde, allele->mkup.anno.str, astr);
            
            ldoc_nde_rm(nde);
            ldoc_nde_free(nde);
        }
    }

    return astr;
}

char* gvf_proc_doc(ldoc_doc_t* doc)
{
    char* attr = gvf_proc_doc_ftr_attrs(doc->rt);
    
    gff_proc_doc_ftr(doc->rt);
    
    if (*attr)
        qk_strcat(";");
    qk_strcat(attr);
    
    // TODO Obsolete?
    // gen_proc_doc_usr(doc->rt);
    
    free(attr);
}
