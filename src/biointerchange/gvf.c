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

#define GVF_SPGM_MAXLBL 1024

// Single-thread buffer for structured pragma lookups in `pgm_idx`:
static char* pgm_idxstr;

// Counter for structured pragma IDs (when not falling under "global"):
static uint32_t spgm_cnt = 0;

// Indices for structured pragma filters:
ldoc_trie_t* pgm_idx;

void gvf_init()
{
    pgm_idx = ldoc_trie_new();
    
    pgm_idxstr = (char*)malloc(GVF_SPGM_MAXLBL);
}

void gvf_cbcks(gen_cbcks_t* cbcks)
{
    cbcks->proc_ln = &gvf_proc_ln;
}

static inline void gvf_idxstr(char* seqid, char* type, char* source, char** str)
{
    // Note: index string needs to match with gvf_idx implementation!
    
    // Reset index string:
    **str = 0;
    
    // Do not check string lengths here. Implement boundary check
    // when reading feature lines!
    uint8_t i = 0;
    char* pth[3] = { seqid, type, source };
    for (; i < 3; i++)
        if (!pth[i] || !*pth[i])
            strcat(*str, ".");
        else
        {
            strcat(*str, " ");
            strcat(*str, pth[i]);
        }
}

static inline void gvf_idx(ldoc_nde_t* seqid, ldoc_nde_t* type, ldoc_nde_t* source, char* ref, uint8_t lvl, char* lbl, ldoc_ent_t* pld)
{
    // Note: index string needs to match with gvf_idxstr implementation!
    
    ldoc_nde_t* nde;
    
    switch (lvl)
    {
        case 0:
            nde = seqid;
            
            lbl = (char*)malloc(GVF_SPGM_MAXLBL);
            
            // TODO Error handling.
            break;
        case 1:
            nde = type;
            break;
        case 2:
            nde = source;
            break;
        default:
            nde = NULL;
            break;
    }
    
    if (lvl == 3)
    {
        // Check if an entry already exists:
        ldoc_trie_nde_t* exst = ldoc_trie_lookup(pgm_idx, lbl, false);
        
        // If a trie entry exists already, then copy over the "payload"
        // (`pld` ldoc_ent_t-contents) and free the memory of the `pld`
        // argument.
        if (exst)
        {
            ldoc_nde_ent_push((ldoc_nde_t*)exst->anno.pld, pld);
            
            return;
        }
        
        // No entry yet: create a new one (ldoc_nde_t that contains
        // the `pld` ldoc_ent_t-type payload):
        
        // Type does not matter -- only the *_ent_* functions are
        // used on this node to iterate through/append information:
        ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_OL);
        
        // TODO Error handling.
        
        ldoc_nde_ent_push(nde, pld);
        
        ldoc_trie_anno_t anno = { 0, nde };
        ldoc_trie_add(pgm_idx, lbl, ASCII, anno);
        
        return;
    }
    
    if (nde)
    {
        size_t lbllen = strlen(lbl);
        
        ldoc_ent_t* ent;
        TAILQ_FOREACH(ent, &(nde->ents), ldoc_ent_entries)
        {
            if (GVF_SPGM_MAXLBL - (long)lbllen - (long)strlen(ent->pld.str) - 2 < 0)
            {
                // TODO Error handling; data labels too long.
            }
            
            // Only add a space for value separation if the value itself
            // is not a wildcard ('.'):
            // Note: Do not need to check for strlen, because a dot would be
            //       a formatting error in GVF.
            if (*ent->pld.str != '.')
                strncat(lbl, " ", GVF_SPGM_MAXLBL - lbllen - 2);
            strncat(lbl, ent->pld.str, GVF_SPGM_MAXLBL - lbllen - strlen(ent->pld.str) - 2);
            
            gvf_idx(seqid, type, source, ref, lvl + 1, lbl, pld);
            
            // Reset label to its original form:
            lbl[lbllen] = 0;
        }
    }
    else
    {
        if (GVF_SPGM_MAXLBL - (long)strlen(lbl) - 2 < 0)
        {
            // TODO Error handling; data labels too long.
        }
        
        strncat(lbl, ".", GVF_SPGM_MAXLBL - strlen(lbl) - 2);
        
        gvf_idx(seqid, type, source, ref, lvl + 1, lbl, pld);
    }
    
}

static inline ldoc_nde_t* gvf_idx_dsc(ldoc_nde_t* prnt, ldoc_res_t* idx, char* mkup)
{
    ldoc_nde_t* dsc;
    
    if (idx)
    {
        dsc = idx->info.nde;
        
        ldoc_res_free(idx);
        
        return dsc;
    }
    
    dsc = ldoc_nde_new(LDOC_NDE_UA);
    
    // TODO Error handling.
    
    dsc->mkup.anno.str = strdup(mkup);
    
    // TODO Error handling.
    
    ldoc_nde_dsc_push(prnt, dsc);
    
    return dsc;
}

static inline ldoc_struct_t gvf_prgm_tpe(char* ky)
{
    char* ky_ = ky;
    
    if (*ky == 'a')
    {
        ky++;
        if (*ky == 't')
        {
            ky++;
            if (*ky == 't')
            {
                ky++;
                if (!strcmp(ky, "ribute-method"))
                    return LDOC_NDE_OO; // attribute-method
            }
        }
    }
    else if (*ky == 'd')
    {
        ky++;
        if (*ky == 'a')
        {
            ky++;
            if (*ky == 't')
            {
                ky++;
                if (!strcmp(ky, "a-source"))
                    return LDOC_NDE_OO; // data-source
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
                if (!strcmp(ky, "e-date"))
                    return LDOC_NDE_UA; // file-date
                else if (!strcmp(ky, "e-version"))
                    return LDOC_NDE_UA; // file-version
            }
        }
    }
    else if (*ky == 'g')
    {
        ky++;
        if (*ky == 'e')
        {
            ky++;
            if (*ky == 'n')
            {
                ky++;
                if (!strcmp(ky, "omic-source"))
                    return LDOC_NDE_UA; // genomic-source
            }
        } else if (*ky == 'v')
        {
            ky++;
            if (*ky == 'f')
            {
                ky++;
                if (!strcmp(ky, "-version"))
                    return LDOC_NDE_UA; // gvf-version
            }
        }
    }
    else if (*ky == 'i')
    {
        ky++;
        if (*ky == 'n')
        {
            ky++;
            if (*ky == 'd')
            {
                ky++;
                if (!strcmp(ky, "ividual-id"))
                    return LDOC_NDE_OL; // individual-id
            }
        }
    }
    else if (*ky == 'm')
    {
        ky++;
        if (*ky == 'u')
        {
            ky++;
            if (*ky == 'l')
            {
                ky++;
                if (!strcmp(ky, "ti-individual"))
                    return LDOC_NDE_OL; // multi-individual
            }
        }
    }
    else if (*ky == 'p')
    {
        ky++;
        if (*ky == 'h')
        {
            ky++;
            if (*ky == 'a')
            {
                ky++;
                if (!strcmp(ky, "sed-genotypes"))
                    return LDOC_NDE_OO; // phased-genotypes
            }
            else if (*ky == 'e')
            {
                ky++;
                if (!strcmp(ky, "notype-description"))
                    return LDOC_NDE_OO; // phenotype-description
            }
        }
        else if (*ky == 'o')
        {
            ky++;
            if (*ky == 'p')
            {
                ky++;
                if (!strcmp(ky, "ulation"))
                    return LDOC_NDE_UA; // population
            }
        }
    }
    else if (*ky == 's')
    {
        ky++;
        if (*ky == 'c')
        {
            ky++;
            if (*ky == 'o')
            {
                ky++;
                if (!strcmp(ky, "re-method"))
                    return LDOC_NDE_OO; // score-method
            }
        }
        else if (*ky == 'e')
        {
            ky++;
            if (*ky == 'x')
            {
                ky++;
                if (!*ky)
                    return LDOC_NDE_UA; // sex
            }
        }
        else if (*ky == 'o')
        {
            ky++;
            if (*ky == 'u')
            {
                ky++;
                if (!strcmp(ky, "rce-method"))
                    return LDOC_NDE_OO; // source-method
            }
        }
    }
    else if (*ky == 't')
    {
        ky++;
        if (*ky == 'e')
        {
            ky++;
            if (*ky == 'c')
            {
                ky++;
                if (!strncmp("hnology-", ky, 8))
                    return LDOC_NDE_OO; // technology-*
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
    uint8_t col = 1;
    off = 0;
    char* coff[9];
    coff[0] = ln;
    for (; col < 9; col++)
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
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_ent_t* ctx = ldoc_ent_new(LDOC_ENT_OR);
    ctx->pld.pair.anno.str = (char*)JSONLD_CTX;
    ctx->pld.pair.dtm.str = (char*)JSONLD_GVF_1;
    ldoc_nde_ent_push(ftr, ctx);
    
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
    
    // Add structured pragma links -- if available:
    // TODO Optimize, so that not all 2^3-1 choices need to be probed:
    uint8_t i = 1;
    for (; i < 8; i++)
    {
        gvf_idxstr(i & 4 ? coff[0] : NULL,
                   i & 2 ? coff[2] : NULL,
                   i & 1 ? coff[1] : NULL,
                   &pgm_idxstr);
        ldoc_trie_nde_t* lnk = ldoc_trie_lookup(pgm_idx, pgm_idxstr, false);
        if (lnk)
        {
            ldoc_nde_t* lnk_idx = (ldoc_nde_t*)lnk->anno.pld;
            ldoc_nde_t* lnk_anno = ldoc_nde_new(LDOC_NDE_OL);
            
            // TODO Error handling.
            
            lnk_anno->mkup.anno.str = (char*)GEN_ANNOTATIONS;
            
            ldoc_ent_t* eidx;
            ldoc_ent_t* ecpy;
            TAILQ_FOREACH(eidx, &(lnk_idx->ents), ldoc_ent_entries)
            {
                ecpy = ldoc_ent_new(eidx->tpe);
                
                // TODO Error handling.
                
                // Note: assumes a type that is accessible as plain string:
                ecpy->pld.str = eidx->pld.str;
                
                ldoc_nde_ent_push(lnk_anno, ecpy);
            }
            
            ldoc_nde_dsc_push(ftr, lnk_anno);
        }
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
    
    // Source, type, score, strand:
    ldoc_nde_ent_push(ftr, src);
    ldoc_nde_ent_push(ftr, tpe);
    ldoc_nde_ent_push(ftr, scr);
    
    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, lm);
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    ldoc_nde_ent_push(lc, strnd);
    
    ldoc_nde_dsc_push(ftr, lc);

    // Reference & variants:
    gen_nde_dsc_opt(ftr, ref, (char*)GEN_REFERENCE);
    gen_nde_dsc_opt(ftr, vars, (char*)GEN_VARIANTS);
    
    // Do not add user defined sub-tree if it is empty:
    if (attrs->dsc_cnt || attrs->ent_cnt)
        ldoc_nde_dsc_push(ftr, attrs);
    
    return doc;
}

inline void gvf_tags(ldoc_nde_t* nde, char* strct, gvf_pgm_t tpe)
{
    // Structure:
    //   Abstract: nde -> sub -> ent
    //   Example: "technology" -> "global" -> "read-length" : 32
    
    ldoc_nde_t* sub = ldoc_nde_new(LDOC_NDE_UA);
    
    // TODO Error handling.
    
    // 32-bit number: 10 digits; plus two letters 'SP'.
    char* uid = (char*)malloc(10 + 2 + 1);
    
    // TODO Error handling.
    
    sprintf(uid, "SP%u", spgm_cnt++);

    sub->mkup.anno.str = uid;
    
    // Supported tags for all GVF structured pragmas: Seqid, Source, Type, Dbxref, Comment
    bool glbl = true;
    bool brk = false;
    char* val;
    char* ky = strct;
    ldoc_nde_t* seqid_nde = NULL;
    ldoc_nde_t* source_nde = NULL;
    ldoc_nde_t* type_nde = NULL;
    do
    {
        if (!*strct || *strct == ';')
        {
            // Note: similar to gen_splt_attrs!
            // If this is the end of the string, then make sure to bail out next:
            if (!*strct)
                brk = true;
            
            // Isolate part of the underlying string:
            *strct = 0;
            
            // Isolate key/value:
            gen_ky(ky, &val);
            
            // Nothing to handle -- break condition 1:
            if (!*ky)
                break;
            
            // Check for key/value pairs that are comma separated:
            if (!strcmp(ky, GEN_DBXREF_GFF3) ||
                !strcmp(ky, GEN_READ_PAIR_SPAN_GVF))
            {
                gen_lwrhyph(ky);
                
                gen_attr_t kwd = { BI_NKW, (char*)NULL };
                gen_csep_dup(sub, kwd, ky, val, true);
            }
            else
            {
                // Determine data type; default: LDOC_ENT_OR
                ldoc_content_t tpe;
                if (!strcmp(ky, "Read_length"))
                    tpe = LDOC_ENT_NR;
                else if (!strcmp(ky, GEN_AVG_COVERAGE_GVF))
                    tpe = LDOC_ENT_NR;
                else
                    tpe = LDOC_ENT_OR;
                
                bool seqid = !strcmp(ky, "Seqid");
                bool source = !strcmp(ky, "Source");
                bool type = !strcmp(ky, "Type");

                // Convert defined keywords to lower case and hyphenate ('_' becomes '-'):
                if (*ky >= 'A' && *ky <= 'Z')
                    gen_lwrhyph(ky);
                
                if (seqid)
                {
                    ky = (char*)GEN_LANDMARKS; // Rename "seqid" to "landmarks"
                    gen_attr_t kwd = { BI_NKW, (char*)NULL };
                    seqid_nde = gen_csep_dup(sub, kwd, ky, val, true);
                    glbl = false;
                }
                else if (source)
                {
                    ky = (char*)GEN_SOURCES; // Rename "source" to "sources"
                    gen_attr_t kwd = { BI_NKW, (char*)NULL };
                    source_nde = gen_csep_dup(sub, kwd, ky, val, true);
                    glbl = false;
                }
                else if (type)
                {
                    ky = (char*)GEN_TYPES; // Rename "type" to "types"
                    gen_attr_t kwd = { BI_NKW, (char*)NULL };
                    type_nde = gen_csep_dup(sub, kwd, ky, val, true);
                    glbl = false;
                }
                else
                {
                    ldoc_ent_t* ent = ldoc_ent_new(tpe);
                    
                    // TODO Error handling.
                    
                    ent->pld.pair.anno.str = strdup(ky);
                    ent->pld.pair.dtm.str = strdup(val);
                    
                    ldoc_nde_ent_push(sub, ent);
                }
            }
            
            // End of string -- break condition 2:
            if (brk)
                break;
            
            // Reposition `ky` to point to the possible next key:
            ky = strct + 1;
        }
        
        strct++;
    } while(true); // See above for break conditions (1 & 2).
    
    if (glbl)
    {
        free(uid);
        
        uid = strdup(GEN_GLOBAL);
        
        // TODO Error handling.
        
        // Overwrite previous custom ID:
        sub->mkup.anno.str = uid;
        
        const char* global_id[] = { uid };
        ldoc_res_t* res_sub = ldoc_find_anno_nde(nde, (char**)global_id, 1);
        
        // Copy over info from `sub`, if such a node already exists in the data:
        if (res_sub)
        {
            ldoc_nde_t* exst = res_sub->info.nde;
            
            ldoc_res_free(res_sub);
            
            ldoc_nde_t* niter;
            TAILQ_FOREACH(niter, &(sub->dscs), ldoc_nde_entries)
            {
                ldoc_nde_rm(niter);
                
                ldoc_nde_dsc_push(exst, niter);
            }
            
            ldoc_ent_t* eiter;
            TAILQ_FOREACH(eiter, &(sub->ents), ldoc_ent_entries)
            {
                ldoc_ent_rm(eiter);
                
                ldoc_nde_ent_push(exst, eiter);
            }
            
            ldoc_nde_free(sub);
        }
        else
            ldoc_nde_dsc_push(nde, sub);
    }
    else
    {
        ldoc_nde_dsc_push(nde, sub);
     
        // Indexing:
        
        // Type does not matter much, since an ldoc_ent_t is only
        // used to conveniently store the information in an ldoc_nde_t
        // later.
        // Choosing LDOC_ENT_TXT, so that it is clear that the information
        // is placed in ent->pld.str:
        ldoc_ent_t* ent = ldoc_ent_new(LDOC_ENT_TXT);
        
        // TODO Error handling.
        
        ent->pld.str = uid;
        
        gvf_idx(seqid_nde, type_nde, source_nde, uid, 0, NULL, ent);
    }
}

static inline void gvf_proc_sctn(ldoc_nde_t* nde, char* sctn, char* ky, char* val)
{
    // Find previous `sctn` section (if it exists):
    ldoc_nde_t* s = NULL;
    ldoc_nde_t* iter;
    TAILQ_FOREACH(iter, &(nde->dscs), ldoc_nde_entries)
    {
        if (!strcmp(iter->mkup.anno.str, sctn))
        {
            s = iter;
            break;
        }
    }
    
    // If no section of name `sctn` exists yet, then create one now:
    if (!s)
    {
        s = ldoc_nde_new(LDOC_NDE_UA);
        s->mkup.anno.str = strdup(sctn);
        
        ldoc_nde_dsc_push(nde, s);
    }
    
    if (!strcmp(ky, GEN_TECHNOLOGY))
        gvf_tags(s, val, GVF_PGM_TECHNOLOGY);
    else if (!strcmp(ky, GEN_DATA_SRC))
        gvf_tags(s, val, GVF_PGM_DATA_SRC);
    else if (!strcmp(ky, GEN_SCORE_MTHD))
        gvf_tags(s, val, GVF_PGM_SCR_MTD);
    else if (!strcmp(ky, GEN_SOURCE_MTHD))
        gvf_tags(s, val, GVF_PGM_SRC_MTD);
    else if (!strcmp(ky, GEN_ATTRIBUTE_MTHD))
        gvf_tags(s, val, GVF_PGM_ATTR_MTD);
    else if (!strcmp(ky, GEN_PHENO_DESCR))
        gvf_tags(s, val, GVF_PGM_PHEN_DESC);
    else if (!strcmp(ky, GEN_PHASED_GENO))
        gvf_tags(s, val, GVF_PGM_PHSD_GT);
    else
    {
        // Find whether a global scope exists already:
        ldoc_nde_t* scp = NULL;
        TAILQ_FOREACH(iter, &(s->dscs), ldoc_nde_entries)
        {
            if (!strcmp(iter->mkup.anno.str, GEN_GLOBAL))
            {
                scp = iter;
                break;
            }
        }
        
        // No scope for the meta information present yet? Create a new node:
        if (!scp)
        {
            scp = ldoc_nde_new(LDOC_NDE_UA);
            scp->mkup.anno.str = strdup(GEN_GLOBAL);
            
            // TODO Error handling.
            
            ldoc_nde_dsc_push(s, scp);
        }
        
        // Determine data type; default: LDOC_ENT_OR
        ldoc_content_t tpe;
        if (!strcmp(ky, "technology-platform-read-length") ||
            !strcmp(ky, "technology-platform-read-pair-span") ||
            !strcmp(ky, "technology-platform-average-coverage"))
        {
            ky += 20;
            tpe = LDOC_ENT_NR;
        }
        else if (!strncmp(ky, "technology-platform-", 20))
        {
            ky += 20;
            tpe = LDOC_ENT_OR;
        }
        else
            tpe = LDOC_ENT_OR;
        
        ldoc_ent_t* ent = ldoc_ent_new(tpe);
        
        // TODO Error handling.
        
        ent->pld.pair.anno.str = strdup(ky);
        ent->pld.pair.dtm.str = strdup(val);
        
        ldoc_nde_ent_push(scp, ent);
    }
}

static inline ldoc_doc_t* gvf_proc_prgm(ldoc_doc_t* doc, char* ln, size_t lnlen, char** cmt)
{
    bool usr_nw;
    ldoc_nde_t* usr = gen_ctx(doc->rt, &usr_nw);

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
    
    ldoc_nde_t* stmt = gen_find_nde(doc->rt, usr, ln);
    
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
            ldoc_nde_t* dst;
            if (tpe != LDOC_NDE_OL ||
                !strcmp(ln, GEN_BUILD) ||
                !strcmp(ln, GEN_INDIVIDUALS_GVF1) ||
                !strcmp(ln, GEN_INDIVIDUALS_GVF2))
                dst = doc->rt;
            else
                dst = usr; // User defined pragmas.
            
            // TODO Use own types, so that this conversion is not necessary:
            stmt = ldoc_nde_new(tpe == LDOC_NDE_OO ? LDOC_NDE_UA : tpe);
            
            if (!strcmp(ln, GEN_SEQUENCE_REGION_GFF3))
                stmt->mkup.anno.str = strdup(GEN_SEQUENCE_REGION);
            else if (!strcmp(ln, GEN_INDIVIDUALS_GVF1) ||
                     !strcmp(ln, GEN_INDIVIDUALS_GVF2))
                stmt->mkup.anno.str = strdup(GEN_INDIVIDUALS);
            else if (!strncmp(ln, "technology-platform-", 20))
                stmt->mkup.anno.str = strdup(ln + 20);
            else
                stmt->mkup.anno.str = strdup(ln);
            
            ldoc_nde_dsc_push(dst, stmt);
        }
    }
    
    // Possible identifier for use in the comment section:
    char* id = NULL;
    
    if (tpe == LDOC_NDE_OO)
    {
        // Fine with string comparisons here, since this case
        // will (hopefully) not be true many times.
        // TODO Can GEN_DATA, SCORE, ATTRIBUTE, etc. be merged
        //      by passing ln instead of a constant?
        if (!strcmp(ln, "sequence-region"))
        {
            ldoc_nde_t* nde = gff_proc_sregion(val, &id);
            
            ldoc_nde_dsc_push(stmt, nde);
        }
        else if (!strcmp(ln, GEN_DATA_SRC))
        {
            gvf_proc_sctn(doc->rt, (char*)GEN_DATA_SRC, ln, val);
        }
        else if (!strcmp(ln, GEN_SCORE_MTHD))
        {
            gvf_proc_sctn(doc->rt, (char*)GEN_SCORE_MTHD, ln, val);
        }
        else if (!strcmp(ln, GEN_SOURCE_MTHD))
        {
            gvf_proc_sctn(doc->rt, (char*)GEN_SOURCE_MTHD, ln, val);
        }
        else if (!strcmp(ln, GEN_ATTRIBUTE_MTHD))
        {
            gvf_proc_sctn(doc->rt, (char*)GEN_ATTRIBUTE_MTHD, ln, val);
        }
        else if (!strcmp(ln, GEN_PHENO_DESCR))
        {
            gvf_proc_sctn(doc->rt, (char*)GEN_PHENO_DESCR, ln, val);
        }
        else if (!strcmp(ln, GEN_PHASED_GENO))
        {
            gvf_proc_sctn(doc->rt, (char*)GEN_PHASED_GENO, ln, val);
        }
        else if (!strncmp(ln, "technology-", 11) && strlen(ln) >= 12) // technology-.+
        {
            gvf_proc_sctn(doc->rt, (char*)GEN_TECHNOLOGY, ln, val);
        }
        else
        {
            // TODO Internal error.
        }
    }
    else if (tpe == LDOC_NDE_OL)
    {
        if (!strcmp(ln, GEN_BUILD))
        {
            gff_proc_gbld(stmt, val);
        }
        else if (!strcmp(ln, GEN_INDIVIDUALS_GVF2))
        {
            gen_csepstr_dup(stmt, val, true);
        }
        else
        {
            // Unknown pragmas -- and GEN_INDIVIDUALS_GVF1 (individual-id):
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
    
    // Add user-defined pragmas -- if those exist:
    if (usr_nw)
        gen_add_nw(doc->rt, usr);
    
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
                else if (ln[2] == '#')
                {
                    // Skip. "###" section separator.
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

// JSON to GVF

static inline void gvf_proc_doc_entlst(ldoc_nde_t* cntnr, char* sep)
{
    ldoc_ent_t* ent;
    TAILQ_FOREACH(ent, &(cntnr->ents), ldoc_ent_entries)
    {
        if (ent != TAILQ_FIRST(&(cntnr->ents)))
            qk_strcat(sep);
        
        qk_strcat(ent->pld.str);
    }
}

static inline void gvf_proc_doc_strct(ldoc_nde_t* prgm, char* ky, char* alt)
{
    char* pth[] = { ky };
    ldoc_res_t* res = ldoc_find_anno_nde(prgm, pth, 1);
    
    if (!res || !res->info.nde->dsc_cnt)
        return;

    ldoc_nde_t* cntnr;
    TAILQ_FOREACH(cntnr, &(res->info.nde->dscs), ldoc_nde_entries)
    {
        if (!qk_heap_empty())
            qk_strcat("\n");

        qk_strcat("##");
        qk_strcat(alt);
        qk_strcat(" ");
        
        bool emty = true;
        ldoc_nde_t* dsc;
        TAILQ_FOREACH(dsc, &(cntnr->dscs), ldoc_nde_entries)
        {
            if (dsc != TAILQ_FIRST(&(cntnr->dscs)))
            {
                emty = false;
                qk_strcat(";");
            }
            
            if (!strcmp(dsc->mkup.anno.str, GEN_LANDMARKS))
            {
                qk_strcat(GEN_LANDMARKS_GVF);
                qk_strcat("=");
                gvf_proc_doc_entlst(dsc, ",");
            }
            else if (!strcmp(dsc->mkup.anno.str, GEN_SOURCES))
            {
                qk_strcat(GEN_SOURCES_GVF);
                qk_strcat("=");
                gvf_proc_doc_entlst(dsc, ",");
            }
            else if (!strcmp(dsc->mkup.anno.str, GEN_TYPES))
            {
                qk_strcat(GEN_TYPES_GVF);
                qk_strcat("=");
                gvf_proc_doc_entlst(dsc, ",");
            }
            else if (!strcmp(dsc->mkup.anno.str, GEN_DBXREF))
            {
                qk_strcat(GEN_DBXREF_GFF3);
                qk_strcat("=");
                gvf_proc_doc_entlst(dsc, ",");
            }
            else
            {
                // TODO Document format error.
            }
        }
        
        ldoc_ent_t* ent;
        TAILQ_FOREACH(ent, &(cntnr->ents), ldoc_ent_entries)
        {
            if (!emty)
                qk_strcat(";");
            else
                emty = false;
            
            if (!strcmp(ent->pld.pair.anno.str, GEN_COMMENT))
            {
                qk_strcat(GEN_COMMENT_GVF);
                qk_strcat("=");
                qk_strcat(ent->pld.pair.dtm.str);
            }
            else
            {
                qk_strcat(ent->pld.pair.anno.str);
                qk_strcat("=");
                qk_strcat(ent->pld.pair.dtm.str);
            }
        }
    }
}

static inline void gvf_proc_doc_prgm(ldoc_nde_t* prgm)
{
    // Version info; gff-version
    gen_proc_doc_prgm_kv(prgm, (char*)GEN_GVFVERSION, (char*)GEN_GVFVERSION_GVF, " ");

    gff_proc_doc_prgm(prgm);

    gen_proc_doc_prgm_kv(prgm, "file-version", "file-version", " ");
    gen_proc_doc_prgm_kv(prgm, "file-date", "file-date", " ");
    
    const char* ind_pth[] = { GEN_INDIVIDUALS };
    ldoc_res_t* ind = ldoc_find_anno_nde(prgm, (char**)ind_pth, 1);
    if (ind)
    {
        if (ind->info.nde->ent_cnt == 1)
        {
            if (!qk_heap_empty())
                qk_strcat("\n");
            
            qk_strcat("##");
            qk_strcat(GEN_INDIVIDUALS_GVF1);
            qk_strcat(" ");
            qk_strcat(TAILQ_FIRST(&(ind->info.nde->ents))->pld.str);
        }
        else if (ind->info.nde->ent_cnt > 1)
        {
            if (!qk_heap_empty())
                qk_strcat("\n");
            
            qk_strcat("##");
            qk_strcat(GEN_INDIVIDUALS_GVF2);
            qk_strcat(" ");
            ldoc_ent_t* ind_ent;
            TAILQ_FOREACH(ind_ent, &(ind->info.nde->ents), ldoc_ent_entries)
            {
                if (ind_ent != TAILQ_FIRST(&(ind->info.nde->ents)))
                    qk_strcat(",");
                
                qk_strcat(ind_ent->pld.str);
            }
        }
    }
    
    gen_proc_doc_prgm_kv(prgm, (char*)GEN_POP, (char*)GEN_POP, " "); // population
    gen_proc_doc_prgm_kv(prgm, (char*)GEN_SEX, (char*)GEN_SEX, " "); // sex
    
    gvf_proc_doc_strct(prgm, (char*)GEN_TECHNOLOGY, (char*)GEN_TECHNOLOGY);
    gvf_proc_doc_strct(prgm, (char*)GEN_DATA_SRC, (char*)GEN_DATA_SRC);
    gvf_proc_doc_strct(prgm, (char*)GEN_SCORE_MTHD, (char*)GEN_SCORE_MTHD);
    gvf_proc_doc_strct(prgm, (char*)GEN_SOURCE_MTHD, (char*)GEN_SOURCE_MTHD);
    gvf_proc_doc_strct(prgm, (char*)GEN_ATTRIBUTE_MTHD, (char*)GEN_ATTRIBUTE_MTHD);
    gvf_proc_doc_strct(prgm, (char*)GEN_PHENO_DESCR, (char*)GEN_PHENO_DESCR);
    gvf_proc_doc_strct(prgm, (char*)GEN_PHASED_GENO, (char*)GEN_PHASED_GENO);
    
    gen_proc_doc_prgm(prgm, " ");
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

char* gvf_proc_doc(ldoc_doc_t* doc, gen_doctype_t tpe)
{
    char* attr;
    
    switch (tpe)
    {
        case GEN_FMT_INF:
            gvf_proc_doc_prgm(doc->rt);
            
            return qk_heap_ptr();
        case GEN_FMT_FTR:
            attr = gvf_proc_doc_ftr_attrs(doc->rt);
            
            gff_proc_doc_ftr(doc->rt);
            
            if (*attr)
                qk_strcat(";");
            qk_strcat(attr);
            
            free(attr);
            
            return qk_heap_ptr();
        default:
            // TODO Internal error.
            return NULL;
    }
    
    // TODO Obsolete?
    // gen_proc_doc_usr(doc->rt);
    

}
