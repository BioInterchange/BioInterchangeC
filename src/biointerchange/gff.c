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

#include "gff.h"

const char* GFF_FA = "\n##FASTA";
const char* GFF_FA_PFX = "##FASTA";

const char* GFF_C1  = "landmark";
const char* GFF_C2  = "source";
const char* GFF_C3  = "type";
const char* GFF_C4  = "start";
const char* GFF_C5  = "end";
const char* GFF_C6  = "score";
const char* GFF_C7  = "strand";
const char* GFF_C8  = "codon-phase";

void gff_cbcks(gen_cbcks_t* cbcks)
{
    cbcks->proc_ln = &gff_proc_ln;
}

static inline void gff_attr()
{
    
}

static inline char* gen_nxt_ws(char* str)
{
    while (*str != ' ' &&
           *str != '\t')
    {
        if (!*str || *str == '\n')
            return NULL;
        
        str++;
    }
    
    // Remember where to place string terminator -- if this function is successful:
    char* sep = str;
    
    while (*str == ' ' ||
           *str == '\t')
    {
        if (!*str || *str == '\n')
            return NULL;
        
        str++;
    }
    
    if (!*str || *str == '\n')
        return NULL;
    
    *sep = 0;
    
    return str;
}

char* gff_proc_cmt(char* ln, size_t lnlen)
{
    // Skip initial character that marks the comment:
    ln++;
    lnlen--;

    // If there is a whitespace, skip over this one too, but
    // leave further whitespace characters as they are:
    if (*ln == ' ' || *ln == '\t')
    {
        ln++;
        lnlen--;
    }
    
    // "Raw" comment, i.e. leading white space removed:
    char* cln = strndup(ln, lnlen);
    
    return cln;
}

/**
 * Optional value handling; duplicates string.
 */
static inline char* gff_proc_optval(char* str)
{
    if (*str == '.' && *(str + 1) == 0)
        return NULL;
    
    // TODO Error handling.
    return strdup(str);
}

/**
 * Optional value handling; keeps string.
 */
static inline char* gff_proc_optvalx(char* str)
{
    if (*str == '.' && *(str + 1) == 0)
        return NULL;
    
    return str;
}

inline void gff_proc_tgt(ldoc_nde_t* nde, char* val)
{
    ldoc_ent_t* id = ldoc_ent_new(LDOC_ENT_OR);
    ldoc_ent_t* st = ldoc_ent_new(LDOC_ENT_NR);
    ldoc_ent_t* en = ldoc_ent_new(LDOC_ENT_NR);
    ldoc_ent_t* strnd = ldoc_ent_new(LDOC_ENT_OR);
    
    // TODO Error handling.

    id->pld.pair.anno.str = (char*)GEN_ID;
    st->pld.pair.anno.str = (char*)GEN_START;
    en->pld.pair.anno.str = (char*)GEN_END;
    strnd->pld.pair.anno.str = (char*)GEN_STRAND;

    char* s = val;
    while (*val && *val != ' ' && *val != ';')
        val++;
    *val = 0;
    
    id->pld.pair.dtm.str = s;
    s = ++val;

    while (*val && *val != ' ' && *val != ';')
        val++;
    *val = 0;
    
    st->pld.pair.dtm.str = s;
    s = ++val;

    while (*val && *val != ' ' && *val != ';')
        val++;
    bool strndd = *val ? true : false;
    *val = 0;
    
    en->pld.pair.dtm.str = s;
    
    if (strndd)
    {
        s = ++val;

        while (*val && *val != ' ' && *val != ';')
            val++;
        *val = 0;
    }
    else
        s = NULL;
    
    strnd->pld.pair.dtm.str = s;
    s = ++val;

    ldoc_nde_ent_push(nde, id);
    ldoc_nde_ent_push(nde, st);
    ldoc_nde_ent_push(nde, en);
    ldoc_nde_ent_push(nde, strnd);
}

void gff_proc_gbld(ldoc_nde_t* nde, char* val)
{
    // Source:
    char* src = val;
    
    // Build:
    char* bld = gen_nxt_ws(val);
    
    if (!bld)
    {
        // TODO Error handling.
    }
    
    // Allocate node and entities upfront, so that strdup's do
    // not need to be freed if any of this fails.
    
    ldoc_nde_t* nde_bld = ldoc_nde_new(LDOC_NDE_UA);
    
    // TODO Error handling.
    
    ldoc_ent_t* ent_bld = ldoc_ent_new(LDOC_ENT_OR);
    
    if (!ent_bld)
    {
        // TODO Error handling.
    }
    
    ent_bld->pld.pair.anno.str = (char*)GEN_SOURCE;
    ent_bld->pld.pair.dtm.str = strdup(src);
    
    ldoc_nde_ent_push(nde_bld, ent_bld);
    
    ent_bld = ldoc_ent_new(LDOC_ENT_OR);
    
    if (!ent_bld)
    {
        // TODO Error handling.
    }
    
    ent_bld->pld.pair.anno.str = (char*)GEN_BUILD_VAL;
    ent_bld->pld.pair.dtm.str = strdup(bld);
    
    ldoc_nde_ent_push(nde_bld, ent_bld);
    
    ldoc_nde_dsc_push(nde, nde_bld);
}

ldoc_nde_t* gff_proc_sregion(char* val, char** cid)
{
    // Landmark identifier:
    char* id = val;

    // Start coordinate:
    char* st = gen_nxt_ws(val);

    if (!st)
    {
        // TODO Error handling.
    }
    
    // End coordinate:
    char* en = gen_nxt_ws(st);
    
    if (!en)
    {
        // TODO Error handling.
    }
    
    // Allocate node and entities upfront, so that strdup's do
    // not need to be freed if any of this fails.
    
    ldoc_nde_t* nde = ldoc_nde_new(LDOC_NDE_UA);

    if (!nde)
    {
        // TODO Error handling.
    }

    ldoc_ent_t* ent_st = ldoc_ent_new(LDOC_ENT_NR);
    
    if (!ent_st)
    {
        // TODO Error handling.
    }
    
    ldoc_ent_t* ent_en = ldoc_ent_new(LDOC_ENT_NR);
    
    if (!ent_en)
    {
        // TODO Error handling.
        
    }
    
    nde->mkup.anno.str = strdup(id);
    
    ent_st->pld.pair.anno.str = strdup(GEN_START);
    ent_st->pld.pair.dtm.str = gff_proc_optval(st);
    
    ent_en->pld.pair.anno.str = strdup(GEN_END);
    ent_en->pld.pair.dtm.str = gff_proc_optval(en);
    
    ldoc_nde_ent_push(nde, ent_st);
    ldoc_nde_ent_push(nde, ent_en);
    
    if (cid)
    {
        *cid = strdup(id);
    }
    
    return nde;
}

ldoc_struct_t gff_prgm_tpe(char* ky)
{
    if (*ky == 'g')
    {
        ky++;
        if (*ky == 'e')
        {
            ky++;
            if (*ky == 'n')
            {
                ky++;
                if (!strcmp(ky, "ome-build")) // genome-build
                    return LDOC_NDE_OL;
            }
        } else if (*ky == 'f')
        {
            ky++;
            if (*ky == 'f')
            {
                ky++;
                if (!strcmp(ky, "-version")) // gff-version
                    return LDOC_NDE_UA;
            }
        }
    }
    else if (*ky == 's')
    {
        ky++;
        if (*ky == 'e')
        {
            ky++;
            if (*ky == 'q')
            {
                ky++;
                if (!strcmp(ky, "uence-region")) // sequence-region
                    return LDOC_NDE_OO;
            }
        }
        else if (*ky == 'p')
        {
            ky++;
            if (*ky == 'e')
            {
                ky++;
                if (!strcmp(ky, "cies")) // species
                    return LDOC_NDE_UA;
            }
        }
    }
    
    return LDOC_NDE_OL;
}

static inline ldoc_doc_t* gff_proc_prgm(ldoc_doc_t* doc, char* ln, size_t lnlen, char** cmt)
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
    
    // Ignore '###' separators:
    if (*ln == '#' && !*(ln + 1))
        return NULL;
    
    ldoc_nde_t* stmt = gen_find_nde(doc->rt, usr, ln);
    ldoc_struct_t tpe = gff_prgm_tpe(ln);
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
                !strcmp(ln, "genome-build"))
                dst = doc->rt;
            else
                dst = usr; // User defined pragmas.
            
            // TODO Use own types, so that this conversion is not necessary:
            stmt = ldoc_nde_new(tpe == LDOC_NDE_OO ? LDOC_NDE_UA : tpe);
            
            if (!strcmp(ln, GEN_SEQUENCE_REGION_GFF3))
                stmt->mkup.anno.str = strdup(GEN_SEQUENCE_REGION);
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
            cent->pld.pair.dtm.str = gen_escstr(*cmt, GEN_FMT_GFF3);
        }
        else
        {
            cent->pld.pair.anno.str = strdup(val);
            cent->pld.pair.dtm.str = gen_escstr(*cmt, GEN_FMT_GFF3);

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

static inline ldoc_doc_t* gff_proc_ftr(int fd, off_t mx, ldoc_trie_t* idx, char* ln, size_t lnlen, char** cmt)
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
    uint8_t col = 1;
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
    
    /* Reimplement later where only used "things" are defined:
    ldoc_nde_t* ctx = ldoc_nde_new(LDOC_NDE_UA);
    ctx->mkup.anno.str = (char*)JSONLD_CTX;
     */
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_ent_t* ctx = ldoc_ent_new(LDOC_ENT_OR);
    ctx->pld.pair.anno.str = (char*)JSONLD_CTX;
    ctx->pld.pair.dtm.str = (char*)JSONLD_GFF3_1;
    ldoc_nde_ent_push(ftr, ctx);

    // User-defined attributes:
    ldoc_nde_t* attrs = ldoc_nde_new(LDOC_NDE_UA);
    attrs->mkup.anno.str = (char*)GEN_ATTRS;
    
    ldoc_ent_t* lm = ldoc_ent_new(LDOC_ENT_OR);
    lm->pld.pair.anno.str = (char*)GFF_C1;
    lm->pld.pair.dtm.str = coff[0];

    ldoc_nde_t* lc = ldoc_nde_new(LDOC_NDE_UA);
    lc->mkup.anno.str = (char*)GEN_LOCUS;

    ldoc_ent_t* src = ldoc_ent_new(LDOC_ENT_OR);
    src->pld.pair.anno.str = (char*)GFF_C2;
    src->pld.pair.dtm.str = coff[1];

    ldoc_ent_t* tpe = ldoc_ent_new(LDOC_ENT_OR);
    tpe->pld.pair.anno.str = (char*)GFF_C3;
    tpe->pld.pair.dtm.str = coff[2];

    ldoc_ent_t* st = ldoc_ent_new(LDOC_ENT_NR);
    st->pld.pair.anno.str = (char*)GFF_C4;
    st->pld.pair.dtm.str = coff[3];
    
    ldoc_ent_t* en = ldoc_ent_new(LDOC_ENT_NR);
    en->pld.pair.anno.str = (char*)GFF_C5;
    en->pld.pair.dtm.str = coff[4];

    ldoc_ent_t* scr = ldoc_ent_new(LDOC_ENT_OR);
    scr->pld.pair.anno.str = (char*)GFF_C6;
    scr->pld.pair.dtm.str = coff[5];

    ldoc_ent_t* strnd = ldoc_ent_new(LDOC_ENT_OR);
    strnd->pld.pair.anno.str = (char*)GFF_C7;
    strnd->pld.pair.dtm.str = coff[6];

    ldoc_ent_t* ph = ldoc_ent_new(LDOC_ENT_NR);
    ph->pld.pair.anno.str = (char*)GFF_C8;
    ph->pld.pair.dtm.str = coff[7];

    // Add comment lines -- if available:
    if (*cmt)
    {
        ldoc_ent_t* c = ldoc_ent_new(LDOC_ENT_OR);
        c->pld.pair.anno.str = (char*)GEN_COMMENT;
        c->pld.pair.dtm.str = gen_escstr(*cmt, GEN_FMT_GFF3);
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
    
    gen_splt_attrs(ftr, attrs, NULL, NULL, coff[8], BI_VAL);

    // Source, type, score, strand, phase:
    ldoc_nde_ent_push(ftr, src);
    ldoc_nde_ent_push(ftr, tpe);
    ldoc_nde_ent_push(ftr, scr);
    ldoc_nde_ent_push(ftr, ph);

    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, lm);
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    ldoc_nde_ent_push(lc, strnd);
    
    ldoc_nde_dsc_push(ftr, lc);
    
    // Do not add user defined sub-tree if it is empty:
    if (attrs->dsc_cnt || attrs->ent_cnt)
        ldoc_nde_dsc_push(ftr, attrs);

    return doc;
}

void gff_proc_fa(char* ln, size_t lnlen)
{
    
}

inline bool gff_srch_chr(char* txt, char c)
{
    if (c == 0)
    {
        return *txt == '\n' || *txt  == '\t' || *txt == ' ';
    }
    
    return *txt == c;
}

static inline off_t gff_fmatch(char* txt, const int len, const char* str, int idx, off_t off, const bool eol_ws)
{
    if (eol_ws)
    {
        if (str[idx])
            return 0;
        
        txt++;
        
        int ws = 1;
        bool seen_nl = false;
        while (off < len) {
            if (*txt == ' ' || *txt == '\t')
                ws++;
            else if (*txt == '\n' || *txt == '\r')
            {
                seen_nl = true;
                ws++;
            }
            else
            {
                if (seen_nl)
                    return ws;
                else
                    return 0;
            }
            
            txt++;
            off++;
        }
        
        return ws;
    }
    else if (!str[idx])
        return 1;
    
    return 0;
}

/**
 * Finds the first position of `str` in `txt`.
 *
 * str needs to be a proper prefix. For example, 'hello' is okay, 'hehello' is not.
 * stat can be NULL, if no statistics should be collected.
 *
 * @return Offset of the found string relative to the beginning of `txt`; -1, if the string
 *         was not found.
 */
static inline off_t gff_ffind(char* txt, const int len, const char* str, gen_fstat* stat, off_t stat_skip, const bool eol_ws, const bool rt_end)
{
    off_t off = 0;
    int idx = 0;
    
    bool last_nl = off == 0 ? true : (*(txt - 1) == '\n' || *(txt - 1) == '\r' ? true : false);
    bool last_hash = false;
    
    // Yes, this handles cases where the input is zero characters!
    while (off < len)
    {
        if (stat)
        {
            if (*txt == '\n')
            {
                last_nl = true;
                last_hash = false;
            }
            else if (last_nl && *txt == '#')
            {
                last_nl = false;
                last_hash = true;
            }
            else if (last_nl)
            {
                last_nl = false;
                if (off >= stat_skip)
                    stat->ftrs++;
            }
            else if (last_hash && *txt == '#')
            {
                last_nl = false;
                last_hash = false;
                if (off >= stat_skip)
                    stat->meta++;
            }
            else if (last_hash)
            {
                last_nl = false;
                last_hash = false;
                if (off >= stat_skip)
                    stat->comms++;
            }
        }
        
        off++;
        
        if (*txt == str[idx])
        {
            idx++;
            
            // Matching complete?
            off_t m = gff_fmatch(txt, len, str, idx, off, eol_ws);
            if (m > 0)
            {
                // Update valid statistics offset:
                if (stat)
                    stat->off = off;
                
                // Match at this position:
                if (rt_end)
                    // Minus 1, because gff_fmatch returns the number of white space
                    // characters + 1:
                    return off + m - 1;
                else
                    return off - idx + 1;
            }
        }
        else
            idx = 0;
        
        txt++;
    }
    
    if (stat)
        stat->off += len;
    
    // Not found:
    return BI_NFOUND;
}

off_t gff_fnd_fa(int fd, gen_fstat* stat, off_t mx)
{
    // 0: no FASTA section found.
    off_t fa = 0;
    
    off_t skip = 0;
    const int pg_size = getpagesize();
    
    // Sanity check: If page size smaller than 3 bytes, then lookback below will cause an error.
    if (pg_size < 3)
    {
        // TODO Error.
    }
    
    fio_mem* mem = NULL;
    const size_t len = pg_size * BI_FA_IDX_MUL;
    while (fa < mx)
    {
        mem = fio_mmap(mem, fd, mx, len, fa);
        
        // TODO Error check.
        
        if (fa > 0)
            skip = pg_size;
        
        // Commence the actual search:
        off_t fa_off = gff_ffind(mem->pg, mem->ln, GFF_FA, stat, skip, true, true);
        
        if (fa_off != BI_NFOUND)
            return fa_off;
        
        // Overlapping search: last/first page searched twice:
        fa += len - pg_size;
    }
    
    if (mem)
        fio_munmap(mem);
    
    // Not found:
    return BI_NFOUND;
}

static inline void gff_idx_fa_pg(ldoc_trie_t* trie, off_t abs, char* txt, size_t len)
{
    off_t off = 0;
    
    char id[BI_FA_ID_MAXLEN];
    id[BI_FA_ID_MAXLEN - 1] = 0;
    
    while (off < len)
    {
        off = gff_ffind(txt + off, len - off, "\n>", NULL, 0, false, true);
        
        // String not found; we are done with this input:
        if (off == BI_NFOUND)
            return;
        
        // Extract sequence ID:
        id[0] = 0;
        uint16_t i;
        uint16_t j;
        for (i = 0; i < BI_FA_ID_MAXLEN - 1; i++)
        {
            // Identifier extends into the next page:
            if (off + i >= len)
                return;
            
            if (*(txt + off + i) == '\n' ||
                *(txt + off + i) == '\r' ||
                *(txt + off + i) == ' ' ||
                *(txt + off + i) == '\t')
            {
                id[i] = 0;
                
                // Bail-out:
                j = i;
                i = BI_FA_ID_MAXLEN;
            }
            else
            {
                id[i] = *(txt + off + i);
            }
        }
        
        // Make sure that empty identifiers are skipped (FASTA format error):
        if (!id[0])
            return;

        // Advance to line ending:
        if (*(txt + off + i) != '\n' &&
            *(txt + off + i) != '\r')
            for (i = j + 1; i < BI_FA_ID_MAXLEN; i++)
            {
                // Reached next page:
                if (off + i >= len)
                    return;
                
                if (*(txt + off + i) == '\n' ||
                    *(txt + off + i) == '\r')
                {
                    // Bail-out:
                    j = i;
                    i = BI_FA_ID_MAXLEN;
                }
            }
        
        // Find sequence line length and type of line ending ('\n', '\r\n', etc.):
        size_t llen = 0; // TODO Could calculate llen based on lbrk and i.
        size_t lbrk = 0;
        for (i = j; i < BI_FA_ID_MAXLEN - 1; i++)
        {
            // Running into the next page.
            if (off + i >= len)
                return;
            
            if (!llen)
            {
                if (*(txt + off + i) == '\n' ||
                    *(txt + off + i) == '\r')
                    lbrk++;
                else
                {
                    // Length of line ending determined, proceed to find length of sequences:
                    j = i;
                    llen++;
                }
            }
            else
            {
                if (*(txt + off + i) != '\n' &&
                    *(txt + off + i) != '\r')
                {
                    llen++;
                }
                else
                {
                    // Bail-out:
                    i = BI_FA_ID_MAXLEN;
                }
            }
        }
        
        ldoc_trie_anno_t anno;
        
        // TODO Error handling.
        
        gen_fa_ntry_t* ntry = (gen_fa_ntry_t*)malloc(sizeof(gen_fa_ntry_t));
        
        // TODO Error handling.
        
        ntry->off = abs + off;
        ntry->llen = llen;
        ntry->lbrk = lbrk;
        ntry->seq = abs + off + j;
        
        anno.cat = BI_FASTA_REF;
        anno.pld = ntry;

        ldoc_trie_add(trie, id, ASCII, anno);
    }
}

ldoc_trie_t* gff_idx_fa(int fd, gen_fstat* stat, off_t mx)
{
    // Find beginning of first FASTA line:
    off_t fa = gff_fnd_fa(fd, stat, mx);
    
    if (fa != BI_NFOUND)
    {
        // "Step back" a character, so that the code for finding identifiers
        // is universal below (no special treatment of any lines; character
        // added due to the `fa--` should be a newline):
        fa--;

        ldoc_trie_t* fa_trie = ldoc_trie_new();
        const size_t pg_size = getpagesize();
        const off_t fa_pg = fa / pg_size;
        
        off_t fa_pg_off = fa - (fa_pg * pg_size);
        
        fio_mem* mem = NULL;
        const size_t len = pg_size * BI_FA_IDX_MUL;
        while (fa < mx)
        {
            mem = fio_mmap(mem, fd, mx, len, fa_pg * pg_size);
            
            // TODO Error check.

            gff_idx_fa_pg(fa_trie, fa, mem->pg + fa_pg_off, mem->ln - fa_pg_off);
            
            fa += len - pg_size;
            fa_pg_off = 0;
        }
        return fa_trie;
    }
    
    return NULL;
}

void gff_idx_release(ldoc_trie_t* trie)
{
    ldoc_trie_free(trie);
}

ldoc_doc_t* gff_proc_ln(int fd, off_t mx, ldoc_doc_t* fdoc, ldoc_trie_t* idx, char* ln, size_t lnlen, gen_prsr_t* st, char** cmt, gen_fstat* stat)
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
                    if (gff_proc_prgm(fdoc, ln, lnlen, cmt))
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
            ldoc = gff_proc_ftr(fd, mx, idx, ln, lnlen, cmt);
            stat->ftrs++;
        }
    }
    
    return ldoc;
}

char* gff_seq(int fd, off_t mx, ldoc_trie_t* idx, const char* id, off_t st, off_t en, bool rv)
{
    ldoc_trie_nde_t* fa = ldoc_trie_lookup(idx, id, false);
    
    // ID could not be found:
    if (!fa)
        return NULL;
    
    // TODO Check annotation type.
    
    // Adjust start coordinate to match 0-based memory indexing:
    st--;
    
    const gen_fa_ntry_t* ntry = fa->anno.pld;
    const size_t pg_size = getpagesize();
    const off_t seq_pg = ntry->seq / pg_size;
    const off_t seq_pg_off = ntry->seq - (seq_pg * pg_size);
    
    fio_mem* mem = NULL;
    size_t lraw = en;
    size_t lbrk = (lraw / ntry->llen) * ntry->lbrk;
    size_t mlen = (((lraw + lbrk + seq_pg_off) / pg_size) + 1) * pg_size;
    mem = fio_mmap(mem, fd, mx, mlen, seq_pg * pg_size);
    
    // TODO Error check.
    
    size_t slen = en - st;
    char* seq = (char*)malloc(slen + 1);
    
    if (!seq)
    {
        // TODO Error handling.
    }
    
    // Stop the string:
    seq[slen] = 0;
    
    // Add in line break characters:
    off_t fa_off = (st / ntry->llen) * ntry->lbrk + st + ntry->seq;
    
    // Set source pointer; advance from here:
    char* src = mem->pg + (fa_off - mem->off);
    char* dst = rv ? seq + (slen - 1) : seq;
    while (slen)
    {
        // Sanity check, in case sequence boundaries were not specified:
        if (mx < (off_t)src - (off_t)mem->pg + mem->off ||
            *src == '>' ||
            *src == ';')
        {
            free(seq);
            fio_munmap(mem);
            
            return NULL;
        }
        
        if (*src == '\n' ||
            *src == '\r')
            src++;
        else
        {
            if (rv)
                *(dst--) = gen_inv(*(src++));
            else
                *(dst++) = *(src++);
            slen--;
        }
    }
    
    fio_munmap(mem);
    
    return seq;
}

//
// JSON to GFF3
//

static inline void gff_join_almnt(ldoc_nde_t* nde, char* attrs)
{
    ldoc_res_t* id = ldoc_find_anno_ent(nde, (char*)GEN_ID);
    ldoc_res_t* st = ldoc_find_anno_ent(nde, (char*)GEN_START);
    ldoc_res_t* en = ldoc_find_anno_ent(nde, (char*)GEN_END);
    ldoc_res_t* strnd = ldoc_find_anno_ent(nde, (char*)GEN_STRAND);
    ldoc_res_t* cgr = ldoc_find_anno_ent(nde, (char*)GEN_CIGAR);
    
    // NOTE Report error, if one is missing?
    if (id && st && en)
    {
        gen_join_attrs_key((char*)GEN_ALIGNMENT_GFF3, NULL, NULL, attrs);
        qk_strcat("=");
        qk_strcat(id->info.ent->pld.pair.dtm.str);
        qk_strcat(" ");
        qk_strcat(st->info.ent->pld.pair.dtm.str);
        qk_strcat(" ");
        qk_strcat(en->info.ent->pld.pair.dtm.str);
        if (strnd)
        {
            if (strnd->info.ent->pld.pair.dtm.str)
            {
                qk_strcat(" ");
                qk_strcat(strnd->info.ent->pld.pair.dtm.str);
            }
            
            ldoc_res_free(strnd);
        }
        
        ldoc_res_free(id);
        ldoc_res_free(st);
        ldoc_res_free(en);
    }
    
    if (cgr)
    {
        gen_join_attrs_key((char*)GEN_CIGAR_GFF3, NULL, NULL, attrs);
        qk_strcat("=");
        gen_qk_revcig(cgr->info.ent->pld.pair.dtm.str);
        
        ldoc_res_free(cgr);
    }
}

inline void gff_proc_doc_prgm(ldoc_nde_t* prgm)
{
    // Version info; gff-version
    gen_proc_doc_prgm_kv(prgm, (char*)GEN_GFFVERSION, (char*)GEN_GFFVERSION_GFF3, " ");
    
    // Genome build; genome-build
    char* gb_pth[] = { (char*)GEN_BUILD };
    ldoc_res_t* gb = ldoc_find_anno_nde(prgm, gb_pth, 1);
    
    if (gb)
    {
        ldoc_nde_t* bld;
        TAILQ_FOREACH(bld, &(gb->info.nde->dscs), ldoc_nde_entries)
        {
            ldoc_res_t* bld_src = ldoc_find_anno_ent(bld, (char*)GEN_SOURCE);
            
            ldoc_res_t* bld_nme = ldoc_find_anno_ent(bld, (char*)GEN_BUILD_VAL);
            
            if (!qk_heap_empty())
                qk_strcat("\n");
            
            qk_strcat("##");
            qk_strcat(GEN_BUILD_GFF3);
            qk_strcat(" ");
            qk_strcat(bld_src->info.ent->pld.pair.dtm.str);
            qk_strcat(" ");
            qk_strcat(bld_nme->info.ent->pld.pair.dtm.str);
        }
    }
    
    // Sequence region; sequence-region (but "contig" in doc)
    char* sr_pth[] = { (char*)GEN_SEQUENCE_REGION };
    ldoc_res_t* sr = ldoc_find_anno_nde(prgm, sr_pth, 1);
    
    if (sr)
    {
        ldoc_nde_t* reg;
        TAILQ_FOREACH(reg, &(sr->info.nde->dscs), ldoc_nde_entries)
        {
            ldoc_res_t* reg_st = ldoc_find_anno_ent(reg, (char*)GEN_START);
            
            ldoc_res_t* reg_en = ldoc_find_anno_ent(reg, (char*)GEN_END);
            
            if (!qk_heap_empty())
                qk_strcat("\n");
            
            qk_strcat("##");
            qk_strcat(GEN_SEQUENCE_REGION_GFF3);
            qk_strcat(" ");
            qk_strcat(reg->mkup.anno.str);
            qk_strcat(" ");
            qk_strcat(reg_st->info.ent->pld.pair.dtm.str);
            qk_strcat(" ");
            qk_strcat(reg_en->info.ent->pld.pair.dtm.str);
        }
    }
    
    gen_proc_doc_prgm(prgm, " ");
}

inline void gff_proc_doc_ftr(ldoc_nde_t* ftr)
{
    // Covered:
    //   - id -> ID
    //   - name -> Name
    //   - dbxref -> Dbxref
    //   - parent -> Parent
    // ! - alignment -> Gap
    // ! - target -> Target
    const char* lc_id[] = { GEN_LOCUS };
    ldoc_res_t* lc = ldoc_find_anno_nde(ftr, (char**)lc_id, 1);

    ldoc_res_t* lm = ldoc_find_anno_ent(lc->info.nde, (char*)GFF_C1);
    ldoc_res_t* src = ldoc_find_anno_ent(ftr, (char*)GFF_C2);
    ldoc_res_t* tpe = ldoc_find_anno_ent(ftr, (char*)GFF_C3);
    ldoc_res_t* st = ldoc_find_anno_ent(lc->info.nde, (char*)GFF_C4);
    ldoc_res_t* en = ldoc_find_anno_ent(lc->info.nde, (char*)GFF_C5);
    ldoc_res_t* scr = ldoc_find_anno_ent(ftr, (char*)GFF_C6);
    ldoc_res_t* strnd = ldoc_find_anno_ent(ftr, (char*)GFF_C7);
    ldoc_res_t* ph = ldoc_find_anno_ent(ftr, (char*)GFF_C8);
    
    qk_strcat(gen_res_req(lm));
    qk_strcat("\t");
    qk_strcat(gen_res_opt(src));
    qk_strcat("\t");
    qk_strcat(gen_res_opt(tpe));
    qk_strcat("\t");
    qk_strcat(gen_res_opt(st));
    qk_strcat("\t");
    qk_strcat(gen_res_opt(en));
    qk_strcat("\t");
    qk_strcat(gen_res_opt(scr));
    qk_strcat("\t");
    qk_strcat(gen_res_opt(strnd));
    qk_strcat("\t");
    qk_strcat(gen_res_opt(ph));
    qk_strcat("\t");
    
    ldoc_res_t* id = ldoc_find_anno_ent(ftr, (char*)GEN_ID);
    ldoc_res_t* nme = ldoc_find_anno_ent(ftr, "name");
    
    const char* dbxref_id[] = { "dbxref" };
    ldoc_res_t* dbxref = ldoc_find_anno_nde(ftr, (char**)dbxref_id, 1);
    
    const char* prnt_pth[] = { "parent" };
    ldoc_res_t* prnt = ldoc_find_anno_nde(ftr, (char**)prnt_pth, 1);
    
    const char* almnt_pth[] = { GEN_ALIGNMENT };
    ldoc_res_t* almnt = ldoc_find_anno_nde(ftr, (char**)almnt_pth, 1);
    
    char* attrs = qk_working_ptr();
    
    if (id && !id->nde)
        gen_join_attrs_ent((char*)GEN_ID_GFF3, id->info.ent, attrs);
    
    if (nme && !nme->nde)
        gen_join_attrs_ent("Name", nme->info.ent, attrs);
    
    if (prnt && prnt->nde)
        gen_join_attrs_nde("Parent", prnt->info.nde, attrs);

    if (dbxref && dbxref->nde)
        gen_join_attrs_nde("Dbxref", dbxref->info.nde, attrs);
    
    if (almnt && almnt->nde)
        gff_join_almnt(almnt->info.nde, attrs);
    
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
}

char* gff_proc_doc(ldoc_doc_t* doc, gen_doctype_t tpe)
{
    // Attributes are joined on the quick heap:
    qk_purge();
    
    // Make sure that the quick heap contains a valid string even
    // in the case where no attributes might be present at all:
    qk_strcat("");
    
    switch (tpe)
    {
        case GEN_FMT_INF:
            gff_proc_doc_prgm(doc->rt);
            
            return qk_heap_ptr();
        case GEN_FMT_FTR:
            gff_proc_doc_ftr(doc->rt);
            
            return qk_heap_ptr();
        default:
            // TODO Internal error.
            return NULL;
    }
    
    //printf("%s\n", qk_heap_ptr());
}
