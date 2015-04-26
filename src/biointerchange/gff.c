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

static const char* GFF_C1  = "seqid";
static const char* GFF_C2  = "source";
static const char* GFF_C3  = "type";
static const char* GFF_C4  = "start";
static const char* GFF_C5  = "end";
static const char* GFF_C6  = "score";
static const char* GFF_C7  = "strand";
static const char* GFF_C8  = "phase";

static const char* GEN_LOCUS = "locus";
static const char* GEN_ATTRS = "user";

static const char* GEN_NULL = "null";
static const char* GEN_TRUE = "true";

static inline void gen_lwr(char* str)
{
    while (*str)
    {
        if (*str >= 'A' && *str <= 'Z')
            *str = *str - 'A' + 'a';
        
        str++;
    }
}

static inline bi_attr gen_kwd(char* str)
{
    if (*str >= 'A' && *str <= 'Z')
        return BI_VAL;
    
    return BI_NKW;
}

static inline char gen_inv(char c)
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
}

static inline void gff_ky(char* attr, char** val)
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

void gff_splt_attrs(ldoc_nde_t* ftr, ldoc_nde_t* usr, char* attrs)
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
            gff_ky(attr, &val);
         
            // Nothing to handle -- break condition 1:
            if (!*attr)
                return;
            
            ldoc_ent_t* kv = ldoc_ent_new(LDOC_ENT_OR);
            
            if (!kv)
            {
                // TODO Error handling.
            }
            
            // Key/value assignment, or, key-only handling:
            ldoc_nde_t* dst;
            bi_attr kind;
            if (kind = gen_kwd(attr))
            {
                dst = ftr;
                
                gen_lwr(attr);
            }
            else
            {
                dst = usr;
                
                if (!val)
                    val = (char*)GEN_TRUE;
            }
            
            kv->pld.pair.anno.str = attr;
            kv->pld.pair.dtm.str = val;
            ldoc_nde_ent_push(dst, kv);
            
            // Attribute EOS reached earlier -- break condition 2:
            if (brk)
                return;
            
            // Next attribute will have to start here:
            attr = attrs + 1;
        }
        
        attrs++;
    } while (true); // See "break condition 1" and "break condition 2".
}

void gff_proc_prgm(char* ln, size_t lnlen)
{

}

static inline ldoc_doc_t* gff_proc_ftr(char* ln, size_t lnlen)
{
    ldoc_doc_t* doc = ldoc_doc_new();
    
    if (!doc)
    {
        // TODO Error handling.
    }
    
    // Remove trailing line breaks:`
    off_t off = lnlen - 1;
    while ((ln[off] == '\n' || ln[off] == '\r') && off > 0) {
        ln[off--] = 0;
    }
    
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
        if (coff[col - 1][0] == '.' && coff[col - 1][1] == 0)
            coff[col - 1] = NULL;
    }
    
    // 0 : landmark
    // 1 : 1
    // 2 : 2
    // 3 : start
    // 4 : end
    // 5 : 5
    // 6 : strand
    // 7 : 7
    // 8 : attr

    ldoc_nde_t* ftr = doc->rt;
    
    /* Reimplement later where only used "things" are defined:
    ldoc_nde_t* ctx = ldoc_nde_new(LDOC_NDE_UA);
    ctx->mkup.anno.str = (char*)JSONLD_CTX;
     */
    
    ldoc_ent_t* ctx = ldoc_ent_new(LDOC_ENT_OR);
    ctx->pld.pair.anno.str = (char*)JSONLD_CTX;
    ctx->pld.pair.dtm.str = (char*)JSONLD_GFF3;

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

    ldoc_ent_t* st = ldoc_ent_new(LDOC_ENT_OR);
    st->pld.pair.anno.str = (char*)GFF_C4;
    st->pld.pair.dtm.str = coff[3];
    
    ldoc_ent_t* en = ldoc_ent_new(LDOC_ENT_OR);
    en->pld.pair.anno.str = (char*)GFF_C5;
    en->pld.pair.dtm.str = coff[4];

    ldoc_ent_t* scr = ldoc_ent_new(LDOC_ENT_OR);
    scr->pld.pair.anno.str = (char*)GFF_C6;
    scr->pld.pair.dtm.str = coff[5];

    ldoc_ent_t* strnd = ldoc_ent_new(LDOC_ENT_OR);
    strnd->pld.pair.anno.str = (char*)GFF_C7;
    strnd->pld.pair.dtm.str = coff[6];

    ldoc_ent_t* ph = ldoc_ent_new(LDOC_ENT_OR);
    ph->pld.pair.anno.str = (char*)GFF_C8;
    ph->pld.pair.dtm.str = coff[7];

    gff_splt_attrs(ftr, attrs, coff[8]);
    
    // JSON-LD context:
    // This needs to be changed when the context is dynamically created.
    ldoc_nde_ent_push(ftr, ctx);

    // Source, type, score, strand, phase:
    ldoc_nde_ent_push(ftr, src);
    ldoc_nde_ent_push(ftr, tpe);
    ldoc_nde_ent_push(ftr, scr);
    ldoc_nde_ent_push(ftr, strnd);
    ldoc_nde_ent_push(ftr, ph);

    // Coordinates, assigned to a locus:
    ldoc_nde_ent_push(lc, st);
    ldoc_nde_ent_push(lc, en);
    
    ldoc_nde_ent_push(ftr, lm);
    ldoc_nde_dsc_push(ftr, lc);
    ldoc_nde_dsc_push(ftr, attrs);

    return doc;
}

void gff_proc_fa(char* ln, size_t lnlen)
{
    
}

void gff_proc_cmt(char* ln, size_t lnlen)
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
    
    // "Step back" a character, so that the code for finding identifiers
    // is universal below (no special treatment of any lines; character
    // added due to the `fa--` should be a newline):
    fa--;
    
    if (fa != BI_NFOUND)
    {
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

static inline void gff_proc_ln(char* ln, size_t lnlen, gen_prsr_t* st)
{
    ldoc_doc_t* doc = NULL;
    
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
            }
            else
            {
                // Comment:
            }
        }
        else if (st->fa_sct)
        {
            
        }
        else
        {
            // Feature:
            doc = gff_proc_ftr(ln, lnlen);
        }
    }
    
    if (doc)
    {
        ldoc_vis_nde_ord_t* vis_nde = ldoc_vis_nde_ord_new();
        vis_nde->vis_setup = ldoc_vis_setup_json;
        vis_nde->vis_teardown = ldoc_vis_teardown_json;
        ldoc_vis_nde_uni(&(vis_nde->pre), ldoc_vis_nde_pre_json);
        ldoc_vis_nde_uni(&(vis_nde->infx), ldoc_vis_nde_infx_json);
        ldoc_vis_nde_uni(&(vis_nde->post), ldoc_vis_nde_post_json);
        
        ldoc_vis_ent_t* vis_ent = ldoc_vis_ent_new();
        ldoc_vis_ent_uni(vis_ent, ldoc_vis_ent_json);

        ldoc_ser_t* ser = ldoc_format(doc, vis_nde, vis_ent);
        
        printf("%s\n", ser->sclr.str);
    }
}

char* gff_rd_ln(fio_mem* mem, off_t mx, size_t llen, char* ln, size_t* ln_len, off_t off)
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

void gff_rd(int fd, off_t mx)
{
    off_t off = 0;
    fio_mem* mem = NULL;
    int incr = 0;
    
    gen_prsr_t st;
    st.fa_sct = false;
    
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
            mem_cpy = gff_rd_ln(mem, mx, lnlen, mem_cpy, &mem_len, off);
            
            gff_proc_ln(mem_cpy, lnlen, &st);
            
            // Next line:
            off += lnlen;
            ln_no++;
            lnlen = fio_lnlen(mem, off);
            
            //printf("%lu: %lu\n", ln_no, lnlen);
            //printf("%s", mem_cpy);
        } while (lnlen);
    }
    
    if (mem)
        fio_munmap(mem);
}

char* gff_seq(int fd, off_t mx, ldoc_trie_t* idx, const char* id, off_t st, off_t en, bool rv)
{
    ldoc_trie_nde_t* fa = ldoc_trie_lookup(idx, id, false);
    
    // ID could not be found:
    if (!fa)
        return NULL;
    
    // TODO Check annotation type.
    
    const gen_fa_ntry_t* ntry = fa->anno.pld;
    const size_t pg_size = getpagesize();
    const off_t seq_pg = ntry->seq / pg_size;
    
    off_t seq_pg_off = ntry->seq - (seq_pg * pg_size);
    
    fio_mem* mem = NULL;
    mem = fio_mmap(mem, fd, mx, (en / pg_size + 1) * pg_size, seq_pg * pg_size);
    
    // TODO Error check.
    
    size_t slen = en - st + 1;
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
    
    return seq;
}
