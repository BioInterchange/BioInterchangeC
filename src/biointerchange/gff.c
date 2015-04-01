//
//  gff.c
//  biointerchange
//
//  Created by Joachim Baran on 2015-03-27.
//
//

#include "gff.h"

const char* GFF_FA = "\n##FASTA";

inline void gff_ky(char* attr, char** val)
{
    while (*attr) {
        if (*attr == '=')
        {
            *attr = 0;
            *val = ++attr;
            
            return;
        }
    }
    
    // No value assignment; just a key.
    *val = NULL;
}

void gff_splt_attrs(char* attrs)
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
    while (*attrs)
    {
        if (*attrs == ';')
        {
            // Isolate attribute assignment by making it an isolated string:
            *attrs = 0;
            
            // Isolate key/value:
            char* val;
            gff_ky(attr, &val);
            
            
            
            // Next attribute will have to start here:
            attr = attrs + 1;
        }
        
        attrs++;
    }
}

void gff_rd_prgm(char* ln)
{

}

void gff_rd_ftr(char* ln)
{
}

void gff_rd_fa(char* ln)
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

/**
 * Finds the first position of `str` in `txt`.
 *
 * str needs to be a proper prefix. For example, 'hello' is okay, 'hehello' is not.
 *
 * Returns -1, if the search string was not found.
 */
inline off_t gff_ffind(char* txt, const int len, const char* str, gen_fstat* stat, off_t stat_skip)
{
    off_t off = 0;
    int idx = 0;
    
    bool last_nl = false;
    bool last_hash = false;
    
    while (idx < len)
    {
        if (off >= stat_skip)
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
                stat->ftrs++;
            }
            else if (last_hash && *txt == '#')
            {
                last_nl = false;
                last_hash = false;
                stat->meta++;
            }
            else if (last_hash)
            {
                last_nl = false;
                last_hash = false;
                stat->comms++;
            }
        }
        
        off++;
        
        if (*txt == str[idx])
        {
            idx++;
            if (!str[idx])
            {
                // Update valid statistics offset:
                stat->off = off;
                
                // Match at this position:
                return txt - idx + 1;
            }
        }
        else
            idx = 0;
    }
    
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
            skip = pg_size - 2;
        // Commence the actual search:
        off_t fa_off = gff_ffind(mem->pg, len, GFF_FA, stat, skip);
        
        if (fa_off != BI_NFOUND)
            return fa_off;
        
        fa += len - pg_size;
    }
    
    if (mem)
        fio_munmap(mem);
    
    // Not found:
    return BI_NFOUND;
}

void* gff_crt_fa_idx(fio_mem* mem, off_t mx, gen_fstat* stat)
{
    // Note: assumes feature indices are at most 63 bytes long.
    char** idx = (char**)malloc(stat->ftrs * BI_FTR_ID_MAX);
    
    return NULL;
}

void gff_idx_fa(fio_mem* mem, off_t mx)
{
    gen_fstat* stat = (gen_fstat*)calloc(1, sizeof(gen_fstat));
    
    off_t fa = gff_fnd_fa(mem, stat, mx);
    
    if (fa != BI_NFOUND)
    {
        gff_crt_fa_idx(mem, mx, stat);
    }
}

void gff_rd(int fd, off_t mx)
{
    off_t off = 0;
    fio_mem* mem = NULL;
    int incr = 0;
    
    while (!mem || mem->off + mem->ln < mx)
    {
    gff_rd_incr_mem:
        mem = fio_mmap(NULL, fd, mx, getpagesize() * (BI_GEN_PG_MUL + incr), off);
        
        // TODO Error checking.
        
        // Figure out if one line can be read (based on current pointer:
        size_t lnlen = fio_lnlen(mem, mem->off);
        
        // Handle case where no line ending is visible:
        if (!lnlen && mem->off + mem->ln < mx)
        {
            incr++;
            goto gff_rd_incr_mem;
        }
        
        // Still no line ending visible? Then read to the end of the buffer (this is implicit; check fio_mmap behavior):
        if (!lnlen)
            lnlen = mem->mx - mem->off;
        
        char* mem_cpy = (char*)malloc(lnlen);
        
        if (!mem_cpy)
        {
            // TODO Error handling.
        }
        
        memcpy(mem_cpy, mem->pg, lnlen);
        
        // TODO DO SOMETHING
        
        // Advance offset for next iteration:
        mem->off += mem->ln;
    }
    
    if (mem)
        fio_munmap(mem);
}
