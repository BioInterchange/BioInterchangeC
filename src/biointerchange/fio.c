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

#include "fio.h"

int fio_opn(const char* path)
{
#ifdef __APPLE__
    int fd = open(path, O_RDONLY | O_NONBLOCK | O_SYMLINK);
#else
    int fd = open(path, O_RDONLY | O_NONBLOCK);
#endif
    
    if (fd == -1)
        gen_err(MAIN_ERR_FLE, path);
    
    return fd;
}

void fio_cls(int fd)
{
    close(fd);
}

size_t fio_len(int fd)
{
    return (size_t)lseek(fd, 0, SEEK_END);
}

fio_mem* fio_mmap(fio_mem* mem, int fd, size_t mx, size_t len, off_t off)
{
    if (!mem)
    {
        // Completely new memory mapping:
        mem = (fio_mem*)malloc(sizeof(fio_mem));
        
        if (!mem)
            gen_err(MAIN_ERR_SYSMALL, "File/memory mapping (new).");
        
        mem->fd = fd;
        mem->mx = mx;
    }
    else
    {
        // Change page size/offset of an existing memory mapping:
        msync(mem->pg, mem->ln, MS_INVALIDATE);
        munmap(mem->pg, mem->ln);
        
        // TODO: Error checking. -- not clear how.
    }
    
    // Adjust len, in case it exceeds the maximum size:
    if (off + len > mem->mx)
        len = mem->mx - off;
    
    mem->pg = mmap(0, len, PROT_READ, MAP_FILE | MAP_SHARED, mem->fd, off);

    if (mem->pg == MAP_FAILED)
        gen_err(MAIN_ERR_SYSMMAP, "File/memory mapping.");
    
    mem->ln = len;
    mem->off = off;
    
    return mem;
}

void fio_munmap(fio_mem* mem)
{
    // TODO Error checking. -- not clear what to check for
    msync(mem->pg, mem->ln, MS_INVALIDATE);
    munmap(mem->pg, mem->ln);
}

inline char* fio_rd(fio_mem* mem, size_t len, off_t off)
{
    // Check conditions under which acquiring a new page mapping becomes necessary:
    if (off < mem->off ||
        off >= mem->off + mem->ln ||
        off + len >= mem->off + mem->ln)
    {
        // TODO Error handling. -- check again
        mem = fio_mmap(mem, 0, 0, len, off);
    }
    
    return (char*)(mem->pg + (off - mem->off));
}

inline size_t fio_lnlen(fio_mem* mem, off_t off)
{
    bool lend = false;
    off_t cur = off;
    
    while (cur < mem->off + mem->ln)
    {
        char* ptr = (char*)(mem->pg + (cur - mem->off));
        
        if (lend && !(*ptr == '\n' || *ptr == '\r'))
            return cur - off;
        if (*ptr == '\n' || *ptr == '\r')
            lend = true;
        else if (cur + 1 >= mem->mx)
            return mem->mx - off;
        
        cur++;
    }
    
    if (cur == mem->mx)
        return cur - off;
    
    return 0;
}
