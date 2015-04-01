//
//  fio.c
//  biointerchange
//
//  Created by Joachim Baran on 2015-03-30.
//
//

#include "fio.h"

int fio_opn(const char* path)
{
    int fd = open(path, O_RDONLY | O_NONBLOCK | O_SYMLINK);
    
    if (fd == -1)
    {
        // TODO Error.
    }
    
    return fd;
}

void fio_cls(int fd)
{
    close(fd);
}

fio_mem* fio_mmap(fio_mem* mem, int fd, size_t mx, size_t len, off_t off)
{
    if (!mem)
    {
        // Completely new memory mapping:
        mem = (fio_mem*)malloc(sizeof(fio_mem));
        
        if (!mem)
        {
            // TODO Error.
        }
        
        mem->fd = fd;
        mem->mx = mx;
    }
    else
    {
        // Change page size/offset of an existing memory mapping:
        munmap(mem->pg, mem->ln);
        
        // TODO Error checking.
    }
    
    // Adjust len, in case it exceeds the maximum size:
    if (len > mem->mx)
        len = mem->mx;
    
    mem->pg = mmap(0, len, PROT_READ, MAP_FILE | MAP_SHARED, mem->fd, off);

    if (mem->pg == MAP_FAILED)
    {
        // TODO Error.
    }
    
    mem->ln = len;
    mem->off = off;
    
    return mem;
}

void fio_munmap(fio_mem* m)
{
    // TODO Error checking.
    munmap(m->pg, m->ln);
}

inline char* fio_rd(fio_mem* mem, size_t len, off_t off)
{
    // Check conditions under which acquiring a new page mapping becomes necessary:
    if (off < mem->off ||
        off >= mem->off + mem->ln ||
        off + len >= mem->off + mem->ln)
    {
        fio_mmap(mem, 0, 0, len, off);
    }
    
    return (char*)(mem->pg + off);
}

inline size_t fio_lnlen(fio_mem* mem, off_t off)
{
    off_t cur = off;
    
    while (cur < mem->ln)
    {
        char* ptr = (char*)(mem->pg + cur);
        
        if (*ptr == '\n' || *ptr == '\r')
            return cur - off;
        else if (cur + 1 == mem->ln && mem->off + mem->ln == mem->mx)
            return mem->ln - cur - 1; // TODO Verify.
        
        cur++;
    }
    
    return 0;
}
