//
//  fio.h
//  biointerchange
//
//  Created by Joachim Baran on 2015-03-30.
//
//

#ifndef biointerchange_fio_h
#define biointerchange_fio_h

#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>

typedef struct fio_mem
{
    int fd;
    size_t mx;
    void* pg;
    size_t ln;
    off_t off;
} fio_mem;

int fio_opn(const char* path);
void fio_cls(int fd);

fio_mem* fio_mmap(fio_mem* mem, int fd, size_t mx, size_t len, off_t off);
void fio_munmap(fio_mem* m);

size_t fio_lnlen(fio_mem* mem, off_t off);

#endif
