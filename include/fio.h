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

#ifndef biointerchange_fio_h
#define biointerchange_fio_h

#include <fcntl.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>

#include <errno.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif
    
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

size_t fio_len(int fd);

char* fio_rd(fio_mem* mem, size_t len, off_t off);
size_t fio_lnlen(fio_mem* mem, off_t off);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
