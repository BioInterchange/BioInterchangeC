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

#include <gtest/gtest.h>

#include <unistd.h>

#include "fio.h"

#define FIO_GFF_SMALL "../test-data/chromosome_BF.gff"
#define FIO_GFF_SMALL_HD "##gff-version   3\n"
#define FIO_GFF_SMALL_NPAGE "ontig GI N"
#define FIO_GFF_SMALL_PAGE6 "rent=DDB02"

TEST(fio, access_testdata)
{
    int fd = fio_opn(FIO_GFF_SMALL);
    
    EXPECT_NE(-1, fd);
    
    size_t mx = fio_len(fd);
    
    EXPECT_EQ(149958, mx);
    
    fio_cls(fd);
}

TEST(fio, mmap)
{
    int fd = fio_opn(FIO_GFF_SMALL);
    size_t mx = fio_len(fd);
    
    fio_mem* mem = fio_mmap(NULL, fd, mx, getpagesize(), 0);
    EXPECT_NE((fio_mem*)NULL, mem);
    
    fio_munmap(mem);
    
    fio_cls(fd);
}

TEST(fio, readline)
{
    int fd = fio_opn(FIO_GFF_SMALL);
    size_t mx = fio_len(fd);
    fio_mem* mem = fio_mmap(NULL, fd, mx, getpagesize(), 0);

    char* hd = fio_rd(mem, strlen(FIO_GFF_SMALL_HD), 0);
    char* hd_chp = strndup(hd, strlen(FIO_GFF_SMALL_HD));
    EXPECT_STREQ(FIO_GFF_SMALL_HD, hd_chp);
    free(hd_chp);
    
    fio_munmap(mem);
    fio_cls(fd);
}

TEST(fio, readline_mmap_adjust)
{
    int fd = fio_opn(FIO_GFF_SMALL);
    size_t mx = fio_len(fd);
    fio_mem* mem = fio_mmap(NULL, fd, mx, getpagesize(), 0);
    
    // Map less than a page size:
    size_t map_len = 10;
    EXPECT_LT(map_len, getpagesize());
    
    char* s = fio_rd(mem, map_len, getpagesize());
    char* s_chp = strndup(s, map_len);
    EXPECT_STREQ(FIO_GFF_SMALL_NPAGE, s_chp);
    free(s_chp);

    // Map a page further in:
    s = fio_rd(mem, map_len, getpagesize() * 5);
    s_chp = strndup(s, map_len);
    EXPECT_STREQ(FIO_GFF_SMALL_PAGE6, s_chp);
    free(s_chp);
    
    fio_munmap(mem);
    fio_cls(fd);
}

TEST(fio, lnlen)
{
    int fd = fio_opn(FIO_GFF_SMALL);
    size_t mx = fio_len(fd);
    fio_mem* mem = fio_mmap(NULL, fd, mx, getpagesize(), 0);

    size_t llen = fio_lnlen(mem, 0);
    EXPECT_EQ(strlen(FIO_GFF_SMALL_HD), llen);
    
    char* s = fio_rd(mem, getpagesize(), getpagesize());
    char* s_chp = strndup(s, getpagesize());
    char* s_ = s_chp;
    while (*s_ != '\n')
        s_++;
    *s_ = 0;
    llen = fio_lnlen(mem, getpagesize());
    EXPECT_EQ(strlen(s_chp) + 1, llen);
    free(s_chp);
    
    s = fio_rd(mem, getpagesize(), getpagesize() * 5);
    s_chp = strndup(s, getpagesize());
    s_ = s_chp;
    while (*s_ != '\n')
        s_++;
    *s_ = 0;
    llen = fio_lnlen(mem, getpagesize() * 5);
    EXPECT_EQ(strlen(s_chp) + 1, llen);
    free(s_chp);
    
    fio_munmap(mem);
    fio_cls(fd);
}
