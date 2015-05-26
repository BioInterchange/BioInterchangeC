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

#include <libgen.h>
#include <pwd.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>

// External (license *will* be external eventually):
#include "document.h"
#include "license.h"

// Internal headers (part of BioInterchange):
#include "gen.h"
#include "fio.h"
#include "gff.h"
#include "gvf.h"
#include "vcf.h"

#define MAIN_SUCCESS  0

// Never change these definitions or it will be very hard to interpret
// error codes across software versions!
#define MAIN_ERR_PGSZ 1
#define MAIN_ERR_HME1 2
#define MAIN_ERR_HME2 3
#define MAIN_ERR_PARA 4
#define MAIN_ERR_LISZ 5
#define MAIN_ERR_LICF 6
#define MAIN_ERR_FNME 7
#define MAIN_ERR_FEXT 8

#define MIN_PAGESIZE 2048

// Relative to the home directory:
#define MAIN_LICPATH "/.biointerchange-license"

// Check size of int without consulting limit.h:
// http://www.pixelbeat.org/programming/gcc/static_assert.html
#define ASSERT_CONCAT_(a, b) a##b
#define ASSERT_CONCAT(a, b) ASSERT_CONCAT_(a, b)
#define ct_assert(e) enum { ASSERT_CONCAT(assert_line_, __LINE__) = 1/(!!(e)) }
ct_assert(sizeof(int)>=4);

static void gff_cbcks(gen_cbcks_t* cbcks)
{
    cbcks->proc_ln = &gff_proc_ln;
}

static void gvf_cbcks(gen_cbcks_t* cbcks)
{
    cbcks->proc_ln = &gvf_proc_ln;
}

static void vcf_cbcks(gen_cbcks_t* cbcks)
{
    cbcks->proc_ln = &vcf_proc_ln;
}

int main(int argc, char* argv[])
{
    char* bname = basename(argv[0]);
    
    // Check whether this system has at least a pagesize of MIN_PAGESIZE:
    if (getpagesize() < MIN_PAGESIZE)
    {
        fprintf(stderr, "Sorry, this system's memory page-size is too small.\n\n");
        fprintf(stderr, "Minimum required page size: %ul\n", MIN_PAGESIZE);
        fprintf(stderr, "Page size of this system  : %ul\n", getpagesize());
        
        exit(MAIN_ERR_PGSZ);
    }
    
    // Get home directory:
    struct passwd *pws;
    pws = getpwuid(geteuid());
    
    if (!pws || !pws->pw_dir)
    {
        fprintf(stderr, "Cannot determine the path of the user's home directory.\n");
        
        if (!pws)
            exit(MAIN_ERR_HME1);
        else
            exit(MAIN_ERR_HME2);
    }
    
    // Get license file path:
    size_t licpath_len = strlen(MAIN_LICPATH) + strlen(pws->pw_dir);
    char* licpath = (char*)malloc(licpath_len + 1);
    snprintf(licpath, licpath_len + 1, "%s%s", pws->pw_dir, MAIN_LICPATH);
    
    // Parameter check; do not check license before the software is called correctly:
    if (argc < 2 || argc > 3)
    {
        fprintf(stderr, "Usage: %s genomicsfile [jsonfile]\n\n", bname);
        fprintf(stderr, "Description:\n    Converts GFF3, GVF and VCF files to JSON-LD.\n");
        fprintf(stderr, "    Output will be written to \"jsonfile\", if given, or printed\n");
        fprintf(stderr, "    on the console otherwise.\n");
        fprintf(stderr, "    The output consists of one JSON document per-line.\n");
        
        exit(MAIN_ERR_PARA);
    }
    
    char* fname = argv[1];
    size_t fname_len = strlen(fname);
    if (fname_len < 5)
    {
        fprintf(stderr, "Filename too short. Needs to be at least one character followed by\n");
        fprintf(stderr, "one of the following extensions: .gff .gtf .gvf .vcf\n");
        
        exit(MAIN_ERR_FNME);
    }
    
    // Determine filetype by extension:
    gen_filetype_t ftype;
    char* ext_s = &fname[fname_len - 4];
    char* ext_l = &fname[fname_len - 5]; // Still works with fname_len above, but means that the basename could be an empty string.
    if (!strcmp(ext_s, ".gff") || !strcmp(ext_l, ".gff3"))
        ftype = GEN_FMT_GFF3;
    else if (!strcmp(ext_s, ".gtf"))
        ftype = GEN_FMT_GTF;
    else if (!strcmp(ext_s, ".gvf"))
        ftype = GEN_FMT_GVF;
    else if (!strcmp(ext_s, ".vcf"))
        ftype = GEN_FMT_VCF;
    else
    {
        fprintf(stderr, "Cannot determine filetype by filename extension.\n");
        fprintf(stderr, "Known filename extensions: .gff .gtf .gvf .vcf\n");
        
        exit(MAIN_ERR_FEXT);
    }
    


    
    // See if the file can be opened; still not check the license:
    int fd = fio_opn(argv[1]);
    
    if (fd == -1)
    {
        // TODO
        exit(1234);
    }

    // Now, everything is in place; go and check whether the software is licensed:
    int lfd = open(licpath, O_RDONLY | O_NONBLOCK | O_SYMLINK);
    
    if (lfd == -1)
    {
        fprintf(stderr, "License file not found.\n\n");
        fprintf(stderr, "Expected license file path: %s\n", licpath);
        
        fio_cls(fd);
        exit(MAIN_ERR_LICF);
    }
    
    // Get license as string (and make sure it is a string -- not binary data):
    size_t asize = (size_t)lseek(lfd, 0, SEEK_END);
    size_t msize = getpagesize();
    if (asize > msize)
    {
        fprintf(stderr, "License file too big; this is not a valid license file.\n");
        
        close(lfd);
        fio_cls(fd);
        
        exit(MAIN_ERR_LISZ);
    }
    char* lid = (char*)malloc(asize + 1);
    char* lidp = lid;
    char* ltxt = mmap(0, msize, PROT_READ, MAP_FILE | MAP_SHARED, lfd, 0);
    char* vrfy = ltxt;
    while (vrfy < ltxt + asize)
    {
        char c = *(vrfy++);
        
        if ((c >= 'A' && c <= 'Z') ||
            (c >= 'a' && c <= 'z') ||
            (c >= '0' && c <= '9'))
        {
            // No problem: save character:
            *(lidp++) = c;
        }
        else if (c == '\n' || c == '\r')
        {
            // Bail out; everything was fine to this point and the string ends here.
            break;
        }
        else
        {
            // This is not good! License file is corrupt:
            fprintf(stderr, "License file contains invalid characters.\n");
            
            close(lfd);
            fio_cls(fd);
            
            exit(MAIN_ERR_LICF);
        }
    }
    *lidp = 0;
    
    munmap(ltxt, msize);
    
    // Index file first to get stats:
    gen_init();
    
    gen_fstat stat = { 0, 0, 0, 0 };
    size_t mx = fio_len(fd);
    
    ldoc_trie_t* idx = NULL;
    switch (ftype)
    {
        case GEN_FMT_GFF3:
        case GEN_FMT_GVF:
            gff_idx_fa(fd, &stat, mx);
            break;
        default:
            // Do nothing.
            break;
    }

    // Check the license now; send stats too:
    lic_status_t status = lic_valid(lid, &stat);
    
    free(lid);
    
    if (status == LICENSE_OK)
    {
        // Everything fine!
    }
    
    gen_cbcks_t cbcks;
    
    switch (ftype)
    {
        case GEN_FMT_GFF3:
            gff_cbcks(&cbcks);
            break;
        case GEN_FMT_GTF:
            // TODO
            gff_cbcks(&cbcks);
            break;
        case GEN_FMT_GVF:
            gvf_cbcks(&cbcks);
            break;
        case GEN_FMT_VCF:
            // TODO
            vcf_cbcks(&cbcks);
            break;
        default:
            // TODO Internal error.
            break;
    }
    
    gen_rd(fd, mx, idx, &cbcks);
    
    fio_cls(fd);
    
    return MAIN_SUCCESS;
}
