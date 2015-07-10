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
#include "ext-python.h"
#include "gen.h"
#include "fio.h"
#include "gff.h"
#include "gvf.h"
#include "vcf.h"

#define MIN_PAGESIZE 2048

// Relative to the home directory:
#define MAIN_LICPATH "/.biointerchange/biointerchange-license"

// Check size of int without consulting limit.h:
// http://www.pixelbeat.org/programming/gcc/static_assert.html
#define ASSERT_CONCAT_(a, b) a##b
#define ASSERT_CONCAT(a, b) ASSERT_CONCAT_(a, b)
#define ct_assert(e) enum { ASSERT_CONCAT(assert_line_, __LINE__) = 1/(!!(e)) }
ct_assert(sizeof(int)>=4);

void usage(char* bname, char* version, int err)
{
    fprintf(stderr, "Usage: %s [parameters] genomicsfile\n\n", bname);
    fprintf(stderr, "Description:\n    Converts GFF3, GVF and VCF files to JSON/JSON-LD.\n");
    fprintf(stderr, "    Outputs one JSON document per line.\n");
    fprintf(stderr, "Parameters:\n");
    fprintf(stderr, "    -o file           : writes output into file (default: STDOUT)\n");
    fprintf(stderr, "    -p package.module : Python API called on package.module\n");
    fprintf(stderr, "    -e                : prints EULA and exits\n");
    fprintf(stderr, "    -v                : prints version number and exits\n");
    fprintf(stderr, "    -h                : this help text\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "BioInterchange %s by CODAMONO\n", version);
    
    exit(err);
}

void eula(char* version)
{
    printf("BioInterchange %s by CODAMONO\n", version);
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
    
    // Allocate quick memory heap (10MB):
    if (!qk_alloc(10*1024*1024))
    {
        fprintf(stderr, "Sorry, there is not enough memory available to fire up the software.\n\n");
        
        exit(MAIN_ERR_PGSZ); // TODO
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
    
    // Execution context -- set defaults:
    gen_ctxt_t ctxt;
    ctxt.fgen = NULL;
    ctxt.fname = NULL;
    ctxt.fout = stdout;
    ctxt.py = false;
    ctxt.pycall = NULL;
    ctxt.usr = NULL;
    ctxt.ver = BIOINTERCHANGE_VERSION;
    
    // Parameter check 1: use getopt to fetch optional parameters:
    int c;
    char* py_opt = NULL;
    
    int err = MAIN_ERR_PARA;
    while ((c = getopt(argc, argv, "o:p:u:vh")) != -1)
    {
        switch (c)
        {
            case 'o':
                ctxt.fname = strdup(optarg);
                break;
            case 'p':
                py_opt = strdup(optarg);
                ctxt.py = true;
                ctxt.pycall = py_opt; // TODO Cleanup.
                break;
            case 'u':
                ctxt.usr = strdup(optarg);
                break;
            case 'e':
                eula(BIOINTERCHANGE_VERSION);
                exit(MAIN_SUCCESS);
            case 'v':
                printf("%s\n", BIOINTERCHANGE_VERSION);
                exit(MAIN_SUCCESS);
            case 'h':
                err = MAIN_SUCCESS;
            default:
                usage(bname, BIOINTERCHANGE_VERSION, err);
        }
    }
    argc -= optind;
    argv += optind;
    
    // Parameter check 2; do not check license before the software is called correctly:
    if (argc != 1)
        usage(bname, BIOINTERCHANGE_VERSION, MAIN_ERR_PARA);
    
    char* fname = argv[0];
    ctxt.fgen = fname; // TODO Cleanup.
    size_t fname_len = strlen(fname);
    if (fname_len < 5)
    {
        fprintf(stderr, "Filename too short. Needs to be at least one character followed by\n");
        fprintf(stderr, "one of the following extensions: .gff .gvf .vcf\n");
        
        exit(MAIN_ERR_FNME);
    }
    
    // Determine filetype by extension:
    char* ext_s = &fname[fname_len - 4];
    char* ext_l = &fname[fname_len - 5]; // Still works with fname_len above, but means that the basename could be an empty string.
    char* ext_xl = fname_len > 7 ? &fname[fname_len - 7] : NULL;
    if (!strcmp(ext_s, ".gff") || !strcmp(ext_l, ".gff3"))
        ctxt.tpe = GEN_FMT_GFF3;
    // else if (!strcmp(ext_s, ".gtf"))
    //    ctxt.tpe = GEN_FMT_GTF;
    else if (!strcmp(ext_s, ".gvf"))
        ctxt.tpe = GEN_FMT_GVF;
    else if (!strcmp(ext_s, ".vcf"))
        ctxt.tpe = GEN_FMT_VCF;
    else if (!strcmp(ext_s, ".ldj") || (ext_xl && !strcmp(ext_xl, ".ldjson")))
        ctxt.tpe = GEN_FMT_LDJ;
    else
    {
        fprintf(stderr, "Cannot determine filetype by filename extension.\n");
        fprintf(stderr, "Known filename extensions: .gff .gvf .ldj .ldjson .vcf\n");
        
        exit(MAIN_ERR_FEXT);
    }
    
    // See if the file can be opened; still not check the license:
    int fd = fio_opn(fname);
    
    if (fd == -1)
    {
        fprintf(stderr, "Cannot access: %s\n", fname);
        
        exit(MAIN_ERR_FACC);
    }
    
    // See if the output should go to stdout, or, a file:
    if (ctxt.fname)
    {
        ctxt.fout = fopen(ctxt.fname, "w");
    }
    
    // Initialize Python -- if needed; not checking license yet:
    if (py_opt)
    {
        py_init(py_opt);
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
    
    gen_init();
    
    // Index file first to get stats:
    gen_fstat stat = { 0, 0, 0, 0 };
    size_t mx = fio_len(fd);
    
    ldoc_trie_t* idx = NULL;
    switch (ctxt.tpe)
    {
        case GEN_FMT_GFF3:
        case GEN_FMT_GVF:
            gvf_init();
            gff_idx_fa(fd, &stat, mx);
            break;
        default:
            // Do nothing.
            break;
    }

    // Check the license now:
    lic_status_t status = lic_valid(lid, &stat);
    
    free(lid);
    
    switch (status)
    {
        case LICENSE_OK:
            break;
        case LICENSE_NET:
            exit(MAIN_ERR_LICN);
        default:
            exit(MAIN_ERR_LICV);
    }
    
    // Set up callback functions based on the file type:
    gen_cbcks_t cbcks;
    switch (ctxt.tpe)
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
        case GEN_FMT_LDJ:
            // Do nothing. Callbacks determined by '@context' in data.
            break;
        default:
            // TODO Internal error.
            break;
    }
    
    if (ctxt.tpe == GEN_FMT_LDJ)
        gen_rd_doc(fd, mx, &ctxt);
    else
        gen_rd(fd, mx, idx, &cbcks, &ctxt);
    
    fio_cls(fd);
    
    qk_free();
    
    if (py_opt)
        py_free();
    
    // Wrap up output -- do not close it if no file name given (stdout used in that case):
    fflush(ctxt.fout);
    if (ctxt.fname)
        fclose(ctxt.fout);
    
    return MAIN_SUCCESS;
}
