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

// Created via: cat biointerchange-l2.txt | tr $'\n' '@' | sed 's/@/\\n/g' | pbcopy
static const char* eula_txt = "BioInterchange Software License\n===============================\n\nTerms and conditions\n--------------------\n\n1.  Preamble: This Agreement governs the relationship between persons or legal entities (hereinafter: Licensee) executing the BioInterchange software version 2.0.0 and later versions, and CODAMONO whose principal place of business is Ontario, Canada (Hereinafter: Licensor). This Agreement sets the terms, rights, restrictions and obligations on using the BioInterchange software (hereinafter: The Software) created and owned by Licensor, as detailed herein.\n\n2.  License Grant: Licensor grants Licensee a non-assignable & non-transferable, commercial, royalty free, without the rights to create derivative works, non-exclusive license, all with accordance with the terms set forth and other legal restrictions set forth in 3rd party software used while running The Software.\n\n    a.  Limited: Licensee may use The Software for the purpose of:\n\n        i.   Running The Software; and\n\n        ii.  Use The Software's output; and\n\n        iii. Distribute verbatim copies of The Software's output.\n\n    b.  Binary Limited:\n\n        i.  Licensee may not distribute nor sublicense The Software as a part of a larger work; and\n\n        ii. Licensee may not reverse engineer, decompile, decode, decrypt, disassemble, or in any way derive source code from The Software.\n\n    c.  Non Assignable & Non-Transferable: Licensee may not assign or transfer his rights and duties under this license.\n\n    d.  Commercial, Royalty Free: Licensee may use The Software for any purpose, including paid-services, without any royalties.\n\n3.  License Restrictions: Licensor restricts Licensee's use rights of The Software per acquired license type.\n\n    a. Trial License: Licensor grants Licensee a personal, non-renewable license, with a termination date set 30 days after Licensee's agreement to the terms and conditions herein.\n\n    b. Reviewer License: Licensor grants Licensee a personal license for the sole purpose of aiding the academic peer-review process, with a termination date set 30 days after Licensee's agreement to the terms and conditions herein.\n\n    c. Individual License: Licensor grants Licensee a personal license for use in personal, non-collaborative, projects in agreement to the terms and conditions herein.\n\n    d. Team License: Licensor grants Licensee the right to simultaneously use The Software by up to five colocated individuals, who are affiliated with the Licensee, in agreement to the terms and conditions herein.\n\n    e. Cloud, Cluster & HIPAA License: Licensor grants Licensee the right to run The Software, within the Licensee's organization, in agreement to the terms and conditions herein.\n\n4.  Term & Termination: The Term of this license shall be until terminated. Licensor may terminate this Agreement, including Licensee's license in the case where Licensee:\n\n    a.  Became insolvent or otherwise entered into any liquidation process; or\n\n    b.  Exported The Software to any jurisdiction where licensor may not enforce his rights under this agreements in; or\n\n    c.  Licensee was in breach of any of this license's terms and conditions and such breach was not cured, immediately upon notification; or\n\n    d.  Licensee in breach of any of the terms of clause 2 or clause 3 to this license; or\n\n    e.  Licensee otherwise entered into any arrangement which caused Licensor to be unable to enforce his rights under this License.\n\n5.  Payment: In consideration of the License granted under clause 2 and clause 3, Licensee shall pay Licensor a FEE. Failure to perform payment shall construe as material breach of this Agreement.\n\n6.  Upgrades, Updates and Fixes: Licensor may provide Licensee, from time to time, with Upgrades, Updates or Fixes, as detailed herein and according to his sole discretion. Licensee warrants to keep The Software up-to-date and will install all relevant updates and fixes, and may, at his sole discretion, purchase upgrades, according to the rates set by Licensor. Licensor shall provide any update or Fix free of charge; however, nothing in this Agreement shall require Licensor to provide Updates or Fixes.\n\n    a.  Upgrades: for the purpose of this license, an Upgrade shall be a material amendment in The Software, which contains new features and or major performance improvements and shall be marked as a new major version number (an increment of X in the version number format X.Y.Z).\n\n    b.  Updates: for the purpose of this license, an update shall be a minor amendment in The Software, which may contain new features or minor improvements and shall be marked as a new minor number (an increment of Y in the version number format X.Y.Z).\n\n    c.  Fix: for the purpose of this license, a fix shall be a minor amendment in The Software, intended to remove bugs or alter minor features which impair The Software's functionality. A fix shall be marked as a new patch number (an increment of Z in the version number format X.Y.Z).\n\n7.  Support: The Software is provided under an AS-IS basis and without any support, updates or maintenance. Nothing in this Agreement shall require Licensor to provide Licensee with support or fixes to any bug, failure, malperformance or other defect in The Software.\n\n    a.  Bug Notification: Licensee may provide Licensor of details regarding any bug, defect or failure in The Software promptly and with no delay from such event; Licensee shall comply with Licensor's request for information regarding bugs, defects or failures and furnish him with information, terminal (shell, command line) output and try to reproduce such bugs, defects or failures.\n\n    b.  Feature Request: Licensee may request additional features in The Software, provided, however, that (i) Licensee shall waive any claim or right in such feature should feature be developed by Licensor; (ii) Licensee shall be prohibited from developing the feature, or disclose such feature request, or feature, to any 3rd party directly competing with Licensor or any 3rd party which may be, following the development of such feature, in direct competition with Licensor; (iii) Licensee warrants that feature does not infringe any 3rd party patent, trademark, trade-secret or any other intellectual property right; and (iv) Licensee developed, envisioned or created the feature solely by himself.\n\n8.  Liability: To the extent permitted under Law, The Software is provided under an as is basis. Licensor shall never, and without any limit, be liable for any damage, cost, expense or any other payment incurred by Licensee as a result of The Software's actions, failure, bugs and/or any other interaction between The Software  and Licensee's end-equipment, computers, other software or any 3rd party, end-equipment, computer or services. Moreover, Licensor shall never be liable for any defect in source code written by Licensee when relying on The Software or using The Software's source code.\n\n9.  Warranty:\n\n    a.  Intellectual Property: Licensor warrants that The Software does not violate or infringe any 3rd party claims in regards to intellectual property, patents and/or trademarks and that to the best of its knowledge no legal action has been taken against it for any infringement or violation of any 3rd party intellectual property rights.\n\n    b.  No-Warranty: The Software is provided without any warranty; Licensor disclaims any warranty that The Software shall be error free, without defects or code that may cause damage to Licensee's computers or to Licensee, and that Software shall be functional. Licensee shall be solely liable to any damage, defect or loss incurred as a result of operating software and undertake the risks contained in running The Software.\n\n    c.  Prior Inspection: Licensee states that he inspected The Software thoroughly and found it satisfactory and adequate to his needs, that it does not interfere with his regular operation and that it does meet the standards and scope of his computer systems and architecture. Licensee found that The Software interacts with his software environment and that it does not infringe any of End User License Agreement of any software Licensee may use in performing his services. Licensee waives any claims regarding The Software's incompatibility, performance, results and features, and warrants that he inspected The Software.\n\n10.  No Refunds: Licensee warrants that he inspected The Software according to clause 9.c and that it is adequate to his needs. Accordingly, as The Software is intangible goods, Licensee shall not be, ever, entitled to any refund, rebate, compensation or restitution for any reason whatsoever, even if The Software contains material flaws.\n\n11. Indemnification: Licensee warrants to hold Licensor harmless and indemnify Licensor for any lawsuit brought against it in regards to Licensee's use of The Software in means that violate, breach or otherwise circumvent this license, Licensor's intellectual property rights or Licensor's title in The Software. Licensor shall promptly notify Licensee in case of such legal action and request Licensee's consent prior to any settlement in relation to such lawsuit or claim.\n\n12. Governing Law, Jurisdiction: Licensee agrees not to initiate class-action lawsuits against Licensor in relation to this license. Licensee agrees to compensate Licensor for any legal fees, cost or attorney fees should any claim brought by Licensee against Licensor be denied, in part or in full.\n\n";

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
    fprintf(stderr, "%s", eula_txt);
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
    
    // Allocate quick memory heap (20MB):
    if (!qk_alloc(20*1024*1024))
    {
        fprintf(stderr, "Sorry, there is not enough memory available to fire up the software.\n\n");
        
        exit(MAIN_ERR_SYSMALL);
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
    while ((c = getopt(argc, argv, "o:p:u:vhe")) != -1)
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
        fprintf(stderr, "one of the following extensions: .gff .gff3 .gvf .ldj .ldjson .vcf\n");
        
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
        fprintf(stderr, "Known filename extensions: .gff .gff3 .gvf .ldj .ldjson .vcf\n");
        
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
#ifdef __APPLE__
    int lfd = open(licpath, O_RDONLY | O_NONBLOCK | O_SYMLINK);
#else
    int lfd = open(licpath, O_RDONLY | O_NONBLOCK);
#endif
    
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
            idx = gff_idx_fa(fd, &stat, mx);
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
            gen_err(MAIN_ERR_LICN, NULL);
        default:
            gen_err(MAIN_ERR_LICV, NULL);
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
