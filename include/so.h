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

#ifndef biointerchange_so_h
#define biointerchange_so_h

#include <document.h>

#ifdef __cplusplus
extern "C" {
#endif

extern const char* ONT[];
    
ldoc_trie_t* ont_trie_new();
    
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
