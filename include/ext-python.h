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

#ifndef biointerchange_python_h
#define biointerchange_python_h

#include <Python.h>

#include <document.h>

#ifdef __cplusplus
extern "C" {
#endif
    
void py_init(char* pypath);

void py_free();
    
ldoc_doc_t* py_setup(ldoc_doc_t* ctxt, ldoc_doc_t* meta);
ldoc_doc_t* py_cleanup(ldoc_doc_t* stats);
ldoc_doc_t* py_process(ldoc_doc_t* ftr);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
