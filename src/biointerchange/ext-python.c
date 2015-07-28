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

#include "ext-python.h"

const char* test_python_fn = "multiply";
const char* test_python_impl = "simplepy.simple";

const char* py_fn_setup = "setup";
const char* py_fn_cleanup = "cleanup";
const char* py_fn_process = "process_feature";

static PyObject* ext_mod;
static PyObject* ext_setup;
static PyObject* ext_cleanup;
static PyObject* ext_process;

static inline ldoc_ser_t* ldoc2py(ldoc_doc_t* doc)
{
    ldoc_vis_nde_ord_t* vis_nde = ldoc_vis_nde_ord_new();
    vis_nde->vis_setup = ldoc_vis_setup_py;
    vis_nde->vis_teardown = ldoc_vis_teardown_py;
    ldoc_vis_nde_uni(&(vis_nde->pre), ldoc_vis_nde_pre_py);
    ldoc_vis_nde_uni(&(vis_nde->infx), ldoc_vis_nde_infx_py);
    ldoc_vis_nde_uni(&(vis_nde->post), ldoc_vis_nde_post_py);
    
    ldoc_vis_ent_t* vis_ent = ldoc_vis_ent_new();
    ldoc_vis_ent_uni(vis_ent, ldoc_vis_ent_py);
    
    ldoc_ser_t* ser = ldoc_format(doc, vis_nde, vis_ent);
    
    ldoc_vis_nde_ord_free(vis_nde);
    ldoc_vis_ent_free(vis_ent);
    
    return ser;
}

static void py_fn(PyObject* mod, char* fn, PyObject** fn_ptr)
{
    PyObject* py_fn = PyObject_GetAttrString(mod, fn);
    
    if (py_fn && PyCallable_Check(py_fn))
    {
        *fn_ptr = py_fn;
    }
    else
    {
        if (PyErr_Occurred())
            PyErr_Print();
        fprintf(stderr, "Cannot find function \"%s\"\n", test_python_fn);
    }
}

void py_init(char* pypath)
{
    PyObject* name;
    
    Py_Initialize();
    
    wchar_t* bname[] = { L"BioInterchange" };
    PySys_SetArgv(1, bname);
    
    name = PyUnicode_FromString(pypath);
    ext_mod = PyImport_Import(name);
    Py_DECREF(name);
    
    if (ext_mod)
    {
        py_fn(ext_mod, (char*)py_fn_setup, &ext_setup);
        py_fn(ext_mod, (char*)py_fn_cleanup, &ext_cleanup);
        py_fn(ext_mod, (char*)py_fn_process, &ext_process);
    }
    else
    {
        PyErr_Print();
        
        gen_err(MAIN_ERR_PYMD, (const char*)pypath);
    }
}

void py_free()
{
    Py_DECREF(ext_setup);
    Py_DECREF(ext_cleanup);
    Py_DECREF(ext_process);
    Py_DECREF(ext_mod);
    
    Py_Finalize();
}

static inline ldoc_doc_t* py_call(PyObject* fn, ldoc_doc_t* d1, ldoc_doc_t* d2)
{
    ldoc_ser_t* ser1 = ldoc2py(d1);
    
    if (!ser1)
        gen_err(MAIN_ERR_SYSMALL, "Python Dict creation for function call (1).");
    
    PyObject* py_d1 = ser1->pld.py.dtm;
    
    PyObject* dict;
    if (!d2)
        dict = PyObject_CallFunctionObjArgs(fn, py_d1, NULL);
    else
    {
        ldoc_ser_t* ser2 = ldoc2py(d2);
        
        if (!ser2)
            gen_err(MAIN_ERR_SYSMALL, "Python Dict creation for function call (2).");
        
        PyObject* py_d2 = ser2->pld.py.dtm;
        
        dict = PyObject_CallFunctionObjArgs(fn, py_d1, py_d2, NULL);
        
        Py_DECREF(py_d2);
        free(ser2);
    }
    
    // TODO Quite frankly, I do not know why Py_DECREF cannot be called
    //      all the time. Checking for None returned works though. Must
    //      be something strange going on when reference counting when
    //      py_d1 is set to None in the external code.
    if (dict != Py_None)
        Py_DECREF(py_d1);
    free(ser1);
    
    if (!dict)
        gen_err(MAIN_ERR_PYMEM, "Function call returned a garbled result.");
    
    ldoc_doc_t* ret = NULL;
    if (dict != Py_None)
    {
        ret = ldoc_pydict2doc(dict);
    
        Py_DECREF(dict);
    }
    
    // Note: DO NOT call the garbage collector by force! Garbage collection
    //       takes an incredible amount of time in Python! It is the worst!
    //PyGC_Collect();
    
    return ret;
}

ldoc_doc_t* py_setup(ldoc_doc_t* ctxt, ldoc_doc_t* meta)
{
    return py_call(ext_setup, ctxt, meta);
}

ldoc_doc_t* py_cleanup(ldoc_doc_t* stats)
{
    return py_call(ext_cleanup, stats, NULL);
}

ldoc_doc_t* py_process(ldoc_doc_t* ftr)
{
    return py_call(ext_process, ftr, NULL);
}