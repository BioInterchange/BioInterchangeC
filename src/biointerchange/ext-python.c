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

#include "ext-python.h"

const char* test_python_fn = "multiply";
const char* test_python_impl = "simplepy.simple";

void python_ftr()
{
    PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue;
    int i;
    
    Py_Initialize();
    
    pName = PyUnicode_FromString(test_python_impl);
    /* Error checking of pName left out */
    
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);
    
    if (pModule != NULL)
    {
        pFunc = PyObject_GetAttrString(pModule, test_python_fn);
        /* pFunc is a new reference */
        
        if (pFunc && PyCallable_Check(pFunc))
        {
            pArgs = PyTuple_New(2);
            for (i = 0; i < 2; ++i)
            {
                pValue = PyLong_FromLong((i + 3) * 2);
                if (!pValue)
                {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    
                    return 1;
                }
                /* pValue reference stolen here: */
                PyTuple_SetItem(pArgs, i, pValue);
            }
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            if (pValue != NULL)
            {
                printf("Result of call: %ld\n", PyLong_AsLong(pValue));
                Py_DECREF(pValue);
            }
            else
            {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else
        {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", test_python_fn);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else
    {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", "TODO");
        
        return 1;
    }
    
    Py_Finalize();
    
    return 0;
}