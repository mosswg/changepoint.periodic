#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
// extern void FreePELT(void *);
// extern void PELTC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sncirc(void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*);

static const R_CMethodDef CEntries[] = {
    // {"FreePELT", (DL_FUNC) &FreePELT,  1},
    // {"PELTC",     (DL_FUNC) &PELTC,     11},
    {"sncirc",     (DL_FUNC) &sncirc,     20},
    {NULL, NULL, 0}
};

void R_init(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

