#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
// extern void FreePELT(void *);
// extern void PELTC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void segneigh (char**, double*, int*, int*, double*, int*, int*, int*, int*, int*, double*);

static const R_CMethodDef CEntries[] = {
    // {"FreePELT", (DL_FUNC) &FreePELT,  1},
    // {"PELTC",     (DL_FUNC) &PELTC,     11},
    {"segneigh",     (DL_FUNC) &segneigh,     11},
    {NULL, NULL, 0}
};

void R_init_changepoint(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

