#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(cindex)(int *n,  int *nevent,
                     double *time, int *event, int *group,
                     double *w, double *cind, double *cwind, double *cnew, double *cwnew);

static const R_FortranMethodDef FortranEntries[] = {
    {"cindex",   (DL_FUNC) &F77_NAME(cindex),    10},
    {NULL, NULL, 0}
};

void R_init_BivRec(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
