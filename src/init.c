#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP kdecopula_eval_beta(SEXP, SEXP, SEXP);
extern SEXP kdecopula_eval_cdf(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kdecopula_eval_hfunc_2d(SEXP, SEXP, SEXP, SEXP);
extern SEXP kdecopula_eval_mr(SEXP, SEXP, SEXP);
extern SEXP kdecopula_eval_t(SEXP, SEXP, SEXP);
extern SEXP kdecopula_interp(SEXP, SEXP, SEXP, SEXP);
extern SEXP kdecopula_interp_2d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kdecopula_inv_hfunc(SEXP, SEXP, SEXP, SEXP);
extern SEXP kdecopula_kern_epan(SEXP, SEXP);
extern SEXP kdecopula_kern_gauss_2d(SEXP, SEXP, SEXP);
extern SEXP kdecopula_renorm(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"kdecopula_eval_beta",     (DL_FUNC) &kdecopula_eval_beta,     3},
    {"kdecopula_eval_cdf",      (DL_FUNC) &kdecopula_eval_cdf,      5},
    {"kdecopula_eval_hfunc_2d", (DL_FUNC) &kdecopula_eval_hfunc_2d, 4},
    {"kdecopula_eval_mr",       (DL_FUNC) &kdecopula_eval_mr,       3},
    {"kdecopula_eval_t",        (DL_FUNC) &kdecopula_eval_t,        3},
    {"kdecopula_interp",        (DL_FUNC) &kdecopula_interp,        4},
    {"kdecopula_interp_2d",     (DL_FUNC) &kdecopula_interp_2d,     5},
    {"kdecopula_inv_hfunc",     (DL_FUNC) &kdecopula_inv_hfunc,     4},
    {"kdecopula_kern_epan",     (DL_FUNC) &kdecopula_kern_epan,     2},
    {"kdecopula_kern_gauss_2d", (DL_FUNC) &kdecopula_kern_gauss_2d, 3},
    {"kdecopula_renorm",        (DL_FUNC) &kdecopula_renorm,        4},
    {NULL, NULL, 0}
};

void R_init_kdecopula(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
