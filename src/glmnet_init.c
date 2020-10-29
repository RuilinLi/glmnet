#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP Getju(SEXP, SEXP);
extern SEXP initialize_plinkmatrix_Xptr(SEXP, SEXP, SEXP);
extern SEXP PlinkMatrix_info(SEXP, SEXP);
extern SEXP PlinkMultiplyv(SEXP, SEXP, SEXP);
extern SEXP PlinkPreMultiplyv(SEXP, SEXP, SEXP);
extern SEXP PlinkSetMean(SEXP, SEXP);
extern SEXP PlinkSetSd(SEXP, SEXP);
extern SEXP storePB(SEXP);
extern SEXP wls_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP wls_plink(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(chg_big)(void *);
extern void F77_NAME(chg_bnorm)(void *, void *);
extern void F77_NAME(chg_dev_max)(void *);
extern void F77_NAME(chg_epsnr)(void *);
extern void F77_NAME(chg_fract_dev)(void *);
extern void F77_NAME(chg_itrace)(void *);
extern void F77_NAME(chg_max_exp)(void *);
extern void F77_NAME(chg_min_flmin)(void *);
extern void F77_NAME(chg_min_lambdas)(void *);
extern void F77_NAME(chg_min_null_prob)(void *);
extern void F77_NAME(chg_mxitnr)(void *);
extern void F77_NAME(coxnet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(elnet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(fishnet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(get_bnorm)(void *, void *);
extern void F77_NAME(get_int_parms)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(get_int_parms2)(void *, void *);
extern void F77_NAME(loglike)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lognet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(multelnet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(multspelnet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(spelnet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(spfishnet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(splognet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(spwls)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"Getju",                       (DL_FUNC) &Getju,                        2},
    {"initialize_plinkmatrix_Xptr", (DL_FUNC) &initialize_plinkmatrix_Xptr,  3},
    {"PlinkMatrix_info",            (DL_FUNC) &PlinkMatrix_info,             2},
    {"PlinkMultiplyv",              (DL_FUNC) &PlinkMultiplyv,               3},
    {"PlinkPreMultiplyv",           (DL_FUNC) &PlinkPreMultiplyv,            3},
    {"PlinkSetMean",                (DL_FUNC) &PlinkSetMean,                 2},
    {"PlinkSetSd",                  (DL_FUNC) &PlinkSetSd,                   2},
    {"storePB",                     (DL_FUNC) &storePB,                      1},
    {"wls_dense",                   (DL_FUNC) &wls_dense,                   27},
    {"wls_plink",                   (DL_FUNC) &wls_plink,                   27},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"chg_big",           (DL_FUNC) &F77_NAME(chg_big),            1},
    {"chg_bnorm",         (DL_FUNC) &F77_NAME(chg_bnorm),          2},
    {"chg_dev_max",       (DL_FUNC) &F77_NAME(chg_dev_max),        1},
    {"chg_epsnr",         (DL_FUNC) &F77_NAME(chg_epsnr),          1},
    {"chg_fract_dev",     (DL_FUNC) &F77_NAME(chg_fract_dev),      1},
    {"chg_itrace",        (DL_FUNC) &F77_NAME(chg_itrace),         1},
    {"chg_max_exp",       (DL_FUNC) &F77_NAME(chg_max_exp),        1},
    {"chg_min_flmin",     (DL_FUNC) &F77_NAME(chg_min_flmin),      1},
    {"chg_min_lambdas",   (DL_FUNC) &F77_NAME(chg_min_lambdas),    1},
    {"chg_min_null_prob", (DL_FUNC) &F77_NAME(chg_min_null_prob),  1},
    {"chg_mxitnr",        (DL_FUNC) &F77_NAME(chg_mxitnr),         1},
    {"coxnet",            (DL_FUNC) &F77_NAME(coxnet),            28},
    {"elnet",             (DL_FUNC) &F77_NAME(elnet),             28},
    {"fishnet",           (DL_FUNC) &F77_NAME(fishnet),           29},
    {"get_bnorm",         (DL_FUNC) &F77_NAME(get_bnorm),          2},
    {"get_int_parms",     (DL_FUNC) &F77_NAME(get_int_parms),      8},
    {"get_int_parms2",    (DL_FUNC) &F77_NAME(get_int_parms2),     2},
    {"loglike",           (DL_FUNC) &F77_NAME(loglike),           11},
    {"lognet",            (DL_FUNC) &F77_NAME(lognet),            30},
    {"multelnet",         (DL_FUNC) &F77_NAME(multelnet),         29},
    {"multspelnet",       (DL_FUNC) &F77_NAME(multspelnet),       31},
    {"spelnet",           (DL_FUNC) &F77_NAME(spelnet),           30},
    {"spfishnet",         (DL_FUNC) &F77_NAME(spfishnet),         31},
    {"splognet",          (DL_FUNC) &F77_NAME(splognet),          32},
    {"spwls",             (DL_FUNC) &F77_NAME(spwls),             31},
    {NULL, NULL, 0}
};

void R_init_glmnet(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
