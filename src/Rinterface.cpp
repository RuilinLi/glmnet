
#include "R.h"
#include "R_ext/Print.h"
#include "Rinternals.h"
#include "glmnetMatrix.h"

void wls_base(double alm0, double almc, double alpha, int m, int no, int ni,
              MatrixGlmnet *X, double *r, const double *v, int intr,
              const int *ju, const double *vp, const double *cl, int nx,
              double thr, int maxit, double *__restrict a, double *aint,
              double *__restrict g, int *__restrict ia, int *__restrict iy,
              int *iz, int *__restrict mm, int *nino, double *rsqc, int *nlp,
              int *jerr);

static void finalizer(SEXP xptr) {
    PlinkMatrix *p = (PlinkMatrix *)R_ExternalPtrAddr(xptr);
    delete p;
}

class Wls_Solver {
   protected:
    double alm0;
    double almc;
    double alpha;
    int m;
    int no;
    int ni;
    double *r;
    double *v;
    int intr;
    int *ju;
    double *vp;
    double *cl;
    int nx;
    double thr;
    int maxit;
    double *a;
    double *aint;
    double *g;
    int *ia;
    int *iy;
    int *iz;
    int *mm;
    int *nino;
    double *rsqc;
    int *nlp;
    int *jerr;
    SEXP result;

   public:
    Wls_Solver(SEXP alm02, SEXP almc2, SEXP alpha2, SEXP m2, SEXP nobs2,
               SEXP nvars2, SEXP x2, SEXP r2, SEXP v2, SEXP intr2, SEXP ju2,
               SEXP vp2, SEXP cl2, SEXP nx2, SEXP thr2, SEXP maxit2, SEXP a2,
               SEXP aint2, SEXP g2, SEXP ia2, SEXP iy2, SEXP iz2, SEXP mm2,
               SEXP nino2, SEXP rsqc2, SEXP nlp2, SEXP jerr2) {
        alm0 = asReal(alm02);
        almc = asReal(almc2);
        alpha = asReal(alpha2);
        m = asInteger(m2);
        no = asInteger(nobs2);
        ni = asInteger(nvars2);
        r = REAL(r2);
        v = REAL(v2);
        intr = asInteger(intr2);
        ju = INTEGER(ju2);
        vp = REAL(vp2);
        cl = REAL(cl2);
        nx = asInteger(nx2);
        thr = asReal(thr2);
        maxit = asInteger(maxit2);
        a = REAL(a2);
        aint = REAL(aint2);
        g = REAL(g2);
        ia = INTEGER(ia2);
        iy = INTEGER(iy2);
        iz = INTEGER(iz2);
        mm = INTEGER(mm2);
        nino = INTEGER(nino2);
        rsqc = REAL(rsqc2);
        nlp = INTEGER(nlp2);
        jerr = INTEGER(jerr2);

        // Setup the result;

        result = PROTECT(allocVector(VECSXP, 20));

        SET_VECTOR_ELT(result, 0, almc2);
        SET_VECTOR_ELT(result, 1, m2);
        SET_VECTOR_ELT(result, 2, nobs2);
        SET_VECTOR_ELT(result, 3, nvars2);
        SET_VECTOR_ELT(result, 4, r2);
        SET_VECTOR_ELT(result, 5, ju2);
        SET_VECTOR_ELT(result, 6, vp2);
        SET_VECTOR_ELT(result, 7, cl2);
        SET_VECTOR_ELT(result, 8, nx2);
        SET_VECTOR_ELT(result, 9, a2);
        SET_VECTOR_ELT(result, 10, aint2);
        SET_VECTOR_ELT(result, 11, g2);
        SET_VECTOR_ELT(result, 12, ia2);
        SET_VECTOR_ELT(result, 13, iy2);
        SET_VECTOR_ELT(result, 14, iz2);
        SET_VECTOR_ELT(result, 15, mm2);
        SET_VECTOR_ELT(result, 16, nino2);
        SET_VECTOR_ELT(result, 17, rsqc2);
        SET_VECTOR_ELT(result, 18, nlp2);
        SET_VECTOR_ELT(result, 19, jerr2);
    }

    SEXP get_result() {
        // Should only be called once
        UNPROTECT(1);
        return result;
    }

    virtual void solve() = 0;
};

class Wls_Solver_Dense : public Wls_Solver {
   private:
    DenseM X;

   public:
    Wls_Solver_Dense(SEXP alm02, SEXP almc2, SEXP alpha2, SEXP m2, SEXP nobs2,
                     SEXP nvars2, SEXP x2, SEXP r2, SEXP v2, SEXP intr2,
                     SEXP ju2, SEXP vp2, SEXP cl2, SEXP nx2, SEXP thr2,
                     SEXP maxit2, SEXP a2, SEXP aint2, SEXP g2, SEXP ia2,
                     SEXP iy2, SEXP iz2, SEXP mm2, SEXP nino2, SEXP rsqc2,
                     SEXP nlp2, SEXP jerr2)
        : Wls_Solver(alm02, almc2, alpha2, m2, nobs2, nvars2, x2, r2, v2, intr2,
                     ju2, vp2, cl2, nx2, thr2, maxit2, a2, aint2, g2, ia2, iy2,
                     iz2, mm2, nino2, rsqc2, nlp2, jerr2),
          X(no, ni, REAL(x2)) {}

    void solve() {
        wls_base(alm0, almc, alpha, m, no, ni, &X, r, v, intr, ju, vp, cl, nx,
                 thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr);
    }
};

class Wls_Solver_Plink : public Wls_Solver {
   private:
    PlinkMatrix *X;
    bool already_solved;

   public:
    Wls_Solver_Plink(SEXP alm02, SEXP almc2, SEXP alpha2, SEXP m2, SEXP nobs2,
                     SEXP nvars2, SEXP x2, SEXP r2, SEXP v2, SEXP intr2,
                     SEXP ju2, SEXP vp2, SEXP cl2, SEXP nx2, SEXP thr2,
                     SEXP maxit2, SEXP a2, SEXP aint2, SEXP g2, SEXP ia2,
                     SEXP iy2, SEXP iz2, SEXP mm2, SEXP nino2, SEXP rsqc2,
                     SEXP nlp2, SEXP jerr2)
        : Wls_Solver(alm02, almc2, alpha2, m2, nobs2, nvars2, x2, r2, v2, intr2,
                     ju2, vp2, cl2, nx2, thr2, maxit2, a2, aint2, g2, ia2, iy2,
                     iz2, mm2, nino2, rsqc2, nlp2, jerr2) {
        X = (PlinkMatrix *)R_ExternalPtrAddr(x2);
        already_solved = false;
    }

    void solve() {
        if (already_solved) {
            // Do nothing
            return;
        }

        wls_base(alm0, almc, alpha, m, no, ni, X, r, v, intr, ju, vp, cl, nx,
                 thr, maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr);
        already_solved = true;
    }
};

#ifdef __cplusplus
extern "C" {
#endif
SEXP initialize_plinkmatrix_Xptr(SEXP fname2, SEXP sample_subset2,
                                 SEXP vsubset2) {
    const char *fname = CHAR(STRING_ELT(fname2, 0));
    int *sample_subset = INTEGER(sample_subset2);
    const uint32_t subset_size = length(sample_subset2);
    int *vsubset = INTEGER(vsubset2);
    const uintptr_t vsubset_size = length(vsubset2);
    PlinkMatrix *p = new PlinkMatrix();
    p->Load(fname, UINT32_MAX, sample_subset, subset_size);
    p->ReadCompact(vsubset, vsubset_size);
    SEXP xptr = R_MakeExternalPtr(p, R_NilValue, R_NilValue);
    PROTECT(xptr);
    R_RegisterCFinalizerEx(xptr, finalizer, TRUE);
    UNPROTECT(1);
    return xptr;
}

SEXP PlinkMatrix_info(SEXP ptr2, SEXP weight2) {
    PlinkMatrix *p = (PlinkMatrix *)R_ExternalPtrAddr(ptr2);
    const double *weight = REAL(weight2);
    const uint32_t subset_size = length(weight2);
    const uintptr_t vsubset_size = p->get_ni();
    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SEXP xm = PROTECT(allocVector(REALSXP, vsubset_size));
    SEXP xs = PROTECT(allocVector(REALSXP, vsubset_size));
    SEXP xmax = PROTECT(allocVector(REALSXP, vsubset_size));
    SEXP xmin = PROTECT(allocVector(REALSXP, vsubset_size));

    for (int j = 0; j < vsubset_size; ++j) {
        double rbuf[4];
        p->get_info(j, weight, subset_size, rbuf);
        REAL(xm)[j] = rbuf[0];
        REAL(xs)[j] = rbuf[1] - rbuf[0] * rbuf[0];  // weighted variance
        REAL(xmax)[j] = rbuf[2];
        REAL(xmin)[j] = rbuf[3];
    }

    SET_VECTOR_ELT(result, 0, xm);
    SET_VECTOR_ELT(result, 1, xs);
    SET_VECTOR_ELT(result, 2, xmax);
    SET_VECTOR_ELT(result, 3, xmin);

    UNPROTECT(5);
    return result;
}

SEXP wls_plink(SEXP alm02, SEXP almc2, SEXP alpha2, SEXP m2, SEXP nobs2,
               SEXP nvars2, SEXP x2, SEXP r2, SEXP v2, SEXP intr2, SEXP ju2,
               SEXP vp2, SEXP cl2, SEXP nx2, SEXP thr2, SEXP maxit2, SEXP a2,
               SEXP aint2, SEXP g2, SEXP ia2, SEXP iy2, SEXP iz2, SEXP mm2,
               SEXP nino2, SEXP rsqc2, SEXP nlp2, SEXP jerr2) {
    Wls_Solver_Plink Solver(alm02, almc2, alpha2, m2, nobs2, nvars2, x2, r2, v2,
                            intr2, ju2, vp2, cl2, nx2, thr2, maxit2, a2, aint2,
                            g2, ia2, iy2, iz2, mm2, nino2, rsqc2, nlp2, jerr2);

    Solver.solve();
    return Solver.get_result();
}

SEXP wls_dense(SEXP alm02, SEXP almc2, SEXP alpha2, SEXP m2, SEXP nobs2,
               SEXP nvars2, SEXP x2, SEXP r2, SEXP v2, SEXP intr2, SEXP ju2,
               SEXP vp2, SEXP cl2, SEXP nx2, SEXP thr2, SEXP maxit2, SEXP a2,
               SEXP aint2, SEXP g2, SEXP ia2, SEXP iy2, SEXP iz2, SEXP mm2,
               SEXP nino2, SEXP rsqc2, SEXP nlp2, SEXP jerr2) {
    Wls_Solver_Dense Solver(alm02, almc2, alpha2, m2, nobs2, nvars2, x2, r2, v2,
                            intr2, ju2, vp2, cl2, nx2, thr2, maxit2, a2, aint2,
                            g2, ia2, iy2, iz2, mm2, nino2, rsqc2, nlp2, jerr2);

    Solver.solve();
    return Solver.get_result();
}

SEXP PlinkSetMean(SEXP ptr2, SEXP xm2) {
    PlinkMatrix *p = (PlinkMatrix *)R_ExternalPtrAddr(ptr2);
    p->setxm(REAL(xm2));
    Rprintf("The last xm is %f", p->xm[199]);
    return R_NilValue;
}

SEXP PlinkSetSd(SEXP ptr2, SEXP xs2) {
    PlinkMatrix *p = (PlinkMatrix *)R_ExternalPtrAddr(ptr2);
    p->setxs(REAL(xs2));
    Rprintf("The last xs is %f", p->xs[199]);
    return R_NilValue;
}

SEXP PlinkMultiplyv(SEXP ptr2, SEXP v2, SEXP r2) {
    PlinkMatrix *p = (PlinkMatrix *)R_ExternalPtrAddr(ptr2);
    p->multiply_vector(REAL(v2), REAL(r2));
    return R_NilValue;
}

SEXP PlinkPreMultiplyv(SEXP ptr2, SEXP v2, SEXP r2) {
    PlinkMatrix *p = (PlinkMatrix *)R_ExternalPtrAddr(ptr2);
    p->pre_multiply_vector(REAL(v2), REAL(r2));
    return R_NilValue;
}
// SEXP wls_dense(SEXP alm02, SEXP almc2, SEXP alpha2, SEXP m2, SEXP nobs2,
//                SEXP nvars2, SEXP x2, SEXP r2, SEXP v2, SEXP intr2, SEXP ju2,
//                SEXP vp2, SEXP cl2, SEXP nx2, SEXP thr2, SEXP maxit2, SEXP a2,
//                SEXP aint2, SEXP g2, SEXP ia2, SEXP iy2, SEXP iz2, SEXP mm2,
//                SEXP nino2, SEXP rsqc2, SEXP nlp2, SEXP jerr2) {
//     double alm0 = asReal(alm02);
//     double almc = asReal(almc2);
//     double alpha = asReal(alpha2);
//     int m = asInteger(m2);
//     int no = asInteger(nobs2);
//     int ni = asInteger(nvars2);
//     DenseM X(no, ni, REAL(x2));
//     double *r = REAL(r2);
//     double *v = REAL(v2);
//     int intr = asInteger(intr2);
//     int *ju = INTEGER(ju2);
//     double *vp = REAL(vp2);
//     double *cl = REAL(cl2);
//     int nx = asInteger(nx2);
//     double thr = asReal(thr2);
//     int maxit = asInteger(maxit2);
//     double *a = REAL(a2);
//     double *aint = REAL(aint2);
//     double *g = REAL(g2);
//     int *ia = INTEGER(ia2);
//     int *iy = INTEGER(iy2);
//     int *iz = INTEGER(iz2);
//     int *mm = INTEGER(mm2);
//     int *nino = INTEGER(nino2);
//     double *rsqc = REAL(rsqc2);
//     int *nlp = INTEGER(nlp2);
//     int *jerr = INTEGER(jerr2);

//     wls_base(alm0, almc, alpha, m, no, ni, &X, r, v, intr, ju, vp, cl, nx,
//     thr,
//              maxit, a, aint, g, ia, iy, iz, mm, nino, rsqc, nlp, jerr);

//     SEXP result;
//     result = PROTECT(allocVector(VECSXP, 1));

//     SET_VECTOR_ELT(result, 0, jerr2);

//     UNPROTECT(1);
//     // Setup the result;
//     // static const char *outnames[] = {
//     //     "almc", "m", "no", "ni", "r",  "ju", "vp",   "cl",   "nx",  "a",
//     //     "aint", "g", "ia", "iy", "iz", "mm", "nino", "rsqc", "nlp",
//     "jerr"};

//     // static const char *outnames[] = {"almc"};
//     // SEXP result;
//     // result = PROTECT(mkNamed(VECSXP, outnames));
//     // SET_VECTOR_ELT(result, 0, almc2);
//     // SET_VECTOR_ELT(result, 0, almc2);
//     // SET_VECTOR_ELT(result, 1, m2);
//     // SET_VECTOR_ELT(result, 2, nobs2);
//     // SET_VECTOR_ELT(result, 3, nvars2);
//     // SET_VECTOR_ELT(result, 4, r2);
//     // SET_VECTOR_ELT(result, 5, ju2);
//     // SET_VECTOR_ELT(result, 6, vp2);
//     // SET_VECTOR_ELT(result, 7, cl2);
//     // SET_VECTOR_ELT(result, 8, nx2);
//     // SET_VECTOR_ELT(result, 9, a2);
//     // SET_VECTOR_ELT(result, 10, aint2);
//     // SET_VECTOR_ELT(result, 11, g2);
//     // SET_VECTOR_ELT(result, 12, ia2);
//     // SET_VECTOR_ELT(result, 13, iy2);
//     // SET_VECTOR_ELT(result, 14, iz2);
//     // SET_VECTOR_ELT(result, 15, mm2);
//     // SET_VECTOR_ELT(result, 16, nino2);
//     // SET_VECTOR_ELT(result, 17, rsqc2);
//     // SET_VECTOR_ELT(result, 18, nlp2);
//     // SET_VECTOR_ELT(result, 19, jerr2);
//     // UNPROTECT(1);

//     return (result);
// }

#ifdef __cplusplus
}
#endif