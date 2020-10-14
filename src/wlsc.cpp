#include <cmath>
#include <cstdlib>

#include "R.h"
#include "R_ext/Print.h"
#include "Rinternals.h"
#include "glmnetMatrix.h"
void wls_base(double alm0, double almc, double alpha, int m, int no, int ni,
              MatrixGlmnet *X, double *__restrict r, const double *v, int intr,
              const int *ju, const double *vp, const double *cl, int nx,
              double thr, int maxit, double *__restrict a, double *aint,
              double *__restrict g, int *__restrict ia, int *__restrict iy,
              int *iz, int *__restrict mm, int *nino, double *rsqc, int *nlp,
              int *jerr) {
    double *__restrict xv = (double *)malloc(sizeof(double) * ni);
    double xmz = MatrixGlmnet::sumv(v, no);
    double ab = almc * alpha;
    double dem = almc * (1.0 - alpha);
    double tlam = alpha * (2.0 * almc - alm0);

    for (int j = 0; j < ni; ++j) {
        if (ju[j]) {
            g[j] = abs(X->dot_product(j, r));
        } else {
            continue;
        }

        if (iy[j]) {
            xv[j] = X->vx2(j, v);
        } else if (g[j] > tlam * vp[j]) {
            iy[j] = 1;
            xv[j] = X->vx2(j, v);
        }
    }

    bool jz = true;

    while (true) {
        if (!((*iz) && jz)) {
            (*nlp)++;
            double dlx = 0.0;
            for (int j = 0; j < ni; ++j) {
                if (!iy[j]) {
                    continue;
                }

                double gj = X->dot_product(j, r);
                double aj = a[j];
                double u = gj + aj * xv[j];
                double au = abs(u) - vp[j] * ab;
                if (au < 0.0) {
                    a[j] = 0.0;
                } else {
                    a[j] = fmax(cl[2 * j],
                                fmin(cl[2 * j + 1],
                                     copysign(au, u) / (xv[j] + vp[j] * dem)));
                }

                if (a[j] == aj) {
                    continue;
                }

                if (mm[j] == 0) {
                    (*nino)++;
                    if ((*nino) > nx) {
                        break;
                    }
                    mm[j] = (*nino);
                    ia[(*nino) - 1] = j;
                }
                double d = a[j] - aj;
                (*rsqc) += d * (2.0 * gj - d * xv[j]);
                X->update_res(j, d, v, r);
                dlx = fmax(xv[j] * d * d, dlx);
            }
            if ((*nino) > nx) {
                break;
            }
            if (intr) {
                double sumr = MatrixGlmnet::sumv(r, no);
                double d = sumr / xmz;
                (*aint) += d;
                (*rsqc) += d * (2.0 * sumr - d * xmz);

                dlx = fmax(dlx, xmz * d * d);

                for (int i = 0; i < no; ++i) {
                    r[i] -= d * v[i];
                }
            }

            // KKT checking here
            if (dlx < thr) {
                bool ixx = false;
                for (int j = 0; j < ni; ++j) {
                    if (iy[j] || (!ju[j])) {
                        continue;
                    }
                    g[j] = abs(X->dot_product(j, r));
                    if (g[j] > ab * vp[j]) {
                        iy[j] = 1;
                        xv[j] = X->vx2(j, v);
                        ixx = true;
                    }
                }

                if (ixx) {
                    continue;
                }
                break;
            }

            if ((*nlp) > maxit) {
                *jerr = -m;
                return;
            }
        }
        (*iz) = 1;

        while (true) {
            
            (*nlp)++;
            double dlx = 0.0;
            for (int l = 0; l < (*nino); ++l) {
                int k = ia[l];
                double gk = X->dot_product(k, r);
                double ak = a[k];
                double u = gk + ak * xv[k];
                double au = abs(u) - vp[k] * ab;
                if (au < 0.0) {
                    a[k] = 0.0;
                } else {
                    a[k] = fmax(cl[2 * k],
                                fmin(cl[2 * k + 1],
                                     copysign(au, u) / (xv[k] + vp[k] * dem)));
                }

                if (ak == a[k]) {
                    continue;
                }
                double d = a[k] - ak;
                (*rsqc) += d * (2.0 * gk - d * xv[k]);
                X->update_res(k, d, v, r);
                dlx = fmax(xv[k] * d * d, dlx);
            }

            if (intr) {
                double sumr = MatrixGlmnet::sumv(r, no);
                double d = sumr / xmz;
                (*aint) += d;
                (*rsqc) += d * (2.0 * sumr - d * xmz);

                dlx = fmax(dlx, xmz * d * d);

                for (int i = 0; i < no; ++i) {
                    r[i] -= d * v[i];
                }
            }

            if (dlx < thr) {
                break;
            }

            if ((*nlp) > maxit) {
                *jerr = -m;
                return;
            }
        }
        jz = false;
    }
    free(xv);
    return;
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

class Wls_Solver_plink {
   private:
    PlinkMatrix X;

   public:
    Wls_Solver_plink(SEXP x2)
           {
               // Now I expect x2 to be a list with
               // x2[0] is the name of the pgen file
               // x2[1] is an integer array specifying the sample subset
               // x2[2] is an integer array specifying the variant subset
               const char * fname = CHAR(VECTOR_ELT(x2, 0));
              Rprintf(fname);
           }

};

#ifdef __cplusplus
extern "C" {
#endif

SEXP wls_plink(SEXP x)
{
    Wls_Solver_plink a(x);
    return R_NilValue;
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