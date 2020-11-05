#include "glmnetMatrix.h"
// #include "mkl.h"
#include <chrono>
#include <math.h>
#ifndef __clang__
#if __GNUC__ < 8
#define _mm256_set_m128d(vh, vl) \
    _mm256_insertf128_pd(_mm256_castpd128_pd256(vl), (vh), 1)
#endif
#endif

#include <Eigen/Core>

MatrixGlmnet::~MatrixGlmnet(){}

double MatrixGlmnet::sumv(const double * v, int len) {
    double result = 0.0;
    for (int i = 0; i < len; ++i) {
        result += v[i];
    }
    return result;
}

DenseM::DenseM(int no, int ni, const double *x) {
    this->no = no;
    this->ni = ni;
    data = x;
}
DenseM::~DenseM() { data = nullptr; }

double DenseM::dot_product(int j, const double *v) {
    // return std::inner_product(data + j * no, data + (j + 1) * no, v, 0.0);

    // double result = 0.0;
    // // If there's no auto vectorization then we can do #pragma clang loop
    // // vectorize(enable) interleave(enable)
    // for (int i = 0; i < no; ++i) {
    //     result += data[j * no + i] * v[i];
    // }
    Eigen::Map<const Eigen::VectorXd> x(data + j * no, no);
    Eigen::Map<const Eigen::VectorXd> y(v, no);
    return x.dot(y);
}

double DenseM::vx2(int j, const double *v) {
    // double result = 0.0;
    // for (int i = 0; i < no; ++i) {
    //     result += data[j * no + i] * data[j * no + i] * v[i];
    // }
    Eigen::Map<const Eigen::ArrayXd> x(data + j * no, no);
    Eigen::Map<const Eigen::ArrayXd> y(v, no);
    return (y * x.square()).sum();
}
double DenseM::column_product(int i, int j) {
    // return std::inner_product(data + i * no, data + (i + 1) * no, data + j * no,
    //                         0.0);
    return 0.0;
}


void DenseM::update_res(int j, double d, const double *v,
                        double *__restrict r) {
    for (int i = 0; i < no; ++i) {
        r[i] -= d * v[i] * data[j * no + i];
    }
}
