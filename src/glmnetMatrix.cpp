#include "glmnetMatrix.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>

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
    return std::inner_product(data + j * no, data + (j + 1) * no, v, 0.0);
}

double DenseM::column_product(int i, int j) {
    return std::inner_product(data + i * no, data + (i + 1) * no, data + j * no,
                              0.0);
}

double DenseM::vx2(int j, const double *v) {
    double result = 0.0;
#pragma omp simd reduction(+ : result)
    for (int i = 0; i < no; ++i) {
        result += data[j * no + i] * data[j * no + i] * v[i];
    }
    return result;
}

void DenseM::update_res(int j, double d, const double *v,
                        double *__restrict r) {
    for (int i = 0; i < no; ++i) {
        r[i] -= d * v[i] * data[j * no + i];
    }
}
