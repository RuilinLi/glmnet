#ifndef GLMNET_MATRIX
#define GLMNET_MATRIX
#include <cstdlib>

#include "pgenlib_ffi_support.h"
#include "pgenlib_read.h"
#include "pvar_ffi_support.h"

class MatrixGlmnet {
   public:
    // Compute the inner product of X[,j] and v
    virtual double dot_product(int j, const double* v) = 0;
    // Compute the inner product of X[,i] and X[,j]
    virtual double column_product(int i, int j) = 0;
    // Compute the inner product of X[,j]^2 and v
    virtual double vx2(int j, const double* v) = 0;

    // Set r = r - d*v*x[,j]
    virtual void update_res(int j, double d, const double* v,
                            double* r) = 0;

    static double sumv(const double* v, int len);

   protected:
    int no;  // Number of rows
    int ni;  // Number of variables
};

class DenseM : public MatrixGlmnet {
   public:
    DenseM(int no, int ni, const double* x);
    ~DenseM();

    double dot_product(int j, const double* v);

    double column_product(int i, int j);

    double vx2(int j, const double* v);

    void update_res(int j, double d, const double* v, double* r);

   private:
    const double* data;
};

class PlinkMatrix : public MatrixGlmnet{
   public:
    void Close();
    PlinkMatrix();

    ~PlinkMatrix();
    void Load(const char* fname, int raw_sample_ct, int* sample_subset,
              const uint32_t subset_size);
    void ReadCompact(int* variant_subset, const uintptr_t vsubset_size);

    uintptr_t** compactM;
    double dot_product(int j, const double* v);
    double vx2(int j, const double* v);
    double column_product(int i, int j);
    void update_res(int j, double d, const double* v, double* r);

   private:
    plink2::PgenFileInfo* _info_ptr;
    plink2::RefcountedWptr* _allele_idx_offsetsp;
    plink2::RefcountedWptr* _nonref_flagsp;
    plink2::PgenReader* _state_ptr;
    uintptr_t* _subset_include_vec;
    uintptr_t* _subset_include_interleaved_vec;
    uint32_t* _subset_cumulative_popcounts;
    plink2::PgrSampleSubsetIndex _subset_index;
    uint32_t _subset_size;
    bool malloc_all;
    uintptr_t _vsubset_size;

    plink2::PgenVariant _pgv;

    plink2::VecW* _transpose_batch_buf;
    // kPglNypTransposeBatch (= 256) variants at a time, and then transpose
    uintptr_t* _multivar_vmaj_geno_buf;
    uintptr_t* _multivar_vmaj_phasepresent_buf;
    uintptr_t* _multivar_vmaj_phaseinfo_buf;
    uintptr_t* _multivar_smaj_geno_batch_buf;
    uintptr_t* _multivar_smaj_phaseinfo_batch_buf;
    uintptr_t* _multivar_smaj_phasepresent_batch_buf;

    void SetSampleSubsetInternal(int* sample_subset,
                                 const uint32_t subset_size);
};

#endif