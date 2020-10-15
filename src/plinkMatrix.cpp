#include <bitset>
#include <iostream>
#include <stdexcept>

#include "glmnetMatrix.h"
#include "pgenlib_ffi_support.h"
#include "pgenlib_read.h"
#include "pvar_ffi_support.h"
// This is mostly adapted from the pgenlibr code
// https://github.com/chrchang/plink-ng/blob/master/2.0/pgenlibr/src/pgenlibr.cpp

static void stop(const char* msg) { throw std::runtime_error(msg); }

PlinkMatrix::PlinkMatrix()
    : _info_ptr(nullptr),
      _allele_idx_offsetsp(nullptr),
      _nonref_flagsp(nullptr),
      _state_ptr(nullptr) {
    malloc_all = false;
}

void PlinkMatrix::Load(const char* fname, int raw_sample_ct, int* sample_subset,
                       const uint32_t subset_size) {
    if (_info_ptr) {
        Close();
    }
    _info_ptr = static_cast<plink2::PgenFileInfo*>(
        malloc(sizeof(plink2::PgenFileInfo)));
    if (!_info_ptr) {
        stop("Out of memory");
    }
    no = (int)subset_size;
    plink2::PreinitPgfi(_info_ptr);
    uint32_t cur_sample_ct = UINT32_MAX;
    cur_sample_ct = raw_sample_ct;
    uint32_t cur_variant_ct = UINT32_MAX;

    plink2::PgenHeaderCtrl header_ctrl;
    uintptr_t pgfi_alloc_cacheline_ct;
    char errstr_buf[plink2::kPglErrstrBufBlen];
    if (PgfiInitPhase1(fname, cur_variant_ct, cur_sample_ct, 0, &header_ctrl,
                       _info_ptr, &pgfi_alloc_cacheline_ct,
                       errstr_buf) != plink2::kPglRetSuccess) {
        stop(&(errstr_buf[7]));
    }
    const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
    if (header_ctrl & 0x30) {
        // no need to zero-initialize this
        _allele_idx_offsetsp = plink2::CreateRefcountedWptr(raw_variant_ct + 1);
        _info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
        // _info_ptr->max_allele_ct updated by PgfiInitPhase2() in this case
    }
    _info_ptr->max_allele_ct = 2;
    if ((header_ctrl & 0xc0) == 0xc0) {
        // todo: load this in pvar, to enable consistency check.  we use a
        // (manually implemented) shared_ptr in preparation for this.
        const uintptr_t raw_variant_ctl =
            plink2::DivUp(raw_variant_ct, plink2::kBitsPerWord);
        // no need to zero-initialize this
        _nonref_flagsp = plink2::CreateRefcountedWptr(raw_variant_ctl + 1);
        _info_ptr->nonref_flags = _nonref_flagsp->p;
    }
    const uint32_t file_sample_ct = _info_ptr->raw_sample_ct;
    unsigned char* pgfi_alloc = nullptr;
    if (plink2::cachealigned_malloc(
            pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
        stop("Out of memory");
    }
    uint32_t max_vrec_width;
    uintptr_t pgr_alloc_cacheline_ct;
    if (PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width,
                       _info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct,
                       errstr_buf)) {
        if (pgfi_alloc && (!_info_ptr->vrtypes)) {
            plink2::aligned_free(pgfi_alloc);
        }
        stop(&(errstr_buf[7]));
    }
    if ((!_allele_idx_offsetsp) && (_info_ptr->gflags & 4)) {
        // Note that it's safe to be ignorant of multiallelic variants when
        // phase and dosage info aren't present; GetAlleleCt() then always
        // returns 2 when that isn't actually true, and all ALTs are treated as
        // if they were ALT1, but otherwise everything works properly.
        stop(
            "Multiallelic variants and phase/dosage info simultaneously "
            "present; pvar required in this case");
    }
    _state_ptr =
        static_cast<plink2::PgenReader*>(malloc(sizeof(plink2::PgenReader)));
    if (!_state_ptr) {
        stop("Out of memory");
    }
    plink2::PreinitPgr(_state_ptr);
    plink2::PgrSetFreadBuf(nullptr, _state_ptr);
    const uintptr_t pgr_alloc_main_byte_ct =
        pgr_alloc_cacheline_ct * plink2::kCacheline;
    const uintptr_t sample_subset_byte_ct =
        plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) *
        plink2::kBytesPerVec;
    const uintptr_t cumulative_popcounts_byte_ct =
        plink2::DivUp(file_sample_ct,
                      plink2::kBitsPerWord * plink2::kInt32PerVec) *
        plink2::kBytesPerVec;
    genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) *
                      plink2::kBytesPerVec;

    const uintptr_t ac_byte_ct = plink2::RoundUpPow2(
        file_sample_ct * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
    const uintptr_t ac2_byte_ct = plink2::RoundUpPow2(
        file_sample_ct * 2 * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
    uintptr_t multiallelic_hc_byte_ct = 0;
    if (_info_ptr->max_allele_ct != 2) {
        multiallelic_hc_byte_ct =
            2 * sample_subset_byte_ct + ac_byte_ct + ac2_byte_ct;
    }
    const uintptr_t dosage_main_byte_ct =
        plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) *
        plink2::kBytesPerVec;
    unsigned char* pgr_alloc;
    if (plink2::cachealigned_malloc(
            pgr_alloc_main_byte_ct +
                (2 * plink2::kPglNypTransposeBatch + 5) *
                    sample_subset_byte_ct +
                cumulative_popcounts_byte_ct +
                (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct +
                multiallelic_hc_byte_ct + dosage_main_byte_ct +
                plink2::kPglBitTransposeBufbytes +
                4 * (plink2::kPglNypTransposeBatch *
                     plink2::kPglNypTransposeBatch / 8),
            &pgr_alloc)) {
        stop("Out of memory");
    }
    plink2::PglErr reterr =
        PgrInit(fname, max_vrec_width, _info_ptr, _state_ptr, pgr_alloc);
    if (reterr != plink2::kPglRetSuccess) {
        if (!plink2::PgrGetFreadBuf(_state_ptr)) {
            plink2::aligned_free(pgr_alloc);
        }
        sprintf(errstr_buf, "PgrInit() error %d", static_cast<int>(reterr));
        stop(errstr_buf);
    }
    unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
    _subset_include_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _subset_include_interleaved_vec =
        reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);

#ifdef USE_AVX2
    _subset_include_interleaved_vec[-3] = 0;
    _subset_include_interleaved_vec[-2] = 0;
#endif
    _subset_include_interleaved_vec[-1] = 0;

    _subset_cumulative_popcounts = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
    _pgv.genovec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
    if (multiallelic_hc_byte_ct) {
        _pgv.patch_01_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _pgv.patch_01_vals =
            reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[ac_byte_ct]);
        _pgv.patch_10_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _pgv.patch_10_vals =
            reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[ac2_byte_ct]);
    } else {
        _pgv.patch_01_set = nullptr;
        _pgv.patch_01_vals = nullptr;
        _pgv.patch_10_set = nullptr;
        _pgv.patch_10_vals = nullptr;
    }
    _pgv.phasepresent = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.phaseinfo = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.dosage_present = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.dosage_main = reinterpret_cast<uint16_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);
    if (sample_subset) {
        SetSampleSubsetInternal(sample_subset, subset_size);
    } else {
        _subset_size = file_sample_ct;
    }
    _transpose_batch_buf = reinterpret_cast<plink2::VecW*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);
    _multivar_vmaj_geno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter =
        &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * genovec_byte_ct]);
    _multivar_vmaj_phasepresent_buf =
        reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(
        pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
    _multivar_vmaj_phaseinfo_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(
        pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
    _multivar_smaj_geno_batch_buf =
        reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch *
                                      plink2::kPglNypTransposeBatch / 4]);
    _multivar_smaj_phaseinfo_batch_buf =
        reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch *
                                      plink2::kPglNypTransposeBatch / 8]);
    _multivar_smaj_phasepresent_batch_buf =
        reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    // pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch *
    // plink2::kPglNypTransposeBatch / 8]);
    if (_subset_size == 0) {
        stop("Subset must be nonnegative");
    }
}

void PlinkMatrix::Close() {
    // don't bother propagating file close errors for now
    if (_info_ptr) {
        CondReleaseRefcountedWptr(&_allele_idx_offsetsp);
        CondReleaseRefcountedWptr(&_nonref_flagsp);
        if (_info_ptr->vrtypes) {
            plink2::aligned_free(_info_ptr->vrtypes);
        }
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgfi(_info_ptr, &reterr);
        free(_info_ptr);
        _info_ptr = nullptr;
    }
    if (_state_ptr) {
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgr(_state_ptr, &reterr);
        if (PgrGetFreadBuf(_state_ptr)) {
            plink2::aligned_free(PgrGetFreadBuf(_state_ptr));
        }
        free(_state_ptr);
        _state_ptr = nullptr;
    }
    if (malloc_all) {
        for (uint32_t i = 0; i != _vsubset_size; ++i) {
            free(compactM[i]);
        }
        free(compactM);
    }
    _subset_size = 0;
}

void PlinkMatrix::SetSampleSubsetInternal(int* sample_subset_1based,
                                          const uint32_t subset_size) {
    const uint32_t raw_sample_ct = _info_ptr->raw_sample_ct;
    const uint32_t raw_sample_ctv =
        plink2::DivUp(raw_sample_ct, plink2::kBitsPerVec);
    const uint32_t raw_sample_ctaw = raw_sample_ctv * plink2::kWordsPerVec;
    uintptr_t* sample_include = _subset_include_vec;
    plink2::ZeroWArr(raw_sample_ctaw, sample_include);

    if (subset_size == 0) {
        stop("Empty sample_subset is not currently permitted");
    }
    uint32_t sample_uidx = sample_subset_1based[0] - 1;
    uint32_t idx = 0;
    uint32_t next_uidx;
    while (1) {
        if (sample_uidx >= raw_sample_ct) {
            char errstr_buf[256];
            sprintf(errstr_buf,
                    "sample number out of range (%d; must be 1..%u)",
                    static_cast<int>(sample_uidx + 1), raw_sample_ct);
            stop(errstr_buf);
        }
        plink2::SetBit(sample_uidx, sample_include);
        if (++idx == subset_size) {
            break;
        }
        next_uidx = sample_subset_1based[idx] - 1;

        // prohibit this since it implies that the caller expects genotypes to
        // be returned in a different order
        if (next_uidx <= sample_uidx) {
            stop("sample_subset is not in strictly increasing order");
        }
        sample_uidx = next_uidx;
    }
    plink2::FillInterleavedMaskVec(sample_include, raw_sample_ctv,
                                   _subset_include_interleaved_vec);
    const uint32_t raw_sample_ctl =
        plink2::DivUp(raw_sample_ct, plink2::kBitsPerWord);
    plink2::FillCumulativePopcounts(sample_include, raw_sample_ctl,
                                    _subset_cumulative_popcounts);
    plink2::PgrSetSampleSubsetIndex(_subset_cumulative_popcounts, _state_ptr,
                                    &_subset_index);
    _subset_size = subset_size;
}
PlinkMatrix::~PlinkMatrix() { Close(); }

void PlinkMatrix::ReadCompact(int* variant_subset,
                              const uintptr_t vsubset_size) {
    if (!_info_ptr) {
        stop("pgen is closed");
    }
    ni = (int)vsubset_size;
    _vsubset_size = vsubset_size;
    // assume that buf has the correct dimensions
    const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;

    compactM = (uintptr_t**)malloc(sizeof(uintptr_t*) * vsubset_size);
    if (!compactM) {
        stop("out of memory\n");
    }
    const uintptr_t byte_ct = genovec_byte_ct;  //(_subset_size + 3) / 4;

    for (uintptr_t col_idx = 0; col_idx != vsubset_size; ++col_idx) {
        compactM[col_idx] = (uintptr_t*)malloc(byte_ct);
        if (!compactM[col_idx]) {
            stop("out of memory\n");
        }

        uint32_t variant_idx = variant_subset[col_idx] - 1;
        if (static_cast<uint32_t>(variant_idx) >= raw_variant_ct) {
            char errstr_buf[256];
            sprintf(errstr_buf,
                    "variant_subset element out of range (%d; must be 1..%u)",
                    variant_idx + 1, raw_variant_ct);
            stop(errstr_buf);
        }
        plink2::PglErr reterr =
            PgrGet(_subset_include_vec, _subset_index, _subset_size,
                   variant_idx, _state_ptr, _pgv.genovec);
        memcpy(compactM[col_idx], _pgv.genovec, byte_ct);
        if (reterr != plink2::kPglRetSuccess) {
            char errstr_buf[256];
            sprintf(errstr_buf, "PgrGet() error %d", static_cast<int>(reterr));
            stop(errstr_buf);
        }
    }
    malloc_all = true;
}

double PlinkMatrix::dot_product(int j, const double* v) {
    if (!malloc_all) {
        stop("Must load the input matrix before calling this function\n");
    }
    if (j > _vsubset_size - 1) {
        stop("Column out of range\n");
    }

    return plink2::LinearCombinationMeanimpute(v, compactM[j], nullptr, nullptr,
                                               _subset_size, 0);
}

double PlinkMatrix::column_product(int i, int j) {
    if ((i > _vsubset_size - 1) || (j > _vsubset_size - 1)) {
        stop("Column out of range\n");
    }
    return plink2::genoarrproduct(compactM[i], compactM[j], _subset_size);
}

double PlinkMatrix::vx2(int j, const double* v) {
    if (!malloc_all) {
        stop("Must load the input matrix before calling this function\n");
    }
    if (j > _vsubset_size - 1) {
        stop("Column out of range\n");
    }
    return plink2::LinearCombinationSquare(v, compactM[j], _subset_size);
}

void PlinkMatrix::update_res(int j, double d, const double* v,
                             double* __restrict r) {
    plink2::update_res_raw(compactM[j], d, v, r, _subset_size);
}

// Wrapper to load the compact matrix
void PlinkMatrix::load_compact_matrix(const char* fname, int raw_sample_ct,
                                      int* sample_subset,
                                      const uint32_t subset_size,
                                      int* variant_subset,
                                      const uintptr_t vsubset_size) {
    Load(fname, raw_sample_ct, sample_subset, subset_size);
    ReadCompact(variant_subset, vsubset_size);
}