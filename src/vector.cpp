#include "amici/vector.h"

#include <algorithm>
#include <functional>

namespace amici {

AmiVector& AmiVector::operator=(AmiVector const& other) {
    vec_ = other.vec_;
    synchroniseNVector(other.get_ctx());
    return *this;
}

realtype* AmiVector::data() { return vec_.data(); }

realtype const* AmiVector::data() const { return vec_.data(); }

N_Vector AmiVector::getNVector() { return nvec_; }

const_N_Vector AmiVector::getNVector() const { return nvec_; }

std::vector<realtype> const& AmiVector::getVector() const { return vec_; }

int AmiVector::getLength() const { return gsl::narrow<int>(vec_.size()); }

void AmiVector::zero() { set(0.0); }

void AmiVector::minus() {
    std::ranges::transform(vec_, vec_.begin(), std::negate<realtype>());
}

void AmiVector::set(realtype const val) { std::ranges::fill(vec_, val); }

realtype& AmiVector::operator[](int const pos) {
    return vec_.at(gsl::narrow<decltype(vec_)::size_type>(pos));
}

realtype const& AmiVector::operator[](int const pos) const {
    return vec_.at(gsl::narrow<decltype(vec_)::size_type>(pos));
}

realtype& AmiVector::at(int const pos) {
    return vec_.at(gsl::narrow<decltype(vec_)::size_type>(pos));
}

realtype const& AmiVector::at(int const pos) const {
    return vec_.at(gsl::narrow<decltype(vec_)::size_type>(pos));
}

void AmiVector::copy(AmiVector const& other) {
    if (getLength() != other.getLength())
        throw AmiException(
            "Dimension of AmiVector (%i) does not "
            "match input dimension (%i)",
            getLength(), other.getLength()
        );
    std::ranges::copy(other.vec_, vec_.begin());
}

void AmiVector::synchroniseNVector(SUNContext const sunctx) {
    if (nvec_)
        N_VDestroy_Serial(nvec_);
    if (sunctx) {
        nvec_ = N_VMake_Serial(
            gsl::narrow<long int>(vec_.size()), vec_.data(), sunctx
        );
    }
}

AmiVector::~AmiVector() {
    if (nvec_)
        N_VDestroy_Serial(nvec_);
}

AmiVectorArray::AmiVectorArray(
    long int const length_inner, long int const length_outer, SUNContext const sunctx
)
    : vec_array_(length_outer, AmiVector(length_inner, sunctx)) {
    nvec_array_.resize(length_outer);
    for (int idx = 0; idx < length_outer; idx++) {
        nvec_array_.at(idx) = vec_array_.at(idx).getNVector();
    }
}

AmiVectorArray& AmiVectorArray::operator=(AmiVectorArray const& other) {
    vec_array_ = other.vec_array_;
    nvec_array_.resize(other.getLength());
    for (int idx = 0; idx < other.getLength(); idx++) {
        nvec_array_.at(idx) = vec_array_.at(idx).getNVector();
    }
    return *this;
}

AmiVectorArray::AmiVectorArray(AmiVectorArray const& vaold)
    : vec_array_(vaold.vec_array_) {
    nvec_array_.resize(vaold.getLength());
    for (int idx = 0; idx < vaold.getLength(); idx++) {
        nvec_array_.at(idx) = vec_array_.at(idx).getNVector();
    }
}

realtype* AmiVectorArray::data(int pos) { return vec_array_.at(pos).data(); }

realtype const* AmiVectorArray::data(int const pos) const {
    return vec_array_.at(pos).data();
}

realtype& AmiVectorArray::at(int const ipos, int const jpos) {
    return vec_array_.at(jpos).at(ipos);
}

realtype const& AmiVectorArray::at(int const ipos, int const jpos) const {
    return vec_array_.at(jpos).at(ipos);
}

N_Vector* AmiVectorArray::getNVectorArray() { return nvec_array_.data(); }

N_Vector AmiVectorArray::getNVector(int const pos) { return nvec_array_.at(pos); }

const_N_Vector AmiVectorArray::getNVector(int const pos) const {
    return nvec_array_.at(pos);
}

AmiVector& AmiVectorArray::operator[](int const pos) { return vec_array_.at(pos); }

AmiVector const& AmiVectorArray::operator[](int const pos) const {
    return vec_array_.at(pos);
}

int AmiVectorArray::getLength() const {
    return gsl::narrow<int>(vec_array_.size());
}

void AmiVectorArray::zero() {
    for (auto& v : vec_array_)
        v.zero();
}

void AmiVectorArray::flatten_to_vector(std::vector<realtype>& vec) const {
    int n_outer = gsl::narrow<int>(vec_array_.size());
    if (n_outer == 0)
        return; // nothing to do ...
    int n_inner = vec_array_.at(0).getLength();

    if (gsl::narrow<int>(vec.size()) != n_inner * n_outer) {
        throw AmiException(
            "Dimension of AmiVectorArray (%ix%i) does not "
            "match target vector dimension (%u)",
            n_inner, n_outer, vec.size()
        );
    }

    for (int outer = 0; outer < n_outer; ++outer) {
        for (int inner = 0; inner < n_inner; ++inner)
            vec.at(inner + outer * n_inner) = this->at(inner, outer);
    }
}

void AmiVectorArray::copy(AmiVectorArray const& other) {
    if (getLength() != other.getLength())
        throw AmiException(
            "Dimension of AmiVectorArray (%i) does not "
            "match input dimension (%i)",
            getLength(), other.getLength()
        );

    for (int iv = 0; iv < getLength(); ++iv) {
        vec_array_.at(iv).copy(other.vec_array_.at(iv));
        nvec_array_[iv] = vec_array_.at(iv).getNVector();
    }
}

} // namespace amici
