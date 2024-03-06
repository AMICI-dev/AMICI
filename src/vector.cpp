#include "amici/vector.h"

#include <algorithm>
#include <functional>

namespace amici {

AmiVector& AmiVector::operator=(AmiVector const& other) {
    vec_ = other.vec_;
    synchroniseNVector();
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
    std::transform(
        vec_.begin(), vec_.end(), vec_.begin(), std::negate<realtype>()
    );
}

void AmiVector::set(realtype val) { std::fill(vec_.begin(), vec_.end(), val); }

realtype& AmiVector::operator[](int pos) {
    return vec_.at(gsl::narrow<decltype(vec_)::size_type>(pos));
}

realtype& AmiVector::at(int pos) {
    return vec_.at(gsl::narrow<decltype(vec_)::size_type>(pos));
}

realtype const& AmiVector::at(int pos) const {
    return vec_.at(gsl::narrow<decltype(vec_)::size_type>(pos));
}

void AmiVector::copy(AmiVector const& other) {
    if (getLength() != other.getLength())
        throw AmiException(
            "Dimension of AmiVector (%i) does not "
            "match input dimension (%i)",
            getLength(), other.getLength()
        );
    std::copy(other.vec_.begin(), other.vec_.end(), vec_.begin());
}

void AmiVector::synchroniseNVector() {
    if (nvec_)
        N_VDestroy_Serial(nvec_);
    nvec_ = N_VMake_Serial(gsl::narrow<long int>(vec_.size()), vec_.data());
}

AmiVector::~AmiVector() {
    if (nvec_)
        N_VDestroy_Serial(nvec_);
}

AmiVectorArray::AmiVectorArray(long int length_inner, long int length_outer)
    : vec_array_(length_outer, AmiVector(length_inner)) {
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

realtype const* AmiVectorArray::data(int pos) const {
    return vec_array_.at(pos).data();
}

realtype& AmiVectorArray::at(int ipos, int jpos) {
    return vec_array_.at(jpos).at(ipos);
}

realtype const& AmiVectorArray::at(int ipos, int jpos) const {
    return vec_array_.at(jpos).at(ipos);
}

N_Vector* AmiVectorArray::getNVectorArray() { return nvec_array_.data(); }

N_Vector AmiVectorArray::getNVector(int pos) { return nvec_array_.at(pos); }

const_N_Vector AmiVectorArray::getNVector(int pos) const {
    return nvec_array_.at(pos);
}

AmiVector& AmiVectorArray::operator[](int pos) { return vec_array_.at(pos); }

AmiVector const& AmiVectorArray::operator[](int pos) const {
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
