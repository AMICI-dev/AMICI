#include "amici/vector.h"

#include <functional>

namespace amici {

AmiVector &AmiVector::operator=(AmiVector const &other) {
    vec = other.vec;
    synchroniseNVector();
    return *this;
}

realtype *AmiVector::data() { return vec.data(); }

const realtype *AmiVector::data() const { return vec.data(); }

N_Vector AmiVector::getNVector() { return nvec; }

const_N_Vector AmiVector::getNVector() const { return nvec; }

std::vector<realtype> const &AmiVector::getVector() const { return vec; }

void AmiVector::clear() { vec.clear(); N_VDestroy_Serial(nvec); }

int AmiVector::getLength() const { return static_cast<int>(vec.size()); }

void AmiVector::reset() { set(0.0); }

void AmiVector::minus() {
    std::transform(vec.begin(), vec.end(),
                   vec.begin(), std::negate<realtype>());
}

void AmiVector::set(realtype val) { std::fill(vec.begin(), vec.end(), val); }

realtype &AmiVector::operator[](int pos) {
    return vec.at(static_cast<decltype(vec)::size_type>(pos));
}

realtype &AmiVector::at(int pos) {
    return vec.at(static_cast<decltype(vec)::size_type>(pos));
}

const realtype &AmiVector::at(int pos) const {
    return vec.at(static_cast<decltype(vec)::size_type>(pos));
}

void AmiVector::copy(const AmiVector &other) {
    if(getLength() != other.getLength())
        throw AmiException("Dimension of AmiVector (%i) does not "
                           "match input dimension (%i)",
                           getLength(), other.getLength());
    std::copy(other.vec.begin(), other.vec.end(), vec.begin());
    synchroniseNVector();
}

void AmiVector::synchroniseNVector() {
    if (nvec)
        N_VDestroy_Serial(nvec);
    nvec = N_VMake_Serial(static_cast<long int>(vec.size()), vec.data());
}

AmiVector::~AmiVector() {
    if (nvec)
        N_VDestroy_Serial(nvec);
}

AmiVectorArray::AmiVectorArray(long int length_inner, long int length_outer)
    : vec_array(length_outer, AmiVector(length_inner)) {
    nvec_array.resize(length_outer);
    for (int idx = 0; idx < length_outer; idx++) {
        nvec_array.at(idx) = vec_array.at(idx).getNVector();
    }
}

AmiVectorArray &AmiVectorArray::operator=(AmiVectorArray const &other) {
    vec_array = other.vec_array;
    nvec_array.resize(other.getLength());
    for (int idx = 0; idx < other.getLength(); idx++) {
        nvec_array.at(idx) = vec_array.at(idx).getNVector();
    }
    return *this;
}

AmiVectorArray::AmiVectorArray(const AmiVectorArray &vaold)
    : vec_array(vaold.vec_array) {
    nvec_array.resize(vaold.getLength());
    for (int idx = 0; idx < vaold.getLength(); idx++) {
        nvec_array.at(idx) = vec_array.at(idx).getNVector();
    }
}

realtype *AmiVectorArray::data(int pos) { return vec_array.at(pos).data(); }

const realtype *AmiVectorArray::data(int pos) const {
    return vec_array.at(pos).data();
}

realtype &AmiVectorArray::at(int ipos, int jpos) {
    return vec_array.at(jpos).at(ipos);
}

const realtype &AmiVectorArray::at(int ipos, int jpos) const {
    return vec_array.at(jpos).at(ipos);
}

N_Vector *AmiVectorArray::getNVectorArray() { return nvec_array.data(); }

N_Vector AmiVectorArray::getNVector(int pos) { return nvec_array.at(pos); }

AmiVector &AmiVectorArray::operator[](int pos) { return vec_array.at(pos); }

const AmiVector &AmiVectorArray::operator[](int pos) const {
    return vec_array.at(pos);
}

int AmiVectorArray::getLength() const {
    return static_cast<int>(vec_array.size());
}

void AmiVectorArray::reset() {
    for (auto &v : vec_array)
        v.reset();
}

void AmiVectorArray::flatten_to_vector(std::vector<realtype> &vec) const {
    int n_outer = static_cast<int>(vec_array.size());
    if (n_outer == 0)
        return; // nothing to do ...
    int n_inner = vec_array.at(0).getLength();

    if (static_cast<int>(vec.size()) != n_inner * n_outer) {
        throw AmiException("Dimension of AmiVectorArray (%ix%i) does not "
                           "match target vector dimension (%u)",
                           n_inner, n_outer, vec.size());
    }

    for (int outer = 0; outer < n_outer; ++outer) {
        for (int inner = 0; inner < n_inner; ++inner)
            vec.at(inner + outer * n_inner) = this->at(inner, outer);
    }
}

void AmiVectorArray::copy(const AmiVectorArray &other) {
    if (getLength() != other.getLength())
        throw AmiException("Dimension of AmiVectorArray (%i) does not "
                           "match input dimension (%i)",
                           getLength(), other.getLength());

    for (int iv = 0; iv < getLength(); ++iv) {
        vec_array.at(iv).copy(other.vec_array.at(iv));
        nvec_array[iv] = vec_array.at(iv).getNVector();
    }
}

} // namespace amici
