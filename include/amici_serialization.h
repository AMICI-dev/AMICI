#ifndef AMICI_SERIALIZATION_H
#define AMICI_SERIALIZATION_H

#include "include/rdata.h"
#include "include/amici_model.h"
#include "include/amici_solver.h"
#include <cassert>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

/* Helper functions and forward declarations for boost::serialization */
namespace boost {
namespace serialization {

template <class Archive, typename T>
void archiveVector(Archive &ar, T **p, int size) {
    if (Archive::is_loading::value) {
        if(*p != nullptr)
            delete[] *p;
        ar &size;
        *p = size ? new T[size] : nullptr;
    } else {
        size = *p == nullptr ? 0 : size;
        ar &size;
    }
    ar &make_array<T>(*p, size);
}

template <class Archive>
void serialize(Archive &ar, amici::Solver &u, const unsigned int version) {
    ar &u.sensi;
    ar &u.atol;
    ar &u.rtol;
    ar &u.maxsteps;
    ar &u.ism;
    ar &u.sensi_meth;
    ar &u.linsol;
    ar &u.interpType;
    ar &u.lmm;
    ar &u.iter;
    ar &u.stldet;
    ar &u.ordering;
}

template <class Archive>
void serialize(Archive &ar, amici::Model &u, const unsigned int version) {
    ar &u.pscale;
    ar &u.nmaxevent;
    ar &u.tstart;
    ar &u.qpositivex;
    ar &u.plist_;
    ar &u.unscaledParameters;
    ar &u.originalParameters;
    ar &u.k_;
    ar &u.ts;
    ar &u.pbar;
    ar &u.x0data;
    ar &u.sx0data;
}


template <class Archive>
void serialize(Archive &ar, amici::ReturnData &r, const unsigned int version) {

    ar &const_cast<int &>(r.np);
    ar &const_cast<int &>(r.nk);
    ar &const_cast<int &>(r.nx);
    ar &const_cast<int &>(r.nxtrue);
    ar &const_cast<int &>(r.ny);
    ar &const_cast<int &>(r.nytrue);
    ar &const_cast<int &>(r.nz);
    ar &const_cast<int &>(r.nztrue);
    ar &const_cast<int &>(r.ne);
    ar &const_cast<int &>(r.nJ);
    ar &const_cast<int &>(r.nplist);
    ar &const_cast<int &>(r.nmaxevent);
    ar &const_cast<int &>(r.nt);
    ar &const_cast<int &>(r.newton_maxsteps);
    ar &const_cast<amici::AMICI_parameter_scaling &>(r.pscale);
    ar &const_cast<amici::AMICI_o2mode &>(r.o2mode);
    ar &const_cast<amici::AMICI_sensi_order &>(r.sensi);
    ar &const_cast<amici::AMICI_sensi_meth &>(r.sensi_meth);

    archiveVector(ar, &r.ts, r.nt);
    archiveVector(ar, &r.xdot, r.nx);
    archiveVector(ar, &r.J, r.nx * r.nx);
    archiveVector(ar, &r.z, r.nmaxevent * r.nz);
    archiveVector(ar, &r.sigmaz, r.nmaxevent * r.nz);
    archiveVector(ar, &r.sz, r.nmaxevent * r.nz * r.nplist);
    archiveVector(ar, &r.ssigmaz, r.nmaxevent * r.nz * r.nplist);
    archiveVector(ar, &r.rz, r.nmaxevent * r.nz);
    archiveVector(ar, &r.srz, r.nmaxevent * r.nz * r.nplist);
    archiveVector(ar, &r.s2rz, r.nmaxevent * r.nz * r.nplist * r.nplist);
    archiveVector(ar, &r.x, r.nt * r.nx);
    archiveVector(ar, &r.sx, r.nt * r.nx * r.nplist);
    archiveVector(ar, &r.y, r.nt * r.ny);
    archiveVector(ar, &r.sigmay, r.nt * r.ny);
    archiveVector(ar, &r.sy, r.nt * r.ny * r.nplist);
    archiveVector(ar, &r.ssigmay, r.nt * r.ny * r.nplist);

    archiveVector(ar, &r.numsteps, r.nt);
    archiveVector(ar, &r.numstepsB, r.nt);
    archiveVector(ar, &r.numrhsevals, r.nt);
    archiveVector(ar, &r.numrhsevalsB, r.nt);
    archiveVector(ar, &r.numerrtestfails, r.nt);
    archiveVector(ar, &r.numerrtestfailsB, r.nt);
    archiveVector(ar, &r.numnonlinsolvconvfails, r.nt);
    archiveVector(ar, &r.numnonlinsolvconvfailsB, r.nt);
    archiveVector(ar, &r.order, r.nt);

    archiveVector(ar, &r.newton_status, r.nt);
    archiveVector(ar, &r.newton_time, r.nt);
    archiveVector(ar, &r.newton_numsteps, r.nt);
    archiveVector(ar, &r.newton_numlinsteps, r.nt);
    archiveVector(ar, &r.x0, r.nx);
    archiveVector(ar, &r.sx0, r.nx * r.nplist);

    archiveVector(ar, &r.llh, 1);
    archiveVector(ar, &r.chi2, 1);
    archiveVector(ar, &r.sllh, r.nplist);
    archiveVector(ar, &r.s2llh, r.nplist * r.nplist);
    archiveVector(ar, &r.status, 1);
}


} // namespace serialization
} // namespace boost

namespace amici {

template <typename T>
char *serializeToChar(const T *data, int *size) {
    std::string serialized;
    ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
    ::boost::iostreams::stream<::boost::iostreams::back_insert_device<std::string>>
        s(inserter);
    ::boost::archive::binary_oarchive oar(s);
    oar << *data;
    s.flush();

    char *charBuffer = new char[serialized.size()];
    memcpy(charBuffer, serialized.data(), serialized.size());

    if (size)
        *size = serialized.size();

    return charBuffer;
}

/**
 * @brief Deserialize object that has been serialized using serializeToChar
 * @param buffer
 * @param size
 * @return The deserialized object
 */

template <typename T>
T deserializeFromChar(const char *buffer, int size) {
    ::boost::iostreams::basic_array_source<char> device(buffer, size);
    ::boost::iostreams::stream<::boost::iostreams::basic_array_source<char>> s(
        device);
    ::boost::archive::binary_iarchive iar(s);
    T data;
    iar >> data;

    return data;
}


template <typename T>
std::string serializeToString(T const& data) {
    std::string serialized;
    ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
    ::boost::iostreams::stream<
        ::boost::iostreams::back_insert_device<std::string>>
        s(inserter);
    ::boost::archive::binary_oarchive oar(s);
    oar << data;
    s.flush();

    return serialized;
}

template <typename T>
std::vector<char> serializeToStdVec(const T *data) {
    std::string serialized;
    ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
    ::boost::iostreams::stream<::boost::iostreams::back_insert_device<std::string>>
        s(inserter);
    ::boost::archive::binary_oarchive oar(s);
    oar << *data;
    s.flush();

    std::vector<char> buf(serialized.begin(), serialized.end());

    return buf;
}

template <typename T>
T deserializeFromString(std::string const& serialized) {
    ::boost::iostreams::basic_array_source<char> device(serialized.data(),
                                                      serialized.size());
    ::boost::iostreams::stream<::boost::iostreams::basic_array_source<char>> s(
        device);
    ::boost::archive::binary_iarchive iar(s);
    T deserialized;
    iar >> deserialized;

    return deserialized;
}


} // namespace amici
#endif // AMICI_SERIALIZATION_H
