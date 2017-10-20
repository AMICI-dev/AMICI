#ifndef AMICI_SERIALIZATION_H
#define AMICI_SERIALIZATION_H

#include "include/udata.h"
#include "include/rdata.h"
#include <cassert>
#include <boost/serialization/array.hpp>
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
void archiveRawArray(Archive &ar, T **p, int size) {
    if (Archive::is_loading::value) {
        assert(*p == nullptr); // ensure it's unset, otherwise would have to deallocate
        ar &size;
        *p = size ? new T[size] : nullptr;
    } else {
        size = *p == nullptr ? 0 : size;
        ar &size;
    }
    ar &make_array<T>(*p, size);
}

template <class Archive>
void serialize(Archive &ar, amici::UserData &u, const unsigned int version) {
    ar &const_cast<int &>(u.np);
    ar &const_cast<int &>(u.nk);
    ar &const_cast<int &>(u.nx);
    ar &u.pscale;
    ar &u.nmaxevent;
    ar &u.nplist;
    ar &u.nt;
    ar &u.tstart;
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

    archiveRawArray(ar, &u.qpositivex, u.nx);
    archiveRawArray(ar, &u.plist, u.nplist);
    archiveRawArray(ar, &u.p, u.np);
    archiveRawArray(ar, &u.k, u.nk);
    archiveRawArray(ar, &u.ts, u.nt);
    archiveRawArray(ar, &u.pbar, u.np);
    archiveRawArray(ar, &u.xbar, u.nx);
    archiveRawArray(ar, &u.x0data, u.nx);
    archiveRawArray(ar, &u.sx0data, u.np * u.nx);
}


template <class Archive>
void serialize(Archive &ar, amici::ReturnData &r, const unsigned int version) {
    ar &r.freeFieldsOnDestruction;

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

    archiveRawArray(ar, &r.ts, r.nt);
    archiveRawArray(ar, &r.xdot, r.nx);
    archiveRawArray(ar, &r.J, r.nx * r.nx);
    archiveRawArray(ar, &r.z, r.nmaxevent * r.nz);
    archiveRawArray(ar, &r.sigmaz, r.nmaxevent * r.nz);
    archiveRawArray(ar, &r.sz, r.nmaxevent * r.nz * r.nplist);
    archiveRawArray(ar, &r.ssigmaz, r.nmaxevent * r.nz * r.nplist);
    archiveRawArray(ar, &r.rz, r.nmaxevent * r.nz);
    archiveRawArray(ar, &r.srz, r.nmaxevent * r.nz * r.nplist);
    archiveRawArray(ar, &r.s2rz, r.nmaxevent * r.nz * r.nplist * r.nplist);
    archiveRawArray(ar, &r.x, r.nt * r.nx);
    archiveRawArray(ar, &r.sx, r.nt * r.nx * r.nplist);
    archiveRawArray(ar, &r.y, r.nt * r.ny);
    archiveRawArray(ar, &r.sigmay, r.nt * r.ny);
    archiveRawArray(ar, &r.sy, r.nt * r.ny * r.nplist);
    archiveRawArray(ar, &r.ssigmay, r.nt * r.ny * r.nplist);

    archiveRawArray(ar, &r.numsteps, r.nt);
    archiveRawArray(ar, &r.numstepsB, r.nt);
    archiveRawArray(ar, &r.numrhsevals, r.nt);
    archiveRawArray(ar, &r.numrhsevalsB, r.nt);
    archiveRawArray(ar, &r.numerrtestfails, r.nt);
    archiveRawArray(ar, &r.numerrtestfailsB, r.nt);
    archiveRawArray(ar, &r.numnonlinsolvconvfails, r.nt);
    archiveRawArray(ar, &r.numnonlinsolvconvfailsB, r.nt);
    archiveRawArray(ar, &r.order, r.nt);

    archiveRawArray(ar, &r.newton_status, r.nt);
    archiveRawArray(ar, &r.newton_time, r.nt);
    archiveRawArray(ar, &r.newton_numsteps, r.nt);
    archiveRawArray(ar, &r.newton_numlinsteps, r.nt);
    archiveRawArray(ar, &r.x0, r.nx);
    archiveRawArray(ar, &r.sx0, r.nx * r.nplist);

    archiveRawArray(ar, &r.llh, 1);
    archiveRawArray(ar, &r.chi2, 1);
    archiveRawArray(ar, &r.sllh, r.nplist);
    archiveRawArray(ar, &r.s2llh, r.nplist * r.nplist);
    archiveRawArray(ar, &r.status, 1);
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
 * @brief Deserialize AMICI::UserData that has been serialized using
 * serializeAmiciUserData
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
