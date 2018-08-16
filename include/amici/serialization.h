#ifndef AMICI_SERIALIZATION_H
#define AMICI_SERIALIZATION_H

#include "amici/rdata.h"
#include "amici/model.h"
#include "amici/solver.h"
#include "amici/solver_cvodes.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>

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
void serialize(Archive &ar, amici::CVodeSolver &u, const unsigned int version) {
    ar & static_cast<amici::Solver&>(u);
}

template <class Archive>
void serialize(Archive &ar, amici::Model &u, const unsigned int version) {
    ar &const_cast<int &>(u.nx);
    ar &const_cast<int &>(u.nxtrue);
    ar &const_cast<int &>(u.ny);
    ar &const_cast<int &>(u.nytrue);
    ar &const_cast<int &>(u.nz);
    ar &const_cast<int &>(u.nztrue);
    ar &const_cast<int &>(u.ne);
    ar &const_cast<int &>(u.nw);
    ar &const_cast<int &>(u.ndwdx);
    ar &const_cast<int &>(u.ndwdp);
    ar &const_cast<int &>(u.nnz);
    ar &const_cast<int &>(u.nJ);
    ar &const_cast<int &>(u.ubw);
    ar &const_cast<int &>(u.lbw);
    ar &const_cast<amici::SecondOrderMode &>(u.o2mode);
    ar &const_cast<std::vector<int> &>(u.z2event);
    ar &const_cast<std::vector<realtype> &>(u.idlist);
    ar &u.h;
    ar &u.unscaledParameters;
    ar &u.originalParameters;
    ar &u.fixedParameters;
    ar &u.plist_;
    ar &u.x0data;
    ar &u.sx0data;
    ar &u.ts;
    ar &u.nmaxevent;
    ar &u.pscale;
    ar &u.tstart;

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
    ar &r.pscale;
    ar &const_cast<amici::SecondOrderMode &>(r.o2mode);
    ar &const_cast<amici::SensitivityOrder &>(r.sensi);
    ar &const_cast<amici::SensitivityMethod &>(r.sensi_meth);

    ar &const_cast<std::vector<realtype> &>(r.ts);
    ar &r.xdot;
    ar &r.J;
    ar &r.z & r.sigmaz;
    ar &r.sz &r.ssigmaz;
    ar &r.rz;
    ar &r.srz;
    ar &r.s2rz;
    ar &r.x;
    ar &r.sx;
    ar &r.y & r.sigmay;
    ar &r.sy & r.ssigmay;

    ar &r.numsteps;
    ar &r.numstepsB;
    ar &r.numrhsevals;
    ar &r.numrhsevalsB;
    ar &r.numerrtestfails;
    ar &r.numerrtestfailsB;
    ar &r.numnonlinsolvconvfails;
    ar &r.numnonlinsolvconvfailsB;
    ar &r.order;

    ar &r.newton_status;
    ar &r.newton_time;
    ar &r.newton_numsteps;
    ar &r.newton_numlinsteps;
    ar &r.x0;
    ar &r.sx0;
    ar &r.llh;
    ar &r.chi2;
    ar &r.sllh;
    ar &r.s2llh;
    ar &r.status;
}


} // namespace serialization
} // namespace boost

namespace amici {

template <typename T>
char *serializeToChar(T const& data, int *size) {
    /**
     * @brief Serialize object to char array
     *
     * @param data input object
     * @param size maximum char length

     * @return The object serialized as char
     */
    std::string serialized;
    ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
    ::boost::iostreams::stream<::boost::iostreams::back_insert_device<std::string>>
        s(inserter);
    ::boost::archive::binary_oarchive oar(s);
    oar << data;
    s.flush();

    char *charBuffer = new char[serialized.size()];
    memcpy(charBuffer, serialized.data(), serialized.size());

    if (size)
        *size = serialized.size();

    return charBuffer;
}



template <typename T>
T deserializeFromChar(const char *buffer, int size) {
    /**
     * @brief Deserialize object that has been serialized using serializeToChar
     *
     * @param buffer serialized object
     * @param size length of buffer
     *
     * @return The deserialized object
     */
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
    /**
     * @brief Serialize object to string
     *
     * @param data input object

     * @return The object serialized as string
     */
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
std::vector<char> serializeToStdVec(T const& data) {
    /**
     * @brief Serialize object to std::vector<char>
     *
     * @param data input object

     * @return The object serialized as std::vector<char>
     */
    std::string serialized;
    ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
    ::boost::iostreams::stream<::boost::iostreams::back_insert_device<std::string>>
        s(inserter);
    ::boost::archive::binary_oarchive oar(s);
    oar << data;
    s.flush();

    std::vector<char> buf(serialized.begin(), serialized.end());

    return buf;
}

template <typename T>
T deserializeFromString(std::string const& serialized) {
    /**
     * @brief Deserialize object that has been serialized using serializeToString
     *
     * @param serialized serialized object
     *
     * @return The deserialized object
     */
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
