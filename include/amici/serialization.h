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
    ar &u.atolB;
    ar &u.rtolB;
    ar &u.atol_fsa;
    ar &u.rtol_fsa;
    ar &u.quad_atol;
    ar &u.quad_rtol;
    ar &u.ss_atol;
    ar &u.ss_rtol;
    ar &u.ss_atol_sensi;
    ar &u.ss_rtol_sensi;
    ar &u.maxsteps;
    ar &u.maxstepsB;
    ar &u.requires_preequilibration;
    ar &u.newton_maxsteps;
    ar &u.newton_maxlinsteps;
    ar &u.newton_damping_factor_mode;
    ar &u.newton_damping_factor_lower_bound;
    ar &u.ism;
    ar &u.sensi_meth;
    ar &u.linsol;
    ar &u.interpType;
    ar &u.lmm;
    ar &u.iter;
    ar &u.stldet;
    ar &u.ordering;
    ar &u.cpu_time;
    ar &u.cpu_timeB;
    ar &u.rdata_mode;
}


template <class Archive>
void serialize(Archive &ar, amici::CVodeSolver &u, const unsigned int version) {
    ar & static_cast<amici::Solver&>(u);
}

template <class Archive>
void serialize(Archive &ar, amici::Model &u, const unsigned int version) {
    ar &u.nx_rdata;
    ar &u.nxtrue_rdata;
    ar &u.nx_solver;
    ar &u.nxtrue_solver;
    ar &u.ny;
    ar &u.nytrue;
    ar &u.nz;
    ar &u.nztrue;
    ar &u.ne;
    ar &u.nw;
    ar &u.ndwdx;
    ar &u.ndwdp;
    ar &u.ndxdotdw;
    ar &u.nnz;
    ar &u.nJ;
    ar &u.ubw;
    ar &u.lbw;
    ar &u.o2mode;
    ar &u.z2event;
    ar &u.idlist;
    ar &u.state.h;
    ar &u.state.unscaledParameters;
    ar &u.originalParameters;
    ar &u.state.fixedParameters;
    ar &u.reinitializeFixedParameterInitialStates;
    ar &u.state.plist;
    ar &u.x0data;
    ar &u.sx0data;
    ar &u.ts;
    ar &u.nmaxevent;
    ar &u.pscale;
    ar &u.tstart;
    ar &u.stateIsNonNegative;
    ar &u.pythonGenerated;
    ar &u.ndxdotdp_explicit;
    ar &u.ndxdotdp_implicit;
}


template <class Archive>
void serialize(Archive &ar, amici::ReturnData &r, const unsigned int version) {
    ar &r.np;
    ar &r.nk;
    ar &r.nx;
    ar &r.nx_solver;
    ar &r.nxtrue;
    ar &r.ny;
    ar &r.nytrue;
    ar &r.nz;
    ar &r.nztrue;
    ar &r.ne;
    ar &r.nJ;
    ar &r.nplist;
    ar &r.nmaxevent;
    ar &r.nt;
    ar &r.newton_maxsteps;
    ar &r.pscale;
    ar &r.o2mode;
    ar &r.sensi;
    ar &r.sensi_meth;

    ar &r.ts;
    ar &r.xdot;
    ar &r.J;
    ar &r.w;
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
    ar &r.cpu_time;
    ar &r.cpu_timeB;
    ar &r.preeq_cpu_time;
    ar &r.preeq_status;
    ar &r.preeq_numsteps;
    ar &r.preeq_numlinsteps;
    ar &r.preeq_wrms;
    ar &r.preeq_t;
    ar &r.posteq_cpu_time;
    ar &r.posteq_status;
    ar &r.posteq_numsteps;
    ar &r.posteq_numlinsteps;
    ar &r.posteq_wrms;
    ar &r.posteq_t;
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
     *
     * @return The object serialized as char
     */
    try {
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
    } catch(boost::archive::archive_exception const& e) {
        throw AmiException("Serialization to char failed: %s", e.what());
    }
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
    try {
        ::boost::iostreams::basic_array_source<char> device(buffer, size);
        ::boost::iostreams::stream<::boost::iostreams::basic_array_source<char>> s(
            device);
        ::boost::archive::binary_iarchive iar(s);
        T data;
        iar >> data;

        return data;
    } catch(::boost::archive::archive_exception const& e) {
        throw AmiException("Deserialization from char failed: %s", e.what());
    }
}


template <typename T>
std::string serializeToString(T const& data) {
    /**
     * @brief Serialize object to string
     *
     * @param data input object
     *
     * @return The object serialized as string
     */
    try {
        std::string serialized;
        ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
        ::boost::iostreams::stream<
            ::boost::iostreams::back_insert_device<std::string>>
            s(inserter);
        ::boost::archive::binary_oarchive oar(s);

        oar << data;
        s.flush();

        return serialized;
    } catch(::boost::archive::archive_exception const& e) {
        throw AmiException("Serialization to string failed: %s", e.what());
    }
}

template <typename T>
std::vector<char> serializeToStdVec(T const& data) {
    /**
     * @brief Serialize object to std::vector<char>
     *
     * @param data input object
     *
     * @return The object serialized as std::vector<char>
     */
    try{
        std::string serialized;
        ::boost::iostreams::back_insert_device<std::string> inserter(serialized);
        ::boost::iostreams::stream<::boost::iostreams::back_insert_device<std::string>>
            s(inserter);
        ::boost::archive::binary_oarchive oar(s);

        oar << data;
        s.flush();

        std::vector<char> buf(serialized.begin(), serialized.end());

        return buf;
    } catch(::boost::archive::archive_exception const& e) {
        throw AmiException("Serialization to StdVec failed: %s", e.what());
    }
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
    try{
        ::boost::iostreams::basic_array_source<char> device(serialized.data(),
                                                          serialized.size());
        ::boost::iostreams::stream<::boost::iostreams::basic_array_source<char>> s(
            device);
        ::boost::archive::binary_iarchive iar(s);
        T deserialized;

        iar >> deserialized;

        return deserialized;
    } catch(::boost::archive::archive_exception const& e) {
        throw AmiException("Deserialization from StdVec failed: %s", e.what());
    }
}


} // namespace amici
#endif // AMICI_SERIALIZATION_H
