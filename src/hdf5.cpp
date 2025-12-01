/**
 * Functions for HDF5 I/O
 *
 * NOTE: Use only `const char*` versions of any HDF5 functions, not the
 * `std::string` version. On many systems, HDF5 libraries are still not compiled
 * with C++11 support, but use the old C++03 ABI, which will lead to linking
 * issues.
 */

#include <amici/hdf5.h>

#include <amici/edata.h>
#include <amici/logging.h>
#include <amici/model.h>
#include <amici/rdata.h>
#include <amici/solver.h>

#include <hdf5_hl.h>

namespace amici::hdf5 {

/**
 * @brief assertMeasurementDimensionsCompatible
 * @param m
 * @param n
 * @param model
 */
void check_measurement_dimensions_compatible(
    hsize_t const m, hsize_t const n, Model const& model
) {
    bool compatible = true;
    // if this is rank 1, n and m can be swapped
    if (n == 1) {
        compatible
            &= (n == (unsigned)model.nt() || n == (unsigned)model.nytrue);
        compatible
            &= (m == (unsigned)model.nytrue || m == (unsigned)model.nt());
        compatible &= (m * n == (unsigned)model.nytrue * model.nt());
    } else {
        compatible &= (n == (unsigned)model.nytrue);
        compatible &= (m == (unsigned)model.nt());
    }

    if (!compatible)
        throw(AmiException(
            "HDF5 measurement data does not match model. "
            "Incompatible dimensions."
        ));
}

/**
 * @brief assertEventDimensionsCompatible
 * @param m
 * @param n
 * @param model
 */
void check_event_dimensions_compatible(
    hsize_t const m, hsize_t const n, Model const& model
) {
    bool compatible = true;

    // if this is rank 1, n and m can be swapped
    if (n == 1) {
        compatible
            &= (n == (unsigned)model.n_max_event()
                || n == (unsigned)model.nztrue);
        compatible
            &= (m == (unsigned)model.nztrue
                || m == (unsigned)model.n_max_event());
        compatible &= (m * n == (unsigned)model.nytrue * model.n_max_event());
    } else {
        compatible &= (n == (unsigned)model.nztrue);
        compatible &= (m == (unsigned)model.n_max_event());
    }

    if (!compatible)
        throw(AmiException(
            "HDF5 event data does not match model. "
            "Incompatible dimensions."
        ));
}

void create_group(
    H5::H5File const& file, std::string const& groupPath, bool recursively
) {
#if H5_VERSION_GE(1, 10, 6)
    H5::LinkCreatPropList const lcpl;
    lcpl.setCreateIntermediateGroup(recursively);
    file.createGroup(groupPath.c_str(), lcpl);
#else
    auto groupCreationPropertyList = H5P_DEFAULT;

    if (recursively) {
        groupCreationPropertyList = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(groupCreationPropertyList, 1);
    }

    hid_t group = H5Gcreate(
        file.getId(), groupPath.c_str(), groupCreationPropertyList, H5P_DEFAULT,
        H5P_DEFAULT
    );
    H5Pclose(groupCreationPropertyList);

    if (group < 0)
        throw(AmiException(
            "Failed to create group in hdf5CreateGroup: %s", groupPath.c_str()
        ));
    H5Gclose(group);
#endif
}

std::unique_ptr<ExpData> read_exp_data_from_hdf5(
    std::string const& hdf5Filename, std::string const& hdf5Root,
    Model const& model
) {
    H5::H5File const file(hdf5Filename.c_str(), H5F_ACC_RDONLY);

    hsize_t m, n;

    auto edata = std::make_unique<ExpData>(model);

    if (attribute_exists(file, hdf5Root, "id")) {
        edata->id = get_string_attribute(file, hdf5Root, "id");
    }

    if (model.ny * model.nt() > 0) {
        if (location_exists(file, hdf5Root + "/Y")) {
            auto const my = get_double_2d_dataset(file, hdf5Root + "/Y", m, n);
            check_measurement_dimensions_compatible(m, n, model);
            edata->set_measurements(my);
        } else {
            throw AmiException(
                "Missing %s/Y in %s", hdf5Root.c_str(), hdf5Filename.c_str()
            );
        }

        if (location_exists(file, hdf5Root + "/Sigma_Y")) {
            auto const sigmay
                = get_double_2d_dataset(file, hdf5Root + "/Sigma_Y", m, n);
            check_measurement_dimensions_compatible(m, n, model);
            edata->set_noise_scales(sigmay);
        } else {
            throw AmiException(
                "Missing %s/Sigma_Y in %s", hdf5Root.c_str(),
                hdf5Filename.c_str()
            );
        }
    }

    if (model.nz * model.n_max_event() > 0) {
        if (location_exists(file, hdf5Root + "/Z")) {
            auto const mz = get_double_2d_dataset(file, hdf5Root + "/Z", m, n);
            check_event_dimensions_compatible(m, n, model);
            edata->set_event_measurements(mz);
        } else {
            throw AmiException(
                "Missing %s/Z in %s", hdf5Root.c_str(), hdf5Filename.c_str()
            );
        }

        if (location_exists(file, hdf5Root + "/Sigma_Z")) {
            auto sigmaz
                = get_double_2d_dataset(file, hdf5Root + "/Sigma_Z", m, n);
            check_event_dimensions_compatible(m, n, model);
            edata->set_event_noise_scales(sigmaz);
        } else {
            throw AmiException(
                "Missing %s/Sigma_Z in %s", hdf5Root.c_str(),
                hdf5Filename.c_str()
            );
        }
    }

    if (location_exists(file, hdf5Root + "/condition")) {
        edata->fixed_parameters
            = get_double_1d_dataset(file, hdf5Root + "/condition");
    }

    if (location_exists(file, hdf5Root + "/conditionPreequilibration")) {
        edata->fixed_parameters_pre_equilibration = get_double_1d_dataset(
            file, hdf5Root + "/conditionPreequilibration"
        );
    }

    if (location_exists(file, hdf5Root + "/conditionPresimulation")) {
        edata->fixed_parameters_presimulation
            = get_double_1d_dataset(file, hdf5Root + "/conditionPresimulation");
    }

    if (attribute_exists(file, hdf5Root, "t_presim")) {
        edata->t_presim
            = get_double_scalar_attribute(file, hdf5Root, "t_presim");
    }

    if (location_exists(file, hdf5Root + "/ts")) {
        edata->set_timepoints(get_double_1d_dataset(file, hdf5Root + "/ts"));
    }

    if (attribute_exists(
            file, hdf5Root, "/reinitializeFixedParameterInitialStates"
        )) {
        edata->reinitialize_fixed_parameter_initial_states
            = static_cast<bool>(get_int_scalar_attribute(
                file, hdf5Root, "/reinitializeFixedParameterInitialStates"
            ));
    }

    if (location_exists(file, hdf5Root + "/parameters")) {
        edata->free_parameters
            = get_double_1d_dataset(file, hdf5Root + "/parameters");
    }

    if (location_exists(file, hdf5Root + "/x0")) {
        edata->x0 = get_double_1d_dataset(file, hdf5Root + "/x0");
    }

    if (location_exists(file, hdf5Root + "/sx0")) {
        edata->sx0 = get_double_1d_dataset(file, hdf5Root + "/sx0");
    }

    if (location_exists(file, hdf5Root + "/pscale")) {
        auto pscaleInt = get_int_1d_dataset(file, hdf5Root + "/pscale");
        edata->pscale.resize(pscaleInt.size());
        for (int i = 0; (unsigned)i < pscaleInt.size(); ++i)
            edata->pscale[i] = static_cast<ParameterScaling>(pscaleInt[i]);
    }

    if (location_exists(file, hdf5Root + "/plist")) {
        edata->plist = get_int_1d_dataset(file, hdf5Root + "/plist");
    }

    if (location_exists(
            file, hdf5Root + "/reinitialization_state_idxs_presim"
        )) {
        edata->reinitialization_state_idxs_presim = get_int_1d_dataset(
            file, hdf5Root + "/reinitialization_state_idxs_presim"
        );
    }

    if (location_exists(file, hdf5Root + "/reinitialization_state_idxs_sim")) {
        edata->reinitialization_state_idxs_sim = get_int_1d_dataset(
            file, hdf5Root + "/reinitialization_state_idxs_sim"
        );
    }

    if (attribute_exists(file, hdf5Root, "tstart")) {
        edata->t_start = get_double_scalar_attribute(file, hdf5Root, "tstart");
    }

    if (attribute_exists(file, hdf5Root, "tstart_preeq")) {
        edata->t_start_preeq
            = get_double_scalar_attribute(file, hdf5Root, "tstart_preeq");
    }

    return edata;
}

void write_exp_data_to_hdf5(
    ExpData const& edata, H5::H5File const& file,
    std::string const& hdf5Location
) {

    if (!location_exists(file, hdf5Location))
        create_group(file, hdf5Location);

    H5LTset_attribute_string(
        file.getId(), hdf5Location.c_str(), "id", edata.id.c_str()
    );

    if (edata.nt())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/ts", edata.get_timepoints()
        );

    if (!edata.fixed_parameters.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/condition", edata.fixed_parameters
        );

    if (!edata.fixed_parameters_pre_equilibration.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/conditionPreequilibration",
            edata.fixed_parameters_pre_equilibration
        );

    if (!edata.fixed_parameters_presimulation.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/conditionPresimulation",
            edata.fixed_parameters_presimulation
        );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "t_presim", &edata.t_presim, 1
    );

    if (!edata.get_measurements().empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/Y", edata.get_measurements(), edata.nt(),
            edata.nytrue()
        );
    if (!edata.get_noise_scales().empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/Sigma_Y", edata.get_noise_scales(),
            edata.nt(), edata.nytrue()
        );
    if (!edata.get_event_measurements().empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/Z", edata.get_event_measurements(),
            edata.nmaxevent(), edata.nztrue()
        );
    if (!edata.get_event_noise_scales().empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/Sigma_Z",
            edata.get_event_noise_scales(), edata.nmaxevent(),
            edata.nztrue()
        );

    int int_attr = edata.reinitialize_fixed_parameter_initial_states;
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(),
        "reinitializeFixedParameterInitialStates", &int_attr, 1
    );

    if (!edata.free_parameters.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/parameters", edata.free_parameters
        );

    if (!edata.x0.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/x0", edata.x0
        );
    if (!edata.sx0.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/sx0", edata.sx0
        );

    std::vector<int> int_buffer;

    if (!edata.pscale.empty()) {
        int_buffer.resize(edata.pscale.size());
        for (int i = 0; (unsigned)i < edata.pscale.size(); i++)
            int_buffer[i] = static_cast<int>(edata.pscale[i]);
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/pscale", int_buffer
        );
    }

    if (!edata.plist.empty()) {
        int_buffer.resize(edata.plist.size());
        for (int i = 0; (unsigned)i < edata.plist.size(); i++)
            int_buffer[i] = static_cast<int>(edata.plist[i]);
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/plist", int_buffer
        );
    }

    if (!edata.reinitialization_state_idxs_presim.empty()) {
        int_buffer.resize(edata.reinitialization_state_idxs_presim.size());
        for (int i = 0;
             (unsigned)i < edata.reinitialization_state_idxs_presim.size(); i++)
            int_buffer[i]
                = static_cast<int>(edata.reinitialization_state_idxs_presim[i]);
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/reinitialization_state_idxs_presim",
            int_buffer
        );
    }

    if (!edata.reinitialization_state_idxs_sim.empty()) {
        int_buffer.resize(edata.reinitialization_state_idxs_sim.size());
        for (int i = 0;
             (unsigned)i < edata.reinitialization_state_idxs_sim.size(); i++)
            int_buffer[i]
                = static_cast<int>(edata.reinitialization_state_idxs_sim[i]);
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/reinitialization_state_idxs_sim", int_buffer
        );
    }

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "tstart", &edata.t_start, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "tstart_preeq",
        &edata.t_start_preeq, 1
    );
}

void write_return_data_to_hdf5(
    ReturnData const& rdata, H5::H5File const& file,
    std::string const& hdf5Location
) {

    if (!location_exists(file, hdf5Location))
        create_group(file, hdf5Location);

    if (!rdata.ts.empty())
        create_and_write_double_1d_dataset(file, hdf5Location + "/t", rdata.ts);

    H5LTset_attribute_string(
        file.getId(), hdf5Location.c_str(), "id", rdata.id.c_str()
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "llh", &rdata.llh, 1
    );
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "chi2", &rdata.chi2, 1
    );

    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "status", &rdata.status, 1
    );

    if (!rdata.sllh.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/sllh", rdata.sllh
        );

    if (!rdata.res.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/res", rdata.res
        );
    if (!rdata.sres.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/sres", rdata.sres, rdata.nt * rdata.nytrue,
            rdata.nplist
        );
    if (!rdata.FIM.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/FIM", rdata.FIM, rdata.nplist, rdata.nplist
        );

    if (!rdata.x0.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/x0", rdata.x0
        );

    if (!rdata.x.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/x", rdata.x, rdata.nt, rdata.nx_rdata
        );

    if (!rdata.y.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/y", rdata.y, rdata.nt, rdata.ny
        );

    if (!rdata.w.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/w", rdata.w, rdata.nt, rdata.nw
        );

    if (!rdata.z.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/z", rdata.z, rdata.nmaxevent, rdata.nz
        );

    if (!rdata.rz.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/rz", rdata.rz, rdata.nmaxevent, rdata.nz
        );

    if (!rdata.sigmay.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/sigmay", rdata.sigmay, rdata.nt, rdata.ny
        );

    if (!rdata.sigmaz.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/sigmaz", rdata.sigmaz, rdata.nmaxevent,
            rdata.nz
        );

    if (!rdata.s2llh.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/s2llh", rdata.s2llh, rdata.nJ - 1,
            rdata.nplist
        );

    if (!rdata.sx0.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/sx0", rdata.sx0, rdata.nplist, rdata.nx_rdata
        );

    if (!rdata.sx.empty())
        create_and_write_double_3d_dataset(
            file, hdf5Location + "/sx", rdata.sx, rdata.nt, rdata.nplist,
            rdata.nx_rdata
        );

    if (!rdata.sy.empty())
        create_and_write_double_3d_dataset(
            file, hdf5Location + "/sy", rdata.sy, rdata.nt, rdata.nplist,
            rdata.ny
        );

    if (!rdata.ssigmay.empty())
        create_and_write_double_3d_dataset(
            file, hdf5Location + "/ssigmay", rdata.ssigmay, rdata.nt,
            rdata.nplist, rdata.ny
        );

    if (!rdata.sz.empty())
        create_and_write_double_3d_dataset(
            file, hdf5Location + "/sz", rdata.sz, rdata.nmaxevent, rdata.nplist,
            rdata.nz
        );

    if (!rdata.srz.empty())
        create_and_write_double_3d_dataset(
            file, hdf5Location + "/srz", rdata.srz, rdata.nmaxevent,
            rdata.nplist, rdata.nz
        );

    if (!rdata.ssigmaz.empty())
        create_and_write_double_3d_dataset(
            file, hdf5Location + "/ssigmaz", rdata.ssigmaz, rdata.nmaxevent,
            rdata.nplist, rdata.nz
        );

    // TODO currently unused
    /*
    if (!rdata.s2rz.empty())
        createAndWriteDouble4DDataset(
            file, hdf5Location + "/s2rz", rdata.s2rz, rdata.nmaxevent,
            rdata.nztrue, rdata.nplist, rdata.nplist
        );
    */

    std::vector<int> int_buffer(1);

    int_buffer[0] = gsl::narrow<int>(rdata.newton_maxsteps);
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "newton_maxsteps",
        int_buffer.data(), 1
    );

    int_buffer[0] = static_cast<int>(rdata.o2mode);
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "o2mode", int_buffer.data(), 1
    );

    int_buffer[0] = static_cast<int>(rdata.sensi);
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "sensi", int_buffer.data(), 1
    );

    int_buffer[0] = static_cast<int>(rdata.sensi_meth);
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "sensi_meth", int_buffer.data(), 1
    );

    int_buffer[0] = static_cast<int>(rdata.rdata_reporting);
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "rdrm", int_buffer.data(), 1
    );

    if (!rdata.pscale.empty()) {
        int_buffer.resize(rdata.pscale.size());
        for (int i = 0; (unsigned)i < rdata.pscale.size(); i++)
            int_buffer[i] = static_cast<int>(rdata.pscale[i]);
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/pscale", int_buffer
        );
    }
    write_log_items_to_hdf5(file, rdata.messages, hdf5Location + "/messages");

    write_return_data_diagnosis(rdata, file, hdf5Location + "/diagnosis");
}

void write_return_data_diagnosis(
    ReturnData const& rdata, H5::H5File const& file,
    std::string const& hdf5Location
) {

    if (!location_exists(file, hdf5Location))
        create_group(file, hdf5Location);

    if (!rdata.xdot.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/xdot", rdata.xdot
        );

    if (!rdata.numsteps.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/numsteps", rdata.numsteps
        );

    if (!rdata.num_rhs_evals.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/numrhsevals", rdata.num_rhs_evals
        );

    if (!rdata.num_err_test_fails.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/numerrtestfails", rdata.num_err_test_fails
        );

    if (!rdata.num_non_lin_solv_conv_fails.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/numnonlinsolvconvfails",
            rdata.num_non_lin_solv_conv_fails
        );

    if (!rdata.order.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/order", rdata.order
        );

    if (!rdata.numsteps_b.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/numstepsB", rdata.numsteps_b
        );

    if (!rdata.num_rhs_evals_b.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/numrhsevalsB", rdata.num_rhs_evals_b
        );

    if (!rdata.num_err_test_fails_b.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/numerrtestfailsB", rdata.num_err_test_fails_b
        );

    if (!rdata.num_non_lin_solv_conv_fails_b.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/numnonlinsolvconvfailsB",
            rdata.num_non_lin_solv_conv_fails_b
        );

    if (!rdata.preeq_status.empty()) {
        std::vector<int> preeq_status_int(rdata.preeq_status.size());
        for (int i = 0; (unsigned)i < rdata.preeq_status.size(); i++)
            preeq_status_int[i] = static_cast<int>(rdata.preeq_status[i]);
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/preeq_status", preeq_status_int
        );
    }

    if (!rdata.preeq_numsteps.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/preeq_numsteps", rdata.preeq_numsteps
        );

    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "preeq_numstepsB",
        &rdata.preeq_numsteps_b, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "preeq_cpu_time",
        &rdata.preeq_cpu_time, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "preeq_cpu_timeB",
        &rdata.preeq_cpu_time_b, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "preeq_t", &rdata.preeq_t, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "preeq_wrms", &rdata.preeq_wrms, 1
    );

    if (!rdata.posteq_status.empty()) {
        std::vector<int> posteq_status_int(rdata.posteq_status.size());
        for (int i = 0; (unsigned)i < rdata.posteq_status.size(); i++)
            posteq_status_int[i] = static_cast<int>(rdata.posteq_status[i]);
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/posteq_status", posteq_status_int
        );
    }

    if (!rdata.posteq_numsteps.empty())
        create_and_write_int_1d_dataset(
            file, hdf5Location + "/posteq_numsteps", rdata.posteq_numsteps
        );

    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "posteq_numstepsB",
        &rdata.posteq_numsteps_b, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "posteq_cpu_time",
        &rdata.posteq_cpu_time, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "posteq_cpu_timeB",
        &rdata.posteq_cpu_time_b, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "posteq_t", &rdata.posteq_t, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "posteq_wrms", &rdata.posteq_wrms, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "cpu_time", &rdata.cpu_time, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "cpu_timeB", &rdata.cpu_time_b, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "cpu_time_total",
        &rdata.cpu_time_total, 1
    );

    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "t_last", &rdata.t_last, 1
    );

    if (!rdata.J.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/J", rdata.J, rdata.nx_rdata, rdata.nx_rdata
        );

    if (!rdata.x_ss.empty())
        create_and_write_double_1d_dataset(
            file, hdf5Location + "/x_ss", rdata.x_ss
        );

    if (!rdata.sx_ss.empty())
        create_and_write_double_2d_dataset(
            file, hdf5Location + "/sx_ss", rdata.sx_ss, rdata.nplist,
            rdata.nx_rdata
        );
}

// work-around for macos segfaults, use struct without std::string
struct LogItemCStr {
    int severity;
    char const* identifier;
    char const* message;
};

void write_log_items_to_hdf5(
    H5::H5File const& file, std::vector<LogItem> const& logItems,
    std::string const& hdf5Location
) {
    if (logItems.empty())
        return;

    try {
        hsize_t dims[1] = {logItems.size()};
        const H5::DataSpace dataspace(1, dims);

        // works on Ubuntu, but segfaults on macos:
        /*
        // Create a compound datatype for the LogItem struct.
        H5::CompType logItemType(sizeof(amici::LogItem));
        logItemType.insertMember(
            "severity", HOFFSET(amici::LogItem, severity),
            H5::PredType::NATIVE_INT
        );
        auto vlstr_type = H5::StrType(H5::PredType::C_S1, H5T_VARIABLE);
        logItemType.insertMember(
            "identifier", HOFFSET(amici::LogItem, identifier), vlstr_type
        );
        logItemType.insertMember(
            "message", HOFFSET(amici::LogItem, message), vlstr_type
        );
        H5::DataSet dataset
            = file.createDataSet(hdf5Location, logItemType, dataspace);

        dataset.write(logItems.data(), logItemType);
        */

        // ... therefore, as a workaround, we use a struct without std::string
        H5::CompType logItemType(sizeof(LogItemCStr));
        logItemType.insertMember(
            "severity", HOFFSET(LogItemCStr, severity), H5::PredType::NATIVE_INT
        );
        auto vlstr_type = H5::StrType(H5::PredType::C_S1, H5T_VARIABLE);
        logItemType.insertMember(
            "identifier", HOFFSET(LogItemCStr, identifier), vlstr_type
        );
        logItemType.insertMember(
            "message", HOFFSET(LogItemCStr, message), vlstr_type
        );
        H5::DataSet dataset
            = file.createDataSet(hdf5Location, logItemType, dataspace);

        // Convert std::vector<LogItem> to std::vector<LogItemCStr>
        std::vector<LogItemCStr> buffer(logItems.size());
        for (size_t i = 0; i < logItems.size(); ++i) {
            buffer[i].severity = static_cast<int>(logItems[i].severity);
            buffer[i].identifier = logItems[i].identifier.c_str();
            buffer[i].message = logItems[i].message.c_str();
        }

        // Write the data to the dataset.
        dataset.write(buffer.data(), logItemType);
    } catch (H5::Exception& e) {
        throw AmiException(e.getCDetailMsg());
    }
}

void write_return_data_to_hdf5(
    ReturnData const& rdata, std::string const& hdf5Filename,
    std::string const& hdf5Location
) {
    auto file = create_or_open_for_writing(hdf5Filename);

    write_return_data_to_hdf5(rdata, file, hdf5Location);
}

std::string get_string_attribute(
    H5::H5File const& file, std::string const& optionsObject,
    std::string const& attributeName
) {
    hsize_t dims;
    H5T_class_t type_class;
    size_t type_size;
    auto status = H5LTget_attribute_info(
        file.getId(), optionsObject.c_str(), attributeName.c_str(), &dims,
        &type_class, &type_size
    );
    if (status < 0) {
        throw AmiException(
            "Could get info for attribute %s for object %s.",
            attributeName.c_str(), optionsObject.c_str()
        );
    }
    std::vector<char> value(type_size);
    status = H5LTget_attribute_string(
        file.getId(), optionsObject.c_str(), attributeName.c_str(), value.data()
    );

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %s\n", attributeName.c_str(), value.data());
#endif

    if (status < 0)
        throw AmiException(
            "Attribute %s not found for object %s.", attributeName.c_str(),
            optionsObject.c_str()
        );

    return {value.data()};
}

double get_double_scalar_attribute(
    H5::H5File const& file, std::string const& optionsObject,
    std::string const& attributeName
) {

    double data = NAN;
    herr_t status = H5LTget_attribute_double(
        file.getId(), optionsObject.c_str(), attributeName.c_str(), &data
    );

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %e\n", attributeName.c_str(), data);
#endif

    if (status < 0)
        throw AmiException(
            "Attribute %s not found for object %s.", attributeName.c_str(),
            optionsObject.c_str()
        );

    return data;
}

int get_int_scalar_attribute(
    H5::H5File const& file, std::string const& optionsObject,
    std::string const& attributeName
) {
    int data = 0;
    herr_t status = H5LTget_attribute_int(
        file.getId(), optionsObject.c_str(), attributeName.c_str(), &data
    );

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d\n", attributeName.c_str(), data);
#endif

    if (status < 0)
        throw AmiException(
            "Attribute %s not found for object %s.", attributeName.c_str(),
            optionsObject.c_str()
        );

    return data;
}

void create_and_write_int_1d_dataset(
    H5::H5File const& file, std::string const& datasetName,
    gsl::span<int const> const buffer
) {
    hsize_t size = buffer.size();
    H5::DataSpace dataspace(1, &size);
    auto const dataset = file.createDataSet(
        datasetName.c_str(), H5::PredType::NATIVE_INT, dataspace
    );
    dataset.write(buffer.data(), H5::PredType::NATIVE_INT);
}

void create_and_write_double_1d_dataset(
    const H5::H5File& file, std::string const& datasetName,
    gsl::span<double const> buffer
) {
    hsize_t const size = buffer.size();
    H5::DataSpace dataspace(1, &size);
    auto const dataset = file.createDataSet(
        datasetName.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace
    );
    dataset.write(buffer.data(), H5::PredType::NATIVE_DOUBLE);
}

void create_and_write_double_2d_dataset(
    const H5::H5File& file, std::string const& datasetName,
    gsl::span<double const> const buffer, hsize_t const m, hsize_t const n
) {
    Expects(buffer.size() == m * n);
    hsize_t const adims[]{m, n};
    H5::DataSpace dataspace(2, adims);
    auto dataset = file.createDataSet(
        datasetName.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace
    );
    dataset.write(buffer.data(), H5::PredType::NATIVE_DOUBLE);
}

void create_and_write_int_2d_dataset(
    H5::H5File const& file, std::string const& datasetName,
    gsl::span<int const> const buffer, hsize_t const m, hsize_t const n
) {
    Expects(buffer.size() == m * n);
    hsize_t const adims[]{m, n};
    H5::DataSpace dataspace(2, adims);
    auto dataset = file.createDataSet(
        datasetName.c_str(), H5::PredType::NATIVE_INT, dataspace
    );
    dataset.write(buffer.data(), H5::PredType::NATIVE_INT);
}

void create_and_write_double_3d_dataset(
    H5::H5File const& file, std::string const& datasetName,
    gsl::span<double const> const buffer, hsize_t const m, hsize_t const n,
    hsize_t const o
) {
    Expects(buffer.size() == m * n * o);
    hsize_t const adims[]{m, n, o};
    H5::DataSpace dataspace(3, adims);
    auto dataset = file.createDataSet(
        datasetName.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace
    );
    dataset.write(buffer.data(), H5::PredType::NATIVE_DOUBLE);
}

bool attribute_exists(
    H5::H5File const& file, std::string const& optionsObject,
    std::string const& attributeName
) {
    AMICI_H5_SAVE_ERROR_HANDLER;
    int result = H5Aexists_by_name(
        file.getId(), optionsObject.c_str(), attributeName.c_str(), H5P_DEFAULT
    );
    AMICI_H5_RESTORE_ERROR_HANDLER;
    return result > 0;
}

bool attribute_exists(
    H5::H5Object const& object, std::string const& attributeName
) {
    AMICI_H5_SAVE_ERROR_HANDLER;
    int result = H5Aexists(object.getId(), attributeName.c_str());
    AMICI_H5_RESTORE_ERROR_HANDLER;
    return result > 0;
}

void write_solver_settings_to_hdf5(
    Solver const& solver, std::string const& hdf5Filename,
    std::string const& hdf5Location
) {
    auto const file = create_or_open_for_writing(hdf5Filename);

    write_solver_settings_to_hdf5(solver, file, hdf5Location);
}

void write_solver_settings_to_hdf5(
    Solver const& solver, H5::H5File const& file,
    std::string const& hdf5Location
) {
    if (!location_exists(file, hdf5Location))
        create_group(file, hdf5Location);

    double dbuffer;
    int ibuffer;

    dbuffer = solver.get_absolute_tolerance();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "atol", &dbuffer, 1
    );

    dbuffer = solver.get_relative_tolerance();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "rtol", &dbuffer, 1
    );

    dbuffer = solver.get_absolute_tolerance_fsa();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "atol_fsa", &dbuffer, 1
    );

    dbuffer = solver.get_relative_tolerance_fsa();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "rtol_fsa", &dbuffer, 1
    );

    dbuffer = solver.get_absolute_tolerance_b();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "atolB", &dbuffer, 1
    );

    dbuffer = solver.get_relative_tolerance_b();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "rtolB", &dbuffer, 1
    );

    dbuffer = solver.get_absolute_tolerance_quadratures();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "quad_atol", &dbuffer, 1
    );

    dbuffer = solver.get_relative_tolerance_quadratures();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "quad_rtol", &dbuffer, 1
    );

    dbuffer = solver.get_steady_state_tolerance_factor();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "ss_tol_factor", &dbuffer, 1
    );

    dbuffer = solver.get_absolute_tolerance_steady_state();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "ss_atol", &dbuffer, 1
    );

    dbuffer = solver.get_relative_tolerance_steady_state();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "ss_rtol", &dbuffer, 1
    );

    dbuffer = solver.get_steady_state_sensi_tolerance_factor();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "ss_tol_sensi_factor", &dbuffer, 1
    );

    dbuffer = solver.get_absolute_tolerance_steady_state_sensi();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "ss_atol_sensi", &dbuffer, 1
    );

    dbuffer = solver.get_relative_tolerance_steady_state_sensi();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "ss_rtol_sensi", &dbuffer, 1
    );

    dbuffer = solver.get_max_time();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "maxtime", &dbuffer, 1
    );

    dbuffer = solver.get_max_step_size();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "max_step_size", &dbuffer, 1
    );

    ibuffer = gsl::narrow<int>(solver.get_max_steps());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "maxsteps", &ibuffer, 1
    );

    ibuffer = gsl::narrow<int>(solver.get_max_steps_backward_problem());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "maxstepsB", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_linear_multistep_method());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "lmm", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_non_linear_solver_iteration());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "iter", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_stability_limit_flag());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "stldet", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_state_ordering());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "ordering", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_interpolation_type());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "interpType", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_sensitivity_method());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "sensi_meth", &ibuffer, 1
    );

    ibuffer
        = static_cast<int>(solver.get_sensitivity_method_pre_equilibration());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "sensi_meth_preeq", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_sensitivity_order());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "sensi", &ibuffer, 1
    );

    ibuffer = gsl::narrow<int>(solver.get_newton_max_steps());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "newton_maxsteps", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_newton_damping_factor_mode());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "newton_damping_factor_mode",
        &ibuffer, 1
    );

    dbuffer = solver.get_newton_damping_factor_lower_bound();
    H5LTset_attribute_double(
        file.getId(), hdf5Location.c_str(), "newton_damping_factor_lower_bound",
        &dbuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_linear_solver());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "linsol", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_internal_sensitivity_method());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "ism", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_return_data_reporting_mode());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "rdrm", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_newton_step_steady_state_check());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "newton_step_steadystate_conv",
        &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_sensi_steady_state_check());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "check_sensi_steadystate_conv",
        &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_max_nonlin_iters());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "max_nonlin_iters", &ibuffer, 1
    );

    ibuffer = static_cast<int>(solver.get_max_conv_fails());
    H5LTset_attribute_int(
        file.getId(), hdf5Location.c_str(), "max_conv_fails", &ibuffer, 1
    );

    create_and_write_double_1d_dataset(
        file, hdf5Location + "/constraints", solver.get_constraints()
    );
}

void read_solver_settings_from_hdf5(
    H5::H5File const& file, Solver& solver, std::string const& datasetPath
) {

    if (attribute_exists(file, datasetPath, "atol")) {
        solver.set_absolute_tolerance(
            get_double_scalar_attribute(file, datasetPath, "atol")
        );
    }

    if (attribute_exists(file, datasetPath, "rtol")) {
        solver.set_relative_tolerance(
            get_double_scalar_attribute(file, datasetPath, "rtol")
        );
    }

    if (attribute_exists(file, datasetPath, "atol_fsa")) {
        solver.set_absolute_tolerance_fsa(
            get_double_scalar_attribute(file, datasetPath, "atol_fsa")
        );
    }

    if (attribute_exists(file, datasetPath, "rtol_fsa")) {
        solver.set_relative_tolerance_fsa(
            get_double_scalar_attribute(file, datasetPath, "rtol_fsa")
        );
    }

    if (attribute_exists(file, datasetPath, "atolB")) {
        solver.set_absolute_tolerance_b(
            get_double_scalar_attribute(file, datasetPath, "atolB")
        );
    }

    if (attribute_exists(file, datasetPath, "rtolB")) {
        solver.set_relative_tolerance_b(
            get_double_scalar_attribute(file, datasetPath, "rtolB")
        );
    }

    if (attribute_exists(file, datasetPath, "quad_atol")) {
        solver.set_absolute_tolerance_quadratures(
            get_double_scalar_attribute(file, datasetPath, "quad_atol")
        );
    }

    if (attribute_exists(file, datasetPath, "quad_rtol")) {
        solver.set_relative_tolerance_quadratures(
            get_double_scalar_attribute(file, datasetPath, "quad_rtol")
        );
    }

    if (attribute_exists(file, datasetPath, "ss_tol_factor")) {
        solver.set_steady_state_tolerance_factor(
            get_double_scalar_attribute(file, datasetPath, "ss_tol_factor")
        );
    }

    if (attribute_exists(file, datasetPath, "ss_atol")) {
        solver.set_absolute_tolerance_steady_state(
            get_double_scalar_attribute(file, datasetPath, "ss_atol")
        );
    }

    if (attribute_exists(file, datasetPath, "ss_rtol")) {
        solver.set_relative_tolerance_steady_state(
            get_double_scalar_attribute(file, datasetPath, "ss_rtol")
        );
    }

    if (attribute_exists(file, datasetPath, "ss_tol_sensi_factor")) {
        solver.set_steady_state_sensi_tolerance_factor(
            get_double_scalar_attribute(
                file, datasetPath, "ss_tol_sensi_factor"
            )
        );
    }

    if (attribute_exists(file, datasetPath, "ss_atol_sensi")) {
        solver.set_absolute_tolerance_steady_state_sensi(
            get_double_scalar_attribute(file, datasetPath, "ss_atol_sensi")
        );
    }

    if (attribute_exists(file, datasetPath, "ss_rtol_sensi")) {
        solver.set_relative_tolerance_steady_state_sensi(
            get_double_scalar_attribute(file, datasetPath, "ss_rtol_sensi")
        );
    }

    if (attribute_exists(file, datasetPath, "maxtime")) {
        solver.set_max_time(
            get_double_scalar_attribute(file, datasetPath, "maxtime")
        );
    }

    if (attribute_exists(file, datasetPath, "max_step_size")) {
        solver.set_max_step_size(
            get_double_scalar_attribute(file, datasetPath, "max_step_size")
        );
    }

    if (attribute_exists(file, datasetPath, "maxsteps")) {
        solver.set_max_steps(
            get_int_scalar_attribute(file, datasetPath, "maxsteps")
        );
    }

    if (attribute_exists(file, datasetPath, "maxstepsB")) {
        solver.set_max_steps_backward_problem(
            get_int_scalar_attribute(file, datasetPath, "maxstepsB")
        );
    }

    if (attribute_exists(file, datasetPath, "lmm")) {
        solver.set_linear_multistep_method(
            static_cast<LinearMultistepMethod>(
                get_int_scalar_attribute(file, datasetPath, "lmm")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "iter")) {
        solver.set_non_linear_solver_iteration(
            static_cast<NonlinearSolverIteration>(
                get_int_scalar_attribute(file, datasetPath, "iter")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "stldet")) {
        solver.set_stability_limit_flag(
            get_int_scalar_attribute(file, datasetPath, "stldet")
        );
    }

    if (attribute_exists(file, datasetPath, "ordering")) {
        solver.set_state_ordering(
            get_int_scalar_attribute(file, datasetPath, "ordering")
        );
    }

    if (attribute_exists(file, datasetPath, "interpType")) {
        solver.set_interpolation_type(
            static_cast<InterpolationType>(
                get_int_scalar_attribute(file, datasetPath, "interpType")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "sensi_meth")) {
        solver.set_sensitivity_method(
            static_cast<SensitivityMethod>(
                get_int_scalar_attribute(file, datasetPath, "sensi_meth")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "sensi_meth_preeq")) {
        solver.set_sensitivity_method_pre_equilibration(
            static_cast<SensitivityMethod>(
                get_int_scalar_attribute(file, datasetPath, "sensi_meth_preeq")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "sensi")) {
        solver.set_sensitivity_order(
            static_cast<SensitivityOrder>(
                get_int_scalar_attribute(file, datasetPath, "sensi")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "newton_maxsteps")) {
        solver.set_newton_max_steps(
            get_int_scalar_attribute(file, datasetPath, "newton_maxsteps")
        );
    }

    if (attribute_exists(file, datasetPath, "newton_damping_factor_mode")) {
        solver.set_newton_damping_factor_mode(
            static_cast<NewtonDampingFactorMode>(get_int_scalar_attribute(
                file, datasetPath, "newton_damping_factor_mode"
            ))
        );
    }

    if (attribute_exists(
            file, datasetPath, "newton_damping_factor_lower_bound"
        )) {
        solver.set_newton_damping_factor_lower_bound(
            get_double_scalar_attribute(
                file, datasetPath, "newton_damping_factor_lower_bound"
            )
        );
    }

    if (attribute_exists(file, datasetPath, "linsol")) {
        solver.set_linear_solver(
            static_cast<LinearSolver>(
                get_int_scalar_attribute(file, datasetPath, "linsol")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "ism")) {
        solver.set_internal_sensitivity_method(
            static_cast<InternalSensitivityMethod>(
                get_int_scalar_attribute(file, datasetPath, "ism")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "rdrm")) {
        solver.set_return_data_reporting_mode(
            static_cast<RDataReporting>(
                get_int_scalar_attribute(file, datasetPath, "rdrm")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "newton_step_steadystate_conv")) {
        solver.set_newton_step_steady_state_check(get_int_scalar_attribute(
            file, datasetPath, "newton_step_steadystate_conv"
        ));
    }

    if (attribute_exists(file, datasetPath, "check_sensi_steadystate_conv")) {
        solver.set_sensi_steady_state_check(get_int_scalar_attribute(
            file, datasetPath, "check_sensi_steadystate_conv"
        ));
    }

    if (attribute_exists(file, datasetPath, "max_nonlin_iters")) {
        solver.set_max_nonlin_iters(
            get_int_scalar_attribute(file, datasetPath, "max_nonlin_iters")
        );
    }

    if (attribute_exists(file, datasetPath, "max_conv_fails")) {
        solver.set_max_conv_fails(
            get_int_scalar_attribute(file, datasetPath, "max_conv_fails")
        );
    }

    if (location_exists(file, datasetPath + "/constraints")) {
        solver.set_constraints(
            get_double_1d_dataset(file, datasetPath + "/constraints")
        );
    }
}

void read_solver_settings_from_hdf5(
    std::string const& hdffile, Solver& solver, std::string const& datasetPath
) {
    H5::H5File file(
        hdffile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT,
        H5::FileAccPropList::DEFAULT
    );

    read_solver_settings_from_hdf5(file, solver, datasetPath);
}

void read_model_data_from_hdf5(
    std::string const& hdffile, Model& model, std::string const& datasetPath
) {
    H5::H5File file(
        hdffile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT,
        H5::FileAccPropList::DEFAULT
    );

    read_model_data_from_hdf5(file, model, datasetPath);
}

void read_model_data_from_hdf5(
    const H5::H5File& file, Model& model, std::string const& datasetPath
) {
    if (attribute_exists(file, datasetPath, "tstart")) {
        model.set_t0(get_double_scalar_attribute(file, datasetPath, "tstart"));
    }

    if (attribute_exists(file, datasetPath, "tstart_preeq")) {
        model.set_t0_preeq(
            get_double_scalar_attribute(file, datasetPath, "tstart_preeq")
        );
    }

    if (location_exists(file, datasetPath + "/pscale")) {
        auto pscaleInt = get_int_1d_dataset(file, datasetPath + "/pscale");
        std::vector<ParameterScaling> pscale(pscaleInt.size());
        for (int i = 0; (unsigned)i < pscaleInt.size(); ++i)
            pscale[i] = static_cast<ParameterScaling>(pscaleInt[i]);
        model.set_parameter_scale(pscale);
    } else if (attribute_exists(file, datasetPath, "pscale")) {
        // if pscale is the same for all parameters,
        // it can be set as scalar attribute for convenience
        model.set_parameter_scale(
            static_cast<ParameterScaling>(
                get_double_scalar_attribute(file, datasetPath, "pscale")
            )
        );
    }

    if (attribute_exists(file, datasetPath, "nmaxevent")) {
        model.set_n_max_event(
            get_int_scalar_attribute(file, datasetPath, "nmaxevent")
        );
    }

    if (attribute_exists(file, datasetPath, "steadyStateComputationMode")) {
        model.set_steady_state_computation_mode(
            static_cast<SteadyStateComputationMode>(get_int_scalar_attribute(
                file, datasetPath, "steadyStateComputationMode"
            ))
        );
    }

    if (attribute_exists(file, datasetPath, "steadyStateSensitivityMode")) {
        model.set_steady_state_sensitivity_mode(
            static_cast<SteadyStateSensitivityMode>(get_int_scalar_attribute(
                file, datasetPath, "steadyStateSensitivityMode"
            ))
        );
    }

    if (location_exists(file, datasetPath + "/theta")) {
        model.set_free_parameters(
            get_double_1d_dataset(file, datasetPath + "/theta")
        );
    }

    if (location_exists(file, datasetPath + "/kappa")) {
        model.set_fixed_parameters(
            get_double_1d_dataset(file, datasetPath + "/kappa")
        );
    }

    if (location_exists(file, datasetPath + "/ts")) {
        model.set_timepoints(get_double_1d_dataset(file, datasetPath + "/ts"));
    }

    if (location_exists(file, datasetPath + "/sens_ind")) {
        auto sensInd = get_int_1d_dataset(file, datasetPath + "/sens_ind");
        model.set_parameter_list(sensInd);
    }

    if (location_exists(file, datasetPath + "/x0")) {
        auto x0 = get_double_1d_dataset(file, datasetPath + "/x0");
        if (!x0.empty())
            model.set_initial_state(x0);
    }

    if (location_exists(file, datasetPath + "/steadystate_mask")) {
        auto mask
            = get_double_1d_dataset(file, datasetPath + "/steadystate_mask");
        if (!mask.empty())
            model.set_steadystate_mask(mask);
    }

    if (location_exists(file, datasetPath + "/sx0")) {
        hsize_t length0 = 0;
        hsize_t length1 = 0;
        auto sx0 = get_double_2d_dataset(
            file, datasetPath + "/sx0", length0, length1
        );
        if (!sx0.empty()) {
            if (length0 != (unsigned)model.nplist()
                && length1 != (unsigned)model.nx_rdata)
                throw(AmiException(
                    "Dimension mismatch when reading sx0. "
                    "Expected %dx%d, got %llu, %llu.",
                    model.nx_rdata, model.nplist(), length0, length1
                ));
            model.set_unscaled_initial_state_sensitivities(sx0);
        }
    }

    if (attribute_exists(file, datasetPath, "sigma_res")) {
        auto sigma_res
            = get_int_scalar_attribute(file, datasetPath, "sigma_res");
        model.set_add_sigma_residuals(static_cast<bool>(sigma_res));
    }

    if (attribute_exists(file, datasetPath, "min_sigma")) {
        auto min_sigma
            = get_double_scalar_attribute(file, datasetPath, "min_sigma");
        model.set_minimum_sigma_residuals(min_sigma);
    }
}

H5::H5File create_or_open_for_writing(std::string const& hdf5filename) {
    AMICI_H5_SAVE_ERROR_HANDLER;
    try {
        H5::H5File file(hdf5filename.c_str(), H5F_ACC_RDWR);
        AMICI_H5_RESTORE_ERROR_HANDLER;
        return file;
    } catch (...) {
        AMICI_H5_RESTORE_ERROR_HANDLER;
        return H5::H5File(hdf5filename.c_str(), H5F_ACC_EXCL);
    }
}

bool location_exists(const H5::H5File& file, std::string const& location) {
    AMICI_H5_SAVE_ERROR_HANDLER;
    auto result = H5Lexists(file.getId(), location.c_str(), H5P_DEFAULT) > 0;
    AMICI_H5_RESTORE_ERROR_HANDLER;
    return result;
}

bool location_exists(std::string const& filename, std::string const& location) {
    H5::H5File file(filename.c_str(), H5F_ACC_RDONLY);
    return location_exists(file, location);
}

std::vector<int>
get_int_1d_dataset(const H5::H5File& file, std::string const& name) {
    auto dataset = file.openDataSet(name.c_str());
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1)
        throw(AmiException("Expected array of rank 1 in %s", name.c_str()));

    hsize_t dim;
    dataspace.getSimpleExtentDims(&dim);
    std::vector<int> result(dim);
    if (!result.empty())
        dataset.read(result.data(), H5::PredType::NATIVE_INT);
    return result;
}

std::vector<double>
get_double_1d_dataset(const H5::H5File& file, std::string const& name) {
    auto dataset = file.openDataSet(name.c_str());
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1)
        throw(AmiException("Expected array of rank 1 in %s", name.c_str()));

    hsize_t dim;
    dataspace.getSimpleExtentDims(&dim);
    std::vector<double> result(dim);
    if (!result.empty())
        dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;
}

std::vector<double> get_double_2d_dataset(
    const H5::H5File& file, std::string const& name, hsize_t& m, hsize_t& n
) {
    m = n = 0;

    auto dataset = file.openDataSet(name.c_str());
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if (rank != 2)
        throw(AmiException("Expected array of rank 2 in %s", name.c_str()));

    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims);
    m = dims[0];
    n = dims[1];

    std::vector<double> result(m * n);
    if (!result.empty())
        dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;
}

std::vector<double> get_double_3d_dataset(
    const H5::H5File& file, std::string const& name, hsize_t& m, hsize_t& n,
    hsize_t& o
) {
    m = n = o = 0;

    auto dataset = file.openDataSet(name.c_str());
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if (rank != 3)
        throw(AmiException("Expected array of rank 3 in %s", name.c_str()));

    hsize_t dims[3];
    dataspace.getSimpleExtentDims(dims);
    m = dims[0];
    n = dims[1];
    o = dims[2];

    std::vector<double> result(m * n * o);
    if (!result.empty())
        dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;
}

void write_exp_data_to_hdf5(
    ExpData const& edata, std::string const& filepath,
    std::string const& hdf5Location
) {
    auto const file = create_or_open_for_writing(filepath);

    write_exp_data_to_hdf5(edata, file, hdf5Location);
}

} // namespace amici::hdf5
