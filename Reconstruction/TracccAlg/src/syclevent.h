/* TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// SYCL include(s)
#include <CL/sycl.hpp>

// io
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/read_geometry.hpp"

// algorithms
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/spacepoint_formation.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "traccc/sycl/clusterization/clusterization_algorithm.hpp"
#include "traccc/sycl/seeding/seeding_algorithm.hpp"
#include "traccc/sycl/seeding/track_params_estimation.hpp"

// performance
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/performance/collection_comparator.hpp"
#include "traccc/performance/container_comparator.hpp"
#include "traccc/performance/timer.hpp"

// options
#include "traccc/options/common_options.hpp"
#include "traccc/options/detector_input_options.hpp"
#include "traccc/options/full_tracking_input_options.hpp"
#include "traccc/options/handle_argument_errors.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/sycl/device_memory_resource.hpp>
#include <vecmem/memory/sycl/host_memory_resource.hpp>
#include <vecmem/memory/sycl/shared_memory_resource.hpp>
#include <vecmem/utils/sycl/async_copy.hpp>

// Project include(s).
#include "traccc/utils/memory_resource.hpp"

// System include(s).
#include <exception>
#include <iomanip>
#include <iostream>

//memory resource
#include <memory_resource>

// csv parser
#include "csv2/writer.hpp"

namespace po = boost::program_options;

struct my_cell {
    uint64_t geometry_id = 0;
    uint64_t hit_id = 0;
    uint32_t channel0 = 0;
    uint32_t channel1 = 0;
    float timestamp = 0.;
    float value = 0.;
    uint64_t particle_id = 0;
};

const auto comp = [](const traccc::cell& c1, const traccc::cell& c2) {
    return c1.channel1 < c2.channel1;
};

namespace traccc{
/// Helper function which finds module from csv::cell in the geometry and
/// digitization config, and initializes the modules limits with the cell's
/// properties
cell_module my_get_module(my_cell c, const geometry* geom,
                          const digitization_config* dconfig) {

    cell_module result;
    result.surface_link = detray::geometry::barcode{c.geometry_id};

    // Find/set the 3D position of the detector module.
    if (geom != nullptr) {

        // Check if the module ID is known.
        if (!geom->contains(result.surface_link.value())) {
            throw std::runtime_error(
                "Could not find placement for geometry ID " +
                std::to_string(result.surface_link.value()));
        }

        // Set the value on the module description.
        result.placement = (*geom)[result.surface_link.value()];
    }

    // Find/set the digitization configuration of the detector module.
    if (dconfig != nullptr) {

        // Check if the module ID is known.
        const digitization_config::Iterator geo_it =
            dconfig->find(result.surface_link.value());
        if (geo_it == dconfig->end()) {
            throw std::runtime_error(
                "Could not find digitization config for geometry ID " +
                std::to_string(result.surface_link.value()));
        }

        // Set the value on the module description.
        const auto& binning_data = geo_it->segmentation.binningData();
        assert(binning_data.size() >= 2);
        result.pixel = {binning_data[0].min, binning_data[1].min,
                        binning_data[0].step, binning_data[1].step};
    }

    return result;
}

void my_read_cells(traccc::io::cell_reader_output& out, std::vector<my_cell> *cells,
                   const geometry* geom, const digitization_config* dconfig) {

    // Create cell counter vector.
    std::vector<unsigned int> cellCounts;
    cellCounts.reserve(5000);

    cell_module_collection_types::host& result_modules = out.modules;
    result_modules.reserve(5000);

    // Create a cell collection, which holds on to a flat list of all the cells
    // and the position of their respective cell counter & module.
    std::vector<std::pair<my_cell, unsigned int>> allCells;
    allCells.reserve(50000);

    for (const auto& iocell: *cells) {

        // Look for current module in cell counter vector.
        auto rit = std::find_if(result_modules.rbegin(), result_modules.rend(),
                                [&iocell](const cell_module& mod) {
                                    return mod.surface_link.value() ==
                                           iocell.geometry_id;
                                });
        if (rit == result_modules.rend()) {
            // Add new cell and new cell counter if a new module is found
            const cell_module mod = my_get_module(iocell, geom, dconfig);
            allCells.push_back({iocell, result_modules.size()});
            result_modules.push_back(mod);
            cellCounts.push_back(1);
        } else {
            // Add a new cell and update cell counter if repeat module is found
            const unsigned int pos =
                std::distance(result_modules.begin(), rit.base()) - 1;
            allCells.push_back({iocell, pos});
            ++(cellCounts[pos]);
        }
    }

    // Transform the cellCounts vector into a prefix sum for accessing
    // positions in the result vector.
    std::partial_sum(cellCounts.begin(), cellCounts.end(), cellCounts.begin());

    // The total number cells.
    const unsigned int totalCells = allCells.size();

    // Construct the result collection.
    cell_collection_types::host& result_cells = out.cells;
    result_cells.resize(totalCells);

    // Member "-1" of the prefix sum vector
    unsigned int nCellsZero = 0;
    // Fill the result object with the read csv cells
    for (unsigned int i = 0; i < totalCells; ++i) {
        const my_cell& c = allCells[i].first;

        // The position of the cell counter this cell belongs to
        const unsigned int& counterPos = allCells[i].second;

        unsigned int& prefix_sum_previous =
            counterPos == 0 ? nCellsZero : cellCounts[counterPos - 1];
        result_cells[prefix_sum_previous++] = traccc::cell{
            c.channel0, c.channel1, c.value, c.timestamp, counterPos, c.particle_id};
    }

    if (cellCounts.size() == 0) {
        return;
    }
    /* This is might look a bit overcomplicated, and could be made simpler by
     * having a copy of the prefix sum vector before incrementing its value when
     * filling the vector. however this seems more efficient, but requires
     * manually setting the 1st & 2nd modules instead of just the 1st.
     */

    // Sort the cells belonging to the first module.
    std::sort(result_cells.begin(), result_cells.begin() + nCellsZero, comp);
    // Sort the cells belonging to the second module.
    std::sort(result_cells.begin() + nCellsZero,
              result_cells.begin() + cellCounts[0], comp);

    // Sort cells belonging to all other modules.
    for (unsigned int i = 1; i < cellCounts.size() - 1; ++i) {
        std::sort(result_cells.begin() + cellCounts[i - 1],
                  result_cells.begin() + cellCounts[i], comp);
    }
}

}

class CEPCAlg{
public:
    virtual void initial(std::string* p_detect,
                         std::string* p_digit){}
    virtual void run(std::vector<my_cell> *cells, int cur_event){}
};

class ACTSTraccc : public CEPCAlg
{
private:
    // queue
    ::sycl::queue q;

    // file path
	std::string detector_file;
    std::string digitization_config_file;
    std::string input_directory;

    // geometry
    traccc::geometry surface_transforms;
    traccc::digitization_config digi_cfg;

    // modules
    // traccc::cell_collection_types::host cells_per_event;
    traccc::cell_module_collection_types::host modules_per_event;

    // config
    traccc::seedfinder_config finder_config;
    traccc::spacepoint_grid_config* grid_config;
    traccc::seedfilter_config filter_config;

    // memory resource
    vecmem::host_memory_resource host_mr;
    vecmem::sycl::host_memory_resource* sycl_host_mr;
    vecmem::sycl::device_memory_resource* device_mr;
    traccc::memory_resource* mr;
    vecmem::sycl::async_copy* copy;

    // algorithm
    traccc::sycl::clusterization_algorithm* ca_sycl;
    traccc::sycl::seeding_algorithm* sa_sycl;
    traccc::sycl::track_params_estimation* tp_sycl;

    // performance
    traccc::seeding_performance_writer* sd_performance_writer;
    traccc::performance::timing_info elapsedTimes;

public:
    void initial(std::string* p_detect,
                 std::string* p_digit);
    void run(std::vector<my_cell> *cells, int cur_event);
};


void ACTSTraccc::initial(std::string* p_detect, std::string* p_digit)
{
    // device
    std::cout << "Running Seeding on device: "
              << q.get_device().get_info<::sycl::info::device::name>() << "\n";

    std::cout << "Initializing ... " << "\n";
    detector_file            =   *p_detect;
    digitization_config_file =   *p_digit;

    // Read the surface transforms & digitization configuration file & input directory
    surface_transforms = traccc::io::read_geometry(detector_file);
    digi_cfg = traccc::io::read_digitization_config(digitization_config_file);

    // grid config
    grid_config  = new traccc::spacepoint_grid_config(finder_config);

    // memory resource
    sycl_host_mr = new vecmem::sycl::host_memory_resource{&q};
    device_mr    = new vecmem::sycl::device_memory_resource{&q};
    mr           = new traccc::memory_resource{*device_mr, sycl_host_mr};
    copy         = new vecmem::sycl::async_copy{&q};

    // algorithms
    ca_sycl = new traccc::sycl::clusterization_algorithm(*mr, *copy, &q, 1024);
    sa_sycl = new traccc::sycl::seeding_algorithm(finder_config, *grid_config,
                                               filter_config, *mr, *copy, &q);
    tp_sycl = new traccc::sycl::track_params_estimation(*mr, *copy, &q);

    // performance
    sd_performance_writer = new traccc::seeding_performance_writer(traccc::seeding_performance_writer::config{});
};

void ACTSTraccc::run(std::vector<my_cell> *cells, int cur_event)
{
    traccc::io::cell_reader_output read_out_per_event((*mr).host);
    traccc::my_read_cells(read_out_per_event, cells,
                          &surface_transforms, &digi_cfg);
    const traccc::cell_collection_types::host& cells_per_event =
        read_out_per_event.cells;
    const traccc::cell_module_collection_types::host&
        modules_per_event = read_out_per_event.modules;

    traccc::clusterization_algorithm::output_type measurements_per_event;
    traccc::spacepoint_formation::output_type spacepoints_per_event;
    traccc::seeding_algorithm::output_type seeds;
    traccc::track_params_estimation::output_type params;

    // Instantiate cuda containers/collections
    traccc::spacepoint_collection_types::buffer spacepoints_sycl_buffer(
        0, *(*mr).host);
    traccc::seed_double_collection_types::buffer seeds_sycl_buffer(0, *(*mr).host);
    traccc::bound_track_parameters_collection_types::buffer
        params_sycl_buffer(0, *(*mr).host);

    {
        traccc::performance::timer wall_t("Wall time", elapsedTimes);

        /*---------------------------------------------------
            Clusterization & Spacepoint formation
        ---------------------------------------------------*/

        // Create device copy of input collections
        traccc::cell_collection_types::buffer cells_buffer(
            cells_per_event.size(), (*mr).main);
        (*copy)(vecmem::get_data(cells_per_event), cells_buffer);
        traccc::cell_module_collection_types::buffer modules_buffer(
            modules_per_event.size(), (*mr).main);
        (*copy)(vecmem::get_data(modules_per_event), modules_buffer);

        {
            traccc::performance::timer t("Clusterization (sycl)",
                                            elapsedTimes);
            // Reconstruct it into spacepoints on the device.
            spacepoints_sycl_buffer =
                (*ca_sycl)(cells_buffer, modules_buffer).first;
            q.wait_and_throw();
        }

        /*----------------------------
            Seeding algorithm
        ----------------------------*/

        {
            traccc::performance::timer t("Seeding (sycl)", elapsedTimes);
            seeds_sycl_buffer = (*sa_sycl)(spacepoints_sycl_buffer);
            q.wait_and_throw();
        }

        /*----------------------------
            Track params estimation
        ----------------------------*/

        {
            traccc::performance::timer t("Track params (sycl)",
                                            elapsedTimes);
            params_sycl_buffer =
                (*tp_sycl)(spacepoints_sycl_buffer, seeds_sycl_buffer,
                        {0.f, 0.f, finder_config.bFieldInZ});
            q.wait_and_throw();
        }
    }

    traccc::spacepoint_collection_types::host spacepoints_per_event_sycl;
    traccc::seed_double_collection_types::host seeds_sycl;
    traccc::bound_track_parameters_collection_types::host params_sycl;
    (*copy)(spacepoints_sycl_buffer, spacepoints_per_event_sycl)->wait();
    (*copy)(seeds_sycl_buffer, seeds_sycl)->wait();
    (*copy)(params_sycl_buffer, params_sycl)->wait();

    // Show which event we are currently presenting the results for.
    std::cout << "===>>> Event " << cur_event << " <<<===" << std::endl;

    // write the spacepoints results to a csv file
    std::cout << "===>>> Writing spacepoints event " << cur_event << " to csv files <<<===" << std::endl;
    std::string spacepoints_filename = "obj/spacepoint/spacepoints_event" + std::to_string(cur_event) + ".csv";
    std::ofstream spacepoints_stream(spacepoints_filename);
    csv2::Writer<csv2::delimiter<','>> spacepoints_writer(spacepoints_stream);
    std::vector<std::string> spacepoints_header = {"x", "y", "z", "radius", "surface_link", "local_x", "local_y"};
    spacepoints_writer.write_row(spacepoints_header);
    const typename traccc::collection_types<traccc::spacepoint>::const_device spacepoints_coll{vecmem::get_data(spacepoints_per_event_sycl)};
    for (const traccc::spacepoint& obj : spacepoints_coll) {
        std::vector<std::string> spacepoints_col = {std::to_string(obj.x()), std::to_string(obj.y()), std::to_string(obj.z()),
                                        std::to_string(obj.radius()), std::to_string(obj.meas.surface_link.value()),
                                        std::to_string(obj.meas.local[0]), std::to_string(obj.meas.local[1])};
        spacepoints_writer.write_row(spacepoints_col);
    }
    spacepoints_stream.close();

    // write the seeds results to a csv file
    std::cout << "===>>> Writing seeds event " << cur_event << " to csv files <<<===" << std::endl;
    std::string seeds_filename = "obj/seed/seeds_event" + std::to_string(cur_event) + ".csv";
    std::ofstream seeds_stream(seeds_filename);
    csv2::Writer<csv2::delimiter<','>> seeds_writer(seeds_stream);
    std::vector<std::string> seeds_header = {"spB_inner", "spB_outer", "spM_inner", "spM_outer", "spT_inner", "spT_outer",
                                             "spB_inner_id", "spB_outer_id", "spM_inner_id", "spM_outer_id", "spT_inner_id", "spT_outer_id"};
    seeds_writer.write_row(seeds_header);
    const typename traccc::collection_types<traccc::seed_double>::const_device seeds_coll{vecmem::get_data(seeds_sycl)};
    traccc::collection_types<traccc::seed_double>::host sex_seeds;
    for (const traccc::seed_double& obj1 : seeds_coll) {
        unsigned int found = 0;
        for (traccc::seed_double& obj2 : sex_seeds){
            if ((obj1.spB_inner == obj2.spB_inner) & (obj1.spB_outer == obj2.spB_outer) &
                (obj1.spT_inner == obj2.spT_inner) & (obj1.spT_outer == obj2.spT_outer)){
                    found = 1;
                    if (spacepoints_coll[obj1.spM_inner].radius() < spacepoints_coll[obj2.spM_inner].radius()){
                        obj2.spM_inner = obj1.spM_inner;
                    } else {
                        obj2.spM_outer = obj1.spM_inner;
                    }
                }
        }
        if (found == 0) sex_seeds.push_back(obj1);
    }
    for (traccc::seed_double& obj : sex_seeds) {
        std::vector<std::string> seeds_col = {std::to_string(obj.spB_inner), std::to_string(obj.spB_outer), std::to_string(obj.spM_inner),
                                              std::to_string(obj.spM_outer), std::to_string(obj.spT_inner), std::to_string(obj.spT_outer),
                                              std::to_string(obj.spB_inner_id), std::to_string(obj.spB_outer_id), std::to_string(obj.spM_inner_id),
                                              std::to_string(obj.spM_outer_id), std::to_string(obj.spT_inner_id), std::to_string(obj.spT_outer_id)};
        seeds_writer.write_row(seeds_col);
    }
    seeds_stream.close();

    // write the track parameters results to a csv file
    std::cout << "===>>> Writing track parameters event " << cur_event << " to csv files <<<===" << std::endl;
    std::string params_filename = "obj/param/params_event" + std::to_string(cur_event) + ".csv";
    std::ofstream params_stream(params_filename);
    csv2::Writer<csv2::delimiter<','>> params_writer(params_stream);
    std::vector<std::string> params_header = {"surface_link", "local_x", "local_y", "phi", "theta", "qoverp"};
    params_writer.write_row(params_header);
    const typename traccc::collection_types<traccc::bound_track_parameters>::const_device params_coll{vecmem::get_data(params_sycl)};
    for (const traccc::bound_track_parameters& obj : params_coll) {
        std::vector<std::string> params_col = {std::to_string(obj.surface_link().value()),
                                        std::to_string(obj.bound_local()[0]), std::to_string(obj.bound_local()[1]),
                                        std::to_string(obj.phi()), std::to_string(obj.theta()), std::to_string(obj.qop())};
        params_writer.write_row(params_col);
    }
    params_stream.close();
};

CEPCAlg* factory(){
    return new ACTSTraccc();
}