/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CH_assert.H"
#include "SmallDataIOReader.hpp"
#include "SurfaceExtractionData.hpp"
#include <chrono>
#include <iostream>

#include "IntegrationMethodSetup.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " "
                  << "</path/to/prefix>" << std::endl;
        return 0;
    }
    /*
    SmallDataIOReader file_reader;

    // open file
    std::cout << "Opening File\n";
    file_reader.open(argv[1]);

    // determine file structure
    std::cout << "Determining File structure\n";
    file_reader.determine_file_structure();

    // Get file structure
    auto file_structure = file_reader.get_file_structure();
    std::cout << "num_blocks = " << file_structure.num_blocks << "\n";
    std::cout << "block_starts = ";
    for (auto start : file_structure.block_starts)
    {
        std::cout << start << " ";
    }
    std::cout << "\nnum_header_rows = ";
    for (auto num : file_structure.num_header_rows)
    {
        std::cout << num << " ";
    }
    std::cout << "\nnum_data_rows = ";
    for (auto num : file_structure.num_data_rows)
    {
        std::cout << num << " ";
    }
    std::cout << "\nnum_coords_columns = ";
    for (auto num : file_structure.num_coords_columns)
    {
        std::cout << num << " ";
    }
    std::cout << "\ncoords_width = " << file_structure.coords_width;
    std::cout << "\nnum_data_columns = ";
    for (auto num : file_structure.num_data_columns)
    {
        std::cout << num << " ";
    }
    std::cout << "\ndata_width = " << file_structure.data_width << std::endl;

    auto header_vals = file_reader.get_data_from_header(0, 1);

    std::cout << "header_vals = ";
    for (auto vals : header_vals)
    {
        std::cout << vals << " ";
    }
    std::cout << std::endl;
    */
    using Clock = std::chrono::steady_clock;

    SurfaceExtractionData extraction_data(argv[1]);
    extraction_data.determine_data_structure();

    auto start_read = Clock::now();
    std::cout << "Reading in data:\n";
    extraction_data.read_data();
    auto end_read = Clock::now();
    std::chrono::duration<double, std::ratio<1>> time_taken =
        end_read - start_read;
    std::cout << "\nTime taken = " << time_taken.count() << "s\n";

    const auto &time_surface_data = extraction_data.get_data()[0][0];
    SurfaceExtractionData::time_surface_data_t integrated_time_surface_data;
    extraction_data.resize_time_surface_data(integrated_time_surface_data);

    auto start_loop = Clock::now();
    std::cout << "\nIntegrating in time:\n";
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(time_surface_data, integrated_time_surface_data, extraction_data)   \
        schedule(guided)
#endif
    for (int istep = 0; istep < integrated_time_surface_data.size(); ++istep)
    {
        // std::cout << "Step " << istep << "/"
        //          << integrated_time_surface_data.size() << "\r";
        extraction_data.integrate_time(integrated_time_surface_data[istep],
                                       time_surface_data, 0, istep);
    }
    auto end_loop = Clock::now();
    time_taken = end_loop - start_loop;
    std::cout << "\nTime taken = " << time_taken.count() << " s\n";
}
