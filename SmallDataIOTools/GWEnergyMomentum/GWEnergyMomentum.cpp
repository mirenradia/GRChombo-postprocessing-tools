/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CH_assert.H"
#include "SmallDataIOReader.hpp"
#include "SurfaceExtractionData.hpp"
#include <chrono>
#include <iomanip>
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

    using Clock = std::chrono::steady_clock;

    SurfaceExtractionData extraction_data(argv[1]);
    extraction_data.determine_data_structure();

    auto start_read = Clock::now();
    // std::cout << "Reading in data:\n";
    extraction_data.read_data();
    auto end_read = Clock::now();
    std::chrono::duration<double, std::ratio<1>> time_taken =
        end_read - start_read;
    std::cout << "\nTime taken = " << time_taken.count() << "s\n";

    const auto &data = extraction_data.get_data();
    SurfaceExtractionData::extracted_data_t integrated_data;
    integrated_data.resize(data.size());

    auto start_loop = Clock::now();
    std::cout << "\nIntegrating in time:\n";
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(data, integrated_data, extraction_data) schedule(guided)
#endif
    for (int istep = 0; istep < integrated_data.size(); ++istep)
    {
        // std::cout << "Step " << istep << "/"
        //          << integrated_time_surface_data.size() << "\r";
        extraction_data.integrate_time(integrated_data[istep], data, 0, istep);
    }
    auto end_loop = Clock::now();
    time_taken = end_loop - start_loop;
    std::cout << "\nTime taken = " << time_taken.count() << " s\n";

    /*
    for (int istep = 0; istep < integrated_data.size(); ++istep)
    {
        std::cout << std::setw(12) << std::fixed << std::setprecision(7)
                  << static_cast<double>(istep) / 999.0;
        std::cout << std::setw(20) << std::scientific << std::setprecision(12)
                  << integrated_data[istep][0][0][0][0] << "\n";
    }
    */
}
