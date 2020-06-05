/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CH_assert.H"
#include "WeylExtractionData.hpp"
#include <chrono>
#include <getopt.h>
#include <limits>
#include <set>
#include <string>
//#include <iomanip>
#include <iostream>

#include "IntegrationMethodSetup.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Options
{
std::string prefix_path = "./Weyl4ExtractionOut_";
std::set<int> radial_indices;
bool compute_energy = true;
bool compute_momentum = true;
bool compute_angular_momentum = true;
std::array<int, 2> step_limits = {0, std::numeric_limits<int>::max()};
} // namespace Options

void print_help()
{
    std::cout
        << "-p, --prefix /path/to/prefix    prefix to ExtractionOut files\n"
        << "                                default: ./Weyl4ExtractionOut_\n"
        << "-r, --radii RADIAL INDICES      which radii to use\n"
        << "                                default: all\n"
        << "-c, --compute QUANTITIES        quantities to compute: \n"
        << "                                energy, linmom, angmom, all\n"
        << "                                default: all\n"
        << "-s, --step-limits MIN MAX       min and max steps to read\n"
        << "-h, --help                      print help\n";
    exit(1);
}

void process_args(int argc, char *argv[])
{
    const char *const short_opts = "p:r:c:s:h";
    const option long_opts[] = {
        {"prefix", required_argument, nullptr, 'p'},
        {"radii", required_argument, nullptr, 'r'},
        {"compute", required_argument, nullptr, 'c'},
        {"step-limits", required_argument, nullptr, 's'},
        {"help", no_argument, nullptr, 'h'}};

    while (true)
    {
        const auto opt =
            getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (opt == -1)
            break;

        switch (opt)
        {
        case 'p':
        {
            Options::prefix_path = optarg;
            break;
        }
        case 'r':
        {
            std::string radii_string = optarg;
            // empty string or "all" means all radii which we choose by making
            // Options::radial_indices be the empty set
            if (!(radii_string.empty() || (radii_string == std::string("all"))))
            {
                std::stringstream radii_ss(radii_string);
                int radial_index;
                while (radii_ss >> radial_index)
                {
                    Options::radial_indices.insert(radial_index);
                }
            }
            break;
        }
        case 'c':
        {
            std::string compute_string = optarg;
            // default case is to compute all anyway
            if (compute_string.find("all") != std::string::npos)
                break;
            else
            {
                Options::compute_energy = false;
                Options::compute_momentum = false;
                Options::compute_angular_momentum = false;
                if (compute_string.find("energy") != std::string::npos)
                    Options::compute_energy = true;
                if (compute_string.find("linmom") != std::string::npos)
                    Options::compute_momentum = true;
                if (compute_string.find("angmom") != std::string::npos)
                    Options::compute_angular_momentum = true;
                break;
            }
        }
        case 's':
        {
            std::string step_limits_str = optarg;
            std::stringstream step_limits_ss(step_limits_str);
            int step_limit;
            int i = 0;
            while (step_limits_ss >> step_limit)
            {
                Options::step_limits[i++] = step_limit;
            }
            break;
        }
        case 'h':
        case '?':
        default:
            print_help();
            break;
        }
    }
}

void print_options()
{
    std::cout << "Prefix: " << Options::prefix_path << "\n";
    std::cout << "Radial indices: ";
    if (Options::radial_indices.empty())
        std::cout << "all";
    else
    {
        for (int radial_index : Options::radial_indices)
        {
            std::cout << radial_index << " ";
        }
    }
    std::cout << "\nComputing: ";
    if (Options::compute_energy)
        std::cout << "energy  ";
    if (Options::compute_momentum)
        std::cout << "linear momentum  ";
    if (Options::compute_angular_momentum)
        std::cout << "angular momentum";
    std::cout << "\nStep Limits: min = " << Options::step_limits[0] << ", "
              << "max = " << Options::step_limits[1];
    std::cout << std::endl;
}

int main(int argc, char *argv[])
{
    process_args(argc, argv);
    print_options();

    using Clock = std::chrono::steady_clock;

    WeylExtractionData weyl_data(Options::prefix_path);

    // determine the layout of the files
    weyl_data.determine_data_structure(Options::radial_indices,
                                       Options::step_limits[0],
                                       Options::step_limits[1]);

    // read in the data
    auto start_timer = Clock::now();
    weyl_data.read_data();
    auto end_timer = Clock::now();
    std::chrono::duration<double, std::ratio<1>> time_taken =
        end_timer - start_timer;
    std::cout << "Time taken = " << time_taken.count() << "s\n";

    if (Options::compute_energy)
    {
        // Calculate power
        std::cout << "\nComputing power: \n";
        start_timer = Clock::now();
        weyl_data.compute_power();
        end_timer = Clock::now();
        time_taken = end_timer - start_timer;
        std::cout << "Time taken = " << time_taken.count() << "s\n";

        // Calculate total energy
        std::cout << "\nComputing total energy:\n";
        weyl_data.compute_energy();
    }

    if (Options::compute_momentum)
    {
        // Calculate momentum flux
        std::cout << "\nComputing momentum flux: \n";
        start_timer = Clock::now();
        weyl_data.compute_momentum_flux();
        end_timer = Clock::now();
        time_taken = end_timer - start_timer;
        std::cout << "Time taken = " << time_taken.count() << "s\n";

        // Calculate total momentum
        std::cout << "\nComputing total momentum:\n";
        weyl_data.compute_momentum();
    }

    if (Options::compute_angular_momentum)
    {
        // Calculate angular_momentum_flux
        std::cout << "\nComputing angular momentum flux: \n";
        start_timer = Clock::now();
        weyl_data.compute_angular_momentum_flux();
        end_timer = Clock::now();
        time_taken = end_timer - start_timer;
        std::cout << "Time taken = " << time_taken.count() << "s\n";

        // Calculate total angular momentum
        std::cout << "\nComputing total angular momentume:\n";
        weyl_data.compute_angular_momentum();
    }
}
