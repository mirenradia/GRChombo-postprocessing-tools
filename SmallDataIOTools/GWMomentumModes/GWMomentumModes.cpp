/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIO.hpp"
#include "TimeData.hpp"
#include "WeylModeData.hpp"
//#include <chrono>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#if __GNUC__ >= 8
#include <filesystem>
#elif __GNUC__ > 5 || (__GNUC__ == 5 && __GNUC_MINOR__ >= 3)
#include <experimental/filesystem>
#else
#error "This code requires GCC v5.3 or greater for the C++17 filesystem library"
#endif

#include "IntegrationMethodSetup.hpp"

namespace Options
{
std::string output_file_prefix;
std::vector<std::pair<int, int>>
    pmodes; // pmode refers to a single (l,m) term in Eq. (8.9.67) in Alcubierre
double extraction_radius = 1.0; // the default is 1.0 because GRChombo output is
                                // already rescaled by extraction radius
} // namespace Options

void print_help(char *argv[])
{
    std::cout
        << "Usage: " << argv[0] << " [OPTION]... INPUT_FILE_PATTERN\n"
        << "INPUT_FILE_PATTERN should have {l} and {m} in place of mode "
           "numbers.\n"
        << "-o, --output FILE_PREFIX output file prefix\n"
        << "                         default: /input/file/path/GW_momentum_\n"
        << "-p, --pmodes MODES       momentum flux pmodes to calculate\n"
        << "                         in form \"l1 m1 l2 m2 ...\" or "
        << "\"all\"\n"
        << "                         default: all\n"
        << "-r, --radius R_EX        output is multiplied by (r_ex)^2\n"
        << "                         default: 1\n"
        << "-h, --help               print help\n";
    exit(1);
}

void process_args(int argc, char *argv[], std::string &input_filename_pattern)
{
    const char *const short_opts = "o:p:r:h";
    const option long_opts[] = {{"output", required_argument, nullptr, 'o'},
                                {"pmodes", required_argument, nullptr, 'p'},
                                {"radius", required_argument, nullptr, 'r'},
                                {"help", no_argument, nullptr, 'h'}};

    while (true)
    {
        const auto opt =
            getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (opt == -1)
            break;

        switch (opt)
        {
        case 'o':
        {
            Options::output_file_prefix = optarg;
            break;
        }
        case 'p':
        {
            std::string pmodes_string = optarg;
            // empty string or "all" means all possible modes are calculated
            if (!(pmodes_string.empty() ||
                  (pmodes_string == std::string("all"))))
            {
                std::stringstream pmodes_ss(pmodes_string);
                int l, m;
                while ((pmodes_ss >> l) && (pmodes_ss >> m))
                {
                    Options::pmodes.emplace_back(l, m);
                }
            }
            break;
        }
        case 'r':
        {
            Options::extraction_radius = std::stod(optarg);
            break;
        }
        case 'h':
        case '?':
        default:
            print_help(argv);
            break;
        }
    }

    if (optind + 1 != argc)
        print_help(argv);

    input_filename_pattern = argv[optind];

#if __GNUC__ >= 8
    namespace fs = std::filesystem;
#elif __GNUC__ > 5 || (__GNUC__ == 5 && __GNUC_MINOR__ >= 3)
    namespace fs = std::experimental::filesystem;
#endif
    fs::path input_path(input_filename_pattern);

    if (Options::output_file_prefix.length() == 0)
    {
        Options::output_file_prefix =
            ((input_path.has_parent_path())
                 ? (input_path.parent_path().string() + "/")
                 : "") +
            "GW_momentum_";
    }
}

int main(int argc, char *argv[])
{
    // using Clock = std::chrono::steady_clock;
    std::string input_filename_pattern;
    process_args(argc, argv, input_filename_pattern);

    WeylModeData weyl4_mode_data(input_filename_pattern);

    std::cout << "\nReading Data\n";
    // auto start_timer = Clock::now();
    weyl4_mode_data.read_mode_data();
    // auto end_timer = Clock::now();
    // std::chrono::duration<double, std::ratio<1>> time_taken =
    //    end_timer - start_timer;
    // std::cout << "Time taken = " << time_taken.count() << "s\n";

    const double dt = weyl4_mode_data.m_time_data.get_dt();
    const double min_time = weyl4_mode_data.m_time_data.get_time_limits().first;

    // the empty set is the whole set i.e. no selected modes means compute all
    // of them
    if (Options::pmodes.empty())
    {
        // all available Psi4 modes are a superset of all computable pmodes
        Options::pmodes = weyl4_mode_data.get_available_modes();
    }

    for (const auto &pmode : Options::pmodes)
    {
        const int l = pmode.first;
        const int m = pmode.second;
        std::cout << "Computing l = " << l << ", m = " << std::setw(2) << m
                  << "\r" << std::flush;

        if (weyl4_mode_data.check_momentum_pmode_computable(l, m))
        {
            std::cout << "Computing l = " << l << ", m = " << std::setw(2) << m
                      << "\r" << std::flush;
            // compute momentum flux for this pmode in each direction
            const auto &pmode_flux =
                weyl4_mode_data.compute_momentum_flux_pmode(
                    l, m, Options::extraction_radius);

            // integrate flux in time to get cumulative radiated momentum for
            // this pmode in each direction
            TimeData::time_multidata_t pmode;
            weyl4_mode_data.m_time_data.integrate_all_time(pmode, pmode_flux);

            // don't bother taking the norm as it isn't very useful
            // take the norm/magnitude of the cumulative radiated momeentum
            // const auto pmode_norm = weyl4_mode_data.m_time_data.norm(pmode);

            // write flux, cumulative and norm data to files
            const int num_columns = pmode_flux.size();
            constexpr int dim = 3;
            const int num_extraction_radii = num_columns / dim;
            const int num_steps = pmode_flux[0].size();
            std::string filename_suffix =
                "l" + std::to_string(l) + "m" + std::to_string(m);
            std::string flux_filename =
                Options::output_file_prefix + "flux_" + filename_suffix;
            std::string cumulative_filename =
                Options::output_file_prefix + filename_suffix;
            // std::string norm_filename =
            //    Options::output_file_prefix + "norm_" + filename_suffix;

            // redundant variables for the SmallDataIO constructor as this is
            // ripped from GRChombo
            constexpr double fake_time = 0.0;
            constexpr double first_step = true;

            // open files for writing
            SmallDataIO flux_file(flux_filename, dt, fake_time, fake_time,
                                  SmallDataIO::APPEND, first_step);
            SmallDataIO cumulative_file(cumulative_filename, dt, fake_time,
                                        fake_time, SmallDataIO::APPEND,
                                        first_step);
            // SmallDataIO norm_file(norm_filename, dt, fake_time, fake_time,
            //                      SmallDataIO::APPEND, first_step);

            // write header in each file
            std::vector<std::string> header_strings(num_columns);
            for (int iradius = 0; iradius < num_extraction_radii; ++iradius)
            {
                header_strings[dim * iradius] = "x";
                header_strings[dim * iradius + 1] = "y";
                header_strings[dim * iradius + 2] = "z";
            }
            // note "time" is automatically in header for the first column
            flux_file.write_header_line(header_strings);
            cumulative_file.write_header_line(header_strings);
            // norm_file.write_header_line({});

            // now write the data
            for (int istep = 0; istep < num_steps; ++istep)
            {
                double time = min_time + istep * dt;
                std::vector<double> flux_data_for_writing(dim *
                                                          num_extraction_radii);
                std::vector<double> cumulative_data_for_writing(
                    dim * num_extraction_radii);
                // std::vector<double>
                // norm_data_for_writing(num_extraction_radii);

                for (int icol = 0; icol < num_columns; ++icol)
                {
                    flux_data_for_writing[icol] = pmode_flux[icol][istep];
                    cumulative_data_for_writing[icol] = pmode[icol][istep];
                }
                // for (int iradius = 0; iradius < num_extraction_radii;
                // ++iradius)
                //{
                //   norm_data_for_writing[iradius] =
                //   pmode_norm[iradius][istep];
                //}
                flux_file.write_data_line(flux_data_for_writing, time);
                cumulative_file.write_data_line(cumulative_data_for_writing,
                                                time);
                // norm_file.write_data_line(norm_data_for_writing, time);
            }
        }
    }
}
