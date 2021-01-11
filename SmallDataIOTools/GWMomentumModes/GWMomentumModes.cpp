/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIO.hpp"
#include "TimeData.hpp"
#include "WeylMomentumModeData.hpp"
//#include <chrono>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
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
bool sum_over_m = false; // output the sum over m (fixed l) of pmodes rather
                         // than invidivual pmodes
std::set<int> l_pmodes;  // values of l to output the sum over m pmodes for
bool sum_over_l =
    false; // also output partial sums over l for selected modes above
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
        << "-l, --l-pmodes[=L_VALS]  values of l to sum up pmodes over m\n"
        << "                         argument omitted means all (default)\n"
        << "-s, --sum-l              also write partial sums over l\n"
        << "                         (implies -l)\n"
        << "-r, --radius R_EX        output is multiplied by (r_ex)^2\n"
        << "                         default: 1\n"
        << "-h, --help               print help\n";
    exit(1);
}

void process_args(int argc, char *argv[], std::string &input_filename_pattern)
{
    const char *const short_opts = "o:p:l::sr:h";
    const option long_opts[] = {{"output", required_argument, nullptr, 'o'},
                                {"pmodes", required_argument, nullptr, 'p'},
                                {"l-pmodes", optional_argument, nullptr, 'l'},
                                {"sum-l", no_argument, nullptr, 's'},
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
            // empty string or "all" means all possible pmodes are calculated
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
        case 'l':
        {
            Options::sum_over_m = true;
            // empty string or "all" means all possible l_pmodes are calculated
            if (optarg != nullptr)
            {
                std::stringstream l_pmodes_ss(optarg);
                int l;
                while (l_pmodes_ss >> l)
                {
                    Options::l_pmodes.insert(l);
                }
            }
            break;
        }
        case 's':
        {
            Options::sum_over_m = true;
            Options::sum_over_l = true;
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

    WeylMomentumModeData weyl4_momentum_mode_data(input_filename_pattern);

    std::cout << "Reading Data\n";
    // auto start_timer = Clock::now();
    weyl4_momentum_mode_data.read_mode_data();
    // auto end_timer = Clock::now();
    // std::chrono::duration<double, std::ratio<1>> time_taken =
    //    end_timer - start_timer;
    // std::cout << "Time taken = " << time_taken.count() << "s\n";

    // only output individual pmodes if not asked to sum over m
    if (!Options::sum_over_m)
    {
        std::cout << "Writing individual pmodes.\n";
        // the empty set is the whole set i.e. no selected modes means compute
        // all of them
        if (Options::pmodes.empty())
        {
            // all available Psi4 modes are a superset of all computable pmodes
            Options::pmodes = weyl4_momentum_mode_data.get_available_modes();
        }

        for (const auto &pmode : Options::pmodes)
        {
            const int l = pmode.first;
            const int m = pmode.second;
            std::cout << "Computing l = " << l << ", m = " << std::setw(2) << m
                      << "\r" << std::flush;

            if (weyl4_momentum_mode_data.check_momentum_pmode_computable(l, m))
            {
                std::cout << "Computing l = " << l << ", m = " << std::setw(2)
                          << m << "\r" << std::flush;
                // compute momentum flux for this pmode in each direction
                const auto &pmode_flux =
                    weyl4_momentum_mode_data.compute_momentum_flux_pmode(
                        l, m, Options::extraction_radius);

                // integrate flux in time to get cumulative radiated momentum
                // for this pmode in each direction
                TimeData::time_multidata_t pmode;
                weyl4_momentum_mode_data.m_time_data.integrate_all_time(
                    pmode, pmode_flux);

                // write flux and cumulative data to files
                std::string filename_suffix =
                    "l" + std::to_string(l) + "m" + std::to_string(m);
                std::string flux_filename_stem =
                    Options::output_file_prefix + "flux_" + filename_suffix;
                std::string cumulative_filename_stem =
                    Options::output_file_prefix + filename_suffix;

                weyl4_momentum_mode_data.m_time_data.write_data(
                    flux_filename_stem, pmode_flux, {"x", "y", "z"});
                weyl4_momentum_mode_data.m_time_data.write_data(
                    cumulative_filename_stem, pmode, {"x", "y", "z"});
            }
        }
    }
    else
    {
        std::cout << "Writing pmodes summed over m\n";
        if (Options::l_pmodes.empty())
        {
            // just put all the l we have available. We'll check if it's
            // actually computable later
            const auto &available_modes =
                weyl4_momentum_mode_data.get_available_modes();
            for (const auto &mode : available_modes)
            {
                Options::l_pmodes.insert(mode.first);
            }
        }

        std::map<int, WeylModeData::mode_data_t> sum_m_pmode_data;
        for (int l : Options::l_pmodes)
        {
            if (weyl4_momentum_mode_data.check_momentum_pmode_computable(l))
            {
                sum_m_pmode_data.emplace(l, WeylModeData::mode_data_t());
                for (int m = -l; m <= l; ++m)
                {
                    std::cout << "Computing l = " << l
                              << ", m = " << std::setw(2) << m << "\r"
                              << std::flush;
                    const auto &pmode_flux =
                        weyl4_momentum_mode_data.compute_momentum_flux_pmode(
                            l, m, Options::extraction_radius);

                    // integrate flux in time to get cumulative radiated
                    // momentum for this pmode in each direction
                    WeylModeData::mode_data_t pmode_data;
                    weyl4_momentum_mode_data.m_time_data.integrate_all_time(
                        pmode_data, pmode_flux);
                    if (m == -l)
                    {
                        sum_m_pmode_data[l] = std::move(pmode_data);
                    }
                    else
                    {
                        const int num_columns = pmode_flux.size();
                        const int num_steps = pmode_flux[0].size();
                        auto &l_sum = sum_m_pmode_data[l];
                        for (int icol = 0; icol < num_columns; ++icol)
                        {
                            for (int istep = 0; istep < num_steps; ++istep)
                            {
                                l_sum[icol][istep] += pmode_data[icol][istep];
                            }
                        }
                    }
                }

                // open file for writing
                std::string sum_m_filename_stem =
                    Options::output_file_prefix + "l" + std::to_string(l);

                weyl4_momentum_mode_data.m_time_data.write_data(
                    sum_m_filename_stem, sum_m_pmode_data[l], {"x", "y", "z"});
            }
        }
        if (Options::sum_over_l)
        {
            // no need to write partial sum for l=2 as we've already done that
            // above
            WeylModeData::mode_data_t sum_l_data =
                sum_m_pmode_data.begin()->second;
            std::string filename_stem = Options::output_file_prefix + "l2";
            for (int l = 3; sum_m_pmode_data.find(l) != sum_m_pmode_data.end();
                 ++l)
            {
                const auto &this_sum_m_data = sum_m_pmode_data[l];
                const int num_columns = sum_l_data.size();
                const int num_steps = sum_l_data[0].size();
                for (int icol = 0; icol < num_columns; ++icol)
                {
                    for (int istep = 0; istep < num_steps; ++istep)
                    {
                        sum_l_data[icol][istep] += this_sum_m_data[icol][istep];
                    }
                }
                filename_stem += "+l" + std::to_string(l);

                weyl4_momentum_mode_data.m_time_data.write_data(
                    filename_stem, sum_l_data, {"x", "y", "z"});
            }
        }
    }
}
