/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIO.hpp"
#include "StrainModeData.hpp"
#include "TimeData.hpp"
#include <chrono>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <limits>
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
std::vector<std::pair<int, int>> modes;
double cutoff_freq_low = 0.02;
bool rescale_cutoff = false;
double cutoff_freq_high = std::numeric_limits<double>::max();
double extraction_radius = 1.0; // the default is 1.0 because GRChombo output is
                                // already rescaled by extraction radius
} // namespace Options

void print_help(char *argv[])
{
    std::cout << "Usage: " << argv[0] << " [OPTION]... INPUT_FILE_PATTERN\n"
              << "INPUT_FILE_PATTERN should have {l} and {m} in place of mode "
                 "numbers.\n"
              << "-o, --output FILE_PREFIX output file prefix\n"
              << "                         default: /input/file/path/strain_\n"
              << "-c, --cutoff-low freq    low cutoff frequency (for l = 2 "
              << "if rescaling)\n"
              << "                         default: 0.02\n"
              << "-d, --cutoff-high freq   high cutoff frequency\n"
              << "                         default: DBL_MAX\n"
              << "-f, --rescale-cutoff     whether to rescale cutoff freqs"
              << "by m * freq / 2\n"
              << "                         default: false\n"
              << "-m, --modes MODES        modes to calculate\n"
              << "                         in form \"l1 m1 l2 m2 ...\" or "
              << "\"all\"\n"
              << "                         default: all\n"
              << "-r, --radius R_EX        output is multiplied by r_ex\n"
              << "                         default: 1\n"
              << "-h, --help               print help\n";
    exit(1);
}

void process_args(int argc, char *argv[], std::string &input_filename_pattern)
{
    const char *const short_opts = "o:c:d:fm:r:h";
    const option long_opts[] = {
        {"output", required_argument, nullptr, 'o'},
        {"cutoff-low", required_argument, nullptr, 'c'},
        {"cutoff-high", required_argument, nullptr, 'd'},
        {"rescale-cutoff", no_argument, nullptr, 'f'},
        {"modes", required_argument, nullptr, 'm'},
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
        case 'm':
        {
            std::string modes_string = optarg;
            // empty string or "all" means all possible pmodes are calculated
            if (!(modes_string.empty() || (modes_string == std::string("all"))))
            {
                std::stringstream modes_ss(modes_string);
                int l, m;
                while ((modes_ss >> l) && (modes_ss >> m))
                {
                    Options::modes.emplace_back(l, m);
                }
            }
            break;
        }
        case 'c':
        {
            Options::cutoff_freq_low = std::stod(optarg);
            std::cout << "cutoff_freq_low = " << Options::cutoff_freq_low
                      << "\n";
            break;
        }
        case 'd':
        {
            Options::cutoff_freq_high = std::stod(optarg);
            std::cout << "cutoff_freq_high = " << Options::cutoff_freq_high
                      << "\n";
            break;
        }
        case 'f':
        {
            Options::rescale_cutoff = true;
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
            "strain_";
    }
}

int main(int argc, char *argv[])
{
    using Clock = std::chrono::steady_clock;
    std::string input_filename_pattern;
    process_args(argc, argv, input_filename_pattern);

    StrainModeData strain_mode_data(input_filename_pattern, Options::modes);

    std::cout << "Reading Data\n";
    strain_mode_data.read_mode_data(Options::extraction_radius);

    if (Options::modes.empty())
    {
        Options::modes = strain_mode_data.get_available_modes();
    }

    auto start_timer = Clock::now();
    bool first_mode = true;
    for (const auto &mode : Options::modes)
    {
        const int l = mode.first;
        const int m = mode.second;
        std::cout << "Computing l = " << l << ", m = " << std::setw(2) << m
                  << "\r" << std::flush;
        double cutoff_freq_low = Options::cutoff_freq_low;
        if (Options::rescale_cutoff)
        {
            cutoff_freq_low *= static_cast<double>(std::abs(m)) / 2.0;
        }
        // make sure we don't divide by zero so force a minimum cutoff frequency
        cutoff_freq_low = std::max(1e-4, cutoff_freq_low);
        // compute this strain
        strain_mode_data.compute_strain(l, m, cutoff_freq_low,
                                        Options::cutoff_freq_high);
    }
    auto end_timer = Clock::now();
    std::chrono::duration<double, std::ratio<1>> time_taken =
        end_timer - start_timer;
    std::cout << "Time taken = " << time_taken.count() << "s\n";

    for (const auto &mode : Options::modes)
    {
        const int l = mode.first;
        const int m = mode.second;
        const auto &strain_data = strain_mode_data.get_strain_data(l, m);

        // write the cutoff frequency in the second header
        std::vector<std::string> header2(strain_data.size(), "");
        header2[0] =
            " = " + std::to_string(strain_mode_data.get_cutoff_freq_low(l, m));
        std::string filename_stem =
            Options::output_file_prefix + std::to_string(l) + std::to_string(m);

        strain_mode_data.m_time_data.write_data(filename_stem, strain_data,
                                                {"h_+", "h_x"}, header2, "time",
                                                "cutoff_freq");
    }
}
