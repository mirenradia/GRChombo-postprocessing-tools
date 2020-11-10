/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIO.hpp"
#include "TimeData.hpp"
#include "WeylModeData.hpp"
#include <cassert>
#include <chrono>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "IntegrationMethodSetup.hpp"

namespace Options
{
std::string output_file_prefix;
std::vector<std::pair<int, int>> modes;
double extraction_radius = 1.0;
} // namespace Options

void print_help(char *argv[])
{
    std::cout
        << "Usage: " << argv[0] << " [OPTION]... INPUT_FILE_PATTERN\n"
        << "INPUT_FILE_PATTERN should have {l} and {m}\n"
        << "in place of mode numbers.\n"
        << "-o, --output FILE_PREFIX output file prefix\n"
        << "                         default: /input/file/path/GW_momentum_\n"
        << "-m, --modes MODES        modes to calculate\n"
        << "                         in form \"l1 m1 l2 m2 ...\"\n"
        << "                         default: all\n"
        << "-r, --radius R_EX        output is multiplied by (r_ex)^2\n"
        << "                         default: 1\n"
        << "-h, --help               print help\n";
    exit(1);
}

void process_args(int argc, char *argv[], std::string &input_filename_pattern)
{
    const char *const short_opts = "o:m:r:";
    const option long_opts[] = {{"output", required_argument, nullptr, 'o'},
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
            // empty string or "all" means all possible modes are calculated
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

    std::filesystem::path input_path(input_filename_pattern);

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
    using Clock = std::chrono::steady_clock;
    std::string input_filename_pattern;
    process_args(argc, argv, input_filename_pattern);
    std::cout << "Options: \n"
              << "output prefix: " << Options::output_file_prefix << "\n"
              << "modes: \n";
    for (auto &mode : Options::modes)
    {
        std::cout << mode.first << " " << mode.second << "\n";
    }
    std::cout << "extraction radius: " << Options::extraction_radius << "\n"
              << "input file pattern: " << input_filename_pattern << "\n";

    WeylModeData weyl4_mode_data(input_filename_pattern);
    std::cout << "(l,m) = (2,2) file: "
              << weyl4_mode_data.get_mode_filepath(2, 2) << "\n";
    std::cout << "index: " << weyl4_mode_data.get_mode_index(2, 2) << "\n";

    std::cout << "\nReading Data\n";
    auto start_timer = Clock::now();
    weyl4_mode_data.read_mode_data();
    auto end_timer = Clock::now();
    std::chrono::duration<double, std::ratio<1>> time_taken =
        end_timer - start_timer;
    std::cout << "Time taken = " << time_taken.count() << "s\n";

    std::cout << "\nTime Integrating\n";
    start_timer = Clock::now();
    weyl4_mode_data.compute_time_integrated_mode_data();
    end_timer = Clock::now();
    time_taken = end_timer - start_timer;
    std::cout << "Time taken = " << time_taken.count() << "s\n";
}
