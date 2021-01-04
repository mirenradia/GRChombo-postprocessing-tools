/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIO.hpp"
#include "TimeData.hpp"
#include <chrono>
#include <cmath>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <set>
#include <string>

#include "IntegrationMethodSetup.hpp"

namespace Options
{
std::string output_filename;
double start_time;
} // namespace Options

void print_help(char *argv[])
{
    std::cout << "Usage: " << argv[0] << " [OPTION]... INPUT_FILE\n"
              << "-o, --output FILE      output filename (w/o extension)\n"
              << "                       default: <input file>_amp_phase.dat\n"
              << "-s, --start-time TIME  start time\n"
              << "                       default: 0.0\n"
              << "-h, --help             print help\n";
    exit(1);
}

void process_args(int argc, char *argv[], std::string &input_filename)
{
    const char *const short_opts = "o:s:";
    const option long_opts[] = {{"output", required_argument, nullptr, 'o'},
                                {"start-time", required_argument, nullptr, 's'},
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
            std::filesystem::path output_path(optarg);
            // remove file extension
            output_path.replace_extension();
            Options::output_filename = output_path.string();
            break;
        }
        case 's':
        {
            std::string start_time_str = optarg;
            Options::start_time = std::stod(start_time_str);
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

    input_filename = argv[optind];

    std::filesystem::path input_path(input_filename);

    if (Options::output_filename.length() == 0)
    {
        std::string parent_path = input_path.parent_path().string();
        if (parent_path.length() > 0)
        {
            parent_path += "/";
        }
        Options::output_filename =
            parent_path + input_path.stem().string() + "_amp_phase";
    }
}

int main(int argc, char *argv[])
{
    std::string input_filename;
    process_args(argc, argv, input_filename);

    using Clock = std::chrono::steady_clock;

    TimeData time_data;
    time_data.read_data(input_filename, Options::start_time);
    double dt = time_data.get_dt();

    const int extraction_radii_header_row = 1;
    const auto &extraction_radii =
        time_data.read_data_from_header(extraction_radii_header_row);
    std::vector<int> column_step_offsets(extraction_radii.size());
    for (int icol = 0; icol < column_step_offsets.size(); ++icol)
    {
        column_step_offsets[icol] = std::ceil(extraction_radii[icol] / dt);
    }

    const auto &data = time_data.get_data();
    const int num_steps = data[0].size();

    auto amp_phase_data =
        time_data.re_im_to_amp_phase(data, column_step_offsets);

    std::vector<std::string> extraction_radii_header_strings(
        extraction_radii.size());
    for (int icol = 0; icol < extraction_radii.size(); ++icol)
    {
        extraction_radii_header_strings[icol] =
            std::to_string(extraction_radii[icol]);
    }
    time_data.write_data(Options::output_filename, amp_phase_data,
                         {"amplitude", "phase"},
                         extraction_radii_header_strings, "time", "r = ");
}
