/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIO.hpp"
#include "TimeData.hpp"
#include <cassert>
#include <chrono>
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
double start_time = 0.0;
} // namespace Options

void print_help(char *argv[])
{
    std::cout << "Usage: " << argv[0] << " [OPTION]... INPUT_FILE\n"
              << "-o, --output FILE      output filename (w/o extension)\n"
              << "                       default: <input_file>_nojunk.dat\n"
              << "-s, --start TIME       start time of nonjunk data\n"
              << "                       default: 0\n"
              << "-h, --help             print help\n";
    exit(1);
}

void process_args(int argc, char *argv[], std::string &input_filename)
{
    const char *const short_opts = "o:s:";
    const option long_opts[] = {{"output", required_argument, nullptr, 'o'},
                                {"start", required_argument, nullptr, 's'},
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
            Options::start_time = std::stod(optarg);
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
        Options::output_filename = input_path.parent_path().string() + "/" +
                                   input_path.stem().string() + "_nojunk";
    }
}

int main(int argc, char *argv[])
{
    std::string input_filename;
    process_args(argc, argv, input_filename);

    using Clock = std::chrono::steady_clock;

    TimeData time_data;
    time_data.read_data(input_filename);

    auto data = time_data.get_data();
    const int num_columns = data.size();
    const int num_steps = data[0].size();
    double dt = time_data.get_dt();

    // assume extraction radii are in column
    constexpr int extraction_radii_header_row = 1;
    const auto extraction_radii =
        time_data.read_data_from_header(extraction_radii_header_row);
    assert(extraction_radii.size() == num_columns);

    // get the value at just before start_time + r_ex
    std::vector<double> junk_subtraction_values(num_columns);

    for (int icol = 0; icol < num_columns; ++icol)
    {
        double column_start_time = extraction_radii[icol] + Options::start_time;
        int column_subtraction_step = static_cast<int>(column_start_time / dt);
        junk_subtraction_values[icol] = -data[icol][column_subtraction_step];
    }

    time_data.add_to_columns(data, junk_subtraction_values);

    const double fake_time = 0.0;
    // this will overwrite any existing files
    const bool first_step = true;
    SmallDataIO output_file(Options::output_filename, dt, fake_time, fake_time,
                            SmallDataIO::APPEND, first_step);
    output_file.write_header_line({});

    std::vector<std::string> header2_strings(num_columns);
    for (int icol = 0; icol < num_columns; ++icol)
    {
        header2_strings[icol] = std::to_string(extraction_radii[icol]);
    }
    output_file.write_header_line(header2_strings, "r = ");

    double min_time = time_data.get_time_limits().first;
    for (int istep = 0; istep < num_steps; ++istep)
    {
        double time = min_time + istep * dt;
        std::vector<double> data_for_writing(num_columns);
        for (int icol = 0; icol < num_columns; ++icol)
        {
            data_for_writing[icol] = data[icol][istep];
        }
        output_file.write_data_line(data_for_writing, time);
    }
}
