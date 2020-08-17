/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIO.hpp"
#include "TimeData.hpp"
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
} // namespace Options

void print_help(char *argv[])
{
    std::cout << "Usage: " << argv[0] << " [OPTION]... INPUT_FILE\n"
              << "-o, --output FILE      output filename (w/o extension)\n"
              << "                       default: Norm_<input file>.dat\n"
              << "-h, --help             print help\n";
    exit(1);
}

void process_args(int argc, char *argv[], std::string &input_filename)
{
    const char *const short_opts = "o:";
    const option long_opts[] = {{"output", required_argument, nullptr, 'o'},
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
            parent_path + "Norm_" + input_path.stem().string();
    }
}

int main(int argc, char *argv[])
{
    std::string input_filename;
    process_args(argc, argv, input_filename);

    using Clock = std::chrono::steady_clock;

    TimeData time_data;
    time_data.read_data(input_filename);

    const auto &data = time_data.get_data();
    const int num_steps = data[0].size();

    auto norm_data = time_data.norm(data, CH_SPACEDIM);
    double dt = time_data.get_dt();

    const double fake_time = 0.0;
    // this will overwrite any existing files
    const bool first_step = true;
    SmallDataIO output_file(Options::output_filename, dt, fake_time, fake_time,
                            SmallDataIO::APPEND, first_step);
    output_file.write_header_line({});

    double min_time = time_data.get_time_limits().first;
    for (int istep = 0; istep < num_steps; ++istep)
    {
        double time = min_time + istep * dt;
        std::vector<double> data_for_writing(norm_data.size());
        for (int idataset = 0; idataset < data_for_writing.size(); ++idataset)
        {
            data_for_writing[idataset] = norm_data[idataset][istep];
        }
        output_file.write_data_line(data_for_writing, time);
    }
}
