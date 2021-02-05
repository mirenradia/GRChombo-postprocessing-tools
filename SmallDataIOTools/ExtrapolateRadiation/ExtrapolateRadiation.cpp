/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "RetardedTimeData.hpp"
#include "SmallDataIO.hpp"
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <set>
#include <string>

#include "IntegrationMethodSetup.hpp"

namespace Options
{
std::string output_file_prefix;
int num_data_per_radius = 1;
bool use_outgoing_ef_coord = false;
double mass;
std::vector<bool> use_radii;
} // namespace Options

void print_help(char *argv[])
{
    std::cout
        << "Usage: " << argv[0] << " [OPTION]... INPUT_FILE\n"
        << "-o, --output PREFIX    output file prefix (w/o extension)\n"
        << "                       default: <input file>\n"
        << "                       writes output_file_prefix_retarded.dat and\n"
        << "                       output_file_prefix_extrapolated.dat\n"
        << "-n, --num-data         number of data columns per radius\n"
        << "                       default: 1\n"
        << "-m, --mass             mass in code units (will use \n"
        << "                       outgoing Eddington-Finkelstein \n"
        << "                       coordinate for retarded time if set)\n"
        << "-r, --radii            which extraction radii to use\n"
        << "                       in form e.g. \"0110...\" where 1 is\n"
        << "                       include and 0 is exclude\n"
        << "                       default: \"111...\"\n"
        << "-h, --help             print help\n";
    exit(1);
}

void process_args(int argc, char *argv[], std::string &input_filename)
{
    const char *const short_opts = "o:n:m:r:h";
    const option long_opts[] = {{"output", required_argument, nullptr, 'o'},
                                {"num-data", required_argument, nullptr, 'n'},
                                {"mass", required_argument, nullptr, 'm'},
                                {"radii", required_argument, nullptr, 'r'},
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
            Options::output_file_prefix = output_path.string();
            break;
        }
        case 'n':
        {
            Options::num_data_per_radius = std::stoi(optarg);
            break;
        }
        case 'm':
        {
            Options::use_outgoing_ef_coord = true;
            Options::mass = std::stod(optarg);
            break;
        }
        case 'r':
        {
            std::string radii_str = optarg;
            Options::use_radii.resize(radii_str.size());
            for (int iradius = 0; iradius < Options::use_radii.size();
                 ++iradius)
            {
                char r_include = radii_str[iradius];
                if (r_include != '1' && r_include != '0')
                    print_help(argv);
                Options::use_radii[iradius] = (r_include == '1') ? true : false;
            }
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

    if (Options::output_file_prefix.length() == 0)
    {
        Options::output_file_prefix = input_path.parent_path().string() + "/" +
                                      input_path.stem().string();
    }
}

int main(int argc, char *argv[])
{
    std::string input_filename;
    process_args(argc, argv, input_filename);

    RetardedTimeData retarded_time_data;
    retarded_time_data.read_data(input_filename);
    retarded_time_data.set_num_data_cols_per_radii(
        Options::num_data_per_radius);
    retarded_time_data.set_radii(Options::use_outgoing_ef_coord, Options::mass);
    retarded_time_data.set_radii_to_use(Options::use_radii);
    retarded_time_data.calculate_time_limits();
    retarded_time_data.retard_data();
    retarded_time_data.extrapolate_data(1);

    const auto &retarded_data = retarded_time_data.get_retarded_data();
    const auto &extrapolated_data = retarded_time_data.get_extrapolated_data();
    const int m_num_retarded_steps = retarded_data[0].size();
    double dt = retarded_time_data.get_dt();
    const auto &retarded_time_limits =
        retarded_time_data.get_retarded_time_limits();

    const double fake_time = 0.0;
    // this will overwrite any existing files
    const bool first_step = true;
    std::string retarded_filename =
        Options::output_file_prefix + std::string("_retarded");
    SmallDataIO retarded_output_file(retarded_filename, dt, fake_time,
                                     fake_time, SmallDataIO::APPEND,
                                     first_step);
    std::string extrapolated_filename =
        Options::output_file_prefix + std::string("_extrapolated");
    SmallDataIO extrapolated_output_file(extrapolated_filename, dt, fake_time,
                                         fake_time, SmallDataIO::APPEND,
                                         first_step);

    retarded_output_file.write_header_line({}, "u");
    extrapolated_output_file.write_header_line({}, "u");
    for (int istep = 0; istep < m_num_retarded_steps; ++istep)
    {
        double time = retarded_time_limits.first + istep * dt;
        std::vector<double> retarded_data_for_writing(retarded_data.size());
        for (int idataset = 0; idataset < retarded_data_for_writing.size();
             ++idataset)
        {
            retarded_data_for_writing[idataset] =
                retarded_data[idataset][istep];
        }
        std::vector<double> extrapolated_data_for_writing(
            extrapolated_data.size());
        for (int idataset = 0; idataset < extrapolated_data_for_writing.size();
             ++idataset)
        {
            extrapolated_data_for_writing[idataset] =
                extrapolated_data[idataset][istep];
        }
        retarded_output_file.write_data_line(retarded_data_for_writing, time);
        extrapolated_output_file.write_data_line(extrapolated_data_for_writing,
                                                 time);
    }
}
