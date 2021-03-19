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
#include <map>
#include <set>
#include <string>

#include "IntegrationMethodSetup.hpp"

namespace Options
{
std::string output_file_prefix;
int num_data_per_radius = 1;
bool use_tortoise_radii_for_retardation = false;
bool use_tortoise_radii_for_extrapolation = false;
double mass;
std::set<int> extrapolation_orders;
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
        << "                       tortoise radial coordinate \n"
        << "                       for retarding time if set)\n"
        << "-t  --tor-extr         use tortoise radius for extrapolation\n"
        << "                       only works if -m set\n"
        << "                       default: false\n"
        << "-r, --radii            which extraction radii to use\n"
        << "                       in form e.g. \"0110...\" where 1 is\n"
        << "                       include and 0 is exclude\n"
        << "                       default: \"111...\"\n"
        << "-e, --order            a space-separated list of the highest\n"
        << "                       order in the polynomial fit of r^{-n}\n"
        << "                       default: \"1 2\"\n"
        << "-h, --help             print help\n";
    exit(1);
}

void process_args(int argc, char *argv[], std::string &input_filename)
{
    const char *const short_opts = "o:n:m:tr:e:h";
    const option long_opts[] = {{"output", required_argument, nullptr, 'o'},
                                {"num-data", required_argument, nullptr, 'n'},
                                {"mass", required_argument, nullptr, 'm'},
                                {"tor-extr", no_argument, nullptr, 't'},
                                {"radii", required_argument, nullptr, 'r'},
                                {"order", required_argument, nullptr, 'e'},
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
            Options::use_tortoise_radii_for_retardation = true;
            Options::mass = std::stod(optarg);
            break;
        }
        case 't':
        {
            Options::use_tortoise_radii_for_extrapolation = true;
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
        case 'e':
        {
            std::string order_str = optarg;
            std::stringstream order_ss(order_str);
            int order;
            while (order_ss >> order)
            {
                Options::extrapolation_orders.insert(order);
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

    if (Options::use_tortoise_radii_for_extrapolation &&
        !Options::use_tortoise_radii_for_retardation)
        print_help(argv);

    input_filename = argv[optind];

    std::filesystem::path input_path(input_filename);

    if (Options::output_file_prefix.length() == 0)
    {
        Options::output_file_prefix = input_path.parent_path().string() + "/" +
                                      input_path.stem().string();
    }

    if (Options::extrapolation_orders.empty())
    {
        Options::extrapolation_orders = {1, 2};
    }
}

int main(int argc, char *argv[])
{
    std::string input_filename;
    process_args(argc, argv, input_filename);

    // Initial set up - set options here
    RetardedTimeData retarded_time_data;
    retarded_time_data.read_data(input_filename);
    retarded_time_data.set_num_data_cols_per_radii(
        Options::num_data_per_radius);
    retarded_time_data.set_radii(Options::use_tortoise_radii_for_retardation,
                                 Options::mass);
    retarded_time_data.set_radii_to_use(Options::use_radii);
    retarded_time_data.calculate_time_limits(
        Options::use_tortoise_radii_for_retardation);

    // retard the time coordinate and get the retarded data
    retarded_time_data.retard_data();
    const auto &retarded_data = retarded_time_data.get_retarded_data();

    // get the names of the column
    constexpr int column_names_row_number = 0;
    std::vector<std::string> column_names =
        retarded_time_data.read_header_strings(column_names_row_number);
    column_names.resize(Options::num_data_per_radius);
    std::cout << "Data column names are: ";
    for (const auto &name : column_names)
    {
        std::cout << name << " ";
    }
    std::cout << "\n";

    // get the chosen radii
    std::vector<double> chosen_radii = retarded_time_data.get_chosen_radii(
        Options::use_tortoise_radii_for_retardation);

    // Now do the extrapolation for all the requested orders
    std::map<int, TimeData::time_multidata_t> extrapolated_order_data_map;
    for (int order : Options::extrapolation_orders)
    {
        retarded_time_data.extrapolate_data(
            order, Options::use_tortoise_radii_for_extrapolation);
        extrapolated_order_data_map[order] =
            retarded_time_data.get_extrapolated_data();
    }

    const int m_num_retarded_steps = retarded_data[0].size();
    double dt = retarded_time_data.get_dt();
    const auto &retarded_time_limits =
        retarded_time_data.get_retarded_time_limits(
            Options::use_tortoise_radii_for_retardation);

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
    std::vector<std::string> retarded_file_header1_strings(
        retarded_data.size());
    std::vector<std::string> retarded_file_header2_strings(
        retarded_data.size());
    std::string retarded_file_pre_header2_string =
        (Options::use_tortoise_radii_for_retardation) ? "r* = " : "r = ";
    for (int icol = 0; icol < retarded_data.size(); ++icol)
    {
        retarded_file_header1_strings[icol] =
            column_names[icol % Options::num_data_per_radius];
        retarded_file_header2_strings[icol] =
            std::to_string(chosen_radii[icol / Options::num_data_per_radius]);
    }
    const int extrapolated_num_cols =
        Options::num_data_per_radius * Options::extrapolation_orders.size();
    std::vector<std::string> extrapolated_file_header1_strings(
        extrapolated_num_cols);
    std::vector<std::string> extrapolated_file_header2_strings(
        extrapolated_num_cols);
    auto order_it = Options::extrapolation_orders.begin();
    for (int icol = 0; icol < extrapolated_num_cols; ++icol)
    {
        extrapolated_file_header1_strings[icol] =
            column_names[icol % Options::num_data_per_radius];
        extrapolated_file_header2_strings[icol] = std::to_string(*order_it);
        if ((icol + 1) % Options::num_data_per_radius == 0)
            ++order_it;
    }

    retarded_output_file.write_header_line(retarded_file_header1_strings, "u");
    retarded_output_file.write_header_line(retarded_file_header2_strings,
                                           retarded_file_pre_header2_string);
    extrapolated_output_file.write_header_line(
        extrapolated_file_header1_strings, "u");
    extrapolated_output_file.write_header_line(
        extrapolated_file_header2_strings, "order = ");
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
            extrapolated_num_cols);
        int iorder = 0;
        for (int order : Options::extrapolation_orders)
        {
            const auto &this_order_data = extrapolated_order_data_map[order];
            for (int idataset = 0; idataset < Options::num_data_per_radius;
                 ++idataset)
            {
                int col = iorder * Options::num_data_per_radius + idataset;
                extrapolated_data_for_writing[col] =
                    this_order_data[idataset][istep];
            }
            ++iorder;
        }
        retarded_output_file.write_data_line(retarded_data_for_writing, time);
        extrapolated_output_file.write_data_line(extrapolated_data_for_writing,
                                                 time);
    }
}
