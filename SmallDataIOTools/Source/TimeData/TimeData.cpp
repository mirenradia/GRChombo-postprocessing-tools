/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "TimeData.hpp"
#include "SmallDataIO.hpp"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

TimeData::TimeData() : m_data_read(false) {}

void TimeData::resize_time_data(time_data_t &a_time_data) const
{
    assert(m_data_read);

    a_time_data.resize(m_num_steps);
}

void TimeData::resize_time_multidata(time_multidata_t &a_time_multidata) const
{
    assert(m_data_read);
    a_time_multidata.resize(m_num_datasets);
    for (auto &time_data : a_time_multidata)
    {
        resize_time_data(time_data);
    }
}

void TimeData::read_data(const std::string &a_filename, const double start_time,
                         const double end_time, int a_block)
{
    m_file_reader.open(a_filename);
    m_file_reader.determine_file_structure();

    const auto &file_structure = m_file_reader.get_file_structure();

    m_num_datasets = file_structure.num_data_columns[a_block];
    const int total_num_steps = file_structure.num_data_rows[a_block];

    const int time_column = 0;
    std::vector<double> times = m_file_reader.get_column(time_column, a_block);
    m_dt = (times[total_num_steps - 1] - times[0]) / (total_num_steps - 1);

    int start_step = 0;
    while (times[start_step] < start_time)
    {
        ++start_step;
    }
    m_time_limits.first = times[start_step];

    int end_step = total_num_steps - 1;
    while (times[end_step] > end_time)
    {
        --end_step;
    }
    m_time_limits.second = times[end_step];

    m_num_steps = end_step - start_step + 1;

    const auto &all_data = m_file_reader.get_all_data_columns(a_block);
    m_data_read = true;
    resize_time_multidata(m_data);
    for (int icol = 0; icol < m_num_datasets; ++icol)
    {
        for (int istep = 0; istep < m_num_steps; ++istep)
        {
            m_data[icol][istep] = all_data[icol][istep + start_step];
        }
    }
}

std::vector<double>
TimeData::read_data_from_header(const int a_header_row_number,
                                const int a_block)
{
    int num_header_rows =
        m_file_reader.get_file_structure().num_header_rows[a_block];
    if (a_header_row_number < num_header_rows)
        return m_file_reader.get_data_from_header(a_header_row_number, a_block);
    else
        return {};
}

std::vector<std::string>
TimeData::read_header_strings(const int a_header_row_number, const int a_block)
{
    int num_header_rows =
        m_file_reader.get_file_structure().num_header_rows[a_block];

    if (a_header_row_number < num_header_rows)
        return m_file_reader.get_header_strings(a_header_row_number, a_block);
    else
        return {};
}

const TimeData::time_multidata_t &TimeData::get_data() const { return m_data; }
double TimeData::get_dt() const { return m_dt; }
const std::pair<double, double> &TimeData::get_time_limits() const
{
    return m_time_limits;
}

// Simpons's 3/8 rule for integration with odd number of points
const IntegrationMethod TimeData::simpson_odd({0.375, 1.125, 1.25});

double TimeData::integrate_time(const time_data_t &in_data, int a_min_step,
                                int a_max_step) const
{
    assert(m_data_read);
    assert(a_min_step <= a_max_step);
    const int num_intervals = a_max_step - a_min_step;
    const int num_steps = num_intervals + 1;
    constexpr bool is_periodic = false;

    double out = 0.0;

    if (num_intervals == 0)
    {
        return out;
    }
    else if (num_steps < 3) // Simpson's rule not possible -> use trapezium rule
    {
        for (int istep = a_min_step; istep <= a_max_step; ++istep)
        {
            out += m_dt *
                   IntegrationMethod::trapezium.weight(istep - a_min_step,
                                                       num_steps, is_periodic) *
                   in_data[istep];
        }
        return out;
    }
    else // Use Simpson's rule
    {
        if (num_intervals % 2 == 0) // even case
        {
            for (int istep = a_min_step; istep <= a_max_step; ++istep)
            {
                out += m_dt *
                       IntegrationMethod::simpson.weight(
                           istep - a_min_step, num_steps, is_periodic) *
                       in_data[istep];
            }
            return out;
        }
        else // odd case
        {
            // do first 3 steps using odd simpson formula
            const int num_odd_intervals = 3;
            const int num_odd_steps = num_odd_intervals + 1;
            const int odd_max_step = a_min_step + num_odd_intervals;
            for (int istep = a_min_step; istep <= odd_max_step; ++istep)
            {
                out += m_dt *
                       TimeData::simpson_odd.weight(
                           istep - a_min_step, num_odd_steps, is_periodic) *
                       in_data[istep];
            }
            // add weight again for last odd step (composition of Newton-Cotes
            // formulae)
            const int even_min_step = odd_max_step;
            const int num_even_steps = a_max_step - even_min_step + 1;
            for (int istep = even_min_step; istep <= a_max_step; ++istep)
            {
                out += m_dt *
                       IntegrationMethod::simpson.weight(
                           istep - even_min_step, num_even_steps, is_periodic) *
                       in_data[istep];
            }
            return out;
        }
    }
}

std::vector<double> TimeData::integrate_time(const time_multidata_t &in_data,
                                             int a_min_step,
                                             int a_max_step) const
{
    std::vector<double> out;
    out.reserve(in_data.size());

    for (const auto &time_data : in_data)
    {
        double integral = integrate_time(time_data, a_min_step, a_max_step);
        out.push_back(integral);
    }

    return out;
}

void TimeData::integrate_all_time(time_multidata_t &out,
                                  const time_multidata_t &in_data) const
{
    assert(m_data_read);
    out.clear();
    out.resize(in_data.size());
    for (auto &time_data : out)
    {
        resize_time_data(time_data);
    }

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(out, in_data, std::cout)
#endif
    for (int istep = 0; istep < m_num_steps; ++istep)
    {
        const auto integrals = integrate_time(in_data, 0, istep);
        for (int idataset = 0; idataset < in_data.size(); ++idataset)
        {
            out[idataset][istep] = integrals[idataset];
        }
    }
}

void TimeData::integrate_all_time(time_data_t &out,
                                  const time_data_t &in_data) const
{
    time_multidata_t out_vect;
    time_multidata_t in_data_vect(1, in_data);

    integrate_all_time(out_vect, in_data_vect);
    out = std::move(out_vect[0]);
}

TimeData::time_multidata_t TimeData::norm(const time_multidata_t &in_data,
                                          const int a_dim)
{
    const int num_columns_in = in_data.size();
    assert((num_columns_in > 0) && (num_columns_in % a_dim == 0));
    const int num_columns_out = num_columns_in / a_dim;
    const int num_steps = in_data[0].size();

    time_multidata_t out(num_columns_out);
    for (auto &time_data : out)
    {
        time_data.resize(num_steps);
    }

    for (int icol = 0; icol < num_columns_out; ++icol)
    {
        for (int istep = 0; istep < num_steps; ++istep)
        {
            double norm2 = 0.;
            for (int idir = 0; idir < a_dim; ++idir)
            {
                int icol_in = a_dim * icol + idir;
                norm2 += in_data[icol_in][istep] * in_data[icol_in][istep];
            }
            double norm = std::sqrt(norm2);
            out[icol][istep] = norm;
        }
    }
    return out;
}

TimeData::time_multidata_t TimeData::re_im_to_amp_phase(
    const time_multidata_t &in_data,
    const std::vector<int> &a_column_step_offsets) const
{
    const int num_columns = in_data.size();
    assert((num_columns > 0) && (num_columns % 2 == 0));
    std::vector<int> column_step_offsets =
        (a_column_step_offsets.size() > 0) ? a_column_step_offsets
                                           : std::vector<int>(num_columns, 0);
    assert(column_step_offsets.size() == num_columns);
    const int num_complex_datasets = num_columns / 2;
    const int num_steps = in_data[0].size();

    time_multidata_t out(num_columns);
    for (auto &time_data : out)
    {
        time_data.resize(num_steps);
    }

    for (int idata = 0; idata < num_complex_datasets; ++idata)
    {
        const int re_col = 2 * idata;
        const int amp_col = re_col;
        const int im_col = 2 * idata + 1;
        const int phase_col = im_col;
        int num_cycles = 0;
        for (int istep = 0; istep < num_steps; ++istep)
        {
            const double &re_data_now = in_data[re_col][istep];
            const double &im_data_now = in_data[im_col][istep];
            out[amp_col][istep] = std::sqrt(re_data_now * re_data_now +
                                            im_data_now * im_data_now);
            const double arg = std::atan2(im_data_now, re_data_now);
            if (istep <= a_column_step_offsets[re_col])
            {
                out[phase_col][istep] = arg;
            }
            else
            {
                const double phase_old = out[phase_col][istep - 1];
                const double phase_minus = arg + (num_cycles - 1) * 2.0 * M_PI;
                const double phase_orig = arg + num_cycles * 2.0 * M_PI;
                const double phase_plus = arg + (num_cycles + 1) * 2.0 * M_PI;
                const double diff_minus = std::abs(phase_minus - phase_old);
                const double diff_orig = std::abs(phase_orig - phase_old);
                const double diff_plus = std::abs(phase_plus - phase_old);
                if (diff_orig < diff_plus && diff_orig < diff_minus)
                {
                    out[phase_col][istep] = phase_orig;
                }
                else if (diff_plus < diff_minus)
                {
                    num_cycles += 1;
                    out[phase_col][istep] = phase_plus;
                }
                else
                {
                    num_cycles -= 1;
                    out[phase_col][istep] = phase_minus;
                }
            }
        }
    }
    return out;
}

void TimeData::write_data(const std::string &a_filename_stem,
                          const time_multidata_t &a_data,
                          const std::vector<std::string> &a_header1_strings,
                          const std::vector<std::string> &a_header2_strings,
                          const std::string &a_pre_header1_string,
                          const std::string &a_pre_header2_string,
                          const int a_data_precision,
                          const int a_time_precision)
{
    // these are redundant variables as we are using the SmallDataIO code from
    // GRChombo
    constexpr double fake_time = 0.0;
    constexpr double first_step = true;

    SmallDataIO file(a_filename_stem, m_dt, fake_time, fake_time,
                     SmallDataIO::APPEND, first_step, ".dat", a_data_precision,
                     a_time_precision);

    const int num_columns = a_data.size();
    assert(num_columns > 0 && a_data[0].size() == m_num_steps);

    const int header1_length = a_header1_strings.size();
    const int header2_length = a_header2_strings.size();
    std::vector<std::string> full_header1_strings(num_columns, "");
    std::vector<std::string> full_header2_strings(num_columns, "");
    if (header1_length > 0)
    {
        for (int icol = 0; icol < num_columns; ++icol)
        {
            full_header1_strings[icol] =
                a_header1_strings[icol % header1_length];
        }
    }
    if (header2_length > 0)
    {
        for (int icol = 0; icol < num_columns; ++icol)
        {
            full_header2_strings[icol] =
                a_header2_strings[icol % header2_length];
        }
    }
    file.write_header_line(full_header1_strings, a_pre_header1_string);
    file.write_header_line(full_header2_strings, a_pre_header2_string);

    for (int istep = 0; istep < m_num_steps; ++istep)
    {
        double time = m_time_limits.first + istep * m_dt;
        std::vector<double> data_for_writing(num_columns);
        for (int icol = 0; icol < num_columns; ++icol)
        {
            data_for_writing[icol] = a_data[icol][istep];
        }
        // can't use write_time_data_line as this will write the time as
        // fake_time
        file.write_data_line(data_for_writing, time);
    }
}

void TimeData::add_to_columns(time_multidata_t &data,
                              const std::vector<double> &a_values_to_add)
{
    const int num_columns = data.size();
    const int num_steps = data[0].size();
    assert(num_columns == a_values_to_add.size());
    for (int icol = 0; icol < num_columns; ++icol)
    {
        for (int istep = 0; istep < num_steps; ++istep)
        {
            data[icol][istep] += a_values_to_add[icol];
        }
    }
}

void TimeData::multiply_columns(time_multidata_t &data,
                                const std::vector<double> &a_factors)
{
    const int num_columns = data.size();
    const int num_steps = data[0].size();
    assert(num_columns == a_factors.size());
    for (int icol = 0; icol < num_columns; ++icol)
    {
        for (int istep = 0; istep < num_steps; ++istep)
        {
            data[icol][istep] *= a_factors[icol];
        }
    }
}

void TimeData::clear()
{
    m_data_read = false;
    m_file_reader.close();
    m_dt = 0.0;
    m_time_limits = {0.0, 0.0};
    m_num_datasets = 0;
    m_num_steps = 0;
    m_data.clear();
}
