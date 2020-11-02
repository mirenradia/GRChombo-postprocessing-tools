/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "TimeData.hpp"
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

void TimeData::read_data(const std::string &a_filename, int a_block)
{
    m_file_reader.open(a_filename);
    m_file_reader.determine_file_structure();

    const auto &file_structure = m_file_reader.get_file_structure();

    m_num_datasets = file_structure.num_data_columns[a_block];
    m_num_steps = file_structure.num_data_rows[a_block];

    const int time_column = 0;
    std::vector<double> times = m_file_reader.get_column(time_column, a_block);
    m_time_limits.first = times[0];
    m_time_limits.second = times[m_num_steps - 1];
    m_dt = (m_time_limits.second - m_time_limits.first) / (m_num_steps - 1);

    m_data = std::move(m_file_reader.get_all_data_columns(a_block));

    m_data_read = true;
}

std::vector<double>
TimeData::read_data_from_header(const int a_header_row_number,
                                const int a_block)
{
    return m_file_reader.get_data_from_header(a_header_row_number, a_block);
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
