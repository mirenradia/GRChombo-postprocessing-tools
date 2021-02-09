/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef RETARDEDTIMEDATA_HPP_
#define RETARDEDTIMEDATA_HPP_

#include "Misc.H"
#include "TimeData.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <iostream>

// A class to retard (by extraction radius) a collection of time series data
class RetardedTimeData : public TimeData
{
  protected:
    int m_num_data_cols_per_radii = 0;
    std::vector<double> m_radii;
    std::vector<double> m_tortoise_radii;
    bool m_radii_read = false;
    std::vector<bool> m_use_radii;
    std::vector<double> m_chosen_radii;
    std::vector<double> m_chosen_tortoise_radii;
    int m_num_radii_to_use;
    bool m_radii_chosen = false;
    std::pair<double, double> m_retarded_time_limits;
    bool m_time_limits_calculated = false;
    bool m_use_tortoise_radii_for_retardation;
    int m_num_retarded_steps;
    time_multidata_t m_retarded_data;
    time_multidata_t m_extrapolated_data;

  public:
    // constructor from TimeData
    using TimeData::TimeData;

    // Set retardation parameters, call before set_radii()
    void set_num_data_cols_per_radii(int a_num_data_cols_per_radius)
    {
        assert(m_data_read);
        m_num_data_cols_per_radii = a_num_data_cols_per_radius;
        std::cout << "Assuming " << m_num_data_cols_per_radii
                  << " data columns per radius which means there are\n"
                  << m_file_reader.get_file_structure().num_data_columns[0] /
                         m_num_data_cols_per_radii
                  << " radii in the file\n";
    }

    // calculate time limits, call after read_data() and
    // set_num_data_cols_per_radii()
    void set_radii(const bool a_use_tortoise_radii = false,
                   const double a_mass = 1)
    {
        assert(m_data_read);
        constexpr int radii_row = 1;
        std::vector<double> repeated_radii =
            m_file_reader.get_data_from_header(radii_row);

        m_radii.resize(repeated_radii.size() / m_num_data_cols_per_radii);
        for (int iradius = 0; iradius < m_radii.size(); ++iradius)
        {
            m_radii[iradius] =
                repeated_radii[m_num_data_cols_per_radii * iradius];
        }

        if (a_use_tortoise_radii)
        {
            std::cout
                << "Using Tortoise radial coordinate for time retardation\n"
                << "with mass = " << a_mass << "\n";
        }
        m_tortoise_radii.resize(m_radii.size());
        for (int iradius = 0; iradius < m_radii.size(); ++iradius)
        {
            double r = m_radii[iradius];
            m_tortoise_radii[iradius] =
                r + 2.0 * a_mass * std::log(r / (2.0 * a_mass) - 1);
        }
        m_radii_read = true;
    }

    // choose which radii to use, call after set_radii()
    void set_radii_to_use(const std::vector<bool> &a_use_radii)
    {
        assert(m_radii_read);
        assert(!m_radii_chosen);
        assert(a_use_radii.size() == m_radii.size() || a_use_radii.size() == 0);
        m_use_radii = (a_use_radii.size() != 0)
                          ? a_use_radii
                          : std::vector<bool>(m_radii.size(), true);
        std::cout << "Using data from coordinate radii = ";
        for (int iradius = 0; iradius < m_radii.size(); ++iradius)
        {
            if (m_use_radii[iradius])
            {
                m_chosen_radii.push_back(m_radii[iradius]);
                m_chosen_tortoise_radii.push_back(m_tortoise_radii[iradius]);
                std::cout << m_radii[iradius] << " ";
            }
        }
        std::cout << "\n";
        m_num_radii_to_use = m_chosen_radii.size();
        m_radii_chosen = true;
    }

    const std::vector<double> &
    get_chosen_radii(bool a_use_tortoise_radii = false)
    {
        assert(m_radii_chosen);
        return (a_use_tortoise_radii) ? m_chosen_tortoise_radii
                                      : m_chosen_radii;
    }

    // calculate time limits, call after set_radii_to_use()
    void
    calculate_time_limits(bool a_use_tortoise_radii_for_retardation = false)
    {
        assert(m_radii_chosen);
        m_use_tortoise_radii_for_retardation =
            a_use_tortoise_radii_for_retardation;
        const auto &chosen_radii = (m_use_tortoise_radii_for_retardation)
                                       ? m_chosen_tortoise_radii
                                       : m_chosen_radii;
        auto minmax_it =
            std::minmax_element(chosen_radii.begin(), chosen_radii.end());
        m_retarded_time_limits.first = m_time_limits.first - *(minmax_it.first);
        m_retarded_time_limits.second =
            m_time_limits.second - *(minmax_it.second);
        m_num_retarded_steps = std::floor((m_retarded_time_limits.second -
                                           m_retarded_time_limits.first) /
                                          m_dt) +
                               1;
        m_time_limits_calculated = true;
        std::cout << "Min retarded time = " << m_retarded_time_limits.first
                  << "\nMax retarded time = " << m_retarded_time_limits.second
                  << "\n";
    }

    const std::pair<double, double> &
    get_retarded_time_limits(bool a_use_tortoise_radii = false)
    {
        assert(m_time_limits_calculated);
        return m_retarded_time_limits;
    }

    // resample data onto common retarded time grid, call after
    // calculate_time_limits()
    void retard_data()
    {
        assert(m_time_limits_calculated);

        std::vector<double> retarded_time(m_num_steps);
        gsl_spline *spline_interp =
            gsl_spline_alloc(gsl_interp_cspline, m_num_steps);
        gsl_interp_accel *interp_acc = gsl_interp_accel_alloc();
        const auto &radii =
            (m_use_tortoise_radii_for_retardation) ? m_tortoise_radii : m_radii;

        for (int iradius = 0; iradius < m_radii.size(); ++iradius)
        {
            if (m_use_radii[iradius])
            {
                for (int istep = 0; istep < m_num_steps; ++istep)
                {
                    retarded_time[istep] =
                        m_time_limits.first + istep * m_dt - radii[iradius];
                }

                for (int icol = 0; icol < m_num_data_cols_per_radii; ++icol)
                {
                    gsl_interp_accel_reset(interp_acc);
                    gsl_spline_init(
                        spline_interp, retarded_time.data(),
                        m_data[iradius * m_num_data_cols_per_radii + icol]
                            .data(),
                        m_num_steps);
                    m_retarded_data.emplace_back(m_num_retarded_steps);
                    std::vector<double> &this_retarded_data_col =
                        m_retarded_data.back();
                    for (int irstep = 0; irstep < m_num_retarded_steps;
                         ++irstep)
                    {
                        double retarded_time =
                            m_retarded_time_limits.first + irstep * m_dt;
                        this_retarded_data_col[irstep] = gsl_spline_eval(
                            spline_interp, retarded_time, interp_acc);
                    }
                }
            }
        }
        gsl_spline_free(spline_interp);
        gsl_interp_accel_free(interp_acc);
    }

    const time_multidata_t &get_retarded_data()
    {
        assert(m_retarded_data.size() > 0);
        return m_retarded_data;
    }

    void extrapolate_data(int a_order, bool a_use_tortoise_radii = false)
    {
        assert(a_order >= 0);
        m_extrapolated_data.resize(m_num_data_cols_per_radii);
        for (auto &data_col : m_extrapolated_data)
        {
            data_col.resize(m_num_retarded_steps);
        }

        gsl_matrix *inv_radii =
            gsl_matrix_alloc(m_num_radii_to_use, a_order + 1);
        const auto &chosen_radii =
            (a_use_tortoise_radii) ? m_chosen_tortoise_radii : m_chosen_radii;
        for (int iradius = 0; iradius < m_num_radii_to_use; ++iradius)
        {
            double inv_radius = 1.0 / chosen_radii[iradius];
            for (int iord = 0; iord <= a_order; ++iord)
            {
                gsl_matrix_set(inv_radii, iradius, iord,
                               ipow(inv_radius, iord));
            }
        }

        gsl_vector *data = gsl_vector_alloc(m_num_radii_to_use);
        gsl_vector *coeff = gsl_vector_alloc(a_order + 1);
        gsl_matrix *covariance_matrix =
            gsl_matrix_alloc(a_order + 1, a_order + 1);
        gsl_multifit_linear_workspace *multifit_workspace =
            gsl_multifit_linear_alloc(m_num_radii_to_use, a_order + 1);

        for (int icol = 0; icol < m_num_data_cols_per_radii; ++icol)
        {
            for (int istep = 0; istep < m_num_retarded_steps; ++istep)
            {
                std::cout << istep << " / " << m_num_retarded_steps << "\r"
                          << std::flush;
                for (int iradius = 0; iradius < m_num_radii_to_use; ++iradius)
                {
                    gsl_vector_set(
                        data, iradius,
                        m_retarded_data[iradius * m_num_data_cols_per_radii +
                                        icol][istep]);
                }
                double chi_sq;
                gsl_multifit_linear(inv_radii, data, coeff, covariance_matrix,
                                    &chi_sq, multifit_workspace);
                m_extrapolated_data[icol][istep] = gsl_vector_get(coeff, 0);
            }
        }

        gsl_matrix_free(inv_radii);
        gsl_vector_free(data);
        gsl_vector_free(coeff);
        gsl_matrix_free(covariance_matrix);
        gsl_multifit_linear_free(multifit_workspace);
    }

    const time_multidata_t &get_extrapolated_data()
    {
        return m_extrapolated_data;
    }
};

#endif /* RETARDEDTIMEDATA_HPP_ */
