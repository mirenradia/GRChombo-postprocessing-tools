/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEYLMODEDATA_HPP_
#define WEYLMODEDATA_HPP_

#include "SmallDataIO.hpp"
#include "TimeData.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#if __GNUC__ >= 8
#include <filesystem>
#elif __GNUC__ > 5 || (__GNUC__ == 5 && __GNUC_MINOR__ >= 3)
#include <experimental/filesystem>
#else
#error "This code requires GCC v5.3 or greater for the C++17 filesystem library"
#endif

#include <iostream>

//! This class reads, stores and manipulates modes of the Weyl scalar, Psi4
class WeylModeData
{
  public:
    using time_multidata_t = TimeData::time_multidata_t;
    using mode_data_t = time_multidata_t;
    TimeData m_time_data; // class to read and manipulate time data

  protected:
    std::string m_mode_file_parent_path;   // path to mode files
    std::string m_mode_file_pattern_begin; // mode filename pattern before l
    std::string
        m_mode_file_pattern_middle; // mode filename pattern between l and m
    std::string m_mode_file_pattern_end; // mode filename pattern after m
    std::vector<std::pair<int, int>>
        m_available_modes;                // modes available in directory
    std::vector<mode_data_t> m_mode_data; // stores the mode data
    mode_data_t
        m_zero_mode; // just an array of zeros to deal with cases such as m > l
    bool m_data_read = false;
    int m_num_modes;
    int m_num_extraction_radii;
    int m_num_steps;
    std::vector<mode_data_t>
        m_time_integrated_mode_data; // stores the time integrated mode data
    std::vector<bool>
        m_mode_data_time_integrated; // stores whether a mode has been time
                                     // integrated or not

    void find_available_modes(const std::string &a_mode_file_pattern)
    {
#if __GNUC__ >= 8
        namespace fs = std::filesystem;
#elif __GNUC__ > 5 || (__GNUC__ == 5 && __GNUC_MINOR__ >= 3)
        namespace fs = std::experimental::filesystem;
#endif

        fs::path pattern_path(a_mode_file_pattern);
        // get just the filename without the path
        std::string pattern = pattern_path.filename();

        // assume l_pos < m_pos
        const auto l_pos_pattern = pattern.find("{l}");
        m_mode_file_pattern_begin = pattern.substr(0, l_pos_pattern);
        const auto m_pos_pattern = pattern.find("{m}");
        m_mode_file_pattern_middle = pattern.substr(
            l_pos_pattern + 3, m_pos_pattern - l_pos_pattern - 3);
        m_mode_file_pattern_end = pattern.substr(m_pos_pattern + 3);

        // set directory to search for mode files
        if (pattern_path.has_parent_path())
        { // if provided pattern has a parent_path use that
            m_mode_file_parent_path = pattern_path.parent_path();
        }
        else
        { // else use the current directory
            m_mode_file_parent_path = fs::current_path();
        }

        // iterate through all files in the directory
        for (auto &file :
             fs::directory_iterator(fs::path(m_mode_file_parent_path)))
        {
            // get the filename of the file that the iterator points to
            const std::string &filename = file.path().filename();
            int l, m;
            int l_pos, m_pos;
            // check if the beginning of the filename matches the mode file
            // pattern begin
            if (filename.compare(0, m_mode_file_pattern_begin.length(),
                                 m_mode_file_pattern_begin) == 0)
            {
                l_pos = m_mode_file_pattern_begin.length();
                // assume l is only 1 digit
                l = std::stoi(filename.substr(l_pos, 1));
                // check if the filename matches the pattern in between l and m
                m_pos = filename.find(m_mode_file_pattern_middle, l_pos + 1) +
                        m_mode_file_pattern_middle.length();
                if (m_pos != std::string::npos)
                {
                    // check if the end of the filename matches the pattern
                    const int m_pos_end =
                        filename.find(m_mode_file_pattern_end, m_pos + 1);

                    if (m_pos_end != std::string::npos)
                    {
                        m = std::stoi(
                            filename.substr(m_pos, m_pos_end - m_pos));
                        m_available_modes.emplace_back(l, m);
                    }
                }
            }
        }
        m_num_modes = m_available_modes.size();
    }

  public:
    WeylModeData(const std::string &a_mode_file_pattern)
    {
        find_available_modes(a_mode_file_pattern);
    }

    // returns the full path to the file for the (l,m) mode data
    std::string get_mode_filepath(const int l, const int m)
    {
        return m_mode_file_parent_path + "/" + m_mode_file_pattern_begin +
               std::to_string(l) + m_mode_file_pattern_middle +
               std::to_string(m) + m_mode_file_pattern_end;
    }

    // returns the index in the m_available_modes vector corresponding to the
    // (l,m) mode or -1 otherwise
    int get_mode_index(const int l, const int m) const
    {
        auto mode_it = std::find(m_available_modes.begin(),
                                 m_available_modes.end(), std::make_pair(l, m));
        if (mode_it == m_available_modes.end())
        {
            return -1;
        }
        return std::distance(m_available_modes.begin(), mode_it);
    }

    // read all modes
    void read_mode_data()
    {
        assert(!m_data_read);
        m_mode_data.resize(m_num_modes);

        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            m_time_data.clear();
            const auto &mode = m_available_modes[imode];
            const std::string filepath =
                get_mode_filepath(mode.first, mode.second);
            m_time_data.read_data(filepath);
            m_mode_data[imode] = std::move(m_time_data.get_data());
        }
        if (m_num_modes > 0)
        {
            m_num_extraction_radii = m_mode_data[0].size() / 2;
            m_num_steps = m_mode_data[0][0].size();
        }
        m_zero_mode.resize(m_num_extraction_radii * 2);
        for (auto &dataset : m_zero_mode)
        {
            dataset = std::move(std::vector<double>(m_num_steps, 0));
        }

        m_time_integrated_mode_data.resize(m_num_modes);
        m_mode_data_time_integrated =
            std::move(std::vector<bool>(m_num_modes, false));
        m_data_read = true;
    }

    // integrates the data for a single mode in time
    // does nothing if already integrated
    void time_integrate_mode(const int a_idx_mode)
    {
        assert(m_data_read);
        if (!m_mode_data_time_integrated[a_idx_mode])
        {
            m_time_data.integrate_all_time(
                m_time_integrated_mode_data[a_idx_mode],
                m_mode_data[a_idx_mode]);
            m_mode_data_time_integrated[a_idx_mode] = true;
        }
    }

    // getters
    const std::vector<mode_data_t> &get_mode_data() const
    {
        assert(m_data_read);
        return m_mode_data;
    }

    const mode_data_t &get_mode_data(const int l, const int m) const
    {
        if (abs(m) > l || l < 2)
        {
            return m_zero_mode;
        }
        assert(m_data_read);
        const int idx_mode = get_mode_index(l, m);
        assert(idx_mode >= 0);
        return m_mode_data[idx_mode];
    }

    const mode_data_t &get_time_integrated_mode_data(const int l, const int m)
    {
        if (abs(m) > l || l < 2)
        {
            return m_zero_mode;
        }
        const int idx_mode = get_mode_index(l, m);
        assert(idx_mode >= 0);
        // will do nothing if data is already time integrated so cheap to call
        // after first time
        time_integrate_mode(idx_mode);
        return m_time_integrated_mode_data[idx_mode];
    }

    const std::vector<std::pair<int, int>> &get_available_modes() const
    {
        return m_available_modes;
    }

    // computes Eq. (8.9.69) in Alcubierre
    inline double momentum_a_coeff(const int l, const int m) const
    {
        return std::sqrt(static_cast<double>((l - m) * (l + m + 1))) /
               static_cast<double>((l * (l + 1)));
    }

    // computes Eq. (8.9.70) in Alcubierre
    inline double momentum_b_coeff(const int l, const int m) const
    {
        double numerator = (l - 2) * (l + 2) * (l + m) * (l + m - 1);
        double denominator = (2 * l - 1) * (2 * l + 1);
        return std::sqrt(numerator / denominator) / static_cast<double>(2 * l);
    }

    // computes Eq. (8.9.71) in Alcubierre
    inline double momentum_c_coeff(const int l, const int m) const
    {
        return static_cast<double>(2 * m) / static_cast<double>(l * (l + 1));
    }

    // computes Eq. (8.9.72) in Alcubierre
    inline double momentum_d_coeff(const int l, const int m) const
    {
        double numerator = (l - 2) * (l + 2) * (l - m) * (l + m);
        double denominator = (2 * l - 1) * (2 * l + 1);
        return std::sqrt(numerator / denominator) / static_cast<double>(l);
    }

    // check that we have the correct modes of Psi4 in order to compute the
    // desired term in Eqs (8.9.67-8.9.68) in Alcubierre
    bool check_momentum_pmode_computable(const int l, const int m)
    {
        bool missing_mode = (get_mode_index(l, m) < 0);
        missing_mode |= (get_mode_index(l, m + 1) < 0) && (abs(m + 1) <= l);
        missing_mode |=
            (get_mode_index(l - 1, m + 1) < 0) && (abs(m + 1) < l && l > 2);
        missing_mode |= (get_mode_index(l + 1, m + 1) < 0);
        missing_mode |= (get_mode_index(l - 1, m) < 0) && (abs(m) < l && l > 2);
        missing_mode |= (get_mode_index(l + 1, m) < 0);
        return !missing_mode;
    }

    // computes Eqs (8.9.67-8.9.68) in Alcubierre
    time_multidata_t compute_momentum_flux_pmode(const int l, const int m,
                                                 const double r_ex)
    {
        assert(l > 1 && abs(m) <= l);
        assert(check_momentum_pmode_computable(l, m));
        // uses the notation from Alcubierre but with _ instead of - and 1
        // instead of +1
        // any non-existent modes with |m|>l or l<2 are just set to zero
        // modes are already time integrated as well
        const mode_data_t &A_lm = get_time_integrated_mode_data(l, m);
        const mode_data_t &A_lm1 = get_time_integrated_mode_data(l, m + 1);
        const mode_data_t &A_l_1m1 =
            get_time_integrated_mode_data(l - 1, m + 1);
        const mode_data_t &A_l1m1 = get_time_integrated_mode_data(l + 1, m + 1);
        const mode_data_t &A_l_1m = get_time_integrated_mode_data(l - 1, m);
        const mode_data_t &A_l1m = get_time_integrated_mode_data(l + 1, m);

        const double a_lm = momentum_a_coeff(l, m);
        const double b_l_m = momentum_b_coeff(l, -m);
        const double b_l1m1 = momentum_b_coeff(l + 1, m + 1);
        const double c_lm = momentum_c_coeff(l, m);
        const double d_lm = momentum_d_coeff(l, m);
        const double d_l1m = momentum_d_coeff(l + 1, m + 1);

        constexpr int dim = 3;
        time_multidata_t out(dim * m_num_extraction_radii,
                             std::vector<double>(m_num_steps));

        for (int iradius = 0; iradius < m_num_extraction_radii; ++iradius)
        {
            // indices for real and imaginary parts in mode data
            int i_re = 2 * iradius;
            int i_im = 2 * iradius + 1;
            // indices for momentum flux direction in output
            int i_x = dim * iradius;
            int i_y = dim * iradius + 1;
            int i_z = dim * iradius + 2;
            for (int istep = 0; istep < m_num_steps; ++istep)
            {
                out[i_x][istep] =
                    r_ex * r_ex *
                    (A_lm[i_re][istep] * (a_lm * A_lm1[i_re][istep] +
                                          b_l_m * A_l_1m1[i_re][istep] -
                                          b_l1m1 * A_l1m1[i_re][istep]) +
                     A_lm[i_im][istep] * (a_lm * A_lm1[i_im][istep] +
                                          b_l_m * A_l_1m1[i_im][istep] -
                                          b_l1m1 * A_l1m1[i_im][istep])) /
                    (8.0 * M_PI);
                out[i_y][istep] =
                    r_ex * r_ex *
                    (A_lm[i_re][istep] * (-a_lm * A_lm1[i_im][istep] -
                                          b_l_m * A_l_1m1[i_im][istep] +
                                          b_l1m1 * A_l1m1[i_im][istep]) +
                     A_lm[i_im][istep] * (a_lm * A_lm1[i_re][istep] +
                                          b_l_m * A_l_1m1[i_re][istep] -
                                          b_l1m1 * A_l1m1[i_re][istep])) /
                    (8.0 * M_PI);
                out[i_z][istep] =
                    r_ex * r_ex *
                    (A_lm[i_re][istep] * (c_lm * A_lm[i_re][istep] +
                                          d_lm * A_l_1m[i_re][istep] +
                                          d_l1m * A_l1m[i_re][istep]) +
                     A_lm[i_im][istep] * (c_lm * A_lm[i_im][istep] +
                                          d_lm * A_l_1m[i_im][istep] +
                                          d_l1m * A_l1m[i_im][istep])) /
                    (16.0 * M_PI);
            }
        }
        return out;
    }
};

#endif /* WEYLMODEDATA_HPP_ */
