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
#include <filesystem>
#include <iostream>

//! This class reads, stores and manipulates modes of the Weyl scalar, Psi4
class WeylModeData
{
  public:
    using mode_data_t = TimeData::time_multidata_t;

  protected:
    std::string m_mode_file_parent_path;
    std::string m_mode_file_pattern_begin;
    std::string m_mode_file_pattern_middle;
    std::string m_mode_file_pattern_end;
    std::vector<std::pair<int, int>> m_available_modes;
    TimeData m_time_data;
    std::vector<mode_data_t> m_mode_data;
    bool m_data_read = false;
    std::vector<mode_data_t> m_time_integrated_mode_data;

    void find_available_modes(const std::string &a_mode_file_pattern)
    {
        namespace fs = std::filesystem;

        fs::path pattern_path(a_mode_file_pattern);
        // get just the filename without the path
        std::string pattern = pattern_path.filename();

        // assume l_pos < m_pos
        const auto l_pos_pattern = pattern.find("{l}");
        m_mode_file_pattern_begin = pattern.substr(0, l_pos_pattern);
        std::cout << "pattern_begin: " << m_mode_file_pattern_begin << "\n";
        const auto m_pos_pattern = pattern.find("{m}");
        m_mode_file_pattern_middle = pattern.substr(
            l_pos_pattern + 3, m_pos_pattern - l_pos_pattern - 3);
        std::cout << "pattern_middle: " << m_mode_file_pattern_middle << "\n ";
        m_mode_file_pattern_end = pattern.substr(m_pos_pattern + 3);
        std::cout << "pattern_end: " << m_mode_file_pattern_end << "\n";

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
                        std::cout << l << " " << m << "\n";
                        m_available_modes.emplace_back(l, m);
                    }
                }
            }
        }
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
    int get_mode_index(const int l, const int m)
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
        m_mode_data.resize(m_available_modes.size());

        for (int imode = 0; imode < m_available_modes.size(); ++imode)
        {
            m_time_data.clear();
            const auto &mode = m_available_modes[imode];
            const std::string filepath =
                get_mode_filepath(mode.first, mode.second);
            m_time_data.read_data(filepath);
            m_mode_data[imode] = std::move(m_time_data.get_data());
        }
        m_data_read = true;
    }

    void compute_time_integrated_mode_data()
    {
        assert(m_data_read);
        m_time_integrated_mode_data.resize(m_available_modes.size());

        for (int imode = 0; imode < m_available_modes.size(); ++imode)
        {
            m_time_data.integrate_all_time(m_time_integrated_mode_data[imode],
                                           m_mode_data[imode]);
        }
    }

    // getters
    const std::vector<mode_data_t> &get_mode_data() const
    {
        assert(m_data_read);
        return m_mode_data;
    }
};

#endif /* WEYLMODEDATA_HPP_ */
