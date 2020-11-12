/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TIMEDATA_HPP
#define TIMEDATA_HPP

#include "IntegrationMethod.hpp"
#include "SmallDataIOReader.hpp"
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

// A class to read, store and manipulate the data extracted by SurfaceExtraction
class TimeData
{
  public:
    // type aliases
    using time_data_t = std::vector<double>;
    using time_multidata_t = std::vector<time_data_t>;

  protected:
    bool m_data_read;
    SmallDataIOReader m_file_reader;
    double m_dt;
    std::pair<double, double> m_time_limits;
    int m_num_datasets;
    int m_num_steps;
    time_multidata_t m_data;

  public:
    // Constructor
    TimeData();

    // resize a time_data_t object to the appropriate length
    void resize_time_data(time_data_t &a_time_data) const;

    // resize a time_multidata_t object to the appropriate size
    void resize_time_multidata(time_multidata_t &a_time_multidata) const;

    // Read data
    void read_data(const std::string &a_filename, int a_block = 0);

    // Read data from header
    std::vector<double> read_data_from_header(const int a_header_row_number,
                                              const int a_block = 0);

    // Getters
    const time_multidata_t &get_data() const;
    double get_dt() const;
    const std::pair<double, double> &get_time_limits() const;

    // Integrate a dataset in time from step a_min_step to a_max_step inclusive.
    double integrate_time(const time_data_t &in_data, int a_min_step,
                          int a_max_step) const;

    // Integrate all datasets in time from step a_min_step to a_max_step
    // inclusive
    std::vector<double> integrate_time(const time_multidata_t &in_data,
                                       int a_min_step, int a_max_step) const;

    // For each step s, the partial integral of a datasets is  computed from the
    // minimum available step to s. This can be expensive so OpenMP is used if
    // available.
    void integrate_all_time(time_data_t &out, const time_data_t &in_data) const;

    // For each step s, the partial integral of a datasets is  computed from the
    // minimum available step to s. This can be expensive so OpenMP is used if
    // available.
    void integrate_all_time(time_multidata_t &out,
                            const time_multidata_t &in_data) const;

    // Take the norm of a_dim dimensional data
    // assumes that columns are of the form,
    // data1[0], ..., data1[a_dim - 1], data2[0], ..., data2[a_dim - 1], ...
    time_multidata_t norm(const time_multidata_t &in_data, const int a_dim = 3);

    // Adds a value to every row in each column. The value to add to column i is
    // given by the ith element of a_values_to_add
    void add_to_columns(time_multidata_t &data,
                        const std::vector<double> &a_values_to_add);

    // Clear all data and file information
    void clear();

    // need Simpson's 3/8 rule IntegrationMethod to do integration with odd
    // number of points
    static const IntegrationMethod simpson_odd;
};

#endif /* TIMEDATA_HPP */
