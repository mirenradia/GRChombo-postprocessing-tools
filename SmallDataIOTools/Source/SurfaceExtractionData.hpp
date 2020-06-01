/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SURFACEEXTRACTIONDATA_HPP
#define SURFACEEXTRACTIONDATA_HPP

#include "CH_assert.H"

#include "IntegrationMethod.hpp"
#include "SmallDataIOReader.hpp"
#include <filesystem>
#include <limits>
#include <map>
#include <string>
#include <vector>

// A class to read and store the data extracted by SurfaceExtraction
class SurfaceExtractionData
{
  public:
    // indices: dataset surface_param_value timestep u v
    using surface_data_t = std::vector<std::vector<double>>;
    using time_surface_data_t = std::vector<surface_data_t>;
    using time_multisurface_data_t = std::vector<time_surface_data_t>;
    using extracted_data_t = std::vector<time_multisurface_data_t>;

  protected:
    bool m_structure_determined;
    bool m_data_read;
    SmallDataIOReader m_file_reader;
    SmallDataIOReader::file_structure_t m_file_structure;
    std::string m_file_prefix;
    std::map<int, std::filesystem::path> m_files;
    double m_dt;
    std::pair<double, double> m_time_limits;
    int m_num_datasets;
    int m_num_surfaces;
    std::vector<double> m_surface_param_values;
    int m_num_steps;
    int m_num_points_u;
    int m_num_points_v;
    extracted_data_t m_data;

  public:
    // Constructor
    SurfaceExtractionData(const std::string &a_file_prefix);

    // Determines the structure of data i.e.
    // how many datasets, number of points in u, number of points in v,
    // number of timesteps, min and max times, dt
    void determine_data_structure(
        int a_min_grchombo_step = 0,
        int a_max_grchombo_step = std::numeric_limits<int>::max());

    // resize a surface_data_t object to the appropriate size
    void resize_surface_data(surface_data_t &a_surface_data) const;

    // resize a time_surface_data_t object to the appropriate size
    void
    resize_time_surface_data(time_surface_data_t &a_time_surface_data) const;

    // resize a time_multisurface_data_t object to the appopriate size
    void resize_time_multisurface_data(
        time_multisurface_data_t &a_time_multisurface_data) const;

    // resize an extracted_data_t object to the appropriate size
    void resize_extracted_data(extracted_data_t &a_extracted_data_t) const;

    // Read data
    void read_data();

    // Get Data
    const extracted_data_t &get_data() const;

    // Integrate a single dataset on a single surface in time from step
    // a_min_step to a_max_step. Note these steps corresponds to the first index
    // in_data
    void integrate_time(surface_data_t &out, const time_surface_data_t &in_data,
                        int a_min_step, int a_max_step) const;

    // Clear all data and structure information
    void clear();

    // need Simpson's 3/8 rule IntegrationMethod to do integration with odd
    // number of points
    static const IntegrationMethod simpson_odd;

  protected:
    // returns the flattened index corresponding to the row in a block of the
    // SurfaceExtraction file
    inline int index(int a_iu, int a_iv) const
    {
        return a_iu * m_num_points_v + a_iv;
    }
};

#endif /* SURFACEEXTRACTIONDATA_HPP */
