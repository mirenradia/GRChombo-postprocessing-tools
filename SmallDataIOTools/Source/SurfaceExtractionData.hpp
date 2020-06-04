/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SURFACEEXTRACTIONDATA_HPP
#define SURFACEEXTRACTIONDATA_HPP

#include "CH_assert.H"

#include "IntegrationMethod.hpp"
#include "SmallDataIOReader.hpp"
#include <algorithm>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

// A class to read, store and manipulate the data extracted by SurfaceExtraction
class SurfaceExtractionData
{
  public:
    // main type alises used for m_data
    // indices: timestep surface_param_value dataset u v
    using surface_data_t = std::vector<std::vector<double>>;
    using surface_multidata_t = std::vector<surface_data_t>;
    using multisurface_multidata_t = std::vector<surface_multidata_t>;
    using extracted_data_t = std::vector<multisurface_multidata_t>;

    // other type aliases which arise from integration
    using multisurface_value_t = std::vector<double>;
    using time_multisurface_value_t = std::vector<multisurface_value_t>;

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
    std::map<int, double>
        m_surfaces; // pairs of {surface_index_in_file, param_value}
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
    // If a_surface_indices is an empty set, all surfaces are extracted
    void determine_data_structure(
        std::set<int> a_surface_indices, int a_min_grchombo_step = 0,
        int a_max_grchombo_step = std::numeric_limits<int>::max());

    // resize a surface_data_t object to the appropriate size
    void resize_surface_data(surface_data_t &a_surface_data) const;

    // resize a time_surface_data_t object to the appropriate size
    void
    resize_surface_multidata(surface_multidata_t &a_surface_multidata) const;

    // resize a time_multisurface_data_t object to the appopriate size
    void resize_multisurface_multidata(
        multisurface_multidata_t &a_multisurface_multidata) const;

    // resize an extracted_data_t object to the appropriate size
    void resize_extracted_data(extracted_data_t &a_extracted_data_t) const;

    // Read data
    void read_data();

    // Get Data
    const extracted_data_t &get_data() const;

    // Integrate all datasets on all surfaces in time from step
    // a_min_step to a_max_step inclusive. Note these steps corresponds to the
    // first index in_data
    void integrate_time(multisurface_multidata_t &out,
                        const extracted_data_t &in_data, int a_min_step,
                        int a_max_step) const;

    // For each step s, the partial integral of all datasets on all surfaces is
    // computed from the minimum available step to s.
    // This can be expensive so OpenMP is used if available.
    void integrate_all_time(extracted_data_t &out,
                            const extracted_data_t &in_data) const;

    // Integrate a value on all surfaces in time from step a_min_step to
    // a_max_step inclusive
    void integrate_time(multisurface_value_t &out,
                        const time_multisurface_value_t &in_data,
                        int a_min_step, int a_max_step) const;

    // this is the type used for integrands on the surface, the vector<double>
    // is a vector of all the datasets at the point on the surface.
    using integrand_t =
        std::function<double(std::vector<double> &, double, double, double)>;

    // Integrate an integrand_t over all surfaces
    template <typename SurfaceGeometry>
    void integrate_surface(
        multisurface_value_t &out, const multisurface_multidata_t &in_data,
        const integrand_t &a_integrand, const SurfaceGeometry &a_geom,
        bool a_use_area_element = true,
        const IntegrationMethod &a_method_u = IntegrationMethod::trapezium,
        const IntegrationMethod &a_method_v =
            IntegrationMethod::trapezium) const;

    // Integrate an integrand_t over all surfaces for each step from a_min_step
    // to a_max_step inclusive
    template <typename SurfaceGeometry>
    void integrate_surface(
        time_multisurface_value_t &out, const extracted_data_t &in_data,
        const integrand_t &a_integrand, const SurfaceGeometry &a_geom,
        int a_min_step, int a_max_step, bool a_use_area_element = true,
        const IntegrationMethod &a_method_u = IntegrationMethod::trapezium,
        const IntegrationMethod &a_method_v =
            IntegrationMethod::trapezium) const;

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

#include "SurfaceExtractionData.impl.hpp"

#endif /* SURFACEEXTRACTIONDATA_HPP */
