/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SurfaceExtractionData.hpp"
#include <cassert>
#include <iterator>
#include <numeric>

// Simpons's 3/8 rule for integration with odd number of points
const IntegrationMethod SurfaceExtractionData::simpson_odd({0.375, 1.125,
                                                            1.25});

// Constructor
SurfaceExtractionData::SurfaceExtractionData(const std::string &a_file_prefix)
    : m_file_prefix(a_file_prefix), m_structure_determined(false),
      m_data_read(false)
{
}

// Determines the structure of data i.e.
// how many datasets, number of points in u, number of points in v,
// number of timesteps, min and max times, dt
void SurfaceExtractionData::determine_data_structure(
    std::set<int> a_surface_indices, int a_min_grchombo_step,
    int a_max_grchombo_step)
{

    namespace fs = std::filesystem;

    // first iterate through all files in directory and find the ones which
    // start with prefix
    fs::path prefix_path(m_file_prefix);
    std::string prefix = prefix_path.stem();
    fs::path parent_path = prefix_path.parent_path();

    for (auto &file : fs::directory_iterator(parent_path))
    {
        const std::string &filename = file.path().filename();
        if (filename.compare(0, prefix.length(), prefix) == 0)
        {
            // remove extension and convert to string
            const std::string &stem_str = file.path().stem();

            // get step number from filename and convert to int
            const std::string &step_str =
                stem_str.substr(prefix.length(), stem_str.length());
            int step = std::stoi(step_str);

            // store in files map
            if (step >= a_min_grchombo_step && step <= a_max_grchombo_step)
                m_files.insert(std::make_pair(step, file.path()));
        }
    }

    m_num_steps = m_files.size();
    // check that the number of steps divides the step number difference
    // between the first and last step - a necessary but not sufficient
    // condition for data equally spaced in time
    assert((m_files.end()->first - m_files.begin()->first) % m_num_steps == 0);

    // next open the first file and get the file structure
    m_file_reader.open(m_files.begin()->second);
    m_file_reader.determine_file_structure();
    m_file_structure = m_file_reader.get_file_structure();
    // the number of datasets will be the same in all blocks so use the first
    // one
    const int first_block = 0;
    m_num_datasets = m_file_structure.num_data_columns[first_block];

    // get the requested surfaces
    int total_num_surfaces = m_file_structure.num_blocks;
    std::set<int> requested_surface_indices = a_surface_indices;

    // if a_surface_indices is empty, request all surfaces
    if (requested_surface_indices.empty())
    {
        for (int isurface = 0; isurface < total_num_surfaces; ++isurface)
        {
            requested_surface_indices.insert(isurface);
        }
    }

    const int first_header = 0;
    for (int isurface : requested_surface_indices)
    {
        assert(0 <= isurface && isurface < total_num_surfaces);
        // the surface param value is second (time is first)
        const int surface_param_header_index = 1;
        double surface_param_value = m_file_reader.get_data_from_header(
            first_header, isurface)[surface_param_header_index];
        m_surfaces.insert(std::make_pair(isurface, surface_param_value));
    }
    m_num_surfaces = m_surfaces.size();

    const int first_column = 0;

    // the u_coords will have blocks of length num_points_v where all the
    // elements are the same.
    std::vector<double> u_coords =
        m_file_reader.get_column(first_column, m_surfaces.begin()->first);
    int iv = 0;
    while (u_coords[++iv] == u_coords[0])
        ;
    m_num_points_v = iv;

    // the number of u points must then be the number of rows divided by the
    // number of v points
    m_num_points_u =
        m_file_structure.num_data_rows[first_block] / m_num_points_v;

    // get time information
    // lower time limit given by first value in first header in first file
    const int time_header_index = 0;
    m_time_limits.first = m_file_reader.get_data_from_header(
        first_header, m_surfaces.begin()->first)[time_header_index];
    m_file_reader.close();

    // upper time limit from last file
    auto filesit = m_files.end();
    // end iterator points past last element so decrement by 1
    std::advance(filesit, -1);
    m_file_reader.open(filesit->second);
    // file structure assumed to be the same for all files
    m_file_reader.set_file_structure(m_file_structure);
    m_time_limits.second =
        m_file_reader.get_data_from_header(first_header)[time_header_index];
    m_file_reader.close();

    // calculate dt assuming equally spaced in time data
    m_dt = (m_time_limits.second - m_time_limits.first) /
           static_cast<double>(m_num_steps - 1);

    m_structure_determined = true;
}

// resize a surface_data_t object to the appropriate size
void SurfaceExtractionData::resize_surface_data(
    surface_data_t &a_surface_data) const
{
    assert(m_structure_determined);
    a_surface_data.resize(m_num_points_u);
    for (auto &fixed_u_data : a_surface_data)
    {
        fixed_u_data.resize(m_num_points_v);
    }
}

// resize a time_surface_data_t object to the appropriate size
void SurfaceExtractionData::resize_surface_multidata(
    surface_multidata_t &a_surface_multidata) const
{
    assert(m_structure_determined);
    a_surface_multidata.resize(m_num_datasets);
    for (auto &fixed_dataset_data : a_surface_multidata)
    {
        resize_surface_data(fixed_dataset_data);
    }
}

// resize a time_multisurface_data_t object to the appopriate size
void SurfaceExtractionData::resize_multisurface_multidata(
    multisurface_multidata_t &a_multisurface_multidata) const
{
    assert(m_structure_determined);
    a_multisurface_multidata.resize(m_num_surfaces);
    for (auto &fixed_surface_data : a_multisurface_multidata)
    {
        resize_surface_multidata(fixed_surface_data);
    }
}

// resize an extracted_data_t object to the appropriate size
void SurfaceExtractionData::resize_extracted_data(
    extracted_data_t &a_extracted_data_t) const
{
    assert(m_structure_determined);
    a_extracted_data_t.resize(m_num_steps);
    for (auto &fixed_step_data : a_extracted_data_t)
    {
        resize_multisurface_multidata(fixed_step_data);
    }
}

void SurfaceExtractionData::read_data()
{
    resize_extracted_data(m_data);

    for (int istep = 0; istep < m_num_steps; ++istep)
    {
        auto filesit = m_files.begin();
        std::advance(filesit, istep);
        m_file_reader.open(filesit->second);
        std::cout << "SurfaceExtractionData::read_data: Reading data from step "
                  << istep << "/" << m_num_steps << "\r";
        m_file_reader.set_file_structure(m_file_structure);
        for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
        {
            auto surface_it = m_surfaces.begin();
            std::advance(surface_it, isurface);
            auto read_data =
                m_file_reader.get_all_data_columns(surface_it->first);
            for (int idataset = 0; idataset < m_num_datasets; ++idataset)
            {
                for (int iu = 0; iu < m_num_points_u; ++iu)
                {
                    for (int iv = 0; iv < m_num_points_v; ++iv)
                    {
                        m_data[istep][isurface][idataset][iu][iv] =
                            read_data[idataset][index(iu, iv)];
                    }
                }
            }
        }
        m_file_reader.close();
    }
    std::cout << "Finished reading data                                    \n";
    long unsigned int memory_count =
        sizeof(m_data[0][0][0][0][0]) * m_num_datasets * m_num_surfaces *
        m_num_steps * m_num_points_u * m_num_points_v;
    long unsigned int memory_count_KiB = memory_count / 1024;
    long unsigned int memory_count_MiB = memory_count_KiB / 1024;
    std::cout << "Memory used = " << memory_count_MiB << " MiB\n";
    m_data_read = true;
}

const SurfaceExtractionData::extracted_data_t &
SurfaceExtractionData::get_data() const
{
    assert(m_data_read);
    return m_data;
}

void SurfaceExtractionData::integrate_time(multisurface_multidata_t &out,
                                           const extracted_data_t &in_data,
                                           int a_min_step, int a_max_step) const
{
    assert(m_structure_determined);
    assert(a_min_step <= a_max_step);
    out.clear();
    resize_multisurface_multidata(out);
    const int num_intervals = a_max_step - a_min_step;
    const int num_steps = num_intervals + 1;
    constexpr bool is_periodic = false;
    for (auto &fixed_surface_data : out)
    {
        for (auto &fixed_dataset_data : fixed_surface_data)
        {
            for (auto &fixed_u_data : fixed_dataset_data)
            {
                for (auto &data : fixed_u_data)
                {
                    data = 0.0;
                }
            }
        }
    }

    if (num_intervals == 0)
    {
        return;
    }
    else if (num_steps < 3) // Simpson's rule not possible -> use trapezium rule
    {
        for (int istep = a_min_step; istep <= a_max_step; ++istep)
        {
            for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
            {
                for (int idataset = 0; idataset < m_num_datasets; ++idataset)
                {
                    for (int iu = 0; iu < m_num_points_u; ++iu)
                    {
                        for (int iv = 0; iv < m_num_points_v; ++iv)
                        {
                            out[isurface][idataset][iu][iv] +=
                                m_dt *
                                IntegrationMethod::trapezium.weight(
                                    istep - a_min_step, num_steps,
                                    is_periodic) *
                                in_data[istep][isurface][idataset][iu][iv];
                        }
                    }
                }
            }
        }
        return;
    }
    else // Use Simpson's rule
    {
        if (num_intervals % 2 == 0) // even case
        {
            for (int istep = a_min_step; istep <= a_max_step; ++istep)
            {
                for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
                {
                    for (int idataset = 0; idataset < m_num_datasets;
                         ++idataset)
                    {
                        for (int iu = 0; iu < m_num_points_u; ++iu)
                        {
                            for (int iv = 0; iv < m_num_points_v; ++iv)
                            {
                                out[isurface][idataset][iu][iv] +=
                                    m_dt *
                                    IntegrationMethod::simpson.weight(
                                        istep - a_min_step, num_steps,
                                        is_periodic) *
                                    in_data[istep][isurface][idataset][iu][iv];
                            }
                        }
                    }
                }
            }
            return;
        }
        else // odd case
        {
            // do first 3 steps using odd simpson formula
            const int num_odd_intervals = 3;
            const int num_odd_steps = num_odd_intervals + 1;
            const int odd_max_step = a_min_step + num_odd_intervals;
            for (int istep = a_min_step; istep <= odd_max_step; ++istep)
            {
                for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
                {
                    for (int idataset = 0; idataset < m_num_datasets;
                         ++idataset)
                    {
                        for (int iu = 0; iu < m_num_points_u; ++iu)
                        {
                            for (int iv = 0; iv < m_num_points_v; ++iv)
                            {
                                out[isurface][idataset][iu][iv] +=
                                    m_dt *
                                    SurfaceExtractionData::simpson_odd.weight(
                                        istep - a_min_step, num_odd_steps,
                                        is_periodic) *
                                    in_data[istep][isurface][idataset][iu][iv];
                            }
                        }
                    }
                }
            }
            // add weight again for last odd step (composition of Newton-Cotes
            // formulae)
            const int even_min_step = odd_max_step;
            const int num_even_steps = a_max_step - even_min_step + 1;
            for (int istep = even_min_step; istep <= a_max_step; ++istep)
            {
                for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
                {
                    for (int idataset = 0; idataset < m_num_datasets;
                         ++idataset)
                    {
                        for (int iu = 0; iu < m_num_points_u; ++iu)
                        {
                            for (int iv = 0; iv < m_num_points_v; ++iv)
                            {
                                out[isurface][idataset][iu][iv] +=
                                    m_dt *
                                    IntegrationMethod::simpson.weight(
                                        istep - even_min_step, num_even_steps,
                                        is_periodic) *
                                    in_data[istep][isurface][idataset][iu][iv];
                            }
                        }
                    }
                }
            }
            return;
        }
    }
}

void SurfaceExtractionData::integrate_all_time(
    extracted_data_t &out, const extracted_data_t &in_data) const
{
    assert(m_structure_determined);
    out.resize(m_num_steps);

    // Use dynamic schedule for OpenMP loop as load is proportional to index
    // (i.e. unbalanced)
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(out, in_data, std::cout)         \
    schedule(dynamic)
#endif /* _OPENMP */
    for (int istep = 0; istep < m_num_steps; ++istep)
    {
        // print every 100 iterations
        if (istep % 100 == 0)
        {
#ifdef _OPENMP
#pragma omp critical(print)
#endif /* _OPENMP */
            {
                std::cout << "SurfaceExtractionData::integrate_all_time: Step "
                          << istep << "/" << m_num_steps - 1 << "\r";
                std::cout << std::flush;
            }
        }
        integrate_time(out[istep], in_data, 0, istep);
    }
}

void SurfaceExtractionData::integrate_time(
    multisurface_value_t &out, const time_multisurface_value_t &in_data,
    int a_min_step, int a_max_step) const
{
    assert(m_structure_determined);
    assert(a_min_step <= a_max_step);
    out.clear();
    out.resize(m_num_surfaces);
    const int num_intervals = a_max_step - a_min_step;
    const int num_steps = num_intervals + 1;
    constexpr bool is_periodic = false;

    std::fill(out.begin(), out.end(), 0.0);

    if (num_intervals == 0)
    {
        return;
    }
    else if (num_steps < 3) // Simpson's rule not possible -> use trapezium rule
    {
        for (int istep = a_min_step; istep <= a_max_step; ++istep)
        {
            for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
            {
                out[isurface] +=
                    m_dt *
                    IntegrationMethod::trapezium.weight(
                        istep - a_min_step, num_steps, is_periodic) *
                    in_data[istep][isurface];
            }
        }
        return;
    }
    else // Use Simpson's rule
    {
        if (num_intervals % 2 == 0) // even case
        {
            for (int istep = a_min_step; istep <= a_max_step; ++istep)
            {
                for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
                {
                    out[isurface] +=
                        m_dt *
                        IntegrationMethod::simpson.weight(
                            istep - a_min_step, num_steps, is_periodic) *
                        in_data[istep][isurface];
                }
            }
            return;
        }
        else // odd case
        {
            // do first 3 steps using odd simpson formula
            const int num_odd_intervals = 3;
            const int num_odd_steps = num_odd_intervals + 1;
            const int odd_max_step = a_min_step + num_odd_intervals;
            for (int istep = a_min_step; istep <= odd_max_step; ++istep)
            {
                for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
                {
                    out[isurface] +=
                        m_dt *
                        SurfaceExtractionData::simpson_odd.weight(
                            istep - a_min_step, num_odd_steps, is_periodic) *
                        in_data[istep][isurface];
                }
            }
            // add weight again for last odd step (composition of Newton-Cotes
            // formulae)
            const int even_min_step = odd_max_step;
            const int num_even_steps = a_max_step - even_min_step + 1;
            for (int istep = even_min_step; istep <= a_max_step; ++istep)
            {
                for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
                {
                    out[isurface] += m_dt *
                                     IntegrationMethod::simpson.weight(
                                         istep - even_min_step, num_even_steps,
                                         is_periodic) *
                                     in_data[istep][isurface];
                }
            }
            return;
        }
    }
}

SurfaceExtractionData::surface_multidata_t SurfaceExtractionData::combine(
    const SurfaceExtractionData::surface_multidata_t &a,
    const SurfaceExtractionData::surface_multidata_t &b)
{
    surface_multidata_t out;
    const int num_datasets_a = a.size();
    const int num_datasets_b = b.size();
    const int num_datasets_out = num_datasets_a + num_datasets_b;
    out.resize(num_datasets_out);
    for (int idataset = 0; idataset < num_datasets_a; ++idataset)
    {
        out[idataset] = a[idataset];
    }
    for (int idataset = 0; idataset < num_datasets_b; ++idataset)
    {
        int idataset_out = idataset + num_datasets_a;
        out[idataset_out] = b[idataset];
    }
    return out;
}

SurfaceExtractionData::multisurface_multidata_t SurfaceExtractionData::combine(
    const SurfaceExtractionData::multisurface_multidata_t &a,
    const SurfaceExtractionData::multisurface_multidata_t &b)
{
    multisurface_multidata_t out;
    out.resize(m_num_surfaces);
    for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
    {
        out[isurface] = combine(a[isurface], b[isurface]);
    }
    return out;
}

SurfaceExtractionData::extracted_data_t
SurfaceExtractionData::combine(const SurfaceExtractionData::extracted_data_t &a,
                               const SurfaceExtractionData::extracted_data_t &b)
{
    extracted_data_t out;
    out.resize(m_num_steps);
    for (int istep = 0; istep < m_num_steps; ++istep)
    {
        out[istep] = combine(a[istep], b[istep]);
    }
    return out;
}

void SurfaceExtractionData::clear()
{
    m_structure_determined = false;
    m_data_read = false;
    m_file_reader.close();
    m_file_structure.clear();
    m_file_prefix.clear();
    m_files.clear();
    m_dt = 0.;
    m_time_limits.first = 0.;
    m_time_limits.second = 0.;
    m_num_datasets = 0;
    m_num_surfaces = 0;
    m_surfaces.clear();
    m_num_steps = 0;
    m_num_points_u = 0;
    m_num_points_v = 0;
    m_data.clear();
}
