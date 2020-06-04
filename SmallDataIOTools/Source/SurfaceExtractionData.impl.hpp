/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SURFACEEXTRACTIONDATA_HPP)
#error "This file should only be included in SurfaceExtractionData.hpp"
#endif

#ifndef SURFACEEXTRACTIONDATA_IMPL_HPP
#define SURFACEEXTRACTIONDATA_IMPL_HPP

#include <cassert>

template <typename SurfaceGeometry>
void SurfaceExtractionData::integrate_surface(
    multisurface_value_t &out, const multisurface_multidata_t &in_data,
    const integrand_t &a_integrand, const SurfaceGeometry &a_geom,
    bool a_use_area_element, const IntegrationMethod &a_method_u,
    const IntegrationMethod &a_method_v) const
{
    double du = a_geom.du(m_num_points_u);
    double dv = a_geom.dv(m_num_points_v);
    out.clear();
    out.resize(m_num_surfaces);
    std::fill(out.begin(), out.end(), 0.0);
    for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
    {
        const auto surface_it = m_surfaces.begin();
        std::advance(surface_it, isurface);
        double surface_param_value = surface_it->second;
        for (int iu = 0; iu < m_num_points_u; ++iu)
        {
            double u = a_geom.u(iu, m_num_points_u);
            double inner_integral = 0.0;
            for (int iv = 0; iv < m_num_points_v; ++iv)
            {
                double v = a_geom.v(iv, m_num_points_v);
                int num_datasets = in_data[isurface].size();
                std::vector<double> data_here(num_datasets);
                for (int idataset = 0; idataset < num_datasets; ++idataset)
                {
                    data_here[idataset] = in_data[isurface][idataset][iu][iv];
                }
                double area_element =
                    a_use_area_element
                        ? a_geom.area_element(surface_param_value, u, v)
                        : 1.0;
                double integrand_with_area_element =
                    a_integrand(data_here, surface_param_value, u, v) *
                    area_element;
                double weight =
                    a_method_v.weight(iv, m_num_points_v, a_geom.is_periodic());
                inner_integral += weight * dv * integrand_with_area_element;
            }
            double weight =
                a_method_u.weight(iu, m_num_points_u, a_geom.is_periodic());
            out[isurface] += weight * du * inner_integral;
        }
    }
}

template <typename SurfaceGeometry>
void SurfaceExtractionData::integrate_surface(
    time_multisurface_value_t &out, const extracted_data_t &in_data,
    const integrand_t &a_integrand, const SurfaceGeometry &a_geom,
    int a_min_step, int a_max_step, bool a_use_area_element,
    const IntegrationMethod &a_method_u,
    const IntegrationMethod &a_method_v) const
{
    assert(a_min_step <= a_max_step);
    const int num_intervals = a_max_step - a_min_step;
    const int num_steps = num_intervals + 1;
    out.clear();
    out.resize(num_steps);
#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(out, in_data, a_integrand, a_geom, a_method_u, a_method_v)
#endif
    for (int istep = a_min_step; istep <= a_max_step; ++istep)
    {
        integrate_surface(out[istep], in_data[istep], a_integrand, a_geom,
                          a_use_area_element, a_method_u, a_method_v);
    }
}

#endif /* SURFACEEXTRACTIONDATA_IMPL_HPP */
