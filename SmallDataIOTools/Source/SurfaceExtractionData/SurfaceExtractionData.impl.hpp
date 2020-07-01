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
#include <iomanip>
#include <iostream>

template <typename SurfaceGeometry>
SurfaceExtractionData::surface_multidata_t
SurfaceExtractionData::compute_surface_derivatives(
    const surface_data_t &in_data, const SurfaceGeometry &a_geom) const
{
    assert(m_num_points_u > 2 && m_num_points_v > 2);
    surface_multidata_t out;
    constexpr int num_surface_dimensions = 2;
    out.resize(num_surface_dimensions);
    for (auto &deriv : out)
    {
        resize_surface_data(deriv);
    }
    constexpr int idu = 0;
    constexpr int idv = 1;
    const bool u_periodic = a_geom.is_u_periodic();
    const bool v_periodic = a_geom.is_v_periodic();
    const double duinv = 1.0 / a_geom.du(m_num_points_u);
    const double dvinv = 1.0 / a_geom.dv(m_num_points_v);
    for (int iu = 0; iu < m_num_points_u; ++iu)
    {
        for (int iv = 0; iv < m_num_points_v; ++iv)
        {
            // first calculate derviatives wrt u
            if (iu == 0)
            {
                if (u_periodic)
                {
                    // use centred stencil in periodic case
                    out[idu][iu][iv] =
                        0.5 * duinv *
                        (in_data[iu + 1][iv] - in_data[m_num_points_u - 1][iv]);
                }
                else
                {
                    // use one-sided stencil otherwise
                    out[idu][iu][iv] =
                        0.5 * duinv *
                        (-3.0 * in_data[iu][iv] + 4.0 * in_data[iu + 1][iv] -
                         in_data[iu + 2][iv]);
                }
            }
            else if (iu == m_num_points_u - 1)
            {
                if (u_periodic)
                {
                    // use centred stencil in periodic case
                    out[idu][iu][iv] =
                        0.5 * duinv * (in_data[0][iv] - in_data[iu - 1][iv]);
                }
                else
                {
                    // use one-sided stencil otherwise
                    out[idu][iu][iv] =
                        0.5 * duinv *
                        (in_data[iu - 2][iv] - 4.0 * in_data[iu - 1][iv] +
                         3.0 * in_data[iu][iv]);
                }
            }
            else
            {
                out[idu][iu][iv] =
                    0.5 * duinv * (in_data[iu + 1][iv] - in_data[iu - 1][iv]);
            }

            // now do derivatives wrt v
            if (iv == 0)
            {
                if (v_periodic)
                {
                    // use centered stencil in periodic case
                    out[idv][iu][iv] =
                        0.5 * dvinv *
                        (in_data[iu][iv + 1] - in_data[iu][m_num_points_v - 1]);
                }
                else
                {
                    // use one-sided stencil otherwise
                    out[idv][iu][iv] =
                        0.5 * dvinv *
                        (-3.0 * in_data[iu][iv] + 4.0 * in_data[iu][iv + 1] -
                         in_data[iu][iv + 2]);
                }
            }
            else if (iv == m_num_points_v - 1)
            {
                if (v_periodic)
                {
                    // use centred stencil in periodic case
                    out[idv][iu][iv] =
                        0.5 * dvinv * (in_data[iu][0] - in_data[iu][iv - 1]);
                }
                else
                {
                    // use one-sided stencil otherwise
                    out[idv][iu][iv] =
                        0.5 * dvinv *
                        (in_data[iu][iv - 2] - 4.0 * in_data[iu][iv - 1] +
                         3.0 * in_data[iu][iv]);
                }
            }
            else
            {
                out[idv][iu][iv] =
                    0.5 * dvinv * (in_data[iu][iv + 1] - in_data[iu][iv - 1]);
            }
        }
    }
    return out;
}

template <typename SurfaceGeometry>
SurfaceExtractionData::surface_multidata_t
SurfaceExtractionData::compute_surface_derivatives(
    const surface_multidata_t &in_data, const SurfaceGeometry &a_geom) const
{
    const int num_datasets = in_data.size();
    constexpr int num_surface_dimensions = 2;
    surface_multidata_t out(num_datasets * num_surface_dimensions);
    for (int idataset = 0; idataset < num_datasets; ++idataset)
    {
        surface_multidata_t derivative_data =
            compute_surface_derivatives(in_data[idataset], a_geom);
        for (int idir = 0; idir < num_surface_dimensions; ++idir)
        {
            out[idir + num_surface_dimensions * idataset] =
                std::move(derivative_data[idir]);
        }
    }
    return out;
}

template <typename SurfaceGeometry>
SurfaceExtractionData::multisurface_multidata_t
SurfaceExtractionData::compute_surface_derivatives(
    const multisurface_multidata_t &in_data,
    const SurfaceGeometry &a_geom) const
{
    multisurface_multidata_t out(m_num_surfaces);
    for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
    {
        out[isurface] = compute_surface_derivatives(in_data[isurface], a_geom);
    }
    return out;
}

template <typename SurfaceGeometry>
SurfaceExtractionData::extracted_data_t
SurfaceExtractionData::compute_surface_derivatives(
    const extracted_data_t &in_data, const SurfaceGeometry &a_geom) const
{
    extracted_data_t out(m_num_steps);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(in_data, a_geom, out, std::cout)
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
                std::cout << "SurfaceExtractionData::compute_surface_"
                             "derivatives: Step "
                          << std::setw(6) << istep << "/" << m_num_steps - 1
                          << "\r";
                std::cout << std::flush;
            }
        }
        out[istep] = compute_surface_derivatives(in_data[istep], a_geom);
    }
    std::cout << "\n";
    return out;
}

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
        auto surface_it = m_surfaces.begin();
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
                double weight = a_method_v.weight(iv, m_num_points_v,
                                                  a_geom.is_v_periodic());
                inner_integral += weight * dv * integrand_with_area_element;
            }
            double weight =
                a_method_u.weight(iu, m_num_points_u, a_geom.is_u_periodic());
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
    shared(out, in_data, a_integrand, a_geom, a_min_step, a_max_step,          \
           num_steps, a_use_area_element, a_method_u, a_method_v, std::cout)
#endif
    for (int istep = a_min_step; istep <= a_max_step; ++istep)
    {
        // print every 100 iterations
        if (istep % 100 == 0)
        {
#ifdef _OPENMP
#pragma omp critical(print)
#endif /* _OPENMP */
            {
                std::cout << "SurfaceExtractionData::integrate_surface: Step "
                          << std::setw(6) << istep << "/" << num_steps - 1
                          << "\r";
                std::cout << std::flush;
            }
        }
        integrate_surface(out[istep], in_data[istep], a_integrand, a_geom,
                          a_use_area_element, a_method_u, a_method_v);
    }
    std::cout << "\n";
}

#endif /* SURFACEEXTRACTIONDATA_IMPL_HPP */
