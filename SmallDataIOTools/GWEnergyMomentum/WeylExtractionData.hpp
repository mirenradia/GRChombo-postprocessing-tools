/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEYLEXTRACTIONDATA_HPP
#define WEYLEXTRACTIONDATA_HPP

#include "SphericalGeometry.hpp"
#include "SurfaceExtractionData.hpp"
#include <array>
#include <cmath>

//! This child class of SurfaceExtractionData reads in the Weyl4ExtractionOut
//! data and uses it to calculate the radiated energy, momentum and angular
//! momentum
class WeylExtractionData : public SurfaceExtractionData
{
  protected:
    SphericalGeometry m_geom = {
        {0.0, 0.0, 0.0}}; // center is irrelevant so just construct with zeros
    extracted_data_t
        m_time_integrated_data; // The once integrated in time Weyl4 data
    bool m_time_integrated_once = false;
    extracted_data_t
        m_time_twice_integrated_data; // The twice integrated in time Weyl4 data
    bool m_time_integrated_twice = false;
    time_multisurface_value_t m_power;
    bool m_power_computed = false;
    multisurface_value_t m_total_energy;
    std::array<time_multisurface_value_t, CH_SPACEDIM> m_momentum_flux;
    bool m_momentum_flux_computed = false;
    std::array<multisurface_value_t, CH_SPACEDIM> m_total_momentum;
    std::array<time_multisurface_value_t, CH_SPACEDIM> m_angular_momentum_flux;
    bool m_angular_momentum_flux_computed = false;
    std::array<multisurface_value_t, CH_SPACEDIM> m_total_angular_momentum;

  public:
    // use SurfaceExtractionData constructor
    using SurfaceExtractionData::SurfaceExtractionData;

  protected:
    void integrate_time_once()
    {
        assert(m_data_read);
        if (!m_time_integrated_once)
        {
            integrate_all_time(m_time_integrated_data, m_data);
        }
        m_time_integrated_once = true;
    }

    void integrate_time_twice()
    {
        assert(m_data_read);
        // do first time integration (does nothing if done already)
        integrate_time_once();
        if (!m_time_integrated_twice)
        {
            integrate_all_time(m_time_twice_integrated_data,
                               m_time_integrated_data);
        }
        m_time_integrated_twice = true;
    }

  public:
    void compute_power()
    {
        if (!m_power_computed)
        {
            // integrate Weyl4 in time (does nothing if already done)
            integrate_time_once();

            // power integrand excluding area element
            auto power_integrand = [](std::vector<double> integrated_Weyl4,
                                      double r, double, double) {
                return (integrated_Weyl4[0] * integrated_Weyl4[0] +
                        integrated_Weyl4[1] * integrated_Weyl4[1]) /
                       (16.0 * M_PI);
            };

            integrate_surface(m_power, m_time_integrated_data, power_integrand,
                              m_geom, 0, m_num_steps - 1,
                              IntegrationMethod::simpson);
        }
        m_power_computed = true;
    }

    void compute_energy()
    {
        // compute power if not done already
        compute_power();

        // integrate in time
        integrate_time(m_total_energy, m_power, 0, m_num_steps - 1);
    }

    void compute_momentum_flux()
    {
        if (!m_momentum_flux_computed)
        {
            // integrate Weyl4 in time (does nothing if already done)
            integrate_time_once();

            // momentum integrands excluding area elements
            auto momentum_integrand_x = [](std::vector<double> integrated_Weyl4,
                                           double r, double theta, double phi) {
                return std::sin(theta) * std::cos(phi) *
                       (integrated_Weyl4[0] * integrated_Weyl4[0] +
                        integrated_Weyl4[1] * integrated_Weyl4[1]) /
                       (16.0 * M_PI);
            };
            auto momentum_integrand_y = [](std::vector<double> integrated_Weyl4,
                                           double r, double theta, double phi) {
                return std::sin(theta) * std::sin(phi) *
                       (integrated_Weyl4[0] * integrated_Weyl4[0] +
                        integrated_Weyl4[1] * integrated_Weyl4[1]) /
                       (16.0 * M_PI);
            };
            auto momentum_integrand_z = [](std::vector<double> integrated_Weyl4,
                                           double r, double theta, double phi) {
                return std::cos(theta) *
                       (integrated_Weyl4[0] * integrated_Weyl4[0] +
                        integrated_Weyl4[1] * integrated_Weyl4[1]) /
                       (16.0 * M_PI);
            };
            const std::array<integrand_t, CH_SPACEDIM> momentum_integrands{
                momentum_integrand_x, momentum_integrand_y,
                momentum_integrand_z};

            for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
                integrate_surface(m_momentum_flux[idir], m_time_integrated_data,
                                  momentum_integrands[idir], m_geom, 0,
                                  m_num_steps - 1, IntegrationMethod::simpson);
            }
        }
        m_momentum_flux_computed = true;
    }

    void compute_momentum()
    {
        // compute momentum flux (does nothing if done already)
        compute_momentum_flux();

        // integrate in time
        for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
            integrate_time(m_total_momentum[idir], m_momentum_flux[idir], 0,
                           m_num_steps - 1);
        }
    }

    void compute_angular_momentum_flux()
    {
        if (!m_angular_momentum_flux_computed)
        {
            // integrate Weyl4 in time twice (does nothing if done already)
            integrate_time_twice();
        }
        m_angular_momentum_flux_computed = true;
    }
};

#endif /* WEYLEXTRACTIONDATA_HPP */
