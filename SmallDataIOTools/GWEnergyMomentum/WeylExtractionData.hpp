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
            m_time_integrated_once = true;
        }
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
            m_time_integrated_twice = true;
        }
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
            bool use_area_element = true;
            integrate_surface(m_power, m_time_integrated_data, power_integrand,
                              m_geom, 0, m_num_steps - 1, use_area_element,
                              IntegrationMethod::simpson);
            m_power_computed = true;
        }
    }

    void print_energy()
    {
        std::cout << "Total energy radiated: \n";
        for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
        {
            auto surface_it = m_surfaces.begin();
            std::advance(surface_it, isurface);
            std::cout << std::fixed << std::setprecision(2)
                      << "at r = " << surface_it->second << ": ";
            std::cout << std::scientific
                      << std::setprecision(m_file_structure.data_width - 10)
                      << m_total_energy[isurface] << "\n";
        }
    }

    void compute_energy()
    {
        // compute power if not done already
        compute_power();

        // integrate in time
        integrate_time(m_total_energy, m_power, 0, m_num_steps - 1);

        // print result
        print_energy();
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

            bool use_area_element = true;
            for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
                integrate_surface(m_momentum_flux[idir], m_time_integrated_data,
                                  momentum_integrands[idir], m_geom, 0,
                                  m_num_steps - 1, use_area_element,
                                  IntegrationMethod::simpson);
            }
            m_momentum_flux_computed = true;
        }
    }

    void print_momentum()
    {
        std::cout << "Total momentum radiated: \n";
        for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
        {
            auto surface_it = m_surfaces.begin();
            std::advance(surface_it, isurface);
            std::cout << std::fixed << std::setprecision(2)
                      << "at r = " << surface_it->second << ": ";
            for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
                std::cout << std::scientific
                          << std::setprecision(m_file_structure.data_width - 10)
                          << m_total_momentum[idir][isurface] << " ";
            }
            std::cout << "\n";
        }
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

        // print momentum
        print_momentum();
    }

    void compute_angular_momentum_flux()
    {
        if (!m_angular_momentum_flux_computed)
        {
            // integrate Weyl4 in time twice (does nothing if done already)
            integrate_time_twice();
            auto d1_time_twice_integrated_data = compute_surface_derivatives(
                m_time_twice_integrated_data, m_geom);
            // this object will have dataset indices:
            // IWeyl4_Re IWeyl4_Im IIWeyl4_Re IIWeyl4_Im
            // du_IIWeyl4_Re dv_IIWeyl4_Re du_IIWeyl4_Im dv_IIWeyl4_Im
            // where I is a time integration and du and dv are derivatives wrt
            // theta and phi respectively
            auto combined_data = combine(
                combine(m_time_integrated_data, m_time_twice_integrated_data),
                d1_time_twice_integrated_data);

            // angular_momentum integrands including area elements
            auto angular_momentum_integrand_x =
                [](std::vector<double> combined_data_here, double r,
                   double theta, double phi) {
                    const double &IWeyl4_Re = combined_data_here[0];
                    const double &IWeyl4_Im = combined_data_here[1];
                    const double &IIWeyl4_Re = combined_data_here[2];
                    const double &IIWeyl4_Im = combined_data_here[3];
                    const double &dtheta_IIWeyl4_Re = combined_data_here[4];
                    const double &dphi_IIWeyl4_Re = combined_data_here[5];
                    const double &dtheta_IIWeyl4_Im = combined_data_here[6];
                    const double &dphi_IIWeyl4_Im = combined_data_here[7];
                    // spin weight
                    constexpr int s = -2;

                    // include sin(theta) from area element to remove
                    // singularities - expressions from Alcubierre p312
                    const double J_x_Re =
                        -std::sin(phi) * std::sin(theta) * dtheta_IIWeyl4_Re -
                        std::cos(phi) * (std::cos(theta) * dphi_IIWeyl4_Re -
                                         s * IIWeyl4_Im);
                    const double J_x_Im =
                        -std::sin(phi) * std::sin(theta) * dtheta_IIWeyl4_Im -
                        std::cos(phi) * (std::cos(theta) * dphi_IIWeyl4_Im -
                                         s * IIWeyl4_Re);

                    // only need real part
                    const double integrand =
                        -r * r * (IWeyl4_Re * J_x_Re + IWeyl4_Im * J_x_Im) /
                        (16.0 * M_PI);
                    return integrand;
                };
            auto angular_momentum_integrand_y =
                [](std::vector<double> combined_data_here, double r,
                   double theta, double phi) {
                    const double &IWeyl4_Re = combined_data_here[0];
                    const double &IWeyl4_Im = combined_data_here[1];
                    const double &IIWeyl4_Re = combined_data_here[2];
                    const double &IIWeyl4_Im = combined_data_here[3];
                    const double &dtheta_IIWeyl4_Re = combined_data_here[4];
                    const double &dphi_IIWeyl4_Re = combined_data_here[5];
                    const double &dtheta_IIWeyl4_Im = combined_data_here[6];
                    const double &dphi_IIWeyl4_Im = combined_data_here[7];
                    // spin weight
                    constexpr int s = -2;

                    // include sin(theta) from area element to remove
                    // singularities - expressions from Alcubierre p312
                    const double J_y_Re =
                        std::cos(phi) * std::sin(theta) * dtheta_IIWeyl4_Re -
                        std::sin(phi) * (std::cos(theta) * dphi_IIWeyl4_Re -
                                         s * IIWeyl4_Im);
                    const double J_y_Im =
                        std::cos(phi) * std::sin(theta) * dtheta_IIWeyl4_Im -
                        std::sin(phi) * (std::cos(theta) * dphi_IIWeyl4_Im -
                                         s * IIWeyl4_Re);

                    // only need real part
                    const double integrand =
                        -r * r * (IWeyl4_Re * J_y_Re + IWeyl4_Im * J_y_Im) /
                        (16.0 * M_PI);
                    return integrand;
                };
            auto angular_momentum_integrand_z =
                [](std::vector<double> combined_data_here, double r,
                   double theta, double phi) {
                    const double &IWeyl4_Re = combined_data_here[0];
                    const double &IWeyl4_Im = combined_data_here[1];
                    const double &dphi_IIWeyl4_Re = combined_data_here[5];
                    const double &dphi_IIWeyl4_Im = combined_data_here[7];
                    // spin weight
                    constexpr int s = -2;

                    // expressions from Alcubierre p312
                    const double J_z_Re = dphi_IIWeyl4_Re;
                    const double J_z_Im = dphi_IIWeyl4_Im;

                    // only need real part
                    // put sin(theta) in now
                    const double integrand =
                        -r * r * std::sin(theta) *
                        (IWeyl4_Re * J_z_Re + IWeyl4_Im * J_z_Im) /
                        (16.0 * M_PI);
                    return integrand;
                };
            const std::array<integrand_t, CH_SPACEDIM>
                angular_momentum_integrands{angular_momentum_integrand_x,
                                            angular_momentum_integrand_y,
                                            angular_momentum_integrand_z};

            bool use_area_element = false;
            for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
                integrate_surface(m_angular_momentum_flux[idir], combined_data,
                                  angular_momentum_integrands[idir], m_geom, 0,
                                  m_num_steps - 1, use_area_element,
                                  IntegrationMethod::simpson);
            }
            m_angular_momentum_flux_computed = true;
        }
    }

    void print_angular_momentum()
    {
        std::cout << "Total angular momentum radiated: \n";
        for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
        {
            auto surface_it = m_surfaces.begin();
            std::advance(surface_it, isurface);
            std::cout << std::fixed << std::setprecision(2)
                      << "at r = " << surface_it->second << ": ";
            for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
                std::cout << std::scientific
                          << std::setprecision(m_file_structure.data_width - 10)
                          << m_total_angular_momentum[idir][isurface] << " ";
            }
            std::cout << "\n";
        }
    }

    void compute_angular_momentum()
    {
        // compute flux (does nothing if done already)
        compute_angular_momentum_flux();

        // integrate in time
        for (int idir = 0; idir < CH_SPACEDIM; ++idir)
        {
            integrate_time(m_total_angular_momentum[idir],
                           m_angular_momentum_flux[idir], 0, m_num_steps - 1);
        }
        print_angular_momentum();
    }
};

#endif /* WEYLEXTRACTIONDATA_HPP */
