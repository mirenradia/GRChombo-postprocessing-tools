/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEYLEXTRACTIONDATA_HPP
#define WEYLEXTRACTIONDATA_HPP

#include "SmallDataIO.hpp"
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

    // write flux and time integrated values to a file
    void write_data(const std::string &a_filename,
                    const std::string &a_time_data_label,
                    const time_multisurface_value_t &a_time_data,
                    const std::string &a_single_data_label = "",
                    const multisurface_value_t &a_single_data = {}) const
    {
        const double fake_time = 0.0;
        // this will overwrite any existing files
        const bool first_step = true;
        SmallDataIO file(a_filename, m_dt, fake_time, fake_time,
                         SmallDataIO::APPEND, true);

        std::vector<std::string> header1_strings(m_num_surfaces,
                                                 a_time_data_label);
        std::vector<std::string> header2_strings(m_num_surfaces);
        for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
        {
            auto surface_it = m_surfaces.begin();
            std::advance(surface_it, isurface);
            header2_strings[isurface] = std::to_string(surface_it->second);
        }
        std::string pre_header2_string = m_geom.param_name() + " = ";

        // write initial header lines
        file.write_header_line(header1_strings);
        file.write_header_line(header2_strings, pre_header2_string);

        for (int istep = 0; istep < m_num_steps; ++istep)
        {
            double time = m_time_limits.first + istep * m_dt;
            file.write_data_line(a_time_data[istep], time);
        }

        // write single data in a comment (header) line
        if (a_single_data.size() == m_num_surfaces)
        {
            file.line_break();
            std::vector<std::string> single_data_strs(m_num_surfaces);
            char single_data_cstr[20];
            std::string format =
                "%." + std::to_string(m_file_structure.data_width - 10) + "e";
            for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
            {
                std::sprintf(single_data_cstr, format.c_str(),
                             a_single_data[isurface]);
                single_data_strs[isurface] = single_data_cstr;
            }
            file.write_header_line(single_data_strs, a_single_data_label);
        }
    }

    void write_data(
        const std::string &a_filename,
        const std::array<std::string, CH_SPACEDIM> &a_time_data_labels,
        const std::array<time_multisurface_value_t, CH_SPACEDIM> &a_time_data,
        const std::string &a_single_data_label = "",
        const std::array<multisurface_value_t, CH_SPACEDIM> &a_single_data = {})
        const
    {
        const double fake_time = 0.0;
        // this will overwrite any existing files
        const bool first_step = true;
        SmallDataIO file(a_filename, m_dt, fake_time, fake_time,
                         SmallDataIO::APPEND, true);

        std::vector<std::string> header1_strings(CH_SPACEDIM * m_num_surfaces);
        std::vector<std::string> header2_strings(CH_SPACEDIM * m_num_surfaces);
        for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
        {
            for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
                const int idx = isurface * CH_SPACEDIM + idir;
                header1_strings[idx] = a_time_data_labels[idir];
                auto surface_it = m_surfaces.begin();
                std::advance(surface_it, isurface);
                header2_strings[idx] = std::to_string(surface_it->second);
            }
        }
        std::string pre_header2_string = m_geom.param_name() + " = ";

        // write initial header lines
        file.write_header_line(header1_strings);
        file.write_header_line(header2_strings, pre_header2_string);

        for (int istep = 0; istep < m_num_steps; ++istep)
        {
            std::vector<double> data_for_writing(CH_SPACEDIM * m_num_surfaces);
            for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
            {
                for (int idir = 0; idir < CH_SPACEDIM; ++idir)
                {
                    const int idx = isurface * CH_SPACEDIM + idir;
                    data_for_writing[idx] = a_time_data[idir][istep][isurface];
                }
            }
            double time = m_time_limits.first + istep * m_dt;
            file.write_data_line(data_for_writing, time);
        }

        // write single data in a comment (header) line
        if (a_single_data[0].size() == m_num_surfaces)
        {
            file.line_break();
            std::vector<std::string> single_data_strs(CH_SPACEDIM *
                                                      m_num_surfaces);
            std::string format =
                "%." + std::to_string(m_file_structure.data_width - 10) + "e\0";

            for (int isurface = 0; isurface < m_num_surfaces; ++isurface)
            {
                for (int idir = 0; idir < CH_SPACEDIM; ++idir)
                {
                    char single_data_cstr[20];
                    std::sprintf(single_data_cstr, "%.10e",
                                 a_single_data[idir][isurface]);
                    single_data_strs[isurface * CH_SPACEDIM + idir] =
                        single_data_cstr;
                }
            }
            file.write_header_line(single_data_strs, a_single_data_label);
            std::cout << "\n";
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

    void compute_energy()
    {
        // compute power if not done already
        compute_power();

        // integrate in time
        integrate_time(m_total_energy, m_power, 0, m_num_steps - 1);

        // print result
        print_energy();

        // write result to file
        write_data("GW_power", "power", m_power, "energy = ", m_total_energy);
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

        // write data
        write_data("GW_momentum_flux", {"x flux", "y flux", "z flux"},
                   m_momentum_flux, "total = ", m_total_momentum);
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

        // write data
        write_data("GW_angular_momentum_flux", {"x flux", "y flux", "z flux"},
                   m_angular_momentum_flux,
                   "total = ", m_total_angular_momentum);
    }
};

#endif /* WEYLEXTRACTIONDATA_HPP */
