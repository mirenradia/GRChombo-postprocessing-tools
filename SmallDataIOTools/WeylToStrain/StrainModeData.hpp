/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef STRAINMODEDATA_HPP_
#define STRAINMODEDATA_HPP_

#include "Misc.H"
#include "WeylModeData.hpp"
// Requires linking with fftw3
#include <algorithm>
#include <cmath>
#include <complex>
#include <fftw3.h>

// Chombo Namespace
#include "UsingNamespace.H"

//! This class reads, stores and manipulates modes of the Weyl scalar, Psi4
class StrainModeData : public WeylModeData
{
  protected:
    std::vector<mode_data_t> m_strain_data;
    std::vector<double> m_cutoff_freq_lows;
    std::vector<bool> m_strain_computed;

  private:
    // The following members are only used by the fftw functions for the FFTs
    // and shouldn't need to be accessed outside of this class
    bool m_set_up_for_fft = false;
    int m_fft_length;
    int m_fft_clength; // since FT data is Hermitian, fftw only stores half of
                       // it
    double *m_weyl_data_ptr;
    std::complex<double> *m_f_weyl_data_ptr;
    std::complex<double> *m_f_strain_data_ptr;
    double *m_strain_data_ptr;
    fftw_plan m_fft_plan;
    fftw_plan m_ifft_plan;
    std::vector<double> m_freq;

  public:
    using WeylModeData::WeylModeData;

    ~StrainModeData()
    {
        if (m_set_up_for_fft)
        {
            fftw_destroy_plan(m_fft_plan);
            fftw_destroy_plan(m_ifft_plan);
            fftw_free(m_weyl_data_ptr);
            fftw_free(m_f_weyl_data_ptr);
            fftw_free(m_f_strain_data_ptr);
            fftw_free(m_strain_data_ptr);
            fftw_cleanup();
        }
    }

    // override this function as we need to resize m_strain_data
    void read_mode_data(const double extraction_radius_multiplier = 1.0)
    {
        // this actually does the reading
        WeylModeData::read_mode_data(extraction_radius_multiplier);
        m_strain_data.resize(m_num_modes);
        m_cutoff_freq_lows.resize(m_num_modes);
        m_strain_computed.resize(m_num_modes, false);
        for (auto &mode : m_strain_data)
        {
            mode.resize(2 * m_num_extraction_radii);
            for (auto &dataset : mode)
            {
                dataset.resize(m_num_steps);
            }
        }
        set_up_for_fft();
    }

    const mode_data_t &get_strain_data(const int l, const int m) const
    {
        const int mode_idx = get_mode_index(l, m);
        assert(m_strain_computed[mode_idx]);
        return m_strain_data[mode_idx];
    }

    double get_cutoff_freq_low(const int l, const int m) const
    {
        const int mode_idx = get_mode_index(l, m);
        assert(m_strain_computed[mode_idx]);
        return m_cutoff_freq_lows[mode_idx];
    }

    void set_up_for_fft()
    {
        assert(m_data_read);
        assert(!m_set_up_for_fft);
        // make the length of the FFT the power of 2 that is at least 4 times
        // as big as the number of steps as this is what is recommended by
        // numerical recipes - we will pad with zeros
        m_fft_length =
            ipow(2, static_cast<int>(std::ceil(std::log2(m_num_steps)) + 2));
        // this comes from the fftw documentation
        m_fft_clength = m_fft_length / 2 + 1;
        std::cout << "num_steps = " << m_num_steps
                  << ", fft_length = " << m_fft_length << "\n";
        m_weyl_data_ptr = fftw_alloc_real(m_fft_length);
        m_f_weyl_data_ptr = reinterpret_cast<std::complex<double> *>(
            fftw_alloc_complex(m_fft_clength));
        m_f_strain_data_ptr = reinterpret_cast<std::complex<double> *>(
            fftw_alloc_complex(m_fft_clength));
        m_strain_data_ptr = fftw_alloc_real(m_fft_length);
        m_freq.resize(m_fft_clength);
        for (int k = 0; k < m_fft_clength; ++k)
        {
            m_freq[k] = 2 * M_PI * k / (m_fft_length * m_time_data.get_dt());
        }

        // plan the FFT
        m_fft_plan = fftw_plan_dft_r2c_1d(
            m_fft_length, m_weyl_data_ptr,
            reinterpret_cast<fftw_complex *>(m_f_weyl_data_ptr), FFTW_ESTIMATE);

        // plan the inverse FFT
        m_ifft_plan = fftw_plan_dft_c2r_1d(
            m_fft_length, reinterpret_cast<fftw_complex *>(m_f_strain_data_ptr),
            m_strain_data_ptr, FFTW_ESTIMATE);

        m_set_up_for_fft = true;
    }

    // computes the strain by FFT, then applying Eq. (27) in arXiv:1006.1632
    // twice (since Psi4 is the second time derivative of h), and then inverse
    // FFT
    void compute_strain(const int l, const int m,
                        const double a_cutoff_freq_low,
                        const double a_cutoff_freq_high)
    {
        assert(m_set_up_for_fft);
        const mode_data_t &mode_data = get_mode_data(l, m);
        const int mode_idx = get_mode_index(l, m);
        m_cutoff_freq_lows[mode_idx] = a_cutoff_freq_low;
        std::vector<double> freq_with_cutoffs = m_freq;

        // For the FFI in arXiv:1006.1632, need to replace frequencies
        // below cutoff frequency with cutoff when integrating
        auto cutoff_low_it =
            std::upper_bound(freq_with_cutoffs.begin(), freq_with_cutoffs.end(),
                             a_cutoff_freq_low);
        std::fill(freq_with_cutoffs.begin(), cutoff_low_it, a_cutoff_freq_low);
        auto cutoff_high_it = std::upper_bound(
            cutoff_low_it, freq_with_cutoffs.end(), a_cutoff_freq_high);
        std::fill(cutoff_high_it, freq_with_cutoffs.end(), a_cutoff_freq_high);

        for (int icol = 0; icol < 2 * m_num_extraction_radii; ++icol)
        {
            for (int istep = 0; istep < m_num_steps; ++istep)
            {
                m_weyl_data_ptr[istep] = mode_data[icol][istep];
            }
            // pad with zeros here
            for (int istep = m_num_steps; istep < m_fft_length; ++istep)
            {
                m_weyl_data_ptr[istep] = 0.0;
            }

            // do the FFT
            fftw_execute(m_fft_plan);

            // now integrate up
            for (int istep = 0; istep < m_fft_clength; ++istep)
            {
                m_f_strain_data_ptr[istep] =
                    -m_f_weyl_data_ptr[istep] /
                    (freq_with_cutoffs[istep] * freq_with_cutoffs[istep]);
            }

            // do the inverse FFT
            fftw_execute(m_ifft_plan);

            // the cross strain is minus the imaginary part
            double sign = (icol % 2 == 0) ? 1 : -1;
            // now put back into m_strain_data
            for (int istep = 0; istep < m_num_steps; ++istep)
            {
                m_strain_data[mode_idx][icol][istep] =
                    sign * m_strain_data_ptr[istep] / m_fft_length;
            }
        }
        m_strain_computed[mode_idx] = true;
    }
};

#endif /* STRAINMODEDATA_HPP_ */
