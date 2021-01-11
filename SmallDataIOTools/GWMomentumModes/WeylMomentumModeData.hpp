/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEYLMOMENTUMMODEDATA_HPP_
#define WEYLMOMENTUMMODEDATA_HPP_

#include "WeylModeData.hpp"

//! This class reads, stores and manipulates modes of the Weyl scalar, Psi4
class WeylMomentumModeData : public WeylModeData
{
  public:
    using WeylModeData::WeylModeData;

    // computes Eq. (8.9.69) in Alcubierre
    inline double momentum_a_coeff(const int l, const int m) const
    {
        return std::sqrt(static_cast<double>((l - m) * (l + m + 1))) /
               static_cast<double>((l * (l + 1)));
    }

    // computes Eq. (8.9.70) in Alcubierre
    inline double momentum_b_coeff(const int l, const int m) const
    {
        double numerator = (l - 2) * (l + 2) * (l + m) * (l + m - 1);
        double denominator = (2 * l - 1) * (2 * l + 1);
        return std::sqrt(numerator / denominator) / static_cast<double>(2 * l);
    }

    // computes Eq. (8.9.71) in Alcubierre
    inline double momentum_c_coeff(const int l, const int m) const
    {
        return static_cast<double>(2 * m) / static_cast<double>(l * (l + 1));
    }

    // computes Eq. (8.9.72) in Alcubierre
    inline double momentum_d_coeff(const int l, const int m) const
    {
        double numerator = (l - 2) * (l + 2) * (l - m) * (l + m);
        double denominator = (2 * l - 1) * (2 * l + 1);
        return std::sqrt(numerator / denominator) / static_cast<double>(l);
    }

    // check that we have the correct modes of Psi4 in order to compute the
    // desired term in Eqs (8.9.67-8.9.68) in Alcubierre
    bool check_momentum_pmode_computable(const int l, const int m)
    {
        bool missing_mode = (get_mode_index(l, m) < 0);
        missing_mode |= (get_mode_index(l, m + 1) < 0) && (abs(m + 1) <= l);
        missing_mode |=
            (get_mode_index(l - 1, m + 1) < 0) && (abs(m + 1) < l && l > 2);
        missing_mode |= (get_mode_index(l + 1, m + 1) < 0);
        missing_mode |= (get_mode_index(l - 1, m) < 0) && (abs(m) < l && l > 2);
        missing_mode |= (get_mode_index(l + 1, m) < 0);
        return !missing_mode;
    }

    // same as above but checks all m for single l
    bool check_momentum_pmode_computable(const int l)
    {
        bool missing_pmode = false;
        for (int m = -l; m <= l; ++m)
        {
            missing_pmode |= !check_momentum_pmode_computable(l, m);
        }
        return !missing_pmode;
    }

    // computes Eqs (8.9.67-8.9.68) in Alcubierre
    time_multidata_t compute_momentum_flux_pmode(const int l, const int m,
                                                 const double r_ex)
    {
        assert(l > 1 && abs(m) <= l);
        assert(check_momentum_pmode_computable(l, m));
        // uses the notation from Alcubierre but with _ instead of - and 1
        // instead of +1
        // any non-existent modes with |m|>l or l<2 are just set to zero
        // modes are already time integrated as well
        const mode_data_t &A_lm = get_time_integrated_mode_data(l, m);
        const mode_data_t &A_lm1 = get_time_integrated_mode_data(l, m + 1);
        const mode_data_t &A_l_1m1 =
            get_time_integrated_mode_data(l - 1, m + 1);
        const mode_data_t &A_l1m1 = get_time_integrated_mode_data(l + 1, m + 1);
        const mode_data_t &A_l_1m = get_time_integrated_mode_data(l - 1, m);
        const mode_data_t &A_l1m = get_time_integrated_mode_data(l + 1, m);

        const double a_lm = momentum_a_coeff(l, m);
        const double b_l_m = momentum_b_coeff(l, -m);
        const double b_l1m1 = momentum_b_coeff(l + 1, m + 1);
        const double c_lm = momentum_c_coeff(l, m);
        const double d_lm = momentum_d_coeff(l, m);
        const double d_l1m = momentum_d_coeff(l + 1, m);

        constexpr int dim = 3;
        time_multidata_t out(dim * m_num_extraction_radii,
                             std::vector<double>(m_num_steps));

        for (int iradius = 0; iradius < m_num_extraction_radii; ++iradius)
        {
            // indices for real and imaginary parts in mode data
            int i_re = 2 * iradius;
            int i_im = 2 * iradius + 1;
            // indices for momentum flux direction in output
            int i_x = dim * iradius;
            int i_y = dim * iradius + 1;
            int i_z = dim * iradius + 2;
            for (int istep = 0; istep < m_num_steps; ++istep)
            {
                out[i_x][istep] =
                    r_ex * r_ex *
                    (A_lm[i_re][istep] * (a_lm * A_lm1[i_re][istep] +
                                          b_l_m * A_l_1m1[i_re][istep] -
                                          b_l1m1 * A_l1m1[i_re][istep]) +
                     A_lm[i_im][istep] * (a_lm * A_lm1[i_im][istep] +
                                          b_l_m * A_l_1m1[i_im][istep] -
                                          b_l1m1 * A_l1m1[i_im][istep])) /
                    (8.0 * M_PI);
                out[i_y][istep] =
                    r_ex * r_ex *
                    (A_lm[i_re][istep] * (-a_lm * A_lm1[i_im][istep] -
                                          b_l_m * A_l_1m1[i_im][istep] +
                                          b_l1m1 * A_l1m1[i_im][istep]) +
                     A_lm[i_im][istep] * (a_lm * A_lm1[i_re][istep] +
                                          b_l_m * A_l_1m1[i_re][istep] -
                                          b_l1m1 * A_l1m1[i_re][istep])) /
                    (8.0 * M_PI);
                out[i_z][istep] =
                    r_ex * r_ex *
                    (A_lm[i_re][istep] * (c_lm * A_lm[i_re][istep] +
                                          d_lm * A_l_1m[i_re][istep] +
                                          d_l1m * A_l1m[i_re][istep]) +
                     A_lm[i_im][istep] * (c_lm * A_lm[i_im][istep] +
                                          d_lm * A_l_1m[i_im][istep] +
                                          d_l1m * A_l1m[i_im][istep])) /
                    (16.0 * M_PI);
            }
        }
        return out;
    }
};

#endif /* WEYLMOMENTUMMODEDATA_HPP_ */
