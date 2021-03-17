/* Copyright 2021 Chan Beom Park <cbpark@gmail.com> */

#include <cmath>  // std::sqrt
#include <iostream>
#include <utility>  // std::pair
#include "lester_mt2_bisect.h"
#include "maos.h"

using std::cout;

/*
 * input:
<event>
  9      1 +7.1353100e-02 1.05830100e+01 7.81860800e-03 1.76468800e-01
        11 -1    0    0    0    0 +0.0000000000e+00 +0.0000000000e+00 +6.9999999897e+00 7.0000000084e+00 5.1100000000e-04 0.0000e+00 -1.0000e+00
       -11 -1    0    0    0    0 -0.0000000000e+00 -0.0000000000e+00 -3.9999999626e+00 3.9999999952e+00 5.1100000000e-04 0.0000e+00 1.0000e+00
        15  2    1    2    0    0 +2.0842344003e+00 +3.7168232237e+00 -1.1871038314e+00 4.7671534922e+00 1.7770000000e+00 0.0000e+00 0.0000e+00
       -15  2    1    2    0    0 -2.0842344003e+00 -3.7168232237e+00 +4.1871038586e+00 6.2328465114e+00 1.7770000000e+00 0.0000e+00 0.0000e+00
        13  1    3    3    0    0 +2.3733444094e+00 +2.9713979118e+00 -1.1149733294e+00 3.9643787403e+00 1.0566000000e-01 0.0000e+00 1.0000e+00
      3000  1    3    3    0    0 -2.8911000907e-01 +7.4542531195e-01 -7.2130502072e-02 8.0277475195e-01 1.0000000000e-03 0.0000e+00 0.0000e+00
       -13  1    4    4    0    0 -1.7146298897e+00 -2.5772167705e+00 +2.6338524472e+00 4.0657526597e+00 1.0566000000e-01 0.0000e+00 1.0000e+00
       -16  1    4    4    0    0 -1.1415344293e-01 +1.1864710228e-01 +1.7841770000e-02 1.6560939638e-01 0.0000000000e+00 0.0000e+00 1.0000e+00
        14  1    4    4    0    0 -2.5545106770e-01 -1.2582535555e+00 +1.5354096414e+00 2.0014844553e+00 0.0000000000e+00 0.0000e+00 -1.0000e+00
</event>

 * output 1 (m_chi1 = 0.0, m_chi2 = 0.0):
MT2: 0.332462
with the MAOS solutions:
  k1x = 1.34024, k1y = 2.06907, k1z(1) = 0.682739, k1z(2) = -2.36249
  k2x = -1.99896, k2y = -2.46325, k2z(1) = 5.5232, k2z(2) = 0.708697

 * output 2 (m_chi1 = 1.0e-3, m_chi2 = 0.921):
MT2: 1.02668
with the MAOS solutions:
  k1x = 14.0665, k1y = 21.716, k1z(1) = -3.71196, k1z(2) = -11.6161
  k2x = -14.7252, k2y = -22.1101, k2z(1) = 28.4877, k2z(2) = 17.2971
 */
int main() {
    // decay topology: Y + Y --> vis_a(p1) chi(k1) + vis_b(p2) chi(k2).
    // {Mass, p1x, p1y, p1z}
    double vis_a[4] =
        {0.10566, 2.3733444094, 2.9713979118, -1.1149733294};
    // {Mass, p2x, p2y, p2z}
    double vis_b[4] =
        {0.10566, -1.7146298897, -2.5772167705, 2.6338524472};
    // (ptmiss_x, ptmiss_y)
    double ptmiss[2] =
        {-2.8911000907e-01 - 1.1415344293e-01 - 2.5545106770e-01,
         7.4542531195e-01 + 1.1864710228e-01 - 1.2582535555e+00};

    // double m_chi1 = 0.0, m_chi2 = 0.0;
    double m_chi1 = 1.0e-3;
    double m_chi2 = 0.921;

    asymm_mt2_lester_bisect::disableCopyrightMessage();
    /*
     * Note that the argumemts of the get_mT2 function are:
     *   asymm_mt2_lester_bisect::get_mT2(
     *       m_vis1, px_vis1, py_vis1,
     *       m_vis2, px_vis2, py_vis2,
     *       px_miss, py_miss, m_inv1, m_inv2);
     */
    double mt2 = asymm_mt2_lester_bisect::get_mT2(
        vis_a[0], vis_a[1], vis_a[2], vis_b[0], vis_b[1], vis_b[2], ptmiss[0],
        ptmiss[1], m_chi1, m_chi2);
    cout << "MT2: " << mt2 << '\n';

    // the scale factor to make the parameters massless.
    // it'll be used for calculating the longitudinal momentum
    double scale =
        (vis_a[0] * vis_a[0] + vis_a[1] * vis_a[1] + vis_a[2] * vis_a[2] +
         vis_b[0] * vis_b[0] + vis_b[1] * vis_b[1] + vis_b[2] * vis_b[2] +
         ptmiss[0] * ptmiss[0] + ptmiss[1] * ptmiss[1] + m_chi1 * m_chi1 +
         m_chi2 * m_chi2) /
        8.0;
    scale = std::sqrt(scale);

    // the parent particle mass
    double m_parent = 1.777;

    cout << "with the MAOS solutions:\n";

    // this is also defined in 'lester_mt2_bisect.h'
    std::pair<double, double> sols =
        ben_findsols(mt2, vis_a[1], vis_a[2], vis_a[0], m_chi1, vis_b[1],
                     vis_b[2], ptmiss[0], ptmiss[1], vis_b[0], m_chi2);

    double k1x = sols.first;
    double k1y = sols.second;

    std::pair<double, double> k1z =
        maos::longitudinalMomentum(vis_a[0], vis_a[1], vis_a[2], vis_a[3], k1x,
                                   k1y, m_parent, m_chi1, scale);

    // note that we have two-fold degenerate solutions for the longitudinal
    // momentum.
    cout << "  k1x = " << k1x << ", k1y = " << k1y << ", k1z(1) = " << k1z.first
         << ", k1z(2) = " << k1z.second << '\n';

    double k2x = ptmiss[0] - k1x;
    double k2y = ptmiss[1] - k1y;
    std::pair<double, double> k2z =
        maos::longitudinalMomentum(vis_b[0], vis_b[1], vis_b[2], vis_b[3], k2x,
                                   k2y, m_parent, m_chi2, scale);

    cout << "  k2x = " << k2x << ", k2y = " << k2y << ", k2z(1) = " << k2z.first
         << ", k2z(2) = " << k2z.second << '\n';
}
