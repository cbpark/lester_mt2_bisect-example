/* Copyright 2021 Chan Beom Park <cbpark@gmail.com> */

#include <cmath>  // std::sqrt
#include <iostream>
#include <utility>  // std::pair
#include "lester_mt2_bisect.h"
#include "maos.h"

using std::cout;

int main() {
    // decay topology: Y + Y --> vis_a(p1) chi(k1) + vis_b(p2) chi(k2).
    double vis_a[4] = {100.0, 410.0, 20.0, -100.0};   // {Mass, p1x, p1y, p1z}
    double vis_b[4] = {150.0, -210.0, -300.0, 50.0};  // {Mass, p2x, p2y, p2z}
    double ptmiss[2] = {-200.0, 280.0};               // {ptmiss_x, ptmiss_y}
    double m_chi = 100.0;                             // invisible particle mass

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
        ptmiss[1], m_chi, m_chi);
    cout << "MT2: " << mt2 << '\n';

    // the scale factor to make the parameters massless.
    // it'll be used for calculating the longitudinal momentum
    double scale =
        (vis_a[0] * vis_a[0] + vis_a[1] * vis_a[1] + vis_a[2] * vis_a[2] +
         vis_b[0] * vis_b[0] + vis_b[1] * vis_b[1] + vis_b[2] * vis_b[2] +
         ptmiss[0] * ptmiss[0] + ptmiss[1] * ptmiss[1] + m_chi * m_chi +
         m_chi * m_chi) /
        8.0;
    scale = std::sqrt(scale);

    // the parent particle mass
    double m_parent = 500.0;

    cout << "with the MAOS solutions:\n";

    // this is also defined in 'lester_mt2_bisect.h'
    std::pair<double, double> sols =
        ben_findsols(mt2, vis_a[1], vis_a[2], vis_a[0], m_chi, vis_b[1],
                     vis_b[2], ptmiss[0], ptmiss[1], vis_b[0], m_chi);

    double k1x = sols.first;
    double k1y = sols.second;

    std::pair<double, double> k1z =
        maos::longitudinalMomentum(vis_a[0], vis_a[1], vis_a[2], vis_a[3], k1x,
                                   k1y, m_parent, m_chi, scale);

    // note that we have two-fold degenerate solutions for the longitudinal
    // momentum.
    cout << "  k1x = " << k1x << ", k1y = " << k1y << ", k1z(1) = " << k1z.first
         << ", k1z(2) = " << k1z.second << '\n';

    double k2x = ptmiss[0] - k1x;
    double k2y = ptmiss[1] - k1y;
    std::pair<double, double> k2z =
        maos::longitudinalMomentum(vis_b[0], vis_b[1], vis_b[2], vis_b[3], k2x,
                                   k2y, m_parent, m_chi, scale);

    cout << "  k2x = " << k2x << ", k2y = " << k2y << ", k2z(1) = " << k2z.first
         << ", k2z(2) = " << k2z.second << '\n';
}
