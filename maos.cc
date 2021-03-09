/* Copyright 2021 Chan Beom Park <cbpark@gmail.com */

#include "maos.h"
#include <cmath>    // std::fabs, std::sqrt
#include <limits>   // std::numerical_limits
#include <utility>  // std::pair
#ifdef DEBUG
#include <iostream>
using std::cout;
#endif

namespace maos {
std::pair<double, double> longitudinalMomentum(double m_vis, double vis_px,
                                               double vis_py, double vis_pz,
                                               double inv_kx, double inv_ky,
                                               double m_parent, double m_inv,
                                               double scale) {
#ifdef DEBUG
    cout << "DEBUG: scale = " << scale << '\n';
#endif
    // set all the dimensionful parameters to be massless.
    m_vis /= scale;
    vis_px /= scale;
    vis_py /= scale;
    vis_pz /= scale;
    inv_kx /= scale;
    inv_ky /= scale;
    m_parent /= scale;
    m_inv /= scale;

    const double m_inv_sq = m_inv * m_inv;
    const double m_vis_sq = m_vis * m_vis;
    const double d = 0.5 * (m_parent * m_parent - m_inv_sq - m_vis_sq) +
                     vis_px * inv_kx + vis_py * inv_ky;

    const double et_vis_sq = m_vis_sq + vis_px * vis_px + vis_py * vis_py;
#ifdef DEBUG
    cout << "DEBUG: et_vis_sq = " << et_vis_sq << '\n';
#endif

    // null momentum?
    if (et_vis_sq < std::numeric_limits<double>::epsilon()) {
        return {1.0e+10, 1.0e+10};
    }

    const double et_inv_sq = m_inv_sq + inv_kx * inv_kx + inv_ky * inv_ky;
#ifdef DEBUG
    cout << "DEBUG: et_inv_sq = " << et_inv_sq << '\n';
#endif

    double disc = d * d - et_vis_sq * et_inv_sq;
#ifdef DEBUG
    cout << "DEBUG: disc = " << disc << '\n';
#endif

    // for numerical stability
    if (std::fabs(disc) < 1.0e-4) { disc = 0.0; }

    // unphysical momentum?
    if (disc < 0.0) { return {-1.0e+10, -1.0e+10}; }

    const double term1 = vis_pz * d;
    const double term2 =
        std::sqrt(et_vis_sq + vis_pz * vis_pz) * std::sqrt(disc);
#ifdef DEUBG
    cout << "DEBUG: term1 = " << term1 << '\n'
         << "DEBUG: term2 = " << term2 << '\n';
#endif

    return {(term1 + term2) / et_vis_sq * scale,
            (term1 - term2) / et_vis_sq * scale};
}
}  // namespace maos
