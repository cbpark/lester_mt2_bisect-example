/* Copyright 2021 Chan Beom Park <cbpark@gmail.com */

#ifndef LESTER_MT2_BISECT_EXAMPLE_MAOS_H_
#define LESTER_MT2_BISECT_EXAMPLE_MAOS_H_

#include <utility>  // std::pair

namespace maos {
/** returns the longitudinal momentum of the invisible particle.
 *
 *  See Eq.(8) in https://arxiv.org/pdf/0810.4853.pdf
 */
std::pair<double, double> longitudinalMomentum(double m_vis, double vis_px,
                                               double vis_py, double vis_pz,
                                               double inv_kx, double inv_ky,
                                               double m_parent, double m_inv,
                                               double scale);
}  // namespace maos

#endif  // LESTER_MT2_BISECT_EXAMPLE_MAOS_H_
