#include <iostream>
#include "lester_mt2_bisect.h"

int main() {
    double vis_a[3] = {100.0, 410.0, 20.0};
    double vis_b[4] = {150.0, -210.0, -300.0};
    double ptmiss[2] = {-200.0, 280.0};
    double m_chi = 100.0;

    asymm_mt2_lester_bisect::disableCopyrightMessage();
    // Note that the argumemts of the get_mT2 function are:
    //   asymm_mt2_lester_bisect::get_mT2(
    //       m_vis1, px_vis1, py_vis1, m_vis2,
    //       px_vis2, py_vis2, px_miss, py_miss,
    //       m_invis1, m_invis2);
    double mt2 = asymm_mt2_lester_bisect::get_mT2(
        vis_a[0], vis_a[1], vis_a[2], vis_b[0], vis_b[1], vis_b[2], ptmiss[0],
        ptmiss[1], m_chi, m_chi);
    std::cout << "MT2: " << mt2 << '\n';
}
