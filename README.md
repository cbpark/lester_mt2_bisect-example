A simple example code for using bisection-based asymmetric MT2 library.

See

* C.G. Lester and B. Nachman, _Bisection-based asymmetric MT2 computation: a higher precision calculator than existing symmetric methods_, [arXiv:1411.4312](https://arxiv.org/abs/1411.4312).
* C.G. Lester, [MT2 / Stransverse Mass / Oxbridge Kinetics Library](https://www.hep.phy.cam.ac.uk/~lester/mt2/).

For the MAOS method,

* W.S. Cho, K. Choi, Y.G. Kim, C.B. Park, _MT2-assisted on-shell reconstruction of missing momenta and its application to spin measurement at the LHC_, [arXiv:0810.4853](https://arxiv.org/abs/0810.4853).

``` no-hightlight
$ make && ./mt2
MT2: 412.628
with the MAOS solutions:
  k1x = 177.061, k1y = 301.832, k1z(1) = 177.774, k1z(2) = -394.728
  k2x = -377.061, k2y = -21.8318, k2z(1) = 362.588, k2z(2) = -238.398
```
