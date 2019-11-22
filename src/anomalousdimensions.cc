//
// ePDF
//

#include "ePDF/anomalousdimensions.h"
#include "ePDF/specialfunctions.h"
#include "ePDF/constants.h"

namespace ePDF
{
  //_________________________________________________________________________________
  std::pair<std::complex<double>, Matrix<std::complex<double>>> andim_lo(std::complex<double> const& N, int const& nf)
  {
    const std::complex<double> Ns = N * N;
    const std::complex<double> N1 = N + 1.;
    const std::complex<double> N2 = N + 2.;
    const std::complex<double> Nm = N - 1.;

    const std::complex<double> s1 = emc + psi(N1);

    const std::complex<double> pqqa = 3. - 4. * s1 + 2. / N / N1;
    const std::complex<double> pqga = 4. * ( Ns + N + 2. ) / N / N1 / N2;
    const std::complex<double> pgqa = 2. * ( Ns + N + 2. ) / N / N1 / Nm;
    const std::complex<double> pgga = 11. / 3. - 4. * s1 + 4. / N / Nm + 4. / N1 / N2;
    const std::complex<double> pggb = - 4. / 3.;

    // Output to the array
    const std::complex<double>         p0ns = CF * pqqa;
    const Matrix<std::complex<double>> p0sg{2, 2, {CF * pqqa, TR * nf * pqga, CF * pgqa, CA * pgga + TR * nf * pggb}};
    return std::make_pair(p0ns, p0sg);
  }

  //_________________________________________________________________________________
  std::pair<std::complex<double>, Matrix<std::complex<double>>> andim_nlo(std::complex<double> const& N, int const& nf)
  {
    const std::complex<double> s1 = emc + psi(N + 1.);
    const std::complex<double> s2 = zeta2 - dpsi(N + 1., 1);

    const std::complex<double> Ns  = N * N;
    const std::complex<double> Nt  = Ns * N;
    const std::complex<double> Nfo = Nt * N;
    const std::complex<double> Nfi = Nfo * N;
    const std::complex<double> Nsi = Nfi * N;
    const std::complex<double> Nse = Nsi * N;
    const std::complex<double> Ne  = Nse * N;
    const std::complex<double> Nn  = Ne * N;

    const std::complex<double> Nm  = N - 1.;
    const std::complex<double> N1  = N + 1.;
    const std::complex<double> N2  = N + 2.;
    const std::complex<double> Nms = Nm * Nm;
    const std::complex<double> N1s = N1 * N1;
    const std::complex<double> N1t = N1s * N1;
    const std::complex<double> N2s = N2 * N2;
    const std::complex<double> N2t = N2s * N2;

    // Analytic continuations of the occuring sums as given in GRV
    // (1990) (with an improved parametrization of the moments of
    // Sp(x)/(1+x).)
    const std::complex<double> N3 = N + 3.;
    const std::complex<double> N4 = N + 4.;
    const std::complex<double> N5 = N + 5.;
    const std::complex<double> N6 = N + 6.;
    const std::complex<double> s11 = s1  + 1. / N1;
    const std::complex<double> s12 = s11 + 1. / N2;
    const std::complex<double> s13 = s12 + 1. / N3;
    const std::complex<double> s14 = s13 + 1. / N4;
    const std::complex<double> s15 = s14 + 1. / N5;
    const std::complex<double> s16 = s15 + 1. / N6;
    const std::complex<double> spmom = 1.0000 * ( zeta2 - s1 / N ) / N  -
                                       0.9992 * ( zeta2 - s11/ N1 ) / N1 +
                                       0.9851 * ( zeta2 - s12/ N2 ) / N2 -
                                       0.9005 * ( zeta2 - s13/ N3 ) / N3 +
                                       0.6621 * ( zeta2 - s14/ N4 ) / N4 -
                                       0.3174 * ( zeta2 - s15/ N5 ) / N5 +
                                       0.0699 * ( zeta2 - s16/ N6 ) / N6;

    const std::complex<double> slc = - 5./8. * zeta3;
    const std::complex<double> slv = - zeta2 / 2.* ( psi(N1 / 2.) - psi(N / 2.) ) + s1 / Ns + spmom;
    const std::complex<double> sschlm = slc - slv;
    const std::complex<double> sstr2m = zeta2 - dpsi(N1 / 2., 1);
    const std::complex<double> sstr3m = 0.5 * dpsi (N1 / 2.,2) + zeta3;

    const std::complex<double> sschlp = slc + slv;
    const std::complex<double> sstr2p = zeta2 - dpsi(N2 / 2., 1);
    const std::complex<double> sstr3p = 0.5 * dpsi(N2 / 2., 2) + zeta3;

    // The contributions to P1NS as given in Gonzalez-Arroyo et
    // al. (1979) (Note that the anomalous dimensions in the
    // literature often differ from these moments of the splitting
    // functions by factors -1 or -2, in addition to possible
    // different normalizations of the coupling)
    const std::complex<double> pnma = ( 16.* s1 * (2.* N + 1.) / (Ns * N1s) +
                                        16.* (2.* s1 - 1./(N * N1)) * ( s2 - sstr2m ) +
                                        64.* sschlm + 24.* s2 - 3. - 8.* sstr3m -
                                        8.* (3.* Nt + Ns -1.) / (Nt * N1t) +
                                        16.* (2.* Ns + 2.* N +1.) / (Nt * N1t) ) * (-0.5);
    const std::complex<double> pnpa = ( 16.* s1 * (2.* N + 1.) / (Ns * N1s) +
                                        16.* (2.* s1 - 1./(N * N1)) * ( s2 - sstr2p ) +
                                        64.* sschlp + 24.* s2 - 3. - 8.* sstr3p -
                                        8.* (3.* Nt + Ns -1.) / (Nt * N1t) -
                                        16.* (2.* Ns + 2.* N +1.)/(Nt * N1t) ) * (-0.5);

    const std::complex<double> pnsb = ( s1 * (536./9. + 8.* (2.* N + 1.) / (Ns * N1s)) -
                                        (16.* s1 + 52./3.- 8./(N * N1)) * s2 - 43./6. -
                                        (151.* Nfo + 263.* Nt + 97.* Ns + 3.* N + 9.) *
                                        4./ (9.* Nt * N1t) ) * (-0.5);
    const std::complex<double> pnsc = ( -160./9.* s1 + 32./3.* s2 + 4./3. +
                                        16.*(11.*Ns+5.*N-3.)/(9.* Ns * N1s))*(-0.5);

    // The contributions to P1SG as given in Floratos et al. (1981)
    // Pure singlet (PS) and QG
    const std::complex<double> ppsa = (5.* Nfi + 32.* Nfo + 49.* Nt+38.* Ns + 28.* N + 8.)
                                      / (Nm * Nt * N1t * N2s) * 2.;

    const std::complex<double> pqga = (-2.* s1 * s1 + 2.* s2 - 2.* sstr2p)
                                      * (Ns + N + 2.) / (N * N1 * N2)
                                      + (8.* s1 * (2.* N + 3.)) / (N1s * N2s)
                                      + 2.* (Nn + 6.* Ne + 15.* Nse + 25.* Nsi + 36.* Nfi
                                             + 85.* Nfo + 128.* Nt + 104.* Ns + 64.* N + 16.)
                                      / (Nm * Nt * N1t * N2t);
    const std::complex<double> pqgb = (2.* s1 * s1 - 2.* s2 + 5.) * (Ns + N + 2.)
                                      / (N * N1 * N2) - 4.* s1 / Ns
                                      + (11.* Nfo + 26.* Nt + 15.* Ns + 8.* N + 4.)
                                      / (Nt * N1t * N2);

    // GQ and GG
    const std::complex<double> pgqa = (- s1 * s1 + 5.* s1 - s2) * (Ns + N + 2.)
                                      / (Nm * N * N1)  -  2.* s1 / N1s
                                      - (12.* Nsi + 30.* Nfi + 43.* Nfo + 28.* Nt - Ns
                                         - 12.* N - 4.) / (2.* Nm * Nt * N1t);
    const std::complex<double> pgqb = (s1*s1 + s2 - sstr2p) * (Ns + N + 2.) / (Nm * N * N1)
                                      - s1 * (17.* Nfo + 41.* Ns - 22.* N - 12.)
                                      / (3.* Nms * Ns * N1);
    + (109.* Nn + 621.* Ne + 1400.* Nse + 1678.* Nsi
      + 695.* Nfi - 1031.* Nfo - 1304.* Nt - 152.* Ns
      + 432.* N + 144.) / (9.* Nms * Nt * N1t * N2s);
    const std::complex<double> pgqc = 4./3.* ( (s1 - 8./3.) * (Ns + N + 2.) / (Nm * N * N1) + 1./ N1s );

    const std::complex<double> pgga = - (2.* Nfi + 5.* Nfo + 8.* Nt + 7.* Ns- 2.* N - 2.)
                                      * 8.* s1 / (Nms * Ns * N1s * N2s) -  67./9.* s1 + 8./3.
                                      - 4.* sstr2p * (Ns + N + 1.) / (Nm * N * N1 * N2)
                                      + 2.* s1 * sstr2p - 4.* sschlp + 0.5 * sstr3p
                                      + (457.* Nn + 2742.* Ne + 6040.* Nse + 6098.* Nsi
                                         + 1567.* Nfi - 2344.* Nfo - 1632.* Nt + 560.* Ns
                                         + 1488.* N + 576.) / (18.* Nms * Nt * N1t * N2t);
    const std::complex<double> pggb = (38.* Nfo + 76.* Nt + 94.* Ns + 56.* N + 12.) *(-2.)
                                      / (9.* Nm * Ns * N1s * N2)  +  20./9.* s1  -  4./3.;
    const std::complex<double> pggc = (2.* Nsi + 4.* Nfi + Nfo - 10.* Nt - 5.* Ns - 4.* N
                                       - 4.) * (-2.) / (Nm * Nt * N1t * N2) - 1.;

    // Valence
    const std::complex<double>         p1ns = CF *((CF-CA/2.)* pnma + CA* pnsb + TR*nf* pnsc);
    const Matrix<std::complex<double>> p1sg{2, 2, {CF *((CF-CA/2.)* pnpa + CA* pnsb + TR*nf* pnsc) + TR*nf*CF*ppsa*4., TR*nf * (CA * pqga + CF * pqgb)*4.,
                                                   (CF*CF*pgqa + CF*CA*pgqb+TR*nf*CF*pgqc)*4., (CA*CA*pgga + TR*nf*(CA*pggb+CF*pggc))*4.
                                                  }};
    return std::make_pair(p1ns, p1sg);
  }
}
