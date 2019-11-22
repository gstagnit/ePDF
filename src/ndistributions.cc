//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "ePDF/ndistributions.h"
#include "ePDF/anomalousdimensions.h"
#include "ePDF/specialfunctions.h"
#include "ePDF/constants.h"

namespace ePDF
{
  //_________________________________________________________________________________
  NDistributions::NDistributions(YAML::Node const& config):
    _aQED(AlphaQED{config}),
    _ipt(_ptmap.at(config["Perturbative order"].as<std::string>())),
    _nl(config["NL"].as<int>()),
    _Qi(config["Initial scale"].as<double>()),
    _ai(_aQED.Evolve(_Qi)),
    _bt0(_aQED.Beta0QED()),
    _bt1(_aQED.Beta1QED()),
    _method(config["PDF"]["evolution"].as<std::string>()),
    _expand(config["PDF"]["expand"] ? config["PDF"]["expand"].as<bool>() : true),
    _scheme(config["Factorisation scheme"].as<std::string>())
  {
  }

  //_________________________________________________________________________
  std::vector<std::complex<double>> NDistributions::Evolve(std::complex<double> const& N, double const& Q) const
  {
    // Call initial-scale distributions
    const auto pdf0 = InitialDistributions(N);

    // If alpha does not run, use a special function that solves the
    // DGLAP exactly.
    if (_aQED.IsFixed())
      {
        const auto EvOp = AlphaFixed(N, Q);
        const auto ns = EvOp.first * pdf0.first;
        const auto sg = EvOp.second * pdf0.second;
        return {sg(0, 0), sg(0, 1), ns};
      }

    // Compute evolution operators
    if (_method == "iterative")
      {
        const auto EvOp = UMatrix(N, Q);
        const auto ns = EvOp.first * pdf0.first;
        const auto sg = EvOp.second * pdf0.second;
        return {sg(0, 0), sg(0, 1), ns};
      }
    else if (_method == "path ordering")
      {
        const auto EvOp = PathOrdering(N, Q);
        const auto ns = EvOp.first * pdf0.first;
        const auto sg = EvOp.second * pdf0.second;
        return {sg(0, 0), sg(0, 1), ns};
      }
    else
      throw std::runtime_error("[NDistributions::Evolve]: Undefined evolution method");
  }

  //_________________________________________________________________________________
  std::pair<std::complex<double>, Matrix<std::complex<double>>> NDistributions::UMatrix(std::complex<double> const& N, double const& Q) const
  {
    // If initial and final scales are equal return the unity
    // operators straight away.
    if (Q == _Qi)
      return std::make_pair(1, Matrix<std::complex<double>> {2, 2, {1, 0, 0, 1}});

    // Number of terms of the iterative solution
    const int nt = 20;

    // b_i coefficients and dglap kernels assignment according to the perturbation order
    double b[nt + 1] = {};
    b[1] = (_ipt > 0 ? _bt1 / _bt0 : 0);

    // Anomalous dimensions
    auto G0 = andim_lo(N, _nl);
    auto G1 = (_ipt > 0 ? andim_nlo(N, _nl) : std::make_pair(0, Matrix<std::complex<double>> {2, 2, {0, 0, 0, 0}}));

    // Computation of r_{1...20} matrices (eq. (2.21) of hep-ph/0408244)
    std::vector<Matrix<std::complex<double>>> Rsg(nt + 1);
    Rsg[0] = G0.second / _bt0;
    Rsg[1] = G1.second / _bt0;
    for (int k = 2; k <= nt; k++)
      Rsg[k] = Matrix<std::complex<double>> {2, 2, {0, 0, 0, 0}};

    std::vector<std::complex<double>> Rns(nt + 1, 0);
    Rns[0] = G0.first / _bt0;
    Rns[1] = G1.first / _bt0;

    for (int k = 1; k <= (_expand ? _ipt : nt); k++)
      for (int l = 1; l <= k; l++)
        Rsg[k] -= Rsg[k-l] * b[l];

    for (int k = 1; k <= nt; k++)
      for (int l = 1; l <= k; l++)
        Rns[k] -= Rns[k-l] * b[l];

    // Computation of r_0 (eq. (2.21) of hep-ph/0408244)
    const std::complex<double> qq0 = Rsg[0](0, 0);
    const std::complex<double> qg0 = Rsg[0](0, 1);
    const std::complex<double> gq0 = Rsg[0](1, 0);
    const std::complex<double> gg0 = Rsg[0](1, 1);

    // Computation of the eigenvalues of r_0 (eq. (2.27) of hep-ph/0408244)
    const std::complex<double> sq = sqrt( pow(qq0 - gg0, 2) + 4. * qg0 * gq0 );
    const std::complex<double> lp = 0.5 * ( qq0 + gg0 + sq );
    const std::complex<double> lm = 0.5 * ( qq0 + gg0 - sq );
    const std::complex<double> ldiff = lm - lp;

    // Computation of the projectors of r_0 (eq. (2.28) of hep-ph/0408244)
    Matrix<std::complex<double>> em{2, 2, {  ( qq0 - lp ) / ldiff,   qg0 / ldiff,   gq0 / ldiff,   ( gg0 - lp ) / ldiff}};
    Matrix<std::complex<double>> ep{2, 2, {- ( qq0 - lm ) / ldiff, - qg0 / ldiff, - gq0 / ldiff, - ( gg0 - lm ) / ldiff}};

    std::vector<Matrix<std::complex<double>>> Usg(nt + 1);
    std::vector<std::complex<double>>         Uns(nt + 1, 0);
    for (int k = 1; k <= nt; k++)
      {
        // Computation of \widetilde{r}_k matrix (eq. (2.25) of hep-ph/0408244)
        Matrix<std::complex<double>> Rtsg = Rsg[k];
        Uns[k] = - Rns[k] / (double) k;
        for (int l = 1; l < k; l++)
          {
            Rtsg += Rsg[l] * Usg[k-l];
            Uns[k] -= Rns[l] * Uns[k-l] / (double) k;
          }

        // Computation of u_k matrix from \widetilde{r}_k (eq. (2.31) of hep-ph/0408244)
        Usg[k] = ( em * Rtsg * em + ep * Rtsg * ep ) / (double) (-k)
                 + em * Rtsg * ep / ( - (double) k - ldiff )
                 + ep * Rtsg * em / ( - (double) k + ldiff );
      }
    // Compute evolution factors
    const double af = _aQED.Evolve(Q);
    const double t  = log(af/_ai);

    // Singlet (eq. (2.29) of hep-ph/0408244)
    const std::complex<double> expm = exp( - lm * t );
    const std::complex<double> expp = exp( - lp * t );
    Matrix<std::complex<double>> Lsg = em * expm + ep * expp;

    // Non-singlet
    std::complex<double> Lns = exp( - Rns[0] * t );

    // Return the LL operators if required
    if (_ipt == 0)
      return std::make_pair(Lns, Lsg);

    const double tmp = (_expand ? af - _ai : log( ( 1 + _bt1 / _bt0 * af ) / ( 1 + _bt1 / _bt0 * _ai ) ) / ( _bt1 / _bt0 ));

    // Non-singlet (eq. (2.34) or (2.35) of hep-ph/0408244)
    const std::complex<double> efnns = Lns * exp( tmp * Uns[1] );

    // Singlet.Calculation of u(a_s) and u(a_0) up to a^20 (which is
    // an approximation of the complete u matrix fot the iterated
    // solution)
    Matrix<std::complex<double>> Usumi{2, 2, {1, 0, 0, 1}};
    Matrix<std::complex<double>> Usumf{2, 2, {1, 0, 0, 1}};
    for (int k = 1; k <= nt; k++)
      {
        Usumi += Usg[k] * pow(_ai, k);
        Usumf += Usg[k] * pow( af, k);
      }

    // Invert Usumi
    Matrix<std::complex<double>> Uinv = Minv(Usumi);

    // Construct singlet evolution matrix
    Matrix<std::complex<double>> efnsg = Usumf * Lsg * Uinv;

    return std::make_pair(efnns, efnsg);
  }

  //_________________________________________________________________________________
  std::pair<std::complex<double>, Matrix<std::complex<double>>> NDistributions::PathOrdering(std::complex<double> const& N, double const& Q) const
  {
    // If initial and final scales are equal return the unity
    // operators straight away.
    if (Q == _Qi)
      return std::make_pair(1, Matrix<std::complex<double>> {2, 2, {1, 0, 0, 1}});

    // Anomalous dimensions
    auto G0 = andim_lo(N, _nl);
    auto G1 = (_ipt > 0 ? andim_nlo(N, _nl) : std::make_pair(0, Matrix<std::complex<double>> {2, 2, {0, 0, 0, 0}}));

    // Scheme-change functions
    const Matrix<std::complex<double>> Jsg = (_scheme == "Delta" ?
                                              Matrix<std::complex<double>> {2, 2, {- fde(N), 0, - fdgm(N), 0}} :
                                              Matrix<std::complex<double>> {2, 2, {0, 0, 0, 0}});
    const std::complex<double> Jns = (_scheme == "Delta" ? - fde(N) : 0.);

    // Initialise singlet evolution operator
    Matrix<std::complex<double>> efnsg{2, 2, {1, 0, 0, 1}};

    // Number of intervals of the path ordered integral
    const int nt = 20;

    const double aff = _aQED.Evolve(Q);
    const double Da  = ( aff - _ai ) / (double) nt;
    double ai = _ai;
    for (int k = 0; k < nt; k++)
      {
        const double af = ai + Da;

        // Singlet
        Matrix<std::complex<double>> intsgi = Integrandsg(G0.second, G1.second, Jsg, ai);
        Matrix<std::complex<double>> intsgf = Integrandsg(G0.second, G1.second, Jsg, af);
        const Matrix<std::complex<double>> spsg = ( intsgi + intsgf ) * Da / 2.;

        // Formulae (15)-(19) on http://mathworld.wolfram.com/matrixexponential.html
        const std::complex<double> delta = sqrt( pow(spsg(0, 0) - spsg(1, 1), 2) + 4. * spsg(0,1) * spsg(1, 0) );
        const std::complex<double> sh    = exp( ( spsg(0, 0) + spsg(1, 1) ) / 2. ) * ( exp( delta / 2. ) - exp( - delta / 2. ) ) / 2.;
        const std::complex<double> ch    = exp( ( spsg(0, 0) + spsg(1, 1) ) / 2. ) * ( exp( delta / 2. ) + exp( - delta / 2. ) ) / 2.;

        const Matrix<std::complex<double>> defnsg{2, 2,
          {
            ch + sh * ( spsg(0, 0) - spsg(1, 1) ) / delta,
            2. * spsg(0,1) * sh / delta,
            2. * spsg(1,0) * sh / delta,
            ch - sh * ( spsg(0, 0) - spsg(1, 1) ) / delta
          }};

        // Apply step evolution
        efnsg = defnsg * efnsg;

        // Increment step
        ai += Da;
      }

    // Non-singlet
    const std::complex<double> intnsi = Integrandns(G0.first, G1.first, Jns, _ai);
    const std::complex<double> intnsf = Integrandns(G0.first, G1.first, Jns, aff);
    std::complex<double> spns = Da * ( intnsi + intnsf ) / 2.;

    double ak = _ai;
    for (int k = 0; k < nt - 1; k++)
      {
        ak += Da;
        spns += Da * Integrandns(G0.first, G1.first, Jns, ak);
      }
    const std::complex<double> efnns = exp(spns);

    return std::make_pair(efnns, efnsg);
  }

  //_________________________________________________________________________________
  std::pair<std::complex<double>, Matrix<std::complex<double>>> NDistributions::AlphaFixed(std::complex<double> const& N, double const& Q) const
  {
    // DGLAP evolution with alpha fixed in not available in the Delta-scheme
    if (_scheme == "Delta")
      throw std::runtime_error("[NDistributions::AlphaFixed]: DGLAP evolution with alpha fixed in the Delta-scheme unavailable");

    // If initial and final scales are equal return the unity
    // operators straight away.
    if (Q == _Qi)
      return std::make_pair(1, Matrix<std::complex<double>> {2, 2, {1, 0, 0, 1}});

    // Log of the scales
    const double ln = 2 * log(Q / _Qi);

    // Anomalous dimensions
    auto G0 = andim_lo(N, _nl);
    auto G1 = (_ipt > 0 ? andim_nlo(N, _nl) : std::make_pair(0, Matrix<std::complex<double>> {2, 2, {0, 0, 0, 0}}));

    // Singlet
    const Matrix<std::complex<double>> spsg = _ai * ( G0.second + _ai * G1.second ) * ln;

    // Formulae (15)-(19) on http://mathworld.wolfram.com/matrixexponential.html
    const std::complex<double> delta = sqrt( pow(spsg(0, 0) - spsg(1, 1), 2) + 4. * spsg(0,1) * spsg(1, 0) );
    const std::complex<double> sh    = exp( ( spsg(0, 0) + spsg(1, 1) ) / 2. ) * ( exp( delta / 2. ) - exp( - delta / 2. ) ) / 2.;
    const std::complex<double> ch    = exp( ( spsg(0, 0) + spsg(1, 1) ) / 2. ) * ( exp( delta / 2. ) + exp( - delta / 2. ) ) / 2.;

    const Matrix<std::complex<double>> efnsg{2, 2,
      {
        ch + sh * ( spsg(0, 0) - spsg(1, 1) ) / delta,
        2. * spsg(0,1) * sh / delta,
        2. * spsg(1,0) * sh / delta,
        ch - sh * ( spsg(0, 0) - spsg(1, 1) ) / delta
      }};

    // Non-singlet
    const std::complex<double> efnns = exp( _ai * ( G0.first + _ai * G1.first ) * ln );

    return std::make_pair(efnns, efnsg);
  }

  //_________________________________________________________________________________
  std::pair<std::complex<double>, Matrix<std::complex<double>>> NDistributions::InitialDistributions(std::complex<double> const& N) const
  {
    const std::complex<double> Nm = N - 1.;
    const std::complex<double> N1 = N + 1.;
    const std::complex<double> s1 = emc + psi(N1);

    std::complex<double> deN  = 0.;
    std::complex<double> dgmN = 0.;
    if (_scheme == "MSbar")
      {
        deN  = fde(N);
        dgmN = fdgm(N);
      }
    const double logfact = 2 * log(_Qi / me);

    const std::complex<double> sg = 1. + (_ipt > 0 ? _ai * ( logfact * 2. * ( 1. / N - 1. / N1 - 2. * s1 + 3. / 2. ) + deN ) : 0);
    const std::complex<double> gm = (_ipt > 0 ? _ai * ( logfact * 2. * ( 2. / Nm - 2. / N + 1. / N1 ) + dgmN ) : 0);
    const std::complex<double> ns = sg;

    return std::make_pair(ns, Matrix<std::complex<double>> {2, 1, {sg, gm}});
  }

  //_________________________________________________________________________________
  std::complex<double> NDistributions::fde(std::complex<double> const& N) const
  {
    const std::complex<double> N1 = N + 1.;
    const std::complex<double> s1 = emc + psi(N1);
    const std::complex<double> s2 = zeta2 - dpsi(N1, 1);
    return 2. * ( - ( 1. / N - 1. / N1 - 2. * s1 + 3. / 2. )
                  - 2. * ( pow(s1, 2) + s2 - s1 / N + s1 / N1 + 1. / pow(N1, 2) - 7. / 4. ) );
  }

  //_________________________________________________________________________________
  std::complex<double> NDistributions::fdgm(std::complex<double> const& N) const
  {
    const std::complex<double> Nm = N - 1.;
    const std::complex<double> N1 = N + 1.;
    return 2. * ( - ( 2. / Nm - 2. / N + 1. / N1 ) - 2. * ( - 2. / pow(Nm, 2) + 2. / pow(N, 2) - 1. / pow(N1, 2) ) );
  }

  //_________________________________________________________________________________
  Matrix<std::complex<double>> NDistributions::Integrandsg(Matrix<std::complex<double>> const& G0,
                                                           Matrix<std::complex<double>> const& G1,
                                                           Matrix<std::complex<double>> const& J,
                                                           double const& a) const
  {
    const Matrix<std::complex<double>> Jd    = Matrix<std::complex<double>> {2, 2, {1, 0, 0, 1}} + J * a;
    const Matrix<std::complex<double>> Jdinv = Minv(Jd);
    Matrix<std::complex<double>> P;
    if (_expand)
      P = (-1.) * ( G0 + a * ( G1 - G0 * (_ipt == 0 ? 0 : _bt1 / _bt0) ) ) / a / _bt0;
    else
      P = ( G0 + G1 * a ) / a / ( - _bt0 - a * (_ipt == 0 ? 0 : _bt1) );
    return Jd * P * Jdinv;
  }

  //_________________________________________________________________________________
  std::complex<double> NDistributions::Integrandns(std::complex<double> const& G0,
                                                   std::complex<double> const& G1,
                                                   std::complex<double> const& J,
                                                   double const& a) const
  {
    const std::complex<double> Jd = 1. + a * J;
    if (_expand)
      return J / Jd - ( G0 + a * ( G1 - G0 * (_ipt == 0 ? 0 : _bt1 / _bt0) ) ) / a / _bt0;
    else
      return J / Jd + ( G0 + a * G1 ) / a / ( - _bt0 - a * (_ipt == 0 ? 0 : _bt1) );
  }

  //_________________________________________________________________________________
  Matrix<std::complex<double>> NDistributions::Minv(Matrix<std::complex<double>> const& M) const
  {
    const std::complex<double> det = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
    return Matrix<std::complex<double>> {2, 2, {M(1, 1) / det, - M(0, 1) / det, - M(1, 0) / det, M(0, 0) / det}};
  }

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, Matrix<std::complex<double>> const& m)
  {
    os << "\n";
    for (int i = 0; i < m.GetLines(); i++)
      {
        os << "| ";
        for (int j = 0; j < m.GetColumns(); j++)
          {
            os << m.GetElement(i, j) << "\t";
          }
        os << "|\n";
      }
    os << "\n";
    return os;
  }
}
