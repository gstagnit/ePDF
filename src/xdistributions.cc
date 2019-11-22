//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "ePDF/xdistributions.h"

#include <math.h>

namespace ePDF
{
  /**
    * Parameters of the dgauss integration with 4, 8, 16, 32 and 64-point integration.
    */
  const std::array<std::vector<double>,5> gq_x =
  {
    {
      std::vector<double>{0.3399810435848562, 0.8611363115940525},
      std::vector<double>{0.1834346424956498, 0.5255324099163289, 0.7966664774136267, 0.9602898564975362},
      std::vector<double>{
        0.0950125098376374, 0.2816035507792589, 0.4580167776572273, 0.6178762444026437,
        0.7554044083550030, 0.8656312023878317, 0.9445750230732325, 0.9894009349916499
      },
      std::vector<double>{
        0.0483076656877383, 0.1444719615827964, 0.2392873622521370, 0.3318686022821276,
        0.4213512761306353, 0.5068999089322293, 0.5877157572407623, 0.6630442669302152,
        0.7321821187402896, 0.7944837959679424, 0.8493676137325699, 0.8963211557660521,
        0.9349060759377396, 0.9647622555875064, 0.9856115115452683, 0.9972638618494815
      },
      std::vector<double>{
        0.0243502926634244, 0.0729931217877990, 0.1214628192961205, 0.1696444204239928,
        0.2174236437400070, 0.2646871622087674, 0.3113228719902109, 0.3572201583376681,
        0.4022701579639916, 0.4463660172534640, 0.4894031457070529, 0.5312794640198945,
        0.5718956462026340, 0.6111553551723932, 0.6489654712546573, 0.6852363130542332,
        0.7198818501716108, 0.7528199072605318, 0.7839723589433414, 0.8132653151227975,
        0.8406292962525803, 0.8659993981540928, 0.8893154459951141, 0.9105221370785028,
        0.9295691721319395, 0.9464113748584028, 0.9610087996520537, 0.9733268277899109,
        0.9833362538846259, 0.9910133714767443, 0.9963401167719552, 0.9993050417357721
      }
    }
  };
  const std::array<std::vector<double>,5> gq_w =
  {
    {
      std::vector<double>{0.6521451548625461, 0.3478548451374538},
      std::vector<double>{0.3626837833783619, 0.3137066458778872, 0.2223810344533744, 0.1012285362903762},
      std::vector<double>{
        0.1894506104550684, 0.1826034150449235, 0.1691565193950025, 0.1495959888165767,
        0.1246289712555338, 0.0951585116824927, 0.0622535239386478, 0.0271524594117540
      },
      std::vector<double>{
        0.0965400885147278, 0.0956387200792748, 0.0938443990808045, 0.0911738786957638,
        0.0876520930044038, 0.0833119242269467, 0.0781938957870703, 0.0723457941088485,
        0.0658222227763618, 0.0586840934785355, 0.0509980592623761, 0.0428358980222266,
        0.0342738629130214, 0.0253920653092620, 0.0162743947309056, 0.0070186100094700
      },
      std::vector<double>{
        0.0486909570091397, 0.0485754674415034, 0.0483447622348029, 0.0479993885964583,
        0.0475401657148303, 0.0469681828162100, 0.0462847965813144, 0.0454916279274181,
        0.0445905581637565, 0.0435837245293234, 0.0424735151236535, 0.0412625632426235,
        0.0399537411327203, 0.0385501531786156, 0.0370551285402400, 0.0354722132568823,
        0.0338051618371416, 0.0320579283548515, 0.0302346570724024, 0.0283396726142594,
        0.0263774697150546, 0.0243527025687108, 0.0222701738083832, 0.0201348231535302,
        0.0179517157756973, 0.0157260304760247, 0.0134630478967186, 0.0111681394601311,
        0.0088467598263639, 0.0065044579689783, 0.0041470332605624, 0.0017832807216964
      }
    }
  };

  /*
   * Function to be fed to GSL for the integration along the Talbot
   * path. The integral has to be between 0 and 1.
   */
  double talbotGSL(double y, void *params)
  {
    struct evol_params *p = (struct evol_params*) params;
    const double t = - log(p->x);
    const int    m = 30;
    const double r = 2 * m / 5 / t;
    const double theta =  M_PI * y;
    const double sigma = theta + ( theta / tan(theta) - 1 ) / tan(theta);
    const std::complex<double> tanm1(1. / tan(theta), 1.);
    const std::complex<double> Nt = r * theta * tanm1 + 1.;
    const std::complex<double> imsigma(1.0, sigma);
    return r * std::real(exp(t * Nt) * imsigma * p->Ndist->Evolve(Nt, p->Q)[p->id]);
  }

  /*
   * Function to be fed to GSL for the integration along the Talbot
   * path. The integral has to be between 0 and 1.
   */
  double straightGSL(double y, void *params)
  {
    struct evol_params *p = (struct evol_params*) params;
    const double t   = - log(p->x);
    const double z   = - log(y) / t;
    const double jac = - 1 / y / t;
    const double c   = 1.9;
    const std::complex<double> ph = exp(3. * std::complex<double>(0, 1) * M_PI / 4.);
    const std::complex<double> N  = c + z * ph;
    const std::complex<double> fact = - p->x * jac * ph * exp( t * N ) / M_PI;
    return std::imag(fact * p->Ndist->Evolve(N, p->Q)[p->id]);
  }

  //_________________________________________________________________________________
  xDistributions::xDistributions(YAML::Node const& config):
    _contour((config["Mellin inverse"]["contour"] ? config["Mellin inverse"]["contour"].as<std::string>() : "talbot")),
    _integrator(config["Mellin inverse"]["integrator"].as<std::string>()),
    _eps(config["Mellin inverse"]["accuracy"].as<double>()),
    _p{0, 0, 0, std::shared_ptr<NDistributions>(new NDistributions{config})},
    _w(gsl_integration_workspace_alloc(10000))
  {
  }

  //_________________________________________________________________________________
  xDistributions::~xDistributions()
  {
    gsl_integration_workspace_free(_w);
  }

  //_________________________________________________________________________
  std::vector<double> xDistributions::gaussGSL(double const& a, double const& b)
  {
    // Allocate GSL function according to the path chosen.
    gsl_function F;
    F.params   = &_p;
    if (_contour == "talbot")
      F.function = &talbotGSL;
    else if (_contour == "straight")
      F.function = &straightGSL;
    else
      throw std::runtime_error("[xDistribution::gaussGSL]: unknown contour");

    double result, error;
    std::vector<double> v(3);
    for (int id = 0; id < 3; id++)
      {
        _p.id = id;
        gsl_integration_qag(&F, a, b, 0, _eps, 6, 10000, _w, &result, &error);
        v[id] = _p.x * result;
      }
    return v;
  }

  //_________________________________________________________________________
  std::vector<double> xDistributions::Evolve(double const& x, double const& Q)
  {
    // Kinematics
    _p.x  = x;
    _p.Q  = Q;
    if (_integrator == "gauss")
      return gauss(0, 1);
    else if (_integrator == "trapezoid")
      return trapezoid(0, 1);
    else if (_integrator == "gauss GSL")
      return gaussGSL(0, 1);
    else
      throw std::runtime_error("[xDistribution::Evolve]: unknown integrator");
  }

  //_________________________________________________________________________
  std::vector<double> xDistributions::integrand(double const& y) const
  {
    if (_contour == "talbot")
      return talbot(y);
    else if (_contour == "straight")
      return straight(y);
    else
      throw std::runtime_error("[xDistribution::integrand]: unknown contour");
  }

  //_________________________________________________________________________
  std::vector<double> xDistributions::talbot(double const& y) const
  {
    const double t = - log(_p.x);
    const int    m = 43;
    const double r = 2 * m / 5 / t;
    const double theta = M_PI * y;
    const double sigma = theta + ( theta / tan(theta) - 1 ) / tan(theta);
    const std::complex<double> tanm1(1. / tan(theta), 1.);
    const std::complex<double> Nt = r * theta * tanm1 + 1.;
    const std::complex<double> imsigma(1.0, sigma);
    const std::complex<double> factor = _p.x * r * exp(t * Nt) * imsigma;
    const std::vector<std::complex<double>> fN = _p.Ndist->Evolve(Nt, _p.Q);
    std::vector<double> res(3);
    for (int j = 0; j < 3; j++)
      res[j] = std::real(factor * fN[j]);
    return res;
  }

  //_________________________________________________________________________
  std::vector<double> xDistributions::straight(double const& y) const
  {
    /*
        const double t = - log(_p.x);
        const double z = - log(y) / t;
        const double jac = 1 / y / t;
        const double c   = 1;
        const std::complex<double> ci(0, 1);
        const std::complex<double> N = - z * ( 1. - ci ) + c;
        const std::complex<double> fact = - jac * ( 1. - ci ) * exp( N * t ) / M_PI / ci;
        const std::vector<std::complex<double>> fN = _p.Ndist->Evolve(N, _p.Q);
        std::vector<double> res(3);
        for (int j = 0; j < 3; j++)
          res[j] = std::real(fact * fN[j]);
        return res;
    */
    const double t   = - log(_p.x);
    const double z   = - log(y) / t;
    const double jac = - 1 / y / t;
    const double c   = 1.9;
    const std::complex<double> p = exp(3. * std::complex<double>(0, 1) * M_PI / 4.);
    const std::complex<double> N = c + z * p;
    const std::complex<double> fact = - _p.x * jac * p * exp( t * N ) / M_PI;
    const std::vector<std::complex<double>> fN = _p.Ndist->Evolve(N, _p.Q);
    std::vector<double> res(3);
    for (int j = 0; j < 3; j++)
      res[j] = std::imag(fact * fN[j]);
    return res;
  }

  //_________________________________________________________________________
  std::vector<double> xDistributions::trapezoid(double const& a, double const& b) const
  {
    const int n = 30;
    const double step = 2 * ( b - a ) / ( 2 * n + 1 );
    std::vector<double> res(3);
    for (int i = 0; i < n; i++)
      {
        const double y = b - step * ( i + 1 );
        const std::vector<double> t = integrand(y);
        for (int j = 0; j < 3; j++)
          res[j] += t[j];
      }
    for (int j = 0; j < 3; j++)
      res[j] *= step;
    return res;
  }

  //_________________________________________________________________________
  std::vector<double> xDistributions::gauss(double const& a, double const& b) const
  {
    const double delta = 1e-25 * std::abs(a - b);
    std::vector<double> dgauss(3, 0.);
    double aa = a;

goto5:
    double y = b - aa;
    if (std::abs(y) <= delta)
      return dgauss;

goto2:
    double bb = aa + y;
    double c1 = 0.5 * ( aa + bb );
    double c2 = c1 - aa;

    std::vector<double> s8(3, 0.);
    for (int i = 0; i < 4; i++)
      {
        double u = gq_x[1][i] * c2;
        const std::vector<double> U = integrand(c1+u);
        const std::vector<double> D = integrand(c1-u);
        for (int j = 0; j < 3; j++)
          s8[j] += gq_w[1][i] * ( U[j] + D[j] );
      }

    std::vector<double> s16(3, 0.);
    for (int i = 0; i < 8; i++)
      {
        double u = gq_x[2][i] * c2;
        const std::vector<double> U = integrand(c1+u);
        const std::vector<double> D = integrand(c1-u);
        for (int j = 0; j < 3; j++)
          s16[j] += gq_w[2][i] * ( U[j] + D[j] );
      }

    for (int j = 0; j < 3; j++)
      {
        s8[j]  *= c2;
        s16[j] *= c2;
        if (std::abs(s16[j] - s8[j]) > _eps * ( 1 + std::abs(s16[j]) ))
          goto goto4;
        dgauss[j] += s16[j];
      }

    aa = bb;
    goto goto5;

goto4:
    y = 0.5 * y;
    if (std::abs(y) > delta)
      goto goto2;

    throw std::runtime_error("[xDistribution::gauss]: too high accuracy required");

    return dgauss;
  }

  //_________________________________________________________________________
  std::vector<std::complex<double>> xDistributions::Moments(std::complex<double> const& N, double const& Q) const
  {
    return _p.Ndist->Evolve(N, Q);
  }
}
