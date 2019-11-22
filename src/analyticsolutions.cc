//
// Author: Giovanni Stagnitto
//

#include "ePDF/analyticsolutions.h"
#include "ePDF/constants.h"
#include "ePDF/specialfunctions.h"
#include "ePDF/numericintegrals.h"

#include <math.h>
#include <iostream>

namespace ePDF
{
  //_________________________________________________________________________________
  AnalyticSolutions::AnalyticSolutions(YAML::Node        const& config,
                                       int               const& orderA,
                                       std::vector<bool> const orderR):
    _aQED(AlphaQED{config}),
    _Qi(config["Initial scale"].as<double>()),
    _nl(config["NL"].as<double>()),
    _b0(_aQED.Beta0QED()/(-4.*M_PI)),
    _b1(_aQED.Beta1QED()/(-16.*Pi2)),
    _a0twopi(_aQED.Evolve(_Qi)*2),
    _L0(log( _Qi * _Qi / ( me * me ) )),
    _orderR(orderR),
    _orderA(orderA),
    _x0photon(1),
    _x1photon(7),
    _pphoton(4),
    _numint(true),
    _rechatON(true)
  {
  }

  //_________________________________________________________________________________
  AnalyticSolutions::AnalyticSolutions(YAML::Node const& config):
    AnalyticSolutions{config,
                      ptmap.at(config["Perturbative order"].as<std::string>()),
                      ptmap.at(config["Perturbative order"].as<std::string>()) ? orderRNLL : orderRLL}
  {
  }

  //_________________________________________________________________________________
  void AnalyticSolutions::SetPhotonMatching(std::vector<double> const& vec)
  {
    _x0photon = vec[0];
    _x1photon = vec[1];
    _pphoton  = vec[2];
  }

  void AnalyticSolutions::SetNumInt(bool const& numint)
  {
    _numint = numint;
  }

  void AnalyticSolutions::SetRecHat(bool const& rechatON)
  {
    _rechatON = rechatON;
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::GetSolution(double const& x, double const& Q, int const& id, std::string const& term)
  {
    Warmup(Q);
    if (term == "A")
      {
        return AsySolution(x, id);
      }
    else if (term == "R")
      {
        return RecSolution(x, id);
      }
    else if (term == "ABAR")
      {
        return AsyBarSolution(x, id);
      }
    else if (term == "RBAR")
      {
        return RecBarSolution(x, id);
      }
    else
      {
        throw std::runtime_error("[AnalyticSolutions::GetSolution]: unknown string");
      }

    return 0.0;
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::TestPhoton(double const& x, double const& Q, std::string test)
  {
    Warmup(Q);
    if (test == "SIMPLIFIED")
      {
        return AsyPhotonSIMPLIFIED(x);
      }
    else
      {
        return AsyPhoton(x);
      }
  }

  //_________________________________________________________________________________
  double Matching(double const& z, double const& x0, double const& x1, double const& p)
  {
    const double x = - log10(1-z);
    if (x < x0)
      return 0.0;
    if (x > x1)
      return 1.0;

    const double xm = ( x - x0 ) / ( x1 - x0 );

    return pow( xm, p ) / ( pow( xm, p ) + pow( ( 1 - xm ), p ) );
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::Evolve(double const& x, double const& Q, int const& id)
  {
    Warmup(Q);

    if (id < 0 || id > 2)
      throw std::runtime_error("[AnalyticSolutions::Evolve]: unknown flavour");

    const double asy    = AsySolution(x, id);
    const double asybar = AsyBarSolution(x, id);
    const double rec    = RecSolution(x, id);
    const double recbar = RecBarSolution(x, id);

    double res(0.0);
    switch (id)
      {
      case 0:
        res = asy + rec - asybar;
        break;
      case 1:
        // TODO : matching with subt = asy or subt = rec ? I would say subt = rec
        res = rec + Matching(x, _x0photon, _x1photon, _pphoton) * ( asy - recbar );
        break;
      case 2:
        res = asy + rec - asybar;
        break;
      }
    return res;
  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::Evolve(double const& x, double const& Q)
  {
    return{Evolve(x, Q, 0), Evolve(x, Q, 1), Evolve(x, Q, 2)};
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::RecSeries(std::vector<double> const& series)  const
  {
    double p  = ( _aQED.IsFixed() ) ? _eta0/2.0 : _t;
    double p2 = p * p;
    double p3 = p2 * p;

    return p * series[0] + (p2/2.) * series[1] + (p3/6.) * series[2]
           + _atwopi * (series[3] + p * series[4] + (p2/2.) * series[5]);
  }

  //_________________________________________________________________________________
  void AnalyticSolutions::Warmup(double const& Q)
  {
    _atwopi = _aQED.Evolve(Q)*2.;
    _t      = 1./(2.*_b0*M_PI)*log(_atwopi/_a0twopi);
    _eta0   = _aQED.Evolve(Q)*4. * log( (Q*Q) / (_Qi*_Qi) );
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::AsySolution(double const& x, int const& id) const
  {
    switch (id)
      {
      case 0:
        return ( _aQED.IsFixed() ) ? AsyEleAF(x) : AsyEleAR(x);
      case 1:
        return AsyPhoton(x);
      case 2:
        return ( _aQED.IsFixed() ) ? AsyEleAF(x) : AsyEleAR(x);
      default:
        throw std::runtime_error("[AnalyticSolutions::AsySolution]: unknown flavour");
      }

    return 0.0;
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::AsyBarSolution(double const& x, int const& id) const
  {
    std::vector<double> asybar(5, 0.0);
    switch (id)
      {
      case 0 :
        asybar = ( _aQED.IsFixed() ) ? AsyEleBarAF(x) : AsyEleBarAR(x);
        break;
      case 1 :
        asybar = AsyPhotonBar(x);
        break;
      case 2 :
        asybar = ( _aQED.IsFixed() ) ? AsyEleBarAF(x) : AsyEleBarAR(x);
        break;
      }
    return RecSeries(asybar);
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::RecBarSolution(double const& x, int const& id) const
  {
    std::vector<double> recbar(5, 0.0);

    // Recursive result
    switch (id)
      {
      case 0:
        recbar = RecBarS(x);
        break;
      case 1:
        recbar = RecBarG(x);
        break;
      case 2:
        recbar = RecBarNS(x);
        break;
      }

    return RecSeries(recbar);
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::RecSolution(double const& x, int const& id) const
  {
    double recbar(0.0), rechat(0.0), recnum(0.0);

    if (_numint && _orderR[5] == true)
      {
        NumericIntegrals numint(id, _nl);
        double tmp = numint.integrate(x);

        // numeric contributions only at order alpha*t^2
        std::vector<double> recnumS = { 0.0, 0.0, 0.0, 0.0, 0.0, tmp };
        recnum = RecSeries(recnumS);
      }

    std::vector<double> rechatseries(5, 0.0);
    if (_rechatON)
      {
        switch (id)
          {
          case 0 :
            rechatseries = RecHatS(x);
            break;
          case 1 :
            rechatseries = RecHatG(x);
            break;
          case 2 :
            rechatseries = RecHatNS(x);
            break;
          }
      }
    rechat = RecSeries(rechatseries);
    recbar = RecBarSolution(x, id);

    return rechat + recbar + recnum;
  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::RecBarNS(double const& z) const
  {
    double JLL1(0.0), JLL2(0.0), JLL3(0.0), JNLL0(0.0), JNLL1(0.0), JNLL2(0.0);

    const double log1mz = log(1-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = pow(log1mz,3);

    double B0(0.0), B02(0.0), B1overB0(0.0);
    if( ! _aQED.IsFixed() )
      {
        B0 = _b0;
        B1overB0 = _b1/_b0;
        B02 = _b0*_b0;
      }

    if (_orderR[0])
      JLL1 = -2 + 2/(1 - z);
    // -----------------------------------------------------------------------
    if (_orderR[1])
      JLL2 = -2 + 6/(1 - z) - 8*log1mz + (8*log1mz)/(1 - z);
    // -----------------------------------------------------------------------
    if (_orderR[2])
      JLL3 = 4.5 + 4*Pi2 + (13.5 - 4*Pi2)/(1 - z) - 12*log1mz + (36*log1mz)/(1 - z)
             - 24*log1mz2 + (24*log1mz2)/(1 - z);
    // -----------------------------------------------------------------------
    if (_orderR[3])
      JNLL0 = 2-2*_L0+(-2+2*_L0)/(1-z)+4*log1mz-(4*log1mz)/(1-z);
    // -----------------------------------------------------------------------
    if (_orderR[4])
      JNLL1 = -2-2*_L0+(32*_nl)/9.+4*B1overB0*M_PI-(4*Pi2)/3.+1/(1-z)+(6*_L0)/(1-z)-(20*_nl)/(9.*(1-z))-(4*B1overB0*M_PI)/(1-z)+(4*Pi2)/(3.*(1-z))+10*log1mz-8*_L0*log1mz-(14*log1mz)/(1-z)+(8*_L0*log1mz)/(1-z)+12*log1mz2-(12*log1mz2)/(1-z)+B0*(-4*M_PI+4*_L0*M_PI+(4*M_PI)/(1-z)-(4*_L0*M_PI)/(1-z)-8*M_PI*log1mz+(8*M_PI*log1mz)/(1-z));
    // -----------------------------------------------------------------------
    if (_orderR[5])
      JNLL2 = -4+(9*_L0)/2.+(22*_nl)/9.+8*B1overB0*M_PI-(10*Pi2)/3.+4*_L0*Pi2+(8*_nl*Pi2)/9.+9/(1-z)+(27*_L0)/(2.*(1-z))-(22*_nl)/(3.*(1-z))-(24*B1overB0*M_PI)/(1-z)+(6*Pi2)/(1-z)-(4*_L0*Pi2)/(1-z)-(8*_nl*Pi2)/(9.*(1-z))-7*log1mz-12*_L0*log1mz+(208*_nl*log1mz)/9.+32*B1overB0*M_PI*log1mz-(40*Pi2*log1mz)/3.-(17*log1mz)/(1-z)+(36*_L0*log1mz)/(1-z)-(160*_nl*log1mz)/(9.*(1-z))-(32*B1overB0*M_PI*log1mz)/(1-z)+(40*Pi2*log1mz)/(3.*(1-z))+36*log1mz2-24*_L0*log1mz2-(60*log1mz2)/(1-z)+(24*_L0*log1mz2)/(1-z)+32*log1mz3-(32*log1mz3)/(1-z)+B02*(8*Pi2-8*_L0*Pi2-(8*Pi2)/(1-z)+(8*_L0*Pi2)/(1-z)+16*Pi2*log1mz-(16*Pi2*log1mz)/(1-z))+B0*(14*M_PI+8*_L0*M_PI-(64*_nl*M_PI)/9.-8*B1overB0*Pi2+(16*pow(M_PI,3))/3.-(24*_L0*M_PI)/(1-z)+(40*_nl*M_PI)/(9.*(1-z))+(8*B1overB0*Pi2)/(1-z)+(-4*M_PI-(16*pow(M_PI,3))/3.)/(1-z)-32*M_PI*log1mz+32*_L0*M_PI*log1mz+(56*M_PI*log1mz)/(1-z)-(32*_L0*M_PI*log1mz)/(1-z)-48*M_PI*log1mz2+(48*M_PI*log1mz2)/(1-z))+40*zeta3-(40*zeta3)/(1-z);

    // std::cout << std:: endl << "test: " << JLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JLL2 << std::endl;
    // std::cout << std:: endl << "test: " << JLL3 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL0 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL2 << std::endl;

    return {JLL1, JLL2, JLL3, JNLL0, JNLL1, JNLL2};
  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::RecHatNS(double const& z) const
  {
    double JLL1(0.0), JLL2(0.0), JLL3(0.0), JNLL0(0.0), JNLL1(0.0), JNLL2(0.0);

    const double z2 = z*z;
    const double z3 = pow(z,3);
    const double z4 = pow(z,4);

    const double log1mz = log(1-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = pow(log1mz,3);

    const double logz = log(z);
    const double logz2 = logz*logz;
    const double logz3 = pow(logz,3);

    const double log1pz = log(1+z);
    const double log1pz2 = log1pz*log1pz;
    const double log1pz3 = pow(log1pz,3);

    const double log2 = log(2.);

    double B0(0.0), B02(0.0), B1overB0(0.0);
    if( ! _aQED.IsFixed() )
      {
        B0 = _b0;
        B1overB0 = _b1/_b0;
        B02 = _b0*_b0;
      }

    if (_orderR[0])
      JLL1 = 1 - z;
    // -----------------------------------------------------------------------
    if (_orderR[1])
      JLL2 = ((1 - z)*(3 + z) - 4*pow(-1 + z,2)*log1mz + log(z) + 3*pow(z,2)*log(z))/(-1 + z);
    // -----------------------------------------------------------------------
    if (_orderR[2])
      JLL3 = ((-1 + z)*(4*Pi2*(-1 + 3*z) - 3*(19 + 5*z)) - 48*pow(-1 + z,2)*pow(log1mz,2) - 2*log(z)*(-3 - 9*pow(z,2) + log(z) + 7*pow(z,2)*log(z)) + 24*log1mz*(3 - pow(z,2) + 2*(1 + pow(z,2))*log(z)) - 48*z*(log(2 - 2*z) - log(2*z)) - 24*(-1 + pow(z,2))*polylog(2,z))/(4.*(-1 + z));
    // -----------------------------------------------------------------------
    if (_orderR[3])
      JNLL0 = (1-z)*(-1+_L0-2*log1mz);
    // -----------------------------------------------------------------------
    if (_orderR[4])
      JNLL1 = (135+54*_L0+68*_nl-36*B0*M_PI+36*B1overB0*M_PI+36*B0*_L0*M_PI+18*Pi2-153*z+18*_L0*z-44*_nl*z+36*B0*M_PI*z-36*B1overB0*M_PI*z-36*B0*_L0*M_PI*z-18*Pi2*z-135*z2-54*_L0*z2-68*_nl*z2+36*B0*M_PI*z2-36*B1overB0*M_PI*z2-36*B0*_L0*M_PI*z2+6*Pi2*z2+153*z3-18*_L0*z3+44*_nl*z3-36*B0*M_PI*z3+36*B1overB0*M_PI*z3+36*B0*_L0*M_PI*z3-6*Pi2*z3-18*log1mz-72*_L0*log1mz-72*B0*M_PI*log1mz-54*z*log1mz+72*_L0*z*log1mz+72*B0*M_PI*z*log1mz+18*z2*log1mz+72*_L0*z2*log1mz+72*B0*M_PI*z2*log1mz+54*z3*log1mz-72*_L0*z3*log1mz-72*B0*M_PI*z3*log1mz+108*log1mz2-108*z*log1mz2-108*z2*log1mz2+108*z3*log1mz2+72*logz+18*_L0*logz+12*_nl*logz+72*z*logz+18*_L0*z*logz+12*_nl*z*logz-90*z2*logz+54*_L0*z2*logz+12*_nl*z2*logz-90*z3*logz+54*_L0*z3*logz+12*_nl*z3*logz-72*z2*log1mz*logz-72*z3*log1mz*logz+27*logz2-9*z*logz2+9*z2*logz2-27*z3*logz2+144*logz*log1pz-144*z*logz*log1pz-72*z2*logz*log1pz+72*z3*logz*log1pz-108*log1pz2+108*z*log1pz2-36*(-1+z)*pow(1+z,2)*polylog(2,1-z)+72*(2-2*z-z2+z3)*polylog(2,-z)-216*polylog(2,1/(1+z))+216*z*polylog(2,1/(1+z)))/(18.*(-1+z2));
    // -----------------------------------------------------------------------
    if (_orderR[5])
      JNLL2 = (24*Pi2*(-1+z)*z+4*_nl*(63+4*Pi2*(-1+z)-19*z)*(-1+z)*z*(1+z)+144*B1overB0*M_PI*(-1+z)*z*(1+z)*(3+z)-12*Pi2*(-1+z)*z2*(3+17*z)-18*(-1+z)*z*(1+z)*(-73+57*z)+9*_L0*(-1+z)*z*(1+z)*(8*Pi2*(-1+z)-3*(19+5*z))+24*log2*(-4*Pi2*z*(-1+2*z+z3)+2*z*(-11+z*(16+z+4*z2))*pow(log2,2)+(-1+z)*(1+z)*(3+z*(7+9*z))*log(8))-144*B02*Pi2*pow(-1+z,2)*z*(1+z)*(-1+_L0-2*log1mz)+576*B1overB0*M_PI*pow(-1+z,2)*z*(1+z)*log1mz-2*(-1+z)*z*(108*_L0*(1+z)*(3+z)+32*_nl*(11+(3-8*z)*z)+3*(81-27*z*(4+7*z)+4*Pi2*(-6+z+9*z2)+120*pow(log2,2)))*log1mz+504*(-1+z)*z*(1+z)*log1mz2-432*_L0*pow(-1+z,2)*z*(1+z)*log1mz2+360*(-1+z)*z2*(1+z)*log1mz2+576*pow(-1+z,2)*z*(1+z)*log1mz3-540*z*logz+54*_L0*z*logz-152*_nl*z*logz-144*B1overB0*M_PI*z*logz-84*Pi2*z*logz+612*z2*logz+486*_L0*z2*logz+136*_nl*z2*logz-144*B1overB0*M_PI*z2*logz+132*Pi2*z2*logz+468*z3*logz+594*_L0*z3*logz-24*_nl*z3*logz-432*B1overB0*M_PI*z3*logz+84*Pi2*z3*logz-684*z4*logz+162*_L0*z4*logz-312*_nl*z4*logz-432*B1overB0*M_PI*z4*logz-36*Pi2*z4*logz-432*log2*logz-288*z*log2*logz+576*z2*log2*logz+288*z3*log2*logz-144*z4*log2*logz+144*z*pow(log2,2)*logz-144*z2*pow(log2,2)*logz-144*z3*pow(log2,2)*logz+144*z4*pow(log2,2)*logz+648*z*(1+z)*log1mz*logz+216*_L0*z*(1+z)*log1mz*logz+96*_nl*z*(1+z)*log1mz*logz-576*z2*(1+z)*log1mz*logz-720*z3*(1+z)*log1mz*logz+648*_L0*z3*(1+z)*log1mz*logz+96*_nl*z3*(1+z)*log1mz*logz-72*z*(1+z)*(1+11*z2)*log1mz2*logz-216*logz2-288*z*logz2-18*_L0*z*logz2-24*_nl*z*logz2+108*z2*logz2-18*_L0*z2*logz2-24*_nl*z2*logz2+288*z3*logz2-126*_L0*z3*logz2-72*_nl*z3*logz2-36*z4*logz2-126*_L0*z4*logz2-72*_nl*z4*logz2+72*(-1+z)*z*(1+z*(5+2*z))*log1mz*logz2-36*z*logz3-12*z2*logz3+36*z3*logz3+12*z4*logz3-576*(-1+z)*z*log1mz*log(z/4.)*log1pz+432*(-1+z)*z*logz*log1pz+576*(-1+z)*z2*logz*log1pz-144*z*(1+z)*(2+z*(-3+2*z))*logz2*log1pz-24*(-1+z)*(Pi2*(z+z3)+12*z*(4+z2)*pow(log2,2)+6*(3+z*(5+z+z2*(-1+log(4)))-2*z*log(2-2*z))*log(2*z))*log1pz+432*(-1+z)*log1pz2+828*(-1+z)*z*log1pz2-36*(-1+z)*z2*(29+14*z)*log1pz2+576*z4*log2*log1pz2-144*z*(-13+15*z)*logz*log1pz2-288*(-1+z)*z*(log(4-4*z)+logz)*log1pz2-288*z3*log(4*z)*log1pz2-96*(-4+z)*(-1+z)*z*(4+z)*log1pz3+288*(-1+z)*z*(log(32)+(-1+z2)*log1mz-(-1+z2)*log1pz)*polylog(2,(1-z)/2.)+792*z*(1+z)*polylog(2,1-z)-216*_L0*z*(1+z)*polylog(2,1-z)-576*z3*(1+z)*polylog(2,1-z)+216*_L0*z3*(1+z)*polylog(2,1-z)+72*z*(1+z)*log(pow(-1+z,8))*polylog(2,1-z)-144*z*(1+z)*logz*polylog(2,1-z)+288*z*(1+z)*log1pz*polylog(2,1-z)+144*z2*(1+z)*((4-8*z)*log1mz+z*logz+2*(-2+z)*log1pz)*polylog(2,1-z)+144*(-1+z2)*(3+z*(2+z*(-1+log(4))-log(4))-2*(-1+z)*z*log1mz)*polylog(2,(-1+z)/(2.*z))+432*(-1+z)*polylog(2,-z)+1440*(-1+z)*z*polylog(2,-z)-72*(-1+z)*z2*(-1+15*z)*polylog(2,-z)-288*z*((-1+z)*log1mz+(1+z3)*logz+(-5+6*z+z2)*log1pz)*polylog(2,-z)-432*(-1+z)*z*polylog(2,1/(1+z))-288*(-1+z)*z2*polylog(2,1/(1+z))+576*(-1+z)*z3*polylog(2,1/(1+z))+864*(-1+z)*z*log1mz*polylog(2,1/(1+z))+288*z*(-4+5*z+z3)*log1pz*polylog(2,1/(1+z))+4*B0*M_PI*z*(-36*B1overB0*M_PI*(-1+z)*(-1+z2)-2*(1-z)*(18*_L0*(1+z)*(3+z)-2*_nl*(1+z)*(-17+11*z)+3*(9+Pi2-15*z+3*(-8+Pi2)*z2))+36*(3-3*(-1+z)*z+4*_L0*pow(-1+z,2)*(1+z))*log1mz-216*pow(-1+z,2)*(1+z)*log1mz2-108*z3*(log1mz-logz)-54*logz-36*_L0*logz-12*_nl*logz-18*z*logz-36*_L0*z*logz-12*_nl*z*logz+144*z2*logz-108*_L0*z2*logz-12*_nl*z2*logz-108*_L0*z3*logz-12*_nl*z3*logz+36*(1+z)*log1mz*logz+180*z2*(1+z)*log1mz*logz-27*logz2+9*z*logz2-9*z2*logz2+27*z3*logz2-72*(-1+z)*(-2+z2)*logz*log1pz-108*(-1+z)*log1pz2-72*(1+z)*polylog(2,1-z)+72*z2*(1+z)*polylog(2,1-z)-72*(-1+z)*(-2+z2)*polylog(2,-z)-216*(-1+z)*polylog(2,1/(1+z)))-72*(-1+z2)*(-6+z*(-8+5*z)+10*(-1+z)*z*log1pz)*polylog(2,z/(1+z))+144*(-1+z2)*(3+z*(2+z*(-1+log(4))-log(4))-2*(-1+z)*z*log1pz)*polylog(2,-1+2/(1+z))-288*pow(-1+z,2)*z*(1+z)*polylog(3,(1-z)/2.)+144*(-1+z)*z*(1+z)*(4+7*z)*polylog(3,1-z)+288*pow(-1+z,2)*z*(1+z)*polylog(3,(-1+z)/(2.*z))-144*z2*(-1+z2)*polylog(3,(-1+z)/z)-144*z*(1+z)*(7+z*(-16+7*z))*polylog(3,-z)+72*(-1+z)*z*(1+z)*(5+7*z)*polylog(3,z)+288*z*(-3+z+z2-3*z3)*polylog(3,1/(1+z))-1008*pow(-1+z,2)*z*(1+z)*polylog(3,z/(1+z))+288*pow(-1+z,2)*z*(1+z)*polylog(3,(2*z)/(1+z))-288*(-1+z)*z*(4+z2)*polylog(3,(1+z)/2.)+18*z*(29+z*(59+z*(-69+61*z)))*zeta3)/(36.*z*(-1+z2));

    // std::cout << std:: endl << "test: " << JLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JLL2 << std::endl;
    // std::cout << std:: endl << "test: " << JLL3 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL0 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL2 << std::endl;

    return {JLL1, JLL2, JLL3, JNLL0, JNLL1, JNLL2};
  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::RecBarS(double const& z) const
  {
    double JLL1(0.0), JLL2(0.0), JLL3(0.0), JNLL0(0.0), JNLL1(0.0), JNLL2(0.0);

    const double log1mz = log(1-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = pow(log1mz,3);

    double B0(0.0), B02(0.0), B1overB0(0.0);
    if( ! _aQED.IsFixed() )
      {
        B0 = _b0;
        B1overB0 = _b1/_b0;
        B02 = _b0*_b0;
      }

    if (_orderR[0])
      JLL1 = -2+2./(1-z);
    // -----------------------------------------------------------------------
    if (_orderR[1])
      JLL2 = -2+6./(1-z)-8*log1mz+(8*log1mz)/(1-z);
    // -----------------------------------------------------------------------
    if (_orderR[2])
      JLL3 = 4.5+4*Pi2+(13.5-4*Pi2)/(1-z)-12*log1mz+(36*log1mz)/(1-z)-24*log1mz2+(24*log1mz2)/(1-z);
    // -----------------------------------------------------------------------
    if (_orderR[3])
      JNLL0 = 2-2*_L0+(-2+2*_L0)/(1-z)+4*log1mz-(4*log1mz)/(1-z);
    // -----------------------------------------------------------------------
    if (_orderR[4])
      JNLL1 = (9-12*_nl-18*z+32*_nl*z-36*B0*M_PI*z+36*B1overB0*M_PI*z-12*Pi2*z+18*_L0*(-2-z+2*B0*M_PI*z)-18*(-2+(-5+4*_L0+4*B0*M_PI)*z)*log1mz+108*z*log1mz2)/(9.*(-1+z));
    // -----------------------------------------------------------------------
    if (_orderR[5])
      JNLL2 = (-9*_L0*(36+16*B02*Pi2*z-(9+8*Pi2)*z-16*B0*M_PI*(2+z))+2*_nl*(44+(22+8*Pi2)*z-8*B0*M_PI*(-3+8*z))+2*(16*_nl*(-3+13*z)+36*_L0*(-6+(-3+8*B0*M_PI)*z)+3*(72+48*B02*Pi2*z+(-21+96*B1overB0*M_PI-40*Pi2)*z-24*B0*M_PI*(3+4*z)))*log1mz-216*(-2+(-3+2*_L0+4*B0*M_PI)*z)*log1mz2+576*z*log1mz3+6*(-15-8*Pi2-12*z-10*Pi2*z+24*B02*Pi2*z+24*B1overB0*M_PI*(2+z)+2*B0*M_PI*(-15+(21-12*B1overB0*M_PI+8*Pi2)*z)+120*z*zeta3))/(18.*(-1+z));

    return {JLL1, JLL2, JLL3, JNLL0, JNLL1, JNLL2};
  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::RecHatS(double const& z) const
  {
    double JLL1(0.0), JLL2(0.0), JLL3(0.0), JNLL0(0.0), JNLL1(0.0), JNLL2(0.0);

    const double z2 = z*z;
    const double z3 = pow(z,3);
    const double z4 = pow(z,4);

    const double log1mz = log(1-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = pow(log1mz,3);

    const double logz = log(z);
    const double logz2 = logz*logz;
    const double logz3 = pow(logz,3);

    const double log1pz = log(1+z);
    const double log1pz2 = log1pz*log1pz;

    const double log2 = log(2.);

    double B0(0.0), B02(0.0), B1overB0(0.0);
    if( ! _aQED.IsFixed() )
      {
        B0 = _b0;
        B1overB0 = _b1/_b0;
        B02 = _b0*_b0;
      }

    if (_orderR[0])
      JLL1 = 1-z;
    // -----------------------------------------------------------------------
    if (_orderR[1])
      JLL2 = -3-z-(2*_nl*(-1+z)*(4+z*(7+4*z)))/(3.*z)-4*(-1+z)*log1mz+((1-4*_nl+(3+4*_nl)*z2)*logz)/(-1+z);
    // -----------------------------------------------------------------------
    if (_orderR[2])
      JLL3 = -(-64*pow(_nl,2)-513*z-552*_nl*z+16*pow(_nl,2)*z-36*Pi2*z+96*_nl*Pi2*z+378*z2+1104*_nl*z2+96*pow(_nl,2)*z2+144*Pi2*z2+135*z3-552*_nl*z3+16*pow(_nl,2)*z3-108*Pi2*z3-96*_nl*Pi2*z3-64*pow(_nl,2)*z4+432*pow(-1+z,2)*z*log1mz2-54*z*logz-96*pow(_nl,2)*z*logz-432*z2*logz+288*_nl*z2*logz-162*z3*logz+96*_nl*z3*logz+96*pow(_nl,2)*z3*logz-384*_nl*z4*logz+18*z*logz2-144*_nl*z*logz2+126*z3*logz2+144*_nl*z3*logz2+24*log1mz*((-1+z)*(9*z*(3+z)+4*_nl*(-4-3*z+3*z2+4*z3))-18*(z+z3)*logz)+72*(3+8*_nl)*z*(-1+z2)*polylog(2,z))/(36.*(-1+z)*z);
    // -----------------------------------------------------------------------
    if (_orderR[3])
      JNLL0 = -1+_L0*(1-z)+z-2*(1-z)*log1mz;
    // -----------------------------------------------------------------------
    if (_orderR[4])
      JNLL1 = -(-72*_nl+48*_L0*_nl+9*z-54*_L0*z-176*_nl*z+36*_L0*_nl*z+36*B0*M_PI*z-36*B1overB0*M_PI*z-36*B0*_L0*M_PI*z+42*Pi2*z+9*z2-18*_L0*z2+368*_nl*z2-84*_L0*_nl*z2-36*B0*M_PI*z2+36*B1overB0*M_PI*z2+36*B0*_L0*M_PI*z2-42*Pi2*z2-9*z3+54*_L0*z3+104*_nl*z3-84*_L0*_nl*z3-36*B0*M_PI*z3+36*B1overB0*M_PI*z3+36*B0*_L0*M_PI*z3-18*Pi2*z3-9*z4+18*_L0*z4-296*_nl*z4+36*_L0*_nl*z4+36*B0*M_PI*z4-36*B1overB0*M_PI*z4-36*B0*_L0*M_PI*z4+18*Pi2*z4+72*_nl*pow(z,5)+48*_L0*_nl*pow(z,5)+18*z*log1mz+72*_L0*z*log1mz+72*B0*M_PI*z*log1mz+54*z2*log1mz-72*_L0*z2*log1mz-72*B0*M_PI*z2*log1mz-18*z3*log1mz-72*_L0*z3*log1mz-72*B0*M_PI*z3*log1mz-54*z4*log1mz+72*_L0*z4*log1mz+72*B0*M_PI*z4*log1mz-108*z*log1mz2+108*z2*log1mz2+108*z3*log1mz2-108*z4*log1mz2-96*_nl*logz-18*_L0*z*logz-192*_nl*z*logz+72*_L0*_nl*z*logz-18*_L0*z2*logz+120*_nl*z2*logz+72*_L0*_nl*z2*logz+18*z3*logz-54*_L0*z3*logz+264*_nl*z3*logz-72*_L0*_nl*z3*logz+18*z4*logz-54*_L0*z4*logz-48*_nl*z4*logz-72*_L0*_nl*z4*logz-96*_nl*pow(z,5)*logz+72*z3*log1mz*logz+72*z4*log1mz*logz+9*z*logz2-108*_nl*z*logz2-27*z2*logz2-108*_nl*z2*logz2+27*z3*logz2+108*_nl*z3*logz2-9*z4*logz2+108*_nl*z4*logz2+144*z*logz*log1pz-144*z2*logz*log1pz-72*z3*logz*log1pz+72*z4*logz*log1pz-108*z*log1pz2+108*z2*log1pz2+36*(-1+z)*z*pow(1+z,2)*polylog(2,1-z)+72*z*(2-2*z-z2+z3)*polylog(2,-z)-216*z*polylog(2,1/(1+z))+216*z2*polylog(2,1/(1+z)))/(18.*z*(-1+z2));
    // -----------------------------------------------------------------------
    if (_orderR[5])
      JNLL2 = (-1152*_nl+32*pow(_nl,2)+192*_L0*pow(_nl,2)+1152*B1overB0*_nl*M_PI+96*_nl*Pi2-1350*z+1539*_L0*z-3660*_nl*z+1656*_L0*_nl*z+720*pow(_nl,2)*z+144*_L0*pow(_nl,2)*z-1296*B1overB0*M_PI*z+864*B1overB0*_nl*M_PI*z+72*Pi2*z+216*_L0*Pi2*z+120*_nl*Pi2*z+486*z2+405*_L0*z2+5004*_nl*z2-1656*_L0*_nl*z2-176*pow(_nl,2)*z2-336*_L0*pow(_nl,2)*z2-432*B1overB0*M_PI*z2-2016*B1overB0*_nl*M_PI*z2+108*Pi2*z2-216*_L0*Pi2*z2-216*_nl*Pi2*z2+1350*z3-1539*_L0*z3+4092*_nl*z3-1656*_L0*_nl*z3-1328*pow(_nl,2)*z3-336*_L0*pow(_nl,2)*z3+1296*B1overB0*M_PI*z3-2016*B1overB0*_nl*M_PI*z3-504*Pi2*z3-216*_L0*Pi2*z3-216*_nl*Pi2*z3-486*z4-405*_L0*z4-3852*_nl*z4+1656*_L0*_nl*z4+144*pow(_nl,2)*z4+144*_L0*pow(_nl,2)*z4+432*B1overB0*M_PI*z4+864*B1overB0*_nl*M_PI*z4+324*Pi2*z4+216*_L0*Pi2*z4+120*_nl*Pi2*z4-432*_nl*pow(z,5)+608*pow(_nl,2)*pow(z,5)+192*_L0*pow(_nl,2)*pow(z,5)+1152*B1overB0*_nl*M_PI*pow(z,5)+96*_nl*Pi2*pow(z,5)+288*Pi2*z*(-1+2*z+z3)*log2-216*(-1+z)*(1+z)*(3+z*(7+9*z))*pow(log2,2)-144*z*(-11+z*(16+z+4*z2))*pow(log2,3)-432*B02*Pi2*pow(-1+z,2)*z*(1+z)*(-1+_L0-2*log1mz)+1728*B1overB0*M_PI*pow(-1+z,2)*z*(1+z)*log1mz+6*(16*pow(_nl,2)*(-3+z)*z*(1+z)*(-1+3*z)+4*_nl*(8+12*_L0*(-4+z*(-3+7*z*(1+z)))+z*(181+z*(-237+z*(-117+229*z))))-3*(-1+z)*z*(36*_L0*(1+z)*(3+z)+Pi2*(-56-4*z+44*z2)+3*(-37+(-36+z)*z-40*pow(log2,2))))*log1mz+1512*(-1+z)*z*(1+z)*log1mz2+1296*_L0*(-1+z)*z*(1+z)*log1mz2+1080*(-1+z)*z2*(1+z)*log1mz2-1296*_L0*(-1+z)*z2*(1+z)*log1mz2+144*_nl*pow(-1+z,2)*(1+z)*(4+7*z+4*z2)*log1mz2+1728*pow(-1+z,2)*z*(1+z)*log1mz3+96*_nl*(4*_nl-9*_L0*z4+4*(-4-3*_L0+_nl)*pow(z,5))*(log1mz-logz)-432*B1overB0*M_PI*z*(1+z)*(1-4*_nl+(3+4*_nl)*z2)*logz-6*(3*_L0*z*(16*pow(_nl,2)*(-1+z)*pow(1+z,2)+16*_nl*z*(3+4*z)-9*(1+z)*(1+z*(8+3*z)))-216*log2+2*z*(8*pow(_nl,2)*(1+z)*(-7+z*(-13+16*z))+3*Pi2*(-13+z*(9+(5-9*z)*z))+2*_nl*(238+18*Pi2*(-1+z)*(1+z)*(2+z)+z*(-83+z*(-152+105*z)))+144*z*log2+9*(1+z)*(3-z2+8*(-1+z)*log2+4*pow(-1+z,2)*pow(log2,2))))*logz+216*z*log1mz*logz-576*pow(_nl,2)*(-1+z)*z*pow(1+z,2)*log1mz*logz-216*z2*(7+10*z)*log1mz*logz+648*_L0*z*(1+z)*(1+3*z2)*log1mz*logz+144*_nl*(1+z)*(8+z*(15+8*z3-9*z*(2+z)+12*_L0*(-1+z2)))*log1mz*logz-216*(1+z)*(z-4*_nl*z+(11+4*_nl)*z3)*log1mz2*logz+648*logz2+432*z*logz2-504*_nl*z*(1+z)*logz2-1440*_nl*z4*(1+z)*logz2+288*pow(_nl,2)*(-1+z)*z*pow(1+z,2)*logz2+72*_nl*z2*(1+z)*(15+8*z)*logz2-54*_L0*z*(1+z)*(1-8*_nl+(7+8*_nl)*z2)*logz2+108*z2*(-13+z*(-8+7*z))*logz2-216*(-1+z)*z*(-5-z*(5+2*z)+2*_nl*(-1+z+2*z2))*log1mz*logz2+36*z*(-1+z2)*(-1+5*z+12*_nl*(2+z))*logz3+72*z*(Pi2*(-1+z)*(1+z2)+(-z3+2*(-1+z)*(4+z2)*log2)*log(64))*log1pz+1728*(-1+z)*z*log1mz*log(z/4.)*log1pz-1296*(-1+z)*z*logz*log1pz-1728*(-1+z)*z2*logz*log1pz+432*z*(1+z)*(2+z*(-3+2*z))*logz2*log1pz+432*(-3+2*z*(-1+z*(2+z))+2*(-1+z)*z3*log2-2*(-1+z)*z*log(2-2*z))*log(2*z)*log1pz-1296*(-1+z)*log1pz2-2484*(-1+z)*z*log1pz2+108*(-1+z)*z2*(29+14*z)*log1pz2-1728*z4*log2*log1pz2+432*z*(-13+15*z)*logz*log1pz2+864*z3*log(4*z)*log1pz2-432*z4*logz*(log((1-z)/2.)+log1pz)+288*(-1+z)*z*log1pz2*(3*(log(4-4*z)+logz)+(-16+z2)*log1pz)+864*(-1+z)*z*(-((-1+z2)*log1mz)+z2*log1pz-log(32*(1+z)))*polylog(2,(1-z)/2.)+648*z*(1+z)*polylog(2,1-z)+648*_L0*z*(1+z)*(-1+z2)*polylog(2,1-z)-576*pow(_nl,2)*z*(1+z)*(-1+z2)*polylog(2,1-z)+288*_nl*(-1+z)*(1+z)*(-2+6*_L0*z*(1+z)+(-3+z)*z*(3+2*z))*polylog(2,1-z)-432*(-1+z)*z*(1+z)*(4*(2+_nl+z+_nl*z)*log1mz+(-1+2*_nl*(-2+z)-z)*logz+2*(-1+z)*log1pz)*polylog(2,1-z)-432*(-1+z2)*(3+z*(2+z*(-1+log(4))-log(4))-2*(-1+z)*z*log1mz)*polylog(2,(-1+z)/(2.*z))-1296*(-1+z)*polylog(2,-z)-4320*(-1+z)*z*polylog(2,-z)+216*(-1+z)*z2*(-1+15*z)*polylog(2,-z)+864*z*((-1+z)*log1mz+(1+z3)*logz+(-5+6*z+z2)*log1pz)*polylog(2,-z)+1296*(-1+z)*z*polylog(2,1/(1+z))+864*(-1+z)*z2*polylog(2,1/(1+z))-1728*(-1+z)*z3*polylog(2,1/(1+z))-2592*(-1+z)*z*log1mz*polylog(2,1/(1+z))-864*z*(-4+5*z+z3)*log1pz*polylog(2,1/(1+z))+12*B0*M_PI*(-224*_nl+96*_L0*_nl+90*z-108*_L0*z-212*_nl*z+72*_L0*_nl*z-36*B1overB0*M_PI*z+54*Pi2*z-36*_L0*z2+556*_nl*z2-168*_L0*_nl*z2+36*B1overB0*M_PI*z2-54*Pi2*z2-90*z3+108*_L0*z3+292*_nl*z3-168*_L0*_nl*z3+36*B1overB0*M_PI*z3-30*Pi2*z3+36*_L0*z4-332*_nl*z4+72*_L0*_nl*z4-36*B1overB0*M_PI*z4+30*Pi2*z4-80*_nl*pow(z,5)+96*_L0*_nl*pow(z,5)+36*z*(-1+z2)*(4*_L0*(-1+z)-3*(1+z))*log1mz+216*(-1+z)*z*(1+z)*log1mz2-216*(-1+z)*z2*(1+z)*log1mz2-6*(1+z)*(3*z*(-1-2*z*(1+z)+_L0*(2+6*z2))+4*_nl*(8+z*(9+4*(-3+z)*z*(1+z)+6*_L0*(-1+z2))))*logz+36*z*(1+z)*(1+5*z2)*log1mz*logz-9*pow(-1+z,3)*z*logz2+180*_nl*z*(1+z)*(-1+z2)*logz2+72*(-1+z)*z*(-2+z2)*logz*log1pz+108*(-1+z)*z*log1pz2+72*z*(1+z)*(-1+z2)*polylog(2,1-z)+72*(-1+z)*z*(-2+z2)*polylog(2,-z)+216*(-1+z)*z*polylog(2,1/(1+z)))+216*(-1+z2)*(-6+z*(-8+5*z)+10*(-1+z)*z*log1pz)*polylog(2,z/(1+z))+432*(-1+z)*(1+z)*(-3+z*(-2+z+log(4)-z*log(4))+2*(-1+z)*z*log1pz)*polylog(2,-1+2/(1+z))+864*pow(-1+z,2)*z*(1+z)*polylog(3,(1-z)/2.)+432*z*(-1+z2)*(8+3*z+4*_nl*(1+z))*polylog(3,1-z)-864*pow(-1+z,2)*z*(1+z)*polylog(3,(-1+z)/(2.*z))-432*z*(6*_nl+z)*(-1+z2)*polylog(3,(-1+z)/z)+432*z*(1+z)*(7+z*(-16+7*z))*polylog(3,-z)+216*z*(9+20*_nl+3*z+8*_nl*z)*(-1+z2)*polylog(3,z)+864*z*(1+z)*(3+z*(-4+3*z))*polylog(3,1/(1+z))+3024*pow(-1+z,2)*z*(1+z)*polylog(3,z/(1+z))-864*pow(-1+z,2)*z*(1+z)*polylog(3,(2*z)/(1+z))+864*(-1+z)*z*(4+z2)*polylog(3,(1+z)/2.)-54*z*(-107+16*_nl*(-1+z)*(1+z)*(5+2*z)+z*(99+z*(67+21*z)))*zeta3)/(108.*z*(-1+z2));

    return {JLL1, JLL2, JLL3, JNLL0, JNLL1, JNLL2};
  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::RecBarG(double const& z) const
  {
    double JLL1(0.0), JLL2(0.0), JLL3(0.0), JNLL0(0.0), JNLL1(0.0), JNLL2(0.0);

    const double log1mz = log(1-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = pow(log1mz,3);

    double B0(0.0), B02(0.0), B1overB0(0.0);
    if( ! _aQED.IsFixed() )
      {
        B0 = _b0;
        B1overB0 = _b1/_b0;
        B02 = _b0*_b0;
      }

    if (_orderR[0])
      JLL1 = 1.0;
    // -----------------------------------------------------------------------
    if (_orderR[1])
      JLL2 = 1.5-(2*_nl)/3.+2*log1mz;
    // -----------------------------------------------------------------------
    if (_orderR[2])
      JLL3 = 2.25-_nl+(4*pow(_nl,2))/9.-(2*Pi2)/3.+(6-(4*_nl)/3.)*log1mz+4*log1mz2;
    // -----------------------------------------------------------------------
    if (_orderR[3])
      JNLL0 = -1+_L0;
    // -----------------------------------------------------------------------
    if (_orderR[4])
      JNLL1 = -4+(3*_L0)/2.-(2*(13+3*_L0)*_nl)/9.-2*B1overB0*M_PI-2*B0*(-1+_L0)*M_PI+(-7+2*_L0-(4*_nl)/3.)*log1mz-3*log1mz2;
    // -----------------------------------------------------------------------
    if (_orderR[5])
      JNLL2 = -5.625-(23*_nl)/6.+(52*pow(_nl,2))/27.+4*B0*M_PI-6*B1overB0*M_PI+(40*B0*_nl*M_PI)/9.+(8*B1overB0*_nl*M_PI)/3.+(11*Pi2)/6.-4*B02*Pi2+4*B0*B1overB0*Pi2+(2*_nl*Pi2)/9.+_L0*(2.25+(4*pow(_nl,2))/9.-6*B0*M_PI-(2*Pi2)/3.+4*B02*Pi2+_nl*(-1+(8*B0*M_PI)/3.))+(-18.5+(8*pow(_nl,2))/9.+18*B0*M_PI-8*B1overB0*M_PI+2*Pi2+_L0*(6-(4*_nl)/3.-8*B0*M_PI)+(4*_nl*(-5+2*B0*M_PI))/3.)*log1mz+(-18.5+4*_L0-(2*_nl)/3.+10*B0*M_PI)*log1mz2-6*log1mz3-6*zeta3;

    // std::cout << std:: endl << "test: " << JLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JLL2 << std::endl;
    // std::cout << std:: endl << "test: " << JLL3 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL0 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL2 << std::endl;

    return {JLL1, JLL2, JLL3, JNLL0, JNLL1, JNLL2};
  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::RecHatG(double const& z) const
  {
    double JLL1(0.0), JLL2(0.0), JLL3(0.0), JNLL0(0.0), JNLL1(0.0), JNLL2(0.0);

    const double z2 = z*z;
    const double z3 = pow(z,3);
    const double z4 = pow(z,4);

    const double log1mz = log(1-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = pow(log1mz,3);

    const double logz = log(z);
    const double logz2 = logz*logz;
    const double logz3 = pow(logz,3);

    const double log1pz = log(1+z);
    const double log1pz2 = log1pz*log1pz;
    const double log1pz3 = pow(log1pz,3);

    const double log2 = log(2.);

    double B0(0.0), B02(0.0), B1overB0(0.0);
    if( ! _aQED.IsFixed() )
      {
        B0 = _b0;
        B1overB0 = _b1/_b0;
        B02 = _b0*_b0;
      }

    if (_orderR[0])
      JLL1 = -3+2/z+z;
    // -----------------------------------------------------------------------
    if (_orderR[1])
      JLL2 = -((-1+z)*(4*_nl*(-2+z)+3*z)-12*(2-3*z+z2)*log1mz+6*(-2+z)*z*logz)/(6.*z);
    // -----------------------------------------------------------------------
    if (_orderR[2])
      JLL3 = (-99+1068*_nl-48*pow(_nl,2)+72*Pi2-(496*_nl)/z+(32*pow(_nl,2))/z+99*z-636*_nl*z+16*pow(_nl,2)*z-24*Pi2*z+64*_nl*z2+(144*(2-3*z+z2)*log1mz2)/z-180*logz-48*_nl*logz-(192*_nl*logz)/z+72*z*logz+240*_nl*z*logz-36*logz2+144*_nl*logz2+18*z*logz2-72*_nl*z*logz2-(24*log1mz*((-1+z)*(2*_nl*(-2+z)+3*z)+6*(2-2*z+z2)*logz))/z-(288*polylog(2,z))/z)/36.;
    // -----------------------------------------------------------------------
    if (_orderR[3])
      JNLL0 = ((-1+_L0)*(2-3*z+z2)-2*(2-2*z+z2)*logz)/z;
    // -----------------------------------------------------------------------
    if (_orderR[4])
      JNLL1 = (45+9*_L0+108*_nl+36*_L0*_nl+108*B1overB0*M_PI-(56*_nl)/z-(24*_L0*_nl)/z-(72*B1overB0*M_PI)/z-45*z-9*_L0*z-52*_nl*z-12*_L0*_nl*z-36*B1overB0*M_PI*z-(54*(2-3*z+z2)*log1mz2)/z+36*_L0*logz-48*_nl*logz+(48*_nl*logz)/z+45*z*logz-18*_L0*z*logz+24*_nl*z*logz-18*logz2+9*z*logz2+(6*log1mz*((-1+z)*(12+8*_nl+6*_L0*(-2+z)-9*z-4*_nl*z)+6*(-2+z)*z*logz))/z-(36*B0*M_PI*((-1+_L0)*(2-3*z+z2)-2*(2-2*z+z2)*logz))/z+36*(-2+z)*polylog(2,1-z))/18.;
    // -----------------------------------------------------------------------
    if (_orderR[5])
      JNLL2 = (-2808+2176*_nl+2976*_L0*_nl-448*pow(_nl,2)-192*_L0*pow(_nl,2)-768*B0*_nl*M_PI-1152*B1overB0*_nl*M_PI-1152*B0*_L0*_nl*M_PI+144*Pi2+1728*B02*Pi2-1728*B0*B1overB0*Pi2+288*_L0*Pi2-1728*B02*_L0*Pi2-96*_nl*Pi2+3645*z+594*_L0*z+12420*_nl*z-6408*_L0*_nl*z+864*pow(_nl,2)*z+288*_L0*pow(_nl,2)*z+648*B0*M_PI*z+432*B1overB0*M_PI*z+432*B0*_L0*M_PI*z+1728*B0*_nl*M_PI*z+1728*B1overB0*_nl*M_PI*z+1728*B0*_L0*_nl*M_PI*z+324*Pi2*z-2592*B02*Pi2*z+2592*B0*B1overB0*Pi2*z-432*_L0*Pi2*z+2592*B02*_L0*Pi2*z+144*_nl*Pi2*z+1971*z2-594*_L0*z2-16324*_nl*z2+840*_L0*_nl*z2+32*pow(_nl,2)*z2+96*_L0*pow(_nl,2)*z2-648*B0*M_PI*z2-432*B1overB0*M_PI*z2-432*B0*_L0*M_PI*z2-192*B0*_nl*M_PI*z2+576*B1overB0*_nl*M_PI*z2+576*B0*_L0*_nl*M_PI*z2-180*Pi2*z2-864*B02*Pi2*z2+864*B0*B1overB0*Pi2*z2-144*_L0*Pi2*z2+864*B02*_L0*Pi2*z2+48*_nl*Pi2*z2-3645*z3-594*_L0*z3-12868*_nl*z3+6024*_L0*_nl*z3-864*pow(_nl,2)*z3-288*_L0*pow(_nl,2)*z3-648*B0*M_PI*z3-432*B1overB0*M_PI*z3-432*B0*_L0*M_PI*z3-1728*B0*_nl*M_PI*z3-1728*B1overB0*_nl*M_PI*z3-1728*B0*_L0*_nl*M_PI*z3-468*Pi2*z3+2592*B02*Pi2*z3-2592*B0*B1overB0*Pi2*z3+432*_L0*Pi2*z3-2592*B02*_L0*Pi2*z3-144*_nl*Pi2*z3+837*z4+594*_L0*z4+14148*_nl*z4-3816*_L0*_nl*z4+416*pow(_nl,2)*z4+96*_L0*pow(_nl,2)*z4+648*B0*M_PI*z4+432*B1overB0*M_PI*z4+432*B0*_L0*M_PI*z4+960*B0*_nl*M_PI*z4+576*B1overB0*_nl*M_PI*z4+576*B0*_L0*_nl*M_PI*z4+180*Pi2*z4-864*B02*Pi2*z4+864*B0*B1overB0*Pi2*z4-144*_L0*Pi2*z4+864*B02*_L0*Pi2*z4+48*_nl*Pi2*z4+448*_nl*pow(z,5)+384*_L0*_nl*pow(z,5)-864*Pi2*z*log2+864*Pi2*z3*log2+2016*z*pow(log2,3)-2016*z3*pow(log2,3)-384*pow(_nl,2)*log(2-2*z)+432*_L0*z2*log(2-2*z)+192*pow(_nl,2)*z2*log(2-2*z)-432*_L0*z4*log(2-2*z)+192*pow(_nl,2)*z4*log(2-2*z)+1080*log1mz+2304*_nl*log1mz+576*_L0*_nl*log1mz-864*B0*M_PI*log1mz+3456*B1overB0*M_PI*log1mz+3456*B0*_L0*M_PI*log1mz-1152*B0*_nl*M_PI*log1mz-576*Pi2*log1mz-3564*z*log1mz-432*_L0*z*log1mz-3744*_nl*z*log1mz-864*_L0*_nl*z*log1mz+576*pow(_nl,2)*z*log1mz+1296*B0*M_PI*z*log1mz-5184*B1overB0*M_PI*z*log1mz-5184*B0*_L0*M_PI*z*log1mz+1728*B0*_nl*M_PI*z*log1mz+1008*Pi2*z*log1mz+1404*z2*log1mz-864*_nl*z2*log1mz-288*_L0*_nl*z2*log1mz+432*B0*M_PI*z2*log1mz-1728*B1overB0*M_PI*z2*log1mz-1728*B0*_L0*M_PI*z2*log1mz+576*B0*_nl*M_PI*z2*log1mz+144*Pi2*z2*log1mz+3564*z3*log1mz+432*_L0*z3*log1mz+3744*_nl*z3*log1mz+864*_L0*_nl*z3*log1mz-576*pow(_nl,2)*z3*log1mz-1296*B0*M_PI*z3*log1mz+5184*B1overB0*M_PI*z3*log1mz+5184*B0*_L0*M_PI*z3*log1mz-1728*B0*_nl*M_PI*z3*log1mz-1008*Pi2*z3*log1mz-2484*z4*log1mz-1440*_nl*z4*log1mz-288*_L0*_nl*z4*log1mz+432*B0*M_PI*z4*log1mz-1728*B1overB0*M_PI*z4*log1mz-1728*B0*_L0*M_PI*z4*log1mz+576*B0*_nl*M_PI*z4*log1mz+432*Pi2*z4*log1mz+1728*log1mz2-1728*_L0*log1mz2+288*_nl*log1mz2-4320*B0*M_PI*log1mz2-2700*z*log1mz2+2592*_L0*z*log1mz2-432*_nl*z*log1mz2+6480*B0*M_PI*z*log1mz2-756*z2*log1mz2+864*_L0*z2*log1mz2-144*_nl*z2*log1mz2+2160*B0*M_PI*z2*log1mz2+2700*z3*log1mz2-2592*_L0*z3*log1mz2+432*_nl*z3*log1mz2-6480*B0*M_PI*z3*log1mz2-972*z4*log1mz2+864*_L0*z4*log1mz2-144*_nl*z4*log1mz2+2160*B0*M_PI*z4*log1mz2+1728*log1mz3-3024*z*log1mz3-432*z2*log1mz3+3024*z3*log1mz3-1296*z4*log1mz3-2880*_nl*logz+1152*_L0*_nl*logz+2304*B0*_nl*M_PI*logz+3456*B02*Pi2*logz-2916*z*logz+1080*_L0*z*logz+18288*_nl*z*logz+288*_L0*_nl*z*logz-384*pow(_nl,2)*z*logz-864*B0*M_PI*z*logz+1728*B1overB0*M_PI*z*logz+1728*B0*_L0*M_PI*z*logz-2304*B0*_nl*M_PI*z*logz+288*Pi2*z*logz-3456*B02*Pi2*z*logz-702*z2*logz+3960*_nl*z2*logz-2592*_L0*_nl*z2*logz+648*B0*M_PI*z2*logz-864*B1overB0*M_PI*z2*logz-864*B0*_L0*M_PI*z2*logz-1152*B0*_nl*M_PI*z2*logz-864*Pi2*z2*logz-1728*B02*Pi2*z2*logz+2916*z3*logz-1080*_L0*z3*logz-17520*_nl*z3*logz-288*_L0*_nl*z3*logz+384*pow(_nl,2)*z3*logz+864*B0*M_PI*z3*logz-1728*B1overB0*M_PI*z3*logz-1728*B0*_L0*M_PI*z3*logz+2304*B0*_nl*M_PI*z3*logz-864*Pi2*z3*logz+3456*B02*Pi2*z3*logz+702*z4*logz-1080*_nl*z4*logz+1440*_L0*_nl*z4*logz-648*B0*M_PI*z4*logz+864*B1overB0*M_PI*z4*logz+864*B0*_L0*M_PI*z4*logz-1152*B0*_nl*M_PI*z4*logz+288*Pi2*z4*logz-1728*B02*Pi2*z4*logz-768*_nl*pow(z,5)*logz-1296*z*log1mz*logz-1728*_L0*z*log1mz*logz-3456*B0*M_PI*z*log1mz*logz-648*z2*log1mz*logz+864*_L0*z2*log1mz*logz+1728*B0*M_PI*z2*log1mz*logz+1296*z3*log1mz*logz+1728*_L0*z3*log1mz*logz+3456*B0*M_PI*z3*log1mz*logz+648*z4*log1mz*logz-864*_L0*z4*log1mz*logz-1728*B0*M_PI*z4*log1mz*logz+3024*z*log1mz2*logz-1512*z2*log1mz2*logz-3024*z3*log1mz2*logz+1512*z4*log1mz2*logz-1152*_nl*logz2-432*z*logz2+216*_L0*z*logz2+1728*_nl*z*logz2-864*_L0*_nl*z*logz2-432*B0*M_PI*z*logz2-432*z2*logz2-108*_L0*z2*logz2+2016*_nl*z2*logz2+432*_L0*_nl*z2*logz2+216*B0*M_PI*z2*logz2+432*z3*logz2-216*_L0*z3*logz2-1728*_nl*z3*logz2+864*_L0*_nl*z3*logz2+432*B0*M_PI*z3*logz2+432*z4*logz2+108*_L0*z4*logz2-864*_nl*z4*logz2-432*_L0*_nl*z4*logz2-216*B0*M_PI*z4*logz2-2592*z*log1mz*logz2+3024*z2*log1mz*logz2+4320*z3*log1mz*logz2-1296*z4*log1mz*logz2-144*z*logz3+1152*_nl*z*logz3-936*z2*logz3-576*_nl*z2*logz3-432*z3*logz3-1152*_nl*z3*logz3+360*z4*logz3+576*_nl*z4*logz3+384*pow(_nl,2)*log(2*z)-432*_L0*z2*log(2*z)-192*pow(_nl,2)*z2*log(2*z)+432*_L0*z4*log(2*z)-192*pow(_nl,2)*z4*log(2*z)+576*Pi2*log1pz+2016*Pi2*z*log1pz-2016*Pi2*z2*log1pz-576*Pi2*z3*log1pz+2592*logz*log1pz+2592*z*logz*log1pz-864*z2*logz*log1pz-3456*z3*logz*log1pz-864*z4*logz*log1pz+864*logz2*log1pz+864*z*logz2*log1pz-432*z2*logz2*log1pz-864*z3*logz2*log1pz-432*z4*logz2*log1pz+1296*log1pz2+3456*z*log1pz2-648*z2*log1pz2-3456*z3*log1pz2-648*z4*log1pz2+864*logz*log1pz2+10368*z*logz*log1pz2-9504*z2*logz*log1pz2-1728*z3*logz*log1pz2-1440*log1pz3-4608*z*log1pz3+3888*z2*log1pz3+288*z3*log1pz3+1872*z4*log1pz3-144*(-1+z2)*(12-12*_L0+4*_nl-15*z-4*_nl*z-24*B0*M_PI*z+6*z2+2*_nl*z2+12*B0*M_PI*z2+(12+18*z-9*z2)*log1mz+3*(-6+z)*z*logz)*polylog(2,1-z)+432*(-1+z)*(-12-30*z-25*z2-3*z3-8*pow(1+z,2)*logz+2*z*(-12+z+3*z2)*log1pz)*polylog(2,-z)+1728*z*polylog(2,1/(1+z))-864*z3*polylog(2,1/(1+z))-864*z4*polylog(2,1/(1+z))-10368*z*log1pz*polylog(2,1/(1+z))+7776*z2*log1pz*polylog(2,1/(1+z))+1728*z3*log1pz*polylog(2,1/(1+z))+864*z4*log1pz*polylog(2,1/(1+z))+2592*polylog(2,z/(1+z))+5184*z*polylog(2,z/(1+z))-1296*z2*polylog(2,z/(1+z))-6048*z3*polylog(2,z/(1+z))-432*z4*polylog(2,z/(1+z))+1728*logz*polylog(2,z/(1+z))-1728*z2*logz*polylog(2,z/(1+z))+3456*z*log1pz*polylog(2,z/(1+z))-3456*z2*log1pz*polylog(2,z/(1+z))-3456*z3*log1pz*polylog(2,z/(1+z))+3456*z4*log1pz*polylog(2,z/(1+z))-4320*z*polylog(3,1-z)+6480*z2*polylog(3,1-z)+7776*z3*polylog(3,1-z)-3024*z4*polylog(3,1-z)+5184*z2*polylog(3,(-1+z)/z)+3456*z3*polylog(3,(-1+z)/z)-1728*z4*polylog(3,(-1+z)/z)+1728*polylog(3,z)-4320*z*polylog(3,z)+5616*z2*polylog(3,z)+7776*z3*polylog(3,z)-3888*z4*polylog(3,z)-432*polylog(3,z2)-432*z*polylog(3,z2)+216*z2*polylog(3,z2)+432*z3*polylog(3,z2)+216*z4*polylog(3,z2)+3456*z*polylog(3,z/(1+z))-1728*z2*polylog(3,z/(1+z))-3456*z3*polylog(3,z/(1+z))+1728*z4*polylog(3,z/(1+z))-3456*zeta3+2592*z*zeta3-864*z2*zeta3-6048*z3*zeta3+864*z4*zeta3)/(216.*z*(-1+z2));

    // std::cout << std:: endl << "test: " << JLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JLL2 << std::endl;
    // std::cout << std:: endl << "test: " << JLL3 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL0 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL2 << std::endl;

    return {JLL1, JLL2, JLL3, JNLL0, JNLL1, JNLL2};
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::AsyEleAF(double const& z) const
  {
    const double lambda0 = 3./4.;

    if ( _orderA == 0 )
      {
        return exp( (lambda0 - emc) * _eta0 ) / tgamma( 1 + _eta0 )
               * _eta0 * pow( 1 - z, -1 + _eta0 );
      }

    const double lambda1 = 3./8. - Pi2/2. + 6.*zeta3 - _nl/18. * (3. + 4.*Pi2);
    const double eta1 = _eta0 * ( 1 - 10./9. * _atwopi * _nl);
    const double etahat1 = _eta0 * ( lambda0 + _atwopi/2. * lambda1 );
    const double factor = exp( etahat1 - emc * eta1 ) / tgamma( 1 + eta1 )
                          * eta1 * pow( 1 - z, -1 + eta1 );
    const double psieta1 = psi(eta1);
    const double A = - emc - psieta1;
    const double B = 0.5 * emc*emc + Pi2/12. + emc*psieta1 + 0.5 * psieta1*psieta1 - 0.5 * dpsi(eta1,1);
    const double curl = 1. + _a0twopi*2. * ( (_L0 - 1)*( A + 3./4. ) - 2*B + 7./4.
                                             + (_L0 - 1 - 2*A)*log(1-z) - log(1-z)*log(1-z) );
    return factor * curl;
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::AsyEleAR(double const& z) const
  {
    if ( _orderA == 0 )
      {
        const double csihat0 = 3./2. * _t;
        const double csi0 = 2.*_t;

        return exp( csihat0 - emc * csi0 ) / tgamma( 1 + csi0 )
               * csi0 * pow( 1 - z, -1 + csi0 );
      }

    const double lambda1 = 3./8. - Pi2/2. + 6.*zeta3 - _nl/18. * (3. + 4.*Pi2);
    const double csi1 = 2. * _t - _atwopi / (2.*M_PI*_b0) * ( 1 - exp(-2*M_PI*_b0*_t) )
                        * (20./9. * _nl + 4.*M_PI*_b1/_b0);
    const double csihat1 = 3./2. * _t + _atwopi / (2.*M_PI*_b0) * ( 1 - exp(-2*M_PI*_b0*_t) )
                           * (lambda1 - 3.*M_PI*_b1/_b0);
    const double factor = exp( csihat1 - emc * csi1 ) / tgamma( 1 + csi1 )
                          * csi1 * pow( 1 - z, -1 + csi1 );
    const double psicsi1 = psi(csi1);
    const double A = - emc - psicsi1;
    const double B = 0.5 * emc*emc + Pi2/12. + emc*psicsi1 + 0.5 * psicsi1*psicsi1 - 0.5 * dpsi(csi1,1);
    const double curl = 1. + _a0twopi*2 * ( (_L0 - 1)*( A + 3./4. ) - 2*B + 7./4.
                                            + (_L0 - 1 - 2*A)*log(1-z) - log(1-z)*log(1-z) );

    return factor*curl;

  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::AsyEleBarAR(double const& z) const
  {
    double JLL1(0.0), JLL2(0.0), JLL3(0.0), JNLL0(0.0), JNLL1(0.0), JNLL2(0.0);

    const double log1mz = log(1-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = pow(log1mz,3);
    const double b02 = _b0*_b0;

    if (_orderR[0])
      JLL1 = 2/(1 - z);
    // -----------------------------------------------------------------------
    if (_orderR[1])
      JLL2 = (2*(-3 - 4*log1mz))/(-1 + z);
    // -----------------------------------------------------------------------
    if (_orderR[2])
      JLL3 = (-27 + 8*Pi2 - 72*log1mz - 48*log1mz2)/(2.*(-1 + z));
    // -----------------------------------------------------------------------
    if (_orderR[3])
      JNLL0 = (2*(-1+_L0))/(1-z)-(4*log1mz)/(1-z);
    // -----------------------------------------------------------------------
    if (_orderR[4])
      JNLL1= (1-(20*_nl)/9.+4*_b0*M_PI-(4*_b1*M_PI)/_b0+(4*Pi2)/3.+_L0*(6-4*_b0*M_PI))/(1-z)-(2*(7-4*_L0-4*_b0*M_PI)*log1mz)/(1-z)-(12*log1mz2)/(1-z);
    // -----------------------------------------------------------------------
    if (_orderR[5])
      JNLL2 = 2*(-(((8.5+(80*_nl)/9.-28*_b0*M_PI+(16*_b1*M_PI)/_b0-(20*Pi2)/3.+8*b02*Pi2-2*_L0*(9-8*_b0*M_PI))*log1mz)/(1-z))-(6*(5-2*_L0-4*_b0*M_PI)*log1mz2)/(1-z)-(16*log1mz3)/(1-z)+(4.5-(12*_b1*M_PI)/_b0-4*b02*Pi2+(3+4*_b1)*Pi2-_nl*(3.6666666666666665-(20*_b0*M_PI)/9.+(4*Pi2)/9.)-(2*_b0*M_PI*(3+4*Pi2))/3.+_L0*(6.75-12*_b0*M_PI-2*Pi2+4*b02*Pi2)-20*zeta3)/(1-z));

    // std::cout << std:: endl << "test: " << JLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JLL2 << std::endl;
    // std::cout << std:: endl << "test: " << JLL3 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL0 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL1 << std::endl;
    // std::cout << std:: endl << "test: " << JNLL2 << std::endl;

    return {JLL1, JLL2, JLL3, JNLL0, JNLL1, JNLL2};
  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::AsyEleBarAF(double const& z) const
  {
    double JLL1(0.0), JLL2(0.0), JLL3(0.0), JNLL0(0.0), JNLL1(0.0), JNLL2(0.0);

    const double log1mz = log(1-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = pow(log1mz,3);

    if (_orderR[0])
      JLL1 = 2/(1 - z);
    // -----------------------------------------------------------------------
    if (_orderR[1])
      JLL2 = (2*(-3 - 4*log1mz))/(-1 + z);
    // -----------------------------------------------------------------------
    if (_orderR[2])
      JLL3 = (-27 + 8*Pi2 - 72*log1mz - 48*log1mz2)/(2.*(-1 + z));
    // -----------------------------------------------------------------------
    if (_orderR[3])
      JNLL0 = (2*(-1+_L0))/(1-z)-(4*log1mz)/(1-z);
    // -----------------------------------------------------------------------
    if (_orderR[4])
      JNLL1= (1+6*_L0-(20*_nl)/9.+(4*Pi2)/3.)/(1-z)-(2*(7-4*_L0)*log1mz)/(1-z)-(12*log1mz2)/(1-z);
    // -----------------------------------------------------------------------
    if (_orderR[5])
      JNLL2 = 2*(-(((8.5-18*_L0+(80*_nl)/9.-(20*Pi2)/3.)*log1mz)/(1-z))-(6*(5-2*_L0)*log1mz2)/(1-z)-(16*log1mz3)/(1-z)+(4.5+3*Pi2+_L0*(6.75-2*Pi2)-_nl*(3.6666666666666665+(4*Pi2)/9.)-20*zeta3)/(1-z));

    return {JLL1, JLL2, JLL3, JNLL0, JNLL1, JNLL2};
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::AsyPhoton(double const& z) const
  {
    double b1overb0(0.0), b0(0.0), t(0.0), chi10(0.0);

    if( ! _aQED.IsFixed() )
      {
        b0 = _b0;
        b1overb0 = _b1/_b0;
        t = _t;
        chi10 = 0.0;
      }
    else
      {
        b0 = 0.0;
        b1overb0 = 0.0;
        t = _eta0/2.0;
        chi10 = _a0twopi * _nl;
      }

    const double lambda1 = 3./8. - Pi2/2. + 6.*zeta3 - _nl/18. * (3. + 4.*Pi2);
    const double csi10 = 2.0 * ( 1.0 - _a0twopi*2.0 * (5.0/9.0 * _nl + M_PI*b1overb0) );
    const double csihat10 = 3.0/2.0 * ( 1.0 + _a0twopi*2.0 * (lambda1/3.0 - M_PI*b1overb0) );
    const double D11 = csi10;
    const double D21 = - (2.0/3.0 *_nl + 2.0*M_PI*b0 + csihat10 + chi10);
    const double C11 = _a0twopi * exp(-D21*t);
    const double C21 = - _a0twopi * ( 5.0 + 4.0/3.0 * _nl ) * exp(-D21*t);
    const double C31 = _a0twopi * (6.0 + Pi2/6.0 + 32.0/9.0 * _nl + 2.0*M_PI*b1overb0) * exp(-D21*t);
    const double D12 = D11;
    const double D22 = D21;
    const double C12 = - _a0twopi;
    const double C22 = _a0twopi*(5.0 + 4.0/3.0*_nl);
    const double C32 = - _a0twopi*(6.0 + Pi2/6.0 + 32.0/9.0*_nl + 2.0*M_PI*b1overb0);
    const double D13 = D11;
    const double D23 = - (2.0/3.0*_nl + csihat10 + chi10);
    const double C13 = 0.0;
    const double C23 = 0.0;
    const double C33 = - exp(-D23*t);
    const double D14 = D13;
    const double D24 = D23;
    const double C14 = 0.0;
    const double C24 = 0.0;
    const double C34 = 1.0;

    const double k1 = csi10*t;
    const double k3 = csi10*t;
    const double k2 = 0.0;
    const double k4 = 0.0;

    const double gamma1 = sumRiMi(C11, C21, C31, D21/D11, D11, z, k1, D11, D21);
    const double gamma2 = sumRiMi(C12, C22, C32, D22/D12, D12, z, k2, D12, D22);
    const double gamma3 = sumRiMi(C13, C23, C33, D23/D13, D13, z, k3, D13, D23);
    const double gamma4 = sumRiMi(C14, C24, C34, D24/D14, D14, z, k4, D14, D24);
    const double gamma5 = _a0twopi * (1.0 + (1-z)*(1-z))/z * (_L0 - 2*log(z) - 1.0);
    const double sumgamma = gamma1 + gamma2 + gamma3 + gamma4 + gamma5;

    return exp( - ( 2.0/3.0 * _nl + chi10 ) * t ) * sumgamma;

  }

  //_________________________________________________________________________________
  double AnalyticSolutions::sumRiMi(double const& C1, double const& C2, double const& C3,
                                    double const& C4, double const& C5,
                                    double const& z, double const& k, double const& M1, double const& M2) const
  {
    const double F0 = 2.0 - Pi2/3.0 + 3.0/2.0 * _L0;
    const double F1 = 2*(1-_L0);
    const double F2 = -2.0;

    const double factor = exp(-emc*k)*pow(1-z,k)/tgamma(1+k);

    const double log1mz = log(1.0-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = log1mz*log1mz2;

    const double den = M2 - M1*log1mz;
    const double den2 = den*den;
    const double den3 = den*den2;

    const double k2 = k*k;

    const double MF1 = 1.0 / den - (Pi2*k - 6*zeta3*k2)*M1/6.0 / den2
                       - (30.0*Pi2 - 360*zeta3*k + Pi2*Pi2*k2)*M1*M1/180.0 / den3;
    const double MF2 = 1.0;
    const double MF3 = - log1mz + Pi2*k/6.0 - zeta3*k2;
    const double MF4 = log1mz2 - Pi2/6.0 + k*(-Pi2/3.0*log1mz + 2*zeta3)
                       + k2*(2*zeta3*log1mz - Pi2*Pi2/180.0);
    const double MF5 = - log1mz3 + Pi2/2.0*log1mz - 2*zeta3
                       + k*(Pi2/2.0*log1mz2 - 6*zeta3*log1mz - Pi2*Pi2/60.0)
                       + k2*(-3*zeta3*log1mz2 + Pi2*Pi2/60.0*log1mz + 3.0/2.0*Pi2*zeta3 - 12*zeta5);

    const double R1 = (C3 - C4*C2 + C4*C4*C1) * (1 + _a0twopi * (F0 - C4*F1 + C4*C4*F2));
    const double R2 = 1/C5 * (C2 - C4*C1)
                      + _a0twopi * (1/C5) * (C2*F0 + C3*F1 - C4*(C1*F0 + C2*F1 + C3*F2)
                                             + C4*C4*(C1*F1 + C2*F2) - C4*C4*C4*C1*F2);
    const double R3 = C1/C5 + _a0twopi * 1/C5 * (C1*F0 + C2*F1 + C3*F2
                                                 - C4*(C1*F1 + C2*F2) + C4*C4*C1*F2);
    const double R4 = _a0twopi * 1/C5 * (C1*F1 + C2*F2 - C4*C1*F2);
    const double R5 = _a0twopi * C1/C5 * F2;

    return factor * (R1*MF1 + R2*MF2 + R3*MF3 + R4*MF4 + R5*MF5);
  }

  //_________________________________________________________________________________
  double AnalyticSolutions::AsyPhotonSIMPLIFIED(double const& z) const
  {
    const double a0 = _a0twopi*2.0*M_PI;
    const double a = _atwopi*2.0*M_PI;
    const double csi10 = 2.0 * ( 1.0 - _a0twopi*2.0 * (5.0/9.0 * _nl + M_PI*_b1/_b0) );

    //    return a0*a0/a * 3.0/(2.0*M_PI*csi10) * log(1-z)
    //     - a0*a0*a0/a * 1.0/(2.0*Pi2*csi10) * pow(log(1-z),3);

    const double gamma5 = _a0twopi * (1.0 + (1-z)*(1-z))/z * (_L0 - 2*log(z) - 1.0);

    const double coeff1 = 3.0/csi10 * log(1-z)
                          + (-31.0 -12*_nl + csi10*(42.0-12.0*_L0+8*_nl))/(6*csi10*csi10);
    const double coeff2 = - 2.0/csi10 * pow(log(1-z),3)
                          + log(1-z)*log(1-z)* (13+csi10*(-36+6*_L0-8*_nl)+4*_nl)/(3.0*csi10*csi10);

    return (a0/a) * ( gamma5 + _a0twopi * coeff1 + _a0twopi*_a0twopi * coeff2 );
  }

  //_________________________________________________________________________________
  std::vector<double> AnalyticSolutions::AsyPhotonBar(double const& z) const
  {
    double JLL1(0.0), JLL2(0.0), JLL3(0.0), JNLL0(0.0), JNLL1(0.0), JNLL2(0.0);

    const double log1mz = log(1-z);
    const double log1mz2 = log1mz*log1mz;
    const double log1mz3 = pow(log1mz,3);

    double B0(0.0), B02(0.0), B1overB0(0.0);
    if( ! _aQED.IsFixed() )
      {
        B0 = _b0;
        B1overB0 = _b1/_b0;
        B02 = _b0*_b0;
      }

    if (_orderR[0])
      JLL1 = 1.0;
    // -----------------------------------------------------------------------
    if (_orderR[1])
      JLL2 = 1.5-(2*_nl)/3.+2*log1mz;
    // -----------------------------------------------------------------------
    if (_orderR[2])
      JLL3 = 2.25-_nl+(4*pow(_nl,2))/9.-(2*Pi2)/3.+(6-(4*_nl)/3.)*log1mz+4*log1mz2;
    // -----------------------------------------------------------------------
    if (_orderR[3])
      JNLL0 = -1+_L0;
    // -----------------------------------------------------------------------
    if (_orderR[4])
      JNLL1 = -4+(3*_L0)/2.-(2*(13+3*_L0)*_nl)/9.-2*B1overB0*M_PI-2*B0*(-1+_L0)*M_PI+(-7+2*_L0-(4*_nl)/3.)*log1mz-3*log1mz2;
    // -----------------------------------------------------------------------
    if (_orderR[5])
      JNLL2 = -5.625-(23*_nl)/6.+(52*pow(_nl,2))/27.+4*B0*M_PI-6*B1overB0*M_PI+(40*B0*_nl*M_PI)/9.+(8*B1overB0*_nl*M_PI)/3.+(11*Pi2)/6.-4*B02*Pi2+4*B0*B1overB0*Pi2+(2*_nl*Pi2)/9.+_L0*(2.25+(4*pow(_nl,2))/9.-6*B0*M_PI-(2*Pi2)/3.+4*B02*Pi2+_nl*(-1+(8*B0*M_PI)/3.))+(-18.5+(8*pow(_nl,2))/9.+18*B0*M_PI-8*B1overB0*M_PI+2*Pi2+_L0*(6-(4*_nl)/3.-8*B0*M_PI)+(4*_nl*(-5+2*B0*M_PI))/3.)*log1mz+(-18.5+4*_L0-(2*_nl)/3.+10*B0*M_PI)*log1mz2-6*log1mz3-6*zeta3;

    return {JLL1, JLL2, JLL3, JNLL0, JNLL1, JNLL2};
  }
}
