//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "ePDF/alphaem.h"
#include "ePDF/constants.h"

namespace ePDF
{
  //_________________________________________________________________________________
  AlphaQED::AlphaQED(YAML::Node const& config):
    AlphaQED(config["Alpha"]["evolution"].as<std::string>(),
             config["Perturbative order"].as<std::string>(),
             config["Alpha"]["aref"].as<double>(),
             config["Alpha"]["Qref"].as<double>(),
             config["NL"].as<int>(), 0)
  {
  }

  //_________________________________________________________________________________
  AlphaQED::AlphaQED(std::string const& method, std::string const& po,
                     double const& aref, double const& Qref,
                     int const& NL, int const& NF):
    m_method(method),
    m_po(po),
    m_aref(aref/4/M_PI),
    m_Qref(Qref),
    m_NL(NL),
    m_NF(NF)
  {
  }

  //_________________________________________________________________________
  double AlphaQED::Evolve(double const& Q) const
  {
    if (m_method == "analytic")
      {
        return EvolveAnalytic(Q);
      }
    else if (m_method == "numerical")
      {
        return EvolveNumerical(Q);
      }
    else if (m_method == "fixed")
      {
        return m_aref;
      }
    else
      {
        throw std::runtime_error("[AlphaQED::Evolve]: Undefined evolution method");
      }
  }

  //_________________________________________________________________________
  double AlphaQED::EvolveAnalytic(double const& Q) const
  {
    double B0 = Beta0QED();
    double B1 = Beta1QED();

    double Q2 = Q * Q;
    double Q02 = m_Qref * m_Qref;
    double LR = log(Q2/Q02);

    if (m_po == "LL")
      {
        return m_aref / ( 1.0 + m_aref * B0 * LR );
      }
    else if (m_po == "NLL")
      {
        double arglog = (1 + m_aref * (B1/B0) + m_aref * B0 * LR)/
                        (1 + m_aref * (B1/B0));
        double den = 1.0 + m_aref * B0 * LR + m_aref * (B1/B0) * log(arglog);
        return m_aref/den;
      }
    else
      {
        throw std::runtime_error("[AlphaQED::EvolveAnalytic]: Undefined logarithmic order");
      }
  }


  //_________________________________________________________________________
  double AlphaQED::EvolveNumerical(double const& Q) const
  {
    double Q2  = Q*Q;
    double Q02 = m_Qref * m_Qref;
    double LR = log(Q2/Q02);

    double dLR = LR / m_nstep;

    // Fourth-order Runge-Kutta
    double alpha = m_aref;
    double xk0, xk1, xk2, xk3;
    for (int i = 0; i < m_nstep; i++)
      {
        xk0 = dLR * FBetaQED(alpha          );
        xk1 = dLR * FBetaQED(alpha + 0.5 * xk0);
        xk2 = dLR * FBetaQED(alpha + 0.5 * xk1);
        xk3 = dLR * FBetaQED(alpha +       xk2);
        alpha = alpha + m_sxth * ( xk0 + 2.0 * xk1 + 2.0 * xk2 + xk3 );
      }

    return alpha;
  }

  //_________________________________________________________________________
  double AlphaQED::Beta0QED() const
  {
    return - 4.0/3.0 * (NC * m_sumch2[m_NF] + m_NL);
  }

  //_________________________________________________________________________
  double AlphaQED::Beta1QED() const
  {
    return - 16.0/4.0 * (NC * m_sumch4[m_NF] + m_NL);
  }

  //_________________________________________________________________________
  double AlphaQED::FBetaQED(double const& alpha) const
  {
    if (m_po == "LL")
      {
        return - alpha*alpha * Beta0QED();
      }
    else if (m_po == "NLL")
      {
        return - alpha*alpha * ( Beta0QED() + alpha * Beta1QED() );
      }
    else
      {
        throw std::runtime_error("[AlphaQED::FBetaQED]: Undefined perturbative order");
      }
  }

  //_________________________________________________________________________
  bool AlphaQED::IsFixed() const
  {
    if (m_method == "fixed")
      return true;
    else
      return false;
  }
}
