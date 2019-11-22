//
// ePDF
//

#pragma once

#include <yaml-cpp/yaml.h>
#include <vector>
#include <string>
#include <math.h>

namespace ePDF
{
  /**
   * @brief The "Alpha QED" class
   */
  class AlphaQED
  {
  public:
    /**
     * @brief The "Alpha QED" constructor.
     * @param config: the YAML:Node with the parameters
     */
    AlphaQED(YAML::Node const& config);

    /**
     * @brief The "Alpha QED" constructor.
     * @param method: the solution method
     * @param po: the perturbative order
     * @param aref: the reference value of alpha
     * @param Qref: the reference scale
     * @param NL: the number of active charged leptons
     * @param NF: the number of active quarks
     */
    AlphaQED(std::string const& method, std::string const& po,
             double const& aref, double const& Qref,
             int const& NL, int const& NF);

    /**
     * @brief Function that return the evolved alpha according to the
     * chose method.
     * @param Q: final scale
     * @note This functions returns alpha / (4pi)
     */
    double Evolve(double const& Q) const;

    /**
     * @brief Operator to return alpha at the scale Q.
     * @param Q: final scale
     */
    double operator()(double const& Q) const { return 4 * M_PI * Evolve(Q); };

    /**
     * @brief Function that return the alpha analytically
     * @param Q: final scale
     */
    double EvolveAnalytic(double const& Q) const;

    /**
     * @brief Function that return the alpha numerically
     * @param Q: final scale
     */
    double EvolveNumerical(double const& Q) const;

    /**
     * @brief The beta function coefficients
     */
    double Beta0QED() const;
    double Beta1QED() const;

    /**
     * @brief Whether alpha is fixed or runs
     */
    bool IsFixed() const;

    /**
     * @brief The beta function
     * @param alpha: the value of alpha to be used to compute the beta function
     */
    double FBetaQED(double const& alpha) const;

  private:
    std::string const m_method;
    std::string const m_po;
    double      const m_aref;
    double      const m_Qref;
    int         const m_NL;
    int         const m_NF;

    // Parameter of numerical integration
    int    const m_nstep = 10;
    double const m_sxth  = 0.166666666666666;

    // Sum of the electric charges squared
    std::vector<double> m_sumch2 = { 0.0, 1.0/9.0, 5.0/9.0, 2.0/3.0,
                                     10.0/9.0, 11.0/9.0, 5.0/3.0
                                   };

    // Sum of the electric charges to the fourth
    std::vector<double> m_sumch4 = { 0.0, 1.0/81.0, 17.0/81.0, 18.0/81.0,
                                     34.0/81.0, 35.0/81.0, 51.0/81.0
                                   };
  };
}
