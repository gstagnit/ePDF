//
// ePDF
//

#pragma once

#include "ePDF/alphaem.h"
#include "ePDF/constants.h"

namespace ePDF
{
  /**
   * @brief Utility definitions
   */
  const std::vector<bool> orderRLL{true, true, true, false, false, false};
  const std::vector<bool> orderRNLL{true, true, true, true, true, true};

  /**
   * @brief The "AnalyticSolutions" class that return the analytic
   * solutions.
   */
  class AnalyticSolutions
  {
  public:

    /**
     * @brief The "AnalyticSolutions" constructor.
     * @param config: the YAML:Node with the parameters
     * @note It reads the perturbative order from config file
     */
    AnalyticSolutions(YAML::Node const& config);

    /**
     * @brief The "AnalyticSolutions" constructor.
     * @param config: the YAML:Node with the parameters
     * @param orderA: 0 == LL, 1 = NLL
     * @param orderR: turn on or off contributions aL,(aL)^2,(aL)^3, a, a^2L, a^3L^2
     * @note It overwrites the perturbative order found in config file
     */
    AnalyticSolutions(YAML::Node const& config, int const& orderA, std::vector<bool> const orderR);

    /**
     * @brief Set the parameters for photon matching
     */
    void SetPhotonMatching(std::vector<double> const& vec);

    /**
     * @brief Turn on or off the numeric integrals (which vanish in
     * the z->1 limit)
     */
    void SetNumInt(bool const& numint);

    /**
     * @brief Turn on or off the terms vanishing in the z->1 limit
     */
    void SetRecHat(bool const& rechatON);

    /**
     * @brief Return only part of the solution
     * @param x: Bjorken x
     * @param Q: the final scale
     * @param term: A == asy, R == rec, ABAR = expansion of asy, RBAR = terms not vanishing in the z->1 limit
     */
    double GetSolution(double const& x, double const& Q, int const& id, std::string const& term);

    /**
     * @brief Function that returns the PDFs in x space.
     * @param x: Bjorken x
     * @param Q: the final scale
     * @param id: PDF "flavour" (0: singlet, 1: photon, 2: nonsinglet)
     */
    double Evolve(double const& x, double const& Q, int const& id);

    /**
     * @brief Function that returns all PDFs in x space as a vector.
     * @param x: Bjorken x
     * @param Q: the final scale
     */
    std::vector<double> Evolve(double const& x, double const& Q);

    /**
     * @brief Function that returns the photon PDF in x space.
     * @param x: Bjorken x
     * @param Q: the final scale
     * @param test: solution type
     */
    double TestPhoton(double const& x, double const& Q, std::string test);

  private:
    AlphaQED _aQED;               //!< Coupling object
    double const _Qi;             //!< Initial scale
    double const _nl;             //!< Active number of

    /**
     * @name Beta-function coefficients
     * Convert Beta0 e Beta1 (relative to alpha/(4 pi) expansion) in
     * b0 e b1 (relative to alpha expansion)
     */
    ///@{
    double const _b0;
    double const _b1;
    ///@}

    /**
     * @nameSave internal values during calculations
     */
    ///@{
    double const _a0twopi;
    double const _L0;
    double _atwopi;
    double _t;
    double _eta0;
    ///@}

    /**
     * @nameselect the accuracy of the perturbative solution
     */
    ///@{
    std::vector<bool> _orderR;
    int               _orderA;
    ///@}

    /**
     * @nameParameters for photon matching
     */
    ///@{
    double _x0photon;
    double _x1photon;
    double _pphoton;
    ///@}

    /**
     * @brief If true, numerical contributions calculated
     */
    bool _numint;

    /**
     * @brief If true, hat terms retained
     */
    bool _rechatON;

    /**
     * @name Utility functions
     */
    ///@{
    /**
     * @brief Build the perturbative series
     */
    double RecSeries(std::vector<double> const& series) const;

    /**
     * @brief Set up the Q dependent terms
     */
    void Warmup(double const& Q);

    /**
     * @brief Returns the different solutions. N.B. only the
     * x-dependent part, should be called after Warmup
     */
    double AsySolution(double const& x, int const& id) const;
    double AsyBarSolution(double const& x, int const& id) const;
    double RecSolution(double const& x, int const& id) const;
    double RecBarSolution(double const& x, int const& id) const;

    /**
     * @brief Returns a vector with the \bar{J}_k(z) for k = LL1, LL2,
     * LL3, NLL0, NLL1, NLL2
     */
    std::vector<double> RecBarNS(double const& z) const;
    std::vector<double> RecBarS(double const& z) const;
    std::vector<double> RecBarG(double const& z) const;

    /**
     * @brief Returns a vector with the \hat{J}_k(z) for k = LL1, LL2,
     * LL3, NLL0, NLL1, NLL2
     */
    std::vector<double> RecHatNS(double const& z) const;
    std::vector<double> RecHatS(double const& z) const;
    std::vector<double> RecHatG(double const& z) const;

    double AsyEleAF(double const& z) const;
    double AsyEleAR(double const& z) const;
    std::vector<double> AsyEleBarAF(double const& z) const;
    std::vector<double> AsyEleBarAR(double const& z) const;

    /**
     * @brief This contains AR and AF
     */
    double AsyPhoton(double const& z) const;

    /**
     * @brief Auxiliary function used in AsyPhoton
     */
    double sumRiMi(double const& C1, double const& C2, double const& C3,
                   double const& C4, double const& C5,
                   double const& z, double const& k, double const& M1, double const& M2) const;

    /**
     * @brief This contains AR and AF
     */
    std::vector<double> AsyPhotonBar(double const& z) const;

    /**
     * @note TO BE TESTED
     */
    double AsyPhotonSIMPLIFIED(double const& z) const;
    ///@}
  };
}
