//
// ePDF
//

#pragma once

#include "ePDF/ndistributions.h"

#include <yaml-cpp/yaml.h>
#include <complex>
#include <vector>
#include <memory>
#include <gsl/gsl_integration.h>
#include <functional>

namespace ePDF
{
  struct evol_params
  {
    double                          x;
    double                          Q;
    int                             id;
    std::shared_ptr<NDistributions> Ndist;
    //_________________________________________________________________________________
    evol_params operator = (evol_params const& p)
    {
      x     = p.x;
      Q     = p.Q;
      id    = p.id;
      Ndist = p.Ndist;
      return *this;
    }
  };

  /**
   * @brief The "xDistribution" class
   */
  class xDistributions
  {
  public:
    /**
     * @brief The "xDistribution" constructor.
     * @param config: the YAML:Node with the parameters
     */
    xDistributions(YAML::Node const& config);

    /**
     * @brief The "xDistribution" destructor.
     */
    ~xDistributions();

    /**
     * @brief Function that returns the PDFs in x space.
     * @param x: Bjorken x
     * @param Q: the final scale
     */
    std::vector<double> Evolve(double const& x, double const& Q);

    /**
     * @brief Interfaces to the integrand functions according to the
     * path chosen to perform the inverse Mellin transform.
     */
    std::vector<double> integrand(double const& y) const;
    std::vector<double> talbot(double const& y) const;
    std::vector<double> straight(double const& y) const;

    /**
     * @brief Integrators functions
     */
    std::vector<double> trapezoid(double const& a, double const& b) const;
    std::vector<double> gauss(double const& a, double const& b) const;
    std::vector<double> gaussGSL(double const& a, double const& b);

    /**
     * @brief Function that returns the N-th (complex) moment of PDFs.
     * @param N: moment
     * @param Q: the final scale
     */
    std::vector<std::complex<double>> Moments(std::complex<double> const& N, double const& Q) const;

    /**
     * @brief Function that sets the parameter structure externally
     * @param p: the input structure
     */
    void SetParameters(evol_params const& p) { _p = p; };

  private:
    std::string               const _contour;
    std::string               const _integrator;
    double                    const _eps;
    evol_params                     _p;
    gsl_integration_workspace      *_w;
  };
}
