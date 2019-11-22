//
// ePDF
//

#pragma once

#include <gsl/gsl_integration.h>

namespace ePDF
{
  struct int_params
  {
    double z;
    int nl;
  };

  /**
   * @brief The "NumericIntegrals" class
   */
  class NumericIntegrals
  {
  public:
    /**
     * @brief The "NumericIntegrals" constructor.
     * @param flv: 0 = singlet, 1 = photon, 2 = nonsinglet
     * @param nl: number of leptons (needed for singlet)
     */
    NumericIntegrals(int const& iflv, int const& nl, double const& epsabs = 0.001, double const& epsrel = 0.001);

    /**
     * @brief The "NumericIntegrals" destructor.
     */
    ~NumericIntegrals();

    /**
     * @brief Perform the integration
     * @param z: momentum fraction
     */
    double integrate(double const& z) const;

  private:
    int const _iflv;
    int const _nl;
    double const _epsabs;
    double const _epsrel;
    gsl_integration_workspace *_w;
  };
}
