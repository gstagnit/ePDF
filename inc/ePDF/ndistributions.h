//
// ePDF
//

#pragma once

#include "ePDF/alphaem.h"
#include "ePDF/matrix.h"

#include <yaml-cpp/yaml.h>
#include <complex>
#include <vector>

namespace ePDF
{
  /**
   * @brief The "NDistribution" class
   */
  class NDistributions
  {
  public:
    /**
     * @brief The "NDistribution" constructor.
     * @param config: the YAML:Node with the parameters
     */
    NDistributions(YAML::Node const& config);

    /**
     * @brief The "NDistribution" destructor.
     */
    virtual ~NDistributions() {};

    /**
     * @brief function that returns the N-th (complex) moment of PDFs.
     * @param N: moment
     * @param Q: the final scale
     */
    virtual std::vector<std::complex<double>> Evolve(std::complex<double> const& N, double const& Q) const;

    /**
     * @brief function that returns the evolution operators in N space
     * using the iterative U matrix method.
     * @param N: moment
     * @param Q: final scale coupling
     */
    std::pair<std::complex<double>, Matrix<std::complex<double>>> UMatrix(std::complex<double> const& N, double const& Q) const;

    /**
     * @brief function that returns the evolution operators in N space
     * using the path-ordering method.
     * @param N: moment
     * @param Q: final scale coupling
     */
    std::pair<std::complex<double>, Matrix<std::complex<double>>> PathOrdering(std::complex<double> const& N, double const& Q) const;

    /**
     * @brief function that returns the evolution operators in N space
     * using for the alpha-fixed case.
     * @param N: moment
     * @param Q: final scale coupling
     */
    std::pair<std::complex<double>, Matrix<std::complex<double>>> AlphaFixed(std::complex<double> const& N, double const& Q) const;

    /**
     * @brief function that returns the initial distributions in N space
     * @param N: moment
     */
    std::pair<std::complex<double>, Matrix<std::complex<double>>> InitialDistributions(std::complex<double> const& N) const;

    /**
     * @brief Scheme functions
     * @param N: moment
     */
    std::complex<double> fde(std::complex<double> const& N) const;
    std::complex<double> fdgm(std::complex<double> const& N) const;

    /**
     * @brief Helper functions
     */
    Matrix<std::complex<double>> Integrandsg(Matrix<std::complex<double>> const& G0,
                                             Matrix<std::complex<double>> const& G1,
                                             Matrix<std::complex<double>> const& J,
                                             double const& a) const;
    std::complex<double> Integrandns(std::complex<double> const& G0,
                                     std::complex<double> const& G1,
                                     std::complex<double> const& J,
                                     double const& a) const;
    Matrix<std::complex<double>> Minv(Matrix<std::complex<double>> const& M) const;

  private:
    std::map<std::string, int> const _ptmap = {{"LL", 0}, {"NLL", 1}};

    AlphaQED    const _aQED;
    int         const _ipt;
    int         const _nl;
    double      const _Qi;
    double      const _ai;
    double      const _bt0;
    double      const _bt1;
    std::string const _method;
    bool        const _expand;
    std::string const _scheme;
  };
}
