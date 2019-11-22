//
// ePDF
//

#pragma once

#include "ePDF/matrix.h"

#include <complex>
#include <array>
#include <utility>

namespace ePDF
{
  /**
   * @brief function for the computation of the leading order splitting functions
   * @param N: complex argument
   * @param nf: the number of active flavours
   * @return The leading order splitting functions function in N
   */
  std::pair<std::complex<double>, Matrix<std::complex<double>>> andim_lo(std::complex<double> const& N, int const& nf);

  /**
   * @brief function for the computation of the next-to-leading order
   * splitting functions
   * @param N: complex argument
   * @param nf: the number of active flavours
   * @return The next-to-leading order splitting functions function in N
   */
  std::pair<std::complex<double>, Matrix<std::complex<double>>> andim_nlo(std::complex<double> const& N, int const& nf);
}
