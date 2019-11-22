//
// Authors: Valerio Bertone
//          Giovanni Stagnitto
//

#include <ePDF/alphaem.h>
#include <ePDF/ndistributions.h>
#include <ePDF/xdistributions.h>
#include <ePDF/analyticsolutions.h>

#include <iostream>
#include <iomanip>
#include <string.h>

int main(int argc, char* argv[])
{
  // Check that the input is correct otherwise stop the code
  if (argc < 2 || strcmp(argv[1], "--help") == 0)
    {
      std::cout << "\nInvalid Parameters:" << std::endl;
      std::cout << "Syntax: ./Evolution <input card>" << std::endl;
      exit(-10);
    }

  // Final scale
  const double Q = 10;
  std::cout << std::scientific << std::setprecision(10) << "\nQ = " << Q << " GeV" << std::endl;

  // YAML config file
  const YAML::Node config = YAML::LoadFile(argv[1]);

  // Allocate AlphaQED object
  const ePDF::AlphaQED a{config};
  std::cout << std::setprecision(10) << "\nAlphaQED(Q) = " << a.Evolve(Q) << std::endl;

  // Allocate N-space PDFs
  ePDF::NDistributions npdfs{config};

  // First and second moment at the initial scale
  const std::vector<std::complex<double>> secmom0 = npdfs.Evolve(std::complex<double>(2.,0.), config["Initial scale"].as<double>());
  const std::vector<std::complex<double>> firmom0 = npdfs.Evolve(std::complex<double>(1.,0.), config["Initial scale"].as<double>());
  std::cout << "\nSum rules at the initial scale..." << std::endl;
  std::cout << "Momentum sum rule: " << secmom0[0] + secmom0[1] << std::endl;
  std::cout << "Valence sum rule: " << firmom0[2] << std::endl;

  // First and second moment at the final scale
  const std::vector<std::complex<double>> secmom = npdfs.Evolve(std::complex<double>(2.,0.), Q);
  const std::vector<std::complex<double>> firmom = npdfs.Evolve(std::complex<double>(1.,0.), Q);
  std::cout << "\nSum rules at the final scale..." << std::endl;
  std::cout << "Momentum sum rule: " << secmom[0] + secmom[1] << std::endl;
  std::cout << "Valence sum rule: " << firmom[2] << std::endl;

  // Allocate x-space PDFs
  ePDF::xDistributions xpdfs{config};

  // Tabulate PDFs
  const std::vector<double> xlha{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999};
  std::cout << "\nNumerical solution:" << std::endl;
  std::cout << "    x    "
            << "   e- + e+  "
            << "   photon   "
            << "   e- - e+  "
            << std::endl;
  for (auto const& x : xlha)
    {
      const std::vector<double> xf = xpdfs.Evolve(x, Q);
      std::cout << std::setprecision(2) << x << "  ";
      std::cout << std::setprecision(4) << xf[0] / x << "  "
                << std::setprecision(4) << xf[1] / x  << "  "
                << std::setprecision(4) << xf[2] / x  << "  "
                << std::endl;
    }
  std::cout << "\n";

  // Allocate AnalyticSolutions object
  std::cout << "Analytic solution:" << std::endl;
  ePDF::AnalyticSolutions pdfsan{config};
  std::cout << "    x    "
            << "   e- + e+  "
            << "   photon   "
            << "   e- - e+  "
            << std::endl;
  for (auto const& x : xlha)
    {
      const std::vector<double> xf = pdfsan.Evolve(x, Q);
      std::cout << std::setprecision(2) << x << "  ";
      std::cout << std::setprecision(4) << xf[0] << "  "
                << std::setprecision(4) << xf[1] << "  "
                << std::setprecision(4) << xf[2] << "  "
                << std::endl;
    }
  std::cout << "\n";

  return 0;
}
