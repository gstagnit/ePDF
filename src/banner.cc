//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "ePDF/banner.h"
#include "ePDF/version.h"

namespace ePDF
{
  //_________________________________________________________________________
  void Banner()
  {
    std::cout << "\033[1;33m\n";
    std::cout << "         ____    ____    _____ \n";
    std::cout << "   ___  |  _ \\  |  _ \\  |  ___|\n";
    std::cout << "  / _ \\ | |_) | | | | | | |_   \n";
    std::cout << " |  __/ |  __/  | |_| | |  _|  \n";
    std::cout << "  \\___| |_|     |____/  |_|    \n";
    std::cout << "\n";
    std::cout << "_____v" << VERSION << ": A code for the computation and the evolution of electron PDFs in QED at NLL\n";
    std::cout << "     webpage: https://github.com/gstagnit/ePDF\n";
    std::cout << "     Authors: V. Bertone, M. Cacciari, S. Frixione, G. Stagnitto\n";
    std::cout << "     If you use this code for a scientific publication, please cite:\n";
    std::cout << "     - S. Frixione, [arXiv:1909.03886]\n";
    std::cout << "     - V. Bertone, M. Cacciari, S. Frixione, G. Stagnitto, [arXiv:1911.12040]\n";
    std::cout << "\033[39m\n";
  }
}
