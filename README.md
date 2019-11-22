# ePDF: A code for the computation and the evolution of electron parton distribution in QED at the next-to-leading logarithmic accuracy

``ePDF`` is a library that implements the evolution in pure QED of the
unpolarised electron parton distribution functions (PDFs) up to
next-to-leading logarithmic (NLL) approximation. The initial
conditions, computed [here](https://arxiv.org/pdf/1909.03886.pdf), can
be evolved either *numerically*, by solving the DGLAP equation through
different numerical algorithms, or *analytically*. The analytical
solutions are obtained by means of an additive formula that matches a
large-*z* solution, that includes all orders in the QED coupling
constant, with a small- and intermediate-*z* solution that includes
terms up to third order in the coupling constant.

## Download

You can obtain ``ePDF`` directly from the github repository:

https://github.com/gstagnit/ePDF/releases

For the last development branch you can clone the master code:

```Shell
git clone https://github.com/gstagnit/ePDF.git
```

## Dependencies

In order to install the code, you need to have installed:

- [``yaml-cpp``](https://github.com/jbeder/yaml-cpp),
- [``GSL``](https://www.gnu.org/software/gsl/doc/html/).

## Installation 

The code can be compiled using the following procedure:

```Shell
cd ePDF
cmake -DCMAKE_INSTALL_PREFIX=/your/installation/path/ .
make && make install
```

By the default, if no prefix specification is given, the program will
be installed in the /usr/local folder. If you want (or need) to use a
different path, remember to export the ``ePDF`` /lib folder into the
LD_LIBRARY_PATH. More configuration options can be accessed through:

```Shell
ccmake .
```

## Documetation

Code documentation generated with Doxygen can be found here:
https://vbertone.github.io/ePDF/html/index.html.

## Relevant references

- S. Frixione, *Initial conditions for electron and photon structure and fragmentation functions*, [arXiv:1909.03886](https://arxiv.org/pdf/1909.03886.pdf).

## Contacts

- Valerio Bertone: valerio.bertone@cern.ch
- Matteo Cacciari: matteo.cacciari@cern.ch
- Stefano Frixione: stefano.frixione@cern.ch
- Giovanni Stagnitto: giovanni.stagnitto@gmail.com
