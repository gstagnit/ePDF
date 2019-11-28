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

## Running the test code

If the installation was successful, the test code ``Evolution`` will be generated in the ``run/`` folder.
This code takes as an input a card containing the evolution settings. An example of such a card is in ``cards/default.yaml``.
To run the code type:

```Shell
cd run/
./Evolution ../cards/default.yaml
```

You should get the following output:

```Shell
Q = 1.0000000000e+01 GeV

AlphaQED(Q) = 5.8974514709e-04

Sum rules at the initial scale...
Momentum sum rule: (1.0000000000e+00,0.0000000000e+00)
Valence sum rule: (1.0000000000e+00,0.0000000000e+00)

Sum rules at the final scale...
Momentum sum rule: (9.9999999995e-01,0.0000000000e+00)
Valence sum rule: (1.0000000000e+00,0.0000000000e+00)

Numerical solution:
    x       e- + e+     photon      e- - e+  
1.00e-01  3.2572e-02  4.8721e-01  2.5949e-02  
2.00e-01  3.2542e-02  2.0718e-01  3.0245e-02  
3.00e-01  3.7671e-02  1.2061e-01  3.6540e-02  
4.00e-01  4.6451e-02  8.0161e-02  4.5805e-02  
5.00e-01  6.0240e-02  5.7539e-02  5.9841e-02  
6.00e-01  8.2542e-02  4.3574e-02  8.2286e-02  
7.00e-01  1.2188e-01  3.4411e-02  1.2171e-01  
8.00e-01  2.0418e-01  2.8134e-02  2.0408e-01  
9.00e-01  4.6076e-01  2.3595e-02  4.6071e-01  
9.50e-01  9.8654e-01  2.1589e-02  9.8652e-01  
9.90e-01  5.2837e+00  1.9248e-02  5.2837e+00  
9.99e-01  5.3974e+01  1.6660e-02  5.3974e+01  

Analytic solution:
    x       e- + e+     photon      e- - e+  
1.00e-01  3.2568e-02  4.8722e-01  2.5948e-02  
2.00e-01  3.2540e-02  2.0719e-01  3.0245e-02  
3.00e-01  3.7671e-02  1.2062e-01  3.6540e-02  
4.00e-01  4.6452e-02  8.0163e-02  4.5807e-02  
5.00e-01  6.0243e-02  5.7540e-02  5.9844e-02  
6.00e-01  8.2548e-02  4.3574e-02  8.2292e-02  
7.00e-01  1.2189e-01  3.4411e-02  1.2172e-01  
8.00e-01  2.0419e-01  2.8133e-02  2.0409e-01  
9.00e-01  4.6077e-01  2.3594e-02  4.6073e-01  
9.50e-01  9.8655e-01  2.1589e-02  9.8653e-01  
9.90e-01  5.2836e+00  1.9254e-02  5.2836e+00  
9.99e-01  5.3974e+01  1.6674e-02  5.3974e+01
```

Different options can be accessed by feeding the executable with the appropriate input card.

## Documetation

Code documentation generated with Doxygen can be found here:
https://vbertone.github.io/ePDF/html/index.html.

## Relevant references

- S. Frixione, *Initial conditions for electron and photon structure and fragmentation functions*, [arXiv:1909.03886](https://arxiv.org/pdf/1909.03886.pdf).
- V. Bertone, M. Cacciari, S. Frixione, G. Stagnitto, *The partonic structure of the electron at the next-to-leading logarithmic accuracy in QED*, [arXiv:1911.12040](https://arxiv.org/pdf/1911.12040.pdf)

## Contacts

- Valerio Bertone: valerio.bertone@cern.ch
- Matteo Cacciari: matteo.cacciari@cern.ch
- Stefano Frixione: stefano.frixione@cern.ch
- Giovanni Stagnitto: giovanni.stagnitto@gmail.com
