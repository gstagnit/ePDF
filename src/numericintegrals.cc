//
// Authors: Valerio Bertone
//          Giovanni Stagnitto
//

#include <ePDF/numericintegrals.h>
#include <ePDF/constants.h>
#include <ePDF/specialfunctions.h>

#include <math.h>
#include <iostream>
#include <vector>
#include <string>


namespace ePDF
{
  //_________________________________________________________________________________
  double lim1(double z)
  {
    double log1mz = log(1-z);
    return 2.*log1mz*log1mz + (2./3)*Pi2 * log1mz + (2./3)*Pi2 - 4.*zeta3;
  }

  //_________________________________________________________________________________
  double fun1(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      ((2 + (-1 + z)*y*(2 + (-1 + z)*y))*
       ((1 + (-1 + z)*y)*log(1 - z)*
        ((-1 + z)*log(1 - z) + 2*z*log(z)) -
        2*z*log(z/(1 + (-1 + z)*y))*
        log(1 - z/(1 + (-1 + z)*y)) +
        (-1 + z)*(-1 + y)*
        pow(log(1 - z/(1 + (-1 + z)*y)),2)))/
      ((-1 + z)*z*y*(1 + (-1 + z)*y)) +
      2*(1 + (-1 + z)*y)*
      (pow(log(1 - z),2)*log(z) -
       log(z/(1 + (-1 + z)*y))*
       pow(log(1 - z/(1 + (-1 + z)*y)),2));

    return result - lim1(z);
  }

  //_________________________________________________________________________________
  double lim2(double z)
  {
    double log1mz = log(1-z);
    return -2.*log1mz*log1mz + (2. - 2./3*Pi2) * log1mz + (Pi2/3) + 4.*zeta3;
  }

  //_________________________________________________________________________________
  double fun2(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = -(((2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*((-1 + z)*((1 + y*(-1 + z))*pow(log(1 - z),2) + (-1 + y)*pow(log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))),2)) - (z + y*(-1 + z)*z)*polylog(2,z) + z*polylog(2,z/(1 + y*(-1 + z)))))/(y*(1 + y*(-1 + z))*(-1 + z)*z)) + 2*(1 + y*(-1 + z))*(log(1 - z)*polylog(2,z) - log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z)))*polylog(2,z/(1 + y*(-1 + z))));

    return result - lim2(z);
  }

  //_________________________________________________________________________________
  double lim3(double z)
  {
    return -2.0;
  }

  //_________________________________________________________________________________
  double fun3(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = ((2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*((1 + y*(-1 + z))*polylog(2,1 - z) - polylog(2,((-1 + y)*(-1 + z))/(1 + y*(-1 + z)))))/(y*(1 + y*(-1 + z))*(-1 + z)) + 2*(1 + y*(-1 + z))*(polylog(3,1 - z) - polylog(3,((-1 + y)*(-1 + z))/(1 + y*(-1 + z))));

    return result - lim3(z);
  }

  //_________________________________________________________________________________
  double lim4(double z)
  {
    return 0.0;
  }

  //_________________________________________________________________________________
  double fun4(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      2*(1 + y*(-1 + z))*(pow(log(z),2)*log(1 + z) - pow(log(z/(1 + y*(-1 + z))),2)*log(1 + z/(1 + y*(-1 + z)))) - ((1 + pow(1 + y*(-1 + z),2))*(-(pow(log(z),2)/(1 + z)) - (2*log(z)*log(1 + z))/z - ((-1 + y)*log(z/(1 + y*(-1 + z)))*(z*log(z/(1 + y*(-1 + z))) + 2*(1 + y*(-1 + z) + z)*log(1 + z/(1 + y*(-1 + z)))))/((1 + y*(-1 + z))*z*(1 + y*(-1 + z) + z))))/y;

    return result - lim4(z);
  }

  //_________________________________________________________________________________
  double lim5(double z)
  {
    return 2*log(2)*log(2);
  }

  //_________________________________________________________________________________
  double fun5(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      -(((1 + pow(1 + y*(-1 + z),2))*((-2*log(z)*log(1 + z))/(1 + z) - pow(log(1 + z),2)/z - (2*(-1 + y)*log(z/(1 + y*(-1 + z)))*log(1 + z/(1 + y*(-1 + z))))/((1 + y*(-1 + z))*(1 + y*(-1 + z) + z)) - ((-1 + y)*pow(log(1 + z/(1 + y*(-1 + z))),2))/(z + y*(-1 + z)*z)))/y) + 2*(1 + y*(-1 + z))*(log(z)*pow(log(1 + z),2) - log(z/(1 + y*(-1 + z)))*pow(log(1 + z/(1 + y*(-1 + z))),2));

    return result - lim5(z);
  }

  //_________________________________________________________________________________
  double lim8(double z)
  {
    return - Pi2/6;
  }

  //_________________________________________________________________________________
  double fun8(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      -(((2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*(log(z)*log(1 + z) - y*log(z)*log(1 + z) + y*z*log(z)*log(1 + z) - log(z/(1 + y*(-1 + z)))*log(1 + z/(1 + y*(-1 + z))) + y*log(z/(1 + y*(-1 + z)))*log(1 + z/(1 + y*(-1 + z))) + (-1 + y - y*z)*polylog(2,-z) - (-1 + y)*polylog(2,-(z/(1 + y*(-1 + z))))))/(y*(1 + y*(-1 + z))*z)) + 2*(1 + y*(-1 + z))*(log(z)*polylog(2,-z) - log(z/(1 + y*(-1 + z)))*polylog(2,-(z/(1 + y*(-1 + z)))));

    return result - lim8(z);
  }

  //_________________________________________________________________________________
  double lim9(double z)
  {
    return -Pi2/12 - 2.*log(2.)*log(2.);
  }

  //_________________________________________________________________________________
  double fun9(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      -(((1 + pow(1 + y*(-1 + z),2))*(pow(log(1 + z),2)/z + ((-1 + y)*pow(log(1 + z/(1 + y*(-1 + z))),2))/(z + y*(-1 + z)*z) - polylog(2,-z)/(1 + z) - ((-1 + y)*polylog(2,-(z/(1 + y*(-1 + z)))))/((1 + y*(-1 + z))*(1 + y*(-1 + z) + z))))/y) + 2*(1 + y*(-1 + z))*(log(1 + z)*polylog(2,-z) - log(1 + z/(1 + y*(-1 + z)))*polylog(2,-(z/(1 + y*(-1 + z)))));

    return result - lim9(z);
  }

  //_________________________________________________________________________________
  double lim10(double z)
  {
    return Pi2/12 - 3./2*log(2.)*log(2.);
  }

  //_________________________________________________________________________________
  double fun10(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      -(((1 + pow(1 + y*(-1 + z),2))*(-((log(z/(1 + z))*log(1 + z))/(1 + z)) - ((-1 + y)*log(z/(1 + y*(-1 + z) + z))*log(1 + z/(1 + y*(-1 + z))))/((1 + y*(-1 + z))*(1 + y*(-1 + z) + z)) - polylog(2,1/(1 + z))/(1 + z) - ((-1 + y)*polylog(2,1/(1 + z/(1 + y*(-1 + z)))))/((1 + y*(-1 + z))*(1 + y*(-1 + z) + z))))/y) + 2*(1 + y*(-1 + z))*(log(1 + z)*polylog(2,1/(1 + z)) - log(1 + z/(1 + y*(-1 + z)))*polylog(2,1/(1 + z/(1 + y*(-1 + z)))));

    return result - lim10(z);
  }

  //_________________________________________________________________________________
  double lim11(double z)
  {
    return -Pi2/6;
  }

  //_________________________________________________________________________________
  double fun11(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      ((2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*((1 + y*(-1 + z))*polylog(2,-z) + (-1 + y)*polylog(2,-(z/(1 + y*(-1 + z))))))/(y*(1 + y*(-1 + z))*z) + 2*(1 + y*(-1 + z))*(polylog(3,-z) - polylog(3,-(z/(1 + y*(-1 + z)))));

    return result - lim11(z);
  }

  //_________________________________________________________________________________
  double lim12(double z)
  {
    return -Pi2/12 + 0.5*log(2.)*log(2.);
  }
  double fun12(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      -(((1 + pow(1 + y*(-1 + z),2))*(polylog(2,1/(1 + z))/(1 + z) + ((-1 + y)*polylog(2,1/(1 + z/(1 + y*(-1 + z)))))/((1 + y*(-1 + z))*(1 + y*(-1 + z) + z))))/y) + 2*(1 + y*(-1 + z))*(polylog(3,1/(1 + z)) - polylog(3,1/(1 + z/(1 + y*(-1 + z)))));

    return result - lim12(z);
  }

  //_________________________________________________________________________________
  double lim13(double z)
  {
    return 2.*log(1-z) - 2. + Pi2/3;
  }

  //_________________________________________________________________________________
  double fun13(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      -(((2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*log((1 - y)/(1 + y*(-1 + z)))*log(1 + y*(-1 + z)))/(y*(-1 + z))) - ((2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*log((1 - y)/(1 + y*(-1 + z)))*log(y - y*z))/(1 + y*(-1 + z)) + ((2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*log(1 + y*(-1 + z))*log(y - y*z))/(1 + y*(-1 + z)) - 2*(1 + y*(-1 + z))*log((1 - y)/(1 + y*(-1 + z)))*log(1 + y*(-1 + z))*log(y - y*z);

    return result - lim13(z);
  }

  //_________________________________________________________________________________
  double lim14(double z)
  {
    return 0.0;
  }

  //_________________________________________________________________________________
  double fun14(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      (log(1 + y*(-1 + z))*(((2 + y*(-1 + z))*(2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*log(1 + y*(-1 + z)))/(1 + y*(-1 + z)) + (2*y*(2 + y*(-1 + z))*(2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*(-1 + z)*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))))/(1 + y*(-1 + z)) + (1 + pow(1 + y*(-1 + z),2))*(2 + y*(-1 + z))*log(1 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))) + y*(1 + pow(1 + y*(-1 + z),2))*(1 - z)*log(1 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))) + 2*y*(1 + y*(-1 + z))*(2 + y*(-1 + z))*(-1 + z)*log(1 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z)))))/pow(2 + y*(-1 + z),2);

    return result - lim14(z);
  }

  //_________________________________________________________________________________
  double lim15(double z)
  {
    return 0.0;
  }
  double fun15(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      (((2 + y*(-1 + z))*(2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*log(1 + y*(-1 + z))*log(2 + y*(-1 + z)))/(1 + y*(-1 + z)) - y*(1 + pow(1 + y*(-1 + z),2))*(1 - z)*log(1 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))) + (y*(2 + y*(-1 + z))*(2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*(-1 + z)*log(2 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))))/(1 + y*(-1 + z)) + (1 + pow(1 + y*(-1 + z),2))*(2 + y*(-1 + z))*log(1 + y*(-1 + z))*log(2 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))) + y*(1 + pow(1 + y*(-1 + z),2))*(1 - z)*log(1 + y*(-1 + z))*log(2 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))) + 2*y*(1 + y*(-1 + z))*(2 + y*(-1 + z))*(-1 + z)*log(1 + y*(-1 + z))*log(2 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))))/pow(2 + y*(-1 + z),2);

    return result - lim15(z);
  }

  //_________________________________________________________________________________
  double lim16(double z)
  {
    return -Pi2/12*log(1-z);
  }

  //_________________________________________________________________________________
  double fun16(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result =
      (-((y*(2 + y*(-1 + z))*(2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*(-1 + z)*log(2 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z))))/(1 + y*(-1 + z))) + ((2 + y*(-1 + z))*(2 + 2*y*(-1 + z) + pow(y,2)*pow(-1 + z,2))*polylog(2,-1 + y - y*z))/(1 + y*(-1 + z)) + (1 + pow(1 + y*(-1 + z),2))*(2 + y*(-1 + z))*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z)))*polylog(2,-1 + y - y*z) + y*(1 + pow(1 + y*(-1 + z),2))*(1 - z)*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z)))*polylog(2,-1 + y - y*z) + 2*y*(1 + y*(-1 + z))*(2 + y*(-1 + z))*(-1 + z)*log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z)))*polylog(2,-1 + y - y*z))/pow(2 + y*(-1 + z),2);

    return result - lim16(z);
  }

  //_________________________________________________________________________________
  double lim20(double z)
  {
    return 0.0;
  }

  //_________________________________________________________________________________
  double fun20(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = polylog(2,((-1 + y)*(-1 + z))/(1 + y*(-1 + z)))/(1 + y*(-1 + z)) + polylog(3,((-1 + y)*(-1 + z))/(1 + y*(-1 + z)));

    return result - lim20(z);
  }

  //_________________________________________________________________________________
  double lim24(double z)
  {
    return 0.0;
  }

  //_________________________________________________________________________________
  double fun24(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = (log(1 + z/(1 + y*(-1 + z)))*(-((-1 + y)*(-1 + z)*(1 + y*(-1 + z) + z)*log(1 + z/(1 + y*(-1 + z)))) + z*log(z/(1 + y*(-1 + z)))*(2*(-1 + y + z - y*z) + (1 + y*(-1 + z) + z)*log(1 + z/(1 + y*(-1 + z))))))/(pow(1 + y*(-1 + z),2)*z*(1 + y*(-1 + z) + z));

    return result - lim24(z);
  }

  //_________________________________________________________________________________
  double lim25(double z)
  {
    return -Pi2/12 * log(2.);
  }

  //_________________________________________________________________________________
  double fun25(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = (((-1 + y)*(-1 + z)*pow(log(1 + z/(1 + y*(-1 + z))),2))/z + ((-1 + y + z - y*z + (1 + y*(-1 + z) + z)*log(1 + z/(1 + y*(-1 + z))))*polylog(2,-(z/(1 + y*(-1 + z)))))/(1 + y*(-1 + z) + z))/pow(1 + y*(-1 + z),2);

    return result - lim25(z);
  }

  //_________________________________________________________________________________
  double lim30(double z)
  {
    return Pi2/12*log(2.) - 0.5*log(2.)*log(2.)*log(2.);
  }

  //_________________________________________________________________________________
  double fun30(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = ((-1 + y + z - y*z)*log(z/(1 + y*(-1 + z) + z))*log(1 + z/(1 + y*(-1 + z))) + (-1 + y + z - y*z + (1 + y*(-1 + z) + z)*log(1 + z/(1 + y*(-1 + z))))*polylog(2,1/(1 + z/(1 + y*(-1 + z)))))/(pow(1 + y*(-1 + z),2)*(1 + y*(-1 + z) + z));

    return result - lim30(z);
  }

  //_________________________________________________________________________________
  double lim32(double z)
  {
    return 0.0;
  }

  //_________________________________________________________________________________
  double fun32(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = (polylog(2,((-1 + y)*(-1 + z))/(1 + y*(-1 + z))) + polylog(3,((-1 + y)*(-1 + z))/(1 + y*(-1 + z))))/pow(1 + y*(-1 + z),2);

    return result - lim32(z);
  }

  //_________________________________________________________________________________
  double lim34(double z)
  {
    return 1./6*log(2.)*log(2.)*log(2.) - Pi2/12*log(2.) + 7./8*zeta3;
  }

  //_________________________________________________________________________________
  double fun34(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = ((-1 + y)*(-1 + z)*polylog(2,1/(1 + z/(1 + y*(-1 + z)))))/((1 + y*(-1 + z))*(1 + y*(-1 + z) + z)) + polylog(3,1/(1 + z/(1 + y*(-1 + z))));

    return result - lim34(z);
  }

  //_________________________________________________________________________________
  double lim35(double z)
  {
    return 1./6*log(2.)*log(2.)*log(2.) - Pi2/12*log(2.) + 7./8*zeta3;
  }

  //_________________________________________________________________________________
  double fun35(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = ((-1 + y)*(-1 + z)*polylog(2,1/(1 + z/(1 + y*(-1 + z)))) + (1 + y*(-1 + z) + z)*polylog(3,1/(1 + z/(1 + y*(-1 + z)))))/(pow(1 + y*(-1 + z),2)*(1 + y*(-1 + z) + z));

    return result - lim35(z);
  }

  //_________________________________________________________________________________
  double lim37(double z)
  {
    return pow(log(1-z),3) - Pi2/3*log(1-z) + 2.*zeta3;
  }

  //_________________________________________________________________________________
  double fun37(double y, void * params)
  {
    double z = *((double *) params);
    double result;

    result = (log(y - y*z)*(log(y - y*z) + log(((-1 + y)*(-1 + z))/(1 + y*(-1 + z)))*(2 + 2*y*(-1 + z) + log(y - y*z))))/pow(1 + y*(-1 + z),2);

    return result - lim37(z);
  }

  //_________________________________________________________________________________
  double limNS(double z)
  {
    return -2./3*Pi2*log(1-z) + 4./3*Pi2 + 10.*log(2.)*log(2.);
  }

  //_________________________________________________________________________________
  double funNS(double y, void * params)
  {
    struct int_params *p = (struct int_params*) params;
    double z = p->z;

    double result;

    result = 4.*fun1(y,&z) + 4.*fun2(y,&z) + 4.*fun3(y,&z)
             + 2.*fun4(y,&z) + 4.*fun5(y,&z)
             + 4.*fun8(y,&z) + 4.*fun9(y,&z) - 4.*fun10(y,&z) - 4.*fun11(y,&z)
             + 8.*fun12(y,&z) - 4.*fun13(y,&z) - 2.*fun14(y,&z) + 8.*fun15(y,&z)
             + 8.*fun16(y,&z);

    return result;
  }

  //_________________________________________________________________________________
  double limS(double z)
  {
    return 2./3*Pi2*log(1-z) + 4.*Pi2 - 10.*log(2.)*log(2.);
  }

  //_________________________________________________________________________________
  double funS(double y, void * params)
  {
    struct int_params *p = (struct int_params*) params;
    double z = p->z;
    double nl = p->nl;

    double result;

    result = 4.*fun1(y,&z) + 4.*fun2(y,&z) + 4.*fun3(y,&z) - 2.*fun4(y,&z) - 4.*fun5(y,&z)
             - 4.*fun8(y,&z) - 4.*fun9(y,&z) + 4.*fun10(y,&z) + 4.*fun11(y,&z)
             - 8.*fun12(y,&z) - 4.*fun13(y,&z) + 2.*fun14(y,&z) - 8.*fun15(y,&z) - 8.*fun16(y,&z)
             - 24.*nl*fun20(y,&z);

    return result;
  }

  //_________________________________________________________________________________
  double limG(double z)
  {
    return -4.*log(1-z)*log(1-z)*log(1-z) + (4./3)*Pi2*log(1-z)
           + (4./3)*Pi2*log(2.) - 4.*log(2.)*log(2.)*log(2.) - 8.*zeta3;
  }

  //_________________________________________________________________________________
  double funG(double y, void * params)
  {
    struct int_params *p = (struct int_params*) params;
    double z = p->z;

    double result;

    result = -8.*fun24(y,&z) - 8.*fun25(y,&z) + 8.*fun30(y,&z)
             + 8.*fun32(y,&z) + 16.*fun34(y,&z) - 16.*fun35(y,&z) - 4.*fun37(y,&z);

    return result;
  }

  //_________________________________________________________________________________
  NumericIntegrals::NumericIntegrals(int const& iflv, int const& nl, double const& epsabs, double const& epsrel):
    _iflv(iflv),
    _nl(nl),
    _epsabs(epsabs),
    _epsrel(epsabs),
    _w(gsl_integration_workspace_alloc(10000))
  {
  }

  //_________________________________________________________________________________
  NumericIntegrals::~NumericIntegrals()
  {
    gsl_integration_workspace_free(_w);
  }

  //_________________________________________________________________________________
  double NumericIntegrals::integrate(double const& z) const
  {
    double (*fun)(double, void*);

    switch(_iflv)
      {
      case 0:
        fun = funS;
        break;
      case 1:
        fun = funG;
        break;
      case 2:
        fun = funNS;
        break;
      default:
        throw std::runtime_error("[NumericIntegrals::integrate]: unknown flavour");
      }

    int_params p { z, _nl };

    // std::cout << "numNS : " << funNS(0.5, &p) << std::endl;
    // std::cout << "numS : " << funS(0.5, &p) << std::endl;

    gsl_function F;
    F.function = fun;
    F.params = &p;

    double result, error;
    gsl_integration_qag(&F, 0.0, 1.0, _epsabs, _epsrel, 10000, 4, _w, &result, &error);

    return result;
  }
}
