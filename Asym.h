// file: Asym.h
//
// Header file to call for asymptotics in scattering
// Claculations
//
// Programmer: Garrett King
//
//**************************************************

#ifndef ASYM_H
#define ASYM_H

#include <complex>
#include <gsl/gsl_sf_coulomb.h> //Gets GL, FL the irregular and regular Coulomb wave functions

class Asym
{

   public:
      Asym(double sommerfeld_parm, double wave_vect, double a, double orb);
      ~Asym();

      std::complex<double> Hplus();
      std::complex<double> Hminus();
      std::complex<double> Hplus_prime();
      std::complex<double> Hminus_prime();
 
   private:
   
      double rmatch;
      double L;
      double eta;
      double k;
};

#endif
