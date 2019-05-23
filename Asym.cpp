// file Asym.cpp
//
// Definitions for the Asym class
//
//**********************************

#include <complex>
#include <cmath>
#include "Asym.h"

#include<gsl/gsl_sf_coulomb.h>

//Returns the asymptotics needed to calculate the S-matrix

Asym::Asym(double sommerfeld_parm, double wave_vect, double a, double orb)
{
  rmatch = a;
  L = orb;
  eta = sommerfeld_parm;
  k = wave_vect;
}

Asym::~Asym()
{ }

std::complex<double> Asym::Hplus()
{
  gsl_sf_result F, G, Fp, Gp;
  double exp_F, exp_G;  

  gsl_sf_coulomb_wave_FG_e(eta,rmatch,L,0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
  
  return std::complex<double> (G.val*exp(exp_G) , F.val*exp(exp_F));

}


std::complex<double> Asym::Hminus()
{
  gsl_sf_result F, G, Fp, Gp;
  double exp_F, exp_G;

  gsl_sf_coulomb_wave_FG_e(eta,rmatch,L,0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
  
  return std::complex<double> (G.val*exp(exp_G) , -1.*F.val*exp(exp_F));

}

std::complex<double> Asym::Hplus_prime()
{
  gsl_sf_result F, G, Fp, Gp;
  double exp_F, exp_G;  

  gsl_sf_coulomb_wave_FG_e(eta,rmatch,L,0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
  
  return std::complex<double> (k*Gp.val*exp(exp_G) , k*Fp.val*exp(exp_F));

}

std::complex<double> Asym::Hminus_prime()
{
  gsl_sf_result F, G, Fp, Gp;
  double exp_F, exp_G;

  gsl_sf_coulomb_wave_FG_e(eta,rmatch,L,0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
  
  return std::complex<double> (k*Gp.val*exp(exp_G), -1.*k*Fp.val*exp(exp_F));

}

