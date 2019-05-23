// file OMP.cpp
// 
// Definitions for the OMP class
//
// Programmer: Garrett King
//
//*****************************************

#include <complex>
#include <cmath>
#include "OMP.h"

//Get optical potential with parameters V, r, a, W, rw, aw for a 
//nucleus with A nucleons

OMP::OMP(double V, double r, double a, double W, double rw, double aw, double A)
{
   re_vol = V;
   re_rad = r;
   re_dif = a;
   im_vol = W;
   im_rad = rw;
   im_dif = aw;
   mass_num = A;

}

OMP::~OMP()
{ }

std::complex<double> OMP::pot(double rad)
{
   double Vv;
   double Wv;

   Vv = -1. * re_vol/(1.+exp( (rad - pow(mass_num, 1./3.)*re_rad) / re_dif ) );
   Wv = -1. * im_vol/(1.+exp( (rad - pow(mass_num, 1./3.)*im_rad) / im_dif ) );
   
   return std::complex<double> (Vv,Wv);
}


