// file: Coulomb.cpp
//
// Function definitions for the Coulomb correction class
//
// Programmer: Garrett King
//
//*********************************************************

#include<complex>
#include "Coulomb.h"


const double esq = 1.4401459854; // e^2 / (4*pi*epsilon_0) in MeV*fm

Coulomb::Coulomb(double charge_t, double charge_p, double rad_c)
{

   Zt=charge_t;
   Zp=charge_p;
   Rc=rad_c;
}

Coulomb::~Coulomb()
{}

std::complex<double> Coulomb::pot(double r)
{

   double result;
   
   if(Rc==0.)
   {
      result = 0.;
   }

   if ((r <= Rc) && (Rc != 0.))
   {
      
       result = Zp*Zt*esq/(2.*Rc)*(3.-(r*r)/(Rc*Rc));

   }

   if ( (r > Rc) && (Rc != 0.))
   {

       result = Zp*Zt*esq/r;  

   }

  return std::complex<double> (result,0.);

}
