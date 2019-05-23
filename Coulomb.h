// file: Coulomb.h
//
// Header file for the Coulomb correction class
//
// Programmer: Garrett King
//
//*********************************************************

#include<complex>

#ifndef COULOMB_H
#define COULOMB_H

class Coulomb
{

   public:
      Coulomb (double charge_t, double charge_p, double rad_c);
      ~Coulomb();

      std::complex<double> pot(double r);

   private:

      double Zt;
      double Zp;
      double Rc;


};



#endif
