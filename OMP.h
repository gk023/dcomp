//  file: OMP.h
//
//  Header file for the Optical Model Potential in C++
//
//  Programmer: Garrett King
//
//*********************************************************
#include<complex>

#ifndef OMP_H
#define OMP_H

class OMP
{

   public:
      OMP (double V, double r, double a, double W, double rw, double aw, double A);
      ~OMP();
      
      std::complex<double> pot(double rad);
   private:

       double re_vol; //Real central volume
       double re_rad; //Real central radius
       double re_dif; //Real central diffuseness
       double im_vol; //Imaginary volume
       double im_rad; //Imaginary radius
       double im_dif; //Imaginary diffuseness
       double mass_num; //Mass number of the nucleus
};

#endif
