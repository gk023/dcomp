//  file: dcomp_main.cpp
//
//  Main program for Dispersive Corrections to Optical Model Parameters
//  (DCOMP) code. 
//  
//  Depends on:
//
//    OMP.h -- Header file to generate optical model parameters
//    
//    Coulomb.h -- Finds the Coulomb correction term for the potential
//
//    gsl/gsl_integration.h -- GSL Integration class to solve dispersive correction
//
//    boost/numer/odeint.hpp -- Boost library used to handle complex differential equation
//
//    Asym.h -- Header file I wrote to get Hankel functions and their derivatives to
//               find S-Matrix elements
//
//   DCOMP takes a card file with a Fixed real central potential and parameters for an
//    imaginary volume potential parameterized like in 'A.J. Koning and J.P. Delaroche
//    Nuc. Phys. A, 713(3-4), (2003).' The card takes in the target and projectile A and Z,
//    the parameters for the fixed real and KD imaginary potential, the matching radius to
//    check asymptotics, and the maximum angular momentum, the number of integration points 
//    for the 4th order Runge-Kutta method employed with odeint, and the errors on the GSL
//    integration for the dispersive corrections. DCOMP uses a python script to run the code
//    and produce plots of the output. The usage is as follows:
//         
//          python run_dcomp.py <path/to/card/file>
//
// *********************************************************************

//Include Headers
#include <cmath>
#include <complex>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <gsl/gsl_integration.h>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

#include "OMP.h" //Optical potential class to find Woods-Saxon
#include "Coulomb.h" //Finds Coulomb correction term
#include "Asym.h" //Class to find Hankel functions and their derivative


typedef std::vector< complex<double> > state_type; //How the wave function will be stored


//**********************************************************************

//Function protoypes

double Wv_E( double E, void *params_ptr );
double doublefactorial(double L);
inline double sqr(double x) {return x*x;};

//***********************************************************************

//Constants and structures

const double pi = 3.141592653589793; // define pi
const double u = 0.02392922805; // atomic mass unit in MeV/(hbar*c)^2
const double esq = 1.4401459854;

//Structure for W_vol parameters
typedef struct
{
   double w1;
   double w2;
   double E_f;
   double n;
}
im_vol_params;

//odeint structure for the rhs of the omp equation
struct omp_rhs
{


   double E; //beam energy
   double L; //Angular momentum
   double mu; //reduced mass
   double V;
   double r;
   double a;
   double W;
   double rw;
   double aw;
   double At;
   double Ap;
   double Zt;
   double Zp;
   double Rc;
   double V_disp;

   omp_rhs(double eng, double orb, double Vr, double rr, double ar, double Wv, double rv, double av, double A_tgt, double A_proj, double Z_tgt, double Z_proj, double rad_c, double V_corr)
   {
    
       E= eng;
       L= orb;
       V = Vr;
       r = rr; 
       a = ar; 
       W = Wv;
       rw = rv;
       aw = av;
       At = A_tgt;
       Ap = A_proj;
       Zt = Z_tgt;
       Zp = Z_proj;
       Rc = rad_c;
       V_disp = V_corr;
       mu = u*(At*Ap)/(At+Ap); //reduced mass per hbar*c
     
    }

   
   void operator() ( const state_type &y, state_type &dydx, const double x) const 
   {
       
       OMP initial_pot(V,r,a,W,rw,aw,At);
       OMP dispersive_pot(V_disp,rw,aw,0.,0.,0.,At);
       Coulomb charge_corr(Zt,Zp,Rc);
      

       dydx[0] = y[1];
       dydx[1] = y[0]*( L*(L+1.)/( sqr(x) ) + 2.*mu*( initial_pot.pot(x) + dispersive_pot.pot(x) + charge_corr.pot(x) - E) );
   }


};


//Will print out wave functions for the runs
struct streaming_observer
{

  std::ostream& m_out;
  streaming_observer ( std::ostream &out) : m_out( out ){}

  template< class State >
  void operator()(const State &x, double t) const
  {

     complex<double> x0 = x[0];
     m_out << t << "          " << x0.real()
           << "          " << x0.imag() << endl;

  }

};

//***********************************************************************

int main()
{
   //Define the RK4 stepper for the complex vector used to store chi_L, chi_L'
   typedef runge_kutta4< state_type > stepper_type;

   cout << " *************************************" << endl;
   cout << " ****** Running DCOMP version 1 ******" << endl;
   cout << " *************************************" << endl;

   double Z_t, A_t, Z_p, A_p, w1, w2, Ef_plus, Ef_minus, n, Ef, rw, aw, Vhf, rhf, ahf, 
          rc, rmax, Lmax, Npts, int_abs_err, int_rel_err;

   //beam energy, wavevector, elastic and reaction cross sections
   double Ebeam, k, cross_section, rxn_cross_section;

   //scattering amplitude scaled by a constant
   complex<double> amp;

   // Read in all the parameters: Z target, A target, Z projectile, 
   // A projectile, w1, w2, Fermi energy, imaginary radius, imaginary diffuseness
   // Hartree-Fock volume, radius, diffuseness, Coulomb radius, matching radius,
   // mesh points, integration absolute error, and integration relative error.
   // Python script provides all the standard input for the code from the card file
   cout << "Running with the following parameters:" << endl;
   cin >> Z_t;
   cout << "Target Z--> " << Z_t << endl;
   cin >> A_t;
   cout << "Target A--> " << A_t << endl;
   cin >> Z_p;
   cout << "Projectile Z--> " << Z_p << endl;
   cin >> A_p;
   cout << "Projectile A--> " << A_p << endl;
   cin >> w1;
   cout << "Imaginary volume w1--> " << w1 << endl;
   cin >> w2;
   cout << "Imaginary volume w2--> " << w2 << endl;
   cin >> Ef_plus;
   cout << "E_F+--> " << Ef_plus << endl;
   cin >> Ef_minus;
   cout << "E_F- --> " << Ef_minus << endl;
   Ef = 0.5*(Ef_plus-Ef_minus);
   cout << "Average Fermi energy--> " << Ef << endl;
   cin >> n;
   cout << "Power of Brown-Rho shape--> " << n << endl;
   cin >> rw;
   cout << "Imaginary volume radius--> " << rw << endl;
   cin >> aw;
   cout << "Imaginary volume diffuseness--> " << aw << endl;
   cin >> Vhf;
   cout << "Hartree-Fock real Volume--> " << Vhf << endl;
   cin >> rhf;
   cout << "Hartree-Fock real radius--> " << rhf << endl;
   cin >> ahf;
   cout << "Hartree-Fock real diffuseness--> " << ahf << endl;
   cin >> rc;
   cout << "Coulomb Radius--> " << rc << endl;
   cin >> rmax;
   cout << "Matching Radius--> " << rmax << endl;
   cin >> Lmax;
   cout << "Max Angular Momentum--> " << Lmax << endl;
   cin >> Npts;
   cout << "Mesh Points--> " << Npts << endl;
   cin >> int_abs_err;
   cout << "GSL absolute integration error--> " << int_abs_err << endl;
   cin >> int_rel_err;
   cout << "GSL relative error--> " << int_rel_err << endl;
   cout << "*********************************************" << endl;
   //Set everthing up to begin calculating the dispersive corrections to the optical potential

   gsl_integration_workspace *work_ptr
      = gsl_integration_workspace_alloc (1000);
    
   double lower_limit = -10e8;
   double upper_limit = 10e8;
   
   //Get everything pointing to the right place for GSL integration
   im_vol_params Wv_E_params;
   Wv_E_params.w1 = w1;
   Wv_E_params.w2 = w2;
   Wv_E_params.E_f = Ef;
   Wv_E_params.n = n;

   gsl_function Wv_E_integrand;
   void *params_ptr = &Wv_E_params;

   Wv_E_integrand.function = &Wv_E;
   Wv_E_integrand.params = params_ptr;
   
   //Start an output file for the cross sections
   ofstream xs_out;
   xs_out.open("cross_sections.out");
   xs_out << "#E     Elastic XS       Rxn XS         Total XS" << endl;

   double mu = (A_t*A_p)/(A_t+A_p)*u;
   double Emin=5;
   double Emax=200;
   double Estep=5;
   for(double E=Emin; E <= Emax; E += Estep)
   {
      Ebeam = E;
      k = sqrt(2*mu*E);
     
      cout << "\n Energy = " << Ebeam;
     
      //Initialize results and errors for dispersive correction integral
      double I1, I2, I3, I4;
      double E1, E2, E3, E4;   
   
      //Integrates each of the four terms of the principal value integral for the 
      //dispersivd correction
             
      gsl_integration_qawc(&Wv_E_integrand, Ef_plus, upper_limit, Ebeam, int_abs_err,
                             int_rel_err, 1000, work_ptr, &I1, &E1);
   
      gsl_integration_qawc(&Wv_E_integrand,Ef_plus, upper_limit, Ef, int_abs_err,
                             int_rel_err, 1000, work_ptr, &I2, &E2);
      
      gsl_integration_qawc(&Wv_E_integrand, lower_limit, Ef_minus, Ef, int_abs_err,
                             int_rel_err, 1000, work_ptr, &I3, &E3);
      
      gsl_integration_qawc(&Wv_E_integrand, lower_limit, Ef_minus, Ebeam, int_abs_err,
                             int_rel_err, 1000, work_ptr, &I4, &E4);
          
      // Add up principal values, get dispersive correction
      double V_disp = ( (-1.0*I1) + I2 + (-1.0*I3) + I4 )/pi;
   
      cout << "\n\n Disperive correction volume term: " << V_disp << endl;
   
      //******************************************************************************************
      //Finding the partial waves starts here
   
      //Two dimensional complex vector, row 0 -> wavefunction, row 1 -> derivative
      state_type y(2);
      
      cross_section = 0.;
      rxn_cross_section = 0.;
   
      complex<double> RL, SL; //R- and S-Matrix elements
      //Loop over all L to get the partial waves
      for(double L=0.; L <= Lmax; L += 1.0)
      {
         //Start and output file for the partial waves

//************************ uncomment these lines to print waves ***************************
         //ostringstream L_wave_stream;
         //L_wave_stream << "wf_E_" << Ebeam << "_" << L << ".out";
         //string L_wave_file = L_wave_stream.str();
         //ofstream L_out;
         //L_out.open(L_wave_file.c_str());
   
         //L_out << "#r       Re(chi)       Im(chi)" << endl;
//*****************************************************************************************
// Must also uncomment lines 371 and 391, comment out line 374 to print wave functions 
    
         //Start a bit off zero so the centrifugal barrier doesn't blow the problem up
         double rmin;
         rmin = (2.*L*rmax)/(Npts+2.*L);
         //Zero will blow it up still, so force to start higher
         if(rmin == 0.)
         { 
         
           rmin=(2.*rmax)/(Npts+2.);    
   
         }
         

         //Since we aren't exactly at R=0, give the wave functions a small value at rmin
         //This choice was taken from "Nuclear Reactions for Astrophysics", I.J. Thompson
         // and F.M. Nunes (2009), Eq. 6.3.1
         double ic = pow( (k*rmin), (L+1.) ) / ( doublefactorial(2.*L+1.) * Npts );
         y[0]= complex<double> (ic,ic);
         y[1]= complex<double> (1.,1.); // Guess something, will be right up to normalization    
      
         double h;
         h=(rmax-rmin)/(Npts); //Step size
      
         // Solve the differential equation for the Lth partial wave

         //integrate with printing wave functions -- (Lmax + 1) * (# of energies) files written 
         //integrate_adaptive( stepper_type(), omp_rhs(Ebeam, L, Vhf, rhf, ahf, Wv_E(Ebeam, &Wv_E_params), rw, aw, A_t, A_p, Z_t, Z_p, rc, V_disp), y, rmin, rmax, h, streaming_observer( L_out ) ); 

         //integrate without printing wave functions
         integrate_adaptive( stepper_type(), omp_rhs(Ebeam, L, Vhf, rhf, ahf, Wv_E(Ebeam, &Wv_E_params), rw, aw, A_t, A_p, Z_t, Z_p, rc, V_disp), y, rmin, rmax, h ); 
   
         
         RL=y[0]/(rmax*y[1]);  //R-matrix element
         cout << "\n R-Matrix element for wave " << L << ": " << RL << endl;
       
         double eta = esq*Z_t*Z_p*k/(2.*E); //Sommerfeld parameter
        
         Asym HL(eta,k,k*rmax,L); //find the asymptotics with arugments eta, k*rmax
         
         // Now get the S-matrix element with the asymptotics
         SL = (HL.Hminus() - rmax*RL*HL.Hminus_prime())/(HL.Hplus() - rmax*RL*HL.Hplus_prime());

         cout << "\n S-Matrix element for wave " << L << ": " << SL << endl;
   
         cout << "\n |S| for wave " << L << ": " << abs(SL) << endl;
   
         //L_out.close(); //uncomment if you are printing wave functions
   
         amp = 1. - SL;
      
         //cross sections, scale by 10 to go from fm^2 to mb//
         cross_section += 10.*pi/(sqr(k))*(2.*L+1.)*norm(amp);
         rxn_cross_section += 10.*pi/(sqr(k))*(2.*L+1.)*(1.-norm(SL));      
   
      }

      cout << "\n The cross section at E=" << Ebeam << ": " << cross_section << endl;
      cout << "The reaction cross section at E=" << Ebeam << ": " << rxn_cross_section << endl;
      cout << "*********************************************************" << endl;

      xs_out << setprecision(10) << Ebeam << "         " << cross_section
             << "         " << rxn_cross_section 
             << "         " << cross_section + rxn_cross_section << endl;
   }
   xs_out.close();

   cout << "\n Wave functions in wf*.out" << endl;
   cout << "Cross sections as a function of beam energy in cross_sections.out" << endl;
  
   return (0);

}

//***********************************************************************

//Energy dependence of the imaginary volume, integrand for GSL QAWC
double Wv_E( double E, void *params_ptr)
{

  im_vol_params *passed_ptr;
  passed_ptr = (im_vol_params *) params_ptr;
  double a = passed_ptr->w1;
  double b = passed_ptr->w2; 
  double c = passed_ptr->E_f;
  double d = passed_ptr->n;

  return a * pow((E - c),d) / ( pow((E - c),d) + pow(b,d) );

}

//Double factorial function
double doublefactorial( double L)
{
   double result=1.;

   for(int n = ((int) L); n > 0; n--)
   {
   
      result *= ((double) n);     

   }

   return result;

}
