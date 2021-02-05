/**
 *	@file		src/include/PreliminaryFunctions.h
 *	@date		22/03/2018
 *	@author	ATroja
 */

#ifndef PRELIMINARYFUNCTIONS_H
#define	PRELIMINARYFUNCTIONS_H

#include <iostream>
#include <complex>

#include "Cosmology.h"

void complexpow (double x, std::complex<double> y, std::complex<double>& z);

// ----- Elementary Functions ----- //

class ElementaryFunctions {

	public:

		/**
		 *	@brief	Default empty Constructor
		 */
		ElementaryFunctions () = default;

		/**
		 *	@brief	Euler Gamma function
		 *	@param	z, the value at which evaluate the gamma
		 *	@return Gamma(z)
		 */
		double GammaC (double z);
		/**
		 *	@brief	Euler Gamma function for complex variables
		 *	@param	z, the value at which evaluate the gamma; result, the complex value of Gamma (z)
		 */
		void GammaC (std::complex<double> z, std::complex<double>& result);

		/**
		 *	@brief	Ratio of Euler Gamma functions
		 *	@param	z1, z2, the values at which evaluate the gammas
		 *	@return Gamma(z1)/Gamma (z2)
		 */
		double GammaRatioC (double z1, double z2);
		/**
		 *	@brief	Ratio of Euler Gamma functions for complex variables
		 *	@param	z1, z2, the values at which evaluate the gammas; result, the complex value of Gamma(z1)/Gamma (z2)
		 */
		void GammaRatioC (std::complex<double> z1, std::complex<double> z2, std::complex<double>& result);

		/**
		 *	@brief	Hypergeometric Functions
		 *	@param	a,b,c,z, parameters required by 2F1
		 *	@return 2F1 (a,b,c;z)
		 */
		double Hyp2F1basic (double a, double b, double c, double z);
		/**
		 *	@brief	Hypergeometric Functions for complex variables
		 *	@param	a,b,c,z, parameters required by 2F1, s the complex value of 2F1 (a,b,c;z)
		 */
		void Hyp2F1basic (std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> z, std::complex<double>& s);

		/**
		 *	@brief	Minimum (a,b)
		 *	@params	a,b variables to compare
		 *	@return The minimum between a and b
		 */
		double MinC (double a, double b);

		/**
		 *	@brief	Maximum (a,b)
		 *	@params	a,b variables to compare
		 *	@return The maximum between a and b
		 */
		double MaxC (double a, double b);


		/**
		 *	@brief	Empty destructor
		 */
		~ElementaryFunctions () {};

};

//  ----- Special Functions ----- //

class SpecialFunctions {

	// ----- Function needed for evaluating I_l (t, nu), eq. (2.19) ----- //
	// ----- See Appendix B for further details ----- //

	public:

/*
		std::complex <double> mickey;
    std::complex <double> goofy;
    std::complex <double> donald;
    std::complex <double> cz;
    std::complex <double> cp;

		std::complex <double> huey;
    std::complex <double> dewey;
    std::complex <double> louie;
    std::complex <double> scrooge;
    std::complex <double> gladstone;
    std::complex<double> carg;
		*/

		SpecialFunctions () = default;

		/**
		 *	@brief	Evaluate the value tminC such that |I_l (nu, tminC)| = 10^-5 |I_l (nu, 1)|
		 *	@param	l, the index of the integral, nu the argument of the integral
		 *	@return	tminC
		 */
		double tminC (int l, std::complex<double> nu);

		/**
		 *	@brief	Evaluation of I_l (nu,t), eq (2.19)
		 *	@param	l, the index of the integral 
		 *	@param	nu, the argument of the integral
		 *	@param	t, the argument of the integral
		 *	@param	Il, the value of the integral
		 */
		void Il (int l, std::complex<double> nu, double t, std::complex <double>& Il);

		~SpecialFunctions () {};

};

// ----- Fourier transform class ----- //
class FourierTransform:public LinearMatterPowerSpectrum {

	private:

		double delta;
		int* n;
		int* m;
		double* kn;
		std::complex <double>* etam;
		std::complex <double>* etan;
		double** Pn;
		std::complex <double>* cn;
		std::complex <double>* cnsym;
		
		int nmax;


	public:

		FourierTransform (std::string pfile):LinearMatterPowerSpectrum (pfile) {};

		void CoeffTransfer (double (&P) (double, double), double b, double cst, int Nmax, double kmin, double kmax, std::complex<double>* x, std::complex<double>* y);

		~FourierTransform () {};
};
#endif
