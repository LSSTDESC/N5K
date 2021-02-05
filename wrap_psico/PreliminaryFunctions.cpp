/**
 *	@file		src/mods/PreliminaryFunctions.cpp
 *	@date		22/03/2018
 *	@author	ATroja
 */

#include "PF.h"
//#include "../include/FileUtils.h"

#include "gsl/gsl_math.h"
#include "gsl/gsl_deriv.h"
#include "gsl/gsl_spline.h"

#include <memory>

void GammaC (std::complex<double> z, std::complex<double>& result) {

	double q0 = 75122.6331530;
	double q1 = 80916.6278952;
	double q2 = 36308.2951477;
	double q3 = 8687.24529705;
	double q4 = 1168.92649479;
	double q5 = 83.8676043424;
	double q6 = 2.50662827511;
	std::complex <double> p1, p2;
	if (z.real () >= 0.) {

		p1 = ( q0 
				 + q1 * z
				 + q2 * z * z
				 + q3 * z * z * z
				 + q4 * z * z * z * z
				 + q5 * z * z * z * z * z
				 + q6 * z * z * z * z * z * z
				 ) /
				 ( z
			   * (z + 1.)
			   * (z + 2.)
			   * (z + 3.)
			   * (z + 4.)
			   * (z + 5.)
			   * (z + 6.)
				 );

		result = p1 * pow ((z + 5.5), (z + 0.5)) * exp (-z-5.5);
	}
	else {

		p1 = ( q0 
				 + q1 * (1. - z)
				 + q2 * (1. - z) * (1. - z)
				 + q3 * (1. - z) * (1. - z) * (1. - z)
				 + q4 * (1. - z) * (1. - z) * (1. - z) * (1. - z)
				 + q5 * (1. - z) * (1. - z) * (1. - z) * (1. - z) * (1. - z)
				 + q6 * (1. - z) * (1. - z) * (1. - z) * (1. - z) * (1. - z) * (1. - z)
				 ) /
				 ( (1. - z)
			   * (2. - z)
			   * (3. - z)
			   * (4. - z)
			   * (5. - z)
			   * (6. - z)
			   * (7. - z)
				 );

		p2 = p1 * pow ((1. - z + 5.5), (1. - z + 0.5)) * exp (-1. + z - 5.5);
		
		result = M_PI / sin (M_PI * z) / p2;
	}
}

void GammaRatioC (std::complex <double> z1, std::complex <double> z2, std::complex <double>& result) {

	result.real (1.);
	result.imag (0.);
	double q0 = 0.0075122633153;
	double q1 = 0.0809166278952;
	double q2 = 0.363082951477;
	double q3 = 0.868724529705;
	double q4 = 1.16892649479;
	double q5 = 0.838676043424;
	double q6 = 0.250662827511;

	for (int i = 0; i < 7; i++) {
		result *= (z2 + i*1.) / (z1 + i*1.);
	}

	result *= exp (z2 - z1);

	std::complex<double> p1, p2;
	p1 = (z1 + 0.5) * log (z1 + 5.5) - (z2 + 0.5) * log (z2 + 5.5);
	result *= exp (p1);

	std::complex <double> z1t = z1 / 10.;
	std::complex <double> z2t = z2 / 10.;

	p1 = q0
		 + q1 * z1t
		 + q2 * z1t * z1t
		 + q3 * z1t * z1t * z1t
		 + q4 * z1t * z1t * z1t * z1t
		 + q5 * z1t * z1t * z1t * z1t * z1t
		 + q6 * z1t * z1t * z1t * z1t * z1t * z1t;
	p2 = q0
		 + q1 * z2t
		 + q2 * z2t * z2t
		 + q3 * z2t * z2t * z2t
		 + q4 * z2t * z2t * z2t * z2t
		 + q5 * z2t * z2t * z2t * z2t * z2t
		 + q6 * z2t * z2t * z2t * z2t * z2t * z2t;

	result *= p1 / p2;

}

void Hyp2F1basic (std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> z, std::complex<double>& s) {

	// ----- Power series representation of 2F1 evaluated around z = 0 ----- //

	s.real (0.);
	s.imag (0.);
	std::complex<double> p (1., 0.);
	double eps = 1.;
	double n = 0.;
	std::complex<double> sold;

	while (eps > 0.000000000000001) {

		sold = s;
		s = s + p;
		p *= (a+n) * (b+n) / (c+n) / (n+1.) * z;
		eps = fabs ((s-sold)/s);
		n++;
	
	}
}

double tminC (int lin, std::complex<double> nu) {

	// ----- I_l (nu, t) is evaluated in the interval [tminC, 1/tminC] and is set to 0 outside this interval ----- //

	double tmin;
	double l = lin * 1.;

	if (l > 60) {

		tmin = ( ( 1. - exp (-0.08405923801793776 * l) ) * pow (l, 1.0388189966482335) )
				 / ( 16.552260860083162 + pow (l, 1.0388189966482335) ) 
				 + ( 1. / 72. )
				 * ( 
					 - ( ( 1. - exp (-0.08405923801793776 * l) ) * pow (l, 1.0388189966482335) )  
					 / ( 16.552260860083162 + pow (l, 1.0388189966482335) )
					 + ( ( 1. - exp (-0.03269491513404876 * l) ) * pow (l, 1.0606484271153198) )
					 / ( 86.60472131391394 + pow (l, 1.0606484271153198) ) 
					 ) 
				 * fabs (nu.imag ());
	}
	else {

		double mickey	= 0.026189266734847335
			   					- 0.04199333649354753 		* l
									+ 0.01813725076906472 		* l * l
				 					- 0.0019512662766506912 	* l * l * l
				 					+ 0.00011476285973931163 	* l * l * l * l
				 					- 4.132495310006262e-6 		* l * l * l * l * l
				 					+ 9.321216021016041e-8 		* l * l * l * l * l * l
				 					- 1.2840836892275476e-9 	* l * l * l * l * l * l * l
				 					+ 9.874351126449866e-12 	* l * l * l * l * l * l * l * l
				 					- 3.247034458438453e-14 	* l * l * l * l * l * l * l * l * l ;
		
		double goofy  = 0.0050534423514964006		
									- 0.004245361441903382		* l
									+ 0.0009644735508629553		* l * l
									- 0.000029194973960888548	* l * l * l
								  - 1.197268126576586e-7		* l * l * l * l
									+ 3.9232441693781885e-8		* l * l * l * l * l
									- 1.3806236786152843e-9		* l * l * l * l * l * l
									+ 2.380296810916152e-11		* l * l * l * l * l * l * l
									- 2.105287890873389e-13		* l * l * l * l * l * l * l * l
									+ 7.627228092016026e-16		* l * l * l * l * l * l * l * l * l;
									

		double donald = 0.014502978209351904		
									- 0.01218174975881159			* l
									+ 0.002817970220966578		* l * l
									- 0.00011942831975390713	* l * l * l
									+ 1.223432213234367e-6		* l * l * l * l
									+ 7.921224590247682e-8		* l * l * l * l * l
									- 3.5781997384264233e-9		* l * l * l * l * l * l
									+ 6.634225862490053e-11		* l * l * l * l * l * l * l
									- 6.057230587166174e-13		* l * l * l * l * l * l * l * l
									+ 2.230575708513619e-15		* l * l * l * l * l * l * l * l * l;

		tmin = mickey
				 + ( 1. / 91512. )
				 * ( 
					 - 4223. * mickey 
					 - 961.  * goofy
					 + 5184. * donald
					 )
				 * fabs (nu.imag ())
				 + ( 1. / 91512. )
				 * (
/*					   41. * ( mickey 
							     + 0.0019512662766506912 * l * l * l 
									 - 2.230575708513619 * 0.0019512662766506912 * l * l * l 
									 )
									 */
						 41. * mickey
					 + 31. * goofy
					 - 72. * donald
					 )
			   * fabs (nu.imag ()) * fabs (nu.imag ());

	}

	return tmin;
}

void complexpow (double x, std::complex<double> y, std::complex<double>& z) {

	double x2y = pow (x, y.real ());
	double lnx = log (x);

	z.real (x2y * cos (y.imag () * lnx));
	z.imag (x2y * sin (y.imag () * lnx));

}

void Il (int lin, std::complex<double> nu, double t, std::complex <double>& result) {

	double l = lin * 1.;
	std::complex<double> res (0., 0.);
	double arg;
	double z1min;
	double z;
	std::complex <double> factor (0., 0.);

	if ( t > 1.) {

		z = 1. / t;
		complexpow (t, -nu, factor);
		
	}
	else {

		z = t;
		factor.real (1.);

	}

	arg = (1. - z*z) * (1. - z*z) / (4. * z * z);

	if (l >= 5) {
		z1min = tminC (l, nu);
	}
	else {
		z1min = 0.;
	}

	if (z < z1min) {

		res.real (0.);
		res.imag (0.);
		
	}
	else {

		std::complex <double> mickey;
		std::complex <double> goofy;
		std::complex <double> donald;
		std::complex <double> cz (z, 0.);
		std::complex <double> cp (0., 0.);
 
		if (z <= 0.7) {

			GammaC ((3.-nu) * 0.5, mickey);
			GammaRatioC ( (l + nu * 0.5), (l + 1.5), goofy);
			Hyp2F1basic ( (nu - 1.)* 0.5, l + nu * 0.5, l + 1.5, cz * cz, donald);

			complexpow (2., nu - 1., cp);

			res = M_PI * M_PI * pow (z, l) * cp / mickey * goofy * donald;

		}


		if (z > 0.7 && z < 1.) {
			
			std::complex <double> huey;
			std::complex <double> dewey;
			std::complex <double> louie;
			std::complex <double> scrooge;
			std::complex <double> gladstone;
			std::complex<double> carg (arg, 0.);

			GammaRatioC ( (l + nu * 0.5), (l + 2. - nu * 0.5), mickey);
			GammaRatioC ( (1. - nu * 0.5), (1.5 - nu * 0.5), goofy);
			Hyp2F1basic ( (l * 0.5 + nu * 0.25), (nu * 0.25 - (l + 1.) * 0.5), nu * 0.5, -carg, donald);
			GammaC (nu - 2., huey);
			dewey = cos (M_PI * nu * 0.5);
			Hyp2F1basic ( (l * 0.5 - nu * 0.25 + 1.), (0.5 - nu * 0.25 - l * 0.5), 2. - nu * 0.5, -carg, louie);
			complexpow (z, -nu * 0.5, scrooge);
			complexpow (arg, 1. - nu * 0.5, gladstone);

			res = M_PI * scrooge
						 * ( sqrt (M_PI) * mickey * goofy * donald 
							 - 2. * huey * dewey * gladstone * louie
							 );

		}

		if (z == 1.) {

			std::complex <double> scrooge;

			GammaC (2. - nu, mickey);
			GammaC ((3. - nu) * 0.5, goofy);
			GammaRatioC ( (l + nu * 0.5), (2. + l - nu * 0.5), donald);
			complexpow (2., nu - 1., scrooge);

			res = M_PI * M_PI * pow (z, l) * scrooge * mickey / goofy / goofy * donald;

		}


	}

	result = res * factor;

}
