<%
cfg['compiler_args'] = ['-std=c++11']
cfg['sources'] = ['PreliminaryFunctions.cpp']
cfg['library_dirs'] = ['/home/iltomo/usr/local/lib/']
cfg['libraries'] = ['fftw3', 'gsl', 'gslcblas', 'm']
setup_pybind11(cfg)
%>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "PF.h"
//#include <iostream>
//#include <complex>
#include <fftw3.h>
#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"

#define REAL 0
#define IMAG 1

namespace py = pybind11;

using namespace pybind11::literals;

// Limber

struct my_f_params {int l; gsl_spline *wg; gsl_spline *plin; gsl_interp_accel *Wacc; gsl_interp_accel *Pacc;};

double rrchi (double chi, void* par) {

  struct my_f_params * params = (struct my_f_params *) par;
  int l = (params->l);
  gsl_spline *wg = (params->wg);
  gsl_spline *plin = (params->plin);
	gsl_interp_accel *Wacc = (params->Wacc);
	gsl_interp_accel *Pacc = (params->Pacc);

	//double radChi = wf.Wgalaxy_noRSD (chi, l);
 	double radChi= gsl_spline_eval (wg, chi, Wacc);
  double arg = (l+0.5)/chi;
  //return radChi * radChi * lmps.Plin (arg)  / chi / chi;
	double Plin = gsl_spline_eval (plin, arg, Pacc);
	
	/*
	if (isnan (radChi)) std::cout << "RadChi " << chi << std::endl;
	if (isnan (arg)) std::cout << "Arg " << arg << std::endl;
	if (isnan (Plin)) std::cout << "Plin " << arg << std::endl;
	*/

  return radChi * radChi * Plin / chi / chi;

}

double Cl_gal_limber_qags (int l, double chimin, double chimax, gsl_spline *spline_wg, gsl_spline *spline_plin){

  double res, err;
	gsl_interp_accel *Wacc = gsl_interp_accel_alloc ();
	gsl_interp_accel *Pacc = gsl_interp_accel_alloc ();

  gsl_function R;
  struct my_f_params alpha = {l, spline_wg, spline_plin, Wacc, Pacc};
  R.function = rrchi;
  R.params = &alpha;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
  gsl_integration_qags (&R, chimin, chimax, 0, 1e-7, 1000, 
w, &res, &err);
  gsl_integration_workspace_free (w);

  return res ;// / Dz (survey.zav) / Dz (survey.zav);

}

// Non-Limber
void herecomplexpow (double x, std::complex<double> y, std::complex<double>& z) {

  double x2y = pow (x, y.real ());
  double lnx = log (x);

  z.real (x2y * cos (y.imag () * lnx));
  z.imag (x2y * sin (y.imag () * lnx));

}


double Cl_auto_gal (int l, double chimin_photoz, double chimax_photoz, double bias, gsl_spline *spline_wg, gsl_spline *spline_plin) {

	gsl_interp_accel *Pacc = gsl_interp_accel_alloc ();
	gsl_interp_accel *Wacc = gsl_interp_accel_alloc ();

  int N = 256; //an O(n log n) algorithm is used even for prime sizes (FFTW Manual, page 7)
  double Nd = (double) N;
	double inv_Nd = 1./ Nd;
  int N2 = N / 2;

  double kmin = 1e-5;//klinsx (1e-3);
  double kmax = 1e3;//klindx (1e-3);

  double dkapa = log (kmax/kmin);
  double dkapa_over_N = dkapa / Nd;

  double k0 = kmin;//kmin;// exp ((log(kmax)+log(kmin))/2.);

  double chi0 = dkapa;

  double res1 = 0.;
  double res2 = 0.;
  std::complex <double> num (0.,0.);

  double expdN = exp (dkapa_over_N);
  double expdNm1 = expdN-1.;

  //double chimin_photoz = Chiz (survey.zmin_photoz);
  //double chimax_photoz = Chiz (survey.zmax_photoz);

  std::complex<double> goofy;

  double dt;
  double acca1 [N], acca2 [N];
  double phip = 0, phipt = 0, phipot = 0;
  double dchi;
  double tq;
  std::complex<double> ilm;

  std::complex <double> nuN2, ilN2, goofyN2, mickey, cN2, c0, donald;

  double *in1, *out2, *out3;
  fftw_complex *out1, *in2, *in3;
  in1 = (double*) fftw_malloc (sizeof (double) * (N));
  out1 = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * (N2+1));
  in2 = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * (N2+1));
  out2 = (double*) fftw_malloc (sizeof (double) * (N));
  in3 = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * (N2+1));
  out3 = (double*) fftw_malloc (sizeof (double) * (N));
 
  fftw_plan pp1, pp2, pp3;
 
  pp1 = fftw_plan_dft_r2c_1d (N, in1, out1, FFTW_ESTIMATE);
  pp2 = fftw_plan_dft_c2r_1d (N, in2, out2, FFTW_ESTIMATE);
  pp3 = fftw_plan_dft_c2r_1d (N, in3, out3, FFTW_ESTIMATE);
 
  double chi [N];
  double k [N];
  for (int n = 0; n < N; n++) {
		chi [n] = chi0 * exp (n * dkapa_over_N);
    k [n] = k0 * exp (n * dkapa_over_N);
  }
 
  double kleft = 5e-5;
  double kright = 5e2;
	double Plin;

  for (int n = 0; n < N; n++) {
		
		Plin = gsl_spline_eval (spline_plin, k[n], Pacc);

    //in1 [n] = k [n]*k[n]*k[n]*Plin (k[n]) * pow (k[n]/ k0, bias);
    //in1 [n] = k [n]*k[n]*k[n] * Plin * pow (k[n]/ k0, bias);														//QUI
    in1 [n] = k [n]*k[n]*k[n] * Plin * pow (k[n], bias);

    if (k[n] < kleft) {
      in1 [n] *= (k [n] - kmin) / (kleft - kmin) - 1./ 2. / M_PI * sin (2.*M_PI * (k[n] - kmin)/(kleft-kmin));
    }
    if (k[n] > kright) {
      in1 [n] *= (kmax - k [n]) / (kmax - kright) - 1./ 2. / M_PI * sin (2.*M_PI * (kmax - k[n])/(kmax-kright));
    }

  }


  fftw_execute (pp1);

  std::complex <double> cn [N2+1];
  std::complex <double> cn0 [N2+1];

  for (int n = 0; n <= N2; n++) {
    cn0 [n] = {out1 [n][REAL], out1 [n][IMAG]};
  }

	/*
	int Nt = 256;
  double Ntd = (double) Nt;
  for (int q=0; q < Nt; ++q){ //THE RANGE IS CORRECT IN ORDER TO GET THE RIGHT INTEGRAL
    dt = 1./Ntd;
    tq = (1.*q + 0.5)/Ntd; //IT IS CORRECT IN ORDER TO GET THE RIGHT INTEGRAL:
		*/
  for (int q=0; q < N; ++q){ //THE RANGE IS CORRECT IN ORDER TO GET THE RIGHT INTEGRAL
    tq = (1.*q + 0.5)*inv_Nd; //IT IS CORRECT IN ORDER TO GET THE RIGHT INTEGRAL:
 
		for (int m = 0; m <= N2; ++m) {
 
			num = {-bias, 2.*M_PI*m/dkapa};
      Il (l, num, tq, ilm);
      herecomplexpow (tq, num-2., goofy);
 
      //herecomplexpow (k0*chi0, -(num+bias), donald);

      cn [m] = cn0 [m] * ilm;// * donald / Nd;
      in2 [m][REAL] = cn [m].real ();
      in2 [m][IMAG] = cn [m].imag ();

      cn [m] *= goofy;
      in3 [m][REAL] = cn [m].real ();
      in3 [m][IMAG] = cn [m].imag ();
		}

		fftw_execute (pp2);
    fftw_execute (pp3);
 
    for (int n = 0; n < N; ++n) {
    	acca1 [n] = out2 [n] * inv_Nd;
      acca2 [n] = out3 [n] * inv_Nd;
    }
 
    res1 = 0.;
    for (int p = 0; p <= N-1; ++p) {

			dchi = chi [p] * expdNm1;

      if (chi [p] > chimin_photoz && chi[p] < chimax_photoz) {
				//phip = wf.Wgalaxy_noRSD (chi [p], l);
  			phip = gsl_spline_eval (spline_wg, chi[p], Wacc);
      }
			else {
				continue;
      }

      if (chi [p] *tq > chimin_photoz && chi[p]*tq < chimax_photoz) {
      	//phipt = wf.Wgalaxy_noRSD (chi [p]*tq, l);
  			phipt = gsl_spline_eval (spline_wg, chi[p]*tq, Wacc);
      }
      else {
				phipt = 0.;
      }
      if (chi [p] /tq > chimin_photoz && chi[p]/tq < chimax_photoz) {
      	//phipot = wf.Wgalaxy_noRSD (chi [p]/tq, l);
  			phipot = gsl_spline_eval (spline_wg, chi[p]/tq, Wacc);
      }
      else {
      	phipot = 0.;
      }

      //res1 += dchi * chi [p] * phip * (phipt * acca1 [p] + phipot * acca2 [p]);
      //res1 += dchi * pow (chi [p], 1.-bias) * 2. * pow (chi0, 2.*bias) * phip * (phipt * acca1 [p] + phipot * acca2 [p]);
      //res1 += dchi * pow (chi [p], 1.+bias) * pow (k0, bias) * phip * (phipt * acca1 [p] + phipot * acca2 [p]);                              // QUI
      res1 += dchi * pow (chi [p], 1.+bias) * phip * (phipt * acca1 [p] + phipot * acca2 [p]);

		}

    //res2 += dt * res1;
    res2 += res1 * inv_Nd;
  }

  return res2 / (2.*M_PI*M_PI);
}

Eigen::MatrixXd Cl_NonLimber (Eigen::MatrixXd larr, Eigen::MatrixXd Chi_in, Eigen::MatrixXd Wg_in, double chimin, double chimax, double bias, Eigen::MatrixXd karr_in, Eigen::MatrixXd Parr_in) {

//ATTENZIONE CHE EIGEN MATRIX LEGGE COLONNA PER COLONNA E NON RIGA PER RIGA

	double *Chi, *Wg;
	Chi = new double [Chi_in.size()];
	Wg  = new double [Chi_in.size()];
	for (int i = 0; i < Chi_in.size (); i++){
		Chi [i] = Chi_in (i);
		Wg [i] = Wg_in (i);
	}
	gsl_spline* spline_Wg = gsl_spline_alloc (gsl_interp_cspline, Chi_in.size ());
	gsl_spline_init (spline_Wg, Chi, Wg, Chi_in.size ());
	double *karr, *Parr;
	karr = new double [karr_in.size()];
	Parr  = new double [karr_in.size()];
	for (int i = 0; i < karr_in.size (); i++){
		karr [i] = karr_in (i);
		Parr [i] = Parr_in (i);
	}
	gsl_spline* spline_Plin = gsl_spline_alloc (gsl_interp_cspline, karr_in.size ());
	gsl_spline_init (spline_Plin, karr, Parr, karr_in.size ());

	for (int i = 0; i < larr.size (); i++)
		larr (i) = Cl_auto_gal (larr (i), chimin, chimax, bias, spline_Wg, spline_Plin);

  return larr;
}

Eigen::MatrixXd Cl_Limber (Eigen::MatrixXd larr, Eigen::MatrixXd Chi_in, Eigen::MatrixXd Wg_in, double chimin, double chimax, Eigen::MatrixXd karr_in, Eigen::MatrixXd Parr_in) {

//ATTENZIONE CHE EIGEN MATRIX LEGGE COLONNA PER COLONNA E NON RIGA PER RIGA

	double *Chi, *Wg;
	Chi = new double [Chi_in.size()];
	Wg  = new double [Chi_in.size()];
	for (int i = 0; i < Chi_in.size (); i++){
		Chi [i] = Chi_in (i);
		Wg [i] = Wg_in (i);
	}
	gsl_spline* spline_Wg = gsl_spline_alloc (gsl_interp_cspline, Chi_in.size ());
	gsl_spline_init (spline_Wg, Chi, Wg, Chi_in.size ());

	double *karr, *Parr;
	karr = new double [karr_in.size()];
	Parr  = new double [karr_in.size()];
	for (int i = 0; i < karr_in.size (); i++){
		karr [i] = karr_in (i);
		Parr [i] = Parr_in (i);
	}
	gsl_spline* spline_Plin = gsl_spline_alloc (gsl_interp_cspline, karr_in.size ());
	gsl_spline_init (spline_Plin, karr, Parr, karr_in.size ());

	for (int i = 0; i < larr.size (); i++){
		larr (i) = Cl_gal_limber_qags (larr (i), chimin, chimax, spline_Wg, spline_Plin);
	}

  return larr;
}

PYBIND11_PLUGIN(wrap) {
  pybind11::module m("wrap", "auto-compiled c++ extension");
  m.def("Cl_NonLimber", &Cl_NonLimber);
  m.def("Cl_Limber", &Cl_Limber);
  return m.ptr();
}

