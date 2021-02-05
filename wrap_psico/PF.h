#ifndef PRELIMINARYFUNCTIONS_H
#define PRELIMINARYFUNCTIONS_H

#include <iostream>
#include <complex>

void GammaC (std::complex<double> z, std::complex<double>& result);

void GammaRatioC (std::complex <double> z1, std::complex <double> z2, std::complex <double>& result);

void Hyp2F1basic (std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> z, std::complex<double>& s);

void Il (int lin, std::complex<double> nu, double t, std::complex <double>& result);

#endif
