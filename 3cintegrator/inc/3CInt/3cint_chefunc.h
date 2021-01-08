#ifndef THREECINT_CHEFUN_SEEN
#define THREECINT_CHEFUN_SEEN
/*
 *  This file is part of Angpow.
 *
 *  Angpow is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Angpow is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Angpow; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  Angpow is being developed at the Linear Accelerateur Laboratory (LAL)
 *  91898 ORSAY CEDEX - FRANCE
 *  main author: J.E Campagne
 *    co-authors: J. Neveu, S. Plaszczynski
 */

#include <math.h> //pow
#include <vector>
#include <algorithm>
#include "angpow_func.h"
#include "angpow_fft.h"

namespace Angpow {

class CheFunc {
 public:
  CheFunc(ClassFunc1D* func, int ordFunc): func_(func), ordFunc_(ordFunc) {
    nOrdFunc_ = pow(2.,(r_8)ordFunc_)+1; //Chebyshev poly degree Forward
    vecDCTFunc_.resize(nOrdFunc_,0.);    //output of DCT forward
    planFunc_ = new FFTPlanning(nOrdFunc_,vecDCTFunc_); //FFTW DCT
  }

  virtual ~CheFunc() {
    if (planFunc_){ delete planFunc_; planFunc_ = 0;}
  }


  //getter
  int orderFunc() {return nOrdFunc_;}
  std::vector<r_8>& vec() {return vecDCTFunc_;}
  std::vector<r_8>& vecFinal() {return vecDCTInvFunc_;}

  //setter
  void updateVec(std::vector<r_8>& vec) {
    vecDCTInvFunc_.resize(vec.size());
    std::copy(vec.begin(),vec.end(),vecDCTInvFunc_.begin());
  }

  // Single value evaluation
  virtual r_8 operator()(r_8 x) const { return (*func_)(x); }

  // Vectorized evaluations
  virtual void getValues(const std::vector<double>& vin, 
			 std::vector<double>& vout) const {
    func_->getValues(vin,vout);
  }

  /*! ChebyshevSampling
    Perform function sampling on the Chebyshev points defined in the range [a,b]
    \input a lower bound of the range
    \input b upper bound of the range
    result of the sampling is in vecDCTFunc
  */
  inline void ChebyshevSampling(r_8 a, r_8 b) {
    int n =  nOrdFunc_  -1;
    r_8 bma = 0.5*(b-a); r_8 bpa = 0.5*(b+a);
    r_8 cte = M_PI/((r_8)n);
    
    std::vector<r_8> kin(n+1);
    for(int k=0;k<=n;k++){
      kin[k] = cos((n-k)*cte)*bma+bpa;
    }
    //call the vectorized f(k_i) = f_i
    func_->getValues(kin,vecDCTFunc_);
    
    //reveser the order to keep the same algorithm
    std::reverse(vecDCTFunc_.begin(), vecDCTFunc_.end());
  }


  /*! ChebyshevCoeffFFT
    Compute the Chebyshev transform using FFT
    use the FFT plan of the Function
    result of the FFT is in vecDCTFunc
  */
   void ChebyshevCoeffFFT(){

     int n =  nOrdFunc_ -1;

      planFunc_->Execute();
    //
    //  Chebyshev Coeff. the noramization of FFTW is propto  = 1/n
    //
     r_8 norm = 1./((r_8)n); 
     std::transform(vecDCTFunc_.begin(), vecDCTFunc_.end(), 
		    vecDCTFunc_.begin(), std::bind1st(std::multiplies<r_8>(),norm));
    
     vecDCTFunc_[0] *= 0.5;
     vecDCTFunc_[n] *= 0.5;
  }


  /*! ChebyshevTransform
    Foward Chebyshev coefficent computation: 
    (1) k-sampling 
    (2) DCT-I throw the FFT
    \input a lower bound of the range
    \input b upper bound of the range
    \output the result is in vecDCTFunc
  */
  void ChebyshevTransform(r_8 a, r_8 b) {
    ChebyshevSampling(a,b);

//     std::cout << "k-Sampling : {";
//     for(size_t i=0;i<vecDCTFunc_.size(); i++) std::cout << vecDCTFunc_[i] << " ";
//     std::cout << "}\n";

    ChebyshevCoeffFFT();
  }


 private:
  ClassFunc1D* func_; //no owner
  int ordFunc_; //!< parameter that define the poly degree
  int nOrdFunc_;  //!<  order of Chebychev polynomial
  std::vector<r_8> vecDCTFunc_; //!< vector used by FFTW for fun
  FFTPlanning*  planFunc_; //!< FFTW plan using  vecDCTFunc_
  std::vector<r_8> vecDCTInvFunc_;
};

}//end namespace

#endif //THREECINT_CHEFUN_SEEN
