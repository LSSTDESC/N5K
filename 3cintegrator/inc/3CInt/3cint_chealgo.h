#ifndef THREECINT_CHEALGO_SEEN
#define THREECINT_CHEALGO_SEEN
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

#include <math.h> //cte numerique
#include <numeric> //inner_product
#include "angpow_fft.h"
#include "3cint_chefunc.h"


namespace Angpow {

class CheAlgo {
public:
  CheAlgo(std::vector<CheFunc*> farr): farr_(farr) {

    //compute the order of the final space
    // 1 + Sum_m 2^{N_m}
    nOrdProd_ = 1;
    for(size_t i=0;i<farr.size();i++){
      nOrdProd_ += farr[i]->orderFunc() - 1;
    }

    //    std::cout << "CheAlgo: nOrdProd = " << nOrdProd_ << std::endl;

    //initialize the FFT plan for the Chebyshev Inverse transform
    vecDCTInv_.resize(nOrdProd_,0.);
    planInv_ = new FFTPlanning(nOrdProd_,vecDCTInv_);
    
    //Initialize the Clenshow-Curtis single quadrature weight
    wCC_.resize(nOrdProd_,0);

    planCC_ = new FFTPlanning(nOrdProd_, wCC_);
    
    //Compute the Clenshow-Curtis single quadrature weight (done once)
    ClenshawCurtisWeightsFast();

//     std::cout << "CC weights: {";
//     for(size_t i=0;i<wCC_.size();i++) std::cout << wCC_[i] << " ";
//     std::cout << "}\n";

    
  }

  virtual ~CheAlgo() {
    if (planInv_){   delete planInv_;   planInv_   = 0;}
    if (planCC_){    delete planCC_;    planCC_    = 0;}
  }


  /*! ClenshawCurtisWeightsFast
    Determine the weights of the Clenshaw-Curtis quadrature using DCT-I algorithm.
    For the normalization one should take into account that
    FFTW transform 
    Y_k = 2 Sum''_j=0^(N-1) X_j Cos[Pi k j/(N-1)]
    while Clenshaw-Curtis weights are
    W_k = 4/(N-1) a_k  Sum''_{j=0, j even}^(N-1) 1/(1-j^2) Cos[Pi k j/(N-1)]
    
    (nb. Sum'' means that the first and last elements of the sum are devided by 2)
    
    \input plan the FFTPlanning object that handle the data and the type of transform
    \ouput w the weights defined for [-1, 1] quadrature
  */
  void ClenshawCurtisWeightsFast() {
    
    //dim of w is n
    fill(wCC_.begin(), wCC_.end(), (r_8)0.0);
    
    int n = wCC_.size();
    for(int k=0;k<n; k +=2){
      wCC_[k] = 1./(1.-(r_8)(k*k));
    }

    //FFTW planning associated to wCC_
    planCC_->Execute();
    
    r_8 norm = 2.0/(r_8)(n-1); //2 * FFTW DCT-1 
    std::transform(wCC_.begin(), wCC_.end(), wCC_.begin(), std::bind1st(std::multiplies<r_8>(),norm));
    wCC_[0] /= (r_8)2; wCC_[n-1] /= (r_8)2;


  }


  /*! InverseChebyshevCoeffFFT
    Compute the sampling using inverse Chebyshev transform coefficient computation
    use the vecDCTInv_ vector and the associated planInv_ FFT plan
   */
  void InverseChebyshevCoeffFFT(){

    vecDCTInv_.front() *= 2.0;
    vecDCTInv_.back()  *= 2.0;
    
    planInv_->Execute();
    
    std::transform(vecDCTInv_.begin(), vecDCTInv_.end(), 
		   vecDCTInv_.begin(), std::bind1st(std::multiplies<r_8>(),0.5)); 
  }

  /*! InverseChebyshevTransform
     perform for each functions the inverse Chebyshev transform in the final dimension space
     update their local sampling vectors
   */
  void InverseChebyshevTransform(){
    for(size_t i=0;i<farr_.size();i++){
      std::fill(vecDCTInv_.begin(), vecDCTInv_.end(), 0.0); //may be not necessary
      std::vector<r_8> vec = farr_[i]->vec();
      std::copy(vec.begin(),vec.end(),vecDCTInv_.begin());
      InverseChebyshevCoeffFFT();
      farr_[i]->updateVec(vecDCTInv_);
    } 
  }


  /*! ComputeIntegralUnscaled
    Compute the non-scale integral computed as 
    h(x) = Prod_i f_i(x)
    integral = Sum_p wCC(p) h(x_p)
    with the h(x_p) from the Chebyshev Transform of each f_(x) 
    and the Expension in the final Chebyshev space
    and the Inverse Chebyshev Transform in the real space.

    Warning: to be properly used mind to multiply by the integral interval length.
   */
  r_8 ComputeIntegralUnscaled() {
    
    // element-by-element product of all the sampling vectors
    std::vector<r_8> prod(nOrdProd_,1.0); //init at 1
    for(size_t i=0;i<farr_.size();i++){
      std::vector<r_8> vec = farr_[i]->vecFinal();
      //      if(vec.size() != prod.size()) printf("HORREUR ComputeIntegralUnscaled !!!! ");
      for(int j=0;j<nOrdProd_;j++){
	prod[j] *= vec[j];
      }
    }
    
    //scalar product Clenshaw-Curtis weights vector and the aboave "prod" vector
    r_8 integral = inner_product(prod.begin(),prod.end(),wCC_.begin(),0.);
    integral *= 0.5; //the 2 division comme from CC quadrature computation throw FFT

    return integral;  
  }

private:

  std::vector<CheFunc*> farr_; //!< array of functions to be processed
  
  int nOrdProd_; //!< size of the vector vecDCTInv_ and wCC_

  std::vector<r_8> vecDCTInv_; //!< vector used by FFTW for the inversion of the product
  FFTPlanning* planInv_;  //!< FFTW plan associated to vecDCTInv_

  std::vector<r_8> wCC_; //!< vector of the Clenshaw-Curtis weights 
  FFTPlanning*  planCC_; //! FFTW plan associated to wCC_

};

}//end namespace
#endif //THREECINT_CHEALGO_SEEN
