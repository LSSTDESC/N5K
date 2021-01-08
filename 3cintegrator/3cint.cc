#include <string.h>
#include <cstdio>
#include <stdlib.h>
#include <iostream>

#include "3CInt/angpow_exceptions.h"
#include "3CInt/angpow_numbers.h"
#include "3CInt/angpow_func.h"
#include "3CInt/3cint_chefunc.h"
#include "3CInt/3cint_chealgo.h"
#include "3CInt/walltimer.h"


namespace Angpow {

struct PARAM {
  int ell;
  r_8 kMin;
  r_8 kMax;
  r_8 R1;
  r_8 R2;
  int chebyshev_order_1;
  int chebyshev_order_2;
  int n_sub_intervals;
} para ;

//smooth common function
class FuncType0: public ClassFunc1D {
public:
  FuncType0() {}
  inline virtual r_8 operator()(r_8 x) const {
    return x * (x-1.)*(x-1.);
  }
  virtual ~FuncType0() {}
private:
  r_8 p_;
  r_8 scale_; 
};

//high oscillatory function
class FuncType1: public ClassFunc1D {

public:
  FuncType1(int ell, r_8 R): ell_(ell), R_(R) {}
  inline virtual r_8 operator()(r_8 x) const {
    return cos(x*R_ -ell_*M_PI*0.5 - M_PI*0.25);
  }
  virtual ~FuncType1() {}
private:
  int ell_;
  r_8 R_; 

};//ProdJBess


  r_8 truth_test0(int l, r_8 R1, r_8 R2, r_8 a, r_8 b) {
  
  r_8 res = 0;
 
  if(R1 != R2) {
    res = ((-3*(-2 + pow(a,2.)*pow(R1 - R2,2.))*cos(a*(R1 - R2)))/pow(R1 - R2,4.) - 
     cos(a*(R1 - R2))/pow(R1 - R2,2.) + (4*a*cos(a*(R1 - R2)))/pow(R1 - R2,2.) + 
     (3*(-2 + pow(b,2.)*pow(R1 - R2,2.))*cos(b*(R1 - R2)))/pow(R1 - R2,4.) + 
     cos(b*(R1 - R2))/pow(R1 - R2,2.) - (4*b*cos(b*(R1 - R2)))/pow(R1 - R2,2.) - 
     (a*(-6 + pow(a,2.)*pow(R1 - R2,2.))*sin(a*(R1 - R2)))/pow(R1 - R2,3.) + 
     (2*(-2 + pow(a,2.)*pow(R1 - R2,2.))*sin(a*(R1 - R2)))/pow(R1 - R2,3.) - 
     (a*sin(a*(R1 - R2)))/(R1 - R2) + 
     (b*(-6 + pow(b,2.)*pow(R1 - R2,2.))*sin(b*(R1 - R2)))/pow(R1 - R2,3.) - 
     (2*(-2 + pow(b,2.)*pow(R1 - R2,2.))*sin(b*(R1 - R2)))/pow(R1 - R2,3.) + 
     (b*sin(b*(R1 - R2)))/(R1 - R2) + 
     (a*(R1 + R2)*cos(l*M_PI - a*(R1 + R2)) + sin(l*M_PI - a*(R1 + R2)))/pow(R1 + R2,2.) - 
     (2*((-2 + pow(a,2.)*pow(R1 + R2,2.))*cos(l*M_PI - a*(R1 + R2)) + 
          2*a*(R1 + R2)*sin(l*M_PI - a*(R1 + R2))))/pow(R1 + R2,3.) + 
     (a*(R1 + R2)*(-6 + pow(a,2.)*pow(R1 + R2,2.))*cos(l*M_PI - a*(R1 + R2)) + 
        3*(-2 + pow(a,2.)*pow(R1 + R2,2.))*sin(l*M_PI - a*(R1 + R2)))/pow(R1 + R2,4.) - 
     (b*(R1 + R2)*cos(l*M_PI - b*(R1 + R2)) + sin(l*M_PI - b*(R1 + R2)))/pow(R1 + R2,2.) + 
     (2*((-2 + pow(b,2.)*pow(R1 + R2,2.))*cos(l*M_PI - b*(R1 + R2)) + 
          2*b*(R1 + R2)*sin(l*M_PI - b*(R1 + R2))))/pow(R1 + R2,3.) - 
     (b*(R1 + R2)*(-6 + pow(b,2.)*pow(R1 + R2,2.))*cos(l*M_PI - b*(R1 + R2)) + 
      3*(-2 + pow(b,2.)*pow(R1 + R2,2.))*sin(l*M_PI - b*(R1 + R2)))/pow(R1 + R2,4.))/2.;

  } else {
    res = (-2*pow(a,2.)*(6 - 8*a + 3*pow(a,2.))*pow(R1,4.) + 
     2*pow(b,2.)*(6 - 8*b + 3*pow(b,2.))*pow(R1,4.) + 
     6*R1*(2 - 4*pow(a,2.)*pow(R1,2.) + 2*pow(a,3.)*pow(R1,2.) + 
        a*(-3 + 2*pow(R1,2.)))*cos(l*M_PI - 2*a*R1) - 
     6*R1*(2 - 4*pow(b,2.)*pow(R1,2.) + 2*pow(b,3.)*pow(R1,2.) + 
        b*(-3 + 2*pow(R1,2.)))*cos(l*M_PI - 2*b*R1) + 
     3*(-3 + (2 - 8*a + 6*pow(a,2.))*pow(R1,2.))*sin(l*M_PI - 2*a*R1) - 
	   3*(-3 + (2 - 8*b + 6*pow(b,2.))*pow(R1,2.))*sin(l*M_PI - 2*b*R1))/(48.*pow(R1,4.)); 
  }
  return res;
}


  /*
    Example of the computation of the integral
    \int_{kmin}^{kmax} f0(k) f1(kR1,ell) f2(kR2,ell) dx
    with 
    f0 a common smooth function
    f1, f2 two higly oscillatory functions
   */


void test0() {
  r_8 kMin = para.kMin;
  r_8 kMax = para.kMax;
  //  int Lmax = para.Lmax; //ell<Lmax
  int nSubInterv = para.n_sub_intervals;
  std::vector<r_8> R(2);
  R[0] = para.R1;
  R[1] = para.R2;
  int ell = para.ell;

  //k-integral bounds
  std::vector<r_8> klp(nSubInterv+1);
  r_8 dK = kMax-kMin;
  for(int i=0; i<= nSubInterv; i++){
    klp[i] = kMin + dK * i/((r_8)nSubInterv);
  }
  printf("ell=%d, Nintervales=%d\n", ell,nSubInterv);



  //Function to be integrated
  FuncType1* f1 = new FuncType1(ell,R[0]);
  FuncType1* f2 = new FuncType1(ell,R[1]);
  FuncType0* f0 = new FuncType0();
  
  // Chebyshev machinery 
  int iOrd0 = 2;
  int iOrd1 = para.chebyshev_order_1;
  int iOrd2 = para.chebyshev_order_2;
  
  std::vector<CheFunc*> farr;
  farr.push_back(new CheFunc(f1, iOrd1));
  farr.push_back(new CheFunc(f2, iOrd2));
  farr.push_back(new CheFunc(f0, iOrd0));
  
  //Initialisation of the Clenshow-Curtis quadrature
  CheAlgo cheAlgo(farr);

  //Integration
  r_8 integral = 0.;

  for(int p = 1; p<= nSubInterv; p++){ //init at p=1

    //get the bounds
    r_8 lowBound = klp[p-1];
    r_8 uppBound = klp[p];
//     std::cout << "current interval: [" << lowBound << ", " << uppBound << "]" 
// 	      << std::endl;
    
    if(lowBound > uppBound)
      throw AngpowError("KIntegrator::Compute uppBound < lowBound Fatal");
    
    //Loop on each function to compute their  Foward Chebyshev coefficents
    for(size_t i=0;i<farr.size();i++){
      farr[i]->ChebyshevTransform(lowBound, uppBound);
    }

    //Compute the sampling of all the functions in the final space dimension
    cheAlgo.InverseChebyshevTransform();

    //Compute the integral thanks to CC quadrature and the function sampling 
    integral += (uppBound - lowBound) * cheAlgo.ComputeIntegralUnscaled();
    
  }//p-loop 

  std::cout << "Approx. Integ = " << integral << std::endl;

  //truth
  r_8 true_int = truth_test0(ell,R[0],R[1],kMin,kMax);
  std::cout << "True Integ = " << true_int 
	    << " diff= " << true_int - integral
	    << std::endl;



  //-------
  // clean
  //-------
  // func
  if(f0) delete f0;
  if(f1) delete f1;
  if(f2) delete f2;
  // Che-stuff
  for(size_t i=0;i<farr.size();i++) delete farr[i];
  farr.clear();

    

  std::cout << "End test0......" << std::endl;

}

}//namespace

//----------------------------------------------
//               Main 
//----------------------------------------------

int main(int narg, char *arg[]) {
  
  using namespace Angpow;

  //  unsigned int maxmemsize = getMemorySize()/1e6;
  //std::cout << "Max Memory size: " <<  maxmemsize << " MBytes" <<  std::endl;

  //The cosmological distance tool
   
  int test=0;
  int ell= 20;
  r_8 R1 = 2000.; //Mpc z=1.0
  r_8 R2 = 2200.; //Mpc z=1.1
  r_8 kMin = 0.;
  r_8 kMax = 1.0; //Mpc^(-1)

  int chebyshev_order_1 = 8;
  int chebyshev_order_2 =  chebyshev_order_1;
  int n_sub_intervals = 5;

  int ka=1;
  while (ka<narg) {
    if (strcmp(arg[ka],"-h")==0) {
      return 0;
    }
    else if (strcmp(arg[ka],"-l")==0) {
      ell=atoi(arg[ka+1]);
      ka+=2;
    }
    else if (strcmp(arg[ka],"-kmin")==0) {
      kMin=atof(arg[ka+1]);
      ka+=2;
    }
    else if (strcmp(arg[ka],"-kmax")==0) {
      kMax=atof(arg[ka+1]);
      ka+=2;
    }
    else if (strcmp(arg[ka],"-r1")==0) {
      R1 =  atof(arg[ka+1]);
      ka+=2;      
    }    
    else if (strcmp(arg[ka],"-r2")==0) {
      R2 =  atof(arg[ka+1]);
      ka+=2;      
    }    
    else if (strcmp(arg[ka],"-che1")==0) {
      chebyshev_order_1 =  atof(arg[ka+1]);
      ka+=2;      
    }    
    else if (strcmp(arg[ka],"-che2")==0) {
      chebyshev_order_2 =  atof(arg[ka+1]);
      ka+=2;      
    }    
    else if (strcmp(arg[ka],"-nroot")==0) {
      n_sub_intervals= atoi(arg[ka+1]);
      ka+=2;      
    }    
    else ka++;
  }//eo while


  para.ell = ell;
  para.kMin = kMin;
  para.kMax = kMax;
  para.R1 = R1;
  para.R2 = R2;

  para.chebyshev_order_1  = chebyshev_order_1;
  para.chebyshev_order_2  = chebyshev_order_2;
  para.n_sub_intervals = n_sub_intervals;
    

  std::cout << "Configuration parameters are set to: " << std::endl;

  int rc=0;
  try {
    
    switch(test) {
    case 0:
      test0();
      break;
    default:
      throw AngpowError("EROOR Test type unkwown");
    }//end of switch


    std::cout << "---/ Fin bloc try ---- " << std::endl;
  }
    
  catch (AngpowError & e) {
    std::cerr << " besssht.cc: Catched Exception (AngpowError)" << (std::string)typeid(e).name() 
	 << " - Msg= " << e.what() << std::endl;
    rc = 99;
  }
  catch (std::exception & e) {
    std::cerr << " 3cint.cc: Catched std::xception "  
	 << " - what()= " << e.what() << std::endl;
    rc = 98;
  }
  catch (...) {
    std::cerr << " 3cint.cc: some other exception (...) was caught ! " << std::endl;
    rc = 97;
  }
  std::cout << " ---- Programme 3cint.cc -  FIN  (Rc=" << rc << ") --- " << std::endl;
  return rc;
}//main
