%module(directors="1") Angpow

class wallTimer {
public:
  wallTimer();
  void start();
  void stop();
};

//////////////////////////////
/// angpow_func.h : //////////
//////////////////////////////

namespace Angpow {

typedef double r_8;
  
class ClassFunc1D {
public:  
  ClassFunc1D();
  virtual ~ClassFunc1D();
public:
  virtual double operator()(double x) const = 0;
};

%feature("director") ClassFunc1D_get_value;
class ClassFunc1D_get_value : public ClassFunc1D {
public:  
  ClassFunc1D_get_value();
  virtual ~ClassFunc1D_get_value();
public:
  virtual double get_value(double x) const;
  virtual double operator()(double x);
};
 
//////////////////////////////
/// 3cint_chefunc.h : ////////
//////////////////////////////

class CheFunc {
public:
  CheFunc(ClassFunc1D* func, int ordFunc);
  virtual ~CheFunc();
public:  
  void ChebyshevTransform(r_8 a, r_8 b);
};

//////////////////////////////
/// 3cint_chealgo.h : ////////
//////////////////////////////

class CheAlgo {
public:
  CheAlgo(std::vector<CheFunc*> farr);
  virtual ~CheAlgo();
public:
  void InverseChebyshevTransform();
  r_8 ComputeIntegralUnscaled();
};

//////////////////////////////
//////////////////////////////
//////////////////////////////

} //end namespace Andpow

%include std_vector.i
%template(std_vector_CheFunc)  std::vector<Angpow::CheFunc*>;
%template(std_vector_double)   std::vector<double>;
