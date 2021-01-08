
#include "inc/3CInt/walltimer.h"
#include "inc/3CInt/angpow_func.h"
#include "inc/3CInt/3cint_chealgo.h"

/*
*/
namespace Angpow {
class ClassFunc1D_get_value : public ClassFunc1D {
  typedef ClassFunc1D parent;
public:
  virtual double get_value(double x) const {
    printf("debug : ClassFunc1D_get_value::get_value %g\n",x);
    return 0;
  }
  virtual double operator()(double x) const {return get_value(x);}
};
}

#include "Angpow_swig_py.icc"
