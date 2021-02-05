<%
cfg['compiler_args'] = ['-std=c++11']
setup_pybind11(cfg)
%>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

using namespace pybind11::literals;

Eigen::MatrixXd Cl_auto_gal (Eigen::MatrixXd larr) {

  return larr * larr;
}

PYBIND11_PLUGIN(wrap, m) {
  pybind11::module m("wrap", "auto-compiled c++ extension");
  m.def("Cl_auto_gal", &Cl_auto_gal);
  return m.ptr();
}

