#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp>

namespace py = pybind11;

// N is the number of zeros we are looking for
// L is an array of orders l
py::array_t<double> bessel_zeros(int N, py::array_t<uint64_t> L) {

  py::buffer_info info = L.request();
  if (info.ndim != 1)
      throw std::runtime_error("Number of dimensions must be one");
  // Number of entries in the L array
  int Nl = info.shape[0];
  // Accessing the array values
  uint64_t *Lptr = static_cast<uint64_t *>(info.ptr);

  // Allocate the qln table and copy over the zeros
  size_t size = Nl*N;
  double *qln = new double[size];

  #pragma omp parallel for schedule(dynamic)
  for(int l=0; l < Nl; l++) {
      std::vector<long double> roots;
      boost::math::cyl_bessel_j_zero((double) (Lptr[l]+0.5), 1, N, std::back_inserter(roots));
      for(int p=0; p <N; p++) {
          qln[N*l + p] = roots[p];
      }
  }
  // Create a Python object that will free the allocated
  // memory when destroyed:
  py::capsule free_when_done(qln, [](void *f) {
      double *qln = reinterpret_cast<double *>(f);
      delete[] qln;
  });

  return py::array_t<double>(
      {Nl, N}, // shape
      {N*8, 8}, // C-style contiguous strides for double
      qln, // the data pointer
      free_when_done); // numpy array references this parent
}

PYBIND11_MODULE(bessel_tools, m) {
  m.doc() = "Module for Bessel stuff";
  m.def("bessel_zeros", &bessel_zeros, "compute Bessel zeros");
}
