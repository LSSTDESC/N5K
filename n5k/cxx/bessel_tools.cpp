#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp>

pybind11::array_t<double> bessel_zeros(int L, int N) {

  // Allocate the qln table and copy over the zeros
  size_t size = (L+1)*N;
  double *qln = new double[size];

  #pragma omp parallel for schedule(guided)
  for(int l=0; l <= L; l++) {
      std::vector<long double> roots;
      boost::math::cyl_bessel_j_zero((double) (l+0.5), 1, N, std::back_inserter(roots));
      for(int p=0; p <N; p++) {
          qln[N*l + p] = roots[p];
      }
  }

  // Create a Python object that will free the allocated
  // memory when destroyed:
  pybind11::capsule free_when_done(qln, [](void *f) {
      double *qln = reinterpret_cast<double *>(f);
      std::cerr << "Element [0] = " << qln[0] << "\n";
      std::cerr << "freeing memory @ " << f << "\n";
      delete[] qln;
  });

  return pybind11::array_t<double>(
      {L+1, N}, // shape
      {N*8, 8}, // C-style contiguous strides for double
      qln, // the data pointer
      free_when_done); // numpy array references this parent
}

PYBIND11_MODULE(bessel_tools, m) {
  m.doc() = "Module for Bessel stuff";
  m.def("bessel_zeros", &bessel_zeros, "compute Bessel zeros");
}
