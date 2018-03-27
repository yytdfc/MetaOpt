#include <pybind11/stl.h>
#include "model/gp.h"
// #include <pybind11/eigen.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

namespace MetaOpt {
template <typename Real>
class PyGp : public Gp<Real>
{
 public:
  PyGp(const int i_kernel = 3, const int i_order = 0, const bool i_norm = true)
      : Gp<Real>(i_kernel, i_order, i_norm){};
  void fit_vec(const std::vector<std::vector<Real>>& x,
               const std::vector<std::vector<Real>>& y) {
    Samples<Real> samples(x.size());
    for (int i = 0; i != x.size(); ++i) {
      samples[i] = Sample<Real>(x[i], y[i]);
    }
    this->fit(move(samples));
  };
  py::tuple predict_vec(const std::vector<std::vector<Real>>& x,
                        const bool return_std = false) {
    std::vector<std::vector<Real>> y(x.size());
    std::vector<std::vector<Real>> st(0);
    std::vector<Real>              t(0);
    if (return_std) {
      st = std::vector<std::vector<Real>>(x.size());
    }
    for (int i = 0; i != x.size(); ++i) {
      if (return_std) {
        t = std::vector<Real>(this->n_y_);
      }
      Sample<Real> s(x[i], std::vector<Real>(this->n_y_));
      this->evaluate(s, t);
      y[i] = move(s.obj());
      if (return_std) {
        st[i] = move(t);
      }
    }
    if (return_std) {
      return py::make_tuple(y, st);
    } else {
      return py::make_tuple(y);
    }
  };
};

PYBIND11_MODULE(gp, m) {
  py::class_<PyGp<double>>(m, "gp")
      .def(py::init<const int, const int, const bool>(),
           "A GaussProcess constructor", py::arg("kernel") = 3,
           py::arg("order") = 1, py::arg("norm") = true)
      .def("fit", &PyGp<double>::fit_vec, "fitting f(x) = y(x)", py::arg("x"),
           py::arg("y"))
      .def("predict", &PyGp<double>::predict_vec, "fitting f(x) = y(x)",
           py::arg("y"), py::arg("return_std") = false);
}

}  // namespace MetaOpt