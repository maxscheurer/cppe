#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <armadillo>

namespace py = pybind11;

// type caster: Matrix <-> NumPy-array
namespace pybind11 { namespace detail {
  template <typename T>
  struct type_caster<arma::Mat<T>> {
    public:

      PYBIND11_TYPE_CASTER(arma::Mat<T>, _("arma::Mat<T>"));

      // Conversion part 1 (Python -> C++)
      bool load(py::handle src, bool convert)
      {
        if (!convert && !py::array_t<T>::check_(src))
          return false;

        auto buf = py::array_t<T, py::array::f_style | py::array::forcecast>::ensure(src);
        if (!buf)
          return false;

        auto dims = buf.ndim();
        if (dims < 1 || dims > 2)
          return false;

        std::vector<size_t> shape(buf.ndim());

        for ( int i=0 ; i<buf.ndim() ; i++ )
          shape[i] = buf.shape()[i];
        
        if (dims == 1) {
            value = arma::Col<T>(buf.data(), shape[0]);
        } else if (dims == 2) {
            value = arma::Mat<T>(buf.data(), shape[0], shape[1]);
        }

        return true;
      }

      //Conversion part 2 (C++ -> Python)
      static py::handle cast(const arma::Mat<T>& src,
        py::return_value_policy policy, py::handle parent) {
        std::vector<size_t> shape(2);
        shape[0] = src.n_rows;
        shape[1] = src.n_cols;
        std::vector<size_t> strides(2);
        std::vector<size_t> ret(2);
        strides[0] = 1;
        strides[1] = shape[0];
        // for later generalization
        for ( size_t i = 0 ; i < 2 ; ++i )
          ret[i] = strides[i] * sizeof(T) ;
        py::array a(std::move(shape), std::move(ret), src.memptr());

        return a.release();
      }
      
  };
  template <typename T>
  struct type_caster<arma::Col<T>> {
    public:
        PYBIND11_TYPE_CASTER(arma::Col<T>, _("arma::Col<T>"));

        // Conversion part 1 (Python -> C++)
        bool load(py::handle src, bool convert) {
          if (!convert && !py::array_t<T>::check_(src))
            return false;

          auto buf = py::array_t<T, py::array::f_style | py::array::forcecast>::ensure(src);
          if (!buf)
            return false;

          auto dims = buf.ndim();
          if (dims < 1 || dims > 2)
            return false;

          std::vector<size_t> shape(buf.ndim());

          for ( int i=0 ; i<buf.ndim() ; i++ )
            shape[i] = buf.shape()[i];
          
          value = arma::Col<T>(buf.data(), shape[0]);

          return true;
        }
        
        static py::handle cast(const arma::Col<T>& src,
          py::return_value_policy policy, py::handle parent) {
          std::vector<size_t> shape(1);
          shape[0] = src.n_rows;
          std::vector<size_t> strides(1);
          std::vector<size_t> ret(1);
          strides[0] = 1;
          // strides[1] = shape[0];
          // for later generalization
          ret[0] = strides[0] * sizeof(T) ;
          py::array a(std::move(shape), std::move(ret), src.memptr());

          return a.release();
        }
  };
}}