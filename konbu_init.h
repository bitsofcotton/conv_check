#include <iostream>
#include <sstream>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>

#if defined(ACC_FLOAT)

typedef float num_t;
#define MANT_DIG FLT_DIG

#elif defined(ACC_DOUBLE)

typedef double num_t;
#define MANT_DIG DBL_DIG

#elif defined(ACC_LDOUBLE)

typedef long double num_t;
#define MANT_DIG LDBL_DIG

#elif defined(ACC_GMP)

#define REAL_ENABLE_CONVERSION_OPERATORS
#include <mpfr.h>
#include <real.hpp>
typedef mpfr::real<BITS> num_t;
//#define MANT_DIG (int)(BITS / 3.33)
#define MANT_DIG DBL_DIG
using namespace mpfr;

#endif

using namespace std;

#if defined(ACC_FLOAT) || defined(ACC_DOUBLE) || defined(ACC_LDOUBLE)
#include <cmath>
#define LIBNAME std
#elif defined (ACC_GMP)
#define LIBNAME mpfr
#endif

namespace Eigen {
  namespace internal {
    num_t max(const num_t& x, const num_t& y) { return x < y ? y : x; }
    num_t min(const num_t& x, const num_t& y) { return x < y ? x : y; }
    num_t conj(const num_t& x) { return x; }
    num_t real(const num_t& x) { return x; }
    num_t imag(const num_t&) { return num_t(0); }
    num_t abs(const num_t& x) { return x < num_t(0) ? - x : x; }
    num_t abs2(const num_t& x) { return x * x; }
    num_t sqrt(const num_t& x) { return LIBNAME::sqrt(x); }
    num_t pow(const num_t& x, num_t y)  { return LIBNAME::pow(x, y); }
    num_t log(const num_t& x) { return LIBNAME::log(x); }
    num_t tan(const num_t& x) { return LIBNAME::tan(x); }
    num_t atan2(const num_t& y, const num_t& x) { return LIBNAME::atan2(y, x); }
    num_t floor(const num_t& x) { return LIBNAME::floor(x); }
    num_t ceil(const num_t& x) { return LIBNAME::ceil(x); }
    bool  isfinite(const num_t& x) { return LIBNAME::isfinite(x); }
  }
}

num_t      m_epsilon(0);

#include <Eigen/Core>
#include <Eigen/Dense>

#if defined(ACC_GMP)
namespace Eigen {
  template<> struct NumTraits<num_t>
  {
    typedef num_t Real;
    typedef num_t NonInteger;
    typedef num_t Nested;
    typedef int   Index;
    
    inline static Real epsilon() {
      return m_epsilon;
    };
    inline static Real dummy_precision() {
      return m_epsilon;
    }
    inline static num_t highest() {
      return 1e100;
    }
    inline static num_t lowest() {
      return -1e100;
    }
    
    enum {
      IsComplex = 0,
      IsInteger = 0,
      IsSigned = 1,
      RequireInitialization = 1,
      ReadCost = 1,
      AddCost = 1,
      MulCost = 1
    };
  };
}
#endif

using namespace Eigen;
typedef Eigen::Matrix<num_t, Eigen::Dynamic, Eigen::Dynamic> Mat;
typedef Eigen::Matrix<num_t, Eigen::Dynamic, 1> Vec;

