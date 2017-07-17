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
typedef mpfr::real<ACC_GMP> num_t;
//#define MANT_DIG (int)(ACC_GMP / 3.33)
#define MANT_DIG DBL_DIG

#endif

#if defined(ACC_FLOAT) || defined(ACC_DOUBLE) || defined(ACC_LDOUBLE)
#include <cmath>
using namespace std;
#elif defined (ACC_GMP)
const num_t& real(const num_t& x) { return x; }
const num_t  imag(const num_t& x) { return num_t(0); }
#endif

#include "konbu.hh"

#if defined(WITHOUT_EIGEN)
typedef SimpleMatrix<num_t> Mat;
typedef SimpleVector<num_t> Vec;
#else
typedef Eigen::Matrix<num_t, Eigen::Dynamic, Eigen::Dynamic> Mat;
typedef Eigen::Matrix<num_t, Eigen::Dynamic, 1> Vec;
#endif

