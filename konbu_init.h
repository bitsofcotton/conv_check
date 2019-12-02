#include <iostream>
#include <sstream>
#include <float.h>
#include <limits>
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
#include <mpreal.h>
typedef mpfr::mpreal num_t;
//#define MANT_DIG (int)(ACC_GMP / 3.33)
#define MANT_DIG DBL_DIG
using mpfr::sqrt;
using mpfr::pow;
using mpfr::log;
using mpfr::isfinite;

#elif defined(ACC_QD_DDOUBLE)

#include <qd/dd_real.h>
#include <qd/dd_inline.h>
typedef dd_real num_t;
#define MANT_DIG DBL_DIG

#elif defined(ACC_QD_QDOUBLE)

#include <qd/qd_real.h>
#include <qd/qd_inline.h>
typedef qd_real num_t;
#define MANT_DIG DBL_DIG

#elif defined(ACC_NO_FLOAT)

#include <cmath>
#include <vector>
#include "ifloat.hh"
typedef SimpleFloat<uint64_t, DUInt<uint64_t, 64>, 64, short> num_t;
#define MANT_DIG DBL_DIG

#endif

#if defined(ACC_FLOAT) || defined(ACC_DOUBLE) || defined(ACC_LDOUBLE) || defined(ACC_NO_FLOAT)
#include <cmath>
using namespace std;
using std::numeric_limits;
template <typename T> T sgn(const T& x) { if(x == 0) return x; else return x > T(0) ? T(1) : - T(1); }
#elif defined(ACC_GMP)
const num_t& real(const num_t& x) { return x; }
const num_t  imag(const num_t& x) { return num_t(0); }
using std::numeric_limits;
#elif defined(ACC_QD_DDOUBLE) || defined(ACC_QD_QDOUBLE)
const num_t& real(const num_t& x) { return x; }
const num_t  imag(const num_t& x) { return num_t(0); }
using std::numeric_limits;
using std::max;
using std::min;
template <typename T> T sgn(const T& x) { if(x == 0) return x; else return x > T(0) ? T(1) : - T(1); }
#endif

#include "konbu.hh"

#if defined(WITHOUT_EIGEN)
typedef SimpleMatrix<num_t> Mat;
typedef SimpleVector<num_t> Vec;
#else
typedef Eigen::Matrix<num_t, Eigen::Dynamic, Eigen::Dynamic> Mat;
typedef Eigen::Matrix<num_t, Eigen::Dynamic, 1> Vec;
#endif

