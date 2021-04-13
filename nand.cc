#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <assert.h>

#if !defined(_FLOAT_BITS_)
  #include <complex>
  #include <cmath>
  using namespace std;
  typedef long double num_t;
#else
  #include "ifloat.hh"
  template <typename T> using complex = Complex<T>;
# if _FLOAT_BITS_ == 8
    typedef uint8_t myuint;
    typedef int8_t  myint;
    typedef SimpleFloat<myuint, uint16_t, 8, int64_t> num_t;
    #define mybits 8
# elif _FLOAT_BITS_ == 16
    typedef uint16_t myuint;
    typedef int16_t  myint;
    typedef SimpleFloat<myuint, uint32_t, 16, int64_t> num_t;
    #define mybits 16
# elif _FLOAT_BITS_ == 32
    typedef uint32_t myuint;
    typedef int32_t  myint;
    typedef SimpleFloat<myuint, uint64_t, 32, int64_t> num_t;
    #define mybits 32
# elif _FLOAT_BITS_ == 64
    typedef uint64_t myuint;
    typedef int64_t  myint;
    typedef SimpleFloat<myuint, DUInt<myuint, 64>, 64, int64_t> num_t;
    #define mybits 64
# elif _FLOAT_BITS_ == 128
    typedef DUInt<uint64_t, 64> uint128_t;
    typedef Signed<uint128_t, 128> int128_t;
    typedef uint128_t myuint;
    typedef int128_t  myint;
    typedef SimpleFloat<myuint, DUInt<myuint, 128>, 128, int64_t> num_t;
    #define mybits 128
# elif _FLOAT_BITS_ == 256
    typedef DUInt<uint64_t, 64> uint128_t;
    typedef DUInt<uint128_t, 128> uint256_t;
    typedef Signed<uint256_t, 128> int256_t;
    typedef uint256_t myuint;
    typedef int256_t  myint;
    typedef SimpleFloat<myuint, DUInt<myuint, 256>, 256, int64_t> num_t;
    #define mybits 256
# else
#   error cannot handle float
# endif
#endif

#include "simplelin.hh"
#include "konbu.hh"

int main(int argc, char* argv[])
{
  std::cout.precision(30);
  std::cerr.precision(30);
  SimpleMatrix<num_t> A(4, 3);
  SimpleVector<num_t> b(A.rows());
  SimpleVector<num_t> one(A.rows());
  const auto computer(- 1);
  num_t MM(0);
  for(int j = 0; j < one.size(); j ++) {
    A(j, 0) = num_t((j & 1) ? 1 : 2);
    A(j, 1) = num_t(((j >> 1) & 1) ? 1 : 2);
    //A(j, 2) = num_t(((j & 1) && ((j >> 1) & 1)) ? 2 : 1);
    A(j, 2) = num_t(1);
    b[j]    = num_t(((j & 1) && ((j >> 1) & 1)) ? 2 : 1);
    one[j]  = num_t(1);
    const auto r(pow(A(j, 0) * A(j, 1) * A(j, 2) * b[j], num_t(computer)));
    for(int i = 0; i < A.cols(); i ++)
      A(j, i) = num_t(1) / A(j, i);
    b[j]      = num_t(1) / b[j];
    A.row(j) *= r;
    b[j]     *= r;
    MM = max(MM, abs(A(j, 0)));
    MM = max(MM, abs(A(j, 1)));
    MM = max(MM, abs(A(j, 2)));
    MM = max(MM, abs(b[j]));
  }
  A  /= MM;
  b  /= MM;
  for(int j = 0; j < one.size(); j ++)
    std::cout << A(j, 0) << ", " << A(j, 1) << ", " << A(j, 2) << ", " << b[j] << std::endl;
  const auto r(pow(num_t(2), - num_t(abs(computer) * 4)));
  const auto err(A * inner<num_t>(A, b - one * r, b + one * r) - b);
        auto M(abs(err[0]));
  for(int j = 1; j < err.size(); j ++)
    M = max(M, abs(err[j]));
  std::cout << M << std::endl;
  // XXX: maxima's simplex result of M limit seems to 2/5 * |delta|.
  // XXX: konbu method result M limit seems to |delta|.
  
  return 0;
}

