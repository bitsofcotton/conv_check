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

int main(int argc, char* argv[])
{
  std::cout.precision(30);
  std::cerr.precision(30);
  SimpleMatrix<num_t> A(4, 4);
  SimpleVector<num_t> one(A.rows());
  const auto computer(- 20);
  for(int j = 0; j < 4; j ++) {
    A(j, 0) = num_t((j & 1) ? 1 : 2);
    A(j, 1) = num_t(((j >> 1) & 1) ? 1 : 2);
    A(j, 2) = num_t(((j & 1) && ((j >> 1) & 1)) ? 2 : 1);
    A(j, 3) = num_t(1);
    for(int i = 0; i < 3; i ++)
      A(j, i) = tan(A(j, i) / num_t(2) * atan2(num_t(1), num_t(1)));
    one[j]  = num_t(1);
    const auto r(pow(abs(A(j, 0) * A(j, 1) * A(j, 2) * A(j, 3)), num_t(computer)));
    A.row(j) *= r;
    A(j, 3)   = num_t(1);
  }
  for(int j = 0; j < one.size(); j ++)
    std::cout << A(j, 0) << ", " << A(j, 1) << ", " << A(j, 2) << ", " << A(j, 3) << std::endl;
  const auto in(A.inner(one, one));
  const auto err(A * in);
        auto M(abs(err[0]));
        auto m(abs(err[0]));
  for(int j = 1; j < err.size(); j ++) {
    M = max(M, abs(err[j]));
    m = min(m, abs(err[j]));
  }
  std::cout << m << ", " << M << ", " << sqrt(in.dot(in)) << std::endl;
  // result of this M, m, ||in|| limit seems to 1.
  return 0;
}

