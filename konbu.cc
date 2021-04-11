#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
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

#define MIN_VAL 1000
#define M_VAL   1000000
int main(int argc, char* argv[])
{
  int N_ELEM = 0, N_INEQ = 0;  
  srandomdev();
  // XXX insecure:
  // srand(time(0));
  N_ELEM = int(ceil(num_t(random()) / num_t(RAND_MAX) * num_t(40)) + num_t(40));
  N_INEQ = int(num_t(N_ELEM) + ceil(num_t(random()) / num_t(RAND_MAX) * num_t(20)));

  SimpleMatrix<num_t> A(N_INEQ, N_ELEM);
  SimpleVector<num_t> b(N_INEQ);
  for(int i = 0; i < N_INEQ; i ++)
    b[i] = abs(ceil((num_t(1) / num_t(2) - num_t(random()) / num_t(RAND_MAX)) * num_t(M_VAL)) / num_t(MIN_VAL));
  for(int i = 0; i < N_ELEM; i ++)
    for(int j = 0; j < N_INEQ; j ++)
      A(j, i) = ceil((num_t(1) / num_t(2) - num_t(random()) / num_t(RAND_MAX)) * num_t(M_VAL)) / num_t(MIN_VAL);
  
  std::cout.precision(30);
  std::cerr.precision(30);
  
  b *= num_t(1000);
  const auto error(A * inner<num_t>(A, b * num_t(0), b) - b);
        auto M(error[0]);
  for(int i = 1; i < error.size(); i ++)
    M = max(M, error[i]);
  std::cout << M << ", " << sqrt(b.dot(b)) << std::endl;
  return 0;
}

