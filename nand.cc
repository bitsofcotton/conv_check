#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>

#include "lieonn.hh"
typedef myfloat num_t;

int main(int argc, char* argv[])
{
  std::cout.precision(30);
  std::cerr.precision(30);
  SimpleMatrix<num_t> A(4, 3);
  SimpleVector<num_t> b(A.rows());
  SimpleVector<num_t> one(A.rows());
  for(int j = 0; j < 4; j ++) {
    A(j, 0) = num_t((j & 1) ? 1 : 2);
    A(j, 1) = num_t(((j >> 1) & 1) ? 1 : 2);
    A(j, 2) = num_t(((j & 1) && ((j >> 1) & 1)) ? 2 : 1);
    b[j]    = num_t(1);
    for(int i = 0; i < A.cols(); i ++)
      A(j, i) = tan(A(j, i) / num_t(4) * atan2(num_t(1), num_t(1)));
    b[j] = tan(b[j] / num_t(4) * atan2(num_t(1), num_t(1)));
    one[j]  = num_t(1);
  }
  const auto eps(pow(num_t(2), - num_t(8)));
  const auto in(A.inner(b - one * eps, b + one * eps));
  const auto err(A * in);
        auto M(abs(err[0] - b[0]));
  for(int i = 1; i < err.size(); i ++)
    M = max(M, abs(err[i] - b[i]));
  std::cout << M << std::endl;
  std::cout << A << b << in << std::endl;
  return 0;
}

