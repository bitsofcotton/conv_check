#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>

#include "ifloat.hh"
typedef myfloat num_t;
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
  const auto in(A.inner(one, one));
  const auto err(A * in);
        auto M(abs(err[0]));
        auto m(abs(err[0]));
  for(int j = 1; j < err.size(); j ++) {
    M = max(M, abs(err[j]));
    m = min(m, abs(err[j]));
  }
  std::cout << A << in << std::endl;
  // result of this M, m, ||in|| limit seems to 1.
  return 0;
}

