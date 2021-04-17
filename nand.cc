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
  SimpleMatrix<num_t> A(5, 4);
  SimpleVector<num_t> one(A.rows());
  assert(1 < argc);
  const auto computer(- abs(std::atoi(argv[1])));
  for(int j = 0; j < 4; j ++) {
    A(j, 0) = num_t((j & 1) ? 1 : 2);
    A(j, 1) = num_t(((j >> 1) & 1) ? 1 : 2);
    A(j, 2) = num_t(((j & 1) && ((j >> 1) & 1)) ? 2 : 1);
    A(j, 3) = num_t(1);
    for(int i = 0; i < 3; i ++)
      A(j, i) = tan(A(j, i) / num_t(4) * atan2(num_t(1), num_t(1)));
    one[j]  = num_t(1);
    const auto r(pow(abs(A(j, 0) * A(j, 1) * A(j, 2) * A(j, 3)), num_t(computer)));
    A.row(j) *= r;
    A(j, 3)   = num_t(1);
  }
  const auto j(4);
  A(j, 0) = A(j, 1) = A(j, 3) = one[j] = num_t(0);
  auto  mone(- one);
  A(j, 2) = mone[j] = num_t(1);
  one[j] = pow(num_t(2), num_t(abs(computer)));
  num_t m(0);
  for(int i = 0; i < A.rows(); i ++)
    for(int j = 0; j < A.cols(); j ++)
      if(m == num_t(0)) m = abs(A(i, j));
      else if(A(i, j) != num_t(0) && abs(A(i, j)) < m) m = abs(A(i, j));
  const auto in(A.inner(mone, one));
  const auto err(A * in);
        auto M(abs(err[0]) - num_t(1));
  for(int j = 1; j < err.size() - 1; j ++)
    M = max(M, abs(err[j]) - num_t(1));
  if(M < num_t(0)) M = num_t(0);
  std::cout << M / m << std::endl;
  std::cout << A << mone << one << in << std::endl;
  return 0;
}

