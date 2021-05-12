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
  SimpleMatrix<num_t> A(5, 4);
  SimpleVector<num_t> left(A.rows());
  SimpleVector<num_t> right(A.rows());
  for(int j = 0; j < 4; j ++) {
    SimpleVector<num_t> work(3);
    work[0]  = num_t(j & 1);
    work[1]  = num_t((j >> 1) & 1);
    work[2]  = num_t(((j & 1) && ((j >> 1) & 1)) ? 0 : 1);
    A.row(j) = makeProgramInvariant<num_t>(work);
    left[j]  = - num_t(100);
    right[j] =   num_t(100);
  }
  A(4, 0) = A(4, 1) = A(4, 3) = num_t(0);
  A(4, 2) = num_t(1);
  left[4] = num_t(1) / num_t(2);
  right[4] = num_t(2);
  const auto in(A.inner(left, right));
  const auto errl(left - A * in);
  const auto errr(A * in - right);
        auto M(max(errl[0], errr[0]));
  for(int i = 1; i < errl.size(); i ++)
    M = max(M, max(errl[i], errr[i]));
  std::cout << M << std::endl;
  std::cout << A << left << right << in << std::endl;
  return 0;
}

