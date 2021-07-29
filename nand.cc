#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>

#include "lieonn.hh"
typedef myfloat num_t;

#if ! defined(_FLOAT_BITS_)
#define _FLOAT_BITS_ 63
#endif
#include <cmath>

int main(int argc, char* argv[]) {
  SimpleMatrix<num_t> A(8 + 1, 6);
  SimpleVector<num_t> left(A.rows());
  SimpleVector<num_t> right(A.rows());
  for(int j = 0; j < 8; j ++) {
    SimpleVector<num_t> work(4);
    work[0]  = num_t(j & 1 ? 1 : 2) / num_t(2);
    work[1]  = num_t((j >> 1) & 1 ? 1 : 2) / num_t(2);
    work[2]  = num_t((j >> 2) & 1 ? 1 : 2) / num_t(2);
    work[3]  = num_t(!((j & 1) && ((j >> 1) & 1)) && ((j >> 2) & 1) ? 2 : 1) / num_t(2);
    A.row(j) = makeProgramInvariant<num_t>(work).first;
    assert(A.cols() == A.row(j).size());
    left[j]  = - (right[j] = sqrt(A.epsilon));
  }
  for(int j = 0; j < A.cols(); j ++)
    A(8, j) = num_t(3 == j ? 1 : 0);
  left[8]  = sqrt(sqrt(A.epsilon));
  right[8] = num_t(1) / sqrt(sqrt(A.epsilon));
  const auto in(A.inner(left, right));
  std::cout << A * in << std::endl;
  return 0;
}

