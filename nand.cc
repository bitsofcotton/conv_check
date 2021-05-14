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

#if ! defined(_FLOAT_BITS_)
#define _FLOAT_BITS_ 63
#endif
#include <cmath>

int main(int argc, char* argv[]) {
  SimpleMatrix<num_t> A(7, 5);
  SimpleVector<num_t> left(A.rows());
  SimpleVector<num_t> right(A.rows());
  for(int j = 0; j < 4; j ++) {
    SimpleVector<num_t> work(4);
    work[0]  = num_t(j & 1);
    work[1]  = num_t((j >> 1) & 1);
    work[2]  = num_t(((j & 1) && ((j >> 1) & 1)) ? 0 : 1);
    work[3]  = (work[0] + work[1] + work[2]) / num_t(2) / num_t(3);
    A.row(j) = makeProgramInvariant<num_t>(work);
    assert(A.cols() == A.row(j).size());
    left[j]  = - (right[j] = sqrt(A.epsilon));
  }
  for(int i = 4; i < A.rows(); i ++) {
    for(int j = 0; j < A.cols(); j ++)
      A(i, j) = num_t(i - 4 == j ? 1 : 0);
    left[i] = (right[i] = num_t(1) / sqrt(sqrt(A.epsilon))) * sqrt(A.epsilon);
  }
  const auto in(A.inner(left, right));
  std::cout << A * in << in << std::endl;
  return 0;
}

