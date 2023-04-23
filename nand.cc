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
  SimpleMatrix<num_t> A(8 + 1, 4 + 1);
  SimpleVector<num_t> left(A.rows());
  SimpleVector<num_t> right(A.rows());
  for(int j = 0; j < 8; j ++) {
    SimpleVector<num_t> work(4);
    work[0]  = num_t(j & 1 ? 0 : 1);
    work[1]  = num_t((j >> 1) & 1 ? 0 : 1);
    work[2]  = num_t((j >> 2) & 1 ? 0 : 1);
    work[3]  = num_t(!((j & 1) && ((j >> 1) & 1)) && ((j >> 2) & 1) ? 1 : 0);
    A.row(j) = makeProgramInvariant<num_t>(work).first;
    left[j]  = - exp(- log(A.epsilon()) / num_t(int(2)) );
    right[j] = - exp(- log(A.epsilon()) / num_t(int(4)) );
  }
  for(int j = 0; j < A.cols(); j ++)
    A(8, j) = atan(num_t(3 == j ? 1 : 0));
  left[8]   = log(num_t(int(3)) / num_t(int(4)));
  right[8]  = log(num_t(int(89)) / num_t(int(90)) );
  auto in(A * A.inner(left, right));
  std::cout << left << right << A << std::endl;
  std::cout << in << in / sqrt(in.dot(in)) << std::endl;
  in = revertProgramInvariant<num_t>(in);
  std::cout << in << in / sqrt(in.dot(in)) << std::endl;
  return 0;
}

