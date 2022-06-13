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
  SimpleVector<num_t> r(A.rows());
  for(int j = 0; j < 8; j ++) {
    SimpleVector<num_t> work(4);
    work[0]  = num_t(j & 1 ? 0 : 1);
    work[1]  = num_t((j >> 1) & 1 ? 0 : 1);
    work[2]  = num_t((j >> 2) & 1 ? 0 : 1);
    work[3]  = num_t(!((j & 1) && ((j >> 1) & 1)) && ((j >> 2) & 1) ? 1 : 0);
    const auto buf(makeProgramInvariant<num_t>(work));
    A.row(j) = buf.first * (r[j] = pow(abs(buf.second), ceil(- log(A.epsilon()) / num_t(int(2)) ) ));
    right[j] = - (left[j]  = - exp(- sqrt(- log(A.epsilon()) )) );
  }
  const auto rr(sqrt(sqrt(sqrt(r[0] * r[1] * r[2] * r[3] * r[4] * r[5] * r[6] * r[7]))));
  A /= rr;
  for(int j = 0; j < A.cols(); j ++)
    A(8, j) = num_t(3 == j ? 1 : 0);
  r[8] = num_t(1);
  left[8]  = sqrt(right[0]);
  right[8] = num_t(1);
  auto in(A * A.inner(left, right));
  std::cout << left << right << A << std::endl;
  std::cout << in << in / sqrt(in.dot(in)) << std::endl;
  for(int i = 0; i < in.size(); i ++)
    if(r[i] != num_t(0))
      in[i] *= rr / r[i];
  std::cout << in << in / sqrt(in.dot(in)) << std::endl;
  return 0;
}

