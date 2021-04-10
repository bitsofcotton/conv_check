#include <cstdio>
#include <cstdlib>
#include "konbu_init.h"

#if defined(WITHOUT_EIGEN)
typedef SimpleMatrix<num_t> Mat;
typedef SimpleVector<num_t> Vec;
#endif

int main(int argc, char* argv[])
{
  cout.precision(MANT_DIG);
  cerr.precision(MANT_DIG);
  Mat A(7, 2);
  Vec b(A.rows());
  for(int j = 0; j < (b.size() + 1) / 2; j ++) {
    A(j, 0) = num_t((j & 1) ? 1 : 2);
    A(j, 1) = num_t(((j >> 1) & 1) ? 1 : 2);
    b[j]    = num_t(((j & 1) && ((j >> 1) & 1)) ? 2 : 1);
    const auto jj(j + (b.size() + 1) / 2);
    if(jj < b.size()) {
      A.row(jj) = - A.row(j);
      b[jj]     = - b[j];
    }
    std::cerr << A(j, 0) << ", " << A(j, 1) << ", " << b[j] << std::endl;
  }
  const auto err(A * Linner<num_t>().inner(A, b) - b);
        auto MM(err[0]);
  for(int j = 1; j < err.size(); j ++)
    MM = max(MM, err[j]);
  std::cout << MM << ", " << sqrt(b.dot(b)) << std::endl;
  
  return 0;
}

