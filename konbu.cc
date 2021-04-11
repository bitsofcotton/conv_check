#include <cstdio>
#include <cstdlib>
#include "konbu_init.h"

#define SZ_NUM_BUF 500
#include "konbu.hh"

#if defined(WITHOUT_EIGEN)
typedef SimpleMatrix<num_t> Mat;
typedef SimpleVector<num_t> Vec;
#endif

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

#if defined(ACC_GMP)
  num_t::set_default_prec(ACC_GMP);
#endif
  
  Mat A(N_INEQ, N_ELEM);
  Vec b(N_INEQ);
  for(int i = 0; i < N_INEQ; i ++)
    b[i] = ceil((num_t(1) / num_t(2) - num_t(random()) / num_t(RAND_MAX)) * num_t(M_VAL)) / num_t(MIN_VAL);
  for(int i = 0; i < N_ELEM; i ++)
    for(int j = 0; j < N_INEQ; j ++)
      A(j, i) = ceil((num_t(1) / num_t(2) - num_t(random()) / num_t(RAND_MAX)) * num_t(M_VAL)) / num_t(MIN_VAL);
  
  cout.precision(MANT_DIG);
  cerr.precision(MANT_DIG);
  
  b *= num_t(1000);
  const auto error(A * Linner<num_t>().inner(A, b) - b);
        auto M(error[0]);
  for(int i = 1; i < error.size(); i ++)
    M = max(M, error[i]);
  std::cout << M << ", " << sqrt(b.dot(b)) << std::endl;
  return 0;
}

