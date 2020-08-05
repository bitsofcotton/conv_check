#include <cstdio>
#include <cstdlib>
#include "konbu_init.h"

#define SZ_NUM_BUF 500
#include "konbu.hh"

#if defined(WITHOUT_EIGEN)
typedef SimpleMatrix<num_t> Mat;
typedef SimpleVector<num_t> Vec;
#endif

#define MIN_VAL 100
#define M_VAL   10000
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
    b[i] = ceil((num_t(1) / num_t(2) - random() / num_t(RAND_MAX)) * M_VAL) / MIN_VAL;
  for(int i = 0; i < N_ELEM; i ++)
    for(int j = 0; j < N_INEQ; j ++)
      A(j, i) = ceil((num_t(1) / num_t(2) - random() / num_t(RAND_MAX)) * M_VAL) / MIN_VAL;
  
  cout.precision(MANT_DIG);
  cerr.precision(MANT_DIG);
  
  Vec  result;
  bool* fix_partial = new bool[A.rows()];
  Linner<num_t> lp;
  bool feas = lp.inner(fix_partial, result, A, b);
  
  if(feas)
    cout << " OPT " << endl;
  
  int n_fixed = 0;
  for(int i = 0; i < A.rows(); i ++)
    if(fix_partial[i]) n_fixed ++;
  cout << "fix_partial ( " << n_fixed << " / " << A.cols() << " ) : ";
  for(int i = 0; i < A.rows(); i ++)
    cout << (const char*)(fix_partial[i] ? "1" : "0");
  cout << endl;
  for(int i = 0; i < result.size(); i ++)
    cout << result[i] << "\t";
  cout << std::endl;
  
  delete[] fix_partial;
  
  return 0;
}

