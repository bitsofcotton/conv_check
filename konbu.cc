#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>
#include <sys/resource.h>

#include "lieonn.hh"
typedef myfloat num_t;

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

  SimpleMatrix<num_t> A(N_INEQ, N_ELEM);
  SimpleVector<num_t> b(N_INEQ);
  for(int i = 0; i < N_INEQ; i ++)
    b[i] = abs(ceil((num_t(1) / num_t(2) - num_t(random()) / num_t(RAND_MAX)) * num_t(M_VAL)) / num_t(MIN_VAL));
  for(int i = 0; i < N_ELEM; i ++)
    for(int j = 0; j < N_INEQ; j ++)
      A(j, i) = ceil((num_t(1) / num_t(2) - num_t(random()) / num_t(RAND_MAX)) * num_t(M_VAL)) / num_t(MIN_VAL);
  
  std::cout.precision(30);
  std::cerr.precision(30);
  struct rusage start, end;
  getrusage(RUSAGE_SELF, &start);
  const auto err(A * A.inner(b * num_t(0), b) - b);
  getrusage(RUSAGE_SELF, &end);
        auto M(err[0]);
  for(int i = 1; i < err.size(); i ++)
    M = max(M, err[i]);
  std::cout << A << b << std::endl;
  std::cout << M << ", " << sqrt(b.dot(b)) << std::endl;
  std::cout << end.ru_utime.tv_sec - start.ru_utime.tv_sec << "[s] and ";
  std::cout << end.ru_utime.tv_usec - start.ru_utime.tv_usec << "[mu s]";
  std::cout << std::endl << std::endl;
  return 0;
}

