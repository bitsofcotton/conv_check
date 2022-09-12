#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <assert.h>
#include <stdlib.h>
#include "lieonn.hh"

int main(int argc, char* argv[]) {
  assert(2 < argc);
  SimpleMatrix<myfloat> A(std::atoi(argv[1]), std::atoi(argv[2]));
  for(int i = 0; i < A.rows(); i ++)
    for(int j = 0; j < A.cols(); j ++)
      A(i, j) = myfloat(int(arc4random())) - myfloat(int(RAND_MAX)) / myfloat(int(2));
  std::cout << A.SVD() * A * A.transpose().SVD().transpose();
  return 0;
}

