#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <assert.h>
#include <stdlib.h>
#include "lieonn.hh"

int main(int argc, char* argv[]) {
  SimpleMatrix<myfloat> A(12, 8);
  for(int i = 0; i < A.rows(); i ++)
    for(int j = 0; j < A.cols(); j ++)
      A(i, j) = myfloat(int(arc4random())) / myfloat(int(arc4random()));
  const auto svd(A.SVD());
        auto rsvd(svd * A);
  for(int i = 0; i < rsvd.rows(); i ++)
    rsvd.row(i) /= sqrt(rsvd.row(i).dot(rsvd.row(i)));
  std::cout <<  svd *  svd.transpose();
  std::cout << rsvd * rsvd.transpose();
  std::cout <<  svd.transpose() * svd;
  std::cout << rsvd.transpose() * rsvd;
  return 0;
}

