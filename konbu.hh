/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/* BSD 3-Clause License:
 * Copyright (c) 2013 - 2017, kazunobu watatsu.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation or other materials provided with the distribution.
 *    Neither the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if !defined(_LINEAR_OPT_)

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>

#if defined(WITHOUT_EIGEN)
#include <assert.h>
#else
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;
#endif

using namespace std;

template <typename T> class SimpleVector {
public:
  SimpleVector();
  SimpleVector(const int& size);
  SimpleVector(const SimpleVector<T>& other);
  ~SimpleVector();
  
        SimpleVector<T>  operator -  () const;
        SimpleVector<T>  operator +  (const SimpleVector<T>& other) const;
  const SimpleVector<T>& operator += (const SimpleVector<T>& other);
        SimpleVector<T>  operator -  (const SimpleVector<T>& other) const;
  const SimpleVector<T>& operator -= (const SimpleVector<T>& other);
        SimpleVector<T>  operator *  (const T& other) const;
  const SimpleVector<T>& operator *= (const T& other);
        SimpleVector<T>  operator /  (const T& other) const;
  const SimpleVector<T>& operator /= (const T& other);
        SimpleVector<T>& operator =  (const SimpleVector<T>& other);
        T                dot         (const SimpleVector<T>& other) const;
        T&               operator [] (const int& idx);
  const T                operator [] (const int& idx) const;
  const int& size() const;
private:
  T*  entity;
  int esize;
};

template <typename T> SimpleVector<T>::SimpleVector() {
  entity = NULL;
  esize  = 0;
}

template <typename T> SimpleVector<T>::SimpleVector(const int& size) {
  assert(size > 0);
  this->entity = new T[size];
  this->esize  = size;
  return;
}

template <typename T> SimpleVector<T>::SimpleVector(const SimpleVector<T>& other) {
  entity = new T[other.esize];
  esize  = other.esize;
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] = other.entity[i];
  return;
}

template <typename T> SimpleVector<T>::~SimpleVector() {
  if(entity)
    delete[] entity;
  return;
}

template <typename T> SimpleVector<T> SimpleVector<T>::operator - () const {
  SimpleVector<T> res(esize);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    res.entity[i] = - entity[i];
  return res;
}

template <typename T> SimpleVector<T> SimpleVector<T>::operator + (const SimpleVector<T>& other) const {
  SimpleVector<T> res(*this);
  return res += other;
}

template <typename T> const SimpleVector<T>& SimpleVector<T>::operator += (const SimpleVector<T>& other) {
  assert(esize == other.esize);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] += other.entity[i];
  return *this;
}

template <typename T> SimpleVector<T> SimpleVector<T>::operator - (const SimpleVector<T>& other) const {
  SimpleVector<T> res(*this);
  return res -= other;
}

template <typename T> const SimpleVector<T>& SimpleVector<T>::operator -= (const SimpleVector<T>& other) {
  return *this += - other;
}

template <typename T> SimpleVector<T> SimpleVector<T>::operator * (const T& other) const {
  SimpleVector<T> res(*this);
  return res *= other;
}

template <typename T> const SimpleVector<T>& SimpleVector<T>::operator *= (const T& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] *= other;
  return *this;
}

template <typename T> SimpleVector<T> SimpleVector<T>::operator / (const T& other) const {
  SimpleVector<T> res(*this);
  return res /= other;
}

template <typename T> SimpleVector<T>& SimpleVector<T>::operator = (const SimpleVector<T>& other) {
  if(entity == other.entity && esize == other.esize)
    return *this;
  if(esize != other.esize) {
    if(entity)
      delete[] entity;
    entity = new T[other.esize];
  }
  esize = other.esize;
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] = other.entity[i];
  return *this;
}

template <typename T> const SimpleVector<T>& SimpleVector<T>::operator /= (const T& other) {
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    entity[i] /= other;
  return *this;
}

template <typename T> T SimpleVector<T>::dot(const SimpleVector<T>& other) const {
  assert(esize == other.esize);
  T res(0);
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    res += entity[i] * other.entity[i];
  return res;
}

template <typename T> T& SimpleVector<T>::operator [] (const int& idx) {
  assert(0 <= idx && idx < esize && entity);
  return entity[idx];
}

template <typename T> const T SimpleVector<T>::operator [] (const int& idx) const {
  assert(0 <= idx && idx < esize && entity);
  return entity[idx];
}

template <typename T> const int& SimpleVector<T>::size() const {
  return esize;
}

template <typename T> class SimpleMatrix {
public:
  SimpleMatrix();
  SimpleMatrix(const int& rows, const int& cols);
  SimpleMatrix(const SimpleMatrix<T>& other);
  ~SimpleMatrix();
  
        SimpleMatrix<T>  operator -  () const;
        SimpleMatrix<T>  operator +  (const SimpleMatrix<T>& other) const;
  const SimpleMatrix<T>& operator += (const SimpleMatrix<T>& other);
        SimpleMatrix<T>  operator -  (const SimpleMatrix<T>& other) const;
  const SimpleMatrix<T>& operator -= (const SimpleMatrix<T>& other);
        SimpleMatrix<T>  operator *  (const T& other) const;
  const SimpleMatrix<T>& operator *= (const T& other);
        SimpleMatrix<T>  operator *  (const SimpleMatrix<T>& other) const;
  const SimpleMatrix<T>& operator *= (const SimpleMatrix<T>& other);
        SimpleVector<T>  operator *  (const SimpleVector<T>& other) const;
        SimpleMatrix<T>  operator /  (const T& other) const;
  const SimpleMatrix<T>& operator /= (const T& other);
        SimpleMatrix<T>& operator =  (const SimpleMatrix<T>& other);
        T&               operator () (const int& y, const int& x);
  const T                operator () (const int& y, const int& x) const;
        SimpleVector<T>& row(const int& y);
  const SimpleVector<T>  row(const int& y) const;
  const SimpleVector<T>  col(const int& x) const;
        void             setCol(const int& x, const SimpleVector<T>& other);
        SimpleMatrix<T>  transpose() const;
        SimpleVector<T>  solve(SimpleVector<T> other) const;
  const int& rows() const;
  const int& cols() const;
private:
  SimpleVector<T>* entity;
  int              erows;
  int              ecols;
};

template <typename T> SimpleMatrix<T>::SimpleMatrix(const int& rows, const int& cols) {
  assert(rows > 0 && cols > 0);
  entity = new SimpleVector<T>[rows];
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < rows; i ++)
    entity[i] = SimpleVector<T>(cols);
  erows = rows;
  ecols = cols;
  return; 
}

template <typename T> SimpleMatrix<T>::SimpleMatrix(const SimpleMatrix<T>& other) {
  erows = other.erows;
  ecols = other.ecols;
  entity = new SimpleVector<T>[other.erows];
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] = other.entity[i];
  return;
}

template <typename T> SimpleMatrix<T>::~SimpleMatrix() {
  if(entity)
    delete[] entity;
  return;
}
  
template <typename T> SimpleMatrix<T> SimpleMatrix<T>::operator - () const {
  SimpleMatrix<T> res(erows, ecols);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++)
    res.entity[i] = - entity[i];
  return res;
}

template <typename T> SimpleMatrix<T> SimpleMatrix<T>::operator + (const SimpleMatrix<T>& other) const {
  SimpleMatrix<T> res(*this);
  return res += other;
}

template <typename T> const SimpleMatrix<T>& SimpleMatrix<T>::operator += (const SimpleMatrix<T>& other) {
  assert(erows = other.erows && ecols == other.ecols);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] += other.entity[i];
  return *this;
}

template <typename T> SimpleMatrix<T> SimpleMatrix<T>::operator - (const SimpleMatrix<T>& other) const {
  SimpleMatrix<T> res(*this);
  return res -= other;
}

template <typename T> const SimpleMatrix<T>& SimpleMatrix<T>::operator -= (const SimpleMatrix<T>& other) {
  return *this += - other;
}

template <typename T> SimpleMatrix<T> SimpleMatrix<T>::operator * (const T& other) const {
  SimpleMatrix<T> res(*this);
  return res *= other;
}

template <typename T> const SimpleMatrix<T>& SimpleMatrix<T>::operator *= (const T& other) {
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] *= other;
  return *this;
}

template <typename T> SimpleMatrix<T> SimpleMatrix<T>::operator * (const SimpleMatrix<T>& other) const {
  assert(ecols == other.erows && entity && other.entity);
  SimpleMatrix<T> derived(other.transpose());
  SimpleMatrix<T> res(erows, other.ecols);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++) {
          SimpleVector<T>& resi(res.entity[i]);
    const SimpleVector<T>& ei(entity[i]);
    for(int j = 0; j < other.ecols; j ++)
      resi[j] = ei.dot(derived.entity[j]);
  }
  return res;

}

template <typename T> const SimpleMatrix<T>& SimpleMatrix<T>::operator *= (const SimpleMatrix<T>& other) {
  return *this = *this * other;
}

template <typename T> SimpleVector<T> SimpleMatrix<T>::operator * (const SimpleVector<T>& other) const {
  assert(ecols == other.size());
  SimpleVector<T> res(erows);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++)
    res[i] = entity[i].dot(other);
  return res;
}

template <typename T> SimpleMatrix<T> SimpleMatrix<T>::operator / (const T& other) const {
  SimpleMatrix<T> res(*this);
  return res /= other;
}

template <typename T> const SimpleMatrix<T>& SimpleMatrix<T>::operator /= (const T& other) {
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] /= other;
  return *this;
}

template <typename T> SimpleMatrix<T>& SimpleMatrix<T>::operator = (const SimpleMatrix<T>& other) {
  if(entity == other.entity && erows == other.erows && ecols == other.ecols)
    return *this;
  if(erows != other.erows || ecols != other.ecols) {
    if(entity)
      delete[] entity;
    entity = new SimpleVector<T>[other.erows];
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for(int i = 0; i < other.erows; i ++)
      entity[i] = SimpleVector<T>(other.ecols);
  }
  erows = other.erows;
  ecols = other.ecols;
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++)
    entity[i] = other.entity[i];
  return *this;

}

template <typename T> T& SimpleMatrix<T>::operator () (const int& y, const int& x) {
  assert(0 <= y && y < erows && entity);
  return entity[y][x];
}

template <typename T> const T SimpleMatrix<T>::operator () (const int& y, const int& x) const {
  assert(0 <= y && y < erows && entity);
  return entity[y][x];
}

template <typename T> SimpleVector<T>& SimpleMatrix<T>::row(const int& y) {
  assert(0 <= y && y < erows && entity);
  return entity[y];
}

template <typename T> const SimpleVector<T> SimpleMatrix<T>::row(const int& y) const {
  assert(0 <= y && y < erows && entity);
  return entity[y];
}

template <typename T> const SimpleVector<T> SimpleMatrix<T>::col(const int& x) const {
  assert(0 <= erows && 0 <= x && x < ecols && entity);
  SimpleVector<T> res(erows);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++)
    res[i] = entity[i][x];
  return res;
}

template <typename T> void SimpleMatrix<T>::setCol(const int& x, const SimpleVector<T>& other) {
  assert(0 <= x && x < ecols && other.size() == erows);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < erows; i ++)
    entity[i][x] = other[i];
  return;
}

template <typename T> SimpleMatrix<T> SimpleMatrix<T>::transpose() const {
  SimpleMatrix<T> res(ecols, erows);
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < ecols; i ++) {
    SimpleVector<T>& resi(res.entity[i]);
    for(int j = 0; j < erows; j ++)
      resi[j] = entity[j][i];
  }
  return res;
}

template <typename T> SimpleVector<T> SimpleMatrix<T>::solve(SimpleVector<T> other) const {
  assert(0 <= erows && 0 <= ecols && erows == ecols && entity && erows == other.size());
  SimpleMatrix<T> work(*this);
  for(int i = 0; i < erows; i ++) {
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for(int j = i + 1; j < erows; j ++) {
      const T ratio(work.entity[j][i] / work.entity[i][i]);
      work.entity[j] -= work.entity[i] * ratio;
      other[j]       -= other[i]       * ratio;
    }
  }
  for(int i = erows - 1; 0 <= i; i --) {
    other[i]         /= work.entity[i][i];
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for(int j = i - 1; 0 <= j; j --)
      other[j]       -= other[i] * work.entity[j][i];
  }
  return other;
}

template <typename T> const int& SimpleMatrix<T>::rows() const {
  return erows;
}

template <typename T> const int& SimpleMatrix<T>::cols() const {
  return ecols;
}


template <typename T> class LP {
public:
#if defined(WITHOUT_EIGEN)
  typedef SimpleMatrix<T> Mat;
  typedef SimpleVector<T> Vec;
#else
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
#endif
  
  LP();
  ~LP();
  
  const T err_norm(const T& err) const;
  
  // <c, x> -> obj, Ax <= b.
  bool minimize(bool* fix_partial, const Mat& A, const Vec& b, const Vec& c) const;
  bool minimize(Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const;
  bool minimize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const;
  bool maximize(bool* fix_partial, const Mat& A, const Vec& b, const Vec& c) const;
  bool maximize(Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const;
  bool maximize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const;
  bool inner(bool* fix_partial, const Mat& A, const Vec& b) const;
  bool inner(Vec& rvec, const Mat& A, const Vec& b) const;
  bool inner(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b) const;
  
  bool isErrorMargin(const Mat& A, const Vec& b, const Vec& x, const bool& disp = true) const;
  
private:
  bool optimize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, Vec c) const;
  bool optimizeNullSpace(bool* fix_partial, Vec& rvec, const Mat& P, Vec& b, const int& cidx) const;
  
  void normalizeRows(Mat& A, Vec& b) const;
  void normalizeCols(Vec& acr, T& bcr, Mat& A, Vec& b) const;
  T    errorCheck(const Mat& A, const Vec& b) const;
  
  bool gainVectors(bool* fix, bool* checked, Vec& rvec, const Mat& P, const Vec& b, const Vec& one, int& n_fixed) const;
  Vec  giantStep(bool* fix, bool* checked, Mat Pverb, Vec mb, int& n_fixed, const Vec& one) const;
  bool checkInner(const Vec& on, const Vec& normalize, const int& max_idx) const;
  int  getMax(const bool* checked, const Vec& on, const Vec& normalize) const;
  
  Mat roughQR(const Mat& A) const;

  T err_raw_epsilon;
  T threshold_feas;
  T threshold_p0;
  T threshold_loop;
  T threshold_inner;
  T err_error;
  T largest_intercept;
  T largest_opt;
};

template <typename T> LP<T>::LP()
{
  err_raw_epsilon   = T(1);
  while(T(1) + err_raw_epsilon != T(1)) err_raw_epsilon /= T(2);
  err_raw_epsilon  *= T(2);
  threshold_feas    = err_raw_epsilon * (- log(err_raw_epsilon) / log(2));
  threshold_p0      = sqrt(threshold_feas);
  threshold_inner   = pow(threshold_p0, T(1) / T(4));
  threshold_loop    = sqrt(threshold_p0 * threshold_inner);
  err_error         = sqrt(threshold_inner);
  if(err_error >= T(1) - err_raw_epsilon)
    cerr << " accuracy not enough in initializing : " << err_error << endl;
  
  // finite maximum 1/abs(cos(theta(x,a_k))) value depends on the problem.
  // If this will larger, feasibility of the problem to solve increases,
  // but if too large, gainVector loop fails because lack of accuracy b.
  largest_intercept = T(2);
  
  // maximum abs(objective_value) ratio, if too large, gainVector fails.
  largest_opt       = pow(threshold_feas, - T(1) / T(4));
  return;
}

template <typename T> LP<T>::~LP()
{
  return;
}

template <typename T> const T LP<T>::err_norm(const T& err) const
{
  // margin needed for stable calculation.
  return err;
}

template <typename T> bool LP<T>::minimize(bool* fix_partial, const Mat& A, const Vec& b, const Vec& c) const
{
  Vec rvec;
  return minimize(fix_partial, rvec, A, b, c);
}

template <typename T> bool LP<T>::minimize(Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const
{
  return minimize(0, rvec, A, b, c);
}

template <typename T> bool LP<T>::minimize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const
{
  return optimize(fix_partial, rvec, A, b, c);
}

template <typename T> bool LP<T>::maximize(bool* fix_partial, const Mat& A, const Vec& b, const Vec& c) const
{
  Vec rvec;
  return maximize(fix_partial, rvec, A, b, c);
}

template <typename T> bool LP<T>::maximize(Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const
{
  return maximize(0, rvec, A, b, c);
}

template <typename T> bool LP<T>::maximize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const
{
  return optimize(fix_partial, rvec, A, b, - c);
}

template <typename T> bool LP<T>::inner(bool* fix_partial, const Mat& A, const Vec& b) const
{
  Vec rvec;
  return inner(fix_partial, rvec, A, b);
}

template <typename T> bool LP<T>::inner(Vec& rvec, const Mat& A, const Vec& b) const
{
  return inner(0, rvec, A, b);
}

template <typename T> bool LP<T>::inner(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b) const
{
  Vec c(A.cols());
  for(int i = 0; i < c.size(); i ++)
    c[i] = T(0);
  return optimize(fix_partial, rvec, A, b, c);
}

template <typename T> bool LP<T>::optimize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, Vec c) const
{
  // cout << A << endl << b << endl << c.transpose() << endl;
  cerr << " (" << A.rows() << ", " << A.cols() << ")";
  rvec = Vec(c.size());
  for(int i = 0; i < A.rows(); i ++)
    if(fix_partial) fix_partial[i] = false;
  if(A.cols() != c.size() || A.rows() != b.size() || c.size() <= 0 || b.size() <= 0) {
    cerr << " Optimize: internal error (form does not match: Ax <= b, <c,x>.)" << endl;
    for(int i = 0; i < rvec.size(); i ++)
      rvec[i] = T(0);
    return false;
  }
  
  Mat AA(A.rows() + 1 + A.cols() * 2, A.cols());
  Vec bb(b.size() + 1 + A.cols() * 2);
  for(int i = 0; i < A.rows(); i ++) {
    AA.row(i) = A.row(i);
    bb[i]     = b[i];
  }
  
  // alpha t - <c, tx> <= 0 <=> alpha <= <c, x> with scaled x.
#if WITHOUT_EIGEN
  AA.row(A.rows()) = c;
#else
  AA.row(A.rows()) = c.transpose();
#endif
  bb[A.rows()]     = T(0);
  for(int i = 0; i < A.cols() * 2; i ++) {
    for(int j = 0; j < A.cols(); j ++)
      AA(A.rows() + 1 + i, j) = T(0);
    bb[A.rows() + 1 + i]      = T(0);
  }
  
  // scale.
  Vec acr(A.cols());
  T   bcr(1);
  for(int i = 0; i < acr.size(); i ++)
    acr[i] = T(1);
  {
    Mat AAA(AA), AAA2(AA);
    Vec bbb(bb), bbb2(bb);
    Vec aacr(acr), aacr2(acr);
    T   bbcr(bcr), bbcr2(bcr);
    normalizeRows(AAA,   bbb);
    normalizeCols(aacr,  bbcr,  AAA,  bbb);
    normalizeCols(aacr2, bbcr2, AAA2, bbb2);
    normalizeRows(AAA2,  bbb2);
    if(errorCheck(AAA,   bbb) < errorCheck(AAA2, bbb2)) {
      AA  = AAA2;
      bb  = bbb2;
      acr = aacr2;
      bcr = bbcr2;
    } else {
      AA  = AAA;
      bb  = bbb;
      acr = aacr;
      bcr = bbcr;
    }
  }
  normalizeRows(AA, bb);
  
  // check accuracy.
  const T errr(errorCheck(AA, bb));
  if(errr < T(1))
    cerr << " err_error?(" << errr << ")";
  
  // guarantee that b is positive.
  const T normbb(sqrt(bb.dot(bb)));
  for(int i = 0; i < A.cols() * 2; i ++) {
    AA(A.rows() + 1 + i, i / 2) = (i % 2 == 0 ? T(1) : - T(1));
    // XXX fixme:
    //bb[A.rows() + 1 + i] = normbb * T(A.cols()) * largest_intercept;
    bb[A.rows() + 1 + i] = normbb * largest_intercept;
  }
  
  bool fflag = false;
  for(int i = 0; i < AA.rows() && !fflag; i ++) {
    if(!isfinite(bb[i])) fflag = true;
    for(int j = 0; j < AA.cols(); j ++)
      if(!isfinite(AA(i, j))) fflag = true;
  }
  if(fflag) {
    cerr << " NaN : infeasible." << endl;
    return false;
  }
  cerr << " NORMALIZE";
  fflush(stderr);
  
  // <c, x> -> obj, Q R x <= b.
  // O(mn * max(m, n))
/*
  Eigen::HouseholderQR<Mat> qr;
  qr.compute(AA);
  Mat Qbuf(qr.householderQ());
  Mat Q(AA.rows(), AA.cols());
  for(int i = 0; i < Q.cols(); i ++)
    Q.col(i) = Qbuf.col(i);
*/
  Mat Q(roughQR(AA));
  Mat R(Q.transpose() * AA);
  cerr << " QR";
  fflush(stderr);

  // optimize <co_c, x'>, co_A [x' 0] <= b. co_A^t co_A = I, in O(mn^2L).
  bool* bfix_partial = new bool[Q.rows()];
  for(int i = 0; i < Q.rows(); i ++)
    bfix_partial[i] = false;
  (void)optimizeNullSpace(bfix_partial, rvec, Q, bb, A.rows());
  if(fix_partial)
    for(int i = 0; i < A.rows(); i ++)
      fix_partial[i] = bfix_partial[i];
  delete[] bfix_partial;
  
  cerr << " NULL_SPACE";
  fflush(stderr);
  
  // pull original solvee.
#if defined(WITHOUT_EIGEN)
  rvec = R.solve(rvec) * bcr;
#else
  rvec = R.inverse() * rvec * bcr;
#endif

  for(int i = 0; i < acr.size(); i ++)
    if(acr[i] != T(0))
      rvec[i] /= acr[i];
  
  cerr << " INVERT" << endl;
  return isErrorMargin(A, b, rvec, true);
}

template <typename T> bool LP<T>::optimizeNullSpace(bool* fix_partial, Vec& rvec, const Mat& P, Vec& b, const int& cidx) const
{
  if(P.rows() != b.size() || P.cols() <= 0 || P.rows() < P.cols()) {
    cerr << " NullSpace: internal error (form does not match: P[x 0] <= b, <c,x>.)" << endl;
    return false;
  }
  
  const T normb(sqrt(b.dot(b)));
  // to avoid stack underflow.
  bool*   checked = new bool[b.size()];
  Vec     one(b.size());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
  for(int i = 0; i < one.size(); i ++)
    one[i] = T(1);
  
  int  n_fixed;
  if(P.row(cidx).dot(P.row(cidx)) <= T(0)) {
    (void)gainVectors(fix_partial, checked, rvec, P, b, one, n_fixed);
    delete[] checked;
    return isErrorMargin(P, b, rvec, true);
  }
  
  // get optimal value.
  //  note: for this step, dynamic accuracy detecting will earn calculation time.
  T mid(1);
  T offset(normb * largest_opt);
  T step(offset * T(2));
  T b_proper(offset);
  int n_steps = 3 + 2 * int(log(- T(2) * log(err_raw_epsilon) / log(T(2))) / log(T(2)));
  for(int i = 0; i <= n_steps; i ++) {
    b[cidx] = - step * mid + offset;
    if(!isfinite(b[cidx]))
      continue;
    if(gainVectors(fix_partial, checked, rvec, P, b, one, n_fixed)) {
      b_proper = b[cidx];
      mid     *= step;
      cerr << ".";
    } else
      cerr << "!";
    step    = sqrt(step);
  }
  cerr << endl;
  
  b[cidx] = b_proper;
  (void)gainVectors(fix_partial, checked, rvec, P, b, one, n_fixed);
  delete[] checked;
  return isErrorMargin(P, b, rvec, true);
}

template <typename T> void LP<T>::normalizeRows(Mat& A, Vec& b) const
{
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < A.rows(); i ++) {
    T min(abs(A(i, 0)));
    T max(min);
    for(int j = 1; j < A.cols(); j ++) {
      const T buf(abs(A(i, j)));
      if(max < buf)
        max = buf;
      if((buf < min || min <= T(0)) && buf != T(0))
        min = buf;
    }
    const T buf(abs(b[i]));
    if(max < buf)
      max = buf;
    if((buf < min || min <= T(0)) && buf != T(0))
      min = buf;
    const T ratio(sqrt(min * max));
    if(ratio > T(0)) {
      A.row(i) /= ratio;
      b[i]     /= ratio;
    }
  }
  return;
}

template <typename T> void LP<T>::normalizeCols(Vec& acr, T& bcr, Mat& A, Vec& b) const
{
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < A.cols(); i ++) {
    T min(abs(A(0, i)));
    T max(min);
    for(int j = 1; j < A.rows(); j ++) {
      const T buf(abs(A(j, i)));
      if(max < buf)
        max = buf;
      if((buf < min || min <= T(0)) && buf != T(0))
        min = buf;
    }
    const T ratio(sqrt(min * max));
    if(ratio > T(0)) {
      acr[i]   *= ratio;
#if defined(WITHOUT_EIGEN)
      A.setCol(i, A.col(i) / ratio);
#else
      A.col(i) /= ratio;
#endif
    }
  }
  {
    T min(abs(b[0]));
    T max(min);
    for(int j = 1; j < A.rows(); j ++) {
      const T buf(abs(b[j]));
      if(max < buf)
        max = buf;
      if((buf < min || min <= T(0)) && buf != T(0))
        min = buf;
    }
    const T ratio(sqrt(min * max));
    if(ratio > T(0)) {
      bcr *= ratio;
      b   /= ratio;
    }
  }
  return;
}

template <typename T> T LP<T>::errorCheck(const Mat& A, const Vec& b) const
{
  // accuracy check.
  T max(1), min(1);
  for(int i = 0; i < A.rows(); i ++) {
    for(int j = 0; j < A.cols(); j ++) {
      const T buf(abs(A(i, j)));
      if(buf > T(0)) {
        if(buf <= min)
          min = buf;
        if(max <= buf)
          max = buf;
      }
    }
    const T buf(abs(b[i]));
    if(buf > T(0)) {
      if(buf <= min)
        min = buf;
      if(max <= buf)
        max = buf;
    }
  }
  return (min / max) / sqrt(err_error);
}

template <typename T> bool LP<T>::gainVectors(bool* fix, bool* checked, Vec& rvec, const Mat& P, const Vec& b, const Vec& one, int& n_fixed) const
{
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
  for(int i = 0; i < P.rows(); i ++) {
    fix[i]     = false;
    checked[i] = false;
  }
  n_fixed = 0;
  
  // set value, orthogonalize, and scale t.
  Vec bb(b - P * (P.transpose() * b));
  if(sqrt(bb.dot(bb)) <= threshold_feas * sqrt(b.dot(b))) {
    for(int i = 0; i < bb.size() - P.cols() * 2 - 1; i ++)
      bb[i] = - T(1);
    for(int i = bb.size() - P.cols() * 2 - 1; i < bb.size(); i ++)
      bb[i] =   b[i];
    const Vec bbb(bb - P * (P.transpose() * bb));
    if(sqrt(bbb.dot(bbb)) <= threshold_feas * sqrt(bb.dot(bb))) {
      rvec  = P.transpose() * (bbb / err_error + b + bb);
      cerr << " Second trivial matrix.";
      return true;
    }
    bb = bbb;
  }
  
  const T normbb(sqrt(bb.dot(bb)));
  bb /= normbb;
  
  rvec = giantStep(fix, checked, P, - bb, n_fixed, one);
  if(n_fixed == P.cols()) {
    std::cerr << " F";
    Mat F(P.cols(), P.cols());
    Vec f(P.cols());
    for(int i = 0, j = 0; i < P.rows() && j < f.size(); i ++)
      if(fix[i]) {
        const T ratio(sqrt(P.row(i).dot(P.row(i))));
        F.row(j) = P.row(i) / ratio;
        f[j]     = bb[i]    / ratio;
        j ++;
      }
#if defined(WITHOUT_EIGEN)
    rvec = F.solve(f);
#else
    rvec = F.inverse() * f;
#endif
  } else
    rvec = P.transpose() * rvec / rvec.dot(- bb); 
  rvec = rvec * normbb + P.transpose() * b;
  return isErrorMargin(P, b, rvec, false);
}

#if defined(WITHOUT_EIGEN)
template <typename T> SimpleVector<T> LP<T>::giantStep(bool* fix, bool* checked, Mat Pverb, Vec mb, int& n_fixed, const Vec& one) const
#else
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> LP<T>::giantStep(bool* fix, bool* checked, Mat Pverb, Vec mb, int& n_fixed, const Vec& one) const
#endif
{
  Vec norm(one.size());
  Vec deltab(one.size());
  Vec on(one.size());
  T   ratiob(1);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
  for(int i = 0; i < deltab.size(); i ++) {
    norm[i]   = T(1);
    deltab[i] = T(0);
    on[i]     = T(0);
  }
  for( ; n_fixed < Pverb.cols(); n_fixed ++) {
    // XXX : we should re implement around mb with accuracy.
    mb *= ratiob;
    mb += deltab;
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
    for(int j = 0; j < Pverb.rows(); j ++) {
      norm[j]     = sqrt(Pverb.row(j).dot(Pverb.row(j)));
      checked[j] |= norm[j] <= threshold_p0;
    }
    if(abs(norm.dot(norm) - T(1)) <= threshold_p0)
      cerr << "?" << std::endl;
    norm   /= sqrt(norm.dot(norm));
    deltab  = Pverb * (Pverb.transpose() * mb);
    mb     -= deltab;
    const T   buf(sqrt(mb.dot(mb)));
    ratiob  = buf;
    mb     /= buf;
    
    // O(mn^2) check for inner or not.
    const Vec work(mb + (norm - Pverb * (Pverb.transpose() * norm)) * threshold_loop);
    on = Pverb * (Pverb.transpose() * (- one)) + work * work.dot(- one) / work.dot(work);
    const int max_idx = getMax(checked, on, norm);
    if(max_idx >= one.size()) {
      cerr << " GiantStep: no more direction.";
      break;
    } else if(checkInner(on, norm, max_idx))
      break;
    
    // O(mn^2) over all in this function.
    // this makes some error after some steps.
#if defined(WITHOUT_EIGEN)
    const Vec orth(Pverb.row(max_idx));
#else
    const Vec orth(Pverb.row(max_idx).transpose());
#endif
    const T   norm2orth(orth.dot(orth));
    const T   mb0(mb[max_idx]);
#if defined(_OPENMP)
#pragma omp for
#endif
    for(int j = 0; j < Pverb.rows(); j ++) {
      const T work(Pverb.row(j).dot(orth) / norm2orth);
#if defined(WITHOUT_EIGEN)
      Pverb.row(j) -= orth * work;
#else
      Pverb.row(j) -= work * orth.transpose();
#endif
      mb[j]        -= mb0  * work;
    }
    fix[max_idx] = 1;
  }
  return on * ratiob + deltab;
}

template <typename T> bool LP<T>::checkInner(const Vec& on, const Vec& normalize, const int& max_idx) const
{
  return max_idx < on.size() && on[max_idx] / normalize[max_idx] <= threshold_inner;
}

template <typename T> int LP<T>::getMax(const bool* checked, const Vec& on, const Vec& normalize) const {
  int result = 0;
  for(; result < on.size(); result ++)
    if(!checked[result])
      break;
  for(int j = result + 1; j < on.size(); j ++)
    if(!checked[j] && on[result] / normalize[result] < on[j] / normalize[j])
      result = j;
  if(checked[result])
    result = on.size();
  return result;
}

template <typename T> bool LP<T>::isErrorMargin(const Mat& A, const Vec& b, const Vec& x, const bool& disp) const
{
  T result(0);
  const Vec err(A * x - b);
  for(int i = 0; i < b.size(); i ++)
    if(result < err[i]) result = err[i];
  if(disp)
    cerr << " errorMargin?(" << sqrt(x.dot(x)) << ", " << sqrt(b.dot(b)) << ", " << result << ")";
  return isfinite(x.dot(x)) && (x.dot(x) > err_norm(err_error) ? result / sqrt(x.dot(x)) : result) <= err_error;
  // XXX: or use this??
  // return result <= err_error;
}

#if defined(WITHOUT_EIGEN)
template <typename T> SimpleMatrix<T> LP<T>::roughQR(const Mat& A) const {
#else
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> LP<T>::roughQR(const Mat& A) const {
#endif
  Mat Q(A.cols(), A.rows());
  for(int i = 0; i < Q.rows(); i ++)
    for(int j = 0; j < Q.cols(); j ++)
      Q(i, j) = T(0);
  Mat work(A.transpose());
  for(int i = 0; i < A.cols(); i ++) {
#if defined(WITHOUT_EIGEN)
    work.row(i) -= Q.transpose() * (Q * work.row(i));
    Q.row(i) = work.row(i) / sqrt(work.row(i).dot(work.row(i)));
#else
    work.row(i) -= (Q.transpose() * (Q * work.row(i).transpose())).transpose();
    Q.row(i) = work.row(i) / sqrt(work.row(i).dot(work.row(i)));
#endif
  }
  return Q.transpose();
}

#define _LINEAR_OPT_
#endif

