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
#include <limits>
#include <assert.h>

#if !defined(WITHOUT_EIGEN)
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::fflush;

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
        void resize(const int& size);
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
  SimpleVector<T> work(other.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < esize; i ++)
    work[i] = entity[i] * other.entity[i];
  for(int i = 0; i < esize; i ++)
    res += work[i];
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

template <typename T> void SimpleVector<T>::resize(const int& size) {
  assert(size > 0);
  if(size != esize) {
    esize = size;
    if(entity)
      delete[] entity;
    entity = new T[esize];
  }
  return;
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
  const SimpleVector<T>& row(const int& y) const;
  const SimpleVector<T>  col(const int& x) const;
        void             setCol(const int& x, const SimpleVector<T>& other);
        SimpleMatrix<T>  transpose() const;
        SimpleVector<T>  solve(SimpleVector<T> other) const;
        SimpleVector<T>  projectionPt(const SimpleVector<T>& other) const;
  const int& rows() const;
  const int& cols() const;
        void resize(const int& rows, const int& cols);
private:
  SimpleVector<T>* entity;
  int              erows;
  int              ecols;
};

template <typename T> SimpleMatrix<T>::SimpleMatrix(const int& rows, const int& cols) {
  assert(rows > 0 && cols > 0);
  entity = new SimpleVector<T>[rows];
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < rows; i ++)
    entity[i].resize(cols);
  erows = rows;
  ecols = cols;
  return; 
}

template <typename T> SimpleMatrix<T>::SimpleMatrix(const SimpleMatrix<T>& other) {
  erows = other.erows;
  ecols = other.ecols;
  entity = new SimpleVector<T>[other.erows];
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
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
#pragma omp parallel for schedule(static, 1)
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
#pragma omp parallel for schedule(static, 1)
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
#pragma omp parallel for schedule(static, 1)
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
#pragma omp parallel for schedule(static, 1)
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
#pragma omp parallel for schedule(static, 1)
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
#pragma omp parallel for schedule(static, 1)
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
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < other.erows; i ++)
      entity[i].resize(other.ecols);
  }
  erows = other.erows;
  ecols = other.ecols;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
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

template <typename T> const SimpleVector<T>& SimpleMatrix<T>::row(const int& y) const {
  assert(0 <= y && y < erows && entity);
  return entity[y];
}

template <typename T> const SimpleVector<T> SimpleMatrix<T>::col(const int& x) const {
  assert(0 <= erows && 0 <= x && x < ecols && entity);
  SimpleVector<T> res(erows);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    res[i] = entity[i][x];
  return res;
}

template <typename T> void SimpleMatrix<T>::setCol(const int& x, const SimpleVector<T>& other) {
  assert(0 <= x && x < ecols && other.size() == erows);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < erows; i ++)
    entity[i][x] = other[i];
  return;
}

template <typename T> SimpleMatrix<T> SimpleMatrix<T>::transpose() const {
  SimpleMatrix<T> res(ecols, erows);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
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
    int xchg = i;
    for(int j = i + 1; j < erows; j ++)
      if(abs(work.entity[j][i]) > abs(work.entity[xchg][i]))
        xchg = j;
    SimpleVector<T> buf(work.entity[i]);
    T               buf2(other[i]);
    work.entity[i]    = work.entity[xchg];
    other[i]          = other[xchg];
    work.entity[xchg] = buf;
    other[xchg]       = buf2;
    const SimpleVector<T>& ei(work.entity[i]);
    const T&               eii(ei[i]);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = i + 1; j < erows; j ++) {
      const T ratio(work.entity[j][i] / eii);
      work.entity[j] -= ei       * ratio;
      other[j]       -= other[i] * ratio;
    }
  }
  for(int i = erows - 1; 0 <= i; i --) {
    const T buf(other[i] / work.entity[i][i]);
    if(!isfinite(buf) || isnan(buf)) {
      assert(!isfinite(work.entity[i][i] / other[i]) || isnan(work.entity[i][i] / other[i]));
      continue;
    }
    other[i]    = buf;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = i - 1; 0 <= j; j --)
      other[j] -= other[i] * work.entity[j][i];
  }
  return other;
}

template <typename T> SimpleVector<T> SimpleMatrix<T>::projectionPt(const SimpleVector<T>& other) const {
  assert(0 < erows && 0 < ecols && ecols == other.size());
  // also needs class or this->transpose() * (*this) == I assertion is needed.
  SimpleMatrix<T> work(erows, ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < work.rows(); i ++)
    work.row(i) = entity[i] * entity[i].dot(other);
  SimpleVector<T> res(ecols);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < other.size(); i ++) {
    res[i] = T(0);
    for(int j = 0; j < erows; j ++)
      res[i] += work(j, i);
  }
  return res;
}

template <typename T> const int& SimpleMatrix<T>::rows() const {
  return erows;
}

template <typename T> const int& SimpleMatrix<T>::cols() const {
  return ecols;
}

template <typename T> void SimpleMatrix<T>::resize(const int& rows, const int& cols) {
  assert(rows > 0 && cols > 0);
  if(rows != erows) {
    erows = rows;
    if(entity)
      delete[] entity;
    entity = new SimpleVector<T>[erows];
  }
  if(cols != ecols) {
    ecols = cols;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < erows; i ++)
      entity[i].resize(ecols);
  }
  return;
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
  bool optimize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const;
  bool optimizeNullSpace(bool* fix_partial, Vec& rvec, const Mat& P, Vec& b, const int& cidx) const;
  
  T    errorCheck(const Mat& A, const Vec& b) const;
  
  bool gainVectors(bool* fix, char* checked, Vec& rvec, const Mat& Pt, Vec b, const Vec& one, int& n_fixed) const;
  Vec  giantStep(bool* fix, char* checked, Mat Pverb, Vec mbb, int& n_fixed, const Vec& one) const;
  bool checkInner(const Vec& on, const Vec& normalize, const int& max_idx) const;
  int  getMax(const char* checked, const Vec& on, const Vec& normalize) const;
  
  Mat  roughQR(const Mat& A) const;

  T threshold_feas;
  T threshold_p0;
  T threshold_loop;
  T threshold_inner;
  T err_error;
  T mr_intercept;
  T largest_intercept;
  T largest_opt;
  int n_opt_steps;
};

template <typename T> LP<T>::LP()
{
  // error rate for orthogonalized b.
  threshold_feas    = pow(numeric_limits<T>::epsilon(), T(3) / T(4));
  
  // if SimpleVector<T>::solve fails, increase this:
  threshold_p0      = sqrt(threshold_feas);
  
  // is we gain inner or not:
  threshold_inner   = sqrt(threshold_p0);
  
  // last stage error rate:
  err_error         = sqrt(threshold_inner);

  // box constraints that Px<=b, not b<=Px:
  largest_intercept = T(1) / sqrt(err_error);
  
  // largest optimized value ratio to be calculated:
  largest_opt       = sqrt(largest_intercept);
  assert(sqrt(largest_opt) > T(1));
  
  // XXX: Please configure me first.
  // if threshold_loop < 0, extends some regions.
  threshold_loop    = - err_error;
  mr_intercept      = T(1e12);
  
  // something bugly with QD library.
  // n_opt_steps       = int(log(largest_opt) / log(T(2))) + 1;
  n_opt_steps       = 18;
  return;
}

template <typename T> LP<T>::~LP()
{
  return;
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

template <typename T> bool LP<T>::optimize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const
{
  assert(A.cols() == c.size() && A.rows() == b.size() && 0 < c.size() && 0 < b.size());
  // cout << A << endl << b << endl << c.transpose() << endl;
  cerr << " (" << A.rows() << ", " << A.cols() << ")";
  rvec = Vec(c.size());
  for(int i = 0; i < A.rows(); i ++)
    if(fix_partial) fix_partial[i] = false;
  
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
  cerr << " err_error(" << errorCheck(AA, bb) << ")";
  
  // guarantee that b is positive.
  T rbb(0);
  for(int i = 0; i < A.rows(); i ++) {
    const T lbb(abs(b[i]) / sqrt(A.row(i).dot(A.row(i))));
    if(isfinite(lbb) && !isnan(lbb))
      rbb = max(rbb, lbb);
  }
  const T intercept(min(largest_intercept, mr_intercept * rbb) / sqrt(T(2 * A.cols())));
  cerr << " intercept(" << rbb / intercept << ")";
  for(int i = 0; i < A.cols() * 2; i ++) {
    AA(A.rows() + 1 + i, i / 2) = (i % 2 == 0 ? T(1) : - T(1));
    bb[A.rows() + 1 + i]        = intercept;
  }
  if(c.dot(c) != T(0))
    bb[A.rows()] = min(largest_opt, intercept);
  
  for(int i = 0; i < AA.rows(); i ++) {
    assert(isfinite(bb[i]) && !isnan(bb[i]));
    for(int j = 0; j < AA.cols(); j ++)
      assert(isfinite(AA(i, j)) && !isnan(AA(i, j)));
  }
  cerr << " .";
  
  // <c, x> -> obj, Q R x <= b.
  // O(mn * max(m, n))
  Mat Q(roughQR(AA));
  cerr << "Q";
  fflush(stderr);
  Mat R(Q.transpose() * AA);
  cerr << "R";
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
  
  cerr << " NULL";
  fflush(stderr);
  
  // pull original solvee.
#if defined(WITHOUT_EIGEN)
  rvec = R.solve(rvec);
#else
  rvec = R.inverse() * rvec;
#endif

  cerr << " INV" << endl;
  return isErrorMargin(A, b, rvec, true);
}

template <typename T> bool LP<T>::optimizeNullSpace(bool* fix_partial, Vec& rvec, const Mat& P, Vec& b, const int& cidx) const
{
  assert(P.rows() == b.size() && 0 < P.cols() && P.cols() < P.rows());
  
  // to avoid stack underflow.
  char* checked = new char[b.size()];
  Vec   one(b.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < one.size(); i ++)
    one[i] = T(1);
  const Mat Pt(P.transpose());
  
  int  n_fixed;
  if(P.row(cidx).dot(P.row(cidx)) <= T(0)) {
    (void)gainVectors(fix_partial, checked, rvec, Pt, b, one, n_fixed);
    delete[] checked;
    return isErrorMargin(P, b, rvec, true);
  }
  
  // get optimal value.
  T offset(b[cidx]);
  T mid(1 / offset);
  T step(offset * offset * T(2));
  for(int i = 0; i < n_opt_steps; i ++) {
    b[cidx] = - step * mid + offset;
    if(!isfinite(b[cidx]) || isnan(b[cidx]))
      continue;
    if(gainVectors(fix_partial, checked, rvec, Pt, b, one, n_fixed)) {
      mid *= step;
      cerr << ".";
    }
    step = sqrt(step);
  }
  cerr << endl;
  
  b[cidx] = - mid + offset;
  (void)gainVectors(fix_partial, checked, rvec, Pt, b, one, n_fixed);
  delete[] checked;
  return isErrorMargin(P, b, rvec, true);
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
  return (max / min) * sqrt(err_error);
}

template <typename T> bool LP<T>::gainVectors(bool* fix, char* checked, Vec& rvec, const Mat& Pt, Vec b, const Vec& one, int& n_fixed) const
{
  Vec norm(Pt.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < Pt.cols(); i ++) {
    fix[i]     = false;
    checked[i] = false;
    norm[i]    = sqrt(Pt.col(i).dot(Pt.col(i)));
  }
  T normb0(0);
  for(int i = 0; i < Pt.cols() - Pt.rows() * 2 - 1; i ++)
    normb0 += b[i] * b[i];
  normb0 = sqrt(normb0);
  
  // extend intercepts with 1 epsilon.
  b -= norm * normb0 * threshold_loop;

  // set value, orthogonalize, and scale t.
#if defined(WITHOUT_EIGEN)
  Vec bb(b - Pt.projectionPt(b));
#else
  Vec bb(b - Pt.transpose() * (Pt * b));
#endif
  if(sqrt(bb.dot(bb)) <= threshold_feas * sqrt(b.dot(b))) {
    cerr << "0";
    fflush(stderr);
    for(int i = 0; i < bb.size() - Pt.rows() * 2 - 1; i ++)
      bb[i] = sqrt(Pt.col(i).dot(Pt.col(i)));
    for(int i = bb.size() - Pt.rows() * 2 - 1; i < bb.size(); i ++)
      bb[i] = b[i];
#if defined(WITHOUT_EIGEN)
    const Vec bbb(bb - Pt.projectionPt(bb));
#else
    const Vec bbb(bb - Pt.transpose() * (Pt * bb));
#endif
    if(sqrt(bbb.dot(bbb)) <= threshold_feas * sqrt(bb.dot(bb))) {
      rvec  = Pt * (bbb / err_error + b + bb);
      cerr << " Second trivial matrix.";
      return true;
    }
    bb  = bbb;
  }
  
  rvec = giantStep(fix, checked, Pt, - bb, n_fixed, one);
  fflush(stderr);
  if(n_fixed == Pt.rows()) {
    cerr << "F";
    Mat F(Pt.rows(), Pt.rows());
    Vec f(Pt.rows());
    for(int i = 0, j = 0; i < Pt.cols() && j < f.size(); i ++)
      if(fix[i]) {
        const T ratio(sqrt(Pt.col(i).dot(Pt.col(i)) + b[i] * b[i]));
        F.row(j) = Pt.col(i) / ratio;
        f[j]     = b[i]      / ratio;
        j ++;
      }
#if defined(WITHOUT_EIGEN)
    rvec = F.solve(f);
#else
    rvec = F.inverse() * f;
#endif
  } else {
    cerr << "g";
    rvec = Pt * (rvec + b);
  }
  fflush(stderr);
  for(int i = 0; i < rvec.size(); i ++)
    if(!isfinite(rvec[i]) || isnan(rvec[i]))
      return false;
  // XXX: Pt.transpose can be cached.
  return isErrorMargin(Pt.transpose(), b, rvec, false);
}

#if defined(WITHOUT_EIGEN)
template <typename T> SimpleVector<T> LP<T>::giantStep(bool* fix, char* checked, Mat Pverb, Vec mbb, int& n_fixed, const Vec& one) const
#else
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> LP<T>::giantStep(bool* fix, char* checked, Mat Pverb, Vec mbb, int& n_fixed, const Vec& one) const
#endif
{
  Vec on(one.size());
  Vec norm(one.size());
  Vec deltab(one.size());
  T   ratiob(1);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < norm.size(); i ++) {
    norm[i]   = T(1);
    deltab[i] = T(0);
    on[i]     = T(0);
  }
  Vec mb(mbb.size());
  Vec orthbuf(mbb.size());
  for(n_fixed = 0 ; n_fixed < Pverb.rows(); n_fixed ++) {
    mb = mbb;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pverb.cols(); j ++) {
      norm[j]     = sqrt(Pverb.col(j).dot(Pverb.col(j)));
      checked[j]  = (fix[j] || norm[j] <= threshold_p0) ? 1 : 0;
    }
#if defined(WITHOUT_EIGEN)
    deltab = Pverb.projectionPt(mb);
#else
    deltab = Pverb.transpose() * (Pverb * mb);
#endif
    mb    -= deltab;
    ratiob = sqrt(mb.dot(mb));
    mb    /= ratiob;
    norm  /= sqrt(norm.dot(norm));
    
    // O(mn^2) check for inner or not.
#if defined(WITHOUT_EIGEN)
    on = Pverb.projectionPt(- one) + mb * mb.dot(- one) / mb.dot(mb);
#else
    on = Pverb.transpose() * (Pverb * (- one)) + mb * mb.dot(- one) / mb.dot(mb);
#endif
    on /= sqrt(on.dot(on));
    const int fidx = getMax(checked, on, norm);
    if(fidx >= one.size()) {
      cerr << " GiantStep: no more direction.";
      break;
    } else if(checkInner(on, norm, fidx)) {
      n_fixed --;
      on /= abs(mb.dot(on));
      break;
    }
    
    // O(mn^2) over all in this function.
    orthbuf = Pverb.col(fidx);
    const Vec& orth(orthbuf);
    const T    norm2orth(orth.dot(orth));
    const T    mbb0(mbb[fidx]);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pverb.cols(); j ++) {
      const T work(Pverb.col(j).dot(orth) / norm2orth);
#if defined(WITHOUT_EIGEN)
      Pverb.setCol(j, Pverb.col(j) - orth * work);
#else
      Pverb.col(j) -= orth * work;
#endif
      mbb[j]       -= mbb0 * work;
    }
    fix[fidx] = true;
  }
  return on * ratiob + deltab;
}

template <typename T> bool LP<T>::checkInner(const Vec& on, const Vec& normalize, const int& max_idx) const
{
  assert(0 <= max_idx && max_idx < on.size());
  return on[max_idx] / normalize[max_idx] <= threshold_inner;
}

template <typename T> int LP<T>::getMax(const char* checked, const Vec& on, const Vec& normalize) const {
  int result = 0;
  for(; result < on.size(); result ++)
    if(!checked[result])
      break;
  for(int j = result + 1; j < on.size(); j ++)
    if(!checked[j] && on[result] / normalize[result] < on[j] / normalize[j])
      result = j;
  if(result < on.size() && checked[result])
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
  return isfinite(result)         && !isnan(result) &&
         isfinite(x.dot(x))       && !isnan(x.dot(x)) &&
         isfinite(sqrt(x.dot(x))) && !isnan(sqrt(x.dot(x))) &&
         (sqrt(x.dot(x)) > err_error ? result / sqrt(x.dot(x)) : result) <=
           err_error &&
         sqrt(x.dot(x)) < sqrt(b.dot(b)) * largest_intercept;
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
  // N.B. .transpose() costs a lot.
  Mat work(A.transpose());
  for(int i = 0; i < A.cols(); i ++) {
#if defined(WITHOUT_EIGEN)
    work.row(i) -= Q.projectionPt(work.row(i));
#else
    work.row(i) -= (Q.transpose() * (Q * work.row(i).transpose())).transpose();
#endif
    // generally, assert norm > error is needed.
    // in this case, not.
    Q.row(i) = work.row(i) / sqrt(work.row(i).dot(work.row(i)));
  }
  return Q.transpose();
}

#define _LINEAR_OPT_
#endif

