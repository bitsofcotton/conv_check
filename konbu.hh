/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/* BSD 3-Clause License:
 * Copyright (c) 2013 - 2019, kazunobu watatsu.
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
using std::flush;
using std::move;

template <typename T> class SimpleVector {
public:
  SimpleVector();
  SimpleVector(const int& size);
  SimpleVector(const SimpleVector<T>& other);
  SimpleVector(SimpleVector<T>&& other);
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
        SimpleVector<T>& operator =  (SimpleVector<T>&& other);
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

template <typename T> SimpleVector<T>::SimpleVector(SimpleVector<T>&& other) {
  esize  = move(other.esize);
  entity = move(other.entity);
  other.entity = 0;
  return;
}

template <typename T> SimpleVector<T>::~SimpleVector() {
  delete[] entity;
  entity = NULL;
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

template <typename T> SimpleVector<T>& SimpleVector<T>::operator = (SimpleVector<T>&& other) {
  if(entity == other.entity && esize == other.esize)
    return *this;
  esize  = move(other.esize);
  delete[] entity;
  entity = move(other.entity);
  other.esize  = 0;
  other.entity = NULL;
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
  SimpleMatrix(SimpleMatrix<T>&& other);
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
        SimpleMatrix<T>& operator =  (SimpleMatrix<T>&& other);
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

template <typename T> SimpleMatrix<T>::SimpleMatrix(SimpleMatrix<T>&& other) {
  erows  = move(other.erows);
  ecols  = move(other.ecols);
  entity = move(other.entity);
  other.entity = NULL;
  return;
}

template <typename T> SimpleMatrix<T>::~SimpleMatrix() {
  delete[] entity;
  entity = NULL;
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
  assert(erows == other.erows && ecols == other.ecols);
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

template <typename T> SimpleMatrix<T>& SimpleMatrix<T>::operator = (SimpleMatrix<T>&& other) {
  if(entity == other.entity && erows == other.erows && ecols == other.ecols)
    return *this;
  erows  = move(other.erows);
  ecols  = move(other.ecols);
  delete[] entity;
  entity = move(other.entity);
  other.erows  = 0;
  other.ecols  = 0;
  other.entity = NULL;
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
      throw "Non full rank matrix SimpleMatrix::solve";
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


template <typename T> class Linner {
public:
#if defined(WITHOUT_EIGEN)
  typedef SimpleMatrix<T> Mat;
  typedef SimpleVector<T> Vec;
#else
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
#endif
  
  Linner();
  ~Linner();
  
  bool inner(bool* fix_partial, const Mat& A, const Vec& b) const;
  bool inner(Vec& rvec, const Mat& A, const Vec& b) const;
  bool inner(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b) const;
  bool isErrorMargin(const Mat& A, const Vec& b, const Vec& x, const bool& disp = true) const;
  
private:
  bool gainVectors(bool* fix, bool* checked, Vec& rvec, const Mat& Pt, const Vec& b) const;
  Mat  roughQR(const Mat& At) const;

  T threshold_feas;
  T threshold_p0;
  T threshold_loop;
  T threshold_inner;
  T err_error;
  T large;
};

template <typename T> Linner<T>::Linner() {
  // error rate for orthogonalized b.
  threshold_feas    = numeric_limits<T>::epsilon() * T(8);
  //threshold_feas    = T(1) >> short(62);
  
  // if SimpleVector<T>::solve fails, increase this:
  threshold_p0      = pow(threshold_feas, T(2) / T(3));
  
  // is we gain inner or not:
  threshold_inner   = pow(threshold_feas, T(1) / T(3));
  
  // last stage error rate:
  err_error         = sqrt(threshold_inner);
  
  // guarantee b is positive.
  large             = T(1) / sqrt(err_error);
  
  assert(sqrt(err_error) < T(1) && T(1) < large);
  
  // XXX: Please configure me first.
  // if threshold_loop < 0, extends some regions.
  threshold_loop    = T(0);
  
  return;
}

template <typename T> Linner<T>::~Linner() {
  return;
}

template <typename T> bool Linner<T>::inner(bool* fix_partial, const Mat& A, const Vec& b) const {
  Vec rvec;
  return inner(fix_partial, rvec, A, b);
}

template <typename T> bool Linner<T>::inner(Vec& rvec, const Mat& A, const Vec& b) const {
  return inner(0, rvec, A, b);
}

template <typename T> bool Linner<T>::inner(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b) const {
  assert(A.rows() == b.size() && 0 < A.cols() && 0 < b.size());
  // cout << A << endl << b << endl << c.transpose() << endl;
  cerr << " (" << A.rows() << ", " << A.cols() << ")";
  Mat AA(A.rows() + A.cols() * 2, A.cols());
  Vec bb(AA.rows());
  for(int i = 0; i < A.rows(); i ++) {
    AA.row(i) = A.row(i);
    bb[i]     = b[i];
  }
  for(int i = 0; i < A.cols() * 2; i ++) {
    for(int j = 0; j < A.cols(); j ++)
      AA(A.rows() + i, j) = i / 2 == j ? (i % 2 ? - T(1) : T(1)) : T(0);
    bb[A.rows()] = T(0);
  }
  
  T M(1), m(1);
  for(int i = 0; i < A.rows(); i ++) {
    for(int j = 0; j < A.cols(); j ++) {
      const T buf(abs(A(i, j)));
      if(buf != T(0)) {
        M = max(M, buf);
        m = min(m, buf);
      }
      assert(isfinite(buf) && !isnan(buf));
    }
    const T buf(abs(b[i]));
    if(buf != T(0)) {
      M = max(M, buf);
      m = min(m, buf);
    }
    assert(isfinite(buf) && !isnan(buf));
  }
  cerr << " err_error(" << (M / m) * sqrt(err_error) << ")" << flush;
  for(int i = A.rows(); i < AA.rows(); i ++)
    bb[i] = max(T(1), sqrt(b.dot(b))) * large;
  
  bool* bfix_partial = new bool[AA.rows()];
  bool* checked = new bool[AA.rows()];
  const Mat Pt(roughQR(AA.transpose()));
  cerr << "Q" << flush;
  const Mat R(Pt * AA);
  cerr << "R" << flush;
  rvec = Vec(A.cols());
  const bool f_inner(gainVectors(bfix_partial, checked, rvec, Pt, bb));
  if(fix_partial)
    for(int i = 0; i < A.rows(); i ++)
      fix_partial[i] = bfix_partial[i];
  delete[] bfix_partial;
  delete[] checked;
#if defined(WITHOUT_EIGEN)
  rvec = R.solve(rvec);
#else
  rvec = R.inverse() * rvec;
#endif
  cerr << "I" << flush;
  return f_inner || isErrorMargin(A, b, rvec, true);
}

template <typename T> bool Linner<T>::gainVectors(bool* fix, bool* checked, Vec& rvec, const Mat& Pt, const Vec& b) const {
  rvec = Vec(Pt.rows());
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < rvec.size(); i ++)
    rvec[i] = T(0);
  int n_fixed;
  Vec one(b.size());
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < Pt.cols(); i ++) {
    one[i]     = T(1);
    fix[i]     = false;
    checked[i] = false;
  }
  
  // set value, orthogonalize, and scale t.
#if defined(WITHOUT_EIGEN)
  Vec bb(b - Pt.projectionPt(b));
#else
  Vec bb(b - Pt.transpose() * (Pt * b));
#endif
  if(sqrt(bb.dot(bb)) <= threshold_feas * sqrt(b.dot(b))) {
    for(int i = 0; i < bb.size() - 1; i ++)
      bb[i] = sqrt(Pt.col(i).dot(Pt.col(i)));
    for(int i = bb.size() - 1; i < bb.size(); i ++)
      bb[i] = b[i];
#if defined(WITHOUT_EIGEN)
    const Vec bbb(bb - Pt.projectionPt(bb));
#else
    const Vec bbb(bb - Pt.transpose() * (Pt * bb));
#endif
    if(sqrt(bbb.dot(bbb)) <= threshold_feas * sqrt(bb.dot(bb))) {
      rvec  = Pt * (b - bb - bbb * large);
      cerr << "t(" << rvec[0] << ")" << flush;
      return isErrorMargin(Pt.transpose(), b, rvec, ! false);
    }
    bb  = bbb;
    std::cerr << "0" << flush;
  }
  Vec mbb(- bb);
  const T normb0(sqrt(mbb.dot(mbb)));
  Mat Pverb(Pt);
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
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pverb.cols(); j ++) {
      norm[j]    = sqrt(Pverb.col(j).dot(Pverb.col(j)));
      checked[j] = fix[j] || norm[j] <= threshold_p0;
    }
    // extend b with threshold_loop. N.B. mbb = - b'.
    mb     = mbb + norm * normb0 * threshold_loop;
#if defined(WITHOUT_EIGEN)
    deltab = Pverb.projectionPt(mb);
#else
    deltab = Pverb.transpose() * (Pverb * mb);
#endif
    mb    -= deltab;
    ratiob = sqrt(mb.dot(mb));
    mb    /= ratiob;
    
    // O(mn^2) check for inner or not.
#if defined(WITHOUT_EIGEN)
    on = Pverb.projectionPt(- one) + mb * mb.dot(- one);
#else
    on = Pverb.transpose() * (Pverb * (- one)) + mb * mb.dot(- one);
#endif
  
    // N.B. This might be rewrited as fix once with sorted fidxs.
    //   When the hypothesis is true, next loop must fix or get inner vector.
    //   Also, it takes O(lg(mn)) time order for this function.
    //   But to do so, theoretical proof is needed and it seems not
    //   because we use mb with little extended and its norm is in error order.
    int fidx(0);
    for( ; fidx < on.size(); fidx ++)
      if(!checked[fidx])
        break;
    for(int j = fidx + 1; j < on.size(); j ++)
      if(!checked[j] && on[fidx] / norm[fidx] < on[j] / norm[j])
        fidx = j;
    if(fidx >= one.size())
      assert(0 && "rank is not full : should not be reached");
    else if(on[fidx] / norm[fidx] <= threshold_inner) {
      n_fixed --;
      on /= norm[fidx];
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
  if(n_fixed == Pt.rows()) {
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
    cerr << "F" << flush;
  } else {
    rvec = Pt * (on * ratiob + deltab + b);
    cerr << "g" << flush;
  }
  return isErrorMargin(Pt.transpose(), b, rvec, ! false);
}

template <typename T> bool Linner<T>::isErrorMargin(const Mat& A, const Vec& b, const Vec& x, const bool& disp) const {
  T result(0);
  const Vec err(A * x - b);
  for(int i = 0; i < b.size(); i ++) {
    //const auto lerr(err[i] / sqrt(A.row(i).dot(A.row(i))));
    const auto& lerr(err[i]);
    if(!isfinite(lerr) || isnan(lerr) ||
       result < lerr) result = lerr;
  }
  if(disp)
    cerr << " errorMargin?(" << result / max(T(1), sqrt(x.dot(x))) << ")";
  return isfinite(result)   && !isnan(result)   &&
         isfinite(x.dot(x)) && !isnan(x.dot(x)) &&
         ((x.dot(x) == T(0) && result == T(0)) ||
          (x.dot(x) != T(0) && result / sqrt(x.dot(x)) <= err_error));
}

template <typename T> typename Linner<T>::Mat Linner<T>::roughQR(const Mat& At) const {
  Mat Q(At.rows(), At.cols());
  for(int i = 0; i < Q.rows(); i ++)
    for(int j = 0; j < Q.cols(); j ++)
      Q(i, j) = T(0);
  for(int i = 0; i < At.rows(); i ++) {
    // N.B. in this case, At is full rank.
#if defined(WITHOUT_EIGEN)
    const Vec work(At.row(i) - Q.projectionPt(At.row(i)));
#else
    const Vec work(At.row(i) - (Q.transpose() * (Q * At.row(i).transpose())).transpose());
#endif
    // generally, assert norm > error is needed.
    // in this case, not.
    Q.row(i) = work / sqrt(work.dot(work));
  }
  return Q;
}

#define _LINEAR_OPT_
#endif

