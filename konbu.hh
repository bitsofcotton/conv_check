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
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

template <typename T, typename S> class LP {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> MatS;
  typedef Eigen::Matrix<S, Eigen::Dynamic, 1> VecS;
  
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
  
  T    steps(bool* fix, bool* checked, Vec& rvec, const Mat& P, const Vec& b, const Vec& bdash, T normbb, const VecS& one) const;
  bool gainVectors(bool& proper, bool* fix, bool* checked, Vec& rvec, const Mat& P, Vec b, const T& orignormbb, const VecS& one, const bool& only_check, int& n_fixed) const;
  bool giantStep(bool* fix, bool* checked, MatS& Pverb, int& n_fixed, const VecS& one) const;
  
  bool checkInner(const VecS& on, const int& max_idx) const;
  int  getMax(const bool* checked, const Vec& on) const;

  bool svdschur(const Mat& A, Mat& U, Vec& w, Mat& V) const;
  
  T err_raw_epsilon;
  T err_singular;
  T err_feas_step;
  T err_fstep_max;
  T err_opt;
  S threshold_feas;
  S threshold_p0;
  S threshold_inner;
  T err_error;
};

template <typename T, typename S> LP<T,S>::LP()
{
  err_raw_epsilon = NumTraits<T>::epsilon();
  err_singular    = internal::sqrt(err_raw_epsilon);
  threshold_feas  = internal::max(S(internal::sqrt(err_singular)), NumTraits<S>::epsilon());
  threshold_p0    = threshold_feas;
  err_feas_step   = internal::pow(err_raw_epsilon, T(1) / T(4));
  err_fstep_max   = T(1000);
  err_opt         = internal::sqrt(err_feas_step);
  threshold_inner = T(0);
  err_error       = internal::sqrt(threshold_feas);
  if(err_error >= T(1) - err_raw_epsilon)
    cerr << " accuracy not enough in initializing : " << err_error << endl;
  return;
}

template <typename T, typename S> LP<T,S>::~LP()
{
  return;
}

template <typename T, typename S> const T LP<T,S>::err_norm(const T& err) const
{
  // margin needed for stable calculation.
  return err;
}

template <typename T, typename S> bool LP<T,S>::minimize(bool* fix_partial, const Mat& A, const Vec& b, const Vec& c) const
{
  Vec rvec;
  return minimize(fix_partial, rvec, A, b, c);
}

template <typename T, typename S> bool LP<T,S>::minimize(Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const
{
  return minimize(0, rvec, A, b, c);
}

template <typename T, typename S> bool LP<T,S>::minimize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const
{
  return optimize(fix_partial, rvec, A, b, c);
}

template <typename T, typename S> bool LP<T,S>::maximize(bool* fix_partial, const Mat& A, const Vec& b, const Vec& c) const
{
  Vec rvec;
  return maximize(fix_partial, rvec, A, b, c);
}

template <typename T, typename S> bool LP<T,S>::maximize(Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const
{
  return maximize(0, rvec, A, b, c);
}

template <typename T, typename S> bool LP<T,S>::maximize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, const Vec& c) const
{
  return optimize(fix_partial, rvec, A, b, - c);
}

template <typename T, typename S> bool LP<T,S>::inner(bool* fix_partial, const Mat& A, const Vec& b) const
{
  Vec rvec;
  return inner(fix_partial, rvec, A, b);
}

template <typename T, typename S> bool LP<T,S>::inner(Vec& rvec, const Mat& A, const Vec& b) const
{
  return inner(0, rvec, A, b);
}

template <typename T, typename S> bool LP<T,S>::inner(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b) const
{
  Vec c(A.cols());
  for(int i = 0; i < c.size(); i ++)
    c[i] = T(0);
  return optimize(fix_partial, rvec, A, b, c);
}

template <typename T, typename S> bool LP<T,S>::optimize(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b, Vec c) const
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
  AA.row(A.rows()) = c.transpose();
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
  
  // guarantee that b is positive.
  for(int i = 0; i < A.cols() * 2; i ++)
    AA(A.rows() + 1 + i, i / 2) = (i % 2 == 0 ? T(1) : - T(1));
  
  // check accuracy.
  const T errr(errorCheck(AA, bb));
  if(errr < T(1))
    cerr << " err_error?(" << errr << ")";
  
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
  
  // <c, x> -> obj, U w V^t x <= b.
  // O(mn * max(m, n))
  Mat U;
  Vec w;
  Mat V;
  if(!svdschur(AA, U, w, V))
    return false;;
  cerr << " SCHUR";
  fflush(stderr);
  
  // optimize <co_c, x'>, co_A [x' 0] <= b. co_A^t co_A = I, in O(mn^2L).
  bool* bfix_partial = new bool[U.rows()];
  for(int i = 0; i < U.rows(); i ++)
    bfix_partial[i] = false;
  (void)optimizeNullSpace(bfix_partial, rvec, U, bb, A.rows());
  if(fix_partial)
    for(int i = 0; i < A.rows(); i ++)
      fix_partial[i] = bfix_partial[i];
  delete[] bfix_partial;
  
  cerr << " NULL_SPACE";
  fflush(stderr);
  
  // pull original solvee.
  const T norm_w(internal::sqrt(w.dot(w)));
  for(int i = 0, idx = 0; i < rvec.size(); i ++)
    if(internal::abs(w[i]) <= norm_w * err_singular)
      rvec[i]  = T(0);
    else
      rvec[i] /= w[i];
  rvec = V * rvec * bcr;
  for(int i = 0; i < acr.size(); i ++)
    if(acr[i] != T(0))
      rvec[i] /= acr[i];
  
  cerr << " INVERT" << endl;
  return isErrorMargin(A, b, rvec, true);
}

template <typename T, typename S> bool LP<T,S>::optimizeNullSpace(bool* fix_partial, Vec& rvec, const Mat& P, Vec& b, const int& cidx) const
{
  if(P.rows() != b.size() || P.cols() <= 0 || P.rows() < P.cols()) {
    cerr << " NullSpace: internal error (form does not match: P[x 0] <= b, <c,x>.)" << endl;
    return false;
  }
  
  const Vec bb(b - P * (P.transpose() * b));
  const T   normbb(internal::sqrt(bb.dot(bb)));
  bool*     checked = new bool[b.size()];
  VecS      one(b.size());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
  for(int i = 0; i < one.size(); i ++)
    one[i] = S(1);
  
  int  n_fixed;
  bool proper = false;
  Vec  bdash(b.size());
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < cidx + 1; i ++)
    bdash[i] = T(0);
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = cidx + 1; i < bdash.size(); i ++)
    bdash[i] = internal::sqrt(P.row(i).dot(P.row(i)));
  fflush(stderr);
  
  if(P.row(cidx).dot(P.row(cidx)) <= T(0)) {
    const T t(steps(fix_partial, checked, rvec, P, b, bdash, normbb, one));
    (void)gainVectors(proper, fix_partial, checked, rvec, P, b + bdash * t, normbb, one, false, n_fixed);
    delete[] checked;
    return isErrorMargin(P, b + bdash * t, rvec, true);
  }
  
  // get optimal value.
  //  note: for this step, dynamic accuracy detecting will earn calculation time.
  bool gained   = false;
  bool feasible = false;
  T    mid(normbb * err_opt);
  T    step(pow(err_opt, - T(2)));
  T    b_proper(mid / step);
  b[cidx] = T(0);
  const T t0(steps(fix_partial, checked, rvec, P, b, bdash, normbb, one));
  if(gainVectors(proper, fix_partial, checked, rvec, P, b + bdash * t0, normbb, one, true, n_fixed) && proper) {
    feasible = true;
    mid     *= - T(1);
    b_proper = T(0);
  }
  int n_steps = 2 + 2 * int(internal::log(- T(2) * internal::log(err_raw_epsilon) / internal::log(T(2))) / internal::log(T(2)));
  for(int i = 0; i <= n_steps; i ++) {
    step    = internal::sqrt(step);
    b[cidx] = mid * step;
    if(!isfinite(b[cidx]))
      continue;
    const T t(steps(fix_partial, checked, rvec, P, b, bdash, normbb, one));
    gained  = gainVectors(proper, fix_partial, checked, rvec, P, b + bdash * t, normbb, one, true, n_fixed);
    if(proper) {
      feasible = true;
      b_proper = b[cidx];
      mid     *= step;
    }
    cerr << (const char*)(proper ? "." : "!");
  }
  cerr << endl;
  
  b[cidx] = b_proper;
  const T t(steps(fix_partial, checked, rvec, P, b, bdash, normbb, one));
  (void)gainVectors(proper, fix_partial, checked, rvec, P, b + bdash * t, normbb, one, false, n_fixed);
  delete[] checked;
  return isErrorMargin(P, b + bdash * t, rvec, true);
}

template <typename T, typename S> void LP<T,S>::normalizeRows(Mat& A, Vec& b) const
{
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < A.rows(); i ++) {
    T min(internal::abs(A(i, 0)));
    T max(min);
    for(int j = 1; j < A.cols(); j ++) {
      const T buf(internal::abs(A(i, j)));
      if(max < buf)
        max = buf;
      if((buf < min || min <= T(0)) && buf != T(0))
        min = buf;
    }
    const T buf(internal::abs(b[i]));
    if(max < buf)
      max = buf;
    if((buf < min || min <= T(0)) && buf != T(0))
      min = buf;
    const T ratio(internal::sqrt(max * min));
    if(ratio > T(0)) {
      A.row(i) /= ratio;
      b[i]     /= ratio;
    }
  }
  return;
}

template <typename T, typename S> void LP<T,S>::normalizeCols(Vec& acr, T& bcr, Mat& A, Vec& b) const
{
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < A.cols(); i ++) {
    T min(internal::abs(A(0, i)));
    T max(min);
    for(int j = 1; j < A.rows(); j ++) {
      const T buf(internal::abs(A(j, i)));
      if(max < buf)
        max = buf;
      if((buf < min || min <= T(0)) && buf != T(0))
        min = buf;
    }
    const T ratio(internal::sqrt(min * max));
    if(ratio > T(0)) {
      acr[i]   *= ratio;
      A.col(i) /= ratio;
    }
  }
  {
    T min(internal::abs(b[0]));
    T max(min);
    for(int j = 1; j < A.rows(); j ++) {
      const T buf(internal::abs(b[j]));
      if(max < buf)
        max = buf;
      if((buf < min || min <= T(0)) && buf != T(0))
        min = buf;
    }
    const T ratio(internal::sqrt(min * max));
    if(ratio > T(0)) {
      bcr *= ratio;
      b   /= ratio;
    }
  }
  return;
}

template <typename T, typename S> T LP<T,S>::errorCheck(const Mat& A, const Vec& b) const
{
  // accuracy check.
  T max(1), min(1);
  for(int i = 0; i < A.rows(); i ++) {
    for(int j = 0; j < A.cols(); j ++) {
      const T buf(internal::abs(A(i, j)));
      if(buf > T(0)) {
        if(buf <= min)
          min = buf;
        if(max <= buf)
          max = buf;
      }
    }
    const T buf(internal::abs(b[i]));
    if(buf > T(0)) {
      if(buf <= min)
        min = buf;
      if(max <= buf)
        max = buf;
    }
  }
  return (min / max) / internal::sqrt(err_error);
}

template <typename T, typename S> T LP<T,S>::steps(bool* fix, bool* checked, Vec& rvec, const Mat& P, const Vec& b, const Vec& bdash, T normbb, const VecS& one) const {
  const T normbd(internal::sqrt(bdash.dot(bdash)));
  int  n_fixed;
  bool proper = false;
  bool gained = false;
  T    mid(max(normbd, normbb / normbd) * err_feas_step * err_feas_step * std::max(T(P.rows()) * T(P.cols()), err_fstep_max));
  T    step(pow(err_feas_step, - T(2)));
  T    b_proper(mid);
  int  n_steps = 2 + 2 * int(internal::log(- internal::log(err_raw_epsilon) / internal::log(T(2))) / internal::log(T(2)));
  for(int i = 0; i < n_steps; i ++) {
    step = internal::sqrt(step);
    const T middle(mid * step);
    if(!isfinite(middle))
      cerr << "inf";
    gained = gainVectors(proper, fix, checked, rvec, P, b + bdash * middle, normbb, one, true, n_fixed);
    if(proper) {
      b_proper = middle;
      mid     *= step;
    }
    cerr << (const char*)(proper ? "." : "!");
  }
  cerr << endl;
  return b_proper;
}

template <typename T, typename S> bool LP<T,S>::gainVectors(bool& proper, bool* fix, bool* checked, Vec& rvec, const Mat& P, Vec b, const T& orignormbb, const VecS& one, const bool& only_check, int& n_fixed) const
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
  if(internal::sqrt(bb.dot(bb)) <= threshold_feas * internal::sqrt(b.dot(b))) {
    for(int i = 0; i < bb.size() - P.cols() * 2 - 1; i ++)
      bb[i] = - T(1);
    for(int i = bb.size() - P.cols() * 2 - 1; i < bb.size(); i ++)
      bb[i] = b[i];
    const Vec bbb(bb - P * (P.transpose() * bb));
    if(internal::sqrt(bbb.dot(bbb)) <= threshold_feas * internal::sqrt(bb.dot(bb))) {
      if(!only_check) {
        rvec  = P.transpose() * (bbb / err_error + b + bb);
        cerr << " Second trivial matrix.";
      }
      return true;
    }
  }
  
  const T normbb(internal::sqrt(bb.dot(bb)));
  bb /= normbb;
  
  MatS Pverb(P.rows(), P.cols() + 1);
  Pverb.col(0) = - bb.template cast<S>();
#if defined(_OPENMP)
#pragma omp for
#endif
  for(int i = 0; i < P.cols(); i ++)
    Pverb.col(i + 1) = P.col(i).template cast<S>();
  
  (void)giantStep(fix, checked, Pverb, n_fixed, one);
  const Vec rvec0((Pverb.transpose() * (- one)).template cast<T>());
  rvec = Pverb.template cast<T>() * rvec0;
  for(int i = 0; i < rvec.size(); i ++)
    if(checked[i])
      rvec[i] = T(0);
  
  const bool proper0(isErrorMargin(Pverb, one * T(0), rvec0, false));
  const T    work(- bb.dot(rvec));
  const T    norm(internal::sqrt(rvec.dot(rvec)));
  rvec = P.transpose() * (rvec / work * normbb + b);
  
  if(!only_check)
    cerr << " giantStep(" << normbb << "," << work / norm << ", " << norm << ")";
  // cerr << endl;
  
  proper = isErrorMargin(P, b, rvec, false);
  return proper || proper0;
}

template <typename T, typename S> bool LP<T,S>::giantStep(bool* fix, bool* checked, MatS& Pverb, int& n_fixed, const VecS& one) const
{
  VecS norm(one.size());
  for( ; n_fixed < Pverb.cols() - 1; n_fixed ++) {
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for
#endif
    for(int j = 0; j < Pverb.rows(); j ++) {
      norm[j]     = internal::sqrt(Pverb.row(j).dot(Pverb.row(j)));
      checked[j] |= norm[j] <= threshold_p0;
    }
    
    // O(mn^2) check for inner or not.
    VecS on(Pverb * (Pverb.transpose() * (- one)));
    if(checkInner(on, getMax(checked, on.template cast<T>())))
      break;
#if defined(_OPENMP)
#pragma omp for
#endif
    for(int j = 0; j < on.size(); j ++)
      if(!checked[j])
        on[j] /= norm[j];
    const int max_idx = getMax(checked, on.template cast<T>());
    if(max_idx >= on.size()) {
      cerr << " GiantStep: no more direction.";
      break;
    }
    
    // O(mn^2) over all in this function.
    // this makes some error after some steps.
    const VecS orth(Pverb.row(max_idx));
    const S    norm2orth(orth.dot(orth));
#if defined(_OPENMP)
#pragma omp for
#endif
    for(int j = 0; j < Pverb.rows(); j ++)
      Pverb.row(j) -= (Pverb.row(j).dot(orth)) * orth / norm2orth;
    fix[max_idx] = Pverb.col(0).dot(Pverb.col(0)) <= err_norm(threshold_p0) ? - 1 : 1;
  }
  
  const VecS on(Pverb * (Pverb.transpose() * (- one)));
  return checkInner(on, getMax(checked, on.template cast<T>()));
}

template <typename T, typename S> bool LP<T,S>::checkInner(const VecS& on, const int& max_idx) const
{
  return max_idx < on.size() && on[max_idx] <= internal::sqrt(on.dot(on)) * threshold_inner;
}

template <typename T, typename S> int LP<T,S>::getMax(const bool* checked, const Vec& on) const
{
  int result = 0;
  for(; result < on.size(); result ++)
    if(!checked[result])
      break;
  for(int j = result + 1; j < on.size(); j ++)
    if(!checked[j] && on[result] < on[j])
      result = j;
  if(checked[result])
    result = on.size();
  return result;
}

template <typename T, typename S> bool LP<T,S>::isErrorMargin(const Mat& A, const Vec& b, const Vec& x, const bool& disp) const
{
  T result(0);
  const Vec err(A * x - b);
  for(int i = 0; i < b.size(); i ++)
    if(result < err[i]) result = err[i];
  if(disp)
    cerr << " errorMargin?(" << internal::sqrt(x.dot(x)) << ", " << internal::sqrt(b.dot(b)) << ", " << result << ")";
  return internal::isfinite(x.dot(x)) && (x.dot(x) > err_norm(err_error) ? result / internal::sqrt(x.dot(x)) : result) <= err_error;
}

template <typename T, typename S> bool LP<T,S>::svdschur(const Mat& A, Mat& U, Vec& w, Mat& V) const
{
  if(A.cols() > A.rows()) {
    cerr << " svdschur : fatal error." << endl;
    return false;
  }
  Eigen::RealSchur<Mat> schur;
  schur.compute(A.transpose() * A, true);
  V = schur.matrixU();
  U = A * V;
  w = Vec(U.cols());
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(int i = 0; i < U.cols(); i ++) {
    w[i] = internal::sqrt(U.col(i).dot(U.col(i)));
    if(w[i] <= internal::sqrt(err_singular)) {
      w[i]      = T(0);
      for(int j = 0; j < U.rows(); j ++)
        U(j, i) = T(0);
      cerr << "Z";
    } else
      U.col(i) /= w[i];
  }
  return true;
}

#define _LINEAR_OPT_
#endif

