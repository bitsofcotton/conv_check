/* You can use one of the both BSD 3-Clause License or GNU Lesser General Public License 3.0 for this source. */
/* BSD 3-Clause License:
 * Copyright (c) 2013 - 2020, kazunobu watatsu.
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

template <typename T> class Linner {
public:
#if defined(_WITHOUT_EIGEN_)
  typedef SimpleMatrix<T> Mat;
  typedef SimpleVector<T> Vec;
#else
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
#endif
  
  inline Linner();
  inline ~Linner();
  
  inline bool inner(bool* fix_partial, const Mat& A, const Vec& b) const;
  inline bool inner(Vec& rvec, const Mat& A, const Vec& b) const;
         bool inner(bool* fix_partial, Vec& rvec, const Mat& A, const Vec& b) const;
  inline bool isErrorMargin(const Mat& A, const Vec& b, const Vec& x, const bool& disp = true) const;
  
private:
         bool gainVectors(bool* fix, bool* checked, Vec& rvec, const Mat& Pt, const Vec& b) const;
  inline Mat  roughQR(const Mat& At) const;

  T threshold_feas;
  T threshold_p0;
  T threshold_loop;
  T threshold_inner;
  T err_error;
  T large;
};

template <typename T> inline Linner<T>::Linner() {
  const auto epsilon(numeric_limits<T>::epsilon());
  // const auto epsilon(T(1) >> short(62));
  
  // error rate for orthogonalized b.
  threshold_feas    = pow(epsilon, T(5) / T(6));
  
  // if SimpleVector<T>::solve fails, increase this:
  threshold_p0      = pow(epsilon, T(4) / T(6));
  
  // is we gain inner or not:
  threshold_inner   = pow(epsilon, T(2) / T(6));
  
  // last stage error rate:
  err_error         = pow(epsilon, T(1) / T(6));
  
  // guarantee b is positive.
  large             = T(1) / sqrt(sqrt(err_error));
  
  // XXX: Please configure me first.
  // if threshold_loop < 0, extends some regions.
  threshold_loop    = T(0);
  
  assert(T(1) < large);
  return;
}

template <typename T> inline Linner<T>::~Linner() {
  return;
}

template <typename T> inline bool Linner<T>::inner(bool* fix_partial, const Mat& A, const Vec& b) const {
  Vec rvec;
  return inner(fix_partial, rvec, A, b);
}

template <typename T> inline bool Linner<T>::inner(Vec& rvec, const Mat& A, const Vec& b) const {
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
  T M(1);
  T m(1);
  for(int i = 0; i < A.rows(); i ++) {
    for(int j = 0; j < AA.cols(); j ++) {
      const T buf(abs(AA(i, j)));
      if(buf != T(0)) m = min(m, buf);
      M = max(M, buf);
      assert(isfinite(buf) && !isnan(buf));
    }
    const T buf(abs(bb[i]));
    if(buf != T(0)) m = min(m, buf);
    M = max(M, buf);
    assert(isfinite(buf) && !isnan(buf));
  }
  cerr << " err_error(" << sqrt(err_error) * M / m << ")" << flush;
  const auto Mm(sqrt(M * m));
  for(int i = A.rows(); i < AA.rows(); i ++) {
    for(int j = 0; j < A.cols(); j ++)
      AA(i, j) = (i - A.rows()) / 2 == j
        ? ((i - A.rows()) % 2 ? Mm : - Mm)
        : T(0);
    bb[i] = large * m;
  }
  bool* bfix_partial = new bool[AA.rows()];
  bool* checked = new bool[AA.rows()];
  const auto Pt(roughQR(AA.transpose()));
  cerr << "Q" << flush;
  const Mat  R(Pt * AA);
  cerr << "R" << flush;
  const bool f_inner(gainVectors(bfix_partial, checked, rvec, Pt, bb));
  if(fix_partial)
    for(int i = 0; i < A.rows(); i ++)
      fix_partial[i] = bfix_partial[i];
  delete[] bfix_partial;
  delete[] checked;
#if defined(_WITHOUT_EIGEN_)
  rvec = R.solve(rvec);
#else
  rvec = R.inverse() * rvec;
#endif
  cerr << "I" << flush;
  return f_inner || isErrorMargin(A, b, rvec, true);
}

template <typename T> bool Linner<T>::gainVectors(bool* fix, bool* checked, Vec& rvec, const Mat& Pt, const Vec& b) const {
  rvec.resize(Pt.rows());
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
#if defined(_WITHOUT_EIGEN_)
  Vec bb(b - Pt.projectionPt(b));
#else
  Vec bb(b - Pt.transpose() * (Pt * b));
#endif
  if(sqrt(bb.dot(bb)) <= threshold_feas * sqrt(b.dot(b))) {
    for(int i = 0; i < Pt.cols() - Pt.rows() * 2; i ++)
      bb[i] = sqrt(Pt.col(i).dot(Pt.col(i)));
    for(int i = Pt.cols() - Pt.rows() * 2; i < Pt.cols(); i ++)
      bb[i] = b[i];
#if defined(_WITHOUT_EIGEN_)
    const Vec bbb(bb - Pt.projectionPt(bb));
#else
    const Vec bbb(bb - Pt.transpose() * (Pt * bb));
#endif
    if(sqrt(bbb.dot(bbb)) <= threshold_feas * sqrt(bb.dot(bb))) {
      rvec = Pt * (b - bb);
      cerr << "t" << flush;
      return isErrorMargin(Pt.transpose(), b, rvec, ! false);
    }
    for(int i = 0; i < bb.size(); i ++)
      bb[i] = T(0);
    std::cerr << "0" << flush;
  }
  Vec mbb(- bb);
  const auto normb00(sqrt(mbb.dot(mbb)));
  const auto normb0(normb00 == T(0) ? T(1) : normb00);
  Mat Pverb(Pt);
  Vec on;
  Vec deltab;
  Vec orth;
  Vec norm(one.size());
  T   ratiob(1);
  for(n_fixed = 0; n_fixed < Pverb.rows(); n_fixed ++) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pverb.cols(); j ++) {
      norm[j]    = sqrt(Pverb.col(j).dot(Pverb.col(j)));
      checked[j] = fix[j] || norm[j] <= threshold_p0;
    }
    // extend b with threshold_loop. N.B. mbb = - b'.
    Vec mb(mbb + norm / sqrt(norm.dot(norm)) * normb0 * threshold_loop);
#if defined(_WITHOUT_EIGEN_)
    mb -= (deltab = Pverb.projectionPt(mb));
#else
    mb -= (deltab = Pverb.transpose() * (Pverb * mb));
#endif
    mb /= (ratiob = sqrt(mb.dot(mb)));
    
    // O(mn^2) check for inner or not.
#if defined(_WITHOUT_EIGEN_)
    on = Pverb.projectionPt(- one) + mb * mb.dot(- one);
#else
    on = Pverb.transpose() * (Pverb * (- one)) + mb * mb.dot(- one);
#endif
  
    // N.B. This might be rewrited as fix once with sorted fidxs.
    //   When the hypothesis is true, next loop must fix or get inner vector.
    //   Also, it takes O(lg(mn)) time order for this function.
    //   But to do so, theoretical proof is needed and it seems not
    //   because we use mb with little extended and its norm is in error order.
    int  fidx(0);
    for( ; fidx < on.size(); fidx ++)
      if(!checked[fidx])
        break;
    for(int j = fidx + 1; j < on.size(); j ++)
      if(!checked[j] && on[fidx] / norm[fidx] < on[j] / norm[j])
        fidx = j;
    if(fidx >= one.size())
      assert(0 && "rank is not full : should not be reached");
    on *= sqrt(norm.dot(norm)) / abs(mb.dot(on)) / norm[fidx];
    if(on[fidx] <= threshold_inner)
      break;
    
    // O(mn^2) over all in this function.
    orth = Pverb.col(fidx);
    const auto norm2orth(orth.dot(orth));
    const auto mbb0(mbb[fidx]);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pverb.cols(); j ++) {
      const auto work(Pverb.col(j).dot(orth) / norm2orth);
#if defined(_WITHOUT_EIGEN_)
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
        const auto ratio(sqrt(Pt.col(i).dot(Pt.col(i)) + b[i] * b[i]));
        F.row(j) = Pt.col(i) / ratio;
        f[j]     = b[i]      / ratio + threshold_loop;
        j ++;
      }
#if defined(_WITHOUT_EIGEN_)
    rvec = F.solve(f);
#else
    rvec = F.inverse() * f;
#endif
    cerr << "F" << flush;
  } else {
    rvec = Pt * (on * ratiob + deltab + b);
    cerr << "G" << flush;
  }
  return isErrorMargin(Pt.transpose(), b, rvec, ! false);
}

template <typename T> inline bool Linner<T>::isErrorMargin(const Mat& A, const Vec& b, const Vec& x, const bool& disp) const {
  T result(0);
  const Vec err(A * x - b);
  for(int i = 0; i < b.size(); i ++) {
    const auto& lerr(err[i]);
    if(!isfinite(lerr) || isnan(lerr) ||
       result < lerr) result = lerr;
  }
  if(disp)
    cerr << " errorMargin?(" << result / max(err_error, sqrt(x.dot(x))) << ")";
  return isfinite(result)   && !isnan(result)   &&
         isfinite(x.dot(x)) && !isnan(x.dot(x)) &&
         ((x.dot(x) == T(0) && result == T(0)) ||
          (x.dot(x) != T(0) && result / sqrt(x.dot(x)) <= err_error));
}

template <typename T> inline typename Linner<T>::Mat Linner<T>::roughQR(const Mat& At) const {
  Mat Q(At.rows(), At.cols());
  for(int i = 0; i < Q.rows(); i ++)
    for(int j = 0; j < Q.cols(); j ++)
      Q(i, j) = T(0);
  for(int i = 0; i < At.rows(); i ++) {
    // N.B. in this case, At is full rank.
#if defined(_WITHOUT_EIGEN_)
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

