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
  
         Vec inner(const Mat& A, const Vec& b) const;
  inline Mat roughQR(const Mat& At) const;
  T epsilon;
};

template <typename T> inline Linner<T>::Linner() {
  epsilon = pow(T(2), - T(16));
  return;
}

template <typename T> inline Linner<T>::~Linner() {
  return;
}

template <typename T> typename Linner<T>::Vec Linner<T>::inner(const Mat& A, const Vec& b) const {
  assert(A.rows() == b.size() && 0 < A.cols() && 0 < b.size());
  // cout << A << endl << b << endl << c.transpose() << endl;
  cerr << " (" << A.rows() << ", " << A.cols() << ")";
  {
    bool trivial(true);
    for(int i = 0; i < b.size() && trivial; i ++)
      trivial = T(0) <= b[i];
    if(trivial) {
      Vec rvec(A.cols());
      for(int i = 0; i < rvec.size(); i ++)
        rvec[i] = T(0);
      return rvec;
    }
  }
  Mat AA(A.rows() + A.cols() * 2 + 1, A.cols() + 1);
  for(int i = 0; i < A.rows(); i ++) {
    for(int j = 0; j < A.cols(); j ++)
      AA(i, j) = A(i, j);
    AA(i, A.cols()) = - b[i];
    AA.row(i) /= sqrt(AA.row(i).dot(AA.row(i)));
  }
  T m(1);
  for(int i = 0; i < A.rows(); i ++)
    for(int j = 0; j < AA.cols(); j ++) {
      const auto buf(abs(AA(i, j)));
      if(buf != T(0)) m = min(m, buf);
      assert(isfinite(buf) && !isnan(buf));
    }
  cerr << " scale(" << m << ")" << flush;
  for(int i = A.rows(); i < AA.rows(); i ++)
    for(int j = 0; j < AA.cols(); j ++)
      AA(i, j) = (i - A.rows()) / 2 == j
        ? T((i - A.rows()) & 1 ? 1 : - 1)
        : T(0);
        auto Pt(roughQR(AA.transpose()));
  cerr << "Q" << flush;
  const Mat  R(Pt * AA);
  cerr << "R" << flush;
  
  Vec one(Pt.cols());
#if defined(_OPENMP)
#pragma omp simd
#endif
  for(int i = 0; i < Pt.cols(); i ++)
    one[i] = T(1);
  for(int n_fixed = 0; n_fixed < Pt.rows() - 1; n_fixed ++) {
#if defined(_WITHOUT_EIGEN_)
    const auto on(Pt.projectionPt(- one));
#else
    const Vec  on(Pt.transpose() * (Pt * (- one)));
#endif
    auto fidx(0);
    for( ; fidx < on.size() && on[fidx] <= T(0); fidx ++) ;
    for(int i = fidx + 1; i < on.size(); i ++)
      if(T(0) < on[i] && on[i] < on[fidx]) fidx = i;
    if(! n_fixed) fidx = Pt.cols() - 1;
    if(on.size() <= fidx || on[fidx] <= T(0)) break;
    
    // O(mn^2) over all in this function.
    const Vec  orth(Pt.col(fidx));
    const auto norm2orth(orth.dot(orth) + epsilon);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
#if defined(_WITHOUT_EIGEN_)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / norm2orth);
#else
      Pt.col(j) -= orth * Pt.col(j).dot(orth) / norm2orth;
#endif
  }
  cerr << "G" << flush;
#if defined(_WITHOUT_EIGEN_)
  auto rvec(- R.solve(Pt * one));
#else
  auto rvec(- R.inverse() * (Pt * one));
#endif
  cerr << "I" << flush;
  Vec rrvec(rvec.size() - 1);
  for(int i = 0; i < rrvec.size(); i ++)
    rrvec[i] = rvec[i] / rvec[rvec.size() - 1];
  return rrvec;
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

